#!/usr/bin/env python3
"""Phase 1C genome-wide methylation feasibility analysis.

Reads local Phase 1/1B derivative files and writes only to phase1C_output.
No local study matrix is transmitted outside the machine.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


CORE_AP = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-DRA",
    "HLA-DRB1",
    "B2M",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "NLRC5",
]
REGIONS = {
    "TSS200": {"TSS200"},
    "TSS1500": {"TSS1500"},
    "promoter_union": {"TSS200", "TSS1500", "5'UTR", "1stExon"},
}
MIN_MATCHED_N = 100
MIN_REGRESSION_N = 100
MIN_RNA_FINITE_FRACTION = 0.80
MIN_RNA_DETECTED_FRACTION = 0.20
BOOTSTRAP_REPS = 200
BOOTSTRAP_NEGATIVE_FRACTION = 0.70
SEED = 20260730


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "t", "1", "yes"})


def bh_fdr(values: pd.Series | np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    result = np.full(array.shape, np.nan, dtype=float)
    valid = np.isfinite(array)
    p = array[valid]
    if p.size == 0:
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def sample_code(sample_id: str) -> str:
    parts = str(sample_id).split("-")
    return parts[3][:2] if len(parts) >= 4 and len(parts[3]) >= 2 else ""


def patient_id(sample_id: str) -> str:
    return "-".join(str(sample_id).split("-")[:3])


def unique_patient_samples(samples: list[str], code: str) -> dict[str, str]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for sample in samples:
        if sample_code(sample) == code:
            grouped[patient_id(sample)].append(sample)
    return {patient: ids[0] for patient, ids in grouped.items() if len(ids) == 1}


def read_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return handle.readline().rstrip("\r\n").split("\t")


def paired_wilcoxon(delta: pd.Series) -> tuple[int, float, float]:
    clean = pd.to_numeric(delta, errors="coerce").dropna()
    if len(clean) < 3:
        return len(clean), float(clean.median()) if len(clean) else np.nan, np.nan
    if np.allclose(clean.to_numpy(), 0):
        return len(clean), float(clean.median()), 1.0
    result = stats.wilcoxon(clean, alternative="two-sided", zero_method="wilcox")
    return len(clean), float(clean.median()), float(result.pvalue)


def safe_spearman(x: pd.Series, y: pd.Series) -> tuple[int, float, float]:
    frame = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if (
        len(frame) < MIN_MATCHED_N
        or frame.iloc[:, 0].nunique() < 3
        or frame.iloc[:, 1].nunique() < 3
    ):
        return len(frame), np.nan, np.nan
    result = stats.spearmanr(frame.iloc[:, 0], frame.iloc[:, 1])
    return len(frame), float(result.statistic), float(result.pvalue)


def ols_methylation_cnv(
    y: pd.Series, methylation: pd.Series, cnv: pd.Series
) -> dict[str, float | int | str]:
    frame = pd.concat([y, methylation, cnv], axis=1)
    frame.columns = ["rna", "methylation", "cnv"]
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    default = {
        "regression_n": len(frame),
        "adjusted_methylation_coefficient": np.nan,
        "adjusted_methylation_se": np.nan,
        "adjusted_methylation_p": np.nan,
        "cnv_coefficient": np.nan,
        "cnv_se": np.nan,
        "cnv_p": np.nan,
        "regression_r2": np.nan,
        "regression_status": "insufficient_data",
    }
    if len(frame) < MIN_REGRESSION_N:
        return default
    x = np.column_stack(
        [
            np.ones(len(frame)),
            frame["methylation"].to_numpy(float),
            frame["cnv"].to_numpy(float),
        ]
    )
    if np.linalg.matrix_rank(x) < x.shape[1]:
        default["regression_status"] = "rank_deficient"
        return default
    yv = frame["rna"].to_numpy(float)
    beta, _, _, _ = np.linalg.lstsq(x, yv, rcond=None)
    residual = yv - x @ beta
    df = len(yv) - x.shape[1]
    sigma2 = float(residual @ residual / df)
    covariance = sigma2 * np.linalg.inv(x.T @ x)
    se = np.sqrt(np.diag(covariance))
    t_values = beta / se
    p_values = 2 * stats.t.sf(np.abs(t_values), df)
    total = float(((yv - yv.mean()) ** 2).sum())
    r2 = 1 - float((residual**2).sum()) / total if total > 0 else np.nan
    return {
        "regression_n": len(frame),
        "adjusted_methylation_coefficient": float(beta[1]),
        "adjusted_methylation_se": float(se[1]),
        "adjusted_methylation_p": float(p_values[1]),
        "cnv_coefficient": float(beta[2]),
        "cnv_se": float(se[2]),
        "cnv_p": float(p_values[2]),
        "regression_r2": r2,
        "regression_status": "complete",
    }


def stable_gene_seed(gene: str, variant: str) -> int:
    token = hashlib.sha256(f"{gene}|{variant}|{SEED}".encode()).hexdigest()[:8]
    return int(token, 16)


def bootstrap_negative_fraction(
    x: pd.Series, y: pd.Series, gene: str, variant: str
) -> tuple[int, float, float, float]:
    frame = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    n = len(frame)
    if n < MIN_MATCHED_N:
        return n, np.nan, np.nan, np.nan
    xr = stats.rankdata(frame.iloc[:, 0].to_numpy(float), method="average")
    yr = stats.rankdata(frame.iloc[:, 1].to_numpy(float), method="average")
    rng = np.random.default_rng(stable_gene_seed(gene, variant))
    indices = rng.integers(0, n, size=(BOOTSTRAP_REPS, n))
    xb = xr[indices]
    yb = yr[indices]
    xb = xb - xb.mean(axis=1, keepdims=True)
    yb = yb - yb.mean(axis=1, keepdims=True)
    denominator = np.sqrt((xb * xb).sum(axis=1) * (yb * yb).sum(axis=1))
    correlations = np.divide(
        (xb * yb).sum(axis=1),
        denominator,
        out=np.full(BOOTSTRAP_REPS, np.nan),
        where=denominator > 0,
    )
    valid = correlations[np.isfinite(correlations)]
    if valid.size == 0:
        return n, np.nan, np.nan, np.nan
    return (
        n,
        float(np.mean(valid < 0)),
        float(np.quantile(valid, 0.025)),
        float(np.quantile(valid, 0.975)),
    )


def load_expression(project: Path, genes: set[str]) -> pd.DataFrame:
    processed = project / "data_processed" / "TCGA_LIHC"
    annotation = pd.read_csv(
        processed / "TCGA_LIHC_gene_annotation.tsv", sep="\t", dtype=str
    )
    annotation = annotation[annotation["gene_name"].isin(genes)]
    id_to_gene = dict(zip(annotation["gene_id"], annotation["gene_name"]))
    path = processed / "TCGA_LIHC_RNA_TPM_matrix.tsv.gz"
    samples = read_header(path)[1:]
    pieces = []
    for chunk in pd.read_csv(path, sep="\t", compression="infer", chunksize=5000):
        selected = chunk[chunk["gene_id"].astype(str).isin(id_to_gene)].copy()
        if selected.empty:
            continue
        values = selected[samples].apply(pd.to_numeric, errors="coerce").astype(
            np.float32
        )
        values.index = selected["gene_id"].map(id_to_gene).to_numpy()
        pieces.append(values)
    if not pieces:
        raise RuntimeError("No promoter-mapped genes were found in TCGA RNA")
    result = pd.concat(pieces).groupby(level=0).median()
    return np.log2(result + 1.0)


def load_cnv(project: Path, genes: set[str]) -> pd.DataFrame:
    path = project / "phase1B" / "TCGA_LIHC_gene_level_CNV.tsv.gz"
    pieces = []
    for chunk in pd.read_csv(path, sep="\t", compression="infer", chunksize=5000):
        selected = chunk[chunk["gene_name"].astype(str).isin(genes)]
        if selected.empty:
            continue
        samples = [
            col for col in selected.columns if col not in {"gene_id", "gene_name"}
        ]
        values = selected[samples].apply(pd.to_numeric, errors="coerce").astype(
            np.float32
        )
        values.index = selected["gene_name"].astype(str).to_numpy()
        pieces.append(values)
    if not pieces:
        return pd.DataFrame()
    return pd.concat(pieces).groupby(level=0).median()


def audit_core_coverage(
    project: Path, out: Path, qc: pd.DataFrame, mapping: pd.DataFrame
) -> pd.DataFrame:
    gene_to_any_probes: dict[str, set[str]] = defaultdict(set)
    for row in qc[["probe_id", "UCSC_RefGene_Name"]].itertuples(index=False):
        for gene in {
            item.strip()
            for item in str(row.UCSC_RefGene_Name).split(";")
            if item.strip() and item.strip().lower() != "nan"
        }:
            gene_to_any_probes[gene].add(str(row.probe_id))
    rows = []
    for gene in CORE_AP:
        gm = mapping[mapping["gene"].eq(gene)].copy()
        any_probes = gene_to_any_probes.get(gene, set())
        promoter_ids = set(gm["probe_id"].astype(str))
        promoter_qc = qc[qc["probe_id"].astype(str).isin(promoter_ids)].copy()
        counts = {}
        for region in ["TSS200", "TSS1500", "5'UTR", "1stExon"]:
            region_map = gm[gm["region"].eq(region)]
            counts[f"{region}_probe_n"] = int(region_map["probe_id"].nunique())
            counts[f"{region}_base_pass_n"] = int(
                region_map.loc[
                    truthy(region_map["base_filter_pass"]), "probe_id"
                ].nunique()
            )
        tss200 = gm[gm["region"].eq("TSS200")]
        if tss200.empty:
            limitation = "no_TSS200_probe_in_HM450_annotation"
        elif not truthy(tss200["base_filter_pass"]).any():
            if truthy(tss200["common_snp_flag"]).all():
                limitation = "TSS200_probes_exist_but_all_common_SNP_flagged"
            elif truthy(tss200["zhou_general_mask"]).all():
                limitation = "TSS200_probes_exist_but_all_crossreactivity_masked"
            else:
                limitation = "TSS200_probes_exist_but_all_fail_prespecified_QC"
        else:
            limitation = "TSS200_covered_after_QC"
        expected_phase1b = counts["TSS200_base_pass_n"] > 0
        rows.append(
            {
                "gene": gene,
                "any_HM450_gene_annotated_probe": bool(any_probes),
                "any_HM450_gene_annotated_probe_n": len(any_probes),
                **counts,
                "promoter_probe_failed_missingness_n": int(
                    promoter_qc.loc[
                        ~truthy(promoter_qc["missingness_pass"]), "probe_id"
                    ].nunique()
                ),
                "promoter_probe_failed_sex_chromosome_n": int(
                    promoter_qc.loc[
                        truthy(promoter_qc["sex_chromosome"]), "probe_id"
                    ].nunique()
                ),
                "promoter_probe_failed_common_SNP_n": int(
                    promoter_qc.loc[
                        truthy(promoter_qc["common_snp_flag"]), "probe_id"
                    ].nunique()
                ),
                "promoter_probe_crossreactivity_flag_n": int(
                    promoter_qc.loc[
                        truthy(promoter_qc["MASK_general"]), "probe_id"
                    ].nunique()
                ),
                "promoter_probe_unmapped_position_n": int(
                    promoter_qc.loc[
                        ~truthy(promoter_qc["mapped_position"]), "probe_id"
                    ].nunique()
                ),
                "multiple_transcript_region_annotation": bool(
                    (pd.to_numeric(gm["annotation_region_count"], errors="coerce") > 1).any()
                    or gm["region"].nunique() > 1
                ),
                "regions_observed": ";".join(sorted(gm["region"].dropna().unique())),
                "primary_limitation": limitation,
                "expected_TSS200_gene_beta_from_audit": expected_phase1b,
                "coverage_interpretation": (
                    "HM450_design_or_prespecified_probe_QC_limitation"
                ),
            }
        )
    audit = pd.DataFrame(rows)
    phase1b_beta = pd.read_csv(
        project / "phase1B" / "gene_promoter_beta_TSS200.tsv.gz",
        sep="\t",
        index_col=0,
    )
    phase1b_present = phase1b_beta.notna().any(axis=1).to_dict()
    audit["phase1B_TSS200_gene_beta_present"] = audit["gene"].map(
        phase1b_present
    ).fillna(False)
    audit["script_mapping_discrepancy"] = (
        audit["expected_TSS200_gene_beta_from_audit"]
        != audit["phase1B_TSS200_gene_beta_present"]
    )
    audit.to_csv(out / "AP_core_gene_probe_coverage_audit.tsv", sep="\t", index=False)
    return audit


def load_selected_beta(path: Path, wanted: set[str]) -> pd.DataFrame:
    pieces = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        index_col=0,
        chunksize=4000,
        low_memory=False,
    ):
        selected = chunk.loc[chunk.index.astype(str).isin(wanted)]
        if selected.empty:
            continue
        selected.index = selected.index.astype(str)
        pieces.append(
            selected.apply(pd.to_numeric, errors="coerce").astype(np.float32)
        )
    if not pieces:
        raise RuntimeError("No eligible promoter probes were loaded")
    result = pd.concat(pieces)
    return result.loc[~result.index.duplicated(keep="first")]


def aggregate_promoters(
    project: Path,
    out: Path,
    mapping: pd.DataFrame,
    expression_genes: set[str],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    eligible = mapping[
        truthy(mapping["base_filter_pass"])
        & mapping["gene"].isin(expression_genes)
        & mapping["region"].isin(REGIONS["promoter_union"])
    ].copy()
    wanted = set(eligible["probe_id"].astype(str))
    beta = load_selected_beta(
        project
        / "data_processed"
        / "TCGA_LIHC"
        / "TCGA_LIHC_methylation_beta_matrix.tsv.gz",
        wanted,
    )
    probe_to_index = {probe: i for i, probe in enumerate(beta.index)}
    beta_values = beta.to_numpy(dtype=np.float32, copy=False)
    matrices: dict[str, pd.DataFrame] = {}
    report_rows = []
    for variant, regions in REGIONS.items():
        selected = eligible[eligible["region"].isin(regions)]
        groups = selected.groupby("gene")["probe_id"].agg(
            lambda values: sorted(set(map(str, values)))
        )
        genes = []
        rows = []
        for gene, ids in groups.items():
            indices = [probe_to_index[p] for p in ids if p in probe_to_index]
            if not indices:
                continue
            values = np.nanmedian(beta_values[indices, :], axis=0)
            genes.append(gene)
            rows.append(values.astype(np.float32))
            report_rows.append(
                {
                    "promoter_definition": variant,
                    "gene": gene,
                    "eligible_probe_n": len(indices),
                    "eligible_probe_ids": ";".join(ids),
                    "aggregation": "median",
                    "sample_n": beta.shape[1],
                    "nonmissing_sample_n": int(np.isfinite(values).sum()),
                }
            )
        matrix = pd.DataFrame(rows, index=genes, columns=beta.columns)
        matrix.index.name = "gene"
        matrices[variant] = matrix
        matrix.to_csv(
            out / f"gene_promoter_beta_genomewide_{variant}.tsv.gz",
            sep="\t",
            compression="gzip",
        )
    report = pd.DataFrame(report_rows)
    report.to_csv(
        out / "genomewide_promoter_aggregation_audit.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    return matrices, report


def series_by_patient(
    matrix: pd.DataFrame, gene: str, samples: dict[str, str]
) -> pd.Series:
    if gene not in matrix.index:
        return pd.Series(dtype=float)
    return pd.Series(
        {
            patient: matrix.at[gene, sample]
            for patient, sample in samples.items()
            if sample in matrix.columns
        },
        dtype=float,
    )


def paired_delta(
    matrix: pd.DataFrame,
    gene: str,
    tumors: dict[str, str],
    normals: dict[str, str],
) -> pd.Series:
    if gene not in matrix.index:
        return pd.Series(dtype=float)
    patients = sorted(set(tumors) & set(normals))
    return pd.Series(
        {
            patient: matrix.at[gene, tumors[patient]]
            - matrix.at[gene, normals[patient]]
            for patient in patients
            if tumors[patient] in matrix.columns
            and normals[patient] in matrix.columns
        },
        dtype=float,
    )


def compute_variant_associations(
    variant: str,
    beta: pd.DataFrame,
    expression: pd.DataFrame,
    cnv: pd.DataFrame,
    probe_counts: dict[str, int],
) -> tuple[pd.DataFrame, dict]:
    meth_tumor = unique_patient_samples(list(beta.columns), "01")
    meth_normal = unique_patient_samples(list(beta.columns), "11")
    rna_tumor = unique_patient_samples(list(expression.columns), "01")
    rna_normal = unique_patient_samples(list(expression.columns), "11")
    cnv_tumor = unique_patient_samples(list(cnv.columns), "01")
    paired_primary = sorted(set(meth_tumor) & set(rna_tumor))
    regression_patients = sorted(set(paired_primary) & set(cnv_tumor))
    rows = []
    genes = sorted(set(beta.index) & set(expression.index))
    for counter, gene in enumerate(genes, start=1):
        meth_primary = series_by_patient(
            beta, gene, {p: meth_tumor[p] for p in paired_primary}
        )
        rna_primary = series_by_patient(
            expression, gene, {p: rna_tumor[p] for p in paired_primary}
        )
        finite_fraction = float(rna_primary.notna().mean())
        detected_fraction = float((rna_primary.fillna(0) > 0).mean())
        expression_eligible = (
            finite_fraction >= MIN_RNA_FINITE_FRACTION
            and detected_fraction >= MIN_RNA_DETECTED_FRACTION
        )
        corr_n, rho, corr_p = (
            safe_spearman(meth_primary, rna_primary)
            if expression_eligible
            else (int(pd.concat([meth_primary, rna_primary], axis=1).dropna().shape[0]), np.nan, np.nan)
        )
        boot_n, boot_fraction, boot_low, boot_high = (
            bootstrap_negative_fraction(
                meth_primary, rna_primary, gene, variant
            )
            if expression_eligible and corr_n >= MIN_MATCHED_N
            else (corr_n, np.nan, np.nan, np.nan)
        )
        meth_delta = paired_delta(beta, gene, meth_tumor, meth_normal)
        meth_pair_n, meth_median, meth_p = paired_wilcoxon(meth_delta)
        rna_delta = paired_delta(expression, gene, rna_tumor, rna_normal)
        rna_pair_n, rna_median, rna_p = paired_wilcoxon(rna_delta)
        if gene in cnv.index:
            model = ols_methylation_cnv(
                series_by_patient(
                    expression,
                    gene,
                    {p: rna_tumor[p] for p in regression_patients},
                ),
                series_by_patient(
                    beta,
                    gene,
                    {p: meth_tumor[p] for p in regression_patients},
                ),
                series_by_patient(
                    cnv,
                    gene,
                    {p: cnv_tumor[p] for p in regression_patients},
                ),
            )
            cnv_available = True
        else:
            model = {
                "regression_n": 0,
                "adjusted_methylation_coefficient": np.nan,
                "adjusted_methylation_se": np.nan,
                "adjusted_methylation_p": np.nan,
                "cnv_coefficient": np.nan,
                "cnv_se": np.nan,
                "cnv_p": np.nan,
                "regression_r2": np.nan,
                "regression_status": "gene_level_CNV_missing",
            }
            cnv_available = False
        rows.append(
            {
                "gene": gene,
                "promoter_definition": variant,
                "eligible_probe_n": int(probe_counts.get(gene, 0)),
                "rna_finite_fraction": finite_fraction,
                "rna_detected_fraction": detected_fraction,
                "rna_detection_eligible": expression_eligible,
                "methylation_rna_n": corr_n,
                "methylation_rna_rho": rho,
                "correlation_p": corr_p,
                "bootstrap_n": boot_n,
                "bootstrap_reps": BOOTSTRAP_REPS,
                "bootstrap_negative_fraction": boot_fraction,
                "bootstrap_rho_q025": boot_low,
                "bootstrap_rho_q975": boot_high,
                "methylation_paired_n": meth_pair_n,
                "tumor_normal_methylation_delta": meth_median,
                "methylation_p": meth_p,
                "rna_paired_n": rna_pair_n,
                "rna_tumor_normal_delta": rna_median,
                "rna_tumor_normal_p": rna_p,
                "cnv_available": cnv_available,
                **model,
                "cnv_value_definition": "GDC ASCAT3 absolute gene copy_number",
                "missing_value_rule": "complete cases per analysis; no imputation",
            }
        )
        if counter % 1000 == 0:
            print(f"{variant}: processed {counter:,}/{len(genes):,} genes", flush=True)
    result = pd.DataFrame(rows)
    for p_col, fdr_col in [
        ("correlation_p", "correlation_fdr"),
        ("adjusted_methylation_p", "adjusted_methylation_fdr"),
        ("cnv_p", "cnv_fdr"),
        ("methylation_p", "methylation_fdr"),
        ("rna_tumor_normal_p", "rna_tumor_normal_fdr"),
    ]:
        result[fdr_col] = bh_fdr(result[p_col])
    result["absolute_rho_ge_0.20"] = result["methylation_rna_rho"].abs() >= 0.20
    result["cnv_adjustment_direction_changed"] = (
        result["methylation_rna_rho"].lt(0)
        & result["adjusted_methylation_coefficient"].ge(0)
    )
    context = {
        "variant": variant,
        "methylation_tumor_unique_patients": len(meth_tumor),
        "rna_tumor_unique_patients": len(rna_tumor),
        "paired_primary_patients": len(paired_primary),
        "regression_patients": len(regression_patients),
        "methylation_tumor_normal_pairs": len(set(meth_tumor) & set(meth_normal)),
        "rna_tumor_normal_pairs": len(set(rna_tumor) & set(rna_normal)),
        "genes_tested": len(result),
    }
    return result, context


def classify_tiers(all_results: pd.DataFrame) -> pd.DataFrame:
    main = all_results[all_results["promoter_definition"].eq("TSS200")].copy()
    sensitivity = all_results[
        all_results["promoter_definition"].isin(["TSS1500", "promoter_union"])
    ][
        [
            "gene",
            "promoter_definition",
            "methylation_rna_rho",
            "correlation_fdr",
            "adjusted_methylation_coefficient",
            "adjusted_methylation_fdr",
        ]
    ].copy()
    sensitivity["negative_direction_both"] = (
        sensitivity["methylation_rna_rho"].lt(0)
        & sensitivity["adjusted_methylation_coefficient"].lt(0)
    )
    sens_summary = sensitivity.groupby("gene").agg(
        sensitivity_definitions_available=("methylation_rna_rho", "count"),
        sensitivity_negative_both_n=("negative_direction_both", "sum"),
    )
    main = main.merge(sens_summary, on="gene", how="left")
    main["sensitivity_definitions_available"] = (
        main["sensitivity_definitions_available"].fillna(0).astype(int)
    )
    main["sensitivity_negative_both_n"] = (
        main["sensitivity_negative_both_n"].fillna(0).astype(int)
    )
    main["sensitivity_any_negative_both"] = (
        main["sensitivity_negative_both_n"] >= 1
    )
    main["sensitivity_all_available_negative_both"] = (
        (main["sensitivity_definitions_available"] > 0)
        & (
            main["sensitivity_negative_both_n"]
            == main["sensitivity_definitions_available"]
        )
    )
    corr_support = main["methylation_rna_rho"].lt(0) & main[
        "correlation_fdr"
    ].lt(0.05)
    adjusted_direction = main["adjusted_methylation_coefficient"].lt(0)
    adjusted_support = adjusted_direction & main[
        "adjusted_methylation_fdr"
    ].lt(0.05)
    bootstrap_support = main["bootstrap_negative_fraction"].ge(
        BOOTSTRAP_NEGATIVE_FRACTION
    )
    sensitivity_support = main["sensitivity_any_negative_both"]
    support_count = (
        corr_support.astype(int)
        + adjusted_direction.astype(int)
        + adjusted_support.astype(int)
        + bootstrap_support.astype(int)
        + sensitivity_support.astype(int)
    )
    stable = (
        corr_support
        & adjusted_support
        & bootstrap_support
        & sensitivity_support
    )
    partial = (~stable) & main["methylation_rna_rho"].lt(0) & support_count.ge(3)
    exploratory = (
        (~stable)
        & (~partial)
        & (
            main["methylation_rna_rho"].lt(0)
            | main["adjusted_methylation_coefficient"].lt(0)
            | (
                main["tumor_normal_methylation_delta"].gt(0)
                & main["rna_tumor_normal_delta"].lt(0)
            )
        )
    )
    main["evidence_support_count_of_5"] = support_count
    main["genomewide_evidence_tier"] = np.select(
        [stable, partial, exploratory],
        ["Stable", "Partially stable", "Exploratory"],
        default="Unsupported",
    )
    main["stable_rule"] = (
        "rho<0; correlation_FDR<0.05; CNV-adjusted_methylation_beta<0; "
        "adjusted_FDR<0.05; bootstrap_negative_fraction>=0.70; "
        "negative direction in >=1 prespecified sensitivity definition"
    )
    return main


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    means = frame.mean(axis=1)
    stds = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(means, axis=0).div(stds, axis=0)


def score_ap_rna(
    out: Path, expression: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tumor_map = unique_patient_samples(list(expression.columns), "01")
    samples = [tumor_map[p] for p in sorted(tumor_map)]
    available = [gene for gene in CORE_AP if gene in expression.index]
    z = zscore_rows(expression.loc[available, samples])
    score = z.mean(axis=0, skipna=True)
    result = pd.DataFrame(
        {
            "patient_id": [patient_id(sample) for sample in samples],
            "sample_id": samples,
            "AP_RNA_output_score": score.reindex(samples).to_numpy(),
            "available_core_gene_n": z.notna().sum(axis=0).reindex(samples).to_numpy(),
            "available_core_genes": ";".join(available),
            "score_direction": "higher = higher antigen-presentation RNA output",
            "standardization": "gene-wise z-score across unique TCGA primary tumors",
        }
    )
    contributions = z.T.copy()
    contributions.index.name = "sample_id"
    contributions.to_csv(
        out / "intermediate" / "AP_RNA_gene_z_contributions.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    result.to_csv(out / "AP_RNA_output_score_TCGA.tsv", sep="\t", index=False)
    return result, z


def load_proteome(project: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = (
        project
        / "data_processed"
        / "CPTAC_HCC"
        / "CPTAC_HCC_proteome_logratio_matrix.tsv.gz"
    )
    header = read_header(path)
    metadata_columns = [
        "feature_type",
        "feature_id",
        "protein_id",
        "gene_symbol",
        "ncbi_gene_id",
        "authority",
        "description",
        "organism",
        "chromosome",
        "locus",
    ]
    samples = header[len(metadata_columns) :]
    pieces = []
    usecols = ["feature_id", "gene_symbol"] + samples
    for chunk in pd.read_csv(
        path, sep="\t", compression="infer", usecols=usecols, chunksize=5000
    ):
        selected = chunk[
            chunk["gene_symbol"].astype(str).isin(CORE_AP)
            & chunk["feature_id"].astype(str).ne("StdDev")
        ]
        if selected.empty:
            continue
        values = selected[samples].apply(pd.to_numeric, errors="coerce")
        values.index = selected["gene_symbol"].astype(str).to_numpy()
        pieces.append(values)
    protein = pd.concat(pieces).groupby(level=0).median()
    clinical = pd.read_csv(
        project
        / "data_processed"
        / "CPTAC_HCC"
        / "CPTAC_HCC_sample_clinical_mapping.tsv",
        sep="\t",
        dtype=str,
    )
    return protein, clinical


def score_ap_protein(
    project: Path, out: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    protein, clinical = load_proteome(project)
    available = [gene for gene in CORE_AP if gene in protein.index]
    z = zscore_rows(protein.loc[available])
    score = z.mean(axis=0, skipna=True)
    meta = clinical.set_index("sample_id")
    rows = []
    for sample in z.columns:
        rows.append(
            {
                "sample_id": sample,
                "case_id": meta.at[sample, "case_id"] if sample in meta.index else "",
                "sample_type": meta.at[sample, "sample_type"]
                if sample in meta.index
                else "",
                "AP_protein_output_score": float(score.get(sample, np.nan)),
                "available_core_gene_n": int(z[sample].notna().sum()),
                "available_core_genes": ";".join(available),
                "score_direction": "higher = higher antigen-presentation protein output",
                "standardization": "gene-wise z-score across CPTAC tumor and paratumor samples",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(out / "AP_protein_output_score_CPTAC.tsv", sep="\t", index=False)
    z.T.to_csv(
        out / "intermediate" / "AP_protein_gene_z_contributions.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    return result, z


def robustness_rows(
    score_name: str,
    full_score: pd.Series,
    z: pd.DataFrame,
    sample_group: pd.Series | None = None,
) -> list[dict]:
    rows = []
    for gene in z.index:
        loo = z.drop(index=gene).mean(axis=0, skipna=True)
        frame = pd.concat([full_score, loo], axis=1).dropna()
        rho = (
            stats.spearmanr(frame.iloc[:, 0], frame.iloc[:, 1]).statistic
            if len(frame) >= 10
            else np.nan
        )
        rows.append(
            {
                "score": score_name,
                "analysis": "leave_one_gene_out",
                "omitted_gene_or_subset": gene,
                "n": len(frame),
                "spearman_vs_full_score": rho,
                "median_score_difference": float(
                    (frame.iloc[:, 1] - frame.iloc[:, 0]).median()
                )
                if len(frame)
                else np.nan,
                "paired_group_delta": np.nan,
                "paired_group_p": np.nan,
            }
        )
    subsets = {
        "MHC_I": ["HLA-A", "HLA-B", "HLA-C", "B2M"],
        "MHC_II": ["HLA-DRA", "HLA-DRB1"],
        "processing_regulation": ["TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"],
    }
    for subset, genes in subsets.items():
        present = [gene for gene in genes if gene in z.index]
        if not present:
            continue
        subset_score = z.loc[present].mean(axis=0, skipna=True)
        frame = pd.concat([full_score, subset_score], axis=1).dropna()
        rho = (
            stats.spearmanr(frame.iloc[:, 0], frame.iloc[:, 1]).statistic
            if len(frame) >= 10
            else np.nan
        )
        rows.append(
            {
                "score": score_name,
                "analysis": "fixed_core_subset",
                "omitted_gene_or_subset": subset,
                "subset_genes_present": ";".join(present),
                "n": len(frame),
                "spearman_vs_full_score": rho,
                "median_score_difference": np.nan,
                "paired_group_delta": np.nan,
                "paired_group_p": np.nan,
            }
        )
    if sample_group is not None:
        joined = pd.concat([full_score.rename("score"), sample_group], axis=1).dropna()
        pairs = []
        for _, group in joined.groupby(sample_group.name):
            if len(group) == 2:
                labels = group.index.astype(str)
                tumor = [x for x in labels if x.startswith("T")]
                para = [x for x in labels if x.startswith("P")]
                if len(tumor) == 1 and len(para) == 1:
                    pairs.append(
                        group.at[tumor[0], "score"] - group.at[para[0], "score"]
                    )
        n, delta, p = paired_wilcoxon(pd.Series(pairs, dtype=float))
        rows.append(
            {
                "score": score_name,
                "analysis": "paired_tumor_minus_paratumor",
                "omitted_gene_or_subset": "none",
                "n": n,
                "spearman_vs_full_score": np.nan,
                "median_score_difference": np.nan,
                "paired_group_delta": delta,
                "paired_group_p": p,
            }
        )
    return rows


def score_robustness(
    out: Path,
    rna_score: pd.DataFrame,
    rna_z: pd.DataFrame,
    protein_score: pd.DataFrame,
    protein_z: pd.DataFrame,
) -> pd.DataFrame:
    rna_full = rna_score.set_index("sample_id")["AP_RNA_output_score"]
    protein_full = protein_score.set_index("sample_id")[
        "AP_protein_output_score"
    ]
    protein_cases = protein_score.set_index("sample_id")["case_id"]
    rows = robustness_rows("AP_RNA_output_score", rna_full, rna_z)
    rows.extend(
        robustness_rows(
            "AP_protein_output_score",
            protein_full,
            protein_z,
            sample_group=protein_cases,
        )
    )
    result = pd.DataFrame(rows)
    result.to_csv(out / "AP_score_robustness.tsv", sep="\t", index=False)
    return result


def write_estimate_input(out: Path, expression: pd.DataFrame) -> None:
    tumor_map = unique_patient_samples(list(expression.columns), "01")
    samples = [tumor_map[p] for p in sorted(tumor_map)]
    selected = expression[samples].copy()
    selected.index.name = "GeneSymbol"
    selected.to_csv(
        out / "intermediate" / "TCGA_primary_log2TPM_gene_symbol.tsv",
        sep="\t",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    project = args.project.resolve()
    out = project / "phase1C_output"
    (out / "intermediate").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    frozen = (
        out
        / "frozen_phase1B"
        / "HCC_multiomics_Phase1B_release_complete.zip"
    )
    if not frozen.exists():
        raise FileNotFoundError("Frozen Phase 1B release is missing")
    frozen_hash = sha256(frozen)
    expected_hash = "c9a464477f7c9c6328529cfc17282fba5e9166cc16a5b4b4d71daa0c922c61d9"
    if frozen_hash != expected_hash:
        raise RuntimeError("Frozen Phase 1B release hash mismatch")

    qc = pd.read_csv(
        project / "phase1B" / "methylation_probe_qc_annotated.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    mapping = pd.read_csv(
        project / "phase1B" / "promoter_probe_gene_mapping.tsv.gz",
        sep="\t",
        low_memory=False,
    )
    coverage = audit_core_coverage(project, out, qc, mapping)
    promoter_genes = set(mapping["gene"].dropna().astype(str))
    print(f"Loading RNA for {len(promoter_genes):,} promoter-mapped genes")
    expression = load_expression(project, promoter_genes)
    matrices, aggregation = aggregate_promoters(
        project, out, mapping, set(expression.index)
    )
    all_genes = set().union(*(set(frame.index) for frame in matrices.values()))
    expression = expression.reindex(sorted(set(expression.index) & all_genes))
    print(f"Loading CNV for {len(expression):,} RNA/promoter genes")
    cnv = load_cnv(project, set(expression.index))

    probe_counts = {
        variant: (
            aggregation[aggregation["promoter_definition"].eq(variant)]
            .set_index("gene")["eligible_probe_n"]
            .to_dict()
        )
        for variant in REGIONS
    }
    results = []
    contexts = []
    for variant in ["TSS200", "TSS1500", "promoter_union"]:
        result, context = compute_variant_associations(
            variant,
            matrices[variant],
            expression,
            cnv,
            probe_counts[variant],
        )
        results.append(result)
        contexts.append(context)
    all_results = pd.concat(results, ignore_index=True)
    all_results.to_csv(
        out / "genomewide_methylation_RNA_results.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    tiers = classify_tiers(all_results)
    tiers.to_csv(
        out / "genomewide_methylation_candidate_tiers.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(contexts).to_csv(
        out / "genomewide_analysis_sample_context.tsv", sep="\t", index=False
    )

    rna_score, rna_z = score_ap_rna(out, expression)
    protein_score, protein_z = score_ap_protein(project, out)
    robustness = score_robustness(
        out, rna_score, rna_z, protein_score, protein_z
    )
    write_estimate_input(out, expression)

    session = {
        "timestamp": datetime.now().astimezone().isoformat(),
        "os": platform.platform(),
        "python": sys.version,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "frozen_phase1B_sha256": frozen_hash,
        "thresholds": {
            "min_matched_n": MIN_MATCHED_N,
            "min_regression_n": MIN_REGRESSION_N,
            "min_rna_finite_fraction": MIN_RNA_FINITE_FRACTION,
            "min_rna_detected_fraction": MIN_RNA_DETECTED_FRACTION,
            "bootstrap_reps": BOOTSTRAP_REPS,
            "bootstrap_negative_fraction": BOOTSTRAP_NEGATIVE_FRACTION,
            "seed": SEED,
        },
        "coverage_script_discrepancies": int(
            coverage["script_mapping_discrepancy"].sum()
        ),
        "tier_counts": tiers["genomewide_evidence_tier"].value_counts().to_dict(),
        "association_rows": len(all_results),
        "rna_score_samples": len(rna_score),
        "protein_score_samples": len(protein_score),
        "robustness_rows": len(robustness),
    }
    (out / "Phase1C_genomewide_session.json").write_text(
        json.dumps(session, indent=2), encoding="utf-8"
    )
    print(json.dumps(session, indent=2))


if __name__ == "__main__":
    main()
