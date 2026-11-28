#!/usr/bin/env python3
"""Phase 2A TCGA purity/composition adjustment under frozen rules."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 20260731
ATTENUATION_LIMIT = 50.0
MIN_N = 80


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def bh(p):
    x = np.asarray(p, dtype=float)
    out = np.full(x.shape, np.nan)
    ok = np.isfinite(x)
    if not ok.any():
        return out
    v = x[ok]
    order = np.argsort(v)
    ranked = v[order] * len(v) / np.arange(1, len(v) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty_like(v)
    tmp[order] = np.minimum(ranked, 1.0)
    out[ok] = tmp
    return out


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def patient_id(sample: str):
    s = str(sample)
    if s.startswith("TCGA-"):
        return s[:12]
    return s


def zscore(s: pd.Series):
    x = pd.to_numeric(s, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else x * np.nan


def fit_ols(y: pd.Series, predictors: dict[str, pd.Series]):
    frame = pd.concat([y.rename("rna")] + [v.rename(k) for k, v in predictors.items()], axis=1)
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    result = {
        "n": len(frame), "status": "insufficient_data", "methylation_beta": np.nan,
        "methylation_se": np.nan, "methylation_p": np.nan,
        "methylation_standardized_beta": np.nan, "r2": np.nan,
    }
    if len(frame) < MIN_N or frame["rna"].nunique() < 3 or frame["methylation"].nunique() < 3:
        return result
    names = list(predictors)
    x = np.column_stack([np.ones(len(frame))] + [frame[n].to_numpy(float) for n in names])
    if np.linalg.matrix_rank(x) < x.shape[1]:
        result["status"] = "rank_deficient"
        return result
    yv = frame["rna"].to_numpy(float)
    beta, _, _, _ = np.linalg.lstsq(x, yv, rcond=None)
    residual = yv - x @ beta
    df = len(yv) - x.shape[1]
    if df <= 0:
        return result
    sigma2 = float(residual @ residual / df)
    covariance = sigma2 * np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    pvals = 2 * stats.t.sf(np.abs(beta / se), df)
    total = float(((yv - yv.mean()) ** 2).sum())
    mi = 1 + names.index("methylation")
    ysd = np.std(yv, ddof=0)
    msd = np.std(frame["methylation"].to_numpy(float), ddof=0)
    result.update({
        "status": "complete",
        "methylation_beta": float(beta[mi]),
        "methylation_se": float(se[mi]),
        "methylation_p": float(pvals[mi]),
        "methylation_standardized_beta": float(beta[mi] * msd / ysd) if ysd > 0 else np.nan,
        "r2": 1 - float(residual @ residual) / total if total > 0 else np.nan,
    })
    for i, name in enumerate(names, start=1):
        result[f"{name}_beta"] = float(beta[i])
        result[f"{name}_p"] = float(pvals[i])
    return result


def vif_table(frame: pd.DataFrame):
    clean = frame.dropna()
    rows = []
    for col in clean:
        y = clean[col].to_numpy(float)
        others = clean.drop(columns=col)
        x = np.column_stack([np.ones(len(clean)), others.to_numpy(float)])
        if np.linalg.matrix_rank(x) < x.shape[1]:
            vif = np.inf
        else:
            beta = np.linalg.lstsq(x, y, rcond=None)[0]
            residual = y - x @ beta
            total = ((y - y.mean()) ** 2).sum()
            r2 = 1 - (residual @ residual) / total if total > 0 else 1.0
            vif = 1 / max(1 - r2, np.finfo(float).eps)
        rows.append({"variable": col, "n_complete": len(clean), "VIF": vif})
    return pd.DataFrame(rows)


def build_covariates(project: Path, out: Path):
    scores = pd.read_csv(
        project / "phase1C_output/intermediate/TCGA_immune_feasibility_scores.tsv.gz",
        sep="\t",
    ).set_index("sample_id")
    selected = {
        "ESTIMATE_TumorPurity": "ESTIMATE_TumorPurity",
        "ESTIMATE_ImmuneScore": "ESTIMATE_ImmuneScore",
        "ESTIMATE_StromalScore": "ESTIMATE_StromalScore",
        "CD8_T_cells": "MCPcounter_CD8_T_cells",
        "NK_cells": "MCPcounter_NK_cells",
        "Monocytic_lineage": "MCPcounter_Monocytic_lineage",
        "B_lineage": "MCPcounter_B_lineage",
        "Endothelial_cells": "MCPcounter_Endothelial_cells",
        "Fibroblasts_CAF_proxy": "MCPcounter_Fibroblasts",
        "cytolytic_activity": "cytolytic_activity",
        "IFN_gamma_score": "IFN_gamma_response_score",
    }
    cov = scores[[v for v in selected.values()]].rename(columns={v: k for k, v in selected.items()})
    cov.index.name = "sample_id"
    absolute_path = out / "resources/TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
    absolute = pd.read_csv(absolute_path, sep="\t")
    absolute["array_key"] = absolute["array"].astype(str).str[:15]
    absolute = absolute.sort_values(["array_key", "call status"]).drop_duplicates("array_key")
    amap = absolute.set_index("array_key")["purity"]
    cov["ABSOLUTE_purity"] = [pd.to_numeric(amap.get(str(i)[:15]), errors="coerce") for i in cov.index]
    cov["patient_id"] = [patient_id(i) for i in cov.index]

    composition = [
        "ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "CD8_T_cells", "NK_cells",
        "Monocytic_lineage", "B_lineage", "Endothelial_cells",
        "Fibroblasts_CAF_proxy", "cytolytic_activity", "IFN_gamma_score",
    ]
    raw = cov[["ESTIMATE_TumorPurity"] + composition]
    corr = raw.corr(method="spearman")
    corr.index.name = "variable"
    corr.to_csv(out / "composition_correlation_matrix.tsv", sep="\t")
    vif = vif_table(raw.apply(zscore))
    vif.to_csv(out / "composition_vif.tsv", sep="\t", index=False)

    purity_z = zscore(cov["ESTIMATE_TumorPurity"])
    residualized = {}
    for name in composition:
        z = zscore(cov[name])
        f = pd.concat([z, purity_z], axis=1).dropna()
        x = np.column_stack([np.ones(len(f)), f.iloc[:, 1].to_numpy(float)])
        b = np.linalg.lstsq(x, f.iloc[:, 0].to_numpy(float), rcond=None)[0]
        r = pd.Series(np.nan, index=cov.index)
        r.loc[f.index] = f.iloc[:, 0] - x @ b
        residualized[name] = r
    residual = pd.DataFrame(residualized).dropna()
    standardized = residual.apply(zscore)
    u, s, vt = np.linalg.svd(standardized.to_numpy(float), full_matrices=False)
    variance = s**2 / np.sum(s**2)
    n_pc = min(5, int(np.searchsorted(np.cumsum(variance), 0.80) + 1))
    pcs = pd.DataFrame(
        u[:, :n_pc] * s[:n_pc],
        index=standardized.index,
        columns=[f"composition_PC{i+1}" for i in range(n_pc)],
    )
    for col in pcs:
        cov[col] = pcs[col]
    q1 = float(cov["ESTIMATE_TumorPurity"].quantile(1 / 3))
    q2 = float(cov["ESTIMATE_TumorPurity"].quantile(2 / 3))
    cov["purity_stratum"] = np.where(
        cov["ESTIMATE_TumorPurity"] >= q2, "high",
        np.where(cov["ESTIMATE_TumorPurity"] >= q1, "medium", "low"),
    )
    cov["purity_stratum_rule"] = f"tertiles fixed before model inspection: low < {q1:.9g}; medium [{q1:.9g}, {q2:.9g}); high >= {q2:.9g}"
    cov["composition_PC_rule"] = (
        f"PCA of 10 standardized composition covariates residualized on ESTIMATE purity; "
        f"minimum PCs reaching >=80% variance, capped at 5; retained={n_pc}; "
        f"cumulative_variance={np.cumsum(variance)[n_pc-1]:.6f}"
    )
    cov.reset_index().to_csv(out / "purity_and_composition_covariates.tsv", sep="\t", index=False)
    pca = pd.DataFrame({
        "component": [f"PC{i+1}" for i in range(len(variance))],
        "variance_fraction": variance,
        "cumulative_variance": np.cumsum(variance),
        "retained": np.arange(len(variance)) < n_pc,
    })
    pca.to_csv(out / "composition_pca_variance.tsv", sep="\t", index=False)
    return cov, [f"composition_PC{i+1}" for i in range(n_pc)], q1, q2


def fit_all_models(project: Path, out: Path, cov: pd.DataFrame, pc_cols: list[str]):
    phase1c = project / "phase1C_output"
    tiers = pd.read_csv(phase1c / "genomewide_methylation_candidate_tiers.tsv", sep="\t")
    genes = set(tiers["gene"].astype(str))
    gw = load_module("phase1c_genomewide", project / "scripts/python/phase1C/phase1C_genomewide.py")
    expression = gw.load_expression(project, genes)
    cnv = gw.load_cnv(project, genes)
    beta = pd.read_csv(
        phase1c / "gene_promoter_beta_genomewide_TSS200.tsv.gz", sep="\t", index_col=0
    )
    beta.index = beta.index.astype(str)

    sample_sets = [set(expression.columns), set(beta.columns), set(cov.index)]
    primary = sorted(set.intersection(*sample_sets))
    cnv_samples = sorted(set(primary) & set(cnv.columns))
    high = [s for s in cnv_samples if cov.loc[s, "purity_stratum"] == "high"]
    medium = [s for s in cnv_samples if cov.loc[s, "purity_stratum"] == "medium"]
    absolute = [s for s in cnv_samples if np.isfinite(cov.loc[s, "ABSOLUTE_purity"])]

    rows = []
    for idx, row in tiers.iterrows():
        gene = str(row["gene"])
        if gene not in expression.index or gene not in beta.index:
            rows.append({"gene": gene, "overall_status": "missing_expression_or_methylation"})
            continue
        y = expression.loc[gene]
        m = beta.loc[gene]
        c = cnv.loc[gene] if gene in cnv.index else pd.Series(dtype=float)
        fits = {
            "unadjusted": fit_ols(y.loc[primary], {"methylation": m.loc[primary]}),
            "cnv": fit_ols(y.loc[cnv_samples], {"methylation": m.loc[cnv_samples], "cnv": c.loc[cnv_samples]}),
            "purity": fit_ols(y.loc[cnv_samples], {
                "methylation": m.loc[cnv_samples], "cnv": c.loc[cnv_samples],
                "ESTIMATE_TumorPurity": cov.loc[cnv_samples, "ESTIMATE_TumorPurity"],
            }),
            "composition": fit_ols(y.loc[cnv_samples], {
                "methylation": m.loc[cnv_samples], "cnv": c.loc[cnv_samples],
                "ESTIMATE_TumorPurity": cov.loc[cnv_samples, "ESTIMATE_TumorPurity"],
                **{pc: cov.loc[cnv_samples, pc] for pc in pc_cols},
            }),
            "high_purity": fit_ols(y.loc[high], {"methylation": m.loc[high], "cnv": c.loc[high]}),
            "medium_purity": fit_ols(y.loc[medium], {"methylation": m.loc[medium], "cnv": c.loc[medium]}),
            "absolute_purity": fit_ols(y.loc[absolute], {
                "methylation": m.loc[absolute], "cnv": c.loc[absolute],
                "ABSOLUTE_purity": cov.loc[absolute, "ABSOLUTE_purity"],
            }),
        }
        rec = {"gene": gene, "overall_status": "complete"}
        for model, fit in fits.items():
            for key, value in fit.items():
                rec[f"{model}_{key}"] = value
        rows.append(rec)
        if (idx + 1) % 1000 == 0:
            print(f"MODEL_PROGRESS {idx+1}/{len(tiers)}", flush=True)
    result = pd.DataFrame(rows)
    for model in ["unadjusted", "cnv", "purity", "composition", "high_purity", "medium_purity", "absolute_purity"]:
        result[f"{model}_methylation_fdr"] = bh(result[f"{model}_methylation_p"])
    result["purity_vs_cnv_attenuation_percent"] = np.where(
        result["cnv_methylation_beta"].abs() > 0,
        (1 - result["purity_methylation_beta"].abs() / result["cnv_methylation_beta"].abs()) * 100,
        np.nan,
    )
    result["composition_vs_cnv_attenuation_percent"] = np.where(
        result["cnv_methylation_beta"].abs() > 0,
        (1 - result["composition_methylation_beta"].abs() / result["cnv_methylation_beta"].abs()) * 100,
        np.nan,
    )
    result["purity_direction_preserved_from_cnv"] = (
        np.sign(result["purity_methylation_beta"]) == np.sign(result["cnv_methylation_beta"])
    )
    result["composition_direction_preserved_from_cnv"] = (
        np.sign(result["composition_methylation_beta"]) == np.sign(result["cnv_methylation_beta"])
    )
    result["high_purity_direction_negative"] = result["high_purity_methylation_beta"] < 0
    result["potential_cell_composition_confounding"] = (
        (result["purity_methylation_fdr"] >= .05)
        | (result["purity_vs_cnv_attenuation_percent"] > ATTENUATION_LIMIT)
        | (~result["purity_direction_preserved_from_cnv"])
        | (result["composition_methylation_fdr"] >= .05)
        | (~result["composition_direction_preserved_from_cnv"])
    )
    result["attenuation_threshold_percent_frozen"] = ATTENUATION_LIMIT
    result["model_missingness_rule"] = "complete cases per model; no imputation"
    result["primary_purity_metric"] = "ESTIMATE_TumorPurity"
    result["absolute_purity_role"] = "secondary sensitivity; PanCanAtlas ABSOLUTE"
    merged = tiers.merge(result, on="gene", how="left", validate="one_to_one")
    merged.to_csv(out / "genomewide_purity_adjusted_models.tsv.gz", sep="\t", index=False, compression="gzip")

    comp_cols = ["gene"] + [c for c in merged.columns if c.startswith("composition_")] + [
        "genomewide_evidence_tier", "potential_cell_composition_confounding",
        "purity_vs_cnv_attenuation_percent", "high_purity_direction_negative",
    ]
    merged[comp_cols].to_csv(out / "immune_composition_adjusted_models.tsv.gz", sep="\t", index=False, compression="gzip")
    hp_cols = ["gene", "genomewide_evidence_tier"] + [c for c in merged.columns if c.startswith("high_purity_") or c.startswith("medium_purity_")]
    merged[hp_cols].to_csv(out / "high_purity_sensitivity_results.tsv.gz", sep="\t", index=False, compression="gzip")

    prior = merged["genomewide_evidence_tier"].isin(["Stable", "Partially stable"])
    technical = prior & (
        merged["purity_status"].ne("complete")
        | merged["composition_status"].ne("complete")
        | merged["high_purity_status"].ne("complete")
    )
    robust = (
        merged["genomewide_evidence_tier"].eq("Stable")
        & merged["purity_methylation_beta"].lt(0)
        & merged["purity_methylation_fdr"].lt(.05)
        & merged["purity_vs_cnv_attenuation_percent"].le(ATTENUATION_LIMIT)
        & merged["high_purity_methylation_beta"].lt(0)
        & merged["bootstrap_negative_fraction"].ge(.70)
        & merged["sensitivity_any_negative_both"].eq(True)
    )
    sensitive = prior & ~technical & ~robust & merged["potential_cell_composition_confounding"].eq(True)
    merged["composition_evidence_class"] = np.select(
        [robust, technical, sensitive],
        ["Composition-robust", "Indeterminate", "Composition-sensitive"],
        default="Not fully robust",
    )
    candidate_cols = [
        "gene", "genomewide_evidence_tier", "composition_evidence_class",
        "methylation_rna_rho", "correlation_fdr", "adjusted_methylation_coefficient",
        "adjusted_methylation_fdr", "bootstrap_negative_fraction",
        "purity_methylation_beta", "purity_methylation_standardized_beta",
        "purity_methylation_p", "purity_methylation_fdr",
        "composition_methylation_beta", "composition_methylation_p",
        "composition_methylation_fdr", "high_purity_methylation_beta",
        "high_purity_methylation_p", "high_purity_methylation_fdr",
        "purity_vs_cnv_attenuation_percent", "composition_vs_cnv_attenuation_percent",
        "potential_cell_composition_confounding", "sensitivity_any_negative_both",
    ]
    merged.loc[prior, candidate_cols].to_csv(
        out / "composition_robust_methylation_candidates.tsv", sep="\t", index=False
    )
    context = {
        "all_gene_n": int(len(merged)),
        "phase1C_stable_n": int(merged["genomewide_evidence_tier"].eq("Stable").sum()),
        "phase1C_partial_n": int(merged["genomewide_evidence_tier"].eq("Partially stable").sum()),
        "composition_robust_stable_n": int((robust & merged["genomewide_evidence_tier"].eq("Stable")).sum()),
        "composition_robust_total_n": int(robust.sum()),
        "composition_sensitive_n": int(sensitive.sum()),
        "indeterminate_n": int(technical.sum()),
        "high_purity_negative_stable_n": int(
            (merged["genomewide_evidence_tier"].eq("Stable") & merged["high_purity_methylation_beta"].lt(0)).sum()
        ),
        "primary_common_sample_n": len(primary),
        "cnv_model_sample_n": len(cnv_samples),
        "high_purity_model_sample_n": len(high),
        "medium_purity_model_sample_n": len(medium),
        "absolute_model_sample_n": len(absolute),
    }
    return merged, context, expression, beta, cov


def expected_direction(frame):
    return (
        ((frame["ranking_metric"] != "rank_hypermethylated_rna_down") & (frame["enrichment_direction"] == "low"))
        | ((frame["ranking_metric"] == "rank_hypermethylated_rna_down") & (frame["enrichment_direction"] == "high"))
    )


def pathway_after_purity(project: Path, out: Path, merged: pd.DataFrame):
    integ = load_module("phase1c_integrate", project / "scripts/python/phase1C/phase1C_integrate.py")
    original = pd.read_csv(project / "phase1C_output/genomewide_ranked_pathway_enrichment.tsv", sep="\t")
    focus = original[
        original["promoter_definition"].eq("TSS200")
        & original["focus_category"].fillna("").str.contains(
            "antigen_processing_presentation|interferon_gamma|inflammatory_response|cytotoxic_lymphocyte"
        )
    ].copy()
    original13 = focus[expected_direction(focus) & focus["fdr"].lt(.05)].copy()

    work = merged.drop_duplicates("gene").set_index("gene")
    tiny = np.finfo(float).tiny
    work["rank_directed_neglog10_fdr"] = (
        -np.log10(work["purity_methylation_fdr"].clip(lower=tiny))
        * np.sign(work["purity_methylation_beta"])
    )
    work["rank_cnv_adjusted_beta"] = work["purity_methylation_beta"]
    work["rank_bootstrap_weighted_effect"] = (
        work["purity_methylation_standardized_beta"]
        * (2 * work["bootstrap_negative_fraction"] - 1)
    )
    md = work["tumor_normal_methylation_delta"]
    pb = work["purity_methylation_standardized_beta"]
    work["rank_hypermethylated_rna_down"] = zscore(md) - zscore(pb)

    rows = []
    for metric in [
        "rank_directed_neglog10_fdr", "rank_cnv_adjusted_beta",
        "rank_bootstrap_weighted_effect", "rank_hypermethylated_rna_down",
    ]:
        values = work[metric].dropna()
        rank_array = stats.rankdata(values.to_numpy(), method="average")
        universe = set(values.index)
        lookup = dict(zip(values.index, rank_array))
        total = len(values)
        _, tie_counts = np.unique(values.to_numpy(), return_counts=True)
        tie_term = float(np.sum(tie_counts**3 - tie_counts))
        for source, name, set_id, members in integ.load_gene_sets(project):
            test = integ.competitive_rank_test(lookup, universe, total, members, tie_term)
            if test is None:
                continue
            n1, n2, effect, p, overlap = test
            rows.append({
                "source": source, "pathway": name, "pathway_id": set_id,
                "ranking_metric": metric, "set_size_tested": n1,
                "background_size": n2, "purity_adjusted_rank_effect": effect,
                "enrichment_direction": "high" if effect > 0 else "low",
                "p_value": p, "member_genes_tested": ";".join(sorted(overlap)),
            })
    all_enrich = pd.DataFrame(rows)
    all_enrich["fdr"] = np.nan
    for _, idx in all_enrich.groupby(["source", "ranking_metric"]).groups.items():
        all_enrich.loc[idx, "fdr"] = bh(all_enrich.loc[idx, "p_value"])

    keys = ["source", "pathway", "ranking_metric"]
    selected = original13.merge(all_enrich, on=keys, how="left", suffixes=("_phase1C", "_purity"))
    selected["expected_direction_after_purity"] = (
        ((selected["ranking_metric"] != "rank_hypermethylated_rna_down") & selected["enrichment_direction_purity"].eq("low"))
        | ((selected["ranking_metric"] == "rank_hypermethylated_rna_down") & selected["enrichment_direction_purity"].eq("high"))
    )
    selected["retained_after_purity"] = selected["expected_direction_after_purity"] & selected["fdr_purity"].lt(.05)
    selected["multiple_testing_scope_after_purity"] = "BH across all tested pathways within source x ranking metric"
    selected.to_csv(out / "pathway_enrichment_after_purity_adjustment.tsv", sep="\t", index=False)
    all_enrich.to_csv(out / "intermediate/pathway_enrichment_after_purity_all_tests.tsv.gz", sep="\t", index=False, compression="gzip")
    return selected, original13


def report(out: Path, context: dict, pathways: pd.DataFrame, cov: pd.DataFrame):
    retained = int(pathways["retained_after_purity"].sum())
    abs_n = int(cov["ABSOLUTE_purity"].notna().sum())
    text = f"""# Cell-composition confounding report

## Frozen model rules

- Main purity covariate: ESTIMATE tumour purity.
- Secondary purity sensitivity: PanCanAtlas ABSOLUTE purity; available for {abs_n} TCGA primary tumour samples in the analysis covariate table.
- Main model: RNA ~ TSS200 promoter methylation + gene-level CNV + ESTIMATE purity.
- Composition sensitivity: the main model plus prespecified residualized composition principal components.
- High-purity threshold: upper ESTIMATE-purity tertile, frozen before gene-model inspection.
- Effect attenuation threshold: >{ATTENUATION_LIMIT:.0f}% relative to the CNV-adjusted methylation coefficient, frozen before result inspection.
- Missing data: complete cases per model; no imputation.

## Results

- Phase 1C Stable genes: {context['phase1C_stable_n']:,}.
- Stable genes meeting all Composition-robust criteria: {context['composition_robust_stable_n']:,}.
- Stable/Partially stable genes flagged Composition-sensitive: {context['composition_sensitive_n']:,}.
- Technically indeterminate genes: {context['indeterminate_n']:,}.
- Prespecified immune pathway tests retained after purity adjustment: {retained}/{len(pathways)}.
- Model sample sizes: primary overlap {context['primary_common_sample_n']}; CNV/purity {context['cnv_model_sample_n']}; high-purity {context['high_purity_model_sample_n']}; medium-purity {context['medium_purity_model_sample_n']}.

These are bulk-tissue associations. A retained coefficient does not establish a tumour-cell-intrinsic or causal methylation effect.
"""
    (out / "cell_composition_confounding_report.md").write_text(text, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()
    project = Path(args.project_root)
    out = project / "phase2A_output"
    (out / "intermediate").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    cov, pc_cols, q1, q2 = build_covariates(project, out)
    merged, context, _, _, cov = fit_all_models(project, out, cov, pc_cols)
    pathways, original13 = pathway_after_purity(project, out, merged)
    context.update({
        "pathway_retained_n": int(pathways["retained_after_purity"].sum()),
        "pathway_original_n": int(len(pathways)),
        "purity_q1": q1, "purity_q2": q2,
        "composition_pc_n": len(pc_cols),
        "script_sha256": sha256(Path(__file__)),
        "seed": SEED,
    })
    (out / "intermediate/tcga_model_context.json").write_text(
        json.dumps(context, indent=2), encoding="utf-8"
    )
    report(out, context, pathways, cov)
    print(json.dumps(context, indent=2), flush=True)


if __name__ == "__main__":
    main()
