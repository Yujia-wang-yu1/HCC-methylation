#!/usr/bin/env python3
"""Lock the Phase 2A immune methylation program before GEO outcome inspection."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 20260731
BOOTSTRAP_REPS = 200


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    assert spec.loader
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def zrows(frame):
    mean = frame.mean(axis=1)
    sd = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(mean, axis=0).div(sd, axis=0)


def expected_direction(frame):
    return (
        ((frame["ranking_metric"] != "rank_hypermethylated_rna_down")
         & frame["enrichment_direction"].eq("low"))
        | ((frame["ranking_metric"] == "rank_hypermethylated_rna_down")
           & frame["enrichment_direction"].eq("high"))
    )


def parse_gpl(path: Path):
    header = None
    mapping = {}
    ended = False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("ID\t") and "Gene Symbol" in line:
                header = line.split("\t")
                continue
            if line == "!platform_table_end":
                ended = True
            if header and line and line[0].isdigit():
                fields = line.split("\t")
                if len(fields) <= max(header.index("ID"), header.index("Gene Symbol")):
                    continue
                pid = fields[header.index("ID")]
                for gene in fields[header.index("Gene Symbol")].split(" /// "):
                    gene = gene.strip()
                    if gene and gene != "---":
                        mapping.setdefault(gene, set()).add(pid)
    if not ended:
        raise RuntimeError("GPL13667 platform table lacks !platform_table_end")
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()
    project = Path(args.project_root)
    out = project / "phase2A_output"
    inter = out / "intermediate"
    inter.mkdir(parents=True, exist_ok=True)

    models = pd.read_csv(out / "composition_robust_methylation_candidates.tsv", sep="\t")
    original = pd.read_csv(project / "phase1C_output/genomewide_ranked_pathway_enrichment.tsv", sep="\t")
    focus = original[
        original["promoter_definition"].eq("TSS200")
        & original["focus_category"].fillna("").str.contains(
            "antigen_processing_presentation|interferon_gamma|inflammatory_response|cytotoxic_lymphocyte"
        )
    ].copy()
    original13 = focus[expected_direction(focus) & focus["fdr"].lt(.05)].copy()
    keys = set(zip(original13["source"], original13["pathway"], original13["ranking_metric"]))
    leading = pd.read_csv(project / "phase1C_output/leading_edge_gene_audit.tsv", sep="\t")
    leading = leading[
        leading.apply(lambda r: (r["source"], r["pathway"], r["ranking_metric"]) in keys, axis=1)
    ].copy()
    provenance = {}
    for _, row in leading.iterrows():
        for gene in str(row["leading_edge_genes"]).split(";"):
            provenance.setdefault(gene, []).append(
                f"{row['source']}|{row['pathway']}|{row['ranking_metric']}"
            )

    gpl = parse_gpl(out / "external/GEO/raw_downloads/GPL13667_full.soft.txt")
    icgc_gene_col = pd.read_csv(
        project / "data_processed/ICGC_LIRI_JP/ICGC_LIRI_JP_expression_matrix.tsv.gz",
        sep="\t", usecols=[0],
    ).iloc[:, 0].astype(str)
    icgc_genes = set(icgc_gene_col)
    candidate = models[
        models["gene"].isin(provenance)
        & models["composition_evidence_class"].eq("Composition-robust")
        & models["sensitivity_any_negative_both"].eq(True)
        & models["gene"].ne("RAB27A")
    ].copy()
    candidate["ICGC_expression_mapped"] = candidate["gene"].isin(icgc_genes)
    candidate["GEO_U219_probe_n"] = candidate["gene"].map(lambda g: len(gpl.get(g, set())))
    candidate["GEO_U219_probe_ids"] = candidate["gene"].map(
        lambda g: ";".join(sorted(gpl.get(g, set())))
    )
    candidate = candidate[
        candidate["ICGC_expression_mapped"] & candidate["GEO_U219_probe_n"].gt(0)
    ].copy()
    candidate["leading_edge_provenance"] = candidate["gene"].map(
        lambda g: ";".join(sorted(set(provenance[g])))
    )
    if candidate.empty:
        raise RuntimeError("Frozen program contains no genes after prespecified mapping filters")
    raw_weights = candidate["purity_methylation_standardized_beta"].abs()
    candidate["MISP_weight"] = raw_weights / raw_weights.sum()
    candidate["weight_rule"] = (
        "absolute standardized methylation coefficient from RNA~TSS200 methylation+CNV+ESTIMATE purity; "
        "normalized to sum 1; no outcome fitting"
    )
    candidate["locked_before_GEO_statistical_validation"] = True
    candidate["lock_seed"] = SEED
    candidate["program_name_status"] = (
        "Methylation-associated Immune Suppression Program (provisional until external validation)"
    )

    columns = [
        "gene", "MISP_weight", "leading_edge_provenance", "genomewide_evidence_tier",
        "composition_evidence_class", "methylation_rna_rho", "correlation_fdr",
        "adjusted_methylation_coefficient", "adjusted_methylation_fdr",
        "purity_methylation_beta", "purity_methylation_standardized_beta",
        "purity_methylation_fdr", "purity_vs_cnv_attenuation_percent",
        "high_purity_methylation_beta", "bootstrap_negative_fraction",
        "sensitivity_any_negative_both", "ICGC_expression_mapped",
        "GEO_U219_probe_n", "GEO_U219_probe_ids", "weight_rule",
        "locked_before_GEO_statistical_validation", "lock_seed", "program_name_status",
    ]
    candidate[columns].sort_values("MISP_weight", ascending=False).to_csv(
        out / "locked_immune_methylation_program_genes.tsv", sep="\t", index=False
    )

    gw = module("phase1c_genomewide_lock", project / "scripts/python/phase1C/phase1C_genomewide.py")
    genes = set(candidate["gene"])
    expr = gw.load_expression(project, genes)
    cnv = gw.load_cnv(project, genes)
    beta200 = pd.read_csv(
        project / "phase1C_output/gene_promoter_beta_genomewide_TSS200.tsv.gz",
        sep="\t", index_col=0,
    )
    cov = pd.read_csv(out / "purity_and_composition_covariates.tsv", sep="\t", index_col="sample_id")
    samples = sorted(set(expr.columns) & set(beta200.columns) & set(cov.index))
    weights = candidate.set_index("gene")["MISP_weight"]
    ordered = [g for g in weights.index if g in expr.index and g in beta200.index]
    weights = weights.loc[ordered]
    methyl_z = zrows(beta200.loc[ordered, samples])
    rna_z = zrows(expr.loc[ordered, samples])
    methyl_score = methyl_z.mul(weights, axis=0).sum(axis=0) / weights.sum()
    rna_suppression = -rna_z.mul(weights, axis=0).sum(axis=0) / weights.sum()
    scores = pd.DataFrame({
        "sample_id": samples,
        "patient_id": [s[:12] for s in samples],
        "MISP_methylation_score": methyl_score.reindex(samples).to_numpy(),
        "MISP_RNA_suppression_score": rna_suppression.reindex(samples).to_numpy(),
        "program_gene_n": len(ordered),
        "program_genes": ";".join(ordered),
        "score_direction": "higher = higher methylation-associated immune-program suppression",
        "weights_frozen_without_outcome": True,
    })
    scores.to_csv(out / "MISP_score_TCGA.tsv", sep="\t", index=False)
    pd.DataFrame({
        "gene": ordered,
        "weight": weights.loc[ordered].to_numpy(),
        "methylation_definition": "TSS200 gene-level median of QC-passing HM450 probes",
        "methylation_standardization": "gene-wise z score across 371 TCGA unique primary tumours",
        "RNA_definition": "negative weighted mean of gene-wise z log2(TPM+1)",
        "RNA_standardization": "gene-wise z score across 371 TCGA unique primary tumours",
        "weight_source": "absolute standardized ESTIMATE-purity-adjusted methylation coefficient",
        "outcome_used": False,
    }).to_csv(out / "MISP_score_definition.tsv", sep="\t", index=False)

    robustness = []
    for gene in ordered:
        remain = [g for g in ordered if g != gene]
        w = weights.loc[remain]
        w = w / w.sum()
        ms = methyl_z.loc[remain].mul(w, axis=0).sum(axis=0)
        rs = -rna_z.loc[remain].mul(w, axis=0).sum(axis=0)
        robustness.append({
            "analysis": "leave_one_gene_out", "omitted": gene, "remaining_gene_n": len(remain),
            "methylation_score_spearman_vs_full": stats.spearmanr(ms, methyl_score).statistic,
            "RNA_score_spearman_vs_full": stats.spearmanr(rs, rna_suppression).statistic,
        })
    pathway_to_genes = {}
    for gene in ordered:
        for item in provenance[gene]:
            source, pathway, metric = item.split("|", 2)
            pathway_to_genes.setdefault(f"{source}|{pathway}", set()).add(gene)
    for pathway, omitted in sorted(pathway_to_genes.items()):
        remain = [g for g in ordered if g not in omitted]
        if len(remain) < 2:
            continue
        w = weights.loc[remain] / weights.loc[remain].sum()
        ms = methyl_z.loc[remain].mul(w, axis=0).sum(axis=0)
        rs = -rna_z.loc[remain].mul(w, axis=0).sum(axis=0)
        robustness.append({
            "analysis": "leave_one_pathway_out", "omitted": pathway,
            "omitted_gene_n": len(omitted), "remaining_gene_n": len(remain),
            "methylation_score_spearman_vs_full": stats.spearmanr(ms, methyl_score).statistic,
            "RNA_score_spearman_vs_full": stats.spearmanr(rs, rna_suppression).statistic,
        })

    rng = np.random.default_rng(SEED)
    fitmod = module("phase2a_tcga_fit", project / "scripts/phase2A/phase2A_tcga_models.py")
    model_samples = sorted(set(samples) & set(cnv.columns))
    for rep in range(BOOTSTRAP_REPS):
        draw = rng.choice(model_samples, len(model_samples), replace=True)
        boot_weights = {}
        for gene in ordered:
            y = pd.Series(expr.loc[gene, draw].to_numpy(), index=np.arange(len(draw)))
            m = pd.Series(beta200.loc[gene, draw].to_numpy(), index=np.arange(len(draw)))
            c = pd.Series(cnv.loc[gene, draw].to_numpy(), index=np.arange(len(draw)))
            p = pd.Series(cov.loc[draw, "ESTIMATE_TumorPurity"].to_numpy(), index=np.arange(len(draw)))
            fit = fitmod.fit_ols(y, {"methylation": m, "cnv": c, "purity": p})
            boot_weights[gene] = abs(fit["methylation_standardized_beta"])
        bw = pd.Series(boot_weights).replace([np.inf, -np.inf], np.nan).dropna()
        if bw.sum() <= 0 or len(bw) < 2:
            continue
        bw /= bw.sum()
        ms = methyl_z.loc[bw.index].mul(bw, axis=0).sum(axis=0)
        rs = -rna_z.loc[bw.index].mul(bw, axis=0).sum(axis=0)
        robustness.append({
            "analysis": "bootstrap_weight_refit", "omitted": f"rep_{rep+1}",
            "remaining_gene_n": len(bw),
            "methylation_score_spearman_vs_full": stats.spearmanr(ms, methyl_score).statistic,
            "RNA_score_spearman_vs_full": stats.spearmanr(rs, rna_suppression).statistic,
        })
    for definition, filename in [
        ("TSS1500", "gene_promoter_beta_genomewide_TSS1500.tsv.gz"),
        ("promoter_union", "gene_promoter_beta_genomewide_promoter_union.tsv.gz"),
    ]:
        alt = pd.read_csv(project / "phase1C_output" / filename, sep="\t", index_col=0)
        present = [g for g in ordered if g in alt.index]
        w = weights.loc[present] / weights.loc[present].sum()
        alt_z = zrows(alt.loc[present, samples])
        alt_score = alt_z.mul(w, axis=0).sum(axis=0)
        robustness.append({
            "analysis": "promoter_definition_sensitivity", "omitted": definition,
            "remaining_gene_n": len(present),
            "methylation_score_spearman_vs_full": stats.spearmanr(alt_score, methyl_score).statistic,
            "RNA_score_spearman_vs_full": np.nan,
        })
    pd.DataFrame(robustness).to_csv(out / "MISP_score_robustness.tsv", sep="\t", index=False)

    program_path = out / "locked_immune_methylation_program_genes.tsv"
    definition_path = out / "MISP_score_definition.tsv"
    lock = {
        "locked_at_stage": "before GSE63898 expression matrix inspection and before GEO statistical validation",
        "gene_n": len(ordered), "genes": ordered,
        "program_table_sha256": sha(program_path),
        "score_definition_sha256": sha(definition_path),
        "seed": SEED, "bootstrap_reps": BOOTSTRAP_REPS,
    }
    (out / "locked_program_hashes.json").write_text(json.dumps(lock, indent=2), encoding="utf-8")
    (out / "locked_program_gene_provenance.md").write_text(
        "# Locked immune methylation program provenance\n\n"
        f"The program contains {len(ordered)} genes. Genes were locked before GEO outcome inspection by intersecting: "
        "Phase 1C prespecified immune-pathway leading edges, Phase 1C Stable status, Phase 2A "
        "Composition-robust status, promoter-sensitivity direction support, RAB27A exclusion, "
        "ICGC expression coverage, and GPL13667 exact gene-symbol probe coverage.\n\n"
        "Weights are absolute standardized ESTIMATE-purity-adjusted methylation coefficients normalized to sum one. "
        "No survival, clinical outcome, or GEO validation result was used.\n",
        encoding="utf-8",
    )
    print(json.dumps(lock, indent=2))


if __name__ == "__main__":
    main()
