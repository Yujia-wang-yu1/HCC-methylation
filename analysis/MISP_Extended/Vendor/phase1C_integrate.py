#!/usr/bin/env python3
"""Phase 1C integration: ranked pathways, immune feasibility, and external audit."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import platform
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


CORE_AP = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
           "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"]
CD8_GENES = ["CD8A", "CD8B", "CD3D", "CD3E", "LCK"]
NK_GENES = ["NCR1", "KLRD1", "NKG7", "GNLY", "PRF1", "XCL1", "XCL2"]
CYTOLYTIC_GENES = ["GZMA", "PRF1"]
CHECKPOINT_GENES = ["PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2"]
MIN_SET = 10
MAX_SET = 500


def bh(values):
    a = np.asarray(values, float)
    out = np.full(a.shape, np.nan)
    ok = np.isfinite(a)
    p = a[ok]
    if not len(p):
        return out
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    restore = np.empty_like(q)
    restore[order] = np.minimum(q, 1)
    out[ok] = restore
    return out


def zscore_rows(frame):
    return frame.sub(frame.mean(axis=1), axis=0).div(
        frame.std(axis=1, ddof=0).replace(0, np.nan), axis=0
    )


def sample_code(x):
    parts = str(x).split("-")
    return parts[3][:2] if len(parts) >= 4 else ""


def patient_id(x):
    return "-".join(str(x).split("-")[:3])


def unique_primary(samples):
    grouped = defaultdict(list)
    for s in samples:
        if sample_code(s) == "01":
            grouped[patient_id(s)].append(s)
    return [v[0] for _, v in sorted(grouped.items()) if len(v) == 1]


def parse_gmt(path, source):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                rows.append((source, fields[0], fields[1], set(fields[2:])))
    return rows


def parse_go(obo, gaf):
    names, namespace, parents = {}, {}, defaultdict(set)
    current = {}
    with open(obo, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]":
                if current.get("id"):
                    names[current["id"]] = current.get("name", current["id"])
                    namespace[current["id"]] = current.get("namespace", "")
                    parents[current["id"]] |= set(current.get("parents", []))
                current = {}
            elif ": " in line:
                key, value = line.split(": ", 1)
                if key in {"id", "name", "namespace"}:
                    current[key] = value
                elif key == "is_a":
                    current.setdefault("parents", []).append(value.split()[0])
                elif key == "relationship" and value.startswith("part_of "):
                    current.setdefault("parents", []).append(value.split()[1])
    if current.get("id"):
        names[current["id"]] = current.get("name", current["id"])
        namespace[current["id"]] = current.get("namespace", "")
        parents[current["id"]] |= set(current.get("parents", []))
    members = defaultdict(set)
    with gzip.open(gaf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or "NOT" in f[3].split("|"):
                continue
            go_id, symbol = f[4], f[2]
            if namespace.get(go_id) == "biological_process":
                members[go_id].add(symbol)
    ancestor_cache = {}
    def ancestors(go_id):
        if go_id in ancestor_cache:
            return ancestor_cache[go_id]
        seen, stack = set(), list(parents.get(go_id, []))
        while stack:
            parent = stack.pop()
            if parent in seen:
                continue
            seen.add(parent)
            stack.extend(parents.get(parent, []))
        ancestor_cache[go_id] = seen
        return seen
    propagated = defaultdict(set)
    for go_id, genes in members.items():
        propagated[go_id] |= genes
        for parent in ancestors(go_id):
            if namespace.get(parent) == "biological_process":
                propagated[parent] |= genes
    return [("GO_BP", f"{go}:{names.get(go, go)}", go, genes)
            for go, genes in propagated.items()]


def parse_kegg(resource):
    names = {}
    for line in (resource / "KEGG_hsa_pathway_list.tsv").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            names[f[0].replace("path:", "")] = re.sub(r" - Homo sapiens.*$", "", f[1])
    gene_symbols = {}
    for line in (resource / "KEGG_hsa_gene_list.tsv").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            symbol = f[1].split(";")[0].split(",")[0].strip()
            gene_symbols[f[0].replace("hsa:", "")] = symbol
    members = defaultdict(set)
    for line in (resource / "KEGG_hsa_pathway_gene_links.tsv").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            pathway = f[0].replace("path:", "")
            gid = f[1].replace("hsa:", "")
            if gid in gene_symbols:
                members[pathway].add(gene_symbols[gid])
    return [("KEGG", f"{pid}:{names.get(pid, pid)}", pid, genes)
            for pid, genes in members.items()]


def load_gene_sets(project):
    public = project / "resources" / "public_annotations"
    resource = public / "phase1C"
    sets = []
    sets += parse_go(public / "go-basic.obo", public / "goa_human.gaf.gz")
    sets += parse_gmt(resource / "reactome_gmt" / "ReactomePathways.gmt", "Reactome")
    sets += parse_kegg(resource)
    sets += parse_gmt(resource / "h.all.v2026.1.Hs.symbols.gmt", "Hallmark")
    return sets


def focus_category(name):
    n = name.lower().replace("_", " ")
    rules = [
        ("antigen_processing_presentation", ["antigen processing", "antigen presentation",
                                             "mhc class i", "mhc class ii"]),
        ("interferon_gamma", ["interferon gamma", "ifn gamma"]),
        ("inflammatory_response", ["inflammatory response", "inflammation"]),
        ("cytotoxic_lymphocyte", ["cytotoxic", "natural killer", "cell killing",
                                  "t cell mediated cytotoxic"]),
        ("vesicle_exocytosis", ["vesicle", "exocytosis", "degranulation"]),
        ("metabolism", ["metabolic process", "metabolism", "glycolysis",
                        "fatty acid", "oxidative phosphorylation", "amino acid"]),
    ]
    return ";".join(label for label, terms in rules if any(t in n for t in terms))


def rank_metrics(frame):
    x = frame.copy()
    tiny = np.finfo(float).tiny
    x["rank_directed_neglog10_fdr"] = (
        -np.log10(x["correlation_fdr"].clip(lower=tiny))
        * np.sign(x["methylation_rna_rho"])
    )
    x["rank_cnv_adjusted_beta"] = x["adjusted_methylation_coefficient"]
    x["rank_bootstrap_weighted_effect"] = (
        x["methylation_rna_rho"] * (2 * x["bootstrap_negative_fraction"] - 1)
    )
    md = x["tumor_normal_methylation_delta"]
    rd = x["rna_tumor_normal_delta"]
    x["rank_hypermethylated_rna_down"] = (
        (md - md.mean()) / md.std(ddof=0) - (rd - rd.mean()) / rd.std(ddof=0)
    )
    return x


def competitive_rank_test(rank_lookup, universe, total, members, tie_term):
    overlap = universe & members
    n1 = len(overlap)
    n2 = total - n1
    if n1 < MIN_SET or n1 > MAX_SET or n2 < MIN_SET:
        return None
    rank_sum = float(sum(rank_lookup[g] for g in overlap))
    u = rank_sum - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2
    variance = n1 * n2 / 12 * ((total + 1) - tie_term / (total * (total - 1)))
    if variance <= 0:
        return None
    z = (u - mean_u) / math.sqrt(variance)
    p = float(2 * stats.norm.sf(abs(z)))
    effect = float((rank_sum / n1) / total - 0.5)
    return n1, n2, effect, p, overlap


def enrich_all(results, gene_sets):
    metrics = [
        "rank_directed_neglog10_fdr", "rank_cnv_adjusted_beta",
        "rank_bootstrap_weighted_effect", "rank_hypermethylated_rna_down"
    ]
    rows = []
    for variant, block in results.groupby("promoter_definition", sort=False):
        ranked = rank_metrics(block).drop_duplicates("gene").set_index("gene")
        for metric in metrics:
            values = ranked[metric].dropna()
            rank_array = stats.rankdata(values.to_numpy(), method="average")
            ranks = pd.Series(rank_array, index=values.index)
            universe = set(values.index)
            rank_lookup = dict(zip(values.index, rank_array))
            total = len(values)
            _, tie_counts = np.unique(values.to_numpy(), return_counts=True)
            tie_term = float(np.sum(tie_counts ** 3 - tie_counts))
            for source, name, set_id, members in gene_sets:
                test = competitive_rank_test(rank_lookup, universe, total, members, tie_term)
                if test is None:
                    continue
                n_in, n_out, effect, p, overlap = test
                rows.append({
                    "source": source, "pathway": name, "pathway_id": set_id,
                    "focus_category": focus_category(name),
                    "promoter_definition": variant, "ranking_metric": metric,
                    "set_size_tested": n_in, "background_size": n_out,
                    "rank_effect": effect,
                    "enrichment_direction": "high" if effect > 0 else "low",
                    "p_value": p,
                    "member_genes_tested": ";".join(sorted(overlap)),
                    "test": "two-sided Wilcoxon/Mann-Whitney competitive rank test; normal approximation with tie correction",
                })
    out = pd.DataFrame(rows)
    out["fdr"] = np.nan
    for _, idx in out.groupby(["source", "promoter_definition", "ranking_metric"]).groups.items():
        out.loc[idx, "fdr"] = bh(out.loc[idx, "p_value"])
    out["correction_scope"] = "BH within source x promoter definition x ranking metric"
    return out


def enrich_without_rab(results, gene_sets):
    rows = []
    for excluded in [False, True]:
        data = results[results["gene"].ne("RAB27A")] if excluded else results
        part = enrich_all(data, gene_sets)
        part = part[part["focus_category"].ne("")].copy()
        part["RAB27A_excluded"] = excluded
        rows.append(part)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def leading_edge(enrichment, results):
    rows = []
    lookup = {
        (v, m): rank_metrics(b).drop_duplicates("gene").set_index("gene")[m]
        for v, b in results.groupby("promoter_definition")
        for m in ["rank_directed_neglog10_fdr", "rank_cnv_adjusted_beta",
                  "rank_bootstrap_weighted_effect", "rank_hypermethylated_rna_down"]
    }
    use = enrichment[(enrichment["fdr"] < .05) | enrichment["focus_category"].ne("")]
    for _, r in use.iterrows():
        values = lookup[(r["promoter_definition"], r["ranking_metric"])]
        genes = [g for g in str(r["member_genes_tested"]).split(";") if g in values.index]
        s = values.loc[genes].dropna().sort_values(
            ascending=r["enrichment_direction"] == "high"
        )
        n = max(1, math.ceil(len(s) * .20))
        edge = s.tail(n) if r["enrichment_direction"] == "high" else s.head(n)
        rows.append({
            "source": r["source"], "pathway": r["pathway"],
            "promoter_definition": r["promoter_definition"],
            "ranking_metric": r["ranking_metric"], "pathway_fdr": r["fdr"],
            "rank_effect": r["rank_effect"], "leading_edge_rule": "most extreme 20% in enrichment direction",
            "leading_edge_n": len(edge),
            "leading_edge_genes": ";".join(edge.index),
            "leading_edge_values": ";".join(f"{x:.6g}" for x in edge.to_numpy()),
            "RAB27A_in_pathway": "RAB27A" in genes,
            "RAB27A_in_leading_edge": "RAB27A" in edge.index,
        })
    return pd.DataFrame(rows)


def targeted_expression(project):
    resource = project / "resources" / "public_annotations" / "phase1C"
    common = pd.read_csv(
        project / "resources" / "tools" / "R-library-4.6" / "estimate" /
        "extdata" / "common_genes.txt", sep="\t", dtype=str
    )["GeneSymbol"].dropna().tolist()
    mcp = pd.read_csv(resource / "MCPcounter_genes.txt", sep="\t")
    hallmark = parse_gmt(resource / "h.all.v2026.1.Hs.symbols.gmt", "Hallmark")
    ifng = set()
    for _, name, _, genes in hallmark:
        if "INTERFERON_GAMMA_RESPONSE" in name:
            ifng |= genes
    targets = set(common) | set(mcp["HUGO symbols"].dropna()) | set(CORE_AP)
    targets |= set(CD8_GENES + NK_GENES + CYTOLYTIC_GENES + CHECKPOINT_GENES) | ifng
    proc = project / "data_processed" / "TCGA_LIHC"
    ann = pd.read_csv(proc / "TCGA_LIHC_gene_annotation.tsv", sep="\t", dtype=str)
    ann = ann[ann["gene_name"].isin(targets)]
    mapping = dict(zip(ann["gene_id"], ann["gene_name"]))
    path = proc / "TCGA_LIHC_RNA_TPM_matrix.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8") as h:
        samples = h.readline().rstrip("\n").split("\t")[1:]
    primary = unique_primary(samples)
    pieces = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=5000):
        sub = chunk[chunk["gene_id"].astype(str).isin(mapping)]
        if sub.empty:
            continue
        vals = sub[primary].apply(pd.to_numeric, errors="coerce").astype(np.float32)
        vals.index = sub["gene_id"].map(mapping).to_numpy()
        pieces.append(vals)
    expr = np.log2(pd.concat(pieces).groupby(level=0).median() + 1)
    return expr, mcp, ifng, targets


def refresh_ap_rna_score(out, expr, update_robustness=False):
    available = [g for g in CORE_AP if g in expr.index]
    if len(available) != len(CORE_AP):
        raise RuntimeError(f"TCGA AP RNA score requires all 11 fixed genes; missing {sorted(set(CORE_AP) - set(available))}")
    z = zscore_rows(expr.loc[available])
    score = z.mean(axis=0)
    result = pd.DataFrame({
        "patient_id": [patient_id(s) for s in expr.columns],
        "sample_id": expr.columns,
        "AP_RNA_output_score": score.reindex(expr.columns).to_numpy(),
        "available_core_gene_n": z.notna().sum(axis=0).reindex(expr.columns).to_numpy(),
        "available_core_genes": ";".join(available),
        "score_direction": "higher = higher antigen-presentation RNA output",
        "standardization": "gene-wise z-score across 371 unique TCGA primary tumors; fixed 11-gene core",
    })
    result.to_csv(out / "AP_RNA_output_score_TCGA.tsv", sep="\t", index=False)
    z.T.rename_axis("sample_id").to_csv(
        out / "intermediate" / "AP_RNA_gene_z_contributions.tsv.gz",
        sep="\t", compression="gzip"
    )
    if update_robustness:
        path = out / "AP_score_robustness.tsv"
        existing = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()
        if len(existing):
            existing = existing[existing["score"].ne("AP_RNA_output_score")]
        rows = []
        full = result.set_index("sample_id")["AP_RNA_output_score"]
        for gene in z.index:
            loo = z.drop(index=gene).mean(axis=0)
            frame = pd.concat([full, loo], axis=1).dropna()
            rows.append({"score": "AP_RNA_output_score", "analysis": "leave_one_gene_out",
                        "omitted_gene_or_subset": gene, "n": len(frame),
                        "spearman_vs_full_score": float(stats.spearmanr(frame.iloc[:, 0], frame.iloc[:, 1]).statistic),
                        "median_score_difference": float((frame.iloc[:, 1] - frame.iloc[:, 0]).median()),
                        "paired_group_delta": np.nan, "paired_group_p": np.nan})
        subsets = {"MHC_I": ["HLA-A", "HLA-B", "HLA-C", "B2M"],
                   "MHC_II": ["HLA-DRA", "HLA-DRB1"],
                   "processing_regulation": ["TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"]}
        for name, genes in subsets.items():
            subset = z.loc[genes].mean(axis=0)
            frame = pd.concat([full, subset], axis=1).dropna()
            rows.append({"score": "AP_RNA_output_score", "analysis": "fixed_core_subset",
                        "omitted_gene_or_subset": name, "subset_genes_present": ";".join(genes),
                        "n": len(frame), "spearman_vs_full_score": float(stats.spearmanr(frame.iloc[:, 0], frame.iloc[:, 1]).statistic),
                        "median_score_difference": np.nan, "paired_group_delta": np.nan,
                        "paired_group_p": np.nan})
        combined = pd.concat([existing, pd.DataFrame(rows)], ignore_index=True, sort=False)
        combined.to_csv(path, sep="\t", index=False)
    return result, z

def prepare_expression(project, out):
    expr, _, _, targets = targeted_expression(project)
    inter = out / "intermediate"
    inter.mkdir(parents=True, exist_ok=True)
    expr.to_csv(inter / "TCGA_primary_targeted_log2TPM.tsv.gz", sep="\t",
                compression="gzip")
    refresh_ap_rna_score(out, expr, update_robustness=False)
    expr.index.name = "GeneSymbol"
    expr.to_csv(inter / "TCGA_primary_log2TPM_gene_symbol.tsv", sep="\t")
    pd.DataFrame({
        "requested_gene_n": [len(targets)], "available_gene_n": [len(expr)],
        "sample_n": [expr.shape[1]], "transformation": ["log2(TPM+1)"],
    }).to_csv(inter / "targeted_expression_audit.tsv", sep="\t", index=False)


def mean_z_score(expr, genes):
    present = [g for g in genes if g in expr.index]
    return zscore_rows(expr.loc[present]).mean(axis=0), present


def immune_associations(project, out):
    expr = pd.read_csv(
        out / "intermediate" / "TCGA_primary_targeted_log2TPM.tsv.gz",
        sep="\t", index_col=0
    )
    mcp = pd.read_csv(
        project / "resources" / "public_annotations" / "phase1C" /
        "MCPcounter_genes.txt", sep="\t"
    )
    hallmark = parse_gmt(
        project / "resources" / "public_annotations" / "phase1C" /
        "h.all.v2026.1.Hs.symbols.gmt", "Hallmark"
    )
    ifng = next(genes for _, name, _, genes in hallmark
                if "INTERFERON_GAMMA_RESPONSE" in name)
    metrics = {}
    definitions = []
    for population, block in mcp.groupby("Cell population"):
        genes = block["HUGO symbols"].dropna().astype(str).tolist()
        present = [g for g in genes if g in expr.index]
        score = expr.loc[present].mean(axis=0)
        name = "MCPcounter_" + re.sub(r"[^A-Za-z0-9]+", "_", population).strip("_")
        metrics[name] = score
        definitions.append((name, "MCPcounter 1.1.0 HUGO-symbol estimator: arithmetic mean of log2(TPM+1) marker expression", ";".join(present)))
    for name, genes in [("CD8_T_cell_score", CD8_GENES), ("NK_cell_score", NK_GENES),
                        ("IFN_gamma_response_score", sorted(ifng)),
                        ("immune_checkpoint_composite", CHECKPOINT_GENES)]:
        score, present = mean_z_score(expr, genes)
        metrics[name] = score
        definitions.append((name, "mean gene-wise z score", ";".join(present)))
    present = [g for g in CYTOLYTIC_GENES if g in expr.index]
    metrics["cytolytic_activity"] = expr.loc[present].mean(axis=0)
    definitions.append(("cytolytic_activity", "mean log2(TPM+1)", ";".join(present)))
    for g in CHECKPOINT_GENES:
        if g in expr.index:
            metrics[f"checkpoint_{g}_log2TPM"] = expr.loc[g]
            definitions.append((f"checkpoint_{g}_log2TPM", "log2(TPM+1)", g))
    estimate = pd.read_csv(out / "intermediate" / "ESTIMATE_scores_TCGA.tsv", sep="\t")
    estimate = estimate.set_index("sample_id")
    estimate.index = estimate.index.astype(str).str.replace(".", "-", regex=False)
    for c in ["ImmuneScore", "StromalScore", "ESTIMATEScore", "TumorPurity"]:
        if c in estimate:
            metrics[f"ESTIMATE_{c}"] = estimate[c]
            definitions.append((f"ESTIMATE_{c}", "ESTIMATE 1.0.13", c))
    matrix = pd.DataFrame(metrics)
    matrix.index.name = "sample_id"
    matrix.to_csv(out / "intermediate" / "TCGA_immune_feasibility_scores.tsv.gz",
                  sep="\t", compression="gzip")
    pd.DataFrame(definitions, columns=["metric", "method", "genes_or_component"]).to_csv(
        out / "intermediate" / "immune_metric_definitions.tsv", sep="\t", index=False
    )
    ap = pd.read_csv(out / "AP_RNA_output_score_TCGA.tsv", sep="\t").set_index("sample_id")
    rows = []
    for metric in matrix:
        f = pd.concat([ap["AP_RNA_output_score"], matrix[metric]], axis=1).dropna()
        if len(f) >= 10 and f.iloc[:, 1].nunique() >= 3:
            s = stats.spearmanr(f.iloc[:, 0], f.iloc[:, 1])
            rho, p = float(s.statistic), float(s.pvalue)
        else:
            rho, p = np.nan, np.nan
        rows.append({"immune_metric": metric, "n": len(f), "spearman_rho": rho,
                     "p_value": p, "association_interpretation":
                     "correlation only; no causal interpretation"})
    result = pd.DataFrame(rows)
    result["fdr"] = bh(result["p_value"])
    result["multiple_testing_scope"] = "BH across all prespecified immune associations"
    result.to_csv(out / "AP_RNA_score_immune_associations.tsv", sep="\t", index=False)
    return result, matrix


def external_audit(project, out):
    proc = project / "data_processed" / "ICGC_LIRI_JP"
    expr_path = proc / "ICGC_LIRI_JP_expression_matrix.tsv.gz"
    expr = pd.read_csv(expr_path, sep="\t", compression="gzip", usecols=["gene_id"])
    genes = set(expr["gene_id"].astype(str))
    clinical = pd.read_csv(proc / "ICGC_LIRI_JP_clinical_survival.tsv", sep="\t",
                           low_memory=False)
    tumor_n = int(clinical["specimen_type"].astype(str).str.contains(
        "Primary tumour", case=False, na=False).sum())
    available = [g for g in CORE_AP if g in genes]
    missing = [g for g in CORE_AP if g not in genes]
    survival_cols = [c for c in ["donor_vital_status", "donor_survival_time",
                                 "donor_interval_of_last_followup"] if c in clinical]
    text = f"""# External validation feasibility report

## ICGC LIRI-JP

- Expression matrix was inspected locally: {len(genes):,} gene rows.
- Fixed AP core coverage: {len(available)}/11 genes available ({'; '.join(available)}).
- Missing core genes: {'none' if not missing else '; '.join(missing)}.
- Primary-tumour sample records in the linked clinical table: {tumor_n}.
- The same `AP_RNA_output_score` formula is computationally feasible: {'yes' if len(available) == 11 else 'only with missing-gene caveat'}.
- Clinical outcome fields present: {', '.join(survival_cols) if survival_cols else 'none identified'}.
- This is a feasibility audit only; no external association or survival validation was run.

## GEO feasibility

- [GSE63898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63898) is an HCC expression cohort; its corresponding HEPTROMIC methylation cohort is [GSE56588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56588). Their sample-level correspondence must be reconstructed and verified before they can be called a paired external methylation-expression validation cohort.
- [GSE64041](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64041) contains paired tumour/non-tumour HCC expression biopsies and is suitable for expression-direction validation, not methylation-expression pairing.
- [GSE54503](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54503) contains paired HCC/non-tumour HM450 methylation samples and can validate methylation direction only unless matched expression is identified.
- [GSE73003](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73003) is a smaller paired HCC/non-tumour 27K methylation cohort with limited promoter coverage.
- No GEO dataset was downloaded or statistically analysed in Phase 1C; therefore external validation is not claimed.

## Feasibility decision

ICGC supports a future same-formula AP RNA score validation. A plausible public paired methylation-expression route exists through the HEPTROMIC GSE63898/GSE56588 study, but sample identity and platform harmonisation remain mandatory unresolved gates.
"""
    (out / "external_validation_feasibility_report.md").write_text(text, encoding="utf-8")
    pd.DataFrame([{
        "ICGC_gene_rows": len(genes), "ICGC_core_available_n": len(available),
        "ICGC_core_missing": ";".join(missing), "ICGC_primary_tumor_records": tumor_n,
        "ICGC_score_feasible": len(available) == 11,
        "ICGC_outcome_fields": ";".join(survival_cols),
    }]).to_csv(out / "intermediate" / "external_validation_audit.tsv",
                sep="\t", index=False)


def add_kegg_manifest(project):
    path = project / "resources" / "public_annotations" / "phase1C" / "KEGG_hsa_gene_list.tsv"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    pd.DataFrame([{
        "name": "KEGG_human_gene_list", "url": "https://rest.kegg.jp/list/hsa",
        "file": path.name, "version": "live KEGG REST snapshot",
        "genome_build": "KEGG hsa identifiers; not coordinate based",
        "license_or_citation": "KEGG REST academic use",
        "download_date": "2026-07-30", "size_bytes": path.stat().st_size,
        "sha256": digest, "local_path": str(path)
    }]).to_csv(project / "phase1C_output" / "intermediate" /
               "additional_public_resource_manifest.tsv", sep="\t", index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project", required=True, type=Path)
    p.add_argument("--mode", choices=["prepare", "final"], required=True)
    args = p.parse_args()
    project = args.project.resolve()
    out = project / "phase1C_output"
    out.mkdir(exist_ok=True)
    (out / "intermediate").mkdir(exist_ok=True)
    if args.mode == "prepare":
        prepare_expression(project, out)
        add_kegg_manifest(project)
        return
    expr = pd.read_csv(out / "intermediate" / "TCGA_primary_targeted_log2TPM.tsv.gz", sep="\t", index_col=0)
    refresh_ap_rna_score(out, expr, update_robustness=True)
    results = pd.read_csv(out / "genomewide_methylation_RNA_results.tsv.gz",
                          sep="\t", compression="gzip")
    gene_sets = load_gene_sets(project)
    enrich = enrich_all(results, gene_sets)
    enrich.to_csv(out / "genomewide_ranked_pathway_enrichment.tsv",
                  sep="\t", index=False)
    sensitivity = enrich_without_rab(results, gene_sets)
    sensitivity.to_csv(out / "immune_pathway_enrichment_sensitivity.tsv",
                       sep="\t", index=False)
    leading_edge(enrich, results).to_csv(
        out / "leading_edge_gene_audit.tsv", sep="\t", index=False
    )
    immune_associations(project, out)
    external_audit(project, out)
    rab = results[results["gene"].eq("RAB27A")].copy()
    rab.to_csv(out / "intermediate" / "RAB27A_genomewide_audit.tsv",
               sep="\t", index=False)
    versions = {
        "timestamp": datetime.now().isoformat(), "python": sys.version,
        "platform": platform.platform(), "pandas": pd.__version__,
        "numpy": np.__version__, "scipy": scipy.__version__,
        "analysis_note": "Phase 1C exploratory extension; thresholds frozen before result inspection",
    }
    (out / "software_versions.json").write_text(
        json.dumps(versions, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
