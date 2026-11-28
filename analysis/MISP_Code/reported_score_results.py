"""Recalculate reported fixed-score bootstraps, immune associations and R-squared results."""

import json

import numpy as np
import pandas as pd

from analysis_core import bh, corr, ols, zrows


def bootstrap(x, y, covariates, repetitions, rng):
    x, y = np.asarray(x, float), np.asarray(y, float)
    covariates = np.asarray(covariates, float)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(covariates).all(axis=1)
    x, y, covariates = x[valid], y[valid], covariates[valid]
    estimates = []
    for draw in rng.integers(0, len(x), size=(repetitions, len(x))):
        estimates.append(corr(x[draw], y[draw], covariates[draw])['rho'])
    lower, upper = np.nanquantile(estimates, (0.025, 0.975))
    return {'lower': lower, 'upper': upper, 'bootstrap_reps': repetitions}


def run(input_dir, output_dir, tcga, geo, icgc, hcc, expression, definitions, genes):
    def save(rows, name):
        pd.DataFrame(rows).to_csv(output_dir / name, sep='\t', index=False)

    tcga_methyl = tcga.MISP_methylation_score.to_numpy()
    tcga_rna = tcga.MISP_RNA_suppression_score.to_numpy()
    geo_methyl = geo.MISP_GEO_methylation_suppression_score.to_numpy()
    geo_rna = geo.MISP_GEO_RNA_suppression_score.to_numpy()
    dna = pd.read_csv(input_dir / 'TCGA_DNA_composition_proxies.tsv', sep='\t', index_col=0).reindex(tcga.index)
    broad = dna[['IC', 'Fib']].to_numpy()
    hierarchical = dna[['Fib'] + [c for c in dna if c.startswith('absolute_proxy_')]].to_numpy()

    fixed_pairs = [('TCGA', tcga_methyl, tcga_rna), ('GEO', geo_methyl, geo_rna)]
    rows = []
    for cohort, x, y in fixed_pairs:
        rows.append({'cohort': cohort, **corr(x, y), **bootstrap(
            x, y, np.empty((len(x), 0)), 5000, np.random.default_rng(20261113))})
    save(rows, 'S13_paired_scores.tsv')

    eligible = np.isfinite(broad).all(axis=1)
    x, y = tcga_methyl[eligible], tcga_rna[eligible]
    rng = np.random.default_rng(20261115)
    rows = []
    for cohort, label, xvalues, yvalues, cov in [
        ('TCGA', 'eligible_unadjusted', x, y, np.empty((len(x), 0))),
        ('TCGA', 'DNA_broad', x, y, broad[eligible]),
        ('GEO', 'eligible_unadjusted', geo_methyl, geo_rna, np.empty((len(geo), 0))),
    ]:
        rows.append({'cohort': cohort, 'adjustment': label, **corr(xvalues, yvalues, cov),
                     **bootstrap(xvalues, yvalues, cov, 2000, rng)})
    save(rows, 'S21a_paired_score_sensitivity.tsv')

    rows = []
    icgc_cols = {'CD8_T_cell_score': 'CD8_expression_score', 'NK_cell_score': 'NK_expression_score'}
    for cohort, use in [('HCC_203', hcc), ('mixed_231', np.ones(len(icgc), bool))]:
        x = icgc.loc[use, 'MISP_RNA_suppression_score'].to_numpy()
        for endpoint in definitions:
            y = icgc.loc[use, icgc_cols.get(endpoint, endpoint)].to_numpy()
            rows.append({'cohort': cohort, 'endpoint': endpoint, **corr(x, y), **bootstrap(
                x, y, np.empty((len(x), 0)), 5000, np.random.default_rng(20261113))})
    output = pd.DataFrame(rows)
    output['FDR_four_immune'] = output.groupby('cohort').p.transform(bh)
    save(output, 'S15_ICGC_immune_correlations.tsv')

    inflamed_genes = pd.read_csv(input_dir / 'Inflamed_gene_list.tsv', sep='\t').gene.tolist()
    excluded = set(genes) | set(inflamed_genes)
    associations, variance = [], []
    for cohort, table, methyl, rna, use in [
        ('TCGA', tcga, tcga_methyl, tcga_rna, np.ones(len(tcga), bool)),
        ('GEO', geo, geo_methyl, geo_rna, np.ones(len(geo), bool)),
        ('ICGC_HCC_203', icgc, None, icgc.MISP_RNA_suppression_score.to_numpy(), hcc),
    ]:
        expr = expression['ICGC' if cohort == 'ICGC_HCC_203' else cohort]
        available = [g for g in inflamed_genes if g in expr.index and expr.loc[g].notna().any()]
        inflamed = zrows(expr.loc[available]).mean(axis=0).to_numpy()
        adjustments = {'unadjusted': np.empty((len(table), 0))}
        if cohort == 'TCGA':
            adjustments.update(DNA_broad=broad, DNA_hierarchical=hierarchical)

        for endpoint, members in definitions.items():
            for rule in ('original', 'exclude_MISP_and_Inflamed'):
                gs = [g for g in members if g in expr.index and (rule == 'original' or g not in excluded)]
                values = expr.loc[gs]
                outcome = (values if endpoint == 'cytolytic_activity' else zrows(values)).mean(axis=0).to_numpy()
                if cohort == 'TCGA' and rule == 'original':
                    outcome = table[endpoint].to_numpy()
                row = {'cohort': cohort, 'outcome': endpoint, 'overlap_rule': rule, 'outcome_genes': len(gs)}

                if methyl is not None:
                    for model, cov in adjustments.items():
                        associations.append({**row, 'adjustment': model, **corr(methyl, outcome, cov)})
                        variance.append({**row, 'adjustment': model, **ols(outcome, methyl, cov)})
                if cohort == 'TCGA':
                    variance.append({**row, 'adjustment': 'DNA_broad_plus_Inflamed',
                                     **ols(outcome, methyl, np.column_stack((broad, inflamed)))})
                if rule != 'original' and cohort != 'GEO':
                    cov = np.column_stack((broad, inflamed)) if cohort == 'TCGA' else inflamed[:, None]
                    variance.append({**row, 'adjustment': 'Inflamed_plus_DNA_if_available_add_MISP_RNA',
                                     **ols(outcome[use], rna[use], cov[use])})
    for name, rows in [('S21b_methylation_immune_associations.tsv', associations),
                       ('S23_incremental_immune_variance.tsv', variance)]:
        output = pd.DataFrame(rows)
        output['FDR_within_family'] = output.groupby(['cohort', 'overlap_rule', 'adjustment']).p.transform(bh)
        save(output, name)
