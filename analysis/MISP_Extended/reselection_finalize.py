"""Finalize v5 recomputed gates, whole-source pathways and three fixed tracks.

Uses only v5 association reruns and frozen public annotation snapshots. Tracks
are reported together and never chosen by their downstream performance.
"""
from __future__ import annotations
import argparse
import json
import math
import pickle
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import reselection_core as F

TRACKS = {
    'T1_literal_crosswindow': {'correct_tail': False, 'TSS200_only': False},
    'T2_correctedtail_crosswindow': {'correct_tail': True, 'TSS200_only': False},
    'T3_correctedtail_TSS200only': {'correct_tail': True, 'TSS200_only': True},
}
IMMUNE_FOCUS = 'antigen_processing_presentation|interferon_gamma|inflammatory_response|cytotoxic_lymphocyte'


def read_associations(spec):
    return pd.concat([pd.read_csv(F.ANALYSIS / spec / f'associations_{v}.tsv.gz', sep='\t') for v in F.VARIANTS], ignore_index=True)


def historical_models():
    destination = F.ANALYSIS / 'historical_baseline'
    data = F.load_inputs()
    all_assoc = read_associations('historical_baseline')
    tiers = F.GW.classify_tiers(all_assoc)
    F.save(tiers, destination / 'recomputed_phase1C_tiers.tsv.gz')
    expected = pd.read_csv(F.PROJECT / 'phase1C_output/genomewide_methylation_candidate_tiers.tsv', sep='\t')
    assert tiers.set_index('gene').genomewide_evidence_tier.sort_index().equals(expected.set_index('gene').genomewide_evidence_tier.sort_index())
    expression, meth, cnv, cov = data['expression'], data['matrices']['TSS200'], data['cnv'], data['old_cov']
    primary = sorted(set(expression.columns) & set(meth.columns) & set(cov.index))
    cs = sorted(set(primary) & set(cnv.columns))
    high = [s for s in cs if cov.loc[s, 'purity_stratum'] == 'high']
    medium = [s for s in cs if cov.loc[s, 'purity_stratum'] == 'medium']
    absolute = [s for s in cs if np.isfinite(cov.loc[s, 'ABSOLUTE_purity'])]
    pcs = [c for c in cov if c.startswith('composition_PC') and c[len('composition_PC'):].isdigit()]
    rows = []
    started = time.monotonic()
    for j, gene in enumerate(tiers.gene):
        y, m, c = expression.loc[gene], meth.loc[gene], cnv.loc[gene]
        definitions = {
            'unadjusted': (primary, []), 'cnv': (cs, ['CNV']),
            'purity': (cs, ['CNV', 'ESTIMATE_TumorPurity']),
            'composition': (cs, ['CNV', 'ESTIMATE_TumorPurity'] + pcs),
            'high_purity': (high, ['CNV']), 'medium_purity': (medium, ['CNV']),
            'absolute_purity': (absolute, ['CNV', 'ABSOLUTE_purity']),
        }
        record = {'gene': gene}
        for name, (samples, cols) in definitions.items():
            predictors = {'methylation': m.loc[samples]}
            predictors.update({col: c.loc[samples] if col == 'CNV' else cov.loc[samples, col] for col in cols})
            fit = F.P2.fit_ols(y.loc[samples], predictors)
            record.update({name + '_' + k: v for k, v in fit.items()})
        rows.append(record)
        if (j + 1) % 1000 == 0:
            print('HISTORICAL_P2', j + 1, len(tiers), round(time.monotonic() - started), flush=True)
    fits = pd.DataFrame(rows)
    for name in ['unadjusted', 'cnv', 'purity', 'composition', 'high_purity', 'medium_purity', 'absolute_purity']:
        fits[name + '_methylation_fdr'] = F.P2.bh(fits[name + '_methylation_p'])
    merged = tiers.merge(fits, on='gene', validate='one_to_one')
    merged['purity_vs_cnv_attenuation_percent'] = 100 * (1 - merged.purity_methylation_beta.abs() / merged.cnv_methylation_beta.abs())
    merged['composition_vs_cnv_attenuation_percent'] = 100 * (1 - merged.composition_methylation_beta.abs() / merged.cnv_methylation_beta.abs())
    robust = merged.genomewide_evidence_tier.eq('Stable') & merged.purity_methylation_beta.lt(0) & merged.purity_methylation_fdr.lt(.05) & merged.purity_vs_cnv_attenuation_percent.le(50) & merged.high_purity_methylation_beta.lt(0) & merged.bootstrap_negative_fraction.ge(.7) & merged.sensitivity_any_negative_both
    merged['Composition_robust_recomputed'] = robust
    F.save(merged, destination / 'historical_phase2A_models.tsv.gz')
    old = pd.read_csv(F.PROJECT / 'phase2A_output/genomewide_purity_adjusted_models.tsv.gz', sep='\t')
    cmp = old.merge(fits, on='gene', suffixes=('_frozen', '_rerun'))
    comparisons = {}
    for name in ['cnv', 'purity', 'composition', 'high_purity', 'absolute_purity']:
        for col in ['methylation_beta', 'methylation_fdr', 'methylation_standardized_beta']:
            a = cmp[f'{name}_{col}_frozen'].to_numpy(float)
            b = cmp[f'{name}_{col}_rerun'].to_numpy(float)
            ok = np.isfinite(a) & np.isfinite(b)
            comparisons[f'{name}_{col}'] = {'max_abs_error': float(abs(a[ok] - b[ok]).max()), 'finite_pattern_mismatch': int((np.isfinite(a) != np.isfinite(b)).sum())}
    context = {'tiers': tiers.genomewide_evidence_tier.value_counts().to_dict(), 'Composition_robust_n': int(robust.sum()), 'comparisons': comparisons, 'elapsed_seconds': time.monotonic() - started}
    (destination / 'historical_phase2A_QA.json').write_text(json.dumps(context, indent=2), encoding='utf-8')
    assert robust.sum() == 1998
    print('HISTORICAL_P2_DONE', json.dumps(context), flush=True)


def new_gate_audit(spec, assoc):
    # Preserve phase1C raw-rho + same-case CNV tier logic. Only the alternative
    # promoter direction check is recomputed using its new joint coefficient.
    tier_input = assoc.copy()
    alt = tier_input.promoter_definition.ne('TSS200')
    tier_input.loc[alt, 'adjusted_methylation_coefficient'] = tier_input.loc[alt, 'joint_beta']
    tiers = F.GW.classify_tiers(tier_input)
    tiers['Composition_robust_recomputed'] = (tiers.genomewide_evidence_tier.eq('Stable') & tiers.joint_beta.lt(0) & tiers.joint_fdr.lt(.05) & tiers.attenuation_vs_same_case_CNV_percent.le(50) & tiers.high_purity_joint_beta.lt(0) & tiers.bootstrap_negative_fraction.ge(.7) & tiers.sensitivity_any_negative_both)
    tiers['high_purity_gate_estimable'] = tiers.high_purity_joint_status.eq('complete')
    tiers['selection_status'] = np.where(~tiers.high_purity_gate_estimable, 'high_purity_gate_not_estimable', np.where(tiers.Composition_robust_recomputed, 'composition_gate_pass', 'composition_gate_not_pass'))
    F.save(tiers, F.ANALYSIS / spec / 'full_gene_gate_audit.tsv.gz')
    return tiers


def ranking_input(assoc, historical):
    if historical:
        return assoc.copy()
    work = assoc.copy()
    # Exact phase2A pathway-after-purity substitutions. Reuse the frozen
    # phase1C rank function through explicit field mapping; legacy metric names
    # are retained to preserve the original pathway key rules.
    work['raw_rho_for_audit'] = work.methylation_rna_rho
    work['raw_correlation_fdr_for_audit'] = work.correlation_fdr
    work['raw_RNA_tumor_normal_delta_for_audit'] = work.rna_tumor_normal_delta
    work['methylation_rna_rho'] = work.joint_std_beta
    work['correlation_fdr'] = work.joint_fdr
    work['adjusted_methylation_coefficient'] = work.joint_beta
    work['rna_tumor_normal_delta'] = work.joint_std_beta
    return work


def leading(enrichment, ranked_input, corrected):
    lookup = {(v, metric): F.INTEG.rank_metrics(block).drop_duplicates('gene').set_index('gene')[metric]
              for v, block in ranked_input.groupby('promoter_definition')
              for metric in ['rank_directed_neglog10_fdr', 'rank_cnv_adjusted_beta', 'rank_bootstrap_weighted_effect', 'rank_hypermethylated_rna_down']}
    rows = []
    use = enrichment[enrichment.fdr.lt(.05) | enrichment.focus_category.ne('')]
    for r in use.itertuples():
        values = lookup[(r.promoter_definition, r.ranking_metric)]
        gs = [g for g in str(r.member_genes_tested).split(';') if g in values.index]
        ordered = values.loc[gs].dropna().sort_values(ascending=r.enrichment_direction == 'high')
        n = max(1, math.ceil(len(ordered) * .2))
        # Literal old implementation: low=descending head (numerically high).
        # Corrected track changes only low head -> low tail; all tie ordering
        # and the high branch remain exactly as in the historical function.
        edge = ordered.tail(n) if corrected or r.enrichment_direction == 'high' else ordered.head(n)
        rows.append({'source': r.source, 'pathway': r.pathway, 'promoter_definition': r.promoter_definition,
                     'ranking_metric': r.ranking_metric, 'window_pathway_fdr': r.fdr,
                     'window_enrichment_direction': r.enrichment_direction, 'rank_effect': r.rank_effect,
                     'leading_edge_n': len(edge), 'leading_edge_genes': ';'.join(edge.index),
                     'leading_edge_values': ';'.join(f'{x:.17g}' for x in edge.to_numpy()),
                     'tail_rule': 'directional_tail_corrected' if corrected else 'literal_historical_tail'})
    return pd.DataFrame(rows)


def pathways_and_selection(spec):
    historical = spec == 'historical_baseline'
    destination = F.ANALYSIS / spec
    assoc = read_associations(spec)
    gate = pd.read_csv(destination / 'historical_phase2A_models.tsv.gz', sep='\t') if historical else new_gate_audit(spec, assoc)
    if historical:
        F.save(gate, destination / 'full_gene_gate_audit.tsv.gz')
    rank_input = ranking_input(assoc, historical)
    ranked = pd.concat([F.INTEG.rank_metrics(block) for _, block in rank_input.groupby('promoter_definition', sort=False)], ignore_index=True)
    cols = ['gene', 'promoter_definition', 'rank_directed_neglog10_fdr', 'rank_cnv_adjusted_beta', 'rank_bootstrap_weighted_effect', 'rank_hypermethylated_rna_down']
    F.save(ranked[cols], destination / 'pathway_ranking_values.tsv.gz')
    with open(F.INPUT / 'frozen_gene_sets.pkl', 'rb') as f:
        sets = pickle.load(f)
    start = time.monotonic()
    enrichment_file = destination / 'pathway_all_tests.tsv.gz'
    if enrichment_file.exists():
        enrichment = pd.read_csv(enrichment_file, sep='\t').fillna({'focus_category': ''})
    else:
        print('PATHWAY_START', spec, len(sets), flush=True)
        enrichment = F.INTEG.enrich_all(rank_input, sets)
        enrichment['ranking_source_model'] = 'historical_phase1C' if historical else spec
        F.save(enrichment, enrichment_file)
    focus = enrichment[enrichment.promoter_definition.eq('TSS200') & enrichment.focus_category.str.contains(IMMUNE_FOCUS)]
    keys = focus[F.LOCK.expected_direction(focus) & focus.fdr.lt(.05)].copy()
    keys = keys.rename(columns={'promoter_definition': 'key_promoter_definition', 'fdr': 'key_TSS200_fdr', 'enrichment_direction': 'key_TSS200_direction'})
    F.save(keys, destination / 'eligible_TSS200_pathway_keys.tsv')
    keycols = ['source', 'pathway', 'ranking_metric']
    keyjoin = keys[keycols + ['key_TSS200_fdr', 'key_TSS200_direction']]
    icgc = set(json.loads((F.PROJECT / 'ICGC_exact_gene_symbols.json').read_text(encoding='utf-8')))
    gpl = json.loads((F.PROJECT / 'GPL13667_exact_symbol_probe_counts.json').read_text(encoding='utf-8'))
    frozen = set(pd.read_csv(F.PROJECT / 'phase2A_output/locked_immune_methylation_program_genes.tsv', sep='\t').gene)
    summaries = []
    members = []
    all_selected = []
    all_track_audit = []
    for track, rule in TRACKS.items():
        le = leading(enrichment, rank_input, rule['correct_tail'])
        join = keyjoin.merge(le, on=keycols, how='left', validate='one_to_many')
        if rule['TSS200_only']:
            join = join[join.promoter_definition.eq('TSS200')]
        F.save(join, destination / f'{track}_pathway_key_window_join.tsv.gz')
        lead_genes = set()
        for s in join.leading_edge_genes.dropna():
            lead_genes.update(str(s).split(';'))
        g = gate.copy()
        g['in_leading_edge_union'] = g.gene.isin(lead_genes)
        g['ICGC_expression_coverage'] = g.gene.isin(icgc)
        g['GPL13667_probe_n'] = g.gene.map(lambda gene: gpl.get(gene, 0))
        g['RAB27A_excluded'] = g.gene.ne('RAB27A')
        g['track'] = track
        g['specification'] = spec
        all_track_audit.append(g[['gene', 'specification', 'track', 'genomewide_evidence_tier', 'Composition_robust_recomputed', 'in_leading_edge_union', 'RAB27A_excluded', 'ICGC_expression_coverage', 'GPL13667_probe_n']])
        steps = [('TSS200 framework', np.ones(len(g), bool)),
                 ('RNA finite/detection eligible', g.rna_detection_eligible),
                 ('negative raw correlation genomewide FDR', g.methylation_rna_rho.lt(0) & g.correlation_fdr.lt(.05)),
                 ('negative same-case CNV genomewide FDR', g.adjusted_methylation_coefficient.lt(0) & g.adjusted_methylation_fdr.lt(.05)),
                 ('200-resample raw-rank negative support >=0.70', g.bootstrap_negative_fraction.ge(.7)),
                 ('alternative-promoter negative directions', g.sensitivity_any_negative_both),
                 ('Stable tier recomputed', g.genomewide_evidence_tier.eq('Stable')),
                 ('composition gate recomputed', g.Composition_robust_recomputed),
                 ('immune pathway leading edge', g.in_leading_edge_union),
                 ('RAB27A exclusion', g.RAB27A_excluded),
                 ('ICGC expression coverage', g.ICGC_expression_coverage),
                 ('GPL13667 exact-symbol coverage', g.GPL13667_probe_n.gt(0))]
        keep = np.ones(len(g), bool)
        previous = len(g)
        for order, (name, flag) in enumerate(steps):
            keep &= np.asarray(flag, bool)
            summaries.append({'specification': spec, 'track': track, 'step_order': order, 'step': name,
                              'retained_n': int(keep.sum()), 'removed_at_step': int(previous - keep.sum())})
            members.extend({'specification': spec, 'track': track, 'step_order': order, 'step': name, 'gene': gene} for gene in g.loc[keep, 'gene'])
            previous = int(keep.sum())
        chosen = g.loc[keep].copy()
        wcol = 'purity_methylation_standardized_beta' if historical else 'joint_std_beta'
        chosen['weight'] = chosen[wcol].abs() / chosen[wcol].abs().sum()
        chosen['in_original_fixed30'] = chosen.gene.isin(frozen)
        chosen['original_fixed30_primary_unchanged'] = True
        chosen['weight_source'] = 'recomputed_historical_ESTIMATE' if historical else 'recomputed_' + spec
        all_selected.append(chosen)
        F.save(chosen, destination / f'{track}_selected_genes.tsv')
        print('SELECTION', spec, track, 'pathway_keys', len(keys), 'leading_union', len(lead_genes), 'selected', len(chosen), 'original_overlap', int(chosen.in_original_fixed30.sum()), flush=True)
    F.save(pd.DataFrame(summaries), destination / 'selection_step_counts.tsv')
    F.save(pd.DataFrame(members), destination / 'selection_step_members.tsv.gz')
    F.save(pd.concat(all_track_audit, ignore_index=True), destination / 'three_track_gene_flags.tsv.gz')
    combined = pd.concat(all_selected, ignore_index=True)
    F.save(combined, destination / 'three_track_selected_genes.tsv')
    context = {'specification': spec, 'pathway_tests': len(enrichment), 'eligible_TSS200_immune_keys': len(keys),
               'track_counts': {track: int(combined.track.eq(track).sum()) if len(combined) else 0 for track in TRACKS},
               'composition_robust_n': int(gate.Composition_robust_recomputed.sum()),
               'stable_n': int(gate.genomewide_evidence_tier.eq('Stable').sum()),
               'elapsed_seconds': time.monotonic() - start,
               'high_purity_gate_status': 'not_estimable_under_original_MIN_N80' if not historical and not gate.high_purity_gate_estimable.any() else 'estimable_for_at_least_some_genes',
               'all_three_predefined_tracks_reported': True}
    if historical:
        legacy = set(combined.loc[combined.track.eq('T1_literal_crosswindow'), 'gene'])
        context['literal_fixed30_exact_match'] = legacy == frozen
        context['literal_missing_vs_frozen'] = sorted(frozen - legacy)
        context['literal_extra_vs_frozen'] = sorted(legacy - frozen)
        old_enrich = pd.read_csv(F.PROJECT / 'phase1C_output/genomewide_ranked_pathway_enrichment.tsv', sep='\t')
        ecmp = old_enrich.merge(enrichment, on=['source', 'pathway', 'promoter_definition', 'ranking_metric'], suffixes=('_frozen', '_rerun'))
        context['historical_comparison_target_scope'] = 'All tier labels; 60 sampled original gene-model rows; predeclared sampled pathway targets; original fixed30 set and 200 weight-refit outputs'
        context['pathway_max_FDR_abs_error'] = float(abs(ecmp.fdr_frozen - ecmp.fdr_rerun).max())
        assert legacy == frozen
    (destination / 'full_selection_context.json').write_text(json.dumps(context, indent=2), encoding='utf-8')
    print('FULL_SELECTION_DONE', json.dumps(context), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['historical_models', 'pathways_selection'], required=True)
    ap.add_argument('--specification', choices=['historical_baseline'] + list(F.SPECS), default='historical_baseline')
    args = ap.parse_args()
    if args.stage == 'historical_models':
        historical_models()
    else:
        pathways_and_selection(args.specification)


if __name__ == '__main__':
    main()
