"""Atomic gate exports and original fixed-selected-set robustness branches.

The 200 raw-rank bootstraps used for gene eligibility have already run across
all genes. This separate 200-refit procedure follows phase2A_lock_program and
conditions on each selected set; it is not bootstrap reselection frequency.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy import stats

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import reselection_core as F
from reselection_finalize import TRACKS


def atomic_gates(spec):
    directory = F.ANALYSIS / spec
    historical = spec == 'historical_baseline'
    gate = pd.read_csv(directory / 'full_gene_gate_audit.tsv.gz', sep='\t')
    flags = pd.read_csv(directory / 'three_track_gene_flags.tsv.gz', sep='\t')
    selected = pd.read_csv(directory / 'three_track_selected_genes.tsv', sep='\t')
    count_rows, member_rows = [], []
    for track in TRACKS:
        g = gate.merge(flags[flags.track.eq(track)][['gene', 'in_leading_edge_union', 'RAB27A_excluded', 'ICGC_expression_coverage', 'GPL13667_probe_n']], on='gene', validate='one_to_one')
        joint_negative = (g.purity_methylation_beta.lt(0) & g.purity_methylation_fdr.lt(.05)) if historical else (g.joint_beta.lt(0) & g.joint_fdr.lt(.05))
        attenuation = g.purity_vs_cnv_attenuation_percent.le(50) if historical else g.attenuation_vs_same_case_CNV_percent.le(50)
        hp_estimable = g.high_purity_status.eq('complete') if historical else g.high_purity_joint_status.eq('complete')
        hp_negative = g.high_purity_methylation_beta.lt(0) if historical else g.high_purity_joint_beta.lt(0)
        steps = [
            ('TSS200 framework', np.ones(len(g), bool)),
            ('RNA finite>=0.80/detected>=0.20', g.rna_detection_eligible),
            ('raw negative rho and genomewide FDR<0.05', g.methylation_rna_rho.lt(0) & g.correlation_fdr.lt(.05)),
            ('same-case CNV negative beta and genomewide FDR<0.05', g.adjusted_methylation_coefficient.lt(0) & g.adjusted_methylation_fdr.lt(.05)),
            ('raw-rank bootstrap negative fraction>=0.70', g.bootstrap_negative_fraction.ge(.7)),
            ('alternative-window negative raw and adjusted directions', g.sensitivity_any_negative_both),
            ('Stable tier recomputed', g.genomewide_evidence_tier.eq('Stable')),
            ('new adjusted-model negative beta and genomewide FDR<0.05' if not historical else 'ESTIMATE model negative beta and genomewide FDR<0.05', joint_negative),
            ('attenuation<=50 percent versus CNV baseline', attenuation),
            ('high-purity model estimable under original n>=80', hp_estimable),
            ('high-purity adjusted beta negative', hp_negative),
            ('Composition-robust recomputed', g.Composition_robust_recomputed),
            ('eligible immune pathway leading-edge union', g.in_leading_edge_union),
            ('RAB27A excluded', g.RAB27A_excluded),
            ('ICGC expression technical coverage', g.ICGC_expression_coverage),
            ('GPL13667 exact-symbol technical coverage', g.GPL13667_probe_n.gt(0)),
        ]
        keep = np.ones(len(g), bool)
        previous = len(g)
        for i, (name, flag) in enumerate(steps):
            keep &= np.asarray(flag, bool)
            count_rows.append({'specification': spec, 'track': track, 'step_order': i, 'step': name, 'retained_n': int(keep.sum()), 'removed_at_step': int(previous - keep.sum()), 'isolated_flag_n': int(np.asarray(flag, bool).sum())})
            member_rows.extend({'specification': spec, 'track': track, 'step_order': i, 'gene': gene} for gene in g.loc[keep, 'gene'])
            previous = int(keep.sum())
        assert set(g.loc[keep, 'gene']) == set(selected.loc[selected.track.eq(track), 'gene'])
    F.save(pd.DataFrame(count_rows), directory / 'selection_atomic_step_counts.tsv')
    F.save(pd.DataFrame(member_rows), directory / 'selection_atomic_step_members.tsv.gz')
    print('ATOMIC_GATES', spec, 'all selected sets match', flush=True)


def coefficient_weight(y, meth, cnv, nuisance):
    ok = np.isfinite(y) & np.isfinite(meth) & np.isfinite(cnv) & np.isfinite(nuisance).all(axis=1)
    y, meth, cnv, nuisance = y[ok], meth[ok], cnv[ok], nuisance[ok]
    if len(y) < 80 or len(np.unique(y)) < 3 or len(np.unique(meth)) < 3:
        return np.nan
    x = np.column_stack([np.ones(len(y)), meth, cnv, nuisance])
    if np.linalg.matrix_rank(x) < x.shape[1]:
        return np.nan
    b = np.linalg.lstsq(x, y, rcond=None)[0]
    return float(abs(b[1] * np.std(meth) / np.std(y)))


def robust_scores(spec):
    historical = spec == 'historical_baseline'
    directory = F.ANALYSIS / spec
    selected = pd.read_csv(directory / 'three_track_selected_genes.tsv', sep='\t')
    if selected.empty:
        (directory / 'selected_score_robustness_context.json').write_text(json.dumps({'specification': spec, 'status': 'not_applicable_no_gate_passing_program', 'reason': 'high-purity gate not estimable at original MIN_N80; not evidence of no biological association'}, indent=2), encoding='utf-8')
        return
    d = F.load_inputs()
    expression, beta, cnv = d['expression'], d['matrices']['TSS200'], d['cnv']
    cv = d['old_cov'] if historical else d['cov']
    primary = sorted(set(expression.columns) & set(beta.columns) & set(cv.index))
    model_samples = sorted(set(primary) & set(cnv.columns))
    nuisance_cols = ['ESTIMATE_TumorPurity'] if historical else F.SPECS[spec]
    model_samples = [s for s in model_samples if np.isfinite(cv.loc[s, nuisance_cols].to_numpy(float)).all()]
    nuisance = cv.loc[model_samples, nuisance_cols].to_numpy(float)
    all_rows, all_weights, score_rows = [], [], []
    started = time.monotonic()
    for track in TRACKS:
        chosen = selected[selected.track.eq(track)].copy()
        if chosen.empty:
            continue
        genes = chosen.gene.tolist()
        fixed_w = chosen.weight.to_numpy(float)
        mz = F.LOCK.zrows(beta.reindex(index=genes, columns=primary)).to_numpy(float)
        rz = F.LOCK.zrows(expression.reindex(index=genes, columns=primary)).to_numpy(float)
        full_m = np.nansum(mz * fixed_w[:, None], axis=0)
        full_r = -np.nansum(rz * fixed_w[:, None], axis=0)
        score_rows.extend({'specification': spec, 'track': track, 'sample_id': s, 'methylation_score': full_m[i], 'RNA_suppression_score': full_r[i], 'calibration': 'original 371 paired primary tumours', 'missing_rule': 'global denominator retained; absent standardized contribution zero'} for i, s in enumerate(primary))
        y_array = expression.reindex(index=genes, columns=model_samples).to_numpy(float)
        m_array = beta.reindex(index=genes, columns=model_samples).to_numpy(float)
        c_array = cnv.reindex(index=genes, columns=model_samples).to_numpy(float)
        # Same original seed and shared patient draws, reset for each fixed set.
        rng = np.random.default_rng(20260731)
        for rep in range(200):
            draw = rng.choice(len(model_samples), len(model_samples), replace=True)
            bw = np.array([coefficient_weight(y_array[j, draw], m_array[j, draw], c_array[j, draw], nuisance[draw]) for j in range(len(genes))])
            valid = np.isfinite(bw)
            if valid.sum() < 2 or bw[valid].sum() <= 0:
                all_rows.append({'specification': spec, 'track': track, 'analysis': 'bootstrap_weight_refit', 'replicate': rep + 1, 'status': 'insufficient_finite_weights', 'remaining_gene_n': int(valid.sum())})
                continue
            bw[~valid] = 0
            bw /= bw.sum()
            ms = np.nansum(mz * bw[:, None], axis=0)
            rs = -np.nansum(rz * bw[:, None], axis=0)
            all_rows.append({'specification': spec, 'track': track, 'analysis': 'bootstrap_weight_refit', 'replicate': rep + 1, 'status': 'complete', 'remaining_gene_n': int(valid.sum()), 'methylation_score_spearman_vs_full': float(stats.spearmanr(ms, full_m).statistic), 'RNA_score_spearman_vs_full': float(stats.spearmanr(rs, full_r).statistic)})
            all_weights.extend({'specification': spec, 'track': track, 'replicate': rep + 1, 'gene': gene, 'weight': float(bw[j]), 'finite_coefficient': bool(valid[j])} for j, gene in enumerate(genes))
        for j, gene in enumerate(genes):
            w = fixed_w.copy()
            w[j] = 0
            w /= w.sum()
            ms, rs = np.nansum(mz * w[:, None], axis=0), -np.nansum(rz * w[:, None], axis=0)
            all_rows.append({'specification': spec, 'track': track, 'analysis': 'leave_one_gene_out', 'omitted': gene, 'status': 'complete', 'remaining_gene_n': len(genes) - 1, 'methylation_score_spearman_vs_full': float(stats.spearmanr(ms, full_m).statistic), 'RNA_score_spearman_vs_full': float(stats.spearmanr(rs, full_r).statistic)})
        provenance = pd.read_csv(directory / f'{track}_pathway_key_window_join.tsv.gz', sep='\t')
        for (source, pathway), group in provenance.groupby(['source', 'pathway']):
            omitted = set()
            for text in group.leading_edge_genes.dropna():
                omitted.update(str(text).split(';'))
            omitted &= set(genes)
            if not omitted:
                continue
            remain = np.array([gene not in omitted for gene in genes])
            if remain.sum() < 2:
                continue
            w = fixed_w.copy()
            w[~remain] = 0
            w /= w.sum()
            ms, rs = np.nansum(mz * w[:, None], axis=0), -np.nansum(rz * w[:, None], axis=0)
            all_rows.append({'specification': spec, 'track': track, 'analysis': 'leave_one_pathway_out', 'omitted': source + '|' + pathway, 'status': 'complete', 'remaining_gene_n': int(remain.sum()), 'omitted_gene_n': len(omitted), 'methylation_score_spearman_vs_full': float(stats.spearmanr(ms, full_m).statistic), 'RNA_score_spearman_vs_full': float(stats.spearmanr(rs, full_r).statistic)})
        for variant in ['TSS1500', 'promoter_union']:
            present = [g for g in genes if g in d['matrices'][variant].index]
            ids = [genes.index(g) for g in present]
            w = fixed_w[ids]
            w /= w.sum()
            az = F.LOCK.zrows(d['matrices'][variant].loc[present, primary]).to_numpy(float)
            ms = np.nansum(az * w[:, None], axis=0)
            all_rows.append({'specification': spec, 'track': track, 'analysis': 'promoter_definition_sensitivity', 'omitted': variant, 'status': 'complete', 'remaining_gene_n': len(present), 'methylation_score_spearman_vs_full': float(stats.spearmanr(ms, full_m).statistic)})
        print('WEIGHT_REFIT', spec, track, len(genes), '200 completed', round(time.monotonic() - started), flush=True)
    result = pd.DataFrame(all_rows)
    F.save(result, directory / 'selected_score_robustness.tsv')
    F.save(pd.DataFrame(all_weights), directory / 'selected_bootstrap_refit_weights.tsv.gz')
    F.save(pd.DataFrame(score_rows), directory / 'selected_scores_TCGA.tsv')
    context = {'specification': spec, 'bootstrap_repetitions_per_fixed_set': 200, 'bootstrap_model_n': len(model_samples), 'bootstrap_seed': 20260731,
               'gene_bootstrap_is_separate': '200 raw-rank resamples were rerun across all genes for eligibility',
               'weight_bootstrap_scope': 'Fixed selected genes, refit coefficients/weights only; not full selection repeated 200 times',
               'elapsed_seconds': time.monotonic() - started, 'complete_bootstrap_by_track': result[result.analysis.eq('bootstrap_weight_refit') & result.status.eq('complete')].groupby('track').size().to_dict()}
    if historical:
        old = pd.read_csv(F.PROJECT / 'phase2A_output/MISP_score_robustness.tsv', sep='\t')
        new = result[result.track.eq('T1_literal_crosswindow') & result.analysis.eq('bootstrap_weight_refit')].sort_values('replicate')
        frozen = old[old.analysis.eq('bootstrap_weight_refit')].copy()
        frozen['replicate'] = frozen.omitted.str.extract(r'(\d+)', expand=False).astype(int)
        frozen = frozen.sort_values('replicate')
        context['historical_bootstrap_methyl_score_max_abs_error'] = float(abs(new.methylation_score_spearman_vs_full.to_numpy() - frozen.methylation_score_spearman_vs_full.to_numpy()).max())
        context['historical_bootstrap_RNA_score_max_abs_error'] = float(abs(new.RNA_score_spearman_vs_full.to_numpy() - frozen.RNA_score_spearman_vs_full.to_numpy()).max())
    (directory / 'selected_score_robustness_context.json').write_text(json.dumps(context, indent=2), encoding='utf-8')
    print('ROBUSTNESS_DONE', json.dumps(context), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--specification', choices=['historical_baseline'] + list(F.SPECS), required=True)
    ap.add_argument('--stage', choices=['atomic', 'robustness', 'all'], default='all')
    args = ap.parse_args()
    if args.stage in ['atomic', 'all']:
        atomic_gates(args.specification)
    if args.stage in ['robustness', 'all']:
        robust_scores(args.specification)


if __name__ == '__main__':
    main()
