"""Transfer all nine predefined selected sets; protect the original fixed30.

No external result selects a preferred track. Original cohort calibration and
global-denominator missing-contribution scoring are retained. All raw missing
standardized entries stay missing in the source arrays and are audited.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy import stats

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_transfer'

def save(frame, path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False)

def sha(path):
    import hashlib
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def zrows(frame):
    mean = frame.mean(axis=1)
    sd = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(mean, axis=0).div(sd, axis=0)

def bh_fdr(values):
    a=np.asarray(values,dtype=float); out=np.full(a.shape,np.nan)
    finite=np.isfinite(a); p=a[finite]
    if not len(p): return out
    order=np.argsort(p); adjusted=p[order]*len(p)/np.arange(1,len(p)+1)
    adjusted=np.minimum.accumulate(adjusted[::-1])[::-1]
    result=np.empty(len(p)); result[order]=np.minimum(adjusted,1)
    out[finite]=result
    return out

SPECS = ['historical_baseline', 'joint_DNA_ABSOLUTE', 'joint_DNA_ABSOLUTE_age_sex']
TRACKS = ['T1_literal_crosswindow', 'T2_correctedtail_crosswindow', 'T3_correctedtail_TSS200only']


def corr(x, y, cov=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    c = np.empty((len(x), 0)) if cov is None else np.asarray(cov, float)
    if c.ndim == 1:
        c = c[:, None]
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(c).all(axis=1)
    x, y, c = x[ok], y[ok], c[ok]
    n, k = len(x), c.shape[1]
    empty = {'n': n, 'rho': np.nan, 'p': np.nan, 'covariate_n': k, 'df': n - k - 2, 'status': 'insufficient_or_constant'}
    if n < 20 or np.std(x) == 0 or np.std(y) == 0:
        return empty
    a, b = stats.rankdata(x), stats.rankdata(y)
    if k:
        c = np.column_stack([np.ones(n)] + [stats.rankdata(c[:, i]) for i in range(k)])
        if np.linalg.matrix_rank(c) < c.shape[1]:
            empty['status'] = 'rank_deficient_covariates'
            return empty
        a -= c @ np.linalg.lstsq(c, a, rcond=None)[0]
        b -= c @ np.linalg.lstsq(c, b, rcond=None)[0]
    r = float(np.corrcoef(a, b)[0, 1])
    df = n - k - 2
    p = float(2 * stats.t.sf(abs(r) * np.sqrt(df / max(1e-16, 1 - r * r)), df))
    return {'n': n, 'rho': r, 'p': p, 'covariate_n': k, 'df': df, 'status': 'complete'}


def weighted(z, weights):
    z = np.where(np.isfinite(z), z, np.nan)
    available = np.isfinite(z)
    available_weight = np.sum(available * weights[:, None], axis=0)
    score = np.nansum(z * weights[:, None], axis=0)
    # No entirely unobserved program is turned into an observed zero score.
    score[available_weight == 0] = np.nan
    return score, available.sum(axis=0), available_weight


def main():
    import argparse
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE/'Inputs/Selected_Transfer')
    parser.add_argument('--selection-dir', type=Path, default=None, help='Directory containing specification/three_track_selected_genes.tsv; default uses supplied complete-selection outputs')
    parser.add_argument('--output-dir', type=Path, default=HERE/'results_transfer')
    args=parser.parse_args()
    inp=args.input_dir.resolve(); selection=(args.selection_dir or inp/'supplied_selection').resolve()
    OUT=args.output_dir.resolve(); OUT.mkdir(parents=True, exist_ok=True)
    started=time.monotonic()
    npz_path=inp/'standardized_program_profiles.npz'; matrices=np.load(npz_path)
    gene_names=list(matrices['genes']); lookup={g:i for i,g in enumerate(gene_names)}
    t=pd.read_csv(inp/'TCGA_metadata_and_fixed30_QA.tsv',sep='\t',index_col=0,float_precision='round_trip')
    ic=pd.read_csv(inp/'ICGC_metadata_and_fixed30_QA.tsv',sep='\t',index_col=0,float_precision='round_trip')
    geo=pd.read_csv(inp/'GEO_metadata_and_fixed30_QA.tsv',sep='\t',float_precision='round_trip')
    assert len(t)==371 and len(geo)==205 and len(ic)==231
    assert t.index.is_unique and ic.index.is_unique and geo.patient_id_normalized.is_unique and geo.GSE63898_GSM.is_unique and geo.GSE56588_sample_title.is_unique
    hcc=ic.tumour_histological_type.astype(str).eq('8170/3').to_numpy(); assert int(hcc.sum())==203
    arrays={key:matrices[key] for key in ['TCGA_RNA','TCGA_methylation','GEO_RNA','GEO_methylation','ICGC_RNA']}
    assert arrays['ICGC_RNA'].shape==(len(gene_names),231)
    with np.load(inp / 'evaluation_covariates.npz') as cv:
        covariates = {key: cv[key] for key in ('TCGA_DNA_broad', 'TCGA_DNA_ABSOLUTE')}
    assert np.isfinite(covariates['TCGA_DNA_broad']).all(axis=1).sum()==354
    assert np.isfinite(covariates['TCGA_DNA_ABSOLUTE']).all(axis=1).sum()==341
    defs=json.loads((inp/'immune_definitions.json').read_text(encoding='utf8'));eps=list(defs)
    ep_npz=np.load(inp/'immune_endpoint_expression.npz')
    raw_ep={co:pd.DataFrame(ep_npz[co+'_expression'],index=ep_npz[co+'_genes']) for co in ['TCGA','GEO','ICGC']}
    core={'TCGA':{ep:t[ep].to_numpy(float) for ep in eps},'ICGC':{'IFN_gamma_response_score':ic.IFN_gamma_response_score.to_numpy(float),'CD8_T_cell_score':ic.CD8_expression_score.to_numpy(float),'NK_cell_score':ic.NK_expression_score.to_numpy(float),'cytolytic_activity':ic.cytolytic_activity.to_numpy(float)},'GEO':{}}
    for ep in eps:
        gs=[g for g in defs[ep] if g in raw_ep['GEO'].index];x=raw_ep['GEO'].loc[gs]
        core['GEO'][ep]=(x.mean(axis=0) if ep=='cytolytic_activity' else zrows(x).mean(axis=0)).to_numpy(float)
    registry = []
    for i, sid in enumerate(t.index):
        registry.append({'cohort': 'TCGA', 'position': i, 'sample_id': sid, 'patient_or_donor_id': t.loc[sid, 'patient_id'], 'primary_analysis': True, 'calibration_n': 371, 'methylation_sample_id': sid, 'RNA_sample_id': sid})
    for i, r in enumerate(geo.itertuples()):
        registry.append({'cohort': 'GEO', 'position': i, 'sample_id': r.patient_id_normalized, 'patient_or_donor_id': r.patient_id_normalized, 'primary_analysis': True, 'calibration_n': 205, 'methylation_sample_id': r.GSE56588_sample_title, 'RNA_sample_id': r.GSE63898_GSM})
    for i, sid in enumerate(ic.index):
        registry.append({'cohort': 'ICGC', 'position': i, 'sample_id': sid, 'patient_or_donor_id': ic.loc[sid, 'icgc_donor_id'], 'primary_analysis': bool(hcc[i]), 'calibration_n': 231, 'RNA_sample_id': sid})
    registry = pd.DataFrame(registry)
    assert not registry.duplicated(['cohort', 'patient_or_donor_id']).any()
    save(registry, OUT / 'sample_registry.tsv')
    save(pd.DataFrame({'sample_id': t.index, 'IC': covariates['TCGA_DNA_ABSOLUTE'][:, 0], 'Fib': covariates['TCGA_DNA_ABSOLUTE'][:, 1], 'ABSOLUTE_purity': covariates['TCGA_DNA_ABSOLUTE'][:, 2]}), OUT / 'TCGA_evaluation_covariates.tsv')
    save(pd.DataFrame([{'cohort': co, 'position': i, 'endpoint': ep, 'value': value, 'original_core_outcome': True} for co in core for ep in eps for i, value in enumerate(core[co][ep])]), OUT / 'original_core_immune_outcomes.tsv.gz')

    paired_rows, immune_rows, gene_coverage, score_rows, overlaps, set_rows = [], [], [], [], [], []
    all_scores, all_weights = {}, []
    checks = {}
    for spec in SPECS:
        selected = pd.read_csv(selection / spec / 'three_track_selected_genes.tsv', sep='\t', float_precision='round_trip')
        for track in TRACKS:
            chosen = selected[selected.track.eq(track)].sort_values('gene')
            genes = chosen.gene.tolist()
            weights = chosen.weight.to_numpy(float, copy=True)
            assert len(genes) and np.isfinite(weights).all() and (weights > 0).all()
            assert abs(weights.sum() - 1) < 1e-12
            weights /= weights.sum()
            if any(g not in lookup for g in genes):
                raise ValueError('A selected gene is outside the supplied transfer profile panel; arbitrary new genes require new profiles and cannot be filled with zero')
            ix = [lookup[g] for g in genes]
            set_id = spec + '__' + track
            set_rows.append({'set_id': set_id, 'specification': spec, 'track': track, 'gene_n': len(genes), 'weight_sum': float(weights.sum()), 'positive_weights': True, 'selection_uses_external_outcomes': False})
            all_weights.extend({'set_id': set_id, 'gene': g, 'weight': float(w)} for g, w in zip(genes, weights))
            scores = {}
            for key, matrix in arrays.items():
                co, omics = key.split('_', 1)
                z = matrix[ix]
                value, available_n, available_weight = weighted(z, weights)
                if omics == 'RNA':
                    value = -value
                scores[key] = value
                ids = t.index if co == 'TCGA' else geo.patient_id_normalized if co == 'GEO' else ic.index
                for j, gene in enumerate(genes):
                    finite_n = int(np.isfinite(z[j]).sum())
                    gene_coverage.append({'set_id': set_id, 'cohort': co, 'omics': omics, 'gene': gene, 'gene_in_source_universe': True, 'weight': float(weights[j]), 'finite_z_n': finite_n, 'unavailable_z_n': z.shape[1] - finite_n, 'all_z_unavailable': finite_n == 0, 'calibration_n': z.shape[1], 'raw_missing_not_filled': True})
                score_rows.extend({'set_id': set_id, 'cohort': co, 'omics': omics, 'position': j, 'sample_id': ids[j], 'score': float(value[j]) if np.isfinite(value[j]) else np.nan,
                                   'available_gene_n': int(available_n[j]), 'available_weight_sum': float(available_weight[j]), 'missing_weight_sum': float(1 - available_weight[j]), 'full_global_weight_denominator': 1.0, 'score_estimable': bool(available_weight[j] > 0)} for j in range(len(value)))
            all_scores[set_id] = scores
            if spec == 'historical_baseline' and track == 'T1_literal_crosswindow':
                for key, target in [('TCGA_RNA', t.MISP_RNA_suppression_score), ('TCGA_methylation', t.MISP_methylation_score), ('GEO_RNA', geo.MISP_GEO_RNA_suppression_score), ('GEO_methylation', geo.MISP_GEO_methylation_suppression_score), ('ICGC_RNA', ic.MISP_RNA_suppression_score)]:
                    error = float(np.max(abs(scores[key] - np.asarray(target))))
                    checks[key + '_max_abs_error'] = error
                    assert error < 5e-13, (key, error)
                float64_rna = -np.nansum(matrices['TCGA_later_float64_RNA_QA'][ix] * weights[:, None], axis=0)
                checks['later_float64_TCGA_RNA_derivative_error_disclosed'] = float(np.max(abs(float64_rna - t.MISP_RNA_suppression_score)))

            outcome_values = {}
            for co in ['TCGA', 'GEO', 'ICGC']:
                for ep in eps:
                    available = [g for g in defs[ep] if g in raw_ep[co].index]
                    shared = sorted(set(available) & set(genes))
                    remaining = [g for g in available if g not in genes]
                    overlaps.append({'set_id': set_id, 'cohort': co, 'endpoint': ep, 'definition_gene_n': len(defs[ep]), 'available_gene_n': len(available), 'overlap_gene_n': len(shared), 'overlap_genes': ';'.join(shared), 'remaining_gene_n': len(remaining), 'remaining_genes': ';'.join(remaining), 'available_genes': ';'.join(available), 'standardization_calibration_n': len(core[co][ep]), 'nonoverlap_outcome_estimable': len(remaining) > 0})
                    outcome_values[(co, ep, 'original_core')] = core[co][ep]
                    if not shared:
                        nonoverlap = core[co][ep].copy()
                    elif not remaining:
                        nonoverlap = np.full(len(core[co][ep]), np.nan)
                    else:
                        x = raw_ep[co].loc[remaining]
                        nonoverlap = (x.mean(axis=0) if ep == 'cytolytic_activity' else zrows(x).mean(axis=0)).to_numpy(float)
                    outcome_values[(co, ep, 'exclude_tested_set_genes')] = nonoverlap

            for co in ['TCGA', 'GEO']:
                n = len(t) if co == 'TCGA' else len(geo)
                setups = [('raw_full', np.ones(n, bool), None)]
                if co == 'TCGA':
                    broad = covariates['TCGA_DNA_broad']
                    eligible_broad = np.isfinite(broad).all(axis=1)
                    setups += [('raw_DNA_broad_same_cases', eligible_broad, None),
                               ('DNA_broad', eligible_broad, broad[eligible_broad])]
                    joint = covariates['TCGA_DNA_ABSOLUTE']
                    eligible_joint = np.isfinite(joint).all(axis=1)
                    setups += [('raw_DNA_ABSOLUTE_same_cases', eligible_joint, None), ('DNA_ABSOLUTE', eligible_joint, joint[eligible_joint])]
                for adjustment, mask, cov in setups:
                    methyl, rna = scores[co + '_methylation'][mask], scores[co + '_RNA'][mask]
                    paired_rows.append({'set_id': set_id, 'specification': spec, 'track': track, 'gene_n': len(genes), 'cohort': co, 'adjustment': adjustment, 'endpoint': 'paired_methylation_RNA', **corr(methyl, rna, cov)})
                    for score_type, x in [('methylation', methyl), ('RNA_suppression', rna)]:
                        for rule in ['original_core', 'exclude_tested_set_genes']:
                            for ep in eps:
                                y = outcome_values[(co, ep, rule)][mask]
                                immune_rows.append({'set_id': set_id, 'specification': spec, 'track': track, 'gene_n': len(genes), 'cohort': co, 'score_type': score_type, 'adjustment': adjustment, 'overlap_rule': rule, 'endpoint': ep, **corr(x, y, cov)})
            for rule in ['original_core', 'exclude_tested_set_genes']:
                for ep in eps:
                    immune_rows.append({'set_id': set_id, 'specification': spec, 'track': track, 'gene_n': len(genes), 'cohort': 'ICGC_HCC_203', 'score_type': 'RNA_suppression', 'adjustment': 'raw_HCC203_original231_calibration', 'overlap_rule': rule, 'endpoint': ep, **corr(scores['ICGC_RNA'][hcc], outcome_values[('ICGC', ep, rule)][hcc])})
            print('TRANSFER_SET', set_id, len(genes), flush=True)

    paired = pd.DataFrame(paired_rows)
    immune = pd.DataFrame(immune_rows)
    family = ['set_id', 'cohort', 'score_type', 'adjustment', 'overlap_rule']
    immune['FDR_four_immune'] = immune.groupby(family).p.transform(bh_fdr)
    immune['planned_family_size'] = 4
    immune['finite_family_size'] = immune.groupby(family).p.transform(lambda p: int(p.notna().sum()))
    save(paired, OUT / 'paired_methylation_RNA_associations.tsv')
    save(immune, OUT / 'score_immune_associations.tsv.gz')
    save(pd.DataFrame(overlaps), OUT / 'set_endpoint_gene_overlap.tsv')
    coverage = pd.DataFrame(gene_coverage)
    save(coverage, OUT / 'gene_coverage_missingness.tsv.gz')
    save(pd.DataFrame(score_rows), OUT / 'sample_scores_available_weights.tsv.gz')
    save(pd.DataFrame(set_rows), OUT / 'set_registry.tsv')
    save(pd.DataFrame(all_weights), OUT / 'selected_gene_weights.tsv')

    # The two reported GEO intervals share paired patient resamples.
    rng = np.random.default_rng(20261115)
    draws = rng.integers(0, 205, size=(2000, 205))
    ci_rows, boot_rows = [], []
    for spec in ['joint_DNA_ABSOLUTE', 'joint_DNA_ABSOLUTE_age_sex']:
        set_id = spec + '__T3_correctedtail_TSS200only'
        scores = all_scores[set_id]
        adjustment = 'raw_full'
        estimates = np.array([corr(scores['GEO_methylation'][ix], scores['GEO_RNA'][ix])['rho'] for ix in draws])
        finite = estimates[np.isfinite(estimates)]
        estimate = paired[paired.set_id.eq(set_id) & paired.cohort.eq('GEO') & paired.adjustment.eq(adjustment)].iloc[0]
        ci_rows.append({'set_id': set_id, 'cohort': 'GEO', 'adjustment': adjustment, 'n': 205, 'rho': estimate.rho, 'p': estimate.p, 'bootstrap_reps': 2000, 'finite_bootstrap_reps': len(finite), 'seed': 20261115, 'lower': float(np.quantile(finite, .025)), 'upper': float(np.quantile(finite, .975)), 'ranks_recomputed_within_each_replicate': True, 'conditional_on_fixed_selection_weights_calibration': True})
        boot_rows.extend({'set_id': set_id, 'adjustment': adjustment, 'replicate': i + 1, 'rho': value} for i, value in enumerate(estimates))
    save(pd.DataFrame(ci_rows), OUT / 'strict_T3_GEO_bootstrap_2000_CI.tsv')
    save(pd.DataFrame(boot_rows), OUT / 'strict_T3_GEO_bootstrap_draw_estimates.tsv.gz')
    np.savez_compressed(OUT / 'strict_T3_GEO_bootstrap_indices.npz', indices=draws)

    # Original primary-model endpoint checks use the untouched legacy cores.
    fixed_id = 'historical_baseline__T1_literal_crosswindow'
    checks['all_original_score_reconstruction_passed_strict_5e13_tolerance'] = True
    checks['missing_values_not_raw_zero_imputed'] = True
    checks['entirely_unavailable_gene_profiles'] = int(coverage.all_z_unavailable.sum())
    checks['entirely_unavailable_sample_scores'] = int((~pd.DataFrame(score_rows).score_estimable).sum())
    checks['elapsed_seconds'] = time.monotonic() - started
    (OUT / 'reconstruction_and_scope_checks.json').write_text(json.dumps(checks, indent=2), encoding='utf-8')
    sources=sorted(p for p in inp.rglob('*') if p.is_file())
    (OUT/'input_provenance.json').write_text(json.dumps([{'path':p.relative_to(inp).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sources],indent=2),encoding='utf8')
    (OUT/'selection_source.json').write_text(json.dumps({'selection_dir_role':'Complete-selection output supplied by the user or accompanying package; weights are inputs to transfer, never fitted from external outcomes','all_nine_sets_reported':True,'profile_gene_n':len(gene_names),'calibration':{'TCGA':371,'GEO':205,'ICGC':231,'ICGC_HCC_inference':203},'input_profile_scope':'Union of the supplied nine selected programs; genes outside this panel cause an error'},indent=2),encoding='utf8')
    print('TRANSFER_DONE', json.dumps(checks), flush=True)
    print(paired[paired.cohort.eq('GEO')][['set_id', 'adjustment', 'n', 'rho', 'p']].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
