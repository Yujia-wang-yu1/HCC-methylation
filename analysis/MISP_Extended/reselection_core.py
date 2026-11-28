"""Local v5 full gate reruns. Source project and v1-v4 are strictly read-only.

Historical functions and thresholds are imported, never invoked against their
original output directory. New-model substitutions are specified in Review/
DNA_reselection_execution_plan.md before results are inspected.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import pickle
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy import stats

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PREPARED = Path(os.environ.get('MISP_REPRO_INPUTS', ROOT.parent)).resolve() / 'MISP_Extended'
PROJECT = PREPARED / 'Vendor/reference_project'
V4 = ROOT
VARIANTS = ['TSS200', 'TSS1500', 'promoter_union']
INPUT = PREPARED / 'Inputs'
ANALYSIS = Path(os.environ.get('MISP_REPRO_OUTPUT', ROOT.parent)).resolve() / 'results_reselection'
REVIEW = ANALYSIS / 'QA'


def module(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sys.modules[name] = m
    sp.loader.exec_module(m)
    return m


GW = module('frozen_phase1c', ROOT / 'Vendor/phase1C_genomewide.py')
INTEG = module('frozen_integrate', ROOT / 'Vendor/phase1C_integrate.py')
P2 = module('frozen_phase2a', ROOT / 'Vendor/phase2A_tcga_models.py')
LOCK = module('frozen_lock', ROOT / 'Vendor/phase2A_lock_program.py')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False)


def prepare():
    raise RuntimeError('This package starts from the supplied full processed matrices; raw-data preparation is outside its scope.')


def load_inputs():
    with open(INPUT / 'full_gene_inputs.pkl', 'rb') as f:
        return pickle.load(f)


def historical_variant(variant):
    d = load_inputs()
    start = time.monotonic()
    result, context = GW.compute_variant_associations(variant, d['matrices'][variant], d['expression'], d['cnv'], d['probe_counts'][variant])
    destination = ANALYSIS / 'historical_baseline'
    save(result, destination / f'associations_{variant}.tsv.gz')
    expected = pd.read_csv(PROJECT / f'phase1C_output/intermediate/association_{variant}.tsv.gz', sep='\t')
    key = ['gene', 'promoter_definition']
    merged = expected.merge(result, on=key, suffixes=('_frozen', '_rerun'), validate='one_to_one')
    comparisons = {}
    for col in ['methylation_rna_rho', 'correlation_p', 'correlation_fdr', 'bootstrap_negative_fraction',
                'bootstrap_rho_q025', 'bootstrap_rho_q975', 'adjusted_methylation_coefficient', 'adjusted_methylation_p', 'adjusted_methylation_fdr']:
        a = merged[col + '_frozen'].to_numpy(float)
        b = merged[col + '_rerun'].to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        comparisons[col] = {'max_abs_error': float(np.max(abs(a[ok] - b[ok]))) if ok.any() else None,
                            'finite_pattern_mismatch': int((np.isfinite(a) != np.isfinite(b)).sum())}
    context.update({'comparison_target_scope': '60 predeclared sampled genes per window, not all-genome target estimates', 'elapsed_seconds': time.monotonic() - start, 'comparisons': comparisons})
    (destination / f'context_{variant}.json').write_text(json.dumps(context, indent=2), encoding='utf-8')
    print('HISTORICAL_DONE', variant, json.dumps(context), flush=True)


def fast_fit(y, meth, cov, min_n=80, cov_names=None, require_three_unique=True):
    """Same OLS definition as frozen code, using NumPy arrays after case selection."""
    y = np.asarray(y, float)
    meth = np.asarray(meth, float)
    cov = np.asarray(cov, float)
    if cov.ndim == 1:
        cov = cov[:, None]
    ok = np.isfinite(y) & np.isfinite(meth) & np.isfinite(cov).all(axis=1)
    y, meth, cov = y[ok], meth[ok], cov[ok]
    n = len(y)
    out = {'n': n, 'status': 'insufficient_data', 'beta': np.nan, 'se': np.nan,
           'p': np.nan, 'std_beta': np.nan, 'lower': np.nan, 'upper': np.nan, 'r2': np.nan}
    if n < min_n or (require_three_unique and (len(np.unique(y)) < 3 or len(np.unique(meth)) < 3)):
        return out
    # As with a categorical factor model, omit unobserved levels and select an
    # observed reference level in that gene's complete cases. Do not remove a
    # constant CNV predictor: historical CNV-rank deficiency remains excluded.
    omitted = []
    if cov_names is not None:
        keep_cov = np.ones(cov.shape[1], dtype=bool)
        for factor in ['stage_numeric', 'grade_numeric']:
            ids = [i for i, name in enumerate(cov_names) if name.startswith(factor + '_level')]
            for i in ids:
                if not np.any(cov[:, i]):
                    keep_cov[i] = False
                    omitted.append(cov_names[i])
            active = [i for i in ids if keep_cov[i]]
            if active and np.all(cov[:, active].sum(axis=1) == 1):
                keep_cov[active[0]] = False
                omitted.append(cov_names[active[0]] + '_observed_reference')
        cov = cov[:, keep_cov]
    x = np.column_stack([np.ones(n), meth, cov])
    if np.linalg.matrix_rank(x) < x.shape[1]:
        out['status'] = 'rank_deficient'
        return out
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    df = n - x.shape[1]
    inverse = np.linalg.pinv(x.T @ x) if require_three_unique else np.linalg.inv(x.T @ x)
    covariance = float(residual @ residual / df) * inverse
    variance = max(covariance[1, 1], 0) if require_three_unique else covariance[1, 1]
    se = float(np.sqrt(variance))
    p = float(2 * stats.t.sf(abs(beta[1] / se), df))
    q = stats.t.ppf(.975, df)
    total = float(((y - y.mean()) ** 2).sum())
    out.update({'status': 'complete', 'beta': float(beta[1]), 'se': se, 'p': p,
                'std_beta': float(beta[1] * np.std(meth) / np.std(y)) if np.std(y) > 0 else np.nan,
                'lower': float(beta[1] - q * se), 'upper': float(beta[1] + q * se),
                'r2': float(1 - residual @ residual / total) if total > 0 else np.nan,
                'omitted_factor_columns': ';'.join(omitted)})
    return out


def raw_bootstrap(x, y, gene, variant):
    """Literal frozen raw-rank bootstrap: rank before, not within, resampling."""
    n = len(x)
    if n < 100:
        return np.nan, np.nan, np.nan
    xr = stats.rankdata(x, method='average')
    yr = stats.rankdata(y, method='average')
    rng = np.random.default_rng(GW.stable_gene_seed(gene, variant))
    indices = rng.integers(0, n, size=(200, n))
    xb, yb = xr[indices], yr[indices]
    xb -= xb.mean(axis=1, keepdims=True)
    yb -= yb.mean(axis=1, keepdims=True)
    den = np.sqrt((xb * xb).sum(axis=1) * (yb * yb).sum(axis=1))
    rho = np.divide((xb * yb).sum(axis=1), den, out=np.full(200, np.nan), where=den > 0)
    finite = rho[np.isfinite(rho)]
    if not len(finite):
        return np.nan, np.nan, np.nan
    return float(np.mean(finite < 0)), float(np.quantile(finite, .025)), float(np.quantile(finite, .975))


SPECS = {
    'joint_DNA_ABSOLUTE': ['IC', 'Fib', 'ABSOLUTE_purity'],
    'joint_DNA_ABSOLUTE_age_sex': ['IC', 'Fib', 'ABSOLUTE_purity', 'age_years', 'male'],
    'joint_DNA_ABSOLUTE_all_clinical': ['IC', 'Fib', 'ABSOLUTE_purity', 'age_years', 'male',
                                     'stage_numeric_level2', 'stage_numeric_level3', 'stage_numeric_level4',
                                     'grade_numeric_level2', 'grade_numeric_level3', 'grade_numeric_level4'],
}


def new_variant(spec, variant):
    d = load_inputs()
    cv = d['cov']
    beta = d['matrices'][variant]
    expression = d['expression']
    cnv = d['cnv']
    primary = sorted(set(beta.columns) & set(expression.columns) & set(cv.index))
    samples = [s for s in primary if s in cnv.columns]
    cov = cv.loc[samples, SPECS[spec]].to_numpy(float)
    high = cv.loc[samples, 'high_original_ESTIMATE'].to_numpy(bool)
    low_missing_cov = np.isfinite(cov).all(axis=1)
    genes = sorted(set(beta.index) & set(expression.index))
    expression_array = expression.reindex(index=genes, columns=samples).to_numpy(float)
    methyl_array = beta.reindex(index=genes, columns=samples).to_numpy(float)
    cnv_array = cnv.reindex(index=genes, columns=samples).to_numpy(float)
    meth_tumor = GW.unique_patient_samples(list(beta.columns), '01')
    meth_normal = GW.unique_patient_samples(list(beta.columns), '11')
    rna_tumor = GW.unique_patient_samples(list(expression.columns), '01')
    rna_normal = GW.unique_patient_samples(list(expression.columns), '11')
    eligible_patients = {str(s)[:12] for s in np.asarray(samples)[low_missing_cov]}
    meth_tumor = {p: s for p, s in meth_tumor.items() if p in eligible_patients}
    rna_tumor = {p: s for p, s in rna_tumor.items() if p in eligible_patients}
    rows = []
    start = time.monotonic()
    for j, gene in enumerate(genes):
        y, m, c = expression_array[j], methyl_array[j], cnv_array[j]
        ok = low_missing_cov & np.isfinite(y) & np.isfinite(m) & np.isfinite(c)
        yv, mv, cn = y[ok], m[ok], c[ok]
        n = len(yv)
        # RNA eligibility is recalculated in the new covariate-eligible cohort,
        # before gene-specific methylation/CNV removal, as in the original rule.
        finite_fraction = float(np.isfinite(y[low_missing_cov]).mean())
        detected_fraction = float((np.nan_to_num(y[low_missing_cov], nan=0) > 0).mean())
        eligible = finite_fraction >= .8 and detected_fraction >= .2
        rho = p = np.nan
        bf = blo = bhi = np.nan
        if eligible and n >= 100:
            if len(np.unique(mv)) >= 3 and len(np.unique(yv)) >= 3:
                cor = stats.spearmanr(mv, yv)
                rho, p = float(cor.statistic), float(cor.pvalue)
            bf, blo, bhi = raw_bootstrap(mv, yv, gene, variant)
        # Phase1C CNV baseline checks n and design rank, but unlike Phase2A
        # joint models does not impose a three-unique-value threshold.
        base = fast_fit(yv, mv, cn, min_n=100, require_three_unique=False)
        cov_names = ['CNV'] + SPECS[spec]
        joint = fast_fit(yv, mv, np.column_stack([cn, cov[ok]]), min_n=80, cov_names=cov_names)
        h = high[ok]
        hp = fast_fit(yv[h], mv[h], np.column_stack([cn[h], cov[ok][h]]), min_n=80, cov_names=cov_names)
        md = GW.paired_delta(beta, gene, meth_tumor, meth_normal)
        rd = GW.paired_delta(expression, gene, rna_tumor, rna_normal)
        record = {'gene': gene, 'promoter_definition': variant, 'specification': spec,
                  'eligible_probe_n': int(d['probe_counts'][variant].get(gene, 0)),
                  'rna_finite_fraction': finite_fraction, 'rna_detected_fraction': detected_fraction,
                  'rna_detection_eligible': eligible, 'methylation_rna_n': n,
                  'methylation_rna_rho': rho, 'correlation_p': p,
                  'bootstrap_n': n, 'bootstrap_reps': 200, 'bootstrap_negative_fraction': bf,
                  'bootstrap_rho_q025': blo, 'bootstrap_rho_q975': bhi,
                  'tumor_normal_methylation_delta': float(md.median()) if len(md) else np.nan,
                  'rna_tumor_normal_delta': float(rd.median()) if len(rd) else np.nan,
                  'methylation_paired_n': len(md), 'rna_paired_n': len(rd),
                  'cnv_available': gene in cnv.index,
                  'case_rule': 'RNA/methylation/CNV/specification covariates joint complete cases; baseline and joint use identical rows'}
        for prefix, fit in [('baseline_CNV', base), ('joint', joint), ('high_purity_joint', hp)]:
            record.update({prefix + '_' + key: val for key, val in fit.items()})
        # Historical phase1C tier fields explicitly refer to the new same-case
        # CNV-only model. The composition gate later uses the new joint model.
        record.update({'regression_n': base['n'], 'adjusted_methylation_coefficient': base['beta'],
                       'adjusted_methylation_p': base['p'], 'adjusted_methylation_se': base['se']})
        rows.append(record)
        if (j + 1) % 1000 == 0:
            print('NEW_PROGRESS', spec, variant, j + 1, len(genes), 'seconds', round(time.monotonic() - start), flush=True)
    result = pd.DataFrame(rows)
    for pcol, qcol in [('correlation_p', 'correlation_fdr'), ('adjusted_methylation_p', 'adjusted_methylation_fdr'),
                       ('baseline_CNV_p', 'baseline_CNV_fdr'), ('joint_p', 'joint_fdr'), ('high_purity_joint_p', 'high_purity_joint_fdr')]:
        result[qcol] = GW.bh_fdr(result[pcol])
    result['attenuation_vs_same_case_CNV_percent'] = 100 * (1 - result.joint_beta.abs() / result.baseline_CNV_beta.abs())
    destination = ANALYSIS / spec
    save(result, destination / f'associations_{variant}.tsv.gz')
    context = {'specification': spec, 'variant': variant, 'genes': len(result),
               'covariates': SPECS[spec], 'covariate_complete_n': int(low_missing_cov.sum()),
               'high_purity_covariate_complete_n': int((low_missing_cov & high).sum()),
               'finite_joint_tests': int(result.joint_p.notna().sum()),
               'finite_high_purity_joint_tests': int(result.high_purity_joint_p.notna().sum()),
               'elapsed_seconds': time.monotonic() - start,
               'bootstrap_definition': 'Original gene/window-seeded 200 resamples of raw ranks, not adjusted-coefficient bootstrap',
               'min_n': {'raw_and_baseline_CNV': 100, 'joint_and_high_purity': 80}}
    (destination / f'context_{variant}.json').write_text(json.dumps(context, indent=2), encoding='utf-8')
    print('NEW_DONE', json.dumps(context), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['prepare', 'historical_variant', 'new_variant'], required=True)
    ap.add_argument('--variant', choices=VARIANTS)
    ap.add_argument('--specification', choices=list(SPECS))
    args = ap.parse_args()
    for p in [INPUT, ANALYSIS, REVIEW]:
        p.mkdir(parents=True, exist_ok=True)
    if args.stage == 'prepare':
        prepare()
    elif args.stage == 'historical_variant':
        historical_variant(args.variant)
    elif args.stage == 'new_variant':
        new_variant(args.specification, args.variant)


if __name__ == '__main__':
    main()
