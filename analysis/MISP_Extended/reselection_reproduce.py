"""Self-contained selection-layer runner using supplied full processed inputs.

No downloads, installations, external source directory or benchmark-as-result
fallback. Run --stage smoke first; --stage all executes the complete workflow.
"""
from pathlib import Path
import argparse
import concurrent.futures
import json
import subprocess
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reselection_core as F


def run_job(command):
    subprocess.run([sys.executable, '-B'] + command, cwd=HERE, check=True)


def run_jobs(jobs, workers):
    if workers == 1:
        for job in jobs:
            run_job(job)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            # Each worker launches an independent Python process. All awaited;
            # input arrays are read-only and each job has a distinct output.
            for _ in pool.map(run_job, jobs):
                pass


def smoke():
    data = F.load_inputs()
    target_path = F.PROJECT / 'smoke_model_targets.tsv'
    targets = pd.read_csv(target_path, sep='\t', float_precision='round_trip')
    rows = []
    for r in targets.itertuples():
        beta = data['matrices'][r.promoter_definition]
        exp, cnv, cv = data['expression'], data['cnv'], data['cov']
        samples = sorted(set(beta.columns) & set(exp.columns) & set(cnv.columns) & set(cv.index))
        y = exp.loc[r.gene, samples].to_numpy(float)
        m = beta.loc[r.gene, samples].to_numpy(float)
        c = cnv.loc[r.gene, samples].to_numpy(float)
        nuisance = cv.loc[samples, F.SPECS[r.specification]].to_numpy(float)
        eligible_cohort = np.isfinite(nuisance).all(axis=1)
        ok = eligible_cohort & np.isfinite(y) & np.isfinite(m) & np.isfinite(c)
        yy, mm, cc = y[ok], m[ok], c[ok]
        base = F.fast_fit(yy, mm, cc, min_n=100, require_three_unique=False)
        joint = F.fast_fit(yy, mm, np.column_stack([cc, nuisance[ok]]), min_n=80, cov_names=['CNV'] + F.SPECS[r.specification])
        high = cv.loc[samples, 'high_original_ESTIMATE'].to_numpy(bool)[ok]
        hp = F.fast_fit(yy[high], mm[high], np.column_stack([cc[high], nuisance[ok][high]]), min_n=80, cov_names=['CNV'] + F.SPECS[r.specification])
        values = {}
        for prefix, fit in [('baseline_CNV', base), ('joint', joint), ('high_purity_joint', hp)]:
            for field in ['n', 'beta', 'p']:
                values[prefix + '_' + field] = fit[field]
        finite_fraction = np.isfinite(y[eligible_cohort]).mean()
        detected = (np.nan_to_num(y[eligible_cohort], nan=0) > 0).mean()
        bf = np.nan
        if finite_fraction >= .8 and detected >= .2 and len(yy) >= 100:
            bf = F.raw_bootstrap(mm, yy, r.gene, r.promoter_definition)[0]
        values['bootstrap_negative_fraction'] = bf
        for field, actual in values.items():
            expected = getattr(r, field)
            if pd.isna(actual) or pd.isna(expected):
                passed = bool(pd.isna(actual) and pd.isna(expected))
                error = np.nan
            else:
                error = float(abs(actual - expected))
                tolerance = 0 if field.endswith('_n') or field == 'bootstrap_negative_fraction' else 1e-10
                passed = error <= tolerance
            rows.append({'specification': r.specification, 'promoter_definition': r.promoter_definition, 'gene': r.gene,
                         'field': field, 'computed': actual, 'expected_QA_target': expected, 'abs_error': error, 'passed': passed})
    fixed = pd.read_csv(F.PROJECT / 'phase2A_output/locked_immune_methylation_program_genes.tsv', sep='\t').set_index('gene').MISP_weight.sort_index()
    frozen = pd.read_csv(F.PROJECT / 'fixed_primary_TCGA_scores.tsv', sep='\t', float_precision='round_trip')
    weights = fixed.to_numpy(float, copy=True)
    weights /= weights.sum()
    mz = F.LOCK.zrows(data['matrices']['TSS200'].loc[fixed.index, frozen.sample_id]).to_numpy(float)
    rz = F.LOCK.zrows(data['expression'].loc[fixed.index, frozen.sample_id]).to_numpy(float)
    for name, computed, expected in [
        ('original_fixed30_methylation_score', np.nansum(mz * weights[:, None], axis=0), frozen.MISP_methylation_score),
        ('original_fixed30_RNA_score', -np.nansum(rz * weights[:, None], axis=0), frozen.MISP_RNA_suppression_score)]:
        error = float(np.max(abs(computed - expected)))
        rows.append({'field': name, 'abs_error': error, 'passed': error < 1e-12})
    checks = pd.DataFrame(rows)
    F.save(checks, F.REVIEW / 'smoke_checks.tsv')
    summary = {'stage': 'smoke', 'passed': bool(checks.passed.all()), 'numerical_checks': len(checks),
               'model_gene_window_specification_rows': len(targets), 'unique_genes': int(targets.gene.nunique()),
               'all_original_full_input_matrices_loaded': True,
               'QA_targets_used_only_for_comparison': True,
               'scope': 'Three specifications, three windows, exact-case CNV/joint/high models, original200 raw-rank bootstrap and fixed30 scoring; full selection requires full stages'}
    (F.REVIEW / 'smoke_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    if not summary['passed']:
        raise AssertionError('Smoke numerical comparison failed; see smoke_checks.tsv')


def historical_finalize():
    run_job(['reselection_finalize.py', '--stage', 'historical_models'])
    run_job(['reselection_finalize.py', '--stage', 'pathways_selection', '--specification', 'historical_baseline'])
    run_job(['reselection_post.py', '--stage', 'all', '--specification', 'historical_baseline'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage', choices=['smoke', 'historical', 'historical_finalize', 'models', 'pathways', 'robustness', 'case_audit', 'all'], default='smoke')
    ap.add_argument('--specification', choices=['historical_baseline'] + list(F.SPECS), default='joint_DNA_ABSOLUTE')
    ap.add_argument('--variant', choices=F.VARIANTS + ['all'], default='all')
    ap.add_argument('--jobs', type=int, default=1)
    args = ap.parse_args()
    if args.jobs < 1:
        ap.error('--jobs must be positive')
    for directory in [F.ANALYSIS, F.REVIEW]:
        directory.mkdir(parents=True, exist_ok=True)
    for name in ['full_gene_inputs.pkl', 'frozen_gene_sets.pkl']:
        if not (F.INPUT / name).exists():
            raise FileNotFoundError(f'Missing prepared input: {F.INPUT / name}')
    variants = F.VARIANTS if args.variant == 'all' else [args.variant]
    if args.stage == 'smoke':
        smoke()
    elif args.stage == 'historical':
        run_jobs([['reselection_core.py', '--stage', 'historical_variant', '--variant', v] for v in variants], args.jobs)
    elif args.stage == 'historical_finalize':
        historical_finalize()
    elif args.stage == 'models':
        if args.specification == 'historical_baseline':
            ap.error('Use --stage historical for the historical baseline')
        run_jobs([['reselection_core.py', '--stage', 'new_variant', '--specification', args.specification, '--variant', v] for v in variants], args.jobs)
    elif args.stage == 'pathways':
        run_job(['reselection_finalize.py', '--stage', 'pathways_selection', '--specification', args.specification])
    elif args.stage == 'robustness':
        run_job(['reselection_post.py', '--stage', 'all', '--specification', args.specification])
    elif args.stage == 'case_audit':
        run_job(['reselection_case_audit.py'])
    else:
        smoke()
        run_jobs([['reselection_core.py', '--stage', 'historical_variant', '--variant', v] for v in F.VARIANTS], args.jobs)
        historical_finalize()
        run_jobs([['reselection_core.py', '--stage', 'new_variant', '--specification', spec, '--variant', v] for spec in F.SPECS for v in F.VARIANTS], args.jobs)
        run_job(['reselection_case_audit.py'])
        run_jobs([['reselection_finalize.py', '--stage', 'pathways_selection', '--specification', spec] for spec in F.SPECS], args.jobs)
        run_jobs([['reselection_post.py', '--stage', 'all', '--specification', spec] for spec in F.SPECS], args.jobs)


if __name__ == '__main__':
    main()
