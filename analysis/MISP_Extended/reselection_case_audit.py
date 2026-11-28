"""Restore exact phase1C CNV-baseline eligibility/inverse in existing v5 runs.

Only the CNV baseline and its BH family change. Raw correlation/bootstrap,
new joint/high-purity models and all original/source versions are preserved.
Also exports the exact same-case bit masks used by both baseline and joint.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import reselection_core as F


def main():
    data = F.load_inputs()
    log = []
    restored = []
    for spec in F.SPECS:
        for variant in F.VARIANTS:
            path = F.ANALYSIS / spec / f'associations_{variant}.tsv.gz'
            old_hash = F.sha(path)
            result = pd.read_csv(path, sep='\t', float_precision='round_trip')
            # Preserve every already-computed non-baseline statistic exactly.
            protected = [c for c in result if not c.startswith('baseline_CNV_') and c not in [
                'adjusted_methylation_coefficient', 'adjusted_methylation_p',
                'adjusted_methylation_se', 'adjusted_methylation_fdr',
                'attenuation_vs_same_case_CNV_percent', 'case_mask_hex']]
            immutable = result[protected].copy()
            beta = data['matrices'][variant]
            expression, cnv, cv = data['expression'], data['cnv'], data['cov']
            samples = sorted(set(beta.columns) & set(expression.columns) & set(cv.index) & set(cnv.columns))
            common_cov = cv.loc[samples, F.SPECS[spec]].to_numpy(float)
            cov_ok = np.isfinite(common_cov).all(axis=1)
            sample_order = pd.DataFrame({'position': np.arange(len(samples)), 'sample_id': samples,
                                         'covariate_complete': cov_ok,
                                         'high_original_ESTIMATE': cv.loc[samples, 'high_original_ESTIMATE'].to_numpy()})
            F.save(sample_order, F.ANALYSIS / spec / f'{variant}_case_sample_order.tsv')
            x = expression.reindex(index=result.gene, columns=samples).to_numpy(float)
            m = beta.reindex(index=result.gene, columns=samples).to_numpy(float)
            c = cnv.reindex(index=result.gene, columns=samples).to_numpy(float)
            before = result.baseline_CNV_p.notna().to_numpy()
            started = time.monotonic()
            fits = []
            masks = []
            for j, gene in enumerate(result.gene):
                ok = cov_ok & np.isfinite(x[j]) & np.isfinite(m[j]) & np.isfinite(c[j])
                fit = F.fast_fit(x[j, ok], m[j, ok], c[j, ok], min_n=100, require_three_unique=False)
                fits.append(fit)
                assert int(ok.sum()) == result.loc[j, 'joint_n']
                masks.append(np.packbits(ok).tobytes().hex())
                if not before[j] and np.isfinite(fit['p']):
                    restored.append({'specification': spec, 'promoter_definition': variant, 'gene': gene,
                                     'n': int(ok.sum()), 'RNA_unique': len(np.unique(x[j, ok])),
                                     'methylation_unique': len(np.unique(m[j, ok])), 'baseline_CNV_p': fit['p']})
            baseline = pd.DataFrame(fits)
            for key in baseline:
                result['baseline_CNV_' + key] = baseline[key]
            result['case_mask_hex'] = masks
            result['baseline_CNV_fdr'] = F.GW.bh_fdr(result.baseline_CNV_p)
            result['adjusted_methylation_coefficient'] = result.baseline_CNV_beta
            result['adjusted_methylation_se'] = result.baseline_CNV_se
            result['adjusted_methylation_p'] = result.baseline_CNV_p
            result['adjusted_methylation_fdr'] = result.baseline_CNV_fdr
            result['attenuation_vs_same_case_CNV_percent'] = 100 * (1 - result.joint_beta.abs() / result.baseline_CNV_beta.abs())
            assert result[protected].equals(immutable)
            F.save(result, path)
            # Same source file is re-read once to ensure stored protected values
            # round-trip; negligible text serialization error is recorded.
            reread = pd.read_csv(path, sep='\t', float_precision='round_trip')
            numerical = immutable.select_dtypes('number').columns
            max_error = float(np.nanmax(abs(reread[numerical].to_numpy() - immutable[numerical].to_numpy())))
            report = {'specification': spec, 'variant': variant, 'rows': len(result),
                      'old_sha256': old_hash, 'new_sha256': F.sha(path),
                      'baseline_finite_before': int(before.sum()), 'baseline_finite_after': int(result.baseline_CNV_p.notna().sum()),
                      'joint_finite_unchanged': int(result.joint_p.notna().sum()),
                      'protected_columns_preserved_in_memory': True, 'protected_max_serialization_error': max_error,
                      'baseline_joint_case_masks_identical': True, 'elapsed_seconds': time.monotonic() - started,
                      'baseline_rule': 'Exact phase1C n>=100/rank criterion, inverse XTX; no >=3-unique requirement',
                      'joint_rule_unchanged': 'Phase2A n>=80 and >=3 unique RNA/methylation values'}
            log.append(report)
            print('BASELINE_REFRESH', json.dumps(report), flush=True)
    F.save(pd.DataFrame(restored), F.REVIEW / 'CNV_baseline_restored_finite_tests.tsv')
    (F.REVIEW / 'same_case_CNV_baseline_correction.json').write_text(json.dumps(log, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
