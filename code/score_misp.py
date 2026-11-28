"""Apply the fixed 30-gene MISP weights from Supplementary Table S3."""

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd


RESULTS = Path(__file__).resolve().parents[1] / 'data' / 'Results.xlsx'


def load_weights():
    table = pd.read_excel(RESULTS, sheet_name='S3', header=1, nrows=30,
                          usecols=['gene', 'MISP weight'])
    weights = table.set_index('gene')['MISP weight'].astype(float)
    if len(weights) != 30 or weights.index.hasnans or not weights.index.is_unique or not (
        np.isfinite(weights).all() and weights.gt(0).all()
    ):
        raise ValueError('Table S3 must contain exactly 30 distinct genes and finite positive weights.')
    return weights


def read_matrix(path):
    with path.open(encoding='utf-8-sig', newline='') as source:
        header = next(csv.reader(source, delimiter='\t'), [])
    samples = header[1:]
    if not samples or any(not name.strip() for name in samples) or len(set(samples)) != len(samples):
        raise ValueError('The TSV header must contain nonempty, unique sample IDs after the gene column.')
    frame = pd.read_csv(path, sep='\t', dtype={header[0]: 'string'}, encoding='utf-8-sig')
    symbols = frame.iloc[:, 0].str.strip()
    if symbols.isna().any() or symbols.eq('').any() or symbols.duplicated().any():
        raise ValueError('Gene symbols must be nonempty and unique.')
    values = frame.iloc[:, 1:].apply(pd.to_numeric, errors='raise')
    values.index = symbols
    return values


def score_misp(frame, weights, assay, standardized=False):
    """Return scores using gene-wise population SD and the full weight denominator."""
    if not frame.index.isin(weights.index).any():
        raise ValueError('The input contains none of the fixed genes in Table S3.')
    values = frame.reindex(weights.index).astype(float)
    values = values.where(np.isfinite(values))
    if not standardized:
        means = values.mean(axis=1)
        deviations = values.std(axis=1, ddof=0).replace(0, np.nan)
        values = values.sub(means, axis=0).div(deviations, axis=0)
    weighted = values.mul(weights, axis=0).sum(axis=0, min_count=1) / weights.sum()
    if assay == 'rna':
        return (-weighted).rename('MISP_RNA_suppression_score')
    return weighted.rename('MISP_methylation_score')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path, help='TSV with gene symbols in rows and sample IDs in columns')
    parser.add_argument('--assay', required=True, choices=('methylation', 'rna'))
    parser.add_argument('--standardized', action='store_true', help='Input already contains cohort-calibrated gene z scores')
    parser.add_argument('--output', required=True, type=Path, help='New output TSV path')
    args = parser.parse_args()
    try:
        scores = score_misp(read_matrix(args.input), load_weights(), args.assay, args.standardized)
        scores.to_csv(args.output, sep='\t', index_label='sample_id', na_rep='NA', mode='x')
    except (OSError, ValueError, TypeError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
