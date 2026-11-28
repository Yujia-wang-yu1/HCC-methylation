"""Run the manuscript analyses from their prepared inputs."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', choices=(
        'fixed', 'matched', 'criteria', 'selection', 'selection-check',
        'transfer', 'clinical-covariates', 'tcga-composition', 'icgc-survival',
    ))
    parser.add_argument('--inputs', type=Path, required=True, help='Prepared input directory described in README')
    parser.add_argument('--output', type=Path, required=True, help='Results directory outside the input directory')
    parser.add_argument('--selection-dir', type=Path, help='Reselection results to transfer instead of the prepared selection weights')
    parser.add_argument('--jobs', type=int, default=1, help='Workers for the selection stages')
    parser.add_argument('--rscript', default='Rscript', help='Path to Rscript for the R analyses')
    args = parser.parse_args()

    inputs = args.inputs.resolve()
    output = args.output.resolve()
    if not inputs.is_dir() or output.is_relative_to(inputs):
        parser.error('Provide an existing input directory and an output directory outside it.')
    if args.selection_dir and args.analysis != 'transfer':
        parser.error('--selection-dir applies to transfer.')
    if args.jobs < 1:
        parser.error('--jobs must be positive.')

    env = os.environ.copy()
    env.update(MISP_REPRO_INPUTS=str(inputs), MISP_REPRO_OUTPUT=str(output))
    fixed = ROOT / 'MISP_Code'
    extended = ROOT / 'MISP_Extended'
    fixed_inputs = inputs / 'MISP_Code/data'
    run_python = [sys.executable, '-B']
    name = args.analysis
    cwd = ROOT

    if name in ('fixed', 'matched', 'criteria'):
        stage = {'fixed': 'scores', 'matched': 'background', 'criteria': 'genes'}[name]
        command = run_python + [str(fixed / 'reproduce.py'), '--stage', stage,
                                '--input-dir', str(fixed_inputs), '--output-dir', str(output / name)]
    elif name in ('selection', 'selection-check'):
        stage = 'all' if name == 'selection' else 'smoke'
        command = run_python + [str(extended / 'reselection_reproduce.py'), '--stage', stage, '--jobs', str(args.jobs)]
    elif name == 'transfer':
        command = run_python + [str(extended / 'reproduce_selected_transfer.py'),
                                '--input-dir', str(inputs / 'MISP_Extended/Inputs/Selected_Transfer'),
                                '--output-dir', str(output / 'transfer')]
        if args.selection_dir:
            command += ['--selection-dir', str(args.selection_dir.resolve())]
    elif name == 'clinical-covariates':
        command = run_python + [str(extended / 'Scripts/clinical_covariate_score_models.py')]
    else:
        filename = {'tcga-composition': 'reproduce_DNA_proxies.R', 'icgc-survival': 'icgc_survival.R'}[name]
        cwd = fixed
        # Some Windows R builds use the active code page for command-line paths.
        try:
            input_arg = os.path.relpath(fixed_inputs, cwd)
            output_arg = os.path.relpath(output / name, cwd)
        except ValueError:  # Different Windows volumes.
            input_arg, output_arg = str(fixed_inputs), str(output / name)
        rscript = str(Path(args.rscript).resolve()) if Path(args.rscript).is_file() else args.rscript
        command = [rscript, filename, input_arg, output_arg]
    if subprocess.run(command, cwd=cwd, env=env).returncode != 0:
        parser.exit(1, f'{name} did not complete. See the error above.\n')


if __name__ == '__main__':
    main()
