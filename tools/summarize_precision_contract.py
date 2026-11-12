#!/usr/bin/env python3
"""Reduce Phase 1-E raw runs into reviewable, reproducible evidence (no vectors).

Run with activated compilers/libraries, --runs /workspace --output NEW_DIRECTORY.
The raw layout is phase1e-{o2,gnu13,bounds,default}-results/report.json,
phase1e-{o2,gnu13}-acceptance/report.json and the regression paths documented in
docs/development/phase1e-validation.md. Never overwrite evidence in place.
"""
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys

BASE = 'b56fe21678e62b8b6216fe0ec742d5dbd62a8df6'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(',', ':')) if isinstance(v, (list, dict)) else v
                             for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    root = args.runs.resolve(); repo = Path(__file__).resolve().parents[1]
    summary = {}; selected = {}
    for profile in ['o2', 'gnu13', 'bounds', 'default']:
        report_path = root/f'phase1e-{profile}-results/report.json'
        report = json.loads(report_path.read_text())
        rows = report['cases']
        reduced = []
        for row in rows:
            row = dict(row)
            assessment = row.pop('candidate_assessment', None)
            row.update({f'assessment_{k}': v for k, v in (assessment or {}).items()})
            reduced.append(row)
        write_csv(out/f'{profile}-cases.csv', reduced)
        info = {k: v for k, v in report.items() if k not in ['cases', 'analytic_controls']}
        info.update(raw_report_sha256=sha(report_path), total=len(rows))
        info['endpoint_comparison'] = {}
        for flag in ['absolute_residual_le_1e_minus10', 'relative_residual_le_1e_minus10',
                     'centered_residual_le_1e_minus10', 'reference_subspace_le_1e_minus5']:
            candidates = [r for r in rows if r.get('candidate_assessment')]
            info['endpoint_comparison'][flag] = dict(
                all_available_candidates=sum(r['candidate_assessment'][flag] for r in candidates),
                successful_candidates=sum(r['candidate_assessment'][flag] for r in candidates if r['legacy_success']))
        summary[profile] = info
        selected[profile] = [r for r in rows if r['id'] in
                             ['0017-omp1', '0018-omp1', '0040-omp2', '0042-omp1', '0043-omp1']]
    write_json(out/'analytic-controls.json', report['analytic_controls'])
    write_json(out/'selected-cases.json', selected)
    write_json(out/'summary.json', summary)
    acceptance = {}
    for profile in ['o2', 'gnu13']:
        path = root/f'phase1e-{profile}-acceptance/report.json'
        report = json.loads(path.read_text())
        write_csv(out/f'{profile}-acceptance.csv', report['cases'])
        acceptance[profile] = {k: v for k, v in report.items() if k != 'cases'}
        acceptance[profile].update(raw_report_sha256=sha(path), total=len(report['cases']))
        failed = [row for row in report['cases'] if row['classification'] == 'solver_failure']
        reasons = [(path.parent/row['id']/'stdout.log').read_text() for row in failed]
        acceptance[profile].update(
            all_failures_are_alg2_tol_1e18_budget_exhaustion=all(row['alg'] == 2 and row['tol'] == 1e-18 and
                                                               'maximum restart iterations without convergence' in text
                                                               for row, text in zip(failed, reasons)),
            failure_eigenvalue_file_count=sum((path.parent/row['id']/'eigenvalues.dat').exists() for row in failed),
            selected_failure=dict(id=failed[0]['id'], stdout_tail=reasons[0][-1600:]) if failed else None)
    write_json(out/'acceptance-summary.json', acceptance)
    regressions = {}
    for profile in ['o2', 'gnu13', 'bounds']:
        regressions[profile] = {}
        for suite in ['stability', 'stability-cli', 'independent', 'fields']:
            path = root/f'phase1e-{profile}-{suite}/report.json'
            if not path.exists():
                regressions[profile][suite] = {'status': 'unrun'}
                continue
            report = json.loads(path.read_text())
            regressions[profile][suite] = {k: v for k, v in report.items()
                                          if k in ['passed', 'failed', 'cli_counts', 'operator_counts', 'permutation_counts']}
            regressions[profile][suite]['report_sha256'] = sha(path)
        for suite in ['lapack', 'standard']:
            path = root/f'phase1e-{profile}-{suite}.log'
            text = path.read_text()
            if suite == 'standard':
                if '# TOTAL: 1' not in text or '# PASS:  1' not in text or '# FAIL:  0' not in text:
                    raise ValueError('unexpected standard test result')
            regressions[profile][suite] = {'log_sha256': sha(path), 'log': text}
    default_log = root/'phase1e-default-standard.log'
    regressions['default'] = {'standard': {'log_sha256': sha(default_log), 'log': default_log.read_text()}}
    write_json(out/'regressions.json', regressions)
    original = repo/'docs/development/phase1d-krylov/gnu13-regression-failures.json'
    instrumentation = {}
    for profile in ['o2', 'gnu13']:
        probe = root/f'phase1e-{profile}-probe'
        instrumentation[profile] = dict(
            provenance=json.loads((probe/'provenance.json').read_text()),
            patch_lines=(probe/'phase1c-instrumentation.patch').read_text().splitlines(),
            build_log=(probe/'phase1c-build.log').read_text())
    write_json(out/'probe-instrumentation.json', instrumentation)
    metadata = dict(numerical_source_base=BASE, python=sys.version, platform=platform.platform(),
                    cpu=next(line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),
                    source_sha256={p.name: sha(p) for p in sorted((repo/'source').glob('*.f90'))},
                    validation_sha256={str(p.relative_to(repo)): sha(p) for p in
                                       [repo/'tests/precision_contract.py', repo/'tests/krylov_stress_driver.f90',
                                        repo/'tests/krylov_acceptance_probe.py', repo/'tests/krylov_benchmark_probe.f90', Path(__file__)]},
                    original_gnu13_failure_sha256=sha(original),
                    original_gnu13_failure=json.loads(original.read_text()),
                    thread_environment=dict(OMP_NUM_THREADS=[1, 2], OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'),
                    profiles={})
    for name, digest in metadata['source_sha256'].items():
        original_bytes = subprocess.check_output(['git', 'show', BASE+':source/'+name], cwd=repo)
        if hashlib.sha256(original_bytes).hexdigest() != digest:
            raise ValueError('production source differs from starting commit: '+name)
    for profile in ['o2', 'gnu13', 'bounds', 'default']:
        build = root/f'ed2-phase1e-{profile}'
        fc = '/workspace/ed2-env/gnu13/usr/bin/gfortran-13' if profile == 'gnu13' else 'gfortran'
        flags = {'o2': '-O2 -fopenmp', 'gnu13': '-O2 -fopenmp',
                 'bounds': '-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra',
                 'default': '-fopenmp -Ofast'}[profile]
        metadata['profiles'][profile] = dict(compiler=subprocess.check_output([fc, '--version'], text=True).splitlines()[0],
                                             flags=flags, driver_sha256=sha(build/'tests/krylov_stress_driver'),
                                             executable_sha256=sha(build/'source/QS3ED2'),
                                             linked_libraries=subprocess.check_output(['ldd', str(build/'source/QS3ED2')], text=True))
        # Verify every production Fortran byte in the uninstrumented build.
        for name, digest in metadata['source_sha256'].items():
            if sha(build/'source'/name) != digest:
                raise ValueError('production source changed: '+str(build/'source'/name))
    replay = root/'phase1e-o2-replay/report.json'
    normal = root/'phase1e-o2-results/report.json'
    if json.loads(replay.read_text()) != json.loads(normal.read_text()):
        raise ValueError('fixed-state replay differs')
    metadata['replay'] = dict(records=136, identical=True, excluded_fields='timings (not collected in contract report)')
    write_json(out/'provenance.json', metadata)
    print(json.dumps({k: v['counts'] for k, v in summary.items()}))


if __name__ == '__main__':
    main()
