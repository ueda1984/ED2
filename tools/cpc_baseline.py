#!/usr/bin/env python3
"""Reproduce the CPC audit in a new, external directory; never change solver code.

Requires GNU Fortran, Autoconf/Automake, make, BLAS/LAPACK, and test-only NumPy.
Activate the cloud toolchain first if necessary. Every failed command is retained;
known numerical defects still cause a nonzero final exit status.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tarfile
import time

import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--long', action='store_true', help='run all nine existing examples')
    parser.add_argument('--debug', action='store_true', help='also run runtime-checked fixtures')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    if out == repo or repo in out.parents:
        parser.error('--output must be outside the source checkout')
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    archive = subprocess.check_output(['git','archive',commit],cwd=repo)
    report = {'schema':1, 'source_commit':commit, 'source_archive_sha256':hashlib.sha256(archive).hexdigest(),
              'harness_sha256':sha(Path(__file__)),
              'numerical_harness_sha256':sha(repo/'tests/numerical_baseline.py'),
              'working_tree_status':subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True),
              'platform':platform.platform(), 'python':sys.version, 'numpy':np.__version__,
              'thread_settings':{k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
              'commands':[], 'reference_comparisons':[]}

    def save():
        (out/'baseline.json').write_text(json.dumps(report,indent=2)+'\n')

    def run(name, command, cwd, required=True):
        start = time.monotonic()
        with (out/(name+'.log')).open('wb') as log:
            result = subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
        item = {'name':name,'command':command,'cwd':str(cwd.relative_to(out)),
                'exit_code':result.returncode,'seconds':time.monotonic()-start,'required':required}
        report['commands'].append(item)
        save()
        print(f'{name}: exit {result.returncode}',flush=True)
        return result.returncode == 0

    def extract(name):
        stage = out/name
        stage.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
            tar.extractall(stage,filter='data')
        return stage

    for tool in ['gfortran','autoconf','automake','make']:
        run(tool+'-version',[tool,'--version'],out)
    run('cpu',['lscpu'],out,required=False)
    stage = extract('release')
    if not run('autoreconf',['autoreconf','-vfi'],stage):
        return 1
    if not run('configure-default',['./configure','FC=gfortran'],stage):
        return 1
    run('build-default',['make','-j2'],stage,required=False)
    flags = '-O2 -fopenmp -fallow-argument-mismatch'
    report['baseline_flags'] = flags
    if not run('configure-baseline',['./configure','FC=gfortran','FCFLAGS='+flags],stage):
        return 1
    if not run('clean',['make','clean'],stage) or not run('build-baseline',['make','-j2'],stage):
        return 1
    exe = stage/'source/QS3ED2'
    report['executable_sha256'] = sha(exe)
    run('linked-libraries',['ldd',str(exe)],out)
    run('check',['make','check'],stage)
    # Verify that the Automake harness actually ran one test.
    text = (stage/'tests/test-suite.log').read_text()
    counts = dict((k,int(v)) for k,v in re.findall(r'^# (TOTAL|PASS|FAIL|ERROR|SKIP):\s*(\d+)',text,re.M))
    report['automake_counts'] = counts
    report['automake_count_check'] = counts == {'TOTAL':1,'PASS':1,'FAIL':0,'ERROR':0,'SKIP':0}
    run('numerical-release',[sys.executable,str(repo/'tests/numerical_baseline.py'),
                            '--exe',str(exe),'--output',str(out/'numerical-release')],out)
    if args.long:
        long_ok = run('check-long',['make','check-long'],stage)
        log = (out/'check-long.log').read_text()
        report['long_passes'] = re.findall(r'^PASS: (test_ex\d.sh)$',log,re.M)
        report['long_count_check'] = long_ok and len(report['long_passes']) == 9
        for refdir in sorted((stage/'examples/ref_data').iterdir()):
            actualdir = stage/'examples'/refdir.name
            item = {'example':refdir.name,'comparisons':[], 'input_sha256':{}}
            for path in sorted(actualdir.glob('*.dat')):
                if path.name.startswith('list_') or path.name in ['input.dat','input_param.dat']:
                    item['input_sha256'][path.name] = sha(path)
            stdout = actualdir/'output.dat'
            if stdout.exists():
                item['dimensions'] = {k:int(v) for k,v in re.findall(r'\b(THS(?:\(k\))?)\s*=\s*(\d+)',stdout.read_text())}
            for filename,indices,atol,rtol in [('eigenvalues.dat',1,1e-9,1e-9),
                                             ('local_mag.dat',1,1e-6,1e-6),
                                             ('two_body_cf_xyz.dat',2,1e-6,1e-6)]:
                comparison = {'file':filename,'absolute_tolerance':atol,'relative_tolerance':rtol}
                try:
                    a = np.loadtxt(actualdir/filename,ndmin=2)
                    b = np.loadtxt(refdir/filename,ndmin=2)
                    assert a.shape == b.shape, 'shape differs'
                    assert np.isfinite(a).all() and np.isfinite(b).all(), 'nonfinite values'
                    np.testing.assert_array_equal(a[:,:indices],b[:,:indices])
                    np.testing.assert_allclose(a[:,indices:],b[:,indices:],atol=atol,rtol=rtol)
                    comparison.update(status='passed',max_absolute_error=float(np.max(np.abs(a[:,indices:]-b[:,indices:]))))
                except (OSError,ValueError,AssertionError) as error:
                    comparison.update(status='failed',error=str(error))
                item['comparisons'].append(comparison)
            report['reference_comparisons'].append(item)
    if args.debug:
        debug = extract('debug')
        debug_flags = '-Og -g -fopenmp -fallow-argument-mismatch -fcheck=all -fbacktrace -Wall -Wextra'
        report['debug_flags'] = debug_flags
        ok = run('debug-autoreconf',['autoreconf','-vfi'],debug)
        ok = ok and run('debug-configure',['./configure','FC=gfortran','FCFLAGS='+debug_flags],debug)
        ok = ok and run('debug-build',['make','-j2'],debug)
        if ok:
            run('numerical-debug',[sys.executable,str(repo/'tests/numerical_baseline.py'),
                                  '--exe',str(debug/'source/QS3ED2'),'--output',str(out/'numerical-debug')],out)
    save()
    failed = any(c['required'] and c['exit_code'] != 0 for c in report['commands'])
    failed |= not report['automake_count_check']
    failed |= args.long and not report.get('long_count_check',False)
    failed |= any(c['status'] != 'passed' for e in report['reference_comparisons'] for c in e['comparisons'])
    return int(failed)


if __name__ == '__main__':
    sys.exit(main())
