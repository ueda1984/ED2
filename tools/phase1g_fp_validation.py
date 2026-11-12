#!/usr/bin/env python3
"""Recheck telemetry expectations in the actual arithmetic mode, without setters.

Only the test driver/harness is refreshed; numerical production files stay fixed.
Initial full-run/Ofast failures remain in their original command manifests/logs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--default-runs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    repo=Path(__file__).resolve().parents[1]
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for key in ['ED2_RNG_SEED','ED2_INITIAL_VECTOR','ED2_DIAGNOSTICS','ED2_RESIDUAL_ATOL','ED2_RESIDUAL_RTOL','ED2_RESIDUAL_SCALE','ED2_ENERGY_ORIGIN']:env.pop(key,None)
    manifest=dict(harness_sha256=sha(Path(__file__)),test_source_sha256={name:sha(repo/'tests'/name) for name in ['krylov_diagnostics_driver.F90','diagnostics_reproducibility.py']},commands=[])
    for profile in ['o2','gnu13','bounds','default']:
        pr=(a.default_runs if profile=='default' else a.runs).resolve()/profile
        dest=out/profile;dest.mkdir();stage=pr/'after';before=pr/'before'
        for name in manifest['test_source_sha256']:
            shutil.copy2(repo/'tests'/name,stage/'tests'/name)
        flags=['-O2','-fopenmp'] if profile in ['o2','gnu13'] else ['-O2','-g','-fopenmp','-fcheck=all','-fbacktrace','-Wall','-Wextra'] if profile=='bounds' else ['-fopenmp','-Ofast']
        fc='/workspace/ed2-env/gnu13/usr/bin/gfortran-13' if profile=='gnu13' else 'gfortran'
        binary=dest/'before-driver'
        commands=[['make','-C',str(stage/'tests'),'krylov_diagnostics_driver'],
         [fc,*flags,'-cpp','-DPHASE1F','-I'+str(before/'source'),str(repo/'tests/krylov_diagnostics_driver.F90'),
          str(before/'source/krylov_solvers.o'),str(before/'source/krylov_observations.o'),str(before/'source/blas_lapack_interfaces.o'),
          '-L/workspace/ed2-env/tools/usr/lib/x86_64-linux-gnu/openblas-pthread','-L/workspace/ed2-env/tools/usr/lib/gcc/x86_64-linux-gnu/14',
          '-llapack','-lopenblas','-o',str(binary)],
         [sys.executable,'-B',str(repo/'tests/diagnostics_reproducibility.py'),'--driver',str(stage/'tests/krylov_diagnostics_driver'),
          '--before-driver',str(binary),'--exe',str(stage/'source/QS3ED2'),'--output',str(dest/'results')]]
        for i,cmd in enumerate(commands):
            log=dest/f'{i}.log'
            with log.open('w') as stream:r=subprocess.run(cmd,cwd=dest,env=env,stdout=stream,stderr=subprocess.STDOUT)
            manifest['commands'].append(dict(profile=profile,args=cmd,exit_code=r.returncode,log_sha256=sha(log)))
            (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            print(profile,i,r.returncode,flush=True)
            if r.returncode:return r.returncode
    return 0


if __name__=='__main__':raise SystemExit(main())
