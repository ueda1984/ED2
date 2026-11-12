#!/usr/bin/env python3
"""Supplement the retained fixed-budget ALG=2 timing with a converged workload."""
import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tests'))
import numpy as np
from diagnostics_reproducibility import run,same,check_observation,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);root=a.runs.resolve()
    n=192;rng=np.random.default_rng(173);u,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    e=np.linspace(-1.,2.,n);h=(u*e)@u.conj().T
    item=dict(n=n,h=h,e=e,u=u,alg=2,tol=1e-10,seed=19,wanted=2,keep=6,capacity=32,budget=400)
    rows=[]
    for profile in ['o2','gnu13','bounds']:
        binary=root/profile/'after/tests/krylov_diagnostics_driver';old=root/profile/'before-driver'
        for repeat in range(7):
            modes=[('before',old,False),('off',binary,False),('on',binary,True)]
            if repeat%2:modes.reverse()
            results={}
            for mode,driver,enabled in modes:
                result=run(driver,out,f'{profile}-r{repeat}-{mode}',item,1,enabled)
                assert result['status']==0,(profile,mode,result)
                values=np.loadtxt(result['work']/'eigenpairs.dat',ndmin=2)[:,1]
                np.testing.assert_allclose(values,e[:2],atol=1e-9,rtol=1e-9)
                row={k:v for k,v in result.items() if k!='work'}
                row.update(profile=profile,repeat=repeat,mode=mode,driver_sha256=sha(driver))
                rows.append(row);results[mode]=result
            same(results['before'],results['off']);same(results['off'],results['on'])
            check_observation(results['on'],item)
    (out/'report.json').write_text(json.dumps(dict(schema=1,harness_sha256=sha(Path(__file__)),
                                                fixture={k:v for k,v in item.items() if k not in ['h','e','u']},
                                                generator='complex Gaussian QR seed173; n192; spectrum linspace(-1,2)',
                                                rows=rows),indent=2,allow_nan=False)+'\n')
    print('63 converged ALG=2 timing runs; before/off/on output bytes agree.')


if __name__=='__main__':main()
