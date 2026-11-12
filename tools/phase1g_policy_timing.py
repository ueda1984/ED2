#!/usr/bin/env python3
"""Matched-start repeated performance/accuracy comparisons, including failures."""
import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tests'))
import numpy as np
from residual_policy_validation import run,reduced,same,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--driver',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);driver=a.driver.resolve(strict=True)
    n=192;rng=np.random.default_rng(173);u,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    e=np.linspace(-1.,2.,n);h=(u*e)@u.conj().T;rows=[]
    for alg in [1,2]:
        t=dict(family='timing-separated',n=n,h=h,e=e,u=u,alg=alg,tol=1e-10,seed=19,wanted=1 if alg==1 else 2,
               keep=6,capacity=32,budget=n if alg==1 else 400)
        for threads in [1,2]:
            for repeat in range(5):
                modes=[('legacy',None),('loose',dict(residual_atol=1e30)),('strict',dict(residual_atol=1e-10))]
                if repeat%2:modes.reverse()
                results={}
                for label,policy in modes:
                    r=run(driver,out,f'a{alg}-t{threads}-r{repeat}-{label}',t,threads,policy)
                    results[label]=r
                    rows.append(dict(alg=alg,threads=threads,repeat=repeat,mode=label,**reduced(r)))
                same(results['legacy'],results['loose'])
                # Do not assume the stronger policy always converges. If it does,
                # independently compare achieved errors and retain all metrics.
                for r in results.values():
                    if r['status']==0:assert max(r['assessment']['lowest_eigenvalue_errors'])<1e-8
    (out/'report.json').write_text(json.dumps(dict(schema=1,n=n,seed=19,energy_tolerance=1e-10,
                                                wanted_by_alg=[1,2],keep=6,capacity=32,alg2_budget=400,
                                                driver_sha256=sha(driver),harness_sha256=sha(Path(__file__)),rows=rows),
                                            indent=2,allow_nan=False)+'\n')
    print('60 paired kernel timing/accuracy runs completed; all statuses retained.')


if __name__=='__main__':main()
