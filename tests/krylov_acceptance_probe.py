#!/usr/bin/env python3
"""Seeded actual-CLI observable stress at unchanged energy-change policies.

Requires the Phase 1-C instrumented binary from krylov_benchmarks.build; reports
all acceptance mismatches, without changing the existing sparse-field harness.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
import numpy as np
from sparse_field_validation import prepare,reference,FIXTURES
from krylov_benchmarks import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--exe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();exe=a.exe.resolve(strict=True);out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records=FIXTURES['single-site4'];h,e,lm,cf=reference(records);results=[]
    for alg in [1,2]:
        for tol in [1e-14,1e-18]:
            for threads in [1,2]:
                for seed in [1000003*s for s in range(1,129)]:
                    name=f'alg{alg}-tol{tol}-omp{threads}-seed{seed}';work=out/name
                    prepare(work,records,alg,-1)
                    text=(work/'input.dat').read_text().replace('wr_wf=0','wr_wf=1').replace('lnc_ene_conv=1d-14',f'lnc_ene_conv={tol}').replace('lnc_ene_conv=1d-13',f'lnc_ene_conv={tol}')
                    (work/'input.dat').write_text(text)
                    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',ED2_BENCH_SEED=str(seed))
                    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:r=subprocess.run([str(exe)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=60)
                    row=dict(id=name,alg=alg,tol=tol,threads=threads,seed=seed,exit_code=r.returncode,input_sha256={p.name:sha(p) for p in work.glob('*.dat') if p.name in ['input.dat','bonds.dat','fields.dat','spins.dat','caps.dat']},rng_state=(work/'benchmark-rng.dat').read_text().splitlines())
                    if r.returncode:row['classification']='solver_failure'
                    else:
                        values=np.loadtxt(work/'eigenvalues.dat',ndmin=2)[:,1];basis=np.loadtxt(work/'benchmark-basis.dat',dtype=int,ndmin=2)
                        idx=[sum(1<<(4-int(s)) for s in b if s!=0) for b in basis];vectors=[]
                        for j in range(len(values)):
                            b=(work/f'wf{j+1}.bin').read_bytes();assert b[:4]==b[-4:] and len(b)==16*16+8
                            q=np.zeros(16,complex);q[idx]=np.frombuffer(b[4:-4],dtype=np.complex128);vectors.append(q)
                        q=np.array(vectors).T;res=np.linalg.norm(h@q-q*values,axis=0)
                        actual_lm=np.loadtxt(work/'local_mag.dat',ndmin=2)[:,1:];actual_cf=np.loadtxt(work/'two_body_cf_xyz.dat',ndmin=2)[:,2:]
                        ok_lm=bool(np.allclose(actual_lm,lm,atol=1e-9,rtol=1e-9));ok_cf=bool(np.allclose(actual_cf,cf,atol=1e-9,rtol=1e-9))
                        row.update(classification='observable_agreement' if ok_lm and ok_cf else 'accepted_observable_failure',eigenvalue_errors=abs(values-e[:len(values)]).tolist(),residual_norms=res.tolist(),max_magnetization_error=float(np.max(abs(actual_lm-lm))),max_correlation_error=float(np.max(abs(actual_cf-cf))),orthogonality_error=float(np.linalg.norm(q.conj().T@q-np.eye(len(values)),2)))
                    results.append(row)
    counts={k:sum(r['classification']==k for r in results) for k in sorted(set(r['classification'] for r in results))}
    (out/'report.json').write_text(json.dumps(dict(executable_sha256=sha(exe),harness_sha256=sha(Path(__file__)),counts=counts,cases=results),indent=2,allow_nan=False)+'\n');print(counts)
    return 0


if __name__=='__main__':raise SystemExit(main())
