#!/usr/bin/env python3
"""Characterize unchanged Krylov policies; numerical findings remain data, not xfails.

Exit 1 means a driver/harness/provenance error. An accepted pair missing the ground
space or a solver nonconvergence remains explicitly counted in report.json.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import tempfile
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixtures():
    for n in [32,96]:
        spectra={
            'separated':np.linspace(-1,2,n),
            'degenerate':np.r_[[-1.]*4,np.linspace(-.5,2,n-4)],
            'clustered':np.r_[[-1.,-1.+1e-10,-1.+2e-10,-1.+3e-10],np.linspace(-.5,2,n-4)],
            'zero-cluster':np.r_[[0.,1e-10,2e-10,3e-10],np.linspace(.5,2,n-4)],
            'zero-centered':np.r_[[-1e-10,0.,1e-10],np.linspace(.5,2,n-3)]}
        rng=np.random.default_rng(8191+n)
        u,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        for family,values in spectra.items():
            h=(u*values)@u.conj().T
            yield family,n,h,values,u


def definitions():
    for family,n,h,e,u in fixtures():
        for alg in [1,2]:
            for tol in [1e-6,1e-12,1e-20]:
                for seed in [19,73,211]:
                    yield dict(family=family,n=n,h=h,e=e,u=u,alg=alg,tol=tol,seed=seed,
                               wanted=1 if alg==1 else 3,keep=6,capacity=12,budget=n if alg==1 else 120)
    family,n,h,e,u=next(x for x in fixtures() if x[0]=='clustered' and x[1]==96)
    for wanted in [1,3,5]:
        for cap in [12,24,48]:
            for budget in [1,5,120]:
                yield dict(family='parameter-grid',n=n,h=h,e=e,u=u,alg=2,tol=1e-12,seed=19,
                           wanted=wanted,keep=max(6,2*wanted),capacity=cap,budget=budget)
    # Exact analytic adversarial starts isolate inclusion from residual validity.
    e=np.array([-1.,-1.+1e-10,0.,2.]);u=np.eye(4,dtype=complex);h=np.diag(e)
    for family,start in [('excluded-ground',np.array([0,1,0,0],complex)),
                         ('weak-ground-overlap',np.array([1e-8,1,0,0],complex))]:
        yield dict(family=family,n=4,h=h,e=e,u=u,alg=1,tol=1e-12,seed=19,wanted=1,keep=2,
                   capacity=4,budget=4,start=start)
    # Nearly rank-deficient ground content under an ordinary random start is also
    # explored by the seeded clustered matrices above; no best seed is selected.


def run(driver,root,t,threads,index):
    work=root/f'{index:04d}-omp{threads}';work.mkdir()
    n=t['n'];h=t['h'];wanted=t['wanted']
    start=t.get('start',np.zeros(n,complex))
    header=f"&problem n={n},alg={t['alg']},wanted={wanted},keep={t['keep']},capacity={t['capacity']},"
    header+=f"maxiter={t['budget']},miniter=1,interval=1,tol={t['tol']},seed={t['seed']},"
    header+=f"use_start={'.true.' if 'start' in t else '.false.'}/\n"
    deck=header+''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in h.flatten())
    deck+=''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in start)
    (work/'input.dat').write_text(deck)
    r={k:v for k,v in t.items() if k not in ['h','e','u','start']}
    r.update(id=work.name,threads=threads,input_sha256=digest(work/'input.dat'),matrix_sha256=hashlib.sha256(h.tobytes()).hexdigest())
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    try:
        wall=time.perf_counter()
        with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
            result=subprocess.run([str(driver)],cwd=work,stdin=inp,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=90)
        r.update(exit_code=result.returncode,wall_seconds=time.perf_counter()-wall)
        lines=(work/'status.dat').read_text().splitlines();status=int(lines[0]);text=(work/'stdout.log').read_text()
        r.update(solver_status=status,reason=lines[1],iterations=int(lines[2]),restart_cycles=len(re.findall('tritr, ene',text)),rng_state=(work/'rng.dat').read_text().splitlines())
        metrics=(work/'metrics.dat').read_text().split();r.update(h_applications=int(metrics[0]),solver_seconds=float(metrics[1]))
        if status!=0:
            if result.returncode==0:raise ValueError('failure status returned success process')
            r['classification']='nonconvergence' if status==1 else 'numerical_failure' if status==2 else 'invalid_parameters'
            estimates=np.loadtxt(work/'last-estimates.dat',ndmin=2)
            r['last_reported_residuals']=estimates[:,1].tolist()
            if status==1 and t['alg']==2:
                data=np.loadtxt(work/'last-vectors.dat',ndmin=2)
                q=(data[:,0]+1j*data[:,1]).reshape(wanted,n).T
                res=np.linalg.norm(h@q-q*estimates[:,0],axis=0)
                err=abs(estimates[:,0]-t['e'][:wanted])
                r.update(last_actual_residuals=res.tolist(),last_eigenvalue_errors=err.tolist(),
                         last_candidate_meets_strict_energy_and_residual=bool(max(err)<1e-10*max(abs(t['e'])) and max(res)<1e-8*max(abs(t['e']))))
        else:
            if result.returncode:raise ValueError('success status returned failed process')
            values=np.loadtxt(work/'eigenpairs.dat',ndmin=2)[:,1]
            data=np.loadtxt(work/'vectors.dat',ndmin=2);q=(data[:,2]+1j*data[:,3]).reshape(wanted,n).T
            residual=np.linalg.norm(h@q-q*values,axis=0);scale=max(abs(t['e']));orth=float(np.linalg.norm(q.conj().T@q-np.eye(wanted),2))
            nearest=np.min(abs(values[:,None]-t['e'][None,:]),axis=1)
            ground_error=abs(values-t['e'][:wanted]);energy_error=float(max(ground_error))
            # Strict accuracy targets are independent assessment criteria, not
            # claims that the production energy-change knob promises these.
            strict_e=energy_error<=1e-10*max(scale,1.)
            idx=np.flatnonzero(t['e']<=t['e'][wanted-1]+1e-12*max(scale,1.))
            v=t['u'][:,idx];leak=float(np.linalg.norm(q-v@(v.conj().T@q),2))
            projector=None
            if len(idx)==wanted:projector=float(np.linalg.norm(q@q.conj().T-v@v.conj().T,2))
            broad_idx=np.flatnonzero(t['e']<=t['e'][wanted-1]+1e-8*max(scale,1.))
            broad=t['u'][:,broad_idx];cluster_leak=float(np.linalg.norm(q-broad@(broad.conj().T@q),2))
            valid=bool(np.isfinite(residual).all() and max(residual)/max(scale,np.finfo(float).tiny)<2e-3 and orth<1e-8)
            r.update(eigenvalues=values.tolist(),reference_eigenvalues=t['e'][:wanted].tolist(),eigenvalue_errors=ground_error.tolist(),nearest_spectral_errors=nearest.tolist(),residual_norms=residual.tolist(),normalized_residuals=(residual/scale).tolist(),orthogonality_error=orth,requested_space_leakage=leak,cluster_space_leakage=cluster_leak,projector_error=projector,strict_energy_target_met=bool(strict_e),requested_space_target_met=bool(leak<1e-5),valid_pair_under_loose_assessment=valid)
            r['classification']='accepted_strict' if strict_e and leak<1e-5 else 'accepted_accuracy_or_inclusion_warning'
            if not np.isfinite(values).all():raise ValueError('nonfinite successful eigenvalues')
    except (OSError,ValueError,subprocess.TimeoutExpired) as error:
        r.update(classification='harness_or_runtime_error',error=str(error))
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--driver',type=Path,required=True);p.add_argument('--output',type=Path)
    args=p.parse_args();driver=args.driver.resolve(strict=True);root=Path(tempfile.mkdtemp(prefix='ed2-krylov-stress-')) if args.output is None else args.output.resolve()
    if args.output is not None:root.mkdir(parents=True,exist_ok=False)
    records=[]
    for i,t in enumerate(definitions()):
        for threads in [1,2]:
            r=run(driver,root,t,threads,i);records.append(r)
        if i%20==0:print(i,records[-1]['classification'],flush=True)
    counts={c:sum(r['classification']==c for r in records) for c in sorted(set(r['classification'] for r in records))}
    report=dict(driver_sha256=digest(driver),harness_sha256=digest(Path(__file__)),numpy=np.__version__,matrix_generator='QR complex Gaussian, seed=8191+n; analytic eigenvalues and Q retained by formula',counts=counts,cases=records)
    (root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(counts,flush=True)
    return int(counts.get('harness_or_runtime_error',0)>0)


if __name__=='__main__':raise SystemExit(main())
