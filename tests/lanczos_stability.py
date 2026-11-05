#!/usr/bin/env python3
"""Residual, invariant-subspace and failure-status checks of production kernels."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np


def rotated(values, seed=3127):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(len(values),len(values))) + 1j*rng.normal(size=(len(values),len(values))))
    return (q*np.array(values))@q.conj().T


def cases():
    out=[]
    def add(name,h,alg,**kwargs):
        n=len(h)
        item=dict(name=name,h=h,alg=alg,wanted=1 if alg==1 else min(2,n),
                  keep=min(4,n),capacity=min(8,n),maxiter=n if alg==1 else 200,
                  miniter=5,interval=1,tol=1e-14,status=0)
        item.update(kwargs);out.append(item)
    for alg in [1,2]:
        add('zero',np.zeros((6,6),complex),alg)
        add('scalar',np.eye(6)*3.25,alg)
        for value in [0.,2.5]:
            add('one-dimensional-'+str(value),np.array([[value]]),alg)
        add('zero-energy',rotated([0.,.5,1.,2.]),alg)
        add('exact-ground-start',np.diag([-2.,-.5,1.,3.]),alg,start=np.array([1,0,0,0],complex))
        add('zero-start',np.eye(4),alg,start=np.zeros(4,complex),status=-1)
        add('nonfinite-start',np.eye(4),alg,start=np.array([complex(float('nan')),0,0,0]),status=-1)
        add('invalid-budget',np.eye(4),alg,maxiter=0,status=-1)
        add('zero-closed-block',np.zeros((16,16),complex),alg,reason='happy breakdown')
        add('scalar-closed-block',np.eye(16)*3.25,alg)
        add('minimal-budget-closed-block',np.zeros((16,16),complex),alg,maxiter=1,reason='happy breakdown')
        add('insufficient-iterations',rotated(np.linspace(-3.,4.,24)),alg,maxiter=1,status=1)
        for scale in [1e-9,1e9]:
            add('scaled-'+str(scale),rotated(np.linspace(-2.,3.,8))*scale,alg)
    add('early-check',rotated([-2.,-.5,1.,3.]),1,miniter=1)
    add('near-breakdown',np.array([[1.,1e-15],[1e-15,2.]]),1,start=np.array([1,0],complex),reason='numerical near-breakdown')
    add('exact-excited-start',np.diag([-2.,-.5,1.,3.]),1,start=np.array([0,0,0,1],complex),expected_value=3.,reason='happy breakdown')
    add('excited-start-restart',np.diag([-2.,-.5,1.,3.]),2,start=np.array([0,0,0,1],complex))
    add('degenerate-full-space',rotated([-2.,-2.,.5,1.,2.,3.]),2,capacity=6,keep=3,projector=True)
    add('degenerate-restarts',rotated([-2.,-2.]+list(np.linspace(-1.,3.,22))),2,projector=True)
    stress=rotated([-2.,-2.]+list(np.linspace(-1.,3.,22)))
    _,stress_vectors=np.linalg.eigh(stress)
    add('excited-start-degenerate-restarts',stress,2,start=stress_vectors[:,-1],projector=True)
    add('diagonal-excited-degenerate-restarts',np.diag([-2.,-2.]+list(np.linspace(-1.,3.,22))),2,
        start=np.eye(24,dtype=complex)[:,-1],projector=True)
    add('zero-energy-restarts',rotated([0.,0.]+list(np.linspace(.5,4.,22))),2,projector=True)
    add('general-restarts',rotated(np.linspace(-3.,4.,24)),2)
    add('minimal-budget-happy',np.zeros((6,6),complex),2,maxiter=1)
    return out


def require(ok,message):
    if not ok:raise ValueError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--driver',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    driver=args.driver.resolve(strict=True)
    if args.output is None:out=Path(tempfile.mkdtemp(prefix='ed2-stability-'))
    else:
        out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    results=[]
    for test in cases():
        for threads in [1,2]:
            h=test['h'];n=len(h);wanted=test['wanted'];name=f"{test['name']}-alg{test['alg']}-omp{threads}"
            work=out/name;work.mkdir()
            start=test.get('start',np.zeros(n,complex))
            header='&problem '+','.join(f'{k}={test[k]}' for k in ['alg','wanted','keep','capacity','maxiter','miniter','interval','tol'])
            header+=f",n={n},use_start={'.true.' if 'start' in test else '.false.'}/\n"
            deck=header+''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in h.flatten())
            deck+=''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in start)
            (work/'input.dat').write_text(deck)
            result={'name':name,'alg':test['alg'],'threads':threads,'expected_status':test['status'],'input_sha256':hashlib.sha256(deck.encode()).hexdigest()}
            try:
                env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
                with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
                    run=subprocess.run([str(driver)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=120)
                lines=(work/'status.dat').read_text().splitlines();status=int(lines[0])
                result.update(exit_code=run.returncode,solver_status=status,reason=lines[1],iterations=int(lines[2]))
                require(status==test['status'],f'unexpected solver status {status}: {lines[1]}')
                if status!=0:
                    require(run.returncode!=0,'failed solve returned success exit')
                    require(not (work/'eigenpairs.dat').exists(),'failed solve wrote eigenpairs')
                else:
                    require(run.returncode==0,'successful solve failed process')
                    pair=np.loadtxt(work/'eigenpairs.dat',ndmin=2);data=np.loadtxt(work/'vectors.dat',ndmin=2)
                    q=(data[:,2]+1j*data[:,3]).reshape(wanted,n).T
                    values=pair[:,1]
                    reference,reference_vectors=np.linalg.eigh(h)
                    expected=reference[:wanted]
                    if 'expected_value' in test:expected=np.array([test['expected_value']])
                    scale=max(np.linalg.norm(h,2),np.finfo(float).tiny)
                    np.testing.assert_allclose(values,expected,atol=1e-9*scale,rtol=1e-9)
                    residual=np.linalg.norm(h@q-q*values,axis=0)
                    relative=residual/scale
                    require(np.isfinite(residual).all(),'nonfinite residual')
                    require(np.max(relative)<2e-7,'actual residual exceeds certification scale')
                    np.testing.assert_allclose(pair[:,2],residual,atol=2e-14*scale,rtol=1e-6)
                    orth=np.linalg.norm(q.conj().T@q-np.eye(wanted),2)
                    require(orth<1e-10,'eigenvectors not orthonormal')
                    # Degeneracies: compare projectors for a complete requested
                    # cluster; for partial scalar clusters validate membership.
                    eigspace=reference_vectors[:,np.abs(reference-values[0])<=1e-8*scale]
                    leakage=np.linalg.norm(q[:,0]-eigspace@(eigspace.conj().T@q[:,0]))
                    require(leakage<1e-5,'vector outside its reference eigenspace')
                    result.update(eigenvalues=values.tolist(),residuals=residual.tolist(),normalized_residuals=relative.tolist(),orthogonality_error=float(orth),eigenspace_leakage=float(leakage))
                    if test.get('projector'):
                        v=reference_vectors[:,:wanted]
                        error=np.linalg.norm(q@q.conj().T-v@v.conj().T,2)
                        require(error<1e-5,'degenerate invariant-subspace projector differs')
                        result['projector_error']=float(error)
                    if 'reason' in test:require(test['reason'] in lines[1],'breakdown misclassified')
                result['status']='passed'
            except (OSError,ValueError,AssertionError,subprocess.TimeoutExpired) as error:
                result.update(status='failed',error=str(error))
            print(result['status'].upper()+': '+name,flush=True)
            results.append(result)
    report={'driver_sha256':hashlib.sha256(driver.read_bytes()).hexdigest(),'numpy':np.__version__,
            'cases':results,'passed':sum(c['status']=='passed' for c in results),'failed':sum(c['status']=='failed' for c in results)}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(report['passed'],'passed;',report['failed'],'failed; evidence:',out)
    return int(report['failed']!=0)


if __name__=='__main__':raise SystemExit(main())
