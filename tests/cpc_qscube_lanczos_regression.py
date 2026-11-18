#!/usr/bin/env python3
"""Independent acceptance/provenance test for the restored QS-Cube two-stage ALG=1.

The two-stage Lanczos *pattern* predates this fork (original published
QS-Cube/ED and QS-Cube/ED2 source/lanczos.f90). This tests the modern residual safety contract, not an
algorithmic novelty claim. The production Krylov kernels are called through
krylov_policy_driver.F90; independently generated Hermitian matrices are used.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def fixtures():
    rng = np.random.default_rng(20261118)
    n=12
    h=np.diag(np.linspace(-2., 2., n)).astype(complex)
    h+=0.18*np.diag(np.ones(n-1),1)+0.18*np.diag(np.ones(n-1),-1)
    yield ('real-tridiagonal',h,None,n,1e-10,None)
    for name,n,ev in [
        ('complex-hermitian',16,np.linspace(-3.,2.,16)),
        ('clustered',20,np.r_[[-1.,-1.+1e-10,-1.+2e-10],np.linspace(-.5,2.,17)]),
        ('degenerate',12,np.r_[[-1.,-1.],np.linspace(-.2,2.,10)]),
        ('near-zero',10,np.r_[[0.,1e-10,2e-10],np.linspace(.5,1.5,7)]),
    ]:
        a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        q,_=np.linalg.qr(a)
        yield (name,(q*ev)@q.conj().T,None,n,1e-10,None)
    yield ('one-dimensional',np.array([[-.7]],complex),np.array([1+0j]),1,1e-10,None)
    yield ('happy-breakdown',np.eye(7,dtype=complex)*-.5,np.arange(1,8,dtype=complex),7,1e-10,None)
    yield ('excluded-ground',np.diag([-2.,-1.,0.,1.]).astype(complex),np.array([0,1,0,0],complex),4,1e-10,None)
    n=20;diag=np.linspace(-2.,2.,n);h=np.diag(diag).astype(complex)
    h+=0.12*np.diag(np.ones(n-1),1)+0.12*np.diag(np.ones(n-1),-1)
    yield ('budget-exhausted',h,None,3,1e-10,None)
    yield ('strict-residual-policy',h,None,n,1e-10,1e-11)
    yield ('impossible-residual-policy',h,None,n,1e-10,1e-300)


def input_text(h,start,budget,tol,policy,alg=1):
    n=h.shape[0]
    controls=(f'&problem n={n},alg={alg},wanted=1,keep=3,capacity=6,maxiter={budget},'
              f'miniter=1,interval=3,tol={tol:.17e},seed=19,use_start=.true.,'
              'enable_diagnostics=.true.,use_seed_option=.false.,')
    if policy is not None:
        controls+=(f'policy_enabled=.true.,scale_provided=.true.,'
                   f'residual_atol={policy:.17e},residual_rtol=0d0,'
                   'residual_scale=1d0,energy_origin=0d0,')
    controls+=' /\n'
    return (controls+''.join(f'{x.real:.17e} {x.imag:.17e}\n' for x in h.flat)+
            ''.join(f'{x.real:.17e} {x.imag:.17e}\n' for x in start))


def run(driver,out,case,threads,mode):
    name,h,start,budget,tol,policy=case
    n=h.shape[0]
    if start is None:
        start=np.array([np.cos(1.73*(i+1))+1j*np.sin(.94*(i+1)) for i in range(n)],complex)
    work=out/f'{name}-t{threads}-{mode}';work.mkdir(parents=True)
    inp=input_text(h,start,budget,tol,policy)
    (work/'input.dat').write_text(inp)
    env={**os.environ,'OMP_NUM_THREADS':str(threads),'OMP_DYNAMIC':'FALSE',
         'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1', 'OMP_PROC_BIND':'close'}
    env.pop('ED2_LANCZOS_TWO_STAGE',None)
    if mode=='two-stage':env['ED2_LANCZOS_TWO_STAGE']='1'
    elif mode=='invalid':env['ED2_LANCZOS_TWO_STAGE']='not-a-mode'
    with (work/'input.dat').open('rb') as fd,(work/'stdout.log').open('wb') as log:
        cp=subprocess.run([str(driver)],cwd=work,stdin=fd,stdout=log,stderr=subprocess.STDOUT,
                          env=env,timeout=90)
    status_lines=(work/'status.dat').read_text().splitlines()
    status=int(status_lines[0]);iterations=int(status_lines[2])
    diagnostics=dict(x.split('=',1) for x in (work/'diagnostics.dat').read_text().splitlines() if '=' in x)
    h_calls=int((work/'metrics.dat').read_text().split()[0])
    assert h_calls==int(diagnostics['h_applications'])
    assert (status==0)==(cp.returncode==0)
    if mode=='invalid':
        assert status==-1 and 'invalid ED2_LANCZOS_TWO_STAGE' in status_lines[1]
        assert not (work/'eigenpairs.dat').exists()
    else:
        assert diagnostics['lowest_state_inclusion'].startswith('not_assessed')
        assert diagnostics.get('lanczos_mode','legacy')==('two_stage' if mode=='two-stage' else 'legacy')
    record=dict(case=name,mode=mode,threads=threads,status=status,reason=status_lines[1],
                iterations=iterations,h_calls=h_calls,exit_code=cp.returncode,
                source_input_sha256=hashlib.sha256(inp.encode()).hexdigest())
    if mode=='two-stage':
        record['reconstructions']=int(diagnostics['vector_reconstructions'])
        if status==0: assert record['reconstructions']>=1
        if diagnostics['centered_variance_available']=='T':
            assert float(diagnostics['centered_variance'])>=0
            record['centered_variance']=float(diagnostics['centered_variance'])
    if status==0:
        values=np.loadtxt(work/'eigenpairs.dat',ndmin=2)
        reported_energy=float(values[0,1]);reported_residual=float(values[0,2])
        a=np.loadtxt(work/'vectors.dat',ndmin=2)
        vector=a[:,2]+1j*a[:,3]
        norm=np.linalg.norm(vector)
        direct=np.linalg.norm(h@vector-reported_energy*vector)
        rayleigh=np.vdot(vector,h@vector).real/(norm**2)
        exact_eig=np.linalg.eigvalsh(h)[0]
        assert abs(norm-1)<1e-11,(name,mode,'norm',norm)
        assert abs(reported_residual-direct)<5e-11,(name,mode,'residual',reported_residual,direct)
        assert direct<1e-7,(name,mode,'direct',direct)
        assert abs(reported_energy-rayleigh)<1e-10,(name,mode,'Rayleigh',reported_energy,rayleigh)
        if policy is not None:assert direct<=policy*1.001
        if name!='excluded-ground':assert abs(exact_eig-reported_energy)<1e-7,(name,mode,'energy',exact_eig,reported_energy)
        else:assert reported_energy>exact_eig+0.5, 'Excluded ground must not be misrepresented as identified'
        record.update(energy=reported_energy,direct_residual=direct,lowest_exact=exact_eig)
        if mode=='two-stage':
            assert abs(direct**2-float(diagnostics['centered_variance']))<1e-10, (name,'variance')
    else:
        assert not (work/'eigenpairs.dat').exists()
    return record



def physical_216(exe,base,threads):
    """Full production ED2 periodic 6^3 three-magnon Hamiltonian check."""
    from cpc_qs3_ed_large_few_magnon import create
    name='216-n3-physical'
    values=[]
    for t in threads:
        for mode in ('legacy','two-stage'):
            work=base/f'{name}-t{t}-{mode}'
            create(work,(6,6,6),3,(0,0,0),False)
            deck=(work/'input.dat').read_text()
            assert 'wk_dim=0,MNTE=-1' in deck
            deck=deck.replace('wk_dim=0,MNTE=-1','wk_dim=-1,MNTE=19')
            (work/'input.dat').write_text(deck)
            env={**os.environ,'OMP_NUM_THREADS':str(t),'OMP_DYNAMIC':'FALSE',
                 'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
                 'ED2_DIAGNOSTICS':'diagnostics.txt','ED2_RNG_SEED':'1984',
                 'ED2_RESIDUAL_ATOL':'1e-8'}
            env.pop('ED2_LANCZOS_TWO_STAGE',None)
            if mode=='two-stage':env['ED2_LANCZOS_TWO_STAGE']='1'
            with (work/'input.dat').open('rb') as fd,(work/'stdout.log').open('wb') as log:
                cp=subprocess.run([str(exe)],cwd=work,stdin=fd,stdout=log,
                                  stderr=subprocess.STDOUT,env=env,timeout=180)
            assert cp.returncode==0,(name,t,mode,(work/'stdout.log').read_text()[-1500:])
            diag=dict(x.split('=',1) for x in (work/'diagnostics.txt').read_text().splitlines() if '=' in x)
            hcalls=int(diag['h_applications']); residual=float(diag['pair_1_residual']) if 'pair_1_residual' in diag else None
            lines=(work/'eigenvalues.dat').read_text().splitlines()
            energy=float(lines[0].split()[-1].replace('D','E'))
            assert abs(energy-148.5648907577122)<2e-8,(mode,t,energy)
            assert diag['lowest_state_inclusion'].startswith('not_assessed')
            assert diag['lanczos_mode'].strip()==('two_stage' if mode=='two-stage' else 'legacy')
            if mode=='two-stage':assert int(diag['vector_reconstructions'])==1
            values.append(dict(case=name,mode=mode,threads=t,energy=energy,h_calls=hcalls,
                               residual=float(diag.get('pair_1_residual',0)) if 'pair_1_residual' in diag else None))
            print('PHYSICAL',name,mode,'OMP',t,'H',hcalls,'E',energy,flush=True)
        old,new=values[-2:]
        assert new['h_calls']<old['h_calls'],(old,new)
        assert abs(new['energy']-old['energy'])<1e-10,(old,new)
    return values

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--driver',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--threads',default='1,4')
    parser.add_argument('--exe',type=Path,help='Optional production ED2 binary for 216-n3 numerical regression')
    args=parser.parse_args()
    driver=args.driver.resolve(strict=True)
    base=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='ed2-qscube-twostage-'))
    if args.output:base.mkdir(parents=True,exist_ok=False)
    threads=[int(x) for x in args.threads.split(',')]
    if not threads or min(threads)<1 or len(set(threads))!=len(threads):parser.error('Invalid threads')
    records=[]
    for case in fixtures():
        for t in threads:
            a=run(driver,base,case,t,'legacy');b=run(driver,base,case,t,'two-stage')
            # Same acceptance classification, but permitted different stopping iterations.
            assert a['status']==b['status'],(case[0],t,a['status'],b['status'])
            if a['status']==0:
                assert abs(a['energy']-b['energy'])<1e-8,(case[0],t)
            records.extend([a,b]);print('PASS',case[0],'OMP',t,a['status'],a['h_calls'],b['h_calls'],flush=True)
    for t in threads:
        bad=run(driver,base,('invalid-environment',np.eye(4,dtype=complex),None,4,1e-10,None),t,'invalid')
        records.append(bad)
    physical = physical_216(args.exe.resolve(strict=True),base,threads) if args.exe else []
    categories={s:sum(x['status']==s for x in records) for s in sorted({x['status'] for x in records})}
    report=dict(schema=1,historical_algorithm='Published QS-Cube/ED and QS-Cube/ED2 source/lanczos.f90: existing two-stage replay',
                driver_sha256=hashlib.sha256(driver.read_bytes()).hexdigest(),
                python_numpy=np.__version__,n_fixtures=len(list(fixtures())),threads=threads,
                counts=categories,records=records,physical=physical,passed=True,
                caveat='Residual accuracy does not certify inclusion of the lowest eigenstate')
    (base/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('QS-CUBE TWO-STAGE REGRESSION',len(records),'DRIVER RUNS;',len(physical),'PHYSICAL RUNS PASS')

if __name__=='__main__':main()
