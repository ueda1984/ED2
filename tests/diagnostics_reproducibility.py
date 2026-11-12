#!/usr/bin/env python3
"""Check telemetry, opt-in initialization and unchanged numerical outcomes.

Optional --before-driver enables exact base-version comparisons. Solver failures
remain failures in characterization counts; assertions check their status and
absence of successful outputs. Python/NumPy remain test-only dependencies.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import numpy as np
import precision_contract as precision
from sparse_field_validation import prepare,FIXTURES,reference


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def telemetry(path):
    fields = {}; pairs = []
    for line in path.read_text().splitlines():
        if '=' in line:
            k,v = line.split('=',1); fields[k] = v.strip()
        elif line and line.split()[0].isdigit():
            p=line.split(); pairs.append(dict(index=int(p[0]),available=p[1]=='T',certified=p[2]=='T',
                                               ritz=float(p[3]),rayleigh=float(p[4]),residual=float(p[5]),
                                               norm=float(p[6]),threshold=float(p[7])))
    return fields,pairs


def environment(threads):
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for name in ['ED2_RNG_SEED','ED2_INITIAL_VECTOR','ED2_DIAGNOSTICS']:
        env.pop(name,None)
    return env


def run(driver,root,name,item,threads,enabled=True,seed_option=False):
    work=root/name; work.mkdir()
    n=item['n']; start=item.get('start',np.zeros(n,complex))
    header=(f"&problem n={n},alg={item['alg']},wanted={item['wanted']},keep={item['keep']},capacity={item['capacity']},"
            f"maxiter={item['budget']},miniter=1,interval=1,tol={item['tol']},seed={item['seed']},"
            f"use_start={'.true.' if 'start' in item else '.false.'},"
            f"enable_diagnostics={'.true.' if enabled else '.false.'},use_seed_option={'.true.' if seed_option else '.false.'}/\n")
    deck=header+''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in item['h'].flatten())
    deck+=''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in start)
    (work/'input.dat').write_text(deck)
    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
        proc=subprocess.run([str(driver)],cwd=work,env=environment(threads),stdin=inp,
                            stdout=log,stderr=subprocess.STDOUT,timeout=120)
    lines=(work/'status.dat').read_text().splitlines();status=int(lines[0])
    assert (status==0)==(proc.returncode==0)
    metrics=(work/'metrics.dat').read_text().split()
    return dict(work=work,status=status,reason=lines[1],iterations=int(lines[2]),h_calls=int(metrics[0]),
                seconds=float(metrics[1]),rss_kib=int(metrics[2]),exit_code=proc.returncode)


def same(a,b):
    for key in ['status','reason','iterations','h_calls','exit_code']:
        assert a[key]==b[key],(key,a,b)
    files=['status.dat','stdout.log','rng.dat','rng-after.dat','eigenpairs.dat','vectors.dat','last-estimates.dat','last-vectors.dat']
    for name in files:
        left=a['work']/name;right=b['work']/name
        assert left.exists()==right.exists(),name
        if left.exists(): assert left.read_bytes()==right.read_bytes(),(name,a['work'],b['work'])


def check_observation(result,item):
    d,pairs=telemetry(result['work']/'diagnostics.dat')
    assert int(d['status'])==result['status'] and d['reason']==result['reason']
    assert int(d['h_applications'])==result['h_calls']
    assert int(d['iterations'])==result['iterations']
    assert int(d['requested_budget'])==item['budget']
    assert d['lowest_state_inclusion'].startswith('not_assessed')
    assert d['roundoff_limit_proven']=='F'
    assert (d['roundoff_risk']=='T')==(0<item['tol']<np.finfo(float).eps)
    assert 'estimate, not a spectral bound' in d['scale_provenance']
    scale=float(d['operator_scale']);eps=np.finfo(float).eps
    np.testing.assert_allclose(float(d['breakdown_threshold']),64*eps*max(scale,np.finfo(float).tiny),atol=0,rtol=2e-15)
    q=None;values=None
    if result['status']==0:
        data=np.loadtxt(result['work']/'vectors.dat',ndmin=2)
        q=(data[:,2]+1j*data[:,3]).reshape(item['wanted'],item['n']).T
        values=np.loadtxt(result['work']/'eigenpairs.dat',ndmin=2)[:,1]
    elif result['status']==1 and item['alg']==2:
        data=np.loadtxt(result['work']/'last-vectors.dat',ndmin=2)
        q=(data[:,0]+1j*data[:,1]).reshape(item['wanted'],item['n']).T
        values=np.loadtxt(result['work']/'last-estimates.dat',ndmin=2)[:,0]
    for j,pair in enumerate(pairs):
        if not pair['available']:continue
        threshold=max(256*eps,np.sqrt(item['tol']))*max(scale,abs(pair['ritz']),np.finfo(float).tiny)
        np.testing.assert_allclose(pair['threshold'],threshold,atol=0,rtol=3e-15)
        assert pair['certified']==(pair['residual']<=pair['threshold'])
        if q is not None:
            v=q[:,j];hv=item['h']@v
            atol=2e-13*max(scale,np.finfo(float).tiny)
            np.testing.assert_allclose(pair['ritz'],values[j],atol=atol,rtol=1e-14)
            np.testing.assert_allclose(pair['norm'],np.linalg.norm(v),atol=3e-14,rtol=0)
            np.testing.assert_allclose(pair['rayleigh'],np.vdot(v,hv).real/np.vdot(v,v).real,atol=atol,rtol=0)
            np.testing.assert_allclose(pair['residual'],np.linalg.norm(hv-pair['ritz']*v),atol=atol,rtol=0)
    if q is not None and d['orthogonality_available']=='T':
        gram=np.max(abs(q.conj().T@q-np.eye(item['wanted'])))
        np.testing.assert_allclose(float(d['max_abs_gram_minus_identity']),gram,atol=3e-13,rtol=0)
    if result['status']!=0:assert 'not a successful output' in d['candidate_role'] or not pairs
    return dict(fields=d,pairs=pairs)


def cli_run(exe,root,name,alg,threads,seed=19,diagnostics=False,initial=None,overrides=None,budget=None):
    work=root/name;prepare(work,FIXTURES['single-site4'],alg,-1)
    text=(work/'input.dat').read_text().replace('wr_wf=0','wr_wf=1')
    if budget is not None:text=text.replace('maxitr=1000',f'maxitr={budget}')
    (work/'input.dat').write_text(text)
    env=environment(threads)
    if seed is not None:env['ED2_RNG_SEED']=str(seed)
    if diagnostics:env['ED2_DIAGNOSTICS']='diagnostics.dat'
    if initial is not None:env['ED2_INITIAL_VECTOR']=str(initial)
    env.update(overrides or {})
    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
        p=subprocess.run([str(exe)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=120)
    return work,p.returncode


def cli_checks(exe,root):
    rows=[]
    for alg in [1,2]:
        for threads in [1,2]:
            for seed in [19,73]:
                name=f'cli-a{alg}-t{threads}-s{seed}'
                a,rc=cli_run(exe,root,name+'-off',alg,threads,seed)
                b,rc2=cli_run(exe,root,name+'-on',alg,threads,seed,True)
                assert rc==rc2==0
                for filename in ['eigenvalues.dat','local_mag.dat','two_body_cf_xyz.dat']+[f'wf{i}.bin' for i in range(1,2 if alg==1 else 3)]:
                    assert (a/filename).read_bytes()==(b/filename).read_bytes(),filename
                def log(p):return '\n'.join(line for line in (p/'stdout.log').read_text().splitlines() if not line.strip().startswith('Time:'))
                assert log(a)==log(b),'legacy normal stdout changed'
                d,pairs=telemetry(b/'diagnostics.dat');assert d['seed_specified']=='T'
                rows.append(dict(name=name,status='passed',sidecar_sha256=sha(b/'diagnostics.dat'),
                                 eigenvalues_sha256=sha(b/'eigenvalues.dat'),h_applications=int(d['h_applications'])))
    start=root/'initial.txt';start.write_text(''.join(f'{i+1}.0 {i/10:.17e}\n' for i in range(16)))
    for alg in [1,2]:
        a,rc=cli_run(exe,root,f'cli-initial-a{alg}-1',alg,2,73,True,start)
        b,rc2=cli_run(exe,root,f'cli-initial-a{alg}-2',alg,2,73,True,start)
        assert rc==rc2==0
        for filename in ['diagnostics.dat','wf1.bin','eigenvalues.dat']:
            assert (a/filename).read_bytes()==(b/filename).read_bytes()
        d,_=telemetry(a/'diagnostics.dat');assert d['explicit_initial']=='T'
        if alg==1:assert int(d['random_vectors'])==0
        rows.append(dict(name=f'explicit-initial-a{alg}',status='passed'))
    for name,seed in [('seed-zero',0),('seed-max',2147483647),('seed-min',-2147483648)]:
        a,rc=cli_run(exe,root,name+'-1',1,2,seed,True)
        b,rc2=cli_run(exe,root,name+'-2',1,2,seed,True)
        assert rc==rc2==0 and (a/'wf1.bin').read_bytes()==(b/'wf1.bin').read_bytes()
        rows.append(dict(name=name,status='passed'))
    short=root/'short.txt';short.write_text('1 0\n')
    zero=root/'zero.txt';zero.write_text('0 0\n'*16)
    nonfinite=root/'nonfinite.txt';nonfinite.write_text('NaN 0\n'+'0 0\n'*15)
    extra=root/'extra.txt';extra.write_text(start.read_text()+'1 0\n')
    existing=root/'existing.txt';existing.write_text('preserve this user file\n')
    bad=[('bad-seed',{'ED2_RNG_SEED':'19x'}),('overflow-seed',{'ED2_RNG_SEED':'9999999999999999'}),
         ('short-vector',{'ED2_INITIAL_VECTOR':str(short)}),('zero-vector',{'ED2_INITIAL_VECTOR':str(zero)}),
         ('nonfinite-vector',{'ED2_INITIAL_VECTOR':str(nonfinite)}),('extra-vector',{'ED2_INITIAL_VECTOR':str(extra)}),
         ('existing-sidecar',{'ED2_DIAGNOSTICS':str(existing)})]
    for name,options in bad:
        work,rc=cli_run(exe,root,name,1,1,19,True,overrides=options)
        assert rc!=0 and not (work/'eigenvalues.dat').exists() and not (work/'wf1.bin').exists()
        rows.append(dict(name=name,status='passed',observed_exit=rc))
    assert existing.read_text()=='preserve this user file\n'
    for alg in [1,2]:
        work,rc=cli_run(exe,root,f'budget-a{alg}',alg,1,19,True,budget=1)
        assert rc!=0 and not (work/'eigenvalues.dat').exists()
        d,_=telemetry(work/'diagnostics.dat');assert d['status']=='1' and d['termination']=='budget_exhausted'
        rows.append(dict(name=f'budget-a{alg}',status='passed',observed_exit=rc))
    return rows


def performance(driver,before,root):
    n=192;rng=np.random.default_rng(173);u,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    e=np.linspace(-1.,2.,n);h=(u*e)@u.conj().T;rows=[]
    for alg in [1,2]:
        item=dict(n=n,h=h,e=e,u=u,alg=alg,tol=1e-10,seed=19,wanted=1 if alg==1 else 2,keep=4,capacity=12,budget=n if alg==1 else 120)
        for repeat in range(7):
            modes=[('off',driver,False),('on',driver,True)]
            if before:modes.insert(0,('before',before,False))
            if repeat%2:modes.reverse()
            results={}
            for label,binary,enabled in modes:
                r=run(binary,root,f'perf-a{alg}-r{repeat}-{label}',item,1,enabled)
                results[label]=r
                rows.append(dict(alg=alg,repeat=repeat,mode=label,seconds=r['seconds'],rss_kib=r['rss_kib'],h_calls=r['h_calls'],status=r['status']))
            same(results['off'],results['on'])
            if before:same(results['before'],results['off'])
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--driver',type=Path,required=True)
    p.add_argument('--before-driver',type=Path);p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--output',type=Path);p.add_argument('--performance',action='store_true')
    a=p.parse_args();driver=a.driver.resolve(strict=True);exe=a.exe.resolve(strict=True)
    before=a.before_driver.resolve(strict=True) if a.before_driver else None
    root=a.output.resolve() if a.output else Path(tempfile.mkdtemp(prefix='ed2-diagnostics-'))
    if a.output:root.mkdir(parents=True,exist_ok=False)
    # Keep the Phase 1-E constructed physical-accuracy counterexample explicit.
    controls=precision.analytic_controls();rows=[];seeds=[];error_paths=[]
    for index,item in enumerate(precision.definitions()):
        for threads in [1,2]:
            name=f'{index:04d}-t{threads}'
            off=run(driver,root,name+'-off',item,threads,False)
            on=run(driver,root,name+'-on',item,threads,True)
            same(off,on);observation=check_observation(on,item)
            if before:
                old=run(before,root,name+'-before',item,threads,False);same(old,off)
            rows.append(dict(id=name,family=item['family'],status=on['status'],reason=on['reason'],
                             h_calls=on['h_calls'],diagnostics=observation,
                             vector_sha256=sha(on['work']/'vectors.dat') if on['status']==0 else None))
            if index in [0,4,16,20,32,34,42,43,44,47,50,53]:
                first=run(driver,root,name+'-seed1',item,threads,True,True)
                second=run(driver,root,name+'-seed2',item,threads,True,True)
                same(first,second)
                assert (first['work']/'diagnostics.dat').read_bytes()==(second['work']/'diagnostics.dat').read_bytes()
                assert (first['work']/'rng.dat').read_bytes()==(first['work']/'rng-after.dat').read_bytes(),'seeded call changed ambient RNG'
                obs=check_observation(first,item)
                if int(obs['fields']['random_vectors'])>0:
                    assert int(obs['fields']['rng_active_threads'])==threads
                    for tid in range(threads):
                        state=list(map(int,obs['fields'][f'rng_initial_state_{tid}'].split()))
                        expected=[(item['seed']+104729*tid+37*k)%2147483646 for k in range(1,len(state)+1)]
                        assert state==expected

                seeds.append(dict(id=name,status=first['status'],h_calls=first['h_calls'],seed=item['seed'],threads=threads))
    h,e,u,_=precision.physical()
    template=dict(n=16,h=h,e=e,u=u,alg=1,tol=1e-12,seed=19,wanted=1,keep=4,capacity=8,budget=16)
    invalid=[('zero-budget',dict(budget=0)),('zero-tolerance',dict(tol=0.)),
             ('nan-tolerance',dict(tol=float('nan'))),('zero-initial',dict(start=np.zeros(16,complex))),
             ('invalid-capacity',dict(alg=2,wanted=2,keep=4,capacity=3))]
    broken=h.copy();broken[0,0]=complex(float('nan'),0)
    invalid.append(('nonfinite-action',dict(h=broken)))
    for name,changes in invalid:
        item=dict(template,**changes)
        off=run(driver,root,name+'-off',item,1,False)
        on=run(driver,root,name+'-on',item,1,True)
        same(off,on)
        if before:
            old=run(before,root,name+'-before',item,1,False);same(old,off)
        assert on['status']==(2 if name=='nonfinite-action' else -1)
        d,_=telemetry(on['work']/'diagnostics.dat')
        assert int(d['h_applications'])==on['h_calls'] and int(d['status'])==on['status']
        error_paths.append(dict(name=name,status=on['status'],termination=d['termination']))
    cli=cli_checks(exe,root)
    perf=performance(driver,before,root) if a.performance else []
    report=dict(schema=1,driver_sha256=sha(driver),before_driver_sha256=sha(before) if before else None,
                executable_sha256=sha(exe),harness_sha256=sha(Path(__file__)),
                controls=controls,comparison_cases=len(rows),before_comparisons=len(rows) if before else 0,
                diagnostic_solver_statuses=dict(Counter(str(r['status']) for r in rows)),
                seeded_replays=seeds,error_paths=error_paths,cli_checks=cli,performance=perf,cases=rows)
    (root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(comparisons=len(rows),before_comparisons=report['before_comparisons'],
                         statuses=report['diagnostic_solver_statuses'],seeded_replays=len(seeds),error_paths=len(error_paths),cli_checks=len(cli),performance_runs=len(perf))))
    return 0


if __name__=='__main__':raise SystemExit(main())
