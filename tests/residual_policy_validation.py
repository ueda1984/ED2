#!/usr/bin/env python3
"""Opt-in residual gates: matched starts, honest failures and independent accuracy.

Characterization includes nonconvergence. Exit zero means assertions/protocols
passed, not that every solve succeeded. --before-driver verifies Phase 1-F bytes.
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

sys.dont_write_bytecode=True
import numpy as np
import precision_contract as precision
from diagnostics_reproducibility import telemetry,same,environment,cli_run,sha
from numerical_baseline import model

CONTROLS=['ED2_RESIDUAL_ATOL','ED2_RESIDUAL_RTOL','ED2_RESIDUAL_SCALE','ED2_ENERGY_ORIGIN']


def run(driver,root,name,t,threads,p=None,diagnostics=True):
    work=root/name;work.mkdir();n=t['n'];start=t.get('start',np.zeros(n,complex))
    header=(f"&problem n={n},alg={t['alg']},wanted={t['wanted']},keep={t['keep']},capacity={t['capacity']},"
            f"maxiter={t['budget']},miniter=1,interval=1,tol={t['tol']},seed={t['seed']},"
            f"use_start={'.true.' if 'start' in t else '.false.'},enable_diagnostics={'.true.' if diagnostics else '.false.'},"
            'use_seed_option=.true.,')
    if p is not None:
        header+='policy_enabled=.true.,'+','.join(f'{k}={str(v).lower() if isinstance(v,str) else v}' for k,v in p.items())+','
    header+='/\n'
    deck=header+''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in t['h'].flatten())
    deck+=''.join(f'{z.real:.17e} {z.imag:.17e}\n' for z in start)
    (work/'input.dat').write_text(deck)
    env=environment(threads)
    for key in CONTROLS:env.pop(key,None)
    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
        proc=subprocess.run([str(driver)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    lines=(work/'status.dat').read_text().splitlines();status=int(lines[0])
    assert (status==0)==(proc.returncode==0)
    met=(work/'metrics.dat').read_text().split()
    r=dict(work=work,status=status,reason=lines[1],iterations=int(lines[2]),h_calls=int(met[0]),
           seconds=float(met[1]),rss_kib=int(met[2]),exit_code=proc.returncode,input_sha256=sha(work/'input.dat'))
    if diagnostics:
        fields,pairs=telemetry(work/'diagnostics.dat');r.update(fields=fields,pairs=pairs)
        assert int(fields['h_applications'])==r['h_calls']
        assert fields['lowest_state_inclusion'].startswith('not_assessed')
        assert fields['roundoff_limit_proven']=='F'
    if status!=0:assert not (work/'eigenpairs.dat').exists() and not (work/'vectors.dat').exists()
    if status==0 or status==1 and t['alg']==2:
        if status==0:
            a=np.loadtxt(work/'vectors.dat',ndmin=2);q=(a[:,2]+1j*a[:,3]).reshape(t['wanted'],n).T
            values=np.loadtxt(work/'eigenpairs.dat',ndmin=2)[:,1]
        else:
            a=np.loadtxt(work/'last-vectors.dat',ndmin=2);q=(a[:,0]+1j*a[:,1]).reshape(t['wanted'],n).T
            values=np.loadtxt(work/'last-estimates.dat',ndmin=2)[:,0]
        assess=precision.diagnostics(t['h'],t['e'],t['u'],q,values)
        if t['family'].startswith('physical'):
            _,_,_,spin=precision.physical();lm,cf=precision.observables(q,spin)
            expected_lm,expected_cf=precision.observables(t['u'][:,:1],spin)
            assess.update(magnetization_error=float(np.max(abs(lm-expected_lm))),correlation_error=float(np.max(abs(cf-expected_cf))),
                          physical_assertions_met=bool(np.allclose(lm,expected_lm,atol=1e-9,rtol=1e-9) and np.allclose(cf,expected_cf,atol=1e-9,rtol=1e-9)))
        if t['family']=='spin-dimer-tiny-gap':
            _,spin=model(2,[],np.zeros((2,3)))
            singlet=np.array([0.,1.,-1.,0.],complex)/np.sqrt(2.)
            mag=np.array([np.vdot(q[:,0],a@q[:,0]).real for site in spin for a in site])
            cf=np.array([np.vdot(q[:,0],spin[0][a]@spin[1][b]@q[:,0]).real for b in range(3) for a in range(3)])
            exact_cf=np.array([np.vdot(singlet,spin[0][a]@spin[1][b]@singlet).real for b in range(3) for a in range(3)])
            assess.update(magnetization_error=float(np.max(abs(mag))),correlation_error=float(np.max(abs(cf-exact_cf))),
                          physical_assertions_met=bool(np.allclose(mag,0,atol=1e-9,rtol=1e-9) and np.allclose(cf,exact_cf,atol=1e-9,rtol=1e-9)))
        r['assessment']=assess
    if p is not None and diagnostics:
        assert fields['residual_policy_enabled']=='T'
        for j,pair in enumerate(pairs,1):
            if not pair['available']:continue
            if f'requested_residual_threshold_{j}' not in fields:
                assert status!=0
                continue
            scale=max(p.get('residual_scale',0.),abs(pair['ritz']-p.get('energy_origin',0.))) if p.get('residual_rtol',0.)>0 else 0.
            threshold=p.get('residual_atol',0.)+p.get('residual_rtol',0.)*scale
            np.testing.assert_allclose(float(fields[f'requested_residual_scale_{j}']),scale,rtol=3e-15,atol=0)
            np.testing.assert_allclose(float(fields[f'requested_residual_threshold_{j}']),threshold,rtol=3e-15,atol=0)
            met=pair['residual']<=threshold
            assert (fields[f'requested_residual_met_{j}']=='T')==met
            assert (fields[f'combined_residual_certified_{j}']=='T')==(met and pair['certified'])
        if status==0:
            assert len(pairs)==t['wanted'] and fields['all_requested_residuals_met']=='T'
            assert all(fields[f'combined_residual_certified_{j}']=='T' for j in range(1,t['wanted']+1))
            assert max(pair['residual']-(p.get('residual_atol',0.)+p.get('residual_rtol',0.)*float(fields[f'requested_residual_scale_{j}'])) for j,pair in enumerate(pairs,1))<=0
    return r


def reduced(r):return {k:v for k,v in r.items() if k!='work'}


def policy_cases():
    yield from precision.definitions()
    # A complete block with an exact zero first pair and inexact later pair.
    e=np.array([0.,1.,2.,3.]);u=np.eye(4,dtype=complex)
    yield dict(family='all-pairs-full-space',n=4,h=np.diag(e),e=e,u=u,alg=2,tol=1e-12,seed=19,
               wanted=2,keep=2,capacity=4,budget=1,start=u[:,0])
    # A genuine two-spin Heisenberg operator with a tiny gap and an identity
    # shift. The analytic singlet supplies gauge-independent physical truth.
    h,_=model(2,[(0,1,[1e-15]*3+[0.]*6)],np.zeros((2,3)));h-=np.eye(4)
    u=np.column_stack([np.array([0,1,-1,0])/np.sqrt(2.),np.array([1,0,0,0]),
                       np.array([0,1,1,0])/np.sqrt(2.),np.array([0,0,0,1])]).astype(complex)
    e=np.real(np.diag(u.conj().T@h@u));order=np.argsort(e);e=e[order];u=u[:,order]
    for alg in [1,2]:
        yield dict(family='spin-dimer-tiny-gap',n=4,h=h,e=e,u=u,alg=alg,tol=1e-12,seed=19,
                   wanted=1 if alg==1 else 2,keep=2,capacity=4,budget=4,origin=-1.,
                   start=np.array([1,1,0,0],complex)/np.sqrt(2.))
    # Exact one-dimensional completion; no epsilon floor may be substituted.
    for alg in [1,2]:
        yield dict(family='exact-zero-scalar',n=1,h=np.zeros((1,1),complex),e=np.zeros(1),u=np.eye(1,dtype=complex),
                   alg=alg,tol=1e-20,seed=19,wanted=1,keep=1,capacity=1,budget=1,start=np.ones(1,complex))


def cli_checks(exe,root):
    rows=[]
    # cli_run inherits process environment; clear new controls for its baseline.
    for key in CONTROLS:os.environ.pop(key,None)
    for alg in [1,2]:
        for threads in [1,2]:
            base,rc=cli_run(exe,root,f'cli-a{alg}-t{threads}-legacy',alg,threads,19,True)
            assert rc==0
            opts={'ED2_RESIDUAL_ATOL':'1e-10'}
            strict,code=cli_run(exe,root,f'cli-a{alg}-t{threads}-strict',alg,threads,19,True,overrides=opts)
            fields,pairs=telemetry(strict/'diagnostics.dat')
            assert (code==0)==(fields['status']=='0')
            if code==0:
                assert fields['all_requested_residuals_met']=='T'
                assert all(p['residual']<=1e-10 for p in pairs)
            else:assert not (strict/'eigenvalues.dat').exists() and not (strict/'wf1.bin').exists()
            loose,code=cli_run(exe,root,f'cli-a{alg}-t{threads}-loose',alg,threads,19,True,overrides={'ED2_RESIDUAL_ATOL':'1e3'})
            assert code==0
            for file in ['eigenvalues.dat','local_mag.dat','two_body_cf_xyz.dat','wf1.bin']:
                assert (base/file).read_bytes()==(loose/file).read_bytes()
            def normal(path):return [line for line in (path/'stdout.log').read_text().splitlines() if not line.strip().startswith('Time:')]
            assert normal(base)==normal(loose)
            # Strict enough to be unattainable in this full numerical solve.
            fail,code=cli_run(exe,root,f'cli-a{alg}-t{threads}-reject',alg,threads,19,True,
                              overrides={'ED2_RESIDUAL_ATOL':'1e-300'},budget=1)
            assert code!=0 and not (fail/'eigenvalues.dat').exists() and not (fail/'wf1.bin').exists()
            f,_=telemetry(fail/'diagnostics.dat');assert f['status']!='0'
            rows.append(dict(alg=alg,threads=threads,strict_exit=int(fields['status']),strict_residuals=[p['residual'] for p in pairs],
                             failure_status=int(f['status']),failure_reason=f['reason']))
    invalid=[{'ED2_RESIDUAL_ATOL':'0'},{'ED2_RESIDUAL_ATOL':'-1'}, {'ED2_RESIDUAL_ATOL':'NaN'},
             {'ED2_RESIDUAL_ATOL':'Inf'}, {'ED2_RESIDUAL_RTOL':'1e-8'}, {'ED2_RESIDUAL_SCALE':'1'},
             {'ED2_ENERGY_ORIGIN':'1'}, {'ED2_RESIDUAL_ATOL':'1e-10 extra'},
             {'ED2_RESIDUAL_RTOL':'1e308','ED2_RESIDUAL_SCALE':'2'},
             {'ED2_RESIDUAL_RTOL':'1e-8','ED2_RESIDUAL_SCALE':'-1'},
             {'ED2_RESIDUAL_RTOL':'1e-8','ED2_RESIDUAL_SCALE':'NaN'},
             {'ED2_RESIDUAL_ATOL':'1e-8','ED2_ENERGY_ORIGIN':'NaN'}]
    for i,opts in enumerate(invalid):
        work,rc=cli_run(exe,root,f'invalid-policy-{i}',1,1,19,True,overrides=opts)
        assert rc!=0 and not (work/'eigenvalues.dat').exists() and not (work/'wf1.bin').exists()
    # A relative-only request with an explicitly supplied zero scale is valid.
    for alg in [1,2]:
        work,rc=cli_run(exe,root,f'cli-relative-a{alg}',alg,2,73,True,overrides={
            'ED2_RESIDUAL_RTOL':'1e-8','ED2_RESIDUAL_SCALE':'0','ED2_ENERGY_ORIGIN':'-1'})
        f,pairs=telemetry(work/'diagnostics.dat');assert f['requested_scale_provided']=='T'
        if rc==0:assert all(f[f'requested_residual_met_{j}']=='T' for j in range(1,len(pairs)+1))
    for alg in [1,2]:
        work,rc=cli_run(exe,root,f'cli-combined-a{alg}',alg,2,19,True,overrides={
            'ED2_RESIDUAL_ATOL':'1e-12','ED2_RESIDUAL_RTOL':'1e-10','ED2_RESIDUAL_SCALE':'2','ED2_ENERGY_ORIGIN':'0'})
        assert rc==0
        f,pairs=telemetry(work/'diagnostics.dat')
        for j,pair in enumerate(pairs,1):
            np.testing.assert_allclose(float(f[f'requested_residual_threshold_{j}']),1e-12+1e-10*max(2.,abs(pair['ritz'])),rtol=3e-15,atol=0)
    return dict(paired=rows,invalid_rejections=len(invalid),relative_cli=2,combined_cli=2)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--driver',type=Path,required=True)
    p.add_argument('--before-driver',type=Path);p.add_argument('--exe',type=Path,required=True)
    p.add_argument('--output',type=Path);a=p.parse_args()
    driver=a.driver.resolve(strict=True);before=a.before_driver.resolve(strict=True) if a.before_driver else None
    root=a.output.resolve() if a.output else Path(tempfile.mkdtemp(prefix='ed2-residual-'))
    if a.output:root.mkdir(parents=True,exist_ok=False)
    controls=precision.analytic_controls();rows=[];equivalence=0
    for index,t in enumerate(policy_cases()):
        for threads in [1,2]:
            name=f'{index:04d}-t{threads}'
            legacy=run(driver,root,name+'-legacy',t,threads)
            if before:
                old=run(before,root,name+'-before',t,threads)
                same(old,legacy);equivalence+=1
            shifted=t['family'].startswith('scale-1.0-shift-1000000')
            origin=t.get('origin',1e6 if shifted else 0.)
            b=float(max(abs(t['e']-origin)))
            policies={'atol-1e-10':dict(residual_atol=1e-10),
                      'rtol-1e-10-centered':dict(residual_rtol=1e-10,residual_scale=b,energy_origin=origin,scale_provided='.true.'),
                      'atol-1e-300':dict(residual_atol=1e-300)}
            # Absolute loose request must not change any legacy trajectory.
            loose=run(driver,root,name+'-loose',t,threads,dict(residual_atol=1e30))
            same(legacy,loose)
            group=dict(id=name,family=t['family'],alg=t['alg'],threads=threads,seed=t['seed'],tol=t['tol'],
                       legacy=reduced(legacy),policies={})
            for label,policy in policies.items():
                r=run(driver,root,name+'-'+label,t,threads,policy)
                group['policies'][label]=reduced(r)
                if t['family']=='all-pairs-full-space' and label=='atol-1e-300':
                    assert legacy['status']==0 and r['status']!=0
                    assert r['fields']['requested_residual_met_1']=='T' and r['fields']['requested_residual_met_2']=='F'
                if t['family']=='excluded-ground' and label=='atol-1e-10':
                    assert r['status']==0 and r['assessment']['projector_error']==1.
                if t['family']=='exact-zero-scalar':assert r['status']==0
            if t['family']=='spin-dimer-tiny-gap' and t['alg']==1:
                accepted=group['policies']['atol-1e-10']
                assert accepted['status']==0 and accepted['assessment']['magnetization_error']>.1
                assert accepted['assessment']['correlation_error']>.1
            if t['family']=='all-pairs-full-space':
                edge=max(pair['residual'] for pair in legacy['pairs'])
                boundary=run(driver,root,name+'-inclusive-boundary',t,threads,dict(residual_atol=edge))
                assert boundary['status']==0
                same(legacy,boundary)
            rows.append(group)
    # Invalid policy objects are rejected before any H action in both kernels.
    template=next(t for t in policy_cases() if t['family']=='physical-random')
    invalid_kernel=[]
    bad=[dict(residual_atol=0.),dict(residual_atol=-1.),dict(residual_atol=float('nan')),
         dict(residual_rtol=1e-8),dict(residual_rtol=1e-8,residual_scale=-1.,scale_provided='.true.'),
         dict(residual_rtol=1e308,residual_scale=2.,scale_provided='.true.')]
    for alg in [1,2]:
        t=dict(template,alg=alg,wanted=1 if alg==1 else 2)
        for i,policy in enumerate(bad):
            r=run(driver,root,f'invalid-kernel-a{alg}-{i}',t,1,policy)
            assert r['status']==-1 and r['h_calls']==0
            invalid_kernel.append(dict(alg=alg,index=i,reason=r['reason']))
    huge=dict(family='overflow-guard',n=1,h=np.array([[1e308]],complex),e=np.array([1e308]),u=np.eye(1,dtype=complex),
              alg=1,tol=1e-12,seed=19,wanted=1,keep=1,capacity=1,budget=1,start=np.ones(1,complex))
    r=run(driver,root,'centered-overflow',huge,1,dict(residual_rtol=1e-10,residual_scale=0.,energy_origin=-1e308,scale_provided='.true.'))
    assert r['status']==2 and 'nonfinite requested' in r['reason']
    r=run(driver,root,'absolute-origin-unused',huge,1,dict(residual_atol=1e-10,energy_origin=-1e308))
    assert r['status']==0
    assert max(controls['constructed_physical_state']['residual_norms'])>1e-10
    cli=cli_checks(a.exe.resolve(strict=True),root)
    summary=dict(legacy=dict(Counter(str(r['legacy']['status']) for r in rows)),
                 policies={label:dict(Counter(str(r['policies'][label]['status']) for r in rows)) for label in policies})
    report=dict(schema=1,harness_sha256=sha(Path(__file__)),driver_sha256=sha(driver),
                before_driver_sha256=sha(before) if before else None,executable_sha256=sha(a.exe.resolve()),
                before_equivalence=equivalence,controls=controls,cli=cli,invalid_kernel=invalid_kernel,overflow_guards=2,summary=summary,cases=rows)
    (root/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(rows),before_equivalence=equivalence,summary=summary,cli=cli)))
    return 0


if __name__=='__main__':raise SystemExit(main())
