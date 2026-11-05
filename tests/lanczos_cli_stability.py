#!/usr/bin/env python3
"""Actual CLI eigenvectors/residuals and failure exits, including physical 1D sectors."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
sys.dont_write_bytecode=True
import numpy as np
from numerical_baseline import model


def fixtures():
    full=[dict(name='zero',n=2,lo=0,hi=2,fields=[[0,0,0]]*2,bonds=[]),
          dict(name='scalar-sector',n=3,lo=1,hi=1,fields=[[0,0,2]]*3,bonds=[]),
          dict(name='zero-energy-sector',n=2,lo=1,hi=1,fields=[[0,0,3]]*2,bonds=[]),
          dict(name='one-dimensional',n=2,lo=0,hi=0,fields=[[0,0,1]]*2,bonds=[]),
          dict(name='zero-one-dimensional',n=1,lo=0,hi=0,fields=[[0,0,0]],bonds=[])]
    c=[1.,.8,.7,.2,.1,.3,.1,.3,-.2]
    general=dict(n=4,lo=0,hi=4,fields=[[.13*(i+1),-.17+.03*i,.21-.07*i] for i in range(4)],
                 bonds=[(0,1,c),(1,2,[.7*x for x in c]),(2,3,[.4*x for x in c])])
    full.append(dict(general,name='early-check',algorithms=[1],minitr=1))
    full.append(dict(general,name='insufficient-budget',maxiter=1,expected_failure=True))
    return full


def require(ok,message):
    if not ok:raise ValueError(message)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe',type=Path,required=True);p.add_argument('--operator-exe',type=Path,required=True)
    p.add_argument('--output',type=Path)
    args=p.parse_args();exe=args.exe.resolve(strict=True);operator=args.operator_exe.resolve(strict=True)
    if args.output is None:out=Path(tempfile.mkdtemp(prefix='ed2-cli-stability-'))
    else:out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    results=[]
    for f in fixtures():
        n=f['n'];h,_=model(n,f['bonds'],f['fields'])
        states=[i for i in range(2**n) if f['lo']<=i.bit_count()<=f['hi']]
        expected_h=h[np.ix_(states,states)];reference=np.linalg.eigvalsh(expected_h);dim=len(states)
        for alg in f.get('algorithms',[1,2]):
            wanted=1 if alg==1 else min(2,dim)
            for threads in [1,2]:
                name=f"{f['name']}-alg{alg}-omp{threads}";work=out/name;work.mkdir()
                (work/'fields.dat').write_text(''.join(f'{i+1} '+' '.join(map(str,values))+'\n' for i,values in enumerate(f['fields'])))
                (work/'bonds.dat').write_text(''.join(f'{i+1} {j+1} '+' '.join(map(str,c))+'\n' for i,j,c in f['bonds']))
                (work/'spins.dat').write_text('0.5\n'*n);(work/'caps.dat').write_text('1\n'*n)
                deck=f'''&input_parameters
NOS={n},NODmin={f['lo']},NODmax={f['hi']},ALG={alg},
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={n},NO_two={len(f['bonds'])},wk_dim=-1,MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=1,re_wf=0,cal_lm=0,cal_cf=0,NO_two_cf=0,FILE_two_cf='unused'
/
&input_lancz
lnc_ene_conv=1d-14,minitr={f.get('minitr',5)},maxitr={f.get('maxiter',1000)},itrint=1
/
&input_TRLan
NOE={wanted},NOK={min(4,dim)},NOM={min(8,dim)},maxitr={f.get('maxiter',1000)},lnc_ene_conv=1d-13,
i_vec_min=1,i_vec_max={wanted}
/
'''
                (work/'input.dat').write_text(deck)
                result={'name':name,'alg':alg,'threads':threads,'expected_failure':f.get('expected_failure',False),
                        'input_sha256':hashlib.sha256(deck.encode()).hexdigest()}
                try:
                    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
                    if f['hi']>0:
                        with (work/'input.dat').open('rb') as inp,(work/'operator.log').open('wb') as log:
                            r=subprocess.run([str(operator)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=120)
                        require(r.returncode==0,'basis exporter failed')
                        basis=np.loadtxt(work/'basis.dat',dtype=int,ndmin=2)[:,1:]
                        indices=[sum(1<<(n-int(site)) for site in row if site!=0) for row in basis]
                    else:indices=[0]
                    actual_h=h[np.ix_(indices,indices)]
                    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
                        run=subprocess.run([str(exe)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=120)
                    result['exit_code']=run.returncode
                    if f.get('expected_failure'):
                        require(run.returncode!=0,'insufficient iteration budget returned success')
                        text=(work/'stdout.log').read_text()
                        require('maximum' in text and 'did not converge' in text,'nonconvergence not diagnosed')
                        require(not (work/'eigenvalues.dat').exists() and not list(work.glob('wf*.bin')),'nonconverged CLI wrote results')
                    else:
                        require(run.returncode==0,'CLI solve failed')
                        values=np.loadtxt(work/'eigenvalues.dat',ndmin=2)[:,1]
                        np.testing.assert_allclose(values,reference[:wanted],atol=1e-9,rtol=1e-9)
                        vectors=[]
                        for index in range(1,wanted+1):
                            data=(work/f'wf{index}.bin').read_bytes()
                            # GNU native unformatted single-record complex(8); this
                            # is a validation reader, not a portable persistence API.
                            length=int.from_bytes(data[:4],sys.byteorder,signed=True)
                            require(length==16*dim and len(data)==length+8,'unexpected wavefunction record length')
                            require(data[:4]==data[-4:],'record markers disagree')
                            vectors.append(np.frombuffer(data[4:-4],dtype=np.complex128).copy())
                        q=np.array(vectors).T
                        residual=np.linalg.norm(actual_h@q-q*values,axis=0)
                        scale=max(np.linalg.norm(actual_h,2),np.finfo(float).tiny)
                        normalized=residual/scale
                        require(np.isfinite(residual).all() and max(normalized)<1e-6,'invalid actual CLI eigenpair residual')
                        orth=float(np.linalg.norm(q.conj().T@q-np.eye(wanted),2))
                        require(orth<1e-10,'CLI wavefunctions not orthonormal')
                        require('NaN' not in (work/'stdout.log').read_text(),'nonfinite CLI quality diagnostic')
                        result.update(eigenvalues=values.tolist(),residuals=residual.tolist(),normalized_residuals=normalized.tolist(),orthogonality_error=orth)
                    result['status']='passed'
                except (OSError,ValueError,AssertionError,subprocess.TimeoutExpired) as error:
                    result.update(status='failed',error=str(error))
                print(result['status'].upper()+': '+name,flush=True);results.append(result)
    report={'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'cases':results,
            'passed':sum(r['status']=='passed' for r in results),'failed':sum(r['status']=='failed' for r in results)}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(report['passed'],'passed;',report['failed'],'failed; evidence:',out)
    return int(report['failed']!=0)


if __name__=='__main__':raise SystemExit(main())
