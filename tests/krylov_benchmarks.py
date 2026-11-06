#!/usr/bin/env python3
"""Serial spin-system benchmarks of Phase 1-B and 1-C, with independent residuals.

Builds instrumented plain copies outside the checkout; immutable numerical sources
are identified by Git/archive hashes and the exact instrumentation diff is saved.
Optional NumPy, compiler, git, BLAS/LAPACK; Linux wait4 measures each child's RSS.
"""
import argparse
import difflib
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import time
import numpy as np

REVISIONS={'phase1b':'2050208999940044c155cce7fb60d42b7228481d',
           'phase1c':'6d3480c5fb8ee7a67f2ac4017156ee4368bdd1d2'}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(repo,output,fc,flags):
    bins={}
    for label,rev in REVISIONS.items():
        stage=output/label;stage.mkdir()
        archive=subprocess.check_output(['git','archive',rev],cwd=repo)
        with tarfile.open(fileobj=io.BytesIO(archive)) as t:t.extractall(stage,filter='data')
        originals={f:(stage/'source'/f).read_text() for f in ['ham2vec.f90','main.f90']}
        s=originals['ham2vec.f90'];s=s.replace('module ham2vec\n','module ham2vec\n  use phase1d_probe, only: probe_h_calls\n',1)
        start=s.index('  subroutine ham_to_vec_wave_vector(');end=s.index('  end subroutine ham_to_vec_wave_vector',start)
        body=s[start:end];needle='    allocate(ni(NODmax),st_list(NODmax))';assert body.count(needle)==1
        body=body.replace(needle,'    probe_h_calls=probe_h_calls+1\n'+needle)
        s=s[:start]+body+s[end:];(stage/'source/ham2vec.f90').write_text(s)
        s=originals['main.f90'].replace('program main\n','program main\n  use phase1d_probe\n',1)
        assert s.count('  call read_ip\n')==1 and s.count('  stop\n')==1
        s=s.replace('  call read_ip\n','  call read_ip\n  call probe_begin()\n')
        # Export the actual basis before timing exits, for independent validation.
        # This disk output is included in all reported whole-process wall times.
        s=s.replace('allocate_lists_omp_SQ\n','allocate_lists_omp_SQ, list_fly\n',1)
        s=s.replace('local_NODmax, THS, ALG','local_NODmax, list_s, THS, ALG',1)
        s=s.replace('  !$ real(8) :: st, en','  integer :: probe_i\n  !$ real(8) :: st, en',1)
        s=s.replace('  stop\n',"  call probe_end()\n  open(93,file='benchmark-basis.dat',status='new')\n  do probe_i=1,THS\n    write(93,*) list_fly(list_s(probe_i),NODmax,NODmin,NOS)\n  end do\n  close(93)\n  stop\n")
        (stage/'source/main.f90').write_text(s)
        patch=''.join(''.join(difflib.unified_diff(originals[f].splitlines(True),(stage/'source'/f).read_text().splitlines(True),fromfile=f'a/source/{f}',tofile=f'b/source/{f}')) for f in originals)
        (output/f'{label}-instrumentation.patch').write_text(patch)
        src=stage/'source';probe=repo/'tests/krylov_benchmark_probe.f90'
        names=['input_param','state_lists','ham2vec','lanczos','eigen_solver','get_expectation_values','main']
        if label=='phase1c':names.insert(0,'krylov_solvers')
        names.insert(0,'blas_lapack_interfaces')
        with (output/f'{label}-build.log').open('w') as log:
            subprocess.run([fc,*flags.split(),'-c',str(probe)],cwd=src,stdout=log,stderr=subprocess.STDOUT,check=True)
            for name in names:subprocess.run([fc,*flags.split(),'-c',name+'.f90'],cwd=src,stdout=log,stderr=subprocess.STDOUT,check=True)
            subprocess.run([fc,*flags.split(),'-o','bench','krylov_benchmark_probe.o',*[n+'.o' for n in names],'-llapack','-lopenblas'],cwd=src,stdout=log,stderr=subprocess.STDOUT,check=True)
        bins[label]=src/'bench'
    return bins


def spin_model(n,kind):
    fields=[[.11+.013*i,-.07+.002*i,.17-.009*i] for i in range(n)]
    couplings=[1.,.83,.71,.09,-.04,.07,.02,-.03,.05]
    if kind=='j1j2':
        bonds=[(i,i+1,couplings) for i in range(n-1)]+[(i,i+2,[.35*c for c in couplings]) for i in range(n-2)]
    elif kind=='heisenberg':
        fields=[[0.,0.,.031*(i+1)] for i in range(n)]
        bonds=[(i,i+1,[1.,1.,1.,0.,0.,0.,0.,0.,0.]) for i in range(n-1)]
    else:bonds=[(i,i+1,couplings) for i in range(n-1)]
    return bonds,fields


def apply_spin(n,bonds,fields,x):
    """Independent Pauli bit actions, up-before-down tensor convention.

    This constructs transitions from spin matrices, not ED2 cache formulas.
    Memory O(2**n), including validation of the large benchmark vectors.
    """
    basis=np.arange(2**n);ops=[np.array([[0,1],[1,0]])/2,np.array([[0,-1j],[1j,0]])/2,np.diag([.5,-.5])]
    def action(v,site,axis):
        bit=1<<(n-site-1);down=(basis&bit)!=0
        if axis==2:return np.where(down,-.5,.5)*v
        # output at b receives matrix element <b|S|b xor bit>.
        return np.where(down,ops[axis][1,0],ops[axis][0,1])*v[basis^bit]
    y=np.zeros_like(x,dtype=complex)
    for i,f in enumerate(fields):
        for a,c in enumerate(f):
            if c:y+=c*action(x,i,a)
    for i,j,c in bonds:
        for a in range(3):y+=c[a]*action(action(x,j,a),i,a)
        for a,b,k in [(1,2,0),(2,0,1),(0,1,2)]:
            y+=c[3+k]*(action(action(x,j,b),i,a)-action(action(x,j,a),i,b))
            y+=c[6+k]*(action(action(x,j,b),i,a)+action(action(x,j,a),i,b))
    return y


def reference(n,bonds,fields):
    if n>8:return None
    eye=np.eye(2**n,dtype=complex)
    h=np.column_stack([apply_spin(n,bonds,fields,eye[:,i]) for i in range(2**n)])
    assert np.max(abs(h-h.conj().T))<1e-13
    return np.linalg.eigh(h)


def process(exe,work,env,timeout):
    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
        started=time.perf_counter();proc=subprocess.Popen([str(exe)],cwd=work,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT)
        while True:
            pid,status,usage=os.wait4(proc.pid,os.WNOHANG)
            if pid:break
            if time.perf_counter()-started>timeout:
                proc.kill();pid,status,usage=os.wait4(proc.pid,0);break
            time.sleep(.002)
        proc.returncode=os.waitstatus_to_exitcode(status)
    return proc.returncode,time.perf_counter()-started,usage.ru_maxrss


def run(exe,root,config,label,rep,ref):
    n=config['n'];alg=config['alg'];wanted=1 if alg==1 else config['wanted'];keep=config['keep'];cap=config['capacity']
    seed=[19,73,211][rep];threads=config['threads'];bonds,fields=spin_model(n,config['model'])
    ident=f"{config['name']}-{label}-omp{threads}-seed{seed}";work=root/ident;work.mkdir()
    (work/'fields.dat').write_text(''.join(f'{i+1} '+' '.join(map(str,f))+'\n' for i,f in enumerate(fields)))
    (work/'bonds.dat').write_text(''.join(f'{i+1} {j+1} '+' '.join(map(str,c))+'\n' for i,j,c in bonds))
    (work/'spins.dat').write_text('0.5\n'*n);(work/'caps.dat').write_text('1\n'*n)
    deck=f'''&input_parameters NOS={n},NODmin=0,NODmax={n},ALG={alg},
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={n},NO_two={len(bonds)},wk_dim={config['cache']},MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=1,re_wf=0,cal_lm=0,cal_cf=0,NO_two_cf=0,FILE_two_cf='unused'/
&input_lancz lnc_ene_conv={config['tol']},minitr=5,itrint=1,maxitr=400/
&input_TRLan NOE={wanted},NOK={keep},NOM={cap},lnc_ene_conv={config['tol']},maxitr=200,
i_vec_min=1,i_vec_max={wanted}/
'''
    (work/'input.dat').write_text(deck)
    r=dict(config,id=ident,revision=REVISIONS[label],label=label,seed=seed,replicate=rep)
    r['input_sha256']={p.name:sha(p) for p in work.glob('*.dat')}
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',ED2_BENCH_SEED=str(seed),OMP_PROC_BIND='close',OMP_PLACES='cores')
    code,wall,rss=process(exe,work,env,120);r.update(exit_code=code,wall_seconds=wall,wait4_peak_rss_kib=rss)
    text=(work/'stdout.log').read_text();r['restart_cycles']=len(re.findall('tritr, ene',text));r['rng_state']=(work/'benchmark-rng.dat').read_text().splitlines()
    if code!=0 or not (work/'eigenvalues.dat').exists():
        r['classification']='solver_failure';r['diagnostic']=text[-1500:];return r
    metrics=(work/'benchmark-metrics.dat').read_text().split();r.update(h_applications=int(metrics[0]),compute_seconds=float(metrics[1]),peak_rss_kib=int(metrics[2]))
    values=np.loadtxt(work/'eigenvalues.dat',ndmin=2)[:,1];basis=np.loadtxt(work/'benchmark-basis.dat',dtype=int,ndmin=2)
    indices=[sum(1<<(n-int(site)) for site in row if site!=0) for row in basis]
    assert len(set(indices))==2**n
    q=[]
    for i in range(wanted):
        b=(work/f'wf{i+1}.bin').read_bytes();assert int.from_bytes(b[:4],sys.byteorder,signed=True)==16*2**n and b[:4]==b[-4:]
        v=np.frombuffer(b[4:-4],dtype=np.complex128);z=np.zeros(2**n,complex);z[indices]=v;q.append(z)
    q=np.array(q).T;hq=np.column_stack([apply_spin(n,bonds,fields,q[:,i]) for i in range(wanted)])
    residuals=np.linalg.norm(hq-q*values,axis=0);orth=float(np.linalg.norm(q.conj().T@q-np.eye(wanted),2))
    scale_bound=sum(sum(abs(x) for x in f)*.5 for f in fields)+sum((sum(map(abs,c[:3]))+2*sum(map(abs,c[3:])))*.25 for _,_,c in bonds)
    r.update(eigenvalues=values.tolist(),residual_norms=residuals.tolist(),operator_norm_upper_bound=scale_bound,orthogonality_error=orth)
    dim=2**n;cached=dim if config['cache']<0 else min(config['cache'],dim);slots=2*n+8*len(bonds)+1
    r['memory_estimate_bytes']={'cache':20*slots*cached,'solver_leading':(192*dim+112*min(400,dim) if alg==1 and label=='phase1c' else 64*dim+112*min(400,dim) if alg==1 else 16*dim*(4*cap+wanted+2)+32*cap**2 if label=='phase1c' else 16*dim*(cap+keep+1)+40*cap**2)}
    r['classification']='valid_residual_only' if max(residuals)/scale_bound<1e-6 and orth<1e-8 else 'accuracy_warning'
    if ref is not None:
        e,u=ref;errors=abs(values-e[:wanted]);leak=np.linalg.norm(q-u[:,:wanted]@(u[:,:wanted].conj().T@q),2)
        r.update(eigenvalue_errors=errors.tolist(),projector_error=float(leak),reference_eigenvalues=e[:wanted].tolist())
        r['classification']='reference_agreement' if max(errors)<1e-9 and max(residuals)<1e-6 and orth<1e-8 else 'accuracy_warning'
    return r


def configs():
    for n in [6,8,10,12,14]:
        for alg in [1,2]:
            yield dict(name=f'chain-n{n}-alg{alg}',model='xyz-chain',n=n,alg=alg,wanted=2,keep=4,capacity=16,cache=-1,tol=1e-12)
    yield dict(name='chain-n14-alg2-c64',model='xyz-chain',n=14,alg=2,wanted=2,keep=4,capacity=64,cache=-1,tol=1e-12)
    for cap,keep,wanted in [(12,4,2),(32,4,2),(32,8,4)]:
        yield dict(name=f'j1j2-c{cap}-k{keep}-w{wanted}',model='j1j2',n=10,alg=2,wanted=wanted,keep=keep,capacity=cap,cache=-1,tol=1e-12)
    for cache in [0,128]:
        yield dict(name=f'heisenberg-cache{cache}',model='heisenberg',n=10,alg=2,wanted=2,keep=4,capacity=16,cache=cache,tol=1e-12)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--fc',default='gfortran');p.add_argument('--flags',default='-O2 -fopenmp');p.add_argument('--repeats',type=int,choices=[1,2,3],default=3)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);bins=build(a.repo.resolve(),out,a.fc,a.flags)
    results=[];refs={}
    metadata={'cpu':Path('/proc/cpuinfo').read_text().split('model name',1)[-1].split('\n',1)[0].strip(), 'load_start':os.getloadavg(), 'cpu_max':Path('/sys/fs/cgroup/cpu.max').read_text().strip() if Path('/sys/fs/cgroup/cpu.max').exists() else None, 'memory_max':Path('/sys/fs/cgroup/memory.max').read_text().strip() if Path('/sys/fs/cgroup/memory.max').exists() else None, 'ldd':subprocess.check_output(['ldd',str(bins['phase1c'])],text=True), 'numpy':np.__version__}
    for config in configs():
        key=(config['n'],config['model'])
        if key not in refs:refs[key]=reference(config['n'],*spin_model(*key))
        for threads in [1,2]:
            for rep in range(a.repeats):
                # Alternate version order to reduce monotonic host-load bias.
                for label in (['phase1b','phase1c'] if rep%2==0 else ['phase1c','phase1b']):
                    r=run(bins[label],out,dict(config,threads=threads),label,rep,refs[key]);results.append(r)
                    print(r['id'],r['classification'],r['wall_seconds'],flush=True)
        (out/'report.json').write_text(json.dumps(dict(compiler=subprocess.check_output([a.fc,'--version'],text=True).splitlines()[0],flags=a.flags,binary_sha256={k:sha(v) for k,v in bins.items()},probe_sha256=sha(a.repo/'tests/krylov_benchmark_probe.f90'),harness_sha256=sha(Path(__file__)),thread_environment=dict(OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores'),affinity=sorted(os.sched_getaffinity(0)),metadata=metadata,load_latest=os.getloadavg(),cases=results),indent=2,allow_nan=False)+'\n')
    return 0


if __name__=='__main__':raise SystemExit(main())
