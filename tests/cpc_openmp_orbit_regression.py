#!/usr/bin/env python3
"""Independent basis-count and cross-thread determinism gate for the OpenMP patch.

The 379 legacy fixtures are tested at every requested thread count. Two known
zero-cap defects are checked against independent product-state enumeration,
not against incorrect legacy output. Larger orbit/bit-word stress fixtures
are repeated to vary thread interleavings. This is not a performance test.
"""
import argparse
import cmath
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import struct
import subprocess

PROBE = '''program basis_probe
  use input_param, only: read_ip,NOS,NODmax,NODmin,local_NODmax,THS,list_s,list_r
  use state_lists, only: get_cum_ltd_rep_combination,cum_ltd_rep_combination,allocate_lists_omp_SQ
  implicit none
  integer :: rawdim,u
  call read_ip
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  rawdim=cum_ltd_rep_combination(NOS,NODmax)
  THS=rawdim
  call allocate_lists_omp_SQ
  open(newunit=u,file='basis.bin',form='unformatted',access='stream',status='new')
  write(u) rawdim,THS
  write(u) list_s
  write(u) list_r
  close(u)
end program
'''


def fixtures():
    out=[]
    def add(label,caps,lo,hi,gens,spins=None):
        out.append((label,caps,spins or [x/2 for x in caps],lo,hi,gens))
    for L in (4,6,8,9,10,12):
        T=[(i+1)%L for i in range(L)]
        for n in sorted(set([0,1,2,3,min(5,L)])):
            for k in range(L):add(f'ring{L}-n{n}-k{k}',[1]*L,n,n,[(T,L,k)])
    for L in (4,6,8):
        T=[(i+2)%L for i in range(L)]
        for lo,hi in ((0,2),(1,3),(2,2),(0,6)):
            for k in range(L//2):add(f'mixed{L}-{lo}-{hi}-k{k}',[1,2]*(L//2),lo,hi,[(T,L//2,k)])
    for nx,ny in ((3,3),(3,4),(4,4)):
        n=nx*ny
        tx=[(i%nx+1)%nx+nx*(i//nx) for i in range(n)]
        ty=[i%nx+nx*((i//nx+1)%ny) for i in range(n)]
        for nd in (2,3):
            for kx in range(nx):
                for ky in range(ny):add(f'plane{nx}x{ny}-n{nd}-k{kx}-{ky}',[1]*n,nd,nd,[(tx,nx,kx),(ty,ny,ky)])
    for L in (6,8,12):
        T=[(i+1)%L for i in range(L)]
        T2=[(i+2)%L for i in range(L)]
        R=[(-i)%L for i in range(L)]
        for k in (0,L//2):
            for p in (0,1):add(f'dihedral{L}-{k}-{p}',[1]*L,2,2,[(T,L,k),(R,2,p)])
        for k in (0,1):add(f'duplicate{L}-{k}',[1]*L,3,3,[(T,L,k),(T2,L//2,2*k%(L//2))])
    for n in (4,6,8,10):add(f'identity{n}',[1,2]*(n//2),0,3*n//2,[])
    for n in (4,6):
        add(f'caps-zero{n}',[0,2]*(n//2),0,4,[([(i+2)%n for i in range(n)],n//2,0)],spins=[.5,1]*(n//2))
    assert len(out)==379
    return out


def stress_fixtures():
    out=[]
    for L in (4,6):
        n=L**3
        tx=[(i%L+1)%L+L*(i//L) for i in range(n)]
        ty=[i%L+L*((i//L+1)%L)+L*L*(i//(L*L)) for i in range(n)]
        tz=[i%(L*L)+L*L*((i//(L*L)+1)%L) for i in range(n)]
        for k in (0,1):
            out.append((f'cube{L}-n3-k{k}',[1]*n,[.5]*n,3,3,[(tx,L,k),(ty,L,0),(tz,L,0)]))
    # Many distinct orbits update different bits in the same bitmap word.
    n=24
    out.append(('shared-word-n24-n4',[1]*n,[.5]*n,4,4,[([(i+12)%n for i in range(n)],2,0)]))
    return out


def expected_dimension(caps,lo,hi,gens):
    n=len(caps);total=0j;group=math.prod(L for _,L,k in gens)
    for exps in itertools.product(*(range(L) for _,L,k in gens)):
        perm=list(range(n));phase=0.
        for (m,L,k),e in zip(gens,exps):
            for _ in range(e):perm=[m[a] for a in perm]
            phase+=k*e/L
        seen=set();poly=[1]+[0]*hi
        for a in range(n):
            if a in seen:continue
            cycle=[];b=a
            while b not in seen:seen.add(b);cycle.append(b);b=perm[b]
            cap=caps[a]
            if any(caps[x]!=cap for x in cycle):raise AssertionError('Noninvariant fixture caps')
            d=len(cycle);new=[0]*(hi+1)
            for j,val in enumerate(poly):
                if val:
                    for q in range(cap+1):
                        if j+d*q<=hi:new[j+d*q]+=val
            poly=new
        total+=sum(poly[lo:hi+1])*cmath.exp(-2j*math.pi*phase)
    value=total/group
    if abs(value.imag)>1e-7 or abs(value.real-round(value.real))>1e-7:
        raise AssertionError('Nonintegral independent dimension')
    return round(value.real)


def zero_cap_reference(caps,lo,hi,gens):
    states=[]
    for occ in itertools.product(*(range(c+1) for c in caps)):
        if lo<=sum(occ)<=hi:
            labels=[0]*(hi-sum(occ))+[i+1 for i,c in enumerate(occ) for _ in range(c)]
            states.append((tuple(reversed(labels)),tuple(occ)))
    states.sort();occ=[x[1] for x in states];idx={s:j+1 for j,s in enumerate(occ)}
    elements=[]
    for exps in itertools.product(*(range(L) for _,L,k in gens)):
        perm=list(range(len(caps)));phase=0.
        for (m,L,k),e in zip(gens,exps):
            for _ in range(e):perm=[m[a] for a in perm]
            phase+=k*e/L
        elements.append((perm,cmath.exp(-2j*math.pi*phase)))
    ranks=[];norms=[]
    for rank,x in enumerate(occ,1):
        images=[];characters=[]
        for perm,phase in elements:
            mapped=[0]*len(caps)
            for old,new in enumerate(perm):mapped[new]=x[old]
            r=idx[tuple(mapped)];images.append(r)
            if r==rank:characters.append(phase)
        if min(images)!=rank:continue
        weight=abs(sum(characters))**2*len(elements)/len(characters)
        if weight>1e-15:ranks.append(rank);norms.append(math.sqrt(weight))
    return tuple(ranks),tuple(norms)


def make_deck(root,case):
    label,caps,spins,lo,hi,gens=case
    root.mkdir(parents=True,exist_ok=False);n=len(caps)
    (root/'caps.dat').write_text(''.join(f'{x}\n' for x in caps))
    (root/'spins.dat').write_text(''.join(f'{x}\n' for x in spins))
    (root/'bonds.dat').write_text('');(root/'fields.dat').write_text('')
    lines=[f'NOS={n},NODmin={lo},NODmax={hi},ALG=3,NO_one=0,NO_two=0,wk_dim=0,MNTE=-1,',
      "FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',FILE_spin='spins.dat',FILE_NODmax='caps.dat',",
      "OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0,"]
    for d in range(6):
        if d<len(gens):
            m,L,k=gens[d]
            (root/f't{d+1}.dat').write_text(''.join(f'{x+1}\n' for x in m))
            lines.append(f"L{d+1}={L},M{d+1}={k},FILE_S{d+1}='t{d+1}.dat',")
        else:lines.append(f'L{d+1}=1,M{d+1}=0,')
    (root/'input.dat').write_text('&input_parameters\n'+'\n'.join(lines)+'\n/\n')


def run(probe,work,case,threads):
    make_deck(work,case)
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1')
    with (work/'input.dat').open('rb') as inp,(work/'stdout.log').open('wb') as log:
        cp=subprocess.run([str(probe)],cwd=work,stdin=inp,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=90)
    text=(work/'stdout.log').read_text()
    if cp.returncode:raise RuntimeError(text[-2000:])
    f=work/'basis.bin'
    if not f.exists():
        if 'THS = 0' not in text:raise AssertionError('Missing basis output')
        return None
    payload=f.read_bytes();raw,dim=struct.unpack_from('=ii',payload)
    if dim<=0 or len(payload)!=8+12*dim:raise AssertionError('Bad binary layout')
    ranks=struct.unpack_from('='+str(dim)+'i',payload,8)
    norms=struct.unpack_from('='+str(dim)+'d',payload,8+4*dim)
    if any(a>=b for a,b in zip(ranks,ranks[1:])):raise AssertionError('Nonordered ranks')
    return {'raw':raw,'dim':dim,'ranks':ranks,'norms':norms,'sha256':hashlib.sha256(payload).hexdigest(),
            'parallel_path':'Basis builder: atomic-orbits' in text}


def compile_probe(source,output,fc,checked):
    output.mkdir(parents=True,exist_ok=False)
    (output/'probe.f90').write_text(PROBE)
    flags=['-O1','-g','-fopenmp','-fcheck=all','-fsanitize=undefined','-fno-sanitize-recover=undefined'] if checked else ['-O2','-fopenmp']
    with (output/'build.log').open('w') as log:
        for name in ('input_param.f90','state_lists.f90'):
            subprocess.run([fc,*flags,'-c',str(source/'source'/name)],cwd=output,check=True,stdout=log,stderr=subprocess.STDOUT)
        subprocess.run([fc,*flags,'probe.f90','input_param.o','state_lists.o','-o','probe'],cwd=output,check=True,stdout=log,stderr=subprocess.STDOUT)
    return output/'probe'


def compare(reference,candidate,expected,case):
    if candidate is None:
        if expected!=0:raise AssertionError('Unexpected empty sector')
        return
    if candidate['dim']!=expected:raise AssertionError('Wrong independent dimension')
    if case[0].startswith('caps-zero'):
        ranks,norms=zero_cap_reference(case[1],case[3],case[4],case[5])
    else:
        if reference is None or reference['raw']!=candidate['raw']:raise AssertionError('Legacy dimension mismatch')
        ranks,norms=reference['ranks'],reference['norms']
    if ranks!=candidate['ranks'] or max((abs(x-y) for x,y in zip(norms,candidate['norms'])),default=0)>1e-11:
        raise AssertionError('Representative ranks or normalization differ')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',required=True,type=Path);p.add_argument('--candidate',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--compiler',default='gfortran')
    p.add_argument('--threads',default='1,2,4');p.add_argument('--stress-repeats',default=4,type=int)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    threads=tuple(int(x) for x in a.threads.split(','))
    if not threads or any(t<1 for t in threads) or len(set(threads))!=len(threads):p.error('Invalid threads')
    baseline=compile_probe(a.baseline.resolve(),root/'build-baseline',a.compiler,False)
    candidate=compile_probe(a.candidate.resolve(),root/'build-candidate',a.compiler,True)
    records=[]
    for stress,cases in ((False,fixtures()),(True,stress_fixtures())):
        for case in cases:
            label=case[0];expect=expected_dimension(case[1],case[3],case[4],case[5])
            ref=run(baseline,root/'baseline'/label,case,1);first=None
            for t in threads:
                for rep in range(a.stress_repeats if stress else 1):
                    got=run(candidate,root/'candidate'/label/f'omp{t}-rep{rep}',case,t)
                    compare(ref,got,expect,case)
                    digest=got['sha256'] if got else 'empty'
                    if first is None:first=digest
                    if digest!=first:raise AssertionError('Cross-thread/repeat basis nondeterminism')
                    records.append({'case':label,'threads':t,'repeat':rep,'stress':stress,'expected':expect,
                                    'sha256':digest,'parallel_path':got['parallel_path'] if got else None})
            if len(records)%100< len(threads):print('Verified',len(records),'candidate executions',flush=True)
    report={'passed':True,'fixtures':379,'thread_counts':threads,'stress_fixtures':len(stress_fixtures()),
            'stress_repeats':a.stress_repeats,'candidate_executions':len(records),'records':records,
            'note':'Oversubscribed thread counts may test correctness but do not establish scaling.'}
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('OPENMP BASIS REGRESSION',len(records),'candidate executions PASS',flush=True)

if __name__=='__main__':main()
