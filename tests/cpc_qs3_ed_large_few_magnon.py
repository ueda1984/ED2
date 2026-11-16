#!/usr/bin/env python3
"""Published QS3 and ED2 identical few-magnon momentum-sector checks.

Independent Burnside sector-count and one-magnon dispersion checks.
Do not interpret shared-runner timings as speedup benchmarks.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess

OLD_SHA = "a3b09934bba5f966bfb564718528de8778c2b46c"
SMOKE = [
 ((12,1,1),2,(0,0,0)), ((24,1,1),3,(0,0,0)),
 ((48,1,1),2,(1,0,0)), ((96,1,1),2,(0,0,0)),
 ((3,3,3),2,(0,0,0)), ((4,4,4),2,(0,0,0)),
 ((10,10,10),1,(5,5,5)),
]
STRESS = [((10,10,10),2,(0,0,0))]

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def sector_dim(lengths,n,k):
    N=math.prod(lengths)
    value=0j
    for x in range(lengths[0]):
      for y in range(lengths[1]):
        for z in range(lengths[2]):
          shifts=(x,y,z)
          order=1
          for d,L in zip(shifts,lengths):
            order=math.lcm(order,L//math.gcd(L,d))
          if n%order:continue
          phase=2*math.pi*sum(d*q/L for d,q,L in zip(shifts,k,lengths))
          value+=math.comb(N//order,n//order)*complex(math.cos(phase),-math.sin(phase))
    v=value/N
    if abs(v.imag)>1e-6 or abs(v.real-round(v.real))>1e-6 or v.real<1:
        raise ValueError(f"Invalid Burnside dimension {v}")
    return round(v.real)

def geometry(lengths):
    nx,ny,nz=lengths
    index=lambda x,y,z:x+nx*(y+ny*z)+1
    bonds=[];translations=[[],[],[]]
    for z in range(nz):
      for y in range(ny):
        for x in range(nx):
          for d,L in enumerate(lengths):
            coords=[x,y,z]
            coords[d]=(coords[d]+1)%L
            t=index(*coords)
            translations[d].append(t)
            if L>1:bonds.append((index(x,y,z),t))
    return bonds,translations

def create(root,lengths,n,k,old):
    root.mkdir(parents=True,exist_ok=False)
    N=math.prod(lengths)
    bonds,translations=geometry(lengths)
    xy,zz=(1.,1.) if lengths[1]>1 else (-1.,.8)
    if old:
        (root/"bonds.dat").write_text("".join(f"{i} {j} {xy:g} {zz:g}\n" for i,j in bonds))
        params=(f"NOS={N},NOD={n},LX={lengths[0]},LY={lengths[1]},LZ={lengths[2]},"
                f"KX={k[0]},KY={k[1]},KZ={k[2]},NOxxz={len(bonds)},"
                f"ALG={3 if n==1 else 1},cal_lm=0,cal_cf=0,cal_dsf=0,"
                "wr_wf=0,re_wf=0,FILExxz='bonds.dat',FILEwf='./',OUTDIR='./'")
        extra=("&input_static\nNOV=1,NOLM=0,NOCF=0,FILElm='unused',FILECF='unused'\n/\n"
               "&input_dynamic\nspsmsz=3,itr_dsf=0,QX=0.,QY=0.,QZ=0.,rfield=0.,FILEpos='unused'\n/\n")
    else:
        (root/"bonds.dat").write_text(
          "".join(f"{i} {j} {xy:g} {xy:g} {zz:g} 0 0 0 0 0 0\n" for i,j in bonds))
        (root/"spins.dat").write_text("0.5\n"*N)
        (root/"caps.dat").write_text("1\n"*N)
        (root/"fields.dat").write_text("".join(f"{i} 0 0 0\n" for i in range(1,N+1)))
        for d,trans in enumerate(translations):
            if lengths[d]>1:
                (root/f"shift{d+1}.dat").write_text("".join(f"{j}\n" for j in trans))
        params=(f"NOS={N},NODmin={n},NODmax={n},"
                f"L1={lengths[0]},L2={lengths[1]},L3={lengths[2]},L4=1,L5=1,L6=1,"
                f"M1={k[0]},M2={k[1]},M3={k[2]},M4=0,M5=0,M6=0,"
                f"ALG={3 if n==1 else 1},NO_one={N},NO_two={len(bonds)},"
                "wk_dim=0,MNTE=-1,FILE_xyz_dm_gamma='bonds.dat',"
                "FILE_hvec='fields.dat',FILE_spin='spins.dat',FILE_NODmax='caps.dat',"
                "FILE_S1='shift1.dat',FILE_S2='shift2.dat',FILE_S3='shift3.dat',"
                "OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0")
        extra=""
    (root/"input.dat").write_text(
       f"&input_parameters\n{params}\n/\n{extra}"
       "&input_lancz\nlnc_ene_conv=1d-10,minitr=5,maxitr=1500,itrint=5\n/\n"
       "&input_TRLan\nNOE=2,NOK=3,NOM=8,maxitr=1500,"
       "lnc_ene_conv=1d-10,i_vec_min=1,i_vec_max=1\n/\n")
    return bonds,xy,zz

def run(exe,root,dim,old,timeout):
    env=dict(os.environ,OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",OPENBLAS_NUM_THREADS="1")
    timing=root/"time.txt"
    with (root/"input.dat").open("rb") as fd:
        try:
            cp=subprocess.run(["/usr/bin/time","-f","%e %M","-o",str(timing),str(exe)],
                              stdin=fd,cwd=root,env=env,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,timeout=timeout,check=False)
        except subprocess.TimeoutExpired as e:
            (root/"stdout.log").write_bytes(e.stdout or b"")
            raise TimeoutError(f"Exceeded {timeout}s")
    (root/"stdout.log").write_bytes(cp.stdout)
    stdout=cp.stdout.decode(errors="replace")
    if cp.returncode:raise RuntimeError(f"Exit {cp.returncode}: {stdout[-1400:]}")
    dims=re.findall(r"THS\(k\)\s*=\s*(\d+)",stdout,re.I)
    if len(dims)!=1 or int(dims[0])!=dim:
        raise RuntimeError(f"Wrong symmetry dimension {dims}, expected {dim}: {stdout[-1000:]}")
    if old and dim==1:
        lines=(root/"full_eigenvalues.dat").read_text().splitlines()
        if len(lines)!=1:raise RuntimeError("Original dense file invalid")
        energy=float(lines[0].split()[-1].replace("D","E"))
    else:
        es=re.findall(r"E_0\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)[EeDd][-+]?\d+)",stdout)
        if len(es)==1:energy=float(es[0].replace("D","E"))
        elif not old:
            lines=(root/"eigenvalues.dat").read_text().splitlines()
            energy=float(lines[0].split()[-1].replace("D","E"))
        else:raise RuntimeError("Missing old QS3 Lanczos energy")
    if not math.isfinite(energy):raise RuntimeError("Nonfinite eigenvalue")
    stat=timing.read_text().split()
    iterations=re.findall(r"total\s+lanczos\s+step:\s*(\d+)",stdout,re.I)
    return dict(energy=energy,momentum_dimension=dim,
                iterations=int(iterations[0]) if len(iterations)==1 else None,
                elapsed_s_uncontrolled=float(stat[0]),peak_rss_kib_uncontrolled=int(stat[1]),
                input_sha256=digest(root/"input.dat"))

def one_energy(lengths,k,bonds,xy,zz):
    n=math.prod(lengths)
    return len(bonds)*zz/4. - len(bonds)*zz/n + xy*sum(
        math.cos(2*math.pi*q/L) for q,L in zip(k,lengths) if L>1)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old",type=Path,required=True)
    p.add_argument("--new",type=Path,required=True)
    p.add_argument("--old-sha",required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--tier",choices=("smoke","stress"),default="smoke")
    p.add_argument("--timeout",type=int,default=120)
    a=p.parse_args()
    if a.old_sha!=OLD_SHA or not 30<=a.timeout<=600:
        p.error("Old source SHA not pinned or timeout invalid")
    old=a.old.resolve(strict=True);new=a.new.resolve(strict=True)
    base=a.output.resolve();base.mkdir(parents=True,exist_ok=False)
    cases=SMOKE if a.tier=="smoke" else STRESS
    accepted=[];failures=[]
    for lengths,n,k in cases:
        label="x".join(map(str,lengths))+f"-n{n}-k"+"-".join(map(str,k))
        work=base/label;work.mkdir(parents=True,exist_ok=False)
        try:
            dim=sector_dim(lengths,n,k)
            bonds,xy,zz=create(work/"old",lengths,n,k,True)
            create(work/"ed2",lengths,n,k,False)
            o=run(old,work/"old",dim,True,a.timeout)
            m=run(new,work/"ed2",dim,False,a.timeout)
            delta=abs(o["energy"]-m["energy"])
            if delta>2e-6:raise AssertionError(f"Old/ED2 energy error {delta}")
            independent=None
            if n==1:
                independent=one_energy(lengths,k,bonds,xy,zz)
                if max(abs(independent-o["energy"]),abs(independent-m["energy"]))>2e-8:
                    raise AssertionError("Analytic one-magnon energy mismatch")
            accepted.append(dict(label=label,sites=math.prod(lengths),n_down=n,k=k,
                                 full_sector_dim=math.comb(math.prod(lengths),n),
                                 momentum_dim=dim,old=o,ed2=m,energy_delta=delta,
                                 analytic_one_magnon=independent))
            print("QS3 LARGE",label,"dim",dim,"delta",delta,"PASS",flush=True)
        except Exception as e:
            err=dict(label=label,type=type(e).__name__,message=str(e))
            failures.append(err);print("QS3 LARGE",err,"FAIL",flush=True)
    passed=len(accepted)==len(cases) and not failures
    (base/"report.json").write_text(json.dumps(dict(
       passed=passed,tier=a.tier,old_source_sha=a.old_sha,
       ed2_source_sha=os.environ.get("GITHUB_SHA","unknown"),
       original_binary_sha256=digest(old),ed2_binary_sha256=digest(new),
       accepted=accepted,failures=failures,
       note="Shared CI timings not controlled; unmatched Lanczos starts; "
            "no 1000-site three-magnon execution or QuSpin speedup implied"),
       indent=2)+"\n")
    print("QS3 LARGE",a.tier,len(accepted),"/",len(cases),"accepted",passed,flush=True)
    return 0 if passed else 1

if __name__=="__main__":raise SystemExit(main())
