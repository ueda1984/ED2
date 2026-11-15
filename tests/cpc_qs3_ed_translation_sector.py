#!/usr/bin/env python3
"""Compare every L=6, Ndown=2, periodic XXZ translation sector:
published original QS-Cube/ED source/QS3 versus generalized ED2,
and a third independently constructed momentum-projector spectrum.
No claim of runtime speedup, binary ABI or legacy DSF support.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import numpy as np
from cpc_qs3_ed_shared_xxz_crosscheck import ORIGINAL_ED_SHA,JXY,JZ,read_energies

L=6
N_DOWN=2
BONDS=[(a,a+1) for a in range(1,L)]+[(1,L)]
TOL=1e-9

def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""):
            h.update(c)
    return h.hexdigest()

def independent_reference():
    states=[x for x in range(1<<L) if x.bit_count()==N_DOWN]
    indices={mask:i for i,mask in enumerate(states)}
    assert len(states)==15
    H=np.zeros((15,15),dtype=complex)
    T=np.zeros((15,15),dtype=complex)
    for j,mask in enumerate(states):
        rolled=((mask<<1)&((1<<L)-1))|(mask>>(L-1))
        T[indices[rolled],j]=1
        for a,b in BONDS:
            ia=(mask>>(a-1))&1
            ib=(mask>>(b-1))&1
            H[j,j]+=JZ*(.5-ia)*(.5-ib)
            if ia!=ib:
                moved=mask^(1<<(a-1))^(1<<(b-1))
                H[indices[moved],j]+=.5*JXY
    if max(np.max(np.abs(H-H.conj().T)),np.max(np.abs(H@T-T@H)))>1e-12:
        raise AssertionError("Independently constructed operator invalid")
    result={}
    for k in range(L):
        P=np.zeros_like(H)
        Tr=np.eye(15,dtype=complex)
        for r in range(L):
            P+=np.exp(-2j*np.pi*k*r/L)*Tr/L
            Tr=Tr@T
        if max(np.max(np.abs(P-P.conj().T)),np.max(np.abs(P@P-P)))>1e-12:
            raise AssertionError("Projection algebra invalid")
        vals,vecs=np.linalg.eigh((P+P.conj().T)/2)
        Q=vecs[:,vals>.5]
        if Q.shape[1]<1:raise AssertionError("Empty momentum sector")
        projected=Q.conj().T@H@Q
        result[k]=np.linalg.eigvalsh((projected+projected.conj().T)/2).tolist()
    flat=sorted(v for spectrum in result.values() for v in spectrum)
    if len(flat)!=15 or max(abs(x-y) for x,y in zip(flat,np.linalg.eigvalsh(H)))>1e-11:
        raise AssertionError("Projected spectrum not complete")
    return result

def prepare(root,k,old):
    root.mkdir(parents=True,exist_ok=False)
    if old:
        (root/"bonds.dat").write_text("".join(
            f"{a} {b} {JXY:.17g} {JZ:.17g}\n" for a,b in BONDS))
        deck=f"""&input_parameters
NOS=6,NOD=2,LX=6,LY=1,LZ=1,KX={k},KY=0,KZ=0,NOxxz=6,
ALG=3,cal_lm=0,cal_cf=0,cal_dsf=0,wr_wf=0,re_wf=0,
FILExxz='bonds.dat',FILEwf='./',OUTDIR='./'
&end
&input_static
NOV=1,NOLM=0,NOCF=0,FILElm='unused',FILECF='unused'
&end
&input_dynamic
spsmsz=3,itr_dsf=0,qx=0.,qy=0.,qz=0.,rfield=0.,FILEpos='unused'
&end
"""
    else:
        (root/"bonds.dat").write_text("".join(
            f"{a} {b} {JXY:.17g} {JXY:.17g} {JZ:.17g} 0 0 0 0 0 0\n" for a,b in BONDS))
        (root/"fields.dat").write_text("".join(f"{i} 0 0 0\n" for i in range(1,7)))
        (root/"spins.dat").write_text("0.5\n"*L)
        (root/"caps.dat").write_text("1\n"*L)
        (root/"translation.dat").write_text("".join(f"{(i%L)+1}\n" for i in range(1,L+1)))
        deck=f"""&input_parameters
NOS=6,NODmin=2,NODmax=2,ALG=3,
L1=6,L2=1,L3=1,L4=1,L5=1,L6=1,
M1={k},M2=0,M3=0,M4=0,M5=0,M6=0,
FILE_S1='translation.dat',
NO_one=6,NO_two=6,wk_dim=0,MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,
cal_lm=0,cal_cf=0
/
&input_lancz
lnc_ene_conv=1d-12,minitr=5,maxitr=2000,itrint=5
/
&input_TRLan
NOE=1,NOK=2,NOM=4,maxitr=2000,lnc_ene_conv=1d-12,
i_vec_min=1,i_vec_max=1
/
"""
    (root/"input.dat").write_text(deck)

def run(exe,work,expected,old,timeout):
    env=dict(os.environ,OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",
        OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as stream:
        r=subprocess.run([str(exe)],stdin=stream,cwd=work,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
            timeout=timeout,env=env,check=False)
    (work/"stdout.log").write_bytes(r.stdout)
    stdout=r.stdout.decode(errors="replace")
    if r.returncode:
        raise AssertionError(f"Solver exit {r.returncode}: {stdout[-1600:]}")
    dims=re.findall(r"THS\(k\)\s*=\s*(\d+)",stdout,re.I)
    if len(dims)!=1 or int(dims[0])!=expected:
        raise AssertionError("Bad reduced dimension "+str(dims)+" "+stdout[-1100:])
    spectrum=read_energies(work/("full_eigenvalues.dat" if old else "eigenvalues.dat"),
                            expected,1 if old else 2)
    return dict(levels=spectrum,dimension=expected,input_sha256=digest(work/"input.dat"))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old-exe",type=Path,required=True)
    p.add_argument("--new-exe",type=Path,required=True)
    p.add_argument("--old-commit",required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--timeout",type=int,default=120)
    a=p.parse_args()
    if a.old_commit!=ORIGINAL_ED_SHA:p.error("Old source SHA not pinned")
    if not 10<=a.timeout<=600:p.error("Invalid timeout")
    old=a.old_exe.resolve(strict=True)
    new=a.new_exe.resolve(strict=True)
    out=a.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    ref=independent_reference()
    rows=[];failures=[]
    for k in range(L):
        root=out/f"k{k}"
        root.mkdir(parents=True,exist_ok=False)
        po=root/"original"
        pn=root/"ed2"
        try:
            prepare(po,k,True)
            prepare(pn,k,False)
            old_values=run(old,po,len(ref[k]),True,a.timeout)
            new_values=run(new,pn,len(ref[k]),False,a.timeout)
            maxerror=dict(
                old_reference=max(abs(x-y) for x,y in zip(old_values["levels"],ref[k])),
                ed2_reference=max(abs(x-y) for x,y in zip(new_values["levels"],ref[k])),
                old_new=max(abs(x-y) for x,y in zip(old_values["levels"],new_values["levels"])))
            if max(maxerror.values())>TOL:
                raise AssertionError("Spectrum disagreement "+str(maxerror))
            rows.append(dict(momentum_index=k,dimension=len(ref[k]),
                independent=ref[k],original=old_values,ed2=new_values,errors=maxerror))
            print("QS3-ED / ED2 MOMENTUM",k,"dim",len(ref[k]),"errors",maxerror,"PASS",flush=True)
        except Exception as ex:
            f=dict(momentum=k,error_type=type(ex).__name__,message=str(ex))
            failures.append(f)
            print("QS3-ED / ED2 MOMENTUM",k,"FAIL",f,flush=True)
    good=len(rows)==6 and not failures
    result=dict(schema=1,pass_complete=good,
        original_source_sha=a.old_commit,
        new_source_sha=os.environ.get("GITHUB_SHA","record source separately"),
        expected_sectors=6,verified_sectors=len(rows),
        expected_levels=15,verified_levels=sum(x["dimension"] for x in rows),
        original_binary_sha256=digest(old),ed2_binary_sha256=digest(new),
        rows=rows,failures=failures,
        caution="Original and ED2 momentum-sector energy continuity, not dynamic output or performance")
    (out/"report.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("QS3->ED2 PERIODIC MOMENTUM",len(rows),"/6,",
        result["verified_levels"],"/15 levels; success",good,flush=True)
    return 0 if good else 1

if __name__=="__main__":
    raise SystemExit(main())
