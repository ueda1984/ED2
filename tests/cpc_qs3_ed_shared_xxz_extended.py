#!/usr/bin/env python3
"""Three pinned common-sector QS3-ED versus ED2 XXZ comparisons.

This *extends* the primary fixed 4-site/Ndown=2 historical-program
test to the legacy near-saturation Ndown=1 and half-filled Ndown=3
sectors of a 6-site open spin-half chain. Both codes are actual Fortran
executables. The separate product-basis matrix gives a third eigenvalue
reference. No altered original QS-Cube/ED source is used.

Tested cases: (size,Ndown) in [(4,2),(6,1),(6,3)];
dimensions: 6, 6, and 20 (32 total eigenvalues).
No near-saturation performance claim or dynamic structure factor.
"""
import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from cpc_qs3_ed_shared_xxz_crosscheck import (
    ORIGINAL_ED_SHA,JXY,JZ,run,digest,read_energies,
)

CASES=((4,2),(6,1),(6,3))
ERROR_TOL=1e-8


def independent(size,nod):
    states=[n for n in range(1<<size) if n.bit_count()==nod]
    index={state:i for i,state in enumerate(states)}
    dim=math.comb(size,nod)
    if len(states)!=dim: raise AssertionError("Independent basis wrong")
    H=np.zeros((dim,dim),dtype=np.float64)
    for col,state in enumerate(states):
        for i in range(size-1):
            a=(state>>i)&1;b=(state>>(i+1))&1
            H[col,col]+=JZ*(.5-a)*(.5-b)
            if a!=b:
                flipped=state^(1<<i)^(1<<(i+1))
                H[col,index[flipped]]+=.5*JXY
    if np.max(np.abs(H-H.T))>1e-14:
        raise AssertionError("Independent XXZ sector not Hermitian")
    return np.linalg.eigvalsh(H).tolist()


def prepare(root,size,nod,old):
    root.mkdir(parents=True,exist_ok=False)
    bonds=[(i,i+1) for i in range(1,size)]
    if old:
        (root/"bonds_old.dat").write_text("".join(
            f"{i} {j} {JXY:.17g} {JZ:.17g}\n" for i,j in bonds))
        deck=f"""&input_parameters
NOS={size}, NOD={nod}, NOxxz={len(bonds)}, ALG=3,
cal_lm=0, cal_cf=0, cal_dsf=0,
wr_wf=0, re_wf=0,
FILExxz='bonds_old.dat', FILEwf='./', OUTDIR='./'
&end
&input_static
NOV=1, NOLM=0, NOCF=0, FILElm='unused-lm.dat', FILECF='unused-cf.dat'
&end
&input_dynamic
spsmsz=3, itr_dsf=0, qx=0., qy=0., qz=0., rfield=0., FILEpos='unused-positions.dat'
&end
"""
    else:
        (root/"bonds_new.dat").write_text("".join(
            f"{i} {j} {JXY:.17g} {JXY:.17g} {JZ:.17g} 0 0 0 0 0 0\n"
            for i,j in bonds))
        (root/"fields.dat").write_text("".join(f"{i} 0 0 0\n" for i in range(1,size+1)))
        (root/"spins.dat").write_text("0.5\n"*size)
        (root/"caps.dat").write_text("1\n"*size)
        deck=f"""&input_parameters
NOS={size}, NODmin={nod}, NODmax={nod}, ALG=3,
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={size}, NO_two={len(bonds)}, wk_dim=0, MNTE=-1,
FILE_xyz_dm_gamma='bonds_new.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=0, cal_cf=0
/
&input_lancz
lnc_ene_conv=1d-12, minitr=4, maxitr=1000, itrint=5
/
&input_TRLan
NOE=2, NOK=2, NOM=4, maxitr=1000, lnc_ene_conv=1d-12,
i_vec_min=1, i_vec_max=1
/
"""
    (root/"input.dat").write_text(deck)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old-exe",type=Path,required=True)
    p.add_argument("--new-exe",type=Path,required=True)
    p.add_argument("--old-commit",required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--timeout",type=int,default=90)
    a=p.parse_args()
    if a.old_commit!=ORIGINAL_ED_SHA or not 5<=a.timeout<=600:
        p.error("Unreviewed legacy SHA or invalid timeout")
    old_exe=a.old_exe.resolve(strict=True)
    new_exe=a.new_exe.resolve(strict=True)
    output=a.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    reports=[]
    errors=[]
    for n,nod in CASES:
        case_root=output/f"size{n}-down{nod}"
        old=case_root/"original-qs3"
        new=case_root/"modern-ed2"
        case_root.mkdir(parents=True,exist_ok=False)
        try:
            prepare(old,n,nod,old=True)
            prepare(new,n,nod,old=False)
            run(old_exe,old,a.timeout)
            run(new_exe,new,a.timeout)
            dim=math.comb(n,nod)
            ref=independent(n,nod)
            legacy=read_energies(old/"full_eigenvalues.dat",dim,1)
            modern=read_energies(new/"eigenvalues.dat",dim,2)
            diffs=dict(old_vs_ref=max(abs(x-y) for x,y in zip(legacy,ref)),
                       ed2_vs_ref=max(abs(x-y) for x,y in zip(modern,ref)),
                       old_vs_ed2=max(abs(x-y) for x,y in zip(legacy,modern)))
            if max(diffs.values())>ERROR_TOL:
                raise AssertionError("Historical solver and ED2 do not reproduce same sector")
            record=dict(size=n,n_down=nod,hilbert_dimension=dim,Jxy=JXY,Jz=JZ,
                        reference_eigenvalues=ref,
                        original_qs3_eigenvalues=legacy,
                        ed2_eigenvalues=modern,
                        errors=diffs,pass_complete=True)
            reports.append(record)
            print("HISTORICAL QS3->ED2 XXZ",f"L={n}",f"Ndown={nod}",
                  f"levels={dim}",diffs,"PASSED",flush=True)
        except Exception as exc:
            error=dict(size=n,n_down=nod,error_type=type(exc).__name__,
                       message=str(exc))
            errors.append(error)
            print("HISTORICAL QS3->ED2 FAILURE",error,flush=True)
    passed=len(reports)==len(CASES) and not errors
    result=dict(schema=1,pass_complete=passed,
                original_repository="https://github.com/QS-Cube/ED",
                old_commit=a.old_commit,
                new_source_sha=os.environ.get("GITHUB_SHA","checkout SHA required"),
                old_exe_sha256=digest(old_exe),new_exe_sha256=digest(new_exe),
                expected_cases=len(CASES),accepted_cases=len(reports),
                expected_levels=sum(math.comb(n,k) for n,k in CASES),
                verified_levels=sum(r["hilbert_dimension"] for r in reports),
                records=reports,failures=errors,
                caution="Common XXZ physical sectors only; does not test old dynamic structure factor or speedup")
    (output/"report.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("HISTORICAL QS3->ED2 SUITE",len(reports),"/",len(CASES),
          "levels",result["verified_levels"],"/",result["expected_levels"],
          "pass",passed,flush=True)
    return 0 if passed else 1


if __name__=="__main__":
    raise SystemExit(main())
