#!/usr/bin/env python3
"""Cross-version numerical equivalence on an exactly shared QS3-ED/ED2 XXZ sector.

Historical reference code: QS-Cube/ED (CPC 277,108369,2022), exact commit
a3b09934bba5f966bfb564718528de8778c2b46c.
Modern code: QS-Cube/ED2 lineage, selected development SHA.

The same FOUR-site open spin-1/2 Hamiltonian is represented in:
  historical QS3_only_u1: NOD=2, Jxy and Jz, 6 basis vectors.
  ED2: NODmin=NODmax=2, Jx=Jy=Jxy, Jz, D=Gamma=h=0, 6 vectors.

Both implementations use their legacy *dense* LAPACK paths (ALG=3),
thus avoiding a confound from different Krylov stopping criteria.
We independently assemble the six-dimensional sector in NumPy.
No claim is made that this tiny test establishes matched performance
or coverage of the original QS3 dynamical spin structure factor.

The old program and its output files are not edited; each run is fresh.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ORIGINAL_ED_SHA="a3b09934bba5f966bfb564718528de8778c2b46c"
JXY=0.8
JZ=1.1
SITES=4
NOD=2
BONDS=((1,2),(2,3),(3,4))
ATOL=3e-9


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for a in iter(lambda:f.read(1024*1024),b""):
            h.update(a)
    return h.hexdigest()


def independent_reference():
    states=[n for n in range(1<<SITES) if n.bit_count()==NOD]
    indices={b:i for i,b in enumerate(states)}
    assert len(states)==6
    H=np.zeros((6,6),dtype=float)
    for k,mask in enumerate(states):
        for i,j in BONDS:
            a=(mask>>(i-1))&1
            b=(mask>>(j-1))&1
            H[k,k]+=JZ*(.5-a)*(.5-b)
            if a!=b:
                flipped=mask^(1<<(i-1))^(1<<(j-1))
                H[k,indices[flipped]]+=JXY*.5
    if np.max(np.abs(H-H.T))>1e-14:
        raise AssertionError("Independent six-dimensional sector non-Hermitian")
    return np.linalg.eigvalsh(H).tolist()


def prepare_old(directory):
    directory.mkdir(parents=True,exist_ok=False)
    (directory/"bonds_old.dat").write_text("".join(
        f"{i} {j} {JXY:.17g} {JZ:.17g}\n" for i,j in BONDS))
    (directory/"input.dat").write_text("""&input_parameters
NOS=4, NOD=2, NOxxz=3, ALG=3,
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
""")
    return directory/"input.dat"


def prepare_new(directory):
    directory.mkdir(parents=True,exist_ok=False)
    (directory/"bonds_new.dat").write_text("".join(
        f"{i} {j} {JXY:.17g} {JXY:.17g} {JZ:.17g} 0 0 0 0 0 0\n" for i,j in BONDS))
    (directory/"fields.dat").write_text("".join(
        f"{i} 0 0 0\n" for i in range(1,5)))
    (directory/"spins.dat").write_text("0.5\n"*4)
    (directory/"caps.dat").write_text("1\n"*4)
    (directory/"input.dat").write_text("""&input_parameters
NOS=4, NODmin=2, NODmax=2, ALG=3,
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one=4, NO_two=3, wk_dim=0, MNTE=-1,
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
""")
    return directory/"input.dat"


def run(executable,dir_name,timeout_seconds):
    env=dict(os.environ,OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    t=time.perf_counter()
    with (dir_name/"input.dat").open("rb") as source:
        proc=subprocess.run([str(executable)],cwd=dir_name,stdin=source,
                            env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=timeout_seconds,check=False)
    (dir_name/"stdout.log").write_bytes(proc.stdout)
    if proc.returncode!=0:
        raise AssertionError("Dense solver exited nonzero: "+str(proc.returncode)+
                             " "+proc.stdout.decode(errors="replace")[-2000:])
    return time.perf_counter()-t


def read_energies(path,number,columns):
    if not path.is_file():raise AssertionError("Missing eigenvalue file "+str(path))
    raw=[]
    for line in path.read_text().splitlines():
        if not line.strip():continue
        values=line.split()
        if len(values)!=columns:raise AssertionError("Wrong eigenvalue file format")
        raw.append(float(values[-1].replace("D","E")))
    if len(raw)!=number or not np.all(np.isfinite(raw)):
        raise AssertionError("Wrong number/nonfinite eigenvalues")
    return raw


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old-exe",type=Path,required=True)
    p.add_argument("--new-exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--old-commit",required=True)
    p.add_argument("--timeout",type=int,default=90)
    a=p.parse_args()
    if a.old_commit!=ORIGINAL_ED_SHA:
        p.error("Old QS3-ED baseline is not pinned to the reviewed CPC-era commit")
    if not 5<=a.timeout<=600:p.error("Invalid timeout")
    old=a.old_exe.resolve(strict=True)
    new=a.new_exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    prepare_old(root/"old-qs3-ed")
    prepare_new(root/"modern-ed2")
    elapsed_old=run(old,root/"old-qs3-ed",a.timeout)
    elapsed_new=run(new,root/"modern-ed2",a.timeout)
    reference=independent_reference()
    legacy=read_energies(root/"old-qs3-ed"/"full_eigenvalues.dat",6,1)
    modern=read_energies(root/"modern-ed2"/"eigenvalues.dat",6,2)
    errors=dict(old_ref=max(abs(a-b) for a,b in zip(legacy,reference)),
                new_ref=max(abs(a-b) for a,b in zip(modern,reference)),
                old_new=max(abs(a-b) for a,b in zip(legacy,modern)))
    good=all(v<=ATOL for v in errors.values())
    result=dict(schema=1,pass_complete=good,
                previous_program="QS-Cube/ED, CPC 277 (2022) 108369",
                old_commit=a.old_commit,
                new_source_sha=os.environ.get("GITHUB_SHA","record tested source"),
                same_physical_model="spin 1/2, open, 4 sites, XXZ Jxy=0.8 Jz=1.1",
                sector="Ndown=2, dimension 6",
                previous_method="QS3_only_u1 dense LAPACK ALG=3",
                new_method="QS3ED2 dense LAPACK ALG=3, NODmin=NODmax=2",
                independent_product_basis_reference=reference,
                old_spectrum=legacy,new_spectrum=modern,
                maximum_absolute_disagreements=errors,
                old_binary_sha256=digest(old),new_binary_sha256=digest(new),
                observed_elapsed_seconds_not_comparable=dict(old=elapsed_old,new=elapsed_new),
                warning="This checks a shared physics limit. It does not show runtime superiority, reproduce legacy DSF, or guarantee all possible old inputs.")
    (root/"report.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("QS3-ED vs ED2 COMMON XXZ SPECTRUM", "6/6",errors,
          "passed",good,flush=True)
    return 0 if good else 1


if __name__=="__main__":
    raise SystemExit(main())
