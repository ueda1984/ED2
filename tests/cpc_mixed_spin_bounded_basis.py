#!/usr/bin/env python3
"""Verify ED2 generalized bounded-occupation basis against independent enumeration.

This *user-facing* numerical test uses the Python standard library only.
The original published QS3 fixed-Ndown basis has size binomial(L,Ndown).
ED2 represents mixed sites with local maximum m_i and an allowed total
occupation window, counted by coefficients of prod_i (1+...+x**m_i).

A diagonal mixed-spin Hamiltonian ensures the complete eigenvalue spectrum
can be enumerated exactly in Python, independent of ED2 basis order.
Three four-site models share identical physical couplings but use
  (NODmin,NODmax)=(0,6),(3,3),(2,4)
to exercise full and restricted representations, dimensions 36,10,26.

This checks the representation and file contract, NOT the accuracy of
projecting arbitrary non-U(1) Hamiltonians to an occupation window.
"""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import subprocess

S=(.5,1.,.5,1.)
CAP=(1,2,1,2)
FIELDS=(.17,-.21,.11,-.09)
JZZ=(.41,-.27,.31)
WINDOWS=((0,6),(3,3),(2,4))
ABS_TOL=5e-12


def integer_counts(caps):
    c=[1]
    for m in caps:
        out=[0]*(len(c)+m)
        for n,v in enumerate(c):
            for z in range(m+1):out[n+z]+=v
        c=out
    return c


def count_spin_half_limit():
    from math import comb
    for sites in (4,6,8):
        count=integer_counts([1]*sites)
        if any(count[n]!=comb(sites,n) for n in range(sites+1)):
            raise AssertionError("Old QS3 binomial limit not recovered")
    if integer_counts(CAP)!=[1,4,8,10,8,4,1]:
        raise AssertionError("Independent bounded-occupation recurrence mismatch")
    return True


def independent_diagonal(low,high):
    spectra=[]
    for state in itertools.product(*(range(m+1) for m in CAP)):
        n=sum(state)
        if low<=n<=high:
            sz=[spin-z for spin,z in zip(S,state)]
            e=sum(h*z for h,z in zip(FIELDS,sz))
            e+=sum(j*sz[k]*sz[k+1] for k,j in enumerate(JZZ))
            spectra.append(e)
    return sorted(spectra)


def prepare(root,low,high):
    root.mkdir(parents=True,exist_ok=False)
    (root/"spins.dat").write_text("".join(f"{x:.1f}\n" for x in S))
    (root/"caps.dat").write_text("".join(f"{x}\n" for x in CAP))
    (root/"fields.dat").write_text(
        "".join(f"{i} 0 0 {h:.17g}\n" for i,h in enumerate(FIELDS,1)))
    (root/"bonds.dat").write_text("".join(
        f"{i} {i+1} 0 0 {JZZ[i-1]:.17g} 0 0 0 0 0 0\n"
        for i in range(1,4)))
    deck=f"""&input_parameters
NOS=4, NODmin={low}, NODmax={high}, ALG=3,
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,
M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one=4,NO_two=3,wk_dim=0,MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,
cal_lm=0,cal_cf=0
/
&input_lancz
lnc_ene_conv=1d-12,minitr=4,maxitr=1000,itrint=5
/
&input_TRLan
NOE=2,NOK=2,NOM=4,maxitr=1000,lnc_ene_conv=1d-12,
i_vec_min=1,i_vec_max=1
/
"""
    (root/"input.dat").write_text(deck)


def sha(path):
    digest=hashlib.sha256()
    with path.open("rb") as file:
        for part in iter(lambda:file.read(1<<20),b""):
            digest.update(part)
    return digest.hexdigest()


def run_one(exe,root,low,high,timeout):
    work=root/f"window-{low}-{high}"
    prepare(work,low,high)
    env=dict(os.environ,OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as input_file:
        p=subprocess.run([str(exe)],stdin=input_file,cwd=work,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                         env=env,timeout=timeout,check=False)
    (work/"stdout.log").write_bytes(p.stdout)
    if p.returncode:
        raise AssertionError("ED2 failed on bounded-occupation test: "+p.stdout.decode(errors="replace")[-1700:])
    text=p.stdout.decode(errors="replace")
    dims=re.findall(r"THS\s*=\s*(\d+)",text)
    ref=independent_diagonal(low,high)
    count=integer_counts(CAP)
    expected=sum(count[low:high+1])
    if len(ref)!=expected or not dims:
        raise AssertionError("Independent basis dimension or ED2 dimension output missing")
    if int(dims[0])!=expected:
        raise AssertionError("Fortran basis count contradicts polynomial coefficients")
    file=work/"eigenvalues.dat"
    if not file.exists():
        raise AssertionError("Dense numerical spectrum missing")
    parsed=[]
    for line in file.read_text().splitlines():
        if not line.strip():continue
        vals=line.split()
        if len(vals)!=2 or int(vals[0])!=len(parsed)+1:
            raise AssertionError("ED2 spectrum row malformed")
        energy=float(vals[1].replace("D","E"))
        if not math.isfinite(energy):raise AssertionError("Non-finite eigenenergy")
        parsed.append(energy)
    if len(parsed)!=expected:
        raise AssertionError("Numerical output dimension doesn't match exact combinatorics")
    error=max(abs(x-y) for x,y in zip(parsed,ref))
    if error>ABS_TOL:
        raise AssertionError("Hamiltonian diagonalization contradicts independent occupation enumeration")
    return dict(min_down=low,max_down=high,
                exact_occupation_coefficients=count,
                expected_dimension=expected,reported_fortran_dimension=int(dims[0]),
                number_of_eigenvalues=len(parsed),
                eigenvalue_error_max=error,reference_spectrum=ref,
                returned_spectrum=parsed,case_path=str(work.relative_to(root)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--timeout",type=int,default=120)
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    if not 5<=args.timeout<=600:p.error("Invalid numerical timeout")
    assert count_spin_half_limit()
    passed=[];issues=[]
    for low,high in WINDOWS:
        try:
            case=run_one(exe,out,low,high,args.timeout)
            passed.append(case)
            print("BOUNDED BASIS PASS",low,high,case["expected_dimension"],
                  case["eigenvalue_error_max"],flush=True)
        except Exception as error:
            issues.append(dict(window=[low,high],error_type=type(error).__name__,
                               message=str(error)))
            print("BOUNDED BASIS ERROR",issues[-1],flush=True)
    good=len(passed)==len(WINDOWS) and not issues
    report=dict(schema=1,pass_complete=good,
                tested_source=os.environ.get("GITHUB_SHA","record checkout SHA"),
                method="independent bounded-occupation polynomial and diagonal Hamiltonian",
                spin_magnitudes=S,local_caps=CAP,fields_z=FIELDS,exchange_z=JZZ,
                exact_qs3_binomial_limit=True,
                expected_cases=len(WINDOWS),completed_cases=len(passed),
                results=passed,failures=issues,executable_sha256=sha(exe),
                caveat="Non-U1 occupation windows are projections, not invariant symmetry sectors")
    (out/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("BOUNDED BASIS",len(passed),"/",len(WINDOWS),
          "expected dimensions 36/10/26, accepted",good,flush=True)
    return 0 if good else 1


if __name__=="__main__":
    raise SystemExit(main())
