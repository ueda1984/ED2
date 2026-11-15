#!/usr/bin/env python3
"""Common-sector iterative Lanczos continuity: published QS3/ED vs ED2.

Build historical QS-Cube/ED at exact source commit
a3b09934bba5f966bfb564718528de8778c2b46c without altering source.
Build ED2 from a clean extracted source package. Execute ALG=1
(conventional Lanczos) in BOTH original executables, using exactly the
same physical open-boundary spin-half XXZ Hamiltonian and Ndown=2
sector at L=12,18,24. Compare each computed lowest energy against an
independently constructed product-basis Hermitian matrix.

The historical program uses an uncontrolled random starting vector.
Do not expect identical trajectories, iteration counts or bitwise
energies. This is a scientific eigenenergy regression, NOT a matched
initial-state or controlled runtime comparison.

NumPy is a TEST-ONLY prerequisite. Fortran ED2 itself does not need it.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time
import numpy as np

# Reuse the original reviewed physical input translator.
from cpc_qs3_ed_shared_xxz_extended import prepare, independent
from cpc_qs3_ed_shared_xxz_crosscheck import ORIGINAL_ED_SHA, JXY, JZ

LENGTHS=(12,18,24)
N_DOWN=2
REQUESTED_CONVERGENCE=1e-11
MAX_ITERATIONS=2400
ACCEPT_ENERGY=2e-8

ENERGY_LINE=re.compile(r"E_0\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)[EeDd][-+]?\d+)")
STEP_LINE=re.compile(r"total\s+lanczos\s+step:\s*(\d+)",re.IGNORECASE)


def sha256(file):
    h=hashlib.sha256()
    with file.open("rb") as s:
        for chunk in iter(lambda:s.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def make_iterative_input(folder,n,old):
    """Translate physical XXZ input into each historical NAMELIST ABI."""
    prepare(folder,n,N_DOWN,old)
    path=folder/"input.dat"
    deck=path.read_text()
    if deck.count("ALG=3")!=1:
        raise AssertionError("Unexpected original input template")
    deck=deck.replace("ALG=3","ALG=1")
    if old:
        deck+=f"""&input_lancz
lnc_ene_conv={REQUESTED_CONVERGENCE:.2e},
minitr=5, maxitr={MAX_ITERATIONS}, itrint=5
&end
"""
    path.write_text(deck)
    return path


def run(exe,folder,timeout):
    env=dict(os.environ,OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             GOTO_NUM_THREADS="1")
    start=time.monotonic()
    with (folder/"input.dat").open("rb") as deck:
        completed=subprocess.run([str(exe)],cwd=folder,stdin=deck,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
            check=False,timeout=timeout,env=env)
    duration=time.monotonic()-start
    (folder/"stdout.log").write_bytes(completed.stdout)
    output=completed.stdout.decode(errors="replace")
    if completed.returncode:
        raise AssertionError(f"ALG1 exited {completed.returncode}: {output[-2000:]}")
    energies=ENERGY_LINE.findall(output)
    steps=STEP_LINE.findall(output)
    if len(energies)!=1 or len(steps)!=1:
        raise AssertionError("Missing/ambiguous solver energy or iterations: "+output[-1900:])
    energy=float(energies[0].replace("D","E"))
    iteration=int(steps[0])
    if not math.isfinite(energy) or not 0<iteration<=MAX_ITERATIONS:
        raise AssertionError("Invalid energy or iteration count")
    return dict(energy=energy,iterations=iteration,
                elapsed_seconds_uncontrolled=duration,
                input_sha256=sha256(folder/"input.dat"))


def validate_ed2_output(folder,reported_energy):
    path=folder/"eigenvalues.dat"
    if not path.is_file():
        raise AssertionError("ED2 eigenvalue output not present")
    lines=[line.split() for line in path.read_text().splitlines() if line.strip()]
    if len(lines)!=1 or len(lines[0])!=2 or int(lines[0][0])!=1:
        raise AssertionError("ED2 energy file does not contain exactly one indexed level")
    value=float(lines[0][1].replace("D","E"))
    if not math.isfinite(value) or abs(value-reported_energy)>5e-11:
        raise AssertionError("ED2 indexed energy differs from STDOUT")
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-exe",type=Path,required=True)
    parser.add_argument("--new-exe",type=Path,required=True)
    parser.add_argument("--old-commit",required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--timeout",type=int,default=120)
    a=parser.parse_args()
    if a.old_commit!=ORIGINAL_ED_SHA:
        parser.error("Published QS3-ED source SHA must match reviewed original")
    if not 10<=a.timeout<=600:
        parser.error("Invalid timeout")
    old_binary=a.old_exe.resolve(strict=True)
    new_binary=a.new_exe.resolve(strict=True)
    out=a.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    records=[];failures=[]
    for n in LENGTHS:
        root=out/f"size{n}-ndown{N_DOWN}"
        root.mkdir(parents=True,exist_ok=False)
        old_dir=root/"old-qs3-ed"
        new_dir=root/"modern-ed2"
        try:
            make_iterative_input(old_dir,n,True)
            make_iterative_input(new_dir,n,False)
            # Exact independent low spectrum; distinct from either solver.
            independent_eigvals=independent(n,N_DOWN)
            exact=independent_eigvals[0]
            next_low=independent_eigvals[1]
            old=run(old_binary,old_dir,a.timeout)
            modern=run(new_binary,new_dir,a.timeout)
            modern_file=validate_ed2_output(new_dir,modern["energy"])
            errors={
                "old_abs_error":abs(old["energy"]-exact),
                "ed2_abs_error":abs(modern["energy"]-exact),
                "old_ed2_abs_difference":abs(old["energy"]-modern["energy"]),
                "ed2_file_stdout_difference":abs(modern_file-modern["energy"])
            }
            if max(errors.values())>ACCEPT_ENERGY:
                raise AssertionError("Energy difference above independent tolerance: "+str(errors))
            r=dict(sites=n,n_down=N_DOWN,dimension=math.comb(n,N_DOWN),
                couplings=dict(Jxy=JXY,Jz=JZ),algorithm="ALG=1",
                requested_energy_convergence=REQUESTED_CONVERGENCE,
                independent_lowest_energy=exact,independent_second_energy=next_low,
                independent_gap=next_low-exact,
                historical_qs3=old,ed2=modern,errors=errors,pass_complete=True)
            records.append(r)
            print("ITERATIVE ORIGINAL QS3/ED -> ED2:",n,
                  "dim",r["dimension"],"old",old["energy"],
                  "new",modern["energy"],"exact",exact,
                  "iterations",old["iterations"],modern["iterations"],
                  "errors",errors,"PASS",flush=True)
        except Exception as ex:
            error=dict(sites=n,error_type=type(ex).__name__,message=str(ex))
            failures.append(error)
            print("ITERATIVE ORIGINAL QS3/ED -> ED2:",error,"FAIL",flush=True)
    passed=len(records)==len(LENGTHS) and not failures
    summary=dict(schema=1,pass_complete=passed,
         published_program="QS3, CPC 277 (2022) 108369, QS-Cube/ED",
         original_source_commit=a.old_commit,
         tested_ed2_commit=os.environ.get("GITHUB_SHA","not set; record source manually"),
         original_executable_sha256=sha256(old_binary),
         ed2_executable_sha256=sha256(new_binary),
         method="same finite-size U(1) physical XXZ Hamiltonian; both ALG1",
         expected_cases=len(LENGTHS),passed_cases=len(records),
         records=records,failures=failures,
         caution="Historical randomized initialization; not matched-start, controlled-node or a timing benchmark")
    (out/"report.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("ITERATIVE QS3 -> ED2 MATCHED XXZ PHYSICS:",len(records),
          "/",len(LENGTHS),"pass",passed,flush=True)
    return 0 if passed else 1


if __name__=="__main__":
    raise SystemExit(main())
