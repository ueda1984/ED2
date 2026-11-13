#!/usr/bin/env python3
"""Construct an actual ED2 excited eigenvector and re-use it as ALG1 start.

This distinguishes residual certification from **global lowest-state
inclusion** without relying on a synthetic operator or invented vector.
First obtain the lowest two eigenvectors from ED2 ALG=2, verify both
energies against independent NumPy dense ED, then initialize ALG=1
from the *first excited eigenvector*, still satisfying the same H.
The start is in ED2's native basis (wavefunction binary), so no
assumed product-state ordering is involved.
"""
import argparse
import json
import os
from pathlib import Path

from mixed_spin_near_degenerate import reference,run_case,load_fortran_vector

DELTA=1e-6
N=4
ATOL=1e-12


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--omp",required=True,type=int,choices=(1,2))
    args=parser.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    ref=reference(N,DELTA)
    source=run_case(exe,root,N,DELTA,2,4,24,ATOL,args.omp,ref,suffix="-source")
    if source["status"]!=0 or len(source["returned_energies"])!=2:
        raise AssertionError("Preparation of accurately verified excited state failed")
    if any(abs(a-b)>5e-9 for a,b in
           zip(source["returned_energies"],ref["eigs"][:2])):
        raise AssertionError("Source solver's requested lowest pair is not independently confirmed")
    dirname=f"n{N}-dp1em06-alg2-k4-m24-a{ATOL:.0e}-omp{args.omp}-source"
    native_excited=load_fortran_vector(root/dirname,2,ref["dimension"])
    result=run_case(exe,root,N,DELTA,1,0,0,ATOL,args.omp,ref,
                    initial_vector=native_excited,suffix="-excited-initial")
    if result["status"]==0:
        e=float(result["returned_energies"][0])
        e0,e1=ref["eigs"][:2]
        is_excited=abs(e-e1)<5e-9 and (e-e0)>.5*DELTA
        observed_sz=float(result["isolated_spin_z"])
        is_excited_observable=abs(observed_sz-.5)<5e-4
        if is_excited and is_excited_observable:
            outcome="accepted_first_excited_not_ground"
        elif abs(e-e0)<5e-9 and abs(observed_sz+.5)<5e-4:
            outcome="converged_to_ground_despite_excited_start"
        else:
            outcome="accepted_mixed_or_other_low_doublet"
        result.update(explicit_excited_state_witness=is_excited and is_excited_observable,
                      excited_energy_error=abs(e-e1),
                      ground_energy_difference=e-e0,
                      outcome=outcome)
    else:
        result.update(explicit_excited_state_witness=False,
                      outcome="failed_or_refused_excited_start")
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record in CI"),
                scenario="isolated spin-half, delta=1e-6, actively non-U1 subsystem",
                target="one lowest eigenpair",requested_residual_atol=ATOL,
                analytic_gap=ref["gap"],reference_low_eigenvalues=ref["eigs"][:2],
                source_solver=source,result=result,omp=args.omp)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("EXCITED-START RESULT:",result["outcome"],
          "status=",result["status"],
          "eigenvalue=",result.get("returned_energies"),
          "independent_ground=",ref["eigs"][0],
          "independent_excited=",ref["eigs"][1],
          "spin0z=",result.get("isolated_spin_z"),
          "residual=",result.get("direct_residuals"),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
