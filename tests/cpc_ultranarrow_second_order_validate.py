#!/usr/bin/env python3
"""Validate a core-resolvent second-order effective Hamiltonian against ED2.

Each of three interacting, non-U(1), eight-site ultraweak edge cases is
independently assembled by SciPy and then treated by a separate 648-site
subproblem constrained Q-space resolvent. Both values are compared to the
existing ED2 two-lowest-eigenpair and observables validation runner.

Important:
  * An accurate gap does NOT establish a high-accuracy individual vector.
  * Difference from a finite-precision eigsh reference is an observed
    disagreement, not a certified error of ED2 alone.
  * The 2x2 second-order reduction omits O(epsilon^3) and delta-dependent
    virtual-excitation terms. Its accuracy is evaluated, not assumed.
  * This is a deliberately modified-edge stress model, not the CPC main
    connected stationary DM-response chain.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np

from cpc_stationary_interacting_sparse_doublet import (
    reference, ULTRA_SCENARIOS, RESTARTS, SPINS,
)
from cpc_ultranarrow_second_order_reference import effective_second_order
from mixed_spin_interacting_avoided_crossing import run_case


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--omp",type=int,choices=(1,2),required=True)
    p.add_argument("--requested-residual",type=float,default=1e-12)
    a=p.parse_args()
    if a.requested_residual not in (1e-12,1e-8):
        p.error("Only validated residual controls 1e-12 and 1e-8 permitted")
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    raw=[]
    refs=[]
    violations=[]
    for eps,delta in ULTRA_SCENARIOS:
        model=reference(SPINS,eps,delta)
        second=effective_second_order(model)
        first_gap=second["first_order_gap"]
        second_gap=second["second_order_gap"]
        sparse_gap=model["ground_gap"]
        gap_first_err=float(abs(first_gap-sparse_gap))
        gap_second_err=float(abs(second_gap-sparse_gap))
        if not (0<sparse_gap<2e-7 and 0<model["outer_gap"]<1):
            raise AssertionError("Reference lost spectral low-doublet isolation")
        if gap_second_err>2e-11:
            raise AssertionError("Second-order gap grossly inconsistent with full Hamiltonian")
        first_z=second["first_order_edge_z"]
        second_z=second["second_order_edge_z"]
        dressed_z=second["second_order_dressed_edge_z"]
        sparse_z=model["edge_sz"]
        refs.append(dict(n_sites=SPINS,epsilon=eps,delta=delta,
                         sparse_ground_edge_z=sparse_z,
                         sparse_low_gap=sparse_gap,
                         sparse_outer_gap=model["outer_gap"],
                         full_sparse_residuals=model["reference_residuals"],
                         first_order_gap_error=gap_first_err,
                         second_order_gap_error=gap_second_err,
                         first_order_edge_z_error_vs_sparse=abs(first_z-sparse_z),
                         second_order_edge_z_error_vs_sparse=abs(second_z-sparse_z),
                         dressed_edge_z_error_vs_sparse=abs(dressed_z-sparse_z),
                         second_order=second))
        print("SECOND-ORDER REFERENCE","eps",eps,"delta",delta,
              "gap_1st",first_gap,"gap_2nd",second_gap,
              "gap_sparse",sparse_gap,
              "gap_error_1st",gap_first_err,
              "gap_error_2nd",gap_second_err,
              "Sz_1st",first_z,"Sz_2nd",second_z,
              "Sz_dressed",dressed_z,"Sz_sparse",sparse_z,
              flush=True)
        for keep,capacity in RESTARTS:
            tag=f"n8-eps{eps:g}-d{delta:+.1e}-k{keep}-m{capacity}-omp{a.omp}"
            try:
                calc=run_case(exe,root,model,alg=2,keep=keep,
                              capacity=capacity,atol=a.requested_residual,
                              omp=a.omp,suffix="-2nd-order")
                if calc["status"]!=0 or not calc.get("inside_low_cluster"):
                    raise AssertionError("ED2 did not accept the independently checked low doublet")
                residuals=calc.get("direct_residuals",[])
                if len(residuals)!=2 or max(residuals)>a.requested_residual*(1+2e-12):
                    raise AssertionError("ED2 direct residual contract violated")
                if calc.get("cluster_ensemble_mag_error") is None or max(
                        calc["cluster_ensemble_mag_error"],
                        calc["cluster_ensemble_cf_error"])>1e-5:
                    raise AssertionError("ED2 native wavefunction trace not independently accepted")
                value=calc["measured_edge_sz"]
                row=dict(tag=tag,OMP=a.omp,epsilon=eps,delta=delta,
                         NOK=keep,NOM=capacity,
                         status=calc["status"],
                         classification=calc["classification"],
                         max_direct_residual=max(residuals),
                         sparse_ground_edge_z=sparse_z,
                         ED2_ground_edge_z=value,
                         ED2_vs_sparse_abs_error=abs(value-sparse_z),
                         ED2_vs_second_order_abs_diff=abs(value-second_z),
                         ED2_vs_dressed_abs_diff=abs(value-dressed_z),
                         ED2_vs_first_order_abs_diff=abs(value-first_z),
                         second_order_vs_sparse_abs_diff=abs(second_z-sparse_z),
                         dressed_vs_sparse_abs_diff=abs(dressed_z-sparse_z),
                         first_order_gap_err=gap_first_err,
                         second_order_gap_err=gap_second_err,
                         outside_gap=model["outer_gap"],
                         H_applications=calc["h_applications"],
                         two_state_trace_mag_error=calc["cluster_ensemble_mag_error"],
                         two_state_trace_corr_error=calc["cluster_ensemble_cf_error"])
                raw.append(row)
                print("SECOND-ORDER ED2 PASS",tag,
                      "ED2_Sz",value,
                      "err_sparse",row["ED2_vs_sparse_abs_error"],
                      "err_2nd",row["ED2_vs_second_order_abs_diff"],
                      "trace",row["two_state_trace_corr_error"],
                      flush=True)
            except Exception as err:
                failure=dict(tag=tag,type=type(err).__name__,message=str(err))
                violations.append(failure)
                print("SECOND-ORDER ED2 FAILURE",failure,flush=True)
    expected=len(ULTRA_SCENARIOS)*len(RESTARTS)
    passed=len(raw)==expected and not violations
    report=dict(schema=1,
                source_sha=os.environ.get("GITHUB_SHA","record exact source SHA"),
                model="8-site alternating spin-1/2,1 XYZ/Gamma/staggered-DM core with ultraweak XY edge",
                caveat="Second-order perturbation and sparse eigsh are separate approximate references, not interval-certified exact wavefunctions",
                OMP=a.omp,requested_direct_residual=a.requested_residual,
                references=refs,raw=raw,
                expected=expected,passed=len(raw),violations=violations,
                pass_complete=passed)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("SECOND-ORDER VALIDATION",len(raw),"/",expected,
          "violations",len(violations),"OMP",a.omp,"pass",passed,flush=True)
    return 0 if passed else 1


if __name__=="__main__":
    raise SystemExit(main())
