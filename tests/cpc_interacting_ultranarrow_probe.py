#!/usr/bin/env python3
"""Classify, rather than hide, near-resolution-limit interacting doublets.

The preceding 20-case interacting sparse test used inner gaps >=1e-4.
This follow-up probes epsilon=5e-7 (anisotropic XY bond) with
longitudinal bias delta=-1e-7,0,+1e-7 and requests absolute residuals
1e-8 and 1e-12, with two restart configurations and OMP 1/2.

A residual-accepted *individual* Ritz vector need not be the unique
ground state when the gap is comparable to residual accuracy. An
accurate low-doublet ensemble does not imply that individual vectors
are uniquely resolved. Honest numerical nonconvergence is counted
separately from a broken runtime/physics validation contract.

The independent reference is a full sparse eight-site non-U(1)
interacting XY-edge Hamiltonian; the first bond of the otherwise
two-site stationary chain has deliberately been weakened.
"""
import argparse
import json
import os
from pathlib import Path
from cpc_stationary_interacting_sparse_doublet import (
    reference,ULTRA_SCENARIOS,RESTARTS,SPINS,
)
from mixed_spin_interacting_avoided_crossing import run_case

RESIDUALS=(1e-8,1e-12)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--omp",type=int,required=True,choices=(1,2))
    a=parser.parse_args()
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    raw=[]
    refs=[]
    violations=[]
    for eps,delta in ULTRA_SCENARIOS:
        model=reference(SPINS,eps,delta)
        refs.append(dict(n_sites=SPINS,epsilon=eps,delta=delta,
                         true_low_energies=model["low_e"],
                         internal_gap=model["ground_gap"],
                         external_gap=model["outer_gap"],
                         sparse_reference_residuals=model["reference_residuals"],
                         edge_component_commutators=model["edge_component_commutators"],
                         reference_edge_linear_entropy=model["edge_linear_entropy"],
                         nonu1_commutator_relative=model["commutator_relative"],
                         effective_low_gap=model["effective_low_gap"],
                         effective_low_gap_error=model["effective_low_gap_error"],
                         effective_firstorder_edge_z=model["effective_edge_z_firstorder"]))
        print("ULTRA REFERENCE",eps,delta,
              "gap",model["ground_gap"],"outer",model["outer_gap"],
              "edge entropy",model["edge_linear_entropy"],
              "effective_gap",model["effective_low_gap"],
              "effective_gap_error",model["effective_low_gap_error"],flush=True)
        for atol in RESIDUALS:
            for keep,capacity in RESTARTS:
                case=f"n{SPINS}-eps{eps:.1g}-delta{delta:+.1g}-k{keep}-m{capacity}-tol{atol:.0e}-omp{a.omp}"
                try:
                    r=run_case(exe,root,model,alg=2,keep=keep,
                               capacity=capacity,atol=atol,omp=a.omp,
                               suffix="-ultranarrow")
                    if r["status"] not in (0,1,2):
                        raise AssertionError("Unexpected status outside numerical protocol")
                    if r["status"]==0:
                        if not r.get("inside_low_cluster"):
                            raise AssertionError("Accepted pair lies outside low doublet")
                        if not r.get("direct_residuals") or max(r["direct_residuals"])>atol*(1+1e-12):
                            raise AssertionError("Residual-accepted pair violates request")
                        if r.get("cluster_ensemble_mag_error") is None:
                            raise AssertionError("Accepted doublet lacks native second-state trace")
                        if max(r["cluster_ensemble_mag_error"],
                               r["cluster_ensemble_cf_error"])>1e-5:
                            raise AssertionError("Accepted low-subspace projector trace disagrees with sparse reference")
                    else:
                        if r.get("classification")!="reported_nonconvergence_or_breakdown":
                            raise AssertionError("Nonconvergence was misclassified")
                    item=dict(case=case,omp=a.omp,epsilon=eps,delta=delta,
                              requested_atol=atol,NOK=keep,NOM=capacity,
                              physical_gap=model["ground_gap"],
                              external_gap=model["outer_gap"],
                              status=r["status"],classification=r["classification"],
                              direct_residuals=r.get("direct_residuals"),
                              individual_edge_z=r.get("measured_edge_sz"),
                              independent_ground_edge_z=model["edge_sz"],
                              individual_edge_z_error=r.get("edge_sz_error"),
                              trace_mag_max_error=r.get("cluster_ensemble_mag_error"),
                              trace_corr_max_error=r.get("cluster_ensemble_cf_error"),
                              h_applications=r["h_applications"])
                    raw.append(item)
                    print("ULTRA CASE",case,"status",item["status"],
                          "classification",item["classification"],
                          "individual_error",item["individual_edge_z_error"],
                          "trace_err",item["trace_corr_max_error"],
                          "Hcalls",item["h_applications"],flush=True)
                except Exception as exc:
                    err=dict(case=case,error=str(exc),error_type=type(exc).__name__)
                    violations.append(err)
                    print("ULTRA HARNESS VIOLATION",err,flush=True)
    expected=len(ULTRA_SCENARIOS)*len(RESIDUALS)*len(RESTARTS)
    accepted=sum(x["status"]==0 for x in raw)
    nonconv=sum(x["status"]!=0 for x in raw)
    individually_unresolved=sum(x["classification"]==
        "accepted_low_doublet_but_individual_ground_unresolved" for x in raw)
    # Do not demand all solver configurations converge; that would
    # conceal the very reliability behavior this test is designed
    # to measure. Every scenario must have *some* accepted solve.
    represented=all(any(x["epsilon"]==e and x["delta"]==d and
                        x["status"]==0 for x in raw)
                    for e,d in ULTRA_SCENARIOS)
    complete=(len(raw)==expected and not violations and represented)
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record source SHA"),
                purpose="ultra-small interacting XY-edge gap, residual versus subspace acceptance",
                warning="Near-resolution-limit characterization, not a universal convergence guarantee",
                omp=a.omp,source_references=refs,raw=raw,
                expected=expected,classified=len(raw),
                accepted=accepted,honest_nonconvergence=nonconv,
                individually_unresolved=individually_unresolved,
                scenarios_with_accepted_pair=represented,
                harness_violations=violations,pass_complete=complete)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("ULTRA SUMMARY",len(raw),"/",expected,
          "accepted",accepted,"nonconverged",nonconv,
          "individually unresolved",individually_unresolved,
          "harness violations",len(violations),flush=True)
    return 0 if complete else 1


if __name__=="__main__":
    raise SystemExit(main())
