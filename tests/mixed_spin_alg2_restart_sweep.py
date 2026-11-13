#!/usr/bin/env python3
"""Exploratory ALG=2 restart workspace sweep at fixed one-pair target.

One explicit starting vector, one requested eigenpair, one residual
threshold, unchanged Hamiltonian, and unchanged ED2 algorithm. Records
normal convergence AND budget exhaustion, without silently changing
tolerances or the legacy stopping gate.
"""
import argparse
import json
import os
from pathlib import Path

from mixed_spin_matched_benchmark import (
    run_one, initial_vector_text, FIELD,
)
from mixed_spin_sparse_scale import reference_model

N=10
RESIDUAL_ATOL=1e-8
WORKSPACES=((2,16),(4,24),(6,32),(8,48),(10,64),(12,96))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    model=reference_model(N,FIELD)
    initial=initial_vector_text(model["dimension"])
    tested=[]
    violations=[]
    baseline=root/"baseline"
    baseline.mkdir()
    try:
        value=run_one(exe,baseline,N,1,RESIDUAL_ATOL,1,model,initial)
        if value["status"]!=0:
            raise AssertionError("Baseline ALG1 did not converge")
        tested.append(value)
        print(f"BASELINE ALG1: Hcalls={value['h_applications']} "
              f"wall={value['wall_seconds']:.3f} RSS={value['peak_rss_kib']}",flush=True)
    except Exception as e:
        violations.append(dict(label="alg1-baseline",error=str(e)))
    for keep,capacity in WORKSPACES:
        block=root/f"keep{keep}-capacity{capacity}"
        block.mkdir()
        try:
            value=run_one(exe,block,N,2,RESIDUAL_ATOL,1,model,initial,
                          keep=keep,capacity=capacity)
            tested.append(value)
            print(f"WORKSPACE keep={keep} cap={capacity}: "
                  f"status={value['status']} reason={value['reason']} "
                  f"Hcalls={value['h_applications']} "
                  f"wall={value['wall_seconds']:.3f} "
                  f"RSS={value['peak_rss_kib']} "
                  f"residual={value.get('residual')}",flush=True)
        except Exception as e:
            violations.append(dict(label=f"keep{keep}-cap{capacity}",error=str(e)))
            print("CONTRACT VIOLATION:",violations[-1],flush=True)
    report=dict(schema=1,kind="ALG2 restart cost at matched one-eigenpair goal",
                git_sha=os.environ.get("GITHUB_SHA","record in CI"),
                n_sites=N,dimension=model["dimension"],field_control=FIELD,
                residual_atol=RESIDUAL_ATOL,explicit_start_shared=True,
                baseline_and_alg2_configurations=tested,
                contract_violations=violations,
                total_recorded=len(tested),
                converged=sum(x["status"]==0 for x in tested),
                nonconverged=sum(x["status"]!=0 for x in tested))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"WORKSPACE SWEEP: {report['total_recorded']} records, "
          f"{report['converged']} accepted, {report['nonconverged']} nonconverged, "
          f"{len(violations)} contract violations",flush=True)
    return 0 if report["total_recorded"]==1+len(WORKSPACES) and not violations else 1


if __name__=="__main__":
    raise SystemExit(main())
