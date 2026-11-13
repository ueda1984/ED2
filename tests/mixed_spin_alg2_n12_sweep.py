#!/usr/bin/env python3
"""Validate tuned ALG2 workspaces at 12 sites against sparse reference.

Validates *one* requested lowest eigenpair, a shared explicit initial
vector and absolute residual target 1e-8. Retains nonconvergent cases
without converting them into successes.
"""
import argparse
import json
import os
from pathlib import Path

from mixed_spin_matched_benchmark import run_one, initial_vector_text, FIELD
from mixed_spin_sparse_scale import reference_model

N=12
RESIDUAL_ATOL=1e-8
WORKSPACES=((2,16),(4,24),(6,32),(10,64))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    model=reference_model(N,FIELD)
    initial=initial_vector_text(model["dimension"])
    records=[]
    violations=[]
    targets=[("ALG1",1,None,None)]+[
        (f"ALG2-k{keep}-m{cap}",2,keep,cap)
        for keep,cap in WORKSPACES
    ]
    for label,alg,keep,cap in targets:
        folder=root/label
        folder.mkdir()
        try:
            if alg==1:
                result=run_one(exe,folder,N,alg,RESIDUAL_ATOL,1,model,initial)
            else:
                result=run_one(exe,folder,N,alg,RESIDUAL_ATOL,1,model,initial,
                               keep=keep,capacity=cap)
            records.append(result)
            print(f"12-SITE WORKSPACE {label}: status={result['status']} "
                  f"Hcalls={result['h_applications']} "
                  f"wall={result['wall_seconds']:.3f} "
                  f"RSS={result['peak_rss_kib']} "
                  f"residual={result.get('residual')} "
                  f"reason={result['reason']}",flush=True)
        except Exception as e:
            violations.append(dict(label=label,error=str(e)))
            print(f"CONTRACT VIOLATION {label}: {e}",flush=True)
    report=dict(schema=1,kind="12-site mixed-spin ALG2 restart sweep",
                source_sha=os.environ.get("GITHUB_SHA","record in CI"),
                n_sites=N,dimension=model["dimension"],
                records=records,contract_violations=violations,
                total_recorded=len(records),
                accepted=sum(x["status"]==0 for x in records),
                nonconverged=sum(x["status"]!=0 for x in records))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"12-SITE SWEEP: {report['total_recorded']} records, "
          f"{report['accepted']} accepted, {report['nonconverged']} nonconverged, "
          f"{len(violations)} contract violations",flush=True)
    return 0 if len(records)==len(targets) and not violations else 1


if __name__=="__main__":
    raise SystemExit(main())
