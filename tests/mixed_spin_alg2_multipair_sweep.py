#!/usr/bin/env python3
"""Test ALG=2 two-eigenpair workspace robustness against independent sparse ED.

All input decks, diagnostics and numerical outputs are preserved, including
unsuccessful solver trials. Accepted results MUST independently match the
lowest two physical eigenvalues and ground-state observables. A correctly
signalled nonconvergence is reported, not reclassified as a successful run.

The default 10/64 workspace is the validated baseline, not a presumed
optimal choice. Configurations are selected before seeing their outcomes.
"""
import argparse
import json
import os
from pathlib import Path

from mixed_spin_sparse_scale import reference_model,run_ed2

WORKSPACES=((2,16),(4,24),(10,64))
MODELS=((8,0.15),(8,0.38),(10,0.15))


def read_fields(path):
    fields={}
    if not path.is_file():
        return fields
    for s in path.read_text().splitlines():
        if "=" in s:
            k,v=s.split("=",1)
            fields[k.strip()]=v.strip()
    return fields


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    results=[]
    failures=[]
    for n,hx in MODELS:
        reference=reference_model(n,hx)
        gap=float(reference["energies"][1]-reference["energies"][0])
        print(f"REFERENCE n{n}, hx={hx}: lowest pair {reference['energies']} "
              f"gap={gap:.7e} independent r={max(reference['eigenpair_residuals']):.3e}",flush=True)
        for keep,capacity in WORKSPACES:
            label=f"mixed-n{n}-hx{hx:.2f}-alg2-omp1"+(
                f"-k{keep}-m{capacity}" if (keep,capacity)!=(10,64) else "")
            try:
                entry=run_ed2(exe,root,n,hx,2,1,reference,keep=keep,capacity=capacity)
                entry["independent_gap"]=gap
                results.append(entry)
                print(f"ACCEPTED {label}: Hcalls={entry['h_applications']} "
                      f"energyErr={entry['energy_max_error']:.2e} "
                      f"r={entry['max_reported_residual']:.2e}",flush=True)
            except Exception as exc:
                folder=root/label
                diag=read_fields(folder/"diagnostics.txt")
                status=diag.get("status","missing")
                reason=diag.get("reason","missing diagnostic status")
                outputs=[n for n in ("eigenvalues.dat","local_mag.dat","two_body_cf_xyz.dat")
                         if (folder/n).exists()]
                if status in ("1","2") and not outputs:
                    entry=dict(case=label,n_sites=n,hx=hx,
                               independent_gap=gap,alg=2,
                               requested_keep=keep,requested_capacity=capacity,
                               status="nonconverged",solver_status=int(status),
                               reason=reason,h_applications=int(diag.get("h_applications","-1")))
                    results.append(entry)
                    print(f"NONCONVERGED {label}: status={status} reason={reason} "
                          f"Hcalls={entry['h_applications']}",flush=True)
                else:
                    entry=dict(case=label,n_sites=n,hx=hx,
                               independent_gap=gap,requested_keep=keep,
                               requested_capacity=capacity,status="test_failure",
                               diagnostics_status=status,diagnostics_reason=reason,
                               outputs_present=outputs,error=str(exc))
                    results.append(entry)
                    failures.append(entry)
                    print(f"TEST FAILURE {label}: {str(exc)[:300]}",flush=True)
    accepted=sum(e["status"]=="passed" for e in results)
    nonconverged=sum(e["status"]=="nonconverged" for e in results)
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record in CI"),
                model="open mixed-spin non-U1 XYZ/DM/Gamma with transverse fields",
                wanted_pairs=2,alg=2,residual_atol=1e-8,omp=1,seed=19,
                workspaces=list(WORKSPACES),models=[list(x) for x in MODELS],
                results=results,total=len(results),accepted=accepted,
                nonconverged=nonconverged,failures=len(failures))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"MULTI-PAIR SWEEP: {len(results)} trials; {accepted} accepted; "
          f"{nonconverged} honest nonconvergences; {len(failures)} test failures",flush=True)
    by_model={(x["n_sites"],x["hx"]):x for x in results
              if x.get("requested_keep")==10 and x.get("requested_capacity")==64}
    if len(results)!=len(WORKSPACES)*len(MODELS) or failures or len(by_model)!=len(MODELS):
        return 1
    if any(x["status"]!="passed" for x in by_model.values()):
        # The prevalidated baseline must work on these tested models.
        return 1
    return 0


if __name__=="__main__":
    raise SystemExit(main())
