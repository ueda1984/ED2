#!/usr/bin/env python3
"""Two-lowest-eigenpair ALG=2 restart study on the stationary non-U(1) DM chain.

A single-eigenpair optimum must not be assumed to be efficient or reliable
for two low eigenpairs. We use precisely the same explicit ED2-basis
initial vector, fixed seed, absolute residual request and observable
checks for all restart workspaces. Independently constructed SciPy
Hamiltonians provide two lowest energy and ground observables, and
every accepted ED2 state must satisfy its direct residual policy.

The first excited state is also verified by energy and residual,
but only the *ground* observables are compared. This does not claim
near-degenerate projector/ensemble validation or an algorithm speedup.
CI-runner resource timings are exploratory rather than journal-grade.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

from mixed_spin_dm_chirality_scan import reference,run_ed2
from mixed_spin_matched_benchmark import initial_vector_text
from cpc_stationary_resource_benchmark import environment_record,stats,parse_workspaces


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="8,10")
    p.add_argument("--D-values",default="0,0.12")
    p.add_argument("--omp",type=int,choices=(1,2),default=1)
    p.add_argument("--workspaces",default="2:16,4:24,10:64")
    p.add_argument("--repeats",type=int,default=1)
    p.add_argument("--warmups",type=int,default=0)
    a=p.parse_args()
    try:
        sizes=[int(t.strip()) for t in a.sizes.split(",")]
        Ds=[float(t.strip()) for t in a.D_values.split(",")]
        wks=parse_workspaces(a.workspaces)
        if not sizes or len(sizes)!=len(set(sizes)) or any(n not in (4,6,8,10,12) for n in sizes):
            raise ValueError("Unique sizes must come from 4,6,8,10,12")
        if not Ds or len(Ds)!=len(set(Ds)) or any(abs(D)>0.5 for D in Ds):
            raise ValueError("Unique DM parameter values within +/-0.5 required")
        if not (1<=a.repeats<=10 and 0<=a.warmups<=5):
            raise ValueError("Invalid measurement count")
    except ValueError as e:
        p.error(str(e))
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    rows=[]
    refs=[]
    violations=[]
    for n in sizes:
        for D in Ds:
            r=reference(n,D,check_dense=(n<=6),profile="stationary")
            initial=initial_vector_text(r["dimension"])
            init_hash=hashlib.sha256(initial.encode()).hexdigest()
            refs.append(dict(sites=n,D=D,dim=r["dimension"],
                             lowest_energies=r["energies"],
                             first_gap=r["gap"],
                             independent_max_residual=max(r["reference_residual"]),
                             initial_sha256=init_hash))
            # Repetition outermost limits drift bias between workspaces.
            for rep in range(1,a.repeats+a.warmups+1):
                phase="warmup" if rep<=a.warmups else "measurement"
                for keep,capacity in wks:
                    name=f"n{n}-D{D:+.3f}-omp{a.omp}-NOE2-k{keep}-m{capacity}-rep{rep}"
                    folder=root/name
                    folder.mkdir(exist_ok=False)
                    try:
                        row=run_ed2(exe,folder,r,alg=2,threads=a.omp,
                                    keep=keep,capacity=capacity,
                                    initial_text=initial)
                        row.update(phase=phase,repeat=rep,
                                   initial_sha256=init_hash)
                        rows.append(row)
                        print("MULTIPAIR",phase,name,"status",row["status"],
                              "Hcalls",row["h_applications"],
                              "max_r",row.get("max_direct_residual"),
                              "seconds",row["wall_seconds"],
                              "rssKiB",row["peak_rss_kib"],flush=True)
                    except Exception as e:
                        violation=dict(case=name,error=str(e),class_name=type(e).__name__)
                        violations.append(violation)
                        print("MULTIPAIR CONTRACT VIOLATION",violation,flush=True)
    summary=[]
    for n in sizes:
        for D in Ds:
            for keep,capacity in wks:
                sample=[x for x in rows
                        if x["n_sites"]==n and x["D"]==D and
                           x["requested_keep"]==keep and x["requested_capacity"]==capacity and
                           x["phase"]=="measurement" and x["status"]==0]
                summary.append(dict(sites=n,D=D,omp=a.omp,
                                    NOE=2,NOK=keep,NOM=capacity,
                                    successful_measurements=len(sample),target=a.repeats,
                                    wall_seconds=stats(x["wall_seconds"] for x in sample),
                                    H_applications=stats(x["h_applications"] for x in sample),
                                    peak_RSS_KiB=stats(x["peak_rss_kib"] for x in sample
                                                       if x["peak_rss_kib"] is not None),
                                    achieved_max_residual=stats(x["max_direct_residual"] for x in sample)))
    expected=len(sizes)*len(Ds)*len(wks)*(a.repeats+a.warmups)
    failed_solver=sum(x["status"]!=0 for x in rows)
    success=(len(rows)==expected and failed_solver==0 and not violations and
             all(x["successful_measurements"]==a.repeats for x in summary))
    report=dict(schema=1,
                scope="same-two-eigenpair-target stationary non-U1 mixed-spin DM test",
                source_sha=os.environ.get("GITHUB_SHA","record from checkout"),
                note="ALG=2 restart workspaces only; no ALG=1 two-state timing comparison",
                timing_warning="GitHub shared runner values are characterizations, not isolated-node results",
                environment=environment_record(),
                settings=dict(sizes=sizes,D_values=Ds,OMP=a.omp,
                              workspaces=wks,NOE=2,absolute_residual=1e-8,
                              repeats=a.repeats,warmups=a.warmups,
                              explicit_initial=True,seed=19,wk_dim=-1,
                              OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1),
                references=refs,raw=rows,summary=summary,
                expected=expected,completed=len(rows),
                accepted=len(rows)-failed_solver,nonconverged=failed_solver,
                contract_violations=violations,pass_complete=success)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    (root/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("STATIONARY MULTIPAIR",len(rows),"/",expected,
          "nonconverged",failed_solver,"violations",len(violations),flush=True)
    return 0 if success else 1


if __name__=="__main__":
    raise SystemExit(main())
