#!/usr/bin/env python3
"""Measure ED2's sparse-row capacity MNTE without changing the Hamiltonian.

At fixed *full* wk_dim caching, compare the automatically sized work arrays
(MNTE=-1) against the smallest capacity actually observed in a full
unrestricted scan (the existing ED2 "Optimal MNTE" counter).
The calibration is made independently for every size and is not assumed
valid for a different Hamiltonian, symmetry sector, or solver basis.

No production numerical routines are modified. Every accepted run must
match the independent stationary non-U(1) mixed-spin reference, both
local magnetization and the nine Cartesian nearest-neighbor correlations,
and the direct residual requirement. A failed calibration, nonconvergence
or unexpected capacity change is recorded and causes a nonzero exit.

CI-hosted timings are characterization only, not a universal speedup.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
from mixed_spin_dm_chirality_scan import reference
from mixed_spin_matched_benchmark import initial_vector_text,run_one
from cpc_stationary_resource_benchmark import environment_record,stats

def unique_csv(value,allow,parse=str):
    xs=[parse(v.strip()) for v in value.split(",")]
    if not xs or len(xs)!=len(set(xs)) or any(x not in allow for x in xs):
        raise ValueError("Invalid/duplicate CSV values: "+value)
    return xs

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="8,10")
    p.add_argument("--algs",default="1")
    p.add_argument("--D",type=float,default=0.12)
    p.add_argument("--omp",type=int,choices=(1,2),default=1)
    p.add_argument("--workspaces",default="2:16")
    p.add_argument("--repeats",type=int,default=2)
    p.add_argument("--warmups",type=int,default=1)
    p.add_argument("--residual-atol",type=float,default=1e-8)
    p.add_argument("--timeout-seconds",type=int,default=300)
    a=p.parse_args()
    try:
        sizes=unique_csv(a.sizes,{4,6,8,10,12},int)
        algs=unique_csv(a.algs,{1,2},int)
        parts=a.workspaces.split(":")
        if len(parts)!=2: raise ValueError("ALG=2 workspace should be NOK:NOM")
        keep,capacity=map(int,parts)
        if not 1<keep<capacity<=512: raise ValueError("Invalid ALG2 workspace")
        if not -0.5<=a.D<=0.5: raise ValueError("D outside tested interval")
        if not 0<a.residual_atol<=1e-5: raise ValueError("Invalid residual tolerance")
        if not 1<=a.repeats<=10 or not 0<=a.warmups<=5:
            raise ValueError("Invalid number of repetitions")
        if not 10<=a.timeout_seconds<=1800: raise ValueError("Invalid timeout")
    except ValueError as e:
        p.error(str(e))

    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(exist_ok=False,parents=True)
    raw=[]
    refs=[]
    capacity_calibration=[]
    violations=[]
    for n in sizes:
        r=reference(n,a.D,check_dense=(n<=6),profile="stationary")
        model=dict(spin=r["spins"],bonds=r["bonds"],fields=r["fields"],
                   pairs=r["pairs"],dimension=r["dimension"],
                   energies=r["energies"],mag=r["lm"],corr=r["corr"])
        initial=initial_vector_text(r["dimension"])
        initial_sha=hashlib.sha256(initial.encode()).hexdigest()
        # In this specific OBC model each S=1/2 site admits one
        # transverse single-site flip, each S=1 site admits up to two,
        # and every mixed-spin NN bond admits up to two joint flips.
        # Including the single diagonal row entry gives
        # 1 + (L/2) + 2*(L/2) + 2*(L-1) = 7L/2 - 1.
        analytical_mnte=7*(n//2)-1
        if (r["spins"]!=[.5,1.]*(n//2) or len(r["bonds"])!=n-1
                or any(i+1!=j for i,j,_ in r["bonds"])):
            raise AssertionError("Analytical MNTE formula requires the OBC alternating-spin chain")
        refs.append(dict(sites=n,dimension=r["dimension"],D=a.D,
                         ground_energy=r["E0"],first_gap=r["gap"],
                         independent_max_residual=max(r["reference_residual"]),
                         explicit_initial_sha256=initial_sha))
        # Always calibrate from the all-representatives cache (not a half
        # prefix). A run with MNTE=-1 covers every generated transition.
        baseline_optimal=None
        for alg in algs:
            for rep in range(1,a.warmups+a.repeats+1):
                phase="warmup" if rep<=a.warmups else "measurement"
                for label in ("default","calibrated"):
                    if label=="calibrated" and baseline_optimal is None:
                        violations.append(dict(n=n,alg=alg,rep=rep,
                                               error="No full-row capacity calibration"))
                        continue
                    assigned_mnte=-1 if label=="default" else baseline_optimal
                    group=root/f"n{n}-D{a.D:+.3f}-omp{a.omp}-alg{alg}-{label}"
                    group.mkdir(exist_ok=True)
                    try:
                        entry=run_one(exe,group,n,alg,a.residual_atol,rep,model,initial,
                                      keep=keep,capacity=capacity,threads=a.omp,
                                      wk_dim=-1,mnte=assigned_mnte,
                                      timeout_seconds=a.timeout_seconds)
                        observed=entry["reported_optimal_mnte"]
                        if observed<=0: raise AssertionError("Zero/negative max transition capacity")
                        if baseline_optimal is None:
                            if observed!=analytical_mnte:
                                raise AssertionError("Observed MNTE violates derived exact bound "+
                                                     str((n,observed,analytical_mnte)))
                            baseline_optimal=observed
                            if baseline_optimal>entry["reported_current_mnte"]:
                                raise AssertionError("Calibration exceeds default capacity")
                            capacity_calibration.append(dict(n=n,dimension=r["dimension"],
                                                             default_mnte=entry["reported_current_mnte"],
                                                             optimal_mnte=baseline_optimal,
                                                             analytical_mnte=analytical_mnte))
                        elif observed!=baseline_optimal:
                            raise AssertionError("Maximum row count changed for identical Hamiltonian "+
                                                 str((n,baseline_optimal,observed)))
                        if label=="calibrated" and entry["reported_current_mnte"]!=baseline_optimal:
                            raise AssertionError("Optimized MNTE did not take effect")
                        entry.update(phase=phase,capacity_policy=label,
                                     mnte_capacity=entry["reported_current_mnte"],
                                     observed_optimal_mnte=observed,
                                     initial_sha256=initial_sha)
                        raw.append(entry)
                        print("MNTE",phase,f"n{n}",f"omp{a.omp}",f"alg{alg}",label,
                              "rep",rep,"MNTE",entry["reported_current_mnte"],
                              "optimal",observed,"status",entry["status"],
                              "seconds",entry["wall_seconds"],
                              "Hcalls",entry["h_applications"],
                              "RSS_KiB",entry["peak_rss_kib"],flush=True)
                    except Exception as e:
                        issue=dict(n=n,alg=alg,rep=rep,policy=label,
                                   error=str(e),kind=type(e).__name__)
                        violations.append(issue)
                        print("MNTE CONTRACT FAILURE",issue,flush=True)
    summary=[]
    for n in sizes:
        calibr=next((x for x in capacity_calibration if x["n"]==n),None)
        for alg in algs:
            for policy in ("default","calibrated"):
                rows=[x for x in raw if x["n_sites"]==n and x["alg"]==alg and
                      x.get("capacity_policy")==policy and x["phase"]=="measurement" and
                      x["status"]==0]
                summary.append(dict(sites=n,dimension=6**(n//2),
                                    omp=a.omp,alg=alg,cache="full",
                                    capacity_policy=policy,
                                    MNTE=calibr["default_mnte" if policy=="default" else "optimal_mnte"]
                                         if calibr else None,
                                    successful_measurements=len(rows),
                                    target=a.repeats,
                                    wall_seconds=stats(x["wall_seconds"] for x in rows),
                                    H_applications=stats(x["h_applications"] for x in rows),
                                    peak_RSS_KiB=stats(x["peak_rss_kib"] for x in rows
                                                       if x["peak_rss_kib"] is not None),
                                    achieved_residual=stats(x["residual"] for x in rows)))
    expected=len(sizes)*len(algs)*2*(a.repeats+a.warmups)
    nonconverged=sum(x["status"]!=0 for x in raw)
    complete=(len(raw)==expected and not nonconverged and not violations and
              len(capacity_calibration)==len(sizes) and
              all(x["successful_measurements"]==a.repeats for x in summary))
    report=dict(schema=1,
                kind="full-cache representative-row MNTE allocation characterization",
                source_sha=os.environ.get("GITHUB_SHA","record checked-out SHA"),
                timing_warning="CI runner measurements are not controlled-node speedup claims",
                environment=environment_record(),
                params=dict(sizes=sizes,algs=algs,D=a.D,OMP=a.omp,
                            NOE=1,ALG2_workspace=[keep,capacity],
                            residual_atol=a.residual_atol,wk_dim=-1,
                            repeats=a.repeats,warmups=a.warmups,seed=19,
                            explicit_initial=True),
                independent_references=refs,
                capacity_calibration=capacity_calibration,raw=raw,summary=summary,
                expected=expected,completed=len(raw),successful=len(raw)-nonconverged,
                numerical_nonconvergence=nonconverged,contract_violations=violations,
                pass_complete=complete)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    (root/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("MNTE CAPACITY",len(raw),"/",expected,
          "nonconverged",nonconverged,"violations",len(violations),
          "complete",complete,flush=True)
    return 0 if complete else 1

if __name__=="__main__":
    raise SystemExit(main())
