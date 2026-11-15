#!/usr/bin/env python3
"""A/B cache-memory frontier for the same stationary non-U(1) spin problem.

Compare three policies *within the same CI job* and identical physical
workload, residual request and explicit solver starting vector:
  half_generic: wk_dim=D//2, MNTE=-1 (generic wide rows)
  full_tight:   wk_dim=D, MNTE=7L/2-1 (derived exact model row bound)
  full_generic: wk_dim=D, MNTE=-1 (generic wide rows)

For even-L open alternating S=1/2,1 chains with NN XYZ/Gamma_xy/Dz
and nonzero transverse fields, the tight bound is 7L/2-1 as justified
in docs/development/cpc-mnte-capacity-evaluation.md. This formula
is NOT a safe general ED2 default and must not be reused blindly.

Every case verifies the ground energy, all local magnetizations,
all nine NN Cartesian correlations and direct residual against an
independent SciPy reference. Status, achieved residual and case outputs
are never silently ignored. These CI timings remain shared-runner
characterization, not dedicated-node speedup claims.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

from mixed_spin_dm_chirality_scan import reference
from mixed_spin_matched_benchmark import initial_vector_text, run_one
from cpc_stationary_resource_benchmark import environment_record,stats

POLICIES=("half_generic","full_tight","full_generic")


def policy_settings(n,dim,policy):
    default=10*n-8
    tight=7*(n//2)-1
    if policy=="half_generic": return dim//2,-1,20*(default+1)*(dim//2)
    if policy=="full_tight": return -1,tight,20*(tight+1)*dim
    if policy=="full_generic": return -1,-1,20*(default+1)*dim
    raise ValueError("Unknown policy")


def comma_ints(value,allow):
    items=[int(s.strip()) for s in value.split(",")]
    if not items or len(set(items))!=len(items) or any(x not in allow for x in items):
        raise ValueError("Invalid repeated/unsupported sizes or solvers")
    return items


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--sizes",default="8,10")
    p.add_argument("--algs",default="1,2")
    p.add_argument("--D",type=float,default=0.12)
    p.add_argument("--omp",type=int,choices=(1,2),default=1)
    p.add_argument("--repeats",type=int,default=2)
    p.add_argument("--warmups",type=int,default=1)
    p.add_argument("--residual-atol",type=float,default=1e-8)
    p.add_argument("--timeout-seconds",type=int,default=350)
    p.add_argument("--keep",type=int,default=2)
    p.add_argument("--capacity",type=int,default=16)
    a=p.parse_args()
    try:
        sizes=comma_ints(a.sizes,{4,6,8,10,12})
        algs=comma_ints(a.algs,{1,2})
        if not -0.5<=a.D<=0.5 or not 0<a.residual_atol<=1e-5:
            raise ValueError("D or residual tolerance outside tested range")
        if not 1<=a.repeats<=10 or not 0<=a.warmups<=5:
            raise ValueError("Invalid repeat controls")
        if not 1<a.keep<a.capacity<=512 or not 10<=a.timeout_seconds<=1800:
            raise ValueError("Invalid workspace or timeout")
    except ValueError as exc:
        p.error(str(exc))
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    entries=[]
    references=[]
    failures=[]
    for n in sizes:
        r=reference(n,a.D,check_dense=(n<=6),profile="stationary")
        dim=r["dimension"]
        if r["spins"]!=[0.5,1.]*(n//2) or len(r["bonds"])!=n-1:
            raise AssertionError("This memory-bound test requires the original OBC chain")
        model=dict(spin=r["spins"],bonds=r["bonds"],fields=r["fields"],
                   pairs=r["pairs"],dimension=dim,
                   energies=r["energies"],mag=r["lm"],corr=r["corr"])
        initial=initial_vector_text(dim)
        initial_sha=hashlib.sha256(initial.encode()).hexdigest()
        budget={name:policy_settings(n,dim,name)[2] for name in POLICIES}
        if not budget["full_tight"]<budget["half_generic"]<budget["full_generic"]:
            raise AssertionError("Expected cache-memory ordering not satisfied")
        references.append(dict(sites=n,dim=dim,D=a.D,E0=r["E0"],
                               independent_max_residual=max(r["reference_residual"]),
                               exact_mnte=7*(n//2)-1,
                               generic_mnte=10*n-8,
                               nominal_cache_bytes=budget,initial_sha256=initial_sha))
        for rep in range(1,a.warmups+a.repeats+1):
            phase="warmup" if rep<=a.warmups else "measurement"
            for alg in algs:
                for policy in POLICIES:
                    wk_dim,mnte,nominal=policy_settings(n,dim,policy)
                    group=root/f"n{n}-D{a.D:+.3f}-omp{a.omp}-alg{alg}-{policy}"
                    group.mkdir(exist_ok=True)
                    try:
                        entry=run_one(exe,group,n,alg,a.residual_atol,rep,
                                      model,initial,keep=a.keep,capacity=a.capacity,
                                      threads=a.omp,wk_dim=wk_dim,mnte=mnte,
                                      timeout_seconds=a.timeout_seconds)
                        if policy in ("full_tight","full_generic"):
                            if entry["reported_optimal_mnte"]!=7*(n//2)-1:
                                raise AssertionError("Source enumerator contradicts analytic row bound")
                        if entry["reported_current_mnte"] != (10*n-8 if mnte<0 else mnte):
                            raise AssertionError("Requested capacity was not applied")
                        entry.update(phase=phase,policy=policy,D=a.D,
                                     nominal_cache_bytes=nominal,
                                     initial_sha256=initial_sha)
                        entries.append(entry)
                        print("PARETO",phase,"n",n,"omp",a.omp,"alg",alg,
                              "policy",policy,"rep",rep,
                              "status",entry["status"],"Hcalls",entry["h_applications"],
                              "r",entry.get("residual"),"wall",entry["wall_seconds"],
                              "rssKiB",entry["peak_rss_kib"],flush=True)
                    except Exception as exc:
                        issue=dict(n=n,omp=a.omp,alg=alg,policy=policy,rep=rep,
                                   error=str(exc),kind=type(exc).__name__)
                        failures.append(issue)
                        print("PARETO CONTRACT FAILURE",issue,flush=True)
    summary=[]
    for n in sizes:
        for alg in algs:
            for policy in POLICIES:
                samples=[x for x in entries if x["n_sites"]==n and x["alg"]==alg
                         and x["policy"]==policy and x["phase"]=="measurement"
                         and x["status"]==0]
                summary.append(dict(n_sites=n,omp=a.omp,alg=alg,policy=policy,
                                    measurements=len(samples),
                                    expected_measurements=a.repeats,
                                    nominal_cache_bytes=policy_settings(n,6**(n//2),policy)[2],
                                    wall_seconds=stats(x["wall_seconds"] for x in samples),
                                    H_applications=stats(x["h_applications"] for x in samples),
                                    max_RSS_KiB=stats(x["peak_rss_kib"] for x in samples
                                                      if x["peak_rss_kib"] is not None),
                                    achieved_residual=stats(x["residual"] for x in samples)))
    expected=len(sizes)*len(algs)*len(POLICIES)*(a.warmups+a.repeats)
    nonconverged=sum(x["status"]!=0 for x in entries)
    passed=len(entries)==expected and not failures and not nonconverged and all(
        x["measurements"]==a.repeats for x in summary)
    report=dict(schema=1,scope="paired same-node-run full-tight versus half-generic cache",
                shared_runner_caution="Not a dedicated-node speedup, measurements may have load variation",
                source_sha=os.environ.get("GITHUB_SHA","record checkout SHA"),
                environment=environment_record(),
                parameters=dict(sizes=sizes,algs=algs,D=a.D,omp=a.omp,
                                NOE=1,residual_atol=a.residual_atol,
                                keep=a.keep,capacity=a.capacity,
                                warmups=a.warmups,repeats=a.repeats,
                                seed=19,explicit_initial=True),
                references=references,raw=entries,summary=summary,
                completed=len(entries),expected=expected,
                accepted=len(entries)-nonconverged,
                numerical_nonconvergence=nonconverged,
                violations=failures,pass_complete=passed)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    (root/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("CACHE PARETO",len(entries),"/",expected,
          "nonconverged",nonconverged,"violations",len(failures),
          "complete",passed,flush=True)
    return 0 if passed else 1

if __name__=="__main__":
    raise SystemExit(main())
