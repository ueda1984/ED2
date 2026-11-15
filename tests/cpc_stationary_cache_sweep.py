#!/usr/bin/env python3
"""Cache-policy crossover for the independently verified stationary mixed-spin ED2 model.

This test changes *only* the Hamiltonian-action caching policy, wk_dim,
for a fixed physical Hamiltonian, eigensolver, explicit initial state,
requested accuracy and observable set. The cache fraction is the number
of cached representative states divided by the full working dimension.

Use small sizes to explore wk_dim=0 (no cache) before extending sizes.
All reported results are independent physical comparisons; timings
on shared GitHub Actions runners are exploratory, not isolated-node
speedup measurements. Every subprocess has its own clean directory.
"""
import argparse
import json
import os
from pathlib import Path
import statistics

from mixed_spin_dm_chirality_scan import reference
from mixed_spin_matched_benchmark import initial_vector_text, run_one
from cpc_stationary_resource_benchmark import environment_record, stats

POLICIES=("none","quarter","half","full")


def cache_count(policy,dimension):
    if policy=="none":
        return 0
    if policy=="quarter":
        return dimension//4
    if policy=="half":
        return dimension//2
    if policy=="full":
        return -1
    raise ValueError("Unknown cache policy "+policy)


def unique_items(text,allowed,name):
    items=[x.strip() for x in text.split(",")]
    if not items or any(not x or x not in allowed for x in items) or len(set(items))!=len(items):
        raise ValueError("Invalid or duplicated "+name+": "+text)
    return items


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--sizes",default="8,10")
    ap.add_argument("--policies",default="none,half,full")
    ap.add_argument("--algs",default="1")
    ap.add_argument("--D",default=0.12,type=float)
    ap.add_argument("--omp",type=int,default=1,choices=(1,2))
    ap.add_argument("--residual-atol",type=float,default=1e-8)
    ap.add_argument("--repeats",type=int,default=2)
    ap.add_argument("--warmups",type=int,default=0)
    ap.add_argument("--alg2-keep",type=int,default=2)
    ap.add_argument("--alg2-capacity",type=int,default=16)
    ap.add_argument("--timeout-seconds",type=int,default=300)
    a=ap.parse_args()
    try:
        sizes=[int(x) for x in unique_items(a.sizes,{"4","6","8","10","12"},"sizes")]
        policies=unique_items(a.policies,set(POLICIES),"policies")
        algs=[int(x) for x in unique_items(a.algs,{"1","2"},"algs")]
    except ValueError as ex:
        ap.error(str(ex))
    if not -0.5<=a.D<=0.5 or not 0<a.residual_atol<=1e-5:
        ap.error("D or residual tolerance is outside accepted control range")
    if not 1<=a.repeats<=20 or not 0<=a.warmups<=5 or not 10<=a.timeout_seconds<=1800:
        ap.error("Invalid number of repeats, warmups or subprocess timeout")
    if not 1<a.alg2_keep<a.alg2_capacity<=512:
        ap.error("Invalid ALG=2 restart workspace")
    if 12 in sizes and "none" in policies:
        ap.error("wk_dim=0 at L=12 is intentionally not part of the controlled pilot; use smaller sizes")
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    entries=[]
    references=[]
    violations=[]
    for n in sizes:
        r=reference(n,a.D,check_dense=(n<=6),profile="stationary")
        model=dict(spin=r["spins"],bonds=r["bonds"],fields=r["fields"],
                   pairs=r["pairs"],dimension=r["dimension"],
                   energies=r["energies"],mag=r["lm"],corr=r["corr"])
        init=initial_vector_text(r["dimension"])
        references.append(dict(sites=n,dim=r["dimension"],D=a.D,E0=r["E0"],
                               first_gap=r["gap"],
                               reference_max_residual=max(r["reference_residual"])))
        # Interleave by run index, then cache policy and solver.
        for rep in range(1,a.warmups+a.repeats+1):
            phase="warmup" if rep<=a.warmups else "measurement"
            for policy in policies:
                wk_dim=cache_count(policy,r["dimension"])
                for alg in algs:
                    group=root/f"n{n}-D{a.D:+.3f}-omp{a.omp}-cache{policy}-alg{alg}"
                    group.mkdir(exist_ok=True)
                    try:
                        row=run_one(exe,group,n,alg,a.residual_atol,rep,model,init,
                                    keep=a.alg2_keep,capacity=a.alg2_capacity,
                                    threads=a.omp,wk_dim=wk_dim,
                                    timeout_seconds=a.timeout_seconds)
                        row.update(D=a.D,cache_policy=policy,phase=phase,
                                   cache_fraction=(r["dimension"] if wk_dim<0 else wk_dim)
                                                  /r["dimension"])
                        entries.append(row)
                        print("CACHE",phase,row["case"],
                              "policy",policy,"alg",alg,"status",row["status"],
                              "seconds",row["wall_seconds"],
                              "Hcalls",row["h_applications"],
                              "RSS KiB",row["peak_rss_kib"],flush=True)
                    except Exception as exc:
                        item=dict(sites=n,policy=policy,alg=alg,rep=rep,phase=phase,
                                  error=str(exc),error_class=type(exc).__name__)
                        violations.append(item)
                        print("CACHE CONTRACT VIOLATION",item,flush=True)
    summary=[]
    for n in sizes:
        for policy in policies:
            for alg in algs:
                group=[x for x in entries if x["n_sites"]==n and x["alg"]==alg and
                       x["cache_policy"]==policy and x["phase"]=="measurement" and
                       x["status"]==0]
                summary.append(dict(sites=n,dimension=next(x["dim"] for x in references
                                                            if x["sites"]==n),
                                    cache_policy=policy,alg=alg,OMP=a.omp,
                                    wk_dim_requested=cache_count(policy,next(
                                        x["dim"] for x in references if x["sites"]==n)),
                                    measurements=len(group),target=a.repeats,
                                    wall_seconds=stats(x["wall_seconds"] for x in group),
                                    H_applications=stats(x["h_applications"] for x in group),
                                    peak_RSS_KiB=stats(x["peak_rss_kib"] for x in group
                                                       if x["peak_rss_kib"] is not None),
                                    achieved_residual=stats(x["residual"] for x in group)))
    expected=len(sizes)*len(policies)*len(algs)*(a.warmups+a.repeats)
    nonconverged=sum(x["status"]!=0 for x in entries)
    complete=(len(entries)==expected and not nonconverged and not violations
              and all(x["measurements"]==a.repeats for x in summary))
    report=dict(schema=1,
                kind="stationary non-U1 mixed-spin DM, matched-target wk_dim cache sweep",
                source_sha=os.environ.get("GITHUB_SHA","record git checkout"),
                timing_limitation="CI timings are characterization only; not a dedicated-node comparison",
                environment=environment_record(),
                params=dict(sizes=sizes,D=a.D,algs=algs,OMP=a.omp,
                            policies=policies,alg2_keep=a.alg2_keep,
                            alg2_capacity=a.alg2_capacity,absolute_residual=a.residual_atol,
                            repeats=a.repeats,warmups=a.warmups,timeout_seconds=a.timeout_seconds,
                            seed=19,explicit_initial=True,openblas_threads=1,mkl_threads=1),
                references=references,raw=entries,summary=summary,
                expected=expected,completed=len(entries),
                accepted=len(entries)-nonconverged,nonconverged=nonconverged,
                contract_violations=violations,pass_complete=complete)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    (root/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("CACHE SWEEP",len(entries),"/",expected,
          "accepted",len(entries)-nonconverged,
          "nonconverged",nonconverged,"violations",len(violations),flush=True)
    return 0 if complete else 1


if __name__=="__main__":
    raise SystemExit(main())
