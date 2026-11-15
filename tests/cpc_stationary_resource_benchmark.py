#!/usr/bin/env python3
"""Repeated matched-target resources for ED2's two-site stationary DM model.

Uses the *same* independently validated non-U(1) mixed-spin model as the
CPC DM/chirality application, rather than silently switching to an easier
conserving chain. Both algorithms request exactly one ground eigenpair,
an identical explicit starting vector, the same absolute residual
threshold and the same physical outputs. The cached operator policy
wk_dim=-1 is held fixed; ALG=2 workspaces are reported, not conflated.

Each case writes fresh input/output directories, raw stdout, solver
diagnostics, wall time, peak child-process RSS and H-application counts.
Wall timings on GitHub-hosted CI runners are smoke/characterization data,
NOT publication-grade isolated-node performance figures.
"""
import argparse
import json
import os
import platform
from pathlib import Path
import shutil
import statistics
import subprocess
import sys

from mixed_spin_dm_chirality_scan import reference
from mixed_spin_matched_benchmark import initial_vector_text,run_one

def parse_workspaces(value):
    out=[]
    for token in value.split(","):
        bits=token.strip().split(":")
        if len(bits)!=2:
            raise ValueError("Expected --workspaces NOK:NOM[,NOK:NOM...]")
        keep,capacity=map(int,bits)
        if not 1<keep<capacity or capacity>512 or (keep,capacity) in out:
            raise ValueError("Invalid or duplicate ALG=2 workspace: "+token)
        out.append((keep,capacity))
    if not out:
        raise ValueError("At least one ALG=2 workspace required")
    return out

def environment_record():
    def probe(cmd):
        try:
            return subprocess.run(cmd,stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT,text=True,
                                  timeout=10,check=False).stdout.strip()
        except (OSError,subprocess.TimeoutExpired):
            return None
    return dict(
        platform=platform.platform(),
        python=sys.version.splitlines()[0],
        machine=platform.machine(),
        cpu_count=os.cpu_count(),
        cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else None,
        cpu_details=probe(["lscpu"]),
        gfortran=probe(["gfortran","--version"]),
        compiler_13=probe(["gfortran-13","--version"]),
        numpy_blas_caveat="BLAS thread counts are controlled in ED2 child only",
        provenance_sha=os.environ.get("GITHUB_SHA"),
        binary=str(shutil.which("gfortran") or "compiler not on PATH"),
    )

def stats(values):
    if not values:
        return None
    vals=[float(x) for x in values]
    med=statistics.median(vals)
    return dict(n=len(vals),median=med,min=min(vals),max=max(vals),
                median_absolute_deviation=statistics.median(abs(x-med) for x in vals))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--sizes",default="8,10")
    p.add_argument("--D",type=float,default=0.12)
    p.add_argument("--omp",type=int,choices=(1,2),default=1)
    p.add_argument("--residual-atol",type=float,default=1e-8)
    p.add_argument("--repeats",type=int,default=3)
    p.add_argument("--warmups",type=int,default=1)
    p.add_argument("--workspaces",default="2:16,4:24,10:64")
    a=p.parse_args()
    sizes=[int(x.strip()) for x in a.sizes.split(",")]
    if not sizes or len(set(sizes))!=len(sizes) or any(x not in (4,6,8,10,12) for x in sizes):
        p.error("Sizes must be distinct members of 4,6,8,10,12")
    if not (-0.5<=a.D<=0.5):
        p.error("D outside tested model range [-0.5,0.5]")
    if a.residual_atol<=0 or a.residual_atol>1e-5:
        p.error("Specify a positive, meaningful absolute residual")
    if not 1<=a.repeats<=20 or not 0<=a.warmups<=5:
        p.error("Invalid measurement/warmup repetition count")
    workspaces=parse_workspaces(a.workspaces)
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    # Distinct from a publication benchmark's isolated-node conditions:
    raw=[]
    violations=[]
    references=[]
    for n in sizes:
        r=reference(n,a.D,check_dense=(n<=6),profile="stationary")
        model=dict(spin=r["spins"],bonds=r["bonds"],fields=r["fields"],
                   pairs=r["pairs"],dimension=r["dimension"],
                   energies=r["energies"],mag=r["lm"],corr=r["corr"])
        initial=initial_vector_text(model["dimension"])
        references.append(dict(n_sites=n,D=a.D,dimension=r["dimension"],
                               E0=r["E0"],first_gap=r["gap"],
                               max_independent_residual=max(r["reference_residual"]),
                               chirality_per_bond=r["staggered_kappa"]/(n-1),
                               nonU1_relative_commutator=r["nonU1_commutator_relative"]))
        configs=[(1,None,None)]+[(2,k,c) for k,c in workspaces]
        # Interleave by repetition to limit systematic drift between solvers.
        for rep in range(1,a.warmups+a.repeats+1):
            phase="warmup" if rep<=a.warmups else "measurement"
            for alg,keep,capacity in configs:
                name=f"n{n}-D{a.D:+.3f}-omp{a.omp}-alg{alg}-"+(
                    f"{keep}x{capacity}" if alg==2 else "lanczos")
                group=root/name
                group.mkdir(exist_ok=True)
                try:
                    row=run_one(exe,group,n,alg,a.residual_atol,rep,
                                model,initial,keep=keep or 10,
                                capacity=capacity or 64,threads=a.omp)
                    row.update(phase=phase,workspace=name,
                               D=a.D,OMP=a.omp)
                    raw.append(row)
                    print("BENCH",phase,name,rep,
                          "status",row["status"],
                          "seconds",row["wall_seconds"],
                          "Hcalls",row["h_applications"],
                          "RSS KiB",row["peak_rss_kib"],flush=True)
                except Exception as exc:
                    entry=dict(case=name,rep=rep,phase=phase,error=str(exc))
                    violations.append(entry)
                    print("BENCH CONTRACT FAILURE",entry,flush=True)
    summary=[]
    for n in sizes:
        for alg,keep,capacity in [(1,None,None)]+[(2,k,c) for k,c in workspaces]:
            group=[x for x in raw if x["n_sites"]==n and
                   x["alg"]==alg and x.get("requested_keep")==keep and
                   x.get("requested_capacity")==capacity and
                   x["phase"]=="measurement" and x["status"]==0]
            summary.append(dict(n_sites=n,D=a.D,omp=a.omp,
                                alg=alg,NOE=1,NOK=keep,NOM=capacity,
                                requested_residual_atol=a.residual_atol,
                                successful_measurements=len(group),
                                wall_seconds=stats(x["wall_seconds"] for x in group),
                                H_applications=stats(x["h_applications"] for x in group),
                                peak_RSS_KiB=stats(x["peak_rss_kib"] for x in group
                                                   if x["peak_rss_kib"] is not None),
                                achieved_residual=stats(x["residual"] for x in group)))
    expected=len(sizes)*(len(workspaces)+1)*(a.repeats+a.warmups)
    nonconverged=sum(x["status"]!=0 for x in raw)
    complete=len(raw)==expected and nonconverged==0 and not violations
    report=dict(schema=1,scope="same-target stationary two-site mixed-spin DM resource study",
                timings_warning="CI runner wall time is NOT an isolated-node performance claim",
                source_sha=os.environ.get("GITHUB_SHA","record revision from git checkout"),
                env=environment_record(),
                ED2_runtime_environment=dict(
                    OMP_NUM_THREADS=a.omp,OMP_DYNAMIC="FALSE",
                    OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,
                    ED2_RNG_SEED=19,explicit_start=True,
                    wk_dim=-1,cal_lm=1,cal_cf=1,NOE=1),
                parameters=dict(sizes=sizes,D=a.D,omp=a.omp,
                                residual_atol=a.residual_atol,
                                repeats=a.repeats,warmups=a.warmups,
                                workspaces=workspaces),
                references=references,raw=raw,summary=summary,
                expected_runs=expected,completed=len(raw),
                accepted=len(raw)-nonconverged,
                honest_nonconvergence=nonconverged,
                contract_violations=violations,pass_complete=complete)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    (root/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("STATIONARY MATCHED BENCHMARK:",len(raw),"/",expected,
          "nonconvergence",nonconverged,"violations",len(violations),
          flush=True)
    return 0 if complete else 1

if __name__=="__main__":
    raise SystemExit(main())
