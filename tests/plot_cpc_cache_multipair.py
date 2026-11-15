#!/usr/bin/env python3
"""Audit frozen CI measurements and make reproducible CPC resource plots.

Two input CSVs are *observations*, not authoritative physics references:
- 52 cache-policy trials (of which 12 small-system measurements per OMP
  plus 8 L=12 measurements per OMP are timed)
- 36 two-lowest-eigenpair restart trials.

The audit guards against missing rows, duplicated trials, warmup leakage
into medians, invalid measurements, and swapped cache/restart columns.
Scientific acceptance lives in the original independent ED2 CI jobs.
The figures are descriptive of shared hosted CI machines, NOT a controlled
speedup or strong-scaling study. No bulk extrapolation is performed.

Usage:
  python3 tests/plot_cpc_cache_multipair.py --validate-only
  python3 tests/plot_cpc_cache_multipair.py --outdir /new/empty/cpc-resource-figs
"""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import statistics

CACHE_FIELDS=(
    "workflow_run_id","job_id","source_sha","n_sites","D","omp",
    "phase","repeat","alg","cache_policy","h_applications",
    "solver_status","wall_seconds","peak_rss_kib"
)
MULTIPAIR_FIELDS=(
    "workflow_run_id","job_id","source_sha","n_sites","D","omp",
    "phase","repeat","NOE","NOK","NOM","solver_status",
    "h_applications","max_direct_residual","wall_seconds","peak_rss_kib"
)
CACHE_RUN_ID="43702889491"
MULTIPAIR_RUN_ID="43703014800"
CACHE_SHA="bdbb4402baeb842a80f5b39779150425ec9226b6"
MULTIPAIR_SHA="c89f018c2adef6af5c1647e83e2b4fed6dbfaa17"
DEFAULT_CACHE=Path(__file__).resolve().parents[1]/"docs/development/cpc-cache-sweep-raw-ci.csv"
DEFAULT_MULTIPAIR=Path(__file__).resolve().parents[1]/"docs/development/cpc-multipair-raw-ci.csv"


def read_csv(path,fields):
    if not path.is_file():
        raise ValueError("Missing input CSV: "+str(path))
    with path.open(newline="",encoding="utf-8") as fp:
        reader=csv.DictReader(fp)
        if reader.fieldnames!=list(fields):
            raise ValueError("Header changed: "+str((reader.fieldnames,fields)))
        rows=list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError("Malformed or empty input CSV: "+str(path))
    return rows


def positive_number(v,context,zero_ok=False):
    try: x=float(v)
    except (ValueError,TypeError): raise ValueError("Non-numeric "+context) from None
    if not math.isfinite(x) or (x<0 if zero_ok else x<=0):
        raise ValueError("Nonfinite/nonpositive "+context+": "+str(v))
    return x


def validate_cache(rows):
    if len(rows)!=52: raise ValueError("Expected 52 cache rows, got "+str(len(rows)))
    grouped=defaultdict(list)
    unique=set()
    for r in rows:
        if r["workflow_run_id"]!=CACHE_RUN_ID or r["source_sha"]!=CACHE_SHA:
            raise ValueError("Cache CI provenance changed")
        n=int(r["n_sites"]); omp=int(r["omp"]); alg=int(r["alg"])
        rep=int(r["repeat"])
        pol=r["cache_policy"]; phase=r["phase"]
        if n not in (8,10,12) or omp not in (1,2) or pol not in ("none","half","full"):
            raise ValueError("Unexpected cache case")
        if alg!=(1 if n<12 else alg) or alg not in ((1,) if n<12 else (1,2)):
            raise ValueError("Invalid solver in cache case")
        if n==12 and pol=="none": raise ValueError("No-cache L=12 not part of pilot")
        if phase not in ("warmup","measurement") or r["solver_status"]!="0":
            raise ValueError("Cache trial did not pass")
        if abs(float(r["D"])-0.12)>1e-12:
            raise ValueError("Cache Hamiltonian changed")
        positive_number(r["wall_seconds"],"cache wall time")
        positive_number(r["peak_rss_kib"],"cache RSS")
        positive_number(r["h_applications"],"cache H calls")
        if not r["job_id"].isdigit() or not r["workflow_run_id"].isdigit():
            raise ValueError("Invalid immutable job reference")
        g=(n,omp,alg,pol)
        key=(g,phase,rep)
        if key in unique: raise ValueError("Duplicate cache trial: "+str(key))
        unique.add(key)
        grouped[g].append(r)
    required=set()
    for n in (8,10):
        for omp in (1,2):
            for pol in ("none","half","full"):
                required.add((n,omp,1,pol))
    for omp in (1,2):
        for alg in (1,2):
            for pol in ("half","full"):
                required.add((12,omp,alg,pol))
    if set(grouped)!=required:
        raise ValueError("Missing or extra cache workload keys "+
                         str((required-set(grouped),set(grouped)-required)))
    for key,items in grouped.items():
        n=key[0]
        phases={(r["phase"],int(r["repeat"])) for r in items}
        expected=({("warmup",1),("measurement",2),("measurement",3)}
                  if n<12 else {("measurement",1),("measurement",2)})
        if phases!=expected: raise ValueError("Cache warmup/repeat mismatch: "+str((key,phases)))
        # All cache measurements at each L/OMP for ALG1 should have identical H calls.
    for n in (8,10):
        for omp in (1,2):
            calls={int(r["h_applications"]) for r in rows if int(r["n_sites"])==n
                   and int(r["omp"])==omp and int(r["alg"])==1}
            if len(calls)!=1:
                raise ValueError("ALG1 cache policy changed H-application count")
    return grouped


def validate_multipair(rows):
    if len(rows)!=36: raise ValueError("Expected 36 multipair rows, got "+str(len(rows)))
    grouped=defaultdict(list)
    seen=set()
    for r in rows:
        if r["workflow_run_id"]!=MULTIPAIR_RUN_ID or r["source_sha"]!=MULTIPAIR_SHA:
            raise ValueError("Multipair CI provenance changed")
        n=int(r["n_sites"]); omp=int(r["omp"])
        D=float(r["D"]); k=int(r["NOK"]); m=int(r["NOM"])
        rep=int(r["repeat"])
        if n not in (8,10,12) or omp not in (1,2) or (k,m) not in ((2,16),(4,24)):
            raise ValueError("Unexpected multipair setting")
        if r["NOE"]!="2" or r["phase"]!="measurement" or r["solver_status"]!="0":
            raise ValueError("Multipair result not accepted")
        if D not in (0.,0.12) or (n==12 and D!=0.12):
            raise ValueError("Unexpected DM parameter in multipair")
        if positive_number(r["max_direct_residual"],"direct residual",zero_ok=True)>1e-8*(1+1e-10):
            raise ValueError("Accepted residual exceeds requested limit")
        positive_number(r["h_applications"],"multipair H calls")
        positive_number(r["wall_seconds"],"multipair wall time")
        positive_number(r["peak_rss_kib"],"multipair RSS")
        if not r["job_id"].isdigit(): raise ValueError("Missing job identifier")
        key=(n,D,omp,k,m)
        trial=(key,rep)
        if trial in seen: raise ValueError("Duplicate multipair trial")
        seen.add(trial)
        grouped[key].append(r)
    required={(n,D,omp,k,m)
              for n in (8,10,12)
              for D in ((0.,0.12) if n<12 else (0.12,))
              for omp in (1,2)
              for k,m in ((2,16),(4,24))}
    if set(grouped)!=required:
        raise ValueError("Missing/extra multipair workloads")
    for key,items in grouped.items():
        expected={1,2} if key[0]<12 else {1}
        if {int(r["repeat"]) for r in items}!=expected:
            raise ValueError("Missing multipair repeat "+str(key))
    return grouped


def median(items,col):
    vals=[float(x[col]) for x in items if x["phase"]=="measurement"]
    if not vals: raise ValueError("No timed replicates for "+col)
    return statistics.median(vals)


def build_summary(cache,multi):
    rows=[]
    for (n,omp,alg,policy),items in sorted(cache.items()):
        count=len([x for x in items if x["phase"]=="measurement"])
        rows.append(dict(kind="cache",n_sites=n,D=0.12,omp=omp,alg=alg,
                         case=policy,wk_fraction={"none":0.,"half":0.5,"full":1.}[policy],
                         measured_repeats=count,
                         wall_median_s=median(items,"wall_seconds"),
                         rss_median_mib=median(items,"peak_rss_kib")/1024,
                         h_median=median(items,"h_applications"),max_residual=""))
    for (n,D,omp,k,m),items in sorted(multi.items()):
        rows.append(dict(kind="multipair",n_sites=n,D=D,omp=omp,alg=2,
                         case=f"NOK={k};NOM={m}",wk_fraction=1.,
                         measured_repeats=len(items),
                         wall_median_s=median(items,"wall_seconds"),
                         rss_median_mib=median(items,"peak_rss_kib")/1024,
                         h_median=median(items,"h_applications"),
                         max_residual=max(float(x["max_direct_residual"]) for x in items)))
    return rows


def make_plots(summary,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=False)
    plots=[]
    for column,suffix,ylabel,logscale in [
        ("wall_median_s","cache_time","ED2 process wall time (s)",True),
        ("rss_median_mib","cache_memory","Peak ED2 RSS (MiB)",False)]:
        fig,ax=plt.subplots(figsize=(6.1,4.0),layout="constrained")
        for n in (8,10,12):
            items=sorted((x for x in summary if x["kind"]=="cache" and
                          x["n_sites"]==n and x["omp"]==1 and x["alg"]==1),
                         key=lambda x:x["wk_fraction"])
            ax.plot([x["wk_fraction"] for x in items],
                    [x[column] for x in items],marker="o",label=f"L={n}")
        ax.set(xlabel="Fraction of representative states cached",ylabel=ylabel)
        ax.set_xticks([0,.5,1])
        if logscale: ax.set_yscale("log")
        ax.legend(title="OMP=1, ALG=1")
        ax.grid(True,alpha=.25)
        for ext in (".svg",".png"):
            path=out/(suffix+ext)
            fig.savefig(path,dpi=230)
            plots.append(path)
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.1,4.0),layout="constrained")
    for n in (8,10,12):
        items=sorted((x for x in summary if x["kind"]=="multipair" and
                      x["n_sites"]==n and x["omp"]==1 and x["D"]==0.12),
                     key=lambda x:x["h_median"])
        ax.plot([x["h_median"] for x in items],
                [x["max_residual"] for x in items],marker="o",label=f"L={n}")
        for x in items:
            ax.annotate(x["case"].replace("NOK=","k").replace(";NOM=","/"),
                        (x["h_median"],x["max_residual"]),
                        xytext=(3,4),textcoords="offset points",fontsize=8)
    ax.axhline(1e-8,ls="--",linewidth=1,label="Requested residual")
    ax.set_yscale("log")
    ax.set(xlabel="Hamiltonian applications (NOE=2, ALG=2)",
           ylabel="Maximum direct eigenpair residual")
    ax.legend(title="OMP=1, D=0.12",fontsize=8)
    ax.grid(True,alpha=.25)
    for ext in (".svg",".png"):
        path=out/("multipair_cost_residual"+ext)
        fig.savefig(path,dpi=230)
        plots.append(path)
    plt.close(fig)
    return plots


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-csv",type=Path,default=DEFAULT_CACHE)
    p.add_argument("--multipair-csv",type=Path,default=DEFAULT_MULTIPAIR)
    p.add_argument("--outdir",type=Path)
    p.add_argument("--validate-only",action="store_true")
    a=p.parse_args()
    cache=validate_cache(read_csv(a.cache_csv,CACHE_FIELDS))
    multi=validate_multipair(read_csv(a.multipair_csv,MULTIPAIR_FIELDS))
    summary=build_summary(cache,multi)
    assert len(summary)==40,(len(summary),len(cache),len(multi))
    print("SOURCE DATA VERIFIED: 52 cache trials, 36 two-eigenpair trials,",
          len(summary),"summarized workload cells",flush=True)
    if a.validate_only:
        if a.outdir: p.error("--validate-only cannot also request --outdir")
        return 0
    if a.outdir is None: p.error("--outdir required to render figures")
    out=a.outdir.resolve()
    paths=make_plots(summary,out)
    with (out/"cpc-resource-summary.csv").open("w",newline="",encoding="utf-8") as fp:
        columns=list(summary[0])
        writer=csv.DictWriter(fp,fieldnames=columns)
        writer.writeheader()
        writer.writerows(summary)
    print("RENDERED",len(paths),"figures:",", ".join(x.name for x in paths),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
