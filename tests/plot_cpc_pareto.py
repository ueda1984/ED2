#!/usr/bin/env python3
"""Audit frozen paired cache-policy CI trials and render a resource-frontier plot.

The same hosted CI job evaluates three policies at identical physical
Hamiltonian, requested ground eigenpair, initial vector, residual target,
local magnetizations and NN Cartesian correlations:
  half_generic (wk_dim=D/2, generic MNTE)
  full_tight   (wk_dim=D, exact model-specific MNTE)
  full_generic (wk_dim=D, generic MNTE)
The original CI scripts, not this plotting tool, check physical correctness.
Shared GitHub runner wall timings are descriptive, NOT a dedicated-node
or general-speedup claim. Warmups are excluded from measured medians.
"""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import statistics

DEFAULT=Path(__file__).resolve().parents[1]/"docs/development/cpc-cache-pareto-raw-ci.csv"
FIELDS=("workflow_run_id","job_id","source_sha","n_sites","D","omp","alg",
        "phase","repeat","policy","status","h_applications",
        "max_direct_residual","wall_seconds","peak_rss_kib")
RUN="43704096102"
SOURCE="08037050aac4cdaa0047c8c977edaaf9fc08ee46"
POLICIES=("half_generic","full_tight","full_generic")


def nominal_cache_bytes(n,pol):
    D=6**(n//2)
    if pol=="half_generic": return 20*(10*n-7)*(D//2)
    if pol=="full_tight": return 20*(7*(n//2))*D
    if pol=="full_generic": return 20*(10*n-7)*D
    raise ValueError("Unknown policy "+str(pol))


def audit(path):
    with path.open(newline="",encoding="utf-8") as fp:
        reader=csv.DictReader(fp)
        if reader.fieldnames!=list(FIELDS):
            raise ValueError("Unexpected source data columns")
        rows=list(reader)
    if len(rows)!=108 or any(None in r for r in rows):
        raise ValueError("Expected 108 complete frozen resource cases")
    groups=defaultdict(list)
    seen=set()
    for row in rows:
        if row["workflow_run_id"]!=RUN or row["source_sha"]!=SOURCE:
            raise ValueError("Invalid source provenance")
        n=int(row["n_sites"]); omp=int(row["omp"]); alg=int(row["alg"])
        rep=int(row["repeat"]); phase=row["phase"]; pol=row["policy"]
        if n not in (8,10,12) or omp not in (1,2) or alg not in (1,2):
            raise ValueError("Unsupported benchmark configuration")
        if pol not in POLICIES or row["status"]!="0":
            raise ValueError("Policy or solver status invalid")
        if phase!=("warmup" if rep==1 else "measurement") or rep not in (1,2,3):
            raise ValueError("Invalid warmup/repetition label")
        if abs(float(row["D"])-.12)>1e-12:
            raise ValueError("Physical Hamiltonian parameter changed")
        for f in ("h_applications","wall_seconds","peak_rss_kib"):
            v=float(row[f])
            if not math.isfinite(v) or v<=0:
                raise ValueError("Invalid measured "+f)
        r=float(row["max_direct_residual"])
        if not math.isfinite(r) or not 0<=r<=1e-8*(1+1e-12):
            raise ValueError("Physical residual policy violation")
        if not row["job_id"].isdigit():
            raise ValueError("Invalid job ID")
        key=(n,omp,alg,pol)
        ident=(key,rep)
        if ident in seen:
            raise ValueError("Duplicate case "+str(ident))
        seen.add(ident)
        groups[key].append(row)
    required={(n,omp,alg,pol)
              for n in (8,10,12) for omp in (1,2)
              for alg in (1,2) for pol in POLICIES}
    if set(groups)!=required: raise ValueError("Missing/extra benchmark configuration")
    for key,items in groups.items():
        if {int(x["repeat"]) for x in items}!={1,2,3}:
            raise ValueError("Missing benchmark repeat")
    for n in (8,10,12):
        if not nominal_cache_bytes(n,"full_tight") < nominal_cache_bytes(n,"half_generic") < nominal_cache_bytes(n,"full_generic"):
            raise ValueError("Analytical cache memory ordering contradicted")
    return groups


def summarize(groups):
    out=[]
    for (n,omp,alg,pol),items in sorted(groups.items()):
        timed=[x for x in items if x["phase"]=="measurement"]
        out.append(dict(n_sites=n,dim=6**(n//2),D=0.12,omp=omp,
                        alg=alg,policy=pol,nominal_cache_mib=nominal_cache_bytes(n,pol)/(1024**2),
                        repeats=len(timed),
                        wall_median_s=statistics.median(float(x["wall_seconds"]) for x in timed),
                        rss_median_mib=statistics.median(float(x["peak_rss_kib"]) for x in timed)/1024,
                        h_median=statistics.median(int(x["h_applications"]) for x in timed),
                        achieved_max_residual=max(float(x["max_direct_residual"]) for x in timed)))
    if len(out)!=36: raise ValueError("Incorrect summary size")
    return out


def plot(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=False)
    assets=[]
    for alg in (1,2):
        fig,ax=plt.subplots(figsize=(6.6,4.1),layout="constrained")
        for n in (8,10,12):
            items=sorted((x for x in rows if x["n_sites"]==n and x["alg"]==alg
                          and x["omp"]==1),
                         key=lambda x:x["nominal_cache_mib"])
            ax.plot([x["rss_median_mib"] for x in items],
                    [x["wall_median_s"] for x in items],
                    marker="o",label=f"L={n}")
            for x in items:
                lab={"half_generic":"half / default",
                     "full_tight":"full / exact MNTE",
                     "full_generic":"full / default"}[x["policy"]]
                ax.annotate(lab,(x["rss_median_mib"],x["wall_median_s"]),
                            xytext=(4,4),textcoords="offset points",fontsize=7)
        ax.set_yscale("log")
        ax.set(xlabel="Peak ED2 RSS (MiB)",
               ylabel="ED2 process wall time (s, shared CI)")
        ax.legend(title=f"ALG={alg}, OMP=1, NOE=1",fontsize=8)
        ax.grid(alpha=.25)
        for ext in (".svg",".png"):
            path=out/f"pareto_alg{alg}{ext}"
            fig.savefig(path,dpi=230)
            assets.append(path)
        plt.close(fig)
    return assets


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,default=DEFAULT)
    p.add_argument("--outdir",type=Path)
    p.add_argument("--validate-only",action="store_true")
    a=p.parse_args()
    rows=summarize(audit(a.csv))
    print("PARETO DATA VALID: 108/108 cases, 36 configurations, two measured repeats each",
          flush=True)
    if a.validate_only:
        if a.outdir is not None: p.error("Only choose one output mode")
        return 0
    if a.outdir is None: p.error("--outdir is required for figures")
    out=a.outdir.resolve()
    files=plot(rows,out)
    with (out/"cpc-pareto-summary.csv").open("w",newline="",encoding="utf-8") as fp:
        cols=list(rows[0])
        w=csv.DictWriter(fp,fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print("PARETO PLOTS:",len(files),"SVG/PNG assets",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
