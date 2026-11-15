#!/usr/bin/env python3
"""Audit archived full-cache MNTE measurements and plot CPC resource figures.

The source data are the 72 independently accepted ED2 CI trials from
workflow 43703623648, with fixed two-site periodic non-U(1) mixed spins.
The "optimal" row width means a measured exact-enumeration maximum
over representative states of *this* Hamiltonian and basis. It is not
a universally safe setting for arbitrary couplings or truncations.

The expected memory reduction is computed from:
  wk_loc default INTEGER(4): 4 bytes per entry
  wk_ele COMPLEX(8): 16 bytes per entry
  cache shape: (MNTE+1, full_dim)
Thus B_saved = 20 * (default_MNTE - measured_MNTE) * full_dim bytes.
Total RSS includes additional overhead and is not predicted by this formula.

Shared CI timings and RSS must not be described as controlled hardware
speedups. Warmups are excluded from all medians.
"""
import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path
import statistics

SOURCE_CSV=Path(__file__).resolve().parents[1]/"docs/development/cpc-mnte-capacity-raw-ci.csv"
FIELDS=("workflow_run_id","job_id","source_sha","n_sites","D","omp",
        "alg","phase","repeat","capacity_policy","reported_mnte",
        "observed_optimal_mnte","solver_status","h_applications",
        "wall_seconds","peak_rss_kib")
RUN_ID="43703623648"
SOURCE_SHA="315bed65b0c21d4db4f3af834d5807d361839d29"
OPTIMUM={8:27,10:34,12:41}


def median(v):
    return statistics.median(v)


def audit(path):
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=list(FIELDS):
            raise ValueError("Unrecognized MNTE CSV columns "+str(reader.fieldnames))
        rows=list(reader)
    if len(rows)!=72 or any(None in r for r in rows):
        raise ValueError("Expected 72 complete MNTE CI records")
    grouped=defaultdict(list)
    seen=set()
    for row in rows:
        if row["workflow_run_id"]!=RUN_ID or row["source_sha"]!=SOURCE_SHA:
            raise ValueError("MNTE provenance has changed")
        n=int(row["n_sites"]); omp=int(row["omp"]); alg=int(row["alg"])
        rep=int(row["repeat"]); phase=row["phase"]
        pol=row["capacity_policy"]; D=float(row["D"])
        mnte=int(row["reported_mnte"]); opt=int(row["observed_optimal_mnte"])
        if n not in (8,10,12) or omp not in (1,2) or alg not in (1,2):
            raise ValueError("Unexpected size/solver")
        if pol not in ("default","calibrated") or rep not in (1,2,3):
            raise ValueError("Unexpected work-array policy/repetition")
        if phase!=("warmup" if rep==1 else "measurement"):
            raise ValueError("Warm-up designation changed")
        if abs(D-.12)>1e-12 or row["solver_status"]!="0":
            raise ValueError("MNTE workload/solver-status changed")
        if mnte!=(10*n-8 if pol=="default" else OPTIMUM[n]) or opt!=OPTIMUM[n]:
            raise ValueError("Capacity calibration does not match frozen ED2 logs")
        if not row["job_id"].isdigit():
            raise ValueError("Invalid job provenance")
        for field in ("wall_seconds","peak_rss_kib","h_applications"):
            v=float(row[field])
            if v<=0 or not math.isfinite(v):
                raise ValueError("Invalid resource value "+field)
        key=(n,omp,alg,pol)
        identity=(key,rep)
        if identity in seen: raise ValueError("Duplicated CI trial "+str(identity))
        seen.add(identity)
        grouped[key].append(row)
    required={(n,omp,alg,pol) for n in (8,10,12)
              for omp in (1,2) for alg in (1,2)
              for pol in ("default","calibrated")}
    if set(grouped)!=required:
        raise ValueError("Missing benchmark configurations")
    for key,rows in grouped.items():
        if {int(r["repeat"]) for r in rows}!={1,2,3}:
            raise ValueError("Incomplete repetitions for "+str(key))
    for n in (8,10,12):
        for omp in (1,2):
            for alg in (1,2):
                num_hcalls={
                    int(r["h_applications"])
                    for policy in ("default","calibrated")
                    for r in grouped[(n,omp,alg,policy)]
                }
                if len(num_hcalls)!=1:
                    raise ValueError("MNTE setting changed Krylov action count: "+
                                     str((n,omp,alg,num_hcalls)))
    return grouped


def summary(grouped):
    results=[]
    for (n,omp,alg,policy),items in sorted(grouped.items()):
        timed=[r for r in items if r["phase"]=="measurement"]
        D=6**(n//2)
        mnte=(10*n-8) if policy=="default" else OPTIMUM[n]
        results.append(dict(n_sites=n,full_dimension=D,omp=omp,alg=alg,
                            capacity_policy=policy,MNTE=mnte,
                            optimum_MNTE=OPTIMUM[n],
                            cached_bytes=20*(mnte+1)*D,
                            measured_repeats=len(timed),
                            wall_median_s=median([float(x["wall_seconds"]) for x in timed]),
                            rss_median_kib=median([float(x["peak_rss_kib"]) for x in timed]),
                            h_applications=median([int(x["h_applications"]) for x in timed])))
    for n in (8,10,12):
        for omp in (1,2):
            for alg in (1,2):
                base=next(x for x in results if (x["n_sites"],x["omp"],x["alg"],x["capacity_policy"])==
                          (n,omp,alg,"default"))
                tuned=next(x for x in results if (x["n_sites"],x["omp"],x["alg"],x["capacity_policy"])==
                           (n,omp,alg,"calibrated"))
                if base["h_applications"]!=tuned["h_applications"]:
                    raise ValueError("Different H-call count in paired capacity comparison")
                measured=(base["rss_median_kib"]-tuned["rss_median_kib"])
                expected=(base["cached_bytes"]-tuned["cached_bytes"])/1024
                if measured<=0: raise ValueError("Optimized capacity increased process RSS")
                if abs(measured-expected)>max(1024,0.1*expected):
                    raise ValueError("Memory allocation and RSS difference inconsistent "+
                                     str((n,omp,alg,measured,expected)))
    return results


def draw(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=False)
    assets=[]
    def save(fig,base):
        for ext in (".svg",".png"):
            name=out/(base+ext)
            fig.savefig(name,dpi=230)
            assets.append(name)
        plt.close(fig)
    # Show both RAM consumption and capacity effect on the same physical task.
    fig,ax=plt.subplots(figsize=(6.0,3.9),layout="constrained")
    for policy in ("default","calibrated"):
        values=[next(x for x in rows if x["n_sites"]==n and x["omp"]==1
                     and x["alg"]==1 and x["capacity_policy"]==policy)
                for n in (8,10,12)]
        ax.plot([x["n_sites"] for x in values],
                [x["rss_median_kib"]/1024 for x in values],marker="o",
                label=f"MNTE {policy}")
    ax.set(xlabel="Mixed-spin chain length L",ylabel="Peak ED2 RSS (MiB)")
    ax.set_xticks([8,10,12])
    ax.grid(alpha=.25)
    ax.legend(title="ALG=1, OMP=1")
    save(fig,"mnte_peak_rss")
    # Compare theoretical work-array delta to the independently observed RSS delta.
    fig,ax=plt.subplots(figsize=(6.0,3.9),layout="constrained")
    theory=[]; observed=[]
    for n in (8,10,12):
        base=next(x for x in rows if x["n_sites"]==n and x["omp"]==1
                  and x["alg"]==1 and x["capacity_policy"]=="default")
        tuned=next(x for x in rows if x["n_sites"]==n and x["omp"]==1
                   and x["alg"]==1 and x["capacity_policy"]=="calibrated")
        theory.append((base["cached_bytes"]-tuned["cached_bytes"])/(1024*1024))
        observed.append((base["rss_median_kib"]-tuned["rss_median_kib"])/1024)
    ax.plot([8,10,12],theory,marker="o",label="Allocated cache reduction")
    ax.plot([8,10,12],observed,marker="s",ls="--",label="Observed RSS reduction")
    ax.set(xlabel="Mixed-spin chain length L",ylabel="Saved memory (MiB)")
    ax.set_xticks([8,10,12])
    ax.grid(alpha=.25)
    ax.legend(title="ALG=1, OMP=1")
    save(fig,"mnte_memory_prediction")
    fig,ax=plt.subplots(figsize=(6.0,3.9),layout="constrained")
    for alg in (1,2):
        ratios=[]
        for n in (8,10,12):
            base=next(x for x in rows if x["n_sites"]==n and x["omp"]==1
                      and x["alg"]==alg and x["capacity_policy"]=="default")
            tuned=next(x for x in rows if x["n_sites"]==n and x["omp"]==1
                       and x["alg"]==alg and x["capacity_policy"]=="calibrated")
            ratios.append(base["wall_median_s"]/tuned["wall_median_s"])
        ax.plot([8,10,12],ratios,marker="o",label=f"ALG={alg}")
    ax.axhline(1.0,ls="--",lw=1)
    ax.set(xlabel="Mixed-spin chain length L",
           ylabel="Default/optimized wall-time ratio on shared CI")
    ax.set_xticks([8,10,12])
    ax.grid(alpha=.25)
    ax.legend(title="OMP=1, one eigenpair")
    save(fig,"mnte_time_ratio")
    return assets


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,default=SOURCE_CSV)
    p.add_argument("--outdir",type=Path)
    p.add_argument("--validate-only",action="store_true")
    a=p.parse_args()
    rows=summary(audit(a.csv))
    if len(rows)!=24: raise ValueError("Expected 24 aggregate cells")
    print("MNTE DATA AUDIT: 72/72 rows, 24 grouped workload cells,",
          "paired memory-allocation check passed",flush=True)
    if a.validate_only:
        if a.outdir is not None: p.error("Use either --validate-only or --outdir")
        return 0
    if a.outdir is None: p.error("--outdir required for plotting")
    path=a.outdir.resolve()
    assets=draw(rows,path)
    with (path/"cpc-mnte-summary.csv").open("w",newline="",encoding="utf-8") as fp:
        fields=list(rows[0])
        writer=csv.DictWriter(fp,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print("MNTE PLOTS:",len(assets),"SVG/PNG files generated",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
