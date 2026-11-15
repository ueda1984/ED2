#!/usr/bin/env python3
"""Validate frozen 12-case eight-site doublet evidence and render figures.

The first/second ED2 eigenvectors within a degenerate doublet have
basis-dependent magnetic moments; their mean is gauge invariant and
matches independently constructed SciPy low-doublet projector traces.
This is a deliberately decoupled edge spin stress test, not a magnetic
phase diagram or direct evidence for the stationary connected chain.

The committed CSV came from successful source-archive CI run 43705159096.
"""
import argparse
import csv
from pathlib import Path
import math

FIELDS=("workflow_run_id","job_id","source_sha","n_sites","delta",
        "NOK","NOM","OMP","gap","outside_gap","edge_z_first",
        "edge_z_second","doublet_mean_z","trace_mag_max_error",
        "trace_corr_max_error")
RUN="43705159096"
SHA="eae33e95cbeb0c25937af4ebc2a8a24c2b3fca58"
DEFAULT=Path(__file__).resolve().parents[1]/"docs/development/cpc-sparse-doublet-8site-raw-ci.csv"
DELTAS=(0.,1e-10,1e-6)
WKS=((2,16),(4,24))


def load_and_audit(path):
    with path.open(encoding="utf-8",newline="") as f:
        r=csv.DictReader(f)
        if r.fieldnames!=list(FIELDS):
            raise ValueError("Unknown source record layout")
        raw=list(r)
    if len(raw)!=12 or any(None in x for x in raw):
        raise ValueError("Expected 12 complete cases")
    rows=[]
    seen=set()
    for item in raw:
        if item["workflow_run_id"]!=RUN or item["source_sha"]!=SHA:
            raise ValueError("Unexpected job/source provenance")
        n=int(item["n_sites"]);d=float(item["delta"])
        k=int(item["NOK"]);m=int(item["NOM"]);omp=int(item["OMP"])
        if n!=8 or d not in DELTAS or (k,m) not in WKS or omp not in (1,2):
            raise ValueError("Unrecognized stress workload")
        expected_job={1:"132784555495",2:"132784555339"}[omp]
        if item["job_id"]!=expected_job:
            raise ValueError("Unexpected source job identifier")
        key=(d,k,m,omp)
        if key in seen:
            raise ValueError("Duplicated test case")
        seen.add(key)
        vals={v:float(item[v]) for v in ("gap","outside_gap","edge_z_first",
               "edge_z_second","doublet_mean_z","trace_mag_max_error",
               "trace_corr_max_error")}
        if any(not math.isfinite(x) for x in vals.values()):
            raise ValueError("Non-finite observable")
        if not math.isclose(vals["gap"],d,abs_tol=1e-14):
            raise ValueError("Wrong low gap for detached edge")
        if not .3<vals["outside_gap"]<.4:
            raise ValueError("Incorrect external gap")
        if any(abs(vals[x])>0.50000001 for x in ("edge_z_first","edge_z_second")):
            raise ValueError("Spin-1/2 physical bound violated")
        if abs(vals["doublet_mean_z"]-(vals["edge_z_first"]+vals["edge_z_second"])/2)>1e-15:
            raise ValueError("Incorrect doublet mean")
        if abs(vals["doublet_mean_z"])>5e-10 or max(
                vals["trace_mag_max_error"],vals["trace_corr_max_error"])>3.7e-10:
            raise ValueError("Trace comparison contradicts accepted CI")
        rows.append(dict(delta=d,omp=omp,keep=k,capacity=m,**vals))
    required={(d,k,m,o) for d in DELTAS for k,m in WKS for o in (1,2)}
    if set(seen)!=required:
        raise ValueError("Missing verification setting")
    return rows


def draw(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=False)
    labels={d:i for i,d in enumerate(DELTAS)}
    pictures=[]
    fig,ax=plt.subplots(figsize=(6.5,3.9),layout="constrained")
    for omp in (1,2):
        for k,m in WKS:
            vals=sorted((r for r in rows if r["omp"]==omp and
                         r["keep"]==k and r["capacity"]==m),
                        key=lambda r:r["delta"])
            # Slight horizontal offset only to show separate numerical samples.
            offset={1:-.07,2:.07}[omp]+({2:-.018,4:.018}[k])
            x=[labels[r["delta"]]+offset for r in vals]
            ax.plot(x,[r["edge_z_first"] for r in vals],marker="o",
                    label=f"state 1, OMP={omp}, {k}/{m}")
            # The second vector is the complementary choice in the doublet:
            ax.scatter(x,[r["edge_z_second"] for r in vals],
                       marker="x",alpha=.7)
    ax.axhline(0.,ls="--",lw=1,label="doublet mean = 0")
    ax.set_xticks(list(range(3)),["0","1e-10","1e-6"])
    ax.set(xlabel="Detached-edge splitting parameter delta",
           ylabel="ED2 edge spin <S_z> (first and second states)")
    ax.grid(alpha=.25)
    ax.legend(fontsize=7,ncol=2)
    for extension in (".svg",".png"):
        p=out/("doublet_individual_vs_trace"+extension)
        fig.savefig(p,dpi=240)
        pictures.append(p)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.5,3.9),layout="constrained")
    for omp in (1,2):
        for k,m in WKS:
            vals=sorted((r for r in rows if r["omp"]==omp and
                         r["keep"]==k and r["capacity"]==m),
                        key=lambda r:r["delta"])
            x=[labels[r["delta"]] for r in vals]
            y=[max(r["trace_mag_max_error"],r["trace_corr_max_error"]) for r in vals]
            ax.plot(x,y,marker="o",label=f"OMP={omp}, {k}/{m}")
    ax.set_yscale("log")
    ax.axhline(3.7e-10,ls="--",lw=1,label="Validated maximum envelope")
    ax.set_xticks(list(range(3)),["0","1e-10","1e-6"])
    ax.set(xlabel="Detached-edge splitting parameter delta",
           ylabel="Largest doublet trace observable error")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8,ncol=2)
    for extension in (".svg",".png"):
        p=out/("doublet_trace_reference_error"+extension)
        fig.savefig(p,dpi=240)
        pictures.append(p)
    plt.close(fig)
    return pictures


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,default=DEFAULT)
    p.add_argument("--output",type=Path)
    p.add_argument("--validate-only",action="store_true")
    a=p.parse_args()
    rows=load_and_audit(a.csv)
    print("LOW-DOUBLET FROZEN EVIDENCE: 12/12 valid; max trace error",
          max(max(r["trace_mag_max_error"],r["trace_corr_max_error"])
              for r in rows),flush=True)
    if a.validate_only:
        if a.output: p.error("--validate-only cannot be combined with --output")
        return 0
    if a.output is None: p.error("--output required for images")
    files=draw(rows,a.output.resolve())
    print("LOW-DOUBLET FIGURES:",len(files),"SVG/PNG generated",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
