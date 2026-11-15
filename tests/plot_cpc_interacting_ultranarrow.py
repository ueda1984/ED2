#!/usr/bin/env python3
"""Strict frozen-CI audit and exploratory plots for 24 ultra-small-gap cases.

This numerical stress set has an interacting, non-U(1) eight-site
Hamiltonian with an inner gap around 1e-7. Individual ground-edge
moments are compared with an independent sparse reference; the
two-state trace is a more robust projector-invariant diagnostic.

The reported 1e-7-scale individual discrepancies are not certified
error bars for the ED2 eigenvector alone: reference eigenvectors at
ultra-small gaps also need extra high-precision checks.
"""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path

DEFAULT=Path(__file__).resolve().parents[1]/"docs/development/cpc-interacting-ultranarrow-raw-ci.csv"
HEADERS=("workflow_run_id","job_id","source_sha","n_sites","epsilon","delta",
         "OMP","NOK","NOM","residual_atol","solver_status","classification",
         "lowest_gap","external_gap","reference_edge_linear_entropy",
         "individual_edge_z_error","trace_correlation_max_error","H_applications")
SHA="cf0b73af4dbc84d902efd6f33526f0a5a8c423a5"
RUN="43706362016"
JOBS={1:"132788202415",2:"132788202597"}
DELTAS=(-1e-7,0.,1e-7)


def audit(path):
    with path.open(newline="",encoding="utf-8") as f:
        parser=csv.DictReader(f)
        if parser.fieldnames!=list(HEADERS):
            raise ValueError("Unrecognized ultra-small-gap CI data")
        raw=list(parser)
    if len(raw)!=24 or any(None in r for r in raw):
        raise ValueError("Expected exactly 24 full records")
    seen=set();out=[];byref=defaultdict(list)
    for r in raw:
        if r["workflow_run_id"]!=RUN or r["source_sha"]!=SHA:
            raise ValueError("Unrecognized source provenance")
        n=int(r["n_sites"]);e=float(r["epsilon"]);d=float(r["delta"])
        omp=int(r["OMP"]);keep=int(r["NOK"]);nom=int(r["NOM"])
        atol=float(r["residual_atol"]);status=int(r["solver_status"])
        if n!=8 or e!=5e-7 or d not in DELTAS or omp not in JOBS:
            raise ValueError("Unexpected Hamiltonian/OMP")
        if r["job_id"]!=JOBS[omp] or (keep,nom) not in ((2,16),(4,24)):
            raise ValueError("Wrong GitHub job or restart configuration")
        if atol not in (1e-8,1e-12) or status!=0:
            raise ValueError("Trial did not pass independently checked residuals")
        if r["classification"]!="accepted_correct_individual_ground":
            raise ValueError("Unexpected physical classification")
        key=(d,omp,keep,nom,atol)
        if key in seen: raise ValueError("Repeated case")
        seen.add(key)
        nums={field:float(r[field]) for field in
              ("lowest_gap","external_gap","reference_edge_linear_entropy",
               "individual_edge_z_error","trace_correlation_max_error",
               "H_applications")}
        if any(not math.isfinite(v) for v in nums.values()):
            raise ValueError("Non-finite performance or physics quantity")
        if not 1e-8<nums["lowest_gap"]<2e-7 or not .3<nums["external_gap"]<.4:
            raise ValueError("Expected separated low doublet absent")
        if not 0<=nums["individual_edge_z_error"]<1e-5:
            raise ValueError("Unaccepted individual magnetic observable")
        if not 0<=nums["trace_correlation_max_error"]<1e-5:
            raise ValueError("Unaccepted subspace trace observable")
        if nums["H_applications"]<=0: raise ValueError("Invalid H action count")
        byref[d].append(nums)
        out.append(dict(delta=d,omp=omp,keep=keep,capacity=nom,
                        residual_atol=atol,**nums))
    expected={(d,omp,k,m,atol) for d in DELTAS
              for omp in (1,2) for k,m in ((2,16),(4,24))
              for atol in (1e-8,1e-12)}
    if seen!=expected:
        raise ValueError("Incomplete accepted run matrix")
    for delta,rs in byref.items():
        if len(rs)!=8: raise ValueError("Missing independent reference cases")
        for name,tol in (("lowest_gap",5e-13),("external_gap",5e-13)):
            if max(x[name] for x in rs)-min(x[name] for x in rs)>tol:
                raise ValueError("Independently computed spectrum changed across solvers")
    return out


def plots(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(exist_ok=False,parents=True)
    assets=[]
    def save(fig,name):
        for extension in (".svg",".png"):
            p=out/(name+extension)
            fig.savefig(p,dpi=230)
            assets.append(p)
        plt.close(fig)
    for column,name,yaxis in (
        ("individual_edge_z_error","ultra_individual_edge_error",
         "Absolute single-ground edge Sz disagreement"),
        ("trace_correlation_max_error","ultra_doublet_trace_error",
         "Maximum doublet-averaged correlation disagreement")):
        fig,ax=plt.subplots(figsize=(6.4,4.),layout="constrained")
        for omp in (1,2):
            for keep in (2,4):
                for atol in (1e-8,1e-12):
                    points=sorted((x for x in rows if
                        x["omp"]==omp and x["keep"]==keep and
                        x["residual_atol"]==atol),key=lambda x:x["delta"])
                    if not points: raise ValueError("Empty plot group")
                    style="-" if atol==1e-8 else "--"
                    ax.plot([x["delta"]/5e-7 for x in points],
                            [x[column] for x in points],ls=style,marker="o",
                            label=f"OMP{omp} {keep}/{16 if keep==2 else 24} req {atol:.0e}")
        ax.set_yscale("log")
        ax.set(xlabel="Edge bias delta / epsilon",ylabel=yaxis)
        ax.grid(alpha=.2)
        ax.legend(fontsize=7,ncol=2)
        save(fig,name)
    fig,ax=plt.subplots(figsize=(6.2,3.8),layout="constrained")
    for omp in (1,2):
        for keep in (2,4):
            mean=[]
            for atol in (1e-8,1e-12):
                vals=[x["H_applications"] for x in rows
                      if x["omp"]==omp and x["keep"]==keep
                      and x["residual_atol"]==atol]
                mean.append(sum(vals)/len(vals))
            ax.plot([1e-8,1e-12],mean,marker="o",
                    label=f"OMP{omp}, {keep}/{16 if keep==2 else 24}")
    ax.invert_xaxis()
    ax.set(xscale="log",xlabel="Requested absolute direct residual",
           ylabel="Mean H applications over three biases")
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
    save(fig,"ultra_hamiltonian_action_cost")
    return assets


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,default=DEFAULT)
    p.add_argument("--outdir",type=Path)
    p.add_argument("--validate-only",action="store_true")
    a=p.parse_args()
    records=audit(a.csv)
    print("ULTRA INTERACTING FROZEN 24/24",
          "max individual error",max(x["individual_edge_z_error"] for x in records),
          "max trace correlation error",
          max(x["trace_correlation_max_error"] for x in records),flush=True)
    if a.validate_only:
        if a.outdir: p.error("Use either --validate-only or --outdir")
        return 0
    if not a.outdir: p.error("--outdir required")
    paths=plots(records,a.outdir.resolve())
    print("ULTRA INTERACTING FIGURES",len(paths),"SVG/PNG",flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
