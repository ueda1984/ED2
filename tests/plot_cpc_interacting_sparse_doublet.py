#!/usr/bin/env python3
"""Audit immutable 20-case interacting eight-site CI evidence and make figures.

The source job independently built the SciPy sparse Hamiltonian with
an interacting XY edge spin, and ED2 separately reconstructed both
native low eigenvectors to compare all Cartesian trace observables.

This script summarizes *accepted numerical tests*, not a bulk phase
diagram, convergence proof, universal solver certification, or CPU
benchmark. The low-subspace leakage bound requires external E2.
"""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path

CSV=Path(__file__).resolve().parents[1]/"docs/development/cpc-interacting-8site-raw-ci.csv"
HEADERS=("workflow_run_id","job_id","source_sha","n_sites","epsilon","delta",
         "omp","NOK","NOM","ground_gap","outer_gap","edge_linear_entropy",
         "trace_mag_max_error","trace_corr_max_error",
         "direct_residual_1","direct_residual_2")
RUN="43705920631"
SHA="539e4d47c98e1eaba8ad80f214f311c468f919ce"
SCENARIOS=((.02,-.002),(.02,0.),(.02,.002),(.0005,0.),(.0005,2e-5))
JOBS={1:"132786843599",2:"132786843765"}


def audit(path):
    with path.open(newline="",encoding="utf-8") as inp:
        reader=csv.DictReader(inp)
        if reader.fieldnames!=list(HEADERS):
            raise ValueError("Unexpected interacting-low-doublet CSV schema")
        raw=list(reader)
    if len(raw)!=20 or any(None in r for r in raw):
        raise ValueError("Expected twenty complete provenance-checked cases")
    seen=set()
    validated=[]
    controls=defaultdict(list)
    for r in raw:
        n=int(r["n_sites"]);eps=float(r["epsilon"]);delta=float(r["delta"])
        omp=int(r["omp"]);keep=int(r["NOK"]);nom=int(r["NOM"])
        if n!=8 or (eps,delta) not in SCENARIOS:
            raise ValueError("Unexpected physical Hamiltonian")
        if (keep,nom) not in ((2,16),(4,24)) or omp not in JOBS:
            raise ValueError("Unexpected solver workspace/OMP")
        if r["source_sha"]!=SHA or r["workflow_run_id"]!=RUN or r["job_id"]!=JOBS[omp]:
            raise ValueError("Original CI code or job provenance changed")
        key=(eps,delta,omp,keep,nom)
        if key in seen:
            raise ValueError("Duplicated model/workspace/thread setting")
        seen.add(key)
        vals={field:float(r[field]) for field in
              ("ground_gap","outer_gap","edge_linear_entropy",
               "trace_mag_max_error","trace_corr_max_error",
               "direct_residual_1","direct_residual_2")}
        if not all(math.isfinite(x) for x in vals.values()):
            raise ValueError("Nonfinite physical or residual observation")
        if not 1e-6<vals["ground_gap"]<.03 or not .3<vals["outer_gap"]<.4:
            raise ValueError("Low doublet isolation requirement changed")
        if not 0<vals["edge_linear_entropy"]<.1:
            raise ValueError("Edge entanglement not detected")
        for f in ("trace_mag_max_error","trace_corr_max_error"):
            if not 0<=vals[f]<=1e-5:
                raise ValueError("Failed doublet trace observable acceptance")
        if any(not 0<=vals[f]<=1e-8*(1+1e-12)
               for f in ("direct_residual_1","direct_residual_2")):
            raise ValueError("Numerical solver residual not accepted")
        # Conditional Frobenius leakage of an orthonormal two-state Ritz
        # subspace: sqrt(sum_j(r_j/(E2-lambda_j))^2).
        # Since E2-lambda_j >= E2-E1, a conservative upper estimate
        # uses the externally computed doublet-to-third-state gap for both.
        leakage_upper=math.hypot(vals["direct_residual_1"],
                                 vals["direct_residual_2"])/vals["outer_gap"]
        # Normalized trace expectation of any operator with norm <= 1
        # differs by at most sqrt(2)*leakage_upper, plus finite-precision
        # reference/eigenvector rounding. This is conditional on correct
        # third eigenvalue and near-orthonormal returned wavefunctions.
        trace_bound=math.sqrt(2.)*leakage_upper
        if max(vals["trace_mag_max_error"],vals["trace_corr_max_error"])>trace_bound+1e-9:
            raise ValueError("Trace error is incompatible with external-gap residual bound")
        entry=dict(epsilon=eps,delta=delta,omp=omp,keep=keep,capacity=nom,
                   **vals,conditional_leakage_upper=leakage_upper,
                   conditional_trace_bound_norm1=trace_bound)
        validated.append(entry)
        controls[(eps,delta)].append(entry)
    required={(e,d,o,k,m) for e,d in SCENARIOS
              for o in (1,2) for k,m in ((2,16),(4,24))}
    if seen!=required:
        raise ValueError("Not all accepted physics configurations present")
    for key,rs in controls.items():
        if len(rs)!=4:
            raise ValueError("Missing model controls")
        if max(x["ground_gap"] for x in rs)-min(x["ground_gap"] for x in rs)>1e-10:
            raise ValueError("Independent reference gap inconsistent across OMP")
        if max(x["edge_linear_entropy"] for x in rs)-min(x["edge_linear_entropy"] for x in rs)>1e-9:
            raise ValueError("Independent reference entropy inconsistent across OMP")
    return validated


def figures(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(exist_ok=False,parents=True)
    assets=[]
    def save(fig,stem):
        for ext in (".svg",".png"):
            path=out/(stem+ext)
            fig.savefig(path,dpi=230)
            assets.append(path)
        plt.close(fig)
    # Physically resolved gap: plot normalized longitudinal bias delta/eps.
    fig,ax=plt.subplots(figsize=(6.2,4.),layout="constrained")
    for eps in (.02,.0005):
        pts=sorted((r for r in rows if r["epsilon"]==eps and
                    r["omp"]==1 and r["keep"]==2),
                   key=lambda x:x["delta"])
        ax.plot([p["delta"]/eps for p in pts],
                [p["ground_gap"] for p in pts],
                marker="o",label=f"epsilon = {eps:g}")
    ax.set(yscale="log",xlabel="Dimensionless edge bias delta / epsilon",
           ylabel="Independent low-doublet gap")
    ax.legend()
    ax.grid(alpha=.25)
    save(fig,"interacting_sparse_gap")
    # Two-state trace accuracy relative to independently generated reference.
    fig,ax=plt.subplots(figsize=(6.2,4.),layout="constrained")
    for omp in (1,2):
        for keep in (2,4):
            pts=sorted((r for r in rows if r["omp"]==omp and r["keep"]==keep),
                       key=lambda x:x["ground_gap"])
            ax.plot([r["ground_gap"] for r in pts],
                    [max(r["trace_mag_max_error"],r["trace_corr_max_error"])
                     for r in pts],
                    marker="o",label=f"OMP={omp}, {keep}/{16 if keep==2 else 24}")
    ax.axhline(1e-8,linewidth=1,ls="--",label="Requested direct residual")
    ax.set(xscale="log",yscale="log",
           xlabel="Independent low doublet gap",
           ylabel="Largest doublet-averaged observable error")
    ax.legend(fontsize=8)
    ax.grid(alpha=.25)
    save(fig,"interacting_sparse_trace_error")
    # Independent entanglement diagnostic is not an observable/output
    # claim about ED2 wavefunctions; it confirms the reference is interacting.
    fig,ax=plt.subplots(figsize=(6.2,4.),layout="constrained")
    for eps in (.02,.0005):
        pts=sorted((r for r in rows if r["epsilon"]==eps and
                    r["omp"]==1 and r["keep"]==2),
                   key=lambda x:x["delta"])
        ax.plot([r["delta"]/eps for r in pts],
                [r["edge_linear_entropy"] for r in pts],
                marker="o",label=f"epsilon = {eps:g}")
    ax.set(yscale="log",xlabel="Dimensionless edge bias delta / epsilon",
           ylabel="Independent edge linear entropy")
    ax.legend()
    ax.grid(alpha=.25)
    save(fig,"interacting_sparse_edge_entropy")
    return assets


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,default=CSV)
    p.add_argument("--outdir",type=Path)
    p.add_argument("--validate-only",action="store_true")
    args=p.parse_args()
    rows=audit(args.csv)
    maxerror=max(max(r["trace_mag_max_error"],r["trace_corr_max_error"])
                 for r in rows)
    print("INTERACTING SPARSE FROZEN EVIDENCE 20/20 verified; max trace error",
          maxerror,flush=True)
    if args.validate_only:
        if args.outdir: p.error("--validate-only and --outdir are exclusive")
        return 0
    if not args.outdir: p.error("--outdir required to render figures")
    assets=figures(rows,args.outdir.resolve())
    with (args.outdir/"interacting-sparse-summary.csv").open(
            "w",newline="",encoding="utf-8") as fp:
        columns=list(rows[0])
        writer=csv.DictWriter(fp,fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print("INTERACTING SPARSE FIGURES",len(assets),"SVG/PNG created",flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
