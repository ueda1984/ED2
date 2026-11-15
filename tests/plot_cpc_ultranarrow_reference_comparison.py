#!/usr/bin/env python3
"""Audit frozen ED2/sparse/dense/second-order low-doublet precision data.

These CSVs are observations from independently accepted CI jobs, not
certified arbitrary-precision ground-state data. This audit checks all
physical configuration cells, source SHA/job provenance, internal
arithmetic of the reported differences, and consistency of the
reference values reused across two independent CI runs.

Plots show the signed difference from the *approximate second-order*
effective spin reference, NOT true unknown exact-solution error bars.
"""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parents[1]/"docs/development"
DENSE_DEFAULT=ROOT/"cpc-ultranarrow-dense-lapack-raw-ci.csv"
ED2_DEFAULT=ROOT/"cpc-ultranarrow-ed2-second-order-raw-ci.csv"

DENSE_FIELDS=("workflow_run_id","job_id","source_sha","n_sites",
              "epsilon","delta","driver","lowest_gap","edge_z",
              "edge_vs_second_order_abs_error",
              "edge_vs_sparse_eigsh_abs_error",
              "trace_correlations_vs_sparse_error")
ED2_FIELDS=("workflow_run_id","job_id","source_sha","delta",
            "OMP","NOK","NOM","residual_atol","ed2_edge_z",
            "sparse_eigsh_edge_z","second_order_edge_z",
            "ed2_vs_sparse_abs_error",
            "ed2_vs_second_order_abs_error","trace_corr_max_error")
DENSE_RUN="43707711743"
DENSE_JOB="132792288717"
DENSE_SHA="a1242ae0d4ea4cfc59d1560cd0d8baae33d0b78f"
ED2_RUN="43707339560"
ED2_SHA="6635e4b6e83eab0c2ccfd9933d27745a8f6a55fd"
ED2_JOBS={1:"132791152148",2:"132791151910"}
DELTAS=(-1e-7,0.,1e-7)


def records(path,fields,count):
    with path.open(newline="",encoding="utf-8") as stream:
        parser=csv.DictReader(stream)
        if parser.fieldnames!=list(fields):
            raise ValueError("Provenance data columns have changed: "+str(path))
        data=list(parser)
    if len(data)!=count or any(None in r for r in data):
        raise ValueError(f"Expected {count} complete records: {path}")
    return data


def finite(v):
    x=float(v)
    if not math.isfinite(x):
        raise ValueError("Nonfinite scientific result")
    return x


def approximately(actual,expected,atol=1e-14):
    if abs(actual-expected)>atol:
        raise ValueError(f"Frozen measurement arithmetic changed: {actual} versus {expected}")


def audit(dense_path,ed2_path):
    dr=records(dense_path,DENSE_FIELDS,9)
    er=records(ed2_path,ED2_FIELDS,24)
    dg=defaultdict(dict)
    eg=defaultdict(list)
    seen=set()
    for r in dr:
        if (r["workflow_run_id"],r["job_id"],r["source_sha"])!=(
                DENSE_RUN,DENSE_JOB,DENSE_SHA):
            raise ValueError("Unknown dense solver source/run/job")
        if r["n_sites"]!="8" or finite(r["epsilon"])!=5e-7:
            raise ValueError("Wrong dense physical Hamiltonian")
        d=finite(r["delta"]);driver=r["driver"]
        if d not in DELTAS or driver not in ("evr","evx","evd"):
            raise ValueError("Unexpected LAPACK driver/bias")
        if driver in dg[d]:
            raise ValueError("Duplicate dense driver")
        datum={k:finite(r[k]) for k in DENSE_FIELDS[7:]}
        if (not 1e-8<datum["lowest_gap"]<2e-7 or
            datum["trace_correlations_vs_sparse_error"]>1e-7 or
            datum["edge_vs_second_order_abs_error"]>7e-8):
            raise ValueError("Dense driver numerical acceptance violation")
        dg[d][driver]=datum
    for r in er:
        if r["workflow_run_id"]!=ED2_RUN or r["source_sha"]!=ED2_SHA:
            raise ValueError("Wrong ED2 source/run")
        d=finite(r["delta"]);omp=int(r["OMP"]);keep=int(r["NOK"])
        nom=int(r["NOM"]);atol=finite(r["residual_atol"])
        if d not in DELTAS or omp not in (1,2) or (keep,nom) not in ((2,16),(4,24)):
            raise ValueError("Unexpected ED2 parameter")
        if atol not in (1e-8,1e-12) or r["job_id"]!=ED2_JOBS[omp]:
            raise ValueError("Missing/altered ED2 job and residual gate")
        key=(d,omp,keep,nom,atol)
        if key in seen:
            raise ValueError("Duplicated ED2 trial")
        seen.add(key)
        item={k:finite(r[k]) for k in ED2_FIELDS[8:]}
        approximately(item["ed2_vs_sparse_abs_error"],
                      abs(item["ed2_edge_z"]-item["sparse_eigsh_edge_z"]))
        approximately(item["ed2_vs_second_order_abs_error"],
                      abs(item["ed2_edge_z"]-item["second_order_edge_z"]))
        if (item["ed2_vs_second_order_abs_error"]>1.1e-8 or
            item["trace_corr_max_error"]>3e-9):
            raise ValueError("Unexpected ED2 physical reference disagreement")
        eg[d].append(item)
    expected={(d,omp,k,m,t) for d in DELTAS for omp in (1,2)
              for k,m in ((2,16),(4,24)) for t in (1e-8,1e-12)}
    if seen!=expected or set(dg)!=set(DELTAS):
        raise ValueError("Provenance table missing workload cases")
    summary=[]
    for d in DELTAS:
        if set(dg[d])!={"evr","evx","evd"} or len(eg[d])!=8:
            raise ValueError("Missing complete solver outputs at one bias")
        sparse_vals={x["sparse_eigsh_edge_z"] for x in eg[d]}
        eff_vals={x["second_order_edge_z"] for x in eg[d]}
        if len(sparse_vals)!=1 or len(eff_vals)!=1:
            raise ValueError("Physical reference differs across requested ED2 workloads")
        sparse=next(iter(sparse_vals))
        eff=next(iter(eff_vals))
        for drv,row in dg[d].items():
            approximately(row["edge_vs_sparse_eigsh_abs_error"],abs(row["edge_z"]-sparse))
            approximately(row["edge_vs_second_order_abs_error"],abs(row["edge_z"]-eff))
        vals=[x["ed2_edge_z"] for x in eg[d]]
        summary.append(dict(
            delta=d,epsilon=5e-7,old_sparse_eigsh_z=sparse,
            dense_evr_z=dg[d]["evr"]["edge_z"],
            dense_evd_z=dg[d]["evd"]["edge_z"],
            second_order_z=eff,
            ed2_edge_z_median=statistics.median(vals),
            ed2_edge_z_min=min(vals),
            ed2_edge_z_max=max(vals),
            sparse_vs_second_abs=abs(sparse-eff),
            evr_vs_second_abs=dg[d]["evr"]["edge_vs_second_order_abs_error"],
            evd_vs_second_abs=dg[d]["evd"]["edge_vs_second_order_abs_error"],
            ed2_vs_second_max=max(x["ed2_vs_second_order_abs_error"] for x in eg[d]),
            ed2_vs_sparse_max=max(x["ed2_vs_sparse_abs_error"] for x in eg[d]),
            independent_dense_evr_gap=dg[d]["evr"]["lowest_gap"],
            independent_dense_evd_gap=dg[d]["evd"]["lowest_gap"],
        ))
    return summary


def charts(summary,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=False)
    assets=[]
    def save(fig,stem):
        for extension in (".svg",".png"):
            dest=out/(stem+extension)
            fig.savefig(dest,dpi=240)
            assets.append(dest)
        plt.close(fig)
    x=[r["delta"]/1e-7 for r in summary]
    fig,ax=plt.subplots(figsize=(6.2,3.9),layout="constrained")
    for k,label in (("old_sparse_eigsh_z","SciPy sparse eigsh"),
                    ("dense_evr_z","Dense LAPACK ZHEEVR"),
                    ("dense_evd_z","Dense LAPACK ZHEEVD"),
                    ("ed2_edge_z_median","ED2, 8 solver settings")):
        ax.plot(x,[r[k]-r["second_order_z"] for r in summary],
                marker="o",label=label)
    ax.axhline(0,ls="--",lw=1,label="Second-order Q resolvent")
    ax.set(xlabel="Edge bias delta / 1e-7",
           ylabel="Signed edge Sz difference from second-order model")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    save(fig,"ultra_edge_reference_offsets")
    fig,ax=plt.subplots(figsize=(6.2,3.9),layout="constrained")
    for k,label in (("sparse_vs_second_abs","SciPy sparse eigsh"),
                    ("evr_vs_second_abs","Dense ZHEEVR"),
                    ("evd_vs_second_abs","Dense ZHEEVD"),
                    ("ed2_vs_second_max","Max ED2, eight settings")):
        ax.plot(x,[r[k] for r in summary],marker="o",label=label)
    ax.set_yscale("log")
    ax.set(xlabel="Edge bias delta / 1e-7",
           ylabel="Absolute edge Sz difference from second-order model")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    save(fig,"ultra_solver_vs_second_order")
    return assets


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dense",type=Path,default=DENSE_DEFAULT)
    parser.add_argument("--ed2",type=Path,default=ED2_DEFAULT)
    parser.add_argument("--validate-only",action="store_true")
    parser.add_argument("--outdir",type=Path)
    a=parser.parse_args()
    summary=audit(a.dense,a.ed2)
    print("FROZEN CROSSREFERENCE AUDIT: 9 dense, 24 ED2 and 3 physical biases",
          "max ED2-to-effective",max(x["ed2_vs_second_max"] for x in summary),
          "max old-sparse-to-effective",max(x["sparse_vs_second_abs"] for x in summary),
          flush=True)
    if a.validate_only:
        if a.outdir is not None:
            parser.error("Select --validate-only OR --outdir")
        return 0
    if a.outdir is None:
        parser.error("--outdir is required to make figures")
    out=a.outdir.resolve()
    generated=charts(summary,out)
    with (out/"ultra_reference_method_summary.csv").open("w",newline="",encoding="utf-8") as fp:
        fields=list(summary[0])
        writer=csv.DictWriter(fp,fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)
    print("FROZEN CROSSREFERENCE FIGURES",len(generated),"SVG/PNG",flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
