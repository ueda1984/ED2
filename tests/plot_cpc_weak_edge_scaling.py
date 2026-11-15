#!/usr/bin/env python3
"""Audit immutable 4-point LAPACK scaling + 12 native ED2 second-order cases.

These files preserve all 3 job IDs and one source SHA for CI #43710019532.
They are numerical test observations on a deliberately modified eight-site
non-U(1) mixed-spin chain, not independent new physics or universal
performance. Warmups/timing are NOT part of this figure.

Usage:
  python3 tests/plot_cpc_weak_edge_scaling.py --validate-only
  python3 tests/plot_cpc_weak_edge_scaling.py --outdir /new/empty/location
"""
import argparse
import csv
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/"docs/development"
SCALING=ROOT/"cpc-weak-edge-scaling-new-ci.csv"
REFERENCE=ROOT/"cpc-second-order-reference-new-ci.csv"
ED2=ROOT/"cpc-second-order-ed2-new-ci.csv"
RUN="43710019532"
SHA="eddc0dafe3f33fca1591460469c8531e508d35e0"
JOBS={1:"132799345686",2:"132799345656"}
SCALES=(5e-7,2e-6,5e-6,1e-5)
DELTAS=(-1e-7,0.,1e-7)


def load(path):
    with path.open(newline="",encoding="utf-8") as fp:
        reader=csv.DictReader(fp)
        if not reader.fieldnames or len(set(reader.fieldnames))!=len(reader.fieldnames):
            raise ValueError("Invalid input header: "+str(path))
        rows=list(reader)
    if not rows or any(None in row for row in rows):
        raise ValueError("Empty/incomplete CSV: "+str(path))
    return rows


def numeric(row,fields):
    out={}
    for k in fields:
        v=float(row[k])
        if not math.isfinite(v): raise ValueError("Non-finite "+k)
        out[k]=v
    return out


def audit(scalingpath,referencepath,ed2path):
    scaling=load(scalingpath)
    reference=load(referencepath)
    ed2=load(ed2path)
    if (len(scaling),len(reference),len(ed2))!=(4,6,12):
        raise ValueError("Expected exactly 4/6/12 frozen cases")
    for collection in (scaling,reference,ed2):
        for row in collection:
            if row["run"]!=RUN or row["sha"]!=SHA:
                raise ValueError("Unrecognized source provenance")
    sx=set();rs=[];coefficients=[]
    for row in scaling:
        if row["job"]!="132799345373":
            raise ValueError("Unexpected independent dense source job")
        v=numeric(row,("epsilon","dense_Sz","Sz_over_epsilon",
                         "derived_coefficient","dense_gap",
                         "second_order_gap_abs_error","second_order_Sz_abs_error"))
        e=v["epsilon"]
        if e not in SCALES or e in sx: raise ValueError("Missing/duplicate epsilon")
        sx.add(e)
        if not 0<v["dense_Sz"]<1e-3 or v["second_order_gap_abs_error"]>1e-12:
            raise ValueError("Failed LAPACK scaling acceptance")
        if abs(v["dense_Sz"]/e-v["Sz_over_epsilon"])>1e-8:
            raise ValueError("Magnetization slope arithmetic inconsistent")
        if abs(v["Sz_over_epsilon"]-v["derived_coefficient"])>.05:
            raise ValueError("Weak coupling coefficient disagreement")
        coefficients.append(v["derived_coefficient"])
        rs.append(v)
    if sx!=set(SCALES) or max(coefficients)-min(coefficients)>1e-9:
        raise ValueError("Source-dependent analytical coefficient")
    if not 1.4<coefficients[0]<1.6: raise ValueError("Coefficient outside expected reference")
    ref={}
    for row in reference:
        o=int(row["OMP"])
        if o not in JOBS or row["job"]!=JOBS[o]:
            raise ValueError("Wrong OMP/reference job")
        v=numeric(row,("eps","delta","gap_1st","gap_2nd","gap_sparse",
                       "err_gap_1st","err_gap_2nd","Sz_1st","Sz_2nd",
                       "Sz_dressed","Sz_sparse"))
        if v["eps"]!=5e-7 or v["delta"] not in DELTAS:
            raise ValueError("Incorrect sparse reference Hamiltonian")
        key=(o,v["delta"])
        if key in ref:raise ValueError("Duplicated independent reference")
        if abs(abs(v["gap_2nd"]-v["gap_sparse"])-v["err_gap_2nd"])>1e-18:
            raise ValueError("Reference gap error arithmetic inconsistent")
        if abs(v["Sz_1st"]+v["delta"]/(2*v["gap_1st"]))>1e-9:
            raise ValueError("First-order two-level edge magnetization inconsistent")
        if v["err_gap_2nd"]>1e-10:
            raise ValueError("Second-order reference gap rejected")
        ref[key]=v
    if set(ref)!={(o,d) for o in JOBS for d in DELTAS}:
        raise ValueError("Incomplete reference matrix")
    seen=set();comparisons=[]
    for row in ed2:
        o=int(row["OMP"]);k=int(row["NOK"]);m=int(row["NOM"])
        if o not in JOBS or row["job"]!=JOBS[o] or (k,m) not in ((2,16),(4,24)):
            raise ValueError("Unexpected ED2 job or solver workspace")
        v=numeric(row,("eps","delta","ED2_Sz","err_sparse","err_2nd","trace_error"))
        if v["eps"]!=5e-7 or v["delta"] not in DELTAS:
            raise ValueError("Unexpected ED2 target Hamiltonian")
        key=(o,v["delta"],k,m)
        if key in seen: raise ValueError("Duplicate ED2 trial")
        seen.add(key)
        q=ref[(o,v["delta"])]
        if abs(abs(v["ED2_Sz"]-q["Sz_2nd"])-v["err_2nd"])>1e-14:
            raise ValueError("ED2/second-order discrepancy inconsistent")
        if abs(abs(v["ED2_Sz"]-q["Sz_sparse"])-v["err_sparse"])>1e-14:
            raise ValueError("ED2/sparse difference inconsistent")
        if not 0<=v["trace_error"]<1e-5 or v["err_2nd"]>1e-7:
            raise ValueError("Unaccepted ED2 individual or trace discrepancy")
        comparisons.append(dict(OMP=o,NOK=k,NOM=m,**v))
    required={(o,d,k,m) for o in JOBS for d in DELTAS for k,m in ((2,16),(4,24))}
    if seen!=required:
        raise ValueError("Missing ED2 configurations")
    return sorted(rs,key=lambda x:x["epsilon"]),ref,comparisons


def draw(scales,refs,comp,outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir.mkdir(parents=True,exist_ok=False)
    files=[]
    def save(fig,base):
        for ext in (".svg",".png"):
            p=outdir/(base+ext)
            fig.savefig(p,dpi=230)
            files.append(p)
        plt.close(fig)
    eps=[v["epsilon"] for v in scales]
    obs=[v["dense_Sz"] for v in scales]
    coeff=scales[0]["derived_coefficient"]
    fig,ax=plt.subplots(figsize=(6.4,4.),layout="constrained")
    ax.plot(eps,obs,"o-",label="Full 8-site LAPACK ZHEEVR")
    ax.plot(eps,[coeff*e for e in eps],"--",label="Core Q-resolvent asymptote C|epsilon|")
    ax.set(xscale="log",yscale="log",
           xlabel="Edge coupling epsilon (>0)",
           ylabel="Ground-edge <Sz> at zero edge bias")
    ax.grid(True,alpha=.25)
    ax.legend(fontsize=8)
    save(fig,"weak_edge_zero_bias_scaling")
    fig,ax=plt.subplots(figsize=(6.4,4.),layout="constrained")
    ax.plot(eps,[v["Sz_over_epsilon"] for v in scales],"o-",
            label="Dense LAPACK response / epsilon")
    ax.axhline(coeff,ls="--",label="Second-order Q-resolvent coefficient")
    ax.set(xscale="log",xlabel="Edge coupling epsilon (>0)",
           ylabel="Ground-edge polarization per unit epsilon")
    ax.grid(True,alpha=.25)
    ax.legend(fontsize=8)
    save(fig,"weak_edge_induced_coefficient")
    fig,ax=plt.subplots(figsize=(6.4,4.),layout="constrained")
    for o in (1,2):
        for k in (2,4):
            pts=sorted((x for x in comp if x["OMP"]==o and x["NOK"]==k),
                       key=lambda x:x["delta"])
            ax.plot([p["delta"] for p in pts],[p["err_2nd"] for p in pts],
                    marker="o",label=f"ED2 OMP{o} NOK{k} vs 2nd order")
    ax.set(yscale="log",xlabel="Edge bias delta (epsilon=5e-7)",
           ylabel="Absolute edge Sz disagreement")
    ax.grid(True,alpha=.25)
    ax.legend(fontsize=8)
    save(fig,"weak_edge_ED2_second_order_disagreement")
    return files


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scaling",type=Path,default=SCALING)
    ap.add_argument("--reference",type=Path,default=REFERENCE)
    ap.add_argument("--ed2",type=Path,default=ED2)
    ap.add_argument("--validate-only",action="store_true")
    ap.add_argument("--outdir",type=Path)
    args=ap.parse_args()
    scales,ref,comp=audit(args.scaling,args.reference,args.ed2)
    print("FROZEN SECOND-ORDER AUDIT: 4/4 dense cases, 6/6 Q-reference",
          "cases, 12/12 ED2-native cases",flush=True)
    if args.validate_only:
        if args.outdir:ap.error("Cannot request both validation-only and figures")
        return 0
    if args.outdir is None:ap.error("--outdir required for figures")
    files=draw(scales,ref,comp,args.outdir.resolve())
    print("SECOND-ORDER FIGURES:",len(files),"reproducible SVG/PNG",flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
