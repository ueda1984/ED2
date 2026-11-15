#!/usr/bin/env python3
"""Make three independent, figure-ready OBC finite-size DM validation plots.

Input: committed cpc-stationary-dm-reference.csv (rounded verified CI log
data) or the full-precision physical-scan.csv from the action artifact.
Output: three individual SVG/PNG figures, no extrapolation/phase claim.
Matplotlib is required only for this optional plotting utility.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIZES=(4,6,8,10,12)
D_GRID=(-.24,-.12,0.,.12,.24)


def read_data(path):
    with path.open(newline="",encoding="utf-8") as fp:
        data=list(csv.DictReader(fp))
    series={}
    for row in data:
        def pick(*fields):
            for name in fields:
                if name in row and row[name] not in ("",None):
                    return float(row[name])
            return None
        n=int(pick("sites"))
        D=pick("D")
        datum=dict(D=D, kappa=pick("staggered_chirality_per_bond","staggered_kappa_per_bond"),
                   gap=pick("first_gap"),energy=pick("independent_E0","E0"),
                   chi=pick("chi_kappa_per_bond","chi_from_kappa_per_bond"))
        if any(datum[k] is None for k in ("kappa","gap","energy")):
            raise ValueError("Missing required scientific quantity in CSV")
        series.setdefault(n,[]).append(datum)
    if set(series)!=set(SIZES):
        raise ValueError(f"Expected sizes {SIZES}, received {sorted(series)}")
    for n in SIZES:
        series[n].sort(key=lambda r:r["D"])
        if len(series[n])!=len(D_GRID) or any(abs(x["D"]-D)>1e-8
                                               for x,D in zip(series[n],D_GRID)):
            raise ValueError("Missing or duplicate DM-scan points at size "+str(n))
        if series[n][2]["chi"] is None or series[n][2]["chi"]<0:
            raise ValueError("Missing or invalid finite-size susceptibility")
        ks=[x["kappa"] for x in series[n]]
        if any(ks[i+1]>ks[i]+2e-7 for i in range(len(ks)-1)):
            raise ValueError("Ground-energy concavity consistency violated")
    return series


def save_plot(path):
    plt.grid(alpha=.2)
    plt.tight_layout()
    plt.savefig(path.with_suffix(".svg"))
    plt.savefig(path.with_suffix(".png"),dpi=190)
    plt.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    data=read_data(args.csv)
    out=args.outdir
    out.mkdir(parents=True,exist_ok=False)

    plt.figure(figsize=(6.6,4.3))
    for n in SIZES:
        a=data[n]
        plt.plot([x["D"] for x in a],[x["kappa"] for x in a],
                 "o-",label=f"L={n}",markersize=4.5)
    plt.xlabel("Oriented staggered DM coefficient D")
    plt.ylabel("Staggered vector chirality per bond")
    plt.axhline(0,linewidth=.7,linestyle=":",color="gray")
    plt.legend(frameon=False,ncol=2)
    save_plot(out/"dm_chirality_vs_D")

    plt.figure(figsize=(6.6,4.3))
    for n in SIZES:
        a=data[n]
        plt.plot([x["D"] for x in a],[x["gap"] for x in a],
                 "o-",label=f"L={n}",markersize=4.5)
    plt.xlabel("Oriented staggered DM coefficient D")
    plt.ylabel("Lowest excitation gap (energy units)")
    plt.legend(frameon=False,ncol=2)
    save_plot(out/"dm_gap_vs_D")

    plt.figure(figsize=(6.6,4.3))
    chi=[data[n][2]["chi"] for n in SIZES]
    plt.plot(list(SIZES),chi,"o-",label="finite open chains")
    plt.xlabel("Number of sites L")
    plt.ylabel("Chiral susceptibility per bond at D=0")
    plt.xticks(list(SIZES))
    plt.legend(frameon=False)
    save_plot(out/"dm_chiral_susceptibility_vs_L")

    print(f"PLOTTING SUCCESS: 3 separate figures, 6 SVG/PNG assets, "
          f"sizes={SIZES}, chi={chi}",flush=True)


if __name__=="__main__":
    raise SystemExit(main())
