#!/usr/bin/env python3
"""QuSpin 1.0.1 comparison: identical few-magnon XXZ and lattice momenta.

This script is independent of the two Fortran executables but accepts their
source-pinned same-physics JSON report. It measures basis construction,
Hamiltonian construction and low energy separately. One Python subprocess
per case isolates peak RSS. Failures are archived; no speedup is implied.
"""
import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from scipy.sparse.linalg import eigsh

from cpc_qs3_ed_large_few_magnon import SMOKE, STRESS, sector_dim, geometry, one_energy


def quspin_case(case_index, tier, reference, output):
    from quspin.basis import spin_basis_general
    from quspin.operators import hamiltonian

    lengths, n_down, k = (SMOKE if tier == "smoke" else STRESS)[case_index]
    n = math.prod(lengths)
    dim = sector_dim(lengths, n_down, k)
    bonds, maps = geometry(lengths)
    jxy, jz = (1.0, 1.0) if lengths[1] > 1 else (-1.0, 0.8)

    kwargs = {}
    if lengths[0] > 1:
        kwargs["kxblock"] = (np.asarray(maps[0], dtype=np.int32) - 1, k[0])
    if lengths[1] > 1:
        kwargs["kyblock"] = (np.asarray(maps[1], dtype=np.int32) - 1, k[1])
    if lengths[2] > 1:
        kwargs["kzblock"] = (np.asarray(maps[2], dtype=np.int32) - 1, k[2])

    t0 = time.perf_counter()
    # QuSpin counts up spins, while QS3 specifies down-spin number.
    basis = spin_basis_general(n, S="1/2", pauli=0, Nup=n-n_down,
                               Ns_block_est=max(128, dim * 2), **kwargs)
    t1 = time.perf_counter()
    if basis.Ns != dim:
        raise AssertionError(f"QuSpin dimension {basis.Ns} not {dim}")
    plusminus = [[0.5*jxy, i-1, j-1] for i,j in bonds]
    zz = [[jz, i-1, j-1] for i,j in bonds]
    static = [["+-", plusminus], ["-+", plusminus], ["zz", zz]]
    H = hamiltonian(static, [], basis=basis, dtype=np.complex128,
                    check_herm=False, check_symm=False, check_pcon=False)
    t2 = time.perf_counter()
    if dim == 1:
        eigvals = H.eigvalsh()
        e0 = float(np.real(eigvals[0]))
        residual = 0.0
    elif dim <= 256:
        e0 = float(np.real(H.eigvalsh()[0]))
        residual = None  # complete dense eigenspectrum reference
    else:
        vals, vecs = eigsh(H.tocsr(), k=1, which="SA", tol=1e-11,
                           maxiter=15000)
        e0 = float(vals[0])
        residual = float(np.linalg.norm(H.dot(vecs[:, 0])-e0*vecs[:, 0]))
        if residual > 5e-8:
            raise AssertionError(f"Independent QuSpin Ritz residual too high {residual}")
    t3 = time.perf_counter()

    label = "x".join(map(str,lengths))+f"-n{n_down}-k"+"-".join(map(str,k))
    ref = next((row for row in reference["accepted"] if row["label"] == label), None)
    if ref is None:
        raise AssertionError(f"Missing source-pinned Fortran reference: {label}")
    error_to_original = abs(e0-ref["old"]["energy"])
    error_to_modern = abs(e0-ref["ed2"]["energy"])
    if max(error_to_original,error_to_modern) > 2e-6:
        raise AssertionError(f"QuSpin/Fortran energy discrepancy: {label} {error_to_original} {error_to_modern}")
    analytic = None
    if n_down==1:
        analytic = one_energy(lengths,k,bonds,jxy,jz)
        if abs(e0-analytic)>2e-8:
            raise AssertionError("QuSpin one-magnon analytic energy mismatch")
    result = dict(case=label,quspin_version=importlib.metadata.version("quspin"),
                  numpy_version=np.__version__,dimension=dim,sites=n,
                  n_down=n_down,k=k,couplings=dict(Jxy=jxy,Jz=jz),
                  QuSpin_E0=e0,original_QS3_E0=ref["old"]["energy"],
                  ED2_E0=ref["ed2"]["energy"],error_to_original=error_to_original,
                  error_to_modern=error_to_modern,analytic_energy=analytic,
                  residual=residual,seconds_basis=t1-t0,
                  seconds_operator=t2-t1,seconds_solver=t3-t2,
                  peak_process_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("QUSPIN MATCHED XXZ",label,"E0",e0,"dim",dim,
          "vs original",error_to_original,"vs ED2",error_to_modern,
          "basis_s",t1-t0,"operator_s",t2-t1,"solve_s",t3-t2,flush=True)


def all_cases(args):
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    cases=SMOKE if args.tier=="smoke" else STRESS
    records=[]
    errors=[]
    for index in range(len(cases)):
        result=out/f"case-{index}.json"
        argv=[sys.executable,__file__,"--case",str(index),
              "--reference",str(args.reference),"--output",str(result),
              "--tier",args.tier]
        try:
            cp=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              timeout=args.timeout,check=False)
            (out/f"case-{index}.log").write_bytes(cp.stdout)
            if cp.returncode:
                raise RuntimeError(f"Subprocess exit {cp.returncode}: {cp.stdout.decode(errors='replace')[-1600:]}")
            records.append(json.loads(result.read_text()))
            print(cp.stdout.decode(errors="replace"),flush=True)
        except Exception as exc:
            errors.append(dict(case=index,type=type(exc).__name__,message=str(exc)))
            print("QUSPIN CASE FAILURE",errors[-1],flush=True)
    ok=len(records)==len(cases)
    (out/"report.json").write_text(json.dumps(dict(
        status="passed" if ok else "incomplete",cases_expected=len(cases),
        cases_verified=len(records),records=records,failures=errors,
        comparison="same bonds, spin-1/2 XXZ, fixed Ndown and lattice momentum",
        caveat="Shared CI hardware and disparate solver initialization. "
        "Raw wall times are NOT controlled performance ratios. "
        "Basis construction and solver work must be evaluated separately."),
        indent=2)+"\n",encoding="utf-8")
    print("QU SPIN MATCHED",len(records),"/",len(cases),"success",ok,flush=True)
    return 0 if ok else 1


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--tier",choices=("smoke","stress"),default="smoke")
    p.add_argument("--case",type=int,default=None)
    p.add_argument("--timeout",type=int,default=180)
    a=p.parse_args()
    reference=json.loads(a.reference.read_text())
    if not reference.get("passed"):
        raise SystemExit("Cannot benchmark against failed/partial Fortran scientific reference")
    if a.case is None:
        raise SystemExit(all_cases(a))
    if a.case<0 or a.case>=len(SMOKE if a.tier=="smoke" else STRESS):
        p.error("Invalid case index")
    quspin_case(a.case,a.tier,reference,a.output)


if __name__=="__main__":
    main()
