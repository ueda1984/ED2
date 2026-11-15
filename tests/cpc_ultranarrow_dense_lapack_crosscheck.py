#!/usr/bin/env python3
"""Independent dense-LAPACK-driver precision audit for the ultraweak 8-site ED test.

The same *independently assembled* SciPy sparse Hamiltonian used in
the original ED2-vs-SciPy physical test is converted to a full
Hermitian 1296-by-1296 matrix and solved with distinct LAPACK
symmetric/Hermitian drivers (ZHEEVR, ZHEEVX, ZHEEVD).
All direct residuals, the first two eigenvalues, edge Sz, and
low-two-state subspace traces are checked against:
(a) SciPy sparse eigsh used by the original reference, and
(b) a separately derived second-order Q-resolvent effective model.

This audit tests the numerical *eigensolvers*, NOT a third independent
Hamiltonian specification. MRRR ('evr') and bisection ('evx') can have
closely related underlying algorithmic steps; the divide-and-conquer
'evd' is additionally included. Near a gap of 1e-7, even 1e-14
residuals do not certify individual eigenvector observables to 1e-14.

Data from these tests must not be called an arbitrary-precision
or rigorous error-bounded ground-state reference.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from cpc_stationary_interacting_sparse_doublet import (
    reference, SPINS, ULTRA_SCENARIOS,
)
from cpc_ultranarrow_second_order_reference import effective_second_order
from mixed_spin_sparse_scale import construct_sparse_hamiltonian

DRIVERS=("evr","evx","evd")
TOL_RES=1e-11
TOL_TRACE=1e-7


def direct_metrics(H,ops,pairs,energies,vectors,full_reference):
    energy=np.asarray(energies[:3],dtype=float)
    vec=np.asarray(vectors[:,:3],dtype=complex)
    norm_error=float(np.max(np.abs(vec.conj().T@vec-np.eye(3))))
    direct_res=[float(np.linalg.norm(H@vec[:,i]-energy[i]*vec[:,i]))
                for i in range(3)]
    if max(direct_res)>TOL_RES or norm_error>1e-8:
        raise AssertionError("Dense reference eigenpairs do not pass direct tests")
    gap=float(energy[1]-energy[0])
    outer=float(energy[2]-energy[1])
    edgez=float(np.vdot(vec[:,0],ops[0][2]@vec[:,0]).real)
    if not 0<gap<1e-5 or outer<.3:
        raise AssertionError("The dense solver did not recover the isolated doublet")
    # Low-cluster traces (basis invariant). Check all local magnetizations
    # and all 252 Cartesian pair correlations, not just edge Sz.
    mag=np.array([[sum(float(np.vdot(vec[:,j],ops[i][a]@vec[:,j]).real)
                       for j in range(2))/2
                   for a in range(3)] for i in range(SPINS)])
    corr=np.array([[sum(float(np.vdot(
                      vec[:,j],(ops[i][a]@ops[k][b])@vec[:,j]).real)
                      for j in range(2))/2
                    for b in range(3) for a in range(3)]
                   for i,k in pairs])
    magerr=float(np.max(np.abs(mag-full_reference["trace_mag"])))
    correrr=float(np.max(np.abs(corr-full_reference["trace_corr"])))
    if max(magerr,correrr)>TOL_TRACE:
        raise AssertionError("Density-matrix spectral projector inconsistent across solvers")
    return dict(ground_energy=float(energy[0]),first_gap=gap,external_gap=outer,
                edge_z=edgez,direct_residuals=direct_res,
                gram_max_error=norm_error,
                trace_mag_eigsh_error=magerr,
                trace_corr_eigsh_error=correrr)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--drivers",default="evr,evx,evd")
    opts=ap.parse_args()
    drivers=[v.strip() for v in opts.drivers.split(",")]
    if not drivers or len(set(drivers))!=len(drivers) or any(x not in DRIVERS for x in drivers):
        ap.error("Drivers must be unique subset of evr,evx,evd")
    root=opts.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    entries=[]
    failures=[]
    for eps,delta in ULTRA_SCENARIOS:
        try:
            model=reference(SPINS,eps,delta)
            eff=effective_second_order(model)
            H,ops=construct_sparse_hamiltonian(
                model["spins"],model["bonds"],model["fields"])
            if H.shape!=(1296,1296):
                raise AssertionError("Unexpected full Hamiltonian dimension")
            A=H.toarray()
            if np.max(np.abs(A-A.conj().T))>3e-13:
                raise AssertionError("Dense matrix is not Hermitian")
            sparse_ground_z=float(model["edge_sz"])
            each=[]
            for method in drivers:
                # LAPACK 'evd' cannot select a subset of eigenvectors.
                energy,vec=eigh(A,driver=method,
                                subset_by_index=None if method=="evd" else (0,2),
                                check_finite=False,overwrite_a=False)
                m=direct_metrics(H,ops,model["pairs"],energy,vec,model)
                m.update(driver=method,
                         second_order_gap_abs_error=abs(m["first_gap"]-eff["second_order_gap"]),
                         second_order_edge_z_abs_error=abs(m["edge_z"]-eff["second_order_dressed_edge_z"]),
                         original_eigsh_edge_z_abs_error=abs(m["edge_z"]-sparse_ground_z),
                         original_eigsh_gap_abs_error=abs(m["first_gap"]-model["ground_gap"]))
                if m["second_order_edge_z_abs_error"]>7e-8:
                    raise AssertionError("Density/reference edge magnetization differs from Q resolvent")
                each.append(m)
                print("DENSE CROSSCHECK",f"delta={delta:+.1e}",method,
                      "gap",m["first_gap"],"edge_z",m["edge_z"],
                      "versus_eff",m["second_order_edge_z_abs_error"],
                      "versus_sparse",m["original_eigsh_edge_z_abs_error"],
                      "trace_err",m["trace_corr_eigsh_error"],flush=True)
            evr=next((x for x in each if x["driver"]=="evr"),None)
            if evr is not None and evr["second_order_edge_z_abs_error"]>2e-8:
                raise AssertionError("EVR edge moment did not agree with independent Q resolvent")
            if "evr" in drivers and "evx" in drivers:
                evx=next(x for x in each if x["driver"]=="evx")
                if abs(evr["edge_z"]-evx["edge_z"])>4e-8:
                    raise AssertionError("Subset LAPACK drivers disagree unexpectedly")
            entries.append(dict(epsilon=eps,delta=delta,
                                original_sparse_eigsh_low_energies=model["low_e"],
                                original_sparse_eigsh_edge_z=sparse_ground_z,
                                second_order_gap=eff["second_order_gap"],
                                second_order_edge_z=eff["second_order_dressed_edge_z"],
                                second_order_up_Q_residual=eff["upper_Q_residual"],
                                second_order_down_Q_residual=eff["lower_Q_residual"],
                                dense=each))
        except Exception as ex:
            case=dict(epsilon=eps,delta=delta,error=str(ex),kind=type(ex).__name__)
            failures.append(case)
            print("DENSE CROSSCHECK FAILURE",case,flush=True)
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","source checkout revision"),
                scope="multi-LAPACK-driver numerical eigenvector precision at interacting 8-site ultraweak gap",
                caution="Floating-point solver cross-validation, not arbitrary-precision ground-state certification",
                drivers=drivers,expected=len(ULTRA_SCENARIOS),
                accepted=len(entries),cases=entries,failures=failures,
                pass_complete=len(entries)==len(ULTRA_SCENARIOS) and not failures)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("DENSE CROSSCHECK SUMMARY",report["accepted"],"/",report["expected"],
          "failures",len(failures),flush=True)
    return 0 if report["pass_complete"] else 1


if __name__=="__main__":
    raise SystemExit(main())
