#!/usr/bin/env python3
"""Zero-bias weak-edge polarization scaling: analytical symmetry + LAPACK.

Stress model only: weak anisotropic XY edge bond connected to the
remaining non-U(1) seven-site XYZ/Gamma/staggered-DM alternating chain.
The PRIMARY CPC physical application keeps all original stationary bonds.

Purpose:
 (1) verify H(-epsilon) = (sigma_z(edge) x I) H(epsilon)
                              (sigma_z(edge) x I), precisely;
 (2) compare independent dense-LAPACK ground Sz / low gap with the
     first- and second-order seven-site core-projected effective model;
 (3) test a robust limiting scaling for epsilon>0, delta=0:
      <S_edge^z> = C * epsilon + O(epsilon**2), epsilon -> 0+,
     where C follows from *asymmetric second-order self-energies*.

At delta=0, the first-order effective Hamiltonian predicts Sz=0,
so the nonzero polarization is a sensitive numerical acceptance test,
not a claim of new bulk physics or a general solver algorithm.

This finite-precision cross-check is not a rigorous interval error bound.
"""
import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from mixed_spin_dm_chirality_scan import create_model
from mixed_spin_sparse_scale import construct_sparse_hamiltonian
from cpc_ultranarrow_second_order_reference import effective_second_order

EPSILONS=(5e-7,2e-6,5e-6,1e-5)
ETA=.6
D=.12
N=8


def model_and_matrix(epsilon):
    spins,bonds,fields=create_model(N,D,profile="stationary")
    bonds[0]=(0,1,[epsilon,ETA*epsilon,0.,0.,0.,0.,0.,0.,0.])
    fields[0]=[0.,0.,0.]
    H,ops=construct_sparse_hamiltonian(spins,bonds,fields)
    model=dict(n_sites=N,spins=spins,bonds=bonds,fields=fields,
               epsilon=epsilon,delta=0.)
    return model,H,ops


def crosscheck(epsilon):
    if epsilon not in EPSILONS:
        raise ValueError("Unsupported epsilon: verify tested precision window")
    model,H,ops=model_and_matrix(epsilon)
    inv_model,H_minus,_=model_and_matrix(-epsilon)
    dim=H.shape[0]
    if dim!=1296 or inv_model["epsilon"]!=-epsilon:
        raise AssertionError("Unexpected mixed-spin eight-site Hilbert dimension")
    core_dim=dim//2
    # Rotation by pi around the edge z spin: Sx -> -Sx, Sy -> -Sy,
    # but Sz and all core operators are invariant.
    edge_rotation=sparse.diags(
        np.r_[np.ones(core_dim),-np.ones(core_dim)],
        format="csr",dtype=complex)
    diff=(H_minus-edge_rotation@H@edge_rotation).tocsr()
    diff.eliminate_zeros()
    sign_symmetry_error=float(np.max(np.abs(diff.data))) if diff.nnz else 0.
    if sign_symmetry_error>1e-12:
        raise AssertionError("Exact weak-bond sign symmetry failed")
    second=effective_second_order(model)
    # A distinct LAPACK eigensolver, not a rerun of the existing
    # SciPy sparse eigsh driver. The physical matrix definition
    # is shared; independence is numerical, not Hamiltonian-level.
    energies,vectors=eigh(H.toarray(),driver="evr",subset_by_index=(0,2),
                          check_finite=False)
    direct=[float(np.linalg.norm(H@vectors[:,j]-energies[j]*vectors[:,j]))
            for j in range(3)]
    gram=float(np.max(np.abs(vectors.conj().T@vectors-np.eye(3))))
    gap=float(energies[1]-energies[0])
    outside=float(energies[2]-energies[1])
    moment=float(np.vdot(vectors[:,0],ops[0][2]@vectors[:,0]).real)
    if max(direct)>1e-10 or gram>1e-10:
        raise AssertionError("Dense independent reference eigenpair residual/Gram failed")
    if not (1e-8<gap<1e-4 and outside>.3):
        raise AssertionError("Unexpected dense weak-edge doublet")
    if moment<=0 or moment>.0001:
        raise AssertionError("Unexpected zero-field chiral-edge moment")
    gaperr1=float(abs(gap-second["first_order_gap"]))
    gaperr2=float(abs(gap-second["second_order_gap"]))
    zerr=float(abs(moment-second["second_order_dressed_edge_z"]))
    if gaperr2>1e-10:
        raise AssertionError("Core-projected second-order gap does not match dense model")
    if zerr>max(2e-8,0.02*abs(moment)):
        raise AssertionError("Core-projected edge magnetization does not match dense model")
    if abs(second["first_order_edge_z"])>1e-12:
        raise AssertionError("First-order zero-bias edge magnetization is not zero")
    # Sigma_up - Sigma_down = epsilon**2 * eta *
    # Im <QXg | (Hcore-eg)^(-1) | QYg>.  The corresponding
    # limiting edge polarization slope is -delta_sigma/
    # (2*gap_first*epsilon), not a new independent fit parameter.
    sigma_diff=second["sigma_difference"]
    slope_core=-sigma_diff/(2*second["first_order_gap"]*epsilon)
    if not math.isfinite(slope_core) or slope_core<=0:
        raise AssertionError("Nonphysical asymptotic polarization coefficient")
    if abs(moment/epsilon-slope_core)>max(.15,0.05*slope_core):
        raise AssertionError("Zero-bias polarization does not show first-order epsilon scaling")
    return dict(epsilon=epsilon,reference_method="LAPACK ZHEEVR",
                independent_dense_gap=gap,external_gap=outside,
                independent_dense_edge_z=moment,
                dense_eigenpair_direct_residuals=direct,
                dense_gram_max_error=gram,exact_sign_symmetry_error=sign_symmetry_error,
                first_order_gap=second["first_order_gap"],
                second_order_gap=second["second_order_gap"],
                second_order_dressed_edge_z=second["second_order_dressed_edge_z"],
                first_order_gap_abs_error=gaperr1,second_order_gap_abs_error=gaperr2,
                second_order_edge_z_abs_error=zerr,
                asymptotic_edge_polarization_coefficient=slope_core,
                observed_edge_z_per_epsilon=moment/epsilon,
                second_order_Q_residual_max=max(second["upper_Q_residual"],
                                                second["lower_Q_residual"]),
                core_gap=second["core_excitation_gap"])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    rows=[];errors=[]
    for epsilon in EPSILONS:
        try:
            row=crosscheck(epsilon)
            rows.append(row)
            print("WEAK EDGE SCALE PASS","epsilon",epsilon,
                  "dense_z",row["independent_dense_edge_z"],
                  "slope",row["observed_edge_z_per_epsilon"],
                  "coefficient",row["asymptotic_edge_polarization_coefficient"],
                  "gap",row["independent_dense_gap"],
                  "gap_2nd_error",row["second_order_gap_abs_error"],
                  "edge_z_2nd_error",row["second_order_edge_z_abs_error"],
                  flush=True)
        except Exception as exc:
            item=dict(epsilon=epsilon,error=str(exc),type=type(exc).__name__)
            errors.append(item)
            print("WEAK EDGE SCALE FAILURE",item,flush=True)
    if rows:
        slopes=[x["observed_edge_z_per_epsilon"] for x in rows]
        relspan=(max(slopes)-min(slopes))/sum(slopes)*len(slopes)
        if relspan>0.05:
            errors.append(dict(type="ScalingViolation",relative_span=relspan,
                               description="Polarization per unit epsilon is not stable"))
    else:
        relspan=None
    pass_complete=len(rows)==len(EPSILONS) and not errors
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record checkout SHA"),
                scope="modified weak-edge interacting 8-site model, not connected main DM chain",
                caution="Finite-precision LAPACK and perturbative results; no interval-certified accuracy",
                expected=len(EPSILONS),passed=len(rows),
                epsilon_positive=EPSILONS,observed_slope_relative_span=relspan,
                results=rows,violations=errors,pass_complete=pass_complete)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("WEAK EDGE SCALING SUMMARY",len(rows),"/",len(EPSILONS),
          "violations",len(errors),flush=True)
    return 0 if pass_complete else 1


if __name__=="__main__":
    raise SystemExit(main())
