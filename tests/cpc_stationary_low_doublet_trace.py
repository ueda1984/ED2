#!/usr/bin/env python3
"""Eight-site sparse-reference low-doublet trace stress test for ED2.

A deliberately decoupled edge spin creates an exact or nearly split low
doublet while the other seven sites keep the *stationary* non-U(1)
alternating mixed-spin XYZ/Gamma/staggered-DM Hamiltonian. The independent
reference uses SciPy sparse Kronecker products and an eigenproblem on the
seven-site active cluster, not ED2's basis enumeration or sparse kernel.

For each converged ED2 ALG=2 run, both low eigenvectors are independently
read back through ED2's native wavefunction I/O. The equal-weight
magnetization/correlation trace is compared with the independently
constructed two-dimensional spectral projector. Individual eigenvector
magnetizations are NOT presumed invariant when the low gap is tiny.

This is an adversarial numerical test, NOT the size-stationary
physical-response model of the CPC main example. Removing the first bond
and field explicitly changes the original physical chain. Keep them
separate in the manuscript.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mixed_spin_dm_chirality_scan import create_model
from mixed_spin_sparse_scale import construct_sparse_hamiltonian, expectation
from mixed_spin_near_degenerate import run_case
from mixed_spin_doublet_trace import readback_second

SIZES=(8,)
DELTAS=(0.0,1e-10,1e-6)
WORKSPACES=((2,16),(4,24))
D=0.12
RESIDUAL_ATOL=1e-8
TRACE_TOL=2e-5
ENERGY_TOL=4e-7


def sparse_reference(n,delta):
    if n not in (8,10) or delta<0 or delta>1e-4:
        raise ValueError("Only the validated small-gap stress profile is supported")
    spins,bonds,fields=create_model(n,D,profile="stationary")
    # Detach the left edge S=1/2; retain the non-U(1) stationary core.
    core_bonds=[(i,j,c) for i,j,c in bonds if i>0 and j>0]
    fields[0]=[0.,0.,delta]
    active_bonds=[(i-1,j-1,c) for i,j,c in core_bonds]
    active_H,active_ops=construct_sparse_hamiltonian(spins[1:],active_bonds,fields[1:])
    dim=active_H.shape[0]
    if dim*2!=6**(n//2):
        raise AssertionError("Unexpected product Hilbert dimension")
    rng=np.random.default_rng(20261115+n)
    initial=rng.standard_normal(dim)+1j*rng.standard_normal(dim)
    initial/=np.linalg.norm(initial)
    vals,vecs=eigsh(active_H,k=2,which="SA",tol=1e-13,
                    ncv=min(dim-1,80),maxiter=40000,v0=initial)
    ix=np.argsort(vals)
    vals=np.asarray(vals[ix],float)
    vecs=vecs[:,ix]
    core_residuals=[float(np.linalg.norm(active_H@vecs[:,j]-vals[j]*vecs[:,j]))
                    for j in range(2)]
    if max(core_residuals)>1e-9:
        raise AssertionError("Independent active-chain eigenpair has large residual")
    gap_active=float(vals[1]-vals[0])
    if not .02<gap_active<2. or delta>=.001*gap_active:
        raise AssertionError("Isolated two-state cluster not sufficiently separated")

    # The exact two-dimensional low-cluster projector is I_edge x |g><g|.
    # Numerically, this avoids an unstable basis choice inside the doublet.
    ground_core=vecs[:,0]
    eye_edge=np.eye(2,dtype=complex)
    Pcols=np.column_stack([np.kron(eye_edge[:,j],ground_core)
                           for j in range(2)])
    full_H,ops=construct_sparse_hamiltonian(spins,core_bonds,fields)
    if full_H.shape!=(2*dim,2*dim):
        raise AssertionError("Unexpected sparse Hamiltonian dimension")
    core_expanded=sparse.kron(sparse.eye(2,format="csr",dtype=complex),
                              active_H,format="csr")
    edge_z=sparse.kron(sparse.diags([.5,-.5],format="csr"),
                       sparse.eye(dim,format="csr",dtype=complex),
                       format="csr")
    factorization=full_H-core_expanded-delta*edge_z
    factor_error=float(np.max(np.abs(factorization.data))) if factorization.nnz else 0.
    if factor_error>1e-12:
        raise AssertionError("Independent full sparse Hamiltonian factorization failed")
    reference_eigs=[float(vals[0]-delta/2),float(vals[0]+delta/2),
                    float(vals[1]-delta/2),float(vals[1]+delta/2)]
    low_pair_gap=delta
    outer_gap=reference_eigs[2]-reference_eigs[1]
    if outer_gap<.01:
        raise AssertionError("Two lowest states not spectrally isolated")
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    mag=np.array([[float(sum(np.vdot(Pcols[:,k],op@Pcols[:,k]).real
                              for k in range(2))/2)
                   for op in ops[i]] for i in range(n)])
    corr=np.array([[float(sum(np.vdot(Pcols[:,k],
                                    (ops[i][a]@ops[j][b])@Pcols[:,k]).real
                               for k in range(2))/2)
                    for b in range(3) for a in range(3)]
                   for i,j in pairs])
    if abs(mag[0,2])>1e-12:
        raise AssertionError("Exact projector must have zero free-edge Sz")
    # Verify the full system is non-U(1); isolated Sz is a *local*
    # accidental conservation law introduced only as a stress probe.
    total_sz=sum((site[2] for site in ops),
                 start=sparse.csr_matrix(full_H.shape,dtype=complex))
    comm=full_H@total_sz-total_sz@full_H
    relative_comm=float(np.linalg.norm(comm.data)/np.linalg.norm(full_H.data))
    if relative_comm<=.01:
        raise AssertionError("The active core must remain non-U(1)")
    gram=Pcols.conj().T@Pcols
    if np.max(np.abs(gram-np.eye(2)))>1e-11:
        raise AssertionError("Independent projector columns are not orthonormal")
    return dict(spins=spins,bonds=core_bonds,fields=fields,dimension=full_H.shape[0],
                pairs=pairs,eigs=reference_eigs,gap=low_pair_gap,outside_gap=outer_gap,
                free_spin_expected=-.5 if delta>0 else None,
                free_spin_ensemble_z=float(mag[0,2]),
                matrix_error=factor_error,commutator_relative=relative_comm,
                ensemble_magnetization=mag,ensemble_correlations=corr,
                core_reference_residuals=core_residuals,
                core_excitation_gap=gap_active)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--omp",type=int,required=True,choices=(1,2))
    ap.add_argument("--sizes",default="8")
    ap.add_argument("--deltas",default="0,1e-10,1e-6")
    ap.add_argument("--workspaces",default="2:16,4:24")
    a=ap.parse_args()
    sizes=[int(x.strip()) for x in a.sizes.split(",")]
    deltas=[float(x.strip()) for x in a.deltas.split(",")]
    wks=[]
    for token in a.workspaces.split(","):
        k,c=map(int,token.split(":"))
        if not 1<k<c<=128 or (k,c) in wks:
            ap.error("Invalid/duplicated restart workspace")
        wks.append((k,c))
    if sizes!=[8] or len(set(deltas))!=len(deltas) or not deltas or any(
            d not in DELTAS for d in deltas):
        ap.error("This controlled test supports only eight-site and listed deltas")
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(exist_ok=False,parents=True)
    results=[]
    failures=[]
    references=[]
    for n in sizes:
        for delta in deltas:
            ref=sparse_reference(n,delta)
            references.append(dict(sites=n,delta=delta,full_dimension=ref["dimension"],
                                   independent_low_eigs=ref["eigs"],
                                   internal_gap=ref["gap"],outer_gap=ref["outside_gap"],
                                   sparse_factor_error=ref["matrix_error"],
                                   relative_nonU1_commutator=ref["commutator_relative"],
                                   active_reference_residuals=ref["core_reference_residuals"]))
            for keep,capacity in wks:
                parent=root/f"n{n}-d{delta:.0e}-k{keep}-m{capacity}"
                parent.mkdir(exist_ok=False)
                case=f"n{n}-delta{delta:.0e}-k{keep}-m{capacity}-omp{a.omp}"
                try:
                    ed=run_case(exe,parent,n,delta,alg=2,keep=keep,
                                capacity=capacity,atol=RESIDUAL_ATOL,
                                omp=a.omp,ref=ref,suffix="-sparse-doublet")
                    if ed["status"]!=0 or not ed.get("inside_low_doublet"):
                        raise AssertionError("ED2 failed low doublet inclusion")
                    err_e=max(abs(x-y) for x,y in
                              zip(ed["returned_energies"],ref["eigs"][:2]))
                    if err_e>ENERGY_TOL:
                        raise AssertionError("Low-doublet energies disagree with independent reference")
                    src=parent/ed["case"]
                    m1=np.loadtxt(src/"local_mag.dat",ndmin=2)
                    c1=np.loadtxt(src/"two_body_cf_xyz.dat",ndmin=2)
                    m2,c2=readback_second(exe,src,parent/"native-second",a.omp)
                    if m1.shape!=(n,4) or m2.shape!=m1.shape:
                        raise AssertionError("Unexpected local magnetization columns")
                    if c1.shape!=(len(ref["pairs"]),11) or c2.shape!=c1.shape:
                        raise AssertionError("Unexpected correlation columns")
                    np.testing.assert_array_equal(m1[:,0],m2[:,0])
                    np.testing.assert_array_equal(c1[:,:2],c2[:,:2])
                    average_mag=(m1[:,1:4]+m2[:,1:4])/2
                    average_corr=(c1[:,2:11]+c2[:,2:11])/2
                    magerr=float(np.max(np.abs(average_mag-ref["ensemble_magnetization"])))
                    correrr=float(np.max(np.abs(average_corr-ref["ensemble_correlations"])))
                    if max(magerr,correrr)>TRACE_TOL:
                        raise AssertionError("Native ED2 doublet ensemble differs from sparse projector "+
                                             str((magerr,correrr,ref["outside_gap"])))
                    # A single vector may be an arbitrary doublet superposition.
                    first_z=float(m1[0,3]);second_z=float(m2[0,3])
                    if abs(first_z)>0.500001 or abs(second_z)>0.500001:
                        raise AssertionError("Free spin-half magnetization out of range")
                    rec=dict(case=case,sites=n,delta=delta,omp=a.omp,
                             keep=keep,capacity=capacity,ed2_status=ed["status"],
                             lowest_eigenvalue_max_error=err_e,
                             max_direct_residual=max(ed["direct_residuals"]),
                             independent_external_gap=ref["outside_gap"],
                             conditional_cluster_leakage_upper=ed["conditional_low_subspace_leakage_upper"],
                             gram_error=ed["orthonormality_max_error"],
                             first_edge_z=first_z,second_edge_z=second_z,
                             doublet_mean_edge_z=float(average_mag[0,2]),
                             independent_trace_mag_max_error=magerr,
                             independent_trace_corr_max_error=correrr,
                             h_applications=ed["h_applications"],status="passed")
                    results.append(rec)
                    print("SPARSE DOUBLET PASS",case,
                          "gap",ref["gap"],"external_gap",ref["outside_gap"],
                          "edgeSz",first_z,second_z,
                          "max_trace_errors",magerr,correrr,flush=True)
                except Exception as error:
                    item=dict(case=case,error=str(error),kind=type(error).__name__)
                    failures.append(item)
                    print("SPARSE DOUBLET FAILURE",item,flush=True)
    expected=len(sizes)*len(deltas)*len(wks)
    good=len(results)==expected and not failures
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record checkout SHA"),
                scope="size-8 non-U1 stationary-core doublet with decoupled edge probe",
                caution="adversarial spectral validation, not thermodynamic phase evidence",
                omp=a.omp,reference=references,results=results,
                expected=expected,passed=len(results),failures=failures,
                pass_complete=good)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("SPARSE DOUBLET SUMMARY",len(results),"/",expected,
          "failures",len(failures),"OMP",a.omp,flush=True)
    return 0 if good else 1


if __name__=="__main__":
    raise SystemExit(main())
