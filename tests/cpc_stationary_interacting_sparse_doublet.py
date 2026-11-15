#!/usr/bin/env python3
"""Independent sparse validation of an *interacting* eight-site low doublet.

The first S=1/2 is weakly coupled to a non-U(1) seven-site mixed-spin
XYZ/Gamma/staggered-DM core by epsilon*(S1x*S2x+eta*S1y*S2y).
A tunable delta*S1z shifts the narrow avoided gap. Unlike the
previous isolated-edge control, the low-doublet projector cannot
be written as a free-spin tensor-product projector.

The complete independent SciPy sparse Hamiltonian supplies E0,E1,E2,
the ground-state observables, and all 2-state trace observables.
ED2 computes two low states using ALG=2. Both native ED2 eigenvectors
are passed through the actual ED2 observable I/O (second-state
re_wf=1 readback), and the average is compared with the independent
low-projector trace. Report nonconvergence or state-specific mismatch
honestly; never infer lowest-state inclusion from a residual alone.

This model modifies the edge bond/field on purpose and MUST NOT be
called the unmodified connected two-site-stationary DM response model.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mixed_spin_dm_chirality_scan import create_model
from mixed_spin_sparse_scale import construct_sparse_hamiltonian
from mixed_spin_interacting_avoided_crossing import run_case

SPINS=8
DM=0.12
EDGE_Y_FRACTION=0.6
SCENARIOS=(
    (.02,-.002),(.02,0.),(.02,+.002),
    (5e-4,0.),(5e-4,2e-5),
)
ULTRA_SCENARIOS=((5e-7,-1e-7),(5e-7,0.),(5e-7,+1e-7))
RESTARTS=((2,16),(4,24))
ATOL=1e-8


def hermitian_mean(matrix,operator):
    v=np.vdot(matrix,operator@matrix)
    if abs(v.imag)>5e-10:
        raise AssertionError("Hermitian expectation acquired imaginary part")
    return float(v.real)


def reference(n,epsilon,delta):
    if n!=SPINS or (epsilon,delta) not in SCENARIOS+ULTRA_SCENARIOS:
        raise ValueError("Only reviewed eight-site interacting stress cases supported")
    spins,bonds,fields=create_model(n,DM,profile="stationary")
    if len(bonds)!=n-1 or bonds[0][:2]!=(0,1):
        raise AssertionError("Unexpected stationary model connectivity")
    # Keep every original core bond, and make the first bond weakly XY
    # interacting (not an exact tensor-product edge).
    bonds[0]=(0,1,[epsilon,EDGE_Y_FRACTION*epsilon,0.,0.,0.,0.,0.,0.,0.])
    fields[0]=[0.,0.,delta]
    H,ops=construct_sparse_hamiltonian(spins,bonds,fields)
    dim=H.shape[0]
    if dim!=6**(n//2):
        raise AssertionError("Wrong eight-site full Hilbert dimension")
    rng=np.random.default_rng(20261115+int(epsilon*1e8)+int((delta+.01)*1e8))
    v0=rng.normal(size=dim)+1j*rng.normal(size=dim)
    v0/=np.linalg.norm(v0)
    vals,vecs=eigsh(H,k=3,which="SA",tol=1e-12,
                    ncv=min(100,dim-1),v0=v0,maxiter=40000)
    order=np.argsort(vals)
    energies=np.asarray(vals[order],float)
    states=vecs[:,order]
    independent_residuals=[float(np.linalg.norm(H@states[:,j]-energies[j]*states[:,j]))
                           for j in range(3)]
    if max(independent_residuals)>5e-9:
        raise AssertionError("Independent SciPy eigenpairs insufficiently accurate")
    if np.max(np.abs(states.conj().T@states-np.eye(3)))>5e-8:
        raise AssertionError("Independent eigenvectors not mutually orthogonal")
    gap=float(energies[1]-energies[0])
    outer=float(energies[2]-energies[1])
    if not 1e-9<gap<.03 or outer<.05 or gap>.1*outer:
        raise AssertionError("Incorrect narrow interacting low-doublet structure "+str((gap,outer)))
    if (epsilon,delta) in SCENARIOS and gap<=1e-6:
        raise AssertionError("Incorrect narrow interacting low-doublet structure "+str((gap,outer)))
    # The edge cannot be factored out. Each of S_edge^{x,y,z}
    # must have a nonzero commutator with the Hamiltonian.
    edge_comm=[]
    for a in range(3):
        comm=H@ops[0][a]-ops[0][a]@H
        edge_comm.append(float(np.linalg.norm(comm.data)/np.linalg.norm(H.data)))
    if min(edge_comm)<(1e-10 if epsilon<1e-5 else 1e-7):
        raise AssertionError("Accidentally conserved edge-spin component "+str(edge_comm))
    total_z=sum((o[2] for o in ops),
                start=sparse.csr_matrix(H.shape,dtype=complex))
    comm=H@total_z-total_z@H
    nonu1=float(np.linalg.norm(comm.data)/np.linalg.norm(H.data))
    if nonu1<.01:
        raise AssertionError("Interacting model unexpectedly conserves total Sz")
    # A separate *first-order low-energy effective Hamiltonian*
    # checks the ultraweak coupling without relying solely on
    # eigenvector conditioning in the full near-doublet:
    # H_eff = e_core*I + (delta*sigmaz + epsilon*mx*sigmax
    #                    + eta*epsilon*my*sigmay)/2.
    # The gap prediction is accurate up to O(epsilon^2/external_gap)
    # for epsilon much smaller than the core excitation gap.
    effective_gap=None
    effective_edge_z=None
    effective_gap_error=None
    if epsilon<1e-5:
        core_spin=spins[1:]
        core_bonds=[(i-1,j-1,c) for i,j,c in bonds if i>0 and j>0]
        H_core,core_ops=construct_sparse_hamiltonian(
            core_spin,core_bonds,fields[1:])
        rng_core=np.random.default_rng(1979)
        core_guess=rng_core.normal(size=H_core.shape[0])+1j*rng_core.normal(size=H_core.shape[0])
        core_guess/=np.linalg.norm(core_guess)
        ecore,vec_core=eigsh(H_core,k=1,which="SA",tol=1e-13,
                             ncv=min(H_core.shape[0]-1,80),
                             maxiter=40000,v0=core_guess)
        phi=vec_core[:,0]
        core_res=float(np.linalg.norm(H_core@phi-ecore[0]*phi))
        if core_res>1e-10:
            raise AssertionError("Effective core ground state is not accurate")
        mx=float(np.vdot(phi,core_ops[0][0]@phi).real)
        my=float(np.vdot(phi,core_ops[0][1]@phi).real)
        effective_gap=float(np.sqrt(delta**2 +
                             epsilon**2*(mx**2+
                                           EDGE_Y_FRACTION**2*my**2)))
        effective_edge_z=float(-delta/(2*effective_gap))
        effective_gap_error=float(abs(gap-effective_gap))
        if effective_gap_error>1e-9:
            raise AssertionError("Ultraweak full sparse gap deviates from independent effective two-state model: "+
                                 str((gap,effective_gap,effective_gap_error)))
    ground=states[:,0]
    cluster=states[:,:2]
    rho_edge=ground.reshape(2,-1)@ground.reshape(2,-1).conj().T
    linear_entropy=float(1.-np.trace(rho_edge@rho_edge).real)
    if epsilon==.02 and linear_entropy<1e-6:
        raise AssertionError("Moderate XY edge bond did not entangle the ground state")

    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    mag=np.array([[hermitian_mean(ground,ops[i][a])
                   for a in range(3)] for i in range(n)])
    corr=np.array([[hermitian_mean(ground,ops[i][a]@ops[j][b])
                    for b in range(3) for a in range(3)]
                   for i,j in pairs])
    trace_mag=np.array([[sum(hermitian_mean(cluster[:,j],ops[i][a])
                             for j in range(2))/2
                         for a in range(3)] for i in range(n)])
    trace_corr=np.array([[sum(hermitian_mean(cluster[:,j],ops[i][a]@ops[k][b])
                              for j in range(2))/2
                          for b in range(3) for a in range(3)]
                         for i,k in pairs])
    # Physical projector invariance under an explicitly constructed
    # complex unitary rotation of the two reference eigenvectors.
    theta=.391
    phase=np.exp(1j*.274)
    U=np.array([[np.cos(theta),-phase.conjugate()*np.sin(theta)],
                [phase*np.sin(theta),np.cos(theta)]],dtype=complex)
    if np.max(np.abs(U.conj().T@U-np.eye(2)))>1e-14:
        raise AssertionError("Invalid independent cluster rotation")
    rotated=cluster@U
    qerror=float(np.max(np.abs(cluster@cluster.conj().T -
                               rotated@rotated.conj().T)))
    if qerror>2e-13:
        raise AssertionError("Independent projector not basis invariant")

    return dict(
        n_sites=n,spins=spins,bonds=bonds,fields=fields,pairs=pairs,
        epsilon=epsilon,delta=delta,edge_y_fraction=EDGE_Y_FRACTION,
        dimension=dim,ground_e=float(energies[0]),
        low_e=energies.tolist(),ground_gap=gap,outer_gap=outer,
        mag=mag,corr=corr,trace_mag=trace_mag,trace_corr=trace_corr,
        edge_sz=float(mag[0,2]),edge_linear_entropy=linear_entropy,
        connected_edge_xx=float(corr[pairs.index((0,1)),0] -
                                mag[0,0]*mag[1,0]),
        commutator_relative=nonu1,edge_component_commutators=edge_comm,
        reference_residuals=independent_residuals,
        projector_rotation_error=qerror,
        effective_low_gap=effective_gap,
        effective_edge_z_firstorder=effective_edge_z,
        effective_low_gap_error=effective_gap_error,
    )


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--omp",required=True,type=int,choices=(1,2))
    a=ap.parse_args()
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    results=[]; failures=[]; refs=[]
    for epsilon,delta in SCENARIOS:
        m=reference(SPINS,epsilon,delta)
        refs.append(dict(sites=SPINS,epsilon=epsilon,delta=delta,
                         low_energies=m["low_e"],inner_gap=m["ground_gap"],
                         outside_gap=m["outer_gap"],
                         edge_linear_entropy=m["edge_linear_entropy"],
                         independent_residuals=m["reference_residuals"],
                         edge_commutators=m["edge_component_commutators"],
                         nonu1_relative=m["commutator_relative"],
                         projector_rotation_error=m["projector_rotation_error"]))
        print("INTERACTING SPARSE REFERENCE",
              "eps",epsilon,"delta",delta,"inner",m["ground_gap"],
              "outer",m["outer_gap"],"entropy",m["edge_linear_entropy"],
              flush=True)
        for keep,capacity in RESTARTS:
            case=f"eps{epsilon:.4g}-d{delta:+.4g}-k{keep}-m{capacity}-omp{a.omp}"
            try:
                ed=run_case(exe,root,m,alg=2,keep=keep,capacity=capacity,
                            atol=ATOL,omp=a.omp,suffix="-independent-sparse")
                if ed["status"]!=0:
                    raise AssertionError("Two-lowest eigenpairs did not converge")
                if not ed.get("inside_low_cluster"):
                    raise AssertionError("Returned pairs not included in independently verified doublet")
                if ed.get("cluster_ensemble_mag_error") is None:
                    raise AssertionError("Missing native second-state observable readback")
                if ed["cluster_ensemble_mag_error"]>1e-5 or ed["cluster_ensemble_cf_error"]>1e-5:
                    raise AssertionError("Doublet trace physical reference disagreement")
                if max(ed["direct_residuals"])>ATOL*(1+1e-12):
                    raise AssertionError("Accepted pair exceeds requested direct residual")
                records=dict(case=case,sites=SPINS,omp=a.omp,epsilon=epsilon,
                             delta=delta,NOK=keep,NOM=capacity,
                             low_gap=m["ground_gap"],outer_gap=m["outer_gap"],
                             reference_edge_entropy=m["edge_linear_entropy"],
                             reference_edge_z=m["edge_sz"],
                             returned_edge_z=ed.get("measured_edge_sz"),
                             ground_classification=ed["classification"],
                             first_energy_error=ed["energy_ground_error"],
                             first_mag_max_error=ed["local_magnetization_error"],
                             first_corr_max_error=ed["two_body_correlation_error"],
                             trace_mag_max_error=ed["cluster_ensemble_mag_error"],
                             trace_corr_max_error=ed["cluster_ensemble_cf_error"],
                             direct_residuals=ed["direct_residuals"],
                             h_applications=ed["h_applications"],status="passed")
                results.append(records)
                print("INTERACTING SPARSE PASS",case,
                      "gap",m["ground_gap"],"traceerrs",
                      records["trace_mag_max_error"],records["trace_corr_max_error"],
                      "res",records["direct_residuals"],flush=True)
            except Exception as e:
                issue=dict(case=case,error=str(e),kind=type(e).__name__)
                failures.append(issue)
                print("INTERACTING SPARSE FAILURE",issue,flush=True)
    expected=len(SCENARIOS)*len(RESTARTS)
    passed=len(results)==expected and not failures
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record checkout SHA"),
                purpose="interacting eight-site sparse ED low-doublet projector/trace validation",
                note="edge bond is modified; this is not the connected stationary chirality model",
                omp=a.omp,reference=refs,results=results,failures=failures,
                expected=expected,passed=len(results),pass_complete=passed)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("INTERACTING SPARSE SUMMARY",len(results),"/",expected,
          "failures",len(failures),"OMP",a.omp,flush=True)
    return 0 if passed else 1


if __name__=="__main__":
    raise SystemExit(main())
