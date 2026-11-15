#!/usr/bin/env python3
"""Validate interacting noncommuting XY avoided crossings at 8/10 sites.

Independent SciPy product-basis sparse Hamiltonian and ARPACK (k=3)
provide reference low levels, ground observables and normalized
two-state-projector traces, without using ED2 matrix construction.

ED2 ALG1/ALG2 outputs are compared with independent numbers and
direct residual checks, including the controlled near-1e-7 gap.
Status-nonzero is recorded as nonconvergence, NOT silently accepted.
The reference calculations are test-only dependencies.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh

from mixed_spin_nonu1_demo import couplings
from mixed_spin_sparse_scale import construct_sparse_hamiltonian,expectation
from mixed_spin_interacting_avoided_crossing import run_case

SIZES=(8,10)
XY_RATIO=.6
LARGE_EPS=.02
TINY_EPS=5e-7
CASES=tuple((n,LARGE_EPS,delta) for n in SIZES for delta in (-.002,0.,.002))+tuple(
    (n,TINY_EPS,0.) for n in SIZES)
SEED=1984


def reference(n,epsilon,delta):
    spins=[.5,1.]*(n//2)
    bonds,fields=couplings(n,.15)
    bonds[0]=(0,1,[epsilon,XY_RATIO*epsilon,0.,0.,0.,0.,0.,0.,0.])
    fields[0]=[0.,0.,delta]
    H,ops=construct_sparse_hamiltonian(spins,bonds,fields)
    k=3
    dim=H.shape[0]
    rng=np.random.default_rng(SEED+n)
    v0=rng.standard_normal(dim)+1j*rng.standard_normal(dim)
    v0=v0/np.linalg.norm(v0)
    values,vecs=eigsh(H,k=k,which="SA",tol=1e-12,ncv=96,
                      v0=v0,maxiter=25000)
    ind=np.argsort(values)
    E=np.array(values[ind],float)
    V=vecs[:,ind]
    residuals=[float(np.linalg.norm(H@V[:,i]-E[i]*V[:,i]))
               for i in range(k)]
    if max(residuals)>1e-8:
        raise AssertionError("Independent ARPACK low-eigenpair residual too large: "+str(residuals))
    gap=float(E[1]-E[0])
    outer=float(E[2]-E[1])
    if not (gap>1e-9 and outer>.1):
        raise AssertionError("Sparse low doublet not separated correctly "+str((gap,outer)))
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    v=V[:,0]
    local=np.array([[expectation(v,ops[i][a]) for a in range(3)]
                    for i in range(n)])
    corr=np.array([
        [expectation(v,ops[i][a]@ops[j][b]) for b in range(3) for a in range(3)]
        for i,j in pairs
    ])
    # Two-state equal-weight invariant independent of arbitrary U(2)
    # rotations in the computed low-doublet eigenbasis.
    trace_mag=np.array([
        [sum(expectation(V[:,k],ops[i][a]) for k in range(2))/2
         for a in range(3)] for i in range(n)
    ])
    trace_corr=np.array([
        [sum(expectation(V[:,k],ops[i][a]@ops[j][b]) for k in range(2))/2
         for b in range(3) for a in range(3)]
        for i,j in pairs
    ])
    rho=v.reshape(2,-1)@v.reshape(2,-1).conj().T
    linent=float(1-np.trace(rho@rho).real)
    if epsilon==LARGE_EPS and delta==0. and linent<1e-5:
        raise AssertionError("Noncommuting XY central ground unexpectedly factorized")
    if epsilon==LARGE_EPS and delta and linent<1e-6:
        raise AssertionError("Visible coupled ground unexpectedly factorized away from centre")
    # Cross-check total-Sz nonconservation without assuming that the
    # applied edge field itself is nonconserving.
    total_sz=sum((o[2] for o in ops),start=H*0.)
    defect=H@total_sz-total_sz@H
    relative_comm=float(np.linalg.norm(defect.data)/np.linalg.norm(H.data))
    if relative_comm<.01:
        raise AssertionError("Chosen core unexpectedly preserves total Sz")
    return dict(
        n_sites=n,epsilon=epsilon,delta=delta,spins=spins,
        bonds=bonds,fields=fields,pairs=pairs,dimension=dim,
        H_nonzeros=int(H.nnz),
        low_e=E.tolist(),ground_e=float(E[0]),
        ground_gap=gap,outer_gap=outer,
        reference_residuals=residuals,
        edge_sz=float(local[0,2]),edge_linear_entropy=linent,
        mag=local,corr=corr,trace_mag=trace_mag,trace_corr=trace_corr,
        nonu1_commutator=relative_comm
    )


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--omp",type=int,choices=(1,2),required=True)
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    records=[]
    models=[]
    violations=[]
    for n,eps,delta in CASES:
        print(f"SPARSE REFERENCE n={n} eps={eps:.3g} "
              f"delta={delta:+.3g}",flush=True)
        model=reference(n,eps,delta)
        models.append(model)
        print(f"SPARSE REFERENCE DONE n={n} eps={eps:.3g} delta={delta:+.3g} "
              f"dimension={model['dimension']} nnz={model['H_nonzeros']} "
              f"E0={model['ground_e']:.14g} gap={model['ground_gap']:.9g} "
              f"outergap={model['outer_gap']:.7g} "
              f"edgeSz={model['edge_sz']:+.9g} "
              f"entropy={model['edge_linear_entropy']:.4g} "
              f"max_reference_residual={max(model['reference_residuals']):.3e}",
              flush=True)
    for n in SIZES:
        for eps in (LARGE_EPS,):
            data={m["delta"]:m for m in models if m["n_sites"]==n and m["epsilon"]==eps}
            if not (data[0.]["ground_gap"]<min(
                    data[-.002]["ground_gap"],data[+.002]["ground_gap"])
                    and data[-.002]["edge_sz"]>.1 and data[+.002]["edge_sz"]<-.1):
                raise AssertionError("Physical large-system avoided crossing not established")
    for model in models:
        configurations=[(1,0,0,1e-8),
                        (2,2,16,1e-8),
                        (2,4,24,1e-8)]
        if model["epsilon"]==TINY_EPS:
            configurations.append((1,0,0,1e-12))
        for alg,keep,cap,atol in configurations:
            try:
                row=run_case(exe,root,model,alg,keep,cap,atol,args.omp)
                records.append(row)
                print("SPARSE CHECK",row["case"],"status",row["status"],
                      "outcome",row["classification"],
                      "max_r",max(row.get("direct_residuals",[float('inf')])),
                      "mag_err",row.get("local_magnetization_error"),
                      "corr_err",row.get("two_body_correlation_error"),
                      "trace_err",row.get("cluster_ensemble_mag_error"),
                      "Hcalls",row["h_applications"],flush=True)
            except Exception as ex:
                record=dict(n_sites=model["n_sites"],epsilon=model["epsilon"],
                            delta=model["delta"],alg=alg,keep=keep,capacity=cap,
                            atol=atol,error=str(ex))
                violations.append(record)
                print("SPARSE CONTRACT FAILURE",record,flush=True)
    summary=dict(schema=1,
                 source_sha=os.environ.get("GITHUB_SHA","record in CI"),
                 kind="interacting anisotropic XY sparse reference for 8/10 sites",
                 omp=args.omp,ref_backend="independent SciPy CSR and ARPACK eigsh(k=3)",
                 independent_models=[
                     dict(n_sites=x["n_sites"],epsilon=x["epsilon"],delta=x["delta"],
                          dimension=x["dimension"],nonzeros=x["H_nonzeros"],
                          energies=x["low_e"],gap=x["ground_gap"],
                          outer_gap=x["outer_gap"],edge_sz=x["edge_sz"],
                          edge_linear_entropy=x["edge_linear_entropy"],
                          max_reference_residual=max(x["reference_residuals"]))
                     for x in models
                 ],cases=records,contract_violations=violations,
                 expected=26,completed=len(records),
                 accepted=sum(x["status"]==0 for x in records),
                 nonconverged=sum(x["status"]!=0 for x in records),
                 individual_ground_unresolved=sum(
                     x["classification"]=="accepted_low_doublet_but_individual_ground_unresolved"
                     for x in records))
    (root/"report.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(f"SPARSE INTERACTING SUMMARY: {summary['completed']}/26 trials; "
          f"{summary['accepted']} accepted; {summary['nonconverged']} nonconverged; "
          f"{summary['individual_ground_unresolved']} individual ground unresolved; "
          f"{len(violations)} contract violations; OMP={args.omp}",flush=True)
    return 0 if summary["completed"]==26 and not violations else 1


if __name__=="__main__":
    raise SystemExit(main())
