#!/usr/bin/env python3
"""CPC physics validation: sparse independent ED for 8/10-site mixed spins.

Builds an *independent* SciPy sparse Hamiltonian in the ordinary product
basis, diagonalizes it with ARPACK (eigsh), and compares ED2 ALG=1/2
results, local magnetization and all Cartesian two-site correlations.
NumPy/SciPy are test-only tools, not ED2 runtime dependencies.

The script refuses nonconvergence and preserves every run's input and logs.
It reports timings and memory from a single environment as observations,
NOT as reproducible throughput/speedup estimates.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

import mixed_spin_nonu1_demo as ref
from mixed_spin_iterative_crosscheck import read_diagnostics, RESIDUAL_ABS

E_ATOL=3e-7
OBS_ATOL=8e-6
SEED=19


def site_ops_sparse(spins):
    local = [[sparse.csr_matrix(x) for x in ref.local_spin_operators(s)] for s in spins]
    dims = [len(local[i][0].indptr)-1 for i in range(len(spins))]
    operators = []
    for i in range(len(spins)):
        site=[]
        for a in range(3):
            mat=sparse.csr_matrix([[1.+0.j]])
            for j,dim in enumerate(dims):
                mat=sparse.kron(mat,local[j][a] if i==j
                                else sparse.identity(dim,format="csr",dtype=complex),
                                format="csr")
            site.append(mat)
        operators.append(site)
    return operators,dims


def construct_sparse_hamiltonian(spins,bonds,fields):
    ops,dims=site_ops_sparse(spins)
    d=int(np.prod(dims))
    H=sparse.csr_matrix((d,d),dtype=complex)
    for i,field in enumerate(fields):
        for a in range(3):
            if field[a]:H+=field[a]*ops[i][a]
    for i,j,c in bonds:
        for a in range(3):
            if c[a]:H+=c[a]*(ops[i][a] @ ops[j][a])
        for k,a,b in ((0,1,2),(1,2,0),(2,0,1)):
            if c[3+k]:H+=c[3+k]*(ops[i][a] @ ops[j][b]
                                 - ops[i][b] @ ops[j][a])
            if c[6+k]:H+=c[6+k]*(ops[i][a] @ ops[j][b]
                                 + ops[i][b] @ ops[j][a])
    H.sum_duplicates()
    H.eliminate_zeros()
    discrepancy=(H-H.getH()).tocoo()
    err=float(np.max(np.abs(discrepancy.data))) if discrepancy.nnz else 0.
    if err>2e-12:
        raise AssertionError("Independent sparse H is non-Hermitian "+str(err))
    return H,ops


def expectation(q,A):
    z=np.vdot(q,A@q)
    if abs(z.imag)>2e-8:
        raise AssertionError("Nonreal Hermitian expectation: "+str(z))
    return float(z.real)


def crosscheck_sparse_and_dense():
    """Ensure SciPy's separate sparse construction equals the prior dense model."""
    records=[]
    for n in (4,6):
        spins=[0.5,1.]*(n//2)
        bonds,fields=ref.couplings(n,0.15)
        Hs,_=construct_sparse_hamiltonian(spins,bonds,fields)
        dense_ops,_=ref.product_operators(spins)
        Hd=ref.dense_hamiltonian(dense_ops,bonds,fields)
        maximum=float(np.max(np.abs(Hs.toarray()-Hd)))
        if maximum>2e-13:
            raise AssertionError("Independent sparse/dense construction differs "+
                                 str((n,maximum)))
        records.append(dict(n_sites=n,dimension=Hs.shape[0],max_matrix_error=maximum))
    return records


def reference_model(n,hx):
    spins=[0.5,1.]*(n//2)
    bonds,fields=ref.couplings(n,hx)
    start=time.perf_counter()
    H,ops=construct_sparse_hamiltonian(spins,bonds,fields)
    dt_build=time.perf_counter()-start
    start=time.perf_counter()
    values,vectors=eigsh(H,k=2,which="SA",tol=1e-12,ncv=64,maxiter=20000)
    order=np.argsort(values)
    energies=np.array(values[order],float)
    vectors=vectors[:,order]
    dt_ref=time.perf_counter()-start
    residuals=[float(np.linalg.norm(H@vectors[:,j]-energies[j]*vectors[:,j]))
               for j in range(2)]
    if max(residuals)>1e-8 or energies[1]-energies[0]<1e-4:
        raise AssertionError("Independent sparse solver did not establish distinct low eigenpairs")
    v=vectors[:,0]
    lm=np.array([[expectation(v,op) for op in site] for site in ops])
    total_sz=sum((site[2] for site in ops),start=sparse.csr_matrix(H.shape,dtype=complex))
    comm=H@total_sz-total_sz@H
    comm_norm=float(np.linalg.norm(comm.data))
    H_norm=float(np.linalg.norm(H.data))
    relative_commutator=comm_norm/H_norm
    if relative_commutator < 1e-5:
        raise AssertionError("Model expected to break U(1), but [H,Sz] is small")
    weights=np.array([0.4*((-1)**i)+0.09*(i+1) for i in range(n)])
    V=sum((weights[i]*ops[i][0] for i in range(n)),
          start=sparse.csr_matrix(H.shape,dtype=complex))
    delta=1e-4
    eval_plus=float(eigsh(H+delta*V,k=1,which="SA",tol=1e-11,
                          ncv=48,maxiter=20000,return_eigenvectors=False)[0])
    eval_minus=float(eigsh(H-delta*V,k=1,which="SA",tol=1e-11,
                           ncv=48,maxiter=20000,return_eigenvectors=False)[0])
    hf_finite_difference=(eval_plus-eval_minus)/(2*delta)
    hf_exact=float(weights@lm[:,0])
    hf_error=abs(hf_finite_difference-hf_exact)
    if hf_error>5e-6:
        raise AssertionError("Independent Hellmann-Feynman identity failed "+
                             str((hf_finite_difference,hf_exact,hf_error)))
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    cf=np.array([[expectation(v,ops[i][a] @ ops[j][b])
                  for b in range(3) for a in range(3)] for i,j in pairs])
    return dict(H=H,spin=spins,bonds=bonds,fields=fields,pairs=pairs,
                dimension=H.shape[0],energies=energies,mag=lm,corr=cf,
                eigenpair_residuals=residuals,
                independent_build_seconds=dt_build,
                independent_eigensolve_seconds=dt_ref,nnz=int(H.nnz),
                relative_sz_commutator_frobenius=relative_commutator,
                hellmann_feynman_finite_difference=hf_finite_difference,
                hellmann_feynman_exact_spin_expectation=hf_exact,
                hellmann_feynman_error=hf_error)


def run_ed2(exe,root,n,hx,alg,omp,reference):
    name=f"mixed-n{n}-hx{hx:.2f}-alg{alg}-omp{omp}"
    work=root/name
    work.mkdir(exist_ok=False)
    bonds,fields,spins,pairs=(reference[x] for x in ("bonds","fields","spin","pairs"))
    (work/"bonds.dat").write_text("".join(
        f"{i+1} {j+1} "+" ".join(f"{q:.17g}" for q in c)+"\n"
        for i,j,c in bonds))
    (work/"fields.dat").write_text("".join(
        f"{i+1} "+" ".join(f"{q:.17g}" for q in field)+"\n"
        for i,field in enumerate(fields)))
    (work/"spins.dat").write_text("".join(f"{s:.1f}\n" for s in spins))
    (work/"caps.dat").write_text("".join(f"{round(2*s)}\n" for s in spins))
    (work/"pairs.dat").write_text("".join(f"{i+1} {j+1}\n" for i,j in pairs))
    cap=sum(round(2*s) for s in spins)
    deck=f"""&input_parameters
NOS={n}, NODmin=0, NODmax={cap}, ALG={alg},
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={n}, NO_two={len(bonds)}, wk_dim=-1, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={len(pairs)}, FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-12, minitr=5, maxitr=700, itrint=10
/
&input_TRLan
NOE=2, NOK=10, NOM=64, maxitr=260, lnc_ene_conv=1d-12,
i_vec_min=1, i_vec_max=2
/
"""
    (work/"input.dat").write_text(deck)
    env=dict(os.environ,OMP_NUM_THREADS=str(omp),OPENBLAS_NUM_THREADS="1",
             MKL_NUM_THREADS="1",ED2_RNG_SEED=str(SEED),
             ED2_DIAGNOSTICS="diagnostics.txt",
             ED2_RESIDUAL_ATOL=str(RESIDUAL_ABS))
    timefile=work/"gnu-time.txt"
    timeexe="/usr/bin/time"
    use_time=Path(timeexe).exists()
    cmd=[str(exe)]
    if use_time:
        cmd=[timeexe,"-f","MAX_RSS_KIB=%M","-o",str(timefile)]+cmd
    tic=time.perf_counter()
    with (work/"input.dat").open("rb") as inp:
        proc=subprocess.run(cmd,cwd=work,stdin=inp,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=240)
    elapsed=time.perf_counter()-tic
    (work/"stdout.log").write_bytes(proc.stdout)
    if proc.returncode!=0:
        raise AssertionError("ED2 failed: "+proc.stdout.decode(errors="replace")[-2500:])
    requested=1 if alg==1 else 2
    diag,pairrows=read_diagnostics(work/"diagnostics.txt",requested)
    if int(diag["dimension"])!=reference["dimension"]:
        raise AssertionError("Incorrect ED2 Hilbert-space dimension")
    eigvals=np.loadtxt(work/"eigenvalues.dat",ndmin=2)
    if eigvals.shape!=(requested,2):raise AssertionError("Wrong eigenvalue output shape")
    np.testing.assert_array_equal(eigvals[:,0],np.arange(1,requested+1))
    energies=eigvals[:,1]
    eerr=float(np.max(np.abs(energies-reference["energies"][:requested])))
    if eerr>E_ATOL:raise AssertionError("Wrong lowest eigenvalue(s): "+str((eerr,energies)))
    m=np.loadtxt(work/"local_mag.dat",ndmin=2)
    cf=np.loadtxt(work/"two_body_cf_xyz.dat",ndmin=2)
    np.testing.assert_array_equal(m[:,0],np.arange(1,n+1))
    np.testing.assert_array_equal(cf[:,:2],
                                  np.array([(i+1,j+1) for i,j in pairs]))
    weights=np.array([0.4*((-1)**i)+0.09*(i+1) for i in range(n)])
    hf_ed2=float(weights @ m[:,1])
    hf_ed2_error=abs(hf_ed2-reference["hellmann_feynman_finite_difference"])
    if hf_ed2_error>5e-6:
        raise AssertionError("ED2 magnetization inconsistent with independent "
                             "energy derivative "+str(hf_ed2_error))
    magerr=float(np.max(np.abs(m[:,1:4]-reference["mag"])))
    cferr=float(np.max(np.abs(cf[:,2:11]-reference["corr"])))
    if max(magerr,cferr)>OBS_ATOL:
        raise AssertionError("Wrong independently referenced observables: "+
                             str((magerr,cferr)))
    peak_rss=None
    if use_time and timefile.exists():
        for line in timefile.read_text().splitlines():
            if line.startswith("MAX_RSS_KIB="):
                peak_rss=int(line.split("=",1)[1])
    return dict(case=name,alg=alg,omp=omp,n_sites=n,dimension=reference["dimension"],
                nnz=reference["nnz"],hx=hx,reference_eigenvalues=reference["energies"].tolist(),
                ed2_eigenvalues=energies.tolist(),energy_max_error=eerr,
                magnetization_max_error=magerr,correlation_max_error=cferr,
                max_reported_residual=max(x["residual"] for x in pairrows),
                relative_sz_commutator_frobenius=reference["relative_sz_commutator_frobenius"],
                hellmann_feynman_ed2_error=hf_ed2_error,
                reference_max_residual=max(reference["eigenpair_residuals"]),
                h_applications=int(diag["h_applications"]),iterations=int(diag["iterations"]),
                cycles=int(diag["cycles"]),wall_seconds=elapsed,peak_rss_kib=peak_rss,
                status="passed")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--threads",default=1,type=int,choices=(1,2))
    p.add_argument("--sizes",default="8,10,12",
                   help="Comma-delimited site counts, currently 8, 10, 12")
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    sizes=tuple(int(s.strip()) for s in args.sizes.split(","))
    if any(n not in (8,10,12) for n in sizes) or len(set(sizes))!=len(sizes):
        raise ValueError("Expected distinct site counts from 8,10,12")
    scenarios=[(n,hx) for n in sizes for hx in ((0.15,0.38) if n==8 else (0.15,))]
    sparse_dense_check=crosscheck_sparse_and_dense()
    print("SPARSE-DENSE INDEPENDENT SELFTEST:",sparse_dense_check,flush=True)
    records=[]
    for n,hx in scenarios:
        print(f"REFERENCE: {n} sites, hx control {hx}",flush=True)
        reference=reference_model(n,hx)
        print(f"REFERENCE DONE: dimension={reference['dimension']} "
              f"nnz={reference['nnz']} "
              f"eigvals={reference['energies'].tolist()} "
              f"residual={max(reference['eigenpair_residuals']):.3e} "
              f"commutator_rel={reference['relative_sz_commutator_frobenius']:.3e} "
              f"HF_error={reference['hellmann_feynman_error']:.3e}",flush=True)
        for alg in (1,2):
            name=f"mixed-n{n}-hx{hx:.2f}-alg{alg}-omp{args.threads}"
            try:
                record=run_ed2(exe,root,n,hx,alg,args.threads,reference)
                records.append(record)
                print(f"PASS {name}: ΔE={record['energy_max_error']:.2e} "
                      f"ΔM={record['magnetization_max_error']:.2e} "
                      f"ΔCF={record['correlation_max_error']:.2e} "
                      f"residual={record['max_reported_residual']:.2e} "
                      f"Hcalls={record['h_applications']} "
                      f"wall={record['wall_seconds']:.3f}s "
                      f"RSS={record['peak_rss_kib']}",flush=True)
            except Exception as error:
                records.append(dict(case=name,status="failed",error=str(error)))
                print(f"FAIL {name}: {error}",flush=True)
    result=dict(schema=1,kind="independent sparse mixed-spin ED2 benchmark",
                git_sha=os.environ.get("GITHUB_SHA","not provided"),
                sparse_dense_matrix_selftest=sparse_dense_check,
                omp=args.threads,sizes=list(sizes),
                total=len(records),passed=sum(x["status"]=="passed" for x in records),
                failed=sum(x["status"]=="failed" for x in records),
                scenarios=records)
    (root/"report.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(f"SPARSE CROSSCHECK: {result['passed']}/{result['total']} passed",flush=True)
    return 0 if result["failed"]==0 and result["total"]==2*len(scenarios) else 1


if __name__=="__main__":
    raise SystemExit(main())
