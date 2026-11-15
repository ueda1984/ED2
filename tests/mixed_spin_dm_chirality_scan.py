#!/usr/bin/env python3
"""DM-conjugate vector-chirality application for generalized mixed-spin ED2.

Model: open S=(1/2,1,1/2,1,...) chain with anisotropic, nonuniform
XYZ exchange, Gamma_xy, site-dependent transverse/longitudinal fields.
The oriented nearest-neighbor Dz coefficient is replaced by D*(-1)^i
(zero-based site i). No artificial fixed-Sz symmetry is imposed.

Physical identity (Hellmann-Feynman):
  dE0/dD = sum_i (-1)^i <S_i^x S_(i+1)^y - S_i^y S_(i+1)^x>.

Independent SciPy construction verifies the spectrum, complete
Cartesian NN correlators, bond chirality, and numerical energy
derivative, then checks original ED2 outputs (whose Cartesian
column order is XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ).
The 4/6-site sparse matrices/eigenvalues are additionally checked
against a separately built *dense* NumPy Hamiltonian.

No claim of a new magnetic phase or algorithmic speedup.
NumPy/SciPy are test-only and not required by ED2.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

import mixed_spin_nonu1_demo as dense
from mixed_spin_sparse_scale import construct_sparse_hamiltonian,expectation
from mixed_spin_near_degenerate import read_diag

D_VALUES=(-.24,-.12,0.,.12,.24)
ALL_SIZES=(4,6,8,10,12)
FD_STEP=2e-4
SEED=19
RES_TOL=1e-8
E_TOL=3e-7
OBS_TOL=8e-6
HF_TOL=3e-5

CSV_FIELDS=("profile","sites","dimension","D","E0","E0_per_site",
            "first_gap","staggered_kappa","staggered_kappa_per_bond",
            "mean_abs_bond_kappa","mean_Sx","mean_Sz",
            "hf_central_difference","hf_error","chi_from_kappa_per_bond",
            "chi_from_curvature_per_bond","chi_spectral_per_bond",
            "chi_selfconsistency_error","max_reference_residual")


def create_model(n,D,profile="diagnostic"):
    spins=[.5,1.]*(n//2)
    if profile=="diagnostic":
        # Historical regression model; J/Gamma/field drift with absolute
        # site index, so DO NOT extrapolate it to the thermodynamic limit.
        bonds,fields=dense.couplings(n,.15)
        prepared=[]
        for i,j,c in bonds:
            if j!=i+1 or len(c)!=9:
                raise AssertionError("Unexpected bond ordering/coupling record")
            cc=list(c)
            cc[5]=float(((-1)**i)*D)
            prepared.append((i,j,cc))
        return spins,prepared,fields
    if profile!="stationary":
        raise ValueError("Unknown physical model profile "+str(profile))

    # Exactly the same two-site-periodic local rules for every size:
    # alternating S=1/2,1 and two alternating oriented bond types.
    # All model coefficients remain bounded as L grows.
    bonds=[]
    for i in range(n-1):
        eta=(-1)**i
        c=[.77+.05*eta,.63+.03*eta,1.09-.04*eta,
           0.,0.,D*eta,
           0.,0.,.065+.015*eta]
        bonds.append((i,i+1,c))
    fields=[]
    for i in range(n):
        eta=(-1)**i
        fields.append([.15*(.4*eta+.1),0.,.17+.035*eta])
    return spins,bonds,fields


def reference(n,D,check_dense=False,profile="diagnostic"):
    spins,bonds,fields=create_model(n,D,profile)
    started=time.perf_counter()
    H,ops=construct_sparse_hamiltonian(spins,bonds,fields)
    build_seconds=time.perf_counter()-started
    dim=H.shape[0]
    pairs=[(i,i+1) for i in range(n-1)]
    rng=np.random.default_rng(1984+n)
    v0=rng.standard_normal(dim)+1j*rng.standard_normal(dim)
    v0=v0/np.linalg.norm(v0)
    started=time.perf_counter()
    vals,vecs=eigsh(H,k=2,which="SA",tol=1e-12,
                    ncv=min(80,dim-1),v0=v0,maxiter=30000)
    order=np.argsort(vals)
    E=np.asarray(vals[order],float)
    V=vecs[:,order]
    solve_seconds=time.perf_counter()-started
    residual=[float(np.linalg.norm(H@V[:,k]-E[k]*V[:,k])) for k in range(2)]
    if max(residual)>5e-9 or E[1]-E[0]<.005:
        raise AssertionError("Independent sparse two-level reference inaccurate "+
                             str((n,D,residual,E)))
    psi=V[:,0]
    lm=np.array([[expectation(psi,op) for op in site] for site in ops])
    corr=np.array([[expectation(psi,ops[i][a]@ops[j][b])
                    for b in range(3) for a in range(3)] for i,j in pairs])
    # Reference correlation is ordered [xx,yx,zx,xy,yy,zy,xz,yz,zz],
    # so kappa^z=xy-yx is index 3 MINUS index 1.
    bond_kappa=np.asarray(corr[:,3]-corr[:,1])
    signs=np.array([(-1)**i for i in range(n-1)],float)
    staggered=float(signs@bond_kappa)
    K=sparse.csr_matrix(H.shape,dtype=complex)
    for i,j in pairs:
        K=K+signs[i]*(ops[i][0]@ops[j][1]-ops[i][1]@ops[j][0])
    full_kappa=expectation(psi,K)
    if abs(full_kappa-staggered)>1e-10:
        raise AssertionError("Component-wise chirality disagrees with direct operator")
    if D==0.0:
        # First derivative from two independent energy solves.
        ep,vp=eigsh(H+FD_STEP*K,k=1,which="SA",v0=v0,tol=1e-12,
                    ncv=min(80,dim-1),maxiter=30000)
        em,vm=eigsh(H-FD_STEP*K,k=1,which="SA",v0=v0,tol=1e-12,
                    ncv=min(80,dim-1),maxiter=30000)
        eplus=float(ep[0])
        eminus=float(em[0])
        derivative=float((eplus-eminus)/(2*FD_STEP))
        hf_err=abs(derivative-staggered)
        if hf_err>HF_TOL:
            raise AssertionError("Independent HF mismatch "+
                                 str((n,derivative,staggered,hf_err)))
        # Second derivative via TWO independent routes:
        # 1. curvature of the ground-state energy,
        # 2. derivative of the DM-conjugate observable.
        kp=expectation(vp[:,0],K)
        km=expectation(vm[:,0],K)
        chi_from_k=float(-(kp-km)/(2*FD_STEP*(n-1)))
        chi_from_energy=float(-(eplus-2*E[0]+eminus)/(FD_STEP**2*(n-1)))
        chi_self_error=abs(chi_from_k-chi_from_energy)
        if chi_from_k< -1e-8 or chi_from_energy< -1e-6 or chi_self_error>2e-4:
            raise AssertionError("Concavity/chiral susceptibility consistency failed "+
                                 str((n,chi_from_k,chi_from_energy,chi_self_error)))
    else:
        derivative=None
        hf_err=None
        chi_from_k=None
        chi_from_energy=None
        chi_self_error=None

    dense_err=None
    dense_eigs_err=None
    chi_from_spectrum=None
    if check_dense:
        full_ops,_=dense.product_operators(spins)
        Hd=dense.dense_hamiltonian(full_ops,bonds,fields)
        dense_err=float(np.max(np.abs(H.toarray()-Hd)))
        dense_eigs,dense_vectors=np.linalg.eigh(Hd)
        dense_eigs_err=float(np.max(np.abs(dense_eigs[:2]-E)))
        if dense_err>5e-13 or dense_eigs_err>5e-10:
            raise AssertionError("Dense and sparse constructions disagree "+
                                 str((dense_err,dense_eigs_err)))
        # Independent finite difference of actual Hamiltonian input
        # coefficients checks the D parameter and its directed sign.
        if D==0.:
            # Third, direct spectral route for the susceptibility,
            # exact within the small-system dense diagonalization.
            matrix_K=K.toarray()
            amplitudes=dense_vectors[:,1:].conj().T @ (matrix_K @ dense_vectors[:,0])
            chi_from_spectrum=float(2*np.sum(
                np.abs(amplitudes)**2/(dense_eigs[1:]-dense_eigs[0]))/(n-1))
            if abs(chi_from_spectrum-chi_from_k)>5e-5:
                raise AssertionError("Dense spectral susceptibility does not match finite difference")
            _,plus_bonds,_=create_model(n,+FD_STEP,profile)
            _,minus_bonds,_=create_model(n,-FD_STEP,profile)
            dense_K=(dense.dense_hamiltonian(full_ops,plus_bonds,fields)-
                     dense.dense_hamiltonian(full_ops,minus_bonds,fields))/(2*FD_STEP)
            operator_err=float(np.max(np.abs(dense_K-K.toarray())))
            if operator_err>1e-10:
                raise AssertionError("Directed DM derivative operator mismatch "+str(operator_err))

    sz=sum((site[2] for site in ops),start=sparse.csr_matrix(H.shape,dtype=complex))
    comm=H@sz-sz@H
    comm_ratio=float(np.linalg.norm(comm.data)/np.linalg.norm(H.data))
    if comm_ratio<.01:
        raise AssertionError("Physical model unexpectedly conserves total Sz")

    return dict(
        profile=profile,n_sites=n,D=float(D),dimension=dim,spins=spins,bonds=bonds,fields=fields,pairs=pairs,
        H_nonzeros=int(H.nnz),max_sparse_row_nonzeros=int(np.max(np.diff(H.indptr))),
        energies=E.tolist(),E0=float(E[0]),
        gap=float(E[1]-E[0]),reference_residual=residual,
        lm=lm,corr=corr,bond_kappa=bond_kappa,staggered_kappa=staggered,
        mean_abs_kappa=float(np.mean(np.abs(bond_kappa))),
        mean_Sx=float(np.mean(lm[:,0])),mean_Sz=float(np.mean(lm[:,2])),
        derivative=derivative,hellmann_feynman_error=hf_err,
        chi_from_kappa_per_bond=chi_from_k,
        chi_from_curvature_per_bond=chi_from_energy,
        chi_spectral_per_bond=chi_from_spectrum,
        chi_selfconsistency_error=chi_self_error,
        reference_build_seconds=build_seconds,reference_solve_seconds=solve_seconds,
        dense_matrix_error=dense_err,dense_eigenvalue_error=dense_eigs_err,
        nonU1_commutator_relative=comm_ratio,
    )


def write_inputs(work,r,alg,keep=4,capacity=24,wk_dim=-1):
    n=r["n_sites"]
    (work/"bonds.dat").write_text(
        "".join(f"{i+1} {j+1} "+" ".join(f"{q:.17g}" for q in c)+"\n"
                for i,j,c in r["bonds"]))
    (work/"fields.dat").write_text(
        "".join(f"{i+1} "+" ".join(f"{q:.17g}" for q in xyz)+"\n"
                for i,xyz in enumerate(r["fields"])))
    (work/"spins.dat").write_text("".join(f"{s:.1f}\n" for s in r["spins"]))
    (work/"caps.dat").write_text("".join(f"{int(round(2*s))}\n" for s in r["spins"]))
    (work/"pairs.dat").write_text("".join(f"{i+1} {j+1}\n" for i,j in r["pairs"]))
    nodmax=sum(int(round(2*s)) for s in r["spins"])
    deck=f"""&input_parameters
NOS={n}, NODmin=0, NODmax={nodmax}, ALG={alg},
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={n}, NO_two={len(r['bonds'])}, wk_dim={wk_dim}, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={len(r['pairs'])},
FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-12, minitr=5, maxitr=1200, itrint=10
/
&input_TRLan
NOE=2, NOK={keep}, NOM={capacity}, maxitr=320, lnc_ene_conv=1d-12,
i_vec_min=1, i_vec_max=2
/
"""
    (work/"input.dat").write_text(deck)


def run_ed2(exe,root,r,alg,threads,keep=4,capacity=24,
            wk_dim=-1,initial_text=None):
    n=r["n_sites"]
    D=r["D"]
    sign="zero" if D==0 else f"{D:+.2f}".replace("+","p").replace("-","m")
    tag=f"{r['profile']}-n{n}-D{sign}-alg{alg}-omp{threads}"
    work=root/tag
    work.mkdir(exist_ok=False)
    write_inputs(work,r,alg,keep=keep,capacity=capacity,wk_dim=wk_dim)
    if initial_text is not None:
        (work/'initial-vector.txt').write_text(initial_text)
    timelog=work/"gnu-time.txt"
    cmd=[str(exe)]
    if Path("/usr/bin/time").is_file():
        cmd=["/usr/bin/time","-f","MAX_RSS_KIB=%M","-o",str(timelog)]+cmd
    env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             ED2_RNG_SEED=str(SEED),ED2_RESIDUAL_ATOL=str(RES_TOL),
             ED2_DIAGNOSTICS="diagnostics.txt")
    if initial_text is not None:
        env['ED2_INITIAL_VECTOR']=str(work/'initial-vector.txt')
    start=time.perf_counter()
    with (work/"input.dat").open("rb") as inp:
        proc=subprocess.run(cmd,stdin=inp,cwd=work,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=320)
    wall=time.perf_counter()-start
    (work/"stdout.log").write_bytes(proc.stdout)
    if not (work/"diagnostics.txt").exists():
        raise AssertionError("ED2 did not write solver diagnostic sidecar "+
                             proc.stdout.decode(errors="replace")[-2200:])
    diag,pairs=read_diag(work/"diagnostics.txt")
    status=int(diag["status"])
    if (status==0)!=(proc.returncode==0):
        raise AssertionError("Solver status and process exit disagree")
    wanted=1 if alg==1 else 2
    if int(diag["wanted"])!=wanted or int(diag["dimension"])!=r["dimension"]:
        raise AssertionError("Unexpected dimension/eigenpair count")
    if diag.get("residual_policy_enabled")!="T" or diag.get("seed_specified")!="T":
        raise AssertionError("Explicit diagnostic contract not activated")
    if initial_text is not None and diag.get('explicit_initial')!='T':
        raise AssertionError('Explicit starting vector not honored')
    if alg==2 and (int(diag['keep'])!=keep or int(diag['capacity'])!=capacity):
        raise AssertionError('Requested ALG=2 restart workspace not honored')
    rss=None
    if timelog.is_file():
        for line in timelog.read_text().splitlines():
            if line.startswith("MAX_RSS_KIB="): rss=int(line.split("=",1)[1])
    row=dict(case=tag,n_sites=n,D=D,dimension=r["dimension"],
             alg=alg,wanted=wanted,threads=threads,
             requested_keep=keep if alg==2 else None,
             requested_capacity=capacity if alg==2 else None,
             requested_wk_dim=wk_dim,
             explicit_initial=initial_text is not None,
             status=status,reason=diag["reason"],
             h_applications=int(diag["h_applications"]),
             iterations=int(diag["iterations"]),cycles=int(diag["cycles"]),
             wall_seconds=wall,peak_rss_kib=rss,
             absolute_residual_request=RES_TOL,
             lowest_state_inclusion=diag.get("lowest_state_inclusion"))
    if status!=0:
        if status not in (1,2):
            raise AssertionError("Unexpected numerical solver exit")
        if any((work/f).exists() for f in
               ("eigenvalues.dat","local_mag.dat","two_body_cf_xyz.dat")):
            raise AssertionError("Unsuccessful solver wrote accepted numerical output")
        row["classification"]="honest_nonconvergence"
        return row
    if len(pairs)!=wanted or any(p["index"]!=j+1 for j,p in enumerate(pairs)):
        raise AssertionError("Successful solver omitted required pair diagnostics")
    if any(p["available"]!="T" or p["certified"]!="T" or
           p["residual"]>RES_TOL*(1+2e-12) for p in pairs):
        raise AssertionError("Direct residual policy contradicted success")
    if diag.get("all_requested_residuals_met")!="T":
        raise AssertionError("Successful solver residual policy not satisfied")
    eig=np.loadtxt(work/"eigenvalues.dat",ndmin=2)
    if eig.shape!=(wanted,2):
        raise AssertionError("Wrong eigenpair count in output")
    np.testing.assert_array_equal(eig[:,0],np.arange(1,wanted+1))
    eigerr=float(np.max(np.abs(eig[:,1]-r["energies"][:wanted])))
    if eigerr>E_TOL:
        raise AssertionError("ED2 low eigenvalue disagrees with independent reference "+
                             str((n,D,alg,eigerr)))
    lm=np.loadtxt(work/"local_mag.dat",ndmin=2)
    corr=np.loadtxt(work/"two_body_cf_xyz.dat",ndmin=2)
    if lm.shape!=(n,4) or corr.shape!=(n-1,11):
        raise AssertionError("Wrong local-magnetization/NN correlation file shape")
    np.testing.assert_array_equal(lm[:,0],np.arange(1,n+1))
    np.testing.assert_array_equal(corr[:,:2],np.array(
        [(i+1,i+2) for i in range(n-1)]))
    magerr=float(np.max(np.abs(lm[:,1:4]-r["lm"])))
    correrr=float(np.max(np.abs(corr[:,2:11]-r["corr"])))
    if max(magerr,correrr)>OBS_TOL:
        raise AssertionError("ED2 physical observables disagree "+
                             str((n,D,alg,magerr,correrr)))
    # File columns are 1-based [site_i,site_j,xx,yx,zx,xy,...].
    # In the zero-based numpy array: xy=5, yx=3.
    kappa_from_ed2=corr[:,5]-corr[:,3]
    kappa_error=float(np.max(np.abs(kappa_from_ed2-r["bond_kappa"])))
    signs=np.array([(-1)**i for i in range(n-1)],float)
    staggered_from_ed2=float(signs@kappa_from_ed2)
    staggered_error=abs(staggered_from_ed2-r["staggered_kappa"])
    if max(kappa_error,staggered_error)>OBS_TOL*(n-1):
        raise AssertionError("ED2 DM conjugate chirality mismatched independent reference")
    if D==0:
        hf_ed2_error=abs(staggered_from_ed2-r["derivative"])
        if hf_ed2_error>HF_TOL:
            raise AssertionError("ED2 DM chirality contradicts independent dE/dD")
    else:
        hf_ed2_error=None
    row.update(classification="passed",
               ed2_low_eigenvalues=eig[:,1].tolist(),
               eigenenergy_max_error=eigerr,
               local_magnetization_max_error=magerr,
               NN_cartesian_correlation_max_error=correrr,
               bond_chirality_max_error=kappa_error,
               staggered_chirality_error=staggered_error,
               staggered_chirality_from_ed2=staggered_from_ed2,
               hellmann_feynman_ed2_error=hf_ed2_error,
               max_direct_residual=max(p["residual"] for p in pairs))
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--sizes",default="4,6,8,10",
                   help="Subset of 4,6,8,10,12 (comma separated)")
    p.add_argument("--omp",default=1,type=int,choices=(1,2))
    p.add_argument("--profile",default="diagnostic",choices=("diagnostic","stationary"))
    args=p.parse_args()
    sizes=tuple(int(s.strip()) for s in args.sizes.split(","))
    if not sizes or any(n not in ALL_SIZES for n in sizes) or len(set(sizes))!=len(sizes):
        raise ValueError("Expected unique even site counts drawn from 4,6,8,10,12")
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    records=[]
    run_rows=[]
    failures=[]
    summary=[]
    for n in sizes:
        for D in D_VALUES:
            ref=reference(n,D,check_dense=(n<=6),profile=args.profile)
            records.append(ref)
            summary.append(dict(
                profile=args.profile,sites=n,dimension=ref["dimension"],D=D,E0=ref["E0"],
                E0_per_site=ref["E0"]/n,first_gap=ref["gap"],
                staggered_kappa=ref["staggered_kappa"],
                staggered_kappa_per_bond=ref["staggered_kappa"]/(n-1),
                mean_abs_bond_kappa=ref["mean_abs_kappa"],
                mean_Sx=ref["mean_Sx"],mean_Sz=ref["mean_Sz"],
                hf_central_difference=ref["derivative"],
                hf_error=ref["hellmann_feynman_error"],
                chi_from_kappa_per_bond=ref["chi_from_kappa_per_bond"],
                chi_from_curvature_per_bond=ref["chi_from_curvature_per_bond"],
                chi_spectral_per_bond=ref["chi_spectral_per_bond"],
                chi_selfconsistency_error=ref["chi_selfconsistency_error"],
                max_reference_residual=max(ref["reference_residual"])
            ))
            print(f"DM REFERENCE n={n} D={D:+.2f} dim={ref['dimension']} "
                  f"E0={ref['E0']:+.12g} gap={ref['gap']:.8g} "
                  f"Kst={ref['staggered_kappa']:+.10g} "
                  f"Kst_per_bond={ref['staggered_kappa']/(n-1):+.10g} "
                  f"HF_error={ref['hellmann_feynman_error']} "
                  f"chi_per_bond={ref['chi_from_kappa_per_bond']} "
                  f"chi_spectrum={ref['chi_spectral_per_bond']} "
                  f"max_ref_r={max(ref['reference_residual']):.3e}",flush=True)
    # Variational concavity: E0(D) is the minimum of functions
    # linear in D, hence it is concave and K(D)=dE0/dD is nonincreasing.
    for n in sizes:
        ordered=[ref["staggered_kappa"] for ref in records if ref["n_sites"]==n]
        if any(ordered[j+1]>ordered[j]+1e-7 for j in range(len(ordered)-1)):
            raise AssertionError("Monotonicity of DM-conjugate ground chirality violated")
    for ref in records:
        n=ref["n_sites"]
        D=ref["D"]
        configurations=(1,2) if D==0. else (1,)
        for alg in configurations:
            try:
                row=run_ed2(exe,root,ref,alg,args.omp)
                run_rows.append(row)
                print("DM VALIDATED",row["case"],
                      "status",row["status"],
                      "classification",row["classification"],
                      "Eerr",row.get("eigenenergy_max_error"),
                      "chirality_err",row.get("staggered_chirality_error"),
                      "HF_err",row.get("hellmann_feynman_ed2_error"),
                      "max_r",row.get("max_direct_residual"),
                      "Hcalls",row["h_applications"],flush=True)
            except Exception as ex:
                violation=dict(n_sites=n,D=D,alg=alg,omp=args.omp,error=str(ex))
                failures.append(violation)
                print("DM CONTRACT FAILURE",violation,flush=True)
    with (root/"physical-scan.csv").open("w",newline="") as fp:
        writer=csv.DictWriter(fp,fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(summary)
    machine=dict(schema=2,source_sha=os.environ.get("GITHUB_SHA","provide from CI"),
                 profile=args.profile,
                 model="non-U1 open alternating S=1/2,1 chain with staggered oriented Dz",
                 DM_values=list(D_VALUES),sizes=sizes,omp=args.omp,
                 calibration="dE0/dD equals staggered vector chirality, full Cartesian NN correlators",
                 independent_physics=summary,
                 cases=run_rows,contract_violations=failures,
                 expected=len(sizes)*6,completed=len(run_rows),
                 successful=sum(x["status"]==0 for x in run_rows),
                 numerical_nonconvergence=sum(x["status"]!=0 for x in run_rows))
    (root/"report.json").write_text(json.dumps(machine,indent=2,allow_nan=False)+"\n")
    print(f"DM SWEEP SUMMARY: {machine['completed']}/{machine['expected']} classified; "
          f"{machine['successful']} accepted; "
          f"{machine['numerical_nonconvergence']} nonconverged; "
          f"{len(failures)} contract violations; OMP={args.omp}",flush=True)
    return 0 if (machine["completed"]==machine["expected"] and
                 machine["numerical_nonconvergence"]==0 and not failures) else 1


if __name__=="__main__":
    raise SystemExit(main())
