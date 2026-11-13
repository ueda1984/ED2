#!/usr/bin/env python3
"""Near-degenerate physical Hamiltonian stress test for ED2 ALG=1/2.

Scientific question: does the returned state lie in the lowest doublet,
and when can its local observables identify the *individual* ground state?

A physically isolated S=1/2 edge spin produces an analytically controlled
doublet splitting abs(delta), while the remaining XYZ/DM/Gamma chain is
non-U(1). Independent NumPy full dense diagonalization verifies the
entire Hamiltonian, spectrum and isolated-spin observables.

Correctly signalled nonconvergence and accurately computed but
individually mixed low-doublet states are CHARACTERIZED, not hidden.
The harness fails for broken I/O, residual contracts or reference
invariants. A status-zero eigenpair is not presumed to be the global
ground state. All runs retain inputs and diagnostic evidence.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from mixed_spin_nonu1_demo import (
    couplings, product_operators, dense_hamiltonian,
)

SEED=19
# Exact spectral doublet, not a machine-generated nearly degenerate QR matrix.
SCENARIOS=((4,1e-2),(4,1e-6),(4,-1e-6),(4,1e-10),(4,0.0),
           (6,1e-6),(6,1e-10),(6,0.0))
BASIC=((1,0,0,1e-8),(2,2,16,1e-8),(2,4,24,1e-8))
TIGHT=((1,0,0,1e-12),(2,4,24,1e-12))


def reference(n,delta):
    spins=[0.5,1.0]*(n//2)
    bonds,fields=couplings(n,0.15)
    # No coupling to site 0, and its only field is +delta Sz_0.
    bonds=[(i,j,c) for i,j,c in bonds if i>0 and j>0]
    fields[0]=[0.0,0.0,delta]
    ops,dims=product_operators(spins)
    h=dense_hamiltonian(ops,bonds,fields)
    es,vs=np.linalg.eigh(h)
    # A second, explicit tensor-factor formula checks the isolation logic.
    active_ops,_=product_operators(spins[1:])
    active_bonds=[(i-1,j-1,c) for i,j,c in bonds]
    h_active=dense_hamiltonian(active_ops,active_bonds,fields[1:])
    hs=np.kron(np.eye(2),h_active)+np.kron(
        np.diag([delta/2,-delta/2]),np.eye(h_active.shape[0]))
    matrix_error=float(np.max(np.abs(h-hs)))
    if matrix_error>2e-13:
        raise AssertionError(f"Independent tensor-factor Hamiltonian mismatch {matrix_error}")
    low_active=float(np.linalg.eigvalsh(h_active)[0])
    expected_low=low_active-abs(delta)/2
    if abs(es[0]-expected_low)>2e-12:
        raise AssertionError("Independent ground energy violates tensor-factor identity")
    gap=float(es[1]-es[0])
    outside=float(es[2]-es[1])
    if abs(gap-abs(delta))>4e-13 or outside<0.15:
        raise AssertionError(f"Doublet and external gap incorrect: {delta=} {gap=} {outside=}")
    # The core is non-U(1), even though spin 0 is exactly decoupled.
    total_sz=sum((ops[i][2] for i in range(n)),start=np.zeros_like(h))
    comm_norm=float(np.linalg.norm(h@total_sz-total_sz@h)/np.linalg.norm(h))
    if comm_norm < 0.01:
        raise AssertionError("Active chain unexpectedly has U(1) symmetry")
    q=vs[:,:2]
    ensemble_z=float(np.trace(q.conj().T@ops[0][2]@q).real/2)
    if abs(ensemble_z)>1e-7:
        raise AssertionError("Low-doublet ensemble must have zero free-spin magnetization")
    return dict(spins=spins,bonds=bonds,fields=fields,dimension=h.shape[0],
                pairs=[(i,j) for i in range(n) for j in range(i+1,n)],
                eigs=es[:4].tolist(),gap=gap,outside_gap=outside,
                free_spin_expected=(-0.5 if delta>0 else 0.5 if delta<0 else None),
                free_spin_ensemble_z=ensemble_z,
                matrix_error=matrix_error,
                commutator_relative=comm_norm)


def label_delta(delta):
    return "zero" if delta==0 else f"{delta:+.0e}".replace("+","p").replace("-","m")


def read_diag(path):
    vals={}
    pairs=[]
    for line in path.read_text().splitlines():
        if "=" in line:
            k,v=line.split("=",1)
            vals[k.strip()]=v.strip()
        else:
            cols=line.split()
            if len(cols)==8 and cols[0].isdigit():
                pairs.append(dict(index=int(cols[0]),available=cols[1],
                                  certified=cols[2],ritz=float(cols[3]),
                                  rayleigh=float(cols[4]),residual=float(cols[5]),
                                  norm=float(cols[6]),legacy_threshold=float(cols[7])))
    if "status" not in vals or "h_applications" not in vals:
        raise AssertionError("Missing Krylov status/metrics")
    return vals,pairs


def load_fortran_vector(work,index,n):
    # GNU gfortran sequential unformatted: one complex128 array record.
    # The source's legacy A100 filename format may insert padding spaces,
    # so locate the exact wfN.bin suffix rather than assume a literal name.
    candidates=[p for p in work.iterdir() if p.is_file() and p.name.endswith(f"wf{index}.bin")]
    if len(candidates)!=1:
        raise AssertionError(f"Expected one wavefunction #{index}, got: {[p.name for p in candidates]}")
    raw=candidates[0].read_bytes()
    expected=n*16
    for marker_len in (4,8):
        if len(raw)!=expected+2*marker_len:
            continue
        before=int.from_bytes(raw[:marker_len],"little",signed=True)
        after=int.from_bytes(raw[-marker_len:],"little",signed=True)
        if before!=expected or after!=expected:
            continue
        vec=np.frombuffer(raw[marker_len:-marker_len],dtype="<c16").copy()
        norm=float(np.linalg.norm(vec))
        if not np.isfinite(norm) or abs(norm-1)>3e-7:
            raise AssertionError(f"Non-unit eigenvector norm: {norm}")
        return vec
    raise AssertionError("Unrecognized Fortran sequential-complex record; "
                         f"file bytes={len(raw)}, expected payload={expected}")


def run_case(exe,root,n,delta,alg,keep,capacity,atol,omp,ref,initial_vector=None,suffix=''):
    name=f"n{n}-d{label_delta(delta)}-alg{alg}-k{keep}-m{capacity}-a{atol:.0e}-omp{omp}{suffix}"
    work=root/name
    work.mkdir(exist_ok=False)
    bonds,fields,spins,pairs=(ref[k] for k in ("bonds","fields","spins","pairs"))
    (work/"bonds.dat").write_text(
        "".join(f"{i+1} {j+1} "+" ".join(f"{x:.17g}" for x in c)+"\n"
                for i,j,c in bonds))
    (work/"fields.dat").write_text(
        "".join(f"{i+1} "+" ".join(f"{x:.17g}" for x in vec)+"\n"
                for i,vec in enumerate(fields)))
    (work/"spins.dat").write_text("".join(f"{x:.1f}\n" for x in spins))
    (work/"caps.dat").write_text("".join(f"{round(2*x)}\n" for x in spins))
    (work/"pairs.dat").write_text(
        "".join(f"{i+1} {j+1}\n" for i,j in pairs))
    cap=sum(round(2*x) for x in spins)
    noe=1 if alg==1 else 2
    deck=f"""&input_parameters
NOS={n}, NODmin=0, NODmax={cap}, ALG={alg},
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={n}, NO_two={len(bonds)}, wk_dim=-1, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=1, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={len(pairs)}, FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-11, minitr=5, maxitr=500, itrint=2
/
&input_TRLan
NOE=2, NOK={keep if alg==2 else 4}, NOM={capacity if alg==2 else 24},
maxitr=180, lnc_ene_conv=1d-11, i_vec_min=1, i_vec_max=2
/
"""
    (work/"input.dat").write_text(deck)
    if initial_vector is not None:
        initial_vector=np.asarray(initial_vector,dtype=complex)
        if initial_vector.shape!=(ref["dimension"],):
            raise AssertionError("Wrong explicit initial-vector shape")
        (work/"initial-vector.txt").write_text(
            "".join(f"{z.real:.17g} {z.imag:.17g}\n" for z in initial_vector))
    env=dict(os.environ,OMP_NUM_THREADS=str(omp),OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             ED2_RNG_SEED=str(SEED),ED2_DIAGNOSTICS="diagnostics.txt",
             ED2_RESIDUAL_ATOL=str(atol))
    if initial_vector is not None:
        env["ED2_INITIAL_VECTOR"]=str(work/"initial-vector.txt")
    start=time.perf_counter()
    with (work/"input.dat").open("rb") as inp:
        proc=subprocess.run([str(exe)],cwd=work,stdin=inp,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=180)
    elapsed=time.perf_counter()-start
    (work/"stdout.log").write_bytes(proc.stdout)
    diag_path=work/"diagnostics.txt"
    if not diag_path.is_file():
        raise AssertionError("Missing solver diagnostic sidecar: "+
                             proc.stdout.decode(errors="replace")[-1200:])
    d,ps=read_diag(diag_path)
    status=int(d["status"])
    if (status==0)!=(proc.returncode==0):
        raise AssertionError("Inconsistent solver status and process exit")
    if int(d["wanted"])!=noe or int(d["dimension"])!=ref["dimension"]:
        raise AssertionError("Unexpected requested eigenpair count or basis dimension")
    if d.get("seed_specified")!="T" or d.get("residual_policy_enabled")!="T":
        raise AssertionError("Seed/residual policy is not recorded")
    result=dict(case=name,site_count=n,delta=delta,split=ref["gap"],
                gap_to_external_states=ref["outside_gap"],alg=alg,omp=omp,
                wanted=noe,keep=keep if alg==2 else None,
                capacity=capacity if alg==2 else None,
                requested_atol=atol,seed=SEED,status=status,
                reason=d["reason"],termination=d.get("termination"),
                inclusion_report=d.get("lowest_state_inclusion"),
                h_applications=int(d["h_applications"]),
                cycles=int(d["cycles"]),wall_seconds=elapsed)
    if status!=0:
        if status not in (1,2):
            raise AssertionError("Unexpected solver status outside tested nonconvergence contract")
        if any((work/p).exists() for p in ("eigenvalues.dat","local_mag.dat","two_body_cf_xyz.dat")):
            raise AssertionError("Nonconvergent solver wrote accepted physical result")
        result["classification"]="reported_nonconvergence_or_breakdown"
        return result
    if len(ps)!=noe or any(p["index"]!=j+1 for j,p in enumerate(ps)):
        raise AssertionError("Successful result lacks requested diagnostic eigenpairs")
    if any(p["available"]!="T" or p["certified"]!="T" or
           p["residual"]>atol*(1+1e-12) for p in ps):
        raise AssertionError("Accepted eigenpair violates requested direct residual")
    if d.get("all_requested_residuals_met")!="T":
        raise AssertionError("Residual gate inconsistent with successful status")
    evals=np.loadtxt(work/"eigenvalues.dat",ndmin=2)
    if evals.shape!=(noe,2):
        raise AssertionError(f"Expected {noe} output eigenvalues, got {evals.shape}")
    np.testing.assert_array_equal(evals[:,0],np.arange(1,noe+1))
    values=evals[:,1]
    vectors=[load_fortran_vector(work,j+1,ref["dimension"]) for j in range(noe)]
    gram=np.array(vectors).conj() @ np.array(vectors).T
    gram_error=float(np.max(np.abs(gram-np.eye(noe))))
    if gram_error>2e-6:
        raise AssertionError("Accepted eigenvectors lost orthonormality")
    if any(abs(p["ritz"]-v)>5e-8 for p,v in zip(ps,values)):
        raise AssertionError("Ritz and output energies disagree")
    mag=np.loadtxt(work/"local_mag.dat",ndmin=2)
    if mag.shape!=(n,4):
        raise AssertionError("Wrong local magnetization output shape")
    np.testing.assert_array_equal(mag[:,0],np.arange(1,n+1))
    m0z=float(mag[0,3])
    if abs(m0z)>0.5000001:
        raise AssertionError("Spin-half magnetization outside physical bound")
    # This threshold distinguishes the isolated lowest doublet from higher
    # active-chain excitations. It does NOT resolve individual states.
    window=min(.05,.20*ref["outside_gap"])
    in_low_cluster=bool(all(abs(float(v)-ref["eigs"][0])<window for v in values))
    error_to_each=[float(abs(float(v)-float(x))) for v,x in zip(values,ref["eigs"])]
    expected=ref["free_spin_expected"]
    mag_err=(None if expected is None else abs(m0z-expected))
    # For exact degeneracy a state-dependent magnetic moment is not unique.
    if not in_low_cluster:
        category="accepted_eigenpair_outside_low_doublet"
    elif delta==0:
        category="accepted_low_doublet_individual_vector_not_unique"
    elif mag_err<.05:
        category="accepted_low_doublet_and_correct_isolated_spin_observable"
    else:
        category="accepted_low_doublet_but_ground_observable_differs"
    rmax=max(p["residual"] for p in ps)
    # Conditional numerical bound using the independent third eigenvalue.
    # For any eigenvector with Rayleigh energy lambda below E2,
    # ||(1-P_low) v|| <= residual / (E2-lambda).
    denom=float(ref["eigs"][2]-max(values))
    cluster_leak_bound=(float(np.sqrt(sum(p["residual"]**2 for p in ps))/denom)
                        if denom>0 else None)
    gap_specific_ratio=(None if delta==0 else float(rmax/abs(delta)))
    result.update(classification=category,returned_energies=values.tolist(),
                  reference_low_energies=ref["eigs"][:noe],
                  error_to_individual_low_levels=error_to_each,
                  isolated_spin_z=m0z,isolated_spin_ground_expected=expected,
                  isolated_spin_error=mag_err,spin0_expectation_within_bound=True,
                  inside_low_doublet=in_low_cluster,
                  orthonormality_max_error=gram_error,
                  direct_residuals=[float(p["residual"]) for p in ps],
                  residual_to_lowest_gap_ratio=gap_specific_ratio,
                  conditional_low_subspace_leakage_upper=cluster_leak_bound)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--omp",type=int,default=1,choices=(1,2))
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    cases=[]
    deviations=[]
    refs=[]
    controls=[]
    for n,delta in SCENARIOS:
        ref=reference(n,delta)
        refs.append(dict(n_sites=n,delta=delta,dimension=ref["dimension"],
                         lowest_four_eigenvalues=ref["eigs"],
                         gap=ref["gap"],outside_gap=ref["outside_gap"],
                         factorization_error=ref["matrix_error"],
                         nonu1_commutator=ref["commutator_relative"],
                         ensemble_free_spin_z=ref["free_spin_ensemble_z"]))
        configs=list(BASIC)
        if delta==0 or abs(delta)<=1e-10:
            configs+=list(TIGHT)
        for alg,keep,capacity,atol in configs:
            try:
                case=run_case(exe,root,n,delta,alg,keep,capacity,atol,args.omp,ref)
                cases.append(case)
                print("MEASURED",case["case"],"status",case["status"],
                      "class",case["classification"],
                      "residuals",case.get("direct_residuals"),
                      "freeSz",case.get("isolated_spin_z"),
                      "gap",case["split"],"Hcalls",case["h_applications"],
                      flush=True)
                if delta==.01 and alg==2 and keep==4:
                    controls.append(case["status"]==0 and case.get("inside_low_doublet"))
            except Exception as error:
                deviations.append(dict(n_sites=n,delta=delta,alg=alg,keep=keep,
                                       capacity=capacity,atol=atol,error=str(error)))
                print("HARNESS VIOLATION",deviations[-1],flush=True)
    summary=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record in workflow"),
                 meaning="controlled physical low-doublet near-degeneracy, characterization not universal guarantee",
                 omp=args.omp,source_model="mixed-spin XYZ/DM/Gamma active chain plus isolated spin-half",
                 reference_cases=refs,results=cases,
                 accepted=sum(c["status"]==0 for c in cases),
                 nonconverged=sum(c["status"]!=0 for c in cases),
                 accepted_ground_magnetization_mismatch=sum(
                     c["classification"]=="accepted_low_doublet_but_ground_observable_differs"
                     for c in cases),
                 harness_violations=deviations,
                 tested=len(cases))
    (root/"report.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    expected=sum(len(BASIC)+(len(TIGHT) if d==0 or abs(d)<=1e-10 else 0)
                 for _,d in SCENARIOS)
    print(f"NEAR-DEGENERATE SUMMARY: {len(cases)}/{expected} classified; "
          f"{summary['accepted']} accepted, {summary['nonconverged']} nonconverged, "
          f"{summary['accepted_ground_magnetization_mismatch']} accepted-observable mismatches, "
          f"{len(deviations)} harness violations; OMP={args.omp}",flush=True)
    return 0 if len(cases)==expected and all(controls) and not deviations else 1


if __name__=="__main__":
    raise SystemExit(main())
