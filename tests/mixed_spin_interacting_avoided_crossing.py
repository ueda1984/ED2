#!/usr/bin/env python3
"""Independent interacting avoided-crossing validation for ED2.

An edge S=1/2 couples by epsilon*S1x*S2x to an active mixed-spin
XYZ/DM/Gamma non-U(1) chain and has a tunable delta*S1z field.
The avoided crossing is real (nonzero finite gap) rather than
two exactly uncoupled copies. Away from delta=0 the edge spin
is measurably entangled with the active subsystem.

This is a scientific characterization, not a claim that accepted
residuals globally certify the lowest state. A near-degenerate
accepted eigenpair with wrong *individual* observables is reported
as such. A two-eigenpair low-cluster projector/trace comparison is
reported and checked independently even when individual states
are ill-conditioned.

Dependencies: NumPy for research tests only, not ED2 runtime.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from mixed_spin_nonu1_demo import (
    couplings, dense_hamiltonian, product_operators,
)
from mixed_spin_near_degenerate import (
    read_diag, load_fortran_vector,
)
from mixed_spin_doublet_trace import readback_second

COUPLINGS=(.02, 5e-7)
FIELDS={
    .02:(-.004,-.002,0.0,.002,.004),
    5e-7:(-1e-7,0.0,1e-7),
}
SIZES=(4,6)
SEED=19
RES_ATOL=1e-8
ABS_ENERGY_TOL=3e-6
OBS_TOL=2e-4
ENSEMBLE_TOL=1e-5


def hermitian_expectation(psi,op):
    value=np.vdot(psi,op@psi)
    if abs(value.imag)>1e-10:
        raise AssertionError("Unexpected complex Hermitian expectation")
    return float(value.real)


def physical_reference(n,delta,epsilon):
    """Full dense matrix independently constructed from physical spin operators."""
    spins=[.5,1.]*(n//2)
    bonds,fields=couplings(n,.15)
    if len(bonds)!=n-1 or bonds[0][:2]!=(0,1):
        raise AssertionError("Unexpected mixed-spin model connection list")
    # Make the edge weakly interacting, with NO artificial local S1x field.
    # All original active-core XYZ/DM/Gamma bonds are preserved.
    bonds[0]=(0,1,[epsilon,0.,0.,0.,0.,0.,0.,0.,0.])
    fields[0]=[0.,0.,delta]
    ops,dims=product_operators(spins)
    H=dense_hamiltonian(ops,bonds,fields)
    if not np.isfinite(H).all():
        raise AssertionError("Nonfinite dense Hamiltonian")
    vals,vecs=np.linalg.eigh(H)
    low=vecs[:,:2]
    psi=low[:,0]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    magnetizations=np.array([
        [hermitian_expectation(psi,ops[i][a]) for a in range(3)]
        for i in range(n)
    ])
    correlations=np.array([
        [hermitian_expectation(psi,ops[i][a]@ops[j][b])
         for b in range(3) for a in range(3)] for i,j in pairs
    ])
    trace_mag=np.array([
        [float(np.trace(low.conj().T@ops[i][a]@low).real/2)
         for a in range(3)] for i in range(n)
    ])
    trace_cf=np.array([
        [float(np.trace(low.conj().T@ops[i][a]@ops[j][b]@low).real/2)
         for b in range(3) for a in range(3)] for i,j in pairs
    ])
    rho=psi.reshape(2,-1)@psi.reshape(2,-1).conj().T
    impurity=float(1-np.trace(rho@rho).real)
    xm=magnetizations[0,0]
    x_core=magnetizations[1,0]
    xx=correlations[pairs.index((0,1)),0]
    connected_xx=float(xx-xm*x_core)
    spinz=sum((site[2] for site in ops),np.zeros_like(H))
    comm_relative=float(np.linalg.norm(H@spinz-spinz@H)/np.linalg.norm(H))
    if comm_relative<.01:
        raise AssertionError("Active chain unexpectedly retains U(1)")
    ground_gap=float(vals[1]-vals[0])
    outer_gap=float(vals[2]-vals[1])
    if ground_gap<=1e-10 or outer_gap<.2:
        raise AssertionError("Wrong physical low-doublet gap construction")
    if epsilon==.02 and delta!=0:
        if impurity<1e-5 or abs(connected_xx)<1e-4:
            raise AssertionError("Weak bond failed to entangle edge with active core")
    if abs(delta)<1e-20 and abs(magnetizations[0,2])>5e-6:
        raise AssertionError("Zero-bias avoided crossing should have zero edge Sz")
    return dict(
        n_sites=n,spins=spins,bonds=bonds,fields=fields,pairs=pairs,
        dimension=H.shape[0],ground_e=float(vals[0]),low_e=vals[:3].tolist(),
        ground_gap=ground_gap,outer_gap=outer_gap,
        mag=magnetizations,corr=correlations,
        trace_mag=trace_mag,trace_corr=trace_cf,
        edge_sz=float(magnetizations[0,2]),
        edge_linear_entropy=impurity,
        connected_edge_xx=connected_xx,
        commutator_relative=comm_relative,
    )


def check_crossing_geometry(references):
    """Check the two branches approach and repel, rather than cross."""
    for n in SIZES:
        refs={(r["epsilon"],r["delta"]):r for r in references if r["n_sites"]==n}
        for epsilon in COUPLINGS:
            delta_values=FIELDS[epsilon]
            center=refs[(epsilon,0.0)]["ground_gap"]
            for d in delta_values:
                r=refs[(epsilon,d)]
                if d!=0 and not r["ground_gap"]>center*1.0001:
                    raise AssertionError("No avoided minimum at zero bias")
                if d>0 and r["edge_sz"]>-.01:
                    raise AssertionError("Edge Sz did not reverse for positive bias")
                if d<0 and r["edge_sz"]<+.01:
                    raise AssertionError("Edge Sz did not reverse for negative bias")
            if epsilon==.02 and not center>2e-3:
                raise AssertionError("Visible avoided crossing gap unexpectedly small")
            if epsilon==5e-7 and not 5e-8<center<2e-7:
                raise AssertionError("Tiny interacting gap not close to planned scale")
    return True


def input_deck(n,model,alg,keep,capacity):
    bonds,fields,spins,pairs=(model[k] for k in ("bonds","fields","spins","pairs"))
    cap=sum(round(2*x) for x in spins)
    return f"""&input_parameters
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
lnc_ene_conv=1d-11, minitr=5, maxitr=600, itrint=2
/
&input_TRLan
NOE=2, NOK={keep}, NOM={capacity}, maxitr=220, lnc_ene_conv=1d-11,
i_vec_min=1, i_vec_max=2
/
"""


def run_case(exe,root,model,alg,keep,capacity,atol,omp):
    n=model["n_sites"]
    epsilon=model["epsilon"]
    delta=model["delta"]
    mark=lambda x: ("zero" if x==0 else f"{x:+.0e}".replace("+","p").replace("-","m"))
    tag=(f"n{n}-epsilon{mark(epsilon)}-delta{mark(delta)}-"
         f"alg{alg}-k{keep}-m{capacity}-a{mark(atol)}-omp{omp}")
    work=root/tag
    work.mkdir(exist_ok=False)
    (work/"bonds.dat").write_text(
        "".join(f"{i+1} {j+1} "+" ".join(f"{x:.17g}" for x in c)+"\n"
                for i,j,c in model["bonds"]))
    (work/"fields.dat").write_text(
        "".join(f"{i+1} "+" ".join(f"{x:.17g}" for x in v)+"\n"
                for i,v in enumerate(model["fields"])))
    (work/"spins.dat").write_text("".join(f"{s:.1f}\n" for s in model["spins"]))
    (work/"caps.dat").write_text("".join(f"{round(2*s)}\n" for s in model["spins"]))
    (work/"pairs.dat").write_text(
        "".join(f"{i+1} {j+1}\n" for i,j in model["pairs"]))
    (work/"input.dat").write_text(input_deck(n,model,alg,keep,capacity))
    noe=1 if alg==1 else 2
    env=dict(os.environ,OMP_NUM_THREADS=str(omp),OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             ED2_RNG_SEED=str(SEED),ED2_DIAGNOSTICS="diagnostics.txt",
             ED2_RESIDUAL_ATOL=str(atol))
    started=time.perf_counter()
    with (work/"input.dat").open("rb") as inp:
        proc=subprocess.run([str(exe)],cwd=work,stdin=inp,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                            timeout=160)
    elapsed=time.perf_counter()-started
    (work/"stdout.log").write_bytes(proc.stdout)
    diag_file=work/"diagnostics.txt"
    if not diag_file.exists():
        raise AssertionError("Solver did not write diagnostics: "+
                             proc.stdout.decode(errors="replace")[-1800:])
    diag,pairs=read_diag(diag_file)
    status=int(diag["status"])
    if (status==0)!=(proc.returncode==0):
        raise AssertionError("Status disagrees with ED2 process exit")
    if int(diag["wanted"])!=noe or int(diag["dimension"])!=model["dimension"]:
        raise AssertionError("Wrong basis dimension or requested pair count")
    if diag.get("seed_specified")!="T" or diag.get("residual_policy_enabled")!="T":
        raise AssertionError("Missing solver provenance or residual policy")
    if alg==2 and (int(diag["keep"])!=keep or int(diag["capacity"])!=capacity):
        raise AssertionError("Solver workspace was changed unexpectedly")
    result=dict(case=tag,n_sites=n,epsilon=epsilon,delta=delta,
                reference_low_e=model["low_e"],reference_ground_gap=model["ground_gap"],
                reference_outer_gap=model["outer_gap"],
                independent_edge_sz=model["edge_sz"],
                independent_edge_linear_entropy=model["edge_linear_entropy"],
                independent_connected_edge_xx=model["connected_edge_xx"],
                alg=alg,keep=keep if alg==2 else None,
                capacity=capacity if alg==2 else None,requested_atol=atol,omp=omp,
                seed=SEED,status=status,reason=diag["reason"],
                h_applications=int(diag["h_applications"]),
                cycles=int(diag["cycles"]),wall_seconds=elapsed,
                lowest_inclusion_diagnostic=diag.get("lowest_state_inclusion"))
    if status!=0:
        if status not in (1,2):
            raise AssertionError("Unexpected solver status")
        if any((work/f).exists() for f in
               ("eigenvalues.dat","local_mag.dat","two_body_cf_xyz.dat")):
            raise AssertionError("Nonconvergent solver emitted accepted outputs")
        result["classification"]="reported_nonconvergence_or_breakdown"
        return result

    if len(pairs)!=noe or any(p["index"]!=j+1 for j,p in enumerate(pairs)):
        raise AssertionError("Missing successful solver eigenpair diagnostics")
    if any(p["available"]!="T" or p["certified"]!="T" or
           p["residual"]>atol*(1+2e-12) for p in pairs):
        raise AssertionError("Accepted residual exceeds opt-in direct criterion")
    if diag.get("all_requested_residuals_met")!="T":
        raise AssertionError("Accepted status without residual-contract satisfaction")
    energies=np.loadtxt(work/"eigenvalues.dat",ndmin=2)
    if energies.shape!=(noe,2):
        raise AssertionError("Unexpected number of output eigenpairs")
    np.testing.assert_array_equal(energies[:,0],np.arange(1,noe+1))
    E=energies[:,1]
    if any(E[j]>E[j+1]+2e-9 for j in range(len(E)-1)):
        raise AssertionError("Accepted eigenvalues are not sorted")
    if any(abs(E[j]-pairs[j]["ritz"])>2e-8 for j in range(noe)):
        raise AssertionError("Output eigenvalues disagree with Ritz diagnostics")
    wavefunctions=[load_fortran_vector(work,j+1,model["dimension"])
                   for j in range(noe)]
    mat=np.stack(wavefunctions,axis=1)
    gram_error=float(np.max(np.abs(mat.conj().T@mat-np.eye(noe))))
    if gram_error>2e-5:
        raise AssertionError("Accepted eigenvectors not approximately orthonormal")
    local_mag=np.loadtxt(work/"local_mag.dat",ndmin=2)
    correlations=np.loadtxt(work/"two_body_cf_xyz.dat",ndmin=2)
    if local_mag.shape!=(n,4) or correlations.shape!=(len(model["pairs"]),11):
        raise AssertionError("Physical observable file shape mismatch")
    np.testing.assert_array_equal(local_mag[:,0],np.arange(1,n+1))
    np.testing.assert_array_equal(correlations[:,:2],
                                  np.array([(i+1,j+1) for i,j in model["pairs"]]))
    mag_err=float(np.max(np.abs(local_mag[:,1:4]-model["mag"])))
    corr_err=float(np.max(np.abs(correlations[:,2:11]-model["corr"])))
    energy_err=float(abs(E[0]-model["ground_e"]))
    direct_max=float(max(p["residual"] for p in pairs))
    low_window=min(.05,.2*model["outer_gap"])
    inside_low_cluster=bool(all(abs(float(v)-model["ground_e"])<low_window for v in E))
    ground_correct=(energy_err<ABS_ENERGY_TOL and
                    mag_err<OBS_TOL and corr_err<OBS_TOL)
    if not inside_low_cluster:
        category="accepted_eigenpair_outside_low_doublet"
    elif ground_correct:
        category="accepted_correct_individual_ground"
    else:
        category="accepted_low_doublet_but_individual_ground_unresolved"
    # For clearly split models, a residual of 1e-8 must produce an
    # accurate individual ground vector; fail that test separately
    # rather than misclassifying it as a benign near-gap effect.
    if epsilon==.02 and not ground_correct:
        raise AssertionError("Accepted state misses well-separated physical ground: "+
                             str((energy_err,mag_err,corr_err)))
    result.update(classification=category,inside_low_cluster=inside_low_cluster,
                  returned_energies=E.tolist(),direct_residuals=[
                      float(p["residual"]) for p in pairs],
                  energy_ground_error=energy_err,local_magnetization_error=mag_err,
                  two_body_correlation_error=corr_err,
                  measured_edge_sz=float(local_mag[0,3]),
                  edge_sz_error=float(abs(local_mag[0,3]-model["edge_sz"])),
                  low_gap_to_residual_ratio=(
                      float(direct_max/model["ground_gap"])),
                  gram_error=gram_error)
    if alg==2 and inside_low_cluster:
        m2,c2=readback_second(exe,work,work/"second-observable",omp)
        np.testing.assert_array_equal(local_mag[:,0],m2[:,0])
        np.testing.assert_array_equal(correlations[:,:2],c2[:,:2])
        trace_mag=(local_mag[:,1:4]+m2[:,1:4])/2
        trace_corr=(correlations[:,2:11]+c2[:,2:11])/2
        tmag=float(np.max(np.abs(trace_mag-model["trace_mag"])))
        tcorr=float(np.max(np.abs(trace_corr-model["trace_corr"])))
        if tmag>ENSEMBLE_TOL or tcorr>ENSEMBLE_TOL:
            raise AssertionError("Low-doublet ensemble observables inaccurate: "+
                                 str((tmag,tcorr)))
        result.update(cluster_ensemble_mag_error=tmag,
                      cluster_ensemble_cf_error=tcorr,
                      cluster_ensemble_edge_sz=float(trace_mag[0,2]))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--omp",type=int,required=True,choices=(1,2))
    args=parser.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    tested=[]
    models=[]
    for n in SIZES:
        for epsilon in COUPLINGS:
            for delta in FIELDS[epsilon]:
                m=physical_reference(n,delta,epsilon)
                m["epsilon"]=epsilon
                m["delta"]=delta
                models.append(m)
                print(f"REFERENCE n{n} eps={epsilon:.2g} delta={delta:+.2g} "
                      f"low-gap={m['ground_gap']:.10g} "
                      f"external-gap={m['outer_gap']:.6g} "
                      f"Szedge={m['edge_sz']:+.7g} "
                      f"linear-entropy={m['edge_linear_entropy']:.5g}",flush=True)
    check_crossing_geometry(models)
    failures=[]
    for model in models:
        jobs=[(1,0,0,RES_ATOL),(2,2,16,RES_ATOL),(2,4,24,RES_ATOL)]
        if model["epsilon"]==5e-7:
            jobs.append((1,0,0,1e-12))
        for alg,keep,capacity,atol in jobs:
            try:
                r=run_case(exe,root,model,alg,keep,capacity,atol,args.omp)
                tested.append(r)
                print("MEASURED",r["case"],"status",r["status"],
                      "classification",r["classification"],
                      "residual",r.get("direct_residuals"),
                      "groundSz",r.get("measured_edge_sz"),
                      "independentSz",r["independent_edge_sz"],
                      "magerr",r.get("local_magnetization_error"),
                      "traceerr",r.get("cluster_ensemble_mag_error"),
                      "Hcalls",r["h_applications"],flush=True)
            except Exception as ex:
                failures.append(dict(
                    site_count=model["n_sites"],epsilon=model["epsilon"],
                    delta=model["delta"],alg=alg,keep=keep,
                    capacity=capacity,atol=atol,error=str(ex)
                ))
                print("CONTRACT VIOLATION",failures[-1],flush=True)
    physical=[dict(n_sites=m["n_sites"],epsilon=m["epsilon"],
                   delta=m["delta"],dimension=m["dimension"],
                   low_energies=m["low_e"],ground_gap=m["ground_gap"],
                   outer_gap=m["outer_gap"],edge_sz=m["edge_sz"],
                   edge_linear_entropy=m["edge_linear_entropy"],
                   connected_edge_xx=m["connected_edge_xx"],
                   nonu1_commutator=m["commutator_relative"]) for m in models]
    count_expected=sum(3+(m["epsilon"]==5e-7) for m in models)
    report=dict(schema=1,
                source_sha=os.environ.get("GITHUB_SHA","record in CI"),
                kind="interacting avoided crossing, gap/ground observables and low-cluster traces",
                omp=args.omp,independent_models=physical,
                trials=tested,violations=failures,expected=count_expected,
                completed=len(tested),
                accepted=sum(t["status"]==0 for t in tested),
                nonconverged=sum(t["status"]!=0 for t in tested),
                individually_unresolved=sum(
                    t["classification"]=="accepted_low_doublet_but_individual_ground_unresolved"
                    for t in tested))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"AVOIDED-CROSSING SUMMARY: {len(tested)}/{count_expected} classified; "
          f"{report['accepted']} accepted; {report['nonconverged']} nonconverged; "
          f"{report['individually_unresolved']} individually unresolved; "
          f"{len(failures)} harness failures; OMP={args.omp}",flush=True)
    return 0 if len(tested)==count_expected and not failures else 1


if __name__=="__main__":
    raise SystemExit(main())
