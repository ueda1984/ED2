#!/usr/bin/env python3
"""Adversarial low-pair test with a genuinely noncommuting XY edge bond.

H=Hactive + delta*S1z + epsilon*(S1x*S2x+0.6*S1y*S2y).
Unlike the initial pure-X coupling, S1x and S1z BOTH fail to commute
with H even at delta=0. The ground state is mixed/entangled at the
avoided-crossing centre; the low gap can be made arbitrarily small.

Run independent NumPy dense ground and doublet-trace checks on ED2.
Also inject an accurately verified *interacting first excited* native
ED2 eigenvector as the ALG1 start and verify the documented distinction
between residual success and global lowest-state inclusion.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np

from mixed_spin_interacting_avoided_crossing import (
    physical_reference, run_case, hermitian_expectation,
)
from mixed_spin_nonu1_demo import couplings, product_operators, dense_hamiltonian
from mixed_spin_near_degenerate import load_fortran_vector

FRACTION_Y=.6
EPSILONS=(.02,5e-7)
SIZES=(4,6)
DELTAS={
    .02:(-.004,-.002,0.,.002,.004),
    5e-7:(-1e-7,0.,1e-7),
}
ABS_RESIDUAL=1e-8


def ensure_interacting(n,delta,epsilon):
    ref=physical_reference(n,delta,epsilon,edge_y_fraction=FRACTION_Y)
    ops,_=product_operators(ref["spins"])
    H=dense_hamiltonian(ops,ref["bonds"],ref["fields"])
    norm=float(np.linalg.norm(H))
    edge_comm=[float(np.linalg.norm(H@ops[0][a]-ops[0][a]@H)/norm)
               for a in (0,1,2)]
    if epsilon==.02 and min(edge_comm)<1e-3:
        raise AssertionError("Supposed noncommuting edge-spin interaction has a conserved component")
    if epsilon==.02 and delta==0 and ref["edge_linear_entropy"]<1e-4:
        raise AssertionError("Avoided-crossing centre unexpectedly factorizes")
    ref["delta"]=delta
    ref["epsilon"]=epsilon
    ref["edge_commutator_norms"]=edge_comm
    return ref


def check_geometry(models):
    for n in SIZES:
        m={(p["epsilon"],p["delta"]):p for p in models if p["n_sites"]==n}
        for epsilon in EPSILONS:
            centre=m[(epsilon,0.)]["ground_gap"]
            points=DELTAS[epsilon]
            far_negative=m[(epsilon,min(points))]
            far_positive=m[(epsilon,max(points))]
            if not (centre>5e-8 and
                    far_negative["ground_gap"]>centre and
                    far_positive["ground_gap"]>centre):
                raise AssertionError("No finite low-energy avoided crossing")
            if not far_negative["edge_sz"]>.1 or not far_positive["edge_sz"]<-.1:
                raise AssertionError("No robust edge magnetization reversal")
            if epsilon==.02 and not (centre>.003 and centre<.005):
                raise AssertionError("Visible noncommuting avoided gap was unexpected")
    return True


def estimate_shifted_crossing():
    """Bracket the avoided-gap minimum using HF derivatives and a parabola.

    This is a reproducible independent finite-size characterization,
    not an exact symbolic solution or a proof of global monotonicity.
    """
    result=[]
    for n in SIZES:
        eps=.02
        a=5e-4
        left=ensure_interacting(n,-a,eps)
        centre=ensure_interacting(n,0.,eps)
        right=ensure_interacting(n,+a,eps)
        gminus=left["ground_gap"]
        gzero=centre["ground_gap"]
        gplus=right["ground_gap"]
        denominator=2*(gminus-2*gzero+gplus)
        if denominator<=0:
            raise AssertionError("Local gap curvature is not positive")
        delta_est=a*(gminus-gplus)/denominator
        minref=ensure_interacting(n,delta_est,eps)
        # d(E1-E0)/d(delta) = <S1z>_1 - <S1z>_0.
        slope0=centre["excited_edge_sz"]-centre["edge_sz"]
        slopeplus=right["excited_edge_sz"]-right["edge_sz"]
        if not (0<delta_est<a and slope0<0 and slopeplus>0
                and minref["ground_gap"]<gzero):
            raise AssertionError("Could not establish a shifted local avoided minimum")
        item=dict(n_sites=n,epsilon=eps,jy_over_jx=FRACTION_Y,
                  left_offset=-a,right_offset=a,
                  gap_at_negative_offset=gminus,gap_at_zero=gzero,
                  gap_at_positive_offset=gplus,
                  slope_at_zero=slope0,slope_at_positive_offset=slopeplus,
                  parabolic_delta_estimate=delta_est,
                  gap_at_parabolic_estimate=minref["ground_gap"])
        result.append(item)
        print("SHIFTED CROSSING ESTIMATE",item,flush=True)
    return result


def first_order_doublet_projection(models):
    """Check weak-bond gaps against independent first-order doublet physics.

    For the decoupled active-chain ground state |g>, projecting
    epsilon*(S1x S2x + eta*S1y S2y) gives
    H_eff=e0 I+(delta*sigma_z+epsilon*mx*sigma_x+
                  epsilon*eta*my*sigma_y)/2,
    where mx=<g|S2x|g>, my=<g|S2y|g>.
    Thus the delta=0 gap is first-order
    epsilon*sqrt(mx^2+eta^2*my^2).
    """
    results=[]
    for n in SIZES:
        spins=[.5,1.]*(n//2)
        bonds,fields=couplings(n,.15)
        if bonds[0][:2]!=(0,1):
            raise AssertionError("Unexpected active bond orientation")
        reduced_spins=spins[1:]
        ob,_=product_operators(reduced_spins)
        Hact=dense_hamiltonian(ob,[
            (i-1,j-1,c) for i,j,c in bonds if i>0 and j>0
        ],fields[1:])
        E,V=np.linalg.eigh(Hact)
        q=V[:,0]
        mx=hermitian_expectation(q,ob[0][0])
        my=hermitian_expectation(q,ob[0][1])
        for eps in EPSILONS:
            g=next(x["ground_gap"] for x in models
                   if x["n_sites"]==n and x["epsilon"]==eps and x["delta"]==0.)
            perturbative=eps*np.sqrt(mx**2+(FRACTION_Y*my)**2)
            relerr=float(abs(g-perturbative)/perturbative)
            threshold=.01 if eps==.02 else 5e-5
            if relerr>threshold:
                raise AssertionError("Independent doublet projection failed "+
                                     str((n,eps,perturbative,g,relerr)))
            item=dict(n_sites=n,epsilon=eps,
                      active_spin2_mx=mx,active_spin2_my=my,
                      predicted_first_order_gap=float(perturbative),
                      observed_exact_gap=float(g),
                      relative_gap_difference=relerr,
                      active_chain_gap=float(E[1]-E[0]))
            results.append(item)
            print("FIRST ORDER DOUBLET PROJECTION",item,flush=True)
    return results


def excited_start_check(exe,root,omp,model):
    work=root/"explicit-interacting-excited"
    work.mkdir(exist_ok=False)
    source=run_case(exe,work,model,2,4,24,1e-12,omp)
    if source["status"]!=0 or not source.get("inside_low_cluster"):
        raise AssertionError("Could not independently prepare both interacting low eigenpairs")
    if any(abs(e-r)>3e-8 for e,r in
           zip(source["returned_energies"],model["low_e"][:2])):
        raise AssertionError("Source's first excited eigenvalue not independently verified")
    source_dir=work/source["case"]
    q=load_fortran_vector(source_dir,2,model["dimension"])
    trial=run_case(exe,work,model,1,0,0,1e-12,omp,initial_vector=q,
                   suffix="-explicit-excited-start")
    if trial["status"]!=0:
        raise AssertionError("Known excited eigenvector was not accepted under low residual")
    E=float(trial["returned_energies"][0])
    E0,E1=model["low_e"][:2]
    gap=float(E1-E0)
    magnetization=float(trial["measured_edge_sz"])
    expected_excited=float(model["excited_edge_sz"])
    if not (abs(E-E1)<3e-8 and E-E0>.9*gap
            and abs(magnetization-expected_excited)<1e-5):
        raise AssertionError("Explicit excited-state inclusion witness not reproduced")
    # This result is NOT a numerical solver bug; the wavefunction is
    # a mathematically correct first excited eigenstate.
    return dict(source=source,trial=trial,first_excited_witness=True,
                ground_energy=float(E0),excited_energy=float(E1),
                returned_energy=E,energy_above_ground=float(E-E0),
                direct_residual=trial["direct_residuals"][0],
                measured_edge_sz=magnetization,
                independent_excited_edge_sz=expected_excited)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--omp",required=True,type=int,choices=(1,2))
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    models=[]
    cases=[]
    failures=[]
    for n in SIZES:
        for epsilon in EPSILONS:
            for delta in DELTAS[epsilon]:
                model=ensure_interacting(n,delta,epsilon)
                models.append(model)
                print(f"XY REFERENCE n={n} epsilon={epsilon:.3g} delta={delta:+.3g} "
                      f"gap={model['ground_gap']:.10g} "
                      f"Szedge={model['edge_sz']:+.8g} "
                      f"edge-entropy={model['edge_linear_entropy']:.4g} "
                      f"edge-comms={model['edge_commutator_norms']}",flush=True)
    check_geometry(models)
    crossing_minima=estimate_shifted_crossing()
    perturbative_checks=first_order_doublet_projection(models)
    for model in models:
        runs=[(1,0,0,ABS_RESIDUAL),
              (2,2,16,ABS_RESIDUAL),
              (2,4,24,ABS_RESIDUAL)]
        if model["epsilon"]==5e-7:
            runs.append((1,0,0,1e-12))
        for alg,keep,capacity,atol in runs:
            try:
                entry=run_case(exe,root,model,alg,keep,capacity,atol,args.omp)
                cases.append(entry)
                print("XY MEASURED",entry["case"],"status",entry["status"],
                      "class",entry["classification"],
                      "edgeSz",entry.get("measured_edge_sz"),
                      "direct-r",entry.get("direct_residuals"),
                      "cluster-trace-error",entry.get("cluster_ensemble_mag_error"),
                      flush=True)
            except Exception as e:
                failures.append(dict(n_sites=model["n_sites"],
                                     epsilon=model["epsilon"],delta=model["delta"],
                                     alg=alg,keep=keep,capacity=capacity,atol=atol,error=str(e)))
                print("XY CONTRACT FAILURE",failures[-1],flush=True)
    try:
        witness_model=next(m for m in models if
                           m["n_sites"]==4 and m["epsilon"]==.02 and m["delta"]==.002)
        witness=excited_start_check(exe,root,args.omp,witness_model)
        print("INTERACTING EXCITED-START WITNESS:",
              "energy_above_ground",witness["energy_above_ground"],
              "residual",witness["direct_residual"],
              "edgeSz",witness["measured_edge_sz"],flush=True)
    except Exception as e:
        witness=None
        failures.append(dict(label="interacting-first-excited-start",error=str(e)))
        print("XY EXCITED-START FAILURE",str(e),flush=True)
    expected=sum(3+(m["epsilon"]==5e-7) for m in models)
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","see CI"),
                model="anisotropic noncommuting XY edge bond plus active XYZ/DM/Gamma",
                edge_y_over_x=FRACTION_Y,
                shifted_gap_minima=crossing_minima,
                first_order_doublet_projection=perturbative_checks,
                omp=args.omp,independent_models=[
                    dict(n_sites=m["n_sites"],epsilon=m["epsilon"],delta=m["delta"],
                         ground_e=m["ground_e"],lowest_three=m["low_e"],
                         ground_gap=m["ground_gap"],external_gap=m["outer_gap"],
                         edge_sz=m["edge_sz"],edge_entanglement_linear_entropy=m["edge_linear_entropy"],
                         edge_commutator_norms=m["edge_commutator_norms"],
                         excited_edge_sz=m["excited_edge_sz"]) for m in models],
                trials=cases,expected=expected,completed=len(cases),
                accepted=sum(c["status"]==0 for c in cases),
                nonconverged=sum(c["status"]!=0 for c in cases),
                individually_unresolved=sum(c["classification"]==
                    "accepted_low_doublet_but_individual_ground_unresolved" for c in cases),
                excited_initial_witness=witness,contract_violations=failures)
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"NONCOMMUTING XY SUMMARY: {len(cases)}/{expected} classified; "
          f"{report['accepted']} accepted, {report['nonconverged']} nonconverged, "
          f"{report['individually_unresolved']} individually unresolved; "
          f"excited-witness={witness is not None}; "
          f"{len(failures)} contract failures; OMP={args.omp}",flush=True)
    return 0 if len(cases)==expected and not failures else 1


if __name__=="__main__":
    raise SystemExit(main())
