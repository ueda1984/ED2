#!/usr/bin/env python3
"""Basis-invariant low-doublet observable averages from native ED2 vectors.

The first and second eigenstates are obtained from ED2 ALG=2 and
their two separate sets of observables are computed by the actual
ED2 observable routines. The second native-basis Fortran wavefunction
is copied into a fresh standard first-state readback path; this
avoids assuming ED2's internal basis ordering.

Then compare the *equal-weight doublet ensemble* with an independently
constructed dense NumPy projector P=|0><0|+|1><1|.

This test does NOT treat any particular vector of an exactly
degenerate doublet as a uniquely defined physical state.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

import mixed_spin_nonu1_demo as dense
from mixed_spin_near_degenerate import reference,run_case

CASES=((4,0.0),(4,1e-10),(6,0.0),(6,1e-10))
WORKSPACES=((2,16),(4,24))
ATOL=1e-8


def expected_ensemble(model):
    spins=model["spins"]
    operators,_=dense.product_operators(spins)
    matrix=dense.dense_hamiltonian(operators,model["bonds"],model["fields"])
    values,vectors=np.linalg.eigh(matrix)
    q=vectors[:,:2]
    sites=len(spins)
    mag=np.array([[float(np.trace(q.conj().T @ operators[i][a] @ q).real)/2
                    for a in range(3)] for i in range(sites)])
    corr=np.array([[float(np.trace(q.conj().T
                      @ operators[i][a] @ operators[j][b] @ q).real)/2
                     for b in range(3) for a in range(3)]
                   for i,j in model["pairs"]])
    return mag,corr,values[:3].tolist()


def readback_second(exe,source,target,omp):
    target.mkdir(exist_ok=False)
    filenames=("bonds.dat","fields.dat","spins.dat","caps.dat","pairs.dat")
    for f in filenames:
        shutil.copy2(source/f,target/f)
    input_text=(source/"input.dat").read_text()
    for before,after in (("ALG=2","ALG=1"),
                         ("wr_wf=1, re_wf=0","wr_wf=0, re_wf=1")):
        if input_text.count(before)!=1:
            raise AssertionError("Unrecognized readback source deck")
        input_text=input_text.replace(before,after)
    (target/"input.dat").write_text(input_text)
    first=[p for p in source.iterdir() if p.is_file() and p.name.endswith("wf1.bin")]
    second=[p for p in source.iterdir() if p.is_file() and p.name.endswith("wf2.bin")]
    energy=[p for p in source.iterdir() if p.is_file() and p.name.endswith("energy.dat")]
    if len(first)!=1 or len(second)!=1 or len(energy)!=1:
        raise AssertionError("Missing native wavefunctions or energy metadata")
    # Preserve exact gfortran filename padding from FILE_wf A100 formatter.
    shutil.copy2(second[0],target/first[0].name)
    shutil.copy2(energy[0],target/energy[0].name)
    env=dict(os.environ,OMP_NUM_THREADS=str(omp),OMP_DYNAMIC="FALSE",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    for key in ("ED2_RESIDUAL_ATOL","ED2_RESIDUAL_RTOL","ED2_DIAGNOSTICS",
                "ED2_INITIAL_VECTOR","ED2_RNG_SEED","ED2_RESIDUAL_SCALE",
                "ED2_ENERGY_ORIGIN"):
        env.pop(key,None)
    with (target/"input.dat").open("rb") as inp:
        result=subprocess.run([str(exe)],cwd=target,stdin=inp,env=env,
                              stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                              timeout=80)
    (target/"stdout.log").write_bytes(result.stdout)
    if result.returncode!=0:
        raise AssertionError("Second-state observable readback failed "+
                             result.stdout.decode(errors="replace")[-1900:])
    return np.loadtxt(target/"local_mag.dat",ndmin=2),np.loadtxt(
        target/"two_body_cf_xyz.dat",ndmin=2)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--omp",type=int,required=True,choices=(1,2))
    args=p.parse_args()
    exe=args.exe.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    records=[]
    errors=[]
    for n,delta in CASES:
        model=reference(n,delta)
        expected_mag,expected_cf,values=expected_ensemble(model)
        if abs(expected_mag[0,2])>1e-7:
            raise AssertionError("Independent doublet free-spin ensemble not neutral")
        for keep,capacity in WORKSPACES:
            sub=root/f"n{n}-delta{delta:.0e}-keep{keep}-capacity{capacity}"
            sub.mkdir(exist_ok=False)
            try:
                first=run_case(exe,sub,n,delta,2,keep,capacity,ATOL,args.omp,
                               model,suffix="-doublet")
                if first["status"]!=0 or not first.get("inside_low_doublet"):
                    raise AssertionError("Requested low-doublet eigenpairs not accepted")
                source=sub/first["case"]
                m1=np.loadtxt(source/"local_mag.dat",ndmin=2)
                c1=np.loadtxt(source/"two_body_cf_xyz.dat",ndmin=2)
                m2,c2=readback_second(exe,source,sub/"second-state-readback",args.omp)
                if m1.shape!=(n,4) or m2.shape!=(n,4):
                    raise AssertionError("Local magnetization columns mismatch")
                if c1.shape!=(len(model["pairs"]),11) or c2.shape!=c1.shape:
                    raise AssertionError("Pair correlation columns mismatch")
                np.testing.assert_array_equal(m1[:,0],m2[:,0])
                np.testing.assert_array_equal(c1[:,:2],c2[:,:2])
                ensemble_m=(m1[:,1:4]+m2[:,1:4])/2
                ensemble_c=(c1[:,2:11]+c2[:,2:11])/2
                mag_err=float(np.max(np.abs(ensemble_m-expected_mag)))
                cf_err=float(np.max(np.abs(ensemble_c-expected_cf)))
                predicted_ensemble_z=float(expected_mag[0,2])
                measured_ensemble_z=float(ensemble_m[0,2])
                if mag_err>3e-6 or cf_err>3e-6:
                    raise AssertionError("Doublet-averaged observable mismatch "+
                                         str((mag_err,cf_err)))
                records.append(dict(case=f"n{n}-delta{delta:.0e}-k{keep}-m{capacity}",
                                    omp=args.omp,n_sites=n,delta=delta,dimension=model["dimension"],
                                    external_gap=model["outside_gap"],
                                    independent_low_eigs=values[:2],
                                    solver_eigs=first["returned_energies"],
                                    direct_residuals=first["direct_residuals"],
                                    numerical_cluster_leakage_bound=first["conditional_low_subspace_leakage_upper"],
                                    solver_gram_error=first["orthonormality_max_error"],
                                    individual_first_z=float(m1[0,3]),
                                    individual_second_z=float(m2[0,3]),
                                    measured_ensemble_z=measured_ensemble_z,
                                    predicted_ensemble_z=predicted_ensemble_z,
                                    ensemble_magnetization_error=mag_err,
                                    ensemble_correlations_error=cf_err,
                                    status="passed"))
                print("DOUBLET TRACE PASS",
                      records[-1]["case"],"OMP",args.omp,
                      "firstSz",records[-1]["individual_first_z"],
                      "secondSz",records[-1]["individual_second_z"],
                      "ensembleSz",measured_ensemble_z,
                      "maxMagError",mag_err,"maxCorrError",cf_err,
                      "clusterBound",records[-1]["numerical_cluster_leakage_bound"],
                      flush=True)
            except Exception as e:
                errors.append(dict(case=f"n{n}-delta{delta:.0e}-k{keep}-m{capacity}",
                                   message=str(e)))
                print("DOUBLET TRACE FAILURE",errors[-1],flush=True)
    report=dict(schema=1,source_sha=os.environ.get("GITHUB_SHA","record from CI"),
                method="native ED2 wavefunction readback; independent dense doublet projector",
                omp=args.omp,records=records,failures=errors,
                passed=len(records),total=len(records)+len(errors))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"DOUBLET TRACE SUMMARY: {len(records)}/{len(records)+len(errors)} "
          f"passed, {len(errors)} failures, OMP={args.omp}",flush=True)
    return 0 if len(records)==len(CASES)*len(WORKSPACES) and not errors else 1


if __name__=="__main__":
    raise SystemExit(main())
