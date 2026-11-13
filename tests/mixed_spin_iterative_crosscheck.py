#!/usr/bin/env python3
"""Independently cross-check ED2 ALG=1/2 on mixed-spin non-U(1) models.

Reuses the separately validated NumPy product-basis reference in
mixed_spin_nonu1_demo.py. The dense reference is *not* an ED2 matrix,
basis reorder, or Krylov implementation. Every solver output is checked
against an independently diagonalized Hamiltonian and its observables.

Expected solver nonconvergence is recorded as a failed scientific
test, rather than being silently skipped. Results are preserved in
new, never-reused directories.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
import mixed_spin_nonu1_demo as ref


SEED = 19
RESIDUAL_ABS = 1e-8
EIGEN_ATOL = 2e-7
OBS_ATOL = 5e-6
HF_ATOL = 4e-6
MODELS = ((4, 0.15, "full"), (4, 0.38, "full"),
          (6, 0.15, "full"), (6, 0.38, "full"),
          (4, 0.15, 2), (4, 0.38, 2))
INPUTS = ("bonds.dat", "fields.dat", "spins.dat", "caps.dat",
          "pairs.dat", "input.dat")


def read_diagnostics(path, expected_pairs):
    data = {}
    pairs = []
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            data[key.strip()] = value.strip()
        else:
            parts = line.split()
            if len(parts) == 8 and parts[0].isdigit():
                pairs.append(dict(index=int(parts[0]), available=parts[1],
                                  certified=parts[2],
                                  ritz=float(parts[3]), rayleigh=float(parts[4]),
                                  residual=float(parts[5]), norm=float(parts[6]),
                                  legacy_threshold=float(parts[7])))
    if len(pairs) != expected_pairs:
        raise AssertionError(f"Expected {expected_pairs} diagnostic pairs, got {pairs}")
    if data.get("status") != "0":
        raise AssertionError("Krylov did not converge: "+str(data.get("reason")))
    if data.get("lowest_state_inclusion") != "not_assessed":
        # Diagnostic inclusion field is advisory, not an independent proof.
        # This test provides external inclusion evidence via dense energies.
        raise AssertionError("Unexpected lowest-state inclusion contract")
    if data.get("seed_specified") != "T":
        raise AssertionError("Random seed not recorded")
    if data.get("all_requested_residuals_met") != "T":
        raise AssertionError("Residual policy not satisfied for all pairs")
    for i, pair in enumerate(pairs, 1):
        if pair["index"] != i or pair["available"] != "T" or pair["certified"] != "T":
            raise AssertionError("Invalid diagnostic eigenpair "+repr(pair))
        if pair["residual"] > RESIDUAL_ABS * (1 + 2e-14):
            raise AssertionError("Reported residual exceeds requested tolerance")
        if data.get(f"requested_residual_met_{i}") != "T":
            raise AssertionError("Per-pair residual policy not met")
    return data, pairs


def hellmann_feynman_error(n, hx, cap, recorded_magnetization):
    """Finite-difference derivative of independent *projected* dense H.

    The same fixed subspace is retained on both sides of the central
    difference. This is not an invalid finite-difference comparison
    between distinct symmetry sectors.
    """
    spins = [0.5, 1.0] * (n // 2)
    full_cap = sum(round(2*s) for s in spins)
    max_occ = full_cap if cap == "full" else int(cap)
    ops, dims = ref.product_operators(spins)
    ids = ref.subspace_indices(dims, max_occ)
    delta = 1e-4
    def exact_energy(x):
        bonds, fields = ref.couplings(n, x)
        h = ref.dense_hamiltonian(ops, bonds, fields)
        return float(np.linalg.eigvalsh(h[np.ix_(ids, ids)])[0])
    numerical = (exact_energy(hx + delta) - exact_energy(hx - delta))/(2*delta)
    weights = np.array([0.4*((-1)**i)+0.09*(i+1) for i in range(n)])
    from_ed2 = float(weights @ recorded_magnetization[:, 0])
    error = abs(numerical - from_ed2)
    if error > HF_ATOL:
        raise AssertionError("Hellmann-Feynman derivative mismatch: "
                             +str((numerical, from_ed2, error)))
    return dict(finite_difference=numerical, magnetic_response=from_ed2,
                absolute_error=error, difference_step=delta)


def run_solver(exe, outroot, baseline, n, hx, cap, alg, omp):
    name = f"mixed-n{n}-hx{hx:.2f}-cap{cap}-alg{alg}-omp{omp}"
    scratch = outroot / name
    scratch.mkdir(exist_ok=False)
    inputs = outroot / "dense-baselines" / f"mixed-n{n}-hx{hx:.2f}-cap{cap}"
    for filename in INPUTS:
        shutil.copyfile(inputs/filename, scratch/filename)

    text = (scratch/"input.dat").read_text()
    if text.count("ALG=3") != 1:
        raise AssertionError("Reference deck ALG selector ambiguous")
    text = text.replace("ALG=3",f"ALG={alg}")
    # Explicitly document a generous algorithmic budget, while the
    # independent residual and physical accuracy tests remain strict.
    text = text.replace("lnc_ene_conv=1d-14, minitr=5, maxitr=1000, itrint=1",
                        "lnc_ene_conv=1d-12, minitr=5, maxitr=400, itrint=1")
    text = text.replace("NOE=2, NOK=4, NOM=8, maxitr=1000, lnc_ene_conv=1d-14",
                        "NOE=2, NOK=8, NOM=48, maxitr=350, lnc_ene_conv=1d-12")
    (scratch/"input.dat").write_text(text)

    env = dict(os.environ,OMP_NUM_THREADS=str(omp),OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", ED2_RNG_SEED=str(SEED),
               ED2_DIAGNOSTICS="diagnostics.txt",
               ED2_RESIDUAL_ATOL=str(RESIDUAL_ABS))
    begin = time.perf_counter()
    with (scratch/"input.dat").open("rb") as inp:
        result = subprocess.run([str(exe)],cwd=scratch,stdin=inp,env=env,
                                stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                timeout=150)
    elapsed = time.perf_counter()-begin
    (scratch/"stdout.log").write_bytes(result.stdout)
    if result.returncode != 0:
        raise AssertionError(
            f"{name} returned {result.returncode}; see stdout and diagnostic:\n"+
            result.stdout.decode(errors="replace")[-3000:])

    num_pairs = 1 if alg == 1 else 2
    diag, pairs = read_diagnostics(scratch/"diagnostics.txt",num_pairs)
    if int(diag["algorithm"])!=alg or int(diag["dimension"])!=baseline["projected_dimension"]:
        raise AssertionError("Incorrect diagnostics algorithm or basis dimension")
    energies = np.loadtxt(scratch/"eigenvalues.dat",ndmin=2)
    reference = np.loadtxt(inputs/"eigenvalues.dat",ndmin=2)
    if energies.shape != (num_pairs,2):
        raise AssertionError("Unexpected eigenvalue count")
    np.testing.assert_array_equal(energies[:,0],np.arange(1,num_pairs+1))
    np.testing.assert_allclose(energies[:,1],reference[:num_pairs,1],
                               atol=EIGEN_ATOL,rtol=0)
    if np.any(np.diff(energies[:,1]) < 0):
        raise AssertionError("Eigenvalues not sorted")
    for j,p in enumerate(pairs):
        np.testing.assert_allclose(p["ritz"], energies[j,1],atol=1e-9,rtol=0)

    actual_mag = np.loadtxt(scratch/"local_mag.dat",ndmin=2)
    exact_mag = np.loadtxt(inputs/"local_mag.dat",ndmin=2)
    actual_cf = np.loadtxt(scratch/"two_body_cf_xyz.dat",ndmin=2)
    exact_cf = np.loadtxt(inputs/"two_body_cf_xyz.dat",ndmin=2)
    np.testing.assert_array_equal(actual_mag[:,0],exact_mag[:,0])
    np.testing.assert_array_equal(actual_cf[:,:2],exact_cf[:,:2])
    mag_err = float(np.max(np.abs(actual_mag[:,1:]-exact_mag[:,1:])))
    cf_err = float(np.max(np.abs(actual_cf[:,2:]-exact_cf[:,2:])))
    if mag_err>OBS_ATOL or cf_err>OBS_ATOL:
        raise AssertionError(f"Observables inaccurate despite success: {mag_err=} {cf_err=}")

    hf = hellmann_feynman_error(n,hx,cap,actual_mag[:,1:])
    return dict(case=name,alg=alg,omp=omp,n_sites=n,hx_control=hx,
                sector_cap=cap,sector_dimension=baseline["projected_dimension"],
                seed=SEED,residual_atol=RESIDUAL_ABS,
                values=energies[:,1].tolist(),reference_values=reference[:num_pairs,1].tolist(),
                energy_max_error=float(np.max(np.abs(energies[:,1]-reference[:num_pairs,1]))),
                max_residual=max(p["residual"] for p in pairs),
                magnetization_max_error=mag_err,correlation_max_error=cf_err,
                hellmann_feynman=hf,h_applications=int(diag["h_applications"]),
                iterations=int(diag["iterations"]),cycles=int(diag["cycles"]),
                reason=diag["reason"],lowest_state_inclusion_verified="compared to independent dense eigenvalues",
                wall_time_seconds=elapsed,status="passed")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=1,choices=(1,2))
    args=parser.parse_args()
    exe=args.exe.resolve(strict=True)
    output=args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    dense=output/"dense-baselines"
    dense.mkdir()
    baselines=[]
    results=[]
    for n,hx,cap in MODELS:
        baselines.append(ref.run_case(exe,dense,n,hx,cap))
    for baseline,(n,hx,cap) in zip(baselines,MODELS):
        for alg in (1,2):
            ident=f"mixed-n{n}-hx{hx:.2f}-cap{cap}-alg{alg}-omp{args.threads}"
            try:
                r=run_solver(exe,output,baseline,n,hx,cap,alg,args.threads)
                results.append(r)
                print(f"PASS {ident}: energy {r['energy_max_error']:.3e} "
                      f"r {r['max_residual']:.3e}, mag {r['magnetization_max_error']:.3e}, "
                      f"cf {r['correlation_max_error']:.3e}",flush=True)
            except Exception as e:
                results.append(dict(case=ident,status="failed",error=str(e)))
                print(f"FAIL {ident}: {e}",flush=True)
    report=dict(schema=1,kind="independent dense-to-Krylov mixed-spin verification",
                source_sha=os.environ.get("GITHUB_SHA","record externally"),
                compiler_omp_threads=args.threads,
                baseline_cases=len(baselines),checks=len(results),
                passed=sum(x["status"]=="passed" for x in results),
                failed=sum(x["status"]=="failed" for x in results),
                cases=results)
    (output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"RESULT: {report['passed']}/{len(results)} passed; "
          f"{report['failed']} failed",flush=True)
    return 0 if report["failed"]==0 and len(results)==12 else 1


if __name__=="__main__":
    raise SystemExit(main())
