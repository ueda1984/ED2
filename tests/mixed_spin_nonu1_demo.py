#!/usr/bin/env python3
"""Independent dense verification of mixed-spin, non-U(1) ED2 Hamiltonians.

This research demonstrator uses only Python/NumPy for testing; it is not
a new dependency of ED2. It deliberately constructs the full product-space
spin operators without calling any ED2 code. Run only against the version
of ED2 whose input field and oriented-bond conventions are documented.

Usage:
  python3 tests/mixed_spin_nonu1_demo.py --exe /path/to/QS3ED2 --output /new/dir
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np


def local_spin_operators(s):
    d = int(round(2 * s)) + 1
    plus = np.zeros((d, d), dtype=complex)
    for k in range(1, d):
        mz = s - k
        plus[k - 1, k] = np.sqrt(s * (s + 1) - mz * (mz + 1))
    minus = plus.conj().T
    return (
        (plus + minus) / 2,
        (plus - minus) / (2j),
        np.diag([s - k for k in range(d)]).astype(complex),
    )


def product_operators(spins):
    dims = [int(round(2 * s)) + 1 for s in spins]
    local = [local_spin_operators(s) for s in spins]
    full = []
    for i in range(len(spins)):
        triple = []
        for a in range(3):
            op = np.ones((1, 1), dtype=complex)
            for site, d in enumerate(dims):
                op = np.kron(op, local[site][a] if i == site else np.eye(d))
            triple.append(op)
        full.append(triple)
    return full, dims


def couplings(n, transverse):
    # Open, alternating-spin chain; nonzero anisotropy and Gamma_xy already
    # break total Sz conservation even for transverse == 0.
    bonds = []
    for i in range(n - 1):
        bonds.append((i, i + 1, [
            0.77 + 0.025 * i, 0.61 + 0.015 * i, 1.09 - 0.02 * i,
            0.0, 0.0, 0.13 * ((-1) ** i),  # oriented Dz
            0.0, 0.0, 0.07 + 0.01 * i,    # symmetric Gamma_xy
        ]))
    fields = [
        [transverse * (0.4 * ((-1) ** i) + 0.09 * (i + 1)), 0.0,
         0.17 + 0.035 * ((-1) ** i)]
        for i in range(n)
    ]
    return bonds, fields


def dense_hamiltonian(ops, bonds, fields):
    n = len(fields)
    h = np.zeros_like(ops[0][0])
    for i, field in enumerate(fields):
        for a in range(3):
            h += field[a] * ops[i][a]
    # Convention: Dx * (Sy_i Sz_j - Sz_i Sy_j), cyclically.
    for i, j, c in bonds:
        for a in range(3):
            h += c[a] * (ops[i][a] @ ops[j][a])
        for k, a, b in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            h += c[3 + k] * (ops[i][a] @ ops[j][b] -
                             ops[i][b] @ ops[j][a])
            h += c[6 + k] * (ops[i][a] @ ops[j][b] +
                             ops[i][b] @ ops[j][a])
    np.testing.assert_allclose(h, h.conj().T, rtol=0, atol=2e-13)
    return h


def subspace_indices(dims, max_occupation):
    return [
        flat
        for flat, occupations in enumerate(np.ndindex(*dims))
        if sum(occupations) <= max_occupation
    ]


def expect(state, operator):
    z = np.vdot(state, operator @ state)
    if abs(z.imag) > 2e-10:
        raise ValueError(f"Observable expected real but has imaginary part {z.imag}")
    return float(z.real)


def run_case(exe, root, n, transverse, cap):
    spins = [0.5, 1.0] * (n // 2)
    full_cap = sum(int(round(2 * s)) for s in spins)
    current_cap = full_cap if cap == "full" else int(cap)
    if current_cap > full_cap:
        raise ValueError("Invalid NOD cutoff")
    bonds, fields = couplings(n, transverse)
    ops, dims = product_operators(spins)
    hfull = dense_hamiltonian(ops, bonds, fields)
    ids = subspace_indices(dims, current_cap)
    projected = hfull[np.ix_(ids, ids)]
    energies, vectors = np.linalg.eigh(projected)
    ground = np.zeros(hfull.shape[0], dtype=complex)
    ground[ids] = vectors[:, 0]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    label = f"mixed-n{n}-hx{transverse:.2f}-cap{cap}"
    work = root / label
    work.mkdir(exist_ok=False)
    (work / "bonds.dat").write_text(
        "".join(f"{i + 1} {j + 1} " +
                " ".join(f"{x:.17g}" for x in c) + "\n"
                for i, j, c in bonds))
    (work / "fields.dat").write_text(
        "".join(f"{i + 1} " + " ".join(f"{x:.17g}" for x in h) + "\n"
                for i, h in enumerate(fields)))
    (work / "spins.dat").write_text("".join(f"{s:.1f}\n" for s in spins))
    (work / "caps.dat").write_text(
        "".join(f"{int(round(2*s))}\n" for s in spins))
    (work / "pairs.dat").write_text(
        "".join(f"{i + 1} {j + 1}\n" for i, j in pairs))
    deck = f"""&input_parameters
NOS={n}, NODmin=0, NODmax={current_cap}, ALG=3,
L1=1, L2=1, L3=1, L4=1, L5=1, L6=1,
M1=0, M2=0, M3=0, M4=0, M5=0, M6=0,
NO_one={n}, NO_two={len(bonds)}, wk_dim=-1, MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat', FILE_hvec='fields.dat',
FILE_spin='spins.dat', FILE_NODmax='caps.dat',
OUTDIR='./', FILE_wf='./', wr_wf=0, re_wf=0,
cal_lm=1, cal_cf=1, NO_two_cf={len(pairs)}, FILE_two_cf='pairs.dat'
/
&input_lancz
lnc_ene_conv=1d-14, minitr=5, maxitr=1000, itrint=1
/
&input_TRLan
NOE=2, NOK=4, NOM=8, maxitr=1000, lnc_ene_conv=1d-14,
i_vec_min=1, i_vec_max=2
/
"""
    (work / "input.dat").write_text(deck)
    env = dict(os.environ, OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (work / "input.dat").open("rb") as stream:
        result = subprocess.run([str(exe)], stdin=stream, cwd=work,
                                env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=120)
    (work / "stdout.log").write_bytes(result.stdout)
    if result.returncode != 0:
        raise AssertionError(f"{label}: ED2 exit {result.returncode}:\n"
                             + result.stdout.decode(errors="replace")[-4000:])
    match = re.search(r"THS\(k\)\s*=\s*(\d+)",result.stdout.decode(errors="replace"))
    if match is None or int(match.group(1)) != len(ids):
        raise AssertionError(f"{label}: incorrect ED2 sector dimension "
                             f"{match.group(1) if match else 'missing'}, expected {len(ids)}")
    file_energies = np.loadtxt(work / "eigenvalues.dat", ndmin=2)
    assert file_energies.shape == (len(energies), 2), (
        label, "spectrum shape", file_energies.shape, len(energies))
    np.testing.assert_array_equal(file_energies[:, 0], np.arange(1,len(ids) + 1))
    np.testing.assert_allclose(file_energies[:, 1], energies, atol=2e-9, rtol=2e-9)
    # ED2's published Cartesian columns are [xx,yx,zx,xy,yy,zy,xz,yz,zz].
    actual_mag = np.loadtxt(work / "local_mag.dat",ndmin=2)
    expected_mag = np.array([[expect(ground, ops[i][a]) for a in range(3)]
                             for i in range(n)])
    np.testing.assert_array_equal(actual_mag[:, 0], np.arange(1,n + 1))
    np.testing.assert_allclose(actual_mag[:, 1:4], expected_mag, atol=2e-9,rtol=2e-9)
    actual_cf = np.loadtxt(work / "two_body_cf_xyz.dat",ndmin=2)
    expected_cf = np.array([
        [expect(ground, ops[i][a] @ ops[j][b])
         for b in range(3) for a in range(3)]
        for i, j in pairs
    ])
    np.testing.assert_array_equal(actual_cf[:, :2],
                                  np.array([(i + 1,j + 1) for i,j in pairs]))
    np.testing.assert_allclose(actual_cf[:, 2:11], expected_cf, atol=2e-9,rtol=2e-9)
    result_record = {
        "case": label,
        "spins": spins,
        "non_u1": True,
        "n_sites": n,
        "full_space_dimension": int(np.prod(dims)),
        "projected_dimension": len(ids),
        "nodmax": current_cap,
        "transverse_field_parameter": transverse,
        "ground_energy": float(energies[0]),
        "first_excitation_gap": float(energies[1] - energies[0]),
        "ed2_ground_energy": float(file_energies[0, 1]),
        "spectrum_max_error": float(np.max(np.abs(file_energies[:, 1]-energies))),
        "magnetization_max_error": float(np.max(np.abs(actual_mag[:,1:4]-expected_mag))),
        "correlation_max_error": float(np.max(np.abs(actual_cf[:,2:11]-expected_cf))),
        "magnetization_z": expected_mag[:,2].tolist(),
    }
    print(json.dumps(result_record,sort_keys=True),flush=True)
    return result_record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args = p.parse_args()
    exe = args.exe.resolve(strict=True)
    args.output.mkdir(parents=True,exist_ok=False)
    records = []
    for n in (4,6):
        for transverse in (0.15,0.38):
            records.append(run_case(exe,args.output,n,transverse,"full"))
    for transverse in (0.15,0.38):
        records.append(run_case(exe,args.output,4,transverse,2))
    for transverse in (0.15,0.38):
        full = next(x for x in records if x["n_sites"]==4
                    and x["nodmax"]==6 and x["transverse_field_parameter"]==transverse)
        cut = next(x for x in records if x["n_sites"]==4
                   and x["nodmax"]==2 and x["transverse_field_parameter"]==transverse)
        if cut["ground_energy"] < full["ground_energy"] - 2e-9:
            raise AssertionError("Projected ground energy below full ground energy")
    report={"schema":1,"kind":"mixed-spin non-U(1) independent ED2 demonstrator",
            "source_ref":"development candidate; record actual Git SHA via CI logs",
            "cases":records,"passed":len(records),"failed":0}
    (args.output/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"Independent mixed-spin non-U(1) demonstrator: {len(records)}/{len(records)} passed")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
