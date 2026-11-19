#!/usr/bin/env python3
"""Diagonalize an independently generated cubic J1-J2 XXZ A1g matrix.

Scientific oracle only: not an inter-library runtime benchmark. Requires
NumPy + SciPy. Supports side=4 and 10 in CI, with deliberate dimension cap.
"""
import argparse
import json
import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import eigsh
from export_cubic_xxz import metadata
from symmetry_dimensions import count_blocks


REFERENCES = {
    4: 40.345908939159166,   # independent full-orbit 44x44 dense oracle
    10: 847.3078940655731,  # previous ED2 PR16 N1000 Gamma/Oh A1g pilot
}


def solve(path, side, tol=1e-12):
    if side not in (4, 10):
        raise ValueError("CI oracle supports side 4 or 10, no performance extrapolation")
    h = mmread(path).tocsr()
    expected = count_blocks((side, side, side))["m3"]["gamma_point_A1"]
    if h.shape != (expected, expected) or expected > 5000:
        raise ValueError("incorrect block size or exceeds dimension safety cap")
    if h.nnz > 300000:
        raise ValueError("matrix exceeds independent oracle nonzero safety cap")
    asym = h - h.T
    symmetry_error = float(np.max(np.abs(asym.data))) if asym.nnz else 0.0
    if symmetry_error > 1e-10:
        raise AssertionError("matrix is not Hermitian")
    v0 = np.linspace(1., 2., h.shape[0], dtype=np.float64)
    v0 /= np.linalg.norm(v0)
    eig, vector = eigsh(h, k=1, which="SA", tol=tol, v0=v0, maxiter=10000)
    energy = float(eig[0])
    residual = float(np.linalg.norm(h @ vector[:, 0] - energy * vector[:, 0]))
    delta = abs(energy - REFERENCES[side])
    if residual > 1e-8 or delta > 1e-8:
        raise AssertionError("reference E/residual outside independent tolerance")
    return dict(metadata(side, 3), status="PASS", symmetry="Gamma/Oh/A1g",
                projected_dimension=expected, matrix_nnz=int(h.nnz),
                matrix_asymmetry_max=symmetry_error, E0=energy,
                direct_residual=residual, reference_energy=REFERENCES[side],
                absolute_energy_difference=delta,
                scientific_scope="physics oracle, no cross-library timing claim")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--side", type=int, choices=[4, 10], required=True)
    p.add_argument("--matrix", required=True)
    p.add_argument("--tol", type=float, default=1e-12)
    args = p.parse_args()
    print(json.dumps(solve(args.matrix, args.side, args.tol), sort_keys=True))
