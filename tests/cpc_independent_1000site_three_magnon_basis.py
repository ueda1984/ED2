#!/usr/bin/env python3
"""Construct the 1000-site / three-magnon Gamma Hamiltonian independently.

This is a *structure* validation: exact 166167 translation-orbit states,
3151137 nonzero real-symmetric matrix entries and a Hermiticity check.
It does NOT run ED2 or claim a converged 1000-site three-magnon ground energy.

Because the 10x10x10 translation group has no element of order three,
no three-magnon configuration has a nontrivial translation stabilizer.
Fixing one occupied site at the origin gives C(999,2) anchored pairs,
three per translation orbit. Canonicalize by re-rooting at each magnon;
the resulting orbit basis has exactly C(1000,3)/1000 = 166167 states.
"""
import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy import sparse

SIDE = 10
VOLUME = SIDE**3
BONDS = 3 * VOLUME
EXPECTED_DIM = 166167
EXPECTED_NNZ = 3151137

def pair(a, b):
    return (a, b) if a < b else (b, a)


def build():
    start = time.perf_counter()
    coords = [(i % SIDE, (i // SIDE) % SIDE, i // (SIDE**2))
              for i in range(VOLUME)]
    difference = np.empty((VOLUME, VOLUME), dtype=np.uint16)
    for a, (ax, ay, az) in enumerate(coords):
        for b, (bx, by, bz) in enumerate(coords):
            difference[b, a] = ((bx - ax) % SIDE) + SIDE * ((by - ay) % SIDE) \
                + SIDE**2 * ((bz - az) % SIDE)

    def canonical(a, b, c):
        return min(pair(int(difference[b, a]), int(difference[c, a])),
                   pair(int(difference[a, b]), int(difference[c, b])),
                   pair(int(difference[a, c]), int(difference[b, c])))

    # No global C(1000,3) enumeration is needed.
    reps = []
    for a in range(1, VOLUME):
        negative_a = int(difference[0, a])
        for b in range(a + 1, VOLUME):
            k = min((a, b),
                    pair(negative_a, int(difference[b, a])),
                    pair(int(difference[0, b]), int(difference[a, b])))
            if k == (a, b):
                reps.append(k)

    if len(reps) != EXPECTED_DIM or len(reps) * VOLUME != math.comb(VOLUME, 3):
        raise AssertionError("Canonical translation-orbit dimension mismatch")
    index = {p: i for i, p in enumerate(reps)}
    if len(index) != len(reps):
        raise AssertionError("Duplicate canonical orbit")

    neighbors = []
    for x, y, z in coords:
        neighbors.append((
            (x + 1) % SIDE + SIDE*y + SIDE**2*z,
            (x - 1) % SIDE + SIDE*y + SIDE**2*z,
            x + SIDE*((y + 1) % SIDE) + SIDE**2*z,
            x + SIDE*((y - 1) % SIDE) + SIDE**2*z,
            x + SIDE*y + SIDE**2*((z + 1) % SIDE),
            x + SIDE*y + SIDE**2*((z - 1) % SIDE),
        ))

    rows, cols, vals = [], [], []
    for col, (a, b) in enumerate(reps):
        sites = (0, a, b)
        occupied = set(sites)
        adjacent_pairs = sum(p in neighbors[q] for q, p in ((0, a), (0, b), (a, b)))
        diagonal = BONDS / 4.0 - 9.0 + adjacent_pairs
        entries = {col: diagonal}
        for site in sites:
            for dest in neighbors[site]:
                if dest in occupied:
                    continue
                hopped = tuple(dest if p == site else p for p in sites)
                row = index[canonical(*hopped)]
                entries[row] = entries.get(row, 0.0) + 0.5
        for row, val in entries.items():
            rows.append(row)
            cols.append(col)
            vals.append(val)

    h = sparse.coo_matrix(
        (vals, (rows, cols)), shape=(EXPECTED_DIM, EXPECTED_DIM),
        dtype=np.float64
    ).tocsr()
    anti = h - h.transpose()
    error = float(np.max(np.abs(anti.data), initial=0.0))
    if h.nnz != EXPECTED_NNZ or error > 1e-12:
        raise AssertionError(f"Unexpected sparse Hamiltonian nnz={h.nnz}, antisym={error}")

    report = {
        "schema": 1,
        "status": "operator_structure_verified_not_diagonalized",
        "model": "1000-site periodic nearest-neighbor S=1/2 isotropic Heisenberg",
        "shape": [SIDE]*3,
        "n_down": 3,
        "momentum": [0, 0, 0],
        "number_of_periodic_bonds": BONDS,
        "Jxy": 1.0,
        "Jz": 1.0,
        "fixed_n_dimension": math.comb(VOLUME, 3),
        "origin_anchored_pairs": math.comb(VOLUME-1, 2),
        "exact_momentum_dimension": len(reps),
        "matrix_nonzeros": int(h.nnz),
        "hermiticity_absolute_error": error,
        "smallest_diagonal_element": float(h.diagonal().min()),
        "largest_diagonal_element": float(h.diagonal().max()),
        "construction_wall_s_uncontrolled": round(time.perf_counter() - start, 6),
        "source_of_hamiltonian": "independent relative-magnon displacements",
        "warning": "A sparse matrix is constructed, but no 1000-site three-magnon eigenvalue "
                   "has been computed by QS3, ED2 or this test."
    }
    return h, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--matrix", type=Path, default=None,
                        help="Optional NPZ for a controlled-node eigensolver run")
    a = parser.parse_args()
    mat, report = build()
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if a.matrix is not None:
        a.matrix.parent.mkdir(parents=True, exist_ok=True)
        sparse.save_npz(a.matrix, mat, compressed=True)
    print("1000-site / 3-magnon independent Hamiltonian structure verified")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
