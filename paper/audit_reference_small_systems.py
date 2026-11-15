#!/usr/bin/env python3
"""Independent small-system physics audit for the ED2 CPC manuscript.

This script never invokes ED2 production matrix elements or parses ED2
wavefunctions. It constructs a separate mixed-spin Hamiltonian in SciPy
from the explicitly printed paper coefficients, then compares the
first two eigenvalues and DM-conjugate observable against the frozen
GitHub reference CSV. L=4,6 are dense; L=8 is sparse.

It complements, but does not replace, the existing L=10,12 Fortran +
independent sparse full-physics CI and the separate figure transcription
gate. It is not evidence of an asymptotic/bulk physical extrapolation.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

LENGTHS = (4, 6, 8)
D_GRID = (-.24, -.12, 0., .12, .24)
FD_STEP = 2e-4


def local_spin(s):
    m = np.arange(s, -s - 1, -1, dtype=float)
    sz = np.diag(m).astype(complex)
    sp = np.zeros_like(sz)
    for i in range(1, len(m)):
        sp[i - 1, i] = math.sqrt(s * (s + 1) - m[i] * (m[i] + 1))
    sm = sp.conj().T
    return tuple(sparse.csr_matrix(x) for x in (
        (sp + sm) / 2, (sp - sm) / (2j), sz))


def kron_product(factors):
    result = sparse.csr_matrix([[1.0 + 0j]])
    for factor in factors:
        result = sparse.kron(result, factor, format="csr")
    return result


def physical_operators(length):
    # Independent site-index/orientation convention: j=0,...,L-1.
    spins = [local_spin(.5 if j % 2 == 0 else 1.) for j in range(length)]
    identity = [sparse.eye(t[0].shape[0], dtype=complex, format="csr")
                for t in spins]
    dim = int(np.prod([x.shape[0] for x in identity]))
    zero = sparse.csr_matrix((dim, dim), dtype=complex)

    def term(j, component, k=None, other=None):
        factors = list(identity)
        factors[j] = spins[j][component]
        if k is not None:
            factors[k] = spins[k][other]
        return kron_product(factors)

    H0 = zero.copy()
    K = zero.copy()
    Sz = zero.copy()
    for j in range(length):
        eta = (-1) ** j
        H0 += .15 * (.4 * eta + .1) * term(j, 0)
        H0 += (.17 + .035 * eta) * term(j, 2)
        Sz += term(j, 2)
    for j in range(length - 1):
        eta = (-1) ** j
        for axis, coupling in (
            (0, .77 + .05 * eta),
            (1, .63 + .03 * eta),
            (2, 1.09 - .04 * eta),
        ):
            H0 += coupling * term(j, axis, j + 1, axis)
        gamma = .065 + .015 * eta
        xy = term(j, 0, j + 1, 1)
        yx = term(j, 1, j + 1, 0)
        H0 += gamma * (xy + yx)
        K += eta * (xy - yx)

    # Catch an accidental U(1)-conserving simplification of the example.
    broken = H0 @ Sz - Sz @ H0
    if not np.linalg.norm(broken.data) > 1e-3:
        raise AssertionError("The intended model must break total Sz conservation")
    if np.linalg.norm((H0 - H0.getH()).data) > 1e-10:
        raise AssertionError("Non-Hermitian reference Hamiltonian")
    if np.linalg.norm((K - K.getH()).data) > 1e-10:
        raise AssertionError("DM-conjugate operator is non-Hermitian")
    return H0.tocsr(), K.tocsr(), dim


def lowest_pair(H):
    if H.shape[0] <= 216:
        values, vectors = np.linalg.eigh(H.toarray())
        energy, gap = float(values[0]), float(values[1] - values[0])
        vector = vectors[:, 0]
    else:
        v0 = np.ones(H.shape[0], dtype=complex)
        values, vectors = eigsh(H, k=2, which="SA", tol=1e-12,
                               maxiter=50000, ncv=40, v0=v0)
        indices = np.argsort(values)
        energy = float(values[indices[0]])
        gap = float(values[indices[1]] - values[indices[0]])
        vector = vectors[:, indices[0]]
    residual = float(np.linalg.norm(H @ vector - energy * vector))
    if residual > 2e-8 or gap <= 0:
        raise AssertionError("Unaccepted independent eigenpair")
    return energy, gap, vector


def load_csv(path):
    with path.open(newline="", encoding="utf-8") as source:
        data = list(csv.DictReader(source))
    selected = {}
    for row in data:
        length = int(row["sites"])
        if length not in LENGTHS:
            continue
        d = float(row["D"])
        matches = [x for x in D_GRID if abs(x - d) < 1e-12]
        if len(matches) != 1:
            raise ValueError("Unexpected finite-size DM parameter")
        key = (length, matches[0])
        if key in selected:
            raise ValueError("Duplicate independent CSV point")
        selected[key] = row
    if len(selected) != 15:
        raise ValueError("Expected exactly 15 reference points for L=4,6,8")
    return selected


def run(csv_path):
    data = load_csv(csv_path)
    report = []
    for length in LENGTHS:
        H0, K, dim = physical_operators(length)
        values = {}
        for d in D_GRID:
            E, gap, vector = lowest_pair(H0 + d * K)
            chirality = float(np.vdot(vector, K @ vector).real)
            ref = data[length, d]
            if int(ref["full_dimension"]) != dim:
                raise AssertionError("Reference dimension mismatch")
            errE = abs(E - float(ref["independent_E0"]))
            errGap = abs(gap - float(ref["first_gap"]))
            errK = abs(chirality - float(ref["staggered_chirality"]))
            if errE > 2e-8 or errGap > 2e-7 or errK > 2e-7:
                raise AssertionError(
                    f"Physics discrepancy at L={length} D={d}: "
                    f"E={errE:.3g}, gap={errGap:.3g}, K={errK:.3g}")
            values[d] = (E, chirality, gap)
            report.append(dict(sites=length, D=d, dimension=dim,
                               E0=E, gap=gap, K=chirality,
                               energy_error=errE, gap_error=errGap,
                               chirality_error=errK))
        if any(values[D_GRID[i + 1]][1] > values[D_GRID[i]][1] + 1e-7
               for i in range(len(D_GRID) - 1)):
            raise AssertionError("DM-conjugate expectation violates concavity")
        e_plus = lowest_pair(H0 + FD_STEP * K)[0]
        e_minus = lowest_pair(H0 - FD_STEP * K)[0]
        e_zero = values[0.][0]
        K_zero = values[0.][1]
        deriv = (e_plus - e_minus) / (2 * FD_STEP)
        chi = -(e_plus - 2 * e_zero + e_minus) / (FD_STEP**2 * (length - 1))
        expected = float(data[length, 0.]["chi_kappa_per_bond"])
        if abs(deriv - K_zero) > 2e-7 or abs(chi - expected) > 5e-5:
            raise AssertionError(
                f"Finite-difference/HF check failed for L={length}: "
                f"derivative={deriv}, chirality={K_zero}, chi={chi}, ref={expected}")
        print(f"L={length}: full dim={dim}, E0(0)={e_zero:.12f}, "
              f"gap(0)={values[0.][2]:.9f}, chi(0)={chi:.9f}, "
              f"HF error={abs(deriv-K_zero):.2e}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path,
                        default=Path("docs/development/cpc-stationary-dm-reference.csv"))
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()
    report = run(args.csv)
    if args.report is not None:
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"CPC INDEPENDENT PHYSICS AUDIT: {len(report)}/15 "
          "L=4,6,8 energies, gaps, and chirality values agree with CSV")
    print("CPC INDEPENDENT PHYSICS AUDIT: 3/3 zero-D HF and susceptibility checks passed")


if __name__ == "__main__":
    main()
