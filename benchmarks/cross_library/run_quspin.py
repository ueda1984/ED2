#!/usr/bin/env python3
"""QuSpin U(1)-only pilot; NOT a spatial-symmetry-projected comparison."""
import argparse
import json
import os
import resource
import time
from reference import bonds_for_shape, case_metadata, parse_shape


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shape", required=True)
    ap.add_argument("--ndown", type=int, required=True)
    ap.add_argument("--max-dim", type=int, default=1000000)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--maxiter", type=int, default=3000)
    ap.add_argument("--matvec-repetitions", type=int, default=25)
    a = ap.parse_args()
    shape = parse_shape(a.shape)
    meta = case_metadata(shape, a.ndown)
    out = dict(meta, implementation="QuSpin/SciPy-CSR", symmetry="U1_only",
               threads_omp=os.environ.get("OMP_NUM_THREADS", "unset"),
               tolerance_requested=a.tol)
    if a.max_dim < 1 or a.matvec_repetitions < 1 or a.tol <= 0:
        raise ValueError("invalid resource/tolerance settings")
    if meta["dimension_u1"] > a.max_dim:
        print(json.dumps(dict(out, status="SKIPPED_DIMENSION_SAFETY_CAP", max_dim=a.max_dim),
                         sort_keys=True))
        return

    import numpy as np
    import scipy
    import quspin
    from scipy.sparse.linalg import eigsh
    from quspin.basis import spin_basis_general
    from quspin.operators import hamiltonian
    out.update(numpy_version=np.__version__, scipy_version=scipy.__version__,
               quspin_version=getattr(quspin, "__version__", "unknown"))
    start = time.perf_counter()
    basis = spin_basis_general(meta["sites"], Nup=meta["sites"]-a.ndown,
                               S="1/2", pauli=0)
    tbasis = time.perf_counter()
    if int(basis.Ns) != meta["dimension_u1"]:
        raise AssertionError("wrong U1 sector")
    bonds = bonds_for_shape(shape)
    static = [
        ["zz", [[1.0, i, j] for i, j in bonds]],
        ["+-", [[0.5, i, j] for i, j in bonds]],
        ["-+", [[0.5, i, j] for i, j in bonds]],
    ]
    H = hamiltonian(static, [], basis=basis, dtype=np.float64,
                    check_herm=False, check_pcon=False, check_symm=False)
    csr = H.tocsr()
    tcsr = time.perf_counter()
    dim = int(csr.shape[0])
    assert dim == meta["dimension_u1"]
    if dim <= 2:
        values, vecs = np.linalg.eigh(csr.toarray())
        e0, vec = float(values[0]), vecs[:, 0]
    else:
        v0 = np.linspace(1.0, 2.0, dim)
        v0 /= np.linalg.norm(v0)
        values, vecs = eigsh(csr, k=1, which="SA", tol=a.tol, v0=v0, maxiter=a.maxiter)
        e0, vec = float(values[0]), vecs[:, 0]
    tsolve = time.perf_counter()
    residual = float(np.linalg.norm(csr @ vec-e0*vec))
    x = np.linspace(.1, 1.1, dim)
    x /= np.linalg.norm(x)
    for _ in range(3):
        y = csr @ x
    t0 = time.perf_counter()
    for _ in range(a.matvec_repetitions):
        y = csr @ x
    t1 = time.perf_counter()
    out.update(status="PASS" if residual <= max(1e-8, 20*a.tol) else "BAD_RESIDUAL",
               dimension_constructed=dim, nnz=int(csr.nnz),
               E0=e0, direct_residual=residual, basis_seconds=tbasis-start,
               csr_seconds=tcsr-tbasis, solve_seconds=tsolve-tcsr,
               time_to_eigenpair_seconds=tsolve-start,
               matvec_seconds=(t1-t0)/a.matvec_repetitions,
               matvec_repetitions=a.matvec_repetitions,
               peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(json.dumps(out, sort_keys=True))
    if out["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
