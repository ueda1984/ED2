#!/usr/bin/env python3
"""QuSpin pilot for large-site U(1), Gamma or Gamma/Oh/A1g sectors.

Every mode is resource-gated by the *raw* fixed-magnetization dimension,
because the basis constructor may scan that space before symmetry reduction.
The Gamma/Oh/A1g sector uses a rigorously factorized set of four point
generators, plus translations; actual QuSpin basis/energy must still be
compared with the independent scientific oracle before claiming equality.
"""
import argparse
import json
import math
import os
import resource
import time
from reference import bonds_for_shape, case_metadata, parse_shape


def translations(shape):
    n = math.prod(shape)
    stride = [math.prod(shape[:i]) for i in range(len(shape))]
    result = {}
    for d, side in enumerate(shape):
        step = stride[d]
        # A generator that translates lattice site i forward along axis d.
        perm = [i + step if (i // step) % side + 1 < side
                else i - (side - 1) * step for i in range(n)]
        assert len(set(perm)) == n
        result["xyz"[d] + "block"] = (perm, 0)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shape", required=True)
    p.add_argument("--ndown", required=True, type=int, choices=[2, 3])
    p.add_argument("--model", choices=["heisenberg-nn", "cubic-j1j2-xxz"],
                   default="heisenberg-nn")
    p.add_argument("--sector", choices=["u1", "gamma", "gamma-a1g"], default="u1")
    p.add_argument("--mode", choices=["dimensions", "construct", "basis", "solve"],
                   default="dimensions")
    p.add_argument("--max-raw-dim", type=int, default=1_000_000)
    p.add_argument("--max-csr-nnz", type=int, default=10_000_000)
    p.add_argument("--prefix-bits", type=int, default=0)
    p.add_argument("--nup-convention", choices=["negative","positive"], default="negative",
                   help="Negative Nup=-m-1 is QuSpin's documented near-saturation shorthand")
    p.add_argument("--tol", type=float, default=1e-10)
    p.add_argument("--maxiter", type=int, default=3000)
    p.add_argument("--matvec-repetitions", type=int, default=25)
    args = p.parse_args()
    shape = parse_shape(args.shape)
    meta = case_metadata(shape, args.ndown)
    if args.model == "cubic-j1j2-xxz":
        if len(shape) != 3 or len(set(shape)) != 1:
            p.error("cubic-j1j2-xxz requires L,L,L shape")
        from export_cubic_xxz import metadata as xxz_metadata
        special = xxz_metadata(shape[0], args.ndown)
        meta.update(model=special["model"], physical_sha256=special["physical_sha256"],
                    bonds=special["bond_count"], all_up_energy=special["all_up_energy"],
                    source=special["source"])
    expected = meta["dimension_u1"]
    if args.sector in ("gamma", "gamma-a1g"):
        from symmetry_dimensions import count_blocks
        if args.sector == "gamma-a1g" and (len(shape) != 3 or len(set(shape)) != 1):
            p.error("Gamma/Oh/A1g is only defined for the cubic L,L,L workload")
        expected_key = "gamma_translation" if args.sector == "gamma" else "gamma_point_A1"
        expected = count_blocks(shape)["m%d" % args.ndown][expected_key]
    out = dict(meta, implementation="QuSpin/spin_basis_general",
               sector=args.sector, mode=args.mode, expected_dimension=expected,
               model_choice=args.model, nup_convention=args.nup_convention,
               threads_omp=os.environ.get("OMP_NUM_THREADS", "unset"))
    if args.mode == "dimensions":
        print(json.dumps(dict(out, status="DIMENSIONS_ONLY"), sort_keys=True))
        return

    # QuSpin may enumerate C(N,m) states even when Gamma reduces to a few
    # thousand states; always apply guard to *raw* U(1) dimension.
    if args.mode != "construct" and meta["dimension_u1"] > args.max_raw_dim:
        print(json.dumps(dict(out, status="SKIPPED_RESOURCE_CAP",
                              max_raw_dim=args.max_raw_dim), sort_keys=True))
        return
    if args.prefix_bits < 0 or args.prefix_bits > 16 or args.tol <= 0 or args.maxiter < 1:
        raise ValueError("invalid settings")
    import quspin
    import numpy as np
    import scipy
    from quspin.basis import spin_basis_general
    from quspin.operators import hamiltonian
    from scipy.sparse.linalg import eigsh
    out.update(quspin_version=getattr(quspin, "__version__", "unknown"),
               numpy_version=np.__version__, scipy_version=scipy.__version__,
               prefix_bits=args.prefix_bits)
    t0 = time.perf_counter()
    kwargs = translations(shape) if args.sector != "u1" else {}
    if args.sector == "gamma-a1g":
        from cubic_oh_generators import verify_factorization, cubic_oh_blocks
        verify_factorization()
        kwargs.update(cubic_oh_blocks(shape[0]))
    nup = (-args.ndown - 1 if args.nup_convention == "negative"
           else meta["sites"]-args.ndown)
    # Disable automatic block-table expansion; use explicit N_p (for dilute
    # systems the default prefix table can be poor). This is a test setting,
    # NOT a demonstrated optimal QuSpin configuration.
    basis = spin_basis_general(meta["sites"], Nup=nup,
                               S="1/2", pauli=0, Ns_block_est=expected,
                               make_basis=False, **kwargs)
    if args.mode == "construct":
        print(json.dumps(dict(out, status="CONSTRUCTOR_ONLY",
                              basis_dtype=str(basis.dtype),
                              constructor_seconds=time.perf_counter()-t0,
                              nup_argument=nup,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                         sort_keys=True))
        return
    basis.make(N_p=args.prefix_bits)
    tbasis = time.perf_counter()
    if int(basis.Ns) != expected:
        raise AssertionError("symmetry dimension mismatch: got %d expected %d" %
                             (basis.Ns, expected))
    if args.mode == "basis":
        print(json.dumps(dict(out, status="BASIS_ONLY",
                              dimension_constructed=int(basis.Ns),
                              basis_seconds=tbasis-t0,
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                         sort_keys=True))
        return

    # Explicit second resource guard for CSR storage. No implied matrix-free
    # result if this is skipped.
    coordination = 18 if args.model == "cubic-j1j2-xxz" else 2*len(shape)
    max_nnz_bound = expected * (1 + args.ndown*coordination)
    if max_nnz_bound > args.max_csr_nnz:
        print(json.dumps(dict(out, status="SKIPPED_CSR_SAFETY_CAP",
                              max_nnz_bound=max_nnz_bound,
                              max_csr_nnz=args.max_csr_nnz,
                              basis_seconds=tbasis-t0), sort_keys=True))
        return
    if args.model == "cubic-j1j2-xxz":
        from export_cubic_xxz import cubic_xxz_bonds
        bonds = cubic_xxz_bonds(shape[0])
    else:
        bonds = [(i, j, 1.0, 1.0) for i, j in bonds_for_shape(shape)]
    static = [
        ["zz", [[jzz, i, j] for i, j, jxy, jzz in bonds]],
        ["+-", [[jxy/2.0, i, j] for i, j, jxy, jzz in bonds]],
        ["-+", [[jxy/2.0, i, j] for i, j, jxy, jzz in bonds]],
    ]
    H = hamiltonian(static, [], basis=basis, dtype=np.float64,
                    check_herm=False, check_pcon=False, check_symm=False)
    csr = H.tocsr()
    tcsr = time.perf_counter()
    d = int(csr.shape[0])
    if d != expected:
        raise AssertionError("CSR dimension mismatch")
    if d <= 2:
        w,v = np.linalg.eigh(csr.toarray())
        e0, vec = float(w[0]),v[:,0]
    else:
        v0 = np.linspace(1.0, 2.0, d)
        v0 /= np.linalg.norm(v0)
        w,v = eigsh(csr, k=1, which="SA", tol=args.tol, v0=v0, maxiter=args.maxiter)
        e0, vec = float(w[0]),v[:,0]
    tsolve = time.perf_counter()
    direct_residual = float(np.linalg.norm(csr @ vec-e0*vec))
    x = np.linspace(.1, 1.1, d)
    x /= np.linalg.norm(x)
    for _ in range(3):
        y = csr @ x
    t1 = time.perf_counter()
    for _ in range(args.matvec_repetitions):
        y = csr @ x
    t2 = time.perf_counter()
    out.update(status="PASS" if direct_residual <= max(1e-8, 20*args.tol)
               else "BAD_RESIDUAL",
               dimension_constructed=d, nnz=int(csr.nnz), E0=e0,
               direct_residual=direct_residual, basis_seconds=tbasis-t0,
               csr_seconds=tcsr-tbasis, solve_seconds=tsolve-tcsr,
               time_to_eigenpair_seconds=tsolve-t0,
               matvec_seconds=(t2-t1)/args.matvec_repetitions,
               matvec_repetitions=args.matvec_repetitions,
               peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(json.dumps(out, sort_keys=True))
    if out["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
