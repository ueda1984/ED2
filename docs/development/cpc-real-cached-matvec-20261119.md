# Opt-in exactly real cached matrix-vector action (experimental PR #22)

**Status:** experimental numerical fast path; default remains unchanged.
This branch is stacked on the fully tested PR #21, which itself
depends on the exact cubic basis and direct-hop improvements. It does not
alter the published-program manuscript branch, upstream, or main.

## Motivation

When a full sparse ED2 Hamiltonian is cached, one Lanczos iteration
performs a column/row gather over `wk_loc(j,i)` and complex products
`wk_ele(j,i)*v1(wk_loc(j,i))`. The matrix elements are stored as
complex even for models with a rigorously **real** representation
(e.g. real XXZ exchange in the four 1D O_h irreps of the cubic test).

A temporary local GNU Fortran 14.2 profile of N=1000/A1g with all
PR19/PR21 optimizations enabled counted 160 native Hamiltonian
applications per two-stage Lanczos run. Matvec wrapper elapsed wall
was typically approximately **85–93% of the Lanczos stage**, so this
repeated operation has better performance leverage than further tiny
one-time basis setup changes. Shared-host timings are noisy and are
not PRESTO performance evidence.

## Narrow correctness contract

The explicit `ED2_REAL_CACHED_MATVEC=1` flag (OFF by default) is
examined only at **the end of native cached-Hamiltonian construction**
and is rejected unless:

1. The full projected matrix is cached (`wk_dim == lv`).
2. Every stored complex matrix element has **exactly zero imaginary
   part**, including zero-filled padding. This check is performed
   once, never by guessing a model family or symmetry.
3. The flag is exactly `0` or `1`.

Once these conditions are met, the cached portion of
`ham_to_vec_wave_vector` uses
`dble(wk_ele(j,i))*v1(wk_loc(j,i))` instead of the full
complex-coefficient multiplication. It retains exactly the
historical per-row summation order and sentinel exit; no basis,
Hamiltonian, solver convergence criterion, phase, input format,
or physical operator is changed. Off path is byte-for-byte the
preexisting code. The distinctive native marker is:

    Hamiltonian action: strictly-real cached coefficients

This is a **general full-cache real-Hamiltonian** optimization, not
limited by code to the cubic 3-magnon 1D-irrep model.

## Evidence

An independent GNU14.2 local build of an experimental copy of
the implementation completed exact full native H-action comparisons:
- **28/28** complete output-vector bitwise SHA256 pairs for L4, L6,
  L10 × four 1D O_h sectors × T1/T4, plus L10 × T8;
  the test input vector includes deterministic **nonzero complex**
  components, so equality is not an artifact of a real input.
- For the actual N1000/A1g two-stage solver, all 8 old/new
  trials per T1 and T4 produced **byte-identical full
  `eigenvalues.dat` files**, with the same 160 H calls.
- An ordinary mixed-spin full-space, non-U(1) Hamiltonian with
  exactly real entries produced identical eigenvalue files;
  an actual complex DM coefficient was rejected on request.
- Partial cache mode and malformed opt-in values were rejected.

The checked-in regression `tests/cpc_real_cached_matvec_regression.py`
repeats the exact native vector oracle, and includes the independent
generic mixed-spin real/nonreal/partial/error-path cases in the
GNU13+GNU14 CI matrix. Prior PR13–21 full H-cache and Ritz proofs
must remain green.

### Exploratory performance (NOT publication evidence)

An interleaved **8 repetitions per mode** on a shared GNU14.2 Linux
machine with source built `-O2 -fopenmp`, fixed N1000/A1g,
160 H calls and all earlier cubic options ON gave:

| Threads | Old action cumulative median | Real-only action median | Old/New |
|---|---:|---:|---:|
| T1 | 91.77 ms | 87.28 ms | 1.051x |
| T4 | 47.04 ms | 42.30 ms | 1.112x |

End-to-end and Lanczos medians also improved in this particular set,
but shared-resource contention caused large run-to-run outliers.
The source-level profiler used to collect action-only wall time was
a **temporary local diagnostic**, not part of the release or a claim
of reproducible operator-attested timing. Re-measure on independently
controlled hardware with a complete immutable source/provenance
package before any performance claim. Never enable this shortcut on
Hamiltonians with nonzero imaginary coefficients.

## Priorities after this PR

The repeated Hamiltonian action is likely to remain important in
Lanczos; however, more invasive optimizations such as operator
compression, coefficient reordering, CSR conversion, fused stages or
parallel vector kernels must be judged against the **complete native
complex-vector bitwise oracle** and independently established physical
residuals. Do not optimize away the independent numerical contracts.

No merge or release tag without explicit maintainer/author approval.
