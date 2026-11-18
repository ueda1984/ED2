# Experimental compact real coefficient cache for repeated ED2 H actions

**Status:** Draft opt-in refinement stacked on PR #22, which in turn
requires the validated PR #21. No change to upstream/main, PR #20
CPC manuscript, final program version or public interface.

## Rationale and contracts

The general `ED2_REAL_CACHED_MATVEC=1` path from PR #22 verifies
that every complex matrix coefficient has **exactly zero imaginary
component** and all projected rows are cached before using real
multiplications. The new, separate
`ED2_PACKED_REAL_CACHE=1` option requires that guard and
copies the coefficient matrix **once** to an allocatable
`real(8)` array with the same per-column/row indices and original
iteration order.

This differs from the prior real-product option, which still
fetches coefficients with a stride-16 complex element. The new
real array has an 8-byte element stride, improving cache utilization
for the repeated matrix-vector loop. It does not change the
`wk_loc` indices, on-the-fly complex Krylov vectors, complex
eigensolver, input/model conventions, basis, matrix construction,
summation order or residual acceptance. Both flags default OFF.

### Safety checks and memory cost

- The packing flag accepts exactly 0/1 and **fails closed**
  unless `ED2_REAL_CACHED_MATVEC=1`.
- Inherited guards reject any nonzero imaginary cache
  coefficient and any partial cache.
- The real array allocates with an explicit status check;
  allocation failure stops rather than silently running an
  invalid or partially packed calculation.
- The original complex cache remains allocated because
  existing public module data and independent bitwise
  matrix proofs expect its availability. The opt-in memory
  **increases** by
  `8 * size(wk_ele,1) * size(wk_ele,2)` bytes.
  At N1000 A1g, `MNTE=64`, `D=4299`, the extra
  memory is **2,235,480 bytes (about 2.13 MiB)**.
- A distinct native marker reports the exact extra bytes:
  `Hamiltonian action: compact-real cache extra_bytes=...`.
- This is a speed/memory tradeoff, not a general
  memory-saving optimization, and must remain user-controlled.

## Numerical proof

A local standalone GNU14.2 copy was built with `-O2 -fopenmp`.
For all 4 1D cubic O_h irreps at L4/L6/L10 T1/T4, plus L10 T8,
the **entire native H*x output** for a deterministic
nonzero-complex test vector had exactly the same SHA256:
**28/28 exact pairs**, packed OFF vs ON, while
`ED2_REAL_CACHED_MATVEC=1` remained ON on both sides.
With N1000/A1g 160-H-call two-stage Lanczos,
the complete `eigenvalues.dat` remained byte-for-byte
identical over eight interleaved runs per setting
at both T1 and T4.

`tests/cpc_packed_real_cache_regression.py` incorporates
the full native vector equality oracle and generic full-space
mixed-spin real-Hamiltonian success test, while explicitly
rejecting nonreal DM, malformed flags, and packing without
prior real-mode validation. All PR13–22 physics tests
continue under the dual gfortran13/14 matrix.

## One shared-host timing sample — NOT PRESTO / journal evidence

Eight interleaved separate old/packed trials per T1/T4,
Linux GNU14.2, pinned N1000/A1g model, 160 H actions,
PR19 and PR21 options ON, all BLAS threads fixed at 1:

| OpenMP threads | Real-only action median | Packed-real action median | Ratio |
|---:|---:|---:|---:|
| 1 | 98.29 ms | 87.12 ms | 1.128x |
| 4 | 55.60 ms | 48.20 ms | 1.153x |

The shared host showed substantial outliers, notably at
T4, and source-level wall instrumentation was a local
development probe only. This is a **directional
experiment**, not a controlled speedup or strong-scaling
claim. Hardware/BLAS/version/CPU binding, memory/RSS,
source archive, raw trial output and operator confirmation
would be needed for formal CPC benchmarking.

## Adoption

Keep this code optional and Draft. Only consider enabling
it by default after:
1. Both GNU13 and GNU14 exact native validation suites pass.
2. Controlled, repeated benchmarks demonstrate a
   performance benefit across compilers and relevant
   real mixed-spin cases.
3. Memory and startup costs are explicitly acceptable for
   the target system; otherwise use PR #22's real-product
   path without the extra allocation or the original
   complex kernel.

No software tag, PR merge or manuscript performance claim
without explicit approval.
