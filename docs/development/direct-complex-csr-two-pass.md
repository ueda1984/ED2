# Experimental direct two-pass complex CSR Hamiltonian construction

**Status: experimental, default OFF, independent general-ED2 optimization.**
This branch starts from main *after* the general complex CSR cache, typed
BLAS/LAPACK, lazy H*x scratch and coefficient-certified capacity changes.
It does not depend on the cubic three-magnon acceleration stack.

## Why the earlier "release after CSR packing" is not enough

The mainline `ED2_CSR_COMPLEX_CACHE=1` mode builds the entire original
padded complex H cache, then allocates CSR values and copies all occupied
records into CSR. If `ED2_CSR_RELEASE_PADDED_CACHE=1` is also set, the
padded cache is deallocated only **after** the copy. Thus startup still
requires both matrices simultaneously, even if the long-lived H action
would only need the CSR representation.

A new explicit opt-in **direct two-pass** construction avoids this full
padded representation in the first place:

```bash
ED2_CSR_COMPLEX_CACHE=1 \
ED2_CSR_RELEASE_PADDED_CACHE=1 \
ED2_CSR_DIRECT_BUILD=1 \
ED2_CACHE_COEFF_BOUND=1 \
./source/QS3ED2 < input.dat
```

The coefficient-bound option is separately optional but helpful for
reducing the temporary batch allocations. Direct build **requires** a full
matrix cache (`wk_dim=lv`) and both CSR and release options. It rejects
partial caching, malformed or oversized flags and missing prerequisites.
Existing default and partial-cache behavior is untouched.

## Construction algorithm: exactly one physical row generator

The original Fortran calculation of each Hamiltonian row is **not copied
into a new, independently maintained implementation**. Instead, the
existing ~200-line one-site/two-site transition builder becomes one
internal `fill_cache_batch(first,last,...)` routine shared by:
- the legacy full padded representation, with first=1, last=wk_dim;
- the new two-pass builder, operating on fixed **64-row batches**.

Within the shared row routine, every original coefficient, ordered row
entry, phase, representative lookup and normalization formula is
unchanged. The only row-storage indexing difference is
`ib=i-first+1` instead of `i`; for the default full padded path
`first=1`, so this is the **same** index.

Direct mode performs:

1. **Count pass.** Allocate one 64-row `wk_loc`/`wk_ele` batch,
   generate its exact physical entries, record the number of
   non-sentinel terms per row in `csr_rowptr`, and immediately
   discard the batch. Repeat until all rows are counted. No
   global padded matrix ever exists.
2. **Exact allocation.** Allocate `csr_target` and `csr_coeff`
   to precisely the total count. Only original target IDs and
   complex numbers are stored; no deduplication, reordering,
   truncation or tolerance pruning.
3. **Fill pass.** Generate the same 64-row batches a second time
   with the **same physical formulas** and copy all occupied entries
   to their CSR row spans in original row order. Count mismatches,
   invalid target IDs, allocation errors or 32-bit row-offset
   overflow fail closed.

After both passes, `wk_loc` and `wk_ele` are **unallocated** and
the existing complex CSR H*x action is reused, exactly as in the
explicit padded-release mode.

**Tradeoff:** the physical row construction is performed twice and may
take approximately twice the one-time builder work. This mode
prioritizes startup cache memory, NOT minimum setup latency, and
should stay an explicit opt-in until controlled hardware testing
demonstrates where the tradeoff is worthwhile. For small problems
where a batch can contain the full matrix, it may save no memory
or have a larger intermediate footprint.

## Precise bounded payload, not measured RSS

Native stdout reports:

    Hamiltonian action: csr-complex cache entries=N, extra_bytes=B
    Hamiltonian action: csr-direct-build passes=2, max_staging_bytes=S, peak_cache_payload_upper_bytes=P

Here `P=B+S`: it is a conservative upper bound for **these
particular Hamiltonian cache arrays only**, not maximum operating
system RSS, total program memory, BLAS buffers, state lists,
OpenMP private data, or the size of the old full-H algorithm's
other allocations. Pass 1 only needs row pointers plus one batch;
pass 2 holds exact CSR plus one batch; therefore B+S bounds the
simultaneously allocated cache-array payload.

Local GNU14.2 N10 alternating S=1/2,S=1, full Hilbert space
D=7776, coefficient-bound ON, one thread:

| Hamiltonian | Existing full padded+CSR maximum payload | Direct CSR cache payload upper bound |
|---|---:|---:|
| Real XYZ or complex DM | 12,731,908 B | **3,940,228 B** |
| U(1) XXZ | 4,074,628 B | **1,144,068 B** |
| Full 9-component tensor | 21,596,548 B | **7,252,228 B** |

For XYZ/complex DM this is a **69% reduction** in the maximum
co-resident cache-array payload relative to existing CSR construction
*when coefficient bounding is already ON*. It is not a 69% peak
RSS reduction. Direct CSR also avoids the unbounded original
padded+CSR peak (18,330,628 B for this case).

## Scientific evidence

The source-level change is unusually invasive because the original
physical row code moves into a shared routine. Therefore the
acceptance test has two parts: default-off equivalence to a
**separately compiled original source** and new-path equivalence
within the candidate.

In independent local GNU Fortran 14.2 runs:
- **128/128** old-source versus refactored-source native
  full H*x bitstreams matched, across N4/6/8/10 × four physical
  model families × T1/T4 × bound OFF/ON × padded/CSR-release,
  both `-O2` **and** the historical `-Ofast`.
- **32/32** direct-CSR physical cases with six storage/cap
  combinations produced exactly identical full complex
  H*x output bits under `-O2`, with two nonreal complex
  test vectors in each case.
- **32/32** physical cases under `-Ofast` passed exact
  coefficient-bound OFF/ON equality *within every layout*
  and the explicitly allowed floating-point tolerance
  for reassociation between the padded and CSR loops.
- N4/N6 independently constructed dense H matrices agreed
  with every H*x output; all five fail-closed invalid-input
  cases passed.
- Unmodified historical `make check` and all `ex1`–`ex8`
  examples passed with GNU Fortran 14.2 `-Ofast`.

The checked-in `tests/cpc_csr_direct_build_regression.py` and
`.github/workflows/main-csr-direct.yml` require native
GNU13/GNU14 success with `-O2` and `-Ofast` before any
adoption. The earlier general complex CSR, LP64 BLAS/LAPACK,
coefficient-bound and release checks remain active in their
own PR CI workflows.

No PRESTO computation, journal-ready performance number,
automatic default enable, upstream QS-Cube main change,
or final release is implied. The default remains the
original padded cache, with no new externally visible
arrays unless this explicit experimental opt-in is used.
