# Experimental release of padded original H cache after complex CSR packing

**Scope:** independent small opt-in improvement of the general complex
CSR backend introduced in fork-main PR26. Requires **both**
environment flags:

    ED2_CSR_COMPLEX_CACHE=1 ED2_CSR_RELEASE_PADDED_CACHE=1 ./source/QS3ED2 < input.dat

Both modes remain OFF unless explicitly specified; existing plain
`ED2_CSR_COMPLEX_CACHE=1` continues to **retain** the original
`wk_ele`/`wk_loc` arrays and keep module-level access compatible.

## Purpose and exact behavior

Historically, a full Hamiltonian cache in ED2 is a padded fixed-capacity
pair of 2D arrays `wk_ele` (complex coefficient) and `wk_loc`
(integer target). PR26's CSR alternative copies *every* non-sentinel
coefficient, including stored zeros, preserving order and numerical
arithmetic but additionally retains the padded originals.

This optional new flag **deallocates the two original arrays** only
after validating and completely constructing the CSR representation
for **all rows**. The repeated H-action calculation accesses only
`csr_rowptr`, `csr_target`, `csr_coeff`, preserving all complex
IEEE input/output bits. The original complex arrays are no longer
available to downstream code after this explicit opt-in.

The program checks `ED2_CSR_RELEASE_PADDED_CACHE` is literally
`0` or `1`, rejects oversized/malformed values and use without
`ED2_CSR_COMPLEX_CACHE=1`, and inherits CSR's strict full-matrix
guard and index/allocation validations. The optional storage
release fails closed on deallocation errors.

## Memory accounting; important limits

The native source prints two machine-readable lines:

    Hamiltonian action: csr-complex cache entries=N, extra_bytes=M
    Hamiltonian action: csr-padded-cache released_bytes=R, net_saved_bytes=S

`R` is the exact allocator payload size of the two original
padded arrays; `M` is the exact payload of the CSR arrays;
`S = R-M` compares **post-construction allocated payload** to
the historical original padded representation. This is NOT a
direct measured peak RSS reduction, not even a guarantee of
lower physical RSS after freeing allocator pages.

Because CSR is first built from the complete padded arrays, the
**initial peak must still hold both representations**. The
optimization only reduces the long-lived cached-storage payload
after the one-time construction. No claim of lowering peak
RAM allocation or preventing initialization OOM is made.

For a representative **ten-site alternating S1/2, S1 OBC**
general XYZ + transverse-field chain (full basis **D=7776**),
real or nonreal DM:

| Data | Bytes |
|---|---:|
| Original padded complex + index arrays released | 14,463,360 |
| CSR complex arrays retained | 3,867,268 |
| Steady-state allocated cache payload saving | **10,596,092** |

Thus this generic example can reduce the long-lived payload
by ~10.60 MB (decimal), but not the construction peak. For a
high-occupancy matrix, CSR could have similar or even greater
payload than the padded representation, so no universal
memory-savings threshold is claimed.

## Scientific equivalence: not just the energy

A native standalone regression now exercises **three**
independent, byte-identical physical inputs for each case:
1. original cached padded `H*x` path;
2. complex CSR with padded arrays retained;
3. complex CSR with padded arrays released.

`tests/cpc_csr_release_regression.py` reads actual compiled
Fortran outputs from a native driver and compares **every
real/imaginary IEEE-754 bit** of **two deterministic complex
input vectors** across all three paths. It additionally
checks the native original arrays' allocation status,
actual nonzero imaginary matrix elements on DM models,
the exact release and retained bytes, and the original
independent dense Hamiltonian reference for N4/N6.
Four invalid configurations must fail closed: missing
complex CSR prerequisite, malformed/oversized release flag,
and incomplete Hamiltonian caching.

On a standalone local GNU Fortran 14.2 computer before
pushing the PR, **16/16 three-way exact bitstream pairs**
passed at N4/N6/N8/N10, both real and genuinely complex
DM, OpenMP 1/4. For each case the native dense-H and
memory accounting validation passed as applicable.
This is a correctness result, not PRESTO or a
journal-grade speed benchmark.

GNU13/14 GitHub science CI must independently pass
against the exact PR head before any default or
merge decision.

## Compiler floating-point optimization policy

ED2's historical default GNU compiler flags include `-Ofast`,
which permits floating-point reassociation/complex arithmetic
transforms. **Do not interpret raw IEEE bit equality between
the original padded multiplication loop and the separately
compiled CSR loop as unconditional under `-Ofast`.**

An independent local GNU14.2 experiment compared N4/N6/N8/N10,
both actual DM=0 and complex DM=0.12, T1/T4,
and two distinct nonreal complex vectors per condition:

- With `-O2 -fopenmp`: original, CSR retained, and CSR
  released were **16/16 identical full output bitstreams**.
- With `-Ofast -fopenmp`: CSR retained versus CSR released
  were **16/16 identical full output bitstreams**.
  Original padded versus CSR stayed bit-identical for
  DM=0 models. For the nonzero complex DM models,
  different compiler reassociation gave tiny per-element
  last-bit differences: maximum absolute difference
  `8.882e-16`, scaled by maximum output amplitude
  `2.816e-16` or less across the tested systems.
- Both compiled modes passed the independent direct
  dense-Hamiltonian reference at N4/N6.

The portable `tests/cpc_csr_release_regression.py`
therefore uses **two explicit policies**:

- Default strict `-O2` testing: *all three* H*x output
  bitstreams must be identical.
- `--allow-fp-reassociation` for `-Ofast` testing:
  **CSR retained versus CSR released must still match
  every IEEE bit**, while original versus CSR must match
  to `1e-12*max(1, max(|H*x|))` or better. The runner
  records both maximum numerical difference and whether
  the original-to-CSR pair was bit-identical.

The dual-compiler CI explicitly runs both optimization
policies on the exact source branch. This tolerance
does NOT justify approximate model/solver changes,
and no stronger reproducibility assertion is made
for all architectures or mathematical inputs.

## Compatibility and adoption limits

**Users with code that reads the public module variables
`ham2vec::wk_ele` or `ham2vec::wk_loc` after
`make_wk_loc_and_ele` must not opt into this release.**
Those allocatable arrays are deliberately unallocated in
release mode. Official native Lanczos matvec and
the existing full-H builder do not reference them in
the release path, but other user-written programs may.

Keep both flags explicitly optional until a separate
module API audit establishes that eliminating the
original arrays by default is safe. This work does NOT
alter BLAS/LAPACK integration, solver convergence,
symmetry, observables or the CPC manuscript.
