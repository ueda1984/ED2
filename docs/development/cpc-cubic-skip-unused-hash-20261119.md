# Elide redundant keyed-orbit hash in exact cubic direct-hop mode (PR #21)

**Status: experimental opt-in; no change to PR #19, #17, #16, #20, or main.**

## Bottleneck and narrow purpose

At N=1000, down spins=3, cubic Gamma/O_h sector A1g on
the validated PR19 basis-pruned path, the optional
`ED2_CUBIC_CACHE_PROFILE=1` timer reported about
44 ms in the representative-key setup and about
97 ms in the total Hamiltonian cache build in one
GNU Fortran 14.2 Linux execution. These are
**shared-host exploratory diagnostics**, not a PRESTO performance
result or an independent timing guarantee.

In `source/cubic_anchored_lookup.f90`,
`setup_cubic_anchored_lookup` has historically:
1. Enumerated canonical orbit keys for each of `lv`
   symmetry representative triples using **48 point actions and
   three translation anchors** (and built an open-addressed
   keyed hash lookup).
2. Separately built `direct_indices` and
   `direct_characters` of length `(NOS-1)(NOS-2)/2`
   with the **same 48 x 3** enumeration.

When `ED2_CUBIC_DIRECT_HOP_TABLE=1`, the runtime resolver
always takes the second, dense **direct** branch; the keyed
hash is never read. Its first setup phase is therefore
redundant in the fully certified cubic direct-hop configuration.

## Implementation

The new `ED2_CUBIC_SKIP_UNUSED_HASH=1` flag is **OFF by default**
and fails closed if the direct-hop table is not active. The
existing anchored-cubic geometry and model guards, point-group
map validation, 1D character checks and **dense direct table's
duplicate representative/character consistency proof** remain.
The opt-in simply skips the separate hash construction and
prevents its unnecessary storage; it still builds the
same dense pair table and sets `active/direct_active` to
the required values.

The dense builder verifies every duplicate mapping and
1D character, so its stabilizer consistency checks remain
in force. The historical hash-builder's own conflict
check is not executed in this mode because the hash
and its corresponding search are not used. All legacy and
general Hamiltonians retain their exact former paths.

A native output marker identifies activation:

    Representative lookup: cubic-direct-only-unused-hash-elided

Invalid values and missing `ED2_CUBIC_DIRECT_HOP_TABLE=1`
are explicitly rejected.

## Exactness proof

`tests/cpc_cubic_skip_unused_hash_regression.py` runs the **real**
compiled ED2 Hamiltonian cache driver with coordinate caching
and PR19 axis pruning both ON and unused-hash elision OFF/ON:
- Four one-dimensional cubic irreps × L4/L6/L10 ×
  OpenMP T1/T4, plus L10 T8, 28 **complete native
  Hamiltonian cache SHA256 pairs per compiler**.
  This checks every target, diagonal and complex coefficient,
  including ordering and all signs; not just eigenvalues.
- L4/L6 × four irreps × T1, eight direct Ritz
  energy/residual/H-application comparisons per compiler.
- Negative cases reject the new flag without direct-hop
  and reject malformed flag values.
- All prior PR13–19 exact cubic basis, branch-and-bound
  and complete H cache regression tests run unchanged.

The change never modifies basis generation, canonical ranks,
norms, eigensolver or physical model. The existing complete
basis proof remains required in the CI matrix.

## Exploratory measurement and next gate

A separate **shared Debian GNU14.2 Linux** sample with
N1000/A1g, fixed flags, six interleaved OFF/ON repetitions
(one initial sample excluded) gave:
- representative lookup setup: median **41.863 vs 22.385 ms**
  (difference 19.478 ms);
- operator cache wall: median **71.418 vs 44.535 ms**;
- end-to-end: median **247.997 vs 236.590 ms**,
  but with high inter-trial variation and other phases not
  stable; this **does not establish** a systematic end-to-end
  performance gain.

Before selecting this optional feature for a release or paper:
1. Both GNU Fortran 13/14 full native checks and negative
   provenance tests must finish successful on the final head.
2. Conduct one same-commit PRESTO OFF/ON performance study with
   PR19 basis pruning ON in both conditions, changing
   ONLY `ED2_CUBIC_SKIP_UNUSED_HASH`.
3. Preserve all native trial/audit records, scheduler placement,
   compiler and source SHA; collect independent operator
   confirmation of exclusive cores for journal-grade timing.
4. If confirmed, consider reducing the next substantial
   common cost in the operator or Lanczos phase,
   not repeatedly micro-optimizing one-time basis setup.

**Do not merge or tag without author approval.**
