# Exact periodic-coordinate cache for cubic three-magnon anchored basis

Status: experimental Draft PR stacked directly on PR #16 (`9692d8b1fd903bcb7b1b67a5f180c57430702510`).
This changes neither `main` nor the pinned PR #16 source.

## Why

PRESTO N1000 three-down-spin Gamma/O_h 4 irreps × T1/T4/T8 independently
reported `ED2_BASIS_PHASE orbit_scan` of 182.6–190.7 ms (99.1–99.3%
of anchored basis time). Startup basis is now approximately 67–71% of
the T8 end-to-end runtime following the verified PR #16 direct-hop path.
These measurements are one-shot exploratory phase records, **not**
publication-grade benchmarks.

The anchored basis examines at most 48 point-group images and 27
axis-wise zero-anchored translations for each affine orbit. The old loop
recalculates three periodic modulo differences for each of three spins
at every translation candidate. It also computes the same differences
for three translation-origin views. For each image, the old path may
evaluate 27×3×3 + 3×3×3 = 270 periodic differences.

## Transformation and exactness

Explicit opt-in `ED2_CUBIC_BASIS_COORD_CACHE=1` (default OFF) computes
`delta_xyz(mu,j,i)=modulo(xyz(mu,j)-xyz(mu,i),L)` for all three sites
and each axis once per O_h point image (27 periodic differences).
The existing anchor-key marking and canonical minimization loops then
compose integer site indices using these stored differences. In both
paths, the 27 candidate enumeration order, insertion sort, strict
colexicographic minimum updates, orbit-key marking, stabilizer operation
order and complex floating-point normalization are untouched.

The flag refuses an absent `ED2_CUBIC_ANCHORED_BASIS=1`, legacy
orbit-force mode, or malformed flag values. This path inherits
the anchored basis canonical cubic Gamma/O_h geometry certificate and
three-down-spin guard. The original code remains executable with the
flag unset/0. A native provenance message distinguishes enabled runs.

## Verification

`tests/cpc_cubic_basis_coordinate_cache_regression.py` runs the real
Fortran basis-only driver on byte-identical native input decks.
It compares the SHA256 digest of every (sorted representative rank,
IEEE-754 norm bits) dump, verifies the independent expected sector
dimension and checks the new provenance marker. The CI matrix exercises
all four real 1D irreps, sides 4/6/8/10, threads 1/4 and side 10
threads 8, plus explicit negative prerequisite rejection. In addition, `tests/cpc_cubic_basis_coordinate_full_cache_regression.py`
constructs the complete native Hamiltonian cache twice with the exact
PR16 direct-hop code and coordinate flag OFF/ON, comparing the SHA256 of
every diagonal/offdiagonal bit, native target and record order for side
4,6 at T1/4 and side10 at T1/4/8, all four irreps. At side4,6 T1 it
also compares complete native Ritz energy and both direct residuals.
Existing original-versus-anchored basis and direct-hop regression
oracles remain in CI. The entire suite runs for gfortran 13 and 14.

**Do not infer PASS from the presence of the test:** wait for both
compiler jobs to finish successfully and inspect all individual
comparisons. Independently repeat PRESTO N1000 unprofiled, same-input
measurements and record basis, cache, end-to-end wall, and RSS before
claiming a measured speedup. Independent operator attestation and
warmups/replicates are still needed for publication-level timing.
No merge into main without explicit review.
