# Validated noncommuting-orbit fast path (development audit, 2026-11-18)

**Scope:** A bounded, verified *implementation-level* acceleration of pre-existing QS-Cube/ED2 site-permutation symmetry reduction. It is **not** the introduction of point-group projection, a new ED algorithm, or a generic multidimensional non-Abelian irrep solver. The current character-based basis is still limited to compatible **one-dimensional representations** of a momentum little group.

## Why the old fast path refused point groups

The OpenMP `try_orbit_basis` path uses one visited bit per basis state, so it can process each representative orbit once rather than testing every raw state against all operations. Previously it required all six cyclic generators to commute. In a square-lattice space group, translations and a rotation generally do not commute, despite a valid normal-form enumeration `T_x^a T_y^b C4^c M^d`. The old code therefore fell back to expensive exhaustive `checkstate_SQ`. Baseline 8x8/10x10 physical measurements are recorded in [cpc-pointgroup-performance-20261118.md](cpc-pointgroup-performance-20261118.md).

## Safety argument of the new bounded validation

**The commuting fast path and its legacy behavior are preserved.** Only for noncommuting generators, and before publishing any visited orbit bits:

1. Reject invalid maps/orders/local spin or occupation-cap mismatch using the existing prechecks.
2. Enumerate the proposed normal form (exactly the order used by `checkstate_SQ`) as site permutations. Reject duplicate permutations by integer fingerprint **plus exact element-wise comparison**. The hash is an index aid, not a probabilistic numerical acceptance test.
3. Check that the last cyclic normal-form element is the identity.
4. For every enumerated permutation `g` and each generator `s`, require `s*g` to appear in the same set. Since all maps are bijective and the set contains identity, generator closure on this finite set verifies the group. The verified group enumeration is exactly the set of orbit images used by the original basis routines.
5. Require the supplied complex phases to obey `chi(s*g)=chi(s)chi(g)` within `1e-10` for each generator and element, i.e. they define a consistent **one-dimensional** character. This rejects e.g. incorrectly treating X momentum as a full-C4v momentum sector.
6. If any check fails or validation resources exceed the explicit limit (8192 group elements or four million site-image entries), **do not** use the accelerated path; fall back to the original calculation, with unchanged numerical behavior. Rejection of acceleration **does not** imply that an invalid physical symmetry input has been repaired or certified.

The optional `ED2_FORCE_LEGACY_ORBITS=1` disables all atomic-orbit paths before basis construction, providing a precise same-input reference for regressions and performance. Its default is unset. This is an engineering opt-in for regression; do not permanently depend on or advertise an unreviewed diagnostic environment variable as a stable user API.

**Important boundary:** The new validation proves group/character consistency and site-label structure, but **does not prove `[H,U_g]=0` for arbitrary user-supplied Hamiltonian tables**. Physical regressions separately inspect each transformed bond table and verify zero leakage of the independently constructed projected Hamiltonian.

## Tests performed locally before push

A frozen source snapshot of the pre-optimization branch (commit `6284bb6ecd6140a0763eeff83bd82a07f17e9d64`) was compiled on Debian 13, GNU Fortran 14 with OpenMP.

- 31 physical translation/point-group cases were reproduced against independent bitmask Hamiltonians and character projectors. For all nonempty compatible 1D point-group sectors, a second **forced-legacy** ED2 run matched the optimized native result. Normal-form validation accepted known ring dihedral, square C4v Gamma/M and C2v X combinations.
- Six separate guard controls include two accepted dihedral characters and four intentionally invalid cases (wrong inversion momentum, full C4v at X, nonclosed normal form, and duplicate normal-form group elements). Every invalid fixture safely declined acceleration.
- The existing bit-exact OpenMP orbit regression compared the frozen baseline and new generator for 379 independent fixtures and five stress fixtures, over 1/2/4 threads with repeated stress schedules: **1167 candidate executions passed** using a bounds/undefined-behavior checked candidate build. Representative ranks and normalization matched the reference, including noncommuting ring cases.
- One local Debian/LAPACK-specific `ZHEEVR` workspace-query failure for a one-dimensional dense eigenspace required skipping the corresponding local dense solver cases; this is **not** counted as a successful local dense test. The untouched pre-optimization local executable showed the same issue, while the prior Ubuntu gfortran-13/14 CI tested those tiny cases successfully. This cross-platform discrepancy is separately actionable, but is not caused by the new orbit algorithm.

## Matched-input exploratory runtime (local Debian GNU Fortran 14; *not* publication evidence)

Three measurements per case, `ED2_LANCZOS_TWO_STAGE=1`, identical exact `input.dat` SHA256 and physical irrep, fixed seed, 1 or 4 OpenMP threads. Times are **medians** of opt-in `ED2_PHASE_TIMING` telemetry. All matched eigenvalues, dimensions and native Ritz residual certifications passed.

| Case | Threads | Legacy basis, ms | Validated fast basis, ms | Basis ratio | Legacy total, ms | Fast total, ms |
| :-- | --: | --: | --: | --: | --: | --: |
| 8x8 Gamma-A1 (D=109) | 1 | 9.038 | 2.193 | 4.12x | 21.98 | 15.55 |
| 8x8 Gamma-A1 (D=109) | 4 | 3.618 | 2.598 | 1.39x | 9.58 | 8.38 |
| 10x10 Gamma-A1 (D=247) | 1 | 42.802 | 6.545 | 6.54x | 84.07 | 50.63 |
| 10x10 Gamma-A1 (D=247) | 4 | 15.503 | 7.095 | 2.19x | 31.60 | 22.88 |
| 10x10 Gamma-B1 (D=218) | 1 | 42.909 | 6.354 | 6.75x | 82.91 | 44.78 |
| 10x10 X-A1 (D=465) | 1 | 47.191 | 4.308 | 10.95x | 91.88 | 52.88 |

For larger thread counts, the up-front single-threaded group validator can limit basis-generation gains. The new code should therefore remain a development experiment until independently reviewed and compared on controlled HPC resources.

## Replication commands

```sh
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j2
make check
make check-pointgroup
python3 tests/cpc_noncommuting_orbit_guard.py --source . --compiler gfortran
python3 tests/cpc_pointgroup_profile.py --exe ./source/QS3ED2 --sizes 8,10 --threads 1 --modes two-stage
ED2_FORCE_LEGACY_ORBITS=1 python3 tests/cpc_pointgroup_profile.py --exe ./source/QS3ED2 --sizes 8,10 --threads 1 --modes two-stage
```

The final comparison requires distinct output directories or the default independent temporary directories. The pinned-baseline OpenMP basis regression is run in CI with the old source reconstructed by commit SHA.

## Review and publication posture

No merge to main, no new default behavior declaration, and no CPC claim of new symmetry methods until independent source review, gfortran-13/14 CI, zero-cap and other stress fixtures, and exclusive PRESTO performance attribution are complete. As in the separate two-stage Lanczos restoration, direct residuals certify returned Ritz pairs but not global lowest-state inclusion. Distinguish demonstration of physical sectors, engineering speedup, and the provenance of published QS³ algorithms.
