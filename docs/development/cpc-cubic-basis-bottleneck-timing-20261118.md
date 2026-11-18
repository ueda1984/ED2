# Cubic 1000-site basis bottleneck: 12-run PRESTO pilot and opt-in phase timing

## Already measured on PRESTO (2026-11-18)

PRESTO job `19673.jpmjpr1911`, pinned source `8608dc28866204aba1d243d3730b11805263be37`, N=1000 spin-1/2, three magnons, J1-J2 real XXZ, Gamma/O_h, and anchored representative cache `MNTE=64`, tested all four 1D irreps at 1, 4, 8 OpenMP threads.

User's independent trial reader reports `12/12 CONFIGURATIONS AND NUMERICAL CHECKS PASS`: unique configuration keys, `energy` within 1e-8 of 1-thread baseline within each irrep, direct Ritz residual ≤1.01e-8, correct dimension/sector, native anchored lookup enabled and `MNTE=64`. This is a one-shot diagnostic study with **no warmups or repeated timing statistics**. The site operator's independent `place=exclhost` attestation is **not supplied**.

| Sector | D | T1 basis | T4 basis | T8 basis | T1 cache | T4 cache | T8 cache | T1 total | T4 total | T8 total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A1g | 4299 | 12.758 | 6.228 | 5.202 | 1.467 | 0.412 | 0.220 | 14.336 | 6.694 | 5.465 |
| A2g | 3624 | 12.744 | 6.230 | 5.194 | 1.239 | 0.348 | 0.186 | 14.073 | 6.626 | 5.420 |
| A1u | 2884 | 12.740 | 6.215 | 5.192 | 0.985 | 0.276 | 0.148 | 13.814 | 6.543 | 5.384 |
| A2u | 3252 | 12.752 | 6.222 | 5.193 | 1.108 | 0.312 | 0.168 | 13.957 | 6.585 | 5.404 |

Across irreps the basis scaling is practically constant despite different projected dimensions. Basis performance improves 2.05x at 4 threads and 2.45x at 8 threads (30.6% basis parallel efficiency at 8 threads). In contrast, the cached operator construction for A1g improves 6.67x at 8 threads. The basis dominates 95%+ of total time at 8 threads.

An **inferred**, not directly measured Amdahl model `T_basis(p) ≈ S + P/p` gives `S≈4.05s` from p=1,4 and `S≈4.17s` from p=4,8 (A1g). A sequential group-validation stage is an obvious hypothesis: `state_lists::validate_noncommuting_action` is executed outside the OpenMP orbit loop and allocates/validates 48,000 affine permutations of 1,001 sites. However, timings alone cannot prove it causes the ~4.1s floor: 8-thread contention, memory bandwidth, sorting and imbalance can mimic serial time.

## Diagnostic-only Fortran instrumentation

New env setting `ED2_BASIS_DETAIL_TIMING=1` produces these lines in native ED2 `stdout.log` for a successful `try_orbit_basis` fast path:

```text
ED2_BASIS_PHASE preflight_other wall_s=<seconds>
ED2_BASIS_PHASE group_validation wall_s=<seconds>
ED2_BASIS_PHASE orbit_scan wall_s=<seconds>
ED2_BASIS_PHASE finalize wall_s=<seconds>
ED2_BASIS_PHASE total wall_s=<seconds>
```

Four mutually exclusive wall time intervals sum to the total. `preflight_other` includes generator checks and pre-enumeration overhead outside full normal-form validation. `group_validation` contains the existing normal-form/character verification. `orbit_scan` contains allocation, visited-bitset initialization and OpenMP orbit scanning. `finalize` contains merging and sorting/normalization. This timing does **not** bypass any validation, alter representatives, numerical formulas or OpenMP scheduling, or change output when disabled. The parent `ED2_PHASE basis` timer remains unchanged; differences on the order of sampling overhead are allowed.

CI regression `tests/cpc_cubic_basis_timing_regression.py` runs matched timed/untimed 4^3 and 6^3 A1g cases at threads 1 and 4, checking phase completeness, nonnegative finite values, wall time decomposition, projected dimension and Ritz energy/direct residual equivalence. Neither branch nor this PR changes the existing ED2 algorithms or PBS scheduler.

## After CI: operator-approved PRESTO diagnostic

Build a new source stage pinned to this PR head (do NOT replace any of the prior stages). On an operator-approved compute node, export `ED2_BASIS_DETAIL_TIMING=1` to the PBS environment; use single-valued `ED2_THREADS` and `ED2_SECTORS` for an inexpensive first N=1000 A1g diagnostic, then compare 1,4,8 threads in the **same** 12-run package if the first looks correct. Use `grep -a 'ED2_BASIS_PHASE'` on `measurements/trials/.../stdout.log`, and preserve the output and run metadata.

The only purpose is phase attribution: **do not yet optimize or remove general noncommuting symmetry validation** without a separate mathematically justified exact certificate and exhaustive regression. Future potential fast basis construction from translation-anchored triples must preserve the original representative ranks and normalization factors across all four Oh irreps.
