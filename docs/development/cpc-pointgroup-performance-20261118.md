# Combined translation / point-group physical profiling (2026-11-18)

**Status:** validated diagnostic benchmark on GitHub-hosted shared hardware. This is not a PRESTO performance result or publication-grade speedup measurement.

## Motivation

The historical QS-Cube/ED2 cyclic-permutation interface can project compatible one-dimensional little-group irreps. The current optimized `try_orbit_basis` builder handles commuting generator collections, but for noncommuting translation + rotation/reflection generators it falls back to an exhaustive representative test. We measured the *whole workflow* as well as individual phases to quantify this tradeoff; restoring two-stage Lanczos is a separate historically established solver-continuity effort.

## Instrumentation and protocol

- Opt-in environment: `ED2_PHASE_TIMING=1`. Otherwise there are no new phase records, and no numerical solver defaults are changed.
- `source/main.f90` emits elapsed wall time from Fortran `system_clock` for input, combinatorics, symmetry-basis construction, operator-cache generation, ALG=1 Lanczos and the total.
- `tests/cpc_pointgroup_profile.py` uses only Python's standard library and independently verifies group-character dimensions with a fixed-magnon cycle-count Burnside projector.
- Fixed square-lattice periodic spin-1/2 J1-J2 XXZ model: NN `Jxy=1,Jz=0.77`; both NNN diagonals `Jxy=0.31,Jz=0.19`; three magnons. Sectors are `Gamma` translation-only, `Gamma A1`, `Gamma B1`, and `X A1`. `C4v` 1D characters are valid at Gamma; `C2v` little-group characters at X.
- One thread, `ED2_LANCZOS_TWO_STAGE=1`, fixed RNG seed, fully cached Hamiltonian and requested direct residual tolerance `1e-8`. Binary, input, and diagnostic SHA256s are archived.
- The formal 31-case physical regression now **requires** actual certified Ritz rows from the diagnostics table, rather than accidentally ignoring residuals due to a non-existent `pair_1_residual` key. Both gfortran-13 and gfortran-14 passed that regression and 16 native profiling runs in [run 44216235216](https://github.com/ueda1984/ED2/actions/runs/44216235216).
- Large-sector measurements below come from [run 44216301856](https://github.com/ueda1984/ED2/actions/runs/44216301856) (gfortran-14, Ubuntu hosted, one unreplicated run per sector). Their archived `report.json` contains exact phase timings, energy, direct residual, memory and Hamiltonian-application counts.

## Large-sector measurements (shared-runner diagnostics, not controlled benchmarks)

| Lattice / sector | Hilbert dimension | Basis, ms | Cache, ms | Lanczos, ms | End to end, ms | Peak RSS, MiB | H applications |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 8x8 / Gamma-T | 651 | 0.655 | 17.450 | 9.233 | 29.396 | 30.4 | 270 |
| 8x8 / Gamma-A1 | 109 | 8.207 | 9.636 | 1.236 | 19.542 | 9.0 | 150 |
| 8x8 / Gamma-B1 | 91 | 8.220 | 8.286 | 0.805 | 17.800 | 8.4 | 110 |
| 8x8 / X-A1 | 200 | 8.173 | 10.424 | 1.489 | 20.564 | 12.6 | 120 |
| 10x10 / Gamma-T | 1,617 | 2.376 | 61.049 | 58.230 | 122.256 | 104.0 | 450 |
| 10x10 / Gamma-A1 | 247 | 38.808 | 32.809 | 3.804 | 76.110 | 19.9 | 250 |
| 10x10 / Gamma-B1 | 218 | 38.701 | 28.898 | 2.522 | 70.818 | 18.2 | 190 |
| 10x10 / X-A1 | 465 | 41.341 | 36.158 | 4.103 | 82.288 | 33.2 | 160 |

Every recorded run passed the native direct-residual certification (residuals approximately `1e-9` to `9e-9`). `Gamma-A1` has the same lowest energy as `Gamma-T` for these specific models, but other irreps and X differ, so their numerical difficulty is not identical. Do not compare Hamiltonian-application counts across unlike sectors as proof of an isolated solver speedup.

**Interpretation:** on the shared CI runner, the point-group reduction lowered *total* runtime and memory for the tested lattices even though noncommuting point-group basis formation grew to roughly 13x (8x8) and 16x (10x10) the translation-only basis time. At 10x10 the `Gamma-A1` sector reduced dimension `1617 -> 247`, but the fallback consumed roughly `39 ms` of a `76 ms` program. This is a concrete reason to evaluate a validated noncommuting orbit builder, not evidence of a new physical algorithm.

## Later experimental implementation

The bounded group-normal-form and character validation described in
[cpc-noncommuting-orbit-validation-20261118.md](cpc-noncommuting-orbit-validation-20261118.md)
has now been implemented on a **separate draft branch** and locally
compared against a forced-legacy basis. The source below remains a
baseline performance record; it does not represent the newer fast path.

## Next engineering experiment

Propose a **separate** opt-in/auditable improvement to `try_orbit_basis` with full enumeration of generator normal forms, exact duplicate-permutation detection, closure under multiplication by each generator, and a verified one-dimensional character homomorphism. **Do not** simply remove the commuting guard without those checks. Retain `ED2_FORCE_LEGACY_ORBITS=1` as a same-Hamiltonian/same-irrep baseline and test both builders against the independent fixed-magnon physical projector (eigenvalues, dimensions, eigenvector residuals) and cross-thread reproducibility.

Publication-level performance evidence still requires a pinned exclusive PRESTO compute node, compiler/affinity/cache configuration, multiple repetitions, stability checks, and independent attestation. No PR should merge to main before the gated numerical tests and source review are satisfactory.
