# ED2 incremental performance audit (2026-11-16)

## Status and decision

Two materially useful candidates remain AFTER the previously validated full/tight Hamiltonian cache: avoid rediscovering symmetry orbits during basis construction, and avoid unused Ritz-vector construction/storage plus premature full-subspace closure checks in ALG2. Working prototypes have been tested. This is an isolated REVIEW PATCH, not a change to main, upstream, or the manuscript branch.

Baseline source: `719da7035e66ce884e9de00376e4ded5b3e5c459`.
Review branch: `codex/ed2-cpc-performance-audit-20261116`.
Candidate: `tests/performance_audit_20261116.patch`, Git blob `8f380d68140d3333d98701d2a00ef64914ad0923`.
It changes only `source/state_lists.f90` and `source/krylov_solvers.f90`.
The test-only timing instrumentation of main.f90 is NOT part of the patch.

## Incremental performance, not a repeat of the cache ON/OFF result

Both versions use FULL, model-specific tight-row cache, identical inputs and explicit initial vectors, complex double precision and absolute residual threshold 1e-8. Five paired/interleaved measured executions per workload, one CPU-affinity core, GNU Fortran 14.2, -O2 -fopenmp, system LAPACK/BLAS, OpenMP=BLAS=1. The machine is shared: these are engineering measurements, not a controlled exclusive-node CPC performance benchmark.

| Workload | Reduced D | Baseline median (s) | Candidate median (s) | Ratio |
| --- | ---: | ---: | ---: | ---: |
| Cubic 1000 sites, two magnons, Gamma, ALG1 | 503 | 0.965905 | 0.057255 | 16.870 |
| Cubic 216 sites, three magnons, Gamma, ALG1 | 7677 | 1.343560 | 0.490311 | 2.740 |
| Non-U1 mixed-spin DM chain, 12 sites, ALG1 | 46656 | 0.958106 | 0.754532 | 1.270 |
| Same DM chain, ALG2 keep=2, cap=32 | 46656 | 6.939522 | 5.359263 | 1.295 |

The corresponding total-time median absolute deviations (baseline/candidate seconds) are 0.007734/0.000639, 0.030341/0.013646, 0.019823/0.013552 and 0.217156/0.138511. All 40 measured runs passed residual acceptance. Within each workload the returned energy, input hashes and Hamiltonian-action counts match exactly between versions. H-action counts are 165, 520, 140 and 363, respectively. Do not multiply these ratios by the older no-cache/full-cache ratios obtained in a different environment.

## Candidate 1: basis construction by marking already visited orbits

The original full-cache 1000-site/two-magnon workload spends about 94% of its total time in basis construction. Median basis time falls from 0.907469 s to 0.004788 s, while cache construction and solver time remain essentially unchanged. For 216 sites/three magnons the corresponding basis times are 0.875328 s and 0.017884 s.

The prototype marks translated raw ranks in a bitset, enumerating each orbit once. It preserves the ascending representative rank ordering, the original phase-sum order and stabilizer normalization. Preconditions check permutations, generator periods, local spin/cap preservation and commutation. Unsupported maps, more than one OpenMP thread, or a visited bitset above 128 MiB use the old path. Identity symmetry bypasses orbit work entirely and also supports multiple threads. This deliberately guarded prototype is NOT a new parallel strong-scaling implementation.

At 1000 sites/three magnons the temporary bitset is approximately 19.81 MiB and is freed before cache/solver setup. This is a bounded setup-memory tradeoff for a substantial speed improvement, not a universal memory saving. In a checked native ED2 basis-only run, the candidate produced raw dimension 166167000 and Gamma dimension 166167 in 5.88 s, with 26908 KiB peak process RSS. No native 1000-site/three-magnon converged eigenenergy is claimed.

## Candidate 2: eliminate unused restarted-solver vectors and checks

The former ALG2 allocates `ritz(n,cap)` and forms all cap Ritz vectors, although only wanted/retained vectors are consumed. The candidate allocates `ritz(n,retained)` and forms only `min(retained,m)` columns. It computes the expensive full-subspace closure test only when Ritz residuals are certified or the iteration budget is exhausted, without weakening acceptance or altering final diagnostics.

For n=46656, retained=2 and cap=32, the removed allocation is 16*46656*(32-2)=22394880 bytes (21.36 MiB). Measured median peak RSS falls from 114044 KiB (111.37 MiB) to 92048 KiB (89.89 MiB), a 19.29% reduction, with total time improved by 1.295x. Tiny ALG1 RSS fluctuations are NOT counted as a memory improvement. No reduced precision, fast-math, relaxed tolerance or approximation is introduced.

## Separate correctness issues found

1. Unused combination-table column: requesting 1000 sites/three magnons computes a needless fourth-magnon column in the original source. C(1000,3)=166167000 fits signed int32, but C(1000,4)=41417124750 does not. A baseline count-only UBSan build reproduced signed integer overflow at state_lists.f90:52. The patch computes only needed columns and uses checked int64 intermediate arithmetic before storing the existing integer index type. It does not claim support for genuinely oversized index spaces or change the BLAS/LAPACK ABI.

2. Local zero occupation caps: plateaus in cumulative counts require the FIRST count >= rank. The old arbitrary-equal binary search produces a four-site normalization error and, in a six-site case, dimension 7 instead of independently enumerated dimension 9. The lower-bound fix agrees with independent product-state enumeration for both. These are intentional corrections, not legacy byte-identical results.

Official compiler context: https://gcc.gnu.org/onlinedocs/gfortran/Behavior-on-integer-overflow.html . Fortran overflow is undefined, not a reliable wraparound strategy.

## Validation

- 379 symmetry/occupation cases at OMP=1 and the SAME 379 at OMP=2. Mixed spins, all ring momenta, occupation windows, redundant commuting maps, noncommuting fallback, empty sectors, identity and zero caps are included. 377 regular cases retain legacy behavior; two zero-cap fixes agree with independently enumerated dimensions and normalization. Do not count these as 758 different Hamiltonians.
- 418 adversarial Krylov cases match all non-timing numerical fields locally: 164 strict accepts, 132 accepted-but-accuracy/inclusion-warning cases, and 122 nonconvergences, unchanged. This does NOT mean 418 successful convergences.
- Independently enumerated bounded-space spectra: 72 levels across three spaces.
- Literature XYZ-DM Hamiltonian: six finite-size cases with independent spectra/observables.
- Stationary non-U1 mixed-spin scans: 18 cases at one thread and the same 18 at two threads, with independent energies, Cartesian nearest-neighbor correlations and chirality.
- GNU Fortran 13/OpenBLAS CI independently verified the exact patched source hashes, release and bounds/UBSan builds, the independent physics above, large-site full-cache references, and the 418-case outcome/iteration classifications. Both jobs succeeded:
  https://github.com/ueda1984/ED2/actions/runs/43789877223 .

## Incorporation recommendation and remaining limits

Recommend incorporating the checked counting/lower-bound fixes, the restarted-solver work reduction, and the guarded orbit path before freezing a performance baseline. Before a production merge, keep the one-thread restriction explicit, review the parallel basis strategy and add all new boundary cases to permanent regression CI. Exclusive-node repeated benchmarks remain necessary for paper-level speedup claims.

Do not add CSR conversion, generic duplicate-transition aggregation, precision reductions or other complexity solely for marginal improvements. The existing generic-vs-tight MNTE allocation difference is already avoidable through the validated configuration and is not a new optimization discovered by this audit.

The original manuscript branch, main, upstream and release metadata remain unchanged. Only review artifacts and a CI that applies the patch inside disposable checkouts were added to this branch.
