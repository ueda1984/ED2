# Phase 1-D: independent Krylov stress and resource characterization

Date: 2026-11-06 (Asia/Tokyo). Required and verified starting HEAD:
`6d3480c5fb8ee7a67f2ac4017156ee4368bdd1d2`.
Branch: `codex/ed2-phase1d-krylov-stress-benchmarks`.

The executor was restored before work. Shell execution, a clean checkout, the
existing remote branch and its exact HEAD were verified before edits. Read
`AGENTS.md`, the [Phase 0 audit](baseline-audit.md), [CPC plan](cpc-library-plan.md),
[Phase 1-A](phase1-lapack-interfaces.md), [Phase 1-B](phase1b-sparse-field-indexing.md)
and [Phase 1-C](phase1c-lanczos-stability.md).

No production Fortran, numerical policy, Hamiltonian, symmetry, field indexing,
LAPACK interface, reference data or existing test assertion is changed. This
phase identifies limitations and supplies a [separate policy proposal](phase1d-policy-proposals.md).
It does not certify the solver as ready for extraction or publication.

## Main findings

1. A small residual is a backward-error statement about an eigenpair. It does not
   certify global ground-state inclusion. A deterministic diagonal example returns
   an exact excited eigenpair with zero residual. Giving that start a nonzero but
   tiny ground component also permits numerical near-breakdown before that component
   is resolved. Both observations reproduce at one and two OpenMP threads.
2. For gaps around 1e-10, the residual allowance proportional to
   `sqrt(lnc_ene_conv)` can accept unresolved mixtures. Ground-energy errors can be
   near 2e-10 while requested individual-state/projector errors approach one.
   The larger isolated cluster can nevertheless be accurately represented. These
   are different scientific questions, not contradictory tests.
3. At tolerances below floating-point resolution, the unchanged energy-change
   gate can reject nearly machine-precision ALG=2 candidates for the entire budget.
   Independent last-candidate residuals are retained; they are **not** relabeled
   successful solves. A small residual alone is not the current ordinary exit policy.
4. The expanded GNU 13 regression run exposed one stochastic sparse-field ALG=1
   observable failure: successful exit, residual 2.644543803251125e-8, and maximum
   magnetization error 1.53227552e-8 against the unchanged 1e-9 absolute/relative
   assertions. Energy agreement and the near-zero variance diagnostic do not
   guarantee observable accuracy. This is not established as compiler-specific.

## Test design and independent references

`tests/krylov_stress_driver.f90` calls the exact production internal kernels via
a dense callback; it counts every H application, records solver time, status,
reason and RNG states, and exports successful vectors. For failed ALG=2 solves it
also exports the last candidate under separately named test files, allowing actual
independent residual checks without changing the production failure contract.
ALG=1 does not expose its failed candidate vector, so its failed residual is only
the reported kernel value, explicitly distinguished in the data.

`tests/krylov_stress.py` runs 418 checks per profile:

- Dimensions 32/96; separated, fourfold degenerate, nearly degenerate, zero-cluster
  and zero-centered spectra. Independent complex Gaussian QR matrices (seed
  `8191+n`) rotate analytic spectra; the analytic eigenvalues and QR columns,
  rather than the ED2 projection/recurrence, define the references.
- Both solvers, energy-change tolerances 1e-6/1e-12/1e-20, three starts and both
  OpenMP thread counts. ALG=2 normally requests three vectors with NOK=6/NOM=12.
- A parameter grid at dimension 96 varies wanted=1/3/5, NOK=6/10, NOM=12/24/48,
  and cycle budgets 1/5/120; ALG=1 uses at most the physical dimension.
  This is a kernel grid: capacity can be less than twice the retained count;
  the existing CLI wrapper would raise NOM to at least twice NOK before solving.
- Two analytic adversarial starts distinguish excluded/weak ground content from
  actual residual validity; both happy and numerical near-breakdown are reported.

Every successful solve is assessed with independently computed `H*Q-Q*Lambda`,
normalization/orthogonality, nearest eigenvalues, lowest requested eigenvalues,
requested-space leakage, and applicable complete-subspace projector errors.
Partial exact degeneracies are assessed by membership, not a unique-vector claim.
For nearly degenerate spectra the requested cutoff uses 1e-12 times spectral scale;
the broader cluster uses 1e-8. The differing cutoffs are explicit in the script.

The `accepted_strict` label means ground-energy error <=1e-10 times
`max(||H||,1)` and requested-space leakage <1e-5. Actual residuals remain separate
columns. `accepted_accuracy_or_inclusion_warning` means the solver returned success
but missed either assessment target, not necessarily that its stated stopping
policy was violated. The loose `valid_pair_under_loose_assessment` flag checks
normalized residual <2e-3 and orthogonality <1e-8; it is not a high-accuracy certificate.
No numerical finding is converted into an expected-failure annotation.

The harness returns nonzero for execution/protocol errors. Legitimate solver
nonconvergence and numerical accuracy warnings are retained in its report and
summary counts even when the characterization command itself exits zero.

## Observed stopping behavior

| Profile | Strict assessment | Accepted warning | Nonconvergence | Total |
| --- | ---: | ---: | ---: | ---: |
| GNU 14 default `-Ofast` | 163 | 134 | 121 | 418 |
| GNU 14 `-O2` | 164 | 135 | 119 | 418 |
| GNU 14 bounds checked | 164 | 135 | 119 | 418 |
| GNU 13 `-O2` | 164 | 135 | 119 | 418 |

No harness/runtime-error or invalid-parameter classification occurred in these
sweeps. The finite-budget grid intentionally includes insufficient budgets; a
nonconvergence count is not an assertion that the method should have succeeded.

Representative O2 cases (IDs index deterministic fixture order):

- `0129-omp1`: dimension-96 clustered ALG=1, tol=1e-12; ordinary energy convergence
  at step 21, ground-energy error about 2.18e-10, ground-vector projector error
  about 0.967. Other seeds/thread counts stop at steps 20/21 and have projector
  errors around 0.92–0.99. An unresolved cluster mixture is not a resolved ground vector.
- `0207-omp1`: diagonal H with eigenvalues `[-1,-1+1e-10,0,2]`, start equal to
  the second basis vector. Exact happy breakdown returns that excited eigenpair;
  residual is zero and ground-projector error is one.
- `0208-omp1`: same problem with start `[1e-8,1,0,0]`; numerical near-breakdown
  returns an almost exact excited state, residual about 1.11e-16, ground-projector
  error one. Nonzero overlap alone gives no practical inclusion guarantee.
- `0015-omp1`: separated ALG=2, tol=1e-20; after 120 cycles, the last three
  candidates have independently checked residuals about 1.7–2.2e-15 and eigenvalue
  errors <=2.3e-16, but the energy gate still reports nonconvergence. Across O2,
  bounds and GNU 13, 63 unsuccessful ALG=2 cases meet the separate last-candidate
  assessment of energy error <1e-10*||H|| and residual <1e-8*||H|| (64 with Ofast).
  That assessment is looser than some requested tolerances; it does not prove
  that every such failure is unnecessary under the existing contract.

Away from zero, the legacy relative energy-change test is retained. Near zero,
Phase 1-C uses observed-operator-scale absolute energy change. The current direct
residual bound is `max(256*epsilon,sqrt(tol))*max(sigma,abs(E),tiny)`.
A square-root bound can be consistent with a second-order energy error while
allowing much larger vector/observable errors; the gap matters. Tolerances such
as 1e-20 cannot promise 1e-20 relative energy accuracy in binary64. This phase
neither tightens the bound nor adds a roundoff floor to the energy gate.

## Representative physical Hamiltonians and resource measurements

`tests/krylov_benchmarks.py` generates open spin-1/2 XYZ+DM+Gamma field chains,
J1–J2 chains, and Heisenberg chains with nonuniform z fields. All use the identity
symmetry sector; no invalid translation symmetry is imposed on nonuniform fields.
Dimensions are 2**N for N=6/8/10/12/14. ALG=1 requests one pair; ALG=2 requests two,
with additional wanted=4/NOK=8 and NOM=12/16/32/64 configurations. Cache-off,
partial and full cases are included.

Independent bit-Pauli actions validate every returned vector without storing a
large dense Hamiltonian. For N<=8, independent dense spectra and subspaces also
check ground-state inclusion. The bit operator was elementwise cross-checked
against the separate existing tensor-Pauli constructor for all three model families
at N=4 (maximum difference zero in this run). For N>8, small residuals plus old/new
agreement are corroboration, **not** independent proof of global lowest states.

The comparison builds exact Phase 1-B and Phase 1-C archives, then instruments
only external copies. The saved patches (losslessly encoded as JSON lines) add seed control, H-call counting and
measurement/basis export. `tests/krylov_benchmark_probe.f90` initializes all OpenMP
thread RNG states and records them. No recurrence, orthogonalization, stopping
policy or operator arithmetic is changed. Both copies use GNU 14 O2 and the same
LP64 OpenBLAS backend. Algorithms and requested accuracy remain different enough
that time ratios alone are not efficiency claims at equal achieved accuracy.

Runs are serial after regression/stress activity finished. Each configuration is
run at OMP=1/2, BLAS=1, `OMP_DYNAMIC=FALSE`, `OMP_PROC_BIND=close`, `OMP_PLACES=cores`.
Three starts (19/73/211) are paired between versions; version order alternates.
The ranges combine start dependence and host/timing variability, not controlled
fixed-work timing confidence intervals. There is no CPU-frequency or exclusive-host
control; this is a cloud characterization, not a publication performance claim.
Thread-specific RNG streams also mean OMP=1 and OMP=2 can do different Krylov
work; their time ratio is not a fixed-start strong-scaling measurement.

Whole-process wall time includes setup, cache, solve, output and basis export.
`compute_seconds` starts after input/seed initialization and ends after output,
before basis export. H-call counts include reconstruction, residual/quality checks
and all actual operator calls, rather than the historical printed step estimate.
Cycle counts are the number of printed `tritr` records and must be interpreted
with each version's initialization/restart semantics. Failed processes retain
status/diagnostics and are excluded from successful-solve timing interpretations.

Peak RSS is the executable's `/proc/self/status` `VmHWM` in KiB at completion,
which retains solver allocation peaks after deallocation. The independent `wait4`
value is also retained, but can include the Python launcher's pre-exec memory
high-water mark and is not used as the primary memory comparison. Probe runs on
small problems demonstrated this distinction (about 5.5 MiB versus 36 MiB).
Tiny-run wall times include process/polling overhead and are not useful speedups.

All **192 benchmark solves exited successfully**. The 48 N<=8 runs agree with independent dense references; the other 144 have valid measured residuals under the stated assessment, without an independent global-inclusion certificate. There were no benchmark solver failures.

Selected medians over the three paired starts (B = Phase 1-B, C = Phase 1-C). Residual columns show the maximum actual norm across those starts and requested pairs.

| Case | OMP | Wall s B/C | H calls B/C | Cycles B/C | RSS MiB B/C | Max residual B/C |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| XYZ N=14 ALG=1 | 1 | 0.491/0.489 | 113/113 | 0/0 | 47.16/49.17 | 1.98e-06/1.98e-06 |
| XYZ N=14 ALG=1 | 2 | 0.292/0.285 | 117/117 | 0/0 | 47.23/49.29 | 2.82e-06/2.82e-06 |
| XYZ N=14 ALG=2 NOM=16 | 1 | 0.578/3.237 | 138/792 | 11/44 | 52.43/60.21 | 8.22e-15/2.65e-11 |
| XYZ N=14 ALG=2 NOM=16 | 2 | 0.340/2.187 | 138/774 | 11/43 | 52.55/60.19 | 1.22e-14/6.05e-15 |
| XYZ N=14 ALG=2 NOM=64 | 1 | 1.012/3.420 | 246/396 | 4/6 | 64.93/96.68 | 8.58e-15/7.61e-15 |
| XYZ N=14 ALG=2 NOM=64 | 2 | 0.630/2.829 | 246/396 | 4/6 | 64.88/96.60 | 5.91e-15/6.82e-15 |
| J1–J2 N=10 wanted=4 | 1 | 0.061/0.584 | 228/1908 | 9/53 | 9.18/9.91 | 5.47e-15/4.88e-15 |
| J1–J2 N=10 wanted=4 | 2 | 0.038/0.406 | 204/1800 | 8/50 | 9.15/9.97 | 5.47e-15/4.53e-15 |
| HB N=10 cache off | 1 | 0.158/0.885 | 114/648 | 9/36 | 5.73/5.91 | 9.10e-15/3.54e-15 |
| HB N=10 cache off | 2 | 0.094/0.483 | 126/648 | 10/36 | 5.98/5.95 | 6.11e-15/5.35e-15 |

At N=14/NOM=16/OMP=1, revised ALG=2 takes about 5.6 times the median wall time and 5.7 times the median H applications of Phase 1-B on this model. This sample does not support a general efficiency improvement. ALG=1 timings/H counts are similar here; its reconstructed-vector peak RSS is about 2 MiB higher. ALG=2 NOM=64 reduces cycle/H counts relative to NOM=16 in the revised solver, but increases projection/orthogonalization cost and RSS; the median whole-process time is slightly larger, not smaller.

Dimension scaling for ALG=2/NOM=16/NOK=4/wanted=2, full cache, OMP=1:

| N / dimension | Wall s B/C | RSS MiB B/C |
| --- | ---: | ---: |
| 6 / 64 | 0.007/0.011 | 5.59/5.42 |
| 8 / 256 | 0.011/0.025 | 5.92/5.88 |
| 10 / 1024 | 0.023/0.105 | 7.67/7.77 |
| 12 / 4096 | 0.115/0.520 | 15.61/17.34 |
| 14 / 16384 | 0.578/3.237 | 52.43/60.21 |

The cache and growing Krylov arrays explain most of the large-dimension RSS increase; tiny problems are dominated by runtime/library pages. NOM/workspace allocations can be partly untouched, so allocation-byte estimates and resident pages need not coincide. See all raw per-case timings, achieved errors and ranges in the evidence CSV and aggregate JSON.

## Memory estimates

Let n be sector dimension, m=NOM and k=NOK. The major complex vector/matrix
allocations (16 bytes per complex value) give these leading estimates:

| Solver/version | Leading vector and projected storage |
| --- | --- |
| ALG=1, Phase 1-B | about `64*n + 112*iterations` bytes |
| ALG=1, Phase 1-C | up to `192*n + 112*iterations` bytes during reconstruction |
| ALG=2, Phase 1-B | about `16*n*(m+k+1) + 40*m*m` bytes |
| ALG=2, Phase 1-C | about `16*n*(4*m+wanted+2) + 32*m*m` bytes |

These are allocation estimates, not measured RSS or rigorous maxima. LAPACK work,
compiler temporaries, allocator/page effects, modules, basis and shared libraries
add overhead. Some buffers have disjoint lifetimes. ALG=2's wrapper output storage,
Q/HQ/Ritz arrays and retained vectors explain its larger dimension-by-NOM term.
ALG=1 keeps linear n storage, despite its larger replay peak; it does not store
an n-by-iteration basis.

With automatic cache slots and 32-bit indices, cache storage is approximately
`20*(2*NO_one + 8*NO_two + 1)*cached_rows` bytes. This can dominate solver workspace.
The benchmark CSV records these leading estimates separately from measured RSS.
Changing cache can change time/RSS without changing the physical operator.
There is no general allocation-failure handling or huge-dimension safety claim.

## Preserved suites and compiler coverage

| Suite | GNU 14 default/O2/bounds (each) | GNU 13 O2 |
| --- | ---: | ---: |
| Clean build, no mismatch workaround | pass | pass |
| Original independent cases | 21/21 | 21/21 |
| Kernel / CLI stability | 82/82; 26/26 | 82/82; 26/26 |
| Sparse-field CLI | 70/70 | **69 passed / 1 failed** |
| Operator / permutation | 10/10; 28/28 | 10/10; 28/28 |
| LAPACK interface/error checks | 11/11 | 11/11 |
| Standard / extended examples | 1/1; 9/9 | 1/1; 9/9 |
| Historical comparisons | 27/27 | 27/27 |

The new optional Automake target was separately tested in a fresh O2 tree, with
418 recorded stress results, 11 LAPACK checks and the standard regression passing.
The target is characterization, not a claim that 418 accuracy checks all passed.

The GNU 13 sparse-field failure remains in the committed failure record and log;
no rerun replaces it. Its unseeded production start was not captured, so the
original run has input/binary/environment provenance but cannot be replayed
bitwise. This is a reproducibility limitation inherited from the CLI.

A separate seeded physical CLI probe (`tests/krylov_acceptance_probe.py`) runs
1,024 cases on that same single-site-four fixture with ALG=1/2, OMP=1/2,
tolerances 1e-14/1e-18 and 128 recorded seeds per combination. In this sample,
768 solves pass the unchanged observable assertions and 256 ALG=2/tol=1e-18 runs
report failure; no accepted observable failure was found in this controlled
sample. This does not erase the original failure or prove its absence elsewhere.
It illustrates the limitations of sampled random starts and excessively tight
energy gates. Twelve selected kernel runs replay exactly in RNG state, status,
reason, H-call count, eigenvalues and residuals; timing fields are excluded.

GNU Fortran 14.2 and 13.3 were tested against the same OpenBLAS 0.3.29 pthread
LP64 backend, Debian 13 x86_64, and NumPy 2.3.5 test environment. GNU 13 packages
were extracted locally from trusted Debian packages, without system changes.
This is two versions of one compiler family, not independent Intel/LLVM/MKL
coverage. No ifx/ifort/flang executable was initially available; MKL/ILP64 and
another frontend remain unvalidated.

## Reproduce and inspect evidence

Activate installed GNU Fortran/Autotools/BLAS and test-only NumPy. Cloud activation
is `. /workspace/.ed2-env/activate.sh`; elsewhere use the system equivalents.
Use fresh external trees and outputs; the tools refuse existing output directories.

```sh
set -eu
repo=/path/to/ED2
stage=/path/to/new-o2-tree
mkdir "$stage"
git -C "$repo" archive HEAD | tar -x -C "$stage"
cd "$stage"
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j2
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
make check-krylov-stress KRYLOV_OUTPUT=/path/to/new-stress-results
make check-lanczos LANCZOS_OUTPUT=/path/to/new-stability-results
make check-sparse-field SPARSE_FIELD_OUTPUT=/path/to/new-field-results
make check-lapack
make check
make check-long
python3 "$repo/tests/numerical_baseline.py" --exe "$stage/source/QS3ED2" \
  --output /path/to/new-independent-results
python3 "$repo/tests/compare_historical.py" --build "$stage" \
  --output /path/to/new-historical.json
```

Repeat in separate trees with the default flags (`./configure FC=gfortran`),
bounds (`-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra`), and GNU 13 O2.
Never reuse compiler module files across versions. No mismatch workaround or
blanket floating-point traps are used. Compile the stress driver through the
optional target or against that profile's matching interface/kernel objects.

After CPU-intensive tests finish, run benchmarks alone from the Git checkout:

```sh
python3 tests/krylov_benchmarks.py --repo "$repo" --output /path/to/new-benchmarks
python3 tests/krylov_acceptance_probe.py \
  --exe /path/to/new-benchmarks/phase1c/source/bench --output /path/to/new-acceptance
```

The benchmark needs both historical Git objects and a GNU-compatible OpenMP
compiler. It deliberately uses the two recorded source revisions, even after the
new test scripts are committed. The source checkout is never built or patched.
Raw input decks, vectors and logs remain in the external run tree. Committed
[evidence](phase1d-krylov/README.md) contains reduced per-case CSVs, provenance,
selected diagnostics and summary tables; no vectors, objects or package binaries.
`tools/krylov_characterization_summary.py` regenerates that reduction from the
profile-named raw layout described by the recorded command manifests.

## Limits and recommendations

There is no exhaustive spectrum, large-scale/mixed-spin performance, reentrancy,
ABI, distribution-packaging, invalid-model or another-backend validation here.
The new tests do not repair existing global ownership or CLI persistence issues.
Benchmarks are moderate identity-sector spin systems, complemented by preserved
historical symmetry/mixed-spin regressions. Their timing cannot establish general
HPC scaling or speedup.

Review the separate convergence contract proposal before extracting a public
solver API. In particular, distinguish energy changes, residual tolerances,
cluster/subspace accuracy and global inclusion, and expose budget exhaustion
without publishing uncertified results. Add explicit start/RNG ownership at the
future API boundary. Keep the present findings and failure data as regression
inputs for a separately authorized numerical-policy change.
