# Phase 1-G: opt-in additive residual convergence

Date: 2026-11-12 (Asia/Tokyo). Existing branch:
`codex/ed2-phase1g-optin-residual-policy`. Fetched/check-out HEAD was verified at
`00702f66d404f9b81fd7db241b7b5da30fb71e2d` before editing. Shell and attached runtime
were operational. Read `AGENTS.md`, Phase 0–1-F reports, especially the
[Phase 1-E contract](phase1e-precision-contract.md) and
[Phase 1-F implementation](phase1f-diagnostics-reproducibility.md).

## Numerical contract and scope

The new optional `residual_policy` is an **additional** necessary condition for
ALG=1/2 success. Absence preserves Phase 1-F behavior. No new nonzero tolerance
is chosen, no energy criterion is replaced, and no Hamiltonian, symmetry, LAPACK
interface, reference data or normal CLI record format changes. ALG=2 projection,
orthogonalization, block initialization, retention and restart algorithms are
unchanged. This is not a public-library API freeze or a numerical optimization.

For each normalized requested pair, the kernels use their existing fresh H action
and direct residual `r_i=||Hq_i-lambda_i q_i||_2`, including after ALG=1 reconstruction.
Ritz values are not replaced by Rayleigh quotients and projected-beta estimates
are not used instead of actual residuals. The additional requirement is

```
r_i <= atol + rtol * S_i
S_i = max(B, abs(lambda_i-c))       when rtol > 0
```

`atol` has energy units; `rtol` is dimensionless. `B` is an explicitly supplied
energy scale for `H-cI`, and `c` is the energy origin. All parameters must be
finite, tolerances/scales nonnegative, and at least one tolerance positive.
A positive `rtol` requires an explicit `B`, including an explicitly specified zero
when appropriate. The caller's scale is labeled **unverified**, not a certified
spectral norm. Its suitability must be justified by the caller if a bound-based
interpretation is needed. An observed maximum action norm is not a spectral bound;
there is deliberately no automatic substitution of Phase 1-F's observed scale.

For absolute-only requests (`rtol=0`), the scale is unused and recorded as zero;
the code does not even form an unnecessary `lambda-c` difference. Zero B and
`lambda=c` give S=0 without inserting an arbitrary unit or epsilon floor. The
origin defaults to zero; it changes only the **new tolerance scale**, not H,
its application, the recurrence, Ritz values or legacy energy checks. Centering
therefore does not recover precision already lost in forming/applying a shifted H.

The threshold is evaluated in binary64. Nonfinite parameters, missing required
scale, zero/zero tolerances and negative values fail validation before iterative
H actions. Overflowing supplied thresholds are rejected; an overflowing candidate
energy difference/threshold produces numerical failure, never an infinite
acceptance threshold. Subnormal representability can limit attainable requests;
no roundoff floor or automatic tolerance relaxation is introduced.

### Integration with legacy termination

| Exit path | Mandatory conditions |
| --- | --- |
| Ordinary ALG=1 | Existing energy criterion and existing residual certificate, plus requested residual |
| Ordinary ALG=2 | Existing two consecutive energy checks on retained values and existing certificates, plus requested residual for **every** wanted pair |
| Existing exact/near closure or full-space completion | Existing structural exit/certification rules, plus requested residual for **every** wanted pair |

The legacy envelope `max(256*epsilon,sqrt(lnc_ene_conv))*max(observed_scale,abs(lambda),tiny)`
and breakdown threshold remain untouched. A loose new request cannot weaken that
legacy certificate. The meaning of `lnc_ene_conv`, the near-zero energy extension,
iteration/cycle budgets and existing structural exceptions remain unchanged.
There is no new residual-only exit: existing closure/full-space exits already
had their own exceptions to ordinary energy-history checks.

ALG=1 continues its existing recurrence/reconstruction checks if a candidate fails
only the added condition. Uncertified breakdown still fails (status 2); budget
exhaustion still fails (status 1). ALG=2 retains its existing bounded cycle flow;
if requested residuals remain unmet, it fails rather than returning the last
candidate as successful. No additional block size, reorthogonalization, special
roundoff success or restart optimization is selected. Both CLI adapters retain
nonzero failure exits and suppress successful eigenvalue/wavefunction files.
Diagnostic last-candidate records are not successful result files.

## User controls and diagnostics

Existing namelists and normal stdout/file formats are preserved. New controls
are environment variables, read alongside the Phase 1-F seed/start controls:

| Variable | Meaning |
| --- | --- |
| `ED2_RESIDUAL_ATOL` | Explicit absolute residual tolerance |
| `ED2_RESIDUAL_RTOL` | Explicit relative residual tolerance |
| `ED2_RESIDUAL_SCALE` | B for `H-cI`, required when rtol is positive |
| `ED2_ENERGY_ORIGIN` | c, default zero |

At least one tolerance variable activates the policy; an omitted additive term
is zero, not an implicitly selected positive tolerance. Scale/origin alone are
rejected to avoid silently ignoring an incomplete request. With all four controls
unset/empty, the optional policy is absent and legacy behavior applies. Controls
apply only to actual ALG=1/2 solves, as with Phase 1-F; ALG=3/readback do not use them.

An absolute-only **example**, not a new recommendation/default:

```sh
OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
ED2_RNG_SEED=19 ED2_RESIDUAL_ATOL=1e-10 ED2_DIAGNOSTICS=solver.txt \
/path/to/source/QS3ED2 < input.dat > output.dat
```

For a caller that has independently established scale B=2 for `H-1e6 I`, a relative
example is `ED2_RESIDUAL_RTOL=1e-10 ED2_RESIDUAL_SCALE=2 ED2_ENERGY_ORIGIN=1e6`.
Both terms may be specified. An internal Fortran caller can append `policy=p`
to either kernel call using `type(residual_policy)` from `krylov_solvers`; old
calls remain source-compatible. A present zero-initialized policy object is
invalid; absence, not a zero/zero request, selects legacy mode.

`krylov_residual_policy.f90` owns parameter validation and the scalar threshold
calculation. The optional argument follows all Phase 1-F arguments. The existing
observations record the requested atol/rtol/B/c, explicit-scale flag/source,
per-pair S and threshold, per-pair requested/combined results, all-pair result and
failure reason. Existing `residual_certified`/pair-column meanings stay **legacy**;
new key/value fields carry the added criterion. With no policy the sidecar stays
in its prior layout. Observed legacy scale and supplied requested scale remain
separate fields. H counts still count every kernel callback and exclude the
existing ALG=1 CLI quality check after return.

A successful residual check neither proves lowest-state inclusion nor establishes
magnetization/correlation accuracy. Inclusion remains `not_assessed`; no gap or
observable certificate is fabricated. Requested residuals may be satisfied while
the unchanged energy gate still fails. Roundoff risk remains diagnostic and
`roundoff_limit_proven` remains false.

## Matched-start experiments and observed failures

`tests/krylov_policy_driver.F90` calls the actual kernels, counts all H actions and
exports vectors only on success (failed ALG=2 diagnostic vectors use distinct
filenames). The optional base build compiles that same driver against the exact
Phase 1-F sources without enabling the new argument. The new test reuses all
136 Phase 1-E definitions and adds full-space/all-pair, tiny-gap spin-dimer and
exact scalar controls, for **146 contexts** per compiler, each at OMP=1/2 where
applicable. It uses Phase 1-F's explicit per-call seed, preserves input matrix and
start bytes, and records/compares RNG state metadata. BLAS threads are one and
OpenMP dynamic teams are disabled.

Every context compares legacy mode with atol=1e-10, centered rtol=1e-10 and
atol=1e-300. Centered tests supply B from their independent reference spectrum and
c=1e6 for the large shifted family (c=-1 for the tiny-gap dimer). These are
**experiment controls**, not newly selected program defaults. A loose atol=1e30
additionally checks identical trajectories/output bytes. Independent calculations
assess actual residuals, Rayleigh values, spectral/subspace errors and physical
observables; failed ALG=1 candidate vectors remain unavailable for independent
assessment and are not assigned invented reference errors.

GNU 14 O2, GNU 13 O2 and GNU 14 bounds each produced exactly:

| Condition | Success (0) | Budget failure (1) | Breakdown/numerical failure (2) |
| --- | ---: | ---: | ---: |
| Legacy, no policy | 116 | 30 | 0 |
| Additional atol=1e-10 | 102 | 44 | 0 |
| Additional centered rtol=1e-10 | 107 | 37 | 2 |
| Additional atol=1e-300 | 11 | 130 | 5 |

All 146 legacy results per profile match Phase 1-F byte-for-byte in serialized
values/vectors, status/reason, legacy iteration count, H calls and RNG states.
Loose-policy results also match that trajectory. The additional absolute
condition turns 14 previous successes into failures; the centered relative
condition turns nine into failures; the ultra-tight condition turns 105 into
failures. No previously failed context becomes successful. Exact zero-residual
cases can legitimately satisfy even the ultra-tight positive request, so it is
not labeled universally unattainable. Failures stay in all tables and CSVs.

The supplementary GNU 14 default `-Ofast` profile records legacy 115/31/0,
absolute 102/44/0, centered relative 106/38/2 and ultra-tight 11/130/5. All 146
no-policy contexts match its own Phase 1-F `-Ofast` baseline; differences from O2
are existing flag-dependent numerical behavior, not default-policy changes.

### Representative cases (GNU 14 O2)

- **Clustered ground values:** `0016-t1`, ALG=1, has legacy residual 3.47e-7 and
  ground-projector error 0.904. Atol=1e-10 lowers the residual to 8.93e-11 and
  raises H calls from 36 to 197, but projector error remains 0.904. The allowance
  is still comparable to the internal 1e-10 gap: residual improvement is not
  individual-ground-vector resolution.
- **Excluded ground:** `0042-t1` returns the exact excited eigenpair with residual
  zero and ground-projector error one, even with atol=1e-300. The weak-overlap
  case `0043-t1` also remains the wrong ground vector at ordinary requested
  tolerances; atol=1e-300 rejects its numerical near-breakdown instead of inventing
  a missing-ground search. Inclusion remains explicitly unassessed.
- **Large shift:** the ALG=1 shifted cases accepted under legacy scaling can
  fail both new strict conditions. The supplied centered scale avoids letting
  1e6 dominate the *new* tolerance, but does not remove action/reconstruction
  cancellation or change the old breakdown rules. The failures are retained,
  not repaired by silently relaxing atol or applying a different H.
- **Every wanted pair:** `0068-t1` completes the full 4D space with requested
  residuals `[0,3.58e-15]`. Atol=1e-300 marks the first pair true and the second
  false and returns failure after its one-cycle budget. A first-pair-only test
  would be wrong. A separate boundary check sets atol to the maximum measured
  residual and succeeds, testing the inclusive `<=` relation.
- **Actual physical-accuracy counterexample:** `0069-t1` is a two-spin Heisenberg
  matrix with J=1e-15 and a -I shift. The analytic singlet independently fixes
  ground magnetization/correlations. ALG=1's residual is only 4.17e-16, yet the
  returned mixture has magnetization error 0.5 and correlation error 0.25. It
  passes both legacy and atol=1e-10; the centered relative request instead fails
  at the unchanged near-breakdown guard. This is an actual kernel outcome,
  distinct from Phase 1-E's retained constructed-state control. Neither success
  nor failure is a new claim of observable precision or global inclusion.

Invalid-policy controls cover zero/negative/nonfinite tolerances, missing or bad
scales, overflowing thresholds, and orphan CLI scale/origin settings. Both kernels
reject six malformed policy objects before any H action. A huge-origin guard
rejects an overflowing candidate difference; an absolute-only request with that
unused origin still succeeds on the exact scalar operator. Combined atol+rtol,
per-pair Ritz-dependent scales, positive rtol with explicit zero scale, exact zero
operators, strict energy requests and normal-output/result-suppression controls
are checked. Scalar zero completion still passes its requested condition without
a numerical floor. No invalid policy or failure is accepted via a legacy fallback.

## Regression results and arithmetic modes

The complete requested matrix ran in fresh before/after source trees. Each of
GNU 14 -O2, GNU 13 -O2 and GNU 14 bounds had these results:

| Suite | Result per build |
| --- | --- |
| LAPACK contracts / standard example | 11 / 1 passed |
| Independent dense physics | 21 passed |
| Lanczos stability / CLI failures | 82 / 26 passed |
| Magnetic fields / operators / permutations | 70 / 10 / 28 passed |
| Extended examples | 9 passed |
| Historical numerical comparisons | 27 passed |
| Phase 1-F telemetry before/after | 136 identical; 106 successful, 30 nonconverged |
| Phase 1-F seeded replay / error / CLI controls | 24 / 6 / 22 passed |
| Existing stress characterization | 164 strict accepts, 135 accuracy/inclusion warnings, 119 nonconvergences (418 total) |
| Existing precision characterization | 92 reference-target accepts, 14 warnings, 30 nonconvergences (136 total) |
| Existing seeded observable acceptance probe | 768 observable agreements, 256 solver failures (1,024 total) |

The last three rows deliberately retain previously documented deficiencies. They
are not claims that every solver run succeeded. In the acceptance probe the 256
failures occur with ALG=2's strict 1e-18 legacy energy request; none produces an
eigenvalue result file. No accepted observable disagreement was observed there.
All required command exit codes are zero. The standalone bounds
`check-residual-policy` and `check-diagnostics` targets were also executed successfully.

An additional GNU 14 default `-Ofast` **quick** run exercised both OMP counts,
policy cases, before/after comparisons, timing and the shorter existing suites.
Extended/historical/1,024-case acceptance were **not** rerun for this additional
profile. Its first Phase 1-F telemetry command failed: the Python checker assumed
gradual underflow while both the old and new executable used flush-to-zero.
For the zero operator, the unchanged breakdown expression becomes 0 instead of
the Python subnormal 3.16e-322. The unchanged legacy residual-limit expression has
the same issue. This was a test expectation problem, not a solver-policy change.

The final test-only driver now **queries**, after the timed solver call, the IEEE
underflow mode; it never sets it. The checker accounts for flushing only when the
reported runtime mode is nongradual and the otherwise expected positive result
is subnormal. All ordinary threshold assertions retain their previous tolerances.
A separate follow-up rebuilt only this driver, linked it against both original
and current solver modules, and reran telemetry/seed/error/CLI comparisons for all
four profiles. All 12 follow-up commands passed; default had 105 successful and
31 nonconverged cases, identical before/after. The three requested profiles report
gradual mode; default reports nongradual. Original command manifests retain the
initial failed default command; the correction has its own manifest and hashes.
Production sources remained identical to the full validation snapshots throughout.

## Accuracy, H-action count, runtime and memory

The reproducible timing workload is a dense complex Hermitian n=192 matrix with
QR-generated eigenvectors (NumPy seed 173) and eigenvalues in [-1,2]. ALG=2 requests
two pairs, capacity 32, retention 6 and budget 400; ALG=1 retains its original
two-pass reconstruction. Energy tolerance is 1e-10. Each compiler/ALG/thread group
runs five repetitions each of legacy, loose absolute 1e30, and strict absolute
1e-10, with the same initial state/seed and alternating order (60 runs/profile).
All 240 timed runs succeeded. Loose controls preserve eigenvalue/vector bytes,
termination and H counts. Times below are solver-only medians in milliseconds;
sidecar and vector export are outside the timer. Each cell is legacy / loose /
strict. No host-exclusive allocation or CPU pinning was used.

| Build | ALG | OMP | H applications | Median ms |
| --- | --- | --- | --- | --- |
| o2 | 1 | 1 | 122 / 122 / 2355 | 5.031 / 4.995 / 79.296 |
| o2 | 1 | 2 | 124 / 124 / 2294 | 4.965 / 5.141 / 77.101 |
| o2 | 2 | 1 | 1020 / 1020 / 1020 | 60.150 / 60.728 / 61.094 |
| o2 | 2 | 2 | 1020 / 1020 / 1020 | 59.282 / 56.538 / 59.148 |
| gnu13 | 1 | 1 | 122 / 122 / 2355 | 4.476 / 4.765 / 70.359 |
| gnu13 | 1 | 2 | 124 / 124 / 2294 | 6.776 / 7.056 / 74.541 |
| gnu13 | 2 | 1 | 1020 / 1020 / 1020 | 65.653 / 61.324 / 61.759 |
| gnu13 | 2 | 2 | 1020 / 1020 / 1020 | 62.375 / 65.210 / 66.157 |
| bounds | 1 | 1 | 122 / 122 / 2355 | 4.618 / 4.615 / 68.445 |
| bounds | 1 | 2 | 124 / 124 / 2294 | 4.665 / 5.410 / 74.678 |
| bounds | 2 | 1 | 1020 / 1020 / 1020 | 59.110 / 62.715 / 60.254 |
| bounds | 2 | 2 | 1020 / 1020 / 1020 | 62.195 / 65.564 / 60.293 |

For GNU 14 OMP=1 ALG=1, independent residual falls from 6.08e-6 to 9.64e-11;
the reference projector error falls from 2.25e-5 to 2.09e-10. This costs 122 to
2,355 H applications and 5.03 to 79.30 ms. OMP=2 changes the floating-point
trajectory (124 to 2,294 actions), while each same-thread replay is deterministic.
ALG=2 already yields approximately 1.2–1.4e-15 residual on this workload; every
policy uses 1,020 H actions in the requested profiles. The new gate adds no H
application to an existing candidate certificate. ALG=1's increased work is
additional existing recurrence/reconstruction attempts, not an ALG=2 optimization.
Loose-gate time differences vary in sign across groups; this small shared-host
sample does not establish a reliable percentage overhead or speedup.

Process VmHWM is sampled after the solve/sidecar and before exporting large vector
files. Across requested profiles and modes it spans 5,128–6,724 KiB; for GNU 14
ALG=1 OMP=1, legacy spans 5,132–5,244 and strict 5,168–5,292 KiB. These include
runtime/BLAS/diagnostic allocations and startup variation, not isolated allocation
measurements. Per-mode ranges and all raw observations are in `timing.csv` and
`timing-summary.json`. The policy adds scalar state plus diagnostic arrays of
O(wanted); it adds no new n-dimensional buffer or Krylov storage. Additional
iterations can still affect runtime/high-water memory. Metadata provenance records
compiler versions, flags, linked library hashes, CPU and the affinity/quota observed
when reducing results; it does not assert stable clock frequency or dedicated CPUs.

## Reproduction and retained evidence

The cloud toolchain setup from Phase 1-F is reused (LP64 OpenBLAS/LAPACK; NumPy
2.3.5). No new mandatory runtime dependency is introduced. The following creates
fresh external source trees and never builds over the working checkout:

```sh
. /workspace/ed2-env/activate.sh
python3 -B tools/phase1g_validation.py --output /workspace/phase1g-validation
python3 -B tools/phase1g_validation.py --profiles default --quick --output /workspace/phase1g-default-validation
python3 -B tools/phase1g_fp_validation.py --runs /workspace/phase1g-validation --default-runs /workspace/phase1g-default-validation --output /workspace/phase1g-fp-final
python3 -B tools/summarize_phase1g.py --runs /workspace/phase1g-validation --default-runs /workspace/phase1g-default-validation --fp-runs /workspace/phase1g-fp-final --output /workspace/phase1g-evidence
```

Output directories must not exist; choose new names for a replay. Paths to GNU 13
and library search directories are configurable in the main runner. The separate
FP follow-up documents the actual test correction; a clean replay already includes
the corrected tests, so the initial default failure is not expected to recur.
The runner records failures and continues independent checks, then exits nonzero
if any required command failed. `--quick` explicitly omits extended, historical and
acceptance suites; it is not equivalent to the complete requested validation.
Optional in-build targets are `make check-residual-policy` and
`make check-diagnostics` (the policy driver needs Python/NumPy).

The [evidence index](phase1g-residual-policy/README.md) describes committed CSV/JSON
results, selected independent assessments, command exit codes, hashes and retained
failure classifications. Full transient builds and raw vectors remain under the
external output directories; generated binaries are not committed. The source
snapshots predate the two-file FP test correction and final report/evidence tools;
production-source hashes are identical. Follow-up test hashes identify the final
checker exactly. Protected Hamiltonian/state/input/LAPACK/main/observable files,
legacy threshold/random/orthogonalization helpers and ordinary CLI write/namelist
statements were checked against the starting commit without changes.

## Remaining decisions and independent follow-ups

- The internal trailing policy type and environment variables are provisional;
  choosing a stable library API or namelist schema belongs to a later phase.
- Caller-supplied scales remain unverified. No general safe operator norm bound is
  inferred from observed action norms. Domain-specific certified scale providers
  need independent derivations and tests.
- Residual-only termination, explicit roundoff-limit acceptance, automatic cluster
  expansion, better near-breakdown handling and shift-centered H application are
  separate changes. None is enabled here. Unattainable requests explicitly fail.
- A residual certifies proximity to some eigenpair, not global ground-state
  inclusion. The exact excited-state and weak-overlap controls remain limitations;
  additional starts or independent inclusion evidence would be required.
- Observable error requires gap/subspace information and operator bounds. Small
  residuals alone do not fix magnetic/correlation errors at tiny gaps. GNU 13's
  prior magnetization discrepancy is not repaired by this optional gate; the
  unchanged legacy probe retains 256 explicit solver failures under strict energy
  tolerance. No reference values were adjusted.
- Cross-compiler/thread bitwise reproducibility is not promised. Seeds control
  the current compiler's RNG sequence, and reduction order still depends on OMP;
  comparisons here fix both seed/state and thread count.
- Broader performance claims require larger physical systems, controlled placement
  and repeated allocations. ALG=2 speed or algorithm redesign is outside this work.

## Changed-file map

- `source/krylov_residual_policy.f90`: internal parameter type, validation and
  overflow-aware tolerance calculation.
- `source/krylov_solvers.f90`: optional argument and additive success gates;
  `source/krylov_observations.f90`: CLI controls and optional diagnostic fields.
- `source/lanczos.f90`, `source/eigen_solver.f90`: CLI adapter forwarding;
  root/source/tests `Makefile.am`: module order, linkage and optional test target.
- `tests/krylov_policy_driver.F90`, `tests/residual_policy_validation.py`: matched
  baseline/policy/adversarial and CLI checks; `tests/krylov_benchmarks.py`: module
  order for the existing external acceptance probe.
- `tests/krylov_diagnostics_driver.F90`, `tests/diagnostics_reproducibility.py`:
  test-only query/check of actual gradual-underflow mode.
- `tools/phase1g_validation.py`, `tools/phase1g_policy_timing.py`,
  `tools/phase1g_fp_validation.py`, `tools/summarize_phase1g.py`: reproducible
  orchestration, timing, test follow-up and evidence reduction.
- This report and `phase1g-residual-policy/`: design, comparisons, failures,
  performance, provenance and remaining decisions. No reference data changes.
