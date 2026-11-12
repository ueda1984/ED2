# Phase 1-F: internal diagnostics and reproducible initialization

Date: 2026-11-12 (Asia/Tokyo). Branch:
`codex/ed2-phase1f-diagnostics-reproducibility`. The existing remote branch was
fetched, checked out, and verified at the required starting HEAD
`b756b0009f08027e8a2ae8a1e349d78d8e81d9c5` before edits. Shell execution worked;
the attached runtime reported running/connected. Read `AGENTS.md`, the Phase 0
through 1-E reports, CPC plan, Phase 1-D proposals, and the
[Phase 1-E precision contract](phase1e-precision-contract.md).

This phase implements that contract's **telemetry and initialization groundwork**.
It does not implement the proposed residual policy, select new defaults, optimize
ALG=2, or make a stable public-library API. Hamiltonian construction, symmetry,
LAPACK interfaces, reference datasets, `lnc_ene_conv`, and every numerical stopping
threshold/formula remain unchanged. No merge, main edit, or PR is authorized.

## Internal implementation and ownership

`source/krylov_observations.f90` defines the per-call `krylov_diagnostics` result,
its initialization/finalization, optional CLI controls, and a separate text writer.
The kernels append optional `diagnostics` and `rng_seed` arguments after the
existing optional `initial` argument; every old call remains valid. Thin wrappers
finalize observations on every normal return, including invalid input, numerical
failure and budget exhaustion. The recurrence/restart implementations remain in
private cores. Diagnostics own only scalar metadata, requested-pair arrays and,
for an explicitly seeded call, a small intrinsic RNG state table. They do not own
or retain an n-by-iteration basis or another copy of the returned vectors.

Example (the internal test seam, not a stable API promise):

```fortran
use krylov_solvers, only: lanczos_ground, krylov_diagnostics
type(krylov_diagnostics) :: observation
! Other arguments retain their Phase 1-C definitions.
call lanczos_ground(apply,n,tol,maxiter,miniter,interval,value,vector, &
                    iterations,status,reason,residual, &
                    diagnostics=observation,rng_seed=19)
```

An existing `initial=` vector is still normalized and processed by the same
kernel. For ALG=2 it remains an **additional direction** alongside the existing
generic block, not a replacement of that block or a new block-start algorithm.
Use a seed as well when reproducibility of those additional random directions
is required. An ALG=1 explicit start needs no random draw.

### Recorded meanings and availability

| Field group | Meaning |
| --- | --- |
| `status`, `reason`, `termination` | Original integer status and text plus an additive classification; never overwrite a failure with success |
| `iterations`, `cycles`, budgets | ALG=1 recurrence steps and effective `min(n,maxiter)` budget; ALG=2 H applications and Rayleigh-Ritz-cycle budget; cycle count separate |
| `h_applications` | Every actual kernel callback, including reconstruction, certification, initialization, expansion and retained-vector refresh |
| `ritz`, `rayleigh`, `residual`, `norm` | Last directly assessed requested pairs; Rayleigh is `real(q^H Hq)/||q||^2`, residual uses the original Ritz value |
| `pair_available`, `candidate_role` | Distinguish available last-candidate measurements from returned successful outputs; a failed solve may have useful candidate diagnostics |
| `residual_certified`, `residual_threshold` | Whether the measured residual meets the **unchanged** current residual envelope; not an observable-accuracy guarantee |
| energy change, gate flag, consecutive checks | Actual legacy energy-change telemetry; ALG=2 still needs two consecutive ordinary energy checks |
| norm/orthogonality | Actual norms and maximum absolute entry of `Q^H Q-I`, explicitly **not** its spectral norm |
| scale/breakdown data | Actual scale used by current stopping logic, availability, threshold, last beta, closure and exact/near-closure event counts |
| inclusion | Always `not_assessed; no independent lowest-state evidence` in this phase |
| roundoff | A below-binary64-epsilon energy-request risk flag and an explicitly unproven assessment; no new successful termination class |
| initialization | Explicit-start/seed flags, seed, random-vector count, configured maximum threads, active seeded threads, compiler/runtime RNG identity and first seeded states |

ALG=1's scale is the maximum observed action norm/absolute recurrence diagonal
in its continuing first pass. ALG=2 updates the action-norm maximum during block
initialization and expansion, as before; retained-vector refresh/certification
calls do not update it. Neither quantity is a certified spectral norm. The
reported threshold is exactly
`max(256*epsilon,sqrt(tol))*max(scale,abs(Ritz),tiny)`, and breakdown remains
`64*epsilon*max(scale,tiny)`. No shift correction or residual absolute/relative
option was added. The H-action count is scoped to the internal solve: ALG=1's existing CLI quality
check after the kernel returns is outside it. No new action is performed to populate diagnostics: Rayleigh
and norm use the vector/action already available at residual certification.
ALG=2's requested-block Gram calculation runs once on return when all its pair
measurements are available. ALG=1's scalar Gram defect follows from its norm.

Availability flags govern values: unavailable scalar placeholders or unmeasured
array entries are not zero-error certificates. ALG=1 can report the last measured
candidate on failure without returning its vector as a successful result. ALG=2
observations reset pair availability at each new cycle; incomplete/error cycles
must not be treated as a complete successful block. The original `status` is the
authoritative success/failure result. A process-level LAPACK error still follows
its pre-existing fatal path; this phase does not convert it to a recoverable API.

### Nonconvergence, breakdown, and roundoff

`budget_exhausted`, exact closure, numerical near-closure and uncertified numerical
breakdown are reported separately from residual certification and inclusion.
Ordinary numerical failures also retain their original reason (for example,
nonfinite action). A below-epsilon request is only an advisory indication of
roundoff risk. `roundoff_limit_proven` remains false because the solver has no
validated action-error enclosure or gap information. A budget-exhausted candidate
with a tiny residual therefore remains nonconverged. This implements the honest
separation in the Phase 1-E specification without silently introducing a new
roundoff stopping policy or changing the requested tolerance.

## Opt-in CLI controls

No namelist field, normal stdout record, eigenvalue/wavefunction/observable file
format was added or changed. New controls are environment variables read only
when ALG=1/2 actually starts a solve:

| Variable | Use |
| --- | --- |
| `ED2_RNG_SEED` | A signed default-integer decimal seed; GNU LP64 builds accept -2147483648 through 2147483647 |
| `ED2_INITIAL_VECTOR` | Path to exactly `THS` text rows, each `real imag`, in the executable's **sector-basis order**, not an assumed product basis or native `wf*.bin` format |
| `ED2_DIAGNOSTICS` | A new sidecar path; existing files are rejected rather than overwritten; choose a path distinct from normal result files |

Example, in an existing input working directory:

```sh
OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
ED2_RNG_SEED=19 ED2_DIAGNOSTICS=solver-diagnostics.txt \
/path/to/source/QS3ED2 < input.dat > output.dat
```

Add `ED2_INITIAL_VECTOR=/absolute/path/start.txt` for an explicit vector. Zero or
nonfinite starts fail through the existing kernel validation; malformed/short/
extra files or invalid seed text fail before solving. Empty variables are treated
as absent. ALG=3 and `re_wf=1` readback do not invoke these controls. The sidecar is
written before the existing adapter checks the solver status, so failed solves
can retain diagnostics while still exiting nonzero and writing no successful
result files. A requested sidecar I/O failure is itself an error; no promise of
extra output is made on fatal runtime/LAPACK errors.

The sidecar is an internal `schema=1` key/value format, with explicit availability
and indexed pair rows. Scientific values use a three-digit exponent field so
binary64 subnormal exponents remain parseable. Pair rows are described by their
`pair_columns` record. It is separate from the historical normal CLI format and
is not a stable library or interchange-format commitment.

## RNG and OpenMP reproducibility

With no seed, the original `random_number` parallel loop and the CLI's existing
initialization remain unchanged. Merely requesting diagnostics does not inspect
or initialize unseeded intrinsic RNG states; doing so could perturb lazy thread
initialization. Such ambient unseeded states are therefore not falsely reported
as captured/replayable.

With an explicit seed, each random-vector draw temporarily initializes each
participating thread's intrinsic generator using

```
state(k) = modulo(seed + 104729*thread_id + 37*k + 1000003*draw_index,
                  huge(default_integer)-1)
```

The arithmetic is performed in int64 before conversion, avoiding overflow for
extreme permitted signed seeds. Draw index starts at zero for the first generic
vector and increments for later generic vectors. Static scheduling specifies
which thread draws each vector entry. Each participating thread's prior intrinsic
state is saved and restored. This supplies per-call independence from its prior
random draws without imposing a new ED2-wide global RNG state. It does not claim
reentrancy of the existing mutable CLI/model modules or recovery of uninitialized
future worker streams in every compiler runtime.

The generator is the compiler/runtime's `RANDOM_NUMBER`/`RANDOM_SEED` implementation,
not an ED2-selected portable RNG algorithm. Diagnostics record `compiler_version`,
thread information and the actual first seeded state vector for each active
thread. The formula, draw count, dimensions and compiler/runtime provenance
specify subsequent seeded draws. Fixed-seed replay requires fixed executable,
compiler/runtime/BLAS, input, thread count and execution conditions. Different
OpenMP thread counts partition streams differently; different compilers/runtimes
need not implement the same generator. Cross-thread-count or cross-compiler
bitwise equivalence is not promised. Disable dynamic teams and nested BLAS
parallelism for the documented reproduction runs.

## Verification and reproduction

`tools/phase1f_validation.py` archives the exact starting revision and snapshots
current maintained files into fresh external trees. It builds separate before/
after trees per compiler profile, records every command/status/log/hash, and
continues independent checks after a failed check. The checkout is never built
or converted into a worktree. Use a new output directory:

```sh
. /workspace/ed2-env/activate.sh
cd /workspace/ED2
python3 -B tools/phase1f_validation.py --output /workspace/new-phase1f-validation
```

Outside this Cloud instance supply the installed GNU 13 compiler via `--gnu13`
and appropriate LP64 search paths via `--ldflags`. GNU 14/13 O2 use
`-O2 -fopenmp`; bounds use `-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra`.
No mismatch workaround, blanket floating-point traps, verification bypass or
new mandatory dependency is introduced. Test-only NumPy is pinned as before.
The optional `make check-diagnostics DIAGNOSTICS_OUTPUT=/new/output` runs the new
checks against a current build; the full runner also builds its matching
before-version driver and records timing/RSS comparisons.

The replay harness reuses the Phase 1-E deterministic definitions and constructed
physical-accuracy counterexample. It compares status/reason, legacy iteration
count, H calls, stdout and serialized eigenpairs/vectors byte-for-byte between
before/after and diagnostics disabled/enabled. It independently recomputes
Rayleigh quotients, norms, actual residuals, Gram defects and the reported
thresholds for available vectors. Failed ALG=1 vectors remain unavailable for
independent recomputation; their direct in-kernel measurements are labeled as
such, not upgraded to independent evidence.

It also checks per-call seeded replay at OMP=1/2, exact restoration of participating
thread RNG states, the recorded seed mapping, explicit starts, extreme signed
seeds, malformed/nonfinite/zero starts, invalid parameters, budget exhaustion,
and CLI result suppression on error. Normal CLI file bytes and stdout (excluding
the already-existing elapsed-time value) match with the sidecar off/on.

### Final results

All three profiles independently reproduced the same observations:

| Check | GNU 14 O2 | GNU 13 O2 | GNU 14 bounds |
| --- | --- | --- | --- |
| Before/off/on kernel comparisons | 136 identical | 136 identical | 136 identical |
| Additional invalid/nonfinite kernel paths | 6 verified | 6 verified | 6 verified |
| Seeded replay / ambient state restoration | 24 pairs | 24 pairs | 24 pairs |
| CLI control/output/error checks | 22 passed | 22 passed | 22 passed |
| Original independent physics suite | 21/21 | 21/21 | 21/21 |
| Sparse-field CLI / operator / permutations | 70/70; 10/10; 28/28 | 70/70; 10/10; 28/28 | 70/70; 10/10; 28/28 |
| Stability kernel / CLI (final format) | 82/82; 26/26 | 82/82; 26/26 | 82/82; 26/26 |
| LAPACK interface/error paths | 11/11 | 11/11 | 11/11 |
| Standard / extended examples | 1/1; 9/9 | 1/1; 9/9 | 1/1; 9/9 |
| Historical full-array comparisons | 27/27 | 27/27 | 27/27 |
| Existing stress characterization | 418 outcomes | 418 outcomes | 418 outcomes |
| Phase 1-E characterization | 136 outcomes | 136 outcomes | 136 outcomes |
| Existing seeded physical CLI probe | 1,024 outcomes | 1,024 outcomes | 1,024 outcomes |

Each new 136-case comparison includes **106 successful solves and 30
nonconvergences**, identically before/off/on; a passing comparison is not a
claim that the nonconverged solves became successful. Each existing stress
characterization records 164 strict assessments, 135 accepted accuracy/inclusion
warnings and 119 nonconvergences. Each existing Phase 1-E run records 92 accepted
reference targets, 14 accuracy/inclusion warnings and 30 nonconvergences. These
numerical findings remain visible, even though the characterization commands
complete normally.

Each 1,024-case physical CLI probe records 768 observable agreements and 256
ALG=2/tolerance=1e-18 budget failures. Every such failure exits nonzero, reports
maximum restart iterations and leaves no eigenvalues result file. No accepted
physical accuracy failure was found in these finite samples. The Phase 1-E
constructed state still demonstrates near-exact energy with magnetization error
about 3.14e-8; diagnostics do not create an accuracy guarantee or erase the
original unseeded GNU 13 finding.

The exact excited-start case still has zero residual and a wrong ground-state
projector; its new inclusion field explicitly says `not_assessed`. Near-degenerate
mixtures and complete degenerate spaces retain their previous behavior. The
below-roundoff request cases still report budget exhaustion rather than a new
successful roundoff exit. Zero-Hamiltonian cases retain their existing guarded
breakdown/zero-energy behavior.

The stability/CLI counts include their existing designed failure controls. The
six new kernel error controls and CLI error checks pass only if the specified
failure is actually observed. Historical comparisons keep the existing 1e-9
energy and 1e-6 observable absolute/relative thresholds; independent physical
assertions remain 1e-9. No reference result or assertion was loosened.

### Diagnosed test serialization failure and correction

The first O2 and GNU 13 stability commands recorded 81 passes/one failure each:
not a solver-status error, but Python failing to parse a legitimate near-zero
finite value written as `-3.420637783583041-308` (O2) or
`-8.914189635110894-309` (GNU 13). Fortran's old test-only `ES25.15` format drops
`E` for three-digit exponents. The first commands therefore did not reach their
CLI suite; those CLI checks were initially **unrun**, not passed.

Only the test driver's serialization was changed to `ES26.15E3`. Numerical
precision, inputs, assertions, production output and algorithms were unchanged.
Fresh corrected runs on all three profiles executed and passed all 82 kernel
and 26 CLI checks. Initial failures, their process/solver statuses, and corrected
report hashes remain in the evidence. The raw command manifest intentionally
retains its two nonzero initial commands; the correction is recorded separately,
not disguised by overwriting statuses or replacing a reference answer. The
committed reproduction runner snapshots the corrected test driver on a fresh run.

## Measured time and memory impact

Measurements are serial, paired before/off/on kernel processes with OMP=1,
BLAS=1, a fixed matrix and initial RNG state; mode order alternates over seven
repetitions. No exclusive host, CPU-frequency or affinity control was applied.
The n=192 complex QR-rotated separated spectrum is fixed by NumPy seed 173.
Solver wall time is measured inside the Fortran driver and includes the optional
kernel diagnostics; text sidecar I/O and vector exports occur after that timer.
Process peak RSS is the executable's own `/proc/self/status` VmHWM, not the Python
parent's inherited `wait4` high-water mark.

For **converged workloads**, median solver milliseconds were:

| Profile / algorithm | Before | Diagnostics off | Diagnostics on | On/off | H calls, all modes |
| --- | ---: | ---: | ---: | ---: | ---: |
| GNU 14 / ALG=1 | 4.383 | 4.474 | 4.595 | 1.027 | 122 |
| GNU 13 / ALG=1 | 4.469 | 4.497 | 4.552 | 1.012 | 122 |
| Bounds / ALG=1 | 4.502 | 4.498 | 4.707 | 1.046 | 122 |
| GNU 14 / ALG=2 | 58.087 | 59.168 | 59.605 | 1.007 | 1,020 |
| GNU 13 / ALG=2 | 57.755 | 59.092 | 58.734 | 0.994 | 1,020 |
| Bounds / ALG=2 | 58.620 | 58.394 | 57.282 | 0.981 | 1,020 |

ALG=1 uses its unchanged recurrence at tolerance 1e-10. The converged ALG=2
workload explicitly requests wanted=2, keep=6, capacity=32, budget=400 and the
same tolerance. All 63 supplemental ALG=2 processes succeeded, matched the
independent lowest two energies at 1e-9, and reproduced before/off/on file bytes.
This is a test-workspace choice, not a new CLI default or ALG=2 optimization.

The originally requested smaller ALG=2 workspace (keep=4, capacity=12,
budget=120) exhausted its budget in **all 63** original timing processes. Those
measurements are retained separately, with identical 1,680 H calls and matching
failed outcomes across modes; they are not used as converged-workload timings.
ALG=1 contributed another 63 original timing processes. Full raw ranges,
statuses and peak RSS are in the [evidence](phase1f-diagnostics/README.md).

On converged workloads, ALG=1 off/on peak RSS ranges were respectively
5116–5256/5120–5260 KiB (GNU 14), 5116–5252/5120–5188 (GNU 13), and
5140–5212/5140–5220 (bounds). ALG=2 ranges were 6296–6504/6288–6500,
6292–6508/6288–6588, and 6312–6568/6304–6484 KiB. The maximum observed on-versus-off
peak increase is 80 KiB; allocator/page/loader variation prevents interpreting
that as a general memory-overhead bound. Scalar/pair telemetry is O(wanted), and
seeded-state storage is O(intrinsic-state-size times configured threads); no new
n-sized vector workspace or H action was introduced.

The small positive and negative timing differences do not establish speedups or
precise overhead percentages. The Gram assessment costs O(n*wanted**2) once on
return, and norm/Rayleigh observations cost O(n*wanted) at each certificate, so
large requested blocks can have more overhead than these wanted=1/2 tests.
Sidecar disk cost, very large problems and scaling remain unmeasured. Disabling
diagnostics avoids those additional numerical assessments; the default wrapper
and presence checks have the bounded before/off observations shown above.

The supplementary successful timing experiment can be reproduced after the full
runner finishes:

```sh
python3 -B tools/phase1f_success_timing.py \
  --runs /workspace/new-phase1f-validation --output /workspace/new-phase1f-success-timing
```

Its raw failed/successful workload distinction is intentional; do not discard
nonconvergence records or compare differing Krylov work as an algorithm speedup.


## Remaining design decisions

- Public result/context types, recoverable LAPACK/action/allocation failures and
  stable serialization remain later library work. This telemetry is internal.
- Certified roundoff/error enclosures, spectral gaps/cluster interpretation,
  observable error guarantees, lowest-state inclusion evidence and new absolute/
  relative residual controls remain the Phase 1-E proposal; none is inferred
  solely from the new fields.
- Portable owned RNG algorithms, complete explicit block starts and deterministic
  behavior across thread counts/compiler families require separate choices.
- Tests cover the recorded GNU versions and LP64 OpenBLAS backend, not an
  independent compiler family, MKL/ILP64, extreme dimensions or concurrent CLI
  instances. Cloud timings are bounded observations, not HPC scaling claims.
