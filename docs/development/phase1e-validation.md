# Phase 1-E: precision experiments and validation report

Date: 2026-11-12 (Asia/Tokyo). Branch:
`codex/ed2-phase1e-precision-contract`. Starting HEAD was fetched from that
existing branch, checked out, and verified as
`b56fe21678e62b8b6216fe0ec742d5dbd62a8df6` before edits. Shell execution worked;
runtime observations reported running/connected with enforced network policy.

Read `AGENTS.md`, the Phase 0 audit, CPC plan, Phase 1-A/B/C/D reports and
`phase1d-policy-proposals.md`. The resulting [precision contract](phase1e-precision-contract.md)
is a proposal for subsequent numerical work, not an implemented stopping rule.
No production Fortran, Hamiltonian, symmetry, numerical default, existing test
assertion or historical reference was changed. Only new validation/documentation
and an optional test target are added. No merge or PR is part of this phase.

## Findings and their limits

1. **Energy stagnation, pair validity and lowest-space accuracy are distinct.**
   For the n=32 near-cluster test, GNU 14 O2 case `0017-omp1` exits successfully
   under current energy/residual checks, with energy error `2.19e-10`, actual
   residual `1.70e-7`, and ground-projector error `0.981`. Leakage from the larger
   four-state cluster is only `1.38e-7`. This is a well-defined cluster membership
   result, but not an accurate individual ground vector. In the corresponding
   ALG=2 four-vector case `0020-omp1`, the complete-cluster projector error is
   `3.10e-15`, with residuals below `1.00e-15`. Those answer different questions.
2. **Residual zero does not prove inclusion.** The existing analytic excluded-start
   case returns the excited energy of `diag(-1,-1+1e-10,0,2)` with residual zero
   and ground-projector error one. The weak-overlap start `[1e-8,1,0,0]` also
   returns projector error one after numerical near-breakdown, with residual
   about `1.11e-16`. Both reproduce at OMP=1/2 in all four build profiles.
3. **Energy origin matters for relative controls.** After adding `1e6 I` to a
   separated n=32 Hamiltonian, case `0040-omp2` still reports ordinary success,
   but its residual is `4.19e-4`, ground-energy error `2.26e-7`, and projector
   error `6.67e-4`. The current energy and residual scales include that offset.
   Centered residuals are an offline alternative here, not a new production
   implementation. Reference eigensystem formation defects are recorded, and
   their action on the assessed block is included in the subspace-bound estimate;
   the shifted reference is not claimed exact below binary64 formation error.
4. **Below-roundoff energy requests can reject excellent candidates.** O2 case
   `0006-omp1` (ALG=2, tol=1e-20, 120 cycles) reports nonconvergence despite
   last-candidate residuals at most `2.71e-15` and energy errors at most `4.44e-16`.
   It remains a failed solve under the existing policy. The proposed roundoff
   diagnostics must not silently relabel it as satisfying `1e-20`.
5. **Physical accuracy needs its own assessment.** A deterministic, independently
   constructed state on the very same single-site-four physical Hamiltonian as
   the GNU 13 failure has energy error `2.22e-16`, residual `1.04e-8`, and maximum
   magnetization error `3.14e-8`. It lies within the current residual envelope at
   energy tolerance `1e-14`, yet fails unchanged `atol=rtol=1e-9` magnetization
   checks. This is a **constructed-state counterexample**, not a claim that a
   production solve returned this particular state.

The six analytic controls also check a stationary mixed-cluster state, its large
shift, exact excited-state residual, unitary rotations of a full degenerate
projector and the absence of a unique partial degenerate projector. Assertions
check the diagnostic distinctions; no failing solver result is hidden as an xfail.

## Experiment design and results

`tests/precision_contract.py` uses the unchanged `krylov_stress_driver.f90` to call
the exact production kernels. Every solve retains a full input deck, kernel status
and reason, process exit, iteration/H-call count, all thread RNG states and output
vectors in a fresh external directory. The test does not implement an alternative
Lanczos recurrence. Independent Python/NumPy calculations evaluate residuals,
orthogonality, projector/membership errors, reference spectral separation and
physical observables. ALG=2 last candidates on failure are labeled separately;
ALG=1 failed vectors are unavailable and are not assessed as if returned.

Per profile, 136 cases cover:

- Complex QR-rotated n=32 spectra with ground-cluster spacings 1e-2, 1e-6,
  1e-10 and exact degeneracy; ALG=1/2; two seeds; tolerances 1e-12/1e-20.
- Zero/scalar operators, factors 1e-9/1e9 and a 1e6 identity shift.
- Exact excluded and weak-overlap starts.
- The independently built n=16 single-site-four physical Hamiltonian, random
  starts at tolerances 1e-14/1e-18, plus explicit starts with ground amplitudes
  1e-2/1e-4/1e-6/1e-8 and sampled orthogonal complements.
- OMP=1/2 for every configuration, with BLAS=1 and dynamic teams disabled.

| Profile | Accepted reference targets | Accepted accuracy/inclusion warning | Nonconvergence | Total |
| --- | ---: | ---: | ---: | ---: |
| GNU 14.2 O2 | 92 | 14 | 30 | 136 |
| GNU 14.2 bounds | 92 | 14 | 30 | 136 |
| GNU 13.3 O2 | 92 | 14 | 30 | 136 |
| GNU 14.2 default Ofast | 92 | 13 | 31 | 136 |

All four characterization commands complete with exit zero and no harness/runtime
error, numerical-failure or invalid-parameter classification. This is **544
recorded solver outcomes**, not 544 successful accuracy checks. The accepted
reference-target label uses lowest-eigenvalue error `<=1e-10*max(||H||,1)` and
membership leakage `<1e-5`; it is not the proposed observable contract and does
not promise unique vectors at a degenerate cut. Thresholds are comparison targets,
not new defaults. The full 136-case O2 run replayed exactly in all recorded fields
(timings are not collected); the final optional Automake target also executed
136 cases in the bounds tree.

For GNU 14 O2, comparing the **same endpoints** gives:

| Offline condition | All available candidates (including failed ALG=2) | Candidates from successful solves |
| --- | ---: | ---: |
| Absolute residual <=1e-10 | 108 | 86 |
| Residual <=1e-10 times full spectral scale | 112 | 90 |
| Residual <=1e-10 times half spectral width | 108 | 86 |
| Complete-space reference bound <=1e-5 | 98 | 76 |

There are 128 available endpoint candidates in this run: 106 successful and 22
failed ALG=2 candidates. Eight failed ALG=1 runs do not expose their candidate
vectors and are excluded from endpoint comparisons. The criteria are neither
equivalent nor interchangeable. They do not retroactively turn the 30 failures into current-policy successes,
establish global inclusion, or predict earlier stop times, cost or speedup under
a new online rule. `last_printed_energy_change` is limited-precision telemetry
for ALG=1; it is not a recomputation of ALG=2's retained-value energy gate.

## GNU 13 magnetization failure: preserved evidence and new probes

The [original failure record](phase1d-krylov/gnu13-regression-failures.json),
[input](phase1d-krylov/gnu13-sparse-failure-input/input.dat), and
[stdout](phase1d-krylov/gnu13-sparse-failure.log) remain unchanged. That run accepted
ALG=1 with residual `2.644543803251125e-8` and magnetization error
`1.53227552e-8`. Its unseeded CLI RNG state was not retained, so its exact initial
vector cannot be recovered; reusing the input is not bitwise replay of that run.

The independent reference ground gap is `0.20806401984119738`. The residual over
this gap is about `1.27e-7`; thus a conservative spin-1/2 magnetization bound is
much larger than 1e-9. The observation is consistent with the permissive residual
threshold; it does not establish a GNU 13 compiler defect. The constructed-state
control above proves the missing implication on this physical fixture without
depending on a lucky stochastic failure.

For each of GNU 14 O2 and GNU 13 O2, reran the **existing unchanged**
`krylov_acceptance_probe.py`: ALG=1/2, tolerances 1e-14/1e-18, OMP=1/2, 128 recorded
seeds per combination, 1,024 actual CLI executions. Only external archived copies
are instrumented using the existing Phase 1-D probe: initialization of every
thread's RNG, H-call counting and basis export. No numerical kernel, stopping
condition or arithmetic expression is patched. Exact patches, compiler/binary
hashes and per-case inputs/RNG states are retained in the evidence.

Each compiler records 768 observable agreements, **zero accepted observable
failures**, and 256 ALG=2/tol=1e-18 failures. Every failure log identifies maximum
restart iterations, exits nonzero, and leaves no eigenvalues result file. These
failures remain failures, not expected successes. The current kernel physical
start experiments and the new unseeded sparse-field regression runs also found
no accepted physical accuracy failure. These finite samples do not erase the
original finding or prove accuracy for all starts. Both compilers share the same
LP64 OpenBLAS/runtime ecosystem; another compiler family/backend remains untested.

## Preserved regression checks and environment

| Check | GNU 14 O2 | GNU 13 O2 | GNU 14 bounds | GNU 14 Ofast |
| --- | --- | --- | --- | --- |
| Clean build without mismatch workaround | pass | pass | pass | pass |
| Kernel / CLI stability | 82/82; 26/26 | 82/82; 26/26 | 82/82; 26/26 | unrun |
| Original independent suite | 21/21 | 21/21 | unrun | unrun |
| Sparse-field CLI / operator / permutation | 70/70; 10/10; 28/28 | 70/70; 10/10; 28/28 | unrun | unrun |
| LAPACK interface/error paths | 11/11 | 11/11 | 11/11 | unrun |
| Standard Automake regression | 1/1 | 1/1 | 1/1 | 1/1 |

The stability counts include their designed invalid-input/nonconvergence controls;
passing means the specified status/behavior was observed. Extended historical
examples, historical-array comparisons and full Phase 1-D benchmarks were not
rerun: this phase changes no production source and makes no new claim for those
suites or for performance. No newly required check is left unexplained.

Debian 13 x86_64; GNU Fortran 14.2.0 and 13.3.0; Autoconf 2.72, Automake 1.17;
OpenBLAS 0.3.29 pthread/LAPACK LP64; NumPy 2.3.5 (test-only). GNU 13 packages were
downloaded from signed Debian metadata and extracted outside the checkout. Both
compiler profiles use the available libgfortran5 runtime. No TLS, signature or
checksum verification was disabled. Detailed source, binary, harness and linked
library hashes and exact flags are in [provenance](phase1e-precision/provenance.json).
All production Fortran bytes were checked against the starting Git revision and
against each uninstrumented build. No worktree was created and main was not edited.

## Reproduce

Use a **fresh** ordinary external staging/output directory per profile. The Cloud
activation path in this task is `/workspace/ed2-env/activate.sh`; other installations
can supply equivalent tools normally. Do not use the older onboarding `install.sh`
for this branch: its historical mismatch workaround is unnecessary after Phase 1-A.

```sh
set -eu
. /workspace/ed2-env/activate.sh
repo=/workspace/ED2
stage=/workspace/new-phase1e-o2
mkdir "$stage"
git -C "$repo" archive HEAD | tar -x -C "$stage"
cd "$stage"
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp' \
  LDFLAGS='-L/workspace/ed2-env/tools/usr/lib/x86_64-linux-gnu/openblas-pthread'
make -j4
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
make check-precision-contract PRECISION_OUTPUT=/workspace/new-phase1e-precision
make check-lanczos LANCZOS_OUTPUT=/workspace/new-phase1e-stability
make check-lapack
make check
make check-sparse-field SPARSE_FIELD_OUTPUT=/workspace/new-phase1e-fields
python3 -B tests/numerical_baseline.py --exe "$stage/source/QS3ED2" \
  --output /workspace/new-phase1e-independent
```

Use `FCFLAGS='-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra'` for bounds;
omit FCFLAGS for the GNU default `-fopenmp -Ofast`. For GNU 13 in this Cloud
instance use `FC=/workspace/ed2-env/gnu13/usr/bin/gfortran-13` with O2 flags and
append `-L/workspace/ed2-env/tools/usr/lib/gcc/x86_64-linux-gnu/14` to LDFLAGS for
the shared available Fortran runtime. Never reuse module files across compilers.
No profile uses `-fallow-argument-mismatch` or blanket floating-point traps.

For the seeded CLI probes, from the checkout (the builder creates new external
copies and verifies patch anchors):

```sh
python3 -B - <<'PY'
from pathlib import Path
import sys, json
sys.path.insert(0, 'tests')
import krylov_benchmarks as kb
kb.REVISIONS = {'phase1c': 'b56fe21678e62b8b6216fe0ec742d5dbd62a8df6'}
out = Path('/workspace/new-phase1e-probe')
out.mkdir()
fc = 'gfortran'
flags = ('-O2 -fopenmp -L/workspace/ed2-env/tools/usr/lib/x86_64-linux-gnu/openblas-pthread '
         '-L/workspace/ed2-env/tools/usr/lib/gcc/x86_64-linux-gnu/14')
bins = kb.build(Path.cwd(), out, fc, flags)
(out/'provenance.json').write_text(json.dumps(dict(
    revision=kb.REVISIONS, compiler=fc, flags=flags,
    executable_sha256=kb.sha(bins['phase1c']),
    probe_sha256=kb.sha(Path('tests/krylov_benchmark_probe.f90'))), indent=2)+'\n')
PY
python3 -B tests/krylov_acceptance_probe.py \
  --exe /workspace/new-phase1e-probe/phase1c/source/bench \
  --output /workspace/new-phase1e-acceptance
```

Repeat with the GNU 13 executable path and fresh directories. The existing
builder label `phase1c` selects the appropriate module list; the actual archived
revision here is the verified Phase 1-E starting commit, not an older algorithm.
The six analytic controls run automatically with each precision characterization.
Re-run with a fresh output and compare report JSON to check fixed-state replay;
timings are excluded by design. A reused output path is rejected.

The [evidence index](phase1e-precision/README.md) identifies the actual run layout
and reduction command. No binaries, full vectors or generated Autotools files
are committed. Independent QR/dense references share NumPy/backend arithmetic;
analytic diagonal/degenerate controls supplement them. Floating-point gap/bound
assessments are not interval-arithmetic certificates. No exhaustive starts,
large-system scalability, independent compiler family, MKL/ILP64 or public-library
ABI/reentrancy guarantee is established. The specification's open choices and
sequenced implementation gates remain for a separately authorized numerical phase.
