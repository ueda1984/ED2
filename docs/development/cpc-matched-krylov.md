# CPC matched one-eigenpair residual and resource study

**Status:** verified numerical characterization, **not** a public release
or a general performance claim.

## Why a matched comparison is necessary

The earlier ED2 mixed-spin scaling validation used ALG=1 to find
**one** low eigenpair, while ALG=2 returned **two**. Comparing their
raw wall times as an algorithmic speedup would be inappropriate.

The separate
[\`tests/mixed_spin_matched_benchmark.py\`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-alg2-restart-sweep/tests/mixed_spin_matched_benchmark.py)
therefore requests exactly **one** lowest eigenpair from **both** solvers
using the same model, sector dimension, explicit starting vector,
random seed and requested absolute residual threshold.

ALG=1 uses the current three-term Lanczos plus candidate reconstruction.
ALG=2 uses the current restarted Hermitian Krylov solver with
\`NOE=1, NOK=10, NOM=64\`. Although the *first* starting vector is
identical, ALG=2 can append additional deterministic seeded random
directions, so the full internal search spaces are not identical.

The test checks success/failure status, direct residual acceptance,
independent SciPy lowest-state energy and physical observables.
Genuine nonconvergence is captured as such rather than treated as a
successful energy/observable result. Each test occurs in a fresh
directory and retains its own \`input.dat\`, initial state, telemetry
and output.

## Independent verified results

[GitHub Actions: 10-site 16-run study](https://github.com/ueda1984/ED2/actions/runs/43397571923)

[GitHub Actions: extended 10-/12-site 20-run study](https://github.com/ueda1984/ED2/actions/runs/43397893641)

GNU Fortran 13, Ubuntu 24.04, OpenBLAS, OMP=1.
\`-O2 -fopenmp\`; both binaries built from a separately unpacked
source distribution.
20/20 trials **accepted with accurate low energy and observables**;
zero internal test-contract violations.
The 10-site measurements were repeated twice at each tolerance.
The 12-site measurements were single observations, not timing
estimates with uncertainty bars.

### 10 sites, 7,776-dimensional full space

| Requested absolute residual | ALG=1 H calls | ALG=1 measured residual | ALG=2 H calls | ALG=2 measured residual |
| ---: | ---: | ---: | ---: | ---: |
| 1e-6 | 140 | 5.157e-9 | 1,625 | 8.226e-15 |
| 1e-8 | 140 | 5.157e-9 | 1,625 | 8.226e-15 |
| 1e-10 | 230 | 8.305e-11 | 1,625 | 8.226e-15 |
| 1e-12 | 330 | 8.342e-13 | 1,625 | 8.226e-15 |

The two runs at each combination had identical solver residuals,
H-application counts and returned observables; timings varied modestly.
For \`atol=1e-8\`, ALG=1 took approximately 0.50 s and ALG=2
approximately 5.22 s on this runner.
Peak ED2 child-process memory was around 22.2 MiB versus 44.8 MiB.

### 12 sites, 46,656-dimensional full space

| Requested absolute residual | ALG=1 H calls | ALG=1 residual | ALG=2 H calls | ALG=2 residual |
| ---: | ---: | ---: | ---: | ---: |
| 1e-8 | 260 | 5.142e-10 | 2,470 | 1.189e-14 |
| 1e-10 | 370 | 1.077e-11 | 2,470 | 1.189e-14 |

At \`atol=1e-8\`, the recorded wall times were 4.259 s (ALG=1)
and 62.838 s (ALG=2). At \`atol=1e-10\` they were 4.947 s and
60.277 s. The readings are single CI observations; changing loads,
memory allocation and hardware can affect them.
The peak RSS for the matched 12-site runs was approximately
117 MiB (ALG=1) and 248 MiB (ALG=2).

## Mathematical and engineering interpretation

1. **Accuracy:** every accepted low eigenpair was independently
   matched to a sparse-reference lowest eigenvalue, not just accepted
   on the basis of a small solver residual.
2. **Additional residual policy:** for the tested ALG=2 case,
   demanding a tighter residual between \`1e-6\` and \`1e-12\`
   did not change its H-application count: the returned residual was
   already approximately \`8e-15\` at 10 sites. This is **consistent
   with** legacy energy convergence/restart closure dominating the
   work, but should not be taken as proof of the only bottleneck.
3. **Fair comparison:** the eigenpair *count*, explicit start and
   requested threshold now match; the methods still have different
   workspaces and internal seeded augmentation. An equal stopping
   threshold does not imply equal achieved residual.
4. **Possible tuning opportunities:** measure \`NOK\`, \`NOM\`,
   reorthogonalization and workspace/memory cost separately. Any
   proposed algorithmic optimization should preserve the existing
   regression suite and compare *achieved* physical accuracy, not
   merely returned status.
5. **No universal claim:** the study concerns one non-U(1)
   alternating-spin model with specified couplings; it is not a
   universal comparison of Lanczos versus thick restart.

## Reproduction

\`\`\`sh
python3 tests/mixed_spin_matched_benchmark.py \
  --exe /absolute/path/to/QS3ED2 \
  --output /new/empty/output-directory \
  --sizes 10,12
\`\`\`

The script requires test-only NumPy and SciPy. It writes a
machine-readable \`report.json\`, per-case inputs, explicit initial
vector, diagnostic sidecars and execution records.
The original GitHub Actions run is the provenance for the numeric
values cited above.

**Release caution:** no change to ED2's default solver policy,
public API, numerical kernels, \`main\`, or the upstream repository
was necessary for this test.
