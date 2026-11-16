# CPC near-saturation cache comparison: explicit physical and memory contract

**Status:** ED2-native source-pinned benchmark and pilot GitHub Actions workflow
are introduced in this branch. A shared CI run is *exploratory only*.
Publication-grade wall-time claims require exclusive-node scheduler/operator
attestation, CPU frequency/load records and replicated measurements.

## Why the row width is as important as cache fraction

The ED2 Hamiltonian cache in `source/ham2vec.f90` consists of two arrays:
`wk_loc(M, wk_dim)`, integer(4), and `wk_ele(M, wk_dim)`,
complex(8). The allocated amount is **20 M wk_dim bytes**, excluding
the solver, basis and runtime allocations.

The generic `MNTE=-1` setting uses `M=2*NO_one+8*NO_two+1`.
For periodic cubic spin-1/2 NN Heisenberg with all longitudinal site
fields zero, `NO_one=L` and `NO_two=3L`, giving `M=26L+1`.
The current benchmark's **model-specific, not general** safe row width
is `M=6*Ndown+2`: there is one diagonal term and at most
six directed NN spin-flip moves per down spin, plus a diagonal term
and one extra slot required by the historical `jd` sentinel index.
This extra slot is essential: CI explicitly observed `Optimal MNTE=13`
for Ndown=2 when `MNTE=12` had been requested, and rejected the run.
Off-diagonal
matrix entries related by symmetry may coalesce and only decrease this
maximum. This bound must NOT be reused for XY anisotropy,
non-U(1) interactions, long-range exchanges or general mixed-spin models
without a fresh derivation and a test of `Optimal MNTE`.

| Lattice | Ndown | Reduced D | full/generic (MiB) | full/tight (MiB) |
| --- | ---: | ---: | ---: | ---: |
| 6x6x6 | 3 | 7677 | 822.48 | 2.93 |
| 10x10x10 | 2 | 503 | 249.5 | 0.134 |

These are **only nominal cache-array allocations**, not observed RSS.

## Modes and safety

The first pilot compares
- `none`: `wk_dim=0, MNTE=6*Ndown`
- `half_tight`: `wk_dim=D//2, MNTE=6*Ndown`
- `full_tight`: `wk_dim=-1, MNTE=6*Ndown`

This isolates the influence of cached fraction *without* conflating
memory consumed by generic row allocation. A separate full/generic
mode is supported for the 1000-site two-magnon test; the 216-site
full/generic case is barred by the runner because it would allocate
over 800 MiB of cache alone.

The benchmark records:
- identical physical input hash and deterministic complex initial vector
  hash for each cache policy; RNG seed 1984;
- completed source-pinned ED2 eigenpairs and independent reference E0;
- `ED2_RESIDUAL_ATOL=1e-8`, direct eigenpair residual, number of
  Hamiltonian applications, iteration count and cache's observed
  `Current/Optimal MNTE` values;
- full subprocess wall time and Linux `/usr/bin/time` maximum RSS;
- an explicit shared/nonpublication flag and a machine/CPU identity.

The numerical acceptance test checks that the residual bound is met,
the independent ground-state energy agrees within 2e-7 and the same
initial state and H application count are obtained across cache settings.
Any failed run invalidates the complete comparison and remains recorded.

## Running at a dedicated compute node

From an immutable checkout of the verified ED2 source, after compiling
`source/QS3ED2`, first confirm the scheduler grants exclusive access
to the selected CPUs. Determine the **actual** permitted CPUs with
`taskset -pc $$`. Do not assume CPU 0 is accessible.

Run the following within the dedicated, attested node allocation:

```bash
python3 -S tests/cpc_few_magnon_cache_benchmark.py \
  --exe "$(pwd)/source/QS3ED2" \
  --output "$PWD/few-magnon-cache-node" \
  --cases 1000-n2,216-n3 \
  --policies none,half_tight,full_tight \
  --warmups 2 --repeats 7 --timeout 900 \
  --residual-atol 1e-8 \
  --exclusive-allocation-verified
```

The flag records an operator attestation, **not automatic proof** of
exclusive scheduler resources. Save scheduler allocation logs alongside
`report.json`. To run on a workstation or shared CI, replace the flag
with `--allow-shared`; no shared-runner result can automatically pass
the publication gate.

The runner randomizes the cache policy ordering with a fixed reproducible
seed independently for each measurement round. One CPU is pinned
from the permitted affinity. `OMP_NUM_THREADS=1`, with BLAS limited
to one thread. For OpenMP scaling, extend the separate controlled-node
protocol; never mix thread counts within this first one-core comparison.

## Interpretation

The scientific question is whether the total practical cost
`T_basis + T_cache + T_solver` decreases when the transitions are
precomputed. It is **not** enough to divide two Krylov loop times
if one mode pays a large setup or allocation cost.

For publication, plot separately:
1. end-to-end time and peak process RSS for no/half/full tight cache,
2. total time divided by Hamiltonian applications (with setup separately
   instrumented before interpreting this as a per-matvec metric),
3. repeated medians with median absolute deviations and all raw runs.

The single-model independent reference cases in the original project are:
- 216-site 3-magnon Gamma E0 ~ 148.5648907577122; D=7677;
- 1000-site 2-magnon Gamma E0 ~ 738.0067742722034; D=503.

The original QS3 and QuSpin are reference engines for numerical
capability, but comparisons of their *caching policies* must be
separate from ED2-internal no/half/full toggles.

Scripts: `tests/cpc_few_magnon_cache_benchmark.py`;
CI: `.github/workflows/cpc-few-magnon-cache-benchmark.yml`.


## First validated shared-runner pilot (2026-11-16)

Source pinned to the [CI workflow run 43784564060](https://github.com/ueda1984/ED2/actions/runs/43784564060),
using GNU Fortran 13, `-O2 -fopenmp`, one OpenMP thread, one
CPU pinned from the runner's allowed affinity and BLAS=1. This is
**one measurement per mode, with no warmup, on a shared runner**,
so all wall-time ratios are preliminary and cannot be promoted to a
journal-grade, controlled-node speedup.

**All six cache configurations passed independently frozen ground
energy and the requested direct Ritz residual tolerance 1e-8.**
Within each physical model, all policies had the same initial-vector
SHA256 and number of Hamiltonian actions:

| Geometry | Ndown | Cache | Wall seconds | Peak RSS (KiB) | H actions | Direct residual |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| 10×10×10 | 2 | none | 8.487879 | 6,928 | 185 | 3.0812e-9 |
| 10×10×10 | 2 | half/tight | 4.796797 | 6,888 | 185 | 3.0813e-9 |
| 10×10×10 | 2 | full/tight | 0.996301 | 6,972 | 185 | 3.0803e-9 |
| 6×6×6 | 3 | none | 125.467950 | 8,176 | 545 | 2.6323e-9 |
| 6×6×6 | 3 | half/tight | 64.081302 | 9,716 | 545 | 2.6323e-9 |
| 6×6×6 | 3 | full/tight | 1.134819 | 11,132 | 545 | 2.6323e-9 |

The observed full-cache vs no-cache end-to-end time ratios were
**8.52** at 1000 sites/two magnons and **110.56** at 216 sites/three
magnons. They have **no uncertainty estimate** from one measurement.
The 216-site ground energy matched
148.5648907577122, and the 1000-site energy matched
738.0067742722034 within 2e-7 in every accepted run.

The **model-specific full/tight cache arrays** consumed
2.929 MiB for 216/3 (7,677 basis states) and 0.134 MiB
for 1000/2 (503 basis states); these are theoretical allocated
array sizes, not measured process RSS. The generic full-cache
allocation would be approximately 822.481 and 249.453 MiB,
respectively, under the standard integer(4)/complex(8) representation.
Cache fraction and cache row width must be optimized separately.

An initial pilot requested `MNTE=6*Ndown` and **correctly failed** the
safety gate because the Fortran-reported optimal row bound was one
greater: for two magnons `Current MNTE=12` but
`Optimal MNTE=13`. The corrected bound
`MNTE=6*Ndown+1` passed all tests. The failure is retained in
[the initial CI run 43784444333](https://github.com/ueda1984/ED2/actions/runs/43784444333)
rather than ignored.

### What the pilot cannot establish

1. Stable speedup factors or confidence intervals: there are no
   independent job allocations and only one observation per mode.
2. Relative speed against QuSpin/old QS3: this is an ED2 *internal*
   cache study, and those programs do not share the same cache-policy
   interface.
3. Pure (H|\psi\rangle) kernel performance: the reported end-to-end
   time includes input/basis/cache construction and solver execution.
   Isolating kernel throughput requires further opt-in timing.
4. Equal performance at arbitrary momentum, DM coupling or mixed
   spins: the narrow row bound is **specific** to periodic isotropic
   NN few-magnon Hamiltonians.

Recommended follow-up: repeat `--warmups 2 --repeats 7` with an
operator-verified exclusive node, record actual CPU topology and
scheduler exclusivity, and run at least two separate allocations.
Inspect outliers, H-call equality, direct residual, independent energy
and RSS before quoting speedup ratios in a submitted paper.
