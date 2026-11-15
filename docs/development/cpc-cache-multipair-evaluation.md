# CPC follow-up: cache policy and two-lowest-eigenpair robustness

**Status:** independently verified on GitHub Actions. These are **test-only**
resource characterizations on shared runners, not dedicated-node performance
claims. All tests retain the original size-stationary, non-U(1), alternating
spin-1/2 and spin-1 XYZ/Gamma/transverse-field chain with oriented staggered
DM interaction; the Hamiltonian is not reduced to XXZ.

Provenance:
- [Cache-policy workflow #43702889491](https://github.com/ueda1984/ED2/actions/runs/43702889491),
  SHA bdbb4402baeb842a80f5b39779150425ec9226b6.
- [Two-lowest-eigenpair workflow #43703014800](https://github.com/ueda1984/ED2/actions/runs/43703014800),
  SHA c89f018c2adef6af5c1647e83e2b4fed6dbfaa17.
- Scripts: `tests/cpc_stationary_cache_sweep.py` and
  `tests/cpc_stationary_multipair_benchmark.py`.
- Raw case inputs, diagnostics, independent-reference comparisons and
  resource telemetry are available in the Actions artifacts for 21 days.

All jobs used Ubuntu 24.04, GNU Fortran 13, `-O2 -fopenmp`,
OpenBLAS restricted to one thread, `OMP_DYNAMIC=FALSE`, full
Hilbert space, explicit ED2-basis starting vector, RNG seed 19,
requested absolute residual 1e-8, and independent SciPy references.
Each accepted calculation checks the lowest state and ground-state
magnetizations and all nearest-neighbor Cartesian correlations.

## Cache-policy dependence

`wk_dim=0` stores no cached rows; `wk_dim=-1` stores all rows.
Half caching stores the **first half** of representative rows, not a
randomly selected half. The physical task and accuracy requirement
are identical. L=8,10 results include one warm-up and two measured
runs per case.

All **36/36 small-system trials** passed at OMP=1/2 (including
warm-ups). The following values are two-run medians.

| L | OMP | Cache | H applications | Wall s | Peak RSS KiB |
| ---: | ---: | --- | ---: | ---: | ---: |
| 8 | 1 | none | 120 | 0.529700 | 6,646 |
| 8 | 1 | half | 120 | 0.278283 | 7,700 |
| 8 | 1 | full | 120 | 0.024038 | 8,604 |
| 10 | 1 | none | 200 | 7.511783 | 8,082 |
| 10 | 1 | half | 200 | 3.916686 | 15,218 |
| 10 | 1 | full | 200 | 0.230234 | 22,278 |
| 8 | 2 | none | 120 | 0.258016 | 6,716 |
| 8 | 2 | half | 120 | 0.143860 | 7,546 |
| 8 | 2 | full | 120 | 0.016639 | 8,442 |
| 10 | 2 | none | 200 | 3.751903 | 8,124 |
| 10 | 2 | half | 200 | 1.946044 | 15,128 |
| 10 | 2 | full | 200 | 0.127878 | 22,238 |

The H-application count is unchanged across cache policies for
ALG=1 at these two sizes. Full versus no caching on this runner
reduces the observed wall time by a factor of about **22.0** at
L=8 and **32.6** at L=10 (OMP=1), while raising memory use.
Neither ratio is a measured portable speedup.

### Implementation-level memory accounting

The Fortran arrays `wk_loc` and `wk_ele` in
`source/ham2vec.f90` have shape
`(M, wk_dim)`, where for default `MNTE=-1`
\[
M=2\,\mathrm{NO\_one}+8\,\mathrm{NO\_two}+1 .
\]
For an open L-site chain, `NO_one=L` and `NO_two=L-1`,
so \(M=10L-7\).
On the verified GNU Fortran build, the two arrays together
allocate \( (4+16)M\,\mathrm{wk\_dim}\) bytes
(default 4-byte integer and 16-byte `complex(8)`).
The nominal **full-cache array** sizes are:

| L | Full dimension | Cache-row width M | Two cache arrays (MiB) |
| ---: | ---: | ---: | ---: |
| 8 | 1,296 | 73 | 1.805 |
| 10 | 7,776 | 93 | 13.793 |
| 12 | 46,656 | 113 | 100.558 |

These are not total RSS estimates. They exclude the solver vectors,
state lists, runtime overhead and allocator effects. The measured
full-minus-none median RSS for L=10, OMP=1 is 14,196 KiB (13.86 MiB),
close to the analytical 13.79 MiB cache allocation.

### L=12 OMP=2 pilot

All **8/8** trials passed independent checks (two repeats of half
and full cache, for ALG=1 and ALG=2 with NOK=2,NOM=16).

| ALG | Cache | H applications | Median wall s | Peak RSS KiB |
| --- | --- | ---: | ---: | ---: |
| 1 | half | 140 | 8.865887 | 67,870 |
| 1 | full | 140 | 1.093300 | 119,376 |
| 2 (2/16) | half | 425 | 27.390648 | 97,454 |
| 2 (2/16) | full | 442 | 3.779526 | 149,056 |

The full-minus-half RSS difference for ALG=1 is 51,506 KiB,
close to the predicted 50.279 MiB cache difference.
ALG=2 has different H-application counts across these
cache policies: equal input and requested residual do
**not** guarantee an identical iterative trajectory.
L=12 OMP=1 evidence must be reviewed separately before adding it.

## Two-lowest-eigenpair ALG=2 workload

To avoid generalizing a single-ground-state workspace optimum,
we request **two lowest eigenpairs** with the same initial vector,
seed, residual criterion, full space and physical observables.
Both low eigenvalues and direct residuals are verified independently,
as are ground-state magnetizations and Cartesian correlations.
This is an ALG=2 workspace comparison, **not** an ALG=1-versus-ALG=2
two-eigenpair timing comparison.

The complete [workflow #43703014800](https://github.com/ueda1984/ED2/actions/runs/43703014800)
passed **36/36** trials: 32 repeated L=8/10, D=0/0.12,
OMP=1/2 measurements, and four single L=12, D=0.12 measurements.
No failed reference check or nonconvergent run was reported.

| L | D | OMP | NOK/NOM | H applications | Largest residual |
| ---: | ---: | ---: | --- | ---: | ---: |
| 8 | 0.12 | 1 | 2/16 | 324 | 4.38e-10 |
| 8 | 0.12 | 1 | 4/24 | 520 | 9.55e-14 |
| 10 | 0.12 | 1 | 2/16 | 360 | 8.62e-09 |
| 10 | 0.12 | 1 | 4/24 | 572 | 2.39e-10 |
| 12 | 0.12 | 1 | 2/16 | 540 | 4.19e-09 |
| 12 | 0.12 | 1 | 4/24 | 754 | 2.52e-11 |
| 12 | 0.12 | 2 | 2/16 | 594 | 4.64e-09 |
| 12 | 0.12 | 2 | 4/24 | 780 | 1.41e-11 |

The 2/16 setting requires fewer H applications in these
particular cases, but its achieved residual approaches the
1e-8 request at L=10/12. The 4/24 setting requires more work
and gives a smaller residual. Both are within their requested
accuracy contracts, including the independent physical checks.
These are not near-degenerate-subspace tests; NOE>2 and
near-degenerate projector/ensemble controls remain necessary.

## Proposed CPC presentation and remaining experiments

The resource figure should show wall time and peak RSS
versus cache fraction with an explicit same-target caption,
alongside a table of H applications, achieved residuals
and the memory-allocation formula. Controlled-node repetitions
(ideally at least five, with warm-up excluded) and pinned CPU/BLAS
affinity are required for credible speed comparisons. Measure
compilers/backends, memory, CPU caches and version/source checksums.
Evaluate the effect of the `MNTE` upper bound separately rather
than conflating it with `wk_dim`.

The present data demonstrate a **real implementation-level
time/memory tradeoff for a nonconserving mixed-spin Hamiltonian**,
not a new physical phase, universal solver superiority or
algorithmic asymptotic improvement.
