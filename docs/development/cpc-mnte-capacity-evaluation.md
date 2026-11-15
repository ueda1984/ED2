# CPC: exact row-capacity bound for the stationary mixed-spin DM example

**Verified numerical resource characterization; not a new eigensolver algorithm.**
The program already supports explicit `MNTE` and prints its observed
`Optimal MNTE`. Here we combine a model-specific upper-bound derivation,
independent physical validation and a measured cache-memory audit.

## Physics and provenance

This experiment uses the same **open**, alternating spin-1/2/spin-1
Hamiltonian as the [main stationary DM response](cpc-stationary-dm-application.md).
It preserves the two-site-periodic non-U(1) XYZ exchange, symmetric
Gamma_xy, transverse/longitudinal fields and oriented staggered DM term
at D=+0.12. The full Hilbert-space dimension is 6^(L/2), no fixed-Sz
projection or symmetry truncation.

- [CI execution 43703623648](https://github.com/ueda1984/ED2/actions/runs/43703623648),
  source SHA `315bed65b0c21d4db4f3af834d5807d361839d29`.
- [Additional independent sparse-row verification #43704146968](https://github.com/ueda1984/ED2/actions/runs/43704146968)
  independently constructs the SciPy CSR Hamiltonian and confirms
  that its maximum nonzero count per row is exactly 27, 34, 41,
  respectively. All four CI jobs passed for this stronger check.
- GitHub Actions checks: **72/72 cases accepted** (L=8,10,12;
  ALG=1/2, OMP=1/2, default/calibrated capacity; one warmup and
  two measured runs for each combination). All solver statuses,
  direct residual requests, independently constructed low energies,
  three local magnetizations per site and nine nearest-neighbor
  Cartesian correlations passed.
- [Frozen raw 72-case CSV](cpc-mnte-capacity-raw-ci.csv) contains
  machine/job IDs, source SHA, capacity, measured maximum,
  H applications, elapsed time and peak ED2 RSS. Warmup rows are
  retained but excluded from all reported medians.
- The script `tests/cpc_stationary_mnte_capacity.py` measures the
  maximum row length with a complete `wk_dim=-1` cache and
  verifies the calibrated value in fresh independent runs.
- GitHub-hosted runner timing is exploratory and **not** a
  portable speedup or an isolated-node benchmarking result.

## 1. Why a tighter MNTE is possible

In the existing `source/ham2vec.f90`, with `MNTE=-1`, the
cached arrays have `(M,wk_dim)` entries, with

\[
M=2\,\mathrm{NO}_{\mathrm{one}}+
  8\,\mathrm{NO}_{\mathrm{two}}+1 .
\]

For this open chain, `NO_one=L` and `NO_two=L-1`.
The **default input capacity** printed by ED2 is therefore

\[
\mathrm{MNTE}_{\rm default}=M-1=10L-8 .
\]

This is a safe allocation bound for more general local interaction
components. It is not the measured maximum number of nonzero
transitions in the present model.

For the **specific** open alternating spin-1/2/spin-1 chain studied
here, all `S^x` transverse fields and the Gamma_z and XY/DM
two-spin-flip terms are nonzero. On a single product state:

1. Every spin-1/2 site admits at most **one** outgoing
   single-site transverse flip: L/2 channels.
2. Every spin-1 site admits at most **two** single-site
   transverse flips (at m=0): L channels.
3. Each of the L-1 nearest-neighbor mixed-spin bonds
   admits at most **two** simultaneous two-site
   transverse flips: 2(L-1) channels.
4. The diagonal contribution is stored in **one** entry.

The complete product basis contains states with all spin-1
sites at m=0. In those states the upper bound is attained
because the relevant transverse and two-spin-flip coefficients
do not vanish. Hence the exact maximum number of generated
entries per representative state for these couplings is

\[
\boxed{\mathrm{MNTE}_{\rm model}=
1+\frac L2+L+2(L-1)=\frac{7L}{2}-1.}
\]

This **derivation is conditional on the model class**:
it is not a globally safe default if one adds other DM
components, Gamma_x/Gamma_y terms, longer-range bonds,
site-dependent spin changes or symmetry reductions.
Different product-basis restrictions also require revalidation.

## 2. Observed maximum agrees exactly with the formula

| L | Full dimension | Default MNTE | Derived MNTE | Observed optimal MNTE |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 1,296 | 72 | 27 | 27 |
| 10 | 7,776 | 92 | 34 | 34 |
| 12 | 46,656 | 112 | 41 | 41 |

All four combinations of ALG=1/2 and OMP=1/2 confirm
these exact row-length maxima at every L, including
the calibrated reruns. The proof addresses the channel count;
the independent SciPy Hamiltonian's maximum CSR-row nonzero count
and ED2's own `Optimal MNTE` agree separately with the formula.
This protects against accidentally validating only two code paths
that share the same transition enumerator.

## 3. Predicted reduction in allocated arrays versus RSS

The current GNU Fortran build uses default INTEGER(4) and
COMPLEX(8) storage for `wk_loc` and `wk_ele`.
Because `wk_dim=-1` caches the complete working space,
changing only MNTE from M0 to M1 saves precisely

\[
\Delta B_{\mathrm{cache}}=
(4+16)(M_0-M_1)\,6^{L/2}\ {\rm bytes}.
\]

Total RSS also includes basis lists, iterative solver vectors,
other arrays and runtime overhead. The following comparison
uses the **ALG=1, OMP=1, two measured-run medians**.

| L | Predicted cache saving (MiB) | Observed peak-RSS reduction (MiB) |
| ---: | ---: | ---: |
| 8 | 1.112 | 1.100 |
| 10 | 8.602 | 8.619 |
| 12 | 63.182 | 63.080 |

The magnitude and scaling agree with the known Fortran
storage layout. At L=12, reducing MNTE from 112 to 41
reduces measured process RSS from 119,452 KiB to
54,858 KiB for ALG=1, and from 148,904 KiB to
84,178 KiB for ALG=2 (OMP=1). The Hamiltonian-action
counts remain 140 and 442 respectively.

### Full versus half cache, with and without tight MNTE

An instructive, **analytical** comparison of the cache arrays
uses the L=12 settings:

- Default full cache: 100.558 MiB in the two cache arrays.
- Default half cache: 50.279 MiB.
- Model-derived full cache: (41+1)/113 × 100.558 = **37.375 MiB**.

Thus a **full cache with a calibrated row bound requires less
array memory** than the default-width *half cache* in this
particular model, while avoiding on-the-fly Hamiltonian
reconstruction for all rows. Any comparison of measured runtimes
between those two policies must still use an appropriately
controlled same-runner or dedicated-node experiment, rather
than combining results from unrelated CI runners.

## 4. Observed timing and numerical invariants

| L | ALG | OMP | Default wall median (s) | Calibrated wall median (s) | H applications |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1 | 1 | 0.034552 | 0.031868 | 120 |
| 8 | 2 | 1 | 0.080853 | 0.076712 | 306 |
| 10 | 1 | 1 | 0.301079 | 0.265619 | 200 |
| 10 | 2 | 1 | 0.613823 | 0.547784 | 340 |
| 12 | 1 | 1 | 1.507066 | 1.300323 | 140 |
| 12 | 2 | 1 | 4.816554 | 3.799433 | 442 |

Every same-size, same-algorithm, same-OMP pair has an unchanged
Hamiltonian-action count and passes the independent physical
acceptance test. The small measured wall-time differences on shared
runners should be reported as observations only. The **memory
reduction is stronger evidence** because it closely matches the
source-level allocation formula.

**Figure prototypes:** use `tests/plot_cpc_mnte.py` with the
committed raw CSV to generate three SVG/PNG files:
`mnte_peak_rss`, `mnte_memory_prediction` and `mnte_time_ratio`.
The plotting script performs a strict frozen-data audit and writes
a generated median-summary CSV.

## 5. Scope and proposed manuscript claims

The result can support a CPC subsection on **physics-aware
resource tuning and reproducible Hamiltonian-action caching**.
It should not be presented as a new general computational
complexity result: the manual `MNTE` control already exists.
The specific contribution is an exact row-capacity derivation
for a representative nonconserving mixed-spin Hamiltonian,
independently checked against the software's enumeration and
observable-level correctness tests.

Recommended next experiment: on one pinned, dedicated node,
compare **default half cache** directly against **tight full
cache** at the same full-space target and residual, and record
median runtime/RSS, uncertainty and achieved eigenpair residual.
Do not change source numerical kernels or default values based
on this limited model alone.
