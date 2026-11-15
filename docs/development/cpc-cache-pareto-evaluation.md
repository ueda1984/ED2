# CPC: full tight cache versus half-width generic cache

**Validation complete: 108/108 independent physics checks accepted.**
The [paired cache CI workflow #43704096102](https://github.com/ueda1984/ED2/actions/runs/43704096102)
completed all four jobs successfully, using source SHA
`08037050aac4cdaa0047c8c977edaaf9fc08ee46`.
This is an exploratory resource study on GitHub-hosted
Ubuntu 24.04 shared runners, not a dedicated-node speedup claim.

The calculation remains the same size-stationary **non-U(1)**
alternating spin-1/2/spin-1 open chain, with XYZ anisotropy,
symmetric off-diagonal Gamma_xy, site-dependent transverse
and longitudinal fields, and oriented staggered DM interaction
at D=+0.12. Full Hilbert-space dimensions are 1296, 7776
and 46656 at L=8,10,12. Every run independently verifies
ground energy, all Cartesian local magnetizations, all nine
nearest-neighbor Cartesian correlations and the requested
absolute direct eigenpair residual threshold of 1e-8.

## 1. Why this head-to-head comparison matters

Previous scans altered the number of representative states
cached (`wk_dim`), while [the exact MNTE study](cpc-mnte-capacity-evaluation.md)
showed that the stored-row capacity can also be made much smaller.
Here we compare **within each same CI job**, at identical physical
target, explicit starting vector, seed, residual policy,
site length, and solver/workspace settings:

| Policy | `wk_dim` | `MNTE` | Meaning |
| --- | --- | --- | --- |
| Half / generic | D/2 | -1 | Cache first half of rows at default broad row width |
| Full / tight | -1 | 7L/2-1 | Cache every row, but at the proven exact row bound |
| Full / generic | -1 | -1 | Cache every row at default broad row width |

This is a **program-configuration** comparison; there are
no changes to the Hamiltonian kernel. The two iterative
algorithms are evaluated separately with one ground
eigenpair requested. ALG=2 uses NOK=2, NOM=16.

## 2. Exact expected ordering of cache-array memory

With `NO_one=L` and `NO_two=L-1`, the broad default
row allocation is \(10L-7\) entries, including one spare
Fortran array index. The model-specific exact row
bound is \(7L/2-1\), which allocates \(7L/2\)
entries per cached row.

On the verified GNU Fortran configuration each cached
entry consists of a 4-byte integer and 16-byte
complex number. With \(D=6^{L/2}\), the nominal
allocation is therefore

\[
\begin{aligned}
B_{\rm half,generic} &=
20(10L-7)\left\lfloor\frac{D}{2}\right\rfloor,\\
B_{\rm full,tight} &=
20\frac{7L}{2} D,\\
B_{\rm full,generic} &=
20(10L-7)D.
\end{aligned}
\]

For all three even lengths,
\(B_{\rm full,tight}<B_{\rm half,generic}<B_{\rm full,generic}\).
This is an algebraic *allocation* comparison and does
not require observing wall-clock runtime.

For example at L=12, the two cache arrays allocate:

| Policy | Cache arrays (MiB) |
| --- | ---: |
| Half / generic | 50.279 |
| Full / tight | 37.375 |
| Full / generic | 100.558 |

Thus full/tight uses **less cache memory than half/generic**,
but all 46,656 representative states are cached.

## 3. Observed L=12 resource results

Each case is run once as a warmup and twice for measurement;
the table gives medians of the two measured runs.
All cases use one BLAS thread, fixed RNG seed, one
requested ground eigenpair, initial vector fixed
across policies, residual request 1e-8 and
GNU Fortran 13 `-O2 -fopenmp`.

| OMP | ALG | Cache policy | H applications | Wall s | Peak RSS KiB |
| ---: | ---: | --- | ---: | ---: | ---: |
| 1 | 1 | half / generic | 140 | 21.8044 | 67,970 |
| 1 | 1 | full / tight | 140 | 1.2973 | 54,886 |
| 1 | 1 | full / generic | 140 | 1.5198 | 119,516 |
| 1 | 2 | half / generic | 459 | 71.1435 | 97,448 |
| 1 | 2 | full / tight | 442 | 3.6522 | 84,238 |
| 1 | 2 | full / generic | 442 | 4.4056 | 148,914 |
| 2 | 1 | half / generic | 140 | 14.2540 | 67,944 |
| 2 | 1 | full / tight | 140 | 0.8920 | 54,944 |
| 2 | 1 | full / generic | 140 | 1.1094 | 119,516 |
| 2 | 2 | half / generic | 425 | 43.8450 | 97,472 |
| 2 | 2 | full / tight | 442 | 3.5102 | 84,322 |
| 2 | 2 | full / generic | 442 | 4.3475 | 149,054 |

**Qualitative result:** in all listed tests the full/tight
case uses **less observed peak process RSS and less
wall time** than half/generic. For ALG=1, the
H-application count is also identical across policies,
supporting a clean measurement of caching rather than
a change in Krylov iteration count. The measured
full/tight advantage over half/generic at L=12,
OMP=1 is about 16.8 in wall time for ALG=1
and 19.5 for ALG=2. These are **two-repeat
same-job CI observations**, not a statistically
established, hardware-independent speedup factor.

For ALG=2, the number of H applications can differ
across policies (459 versus 442 at OMP=1;
425 versus 442 at OMP=2). Achieved residuals
and physical correctness are checked separately.
Even matched requested convergence does not
force identical iterative trajectories.

### Independent memory accounting

At L=12, ALG=1, OMP=1, the measured RSS difference
between half/generic and full/tight is 13,084 KiB,
or approximately **12.777 MiB**. The algebraic
difference between the two cache array allocations
is approximately **12.904 MiB**. The close
agreement independently supports the explanation
based on the actual Fortran cache dimensions.

## 4. Complete workload matrix and provenance

All **108/108** cases passed, including 36 warmups
and 72 independent timed runs: L=8,10,12;
OMP=1/2, ALG=1/2, three cache policies,
one warmup and two measurements per workload.
Zero numerical nonconvergence and zero contract violations
were reported. Individual inputs, solver logs, output
observables and independent reference results are
stored in GitHub Actions artifacts with 21-day retention.
To avoid losing the trial measurements when artifacts
expire, the [full 108-row CSV](cpc-cache-pareto-raw-ci.csv)
is version-controlled, including run/job IDs,
source SHA, H applications, attained residual,
wall time and ED2 child maximum RSS.

The `tests/plot_cpc_pareto.py` script performs
a strict 108-row audit, excludes warmups,
generates machine-readable medians and
produces separate ALG=1/2 SVG/PNG exploratory
plots of wall time against peak RSS.
Do not use these plots as controlled-node
speedup figures without repeated pinned-node
measurements.

## 5. Interpretation for the manuscript

This result suggests a **physics-aware cache allocation
strategy** for generalized quantum-spin ED:
derive the maximum needed row capacity from the
allowed local transition channels, then cache all rows
when memory permits. The present case combines
a mathematically exact per-row bound with independent
matrix/observable verification. It does **not**
establish a new general algorithm or change ED2's
pre-existing `wk_dim`/`MNTE` controls.

The next production benchmark should repeat the
three policies under fixed, isolated CPU affinity
and fixed BLAS, record hardware/cache/compiler
information and multiple independent timing runs.
It should also test a Hamiltonian with a different
transition structure before suggesting an
automatic default. The current tests leave
`source/`, `main` and the upstream repository
untouched.
