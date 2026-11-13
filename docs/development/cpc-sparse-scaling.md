# CPC: independent sparse reference and preliminary resource scaling

**Status:** independently verified scientific scaling **on one GitHub runner**;
not a tuned performance benchmark or a release claim.

## Model and exact-space meaning

The source code [\`tests/mixed_spin_sparse_scale.py\`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-sparse-scaling/tests/mixed_spin_sparse_scale.py)
constructs an independent complex Hermitian Hamiltonian using SciPy
sparse matrices in the ordinary spin-product basis and computes its
lowest two eigenpairs with ARPACK \`eigsh(which="SA")\`.

The model is the same open alternating \(S=1/2,1\) chain as in
[the small-system dense ED study](cpc-mixed-spin-demo.md), with
nonuniform anisotropic exchange, oriented \(D_z\), symmetric
\(\Gamma_{xy}\), and site-dependent transverse/longitudinal fields.
The magnetic term uses **\( + \mathbf h_i\cdot\mathbf S_i\)**.

The complete local product space is represented:
\(D=\prod_i(2S_i+1)=6^{N/2}\), with
\(\mathrm{NODmin}=0\) and
\(\mathrm{NODmax}=\sum_i 2S_i\).
No symmetry reduction or projected truncation is used in this scaling
study. The non-U(1) Hamiltonian is *not* an invariant fixed-magnetization
block.

The comparison validates the lowest energies (ALG=1: one, ALG=2: two),
local magnetization in all Cartesian directions, all nine Cartesian
pair correlations for each distinct pair, and the optional direct
residual criterion \(\|H\psi-\lambda\psi\|_2\le10^{-8}\).
The independent sparse eigenpair residuals are checked directly, and
the lowest two energies have a nonzero finite-size gap.

**Important distinction:** ARPACK gives a reproducible independently
constructed numerical reference, not a proof of global ground-state
inclusion for every hypothetical Hamiltonian. The 4/6-site cases were
separately compared against a full dense spectrum in
[the small-system study](cpc-mixed-spin-demo.md).

## Verified OMP=1 results

[Primary GitHub Actions execution](https://github.com/ueda1984/ED2/actions/runs/43396341796)
(13 November 2026; Ubuntu 24.04, GNU Fortran 13, \`-O2 -fopenmp\`,
OpenBLAS one thread, OpenMP one thread).
An **unrelated extracted source archive** was compiled for ED2;
the physics reference uses SciPy, not ED2's Hamiltonian construction.

| Sites | Hilbert \(D\) | Sparse \(H\) nonzeros | Independent \(E_0\) | Independent \(\Delta=E_1-E_0\) |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 1,296 | 25,480 | -4.922896708360970 | 0.210915273554495 |
| 10 | 7,776 | 191,808 | -6.284268018540742 | 0.180644980189410 |
| 12 | 46,656 | 1,384,044 | -7.660800376520366 | 0.136192614266859 |

These rows use the transverse-field control value 0.15. A second
8-site control value 0.38 was also verified.
**All eight ALG=1/2 cases passed** independent energy,
magnetization, correlation and residual assertions.

| Sites | ALG | Requested eigenpairs | Max energy error | Max residual | \(H\) applications | Wall time (s) | Max RSS (KiB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1 | 1 | 2.04e-14 | 1.93e-09 | 120 | 0.068 | 8,620 |
| 8 | 2 | 2 | 2.22e-14 | 5.86e-15 | 792 | 0.264 | 13,764 |
| 10 | 1 | 1 | 2.40e-14 | 9.49e-09 | 140 | 0.297 | 22,240 |
| 10 | 2 | 2 | 2.31e-14 | 7.09e-15 | 1,782 | 3.655 | 45,996 |
| 12 | 1 | 1 | 8.08e-14 | 4.39e-10 | 260 | 3.888 | 118,848 |
| 12 | 2 | 2 | 7.99e-14 | 1.59e-14 | 2,640 | 56.452 | 254,100 |

**At 12 sites**, the maximum magnetization discrepancy was
\(\le 7.53\times10^{-11}\) for ALG=1 and \(\le 2.34\times10^{-14}\)
for ALG=2; the correlation discrepancies were
\(\le1.05\times10^{-10}\) and
\(\le1.82\times10^{-14}\), respectively.

Memory is the peak resident set of the **ED2 child process**
measured by GNU \`/usr/bin/time\`. Wall-clock includes process startup,
the eigenproblem and observable output. Values are individual samples
on shared CI hardware, **without repeated timings, uncertainty bars
or guaranteed dedicated CPUs**.

### Why these numbers do *not* justify an ALG1-versus-ALG2 speed claim

ALG=1 requests one eigenpair; ALG=2 requests two. Both are asked for
absolute residual \(\le10^{-8}\), but ALG=2's actual residual is
orders of magnitude smaller. Their convergence histories, Krylov
subspace dimensions and memory requirements are different.

Therefore:
- do **not** describe the time or memory ratios as a fair same-workload
  speedup or algorithmic superiority;
- compare **matched requested eigenpair counts** and tolerances before
  publishing any numerical advantage claim;
- separate executable wall time from the H-application counts and
  independent reference calculation resources.

The present observations do, however, show clearly why
**restart/workspace cost and the number of requested eigenpairs**
are first-class benchmark parameters.

## Current evidence and next scientific questions

- Small full-space 4- and 6-site cases: all eigenvalues and observables
  cross-checked against independent dense NumPy ED.
- Extended 8-, 10-, 12-site cases: two lowest sparse-reference energies
  and ground-state local/two-site observables independently checked.
- Fixed seed and optional residual gate used for ALG=1/2; status zero
  was not accepted without external eigenvalue/observable comparisons.
- Same code path can be exercised with OpenMP 2 threads; correctness
  is a necessary but insufficient parallel-performance metric.
- A controlled **one-eigenpair versus one-eigenpair** study with
  matched requested and achieved precision, repeated timings and
  memory evidence is needed before drawing performance conclusions.
- A broader physical application or parameter sweep with
  interpretable size-dependence remains necessary before a strong CPC
  physics-performance claim.

The benchmark does **not** show that generic ED2 code is asymptotically
optimal, that every requested state is globally lowest for arbitrary
starts, or that a public stable Fortran library API exists.
