# CPC follow-up: published XYZ-DM Hamiltonian and matched stationary-resource characterization

**Status: GitHub Actions verified, 2026-11-15.** [Run 43702189530](https://github.com/ueda1984/ED2/actions/runs/43702189530), source SHA
97101cfb8b47dfe04b0bbdf5d7e1a329ea7a4a61 at the time of the run.
The independent physical test and the resource benchmark ran from a
**separately unpacked source archive**, with GNU Fortran 13,
-O2 -fopenmp, Ubuntu 24.04 and single-thread BLAS.
Individual inputs, ED2 outputs, solver telemetry, comparison reports
and selected timing metadata are retained as GitHub Actions artifacts
for that run. There are **no changes to source/** in this study.

## A. Published Hamiltonian: what is and is not reproduced

Yi, You, Wu and Oles, *Physical Review B* **100**, 024423 (2019),
https://doi.org/10.1103/PhysRevB.100.024423, Eq. (1), use
**Pauli** matrices, periodic boundaries, an XYZ exchange anisotropy,
a uniform directed z-component DM term and a negative Zeeman field:

\[
H=J\sum_i\left[
\frac{1+\gamma}{2}\sigma_i^x\sigma_{i+1}^x+
\frac{1-\gamma}{2}\sigma_i^y\sigma_{i+1}^y+
\Delta\sigma_i^z\sigma_{i+1}^z\right]
+D\sum_i(\sigma_i^x\sigma_{i+1}^y-\sigma_i^y\sigma_{i+1}^x)
-h\sum_i\sigma_i^z .
\]

For ED2's S=sigma/2 and +hvec dot S convention, the exact mapping is
Jx=2J(1+gamma), Jy=2J(1-gamma), Jz=4J Delta, Dz=4D, hz=-2h.
The closing periodic bond (L,1) is included as an **oriented** bond.
Unlike open-chain tests, this probes the sign correction applied to
a reversed input bond and the reversed output orientation of xy and yx.

The deterministic test points are J=1, gamma=0.37, Delta=0.41,
h=0.58, D=0 and D=0.19, and L=4,6,8. They were **not taken from
figure coordinates in the publication**. The calculation is a
**Hamiltonian-level finite-size reproduction and independent ED2
verification**, not a published phase-diagram or numerical-dataset
reproduction. The original 4–12-site alternating-spin stationary
DM demonstrator remains the **primary** scientific example.

The independent NumPy Pauli Hamiltonian is verified against a second
spin-operator representation with the mapped ED2 coefficients.
All eigenvalues, local magnetizations, all nine Cartesian correlations,
the oriented DM-conjugate chirality, a finite-difference
Hellmann–Feynman check and non-U(1) commutator are examined.

| L | D | Lowest energy from independent Pauli H | First gap |
| ---: | ---: | ---: | ---: |
| 4 | 0.00 | -3.854733105562 | 0.478797061 |
| 4 | 0.19 | -3.854733105562 | 0.478797061 |
| 6 | 0.00 | -5.466704725748 | 0.006836479 |
| 6 | 0.19 | -5.468567485505 | 0.008699239 |
| 8 | 0.00 | -7.368761434736 | 0.172597665 |
| 8 | 0.19 | -7.369394757989 | 0.170241793 |

**Result:** 6/6 cases for OMP=1 and 6/6 for OMP=2 accepted;
no reported validation violations. The largest full-spectrum
difference from ED2 was below 3e-14, while the local magnetization,
Cartesian correlation and DM-conjugate expectation errors were
below 1e-14. The independently verified mapping difference is zero
at numerical precision. This validates Hamiltonian conventions,
not a previously unknown physical phase.

## B. Same-workload resource measurement on the primary mixed-spin model

The second test uses the stationary two-site-rule alternating
spin-1/2,spin-1 chain already reported in
[cpc-stationary-dm-application.md](cpc-stationary-dm-application.md).
It retains non-U(1) XYZ exchange, Gamma_xy, transverse/longitudinal
two-site-periodic fields and staggered oriented DM, at **D=0.12**.
Thus we do **not** replace the primary science case by a
spin-conserving or simpler computational model.

ALG=1 and ALG=2 both request the **one lowest eigenpair** with
an identical explicit complex initial state, supplementary RNG seed 19,
an absolute residual requirement of 1e-8, and the same complete
local magnetization and nearest-neighbor two-body correlation outputs.
ALG=2's NOK/NOM workspace is varied as a separate variable.
The operator cache is held at wk_dim=-1.
For each parameter point one warmup is excluded, followed by
two interleaved measured runs, all in fresh output directories.

**Numerical contract:** all 18 runs per OMP configuration were
successful (36/36 across OMP=1/2), zero nonconvergent cases and
zero contract violations. Each accepted run independently checked
the lowest energy and complete observable arrays against the
SciPy reference, as well as direct residual acceptance. The
table uses medians of the two measurements; maximum RSS is measured
for the ED2 subprocess using GNU time, in KiB.

| L | OMP | Solver (NOK/NOM) | Median H applications | Median wall s | Median RSS KiB |
| ---: | ---: | --- | ---: | ---: | ---: |
| 8 | 1 | ALG=1 | 120 | 0.030616 | 8,454 |
| 8 | 1 | ALG=2 2/16 | 306 | 0.075426 | 10,160 |
| 8 | 1 | ALG=2 4/24 | 475 | 0.123207 | 10,716 |
| 10 | 1 | ALG=1 | 200 | 0.285188 | 22,128 |
| 10 | 1 | ALG=2 2/16 | 340 | 0.547111 | 27,850 |
| 10 | 1 | ALG=2 4/24 | 550 | 0.930786 | 30,800 |
| 8 | 2 | ALG=1 | 120 | 0.017098 | 8,582 |
| 8 | 2 | ALG=2 2/16 | 306 | 0.048065 | 10,238 |
| 8 | 2 | ALG=2 4/24 | 525 | 0.091589 | 10,708 |
| 10 | 2 | ALG=1 | 200 | 0.135527 | 22,272 |
| 10 | 2 | ALG=2 2/16 | 340 | 0.358266 | 28,054 |
| 10 | 2 | ALG=2 4/24 | 550 | 0.703913 | 30,976 |

The H-application count at L=8, ALG=2 4/24 is 475 at OMP=1
but 525 at OMP=2; the comparison must **not claim identical Krylov
trajectories across thread counts**. Achieved residuals and
physical correctness are assessed independently in each run.
Two repeats are insufficient to estimate reliable timing error
bars or strong-scaling efficiency.

The relative solver times are **conditional on this single physical
model and selected restart workspaces**. For the tested
single-ground-state task, ALG=1 used fewer H applications and less
RSS than ALG=2, but no general algorithm superiority follows.
The current ALG=2 implementation can use deterministic supplementary
random directions beyond the common initial vector. This is a fair
matched **target**, not a promise of identical internal subspaces.

## C. What the paper can cite now, and what remains

The literature-model test and all 36 resource smoke trials now pass.
It is reasonable to report the finite-size Hamiltonian agreement
and the independent-residual and observable validation.
**Do not** claim strong parallel scaling, portability, universal
workspace optima, an extrapolated chirality susceptibility or a
CPC-ready performance evaluation from this single CI-hosted study.

For a journal-ready benchmark, execute the same script repeatedly
on a controlled node, pin affinity and BLAS threads, collect source SHA,
compiler and library versions, host CPU/cache and memory, compare
different cache policies and a second spectral workload, and extend
to L=12 where feasible. Retain raw output and report medians,
variability and honest nonconvergence rather than only successes.

## Reproduction

**Literature calculation (small dense):**

    python3 tests/cpc_literature_yi_xyz_dm.py --exe /path/to/QS3ED2 \
      --sizes 4,6,8 --omp 1 --output /new/empty/yi-omp1

**Matched stationary resource calculation (SciPy reference):**

    python3 tests/cpc_stationary_resource_benchmark.py --exe /path/to/QS3ED2 \
      --sizes 8,10 --D 0.12 --omp 1 --residual-atol 1e-8 \
      --repeats 2 --warmups 1 --workspaces 2:16,4:24 \
      --output /new/empty/matched-omp1

CI provenance: https://github.com/ueda1984/ED2/actions/runs/43702189530 .
