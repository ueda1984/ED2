# Legacy QS³ versus ED2 at large site counts: experiment

The historical QS-Cube/ED source is pinned at
`a3b09934bba5f966bfb564718528de8778c2b46c`.
The modern ED2 version is the GitHub Actions checkout SHA.
Source is unmodified; both use fixed Ndown, identical XXZ
couplings, the same momentum, and native symmetry-adapted bases.

Default CI (seven source-pinned regressions): 1D periodic rings
L=12,Ndown=2,k=0; L=24,Ndown=3,k=0; L=48,Ndown=2,k=1;
L=96,Ndown=2,k=0; cubic 3×3×3 and 4×4×4 Ndown=2 at Gamma;
cubic 10×10×10 Ndown=1 at R. Shared-runner wall time and RSS
are archived but are **not** an HPC/node-controlled performance claim.

A separate 10×10×10, 1000-site, Ndown=2, Gamma-point stress
with 330-second per-executable timeout is explicitly **non-gating**.
Failure, timeout and asymmetry between executables must be reported
honestly; neither completion nor three-magnon continuity is assumed.

The independent sector dimension is computed by character-weighted
Burnside counting. A translation of order t has C(L/t,Ndown/t)
fixed bitstrings if t divides Ndown, or zero otherwise.
At 10×10×10 Gamma, the expected n=1,2,3 sector dimensions
are 1, 503 and 166167 respectively. The last number is a
**combinatorial prediction, not a validated ED2 run**.
For n=1 the energy is checked against the analytic lattice dispersion.

Original QS³ and modern ED2 use independent random Lanczos
initialization. Energy agreement does not imply matched
iteration histories; observed times are exploratory and not speedups.
A QuSpin comparison must use the same Hamiltonian, conservation
sector, momentum, hardware and accuracy criteria.


## First verified results (2026-11-16)

Source-pinned experiment:
[run 43778530930](https://github.com/ueda1984/ED2/actions/runs/43778530930),
full script and archived inputs in the Actions artifact. The seven
regular test cases passed all momentum-sector size and energy checks
between **two actual Fortran executables**.

| Lattice | Ndown | Momentum | Independent dimension | Original and ED2 agree |
| --- | ---: | --- | ---: | --- |
| 12-site ring | 2 | 0 | 6 | yes |
| 24-site ring | 3 | 0 | 85 | yes |
| 48-site ring | 2 | 2pi/48 | 23 | yes |
| 96-site ring | 2 | 0 | 48 | yes |
| 3x3x3 cube | 2 | Gamma | 13 | yes |
| 4x4x4 cube | 2 | Gamma | 35 | yes |
| 10x10x10 cube | 1 | R=(pi,pi,pi) | 1 | yes, analytic E=744 |

For these small momentum-reduced sectors, we use ALG=3 in both codes.
The original unmodified ALG=1 sometimes segfaults when the number of
Lanczos steps reaches the reduced dimension without satisfying
the convergence exit condition (see next section).
We do **not** silently equate the small-sector dense solver checks
with Lanczos performance. The original and ED2 model inputs and
Fortran compiler provenance are separately archived.

In a distinct **10x10x10, 1000-site, Ndown=2, Gamma** run both
unmodified QS3 and ED2 executed **ALG=1** successfully.
The full fixed-Ndown dimension is 499500; Burnside gives 503
momentum representatives. Energies in the run were

- original QS3: E=738.006774272212,
- ED2: E=738.006774272203 (observed difference approximately 1e-11).

This is corroborated by a third, physically equivalent but
algorithmically independent construction using *relative magnon
displacements modulo inversion*, not the QS3/ED2 basis routines.
Its 503x503 sparse real-symmetric Hamiltonian has 3497 nonzeros,
ground energy **738.0067742722034**, first gap
**0.37531313226088514**, and residual
**9.42e-13**. Verified independently in
[run 43778483902](https://github.com/ueda1984/ED2/actions/runs/43778483902).
This validates both the symmetry basis size and a nontrivial
two-magnon interaction energy at 1000 sites.

**QuSpin 1.0.1 is not absent from this physics regime.**
Using `spin_basis_general` with `Nup=L-Ndown`, explicit
translation maps, `pauli=0` and matching XXZ operators,
all **7/7** small-sector benchmark cases reproduced
the original QS3 and ED2 energies. Including a 1000-site
one-magnon R-sector calculation, the maximum absolute
QuSpin-versus-Fortran energy discrepancy was below 6e-14.
This explicitly rejects any claim that QS3 is the only
software supporting very large site counts or momentum sectors.
A fair *numerical performance* comparison needs CPU-pinned
exclusive runs, separate basis-construction timing and matched
solver accuracy; these shared-runner wall times cannot be used
as journal speedup evidence.

### Historical Lanczos boundary (not an ED2 regression)

The first smoke attempt at [run 43778009408](https://github.com/ueda1984/ED2/actions/runs/43778009408)
triggered three original QS3 ALG=1 SIGSEGV failures on
small reduced sectors. Source audit shows the historical loop
sets `mitr=min(maxitr,dim)` and may finish at `mitr` without
the ordinary convergence exit. An eigenvector coefficient
array `z` can then be deallocated before it is accessed.
The reproducible workaround for the scientific overlap is
dense diagonalization (ALG=3) in both codes for reduced
dimensions <=100; we have **not patched** the original
source. This is a software-correctness observation, not a
new physical phase or a generic Lanczos convergence theorem.

### Three-way 216-site / 1000-site stress results

**Verified [run 43778708860](https://github.com/ueda1984/ED2/actions/runs/43778708860)**:
both source-pinned Fortran executables and independent QuSpin 1.0.1
accepted all seven smaller momentum sectors, then two more
conventional-Lanczos stress cases. The stress results are:

| Lattice | Ndown | Gamma dimension | Original QS3 E0 | ED2 E0 | QuSpin E0 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6x6x6 | 3 | 7,677 | 148.5648907577822 | 148.5648907577438 | 148.5648907577122 |
| 10x10x10 | 2 | 503 | 738.0067742722275 | 738.0067742723085 | 738.0067742722047 |

The 216-site three-magnon unsymmetrized Hilbert space contains
**1,656,360** states, and the maximum tested old-to-QuSpin
energy discrepancy is **7.0e-11**. The 1000-site two-magnon
sector is independently validated by the relative-displacement
Hamiltonian described above, with maximum three-program
discrepancy in this stress run **1.1e-10**.
All these are *finite-size, same-Hamiltonian checks*, not
evidence that the two software stacks have identical runtime costs.

The exploratory metrics from one **shared CI host** for the
stress cases are listed only to aid future controlled benchmarking;
they cannot support speedup claims:

| Case | Engine | Basis and Hamiltonian construction | Solver plus setup / execution | Peak process RSS |
| --- | --- | --- | --- | --- |
| 216, 3 magnons | Original QS3 | included in total | 40.33 s total | 7,204 KiB |
| 216, 3 magnons | ED2 | included in total | 43.82 s total | 8,044 KiB |
| 216, 3 magnons | QuSpin | 5.83 s basis + 28.28 s operator | 0.067 s eigensolver | 232,600 KiB |
| 1000, 2 magnons | Original QS3 | included in total | 4.22 s total | 6,724 KiB |
| 1000, 2 magnons | ED2 | included in total | 3.78 s total | 6,868 KiB |
| 1000, 2 magnons | QuSpin | 2.90 s basis + 9.37 s operator | 0.003 s eigensolver | 227,816 KiB |

The underlying measurements are **not methodologically symmetric**:
the Fortran timings include native setup and eigenvalue work, while
QuSpin times are segmented Python calls and exclude import/startup
overheads. QuSpin's RSS includes Python and its imported scientific
stack, while the Fortran process has none of those costs. All cases
have unpinned CPU affinity, no exclusive-node attestation, one
measurement per engine, and unmatched Lanczos initial states.
A publication-grade speed ratio cannot be extracted from this table.

### Pending

- **1000-site three-magnon** QS3-to-ED2 numerical reproduction
  is not verified; the 166,167 Gamma-space dimension is only
  independently predicted by symmetry counting.
- Repeated pinned-node benchmark with controlled environment,
  explicit BLAS/OpenMP affinity, numerical tolerances and
  setup/solver separation is not complete.
- No software release/tag, official CPC Program Library archive
  or assigned DOI is implied by these CI validations.

## Independent 1000-site three-magnon matrix construction

**Verified [CI #43779705931](https://github.com/ueda1984/ED2/actions/runs/43779705931)**.
The new
[`tests/cpc_independent_1000site_three_magnon_basis.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-stationary-dm-scaling/tests/cpc_independent_1000site_three_magnon_basis.py)
does not run ED2 and does not compute eigenvalues. It constructs the
periodic 10x10x10, spin-1/2, isotropic nearest-neighbor Heisenberg
Hamiltonian at Ndown=3 and Gamma **directly from translation-equivalent
relative coordinates**. Its results are independently checked:

- Full fixed-Ndown space: **C(1000,3) = 166,167,000**.
- Origin-anchored distinct relative pairs: **C(999,2) = 498,501**.
- Canonical translation representatives: **166,167**, exactly one-third
  of the anchored pairs, with no order-three translation stabilizers.
- Resulting 166,167-by-166,167 real-symmetric sparse Hamiltonian:
  **3,151,137 nonzero entries**; antisymmetric-part error **0**.
- A compressed sparse matrix (.npz) and machine-readable report are
  retained as GitHub Actions artifacts. One shared-runner matrix
  construction took about 9.35 seconds, **not a controlled benchmark**.

This is a constructive validation of the `Ndown=3`, 1000-site
symmetry-reduced Hamiltonian *representation*. The eigensolver's
convergence on this large matrix was **not established** in the tested
workflow: an exploratory local `eigsh(k=2)` attempt did not finish
within its execution cutoff. No source-pinned ED2 1000-site
three-magnon energy has been reported or inferred. The matrix may
serve as a future independent reference once a controlled numerical
calculation is completed.

### Large-numerical-space engineering limit

The current Fortran basis implementation uses default-width integers
for cumulative bounded-combination counts. Even though
C(1000,3)=166,167,000 remains below signed 32-bit maximum,
C(1000,4)=41,417,124,750 exceeds it substantially.
A claimed generalization to four or more dilute magnons
at 1000 sites would therefore require an explicit checked
64-bit counting/indexing strategy and independent validation.
