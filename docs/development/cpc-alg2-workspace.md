# CPC: ALG=2 restart workspace sensitivity

**Status:** validated parameter study for the current ED2 solver,
**not** a universal algorithmic speedup claim.
Numerical kernels, defaults and the published upstream are unchanged.

## Purpose and controlled conditions

A prior [matched single-eigenpair comparison](cpc-matched-krylov.md)
showed that the current restarted Hermitian solver (ALG=2) uses more
Hamiltonian applications than the three-term Lanczos solver (ALG=1)
for a representative mixed-spin, non-U(1) chain.

This study varies **only** ALG=2's retained-subspace size `NOK`
and total workspace `NOM`, keeping these fixed:

- The same open alternating `S=1/2,1` non-U(1) spin chain,
  generalized XYZ/DM/Gamma couplings and magnetic field.
- Full Hilbert space, no symmetry restriction or NOD projection.
- **One** lowest eigenpair, `NOE=1`, requested direct residual
  `ED2_RESIDUAL_ATOL=1e-8`.
- Exactly the same explicit complex starting vector (with
  `ED2_INITIAL_VECTOR`) and fixed supplementary random seed.
- GNU Fortran 13, `-O2 -fopenmp`, one OpenMP and BLAS thread.

Every accepted answer must independently match the lowest eigenvalue,
ground-state local magnetization and all Cartesian pair correlations
from a separately constructed SciPy sparse Hamiltonian. A small
Krylov residual alone is **not** accepted as proof of the lowest state.

## 10-site, 7,776-dimensional numerical results

[Validation execution](https://github.com/ueda1984/ED2/actions/runs/43398482506).
All seven runs (one ALG=1 baseline plus six ALG=2 workspaces)
were accepted, with **zero** contract violations.

| Solver and restart configuration | H applications | Wall time (s) | Peak RSS (KiB) | Measured residual |
| --- | ---: | ---: | ---: | ---: |
| ALG=1, one eigenpair | 140 | 0.529 | 22,304 | Within requested 1e-8 |
| ALG=2: `NOK=2, NOM=16` | 289 | 0.828 | 27,848 | 2.93e-12 |
| ALG=2: `NOK=4, NOM=24` | 775 | 1.680 | 30,840 | 4.54e-15 |
| ALG=2: `NOK=6, NOM=32` | 1,254 | 2.798 | 33,780 | 6.23e-15 |
| ALG=2: `NOK=8, NOM=48` | 1,225 | 3.431 | 39,884 | 6.38e-15 |
| ALG=2: `NOK=10, NOM=64` | 1,625 | 5.575 | 45,824 | 8.23e-15 |
| ALG=2: `NOK=12, NOM=96` | 1,261 | 6.124 | 57,760 | 9.75e-15 |

On this model, the `NOK=2, NOM=16` setting used **82% fewer**
Hamiltonian applications than the original `10/64` setting,
and the observed wall time decreased by approximately **85%**.
The measured residual `2.93e-12` remains well inside the
requested `1e-8` limit.

This shows that **workspace/restart policy strongly affects the
current ALG=2 implementation** and makes a speedup claim based on
one arbitrary `NOK/NOM` choice scientifically weak.

## 12-site, 46,656-dimensional follow-up

A second, independent
[GitHub Actions validation](https://github.com/ueda1984/ED2/actions/runs/43398889217)
on the *same physical model* extended the single-pair fixed-start
and requested absolute residual `1e-8` test to 12 sites.
Every workspace passed the independent energy, magnetization and
two-body correlation checks. The five configurations
(one ALG=1 baseline, four ALG=2 workspaces) were **all accepted**;
zero contract violations and no nonconvergent trials were observed.

| Solver configuration | H applications | Wall time (s) | Peak RSS (KiB) | Measured residual |
| --- | ---: | ---: | ---: | ---: |
| ALG=1 | 260 | 6.746 | 119,408 | 5.142e-10 |
| ALG=2: `NOK=2, NOM=16` | 323 | 9.059 | 148,924 | 5.709e-11 |
| ALG=2: `NOK=4, NOM=24` | 825 | 19.585 | 166,396 | 7.983e-15 |
| ALG=2: `NOK=6, NOM=32` | 1,287 | 31.194 | 183,912 | 1.183e-14 |
| ALG=2: `NOK=10, NOM=64` | 2,470 | 82.157 | 254,316 | 1.189e-14 |

For this one-instance 12-site comparison, `2/16` used about
**86.9% fewer Hamiltonian actions** and **89% less observed wall
time** than `10/64`. It used about **41% less peak RSS**.
At the same time, `2/16` remained more expensive than ALG=1
in this single lowest-eigenpair workload (323 versus 260 H actions).
All timings are *single shared-runner observations*, not independent
replicate samples or speedup confidence intervals.

The result strengthens the narrow conclusion that the current
ALG=2 implementation is **highly sensitive to restart workspace
choices**; it does **not** validate `2/16` as a universal default.
In particular, `NOK=2` cannot retain two independent directions
per low-energy cluster when multiple degenerate eigenpairs are wanted.
A separate several-eigenpair and near-degenerate workload matrix
must be passed before recommending an automatic workspace policy.

## Interpretation and limitations

- All settings were tested on the same source and comparable physical
  accuracy: the direct residual test and independent physical
  observables both passed.
- These are **single measurements** on a CI runner, not a full
  performance distribution or error bars. Hamiltonian application
  counts and reported solver statuses are more reproducible than
  the wall-clock times.
- ALG=2 can add supplementary deterministic seeded random directions
  even when its first start is explicitly fixed; the exact internal
  search trajectory is algorithm-dependent.
- A smaller workspace **may fail or slow down** on clustered low
  spectra, requested multiple eigenpairs, nearly degenerate states
  or other coupling patterns. It is not suitable to change a public
  default based on this one case.
- The current implementations have different reorthogonalization,
  reconstruction and restart costs, so `H` counts do not capture
  every source of work.
- Any optimized default must be validated across a **diverse model
  and spectral gap matrix** and compared against the independent
  accuracy/diagnostic regression suite, not merely this example.

## What to do before a CPC performance figure

1. **12-site OMP=1 follow-up passed** (table above). Repeat tuned
   workspace comparisons at OMP=2, and track nonconverged settings
   explicitly on multiple models and spectra.
2. Repeat timings under controlled node conditions (not only CI),
   with RSS and H-application counts, and report medians and spread.
3. Include both `NOE=1` and `NOE>1`, a near-degenerate spectrum,
   and at least one nontrivial symmetry-reduced test.
4. Only after these checks consider adjusting recommended `NOK/NOM`;
   the solver's numerical algorithm and default behavior must remain
   unchanged until a separately reviewed optimization.

## Reproducibility

```sh
python3 tests/mixed_spin_alg2_restart_sweep.py \
  --exe /absolute/path/to/QS3ED2 \
  --output /new/empty/output-directory
```

Requires test-only NumPy and SciPy. All trial-specific inputs,
diagnostic sidecars and the machine-readable `report.json` are
saved under the requested directory. The CI run linked above records
the exact tested source commit.
