# Native CPC physical-response reproduction on generic ED2 main

**Provisional review evidence — 19 November 2026.**
Do not mistake this for a tagged CPC Program Library release or
direct residual certification.

## Exact tested source and tools

- Source head at reproduction: `c40c904cd79c5e348d85ef876076d4cda58dea67`,
  [GitHub Actions proof #44445263717](
  https://github.com/ueda1984/ED2/actions/runs/44445263717).
  The `source/ham2vec.f90` Git blob SHA
  `721e6ee5dffc0d999aaefa6d11eae5fc4038da3e`
  is identical in the current fork `main` and this paper branch.
- GNU Fortran 14, `-O2 -fopenmp`, system OpenBLAS/LAPACK,
  `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`.
- Actual binary SHA256
  `4d772caa1e37b59f731587a2d818850812060a96e7e4e717d1d49e49ab5e3cd9`.
- Regenerator: `paper/tools/cpc_portable_dm_quickstart.py`,
  same source for all (D=-0.0002,0,+0.0002).
- Native input: genuine non-`U(1)` OBC alternating
  spin-1/2 and spin-1 XYZ + symmetric off-diagonal
  Gamma + staggered oriented DM and transverse field.
  Full local Hilbert spaces: 36 states (4 sites) and
  216 states (6 sites), **not** a falsely conserved
  fixed-Sz sector.
- Full raw archive: `cpc-ed2-current-main-native-dm`
  workflow artifact #14492907330, containing input decks,
  eigenvalues, local magnetizations, all nine Cartesian
  correlation channels, stdout and structured `report.json`
  for each strength and both system sizes.

## Results

| Physical quantity | Four sites | Six sites |
|---|---:|---:|
| D=0 E0 | -2.230937250601721 | -3.583863495866078 |
| D=0 staggered chirality $K$ | +0.0011502476661632913 | +0.0010900534820340556 |
| Central finite-difference $\partial_D E_0$ | +0.001150247719294839 | +0.001090053619590137 |
| $|\partial_D E_0-K|$ | 5.3132e-11 | 1.3756e-10 |
| Susceptibility $\chi_K/(L-1)$ | 0.6343111813578161 | 0.5687329207582081 |
| Independent source reference susceptibility | 0.6343112017685585 | 0.5687329253277751 |
| Complete native runs | 3/3 PASS | 3/3 PASS |

The independent 25-case small-spin reference values were prepared
on the earlier physics branch and checked again in the paper CI
using a separately constructed NumPy/SciPy matrix for 15
L4/L6/L8 points. These source and native comparisons are
**complementary**; the native values are not produced by
replaying the independent matrix builder.

### Solver diagnostic caveat

The generic `main` has no native
`ED2_DIAGNOSTICS` direct-residual sidecar.
The portable paper quickstart therefore checks native
thick-restart convergence flags, the final reported
restart indicator, the two native `<Phi|H|Phi>` energy
expectations, and their consistency with the two
eigenvalues. The largest discrepancy among the native
expectation checks was below $1.4\times10^{-14}$.

Those indicators **do not mathematically equal** the direct
residual norm $\|H\psi-E\psi\|_2$.
The exact-source program cannot claim an independent
direct-residual certificate on the basis of this check.
Each structured `report.json` records
`direct_residual_certified=false`.
The `ED2_RNG_SEED` variable, retained for compatibility
with experimental branches, is not proof of a
seed-controlled original-main solver.

The full native cache `H|\psi\rangle` bits and
independent dense-H construction are separately
regressed on the generic mainline using the
`tests/cpc_csr_main_regression.py` suite. This is a
check of operator implementation, not of residual
norms for every native solver eigenpair.

## Repeat locally

```bash
autoreconf -vfi
./configure FC=gfortran-14 FCFLAGS="-O2 -fopenmp"
make -j2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -S paper/tools/cpc_portable_dm_quickstart.py \
  --exe "$PWD/source/QS3ED2" --size 4 --omp 1 \
  --output "$PWD/cpc-demo-l4"
python3 -S paper/tools/cpc_portable_dm_quickstart.py \
  --exe "$PWD/source/QS3ED2" --size 6 --omp 1 \
  --output "$PWD/cpc-demo-l6"
```

Use fresh output directories. Check the `report.json` and
retain the source SHA. These runs establish a
valid finite-size demonstrator, **not** new bulk phase
physics, a speedup factor, or the final publisher archive.
