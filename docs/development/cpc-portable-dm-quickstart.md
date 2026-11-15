# CPC worked example: mixed-spin DM/chirality response without Python numerical packages

**Status:** source-archive execution verified for **4 site/thread combinations**, each running **three independently checked DM strengths**, in [CI workflow #43711037208](https://github.com/ueda1984/ED2/actions/runs/43711037208). This example is provided for *end users* and uses only the **Python 3 standard library** after ED2 has been built. Unlike the full independent-reference validation suite, it does **not** import NumPy or SciPy. The Fortran executable still requires a compiler to build and BLAS/LAPACK at run time.

This is a compact reproduction of a *finite-size software example*, not a novel chiral phase, a thermodynamic extrapolation, a speedup benchmark or an installed Fortran library interface.

## What it computes

Use a connected, open, two-site-periodic alternating \(S=\tfrac12,1\) chain,

\[
H(D)=H_0+D\sum_{i=1}^{L-1}(-1)^{i-1}
(S_i^xS_{i+1}^y-S_i^yS_{i+1}^x),
\]

where \(H_0\) retains distinct \(J_x,J_y,J_z\), symmetric off-diagonal \(\Gamma_z\), and two-site-periodic transverse and longitudinal fields. **No fixed-\(S^z\) sector** is imposed. The sizes are \(L=4,6\), with Hilbert dimensions \(36\) and \(216\).

The script produces a three-point scan at \(D=-2\times10^{-4},0,+2\times10^{-4}\). From the ED2 *native* correlation file it calculates

\[
K(D)=\sum_{i=1}^{L-1}(-1)^{i-1}
\left[\langle S_i^xS_{i+1}^y\rangle
-\langle S_i^yS_{i+1}^x\rangle\right].
\]

**Important:** the ED2 Cartesian correlation output columns are **xx, yx, zx, xy, yy, zy, xz, yz, zz**, not the textbook row-major order. The program deliberately subtracts `xy − yx`, so a component-order swap fails the Hellmann–Feynman acceptance.

It verifies

\[
\frac{E_0(+h)-E_0(-h)}{2h}\approx K(0),
\quad
\chi_\kappa(0)\approx
-\frac{E_0(+h)-2E_0(0)+E_0(-h)}{h^2(L-1)}
\]

against the independently generated frozen reference rows in
[cpc-stationary-dm-reference.csv](cpc-stationary-dm-reference.csv).
As with any finite-difference approximation, these equations have
nonzero truncation and floating-point errors.

## Reproduce on Linux after building the branch

For a fresh clone of the **development branch** (not the upstream published release):

```bash
git clone --branch codex/ed2-cpc-stationary-dm-scaling \
  https://github.com/ueda1984/ED2.git ed2-cpc
cd ed2-cpc
autoreconf -vfi
./configure FC=gfortran FCFLAGS="-O2 -fopenmp"
make -j2

export OPENBLAS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE

python3 -S tests/cpc_portable_dm_quickstart.py \
  --exe "$(pwd)/source/QS3ED2" \
  --size 4 --omp 1 \
  --output "$(pwd)/cpc-portable-l4"
```

The output directory **must not already exist**; this prevents accidentally reading stale solver results. For a larger example use `--size 6` and a different output directory. To test 2 OpenMP threads, set `--omp 2`. The script pins its finite-difference step to the reviewed \(h=2\times10^{-4}\) reference, and sets a fixed RNG seed and requested residual threshold of \(10^{-8}\).

**Outputs:** each directory contains complete input records (`input.dat`, `spins.dat`, `caps.dat`, `bonds.dat`, `fields.dat`, `pairs.dat`), ED2 low eigenvalues, local magnetizations, nine Cartesian nearest-neighbor correlations, solver diagnostics and captured stdout. The parent `report.json` includes the energies, DM-conjugate chirality, numerical derivative, susceptibility, errors, binary SHA256 and convergence status.

The script returns nonzero when a requested eigenpair fails to converge, residual checks fail, an expected output file is missing/malformed, the Hilbert dimension or output component order is inconsistent, or physics identities/reference values exceed the acceptance tolerances. It records failures rather than silently accepting them.

## Reference values at zero DM coupling

| \(L\) | Full dimension | \(E_0(0)\) | \(K(0)\) | \(\chi_\kappa(0)\) per bond |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 36 | -2.2309372506 | +0.001150247666 | 0.6343112018 |
| 6 | 216 | -3.5838634959 | +0.001090053482 | 0.5687329253 |

The values were computed and independently checked in the
[60-case stationary DM study](cpc-stationary-dm-application.md), not invented by this new quickstart. The quickstart is a **reproducer and acceptance check** using frozen independent reference values; it does not recompute an independent Hilbert-space Hamiltonian, which is why no SciPy dependency is needed.

## Limits

- Only the two reviewed sizes \(L=4\) and 6, the frozen model coefficients and fixed \(h=2\times10^{-4}\) have independent acceptance thresholds in the standalone script. For \(L=8,10,12\) or larger DM scans, use `tests/mixed_spin_dm_chirality_scan.py` with NumPy/SciPy and the original reference dataset.
- A direct residual below \(10^{-8}\) does **not** alone certify global ground-state inclusion in arbitrary near-degenerate models; this example also checks reference energies and requires a separated low gap.
- The script does not assert performance superiority, installed library API stability, MPI support or quantum critical behavior.
- Source hashes, archive/build versions, compiler/flags and external BLAS library should be recorded for any use in a formal program release. `configure.ac` and `CITATION.cff` still contain inconsistent version metadata awaiting maintainer resolution.

The exact test source is [cpc_portable_dm_quickstart.py](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-stationary-dm-scaling/tests/cpc_portable_dm_quickstart.py). The [4-job CI workflow](https://github.com/ueda1984/ED2/actions/workflows/cpc-portable-dm-quickstart.yml) builds the extracted source distribution and uses `python3 -S` to demonstrate the absence of optional Python numerical dependencies.