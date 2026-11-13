# Installation Guide

This page provides a **practical, command‑oriented installation guide** for QS³‑ED2.
Conceptual overview, design philosophy, and citation information are documented on the main documentation home page.

This guide focuses strictly on:

- Building the code
- Verifying the build
- Running examples
- Ensuring reproducibility on HPC systems

---

## 1. Requirements

### Required

- A compatible Fortran compiler. Phase 1-G development validation covered GNU Fortran 13 and 14 on Linux; equivalent validation of Intel compilers has not yet been recorded.
- BLAS / LAPACK
  - OpenBLAS (preferred)
  - Intel MKL (configure attempts an automatic MKL fallback only with Intel `ifort`)
- GNU Make
- POSIX shell

### Optional (recommended for HPC)

- Intel MKL
- Environment modules

---

## 2. Obtain the Source

```bash
git clone https://github.com/QS-Cube/ED2.git
cd ED2
```

If building from a fresh Git checkout:

```bash
autoreconf -vfi
```

---

## 3. Configure

### Basic

```bash
./configure
```

or explicitly specify compiler:

```bash
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
# Other toolchains require their own compatibility validation.
```

### BLAS / LAPACK detection

The configure script attempts, in order:

1. OpenBLAS
2. Generic BLAS / LAPACK
3. Intel MKL (`-qmkl=parallel`)

Manual override examples:

```bash
./configure LIBS="-lopenblas"
```

```bash
./configure FC=ifort --with-lapack=mkl
```

The explicit MKL configuration above is supported only with Intel `ifort` in the current `configure.ac`. It is **not** a tested GNU Fortran configuration.

---

## 4. Build

```bash
make -j
```

This produces:

```
source/QS3ED2
```

Important design choices:

- A normal build leaves the executable in `source/`; optional `make install` installs it under the configured prefix.
- Example scripts locate the binary explicitly
- No recompilation occurs inside example runs

This avoids environment‑dependent execution failures on HPC systems.

---

## 5. Regression Tests

### Lightweight smoke test

```bash
make check
```

Runs Example 1 only.
Intended for:

- CI
- quick validation

### Extended regression suite

```bash
make check-long
```

Runs Examples 1–9 sequentially.

Recommended for:

- new clusters
- compiler changes
- BLAS / MKL changes

---

## 6. Running Examples Manually

Example:

```bash
cd script_chain
./run.sh
```

All example scripts:

- use the already built executable
- never recompile sources
- never rely on PATH
- accept explicit override via `ED2_EXE`

Example:

```bash
ED2_EXE=../source/QS3ED2 ./run.sh
```

Helper programs (`mk_input`) are built at `make` time and executed directly by scripts.

---

## 7. Thread Control (IMPORTANT)

QS³‑ED2 uses:

- OpenMP in the main solver
- threaded BLAS / LAPACK

To avoid nested parallelism:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Users may tune these for benchmarking.

---

## 8. Reproducibility Checklist

For numerical reproducibility, record:

- compiler name and version
- BLAS / LAPACK implementation and version
- OpenMP thread count
- BLAS thread count
- success of `make check` and `make check-long`

Floating‑point roundoff may vary across systems, but small residuals alone do not guarantee accurate individual eigenvectors or
observables in near-degenerate problems. See
[Solver accuracy and reproducibility](../user-guide/solver-accuracy.md).

---

## 9. Cleaning

```bash
make clean
```

Remove configuration:

```bash
make distclean
```

---

## 10. Installation (Optional)

Installation is not required for normal usage.

If desired:

```bash
make install
```

Prefix can be set via:

```bash
./configure --prefix=/path/to/install
```

---

## 11. Troubleshooting

### BLAS not detected

```bash
./configure LIBS="-lopenblas"
```

or

```bash
./configure LIBS="-qmkl=parallel"
```

### OpenMP

Ensure compiler supports OpenMP:

- ifort: use the OpenMP settings appropriate to the selected toolchain
- gfortran: may require `-fopenmp`

### Test failures

```bash
make distclean
autoreconf -vfi
./configure
make -j
make check-long
```

---

These instructions are a starting point, not a guarantee for all platforms.
Source distribution and `make distcheck` are under separate audit before release.
