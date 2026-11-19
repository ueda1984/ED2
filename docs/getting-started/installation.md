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

- Fortran compiler (tested with `ifort`, `gfortran`)
- BLAS / LAPACK
  - OpenBLAS (preferred)
  - Intel MKL (automatic fallback supported)
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
./configure FC=gfortran
# or
./configure FC=ifort
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
./configure LIBS="-qmkl=parallel"
```

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

- The executable is **not installed into PATH**
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

!!! warning "Heavy regression case"

    `make check-long` runs the historical ex1–ex9
    examples. The final ex9 is a large 1000-site
    three-magnon workload, intended for systems with
    sufficient CPU resources. On an ordinary workstation,
    first run `make check` and the ex1–ex8 subset, then
    schedule ex9 on suitable resources if needed.

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

Floating‑point roundoff may vary across systems, but physical observables should remain consistent within numerical precision.

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

- ifort: enabled by default
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

This installation procedure is designed for portability, reproducibility, and HPC robustness.
