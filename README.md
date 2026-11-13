# ED2 (QS³-ED2)
**Exact Diagonalization for Quantum Spin Systems**

ED2 (QS³-ED2) is a **Fortran-based Exact Diagonalization (ED) code** for numerical studies of
**quantum spin lattice models**. The project targets

- reproducible research,
- shared-memory performance via OpenMP,
- minimal external dependencies.

---

## Key features

- **General quantum spin Hamiltonians**
  - Heisenberg / XYZ-type interactions
  - External magnetic fields
- **Controlled Hilbert-space truncation**
  - Restriction by excitation number (e.g. number of spin-down states)
  - Efficient access to dilute or constrained sectors
- **Numerical solvers**
  - Full diagonalization for very small systems
  - Lanczos-type iterative solvers for larger Hilbert spaces
- **Parallelization**
  - Shared-memory parallelism via **OpenMP**
  - No MPI dependency (single-node execution model)
- **Linear algebra backends**
  - BLAS/LAPACK (OpenBLAS, Intel MKL, or vendor implementations)

---

## Documentation

Upstream documentation (for `QS-Cube/ED2`):

https://ed2.readthedocs.io/en/latest/

Development-branch documentation may be newer than the published upstream site.
Verify the repository and Git SHA before using newly introduced solver controls.

Includes:

- Installation instructions
- Quickstart examples
- Input / output specifications
- Reproducible reference examples
- Validation benchmarks
- Algorithmic background

---

## Requirements

- Linux (x86_64)
- Fortran compiler
  - GNU Fortran (`gfortran`; source-archive verification recorded with versions 11 and 13)
  - Intel oneAPI Fortran (`ifort` / `ifx`; compatibility requires separate validation)
- BLAS/LAPACK (OpenBLAS or Intel MKL)
- OpenMP-capable compiler

---

## Build & Installation (Autotools)

ED2 uses GNU Autotools.

```bash
git clone https://github.com/QS-Cube/ED2.git
cd ED2
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j
```

BLAS/LAPACK:

- OpenBLAS is detected automatically when available.
- Configure checks other BLAS/LAPACK libraries when OpenBLAS is unavailable. The automatic Intel MKL fallback applies **only to `ifort`**, not to GNU Fortran or `ifx`.

No PATH dependency is assumed for the ED2 executable.

---

## Testing

Two test levels are provided.

### Lightweight regression test

```bash
make check
```

Runs a minimal smoke test (ex1 only).
This is intended to be fast and suitable for CI.

### Extended regression suite

```bash
make check-long
```

Runs examples ex1–ex9 sequentially.
This performs full physical workflows and may take longer.

`check-long` is intentionally implemented outside Automake’s test harness to improve robustness on HPC systems.

---

## Directory layout

```
source/     Fortran source code
examples/   Reference physical examples
tests/      Regression tests (ex1–ex9)
tools/      Auxiliary scripts
docs/       Documentation sources
```

---

## Reproducibility

For published results, please record:

- QS³-ED2 version or Git commit hash
- Compiler name and version
- BLAS/LAPACK backend
- OpenMP environment (OMP_NUM_THREADS etc.)
- Input files

ED2 is designed so all examples can be reproduced from a clean `git clone`.

---

## Citation

Before DOI assignment:

```
ED2 (QS³-ED2), QS-Cube/ED2, GitHub repository,
version <tag> or commit <hash>.
```

After DOI release, please use the Zenodo DOI.

A `CITATION.cff` file is provided. Its release version, date and author metadata
must be reconciled with the release manifest before the next formal publication;
for a development version, record and cite the full Git SHA.

---

## License

MIT License. See `LICENSE`.

---

## Authors and maintainers

- Hiroshi Ueda
- Daisuke Yamamoto

---

## Contributing

Bug reports and pull requests are welcome.

For scientific usage, please ensure sufficient information is provided to guarantee reproducibility.
