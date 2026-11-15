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

## CPC development and reproducible worked example

**Publication status:** this is the first proposed ED2 software-program
journal paper. ED2 source was publicly available before this manuscript,
but ED2 has **not** previously appeared as a CPC program paper.
The already-published predecessor is QS³
([QS-Cube/ED](https://github.com/QS-Cube/ED),
*CPC* **277**, 108369, 2022), not this ED2 repository.
A [source-level and capability comparison](docs/development/cpc-qs3-ed-to-ed2-scope-matrix.md)
documents the spin-model generalization and records the
historical dynamical structure-factor functions not
currently demonstrated by the new ED2 executable.
The [pinned original-program XXZ cross-check](docs/development/cpc-qs3-ed-to-ed2-scope-matrix.md)
now reproduces **32 levels across three four-/six-site U(1)
sectors**, including the original near-saturation \(N_\downarrow=1\)
regime. The [frozen 32-level spectrum](docs/development/cpc-qs3-ed-ed2-32level-extended-ci.csv)
shows a maximum original-to-ED2 difference below 3e-15.

The [CPC submission-focused draft](docs/development/cpc-submission-core-draft.md)
and [upstream-versus-development contribution audit](docs/development/cpc-software-provenance-contribution-audit.md)
are **development documents, not a submitted article or official release**.
For a quick physics reproduction after building this branch,
run the [mixed-spin DM/chirality quickstart](docs/development/cpc-portable-dm-quickstart.md):

```bash
python3 -S tests/cpc_portable_dm_quickstart.py \
  --exe "$(pwd)/source/QS3ED2" \
  --size 4 --omp 1 --output "$(pwd)/cpc-portable-example"
```

This example uses only the Python standard library and the compiled ED2
binary. It verifies an open, two-site-periodic mixed-spin XYZ/Gamma/DM
Hamiltonian with **nonconserved total Sz**, its Hellmann–Feynman
vector-chirality response and independent frozen finite-size results.
The source-archive CI validates L=4/6 at OMP=1/2. The exhaustive
independent reference scripts remain optional and require NumPy/SciPy.
No controlled-node speedup, stable installed Fortran API or software DOI
is claimed by this development fork.

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
