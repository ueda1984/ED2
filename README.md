# ED2 (QS³-ED2)

**Exact Diagonalization for Quantum Spin Systems**

ED2 is a Fortran code for quantum spin lattice calculations with site-dependent
spins, anisotropic exchange, Dzyaloshinskii–Moriya interactions, fields, and
explicit spatial-symmetry inputs. It provides dense diagonalization and
Lanczos-type solvers, with shared-memory OpenMP parallelism and LP64
BLAS/LAPACK. It has no MPI dependency.

## Documentation and source scope

The official documentation is [ED2 on Read the Docs](https://ed2.readthedocs.io/en/latest/).
Start with the [installation guide](docs/getting-started/installation.md) and
[quickstart](docs/getting-started/quickstart.md), then consult the
[input](docs/user-guide/input-format.md) and [output](docs/user-guide/outputs.md)
references.

### Historical manual and retired setup script

The [2023 Japanese manual](https://github.com/QS-Cube/ED2/blob/main/user_manual_jpn.pdf)
is historical material. Its
field-sign convention and build instructions differ from the current code;
use the current [Read the Docs installation guide](https://ed2.readthedocs.io/en/latest/getting-started/installation/)
and input reference above for current behavior. The legacy `setup.sh` entry
point is retired: it exits with an error and prints the current build commands.
The archived PDF is preserved unchanged.

The general complex-CSR cache, coefficient-certified capacity, padded-cache
release, and two-pass construction were merged into official `main` in
[PR #1](https://github.com/QS-Cube/ED2/pull/1) on 19 November 2026. All four
controls remain **off by default**. See
[performance and memory](docs/user-guide/performance-and-memory.md) for their
prerequisites and allocation accounting. Separate cubic-specific experimental
fork features are not implied by this merge. Record the actual source commit;
`main` and Read the Docs `latest` are moving references, not release identifiers.

## Build with GNU Fortran

A Git checkout requires Git, a C compiler, GNU Make, Autoconf (2.69 or later),
Automake, GNU Fortran with OpenMP, and LP64 BLAS/LAPACK development libraries.
On Debian/Ubuntu:

```bash
sudo apt-get update
sudo apt-get install git build-essential gfortran autoconf automake libopenblas-dev liblapack-dev
git clone https://github.com/QS-Cube/ED2.git
cd ED2
git rev-parse HEAD
autoreconf -vfi
./configure FC=gfortran FCFLAGS="-O2 -fopenmp"
make -j2
```

This builds `source/QS3ED2` and the `script_*/mk_input` helpers. No executable
installation is required to run the examples. Optional `make install` installs
`QS3ED2` under the configured `prefix/bin`; see the installation guide for a
writable local prefix.

The default BLAS search is OpenBLAS, BLAS, then f77blas, with a separate LAPACK
check. The legacy `ifort` option remains supported as a build route; `ifx` with
Intel oneAPI LP64 oneMKL is available via `FC=ifx --with-lapack=mkl`
and is verified separately from GNU. Use conservative
`FCFLAGS="-O2 -qopenmp -fp-model precise"`; see the
[installation guide](docs/getting-started/installation.md) for the exact
commands. ILP64 libraries are incompatible with the declared LP64 interfaces.

The audit used Debian 13/GNU Fortran 14.2 and system BLAS/LAPACK. The repository
CI also exercises GNU Fortran 13 and 14 on Ubuntu 24.04; this is not a blanket
claim of portability to all compilers and platforms.

## Run and inspect an example

Set thread counts before validation:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
make check
(cd script_chain && ./run.sh)
head -n 5 examples/chain/eigenvalues.dat
```

`make check` runs ex1 in a disposable directory, checks nonempty outputs, and
requires a positive state index and finite-number energy syntax on every row
of `eigenvalues.dat`.
The explicit script invocation leaves inputs, logs, energies and observables
under `examples/chain/`. Invoke each script from its own directory because the
input generator reads `fort.30` from the working directory.

In `eigenvalues.dat`, column 1 is the state index and column 2 is the energy.
Direct runs append energy and observable records; the example script clears
only its energy table. Use fresh output directories when preserving separate
calculations. A nonempty file or zero process exit status alone is not a
numerical-accuracy check.

### Additional validation

The historical ex1–ex8 execution suite is:

```bash
make -C tests check-long \
  LONG_TESTS="test_ex1.sh test_ex2.sh test_ex3.sh test_ex4.sh test_ex5.sh test_ex6.sh test_ex7.sh test_ex8.sh"
```

Both wrappers reject malformed energy rows, including `NaN`/`Inf` output;
they do not compare reference energies. For independent small-system
physics checks, install Python 3 and NumPy, then run:

```bash
python3 tests/audit_physics_regression.py \
  --exe "$PWD/source/QS3ED2" \
  --output "$PWD/physics-validation" \
  --with-iterative --field-cache-variants
```

Use a new output directory and retain `summary.json` and the case logs. These
checks cover the tested cases, not all possible inputs. The unqualified
`make check-long` additionally runs ex9, a 1000-site cubic three-magnon
calculation that requires separate resource planning. The archived output
format can be validated quickly with `python3 tests/verify_ex9.py`; this
**does not run ex9**. To check an independently completed HPC run, use
`python3 tests/verify_ex9.py --actual examples/cubic_sp_HB --legacy-ex9-symmetry-projection` and retain the full runtime log.
The historical ladder archive contains non-Hermitian +-/-+ entries because
reflection generators were omitted in the old correlation average; the
explicit projection applies **only** to these ex9 legacy channels, and raw
differences remain reported. Inside an authorized compute-node allocation, use
`sh tests/run_ex9_reproducibility.sh` from a fresh checkout to capture logs,
verify both correlation files, and archive input/output checksums.
See [release qualification](docs/reference/release-qualification.md).

For bounds-checked wavefunction I/O, site-magnetization output for multiple
vectors, and basis-allocation checks, a separate runner needs only Python's
standard library and GNU Fortran:

```bash
python3 tests/audit_runtime_regression.py \
  --fc gfortran --output "$PWD/runtime-validation"
```

It builds isolated test drivers with runtime bounds checking and compares
wavefunction round-trips and representative counts against independent
references. Use a new output directory for this runner too.

## Physical conventions

Input fields enter the Hamiltonian as **`+h.S`**. For a physical Zeeman term
`-B.S`, set input `h=-B`. Bond orientation matters for the DM terms; see the
input reference before reversing site pairs.

An excitation restriction is an exact symmetry sector only when the
Hamiltonian preserves it. Otherwise the code calculates the projected
operator `P H P`. Spatial permutation inputs must preserve the full model,
including fields, anisotropies, and local spins. Users must establish these
conditions and verify solver convergence for their calculation.

## Source tree

| Directory | Contents |
|---|---|
| `source/` | Numerical solver and Fortran modules |
| `input/` | Common input template |
| `script_*/` | Example generators and run scripts |
| `examples/` | Generated example inputs and outputs after running scripts |
| `tests/` | Smoke tests, native cache/ABI regressions, and independent physics checks |
| `docs/` | MkDocs documentation sources |

## Reproducibility and citation

Retain the source commit and any local changes, compiler/version/flags,
BLAS/LAPACK backend, thread settings, complete input files, solver diagnostics,
and reference-check results. Consult
[citation and versioning](docs/reference/citation-versioning.md),
[release qualification](docs/reference/release-qualification.md), and
[CITATION.cff](CITATION.cff). Until a versioned archival DOI is available, cite
the repository with the full commit hash used for the calculation. No GitHub
release is implied by the historical `AC_INIT` package number. Add the CFF
version and actual planned release date only in a reviewed release candidate,
before validating that candidate and creating its tag. A passing metadata
check does not publish a tag, release, or DOI.

## License and maintainers

The repository uses the [MIT License](LICENSE).
Maintainers listed here are Hiroshi Ueda and Daisuke Yamamoto; contributor
citation metadata is provided in `CITATION.cff`.

Report reproducible problems through [GitHub Issues](https://github.com/QS-Cube/ED2/issues).
Include the source revision, build configuration, smallest failing input, and
relevant logs. Contributions can be proposed through pull requests.
