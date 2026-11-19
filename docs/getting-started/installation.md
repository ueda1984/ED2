# Installation Guide

Build ED2 from a clean source checkout, then validate the executable on small
inputs before using a larger model. The commands below use GNU Fortran and
LP64 BLAS/LAPACK on Linux.

## 1. Requirements

A Git checkout needs:

- Git and a POSIX shell;
- a C compiler, GNU Make, GNU Autoconf (2.69 or later), and Automake;
- GNU Fortran with OpenMP support;
- development libraries for LP64 BLAS and LAPACK.

The C compiler is required by `configure` even though the numerical solver is
Fortran. On Debian/Ubuntu, the following packages supply these prerequisites:

```bash
sudo apt-get update
sudo apt-get install git build-essential gfortran autoconf automake libopenblas-dev liblapack-dev
```

On a managed cluster, use the corresponding compiler and BLAS/LAPACK modules
provided by the administrator. Record the selected module versions.

The 19 November 2026 audit built a clean checkout on Debian 13 with GNU
Fortran 14.2 and the system LP64 BLAS/LAPACK libraries. The existing CI also
exercises GNU Fortran 13 and 14 on Ubuntu 24.04. This does not establish that
every compiler, operating system, or BLAS variant is supported.

The configure script supports GNU `gfortran`, classic Intel `ifort`, and
Intel oneAPI `ifx`. Select `FC=gfortran` for the GNU recipe below or `FC=ifx`
for the Intel recipe. If `FC` is not set, automatic compiler discovery tries
`ifx`, then `ifort`, then `gfortran`; an available Intel compiler can therefore
be selected automatically. Intel runtime/compiler verification is tracked
separately from the GNU 13/14 baseline.
The explicit Fortran interfaces require **LP64** integer arguments;
ILP64 BLAS/LAPACK is not supported by these interfaces.

## 2. Obtain and build the source

```bash
git clone https://github.com/QS-Cube/ED2.git
cd ED2
git rev-parse HEAD
autoreconf -vfi
./configure FC=gfortran FCFLAGS="-O2 -fopenmp"
make -j2
```

Record the printed commit hash. `main` can change after the checkout was
created. These commands build `source/QS3ED2` and the `mk_input` helper in each
`script_*` directory. The helpers are needed to generate the supplied examples.
Running `make` alone does not install the executable into a directory on `PATH`.

The explicit `-O2 -fopenmp` flags are the reproducibility-oriented starting
point used by the strict cache-equivalence tests. Omitting `FCFLAGS` with a
recognized GNU compiler selects the historical `-Ofast -fopenmp` defaults.
`-Ofast` allows floating-point reassociation, so bitwise agreement between
cache layouts is not guaranteed. See [Performance and memory](../user-guide/performance-and-memory.md).

### BLAS/LAPACK selection

For the default `--with-lapack=lapack` configuration, `configure` searches for
BLAS in this order: `-lopenblas`, `-lblas`, then `-lf77blas`. It separately
checks `-llapack`. If the BLAS checks fail with an Intel compiler, it attempts
the MKL fallback `-qmkl=parallel` for `ifort` or `-qmkl=sequential` for `ifx`.
If the separate LAPACK check fails while either MKL option is active, the
configuration assumes MKL supplies LAPACK; the subsequent link and runtime
checks must still succeed.

The explicit `--with-lapack=mkl` path accepts `FC=ifort` (historical
parallel MKL) or `FC=ifx` (sequential LP64 MKL).
Installing a library or setting `LIBS` does not prove that its ABI or all
required symbols match. Check `config.log`, the final link command, and the
runtime tests. For libraries outside standard search paths, supply the
appropriate `LDFLAGS` and module environment and retain the resulting log.


### Intel oneAPI ifx and LP64 oneMKL (optional)

Install the oneAPI Fortran compiler and oneMKL development package from
Intel's authenticated APT repository, or load matching oneAPI modules on
your cluster. Use a separate fresh checkout for this build, at its repository
root. If reusing an existing build tree, run `make distclean` before changing
compilers. The following Bash commands use Intel's standard installation
location; with cluster modules, load the matching environment instead of
the `source` command:

```bash
source /opt/intel/oneapi/setvars.sh
ifx --version
echo "$MKLROOT"
autoreconf -vfi
./configure FC=ifx FCFLAGS="-O2 -qopenmp -fp-model precise" --with-lapack=mkl
make -j2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
make check
```

The explicit `--with-lapack=mkl` route selects Intel `-qmkl=sequential`
for `ifx`, which uses the **LP64** interface. Start correctness checks with
the single-thread settings above; later exercise 4 OpenMP threads independently.
Review the final link command and `ldd source/QS3ED2` to confirm actual oneMKL linkage.
Do not infer that a successful build alone proves physical correctness.

The repository's
[oneAPI CI workflow](https://github.com/QS-Cube/ED2/actions/workflows/intel-oneapi-ifx-mkl.yml)
runs automatically for its configured source/test/build-path changes on pushes
and pull requests, and also supports manual dispatch. It tests the Intel
backend with independent spin-matrix cases, typed complex BLAS/LAPACK and cache
matrix-action regressions. It does not replace the full 1000-site ex9 benchmark
or imply ILP64 support.


## 3. Set threads before running tests

Use one OpenMP thread and one BLAS thread for the first validation:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

For later timing experiments, change `OMP_NUM_THREADS` deliberately and keep
the BLAS thread count under control. Record both settings; nested threading
can change runtime substantially.

## 4. Run the smoke and example checks

From the repository root:

```bash
make check
```

Expected harness output includes `PASS: test_ex1.sh`. This runs the chain
example in a temporary directory, checks that the requested output files are
nonempty, and validates every nonblank row of `eigenvalues.dat`: exactly two
columns, a positive integer state index, and finite-number energy syntax.
Malformed rows, including `NaN`, `Inf`, and overflow asterisks, fail this
check. The temporary results are removed afterward. Despite the helper
filename `test_common_strict.sh`, this is **not a comparison against reference
energies or a proof of solver convergence**.

Run the CPU-feasible historical examples separately:

```bash
make -C tests check-long \
  LONG_TESTS="test_ex1.sh test_ex2.sh test_ex3.sh test_ex4.sh test_ex5.sh test_ex6.sh test_ex7.sh test_ex8.sh"
```

These scripts populate the corresponding `examples/` directories. The
long-test wrapper performs the same energy-row syntax checks and reports
`E0` from the second column of the first row. Column 1 is the state index;
column 2 is the energy. These execution and format checks do not compare
physical reference values. Inspect the actual energies, solver diagnostics,
and applicable independent references.

!!! warning "The complete suite includes an HPC workload"

    `make check-long` includes ex9 (`cubic_sp_HB`): a 1000-site cubic model
    in a three-magnon sector. Treat it as a separate resource-dependent
    calculation, not an installation smoke test. Read the
    [ex9 page](../examples/cubic_sp_HB.md) before scheduling it.

### Independent small-system physics checks

The additional audit regression builds spin matrices and tensor-product
Hamiltonians independently with NumPy, then compares ED2 energies,
wavefunctions/residuals, and observables on small models. Python and NumPy are
only needed for this additional test, not for building the solver:

```bash
sudo apt-get install python3 python3-numpy
python3 tests/audit_physics_regression.py \
  --exe "$PWD/source/QS3ED2" \
  --output "$PWD/physics-validation" \
  --with-iterative --field-cache-variants
```

The output directory must not already exist; choose a new name on a later run.
Read `physics-validation/summary.json` and preserve the per-case logs. Passing
these cases validates their tested models and settings; it does not establish
convergence or symmetry validity for every user input.

### Bounds-checked I/O and basis checks

A separate regression uses only Python's standard library and GNU Fortran:

```bash
python3 tests/audit_runtime_regression.py \
  --fc gfortran --output "$PWD/runtime-validation"
```

It builds isolated test drivers with `-fcheck=all`, checks wavefunction
round-trips with three and twelve saved vectors, verifies that two supplied
states produce independent site-magnetization results with one/four threads,
and compares representative counts against an analytical symmetry count with
one/four threads and redundant generators. This comprises eight cases. Use an output directory that does not
already exist, and preserve its `summary.json` and logs.

## 5. Run a supplied example manually

Keep the parent shell at the repository root:

```bash
(cd script_chain && ./run.sh)
head -n 5 examples/chain/eigenvalues.dat
```

The subshell changes directory because `mk_input` reads `fort.30` from its
working directory. Do not replace this with `./script_chain/run.sh` from an
arbitrary directory. The script uses the already compiled solver and helper;
it does not recompile them.

Generated files are placed in `examples/chain/`, including `input.dat`,
interaction lists, `output.dat`, `eigenvalues.dat`, and requested observable
files. The script uses `source/QS3ED2` by default. To select another executable,
pass an absolute path through `ED2_EXE` when invoking the script.

For direct execution of a prepared input, the path from that example directory
is two levels back to the source executable:

```bash
(cd examples/chain && ../../source/QS3ED2 < input.dat > output-rerun.dat)
```

Direct solver runs append energy and observable tables. The example script
clears `eigenvalues.dat` but does not clear every observable file. Archive
previous results and use a fresh output directory when independent runs are
required; otherwise records from multiple runs can be concatenated. See
[Outputs](../user-guide/outputs.md).

## 6. Optional installation

Running from the build tree is sufficient for the supplied examples. If an
installed executable is wanted, select a writable prefix before building:

```bash
./configure FC=gfortran FCFLAGS="-O2 -fopenmp" --prefix="$PWD/local-install"
make -j2
make install
```

This installs `local-install/bin/QS3ED2`. It does not automatically add that
directory to `PATH`. The example generators are `noinst_PROGRAMS`, so retain
the source tree and built helpers when using the example scripts. Choose a
new installation prefix if an earlier installation must be preserved.

## 7. Troubleshooting and rebuilding

| Symptom | Check |
|---|---|
| `autoreconf: command not found` | Install Autoconf and Automake before the first build. |
| No acceptable C compiler | Install the C toolchain; `configure` requires it. |
| BLAS or LAPACK not found | Install the development libraries or load the correct modules; inspect `config.log`. |
| Missing `source/QS3ED2` or `mk_input` | Complete `make -j2` at the repository root. |
| Missing `fort.30` | Run the script from its own `script_*` directory as shown above. |
| Missing `examples/chain` after `make check` | The smoke test uses a disposable directory; run the example script explicitly. |
| Outputs exist but energies look wrong | Check model signs, projection/symmetries, solver diagnostics, and the independent small-system tests. |

Use `make clean` to remove compiled products. Before switching compiler,
compiler flags, or library ABI, use `make distclean`, then repeat the configure
and build steps. Preserve input files and results separately; a clean build is
not a substitute for checking numerical convergence.
