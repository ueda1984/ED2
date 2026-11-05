# Agent instructions for ED2 (QS3-ED2)

## Mission and scope
This repository is a fork of https://github.com/QS-Cube/ED2. The development objective is a maintainable, reusable Fortran software library and a reproducible software paper targeting *Computer Physics Communications* (CPC). Start by understanding and preserving the scientific behavior of the existing ED2 executable; do not treat a successful build as sufficient numerical validation.

The upstream repository is the historical reference. Work in feature branches and submit focused pull requests to this fork. Do not change `main` directly, force-push, or create release tags without explicit authorization.

## Existing project structure
- `source/`: Fortran implementation; `source/main.f90` is the current executable entry point.
- `source/Makefile.am`: Autotools program target `QS3ED2` and explicit Fortran module dependency order.
- `configure.ac`, `Makefile.am`: Autoconf/Automake configuration.
- `tests/`: regression scripts. Currently `make check` runs example ex1, while `make check-long` runs ex1 through ex9.
- `script_*/`, `input/`, `examples/`: existing model generators, inputs, and reference workflows.
- `docs/`: user and developer documentation.

Check the repository contents rather than relying on this summary if filenames or targets change.

## First task: baseline audit (before refactoring)
1. Record the starting Git commit and the actual compiler/build environment.
2. Inspect module dependencies, global state, inputs/outputs, symmetry handling, matrix-vector multiplication, and eigensolver interfaces.
3. On a clean Linux build environment with GNU Fortran, Autotools, and BLAS/LAPACK available, try:
   ```sh
   autoreconf -vfi
   ./configure FC=gfortran
   make -j2
   make check
   ```
   Run `make check-long` when resources permit. Report any unsupported or failed step explicitly, with relevant logs. Do not claim tests passed unless they were executed.
4. Capture numerical reference results (energies, eigenvalues, observables, sector dimensions where available) with stated tolerances and reproducible inputs. Include edge cases such as zero interactions, small systems, and nontrivial symmetry sectors where feasible.
5. Propose a modularization plan before making invasive API or solver changes.

## Scientific invariants and validation
- Preserve Hamiltonian conventions, site/bond indexing, signs, normalization, basis ordering, symmetry sector definitions, truncation semantics, output interpretation, and numerical accuracy.
- Preserve the original command-line executable and its input format while adding reusable library interfaces, unless an explicit change is approved.
- Use deterministic, small-system checks against dense exact diagonalization or independently verified results whenever practical; distinguish absolute and relative error tolerances.
- Validate both single-thread and OpenMP execution; control nested BLAS parallelism (for example `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` for baselines).
- Avoid broad algorithmic rewrites, numerical optimizations, or silent changes to floating-point behavior during initial library extraction.
- Do not edit or regenerate reference output solely to make failing tests pass.

## Build and library design
- Keep existing Autotools build and `QS3ED2` operational during migration.
- Prefer explicit module ownership and public/private interfaces, well-defined kinds, resource cleanup, error reporting, and separate I/O from computational kernels.
- Minimize mutable global module state where possible; plan for independent solver instances, library reuse, and OpenMP safety.
- Document library API consumers, compiler requirements, ABI/module compatibility constraints, and installation/linking usage. Do not promise a stable binary ABI without design and tests.
- Avoid introducing mandatory external dependencies without a demonstrated need.
- Ensure any new tests are integrated into a reproducible test command and suitable for CI.

## Publication and provenance
- Retain the upstream attribution, contributors, MIT license, and existing citations. Treat software version and DOI updates as a deliberate release step, not an automatic cleanup.
- Document scientifically meaningful capabilities separately from engineering improvements, with honest novelty claims.
- Capture benchmark machine/compiler/BLAS/OpenMP information and input datasets for CPC reproducibility.
- Do not assert acceptance, DOI assignment, or journal suitability as a result of code changes.

## Work reporting
For each feature branch or PR, summarize:
- Purpose and changed files;
- Physical/numerical invariants affected;
- Build and tests actually run, with pass/fail/skipped status;
- Known risks and next validation steps.

Never merge to `main` or submit an upstream PR without an explicit request.
