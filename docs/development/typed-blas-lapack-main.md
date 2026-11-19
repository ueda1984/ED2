# Mainline ED2 explicit BLAS/LAPACK interfaces — compiler portability

**Scope:** General ED2 implementation, independent of cubic point-group
optimizations and of the optional complex CSR cache. This change is not
a solver rewrite. It makes the intended Fortran LP64 ABI explicit and
ensures invalid workspace queries are rejected rather than used as array
lengths.

## Changed internals

- `source/blas_lapack_interfaces.f90` supplies explicit interfaces for
  the ED2 BLAS/LAPACK routines used by its classic Lanczos and
  thick-restarted Lanczos/dense diagonalization. It intentionally uses
  **Fortran calling conventions** (no `bind(C)`), especially for
  the complex-return `ZDOTC` function.
- `source/eigen_solver.f90` imports the typed `ZHEEVR`, `DSTEVX`
  and BLAS routines. The full-Hermitian `ZHEEVR` arguments `VL` and
  `VU` are now `real(8)` rather than incorrectly declared integers,
  and the three workspace-query buffers are rank-one arrays
  with exactly one element. Both iterative and dense
  `ZHEEVR` paths use bounded, finite workspace-length checks.
- `source/lanczos.f90` imports explicit BLAS/LAPACK
  function/subroutine interfaces instead of hand-declaring
  external complex function results.
- `source/Makefile.am` **and the checked-in generated
  `source/Makefile.in`** enumerate the module as a first
  Fortran source file and declare the cross-module prerequisites
  for parallel compilation.

No input parameters, Hamiltonian construction, wavefunction format,
algorithm selection, numerical convergence thresholds, default
feature flags or source-level user-facing API are changed.

## Existing base issue and resolved behavior

The inherited `eigen_solver.f90` had multiple implicit calls to
`ZHEEVR` using both scalar workspace-query temporaries and
array workspaces. Because several calls were to the same external
symbol with incompatible assumed argument ranks, GNU13 and GNU14
rejected compilation without `-fallow-argument-mismatch`.

The modernization declares all calls through one consistent,
typed interface, fixes the incorrectly typed `VL`/`VU`,
and checks the returned workspace size before allocation.

**Verification must build with `-O2 -fopenmp` and WITHOUT
`-fallow-argument-mismatch`.** This is part of the CI contract.
Fortran LP64 `integer` and `real64` are required; ILP64 BLAS
backends are not claimed to be supported.

## Native independent checks

The test-only `tests/cpc_blas_lapack_contract_driver.f90`
uses the actual ED2 `fulldia_routines_with_momentum` with
a genuinely complex 2x2 Hermitian matrix whose exact
eigenvalues are `3 +/- sqrt(3)`. It separately checks
eigenvalue ordering, complex eigenvector residuals,
normalization and orthogonality for a nontrivial 3x3
Hermitian matrix. The code then confirms the BLAS
`ZDOTC` complex return ABI by comparison with explicit
`sum(conjg(x)*y)`, `DZNRM2` by a pure Fortran norm,
and `ZCOPY` by exact copied coefficients.

The existing `.github/workflows/main-generic-csr.yml`
checks this driver in both GNU Fortran 13 and 14, as
well as native main-source `make check`, historical
ex1–ex8 and full bitwise complex `H*x` comparisons
on physically real/nonreal DM mixed-spin models.

Independent GNU Fortran 14.2 local smoke results
before the main-branch change:
- Built after `autoreconf -vfi`, `./configure FC=gfortran
  FCFLAGS='-O2 -fopenmp'` and `make -j2`, no compatibility
  compiler flag — PASS.
- `make check` and historical `ex1`–`ex8` — PASS.
- Native `ZHEEVR` analytic Hermitian and complex BLAS
  ABI checks — PASS.
- Original generalized CSR full `H*x` proof across
  N4/N6/N8/N10, mixed spin S1/2,S1, both exactly real and
  genuinely complex DM, OpenMP 1/4 — **16/16 bitwise
  input/output pairs**, 8 dense-H crosschecks and 3 negative
  configuration cases — PASS.

The original N1000 cubic workload `ex9` remains an
HPC-sized test and is not run as part of the quick dual
compiler CI. PRESTO exclusive-node performance
measurements are not claimed. No new performance result
is asserted by this interface modernization.

The separate cubic speedup and CPC draft PRs
remain untouched.
