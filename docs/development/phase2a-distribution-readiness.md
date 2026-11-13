# Phase 2-A distribution validation and release-readiness checklist

Status: **Development validation; not a public release.**
Source basis: Phase 1-G `d720a3538d32ceadb598c956fae7564fb5b85460`.
Distribution branch: `codex/ed2-phase2-distribution-validation`.
Combined documentation/distribution branch: `codex/ed2-docs-distribution-integration`.

## Why this work was necessary

The original in-tree build and `make check` did not demonstrate that
a source archive could independently reproduce the executable and its
regressions. The initial GitHub Actions distribution run confirmed that
a tarball could be built, but out-of-tree `make check` failed.
Subsequent independently recorded failures identified:

- Missing example helper scripts, templates and baseline fixtures in the archive.
- Test scripts assuming the source tree, build tree and output tree coincided.
- Read-only source permissions propagated to private test sandboxes by `cp -R`.
- Fortran compiler-generated `.mod` files left after `make distclean`.

Each was repaired in Automake declarations or test harness utilities;
no physics routines, Hamiltonian definitions or reference numerical values
were changed.

## Reproducible validation

Run on Ubuntu 22.04 GitHub Actions using:

- GNU Fortran 11.4.0, GNU Autoconf 2.71, Automake 1.16.5,
  Ubuntu OpenBLAS 0.3.20 and LAPACK 3.10.0.
- `FCFLAGS='-O2 -fopenmp'`; `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`.

The source-distribution checks were:

```sh
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j2
make check
make dist
make distcheck DISTCHECK_CONFIGURE_FLAGS='FC=gfortran FCFLAGS=-O2\ -fopenmp'
```

All of these **passed** in the recorded
[GitHub Actions run](https://github.com/ueda1984/ED2/actions/runs/43350193115).
This includes out-of-tree configure/compile/test, installation,
uninstallation and distclean verification.

The additional archive-only validation has now **passed**:
`make check-long` ran the nine extended model examples (**9/9 passed**),
and `tests/compare_historical.py` reported **27/27 historical reference
comparisons passed** using the unchanged reference data. These comparisons
use the project's documented numerical acceptance criteria, not necessarily
bitwise file identity. The CI log provides reproducible evidence:
[full source-distribution and numerical validation](https://github.com/ueda1984/ED2/actions/runs/43350378812).
The validated source commit was `33fdefff12fab1ae393a46956a7afcc365a25554`.

## Build and documentation integration

The combined branch adds the reviewed documentation changes from
`codex/ed2-docs-cpc-preparation` to the Phase 2-A source.

- The MkDocs strict build has passed on the combined branch.
- A static HTML preview can be obtained from
  [this Actions run](https://github.com/ueda1984/ED2/actions/runs/43350867253)
  under `ed2-docs-html-preview`.
- The public `QS-Cube/ED2` main branch and
  `ed2.readthedocs.io/en/latest/` have **not** been changed.

## Remaining publication gates

- Completed: distribution-archive extended models (9/9) and historical
  comparisons (27/27) passed in the recorded CI run; retain this evidence
  when cutting the release tag.
- Run additional numerical validation from the tarball, including
  the independent Hamiltonian fixtures and the optional solver policy tests.
- Choose the intended public version: `configure.ac` currently uses
  `1.0.0` while `CITATION.cff` says `0.1.0` with an earlier release date.
  These metadata values must not be changed without a release decision.
- Decide whether the CPC paper is an executable software article or
  commits to a separately documented, installed Fortran library API.
  The current source distribution primarily installs the executable.
- Review the final authorship, software citation, DOI strategy and
  publication claims.
- Recheck the user documentation against the exact released code SHA.
- Obtain explicit approval before upstream integration, publishing to
  Read the Docs, creating a PR, merging `main` or making a release tag.

A successful `make distcheck` is an important engineering gate, not a
proof of numerical accuracy for arbitrary physical models.
