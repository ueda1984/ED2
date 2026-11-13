# Phase 2-B: source-archive numerical validation and release metadata audit

**Status:** validated development work, **not** a released or accepted CPC program.
**Source branch:** `codex/ed2-phase2b-release-audit`.
**Initial parent:** `f8dd63a2ef5eddf9b8a1a0cc5534769ffdda85e3` (Phase 2-A plus
the reviewed MkDocs documentation).

## Archive-only numerical validation

A new [GitHub Actions matrix run](https://github.com/ueda1984/ED2/actions/runs/43392314521)
builds a source distribution, extracts it into an unrelated temporary tree,
then compiles and tests **only the extracted source**.

Profiles: Ubuntu 22.04 with GNU Fortran 11, and Ubuntu 24.04 with
GNU Fortran 13; `-O2 -fopenmp`, OpenBLAS/LAPACK, with the process
default of `OMP_NUM_THREADS=1` and single-threaded BLAS. Several
individual tests explicitly exercise `OMP_NUM_THREADS=2`.

| Test from a fresh source archive | GNU 11 | GNU 13 |
| --- | --- | --- |
| Independent Hamiltonian/basis/observable comparisons | 21/21 | 21/21 |
| LAPACK numerical/interface/error tests | 11 checks passed | 11 checks passed |
| Sparse-field CLI | 70/70 | 70/70 |
| Sparse-field operator/Hermiticity | 10/10 | 10/10 |
| Sparse-field permutation checks | 28/28 | 28/28 |
| Lanczos kernel tests | 82/82 | 82/82 |
| Lanczos CLI tests | 26/26 | 26/26 |
| Phase 1-F diagnostic/replay test cases | 136 executed, test passed | 136 executed, test passed |
| Fixed-seed replay tests | 24 | 24 |
| Phase 1-G residual-policy test cases | 146 executed, test passed | 146 executed, test passed |

The final two suites check that **expected** nonconvergence and rejection
statuses are reported; they are **not** claiming 136 or 146 physical
solves succeeded.

The GNU 11 Phase 1-F tests classified 105 accepted / 31 budget-limited,
versus GNU 13's 106 / 30. The Phase 1-G legacy-policy group similarly
had 115 / 31 versus 116 / 30. Such differences demonstrate why exact
cross-compiler iteration or acceptance-count identity should not be
promised. The test assertions and acceptance contract passed in both
toolchains; no failing numerical test was silently waived.

The earlier Phase 2-A
[archive build and reference comparison](https://github.com/ueda1984/ED2/actions/runs/43350378812)
also passed `make distcheck`, the nine extended physical model
examples (9/9), and 27 reference table comparisons (27/27).
These are distinct but complementary verification layers.

## Release metadata: confirmed inconsistencies

| Source | Recorded value | Assessment |
| --- | --- | --- |
| `configure.ac` | `AC_INIT([QS3-ED2], [1.0.0], ...)` | Current Autotools package label; **not independently a released tag** |
| `CITATION.cff` | `version: "0.1.0"` | Different from Autotools; publication decision required |
| `CITATION.cff` | `date-released: "2026-01-01"` | The corresponding release event has not been verified; do not promote as proven |
| `README.md` | Hiroshi Ueda and Daisuke Yamamoto under maintainers | A role list, not necessarily the complete scientific authorship |
| `CITATION.cff` | Hiroshi Ueda, Daisuke Yamamoto, Tokuro Shimokawa under authors | Preserve existing attribution; confirm citation order and contributions before release |
| `LICENSE` | MIT; Copyright (c) 2023 QS-Cube | Preserve upstream license and notice |
| `source/Makefile.am` | `bin_PROGRAMS = QS3ED2` | Installs an executable; **no installed general-purpose Fortran library API** |

No official Git release tag or software DOI was established during this
audit. Existing numeric versions and `CITATION.cff` authors were
**not changed**. A CPC article DOI, if/when assigned by the journal,
is distinct from any Zenodo software archive DOI.

## Required release decisions

1. **Publication scope:** submit the current ED2 executable as a
   documented scientific program, or implement and validate a separately
   installed public Fortran library API before advertising a library.
2. **Version authority:** approve one exact release version, synchronize
   `AC_INIT`, `CITATION.cff`, documentation and the release tag, and
   verify any `date-released` against a real release event.
3. **Credit:** confirm the intended software authors, contributors,
   historical QS³ attribution and maintainers with the project owners.
4. **Archiving:** confirm upstream integration, review a tagged source
   tarball, select archival DOI policy, and reconcile software and
   prospective article citations separately.
5. **Scientific manuscript:** compare ED2's supported Hamiltonians,
   symmetry sectors, truncation and numerical solvers explicitly with
   existing QS³ and other public ED software. Add a physics-relevant
   reproducible application with machine/resource metadata.

## Claims that remain unsupported

The finite suite does not certify the global lowest state for arbitrary
initial vectors, observable accuracy near degeneracies, Intel `ifort`
or `ifx` builds, GPU or MPI use, or a stable interoperable library ABI.
ALG=2 efficiency and near-degenerate inclusion remain recognized
limitations. The original isolated GNU 13 magnetization anomaly was
not deterministically reconstructed even though the accuracy
counterexamples are documented.

## Publication boundary

This work changes only development branches and validation/documentation.
**Do not** merge `main`, create or announce a release, publish Read the
Docs `latest`, add a DOI, or open a PR without a separate review
and explicit approval.
