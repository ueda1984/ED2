# ED2 (QS³-ED2) — v1.0.0 release notes

**First versioned software release.** These notes accompany the reviewed
`v1.0.0` source, once its tag and GitHub Release are published. Before
publication, cite the exact Git commit instead of treating these notes
as proof of a release.

## Software and scope

ED2 is an OpenMP Fortran exact-diagonalization program for quantum-spin
lattice models, supporting site-dependent spins, anisotropic exchange,
Dzyaloshinskii–Moriya interactions, external fields, user-defined spatial
symmetries, Lanczos and dense Hermitian eigensolvers, and LP64 BLAS/LAPACK.
It provides example input generators, documented inputs/outputs, and
optional complex-CSR/partial-cache Hamiltonian actions.

## Qualification evidence

- Independent small-system spin-operator and dense-H comparisons, complex
  BLAS/LAPACK ABI tests and native runtime/invalid-input regressions.
- GNU Fortran 13/14 with LP64 BLAS/LAPACK; Intel oneAPI ifx 2026.2.0
  with oneMKL (sequential LP64), as recorded in the CI logs.
- Unpacked `make distcheck` and SHA-256 checks for all author-maintained
  release inputs and generated build files.
- A separate full 1000-site, three-magnon, specified-symmetry-sector
  `ex9` reproduction with one energy, 1000 magnetizations and nine
  two-body correlations in each output format. The Ubuntu/Debian full
  runs use the same algorithm, not independent numerical implementations.

See [compiler validation](docs/reference/backend-validation.md) and
[release qualification](docs/reference/release-qualification.md) for exact
source revisions, CI runs, solver residuals, benchmark limits and raw
vs historically projected ladder-correlation comparisons.

## Important limitations

- An occupation restriction describes a true symmetry sector only when
  the chosen Hamiltonian preserves it; otherwise ED2 diagonalizes the
  projected operator `PHP`.
- Spatial symmetry labels must define a consistent one-dimensional
  representation of the generated symmetry group.
- Lanczos convergence is not automatically certified for all parameter
  choices. Unresolved near-breakdown cases may be rejected deliberately;
  record independent eigenvector residual checks for scientific results.
- The legacy 1000-site reference has a known ladder-correlation discrepancy.
  The opt-in comparison projection applies only to the four identified
  historical reference tables, not arbitrary correlators or models.
- No claim of MPI, GPU offload, ILP64 compatibility, comprehensive OS
  portability, or a completed CPC publication is made. Intel ifx full
  1000-site `ex9` was not independently rerun.

## Reproduction and release artifacts

Record the exact source commit, configuration, compiler flags, linked
libraries, hardware, inputs, logs, numerical tolerances and output files.
For a checked release source tree:

```bash
python3 -B tests/check_release_metadata.py --release-tag v1.0.0
autoreconf -vfi
./configure FC=gfortran FCFLAGS="-O2 -fopenmp"
make -j2
make check
make distcheck
make dist
python3 -B tests/release_archive_preflight.py \
  --archive qs3-ed2-1.0.0.tar.gz \
  --expected-commit "$(git rev-parse HEAD)" \
  --release-tag v1.0.0 \
  --json-out release-source-manifest.json
```

The versioned source archive and SHA-256 manifest must match the *final
tagged commit*. GitHub automatic source snapshots are not substitutes
for the separately tested Autotools source archive. See
[CITATION.cff](CITATION.cff) for the approved version/date; use a Zenodo
DOI or CPC article DOI **only after it is actually issued**.

## Historical program materials

The 2023 Japanese PDF manual and archived ex9 reference data are retained
as historical materials; see the current online installation/input
documentation for actual semantics and build commands.
