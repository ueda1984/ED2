# Citation and Versioning

Record the exact source revision used for a calculation. A branch name such
as `main`, the Read the Docs label `latest`, and a successful documentation
build each identify different things; none alone fixes the executable used
to produce a scientific result.

## Current citation

Until an archival software DOI or a published ED2 article is available,
cite the repository and the complete Git commit hash:

> ED2 (QS³-ED2): Exact Diagonalization for Quantum Spin Systems,
> QS-Cube/ED2, GitHub repository, commit `<full-commit-hash>`.

Obtain the revision from the source checkout used to build the executable:

```bash
git rev-parse HEAD
git status --short
```

If there are local changes, also archive the patch and the modified source
files; a commit hash alone does not describe a modified executable. Keep
the build log and compiler flags with the calculation results.

The repository's [CITATION.cff](https://github.com/QS-Cube/ED2/blob/main/CITATION.cff)
contains the software title and contributor metadata for citation tools.
The approved first software version is **`v1.0.0`**. The release-source
`CITATION.cff` records `version: 1.0.0` and
`date-released: 2026-11-20`, matching the planned initial publication
date. If publication is rescheduled, update that date and rerun the complete
release checks **before** tagging. These metadata fields alone do not
establish the existence of a release: first verify that an actual
`v1.0.0` tag and GitHub Release exist and point to the qualified source.
Until then, cite the full Git commit rather than a release not yet issued.

The older CFF historically contained an unissued `0.1.0` version and an
unverified `2026-01-01` release date. The Autotools `AC_INIT` package
identifier already reads `1.0.0` and now matches the approved software
version; it is not by itself proof of publication.

## Tags, releases, and documentation versions

At the 19 November 2026 documentation audit, the official repository
published no Git tags or GitHub releases. Check the current
[tags](https://github.com/QS-Cube/ED2/tags) and
[releases](https://github.com/QS-Cube/ED2/releases) before citing a release.
The `latest` documentation follows the configured development branch;
it is not an immutable software version.

For a reproducible documentation/source comparison, record:

- the source commit used to build ED2;
- the Read the Docs build number and its checked-out commit;
- the documentation URL and the date consulted;
- any local source or documentation changes.

The [Read the Docs build history](https://app.readthedocs.org/projects/ed2/builds/)
provides the build's commit and command output. A successful build proves
that documentation was generated; numerical validation requires the
tests and reference calculations described in the user guide.

## Intended release policy

The intended naming convention for future releases is
`vMAJOR.MINOR.PATCH`:

| Component | Intended meaning |
|---|---|
| MAJOR | Changes incompatible with earlier documented input or behavior |
| MINOR | New features compatible with the previous public interface |
| PATCH | Compatible bug fixes, performance changes, or documentation fixes |

The first software version `v1.0.0` follows this scheme. Its tag, when
published, must identify the frozen, tested source commit with corresponding
documentation, verification evidence, and archive. A merge to `main`
alone does not complete that process.

For the full release gate, including the expensive cubic benchmark, see
[Release qualification and ex9 validation](release-qualification.md).
Check metadata using `python3 tests/check_release_metadata.py`. Development
may omit both fields; a version-only candidate returns `release_candidate`.
The final versioned source includes both `version` and the actual
`date-released`. Run `python3 tests/check_release_metadata.py --release-tag
v1.0.0` against the frozen final commit, and rerun all GNU/Intel/archive
checks on that *same* source tree. The tag-specific check rejects missing
release dates or mismatched versions. A passing check proves metadata
consistency, not that publication has actually occurred.

## Archival software DOI and article citation

Zenodo integration is a planned way to archive versioned source releases.
Once an archive exists, cite the DOI for the specific software version
used and retain the source commit in the reproducibility record.

A software-archive DOI and a journal-article DOI identify different works.
If an ED2 article is published, use its actual bibliographic metadata in
addition to the versioned software citation. This page does not supply a
placeholder DOI or assert that a CPC article has been published.

## Reproducibility checklist

For published calculations, retain:

- source revision, local patch if any, and executable checksum;
- compiler version, compilation flags, BLAS/LAPACK implementation and ABI;
- operating system, processor, OpenMP and BLAS thread settings;
- all input and interaction files, including solver and symmetry parameters;
- enabled `ED2_*` options and the Hamiltonian-cache diagnostics;
- complete solver logs, convergence/residual checks, and output files;
- commands and results for the regression tests and independent reference
  calculations relevant to the model.

Historical timings and figures require their own source revision,
hardware, inputs, scripts, and raw data. A current successful build does
not retroactively establish that provenance.
