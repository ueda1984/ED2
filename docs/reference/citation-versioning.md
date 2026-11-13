# Citation and Versioning

This page explains how to **cite ED2** and describes the **versioning policy**
used to ensure long-term reproducibility.

Clear citation and versioning are essential requirements for a
*Computer Physics Communications (CPC)* software article.

---

## How to cite ED2

If you use ED2 in academic work, please cite the software using one of the
methods below.

### Preferred citation (after DOI release)

Once a DOI is available, cite ED2 as:

> ED2 (QS³-ED2): Exact Diagonalization for Quantum Spin Systems,  
> Hiroshi Ueda *et al.*, Computer Physics Communications,  
> DOI: `<DOI>`.

**Two DOI types must not be conflated.** The journal assigns a DOI to an
accepted CPC article; Zenodo can assign a **different software-archive DOI**.
Do not cite a prospective article as published until publication is confirmed.

---

### Temporary citation (before DOI release)

Until the CPC article and DOI are published, please cite ED2 using the GitHub
repository and a specific version identifier:

> ED2 (QS³-ED2), QS-Cube/ED2, GitHub repository,  
> version `<tag>` or commit `<hash>`.

Always include the exact version or commit hash used in your calculations.

---

## Versioning policy

ED2 **proposes** a semantic versioning–inspired scheme for future releases:

```
vMAJOR.MINOR.PATCH
```

### Version number meaning

- **MAJOR**  
  Incremented when backward-incompatible changes are introduced
  (e.g., input format changes).

- **MINOR**  
  Incremented when new features or models are added in a backward-compatible way.

- **PATCH**  
  Incremented for bug fixes, performance improvements, or documentation updates
  that do not change user-facing behavior.

---

## Release tags

Proposed official releases would use **Git tags** in the GitHub repository (illustrative names, **not assertions that these tags exist**):

```
v1.0.0
v1.1.0
v1.1.1
```

Once the release workflow has been approved, each release tag should correspond to:
- a frozen code state,
- matching documentation and automated test evidence,
- an archival software DOI if enabled through Zenodo.

---

## Development versions

Commits between tagged releases are considered **development versions**.

For reproducibility:
- always record the full Git commit hash,
- avoid citing branch names such as `main` or `develop`.

---

## Zenodo integration

ED2 is intended to be archived on **Zenodo** to obtain a DOI.

Recommended procedure:
1. Connect the GitHub repository to Zenodo.
2. Enable DOI creation for new releases.
3. Create a Git tag for each official release.
4. Include the DOI in the documentation and CPC manuscript.

---

## CITATION.cff file

The repository includes a `CITATION.cff` file to support automated citation
by tools such as GitHub and reference managers.

!!! warning "Metadata reconciliation required before a tagged release"
    `configure.ac` currently declares `1.0.0`, while `CITATION.cff`
    declares `0.1.0` and `date-released: 2026-01-01`. These values have
    not been established as one coherent released version. The README
    lists two maintainers, while the citation file names three authors;
    the differing roles must be confirmed, without deleting attribution.
    Until a release decision is recorded, cite the exact source commit.

Users can obtain citation information directly from the repository page.

---

## Reproducibility checklist

For published results using ED2, please record:

- ED2 version or Git commit hash,
- compiler name and version,
- BLAS/LAPACK backend,
- OpenMP settings,
- input files and solver parameters.

---

## Next steps

- Add the Zenodo DOI once available.
- Update the CPC manuscript with the final citation.
- Maintain this page as new releases are published.
