# ED2 — first CPC Computer Programs in Physics manuscript

**Working author-review manuscript. Not submitted, accepted, tagged or archived.**

This directory is a **paper-only port** of the earlier publication
draft (historical Draft PR #20), rebased conceptually on the
general-purpose `ueda1984/ED2` mainline as of 19 November 2026.
It does not import the large experimental cubic-kernel branch stack.

The intended article category is *Computer Programs in Physics*
(CPiP) of [Computer Physics Communications](
https://shop.elsevier.com/journals/computer-physics-communications/0010-4655).
The proposed novelty is **ED2's validated generalization** of the
previously published QS³ (CPC 277, 108369, 2022) to mixed local
spin magnitudes, general real/complex anisotropic spin interactions,
nonconserved (S^z), configurable compatible symmetries and complete
Cartesian static observables. Standard Lanczos or CSR are **not**
presented as newly invented numerical algorithms.

## Manuscript files and proof

- `cpc_ed2_main.tex`: complete Elsevier `elsarticle` review draft,
  with Program Summary and a clearly provisional author list.
- `data/cpc-stationary-dm-reference.csv`: archived independent
  scientific reference for 25 chirality/gap points, five susceptibility
  points and rounded table cells; this is previous scientific
  evidence and must be revalidated on the final source revision.
- `check_dm_figure.py`: checks every plotted/tabled physics value,
  including deliberate corruption self-tests.
- `audit_reference_small_systems.py`: independent NumPy/SciPy
  Hamiltonian and response reconstruction at N4, N6, N8
  (not an ED2 operator implementation).
- `tools/cpc_portable_dm_quickstart.py`: end-user reproduction
  script using only Python standard library and the compiled ED2
  program; it checks the chirality/energy-response identity and
  reference values for N4/N6.
- `check_current_main_native_dm.py`: matches the mainline four/six-site native DM table to the archived CI-produced source/binary-provenance JSON, with a deliberately corrupted value rejected.
- `MAINLINE_NATIVE_DM_REPRODUCTION.md`: native current-main four/six-site energy, staggered DM chirality, finite-difference Hellmann--Feynman, source/binary checksums and **no direct residual sidecar** caveat.
- `check_mainline_cache_claims.py`: checks documented generic
  mainline opt-ins and exactly recomputes the new cache-allocation
  table from the recorded spin-model dimensions/record count. It
  **does not** claim independent physical generation of Nrec=191808.
- `data/pr17_presto_coordinate_cache_pilot.csv` and
  `check_pr17_presto_table.py`: **unmerged experimental branch**
  provenance. Its source was pinned to PR17, **not** the current
  mainline. This table now appears only in the manuscript appendix,
  not as a program-wide speedup or a feature of the official release.
- `CPC_SUBMISSION_GATES.md`: authoritative list of remaining
  decisions before journal submission.
- `CPC_PROGRAM_PACKAGE_PLAN.md`: intended CPiP source archive and
  example distribution, pending approval.
- `CPC_COVER_LETTER_DRAFT.md`: editorial proposal for coauthor review.

## Local editorial checks

From the repository root:

```bash
python3 -S paper/check_dm_figure.py --self-test
python3 -S paper/check_pr17_presto_table.py
python3 -S paper/check_mainline_cache_claims.py --self-test
python3 -S paper/audit_cpc_submission.py --output paper-readiness.json
```

The independent numerical physics reference additionally requires
NumPy and SciPy:

```bash
python3 paper/audit_reference_small_systems.py \
  --csv paper/data/cpc-stationary-dm-reference.csv
```

The portable end-user demo, after a clean ED2 build, is:

```bash
python3 -S paper/tools/cpc_portable_dm_quickstart.py \
  --exe "$(pwd)/source/QS3ED2" --size 4 --omp 1 \
  --output "$(pwd)/paper-dm-example-l4"
```

Use a **new** output directory. Tests check physical response and
sign conventions; if the historical TRLan solver fails to converge
on a particular compiler/start, do not lower accuracy requirements
to force a pass. Report the issue and pin the source/compiler/seed.

Compile and inspect the manuscript with a complete LaTeX installation:

```bash
mkdir -p paper/build
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory=paper/build paper/cpc_ed2_main.tex
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory=paper/build paper/cpc_ed2_main.tex
```

The paper CI records source/PDF hashes and preserves a preview.
A successful build proves only typesetting integrity and automated
numerical transcription, **not** publication readiness.

## Current mainline cache options in the paper

General source code supports `ED2_CACHE_COEFF_BOUND=1`,
`ED2_CSR_COMPLEX_CACHE=1`,
`ED2_CSR_RELEASE_PADDED_CACHE=1`, and
`ED2_CSR_DIRECT_BUILD=1`. All are default OFF.
The new CPC table documents **cache-array payload**, not
measured peak RSS or a portable timing speedup. Full/partial
cache eligibility and the more exact user instructions are in
`docs/user-guide/performance-and-memory.md`.

The special N1000 cubic (\Gamma/O_h) coordinate-cache pilot was
an audited numerical comparison between two PRESTO allocations,
but hardware-exclusivity confirmation remains unavailable.
It cannot be used to claim that the final mainline release
achieves the same timing ratios.

## Outstanding external decisions

The complete author list/order/affiliations and CRediT, official
QS-Cube maintainer approval, a frozen program source/tag,
author-approved release metadata, exact publication category,
CPC Program Library/Mendeley Data accession and DOI, and final
editorial declarations are **not yet resolved**. The repository
currently has a metadata discrepancy (`configure.ac` **1.0.0**
versus `CITATION.cff` **0.1.0** with a date that is not independently
verified). The manuscript intentionally does not conceal it.

The [published Read the Docs site](
https://ed2.readthedocs.io/en/latest/) follows
`QS-Cube/ED2`, not automatically the newer
`ueda1984/ED2` development fork. Even merged fork
documentation requires an upstream-authorized merge and a
fresh Read the Docs build to reach users.

**Do not publish, deposit, tag, assign a DOI, or claim final
authorship from this draft.**
