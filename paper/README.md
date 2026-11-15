# ED2 first CPC program-paper source (provisional)

This directory contains a **journal-formatted, author-review LaTeX draft**, not
an accepted, submitted or finalized article.

- [Main manuscript](cpc_ed2_main.tex) uses Elsevier's `elsarticle` class.
- [Paper compilation CI](https://github.com/ueda1984/ED2/actions/workflows/cpc-manuscript-latex.yml)
  first verifies that **all 25 plotted DM/chirality points** match the
  committed independent stationary-reference CSV via
  `paper/check_dm_figure.py`. It then compiles the PDF twice,
  checks it with `pdfinfo`/`pdftotext`, and archives a proof artifact
  for visual review. The workflow does not publish the PDF.
- The [compact manuscript working draft](../docs/development/cpc-submission-core-draft.md)
  contains the longer scientific argument and release TODOs.
- The [detailed technical draft](../docs/development/cpc-manuscript-draft.md)
  retains near-degenerate and engineering failure analyses that may be
  moved to supplementary material.
- The [published-QS³ to ED2 source and feature audit](../docs/development/cpc-qs3-ed-to-ed2-scope-matrix.md)
  establishes what is new **relative to the published 2022 QS³ program**.

## Reproduce locally

On a system with TeX Live, including `elsarticle`, `amsmath`, `bm`,
`hyperref` and `booktabs`:

```sh
mkdir -p paper/build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper/build paper/cpc_ed2_main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper/build paper/cpc_ed2_main.tex
```

The compiled PDF is a draft proof. All numerical values have independent
source-linked evidence; the PDF is not a substitute for the tested code
archive or the data files.

## Scientific evidence matrix

| Manuscript part | Regeneration / verification evidence |
| --- | --- |
| Published QS³ fixed-down-spin XXZ continuity | [32-level three-sector CSV](../docs/development/cpc-qs3-ed-ed2-32level-extended-ci.csv) |
| Published QS³ momentum continuity | [15-level six-sector CSV](../docs/development/cpc-qs3-ed-ed2-momentum-15level-ci.csv) |
| Conventional-Lanczos continuity, 12/18/24 sites | [Both original executable runs and independent eigenenergy checks](../docs/development/cpc-original-qs3-ed2-iterative-continuity.md) |
| Bounded-occupation formula | [Full/bounded mixed-spin basis tests, 72 levels](../docs/development/cpc-qs3-ed-to-ed2-basis-generalization.md) |
| Non-U(1) alternating-spin staggered DM response | [Physical model and 60-case verification](../docs/development/cpc-stationary-dm-application.md) |
| Literature spin-convention check | [Independent finite-size XYZ-DM verification](../docs/development/cpc-literature-matched-benchmark.md) |
| Cache capacity and resource caveats | [Cache allocation derivation](../docs/development/cpc-cache-pareto-evaluation.md) |
| Controlled performance measurements (not yet complete) | [Dedicated-node gated protocol](../docs/development/cpc-controlled-node-protocol.md) |

## Pre-submission blockers

The authorized authors and order, affiliation, official code release/tag,
program archive DOI, CPC Program Library accession and final manuscript
category are **not confirmed**. The original QS³ code contains a
dynamical-structure-factor route which has not been shown to be
natively reproduced by the current ED2 executable. Do not describe
ED2 as automatically backward-compatible or the whole output feature
set as a superset.

The first physical figure (five system sizes, 25 finite-size DM-response
points) is now embedded directly as a `pgfplots` figure. Its coordinates
must pass the independent CSV-integrity gate in the build workflow.
The full scientific figure set (in particular the excitation-gap and
finite-size susceptibility visualizations) still needs author review
for main-text versus supplementary placement. The present resource
numbers are **not** a controlled-node performance study.
Scientific results and references require final author review before
external submission.
