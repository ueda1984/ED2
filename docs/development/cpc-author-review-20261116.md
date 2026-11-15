# ED2 first CPC paper: author-review checkpoint (2026-11-16)

**Manuscript basis:** `paper/cpc_ed2_main.tex` at development commit
`23312f880db98299defa0a0510f53f71cb552ee5`. Later commits may
change documentation or build status; this is **not** an immutable release
or a formal editor's verdict.

**Purpose:** decide whether the manuscript is mature enough for PI/author
scientific review. This rubric is an *internal editorial self-assessment*,
not a prediction of journal acceptance, a peer-review decision, or a
submission authorization.

## Multi-angle CPC CPiP readiness assessment

| Dimension | Weight | Internal score | Evidence and explicit reservation |
| --- | ---: | ---: | --- |
| Scientific/software contribution and novelty positioning | 20 | 17 | A first ED2 journal-program description relative to CPC-published QS³/ED; audited pre-existing upstream ED2 features are credited; HΦ, ALPS and QuSpin are acknowledged. Novelty is program-specific, not a new Lanczos algorithm. |
| Physics correctness and independent validation | 20 | 19 | 32 prior-code XXZ levels, 15 momentum levels, 72 bounded-basis levels, 60 stationary mixed-spin DM trials, 12 literature-Hamiltonian checks, and near-degenerate projector tests; cases overlap in properties and must not be summed as independent problems. |
| Reproduction and numerical software reliability | 20 | 18 | Small standard-library-only user driver, independent SciPy matrices, physical/figure integrity gates, residual contracts, archived CI outputs; model-specific memory counting. The 12-site scientific data are independently verified in the earlier scientific workflow, not by the new 4--8-site audit alone. |
| Manuscript reasoning, figures and scholarly positioning | 20 | 18 | Four figures (workflow, response, gap, susceptibility), program-feature and validation matrices, model definition, Pauli-to-spin literature mapping, empirical cache/RSS caveats, full references and a reproduction recipe. Figure placement, English and numerical prose still need author review. |
| Release, benchmarking and program-library packaging | 20 | 10 | Provisional Program Summary exists, but authorship, official merged source, release version/tag, DOI, program archive, CPC article category and dedicated-node benchmarks remain unresolved. |
| **Total (draft review readiness only)** | **100** | **82** | **The manuscript has crossed an internal 80-point threshold for author scientific review. It is not submission-ready.** |

## Three distinct validation contracts

1. **Physics truth within the tested finite models.** Independent spin
   matrices and references compare eigenvalues, gaps, observables and
   Hellmann--Feynman identities against ED2. For the primary scan,
   L=4--12 are covered by the prior independent scientific workflow
   (60 Fortran trials, 25 numerical reference parameter points).
2. **Fresh small-system source-data audit.** The separate
   `paper/audit_reference_small_systems.py` reconstructs and checks
   L=4,6,8 at five DM couplings (15 source points) and three zero-D
   energy-derivative/susceptibility identities. This overlaps the
   primary scan and is deliberately **not** counted as 15 novel
   scientific cases.
3. **Typesetting integrity.** `paper/check_dm_figure.py` checks
   25 chirality + 25 first-gap + five susceptibility plotted ordinates
   and five rounded table entries. Four deliberately tampered values
   must be rejected. A successful figure gate does not independently
   verify physics.

## Likely referee/editor objections and current answers

**Why publish ED2 if mixed-spin ED already exists?**
ED2's legitimate paper is its *first formal program description*,
relative to the already-CPC-published QS³/ED code; most generalized
operators already existed in the public upstream ED2 before the present
fork's correctness work. Source-level engineering and the test/release
infrastructure must not be misrepresented as initial invention.
Wider ED packages (HΦ, ALPS, QuSpin) can already solve overlapping
physical models. The manuscript emphasizes a specific Fortran program,
its input/basis formalism and independently verified numerical scope.

**Is the scientific model physically novel?**
No. The two-site-periodic, explicitly non-U(1) mixed-spin chain is a
repeatable software-validation application. Chirality, susceptibility
and gaps are finite-size observables; no new magnetic phase or bulk
extrapolation is claimed. Yi et al.'s model is reproduced at the
Hamiltonian level with a Pauli/spin conversion, not by redrawing its
published phase diagram.

**Are the timing results reproducible speedup claims?**
No. The source-archived CI runner is shared and the measurements have
only two measured repeats. Nominal cache allocations are derived,
and selected RSS values are *observed* as an exploratory check.
Any journal-level speedup/scaling claim requires a verified exclusive
node, pinned threads, provenance, multiple repeats and error bars.

**Do small residuals prove correct ground-state observables?**
No. The paper distinguishes residual accuracy from ground-state
inclusion and uses independent low-energy references and rank-two
projector traces for near-degenerate doublets.

**Is ED2 a functional superset of original QS³?**
Not established. The published original supports a native dynamical
structure factor route that has not been demonstrated in ED2;
continuity claims are explicitly restricted to shared spectra and
static observables.

## Submission blockers and decisions reserved for authors

- Agree author list, order, affiliations, contributor credits,
  provenance and copyright/license attribution.
- Set one consistent version and release date: the observed
  `configure.ac` and `CITATION.cff` metadata are inconsistent.
- Decide official upstream integration, frozen source SHA/tag,
  archived source package and CPC Program Library/DOI workflow.
  Do not invent an accession or assign a DOI before it exists.
- Complete controlled-node performance characterization with
  independent verification of machine exclusivity and a second
  spectral workload. Shared-runner timings cannot fill this role.
- Read and check the compiled proof page-by-page, choose which
  finite-size figures belong in main text versus supplement, and
  recheck the latest CPC author-facing requirements at submission.

**Recommended current action:** PI review of the physics,
the claim of novelty relative to the already-public ED2 implementation,
and the editorial usefulness of the three physical plots. Metadata,
release and controlled-machine benchmarking can then be frozen for a
genuine submission candidate.
