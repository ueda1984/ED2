# ED2 CPC CPiP submission gates — working plan (2026-11-19)

**State: manuscript review in progress; NOT submission-ready.**
This is the first proposed ED2 journal-program article, not a claim that
ED2 was previously unpublished as *source code*. The relevant published
software predecessor is QS^3, *Computer Physics Communications* **277**,
108369 (2022), DOI [10.1016/j.cpc.2022.108369](https://doi.org/10.1016/j.cpc.2022.108369);
its code is [QS-Cube/ED](https://github.com/QS-Cube/ED).
The generalized ED2 source is already public at
[QS-Cube/ED2](https://github.com/QS-Cube/ED2).
The publication-scoped working source here is a separate fork branch,
based on experimentally audited PR #17. The unmerged PRs are **not**
an official or licensed-by-the-project release.

## Journal route and editorial thesis

**Preferred route: Computer Programs in Physics (CPiP),**
subject to coauthor/editor approval. Elsevier describes CPiP as
significant open-source computer programs archived in the
CPC Program Library (Mendeley Data); CP research papers emphasize new
computational methods and implementation. The available ED2
contribution is more naturally the first complete release paper of
the generalized `QS3ED2` executable, **not** a new Lanczos algorithm
or a generic installed Fortran library.
Official publisher scope:
https://shop.elsevier.com/journals/computer-physics-communications/0010-4655

**Central claim:** ED2 generalizes the previously published QS^3/ED
program from spin-1/2 near-saturation XXZ models to heterogeneous
local spins, non-U(1) XYZ, directed DM and symmetric off-diagonal
exchange, Cartesian fields/observables and compatible user-specified
symmetry sectors. It documents exact versus projected Hilbert
subspaces, physical input and sign conventions, numerical acceptance
and independent reproducibility. These model features were
present in the earlier publicly released ED2 source: do not attribute
their invention to the current fork's numerical-reliability PRs.
Do not claim dynamic structure factor feature parity, non-Abelian
irreps in general or a new solver.

## Staged gates and artifacts

| Gate | Criterion of completion | Existing evidence / action | Status |
| --- | --- | --- | --- |
| G0. Article category and novelty | Coauthors approve CPiP and the scope relative to 2022 QS^3; explicit comparison with QuSpin/H-Phi/ALPS | `paper/cpc_ed2_main.tex`; `docs/development/cpc-qs3-ed-to-ed2-scope-matrix.md` | **Draft; author sign-off needed** |
| G1. Source/physics correctness | Clean source distribution, repeatable independent spectra, observable identities and residual acceptance on the final frozen code | 32 old/new XXZ levels, 15 momentum levels, mixed-spin 72 levels, 60 physical DM runs, 12 literature-Hamiltonian tests, existing tests and CI | **Strong evidence; final release rerun needed** |
| G2. Large-site noncommuting symmetry example | Independent sector dimensions, native numerical residuals, exact old/new basis and H-cache bits, realistic N1000 sample | PR17 L10 Ndown3: 4 irreps D=(4299,3624,2884,3252), audited 48/48 cross-run energy/residual; 36/36 basis, 28/28 full H-cache pairs per compiler in CI | **Numerically verified; no independent N1000/3 spectrum** |
| G3. Dedicated-node performance | Freeze hardware/compiler/BLAS, source and input SHA, independent exclusivity evidence, full raw trials with multi-allocation variability; at least 2 model families | PRESTO jobs 22282 OFF and 22322 ON, 1 warmup + 3 repeats per 12 configurations, 48/48 physics PASS, basis 2.26–2.34x, total 1.44–1.71x; both operator review eligible=False | **BLOCKED for journal-speedup claims** |
| G4. Program packaging and examples | Final upstream-reviewed snapshot, clean `make distcheck` and source-archive tests, three representative model input/output fixtures, checksums, user instructions | Autotools, quickstart example, tests and CI source distribution already exist | **Needs final frozen release bundle** |
| G5. Credit, license and version | Complete author order, affiliations, CRediT, QS^3 historical attribution, MIT copyright notices, version/tag/date consistency | `configure.ac` states 1.0.0, `CITATION.cff` states 0.1.0 with date 2026-01-01; author placeholder in TeX | **BLOCKED: author/maintainer decisions** |
| G6. CPC code archive and DOI | Final approved code deposited/archived as CPiP requires, Mendeley Data/CPC Program Library accession linked to paper; optional Zenodo DOI separately tracked | Public development GitHub and MIT license | **BLOCKED: no official tag, accession or DOI** |
| G7. Manuscript technical QA | `elsarticle` builds twice, full independent DM figure checker, PR17 pilot CSV/table checker, no broken citations, manual visual proof and reference audit | `paper/cpc_ed2_main.tex`, `paper/check_dm_figure.py`, `paper/check_pr17_presto_table.py`, dedicated manuscript CI | **Automated checks in progress; human proof pending** |
| G8. Submission package | Program Summary, source distribution + examples, cover letter, data/code statement, funding/competing-interest disclosures, approved authorship and formal author consent | Provisional Program Summary exists | **BLOCKED until G4–G7** |

The **software must not be described as an installed library** unless a
public installed/ABI-stable Fortran library target is implemented and
validated. Existing Automake builds the `QS3ED2` executable.

### What is already adequate for the CPC scientific story

- **Generalized Hamiltonian and basis** relative to the *published*
  earlier program, with explicit oriented DM and nonconserving sectors.
- **Physical example with independent truth:** alternating S=1/2,S=1
  chain, 2-site-stationary XYZ/Gamma/staggered-DM interactions and
  field; verified vector chirality/Hellmann--Feynman/susceptibility.
- **Reproducibility and continuity:** direct old-QS^3 code comparisons
  plus independently assembled spin Hamiltonians; near-degenerate
  eigensubspace caveats and explicit residual acceptance.
- **Practical implementation characteristics:** dense/Krylov,
  OpenMP, mixed-spin full-space and optional sparse transition caches.
- **Noncentral yet useful demonstration:** four 1D cubic O_h sectors
  at N1000/Ndown3; PR17 paired physical equivalence and exploratory
  median timings. See `paper/data/README.md`.

### What can stay out of the critical path

- **PR #18** (opt-in 3-element sorter): exact but compiler-dependent
  and potentially slower; omit from the proposed final release
  feature claims until portable improvement is established.
- **PR #19** (exact canonical pruning): promising experiment; finish
  standalone CI/real PRESTO OFF/ON audit before considering inclusion.
  A CPC paper can proceed using validated **PR #17** as its
  feature baseline. Never merge PR #19 merely to add one extra speedup.
- New magnetic phases, full thermodynamic extrapolations, an
  independent N1000 three-magnon full spectrum, native DSF feature
  parity or a newly invented Lanczos solver **are not required claims**.

## Suggested execution sequence

1. **Now:** maintain this dedicated publication PR based on the
   PR17-reviewed source; update the TeX manuscript and capture the
   12-row PRESTO pilot with provenance, while retaining the warning
   that the full two-allocation speedup is preliminary.
2. **Next scientific freeze:** choose reviewed source features and
   exact release candidate; avoid integrating experimental sort3
   or canonical pruning automatically. Re-run clean source-archive
   tests and independent spectra/DM/Hellmann--Feynman references.
3. **Parallel human decisions:** agree on software/primary authors,
   complete affiliations and author contributions, upstream owner
   approval and exact release tag; reconcile 1.0.0 vs 0.1.0 and the
   CFF release date against a documented real release.
4. **Performance evidence:** obtain independent operator proof of
   exclusive hardware; repeat frozen workload(s) in multiple
   allocations with >=2 warmups and >=7 independent measurements
   per condition or a preregistered equivalent design; include at
   least one non-U(1) mixed-spin case. Report medians and MADs,
   raw observations, RSS, cache capacity and H-application counts,
   not just favorable ratios.
5. **Freeze/archiving:** prepare a complete *reproducible* source
   archive (program source, build scripts, inputs, selected outputs,
   expected checksums and short user guide), publish only after
   owner/author approvals and deposit according to CPC instructions.
   Do not invent an accession or DOI; a CPC Program Library record
   and optional Zenodo software DOI are different metadata objects.
6. **Submission:** compile and visually inspect each manuscript
   page, check captions/data, Program Summary and references,
   complete disclosures/cover letter, approve the submission
   package together with all coauthors, then submit as CPiP.

## Formal gate audit

Run in the source root:

    python3 -S paper/check_pr17_presto_table.py
    python3 -S paper/check_dm_figure.py --csv docs/development/cpc-stationary-dm-reference.csv --tex paper/cpc_ed2_main.tex --self-test
    python3 -S paper/audit_cpc_submission.py --output /path/to/review/report.json

The last command **reports** outstanding decisions but intentionally
does not modify source, release a package, or claim readiness.
For a genuinely frozen release candidate, run it with
`--require-ready`; this returns nonzero until the hard release gates
have been satisfied. A green LaTeX CI is not a substitute for these
release/consent checks.

**Human approval required:** publication route, final contributions,
authorship order, upstream release plan and final merge/tag. We do
not take any of these irreversible steps automatically.
