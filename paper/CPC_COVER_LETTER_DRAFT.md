# CPC CPiP cover letter — author-review draft (NOT submitted)

**Manuscript title:** ED2: Generalized exact diagonalization of mixed-spin
quantum Hamiltonians with anisotropic and nonconserving interactions

**Proposed manuscript class:** Computer Programs in Physics (CPiP).

**To:** Editors, Computer Physics Communications

Dear Editors,

We wish to submit the manuscript entitled "ED2: Generalized exact
diagonalization of mixed-spin quantum Hamiltonians with anisotropic
and nonconserving interactions" for consideration as a
*Computer Programs in Physics* article in
*Computer Physics Communications*.

ED2 is a Fortran-based exact-diagonalization program for quantum-spin
lattice Hamiltonians with heterogeneous local spin magnitudes,
three-axis exchange anisotropy, directed Dzyaloshinskii--Moriya
interactions, symmetric off-diagonal exchange and site-dependent
magnetic fields. Its Hilbert-space representation supports full
nonconserving calculations as well as occupation restrictions and
compatible spatial-symmetry sectors. The paper provides an explicit
description of these model and input conventions and demonstrates
independent verification of eigenenergies and Cartesian spin
observables.

**Relationship to previously published software.** This is the
first proposed journal-program article describing the generalized
ED2 executable, whose source has already been publicly available
in the separate QS-Cube/ED2 repository. Our earlier paper,
H. Ueda, S. Yunoki and T. Shimokawa,
"Quantum spin solver near saturation: QS3",
*Computer Physics Communications* **277**, 108369 (2022),
DOI 10.1016/j.cpc.2022.108369,
described the distinct QS-Cube/ED code for spin-1/2 XXZ-type
systems near saturation. We identify what ED2 adds relative
to that published program without claiming the development of
a new Lanczos algorithm, general superiority to QuSpin/ALPS/H-Phi,
or feature parity with the original dynamical-response pathway.

**Numerical and physical validation.** The manuscript directly
checks the common XXZ limits against the published QS3 executable
and independently constructed Hamiltonians, including
32 level-by-level and 15 momentum-sector eigenvalue comparisons.
For a genuinely non-U(1) alternating spin-1/2/spin-1 chain,
we compare finite-size spectra, correlators, the
Dzyaloshinskii--Moriya-conjugate vector chirality and its
Hellmann--Feynman derivative against independent calculations.
We additionally include a source-pinned, four-irrep
1000-site/three-magnon cubic symmetry stress test in which
48 paired production solver trials preserve the numerical
results under a carefully isolated, optional basis-coordinate
optimization. The performance data are explicitly marked as
exploratory rather than operator-attested benchmark estimates.

We believe this combination of expanded scientific model
coverage, explicit symmetry and truncation contracts,
reproducibility tools and independent physical verification
will be useful to researchers investigating finite quantum-spin
systems with complex anisotropic interactions.

**Software availability and reproducibility.**
The current development source is available under the MIT license;
the final author-approved version, program package, representative
input/output files, citation metadata and appropriate
CPC Program Library archive record will accompany the
submission in accordance with the journal's requirements.

Thank you for considering our manuscript.

Sincerely,

[Corresponding author and coauthors: approval pending]

---

## Required before using this letter

- Final author order, corresponding author, institutional affiliation
  and agreement by all listed authors.
- Final manuscript version/title and actual program release/tag.
- Confirm the scope, licensing/copyright and official
  QS-Cube/ED2 provenance with upstream maintainers.
- Confirm CPC CPiP submission procedure and required data/code
  availability declarations at submission time.
- Freeze source and independent reproducibility package;
  link the actual CPC library/Mendeley archive DOI only
  when assigned; do not invent a DOI.
- Recheck that the 1000-site performance paragraph does not
  accidentally imply a scheduler-verified publication speedup.
