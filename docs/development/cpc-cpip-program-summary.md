# ED2: draft CPC program summary and submission-position decision

> **DRAFT FOR AUTHOR REVIEW ONLY.** Not a release, software DOI,
> submitted article, accepted CPiP designation or completed journal
> metadata form. The journal category, author list and source version
> require the project maintainers' approval.

## CPC journal positioning

The publisher describes two paper types:

- **Computer Programs in Physics (CPiP):** significant executable
  programs intended for the CPC Program Library on Mendeley Data,
  with an approved open-source license.
- **Computational Physics Papers (CP):** methodological or
  computational-physics research, normally with software and
  performance evidence.

Official CPC scope:
https://shop.elsevier.com/journals/computer-physics-communications/0010-4655 .
Official CPC Program Library collection:
https://data.mendeley.com/journal/00104655 .

**Working recommendation:** position ED2 as a **CPiP
software-program article**, subject to editorial review.
The evidence here consists mainly of expanded mixed-spin and
nonconserving-Hamiltonian model coverage, observables and
independent software validation, not a new general iterative
eigensolver or asymptotic method. A CP paper would need a
separately articulated and demonstrable methodological novelty.
The current installed deliverable is an executable, not a
public reusable Fortran library.

The CPiP claim must distinguish (i) the earlier
CPC 2022 near-saturation QS³ publication and
(ii) the already-public QS-Cube/ED2 upstream program
from (iii) subsequent development, tests and documentation
on the fork. A new scientific result about Hellmann–Feynman
or a universally optimal `wk_dim`/`MNTE` setting is **not**
being claimed.

## Draft Program Summary fields

| Item | Provisional content | Release review |
| --- | --- | --- |
| Program title | ED2 (QS³-ED2): Generalized quantum spin exact diagonalization | Confirm official release title |
| Program version | **Unassigned** development candidate | Reconcile `configure.ac`, citation version and tag |
| Code repository | https://github.com/QS-Cube/ED2 | Confirm upstream release integration; the current work is in the fork https://github.com/ueda1984/ED2 |
| CPC Program Library accession | Not assigned | Do not invent an accession or DOI |
| Software DOI | Not assigned | Decide software archive and versioning |
| License | MIT; © 2023 QS-Cube | Preserve the notice and upstream credit |
| Programming language | Fortran | Validate claimed standard/compiler portability |
| Program/build tools | GNU Autotools, Fortran compiler, BLAS and LAPACK; optional OpenMP shared-memory execution | Confirm installation instructions |
| Executable interface | `QS3ED2` binary reading structured input files | No stable installed library API is currently claimed |
| Independent test dependencies | NumPy/SciPy and optionally Matplotlib for test/plot scripts | Not mandatory for the program's runtime |
| Nature of problem | Lowest spectral levels and static observables of lattice spin Hamiltonians with unequal on-site spins, nonconserving exchange and fields | Narrow to genuinely supported input operators |
| Solution method | Product-state or selected Hilbert representations; Hermitian dense or Krylov diagonalization; computation of local magnetizations and Cartesian two-spin correlations | Distinguish invariant sectors and projected truncations |
| Unusual features | Directed DM and symmetric off-diagonal exchange; user-defined spin magnitudes; optional cached Hamiltonian action, numerical acceptance controls and reproducibility sidecars | No claim of universal optimal caching |
| Typical workload | Full-space open spin-1/2/spin-1 chain through L=12 (dimension 46,656) and independent physical response checks | Controlled-node benchmark to follow |
| Known limitations | Exponential full-space growth, observable sensitivity to small gaps, no verified Intel portability in this workflow, no MPI or stable installed library ABI | Confirm final constraints before submission |

## Proposed short prose (not yet manuscript metadata)

**Nature of problem:** Standard exact diagonalization in fixed
magnetization sectors does not directly represent generic
quantum-spin Hamiltonians that combine unequal local spin
magnitudes, anisotropic exchange, antisymmetric DM terms,
symmetric off-diagonal exchange and transverse fields.
These nonconserving interactions require careful basis
construction, spectrum validation and observable checks.

**Solution method:** ED2 constructs a finite spin Hilbert
representation from on-site spin values and user-defined
interaction lists and applies dense LAPACK or iterative
Hermitian eigensolvers. Optional symmetry representations
and occupation windows reduce the working problem;
in noninvariant cases the latter represent a projection
`PHP`, not an exact symmetry sector. Local magnetization
and two-site Cartesian correlation files provide physical
outputs. Source-archive builds, independently constructed
dense and sparse Hamiltonians, solver residual diagnostics
and reproducible starting-vector controls support numerical
acceptance.

**Scientific example:** The main reproducible example
is a size-stationary two-site-periodic, alternating
spin-1/2–spin-1 chain with XYZ, Gamma_xy, transverse/longitudinal
fields and oriented staggered DM. The conjugate vector
chirality is independently compared with the energy
derivative and a finite-size susceptibility. An additional
periodic spin-1/2 XYZ-DM Hamiltonian checks published
spin conventions and oriented boundary bonds.

**Memory and timing:** The original ED2 program already
supports `wk_dim` caching and `MNTE` row capacity.
For this physical example, an exact model-specific
maximum `MNTE=7L/2-1` follows from counting allowed
transverse and two-spin-flip channels and is verified
independently through L=12. Measured memory savings
match the cache-array allocation formula; shared-runner
timings are characterized, but dedicated-node performance
claims remain pending.

## Verified software status and open release decisions

At the time of this draft:
- A reproducible development workflow validates a 60-case
  stationary DM/chirality response, 12 literature-Hamiltonian
  cases, 36 initial matched-resource trials, 52 cache-policy
  trials, 36 two-eigenpair restart trials and 72 full-cache
  row-capacity trials. **Do not add these counts blindly:**
  some categories reuse models or runs but test different
  properties and source commits.
- Archived GNU Fortran 11 and 13 source-package verification
  and the nine historical examples were completed during
  Phase 2-A/2-B; current physics benchmark CI specifically
  uses GNU Fortran 13 on Ubuntu.
- `source/` and the public fork's `main` are not being changed
  for the resource test suite; no official release tag is
  being created.

The metadata mismatch remains unresolved:
`configure.ac` declares `QS3-ED2 1.0.0`,
whereas `CITATION.cff` declares `0.1.0` with
`date-released: 2026-01-01`, without verification
of that release. The MIT license attributes © 2023
QS-Cube; the existing citation metadata lists Hiroshi Ueda,
Daisuke Yamamoto and Tokuro Shimokawa. **Do not rewrite
authorship, dates or license metadata automatically.**
Confirm them with the maintainers before packaging the
CPiP program archive.

## Submission-ready checklist

- Confirm the CPC track with the handling editor/author team;
  verify latest author instructions and Program Summary format
  at actual submission.
- Agree exact scientific/software contribution and differences
  from the published QS³ program and public ED2 upstream.
- Confirm author order and contributions, project ownership,
  appropriate citations and the existing MIT notice.
- Reconcile single authoritative software version, source
  tag, release date, citation file and archived source digest.
- Freeze a source tarball, input decks, selected expected
  outputs, test scripts and same-source validation logs.
- Obtain repeated pinned-node timing/RSS/H-application
  benchmarks with a second gap/eigenpair workload and
  honest nonconvergence reporting.
- Assemble the CPC Program Summary, statement of software
  availability, reproducibility, dependencies and restrictions.
- Never claim a Mendeley Data DOI or a journal acceptance
  before officially assigned.
