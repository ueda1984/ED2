# CPC manuscript plan: ED2 versus the published QS³ program

> **Working outline, not a submission or a verified novelty claim.**
> Prepared during Phase 2-B from the audited ED2 development source.
> Authors, title, exact numerical tables and the choice of CPC article
> category remain subject to review. This draft is **not** a released
> program or a peer-reviewed manuscript.

## Scientifically defensible central claim

**Proposed thesis:** ED2 extends the symmetry-adapted, near-saturation
quantum spin solver QS³ toward exact diagonalization of more general
quantum-spin Hamiltonians, including site-dependent spin magnitudes,
anisotropic/DM/symmetric off-diagonal exchange and transverse fields,
while explicitly distinguishing invariant symmetry sectors from
nonconserving excitation-number truncations and recording independent
numerical verification.

The main contribution to a *CPC* submission should not be described
as merely “the Lanczos algorithm” or “OpenMP parallelization”: the
published QS³ program already has Lanczos/thick-restart Lanczos,
one- and two-point observables, static/dynamic structure factors and
OpenMP performance analysis.

Prior primary reference: **H. Ueda, S. Yunoki and T. Shimokawa**,
“Quantum spin solver near saturation: QS³,”
*Computer Physics Communications* **277**, 108369 (2022),
[doi:10.1016/j.cpc.2022.108369](https://doi.org/10.1016/j.cpc.2022.108369).
The publication's reference implementation is associated with
[QS-Cube/ED](https://github.com/QS-Cube/ED), while this project's
upstream is [QS-Cube/ED2](https://github.com/QS-Cube/ED2).
These are **distinct repositories**.

The journal identifies **Computer Programs in Physics (CPiP)** for
significant open-source programs archived in its Mendeley Data program
library and **Computational Physics Papers (CP)** for methods-oriented
research. The intended category requires editorial and author review:
[CPC aims and scope](https://shop.elsevier.com/journals/computer-physics-communications/0010-4655).

## Distinguish implementation, validation and novelty

| Capability | 2022 QS³ paper | Current ED2 development source | What remains to establish |
| --- | --- | --- | --- |
| Local spin values | Spin-1/2 framework | Site-specific spin values (including mixed spin) | Explicit small-system and scaling case |
| Exchange models | XXZ-type near saturation | XYZ, DM and symmetric off-diagonal Gamma terms | Reproducible physically motivated example |
| Excitation-number conservation | U(1) sector is essential to published scope | Full sectors and optionally restricted NOD interval | Explain projected vs invariant sector for non-U(1) terms |
| Lattice symmetries | Momentum sectors for translation | User-provided site permutations/characters | Check compatible one-dimensional sectors; no arbitrary non-Abelian claim |
| Eigensolvers | Lanczos and thick-restart Lanczos already published | ALG=1/2 and dense ALG=3 with tested residual/status handling | Avoid claiming solver families themselves are new |
| Observables | Local, two-point and static/dynamic structure factors | Local magnetization and static Cartesian/ladder two-body correlations verified | Do **not** assert dynamic structure factors exist in ED2 without code evidence |
| Software quality | Program published in CPC | Source distribution, diagnostics, replay and numerical acceptance contract | Engineering value versus *scientific* novelty must be argued separately |
| Installed Fortran library API | Not the claim here | Current Automake target installs \`QS3ED2\` executable | No public library API is currently shipped |

**Priority:** independently demonstrate an important model that the
published QS³ formulation does not directly cover, rather than
reproducing its original near-saturation benchmarks alone.

## Candidate title

**ED2: Exact diagonalization of generalized quantum spin Hamiltonians
with symmetry-adapted and truncated Hilbert spaces**

Avoid the term “Fortran library” in the final title until an installed,
documented, externally exercised reusable API exists.

## Abstract scaffold (rewrite after completing physics results)

> We present ED2, a Fortran exact-diagonalization program for quantum
> spin lattice models that extends the near-saturation QS³ framework to
> [precisely describe verified local-spin and interaction generality].
> The program constructs [supported basis and symmetry sectors] and
> provides [verified numerical solvers], with [specified observable
> outputs]. A documented source archive, automated reference
> comparisons and optional solver residual diagnostics support
> reproducible calculations. We demonstrate the program using [one
> physically relevant mixed-spin/non-U(1) model with parameter values]
> and benchmark [system sizes, compiler, OpenMP, memory and runtime].
> Our results establish [measured capability, not an unverified
> advantage or universal precision guarantee].

## Manuscript sections

### 1. Introduction

- Explain why general quantum spin models with mixed local spins,
  anisotropies, DM/Gamma terms and spatially varying fields remain
  important, alongside the Hilbert-space scaling challenge.
- Cite the 2022 QS³ paper precisely and state what it already provides.
- Explain **where ED2 adds scientifically usable model coverage**,
  not just improved build/CI infrastructure.
- Explain the intended users and boundary between exact full-space
  diagonalization, symmetry sectors and projected truncations.

### 2. Hamiltonians and Hilbert-space representation

Use the verified code sign convention (real couplings shown only for
illustration; consult the parser for ordering):

\[
H = \sum_{\langle i,j \rangle}
\left[
\sum_{\mu=x,y,z} J_{ij}^{\mu} S_i^{\mu} S_j^{\mu}
+\mathbf D_{ij}\cdot(\mathbf S_i\times\mathbf S_j)
+\sum_{\mu=x,y,z}\Gamma_{ij}^{\mu}
   (S_i^{\nu_\mu}S_j^{\rho_\mu}+S_i^{\rho_\mu}S_j^{\nu_\mu})
\right]
+\sum_i\mathbf h_i\cdot\mathbf S_i.
\]

Here \((\nu_x,\rho_x)=(y,z)\),
\((\nu_y,\rho_y)=(z,x)\) and
\((\nu_z,\rho_z)=(x,y)\).
**The field term is \(+\mathbf h\cdot\mathbf S\)** in the ED2 input
convention; a physical Zeeman term \(-\mathbf B\cdot\mathbf S\)
requires \(\mathbf h=-\mathbf B\). DM signs depend on directed
bonds; reversing input bonds changes the oriented DM coefficients.

For site-specific \(S_i\), the full product basis has dimension
\(\prod_i(2S_i+1)\). Explain how local occupation bounds and total
\`NODmin\`/\`NODmax\` define the represented space, and why
\(\,P H P\,\) for an NOD-restricted, non-U(1) Hamiltonian is **not**
a dynamically invariant sector in general.

Clarify supported site-permutation symmetries and one-dimensional
characters; user-provided permutations are not automatically checked
against all Hamiltonian terms or arbitrary noncommuting group
relations.

### 3. Implementation and diagonalization

- Fortran code architecture, basis enumeration and cache operation.
- ALG=1 reconstructed-vector Lanczos; ALG=2 restarted Hermitian
  solver; ALG=3 dense LAPACK.
- Optional diagnostics, controlled initial vector/seeding and
  independently specified absolute/relative residual acceptance.
- Scope limitations: a tiny residual does not certify global
  ground-state inclusion or individual observables at small gap;
  no stable library ABI/installed module interface is claimed.
- Separate **numerical algorithms** from **recent maintenance
  improvements** to avoid misrepresenting their provenance.

### 4. Independent verification and distribution

- Independent 21-case Hamiltonian/product-basis and symmetry checks.
- Sparse field checks (70 CLI, 10 operator, 28 permutation).
- Lanczos 82 kernel plus 26 CLI tests, and Phase 1-F/1-G
  diagnostic/policy acceptance cases. Distinguish *test passed* from
  *physical solver converged*.
- Archived source builds on GNU Fortran 11/13 with separate build
  directories; \`make distcheck\`, extended nine examples and
  27 historical reference comparisons (use archived CI links).
- Show examples of two numerical pitfalls: exact excited-vector
  starting state and low-gap observable errors.
- Do **not** claim identical status counts or random bitstreams
  across compilers.

### 5. Scientific application and performance (requires new data)

**Recommended compact demonstrator:**
an alternating mixed-spin chain
(\(S_{2i-1}=1/2,\ S_{2i}=1\)) with inequivalent
\(J_x,J_y,J_z\), a DM term, a nonzero \(\Gamma_{xy}\) or transverse
field, and open boundary conditions. This simultaneously exercises
general spin, interaction anisotropy, nonconserved \(S^z\) and local
observables; it is outside the original spin-1/2 U(1)-restricted
QS³ use case.

First run small \(L\) with the **complete** occupation range
\(\text{NODmin}=0,\ \text{NODmax}=\sum_i 2S_i\).
Use an independently constructed dense Hamiltonian to compare
energies, residuals, local magnetization and correlation functions.
Keep site/bond orientation and correlation-column ordering explicit.
Only then scale to larger \(L\), measuring runtime, RSS memory,
Hamiltonian actions and the effect of \`wk_dim\`, and compare
OpenMP 1 versus 2 threads.

Choose a physically motivated parameter cut and interpret at least
one nontrivial observable or avoided crossing rather than reporting
only a performance microbenchmark. Record parameters and figure
data and retain all runs, including nonconvergent trials.

A second optional example could compare an invariant U(1) sector
with a **projected** excitation-number truncation after transverse
terms are introduced, showing the approximation error rather than
calling the restricted space an exact symmetry block.

### 6. Conclusion and limitations

The conclusion should explicitly state model coverage, evidence from
the worked physics example, reproduction procedure, and unresolved
limits (very small spectral gaps, global ground-state inclusion,
compiler portability, algorithmic cost, API scope).

## Proposed figure/table inventory

1. Diagram of basis/sector construction and ED2 execution workflow.
2. Independent dense versus ED2 energies and selected observables for
   the non-U(1) mixed-spin demonstrator (with error versus system size).
3. Performance/peak-memory versus dimension, caching and OpenMP
   threads for a fixed scientifically motivated Hamiltonian.
4. Numerical accuracy demonstration separating residual, gap and
   observable errors.
5. Concise comparison table: published QS³ capabilities, verified
   ED2 capabilities, and explicitly unsupported/unverified items.

## Preconditions before a full first submission draft

- Run and archive the **scientific demonstrator** above; without
  it, novelty is a proposal rather than demonstrated physics.
- Agree author/order and distinguish the **2022 QS³** authors
  from the current software's attribution.
- Choose CPC CPiP or CP route and associated software-archive
  obligations; check current author instructions at submission.
- Decide whether any Fortran library interface is part of the paper.
- Reconcile \`AC_INIT\`, \`CITATION.cff\`, release date, tag and
  software DOI without confusing the CPC article DOI.
- Freeze input decks, result files, compiler/backend metadata,
  provenance and program summary against one source SHA.

## Existing reproducibility evidence

- [Phase 2-A archive build, nine models and historical references](https://github.com/ueda1984/ED2/actions/runs/43350378812).
- [Phase 2-B independent GNU 11/13 archive-based validation](https://github.com/ueda1984/ED2/actions/runs/43392314521).
- [Release/attribution audit](phase2b-release-audit.md).

**No manuscript authorship, official source tag, DOI or journal
publication is assigned by this planning document.**
