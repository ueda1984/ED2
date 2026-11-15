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
| Installed Fortran library API | Not the claim here | Current Automake target installs `QS3ED2` executable | No public library API is currently shipped |

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
`NODmin`/`NODmax` define the represented space, and why
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
  directories; `make distcheck`, extended nine examples and
  27 historical reference comparisons (use archived CI links).
- Show examples of two numerical pitfalls: exact excited-vector
  starting state and low-gap observable errors.
- Do **not** claim identical status counts or random bitstreams
  across compilers.

### 5. Scientific application and performance (physics demonstration verified; publication benchmark pending)

**Recommended compact demonstrator:**
an alternating mixed-spin chain
(\(S_{2i-1}=1/2,\ S_{2i}=1\)) with inequivalent
\(J_x,J_y,J_z\), a DM term, a nonzero \(\Gamma_{xy}\) or transverse
field, and open boundary conditions. This simultaneously exercises
general spin, interaction anisotropy, nonconserved \(S^z\) and local
observables; it is outside the original spin-1/2 U(1)-restricted
QS³ use case.

The initial independent [mixed-spin demonstration](cpc-mixed-spin-demo.md)
verified full-space four- and six-site chains and a four-site
projected non-U(1) space against all dense eigenvalues and
local/two-site observables (six cases). The follow-up
[independent sparse study](cpc-sparse-scaling.md)
extended fully nonconserving calculations to eight, ten and
twelve sites (Hilbert dimension up to 46,656) with lowest-pair,
magnetization and correlation checks at OMP=1 and 2.
The [one-eigenpair matched-start study](cpc-matched-krylov.md)
and [ALG2 workspace scan](cpc-alg2-workspace.md) then measured
restart parameters at fixed requested accuracy, including
two-low-eigenpair runs with independent physical references.

These are **verified implementation and resource-characterization
results**, not yet a claim of new physics or a universal
performance advantage. Separate low-gap/ground-inclusion studies now
test near-degeneracy explicitly; do not confuse those controlled
numerical stress cases with evidence of a phase transition.

A complementary **size-stationary science application** now scans
a staggered oriented \(D_z\) coupling in the same alternating
mixed-spin, non-U(1) chain. The historical regression Hamiltonian's
fields and bond coefficients varied with absolute site index, so
those runs were **not** suitable for a thermodynamic size
extrapolation. The new
[stationary DM/chirality response study](cpc-stationary-dm-application.md)
instead fixes a two-site repeat rule for \(J^{x,y,z}\),
\(\Gamma_{xy}\), \(h_x\) and \(h_z\) for every length.

The conjugate chirality
\(\sum_j(-1)^j\langle
S_j^x S_{j+1}^y-S_j^y S_{j+1}^x\rangle\)
is checked against \(\partial E_0/\partial D\)
using both an independent sparse Hamiltonian and
actual ED2 NN correlations. The positive finite-size
susceptibility \(-\partial^2 E_0/\partial D^2/(L-1)\)
is checked by energy curvature and direct chiral response;
four-/six-site models also use the complete spectral sum.
**60/60 scientific trials passed**, from \(L=4\)
to \(12\), including OMP=1 and OMP=2 at **every tested length
through \(L=12\)**.

For the CPC draft, this establishes a **reproducible
non-U(1) response example with size-stationary local
rules**, not a new chiral phase, thermodynamic
extrapolation, performance advantage or material
parameterization. The [reference CSV](cpc-stationary-dm-reference.csv)
is committed for figure reproduction.
A [Hamiltonian-level reproduction of Yi et al. (2019)](cpc-literature-matched-benchmark.md)
now validates a published XYZ/DM/field model in finite periodic chains,
with 12/12 independent ED2 comparisons successful at OMP=1/2.
This is **not** a reproduction of the article's phase diagram or its
numerical plotted datasets. The main scientific example remains the
size-stationary mixed-spin chain, and controlled-node resource
benchmarks remain a publication priority.

**Literature position (not a parameter-source claim):**
Lahiri and Pati's alternating spin-\(1/2\)/spin-1
chain with NN/NNN exchanges and \(D_z\)
([*J. Phys.: Condens. Matter*, 2022](https://doi.org/10.1088/1361-648X/ac4aab))
and Yi *et al.*'s anisotropic spin-\(1/2\) XYZ/DM/field
example ([*Phys. Rev. B* **100**, 024423 (2019)](https://doi.org/10.1103/PhysRevB.100.024423))
motivate the importance of these interaction classes.
**Neither paper supplied the present stationary
demonstrator's exact coefficients**: they should not be
cited as if we reproduced their numerical phases.
The 2022 QS³ scope comparison belongs in the
Introduction, not in a claim of new underlying chiral physics.

For any additional model, begin at small \(L\) with the **complete**
occupation range \(\text{NODmin}=0,\ \text{NODmax}=\sum_i 2S_i\).
Cross-check energies, residuals, magnetization and two-body correlation
functions against an independent dense or sparse Hamiltonian. The
existing scaling study has already recorded wall time, peak RSS and
H-action counts at OMP=1/2 for a single couplings family; it does
**not** replace a controlled repeated-node benchmark, a cache-policy
scan (`wk_dim`), or a physically diverse spectral study.

The DM-conjugate response is now a **validated illustrative
parameter cut**; retain its input files, source SHA, figure
CSV, and honest records of any nonconvergence. To argue
broader scientific relevance, identify a literature-grounded
Hamiltonian and explain why its parameter choices matter.
Do not infer a critical point from five DM values or
five open-chain sizes.

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
2. **DM/chirality figure:** staggered NN vector chirality per bond
   versus oriented DM coefficient, (L=4,6,8,10,12), plus the
   finite-size (chi_kappa(L)). Include independently checked
   (dE_0/dD), not just solver agreement.
3. Performance/peak-memory versus dimension, caching and OpenMP
   threads for a **fixed workload and matched eigenpair target**;
   obtain replicated dedicated-node timings before making speed claims.
4. Numerical accuracy demonstration separating residual, gap and
   observable errors.
5. Concise comparison table: published QS³ capabilities, verified
   ED2 capabilities, and explicitly unsupported/unverified items.

## Preconditions before a submission-ready manuscript

- **Illustrative physics application completed:** the 4–12-site
  two-site-periodic mixed-spin DM/chirality response is verified
  against independent finite-difference, sparse/dense and full
  correlation checks. An additional published-Hamiltonian
  finite-size check and a preliminary matched-target resource
  study have been verified; see the [follow-up report](cpc-literature-matched-benchmark.md).
  For submission, verify a larger length/window where feasible
  and obtain repeated isolated-node resource benchmarks. A new
  [cache and two-eigenpair study](cpc-cache-multipair-evaluation.md)
  verifies 52/52 cache-policy trials and 36/36 two-lowest-state
  trials with independent physical checks, but shared-runner times
  do not establish publication-grade speedups. No new phase or
  asymptotic numerical advantage has been established.
- Agree author/order and distinguish the **2022 QS³** authors
  from the current software's attribution.
- Choose CPC CPiP or CP route and associated software-archive
  obligations; check current author instructions at submission.
- Decide whether any Fortran library interface is part of the paper.
- Reconcile `AC_INIT`, `CITATION.cff`, release date, tag and
  software DOI without confusing the CPC article DOI.
- Freeze input decks, result files, compiler/backend metadata,
  provenance and program summary against one source SHA.

### Near-degenerate accuracy interpretation

A new [physical near-degenerate study](cpc-near-degenerate-gaps.md)
verified an **exactly twofold** or weakly split low-energy doublet
with an independent dense reference. For \(\Delta=10^{-10}\),
an ALG=1 result can pass a direct residual `1e-8`
while the unique ground-state edge-spin observable is inaccurate;
tightening the residual resolves these particular starts but
increases Hamiltonian applications.

A controlled first-excited initial vector gives a stronger general
counterexample: for \(\Delta=10^{-6}\), ALG=1 returned
a **correct excited eigenpair**, residual \(10^{-15}\), success
status zero, and the **wrong state to call the global ground**.
This is mathematical Krylov behavior, not a regression failure.
The manuscript should explicitly distinguish residual certification,
low-doublet inclusion, individual-state observable fidelity and
global ground-state identification. Do **not** advertise an
unimplemented spectral inclusion guarantee.

### Interacting avoided crossing with independently validated observables

The new [interacting avoided-crossing study](cpc-interacting-avoided-crossing.md)
extends the prior exactly decoupled low-doublet model with
\(\varepsilon(S_1^xS_2^x+\eta S_1^yS_2^y)\)
for \(\eta=0\) and \(0.6\). For \(\eta=0.6\),
the central ground state is **entangled** and no single edge-spin
component commutes with the full Hamiltonian, avoiding the
special product-state structure of the pure XX coupling.

Independent full-matrix validation confirmed the avoided gap,
edge magnetization reversal, local/correlation observables and
two-state projector/trace quantities on four and six sites,
including gaps around \(10^{-7}\). Across the two coupling
families, OMP=1/2 and several solver configurations, there are
**216 successfully classified/validated ordinary runs**.
A separate explicit interacting first-excited starting vector
also reproduced an ALG=1 residual near \(10^{-15}\)
with a first-excited rather than globally lowest energy.

A subsequent [8-/10-site interacting sparse study](cpc-interacting-sparse810.md)
extends the **same** noncommuting XY avoided-crossing model to
dimensions 1,296 and 7,776, including independent ARPACK
references for a ground doublet split by approximately \(10^{-7}\).
Across OMP=1/2, **52/52 additional trials passed**.
Near this tiny gap, the individual ED2/reference observable
difference can reach \(1.4\times10^{-7}\), while the
independent low-doublet ensemble trace differs by less
than \(6\times10^{-9}\); the reference's own residual-to-gap
sensitivity must be discussed explicitly.

These data are appropriate for **scientific verification figures**,
not evidence of a newly discovered phase transition or general
solver superiority. The existing 8--12-site sparse study validates
a different model, while the present 8/10-site test uses the
same interacting avoided-crossing Hamiltonian. More
physics-motivated application and repeated-node performance
benchmarks remain necessary.

## Existing reproducibility evidence

- [Phase 2-A archive build, nine models and historical references](https://github.com/ueda1984/ED2/actions/runs/43350378812).
- [Phase 2-B independent GNU 11/13 archive-based validation](https://github.com/ueda1984/ED2/actions/runs/43392314521).
- [Release/attribution audit](phase2b-release-audit.md).
- [Independent sparse scaling with physical tests](cpc-sparse-scaling.md).
- [Matched one-eigenpair cost and residual](cpc-matched-krylov.md).
- [ALG2 workspace and two-pair robustness study](cpc-alg2-workspace.md).
- [Physical near-degenerate and excited-start controls](cpc-near-degenerate-gaps.md).
- [Interacting avoided crossings and noncommuting bond physics](cpc-interacting-avoided-crossing.md).
- [Interacting 8-/10-site sparse spectral and observable checks](cpc-interacting-sparse810.md).
- [Size-stationary DM/chirality response, susceptibility and figure data](cpc-stationary-dm-application.md).
- [Published XYZ-DM Hamiltonian and matched-target resource verification (12+36 cases)](cpc-literature-matched-benchmark.md).
- [Cache-policy sweep and ALG=2 two-eigenpair robustness (52+36 trials)](cpc-cache-multipair-evaluation.md).
- [Exact MNTE row-capacity derivation, independent CSR sparsity and resource measurements](cpc-mnte-capacity-evaluation.md).
- [72-run frozen full-cache MNTE CI measurements](cpc-mnte-capacity-raw-ci.csv).
- [CPC CPiP positioning and draft program-summary fields](cpc-cpip-program-summary.md).
- [Eight-site sparse-reference low-doublet native-readback validation (12/12)](cpc-stationary-sparse-doublet.md).
- [Interacting eight-site XY-edge sparse low-doublet and projector verification (20/20)](cpc-interacting-8site-sparse-doublet.md).
- [Ultra-narrow interacting doublet classification (24/24), with reference-precision cautions](cpc-interacting-ultranarrow.md).
- [Independent second-order core resolvent and three dense-LAPACK drivers: diagnosing the 1e-7 sparse-reference edge-spin discrepancy](cpc-ultranarrow-reference-diagnosis.md).
- [Exact coupling-sign unitary, induced edge polarization coefficient and four new dense plus twelve ED2 cases](cpc-weak-edge-polarization-scaling.md).
- [Dedicated-node CPU affinity, source checksum and statistical measurement protocol](cpc-controlled-node-protocol.md).
- [Working CPC manuscript prose draft](cpc-manuscript-draft.md).

**No manuscript authorship, official source tag, DOI or journal
publication is assigned by this planning document.**
