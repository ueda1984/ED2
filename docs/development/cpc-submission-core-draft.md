# CPC submission-focused manuscript core — author-review draft

**Working document — not submitted or accepted.** Prepared 15 November
2026. This is the **first proposed journal program paper for ED2**, whose
source was already publicly available, and not a post-publication ED2
update. The relevant *published software* comparator is the original
2022 QS³ code [QS-Cube/ED](https://github.com/QS-Cube/ED).
Authors/order, software version, upstream/fork integration,
program archive and CPC article category require confirmation.
This compact version intentionally separates the *main software
contribution* from exploratory numerical stress tests in the
[detailed technical draft](cpc-manuscript-draft.md).
A formally journal-formatted manuscript is not yet assembled.

## Proposed title

**ED2: Generalized exact diagonalization of mixed-spin quantum
Hamiltonians with anisotropic and nonconserving interactions**

## Abstract

We present the first journal program description of ED2, a Fortran
exact-diagonalization code extending the published QS³ solver
([QS-Cube/ED](https://github.com/QS-Cube/ED)) to Hamiltonians with unequal
on-site spin magnitudes, anisotropic exchange, directed
Dzyaloshinskii–Moriya interactions, symmetric off-diagonal
exchange and site-dependent magnetic fields. We document
the generalized basis and input conventions, dense and iterative
Hermitian eigensolvers, local and Cartesian two-spin
observables, and numerical acceptance diagnostics.
In the common XXZ limit, we directly compare the original
published QS³ executable with ED2 in three exactly
matching invariant sectors of four- and six-site spin-1/2
chains, verifying **all 32 eigenvalues** against a
separately assembled matrix. The maximum cross-code
difference is \(2.89\times10^{-15}\).
Independent dense and sparse constructions validate
fully nonconserving alternating spin-1/2–spin-1 chains
of up to 12 sites (Hilbert-space dimension 46,656).
In a reproducible two-site-periodic example, the
response of staggered vector chirality to a directed
DM coupling is checked against a ground-energy
Hellmann–Feynman derivative and independently
evaluated finite-size susceptibilities. Additional
tests verify a literature-defined XYZ–DM Hamiltonian
and assess eigenpair residuals and low-subspace
observables under near-degeneracy. We quantify
model-dependent cache storage from an analytically
verified maximum row capacity and independently
test numerical equivalence across cache policies.
A standalone, standard-library-only Python driver
reproduces the small-system physical example using
the compiled ED2 executable. The results establish
model coverage, reproducibility and numerical
validation, without proposing a new eigensolver
or claiming a novel magnetic phase.

*Editorial note:* do not add a universal speedup, dedicated-node
efficiency or DOI claim until controlled benchmarks and the
formal release are complete.

## 1. Introduction

Exact diagonalization supplies finite-size reference
information for interacting quantum-spin systems.
The exponential growth of
\(\mathcal D=\prod_i(2S_i+1)\)
restricts accessible sizes, so efficient representation
of physical Hilbert spaces and reliable treatment of
their eigenspectra are central software concerns.
Interactions with spin-orbit character often
combine unequal spin magnitudes, anisotropic
exchange, antisymmetric DM couplings and
transverse fields, and hence cannot generically
be described by a conserved total \(S^z\) block.

The published **QS³** program [1] and its
[CPC-era source `QS-Cube/ED`](https://github.com/QS-Cube/ED)
focused on large, low-down-spin-density,
spin-1/2 XXZ-type models in a fixed
\(N_\downarrow\) sector, exploiting momentum
blocks when appropriate. In addition to
Lanczos and static spin observables, the
original code included a dynamical
spin-structure-factor pathway. ED2
**is not a new Lanczos algorithm** and
does not currently show the same native
dynamic-response path: it should not
be described as an unqualified functional
superset of QS³.

**This is ED2's first proposed journal
program description**, despite the
generalized code having been publicly
accessible in [`QS-Cube/ED2`](https://github.com/QS-Cube/ED2)
before this work. Relative to the
*published* QS³, ED2's
central program advance is its
new physical operator and basis
scope: mixed on-site spin values,
anisotropic Cartesian exchange,
directed DM, symmetric Gamma terms,
site-dependent vector fields and
nonconserving full-space computations.
The later development fork also
repaired numerical-library interfaces
and field indexing and added explicit
accuracy/failure and reproducibility
controls. Those maintenance changes
must be credited, but are not the
sole new contribution of the first
ED2 software article. The
[primary published-program comparison](cpc-qs3-ed-to-ed2-scope-matrix.md)
identifies these boundaries precisely.

Other packages, including **QuSpin** [4],
can express a broad range of quantum-spin
Hamiltonians. Thus our argument is
**not exclusivity of the physical models**,
but documented use of a compact Fortran
executable with explicit physical input
records, verifiable numerical contracts,
and a controlled comparison against
independent Hamiltonian constructors.
The following sections distinguish the
supported physics, correctness evidence,
and resource behavior rather than
claiming unprecedented algorithmic complexity.

### 1.1 Relation to the previously published QS³ software

The scope difference is summarized below for the
**actual CPC-published original**
[QS-Cube/ED](https://github.com/QS-Cube/ED),
rather than the already-public but as yet
**unpublished-in-journal**
[QS-Cube/ED2](https://github.com/QS-Cube/ED2)
repository. The individual program inputs
were checked in their Fortran source;
see the [complete capability audit](cpc-qs3-ed-to-ed2-scope-matrix.md).

| Capability | Published QS³/ED | This first ED2 program |
| --- | --- | --- |
| Local Hilbert space | Fixed spin-1/2 | Site-specific \(S_i\), including mixed spins |
| Bilinear bonds | XXZ-type \(J_{xy},J_z\) | XYZ, all DM vector components, all symmetric Gamma components |
| Site fields | No general three-axis site-record Hamiltonian input | Individually specified Cartesian \(h_i^{x,y,z}\) |
| Excitation sectors | Fixed \(N_\downarrow\), optionally momentum-reduced | Full/product space, or local/global occupation windows and compatible user-specified symmetry |
| Static observables | Longitudinal/transverse channel correlations | All local Cartesian moments and nine two-site Cartesian components |
| Dynamic structure factor | Original DSFLan route available | Equivalent native ED2 route not demonstrated |
| Eigensolvers/OpenMP | Already implemented | Mathematical solver families inherited and numerically revised, not invented |

We separately verify the **same XXZ Hamiltonian**
in both executable codes for three invariant
magnetization sectors. All **32 eigenvalues**
of the combined \(L=4\) and \(L=6\) tests
match an independent Hamiltonian;
the maximum original-to-ED2 discrepancy is
\(2.89\times10^{-15}\).
This tests backward *physical-model*
continuity, not binary/file compatibility
or legacy dynamical features.

## 2. Physical model and input conventions

ED2 represents a lattice of site-dependent
spin magnitudes \(S_i\) and input-defined
oriented bonds \(i,j\), with a general
two-body bilinear Hamiltonian

\[
H=\sum_{(i,j)}\left[
\sum_{\alpha=x,y,z}J_{ij}^\alpha S_i^\alpha S_j^\alpha
+\boldsymbol D_{ij}\cdot(\boldsymbol S_i\times\boldsymbol S_j)
+\sum_\alpha\Gamma_{ij}^\alpha
(S_i^\beta S_j^\gamma+S_i^\gamma S_j^\beta)
\right]+\sum_i\boldsymbol h_i\cdot\boldsymbol S_i ,
\]

where \((\alpha,\beta,\gamma)\) are the Cartesian
cyclic permutations. ED2's physical input
convention is \(+\boldsymbol h\cdot\boldsymbol S\);
to encode the usual Zeeman energy
\(-\boldsymbol B\cdot\boldsymbol S\),
the user supplies \(\boldsymbol h=-\boldsymbol B\).
The sign of a DM coefficient depends on
the orientation of the directed input bond.

The full local product Hilbert space has
dimension \(\mathcal D=\prod_i(2S_i+1)\).
Selected basis representations may exploit
compatible symmetries or restrict local/total
occupation ranges. **A restricted excitation
window in a non-U(1) model is generally a
projection \(PHP\), not an exact invariant
sector**. Symmetry operations are user-provided;
no unimplemented universal
commutation-checking algorithm is claimed.
The validated model instances below use
the full product space, identity spatial
sector and no fixed-\(S^z\) restriction.

## 3. Numerical implementation and acceptance

The executable provides dense Hermitian
LAPACK diagonalization and two iterative
modes: a reconstructed-vector Lanczos
path and restarted Hermitian iteration.
The Fortran code relies on BLAS/LAPACK
and optional OpenMP on a shared-memory
node; it is **not** a distributed MPI
program or a stable installed library API.

Source maintenance on the development
line covers explicit LP64 BLAS/LAPACK
interfaces, error-path checks, an indexing
repair for sparse/site-reordered fields,
Lanczos breakdown handling, reproducible
random initialization and opt-in
direct-residual policies. These are
engineering correctness changes, not
new eigenvalue algorithms. Their
exact upstream/source comparison is
recorded in the [contribution audit](cpc-software-provenance-contribution-audit.md).

For a returned Ritz pair
\((E_j,\psi_j)\), define its
direct residual

\[
r_j=\|H\psi_j-E_j\psi_j\|_2 .
\]

Small \(r_j\) tests the eigenvalue
equation but does **not** itself
prove that the globally lowest
state was obtained. In a near
degenerate spectrum, it also
does not guarantee accurate
state-specific observables.
We therefore test direct residuals,
spectrum inclusion against independent
references and physical observables
as **separate conditions**.
Numerical failures and output
contract violations are reported
rather than silently converted
into successful runs.

## 4. Independent validation strategy

To verify continuity with the **actual published**
QS³ implementation, we compiled the original
`QS-Cube/ED` source at the pinned CPC-era SHA
`a3b09934bba5f966bfb564718528de8778c2b46c`
and the ED2 source distribution at SHA
`88c72032c26e52e507bc1217ec030479ada486fc`.
On a four-site open spin-1/2 XXZ chain
with \(J_{xy}=0.8\), \(J_z=1.1\) and
\(N_\downarrow=2\), both executable
dense solvers represent the same
six-dimensional invariant sector.
All six eigenvalues coincide with a
third independently assembled
product-state Hamiltonian. The largest
old-to-new absolute difference is
\(2.89\times10^{-15}\).
The expanded source-pinned test
further verifies \(L=6,N_\downarrow=1\) and
\(L=6,N_\downarrow=3\) at the same physical
XXZ couplings. Together, the three invariant
sectors contain **32 independently verified
eigenvalues**, with maximum original-to-ED2
difference \(2.89\times10^{-15}\) and maximum
ED2-to-independent difference
\(3.44\times10^{-15}\).
The [32-row frozen spectrum](cpc-qs3-ed-ed2-32level-extended-ci.csv)
and the [full original/modern source comparison](cpc-qs3-ed-to-ed2-scope-matrix.md)
document both the original
near-saturation regime and the
more populated conserved sector.

This is a direct continuity test for shared
physics, not a demonstration of new
solver performance, fully backward-compatible
input format or the old dynamical
structure-factor pathway.

Small-system dense Hamiltonians
and independent sparse product-basis
constructors validate energies,
symmetry-compatible representations,
field indexing, and selected wavefunction
observables. The verification covers
the three diagonal exchange terms,
directed DM, symmetric exchange and
site-dependent fields, without
assuming total \(S^z\) conservation.

A low-energy doublet stress model
deliberately weakens a terminal
anisotropic XY bond while retaining
the fully interacting non-U(1)
mixed-spin core. An independent
1296-dimensional sparse problem
computes its first three levels.
Native ED2 first- and second-state
wavefunctions are re-read through
the observable pipeline and their
two-state ensemble compared with the
independent spectral projector.
This detects numerical mistakes
not visible from energy values alone.
A separate dense-LAPACK and
core-projected effective-model
cross-check explains why individual
eigenvectors can be less stable than
the isolated **two-dimensional
low-energy projector**.
These are adversarial numerical
tests; their modified edge Hamiltonian
is not the main connected physical
example and does not imply
a bulk phase transition.

## 5. Scientific application: staggered DM response

Consider an OBC chain with alternating
spin-1/2 and spin-1 sites.
The on-site fields and nearest-neighbor
\(J_x,J_y,J_z,\Gamma_z\) couplings
obey the same two-site-periodic
rules for every system length
\(L=4,6,8,10,12\). The oriented
DM term is

\[
H_{\mathrm{DM}}=D\sum_{i=1}^{L-1}
(-1)^{i-1}
(S_i^xS_{i+1}^y-S_i^yS_{i+1}^x).
\]

This preserves nonzero transverse
fields, XYZ anisotropy and symmetric
off-diagonal exchange; it is *not*
simplified to an XXZ/U(1) model.
The coefficient set is explicitly
defined in the
[worked physical example](cpc-stationary-dm-application.md),
rather than fitted to a specific material.

Defining the conjugate
staggered chirality operator \(K\),
the nondegenerate ground state obeys

\[
\frac{\partial E_0}{\partial D}=\langle K\rangle,
\quad
\chi_\kappa(L)=-\frac{1}{L-1}
\left.\frac{\partial^2E_0}{\partial D^2}\right|_{D=0}.
\]

For \(4\le L\le12\), **60/60**
source-archived ED2 trials
(OMP=1/2) pass independent
energy and full Cartesian
nearest-neighbor correlation
checks. The vector chirality
obtained from ED2's unusual
`xx,yx,zx,xy,yy,zy,xz,yz,zz`
output order agrees with the
independent energy derivative.
The susceptibility is independently
checked through energy curvature
and a derivative of the chiral
observable. For 4 and 6 sites
it is additionally compared to
a full spectral sum.

The calculated
\(\chi_\kappa(0;L)\) is approximately
0.6343, 0.5687, 0.5374, 0.5197
and 0.5085 for \(L=4,6,8,10,12\),
respectively. These are **finite
sizes only**; no thermodynamic
extrapolation or transition
is claimed.

As a distinct limiting test,
the program reproduces the
properly converted finite-size
Hamiltonian of a published
XYZ–DM spin-1/2 model [3],
including boundary orientation.
This does not reproduce the
published phase diagram.

## 6. Storage and runtime characterization

The existing program controls
`wk_dim` (number of cached
basis rows) and `MNTE` (cached-row
entry capacity) provide a
model-dependent memory–work
tradeoff. For the selected
fully interacting alternating
OBC Hamiltonian, counting the
allowed single-spin transverse
and two-site transition channels
gives an exact per-row capacity

\[
\boxed{\mathrm{MNTE}_{\rm model}=7L/2-1}.
\]

This bound was independently
verified against sparse
row enumerations for tested
system sizes. On the current
GNU LP64 layout, the two
cache work arrays require
20 bytes per stored entry.
For \(L=12\), caching all
\(\mathcal D=46{,}656\) rows
at this tight capacity consumes
approximately **37.38 MiB**,
less than the approximately
**50.28 MiB** allocated for
a default-width cache of
half the rows. These values
describe **cache arrays**,
not total process RSS.

The [108-case same-CI-job
resource study](cpc-cache-pareto-evaluation.md)
confirms numerical equivalence
across three cache policies
and records process peak RSS,
wall time and Hamiltonian
application counts for matched
physical targets. Although
the timing differences are
large for some configurations,
these measurements used
shared hosted runners, with
only two measured repetitions.
They are **not** presented as
hardware-independent speedups,
strong-scaling results or
proof of a universally optimal
cache setting.

A CPU-pinned, source-/input-hashed
[dedicated-node protocol](cpc-controlled-node-protocol.md)
and strict [postprocessor](cpc-controlled-node-protocol.md#d-relevance-to-the-cpc-paper)
are available; fully repeated,
operator-verified exclusive-node
performance measurements are
still pending.

## 7. Reproduction and software availability

The scientific application
requires only the compiled
ED2 executable, its
Fortran BLAS/LAPACK runtime,
and input records. We supply
a [portable quickstart](cpc-portable-dm-quickstart.md)
that builds the complete
\(L=4,6\) physical models,
runs \(D=-2\times10^{-4},0,
+2\times10^{-4}\), and checks
the energy derivative, chirality,
susceptibility and frozen
independent reference.
The Python driver uses the
**standard library only**;
NumPy/SciPy are confined to
the optional independent
verification tools.

ED2 has not previously been reported in a
journal program article, although its
source was public before the present work.
The earlier published QS³ program and its
original [`QS-Cube/ED` source](https://github.com/QS-Cube/ED)
are the paper's principal software lineage.
The [direct feature and XXZ-spectrum comparison](cpc-qs3-ed-to-ed2-scope-matrix.md)
is part of this first ED2 submission.

The current development work
is publicly visible in the
[ED2 fork](https://github.com/ueda1984/ED2/tree/codex/ed2-cpc-stationary-dm-scaling).
The [upstream source repository](https://github.com/QS-Cube/ED2)
and MIT notice remain the
underlying program lineage.
The formal CPC source archive
version, repository integration,
citation metadata, Mendeley
Data accession and DOI **have
not yet been agreed or assigned**.
The program is a CLI executable;
no stable Fortran/Python library
ABI or MPI support is asserted.

## 8. Conclusions and outstanding limitations

ED2 provides a documented
generalized quantum-spin
exact-diagonalization executable
beyond the model scope of the
published QS³ near-saturation
program. The source-maintenance
and independent physical checks
make its use for mixed-spin,
nonconserving Hamiltonians
reproducible and expose solver
failures and spectral limitations
that can otherwise be
misinterpreted. The staggered
DM response illustrates a
scientifically interpretable
observable without claiming
new chiral physics.

Before submission, the
maintainers must confirm
authorship, code-version metadata,
official source archive,
program accessibility and
the manuscript category.
Any publication-quality
timing comparison should come
from repeatedly verified,
exclusive compute-node
measurements and should
retain an explicit physical
workload and uncertainty.
The limitations of exponentially
growing Hilbert spaces,
near-degenerate ground-state
identification and the absence
of a stable installed library
interface must remain explicit.

## References (provisional)

1. H. Ueda, S. Yunoki and T. Shimokawa, QS³ near-saturation
   spin solver, *Computer Physics Communications* **277**,
   108369 (2022).
   https://doi.org/10.1016/j.cpc.2022.108369
2. A. Lahiri and S. K. Pati,
   mixed-spin chain with DM couplings,
   *Journal of Physics: Condensed Matter* (2022).
   https://doi.org/10.1088/1361-648X/ac4aab
3. Yi *et al.*, XYZ–DM spin-1/2 chain,
   *Physical Review B* **100**, 024423 (2019).
   https://doi.org/10.1103/PhysRevB.100.024423
4. P. Weinberg and M. Bukov, QuSpin numerical
   exact-diagonalization framework.
   [Reference metadata to be verified for final bibliography.]
