# CPC manuscript: scientific first draft (for author review)

> **Working draft, 15 November 2026.** Not submitted, accepted or associated
> with a software DOI. Author list, article category, program summary and
> release/version information are to be approved. The main 60/60 DM response
> results and the auxiliary 12/12 literature-Hamiltonian plus 36/36
> resource-characterization trials, 52/52 cache-policy trials,
> 36/36 two-eigenpair trials and 72/72 MNTE-capacity trials have
> archived independent validation.
> Repeated controlled-node performance and a public installed library
> interface have **not** been established.

## Provisional title

**ED2: Exact diagonalization of generalized quantum spin Hamiltonians with mixed local spins and anisotropic nonconserving interactions**

## Abstract

We present ED2, a Fortran exact-diagonalization program for lattice
quantum-spin Hamiltonians with site-dependent spin magnitudes,
anisotropic exchange, Dzyaloshinskii–Moriya interactions, symmetric
off-diagonal exchange and magnetic fields. Extending the model
coverage of the published near-saturation QS³ solver, ED2 provides
full Hilbert-space calculations, user-selected symmetry representations
and optional occupation-number restrictions. We distinguish exact
symmetry reductions from projections that need not commute with
the Hamiltonian. Dense and iterative eigensolvers, local magnetic
moments and Cartesian two-spin correlations are verified using
independently constructed Hamiltonians, with numerical diagnostics
that distinguish convergence from ground-state inclusion. As a
reproducible scientific application we compute the response to a
staggered directed DM interaction in an alternating spin-1/2 and
spin-1 chain containing XYZ, symmetric off-diagonal and transverse
field terms. For 4–12 sites, ground energies and all nearest-neighbor
correlations agree with independent references, and the DM-conjugate
vector chirality satisfies Hellmann–Feynman and finite-size
susceptibility identities. A separate test reproduces the finite-size
Hamiltonian of a published XYZ-DM spin-1/2 chain, including periodic
bond orientation and Pauli-to-spin conversion. We characterize the
memory and Hamiltonian-application cost of matched one-eigenpair
calculations, quantify cache storage against source-level allocation,
and independently verify two-lowest-eigenpair restart workspaces.
For the selected open alternating-spin model, a model-specific
exact cached-row capacity of \(7L/2-1\) is derived and verified,
reducing allocated memory without changing the
Hamiltonian-application count. Definitive timing and parallel-speed
comparisons are reserved for controlled-node measurements.

**Abstract revision required before submission:** replace the
general last sentence with a quantitative controlled-node result
once it is available. Avoid an unsupported claim of performance
superiority, universal accuracy or new chiral phases.

## 1. Introduction

Exact diagonalization remains a reference method for studying finite
quantum many-body systems. In quantum-spin models, exploiting exact
conservation laws and spatial symmetries can reduce the working
Hilbert space, but interactions with spin-orbit character, unequal
local spin magnitudes and transverse magnetic fields challenge
special-purpose sector decompositions. These features are especially
relevant when the desired outputs are physical observables sensitive
to the wave function, not only ground-state energies.

The QS³ program, published in *Computer Physics Communications*
in 2022 [1], addresses near-saturation quantum-spin calculations
and already provides iterative eigensolvers, observables and OpenMP
parallelism. ED2 must therefore be positioned as an extension
in the supported local Hilbert spaces and interaction structures,
rather than as the invention of Lanczos, restarted Krylov iteration
or parallel exact diagonalization. In particular, ED2's generalized
exchange and field inputs can represent Hamiltonians in which
total \(S^z\) is not conserved.

Other exact-diagonalization frameworks, including QuSpin [4],
already offer general spin-model Hamiltonians and
user-programmable symmetry representations. Therefore
our claim is **not** that nonconserving XYZ or DM models
are inaccessible to all existing software. Rather,
we seek to document and independently validate the
specific native mixed-local-spin interaction inputs,
numerical acceptance contracts and resource tradeoffs
of ED2 relative to its previously published QS³ lineage.

A challenge in assessing generalized spin software is to establish
that independent definitions of the Hamiltonian, symmetry,
spectrum and observables agree. A directed
Dzyaloshinskii–Moriya (DM) interaction supplies an exact and
sensitive diagnostic: differentiation of the ground-state energy
with respect to the DM coefficient gives a vector-chirality
expectation value. This identity tests both an oriented
spin operator and the component order of the program's
two-body correlation output.

The spin-chain literature provides relevant physical context.
Lahiri and Pati [2] discuss alternating spin-1/2 and spin-1
chains with competing exchanges and DM couplings; Yi and
colleagues [3] examine XYZ anisotropy, magnetic field and DM
exchange for spin-1/2 chains. Our main mixed-spin demonstrator
uses newly selected, size-stationary coefficients and is not a
reproduction of their phases or data. Separately, we verify
ED2 against the precisely defined periodic XYZ-DM Hamiltonian
of Ref. [3]. The work is thus a software-model and
observable-validation study, not a claim of new chiral physics.

## 2. Supported models and basis representation

ED2 accepts local spin values \(S_i\) and two-site records with
Cartesian exchange \(J_{ij}^{x,y,z}\), DM exchange
\(\boldsymbol D_{ij}\), and symmetric off-diagonal
\(\Gamma^\mu_{ij}\) terms. One representative Hamiltonian is

\[
H=\sum_{(i,j)}\Big[
\sum_{\mu=x,y,z}J_{ij}^\mu S_i^\mu S_j^\mu
+\boldsymbol D_{ij}\cdot(\boldsymbol S_i\times\boldsymbol S_j)
+\Gamma^x_{ij}(S_i^yS_j^z+S_i^zS_j^y)
+\Gamma^y_{ij}(S_i^zS_j^x+S_i^xS_j^z)
+\Gamma^z_{ij}(S_i^xS_j^y+S_i^yS_j^x)
\Big]+\sum_i\boldsymbol h_i\cdot\boldsymbol S_i.
\]

The input field term is **positive**: the conventional
physical Zeeman term \(-\boldsymbol B_i\cdot\boldsymbol S_i\)
requires \(\boldsymbol h_i=-\boldsymbol B_i\).
The DM term requires a directed pair; reversing its
orientation also reverses the oriented DM coefficient.

With local spin \(S_i\), the complete product-space dimension is
\(\mathcal D=\prod_i(2S_i+1)\). Occupation-number limits
restrict the represented basis, but such a restriction need
not be a symmetry. For a Hamiltonian without the required
\(U(1)\) invariance, it defines \(PHP\), generally an
approximation to the unrestricted problem rather than a true
invariant sector. This distinction is explicit in the
validation and user-facing numerical documentation.
Likewise, arbitrary site-permutation inputs do not establish
that the Hamiltonian commutes with the requested operation:
compatible symmetry representations require checking the
full set of local and bond coefficients.

## 3. Implementation and numerical acceptance

The existing executable provides reconstructed-vector
Lanczos (ALG=1), restarted Hermitian iteration (ALG=2),
and dense LAPACK diagonalization (ALG=3). Fortran
source-distribution validation, residual reporting,
explicit starting vectors, deterministic-seed controls and
diagnostic sidecars support reproducible scientific
verification. These engineering features are not presented
as distinct new eigensolver theories.

For an approximate eigenpair \((E_n,\psi_n)\),
the direct residual \(r_n=\|H\psi_n-E_n\psi_n\|\)
measures satisfaction of the eigenvalue equation.
A small \(r_n\) does not alone certify that the global
ground state has been found, and it does not ensure
accurate individual observables when the low-energy
spectrum is nearly degenerate. Therefore the acceptance
tests compare the target eigenvalues and physical
observables with independently assembled dense or
sparse reference Hamiltonians. Actual nonconvergence
and invalid solver-output states are recorded separately
from successful numerical validations.

### 3.1 Near-degenerate spectral subspaces

We further test a deliberately adversarial eight-site Hamiltonian
formed by detaching the first spin-1/2 from a seven-site,
non-\(U(1)\), XYZ/Gamma/DM interacting core. The detached spin
is subjected to a small longitudinal bias \(\delta\). This
**stress Hamiltonian is not** the connected size-stationary
DM-response application of Sec. 4. It admits the independently
checkable tensor structure
\[
H_{\rm stress}=I_2\otimes H_{\rm core}
+\delta S^z_{\rm edge}\otimes I_{\rm core}.
\]
The first two levels are separated by \(|\delta|\) for
\(\delta=0,10^{-10},10^{-6}\), while their separation
from the third level remains approximately 0.339. An
independent sparse calculation constructs the two-dimensional
projector \(P_{\mathcal C}\) using the active-core ground
vector and verifies the tensor factorization. Within a
doublet, an individual eigenvector may change under an
arbitrary unitary rotation even if its direct residual
is tiny; for an observable \(O\) the invariant comparison is
\[
\langle O\rangle_{\mathcal C}=\tfrac12
\operatorname{Tr}(P_{\mathcal C}O).
\]
We independently read the two ED2 native Fortran
wavefunctions and evaluate the local moments and all nine
Cartesian correlations for all 28 site pairs. All **12/12**
OMP=1/2 and two-workspace doublet trials passed; the
largest trace-averaged correlation or magnetization
error was below \(3.7\times10^{-10}\). At exact degeneracy,
the first-state edge magnetization ranged, for example,
from approximately \(-0.025\) to \(-0.400\) between two
OMP runs of NOK/NOM=2/16, while the two-state mean
remained zero. This verifies **subspace-invariant
observables**; it neither makes each state uniquely
physical nor certifies global ground-state inclusion
from a small residual. The complete setup and outputs
are archived in the [eight-site doublet report](cpc-stationary-sparse-doublet.md).

### 3.2 Interacting near-degenerate eight-site validation

The preceding rank-two test admits an exactly
factorized detached-edge Hamiltonian, so it is not
sufficient evidence that low-subspace observables
remain accurate for *interacting* near-degenerate
states. We therefore retain the same seven-site
non-\(U(1)\) XYZ/Gamma/staggered-DM core but couple
its edge site to an additional spin-1/2 through

\[
H_{\rm edge}=\delta S_1^z+
\varepsilon(S_1^x S_2^x+0.6 S_1^yS_2^y).
\]

This is explicitly a modified *stress Hamiltonian*,
not the connected, uniform-rule physical DM
demonstrator of Sec. 4. Independent SciPy sparse
diagonalization of the full 1296-dimensional
matrix computes the three lowest energies,
their direct eigenpair residuals, edge-spin
entanglement and the two-dimensional low
spectral projector. The XY bond generates
a genuinely entangled ground state; at
\(\varepsilon=0.02,\delta=0\), its edge
linear entropy is approximately
\(2.19\times10^{-4}\).

Five parameter combinations yield a first
gap between \(1.07\times10^{-4}\) and
\(4.84\times10^{-3}\), while the external
separation \(E_2-E_1\) stays between
approximately 0.333 and 0.339.
For each Hamiltonian the two ED2 low
native eigenvectors are verified with
NOK/NOM=2/16 and 4/24 at OMP=1/2.
Both vectors are passed through the actual
ED2 observable output path, using
native-basis second-state readback for
the latter. The equal-weight local
magnetizations (24 components) and
Cartesian correlations (252 components)
are compared against the independently
constructed interacting rank-two projector.
**20/20 tests passed**, with the largest
trace-observable difference
\(3.7433\times10^{-10}\) and maximum
direct residual \(2.676\times10^{-9}\)
under a \(10^{-8}\) request.

This supplies an interacting complement
to the detached-edge doublet experiment,
rather than relying on a tensor-product
identity alone. A residual-controlled
Ritz pair does not, by itself, certify
global lowest-state inclusion. If
the exact external third level \(E_2\)
is independently known, the low-doublet
projector admits a conditional leakage
estimate:
for orthonormal Ritz vectors \(u_j\)
with \(\lambda_j<E_2\) and
\(r_j=\|(H-\lambda_j)u_j\|\),

\[
\|(I-P_{\mathcal C})u_j\|
\le\frac{r_j}{E_2-\lambda_j}.
\]

It follows from the spectral decomposition
and the rank-two projector norm inequalities
that for \(\|O\|_{\rm op}\le1\),

\[
\left|\frac12\operatorname{Tr}[(Q-P_{\mathcal C})O]\right|
\le\sqrt2\left[
\sum_{j=0}^1
\frac{r_j^2}{(E_2-\lambda_j)^2}
\right]^{1/2},
\]

where \(Q\) is the returned orthonormal
rank-two projector. This bound uses the
**external** gap, not the small inner
doublet gap, and does not constitute
a self-contained ED2 ground-state
certificate. Finite floating-point
orthogonality and reference-level
accuracy must be controlled separately.
The [interacting sparse validation report](cpc-interacting-8site-sparse-doublet.md)
and [frozen 20-case source data](cpc-interacting-8site-raw-ci.csv)
document the observed results and
qualification.

## 4. Scientific application: staggered DM vector chirality

Consider an open chain with
\(S_i=1/2,1,1/2,1,\ldots\), with zero-based
site and bond-start indices and \(\eta_i=(-1)^i\).
We examine

\[
\begin{aligned}
H(D)=&\sum_{i=0}^{L-2}\big[
J_i^x S_i^xS_{i+1}^x+
J_i^y S_i^yS_{i+1}^y+
J_i^z S_i^zS_{i+1}^z\\
&\quad+\Gamma_i^z(S_i^xS_{i+1}^y+S_i^yS_{i+1}^x)
+D\eta_i(S_i^xS_{i+1}^y-S_i^yS_{i+1}^x)\big]\\
&+\sum_{i=0}^{L-1}(h_i^x S_i^x+h_i^z S_i^z).
\end{aligned}
\]

The local coefficients are fixed by a two-site rule,
independently of \(L\):

\[
\begin{aligned}
J_i^x&=0.77+0.05\eta_i, & J_i^y&=0.63+0.03\eta_i,\\
J_i^z&=1.09-0.04\eta_i, & \Gamma_i^z&=0.065+0.015\eta_i,\\
h_i^x&=0.15(0.4\eta_i+0.1), &
h_i^z&=0.17+0.035\eta_i .
\end{aligned}
\]

The presence of unequal XYZ exchange, nonzero
\(\Gamma^z\) and transverse fields makes this
a genuine non-\(U(1)\) calculation. The full
Hilbert space has dimension \(6^{L/2}\).
The parameters specify a reproducible illustrative
test model, not a material Hamiltonian inferred
from experiment.

Defining
\[
\widehat K=\sum_{i=0}^{L-2}\eta_i
(S_i^xS_{i+1}^y-S_i^yS_{i+1}^x),
\]
the Hellmann–Feynman identity for a differentiable
nondegenerate ground state is
\[
\frac{\partial E_0}{\partial D}
=\langle\widehat K\rangle,\qquad
\bar\kappa_{\rm stag}=\langle\widehat K\rangle/(L-1).
\]
The corresponding finite-size susceptibility is
\[
\chi_\kappa(L)=
-\frac{1}{L-1}\left.\frac{\partial^2E_0}{\partial D^2}\right|_{D=0}
=-\frac{1}{L-1}\left.
\frac{\partial\langle\widehat K\rangle}{\partial D}\right|_{D=0}.
\]
For the finite nondegenerate ground state it can
also be written as a spectral sum,
\[
\chi_\kappa(L)=\frac{2}{L-1}\sum_{m>0}
\frac{|\langle m|\widehat K|0\rangle|^2}{E_m-E_0}.
\]
This third route is evaluated for \(L=4,6\).
The ground energy, as the minimum of affine
functions of \(D\), is concave, so the conjugate
chirality must be nonincreasing with \(D\) where
the derivative exists.

### 4.1 Results and interpretation

| L | Full dimension | Ground-to-first gap at D=0 | Chiral mean at D=0 | Chiral susceptibility |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 36 | 0.17235306 | +0.000383416 | 0.63431120 |
| 6 | 216 | 0.22526115 | +0.000218011 | 0.56873293 |
| 8 | 1,296 | 0.28340937 | +0.000156291 | 0.53735526 |
| 10 | 7,776 | 0.32938732 | +0.000133476 | 0.51970132 |
| 12 | 46,656 | 0.35158094 | +0.000123044 | 0.50850099 |

For \(D=-0.24,-0.12,0,+0.12,+0.24\),
the ED2 low spectrum, local magnetization and
all nine Cartesian nearest-neighbor correlations
are checked against an independent sparse
Hamiltonian; dense references are also used at
four and six sites. All 60 trials (five lengths,
two OMP settings, five DM values plus a
two-eigenpair check at \(D=0\)) passed,
with the greatest eigenenergy error below
\(1.4\times10^{-13}\) and the total staggered
chirality difference below \(1.5\times10^{-9}\).
The energy finite-difference derivative agreed
with the independent chirality expectation
to within \(2.5\times10^{-10}\) at \(D=0\).

The unusual Cartesian output order after the
two site labels is
\((xx,yx,zx,xy,yy,zy,xz,yz,zz)\).
Explicit verification of both oriented
off-diagonal components is therefore crucial.
The mean chirality varies smoothly along the
tested finite-size parameter cut, and the
susceptibility decreases as the chain length
grows through twelve sites. Neither observation
establishes a thermodynamic extrapolation,
spontaneous order or a phase transition. In
particular, small nonzero \(D=0\) chirality
is compatible with the explicit anisotropy,
site parity and finite boundaries of this model.

## 5. Independent comparison with a published XYZ-DM Hamiltonian

As a complementary model-level test we implement
the periodic spin-1/2 Hamiltonian of Yi *et al.*
[3], with Pauli matrices \(\sigma^\mu\):
\[
H_{\rm Yi}=J\sum_i\left[
\frac{1+\gamma}{2}\sigma_i^x\sigma_{i+1}^x+
\frac{1-\gamma}{2}\sigma_i^y\sigma_{i+1}^y+
\Delta\sigma_i^z\sigma_{i+1}^z\right]
+D\sum_i(\sigma_i^x\sigma_{i+1}^y-\sigma_i^y\sigma_{i+1}^x)
-h\sum_i\sigma_i^z.
\]
The exact conversion to ED2's \(S=\sigma/2\)
and positive input-field conventions is
\[
J_x=2J(1+\gamma),\quad
J_y=2J(1-\gamma),\quad
J_z=4J\Delta,\quad D_z=4D,\quad h_z=-2h .
\]
For independently chosen control parameters
\(J=1,\gamma=0.37,\Delta=0.41,h=0.58\),
\(D=0,0.19\), and \(L=4,6,8\), the
independently assembled Pauli matrix gives
the full spectrum and physically oriented
two-site correlations. All 12 ED2 comparisons,
including both OMP=1/2, passed. Maximum
full-spectrum differences were below
\(3\times10^{-14}\), while magnetization
and correlation differences were below
\(10^{-14}\). The periodic closing bond
tests the directed DM sign and the output
correlation-order convention.

This is a reproduction of a **published
Hamiltonian definition** at selected finite
sizes, not a numerical reproduction of its
reported bulk phase diagram or Wigner–Yanase
coherence quantities. The main scientific
application remains the mixed-spin
nonconserving model in Sec. 4.

## 6. Matched-target resource characterization

A comparison of two algorithms can be misleading
when one requests one state while the other
requests two, or when their numerical acceptance
criteria differ. We therefore request one
lowest eigenpair from both ALG=1 and ALG=2,
with an identical explicit complex initial
vector, the same fixed RNG seed, full
representation, physical observables and
absolute residual requirement \(10^{-8}\).
The ALG=2 restart dimensions (NOK,NOM) are
recorded as separate choices. We use
one uncounted warm-up and two interleaved
measurement repetitions for each condition.

On the stationary alternating-spin model
at \(D=0.12\), the following data were
measured on Ubuntu 24.04 GitHub-hosted
CI with GNU Fortran 13, -O2 -fopenmp and
single-thread BLAS. The full 36-run
acceptance suite passed independently
for both OMP=1 and OMP=2.

| L | OMP | Solver | Median H applications | Median wall time (s) | Median RSS (KiB) |
| ---: | ---: | --- | ---: | ---: | ---: |
| 8 | 1 | ALG=1 | 120 | 0.030616 | 8,454 |
| 8 | 1 | ALG=2 NOK=2,NOM=16 | 306 | 0.075426 | 10,160 |
| 8 | 1 | ALG=2 NOK=4,NOM=24 | 475 | 0.123207 | 10,716 |
| 10 | 1 | ALG=1 | 200 | 0.285188 | 22,128 |
| 10 | 1 | ALG=2 NOK=2,NOM=16 | 340 | 0.547111 | 27,850 |
| 10 | 1 | ALG=2 NOK=4,NOM=24 | 550 | 0.930786 | 30,800 |

The ALG=2 workspace strongly influences the
observed resource cost; the listed ALG=1
measurements use fewer Hamiltonian applications
for this one-ground-state task. We draw no
general algorithmic ranking or strong-scaling
conclusion from these shared-runner,
two-repeat measurements. The performance
section must be augmented with repeated
isolated-node timings and pinned compiler/BLAS,
CPU/cache and affinity records before submission.
A separate cache-policy scan and two-lowest-state
workspace study have now been independently
verified and are summarized below.

### 6.1 Cache policy and the memory-time tradeoff

The ED2 Hamiltonian-action kernel can store matrix elements
for a chosen prefix of representative states. Holding the
physical model, eigenpair target, explicit starting vector
and requested residual fixed, we compare no cached rows
(`wk_dim=0`), half of the rows and all rows
(`wk_dim=-1`). At \(D=0.12\), the 8- and 10-site
full-space ALG=1 runs have the same H-application
counts across cache settings, providing a controlled
test of caching rather than one that conflates
performance with the number of Krylov iterations.

| L | OMP | Cache | H applications | Wall time (s) | Peak RSS (KiB) |
| ---: | ---: | --- | ---: | ---: | ---: |
| 8 | 1 | none | 120 | 0.529700 | 6,646 |
| 8 | 1 | half | 120 | 0.278283 | 7,700 |
| 8 | 1 | full | 120 | 0.024038 | 8,604 |
| 10 | 1 | none | 200 | 7.511783 | 8,082 |
| 10 | 1 | half | 200 | 3.916686 | 15,218 |
| 10 | 1 | full | 200 | 0.230234 | 22,278 |

All 36 small-system cache trials across OMP=1/2
passed independent lowest-energy, ground-observable
and direct-residual checks, including the warm-ups.
At L=12 we compare half and full caching for
both ALG=1 and ALG=2 (NOK=2,NOM=16), again
at OMP=1 and 2. All 16 additional trials
also passed the physical acceptance criteria.
For example, at L=12 and OMP=1, ALG=1 takes
22.146 s with half caching versus 1.954 s
with full caching, with median peak RSS of
67,980 versus 119,468 KiB. These ratios
are **observations on shared CI hosts**, not
portably reproducible speedups.

The expected cache memory follows directly
from the two Fortran work arrays of shape
\((M,\mathrm{wk\_dim})\), with
\[
M=2\,\mathrm{NO\_one}+8\,\mathrm{NO\_two}+1.
\]
For the open chain this is \(M=10L-7\).
With a 4-byte default integer and 16-byte
`complex(8)` on the tested compiler, the two
work arrays have nominal allocation
\[
B_{\rm cache}=20(10L-7)\,\mathrm{wk\_dim}
\quad\text{bytes}.
\]
At L=12, full cache alone requires
100.558 MiB, compared with 50.279 MiB
for half cache. The observed OMP=1 ALG=1
median peak RSS difference, 51,488 KiB,
is in close agreement with this expected
allocation difference. The overall RSS
includes additional solver vectors, basis
tables and runtime overhead.

The combined **52/52 cache-policy validations**
therefore provide both independent physical
correctness and a mechanistic accounting
of memory cost. The full raw CI measurements,
including warm-ups, are preserved in
[the cache CSV](cpc-cache-sweep-raw-ci.csv)
and interpreted in the
[resource report](cpc-cache-multipair-evaluation.md).
The cache fraction refers to a prefix of
representatives rather than a random subset,
so its intermediate performance is not
necessarily a universal linear function.

### 6.2 Two-lowest-state restart robustness

A restart policy that is efficient for one
ground state is not automatically appropriate
for multiple eigenpairs. We separately
request the two lowest eigenpairs with
ALG=2, an identical explicit ED2-basis
starting vector, fixed seed, full sector
and the same residual threshold.
For L=8 and 10 at D=0 and 0.12,
two repetitions were made for OMP=1/2.
L=12 at D=0.12 was also examined at
both thread counts. The **36/36 trials**
passed independent checks of both lowest
energies, both reported direct residuals,
and ground-state local/NN observables.

At L=12, D=0.12, OMP=1, the
NOK/NOM=2/16 restart uses 540 H
applications with its largest direct
residual \(4.19\times10^{-9}\),
whereas 4/24 uses 754 H applications
with a residual \(2.52\times10^{-11}\).
Both satisfy the \(10^{-8}\) requirement.
At L=10 the smaller workspace produces
a residual as high as \(8.62\times10^{-9}\),
close to the acceptance boundary.
Consequently, the lowest H-application
count need not imply a desirable
numerical safety margin across more
demanding spectra. The
[two-eigenpair source CSV](cpc-multipair-raw-ci.csv)
preserves per-run measured costs and
residuals. These finite-gap data
do **not** validate arbitrary near-degenerate
multiplets or NOE>2.

**Figure planned:** a two-part performance
figure showing wall time and RSS versus
cached fraction with markers for L and
OMP, and a separate restart-workspace
cost-versus-achieved-residual plot for
two lowest eigenpairs. Plot data and
captions must distinguish the shared
CI runner from future controlled-node
timing evidence. No universal
performance comparison with QS³ follows
from these observations.

### 6.3 Exact cache-row capacity for the selected Hamiltonian

The preceding cache study varies the fraction of basis states
stored in memory. A complementary option is to reduce the
**row width** of the cache itself while retaining *all*
basis states in the cached representation. The existing ED2
`MNTE` input specifies a conservative number of matrix
elements per representative state. With `NO_one=L` and
`NO_two=L-1`, the generic default `MNTE=-1` expands to a
capacity of \(10L-8\), which exceeds the transitions
actually present in the selected alternating-spin model.

The exact maximum for this OBC Hamiltonian can be derived in
the local spin-product basis. Each spin-1/2 site supports
at most one outgoing single-site transverse-field flip,
and each spin-1 site supports at most two. A nearest-neighbor
spin-1/2/spin-1 bond allows at most two simultaneous
two-site flips. Including one diagonal entry gives

\[
\mathrm{MNTE}_{\rm model} =
1+\frac L2+L+2(L-1) = \frac{7L}{2}-1 .
\]

The maximum is attained when all spin-1 sites occupy
\(S^z=0\), because the relevant transverse and
two-spin-flip coefficients are nonzero. The observed
program maxima at \(L=8,10,12\) are
27, 34 and 41, in exact agreement with the derived
bound. The formula is **model- and representation-specific**
and cannot be reused for Hamiltonians with other
DM/Gamma components or longer-range interactions
without a new derivation.

The two cache arrays of the current GNU Fortran build
occupy \(20(\mathrm{MNTE}+1)\mathcal D\) bytes,
where \(\mathcal D=6^{L/2}\). Thus the model-specific
allocation reduction relative to the generic default
is \(20[(10L-8)-(7L/2-1)]6^{L/2}\) bytes.
The predicted cache saving and measured peak-RSS
reduction agree within about 1.2% for ALG=1,
OMP=1 across the three tested lengths:

| L | Predicted cache saving (MiB) | Measured RSS saving (MiB) |
| ---: | ---: | ---: |
| 8 | 1.112 | 1.100 |
| 10 | 8.602 | 8.619 |
| 12 | 63.182 | 63.080 |

At \(L=12\), ALG=1, OMP=1, the two-run median
peak RSS decreases from 119,452 to 54,858 KiB
with the Hamiltonian-action count unchanged at 140.
For ALG=2 (NOK=2,NOM=16), the RSS decreases
from 148,904 to 84,178 KiB with 442
Hamiltonian applications in both cases.
**72/72 trials** across three lengths, two
solvers and two thread counts pass independent
low-energy, magnetization, nearest-neighbor
correlation and direct-residual verification.

At \(L=12\), a complete cache with calibrated
row width allocates approximately **37.38 MiB**
in the two work arrays, less than the
**50.28 MiB** needed to cache only half
the rows at the conservative default width.
This is a structural memory result, not a
dedicated-node timing comparison between
those two alternatives. The pre-existing
`MNTE` control is not claimed as a new
algorithm. The
[capacity derivation and resource analysis](cpc-mnte-capacity-evaluation.md),
[frozen 72-case measurements](cpc-mnte-capacity-raw-ci.csv)
and plotting script `tests/plot_cpc_mnte.py`
provide independent reproducibility.

### 6.4 Paired memory-time frontier at the same numerical target

The two resource controls can be combined. Rather than
allocating the generic row width and caching only the first
half of the representative states, we may cache **all**
representative states at the proven model-specific
row capacity. The nominal cache allocation is

\[
B_{\rm half,generic}=20(10L-7)\mathcal D/2,\qquad
B_{\rm full,tight}=20(7L/2)\mathcal D,
\]

for the even open-chain lengths and
\(\mathcal D=6^{L/2}\). At \(L=12\),
these cache arrays require 50.279 and
37.375 MiB, respectively. Thus the
tight full cache uses **less memory than
the default half cache** despite storing
the Hamiltonian actions on every basis state.

A separate same-job comparison with
one warm-up and two measured repetitions
per condition now independently validates
**108/108 cases** across \(L=8,10,12\),
OMP=1/2, ALG=1/2 and three cache policies.
The following \(L=12\) numbers are measured
two-run medians, on GitHub-hosted CI,
for one lowest eigenpair and the same
requested residual accuracy.

| OMP | ALG | Policy | H actions | Wall (s) | Peak RSS (KiB) |
| ---: | ---: | --- | ---: | ---: | ---: |
| 1 | 1 | half/generic | 140 | 21.804 | 67,970 |
| 1 | 1 | full/tight | 140 | 1.297 | 54,886 |
| 1 | 2 | half/generic | 459 | 71.144 | 97,448 |
| 1 | 2 | full/tight | 442 | 3.652 | 84,238 |
| 2 | 1 | half/generic | 140 | 14.254 | 67,944 |
| 2 | 1 | full/tight | 140 | 0.892 | 54,944 |
| 2 | 2 | half/generic | 425 | 43.845 | 97,472 |
| 2 | 2 | full/tight | 442 | 3.510 | 84,322 |

For ALG=1 at fixed size the
Hamiltonian-application count is identical,
so the observed wall-time difference is not
explained by more Krylov iterations. For
ALG=2 the action count differs between
policies and must be reported. In every
listed case the model-specific full cache
uses less measured process RSS and less
observed wall time than the generic half cache.
These **shared-hosted-runner** results
motivate but cannot replace controlled
CPU-pinned, statistically repeated benchmarks.

The [full resource analysis](cpc-cache-pareto-evaluation.md),
[108-row trial CSV](cpc-cache-pareto-raw-ci.csv)
and `tests/plot_cpc_pareto.py` reproduce
the calculation, raw timings and SVG/PNG
figure prototypes. This represents a
**model-dependent configuration result
within an existing cache implementation**,
not a new general solver algorithm.

### 6.5 Dedicated-node measurement protocol

We implemented a CPU-affinity pinned benchmark wrapper
([protocol](cpc-controlled-node-protocol.md)) that runs the
independently validated paired-cache cases while recording
source-archive and binary SHA256 hashes, input-file digests,
CPU topology, compiler and linked-library information,
OpenMP/BLAS thread controls and all trial-level output.
The launcher requires explicit operator attestation before
a run is labelled as taking place in an exclusive
allocation; a shared-runner smoke run is labelled
exploratory. This **reproducibility protocol is not itself
evidence of a completed dedicated-node benchmark**.
No isolated-node timing or strong-scaling efficiency
is claimed until the requested physical workloads have
been repeated under verified exclusive CPU allocation
with statistical uncertainty assessment.

## 7. Conclusion and limitations

ED2 provides a generalized exact-diagonalization
setting for unequal spins and anisotropic,
nonconserving interactions beyond the
published near-saturation QS³ model scope.
Its physically meaningful validation
combines low spectral levels with local and
Cartesian correlations, a parameter-conjugate
operator identity and multiple independent
susceptibility formulations. The literature
Hamiltonian check adds an independently
verified spin-convention and oriented
periodic-boundary test.

These finite-size results demonstrate
software correctness and model coverage,
not novel chiral order, exact bulk
properties or a universally advantageous
eigensolver. The full Hilbert-space
dimension remains exponential; finite
small gaps, truncation conventions,
public API scope, compiler portability
and controlled-node performance require
careful documentation.

## References and evidence links

[1] H. Ueda, S. Yunoki, T. Shimokawa,
*Computer Physics Communications* **277**,
108369 (2022). DOI: 10.1016/j.cpc.2022.108369.

[2] A. Lahiri and S. K. Pati,
*Journal of Physics: Condensed Matter* (2022).
DOI: 10.1088/1361-648X/ac4aab.

[3] T.-C. Yi, W.-L. You, N. Wu and A. M. Oleś,
*Physical Review B* **100**, 024423 (2019).
DOI: 10.1103/PhysRevB.100.024423.

[4] P. Weinberg and M. Bukov, *QuSpin: a Python package
for dynamics and exact diagonalisation of quantum many body
systems. Part I: spin chains*, *SciPost Physics* **2**, 003 (2017).
Project documentation: https://quspin.github.io/QuSpin/ .

Primary independent DM response:
https://github.com/ueda1984/ED2/actions/runs/43690041312

Literature and matched-resource validation:
https://github.com/ueda1984/ED2/actions/runs/43702189530

Related development reports:
docs/development/cpc-stationary-dm-application.md
and docs/development/cpc-literature-matched-benchmark.md.
The cache/multipair evaluation is documented in
`docs/development/cpc-cache-multipair-evaluation.md`.


**Unresolved editorial items:** agree author list,
CPC CPiP versus CP category, program release,
source DOI/archive and licensing provenance;
run isolated-node benchmarks and perform
full technical copy-edit before submission.
