# CPC numerical verification: eight-site sparse-reference low doublet

**Status: validated 12/12 independent comparisons** on
[GitHub Actions #43705159096](https://github.com/ueda1984/ED2/actions/runs/43705159096),
using OMP=1 and OMP=2 and extracting/building from the source
distribution. The test is at source SHA
`eae33e95cbeb0c25937af4ebc2a8a24c2b3fca58`.

This is a **numerical stress test**, not a thermodynamic
or physical-application result for the unmodified stationary
Hamiltonian. The principal CPC scientific application
[here](cpc-stationary-dm-application.md) uses a connected chain.

## 1. Construction and what is tested

Start with the same two-site-periodic alternating spin
\(S=1/2,1,\ldots\) XYZ/Gamma_xy/oriented staggered-DM and
transverse/longitudinal-field chain as the primary scientific
application, at \(D=0.12\). **For this stress test only**, delete
the bond to the first spin \(S_1=1/2\), and replace its field by
a small positive longitudinal bias \(\delta\).

The eight-site Hilbert space has dimension \(6^4=1296\),
and the remaining **seven-site active core** remains fully
interacting and non-\(U(1)\), with the original stationary local
couplings among its sites. For a suitable ordering of local
tensor factors, the independently constructed sparse Hamiltonian
obeys

\[
H_{\rm stress}(\delta)=
I_2\otimes H_{\rm core}+
\delta S^z_{\rm edge}\otimes I_{\rm core}.
\]

Let \(e_0<e_1\) denote the two lowest eigenvalues of the
active core. For \(|\delta|\ll e_1-e_0\), the lowest
doublet is

\[
E_0=e_0-\frac{|\delta|}{2},\qquad
E_1=e_0+\frac{|\delta|}{2},\qquad
E_2=e_1-\frac{|\delta|}{2}.
\]

The inner gap is \(|\delta|\), whereas the spectral
separation between the doublet and the third level is
\(e_1-e_0-|\delta|\).

SciPy independently assembles both the seven-site
core and full eight-site Hamiltonian in the product
basis using sparse tensor products, verifies the
tensor-factor identity, checks the nonzero commutator
with total \(S^z\), and computes the core's low
eigenpair with a direct residual check.

This formulation provides an independent and unambiguous
**low-subspace projector**:

\[
P_{\mathcal C}=I_2\otimes|g_{\rm core}\rangle
\langle g_{\rm core}|.
\]

For any physical observable \(O\), the doublet
ensemble expectation is

\[
\langle O\rangle_{\mathcal C}
=\frac12\operatorname{Tr}(P_{\mathcal C}O).
\]

Unlike an individual eigenvector in an exactly or nearly
degenerate doublet, this quantity is invariant under any
unitary rotation of the chosen basis within the doublet.

## 2. Native ED2 readback and independent acceptance

For \(\delta=0,10^{-10},10^{-6}\) we test
ALG=2 with NOK/NOM=2/16 and 4/24,
separately at OMP=1 and OMP=2.
**All twelve solver runs passed** the requested
absolute direct residual threshold \(10^{-8}\),
their two-state inclusion checks and the
independent doublet-observable comparisons.

The first vector's local moments and every Cartesian
two-body correlation are read from the ordinary ED2
outputs. The second **native Fortran wavefunction record**
is copied to a new input directory and evaluated through
ED2's own `re_wf=1` readback path.
The averaged ED2 output is compared with the independently
constructed sparse low-doublet projector for **all
8×3 local moments** and **28×9 Cartesian two-spin
correlations**. This checks ED2's observable
implementation without assuming its internal basis
layout matches the SciPy product basis.

The external gap remains approximately
\(\Delta_{\rm outer}=0.33909\) at \(\delta=0\);
the largest error in any doublet-averaged
magnetization or correlation over the test
matrix is **below \(3.7\times10^{-10}\)**.

### Representative exactly degenerate results

At \(\delta=0\), the individual isolated-spin
\(S_1^z\) expectation can vary substantially
with the solver workspace and OMP profile:

| OMP | NOK/NOM | State 1 \(\langle S_1^z\rangle\) | State 2 \(\langle S_1^z\rangle\) | Doublet mean |
| ---: | --- | ---: | ---: | ---: |
| 1 | 2/16 | -0.0250678223 | +0.0250678223 | \(\approx0\) |
| 1 | 4/24 | +0.0326556298 | -0.0326556298 | \(\approx0\) |
| 2 | 2/16 | -0.4004892609 | +0.4004892609 | \(\approx0\) |
| 2 | 4/24 | +0.0274481190 | -0.0274481190 | \(\approx0\) |

These differences are **not** physically distinct
ground-state order parameters. The two-dimensional
ground subspace is exactly degenerate, and its
individual orthonormal basis vectors are not unique.
The trace is physical and basis independent.

At the extremely small but nonzero
\(\delta=10^{-10}\), these particular
runs selected edge-spin expectations close to
\(-1/2\) and \(+1/2\); such behavior should
**not** be presented as an automatic guarantee
at small gaps. Earlier [adversarial tests](cpc-near-degenerate-gaps.md)
show that even residual-accepted states can
return inaccurate individual ground observables
or entirely miss the true lowest state when
the starting subspace is unsuitable.

## 3. Interpretation and reproducibility

For a normalized vector with Ritz value \(\lambda<E_2\)
and residual \(r=\|(H-\lambda)\psi\|\), the
reference-based conditional bound on leakage
outside the doublet is

\[
\|(I-P_{\mathcal C})\psi\|
\le\frac{r}{E_2-\lambda}.
\]

This requires *external knowledge* of the
correct third eigenvalue \(E_2\). ED2's internal
residual gate cannot independently certify
global ground-state inclusion.

The test establishes that low-energy
subspace-invariant observables are reliable
even when individual basis vectors change
substantially under restart workspaces or
thread-count changes. It complements the main
chiral-response verification and the
resource characterization on the **connected**
size-stationary Hamiltonian, rather than
replacing either one.

Reproduce with:

    python3 tests/cpc_stationary_low_doublet_trace.py \
      --exe /path/to/QS3ED2 --omp 1 \
      --sizes 8 --deltas 0,1e-10,1e-6 \
      --workspaces 2:16,4:24 \
      --output /new/uncreated/sparse-doublet-omp1

Run separately with `--omp 2`. Both runs
preserve full input decks, diagnostics, binary
wavefunctions, ED2 correlation outputs, native
second-state readback and independently computed
`report.json`, available as short-retention
CI artifacts in the linked workflow.

**Limits:** this is an exactly factorized edge
stress construction; a physically interacting
near-degenerate cluster with 8 or more sites,
a genuine three-state cluster, and associated
projector/individual-observable checks are
further tests rather than implied results.
