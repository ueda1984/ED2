# CPC: diagnosing individual-state reference precision at an ultra-small gap

**Verified numerical result:** the apparent \(\sim10^{-7}\)
difference between ED2's ground-state edge magnetization and
the original SciPy sparse-`eigsh` reference is **not sufficient
evidence of an ED2 error**. A separate core-projected
second-order effective model and three dense LAPACK drivers
agree much more closely with ED2. The raw comparison now
makes the numerical-method dependence explicit.

This is a **finite-size numerical reliability study** for the
modified weak-edge, eight-site alternating-spin XYZ/Gamma/DM
Hamiltonian, not a newly discovered bulk phase or a revised
DM susceptibility of the primary connected chain.

## 1. Evidence sources and scope

Three independently archived sources are used:

- [24 ED2 doublet cases versus a Q-resolvent second-order model](cpc-ultranarrow-ed2-second-order-raw-ci.csv)
  from [CI run #43707339560](https://github.com/ueda1984/ED2/actions/runs/43707339560),
  source SHA `6635e4b6e83eab0c2ccfd9933d27745a8f6a55fd`,
  with OMP=1/2, NOK/NOM=2/16 and 4/24,
  and requested direct residual \(10^{-8}\) or \(10^{-12}\).
- [Nine independent dense-LAPACK results](cpc-ultranarrow-dense-lapack-raw-ci.csv)
  from [CI run #43707711743](https://github.com/ueda1984/ED2/actions/runs/43707711743),
  source SHA `a1242ae0d4ea4cfc59d1560cd0d8baae33d0b78f`.
  This computes all three low levels and individual edge moments
  with LAPACK `ZHEEVR`, `ZHEEVX` and `ZHEEVD`.
- The [original 24-case ultra-small-gap study](cpc-interacting-ultranarrow.md)
  and the [independent core-effective first-order
  gap data](cpc-interacting-effective-gap-raw-ci.csv).

The dense reference uses the **same independently constructed
SciPy physical Hamiltonian** as the sparse `eigsh` reference,
but it employs mathematically different eigensolver algorithms.
It is therefore a *solver-level numerical cross-check*, not
a third independently specified Hamiltonian. The projected
effective reference starts from a separate seven-site core
calculation and explicitly resolves virtual excitations
through a constrained sparse linear system.

All source data and plotting scripts remain in the
CPC development branch. Neither the ED2 numerical kernel
nor the public fork's main branch has been modified.

## 2. Second-order effective model

The first bond is replaced by

\[
V_{\rm edge}=
\varepsilon S_1^x S_2^x
+0.6\varepsilon S_1^y S_2^y,
\qquad \varepsilon=5\times10^{-7},
\]

while the active seven-site core retains the XYZ,
Gamma_xy, staggered DM and site-field couplings.

Let \(|g\rangle\) denote the core ground state
with energy \(e_g\), and let
\(Q=I-|g\rangle\langle g|\).
Define \(R=Q(H_{\rm core}-e_g)^{-1}Q\),
interpreted as an inverse on the *orthogonal
complement* only. The Python numerical implementation
uses an augmented sparse linear system instead of
attempting to invert the singular full matrix.

For

\[
w_\pm=Q(S_2^x\pm0.6iS_2^y)|g\rangle
\]

the two diagonal self-energies are

\[
\Sigma_\uparrow=
-\frac{\varepsilon^2}{4}
\langle w_+|R|w_+\rangle,\qquad
\Sigma_\downarrow=
-\frac{\varepsilon^2}{4}
\langle w_-|R|w_-\rangle.
\]

The second-order effective Hamiltonian in
the edge-spin subspace is

\[
H_{\rm eff}^{(2)}=
e_gI+
\begin{pmatrix}
\delta/2+\Sigma_\uparrow &
\varepsilon(m_x-0.6im_y)/2\\
\varepsilon(m_x+0.6im_y)/2 &
-\delta/2+\Sigma_\downarrow
\end{pmatrix}.
\]

Its low-state eigenvector yields the
effective edge \(S_1^z\) expectation.
First-order correction from the
virtual excited-core amplitudes is
also included when calculating the
normalized *dressed* spin expectation.

At this extremely weak coupling,
the omitted energy-dependent and
higher-order terms scale at least
quadratically in the small parameter
\(\varepsilon/\Delta_{\rm core}\),
but no strict arbitrary-precision
error enclosure has been established.

## 3. What the different methods actually report

| Edge bias \(\delta\) | Original sparse `eigsh` edge \(S_1^z\) | Dense `ZHEEVR` | Second-order core reduction |
| ---: | ---: | ---: | ---: |
| \(-10^{-7}\) | 0.3409440783 | 0.3409440978 | 0.3409441010 |
| 0 | \(9.1074011\times10^{-7}\) | \(7.3714597\times10^{-7}\) | \(7.3738945\times10^{-7}\) |
| \(+10^{-7}\) | -0.3409435354 | -0.3409435237 | -0.3409435238 |

The largest ED2/second-order edge moment
disagreement among all 24 solver trials is
\(9.95\times10^{-9}\). By comparison,
the largest ED2/original-sparse-`eigsh`
disagreement in the same trials is
\(1.75\times10^{-7}\), and the
largest dense-`ZHEEVR`/original-sparse
disagreement is \(1.80\times10^{-7}\).

At \(\delta=0\), for example:

- Original sparse `eigsh` vs core reduction:
  absolute difference \(1.73\times10^{-7}\).
- Dense `ZHEEVR` vs core reduction:
  absolute difference \(2.43\times10^{-10}\).
- Dense `ZHEEVD` vs core reduction:
  absolute difference \(9.02\times10^{-9}\).
- ED2 across all eight OMP/restart/residual
  settings vs core reduction: absolute
  differences ranging from about
  \(1.82\times10^{-9}\) to \(8.95\times10^{-9}\).

This pattern is **consistent with an eigenvector
conditioning problem in the original sparse
reference**. It does not support labelling
the \(10^{-7}\)-scale ED2/sparse difference
as an actual physical ED2 ground-state error.

## 4. Gap cross-check is more stable than individual states

The internal low gap is
\(E_1-E_0\approx1.07\times10^{-7}\)
at zero field, and the external doublet
separation remains approximately 0.339.

The dense `ZHEEVR` zero-bias gap is
\(1.072694271186947\times10^{-7}\);
the second-order reduction gives
\(1.072694270570897\times10^{-7}\).
Their double-precision difference is
approximately \(6.2\times10^{-17}\).

At the two nonzero biases, their
differences are approximately
\(1.7\times10^{-15}\) and
\(8.5\times10^{-16}\).
These small differences are consistent
with the perturbative estimate and
numerical floating-point limitations.
They are **not** rigorous proof of
15–17-digit ground-state eigenvector accuracy.

## 5. Why residual-accepted eigenvectors can disagree

Suppose a normalized approximate ground
vector has Rayleigh value \(\lambda<E_1\)
and residual \(r=\|(H-\lambda)u\|\).
Given the *external* true first-excited
energy \(E_1\), the spectral theorem gives

\[
\|(I-|g\rangle\langle g|)u\|
\le\frac{r}{E_1-\lambda}.
\]

Thus a direct residual of \(10^{-9}\)
with an internal gap near \(10^{-7}\)
only gives an eigenvector-angle bound
of order \(10^{-2}\), even if the actual
error happens to be much smaller.
This explains why an excellent
eigenpair residual does not numerically
certify state-specific observables
to \(10^{-9}\).

For the rank-two low cluster, the analogous
leakage bound instead uses the much larger
**external** separation \(E_2-E_1\approx0.339\).
It is therefore natural that the two-state
Cartesian correlation traces remain stable
at about \(10^{-9}\) even when individual
ground eigenvectors shift noticeably
between numerical eigensolvers.

These bounds are conditional on the actual
spectral ordering and independently known
separation. They do not allow ED2's
residual alone to certify that the
correct lowest eigenpair was returned.

## 6. Consequences for the CPC paper

The earlier statement that ED2 exhibits
an individual observable "error" of
approximately \(1.90\times10^{-7}\)
against the sparse reference is
**over-strong**. The defensible
replacement is that there was an
"ED2–sparse-reference disagreement",
and that a dense-LAPACK plus second-order
effective cross-check strongly points to
the sparse iterative reference as the
dominant source of that discrepancy.

The CPC manuscript should lead with:
1. independent reproduction of physical
   Hamiltonians and chirality;
2. solver residual checks and measured
   low-eigenpair inclusion;
3. robust low-cluster projector traces;
4. a careful example of *why reference
   eigensolvers need numerical validation*
   in very small gaps.

We **do not** claim new physics, arbitrary
precision, a universal ED2 convergence
guarantee or a certified high-precision
edge magnetization.

Reproduce the dense numerical reference
with `tests/cpc_ultranarrow_dense_lapack_crosscheck.py`
and the ED2/core reference with
`tests/cpc_interacting_ultranarrow_probe.py`
and `tests/cpc_ultranarrow_second_order_reference.py`.
The [frozen-data audit and two SVG/PNG
comparison figures](../../tests/plot_cpc_ultranarrow_reference_comparison.py)
regenerate the method comparison without
rerunning the full Hamiltonian tests.
