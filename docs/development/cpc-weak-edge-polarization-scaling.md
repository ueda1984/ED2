# CPC: zero-bias weak-edge polarization and a second-order numerical reference

**Validated numerical stress experiment:** [GitHub Actions #43710019532](https://github.com/ueda1984/ED2/actions/runs/43710019532), source SHA `eddc0dafe3f33fca1591460469c8531e508d35e0`. All three jobs passed: 4/4 independent dense-reference cases and 12/12 ED2 two-eigenpair comparisons, with independent observable and direct-residual checks. No Fortran numerical kernel was changed.

**Scope:** this is a deliberately modified first bond in an eight-site alternating spin-1/2/spin-1, non-U(1) XYZ/Gamma/staggered-DM chain. It is a *numerical stress benchmark*, not the connected two-site-stationary DM susceptibility model or evidence of a bulk phase transition.

## 1. Exact coupling-sign symmetry

Keep the full seven-site interacting core and replace the first bond by

\[
H_{\rm edge}(\varepsilon,\delta)
=\varepsilon(S_1^x S_2^x+0.6S_1^y S_2^y)+\delta S_1^z.
\]

A rotation by \(\pi\) about the edge spin's \(z\) axis acts as \(U=2S_1^z\), up to an irrelevant phase. It reverses both transverse edge operators while preserving the remaining core and longitudinal bias. Thus

\[
\boxed{H(-\varepsilon,\delta)=U H(\varepsilon,\delta) U^\dagger.}
\]

The independently assembled 1296-dimensional sparse Hamiltonians pass this identity to \(10^{-12}\) in all four coupling cases. Away from the exact \(\varepsilon=0\) degeneracy, the lowest-state edge \(S_1^z\) expectation must be **even** in \(\varepsilon\).

## 2. Second-order core projection and a nonzero edge moment

Let \(|g\rangle\) and \(e_g\) be the core's ground state and energy; let \(Q=I-|g\rangle\langle g|\). The positive complement resolvent is \(R=Q(H_c-e_g)^{-1}Q\), where the inverse is taken **only in Q space** using an augmented sparse system.

Define \(m_\alpha=\langle g|S_2^\alpha|g\rangle\),
\(w_\pm=Q(S_2^x\pm0.6iS_2^y)|g\rangle\), and

\[
\Sigma_\uparrow=-\frac{\varepsilon^2}{4}\langle w_+|R|w_+\rangle,
\qquad
\Sigma_\downarrow=-\frac{\varepsilon^2}{4}\langle w_-|R|w_-\rangle.
\]

The core-projected effective Hamiltonian is

\[
H_{\rm eff}^{(2)}
=e_gI+
\begin{pmatrix}
\delta/2+\Sigma_\uparrow &
\varepsilon(m_x-0.6im_y)/2\\
\varepsilon(m_x+0.6im_y)/2 &
-\delta/2+\Sigma_\downarrow
\end{pmatrix}.
\]

With \(a=QS_2^x|g\rangle\) and \(b=QS_2^y|g\rangle\), its diagonal shift difference obeys

\[
\Sigma_\uparrow-\Sigma_\downarrow=
0.6\varepsilon^2\,\operatorname{Im}\langle a|R|b\rangle.
\]

At zero edge bias, the first-order two-level Hamiltonian predicts \(\langle S_1^z\rangle=0\). The second-order term instead generates a leading moment

\[
\boxed{\langle S_1^z\rangle=
C|\varepsilon|+O(\varepsilon^2),\quad
C=-\frac{0.6\,\operatorname{Im}\langle a|R|b\rangle}
{2\sqrt{m_x^2+0.36m_y^2}}.}
\]

The absolute value follows from the **exact sign-flip unitary**, not from an assumed analyticity across the degeneracy at zero coupling. Its cusp is a finite-system level-selection effect, **not** thermodynamic phase-transition evidence. The constrained core calculation yields \(C=1.4747790903698443\) for the selected Hamiltonian.

## 3. Four independent LAPACK measurements

The same independently assembled physical sparse Hamiltonian is converted to a full dense Hermitian matrix and diagonalized with LAPACK `ZHEEVR`. Every eigenpair passes direct residual and orthogonality checks. This is *eigensolver-level* independence, not an additional independently specified Hamiltonian.

| \(\varepsilon\) | Dense \(\langle S_1^z\rangle\) | Dense \(\langle S_1^z\rangle/\varepsilon\) | Dense versus second-order gap error |
| ---: | ---: | ---: | ---: |
| 5e-7 | 7.371459725888e-7 | 1.474291945 | 6.16e-17 |
| 2e-6 | 2.948309088004e-6 | 1.474154544 | 6.49e-16 |
| 5e-6 | 7.375266809054e-6 | 1.475053362 | 5.01e-16 |
| 1e-5 | 1.474760547388e-5 | 1.474760547 | 1.48e-15 |

The near-constant ratio over a factor of 20 in coupling agrees with the independently derived asymptotic coefficient. Gap differences this small in double precision do **not** certify state-specific expectations to 15 digits.

[Four-point source/provenance CSV](cpc-weak-edge-scaling-new-ci.csv).

## 4. ED2 two-state validation

ALG=2 computes two states for three biases \(\delta=-10^{-7},0,+10^{-7}\), two workspaces NOK/NOM=2/16 and 4/24, and OMP=1/2, with requested absolute residual \(10^{-12}\). The complete **12/12** matrix passes the native ED2 observable readback and independently assembled sparse low-projector comparisons.

The largest ED2 versus second-order *two-level* ground edge moment difference is \(1.17\times10^{-8}\), (the dressed result differs from this two-level estimate by less than \(10^{-12}\) in these ultraweak cases), and the maximum doublet-averaged Cartesian correlation difference from the independent low projector is approximately \(2.26\times10^{-9}\). A small residual does not, in isolation, certify globally correct lowest-state inclusion.

The [12-row ED2 comparison](cpc-second-order-ed2-new-ci.csv) and [6-row model reference comparison](cpc-second-order-reference-new-ci.csv) preserve the job IDs and source SHA; do not merge them with earlier runs as if they were one measurement session.

## 5. Publication limitations and reproduction

This gives a sensitive **numerical correctness and conditioning** example: a vanishing first-order moment, a second-order induced edge moment, an exact sign-flip symmetry, and quantitative comparisons among ED2, sparse and dense reference solvers. The second-order model omits higher-order core excitations and finite-bias corrections. Neither it nor dense LAPACK supplies rigorous interval-certified eigenvectors.

Keep the connected staggered-DM/chirality response as the primary scientific application, and the matched-target cache memory/time results as the primary software-resource demonstration. The weak-edge construction is an additional adversarial verification, not a claim of newly discovered bulk physics.

Reproduce with `tests/cpc_weak_edge_polarization_scaling.py` and `tests/cpc_ultranarrow_second_order_validate.py`. The [CI workflow](https://github.com/ueda1984/ED2/actions/workflows/cpc-ultranarrow-second-order-extension.yml) preserves source-archive provenance and native output logs.