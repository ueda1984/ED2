# CPC: independently verified interacting eight-site low-doublet test

**Status:** all **20/20** numerical tests passed, across two
OpenMP thread counts and two restarted-Krylov workspaces.
[Source-archive GitHub Actions run #43705920631](https://github.com/ueda1984/ED2/actions/runs/43705920631)
executed source SHA
`539e4d47c98e1eaba8ad80f214f311c468f919ce`.
The two numerical jobs both concluded successfully.
The complete run inputs, diagnostics, native Fortran
wavefunctions and second-state observable readback are
retained as temporary GitHub Actions artifacts. A
[frozen 20-case CSV](cpc-interacting-8site-raw-ci.csv)
retains the numeric evidence and job IDs beyond artifact expiration.

**Scientific boundary:** this is an *adversarial numerical test*,
not the actual connected size-stationary DM-response Hamiltonian
used for the main CPC scientific application. We modify its first
bond and field deliberately, while retaining the rest of the
seven-site non-U(1) alternating-spin core.

## Model and independent reference

Use the eight-site OBC chain with local spins
`[1/2,1,1/2,1,1/2,1,1/2,1]` and the original
two-site-periodic XYZ/Gamma_xy/staggered-DM
couplings among sites 2–8, with `D=0.12`.
The first spin is connected to site 2 via

\[
H_{\rm edge}(\delta,\varepsilon)=
\delta S_1^z+
\varepsilon (S_1^x S_2^x +0.6\,S_1^yS_2^y).
\]

The edge is now **genuinely interacting**, unlike the
earlier [exactly decoupled control](cpc-stationary-sparse-doublet.md).
The two transverse couplings fail to commute with each
Cartesian component of the edge spin (independently
checked using sparse commutators), and the full
Hamiltonian remains non-\(U(1)\).

The independently implemented SciPy sparse constructor
computes the lowest **three** eigenpairs, with strict
direct residual and mutual orthogonality checks. The
third eigenvalue is essential: it distinguishes the
small first gap \(E_1-E_0\) from the isolated
two-state cluster's external gap \(E_2-E_1\).

For the independent ground vector, the
edge-spin reduced density matrix yields the
entanglement diagnostic

\[
S_L=1-\operatorname{Tr}(\rho_{\rm edge}^2).
\]

At \(\varepsilon=0.02,\delta=0\), the independent
reference has \(S_L\simeq2.19\times10^{-4}>0\).
This distinguishes the test from an exactly
factorized doublet. Smaller but nonzero
entanglement persists at other sampled biases.

## Independently computed spectrum

All energies are in the input Hamiltonian units.

| Edge coupling \(\varepsilon\) | Bias \(\delta\) | \(E_1-E_0\) | \(E_2-E_1\) | Edge linear entropy |
| ---: | ---: | ---: | ---: | ---: |
| 0.02 | -0.002 | 0.004844376 | 0.333038808 | 5.776e-4 |
| 0.02 | 0 | 0.004296818 | 0.333507930 | 2.188e-4 |
| 0.02 | +0.002 | 0.004630955 | 0.333252813 | 2.526e-5 |
| 0.0005 | 0 | 0.000107270 | 0.338950016 | 1.152e-7 |
| 0.0005 | +0.00002 | 0.000109089 | 0.338948539 | 5.494e-8 |

The moderate-coupling avoided-crossing gap is
**not exactly symmetric** in the longitudinal bias.
No claim is made that the sampled bias points include
the true exact location of its minimum. The table
is a finite numerical stress profile, not a phase
diagram or a thermodynamic prediction.

## ED2 numerical acceptance and projector traces

For each of the five reference Hamiltonians, ED2
ALG=2 requests the lowest **two eigenpairs**.
Both NOK/NOM workspaces `2/16` and `4/24`
are tested at OMP=1 and OMP=2, producing
**5 × 2 × 2 = 20 accepted runs**.

Each accepted case verifies:

1. Solver exit status, requested direct residual
   \(10^{-8}\), requested eigenpair count,
   numerical energy inclusion in the independent
   low-doublet reference, and native eigenvector
   orthonormality.
2. Ground-state local magnetizations and every
   Cartesian two-body correlation against SciPy.
   With near-degeneracy, any unresolved **individual**
   ground observable is classified separately;
   a small residual alone is never claimed to
   certify globally lowest-state inclusion.
3. First-state observables are read from the
   regular ED2 files; the second Fortran
   native-basis eigenvector is independently
   passed through ED2's `re_wf=1` observable
   implementation. Their equal-weight trace
   is then compared with the independently
   constructed sparse two-state spectral projector
   for **24 local moments and 252 pair correlations**.
4. The independent rank-two projector is also
   invariant under an explicit nontrivial
   complex unitary rotation of its two
   eigenvector columns.

The maximum observed error in any doublet
mean local moment or Cartesian pair correlation
across all twenty trials is

\[
\boxed{\epsilon_{\rm trace,max}
=3.7433\times10^{-10}.}
\]

The largest requested-pair direct residual
reported by ED2 is
\(2.676\times10^{-9}\), below
the required \(10^{-8}\) threshold.
Both OMP profiles and both restart workspaces
pass the independent physical checks.

In contrast to the isolated-edge control,
the rank-two projector **cannot** be written as
\(I_2\otimes|g_{\rm core}\rangle
\langle g_{\rm core}|\) when \(\varepsilon\ne0\).
The reference uses the actual two interacting
eigenvectors obtained from the full
1296-dimensional sparse problem.

## An external-gap-based conditional error estimate

Let \(P\) be the *exact* spectral projector
onto the two lowest eigenvectors, with
\(E_2>E_1\), and let \(u_0,u_1\) denote
normalized, orthonormal approximate
eigenvectors with Ritz values
\(\lambda_j<E_2\) and direct residuals
\(r_j=\|(H-\lambda_j)u_j\|\).

The spectral theorem gives the separate
conditional bounds

\[
\|(I-P)u_j\|\le
\frac{r_j}{E_2-\lambda_j}.
\]

For the rank-two returned projector
\(Q=\sum_{j=0}^1|u_j\rangle\langle u_j|\),
we have

\[
\|P-Q\|_{\rm F}
=\sqrt2\|(I-P)U\|_{\rm F}
\le\sqrt{2\sum_{j=0}^1
\frac{r_j^2}{(E_2-\lambda_j)^2}} .
\]

Using \(\operatorname{rank}(P-Q)\le4\)
and the Schatten norm inequalities,
for a Hermitian observable \(O\),

\[
\left|\frac12\operatorname{Tr}[(Q-P)O]\right|
\le\sqrt2\,\|O\|_{\rm op}\,
\sqrt{\sum_{j=0}^1
\frac{r_j^2}{(E_2-\lambda_j)^2}} .
\]

Here local spin components and Cartesian
two-site products have
\(\|O\|_{\rm op}\le1\) for the
selected \(S=1/2,1\) spins.

This is an **external conditional estimate**:
the independent SciPy \(E_2\) supplies the
unknown spectral gap; ED2's residual
alone does not certify it or the global
ground-state inclusion. The finite-precision
Gram matrix and independently computed
reference eigenvalues also introduce
numerical rounding, which the frozen CSV
audit handles separately from the exact
mathematical statement. This bound becomes
useful precisely because the external gap
is about 0.33 while the *inner* doublet
gap may be only \(10^{-4}\).

## Reproduction and figure data

The independent numerical test lives in
`tests/cpc_stationary_interacting_sparse_doublet.py`.
The source archive is extracted afresh and
built with GNU Fortran 13 on Ubuntu 24.04
in [CI #43705920631](https://github.com/ueda1984/ED2/actions/runs/43705920631).
A direct reproduction command is:

    python3 tests/cpc_stationary_interacting_sparse_doublet.py \
      --exe /absolute/path/to/QS3ED2 \
      --omp 1 \
      --output /new/interacting-omp1

Repeat at OMP=2 in a fresh directory.
Raw source logs, diagnostic sidecars,
per-case independent SciPy references and
native wavefunction readback files are
retained in the action artifacts.

The committed
[frozen CSV](cpc-interacting-8site-raw-ci.csv)
contains one row per accepted case,
source SHA, CI job ID, gap, entropy,
OMP/restart settings, residuals and
all trace-observable maximum errors.
`tests/plot_cpc_interacting_sparse_doublet.py`
strictly validates every one of the twenty
rows, checks the external-gap consistency
and derives the conditional residual
bound before producing SVG/PNG
figure prototypes for (i) avoided
gaps, (ii) doublet trace errors and
(iii) independent edge linear entropy.

**Manuscript interpretation:** the ED2
code reliably reproduces basis-invariant
low-subspace observables for this
interacting, nonconserving, near-degenerate
system. The test does not establish
a universally reliable ground-state
identifier, a new bulk quantum phase,
or a software speedup. Source numerical
kernels are unchanged.
