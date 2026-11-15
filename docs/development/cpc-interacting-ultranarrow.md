# CPC: interacting eight-site doublets near the residual-resolution limit

**Status: 24/24 numerical trials accepted**, two OpenMP settings,
three weak longitudinal biases, two restarted-Krylov workspaces
and two requested absolute residual thresholds.

- [Successful source-archive CI run #43706362016](https://github.com/ueda1984/ED2/actions/runs/43706362016),
  SHA `cf0b73af4dbc84d902efd6f33526f0a5a8c423a5`.
- [Frozen complete 24-case classification CSV](cpc-interacting-ultranarrow-raw-ci.csv):
  source/job identities, gap, requested tolerance, solver status,
  individual edge-magnetization disagreement, doublet correlation
  error and number of Hamiltonian actions.
- Independent SciPy low-energy references, native ED2 outputs
  and per-case diagnostic sidecars reside in the GitHub Actions
  artifacts (temporarily retained for 21 days).

## Why this is a different test from the 20-case interacting baseline

The [interacting eight-site benchmark](cpc-interacting-8site-sparse-doublet.md)
already checked a connected anisotropic XY edge weakly coupled
to the non-\(U(1)\) alternating-spin XYZ/Gamma/DM core,
with first gaps between about \(10^{-4}\) and \(5\times10^{-3}\).
The new test sets

\[
\varepsilon=5\times10^{-7}, \qquad
H_{\rm edge}=\delta S_1^z+\varepsilon
(S_1^x S_2^x+0.6 S_1^yS_2^y).
\]

The three longitudinal biases are
\(\delta=-10^{-7},0,+10^{-7}\).
All other seven-site core couplings are unchanged
from the two-site-periodic stationary XYZ/Gamma/DM model
at \(D=0.12\), but the first bond remains deliberately
modified from that main scientific application.
The full dimension is \(6^4=1296\), without a
false fixed-\(S^z\) symmetry sector.

**The bond is nonzero and the edge-spin Cartesian
component commutators are nonzero, so this is not
the factorized-edge Hamiltonian.** At this extremely
weak coupling, however, the independent edge
linear entropies are approximately \(10^{-13}\) or
smaller: they are close to numerical floating-point
resolution. We therefore do **not** make a quantitative
entanglement claim from those tiny values.

## Spectrum and external separation

| \(\delta\) | First gap \(E_1-E_0\) | External gap \(E_2-E_1\) |
| ---: | ---: | ---: |
| \(-1\times10^{-7}\) | 1.46651847e-7 | 0.339089157 |
| 0 | 1.07269404e-7 | 0.339089190 |
| \(+1\times10^{-7}\) | 1.46651620e-7 | 0.339089157 |

The *internal* gap is between one and two orders
of magnitude above the loose requested residual
\(10^{-8}\), but the *external* gap stays
more than six orders of magnitude larger
than the smallest low-gap scale.

For the independent reference, SciPy constructs
the full sparse Hamiltonian, checks its
Hermiticity and the breaking of total \(S^z\),
computes three low eigenpairs and verifies
their direct residuals and orthogonality.
The computed \(E_2\) is external numerical
reference information, **not** a guarantee
from ED2's own residual sidecar.

## How the 24 trials were classified

Each model is solved at OMP=1 and 2
with ALG=2 (two requested low states),
NOK/NOM=2/16 and 4/24, and absolute
residual requests \(10^{-8}\) and \(10^{-12}\).
A nonconvergent run is allowed in the
test design but must be explicitly reported;
the final CI evidence was:

| Outcome | Cases |
| --- | ---: |
| Classified total | 24 |
| Accepted direct-residual policy | 24 |
| Honesty-classified nonconvergence | 0 |
| Accepted-but-individually-unresolved under the existing \(2\times10^{-4}\) observable threshold | 0 |
| Physical/harness violations | 0 |

The last two counts are **threshold-dependent**
software-test outcomes. They should not be
interpreted as a precise bound on the
individual eigenvector error.

The largest absolute difference between
the ED2 first state's edge \(S^z\) and the
independent reference was approximately

\[
1.90\times10^{-7},
\]

whereas the maximum doublet-averaged
Cartesian correlation disagreement was

\[
\boxed{2.255\times10^{-9}.}
\]

The latter is significantly smaller
than the edge magnetization's individual
absolute difference and remains stable
across solver configurations. This
is consistent with the cluster's large
external spectral gap.

Importantly, tightening the requested
residual from \(10^{-8}\) to \(10^{-12}\)
does **not** monotonically reduce the
measured individual edge-observable error
in this data set. The independent SciPy
reference eigenvector accuracy itself can
be sensitive to the \(10^{-7}\) splitting.
A meaningful higher-precision conclusion
about individual eigenvectors would require
a separately validated reference, preferably
with extended precision or a controlled
effective-doublet construction.
We do **not** attribute this apparent
error floor entirely to ED2.

## An independently projected effective two-level cross-check

To avoid relying only on a potentially ill-conditioned
full eight-site eigenvector, the active **seven-site**
non-\(U(1)\) core is diagonalized independently with
SciPy. For its lowest vector \(|g_{\rm core}\rangle\),
define

\[
m_x=\langle g_{\rm core}|S_2^x|g_{\rm core}\rangle,
\qquad
m_y=\langle g_{\rm core}|S_2^y|g_{\rm core}\rangle .
\]

Projecting the weak edge coupling onto the two
free-edge-spin orientations gives, to first
order in \(\varepsilon\),

\[
H_{\rm eff}=e_{\rm core}I_2+
\frac12\left(
\delta\sigma_z+
\varepsilon m_x\sigma_x+
0.6\varepsilon m_y\sigma_y
\right),
\]

and hence the first-order avoided gap

\[
\boxed{\Delta_{\rm eff}=
\sqrt{\delta^2+\varepsilon^2(m_x^2+0.36m_y^2)}.}
\]

The neglected virtual core excitations introduce
corrections of order \(\varepsilon^2/\Delta_{\rm core}\)
to this effective description. It is **not**
an exact identity for the interacting full
Hamiltonian. Nonetheless, at
\(\varepsilon=5\times10^{-7}\)
the agreement is extremely close:

| Edge bias \(\delta\) | Full independent sparse gap | First-order gap prediction | Absolute difference |
| ---: | ---: | ---: | ---: |
| -1e-7 | 1.4665184711e-7 | 1.4665173024e-7 | 1.17e-13 |
| 0 | 1.0726940403e-7 | 1.0726942706e-7 | 2.30e-14 |
| +1e-7 | 1.4665161974e-7 | 1.4665173024e-7 | 1.10e-13 |

The independent core-to-effective reduction was
checked by both OpenMP CI profiles in
[workflow #43706671626](https://github.com/ueda1984/ED2/actions/runs/43706671626)
with all three jobs successful. A separate
[6-row frozen provenance CSV](cpc-interacting-effective-gap-raw-ci.csv)
retains the source SHA, job IDs, independent
sparse gaps and effective predictions.

This strengthens the interpretation of the
*energy splitting*, but does not by itself
provide a \(10^{-13}\)-accurate state-specific
magnetization reference. Eigenvector
orientation in a nearly degenerate doublet
requires separately controlled precision.

## Krylov work and numerical safety

With NOK/NOM=2/16 in these cases,
typical H applications were 180
at the loose threshold and 252
at the tight threshold. NOK/NOM=4/24
used 520–546 H applications in these
particular trials. All corresponding
physical reference and direct-residual
tests passed.

This is a useful reminder that
a solver reaching a requested residual
does not automatically establish a
uniquely accurate global ground vector,
nor does a smaller workspace guarantee
greater or worse numerical accuracy
in every model. We retain the independent
energy and observable checks as part of
the numerical acceptance criterion.

## Manuscript use and limitations

The ultra-small-gap experiment supplements
the regular scientific DM-response curve,
the non-factorized interacting doublet
tests and the analytically controlled
detached-edge cases. Its defensible
contribution is a **transparent,
independently checked numerical
characterization near the residual
resolution scale**.

Avoid claiming a new avoided-crossing
physical phenomenon, universal solver
ground-state certification or a
high-precision individual observable
error bound from the present data.
The existing numerical kernels, main
branch and upstream repository have
not been modified by this test.
