# CPC numerical accuracy: a physical near-degenerate low doublet

**Status:** independently verified physical numerical characterization;
**not** a change to ED2's solver algorithm, numerical defaults or
public library API.

## Why this experiment is necessary

ED2 Phase 1-G's optional residual threshold certifies that a
**returned Ritz pair** has a sufficiently small direct residual.
It does **not** independently prove that the returned pair is the
lowest eigenpair, nor does it guarantee a unique eigenvector or
state-specific observable near degeneracy.

Earlier matched-start ALG=1/2 studies used well-separated eigenvalues.
This experiment deliberately controls the **lowest gap**, while
keeping the gap from the lowest doublet to the rest of the spectrum
macroscopically separated.

## Constructed physical Hamiltonian

Take the existing generalized mixed-spin chain with an alternating
`S=[1/2,1,1/2,1,...]`, XYZ exchange, oriented DM,
`Gamma_xy` and transverse fields. Remove the bond(s) to the first
`S=1/2` spin, leaving an independently interacting and non-U(1)
active chain. On the isolated spin apply only `+delta * S_1^z`:

\[
H(\delta)=I_2\otimes H_{\mathrm{active}}+
\delta\,S_1^z\otimes I_{\mathrm{active}}.
\]

The low-energy spectrum is therefore an *exact* shifted copy of the
active-chain spectrum. For `|delta|` well below the active gap,

\[
E_0=e_0-|\delta|/2,\quad
E_1=e_0+|\delta|/2,\quad
E_1-E_0=|\delta|.
\]

For \(\delta>0\), the **unique** ground state satisfies
\(\langle S_1^z\rangle=-1/2\); for \(\delta<0\) it is
\(+1/2\). For \(\delta=0\), the ground space is exactly
two-dimensional, individual state vectors are not unique, and the
**equal-weight doublet ensemble** has
\(\langle S_1^z\rangle_{\mathrm{ensemble}}=0\).

The independent Python reference verified the **entire dense matrix**
against the explicit tensor-product factorization, and confirmed that
\([H,S^z_{\mathrm{tot}}]\ne 0\). The low-doublet external gap is
about 0.24 for four sites and 0.26 for six sites. These are
adversarial *numerical-verification* Hamiltonians, not new interacting
phase examples.

## Reproducibility and test design

Test script:
[`tests/mixed_spin_near_degenerate.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-near-degenerate-physics/tests/mixed_spin_near_degenerate.py).

Validated [GNU 13, OMP 1/2 execution](https://github.com/ueda1984/ED2/actions/runs/43402777979).
The ED2 binary is compiled **only from a newly extracted release
source archive**; NumPy constructs the reference matrix without
using ED2's operator implementations.

The experiment covered \(N=4,6\), \(\delta=10^{-2},10^{-6},
-10^{-6},10^{-10},0\) in the explicitly selected combinations,
ALG=1 and ALG=2 with distinct `NOK/NOM`, residual thresholds
`1e-8` and selected `1e-12`, and OMP=1 and OMP=2.

**64/64 trials** were classified with zero harness violations.
All 64 solves returned `status=0` and satisfied their requested
direct residual policy. However, this is **not** equivalent to
64/64 individually accurate ground-state observables.

## Key accepted-but-inaccurate ground-observable examples

For \(\delta=10^{-10}\), the unique ground-state free-spin
magnetization must be \(-1/2\). ALG=1 returned these accepted
solutions:

| Sites | OMP | Residual request | Direct residual | Measured \(\langle S_1^z\rangle\) | \(H\) applications |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 1 | 1e-8 | 6.377e-11 | +0.1953038858 | 36 |
| 4 | 1 | 1e-12 | 2.169e-13 | -0.499999999999 | 214 |
| 6 | 1 | 1e-8 | 3.409e-9 | -0.2722090404 | 250 |
| 6 | 1 | 1e-12 | 2.594e-13 | -0.500000000000 | 1,348 |
| 4 | 2 | 1e-8 | 5.761e-11 | +0.2506682114 | 36 |
| 4 | 2 | 1e-12 | 2.286e-13 | -0.499999999995 | 214 |
| 6 | 2 | 1e-8 | 3.313e-9 | -0.1031745503 | 250 |
| 6 | 2 | 1e-12 | 2.362e-13 | -0.500000000000 | 1,348 |

The four loose-threshold cases above are **not** residual-policy
violations. They are accepted *low-doublet* approximate eigenpairs,
but they fail to identify the unique lowest state well enough to
reproduce its **state-specific** observable.

Making the residual stricter resolved these particular examples,
at the cost of more Hamiltonian applications. It is **not a general
solution** to global-ground-state certification: an exact excited
eigenvector would have zero residual yet remain excited.

For the **exactly degenerate** cases \(\delta=0\), the accepted
ALG=1 and ALG=2 vectors gave a range of isolated-site
\(\langle S_1^z\rangle\) values. This does **not** indicate
physical inaccuracy: one must instead compare the *projector* onto
the whole ground doublet, or basis-invariant doublet averages.

## A direct excited-state-start counterexample

To separate **residual accuracy** from **ground-state inclusion**
even more sharply, the follow-up test prepares the second eigenvector
with ED2 ALG=2, verifies that eigenvector's eigenvalue against the
independent dense reference, and passes the resulting **native-basis**
vector to ALG=1 via `ED2_INITIAL_VECTOR`. The source eigenvector
is not constructed by guessing ED2's internal basis ordering.

This is still a *physical* four-site mixed-spin model, with the
isolated edge spin splitting \(\delta=10^{-6}\). The true ground
state has \(\langle S_1^z\rangle=-1/2\) and the first excited
state has \(+1/2\).

[Independent OMP=1/2 witness execution](https://github.com/ueda1984/ED2/actions/runs/43403308571):

| Quantity | OMP=1 | OMP=2 |
| --- | ---: | ---: |
| Independently verified \(E_0\) | -1.533921122613122 | -1.533921122613122 |
| Independently verified \(E_1\) | -1.533920122613121 | -1.533920122613121 |
| ALG=1 returned energy | -1.533920122613122 | -1.533920122613122 |
| ALG=1 direct residual | 1.279e-15 | 5.960e-16 |
| Requested absolute residual | 1e-12 | 1e-12 |
| ALG=1 success status | 0 | 0 |
| Returned \(\langle S_1^z\rangle\) | +0.5 | +0.5 |

Thus ALG=1 returned a **mathematically correct excited eigenpair**
under the specified residual policy, but **not the lowest eigenpair**.
This remains true even though the residual is three orders of
magnitude smaller than the requested tolerance. It illustrates a
fundamental Krylov starting-space limitation, **not** an ED2-specific
bug or a failure to enforce the documented residual contract.

Reproduce with the separate
[`tests/mixed_spin_excited_start_probe.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-near-degenerate-physics/tests/mixed_spin_excited_start_probe.py)
and retain both the source and candidate-wavefunction records.

## What can be certified from residual plus a *known external gap*

Let \(P_{\mathcal C}\) be the exact spectral projector on the two
lowest eigenvalues \(E_0,E_1\), and let \(E_2\) be the next
eigenvalue. For a normalized Ritz vector \(\psi\), Ritz value
\(\lambda<E_2\), and direct residual
\(r=\|(H-\lambda)\psi\|\), Hermitian spectral decomposition gives

\[
\|(I-P_{\mathcal C})\psi\|\le
\frac{r}{E_2-\lambda}.
\]

This is a **conditional a posteriori** bound requiring knowledge
of \(E_2\) for the actual Hamiltonian. It is **not** an ED2
internal lowest-state guarantee. For two orthonormal returned
vectors, the combined Frobenius leakage is bounded by the
root-sum-square of the two individual bounds.

When \(\delta=10^{-10}\), the two-dimensional cluster stays well
separated (\(E_2-E_1\approx0.25\)) even though its individual
members are not well separated. Therefore \(r/(E_2-\lambda)\)
may be tiny while \(r/(E_1-E_0)\) is of order one or much larger.
This distinction explains the observed outcome:

- **Low-energy doublet membership can be reliable.**
- **Individual ground-state magnetization can be unreliable.**

An internal diagnostic should not use the word "certified ground
state" solely because the residual passes. A reliable package may
report `residual_accepted`, `lowest_inclusion_not_assessed`,
and `gap_sensitivity_warning` as *separate* notions.

## Verified doublet-invariant observables from native ED2 states

In the additional [16-case trace validation](https://github.com/ueda1984/ED2/actions/runs/43404082525),
ED2 ALG=2 returned the lowest two eigenvectors for 4- and 6-site
chains at exact degeneracy (\(\delta=0\)) or splitting
\(\delta=10^{-10}\), with `NOK/NOM=2/16` and
`4/24`, OMP=1 and OMP=2.

The code loads the second eigenvector in its **native ED2 basis**
through the existing `re_wf=1` path (after copying its Fortran
binary record into an isolated first-state readback directory).
It separately computes magnetic and all Cartesian two-site
observables of the first and second eigenvectors, then compares
their equal-weight average with
\[
\langle O\rangle_{\mathcal C}
=\tfrac12\operatorname{Tr}(P_{\mathcal C}O),
\qquad
P_{\mathcal C}=\sum_{j=0}^{1}|\phi_j\rangle\langle\phi_j|,
\]
where the independent NumPy dense eigensolver supplies
the reference low-doublet eigenvectors \(\phi_j\).

**16/16 comparisons passed**: over all reported local
magnetizations and nine Cartesian pair correlations, the maximum
discrepancy was \(<3\times10^{-9}\). The observed two-state
ensemble free-spin magnetization was zero (within floating-point
roundoff), independent of the gauge or rotation inside the
low doublet.

Representative 6-site **exact-doublet** results:

| OMP | `NOK/NOM` | First \(\langle S_1^z\rangle\) | Second \(\langle S_1^z\rangle\) | Mean | Largest mag / correlation error |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | 2 / 16 | +0.122891 | -0.122891 | 0 | 2.43e-10 / 2.31e-10 |
| 1 | 4 / 24 | +0.484167 | -0.484167 | 0 | 2.88e-14 / 2.88e-14 |
| 2 | 2 / 16 | -0.251835 | +0.251835 | 0 | 4.84e-10 / 6.47e-10 |
| 2 | 4 / 24 | +0.490023 | -0.490023 | 0 | 3.64e-15 / 4.05e-15 |

The test also verifies that both native eigenvectors are
approximately orthonormal and records a *conditional* cluster
leakage estimate based on independent knowledge of the external
spectral gap. It is not a universal lowest-subspace guarantee.
The equal-weight ensemble is the **normalized projector**; it
does not assert that every experimental preparation is mixed.

Reproduction script:
[`tests/mixed_spin_doublet_trace.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-near-degenerate-physics/tests/mixed_spin_doublet_trace.py).
The existing ED2 numerical algorithms were not modified.

## Remaining research tasks

1. **Completed:** controlled excited-start tests independently
   verified that a correct excited state can pass even `1e-12`
   absolute residual, in both OMP profiles.
2. **Completed:** the equal-weight doublet projector/trace
   observables were independently verified in **16/16** native-basis
   ED2 readback tests, despite large differences among individual states.
3. Study genuinely *interacting* avoided crossings and
   near-degenerate clusters (this intentionally decoupled
   spin is a transparent numerical counterexample).
4. Do **not** change ALG=1/2 defaults or promise an automatic
   smallest-eigenvalue certificate without a separately reviewed
   algorithmic design and independent evidence.

## Reproduce the physical test

```sh
python3 tests/mixed_spin_near_degenerate.py \
  --exe /absolute/path/to/QS3ED2 \
  --omp 1 \
  --output /new/nonexistent/results-omp1
```

Repeat with `--omp 2`. The output contains each case's input
deck, full stdout, diagnostics, native-binary eigenvectors and
machine-readable `report.json`; the linked GitHub Actions jobs establish completion of the physical
64-case matrix and the two additional excited-start probes.
