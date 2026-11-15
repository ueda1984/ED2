# CPC: 8- and 10-site interacting avoided-crossing validation

**Status:** 52/52 independently validated numerical trials,
including 8- and 10-site mixed-spin non-U(1) Hamiltonians.
This is a physical correctness and spectral-sensitivity study,
not a general performance comparison or a certified eigensolver.

## Why an independent sparse reference?

The [four-/six-site study](cpc-interacting-avoided-crossing.md)
used full dense NumPy diagonalization to verify the
anisotropic noncommuting XY bond

\[
H=H_{\mathrm{active}}+\delta S_1^z+
\varepsilon(S_1^xS_2^x+0.6 S_1^yS_2^y).
\]

Here the alternating-spin pattern \(S=(1/2,1,1/2,\ldots)\)
and the active-core XYZ, DM, Gamma and transverse-field
coefficients are **unchanged**, but the number of sites
is increased to 8 and 10, with dimensions
\(6^{N/2}=1,296\) and \(7,776\).

Instead of attempting an impractical full dense diagonalization,
the independent
[`tests/mixed_spin_interacting_sparse810.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-interacting-sparse-810/tests/mixed_spin_interacting_sparse810.py)
uses the separately constructed SciPy CSR Hamiltonian
and ARPACK `eigsh(k=3, which="SA")` to provide
the lowest three eigenpairs, directly assessed residuals,
and the gap to the remainder of the spectrum.
The sparse construction was itself previously checked
against a fully dense operator, element by element, at
four and six sites.

The script verifies full Hilbert-space coverage (no
artificial fixed-\(S^z\) sector, no projected NOD cut),
nonconservation of \(S^z_{\mathrm{tot}}\), a
nonzero low-energy gap and magnetization reversal.
It compares ALG=1's single returned ground candidate,
ALG=2's two returned low pairs (`NOK/NOM=2/16`
or `4/24`), all local Cartesian magnetic components,
nine pair correlations per pair, and the **two-state
equal-weight projector ensemble**. The second native
ED2 wavefunction is passed through the actual
`re_wf=1` observable readback path.

For each size:
- Three biases \(\delta=-0.002,0,+0.002\) are
  tested at \(\varepsilon=0.02\);
- One near-degenerate bias \(\delta=0\) is
  tested at \(\varepsilon=5\times10^{-7}\);
- The near-gap ALG=1 runs also include an
  optional absolute residual requirement of `1e-12`.

Every test is run at OMP=1 and OMP=2, with
one BLAS thread and fixed solver seed. All ED2 executables
are built from a separately extracted source tarball.

## Independent physical-reference results

For \(\varepsilon=0.02\):

| Sites | Full dimension | \(\delta\) | \(E_1-E_0\) | \(E_2-E_1\) | \(\langle S_1^z\rangle\) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1,296 | -0.002 | 0.005021923 | 0.230588 | +0.226737 |
| 8 | 1,296 | 0 | 0.004481714 | 0.231103 | +0.031433 |
| 8 | 1,296 | +0.002 | 0.004789102 | 0.230847 | -0.178997 |
| 10 | 7,776 | -0.002 | 0.005207850 | 0.188469 | +0.218689 |
| 10 | 7,776 | 0 | 0.004688795 | 0.188977 | +0.030103 |
| 10 | 7,776 | +0.002 | 0.004983098 | 0.188706 | -0.171965 |

The central state is **entangled** even at zero bias:
the edge-spin linear entropies at \(\delta=0\)
are about \(2.88\times10^{-4}\) for 8 sites
and \(2.90\times10^{-4}\) for 10 sites.

At \(\varepsilon=5\times10^{-7}\), \(\delta=0\):

| Sites | \(E_1-E_0\) | \(E_2-E_1\) | Independent sparse reference's maximum low-pair residual |
| ---: | ---: | ---: | ---: |
| 8 | \(1.11867022\times10^{-7}\) | 0.236111 | below \(3\times10^{-13}\) |
| 10 | \(1.17055513\times10^{-7}\) | 0.193902 | below \(5\times10^{-14}\) |

The separated two-dimensional cluster is therefore
numerically much better conditioned than the individual
very closely spaced eigenvectors.

## ED2 validation and numerical precision

[Successful CI execution (GNU Fortran 13, OMP=1 and 2)](https://github.com/ueda1984/ED2/actions/runs/43687903591).

**52/52 trials returned success and passed the independent
numerical and observable assertions**, with zero recorded
contract failures and zero reported nonconvergences.

Across **all** tested scenarios and both thread profiles:

| Quantity | Largest observed difference |
| --- | ---: |
| Individual ground-state local magnetization | \(1.42\times10^{-7}\) |
| Individual ground-state Cartesian pair correlations | \(1.25\times10^{-7}\) |
| Two-state mean magnetization versus independent projector | \(5.61\times10^{-9}\) |

The largest *individual* differences arise in
the **smallest-gap** models. In the more strongly
split \(\varepsilon=0.02\) models, individual
magnetization/correlation comparisons were instead
at approximately \(10^{-9}\) or better for the
tested input, solver configurations and compiler.

### Numerical subtlety: the reference eigenvectors have their own error

When \(\Delta_{01}\sim10^{-7}\), a nominally
well-converged independent sparse Ritz vector with
residual \(\sim10^{-13}\) can still be
ill-conditioned at the **individual eigenvector** level.
The rough sensitivity scale \(\|r_{\mathrm{ref}}\|/
\Delta_{01}\) can be \(10^{-6}\) or \(10^{-7}\)
in these test cases. This is **not** a rigorous
standalone ground-inclusion certificate; it is a
warning that residual errors and very small
spectral separations affect both ED2 and the
reference diagonalizer.

Thus it would be scientifically incorrect to interpret
every \(\sim10^{-7}\) individual-magnetization
difference as an ED2 error.
The independent **two-state trace**, with a
third-state separation of order 0.19--0.24,
is more stable and shows disagreement below
\(6\times10^{-9}\) in the measured observables.

These are double-precision reference computations;
high-precision or separately certified spectral
inclusion would be needed for stronger
individual-eigenvector error claims.

## Hamiltonian action counts: characterize, do not over-interpret

On the 10-site \(\varepsilon=0.02\) family,
the single-state ALG=1 required approximately
768--856 Hamiltonian applications, whereas
ALG=2 with `NOK/NOM=2/16` used 288 and
`4/24` used 468, all passing the
independent physical checks. The small-gap runs
have a different convergence pattern and
can require more Lanczos applications.

**This is not a speedup claim:** ALG=1
requested one eigenpair and ALG=2 requested
two, and achieved residual accuracies differed.
Timing was not replicated under controlled
dedicated-node conditions, and reorthogonalization
or memory costs are not represented by H-action
counts alone.

## Development and paper implications

The main contribution remains **model coverage plus
independently verified scientific output**:
the generalized ED2 Hamiltonian supports
mixed spins, anisotropic exchange, complex
oriented DM terms, Gamma couplings, transverse
fields, interacting avoided crossings and
their one-/two-site observables.

The present larger reference uses a
*synthetically controlled low-energy crossing*
for verification. It does not establish
a new magnetic phase, universal eigenpair
inclusion guarantee, or performance superiority
over the published QS³ program. A physically
motivated application with scaling/correlation
trends and repeated performance benchmarks
is still desirable before CPC submission.

## Reproduce

```sh
python3 tests/mixed_spin_interacting_sparse810.py \
  --exe /absolute/path/to/QS3ED2 \
  --omp 1 \
  --output /new/nonexistent/sparse810-omp1
```

Repeat with `--omp 2` and a different
new output path. NumPy/SciPy are used only in
this independent test. The output includes
input decks, ED2 solver and observable records,
native wavefunctions and a per-trial
machine-readable `report.json`.
