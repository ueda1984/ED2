# CPC application: DM-conjugate vector chirality in a size-stationary mixed-spin chain

**Status:** independently verified, open-boundary, finite-size scientific
software demonstration. This is **not** a discovery of a new phase
or evidence for a thermodynamic critical point.

## Motivation and correction to the earlier diagnostic model

The historical mixed-spin regression model in
[`mixed_spin_nonu1_demo.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-stationary-dm-scaling/tests/mixed_spin_nonu1_demo.py)
includes coefficients that increase with the **absolute bond/site
index**. Those regression results correctly verify ED2, but changing
\(L\) changes the local Hamiltonian in a way that prevents an
ordinary finite-size extrapolation to a single stationary bulk model.

This follow-up keeps **all local couplings and fields on a fixed
two-site repeat rule** as \(L=4,6,8,10,12\) changes, while preserving
site-dependent spin magnitudes, non-U(1) anisotropy, off-diagonal
symmetric exchange and directed DM interactions.

The earlier 54-trial nonstationary diagnostic is retained separately.
The size-stationary study forms a more defensible **physics application
figure**, not a retroactive replacement of prior test evidence.

## Hamiltonian and sign conventions

For an open alternating chain \(S_i=\tfrac12,1,\tfrac12,1,\ldots\),
write the zero-based bond index as \(j=0,\ldots,L-2\) and
\(\eta_j=(-1)^j\). The input contains

\[
H(D)=
\sum_{j=0}^{L-2}
\Big[
J_j^x S_j^xS_{j+1}^x+
J_j^y S_j^yS_{j+1}^y+
J_j^z S_j^zS_{j+1}^z+
\Gamma_j^z(S_j^x S_{j+1}^y+S_j^y S_{j+1}^x)
+
D\eta_j(S_j^xS_{j+1}^y-S_j^yS_{j+1}^x)
\Big]
+\sum_{j=0}^{L-1}(h_j^x S_j^x+h_j^z S_j^z).
\]

All coefficients use the same energy unit, with
**ED2's positive input-field convention**.
The two-site-repeat coefficients are

\[
\begin{aligned}
J_j^x&=0.77+0.05\eta_j,&
J_j^y&=0.63+0.03\eta_j,\\
J_j^z&=1.09-0.04\eta_j,&
\Gamma_j^z&=0.065+0.015\eta_j,\\
h_j^x&=0.15(0.4\eta_j+0.1),&
h_j^z&=0.17+0.035\eta_j .
\end{aligned}
\]

Here the \(\eta_j\) in the field means site parity,
and \(\eta_j\) in the bond means bond-start parity.
All other DM components and symmetric off-diagonal components
are zero. The full Hilbert space of dimension \(6^{L/2}\)
is represented with no fixed-\(S^z\) restriction or NOD projection.
Because of the \(J^x\ne J^y\), \(\Gamma^z\) and transverse
field terms, this is genuinely **non-U(1)**.

The choice of staggered \(D_z\) is deliberate: the *oriented*
bond \((j,j+1)\) has coefficient \(D_j^z=D\eta_j\).
Reversing the bond orientation without simultaneously transforming
the DM sign would describe a different Hamiltonian.

## Observable and three independent verification routes

Define the **z vector chirality** of an oriented bond as

\[
\kappa_j^z=
\langle S_j^xS_{j+1}^y-S_j^yS_{j+1}^x\rangle,
\qquad
K=\sum_{j=0}^{L-2}\eta_j \kappa_j^z.
\]

The script independently constructs the operator
\(\widehat K=\partial H/\partial D\) using SciPy's
spin-product matrices. It then checks

\[
\frac{\partial E_0}{\partial D}
=\langle\widehat K\rangle
=K,
\qquad
\bar\kappa_{\rm stag}=\frac{K}{L-1}.
\]

At \(D=0\), three independent numerical routes are checked:

1. **Direct operator:** \(\langle\widehat K\rangle\),
   computed from the SciPy ground vector.
2. **Energy finite difference:**
   \([E_0(h)-E_0(-h)]/(2h)\) with \(h=2\times10^{-4}\).
3. **ED2 correlation output:** NN `xy-yx` from the
   actual `two_body_cf_xyz.dat` columns. **For this file
   `xy` is column 6 and `yx` is column 4 (1-based),
   including the two site-label columns**.

The **bond-normalized chiral susceptibility**
is defined as

\[
\chi_\kappa(L)=
-\frac1{L-1}\left.
\frac{\partial^2E_0}{\partial D^2}\right|_{D=0}
=-\frac1{L-1}\left.
\frac{\partial K}{\partial D}\right|_{D=0}.
\]

The numerical calculation verifies both finite-difference forms.
For \(L=4,6\), an **additional fully dense independent
spectral calculation** checks

\[
\chi_\kappa(L)=
\frac{2}{L-1}\sum_{m>0}
\frac{|\langle m|\widehat K|0\rangle|^2}{E_m-E_0}.
\]

This formula presumes a **nondegenerate** finite-system ground state.
Ground-energy concavity as a function of the linearly coupled \(D\)
requires \(K(D)\) to be nonincreasing; every sampled
size and DM sweep passed this mathematical consistency check.

The value at \(D=0\) is **not** identified as a
spontaneous chiral order parameter. Nonzero staggered chirality
can arise from the other explicitly anisotropic couplings,
finite boundaries and site-dependent fields.

## Reproducible numerical results

[GNU Fortran 13 / OpenBLAS / OMP=1,2 scientific validation](https://github.com/ueda1984/ED2/actions/runs/43689548834).

The full-size scans were built from a **separately extracted
source distribution**. The reference uses independent SciPy
CSR matrices and ARPACK; four- and six-site full Hamiltonian
matrices were also compared element by element against
a second, dense NumPy construction.

For \(D=0\):

| Sites | Hilbert dimension | \(E_1-E_0\) | \(\bar\kappa_{\rm stag}\) | \(\chi_\kappa(L)\) |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 36 | 0.17235306 | +0.000383416 | 0.63431120 |
| 6 | 216 | 0.22526115 | +0.000218011 | 0.56873293 |
| 8 | 1,296 | 0.28340937 | +0.000156291 | 0.53735526 |
| 10 | 7,776 | 0.32938732 | +0.000133476 | 0.51970132 |
| 12 | 46,656 | 0.35158094 | +0.000123044 | 0.50850099 |

Over the full DM scan \(-0.24,-0.12,0,+0.12,+0.24\),
the bond-normalized staggered chirality is:

| Sites | \(D=-0.24\) | \(D=-0.12\) | \(D=+0.12\) | \(D=+0.24\) |
| ---: | ---: | ---: | ---: | ---: |
| 4 | +0.149067 | +0.076062 | -0.075254 | -0.148151 |
| 6 | +0.134737 | +0.068232 | -0.067733 | -0.134045 |
| 8 | +0.127668 | +0.064465 | -0.064099 | -0.127120 |
| 10 | +0.123635 | +0.062350 | -0.062035 | -0.123157 |
| 12 | +0.121066 | +0.061010 | -0.060720 | -0.120624 |

The response is smooth, sign-changing and visibly
size-dependent over these **finite** open chains.
The monotonic reduction of \(K(D)\) with increasing
\(D\) is consistent with the general variational
concavity of the ground-state energy.

**54/54 ED2 trials passed** (24 for OMP=1, sizes
4--10; 24 for OMP=2, sizes 4--10; six additional
OMP=1 trials at 12 sites): no solver nonconvergence,
no residual-contract violations, and independent
energy, magnetic, all-nine-NN-correlation and chiral
observable comparisons passed. The largest detected
energy mismatch was below \(1.4\times10^{-13}\);
the largest staggered-chirality-sum difference was
below \(1.5\times10^{-9}\) and the \(D=0\) energy
derivative/chirality discrepancy was below
\(2.5\times10^{-10}\) across **all five independently
referenced sizes**.
The exact case-by-case errors and compiler provenance
are retained in the linked workflow artifacts.

At \(D=0\), the two susceptibility finite-difference
routes agree within a stated \(2\times10^{-4}\)
acceptance threshold. For the 4/6-site calculations,
the direct spectral formula agrees closely with
the response estimate, providing a separate
check beyond ED2 and the finite difference.

Data suitable for plotting:
[cpc-stationary-dm-reference.csv](cpc-stationary-dm-reference.csv).
The committed CSV includes values **rounded to the precision
printed in the verified CI logs**; the
`physical-scan.csv` file in each CI artifact
contains the full numerical values and a
case-by-case `report.json`.

### Limitations and what not to infer

- **No thermodynamic extrapolation yet.** The coefficients
  are now stationary with a two-site unit cell, but
  open-boundary \(L=4\ldots12\) data are too limited
  to claim a converged bulk value or a critical exponent.
- The near-zero \(D=0\) chiral signal is not evidence of
  spontaneous chirality.
- \(L=12\) is tested at OMP=1; OMP=2 is currently
  verified up through \(L=10\).
- ALG=1 requests one eigenpair, ALG=2 requests two at
  \(D=0\), and achieved residuals are not matched.
  These trials are **not** a fair timing comparison.
- This is a carefully controlled **physics application
  and software-validation example**, not a novelty
  claim about chiral susceptibility or Hellmann--Feynman.
- A CPC paper will benefit from selected real
  materials or literature-grounded Hamiltonians,
  careful convergence checks, and repeated performance
  measurements on controlled nodes before submission.

## Reproduce

The implemented script supports two separate model profiles:
`diagnostic` for the historical nonstationary
regression Hamiltonian and **`stationary`** for
the fixed two-site repeat rule used in the
scientific data above.

```sh
python3 tests/mixed_spin_dm_chirality_scan.py \
  --exe /absolute/path/to/QS3ED2 \
  --profile stationary \
  --sizes 4,6,8,10 \
  --omp 1 \
  --output /new/nonexistent/stationary-dm-omp1

python3 tests/mixed_spin_dm_chirality_scan.py \
  --exe /absolute/path/to/QS3ED2 \
  --profile stationary \
  --sizes 12 \
  --omp 1 \
  --output /new/nonexistent/stationary-dm-n12
```

To create reproducible **three separate** SVG/PNG
figure prototypes from the committed CSV, install
Matplotlib in the analysis environment and run:

```sh
python3 tests/plot_cpc_dm_chirality.py \
  --csv docs/development/cpc-stationary-dm-reference.csv \
  --outdir /new/nonexistent/cpc-dm-figures
```

The plotting script produces `dm_chirality_vs_D`,
`dm_gap_vs_D`, and `dm_chiral_susceptibility_vs_L`
as **individual** SVG and PNG files without
assuming or fitting a thermodynamic extrapolation.
These figures use values rounded in the GitHub
validation logs; use the full-precision numerical
CSV in the Actions artifact for manuscript-quality
quantitative reanalysis.

All inputs, wavefunction-independent observable
files, solver diagnostic sidecars, independent
physics CSV and machine-readable
`report.json` are retained per run.
Do not reuse an output directory because the
legacy ED2 programs can append to output files.
