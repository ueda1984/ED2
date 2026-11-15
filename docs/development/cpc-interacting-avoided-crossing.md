# CPC scientific validation: interacting avoided crossings

**Status:** independently validated numerical evidence for the ED2
development executable; not a physical discovery, a released API,
or a universal ground-state certificate.

## Why extend the decoupled-edge-spin example?

The [previous low-doublet stress test](cpc-near-degenerate-gaps.md)
introduced `delta*S1z` on an exactly isolated `S=1/2` spin.
That example proves the distinction between a residual and
ground-state identification, but its exact factorization limits
the complexity of its physics.

The present study adds a small, **nonzero intersite coupling**:

\[
H(\delta,\varepsilon,\eta)=H_{\mathrm{active}}
+\delta S_1^z
+\varepsilon(S_1^xS_2^x+\eta S_1^yS_2^y).
\]

The active subsystem is an open alternating-spin chain with
`[S_1,S_2,S_3,...]=[1/2,1,1/2,...]`.
Sites 2 through \(N\) retain the earlier site-dependent XYZ,
oriented DM, symmetric off-diagonal Gamma interactions, and
transverse/longitudinal magnetic fields, all with the ED2
\(+\mathbf h\cdot\mathbf S\) sign convention.
The explicit numerical couplings are in the two test scripts
linked below.

Two nested tests are considered:

- **\(\eta=0\):** a single \(S_1^xS_2^x\) coupling. At
  \(\delta=0\), \([H,S_1^x]=0\), so an individual eigenstate
  can be a product at the central gap. For \(\delta\ne0\),
  the edge and core are genuinely entangled. This special
  model permits transparent level-repulsion and observable checks.
- **\(\eta=0.6\):** anisotropic XY coupling in which all three
  edge-spin component commutators are nonzero even at \(\delta=0\).
  The central ground state is entangled, and the minimum of the
  avoided crossing is no longer exactly at zero bias.

The full product Hilbert spaces are dimension 36 for \(N=4\)
and 216 for \(N=6\). The active core violates
total \(S^z\) conservation; no spurious fixed-magnetization
sector or NOD projection is applied.
Both \(\varepsilon=0.02\) and \(5\times10^{-7}\) are tested,
allowing the first gap to vary from order \(10^{-3}\) to
order \(10^{-7}\) while remaining separated from the third
state by an energy of order 0.24--0.26.

## Independent numerical methodology

NumPy constructs the **full Hermitian physical spin Hamiltonian**
from separately implemented tensor products of local spin matrices.
Neither ED2's native Hamiltonian storage nor its Krylov solvers
are used in the independent reference.
The dense reference computes \(E_0,E_1,E_2\), all local
\(\langle S_i^{x,y,z}\rangle\), and all nine Cartesian
two-body correlations for each pair.

Each ED2 run has a fresh input/output directory and opt-in
residual diagnostics. The following are evaluated **separately**:

1. Correct Hamiltonian/spectrum and model properties:
   non-U(1), nonzero avoided gap, isolation of the doublet,
   magnetization reversal, and an interacting/entangled
   edge spin at nonzero bias.
2. ALG=1 and ALG=2 returned eigenpair residual acceptance,
   statuses, numerically ordered energies and native
   wavefunction orthonormality.
3. Independent **individual ground-state** energy, local
   magnetization and two-point correlations.
4. **Low-doublet ensemble** observables computed by reading
   back *both* ED2 native eigenvectors through `re_wf=1`
   and checking their equal-weight average against the
   independently constructed dense doublet projector.
5. The effect of a **verified first-excited initial state**
   on a nominal ground-state Lanczos calculation.

ALG=1 is tested at absolute residual `1e-8`, with additional
`1e-12` runs for the extremely small-gap family.
ALG=2 computes **two** requested low states using both
`NOK/NOM=2/16` and `4/24`, at `1e-8`.

All reported SciPy/NumPy/Fortran toolchain dependencies are
**test-only**; ED2's production dependency set is unchanged.
The tested ED2 binary is built from a separately extracted
source archive (not from a stale build tree).

### Reproducible CI evidence

- [Pure XX edge coupling: OMP=1 and OMP=2](https://github.com/ueda1984/ED2/actions/runs/43687308046).
- [Noncommuting XY edge coupling, OMP=1/2](https://github.com/ueda1984/ED2/actions/runs/43687308058).
- [Shifted-gap minimum and finite-size reference](https://github.com/ueda1984/ED2/actions/runs/43687467478).

**108/108 accepted runs for each coupling family**,
for **216/216 ordinary validation trials total**, across
both OMP profiles. These are separate from the two
explicit excited-state-start counterexamples described below.
Neither run matrix generated a numerical/contract failure or
an accepted but individually unresolved ground observable in
its naturally seeded tests.

The machine-readable plotting grid contains the independent
reference's 32 points:
[cpc-avoided-crossing-reference.csv](cpc-avoided-crossing-reference.csv).
The CI artifacts retain the raw input decks, native
wavefunctions, solver residuals, full observables and
per-trial JSON results. The CSV records OMP=1
reference values rounded in the CI log output.

## Quantitative low-energy structure

For \(\varepsilon=0.02\), the **reference** smallest
two-level gaps and ground-state edge-spin magnetizations are:

| Coupling | Sites | \(\delta\) | \(E_1-E_0\) | \(\langle S_1^z\rangle\) |
| --- | ---: | ---: | ---: | ---: |
| XX | 4 | -0.002 | 0.004346338 | +0.229635 |
| XX | 4 | 0 | 0.003859856 | 0 |
| XX | 4 | +0.002 | 0.004346338 | -0.229635 |
| XX | 6 | -0.002 | 0.004595231 | +0.217278 |
| XX | 6 | 0 | 0.004137931 | 0 |
| XX | 6 | +0.002 | 0.004595231 | -0.217278 |
| XY, \(\eta=0.6\) | 4 | -0.002 | 0.004524011 | +0.254117 |
| XY, \(\eta=0.6\) | 4 | 0 | 0.003904007 | +0.039040 |
| XY, \(\eta=0.6\) | 4 | +0.002 | 0.004242130 | -0.199253 |
| XY, \(\eta=0.6\) | 6 | -0.002 | 0.004764821 | +0.239029 |
| XY, \(\eta=0.6\) | 6 | 0 | 0.004191426 | +0.033683 |
| XY, \(\eta=0.6\) | 6 | +0.002 | 0.004518422 | -0.189649 |

For XX at \(\delta=0\), \([H,S_1^x]=0\) and the
independent ground state has zero edge/core linear entropy
within numerical roundoff. At \(\delta=\pm0.004\),
the calculated edge-spin linear entropies
\(1-\mathrm{Tr}\rho_1^2\) are approximately
\(4.86\times10^{-4}\) (four sites) and
\(3.67\times10^{-4}\) (six sites).

For the more general XY coupling at **\(\delta=0\)**,
the computed linear entropies are already
\(3.59\times10^{-4}\) and \(2.92\times10^{-4}\)
for four and six sites. The relative commutator norms
\(\|[H,S_1^{x,y,z}]\|_F/\|H\|_F\) are nonzero,
which removes the accidental edge-spin conservation of
the pure XX construction.

For the interacting **small-gap** variant
\(\varepsilon=5\times10^{-7}\), the \(\delta=0\) gaps are:

| Edge coupling | Four sites | Six sites |
| --- | ---: | ---: |
| XX | \(9.65334\times10^{-8}\) | \(1.03485\times10^{-7}\) |
| XY, \(\eta=0.6\) | \(9.73459\times10^{-8}\) | \(1.04590\times10^{-7}\) |

The third state stays separated from the lowest doublet
by roughly \(0.24\) or more, so the individual lowest
gap can be much smaller than the isolation of the
two-state cluster.

### The XY gap minimum is shifted, not pinned to \(\delta=0\)

For the XY model with \(\varepsilon=0.02\), a
central three-point parabolic estimate using biases
\(\delta=\{-5\times10^{-4},0,5\times10^{-4}\}\)
locates a **nearby local gap minimum** at approximately

- \(\delta_*=3.10\times10^{-4}\) for four sites,
  with gap around \(0.003892\);
- \(\delta_*=2.87\times10^{-4}\) for six sites,
  with gap around \(0.004182\).

The independent Hellmann--Feynman derivative
\(d(E_1-E_0)/d\delta=
\langle S_1^z\rangle_1-\langle S_1^z\rangle_0\)
is negative at \(\delta=0\) and positive at
\(\delta=5\times10^{-4}\), corroborating
the displaced minimum. **These are parabolic estimates,
not globally optimized or mathematically exact minimizers.**

### Small-coupling two-level consistency

To first order in \(\varepsilon\), let
\(m_x=\langle g_{\mathrm{active}}|S_2^x|g_{\mathrm{active}}\rangle\)
and \(m_y=\langle g_{\mathrm{active}}|S_2^y|g_{\mathrm{active}}\rangle\).
Projecting onto the two decoupled edge-spin orientations
predicts an avoided gap

\[
\Delta_{\mathrm{eff}}(\delta)
=\sqrt{\delta^2+\varepsilon^2
\left(m_x^2+\eta^2m_y^2\right)}
\quad\text{(leading order only)}.
\]

The full physical diagonalization and the independent
effective model are compared programmatically.
Second-order shifts explain why the XY central gap
need not be exactly minimized at zero longitudinal bias.
This approximation is **not substituted** for the exact
dense reference when validating ED2.

## Explicit interacting excited-state-start counterexample

Even in this model with nonzero edge/core interaction,
the direct residual does not prove lowest-state inclusion.

For four sites, \(\varepsilon=0.02,\eta=0.6\),
\(\delta=+0.002\), a validated ALG=2 **first-excited**
native eigenvector is used as the explicit initial
vector for ALG=1. The ALG=1 run reports normal success
at an absolute residual requirement of `1e-12`,
but returns the **first-excited** energy.

- Verified energy above the true ground:
  \(E_1-E_0=0.004242130472\).
- Actual ALG=1 direct residual:
  approximately \(1.54\times10^{-15}\) (OMP=1) and
  below \(10^{-15}\) (OMP=2).
- Returned edge magnetization:
  about \(+0.198199\); the independent unique
  ground-state value is approximately \(-0.199253\).

This is a **mathematically correct excited eigenpair**,
not an ED2 bug. The experiment makes the distinction
between `residual_accepted` and
`global_ground_inclusion_not_assessed` concrete
on a genuinely interacting physical Hamiltonian.

## Interpretation for a CPC paper

The defensible software claim is that ED2
**accurately supports generalized mixed-spin, XYZ/DM/Gamma,
non-U(1) Hamiltonians and their verified observables**,
not that a new avoided-crossing physical phenomenon
was discovered or that a Krylov residual universally
certifies ground-state inclusion.

Possible manuscript figures:

1. \(\Delta_{01}\) and edge \(\langle S_1^z\rangle\)
   versus \(\delta\), for XX and anisotropic XY,
   at \(N=4,6\), with a zoom of the tiny-gap cases.
2. The edge/core linear entropy, showing the
   noncommuting XY central state is entangled while
   the XX model has a conserved edge component at
   the central point.
3. Absolute ground-observable errors and conditional
   low-doublet subspace or trace errors versus the
   achieved direct residual, with the separate
   excited-state-start counterexample annotated.

The primary **fully dense** independent validations have
36 and 216 basis states. The separately completed
[8-/10-site interacting sparse reference study](cpc-interacting-sparse810.md)
extends this very same XY avoided-crossing Hamiltonian
to **1,296 and 7,776 states**, with 52/52 successful
OMP=1/2 trials and additional near-(10^{-7}) gaps.
The largest individual-state differences at tiny gaps are
about (10^{-7}), while the low-doublet trace comparisons
are more stable (below (6	imes10^{-9})).
At tiny gaps, the **independent reference eigenvector**
also has finite residual and must not be treated as
mathematically exact.

The original [8--12-site sparse scaling study](cpc-sparse-scaling.md)
remains a complementary, differently parameterized spin-chain
benchmark. Neither study establishes a new phase diagram,
universal ground-state certification, or general solver speedup.

## Run the two independent tests

From a source archive built into ED2 `source/QS3ED2`:

```sh
python3 tests/mixed_spin_interacting_avoided_crossing.py \
  --exe /absolute/path/to/QS3ED2 \
  --omp 1 \
  --output /new/nonexistent/xx-results

python3 tests/mixed_spin_noncommuting_xy.py \
  --exe /absolute/path/to/QS3ED2 \
  --omp 1 \
  --output /new/nonexistent/xy-results
```

Each output directory must be new. Repeat with
`--omp 2` using distinct destinations.
The full scripts, actual CI logs, native wavefunctions
and independent numerical records are essential to the
reproducibility claims above.
