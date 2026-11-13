# Solver accuracy and reproducibility

!!! warning "Development-branch controls"
    This page documents Phase 1-F/1-G features of `ueda1984/ED2`
    through commit `d720a3538d32ceadb598c956fae7564fb5b85460`.
    They were not in the upstream `QS-Cube/ED2` `main` at the
    time this page was prepared. Check your actual executable SHA
    before using them; the hosted `latest` site may document an
    older source.

## Available solvers

- `ALG=1`: Lanczos ground-candidate calculation.
- `ALG=2`: restarted Hermitian solver for several requested pairs.
- `ALG=3`: dense LAPACK diagonalization for small sectors.

For ALG=1/2, the legacy `lnc_ene_conv` primarily controls an
**energy-change stopping criterion**. It is not itself an
absolute error bound on individual observables or a guarantee
that the lowest eigenspace was found.

## Additional opt-in environment controls

The following settings are **environment variables**, not Fortran
NAMELIST keys.

| Variable | Meaning |
| --- | --- |
| `ED2_RNG_SEED` | Integer seed for iterative random starts |
| `ED2_INITIAL_VECTOR` | Path to an explicit complex starting vector |
| `ED2_DIAGNOSTICS` | Path for a new diagnostic text sidecar |
| `ED2_RESIDUAL_ATOL` | Absolute residual tolerance, energy units |
| `ED2_RESIDUAL_RTOL` | Dimensionless relative residual tolerance |
| `ED2_RESIDUAL_SCALE` | Caller-supplied scale `B` for relative tolerance |
| `ED2_ENERGY_ORIGIN` | Caller-supplied energy origin `c`, default zero |

Without the new settings, the existing random initialization and
solver stopping behavior remain unchanged.

An explicit starting-vector file has exactly `n` lines with one
`real imaginary` pair per line. Here `n` is the constructed
**symmetry-sector dimension**, not necessarily the full spin Hilbert
space dimension. A zero/nonfinite vector is rejected. ALG=2 can
generate additional random directions even with an explicit start.

Seeding provides replayability under the same compiler/runtime,
BLAS and OpenMP configuration. It does not promise bitwise results
across different toolchains or thread counts.

`ED2_DIAGNOSTICS` must name a **nonexistent** file. Existing files
are not overwritten. The sidecar is currently an internal, evolving
diagnostic format rather than a stable library API.

## Residual acceptance

For each normalized requested pair `(lambda_i,psi_i)`, the code
computes the **direct** residual

$$
r_i=\|H\psi_i-\lambda_i\psi_i\|_2.
$$

With opt-in residual controls, successful ALG=1/2 termination also
requires

$$
r_i \le r_{\mathrm{atol}}+r_{\mathrm{rtol}} S_i,
\qquad
S_i=\max\{B,|\lambda_i-c|\}.
$$

`B` and `c` are explicit **user choices**; `B` is **not**
certified by ED2 as an operator-norm bound. Positive relative
tolerance requires explicitly specifying `ED2_RESIDUAL_SCALE`.
For absolute-only requests the scale is unnecessary. At least
one tolerance must be positive.

This is an **additional necessary success condition**. It does
not replace the legacy energy gate or change the default tolerance;
there is no silent epsilon floor or automatic threshold relaxation.
ALG=2 tests **all requested pairs**. Unattainably strict settings
can lead to budget exhaustion or numerical failure, not success.

### Example: absolute threshold and diagnostic sidecar

After building the development executable and generating the
`examples/chain/` inputs:

```bash
cd examples/chain
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
ED2_RNG_SEED=19 ED2_DIAGNOSTICS=solver-diagnostics.txt \
ED2_RESIDUAL_ATOL=1e-10 \
../../source/QS3ED2 < input.dat > output.dat
```

Use an empty output directory and unique sidecar path for each run.
This threshold is a **syntax example**, not a universally recommended
scientific accuracy setting.

A relative-only request needs both `ED2_RESIDUAL_RTOL` and
`ED2_RESIDUAL_SCALE`. Setting `ED2_ENERGY_ORIGIN` changes
the **tolerance scale** but does not shift the Hamiltonian's
actual matrix-vector action or prevent numerical cancellation.

## Diagnostic status and limits

The sidecar records reason and status, Ritz/Rayleigh values, direct
residuals, thresholds, iteration and Hamiltonian-action counts,
seed information, and (if requested) per-pair residual acceptance.

| Status | Meaning |
| --- | --- |
| 0 | Existing acceptance and optional added condition satisfied |
| 1 | Iteration/restart budget exhausted |
| 2 | Numerical failure or breakdown |
| -1 | Invalid input or parameters |

The `lowest_state_inclusion` field remains unassessed without
independent evidence. The observed `operator_scale` is an
**estimate from computed actions**, not a certified spectral norm.

## Accuracy and cost: why the distinction matters

Phase 1-G reports a matched-start matrix calculation where ALG=1's
residual improved from approximately `6.1e-6` to `9.6e-11`
under stricter acceptance, while H applications increased
from `122` to `2355`. These are workload-specific observations,
not universal performance guarantees.

Its tiny-gap two-spin Heisenberg example retains an individual
magnetization error of `0.5` despite a residual around
`4.2e-16`. Small residuals alone do not establish ground-state
inclusion, near-degenerate eigenvector fidelity or the accuracy
of physical observables.

### A concrete counterexample: an accurate excited eigenpair is not the ground state

The [CPC near-degeneracy validation](../development/cpc-near-degenerate-gaps.md)
uses an actual mixed-spin XYZ/DM/Gamma Hamiltonian with a
controlled, nonzero lowest gap \(\Delta=10^{-6}\).

A validated first-excited eigenvector, supplied using
`ED2_INITIAL_VECTOR`, makes ALG=1 return `status=0`
with direct residual around \(10^{-15}\), although the returned
energy exceeds the independently verified ground energy by
\(10^{-6}\). The returned edge-spin magnetization is \(+1/2\),
while the unique ground state has \(-1/2\).

For a smaller gap \(\Delta=10^{-10}\) and an ALG=1 random
start, `ED2_RESIDUAL_ATOL=1e-8` also permits a low-doublet
state whose individual ground-state magnetization is markedly
wrong. Tightening the residual to `1e-12` fixes **those specific
runs**, but does not solve the excited-start counterexample.

This is the standard mathematical distinction between
**eigenpair residual**, **membership in a low-energy cluster**,
and **identification of the lowest individual eigenvector**.
At exact degeneracy, a state-specific magnetization is not
unique at all; use cluster projectors, traces, or invariant
quantities.

The current `lowest_state_inclusion` diagnostic remains
`not_assessed` for precisely this reason. No output flag should
be interpreted as a universal ground-state certificate.

For published calculations record the full Git SHA, all model inputs,
compiler/BLAS versions, OpenMP thread settings, initial state/seed,
the requested and achieved tolerances, and unsuccessful runs.

## Development validation

- [Phase 1-E precision contract](../development/phase1e-precision-contract.md)
- [Phase 1-E validation](../development/phase1e-validation.md)
- [Phase 1-F diagnostics](../development/phase1f-diagnostics-reproducibility.md)
- [Phase 1-G residual policy](../development/phase1g-optin-residual-policy.md)
- [CPC physical near-degenerate and excited-start tests](../development/cpc-near-degenerate-gaps.md)
