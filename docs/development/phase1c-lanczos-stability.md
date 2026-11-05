# Phase 1-C: Lanczos breakdown, convergence and eigenpair certification

Base: `2050208999940044c155cce7fb60d42b7228481d`.
Branch: `codex/ed2-phase1c-lanczos-stability`.
Date: 2026-11-05 (UTC).

Read the [baseline audit](baseline-audit.md), [Phase 1-A](phase1-lapack-interfaces.md),
[Phase 1-B](phase1b-sparse-field-indexing.md) and [CPC plan](cpc-library-plan.md).
This phase intentionally repairs numerical solver behavior. It does not change
the Hamiltonian, symmetry, field indexing, LAPACK interfaces or historical answers.

## Before changing the solvers

The Phase 1-B bounds executable was checked against its committed executable
SHA-256; both solver source files matched the base. A fresh run of the unchanged
`tests/numerical_baseline.py` reproduced 18 passes and the three known failures:

- `early-check-lanczos`: dimension 16 returned its first Rayleigh quotient
  (about -0.163965 in this run) instead of E0 ≈ -1.1509634025303355;
- `zero-alg1`: nonfinite arithmetic/quality diagnostics;
- `zero-alg2`: nonfinite arithmetic/outputs.

The original RNG is not fixed; the wrong early-return value is an observation,
not a deterministic reference. Pre-fix inputs are reproducible with the unchanged
baseline generator; representative stdout and the old report's failure details
are retained in [the evidence directory](phase1c-lanczos/README.md).

## Numerical decisions and explicit stopping policy

Let `sigma` be the largest observed norm of H acting on a normalized basis vector,
and epsilon be double-precision machine epsilon. Two separate tolerances are used:

- Near-breakdown threshold: `64*epsilon*max(sigma,tiny)`.
- Actual eigenpair acceptance threshold:
  `max(256*epsilon,sqrt(lnc_ene_conv))*max(sigma,abs(lambda),tiny)`.

The second threshold applies to the **directly recomputed** `||H*psi-lambda*psi||2`
for a normalized output vector, not just a projected residual estimate or a
relative-variance expression. The square-root relation is motivated by the usual
second-order energy error versus vector error; it is not a universal eigenvalue
error bound without spectral-gap assumptions. Independent tests also evaluate
residuals using their own Hamiltonian matrices and exact spectral norms.

Away from zero (`abs(Enew)>sqrt(epsilon)*max(sigma,tiny)`), the energy-change
criterion remains exactly the legacy `abs(Eold/Enew-1) < lnc_ene_conv`. Near zero
it is explicitly extended to `abs(Eold-Enew)/max(sigma,tiny) < lnc_ene_conv`.
This is a scale-relative absolute-energy test at/near zero, not relative accuracy
with respect to an eigenvalue that vanishes. No unseen previous energy is used
for the first comparison. The tolerance is not silently loosened when an
iteration budget is exhausted.

Ordinary ALG=1 energy convergence still obeys `minitr` and `itrint`. Ordinary
ALG=2 convergence still requires **two consecutive** energy-converged checks of
the retained `NOK` Ritz values, now comparing actual successive Ritz values rather
than an initial basis diagonal. Both ordinary exits additionally require actual
residual certification. Full-space completion or an exact/numerically closed
invariant subspace can exit without two energy checks, because a direct residual
certificate is available and another Krylov step can be undefined.

Termination classes remain distinct in kernel status/reason and CLI diagnostics:

| Condition | Behavior |
| --- | --- |
| Exactly zero beta / exactly closed computed subspace | Happy breakdown only after eigenpair residual certification |
| Positive beta below the roundoff-scaled threshold | Numerical near-breakdown; never divide by it; certify the pair or report numerical failure |
| Ordinary energy and residual criteria satisfied | Converged |
| Full finite-dimensional space and residuals certified | Completed |
| Budget exhausted without a permitted certified exit | Status 1 / CLI nonzero exit, no result files |
| Nonfinite action, invalid reconstructed vector or uncertified numerical breakdown | Status 2 / CLI nonzero exit |
| Invalid parameters or zero/nonfinite supplied initial vector | Status -1 / CLI nonzero exit |

A happy breakdown certifies eigenpairs **inside the starting invariant/Krylov
subspace**, not global extremality over an unvisited orthogonal complement.
In particular, ALG=1 started in an exact excited eigenvector returns that eigenpair;
it cannot infer an excluded ground state from residual zero. This behavior is
explicitly tested. Normal CLI starts remain generic random starts. No finite set
of residual checks proves global optimality if the starting subspace misses a
lower invariant block.

## ALG=1: guarded three-term recurrence and two-pass reconstruction

The conventional method remains a three-term Lanczos recurrence with O(n) vector
storage, tridiagonal DSTEVX and a second pass for the Ritz vector. It does not
start storing an n-by-iteration Krylov basis or introduce full reorthogonalization.

- Iteration `k==1` is no longer treated as a one-dimensional physical sector.
  A one-element projected eigensystem has coefficient 1; reconstruction multiplies
  the actual normalized start vector rather than replacing its first component.
- Beta is classified before normalizing the next vector. No normalization occurs
  after exact/near breakdown or after the final iteration.
- The second pass replays only previously accepted, non-breakdown beta values.
  Its work arrays are separate from the continuing recurrence, so a failed
  residual certificate can continue iteration instead of corrupting the first pass.
- Every successful result is normalized and checked with a fresh H action. The
  iteration count is initialized on every path, including 1D and invalid inputs.
- At `maxitr`, an unconverged candidate is not returned as a successful result.

The quality report now obtains `<H^2>` from `||H*psi||2^2` for a Hermitian H,
calculating the old ratio as `(mean(H)/||H*psi||2)^2` to avoid squaring large/small
energies before division. For exact `H*psi=0` the displayed legacy ratio is defined
as zero by convention. The actual residual norm is reported independently and
is authoritative. Near-zero-energy output warns that the relative-variance ratio
is ill-conditioned; it is not used for convergence.

## ALG=2: stabilized block thick restart

This is an explicit extension of the old single-start thick restart, not merely
an IF statement around beta division. A single-vector Krylov space contains at
most one independent direction from each exact eigenspace; guarding division
alone cannot recover multiplicities or multiple requested vectors of H=0.

The internal kernel now:

1. Initializes a block of `NOE` independent generic starts with two-pass modified
   Gram-Schmidt. An explicitly supplied test start is an additional direction
   when capacity permits, so an exact excited start does not replace one of the
   generic directions needed to resolve a degenerate ground space.
2. Expands a Hermitian block Krylov space using stored H actions and two-pass
   orthogonalization. Exact/near closures are logged separately; a deterministic
   coordinate-complement direction continues exploration instead of dividing by
   a zero/tiny residual or duplicating an eigenvector.
3. Forms the actual Hermitian Rayleigh-Ritz projection `Q^H H Q`, solved by the
   existing checked ZHEEVR interface. This avoids fragile hand-maintained complex
   arrowhead couplings and handles multiple closed blocks. It diagonalizes only
   the bounded projected problem, not a new dense n-by-n global Hamiltonian.
4. Retains the lowest `NOK` Ritz vectors, reorthogonalizes them and recomputes their
   H actions for the next expansion. All requested output pairs get direct residual
   checks; a closed subspace is tested as a whole, not inferred from one small beta.

`NOM`/`NOK` workspace sizes are capped by the sector dimension, but requested `NOE`
is never silently reduced. A 1x1 projected problem is solved as the scalar it is,
without a LAPACK workspace query. There must be expansion room unless the full
sector fits. ALG=2's `maxitr` is a budget of Rayleigh-Ritz cycles (including the
initial block solve); reaching it does not authorize publication of the last
unconverged iterate. The printed total step count now records actual H applications
instead of the old `restarts*NOM` estimate. ALG=1 still reports recurrence steps,
excluding reconstruction/certification actions.

Two-pass orthogonalization, stored H images, exact projected matrices, extra generic
starts for explicit initial vectors and residual checks change work/memory and
iteration histories. These are documented stability choices. `lnc_ene_conv` remains
an energy-change control with the stated zero-energy extension, not a newly
renamed residual tolerance. No performance improvement is claimed.

## Implementation boundaries

`source/krylov_solvers.f90` contains internal, callback-based numerical kernels
and checked normalization. This is a test seam, not a promised stable library API.
`lanczos.f90` and `eigen_solver.f90` retain namelist/CLI adapters and process-level
error handling. Module-level private H adapters avoid GNU executable-stack
trampolines; test matrices live only in the test driver. The existing dense
ALG=3 implementation and the BLAS/LAPACK interface file are unchanged.

`ham2vec.f90`, `state_lists.f90`, `input_param.f90`, `main.f90`, field indexing,
all reference datasets and existing independent/sparse-field/LAPACK assertions
are unchanged. Only necessary Automake module/link dependencies and an optional
`check-lanczos` target were added. There are no new mandatory dependencies.

## Validation and reproduction

Use a new external tree for each profile, because historical examples append some
observable files. GNU Fortran 14.2.0 / Debian 13 x86_64, existing LP64 OpenBLAS
0.3.29 pthread/LAPACK and test-only NumPy 2.3.5 were used. Example:

```sh
set -eu
. /workspace/.ed2-env/activate.sh  # Cloud activation; use installed tools elsewhere.
repo=/workspace/ED2
stage=/workspace/ed2-phase1c-reproduction  # Must not already exist.
mkdir "$stage"
git -C "$repo" archive HEAD | tar -x -C "$stage"
cd "$stage"
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
make check-lanczos LANCZOS_OUTPUT=/workspace/phase1c-kernels
make check-sparse-field SPARSE_FIELD_OUTPUT=/workspace/phase1c-fields
make check-lapack
make check
make check-long
python3 "$repo/tests/numerical_baseline.py" --exe "$stage/source/QS3ED2" \
  --output /workspace/phase1c-independent
python3 "$repo/tests/compare_historical.py" --build "$stage" \
  --output /workspace/phase1c-historical.json
```

GNU defaults use `./configure FC=gfortran` (`-fopenmp -Ofast`); bounds use
`FCFLAGS='-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra'`. No profile uses
`-fallow-argument-mismatch`. The unchanged 21-case suite now returns **0**, with
21 actual successes and no failed/expected-failure/skipped cases.

`check-lanczos` runs 82 internal-kernel checks and 26 actual CLI checks, each at
one/two OpenMP threads and single-threaded BLAS. The internal checks include 16
intended invalid-input/nonconvergence controls; the CLI checks include four
intended insufficient-budget failures. A check passes only if its specified
status/exit behavior is observed. These controls do not mask the formerly failing
21-case suite, which must genuinely pass.

Coverage includes zero/scalar operators, 1D sectors, exact ground/excited starts,
positive near-breakdown beta, zero-energy states, scaled operators (1e-9/1e9),
block/restart degeneracy, complete degenerate-cluster projectors, orthonormality,
closed-block termination at minimal budget, zero/nonfinite starts, and invalid or
insufficient iteration budgets. Test matrices/spin operators are independently
constructed. Complete degenerate clusters are compared via spectral projectors;
partial scalar eigenspaces use membership, residual and orthonormality checks,
not component-wise eigenvector comparisons.

The CLI tests read the actual GNU native unformatted wavefunction records and
recompute residuals independently, covering reconstruction and the CLI adapter,
not only the kernel callback. The reader is explicitly a native-format validation
tool, not a new persistence format. Nonconverged CLI cases must exit nonzero and
create neither eigenvalue nor wavefunction result files.

All three final profiles also pass the unchanged 70 sparse-field CLI cases,
10 operator cases, 28 record permutations, 11 LAPACK checks, one standard regression,
nine extended examples and all 27 historical comparisons. Eleven focused ALG=1
kernel/CLI controls additionally passed with traps for invalid/zero-division/
overflow. Broad trap-enabled ZHEEVR testing is not claimed: its legitimate internal
IEEE probes were already diagnosed in Phase 1-A.

The [evidence index](phase1c-lanczos/README.md) and
[validation record](phase1c-lanczos/validation.json) include source/compiler/library
hashes, exact flags, statuses, residuals, orthogonality/projector errors and
historical errors. Observed maximum normalized residuals were about 4.00e-15 for
the kernel checks and 2.06e-14 for CLI wavefunctions; maximum degenerate projector
error was about 3.01e-15. Historical maximum absolute differences across profiles
were 2.28e-13 (energies), 8.22e-12 (magnetization) and 1.67e-10 (XYZ correlations).
These are observed test errors, not universal accuracy bounds.

## Remaining limitations and review risks

- This is a numerical-method change for ALG=2; independent review of block
  initialization, projection, restart and convergence semantics remains valuable.
- ALG=1 retains the memory-saving three-term method without full reorthogonalization.
  Difficult clustered/large problems can lose orthogonality or fail certification;
  they now report failure rather than silently accepting an invalid pair. The
  certificate does not prove global ground-state inclusion for an excluded start.
- ALG=2 also assumes the sampled block has adequate overlap with wanted invariant
  spaces. Random-block rank is generic, not a deterministic proof of global
  extremality for every problem. Finite budgets can legitimately fail.
- Stored basis/H images/projected Ritz vectors increase ALG=2 peak memory and
  additional H actions/orthogonalization increase work; ALG=1 replay/certification
  adds work and O(n) buffers. No controlled scalability/performance benchmark was run.
- Energy-change thresholds below floating-point resolution cannot guarantee
  meaningful relative energy accuracy. The residual roundoff floor and explicit
  near-zero energy policy must be considered when choosing tolerances.
- IEEE binary64 / native complex / LP64 assumptions remain; other compilers,
  MKL/ILP64, extreme subnormal/overflow ranges and exhaustive spectra were not
  validated. Checks at scales 1e-9 and 1e9 are not a proof for all exponents.
- The new 1D handling is for ALG=1/2. ALG=3's caller/legacy workspace-query guards
  and unrelated wavefunction readback/packaging/global-state issues were not repaired.
  No public reentrancy or stable ABI guarantee is introduced by the internal test seam.
- Production wrappers still use process-level errors; a recoverable user-facing
  library error API and explicit owned contexts remain later CPC plan phases.
- Historical output/log differences in iteration counts and diagnostics are
  intentional; input/output file formats, Hamiltonian conventions and reference
  answers are unchanged. Concurrent cloud timings are not publication benchmarks.
