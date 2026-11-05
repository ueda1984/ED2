# Phase 1: explicit BLAS/LAPACK interfaces

Base commit: `bfe96f0754c662d40d0207e7dc3d1ec6863486dc`.
Feature branch: `codex/ed2-phase1-lapack-interfaces`.
Date: 2026-11-05 (Asia/Tokyo).

This change repairs numerical-library calling contracts. It is a subset of
[phase 1 of the CPC plan](cpc-library-plan.md), not a solver modernization or a
repair of the [three baseline numerical failures](baseline-audit.md).

## Changes and preserved behavior

`source/blas_lapack_interfaces.f90` declares native Fortran interfaces for the ten
external entry points actually used: ZHEEVR, DSTEVX, DLAMCH, DZNRM2, ZDOTC, ZGEMM,
ZCOPY, ZAXPY, ZSCAL and ZDSCAL. Assumed-size arrays retain normal Fortran sequence
association for vector elements, columns and flattened matrices. ZDOTC retains
the native Fortran complex-result ABI; no C binding, wrapper library or dependency
is introduced. The interfaces explicitly select real64 and 32-bit BLAS integers,
matching the existing LP64 build; alternate default-kind/ILP64 builds are not
silently accommodated.

Both solver modules import these interfaces instead of locally declaring external
function result types. All production BLAS/LAPACK callers compile with
`-Werror=implicit-interface`. The Automake source list and module dependency order
include the new module, including parallel builds.

Dense ZHEEVR now receives real(8) VL/VU and one-element **arrays** for all three
workspace-query outputs. TRLan also uses rank-one query buffers and its actual
ISUPPZ buffer. Queries are checked for nonzero INFO before sizes are read. Query
lengths must be finite, meet documented minima, and fit LP64 integers; fractional
real recommendations round upward before allocation. NaN/Inf rejection uses IEEE
binary64 exponent bits, because GNU `-Ofast` can eliminate IEEE_IS_FINITE checks.
The default fast build's negative tests exposed this during development, and
all final profiles exercise the replacement guard. Binary64 representation is an
explicit platform assumption of this validation, not a claim of non-IEEE support.

TRLan's initial DSTEVX call uses the already initialized `abstol = 2*DLAMCH('S')`
instead of the uninitialized local EPS. This defines a previously undefined
LAPACK argument using the value computed immediately before the call. It does
not change the namelist energy-convergence threshold, restart criteria or any
Krylov recurrence. Existing ZHEEVR absolute tolerances remain unchanged.

Every ZHEEVR and DSTEVX query/solve checks INFO before consuming its results.
Failures now report the routine/context and INFO to stderr and terminate with a
nonzero status, instead of printing and potentially continuing/writing output or
using a bare STOP that can exit successfully. Successful CLI input/output formats
are unchanged; error diagnostics/status intentionally become reliable.

Hamiltonian formulas, symmetry permutations/phases and normalization, basis and
site/bond indexing, truncation/cache semantics, eigenvalue algorithms, random
initialization, observables, convergence policies and reference datasets are
unchanged. Neither `ham2vec.f90`, `state_lists.f90`, `input_param.f90`, `main.f90`
nor the historical reference data is modified. No compatibility flag is needed
for any validated build.

## Tests and reproduction

The machine/compiler/backend match the original audit: Debian 13 x86_64,
GNU Fortran 14.2.0, Autoconf 2.72, Automake 1.17, GNU Make 4.4.1, and
OpenBLAS 0.3.29 pthread with LAPACK. NumPy 2.3.5 is used only by existing independent
validation and the new historical-output comparator; it remains a test-only
optional installation requirement, not a mandatory ED2 dependency.

Use a **fresh source tree per profile**, because example programs append some
observable outputs. From the committed feature branch, for example:

```sh
set -eu
# Cloud-specific tool activation; use your installed toolchain elsewhere.
. /workspace/.ed2-env/activate.sh
repo=/workspace/ED2
stage=/workspace/ed2-phase1-reproduction  # Must not already exist.
mkdir "$stage"
git -C "$repo" archive HEAD | tar -x -C "$stage"
cd "$stage"
autoreconf -vfi
# Conservative profile:
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
make check-lapack
make check
make check-long
python3 "$repo/tests/compare_historical.py" --build "$stage" \
  --output /workspace/phase1-historical-report.json
```

For the default profile use `./configure FC=gfortran` (the unchanged default
`-fopenmp -Ofast`). For the bounds profile use:

```sh
./configure FC=gfortran \
  FCFLAGS='-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra'
```

Neither profile adds `-fallow-argument-mismatch`. Bounds checks use `-O2` with
runtime checks, rather than the baseline audit's `-Og` setting. Blanket floating
point traps are not used: the baseline standalone probe demonstrated that the
linked ZHEEVR can itself trap during valid calls. Historical runs use one OpenMP
thread, while the independent suite covers one and two threads with BLAS nesting
controlled.

Run the unchanged independent suite against each profile:

```sh
python3 /workspace/ED2/tests/numerical_baseline.py \
  --exe /workspace/ed2-phase1-reproduction/source/QS3ED2 \
  --output /workspace/phase1-independent-results
```

This command intentionally returns **1**, with 18 passes and these three failures
still visible: `early-check-lanczos`, `zero-alg1`, `zero-alg2`. Do not describe the
21-case suite as fully passing or change its assertions/reference answers. The
new INFO handling can change how a failing zero-Hamiltonian case exits; that is
not a repair of Krylov breakdown.

`make check-lapack` runs 11 checks without NumPy: complex Hermitian ZHEEVR
spectra/residuals, eigenvalues-only/value-range ZHEEVR, full/index-selected DSTEVX,
workspace rounding, eight rejection/error modes, and two dense-caller injected
failures. The first positive group counts as one check. Its test-only ZHEEVR stub
is linked only into `lapack_failure_check`, never into QS3ED2. Query failure leaves
length buffers undefined deliberately; the caller must stop before reading them.
The solve-failure check verifies that no eigenvalues file is written. The target
is separate from historical `make check`, preserving its one-test Automake count.

`tests/compare_historical.py` checks all 27 full-array historical comparisons,
including output shapes, exact index columns, finite numerical columns, energies
at `atol=rtol=1e-9` and magnetization/XYZ correlations at `atol=rtol=1e-6`. It records
input hashes, dimensions, energies and errors, and refuses to overwrite an existing
report. It must be run after check-long in a fresh tree: it does not independently
prove output freshness for arbitrarily reused directories.

The original `tools/cpc_baseline.py` and baseline evidence are retained unchanged
as historical audit tools/records; that runner still requests the old compatibility
flag. Use the commands above for validation of this feature.

## Results and remaining risks

See [the machine-readable validation record](phase1-lapack/validation.json) for
source hashes, flags, counts, executable hashes and the per-case failures, and
[the evidence index](phase1-lapack/README.md) for selected logs.

The default, conservative and bounds-checked profiles build successfully, pass
11 interface/error-path checks and the standard regression (TOTAL=1, PASS=1), and
each independent suite has 18 passes/3 diagnosed failures. The conservative and
bounds profiles additionally run the full nine-example regression and historical
comparisons: all nine examples and all 27 comparisons pass in each profile. The default
fast profile is not used for historical accuracy claims or an extended-suite claim.

Residual risks/limits:

- The early-check Lanczos return and both zero-Hamiltonian failures remain;
  sparse-field indexing, RNG control and other baseline defects are out of scope.
- DLAMCH-based tolerance initialization corrects undefined TRLan input; only the
  tested cases establish numerical agreement, not every spectrum/backend.
- LP64 and native Fortran complex ABI compatibility must match the linked library.
  Intel/MKL, other compilers, ILP64 and non-IEEE platforms were not validated.
- INFO checks and workspace guards provide process-level failures, not yet a
  recoverable library status API. Allocation exhaustion is still handled by the
  Fortran runtime; global dimension/product overflow outside the guarded queries
  is not fixed here.
- Existing unrelated compiler warnings remain. Stub routines intentionally have
  unused dummy arguments; absence of argument-mismatch errors is not a claim
  that the whole repository is warning-free or numerically correct.
- The source-distribution omissions recorded in the baseline audit are not repaired
  by this change. Do not claim a complete distcheck/release validation.
- Concurrent cloud run timings are not controlled performance benchmarks.

Development-only runs before final validation were superseded after strengthening
NaN/Inf guards. A staging command also fell through after a failed directory setup;
its job was stopped, and only generated checkout changes were restored. Final
validation uses new external trees with stop-on-error commands; generated checkout
files/build products are not part of the commit.
