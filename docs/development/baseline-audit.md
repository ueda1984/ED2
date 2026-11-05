# ED2 executable and library-readiness audit

Audit date: 2026-11-05 (Asia/Tokyo). Branch: `codex/cpc-library-preparation`.
Starting revision: `ca812ca0ca402bed2525556bd797cec5e0a002e8`.
The parent scientific implementation is `49941dc4b437a8fb6c2b5c8dac7d621866731d7c`;
the starting branch adds only `AGENTS.md`. Those instructions were read before
the audit. No files under `source/`, numerical algorithms, reference datasets,
existing build declarations, or existing regression scripts were changed.

This is a working development baseline with diagnosed defects, not a claim that
all scientific inputs are validated or that a reusable library already exists.
See the [phased CPC plan](cpc-library-plan.md).

## Reproduce

Prerequisites: Linux x86_64, GNU Fortran 14.2.0 (other compilers not validated in
this audit), GCC/binutils, GNU Make 4.4.1, Autoconf 2.72, Automake 1.17, GNU m4,
BLAS/LAPACK, Python 3.11 or newer and the test-only dependency in
`tests/requirements-baseline.txt` (NumPy 2.3.5). NumPy is not an ED2 runtime or
build dependency. The measured cloud machine uses Debian 13 and OpenBLAS 0.3.29
pthread; its CPU is reported as Intel Xeon Platinum 8573C with five visible CPUs.
The snapshot's local toolchain and trusted Debian packages are retained outside
the checkout under `/workspace/.ed2-env`; no credentials are required by tests.

```sh
# In the prepared cloud environment, activate relocated tools and library paths:
. /workspace/.ed2-env/activate.sh
cd /workspace/ED2
# On another machine install the system prerequisites and test-only NumPy first.
# Use a NEW output directory outside the checkout; do not reuse stale artifacts.
python3 tools/cpc_baseline.py --output /workspace/ed2-reproduction --long --debug
```

The runner archives the exact local `HEAD` with `git archive`, uses plain external
staging directories (not Git worktrees), and preserves the checkout. It records
source revision, archive/harness/executable hashes, compiler/tool versions, CPU,
linked libraries, flags, thread settings, commands, exit statuses, and results.
Uncommitted solver edits are intentionally NOT part of the archived source;
working-tree status and the separately executed harness hashes are recorded.
Every generated small-system input and stdout log is retained. The command
returns nonzero for the known failures; it does not convert them to passes.

The documented build sequence was tried literally in a clean archive:

```sh
autoreconf -vfi
./configure FC=gfortran
make -j2
make check
```

`autoreconf` and configure passed; compilation failed (exit 2) on incompatible
implicit LAPACK arguments. Therefore `make check` could not run for that default
profile. The tested compatibility profile is:

```sh
./configure FC=gfortran FCFLAGS='-O2 -fopenmp -fallow-argument-mismatch'
make clean
make -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 make check
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 make check-long
```

The flag is a documented compatibility workaround, not a repair of the LAPACK
interfaces. `-O2` deliberately avoids the repository default `-Ofast` for numerical
validation. The default and conservative flags were both recorded; numerical
results are from the conservative build. A second build uses `-Og -g -fopenmp
-fallow-argument-mismatch -fcheck=all -fbacktrace -Wall -Wextra`.

A blanket `-ffpe-trap=invalid,zero,overflow` exploratory build caused failures
inside ZHEEVR even for an independent correctly typed 16x16 diagonal matrix.
That standalone probe succeeds without traps. Such failures are not evidence
that ED2's dense spectrum is wrong; the portable checked baseline uses bounds
checks without blanket traps. Separately, trapping a one-dimensional Lanczos
case points to division by zero in `lanczos.f90:88`, before the special-case return.

## Architecture and numerical contracts

The implementation comprises 2,414 Fortran lines in six modules and one program.
The build's dependency order is broadly correct but semantic coupling is wider
than explicit `use` declarations: kernels import `input_param` inside routines.

| File/module | Current responsibility and coupling |
| --- | --- |
| `input_param` | namelists and disk I/O; mutable model, sector, basis, wavefunction, energies, phase factors and observable lists |
| `state_lists` | bounded-occupation combinatorics, ranking/unranking, representative search, normalization, sorted basis; consumes global sector/permutation state |
| `ham2vec` | transition cache, dense Hamiltonian, cached/on-demand H-vector application, local/correlation kernels; consumes model/basis globals |
| `lanczos` | conventional three-vector Lanczos, eigenvector reconstruction, module-level iteration parameters/dimension |
| `eigen_solver` | thick-restarted Lanczos and dense diagonalization, module-level solver parameters/dimension |
| `get_expectation_values` | symmetry-equivalent-site aggregation, observable output and component ordering |
| `main` | read → basis → solver selection → wavefunction/energy I/O → ground-state observables |

Almost all module entities are public by default. Repeated solves cannot safely
be assumed independent: persistent allocatables and module dimensions/settings
have no complete reset/destroy API. Fixed units, stdin namelist reads, filenames
and `stop` paths tie computations to a process. OpenMP loops are present, but
OpenMP within one solve does not establish concurrent-instance reentrancy.

Local spins are represented by bounded excitation occupations, with total
`NODmin:NODmax`; omitted spin sectors mean a projected Hamiltonian, not generally
a conserved magnetization block. Combinatorial ranks, sorting, local ladder
factors, phase signs, and representative-normalization ratios are scientific
contracts. Up to six supplied site permutations and character phases are used
for sector projection, assuming appropriate finite commuting symmetry operations.
The code does not establish their validity, mutual compatibility or invariance of
the supplied Hamiltonian. Geometrical bond inversion in the cubic example is not
spin inversion; the generic documentation's spin-inversion capability is not
supported by an identified spin-flip implementation in this audit.

`read_ip` sorts each bond's sites and reverses DM signs when needed. The nine bond
couplings are XYZ, DM xyz, and symmetric off-diagonal Gamma xyz. Fields enter with
a **positive** `h·S` convention. Local excitation caps, basis ordering and the
nonstandard text correlation order (xx,yx,zx,xy,yy,zy,xz,yz,zz) must be preserved.
The independent fixtures construct ordinary spin matrices to validate these signs
and observables, rather than copying the Hamiltonian implementation.

The operator has three substantially duplicated formula paths: full dense
construction, precomputed transitions, and the remainder evaluated during H-vector
application. `wk_dim` controls cached rows, not physics; `MNTE` bounds slots.
With default flags, default integers occupy four bytes and complex(8) sixteen;
cache storage is approximately `(sizeof(integer)+sizeof(complex))*slots*cached_rows`,
plus basis and temporary storage. Dense eigenvectors/matrices are quadratic in
sector dimension. Off-diagonal contributions below `1e-15` are omitted. Moving
that threshold or changing summation order is not merely an API change.

ALG=1 uses a three-term Lanczos recurrence without full reorthogonalization and
reconstructs the eigenvector in a second pass; convergence compares successive
energies, not an explicit residual criterion. ALG=2 retains a subspace, performs
reorthogonalization, and uses `DSTEVX`, `ZHEEVR` and `ZGEMM`; convergence again uses
relative energy changes. ALG=3 builds the dense Hermitian operator and calls
`ZHEEVR`. BLAS complex norms/dot products/AXPY/scaling/copies are shared. Existing
namelist random initialization calls `random_seed` without a user-controlled
seed; iterative starts, convergence histories and timings are not bitwise
reproducible. Baseline assertions compare toleranced spectra and gauge-invariant
observables, not raw vectors or exact iteration counts.

## Observed test coverage and outcomes

The compatibility build passed `make check` (TOTAL=1, PASS=1, FAIL=0, ERROR=0,
SKIP=0) and `make check-long` (all nine scripts reported PASS). All nine historical
examples select ALG=1. They do not cover dense/TRLan, breakdown, invalid inputs,
nonzero momenta comprehensively, or repeated library instances.

The so-called strict harness only checks nonempty output files. The relaxed
harness checks the first token of `eigenvalues.dat`, which is the **index**, not
the energy; neither compares against reference results or rejects all NaN cases.
The new runner compares eigenvalues at `atol=rtol=1e-9`, and local magnetization
and XYZ correlations at `atol=rtol=1e-6` against the historical files, checking
shapes, indices, finite values and current-run input hashes. These looser
observable tolerances allow solver stopping error and random-start differences;
independent small-system observable assertions use `atol=rtol=1e-9`.

`tests/numerical_baseline.py` runs 21 black-box cases, individually timed out at
120 seconds, with fresh directories and no solver modification:

- Ten general complex XYZ+DM+Gamma/field cases: dense, thick-restarted, and
  conventional Lanczos with full/off/partial cache, at one and two OpenMP threads.
  All pass independent dense spectra and ground-state observables at 1e-9.
- The Heisenberg dimer: full spectrum `[-0.75,0.25,0.25,0.25]`, zero local
  magnetization and diagonal singlet correlations `-0.25`; passes.
- Zero Hamiltonian: dense passes; both iterative solvers fail finite-result checks.
- Vacuum, one-dimensional conventional Lanczos: spectrum passes in the normal
  and bounds-checked profiles; floating-point trapping exposes an intermediate
  zero-norm division, so spectrum agreement alone does not validate arithmetic.
- Six-site ring, exactly two down spins, momentum 1 and 3, one/two threads:
  projected dense spectra/dimensions pass independent translation-projector
  calculations (complex phase included). Sector observables are not checked.
- A nonconserved excitation-truncated general Hamiltonian: projected dense
  spectrum and ground-state observables pass.
- Early convergence checking (`minitr=1,itrint=1`) for conventional Lanczos:
  fails the independent ground-energy comparison.

Both normal and bounds-checked 21-case runs have **18 passes and 3 failures**.
Those failures are left visible with nonzero harness exit status. They are not
marked expected-failure or silently excluded. Multithreaded checks here cover
small fixtures, not the full nine-example suite or all scheduling patterns.

## Findings requiring separate fixes

| Priority/evidence | Finding | Consequence and next step |
| --- | --- | --- |
| High; default build failed | Implicit LAPACK calls mix scalar/array work arguments; dense `vl,vu` are declared INTEGER although ZHEEVR expects DOUBLE PRECISION (`eigen_solver.f90:270,305`) | Add explicit, tested interfaces with correct argument kinds/ranks and workspace queries; the compatibility flag masks diagnostics, not ABI errors |
| High; numerical fixture failed | `lanczos.f90:95` tests iteration `i==1` as if the physical dimension were one | With `minitr=1,itrint=1`, a dimension-16 problem returns its random first Rayleigh quotient; validate the real dimension and preserve a failing regression before fixing |
| High; two numerical fixtures failed | No robust exact-breakdown guard before normalization, relative energy convergence divides by energy, and quality diagnostic divides by `<H²>` | Zero Hamiltonian causes nonfinite outputs/diagnostics in ALG=1/2; dimension-one arithmetic can also divide by zero even if the final energy is correct |
| High; bounds-checked targeted reproducer | Field arrays have shape `(3,NO_one)` but Hamiltonian paths address `hvec(:,p_one(j))` rather than the field record (`ham2vec.f90:56,268,468`) | `NO_one=1` with its field at site 2 raises an index-2/upper-bound-1 error in dense construction; all standard fixtures use full site-ordered field records and miss this |
| High; static source inspection | TRLan declares `eps` and passes it as `DSTEVX` ABSTOL (`eigen_solver.f90:130`) without assigning it; a separate `abstol` is initialized | Normal fixtures pass, but this does not validate the uninitialized argument; correct and test in a separate numerical-interface change |
| Medium; generated archive tested | `make dist` succeeds but omits `tests/test_common_strict.sh`, `script_*/run.sh` and model data, and `input/input.dat` | Fresh distribution configure/build succeed under compatibility flags; `make check` fails, TOTAL=1/PASS=0/FAIL=1, because the strict helper is missing (test exit 127) |
| Medium; static source inspection | Default-integer combinatorial ranks/dimensions and unchecked products can overflow before the negative-THS check | Current tested dimensions fit 32-bit integers; the diagnostic claiming an eight-byte limit is misleading. Plan checked dimension arithmetic and explicit BLAS integer limits |
| Medium; static source inspection | I/O errors mostly rely on runtime/STOP, several STOP paths can return success, dense LAPACK failure prints without a structured error, and wavefunctions use compiler-specific unformatted records and `I1` vector numbering | Library callers cannot recover reliably; define status/cleanup and persistence contracts while preserving the CLI format |
| Medium; coverage limitation | Numerical tests cover only a subset of spin/symmetry/truncation space; all historical examples use ALG=1 | Mixed-spin independent reference, reflection validation, repeated contexts, residual/orthogonality and restart edge cases remain necessary |
| Medium; static/source packaging evidence | Generated Makefiles/configuration, caches, historical `.mod` and `a.out` files are tracked; `setup.sh` hardcodes gfortran+MKL although configure only accepts MKL with ifort | Use pristine source staging for this audit; handle source hygiene and setup defaults separately |
| Release planning; metadata inspection | `configure.ac` identifies 1.0.0, tracked generated Makefile identifies 2.0.0, CITATION.cff identifies 0.1.0 | Maintain attribution; resolve version/release provenance deliberately rather than changing citation metadata during an audit |

No finding above was repaired in the numerical implementation. Compiler warnings
about potentially uninitialized indices and an apparent subscript-zero reference
also require review; they are **warnings**, not independently demonstrated failures
in every reported path. In particular, the subscript-zero warning is inside an
`i>1` branch. Do not count all warnings as confirmed runtime bugs.

## Recorded evidence and limits

The [baseline evidence directory](baseline/README.md) contains machine records and selected logs/inputs.
Captured log files have trailing whitespace and excess blank lines at EOF trimmed
for repository review; input decks, numerical JSON records and reported hashes
are unchanged. Raw logs remain in the external run directories. The full workflow record includes an exploratory
trapping debug run (6 passes/15 failures, many backend-induced), and a separately
repeated bounds-only workflow supplies the valid checked counts (18/3). These
records intentionally retain the diagnostic change; the current reproducer's
`--debug` profile is bounds-only. Both workflow invocations returned exit 1 because
the three real numerical failures remain. The exact committed harness's final
small-system reruns and fresh-output refusal are recorded separately.

The nine historical spectra, local-magnetization arrays and XYZ-correlation arrays
all match the reference files: 27 comparisons passed. In this run, maximum absolute
errors are approximately 4.55e-13 for energies, 1.58e-11 for local magnetization and
1.93e-10 for XYZ correlations. These are observed errors, not universally proven
bounds. The chain baseline has THS=166751, sector dimension 1669 and
E0≈-12.99620173300453; the 1000-site cubic inversion example has THS=166167000,
sector dimension 23719 and E0≈736.7224213680655. Full values and input hashes are
recorded in the JSON evidence.

Historical outputs are corroborating regressions, not independent physical truth.
The small NumPy reference is independent operator construction but shares the
installed BLAS/LAPACK ecosystem; the analytic singlet spectrum/observables and
simple ring-sector values provide additional checks. Dense projected sectors
cover momenta 1 and 3 with dimension 2 and spectra [-0.8,1.2] and [-0.3,0.7],
respectively. Degenerate zero-Hamiltonian eigenvectors/observables are intentionally
not compared. Checked compilation is not a memory-safety proof.

No full second-compiler/MKL/ILP64 matrix, exhaustive symmetry audit, bitwise RNG
reproducibility, public API, installable library, benchmark scaling study, or
fresh-cloud restoration/publication has been validated. Concurrent exploratory
runs make recorded timings unsuitable for CPC performance claims. `make distcheck`
as a complete target was not run; the generated archive was extracted, configured,
built and tested directly, exposing its test-packaging failure.
