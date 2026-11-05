# Phased library preparation for Computer Physics Communications

This is a proposal, not an implemented library or a submission-readiness claim.
Read the [baseline audit](baseline-audit.md) and its recorded failures before
changing APIs or numerical behavior. The historical QS-Cube/ED2 repository,
contributors, MIT license, and scientific conventions remain the reference.

## Phase 0 — establish and preserve the executable baseline

Deliver the immutable source revision, compiler/backend/thread information,
reproduction commands, input hashes, independent small-system checks, and the
existing nine workflows. Keep failed cases visible and retain their inputs/logs.
This branch implements this audit and baseline harness only.

Gate: every observed failure has a reproducer and a diagnosis; successful builds
and file-presence smoke tests are not treated as numerical validation. The present
baseline meets the diagnosis goal but has unresolved defects, so it is not a
scientifically validated release baseline for every supported input.

## Phase 1 — build, test, and numerical-interface hygiene

Use separate focused changes for each defect. First add regression assertions,
then repair explicit BLAS/LAPACK interfaces and correct argument types/workspace
queries without changing Hamiltonian construction. Add checked LP64 linkage;
do not silently switch to ILP64 or widen all default integers. Distinguish solver
failure, maximum-iteration termination, empty sector, invalid model, and successful
convergence using statuses while retaining CLI diagnostics.

Handle zero-norm Krylov breakdown and zero-energy convergence only in explicitly
reviewed numerical fixes. Correct the early-check Lanczos branch and sparse-field
indexing with before/after tests. Preserve legacy output interpretation; do not
regenerate historical results merely to pass a test. Investigate uninitialized
TRLan `eps` and all LAPACK `info` paths. Clarify documented symmetry capabilities.

Repair source distribution contents, add numerical targets to Automake/CI, and
separate generated files and binaries from maintained source in a dedicated
repository-hygiene change. Keep Autotools and `QS3ED2` operational. Prefer an
explicit conservative validation profile over silently inheriting `-Ofast`.

Gate: the three currently failing numerical fixtures and the sparse-field
reproducer pass; normal and bounds-checked builds pass, the distribution builds
and runs its tests, and a compiler/backend matrix distinguishes real defects from
backend-specific diagnostic behavior. IEEE-trapping checks must account for
LAPACK's own IEEE arithmetic probes; blanket traps are not a portable test gate.

## Phase 2 — specify ownership and a small Fortran API

Initial consumers: a Fortran research program making many independent ED calls,
a parameter-sweep driver, and the existing CLI. Proposed logical boundaries:

| Boundary | Owns | Contract |
| --- | --- | --- |
| `ed2_kinds`, linear-algebra adapter | real/complex kinds, checked BLAS integer conversion | compile-time precision/backend contract; explicit interfaces |
| model | spins, fields, directed bonds, local excitation caps | validate indices and conventions; immutable after construction |
| sector/basis | excitation interval, permutation generators/characters, ordered representatives and norms | validated sector, deterministic ordering, dimension/overflow checks |
| operator | model + basis reference and transition cache | apply `H` to a vector; cache policy must not change mathematical results |
| solver context/options | Krylov workspace, seed/state, iteration settings | result with eigenvalues, vectors, residuals, convergence/status |
| observables | expectation-value kernels | return arrays; no implicit files or global wavefunction |
| CLI and persistence adapters | namelist parsing, file layout, legacy binary I/O | translate legacy input/output to the API; preserve existing CLI |

Proposed call sequence: construct/validate model → construct sector/basis → create
operator → solve → compute observables → release resources. The exact signatures
need a design review; these names are not a stability promise. Start with a
Fortran module API. `.mod` files and compiler-specific complex/unformatted-file
layouts do not constitute a portable ABI. A C binding or Python wrapper is a
later optional consumer, not a prerequisite or an assumed stable binary ABI.

Specify error ownership, allocation/cleanup rules, instance lifetime, zero-length
outputs, and cancellation before implementing the API. Keep one documented spin,
site/bond indexing, normalization, DM orientation, and output convention. Do not
add a mandatory external numerical dependency or replace the solvers as part of
API extraction.

Gate: a reviewed API specification and two independent-instance use cases, with
an explicit map from every mutable module variable to an owner.

## Phase 3 — extract without changing the mathematics

Move I/O to adapters, initially passing the current arrays explicitly. Then move
persistent model, basis, cache, and Krylov state into owned contexts; make module
visibility private by default with deliberate public entry points. Keep the old
CLI as a thin adapter and compare it with direct library calls on identical input.
Preserve basis ordering, representative normalization, phase factors, excitation
truncation, cache thresholds, Hamiltonian formulas, and solver stopping settings.

Do not combine extraction with vectorization, new reorthogonalization, new
sparse formats, a replacement eigensolver, or a new default tolerance. Expose
legacy behavior first; stage any needed scientific correction separately with
its own reference comparisons. Introduce deterministic seed control as an
explicit option, not a silent change to the executable's current random starts.

Gate: two solver instances run sequentially without leaked state; concurrent
instances are tested before claiming reentrancy/thread safety; the CLI and API
agree in energies, residuals, sector dimensions, and gauge-invariant observables.
Benchmark cache-off/partial/full modes, and validate one and multiple OpenMP
threads with BLAS nesting controlled.

## Phase 4 — broaden scientific validation and packaging

Add small independent tests for mixed spin, nonuniform local caps, bond reversal
with DM sign changes, fields on subsets/permuted sites, nonzero momentum and
reflection sectors, invalid/noncommuting symmetry generators, excluded sectors,
NODmin=NODmax, vacuum/fully polarized states, degeneracies, and exact Krylov
breakdown. Validate cached application against dense multiplication and inner
product Hermiticity. Check eigenpair residual norms, normalization, orthogonality,
and invariant subspaces for degeneracies rather than comparing eigenvectors
component by component.

Test all three solvers and restart convergence; exercise exhaustion/error exits,
wavefunction round trips (including >9 vectors and tiny sectors), and observable
index/order conventions. Use absolute-plus-relative tolerances with justification
and tighter operator-level checks than end-to-end stopping-error tolerances.
Never assume agreement with a historical output proves physical correctness.

Package an installable static library first if sufficient, export only reviewed
modules, and provide an external consumer compile/link example. Test clean
installation, `make distcheck`, and execution from the source archive. Document
compiler/module compatibility and BLAS LP64 requirements. Benchmark release
flags separately from validation flags; retain unoptimized checked runs.

Gate: clean downstream build/link/run, clean source distribution, documented
resource requirements, and a validation matrix with no unexplained required
failures. Performance claims require controlled, repeated measurements and error
bars; concurrent cloud audit timings are not publication benchmarks.

## Phase 5 — prepare the CPC software and reproducibility package

Prepare a manuscript explaining the physics problem, actual symmetry/truncation
capabilities, basis and operator representation, all solver choices and limitations,
validation against independent calculations, and bounded performance/resource
claims. Separate scientific capabilities from library engineering improvements;
assess novelty against relevant exact-diagonalization software and the historical
upstream code without claiming novelty from modularization alone.

Prepare a program summary, supported platforms, dependencies, licensing and
attribution, installation and consumer examples, input/output specification,
algorithm description, test-data provenance, and a complete reproducibility
bundle. Include exact code revision, build commands, dependency versions,
compiler flags, hardware, threads, random initialization policy, input datasets,
and checksums for each reported result. Supply scripts to regenerate tables and
figures plus archived results and known limitations.

Resolve current version/citation inconsistencies deliberately with maintainers.
Select release/version/DOI/archive steps only with authorization; retain upstream
attribution and contributors. Check the journal's current author and software/data
requirements at submission time: this plan is not a verified checklist of current
CPC policy and makes no assertion of acceptance or DOI assignment.

Gate: an independent fresh-environment reproduction, reviewed manuscript and
software archive, and explicit maintainer authorization for release/submission.
No phase authorizes merging this branch into `main`.
