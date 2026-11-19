# Cross-library ED benchmark: protocol (pilot, 2026-11-19)

**Status:** experimental protocol; no external-library timing or superiority claim is made.
This work is separated from the CPC draft and from the unmerged cubic-specific ED2
kernels. Do not merge or copy benchmark numbers into the paper without review.

## Scientific question

Measure (a) performance of the general-purpose public ED2 operator/basis,
(b) performance of validated generic cache/CSR options, and (c) the performance
ceiling of explicitly *specialized* dilute-magnon algorithms. Compare with
QuSpin, XDiag, and optionally lattice-symmetries and HΦ, each at a pinned
source version. The three tiers must not be conflated:

| Tier | ED2 | External solvers | Permitted claims |
| --- | --- | --- | --- |
| G: generic, matched U(1) sector | Official upstream or pinned fork default | Pinned QuSpin/XDiag general fixed-magnetization block | End-to-end / setup / matvec / eigenpair on identical physical subspace |
| T: optimized but general | Mainline coefficient bound, CSR and memory controls if physically applicable | Competitors' best documented generic U(1) implementations (including 1D-specialized general library APIs) | Best generic implementation per model, clearly reporting configuration |
| S: explicitly specialized | Unmerged 3-magnon cubic Gamma/Oh development branches, with independent physics proof | Competitors' best *equivalent symmetry sector* if supported; otherwise do not assert direct speedup | Specialization advantage **within ED2**; cross-library feasibility or time-to-equivalent-physics, qualified |

A competitor running only U(1) must NOT be described as using the same
Hamiltonian block as ED2 in a Gamma/Oh A1g block. The U(1)-only dimension is
not the sector dimension after space-group projection. Numerical agreement
of the lowest energy of distinct blocks must not be presumed.

## Canonical interoperable physical workload (G and T)

Spin-1/2 periodic nearest-neighbor isotropic Heisenberg model:

    H = sum_(i,j in undirected positive-axis NN bonds) S_i dot S_j

Set J=1 and all fields zero. Work at exactly m down spins, m=2 or 3.
This is a true conserved U(1) sector, not a truncated nonconserving model.
For bond (i,j), the Sz product contributes +1/4 (parallel) or -1/4
(antiparallel), and each opposite-spin exchange contributes +1/2 to the
spin-flipped state. NO arbitrary constant-energy shift is permitted.
The site order is x-fastest:

    i = x + Lx*y + Lx*Ly*z

Generate one undirected positive-axis bond per site per axis; every periodic
axis must have length >=3 to avoid double counting short periodic bonds.
For QuSpin use pauli=0, fixed Nup=N-m; for XDiag use Spinhalf(N,N-m).
Use explicit matching physical bond lists, not the pre-existing anisotropic
cubic_sp_HB example or its energy.

The standalone reference.py prints a physical SHA256 derived from the shape,
m, precise bond list, site order, J and field. Each runner must report that
identifier, the actual U(1) dimension C(N,m), and an energy / direct residual.
A different hash means different physics and disqualifies a comparison.

### Primary 1000-site-class benchmark grid (supersedes the earlier small-site pilot)

**Study variable:** spatial site count N approximately 10^3, with *few*
down spins. Do **not** replace this with 12--64 site performance graphs.
Small ED instances only validate the physical model, symmetry projectors
and measured code paths; they are not the claimed benchmark target.

| Main workload | Geometry | N | m | U(1)-only dimension | Mandatory symmetry comparisons |
| --- | --- | ---: | ---: | ---: | --- |
| 1D-1000 | PBC chain 1000 | 1000 | 2 / 3 | 499500 / 166167000 | U(1), Gamma, Gamma + inversion-even |
| 2D-1024 | PBC square 32 x 32 | 1024 | 2 / 3 | 523776 / 178433024 | U(1), Gamma, Gamma + C4v A1 |
| 3D-1000 | PBC cube 10 x 10 x 10 | 1000 | 2 / 3 | 499500 / 166167000 | U(1), Gamma, Gamma + Oh A1g |

**Primary scientific target:** N=1000, m=3, cubic Gamma/Oh A1g,
including the experimentally validated ED2 specialized opt-ins, on
**the exact same physical Hamiltonian**. A 1000-site chain and
1024-site square form complementary geometry/overhead controls.
m=2 is the lower-memory U(1) control; m=4 or different irreps are
optional challenges conditional on measured memory/time feasibility.

For 1000-site cubic m=3, dimension is **166167000** at U(1) only,
**166167** at Gamma using just translation symmetries, and **4299**
at Gamma/Oh A1g for the validated cubic sector. The latter is the
development PR's dimension for its *compatible isotropic/XXZ*
space-group representation; equal Hilbert dimension does not alone
prove equality of Hamiltonians.

The first shared comparison *at the stated N* is a **feasibility matrix**
for each library/setting:
(i) whether 1000/1024-site spin representations are accepted,
(ii) whether correct fixed-m sectors are constructible,
(iii) whether Gamma and full point-group A1 projections are
supported and agree with independent dimension checks,
(iv) whether the target eigenpair converges within the declared
resource budget. Only **after** those checks is a runtime ratio
meaningful.

For Gamma/point-group benchmarks construct the *identical projected
block* in each library and independently check dimensions and energy.
Noncommuting space-group generators may need complete group projectors,
not merely a blind combination of commuting-generator APIs.
A solver limited to U(1) must **not** be directly compared to
ED2's 4299-dimensional A1g block as a speed victory. Report its
U(1)-only ability and resource status, and additionally compare against
ED2 U(1) on the same physical Hamiltonian if feasible.

QuSpin supports custom large-integer representations above 64 sites,
but that is not evidence that it can efficiently enumerate m=3 at
N=1000; the current U(1) pilot's cap intentionally refuses 166 million
states. QuSpin's docs specifically discuss large bit types and
noncommuting symmetries. XDiag's documentation lists Spinhalf and
generic lattice symmetries but does not itself establish a 1000-site
successful benchmark: test and record exact version, compile-time /
runtime feasibility and errors before drawing conclusions.

**Memory and time guardrails:** do not attempt a 166/178 million
state CSR construction by default, nor infer resource needs solely
from the final reduced block dimension (basis enumeration and orbit
canonicalization can still visit the full fixed-m sector). Limit
memory / wall time for feasibility tests. Use a dedicated approved HPC
allocation for any large attempt; failed builds, unsupported N or
resource exhaustion are distinct statuses and *not* measured
speedups. The run_quspin.py prototype is U(1)-only and cannot
benchmark the m=3 projected target until a separately verified
large-site space-group adapter exists.

### Independent validation (not performance data)

Dense-oracle prechecks at N=6, 12, 16, 36 and 64 when manageable;
tiny cases do not appear as headline runtime performance.
Run a 1000-site metadata/symmetry-dimension check before any solver
allocation, and save its physical input hash. Verify all large-case
Hamiltonians are truly U(1) conserving and exactly obey specified
space-group symmetries.

## External library source-level feasibility audit (2026-11-19)

These are *code-path observations*, not executed N=1000 benchmarks.

- **QuSpin 1.0.1:** the documented spin_basis_general constructor supports
  large custom integer dtypes including uint1024, and custom translation
  maps. Its make(N_p=0) option avoids allocating an unnecessarily large
  prefix table for dilute bases. These features do not establish practical
  runtime, peak RSS, or correctness for 1000/1024-site m=3. The current
  run_quspin.py probes exactly U(1) or commuting-translation Gamma;
  noncommuting point-group symmetry projections require separate proof.
- **XDiag 0.6.0:** source xdiag/blocks/spinhalf.cpp (at upstream commit
  27065537a743c642c35efa0e1b5ae4f326fd984c) dispatches fixed-nup
  cases with nsites>128 to Combinations<BitsetDynamic>. Thus it is false
  to claim a universal 64-site limit for fixed magnetization. However
  its fast <k>sublattice symmetry backend explicitly rejects nsites>64.
  xdiag/blocks/distributed/spinhalf_distributed.cpp independently rejects
  nsites>64 and distributed permutation symmetries. Whether the generic
  BitsetDynamic+BasisSymmetric implementation finishes N1000/m3 efficiently
  is a resource/benchmark question, not resolved by the dispatch alone.
- **XDiag source probe:** xdiag_preflight.cpp is a preliminary single-node
  C++17 U(1)-only capacity/eigenvalue adapter, based on documented
  Spinhalf(n,nup), OpSum/SdotS and eigval0 APIs. The code has not yet been
  compiled against the external XDiag library. Its optional solve path
  records the eigenvalue but does **not** compute the direct residual;
  consequently it is NOT eligible for a CPC performance table.
- **Scope warning:** an U(1)-only success at 1000 sites is not the same
  calculation as a Gamma+Oh A1g success (dimension 4299). Separately report
  XDiag generic symmetries, sublattice optimized symmetries and MPI support.

Sources:
- https://quspin.github.io/QuSpin/generated/quspin.basis.spin_basis_general.html
- https://quspin.github.io/QuSpin/basis.html
- https://github.com/awietek/xdiag/blob/27065537a743c642c35efa0e1b5ae4f326fd984c/xdiag/blocks/spinhalf.cpp
- https://github.com/awietek/xdiag/blob/27065537a743c642c35efa0e1b5ae4f326fd984c/xdiag/blocks/distributed/spinhalf_distributed.cpp

### XDiag full cubic space-group API: a specific source-level memory bottleneck

The pinned XDiag 0.6.0 `PermutationGroup` implementation explicitly stores
an `arma::Mat<int64_t>` of **all group multiplication products**
(`G x G`) and another `nsites x G` int64 matrix of site permutations.
For a cubic (L^3=N=1000) site system, Gamma/Oh A1g requires a
space group of (G=1000\times48=48000) elements, giving at least:

- Multiplication table: **18,432,000,000 bytes (17.17 GiB)**.
- Site-permutation matrix: **384,000,000 bytes (366.2 MiB)**.
- Inverse indices: 384,000 bytes.
- Sum of these **allocated payload arrays alone**:
  **18,816,384,000 bytes (17.52 GiB)**, excluding temporary allocations,
  Hilbert-space basis, Lanczos and runtime overhead.

Moreover, the current `compute_multiply` implementation loops over
every ordered pair of group elements and searches for the product in the
permutation list, so the construction cost is a substantial additional
scaling concern. `RepresentationSet` accepts only one SitePermutation
group type; putting translations and Oh into two independent groups is
not a supported substitute in this native API.

This is a **structural bound for the present full-group API** from
verified source, **not measured XDiag RSS or a proof that all potential
XDiag extensions/algorithms are intrinsically unsuitable**. In
particular, QuSpin's generator factorization stores small generators
rather than this multiplication table. The comparison therefore
highlights a valid research question: whether a library exploits the
semidirect product / few-magnon structure or materializes the group.

The independent `group_memory_floor.py --self-test` automates the
accounting without constructing any group or calling XDiag.

Source:
- https://github.com/awietek/xdiag/blob/27065537a743c642c35efa0e1b5ae4f326fd984c/xdiag/symmetries/permutation_group.cpp
- https://github.com/awietek/xdiag/blob/27065537a743c642c35efa0e1b5ae4f326fd984c/xdiag/symmetries/representation_set.cpp

## Specialized cubic lane (S) is a distinct physical workload

Existing ED2 draft PRs #13--#21 contain exact but highly scoped cubic
three-magnon algorithms (anchored-orbit basis, local-bond cache, direct-hop
pair table, coordinate caching and canonical pruning). The direct pair table
requires a validated 3-magnon S=1/2 Gamma / one-dimensional Oh
representation and a specific compatible XXZ bond/model certificate.
This is not a generic fast path. Do NOT infer that its runtime applies
to anisotropic DM, nonconserved Ndown, arbitrary momentum or arbitrary
point-group representations.

Reuse the **exact source commit and archived physical input** from that
experimental stack, preserving every model parameter, symmetry and SHA.
This J1-J2/XXZ comparison is not interchangeable with the canonical
nearest-neighbor Heisenberg workload above. Document each baseline and
fast opt-in separately, and run full equality checks on basis rank/norm,
Hamiltonian cache (when possible), Ritz E and direct residual.

The recorded N1000 PR #14 vs PR #16 pilot samples were 1.546690 vs
0.364112 s at one thread, 0.590361 vs 0.296340 s at four, and
0.417922 vs 0.275544 s at eight in A1g. They were from separate
allocations, without complete publication-level repeat / independent
node-exclusivity attestation. They must NOT be used as publication-quality
speedup figures or cross-library results.

## Exact source-equivalent cubic J1-J2 XXZ model: primary Tier S

**This is the actual model in ED2 cubic PR16**, distinct from the
nearest-neighbor isotropic Heisenberg control:
`H = sum_b [ Jxy(b) * (Sx_i Sx_j + Sy_i Sy_j) + Jzz(b) * Sz_i Sz_j ]`.

| Shell | Displacements per site (positive representatives) | Bonds at N1000 | Jxy | Jzz |
| --- | --- | ---: | ---: | ---: |
| J1 axial | 3 | 3000 | 1.0 | 0.77 |
| J2 face diagonal | 6 | 6000 | 0.31 | 0.19 |

This 9N-bond Hamiltonian conserves Ndown and is invariant under the
cubic full point group. All-up energy at N=1000 equals
`E_all_up = (3000*.77+6000*.19)/4 = 862.5` exactly in the mathematical
model. For three magnons, the U(1), Gamma and Gamma/Oh A1g dimensions
are 166167000, 166167 and 4299, respectively. The existing ED2 PR16
PRESTO A1g run reported Ritz energy about 847.307894065573,
but that is a prior experimental ED2 result, **not** a new external
comparison or an independent large-N spectrum.

- `export_cubic_xxz.py` copies the exact bond-direction list and
  real XXZ constants from PR16 pinned tests/cpc_cubic_pointgroup.py.
  It writes 0-based, 4-column `i j Jxy Jzz` bonds and a JSON manifest
  recording source, all-up energy, bond-file SHA256 and Hamiltonian
  physical SHA256. The axial/face bond counts, all-up energy and
  invariance under independent translational/O_h generators are tested.
- `validate_native_xxz_input.py` checks **only bonds** of an
  original ED2 11-column `bonds.dat` against those expected couplings.
  It does NOT establish physical equivalence of the spin table, fields,
  NOD sector, eigenstate irrep, solver settings or bytewise input files.
- `xxz_projected_oracle.py` independently builds the full 4^3-site,
  three-magnon space, explicitly enumerates all 44 Gamma/Oh A1g
  orbits and computes the normalized projected Hamiltonian.
  Local verification yielded E0=40.345908939159166,
  direct eigen-residual about 3.35e-14 and matrix asymmetry below 1e-15.
  This is a correctness-only oracle, not an N1000 performance datum.
- `run_quspin.py --model cubic-j1j2-xxz` now builds exactly the
  Jzz SzSz and Jxy/2 S+S- plus S-S+ terms, with `pauli=0`.
  N1000/m3 remains resource-gated by raw C(1000,3), not by its
  small post-projection dimension.
- `xdiag_preflight.cpp` accepts a 4-column XXZ bond file as
  its optional final argument `cubic-j1j2-xxz`.
  Its source-level implementation uses `Jxy*Exchange+Jzz*SzSz`,
  consistent with XDiag operator definitions. XDiag compilation,
  actual runtime and direct residual **remain unverified**.

Example canonical inputs (no solver allocation):

~~~bash
python3 benchmarks/cross_library/export_cubic_xxz.py --self-test
python3 benchmarks/cross_library/xxz_projected_oracle.py --self-test
python3 benchmarks/cross_library/export_cubic_xxz.py --side 10 --ndown 3 --output-prefix /tmp/cubic_j1j2_n1000
python3 benchmarks/cross_library/run_quspin.py --shape 10,10,10 --ndown 3 --model cubic-j1j2-xxz --sector gamma --mode dimensions
# Source-only XDiag demonstration after installation and separately verified build:
./xdiag_preflight 1000 3 dimension 1000000 /tmp/cubic_j1j2_n1000.xxz_bonds cubic-j1j2-xxz
~~~

Verify native historical ED2 physical bonds with
`python3 benchmarks/cross_library/validate_native_xxz_input.py --side 10
--ed2-bonds /path/to/original/bonds.dat` (one command on one line).
Make sure inputs/symmetry/energy are also independently verified
before any speed ratio is attempted.

Reference: https://github.com/ueda1984/ED2/blob/9692d8b1fd903bcb7b1b67a5f180c57430702510/tests/cpc_cubic_pointgroup.py
and https://awietek.github.io/xdiag/examples/entanglement_ground_state/
(XDiag `Exchange` + `SzSz` example).

## Independently recovered 1000-site cubic physical eigenvalue

A **standalone C++17** implementation (independent_cubic_a1g.cpp)
constructs the complete Gamma/Oh A1g projected Hamiltonian of the
certified cubic J1-J2 XXZ model, without calling any ED2 Fortran,
QuSpin or XDiag routines. It enumerates anchored three-site shapes
(C(N-1,2), not C(N,3)), closes each spatial orbit under all signed
axis permutations and three anchor translations, and uses full raw
orbit cardinalities for the correctly normalized symmetric basis.

The code asserts both (i) orbit-size sum equals C(N,3) and
(ii) max_abs(H-H^T)<=1e-9. It emits a Matrix Market sparse matrix;
independent_cubic_solver.py then solves for the smallest eigenvalue
with SciPy eigsh and checks the direct eigenvector residual.

For **N=1000** the independent local calculation yielded:
- 166167000 total fixed-three-down states
- 4299 Gamma/Oh A1g projected basis states
- 193791 stored matrix elements (including diagonal)
- eigenvalue E0=847.3078940655738 (local sample),
  direct residual=6.04e-10 at eigsh tol=1e-12,
  matrix asymmetry <=4.45e-16
- previous ED2 PRESTO PR16 E0=847.3078940655731, within 1e-12
  (floating-point environment dependent).

At side=4, D=44, standalone C++ construction returned
E0=40.34590893916067, consistent within 2e-12 with the independent
fully enumerated 41664-state Python-orbit-projected Hamiltonian
E0=40.345908939159166. This double independent check makes a
model/normalization mistake much less likely.

Build and verify (all calculations are *physics checks*, not timing
benchmarks; requires a C++17 compiler, NumPy and SciPy):

~~~bash
g++ -std=c++17 -O2 -Wall -Wextra benchmarks/cross_library/independent_cubic_a1g.cpp -o /tmp/independent_cubic_a1g
/tmp/independent_cubic_a1g 4 /tmp/independent_A1g_N64.mtx
OPENBLAS_NUM_THREADS=1 python3 benchmarks/cross_library/independent_cubic_solver.py --side 4 --matrix /tmp/independent_A1g_N64.mtx
/tmp/independent_cubic_a1g 10 /tmp/independent_A1g_N1000.mtx
OPENBLAS_NUM_THREADS=1 python3 benchmarks/cross_library/independent_cubic_solver.py --side 10 --matrix /tmp/independent_A1g_N1000.mtx
~~~

An independent GitHub workflow compiled and verified both sides
successfully on the experiment branch:
https://github.com/ueda1984/ED2/actions/runs/44468545819

**Limitations:** This is an independent code implementation of the
same mathematical orbit-reduction concept, not a competing software
benchmark. It does not prove QuSpin/XDiag symmetry compatibility
or timing, and cannot replace independent scientific and
exclusive-resource review of PRESTO performance claims.

## Required measurements for every individual trial

- Exact repo URL, commit, compiler/version/flags, library + BLAS version,
  job/allocation/host CPU model, NUMA/affinity, OpenMP and BLAS thread
  settings, input and binary checksums, physical SHA.
- Model, shape, m, exact conserved symmetries, irrep quantum numbers,
  full dimension, projected dimension (if used), and matrix representation
  (matrix-free, padded, CSR, cached, direct CSR).
- Separate basis construction, H setup, single H*vector (warmup and
  repeats), Lanczos/eigensolver, and **basis-to-converged-eigenpair**
  wall times; whole-process wall / peak RSS. Do not label cache-array
  payload bytes as peak RSS.
- E0, explicit direct residual norm ||H psi-E psi|| for normalized psi,
  actual iterations and H-application count, target convergence settings.
  Internal relative/squared criteria are NOT substitutes for the direct
  residual.
- Status: PASS, PHYSICS_MISMATCH, UNSUPPORTED_SECTOR,
  RESOURCE_LIMIT, TIMED_OUT, BUILD_FAILED or NOT_RUN.
  Resource refusals are not zero runtimes and not proof that a library
  categorically lacks the capability.

For final measurements use dedicated/independently attested resources,
fixed CPU affinity, interleaved configurations, warmup and >=5 repeated
measurements. Publish per-trial raw JSON/CSV, median, spread and preferably
a confidence interval. Never conflate software initialization, Python
import time, matrix setup, H application, and eigensolver iteration.
Differences between SciPy ARPACK and ED2 Lanczos must be stated;
both end-to-end and per-matvec numbers are required.

### Correctness gate

1. Verify physical SHA, exact sector definitions and dimensions.
2. Independently compare small (12 or fewer sites, bounded dimension)
   spectra against reference.py; require |Delta E0| <=1e-8 and a direct
   residual norm <=1e-8 for the target eigenvector.
3. For different space-group representations require an independent
   symmetry projector / rank check and a physical eigenpair match.
4. Only after the above run timed cases and record resource failures.
5. Independent review of PRESTO allocation exclusivity and archival
   provenance is required before any paper performance claim.

## Tools currently provided

- reference.py: NumPy dense small-sector oracle (NO performance inference),
  bond generator and model SHA256. Use --metadata-only at N1000.
- symmetry_dimensions.py: independent Burnside group-orbit dimension
  calculations at large N, with N1000/1024 dimensions verified locally.
- verify_orbits.py: independent explicit small-geometry orbit traversal.
  Six small cases were locally verified against Burnside, including
  (4,4,4), three magnons: 44 space-group A1 orbits.
- export_case.py: deterministic same-physics .bonds and .json manifest,
  with physical model and bond SHA256 for the other language adapters.
- run_quspin.py: gated U(1)-only or translation-Gamma QuSpin prototype,
  staged --mode dimensions|basis|solve. Default dimensions stage
  requires NO QuSpin and allocates no many-body states; basis/solve
  require an installed QuSpin. The source has not yet been executed with
  QuSpin (only NumPy preflight checks are covered by CI).
- xdiag_preflight.cpp: XDiag 0.6.0 API-based preliminary C++ U(1)
  constructor and energy probe, for Heisenberg NN or exact cubic J1-J2 XXZ.
  Requires an external XDiag build; not yet compiler-tested against it.
- .github/workflows/cross-library-science.yml: independent model,
  orbit-count, export, guard and syntax checks; does not install or run
  QuSpin/XDiag, and does not produce performance claims.

### Non-allocating feasibility commands

~~~bash
python3 benchmarks/cross_library/reference.py --self-test
OPENBLAS_NUM_THREADS=1 python3 benchmarks/cross_library/symmetry_dimensions.py --self-test
OPENBLAS_NUM_THREADS=1 python3 benchmarks/cross_library/verify_orbits.py --self-test
python3 benchmarks/cross_library/reference.py --shape 10,10,10 --ndown 3 --metadata-only
python3 benchmarks/cross_library/export_case.py --shape 10,10,10 --ndown 3 --output-prefix /tmp/cubic_n1000_m3
OPENBLAS_NUM_THREADS=1 python3 benchmarks/cross_library/run_quspin.py --shape 10,10,10 --ndown 3 --sector gamma --mode dimensions
# Intentionally SKIPPED_RESOURCE_CAP; never construct C(1000,3) by default:
OPENBLAS_NUM_THREADS=1 python3 benchmarks/cross_library/run_quspin.py --shape 10,10,10 --ndown 3 --sector gamma --mode basis --max-raw-dim 1000000
~~~

### Controlled, conditional probe on an installed external library

QuSpin 1000-site **two-magnon** basis-only evaluation may be attempted,
only with CPU/memory/time limits approved by the host operator:

~~~bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 120s \
  python3 benchmarks/cross_library/run_quspin.py --shape 10,10,10 --ndown 2 --sector gamma --mode basis --max-raw-dim 600000
~~~

The default m3 safety cap is not a claim that QuSpin cannot run m3.
Record skipped resource gates distinctly from runtime failures.

After the user separately builds xdiag_preflight.cpp against XDiag v0.6.0
and pins the installed source SHA, the equivalent source-only U(1) probe
takes arguments: N, m, dimension|basis|solve, maximum raw U(1)
dimension, canonical .bonds path. For example:

~~~bash
./xdiag_preflight 1000 3 dimension 1000000 /tmp/cubic_n1000_m3.bonds
./xdiag_preflight 1000 3 basis 1000000 /tmp/cubic_n1000_m3.bonds
~~~

The latter must return SKIPPED_RESOURCE_CAP without constructing
C(1000,3). Actual 1000-site XDiag timing remains unmeasured.

## Specific audit warning about existing documentation

The docs/examples/cubic_sp_HB.md input describes generic anisotropies and
transverse fields while fixing NOD=3. A fixed-NOD window is not a
conserved subspace under such terms. Separately, vp/wk_dim_scaling.md
discusses a Heisenberg benchmark 'except observables'. Do not assume the
two narratives define the identical Hamiltonian without auditing their
exact input decks. A projected-Hamiltonian calculation must be explicitly
labeled as such and must never be compared as an unrestricted eigenproblem.

## Upstream documentation / existing development

- QuSpin general basis: https://quspin.github.io/QuSpin/generated/quspin.basis.spin_basis_general.html
- QuSpin large integer support: https://quspin.github.io/QuSpin/basis.html
- XDiag documentation / release: https://awietek.github.io/xdiag/
- XDiag benchmarks: https://awietek.github.io/xdiag/benchmarks/
- ED2 cubic development PR16: https://github.com/ueda1984/ED2/pull/16
- ED2 mainline CPC draft PR35: https://github.com/ueda1984/ED2/pull/35
