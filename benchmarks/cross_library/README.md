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

### Proposed pilot grid

| Geometry | Shape | Down spins | U(1) dimension (formula) | Purpose |
| --- | --- | --- | --- | --- |
| PBC chain | (12,) and (24,) | 2, 3 | C(N,m) | Generic pilot, check 1D-specific optimized modes |
| PBC square | (4,4) and (6,6) | 2, 3 | C(N,m) | Two-dimensional translation-free baseline |
| PBC cubic | (4,4,4) | 2, 3 | C(64,m) | Local scaling and early memory behavior |
| PBC cubic | (6,6,6) | 2, 3 | C(216,m) | Moderate size: enforce memory budget |
| PBC cubic | (10,10,10) | 3 | C(1000,3)=166167000 | **Not** a default allocation: only with an approved memory/time budget and appropriate spatial-symmetry approach |

For the first stage compare *U(1) only* (no translations, reflection, or
Oh projection) to ensure same Hilbert-space block. In a second stage use
matched translational sectors (Gamma) and identical admissible point
representations, with independent checks of noncommuting generator handling.
For generic non-U(1) models, use the full physical sector and separately
validate terms with an independent small-Hilbert-space matrix oracle.

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

- reference.py: NumPy dense small-sector oracle, physical case metadata,
  bond construction, independent self-tests. Safe metadata-only mode for
  large cases.
- run_quspin.py: executable *U(1)-only* QuSpin/SciPy CSR pilot with a hard
  dimension safety cap; prints JSON and records actual direct residual.
  It is intentionally not a claim of the best possible QuSpin performance.
  Tests against an installed QuSpin runtime remain to be performed.
- XDiag: next implementation step, starting with Spinhalf(N,N-m) and
  matrix-free versus CSR modes (C++ recommended). XDiag release and
  source SHA must be pinned. The XDiag benchmark guide explicitly compares
  single Lanczos iteration, so record end-to-end timing independently.

Suggested local smoke (the second command needs NumPy, the third requires
QuSpin, NumPy, SciPy):

~~~bash
python3 benchmarks/cross_library/reference.py --self-test
python3 benchmarks/cross_library/reference.py --shape 6 --ndown 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 benchmarks/cross_library/run_quspin.py --shape 6 --ndown 2
~~~

On the current development host, QuSpin and XDiag are unavailable. Only
the standalone reference and Python syntax checks have been executed.

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
