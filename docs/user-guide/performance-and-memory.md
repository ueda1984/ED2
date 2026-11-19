# Hamiltonian caching, performance, and memory

This guide describes **general features of the current development
`ueda1984/ED2` main branch**. The public Read the Docs site is built
from [QS-Cube/ED2](https://github.com/QS-Cube/ED2); changes available
only in a fork are **not available in the official executable** until
the upstream project merges and publishes them. Always check your
executable's source commit before enabling an option.

No optimization listed below changes the physical Hamiltonian.
Most controls are intentionally **off by default**; explicit opt-in
and a suitable verification run are recommended. A valid spatial
symmetry must commute with the entire Hamiltonian, including all
oriented DM couplings and fields.

## Check the basic run first

Start from a known working input deck:

```bash
git rev-parse HEAD
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
./source/QS3ED2 < input.dat > baseline.log
```

Keep the input files, executable commit, compiler, compilation
flags, BLAS/LAPACK backend, OpenMP/BLAS thread counts, and computed
energies/observables. Test an optimized mode on a **new output
directory**, not over the baseline result. A run that exits normally
is not itself proof that a requested symmetry sector is physical.

### `wk_dim` and `MNTE`

These are **Fortran NAMELIST input parameters**, not environment
variables:

| Setting | Meaning | Caution |
| --- | --- | --- |
| `wk_dim=-1` | Cache the entire projected Hamiltonian | May require considerable memory |
| `wk_dim=1` or another smaller positive integer | Cache only the specified initial rows and calculate the remainder on demand | CSR options requiring a full cache are incompatible |
| `MNTE=-1` | Automatically choose the historical conservative row capacity | Often overallocates when most transition coefficients vanish |
| `MNTE>=0` | Explicit number of non-diagonal transition slots per row | Too small can terminate with `MNTE is too small` |

The exact projected dimension depends on the occupation window and
the chosen symmetry sector. The number of sites by itself is **not**
the dimension of the cached matrix.

## Generic optional cache controls

The following environment variables are available in the **fork
development main as of 19 November 2026**, not necessarily in
the official QS-Cube release.

| Environment variable | Effect | Conditions and tradeoff |
| --- | --- | --- |
| `ED2_CACHE_COEFF_BOUND=1` | Certifies a conservative row-capacity bound from nonzero field/exchange coefficients **before** allocation | Applies when `MNTE<0`; does **not** override explicit `MNTE`. No compression if all channels are present |
| `ED2_CSR_COMPLEX_CACHE=1` | Keeps every occupied cached matrix record in ordered complex CSR form | Requires the full projected matrix cached. Original padded arrays remain allocated unless separately released |
| `ED2_CSR_RELEASE_PADDED_CACHE=1` | Releases the original padded arrays **after** successful CSR construction | Requires complex CSR. Reduces long-lived cache payload; the old arrays cannot be read by external clients afterwards. Does **not** reduce the initial construction peak |
| `ED2_CSR_DIRECT_BUILD=1` | Builds exact ordered complex CSR directly in **two physical row-generation passes**, using 64-row staging batches rather than a complete padded matrix | Requires both CSR and release flags, and a full cache. Uses less simultaneous cache-array payload but computes each row twice; may be slower during initialization |

All flags are omitted/zero by default. Invalid values or incompatible
prerequisites are rejected. **`ED2_CSR_DIRECT_BATCH_ROWS` is not
included in this main-branch interface**; batch-size tuning is a
separate experimental PR and must not be documented as an available
release feature.

For memory-constrained, fully cached generic models, an example is:

```bash
export ED2_CACHE_COEFF_BOUND=1
export ED2_CSR_COMPLEX_CACHE=1
export ED2_CSR_RELEASE_PADDED_CACHE=1
export ED2_CSR_DIRECT_BUILD=1
./source/QS3ED2 < input.dat > direct-csr.log
```

For minimum one-time setup latency, **do not assume two-pass
construction is faster**. First compare the original cache and
ordinary CSR with the same model, solver tolerances and threads.
For a partially cached run, leave all CSR options off and use the
original on-the-fly operator path.

## What the memory numbers mean

With automatic `MNTE`, the historical padded row capacity is

```text
width = 1 + 2*NO_one + 8*NO_two
```

This is a capacity bound on *possible transition records*, not the
number of physically nonzero matrix elements. The coefficient-based
bound reduces the row allocation based on exactly zero operator
channels. It retains the original coefficients and transition order.

The compressed CSR form stores a one-based row pointer array plus
ordered target indices and complex coefficients. Stored zero
coefficients are retained: the optimization **does not apply a
numerical threshold** or reorder the Hamiltonian action.

Direct two-pass CSR prints a marker like

```text
Hamiltonian action: csr-direct-build passes=2, max_staging_bytes=..., peak_cache_payload_upper_bytes=...
```

This upper bound applies only to **simultaneously allocated
Hamiltonian-cache arrays**. It is **not** the maximum process resident
set size (RSS), and excludes basis arrays, OpenMP data, solver
workspaces, BLAS buffers and allocator overhead.

For an alternating spin-1/2 / spin-1 chain with ten sites
(projected dimension 7,776), a representative local check of
XYZ/complex-DM Hamiltonians found:

| Cache construction | Cache-array payload estimate |
| --- | ---: |
| Coefficient-bounded padded plus CSR (simultaneous) | 12,731,908 bytes |
| Direct two-pass CSR with a 64-row staging batch | 3,940,228 bytes |

These are **native allocation-accounting results for a specific
model**, not a claim of a percentage reduction in process peak RSS or
general runtime improvement. Different exchange terms and projected
dimensions yield different results.

## Numerical validation policy

ED2's tested **`-O2 -fopenmp`** profiles preserve all complex
`H|\psi\rangle` output bits between ordinary cached and CSR
representations for the validated real and complex DM examples.
Under **`-Ofast`**, compiler floating-point reassociation can
produce last-bit differences between padded and CSR multiplication
even with identical mathematical terms. The correctness criterion is
then an explicitly stated numerical tolerance, **not** unconditional
bitwise reproducibility. The direct CSR and CSR-release variants
are tested separately to ensure no terms were omitted.

Relevant reproducibility tests in the development fork include:

```bash
make check
# Time permitting, on a machine suitable for the large cubic ex9:
make check-long
```

`make check-long` includes **ex9, a large three-magnon cubic
workload**. On ordinary CI machines, the project runs the lighter
historical ex1–ex8 cases independently; do not describe the complete
`make check-long` suite as a short smoke test.

For scientific work, inspect eigensolver residuals and compare
independent small-system spectra, including nonzero DM models.
For numerical proof, the fork tests check **two independent complex
input vectors and the complete output bit patterns**, plus
independently constructed dense Hamiltonians for small systems.

## Which option should I select?

| Scenario | Initial choice |
| --- | --- |
| Small ED example / first-time user | Historical defaults, `wk_dim=-1`, `MNTE=-1` |
| Fixed U(1)-conserving XXZ sector | Select only compatible symmetries; try `ED2_CACHE_COEFF_BOUND=1` |
| Large real or complex DM full-cache run with many repeated H actions | Test complex CSR against the same baseline |
| Long-lived CSR storage is too large | Consider padded release if external code does not access `wk_ele`/`wk_loc` |
| **Initial Hamiltonian-cache allocation** causes memory pressure | Consider direct two-pass CSR; account for the extra row-generation pass |
| Full cache cannot be allocated at all | Use partial `wk_dim` and the ordinary on-demand Hamiltonian action; direct CSR is ineligible |

The separate [lattice and model tips](../examples/optimization-tips.md)
cover physical symmetries and examples, without promoting
unmerged cubic-specific kernels to stable public features.
