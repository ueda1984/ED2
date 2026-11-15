# CPC performance: reproducible pinned-node experiment protocol

**Status:** script and a GitHub-hosted *shared-runner smoke test*.
A controlled dedicated compute-node benchmark has **not yet been
executed** as part of this work. The present wall-time numbers
in the CPC manuscript are exploratory CI observations.

## Why this is needed

The 8–12-site memory-time studies on the alternating
spin-1/2/spin-1 XYZ/Gamma/staggered-DM chain passed independent
numerical and physical checks. They demonstrate strong
configuration dependence, but wall time on GitHub-hosted
shared runners is sensitive to CPU placement and other
uncontrolled loads. It is premature to publish a universal
speedup or parallel-efficiency result from those timings.

The [paired-cache study](cpc-cache-pareto-evaluation.md)
shows that *full caching at the independently proven
row capacity* can use less nominal cache memory
than default-width half caching on the specific model.
This motivates a repeatable **dedicated-node** comparison
using the very same model, starting vector, eigenpair
target, residual request, observables and solver workspaces.

## A. Benchmark script and provenance

The runner is
`tests/cpc_controlled_node_benchmark.py`.
It delegates every numerical case to the existing
`tests/cpc_stationary_cache_pareto.py`, which
has independent SciPy reference checks and does not
modify ED2's numerical kernels.

The wrapper:
- sets and verifies Linux CPU affinity through
  `os.sched_setaffinity`, inherited by child ED2 processes;
- pins OpenMP via `OMP_DYNAMIC=FALSE`,
  `OMP_PROC_BIND=TRUE`, `OMP_PLACES=threads`,
  and passes 1, 2 (or more) threads explicitly;
- forces OpenBLAS/MKL/Goto/NumExpr thread count to one
  to avoid unintended oversubscription;
- records CPU model and topology, allowed and
  selected CPU affinities, physical memory,
  compiler, linked library information, platform,
  selected scheduler job metadata and the full
  execution environment;
- hashes the ED2 binary, benchmark and wrapper
  scripts, optional exact source archive and
  **every generated input deck** using SHA256;
- performs requested warmups separately from
  measured repetitions, retains all raw solver
  outputs and numerical diagnostic sidecars,
  and writes a per-OMP `report.json` plus
  `environment.json`, `input-manifest.json`
  and top-level `summary.json`.

The machine-report object is evidence of *what
was observed*. The source SHA identifies source
provenance; it does not prove the binary was
built from that source without a clean build or
archived build record.

### Operator attestation, not inferred exclusivity

The script **cannot automatically prove**
that no other job shares the CPU cores
or that the scheduler allocated an exclusive
node. A scheduler job ID alone is insufficient.

For measurement intended for the manuscript,
an operator must first verify an exclusive
allocation, then explicitly choose the
`--exclusive-allocation-verified` flag.
It is a human attestation recorded in JSON,
not an automatic certification. The script
also requires a full source SHA or exact
source archive in this mode.

For exploratory checking on CI or a workstation,
use `--allow-shared` instead. Results are
explicitly labelled
`shared_or_exclusivity_unverified` and
**must not** be promoted as controlled-node
speedups.

## B. Procedure on an already allocated dedicated node

**1. Verify the machine and node reservation.**
Confirm the batch scheduler has granted
exclusive use of the intended CPU cores, check
CPU topology and allowed affinity, and ensure
other system load is absent or documented:

    hostname
    lscpu -e=CPU,CORE,SOCKET,NODE
    taskset -pc $$
    free -h

Choose a CPU-list entirely within the allowed
affinity. Avoid comparing different NUMA nodes
or hyperthread siblings unintentionally.
Do not infer or hardcode CPU indices from a
different machine.

**2. Build a fixed source revision.**
Prefer a reproducibly extracted source tarball
built on that same machine with recorded
Fortran flags and BLAS/LAPACK linkage.
Keep the exact tarball and full source SHA;
record toolchain versions.

**3. Run at least seven measured repetitions**
with two separate warmups (adjust based on
available resources). An example for an
operator-verified exclusive allocation is:

    python3 tests/cpc_controlled_node_benchmark.py \
      --exe /absolute/path/to/QS3ED2 \
      --output /new/nonexistent/cpc-node-results \
      --cpu-list 0-7 \
      --threads 1,2,4,8 \
      --sizes 8,10,12 --algs 1,2 \
      --D 0.12 --keep 2 --capacity 16 \
      --residual-atol 1e-8 \
      --warmups 2 --repeats 7 \
      --source-revision <full-40-character-git-sha> \
      --source-archive /absolute/path/to/source.tar.gz \
      --exclusive-allocation-verified

**Replace** `0-7` and the source identifiers
with the ones actually allocated. This is a
template, not evidence of a run. Never pass
the exclusivity flag on shared GitHub CI.

**4. Check acceptance before using the timings.**
Each `paired-ompN/report.json` must show
`pass_complete=true` with all cases accepted
and no numerical or physical reference violations.
Audit individual residuals, low-state energy
inclusion, H application counts, local
moments and Cartesian two-spin correlations.
If a solver does not converge, retain the
failure and modify the experimental design
openly instead of silently dropping the case.

**5. Report environment and statistics.**
Use the raw measurements rather than only
summary medians. Report median and MAD
or an appropriate uncertainty assessment,
CPU model/cache hierarchy, compiler and
flags, BLAS/LAPACK, affinity, working
precision, cache allocation, and source
plus binary checksums.

### Strict post-run numerical and resource audit

After the launcher finishes, use the independent
**standard-library-only** postprocessor:

```bash
python3 tests/cpc_controlled_node_analyze.py \
  --root /path/to/cpc-node-results \
  --output /new/nonexistent/cpc-node-audit \
  --require-publication
```

The `--require-publication` flag **refuses** shared/exclusivity-unverified
runs and fewer than five measured repetitions. The checker rehashes
every original input file, validates the source/binary provenance,
checks that each requested OMP, size, algorithm, cache policy,
measurement repeat and warmup exists exactly once, and re-evaluates
numerical acceptance (direct residual, independent energy,
all local magnetic moments and all Cartesian nearest-neighbor
correlations). It recomputes cache allocation formulas and confirms
the same explicit starting vector, seed, requested numerical target
and cache policy along each matched comparison.

Successful output is `audit.json` and
`paired-policy-stats.csv` with trial counts,
per-policy **wall-time medians and median absolute deviations**,
peak process RSS medians/MAD, Hamiltonian-application counts,
nominal cache array sizes and matched half-generic/full-tight
ratios. Do not quote only the ratios: retain all
observations, source/inputs and the full CPU/environment log.

For **shared CI testing**, omit `--require-publication`.
The report is labelled
`exploratory_shared_or_insufficient_repeats`;
the gate for publication-quality results remains **false**.
The [smoke workflow](https://github.com/ueda1984/ED2/actions/workflows/cpc-controlled-node-protocol.yml)
explicitly checks that its shared-runner input
**is rejected** if `--require-publication` is requested,
while the exploratory audit still succeeds.

**Human review remains mandatory.** The recorded operator attestation
is not machine-verifiable proof of exclusive scheduling;
the postprocessor also cannot independently establish
stable CPU frequency, NUMA placement, thermal/load conditions,
or statistical representativeness of a single node allocation.
Its `publication_gate_passed` Boolean represents a
**minimum data-completeness screening**, not editorial approval
or a hardware-independent performance guarantee.

### Important limitations

The current launcher runs each thread count
sequentially, with configuration policies
interleaved over repeats but in a fixed order.
A journal-ready strong-scaling claim would
also benefit from repeated independent job
allocations or randomized ordering, thread
binding checks, and per-core power/frequency
observations where available. Even when
CPU affinity is fixed, hardware turbo
modes and NUMA placement can affect
timing variability.

The program's current one-ground-state
study does **not** establish behavior for
all low-energy spectral structures. The
two-eigenpair restart-workspace and
[near-degenerate projector](cpc-stationary-sparse-doublet.md)
tests must be reported separately. The
decoupled-edge model is deliberately
adversarial and must not be presented as
the connected stationary DM physical model.

## C. Shared-runner functional CI

[Completed shared-runner protocol smoke #43705348634](https://github.com/ueda1984/ED2/actions/runs/43705348634)
passed, including both OMP=1 and OMP=2, the
source SHA, source-archive SHA256, generated-input
checksum manifest and independent numerical acceptance.
The [CI workflow](https://github.com/ueda1984/ED2/actions/workflows/cpc-controlled-node-protocol.yml)
builds a fresh, separately extracted source
distribution, sets a small allowed CPU
affinity, and exercises OMP=1/2 on the
8-site model in explicitly **shared**
mode. It checks the presence of a full
source revision, source-archive digest,
input-file manifest and complete raw
per-thread reports. Success in this
CI run verifies **the protocol mechanics**,
not isolated-node or CPU-scaling performance.

## D. Relevance to the CPC paper

With actual controlled-node measurements,
the manuscript can make reproducible
limited-scope statements about the
memory/runtime frontier and OpenMP
efficiency at the *specified* physical
workload. Until then, keep the analytic
cache-size formulas and independent
sparsity-row bound as strongest evidence,
and present shared-runner elapsed times
as exploratory only.
