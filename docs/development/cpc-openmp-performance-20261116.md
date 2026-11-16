# ED2: OpenMP-compatible incremental improvements (2026-11-16)

## Scope and source status

The earlier single-thread-only orbit prototype is replaced by a shared-memory OpenMP implementation. The optimized path now runs at the requested team size; it does not set the global OpenMP thread count to one. The existing Hamiltonian cache, Hamiltonian-vector parallelism and eigensolver precision/acceptance policy are retained.

This is a tested combined REVIEW PATCH on `codex/ed2-cpc-performance-audit-20261116`, not a production merge. Main, upstream and the manuscript branch remain unchanged. `tests/performance_openmp_20261116.patch` is standalone against baseline `719da7035e66ce884e9de00376e4ded5b3e5c459`; it replaces, and must NOT be applied on top of, the earlier `performance_audit_20261116.patch`.

Only two native source files change:
- state_lists.f90 SHA256: `52d2484a8e878edb4905058ab2bcc74b2164722b8f728f827a85b6ba2a24eabb`.
- krylov_solvers.f90 SHA256: `5f2b91a24c8cfa33d58b79a8204ce2cf6b889f5d02554b2dc80f3821c6503356`.

The downloadable reproduction bundle includes native sources, the exact baseline archive, the patch, reproducible build/test scripts, all paired timing records, inputs and logs, and independent CI summaries. Temporary main.f90 stage timers used for measurements are test-only and absent from the production patch.

## Race-avoidance design

A single shared visited bitmap is read and updated using OpenMP atomic operations, including when workers access different bits in the same 64-bit word. Each worker owns its occupation, translated-rank and output buffers. A worker must finish the global-minimality check BEFORE marking an orbit. A nonminimal candidate never publishes visited bits, so it cannot suppress a lower-ranked representative that another worker has not processed. A minimal orbit with vanishing stabilizer character can be marked because it contributes no state to the requested sector.

Outputs are joined and sorted by the original representative rank after the parallel barrier. Phase sums use the original generator order. Thus the implementation preserves the original basis ordering and normalizations regardless of thread scheduling. Dynamic chunking balances work without a global per-representative insertion lock. Identity symmetry directly constructs the unprojected rank list in parallel.

Permutation validity, generator periods, local spin/cap invariance and commutation are checked before the optimized path. The shared bitmap is capped at 128 MiB; total worker rank workspace is separately capped at 64 MiB. Unsupported cases retain the legacy path. These guards do NOT fall back solely because more than one thread is requested. The final sort remains serial, and no claim of ideal scaling at arbitrary core counts is made.

OpenMP specification context: atomic IOR/update and read semantics are documented at https://www.openmp.org/spec-html/5.1/openmpsu105.html . Bounds/undefined-behavior sanitizers are not race detectors; race avoidance follows the ownership/atomic protocol and is additionally exercised by repeated cross-thread tests.

## Paired total-runtime measurements

Both baseline and candidate use full, model-specific tight-row Hamiltonian cache. These are improvements BEYOND enabling the cache, not a repetition of the earlier no-cache/full-cache experiment. Same physical input, explicit starting vector, RNG seed 1984 and absolute residual target 1e-8 are used per pair. GNU Fortran 14.2, -O2 -fopenmp, system LAPACK/BLAS, BLAS=1. Each version runs one warmup and five interleaved measurements per workload/thread count. CPU masks use 1/2/4 allowed CPUs. The host is shared and has a four-CPU cgroup quota; this is engineering evidence, not an exclusive-node publication benchmark.

| Workload | OpenMP threads | Baseline median s | Candidate median s | Incremental ratio |
| --- | ---: | ---: | ---: | ---: |
| 1000-site, 2 magnons, Gamma, ALG1 | 1 | 0.949629 | 0.059183 | 16.05 |
| Same | 2 | 0.535701 | 0.044078 | 12.15 |
| Same | 4 | 0.275076 | 0.032851 | 8.37 |
| 216-site, 3 magnons, Gamma, ALG1 | 1 | 1.368576 | 0.494245 | 2.77 |
| Same | 2 | 0.742673 | 0.284480 | 2.61 |
| Same | 4 | 0.451708 | 0.195556 | 2.31 |
| 12-site non-U1 mixed-spin DM, ALG1 | 4 | 0.393266 | 0.311710 | 1.26 |
| Same DM, ALG2 keep=2/cap=32 | 4 | 6.207161 | 4.349294 | 1.43 |

All 96 completed runs (80 measurements plus 16 warmups) satisfy the residual bound; maximum observed residual is 8.7913e-9. All 48 baseline/candidate pairs have identical energies, H-application counts, physical input and initial-vector hashes in this local experiment. One externally interrupted measurement attempt was retained and rerun; it is not classified as a numerical solver failure. The summary CSV includes medians, MAD, RSS and H calls, not just speed ratios.

For the four-thread restarted solver, median peak process RSS falls from 113824 KiB (111.16 MiB) to 91848 KiB (89.70 MiB), a 19.31% reduction, while total time improves by 1.43x. This comes from constructing/storing only the Ritz vectors actually consumed and postponing the expensive subspace-closure test until it can affect acceptance. No precision reduction, relaxed residual threshold or fast-math is used. Small RSS changes in ALG1 are not counted as memory improvements.

The timed source differs from the final source only by clarifying comments. Both exact source snapshots and binary hashes are retained in source-comment-provenance.json; the binaries are not claimed byte-identical because embedded diagnostic line numbers can differ. The final source bytes above are independently verified in CI.

## Correctness and multithread validation

Local GNU14 checked runs tested the same 379 symmetry/occupation fixtures at 1,2,3,4,8 threads. Every candidate basis was byte-identical across these thread counts. Eight threads oversubscribe this host and are used for correctness only. Two existing zero-cap rank/unrank defects are deliberately corrected and checked against independent product-state enumeration; ordinary cases preserve legacy ranks and normalization.

Independent GNU13/OpenBLAS CI completed successfully:
https://github.com/ueda1984/ED2/actions/runs/43793344388

The CI applies the standalone patch to the immutable baseline and checks both source hashes before compilation. It verifies:
- 379 fixtures at 1/2/4 threads plus five larger/shared-bit-word fixtures repeated four times: 1197 candidate basis executions, deterministic rank and normalization output.
- Native full-cache 216-site/three-magnon and 1000-site/two-magnon eigenpairs at 2/4 threads, with the multithread optimized path explicitly observed in logs.
- 18 independent non-U1 mixed-spin energy/Cartesian-observable checks at four threads.
- 209 adversarial Krylov definitions at 1/2/4 threads: 627 paired comparisons with identical classifications and independent acceptance results. There are 251 strict accepts, 200 accuracy/inclusion warnings and 176 nonconvergences in EACH implementation, not 627 converged cases. With matched checked-build flags, iteration histories also agree in all 627 comparisons.

An earlier CI attempt compared -O2 baseline to -O1 checked candidate and incorrectly demanded identical iteration counts. One clustered-eigenvalue case converged accurately in both versions but used 54 versus 53 restarts, causing that test gate to fail. The revised gate retains all trajectory differences, compares independently assessed numerical accuracy, and uses matching compiler flags. No solver or independent accuracy tolerance was loosened. The failed first run is retained at https://github.com/ueda1984/ED2/actions/runs/43792608076 .

## Native 1000-site / three-magnon basis only

The checked native candidate constructs raw dimension 166167000 and Gamma dimension 166167 with exactly identical basis bytes at 1/2/4 threads. Observed basis times were 7.441, 3.879 and 2.150 seconds respectively, one run each on the shared host. Peak RSS remained about 26.4 MiB. All norms match sqrt(1000), ranks are strictly ordered, and 31 sampled configurations were independently checked against every translation. This validates the native BASIS path only; no converged 1000-site/three-magnon ED2 eigenenergy is claimed.

## Adoption boundary

This candidate meets the intended single-node multithread design: the new basis builder uses the requested OpenMP team, and the restarted-solver work/memory reduction does not remove existing Hamiltonian or vector parallelism. Before freezing CPC performance claims, test on the actual exclusive PRESTO node with physical-core affinity, controlled BLAS threading and repeated allocations. Performance at 16/32 or more cores, NUMA cross-socket effects and other Fortran compilers remain unmeasured.

No global OMP_NUM_THREADS override, thread-count reduction inside production code, native algorithmic approximation, source release/tag or publication-ready speedup claim is introduced.
