# Opt-in cubic local-bond Hamiltonian cache phase profiling (PR #15)

Status: Draft stacked on PR #14, disabled by default. This is a strictly observational experiment, **not** another Hamiltonian optimization. Neither the older PRs nor main are altered.

## Complete physics and performance context

The pinned N1000 PR13/PR14 PRESTO jobs have independent A1g and A2g/A1u/A2u trial audits with **12/12 matching source-input SHA256, dimensions and eigenvalues within 2.274e-13**, and both native direct residuals below the requested 1e-8 threshold. One of the twelve, A1u at 4 OpenMP threads, converges at Lanczos k=105 (PR14) instead of k=110 (PR13), yielding H applications 210 vs 220, but both runs reach the same Ritz eigenvalue within 1.137e-13 and pass direct residual (9.71160454e-9 vs 1.60351564e-9). Both have one reconstructed vector and seed=1984. This is consistent with differing convergence gates after a small floating-point diagonal summation-order change, **not proof of the cause**. Never require bitwise identical H counts in a physical equivalence check.

The PR14 local sparse cache scans <=54 incident bonds per 3-magnon representative instead of 9*N=9000 bonds per state, but the PRESTO 1000-site A1g cache phase only improved from 1.5553s to 1.2503s (T1), 0.4369s to 0.3535s (T4), 0.2328s to 0.1910s (T8), corresponding to ~1.24x/1.24x/1.22x. It is therefore important to attribute the remaining cache cost before making another change.

## What gets measured

Set `ED2_CUBIC_CACHE_PROFILE=1` **in addition to** all four PR14 fast-path opt-ins, including `ED2_CUBIC_LOCAL_BOND_CACHE=1`. The native builder rejects enabling the profiler without the local bond cache. It first times the existing anchored representative-key setup (wall clock, serial), then native incidence-list setup (wall clock, serial), then the complete OpenMP column-processing region (wall clock, parallel).

Within the OpenMP column loop the profile counts columns, unique incident input bond IDs, active exchange bonds, resolver invocations, sampled resolver calls and stored nonzero hops. It logs aggregate per-thread elapsed-time sums separately for (a) projected state initialization and spin decoding, (b) incidence gathering and bond-ID sorting, (c) active-bond mathematical operations including lookup and storage. The latter already contains the resolver cost; it must not be added to the sampled resolver contribution.

To avoid instrumenting every hot-path lookup, only pairs satisfying `modulo(column_ID+input_bond_ID,31)==0` are timed around `resolve_cached_hop` (1-in-31 deterministic stratified bond-ID samples). The actual total invocation and sampled invocation counts are logged. The sampled resolver aggregate is extrapolated by total_calls/sampled_calls as a **rough diagnostic**, not an unbiased estimate of true hot-path wall time.

All per-thread phase sums are *sums of measured thread elapsed seconds*, and the parallel wall time is a distinct measurement. Do NOT add the two units together or report the extrapolated resolver time as a direct 1/4/8-thread speedup. Probe overhead and sampling bias make profiled total times inappropriate for publication benchmarking. Always compare against an unprofiled control from the exact pinned source.

The same native basis, local-bond cache and offdiagonal/diagonal code execute with the profiler disabled (default) or enabled. The profiler never modifies the cached matrix, lookup table, Lanczos thresholds or random number seed.

## Independent proof that sampling cannot perturb cache bits

The CI test `tests/cpc_cubic_cache_phase_profile_regression.py` runs the production `tests/cubic_cache_driver` twice from byte-identical generated inputs: once with profiling disabled, once enabled. The `ED2_CUBIC_CACHE_PROOF_DUMP=1` test-only output contains every nonzero and diagonal complex cache value by its exact IEEE-754 bits. The proof requires byte-for-byte matching complete native cache dumps, checks the native flags/printed timer/counter fields, and enforces a negative unsupported-profile-without-local-bond-cache test. Runs for cubic sides 4 and 10 at OpenMP 1/4; gfortran13 and 14 matrix required.

All previous CI source-provenance, exact offdiagonal PR14 reference and physical Lanczos checks must remain successful. The independent PRESTO benchmark/audit JSON records `cache_profile` separately and requires native profile markers only for the fast opt-in path, not the side-4 forced-legacy baseline.

## PRESTO protocol

Use a new staged immutable commit and a new PBS job directory. Pilot N1000 Gamma/O_h A1g (1/4/8 threads) with profile explicitly enabled; collect output lines starting `Cubic cache profile:` from each native trial. Compare raw input SHA256 and native physics to the previous PR14 unprofiled A1g stage; do NOT compare profiled execution times as speedup numbers. If needed, run an unprofiled PR15 control separately to show instrumentation overhead. Continue all four sectors after the A1g pilot.

The `ED2_CUBIC_CACHE_PROFILE=1` condition must be removed from any normal/generic ED2 `make check` tests and from publication-grade timing jobs. PRESTO exclusive-node operator review and appropriately repeated unprofiled measurements remain pending.