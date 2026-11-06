# Phase 1-D proposals: convergence and inclusion contracts

These are proposals only. Phase 1-D changes no production solver or policy.
Read [the stress report](phase1d-krylov-stress-benchmarks.md) and its failures first.

Before library extraction, agree on separate caller-visible controls/results:

1. Preserve `lnc_ene_conv` as the legacy energy-change control for compatibility.
   Add independently specified absolute/relative residual tolerances and record
   the scale used. A `sqrt(energy_tolerance)` default is not an observable or
   individual-eigenvector accuracy guarantee, especially across weak gaps.
   Changing the default requires explicit review and historical comparisons.
2. Define binary64 behavior for requested energy changes below roundoff. Options
   include a diagnostic/rejection of unattainable requests or an explicitly
   reported roundoff-limited stop after a residual/energy-stagnation assessment.
   Do not silently loosen the user's tolerance or convert budget exhaustion into
   success. Test zero/offset/scaled Hamiltonians and large relative cancellations.
3. Make statuses distinguish eigenpair validity, convergence by the requested
   policy, and evidence of lowest-state inclusion. Residual zero in an excluded
   starting invariant subspace cannot establish global optimality. Expose initial
   vectors/blocks and owned RNG state, with deterministic multi-thread tests.
   Independent starts and block-rank diagnostics can increase confidence, not
   provide an unconditional global proof.
4. Treat an isolated degenerate/near-degenerate cluster as an invariant subspace.
   Report Ritz values, actual residuals, orthogonality and a posteriori spectral
   separation estimates where available. Do not promise individual vectors when
   the residual allowance is larger than the resolving gap. Wanted cuts through
   clusters need a documented contract and potentially an enlarged requested block.
5. For physical observables, require an accuracy assessment consistent with the
   gap and operator norm. The observed successful energy/variance diagnostics
   with 1.5e-8 magnetization error show why scalar energy agreement is insufficient.
   Keep the unchanged 1e-9 physical assertions; use saved counterexamples to test
   a separately proposed residual policy, rather than relaxing reference checks.
6. Evaluate additional orthogonality/reorthogonalization strategies for ALG=1
   only as a separately reviewed method/performance tradeoff. Do not combine that
   change with ownership/API extraction. Profile reconstruction retry work and
   cluster failure cases before choosing a policy.
7. Reduce ALG=2 duplicate wrapper/workspace storage only in an engineering change
   with allocation-lifetime and result-ownership tests. Investigate work from
   retained H images and restart closure diagnostics with independent comparisons;
   these are not authorized optimizations in this phase.

Suggested gate for a subsequent numerical phase: explicitly approved policy;
seeded physical observable regression cases; unresolved-cluster and excluded-start
cases with honest statuses; zero/shifted/scaled operators; residual/energy targets
above and below attainable precision; full parameter/budget sweeps; at least one
independent compiler family/backend where available; and accuracy-matched timing
comparisons. Retain all nonconvergence and random-start findings in the report.
