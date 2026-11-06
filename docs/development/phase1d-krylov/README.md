# Phase 1-D evidence index

Read [the characterization report](../phase1d-krylov-stress-benchmarks.md) and
[the unimplemented policy proposals](../phase1d-policy-proposals.md).
Base: `6d3480c5fb8ee7a67f2ac4017156ee4368bdd1d2`.

- `summary.json`: every profile's numerical classification and regression counts.
- `default-stress.csv`, `o2-stress.csv`, `bounds-stress.csv`, `gnu13-stress.csv`:
  all 418 kernel runs per profile, including nonconvergence, input/matrix hashes,
  requested settings, actual residuals, eigenvalue/projector/orthogonality errors,
  H applications and last-candidate diagnostics. None are filtered to successful
  solves. `rng_key` resolves to `rng-states.json`; blank fields are inapplicable.
- `selected-stress-cases.json`: representative early acceptance, excluded-ground,
  weak-overlap and very-tight-energy-gate cases, with full provenance.
- `replay-check.json`: twelve selected seeded runs repeated exactly, excluding times.
- `acceptance.csv`: all 1,024 seeded physical CLI checks. The 256 failures remain
  failures; the other 768 pass the unchanged physical observable assertions.
  `acceptance-input-provenance.json` shows the physical fixture input hashes;
  per-case settings plus the committed preparation/probe scripts reproduce decks
  and thread RNG states. The probe does not erase the unseeded regression failure.
- `benchmarks.csv`: all 192 spin benchmark runs, wall/compute times, executable
  RSS high-water marks, separate wait4 RSS, H calls/cycles, memory estimates,
  actual residuals and available independent spectral/subspace errors.
- `benchmark-aggregates.json`: medians/ranges over the three starts for every
  version/model/settings/thread group. Ranges include random-start variation.
- `benchmark-provenance.json`, `benchmark-input-hashes.json`: compiler, flags,
  CPU/cgroup/affinity/backend/load information, executable/probe/harness hashes
  and input hashes. The seed formula and seed/thread IDs reconstruct RNG state.
- `benchmark-*-instrumentation.json`: lossless unified-diff lines of the external
  benchmark-only source changes. `''.join(record['diff_lines'])` reconstructs the
  exact raw diff identified by its SHA-256. JSON avoids altering Fortran source
  whitespace or adding patch-context whitespace errors to the repository.
- `*-historical.json`: every preserved historical comparison, tolerances, errors,
  input hashes and physical metadata. No reference output was regenerated.
- `regressions.json`, `regressions-gnu13.json` and `*-*.log`: exact commands/statuses
  and build/test evidence. GNU 13's sparse command exits 2; unaffected checks were
  explicitly continued and their statuses retained. Logs have only trailing
  whitespace and excess trailing blank lines stripped for review.
- `gnu13-regression-failures.json`, `gnu13-sparse-failure.log`, and
  `gnu13-sparse-failure-input/`: the original successful-exit/failed-observable
  regression, preserved without rerunning it away. Its original CLI RNG state is
  unavailable; input/binary/environment/log provenance is retained, not a claim
  of bitwise replay of that unseeded start.
- `validation.json`: protected source/oracle and new harness hashes, compiler/
  backend/package details and final integration/review evidence.

The characterization runner's zero exit means its protocol completed; it does
not mean every returned eigenpair meets the stricter assessment or every budget
converged. Use the CSV classification/status columns and summary counts.
`accepted_strict` assesses the stated ground-energy and space targets; actual
residuals are separate columns. Large-system benchmark classification
`valid_residual_only` explicitly does not certify global ground-state inclusion.

No executables, package archives, object/module files, full vectors/dense matrices,
example output trees, caches, credentials or temporary build files are committed.
The raw run tree was `/workspace/ed2-phase1d`. Recreate its profile-named layout
using the report's commands and run `tools/krylov_characterization_summary.py
--raw RAW --output NEW_REDUCTION` to regenerate the compact tables. The evidence
is about 1.6 MB; four full per-case stress CSVs preserve compiler/profile effects.

Development preflights corrected a Python list-norm expression and discovered
wait4 RSS contamination by the launcher; no production numerical fix was made.
Earlier exploratory/partial runs remain external and are not used for final
accuracy, regression or benchmark claims. The final source hashes identify the
submitted harness versions used for the final datasets.
