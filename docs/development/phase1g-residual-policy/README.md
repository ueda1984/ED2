# Phase 1-G evidence index

See [the report](../phase1g-optin-residual-policy.md) for policy semantics, limitations
and commands. These files are reductions of actually executed tests, not updated
reference data. Compiler/profile names are `o2` (GNU 14), `gnu13`, `bounds` (GNU 14),
and additional `default` (GNU 14 -Ofast, quick only).

- `summary.json`: per-profile status transitions, unchanged-start checks, CLI and
  malformed controls, regression counts, final FP-aware checks and report hashes.
- `*-cases.csv`: 146 contexts times legacy and three policies (584 rows/profile).
  Inputs, seeds, budgets, actions, candidate diagnostics and independent dense
  errors are retained. Vector/assessment availability is explicit; absent errors
  on failed ALG=1 calls must not be read as zero. Arrays/objects are JSON cells.
- `selected-cases.json`: detailed shifted, close-gap, excluded/weak-ground,
  tiny-gap spin-dimer, scalar-zero and all-pairs full-space evidence.
- `analytic-controls.json`: independently constructed physical-accuracy controls.
- `*-acceptance.csv`: all 1,024 legacy seeded probe outcomes for each requested
  profile; explicit solver failures are retained (not observable agreements).
- `*-historical.json`: all 27 numerical historical comparisons per requested build.
- `timing.csv`, `timing-summary.json`: 240 measurements with actual residuals,
  reference subspace assessment, H actions, seconds and process VmHWM KiB.
- `required-commands.json`: full three-profile commands, source archive hashes,
  compiler/library metadata, log hashes and exit codes; all exit zero.
- `default-commands.json`: supplementary quick commands. The initial telemetry
  checker failure is deliberately retained. It assumed gradual underflow under
  -Ofast; production expressions matched the baseline. This profile omits
  extended/historical/acceptance suites.
- `fp_followup-commands.json`: final test-only arithmetic-mode correction,
  recompiled before/after drivers and telemetry/seed/error/CLI checks; all 12
  commands exit zero. No solver source or floating-point mode was changed.
- `provenance.json`: protected-file/helper checks, current source/harness hashes,
  linked library hashes, NumPy/Python and machine observations.

Full build trees, vector files and logs were transient external artifacts under
`/workspace/phase1g-validation`, `/workspace/phase1g-default-validation` and
`/workspace/phase1g-fp-final`. Their commands and hashes are retained here; rebuild
with new output names using the maintained runner. Initial archives predate the
final test-only underflow checker correction; production hashes are unchanged.
Values are binary64 measurements; observed scales and dense error estimates are
not certified spectral bounds. Warnings and nonconvergence are scientific results,
not successful precision guarantees. No failure result file is promoted to success.
