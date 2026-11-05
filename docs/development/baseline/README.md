# Recorded baseline evidence

Captured `.log` files have trailing whitespace and excess blank lines at EOF
trimmed; numerical records and input decks are unchanged. Raw logs remain in the
external run directories identified in the workflow records.

These files preserve observations from source revision
`ca812ca0ca402bed2525556bd797cec5e0a002e8` on 2026-11-05. They are not regenerated
reference answers for the historical example suite. Read [the audit](../baseline-audit.md)
for interpretation and [the plan](../cpc-library-plan.md) for subsequent changes.

- `full-workflow.json`: archived-source build, standard test, all nine long tests,
  27 historical numerical comparisons and an exploratory debug profile with
  floating-point traps. The trapping profile's 6 passes/15 failures are diagnostic;
  ZHEEVR also traps on the standalone correctly typed probe below. This file's
  harness hash precedes the removal of blanket traps from the default debug profile.
- `checked-workflow.json`: final reproducer run after correcting that diagnostic
  profile; includes default build failure, compatibility build, actual one-test
  counts and bounds-only debug build. This run did not repeat the long examples;
  the full workflow above records that unchanged path separately. Both runners
  returned exit 1 because of the numerical failures.
- `numerical-release.json`, `numerical-checked.json`: all 21 cases of the final
  numerical harness, each with expected/actual energies and dimensions, checked
  observables where applicable, input hashes and precise failures. Each has
  18 passes/3 failures. `numerical-*.log` provides the case count/status listing.
- `reference-results.json`: all nine historical energies/dimensions, input and
  output hashes, local-magnetization summaries and XYZ correlation values. Full
  generated inputs and outputs remain in the external staging directory recorded
  by `full-workflow.json`; the reproducer regenerates them. It compares the full
  arrays, not just the summaries stored here.
- `build-default.log`, `configure-baseline.log`, `check.log`, `check-long.log` and
  `checked-build.log`: relevant build/test logs; version/CPU/link logs record the
  actual environment. Paths describe this machine, not portable installation paths.
- `exploratory-traps.json`: full 21-case trapping-build outcomes; not the accepted
  bounds-only debug result or an additional set of 15 confirmed ED2 bugs.
- `failures/`: input lists, namelist and stdout for the three failing normal-build
  cases. From a NEW copy of one case directory, run the built executable with
  `QS3ED2 < input.dat > new-stdout.log`; delete prior generated numerical outputs
  in that new copy first since the executable appends files. The output energy in
  the early-check case depends on the executable's uncontrolled random start.
- `sparse-field/`: separate bounds-checked reproducer (not part of the 21-case
  totals); the input has one field record on site 2. Its `sparse.log` records the
  out-of-bounds access. Reproduce in a fresh copy using the checked executable.
- `distribution-test-suite.log`: actual fresh source-distribution test failure;
  its sole test exits 127 because `test_common_strict.sh` is absent. The archive's
  configure and compatibility build succeeded before this check.
- `lapack-probe.f90`, `lapack-probe-traps.log`: standalone diagonal-matrix ZHEEVR
  probe. With the activated toolchain, compile with `gfortran lapack-probe.f90
  -lopenblas -o /tmp/lapack-probe` and run it; it returns INFO=0 and endpoints 1,16.
  Adding `-ffpe-trap=invalid,zero,overflow` reproduces the backend trapping behavior
  independently of ED2. It does not establish that every trapping error has the
  same cause; the Lanczos vacuum division was separately localized to ED2.
- `vacuum-traps.log`: separate floating-point-trapping evidence for the
  one-dimensional Lanczos intermediate division at `lanczos.f90:88`.
- `stale-output-check.txt`: a reused numerical-output directory was rejected and
  its previous report checksum remained unchanged.

The records include failures, not just successful results. The current runner
uses one OpenMP thread for historical examples and one/two threads for independent
small fixtures. These elapsed times are not controlled performance benchmarks.
Dependencies, random-start policy, unrun validations, and tolerance rationale are
in the audit. No credentials, reference-data edits, or solver patches are included.
