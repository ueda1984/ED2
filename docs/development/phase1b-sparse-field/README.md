# Phase 1-B validation evidence

The [change report](../phase1b-sparse-field-indexing.md) explains the fix and
reproduction commands. `validation.json` records the exact base, tested source
hashes, compiler/backend/library hashes, flags, counts and maximum errors. The
three profile JSON files contain individual solver, operator and record-order
comparisons. Full generated outputs remain outside the checkout.

| Profile | Sparse-field CLI | Operator cases | Record permutations | LAPACK checks | Standard | Extended | Historical | Unchanged baseline |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GNU default `-Ofast` | 70 pass | 10 pass | 28 pass | 11 pass | 1 pass | not run | not run | 18 pass / 3 fail |
| Normal `-O2` | 70 pass | 10 pass | 28 pass | 11 pass | 1 pass | 9 pass | 27 pass | 18 pass / 3 fail |
| Bounds `-O2 -fcheck=all` | 70 pass | 10 pass | 28 pass | 11 pass | 1 pass | 9 pass | 27 pass | 18 pass / 3 fail |

`pre-fix.json` records the old implementation's failure counts, three affected
paths' bounds diagnostics and the silently wrong but Hermitian full-field
permutation example. `pre-fix-{cached,on-demand,dense}.log` preserves the path
failures, and `pre-fix-input/` contains the single-site4 fixture's exact input
files. The old executable hash matches the Phase 1-A bounds validation record.
The pre-fix suite returned exit 1 (14 CLI passes/56 failures, two operator
passes/eight failures, and 28 failed order comparisons).

Each operator case checks the full dense matrix and application to the complete
basis under all three cache policies; the count 10 is fixture/thread combinations,
not a claim of only ten matrix elements or only one cache setting. The three
profile JSON files retain every operator error and each basis permutation.
Solver records retain spectrum/observable errors and input-deck hashes; the
harness/reference constructor hashes fix the reproducible fixture generation.

`*-sparse.log`, `*-lapack.log`, `*-check.log`, and `*-check-long.log` show target
completion; `*-build.log` shows compiler commands without the old mismatch
workaround. The unchanged baseline suite separately returns exit 1 for all
profiles, with the same three diagnosed failures listed in `validation.json`.
Those failures are not converted into passing or expected-failure results.

Automatic temporary-output rerun and refusal to overwrite an existing report
were also exercised; their statuses are in `target_repeatability`. Normal and
bounds historical comparisons check all 27 complete arrays and retain observed
errors, eigenvalues, dimensions and generated-input hashes in the validation
record. Historical files themselves are unchanged.

Captured `.log` files have trailing whitespace and excess blank lines at EOF
trimmed. Input and JSON hashes/values are not normalized. The full raw results
and source/build snapshots remain under `/workspace/ed2-phase1b-pre-fix/` and
`/workspace/ed2-phase1b-final/` on the validation machine; no binaries, objects,
module files or temporary output directories are committed. Concurrent cloud
timings are not performance benchmarks.
