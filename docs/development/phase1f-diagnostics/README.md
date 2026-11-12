# Phase 1-F evidence

Base: `b756b0009f08027e8a2ae8a1e349d78d8e81d9c5`.
Read [the implementation/validation report](../phase1f-diagnostics-reproducibility.md).

| File | Meaning |
| --- | --- |
| `command-manifest.json` | Actual before/after builds, compiler flags/libraries/source hashes, every original command/status/log hash and elapsed time |
| `summary.json` | New comparison/control counts, existing suite outcomes and report hashes, including separately labeled initial versus corrected stability runs |
| `kernel-comparisons.csv` | All 408 before/off/on comparisons; original solver statuses, exact vector hashes and complete observed field/pair data |
| `selected-observations.json` | Near-cluster, shifted, excluded/weak-ground examples showing residual certification separately from inclusion |
| `analytic-controls.json` | Phase 1-E deterministic controls, including the **constructed**, non-solver physical-accuracy counterexample |
| `*-stability-format-fixed*.json` | Actual final 82 kernel / 26 CLI results per profile, with unchanged scientific assertions |
| `format-correction.json` | Initial parse failures retained, original solver status zero, exact test-only exponent correction and corrected source hash |
| `*-historical.json` | All 27 historical comparisons per profile, input hashes, dimensions, energies and measured errors |
| `*-acceptance.csv` | All 3,072 seeded physical CLI outcomes, including 768 failed processes across the three profiles |
| `*-probe-provenance.json` | External-copy instrumentation patches and executable/archive provenance for the unchanged existing CLI probe |
| `performance.csv` | All 126 original before/off/on timings, RSS, H calls **and actual solver statuses**; the ALG=2 subset is budget-exhausted |
| `converged-alg2-timing.json` | All 63 additional successful ALG=2 timings at explicitly documented larger test workspace |
| `performance-summary.json` | Paired workload summaries, status sets, timing ranges, RSS ranges and common H-call counts |
| `provenance.json` | Protected-source/helper/stdout checks, maintained source hashes, linked-library hashes and observed machine/runtime metadata |

JSON-valued CSV cells retain their arrays/mappings. Unavailable diagnostics have
explicit flags; a numeric placeholder is not an accuracy assertion. Numerical
findings in characterization commands remain warnings/nonconvergences even when
the command itself returns zero. A passing failure-path check verifies failure;
it does not make the attempted solve successful.

The two initial nonzero stability commands remain unchanged in the raw command
manifest. Fresh, corrected kernel/CLI runs resolve their test-format failure.
The initial O2/GNU 13 CLI stability suites were not reached; `summary.json` marks
them unrun and lists the later actual 26-case results separately. No production
source was changed during those reruns; only the test driver's exponent width was
corrected. The source archive predates that test-only correction and later report/
reduction tools; hashes of the corrected test driver and final maintained code
are recorded separately. Future clean runs of the committed code include the fix.

Current-instance raw evidence is outside the checkout:

```
/workspace/phase1f-validation/manifest.json
/workspace/phase1f-validation/{before,after}-source.tar
/workspace/phase1f-validation/{o2,gnu13,bounds}/before/
/workspace/phase1f-validation/{o2,gnu13,bounds}/after/
/workspace/phase1f-validation/{o2,gnu13,bounds}/diagnostics/
/workspace/phase1f-validation/{o2,gnu13,bounds}/stability-format-fixed{,-cli}/
/workspace/phase1f-validation/{o2,gnu13,bounds}/{fields,stress,precision,independent,acceptance}/
/workspace/phase1f-validation/{o2,gnu13,bounds}/historical.json
/workspace/phase1f-success-timing/report.json
```

The full runner and supplementary timing command are documented in the report.
For the retained run's focused serialization correction, the exact follow-up was
performed for each profile with a **fresh** output path:

```sh
. /workspace/ed2-env/activate.sh
profile=o2  # repeat for gnu13 and bounds
cp /workspace/ED2/tests/lanczos_stability_driver.f90 \
  /workspace/phase1f-validation/$profile/after/tests/
make -C /workspace/phase1f-validation/$profile/after check-lanczos \
  LANCZOS_OUTPUT=/workspace/phase1f-validation/$profile/stability-format-fixed
```

Regenerate this reduction only after the raw runs/follow-ups complete, into a
new directory (the reducer refuses an existing output):

```sh
python3 -B tools/summarize_phase1f.py \
  --runs /workspace/phase1f-validation \
  --success-timing /workspace/phase1f-success-timing \
  --output /workspace/new-phase1f-evidence
```

No binaries, package files, native wavefunctions, object/module files, full
vectors or generated Autotools files are committed. Raw directories are not a
promise of remote persistence: use the committed generation code, seeds,
compiler/flags, base revision and file/library hashes for reproduction. Machine
metadata was observed after validation; it is not proof of an exclusive timing
host or a snapshot publication. No Cloud publication or cross-compiler-family
certification is claimed.
