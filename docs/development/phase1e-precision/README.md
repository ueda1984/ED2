# Phase 1-E evidence index

Source base: `b56fe21678e62b8b6216fe0ec742d5dbd62a8df6`.
See [the validation report](../phase1e-validation.md) and
[proposed precision contract](../phase1e-precision-contract.md).

| File | Evidence |
| --- | --- |
| `summary.json` | Four profiles, 136 outcomes each; criterion comparisons and exact harness/driver hashes |
| `*-cases.csv` | Every kernel outcome, input hash, RNG states, status/reason, H calls and independent endpoint metrics |
| `analytic-controls.json` | Constructed states illustrating distinctions; these are not solver outputs |
| `selected-cases.json` | Full representative near-cluster, nonconvergence, shifted, excluded and weak-start records |
| `*-acceptance.csv` | All 1,024 seeded physical CLI probes per compiler, including 256 failures per profile |
| `acceptance-summary.json` | Probe counts, executable/input provenance, failed-budget diagnosis and a representative log tail |
| `probe-instrumentation.json` | Exact patches to external copies only, probe build logs and binary provenance |
| `regressions.json` | Counts/hashes from actually run existing suites and LAPACK/standard logs; unrun suites explicit |
| `provenance.json` | Compiler/flags, linked libraries, production-source hashes verified against the starting Git revision, harness hashes and fixed-state replay |

CSV cells containing arrays or nested mappings are JSON-encoded; booleans outside
those cells are written as `True`/`False`. An empty metric denotes unavailable
evidence (not zero). In particular, failed ALG=1 vectors are unavailable. Reference
spectral scales, gaps and subspace bounds use independently known eigensystems;
they are not quantities the production solver has learned or certified online.

Current raw layout (external to Git):

```
/workspace/ed2-phase1e-{o2,gnu13,bounds,default}/  # builds
/workspace/phase1e-{o2,gnu13,bounds,default}-results/  # final contract runs
/workspace/phase1e-o2-replay/                   # identical fixed-state repeat
/workspace/phase1e-target-final/                # final optional Automake target
/workspace/phase1e-{o2,gnu13}-probe/             # instrumented archived copies
/workspace/phase1e-{o2,gnu13}-acceptance/        # seeded CLI probes
/workspace/phase1e-{o2,gnu13,bounds}-stability{,-cli}/
/workspace/phase1e-{o2,gnu13}-{independent,fields}/
/workspace/phase1e-{o2,gnu13,bounds}-{lapack,standard}.log
/workspace/phase1e-default-standard.log
```

Fresh-output safeguards intentionally reject existing result paths. Earlier
development runs (`*-final`, `*-pre-reference-correction`) are not the evidence
selected here. The final harness explicitly includes reference-eigensystem
formation defects in its numerical bound assessment. No failed solver outcome
was deleted or reclassified as success by the evidence reduction.

After reproducing the documented commands with that raw layout, regenerate this
reduction into a **new** directory:

```sh
. /workspace/ed2-env/activate.sh
cd /workspace/ED2
python3 -B tools/summarize_precision_contract.py \
  --runs /workspace --output /workspace/new-phase1e-evidence
```

The probe directories additionally contain `provenance.json` with revision,
compiler, flags and executable/probe hashes (see the probe reproduction example
in the report). The reducer checks unchanged production source in all four builds
and equality of both fixed-state reports before writing provenance. It does not
run tests or turn characterization exit zero into a claim that all solves passed.

No binaries, vectors, object/module files, native wavefunction files or dependency
packages are committed. Raw directories are current-instance evidence, not a
promise of remote persistence; formulas, seeds, committed inputs from prior phases,
probe code and hashes support reconstruction. Numerical runtime validation is
distinct from Cloud snapshot publication; no new publication was performed here.
