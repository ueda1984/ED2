# Phase 1 LAPACK validation evidence

Read [the change report](../phase1-lapack-interfaces.md) for scope, reproduction,
interpretation and limitations. This evidence concerns the feature based on
`bfe96f0754c662d40d0207e7dc3d1ec6863486dc`; the exact tested source bytes are
identified in `validation.json` by file SHA-256, independently of the later
commit containing this documentation.

| Profile | Build, no mismatch workaround | Interface checks | Standard regression | Extended regression | Historical arrays | Independent cases |
| --- | --- | --- | --- | --- | --- | --- |
| GNU defaults (`-Ofast`) | pass | 11 pass | 1 pass | not run | not run | 18 pass, 3 known failures |
| Conservative (`-O2`) | pass | 11 pass | 1 pass | 9 pass | 27 pass | 18 pass, 3 known failures |
| Bounds (`-O2 -fcheck=all`) | pass | 11 pass | 1 pass | 9 pass | 27 pass | 18 pass, 3 known failures |

- `validation.json` records compiler/backend/platform, flags, executable/source
  hashes, counts, per-case independent outcomes, all 54 historical comparisons
  across the two accuracy profiles, current-run input hashes, and strict
  implicit-interface compilation results.
- `*-build.log` records actual compiler/link invocations with no
  `-fallow-argument-mismatch`. The bounds log retains unrelated warnings; this is
  not a warning-free-build claim.
- `*-lapack.log` records the explicit interface/error-path target and all 11
  checks passing, including the default fast-math profile.
- `*-check.log` and `*-check-long.log` record the standard and extended regression
  commands' successful completion. The independent numerical suite separately
  returns exit 1 for each profile; it is not hidden by these successful tests.
- `early-check-lanczos.log`, `zero-alg1.log`, `zero-alg2.log` retain representative
  conservative-build stdout for the three known failures. Their deterministic
  input definitions are unchanged in `tests/numerical_baseline.py`; corresponding
  baseline input decks remain in `../baseline/failures/`.
- A comparator negative control, also recorded in `validation.json`, replaced
  one energy with NaN **in a temporary copy outside the checkout**. It returned
  exit 1 with 26 passing comparisons and the one intended failure. No historical
  reference file was changed.

Captured `.log` files have trailing whitespace and excess blank lines at EOF
trimmed for review. Numerical JSON/input hashes are not normalized. Raw results,
full independent reports, all generated inputs and outputs remain under
`/workspace/ed2-phase1-final/` on the validation machine; they are reproducible with
the documented commands and are not committed as build artifacts. Timings from
concurrent cloud processes are not performance benchmarks. Superseded preliminary
runs are excluded from this final matrix.
