# Phase 1-C evidence

See [the numerical change report](../phase1c-lanczos-stability.md) for the
algorithms, stopping policy, reproduction commands and limitations.
Base commit: `2050208999940044c155cce7fb60d42b7228481d`.

Each profile was built from a fresh external source copy. Recorded hashes verify
that the tested production, test and Automake sources equal the submitted sources.
The unchanged independent tests were run before and after the fixes; neither
reference data nor their assertions were adjusted.

| Check | Default (`-Ofast`) | Conservative (`-O2`) | Bounds (`-O2 -fcheck=all`) |
| --- | ---: | ---: | ---: |
| Clean GNU Fortran build, no argument-mismatch workaround | pass | pass | pass |
| Original independent cases | 21/21 | 21/21 | 21/21 |
| Kernel stability checks | 82/82 | 82/82 | 82/82 |
| CLI stability checks | 26/26 | 26/26 | 26/26 |
| Sparse-field CLI | 70/70 | 70/70 | 70/70 |
| Hamiltonian operator checks | 10/10 | 10/10 | 10/10 |
| Field-record permutations | 28/28 | 28/28 | 28/28 |
| LAPACK interface/error paths | 11/11 | 11/11 | 11/11 |
| Standard regression | 1/1 | 1/1 | 1/1 |
| Extended examples | 9/9 | 9/9 | 9/9 |
| Historical numerical comparisons | 27/27 | 27/27 | 27/27 |

Kernel and CLI stability cases run at one and two OpenMP threads, with BLAS
single-threaded. The preserved independent and sparse suites also cover these
thread counts. Standard/extended examples use one OpenMP thread. Exact flags,
compiler/runtime hashes and error maxima are in [validation.json](validation.json).
`normal` in filenames denotes the conservative profile.

The 82 kernel checks include 16 intentional invalid-input/nonconvergence controls;
the 26 CLI checks include four intentional insufficient-iteration controls. These
pass only when the expected nonzero status/exit and absence of result files are
observed. They are separate from the 21 formerly failing independent cases, which
now all succeed without expected failures. Eleven additional focused ALG=1
floating-point-trap controls pass; no general trap-safe LAPACK claim is made.

Files retained:

- `pre-fix.json` and `pre-fix-*.log`: the three original failures, before changes;
  all-case counts were 18 passed / 3 failed. The old executable hash matches
  Phase 1-B's committed bounds-build evidence.
- `*-kernel.json`, `*-cli.json`: individual case input hashes, eigenvalues,
  independent actual residuals, status/reasons, orthogonality and applicable
  degenerate-space projector errors.
- `*-independent.json`: all 21 original cases, including eigenvalues and input
  hashes; bulky magnetization/correlation arrays omitted.
- `*-sparse.json`: preserved cases, hashes, tolerances and error summaries;
  bulky observable arrays omitted.
- `*-historical.json`: all 27 comparisons with tolerances, observed differences,
  reference hashes, output hashes and physical metadata from the existing harness.
- `*-autoreconf.log`, `*-configure.log`, `*-build.log`, `*-stability.log`,
  `*-sparse.log`, `*-lapack.log`, `*-check.log`, `*-check-long.log`,
  `*-baseline.log`: build and test evidence. Trailing whitespace/blank lines were
  stripped for repository whitespace checks; substantive messages are unchanged.
- `trap-controls.json`: eleven additional zero/scalar/1D/exact-start/near-breakdown
  ALG=1 kernel and CLI runs using `-ffpe-trap=invalid,zero,overflow`.

Only text evidence is committed. Executables, objects, module files, full example
outputs, temporary trees and toolchain packages remain outside the repository.
Recorded absolute paths identify this cloud run; they are not required by the
solver or the portable reproduction procedure. RNG starts remain unfixed, so
reruns should satisfy the numerical checks rather than reproduce individual
iteration counts or last-bit errors.
