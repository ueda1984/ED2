# PRESTO 24608: cubic 3-magnon size/sector scaling (aggregate evidence)

**Status:** audit PASS, 360 distinct trials, 32 fast configurations; independent operator review **PENDING**.
These numbers are a **user-transcribed summary table**, not the original raw
`trials.jsonl` and `audit.json`. For original source and compiled binary
identities, PBS allocation/affinity details and per-trial uncertainties, consult
`$HOME/ED2_CPC_AUTO_MNTE_RESULTS/24608.jpmjpr1911/` on PRESTO; the job
archive and its hash must be preserved before CPC submission.

## Exact experiment definition

- Source: experimental cubic auto-MNTE ED2 code supplied in the previous
  `ED2_CPC_AUTO_MNTE_PRESTO_20261119.zip`. Source commit must be verified
  from `measurements/summary.json` and `cubic_auto_mnte/source-commit.txt`;
  this note is newer than the transferred experiment binary.
- Physics: periodic simple cubic spin-1/2, 3 magnons; 9N J1/J2 XXZ bonds;
  full translations at Gamma and four one-dimensional O_h irreps
  A1g, A2g, A1u, A2u; dimension matched to independent Burnside oracle.
- N=64, 216, 512, 1000; OpenMP T=1, 16.
- `ED2_MNTE=-1` autoselects 64 transitions on the local cubic opt-in path.
  Coordinate-cache, axis-pruned canonicalization, direct-hop and unused-hash
  elision enabled. Phase time excludes PBS build/setup.
- 32 fast conditions times 9 trials = 288 fast; extra 8 N64
  fast/forced-legacy paired cases times 9 = 72 legacy; 360 in total.
  2 warmups + 7 measured repeats for each distinct configuration.
- `AUDIT PASS: 360 distinct trials, 32 configurations; operator review eligible=False`.

**Frozen transcribed 32-row input:**
[`data/presto-24608-cubic-size-sector-aggregates.csv`](data/presto-24608-cubic-size-sector-aggregates.csv).

For every row the source table provides the median native phase timing,
median maximum RSS and the maximum (over measured repeats) direct
Ritz residual. The per-trial median absolute deviations (MAD), energies,
Hamiltonian application counts and independent node allocation evidence
were **not** in the copied table; do not synthesize them.

## Outcome: size, symmetry and OpenMP

| N | 1-thread A1g total (s) | 16-thread A1g total (s) | A1g speedup | 16-thread A1g RSS (MiB) |
|---:|---:|---:|---:|---:|
| 64 | 0.008320 | 0.009316 | 0.893x | 6.54 |
| 216 | 0.021005 | 0.019020 | 1.104x | 7.27 |
| 512 | 0.068522 | 0.046147 | 1.485x | 10.07 |
| 1000 | 0.219021 | 0.124789 | 1.755x | 18.05 |

| N | A1g t1/t16 | A2g t1/t16 | A1u t1/t16 | A2u t1/t16 |
|---:|---:|---:|---:|---:|
| 64 | 0.893 | 0.909 | 0.942 | 0.922 |
| 216 | 1.104 | 1.065 | 0.968 | 1.045 |
| 512 | 1.485 | 1.305 | 1.233 | 1.322 |
| 1000 | 1.755 | 1.595 | 1.550 | 1.590 |

### Main algorithmic interpretation

At N=1000, the anchored-cubic basis construction takes 0.0636-0.0639 s
across all four O_h irreps, **independently of their different projected
dimensions** D=4299,3624,2884,3252. For T=16 this is 51.1-53.8% of
the median native total. The code iterates over at most
`choose(N-1,2)` translation-anchored pair keys (N1000: 498,501),
and uses 48 O_h operations to canonicalize previously unseen orbits.
The equality of timings is consistent with a cost paid before
irrep-specific output sparsity is realized. It is not by itself
a formal empirical estimate of asymptotic complexity, and it should
not be interpreted as a claim of perfect thread-independent work.

At N1000, 1-to-16 speedup is approximately 1.55-1.76x across sectors.
At N64 all four are **slower** using 16 threads, as expected when
overhead dominates sub-10ms tasks. The maximum reported direct Ritz
residual in these 32 rows is **9.496e-9** (N1000/A2u/T1), within the
1.01e-8 audit threshold but relatively close to it. Avoid comparing
raw Lanczos walltimes across sectors without recording each case's
Hamiltonian applications and stopping behavior.

The native phase-total is not PBS turnaround and not necessarily the
sum of individually selected phase medians. If presenting a stacked bar
from medians, label its residual as **approximate**
`median(total) - sum(median(phases))`, not a measured per-run phase.

### Publication figures and reproducibility

Four draft figures were generated from the frozen 32-line aggregate
CSV in a separate user-transfer ZIP:
1. total median phase time vs N for each irrep, T=1 and T=16;
2. approximate A1g phase decomposition vs N;
3. 1-to-16 median-total speedup by sector vs N;
4. process RSS vs N by sector and thread count.

The ZIP contains PNG (300 dpi), SVG, Python plotting source, source CSV,
and derived-statistics JSON. These figures have **no error bars** until
the original per-trial data provide MAD or other variability estimates.

### Priority follow-up

1. Archive and SHA256 the original `24608.jpmjpr1911` job folder, with
   `audit/audit.json`, `measurements/trials.jsonl`, PBS pre/postflight,
   source and binary identities, and all N64 forced-legacy controls.
2. Export per-condition per-trial total, phase timings, residual and
   `h_applications` so valid MAD/error bars and solver-work normalization
   can be plotted without re-running PRESTO.
3. Obtain site-operator attestation for node exclusivity. The current
   `review_eligible=false` must not be manually overridden.
4. Matched cubic N64/N216 QuSpin (same J1/J2 XXZ and Gamma/O_h A1g,
   same energy/residual/physical sector) before attempting N1000. If
   QuSpin exceeds a declared budget, report timeout/right-censoring,
   not an invented runtime or a speedup.
5. Preserve general ED2 (noncubic, mixed-spin/non-U(1)/DM) validation;
   the cubic kernel is a performance case study, not the entire CPC
   program-library contribution.

This experimental branch is stacked on unmerged PR #21. Do not merge it
into main solely because the reported PRESTO numerical audit succeeded.
