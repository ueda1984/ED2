# ED2 CPC / PRESTO 2026-11-19: automatic cubic MNTE and reproducible measurements

This package is a **manual PRESTO transfer bundle**, not a submitted job.
The experiments on this branch are **not merged into main**. The source
commit is fixed in `cubic_auto_mnte/source-commit.txt`; main and QuSpin
comparator code are separate pinned archives with their own commit IDs.
Nothing in the ZIP writes to the existing `$HOME/ED2_EXPERIMENT_RESULTS`.

## 1. Scientific result so far (user-provided PRESTO logs; not archived raw data here)

The user's previous PR #21 commit `6da9ba96fc708de18a7b43d1db4ee818f8a8e7f4`
showed on PRESTO for periodic 10×10×10 spin-½ J1-J2 XXZ,
three down-spins, Gamma/O_h A1g, D=4299, E0=847.3078940655731:

| Configuration | Threads | Basis median (s) | Cache median (s) | Total median (s) | RSS (GiB) | Repeats |
|---|---:|---:|---:|---:|---:|---:|
| PR21 auto legacy MNTE=-1, unoptimized basis | 1 | .206993 | 1.580768 | 1.923757 | 5.778 | 7 |
| PR21 same code, unoptimized basis | 4 | .206542 | .433169 | .697635 | 5.778 | 7 |
| PR21 same code, unoptimized basis | 8 | .209264 | .240903 | .494663 | 5.778 | 7 |
| PR21 same code, unoptimized basis | 16 | .207710 | .185130 | .433110 | 5.778 | 7 |
| PR21 same code, unoptimized basis | 32 | .206488 | .155882 | .407305 | 5.778 | 7 |
| PR21 coordinate cache + canonical pruning | 32 | .063892 | .155910 | .266326 | 5.778 | 7 |
| PR21 coordinate cache + pruning + explicit MNTE=64 | 32 | .064493 | .021462 | .131079 | .017872 | 7 |

All reported runs had `AUDIT PASS`, same expected energy, and the
one-thread/32-thread Ritz residual checks were below the `1.01e-8`
acceptance threshold. **Their operator attestations are missing**.
The source logs above are user-reported, *not copied into this package*;
retain and independently verify the original PRESTO job directories.

The automatic-capacity change on THIS NEW BRANCH is designed to produce
the 64-transition storage whenever `MNTE=-1` AND the experimental
`ED2_CUBIC_LOCAL_BOND_CACHE=1` fixed-three-magnon path is selected.
The default generic ED2 path remains unchanged.

## 1a. Completed paired PRESTO N1000 run (2026-11-19; reported aggregates)

**Provenance:** User-provided native logs and summary for PBS
`24582.jpmjpr1911`, not the raw job archive. The candidate source
was from the **already transferred** 2026-11-19 ZIP; do not infer that
later documentation-only repository commits were deployed to PRESTO.
All 5 thread configurations had 2 warmups and 7 measured trials per
capacity mode, i.e. 45 distinct trials per mode. Both independent
producer-data audits reported `AUDIT PASS`; paired comparison reported
`CPC_AUTO_MNTE_PAIRED_PASS job=24582.jpmjpr1911 trials=45 cases=5`.
Independent operator attestation is still missing (`review_eligible=false`).

For periodic N=1000, 3 magnons, Gamma/O_h A1g D=4299,
with basis coordinate-cache + canonical pruning + direct-hop table,
the **automatic MNTE=-1** measured native phase medians were:

| Threads | Basis (s) | Cache (s) | Lanczos (s) | Total (s) | RSS (MiB) | Max direct Ritz residual |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.063686 | 0.043116 | 0.086935 | 0.216955 | 17.820 | 2.244e-9 |
| 4 | 0.063730 | 0.025627 | 0.028661 | 0.140868 | 17.816 | 2.837e-9 |
| 8 | 0.063865 | 0.022498 | 0.018013 | 0.128415 | 17.867 | 9.568e-9 |
| 16 | 0.063476 | 0.021243 | 0.015358 | 0.124013 | 18.090 | 4.918e-9 |
| 32 | 0.063512 | 0.020985 | 0.017671 | 0.129892 | 18.371 | 4.640e-9 |

Paired **explicit MNTE=64** median native total times for the same
PBS job were T1=0.220445s, T4=0.143192s, T8=0.129012s,
T16=0.125475s, T32=0.130942s. Across these thread counts the
automatic/explicit total median ratio was 0.9838-0.9954. Do not claim
automatic MNTE is faster without a paired significance/ordering
analysis; the scientific claim is **equal physical results with
automatic safe capacity**.

At this size 16 threads minimized the *native ED2 phase-total*
(0.124013s); native 1-to-16 speedup 1.75x (not 16x) because the
near-constant ~0.064s basis phase accounts for 51.2% of the T16 total.
At T8 the maximum reported Ritz residual 9.568e-9 is below the
1.01e-8 audit gate but close enough to disclose explicit thresholds
and sensitivity. The number of Hamiltonian applications may vary with
thread count; normalize solver-cost interpretation accordingly.
RSS remains 17.8-18.4 MiB. These are Fortran phase-total wall times,
**not** whole PBS turnaround times.

**Publication gate still open:** retain original
`$HOME/ED2_CPC_AUTO_MNTE_RESULTS/24582.jpmjpr1911/` including
`measurements/trials.jsonl`, `measurements-explicit64/trials.jsonl`,
both `audit.json` files, `paired-comparison.json`, PBS node
preflight/postflight, compiler/process affinity and binary/source SHA
records. Archive it with a SHA256 manifest. Independently obtain
site-operator confirmation of scheduler exclusion before publication.

## 2. Unpack and verify (no scheduler action)

Transfer this ZIP to PRESTO via your normal method and unzip it into a NEW
location. For example, let it unpack to

```bash
BUNDLE="$HOME/ED2_CPC_AUTO_MNTE_PRESTO_20261119"
cd "$BUNDLE"
sha256sum -c DELIVERY_SHA256SUMS
(cd cubic_auto_mnte && sha256sum -c ed2-source.sha256)
cat cubic_auto_mnte/source-commit.txt
cat latest_main_snapshot/main-commit.txt
cat cross_library_comparator/comparator-commit.txt
```

Check the candidate source includes the NEW automatic-capacity marker:

```bash
tar -xOzf "$BUNDLE/cubic_auto_mnte/ed2-source.tar.gz" source/ham2vec.f90 \
  | grep 'Cubic cache capacity: auto-MNTE=64'
```

Do not describe the main archive as including unmerged cubic features.

## 3. Next priority: paired auto/explicit64 + 1/4/8/16/32 threads, SAME node

Before submission an independent PRESTO operator must approve the
`standard` queue's `place=exclhost` contract and verify the operator
allocation record. **Setting ED2_EXCLUSIVE_CONFIRMED=1 is a user's
acknowledgment, NOT independent attestation**.

This command *manually submits one PBS job*. The auto and explicit64
comparators run sequentially within the SAME node, using the exact SAME
binary. Each mode uses 2 warmups and 7 measured repeats at each of five
OpenMP thread counts; individual configurations are seeded and shuffled.

Because `qsub -v` splits comma-delimited environment variables,
**use plus signs in `ED2_THREADS=1+4+8+16+32`**. The new PBS script
converts the plus signs back to the benchmark's comma-separated format.

```bash
BUNDLE="$HOME/ED2_CPC_AUTO_MNTE_PRESTO_20261119"
RESULTS="$HOME/ED2_CPC_AUTO_MNTE_RESULTS"
mkdir -p "$RESULTS"

JOBID=$(qsub -v \
"ED2_STAGE=$BUNDLE/cubic_auto_mnte,ED2_RESULTS_ROOT=$RESULTS,ED2_PRESET=cubic-main,ED2_SIZES=10,ED2_SECTORS=oh-a1g,ED2_THREADS=1+4+8+16+32,ED2_WARMUPS=2,ED2_REPEATS=7,ED2_TIMEOUT=1800,ED2_MNTE=-1,ED2_CPC_PAIRED_MNTE=1,ED2_CUBIC_FAST_GROUP_VALIDATION=1,ED2_CUBIC_ANCHORED_REP=1,ED2_CUBIC_ANCHORED_BASIS=1,ED2_CUBIC_BASIS_COORD_CACHE=1,ED2_CUBIC_BASIS_PRUNED_CANON=1,ED2_CUBIC_LOCAL_BOND_CACHE=1,ED2_CUBIC_DIRECT_HOP_TABLE=1,ED2_CUBIC_SKIP_UNUSED_HASH=1,ED2_EXCLUSIVE_CONFIRMED=1" \
"$BUNDLE/cubic_auto_mnte/presto_cubic_exclusive.pbs")

echo "PRESTO job: $JOBID"
qstat -f "$JOBID"
```

The scheduler script refuses missing PBS variables and refuses the
`ED2_EXCLUSIVE_CONFIRMED` guard if not explicitly set. The archive
checksum is checked, production source rebuilt, and `make check` runs
on the allocated node before any benchmark. Other jobs are not touched.

### Result extraction

```bash
RUN="$RESULTS/$JOBID"
tail -n 20 "$RUN/benchmark-console.log"
tail -n 20 "$RUN/audit-console.log"
tail -n 20 "$RUN/benchmark-explicit64-console.log"
tail -n 20 "$RUN/audit-explicit64-console.log"
cat "$RUN/paired-comparison-console.log"
python3 -m json.tool "$RUN/paired-comparison.json" | head -n 100
grep -ahE "Cubic cache capacity: auto-MNTE|Current MNTE|Optimal MNTE" \
 "$RUN"/measurements/trials/*/stdout.log | head -n 12
```

Both audits should PASS and `paired-comparison.json` should say
`paired_audit=PASS`. `review_eligible=false` is EXPECTED until
independent operator attestation is supplied; do not change it manually.

The original PR21 `MNTE=-1` run used a 72,001-entry-per-column
generic capacity. The new automatic policy uses 65 slots per column.
Auto and explicit64 at the NEW commit are primarily a **correctness
and reproducibility control**, not a direct proof of speedup over old code.
For a publication-quality speedup claim, compare against the pinned
OLD PR21 source on the same node under a controlled experimental design.

## 4. CPC / Computer Physics Communications calculation programme

### Tier 0: scientific correctness

* Exact 3D geometry and 9000 J1/J2 XXZ bonds; use the SAME source-input
  SHA and independent physical model SHA across libraries.
* Burnside independent symmetry dimensions, 48 point operations,
  Gamma/O_h A1g, A2g, A1u, A2u for N=64,216,512,1000.
* Full native bit-for-bit projected Hamiltonian cache
  (auto-MNTE=-1 vs explicit64) for all irreps and representative sizes.
* Certified two-stage Ritz pairs: reference energy, direct Hψ-Eψ
  residual <=1.01e-8, fixed seed, H application counts. Negative
  tests for other magnon count/non-XXZ/invalid degree. General spin
  and non-U(1) ED2 regression must remain unchanged.

### Tier 1: complete run vs separate phase costs

* N=1000 A1g, 1/4/8/16/32 threads **in one allocation**; 2 warmups,
  7 repetitions, medians and MAD. Record `basis`, `operator_cache`,
  `lanczos`, `total`, external RSS, direct residual, solver H calls.
* Ablation: original PR21 default capacity vs explicit64;
  fixed basis options A(no coord/no prune), B(coord), C(coord+prune);
  direct pair table with/without unused hash. Pin individual commits
  when results come from different versions.
* Record row-cache memory law: `20*(MNTE+1)*D` bytes for int32 + complex128
  (plus real program overhead), and check against actual RSS.
* Whole-job `make check` overhead is *not included* in native
  `ED2_PHASE total`; do not call native timing whole PBS turnaround.

### Tier 2: volume, symmetry, and competitors

* N=64,216,512,1000 × four 1D Oh sectors, thread count 1 and 8
  where feasible; stay within memory/timeout constraints, recording
  exact symmetry-reduced dimension for each case.
* QuSpin 1.0.1: first prove matched A1g Hamiltonian, basis, energy,
  and numerical residual at N=64/216; attempt N=1000 only if it
  completes within the declared memory/time budget. Do not imply
  an exact N=1000 result already exists; earlier CI timed out.
* XDiag: use comparable U(1) or Gamma-only sectors if exact A1g
  capability is absent, but label these **distinct physics sectors**.
  Never interpret cross-sector runtime ratios as speedup of ED2.
* On the same node, capture compiler/OpenMP/BLAS, CPU affinity,
  peak RSS, version/commit, all logs and scheduler evidence.

### Tier 3: generality expected of a CPC program-library manuscript

The cubic specialization is a *case study*, not the whole contribution.
Retain and benchmark the generic ED2 path on independently soluble spin
models: spin-1/2 U(1) chain/ladder, mixed-spin bounded bases,
non-conserving fields, real/complex anisotropy and DM interactions,
Krylov degeneracy/restart, MPI/OpenMP if actually supported. This
supports a program-library claim rather than only one engineered lattice.

For any reported speedup, specify benchmark denominator and exact
equivalent physical task; provide full raw per-run data, not just medians.

## 5. Publication and scheduler safeguards

* GitHub Actions shared runners **cannot** certify PRESTO performance.
* `place=exclhost` is a scheduler request, not independent evidence.
* Native direct residuals are physically meaningful and essential.
* CI cache-proof SHA, finite-group dimensions, native source SHA and
  node PBS provenance are separately retained.
* A draft manuscript is not publication-ready until software license,
  references, install tests, program summary, frozen reproducibility
  scripts and comparisons are editorially checked.
* **Do not merge any experiment into main until all correctness gates pass**.
