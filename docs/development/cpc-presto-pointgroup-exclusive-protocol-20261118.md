# PRESTO single-node exclusive ED2 point-group benchmark

**Status (2026-11-18): prepared for operator review and manual submission;
NO PRESTO job has been submitted by this workflow.**
Baseline is the development branch
`codex/ed2-cpc-noncommuting-orbits`, not `main`.
Results are engineering evidence until allocation exclusivity is independently
confirmed. No scheduler/queue changes, no forced pushes, no merge.

## Scientific question and comparison matrix

Use the J1-J2 XXZ spin-1/2 square-lattice model with three flipped spins,
periodic boundaries, fixed identical bonds, a deterministic initial random
seed (1984), the same Hamiltonian cache (`wk_dim=-1`, `MNTE=-1`),
`ALG=1` and opt-in historical two-stage reconstruction,
`ED2_LANCZOS_TWO_STAGE=1`. Direct Ritz residual <= 1e-8.

These are **two different scientific comparisons**:

1. Same sector, same \(H\), same thread count and input bytes:
   **forced historical orbit enumeration** (`ED2_FORCE_LEGACY_ORBITS=1`)
   vs **validated fast orbit construction** (unset). Compare the
   `basis`, `operator_cache`, `lanczos` and `total` stages; all
   input SHA256, sector dimensions, direct residual certification and
   energies must agree. This is the primary **algorithmic speedup** comparison.
2. \(\Gamma\)-translation-only vs \(\Gamma/A_1\), \(\Gamma/B_1\) and
   \(X/A_1\) compatible **1D** point-group sectors. Compare absolute
   runtime, phase and memory demand versus reduced Hilbert dimension.
   **Do not call this a same-problem isolated solver speedup.**
   \(X\) uses the \(C_{2v}\) little group, not \(C_{4v}\). The
   two-dimensional \(C_{4v}\) E irrep is still unsupported by the
   scalar-character input.

| Tier | Square lattice | OMP threads | Two warmups | Recorded repeats | Runs |
| --- | --- | --- | ---: | ---: | ---: |
| Small smoke (shared/debug, not timing evidence) | 4x4 only | 1 | 0 | 1 | 7 |
| Main exclusive | 12x12 and 16x16 | 1,2,4,8,16,32 | 2 | 7 | 756 |
| Extended exclusive, separate reviewed job | 20x20 | 1,4,8,16,32 | 2 | 7 | 315 |

In every square-lattice/thread/repetition combination there is one
translation-only fast run and three point-group sectors with two variants
each: **7 native ED2 runs**. Warmups are retained but excluded from
statistics. Within each matched pair, fast/legacy order alternates;
sector/thread order is shuffled with a fixed seed per repetition to
reduce thermal/dynamic-order bias. The initial model and solver inputs
remain identical.

Start with the main tier. Submit an extended tier only after reviewing
memory, job duration, residual convergence and sector dimensions from
the main run; the 12-hour PBS limit is a hard ceiling, not a target.
The runner refuses to overwrite output paths and stops on any failed
numerical check or elapsed-budget guard, retaining failed raw logs.

## Site allocation request

The PBS file `tools/presto_pointgroup_exclusive.pbs` requests:

```text
#PBS -q standard
#PBS -l select=1:ncpus=96:mem=64gb
#PBS -l place=exclhost
#PBS -l walltime=12:00:00
```

These values follow the prior PRESTO dedicated-node experiment. Do not
confuse requested `ncpus=96` with proof of exclusivity. Before
submission, a PRESTO operator must verify that `exclhost` is honored
for this queue, that a single 96-core node/64 GB is permitted, and
that available CPU affinity and `exec_vnode` match policy. The user
must deliberately pass `ED2_EXCLUSIVE_CONFIRMED=1` on submission
**only after** this site-policy review. A runtime `qstat -f` and
`PBS_NODEFILE` check then refuses obvious nonexclusive/multinode
placements. That still does not prove no other workloads share the node.

The script will not run on a login node: it requires `PBS_JOBID`,
`PBS_NODEFILE`, `qstat -f` containing `place=exclhost` and
a one-node `Resource_List.select`, the queue name `standard`,
and `exec_vnode`. It captures `qstat -f`, `pbsnodes -av`,
CPU topology, taskset, source hash, compiler, BLAS limits and
binary hash. These are read-only actions; it does not change scheduler
settings, terminate or cancel jobs.

The runner chooses physically distinct CPUs from the allocation's
allowed Linux affinity mask, pins the parent (inherited by children),
and uses `OMP_PLACES=cores`, `OMP_PROC_BIND=close`,
`OMP_DYNAMIC=FALSE`, and one BLAS thread. An operator may provide
`ED2_CPU_LIST` to explicitly choose CPU IDs when NUMA/topology
requires it, but the runner rejects CPUs outside the allocation and
hyperthread sibling duplication. Inspect `lscpu -e` first.

## Preparatory commands on PRESTO login node (NO job submission)

Use a reviewed checkout of the draft benchmark branch, with a clean
tracked worktree. The branch currently sits on top of
PR #4 `codex/ed2-cpc-noncommuting-orbits`, itself stacked on PR #3,
PR #2 and PR #1.

```sh
cd /path/to/ED2
git fetch origin codex/ed2-cpc-presto-pointgroup-benchmark
git checkout codex/ed2-cpc-presto-pointgroup-benchmark
git pull --ff-only origin codex/ed2-cpc-presto-pointgroup-benchmark
git rev-parse HEAD
bash tools/presto_pointgroup_stage.sh "$PWD" "$HOME/ed2-presto-stage-20261118"
(cd "$HOME/ed2-presto-stage-20261118" && sha256sum -c ed2-source.sha256)
```

The staging directory must not exist beforehand. The shell script
refuses dirty tracked files; `git archive` preserves only the exact
committed source. Confirm the printed source SHA matches the reviewed
pull request. Make sure the staging archive and results directory are
readable/writable from PRESTO compute nodes.

One can inspect the **planned** run count without launching any ED2
simulation using the command below; `--allow-shared` is a label required
by the argument parser and does not request an allocation:

```sh
python3 tools/cpc_presto_pointgroup_benchmark.py \
  --output /tmp/unused-plan-path \
  --source-revision "$(git rev-parse HEAD)" \
  --source-archive "$HOME/ed2-presto-stage-20261118/ed2-source.tar.gz" \
  --preset main --allow-shared --plan-only
```

## Manual submission only after site approval

The following is a template **not executed by ChatGPT**. Use absolute,
site-valid paths without commas or whitespace in the `qsub -v` values.
Do not submit until the site-policy review above is complete.

```sh
mkdir -p "$HOME/ED2_POINTGROUP_RESULTS"
qsub -v "ED2_STAGE=$HOME/ed2-presto-stage-20261118,ED2_RESULTS_ROOT=$HOME/ED2_POINTGROUP_RESULTS,ED2_PRESET=main,ED2_EXCLUSIVE_CONFIRMED=1" \
  tools/presto_pointgroup_exclusive.pbs
```

The PBS launcher checks source checksum, creates a never-overwritten
job-ID directory, extracts the source on the allocated node,
runs `autoreconf`, `./configure FC=gfortran FCFLAGS="-O2 -fopenmp"`,
`make -j4`, and `make check`, recording build logs. It then runs
the benchmark and an **independent offline audit**. It does not auto
requeue, submit another job, modify PBS/queues or merge the code.

If a main-tier run uses too much time or memory, retain its
`trials.jsonl`, raw logs and `summary.json` with
`status=incomplete`; do not present only its completed subset as a
seven-repeat study.

## Data audit and operator attestation

Every run stores a unique trial directory containing the *exact* ED2
input deck and bonds, stdout, eigenvalues, diagnostics, resource
measurements and SHA256 metadata. Output also includes:

- `allocation-preflight.txt`, `allocation-postflight.txt`, compiler/build log
- `measurements/environment.json`, `measurements/trials.jsonl`,
  `measurements/summary.json`, `measurements/trials/*`
- `audit/audit.json`: completeness checks, rehashed inputs, direct Ritz
  residuals, numerical pairwise equality, and median/MAD for each
  phase and peak RSS. Missing/duplicate trials are rejected.

The audit deliberately returns
`review_eligible=false` until a named site operator separately
reviews the *actual* allocation, competing jobs, topology, and
exclusive use. An optional later operator-attestation JSON can contain:

```json
{
  "job_id": "ACTUAL_JOB_ID",
  "host": "ACTUAL_COMPUTE_HOST",
  "verified_by": "NAMED_SITE_OPERATOR",
  "exclusive_confirmed": true,
  "evidence_reviewed": "Actual qstat/pbsnodes records and overlapping jobs checked"
}
```

These are **illustrative placeholders**, not an attestation; do not
submit or manufacture this document without actual operator review.
After human verification, rerun:

```sh
python3 tools/cpc_presto_pointgroup_audit.py \
  --root /actual/job/measurements --output /new/audit-reviewed \
  --source-archive /actual/stage/ed2-source.tar.gz \
  --exe /actual/job/source/source/QS3ED2 \
  --operator-attestation /actual/operator-attestation.json \
  --require-publication
```

Passing the gate establishes minimum input/completeness/attestation
conditions only; it is not a hardware-independent speedup guarantee
or manuscript approval. Report all individual trials, median and MAD,
CPU model/topology, total and phase times, memory, H applications,
compiler/LAPACK/BLAS, and whether different sectors target different
eigenvalues. No source/numerical policy change is part of these scripts.

## CI and external dependencies

The optional workflow
`.github/workflows/cpc-presto-pointgroup-smoke.yml` runs a
**shared-runner functional smoke**, not a PRESTO benchmark, and
checks that shared data fail the publication prerequisite. Only
Python's standard library is needed for benchmarking and audit;
`gfortran`, Autotools, BLAS/LAPACK and `/usr/bin/time` are build
and measurement dependencies. Optional NumPy/SciPy tests remain
separate. Do not interpret CI success as an exclusive-node result.
