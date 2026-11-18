# PR #17 PRESTO 1000-site pilot: source and limitation notes

Data table: `pr17_presto_coordinate_cache_pilot.csv`.

**Scope:** periodic 10 x 10 x 10 simple-cubic S=1/2
nearest-neighbor Heisenberg model, three down spins,
Gamma momentum and one of four **one-dimensional** O_h irreps:
A1g (D=4299), A2g (D=3624), A1u (D=2884), A2u (D=3252).
OpenMP threads 1, 4, 8; source-pinned local transition cache
and direct-hop lookup ON; MNTE=64.

**Provenance:** PRESTO group `jpmjpr1911`, two distinct job
allocations at the same source commit
`3fbe2d7dbb3314bbf57391565937b271c259ecb8`,
source-archive SHA256
`9b20299486c03a8ef5096e72bf1d9da80a718d474db21db320358554fbf173cc`.
- OFF: job `22282.jpmjpr1911`, `ED2_CUBIC_BASIS_COORD_CACHE=0`.
- ON: job `22322.jpmjpr1911`, `ED2_CUBIC_BASIS_COORD_CACHE=1`.
- All other native exact cubic optimization flags and case geometry
  are matched; each job ran one warmup plus three measured
  repeats per irreducible representation and OpenMP count.
- Both job-local independent audit programs reported
  `audited=True`, `status=complete` and `raw_trials=48`.
  They independently checked native diagnostic sidecars,
  physical residuals, phases, input and source provenance.
- The *paired user-executed* cross-job comparison asserted
  equal source/archive hashes, all 48 native inputs and
  generated input SHA256s, sector dimensions, energy
  differences <1e-8, direct residual <=1.01e-8 and
  Hamiltonian action counts. Its output:
  `PHYSICS CHECKS PASS: 48/48`;
  maximum `|Delta E|=0.00000e+00` at recorded precision;
  maximum Ritz residual `9.71161e-09`;
  differing Hamiltonian-application counts `0`.
- Every row in this CSV is a **median of the three measured
  trials only**, excluding the warmup. Units are milliseconds.
  Speedup is OFF time divided by ON time. The exact T8
  manuscript rows are guarded by
  `paper/check_pr17_presto_table.py`.

**Evidence boundaries:**
1. The raw PRESTO `measurements/trials.jsonl`, all 96 trial
   directories, `summary.json`, `audit/audit.json`,
   allocation preflight, build logs and operator attestation
   **have not been deposited in this repository**. This
   CSV is a transcription of the supplied output, not
   an independently reprocessed raw trial archive.
2. Independent scheduler/operator confirmation of exclusive
   resource use was **not provided** for either job;
   both `operator review eligible=False`.
3. Three measured repetitions and two allocations do not
   establish a publication-quality hardware speedup,
   statistical confidence interval or strong-scaling law.
4. The 48/48 equivalence validates the two coordinate-cache
   modes of **the same code release**. It is not an independent
   spectral reference for the full 1000-site three-magnon
   Hamiltonian. Independent orbit counting and a separate
   sparse Hermiticity/structural check for the translation-only
   Gamma sector are distinct evidence, not a replacement
   for an independent four-irrep energy calculation.
5. Do not label this table as a CPC-ready result or report
   the pilot speedups without the caveats above.

**Publication-evidence next actions:** preserve full raw
outputs outside the source repository until the approved
archive destination is chosen; verify allocations independently;
collect at least two fresh resource-attested allocations
with >=2 warmups and >=7 measurements per case (or
pre-register a different sample design), include a second
physical benchmark (mixed-spin anisotropic non-U(1) case),
and regenerate an immutable raw-to-table pipeline. Do
not create an official release, DOI or upstream merge
without author/maintainer approval.
