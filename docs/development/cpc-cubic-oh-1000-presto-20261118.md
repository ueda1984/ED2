# PRESTO cubic (O_h) 1000-site benchmark — staged protocol (2026-11-18)

**Status: prepared, pending GitHub CI and independent PRESTO operator approval.
No PRESTO job has been submitted.** This is a separate stacked draft change on
top of the PR #5 square-lattice protocol. It does not alter ED2's default
Lanczos mode, normal-form limits, queue configuration, or `main`.

## Why a 3D simple-cubic lattice, not a 1000-site ring?

The previous 1000-site ring is a useful established performance and solver
baseline, but its nontrivial point operation is essentially inversion. It
does not exercise three translations and the noncommuting cubic (O_h)
point group. A periodic `10 x 10 x 10` spin-1/2 lattice provides **both**
a 1000-site scale and a meaningful 48-element cubic point group.

The current ED2 input has exactly six cyclic generator slots. A checked
normal form for all 48 signed coordinate permutations uses three point
generators in the remaining slots:

- `Tx, Ty, Tz`, order `(L,L,L)`.
- `A = I*C4x`, mapping `(x,y,z)->(-x,z,-y)`, order 4.
- `B = C4z^-1`, mapping `(x,y,z)->(y,-x,z)`, order 4.
- `C = C3[111]`, mapping `(x,y,z)->(z,x,y)`, order 3.

Explicitly, all products `A^a B^b C^c` for
`0<=a,b<4`, `0<=c<3` give **48 different** point maps.
Together with translations they enumerate `L^3*48` affine elements,
including `48,000` elements at `L=10`.
The code verifies uniqueness, point group closure and every 1D character
against multiplication, rather than assuming the normal form is valid.

Only the four 1D representations `A1g,A2g,A1u,A2u` can currently be
specified via ED2's scalar phases. Their `(M4,M5,M6)` are respectively
`(0,0,0), (2,2,0), (2,0,0), (0,2,0)`. This does **not**
support the two- or three-dimensional `O_h` irreps; no claim of a
general non-Abelian representation engine is made.

## Physical model and independently computed dimensions

Periodic 3D simple cubic, spin-1/2, exactly three down spins:
`H = sum_b Jxy_b (Sx_i Sx_j+Sy_i Sy_j)+Jz_b Sz_i Sz_j`.
Every positive axial direction uses `(Jxy,Jz)=(1,0.77)`.
All six distinct face-diagonal directions use `(0.31,0.19)`.
This is `3N + 6N = 9N` unique physical bonds and is invariant under
every cubic point operation (including inversion) and translation.

| Side | Sites | Raw 3-down dimension | Gamma translation | Gamma A1g | Gamma A2g | Gamma A1u | Gamma A2u |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 64 | 41,664 | 651 | 44 | 16 | 4 | 9 |
| 6 | 216 | 1,656,360 | 7,677 | 289 | 177 | 95 | 134 |
| 8 | 512 | 22,238,720 | 43,435 | 1,271 | 967 | 679 | 820 |
| 10 | **1000** | **166,167,000** | **166,167** | **4,299** | **3,624** | **2,884** | **3,252** |

`tests/cpc_cubic_pointgroup.py` obtains these dimensions independently
using the fixed-three-subset character projector. For each affine group
element `g`, its fixed 3-spin configurations are counted from the
numbers of fixed 1-, 2- and 3-cycles:
`C(f1,3) + f1*f2 + f3`. The sums are divided by group order and
checked for integrality; no raw 166-million-state enumeration is
needed for the dimension oracle. All transformed bond tables are
checked to remain exactly equal. The 1000-site figures are **group
dimension predictions**, not a claim that native ED2 on PRESTO has run.

### Numerical evidence before PRESTO

Local GNU Fortran 14 reference:
- `4^3` all four 1D sectors: independent sparse-bitmask Hamiltonian
  and character-projected full spectra, with all native ground-state
  Lanczos energies and direct residuals matching fast and forced
  historical builders. These are verified in
  `tests/cpc_cubic_pointgroup_regression.py` using NumPy for the
  small reduced eigenproblem.
- `6^3` A1g and A2g native three-magnon runs: with opt-in large
  group validation, both passed direct residuals `<1e-8`.
- `8^3` and `10^3` full group, J1-J2 bonds and sector dimensions
  have independent algebraic tests; full native Lanczos on those
  geometries remains a PRESTO objective.

## Bounded group validator — necessary 1000-site opt-in

The previously reviewed noncommuting fast orbit validator conservatively
limits groups to 8192 elements and four million stored site images.
For a `10^3` cubic `O_h` Gamma sector the normal form has
`1000*48=48000` elements and `48000*(1000+1)=48048000`
integer site images (about **183 MiB** in 32-bit storage).

The separate development change allows an explicit
`ED2_ALLOW_LARGE_GROUP_VALIDATION=1` opt-in. The extended hard
bounds are **65536 group elements / 64 million site images**;
they do not bypass permutation uniqueness, group closure or character
homomorphism checks. Default limits remain unchanged. The cubic
benchmark driver requests the opt-in, and checks ED2 output for a
genuine validated fast-basis path. If the group fails validation or
resource ceilings, it is **not silently benchmarked as fast**.

## Exclusive-node measurement tiers

The new driver is `tools/cpc_presto_cubic_benchmark.py`. It reuses the
source-provenance, affinity, qstat and direct-residual acceptance
protocol of PR #5 while checking the cubic model independently.

| Tier | Geometry | Default threads | Warmups | Measured per case | Trials |
| --- | --- | --- | ---: | ---: | ---: |
| `cubic-smoke` | 4³ | 1 | 0 | 1 | 9 |
| `cubic-pilot` | 4³ + 6³ | 1,4 | 0 | 1 | 28 |
| `cubic-main` | **10³** | 1,4,8,16,32 | 1 | 3 | 100 |
| `cubic-publication` | **10³** | 1,8,32 | 2 | 7 | 135 |

The publication tier can be scoped with
`ED2_SECTORS=gamma-t,oh-a1g` (then 54 trials).
**First submit the pilot**, independently inspect allocations and
runtime, then submit a short single-configuration 10³ run before the
full main/publication tier. The 12-hour queue limit is not a license
to launch an unmeasured all-pairs grid.

The prior legacy basis enumerator is paired against the accelerated
basis for the *same physical sector* only at `L=4`. For 216–1000
sites, the historic fallback is potentially prohibitively expensive;
those tiers intentionally run **fast only** and compare absolute costs
of different symmetry sectors. They **do not** establish a same-sector
algorithmic speedup at 1000 sites. Dedicated larger legacy comparisons
would require a separate bounded pilot and honest timeout reporting.

Report exact phase wall times for basis, operator cache, Lanczos,
and total; peak RSS and H-application counts, binary/source SHA256,
all exact physical inputs, per-run Ritz energies/residuals and
`median +/- MAD`. Store all warmups and every failed trial, too.
Direct residual certification does not prove global lowest-state
inclusion.

## Reproduction and manual PRESTO deployment

The old source archive and 2D protocol are **not suitable** for cubic
benchmarking; stage the current full reviewed cubic branch instead.
The staging script is already in PR #5:

```sh
cd /path/to/ED2
git fetch origin codex/ed2-cpc-cubic-oh-1000-presto
git checkout codex/ed2-cpc-cubic-oh-1000-presto
git pull --ff-only origin codex/ed2-cpc-cubic-oh-1000-presto
git rev-parse HEAD
bash tools/presto_pointgroup_stage.sh "$PWD" "$HOME/ED2_CUBIC_OH_stage"
(cd "$HOME/ED2_CUBIC_OH_stage" && sha256sum -c ed2-source.sha256)
```

The stage directory must not already exist; archive files are
write-once. Do not run the benchmark on a login node.

**Site policy and exclusivity must be verified before submitting:**
`standard`, `select=1:ncpus=96:mem=64gb`, `place=exclhost`,
`walltime=12:00:00` are requests, not proof of actual allocation.
Only after a PRESTO operator has confirmed node-wide exclusivity,
CPU mask, fair-share policy and other overlapping jobs should the
user manually execute:

```sh
cd /path/to/ED2
mkdir -p "$HOME/ED2_CUBIC_OH_results"
qsub -v "ED2_STAGE=$HOME/ED2_CUBIC_OH_stage,ED2_RESULTS_ROOT=$HOME/ED2_CUBIC_OH_results,ED2_PRESET=cubic-pilot,ED2_EXCLUSIVE_CONFIRMED=1" \
  tools/presto_cubic_exclusive.pbs
```

After the pilot has been audited, set `ED2_PRESET=cubic-main`
for 1000 sites. Override `ED2_THREADS`, `ED2_SIZES`, or
`ED2_SECTORS` to narrow the first large-node smoke.
The runner fails closed on missing scheduler records,
mislabelled basis paths, unequal paired results, incomplete
repetitions, nonconvergence and an elapsed budget guard.
Manual operator attestation remains required after the job.

The independently implemented offline auditor is
`tools/cpc_presto_cubic_audit.py`. It rechecks all exact
physical input hashes, dimensions, all native Ritz-pair residual rows
and failed/missing trial sets and explicitly refuses shared data for
publication qualification. It does not assert exclusive use by itself.

The GitHub Actions test
`.github/workflows/cpc-cubic-oh-presto.yml` uses shared runners
**solely for functional and mathematical verification**. It cannot
supply PRESTO timing evidence.

## Actual N=1000 native group-path smoke (one magnon, NOT 3 magnons)

The exact 1000-site cube with all 48,000 affine translations and O_h group
elements has already been exercised locally in the native ED2 executable
using **only one flipped spin**. The bounded extended validator accepted the
full normal form; the accelerated atomic-orbit basis was selected; the sector
dimension was exactly one. This analytically solvable state has

`E = (3N*0.77+6N*0.19)/4 - (6*0.77+12*0.19)/2 + (6+12*0.31)/2 = 863.91`

for `N=1000`, matching ED2 to `~6e-11`, with a directly observed
Ritz residual `3.43e-13`. Local Fortran phase timing for basis formation
was `3.46 s` on a shared, non-PRESTO one-thread machine. This **only**
demonstrates the 1000-site group-validation path and physical labeling.
It is **not** evidence for 1000-site three-magnon memory, solver scaling,
throughput or PRESTO performance. The native smoke is reproducible with
`tests/cpc_cubic_1000_one_magnon_smoke.py`.

## First full 1000-site three-magnon job (one configuration only)

After the small `cubic-pilot` allocation has completed successfully and
an operator has independently verified the exclusive-node policy, submit
**one** 1000-site A1g three-magnon calculation (not all 100 main trials),
using the exact staged cubic source:

```sh
qsub -v "ED2_STAGE=$HOME/ED2_CUBIC_OH_stage,ED2_RESULTS_ROOT=$HOME/ED2_CUBIC_OH_results,ED2_PRESET=cubic-main,ED2_SIZES=10,ED2_THREADS=1,ED2_SECTORS=oh-a1g,ED2_WARMUPS=0,ED2_REPEATS=1,ED2_TIMEOUT=1800,ED2_EXCLUSIVE_CONFIRMED=1" \
  tools/presto_cubic_exclusive.pbs
```

That job creates the real `D=4299` cubic three-magnon basis, performs
Lanczos, records the direct eigenpair residual and phase timings, and
passes its raw records through the independent audit (without an
unearned publication-ready status). The per-trial timeout is 1800 s;
nonconvergence must be retained as a failed trial, not dropped.
Only if this succeeds should the larger multi-OMP measurement be
scheduled. No command in this document has been executed on PRESTO.


## PRESTO runtime compatibility follow-up (2026-11-18)

The first PRESTO allocations ended **before** any cubic three-magnon
Lanczos timing was recorded. Keep both failed jobs as audit evidence:

- `19307.jpmjpr1911`: `configure` failed to find the unversioned
  BLAS/LAPACK development linker names even though an OpenBLAS runtime
  was available. The production Fortran sources were not at fault.
- `19353.jpmjpr1911`: with a private
  `libopenblas.so -> libopenblas.so.0` and
  `liblapack.so -> libopenblas.so.0` link workaround, `make check`
  reported **PASS=1, FAIL=0, ERROR=0**; the following independent
  Python geometry check then failed at `math.comb`, which was
  introduced in Python 3.8. The actual PRESTO minor Python version is
  not yet recorded here; query it with `python3 --version`.

The updated job script embeds the PRESTO private-link fallback (without
changing system libraries). The cubic geometry oracle now computes
integer `C(n,3)` directly and its shared profile helper computes an
integer binomial without `math.comb`; the auditor and 1000-site
one-magnon smoke no longer use Python 3.8 assignment expressions.
The runner uses pre-3.7 portable `subprocess.run` keyword arguments.
GitHub CI checks the Python **3.6 syntax grammar** and deletes
`math.comb` during oracle tests before repeating the physical
regressions on gfortran-13/14.

The source archive SHA and Git revision MUST be refreshed when using
this compatibility fix. Do not patch the original stage in place and
continue to identify the modified source by the old pinned commit.
Create a new stage/output path so neither the historical
`19307` nor `19353` evidence is overwritten. Re-run
`cubic-pilot`, then independently inspect native direct residuals,
actual allocation, logs, and audit before the first 1000-site
three-magnon job.
