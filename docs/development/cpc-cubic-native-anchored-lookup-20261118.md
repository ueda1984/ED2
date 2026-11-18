# Three-magnon cubic anchored-orbit representative lookup: native experiment (2026-11-18)

## Why this matters

Measured PRESTO 10×10×10, spin-1/2, J1–J2 XXZ, 3 magnons, Gamma/O_h A1g:
D=4299, E=847.3078940655731, direct residual 2.2401986516701506e-9,
170 Hamiltonian applications, and 78.233475885 s total with
`MNTE=64` on the user's PRESTO node. Of this, 65.595136833 s
(83.84%) is the `ham2vec::make_wk_loc_and_ele` phase.
RSS was 196420 KiB. That is a **single PRESTO run**, not replicated.

For each of 4299 representatives, the original cache loops over all
9000 physical bonds, with up to 54 allowed directed magnon hops.
The old `state_lists::representative_SQ` examines 10^3×48=48000
affine operations **per hop**, even for exactly 3 occupied sites at
Gamma. The point-group phase has to be preserved, not just the orbit.

## Implemented algorithm (opt-in only)

A distinct module `source/cubic_anchored_lookup.f90` is activated by
`ED2_CUBIC_ANCHORED_REP=1`. There are **no default behavior changes**
and no modifications to `state_lists::representative_SQ`,
`state_lists::try_orbit_basis`, the Lanczos solver or the Hamiltonian
coefficients.

1. Fail closed unless the input is the canonical 3D even cubic
   site-numbering and exact six transformations (Tx,Ty,Tz,
   inversion*C4x, C4z^-1, C3[111]), with Gamma translations,
   four real 1D O_h character choices, exactly three spin-1/2
   down spins, 9*N J1–J2 bonds, and real number-conserving XXZ.
   Generator permutations are checked **on every site**; the
   independent normal-form validation remains mandatory.
2. Precompute 48 point-group operations on physical sites and their
   irreducible characters, without storing 48000 affine site maps.
3. For any triple of occupied sites, apply the 48 signed point maps;
   for each transformed triple, translate one of its 3 sites to the
   origin. Lexicographically minimize the resulting sorted triple:
   **at most 3×48=144 candidates**. Record its real 1D character
   (±1). Character-inconsistent stabilizers signal forbidden sectors.
4. During setup, hash each key to the **existing** ED2 projected
   representative index `list_s(i)` and its sign relative to the
   key. The anchored triple is **not** used directly as an ED2
   basis representative; PR #9 gives a counterexample.
5. During cache building, compute the hop key, look up the
   preexisting representative index, and recover the phase via
   the ratio of the two real 1D characters. Invalid projected
   states remain absent. Only the cached ALG=1 hot path uses
   this optimization; all other operator paths are unchanged.

`ED2_CUBIC_ANCHORED_VERIFY=1` is a **test-only** setting: for
**every** generated cache hop, also invoke original
`representative_SQ` and `findstate` and require the same
representative index and the same complex character within 1e-10.
Do not set this during benchmark timings.

## Scientific checks (actually executed locally)

Using gfortran with `-O2 -fopenmp -fcheck=all -fbacktrace`, the
local build and `make check` passed 1/1. Optional LAPACK interface
and sparse-field operator targets also passed.

At 4×4×4 and 6×6×6, **all four O_h one-dimensional representations**
were compared with the default ED2 cache, with identical input deck
hash, projected dimension, ground energy within 1e-8 and certified
Ritz residual ≤1.01e-8. At one OpenMP thread, the instrumented
per-hop reference mode verified the **original representative index
and character for every generated cache element**. At four
threads the uninstrumented native path passed the same
energy/residual tests. Krylov H-application counts can differ
at four threads due to termination near the prescribed residual
threshold, so they are recorded but not forced equal there.

The original independent 4×4×4 dense projected-Hamiltonian
spectra (previous CI) agree with the legacy path. The PR #9
Python oracle exhaustively classified all 41664 three-down states,
recovering 44 affine orbits and four sector sizes (44,16,4,9).

A separate **local, non-exclusive computer** executed 1000-site A1g,
MNTE=64, one thread, with the accelerated cache:
- E=847.3078940655731 (matches original PRESTO energy exactly as printed)
- direct residual=7.64745229990986e-9 (passes 1e-8)
- H applications=160 (different compiler/system, not a mismatch gate)
- peak RSS=193004 KiB
- basis=11.765096665 s
- operator_cache=1.367447997 s
- total=13.295097245 s

**DO NOT divide the local 1.37 seconds by the PRESTO 65.60 seconds
and call it a measured machine-controlled speedup.** The local
legacy N=1000 run did not finish within the tool's execution limit.
The direct 1-thread matched-machine speedup must be measured on PRESTO.

## PRESTO controlled first pilot (requires verified exclhost policy)

Stage a **new immutable Git revision** and SHA256 for this PR; the
archive from PR #8 does **not** contain this Fortran module.

Use a new stage/source directory and preserve PRESTO jobs 19383 and
19447. From the source directory of this new stage, after independent
operator approval of `place=exclhost`, submit one 10³ A1g case:

```sh
JOB=$(qsub -v "ED2_STAGE=$ED2_STAGE,ED2_RESULTS_ROOT=$ED2_RESULTS,ED2_PRESET=cubic-main,ED2_SIZES=10,ED2_THREADS=1,ED2_SECTORS=oh-a1g,ED2_WARMUPS=0,ED2_REPEATS=1,ED2_TIMEOUT=1800,ED2_MNTE=64,ED2_CUBIC_ANCHORED_REP=1,ED2_EXCLUSIVE_CONFIRMED=1" tools/presto_cubic_exclusive.pbs)
echo "$JOB"
```

The audit checks the provenance mode against the **native stdout log**
as well as the original energy, binary/source/input hashes,
dimensions and certified direct residual. Confirm the residual
≤1e-8 and 4299-dimensional sector, then compare cache build time,
peak RSS and full elapsed time to a matched PRESTO legacy run.

The current pilot has zero warmups/one repetition and is
**not publication eligible**, even if independently audited.
Publication comparisons require at least 2 warmups, 7 repetitions,
controlled physical cores, environment attestation and independent
review. No qsub or queue policy changes are performed by this PR.

## Limits and follow-up

- Only gamma-point and four real 1D O_h sectors, 3 down spins and
  the canonical even cubic test geometry; no higher-dimensional
  point irreps, nonzero momenta, magnetic fields, multiple-magnon
  sectors or arbitrary input symmetry groups.
- The retained 48000-site-map validation and original basis
  construction still dominate memory (~190 MiB) and ~12 s of setup.
- Cache setup still scans all 9,000 bonds per projected state.
  A separate guarded adjacency-based builder is the next candidate,
  but requires additional independent matrix-element validation.
- The module caches lookup arrays once per executable process;
  reusable library contexts are not yet supported. This is an ED2
  executable pilot, not the final CPC library interface.
- PR remains **Draft**, stacked on PR #9; do not merge or report
  an accepted CPC improvement without independent review.
