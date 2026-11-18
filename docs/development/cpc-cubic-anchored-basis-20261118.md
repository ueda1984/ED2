# Exact three-magnon cubic Oh basis construction by translation anchoring

**Status:** experimental opt-in Draft PR. General-purpose ED2 remains unchanged when disabled.

## Measured motivation

PRESTO PR #12, N=1000 spin-1/2 J1-J2 XXZ Gamma/Oh, 8-thread basis phase:
- certified point group: 0.000220 seconds
- complete raw orbit scan: 2.102832 seconds
- basis total: 2.103429 seconds

The expensive scan is over C(1000,3)=166167000 raw spin configurations.

## Exact mathematical construction

A translation orbit of any triple of distinct occupied sites can be encoded by translating any one of the spins to site zero. The remaining two form a sorted pair (a,b), with 0<a<b<N. Therefore the candidate quotient size is C(N-1,2)=498501 for N=1000. A triangular visited array marks the 48 point transforms, each with 3 possible choices of zero anchor, as a complete affine orbit.

**Preserve the historical ED2 representative.** The anchored key is not necessarily the old minimum-ranked representative. For each of 48 point transforms it is sufficient to test at most 3^3=27 translations, each placing the minimum occupied x, y, and z coordinate at zero. Proof: if all occupied x coordinates were positive, translating all sites by minus one x step would decrease all sorted numeric site indices and thus the combinatorial rank without wrapping. The same argument holds separately for y and z. The old colex minimum must therefore have an occupied site at coordinate zero in each axis; the minimizer appears among these 27 shifts.

**Preserve the historical normalization bits.** Given the old representative, for each point transform find its stabilizing translations by mapping its first transformed occupied site onto each of the original three sites, and checking the full set. There are at most three candidates per point transform. At Gamma their translation characters are exactly one, so add the original ED2 complex character term once per stabilizing translation, in the *same* ordered point normal form. Calculate the old expression ro=(abs(coef)**2)*|G|/count and sqrt(ro) in the same order. This reproduces the 64-bit IEEE representation of the old list_r, as well as the sorted integer list_s.

## Scientific scope and fail-closed conditions

Only explicitly enabled through ED2_CUBIC_ANCHORED_BASIS=1, which also requires ED2_CUBIC_FAST_GROUP_VALIDATION=1, with independent all-site generator checks and 48-element point character multiplication. Restricted to the canonical cubic 4,6,8,10 side, spin-1/2, exactly three down spins, Gamma, and real 1D A1g/A2g/A1u/A2u sectors. Any incompatible input is rejected; a forced legacy orbit builder bypasses this optimization.

The original generic builder is unchanged and is still the default. The separate ED2_CUBIC_ANCHORED_REP flag controls Hamiltonian-cache representative lookup, not basis generation.

## Independent initial checks (GNU Fortran 14, non-PRESTO)

- Native make check PASS 1/1 with optimization disabled.
- **16/16** distinct sector/size combinations: side=4,6,8,10 and all four irreps, one thread, complete byte-for-byte equivalence of every legacy list_s rank and floating-point list_r bit pattern. More 4-thread cases at sides 4,6,8 and an 8-thread side=10 A1g case matched.
- Full Lanczos physics matched at sides 4,6,8, all four irreps (energy, direct residual, H applications at one thread) and also at N=1000 A1g/A2g; N=1000 A1u/A2u full-Lanczos comparison is not yet complete, although their complete basis rank/norm files match exactly.
- For the local N=1000, one-thread base builder orbit scan was about 7–10 seconds, and the anchored builder about 0.14–0.22 seconds on the same non-exclusive local machine. These are noisy exploratory measurements and **NOT PRESTO speedup claims**.
- A separate Fortran test-only cubic_basis_driver avoids slow Lanczos runs when exhaustively comparing N=1000 basis ranks and norms; test generator, same input SHA, exact basis bytes and dimensions are checked independently.

## PRESTO acceptance plan

After the GitHub CI (gfortran-13 and gfortran-14) and independent operator allocation review, stage an immutable archive under a *new* directory. First run N=1000 Gamma/Oh A1g at 1,4,8 threads, MNTE64, using the two validated group/cache opt-ins plus ED2_CUBIC_ANCHORED_BASIS=1. Enable detailed basis timing for the initial pilot only. Check 4299-dimensional basis, exact pinned input and output hashes, native energy, direct residual, Hamiltonian application count, RSS and audit before comparing with PR #12. Do not overwrite older source stages, alter existing running jobs or merge any Draft PR.

One-shot PRESTO timings are exploratory, not publication eligible; warmup/repeats and actual site-operator exclusive-node attestation remain necessary.
