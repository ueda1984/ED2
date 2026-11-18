# Cubic three-magnon exact dense direct-hop representative table (Draft PR #16)

Status: experimental explicitly opt-in direct-table optimization stacked on PR #15; the general ED2 Hamiltonian, eigen solver, input, original anchored-hash resolver and all other cubic fast paths are preserved unchanged by default. No automatic main merge.

## Why optimize here?

PRESTO N1000 A1g local cache profile (PR15) reported independent audit PASS for all 1/4/8 threads. It processed 231566 active spin exchanges for D=4299 symmetry representatives; every exchange called anchored_key_for_triple, which evaluates all 48 point operations and all 3 translation anchors (up to 144 candidate keys) and performs table lookup. The 7439 sampled calls accounted for 0.0393, 0.0421 and 0.0428 total sampled thread seconds at 1/4/8 threads. Naive extrapolated 1.224,1.312,1.333 *thread seconds* is close to the 1.227,1.322,1.341 thread-second sums of the bond-math loop. Sampling bias, timer overhead and thread-vs-wall distinctions mean these are qualitative hot-path clues, **not an exact timing decomposition**.

## Mathematical construction

At Gamma, translations act with character +1. A three-site configuration can be translated so that any one of its three occupied sites becomes the origin, leaving two DISTINCT non-origin sites a,b from 1..N-1. Thus only choose(N-1,2) anchored pair keys exist. For N=1000 this is choose(999,2)=498501, using two int32 arrays (representative ID and sign), 3,988,008 bytes of additional table memory (~3.8 MiB).

Use the triangular bijection slot=(max(a,b)-1)*(max(a,b)-2)/2 + min(a,b), with 1<=a<b<=N-1. For each projected basis representative from the unchanged list_s/normalization, apply the independently certified 48 O_h permutations; for each transformed triple, choose each of its three occupied sites as translation anchor. Insert the resulting pair key, the ORIGINAL projected representative index i, and the real one-dimensional character chi(g)=+/-1 of that point action. In a surviving orbit, collisions must agree in BOTH representative index and character, otherwise fail closed on an invalid stabilizer/orbit mapping. Unpopulated slots represent states that vanish in that irrep.

With 4299 A1g representatives at N1000, precomputation attempts 4299*48*3=619056 insertions, versus 231566*48*3=33,345,504 anchor candidate evaluations incurred by the old per-hop resolver for that single cache build. The new hot path selects the FIRST occupied spin as translation anchor, computes two translated site positions, performs one integer triangular-array lookup, and emits the original representative ID and exact chi(g) phase. Any other anchor yields the same character in the surviving 1D representation. No change is made to the physical exchange amplitude or Lanczos.

Original anchored-key lookup can still be used by leaving ED2_CUBIC_DIRECT_HOP_TABLE unset (default=0), and ED2_CUBIC_ANCHORED_VERIFY=1 can independently cross-check direct lookup results against the legacy full group representative and character.

## Fail-closed activation

ED2_CUBIC_DIRECT_HOP_TABLE=1 requires ED2_CUBIC_ANCHORED_REP=1, ED2_CUBIC_FAST_GROUP_VALIDATION=1, ED2_CUBIC_ANCHORED_BASIS=1, ED2_CUBIC_LOCAL_BOND_CACHE=1. Underlying certified setup further requires canonical 3D simple cubic periodic geometry, three spin-half down spins, Gamma, real one-dimensional O_h character, exactly 9N conserving real XXZ bonds with local degree=18. Unsupported models, missing flags, character collisions, table index overflow and incomplete orbit enumeration fail closed.

## Test and benchmarking plan

Exact **complete** cached Hamiltonian bit-for-bit proof: `tests/cpc_cubic_direct_hop_table_regression.py` runs native production basis+local cache twice with byte-identical input and forces `ED2_CUBIC_DIRECT_HOP_TABLE=0/1`. It compares the SHA256 of the full `cubic_cache_proof_dump.txt` containing every diagonal and every nonzero offdiagonal representative target and *both* IEEE-754 complex components, and checks origin/character/flags/census. N64/N216 all four irreps T1/T4, N1000 all four irreps T1/T4/T8. Small volumes also compare certified two-stage native Ritz energies and direct residuals. Explicit negative tests ensure direct-hop without the prerequisite local cache fails. Existing PR14 and PR15 Hamiltonian, group, basis, legacy chain and independent PRESTO audit tests must also continue passing under both gfortran13 and gfortran14.

The independent PRESTO benchmark/audit records the exact `direct_hop_table` boolean and native message `Representative lookup: cubic-direct-hop-pairs` on fast trials, with the flag cleared for forced-legacy controls. PBS wrapper validates flag 0/1 and ensures it does not leak into the generic baseline make check.

After both CI compiler matrices pass, stage a pinned immutable source tar and run a new PRESTO job at N1000 A1g T1/T4/T8 (MNTE64) with direct_hop_table=1 and **cache_profile=0**. Compare input SHA256 to PR14/PR15 corresponding unprofiled baseline, exact result energy and direct residual, Hamiltonian applications (not necessarily identical), phase/basis/total wall and peak memory. Record the new upfront table construction cost separately via optional profiler run; benchmark figures must come from unprofiled runs. Do not assert a performance gain in advance, because upfront table construction could dominate small workloads. Independent operator node-exclusivity approval and warmup/repeated samples remain required for publication claims.