# Experimental cubic three-magnon local-bond Hamiltonian cache

Status: independent, opt-in, fail-closed Draft branch stacked on PR #13. Existing general ED2 models and default algorithms are preserved.

## Empirical motivation

Previous N=1000 PRESTO PR13 work has 12/12 independently audited 1D Oh symmetry-sector/1-4-8-thread runs, with old/new input SHA256 equality, matched energies and certified Ritz residuals. At eight threads, wall times are A1g=0.4606s, A2g=0.4212s, A1u=0.3828s and A2u=0.4057s. For A1g, Hamiltonian cache construction is 0.2328s and anchored basis construction 0.1860s. At one thread the cache construction consumes 1.5553s of the 1.8519s job.

The generic make_wk_loc_and_ele visits every one of 9000 interaction records for every projected state (38.7 million visits at D_A1g=4299), but exactly three down spins can participate in offdiagonal real XXZ exchange. Each physical site has 18 incident bonds: six nearest-neighbor and twelve face-diagonal neighbors.

## Mathematically exact sparse reduction

For spin-half and occupation n_i in {0,1}, Sz_i Sz_j=(1/2-n_i)(1/2-n_j)=1/4-(n_i+n_j)/2+n_i*n_j. Consequently E_diag = E_all_up - (1/2) SUM_over_bonds_with_one_down_endpoint Jz_b, where E_all_up=(1/4) SUM_all_bonds Jz_b.

A single setup scan over the actual p_two/Jint arrays constructs each site's sorted incident bond IDs, requiring degree=18 per site. Per representative, union/uniquify at most 3*18=54 bond IDs, sorted by original input bond index. For an active bond with exactly one occupied endpoint, apply the identical old j_flip_ni and resolve_cached_hop functions, and the *same* coefficient and list_r arithmetic order as the original XXZ branch. All nonzero offdiagonal cache indices, entries and their order are therefore intended to match the old code byte-for-byte.

Numerical qualification: regrouping the 9000 diagonal contributions into the once-computed all-up baseline plus at most 54 corrections changes the floating-point summation order. Therefore the diagonal is *mathematically exact* but need not have the same IEEE-754 bits. Independent native cache checks require absolute diagonal differences <=1e-8 per projected state and identical complex bits for all offdiagonal entries; independent Ritz energy/direct-residual acceptance is also mandatory. We do not claim bitwise diagonal equality.

## Scope and guards

ED2_CUBIC_LOCAL_BOND_CACHE=1 explicitly enables the new path. It requires all three existing opt-ins ED2_CUBIC_ANCHORED_REP=1, ED2_CUBIC_ANCHORED_BASIS=1 and ED2_CUBIC_FAST_GROUP_VALIDATION=1, which certify the canonical even cubic Gamma/O_h geometry, its exact 48-point semidirect character and the correct projected basis.

Additional fail-closed preconditions: exactly three down spins, spin-half, no fields, NO_two=9*NOS, degree 18 on every site, sorted distinct endpoints, exactly Jx==Jy (real XXZ), all Jint(4:9)=0, and complete projected wk_dim==lv. This excludes any general XXZ ladders, fields, anisotropic transverse exchanges, DM or other interactions that could alter the Hamiltonian. Feature default is disabled; original generic Hamiltonian and matrix-free path remain unchanged. PRESTO PBS preflight logs the flag, and only the generic make-check build subshell unsets it.

## Independent science checks

The new test-only native cubic_cache_driver constructs an actual production basis and cache without Lanczos. Opt-in ED2_CUBIC_CACHE_PROOF_DUMP=1 emits all diagonal and nonzero offdiagonal cache records with their full IEEE real64 bit patterns. Independent Python regression runs original full-bond scan and new local-bond scan on precisely the same input and checks:

1. Input SHA256, native output mode, projected dimension, representative resolver and anchored basis all match.
2. Every offdiagonal target state index and real/imag IEEE bits match exactly, in original input-bond order, with no duplicate or missing entry.
3. Every diagonal entry is within absolute 1e-8; maximum deviation is reported.
4. Real Lanczos matching energies within 1e-8 and direct residual <=1.01e-8 at sides 4 and 6 in all four 1D sectors.
5. Negative tests on unrelated spin chains or invalid Jx!=Jy bonds fail closed.

The PRESTO benchmark/audit schema records the local_bond_cache flag and requires the native local-bond marker for its fast variant while retaining the original full-bond-scan legacy comparison.

## PRESTO next step

After gfortran13/14 CI and exact source pinning, make a new immutable source stage. On an operator-approved exclusive PBS allocation, run N1000 A1g 1/4/8 threads with MNTE64 and the three earlier opt-ins plus local_bond_cache=1, compare every input SHA256 against PR13, and report Ritz E, direct residual, H application count, cache and total wall times, and RSS. The cache proof dump must stay off for timings. Then test all four Oh irreps and reproducibility with repeated samples and operator attestation.

Neither this Draft PR nor any performance investigation warrants changing main or claiming a publication-grade speedup before those independent checks.