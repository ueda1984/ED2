# Experimental exact three-site sorting network (stacked on PR #17)

Status: **experimental opt-in** / no claim of measured PRESTO improvement.

## Motivation and scientific invariant

On PRESTO, PR17 ON/OFF 48/48 same-input native physics checks passed,
no H-application count changes, and the anchored basis wall time fell
from 207–212 ms to 90–92 ms (2.26–2.34x; warmup-excluded median over
three measured repeats for each of 4 O_h irreps and 1/4/8 threads).
At T8 the PR17 basis still takes approximately 51–55% of overall
two-stage solver wall time (operator review eligibility remains false).

Within the basis orbit scan, every one of up to 27 coordinate translations
for each of up to 48 point operations sorts a triple of integer site
indices, then compares the sorted triple against the current canonical
colexicographic best. For N1000/A1g there can be approximately
4299 × 48 × 27 = 5.57 million sorting calls. This is a plausible
hot path; its actual contribution must be measured, not assumed.

## Exact opt-in

`ED2_CUBIC_BASIS_SORT3_NETWORK=1` (default 0) is **only allowed when
both** `ED2_CUBIC_ANCHORED_BASIS=1` and
`ED2_CUBIC_BASIS_COORD_CACHE=1`; all existing cubic Oh Gamma
geometry and three-down-spin guards remain in force.

The previous generic `insertion_sort(tmp,3)` is replaced **only** in
the 27-candidate canonical-rank search loop. For the three integers
`a=tmp(1), b=tmp(2), c=tmp(3)`, let
`lo=min(a,b), hi=max(a,b)`, then assign
`tmp(1)=min(lo,c)`, `tmp(2)=max(lo,min(hi,c))`, and
`tmp(3)=max(hi,c)`. This is a branch-free exact sorting network with
the same ordered three integer values for all `a,b,c`, including ties,
without integer addition or overflow concerns. Existing candidate
enumeration, colex compare/update, orbit marking, stabilizer complex
summation and normalizations are untouched. Disabled mode retains
the exact original insertion_sort.

A native marker `Basis canonical sorter: exact-three-site-network`
makes the opt-in observable.

## Validation

`tests/cpc_cubic_basis_sort3_regression.py` compares the complete
native production basis dump (all historic representative ranks and
`list_r` IEEE-754 bits), not sampled states, with ON vs OFF sort3
but coordinate optimization **ON on both sides**. It checks 4 1D Oh
irreps and side lengths 4,6,8,10 at threads 1/4, plus all 4 sectors
at side10 threads8: 36 basis pairs per compiler. For side4/6
threads1/4 and side10 threads1/4/8, it also compares SHA256 of the
**entire native full Hamiltonian cache dump** (28 pairs per compiler),
and side4/6 threads1 compare real native Ritz energy, direct residual
tolerance and H-application count (8 pairs per compiler). Missing
coordinate prerequisite is rejected. Previous PR16/17 tests remain.
CI targets gfortran13+14.

No speedup is implied merely from compilation or exact numerical
equivalence. After both CI compiler jobs succeed, run a controlled
PRESTO sort3-OFF/ON comparison on **the same PR head commit**,
with coordinate cache ON in both, all previous fast flags ON, MNTE64,
all four irreps, thread counts 1/4/8, warmup1+repeat3 minimum,
and exact trial-level physics/audit checks. Keep PR18 Draft; never merge
into PR17 or main without separate review.
