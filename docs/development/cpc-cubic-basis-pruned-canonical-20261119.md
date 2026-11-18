# Experimental exact axis-pruned cubic canonical enumeration (PR #19)

Status: Draft experimental optimization stacked on the **CI + PRESTO-validated**
PR #17 head `3fbe2d7dbb3314bbf57391565937b271c259ecb8`.
No change to PR #17, PR #16 or `main`. This branch has not yet
undergone native PRESTO testing.

## Motivation

PRESTO PR17 paired OFF job 22282 versus ON job 22322 yielded 48/48
matching native energies, direct Ritz residuals, input hashes and
Hamiltonian-application counts, and basis speedup of 2.26–2.34x
with coordinate caching ON. Four one-dimensional O_h sectors ×
OpenMP threads 1/4/8, with one warmup and three measured repeats each.
Yet the basis remains 51.0–55.1% of total at T8 (operator allocation
attestation still pending).

One previous micro-optimization (PR #18: 3-element sorting network)
is not portable: gfortran13 improved, while gfortran14 became slower,
despite exact physics. A more algorithmic reduction of the
canonical translation search is therefore preferable.

## Mathematical reduction (exact at all cubic sides)

For a point-transformed three-site state, every historical canonical
candidate has integer site indices
`r_j=1+x_j+L*(y_j+L*z_j)`, with each coordinate in
`0,...,L-1`, translated by independently choosing one of the three
sites to zero along x, y and z. There are at most 27 candidate triples.
Sort their site indices ascending, and choose the triple
lexicographically minimal when examined from largest to smallest
(the historical colex rank).

**Observation:** The largest site index `max_j r_j` can be minimized
with strict coordinate-wise bounds. For any candidate z-anchor,
let `Z=max_j z_j`. Because the lower axes contribute at most
`L^2-1`, any anchor with strictly larger `Z` must have strictly
larger `max_j r_j`; eliminate those anchors. For a retained
z-anchor, among sites achieving its maximum z, let
`Y=max y_j`. Any y-anchor with strictly larger Y is worse in
largest site index for this z-anchor, irrespective of x; eliminate
those. For a retained z/y-anchor, among sites achieving both
maxima, let `X=max x_j`, and eliminate x-anchors with larger X.
Critically, **retain every tie at all three stages**: competing
translations can yield the same largest site index but different
second/third values. Finally construct exactly the surviving
site-index triples, use the unchanged `insertion_sort` and original
colex update to find the historical canonical state.

The discarded candidates have provably worse largest site index
than another retained candidate, so none can minimize the complete
sorted colex tuple. This preserves the global result for every
configuration, including coincident coordinate components and ties.
No floating-point calculations are modified.

## Implementation and safety

`ED2_CUBIC_BASIS_PRUNED_CANON=1` (OFF by default) requires
`ED2_CUBIC_BASIS_COORD_CACHE=1` and the previously certified
three-down-spin anchored cubic Gamma/O_h group certificate.
Missing prerequisites/invalid flag values must fail closed.
The original 27 candidates remain executable when the flag is off.
Distinct native provenance line:
`Basis canonical enumeration: exact-axis-pruned`.

## Validation gates

The regression runs the actual native ED2 basis, Hamiltonian cache and
two-stage solver on **identical input decks** with coordinate caching
ON on both sides and pruning OFF/ON:

- All 36 complete basis pairs per compiler (L=4,6,8,10; four
  1D irreps; T1/T4 plus L10 T8), byte-for-byte of the entire
  sorted historical representative and norm dump
- All 28 complete native Hamiltonian cache pairs per compiler
  (L4/L6 T1/T4, L10 T1/T4/T8), full SHA256 of every target,
  coefficient and record order
- Native Ritz energy, direct residual and H-application count for
  L4/L6 T1 all four irreps (8 pairs per compiler)
- Negative opt-in missing prerequisite rejection
- All historical PR13–17 group/geometry/cache/baseline safety tests
- Both gfortran13 and gfortran14 CI matrix jobs must PASS

An independent seeded Python integer oracle is now checked into
`tests/cpc_cubic_axis_pruning_oracle.py`, performing **exhaustive**
comparison of all 41,664 distinct 3-site configurations at L=4,
plus 6,000 seeded random configurations each at L=6,8,10
against the complete original 27-candidate search.
The earlier standalone 45,000-case random check also passed.
This independent oracle is auxiliary evidence, not a substitute
for the complete native ED2 basis/cache/Ritz tests.

Only if numerical validation and a meaningful speed improvement both
hold should the PR advance to PRESTO same-source sort-pruning OFF/ON
benchmarks, independent audit and operator resource attestation.
Never merge to main without explicit user decision.
