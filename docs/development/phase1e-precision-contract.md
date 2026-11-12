# Phase 1-E: proposed Lanczos precision and convergence contract

Date: 2026-11-12 (Asia/Tokyo). Starting revision:
`b56fe21678e62b8b6216fe0ec742d5dbd62a8df6`.
Branch: `codex/ed2-phase1e-precision-contract`.

This is a specification proposal and validation phase. It implements **no new
solver policy**, changes no numerical production source, and changes neither
the default `lnc_ene_conv` nor its meaning. Read the [Phase 1-D proposals](phase1d-policy-proposals.md),
[stress report](phase1d-krylov-stress-benchmarks.md), and [Phase 1-C policy](phase1c-lanczos-stability.md)
for the existing behavior. The [Phase 1-E experiments](phase1e-validation.md)
distinguish results actually obtained from requirements for a later implementation.

The proposed new option/status names below are conceptual fields, **not accepted
CLI keywords or a stable API**. They require review before implementation. In
particular, this document does not choose new numerical defaults on the user's behalf.

## 1. Scientific questions must have separate answers

A successful solve should eventually return independently assessable statements:

| Question | Evidence/result | What it cannot establish alone |
| --- | --- | --- |
| Did iteration stop under the requested policy? | termination reason, budgets, energy checks, residual checks | Accuracy of a requested observable |
| Is the returned pair an approximate eigenpair? | direct residual for a normalized vector, finite arithmetic, norm error | That it is the lowest pair in the physical sector |
| Was the requested lowest space included? | independent spectrum/enclosures or explicitly qualified search evidence | A global proof from a random start or zero residual |
| Is a state or cluster resolved? | residual, orthogonality, boundary separation, dimension/rank | Unique vectors inside an exactly degenerate eigenspace |
| Is an observable accurate? | operator norm and applicable state/subspace bound, or independent reference | Accuracy inferred from energy changes or variance cancellation |

Do not collapse these into one `converged` boolean. Status zero in the current
kernel means its **current** exit/certification rules were satisfied. It is not
the proposed full result contract and does not establish lowest-state inclusion.
An exact excited start is a required, legitimate counterexample.

## 2. Energy, residual and scale definitions

For normalized `q`, retain the returned Ritz value `lambda`, independently
computed Rayleigh value `rho = real(q^H H q)`, norm defect, and
`r = ||H q - lambda q||_2`. Do not silently replace the reported Ritz value by
`rho`; report their difference. Use a fresh operator application for final
validation, including after ALG=1 reconstruction and ALG=2 restart. A projected
recurrence estimate is useful telemetry, not a substitute for that validation.

For Hermitian H in exact arithmetic, `min_j |lambda-E_j| <= r`; this locates
**some** eigenvalue. A small residual does not identify its index. A normalized
Rayleigh quotient obeys the ground-state variational bound, but that alone is
not an error bound on the requested ground energy.

Proposed pair criterion:

```
r <= residual_atol + residual_rtol * S
S = max(operator_scale, abs(lambda-c))
```

`residual_atol` has energy units; `residual_rtol` is dimensionless. Require finite,
nonnegative values with at least one positive control. `c` is a recorded energy
origin (zero by default for compatibility); use the same origin for the operator
scale. Record whether that scale is a caller-supplied bound, a certified bound,
or an observed estimate. A maximum observed `||Hq||` is not `||H||` and must not
be described as such. For H=0 and S=0, the absolute term and explicitly defined
zero-residual behavior apply; do not divide by zero or insert an arbitrary unit
energy into a relative test. Stable norm computation and finite-value checks
remain necessary across scales.

Centering by `c` can avoid letting a large identity shift dominate a relative
tolerance. It does not recover bits lost when forming/applying `H-cI`. An
unvalidated midpoint or estimated spectral width is not a certified operator
bound. This phase uses the **known reference** half spectral width only for an
offline comparison; it supplies no production scale estimator.

Keep `lnc_ene_conv` as the legacy energy-change control. Current behavior is
`abs(Eold/Enew-1)` away from zero and an observed-scale absolute change near zero;
ALG=2 requires two successive successful energy checks on retained `NOK` values.
Current exact/numerical closure and full-space exits have separate conditions.
Preserve those semantics until a numerical-policy change is authorized. Future
energy telemetry should include absolute change, relative/scale-normalized change,
check interval and consecutive-check count. Small change is a stagnation or
progress diagnostic, not a universal energy/vector-accuracy certificate.

## 3. State and invariant-subspace accuracy

For a known simple target eigenvector `u0`, let `delta` be a **positive lower
bound** on the distance from `lambda` to every unwanted eigenvalue. Then

```
sin(angle(q,u0)) <= r / delta
|<q,Aq> - <u0,A u0>| <= 2 ||A||_2 sin(angle(q,u0))
```

The first relation follows by expansion in H's orthonormal eigenbasis. Without
target association and spectral-separation evidence, it is not a certificate of
the ground state. A neighboring Ritz-value gap need not bound the actual gap:
an omitted lower state or unresolved cluster can make that interpretation wrong.
Label such a gap and resulting accuracy estimate **heuristic**.

For an orthonormal block Q with diagonal Ritz matrix Lambda, record the block
residual `R=HQ-Q Lambda`, individual residuals and `||Q^H Q-I||`. If a complete
target invariant space of the same dimension is isolated from the unwanted
spectrum by the corresponding `delta`, a safe reference assessment is
`||sin Theta||_F <= ||R||_F/delta`; it also bounds spectral-norm leakage.
Do not assume a general spectral-norm Sylvester bound with constant one without
the additional separation hypotheses. Check the projector distance when the
reference space and requested space have equal dimension.

For exact degeneracy, only the full eigenspace is unique. A partial selection may
be a valid orthonormal subset, but its projector is not uniquely prescribed.
For a near cluster, small residual compared with the **external** gap can resolve
cluster membership while leaving individual states unresolved at internal gaps.
Return `unresolved_cluster` when the requested cutoff bisects an unresolved
cluster. Offer a separately controlled enlargement of the wanted block or an
honest partial-space result; never silently increase memory or promise a unique
ground vector. Exact degeneracy also makes arbitrary pure-state magnetization
nonunique. A cluster-average observable such as `Tr(P A)/rank(P)` has a different,
gauge-invariant meaning and must be explicitly requested.

For the spin-1/2 tests, `||S_i^a||=1/2` and, for distinct sites,
`||S_i^a S_j^b||=1/4`. A sufficient uniform absolute magnetization target `1e-9`
thus needs a state-angle bound at most `1e-9`, not merely `1e-5` or a small energy
change. Given suitable separation, that implies `r <= delta*1e-9`. This is a
sufficient conservative condition, not a necessary condition or a new default.
The existing elementwise comparisons remain `atol=rtol=1e-9`; they are not relaxed.

## 4. Lowest-state inclusion and initialization

Describe inclusion separately as `not_assessed`, `heuristic_search_evidence`, or
`certified_by_external_evidence`, with method/provenance. Residual success,
happy breakdown, full *visited* Krylov dimension, repeat agreement and an initial
nonzero ground component do not by themselves prove global inclusion.

For future contexts, expose explicit vectors/blocks and owned RNG state. Record
generator/version, every thread's state, scheduling/thread count and initial
block rank. A scalar start samples at most one direction per exact eigenspace;
multiple generic starts improve coverage but supply no unconditional proof.
Check initial/retained block rank, independent starts and disagreement, without
silently rerunning until a favored answer appears. A deterministic full-sector
reference or validated spectral enclosure is stronger external evidence; imposing
it on every large solve would defeat the purpose of the sparse method.

## 5. Roundoff, budgets and failure semantics

Binary64 cannot promise relative energy changes below its effective arithmetic
resolution. Identity shifts, near-zero values, cancellation in H applications and
operator/model error affect attainability. Record requested tolerances and any
estimated arithmetic/operator uncertainty separately. Do not silently raise a
requested tolerance to `epsilon` or `sqrt(epsilon)` and call that request satisfied.

Proposed termination reasons include `requested_policy_satisfied`,
`exact_subspace_closed`, `numerical_subspace_closed`, `budget_exhausted`,
`roundoff_limited`, `invalid_input` and `nonfinite_or_failed_action`. Keep these
orthogonal to pair/inclusion/cluster assessments. A future roundoff-limited
assessment requires stagnation and residual evidence plus an explicit uncertainty
estimate; it does not turn budget exhaustion into convergence. A small residual
may make a failed candidate scientifically interesting without authorizing its
publication as a successful current CLI solve.

Retain current CLI nonzero failure exits and result-file suppression. A future
library may return last candidates in a clearly separate diagnostic object with
validity flags, never indistinguishably in its successful-result fields. The
Phase 1-D test driver already exports failed ALG=2 candidates under distinct
filenames for **tests only**; ALG=1 has no equivalent failed-vector export.

## 6. Proposed implementation sequence and approval gates

1. **Telemetry and ownership first.** Add additive internal result fields for
   tolerances, scales/origins, Rayleigh/Ritz values, direct residuals, norm/rank,
   H applications, budgets and all independent assessments. Introduce explicit
   RNG/start ownership in a separate change. Preserve old CLI defaults/results.
   Test deterministic replay at fixed thread count, cleanup and failure outputs.
2. **Opt-in residual policy.** After review of option semantics, add independent
   absolute/relative controls with explicit zero/shift/scale behavior. Compare
   with legacy stopping on identical starts and achieved accuracy. Do not infer
   online stop counts or speedups from this phase's endpoint reclassification.
3. **Cluster and observable assessment.** Add rank/closure diagnostics, uncertainty
   labels and caller-supplied or independently justified separation information.
   Test wanted cuts, exact multiplicity, small internal/large external gaps, and
   observable error bounds. Decide explicitly whether a request needs individual
   vectors, full-cluster projectors, or cluster-average observables.
4. **Numerical methods only after a separate decision.** Evaluate ALG=1 selective
   or full reorthogonalization versus the two-pass memory contract; evaluate
   ALG=2 restart closure/retention and redundant H-image/workspace work separately.
   Measure memory, H calls and time at matched achieved accuracy. No such method
   or optimization is implemented or selected here.
5. **Defaults and API extraction last.** Require approved tolerances/status
   semantics, unchanged physical regressions, excluded/weak starts, zero/shifted/
   scaled H, complete/partial clusters, budget/roundoff tests, and fixed-state
   one/two-thread checks. Obtain another compiler family/backend before broad
   portability claims. Preserve historical failed findings and require explicit
   review before changing a default or claiming CPC submission readiness.

Open decisions: default residual controls and scale estimator; whether legacy
energy checks remain mandatory under a new opt-in policy; certified versus
estimated gap inputs; cluster-enlargement budget limits; last-candidate ownership;
and RNG compatibility across compilers. These remain proposals for the next
authorized numerical phase, not unfinished production edits in Phase 1-E.
