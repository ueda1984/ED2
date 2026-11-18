# Canonical cubic \(O_h\) group certificate: experimental exact accelerator

**Status:** opt-in Draft code experiment. Numerical algorithms, accepted solver defaults and general group validator are unchanged. Do **not** claim publication-grade performance until a replicated, operator-attested PRESTO measurement.

## Empirical motivation

From PRESTO job \`19673.jpmjpr1911\` (N=1000, three magnons, Gamma/O_h, MNTE64, anchored cache), the four irreps independently passed 12/12 native Ritz numerical checks. The 8-thread total elapsed time is 5.38–5.47 s and the 8-thread basis time is ~5.20 s. We instrumented that basis stage separately in follow-up PRESTO output:

| Basis wall time (s) | T1 | T4 | T8 |
|---|---:|---:|---:|
| Preflight other | 0.000068 | 0.000069 | 0.000069 |
| Normal-form validation | 3.082175 | 3.087042 | 3.086126 |
| Orbit scan | 9.679079 | 3.141126 | 2.111575 |
| Finalization | 0.000160 | 0.000342 | 0.000275 |
| Total | 12.761482 | 6.228579 | 5.198046 |

The repeated *3.08 s* validation across thread counts demonstrates a serial bottleneck. Orbit scan itself scales 4.58× at 8 threads. The original validator allocates 48000×1001 32-bit site-image entries, ~183.3 MiB of temporary affine maps. The proposed scheme avoids that table.

## Exact reduction for the canonical cubic subgroup

Only when \`ED2_CUBIC_FAST_GROUP_VALIDATION=1\`:
- Require a **canonical even 4×4×4, 6×6×6, 8×8×8 or 10×10×10** spin-1/2 three-magnon lattice (\`NODmin=NODmax=3\`) at Gamma.
- Require normal-form periods \`(L,L,L,4,4,3)\`; \`M1=M2=M3=M6=0\`; \`M4,M5∈{0,2}\`; and spin/site occupancy constraints.
- Check the ACTUAL input six site permutations on every site against:
  \`Tx(x,y,z)=(x+1,y,z)\`, \`Ty=(x,y+1,z)\`, \`Tz=(x,y,z+1)\`, \`A=(-x,z,-y)\`, \`B=(y,-x,z)\`, and \`C=(z,x,y)\`, all modulo L.
- A,B,C are signed 3-axis permutations fixing the origin and have respective orders 4,4,3. Enumerate \`A^a B^b C^c\` for \`a=1..4,b=1..4,c=1..3\` in *the legacy ED2 order*; explicitly check all 48 maps are distinct and the last is identity.
- Signed 3-axis permutations have at most \`3! × 2^3 = 48\` elements, so 48 distinct valid products exhaust \`O_h\`. Multiplication for all \`48²\` pairs is checked exactly at the level of signed generators; the corresponding **runtime ED2 complex phases** must satisfy \`χ(PQ)=χ(P)χ(Q)\` within \`1e-10\`. All phases must be real ±1.
- The translations form \`(Z_L)^3\` and are checked on **every** site. Since point operations fix the origin, affine products \`T^tP\` are distinct for all \`L³×48\` normal forms. The point transforms preserve the translation lattice. At Gamma translations all have \`χ=1\`; thus the point representation extends to a full-group 1D character.

This is an exact finite-group certificate in the restricted geometry—not skipping the necessary conditions. If the explicit flag requests the certificate for an incompatible group, the code aborts without substituting results. If the flag is absent, the original \`validate_noncommuting_action\` runs unchanged, with its existing memory/time limits and fallbacks.

**Mathematical scope:** the certificate is specialized to the exact generators, coordinate numbering and one-dimensional real Gamma characters above. It is NOT a general replacement for noncommuting groups, nonzero momenta, other spin sectors, nor arbitrary crystallographic actions. The existing orbit enumeration, normalization, Hamiltonian cache and solver are identical after successful certification.

## Exact-basis reference checks

\`tests/cpc_cubic_fast_group_certificate_regression.py\` executes both validators against the exact **same input deck**, independent original anchored cache and same Lanczos settings, at L=4,6 for all four 1D Oh sectors and T=1,4, and at L=10 A1g T=1. A test-only \`ED2_CUBIC_BASIS_PROOF_DUMP=1\` writes every legacy \`list_s\` integer rank and 64-bit IEEE-754 bit pattern of \`list_r\` in sorted basis order. The file bytes must match **exactly**, not merely energies, between general and specialized validators. Energies, direct residuals, dimensions and, for T=1, Hamiltonian-action counts also must match.

The optional rank/norm dump modifies no computation, but is intended only for numerical regression: it writes a file in the per-trial working directory. It is **not** part of performance measurements. \`ED2_BASIS_DETAIL_TIMING=1\` is similarly optional and intended to measure the validation subphase. The preset PRESTO benchmark/audit records \`fast_group_validation\` per trial, checks the native certificate marker, and fails closed on translation-only sectors.

The PRESTO PBS script logs \`FAST_GROUP_VALIDATION\` in its allocation-preflight record, validates it is 0 or 1, and unsets it only inside the generic \`make check\` build subshell (the 100-site anisotropic chain regression is deliberately incompatible with cubic-only settings). The actual 1000-site benchmark retains the parent's explicit flag.

## PRESTO safety and next check

Create a *new immutable source stage*, pinned to the reviewed Draft PR head. Verify the \`ed2-source.sha256\` before extracting, and do **not** reuse prior PR #10/#11 source stages. Only after the operator's independent \`place=exclhost\` attestation, run a small one-sector N=1000 A1g 1/4/8-thread pilot. The PBS \`qsub -v\` argument must list **variable names only**, not inline values containing commas. Enable:
- \`ED2_CUBIC_ANCHORED_REP=1\`;
- \`ED2_CUBIC_FAST_GROUP_VALIDATION=1\`;
- \`ED2_BASIS_DETAIL_TIMING=1\` for a phase-specific first pilot (turn off for performance-quality measurements);
- \`ED2_MNTE=64\`.

Audit exact state dimension, E, Ritz residual, H applications, bit-for-bit input SHA, all timings and peak RSS against the old PRESTO job \`19673.jpmjpr1911\`. A single non-replicated pilot does NOT prove reproducible speedup or operator-independent exclusivity.

If validation falls below 0.1 s, A1g at 8 threads would have about \`2.11s\` basis scan plus small extra setup, **not a guaranteed wall-clock target**. A remaining roughly 1-second non-scaling component is inferred from the orbit scan scaling and must be instrumented separately before more code changes.
