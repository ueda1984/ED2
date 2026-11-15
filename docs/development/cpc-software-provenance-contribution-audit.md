# CPC contribution audit: published QS³, upstream ED2 and the validated fork

**Purpose:** state exactly what would be new to a CPC reader, without
renaming an existing solver algorithm, misattributing upstream work,
confusing maintenance with scientific model coverage, or claiming that
the test-only SciPy programs form part of the ED2 production executable.

**Snapshot:** 15 November 2026. The public upstream
[QS-Cube/ED2](https://github.com/QS-Cube/ED2) and the fork's
[main](https://github.com/ueda1984/ED2/tree/main) have the same
observed commit `49941dc4b437a8fb6c2b5c8dac7d621866731d7c`.
The tested development line is
[`codex/ed2-cpc-stationary-dm-scaling`](https://github.com/ueda1984/ED2/tree/codex/ed2-cpc-stationary-dm-scaling).
Always replace the branch with an immutable full SHA for an actual archive release.

## The correct published-program baseline (first ED2 program paper)

**Editorial correction:** ED2 itself has **not yet
been reported in a peer-reviewed CPC software article**.
The already-public
[QS-Cube/ED2](https://github.com/QS-Cube/ED2)
repository must not be conflated with
an earlier *published ED2* paper.
The previous published CPC software is
**QS³** [CPC **277**, 108369 (2022)],
and its [original Fortran repository is
QS-Cube/ED](https://github.com/QS-Cube/ED).

The primary *research contribution of a first
ED2 program article* is the new computational
model coverage relative to **this original
published QS³/ED code**, not merely
the Phase 1 code-hygiene changes made
to ED2 after it had appeared on GitHub.

The [old-to-new source and feature matrix]
(cpc-qs3-ed-to-ed2-scope-matrix.md)
directly checks the two repositories.
A [six-level same-physics test]
(cpc-qs3-ed-ed2-shared-xxz-spectrum-ci.csv)
also builds the original and modern
executables and compares the identical
fixed-magnetization XXZ sector with
an independent 6x6 Hamiltonian:
the largest cross-program difference
is \(2.89\times10^{-15}\).
This makes the progression from
published QS³ to generalized ED2 a
**checked mathematical and source-level
transition**, not just a narrative claim.

The old QS³ implementation also provides
a dynamical spin-structure-factor
pathway. The current ED2 executable
has not been shown to include a
matching generalized dynamic path.
Do **not** market ED2 as a strict
backward-compatible superset of all
QS³ outputs. Separate absent dynamic
capabilities from genuine improvements
in mixed-spin, nonconserving
Hamiltonian representations.

## 1. Compare against two *different* prior baselines

The [published 2022 QS³ article](https://doi.org/10.1016/j.cpc.2022.108369)
describes a near-saturation spin-1/2-focused quantum-spin solver.
ED2 broadens the **scientifically accessible model domain** to
non-identical spin magnitudes, anisotropic XYZ exchange, directed DM
and symmetric off-diagonal exchange, together with site-dependent
magnetic fields and nonconserving Hamiltonians. This broadening is
the paper's potentially significant computational-physics scope.

But the [upstream ED2 repository](https://github.com/QS-Cube/ED2)
**already contains** these generalized model inputs, its
Fortran eigensolver kernels, the existing `wk_dim` and `MNTE`
cache parameters, and the output facilities. **Do not claim**
that our development branch newly invented the generalized
Hamiltonian specification, Krylov iteration, OpenMP or cache controls.

The contribution of the *current fork development* is instead
numerical-interface repair, reproducibility and reliability
contracts, portable source distribution, independent
physical reference testing, and documented performance
properties for a concrete mixed-spin nonconserving model.

The paper must make clear whether it describes the
generalized ED2 **as a program**, including upstream work,
or claims specific source-level engineering advances on top
of the already-public upstream ED2. Both are possible,
but neither justifies citing the fork as a wholly
independent new program.

## 2. Actual checked production source delta against fork main

We compared the **GitHub content hashes for the entire `source/`
directory** between fork `main` and the development branch.
The numerical code is **not unchanged across that full comparison**.

| Path | Changed in development fork? | Grounded purpose |
| --- | :---: | --- |
| `source/blas_lapack_interfaces.f90` | Added | Explicit LP64 BLAS/LAPACK signatures; compiler calling contract |
| `source/ham2vec.f90` | Yes | Site/field record-index correction in cached, dense and on-demand Hamiltonian paths |
| `source/lanczos.f90` | Yes | Previously diagnosed early-check and zero-Hamiltonian Lanczos behavior |
| `source/krylov_solvers.f90` | Added | Isolated/reorganized Krylov implementation with numerical checks |
| `source/krylov_observations.f90` | Added | External-facing numerical diagnostic and reproducibility facilities |
| `source/krylov_residual_policy.f90` | Added | Optional direct-eigenpair residual contract |
| `source/eigen_solver.f90` | Yes | Solver call paths and checked linear-algebra interfaces |
| `source/Makefile.am` | Yes | Explicit build/distribution of new and reorganized Fortran modules |

No files were removed from the listed production `source/`
directory in the comparison. Other unchanged production
paths inherit the upstream code.

The key **time distinction** is that these Fortran changes
were introduced by the *earlier Phase 1 maintenance*
work. The recent stationary DM, nearly-degenerate
reference, Q-resolvent, cached-row and publication
documentation additions did **not** independently add
new numerical Fortran kernels. Therefore both
statements can be true: the recent scientific test phase
leaves `source/` unchanged, and the full development
branch differs from upstream in eight `source/` entries.

Selected evidence and failure-before/fix-after documentation:

- [BLAS/LAPACK typed-interface and INFO-path audit](phase1-lapack-interfaces.md)
- [Sparse magnetic-field site/record-index bug reproduction and repair](phase1b-sparse-field-indexing.md)
- [Lanczos early-exit/breakdown repair](phase1c-lanczos-stability.md)
- [Krylov stress and reproducibility](phase1d-krylov-stress-benchmarks.md)
- [Precision/convergence contract](phase1e-precision-contract.md)
- [Diagnostic/reproducibility and opt-in residual tests](phase1f-diagnostics-reproducibility.md), [phase1g](phase1g-optin-residual-policy.md)

This is a path/content comparison, not an independent
license, authorship or API-stability certification.

## 3. Claims that the present evidence *does* establish

| Manuscript claim | Evidence | Qualified limit |
| --- | --- | --- |
| Mixed-spin, fully non-U(1) XYZ/Gamma/DM input and output in the executable | [stationary DM response 60/60](cpc-stationary-dm-application.md) | Already present in upstream ED2, not a newly invented solver |
| Oriented DM sign and unusual xy/yx output order are handled correctly | Independent correlations, energy derivative, exact small-size spectrum | Known Hellmann–Feynman identity, not new physics |
| A published spin-1/2 XYZ-DM Hamiltonian is represented with correct spin conventions | [12 independent literary-model checks](cpc-literature-matched-benchmark.md) | Finite-size Hamiltonian comparison, not complete reproduction of the article's phase diagram |
| Numerical acceptance is more than a solver status flag | Native wavefunctions, independent dense/sparse Hamiltonians, projector trace and direct residual | Global lowest-state inclusion is not certified by residual alone |
| Upstream `wk_dim` and `MNTE` resource controls have analytically quantifiable consequences | [72 MNTE cases](cpc-mnte-capacity-evaluation.md), [108 paired trials](cpc-cache-pareto-evaluation.md) | Exact row bound \(7L/2-1\) is **model-specific**, not a safe default for arbitrary Hamiltonians |
| A small worked physical example needs no Python numerical stack | [portable 4/6-site three-D quickstart](cpc-portable-dm-quickstart.md) | Reference answer is frozen from independent SciPy/dense tests; script is an example, not library ABI |
| Repeated CPU-pinned measurements can be reproduced with source/input hashes | [protocol](cpc-controlled-node-protocol.md) and [postprocessor](cpc-controlled-node-protocol.md#d-relevance-to-the-cpc-paper) | Only shared-runner functional checks are complete; no independent dedicated-node performance evidence yet |

The significance argument should be **usable physical
coverage + correctness + distribution**, not a claim of
unprecedented asymptotic speed, universal precision,
new magnetic phases or automatic symmetry certification.

## 4. Position relative to other available ED tools

General exact-diagonalization programs such as
[QuSpin](https://quspin.github.io/QuSpin/) can already
represent many spin Hamiltonians and custom symmetry
representations. ED2's relevant user-facing value is
a compact Fortran executable with a documented
file-based interface and verified mixed-local-spin
operators, not proof that competing packages are
incapable of solving equivalent finite-dimensional
physics. Claims about relative speed, memory or
API simplicity would require **matched, controlled,
explicitly defined cross-program benchmarks**,
which have not been performed.

## 5. Release and publication blockers to resolve with maintainers

1. **Version conflict:** `configure.ac` declares `QS3-ED2 1.0.0`,
   while `CITATION.cff` declares `0.1.0` and an unverified
   release date. These must not be silently rewritten;
   an official release/tag must be chosen.
2. **Authorship/credit:** `README.md` and `CITATION.cff`
   currently differ in the contributors displayed;
   2022 QS³ authors and public ED2 maintainers must
   be credited correctly, with the actual CPC author
   list confirmed by the responsible team.
3. **Code archive:** No Mendeley Data CPC program
   accession or software DOI is assigned; don't
   invent one. Prepare a buildable source distribution
   at one immutable source revision.
4. **Performance:** obtain a repeatable, isolated
   and CPU-pinned run with at least several timed
   repetitions; report medians, MAD/uncertainty,
   source/binary hashes, compiler/BLAS, CPU topology,
   RSS and Hamiltonian applications. Existing
   GitHub-hosted wall time is exploratory only.
5. **Program interface:** the production deliverable is
   currently an executable, **not** an installed
   reentrant Fortran library or stable ABI. A CPiP
   article can describe an executable program;
   do not publish an unimplemented API.
6. **Publication track:** confirm CPiP versus CP and
   requirements at submission. CPC's [official
   journal scope](https://shop.elsevier.com/journals/computer-physics-communications/0010-4655)
   calls for a significant open-source program,
   a substantive contemporary physics context and
   a clear explanation of novelty and importance.
   Evaluation by editors is not guaranteed.

## 6. Proposed defensible thesis (for author review)

> ED2 extends the near-saturation QS³ computational framework
> to input-defined quantum-spin Hamiltonians with mixed local
> spin magnitudes and symmetry-breaking anisotropic interactions.
> We document a tested source-distributed executable,
> independent Hamiltonian/observable validation, explicit
> numerical failure and accuracy diagnostics, and model-aware
> cache resource costs. The reproducer evaluates a staggered
> DM response in a nonconserving mixed-spin chain, rather
> than claiming a new eigensolver or a previously unknown
> chiral phase.

Before submission, verify that the final public program
archive and the authorship accurately cover the generalized
capabilities and code contributions being described.