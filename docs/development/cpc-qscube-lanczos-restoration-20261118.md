# Original QS-Cube two-stage Lanczos restoration: provenance and numerical validation

**Status:** development branch / opt-in implementation, 2026-11-18. This is an engineering and numerical-robustness change, **not the invention of a new Lanczos algorithm**.

## Historical evidence and what was actually changed

The original **published QS³** [QS-Cube/ED `source/lanczos.f90`](https://github.com/QS-Cube/ED/blob/main/source/lanczos.f90) already used the energy-first/replay-once pattern; the subsequent [QS-Cube/ED2 `source/lanczos.f90`](https://github.com/QS-Cube/ED2/blob/main/source/lanczos.f90) retained that design. Both sources already implemented the following *conventional low-memory two-stage Lanczos workflow* in `lancz_eigen_val`:

1. In the first loop, construct the tridiagonal Lanczos coefficients `alpha(i)`, `beta(i)` and monitor the lowest `DSTEVX` eigenvalue (`conv`, `lnc_ene_conv`). No full-space Ritz vector is reconstructed at each convergence check.
2. After the loop, restart from the saved initial vector `psi_io`, replay the Lanczos recurrence once, and accumulate its full-space Ritz vector with the tridiagonal eigenvector coefficients `z(i,1)`.
3. The caller computes `<GS|H|GS>` and `<GS|H^2|GS>` for a quality indicator. The original version **does not use a directly measured residual as its acceptance criterion** and **does not implement the inverse-iteration refinement** associated with some traditional TITPACK workflows.

The later [ED2 development source, pinned at 673a7faa](https://github.com/ueda1984/ED2/blob/673a7faa075b6f8a7839fb00bbcf193ee920790d/source/krylov_solvers.f90) introduced a safety-focused numerical policy: after apparent energy convergence it reconstructed a full-space eigenvector and independently computed `||H psi - lambda psi||`. If that check failed, it continued iterating and re-created the full-space vector at subsequent energy-stagnation checks. The stronger acceptance condition was useful, but reconstructing the vector many times could dominate runtime.

The opt-in `ED2_LANCZOS_TWO_STAGE=1` branch **restores the original economical replay pattern** while adding properties not present in the original QS-Cube source:

- A cheap tridiagonal projected residual proxy `abs(beta_k*z_k)` as an **intermediate gate**, never treated as the final numerical certificate.
- A final **direct** Ritz residual `||H psi-lambda psi||`, together with the optional pre-existing absolute/relative residual policy and the pre-existing energy-change gate.
- Controlled retry after a failed direct check, a conservative breakdown/iteration-budget exit and finite-value checking.
- Centered energy variance `||H psi-<H>psi||^2` measured without subtracting two large floating-point numbers.
- An explicit opt-in switch, matched-binary experiment and independent diagnostics and regression verification.

Preserve the published QS³ attribution (CPC 277, 108369 (2022), DOI 10.1016/j.cpc.2022.108369), and existing QS-Cube/ED and QS-Cube/ED2 contributor/license provenance. The source history is evidence of algorithmic continuity, not evidence of authorship of the present experimental patch. Avoid describing the two-stage algorithm as novel, invented, or newly introduced by this project. The narrow defensible claim is a **restoration/reconciliation of a historically existing efficient workflow with modern residual-based accuracy checks**.

## Controlled PRESTO demonstration (job 18730.jpmjpr1911)

1000-site periodic cubic `S=1/2` Heisenberg model with three magnons, translation-only Gamma sector dimension `D=166,167`; `standard`, `node03`, PBS `place=exclhost`, same executable hash and fixed input, one warmup and three measurements per mode/thread count. Both modes used the same compiler, real-cache setting and accuracy policy. Reported values are medians across the three measured repeats.

| OpenMP threads | Historical-style repeated-residual checking, s | Restored two-stage + direct residual, s | Observed wall-time ratio |
|---:|---:|---:|---:|
| 1 | 197.0266 | 35.0962 | 5.61x |
| 16 | 49.2224 | 4.5382 | 10.85x |
| 32 | 50.1983 | 3.8098 | 13.18x |

Source: preserved PRESTO job 18730 archive and its independent review `ED2_PRESTO_twostage_18730_verified_report.md` (separate evidence artifact, not part of the source repository). In the solver, the ALG=1 Hamiltonian-application count fell from `23,170` to `990` (excluding one post-solver quality `H` application in each case). Across modes, `495` Lanczos steps, `E=736.7224213680653`, and direct Ritz residual `5.107622778767303e-9` matched. The centered energy variance was `2.6087810890994963e-17` for the two-stage run.

These numbers demonstrate an implementation/runtime improvement **relative to the current repeated-reconstruction path**. They are **not** evidence of a new eigensolver algorithm or of faster runtime than the original QS-Cube software under matched hardware and numerical controls.

The independent audit checked 24/24 records, 24/24 diagnostic hashes and nine physically present deduplicated input files. `publication_ready=false`: operator-side exclusivity attestation, hardware stability and stronger statistical evidence remain to be established. A small direct residual certifies an approximate eigenpair, **not inclusion of the global lowest state**.

## Integrated regression and acceptance boundary

The regression script `tests/cpc_qscube_lanczos_regression.py` exercises the same compiled kernel with the mode switch off/on, explicit starts, a separate NumPy Hermitian eigenspectrum reference, complex noncommuting matrices, degenerate/clustered energies, zero-energy states, a one-dimensional space, numerical breakdown, a ground-excluding start, an insufficient iteration budget, satisfiable and unattainable strict residual policies, and invalid switch values. Optional production-model testing uses the periodic 216-site/three-magnon Heisenberg sector and verifies a reduction in H applications at matched physical parameters.

The development default remains the existing path until a focused code review, compiler-independent regression, and provenance review of this change are complete. The separate inverse-iteration (e.g. MINRES with appropriate shifts/preconditioning) is **not** part of this change and should be assessed only on cases where direct residual accuracy proves inadequate.

## Attribution and CPC scope

The central CPC submission must distinguish the historical QS-Cube algorithms and already-public physical-model capabilities from modern software engineering: broader model handling (with provenance), robust numerical acceptance and diagnostic contracts, reproducible tests, source distribution and measured performance of implementation repairs. Do not claim that 5.61-13.18x results reflect an invented Lanczos algorithm; present them as the impact of removing repeated full-space reconstruction while preserving the stronger direct-residual criterion.
