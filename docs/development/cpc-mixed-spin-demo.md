# CPC science demonstrator: mixed spins and nonconserving interactions

**Status: independent small-system validation completed.**
This is a testable physics example, not a finite-size extrapolation,
phase transition, new physical discovery or large-scale performance claim.

## Model and conventions

An open alternating spin chain with
`S = [1/2, 1, 1/2, 1, ...]`,
inequivalent `Jx, Jy, Jz`, oriented `Dz`, symmetric `Gamma_xy`,
site-dependent `hz` and an `hx` control parameter. The complete
coupling functions are recorded in
[`tests/mixed_spin_nonu1_demo.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-mixed-spin-demo/tests/mixed_spin_nonu1_demo.py).
Both exchange anisotropy and Gamma break total `Sz` conservation;
the transverse field further breaks it.

The independently constructed Hamiltonian uses physical spin matrices
with dimensions `2S_i+1` and the **same declared source conventions**:
`+ h_i · S_i`; antisymmetric oriented DM term
`D_ij · (S_i × S_j)`; and symmetric Gamma cross terms.
All correlation columns are compared in ED2's actual
`xx,yx,zx,xy,yy,zy,xz,yz,zz` order.

Two spaces are compared for `N=4`:

- Full product Hilbert space: `NODmax=6`, dimension 36.
- Projected nonconserving space: `NODmax=2`, dimension 13.
  This is `P H P`, **not** an exact symmetry sector.

For `N=6`, the full space is `NODmax=9`, dimension 216.
No user-supplied lattice permutation symmetry is imposed.

## Independently measured results

**Source commit:** `86768728b3191c20e083062762832a9dafd80107`.
**CI validation:** [GitHub Actions run](https://github.com/ueda1984/ED2/actions/runs/43393654799),
GNU Fortran 13, Ubuntu 24.04, OpenBLAS, one OpenMP thread.
Six cases, all passed; no reference data modified to fit output.

| Sites | `hx` control | `NODmax` | Dimension | `E0` | Max spectral error | Max local-mag error | Max Cartesian CF error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 0.15 | 6 | 36 | -2.2184277260 | 5.77e-15 | 2.54e-15 | 2.44e-15 |
| 4 | 0.38 | 6 | 36 | -2.2900665380 | 4.44e-15 | 3.55e-15 | 3.44e-15 |
| 6 | 0.15 | 9 | 216 | -3.5693260311 | 3.11e-15 | 8.38e-15 | 5.55e-15 |
| 6 | 0.38 | 9 | 216 | -3.6640192356 | 6.22e-15 | 6.27e-15 | 4.72e-15 |
| 4 | 0.15 | 2 | 13 | -2.0192793301 | 2.22e-15 | 1.44e-15 | 9.99e-16 |
| 4 | 0.38 | 2 | 13 | -2.0351796775 | 3.11e-15 | 7.77e-16 | 6.66e-16 |

The errors compare ED2 output against a **separately constructed,
NumPy diagonalized product-space matrix**, not against another ED2
solver or the same implementation's matrix storage. The maximum
spectral error covers **all** eigenvalues in the represented space.

For four sites at transverse-field control 0.15,
restricting to `NODmax=2` increases the energy by
`0.1991483959`.
Site 1's `<Sz>` changes from
`0.232079` (full) to
`-0.316994` (projected);
site 2's `<Sz>` changes from
`-0.608049` to
`0.753866`.
These are qualitative differences caused by the approximation, **not
disagreement between ED2 and independent diagonalization**.

At transverse-field control 0.38, the corresponding projected-space
energy difference is `0.2548868605`.

The full machine-readable records, including gaps and local `<Sz>`,
are in [cpc-mixed-spin-demo-data.json](cpc-mixed-spin-demo-data.json).

## How to reproduce

From a source archive containing this test script (test-only NumPy is
required), build `source/QS3ED2`, then run:

```sh
python3 tests/mixed_spin_nonu1_demo.py \
  --exe /absolute/path/to/QS3ED2 \
  --output /new/empty/parent/demo-results
```

The parent output path must **not** already exist. The script writes
the exact NAMELIST input, site spins, local caps, oriented bonds,
fields, pair lists, STDOUT, ED2 outputs and a JSON summary under it.

The script deliberately compares all 36 or 216 dense eigenvalues,
the ground-state local magnetizations and all Cartesian pair
correlations. This is **not** evidence that large-scale iterative
`ALG=1/2` ground-state inclusion is certified.

## Next physics milestones

1. Verify **iterative ALG=1 and ALG=2** on the same nonconserving
   mixed-spin Hamiltonian with independently checked residuals and
   energies (do not assume a small residual guarantees the ground state).
2. Add a physically motivated observable or parameter sweep at
   larger size, with finite-size and convergence diagnostics.
3. Benchmark controlled caching and OpenMP settings, recording wall
   time, peak RSS and Hamiltonian applications, **without** making
   unmeasured performance comparisons with QS³.
4. Use this validated code path as the basis for the CPC physics
   application and numerical methods sections. The published
   QS³ code already supported Lanczos, correlation functions
   and OpenMP, so those alone are not novelty claims.
