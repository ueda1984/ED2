# Theory and numerical algorithms

This page describes features implemented in the source, rather than
promising support for every possible symmetry or representation.

## Basis construction and truncation

ED2 constructs states from the fully polarized reference state by
applying spin-lowering operators. A site may have its own spin `S_r`
and local occupation cap. Without restrictions, the dimension is

$$
\dim\mathcal H_{\mathrm{full}}=\prod_{r=1}^{N}(2S_r+1).
$$

`NODmin` and `NODmax` restrict the total number of lowering
excitations; `FILE_NODmax` additionally restricts each site.

!!! important "Truncation is not always a symmetry sector"
    If transverse fields or other terms change excitation number,
    a restricted NOD calculation is a **projected Hamiltonian** on
    that space, not generally a conserved block of the full model.

## Hamiltonian convention

The supported interactions include XYZ exchange, DM exchange,
symmetric off-diagonal Gamma terms, and site-dependent fields.
For the field terms, the implementation uses

$$
H_{\mathrm{field}} =
+\sum_r (h_r^x S_r^x+h_r^y S_r^y+h_r^z S_r^z).
$$

To model a physical Zeeman term `-B · S`, input the negative
coefficients `h = -B`. DM signs follow directed bonds; the
parser reverses the DM coefficients if it swaps the two site indices.
See [Input format](../user-guide/input-format.md).

## Symmetry reduction

ED2 accepts up to six user-defined site permutations and their
orders and one-dimensional character labels. These may represent
translations and **compatible** point-group operations.

It is the user's responsibility to ensure that the supplied operations
are genuine symmetries of the Hamiltonian and mutually compatible.
The input parser does not generally validate noncommuting group
relations or construct arbitrary point-group irreducible
representations.

A lattice reflection/inversion permutation is **not**
the global *spin-inversion* map `S_i^z -> -S_i^z`.
A separate general spin-inversion sector implementation has not
been established in the current source.

## Hamiltonian-vector product

ED2 applies the Hamiltonian using a transition cache and
on-demand matrix elements. `wk_dim` controls the extent of caching.
The represented Hamiltonian must not depend on the cache policy.
Independent small-system and cache-on/off tests are included in
the development validation suite.

## Eigensolvers

| `ALG` | Method | Use |
| --- | --- | --- |
| 1 | Three-term Lanczos, with candidate-vector reconstruction | Single low-lying eigenpair |
| 2 | Restarted Hermitian Krylov solver, multiple Ritz pairs | Several low-lying eigenpairs |
| 3 | Dense Hermitian diagonalization via LAPACK | Small sectors |

ALG=1/2 measure direct `||H psi - lambda psi||` residuals.
In the Phase 1-G development branch, an optional user-supplied
residual acceptance condition is added to their legacy
energy-convergence rules. See
[Solver accuracy and reproducibility](../user-guide/solver-accuracy.md).

## Numerical interpretation

A small residual for a normalized vector means that it nearly solves
the eigenvalue equation for *some* eigenvalue of the represented
Hermitian operator. It does **not** establish that the true lowest
eigenspace was found. If levels are nearly degenerate, individual
eigenvector components and physical expectation values may be much
less accurate than the energies or residuals suggest.

Validate observables and ground-space inclusion separately whenever
these distinctions affect the scientific result. The
[Phase 1-G report](../development/phase1g-optin-residual-policy.md)
records reproducible counterexamples and the limitations of the
optional tolerance policy.
