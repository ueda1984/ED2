# Theory and algorithms

ED2 constructs a quantum-spin Hamiltonian in a selected basis and solves its
Hermitian eigenvalue problem. The result is exact up to numerical error **for
the retained space**. A truncation that is not an invariant subspace of the
physical Hamiltonian is an approximation to the full problem.

## Local spins and occupation basis

Site `i` has spin `S_i` and magnetic quantum numbers
`m_i=-S_i,-S_i+1,...,S_i`. Different sites may have different spins. Define
`n_i=S_i-m_i`, measured from the fully polarized state. Before restrictions,

$$
\dim\mathcal H_{\mathrm{full}}=\prod_i(2S_i+1).
$$

ED2 enumerates states satisfying the local caps in `FILE_NODmax` and
`NODmin <= sum_i n_i <= NODmax`. For equal spins this reduces to the familiar
full-space dimension `(2S+1)^NOS` only when no states are excluded.

When `[H, S^z_total]=0`, a fixed occupation selects an exact magnetization
sector. If transverse fields, unequal `Jx,Jy`, or other nonconserving terms
couple that sector to excluded states, the calculation instead solves `P H P`,
where `P` projects onto the retained states. Convergence with respect to the
occupation window must then be checked independently of eigensolver convergence.

## Spatial symmetry projection

Up to six finite-order site permutations are supplied explicitly, each with
an order `L_m`, sector `M_m`, and site-image list. The code finds representative
product states, computes their orbit normalization, and assembles the matrix
in the resulting symmetry-adapted basis. Let `T_m` actively move the spin
occupation at each site to the site specified by its image list. For that
active convention the sector eigenvalue is

$$
T_m|\Psi\rangle=e^{+2\pi iM_m/L_m}|\Psi\rangle.
$$

Translations and compatible point-group operations can be used when they
preserve the full Hamiltonian and the retained basis. Site-dependent spins,
fields, oriented DM vectors, and the local occupation caps all enter this
condition. Geometry alone is insufficient, and ED2 does not automatically
verify commutation with `H`.

The ordered products enumerated by ED2 must cover the generated finite group
uniformly, and the assigned scalar phases must define a consistent
one-dimensional character, including all relations between generators.
Commuting generators with compatible characters are a sufficient case;
commutativity itself is not required. ED2 does not automatically validate
these group and character compatibility conditions.

For a translation `T` and reflection `R`, the relation `R T R^-1 = T^-1`
requires a scalar translation character equal to its inverse: momentum zero,
or pi when available. The zero-momentum/even-reflection ex9 sector is
therefore compatible even though its translation and reflection permutations
do not commute. Generic nonzero momentum combined with reflection requires
a multidimensional representation, which these scalar labels do not
implement. Defining the translation as the inverse/passive permutation would
instead give the negative sign in the phase convention above.

There is also no separate spin-inversion parity input in this version.
Symmetry reduction factors depend on orbit sizes and the chosen character;
division by the number of sites is only a rough estimate in suitable
translation sectors.

## Hamiltonian conventions

The implemented model has XYZ exchange, directed DM coupling, symmetric
off-diagonal exchange, and a **positive** field coefficient `+h_i.S_i`. A
physical Zeeman coupling `-B_i.S_i` therefore requires `h_i=-B_i`. The complete
bond-column order and orientation rules are in the
[input reference](../user-guide/input-format.md#hamiltonian-lists-and-sign-conventions).

Off-diagonal channel coefficients are tested against the implementation's
absolute threshold `1e-15`. Terms below that threshold are omitted during
physical row generation. Optional CSR conversion adds no further drop
tolerance and preserves occupied records and their order.

## Matrix storage and multiplication

For `ALG=1` and `ALG=2`, `wk_dim` controls how many initial matrix rows are
cached. The remaining rows are generated each time the Hamiltonian acts on a
vector. `MNTE` sets the padded row capacity, not the Hilbert-space dimension.

The optional coefficient bound, complex CSR cache, padded-cache release, and
two-pass direct CSR build change storage and construction. They do not repair
an invalid physical input or symmetry. Direct CSR requires a full cache and
performs two physical row-generation passes. See
[performance and memory](../user-guide/performance-and-memory.md) for the exact
prerequisites and the distinction between cache-array payload and process RSS.

For `ALG=3`, ED2 explicitly constructs the dense complex Hermitian matrix and
calls LAPACK `ZHEEVR`. Dense storage grows quadratically with the projected
dimension and includes eigensolver work arrays and eigenvectors.

## Iterative solvers and their limits

`ALG=1` is a conventional Lanczos implementation for one lowest state. It
stores a small number of full vectors, diagonalizes the tridiagonal Krylov
matrix, and reconstructs the eigenvector in another pass. Its stopping test
uses relative changes between successive lowest Ritz estimates, not an
explicit residual criterion. The first Ritz estimate is recorded but cannot
itself satisfy convergence because there is no previous estimate. Exact zero
Krylov residuals are handled without normalizing a zero vector; an exactly
one-dimensional invariant subspace retains its complete normalized vector.
The implementation does not perform full basis reorthogonalization in this
path.

A small nonzero Krylov residual relative to `||H q||` alone can be misleading
when a large constant energy offset hides the spectral width. Before claiming
early closure, `ALG=1` also requires that residual to be small relative to the
variation already observed in the diagonal Lanczos coefficients and previous
off-diagonal coefficients. A positive near-breakdown residual in the first
iteration has no such observed scale and causes a nonzero error exit. Exact
zero breakdown and reaching the complete projected dimension retain their
normal completion paths. This conservative check prevents an identified
false-success case; it is not a general residual-based convergence guarantee.

Even a scalar Hamiltonian may therefore stop if roundoff leaves a nonzero
first-step residual. Remove an avoidable constant energy offset, validate
against an independent solver, or use `ALG=3` when the dense matrix fits in
memory. Do not interpret this diagnostic as a valid computed ground state.

`ALG=2` implements thick-restart Lanczos. `NOM` is the Krylov-space size and
`NOK` is the number of retained Ritz vectors. The implementation performs
reorthogonalization and checks relative changes in the retained Ritz energies.
It computes `NOE` reported eigenpairs, subject to the workspace and projected
dimension constraints listed in the
[input reference](../user-guide/input-format.md#input_trlan-alg2).

Both methods use random initial vectors, with no user-facing seed parameter
in the public NAMELIST. Their relative energy-change denominators are bounded
below by `tiny(1d0)` to avoid division by zero. Both print explicit residual
norms `||Hψ-Eψ||`, which must be assessed separately from the stopping test.
Finite precision, loss of orthogonality, nearly zero energies, and iteration
limits still require attention.

The ordinary relative Ritz-energy stopping test remains sensitive to a large
constant energy shift. The additional near-breakdown guard above addresses
one false-closure failure; it does not make all Lanczos convergence tests
shift-invariant. Assess residuals against the energy scale and accuracy
relevant to the observable, rather than against a large constant offset.

Exhausting the configured iteration or restart limit before convergence
causes a nonzero error exit. A complete Lanczos basis or a sufficiently
resolved closed Krylov subspace can terminate successfully without an
energy-change check.

The TRLan implementation stops with a nonzero status if it encounters a
numerically vanishing Krylov direction. It does not construct additional
independent vectors to complete a degenerate invariant subspace. Use dense
`ALG=3` for small problems in this regime. An unchanged energy can coexist
with a poor eigenvector: check normalization and residuals, compare with
independently constructed small-system spectra, and distinguish a single
vector from the whole degenerate eigenspace.
The [output guide](../user-guide/outputs.md#standard-output-and-convergence)
explains what the printed diagnostics actually measure.

## OpenMP and numerical reproducibility

ED2 uses shared-memory OpenMP and has no MPI decomposition. BLAS may have its
own threads; control those separately to avoid nested oversubscription. Record
the compiler, optimization flags, BLAS backend, and thread counts. Mathematical
agreement does not imply bitwise reproducibility across compilers, `-Ofast`
reassociation, or parallel reductions.

The bundled example smoke tests are useful for detecting execution failures.
A scientific validation also needs numerical comparisons: spectra, residuals,
operator conventions, and observables for suitable independent reference
problems. The implementation and the archived source revision remain the
reference for the available behavior; a manuscript under preparation does not
establish that every proposed method is present in the public code.
