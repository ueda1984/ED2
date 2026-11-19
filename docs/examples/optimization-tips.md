# Lattice- and model-specific optimization tips

These tips complement the [general cache and memory controls](
../user-guide/performance-and-memory.md). They concern **physical
model choice and valid symmetries**, not universal speedup promises.
Always recheck `[H,G]=0` for every requested generator `G`;
including a permutation that does not commute with the model can
produce a different projected operator.

## One-dimensional spin-1/2 chain

- For an XXZ Hamiltonian with a field parallel to $z$ and no
  spin-nonconserving terms, work at fixed
  $N_{\downarrow}$ using `NODmin=NODmax`. The sector
  dimension for spin-1/2 is $\binom{L}{N_{\downarrow}}$.
  This is a **conserved sector**, not an approximation.
- Use lattice translation only if bond couplings and fields
  are periodic with the **same** translation. For a
  dimerized chain, translation by one site is generally
  not a symmetry; translation by the larger unit cell
  may be.
- `ED2_CACHE_COEFF_BOUND=1` can reduce the padded
  transition capacity for XXZ exchange when `MNTE<0`,
  because many anisotropic transition channels have
  exactly zero coefficients.

## Alternating mixed-spin chain

For an alternating spin-1/2/spin-1 chain with (L) unit
cells, the full Hilbert space has dimension $6^L$.
Translation by **one unit cell** preserves the spin
magnitudes; a one-site shift generally does not.

An XYZ model with unequal (J_x,J_y), transverse
field, or oriented DM interaction generally breaks
conservation of total (S^z). In that case:

- Use the **full physical occupation window** if the
  goal is the unrestricted spectrum. A smaller
  `NODmin:NODmax` window represents a **projected
  Hamiltonian**, not a symmetry-conserved sector.
- Verify the orientation of every DM bond. Reversing
  the directed bond changes the sign of the
  antisymmetric exchange.
- For modest Hilbert dimensions, compare native
  energies and $H|\psi\rangle$ with an independently
  assembled dense Hamiltonian.
- For a fully cached sparse operator, compare regular
  CSR and direct two-pass CSR, including setup time and
  memory use. They have different tradeoffs.

See the [mixed-spin reference example](mixed_spin_chain.md).
The example's parameters and small local numerical tests should not
be conflated with the full physical Hilbert-space tests used in
the development CI.

## Two-dimensional frustrated lattices

Square, triangular, kagome and honeycomb models may exploit
translations and compatible point operations when the
**entire anisotropic interaction tensor** respects them.
A lattice that is geometrically symmetric is not
necessarily a symmetry of its Hamiltonian.

- Start with translations preserved by the coupling pattern.
  Then test whether additional point operations commute with
  the model and with the chosen symmetry-sector projectors.
- Be especially careful with bond-dependent anisotropy,
  oriented DM vectors and nonuniform fields.
- At fixed small magnon number, a few-particle basis
  can be vastly smaller than the full spin-1/2 space.
  This is exact only when magnon number is conserved.
- Keep model-specific input/reference fixtures and use
  an independent small-cluster calculation to validate
  symmetry and sign conventions.

## Simple-cubic dilute-magnon example

For a periodic simple-cubic lattice with $N=L^3$
spin-1/2 sites and **exactly three conserved magnons**,
the unprojected fixed-magnon dimension is
$\binom{N}{3}$. A 10×10×10 lattice has
$\binom{1000}{3}=166{,}167{,}000$ such basis states
before lattice-symmetry reduction.

A development-only experiment on the $\Gamma$
point and four one-dimensional $O_h$ sectors found
projected dimensions 4,299; 3,624; 2,884; and 3,252.
These are **specific to the verified isotropic
few-magnon setup**, not generic input/output
guarantees for a field-anisotropic cubic example.
Noncommuting point-group generators must be handled
with an explicitly verified representation, not
assumed to be available in the stable generic
commuting-generator interface.

The audited PRESTO coordinate-cache comparison was
performed on a separate development stack,
not on the current generic `main` release. Its
numeric equivalence is useful evidence, but the
exclusive-resource attestation and the precise
software release target remain unresolved.
Accordingly, its timing ratios are **not** recommended
as general CPC performance figures.

The stable [cubic example](cubic.md) includes
an anisotropic XYZ/DM model. Its symmetries should
be checked against its actual couplings; do
**not** transfer the isotropic $O_h$ example's
sector dimensions or benchmark settings to it.

## Benchmarking checklist

When reporting a lattice-specific optimization,
record the following in the same run package:

1. Exact repository and commit, input files and checksums.
2. The model Hamiltonian, spin lengths, boundary
   conditions, physical conservation laws, symmetry
   generators and the **actual projected dimension**.
3. Compiler and flags, BLAS/LAPACK version, thread
   bindings and number of independent solver starts.
4. Memory representation (`wk_dim`, `MNTE`, environment
   cache flags), initial cache-setup time, H-action count,
   Lanczos time, total time, and peak RSS if measured.
5. Energy, direct residual, and an independent physical
   reference when available.
6. Hardware allocation/affinity and raw per-run results.
   Shared CI runner times are useful for exploratory
   comparisons but not equivalent to an attested
   dedicated-node benchmark.

For reproducible general-model performance guidance,
use [Hamiltonian caching, performance, and memory](
../user-guide/performance-and-memory.md).
