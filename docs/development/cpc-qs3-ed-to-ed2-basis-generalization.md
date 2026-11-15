# From QS³'s binomial basis to ED2's bounded-occupation representation

**Source-audited, analytical note for the first ED2 program article.**
The published QS³ code is
[QS-Cube/ED](https://github.com/QS-Cube/ED),
not the publicly hosted but journal-unpublished ED2.
This derivation is about **basis representation**,
not a new eigensolver, numerical speedup or a proof
that every truncated nonconserving sector is invariant.

## 1. Exact common limit: binomial counting

In the old `QS-Cube/ED/source_only_u1/state_lists.f90`,
the basis of \(L\) spin-1/2 sites with a
fixed number \(N_\downarrow=n\) of down spins
has

\[
\mathcal D_{\mathrm{QS^3}}(L,n)=\binom{L}{n}.
\]

The original source constructs binomial-coefficient
tables (`get_combination`) and a rank/unrank mapping
over the sorted down-spin site list
(`inv_list`, `list_fly`). No general-purpose
bit-string enumeration is required.

This limit is **directly tested** against
the actual old and new Fortran executables
in the [32-level matched XXZ crosscheck](cpc-qs3-ed-to-ed2-scope-matrix.md#31-extension-into-the-original-near-saturation-sector-and-larger-conserved-sector).

## 2. Generalization to different site occupations

For ED2, let the local spin be \(S_i\) and
the allowed local down-spin occupation be
\(n_i=0,\ldots,m_i\), where
\(0\le m_i\le 2S_i\).
With no site-specific truncation,
\(m_i=2S_i\) and the local dimension is
\(2S_i+1\). Define

\[
P_L(x)=\prod_{i=1}^{L}\left(\sum_{k=0}^{m_i}x^k\right),
\qquad
a_L(n)=[x^n]P_L(x).
\]

Then the **pre-spatial-symmetry** basis
size for an allowed global occupation window
\(N_{\min}\le n\le N_{\max}\) is

\[
\boxed{
\mathcal D_{\mathrm{ED2}}=
\sum_{n=N_{\min}}^{N_{\max}}
a_L(n).
}
\]

The exact coefficient recurrence is

\[
a_i(n)=\sum_{k=0}^{m_i} a_{i-1}(n-k),
\quad a_0(0)=1,\quad a_0(n\ne0)=0.
\]

It has a sliding-window implementation
with constant-time updates per coefficient
after preprocessing. Inspection of
[`QS-Cube/ED2/source/state_lists.f90`](https://github.com/QS-Cube/ED2/blob/main/source/state_lists.f90)
confirms a bounded-repetition (`ltd_rep_combination`)
table and a cumulative/sliding-window
(`cum_ltd_rep_combination`) table. Those
tables support state ranking, unranking
and selection of a global occupation window
without enumerating the entire tensor-product
basis as explicit bit strings.

For full local occupation and the entire
global window this expression reduces to

\[
\mathcal D_{\mathrm{full}}
=P_L(1)=\prod_{i=1}^L(2S_i+1).
\]

For all \(m_i=1\) it reduces to
\(P_L(x)=(1+x)^L\), hence
\(a_L(n)=\binom{L}{n}\): precisely
the old QS³ fixed-occupation count.
That is an **exact mathematical
specialization**, independently
of any solver convergence tolerance.

## 3. Explicit coefficient checks

The following integer coefficients were
independently obtained from the above
polynomial product (no diagonalization or
ED2 solver residual is involved).

| Representation | \(P_L(x)\) occupation caps | Selected exact counts | Full product dimension |
| --- | --- | --- | ---: |
| Old 4-site spin-1/2 | (1,1,1,1) | \(a_4(2)=6\) | 16 |
| ED2 alternating 4-site spin-1/2,1 | (1,2,1,2) | \(a_4(3)=10\) | 36 |
| ED2 alternating 6-site spin-1/2,1 | (1,2,1,2,1,2) | \(a_6(4)=48\) | 216 |
| ED2 alternating 12-site spin-1/2,1 | (1,2) repeated six times | \(a_{12}(9)=7780\) | 46,656 |

The mixed four-site polynomial has coefficients
\((1,4,8,10,8,4,1)\);
for six sites it is
\((1,6,18,35,48,48,35,18,6,1)\).
These are exact integer checks, **not**
measurements of computation speed.
The 4- and 6-site full-space dimensions
are also independently exercised by the
[portable ED2 DM quickstart](cpc-portable-dm-quickstart.md).

## 4. Crucial physical boundary: invariant versus projected

The global total occupation operator
\(N_\downarrow=\sum_i(S_i-S_i^z)\)
is conserved in the original
spin-1/2 XXZ Hamiltonian,
so the old fixed-\(N_\downarrow\)
subspace is **exactly invariant**.

With nonzero transverse fields,
\(J_x\ne J_y\), DM terms not
aligned with \(z\), or symmetric
off-diagonal interactions,
the full general ED2 Hamiltonian
can connect distinct \(N_\downarrow\)
values. Restricting a window then
usually computes a **projected Hamiltonian
\(PHP\)**. Its spectrum need not
equal a symmetry sector of the
unrestricted \(H\), and it is
incorrect to market that reduction
as exact global \(S^z\) conservation.

Similarly, site-permutation character
projections only define valid
invariant symmetry sectors when the
specified operations actually
commute with the Hamiltonian.
The present software does not
certify arbitrary user-supplied
permutations automatically.

For the paper's main mixed-spin
XYZ/Gamma/DM response,
we therefore use the **full**
\(N_\downarrow=0,\ldots,\sum_i 2S_i\)
space and identity spatial sector.
The demonstrated 12-site dimension
is 46,656, not the artificially
smaller 7,780-dimensional
fixed-\(N_\downarrow=9\) projection.

## 5. CPC Methods use

The generalized occupation-counting
method and its exact QS³ binomial limit
should be a distinct part of the ED2
implementation description, rather than
only listing supported spin magnitudes.
It explains the representation change
that permits mixed \(S_i\) and
nonconserving interactions. Cite the
2022 QS³ publication for the original
binomial/momentum framework; identify
the generalized bounded-repetition
recurrence as an ED2 implementation
feature already present in the
public ED2 code, not an invention
of later numerical-maintenance tests.

Do not compare maximum lattice site
counts at fixed \(N_\downarrow\)
and unrestricted mixed-spin dimensions
as if the computational targets
were equivalent.
