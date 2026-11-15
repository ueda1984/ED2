# Direct iterative and translation-sector continuity: published QS³/ED -> ED2

**Scope and status:** original published source
[QS-Cube/ED](https://github.com/QS-Cube/ED), exact old
SHA `a3b09934bba5f966bfb564718528de8778c2b46c`,
versus the first proposed ED2 program article, based on
[QS-Cube/ED2](https://github.com/QS-Cube/ED2)
with the reviewed development tests.

The article's primary comparison is to **published QS³/ED**,
not merely to ED2's pre-existing GitHub implementation.
This is continuity of shared *physical problems*,
not a claim of identical input ABI, old dynamic-spectrum
output or new asymptotic solver theory.

## 1. What the previously verified full spectra showed

The [existing dense test](cpc-qs3-ed-to-ed2-scope-matrix.md)
compiles the genuine 2022 source and the current ED2 distribution,
and compares a shared spin-1/2 OBC XXZ model
\(J_{xy}=0.8,J_z=1.1\) in three
fixed-down-spin sectors. With ALG=3
dense LAPACK on both programs, all
**32/32 eigenvalues** agree against
an independent full product-sector matrix;
the maximum old-to-new discrepancy is
\(2.89\times10^{-15}\).
See the [frozen data](cpc-qs3-ed-ed2-32level-extended-ci.csv)
and [passing CI](https://github.com/ueda1984/ED2/actions/runs/43712595377).

## 2. NEW: both actual conventional Lanczos executables

The exact old QS³ `source_only_u1/QS3_only_u1`
and ED2 `source/QS3ED2` now run **ALG=1**,
each independently finding the lowest
energy in the same OBC Hamiltonian
at \(N_\downarrow=2\).
We also diagonalize an independent,
product-bit-basis \(N_\downarrow=2\)
Hamiltonian to avoid circular
agreement between the two Fortran codes.

**Verified [GitHub Actions run #43768984843](https://github.com/ueda1984/ED2/actions/runs/43768984843)**,
GNU Fortran 13, Ubuntu 24.04, OpenBLAS;
old code SHA
`a3b09934bba5f966bfb564718528de8778c2b46c`,
ED2 tested SHA
`1e9df6b5c357daea85e59cd3501c9a769f291356`,
OpenMP=1. Both executables
were built without modifying
the old Fortran source code.
The test and data are reproducible with
[`tests/cpc_qs3_ed_iterative_xxz.py`](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-stationary-dm-scaling/tests/cpc_qs3_ed_iterative_xxz.py).

| \(L\) | \(\dim\binom{L}{2}\) | Independent \(E_0\) | Old QS³ \(E_0\) | ED2 \(E_0\) | Absolute difference | Reported Lanczos iterations (old/new) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 66 | -0.6225631337454959 | -0.6225631337454958 | -0.6225631337454959 | 1.11e-16 | 35 / 35 |
| 18 | 153 | +0.9406443730741343 | +0.9406443730741542 | +0.9406443730741364 | 1.78e-14 | 45 / 55 |
| 24 | 276 | +2.5612930570490877 | +2.5612930570490990 | +2.5612930570490870 | 1.20e-14 | 60 / 60 |

**3/3 original-vs-ED2 iterative cases passed**.
The largest old-to-independent
difference is \(1.99\times10^{-14}\);
the largest ED2-to-independent
difference is \(2.11\times10^{-15}\).

The original QS³ algorithm uses
uncontrolled random starts, while
the modern solver has different
numerical handling and convergence
checks. Therefore the iteration
counts are **recorded context**
and are *not* performance ratios
or measurements under identical
initial vectors. The exact
Hamiltonian is the controlled
scientific target; Krylov
trajectories are not matched.
The purpose is **scientific
continuity in the iteration
mode previously important to
QS³**, not proof of a speedup.

## 3. Mathematically fair comparison boundaries

- The old fixed-\(N_\downarrow\)
  basis is a genuine invariant
  subspace of the shared XXZ
  Hamiltonian. ED2's setting
  \(N_{\min}=N_{\max}=N_\downarrow\)
  recovers that exact sector,
  rather than projecting out
  physically coupled states.
- The models have identical
  \(J_x=J_y=J_{xy}\),
  \(J_z\), open boundaries and
  zero DM, Gamma and site fields.
  These restrictions are chosen
  solely to establish the
  historical overlap.
- The central ED2 scientific
  example remains the
  **non-U(1), alternating
  spin-1/2/spin-1 XYZ/Gamma/DM
  chain**, verified separately
  through 12 sites. That
  model is beyond the published
  QS³ input Hamiltonian, and
  cannot be forced into a
  misleading old-code
  performance comparison.
- Timing or peak-RSS
  superiority cannot be
  established by these
  source-archived shared
  hosted-runner jobs. The
  controlled-node benchmark
  gate remains described in
  [the dedicated-node protocol]
  (cpc-controlled-node-protocol.md).

## 4. Translation sectors under review

The original QS³ also exploits
crystal momentum. A separate
[pinned-source translation test](https://github.com/ueda1984/ED2/blob/codex/ed2-cpc-stationary-dm-scaling/tests/cpc_qs3_ed_translation_sector.py)
constructs every one-site
translation sector of an
\(L=6\), \(N_\downarrow=2\)
periodic XXZ ring with \(k=0,\ldots,5\)
and independently constructs
the projector

\[
P_k = \frac{1}{6}\sum_{r=0}^{5}
e^{-2\pi i k r/6}T^r .
\]

This is designed to verify
all 15 basis energies after
momentum resolution. **Do not
quote a numerical success
until its own CI has passed.**
The predecessor's dynamical
spin-structure-factor path
is not included in this
static spectral equivalence
check.
