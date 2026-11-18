# Experimental general complex CSR-like cache for ED2 H-vector products

**Status:** default OFF; separate Draft experiment stacked on
PR #22's complete native H-vector proof machinery. It is
an **alternative** to both PR #23's dense real coefficient copy
and PR #24's compressed real coefficient copy. No production
merge, tag, PRESTO allocation or CPC manuscript claim.

## Broad model support

Unlike the real-only experiments, this optional compressed
complex cache can handle **every fully cached ED2 Hamiltonian**,
including models with nonzero imaginary matrix elements from
oriented Dzyaloshinskii--Moriya exchange and y magnetic fields.
It does not attempt to change complex arithmetic, point-group
symmetry, truncated Hilbert-space rules, convergence checks,
observable definitions or Krylov solver.

Opt in explicitly with `ED2_CSR_COMPLEX_CACHE=1`. It requires
`wk_dim=lv`, i.e. a complete cached projected matrix, and
rejects invalid/oversized values. To prevent ambiguous overlapping
experimental kernels, it also rejects combining with
`ED2_REAL_CACHED_MATVEC=1` in this Draft. Both remain OFF by
default. The preexisting `wk_ele`/ `wk_loc` arrays are kept,
so **the new arrays add memory** to the original implementation.

Ordered CSR-like row storage:
- `complex_csr_rowptr(1:lv+1)`, one-based cumulative row offsets;
- `complex_csr_target(1:nentries)`, every original row target;
- `complex_csr_coeff(1:nentries)`, exactly copied complex numbers.

Every original non-sentinel record is retained in order,
**including explicitly stored complex zeros and both signed
zero components**. No coefficient is removed due to magnitude.
Each H action starts with the historical complex zero and
accumulates exact original `coeff*v1(target)` terms in original
order, under the same parallel row ownership.
Therefore no matrix or numerical physics is intentionally
approximated.

The program fails closed on a partial cache, any invalid
target, 32-bit row offset overflow and allocation failure.
Its provenance marker reports the exact new storage cost:

    Hamiltonian action: csr-complex cache entries=<n>, extra_bytes=<m>

For standard default 32-bit `integer`, complex(8) 16 bytes:
extra bytes `4*(lv+1)+20*nentries`.

## Scientific proof (before a merge)

`tests/cpc_csr_complex_cache_regression.py` runs actual
compiled ED2, with full SHA256 of the **entire complex
H*psi output vector** for deterministic nonreal-complex
trial vectors. Both old and new use byte-identical input.

- Cubic exact Gamma/O_h, all four 1D irreps,
  L4/L6/L10 T1/T4 and L10 T8: **28 complete
  output-vector bitwise pairs/compiler**.
- General alternating S1/2,S1 full-space model, both
  exactly real and genuinely complex (DM=0.12),
  L4/L6 at T1/T4: **8 native complex-vector bitwise
  pairs/compiler**.
- Native two-stage spectra for 6-site real and
  complex DM models must produce identical
  `eigenvalues.dat` SHA256, with certified residual.
- Negative tests reject partial matrix cache,
  malformed/oversized flags and conflicting real mode.
- Original PR13–22 native basis, H matrix, Ritz, DM,
  symmetry and conservation correctness tests remain
  required in both GNU13/GNU14 CI jobs.

An independent local GNU14.2 prototype completed 40
separate two-stage generic **complex** DM/full-space
solver runs and another 40 genuinely non-U(1)
**real** XYZ/transverse-field runs. It compared
OFF and ON for each identical input/seed and
confirmed all 80 eigenvalue-file SHA256s and
Hamiltonian-application counts matched. The
PR-head dual-compiler CI is the separate acceptance
gate; never infer completion from local prototypes.

## Shared GNU14.2 pilot: preliminary and un-attested

Five interleaved measurements/mode, T1/T4, 8- or
10-site OBC alternating mixed-spin chains. Both real
and nonreal DM examples are physically non-U(1);
all timing numbers below are native wall-milliseconds.
The H-action timer is a **temporary local probe**
not committed to the PR. No operator-exclusive
resource claim, publication speedup, or predefined
dimension threshold is implied.

| Physical model | Sites | T | Original H-action | Complex CSR H-action | Original total | CSR total |
|---|---:|---:|---:|---:|---:|---:|
| Complex y-field non-U(1) | 8 | 1 | 7.408 | 5.367 | 14.093 | 12.284 |
| Complex y-field non-U(1) | 8 | 4 | 7.369 | 4.154 | 11.831 | 9.284 |
| Complex y-field non-U(1) | 10 | 1 | 83.098 | 39.577 | 130.011 | 91.512 |
| Complex y-field non-U(1) | 10 | 4 | 46.167 | 32.917 | 70.323 | 62.784 |
| Real XYZ + hx non-U(1) | 8 | 1 | 7.316 | 5.067 | 14.380 | 11.510 |
| Real XYZ + hx non-U(1) | 8 | 4 | 5.159 | 4.944 | 10.119 | 9.564 |
| Real XYZ + hx non-U(1) | 10 | 1 | 91.375 | 45.142 | 152.081 | 112.970 |
| Real XYZ + hx non-U(1) | 10 | 4 | 50.007 | 27.405 | 77.050 | 51.905 |

The 10-site mixed-spin Hilbert-space dimension is
7776. The CSR representation contains 191,808
original cached row entries: extra storage
`20*191808+4*7777 = 3,867,268 B` (~3.69 MiB).
The alternate real-only CSR storage in PR24 uses
2,332,804 B (~2.22 MiB) on the same real
Hamiltonian; on a complex Hamiltonian, however,
PR24's exact-real-only guard rightly rejects
and PR25's complex CSR is usable. Both keep
the original padded complex arrays allocated.

In the 10-site complex model at T1, the provisional
H-action ratio is ~2.10x and the full runtime ratio
~1.42x. But shared-host variability, model size,
cache setup and thread count can change the outcome
significantly. A re-run on independent hardware
would be necessary before recommending defaults.

## Decision framework

Prefer retaining the generic complex CSR implementation
as an **experimental, user-opt-in** candidate if dual
compiler and general real/complex DM tests pass. It
provides the broadest scientific model coverage among
these cache-layout experiments. Keep PR24's real CSR
as a lower-extra-memory alternative for real H, and
PR23's dense copy as a potentially useful option on
high occupancy matrices. Do not combine incompatible
flags, merge the chain or promote benchmark figures
to the CPC program-paper claims without further work.

Possible longer-term improvement: after a separate
API/output compatibility audit, allow the compressed
representation to **replace** rather than duplicate the
padded complex cache to reduce total memory; do NOT
make that destructive change without preserving all
current source/input/proof and expectation-value interfaces.
