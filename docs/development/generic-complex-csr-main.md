# Complex CSR for general ED2 — original-main integration

Status: independent main-based integration branch, not the older
cubic-specific PR13–21 stack. The existing ED2 default path is unchanged.

## User interface and scope

Use the exact opt-in environment variable:

    ED2_CSR_COMPLEX_CACHE=1 ./source/QS3ED2 < input.dat

The option is OFF by default (`unset` or `0`). The program
requires `wk_dim == THS`, a **fully cached** Hamiltonian.
It fails closed if a partial cache is requested, the setting is
not literally `0` or `1`, a target is outside the basis,
the row offset would overflow the default 32-bit integer,
or an allocation fails.

This option applies to both **exactly real** and **genuinely
complex** Hamiltonians: Dzyaloshinskii--Moriya interactions,
Cartesian fields and arbitrary model-compatible projected
Hilbert sectors are not excluded by design.

The original padded arrays `wk_loc` and `wk_ele` remain intact
and allocated; a **separate** CSR-like array copies every
non-sentinel row entry in original order, *including stored zero
coefficients*. We neither fuse duplicate indices nor reorder
floating-point operations. The existing complex Krylov vector
path is unchanged.

Therefore this is an extra-memory versus repeated-matvec
tradeoff, **not** a reduction in the program's total memory
allocation. The native startup output reports:

    Hamiltonian action: csr-complex cache entries=N, extra_bytes=B

The additional storage is `4*(THS+1)+20*N` bytes when default
`integer` is 32 bits and `complex(8)` is 128 bits. The source
computes actual storage bytes with `storage_size` rather than
hardcoding these sizes in the executable.

## Scientific integration proof

Unlike the prototype PR25, this patch is based **directly
on the original `ueda1984/ED2` main commit**
`49941dc4b437a8fb6c2b5c8dac7d621866731d7c`
and only introduces:
- `source/ham2vec.f90`: guarded one-time CSR construction
  after original Hamiltonian cache creation, plus a
  multiply loop that preserves the exact old complex
  accumulation order. Old branch used for opt-in OFF.
- `tests/cpc_csr_matvec_driver.f90`: independently compiled
  native driver that builds the **original** ED2 basis and
  Hamiltonian, applies H to *two deterministic complex vectors*,
  and records the entire real/imaginary IEEE bit pattern.
- `tests/cpc_csr_main_regression.py`: standard-library
  input generator for the alternating mixed spin-1/2/spin-1
  OBC non-U(1) XYZ+hx model at N4/6/8,
  with both real and truly complex DM terms. It verifies
  identical generated input SHA256, identical full H*x
  SHA256, exact output vector length and full-cache
  provenance, and checks old standalone TRLan spectra
  against a separately constructed dense Hamiltonian H for
  N4/N6, with an absolute complex-vector tolerance of 1e-10.
  The original TRLan can stop before accurate convergence even
  when started from the same seed, so *its* final eigenvalues
  are deliberately not used as the opt-in H*x proof.
  Every individual native H*x bit is still required to match
  the legacy cached matvec exactly; the dense reference uses
  a distinct ED2 full-H construction path.
  Partial cache and malformed/oversized flags are
  independently checked to fail closed.
- `.github/workflows/main-generic-csr.yml`:
  two-compiler original-main build, old `make check`
  and `make check-long`, full native H*x reference
  comparisons and archived run logs.

The branch does **not** import PR17–21 cubic source
modifications, newer two-stage Krylov enhancements,
permuted basis builders, program-format changes,
new public ABI or changes to the CPC manuscript.

## Relative priority of other experiments

- PR22 direct-real multiply is separate and may be useful
  with a rigorously real fully cached Hamiltonian but
  does not handle complex DM generically.
- PR23 padded dense-real copy, PR24 CSR-real and PR25
  prototype CSR-complex remain valuable **comparative
  experiments** and independent evidence. Do not merge
  their stacked branches wholesale into main.
- The generic CSR feature is intentionally opt-in
  until repeated controlled, hardware-bound performance
  studies and memory-budget reviews justify a
  default-enabled adaptive policy.

No PRESTO access is required for the correctness tests.
The earlier PRESTO operator exclusivity attestation is
unrelated to this scientific equivalence proof.
