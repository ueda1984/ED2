# Experimental CSR-like row storage for exactly real cached ED2 Hamiltonians

**Status: opt-in Draft, not a merged or released optimization.**
This branch is independently stacked on the general real-operator
PR #22. It does **not** require the separate dense packed-real
PR #23 and does not touch CPC manuscript PR #20, upstream, or main.

## Why this representation is different

ED2's historical cached Hamiltonian uses padded 2D arrays,
`wk_loc(j,i)` and `wk_ele(j,i)`, retaining an often significantly
larger fixed column capacity than the actual number of stored row
terms. The faster dense real cache in experimental PR #23
copies the **entire** padded matrix into a separate real array,
trading extra bytes and construction time for better strides.

For general large alternating-spin non-U(1) Hamiltonians,
a separate **CSR-like ordered row representation** can provide
an even smaller working set without changing the numerical
sum. The CSR-like data consist of:
- `csr_rowptr(1:lv+1)`: one-based offsets for each cached row,
- `csr_target(1:nentries)`: all occupied original target indices,
- `csr_coeff(1:nentries)`: the exactly real coefficient of
  each occupied term.

**All terms, including any explicit numerically zero coefficients,
are retained in the original row order.** The transformation only
removes unoccupied sentinel/padding slots. This is important:
ignoring terms based on numerical tolerance could change IEEE
sum results, and is deliberately not done.

The original `wk_ele` and `wk_loc` remain allocated for
compatibility with other ED2 functionality and the independent
full-matrix proof. Thus CSR is an optional **extra allocation**,
not a reduction of the program's total memory footprint.

## Exact eligibility and fail-closed behavior

Set `ED2_REAL_CACHED_MATVEC=1` and
`ED2_CSR_REAL_CACHE=1` explicitly. Both default OFF.
The new CSR flag rejects missing real-mode prerequisite,
invalid/oversized values, partially cached projected bases
(`wk_dim != lv`) and **any nonzero imaginary part** in
the full cached operator. The original real-mode guard is
performed before CSR construction.

The code checks all target index ranges, 32-bit row-offset
capacity, array allocation status and each packed row count.
The exact additional memory is printed by the native source:

    Hamiltonian action: csr-real cache entries=<n>, extra_bytes=<m>

Additional bytes = `4*(lv+1)+12*nentries` with the current
32-bit integer and 64-bit real build. This is platform-dependent
if default integer kind were changed.

At every H action, each row is recomputed in its original
sequence via `csr_coeff(j)*v1(csr_target(j))`; the vectors
remain **complex**, the arithmetic remains exact real-times-complex
and the native two-stage Krylov solver is unmodified. This is
a storage/lookup optimization, not a new numerical algorithm.

## Numerical acceptance

Checked-in `tests/cpc_csr_real_cache_regression.py` runs the
actual native ED2 cubic cache driver and compares **every
complex output IEEE-754 bit** for a fixed, nonreal complex input
vector across:
- L=4,6,10 × four cubic 1D O_h irreps × T=1,4,
- L=10 × the same four irreps × T=8:
  **28 complete H-vector pairs per compiler**.

The test independently checks native source mode markers, input
SHA256, expected projected dimension, the exact count of
bit-pattern lines and the reported CSR memory formula.
Additional generic mixed-spin full-space real-Hamiltonian
tests require eigenvalue-file SHA256 equality; complex DM,
partial cache, missing real mode, invalid values and
oversized flags must all be rejected.
All earlier PR13–22 basis, Hamiltonian cache, Lanczos,
symmetry and independent model tests remain in the
two-compiler GitHub Actions suite.

Independent temporary local GNU14.2 prototype already
passed 28/28 complete complex-vector SHA256 comparisons
and 80 generic full native ED2 trial runs across N8/N10
at T1/T4 with all four operator forms (original complex,
PR22 real multiply, PR23 dense packed real, CSR real).
The full eigenvalue-file SHA256 and H-application counts
matched within every native case. The exact PR-head CI
still needs to complete independently.

## Performance *direction only*: shared Linux, not PRESTO

A short GNU14.2 shared-host, 5-interleaved-repeat pilot was
run with the same input, seed and `-O2 -fopenmp` binary
on the **genuinely non-U(1)** open alternating
spin-1/2/spin-1 XYZ + transverse hx model. The Hamiltonian
is real in the Sz basis and physical U(1) symmetry is broken.

Warmups/repetition variability on shared CPUs remain
uncontrolled. Values are wall-time medians in milliseconds
from a **temporary local** source instrumentation around
all native Hamiltonian applications and the standard ED2
phase timer:

| D (N) | T | Complex action | Real-only action | Dense packed action | CSR action | Complex full | CSR full |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1296 (8) | 1 | 6.972 | 6.470 | 5.752 | 5.259 | 12.772 | 11.470 |
| 1296 (8) | 4 | 4.185 | 4.607 | 4.751 | 4.532 | 9.564 | 9.299 |
| 7776 (10) | 1 | 94.259 | 94.673 | 74.807 | 45.805 | 137.842 | 95.280 |
| 7776 (10) | 4 | 46.667 | 49.464 | 32.255 | 27.025 | 71.272 | 67.315 |

For N10/T1 the measured *action itself* is ~2.06x
faster and full run ~1.45x faster than original complex
storage, but this **is not hardware-exclusive,
publication-grade or a general guarantee**. In N8/T4,
CSR's H-action is slightly slower than the original
(4.532 vs 4.185 ms). T4 full-time effects can
be much smaller than H-action improvements.

CSR stores **191,808** row entries at N10/D7776.
Exact extra memory: `12*191808 + 4*7777 =
2,332,804 bytes` (~2.22 MiB), less than the
padded dense packed real array's
`8*93*7776 = 5,785,344 bytes` (~5.52 MiB).

A separate N1000/A1g cubic Gamma/O_h prototype
with 160 H actions gave an exploratory median
T1 action **104.180 -> 58.285 ms** and total
**232.763 -> 179.902 ms**; at T4 action
**45.302 -> 36.556 ms**, total **154.836 ->
151.665 ms**. There the dense packed real
path was about as fast at T4 and uses less
extra memory: CSR had 235,865 occupied entries,
requiring `2,847,580` extra bytes versus
`2,235,480` for dense real copy at MNTE=64.

The correct choice is **model-, thread- and
memory-dependent**. No automatic threshold, statistical
speedup, HPC reproducibility or release default is
claimed. PR23 and this alternative PR are intentionally
kept separate until enough controlled evidence
exists to choose one.

## Next checks

1. Both GNU13/14 CI full complex-output and general
   mixed-spin negative-case suites pass on immutable PR head.
2. Evaluate the CSR representation in a controlled,
   **non-PRESTO** physical benchmark environment with
   thread affinity, core exclusivity and raw trial retention.
   Do not suggest PRESTO job submission while access is
   unavailable.
3. Compare memory use, cache initialization and solver
   repeated H actions for several model families, not
   just one cubic three-magnon example.
4. Consider a future **shared representation** that frees
   the original complex storage only after all public
   module contracts and independent H-cache/dump
   interfaces are proven compatible. That would be a
   new, separately reviewed engineering change.

No release, merge or manuscript performance claims.
