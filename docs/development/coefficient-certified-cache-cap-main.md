# Coefficient-certified padded Hamiltonian row capacity (general ED2)

**Status:** opt-in general peak-construction *allocation* improvement, independent
of the experimental cubic-symmetry stack and applicable before the mainline
complex CSR conversion. Default behavior is unchanged.

## Purpose and exact mechanism

When input `MNTE < 0`, the original ED2 `make_wk_loc_and_ele` allocates
`wk_loc` and `wk_ele` as dense padded arrays with each row dimension

\[
w_{\mathrm{old}} = 1 + 2N_{\mathrm{one}} + 8N_{\mathrm{two}}.
\]

This bounds the **formal number** of field and exchange transition types,
including transitions whose coefficient is identically zero. It is often
very wasteful for an XXZ or real XYZ spin Hamiltonian, even if the physical
nonzero matrix has many possible basis transitions.

Opt in with:

    ED2_CACHE_COEFF_BOUND=1 ./source/QS3ED2 < input.dat

An input with `MNTE < 0` computes a new conservative row bound *before*
the padded arrays are allocated: one guaranteed diagonal record plus one
possible row record for **each nonzero one-site (2 per field) and two-site
(8 per exchange tensor) transition coefficient**. Every presence test uses
the exact same complex coefficient formula and `abs(c)>eps` test as the
physical row generation below, so no physically allowed transition can
be omitted on account of the new upper bound. The actual matrix generation,
operation order, phase, normalization, basis enumeration, and native H*x
arithmetic are completely unchanged. A symmetry sector can only suppress
further transitions, so this bound remains safe with the general ED2
projection machinery as well.

**The upper bound need not be tight**. If all nine exchange tensor terms
are present, all eight offdiagonal transition channels can occur and
no reduction may be possible. The source detects an impossible bound
before allocation and preserves historical `MNTE too small` checks for
each row. A 64-bit intermediate protects its coefficient count.

Explicit user `MNTE>=0` is intentionally **not** overridden. A machine
marker states that the manual cap was retained. Malformed/oversized
`ED2_CACHE_COEFF_BOUND` values reject rather than silently change
semantics. Setting `0` or leaving the variable unset is identical to
the historical allocation policy.

## Native memory accounting and scope

The source reports:

    Hamiltonian cache: coefficient-bound width=W, nominal_width=N, initial_bytes_saved=B

The saving is based on the allocator payload of `wk_loc` and `wk_ele`,
using the source's actual `storage_size` element sizes and complete
`wk_dim` rows. It is an **initial padded matrix allocation saving**,
which also reduces the memory that coexists with the newly constructed
CSR arrays when `ED2_CSR_COMPLEX_CACHE=1`. It is *not a measured peak RSS*
or virtual address space bound, and it does not eliminate the padded
arrays entirely. The independently optional
`ED2_CSR_RELEASE_PADDED_CACHE=1` can subsequently remove their
long-lived storage after packing, but cannot eliminate CSR construction
overlap by itself.

Local GNU Fortran 14.2 `-O2` example, N10 alternating S=1/2,S=1,
full non-U1 Hilbert space D=7776, OBC:

| Model | Historical row width | Certified width | Initial cache bytes avoided |
|---|---:|---:|---:|
| Real anisotropic XYZ + transverse field | 93 | 57 | **5,598,720 B** (38.7%) |
| Complex DM=0.12 with XYZ/field | 93 | 57 | **5,598,720 B** (38.7%) |
| U(1)-conserving XXZ (zero transverse field) | 93 | 19 | **11,508,480 B** (79.6%) |
| Fully populated 9-component exchange tensor, hx, hy | 93 | 93 | **0 B** |

The nominal original padded matrix payload is 14,463,360 bytes at
this N10/D7776; the actual process RSS includes much more, and these
numbers are **not empirical RSS observations**. For the complex CSR
release mode, a smaller initial padded array lowers the peak
*cache-array payload*, but a measured reduction in peak resident
RAM must be independently audited.

## Exact numerical evidence and compiler policies

The checked-in `tests/cpc_coeff_cache_bound_regression.py`
constructs genuine native ED2 mixed-spin alternating S1/2,S1
Hamiltonians in **four distinct families**: real anisotropic XYZ,
genuinely complex DM, U1 XXZ, and full exchange tensors with
all nine components nonzero. It compares byte-identical input
and both **full complex H*x bitstreams** (all real and imaginary
IEEE-754 bits) for:

- original padded cache with/without coefficient bound;
- complex CSR retained with/without coefficient bound;
- complex CSR plus padded release with/without coefficient bound;
- small-system independent native full dense-H reference;
- an explicit `MNTE` override and two malformed flag cases.

A separate local GNU14.2 copy of main with this change has:
- `-O2 -fopenmp`: **32/32 model/size/thread cases** passed
  all six modes with *identical full vector bit hashes*,
  N4/N6/N8/N10, T1/T4, all four physical families.
- `-Ofast -fopenmp`: **32/32** passed strict bitwise equality
  between bounded/unbounded within each layout, and strict
  CSR retained versus released bit equality. The historical
  original padded versus CSR layout is compared only under
  explicit `1e-12*max(1,max(|H*x|))` numeric tolerance because
  this compiler option permits reassociation of complex products.
- `make check` and historical CPU-feasible `ex1`–`ex8`:
  PASS under `-O2`.

GNU13/GNU14 CI independently checks the exact PR source with
the same original/CSR/release, real/complex/tensor, `-O2` versus
`-Ofast` numerical policies. No PRESTO node, hardware-exclusive
performance figure, CPC manuscript change or release tag is claimed.

## Long-term engineering option

A true *direct-to-CSR streaming constructor* that never holds
a global padded matrix could further lower the overlap of
the two representations. That would require refactoring
the single physical row builder to preserve all source tensor
terms and their original summation order, with independent
full-H and physics checks. It is a separate, larger change;
do not simulate an unavailable peak-RSS improvement by merely
relabeling this bounded padded allocation.
