# Native ED2 cached H*x scratch allocation elision

A general, non-cubic source-level optimization that eliminates
two unused per-application allocations in the ED2
`ham_to_vec_wave_vector` routine when `wk_dim == lv`
(all Hamiltonian rows are already cached). It is independent
of the type of magnetic interactions, real/complex DM,
point-group symmetry and the CSR storage choices.

## Source change

Previously each H*x call unconditionally allocated the
`ni(:)` and `st_list(:)` temporary state-transition buffers
before computing cached rows, then considered the on-the-fly
loop over `i=wk_dim+1,...,lv`. When `wk_dim==lv`, this
second loop has no iterations and **never reads the buffers**.

Allocate only when `wk_dim<lv`:

    if(wk_dim<lv)allocate(ni(NODmax),st_list(NODmax))

The old Fortran/OpenMP `private(ni,st_list)` behavior
remains unchanged in **partially cached** cases, because
both buffers are allocated as before. No mathematical
term, coefficient, summation order or default flag is changed.

The optimization applies to the historical padded cache,
the complex CSR cache and the optional CSR padded-cache
release introduced by PR29. It avoids unnecessary repeated
memory-management operations for **any** full cache.

## Independent source-to-source correctness

A separate original-main local GNU14.2 benchmark built both
the prechange and modified `ham2vec` from otherwise identical
sources with `-O2 -fopenmp`, using the same native complete
complex H*x proof driver and byte-identical input files.
There were **48/48 fully bitwise identical pairs**:
- N4/N6/N8/N10 full alternating S1/2,S1 OBC
  non-U1 XYZ/transverse-field models,
- both strictly real and actually nonreal DM,
- T1/T4,
- three cache layouts per physical case: full padded,
  fully packed complex CSR, and partial `wk_dim=1`
  padded/on-the-fly Hamiltonian.

For each input the two native drivers used identical
complex input vectors, produced output arrays of identical
dimensions and **identical SHA256 of every complex vector
output bit**. The partially cached variant demonstrates
that the scratch is still correctly allocated and used
when necessary. No PRESTO measurements were used.

The checked-in `tests/cpc_partial_cache_allocation_regression.py`
independently verifies that the *modified* main full and
partially cached operators agree numerically on physical
inputs, with the separately constructed dense Hamiltonian
as another reference for N4/N6. GNU Fortran 13/14
CI retains the previous exact three-mode CSR tests,
BLAS/LAPACK ABI proof and historical ex1–ex8 examples.

## Performance claim boundary

This is a proven elimination of two allocations at each
fully cached H*x call, **not** a measured or
publication-quality end-to-end speedup. The actual
effect depends on system size, threading, allocator,
cache representation and the total H-action count.
The improvement does not change memory at the peak
operator construction stage. No default flag,
public symbol or CPC manuscript result is changed.
