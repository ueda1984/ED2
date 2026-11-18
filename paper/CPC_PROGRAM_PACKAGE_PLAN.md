# CPC Program Library source-and-test package: concrete assembly plan

**Status: proposed structure only, not an official release or CPC submission.**
The CPC target is a CPiP article for the `QS3ED2` Fortran executable,
not an installed reusable Fortran library. The primary published predecessor
is QS^3/ED (CPC 277, 108369, 2022); the intended official generalized ED2
source is QS-Cube/ED2, subject to upstream/project approval.

## Final archive contents (freeze only after release approval)

| Directory | Contents | Acceptance |
| --- | --- | --- |
| `source/`, `Makefile.am`, `configure.ac` | Reviewed original program sources and Autotools build recipe | A clean extracted archive builds without untracked files |
| `LICENSE`, `CITATION.cff`, `README.md` | MIT notice, original contributor attribution, version and citation metadata | Version, date and tag match the approved release |
| `tests/`, `examples/` | Standard reference examples, regression and scientific verification scripts | `make check`, `make check-long`, and packaged example tests pass |
| `program_examples/` | Three **reviewed** sample input/output sets with expected diagnostics | Every input hash, energy, residual and output checksum documented |
| `docs/`, `paper/` | Model conventions, executable usage, reproducibility appendix and manuscript | README instructions tested from unpacked archive |
| `provenance/` | Compiler, numerical libraries, source/binary SHA256, relevant CI and independent reference records | Exact build/release provenance and separate numerical evidence |
| `performance_pilot/` (optional) | Full read-only raw PRESTO input/log/trace/audit records, subject to redistribution clearance | Both original jobs and operator reservation evidence present; no partial-only favorable trial selection |

**Do not manufacture a software DOI, CPC accession, release date or
author list.** The approved CPiP/Mendeley Data deposit workflow may supply
a permanent Program Library record; a separate Zenodo DOI, if chosen,
requires its own metadata and version policy.

## Three practical example families for the program package

### Example A: basic conserved spin model / regression continuity

From a clean archive on Linux:

    autoreconf -vfi
    ./configure FC=gfortran FCFLAGS="-O2 -fopenmp"
    make -j2
    make check
    make check-long

This checks the existing nine native example families and the basic
build/legacy input contract. Select at least one historical
XXZ test input and its *actual native expected output* from the
clean run for the final `program_examples/01-xxz/` folder.
Do **not** claim independent original QS3/QS-Cube/ED comparison
unless the separate old executable and source hash are included.

### Example B: non-U(1) alternating mixed-spin oriented-DM response

After building the executable, no Python scientific stack is needed:

    export OPENBLAS_NUM_THREADS=1
    python3 -S tests/cpc_portable_dm_quickstart.py       --exe "$(pwd)/source/QS3ED2" --size 4 --omp 1       --output "$(pwd)/cpc-example-dm-l4"

This verifies a fully specified spin-1/2/spin-1 Hamiltonian
at neighboring DM strengths, native vector-chirality output and
an independent finite-size reference/response identity.
Select actual source-generated `input.dat`, native outputs,
diagnostics and independent expected checks into
`program_examples/02-mixed-spin-dm/` after final freeze.

### Example C: literature Hamiltonian convention and observable test

Requires NumPy/SciPy for the *independent reference program*
(the ED2 executable itself only needs Fortran BLAS/LAPACK):

    python3 tests/cpc_literature_yi_xyz_dm.py       --exe "$(pwd)/source/QS3ED2"       --sizes 4,6,8 --omp 1       --output "$(pwd)/cpc-example-xyz-dm"

The reference validates the finite-size Hamiltonian from
Yi et al., Phys. Rev. B 100, 024423 (2019), including
Pauli-versus-spin conventions, oriented periodic DM bonds,
spectra, local moments and Cartesian correlators.
Package source-generated cases plus reference outputs under
`program_examples/03-xyz-dm-literature/`.

### Advanced optional example: 1000-site cubic Gamma/O_h

PR17's checked N1000/Ndown3 1D-irrep benchmark is a useful
large-site demonstration, but it has dedicated HPC resource
requirements and is **not** an appropriate sole beginner
test. The 12-row medians and scope are in
`paper/data/pr17_presto_coordinate_cache_pilot.csv` and
`paper/data/README.md`. Formal article-level timing use
still requires allocation attestation, raw trial preservation
and further controlled repetitions. The CPiP package can
include this as an optional workload once independent
archiving rights and final source are established.

## Release acceptance matrix

Before labeling the package ready:

- [ ] Coauthors/maintainers choose the official QS-Cube/ED2 release SHA,
      tag, date, Program Title, version and CRediT metadata.
- [ ] `configure.ac` and `CITATION.cff` match; the current
      **1.0.0 vs 0.1.0** mismatch is not silently changed.
- [ ] Source archive is recreated reproducibly, hashed and independently
      extracted/compiled on at least two supported compilers.
- [ ] Numerical regression cases and all three example families pass
      with archived inputs, diagnostics and reference outputs.
- [ ] Program Summary and README list actual runtime dependencies and
      restrictions, including nonconserving occupation-window
      projections and no installed library ABI.
- [ ] Author-approved software archive is deposited through the
      journal's required program-library route; DOI/accession
      is copied only from the real assigned record.
- [ ] Final article, figures, declarations and cover letter have
      unanimous sign-off. No submission before that point.

The existing source-pinned CI demonstrations establish much of the
physics case already, but the final program library bundle must
refer to an immutable **approved** release rather than a moving
development branch. Refer to `paper/CPC_SUBMISSION_GATES.md`.
