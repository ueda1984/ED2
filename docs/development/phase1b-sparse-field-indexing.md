# Phase 1-B: magnetic-field record indexing

Base commit: `233b78551aac7dcef9f106bf1fe2d7ed3491b228`.
Branch: `codex/ed2-phase1-sparse-field-indexing` (the existing branch was fetched
and verified at that base before work). Date: 2026-11-05 (Asia/Tokyo).

Read the [baseline audit](baseline-audit.md), [CPC plan](cpc-library-plan.md) and
[Phase 1-A report](phase1-lapack-interfaces.md) for the surrounding work. This
change does not implement the later Lanczos breakdown/convergence repairs.

## Root cause and exact production change

`input_param` allocates `hvec(3,NO_one)` and reads record `j` into `hvec(:,j)`
while reading its target site into `p_one(j)`. This storage/read convention is
correct and unchanged. An audit of all Fortran uses of `p_one` and `hvec` found
that three Hamiltonian paths incorrectly used the **site** as the field-array
column. Site indices and field-record indices coincide only for a complete,
site-ordered field list; completeness alone is not sufficient.

In `source/ham2vec.f90`, the field reads now use `hvec(1:3,j)` in:

1. `make_wk_loc_and_ele`, the precomputed transition cache;
2. `make_full_hamiltonian`, the dense Hamiltonian;
3. `ham_to_vec_wave_vector`, the uncached/on-demand rows.

There are fifteen corrected component reads in nine assignment statements, plus
three comments distinguishing record `j` from site `i1=p_one(j)`. The site variable
still selects the occupation, local spin/cap, and ladder operation. The positive
`h·S` convention, the `hx ± i*hy` coefficients, the factor 1/2, local ladder factors,
`Sz=S-n`, truncation conditions, basis order, symmetry factors, and the rest of
the Hamiltonian are unchanged. No field-array redistribution or parser change
was introduced.

Lanczos recurrences/breakdown/convergence, LAPACK interfaces, reference data,
existing independent tests, CLI formats and compiler/linker settings are unchanged.
The only build-file additions register optional tests and their Fortran exporter;
Python/NumPy remain test-only dependencies, not mandatory ED2 build/runtime
requirements.

## Failure reproduced before the fix

The pre-fix run used the Phase 1-A bounds-checked executable; its SHA-256 was
verified against the committed Phase 1-A validation record, and its Hamiltonian
source bytes were checked against the branch base before editing. The new
operator-export driver was linked against those same pre-fix objects.

A single field at site 4 with `NO_one=1` exited with a bounds error: field-array
column 4 exceeds upper bound 1. This was reproduced separately in cached ALG=1,
cache-disabled ALG=1 and dense ALG=3, locating the three affected paths. The full
new suite before the fix recorded:

| Group | Passed | Failed |
| --- | ---: | ---: |
| CLI solver cases | 14 | 56 |
| Operator fixtures/thread combinations | 2 | 8 |
| Record-order comparisons | 0 | 28 |

Only the site-ordered full-field control passed all its solver configurations.
The full-field **permuted** dense case returned exit 0 but had maximum spectral
error about `0.1094178576` and matrix error `0.21` against the independent reference.
That incorrect matrix was still Hermitian (measured error zero). Hermiticity alone
therefore cannot detect the field-to-site reassignment. The pre-fix summary,
representative inputs and bounds diagnostics are retained in the evidence folder.

## Independent tests and tolerances

`tests/sparse_field_validation.py` reuses only the independently constructed
spin-1/2 product-basis operator reference `model` from `numerical_baseline.py`;
it does not reuse ED2 Hamiltonian formulas or eigenvalues. Sparse records are
accumulated into reference fields by their explicit site identifiers. The fixed
four-site model has anisotropic XYZ, DM and Gamma bonds and a nondegenerate ground
state, making ground-state magnetization/correlations meaningful. The fixtures include
nonzero x, y and z field components.

Five fixtures cover a single field at site 4; fields at nonconsecutive sites 1/3;
those sparse records reversed; four site-ordered records; and all four records
permuted as sites 3/1/4/2. For each fixture and one/two OpenMP threads:

- ALG=1 and ALG=2 run with `wk_dim=0`, `3`, and `-1` (no, partial, full cache).
- ALG=3 runs the dense spectrum once; it does not use the cache parameter.
- Spectra and all four sites' x/y/z magnetizations and six pairs' nine XYZ
  correlations are compared with independent dense diagonalization. Output index
  columns and array shapes are also checked. Absolute and relative
  tolerances are both `1e-9`. All inputs and stdout are retained in fresh directories.

This gives **70 CLI cases**, plus **28 explicit order-invariance comparisons**
for sparse/full lists under matching solver/cache/thread configurations.

`tests/sparse_field_operator_check.f90` is a test-only exporter linked to the
unchanged production module entry points. It reads the same model, exports ED2's
basis occupations and dense Hamiltonian, then applies H to every unit basis vector
with no/partial/full cache. Python maps the exported occupations to the independent
product basis instead of assuming ED2 basis order. Every element of all four
16x16 matrices is checked against the independent Hamiltonian; Hermiticity and
agreement between dense and applied matrices are checked at absolute tolerance
`1e-12`, relative tolerance zero. **10 operator fixtures/thread combinations**
cover all five models at both thread counts, each checking the three cache modes.
The exporter does not instrument or alter the production executable.

All cases use the identity symmetry sector: nonuniform sparse fields need not
respect translations. Nontrivial symmetry behavior is checked for regressions by
the existing historical workflows, not by imposing an invalid symmetry on these
new fixtures. OpenMP dynamic teams are disabled for the new cases, and nested
BLAS threads are fixed at one.

## Reproduction

Use a fresh external source tree per build profile. Example from this feature
branch after committing, in the prepared cloud environment:

```sh
set -eu
. /workspace/.ed2-env/activate.sh
repo=/workspace/ED2
stage=/workspace/ed2-phase1b-reproduction  # Must not already exist.
mkdir "$stage"
git -C "$repo" archive HEAD | tar -x -C "$stage"
cd "$stage"
autoreconf -vfi
./configure FC=gfortran FCFLAGS='-O2 -fopenmp'
make -j2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
make check-sparse-field SPARSE_FIELD_OUTPUT=/workspace/phase1b-new-results
make check-lapack
make check
make check-long
python3 "$repo/tests/compare_historical.py" --build "$stage" \
  --output /workspace/phase1b-historical.json
```

The bounds profile uses `FCFLAGS='-O2 -g -fopenmp -fcheck=all -fbacktrace -Wall
-Wextra'`. The GNU default profile uses `./configure FC=gfortran` and retains the
existing `-fopenmp -Ofast` defaults. No argument-mismatch workaround is used.
When `SPARSE_FIELD_OUTPUT` is omitted, the optional target creates a new temporary
results directory. Explicit output directories must not exist already.

The unchanged baseline suite is also run against each executable:

```sh
python3 /workspace/ED2/tests/numerical_baseline.py \
  --exe /workspace/ed2-phase1b-reproduction/source/QS3ED2 \
  --output /workspace/phase1b-baseline-results
```

Its exit **1**, 18 passes and three failures (`early-check-lanczos`, `zero-alg1`,
`zero-alg2`) remain intentional observations, not expected-failure annotations or
suppressed failures. Inspect the JSON counts and names; an arbitrary nonzero exit
is not validation.

To reproduce the original failure, archive commit `233b785` instead, configure
that tree with bounds checks and build it. Link the feature branch's exporter
against the old objects, for example:

```sh
gfortran -O2 -g -fopenmp -fcheck=all -fbacktrace -I "$stage/source" \
  "$repo/tests/sparse_field_operator_check.f90" "$stage/source/ham2vec.o" \
  "$stage/source/state_lists.o" "$stage/source/input_param.o" \
  -o /workspace/phase1b-old-operator
python3 "$repo/tests/sparse_field_validation.py" --exe "$stage/source/QS3ED2" \
  --operator-exe /workspace/phase1b-old-operator --output /workspace/phase1b-old-results
```

This last Python command returns 1 with the pre-fix failures listed above. Use a
separate old tree and fresh results path; do not overwrite historical references.

## Results and remaining limitations

The [evidence index](phase1b-sparse-field/README.md) and
[validation record](phase1b-sparse-field/validation.json) identify exact tested
source hashes, binaries, flags, counts and numerical errors. GNU Fortran 14.2.0,
OpenBLAS 0.3.29 pthread/LAPACK LP64 and the machine match Phase 1-A.

Each of the default, normal (`-O2`) and bounds builds passes all 70 CLI cases,
10 operator tests, 28 permutation comparisons, 11 LAPACK checks and the standard
regression test. Each also preserves the 18/3 baseline outcome. Across the new
cases, maximum observed absolute errors are approximately `2.78e-14` (energies),
`1.18e-14` (magnetization), `4.48e-15` (correlations), and `1.11e-16` (operator checks).
These measured errors are not general error bounds or timing claims.

The normal and bounds builds additionally run the full nine historical examples
and the 27 historical numerical comparisons: all nine examples and all 27
comparisons pass in each profile. The default fast profile is not used for
historical accuracy claims.

The three unrelated Lanczos defects remain for Phase 1-C. These new tests cover
four spin-1/2 sites and distinct, valid field sites; they are not exhaustive tests
of mixed spin, duplicate field records, invalid sites, large-system errors, other
compilers or MKL/ILP64. Historical mixed-spin/symmetry examples supplement, rather
than replace, independent small-system checks. Global-state/reentrancy, source
packaging omissions and other baseline issues remain outside this change.
