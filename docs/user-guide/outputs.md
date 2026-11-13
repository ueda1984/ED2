# Output files and column conventions

ED2 writes eigenvalues, observables, and optional wavefunctions to the locations
specified by `OUTDIR` and `FILE_wf` in `&input_parameters`. These files do
**not** all use the same text format.

!!! warning "Start with an empty output directory"
    Several output routines open result files with `position='append'`.
    Re-running an example in the same directory can accumulate stale rows.
    Use a fresh run directory or intentionally remove old results first.

## Eigenvalues: `eigenvalues.dat`

A new diagonalization writes a text file under `OUTDIR`:

| Column | Meaning |
| --- | --- |
| 1 | Eigenstate index (1-based) |
| 2 | Energy in the units of the Hamiltonian coefficients |

ALG=1 normally writes one low-energy eigenvalue; ALG=2 writes its requested
eigenvalues; ALG=3 writes the full dense spectrum. The first column is an
**index**, not the energy. A generated file does not, by itself, prove
ground-state inclusion or sufficient observable accuracy.

## Local magnetization: `local_mag.dat`

When `cal_lm = 1`, the current main executable writes values for one
selected wavefunction:

| Column | Quantity |
| --- | --- |
| 1 | Site `i` |
| 2 | `<S_i^x>` |
| 3 | `<S_i^y>` |
| 4 | `<S_i^z>` |

## Two-body Cartesian correlations: `two_body_cf_xyz.dat`

When `cal_cf = 1`, each row contains the original requested site labels
and nine **real-valued** correlations. The order is **not** conventional
row-major `xx, xy, xz, yx, ...` order. The source
(`source/get_expectation_values.f90`) uses:

| Column | Quantity |
| --- | --- |
| 1 | Site `i` |
| 2 | Site `j` |
| 3 | `<S_i^x S_j^x>` |
| 4 | `<S_i^y S_j^x>` |
| 5 | `<S_i^z S_j^x>` |
| 6 | `<S_i^x S_j^y>` |
| 7 | `<S_i^y S_j^y>` |
| 8 | `<S_i^z S_j^y>` |
| 9 | `<S_i^x S_j^z>` |
| 10 | `<S_i^y S_j^z>` |
| 11 | `<S_i^z S_j^z>` |

The program handles reversed requested pairs using its internal swap list;
retain the file's displayed site order when interpreting the results.

## Ladder-operator correlations: `two_body_cf_z+-.dat`

This companion file is also written when `cal_cf = 1`. It contains
site labels followed by **complex** correlations in the
`z,+,-` operator basis. Its complex text representation and
component order are distinct from those of `two_body_cf_xyz.dat`.
Check the source's `expe_cf` ordering before using this file in
post-processing.

## Optional wavefunctions and energies

With `wr_wf = 1`, ED2 writes wavefunctions to `FILE_wf` using names
such as `wf1.bin`, `wf2.bin`, etc., plus an `energy.dat` text
file. The actual file format is **Fortran sequential unformatted
binary** with complex amplitudes in the constructed **symmetry-sector
basis**. It is **not** a readable `eigenvec_XXXX.dat` text table.

- `energy.dat` stores the state index and energy and supports
  the corresponding readback workflow (`re_wf = 1`).
- Compiler-specific record markers and complex representation mean
  that a binary file may not be portable across compilers.
- The current legacy filename formatter uses an `I1` integer field;
  more than nine indexed wavefunction filenames are **not validated
  or generally supported** by this output path.
- The `FILE_wf` string is concatenated with the output filename.
  Supply a path with the appropriate directory separator.

## Logs and optional solver diagnostics

Stdout can be redirected with `> output.dat` from the directory
containing a valid `input.dat` and its referenced input files.

In the **Phase 1-F/1-G development branch**, `ED2_DIAGNOSTICS`
optionally creates a separate text sidecar with solver status,
residuals, iterations, and provenance. It is not written by default.
See [Solver accuracy and reproducibility](solver-accuracy.md).

Archive the full source SHA, compiler/BLAS versions, input and list
files, thread counts, seed or explicit initial vector, and diagnostic
sidecar (if enabled) for reproducible numerical work.
