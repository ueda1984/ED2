# Release qualification and large-example validation

ED2 development `main` and Read the Docs `latest` are moving targets. A
successful merge or documentation build does **not** constitute a tagged
software release. `v1.0.0` is the approved first release **candidate**, not
an issued release. Until the tag/Release are published, cite the full Git
commit used for each calculation.

## Before publishing the approved `v1.0.0` candidate

1. Review a candidate branch and ensure the source tree is clean. Record its
   SHA, compiler/flags, BLAS/LAPACK implementation, operating system, and
   `ED2_*` environment options.
2. Check GNU Fortran 13 and 14 CI, independent small-system physics checks,
   bounds-checked runtime tests, and `make distcheck`. Keep logs and artifacts.
3. Run a fresh, resource-controlled `ex9` calculation, **including both
   two-body correlation files**, or document that the expensive example is not
   yet independently reproduced. `make check` and the regular CI do not
   execute this complete 1000-site workload.
4. Review solver residuals and numerical limitations. Treat GNU 13/14
   results and any Intel oneAPI ifx/LP64 oneMKL run as **separate recorded
   configurations**, with compiler versions, link libraries, and audit
   artifacts. An Intel compiler build alone is not a validation of all
   scientific input families; a 1000-site ifx ex9 run is an independent
   optional benchmark.
5. The initial candidate version is approved as **`v1.0.0`**; `AC_INIT` already
   says `1.0.0`. Set `version: 1.0.0` in `CITATION.cff` for candidate CI,
   but **do not assert a date of publication before publication is scheduled**.
   In this state omit `date-released`. Review release notes and all limitations.
6. First verify the candidate with `python3 tests/check_release_metadata.py`
   (mode `release_candidate`), GNU/Intel physics, source-distribution and
   archived ex9 checks. Only after approving the *actual* publication date,
   add `date-released: YYYY-MM-DD` to `CITATION.cff`, commit that change, and
   freeze the **new final SHA**. Run
   `python3 tests/check_release_metadata.py --release-tag v1.0.0` and
   `tests/release_archive_preflight.py --archive <tar.gz>
   --expected-commit <full-final-SHA> --release-tag v1.0.0`. Rerun all
   GNU/Intel physics and source-distribution checks on the exact **final**
   commit, retain the archive and SHA-256 manifest, then create the tag.
7. With the appropriate repository permissions, create and publish the
   reviewed tag/release. Confirm that the tag points to the tested commit.
   Build the versioned Read the Docs pages and archive the matching source and
   checksums. Add an actual archival DOI only after issuance; do not cite a
   placeholder Zenodo or journal DOI.

The default `python3 tests/check_release_metadata.py` accepts (1) development
with both CFF fields absent, (2) `release_candidate` with `version` present
and `date-released` absent, and (3) complete versioned metadata with both
fields present and consistent with `AC_INIT`. The `--release-tag` mode
additionally requires a valid actual release date and exact tag version.
A date without a version, version mismatches, malformed dates and missing
date on tag builds all fail. Passing a check never creates a release or DOI.

The archive preflight requires every tracked release input: numerical
Fortran modules, test scripts/drivers, example generators and inputs,
historical reference tables, documentation and build templates. It compares
those bytes with the frozen Git commit. Regenerated Autotools files are
listed and checked separately against the build tree, because they depend
on the toolchain. Missing, changed or unclassified archive files fail. These
checks complement a successful unpacked `make distcheck`; neither replaces
scientific tests. Run `python3 -B tests/test_release_guards.py` for the negative
cases covering omitted/tampered solver sources and metadata transitions.

## How to complete the 1000-site ex9 verification

The `script_cubic_sp_HB/run.sh` example uses a periodic \(10\times10\times10\)
spin-1/2 Heisenberg lattice, `NODmin=NODmax=3`, three translations and three
reflections, with 9 requested site pairs. The unreduced sector contains
166,167,000 spin configurations; the example's saved reference reports
23,719 representatives in the requested symmetry sector. The expensive
pair-correlation stage is **not** covered by lightweight CI.

Run it in a **fresh source checkout or clean scratch copy** on a compute node,
not on a directory containing valuable past output. Its current example
script clears `eigenvalues.dat` but observable files are opened for append,
so preexisting outputs can invalidate comparisons. From the source root:

```bash
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

git rev-parse HEAD
git status --short
sh tests/run_ex9_reproducibility.sh

# If you already have a completed full run, compare it independently:
python3 tests/verify_ex9.py \
  --actual examples/cubic_sp_HB \
  --legacy-ex9-symmetry-projection \
  --json-out ex9-comparison.json
```

The compute is CPU/memory intensive. The wrapper does not submit a job or
request scheduler resources. It refuses pre-existing outputs to prevent
append contamination, then preserves complete inputs, output files, checksums,
and a run log under a new `ex9-audit-*` directory. It records the configured
compiler/flags/libraries from the generated Makefile, `config.log` when
available, executable hashes and dynamic linkage. A different compiler on
the shell's default PATH does not identify the compiler that built ED2.
Select the node class,
thread count, and job walltime under the local scheduler's policy; retain scheduler peak RSS,
walltime, stdout/stderr, complete input files, and the resulting
`eigenvalues.dat`, `local_mag.dat`, `two_body_cf_xyz.dat`, and
`two_body_cf_z+-.dat`. A killed, partial, or appended observable file is
**not** a successful result.

The comparison verifies all of the following:

- exactly 1 finite energy, 1000 complete local-magnetization rows, and 9
  rows of each two-body correlation format;
- ordered site/pair identifiers, spin-1/2 operator bounds, vanishing
  transverse one-site magnetizations, and total
  \(\sum_i\langle S_i^z\rangle=1000/2-3=497\);
- agreement between the Cartesian and ladder correlator output formats;
- absolute differences from archived reference values (energy tolerance
  `1e-7`, observable tolerance `1e-8`). The legacy ladder +-/-+ entries
  require the explicitly enabled real/even-reflection projection, while
  both the raw and symmetry-projected differences remain reported.

Use `python3 tests/verify_ex9.py` without `--actual` to check **only the
archived reference fixture**. That fast check is suitable for CI, but it
must **not** be described as a new 1000-site calculation or as independent
reproduction of the archived correlation functions. Even a full comparison
to archived outputs needs solver convergence and model validation to support
scientific claims.

If the full calculation is too slow, record the timed-out stage separately.
Optimize or parallelize the pair-correlation measurement only with
corresponding small-system independent-operator tests; do not change the
scientific benchmark or skip `cal_cf=1` merely to obtain a green check.

## Completed two-environment ex9 audit (19 November 2026)

Full ex9 calculations were completed with unchanged numerical source from
official main commit `63f81c198df6adfa992932e1b7c7fa5c359c8c43`.
The GitHub Actions audit checkout `035d81fa225e2cf362199f3e4b8d1a10a3da4ef2`
differs from main **only** by the audit workflow. Both executions used
GNU Fortran 14, `-O2 -fopenmp`, two OpenMP threads, one BLAS thread and
the supplied 1000-site/three-magnon input:

- Ubuntu 24.04, LP64 OpenBLAS/LAPACK:
  [full ex9 evidence run #44514651378](https://github.com/ueda1984/ED2/actions/runs/44514651378),
  6:50.75 wall-clock elapsed, 140 Lanczos steps;
- Debian/GNU Fortran 14.2, system LP64 BLAS/LAPACK:
  9:45.04 wall-clock elapsed; independently rebuilt from the same archived
  source snapshot.

The Ubuntu build job explicitly used `FC=gfortran-14` and linked
`-llapack -lopenblas`. Its older `environment.txt` collector printed the
default `gfortran` (13.3), which was not the selected build compiler. Use
the actual compile/link log for that audit. Generated build files in the
separately exported source snapshot are historical checkout files; the
Debian reproduction regenerated and reconfigured them before building.

The first run produced all four numerical output files but its historical
**raw** ladder comparison failed. This is a genuine archived-fixture
inconsistency, **not** an incomplete ex9 calculation. The previous
`calcu_cf` averaged only translation generators 1--3, omitting reflection
generators 4--6. The historical `two_body_cf_z+-.dat` consequently violates
the operator-adjoint identity
`<S^-_i S^+_j> = conj(<S^+_i S^-_j>)`
by as much as `2.7829938541e-5`. For the real, zero-momentum and
even-reflection ex9 state, these values must moreover be equal and real.

The corrected implementation produces Hermitian, real, pair-exchange
symmetric ladder correlations. The historical two channels are projected
**for comparison only** to their common real arithmetic mean, with
imaginary parts set to zero; all historical archive files remain untouched.
The explicit CLI flag
`--legacy-ex9-symmetry-projection` enables this *ex9-specific* projection,
and always records the unsymmetrized raw discrepancy.

This average is justified for the archived collinear pairs `(1,2)` through
`(1,10)` of this particular isotropic, zero-field, zero-momentum,
even-reflection ex9 input. Reflection along that line exchanges the two
oriented ladder channels; the transverse reflections preserve the pair
displacement. It is not a general correction for arbitrary real or
reflection-symmetric models: independent small-cube checks exhibit
off-axis pairs for which the two-channel mean is wrong. Projection mode
therefore requires the exact four historical numerical reference tables
by SHA-256; ordinary strict comparisons remain available for other data.
The computed output is never symmetrized or modified by the verifier.

| Full output comparison | Maximum absolute discrepancy |
| --- | ---: |
| Ubuntu vs Debian: energy | `6.82e-13` |
| Ubuntu vs Debian: all 1000 local magnetizations | `2.07e-12` |
| Ubuntu vs Debian: Cartesian correlations (9 pairs) | `4.66e-12` |
| Ubuntu vs Debian: ladder correlations (9 pairs) | `9.32e-12` |
| Ubuntu vs raw historical ladder archive | `1.39150e-5` |
| Ubuntu vs **symmetry-projected** historical ladder | `1.89112e-11` |
| Debian vs symmetry-projected historical ladder | `9.59420e-12` |

The matched energy is `736.7224213680655` (sector-specific);
`sum_i <S^z_i>=497`. The new runs are independent builds/environments of
the same code, **not** independent algorithms or independent verification
of every Hamiltonian term. They resolve the ex9 benchmark reproduction
question with the historical reference discrepancy explicitly documented.

The former `examples/ref_data/cubic_sp_HB` directory remains a historical
fixture. Do not overwrite or silently loosen comparisons to make it
appear to match. A future immutable release should separately archive
the newly computed outputs, build logs and source checksums.
