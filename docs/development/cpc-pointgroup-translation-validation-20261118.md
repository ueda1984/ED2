# Translation and point-group physical verification (2026-11-18)

**CI-verified numerical regression:** the 31-case standard physical suite passed on both gfortran-13 and gfortran-14 at [GitHub Actions run 44215282646](https://github.com/ueda1984/ED2/actions/runs/44215282646) (2026-11-18). This is not a claim that the optional 6x6 three-magnon stress tier or exclusive-node PRESTO performance tests have passed. The tracked script is `tests/cpc_pointgroup_translation_crosscheck.py`, invoked using `make check-pointgroup`. As a research-only regression, it requires Python/NumPy in addition to the normal ED2 build and is not run by the default `make check`.

## Motivation and provenance

The current ED2 source accepts up to six site permutations `FILE_S1..FILE_S6` with cyclic periods `L1..L6` and character indices `M1..M6`. This already existed in the QS-Cube-derived implementation. The optimized `try_orbit_basis` only accepts mutually commuting generators; for a noncommuting point-group/translation semidirect product it conservatively returns to the older `checkstate_SQ` basis builder. The test deliberately checks that behavior and **does not** advertise a new non-Abelian irreducible-representation engine.

The original QS³/ED and QS-Cube/ED2 were capable of some permutation-based sector projections; modern claims must be limited to accurately checked pre-existing capability and newly added **independent numerical validation**.

## Physical models and group

The finite periodic spin-1/2 models are
```
H = sum_<ij> Jxy_ij (Sx_i Sx_j + Sy_i Sy_j) + Jz_ij Sz_i Sz_j
```
with nearest-neighbor `(Jxy,Jz)=(1,0.77)` and next-neighbor `(0.31,0.19)`. The ring uses distance-two neighbors; square lattices use both diagonals. The anisotropic but spin-U(1)-preserving exchange checks nontrivial off-diagonal terms and physically exact fixed-magnon sectors.

- Periodic 12-site ring: translations `T` and inversion `I` at `k=0,pi`, even/odd parity.
- Periodic 4x4 square lattice: `T_x,T_y,C4,mirror` at Gamma and M. The one-dimensional `C4v` irreps `A1,A2,B1,B2` are represented by characters `C4 = +/-1`, `mirror = +/-1`.
- Periodic 4x4 square lattice at X=(pi,0): `T_x,T_y,C2,mirror`. Its little group is `C2v`, **not** the full `C4v`, as 90-degree rotation exchanges X and Y. The code rejects unsupported generic momentum point-group combinations.
- Periodic 6x6 square lattice: two-magnon Gamma sectors and, on demand, three-magnon stress cases.

For Gamma/M with full `C4v`, the two-dimensional `E` irrep is **not** represented by the current scalar-character input. Its absence is measured as missing dimension; the test does not claim that all symmetry sectors are included.

## Independent reference / acceptance criteria

1. Explicitly enumerate bitmask spin states and every group element, checking the group elements are distinct and each generator permutes an invariant bond set. This is an independent guard, not a production ED2 input validator.
2. Compute `dim = sum_g chi(g)^* Tr(U_g)/|G|` from full Hilbert-space permutations; ensure an integer dimension. Explicitly construct character-projected orbit vectors and compare the count with the character-trace dimension. Verify orthonormality and each generator's eigenvalue.
3. Apply an independent XXZ Hamiltonian on the full fixed-magnon bitmask basis, project by the independently constructed vectors, and verify Hermiticity and **no leakage** out of the irrep: `||(1-P) H Q|| < 2e-8`.
4. Compare **all** eigenvalues from native ED2 dense diagonalization with the independently projected matrix to absolute `5e-9`. Check ring parity and X-point `C2v` spectrum decomposition. For Gamma/M track the unrepresented `E` sector rather than asserting false completeness.
5. For selected noncommuting point-group sectors run **both** modern ALG=1 modes against the independently established lowest eigenvalue and check direct-residual telemetry. This makes the exercise relevant to the revived historical two-stage Lanczos PR.
6. Record compiled-binary SHA256, exact per-case input hashes, projected levels, and the employed compiler in archived CI artifacts.

## Verified standard-suite sector dimensions (CI 44215282646)

| System / momentum | Translation-only dim | 1D irrep dimensions (A1,A2,B1,B2 for square) | Unrepresented dimension |
| --- | ---: | --- | ---: |
| Ring 12, 3-down, Gamma (even/odd) | 19 | 12, 7 | 0 |
| Ring 12, 3-down, pi (even/odd) | 19 | 12, 7 | 0 |
| Square 4x4, 3-down, Gamma (C4v) | 35 | 10, 1, 6, 4 | 14 (E irrep) |
| Square 4x4, 3-down, M (C4v) | 35 | 10, 1, 6, 4 | 14 (E irrep) |
| Square 4x4, 3-down, X (C2v) | 35 | 16, 5, 7, 7 | 0 |
| Square 4x4, 2-down, Gamma (C4v) | 9 | 5, 0, 3, 1 | 0 |
| Square 6x6, 2-down, Gamma (C4v) | 19 | 9, 1, 6, 3 | 0 |

The empty square-4x4 two-down Gamma A2 sector is counted as a mathematically valid *empty sector*; no ED2 eigenvalue run is attempted for it. The remaining nonempty sectors were checked against native ED2 dense spectra. The standard suite additionally runs both ALG=1 Lanczos modes for selected noncommuting sectors.

## Source-derived caveats

- Supplied maps and phases are generally **not** independently checked for Hamiltonian invariance or legitimate little-group characters by the ED2 production reader. Incorrect input could yield a formally projected but physically invalid spectrum. The test suite supplies a separate oracle for its specific fixtures.
- Because `C4v` is non-Abelian, one cannot infer a general non-Abelian implementation from its four one-dimensional representations.
- ED2's optimized atomic-orbit construction currently does not handle noncommuting generators directly. Benchmarking speedup and solver correctness must not be conflated.
- No complete real-world performance claim is permitted from ordinary GitHub hosted CI nodes. A separate controlled PRESTO measurement must compare translation-only versus translation+point-group sector costs under pinned workloads, thread counts, cache settings, and exclusive-node attestation.
- Full physical certification additionally needs multiple momentum/sector controls, an independent lowest-state-inclusion argument, and review of all null-sector and degeneracy edge cases.

## Reproduction

```sh
autoreconf -vfi
./configure FC=gfortran-13 FCFLAGS="-O2 -fopenmp"
make -j2
make check
make check-pointgroup POINTGROUP_OUTPUT=/tmp/ed2-pointgroup
python3 tests/cpc_pointgroup_translation_crosscheck.py --exe ./source/QS3ED2 --stress
```

The second invocation chooses an automatic temporary output directory. The script writes `report.json`, each case's complete input/STDOUT, and independently computed eigenvalues. No modification to physics production source is needed for this verification step.

## CPC presentation

The defensible methods statement is: *ED2's pre-existing site-permutation symmetry interface is exercised and validated for combined translation and compatible one-dimensional point-group characters by independent projection, group-theoretic sector counting, and fully resolved spectra.* This is a **test coverage and verification contribution**, not invention of point-group reduction or low-memory Lanczos.
