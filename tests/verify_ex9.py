#!/usr/bin/env python3
"""Validate archived or freshly computed ED2 N=1000 cubic three-magnon outputs.

This checker does not run ex9. An archive-only pass is NOT a fresh reproduction.
Uses Python's standard library; fails on partial/appended/NaN/Inf output.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

TABLES = {
    "energy": ("eigenvalues.dat", 1, 2),
    "mag": ("local_mag.dat", 1000, 4),
    "xyz": ("two_body_cf_xyz.dat", 9, 11),
    "zpm": ("two_body_cf_z+-.dat", 9, 20),
}

# This correction is justified for one known archive, not for arbitrary
# correlation data. Pair IDs alone cannot establish the lattice geometry or
# the symmetry sector. Keep byte-exact copies of all four original tables;
# changing even whitespace deliberately requires a new provenance review.
LEGACY_EX9_REFERENCE_SHA256 = {
    "eigenvalues.dat": "0f3824d777fd4e21336d6a48bb3a3c84c0be8a4afdc05f72710f6741da17c6bf",
    "local_mag.dat": "14f6aa1bded41068080bb42b3683d27989a580fadc563455bc61df530ddf61f3",
    "two_body_cf_xyz.dat": "a0ae1b668bace46ea22ee292c50ce20a6dfc14c0dfac8f1058cc1b27814f6a92",
    "two_body_cf_z+-.dat": "d764764d76b21944bd6ddd3632e3fc94ec3d5f2032ca8eb295afe213b219e03b",
}


def check_legacy_reference_identity(folder):
    for name, expected in LEGACY_EX9_REFERENCE_SHA256.items():
        observed = hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if observed != expected:
            raise ValueError(
                f"{name}: not the pinned historical ex9 reference; "
                "symmetry projection requires byte-exact original tables")
    return "SHA256 matches all four original historical ex9 tables"


def read_table(folder, file_name, count, width):
    path = folder / file_name
    rows = []
    with path.open(encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != width:
                raise ValueError(f"{path}:{line_no}: expected {width} columns, got {len(fields)}")
            try:
                row = [float(v.replace("D", "E").replace("d", "e")) for v in fields]
            except ValueError as error:
                raise ValueError(f"{path}:{line_no}: invalid numeric value") from error
            if not all(math.isfinite(v) for v in row):
                raise ValueError(f"{path}:{line_no}: NaN or Inf")
            rows.append(row)
    if len(rows) != count:
        raise ValueError(f"{path}: expected exactly {count} rows, got {len(rows)} (partial/append?)")
    return rows


def read_results(folder):
    return {key: read_table(folder, *spec) for key, spec in TABLES.items()}


def cartesian_from_ladder(row):
    # Fortran's nine column-major complex <A_i B_j>, A,B ordered (z,-,+).
    z = [complex(row[2 + 2*k], row[3 + 2*k]) for k in range(9)]
    return [
        (z[8] + z[7] + z[5] + z[4]).real * .25,        # xx
        ((z[8] - z[7] + z[5] - z[4]) * (-.25j)).real,  # yx
        ((z[6] + z[3]) * .5).real,                      # zx
        ((z[8] + z[7] - z[5] - z[4]) * (-.25j)).real,  # xy
        ((z[8] - z[7] - z[5] + z[4]) * (-.25)).real,   # yy
        ((z[6] - z[3]) * (-.5j)).real,                 # zy
        ((z[2] + z[1]) * .5).real,                      # xz
        ((z[2] - z[1]) * (-.5j)).real,                 # yz
        z[0].real,                                     # zz
    ]


def zpm_hermiticity_gap(rows):
    # For distinct sites, (S^-_i S^+_j)^dag = S^+_i S^-_j.
    return max(abs(complex(row[12], row[13]) -
                   complex(row[16], -row[17])) for row in rows)


def symmetrize_historical_zpm(row):
    # EX9 ONLY: site-exchange symmetry and operator adjoints force equal real
    # +- and -+ channels in this k=0, even-reflection sector. The old archive
    # included translations but omitted reflections in the correlation average.
    # For its nine collinear pairs, full reflection averaging exchanges these
    # two oriented channels. This mean is not valid for general off-axis pairs.
    # Preserve the unmodified archive and project only these channels.
    result = list(row)
    symmetric = 0.5 * (row[12] + row[16])
    result[12] = symmetric
    result[16] = symmetric
    result[13] = 0.0
    result[17] = 0.0
    return result


def check_physics(data, tolerance=2e-7):
    energy, mag, xyz, zpm = (data[k] for k in ("energy", "mag", "xyz", "zpm"))
    if energy[0][0] != 1:
        raise ValueError("Expected one ground-state energy with state index 1")
    for i, row in enumerate(mag, 1):
        if row[0] != i:
            raise ValueError(f"Local-magnetization site index error at {i}")
        if any(abs(x) > .5 + tolerance for x in row[1:]):
            raise ValueError(f"Spin-1/2 local-magnetization bound violated at site {i}")
        if abs(row[1]) > tolerance or abs(row[2]) > tolerance:
            raise ValueError(f"Fixed Sz sector has transverse magnetization at {i}")
    sz = math.fsum(row[3] for row in mag)
    if abs(sz - 497.0) > 1000.0 * tolerance:
        raise ValueError(f"Total Sz must be 1000/2-3=497, got {sz}")
    largest_identity_error = 0.0
    for i, (car, lad) in enumerate(zip(xyz, zpm), 2):
        if (car[:2] != [1.0, float(i)] or lad[:2] != [1.0, float(i)]):
            raise ValueError(f"Pair {i} must be (1,{i}) in both outputs")
        if any(abs(v) > .25 + tolerance for v in car[2:]):
            raise ValueError(f"Spin-1/2 Cartesian correlator bound violated at pair {i}")
        predicted = cartesian_from_ladder(lad)
        largest_identity_error = max(
            largest_identity_error,
            max(abs(x - y) for x, y in zip(car[2:], predicted)))
    if largest_identity_error > tolerance:
        raise ValueError(f"Ladder/Cartesian consistency discrepancy {largest_identity_error}")
    return {"ground_energy": energy[0][1], "total_Sz": sz,
            "max_ladder_cartesian_error": largest_identity_error,
            "max_ladder_hermiticity_gap": zpm_hermiticity_gap(zpm)}


def compare(actual, reference, energy_atol, observable_atol,
            legacy_ex9_projection=False):
    differences = {}
    for name in TABLES:
        start = 1 if name in ("energy", "mag") else 2
        raw_error = max(abs(a[k] - b[k])
                        for a, b in zip(actual[name], reference[name])
                        for k in range(start, len(a)))
        comparison_ref = (list(map(symmetrize_historical_zpm, reference[name]))
                          if name == "zpm" and legacy_ex9_projection
                          else reference[name])
        error = max(abs(a[k] - b[k])
                    for a, b in zip(actual[name], comparison_ref)
                    for k in range(start, len(a)))
        limit = energy_atol if name == "energy" else observable_atol
        differences[name] = {"max_abs_error": error, "raw_max_abs_error": raw_error,
                             "atol": limit,
                             "historical_symmetry_projection":
                             name == "zpm" and legacy_ex9_projection}
        if error > limit:
            raise ValueError(f"{name}: computed/reference difference {error} exceeds {limit}"
                             + (" after historical symmetry projection"
                                if name == "zpm" and legacy_ex9_projection else ""))
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path,
                        default=Path(__file__).resolve().parents[1] /
                        "examples/ref_data/cubic_sp_HB")
    parser.add_argument("--actual", type=Path,
                        help="Completed ex9 output directory; omit for archive-only validation")
    parser.add_argument("--energy-atol", type=float, default=1e-7)
    parser.add_argument("--observable-atol", type=float, default=1e-8)
    parser.add_argument("--legacy-ex9-symmetry-projection", action="store_true",
                        help="Compare to the historical ladder archive after the "
                             "exact real even-reflection projection; requires byte-exact "
                             "original reference tables and reports raw differences")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    report = {"mode": "archive_only", "reference": str(args.reference.resolve())}
    try:
        if not (math.isfinite(args.energy_atol) and args.energy_atol > 0
                and math.isfinite(args.observable_atol) and args.observable_atol > 0):
            raise ValueError("Positive, finite tolerances required")
        if args.legacy_ex9_symmetry_projection:
            report["legacy_reference_identity"] = check_legacy_reference_identity(args.reference)
        ref = read_results(args.reference)
        report["reference_checks"] = check_physics(ref)
        if report["reference_checks"]["max_ladder_hermiticity_gap"] > 1e-8:
            report["historical_archive_warning"] = (
                "Legacy ex9 ladder data violate operator Hermiticity; "
                "archive-only PASS checks structural data, not this physical identity. "
                "Fresh runs require an explicit legacy symmetry projection.")
        if args.actual is not None:
            observed = read_results(args.actual)
            report["mode"] = "computed_vs_archive"
            report["actual"] = str(args.actual.resolve())
            report["observed_checks"] = check_physics(observed)
            if report["observed_checks"]["max_ladder_hermiticity_gap"] > 1e-8:
                raise ValueError("Fresh ex9 violates ladder-operator Hermiticity")
            if args.legacy_ex9_symmetry_projection:
                # Pair-exchange symmetry plus Hermiticity makes these channels real.
                if max(abs(row[k]) for row in observed["zpm"]
                       for k in (13, 17)) > 1e-8:
                    raise ValueError("Fresh ex9 violates real even-reflection sector")
            report["comparison"] = compare(
                observed, ref, args.energy_atol, args.observable_atol,
                legacy_ex9_projection=args.legacy_ex9_symmetry_projection)
        report["status"] = "PASS"
    except (OSError, ValueError) as error:
        report.update(status="FAIL", error=str(error))
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
