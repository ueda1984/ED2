#!/usr/bin/env python3
"""Fail closed if the 25 DM response points in CPC LaTeX differ from source CSV.

The paper has deliberately plotted the independent reference data, not
a smoothed fit. This checks all five sizes, five D values, all ordinates
and legend labels against the version-controlled source CSV. It uses
only Python's standard library and is intentionally separate from
Fortran solver tests: it certifies figure transcription, not physical
truth of the reference calculation.
"""
import argparse
import csv
import math
from pathlib import Path
import re

PLOT = re.compile(
    r"\\addplot\+\[[^\]]+\]\s+coordinates\s*\{(?P<points>[^}]+)\};"
    r"\s*\\addlegendentry\{\$L=(?P<size>\d+)\$\}",
    re.S,
)
POINT = re.compile(
    r"\(\s*(?P<x>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
    r"\s*,\s*(?P<y>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*\)"
)
EXPECTED_SIZES = (4, 6, 8, 10, 12)
EXPECTED_D = (-0.24, -0.12, 0., 0.12, 0.24)


def read_reference(path):
    with path.open(newline="", encoding="utf-8") as stream:
        records = list(csv.DictReader(stream))
    if len(records) != len(EXPECTED_SIZES) * len(EXPECTED_D):
        raise ValueError("Unexpected physical reference table size")
    expected = {}
    for row in records:
        l = int(row["sites"])
        d = float(row["D"])
        y = float(row["staggered_chirality_per_bond"])
        if l not in EXPECTED_SIZES or not math.isfinite(y):
            raise ValueError("Invalid reference model entry")
        matches = [v for v in EXPECTED_D if abs(v-d) < 1e-12]
        if len(matches) != 1:
            raise ValueError("Unexpected D reference grid")
        key = (l, matches[0])
        if key in expected:
            raise ValueError("Duplicate physical data point")
        expected[key] = y
    return expected


def check_figure(path, expected):
    text = path.read_text(encoding="utf-8")
    start = text.find(r"\label{fig:dm-response}")
    if start < 0:
        raise ValueError("Figure label missing")
    parsed = PLOT.findall(text)
    if len(parsed) != len(EXPECTED_SIZES):
        raise ValueError("Expected exactly five plotted finite-size curves")
    actual = {}
    for i, (coordinates, size) in enumerate(parsed):
        l = int(size)
        if l != EXPECTED_SIZES[i]:
            raise ValueError("Figure curve ordering/legend mismatch")
        points = list(POINT.finditer(coordinates))
        if len(points) != len(EXPECTED_D):
            raise ValueError("Each L must have five physical D points")
        for j, match in enumerate(points):
            x, y = float(match.group("x")), float(match.group("y"))
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError("Nonfinite plotted physical data")
            d = EXPECTED_D[j]
            if abs(x-d) > 1e-12:
                raise ValueError("Plot x-grid or ordering mismatch")
            if abs(y-expected[(l, d)]) > 5e-12:
                raise ValueError(
                    f"DM response figure differs from independent CSV: "
                    f"L={l}, D={d}, CSV={expected[(l,d)]}, plotted={y}"
                )
            if (l, d) in actual:
                raise ValueError("Duplicate plotted point")
            actual[(l, d)] = y
        remaining = POINT.sub("", coordinates)
        if remaining.strip():
            raise ValueError("Unexpected/unparsed coordinate data")
    if set(actual) != set(expected):
        raise ValueError("Missing reference or plotted values")
    return len(actual)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path, default=Path("docs/development/cpc-stationary-dm-reference.csv"))
    p.add_argument("--tex", type=Path, default=Path("paper/cpc_ed2_main.tex"))
    args = p.parse_args()
    count = check_figure(args.tex, read_reference(args.csv))
    print(f"CPC DM FIGURE: {count}/25 independent CSV points match LaTeX, 5/5 size curves")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
