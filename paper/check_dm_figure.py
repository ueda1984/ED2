#!/usr/bin/env python3
"""Fail closed if the scientific manuscript data do not match independent CSV.

Checks the 25 chirality points, 25 first gaps, five D=0 susceptibilities
and the five corresponding rounded table entries. This checks *transcription*,
not the scientific validity of the independent solver calculations.
"""
import argparse
import csv
import math
import re
from pathlib import Path

SIZES = (4, 6, 8, 10, 12)
DM = (-.24, -.12, 0., .12, .24)
FIGURE = re.compile(r"\\begin\{figure\}(?:\[[^\]]*\])?(.*?)\\end\{figure\}", re.S)
TABLE = re.compile(r"\\begin\{table\}(?:\[[^\]]*\])?(.*?)\\end\{table\}", re.S)
PLOT = re.compile(r"\\addplot\+\[[^\]]+\]\s*coordinates\s*\{([^}]+)\};", re.S)
LEGEND = re.compile(r"\\addlegendentry\{\$L=(\d+)\$\}")
POINT = re.compile(
    r"\(\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
    r"\s*,\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*\)"
)
ROW = re.compile(r"(?m)^\s*(\d+)\s*&\s*([+-]?\d+(?:\.\d+)?)\s*\\\\\s*$")


def reference(path):
    with path.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 25:
        raise ValueError("Expected 25 independent source rows")
    chirality, gaps, chi = {}, {}, {}
    for row in rows:
        size, d = int(row["sites"]), float(row["D"])
        candidates = [v for v in DM if abs(d - v) < 1e-12]
        if size not in SIZES or len(candidates) != 1:
            raise ValueError("Unexpected finite-size parameter grid")
        key = (size, candidates[0])
        if key in chirality:
            raise ValueError("Duplicate reference point")
        chirality[key] = float(row["staggered_chirality_per_bond"])
        gaps[key] = float(row["first_gap"])
        if not math.isfinite(chirality[key]) or not math.isfinite(gaps[key]) or gaps[key] <= 0:
            raise ValueError("Invalid physical reference number")
        if candidates[0] == 0:
            chi[size] = float(row["chi_kappa_per_bond"])
            if not math.isfinite(chi[size]):
                raise ValueError("Invalid susceptibility reference")
    if len(chirality) != 25 or len(gaps) != 25 or len(chi) != 5:
        raise ValueError("Incomplete reference grid")
    return chirality, gaps, chi


def labeled_block(text, regex, label):
    matches = [block.group(1) for block in regex.finditer(text)
               if r"\label{" + label + "}" in block.group(1)]
    if len(matches) != 1:
        raise ValueError("Expected exactly one labeled block: " + label)
    return matches[0]


def points(text, required):
    found = POINT.findall(text)
    if len(found) != required or POINT.sub("", text).strip():
        raise ValueError("Unexpected/ambiguous coordinate syntax")
    coords = [(float(x), float(y)) for x, y in found]
    if not all(math.isfinite(x) and math.isfinite(y) for x, y in coords):
        raise ValueError("Nonfinite plot coordinates")
    return coords


def five_curves(text, figure_label, expected):
    body = labeled_block(text, FIGURE, figure_label)
    coords = PLOT.findall(body)
    labels = list(map(int, LEGEND.findall(body)))
    if len(coords) != 5 or labels != list(SIZES):
        raise ValueError("Incorrect curves or legend order in " + figure_label)
    for size, curve in zip(SIZES, coords):
        for (x, y), d in zip(points(curve, 5), DM):
            if abs(x - d) > 1e-12 or abs(y - expected[(size, d)]) > 5e-12:
                raise ValueError(f"{figure_label}: L={size}, D={d} mismatches CSV")


def susceptibilities(text, expected):
    body = labeled_block(text, FIGURE, "fig:dm-susceptibility")
    curves = PLOT.findall(body)
    if len(curves) != 1:
        raise ValueError("Expected one susceptibility plot")
    for (x, y), size in zip(points(curves[0], 5), SIZES):
        if abs(x - size) > 1e-12 or abs(y - expected[size]) > 5e-12:
            raise ValueError(f"Susceptibility graph disagrees at L={size}")
    body = labeled_block(text, TABLE, "tab:susc")
    entries = ROW.findall(body)
    if len(entries) != 5:
        raise ValueError("Expected five susceptibility table entries")
    for (size, printed), expected_size in zip(entries, SIZES):
        if int(size) != expected_size or abs(float(printed) - expected[expected_size]) > 5.01e-7:
            raise ValueError(f"Susceptibility table mismatch at L={expected_size}")


def check(text, src):
    chirality, gaps, chi = src
    five_curves(text, "fig:dm-response", chirality)
    five_curves(text, "fig:dm-gap", gaps)
    susceptibilities(text, chi)


def negative_tests(text, src):
    probes = (
        ("response sign", "(0.24,-0.1481511796)", "(0.24,0.1481511796)"),
        ("gap", "(0.24,0.17255018)", "(0.24,0.17255019)"),
        ("chi plot", "(12,0.5085009862217026)", "(12,0.5085019862217026)"),
        ("chi table", r"12 & 0.508501\\", r"12 & 0.608501\\"),
    )
    for name, original, corrupt in probes:
        if text.count(original) != 1:
            raise ValueError(f"Missing unique self-test coordinate: {name}")
        try:
            check(text.replace(original, corrupt, 1), src)
        except ValueError:
            continue
        raise AssertionError(f"Checker incorrectly accepted corrupted {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=Path("docs/development/cpc-stationary-dm-reference.csv"))
    parser.add_argument("--tex", type=Path, default=Path("paper/cpc_ed2_main.tex"))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    src = reference(args.csv)
    text = args.tex.read_text(encoding="utf-8")
    check(text, src)
    print("CPC FIGURE GATE: 25/25 chirality + 25/25 gaps + "
          "5/5 susceptibility points + 5/5 rounded table rows match the CSV")
    if args.self_test:
        negative_tests(text, src)
        print("CPC FIGURE GATE: 4/4 deliberate corruptions rejected")


if __name__ == "__main__":
    main()
