#!/usr/bin/env python3
"""Validate CPC PR17 1000-site three-magnon pilot CSV and manuscript table.

This is an editorial integrity gate, NOT independent reconstruction of the
original PRESTO trial logs and NOT operator exclusivity attestation.
Requires only Python standard library (Python >= 3.6).
"""
import argparse
import csv
import math
from pathlib import Path

SECTORS = ("oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u")
THREADS = (1, 4, 8)
DIMENSIONS = {
    "oh-a1g": 4299,
    "oh-a2g": 3624,
    "oh-a1u": 2884,
    "oh-a2u": 3252,
}
LABELS = {
    "oh-a1g": "A_{1g}",
    "oh-a2g": "A_{2g}",
    "oh-a1u": "A_{1u}",
    "oh-a2u": "A_{2u}",
}
COLUMNS = (
    "sector", "threads", "basis_off_ms", "basis_on_ms",
    "total_off_ms", "total_on_ms",
)


def verify(csv_path, tex_path):
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert tuple(reader.fieldnames or ()) == COLUMNS, "CSV schema changed"
        raw = list(reader)
    assert len(raw) == 12, "Expected 12 sectors-by-threads median records"

    items = {}
    bs_ratios, total_ratios = [], []
    for row in raw:
        sector, threads = row["sector"], int(row["threads"])
        assert sector in SECTORS and threads in THREADS
        key = (sector, threads)
        assert key not in items, ("Duplicate case", key)
        measured = tuple(float(row[name]) for name in COLUMNS[2:])
        assert all(math.isfinite(v) and v > 0 for v in measured)
        b0, b1, t0, t1 = measured
        assert b0 < t0 and b1 < t1, "Basis cannot exceed entire run"
        bs_ratio, total_ratio = b0 / b1, t0 / t1
        assert bs_ratio > 1 and total_ratio > 1, ("No speedup", key)
        assert 2.25 < bs_ratio < 2.35, ("Unexpected PR17 basis ratio", key)
        assert 1.43 < total_ratio < 1.72, ("Unexpected PR17 total ratio", key)
        items[key] = measured
        bs_ratios.append(bs_ratio)
        total_ratios.append(total_ratio)
    assert set(items) == {(s, t) for s in SECTORS for t in THREADS}

    tex = tex_path.read_text(encoding="utf-8")
    assert tex.count("% PR17_PRESTO_T8_BEGIN:") == 1
    assert tex.count("% PR17_PRESTO_T8_END") == 1
    table = tex.split("% PR17_PRESTO_T8_BEGIN:", 1)[1].split(
        "% PR17_PRESTO_T8_END", 1)[0]
    assert r"\label{tab:pr17-pilot}" in tex
    for sector in SECTORS:
        b0, b1, t0, t1 = items[(sector, 8)]
        expected = (
            "$" + LABELS[sector] + "$ & {} & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f}"
            .format(DIMENSIONS[sector], b0, b1, t0, t1, t0 / t1)
        )
        assert expected in table, ("LaTeX/CPC CSV mismatch", sector, expected)

    print("CPC_PR17_PILOT_TABLE_PASS: 12 CSV cases; 4 exact T8 manuscript rows; "
          "basis ratio range {:.3f}-{:.3f}; total ratio range {:.3f}-{:.3f}"
          .format(min(bs_ratios), max(bs_ratios),
                  min(total_ratios), max(total_ratios)), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path,
                   default=Path("paper/data/pr17_presto_coordinate_cache_pilot.csv"))
    p.add_argument("--tex", type=Path, default=Path("paper/cpc_ed2_main.tex"))
    a = p.parse_args()
    verify(a.csv, a.tex)


if __name__ == "__main__":
    main()
