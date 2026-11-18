#!/usr/bin/env python3
"""Ensure opt-in cubic basis timing is complete and does not change the Ritz pair.

Runs a matched pair with and without ED2_BASIS_DETAIL_TIMING.  Python >=3.6.
This is a wall-time attribution test, not a performance acceptance benchmark.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import tempfile

from cpc_cubic_pointgroup import canonical_bonds, sector_dimensions, symmetry_generators, write_input
from cpc_pointgroup_profile import run_case

PHASES = ("preflight_other", "group_validation", "orbit_scan", "finalize", "total")


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sides", default="4,6")
    parser.add_argument("--threads", default="1,4")
    args = parser.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-cubic-basis-timing-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    exe = args.exe.resolve(strict=True)
    report = []
    for side in (int(v) for v in args.sides.split(",")):
        sector = "oh-a1g"
        dimension = sector_dimensions(side)[sector]
        bonds = canonical_bonds(side)
        for threads in (int(v) for v in args.threads.split(",")):
            paired = []
            for timed in (False, True):
                work = root / "s{}-t{}-{}".format(side, threads, "timed" if timed else "default")
                write_input(work, side, 3, symmetry_generators(side, sector), bonds, mnte=64)
                overrides = {
                    "ED2_CUBIC_ANCHORED_REP": "1",
                    "ED2_ALLOW_LARGE_GROUP_VALIDATION": "1",
                    "ED2_BASIS_DETAIL_TIMING": "1" if timed else None,
                }
                previous = {k: os.environ.get(k) for k in overrides}
                try:
                    for key, value in overrides.items():
                        if value is None:
                            os.environ.pop(key, None)
                        else:
                            os.environ[key] = value
                    result = run_case(exe, work, side ** 3, 3, sector, dimension,
                                      threads, "two-stage", 600)
                finally:
                    for key, value in previous.items():
                        if value is None:
                            os.environ.pop(key, None)
                        else:
                            os.environ[key] = value
                output = (work / "stdout.log").read_text(errors="replace")
                pairs = re.findall(
                    r"ED2_BASIS_PHASE ([a-z_]+) wall_s=\s*([+\-0-9.eE]+)",
                    output
                )
                phases = {}
                for key, value in pairs:
                    assert key not in phases, ("duplicate basis timing phase", key)
                    phases[key] = float(value)
                if timed:
                    assert set(phases) == set(PHASES), (phases, PHASES)
                    assert all(math.isfinite(x) and x >= 0 for x in phases.values())
                    assert phases["group_validation"] > 0
                    total = sum(phases[k] for k in PHASES if k != "total")
                    assert abs(total - phases["total"]) < 1e-7, phases
                    assert abs(phases["total"] - result["phase_wall_s"]["basis"]) < 0.05, (
                        phases["total"], result["phase_wall_s"]["basis"])
                    print("BASIS_DETAIL_TIMING_PASS side={} threads={} preflight={:.6f}s "
                          "validation={:.6f}s orbit={:.6f}s finalize={:.6f}s".format(
                              side, threads, phases["preflight_other"],
                              phases["group_validation"], phases["orbit_scan"],
                              phases["finalize"]), flush=True)
                else:
                    assert not phases, "disabled basis detail timing unexpectedly logged"
                paired.append((result, phases))
            a, b = paired[0][0], paired[1][0]
            assert a["input_sha256"] == b["input_sha256"]
            assert a["dimension"] == b["dimension"] == dimension
            assert abs(a["energy"]-b["energy"]) <= 1e-8
            assert a["residual"] <= 1.01e-8 and b["residual"] <= 1.01e-8
            report.append(dict(side=side,threads=threads,energy=b["energy"],
                               basis_breakdown_s=paired[1][1]))
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print("CUBIC_BASIS_DETAIL_TIMING_REGRESSION_PASS {} matched pairs".format(len(report)),
          flush=True)


if __name__ == "__main__":
    run()
