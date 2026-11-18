#!/usr/bin/env python3
"""Prove complete native Hamiltonian cache identity with basis coordinate opt-in.

Compare the original PR16 exact direct-hop path against the SAME path with
only ED2_CUBIC_BASIS_COORD_CACHE=1. Reuse the established production cache
and Ritz drivers, rather than an independently reconstructed Hamiltonian.
Python >= 3.6; standard library only.
"""
import argparse
from contextlib import contextmanager
import itertools
import json
import os
from pathlib import Path
import tempfile

from cpc_cubic_direct_hop_table_regression import cache_case, ritz_case, sha
from cpc_cubic_pointgroup import sector_dimensions


@contextmanager
def coord_mode(enabled):
    key = "ED2_CUBIC_BASIS_COORD_CACHE"
    saved = os.environ.get(key)
    try:
        if enabled:
            os.environ[key] = "1"
        else:
            os.environ.pop(key, None)
        yield
    finally:
        if saved is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = saved


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--sides", default="4,6,10")
    p.add_argument("--threads", default="1,4")
    p.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    p.add_argument("--timeout", type=int, default=1200)
    p.add_argument("--skip-lanczos", action="store_true")
    p.add_argument("--output", type=Path)
    args = p.parse_args()

    driver = args.driver.resolve(strict=True)
    exe = args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-cubic-coordinate-fullcache-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)

    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert sides and all(x in (4, 6, 8, 10) for x in sides)
    assert threads and all(x >= 1 for x in threads)
    assert set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}

    rows = []
    for side, sector, thread in itertools.product(sides, sectors, threads):
        hashes = []
        dims = sector_dimensions(side)
        for enabled in (False, True):
            dest = root / ("new" if enabled else "old")
            with coord_mode(enabled):
                work = cache_case(dest, driver, side, sector, thread, True, args.timeout)
            log = (work / "driver.log").read_text(errors="replace")
            assert ("Basis orbit coordinates: exact precomputed periodic deltas" in log) == enabled
            hashes.append((sha(work / "input.dat"),
                           sha(work / "cubic_cache_proof_dump.txt")))
        assert hashes[0][0] == hashes[1][0], ("INPUT SHA MISMATCH", side, sector, thread)
        assert hashes[0][1] == hashes[1][1], ("COMPLETE CACHE BITS MISMATCH",
                                            side, sector, thread)
        rows.append(dict(side=side, sector=sector, threads=thread,
                         dimension=dims[sector], input_sha256=hashes[0][0],
                         cache_sha256=hashes[0][1]))
        print("CUBIC_COORD_FULL_CACHE_BITS_PASS side={} sector={} threads={} "
              "D={}".format(side, sector, thread, dims[sector]), flush=True)

        if side in (4, 6) and thread == 1 and not args.skip_lanczos:
            outcomes = []
            for enabled in (False, True):
                with coord_mode(enabled):
                    outcomes.append(ritz_case(root / ("new" if enabled else "old"),
                                              exe, side, sector, True, args.timeout))
            old, new = outcomes
            assert old["input_sha256"] == new["input_sha256"]
            assert abs(old["energy"] - new["energy"]) < 1e-8
            assert old["residual"] <= 1.01e-8 and new["residual"] <= 1.01e-8
            print("CUBIC_COORD_RITZ_PASS side={} sector={} deltaE={:+.3e} "
                  "residuals={:.3e}/{:.3e}".format(side, sector,
                    new["energy"] - old["energy"], old["residual"],
                    new["residual"]), flush=True)

    (root / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("CUBIC_COORD_FULL_CACHE_REGRESSION_PASS {} full-cache exact pairs"
          .format(len(rows)), flush=True)


if __name__ == "__main__":
    main()
