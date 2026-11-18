#!/usr/bin/env python3
"""Compare both exact cubic anchored-basis coordinate paths bit for bit.

The original translation-anchored basis is used on both sides: only the
integer periodic-coordinate difference cache is toggled. Python >=3.6.
"""
import argparse
import hashlib
import itertools
import os
from pathlib import Path
import re
import subprocess
import tempfile

from cpc_cubic_pointgroup import (canonical_bonds, sector_dimensions,
                                  symmetry_generators, write_input)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def run(driver, work, side, sector, threads, enabled, timeout):
    write_input(work, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=64)
    env = dict(os.environ)
    for key in ("ED2_CUBIC_BASIS_COORD_CACHE", "ED2_FORCE_LEGACY_ORBITS",
                "ED2_CUBIC_BASIS_PROOF_DUMP", "ED2_CUBIC_DIRECT_HOP_TABLE"):
        env.pop(key, None)
    env.update(ED2_CUBIC_FAST_GROUP_VALIDATION="1",
               ED2_CUBIC_ANCHORED_BASIS="1",
               ED2_CUBIC_ANCHORED_REP="1",
               ED2_CUBIC_BASIS_PROOF_DUMP="1",
               ED2_BASIS_DETAIL_TIMING="1",
               OPENBLAS_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE",
               OMP_NUM_THREADS=str(threads))
    if enabled:
        env["ED2_CUBIC_BASIS_COORD_CACHE"] = "1"
    with (work / "input.dat").open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, cwd=str(work), env=env,
                            timeout=timeout, check=False)
    log = cp.stdout.decode("utf-8", errors="replace")
    (work / "driver.log").write_text(log)
    assert cp.returncode == 0, (work, cp.returncode, log[-3000:])
    assert "Basis builder: anchored-cubic" in log
    assert ("Basis orbit coordinates: exact precomputed periodic deltas" in log) == enabled
    dim = sector_dimensions(side)[sector]
    matches = re.findall(r"ED2_BASIS_DRIVER_DIM=(\\d+)", log)
    assert len(matches) == 1 and int(matches[0]) == dim
    proof = work / "cubic_basis_proof_dump.txt"
    assert len(proof.read_bytes().splitlines()) == dim
    timing = re.findall(r"ED2_BASIS_PHASE orbit_scan wall_s=\\s*([+\\-0-9.eE]+)", log)
    assert len(timing) == 1 and float(timing[0]) >= 0
    return sha256(work / "input.dat"), sha256(proof), float(timing[0])


def reject_bad_flag(driver, root):
    side, sector = 4, "oh-a1g"
    work = root / "reject-missing-anchored-basis"
    write_input(work, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=64)
    env = dict(os.environ)
    env.update(ED2_CUBIC_BASIS_COORD_CACHE="1",
               ED2_CUBIC_FAST_GROUP_VALIDATION="1",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    env.pop("ED2_CUBIC_ANCHORED_BASIS", None)
    with (work / "input.dat").open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, cwd=str(work), env=env,
                            timeout=60, check=False)
    assert cp.returncode != 0
    assert b"Basis coordinate cache requires ED2_CUBIC_ANCHORED_BASIS=1" in cp.stdout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver", required=True, type=Path)
    p.add_argument("--sides", default="4,6,8,10")
    p.add_argument("--threads", default="1,4,8")
    p.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    driver = args.driver.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-cubic-coord-cache-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    rows = []
    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}
    for side, sector, thread in itertools.product(sides, sectors, threads):
        times = []
        for enabled in (False, True):
            dest = root / "s{}-{}-t{}-{}".format(
                side, sector, thread, "optimized" if enabled else "original")
            i_hash, b_hash, timing = run(driver, dest, side, sector, thread,
                                         enabled, args.timeout)
            times.append((i_hash, b_hash, timing))
        assert times[0][0] == times[1][0], ("INPUT DIFFERENCE", side, sector, thread)
        assert times[0][1] == times[1][1], ("BASIS BITS DIFFER", side, sector, thread)
        rows.append((side, sector, thread, times[0][2], times[1][2]))
        print("CUBIC_COORD_CACHE_BITS_PASS side={} sector={} threads={} "
              "old={:.6f}s new={:.6f}s".format(*rows[-1]), flush=True)
    reject_bad_flag(driver, root)
    print("CUBIC_COORD_CACHE_REGRESSION_PASS {} exact basis pairs, "
          "missing-prerequisite rejected".format(len(rows)), flush=True)


if __name__ == "__main__":
    main()
