#!/usr/bin/env python3
"""Exact regression for skipping UNUSED cubic orbit-key hash when dense-hop is enabled.

The direct-pair table contains every actual lookup for the opt-in
cubic 3-magnon Gamma/O_h Hamiltonian cache. This test deliberately runs
the real native ED2 Fortran code, not a Python reimplementation, and
compares the COMPLETE target/complex-bit Hamiltonian cache dump.
Both paths share the same PR19 basis pruning and PR17 coordinate cache.
Python >=3.6 standard library only.
"""
import argparse
from contextlib import contextmanager
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

from cpc_cubic_direct_hop_table_regression import cache_case, make_env, ritz_case, sha
from cpc_cubic_pointgroup import canonical_bonds, symmetry_generators, write_input

OPT_FLAG = "ED2_CUBIC_SKIP_UNUSED_HASH"
MARKER = "Representative lookup: cubic-direct-only-unused-hash-elided"
COORD_FLAG = "ED2_CUBIC_BASIS_COORD_CACHE"
PRUNE_FLAG = "ED2_CUBIC_BASIS_PRUNED_CANON"


@contextmanager
def mode(enabled):
    keys = (OPT_FLAG, COORD_FLAG, PRUNE_FLAG)
    previous = {key: os.environ.get(key) for key in keys}
    try:
        os.environ[COORD_FLAG] = "1"
        os.environ[PRUNE_FLAG] = "1"
        if enabled:
            os.environ[OPT_FLAG] = "1"
        else:
            os.environ.pop(OPT_FLAG, None)
        yield
    finally:
        for key in keys:
            if previous[key] is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = previous[key]


def assert_marker(path, enabled):
    txt = (path / "driver.log").read_text(errors="replace")
    assert (MARKER in txt) == enabled, ("Direct-only provenance mismatch", path)
    assert "Representative lookup: cubic-direct-hop-pairs" in txt
    assert "Basis orbit coordinates: exact precomputed periodic deltas" in txt
    assert "Basis canonical enumeration: exact-axis-pruned" in txt


def reject_bad_flags(driver, root):
    for name, direct, value, expected in (
        ("no-direct", False, "1",
         "skipping unused hash requires ED2_CUBIC_DIRECT_HOP_TABLE=1"),
        ("invalid-value", True, "2",
         "ED2_CUBIC_SKIP_UNUSED_HASH must be 0 or 1"),
    ):
        dest = root / ("reject-" + name)
        write_input(dest, 4, 3, symmetry_generators(4, "oh-a1g"),
                    canonical_bonds(4), mnte=64)
        with mode(False):
            env = make_env(1, direct)
        env[OPT_FLAG] = value
        env[COORD_FLAG] = "1"
        env[PRUNE_FLAG] = "1"
        with (dest / "input.dat").open("rb") as inp:
            p = subprocess.run([str(driver)], stdin=inp, cwd=str(dest),
                               env=env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=60,
                               check=False)
        assert p.returncode != 0, ("Unsafe hash bypass accepted", name)
        assert expected.encode("utf-8") in p.stdout, (name, p.stdout[-1000:])


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
    driver, exe = args.driver.resolve(strict=True), args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-skip-unused-hash-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert sides and all(s in (4,6,8,10) for s in sides)
    assert threads and all(t >= 1 for t in threads)
    assert sectors and set(sectors) <= {"oh-a1g","oh-a2g","oh-a1u","oh-a2u"}
    count_cache = count_ritz = 0
    report = []
    for side, sector, t in itertools.product(sides, sectors, threads):
        cache = []
        for enabled in (False, True):
            with mode(enabled):
                dest = cache_case(root / ("hash" + str(int(enabled))),
                                  driver, side, sector, t, True, args.timeout)
            assert_marker(dest, enabled)
            cache.append((sha(dest / "input.dat"),
                          sha(dest / "cubic_cache_proof_dump.txt")))
        assert cache[0] == cache[1], (
            "Skipping unused hash changed native Hamiltonian bits",
            side, sector, t)
        count_cache += 1
        report.append(dict(side=side, sector=sector, threads=t,
                           cache_sha256=cache[0][1], input_sha256=cache[0][0]))
        print("CUBIC_SKIP_HASH_FULL_CACHE_BITS_PASS side={} sector={} t={}".format(
            side, sector, t), flush=True)
        if side in (4,6) and t == 1 and not args.skip_lanczos:
            obs = []
            for enabled in (False, True):
                with mode(enabled):
                    obs.append(ritz_case(root / ("ritz" + str(int(enabled))),
                                         exe, side, sector, True, args.timeout))
            a,b = obs
            assert a["input_sha256"] == b["input_sha256"]
            assert abs(a["energy"] - b["energy"]) < 1e-8
            assert a["residual"] <= 1.01e-8 and b["residual"] <= 1.01e-8
            assert a["h_applications"] == b["h_applications"]
            count_ritz += 1
            print("CUBIC_SKIP_HASH_RITZ_PASS side={} sector={} dE={:+.3e}".format(
                side, sector, b["energy"] - a["energy"]), flush=True)
    reject_bad_flags(driver, root)
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print("CUBIC_SKIP_HASH_REGRESSION_PASS cache_pairs={} ritz_pairs={} "
          "bad_flags_rejected=2".format(count_cache, count_ritz), flush=True)


if __name__ == "__main__":
    main()
