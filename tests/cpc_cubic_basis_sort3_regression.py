#!/usr/bin/env python3
"""Bit-for-bit scientific proof for opt-in cubic three-site sorting network.

Compare two native production paths: both use the PR17 periodic coordinate
cache, while only ED2_CUBIC_BASIS_SORT3_NETWORK is toggled. Every historical
basis representative, norm bit, entire native local H-cache record and small
native Ritz/residual must agree. Python 3.6+, stdlib only.
"""
import argparse
from contextlib import contextmanager
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

from cpc_cubic_basis_coordinate_cache_regression import run as basis_case
from cpc_cubic_direct_hop_table_regression import cache_case, ritz_case, sha
from cpc_cubic_pointgroup import (canonical_bonds, sector_dimensions,
                                  symmetry_generators, write_input)

MARKER = "Basis canonical sorter: exact-three-site-network"


@contextmanager
def sorter_flag(enabled):
    keys = ("ED2_CUBIC_BASIS_COORD_CACHE", "ED2_CUBIC_BASIS_SORT3_NETWORK")
    saved = dict((key, os.environ.get(key)) for key in keys)
    try:
        os.environ["ED2_CUBIC_BASIS_COORD_CACHE"] = "1"
        if enabled:
            os.environ["ED2_CUBIC_BASIS_SORT3_NETWORK"] = "1"
        else:
            os.environ.pop("ED2_CUBIC_BASIS_SORT3_NETWORK", None)
        yield
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def check_native_marker(work, enabled, logname):
    text = (work / logname).read_text(errors="replace")
    assert (MARKER in text) == enabled, ("SORT3 PROVENANCE", work)
    assert "Basis orbit coordinates: exact precomputed periodic deltas" in text


def reject_missing_coordinate(driver, root):
    dest = root / "reject-no-coordinates"
    side, sector = 4, "oh-a1g"
    write_input(dest, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=64)
    env = dict(os.environ)
    env.update(ED2_CUBIC_ANCHORED_BASIS="1",
               ED2_CUBIC_FAST_GROUP_VALIDATION="1",
               ED2_CUBIC_BASIS_SORT3_NETWORK="1",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    env.pop("ED2_CUBIC_BASIS_COORD_CACHE", None)
    env.pop("ED2_FORCE_LEGACY_ORBITS", None)
    with (dest / "input.dat").open("rb") as inp:
        proc = subprocess.run([str(driver)], stdin=inp, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, cwd=str(dest), env=env,
                              timeout=60, check=False)
    assert proc.returncode != 0
    assert b"Basis three-site sorter requires ED2_CUBIC_BASIS_COORD_CACHE=1" in proc.stdout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--basis-driver", type=Path, required=True)
    p.add_argument("--cache-driver", type=Path, required=True)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--sides", default="4,6,8,10")
    p.add_argument("--threads", default="1,4")
    p.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    p.add_argument("--timeout", type=int, default=1200)
    p.add_argument("--skip-lanczos", action="store_true")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    basis_driver = args.basis_driver.resolve(strict=True)
    cache_driver = args.cache_driver.resolve(strict=True)
    exe = args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-cubic-sort3-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)

    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert sides and all(x in (4, 6, 8, 10) for x in sides)
    assert threads and all(x >= 1 for x in threads)
    assert sectors and set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}
    rows = []
    basis_pairs = cache_pairs = ritz_pairs = 0

    for side, sector, thread in itertools.product(sides, sectors, threads):
        native = []
        for enabled in (False, True):
            label = "on" if enabled else "off"
            work = root / "basis-{}-n{}-{}-t{}".format(label, side, sector, thread)
            with sorter_flag(enabled):
                input_sha, basis_sha, elapsed = basis_case(
                    basis_driver, work, side, sector, thread, True, args.timeout)
            check_native_marker(work, enabled, "driver.log")
            native.append((input_sha, basis_sha, elapsed))
        assert native[0][:2] == native[1][:2], (
            "SORT3 CHANGED HISTORICAL BASIS BITS", side, sector, thread)
        basis_pairs += 1
        print("CUBIC_SORT3_BASIS_BITS_PASS side={} sector={} t={} "
              "old={:.6f}s new={:.6f}s".format(side, sector, thread,
                                               native[0][2], native[1][2]), flush=True)

        if side in (4, 6, 10):
            hashes = []
            for enabled in (False, True):
                label = "on" if enabled else "off"
                with sorter_flag(enabled):
                    work = cache_case(root / ("cache-" + label), cache_driver,
                                      side, sector, thread, True, args.timeout)
                check_native_marker(work, enabled, "driver.log")
                hashes.append((sha(work / "input.dat"),
                               sha(work / "cubic_cache_proof_dump.txt")))
            assert hashes[0] == hashes[1], (
                "SORT3 CHANGED COMPLETE HAMILTONIAN CACHE BITS", side, sector, thread)
            cache_pairs += 1
            print("CUBIC_SORT3_FULL_CACHE_BITS_PASS side={} sector={} t={} D={}".format(
                side, sector, thread, sector_dimensions(side)[sector]), flush=True)

            if side in (4, 6) and thread == 1 and not args.skip_lanczos:
                observed = []
                for enabled in (False, True):
                    with sorter_flag(enabled):
                        obs = ritz_case(root / ("ritz-" + ("on" if enabled else "off")),
                                        exe, side, sector, True, args.timeout)
                    observed.append(obs)
                x, y = observed
                assert x["input_sha256"] == y["input_sha256"]
                assert abs(x["energy"] - y["energy"]) < 1e-8
                assert x["residual"] <= 1.01e-8 and y["residual"] <= 1.01e-8
                assert x["h_applications"] == y["h_applications"]
                ritz_pairs += 1
                print("CUBIC_SORT3_RITZ_PASS side={} sector={} dE={:+.3e}".format(
                    side, sector, y["energy"] - x["energy"]), flush=True)
        rows.append(dict(side=side, sector=sector, threads=thread,
                         basis_sha256=native[0][1], old_orbit_s=native[0][2],
                         new_orbit_s=native[1][2]))
    reject_missing_coordinate(basis_driver, root)
    (root / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("CUBIC_SORT3_REGRESSION_PASS basis_pairs={} cache_pairs={} "
          "ritz_pairs={} bad_flag_rejected=1".format(
              basis_pairs, cache_pairs, ritz_pairs), flush=True)


if __name__ == "__main__":
    main()
