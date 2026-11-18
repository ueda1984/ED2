#!/usr/bin/env python3
"""Full native bitwise oracle for compact real H-cache storage.

Compare ED2_REAL_CACHED_MATVEC=1 with packing OFF/ON using complex
deterministic vector entries. All historical phase, basis and matrix
construction options remain identical; only physical coefficient storage
changes. Includes independent real mixed-spin pass and fail-closed tests.
Python >=3.6 / stdlib only.
"""
import argparse
import itertools
import os
from pathlib import Path
import tempfile

from cpc_cubic_pointgroup import (
    canonical_bonds, sector_dimensions, symmetry_generators, write_input)
from cpc_real_cached_matvec_regression import (
    MARKER as REAL_MARKER, cubic_env, generic_env, run, sha)
from cpc_portable_dm_quickstart import write_input as write_dm_input

FLAG = "ED2_PACKED_REAL_CACHE"
MARKER = "Hamiltonian action: compact-real cache extra_bytes="


def env_cubic(threads, packed):
    env = cubic_env(threads, True)
    env.pop(FLAG, None)
    if packed:
        env[FLAG] = "1"
    return env


def generic(root, exe, name, packed, timeout, real=True):
    dest = root / ("generic-{}-packed{}".format(name, int(packed)))
    write_dm_input(dest, 4, 0.0 if real else 0.12)
    if real:
        bonds = []
        for line in (dest / "bonds.dat").read_text().splitlines():
            parts = line.split()
            parts[-1] = "0"
            bonds.append(" ".join(parts))
        (dest / "bonds.dat").write_text("\n".join(bonds) + "\n")
    env = generic_env("1")
    env.pop(FLAG, None)
    if packed:
        env[FLAG] = "1"
    return dest, *run(exe, dest, env, timeout)


def reject_flags(root, exe, timeout):
    src = root / "reject-no-real"
    write_dm_input(src, 4, 0.0)
    env = generic_env(None)
    env[FLAG] = "1"
    code, output = run(exe, src, env, timeout)
    assert code != 0
    assert "ED2_PACKED_REAL_CACHE requires ED2_REAL_CACHED_MATVEC=1" in output
    print("PACKED_REAL_MISSING_REAL_MODE_REJECTED", flush=True)

    src = root / "reject-bad-packed"
    write_dm_input(src, 4, 0.0)
    env = generic_env("1")
    env[FLAG] = "2"
    code, output = run(exe, src, env, timeout)
    assert code != 0
    assert "ED2_PACKED_REAL_CACHE must be 0 or 1" in output
    print("PACKED_REAL_INVALID_FLAG_REJECTED", flush=True)

    src = root / "reject-nonreal"
    write_dm_input(src, 4, 0.12)
    env = generic_env("1")
    env[FLAG] = "1"
    code, output = run(exe, src, env, timeout)
    assert code != 0
    assert "ED2_REAL_CACHED_MATVEC requires exactly real matrix elements" in output
    print("PACKED_REAL_NONREAL_MATRIX_REJECTED", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--sides", default="4,6,10")
    p.add_argument("--threads", default="1,4")
    p.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--skip-generic", action="store_true")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    exe, driver = args.exe.resolve(strict=True), args.driver.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-packed-real-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)

    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert sides and all(s in (4, 6, 8, 10) for s in sides)
    assert threads and all(t >= 1 for t in threads)
    assert set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"} and sectors

    count = 0
    for side, sector, t in itertools.product(sides, sectors, threads):
        proofs = []
        for packed in (False, True):
            dest = root / "n{}-{}-t{}-packed{}".format(
                side, sector, t, int(packed))
            write_input(dest, side, 3, symmetry_generators(side, sector),
                        canonical_bonds(side), mnte=64)
            rc, log = run(driver, dest, env_cubic(t, packed), args.timeout)
            assert rc == 0, (dest, log[-3000:])
            assert REAL_MARKER in log
            assert (MARKER in log) == packed
            assert "ED2_CUBIC_MATVEC_PROOF_ROWS={}".format(
                sector_dimensions(side)[sector]) in log
            proof = dest / "cubic_matvec_proof_dump.txt"
            assert len(proof.read_bytes().splitlines()) == sector_dimensions(side)[sector]
            proofs.append((sha(dest / "input.dat"), sha(proof)))
        assert proofs[0] == proofs[1], (
            "Compact real coefficient cache altered native action bits", side, sector, t)
        print("PACKED_REAL_MATVEC_EXACT_BITS_PASS side={} sector={} threads={}".format(
            side, sector, t), flush=True)
        count += 1

    if not args.skip_generic:
        off, rc, log = generic(root, exe, "mixed-real", False, args.timeout)
        assert rc == 0
        on, rc, log = generic(root, exe, "mixed-real", True, args.timeout)
        assert rc == 0, log[-3000:]
        assert MARKER in log
        assert sha(off / "eigenvalues.dat") == sha(on / "eigenvalues.dat")
        print("PACKED_REAL_GENERIC_MIXED_SPIN_BITS_PASS", flush=True)
        reject_flags(root, exe, args.timeout)

    print("PACKED_REAL_CACHE_REGRESSION_PASS exact_pairs={} generic_tested={}".format(
        count, not args.skip_generic), flush=True)


if __name__ == "__main__":
    main()
