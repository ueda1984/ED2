#!/usr/bin/env python3
"""Exact native matrix-vector equivalence for opt-in real cached coefficients.

The production operator is evaluated on a fixed *complex* test vector, with
the existing complex kernel and with ED2_REAL_CACHED_MATVEC=1. Every output
real/imaginary IEEE-754 bit is compared across all cubic 1D O_h sectors.
Extra generic mixed-spin tests prove that nonreal matrix elements reject the
optimization while a real full-cache Hamiltonian remains supported.
Python >= 3.6, standard library only.
"""
import argparse
import hashlib
import itertools
import os
from pathlib import Path
import re
import subprocess
import tempfile

from cpc_cubic_pointgroup import (
    canonical_bonds, sector_dimensions, symmetry_generators, write_input)
from cpc_portable_dm_quickstart import write_input as write_dm_input

MARKER = "Hamiltonian action: strictly-real cached coefficients"
FLAG = "ED2_REAL_CACHED_MATVEC"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def run(driver, dest, env, timeout):
    with (dest / "input.dat").open("rb") as inp:
        completed = subprocess.run(
            [str(driver)], cwd=str(dest), stdin=inp, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=timeout, check=False)
    output = completed.stdout.decode("utf-8", errors="replace")
    (dest / "driver.log").write_text(output)
    return completed.returncode, output


def cubic_env(thread, enabled):
    env = dict(os.environ)
    for key in tuple(env):
        if key.startswith("ED2_CUBIC_") or key == FLAG or key == "ED2_FORCE_LEGACY_ORBITS":
            env.pop(key, None)
    env.update(
        ED2_CUBIC_ANCHORED_REP="1",
        ED2_CUBIC_FAST_GROUP_VALIDATION="1",
        ED2_CUBIC_ANCHORED_BASIS="1",
        ED2_CUBIC_LOCAL_BOND_CACHE="1",
        ED2_CUBIC_DIRECT_HOP_TABLE="1",
        ED2_CUBIC_BASIS_COORD_CACHE="1",
        ED2_CUBIC_BASIS_PRUNED_CANON="1",
        ED2_CUBIC_SKIP_UNUSED_HASH="1",
        ED2_CUBIC_MATVEC_PROOF_DUMP="1",
        OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OMP_NUM_THREADS=str(thread))
    if enabled:
        env[FLAG] = "1"
    return env


def generic_env(value):
    env = dict(os.environ)
    for key in tuple(env):
        if key.startswith("ED2_CUBIC_") or key.startswith("ED2_FORCE_LEGACY"):
            env.pop(key, None)
    env.update(ED2_RNG_SEED="19", OPENBLAS_NUM_THREADS="1",
               OMP_DYNAMIC="FALSE", OMP_NUM_THREADS="1")
    if value is None:
        env.pop(FLAG, None)
    else:
        env[FLAG] = value
    return env


def generic_case(exe, root, kind, flag, timeout):
    work = root / ("generic-{}-{}".format(kind, flag or "off"))
    write_dm_input(work, 4, 0.12 if kind == "nonreal-dm" else 0.0)
    if kind in ("real-mixed", "partial-real"):
        records = []
        for line in (work / "bonds.dat").read_text().splitlines():
            parts = line.split()
            parts[-1] = "0"  # Remove the SxSy+SySx term.
            records.append(" ".join(parts))
        (work / "bonds.dat").write_text("\n".join(records) + "\n")
    if kind == "partial-real":
        deck = (work / "input.dat").read_text()
        assert deck.count("wk_dim=-1") == 1
        (work / "input.dat").write_text(deck.replace("wk_dim=-1", "wk_dim=1"))
    code, output = run(exe, work, generic_env(flag), timeout)
    return work, code, output


def check_generic(exe, root, timeout):
    a, code, output = generic_case(exe, root, "real-mixed", None, timeout)
    assert code == 0, output[-2500:]
    b, code, output = generic_case(exe, root, "real-mixed", "1", timeout)
    assert code == 0, output[-2500:]
    assert MARKER in output
    assert sha(a / "eigenvalues.dat") == sha(b / "eigenvalues.dat")
    print("REAL_MATVEC_GENERIC_REAL_MODEL_PASS", flush=True)

    _, code, output = generic_case(exe, root, "nonreal-dm", "1", timeout)
    assert code != 0
    assert "ED2_REAL_CACHED_MATVEC requires exactly real matrix elements" in output
    print("REAL_MATVEC_NONREAL_MODEL_REJECTED", flush=True)

    _, code, output = generic_case(exe, root, "partial-real", "1", timeout)
    assert code != 0
    assert "ED2_REAL_CACHED_MATVEC requires fully cached matrix" in output
    print("REAL_MATVEC_PARTIAL_CACHE_REJECTED", flush=True)

    _, code, output = generic_case(exe, root, "real-mixed", "2", timeout)
    assert code != 0
    assert "ED2_REAL_CACHED_MATVEC must be 0 or 1" in output
    print("REAL_MATVEC_INVALID_SETTING_REJECTED", flush=True)


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
    driver = args.driver.resolve(strict=True)
    exe = args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-real-matvec-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)

    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert all(s in (4, 6, 8, 10) for s in sides) and sides
    assert threads and all(t >= 1 for t in threads)
    assert sectors and set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}
    total = 0
    for side, sector, threads_value in itertools.product(sides, sectors, threads):
        results = []
        for enabled in (False, True):
            dest = root / "cubic-s{}-{}-t{}-real{}".format(
                side, sector, threads_value, int(enabled))
            write_input(dest, side, 3, symmetry_generators(side, sector),
                        canonical_bonds(side), mnte=64)
            code, output = run(driver, dest, cubic_env(threads_value, enabled),
                               args.timeout)
            assert code == 0, (dest, output[-3000:])
            assert (MARKER in output) == enabled, (dest, "native opt-in marker")
            assert "Basis canonical enumeration: exact-axis-pruned" in output
            assert "Representative lookup: cubic-direct-only-unused-hash-elided" in output
            assert "ED2_CUBIC_MATVEC_PROOF_ROWS={}".format(
                sector_dimensions(side)[sector]) in output
            proof = dest / "cubic_matvec_proof_dump.txt"
            assert len(proof.read_bytes().splitlines()) == sector_dimensions(side)[sector]
            results.append((sha(dest / "input.dat"), sha(proof)))
        assert results[0] == results[1], (
            "Cached real H-application changed exact output bits", side,
            sector, threads_value)
        total += 1
        print("REAL_MATVEC_EXACT_BITS_PASS side={} sector={} threads={}".format(
            side, sector, threads_value), flush=True)

    if not args.skip_generic:
        check_generic(exe, root, args.timeout)
    print("REAL_MATVEC_REGRESSION_PASS exact_pairs={} generic_tested={}".format(
        total, not args.skip_generic), flush=True)


if __name__ == "__main__":
    main()
