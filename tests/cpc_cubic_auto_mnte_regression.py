#!/usr/bin/env python3
"""Exact native automatic MNTE test for the opt-in three-magnon cubic cache.

Compare the complete projected Hamiltonian cache, byte for byte, when native
ED2 selects the safe local bound automatically (MNTE=-1 in the input) and when
the SAME geometry uses a manually validated MNTE=64. Compare two-stage Ritz
physics for selected cases. The input decks are identical except for MNTE.
No performance result is asserted on shared CI runners; no PRESTO submission.
Python >= 3.6, no NumPy requirement.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from cpc_cubic_pointgroup import (
    canonical_bonds, check_hamiltonian, sector_dimensions,
    symmetry_generators, write_input)
from cpc_pointgroup_profile import run_case

FLAG_NAMES = (
    "ED2_CUBIC_ANCHORED_REP", "ED2_CUBIC_FAST_GROUP_VALIDATION",
    "ED2_CUBIC_ANCHORED_BASIS", "ED2_CUBIC_LOCAL_BOND_CACHE",
    "ED2_CUBIC_DIRECT_HOP_TABLE", "ED2_CUBIC_SKIP_UNUSED_HASH",
    "ED2_CUBIC_BASIS_COORD_CACHE", "ED2_CUBIC_BASIS_PRUNED_CANON",
)
MARKER = "Cubic cache capacity: auto-MNTE=64, conservative 3x18 degree bound"


def sha(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as src:
        for block in iter(lambda: src.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def environment(threads):
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("ED2_CUBIC_") or key == "ED2_FORCE_LEGACY_ORBITS":
            env.pop(key)
    for name in FLAG_NAMES:
        env[name] = "1"
    env.update(ED2_CUBIC_CACHE_PROOF_DUMP="1", OMP_NUM_THREADS=str(threads),
               OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               ED2_RNG_SEED="1984")
    return env


def normalized_input(raw):
    fixed = re.sub(rb"\bMNTE\s*=\s*(?:-1|64)\b",
                   b"MNTE=SAFE_AUTO_VS_EXPLICIT", raw)
    assert fixed != raw, "Input is missing MNTE selector"
    return fixed


def run_cache(driver, root, side, sector, threads, mnte, timeout):
    label = "auto" if mnte == -1 else "explicit64"
    work = root / "cache-n{}-{}-t{}-{}".format(side, sector, threads, label)
    write_input(work, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=mnte)
    env = environment(threads)
    with (work / "input.dat").open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, cwd=str(work), env=env,
                            timeout=timeout, check=False)
    log = cp.stdout.decode("utf8", errors="replace")
    (work / "driver.log").write_text(log)
    assert cp.returncode == 0, (work, cp.returncode, log[-3000:])
    assert (MARKER in log) == (mnte == -1), (work, "Auto-policy marker mismatch")
    assert re.search(r"Current MNTE\s*=\s*64\b", log), (work, "Capacity not 64")
    assert "Cubic cache builder: local-bond-sparse" in log
    assert "Representative lookup: cubic-direct-hop-pairs" in log
    assert "Basis builder: anchored-cubic" in log
    assert "Basis canonical enumeration: exact-axis-pruned" in log
    assert "Cubic cache proof dump: projected rows={}".format(
        sector_dimensions(side)[sector]) in log
    return work


def ritz(exe, root, side, sector, mnte, timeout):
    label = "auto" if mnte == -1 else "explicit64"
    work = root / "ritz-n{}-{}-{}".format(side, sector, label)
    write_input(work, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=mnte)
    env = environment(1)
    old_env = dict(os.environ)
    try:
        for k, v in env.items():
            if k.startswith("ED2_CUBIC_") or k in (
                    "ED2_RNG_SEED", "OMP_NUM_THREADS", "OMP_DYNAMIC",
                    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
                os.environ[k] = v
        os.environ.pop("ED2_CUBIC_CACHE_PROOF_DUMP", None)
        obs = run_case(exe, work, side ** 3, 3, sector,
                       sector_dimensions(side)[sector],
                       1, "two-stage", timeout)
    finally:
        os.environ.clear()
        os.environ.update(old_env)
    assert obs["residual"] <= 1.01e-8, obs
    return obs


def reject_invalid_geometry(driver, root, timeout):
    # The CPU cheap N64 case is cubic in symmetry but not the required
    # local 18-degree, 9N-bond Hamiltonian if one interaction is dropped.
    work = root / "reject-incorrect-9n-bond-count"
    write_input(work, 4, 3, symmetry_generators(4, "oh-a1g"),
                canonical_bonds(4), mnte=-1)
    p = work / "input.dat"
    original = p.read_text()
    assert "NO_two=576" in original
    p.write_text(original.replace("NO_two=576", "NO_two=575"))
    env = environment(1)
    with p.open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            cwd=str(work), env=env, timeout=timeout,
                            check=False)
    assert cp.returncode != 0, "Auto capacity accepted invalid 9N model"
    assert b"Auto cubic MNTE requires full 9*N cubic bonds and basis" in cp.stdout, (
        cp.stdout[-1800:])

    # Reject an incorrect particle sector, rather than silently shrinking a
    # general Hamiltonian using a three-magnon-only capacity formula.
    work = root / "reject-not-three-magnons"
    write_input(work, 4, 3, symmetry_generators(4, "oh-a1g"),
                canonical_bonds(4), mnte=-1)
    p = work / "input.dat"
    orig = p.read_text()
    assert "NODmin=3,NODmax=3" in orig
    p.write_text(orig.replace("NODmin=3,NODmax=3", "NODmin=2,NODmax=2"))
    with p.open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, cwd=str(work),
                            env=env, timeout=timeout, check=False)
    assert cp.returncode != 0, "Auto capacity accepted other magnon count"
    # The exact-cubic anchored basis guard may reject before the cache guard.
    assert b"error" in cp.stdout.lower() or b"requires" in cp.stdout.lower(), (
        cp.stdout[-1800:])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--driver", required=True, type=Path)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--sides", default="4,6,10")
    ap.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--timeout", type=int, default=1200)
    ap.add_argument("--skip-lanczos", action="store_true")
    ap.add_argument("--output", type=Path)
    a = ap.parse_args()
    driver, exe = a.driver.resolve(strict=True), a.exe.resolve(strict=True)
    root = a.output or Path(tempfile.mkdtemp(prefix="ed2-auto-mnte-"))
    if a.output:
        root.mkdir(parents=True, exist_ok=False)
    sides = [int(x) for x in a.sides.split(",")]
    sectors = a.sectors.split(",")
    threads = [int(x) for x in a.threads.split(",")]
    assert all(s in (4, 6, 8, 10) for s in sides)
    assert sectors and set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}
    assert threads and min(threads) > 0
    rows = []
    for side, sector, nthreads in itertools.product(sides, sectors, threads):
        assert check_hamiltonian(side)
        auto = run_cache(driver, root, side, sector, nthreads, -1, a.timeout)
        explicit = run_cache(driver, root, side, sector, nthreads, 64, a.timeout)
        assert normalized_input((auto / "input.dat").read_bytes()) == \
               normalized_input((explicit / "input.dat").read_bytes())
        assert sha(auto / "bonds.dat") == sha(explicit / "bonds.dat")
        assert sha(auto / "cubic_cache_proof_dump.txt") == \
               sha(explicit / "cubic_cache_proof_dump.txt"), (
                   "Auto capacity changed native Hamiltonian cache bits",
                   side, sector, nthreads)
        row = dict(side=side, sector=sector, threads=nthreads,
                   dimension=sector_dimensions(side)[sector],
                   cache_sha256=sha(auto / "cubic_cache_proof_dump.txt"))
        rows.append(row)
        print("CUBIC_AUTO_MNTE_FULL_CACHE_BITS_PASS side={} sector={} t={} D={}".format(
            side, sector, nthreads, row["dimension"]), flush=True)
        if not a.skip_lanczos and nthreads == 1 and side in (4, 6, 10) and sector == "oh-a1g":
            ra = ritz(exe, root, side, sector, -1, a.timeout)
            rb = ritz(exe, root, side, sector, 64, a.timeout)
            assert abs(ra["energy"] - rb["energy"]) < 1e-10
            assert ra["h_applications"] == rb["h_applications"]
            print("CUBIC_AUTO_MNTE_RITZ_PASS side={} E={:.13f} residual={:.3e}".format(
                side, ra["energy"], ra["residual"]), flush=True)
    reject_invalid_geometry(driver, root, min(a.timeout, 60))
    (root / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("CUBIC_AUTO_MNTE_REGRESSION_PASS cache_pairs={} negative=2".format(
        len(rows)), flush=True)


if __name__ == "__main__":
    main()
