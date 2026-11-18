#!/usr/bin/env python3
"""Exact full native complex-vector H-action proof for CSR real cache.

Compares the PR22 direct-real multiply to a new CSR row list, retaining
all matrix terms (including explicit stored zeros) in original order.
The native Fortran cubic_cache_driver emits EVERY output IEEE-754 bit for
a deterministic nonreal complex input vector; no Python H matrix rewrite.
Python >=3.6 / stdlib only.
"""
import argparse
import itertools
import os
from pathlib import Path
import re
import tempfile

from cpc_cubic_pointgroup import (
    canonical_bonds, sector_dimensions, symmetry_generators, write_input)
from cpc_portable_dm_quickstart import write_input as write_dm_input
from cpc_real_cached_matvec_regression import (
    MARKER as REAL_MARKER, cubic_env, generic_env, run, sha)

FLAG = "ED2_CSR_REAL_CACHE"
MARKER = "Hamiltonian action: csr-real cache entries="


def cubic_mode(threads, enabled):
    env = cubic_env(threads, True)
    env.pop(FLAG, None)
    if enabled:
        env[FLAG] = "1"
    return env


def dm_case(root, exe, name, flag, real_matrix, partial, timeout):
    dest = root / name
    write_dm_input(dest, 4, 0.0 if real_matrix else 0.12)
    if real_matrix:
        lines = []
        for line in (dest / "bonds.dat").read_text().splitlines():
            fields = line.split()
            fields[-1] = "0"
            lines.append(" ".join(fields))
        (dest / "bonds.dat").write_text("\n".join(lines) + "\n")
    if partial:
        deck = (dest / "input.dat").read_text()
        assert deck.count("wk_dim=-1") == 1
        (dest / "input.dat").write_text(deck.replace("wk_dim=-1", "wk_dim=1"))
    env = generic_env("1" if flag is not None else None)
    env.pop(FLAG, None)
    if flag is not None:
        env[FLAG] = flag
    code, log = run(exe, dest, env, timeout)
    return dest, code, log


def test_generic(root, exe, timeout):
    before, code, log = dm_case(root, exe, "generic-real-off", "0", True, False, timeout)
    assert code == 0, log[-2200:]
    after, code, log = dm_case(root, exe, "generic-real-on", "1", True, False, timeout)
    assert code == 0, log[-2200:]
    assert REAL_MARKER in log and MARKER in log
    assert sha(before / "eigenvalues.dat") == sha(after / "eigenvalues.dat")
    print("CSR_REAL_GENERIC_MIXED_SPIN_EIGEN_BITS_PASS", flush=True)

    _, code, log = dm_case(root, exe, "reject-no-real-mode", None, True, False, timeout)
    assert code == 0 and MARKER not in log  # controls have no CSR opt-in
    dest = root / "reject-no-real-mode-request"
    write_dm_input(dest, 4, 0.0)
    env = generic_env(None)
    env[FLAG] = "1"
    code, log = run(exe, dest, env, timeout)
    assert code != 0 and "ED2_CSR_REAL_CACHE requires ED2_REAL_CACHED_MATVEC=1" in log
    print("CSR_REAL_MISSING_REAL_MODE_REJECTED", flush=True)

    _, code, log = dm_case(root, exe, "reject-nonreal-dm", "1", False, False, timeout)
    assert code != 0 and "ED2_REAL_CACHED_MATVEC requires exactly real matrix elements" in log
    print("CSR_REAL_NONREAL_DM_REJECTED", flush=True)

    _, code, log = dm_case(root, exe, "reject-partial-cache", "1", True, True, timeout)
    assert code != 0 and "ED2_REAL_CACHED_MATVEC requires fully cached matrix" in log
    print("CSR_REAL_PARTIAL_CACHE_REJECTED", flush=True)

    _, code, log = dm_case(root, exe, "reject-invalid-flag", "2", True, False, timeout)
    assert code != 0 and "ED2_CSR_REAL_CACHE must be 0 or 1" in log
    print("CSR_REAL_INVALID_FLAG_REJECTED", flush=True)

    _, code, log = dm_case(root, exe, "reject-oversized-flag", "11111111111111111111",
                           True, False, timeout)
    assert code != 0 and "ED2_CSR_REAL_CACHE flag too long or invalid" in log
    print("CSR_REAL_OVERSIZED_FLAG_REJECTED", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver", required=True, type=Path)
    p.add_argument("--exe", required=True, type=Path)
    p.add_argument("--sides", default="4,6,10")
    p.add_argument("--threads", default="1,4")
    p.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    p.add_argument("--timeout", default=600, type=int)
    p.add_argument("--skip-generic", action="store_true")
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    driver, exe = a.driver.resolve(strict=True), a.exe.resolve(strict=True)
    root = a.output or Path(tempfile.mkdtemp(prefix="ed2-csr-real-"))
    if a.output:
        root.mkdir(parents=True, exist_ok=False)
    sides = [int(x) for x in a.sides.split(",")]
    threads = [int(x) for x in a.threads.split(",")]
    sectors = a.sectors.split(",")
    assert sides and all(x in (4, 6, 8, 10) for x in sides)
    assert threads and all(x >= 1 for x in threads)
    assert sectors and set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}

    count = 0
    for side, sector, t in itertools.product(sides, sectors, threads):
        hashes = []
        for mode in (False, True):
            dest = root / "s{}-{}-t{}-csr{}".format(side, sector, t, int(mode))
            write_input(dest, side, 3, symmetry_generators(side, sector),
                        canonical_bonds(side), mnte=64)
            code, log = run(driver, dest, cubic_mode(t, mode), a.timeout)
            assert code == 0, (dest, log[-3000:])
            assert REAL_MARKER in log
            assert (MARKER in log) == mode
            assert "ED2_CUBIC_MATVEC_PROOF_ROWS={}".format(
                sector_dimensions(side)[sector]) in log
            proof = dest / "cubic_matvec_proof_dump.txt"
            assert len(proof.read_bytes().splitlines()) == sector_dimensions(side)[sector]
            hashes.append((sha(dest / "input.dat"), sha(proof)))
            if mode:
                line = re.search(r"csr-real cache entries=(\d+), extra_bytes=(\d+)", log)
                assert line, (dest, "missing exact memory report")
                nonzero, memory = map(int, line.groups())
                assert nonzero >= sector_dimensions(side)[sector]
                assert memory == 12*nonzero + 4*(sector_dimensions(side)[sector]+1)
        assert hashes[0] == hashes[1], ("CSR REAL H ACTION BIT DIFFERENCE", side, sector, t)
        count += 1
        print("CSR_REAL_NATIVE_MATVEC_BITS_PASS L={} irrep={} t={}".format(
            side, sector, t), flush=True)

    if not a.skip_generic:
        test_generic(root, exe, a.timeout)
    print("CSR_REAL_CACHE_REGRESSION_PASS exact_pairs={} generic_tested={}"
          .format(count, not a.skip_generic), flush=True)


if __name__ == "__main__":
    main()
