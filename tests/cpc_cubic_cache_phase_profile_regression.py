#!/usr/bin/env python3
"""Prove optional sampled cubic-cache profiling changes no cache bit pattern.

Only a test-only native cache driver is invoked.  The production cache
is created twice from byte-identical inputs, and every source-column
diagonal/offdiagonal/zero entry is checked with a full SHA256 equality of
the native IEEE-754 bit dump. Timing attribution is *not* tested against
hardcoded wall times on a shared CI machine.
"""
import argparse
import hashlib
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from cpc_cubic_pointgroup import (
    canonical_bonds, sector_dimensions, symmetry_generators, write_input)


def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as src:
        for block in iter(lambda: src.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def produce(driver, root, side, sector, threads, profiling, timeout):
    dest = root / ("n{}-{}-t{}-{}".format(
        side, sector, threads, "profiled" if profiling else "plain"))
    write_input(dest, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=64)
    env = dict(os.environ)
    env.update(ED2_CUBIC_ANCHORED_REP="1",
               ED2_CUBIC_FAST_GROUP_VALIDATION="1",
               ED2_CUBIC_ANCHORED_BASIS="1",
               ED2_CUBIC_LOCAL_BOND_CACHE="1",
               ED2_CUBIC_CACHE_PROOF_DUMP="1",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS=str(threads),
               OMP_DYNAMIC="FALSE")
    env.pop("ED2_CUBIC_CACHE_PROFILE", None)
    env.pop("ED2_FORCE_LEGACY_ORBITS", None)
    if profiling:
        env["ED2_CUBIC_CACHE_PROFILE"] = "1"
    with (dest / "input.dat").open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            cwd=str(dest), env=env, timeout=timeout,
                            check=False)
    log = cp.stdout.decode("utf-8", errors="replace")
    (dest / "driver.log").write_text(log)
    assert cp.returncode == 0, (dest, cp.returncode, log[-3000:])
    assert "Cubic cache builder: local-bond-sparse" in log, (dest, log[-1000:])
    assert ("Cubic cache profile: representative_key_setup_wall_s=" in log) is profiling
    assert ("Cubic cache profile: resolver_calls=" in log) is profiling
    assert "Basis builder: anchored-cubic" in log
    dim = sector_dimensions(side)[sector]
    assert "Cubic cache proof dump: projected rows={}".format(dim) in log

    if profiling:
        def metric(key):
            prefix = (r", " if key == "omp_max_threads" else
                      r"Cubic cache profile: ")
            r = re.search(
                prefix + re.escape(key) +
                r"=\s*([+-.0-9eEdD]+)", log)
            assert r, (key, dest)
            x = float(r.group(1).replace("D", "E"))
            assert math.isfinite(x) and x >= 0, (key, x)
            return x
        assert metric("columns") == dim
        assert metric("omp_max_threads") >= threads
        assert 0 < metric("unique_incident_bonds") <= 54 * dim
        assert 0 < metric("active_exchange_bonds") <= metric("unique_incident_bonds")
        assert 0 < metric("stored_hops") <= metric("resolver_calls")
        assert 0 < metric("resolver_sampled_calls") <= metric("resolver_calls")
        assert metric("adjacency_setup_wall_s") >= 0
        assert metric("representative_key_setup_wall_s") >= 0
        assert metric("state_decode_thread_elapsed_sum_s") >= 0
        assert metric("bond_union_sort_thread_elapsed_sum_s") >= 0
        assert metric("bond_math_thread_elapsed_sum_s") >= 0
        assert metric("resolver_sampled_thread_elapsed_sum_s") >= 0
        assert metric("resolver_extrapolated_thread_elapsed_sum_s") >= 0
        assert "sampling bias prohibit publication timing" in log
    return dest


def reject_invalid_profile(driver, root):
    side = 4
    sector = "oh-a1g"
    dest = root / "reject-profile-no-local-cache"
    write_input(dest, side, 3, symmetry_generators(side, sector),
                canonical_bonds(side), mnte=64)
    env = dict(os.environ)
    env.update(ED2_CUBIC_ANCHORED_REP="1",
               ED2_CUBIC_FAST_GROUP_VALIDATION="1",
               ED2_CUBIC_ANCHORED_BASIS="1",
               ED2_CUBIC_CACHE_PROFILE="1",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    env.pop("ED2_CUBIC_LOCAL_BOND_CACHE", None)
    with (dest / "input.dat").open("rb") as inp:
        cp = subprocess.run([str(driver)], stdin=inp,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            cwd=str(dest), env=env, timeout=60, check=False)
    assert cp.returncode != 0
    assert b"Cubic cache profile requires ED2_CUBIC_LOCAL_BOND_CACHE=1" in cp.stdout
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--sides", default="4,6,10")
    p.add_argument("--threads", default="1,4")
    p.add_argument("--sectors", default="oh-a1g")
    p.add_argument("--output", type=Path)
    p.add_argument("--timeout", type=int, default=900)
    a = p.parse_args()
    driver = a.driver.resolve(strict=True)
    root = a.output or Path(tempfile.mkdtemp(prefix="ed2-cache-profile-"))
    if a.output:
        root.mkdir(parents=True, exist_ok=False)
    sectors = tuple(a.sectors.split(","))
    assert set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}
    count = 0
    for side in [int(z) for z in a.sides.split(",")]:
        for sector in sectors:
            for threads in [int(z) for z in a.threads.split(",")]:
                plain = produce(driver, root, side, sector, threads, False, a.timeout)
                prof = produce(driver, root, side, sector, threads, True, a.timeout)
                assert sha256(plain / "input.dat") == sha256(prof / "input.dat")
                assert sha256(plain / "cubic_cache_proof_dump.txt") == sha256(
                    prof / "cubic_cache_proof_dump.txt"), (
                    "Profiling changed native Hamiltonian cache bits",
                    side, sector, threads)
                count += 1
                print("CUBIC_CACHE_PHASE_PROFILE_BITS_PASS side={} sector={} threads={} "
                      "D={}".format(side, sector, threads,
                                    sector_dimensions(side)[sector]), flush=True)
    reject_invalid_profile(driver, root)
    print("CUBIC_CACHE_PHASE_PROFILE_REGRESSION_PASS {} byte-identical "
          "profiled/unprofiled complete cache pairs; fail-closed negative PASS".format(count))


if __name__ == "__main__":
    main()
