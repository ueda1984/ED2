#!/usr/bin/env python3
"""Native complex CSR cache: exact H*x bit proof and true-complex DM tests.

The three-magnon 1D Oh cases verify all 28 output vector bit hashes for a
deterministic complex trial vector. Generic alternating mixed-spin
Hamiltonians with strictly real and genuinely complex DM entries are also
checked from the actual Fortran matrix builder, then full native Lanczos.
The optional CSR mode retains every cached row term in its historical order.
Python >=3.6, standard library only.
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
    cubic_env, generic_env, run, sha)

FLAG = "ED2_CSR_COMPLEX_CACHE"
MARKER = "Hamiltonian action: csr-complex cache entries="


def cubic_mode(threads, enabled):
    env = cubic_env(threads, False)
    env.pop(FLAG, None)
    env.pop("ED2_REAL_CACHED_MATVEC", None)
    env["ED2_CUBIC_MATVEC_PROOF_DUMP"] = "1"
    if enabled:
        env[FLAG] = "1"
    return env


def generic_mode(threads, enabled, proof):
    env = generic_env(None)
    env.pop(FLAG, None)
    env.pop("ED2_REAL_CACHED_MATVEC", None)
    env.update(OMP_NUM_THREADS=str(threads), ED2_LANCZOS_TWO_STAGE="1",
               ED2_RNG_SEED="1984", ED2_RESIDUAL_ATOL="1e-8")
    if proof:
        env["ED2_CUBIC_MATVEC_PROOF_DUMP"] = "1"
    if enabled:
        env[FLAG] = "1"
    return env


def make_generic(work, sites, complex_dm, partial=False):
    write_dm_input(work, sites, 0.12 if complex_dm else 0.0)
    if not complex_dm:
        lines = []
        for line in (work / "bonds.dat").read_text().splitlines():
            fields = line.split()
            fields[-1] = "0"
            lines.append(" ".join(fields))
        (work / "bonds.dat").write_text("\n".join(lines) + "\n")
    if partial:
        deck = (work / "input.dat").read_text()
        assert deck.count("wk_dim=-1") == 1
        (work / "input.dat").write_text(deck.replace("wk_dim=-1", "wk_dim=1"))


def reject(exe, root, name, opts, expect, sites=4, complex_dm=False, partial=False):
    work = root / ("reject-" + name)
    make_generic(work, sites, complex_dm, partial=partial)
    env = generic_mode(1, False, False)
    env.update(opts)
    code, log = run(exe, work, env, 90)
    assert code != 0 and expect in log, ("Missing fail-closed check", name, log[-2000:])
    print("CSR_COMPLEX_BAD_CONFIGURATION_REJECTED " + name, flush=True)


def generic_tests(driver, exe, root, timeout):
    count = 0
    for sites, complex_dm, threads in itertools.product((4, 6), (False, True), (1, 4)):
        hashes = []
        for enabled in (False, True):
            name = "s{}-dm{}-t{}-csr{}".format(sites, int(complex_dm),
                                               threads, int(enabled))
            work = root / ("generic-" + name)
            make_generic(work, sites, complex_dm)
            code, log = run(driver, work, generic_mode(threads, enabled, True), timeout)
            assert code == 0, (name, log[-2300:])
            assert (MARKER in log) is enabled
            dim = re.findall(r"ED2_CUBIC_MATVEC_PROOF_ROWS=(\d+)", log)
            assert len(dim) == 1
            proof = work / "cubic_matvec_proof_dump.txt"
            assert len(proof.read_bytes().splitlines()) == int(dim[0])
            hashes.append((sha(work / "input.dat"), sha(proof)))
        assert hashes[0] == hashes[1], ("Generic mixed-spin DM H*x bits differ", sites,
                                       complex_dm, threads)
        count += 1
        print("CSR_COMPLEX_GENERIC_NATIVE_BITS_PASS sites={} complex_dm={} threads={}"
              .format(sites, complex_dm, threads), flush=True)

    for dm in (False, True):
        records = []
        for enabled in (False, True):
            label = "dm{}-csr{}".format(int(dm), int(enabled))
            work = root / ("generic-eigen-" + label)
            make_generic(work, 6, dm)
            code, log = run(exe, work, generic_mode(1, enabled, False), timeout)
            assert code == 0, (label, log[-2300:])
            assert (MARKER in log) is enabled
            # This generic DM input selects the existing thick-restarted
            # TRLan solver, not the cubic two-stage Lanczos pathway.
            assert "TRLan status:" in log and "residuals certified" in log, (
                label, "Native thick-restart solver did not certify residuals")
            matches = re.findall(
                r"Eigenpair residual norms:\s*([0-9.Ee+\-]+)\s+([0-9.Ee+\-]+)", log)
            assert len(matches) == 1 and all(
                float(val) <= 1.01e-8 for val in matches[0])
            records.append((sha(work / "input.dat"),
                            sha(work / "eigenvalues.dat")))
        assert records[0] == records[1], ("Mixed-spin energy bits changed", dm)
    print("CSR_COMPLEX_GENERIC_RESTARTED_RITZ_BITS_PASS 2 complex/real models", flush=True)

    reject(exe, root, "partial", {FLAG: "1"},
           "ED2_CSR_COMPLEX_CACHE requires fully cached matrix", partial=True)
    reject(exe, root, "invalid", {FLAG: "2"},
           "ED2_CSR_COMPLEX_CACHE must be 0 or 1")
    reject(exe, root, "oversized", {FLAG: "11111111111111111111"},
           "ED2_CSR_COMPLEX_CACHE flag too long or invalid")
    reject(exe, root, "conflict-real", {FLAG: "1", "ED2_REAL_CACHED_MATVEC": "1"},
           "ED2_CSR_COMPLEX_CACHE incompatible with ED2_REAL_CACHED_MATVEC=1")
    print("CSR_COMPLEX_GENERIC_REGRESSION_PASS pairs={} negative_cases=4"
          .format(count), flush=True)


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
    driver, exe = args.driver.resolve(strict=True), args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-csr-complex-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    sides = [int(x) for x in args.sides.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert sides and all(x in (4, 6, 8, 10) for x in sides)
    assert threads and all(t >= 1 for t in threads)
    assert sectors and set(sectors) <= {"oh-a1g", "oh-a2g", "oh-a1u", "oh-a2u"}
    count = 0
    for side, sector, t in itertools.product(sides, sectors, threads):
        hashes = []
        for enabled in (False, True):
            dest = root / "s{}-{}-t{}-csr{}".format(side, sector, t, int(enabled))
            write_input(dest, side, 3, symmetry_generators(side, sector),
                        canonical_bonds(side), mnte=64)
            code, log = run(driver, dest, cubic_mode(t, enabled), args.timeout)
            assert code == 0, (dest, log[-2200:])
            assert (MARKER in log) == enabled
            assert "ED2_CUBIC_MATVEC_PROOF_ROWS={}".format(
                sector_dimensions(side)[sector]) in log
            dump = dest / "cubic_matvec_proof_dump.txt"
            assert len(dump.read_bytes().splitlines()) == sector_dimensions(side)[sector]
            hashes.append((sha(dest / "input.dat"), sha(dump)))
            if enabled:
                found = re.search(r"csr-complex cache entries=(\d+), extra_bytes=(\d+)", log)
                assert found, (dest, "missing numeric memory provenance")
                entries, bytes_extra = map(int, found.groups())
                assert entries >= sector_dimensions(side)[sector]
                assert bytes_extra == 20*entries + 4*(sector_dimensions(side)[sector]+1)
        assert hashes[0] == hashes[1], ("Cubic exact complex H*x mismatch", side, sector, t)
        count += 1
        print("CSR_COMPLEX_CUBIC_MATVEC_BITS_PASS side={} sector={} t={}"
              .format(side, sector, t), flush=True)

    if not args.skip_generic:
        generic_tests(driver, exe, root, args.timeout)
    print("CSR_COMPLEX_CACHE_REGRESSION_PASS cubic_pairs={} generic_tested={}"
          .format(count, not args.skip_generic), flush=True)


if __name__ == "__main__":
    main()
