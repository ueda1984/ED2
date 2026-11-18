#!/usr/bin/env python3
"""Compare *exact* list_s/list_r bits under generic and cubic group certificates.

The fast certificate uses the known Z_L^3 semidirect O_h group only when
explicitly opted in. A single native input and the identical solver/basis
builder must produce the SAME sorted representative rank and SAME floating
point bits for each normalization. This is stronger than just matching
energy/residual; every basis member is compared byte-for-byte.

Python 3.6+, stdlib. Independent of the fast validator's own implementation.
"""
import argparse
import json
import os
from pathlib import Path
import re
import tempfile

from cpc_cubic_pointgroup import (SECTORS, canonical_bonds,
                                  sector_dimensions, symmetry_generators, write_input)
from cpc_pointgroup_profile import run_case

FAST_FLAG = "ED2_CUBIC_FAST_GROUP_VALIDATION"
DUMP_FLAG = "ED2_CUBIC_BASIS_PROOF_DUMP"
CERT_LINE = "Symmetry certificate: canonical cubic O_h semidirect product verified"
PHASE_RE = re.compile(r"ED2_BASIS_PHASE ([a-z_]+) wall_s=\s*([+\-0-9.eE]+)")


def compare():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sides", default="4,6")
    parser.add_argument("--threads", default="1,4")
    parser.add_argument("--sectors", default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    sides = [int(s) for s in args.sides.split(",")]
    threads = [int(s) for s in args.threads.split(",")]
    sectors = args.sectors.split(",")
    assert sides and threads and sectors
    assert set(sectors) <= set(SECTORS[1:])
    exe = args.exe.resolve(strict=True)
    root = args.output or Path(tempfile.mkdtemp(prefix="ed2-cubic-group-certificate-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    report = []
    for side in sides:
        bonds = canonical_bonds(side)
        dimensions = sector_dimensions(side)
        for sector in sectors:
            for nt in threads:
                runs = []
                for fast in (False, True):
                    work = root / "s{}-{}-t{}-{}".format(
                        side, sector, nt, "new" if fast else "old")
                    write_input(work, side, 3,
                                symmetry_generators(side, sector), bonds, mnte=64)
                    changes = {
                        "ED2_CUBIC_ANCHORED_REP": "1",
                        "ED2_ALLOW_LARGE_GROUP_VALIDATION": "1",
                        "ED2_BASIS_DETAIL_TIMING": "1",
                        DUMP_FLAG: "1",
                        FAST_FLAG: "1" if fast else None,
                    }
                    originals = {key: os.environ.get(key) for key in changes}
                    try:
                        for key, value in changes.items():
                            if value is None:
                                os.environ.pop(key, None)
                            else:
                                os.environ[key] = value
                        obs = run_case(exe, work, side ** 3, 3, sector,
                                       dimensions[sector], nt, "two-stage", args.timeout)
                    finally:
                        for key, value in originals.items():
                            if value is None:
                                os.environ.pop(key, None)
                            else:
                                os.environ[key] = value
                    output = (work / "stdout.log").read_text(errors="replace")
                    assert (CERT_LINE in output) is fast, (side, sector, nt, fast)
                    assert "Symmetry action: noncommuting normal form and character validated" in output
                    assert "Basis proof dump: exactly ranked states and bitwise norms=" in output
                    assert "Representative lookup: cubic-anchored-orbits" in output
                    dump_path = work / "cubic_basis_proof_dump.txt"
                    proof = dump_path.read_bytes()
                    assert len(proof.splitlines()) == dimensions[sector], (work, dimensions[sector])
                    timings = dict((key, float(value)) for key, value in PHASE_RE.findall(output))
                    assert "group_validation" in timings
                    runs.append((obs, proof, timings))
                old, new = runs
                assert old[0]["input_sha256"] == new[0]["input_sha256"]
                assert old[0]["dimension"] == new[0]["dimension"] == dimensions[sector]
                assert old[1] == new[1], (
                    "EXACT list_s/list_r BITWISE MISMATCH", side, sector, nt)
                assert abs(old[0]["energy"]-new[0]["energy"]) < 1e-8
                assert old[0]["residual"] <= 1.01e-8 and new[0]["residual"] <= 1.01e-8
                if nt == 1:
                    assert old[0]["h_applications"] == new[0]["h_applications"]
                data = dict(side=side, sector=sector, threads=nt,
                            dimension=dimensions[sector],
                            energy=new[0]["energy"],
                            residual=new[0]["residual"],
                            old_validation_s=old[2]["group_validation"],
                            new_validation_s=new[2]["group_validation"],
                            exact_basis_and_norm_bits=True)
                report.append(data)
                print("EXACT_CUBIC_CERT_PASS side={} sector={} threads={} "
                      "D={} validation_old={:.6f}s validation_new={:.6f}s".format(
                          side, sector, nt, dimensions[sector],
                          data["old_validation_s"], data["new_validation_s"]), flush=True)
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print("CUBIC_OH_FAST_CERTIFICATE_REGRESSION_PASS {} bitwise matched pairs".format(
        len(report)), flush=True)


if __name__ == "__main__":
    compare()
