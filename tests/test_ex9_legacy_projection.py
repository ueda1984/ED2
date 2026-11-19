#!/usr/bin/env python3
"""Unit-test the historical ex9 ladder-symmetry projection without new compute.

The mock 'actual' tables are derived from the archived tables and must never
be confused with a genuine ex9 recomputation. Full-run evidence is separate.
"""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
original = root / "examples/ref_data/cubic_sp_HB"
validator = root / "tests/verify_ex9.py"
files = ("eigenvalues.dat", "local_mag.dat", "two_body_cf_xyz.dat",
         "two_body_cf_z+-.dat")


def run(folder, projected, reference=None):
    command = [sys.executable, str(validator), "--actual", str(folder)]
    if projected:
        command.append("--legacy-ex9-symmetry-projection")
    if reference is not None:
        command.extend(["--reference", str(reference)])
    return subprocess.run(command, capture_output=True, text=True)


def copy_tables(source, destination):
    destination.mkdir()
    for name in files:
        shutil.copyfile(source / name, destination / name)


def project_mock(folder):
    ladder = folder / "two_body_cf_z+-.dat"
    out = []
    for line in ladder.read_text().splitlines():
        columns = line.split()
        real = (float(columns[12]) + float(columns[16])) / 2
        columns[12] = str(real)
        columns[16] = str(real)
        columns[13] = "0.0"
        columns[17] = "0.0"
        out.append(" ".join(columns))
    ladder.write_text("\n".join(out) + "\n")


def expect_identity_rejection(result):
    assert result.returncode != 0, "unrecognized projection reference accepted"
    assert "not the pinned historical ex9 reference" in result.stdout, result.stdout + result.stderr


with tempfile.TemporaryDirectory(prefix="ed2_ex9_projection_test_") as temp:
    target = Path(temp)
    for name in files:
        shutil.copyfile(original / name, target / name)
    ladder = target / "two_body_cf_z+-.dat"
    project_mock(target)
    assert run(target, projected=False).returncode != 0, "raw archive incorrectly accepted"
    result = run(target, projected=True)
    assert result.returncode == 0, result.stdout + result.stderr

    # Identity is based on contents, so an unchanged copy in another folder works.
    copied = target / "copied-reference"
    copy_tables(original, copied)
    result = run(target, projected=True, reference=copied)
    assert result.returncode == 0, result.stdout + result.stderr

    # A self-consistent but different model with the same pair IDs must not be
    # eligible for this archive-specific correction. Adjust both representations
    # together so ladder/Cartesian consistency alone cannot reject the fixture.
    changed = target / "different-reference"
    copy_tables(original, changed)
    for name, indices, delta in (("two_body_cf_z+-.dat", (12, 16), 1e-4),
                                 ("two_body_cf_xyz.dat", (2, 6), 5e-5)):
        path = changed / name
        rows = [line.split() for line in path.read_text().splitlines()]
        for index in indices:
            rows[0][index] = str(float(rows[0][index]) + delta)
        path.write_text("\n".join(" ".join(row) for row in rows) + "\n")
    changed_actual = target / "different-mock-actual"
    copy_tables(changed, changed_actual)
    project_mock(changed_actual)
    expect_identity_rejection(run(changed_actual, projected=True, reference=changed))

    # Site 12 is diagonally displaced in the documented 10x10x10 numbering.
    off_axis = target / "off-axis-reference"
    copy_tables(original, off_axis)
    for name in ("two_body_cf_xyz.dat", "two_body_cf_z+-.dat"):
        path = off_axis / name
        rows = [line.split() for line in path.read_text().splitlines()]
        rows[0][1] = "12"
        path.write_text("\n".join(" ".join(row) for row in rows) + "\n")
    expect_identity_rejection(run(target, projected=True, reference=off_axis))

    # Byte identity is deliberate, including otherwise harmless formatting.
    with (copied / "eigenvalues.dat").open("a") as stream:
        stream.write("\n")
    expect_identity_rejection(run(target, projected=True, reference=copied))

    # Perturb a physical ladder correlator without violating +-/-+ symmetry.
    bad = [line.split() for line in ladder.read_text().splitlines()]
    for index in (12, 16):
        bad[0][index] = str(float(bad[0][index]) + 1e-4)
    ladder.write_text("\n".join(" ".join(line) for line in bad) + "\n")
    assert run(target, projected=True).returncode != 0, "tampered data accepted"
print("EX9_LEGACY_HERMITIAN_PROJECTION_TEST_PASS")
