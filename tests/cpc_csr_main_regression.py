#!/usr/bin/env python3
"""Standalone regression for opt-in complex CSR H*x in ORIGINAL ED2 main.

No cubic-specific modules, newer Lanczos features, or scientific Python
packages are required. Compares COMPLETE native H*x bitstreams for two
nonreal-complex trial vectors across distinct model sizes, OpenMP counts,
and actual real/complex non-U(1) mixed-spin Hamiltonians.
Python 3.8+, standard library only.
"""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import subprocess

FLAG = "ED2_CSR_COMPLEX_CACHE"
MARKER = "Hamiltonian action: csr-complex cache entries="
COLUMNS = ("input.dat", "spins.dat", "caps.dat", "bonds.dat", "fields.dat")


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def deck(work, n, dm, partial=False):
    """Physical OBC alternating S=1/2,S=1, XYZ+hx, full Hilbert space."""
    assert n in (4, 6, 8, 10) and n % 2 == 0
    work.mkdir(parents=True, exist_ok=False)
    spins = [0.5, 1.0] * (n // 2)
    (work / "spins.dat").write_text(
        "".join("{:.1f}\n".format(s) for s in spins), encoding="utf-8")
    (work / "caps.dat").write_text(
        "".join("{}\n".format(int(2 * s)) for s in spins), encoding="utf-8")
    bonds = []
    for j in range(n - 1):
        sign = -1 if j % 2 else 1
        values = [
            0.77 + 0.05 * sign, 0.63 + 0.03 * sign,
            1.09 - 0.04 * sign,
            0.0, 0.0, dm * sign, 0.0, 0.0, 0.0
        ]
        bonds.append("{} {} {}\n".format(
            j + 1, j + 2, " ".join("{:.17g}".format(x) for x in values)))
    (work / "bonds.dat").write_text("".join(bonds), encoding="utf-8")
    fields = []
    for j in range(n):
        sign = -1 if j % 2 else 1
        # hx and Jxx!=Jyy ensure no U(1) conservation even for dm=0.
        vec = [0.15 * (0.4 * sign + 0.1), 0.0, 0.17 + 0.035 * sign]
        fields.append("{} {}\n".format(
            j + 1, " ".join("{:.17g}".format(x) for x in vec)))
    (work / "fields.dat").write_text("".join(fields), encoding="utf-8")
    caps = sum(int(2 * x) for x in spins)
    # ALG=2 uses the repository's original thick-restart Lanczos for the
    # optional whole-solver comparison. The H*x probe skips this solver.
    text = """&input_parameters
NOS={n}, NODmin=0, NODmax={caps}, ALG=2,
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,
M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={n},NO_two={num_bonds},wk_dim={wk_dim},MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0
/
&input_lancz
lnc_ene_conv=1d-12,minitr=5,maxitr=320,itrint=10
/
&input_TRLan
NOE=2,NOK=4,NOM=24,maxitr=320,lnc_ene_conv=1d-12,
i_vec_min=1,i_vec_max=2
/
""".format(n=n, caps=caps, num_bonds=n - 1,
           wk_dim=1 if partial else -1)
    (work / "input.dat").write_text(text, encoding="utf-8")


def env_for(t, flag):
    env = dict(os.environ)
    for name in list(env):
        if name.startswith("ED2_CUBIC_") or name in (
                FLAG, "ED2_REAL_CACHED_MATVEC",
                "ED2_CSR_REAL_CACHE", "ED2_PACKED_REAL_CACHE"):
            env.pop(name, None)
    env["OMP_NUM_THREADS"] = str(t)
    env["OMP_DYNAMIC"] = "FALSE"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    if flag is not None:
        env[FLAG] = flag
    return env


def native(exe, work, t, flag, timeout):
    with (work / "input.dat").open("rb") as inp:
        run = subprocess.run(
            [str(exe)], stdin=inp, cwd=str(work),
            env=env_for(t, flag), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=timeout, check=False)
    log = run.stdout.decode("utf-8", errors="replace")
    (work / "stdout.log").write_text(log, encoding="utf-8")
    return run.returncode, log


def matching_inputs(work):
    return tuple(sha256(work / name) for name in COLUMNS)


def check_vectors(driver, root, n, dm, t, timeout):
    cases = []
    for enabled in (False, True):
        work = root / "s{}-dm{}-t{}-csr{}".format(n, int(dm != 0), t, int(enabled))
        deck(work, n, dm)
        status, out = native(driver, work, t, "1" if enabled else "0", timeout)
        assert status == 0, (work, "native failure", out[-2500:])
        assert (MARKER in out) == enabled, (work, "native mode mismatch")
        expected_dim = 6 ** (n // 2)
        assert "CPC_CSR_H_VECTOR_DIM={}".format(expected_dim) in out
        assert "CPC_CSR_H_VECTOR_COUNT=2" in out
        assert out.count("CPC_CSR_DENSE_REFERENCE_PASS") == (
            2 if n<=6 else 0), (work, "dense H reference not certified")
        imag_match = re.search(r"CPC_CSR_NATIVE_IMAG_COUNT=(\d+)", out)
        assert imag_match is not None
        imag_count = int(imag_match.group(1))
        assert (imag_count > 0) == (dm != 0), (
            "Test must cover physically complex Hamiltonian", work, imag_count)
        outfile = work / "csr_matvec_bits.txt"
        assert len(outfile.read_bytes().splitlines()) == 2 * expected_dim
        if enabled:
            found = re.search(
                r"csr-complex cache entries=(\d+), extra_bytes=(\d+)", out)
            assert found, (work, "missing exact memory report")
            entries, nbytes = map(int, found.groups())
            assert entries >= expected_dim
            assert nbytes == 20 * entries + 4 * (expected_dim + 1)
        cases.append((matching_inputs(work), sha256(outfile)))
    assert cases[0] == cases[1], (
        "Full complex H-vector bitwise mismatch", n, dm, t, cases)
    print("MAIN_COMPLEX_CSR_FULL_HV_BITS_PASS size={} dm={} threads={}".format(
        n, dm, t), flush=True)



def reject(driver, root, name, n, dm, flag, partial, substring):
    work = root / ("reject-" + name)
    deck(work, n, dm, partial=partial)
    code, stdout = native(driver, work, 1, flag, 60)
    assert code != 0 and substring in stdout, (
        "Expected fail-closed check did not trigger", name, code, stdout[-2500:])
    print("MAIN_COMPLEX_CSR_INVALID_SETTING_REJECTED " + name, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", type=Path, required=True)

    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", default="4,6,8")
    parser.add_argument("--threads", default="1,4")
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()
    driver = args.driver.resolve(strict=True)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    sizes = [int(x) for x in args.sizes.split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    assert sizes and all(x in (4, 6, 8, 10) for x in sizes)
    assert threads and all(x > 0 for x in threads)

    for n, dm, t in itertools.product(sizes, (0.0, 0.12), threads):
        check_vectors(driver,root,n,dm,t,args.timeout)
    # Full H*x proof for N4/N6 is also independently checked against
    # the separately built full dense Hamiltonian in the native driver.
    reject(driver,root,"partial",4,0.12,"1",True,
           "ED2_CSR_COMPLEX_CACHE requires fully cached matrix")
    reject(driver,root,"invalid",4,0.12,"2",False,
           "ED2_CSR_COMPLEX_CACHE must be 0 or 1")
    reject(driver,root,"oversized",4,0.12,"11111111111111111111",
           False,"ED2_CSR_COMPLEX_CACHE flag too long or invalid")
    report = {
        "source": "standalone original-main branch independent scientific test",
        "vector_pairs": len(sizes)*2*len(threads),
        "dense_reference_model_cases": 2*len(threads)*sum(n<=6 for n in sizes),
        "invalid_cases": 3,
        "no_presto": True,
        "mode": "default-off / ED2_CSR_COMPLEX_CACHE=1",
        "status": "PASS"
    }
    (root / "summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print("MAIN_COMPLEX_CSR_SCIENCE_PASS "
          "vector_pairs={} dense_reference_cases={} rejections=3".format(
              report["vector_pairs"],report["dense_reference_model_cases"]),flush=True)


if __name__ == "__main__":
    main()
