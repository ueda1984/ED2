#!/usr/bin/env python3
"""Verify optional CSR padded-cache release against full native H*x bits.

Requires only Python stdlib, original ED2 modules, and the native proof
driver. All three variants run byte-identical physical inputs:
  original: CSR off, padded cache retained;
  csr:      CSR on, padded cache retained;
  release:  CSR on, padded cache deallocated after copying.
For every model, compare both entire complex output vectors bit-for-bit.
"""
import argparse
import itertools
import json
import math
import struct
import os
from pathlib import Path
import re
import subprocess

from cpc_csr_main_regression import deck, matching_inputs, sha256

MODES = {
    "original": ("0", "0"),
    "csr": ("1", "0"),
    "release": ("1", "1"),
}
CSR_MARKER = "Hamiltonian action: csr-complex cache entries="
RELEASE_MARKER = "Hamiltonian action: csr-padded-cache released_bytes="


def run(driver, work, threads, csr, release, timeout):
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("ED2_CUBIC_") or key in (
            "ED2_CSR_COMPLEX_CACHE", "ED2_CSR_RELEASE_PADDED_CACHE",
            "ED2_CSR_REAL_CACHE", "ED2_PACKED_REAL_CACHE",
            "ED2_REAL_CACHED_MATVEC"):
            env.pop(key, None)
    env.update(
        ED2_CSR_COMPLEX_CACHE=csr,
        ED2_CSR_RELEASE_PADDED_CACHE=release,
        OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (work / "input.dat").open("rb") as source:
        proc = subprocess.run([str(driver)], cwd=work, stdin=source,
                              env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=timeout)
    log = proc.stdout.decode("utf-8", errors="replace")
    (work / "stdout.log").write_text(log, encoding="utf-8")
    return proc.returncode, log


def parse_complex_dump(data):
    out=[]
    for raw in data.decode("ascii").splitlines():
        fields=raw.split()
        assert len(fields)==4
        bits_re=int(fields[2],16).to_bytes(8,"big")
        bits_im=int(fields[3],16).to_bytes(8,"big")
        re_val=struct.unpack(">d",bits_re)[0]
        im_val=struct.unpack(">d",bits_im)[0]
        assert math.isfinite(re_val) and math.isfinite(im_val)
        out.append(complex(re_val,im_val))
    return out


def check_case(root, driver, n, dm, threads, timeout, allow_fp_reassociation=False):
    results = []
    bitstreams=[]
    rows = 6 ** (n//2)
    for mode, (csr, release) in MODES.items():
        work = root / ("s{}-dm{}-t{}-{}".format(
            n, int(dm != 0), threads, mode))
        deck(work, n, dm)
        code, log = run(driver, work, threads, csr, release, timeout)
        assert code == 0, (work, code, log[-3000:])
        assert (CSR_MARKER in log) == (mode != "original")
        assert (RELEASE_MARKER in log) == (mode == "release")
        assert ("CPC_CSR_PADDED_RELEASE_PASS" in log) == (mode == "release")
        assert ("CPC_CSR_PADDED_RETAINED_PASS" in log) == (mode != "release")
        assert out_rows(log) == rows
        assert log.count("CPC_CSR_DENSE_REFERENCE_PASS") == (2 if n<=6 else 0)
        imag = re.search(r"CPC_CSR_NATIVE_IMAG_COUNT=(\d+)", log)
        assert imag and (int(imag.group(1)) > 0) == (dm != 0)
        dump = work / "csr_matvec_bits.txt"
        bits=dump.read_bytes()
        assert len(bits.splitlines()) == 2*rows
        bitstreams.append(bits)
        data = {
            "mode": mode, "sha256": sha256(dump), "input_sha256": matching_inputs(work)
        }
        if mode != "original":
            m = re.search(
                r"csr-complex cache entries=(\d+), extra_bytes=(\d+)", log)
            assert m
            entries, extra = map(int, m.groups())
            assert entries >= rows
            assert extra == 20*entries+4*(rows+1)
            data.update(entries=entries, csr_extra_bytes=extra)
        if mode == "release":
            m = re.search(
                r"csr-padded-cache released_bytes=(\d+), net_saved_bytes=(-?\d+)", log)
            assert m
            released, saved = map(int, m.groups())
            assert saved == released-data["csr_extra_bytes"]
            assert released > 0
            data.update(released_bytes=released, net_saved_bytes=saved)
        results.append(data)
    assert len({x["input_sha256"] for x in results}) == 1
    # Both CSR storage variants run the identical multiplication loop even
    # under -Ofast and MUST retain exact original IEEE bit patterns.
    assert bitstreams[1] == bitstreams[2], (
        "Releasing padded data changed CSR H-vector bits",n,dm,threads)
    baseline_exact=(bitstreams[0] == bitstreams[1])
    if allow_fp_reassociation:
        old=parse_complex_dump(bitstreams[0])
        csr=parse_complex_dump(bitstreams[1])
        assert len(old)==len(csr)==2*rows
        max_delta=max(abs(x-y) for x,y in zip(old,csr))
        scale=max(1.0,max(abs(x) for x in old))
        assert max_delta <= 1.0e-12*scale, (
            "Fast-math changed H*x beyond double roundoff",n,dm,threads,
            max_delta,scale)
        results[0]["original_vs_csr_max_abs_delta"]=max_delta
        results[0]["original_vs_csr_max_scaled_error"]=max_delta/scale
        print("CSR_PADDED_RELEASE_FAST_MATH_NUMERIC_PASS n={} dm={} t={} max_delta={:.3e}".format(
            n,dm,threads,max_delta),flush=True)
    else:
        assert baseline_exact, (
            "Strict -O2 full complex H-vector bits changed",n,dm,threads,results)
        print("CSR_PADDED_RELEASE_ALL_THREE_NATIVE_BITS_PASS n={} dm={} threads={}".format(
            n, dm, threads), flush=True)
    results[0]["original_vs_csr_exact_bits"]=baseline_exact
    return results


def out_rows(log):
    m = re.search(r"CPC_CSR_H_VECTOR_DIM=(\d+)",log)
    assert m and "CPC_CSR_H_VECTOR_COUNT=2" in log
    return int(m.group(1))


def reject(root, driver, name, csr, release, partial, message):
    work = root / ("reject-"+name)
    deck(work, 4, 0.12, partial=partial)
    code, log = run(driver, work, 1, csr, release, 60)
    assert code != 0 and message in log, (name, code, log[-2300:])
    print("CSR_PADDED_RELEASE_INVALID_SETTING_REJECTED "+name, flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sizes", default="4,6,8")
    parser.add_argument("--threads", default="1,4")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--allow-fp-reassociation", action="store_true",
                        help="For -Ofast only: require numeric original/CSR equality; still require exact CSR/release bits")
    args=parser.parse_args()
    driver=args.driver.resolve(strict=True)
    sizes=[int(x) for x in args.sizes.split(",")]
    threads=[int(x) for x in args.threads.split(",")]
    assert sizes and all(n in (4,6,8,10) for n in sizes)
    assert threads and all(t>0 for t in threads)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    trials={}
    for n, dm, t in itertools.product(sizes, (0.0,0.12), threads):
        trials["n{}-dm{}-t{}".format(n,dm,t)]=check_case(
            root,driver,n,dm,t,args.timeout,args.allow_fp_reassociation)
    reject(root,driver,"missing-csr","0","1",False,
           "ED2_CSR_RELEASE_PADDED_CACHE requires ED2_CSR_COMPLEX_CACHE=1")
    reject(root,driver,"bad-flag","1","2",False,
           "ED2_CSR_RELEASE_PADDED_CACHE must be 0 or 1")
    reject(root,driver,"long-flag","1","11111111111111111111",False,
           "ED2_CSR_RELEASE_PADDED_CACHE flag too long or invalid")
    reject(root,driver,"partial-cache","1","1",True,
           "ED2_CSR_COMPLEX_CACHE requires fully cached matrix")
    (root/"summary.json").write_text(json.dumps({
        "status":"PASS", "model_cases":len(trials),
        "per_case_modes":list(MODES), "rejected_conditions":4,
        "no_presto":True,
        "fp_policy":"fast_math_numeric_baseline" if args.allow_fp_reassociation else "strict_bits",
        "memory_note":"steady-state cache accounting, not peak RSS",
        "trials":trials
    }, indent=2)+"\n", encoding="utf-8")
    print("CSR_PADDED_RELEASE_REGRESSION_PASS cases={} invalid_conditions=4".format(
        len(trials)),flush=True)


if __name__=="__main__":
    main()
