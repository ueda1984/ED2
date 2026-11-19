#!/usr/bin/env python3
"""Exercise uncached transitions after avoiding unnecessary full-cache scratch.

The original ED2 on-the-fly rows use temporary state arrays. This native
test checks that fully cached and wk_dim=1 partial caching both compute
the same complex H*x to strict numerical tolerance on physical mixed-spin
non-U1 real and nonreal DM models. Each mode independently matches the
separate full dense-H builder at small N. No HPC or scientific wheels.
"""
import argparse
import itertools
import json
from pathlib import Path
import math

from cpc_csr_main_regression import deck, matching_inputs, COLUMNS
from cpc_csr_release_regression import run, parse_complex_dump


def one_case(root, driver, n, dm, threads, timeout):
    records=[]
    for partial in (False,True):
        work=root/"n{}-dm{}-t{}-partial{}".format(n,int(dm!=0),threads,int(partial))
        deck(work,n,dm,partial=partial)
        rc,log=run(driver,work,threads,"0","0",timeout)
        assert rc==0,(work,log[-2000:])
        assert "CPC_CSR_PADDED_RETAINED_PASS" in log
        assert "CPC_CSR_PADDED_RELEASE_PASS" not in log
        assert log.count("CPC_CSR_DENSE_REFERENCE_PASS")==(
            2 if n<=6 else 0), (work, "independent dense Hamiltonian failed")
        dump=work/"csr_matvec_bits.txt"
        rows=6**(n//2)
        assert len(dump.read_bytes().splitlines())==2*rows
        payload=parse_complex_dump(dump.read_bytes())
        assert len(payload)==2*rows and all(
            math.isfinite(z.real) and math.isfinite(z.imag) for z in payload)
        records.append({
            "partial":partial,
            "files":{name:work.joinpath(name).read_bytes() for name in COLUMNS},
            "vector":payload,
        })
    for name in COLUMNS:
        if name=="input.dat":
            old=records[0]["files"][name]
            new=records[1]["files"][name]
            assert old.count(b"wk_dim=-1")==1 and new.count(b"wk_dim=1")==1
            assert old.replace(b"wk_dim=-1",b"wk_dim=1")==new
        else:
            assert records[0]["files"][name]==records[1]["files"][name]
    full,partial=records[0]["vector"],records[1]["vector"]
    max_delta=max(abs(x-y) for x,y in zip(full,partial))
    scale=max(1.0,max(abs(x) for x in full))
    assert max_delta <= 1e-10*scale, (
        "Full vs partially cached H*x disagreement",n,dm,threads,max_delta,scale)
    print("MAIN_PARTIAL_CACHE_NATIVE_HV_PASS n={} dm={} t={} max_delta={:.3e}".format(
        n,dm,threads,max_delta),flush=True)
    return {
        "n":n,"dm":dm,"threads":threads,
        "full_vs_partial_max_delta":max_delta,
        "scale":scale,"dense_reference_per_mode":n<=6,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="4,6,8")
    p.add_argument("--threads",default="1,4")
    p.add_argument("--timeout",type=int,default=180)
    args=p.parse_args()
    sizes=[int(x) for x in args.sizes.split(",")]
    threads=[int(x) for x in args.threads.split(",")]
    assert sizes and all(n in (4,6,8,10) for n in sizes)
    assert threads and all(t>0 for t in threads)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    driver=args.driver.resolve(strict=True)
    records=[
        one_case(root,driver,n,dm,t,args.timeout)
        for n,dm,t in itertools.product(sizes,(0.0,0.12),threads)
    ]
    (root/"summary.json").write_text(
        json.dumps({"status":"PASS","cases":len(records),"records":records},indent=2)+"\n")
    print("MAIN_PARTIAL_CACHE_REGRESSION_PASS cases={}".format(len(records)),flush=True)


if __name__=="__main__":
    main()
