#!/usr/bin/env python3
"""Prove direct two-pass complex CSR construction matches the native padded H.

Physical alternating mixed S=1/2,1 XYZ/DM/XXZ/full-tensor Hamiltonians,
N=4/6/8/10, one- and multi-threaded, two nonreal complex H*x vectors.
The single Fortran row generator is used with a full padded layout or
64-row staging buffers. No PRESTO access or scientific Python required.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess

from cpc_csr_main_regression import deck
from cpc_coeff_cache_bound_regression import model_inputs, parse_values

MODES={
    "padded":("0","0","0"),
    "release":("1","1","0"),
    "direct":("1","1","1"),
}
MARKER="Hamiltonian action: csr-direct-build passes="


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(driver,work,n,model,threads,mode,bounded,timeout):
    deck(work,n,0.12 if model=="complex_dm" else 0.)
    model_inputs(work,n,model)
    csr,release,direct=MODES[mode]
    env=dict(os.environ)
    for key in list(env):
        if key.startswith("ED2_CUBIC_") or key in (
                "ED2_CSR_COMPLEX_CACHE","ED2_CSR_RELEASE_PADDED_CACHE",
                "ED2_CSR_DIRECT_BUILD","ED2_CACHE_COEFF_BOUND",
                "ED2_REAL_CACHED_MATVEC"):
            env.pop(key,None)
    env.update(ED2_CSR_COMPLEX_CACHE=csr,
               ED2_CSR_RELEASE_PADDED_CACHE=release,
               ED2_CSR_DIRECT_BUILD=direct,
               ED2_CACHE_COEFF_BOUND=str(int(bounded)),
               OMP_NUM_THREADS=str(threads),OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as f:
        p=subprocess.run([str(driver)],cwd=work,stdin=f,
                         env=env,stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,timeout=timeout)
    log=p.stdout.decode(errors="replace")
    (work/"stdout.log").write_text(log)
    assert p.returncode==0,(work,p.returncode,log[-2500:])
    assert (MARKER in log)==(mode=="direct")
    assert ("CPC_CSR_PADDED_RELEASE_PASS" in log)==(mode!="padded")
    assert ("CPC_CSR_PADDED_RETAINED_PASS" in log)==(mode=="padded")
    assert log.count("CPC_CSR_DENSE_REFERENCE_PASS")==(
        2 if n<=6 else 0)
    dimension=6**(n//2)
    assert "CPC_CSR_H_VECTOR_DIM={}".format(dimension) in log
    proof=work/"csr_matvec_bits.txt"
    assert len(proof.read_bytes().splitlines())==2*dimension
    if mode=="direct":
        stage=re.search(
            r"csr-direct-build passes=(\d+), max_staging_bytes=(\d+), peak_cache_payload_upper_bytes=(\d+)",
            log)
        csr_m=re.search(
            r"csr-complex cache entries=(\d+), extra_bytes=(\d+)",log)
        assert stage and csr_m,(work,"missing native CSR memory provenance")
        passes,stage_bytes,peak=map(int,stage.groups())
        entries,csr_bytes=map(int,csr_m.groups())
        assert passes==2 and entries>=dimension
        assert peak==csr_bytes+stage_bytes and stage_bytes>0
        width=int(re.search(r"Current MNTE =\s*(\d+)",log).group(1))+1
        assert stage_bytes==20*width*min(64,dimension)
        memory={"peak_cache_payload_upper_bytes":peak,
                "max_staging_bytes":stage_bytes,
                "csr_bytes":csr_bytes}
    else:
        memory={}
    return {"hash":sha(proof),"proof":proof,
            "input":tuple(sha(work/name) for name in (
                "input.dat","spins.dat","caps.dat","bonds.dat","fields.dat")),
            **memory}


def compare(root,driver,n,model,t,timeout,fast_math):
    data={}
    for mode,bound in itertools.product(MODES,(False,True)):
        work=root/"s{}-{}-t{}-{}-bound{}".format(
            n,model,t,mode,int(bound))
        data[mode,bound]=run(driver,work,n,model,t,mode,bound,timeout)
    assert len({v["input"] for v in data.values()})==1
    for mode in MODES:
        assert data[mode,False]["hash"]==data[mode,True]["hash"],(
            "Coefficient width changed full H*x IEEE bits",n,model,t,mode)
    if fast_math:
        base=parse_values(data["padded",False]["proof"])
        scale=max(1.0,max(abs(z) for z in base))
        for item in data.values():
            alternative=parse_values(item["proof"])
            assert len(base)==len(alternative)
            assert max(abs(a-b) for a,b in zip(base,alternative))<=1e-12*scale
    else:
        assert len({v["hash"] for v in data.values()})==1,(
            "Strict full native complex H*x bit difference",n,model,t)
    chosen=data["direct",True]
    print("CSR_DIRECT_PHYSICS_PASS sites={} model={} threads={}"
          .format(n,model,t),flush=True)
    return {"sites":n,"model":model,"threads":t,
            "direct_peak_cache_payload_upper":chosen["peak_cache_payload_upper_bytes"],
            "csr_bytes":chosen["csr_bytes"],
            "maximum_staging_bytes":chosen["max_staging_bytes"]}


def reject(root,driver,name,overrides,expected,partial=False):
    work=root/("reject-"+name)
    deck(work,4,0.12,partial=partial)
    env=dict(os.environ)
    env.update(ED2_CSR_COMPLEX_CACHE="1",
               ED2_CSR_RELEASE_PADDED_CACHE="1",
               ED2_CSR_DIRECT_BUILD="1",
               ED2_CACHE_COEFF_BOUND="1",
               OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    env.update(overrides)
    with (work/"input.dat").open("rb") as f:
        p=subprocess.run([str(driver)],cwd=work,stdin=f,
                         env=env,stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT,timeout=60)
    output=p.stdout.decode(errors="replace")
    assert p.returncode!=0 and expected in output,(
        name,p.returncode,output[-1500:])
    print("CSR_DIRECT_INVALID_INPUT_REJECTED "+name,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="4,6,8")
    p.add_argument("--threads",default="1,4")
    p.add_argument("--models",
                   default="real_xyz,complex_dm,u1_xxz,full_tensor")
    p.add_argument("--timeout",type=int,default=180)
    p.add_argument("--allow-fp-reassociation",action="store_true")
    args=p.parse_args()
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    driver=args.driver.resolve(strict=True)
    sizes=[int(n) for n in args.sizes.split(",")]
    threads=[int(t) for t in args.threads.split(",")]
    models=args.models.split(",")
    assert sizes and all(n in (4,6,8,10) for n in sizes)
    assert threads and all(t>0 for t in threads)
    valid=("real_xyz","complex_dm","u1_xxz","full_tensor")
    assert models and all(m in valid for m in models)
    rows=[compare(root,driver,n,model,t,args.timeout,
                  args.allow_fp_reassociation)
          for n,model,t in itertools.product(sizes,models,threads)]
    reject(root,driver,"missing-csr",
           {"ED2_CSR_COMPLEX_CACHE":"0"},
           "ED2_CSR_DIRECT_BUILD requires ED2_CSR_COMPLEX_CACHE=1")
    reject(root,driver,"missing-release",
           {"ED2_CSR_RELEASE_PADDED_CACHE":"0"},
           "ED2_CSR_DIRECT_BUILD requires ED2_CSR_RELEASE_PADDED_CACHE=1")
    reject(root,driver,"invalid-flag",
           {"ED2_CSR_DIRECT_BUILD":"2"},
           "ED2_CSR_DIRECT_BUILD must be 0 or 1")
    reject(root,driver,"oversized-flag",
           {"ED2_CSR_DIRECT_BUILD":"11111111111111111111"},
           "ED2_CSR_DIRECT_BUILD flag too long or invalid")
    reject(root,driver,"partial-cache",{},
           "ED2_CSR_DIRECT_BUILD requires fully cached matrix",
           partial=True)
    (root/"summary.json").write_text(json.dumps({
        "status":"PASS","cases":len(rows),"results":rows,
        "fast_math":args.allow_fp_reassociation,"invalid":5,
        "memory_note":"cache-array payload upper bound, NOT physical peak RSS"
    },indent=2)+"\n")
    print("CSR_DIRECT_REGRESSION_PASS cases={} invalid=5".format(len(rows)),
          flush=True)


if __name__=="__main__":
    main()
