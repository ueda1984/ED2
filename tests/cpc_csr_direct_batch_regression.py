#!/usr/bin/env python3
"""Bitwise native H*x regression for optional direct two-pass CSR batch widths.

All full native complex output vector bits must remain identical when only
the CSR 2-pass staging row size changes. Real, DM, XXZ, and fully populated
exchange tensors, independent dense H reference for N<=6, and negative
flag validation. No scientific Python packages or PRESTO.
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

from cpc_csr_main_regression import deck
from cpc_coeff_cache_bound_regression import model_inputs

MARKER = "Hamiltonian action: csr-direct-batch-rows="
MODELS = ("real_xyz", "complex_dm", "u1_xxz", "full_tensor")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def start(driver,work,n,model,t,rows,timeout):
    deck(work,n,0.12 if model=="complex_dm" else 0.0)
    model_inputs(work,n,model)
    env=dict(os.environ)
    for key in tuple(env):
        if key.startswith("ED2_CUBIC_") or key in (
            "ED2_CSR_COMPLEX_CACHE","ED2_CSR_RELEASE_PADDED_CACHE",
            "ED2_CSR_DIRECT_BUILD","ED2_CACHE_COEFF_BOUND",
            "ED2_CSR_DIRECT_BATCH_ROWS"):
            env.pop(key,None)
    env.update(ED2_CSR_COMPLEX_CACHE="1",
               ED2_CSR_RELEASE_PADDED_CACHE="1",
               ED2_CSR_DIRECT_BUILD="1",
               ED2_CACHE_COEFF_BOUND="1",
               OMP_NUM_THREADS=str(t),OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    if rows is not None:
        env["ED2_CSR_DIRECT_BATCH_ROWS"]=str(rows)
    with (work/"input.dat").open("rb") as deck_file:
        proc=subprocess.run([str(driver)],cwd=work,stdin=deck_file,
                            env=env,stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,timeout=timeout)
    text=proc.stdout.decode(errors="replace")
    (work/"stdout.log").write_text(text)
    assert proc.returncode==0,(work,proc.returncode,text[-2200:])
    chosen=64 if rows is None else rows
    assert MARKER+str(chosen) in text,(work,"wrong selected batch provenance")
    assert "CPC_CSR_PADDED_RELEASE_PASS" in text
    assert text.count("CPC_CSR_DENSE_REFERENCE_PASS")==(
        2 if n<=6 else 0)
    dimension=6**(n//2)
    assert "CPC_CSR_H_VECTOR_DIM={}".format(dimension) in text
    bits=work/"csr_matvec_bits.txt"
    assert len(bits.read_bytes().splitlines())==2*dimension
    info=re.search(
        r"csr-direct-build passes=(\d+), max_staging_bytes=(\d+), peak_cache_payload_upper_bytes=(\d+)",text)
    values=re.search(r"csr-complex cache entries=(\d+), extra_bytes=(\d+)",text)
    width=re.search(r"Current MNTE =\s*(\d+)",text)
    assert info and values and width
    passes,stage_bytes,peak=map(int,info.groups())
    entries,csr_bytes=map(int,values.groups())
    row_width=int(width.group(1))+1
    assert passes==2 and entries>=dimension
    assert peak==csr_bytes+stage_bytes
    assert stage_bytes==20*row_width*min(chosen,dimension),(
        "Unaccounted source batch peak",chosen,stage_bytes,row_width,dimension)
    return {
        "sha256":sha(bits),
        "input_sha256":tuple(sha(work/name) for name in (
            "input.dat","spins.dat","caps.dat","bonds.dat","fields.dat")),
        "peak_cache_bytes_upper":peak,
        "stage_bytes":stage_bytes,
        "csr_bytes":csr_bytes,
    }


def reject(driver,root,name,direct,batch,expect,timeout):
    work=root/("reject-"+name)
    deck(work,4,0.12)
    env=dict(os.environ)
    for key in tuple(env):
        if key.startswith("ED2_CUBIC_") or key in (
            "ED2_CSR_COMPLEX_CACHE","ED2_CSR_RELEASE_PADDED_CACHE",
            "ED2_CSR_DIRECT_BUILD","ED2_CACHE_COEFF_BOUND",
            "ED2_CSR_DIRECT_BATCH_ROWS"):
            env.pop(key,None)
    env.update(ED2_CSR_COMPLEX_CACHE="1",
               ED2_CSR_RELEASE_PADDED_CACHE="1",
               ED2_CSR_DIRECT_BUILD=str(direct),
               ED2_CSR_DIRECT_BATCH_ROWS=batch,
               OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as inp:
        proc=subprocess.run([str(driver)],cwd=work,stdin=inp,
                            env=env,stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,timeout=timeout)
    log=proc.stdout.decode(errors="replace")
    assert proc.returncode!=0 and expect in log,(name,proc.returncode,log[-2000:])
    print("CSR_BATCH_INVALID_FLAG_REJECTED "+name,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="4,6,10")
    p.add_argument("--threads",default="1,4")
    p.add_argument("--models",default="real_xyz,complex_dm,u1_xxz,full_tensor")
    p.add_argument("--batches",default="1,64,256,512,1024")
    p.add_argument("--timeout",type=int,default=180)
    args=p.parse_args()
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    driver=args.driver.resolve(strict=True)
    ns=[int(x) for x in args.sizes.split(",")]
    ts=[int(x) for x in args.threads.split(",")]
    models=args.models.split(",")
    batches=[int(x) for x in args.batches.split(",")]
    assert ns and all(n in (4,6,8,10) for n in ns)
    assert ts and all(t>=1 for t in ts)
    assert models and all(x in MODELS for x in models)
    assert batches and all(1<=x<=1024 for x in batches)
    cases=[]
    for n,model,t in itertools.product(ns,models,ts):
        base=start(driver,root/("n{}-{}-t{}-default".format(n,model,t)),
                   n,model,t,None,args.timeout)
        for batch in batches:
            other=start(driver,root/("n{}-{}-t{}-b{}".format(n,model,t,batch)),
                        n,model,t,batch,args.timeout)
            assert base["sha256"]==other["sha256"],(
                "Batch size changed full IEEE complex output",n,model,t,batch)
            assert base["input_sha256"]==other["input_sha256"]
            assert base["csr_bytes"]==other["csr_bytes"]
            cases.append({"n":n,"model":model,"threads":t,"batch":batch,
                          "stage_bytes":other["stage_bytes"],
                          "peak_cache_bytes_upper":other["peak_cache_bytes_upper"]})
        print("CSR_BATCH_NATIVE_BITWISE_PASS n={} model={} t={} batches={}".format(
            n,model,t,len(batches)),flush=True)

    for name,direct,batch,expect in [
        ("zero","1","0","out of range [1,1024]"),
        ("too-large","1","1025","out of range [1,1024]"),
        ("negative","1","-2","must be an integer in [1,1024]"),
        ("alpha","1","five","must be an integer in [1,1024]"),
        ("blank","1","  ","must be an integer in [1,1024]"),
        ("overflow","1","9999999999999999999999","integer parsing failed"),
        ("oversized","1","1"*40,"flag too long or invalid"),
        ("missing-direct","0","256","requires ED2_CSR_DIRECT_BUILD=1")]:
        reject(driver,root,name,direct,batch,
               "ED2_CSR_DIRECT_BATCH_ROWS "+expect,args.timeout)
    (root/"summary.json").write_text(json.dumps({
        "status":"PASS","settings":len(cases),"physical_cases":len(cases)//len(batches),
        "invalid":8,"mode":"strict same-CSR IEEE bits under each compiler",
        "cases":cases},indent=2)+"\n")
    print("CSR_BATCH_REGRESSION_PASS settings={} invalid=8".format(len(cases)),
          flush=True)


if __name__=="__main__":
    main()
