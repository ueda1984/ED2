#!/usr/bin/env python3
"""Native exact-H*x proof of ED2_CACHE_COEFF_BOUND across general spin models.

Tests real XYZ, complex DM, U(1)-conserving XXZ, and full nine-component
exchange tensors, with padded / CSR retained / CSR released layouts. All
bounded-vs-unbounded outputs within a layout must have identical IEEE bits.
The historic -Ofast case permits only the known original-vs-CSR last-bit
reassociation, not any change caused by this coefficient bound.
"""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess

from cpc_csr_main_regression import deck

LAYOUTS={"padded":("0","0"),"csr":("1","0"),"released":("1","1")}
MODELS=("real_xyz","complex_dm","u1_xxz","full_tensor")
MARKER="Hamiltonian cache: coefficient-bound width="


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):
            h.update(block)
    return h.hexdigest()


def model_inputs(work,n,model):
    if model in ("real_xyz","complex_dm"):
        return
    bonds=[]
    for j,line in enumerate((work/"bonds.dat").read_text().splitlines()):
        fields=line.split()
        assert len(fields)==11
        if model=="u1_xxz":
            v=[0.87,0.87,1.09,0,0,0,0,0,0]
        elif model=="full_tensor":
            sign=(-1)**j
            v=[0.87,0.59,1.09,0.13*sign,0.09,0.17*sign,
               0.05,0.07*sign,0.11]
        else:
            raise AssertionError(model)
        bonds.append(" ".join(fields[:2]+["{:.17g}".format(x) for x in v]))
    (work/"bonds.dat").write_text("\n".join(bonds)+"\n")
    fields=[]
    for j,line in enumerate((work/"fields.dat").read_text().splitlines()):
        index=line.split()[0]
        if model=="u1_xxz":
            v=[0.0,0.0,0.13]
        else:
            v=[0.17+0.015*((-1)**j),0.08*((-1)**j),0.13]
        fields.append(" ".join([index]+["{:.17g}".format(x) for x in v]))
    (work/"fields.dat").write_text("\n".join(fields)+"\n")


def run(driver,work,n,model,threads,layout,bounded,timeout,explicit=False):
    deck(work,n,0.12 if model=="complex_dm" else 0.0)
    model_inputs(work,n,model)
    if explicit:
        file=work/"input.dat"
        text=file.read_text()
        assert text.count("MNTE=-1")==1
        file.write_text(text.replace("MNTE=-1","MNTE=100"))
    env=dict(os.environ)
    for key in list(env):
        if key.startswith("ED2_CUBIC_") or key in (
                "ED2_CACHE_COEFF_BOUND","ED2_CSR_COMPLEX_CACHE",
                "ED2_CSR_RELEASE_PADDED_CACHE","ED2_REAL_CACHED_MATVEC"):
            env.pop(key,None)
    csr,release=LAYOUTS[layout]
    env.update(ED2_CACHE_COEFF_BOUND=str(int(bounded)),
               ED2_CSR_COMPLEX_CACHE=csr,ED2_CSR_RELEASE_PADDED_CACHE=release,
               OMP_NUM_THREADS=str(threads),OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as inp:
        p=subprocess.run([str(driver)],cwd=work,stdin=inp,env=env,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                         check=False,timeout=timeout)
    out=p.stdout.decode("utf-8",errors="replace")
    (work/"stdout.log").write_text(out)
    assert p.returncode==0,(work,p.returncode,out[-2500:])
    assert (MARKER in out)==(bounded and not explicit)
    assert ("Hamiltonian cache: explicit MNTE override preserved" in out)==(
        bounded and explicit)
    assert ("CPC_CSR_PADDED_RELEASE_PASS" in out)==(layout=="released")
    assert ("CPC_CSR_PADDED_RETAINED_PASS" in out)==(layout!="released")
    assert out.count("CPC_CSR_DENSE_REFERENCE_PASS")==(
        2 if n<=6 else 0)
    dimension=6**(n//2)
    assert "CPC_CSR_H_VECTOR_DIM={}".format(dimension) in out
    proof=work/"csr_matvec_bits.txt"
    assert len(proof.read_bytes().splitlines())==2*dimension
    found=re.search(r"Current MNTE =\s*(\d+)",out)
    assert found,(work,"missing native capacity provenance")
    width=int(found.group(1))+1
    if bounded and not explicit:
        found=re.search(
            r"coefficient-bound width=(\d+), nominal_width=(\d+), initial_bytes_saved=(\d+)",out)
        assert found,(work,"missing exact capacity accounting")
        actual,nominal,saved=map(int,found.groups())
        assert actual==width and nominal==1+2*n+8*(n-1)
        assert 1<=actual<=nominal
        assert saved==20*(nominal-actual)*dimension
    else:
        saved=0
    return {"proof":proof,"digest":sha(proof),"width":width,"saved_bytes":saved,
            "inputs":tuple(sha(work/name) for name in (
                "input.dat","spins.dat","caps.dat","bonds.dat","fields.dat"))}


def parse_values(proof):
    out=[]
    for line in proof.read_text().splitlines():
        v=line.split()
        assert len(v)==4
        r=struct.unpack(">d",int(v[2],16).to_bytes(8,"big"))[0]
        i=struct.unpack(">d",int(v[3],16).to_bytes(8,"big"))[0]
        assert math.isfinite(r) and math.isfinite(i)
        out.append(complex(r,i))
    return out


def compare_case(root,driver,n,model,threads,timeout,fast_math):
    data={}
    for layout,bound in itertools.product(LAYOUTS,(False,True)):
        work=root/"n{}-{}-t{}-{}-bounded{}".format(
            n,model,threads,layout,int(bound))
        data[layout,bound]=run(driver,work,n,model,threads,
                               layout,bound,timeout)
    assert len({r["inputs"] for r in data.values()})==1
    for layout in LAYOUTS:
        old,new=data[layout,False],data[layout,True]
        assert old["digest"]==new["digest"],(
            "Coefficient bounding changed IEEE output",n,model,threads,layout)
        assert 1<=new["width"]<=old["width"]
    assert data["csr",True]["digest"]==data["released",True]["digest"]
    if fast_math:
        reference=parse_values(data["padded",False]["proof"])
        scale=max(1.0,max(abs(z) for z in reference))
        for item in data.values():
            alternative=parse_values(item["proof"])
            assert len(reference)==len(alternative)
            maximum=max(abs(x-y) for x,y in zip(reference,alternative))
            assert maximum<=1e-12*scale,(
                "Unexpected fast-math numerical difference",n,model,threads,maximum)
    else:
        assert len({v["digest"] for v in data.values()})==1,(
            "Strict O2 exact full-vector mismatch",n,model,threads)
    saved=data["padded",True]["saved_bytes"]
    print("ED2_COEFF_BOUND_NATIVE_PASS n={} model={} threads={} saved_bytes={}"
          .format(n,model,threads,saved),flush=True)
    return {"sites":n,"model":model,"threads":threads,
            "nominal_width":data["padded",False]["width"],
            "certified_width":data["padded",True]["width"],"saved_bytes":saved}


def reject(driver,root,value,expected):
    work=root/("reject-{}".format(len(list(root.glob("reject-*")))))
    deck(work,4,0.12)
    env=dict(os.environ)
    env.update(ED2_CACHE_COEFF_BOUND=value,
               ED2_CSR_COMPLEX_CACHE="0",
               ED2_CSR_RELEASE_PADDED_CACHE="0",
               OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    with (work/"input.dat").open("rb") as inp:
        p=subprocess.run([str(driver)],cwd=work,stdin=inp,env=env,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                         timeout=60)
    assert p.returncode!=0 and expected in p.stdout.decode(errors="replace")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--sizes",default="4,6,8")
    p.add_argument("--threads",default="1,4")
    p.add_argument("--models",default=",".join(MODELS))
    p.add_argument("--timeout",type=int,default=180)
    p.add_argument("--allow-fp-reassociation",action="store_true")
    args=p.parse_args()
    sizes=[int(x) for x in args.sizes.split(",")]
    threads=[int(x) for x in args.threads.split(",")]
    models=args.models.split(",")
    assert sizes and all(n in (4,6,8,10) for n in sizes)
    assert threads and all(t>=1 for t in threads)
    assert models and all(m in MODELS for m in models)
    driver=args.driver.resolve(strict=True)
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    records=[compare_case(root,driver,n,mode,t,args.timeout,args.allow_fp_reassociation)
             for n,mode,t in itertools.product(sizes,models,threads)]
    extra=run(driver,root/"explicit-mnte",4,"complex_dm",1,
              "padded",True,args.timeout,explicit=True)
    assert extra["width"]==101 and extra["saved_bytes"]==0
    reject(driver,root,"2","ED2_CACHE_COEFF_BOUND must be 0 or 1")
    reject(driver,root,"1111111111111111111111",
           "ED2_CACHE_COEFF_BOUND flag too long or invalid")
    summary={"status":"PASS","cases":len(records),"records":records,
             "explicit_MNTE_preserved":True,"invalid_flags":2,
             "policy":"fast_math_numeric" if args.allow_fp_reassociation else "strict_bits"}
    (root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("ED2_COEFF_BOUND_REGRESSION_PASS cases={} rejections=2"
          .format(len(records)),flush=True)


if __name__=="__main__":
    main()
