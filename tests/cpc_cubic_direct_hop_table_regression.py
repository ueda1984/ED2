#!/usr/bin/env python3
"""Exact full-cache, Ritz and group-character proof for dense direct-hop table.

Build identical production local-XXZ Hamiltonian caches in three-magnon
Gamma/O_h with and without ED2_CUBIC_DIRECT_HOP_TABLE=1.  Verify *all*
diagonal/offdiagonal IEEE-754 bytes and cache record targets agree exactly,
including every 1D character.  For small volumes also run certified Lanczos
with both paths and compare energies/direct residuals.  Negative tests must
fail closed if all exact-cubic prerequisite switches are not enabled.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

from cpc_cubic_pointgroup import (canonical_bonds, sector_dimensions,
                                  symmetry_generators, write_input)
from cpc_pointgroup_profile import run_case


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as src:
        for block in iter(lambda: src.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def make_env(threads, direct):
    env=dict(os.environ)
    env.update(ED2_CUBIC_ANCHORED_REP="1",
               ED2_CUBIC_FAST_GROUP_VALIDATION="1",
               ED2_CUBIC_ANCHORED_BASIS="1",
               ED2_CUBIC_LOCAL_BOND_CACHE="1",
               ED2_CUBIC_CACHE_PROOF_DUMP="1",
               ED2_RNG_SEED="1984", OMP_NUM_THREADS=str(threads),
               OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1")
    for key in ("ED2_CUBIC_DIRECT_HOP_TABLE", "ED2_CUBIC_CACHE_PROFILE",
                "ED2_FORCE_LEGACY_ORBITS", "ED2_CUBIC_ANCHORED_VERIFY"):
        env.pop(key, None)
    if direct:
        env["ED2_CUBIC_DIRECT_HOP_TABLE"]="1"
    return env


def cache_case(root, driver, side, sector, threads, direct, timeout):
    dest=root/"n{}-{}-t{}-{}".format(side,sector,threads,
                                     "direct" if direct else "reference")
    write_input(dest,side,3,symmetry_generators(side,sector),
                canonical_bonds(side),mnte=64)
    env=make_env(threads,direct)
    with (dest/"input.dat").open("rb") as inp:
        run=subprocess.run([str(driver)],stdin=inp,stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,cwd=str(dest),
                           timeout=timeout,env=env,check=False)
    log=run.stdout.decode("utf-8",errors="replace")
    (dest/"driver.log").write_text(log)
    assert run.returncode==0,(dest,run.returncode,log[-3000:])
    assert ("Representative lookup: cubic-direct-hop-pairs" in log) is direct
    assert "Representative lookup: cubic-anchored-orbits" in log
    assert "Cubic cache builder: local-bond-sparse" in log
    assert "Basis builder: anchored-cubic" in log
    dim=sector_dimensions(side)[sector]
    assert "Cubic cache proof dump: projected rows={}".format(dim) in log
    return dest


def ritz_case(root, exe, side, sector, direct, timeout):
    dest=root/"physical-n{}-{}-{}".format(
        side,sector,"direct" if direct else "reference")
    write_input(dest,side,3,symmetry_generators(side,sector),
                canonical_bonds(side),mnte=64)
    env=make_env(1,direct)
    for key in ("ED2_CUBIC_CACHE_PROOF_DUMP",):
        env.pop(key,None)
    knobs=tuple(x for x in env if x.startswith("ED2_CUBIC_"))
    prev={key:os.environ.get(key) for key in knobs}
    try:
        for key in knobs:
            os.environ[key]=env[key]
        for key in ("ED2_CUBIC_DIRECT_HOP_TABLE","ED2_CUBIC_CACHE_PROFILE",
                    "ED2_CUBIC_CACHE_PROOF_DUMP","ED2_FORCE_LEGACY_ORBITS",
                    "ED2_CUBIC_ANCHORED_VERIFY"):
            if key not in env:
                os.environ.pop(key,None)
        obs=run_case(exe,dest,side**3,3,sector,sector_dimensions(side)[sector],
                     1,"two-stage",timeout)
    finally:
        for key,val in prev.items():
            if val is None:
                os.environ.pop(key,None)
            else:
                os.environ[key]=val
        for key in ("ED2_CUBIC_DIRECT_HOP_TABLE",):
            if key not in prev:
                os.environ.pop(key,None)
    assert obs["residual"]<=1.01e-8
    assert obs["dimension"]==sector_dimensions(side)[sector]
    return obs


def reject_bad_flag(driver, root):
    side=4
    sector="oh-a1g"
    dest=root/"reject-direct-without-local-bond-cache"
    write_input(dest,side,3,symmetry_generators(side,sector),
                canonical_bonds(side),mnte=64)
    env=make_env(1,True)
    env.pop("ED2_CUBIC_LOCAL_BOND_CACHE",None)
    with (dest/"input.dat").open("rb") as inp:
        res=subprocess.run([str(driver)],stdin=inp,stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,cwd=str(dest),env=env,
                           timeout=60,check=False)
    assert res.returncode!=0
    assert b"Cubic direct hop missing required flag: ED2_CUBIC_LOCAL_BOND_CACHE" in res.stdout
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--driver",type=Path,required=True)
    p.add_argument("--exe",type=Path,required=True)
    p.add_argument("--sides",default="4,6,10")
    p.add_argument("--threads",default="1,4")
    p.add_argument("--sectors",default="oh-a1g,oh-a2g,oh-a1u,oh-a2u")
    p.add_argument("--timeout",type=int,default=1200)
    p.add_argument("--skip-lanczos",action="store_true")
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    driver=args.driver.resolve(strict=True)
    exe=args.exe.resolve(strict=True)
    root=args.output or Path(tempfile.mkdtemp(prefix="ed2-direct-hop-table-"))
    if args.output:
        root.mkdir(parents=True,exist_ok=False)
    sides=tuple(int(x) for x in args.sides.split(","))
    threads=tuple(int(x) for x in args.threads.split(","))
    sectors=tuple(args.sectors.split(","))
    assert set(sectors)<=set(("oh-a1g","oh-a2g","oh-a1u","oh-a2u"))
    report=[]
    for side,sector,t in itertools.product(sides,sectors,threads):
        old=cache_case(root,driver,side,sector,t,False,args.timeout)
        new=cache_case(root,driver,side,sector,t,True,args.timeout)
        assert sha(old/"input.dat")==sha(new/"input.dat"),(side,sector,t)
        old_hash=sha(old/"cubic_cache_proof_dump.txt")
        new_hash=sha(new/"cubic_cache_proof_dump.txt")
        assert old_hash==new_hash,("Direct hop altered cache bits or order",side,sector,t)
        report.append(dict(side=side,sector=sector,threads=t,
                           dimension=sector_dimensions(side)[sector],
                           input_sha256=sha(new/"input.dat"),
                           full_cache_sha256=new_hash))
        print("CUBIC_DIRECT_HOP_BITS_PASS side={} sector={} threads={} dim={}".format(
            side,sector,t,sector_dimensions(side)[sector]),flush=True)
        if side in (4,6) and t==1 and not args.skip_lanczos:
            a=ritz_case(root,exe,side,sector,False,args.timeout)
            b=ritz_case(root,exe,side,sector,True,args.timeout)
            assert a["input_sha256"]==b["input_sha256"]
            assert abs(a["energy"]-b["energy"])<1e-8
            assert a["residual"]<=1.01e-8 and b["residual"]<=1.01e-8
            print("CUBIC_DIRECT_HOP_RITZ_PASS side={} sector={} deltaE={:+.3e}".format(
                side,sector,b["energy"]-a["energy"]),flush=True)
    reject_bad_flag(driver,root)
    (root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print("CUBIC_DIRECT_HOP_REGRESSION_PASS {} exact full-cache pairs "
          "and negative model check".format(len(report)),flush=True)


if __name__=="__main__":
    main()
