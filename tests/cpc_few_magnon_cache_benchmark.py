#!/usr/bin/env python3
"""Controlled-design ED2 few-magnon cache experiment, pilot-capable.

Compares no/partial/full cache at fixed physics, fixed initial vector, fixed
RNG seed, fixed residual contract and identical source/binary. Explicitly
labels shared GitHub CI as exploratory; the source is not instrumented.
The Hamiltonian is the exactly validated Gamma-sector periodic cubic XXZ
model. Do not use MNTE=6*n_down for a general XYZ/DM/mixed-spin model.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import subprocess
import time

from cpc_qs3_ed_large_few_magnon import create, sector_dim

CASES={
    "216-n3": ((6,6,6),3,(0,0,0),148.56489075771222),
    "1000-n2": ((10,10,10),2,(0,0,0),738.0067742722034),
}
POLICIES=("none","half_tight","full_tight","full_generic")
SEED=1984

def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()

def initial_vector(dim):
    # No NumPy needed, explicitly identical in every test and repeat.
    return "".join(f"{math.sin((i+1)*0.113):.17g} {math.cos((i+1)*0.197):.17g}\n"
                   for i in range(dim))

def settings(policy,dim,nsites,n,nbonds):
    tight=6*n # One magnon has at most six directed nearest-neighbor hops.
    broad=2*nsites+8*nbonds
    if policy=="none":return (0,tight)
    if policy=="half_tight":return (dim//2,tight)
    if policy=="full_tight":return (-1,tight)
    if policy=="full_generic":return (-1,-1)
    raise ValueError("Unknown policy")

def diag_parse(path):
    if not path.is_file():raise AssertionError("Missing ED2 Krylov diagnostics")
    vals={}
    for line in path.read_text().splitlines():
        if "=" in line:
            k,v=line.split("=",1)
            vals[k.strip()]=v.strip()
    if vals.get("status") is None or vals.get("h_applications") is None:
        raise AssertionError("Diagnostics are incomplete")
    return vals

def prepare(work,lengths,n,k,policy,dim,initial):
    bonds,xy,zz=create(work,lengths,n,k,False)
    assert xy==zz==1.
    wk,mnte=settings(policy,dim,math.prod(lengths),n,len(bonds))
    deck=work/"input.dat"
    old=deck.read_text()
    needle="wk_dim=0,MNTE=-1"
    if old.count(needle)!=1:raise AssertionError("Unique cache settings not found")
    deck.write_text(old.replace(needle,f"wk_dim={wk},MNTE={mnte}"))
    (work/"initial.txt").write_text(initial)
    physics_files=["bonds.dat","spins.dat","caps.dat","fields.dat",
                   "shift1.dat","shift2.dat","shift3.dat"]
    hash_physics=hashlib.sha256()
    for file in physics_files:
        data=(work/file).read_bytes()
        hash_physics.update(file.encode()+b"\0"+data)
    return dict(wk_dim=wk,MNTE=mnte,NBonds=len(bonds),
                physical_input_sha256=hash_physics.hexdigest(),
                input_sha256=sha_file(deck),initial_sha256=sha_file(work/"initial.txt"))

def execute(exe,work,meta,dim,target,tol,timeout,cpus,phase,rep,policy,case):
    env=dict(os.environ,OMP_NUM_THREADS="1",OMP_DYNAMIC="FALSE",
             OMP_PROC_BIND="TRUE",OMP_PLACES="cores",
             OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
             GOTO_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1",
             ED2_RNG_SEED=str(SEED),
             ED2_INITIAL_VECTOR=str((work/"initial.txt").resolve()),
             ED2_DIAGNOSTICS="diagnostics.txt",ED2_RESIDUAL_ATOL=str(tol))
    timing=work/"gnu-time.txt"
    cmd=["/usr/bin/time","-f","RSS_KIB=%M\nELAPSED_S=%e\nCPU_USER_S=%U\nCPU_SYS_S=%S",
         "-o",str(timing),str(exe)]
    t0=time.perf_counter()
    with (work/"input.dat").open("rb") as input_stream:
        try:
            run=subprocess.run(cmd,cwd=work,stdin=input_stream,env=env,
                               stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                               timeout=timeout,check=False)
        except subprocess.TimeoutExpired as e:
            (work/"stdout.log").write_bytes(e.stdout or b"")
            raise TimeoutError(f"Timeout {timeout}s: {case}/{policy}/{rep}")
    runtime=time.perf_counter()-t0
    (work/"stdout.log").write_bytes(run.stdout)
    out=run.stdout.decode(errors="replace")
    data=dict(case=case,phase=phase,repeat=rep,policy=policy,dimension=dim,
              source="ED2 native ALG=1",start_seed=SEED,
              selected_cpu_ids=cpus,wall_seconds=runtime,
              process_returncode=run.returncode,**meta)
    for line in timing.read_text().splitlines():
        if "=" in line:
            key,val=line.split("=",1)
            if key=="RSS_KIB":data["peak_rss_kib"]=int(val)
            if key=="ELAPSED_S":data["gnu_wall_seconds"]=float(val)
            if key=="CPU_USER_S":data["user_cpu_seconds"]=float(val)
            if key=="CPU_SYS_S":data["system_cpu_seconds"]=float(val)
    if run.returncode:
        raise RuntimeError(f"ED2 failed ({run.returncode}): {case}/{policy} {out[-1300:]}")
    pattern=lambda p:re.findall(p,out,re.I)
    d=pattern(r"THS\(k\)\s*=\s*(\d+)")
    if len(d)!=1 or int(d[0])!=dim:
        raise AssertionError(f"Unexpected momentum dimension: {d}; expected {dim}")
    e=pattern(r"E_0\s*=\s*([-+]?\d+(?:\.\d+)?[EDed][-+]?\d+)")
    if len(e)!=1:raise AssertionError("ED2 did not print exactly one E_0")
    energy=float(e[0].replace("D","E").replace("d","e"))
    if abs(energy-target)>2e-7:
        raise AssertionError(f"Independent reference mismatch {energy} vs {target}")
    current=pattern(r"Current MNTE\s*=\s*(\d+)")
    optimal=pattern(r"Optimal MNTE\s*=\s*(\d+)")
    if len(current)!=1 or len(optimal)!=1:
        raise AssertionError("Cache row capacity diagnostics not recorded")
    expected=(2*math.prod(CASES[case][0])+8*meta["NBonds"]
              if meta["MNTE"]<0 else meta["MNTE"])
    if int(current[0])!=expected or int(optimal[0])>expected:
        raise AssertionError(f"Cache row bounds violated: {current} {optimal} {expected}")
    diag=diag_parse(work/"diagnostics.txt")
    if diag["status"]!="0" or diag.get("explicit_initial")!="T" or diag.get("seed_specified")!="T":
        raise AssertionError(f"Unaccepted numerical run: {diag}")
    if diag.get("requested_residual_met_1")!="T":
        raise AssertionError("User-specified residual policy not met")
    pairs=[line.split() for line in (work/"diagnostics.txt").read_text().splitlines()
           if len(line.split())==8 and line.split()[0]=="1"]
    if len(pairs)!=1:raise AssertionError("Expected one accepted Ritz pair")
    residual=float(pairs[0][5])
    if residual>tol*(1+3e-13):
        raise AssertionError(f"Direct residual {residual} exceeds {tol}")
    actual_dim=dim if meta["wk_dim"]<0 else meta["wk_dim"]
    allocated_entries=1+expected
    nominal_bytes=20*allocated_entries*actual_dim
    data.update(energy=energy,reference_E0=target,
                independent_energy_error=abs(energy-target),
                direct_residual=residual,
                h_applications=int(diag["h_applications"]),
                iterations=int(diag["iterations"]),
                reported_optimal_mnte=int(optimal[0]),
                reported_current_mnte=int(current[0]),
                cached_rows=actual_dim,
                nominal_cache_bytes=nominal_bytes,
                diagnostics_sha256=sha_file(work/"diagnostics.txt"),
                stdout_sha256=sha_file(work/"stdout.log"))
    return data

def stats(vals):
    if not vals:return None
    med=statistics.median(vals)
    return dict(median=med,mad=statistics.median(abs(x-med) for x in vals),
                minimum=min(vals),maximum=max(vals),n=len(vals))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--cases",default="1000-n2,216-n3")
    p.add_argument("--policies",default="none,half_tight,full_tight")
    p.add_argument("--warmups",type=int,default=1)
    p.add_argument("--repeats",type=int,default=2)
    p.add_argument("--timeout",type=int,default=220)
    p.add_argument("--residual-atol",type=float,default=1e-8)
    p.add_argument("--cpu",type=int,default=None,
                   help="One CPU in the existing process's allowed affinity")
    p.add_argument("--allow-shared",action="store_true",
                   help="Explicitly label result as non-publication exploratory")
    p.add_argument("--exclusive-allocation-verified",action="store_true",
                   help="Require operator's independent scheduler attestation")
    a=p.parse_args()
    cases=a.cases.split(",")
    policies=a.policies.split(",")
    if len(set(cases))!=len(cases) or any(x not in CASES for x in cases):p.error("Invalid case set")
    if len(set(policies))!=len(policies) or any(x not in POLICIES for x in policies):
        p.error("Invalid cache policies")
    if not a.allow_shared and not a.exclusive_allocation_verified:
        p.error("Specify exactly one of --allow-shared or --exclusive-allocation-verified")
    if a.allow_shared and a.exclusive_allocation_verified:
        p.error("Contradictory environment attestation flags")
    if not 0<=a.warmups<=4 or not 1<=a.repeats<=20 or not 30<=a.timeout<=1800:
        p.error("Invalid timing controls")
    if not 0<a.residual_atol<=1e-6:p.error("Invalid residual tolerance")
    # Do not accidentally allocate nearly a GiB in 216-site full-generic.
    if "full_generic" in policies and "216-n3" in cases:
        p.error("216-site full generic cache is excluded from this benchmark (estimated 824 MiB)")
    if not hasattr(os,"sched_getaffinity"):
        p.error("Linux CPU affinity support required")
    allowed=sorted(os.sched_getaffinity(0))
    cpu=a.cpu if a.cpu is not None else allowed[0]
    if cpu not in allowed:p.error("CPU not in current allowed affinity")
    os.sched_setaffinity(0,{cpu})
    exe=a.exe.resolve(strict=True)
    root=a.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    allrows=[];violations=[];physical_hashes={}
    for name in cases:
        lengths,n,k,reference=CASES[name]
        dim=sector_dim(lengths,n,k)
        initial=initial_vector(dim)
        for rep in range(a.warmups+a.repeats):
            phase="warmup" if rep<a.warmups else "measurement"
            ordered=policies.copy()
            random.Random(20261116+100*rep+dim).shuffle(ordered)
            for policy in ordered:
                work=root/name/phase/f"rep-{rep+1:02d}"/policy
                try:
                    meta=prepare(work,lengths,n,k,policy,dim,initial)
                    if name in physical_hashes and physical_hashes[name]!=meta["physical_input_sha256"]:
                        raise AssertionError("Physical input changed between policies")
                    physical_hashes[name]=meta["physical_input_sha256"]
                    row=execute(exe,work,meta,dim,reference,a.residual_atol,
                                a.timeout,[cpu],phase,rep+1,policy,name)
                    allrows.append(row)
                    print("CACHE FEW MAGNON",name,phase,policy,"dim",dim,
                          "wall",round(row["wall_seconds"],5),"RSS",row["peak_rss_kib"],
                          "Hcalls",row["h_applications"],"resid",row["direct_residual"],
                          "row_bound",row["reported_optimal_mnte"],"PASS",flush=True)
                except Exception as e:
                    failure=dict(case=name,phase=phase,repeat=rep+1,policy=policy,
                                 error=str(e),kind=type(e).__name__)
                    violations.append(failure)
                    print("CACHE FEW MAGNON FAILURE",failure,flush=True)
    matrix=[];comparison=[]
    for name in cases:
        for policy in policies:
            vals=[x for x in allrows if x["case"]==name and x["policy"]==policy
                  and x["phase"]=="measurement"]
            matrix.append(dict(case=name,policy=policy,
                count=len(vals),expected=a.repeats,
                wall_s=stats([x["wall_seconds"] for x in vals]),
                maxrss_kib=stats([x["peak_rss_kib"] for x in vals]),
                h_calls=stats([x["h_applications"] for x in vals]),
                cache_bytes=vals[0]["nominal_cache_bytes"] if vals else None))
        sample=[x for x in allrows if x["case"]==name and x["phase"]=="measurement"]
        counts={x["h_applications"] for x in sample}
        energies=[x["energy"] for x in sample]
        initials={x["initial_sha256"] for x in sample}
        consistency=(len(counts)==1 and len(initials)==1 and
                     (max(energies)-min(energies) <=2e-7 if energies else False))
        if not consistency:
            violations.append(dict(case=name,error="Matched initial/H-actions/energy contract differs",
                                   H_application_counts=sorted(counts),energy_range=
                                   (max(energies)-min(energies)) if energies else None))
        no=[x["wall_seconds"] for x in sample if x["policy"]=="none"]
        full=[x["wall_seconds"] for x in sample if x["policy"]=="full_tight"]
        speedup=statistics.median(no)/statistics.median(full) if no and full else None
        comparison.append(dict(case=name,matched_h_applications=sorted(counts),
                               consistent=consistency,ratio_no_to_full_exploratory=speedup))
    expected=len(cases)*len(policies)*(a.warmups+a.repeats)
    done=(len(allrows)==expected and not violations)
    env=dict(hostname=os.uname().nodename,linux_release=os.uname().release,
             selected_cpu=cpu,original_affinity=allowed,
             current_affinity=sorted(os.sched_getaffinity(0)),
             configured_threads=1,OMP_DYNAMIC=False,
             source_sha=os.environ.get("GITHUB_SHA","not-available"),
             executable_sha256=sha_file(exe),
             mode=("shared_or_exclusivity_unverified" if a.allow_shared else
                   "operator_attested_exclusive_not_machine_verified"))
    report=dict(schema=1,scope="identical ED2 XXZ near-saturation cache ON/OFF",
                parameters=dict(cases=cases,policies=policies,warmups=a.warmups,
                   repeats=a.repeats,timeout=a.timeout,residual_atol=a.residual_atol,
                   seed=SEED),environment=env,physical_input_sha256=physical_hashes,
                observations=allrows,summary=matrix,paired_comparisons=comparison,
                failures=violations,expected=expected,completed=len(allrows),
                all_numerical_and_coverage_checks_passed=done,
                publication_ready=bool(done and a.exclusive_allocation_verified
                    and a.repeats>=7),
                caution=("Shared-runner data are NOT controlled-node speedups; "
                         "independent exclusive host verification and repeat-node "
                         "variation remain necessary. Cache row bound applies "
                         "ONLY to this NN isotropic few-magnon model."))
    (root/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(f"CACHE FEW MAGNON COMPLETE {len(allrows)}/{expected} "
          f"violations={len(violations)} passed={done} publication_eligible="
          f"{report['publication_ready']}",flush=True)
    return 0 if done else 1

if __name__=="__main__":
    raise SystemExit(main())
