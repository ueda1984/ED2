#!/usr/bin/env python3
"""PRESTO N1000/m3 exact-A1g paired ED2/QuSpin driver (Python >=3.6).

Never runs on a login node without an actual PBS allocation. An execution
requires the complete group, numerical residual, model hash and a pinned
build of latest main. Every raw attempt is retained. XDiag is excluded
until an exact Gamma/Oh A1g adapter is validated; U1 is NOT the same block.
"""
import argparse
from datetime import datetime
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import random
import re
import resource
import shutil
import signal
import statistics
import subprocess
import sys
import time

MODEL_SHA = "a774358f0c69f81982d3495707f45433acd223b3a40aacd109304fb371903fd7"
REFERENCE_E0 = 847.3078940655731
SIDE, N, DOWN, DIM = 10, 1000, 3, 4299
POINTS = ((-1,-3,2), (-2,1,3), (2,3,1))  # ED2 PR16 checked 4,4,3 Oh factorization
POINT_ORDERS = (4,4,3)
SUPPORTED = ("ed2-cubic", "ed2-legacy", "quspin")


def sha(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for piece in iter(lambda: f.read(1024*1024), b""):
            h.update(piece)
    return h.hexdigest()


def write_json(path, obj):
    path.write_text(json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n")


def site(x,y,z):
    return x%SIDE + SIDE * ((y%SIDE) + SIDE*(z%SIDE))


def site_xyz(i):
    return (i%SIDE,(i//SIDE)%SIDE,i//(SIDE*SIDE))


def signed(mat, coords):
    transformed = [0,0,0]
    for axis, sign_axis in enumerate(mat):
        transformed[abs(sign_axis)-1] = (1 if sign_axis>0 else -1)*coords[axis]
    return transformed


def generators():
    # No spin reflection or approximate projection: exact spatial affine maps.
    gx = []
    gy = []
    gz = []
    for i in range(N):
        x,y,z=site_xyz(i)
        gx.append(site(x+1,y,z))
        gy.append(site(x,y+1,z))
        gz.append(site(x,y,z+1))
    result=[(gx,SIDE,0),(gy,SIDE,0),(gz,SIDE,0)]
    for mat,order in zip(POINTS,POINT_ORDERS):
        result.append(([site(*signed(mat,site_xyz(i))) for i in range(N)],
                      order,0))
    for p,order,_ in result:
        if sorted(p)!=list(range(N)):
            raise AssertionError("Invalid symmetry site permutation")
        current=list(range(N))
        for _ in range(order):
            current=[p[i] for i in current]
        if current!=list(range(N)):
            raise AssertionError("Incorrect symmetry generator order")
    return result


def model_bonds():
    axial=((1,0,0),(0,1,0),(0,0,1))
    faces=((1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1))
    ans=[]
    for z in range(SIDE):
        for y in range(SIDE):
            for x in range(SIDE):
                i=site(x,y,z)
                for dx,dy,dz in axial:
                    ans.append((i,site(x+dx,y+dy,z+dz),1.0,.77))
                for dx,dy,dz in faces:
                    ans.append((i,site(x+dx,y+dy,z+dz),.31,.19))
    if len(ans)!=9000 or len(set((min(i,j),max(i,j)) for i,j,_,_ in ans))!=9000:
        raise AssertionError("Invalid canonical bond list")
    return ans


def physical_sha(bonds):
    canonical=sorted((min(i,j),max(i,j),xy,zz) for i,j,xy,zz in bonds)
    data=dict(model="periodic_spin_half_cubic_J1J2_XXZ",
              shape=[SIDE]*3,ndown=DOWN,site_order="x-fastest",field=0.0,
              bonds=[list(x) for x in canonical],
              spin_operator_convention="SdotS uses spin (not Pauli) normalization")
    calculated=hashlib.sha256(json.dumps(
        data,sort_keys=True,separators=(",",":")).encode("utf8")).hexdigest()
    if calculated!=MODEL_SHA:
        raise AssertionError("Physical model SHA mismatches independent PR36 export")
    return calculated


def prepare_native(work):
    work.mkdir(parents=True,exist_ok=False)
    bonds=model_bonds()
    phys=physical_sha(bonds)
    gens=generators()
    for p,_,_ in gens:
        for i in (0,1,234,999):
            assert 0<=p[i]<N
    (work/"spins.dat").write_text("0.5\n"*N)
    (work/"caps.dat").write_text("1\n"*N)
    (work/"fields.dat").write_text("")
    (work/"bonds.dat").write_text("".join(
        "%d %d %.17g %.17g %.17g 0 0 0 0 0 0\n" %
        (i+1,j+1,xy,xy,zz) for i,j,xy,zz in bonds))
    lines=[
        "NOS=1000,NODmin=3,NODmax=3,ALG=1,",
        "NO_one=0,NO_two=9000,wk_dim=-1,MNTE=64,",
        "FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',",
        "FILE_spin='spins.dat',FILE_NODmax='caps.dat',",
        "OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0,"]
    for idx,(perm,order,phase) in enumerate(gens):
        name="symmetry%d.dat"%(idx+1)
        (work/name).write_text("".join("%d\n"%(k+1) for k in perm))
        lines.append("L%d=%d,M%d=%d,FILE_S%d='%s'," %
                     (idx+1,order,idx+1,phase,idx+1,name))
    (work/"input.dat").write_text(
        "&input_parameters\n" + "\n".join(lines) + "\n/\n"
        "&input_lancz\nlnc_ene_conv=1d-11,minitr=5,maxitr=2000,itrint=5\n/\n"
        "&input_TRLan\nNOE=2,NOK=3,NOM=8,maxitr=2000,"
        "lnc_ene_conv=1d-11,i_vec_min=1,i_vec_max=1\n/\n")
    return dict(physical_sha256=phys,
                input_sha256={p.name:sha(p) for p in sorted(work.iterdir())
                              if p.is_file()},expected_dimension=DIM,
                bonds=len(bonds))


def physical_cpus(nthreads):
    allowed=sorted(os.sched_getaffinity(0))
    by_core={}
    for cpu in allowed:
        base=Path("/sys/devices/system/cpu/cpu%d/topology"%cpu)
        key=(int((base/"physical_package_id").read_text()),
             int((base/"core_id").read_text()))
        by_core.setdefault(key,[]).append(cpu)
    selection=[min(by_core[k]) for k in sorted(by_core)][:nthreads]
    if len(selection)!=nthreads:
        raise RuntimeError("Not enough allowed distinct physical CPU cores")
    return selection


def launch(command,cwd,env,output,timeout,cpus,vmem_gib):
    if output.exists():
        raise RuntimeError("Do not overwrite a trial: "+str(output))
    timepath=cwd/"process-resource.txt"
    cmd=list(command)
    if Path("/usr/bin/time").is_file():
        cmd=["/usr/bin/time","-f","%e %M %U %S",
             "-o",str(timepath)] + cmd

    def private_setup():
        os.setsid()
        os.sched_setaffinity(0,set(cpus))
        mem=int(vmem_gib*(1024**3))
        resource.setrlimit(resource.RLIMIT_AS,(mem,mem))
    t0=time.time()
    timed_out=False
    with output.open("wb") as f:
        with (cwd/"input.dat").open("rb") if (cwd/"input.dat").exists() else open(os.devnull,"rb") as stdin:
            proc=subprocess.Popen(cmd,cwd=str(cwd),stdin=stdin,
                                  stdout=f,stderr=subprocess.STDOUT,
                                  env=env,preexec_fn=private_setup)
            try:
                code=proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out=True
                os.killpg(proc.pid,signal.SIGTERM)
                try: proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL)
                    proc.wait()
                code=124
    return dict(returncode=code,timeout=timed_out,
                full_process_wall_seconds=time.time()-t0,
                log_sha256=sha(output),
                resource_text=timepath.read_text().strip()
                              if timepath.exists() else None)


def parse_native(work):
    stdout=(work/"stdout.log").read_text(errors="replace")
    found=re.findall(r"THS\(k\)\s*=\s*(\d+)",stdout)
    if len(found)!=1 or int(found[0])!=DIM:
        raise RuntimeError("ED2 sector dimension missing/mismatch")
    eigen=(work/"eigenvalues.dat")
    if not eigen.exists(): raise RuntimeError("ED2 eigenvalues.dat absent")
    energies=[float(x.split()[-1].replace("D","E"))
              for x in eigen.read_text().splitlines() if x.strip()]
    if len(energies)!=1 or not math.isfinite(energies[0]):
        raise RuntimeError("ED2 eigenvalue missing or nonfinite")
    lines=(work/"diagnostics.txt").read_text().splitlines()
    diag=dict(x.split("=",1) for x in lines if "=" in x)
    if diag.get("lanczos_mode","").strip()!="two_stage":
        raise RuntimeError("ED2 two-stage certified Lanczos did not run")
    if sum(x.startswith("pair_columns=") for x in lines)!=1:
        raise RuntimeError("ED2 residual schema marker absent")
    pairs=[x.split() for x in lines if re.match(r"^\s*1\s+[TF]\s+[TF]\s+",x)]
    if len(pairs)!=1 or len(pairs[0])!=8 or pairs[0][1:3]!=["T","T"]:
        raise RuntimeError("ED2 direct Ritz pair not certified")
    energy=energies[0]; residual=float(pairs[0][5])
    if abs(energy-REFERENCE_E0)>1e-8 or not math.isfinite(residual) or residual>1e-8:
        raise RuntimeError("ED2 E0/direct residual differs from independent oracle")
    stages={}
    for name,val in re.findall(
            r"ED2_PHASE ([a-z_]+) wall_s=\s*([+\-0-9.eE]+)",stdout):
        if name in stages: raise RuntimeError("Duplicate ED2 phase")
        stages[name]=float(val)
    needed={"input","combinatorics","basis","operator_cache","lanczos",
            "observables_and_output","total"}
    if set(stages)!=needed or not all(math.isfinite(x) and x>=0
                                       for x in stages.values()):
        raise RuntimeError("Incomplete ED2 phase timing")
    rss=None
    line=(work/"process-resource.txt")
    if line.is_file():
        words=line.read_text().split()
        if len(words)>=2: rss=int(words[1])
    return dict(status="PASS",dimension=DIM,energy=energy,residual=residual,
                h_applications=int(diag["h_applications"]),phases=stages,
                max_rss_kib=rss,diagnostic_sha256=sha(work/"diagnostics.txt"))


def parse_quspin(work):
    log=(work/"stdout.log").read_text()
    data=json.loads(log.strip().splitlines()[-1])
    if (data.get("status")!="PASS" or
        data.get("physical_sha256")!=MODEL_SHA or
        data.get("dimension_constructed")!=DIM or
        data.get("sector")!="gamma-a1g" or
        data.get("model_choice")!="cubic-j1j2-xxz"):
        raise RuntimeError("QuSpin input/sector/eigensolver status mismatch")
    energy=float(data["E0"]); residual=float(data["direct_residual"])
    if (abs(energy-REFERENCE_E0)>1e-8 or
        not math.isfinite(residual) or residual>1e-8):
        raise RuntimeError("QuSpin E0/direct residual outside physics tolerance")
    return dict(status="PASS",dimension=DIM,energy=energy,residual=residual,
                phases=dict(basis=data["basis_seconds"],
                            operator_cache=data["csr_seconds"],
                            lanczos=data["solve_seconds"],
                            total=data["time_to_eigenpair_seconds"]),
                max_rss_kib=data.get("peak_rss_kib"),
                raw_quspin=data)


def append(path,obj):
    with path.open("a") as f:
        f.write(json.dumps(obj,allow_nan=False,sort_keys=True)+"\n")
        f.flush()
        os.fsync(f.fileno())


def run(args):
    if not os.environ.get("PBS_JOBID") or os.environ.get("ED2_EXCLUSIVE_CONFIRMED")!="1":
        raise RuntimeError("Execution requires approved real PBS allocation")
    if not args.out.is_dir(): raise RuntimeError("Output directory missing")
    programs=args.programs.split(",")
    if not programs or len(set(programs))!=len(programs) or any(
            p not in SUPPORTED for p in programs):
        raise ValueError("Only "+",".join(SUPPORTED)+" allowed")
    if args.warmups<0 or args.repeats<1 or args.repeats>10:
        raise ValueError("Invalid repeats/warmups")
    if args.threads not in (1,4,8):
        raise ValueError("Pilot thread count restricted to 1/4/8")
    if "quspin" in programs:
        if args.quspin_python is None or not args.quspin_python.is_file():
            raise RuntimeError("QuSpin v1.0.1 python environment not configured")
    if any(x.startswith("ed2") for x in programs) and (
            args.exe is None or not args.exe.is_file()):
        raise RuntimeError("Latest ED2 executable not configured")
    env=dict(os.environ)
    env.update(OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close", OMP_PLACES="cores",
               OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1",BLIS_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1",
               ED2_RNG_SEED="1984")
    cpus=physical_cpus(args.threads)
    manifest=dict(pbs_jobid=os.environ["PBS_JOBID"],threads=args.threads,
                  cpu_ids=cpus,programs=programs,
                  expected_energy=REFERENCE_E0,expected_dimension=DIM,
                  physical_sha256=MODEL_SHA,
                  ed2_exe_sha256=sha(args.exe) if args.exe else None,
                  quspin_python=str(args.quspin_python) if args.quspin_python else None,
                  benchmark_py_sha256=sha(args.bench/"run_quspin.py"),
                  warmups=args.warmups,repeats=args.repeats,
                  timeout_seconds=args.timeout_seconds,
                  virtual_memory_gib=args.max_virtual_gib,
                  requires_independent_exclusivity_attestation=True,
                  publication_ready=False)
    write_json(args.out/"plan.json",manifest)
    trial_file=args.out/"trials.jsonl"
    result=[]
    for n in range(args.warmups+args.repeats):
        period="warmup" if n<args.warmups else "measured"
        order=list(programs)
        random.Random(1984+n).shuffle(order)
        for kind in order:
            key="%s-r%d-%s"%(period,n,kind)
            work=args.out/key
            meta=prepare_native(work) if kind.startswith("ed2") else None
            if meta is None:
                work.mkdir(parents=True,exist_ok=False)
            env_i=dict(env)
            if kind=="ed2-cubic":
                env_i.update(ED2_CUBIC_FAST_GROUP_VALIDATION="1",
                             ED2_CUBIC_ANCHORED_REP="1",
                             ED2_CUBIC_ANCHORED_BASIS="1",
                             ED2_CUBIC_LOCAL_BOND_CACHE="1",
                             ED2_CUBIC_DIRECT_HOP_TABLE="1")
            elif kind=="ed2-legacy":
                env_i.pop("ED2_CUBIC_FAST_GROUP_VALIDATION",None)
                env_i.pop("ED2_CUBIC_ANCHORED_REP",None)
                env_i.pop("ED2_CUBIC_ANCHORED_BASIS",None)
                env_i.pop("ED2_CUBIC_LOCAL_BOND_CACHE",None)
                env_i.pop("ED2_CUBIC_DIRECT_HOP_TABLE",None)
                env_i["ED2_ALLOW_LARGE_GROUP_VALIDATION"]="1"
            if kind.startswith("ed2"):
                env_i.update(ED2_LANCZOS_TWO_STAGE="1",ED2_PHASE_TIMING="1",
                             ED2_DIAGNOSTICS="diagnostics.txt",
                             ED2_RESIDUAL_ATOL="1e-8")
                cmd=[str(args.exe)]
            else:
                cmd=[str(args.quspin_python),str(args.bench/"run_quspin.py"),
                     "--shape","10,10,10","--ndown","3",
                     "--model","cubic-j1j2-xxz","--sector","gamma-a1g",
                     "--mode","solve","--max-raw-dim","200000000",
                     "--max-csr-nnz","1000000","--nup-convention","negative",
                     "--prefix-bits","0","--tol","1e-10"]
            rec=dict(key=key,kind=kind,period=period,repeat=n,
                     physical_sha256=MODEL_SHA,
                     input_files=meta["input_sha256"] if meta else None,
                     start_utc=datetime.utcnow().isoformat()+"Z")
            try:
                runtime=launch(cmd,work,env_i,work/"stdout.log",
                               args.timeout_seconds,cpus,args.max_virtual_gib)
                rec.update(runtime)
                if runtime["returncode"]!=0:
                    raise RuntimeError("Exit=%d timeout=%s"%(runtime["returncode"],
                                                            runtime["timeout"]))
                result_phys=parse_native(work) if kind.startswith("ed2") else parse_quspin(work)
                rec.update(result_phys)
            except Exception as exc:
                rec.update(status="FAILED",error=str(exc))
                append(trial_file,rec)
                raise
            append(trial_file,rec)
            result.append(rec)
            print("PASS %s E=%.13f residual=%g wall=%.3fs" %
                  (kind,rec["energy"],rec["residual"],
                   rec["full_process_wall_seconds"]),flush=True)
    summary={}
    for kind in programs:
        x=[r["full_process_wall_seconds"] for r in result
           if r["kind"]==kind and r["period"]=="measured"]
        med=statistics.median(x)
        summary[kind]=dict(n=len(x),median_wall_s=med,
                           median_abs_deviation_s=statistics.median(
                               [abs(y-med) for y in x]))
    write_json(args.out/"summary.json",dict(
        complete=True,pbs_jobid=os.environ["PBS_JOBID"],
        source_sha=os.environ.get("ED2_MAIN_COMMIT",""),
        experiment_summary=summary,
        warnings=["Not publication-ready without independent PRESTO resource attestation.",
                  "Cross-library time uses same hardware, but solver algorithms differ.",
                  "No XDiag A1g result because equivalent-symmetry adapter is unavailable."]))
    return 0


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--plan-only",action="store_true")
    p.add_argument("--exe",type=Path)
    p.add_argument("--quspin-python",type=Path)
    p.add_argument("--bench",type=Path,default=Path("benchmarks/cross_library"))
    p.add_argument("--out",type=Path,default=Path("."))
    p.add_argument("--programs",default="ed2-cubic,quspin")
    p.add_argument("--threads",type=int,default=1)
    p.add_argument("--warmups",type=int,default=0)
    p.add_argument("--repeats",type=int,default=1)
    p.add_argument("--timeout-seconds",type=int,default=1800)
    p.add_argument("--max-virtual-gib",type=int,default=48)
    a=p.parse_args()
    if a.plan_only:
        print(json.dumps(dict(status="PLAN_ONLY",required_nsites=N,
                              magnons=DOWN,dimension=DIM,
                              physical_sha256=physical_sha(model_bonds()),
                              programs=a.programs,
                              required_latest_main_flags=[
                               "ED2_CUBIC_FAST_GROUP_VALIDATION",
                               "ED2_CUBIC_ANCHORED_REP",
                               "ED2_CUBIC_ANCHORED_BASIS",
                               "ED2_CUBIC_LOCAL_BOND_CACHE",
                               "ED2_CUBIC_DIRECT_HOP_TABLE",
                               "ED2_LANCZOS_TWO_STAGE",
                               "ED2_PHASE_TIMING","ED2_DIAGNOSTICS"],
                              note="No native/QuSpin execution performed"),
                         sort_keys=True))
        return 0
    return run(a)


if __name__=="__main__":
    sys.exit(main())
