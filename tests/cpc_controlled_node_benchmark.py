#!/usr/bin/env python3
"""Single-node ED2 reproducibility orchestrator for the paired-cache study.

This wrapper uses the *existing independently referenced*
cpc_stationary_cache_pareto.py without changing the ED2 solver kernel.
It pins its own CPU affinity (inherited by ED2 children), disables BLAS
oversubscription, executes identical workloads at requested OMP counts,
and records all machine/source/binary/provenance and raw case outputs.

Exclusivity cannot be certified programmatically in a portable manner:
the operator must explicitly attest --exclusive-allocation-verified.
For functional CI smoke only, use --allow-shared, which labels the
results as exploratory and NOT publication-grade speedup evidence.
Never infer exclusivity just from a PBS/Slurm job ID.

Example dedicated node after operator verifies an exclusive allocation:
  python3 tests/cpc_controlled_node_benchmark.py \
    --exe /absolute/path/to/QS3ED2 --output /new/results \
    --cpu-list 0-7 --threads 1,2,4,8 \
    --sizes 8,10,12 --algs 1,2 \
    --repeats 7 --warmups 2 --exclusive-allocation-verified
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
BENCH=HERE/"cpc_stationary_cache_pareto.py"
SAFE_SCHEDULER_KEYS=(
    "SLURM_JOB_ID","SLURM_JOB_NODELIST","SLURM_JOB_NUM_NODES",
    "SLURM_CPUS_ON_NODE","SLURM_JOB_CPUS_PER_NODE",
    "PBS_JOBID","PBS_QUEUE","LSB_JOBID","LSB_QUEUE",
)


def cpu_list(value):
    ids=[]
    for token in value.split(","):
        token=token.strip()
        if re.fullmatch(r"\d+",token):
            ids.append(int(token))
        elif re.fullmatch(r"\d+-\d+",token):
            lo,hi=map(int,token.split("-"))
            if lo>hi or hi-lo>2048:
                raise ValueError("Invalid CPU index range")
            ids.extend(range(lo,hi+1))
        else:
            raise ValueError("Invalid CPU-list token "+token)
    if not ids or len(set(ids))!=len(ids):
        raise ValueError("Empty or duplicate CPU-list entries")
    return sorted(ids)


def integer_csv(text,allowed,field):
    out=[int(x.strip()) for x in text.split(",")]
    if not out or len(set(out))!=len(out) or any(x not in allowed for x in out):
        raise ValueError("Invalid or repeated "+field)
    return out


def sha256(path):
    digest=hashlib.sha256()
    with path.open("rb") as src:
        while True:
            chunk=src.read(1024*1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def command_output(cmd):
    try:
        out=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,timeout=15,check=False)
        return dict(exit_code=out.returncode,output=out.stdout[:30000])
    except (OSError,subprocess.TimeoutExpired) as err:
        return dict(error=str(err))


def record_environment(binary,attestation,selected,allowed,threads,params,
                       source_revision=None,source_archive=None):
    return dict(
        schema=1,utc_epoch_time=time.time(),
        python=sys.version,platform=platform.platform(),
        uname=list(platform.uname()),hostname=platform.node(),
        pid=os.getpid(),euid=os.geteuid() if hasattr(os,"geteuid") else None,
        benchmark_mode=("operator_verified_exclusive" if attestation else
                        "shared_or_exclusivity_unverified"),
        exclusive_operator_attestation=attestation,
        allowed_affinity=allowed,selected_affinity=selected,
        process_affinity=(sorted(os.sched_getaffinity(0))
                          if hasattr(os,"sched_getaffinity") else None),
        requested_omp_counts=threads,
        known_scheduler_fields={key:os.environ[key]
                                for key in SAFE_SCHEDULER_KEYS if key in os.environ},
        relevant_env={key:os.environ.get(key)
                      for key in ("OMP_NUM_THREADS","OMP_DYNAMIC","OMP_PROC_BIND",
                                  "OMP_PLACES","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
                                  "GOTO_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS")},
        physical_cpu=command_output(["lscpu"]),
        memory=command_output(["free","-b"]),
        compiler=command_output(["gfortran","--version"]),
        compiler13=command_output(["gfortran-13","--version"]),
        linked_libraries=command_output(["ldd",str(binary)]),
        numactl=command_output(["numactl","--hardware"]),
        python_numpy_scipy=command_output([sys.executable,"-c",
                                          "import numpy,scipy; print(numpy.__version__,scipy.__version__)"]),
        git_revision=command_output(["git","-C",str(HERE.parent),"rev-parse","HEAD"]),
        executable_path=str(binary),executable_sha256=sha256(binary),
        benchmark_script_sha256=sha256(BENCH),
        source_revision=source_revision,
        source_archive_path=str(source_archive) if source_archive else None,
        source_archive_sha256=sha256(source_archive) if source_archive else None,
        wrapper_sha256=sha256(Path(__file__).resolve()),
        parameters=params,
        caution=("Operator exclusivity is self-attested, not verified by code. "
                 "Other workloads, turbo/governor, NUMA or thermal effects require "
                 "separate assessment."))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--cpu-list",required=True,
                   help="Explicit CPU affinity, e.g. 0-7,16-23 (not CPU count)")
    p.add_argument("--threads",default="1,2")
    p.add_argument("--sizes",default="8,10,12")
    p.add_argument("--algs",default="1,2")
    p.add_argument("--D",type=float,default=0.12)
    p.add_argument("--repeats",type=int,default=7)
    p.add_argument("--warmups",type=int,default=2)
    p.add_argument("--keep",type=int,default=2)
    p.add_argument("--capacity",type=int,default=16)
    p.add_argument("--residual-atol",type=float,default=1e-8)
    p.add_argument("--source-revision",help="Full 40- or 64-digit source revision SHA")
    p.add_argument("--source-archive",type=Path,
                   help="Exact source tarball used for the executable build")
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--exclusive-allocation-verified",action="store_true",
                      help="Operator attests to a real exclusive node reservation")
    mode.add_argument("--allow-shared",action="store_true",
                      help="Only for CI smoke/exploration; not publication timing")
    a=p.parse_args()
    try:
        threads=integer_csv(a.threads,set(range(1,129)),"OMP thread counts")
        sizes=integer_csv(a.sizes,{4,6,8,10,12},"system sizes")
        algs=integer_csv(a.algs,{1,2},"solver algorithms")
        cpus=cpu_list(a.cpu_list)
        if len(cpus)<max(threads):
            raise ValueError("Selected CPU affinity smaller than requested OMP count")
        if not -0.5<=a.D<=0.5 or not 0<a.residual_atol<=1e-5:
            raise ValueError("DM parameter/residual tolerance outside validated range")
        if not 1<=a.repeats<=20 or not 0<=a.warmups<=5:
            raise ValueError("Invalid number of measurements or warmups")
        if not 1<a.keep<a.capacity<=512:
            raise ValueError("Invalid ALG2 restart workspace")
        if a.source_revision and not re.fullmatch(r"[0-9a-fA-F]{40}([0-9a-fA-F]{24})?",a.source_revision):
            raise ValueError("Source revision must be a full 40/64-hex SHA")
        if a.exclusive_allocation_verified and not (a.source_revision or a.source_archive):
            raise ValueError("Publication-mode runs require a source SHA or source archive")
        source_archive=(a.source_archive.resolve(strict=True) if a.source_archive else None)
        if source_archive is not None and not source_archive.is_file():
            raise ValueError("Source archive not a regular file")
        if not hasattr(os,"sched_getaffinity") or not hasattr(os,"sched_setaffinity"):
            raise ValueError("Linux sched_getaffinity/setaffinity required for controlled CPU pinning")
        allowed=sorted(os.sched_getaffinity(0))
        if not set(cpus).issubset(allowed):
            raise ValueError("CPU set not allowed by scheduler/OS; allowed "+str(allowed))
    except ValueError as err:
        p.error(str(err))
    exe=a.exe.resolve(strict=True)
    if not exe.is_file():
        p.error("ED2 executable is not a regular file")
    if not BENCH.is_file():
        p.error("Missing paired-cache benchmark source")
    output=a.output.resolve()
    output.mkdir(exist_ok=False,parents=True)
    os.sched_setaffinity(0,set(cpus))
    os.environ.update(
        OMP_DYNAMIC="FALSE",OMP_PROC_BIND="TRUE",OMP_PLACES="threads",
        OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
        GOTO_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1")
    params=dict(sizes=sizes,algs=algs,D=a.D,threads=threads,
                repeats=a.repeats,warmups=a.warmups,
                keep=a.keep,capacity=a.capacity,
                residual_atol=a.residual_atol)
    metadata=record_environment(exe,a.exclusive_allocation_verified,
                                cpus,allowed,threads,params,
                                source_revision=a.source_revision,
                                source_archive=source_archive)
    (output/"environment.json").write_text(json.dumps(metadata,indent=2)+"\n")
    tasks=[]
    for omp in threads:
        childroot=output/f"paired-omp{omp}"
        args=[sys.executable,str(BENCH),"--exe",str(exe),
              "--output",str(childroot),"--sizes",",".join(map(str,sizes)),
              "--algs",",".join(map(str,algs)),"--D",str(a.D),
              "--omp",str(omp),"--repeats",str(a.repeats),
              "--warmups",str(a.warmups),"--residual-atol",str(a.residual_atol),
              "--keep",str(a.keep),"--capacity",str(a.capacity)]
        start=time.perf_counter()
        try:
            with (output/f"run-omp{omp}.log").open("w") as log:
                run=subprocess.run(args,cwd=HERE,env=dict(os.environ),
                                   stdout=log,stderr=subprocess.STDOUT,
                                   check=False)
            exitcode=run.returncode
            reportpath=childroot/"report.json"
            report=json.loads(reportpath.read_text()) if reportpath.exists() else None
            tasks.append(dict(omp=omp,exit_code=exitcode,
                              wall_seconds=time.perf_counter()-start,
                              report_path=str(reportpath),
                              completed=(report or {}).get("completed"),
                              expected=(report or {}).get("expected"),
                              pass_complete=(report or {}).get("pass_complete",False),
                              violations=(report or {}).get("violations")))
            print("NODE BENCH",omp,"threads","exit",exitcode,
                  "complete",tasks[-1]["pass_complete"],flush=True)
        except Exception as exc:
            tasks.append(dict(omp=omp,error=str(exc),exit_code=None,pass_complete=False))
            print("NODE BENCH FAILED",omp,str(exc),flush=True)
    good=len(tasks)==len(threads) and all(
        x.get("exit_code")==0 and x.get("pass_complete") for x in tasks)
    manifest=[]
    for folder in sorted(output.glob("paired-omp*")):
        for name in ("input.dat","bonds.dat","fields.dat","spins.dat",
                     "caps.dat","pairs.dat","initial-vector.txt"):
            for item in sorted(folder.rglob(name)):
                manifest.append(dict(relative_path=str(item.relative_to(output)),
                                     sha256=sha256(item),size_bytes=item.stat().st_size))
    (output/"input-manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    if good and not manifest:
        good=False
    summary=dict(schema=1,metadata_sha256=sha256(output/"environment.json"),
                 input_manifest_sha256=sha256(output/"input-manifest.json"),
                 input_file_count=len(manifest),
                 mode=metadata["benchmark_mode"],
                 raw_measurements="paired-omp*/report.json contains all raw cases",
                 tasks=tasks,pass_complete=good)
    (output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("NODE BENCHMARK SUMMARY",sum(x.get("pass_complete",False) for x in tasks),
          "/",len(threads),"thread configurations accepted",flush=True)
    return 0 if good else 1

if __name__=="__main__":
    raise SystemExit(main())
