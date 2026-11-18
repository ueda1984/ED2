#!/usr/bin/env python3
"""PRESTO *cubic O_h* performance run: pinned physical cores, immutable trials.

No scheduler mutations or job submission. A PBS exclhost request does not itself
prove absence of another job; publication attestation remains pending.
Python standard library only; imports the existing independent ED2 test oracle.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import socket
import statistics
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import cpc_pointgroup_profile as physics
import cpc_cubic_pointgroup as cubic

SECTORS = cubic.SECTORS
PHASES = ("basis", "operator_cache", "lanczos", "total")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def capture(command, timeout=15):
    try:
        cp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, timeout=timeout)
        return dict(command=command, returncode=cp.returncode,
                    stdout=cp.stdout[-200000:], stderr=cp.stderr[-20000:])
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        return dict(command=command, error=str(exc))


def core_map(allowed):
    cores = {}
    for cpu in sorted(allowed):
        base = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        try:
            package = int((base / "physical_package_id").read_text())
            core = int((base / "core_id").read_text())
        except (FileNotFoundError, ValueError):
            raise RuntimeError(f"Cannot verify physical core topology for CPU {cpu}")
        cores.setdefault((package, core), []).append(cpu)
    return cores


def parse_cpus(text):
    cpus = set()
    for word in text.split(","):
        word = word.strip()
        if "-" in word:
            lo, hi = map(int, word.split("-", 1))
            if lo > hi or hi - lo > 4096:
                raise ValueError("Invalid CPU range")
            cpus.update(range(lo, hi + 1))
        else:
            cpus.add(int(word))
    return sorted(cpus)


def allocation(shared, selected_threads, cpu_list):
    allowed = set(os.sched_getaffinity(0))
    cores = core_map(allowed)
    if cpu_list:
        selected = parse_cpus(cpu_list)
        if not set(selected) <= allowed:
            raise RuntimeError("Requested CPU list is not contained in allowed affinity")
        chosen_cores = {(k[0], k[1]) for k, xs in cores.items() for x in xs if x in selected}
        if len(selected) != len(chosen_cores):
            raise RuntimeError("CPU list includes logical siblings of a physical core")
    else:
        selected = [sorted(xs)[0] for _, xs in sorted(cores.items())][:max(selected_threads)]
    if len(selected) < max(selected_threads):
        raise RuntimeError("Insufficient distinct physical cores for requested threads")
    # All ED2 children inherit an identical physical-core affinity. OpenMP
    # chooses its first N cores from that same fixed allowed set.
    os.sched_setaffinity(0, selected)
    if set(os.sched_getaffinity(0)) != set(selected):
        raise RuntimeError("Failed to pin process affinity")
    pbs_job = os.environ.get("PBS_JOBID", "")
    nodefile = os.environ.get("PBS_NODEFILE", "")
    qstat = capture(["qstat", "-f", pbs_job]) if pbs_job else {}
    pbsnodes = capture(["pbsnodes", "-av"]) if pbs_job else {}
    node_entries = []
    if nodefile:
        p = Path(nodefile)
        if not p.is_file():
            raise RuntimeError("PBS_NODEFILE is absent")
        node_entries = [v.strip() for v in p.read_text().splitlines() if v.strip()]
    if not shared:
        if not pbs_job or not node_entries or qstat.get("returncode") != 0:
            raise RuntimeError("Exclusive mode requires a readable real PBS allocation and qstat -f")
        hosts = {h.split(".")[0] for h in node_entries}
        if len(hosts) != 1 or socket.gethostname().split(".")[0] not in hosts:
            raise RuntimeError("PBS_NODEFILE does not describe this one exclusive node")
        qs = qstat["stdout"]
        if not re.search(r"Resource_List\.place\s*=\s*exclhost", qs):
            raise RuntimeError("PBS scheduler did not report place=exclhost")
        if not re.search(r"Resource_List\.select\s*=\s*1:ncpus=96\b", qs):
            raise RuntimeError("Expected a one-node CPU selection")
        if not re.search(r"\bqueue\s*=\s*standard\b", qs):
            raise RuntimeError("Not running in the expected standard queue")
        if "exec_vnode" not in qs:
            raise RuntimeError("Missing actual exec_vnode allocation record")
    return dict(mode="shared_smoke" if shared else "pbs_exclhost_requested_attestation_pending",
                host=socket.gethostname(), pbs_jobid=pbs_job,
                requested_nodefile=nodefile, nodefile_hosts=sorted(set(node_entries)),
                qstat=qstat, pbsnodes=pbsnodes,
                allowed_cpu_ids=sorted(allowed), selected_cpu_ids=selected,
                selected_core_topology=[dict(cpu=cpu, package=package, physical_core=core)
                  for (package, core), xs in sorted(cores.items())
                  for cpu in xs if cpu in selected],
                max_threads=max(selected_threads),
                commands={x: capture([x, "--version"]) for x in ("gfortran", "python3")},
                cpu_info={x: capture(x.split()) for x in
                          ("lscpu", "lscpu -e=CPU,CORE,SOCKET,NODE", "free -h", "uptime")},
                attestation="PENDING: a site operator must independently attest exclusive ownership")


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def medmad(data):
    m = statistics.median(data)
    return dict(n=len(data), median=m,
                mad=statistics.median([abs(v - m) for v in data]),
                minimum=min(data), maximum=max(data))


def append_trial(path, rec):
    with open(path, "a", encoding="utf8") as f:
        f.write(json.dumps(rec, allow_nan=False) + "\n")
        f.flush()
        os.fsync(f.fileno())


def execute(exe, out, spec, allocation_info, timeout, mnte):
    side, sector, threads, repeat, period, builder = spec
    key = f"s{side}-{sector}-t{threads}-{period}{repeat}-{builder}"
    work = out / "trials" / key
    generators = cubic.symmetry_generators(side, sector)
    bonds = cubic.canonical_bonds(side)
    assert cubic.check_hamiltonian(side)
    dimension = cubic.sector_dimensions(side)[sector]
    cubic.write_input(work, side, 3, generators, bonds, 2000, -1, mnte=mnte)
    before = os.environ.get("ED2_FORCE_LEGACY_ORBITS")
    before_local = os.environ.get("ED2_CUBIC_LOCAL_BOND_CACHE")
    before_profile = os.environ.get("ED2_CUBIC_CACHE_PROFILE")
    before_direct = os.environ.get("ED2_CUBIC_DIRECT_HOP_TABLE")
    try:
        if builder == "forced-legacy":
            os.environ["ED2_FORCE_LEGACY_ORBITS"] = "1"
            os.environ.pop("ED2_CUBIC_LOCAL_BOND_CACHE", None)
            os.environ.pop("ED2_CUBIC_CACHE_PROFILE", None)
            os.environ.pop("ED2_CUBIC_DIRECT_HOP_TABLE", None)
        else:
            os.environ.pop("ED2_FORCE_LEGACY_ORBITS", None)
        obs = physics.run_case(exe, work, side ** 3, 3, sector, dimension,
                               threads, "two-stage", timeout)
    finally:
        if before_direct is None:
            os.environ.pop("ED2_CUBIC_DIRECT_HOP_TABLE", None)
        else:
            os.environ["ED2_CUBIC_DIRECT_HOP_TABLE"] = before_direct
        if before_profile is None:
            os.environ.pop("ED2_CUBIC_CACHE_PROFILE", None)
        else:
            os.environ["ED2_CUBIC_CACHE_PROFILE"] = before_profile
        if before_local is None:
            os.environ.pop("ED2_CUBIC_LOCAL_BOND_CACHE", None)
        else:
            os.environ["ED2_CUBIC_LOCAL_BOND_CACHE"] = before_local
        if before is None:
            os.environ.pop("ED2_FORCE_LEGACY_ORBITS", None)
        else:
            os.environ["ED2_FORCE_LEGACY_ORBITS"] = before
    output = (work / "stdout.log").read_text(errors="replace")
    anchored = os.environ.get("ED2_CUBIC_ANCHORED_REP", "") == "1"
    anchor_verified = os.environ.get("ED2_CUBIC_ANCHORED_VERIFY", "") == "1"
    fast_group = os.environ.get("ED2_CUBIC_FAST_GROUP_VALIDATION", "") == "1"
    fast_group_reported = "Symmetry certificate: canonical cubic O_h semidirect product verified" in output
    if fast_group_reported != (fast_group and builder == "fast"):
        raise AssertionError("Fast cubic group certificate mismatch for selected native basis builder")
    anchored_reported = "Representative lookup: cubic-anchored-orbits" in output
    if anchored_reported != anchored:
        raise AssertionError("Requested anchored lookup disagrees with native stdout")
    if anchored and anchor_verified and "per-hop legacy verification=T" not in output:
        raise AssertionError("Missing per-hop legacy verification confirmation")
    anchored_basis_requested = os.environ.get("ED2_CUBIC_ANCHORED_BASIS", "") == "1"
    anchored_basis_used = "Basis builder: anchored-cubic" in output
    atomic = "Basis builder: atomic-orbits" in output
    if atomic and anchored_basis_used:
        raise AssertionError("Two incompatible basis builders were logged")
    if anchored_basis_used != (builder == "fast" and anchored_basis_requested):
        raise AssertionError("Requested anchored cubic basis not used exactly for fast variant")
    local_bond_cache = os.environ.get("ED2_CUBIC_LOCAL_BOND_CACHE", "") == "1"
    local_cache_reported = "Cubic cache builder: local-bond-sparse" in output
    if local_cache_reported != (local_bond_cache and builder == "fast"):
        raise AssertionError("Native local-bond cache builder disagrees with requested fast variant")
    cache_profile = os.environ.get("ED2_CUBIC_CACHE_PROFILE", "") == "1"
    direct_hop_table = os.environ.get("ED2_CUBIC_DIRECT_HOP_TABLE", "") == "1"
    direct_hop_reported = "Representative lookup: cubic-direct-hop-pairs" in output
    if direct_hop_reported != (direct_hop_table and builder == "fast"):
        raise AssertionError("Native direct-hop table disagrees with requested fast path")
    profile_reported = "Cubic cache profile: representative_key_setup_wall_s=" in output
    local_profile_reported = "Cubic cache profile: resolver_calls=" in output
    if profile_reported != (cache_profile and builder == "fast"):
        raise AssertionError("Native cache lookup profile marker disagreement")
    if local_profile_reported != (cache_profile and builder == "fast"):
        raise AssertionError("Native local cache profile marker disagreement")
    normal_form = "Symmetry action: noncommuting normal form and character validated" in output
    if builder == "forced-legacy" and atomic:
        raise AssertionError("Forced-legacy run unexpectedly used atomic orbit basis")
    if builder == "fast" and not (atomic or anchored_basis_used):
        raise AssertionError("Fast run unexpectedly bypassed both validated basis builders")
    if builder == "fast" and sector != "gamma-t" and not normal_form:
        raise AssertionError("Noncommuting group fast path lacked validation record")
    obs.update(dict(key=key, side=side, sector=sector, threads=threads,
                    repetition=repeat, period=period, builder=builder,
                    atomic_orbits=atomic, anchored_basis=anchored_basis_requested,
                    native_basis_builder=("anchored-cubic" if anchored_basis_used else
                                          "atomic-orbits" if atomic else "forced-legacy"),
                    validated_noncommuting_action=normal_form,
                    local_bond_cache=local_bond_cache,
                   cache_profile=cache_profile, direct_hop_table=direct_hop_table,
                    native_cache_builder=("local-bond-sparse" if local_cache_reported else "full-bond-scan"),
                    mnte=mnte, anchored_lookup=anchored,
                     fast_group_validation=fast_group,
                     anchored_verify=anchor_verified, sample_utc=now(), status="PASS",
                    all_input_sha256={f.name: sha(f) for f in work.iterdir()
                                      if f.suffix == ".dat" and f.name not in
                                      ("eigenvalues.dat", "diagnostics.dat")}))
    return obs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, help="Absolute prebuilt, tested ED2 binary")
    p.add_argument("--output", type=Path, required=True, help="Must NOT exist")
    p.add_argument("--source-revision", required=True, help="Pinned 40-digit Git commit")
    p.add_argument("--source-archive", type=Path, required=True)
    p.add_argument("--preset", choices=("cubic-smoke", "cubic-pilot", "cubic-main", "cubic-publication"), default="cubic-smoke")
    p.add_argument("--sectors", help="Comma-separated Gamma+O_h sectors; defaults all 5")
    p.add_argument("--sizes", help="Even side lengths; overrides preset")
    p.add_argument("--threads", help="OpenMP counts; overrides preset")
    p.add_argument("--cpu-list", help="Optional explicitly verified physical logical CPU IDs")
    p.add_argument("--warmups", type=int)
    p.add_argument("--repeats", type=int)
    p.add_argument("--mnte", type=int, default=-1,
                   help="Cubic fixed-three-magnon cache capacity: -1 original auto, >=54 conservative degree bound")
    p.add_argument("--timeout", type=int, default=1200)
    p.add_argument("--time-budget", type=int, default=36000)
    p.add_argument("--seed", type=int, default=1984)
    p.add_argument("--plan-only", action="store_true")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--require-pbs-exclusive", action="store_true")
    group.add_argument("--allow-shared", action="store_true")
    a = p.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", a.source_revision):
        p.error("--source-revision must be a full 40-digit Git SHA")
    defaults = {
        "cubic-smoke": dict(sizes="4", threads="1", warmups=0, repeats=1),
        "cubic-pilot": dict(sizes="4,6", threads="1,4", warmups=0, repeats=1),
        "cubic-main": dict(sizes="10", threads="1,4,8,16,32", warmups=1, repeats=3),
        "cubic-publication": dict(sizes="10", threads="1,8,32", warmups=2, repeats=7),
    }[a.preset]
    sizes = [int(x) for x in (a.sizes or defaults["sizes"]).split(",")]
    threads = [int(x) for x in (a.threads or defaults["threads"]).split(",")]
    sectors = tuple((a.sectors or ",".join(SECTORS)).split(","))
    if not sectors or len(sectors) != len(set(sectors)) or not set(sectors) <= set(SECTORS):
        p.error("Unknown or duplicate cubic Gamma/O_h sector")
    warmups = a.warmups if a.warmups is not None else defaults["warmups"]
    repeats = a.repeats if a.repeats is not None else defaults["repeats"]
    if not sizes or len(sizes) != len(set(sizes)) or any(n < 4 or n % 2 for n in sizes):
        p.error("Sizes must be unique even numbers >= 4")
    if not threads or len(threads) != len(set(threads)) or min(threads) < 1:
        p.error("Threads must be distinct positive integers")
    if not 0 <= warmups <= 5 or not 1 <= repeats <= 20:
        p.error("Invalid warmups/repeats")
    anchored = os.environ.get("ED2_CUBIC_ANCHORED_REP", "") == "1"
    anchored_verify = os.environ.get("ED2_CUBIC_ANCHORED_VERIFY", "") == "1"
    if anchored_verify and not anchored:
        p.error("ED2_CUBIC_ANCHORED_VERIFY requires ED2_CUBIC_ANCHORED_REP=1")
    fast_group = os.environ.get("ED2_CUBIC_FAST_GROUP_VALIDATION", "") == "1"
    anchored_basis = os.environ.get("ED2_CUBIC_ANCHORED_BASIS", "") == "1"
    local_bond_cache = os.environ.get("ED2_CUBIC_LOCAL_BOND_CACHE", "") == "1"
    profile_setting = os.environ.get("ED2_CUBIC_CACHE_PROFILE", "0")
    if profile_setting not in ("0", "1"):
        p.error("ED2_CUBIC_CACHE_PROFILE must be 0 or 1")
    cache_profile = profile_setting == "1"
    direct_setting = os.environ.get("ED2_CUBIC_DIRECT_HOP_TABLE", "0")
    if direct_setting not in ("0", "1"):
        p.error("ED2_CUBIC_DIRECT_HOP_TABLE must be 0 or 1")
    direct_hop_table = direct_setting == "1"
    if direct_hop_table and not (local_bond_cache and anchored and
                                 fast_group and anchored_basis):
        p.error("Direct cubic hop table requires all four previous cubic opt-ins")
    if cache_profile and not local_bond_cache:
        p.error("Cache profiling requires ED2_CUBIC_LOCAL_BOND_CACHE=1")
    if local_bond_cache and not (anchored and fast_group and anchored_basis):
        p.error("Local cubic cache requires anchored representative, fast group and anchored basis")
    if local_bond_cache and "gamma-t" in sectors:
        p.error("Local cubic cache requires a 1D O_h sector")
    if local_bond_cache and any(side not in (4, 6, 8, 10) for side in sizes):
        p.error("Local cubic cache requires cubic side 4, 6, 8 or 10")
    if anchored_basis and not fast_group:
        p.error("Anchored cubic basis requires ED2_CUBIC_FAST_GROUP_VALIDATION=1")
    if anchored_basis and "gamma-t" in sectors:
        p.error("Anchored cubic basis requires a 1D O_h sector")
    if anchored_basis and any(side not in (4, 6, 8, 10) for side in sizes):
        p.error("Anchored cubic basis currently requires cubic side 4, 6, 8 or 10")
    if fast_group and "gamma-t" in sectors:
        p.error("Fast cubic group certificate requires a noncommuting 1D O_h sector")
    if anchored and "gamma-t" in sectors:
        p.error("Anchored cubic representatives require a 1D O_h sector")
    if a.mnte != -1 and a.mnte < 54:
        p.error("MNTE must be -1 (original) or >=54 (3 magnons times 18 neighbors)")
    if a.timeout < 5 or a.time_budget < a.timeout:
        p.error("Invalid per-run timeout / total time budget")
    if a.time_budget > 40000:
        p.error("Cubic profile must fit within safe 12-hour PBS walltime")
    root = a.output.resolve()
    archive = a.source_archive.resolve(strict=True)
    if a.plan_only:
        count = sum(1 if sec == "gamma-t" or side > 4 else 2
                    for side in sizes for sec in sectors for _ in threads) * (warmups + repeats)
        print(json.dumps(dict(sizes=sizes, threads=threads, sectors=sectors, warmups=warmups, repeats=repeats,
                              trial_count=count, mnte=a.mnte,
                               anchored_lookup=anchored, anchored_verify=anchored_verify,
                               fast_group_validation=fast_group, anchored_basis=anchored_basis,
                   local_bond_cache=local_bond_cache,
                   cache_profile=cache_profile, direct_hop_table=direct_hop_table,
                              source_revision=a.source_revision, archive_sha256=sha(archive)),
                         indent=2))
        return
    exe = a.exe.resolve(strict=True) if a.exe else p.error("--exe is required for execution")
    if root.exists():
        p.error("Output directory already exists; refusal prevents measurement overwrite")
    root.mkdir(parents=True)
    (root / "trials").mkdir()
    started = time.monotonic()
    original_affinity = os.sched_getaffinity(0)
    record = dict(schema=1, created_utc=now(),
                  source_revision=a.source_revision,
                  archive_sha256=sha(archive), archive_path=str(archive),
                  binary_path=str(exe), binary_sha256=sha(exe),
                  bench_sha256=sha(__file__),
                  profile_oracle_sha256=sha(physics.__file__),
                  cubic_oracle_sha256=sha(cubic.__file__), geometry="3D_simple_cubic",
                  symmetry_scope="Gamma translations and four 1D O_h irreps",
                  legacy_rule="force old basis only for side=4, no infeasible large-side legacy",
                  preset=a.preset, sizes=sizes, threads=threads, sectors=sectors, warmups=warmups,
                  repeats=repeats, mnte=a.mnte, anchored_lookup=anchored,
                   anchored_verify=anchored_verify,
                   fast_group_validation=fast_group, anchored_basis=anchored_basis,
                   local_bond_cache=local_bond_cache,
                   cache_profile=cache_profile, direct_hop_table=direct_hop_table,
                   timeout=a.timeout, time_budget=a.time_budget,
                  seed=a.seed, protocol="preset-specific measured repeats and warmups",
                  publication_ready=False,
                  reason="Operator confirmation of exclusivity and independent review required")
    try:
        allocation_info = allocation(a.allow_shared, threads, a.cpu_list)
        record["environment"] = allocation_info
        write_json(root / "environment.json", record)
        os.environ["OMP_DYNAMIC"] = "FALSE"
        os.environ["OMP_PLACES"] = "cores"
        os.environ["OMP_PROC_BIND"] = "close"
        os.environ["ED2_RNG_SEED"] = str(a.seed)
        os.environ["ED2_PHASE_TIMING"] = "1"
        os.environ["ED2_RESIDUAL_ATOL"] = "1e-8"
        os.environ["ED2_ALLOW_LARGE_GROUP_VALIDATION"] = "1"
        for var in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "GOTO_NUM_THREADS"):
            os.environ[var] = "1"
        os.environ.pop("GOMP_CPU_AFFINITY", None)
        jobs = []
        for side in sizes:
            for sector in sectors:
                for nt in threads:
                    jobs.append((side, sector, nt))
        journal = root / "trials.jsonl"
        all_records = []
        for period, count in (("warmup", warmups), ("measure", repeats)):
            for repeat in range(count):
                shuffled = list(jobs)
                random.Random(a.seed + repeat + (0 if period == "warmup" else 100000)).shuffle(shuffled)
                for side, sector, nt in shuffled:
                    variants = (["fast", "forced-legacy"] if sector != "gamma-t" and side == 4
                                else ["fast"])
                    if repeat % 2 == 1:
                        variants.reverse()
                    pair = []
                    for variant in variants:
                        if time.monotonic() - started + a.timeout + 60 > a.time_budget:
                            raise TimeoutError("Time budget reached before next trial; results incomplete")
                        spec = (side, sector, nt, repeat, period, variant)
                        try:
                            obs = execute(exe, root, spec, allocation_info, a.timeout, a.mnte)
                            obs["source_revision"] = a.source_revision
                            obs["archive_sha256"] = record["archive_sha256"]
                            append_trial(journal, obs)
                            all_records.append(obs)
                            pair.append(obs)
                            print(f"PASS {obs['key']} D={obs['dimension']} "
                                  f"basis={obs['phase_wall_s']['basis']:.6g}s "
                                  f"total={obs['phase_wall_s']['total']:.6g}s",
                                  flush=True)
                        except Exception as exc:
                            err = dict(status="FAIL", spec=spec, error=repr(exc), utc=now())
                            append_trial(journal, err)
                            raise
                    if len(pair) == 2:
                        x, y = pair
                        if x["dimension"] != y["dimension"] or abs(x["energy"]-y["energy"]) > 1e-8:
                            append_trial(journal, dict(status="FAIL", error="Matched Ritz energy/dimension disagreement",
                                                       first=x["key"], second=y["key"], utc=now()))
                            raise AssertionError("Same-sector fast/legacy physical mismatch")
                        if x["all_input_sha256"] != y["all_input_sha256"]:
                            raise AssertionError("Same-sector generated input files are not identical")
        record["finished_utc"] = now()
        record["status"] = "complete"
        record["pass_count"] = len(all_records)
        write_json(root / "summary.json", record)
        print(f"SUCCESS: {len(all_records)} runs; independent publication attestation still pending")
    except Exception as exc:
        record["finished_utc"] = now()
        record["status"] = "incomplete"
        record["failure"] = repr(exc)
        write_json(root / "summary.json", record)
        print("INCOMPLETE:", repr(exc), file=sys.stderr)
        raise
    finally:
        os.sched_setaffinity(0, original_affinity)


if __name__ == "__main__":
    main()
