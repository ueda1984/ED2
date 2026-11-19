#!/usr/bin/env python3
"""Read-only PRESTO staging and latest-main feature gate (Python >=3.6).

The scientific A1g comparison requires ED2's *noncommuting* cubic group
validation AND direct Ritz-pair diagnostics. An older generic main can be
archived but must NOT be presented as an equivalent A1g timing comparison.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile


NECESSARY = (
    "ED2_CUBIC_FAST_GROUP_VALIDATION",
    "ED2_CUBIC_ANCHORED_REP",
    "ED2_CUBIC_ANCHORED_BASIS",
    "ED2_CUBIC_LOCAL_BOND_CACHE",
    "ED2_CUBIC_DIRECT_HOP_TABLE",
    "ED2_LANCZOS_TWO_STAGE",
    "ED2_PHASE_TIMING",
    "ED2_DIAGNOSTICS",
)
ARCHIVES = ("ed2-main.tar.gz", "benchmark-tools.tar.gz", "launch-tools.tar.gz")
REQUIRED_BENCH = (
    "benchmarks/cross_library/run_quspin.py",
    "benchmarks/cross_library/export_cubic_xxz.py",
    "benchmarks/cross_library/cubic_oh_generators.py",
)
REQUIRED_LAUNCH = (
    "tools/presto_cross_library/launch.pbs",
    "tools/presto_cross_library/runner.py",
    "tools/presto_cross_library/preflight.py",
)


def digest(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_archive(stage, filename):
    path = stage / filename
    if not path.is_file():
        raise RuntimeError("Missing immutable source archive " + filename)
    with tarfile.open(str(path), "r:gz") as tf:
        names = set(m.name.lstrip("./") for m in tf.getmembers())
        source = {}
        for member in tf.getmembers():
            if member.isfile() and member.name.startswith("source/") and (
                    member.name.endswith(".f90") or member.name.endswith(".F90")):
                f = tf.extractfile(member)
                source[member.name] = f.read().decode("utf8", "replace")
    return names, source


def check_stage(stage):
    manifest = stage / "SHA256SUMS"
    if not manifest.is_file():
        raise RuntimeError("No SHA256SUMS for immutable PRESTO stage")
    actual = {}
    for record in manifest.read_text().splitlines():
        m = re.fullmatch(r"([0-9a-f]{64})  ([A-Za-z0-9._-]+)", record.strip())
        if m is None:
            raise RuntimeError("Malformed SHA256SUMS")
        name = m.group(2)
        if name in actual or name not in ARCHIVES:
            raise RuntimeError("Duplicate/unexpected archive " + name)
        have = digest(stage / name)
        if have != m.group(1):
            raise RuntimeError("Archive hash mismatch for " + name)
        actual[name] = have
    if set(actual) != set(ARCHIVES):
        raise RuntimeError("Incomplete archive manifest")
    commits = {}
    for kind in ("main", "benchmark", "launcher"):
        sha = (stage / (kind + "-commit.txt")).read_text().strip()
        if not re.fullmatch("[0-9a-f]{40}", sha):
            raise RuntimeError("Bad pinned " + kind + " commit")
        commits[kind] = sha
    _, src = read_archive(stage, ARCHIVES[0])
    if not src:
        raise RuntimeError("ED2 source archive contains no Fortran source")
    bench_names, _ = read_archive(stage, ARCHIVES[1])
    launch_names, _ = read_archive(stage, ARCHIVES[2])
    for target in REQUIRED_BENCH:
        if target not in bench_names:
            raise RuntimeError("Benchmark snapshot missing " + target)
    for target in REQUIRED_LAUNCH:
        if target not in launch_names:
            raise RuntimeError("Launch snapshot missing " + target)
    joined = "\n".join(src.values())
    capabilities = {x: x in joined for x in NECESSARY}
    # The specific fast cubic validation marker implies a compiled
    # noncommuting space-group certificate. Never infer from support of six
    # commuting cyclic generators alone.
    ready = all(capabilities.values())
    return dict(status="READY" if ready else "SOURCE_NOT_READY",
                commits=commits, archive_sha256=actual,
                capabilities=capabilities,
                source_file_count=len(src),
                meaning="The full N1000 A1g ED2 baseline and specialized "
                        "comparison requires these compiled features; "
                        "staging does not imply execution or correctness.")


def latest_main(repo):
    cp = subprocess.run(["git", "-C", str(repo), "ls-remote", "origin",
                         "refs/heads/main"], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE, universal_newlines=True,
                        timeout=30)
    if cp.returncode != 0 or not cp.stdout.strip():
        raise RuntimeError("Cannot verify remote main HEAD before submission")
    sha = cp.stdout.split()[0]
    if not re.fullmatch("[0-9a-f]{40}", sha):
        raise RuntimeError("Invalid remote main SHA")
    return sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, type=Path)
    ap.add_argument("--mode", choices=("stage", "submit"), default="stage")
    ap.add_argument("--repo", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    report = check_stage(args.stage)
    report["mode"] = args.mode
    if args.mode == "submit":
        if not args.repo:
            ap.error("--repo required to prove current main before submission")
        fresh = latest_main(args.repo)
        report["current_remote_main_sha"] = fresh
        if fresh != report["commits"]["main"]:
            report["status"] = "STALE_MAIN_RESTAGE_REQUIRED"
        elif report["status"] != "READY":
            report["status"] = "BLOCKED_MISSING_MAIN_CAPABILITIES"
    blob = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        if args.output.exists():
            raise RuntimeError("Refusing overwrite of existing audit report")
        args.output.write_text(blob)
    sys.stdout.write(blob)
    if args.mode == "submit" and report["status"] != "READY":
        sys.stderr.write("REFUSED: do not submit an unready/stale ED2 benchmark\n")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
