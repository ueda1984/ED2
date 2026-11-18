#!/usr/bin/env python3
"""Reproducible shared-host pilot: alternating mixed-spin real H, 3 H-action modes.

This is a *scientific operator equivalence and phase-observation harness*,
not a journal-grade performance benchmark. Real mixed S=1/2, S=1 chain,
genuinely non-U(1) (Jxx != Jyy plus transverse hx), OBC, no DM.
Only stdlib Python >=3.6 is required; ED2 needs normal Fortran BLAS/LAPACK.

Three otherwise identical executions per case:
 complex: ED2_REAL_CACHED_MATVEC=0, ED2_PACKED_REAL_CACHE=0
 real:    ED2_REAL_CACHED_MATVEC=1, ED2_PACKED_REAL_CACHE=0
 packed:  ED2_REAL_CACHED_MATVEC=1, ED2_PACKED_REAL_CACHE=1
"""
import argparse
import datetime
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import random
import re
import statistics
import subprocess


MODES = {
    "complex": ("0", "0"),
    "real": ("1", "0"),
    "packed": ("1", "1"),
}
FILES = ("input.dat", "spins.dat", "caps.dat", "bonds.dat", "fields.dat")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for buf in iter(lambda: handle.read(1048576), b""):
            h.update(buf)
    return h.hexdigest()


def median_mad(values):
    mid = statistics.median(values)
    return {
        "n": len(values),
        "median_ms": mid,
        "mad_ms": statistics.median(abs(v - mid) for v in values),
        "minimum_ms": min(values),
        "maximum_ms": max(values),
    }


def make_case(dest, n):
    if n not in (4, 6, 8, 10):
        raise ValueError("Only 4,6,8,10-site mixed-spin pilot cases are allowed")
    dest.mkdir(parents=True, exist_ok=False)
    spins = [0.5, 1.0] * (n // 2)
    (dest / "spins.dat").write_text("".join("{:.1f}\n".format(v) for v in spins))
    (dest / "caps.dat").write_text("".join("{}\n".format(int(2 * v)) for v in spins))
    bonds = []
    for j in range(n - 1):
        sign = (-1) ** j
        # Jxx != Jyy breaks U(1); Jxy/DM and y fields remain zero.
        val = [0.77 + 0.05 * sign, 0.63 + 0.03 * sign,
               1.09 - 0.04 * sign, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        bonds.append("{} {} {}\n".format(
            j + 1, j + 2, " ".join("{:.17g}".format(x) for x in val)))
    (dest / "bonds.dat").write_text("".join(bonds))
    fields = []
    for j in range(n):
        sign = (-1) ** j
        val = [0.15 * (0.4 * sign + 0.1), 0.0, 0.17 + 0.035 * sign]
        fields.append("{} {}\n".format(
            j + 1, " ".join("{:.17g}".format(x) for x in val)))
    (dest / "fields.dat").write_text("".join(fields))
    cap = sum(int(2 * x) for x in spins)
    deck = """&input_parameters
NOS={n},NODmin=0,NODmax={cap},ALG=1,
L1=1,L2=1,L3=1,L4=1,L5=1,L6=1,
M1=0,M2=0,M3=0,M4=0,M5=0,M6=0,
NO_one={n},NO_two={bonds},wk_dim=-1,MNTE=-1,
FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',
FILE_spin='spins.dat',FILE_NODmax='caps.dat',
OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0,
/
&input_lancz
lnc_ene_conv=1d-12,minitr=5,maxitr=320,itrint=10
/
&input_TRLan
NOE=2,NOK=4,NOM=24,maxitr=320,lnc_ene_conv=1d-12,
i_vec_min=1,i_vec_max=2
/
""".format(n=n, cap=cap, bonds=n - 1)
    (dest / "input.dat").write_text(deck)


def run_case(exe, work, n, threads, repetition, period, mode, timeout, seed):
    make_case(work, n)
    env = dict(os.environ)
    for key in tuple(env):
        if key.startswith("ED2_CUBIC_") or key == "ED2_FORCE_LEGACY_ORBITS":
            env.pop(key, None)
    real, packed = MODES[mode]
    env.update(
        ED2_REAL_CACHED_MATVEC=real,
        ED2_PACKED_REAL_CACHE=packed,
        ED2_LANCZOS_TWO_STAGE="1",
        ED2_RESIDUAL_ATOL="1e-8",
        ED2_RNG_SEED=str(seed),
        ED2_PHASE_TIMING="1",
        OMP_DYNAMIC="FALSE",
        OMP_NUM_THREADS=str(threads),
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        BLIS_NUM_THREADS="1")
    with (work / "input.dat").open("rb") as inp, (work / "stdout.log").open("wb") as out:
        proc = subprocess.run([str(exe)], stdin=inp, stdout=out, stderr=subprocess.STDOUT,
                              env=env, cwd=str(work), timeout=timeout, check=False)
    stdout = (work / "stdout.log").read_bytes().decode("utf8", errors="replace")
    assert proc.returncode == 0, ("native ED2 run failed", work, stdout[-2500:])
    expected_real = mode != "complex"
    expected_packed = mode == "packed"
    assert ("Hamiltonian action: strictly-real cached coefficients" in stdout) == expected_real
    assert ("Hamiltonian action: compact-real cache extra_bytes=" in stdout) == expected_packed
    lanczos_status = re.findall(r"Lanczos status: ([^\\n]+)", stdout)
    assert len(lanczos_status) == 1 and "residual certified" in lanczos_status[0], (
        "Solver did not certify its final native residual", work, stdout[-2500:])
    pattern = r"ED2_PHASE {} wall_s=\s*([0-9.Ee+\-]+)"
    phases = {}
    for phase in ("basis", "operator_cache", "lanczos", "total"):
        m = re.findall(pattern.format(phase), stdout)
        assert len(m) == 1, (work, phase)
        phases[phase] = float(m[0]) * 1000
        assert math.isfinite(phases[phase]) and phases[phase] >= 0

    steps = re.findall(r"Two-stage replay count=\d+, H-applications=(\d+)", stdout)
    assert len(steps) == 1
    residuals = re.findall(r"Eigenpair residual norm:\s*([0-9.Ee+\-]+)", stdout)
    assert len(residuals) == 1 and float(residuals[0]) <= 1.01e-8
    return {
        "n": n, "threads": threads, "repetition": repetition,
        "period": period, "mode": mode,
        "input_sha256": {name: digest(work / name) for name in FILES},
        "eigen_sha256": digest(work / "eigenvalues.dat"),
        "h_applications": int(steps[0]), "residual": float(residuals[0]),
        "phase_ms": phases, "status": "PASS",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True, help="Must not exist")
    p.add_argument("--sizes", default="8,10")
    p.add_argument("--threads", default="1,4")
    p.add_argument("--warmups", type=int, default=1)
    p.add_argument("--repeats", type=int, default=5)
    p.add_argument("--timeout", type=int, default=120)
    p.add_argument("--seed", type=int, default=1984)
    a = p.parse_args()
    exe = a.exe.resolve(strict=True)
    sizes = [int(s) for s in a.sizes.split(",")]
    threads = [int(s) for s in a.threads.split(",")]
    assert sizes and all(s in (4, 6, 8, 10) for s in sizes)
    assert threads and all(1 <= t <= 8 for t in threads)
    assert a.warmups >= 0 and 1 <= a.repeats <= 20
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)

    records = []
    index = 0
    for n, t in itertools.product(sizes, threads):
        for period, count in (("warmup", a.warmups), ("measure", a.repeats)):
            for rep in range(count):
                modes = list(MODES)
                random.Random(a.seed + n * 100 + t * 10 + rep).shuffle(modes)
                group = []
                for mode in modes:
                    work = root / "n{}-t{}-{}{}-{}".format(n, t, period, rep, mode)
                    obs = run_case(exe, work, n, t, rep, period, mode,
                                   a.timeout, a.seed)
                    records.append(obs)
                    group.append(obs)
                    with (root / "trials.jsonl").open("a") as out:
                        out.write(json.dumps(obs, sort_keys=True) + "\n")
                    index += 1
                    print("GENERIC_REAL_CACHE_NATIVE_PASS n={} t={} {}{} {} {}/{}"
                          .format(n, t, period, rep, mode, index,
                                  len(sizes) * len(threads) *
                                  (a.warmups + a.repeats) * 3),
                          flush=True)
                assert len({json.dumps(o["input_sha256"], sort_keys=True)
                            for o in group}) == 1
                assert len({o["eigen_sha256"] for o in group}) == 1
                assert len({o["h_applications"] for o in group}) == 1
    summary = {}
    for n, t in itertools.product(sizes, threads):
        per_mode = {}
        for mode in MODES:
            r = [v for v in records if v["n"] == n and v["threads"] == t
                 and v["period"] == "measure" and v["mode"] == mode]
            assert len(r) == a.repeats
            per_mode[mode] = {
                phase: median_mad([v["phase_ms"][phase] for v in r])
                for phase in ("operator_cache", "lanczos", "total")
            }
        summary["n{}-t{}".format(n, t)] = per_mode
    report = {
        "status": "complete_shared_unattested_pilot",
        "source_type": "native ED2 executable",
        "binary_sha256": digest(exe),
        "source_sha256_or_git_revision": "record separately; not inferred from binary",
        "model": "open mixed S=1/2,S=1 XYZ real non-U1 chain",
        "sizes": sizes, "threads": threads,
        "warmups": a.warmups, "repeats": a.repeats,
        "operator_attested": False,
        "results": summary,
    }
    (root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print("GENERIC_REAL_CACHE_PILOT_PASS matched_trials={} output_sha_exact=True"
          .format(len(records)), flush=True)
    print("Shared-host exploratory only; operator exclusivity NOT attested.", flush=True)


if __name__ == "__main__":
    main()
