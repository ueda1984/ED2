#!/usr/bin/env python3
"""Reproducible, provenance-rich ED2 translation vs point-group performance profile.

This is *not* an apples-to-apples eigensolver speedup claim: symmetry
sectors have different dimensions. No numeric timing threshold is a test gate.
Only the Python standard library is required (including on PRESTO).
"""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


def binomial(n, k):
    """Exact C(n,k) without math.comb, available only in Python >=3.8."""
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    result = 1
    for i in range(1, k + 1):
        result = result * (n - i + 1) // i
    return result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_bonds(side):
    def p(x, y):
        return (x % side) + side * (y % side)
    bonds = []
    for y in range(side):
        for x in range(side):
            for dx, dy, xy, zz in ((1, 0, 1., .77), (0, 1, 1., .77),
                                   (1, 1, .31, .19), (1, -1, .31, .19)):
                bonds.append((p(x, y), p(x + dx, y + dy), xy, zz))
    assert len(set((min(a, b), max(a, b), xy, zz) for a, b, xy, zz in bonds)) == len(bonds)
    return bonds


def symmetry_generators(side, symmetry):
    n = side * side
    def permutation(f):
        return tuple(f(i % side, i // side) for i in range(n))
    tx = permutation(lambda x, y: ((x + 1) % side) + side * y)
    ty = permutation(lambda x, y: x + side * ((y + 1) % side))
    r = permutation(lambda x, y: ((-y) % side) + side * x)
    mirror = permutation(lambda x, y: x + side * ((-y) % side))
    c2 = permutation(lambda x, y: ((-x) % side) + side * ((-y) % side))
    if symmetry == "gamma-t":
        return [(tx, side, 0, "Tx"), (ty, side, 0, "Ty")]
    if symmetry == "gamma-a1":
        return [(tx, side, 0, "Tx"), (ty, side, 0, "Ty"),
                (r, 4, 0, "C4"), (mirror, 2, 0, "mirror")]
    if symmetry == "gamma-b1":
        return [(tx, side, 0, "Tx"), (ty, side, 0, "Ty"),
                (r, 4, 2, "C4"), (mirror, 2, 0, "mirror")]
    if symmetry == "x-a1":
        assert side % 2 == 0
        return [(tx, side, side // 2, "Tx"), (ty, side, 0, "Ty"),
                (c2, 2, 0, "C2"), (mirror, 2, 0, "mirror")]
    raise ValueError("Invalid symmetry")


def check_bonds(bonds, gens):
    edges = sorted((min(a, b), max(a, b), jxy, jz) for a, b, jxy, jz in bonds)
    for p, period, m, name in gens:
        assert sorted((min(p[a], p[b]), max(p[a], p[b]), jxy, jz)
                      for a, b, jxy, jz in bonds) == edges, ("Noninvariant Hamiltonian", name)
        assert sorted(p) == list(range(len(p))), ("Nonpermutation", name)


def independent_dimension(nsites, down, gens):
    """Trace-character projector with fixed-magnon Burnside counting."""
    order = math.prod(g[1] for g in gens)
    seen = set()
    trace_sum = 0
    for exps in itertools.product(*(range(g[1]) for g in gens)):
        perm = list(range(nsites))
        for (g, _, _, _), exp in zip(reversed(gens), reversed(exps)):
            for _ in range(exp):
                perm = [g[i] for i in perm]
        perm = tuple(perm)
        assert perm not in seen, "Generator normal form contains duplicate group elements"
        seen.add(perm)
        phase = sum(m * e / period for (_, period, m, _), e in zip(gens, exps))
        chi = math.cos(2 * math.pi * phase)
        assert abs(abs(chi) - 1) < 1e-8, "Not a real one-dimensional character"
        character = round(chi)
        visited = [False] * nsites
        poly = [1] + [0] * down
        for i in range(nsites):
            if visited[i]:
                continue
            cur = i
            size = 0
            while not visited[cur]:
                visited[cur] = True
                size += 1
                cur = perm[cur]
            # Only occupations 0 or 1 allowed at each spin-1/2 site.
            for j in range(down, size - 1, -1):
                poly[j] += poly[j - size]
        trace_sum += character * poly[down]
    assert len(seen) == order and trace_sum % order == 0, "Nonintegral character projector"
    dimension = trace_sum // order
    assert 0 < dimension <= binomial(nsites, down), "Invalid physical symmetry dimension"
    return dimension


def write_input(work, side, down, gens, bonds, maxiter, wk_dim):
    work.mkdir(parents=True, exist_ok=False)
    n = side * side
    (work / "spins.dat").write_text("0.5\n" * n)
    (work / "caps.dat").write_text("1\n" * n)
    (work / "fields.dat").write_text("")
    (work / "bonds.dat").write_text("".join(
        f"{i+1} {j+1} {xy:g} {xy:g} {zz:g} 0 0 0 0 0 0\n"
        for i, j, xy, zz in bonds))
    lines = [
        f"NOS={n},NODmin={down},NODmax={down},ALG=1,",
        f"NO_one=0,NO_two={len(bonds)},wk_dim={wk_dim},MNTE=-1,",
        "FILE_xyz_dm_gamma='bonds.dat',FILE_hvec='fields.dat',",
        "FILE_spin='spins.dat',FILE_NODmax='caps.dat',",
        "OUTDIR='./',FILE_wf='./',wr_wf=0,re_wf=0,cal_lm=0,cal_cf=0,"]
    for k in range(6):
        if k < len(gens):
            g, order, char, label = gens[k]
            f = f"symmetry{k+1}.dat"
            (work / f).write_text("".join(f"{index+1}\n" for index in g))
            lines.append(f"L{k+1}={order},M{k+1}={char},FILE_S{k+1}='{f}',")
        else:
            lines.append(f"L{k+1}=1,M{k+1}=0,")
    (work / "input.dat").write_text(
        "&input_parameters\n" + "\n".join(lines) + "\n/\n"
        f"&input_lancz\nlnc_ene_conv=1d-11,minitr=5,maxitr={maxiter},itrint=5\n/\n"
        f"&input_TRLan\nNOE=2,NOK=3,NOM=8,maxitr={maxiter},"
        "lnc_ene_conv=1d-11,i_vec_min=1,i_vec_max=1\n/\n")


def run_case(exe, work, nsites, down, symmetry, dimension, threads, mode, timeout):
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", OMP_PROC_BIND="close",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               ED2_RNG_SEED="1984", ED2_PHASE_TIMING="1",
               ED2_DIAGNOSTICS="diagnostics.txt", ED2_RESIDUAL_ATOL="1e-8")
    env.pop("ED2_LANCZOS_TWO_STAGE", None)
    if mode == "two-stage":
        env["ED2_LANCZOS_TWO_STAGE"] = "1"
    command = [str(exe)]
    util = shutil.which("/usr/bin/time")
    if util:
        command = [util, "-f", "%e %M %U %S", "-o", str(work / "resource.txt"), str(exe)]
    start = time.perf_counter()
    with (work / "input.dat").open("rb") as infile, (work / "stdout.log").open("wb") as log:
        cp = subprocess.run(command, cwd=work, stdin=infile, stdout=log,
                            stderr=subprocess.STDOUT, env=env, timeout=timeout, check=False)
    elapsed = time.perf_counter() - start
    output = (work / "stdout.log").read_text(errors="replace")
    assert cp.returncode == 0, (work, cp.returncode, output[-3000:])
    found = re.findall(r"THS\(k\)\s*=\s*(\d+)", output)
    assert len(found) == 1 and int(found[0]) == dimension, (symmetry, found, dimension)
    energies = [float(line.split()[-1].replace("D", "E"))
                for line in (work / "eigenvalues.dat").read_text().splitlines() if line.strip()]
    assert len(energies) == 1 and math.isfinite(energies[0])
    diag_lines = (work / "diagnostics.txt").read_text().splitlines()
    diag = dict(line.split("=", 1) for line in diag_lines if "=" in line)
    assert diag.get("lanczos_mode", "").strip() == ("two_stage" if mode == "two-stage" else "legacy")
    # ED2's schema=1 writes a *tabular* Ritz-pair row after pair_columns=... .
    # The old 'pair_1_residual' dictionary key does not exist: never skip
    # checking actual residuals because that optional key is absent.
    assert sum(line.startswith("pair_columns=") for line in diag_lines) == 1
    pairs = [line.split() for line in diag_lines if re.match(r"^\s*1\s+[TF]\s+[TF]\s+", line)]
    assert len(pairs) == 1 and len(pairs[0]) == 8, ("Missing native Ritz-pair diagnostic", work)
    pair = pairs[0]
    assert pair[1] == "T" and pair[2] == "T", ("Uncertified eigenpair", work, pair)
    assert abs(float(pair[3]) - energies[0]) < 1e-8
    assert abs(float(pair[6]) - 1.) < 1e-8
    residual = float(pair[5])
    threshold = float(pair[7])
    assert math.isfinite(residual) and residual <= min(1.01e-8, threshold * 1.01), (symmetry, residual, threshold)
    h_calls = int(diag["h_applications"])
    assert h_calls > 0
    phases = {}
    for p, seconds in re.findall(r"ED2_PHASE ([a-z_]+) wall_s=\s*([+\-0-9.eE]+)", output):
        assert p not in phases, ("Duplicate timing phase", p)
        phases[p] = float(seconds)
    required = {"input", "combinatorics", "basis", "operator_cache", "lanczos",
                "observables_and_output", "total"}
    assert required == set(phases), (symmetry, required - set(phases), set(phases) - required)
    assert all(math.isfinite(v) and v >= 0 for v in phases.values())
    assert sum(v for k, v in phases.items() if k != "total") <= phases["total"] + 1e-3
    rss = None
    cpu_s = None
    if (work / "resource.txt").exists():
        r = (work / "resource.txt").read_text().strip().split()
        if len(r) == 4:
            rss, cpu_s = int(r[1]), float(r[2]) + float(r[3])
    return {
        "n_sites": nsites, "down_spins": down, "symmetry": symmetry,
        "dimension": dimension, "mode": mode, "threads": threads,
        "energy": energies[0], "residual": residual,
        "h_applications": h_calls,
        "reconstructions": int(diag.get("vector_reconstructions", 0)),
        "phase_wall_s": phases, "process_elapsed_s": elapsed,
        "peak_rss_kib": rss, "cpu_s": cpu_s,
        "input_sha256": digest(work / "input.dat"),
        "diagnostic_sha256": digest(work / "diagnostics.txt"),
        "binary_sha256": digest(exe)
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tier", choices=("smoke", "stress"), default="smoke")
    parser.add_argument("--sizes", default=None, help="Comma-separated square side lengths")
    parser.add_argument("--threads", default="1,4", help="Comma-separated OpenMP thread counts")
    parser.add_argument("--modes", default="legacy,two-stage")
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--wk-dim", type=int, default=-1)
    args = parser.parse_args()
    exe = args.exe.resolve(strict=True)
    sizes = [int(x) for x in (args.sizes or ("4,6" if args.tier == "smoke" else "8,10")).split(",")]
    threads = [int(x) for x in args.threads.split(",")]
    modes = args.modes.split(",")
    assert sizes and all(x >= 4 and x % 2 == 0 for x in sizes)
    assert threads and all(t >= 1 for t in threads) and len(threads) == len(set(threads))
    assert modes and set(modes).issubset({"legacy", "two-stage"})
    assert args.repeat >= 1 and args.timeout >= 5
    root = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix="ed2-pg-benchmark-"))
    if args.output:
        root.mkdir(parents=True, exist_ok=False)
    results = []
    meta = {}
    for side in sizes:
        bonds = canonical_bonds(side)
        for symmetry in ("gamma-t", "gamma-a1", "gamma-b1", "x-a1"):
            gens = symmetry_generators(side, symmetry)
            check_bonds(bonds, gens)
            dim = independent_dimension(side * side, 3, gens)
            key = f"square{side}-{symmetry}"
            meta[key] = dict(n_sites=side * side, sector=symmetry, dimension=dim,
                             group_order=math.prod(g[1] for g in gens), magnons=3)
            for t in threads:
                for repetition in range(args.repeat):
                    pair = []
                    for mode in modes:
                        work = root / f"{key}-t{t}-r{repetition}-{mode}"
                        write_input(work, side, 3, gens, bonds, 2000, args.wk_dim)
                        r = run_case(exe, work, side * side, 3, symmetry, dim,
                                     t, mode, args.timeout)
                        r["repetition"] = repetition
                        results.append(r)
                        pair.append(r)
                        print(f"PASS {key} OMP={t} repeat={repetition} {mode} D={dim} "
                              f"basis_s={r['phase_wall_s']['basis']:.6g} "
                              f"H={r['h_applications']} E={r['energy']:.12g}", flush=True)
                    if len(pair) == 2:
                        assert pair[0]["input_sha256"] == pair[1]["input_sha256"]
                        assert abs(pair[0]["energy"] - pair[1]["energy"]) < 1e-8, (key, pair)
    report = dict(schema=1, passed=True, tier=args.tier, sizes=sizes, threads=threads,
                  modes=modes, repeat=args.repeat, wk_dim=args.wk_dim,
                  binary_sha256=digest(exe), fixtures=meta, results=results,
                  warning="CI/shared-node timings are diagnostic only, not publication benchmarks. "
                          "Translation and point-group runs have different symmetry dimensions; "
                          "phase and full-program times cannot be interpreted as isolated solver speedups.",
                  opt_in="ED2_PHASE_TIMING=1; no profiler output or source physics behavior changes by default")
    (root / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"POINTGROUP BENCHMARK PASS {len(results)} native runs; output={root}", flush=True)


if __name__ == "__main__":
    main()
