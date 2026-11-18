#!/usr/bin/env python3
"""Independent offline audit for PRESTO cubic O_h ED2 measurements.

Re-reads every physical log/input and all trials, rather than trusting the
producer's summary. Operator attestations are declarations, not automated
proof of exclusive use or independent editorial approval.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import cpc_cubic_pointgroup as cubic
import statistics
import sys

PHASES = ("basis", "operator_cache", "lanczos", "total")
SECTORS = cubic.SECTORS


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stat(data):
    m = statistics.median(data)
    return dict(n=len(data), median=m,
                mad=statistics.median(abs(v-m) for v in data),
                min=min(data), max=max(data))


def checked_pair(lines):
    assert sum(line.startswith("pair_columns=") for line in lines) == 1
    rows = []
    for line in lines:
        fields = line.split()
        if len(fields) == 8 and fields[0] == "1" and fields[1] in ("T", "F"):
            rows.append(fields)
    assert len(rows) == 1, "Missing native Ritz pair"
    fields = rows[0]
    assert fields[1:3] == ["T", "T"], "Uncertified Ritz pair"
    residual = float(fields[5])
    norm = float(fields[6])
    tol = float(fields[7])
    assert math.isfinite(residual) and residual <= min(1.01e-8, tol*1.01)
    assert math.isfinite(norm) and abs(norm - 1) < 1e-8
    return float(fields[3]), residual


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path, help="Must not exist")
    p.add_argument("--source-archive", type=Path)
    p.add_argument("--exe", type=Path)
    p.add_argument("--operator-attestation", type=Path,
                   help="External operator-signed JSON: job_id, host, verified_by, exclusive_confirmed")
    p.add_argument("--require-publication", action="store_true")
    a = p.parse_args()
    root = a.root.resolve(strict=True)
    output = a.output.resolve()
    if output.exists():
        p.error("Refusing to overwrite existing audit")
    env = load(root / "environment.json")
    summary = load(root / "summary.json")
    assert summary["status"] == "complete"
    assert summary["geometry"] == "3D_simple_cubic"
    assert isinstance(summary["mnte"], int)
    assert isinstance(summary["anchored_lookup"], bool)
    assert isinstance(summary["anchored_verify"], bool)
    assert isinstance(summary["fast_group_validation"], bool)
    assert isinstance(summary["anchored_basis"], bool)
    assert isinstance(summary["local_bond_cache"], bool)
    assert isinstance(summary["cache_profile"], bool)
    assert not summary["cache_profile"] or summary["local_bond_cache"]
    if summary["local_bond_cache"]:
        assert summary["anchored_lookup"] and summary["anchored_basis"] and summary["fast_group_validation"]
        assert "gamma-t" not in summary["sectors"]
    assert not summary["anchored_basis"] or summary["fast_group_validation"]
    if summary["anchored_basis"]:
        assert "gamma-t" not in summary["sectors"]
    if summary["fast_group_validation"]:
        assert "gamma-t" not in summary["sectors"]
    assert not summary["anchored_verify"] or summary["anchored_lookup"]
    if summary["anchored_lookup"]:
        assert "gamma-t" not in summary["sectors"]
    assert summary["mnte"] == -1 or summary["mnte"] >= 54
    assert summary["cubic_oracle_sha256"] == sha(cubic.__file__)
    assert env["source_revision"] == summary["source_revision"]
    assert env["binary_sha256"] == summary["binary_sha256"]
    assert re.fullmatch(r"[0-9a-f]{40}", summary["source_revision"])
    archive = a.source_archive or Path(summary["archive_path"])
    exe = a.exe or Path(summary["binary_path"])
    assert archive.is_file() and sha(archive) == summary["archive_sha256"]
    assert exe.is_file() and sha(exe) == summary["binary_sha256"]
    assert env["environment"]["selected_cpu_ids"]
    assert len(env["environment"]["selected_core_topology"]) >= max(summary["threads"])
    raw = [(i+1, json.loads(s)) for i, s in
           enumerate((root / "trials.jsonl").read_text().splitlines()) if s.strip()]
    assert raw and all(rec.get("status") == "PASS" for _, rec in raw), "Failed trial in journal"
    sectors = tuple(summary["sectors"])
    assert len(sectors) == len(set(sectors)) and set(sectors) <= set(SECTORS)
    for side in summary["sizes"]:
        assert cubic.check_hamiltonian(side)
        dims = cubic.sector_dimensions(side)
        assert all(dims[sec] > 0 for sec in sectors)
    expected = {(side, sector, threads, period, i, builder)
                for side in summary["sizes"] for sector in sectors
                for threads in summary["threads"]
                for period, count in (("warmup", summary["warmups"]), ("measure", summary["repeats"]))
                for i in range(count)
                for builder in (("fast", "forced-legacy") if sector != "gamma-t" and side == 4
                                else ("fast",))}
    seen = {}
    for index, record in raw:
        key = (record["side"], record["sector"], record["threads"],
               record["period"], record["repetition"], record["builder"])
        assert key not in seen, ("Duplicate trial", key)
        seen[key] = record
        assert record["source_revision"] == summary["source_revision"]
        assert record["archive_sha256"] == summary["archive_sha256"]
        assert record["binary_sha256"] == summary["binary_sha256"]
        assert record["mode"] == "two-stage"
        assert record["mnte"] == summary["mnte"]
        assert record["anchored_lookup"] is summary["anchored_lookup"]
        assert record["anchored_verify"] is summary["anchored_verify"]
        assert record["fast_group_validation"] is summary["fast_group_validation"]
        assert record["anchored_basis"] is summary["anchored_basis"]
        assert record["local_bond_cache"] is summary["local_bond_cache"]
        assert record["cache_profile"] is summary["cache_profile"]
        raw_input = (root / "trials" / record["key"] / "input.dat").read_text()
        match = re.search(r"\bMNTE\s*=\s*(-?\d+)", raw_input, flags=re.IGNORECASE)
        assert match and int(match.group(1)) == summary["mnte"]
        assert record["n_sites"] == record["side"] ** 3
        assert record["dimension"] == cubic.sector_dimensions(record["side"])[record["sector"]]
        assert record["residual"] <= 1.01e-8
        path = root / "trials" / record["key"]
        assert path.is_dir(), ("Missing trial directory", path)
        assert sha(path / "input.dat") == record["input_sha256"]
        for filename, checksum in record["all_input_sha256"].items():
            assert sha(path / filename) == checksum, (path, filename, "source input changed")
        assert sha(path / "diagnostics.txt") == record["diagnostic_sha256"]
        lines = (path / "diagnostics.txt").read_text().splitlines()
        diagnostic = dict(x.split("=", 1) for x in lines if "=" in x)
        assert diagnostic["lanczos_mode"].strip() == "two_stage"
        assert int(diagnostic["h_applications"]) == record["h_applications"]
        energy, residual = checked_pair(lines)
        assert abs(energy - record["energy"]) < 1e-8
        assert abs(residual - record["residual"]) < 1e-12
        values = [float(v.split()[-1].replace("D", "E")) for v in
                  (path / "eigenvalues.dat").read_text().splitlines() if v.strip()]
        assert len(values) == 1 and abs(values[0] - energy) < 1e-8
        log = (path / "stdout.log").read_text(errors="replace")
        assert ("Representative lookup: cubic-anchored-orbits" in log) is summary["anchored_lookup"]
        expected_certificate = summary["fast_group_validation"] and record["builder"] == "fast"
        assert ("Symmetry certificate: canonical cubic O_h semidirect product verified" in log) is expected_certificate
        if summary["anchored_verify"]:
            assert "per-hop legacy verification=T" in log
        atomic = "Basis builder: atomic-orbits" in log
        anchored_basis_used = "Basis builder: anchored-cubic" in log
        expected_anchored = record["builder"] == "fast" and summary["anchored_basis"]
        assert anchored_basis_used == expected_anchored
        assert atomic == record["atomic_orbits"]
        assert atomic == (record["builder"] == "fast" and not summary["anchored_basis"])
        expected_sparse_cache = summary["local_bond_cache"] and record["builder"] == "fast"
        assert ("Cubic cache builder: local-bond-sparse" in log) is expected_sparse_cache
        expected_cache = "local-bond-sparse" if expected_sparse_cache else "full-bond-scan"
        assert record["native_cache_builder"] == expected_cache
        expect_profile = summary["cache_profile"] and record["builder"] == "fast"
        assert ("Cubic cache profile: representative_key_setup_wall_s=" in log) is expect_profile
        assert ("Cubic cache profile: resolver_calls=" in log) is expect_profile
        expected_builder = ("anchored-cubic" if expected_anchored else
                            "atomic-orbits" if atomic else "forced-legacy")
        assert record["native_basis_builder"] == expected_builder
        if record["builder"] == "fast" and record["sector"] != "gamma-t":
            assert "Symmetry action: noncommuting normal form and character validated" in log
        phases = dict((phase, float(value)) for phase, value in
                      re.findall(r"ED2_PHASE ([a-z_]+) wall_s=\s*([+\-0-9.eE]+)", log))
        assert all(phase in phases and math.isfinite(phases[phase]) and phases[phase] >= 0
                   for phase in PHASES)
        for phase in PHASES:
            assert abs(phases[phase] - record["phase_wall_s"][phase]) < 5e-10, (path, phase)
    assert set(seen) == expected, ("Incomplete trial factorial",
                                   len(expected - set(seen)), len(set(seen) - expected))
    # Audit exact same input and physics per matched fast/legacy pair. Compare
    # actual Ritz-pair energies and all generated physical input checksums.
    for side in summary["sizes"]:
        if side != 4:
            continue
        for sector in sectors:
            if sector == "gamma-t":
                continue
            for threads in summary["threads"]:
                for period, count in (("warmup", summary["warmups"]),
                                      ("measure", summary["repeats"])):
                    for i in range(count):
                        left = seen[side, sector, threads, period, i, "fast"]
                        right = seen[side, sector, threads, period, i, "forced-legacy"]
                        assert left["all_input_sha256"] == right["all_input_sha256"]
                        assert left["dimension"] == right["dimension"]
                        assert abs(left["energy"]-right["energy"]) <= 1e-8
    stats = {}
    for side in summary["sizes"]:
        for sector in sectors:
            for threads in summary["threads"]:
                records = {builder: [seen[side, sector, threads, "measure", i, builder]
                                    for i in range(summary["repeats"])]
                           for builder in (("fast", "forced-legacy") if sector != "gamma-t" and side == 4
                                           else ("fast",))}
                row = {}
                for builder, trials in records.items():
                    row[builder] = {
                        "timings_s": {phase: stat([rec["phase_wall_s"][phase]
                                                   for rec in trials]) for phase in PHASES},
                        "rss_kib": stat([rec["peak_rss_kib"] for rec in trials
                                         if rec["peak_rss_kib"] is not None]) if all(
                            rec["peak_rss_kib"] is not None for rec in trials) else None,
                        "h_applications": sorted({rec["h_applications"] for rec in trials}),
                        "dimension": trials[0]["dimension"],
                        "energy": trials[0]["energy"],
                    }
                if "forced-legacy" in row:
                    row["same_sector_legacy_over_fast_ratio"] = {
                        phase: row["forced-legacy"]["timings_s"][phase]["median"] /
                               max(1e-14, row["fast"]["timings_s"][phase]["median"])
                        for phase in PHASES}
                stats[f"s{side}-{sector}-t{threads}"] = row
    attested = False
    attestation = None
    if a.operator_attestation:
        attestation = load(a.operator_attestation)
        assert attestation["job_id"] == env["environment"]["pbs_jobid"]
        assert attestation["host"].split(".")[0] == env["environment"]["host"].split(".")[0]
        assert attestation["verified_by"].strip() and attestation["exclusive_confirmed"] is True
        assert attestation.get("evidence_reviewed", "").strip()
        attested = True
    eligible = (env["environment"]["mode"] == "pbs_exclhost_requested_attestation_pending"
                and attested and summary["warmups"] >= 2 and summary["repeats"] >= 7
                 and not summary["anchored_verify"])
    if a.require_publication and not eligible:
        raise RuntimeError("Publication prerequisite gate fails: operator evidence, PBS or replicates insufficient")
    report = dict(schema=1, audited=True, raw_trials=len(raw),
                  warmups=summary["warmups"], repeats=summary["repeats"],
                  source_revision=summary["source_revision"],
                  binary_sha256=summary["binary_sha256"],
                  archive_sha256=summary["archive_sha256"],
                  job_id=env["environment"]["pbs_jobid"],
                  host=env["environment"]["host"],
                  pbs_mode=env["environment"]["mode"],
                  geometry=summary["geometry"], sectors=sectors, mnte=summary["mnte"],
                  anchored_lookup=summary["anchored_lookup"],
                  anchored_verify=summary["anchored_verify"],
                  fast_group_validation=summary["fast_group_validation"],
                  anchored_basis=summary["anchored_basis"],
                  local_bond_cache=summary["local_bond_cache"],
                  cache_profile=summary["cache_profile"],
                  operator_attestation_present=attested,
                  review_eligible=eligible,
                  caution="Operator claim is not independent scheduler proof, and review gate is not publication approval.",
                  attestation=attestation, statistics=stats)
    output.mkdir(parents=True)
    (output / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(f"AUDIT PASS: {len(raw)} distinct trials, {len(stats)} configurations; "
          f"operator review eligible={eligible}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
