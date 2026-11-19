#!/usr/bin/env python3
"""Offline paired PRESTO auto-MNTE vs explicit-MNTE64 scientific comparison.

Only accept independently PASS-audited producer output from the SAME PBS
allocation and pinned native binary. Compare exact generated physical inputs
(except the declared MNTE token), Ritz E/residual, Lanczos H applications and
paired runtime/RSS metrics. Does not confer independent scheduler attestation.
No NumPy; Python >=3.6.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def norm_input(raw):
    assert len(re.findall(rb"\bMNTE\s*=\s*(?:-1|64)\b", raw)) == 1
    return re.sub(rb"\bMNTE\s*=\s*(?:-1|64)\b", b"MNTE=PAIRED", raw)


def records(root):
    out = {}
    for line in (root / "trials.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        assert rec["status"] == "PASS", ("Incomplete trial", rec)
        key = rec["key"]
        assert key not in out, ("Duplicate trial", key)
        out[key] = rec
    assert out
    return out


def medmad(values):
    median = statistics.median(values)
    return {"n":len(values), "median":median,
            "mad":statistics.median(abs(x-median) for x in values),
            "min":min(values), "max":max(values)}


def compare(auto_dir, explicit_dir):
    auto_dir, explicit_dir = Path(auto_dir), Path(explicit_dir)
    aa = load(auto_dir / "summary.json")
    bb = load(explicit_dir / "summary.json")
    assert aa["status"] == bb["status"] == "complete"
    assert aa["mnte"] == -1 and bb["mnte"] == 64
    assert aa["geometry"] == bb["geometry"] == "3D_simple_cubic"
    for k in ("source_revision","archive_sha256","binary_sha256",
              "sizes","sectors","threads","warmups","repeats",
              "fast_group_validation","anchored_basis","anchored_lookup",
              "local_bond_cache","direct_hop_table","basis_coord_cache",
              "basis_pruned_canon","skip_unused_hash","cache_profile"):
        assert aa[k] == bb[k], ("Control variable differs", k)
    assert aa["local_bond_cache"] and aa["basis_pruned_canon"]
    assert aa["cache_profile"] is False
    env_a, env_b = load(auto_dir / "environment.json"), load(explicit_dir / "environment.json")
    assert env_a["environment"]["pbs_jobid"] == env_b["environment"]["pbs_jobid"]
    assert env_a["environment"]["host"] == env_b["environment"]["host"]
    assert env_a["environment"]["mode"] == env_b["environment"]["mode"]
    assert env_a["binary_sha256"] == env_b["binary_sha256"]
    ar, br = records(auto_dir), records(explicit_dir)
    assert ar.keys() == br.keys()
    results = []
    for key in sorted(ar):
        a, b = ar[key], br[key]
        assert a["mnte"] == -1 and b["mnte"] == 64
        for k in ("side","sector","threads","period","repetition",
                  "builder","dimension","source_revision","binary_sha256",
                  "archive_sha256"):
            assert a[k] == b[k], (key, k)
        pa, pb = auto_dir / "trials" / key, explicit_dir / "trials" / key
        left = (pa / "input.dat").read_bytes()
        right = (pb / "input.dat").read_bytes()
        assert b"MNTE=-1" in left and b"MNTE=64" in right
        assert norm_input(left) == norm_input(right), ("Different physical input", key)
        assert a["all_input_sha256"].keys() == b["all_input_sha256"].keys()
        for fname in a["all_input_sha256"]:
            if fname != "input.dat":
                assert a["all_input_sha256"][fname] == b["all_input_sha256"][fname], (
                    "Changed physical input file", key, fname)
        assert a["dimension"] == b["dimension"]
        assert abs(a["energy"]-b["energy"]) <= 1e-10, ("Energy mismatch", key)
        assert 0 <= a["residual"] <= 1.01e-8
        assert 0 <= b["residual"] <= 1.01e-8
        assert a["h_applications"] == b["h_applications"], ("Lanczos operation count", key)
        results.append((a,b))
    phases = ("basis", "operator_cache", "lanczos", "total")
    cases = {}
    combos = sorted({(a["side"],a["sector"],a["threads"],a["builder"])
                    for a,b in results})
    for side, sector, threads, builder in combos:
        subset = [(a,b) for a,b in results if
                  a["side"] == side and a["sector"] == sector and
                  a["threads"] == threads and a["builder"] == builder and
                  a["period"] == "measure"]
        assert len(subset) == aa["repeats"]
        key = "s{}-{}-t{}-{}".format(side,sector,threads,builder)
        rec = {
            "dimension": subset[0][0]["dimension"],
            "energy": subset[0][0]["energy"],
            "h_applications": subset[0][0]["h_applications"],
            "phases_s": {}, "rss_kib": {}
        }
        for phase in phases:
            xa=medmad([a["phase_wall_s"][phase] for a,b in subset])
            xb=medmad([b["phase_wall_s"][phase] for a,b in subset])
            rec["phases_s"][phase] = {"automatic":xa, "explicit64":xb,
                "auto_over_explicit_median":xa["median"]/max(xb["median"],1e-12)}
        if all(a["peak_rss_kib"] is not None and b["peak_rss_kib"] is not None
               for a,b in subset):
            rec["rss_kib"] = {
                "automatic":medmad([a["peak_rss_kib"] for a,b in subset]),
                "explicit64":medmad([b["peak_rss_kib"] for a,b in subset])}
        cases[key] = rec
    return {"schema":1,"paired_audit":"PASS","publication_ready":False,
            "caution":"Same PBS allocation and auditable input pairing are not independent site operator attestation.",
            "job_id":env_a["environment"]["pbs_jobid"],
            "host":env_a["environment"]["host"],
            "source_revision":aa["source_revision"],"binary_sha256":aa["binary_sha256"],
            "warmups":aa["warmups"],"repeats":aa["repeats"],
            "validated_trials":len(results),
            "model":"3D cubic fixed-three-magnon Gamma/Oh exact paired input",
            "cases":cases}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--auto",type=Path,required=True)
    p.add_argument("--explicit64",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    assert not a.output.exists(), "Refusing to overwrite paired comparison"
    for root in (a.auto,a.explicit64):
        audit=load(root.parent / ("audit" if root.name == "measurements" else
                                    "audit-explicit64") / "audit.json")
        assert audit["audited"] is True, ("Audit not PASS", root)
        assert audit["source_revision"] == load(root / "summary.json")["source_revision"]
    comparison=compare(a.auto,a.explicit64)
    a.output.write_text(json.dumps(comparison,sort_keys=True,indent=2)+"\n")
    print("CPC_AUTO_MNTE_PAIRED_PASS job={} trials={} cases={}".format(
        comparison["job_id"],comparison["validated_trials"],len(comparison["cases"])))
    for key in sorted(comparison["cases"]):
        c=comparison["cases"][key]
        aT=c["phases_s"]["total"]["automatic"]["median"]
        eT=c["phases_s"]["total"]["explicit64"]["median"]
        print("{} auto={:.6f}s explicit64={:.6f}s ratio={:.4f}".format(
            key,aT,eT,aT/eT))


if __name__=="__main__":
    main()
