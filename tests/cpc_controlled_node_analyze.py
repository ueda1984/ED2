#!/usr/bin/env python3
"""Numerically strict post-processing for a complete CPU-pinned ED2 cache study.

Inputs are the *existing* cpc_controlled_node_benchmark.py result directory,
which contains environment.json, input-manifest.json, summary.json and one
paired-ompN/report.json per thread setting. Independently audit the factorial
experiment, numerical residual/physical acceptance, input SHA256 and
source/binary provenance before calculating matched-policy robust statistics.

A paired cache policy speed ratio is REPORTED, but is NOT presented as a
publication-quality speedup unless operator-attested exclusivity, >=5 measured
repeats, a source revision/archive and a complete run have been verified.
Even in that case, a single node allocation cannot rule out drift or establish
universal speedup; publish uncertainty and hardware configuration.

No NumPy or SciPy dependency is needed for this postprocessor.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

POLICIES=("half_generic","full_tight","full_generic")
COLUMNS=("n_sites","omp","alg","policy","count","wall_median_s",
         "wall_mad_s","rss_median_mib","rss_mad_mib","h_actions_median",
         "nominal_cache_mib","wall_ratio_half_over_tight",
         "rss_ratio_half_over_tight","measurement_class")


def load(path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as src:
        for piece in iter(lambda:src.read(1024*1024),b""):
            h.update(piece)
    return h.hexdigest()


def finite_number(value,positive=False):
    if not isinstance(value,(int,float)):
        raise ValueError("Expected a numeric value")
    x=float(value)
    if not math.isfinite(x) or (positive and x<=0):
        raise ValueError("Non-finite/nonpositive numerical observation")
    return x


def median_mad(xs):
    values=[finite_number(x,positive=True) for x in xs]
    if not values:raise ValueError("No measurements to summarize")
    med=statistics.median(values)
    return med,statistics.median(abs(y-med) for y in values)


def expected_cache_bytes(n,policy):
    dim=6**(n//2)
    default=10*n-7
    tight=7*n//2
    if policy=="half_generic":return 20*default*(dim//2)
    if policy=="full_tight":return 20*tight*dim
    if policy=="full_generic":return 20*default*dim
    raise ValueError("Unexpected cache policy")


def audit(root,publication_requested=False):
    environment=load(root/"environment.json")
    summary=load(root/"summary.json")
    manifest=load(root/"input-manifest.json")
    if environment.get("schema")!=1 or summary.get("schema")!=1:
        raise ValueError("Unrecognized wrapper record schema")
    if summary.get("mode")!=environment.get("benchmark_mode"):
        raise ValueError("Contradictory shared/exclusive provenance")
    if not summary.get("pass_complete"):
        raise ValueError("Source benchmark reported numerical or harness failure")
    if summary["metadata_sha256"]!=sha(root/"environment.json"):
        raise ValueError("Source environment SHA256 mismatch")
    if summary["input_manifest_sha256"]!=sha(root/"input-manifest.json"):
        raise ValueError("Input manifest itself was mutated")
    if not isinstance(manifest,list) or len(manifest)!=summary["input_file_count"] or not manifest:
        raise ValueError("Incomplete source input manifest")
    root_real=root.resolve()
    counted=set()
    for item in manifest:
        rel=Path(item["relative_path"])
        if rel.is_absolute() or ".." in rel.parts or not str(rel).startswith("paired-omp"):
            raise ValueError("Invalid input-manifest path")
        if str(rel) in counted:raise ValueError("Duplicate input-manifest entry")
        counted.add(str(rel))
        file=(root/rel).resolve(strict=True)
        if not file.is_relative_to(root_real):
            raise ValueError("Input file escapes archived benchmark tree")
        if file.stat().st_size!=item["size_bytes"] or sha(file)!=item["sha256"]:
            raise ValueError("Input file contents failed SHA256 verification: "+str(rel))
    allowed=set(environment["allowed_affinity"])
    selected=set(environment["selected_affinity"])
    if not selected or not selected.issubset(allowed):
        raise ValueError("Invalid selected CPU affinity")
    if set(environment["process_affinity"])!=selected:
        raise ValueError("Actual CPU affinity not equal to requested set")
    modes=("operator_verified_exclusive","shared_or_exclusivity_unverified")
    mode=environment["benchmark_mode"]
    if mode not in modes:raise ValueError("Unrecognized operator-attestation mode")
    if bool(environment["exclusive_operator_attestation"])!=(mode==modes[0]):
        raise ValueError("Exclusive label/attestation inconsistent")
    if publication_requested and mode!=modes[0]:
        raise ValueError("Publication timing requires verified dedicated-node operator attestation")
    if publication_requested and not (environment.get("source_revision") or
                                      environment.get("source_archive_sha256")):
        raise ValueError("Publication-mode benchmark source revision is unrecorded")
    if len(environment.get("executable_sha256",""))!=64:
        raise ValueError("Missing executable SHA256")
    params=environment["parameters"]
    sizes=params["sizes"];algs=params["algs"];omps=params["threads"]
    repeats=int(params["repeats"]);warmups=int(params["warmups"])
    if not sizes or not algs or not omps or len(set(omps))!=len(omps):
        raise ValueError("Missing/duplicate experiment axes")
    if sorted(set(omps))!=sorted(x["omp"] for x in summary["tasks"]):
        raise ValueError("Per-OMP tasks do not match requested threads")
    if max(omps)>len(selected):
        raise ValueError("Number of requested OpenMP threads exceeds bound cores")
    if publication_requested and repeats<5:
        raise ValueError("Fewer than five measured repeats do not justify paper timing")
    accepted=[]
    groups=defaultdict(list)
    for task in summary["tasks"]:
        omp=int(task["omp"])
        if task.get("exit_code")!=0 or not task.get("pass_complete"):
            raise ValueError("Source numerical job failed at OMP "+str(omp))
        report_path=root/f"paired-omp{omp}"/"report.json"
        data=load(report_path)
        p=data["parameters"]
        if not data["pass_complete"] or data["violations"] or data["numerical_nonconvergence"]!=0:
            raise ValueError("Underlying independent physics checks failed")
        if p["sizes"]!=sizes or p["algs"]!=algs or int(p["omp"])!=omp:
            raise ValueError("Underlying physical target matrix changed")
        if (p["warmups"],p["repeats"],p["NOE"],p["D"])!=(warmups,repeats,1,params["D"]):
            raise ValueError("Experiment/warmup/one-eigenpair mismatch")
        if p["residual_atol"]!=params["residual_atol"] or not p["explicit_initial"]:
            raise ValueError("Residual request/initial-vector contract changed")
        records=data["raw"]
        wanted=len(sizes)*len(algs)*len(POLICIES)*(warmups+repeats)
        if len(records)!=wanted or data["completed"]!=wanted or data["expected"]!=wanted:
            raise ValueError("Missing/extra numerical trials at OMP "+str(omp))
        keys=set()
        for x in records:
            n=int(x["n_sites"]);alg=int(x["alg"]);policy=x["policy"]
            rep=int(x["repeat"]);phase=x["phase"]
            if n not in sizes or alg not in algs or policy not in POLICIES:
                raise ValueError("Unexpected model/solver/cache condition")
            if phase!=("warmup" if rep<=warmups else "measurement") or not 1<=rep<=warmups+repeats:
                raise ValueError("Warmup was mislabeled or repetition missing")
            key=(n,omp,alg,policy,rep)
            if key in keys: raise ValueError("Duplicate measured/repeated trial")
            keys.add(key)
            if x["status"]!=0 or x["exit_code"]!=0:
                raise ValueError("Nonconvergent ED2 trial cannot appear in accepted table")
            if x["threads"]!=omp or x["wanted"]!=1:
                raise ValueError("Thread count/number of requested states differs")
            if x["nominal_cache_bytes"]!=expected_cache_bytes(n,policy):
                raise ValueError("Unexpected cache allocation")
            if int(x["reported_optimal_mnte"])!=7*n//2-1:
                raise ValueError("Source-calibrated tight cache capacity inconsistent")
            if x["seed"]!=19 or x["initial_sha256"]!=data["references"][sizes.index(n)]["initial_sha256"]:
                raise ValueError("Seed/starting vector modified during cache comparison")
            if x["residual_atol"]!=params["residual_atol"]:
                raise ValueError("Absolute residual request differs across cache policies")
            if finite_number(x["residual"])>params["residual_atol"]*(1+1e-10):
                raise ValueError("Residual gate exceeded")
            if finite_number(x["energy_error"])>3e-7 or \
                finite_number(x["magnetization_max_error"])>8e-6 or \
                finite_number(x["correlation_max_error"])>8e-6:
                raise ValueError("Independent physical comparison no longer accepted")
            for field in ("wall_seconds","h_applications","peak_rss_kib"):
                finite_number(x[field],positive=True)
            if phase=="measurement":
                groups[(n,omp,alg,policy)].append(x)
            accepted.append(x)
        required={(n,omp,alg,pol,rep) for n in sizes for alg in algs
                  for pol in POLICIES for rep in range(1,warmups+repeats+1)}
        if keys!=required:raise ValueError("Factorial measurement matrix incomplete")
    expected_group={(n,omp,alg,pol) for n in sizes for omp in omps
                    for alg in algs for pol in POLICIES}
    if set(groups)!=expected_group or any(len(g)!=repeats for g in groups.values()):
        raise ValueError("Incomplete policy-solver-size measured matrix")
    quality=("operator_attested_node_replicated" if mode==modes[0] and
             repeats>=5 else "exploratory_shared_or_insufficient_repeats")
    rows=[]
    stats={}
    for (n,omp,alg,policy),samples in sorted(groups.items()):
        wall,wallmad=median_mad(x["wall_seconds"] for x in samples)
        rss,rssmad=median_mad(x["peak_rss_kib"] for x in samples)
        actions,_=median_mad(x["h_applications"] for x in samples)
        stats[(n,omp,alg,policy)]=(wall,rss)
        rows.append(dict(n_sites=n,omp=omp,alg=alg,policy=policy,count=len(samples),
                         wall_median_s=wall,wall_mad_s=wallmad,
                         rss_median_mib=rss/1024,rss_mad_mib=rssmad/1024,
                         h_actions_median=actions,
                         nominal_cache_mib=expected_cache_bytes(n,policy)/(1024**2),
                         wall_ratio_half_over_tight=None,rss_ratio_half_over_tight=None,
                         measurement_class=quality))
    for row in rows:
        n=row["n_sites"];omp=row["omp"];alg=row["alg"]
        half=stats[(n,omp,alg,"half_generic")]
        tight=stats[(n,omp,alg,"full_tight")]
        row["wall_ratio_half_over_tight"]=half[0]/tight[0]
        row["rss_ratio_half_over_tight"]=half[1]/tight[1]
    return dict(schema=1,source_mode=mode,measurement_class=quality,
                author_attested_exclusivity=(mode==modes[0]),
                publication_gate_passed=(quality=="operator_attested_node_replicated"),
                caution=("Operator attestation does not guarantee CPU clock, "
                         "NUMA or contention control. Data from a single "
                         "allocation do not certify a universal speedup."),
                sizes=sizes,algs=algs,threads=omps,warmups=warmups,
                repeats=repeats,observed_total_records=len(accepted),
                input_manifest_entries=len(manifest),
                binary_sha256=environment["executable_sha256"],
                source_revision=environment.get("source_revision"),
                cache_formula="two cache arrays, integer4 + complex16 per stored entry",
                rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,required=True,
                   help="Previously finished cpc_controlled_node_benchmark.py output root")
    p.add_argument("--output",type=Path,required=True,
                   help="New, nonexistent output directory for audited paper tables")
    p.add_argument("--require-publication",action="store_true",
                   help="Refuse shared runs or <5 measured repeats")
    args=p.parse_args()
    root=args.root.resolve(strict=True)
    report=audit(root,publication_requested=args.require_publication)
    out=args.output.resolve()
    out.mkdir(exist_ok=False,parents=True)
    (out/"audit.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    with (out/"paired-policy-stats.csv").open("w",newline="",encoding="utf-8") as fp:
        writer=csv.DictWriter(fp,fieldnames=COLUMNS)
        writer.writeheader()
        for row in report["rows"]:
            writer.writerow(row)
    print("CONTROLLED NODE POSTPROCESS:",report["observed_total_records"],"all cases audited;",
          report["measurement_class"],"; rows",len(report["rows"]),flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
