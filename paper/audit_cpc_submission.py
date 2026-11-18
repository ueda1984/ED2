#!/usr/bin/env python3
"""Produce a conservative CPC CPiP pre-submission gate report for ED2.

This is a source/metadata audit, not journal acceptance prediction.
In default reporting mode it succeeds even when release decisions remain
open. Use --require-ready only for a consciously frozen release.
Python >= 3.6, standard library only.
"""
import argparse
import json
from pathlib import Path
import re
import sys


def parse_version(text, expression, origin):
    found = re.search(expression, text, re.MULTILINE)
    if not found:
        raise ValueError("Cannot parse version from {}".format(origin))
    return found.group(1)


def examine(repo):
    read = lambda name: (repo / name).read_text(encoding="utf-8")
    tex = read("paper/cpc_ed2_main.tex")
    conf = read("configure.ac")
    citation = read("CITATION.cff")
    license_file = read("LICENSE")
    makefile = read("source/Makefile.am")

    autotools = parse_version(
        conf, r"AC_INIT\(\[QS3-ED2\],\s*\[([^\]]+)\]", "configure.ac")
    cff = parse_version(
        citation, r"^version:\s*[\x22\x27]?([0-9A-Za-z._-]+)",
        "CITATION.cff")
    released_date = re.search(r"^date-released:\s*(.+)$",
                              citation, re.MULTILINE)

    authors_approved = (
        "Author list pending project approval" not in tex
        and "Affiliations to be confirmed before submission" not in tex
    )
    version_approved = (
        autotools == cff
        and "To be assigned after integration and release review" not in tex
    )
    accession_assigned = (
        "Program Library accession/DOI:] Not yet assigned." not in tex
        and "CPC Program Library accession, software DOI and complete" not in tex
    )

    report = {
        "purpose": "CPC CPiP ED2 submission gate audit",
        "article_type_proposed": "Computer Programs in Physics (CPiP)",
        "publication_ready": False,
        "project": "ED2 (QS3-ED2)",
        "source_version": {
            "autotools": autotools,
            "citation_cff": cff,
            "citation_release_date_claim": released_date.group(1).strip()
            if released_date else None,
            "versions_consistent": autotools == cff,
            "release_date_event_verified": False,
        },
        "gates": {
            "source_open_license": {
                "pass": license_file.startswith("MIT License"),
                "detail": "OSI-compatible MIT license source exists; full contributor/copyright review is manual.",
            },
            "executable_not_public_library_abi": {
                "pass": "bin_PROGRAMS" in makefile,
                "detail": "ED2 ships an executable. Do not advertise an installed Fortran library interface.",
            },
            "scientific_manuscript_drafted": {
                "pass": all(x in tex for x in (
                    r"\begin{abstract}", r"\section{Conclusions}",
                    r"\section*{Program summary", r"\label{tab:pr17-pilot}"
                )),
                "detail": "A complete manuscript draft exists; scientific/visual editorial review remains open.",
            },
            "pilot_table_and_checker_present": {
                "pass": all((repo / x).is_file() for x in (
                    "paper/data/pr17_presto_coordinate_cache_pilot.csv",
                    "paper/check_pr17_presto_table.py",
                    "paper/data/README.md",
                )),
                "detail": "Only published pilot medians are staged, not full raw PRESTO trial evidence.",
            },
            "author_credit_approved": {
                "pass": authors_approved,
                "detail": "Need approved complete author/order/affiliation/CRediT and historical software contributor credits.",
            },
            "release_version_consistent_and_approved": {
                "pass": version_approved,
                "detail": "Synchronize configure.ac, CITATION.cff, release tag and manuscript; check actual release date.",
            },
            "official_upstream_release_frozen": {
                "pass": False,
                "detail": "Requires confirmed QS-Cube/ED2 integration, tagged source archive and checksums (no unapproved merge/tag).",
            },
            "cpc_program_library_and_software_archive": {
                "pass": accession_assigned,
                "detail": "Do not invent a CPC Program Library accession or a Mendeley Data / software DOI.",
            },
            "raw_trial_provenance_and_operator_attestation": {
                "pass": False,
                "detail": "PR17 jobs 22282/22322 passed 48-pair physics checks, but exclusivity operator review is false and raw trials are not archived with the paper.",
            },
            "representative_input_output_and_clean_build": {
                "pass": False,
                "detail": "Freeze at least three example input/output cases and verify build, make check, make distcheck on immutable release.",
            },
            "paper_proof_and_submission_declarations": {
                "pass": False,
                "detail": "Requires final human page-wise PDF/figures proof, cover letter, CRediT, conflicts, funding and data statement.",
            },
        },
    }
    report["blockers"] = [
        key for key, entry in report["gates"].items()
        if not entry["pass"]
    ]
    report["publication_ready"] = not report["blockers"]
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path("."))
    p.add_argument("--output", type=Path)
    p.add_argument("--require-ready", action="store_true",
                   help="Refuse if any gate remains blocked")
    args = p.parse_args()
    report = examine(args.repo.resolve(strict=True))
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print("CPC SUBMISSION AUDIT: {} unfulfilled gates".format(
          len(report["blockers"])))
    for key, entry in report["gates"].items():
        print("  {} {}".format("PASS" if entry["pass"] else "BLOCK", key))
    print("Versions: configure.ac={} / CITATION.cff={}".format(
          report["source_version"]["autotools"],
          report["source_version"]["citation_cff"]))
    print("Submission status: {}".format(
          "READY" if report["publication_ready"] else
          "NOT READY (author/release decisions outstanding)"))
    if args.require_ready and not report["publication_ready"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
