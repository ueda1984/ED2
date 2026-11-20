#!/usr/bin/env python3
"""Validate development, candidate, or tagged metadata without publishing.

Without a tag: absent CFF fields mean development; a version without an
actual release date means an approved, UNPUBLISHED candidate; version+date
is valid complete metadata after scheduling/final approval.
--release-tag requires both fields and an exact vMAJOR.MINOR.PATCH match.
No check creates or proves publication of a tag, Release, or DOI.
"""
import argparse
from datetime import date
import json
from pathlib import Path
import re


def field(text, name):
    values = re.findall("^" + re.escape(name) + r":[ \t]*(.*?)[ \t]*$", text, re.M)
    if len(values) > 1:
        raise ValueError(f"Duplicate {name}")
    if not values:
        return None
    value = values[0]
    # Only simple scalar fields are supported here, not a full CFF/YAML
    # schema validator. Never turn an unmatched quote into valid metadata.
    if value.startswith(('"', "'")):
        quote = value[0]
        match = re.fullmatch(re.escape(quote) + "([^" + quote + r"\n]*)" +
                             re.escape(quote) + r"(?:[ \t]+#.*)?", value)
        if not match:
            raise ValueError(f"Invalid quoted scalar for {name}")
        return match.group(1)
    value = re.split(r"[ \t]+#", value, maxsplit=1)[0].strip()
    if not value or '"' in value or "'" in value:
        raise ValueError(f"Invalid scalar for {name}")
    return value


def validate(root, tag):
    configure = (root / "configure.ac").read_text()
    citation = (root / "CITATION.cff").read_text()
    match = re.search(r"^AC_INIT\(\s*\[([^]]+)\]\s*,\s*\[([^]]+)\]", configure, re.M)
    if not match or match.group(1) != "QS3-ED2":
        raise ValueError("Cannot identify expected Autotools package/version")
    package_version = match.group(2)
    if field(citation, "cff-version") != "1.2.0":
        raise ValueError("Expected CFF format 1.2.0")
    version, released = field(citation, "version"), field(citation, "date-released")
    if tag and not re.fullmatch(r"v[0-9]+\.[0-9]+\.[0-9]+", tag):
        raise ValueError("Release tag must be vMAJOR.MINOR.PATCH")
    if version is None and released is None:
        if tag:
            raise ValueError("Tagged release requires CFF version and date-released")
        return {"mode": "unreleased_development",
                "autotools_package_version": package_version,
                "citation": "Use the full Git commit, not this package number"}
    if version is None:
        raise ValueError("CFF date-released requires version")
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        raise ValueError("CFF version must be MAJOR.MINOR.PATCH")
    if package_version != version or (tag and version != tag[1:]):
        raise ValueError(f"Version mismatch: tag={tag}, AC_INIT={package_version}, CFF={version}")
    if released is None:
        if tag:
            raise ValueError("Tagged release requires an actual CFF date-released")
        return {"mode": "release_candidate", "version": version,
                "date_released": None,
                "publication": "Not published; cite the full Git commit"}
    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", released):
        raise ValueError("CFF date-released must use YYYY-MM-DD")
    try:
        date.fromisoformat(released)
    except ValueError as error:
        raise ValueError(f"Invalid CFF release date: {released}") from error
    return {"mode": "tag_metadata" if tag else "versioned_metadata",
            "version": version, "date_released": released,
            "publication": "Not verified by this metadata check"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--release-tag", help="Prospective vX.Y.Z tag, not an action")
    args = parser.parse_args()
    try:
        result = {"status": "PASS", **validate(args.repo, args.release_tag)}
    except (OSError, ValueError) as error:
        result = {"status": "FAIL", "error": str(error)}
    print(json.dumps(result, indent=2, sort_keys=True))
    return int(result["status"] != "PASS")


if __name__ == "__main__":
    raise SystemExit(main())
