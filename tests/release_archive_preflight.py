#!/usr/bin/env python3
"""Verify release sources against Git and a complete make-dist archive."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile

ROOT_SOURCES = (
    "README.md", "LICENSE", "CITATION.cff", "configure.ac",
    "Makefile.am", "mkdocs.yml", ".readthedocs.yaml",
)
GENERATED_DIST = (
    "configure", "aclocal.m4", "config.h.in", "compile",
    "depcomp", "install-sh", "missing", "test-driver",
)
# Historical generated files are tracked in this repository. Autotools/compiler
# output may change, but authored files outside this allowlist must stay clean.
GENERATED_TRACKED = set(GENERATED_DIST) | {
    "config.h", "config.log", "config.status", "stamp-h1",
    "script_cubic_sp_HB/functions.mod",
}
SOURCE_PREFIXES = ("docs/", "tests/", "examples/ref_data/", "input/")


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def git_paths(root, *args):
    data = subprocess.check_output(["git", "-C", str(root), *args])
    return [p.decode("utf-8") for p in data.split(b"\0") if p]


def committed_executable_modes(root):
    modes = {}
    for record in subprocess.check_output(
            ["git", "-C", str(root), "ls-tree", "-r", "-z", "HEAD"]).split(b"\0"):
        if record:
            metadata, name = record.split(b"\t", 1)
            modes[name.decode("utf-8")] = bool(int(metadata.split()[0], 8) & 0o100)
    return modes


def generated_path(name):
    path = PurePosixPath(name)
    return (name in GENERATED_TRACKED or
            name.startswith("autom4te.cache/") or
            path.name in ("Makefile", "Makefile.in"))


def source_paths(tracked):
    """Release inputs, determined from Git rather than archive contents."""
    selected = set(ROOT_SOURCES)
    for name in tracked:
        path = PurePosixPath(name)
        if generated_path(name):
            continue
        if name.startswith(SOURCE_PREFIXES):
            selected.add(name)
        elif name.startswith("source/") and path.suffix in (".f90", ".am"):
            selected.add(name)
        elif name.startswith("script_") and (
                path.name in ("Makefile.am", "mk_input.f90", "run.sh") or
                path.suffix == ".dat"):
            selected.add(name)
    return sorted(selected)


def regular_bytes(root, name):
    path = root / name
    if not path.is_file() or path.is_symlink():
        raise ValueError("Missing regular source: " + name)
    return path.read_bytes()


def audit(root, tarball, expected_sha, release_tag):
    sha = git(root, "rev-parse", "HEAD")
    if expected_sha and sha != expected_sha:
        raise ValueError("HEAD and expected full Git SHA differ")
    tracked = set(git_paths(root, "ls-files", "-z"))
    modified = git_paths(root, "diff", "HEAD", "--name-only", "-z")
    authored_changes = [name for name in modified if not generated_path(name)]
    if authored_changes:
        raise ValueError("Authored tracked files modified after checkout: " +
                         ", ".join(authored_changes))
    expected = {}
    source_modes = committed_executable_modes(root)
    for name in source_paths(tracked):
        if name not in tracked:
            raise ValueError("Untracked required source: " + name)
        commit_bytes = subprocess.check_output(
            ["git", "-C", str(root), "show", "HEAD:" + name])
        if regular_bytes(root, name) != commit_bytes:
            raise ValueError("Tracked source modified after checkout: " + name)
        expected[name] = sha256(commit_bytes)

    # Compare the checker to HEAD before executing it.
    check = [sys.executable, str(root / "tests/check_release_metadata.py"),
             "--repo", str(root)]
    if release_tag:
        check += ["--release-tag", release_tag]
    done = subprocess.run(check, capture_output=True, text=True)
    if done.returncode:
        raise ValueError("Citation/version validation failed: " +
                         done.stdout.strip() + done.stderr.strip())
    metadata = json.loads(done.stdout)

    # Generated build inputs are necessarily toolchain-dependent. Require them
    # in the archive and compare them with this verified build tree, separately
    # from the authored-source hashes tied to the commit.
    generated_names = set(GENERATED_DIST)
    generated_names.update(
        str(PurePosixPath(name).with_name("Makefile.in"))
        for name in tracked if PurePosixPath(name).name == "Makefile.am")
    generated = {name: sha256(regular_bytes(root, name))
                 for name in sorted(generated_names)}
    executable_modes = {name: source_modes[name] for name in expected}
    executable_modes.update({name: bool((root / name).stat().st_mode & 0o100)
                             for name in generated})
    archive_hashes = {}
    with tarfile.open(tarball, "r:*") as tar:
        members = tar.getmembers()
        if not members:
            raise ValueError("Empty archive")
        parsed = [PurePosixPath(m.name) for m in members]
        if any(p.is_absolute() or ".." in p.parts for p in parsed):
            raise ValueError("Unsafe archive path")
        roots = {p.parts[0] for p in parsed if p.parts}
        if len(roots) != 1:
            raise ValueError("Source archive must have one package root")
        root_dir = next(iter(roots))
        index = {}
        for m, path in zip(members, parsed):
            key = str(path)
            if key in index:
                raise ValueError("Duplicate archive member " + key)
            if not (m.isfile() or m.isdir()):
                raise ValueError("Unsupported archive member: " + key)
            index[key] = m
        for name, wanted in {**expected, **generated}.items():
            m = index.get(root_dir + "/" + name)
            if m is None or not m.isfile():
                raise ValueError("Missing regular archive member: " + name)
            if not m.mode & 0o400 or bool(m.mode & 0o100) != executable_modes[name]:
                raise ValueError("Archive read/execute mode mismatch for " + name)
            stream = tar.extractfile(m)
            if stream is None or sha256(stream.read()) != wanted:
                raise ValueError("Archive byte mismatch for " + name)
        for key, member in index.items():
            if not member.isfile():
                continue
            relative = str(PurePosixPath(key).relative_to(root_dir))
            if relative not in expected and relative not in generated:
                raise ValueError("Unclassified archive file: " + relative)
            archive_hashes[relative] = expected.get(relative, generated.get(relative))
    return {
        "status": "PASS",
        "scope": "all declared release inputs/archive byte identity, NOT physics",
        "git_sha": sha, "citation": metadata, "archive_name": tarball.name,
        "archive_sha256": sha256(tarball.read_bytes()), "package_root": root_dir,
        "verified_files": len(expected), "source_sha256": expected,
        "generated_files": len(generated), "generated_sha256": generated,
        "executable_paths": sorted(name for name, mode in executable_modes.items() if mode),
        "archive_files": len(archive_hashes),
        "postbuild_modified_tracked_paths": modified,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--expected-commit", help="Full 40-character Git SHA")
    ap.add_argument("--release-tag", help="Validate vX.Y.Z metadata, never publish a tag")
    ap.add_argument("--json-out", type=Path)
    args = ap.parse_args()
    try:
        if args.expected_commit and (len(args.expected_commit) != 40 or
             any(c not in "0123456789abcdef" for c in args.expected_commit)):
            raise ValueError("Expected commit must be a full lowercase Git SHA")
        result = audit(args.repo.resolve(), args.archive.resolve(),
                       args.expected_commit, args.release_tag)
    except (OSError, ValueError, subprocess.CalledProcessError, tarfile.TarError,
            json.JSONDecodeError) as exc:
        result = {"status": "FAIL", "error": str(exc)}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(encoded, end="")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(encoded)
    return int(result["status"] != "PASS")


if __name__ == "__main__":
    raise SystemExit(main())
