#!/usr/bin/env python3
"""Negative release tests: candidate transitions and complete source archives."""
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


metadata = load("check_release_metadata")
preflight = load("release_archive_preflight")


class ReleaseMetadataTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "configure.ac").write_text("AC_INIT([QS3-ED2], [1.0.0])\n")
        self.citation()

    def citation(self, extra=""):
        (self.root / "CITATION.cff").write_text("cff-version: 1.2.0\n" + extra)

    def test_unreleased_development(self):
        self.assertEqual(metadata.validate(self.root, None)["mode"], "unreleased_development")

    def test_checked_out_release_tag_metadata_if_versioned(self):
        # Tag-mode qualification must be tested against the actual checkout,
        # not only synthetic CFF fixtures. Development/candidate branches
        # without a release date are intentionally exempt from tag mode.
        checkout = HERE.parent
        current = metadata.validate(checkout, None)
        if current["mode"] != "versioned_metadata":
            self.skipTest("Full version/date metadata not yet present")
        tag = "v" + current["version"]
        tagged = metadata.validate(checkout, tag)
        self.assertEqual(tagged["mode"], "tag_metadata")
        self.assertEqual(tagged["version"], current["version"])
        self.assertEqual(tagged["date_released"], current["date_released"])

    def test_approved_candidate_without_pretend_publication_date(self):
        self.citation("version: 1.0.0\n")
        candidate = metadata.validate(self.root, None)
        self.assertEqual(candidate["mode"], "release_candidate")
        self.assertEqual(candidate["version"], "1.0.0")
        self.assertIsNone(candidate["date_released"])
        with self.assertRaisesRegex(ValueError, "requires an actual CFF date-released"):
            metadata.validate(self.root, "v1.0.0")

    def test_candidate_version_must_equal_package_version(self):
        self.citation("version: 1.0.1\n")
        with self.assertRaisesRegex(ValueError, "mismatch"):
            metadata.validate(self.root, None)

    def test_candidate_then_tag_same_metadata(self):
        self.citation("version: 1.0.0\ndate-released: 2026-11-20\n")
        self.assertEqual(metadata.validate(self.root, None)["mode"], "versioned_metadata")
        self.assertEqual(metadata.validate(self.root, "v1.0.0")["mode"], "tag_metadata")

    def test_date_without_version_is_rejected(self):
        self.citation("date-released: 2026-11-20\n")
        with self.assertRaisesRegex(ValueError, "requires version"):
            metadata.validate(self.root, None)

    def test_tag_requires_metadata(self):
        with self.assertRaises(ValueError):
            metadata.validate(self.root, "v1.0.0")

    def test_wrong_version_or_tag(self):
        self.citation("version: 0.1.0\ndate-released: 2026-11-20\n")
        with self.assertRaisesRegex(ValueError, "mismatch"):
            metadata.validate(self.root, None)
        self.citation("version: 1.0.0\ndate-released: 2026-11-20\n")
        with self.assertRaisesRegex(ValueError, "mismatch"):
            metadata.validate(self.root, "v1.0.1")

    def test_invalid_date(self):
        for date in ("2026-02-30", "20261120", "tomorrow"):
            with self.subTest(date=date):
                self.citation("version: 1.0.0\ndate-released: " + date + "\n")
                with self.assertRaises(ValueError):
                    metadata.validate(self.root, None)

    def test_balanced_quoted_scalars_and_comments(self):
        self.citation('version: "1.0.0" # approved candidate\ndate-released: \'2026-11-20\'\n')
        self.assertEqual(metadata.validate(self.root, "v1.0.0")["version"], "1.0.0")

    def test_unmatched_quotes_rejected(self):
        for value in ('"1.0.0', "'1.0.0", '"1.0.0\''):
            with self.subTest(value=value):
                self.citation('version: ' + value + '\ndate-released: 2026-11-20\n')
                with self.assertRaisesRegex(ValueError, "quoted scalar"):
                    metadata.validate(self.root, "v1.0.0")


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "source"
        self.root.mkdir()
        files = {
            "README.md": "Release test fixture\n",
            "RELEASE_NOTES_v1.0.0.md": "v1.0.0 notes fixture\n",
            "LICENSE": "MIT\n",
            "CITATION.cff": "cff-version: 1.2.0\n",
            "configure.ac": "AC_INIT([QS3-ED2], [1.0.0])\n",
            "Makefile.am": "SUBDIRS = source tests script_chain\n",
            "mkdocs.yml": "site_name: Test\n", ".readthedocs.yaml": "version: 2\n",
            "source/Makefile.am": "bin_PROGRAMS = QS3ED2\n",
            "source/main.f90": "! main fixture\n",
            "source/ham2vec.f90": "! essential numerical source\n",
            "source/eigen_solver.f90": "! eigensolver fixture\n",
            "source/new_scientific_module.f90": "! future modules must be covered\n",
            "tests/Makefile.am": "EXTRA_DIST = check_release_metadata.py\n",
            "tests/check_release_metadata.py": (HERE / "check_release_metadata.py").read_text(),
            "tests/release_archive_preflight.py": (HERE / "release_archive_preflight.py").read_text(),
            "docs/getting-started/installation.md": "# Installation\n",
            "script_chain/Makefile.am": "noinst_PROGRAMS = mk_input\n",
            "script_chain/mk_input.f90": "! generator fixture\n",
            "script_chain/input_param.dat": "input fixture\n",
            "script_chain/run.sh": "#!/bin/sh\n",
            "examples/ref_data/chain/eigenvalues.dat": "1 -0.75\n",
            "input/input.dat": "input fixture\n",
            "configure": "# generated old configure\n",
            "config.log": "old generated log\n",
        }
        for name, content in files.items():
            self.write(name, content)
        (self.root / "configure").chmod(0o755)
        (self.root / "script_chain/run.sh").chmod(0o755)
        self.git("init", "-q")
        self.git("add", ".")
        self.git("-c", "user.name=ED2 test fixture", "-c",
                 "user.email=ed2-test@example.invalid", "commit", "-qm", "fixture")
        self.sha = self.git("rev-parse", "HEAD").strip()
        # A build updates generated files, as real autoreconf/configure does.
        generated = set(preflight.GENERATED_DIST) | {
            "Makefile.in", "source/Makefile.in", "tests/Makefile.in",
            "script_chain/Makefile.in",
        }
        for name in generated:
            self.write(name, "# generated current build\n")
        self.write("config.log", "new generated log\n")
        self.names = sorted(set(files) - {"config.log"} | generated)
        self.archive = Path(self.tmp.name) / "source.tar.gz"

    def git(self, *args):
        return subprocess.check_output(
            ["git", "-C", str(self.root), *args], text=True, stderr=subprocess.DEVNULL)

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)

    def pack(self, omit=(), corrupt=None, extra=None, strip_execute=False):
        with tarfile.open(self.archive, "w:gz") as tar:
            for name in self.names:
                if name in omit:
                    continue
                data = (self.root / name).read_bytes()
                if name == corrupt:
                    data += b"! unexpected bytes\n"
                info = tarfile.TarInfo("qs3-ed2-1.0.0/" + name)
                info.size = len(data)
                info.mode = (self.root / name).stat().st_mode & 0o777
                if strip_execute:
                    info.mode &= ~0o111
                tar.addfile(info, io.BytesIO(data))
            if extra:
                info = tarfile.TarInfo("qs3-ed2-1.0.0/" + extra)
                info.size = 1
                tar.addfile(info, io.BytesIO(b"x"))

    def audit(self):
        return preflight.audit(self.root, self.archive, self.sha, None)

    def test_complete_archive_and_generated_changes(self):
        self.pack()
        result = self.audit()
        self.assertEqual(result["status"], "PASS")
        self.assertIn("source/new_scientific_module.f90", result["source_sha256"])
        self.assertIn("tests/check_release_metadata.py", result["source_sha256"])
        self.assertIn("configure", result["generated_sha256"])

    def test_missing_versioned_release_notes_rejected(self):
        self.pack(omit=("RELEASE_NOTES_v1.0.0.md",))
        with self.assertRaisesRegex(ValueError, "RELEASE_NOTES_v1.0.0.md"):
            self.audit()

    def test_tampered_versioned_release_notes_rejected(self):
        self.pack(corrupt="RELEASE_NOTES_v1.0.0.md")
        with self.assertRaisesRegex(ValueError, "byte mismatch"):
            self.audit()

    def test_missing_solver_rejected(self):
        self.pack(omit=("source/ham2vec.f90",))
        with self.assertRaisesRegex(ValueError, "source/ham2vec.f90"):
            self.audit()

    def test_tampered_solver_archive_rejected(self):
        self.pack(corrupt="source/eigen_solver.f90")
        with self.assertRaisesRegex(ValueError, "byte mismatch"):
            self.audit()

    def test_modified_scientific_worktree_rejected(self):
        self.write("source/ham2vec.f90", "! uncommitted scientific change\n")
        self.pack()
        with self.assertRaisesRegex(ValueError, "modified after checkout"):
            self.audit()

    def test_modified_checker_rejected_before_execution(self):
        self.write("tests/check_release_metadata.py", "raise RuntimeError('must not run')\n")
        self.pack()
        with self.assertRaisesRegex(ValueError, "modified after checkout"):
            self.audit()

    def test_stale_generated_archive_rejected(self):
        self.pack(corrupt="configure")
        with self.assertRaisesRegex(ValueError, "byte mismatch"):
            self.audit()

    def test_missing_generated_build_input_rejected(self):
        self.pack(omit=("source/Makefile.in",))
        with self.assertRaisesRegex(ValueError, "source/Makefile.in"):
            self.audit()

    def test_wrong_commit_rejected(self):
        self.pack()
        with self.assertRaisesRegex(ValueError, "HEAD"):
            preflight.audit(self.root, self.archive, "0" * 40, None)

    def test_unclassified_member_rejected(self):
        self.pack(extra="source/unverified.f90")
        with self.assertRaisesRegex(ValueError, "Unclassified"):
            self.audit()

    def test_archive_without_executable_bits_rejected(self):
        self.pack(strip_execute=True)
        with self.assertRaisesRegex(ValueError, "mode mismatch"):
            self.audit()


if __name__ == "__main__":
    unittest.main(verbosity=2)

