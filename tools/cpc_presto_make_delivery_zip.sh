#!/usr/bin/env bash
# Standalone immutable source ZIP: candidate cubic auto-MNTE + frozen original
# PR21 + latest-main snapshot + cross-library comparison tools.
# No merges, no scheduler, no PRESTO connectivity; never overwrite output.
set -euo pipefail
umask 077
if [[ $# -ne 2 ]]; then
  echo "Usage: bash tools/cpc_presto_make_delivery_zip.sh /path/to/ED2 /absolute/output.zip" >&2
  exit 2
fi
repo=$(cd -- "$1" && pwd -P)
target=$2
[[ "$target" == /* ]] || { echo "Output path must be absolute" >&2; exit 2; }
[[ ! -e "$target" ]] || { echo "Refusing existing ZIP: $target" >&2; exit 2; }
git -C "$repo" rev-parse --is-inside-work-tree | grep -qx true
[[ -z $(git -C "$repo" status --porcelain --untracked-files=no) ]] || {
  echo "Refusing dirty source tree" >&2; exit 2;
}
source_sha=$(git -C "$repo" rev-parse 'HEAD^{commit}')
[[ "$source_sha" =~ ^[0-9a-f]{40}$ ]] || exit 2
git -C "$repo" cat-file -e "$source_sha:tests/cpc_cubic_auto_mnte_regression.py"
git -C "$repo" cat-file -e "$source_sha:tools/cpc_presto_compare_mnte.py"
git -C "$repo" cat-file -e "$source_sha:tools/presto_cubic_exclusive.pbs"
git -C "$repo" cat-file -e "$source_sha:docs/development/cpc-presto-auto-mnte-plan-20261119.md"

tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
bundle="$tmp/ED2_CPC_AUTO_MNTE_PRESTO_20261119"
mkdir -- "$bundle"
mkdir -- "$bundle/cubic_auto_mnte" "$bundle/experimental_cubic_previous_PR21" \
        "$bundle/latest_main_snapshot" "$bundle/cross_library_comparator"
cp -- "$repo/docs/development/cpc-presto-auto-mnte-plan-20261119.md" \
   "$bundle/README_PRESTO_JP.md"
git -C "$repo" archive --format=tar "$source_sha" | gzip -n \
  > "$bundle/cubic_auto_mnte/ed2-source.tar.gz"
printf '%s\n' "$source_sha" > "$bundle/cubic_auto_mnte/source-commit.txt"
git -C "$repo" show "$source_sha:tools/presto_cubic_exclusive.pbs" \
  > "$bundle/cubic_auto_mnte/presto_cubic_exclusive.pbs"
(cd "$bundle/cubic_auto_mnte" && sha256sum ed2-source.tar.gz > ed2-source.sha256)
bash -n "$bundle/cubic_auto_mnte/presto_cubic_exclusive.pbs"

# Copy the earlier PRESTO source EXACTLY: never assume the experimental PR
# branch still points at the snapshot that produced the user's results.
old_sha="6da9ba96fc708de18a7b43d1db4ee818f8a8e7f4"
git -C "$repo" cat-file -e "$old_sha^{commit}"
git -C "$repo" archive --format=tar "$old_sha" | gzip -n \
  > "$bundle/experimental_cubic_previous_PR21/ed2-source.tar.gz"
printf '%s\n' "$old_sha" > "$bundle/experimental_cubic_previous_PR21/source-commit.txt"
git -C "$repo" show "$old_sha:tools/presto_cubic_exclusive.pbs" \
  > "$bundle/experimental_cubic_previous_PR21/presto_cubic_exclusive.pbs"
(cd "$bundle/experimental_cubic_previous_PR21" && sha256sum ed2-source.tar.gz > ed2-source.sha256)

git -C "$repo" fetch --no-tags origin refs/heads/main
main_sha=$(git -C "$repo" rev-parse 'FETCH_HEAD^{commit}')
git -C "$repo" archive --format=tar "$main_sha" | gzip -n \
  > "$bundle/latest_main_snapshot/ed2-main.tar.gz"
printf '%s\n' "$main_sha" > "$bundle/latest_main_snapshot/main-commit.txt"
(cd "$bundle/latest_main_snapshot" && sha256sum ed2-main.tar.gz > ed2-main.sha256)

comparison_ref="codex/ed2-cross-library-benchmark-protocol-20261119"
git -C "$repo" fetch --no-tags origin "refs/heads/$comparison_ref"
comparison_sha=$(git -C "$repo" rev-parse 'FETCH_HEAD^{commit}')
git -C "$repo" archive --format=tar "$comparison_sha" | gzip -n \
  > "$bundle/cross_library_comparator/comparator-source.tar.gz"
printf '%s\n' "$comparison_sha" > "$bundle/cross_library_comparator/comparator-commit.txt"
(cd "$bundle/cross_library_comparator" && sha256sum comparator-source.tar.gz > comparator-source.sha256)

python3 - "$bundle" "$source_sha" "$old_sha" "$main_sha" "$comparison_sha" <<'PY'
from pathlib import Path
import json,sys
root=Path(sys.argv[1])
source_sha,old_sha,main_sha,comparison_sha=sys.argv[2:]
assert all(len(s)==40 for s in (source_sha,old_sha,main_sha,comparison_sha))
data={
 "format":"ED2_CPC_AUTO_MNTE_PRESTO_20261119",
 "purpose":"manual source transfer only; no PRESTO job submitted",
 "candidate_commit":source_sha,
 "old_measured_PR21_commit":old_sha,
 "latest_remote_main_at_stage":main_sha,
 "cross_library_comparator_commit":comparison_sha,
 "candidate_contains_unmerged_cubic_features":True,
 "main_source_is_separate":True,
 "user_reported_presto_measurements_are_not_archived_raw_logs":True,
 "independent_operator_attestation":"PENDING",
 "citation_ready":False,
 "physical_model_sha256":"a774358f0c69f81982d3495707f45433acd223b3a40aacd109304fb371903fd7"
}
(root/"DELIVERY_MANIFEST.json").write_text(
  json.dumps(data,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
PY
(
  cd "$bundle"
  find . -type f ! -name DELIVERY_SHA256SUMS -print0 \
    | LC_ALL=C sort -z | xargs -0 sha256sum > DELIVERY_SHA256SUMS
  sha256sum -c DELIVERY_SHA256SUMS
  (cd cubic_auto_mnte && sha256sum -c ed2-source.sha256)
  (cd experimental_cubic_previous_PR21 && sha256sum -c ed2-source.sha256)
  (cd latest_main_snapshot && sha256sum -c ed2-main.sha256)
  (cd cross_library_comparator && sha256sum -c comparator-source.sha256)
  tar -tzf cubic_auto_mnte/ed2-source.tar.gz | grep -qx source/ham2vec.f90
  tar -tzf cubic_auto_mnte/ed2-source.tar.gz | grep -qx tools/cpc_presto_compare_mnte.py
  tar -tzf cubic_auto_mnte/ed2-source.tar.gz | grep -qx tests/cpc_cubic_auto_mnte_regression.py
)

# A moving main can invalidate the phrase 'latest-main snapshot'. Fail
# rather than silently delivering a version stale during creation.
fresh_main=$(git -C "$repo" ls-remote origin refs/heads/main | awk '{print $1}')
[[ "$main_sha" == "$fresh_main" ]] || {
  echo "Aborting: origin/main changed during source ZIP assembly" >&2; exit 3;
}
mkdir -p -- "$(dirname "$target")"
(cd "$tmp" && zip -q -r "$target" ED2_CPC_AUTO_MNTE_PRESTO_20261119)
unzip -tqq "$target"
echo "CPC_PRESTO_ZIP_OK: $target"
echo "CANDIDATE=$source_sha"
echo "OLD_PR21=$old_sha"
echo "MAIN=$main_sha"
echo "COMPARATOR=$comparison_sha"
sha256sum "$target"
