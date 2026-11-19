#!/usr/bin/env bash
# Produce the user's standalone PRESTO delivery ZIP on an authorized machine.
# It includes immutable current remote main, cross-library tools, launcher,
# and a CLEARLY SEPARATE unmerged experimental cubic source.
# NO qsub, no PRESTO connection, no merge, and no modifications to working tree.
set -euo pipefail
umask 077
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 /path/to/ueda1984-ED2 /absolute/output.zip" >&2; exit 2
fi
repo=$(cd -- "$1" && pwd -P)
target=$2
[[ "$target" == /* ]] || { echo "Output ZIP must be an absolute path" >&2; exit 2; }
[[ ! -e "$target" ]] || { echo "Refusing overwrite: $target" >&2; exit 2; }
root=$(mktemp -d)
trap 'rm -rf -- "$root"' EXIT
bundle="$root/ED2_PRESTO_DELIVERY_20261119"
mkdir -- "$bundle"
bash "$repo/tools/presto_cross_library/stage_latest.sh" "$repo"   "$bundle/latest_main_stage" > "$bundle/STAGING_LOG.txt"
# Capture optional unmerged special-purpose source with its own complete
# provenance. This source MUST NOT be labelled as 'latest ED2 main'.
special_ref="codex/ed2-cubic-direct-hop-elide-unused-hash"
git -C "$repo" fetch --no-tags origin "refs/heads/$special_ref"
special_sha=$(git -C "$repo" rev-parse "FETCH_HEAD^{commit}")
[[ "$special_sha" =~ ^[0-9a-f]{40}$ ]] || exit 2
git -C "$repo" cat-file -e "$special_sha:tools/presto_cubic_exclusive.pbs"
special="$bundle/experimental_cubic_PR21"
mkdir -- "$special"
git -C "$repo" archive --format=tar "$special_sha" | gzip -n   > "$special/ed2-source.tar.gz"
printf '%s\n' "$special_sha" > "$special/source-commit.txt"
printf '%s\n' "$special_ref" > "$special/source-ref.txt"
(cd "$special" && sha256sum ed2-source.tar.gz > ed2-source.sha256)
git -C "$repo" show "$special_sha:tools/presto_cubic_exclusive.pbs"   > "$special/presto_cubic_exclusive.pbs"
bash -n "$special/presto_cubic_exclusive.pbs"

launcher_sha=$(cat "$bundle/latest_main_stage/launcher-commit.txt")
git -C "$repo" show "$launcher_sha:tools/presto_cross_library/ZIP_DEPLOY_JP.md"   > "$bundle/README_PRESTO_JP.md"

# Distinguish four source commits in a machine-readable delivery manifest.
python3 - "$bundle" "$special_sha" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1]); special_sha=sys.argv[2]
stage=root/"latest_main_stage"
shas={k:(stage/(k+"-commit.txt")).read_text().strip()
      for k in ("main","benchmark","launcher")}
shas["experimental_cubic_pr21"]=special_sha
physical_hash="a774358f0c69f81982d3495707f45433acd223b3a40aacd109304fb371903fd7"
data={"package":"ED2 PRESTO user-transfer bundle",
      "purpose":"manual PRESTO deployment, no jobs submitted",
      "main_is_actual_remote_ref_at_staging":"origin/main",
      "pinned_commits":shas,
      "model":"N=1000 m=3 simple cubic J1-J2 XXZ Gamma/Oh A1g D=4299",
      "physical_model_sha256":physical_hash,
      "experimental_branch_is_NOT_main":True,
      "publication_ready":False,
      "must_verify_remote_main_before_1000_A1g_run":True}
(root/"DELIVERY_MANIFEST.json").write_text(
    json.dumps(data,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
PY
(
  cd "$bundle"
  sha256sum README_PRESTO_JP.md STAGING_LOG.txt DELIVERY_MANIFEST.json     latest_main_stage/* experimental_cubic_PR21/* > DELIVERY_SHA256SUMS
  sha256sum -c DELIVERY_SHA256SUMS
  (cd latest_main_stage && sha256sum -c SHA256SUMS)
  (cd experimental_cubic_PR21 && sha256sum -c ed2-source.sha256)
  python3 latest_main_stage/preflight.py     --stage latest_main_stage --mode stage > /dev/null
)
# Verify the source main did not change during the long archive build.
main_now=$(git -C "$repo" ls-remote origin refs/heads/main | awk '{print $1}')
main_staged=$(cat "$bundle/latest_main_stage/main-commit.txt")
if [[ "$main_now" != "$main_staged" ]]; then
  echo "ABORT delivery: main advanced while building ZIP." >&2
  exit 3
fi
mkdir -p -- "$(dirname "$target")"
(cd "$root" && zip -q -r "$target" ED2_PRESTO_DELIVERY_20261119)
unzip -tqq "$target"
printf 'CREATED %s\n' "$target"
printf 'ED2 MAIN %s\n' "$main_staged"
printf 'CUBIC PR21 (EXPERIMENTAL) %s\n' "$special_sha"
sha256sum "$target"
