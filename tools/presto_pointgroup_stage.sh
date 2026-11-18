#!/usr/bin/env bash
# Stage a clean, immutable ED2 source archive. Does not build, submit or
# modify PRESTO queues. Usage: bash tools/presto_pointgroup_stage.sh REPO STAGE
set -euo pipefail
if [[ $# -ne 2 ]]; then
  printf 'Usage: %s /path/to/ED2 /new/nonexistent/stage-dir\n' "$0" >&2
  exit 2
fi
repo=$(cd -- "$1" && pwd -P)
stage=$2
[[ ! -e "$stage" ]] || { echo "Refusing existing staging directory: $stage" >&2; exit 2; }
git -C "$repo" rev-parse --is-inside-work-tree | grep -qx true
if [[ -n $(git -C "$repo" status --porcelain --untracked-files=no) ]]; then
  echo "Tracked worktree is dirty. Commit and review the exact source first." >&2
  exit 2
fi
revision=$(git -C "$repo" rev-parse --verify HEAD)
[[ "$revision" =~ ^[0-9a-f]{40}$ ]] || exit 2
if [[ -v ED2_EXPECTED_COMMIT && "$revision" != "$ED2_EXPECTED_COMMIT" ]]; then
  echo "Source revision does not match ED2_EXPECTED_COMMIT" >&2
  exit 2
fi
mkdir -- "$stage"
stage=$(cd -- "$stage" && pwd -P)
printf '%s\n' "$revision" > "$stage/source-commit.txt"
git -C "$repo" archive --format=tar HEAD | gzip -n > "$stage/ed2-source.tar.gz"
(cd -- "$stage" && sha256sum ed2-source.tar.gz > ed2-source.sha256)
printf 'READY: commit %s\nSTAGE: %s\n' "$revision" "$stage"
printf 'VERIFY:\n  (cd "%s" && sha256sum -c ed2-source.sha256)\n' "$stage"
echo "No PBS jobs were submitted."
