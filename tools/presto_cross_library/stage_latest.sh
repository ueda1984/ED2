#!/usr/bin/env bash
# PRESTO cross-library: create immutable archives of the CURRENT remote main,
# a separate benchmark implementation ref, and this launch machinery.
# Usage: bash tools/presto_cross_library/stage_latest.sh /path/to/ED2 STAGE_NEW
# Never checks out branches, merges, pushes, modifies source, or submits PBS.
set -euo pipefail
umask 077
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 /path/to/ED2 /new/nonexistent/STAGE" >&2
  exit 2
fi
repo=$(cd -- "$1" && pwd -P)
dest=$2
[[ ! -e "$dest" ]] || { echo "Stage exists; choose a new unique path: $dest" >&2; exit 2; }
git -C "$repo" rev-parse --is-inside-work-tree | grep -qx true
remote=$(git -C "$repo" remote get-url origin)
case "$remote" in
  *github.com*/ueda1984/ED2*|*github.com:ueda1984/ED2*) ;;
  *) echo "Refusing unexpected origin: $remote" >&2; exit 2 ;;
esac
benchmark_ref=${ED2_BENCHMARK_REF:-codex/ed2-cross-library-benchmark-protocol-20261119}
launcher_ref=${ED2_LAUNCHER_REF:-codex/ed2-presto-cross-library-latest-main-20261119}
case "$benchmark_ref $launcher_ref" in
  *" "*)
    [[ "$benchmark_ref" =~ ^[A-Za-z0-9._/-]+$ && "$launcher_ref" =~ ^[A-Za-z0-9._/-]+$ ]] ||
      { echo "Invalid branch refs" >&2; exit 2; } ;;
esac

# Refresh remote main immediately before staging. SAVE each fetched SHA
# before the next git fetch changes FETCH_HEAD. Do not depend on stale HEAD.
git -C "$repo" fetch --no-tags origin refs/heads/main
main_sha=$(git -C "$repo" rev-parse "FETCH_HEAD^{commit}")
[[ "$main_sha" =~ ^[0-9a-f]{40}$ ]] || exit 2
if [[ -n ${ED2_EXPECTED_MAIN_SHA:-} && "$main_sha" != "$ED2_EXPECTED_MAIN_SHA" ]]; then
  echo "Main moved unexpectedly: $main_sha != $ED2_EXPECTED_MAIN_SHA" >&2
  exit 2
fi
git -C "$repo" fetch --no-tags origin "refs/heads/$benchmark_ref"
benchmark_sha=$(git -C "$repo" rev-parse "FETCH_HEAD^{commit}")
git -C "$repo" fetch --no-tags origin "refs/heads/$launcher_ref"
launcher_sha=$(git -C "$repo" rev-parse "FETCH_HEAD^{commit}")
for rev in "$main_sha" "$benchmark_sha" "$launcher_sha"; do
  [[ "$rev" =~ ^[0-9a-f]{40}$ ]] || exit 2
done
# Verify required benchmark/launch files exist in their own separate refs.
git -C "$repo" cat-file -e "$benchmark_sha:benchmarks/cross_library/run_quspin.py"
git -C "$repo" cat-file -e "$benchmark_sha:benchmarks/cross_library/export_cubic_xxz.py"
git -C "$repo" cat-file -e "$launcher_sha:tools/presto_cross_library/launch.pbs"
git -C "$repo" cat-file -e "$launcher_sha:tools/presto_cross_library/preflight.py"

mkdir -- "$dest"
stage=$(cd -- "$dest" && pwd -P)
git -C "$repo" archive --format=tar "$main_sha" | gzip -n > "$stage/ed2-main.tar.gz"
git -C "$repo" archive --format=tar "$benchmark_sha" benchmarks/cross_library |
  gzip -n > "$stage/benchmark-tools.tar.gz"
git -C "$repo" archive --format=tar "$launcher_sha" tools/presto_cross_library |
  gzip -n > "$stage/launch-tools.tar.gz"
printf '%s\n' "$main_sha" > "$stage/main-commit.txt"
printf '%s\n' "$benchmark_sha" > "$stage/benchmark-commit.txt"
printf '%s\n' "$launcher_sha" > "$stage/launcher-commit.txt"
printf '%s\n' "$benchmark_ref" > "$stage/benchmark-ref.txt"
printf '%s\n' "$launcher_ref" > "$stage/launcher-ref.txt"
(cd "$stage" && sha256sum ed2-main.tar.gz benchmark-tools.tar.gz launch-tools.tar.gz   > SHA256SUMS && sha256sum -c SHA256SUMS)
tar -xOzf "$stage/launch-tools.tar.gz"   tools/presto_cross_library/launch.pbs > "$stage/launch.pbs"
tar -xOzf "$stage/launch-tools.tar.gz"   tools/presto_cross_library/preflight.py > "$stage/preflight.py"
bash -n "$stage/launch.pbs"

# Explicit last-minute race detection. A new main merge invalidates only
# "latest main" labeling, not source integrity. Do not silently restage.
current_main=$(git -C "$repo" ls-remote origin refs/heads/main | awk '{print $1}')
if [[ "$current_main" != "$main_sha" ]]; then
  echo "Remote main advanced DURING staging ($main_sha -> $current_main)." >&2
  echo "Keep existing archive for provenance; make a NEW stage before submit." >&2
  exit 3
fi
python3 "$stage/preflight.py" --stage "$stage" --mode stage
echo "READY (no PBS submission)"
echo "STAGE=$stage"
echo "ED2_CURRENT_MAIN_COMMIT=$main_sha"
echo "BENCHMARK_COMMIT=$benchmark_sha"
echo "LAUNCHER_COMMIT=$launcher_sha"
echo "Important: run preflight --mode submit immediately before manual qsub."
