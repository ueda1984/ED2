#!/bin/sh
set -eu
# Usage: test_common.sh script_chain examples/chain/output.dat ...
# All numerical output goes into a writable sandbox. Build binaries come
# from the build tree, while scripts and model templates come from source.

model_dir=$1
shift
source_root=${ED2_TEST_SOURCE_ROOT:-$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)}
build_root=${ED2_TEST_BUILD_ROOT:-$source_root}
case "$source_root" in
  /*) ;;
  *) source_root=$(CDPATH= cd -- "$source_root" && pwd) ;;
esac
case "$build_root" in
  /*) ;;
  *) build_root=$(CDPATH= cd -- "$build_root" && pwd) ;;
esac

exe=${ED2_EXE:-$build_root/source/QS3ED2}
if [ ! -x "$exe" ]; then
  echo "ERROR: missing ED2 executable: $exe" >&2
  exit 1
fi
generator="$build_root/$model_dir/mk_input"
if [ ! -x "$generator" ]; then
  echo "ERROR: missing model generator: $generator" >&2
  exit 1
fi
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export FC=${FC:-gfortran}

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT HUP INT TERM
sandbox="$scratch/ED2"
mkdir -p "$sandbox"
cp -R "$source_root/$model_dir" "$sandbox/"
cp -R "$source_root/input" "$sandbox/"
# Automake distcheck intentionally makes the unpacked source tree read-only.
# cp -R propagates those modes, so normalize permissions in the private copy.
chmod -R u+w "$sandbox"
cp "$generator" "$sandbox/$model_dir/mk_input"
(
  cd "$sandbox/$model_dir"
  TOPDIR="$sandbox" ED2_EXE="$exe" ./run.sh > test.log 2>&1
) || {
  echo "ERROR: $model_dir/run.sh failed (isolated test):" >&2
  tail -n 200 "$sandbox/$model_dir/test.log" >&2 || true
  exit 1
}

for output in "$@"; do
  if [ ! -s "$sandbox/$output" ]; then
    echo "ERROR: required output missing or empty: $output" >&2
    tail -n 200 "$sandbox/$model_dir/test.log" >&2 || true
    exit 1
  fi
done

if [ "${LONGTEST:-0}" = 1 ]; then
  # Preserve the historical extended-test eigenvalue-token sanity check.
  name=${model_dir#script_}
  eigen="$sandbox/examples/$name/eigenvalues.dat"
  if [ -s "$eigen" ]; then
    e0=$(awk 'NR==1 { print $1; exit }' "$eigen")
    case "$e0" in
      ""|*[!0-9eE+.-]*) echo "ERROR: nonnumeric first eigenvalue token: '$e0'" >&2; exit 1 ;;
    esac
  fi
fi
# Optional capture of the exact generated outputs for independent historical
# comparisons. The default sandbox remains ephemeral and self-cleaning.
if [ -n "${ED2_TEST_CAPTURE_ROOT:-}" ]; then
  name=${model_dir#script_}
  mkdir -p "$ED2_TEST_CAPTURE_ROOT"
  if [ -e "$ED2_TEST_CAPTURE_ROOT/$name" ]; then
    echo "ERROR: capture destination exists: $ED2_TEST_CAPTURE_ROOT/$name" >&2
    exit 1
  fi
  cp -R "$sandbox/examples/$name" "$ED2_TEST_CAPTURE_ROOT/$name"
fi

echo "OK: $model_dir passed in isolated sandbox."
