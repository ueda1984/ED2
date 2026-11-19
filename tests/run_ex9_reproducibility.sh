#!/bin/sh
# Full ex9; run only inside a separately approved compute-node allocation.
# This wrapper does not select resources or submit a scheduler job.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out="$root/examples/cubic_sp_HB"
if [ ! -x "$root/source/QS3ED2" ] || [ ! -x "$root/script_cubic_sp_HB/mk_input" ]; then
  echo "ERROR: build ED2 and mk_input first (make)." >&2
  exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: Python 3 required for verification." >&2
  exit 2
fi
# Original example scripts append observables; never mix old and new rows.
for file in eigenvalues.dat local_mag.dat two_body_cf_xyz.dat two_body_cf_z+-.dat; do
  if [ -e "$out/$file" ]; then
    echo "ERROR: existing $out/$file; use a fresh checkout or scratch copy." >&2
    exit 2
  fi
done
: "${OMP_NUM_THREADS:=1}"
: "${OPENBLAS_NUM_THREADS:=1}"
: "${MKL_NUM_THREADS:=1}"
export OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS
export OMP_DYNAMIC=FALSE
audit="${ED2_EX9_AUDIT_DIR:-$root/ex9-audit-$(date -u +%Y%m%dT%H%M%SZ)}"
if [ -e "$audit" ]; then
  echo "ERROR: audit directory already exists: $audit" >&2
  exit 2
fi
mkdir -p "$audit"
{
  date -u
  git -C "$root" rev-parse HEAD 2>/dev/null || true
  git -C "$root" status --short 2>/dev/null || true
  if [ -f "$root/Makefile" ]; then
    echo 'Configured build variables (not the default compiler on PATH):'
    sed -n '/^FC = /p; /^FCFLAGS = /p; /^LDFLAGS = /p; /^LIBS = /p' "$root/Makefile"
  fi
  printf 'OMP_NUM_THREADS=%s\nOPENBLAS_NUM_THREADS=%s\nMKL_NUM_THREADS=%s\n' \
    "$OMP_NUM_THREADS" "$OPENBLAS_NUM_THREADS" "$MKL_NUM_THREADS"
} > "$audit/environment.txt"
if [ -f "$root/config.log" ]; then
  cp "$root/config.log" "$audit/build-config.log"
fi
sha256sum "$root/source/QS3ED2" "$root/script_cubic_sp_HB/mk_input" \
  > "$audit/executables.sha256"
if command -v ldd >/dev/null 2>&1; then
  ldd "$root/source/QS3ED2" > "$audit/linkage.txt" 2>&1 || true
fi

echo "Running full ex9; log: $audit/run.log"
if ! (cd "$root/script_cubic_sp_HB" && TOPDIR="$root" \
        ED2_EXE="$root/source/QS3ED2" ./run.sh) > "$audit/run.log" 2>&1; then
  echo "ERROR: ex9 failed/interrupted; no reproduction claim is supported." >&2
  exit 1
fi
if ! python3 "$root/tests/verify_ex9.py" \
      --reference "$root/examples/ref_data/cubic_sp_HB" \
      --actual "$out" --legacy-ex9-symmetry-projection \
      --json-out "$audit/comparison.json" \
      > "$audit/verification.log" 2>&1; then
  echo "ERROR: ex9 outputs failed verification; see $audit/comparison.json." >&2
  exit 1
fi
for file in input.dat input_param.dat list_xyz_dm_gamma.dat list_hvec.dat \
    list_spin.dat list_NODmax.dat list_ij_cf.dat \
    list_p1.dat list_p2.dat list_p3.dat list_p4.dat list_p5.dat list_p6.dat \
    eigenvalues.dat local_mag.dat two_body_cf_xyz.dat two_body_cf_z+-.dat; do
  if [ -f "$out/$file" ]; then
    cp "$out/$file" "$audit/$file"
  fi
done
( cd "$audit" && sha256sum ./*.dat ./*.json ./*.txt ./*.log ./*.sha256 > SHA256SUMS )
echo "PASS: full ex9 completed and verified. Results: $audit"
