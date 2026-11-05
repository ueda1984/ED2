#!/bin/sh
set -eu
check=$1
failure=$2
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
"$check"
count=1
for mode in info nan inf negative overflow small integer dimension; do
  if "$check" "$mode" > "$work/$mode.log" 2>&1; then
    echo "FAIL: $mode was accepted" >&2
    exit 1
  fi
  case "$mode" in
    info) expected='injected query: LAPACK INFO = -18' ;;
    *) expected='ZHEEVR:' ;;
  esac
  grep -F "$expected" "$work/$mode.log" > /dev/null
  count=$((count+1))
done
# No numerical output may be written after either injected failure.
case "$failure" in /*) ;; *) failure="$(pwd)/$failure" ;; esac
for mode in query solve; do
  mkdir "$work/$mode"
  if (cd "$work/$mode" && "$failure" "$mode") > "$work/$mode.log" 2>&1; then
    echo "FAIL: dense $mode ignored LAPACK INFO" >&2
    exit 1
  fi
  case "$mode" in
    query) expected='ZHEEVR (dense workspace query): LAPACK INFO = -18' ;;
    solve) expected='ZHEEVR (dense solve): LAPACK INFO = 2' ;;
  esac
  grep -F "$expected" "$work/$mode.log" > /dev/null
  test ! -e "$work/$mode/eigenvalues.dat"
  count=$((count+1))
done
echo "PASS: $count interface/error-path checks"
