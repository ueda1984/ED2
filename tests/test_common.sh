#!/bin/sh
set -eu
# Common runner for the historical nine example entrypoints.
# Isolate writable data for both in-tree and out-of-tree builds.
exec "$(dirname "$0")/test_common_strict.sh" "$@"
