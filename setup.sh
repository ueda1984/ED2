#!/bin/sh
# Preserve the historical entry point only to explain its replacement.
cat >&2 <<'EOF'
setup.sh has been retired because its historical build/example workflow is obsolete.

Follow the current installation guide:
https://ed2.readthedocs.io/en/latest/getting-started/installation/

After installing the documented GNU Fortran, Autotools and LP64 BLAS/LAPACK
prerequisites, use a fresh checkout and run from its repository root:

  autoreconf -vfi
  ./configure FC=gfortran FCFLAGS="-O2 -fopenmp"
  make -j2

Then set thread counts and run the smoke test:

  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  make check

Use the current quickstart to run individual examples:
https://ed2.readthedocs.io/en/latest/getting-started/quickstart/
EOF
exit 2
