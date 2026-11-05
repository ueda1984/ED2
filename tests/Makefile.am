# tests/Makefile.am

# make check は軽量（ex1のみ）→ これは Automake ハーネスでOK
TESTS = test_ex1.sh
EXTRA_DIST = test_ex1.sh test_common.sh \
             test_ex2.sh test_ex3.sh test_ex4.sh test_ex5.sh \
             test_ex6.sh test_ex7.sh test_ex8.sh test_ex9.sh

# check-long は“任意”なので自前ループで堅牢に実行
LONG_TESTS = \
  test_ex1.sh \
  test_ex2.sh \
  test_ex3.sh \
  test_ex4.sh \
  test_ex5.sh \
  test_ex6.sh \
  test_ex7.sh \
  test_ex8.sh \
  test_ex9.sh

.PHONY: check-long
check-long:
	@set -e; \
	echo "== Running long tests: $(LONG_TESTS) =="; \
	for t in $(LONG_TESTS); do \
	  echo "---- $$t ----"; \
	  rm -f "$$t.long.log"; \
	  TOPDIR="$(top_srcdir)" \
	  ED2_EXE="$(top_builddir)/source/QS3ED2" \
	  LONGTEST=1 $(SHELL) "$(srcdir)/$$t" > "$$t.long.log" 2>&1 || { \
	    echo "FAIL: $$t"; \
	    echo "===== last 200 lines of $$t.long.log ====="; \
	    tail -n 200 "$$t.long.log" || true; \
	    exit 1; \
	  }; \
	  echo "PASS: $$t"; \
	done; \
	echo "OK: check-long finished."

# Separate target preserves the historical one-test Automake baseline.
check_PROGRAMS = lapack_interface_check lapack_failure_check
AM_FCFLAGS = -I$(abs_top_builddir)/source
lapack_interface_check_SOURCES = lapack_interface_check.f90
lapack_interface_check_LDADD = $(top_builddir)/source/blas_lapack_interfaces.o
lapack_failure_check_SOURCES = lapack_failure_check.f90
lapack_failure_check_LDADD = $(top_builddir)/source/eigen_solver.o \
  $(top_builddir)/source/ham2vec.o $(top_builddir)/source/state_lists.o \
  $(top_builddir)/source/input_param.o $(top_builddir)/source/blas_lapack_interfaces.o
EXTRA_DIST += test_lapack_interfaces.sh
.PHONY: check-lapack
check-lapack: $(check_PROGRAMS)
	$(SHELL) $(srcdir)/test_lapack_interfaces.sh ./lapack_interface_check ./lapack_failure_check

# Optional numerical validation; Python/NumPy are not runtime/build dependencies.
check_PROGRAMS += sparse_field_operator_check
sparse_field_operator_check_SOURCES = sparse_field_operator_check.f90
sparse_field_operator_check_LDADD = $(top_builddir)/source/ham2vec.o \
  $(top_builddir)/source/state_lists.o $(top_builddir)/source/input_param.o
EXTRA_DIST += sparse_field_validation.py numerical_baseline.py
PYTHON ?= python3
.PHONY: check-sparse-field
check-sparse-field: sparse_field_operator_check
	$(PYTHON) $(srcdir)/sparse_field_validation.py \
	  --exe "$(abs_top_builddir)/source/QS3ED2" --operator-exe "$(abs_builddir)/sparse_field_operator_check" \
	  $(if $(SPARSE_FIELD_OUTPUT),--output "$(SPARSE_FIELD_OUTPUT)",)
