# tests/Makefile.am

# make check は軽量（ex1のみ）→ これは Automake ハーネスでOK
TESTS = test_ex1.sh
EXTRA_DIST = test_ex1.sh test_common.sh test_common_strict.sh compare_historical.py requirements-baseline.txt \
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
	  ED2_TEST_SOURCE_ROOT="$(abs_top_srcdir)" \
	  ED2_TEST_BUILD_ROOT="$(abs_top_builddir)" \
	  ED2_EXE="$(abs_top_builddir)/source/QS3ED2" \
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
  $(top_builddir)/source/krylov_solvers.o $(top_builddir)/source/krylov_observations.o $(top_builddir)/source/krylov_residual_policy.o \
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


# Numerical kernel tests require only the existing optional NumPy test dependency.
check_PROGRAMS += lanczos_stability_driver
lanczos_stability_driver_SOURCES = lanczos_stability_driver.f90
lanczos_stability_driver_LDADD = $(top_builddir)/source/krylov_solvers.o $(top_builddir)/source/krylov_observations.o $(top_builddir)/source/krylov_residual_policy.o \
  $(top_builddir)/source/blas_lapack_interfaces.o
EXTRA_DIST += lanczos_stability.py lanczos_cli_stability.py
.PHONY: check-lanczos
check-lanczos: lanczos_stability_driver sparse_field_operator_check
	$(PYTHON) $(srcdir)/lanczos_stability.py --driver "$(abs_builddir)/lanczos_stability_driver" \
	  $(if $(LANCZOS_OUTPUT),--output "$(LANCZOS_OUTPUT)",)

	$(PYTHON) $(srcdir)/lanczos_cli_stability.py --exe "$(abs_top_builddir)/source/QS3ED2" \
	  --operator-exe "$(abs_builddir)/sparse_field_operator_check" \
	  $(if $(LANCZOS_OUTPUT),--output "$(LANCZOS_OUTPUT)-cli",)

# Characterization is optional and records numerical findings without hiding them.
EXTRA_PROGRAMS = krylov_stress_driver
krylov_stress_driver_SOURCES = krylov_stress_driver.f90
krylov_stress_driver_LDADD = $(top_builddir)/source/krylov_solvers.o $(top_builddir)/source/krylov_observations.o $(top_builddir)/source/krylov_residual_policy.o \
  $(top_builddir)/source/blas_lapack_interfaces.o
EXTRA_DIST += krylov_stress.py krylov_benchmarks.py krylov_benchmark_probe.f90 krylov_acceptance_probe.py
.PHONY: check-krylov-stress
check-krylov-stress: krylov_stress_driver
	$(PYTHON) $(srcdir)/krylov_stress.py --driver "$(abs_builddir)/krylov_stress_driver" \
	  $(if $(KRYLOV_OUTPUT),--output "$(KRYLOV_OUTPUT)",)

# Phase 1-E: diagnostic/contract experiments, no production policy changes.
EXTRA_DIST += precision_contract.py
.PHONY: check-precision-contract
check-precision-contract: krylov_stress_driver
	$(PYTHON) $(srcdir)/precision_contract.py --driver "$(abs_builddir)/krylov_stress_driver" \
	  $(if $(PRECISION_OUTPUT),--output "$(PRECISION_OUTPUT)",)

# Additive internal telemetry/replay checks. A base-version driver is optional.
EXTRA_PROGRAMS += krylov_diagnostics_driver
krylov_diagnostics_driver_SOURCES = krylov_diagnostics_driver.F90
krylov_diagnostics_driver_FCFLAGS = $(AM_FCFLAGS) $(FCFLAGS) -DPHASE1F
krylov_diagnostics_driver_LDADD = $(top_builddir)/source/krylov_solvers.o \
  $(top_builddir)/source/krylov_observations.o $(top_builddir)/source/krylov_residual_policy.o $(top_builddir)/source/blas_lapack_interfaces.o
EXTRA_DIST += diagnostics_reproducibility.py
.PHONY: check-diagnostics
check-diagnostics: krylov_diagnostics_driver
	$(PYTHON) $(srcdir)/diagnostics_reproducibility.py --driver "$(abs_builddir)/krylov_diagnostics_driver" \
	  --exe "$(abs_top_builddir)/source/QS3ED2" \
	  $(if $(DIAGNOSTICS_OUTPUT),--output "$(DIAGNOSTICS_OUTPUT)",)

# Optional residual policy; default Automake TESTS and stopping rules stay intact.
EXTRA_PROGRAMS += krylov_policy_driver
krylov_policy_driver_SOURCES = krylov_policy_driver.F90
krylov_policy_driver_FCFLAGS = $(AM_FCFLAGS) $(FCFLAGS) -DPHASE1F -DPHASE1G
krylov_policy_driver_LDADD = $(top_builddir)/source/krylov_solvers.o \
  $(top_builddir)/source/krylov_observations.o $(top_builddir)/source/krylov_residual_policy.o \
  $(top_builddir)/source/blas_lapack_interfaces.o
EXTRA_DIST += residual_policy_validation.py
.PHONY: check-residual-policy
check-residual-policy: krylov_policy_driver
	$(PYTHON) $(srcdir)/residual_policy_validation.py --driver "$(abs_builddir)/krylov_policy_driver" \
	  --exe "$(abs_top_builddir)/source/QS3ED2" \
	  $(if $(RESIDUAL_OUTPUT),--output "$(RESIDUAL_OUTPUT)",)

# Source tree can be read-only and distinct from the build tree in distcheck.
AM_TESTS_ENVIRONMENT = ED2_TEST_SOURCE_ROOT='$(abs_top_srcdir)'; export ED2_TEST_SOURCE_ROOT; \
                       ED2_TEST_BUILD_ROOT='$(abs_top_builddir)'; export ED2_TEST_BUILD_ROOT; \
                       FC='$(FC)'; export FC;
