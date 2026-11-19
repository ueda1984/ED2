# tests/Makefile.am

# make check は軽量（ex1のみ）→ これは Automake ハーネスでOK
TESTS = test_ex1.sh
EXTRA_DIST = test_ex1.sh test_common.sh \
             test_ex2.sh test_ex3.sh test_ex4.sh test_ex5.sh \
             test_ex6.sh test_ex7.sh test_ex8.sh test_ex9.sh \
             cpc_csr_matvec_driver.f90 cpc_csr_main_regression.py \
             cpc_blas_lapack_contract_driver.f90 \
             cpc_csr_release_regression.py \
             cpc_partial_cache_allocation_regression.py \
             cpc_coeff_cache_bound_regression.py \
             cpc_csr_direct_build_regression.py \
             cpc_csr_direct_legacy_crossbuild.py

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
