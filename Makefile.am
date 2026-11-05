SUBDIRS = \
  source \
  tests \
  script_chain \
  script_cubic \
  script_cubic_sp_HB \
  script_honeycomb \
  script_j1j2_chain \
  script_kagome \
  script_mixed_spin_chain \
  script_square \
  script_triangular

.PHONY: check-long
check-long:
	$(AM_V_at)$(MAKE) -C tests check-long

.PHONY: check-lapack
check-lapack:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-lapack
