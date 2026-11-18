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

.PHONY: check-sparse-field
check-sparse-field:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-sparse-field SPARSE_FIELD_OUTPUT="$(SPARSE_FIELD_OUTPUT)"

.PHONY: check-lanczos
check-lanczos:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-lanczos LANCZOS_OUTPUT="$(LANCZOS_OUTPUT)"

.PHONY: check-krylov-stress
check-krylov-stress:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-krylov-stress KRYLOV_OUTPUT="$(KRYLOV_OUTPUT)"

.PHONY: check-precision-contract
check-precision-contract:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-precision-contract PRECISION_OUTPUT="$(PRECISION_OUTPUT)"

.PHONY: check-diagnostics
check-diagnostics:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-diagnostics DIAGNOSTICS_OUTPUT="$(DIAGNOSTICS_OUTPUT)"

.PHONY: check-residual-policy
check-residual-policy:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-residual-policy RESIDUAL_OUTPUT="$(RESIDUAL_OUTPUT)"

# Opt-in numerical regression: restored historical QS-Cube two-stage ALG=1.
.PHONY: check-qscube-lanczos
check-qscube-lanczos:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-qscube-lanczos QSCUBE_OUTPUT="$(QSCUBE_OUTPUT)"

# Optional independent physical test of combined translations and point-group irreps.
.PHONY: check-pointgroup
check-pointgroup:
	$(AM_V_at)$(MAKE) -C source all
	$(AM_V_at)$(MAKE) -C tests check-pointgroup POINTGROUP_OUTPUT="$(POINTGROUP_OUTPUT)"

# Include reference fixtures, example templates and user documentation in
# source archives used outside the Git repository.
EXTRA_DIST = input/input.dat examples/ref_data docs tools CITATION.cff LICENSE mkdocs.yml .readthedocs.yaml
