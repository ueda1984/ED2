bin_PROGRAMS = QS3ED2
QS3ED2_SOURCES = blas_lapack_interfaces.f90 \
              input_param.f90 \
              state_lists.f90 \
              ham2vec.f90 \
              lanczos.f90 \
              eigen_solver.f90 \
              get_expectation_values.f90 \
              main.f90

# ------------------------------------------------------------
# Fortran module dependency order (required for parallel build)
# ------------------------------------------------------------
# input_param.f90 -> input_param.mod
# state_lists.f90 uses input_param
state_lists.o: input_param.o

# ham2vec uses state_lists (+ likely input_param indirectly)
ham2vec.o: input_param.o state_lists.o

# lanczos uses ham2vec and input_param
lanczos.o: input_param.o ham2vec.o blas_lapack_interfaces.o

# eigen_solver uses ham2vec (and possibly input_param)
eigen_solver.o: input_param.o ham2vec.o blas_lapack_interfaces.o

# get_expectation_values uses input_param and ham2vec
get_expectation_values.o: input_param.o ham2vec.o

# main uses state_lists and input_param (+ calls others)
main.o: input_param.o state_lists.o ham2vec.o lanczos.o eigen_solver.o get_expectation_values.o
