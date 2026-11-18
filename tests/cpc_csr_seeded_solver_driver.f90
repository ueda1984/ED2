! Test-only reproduction of ORIGINAL ED2 thick-restarted Lanczos with fixed seed.
! This program does not alter any production RNG or eigensolver behavior.
program cpc_csr_seeded_solver_driver
  use input_param, only: read_ip,NOS,NODmax,NODmin,local_NODmax,THS, &
    write_index_ene
  use state_lists, only: get_cum_ltd_rep_combination, &
    cum_ltd_rep_combination,allocate_lists_omp_SQ
  use ham2vec, only: make_wk_loc_and_ele
  use eigen_solver, only: my_trlanczos_routines_complex,NOE
  implicit none
  integer :: nseed
  integer, allocatable :: seed(:)
  real(8), allocatable :: energies(:)
  complex(8), allocatable :: vectors(:,:)

  call random_seed(size=nseed)
  allocate(seed(nseed))
  seed=1984
  call random_seed(put=seed)
  call read_ip
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  if(NODmax>0)then
    THS=cum_ltd_rep_combination(NOS,NODmax)
  else
    THS=1
  end if
  if(THS<=0)error stop 'CSR test: invalid Lanczos dimension'
  call allocate_lists_omp_SQ
  call make_wk_loc_and_ele(THS)
  call my_trlanczos_routines_complex(energies,vectors,THS)
  call write_index_ene(energies,NOE)
  write(*,'(a)') 'CPC_CSR_OLD_TRLAN_FIXED_SEED=1984'
end program cpc_csr_seeded_solver_driver
