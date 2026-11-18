! Test-only driver for ED2 symmetry-basis construction, no Hamiltonian/Lanczos.
program cubic_basis_driver
  use input_param, only: read_ip,NOS,NODmax,NODmin,local_NODmax,THS
  use state_lists, only: get_cum_ltd_rep_combination,cum_ltd_rep_combination,allocate_lists_omp_SQ
  implicit none
  call read_ip
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  if(NODmax>0)then
    THS=cum_ltd_rep_combination(NOS,NODmax)
  else
    THS=1
  end if
  call allocate_lists_omp_SQ
  write(*,'(a,i0)') 'ED2_BASIS_DRIVER_DIM=',THS
end program cubic_basis_driver
