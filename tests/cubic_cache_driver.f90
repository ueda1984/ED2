! Test-only native cache driver for full-column old/new matrix element proofs.
program cubic_cache_driver
  use input_param, only: read_ip,NOS,NODmax,NODmin,local_NODmax,THS
  use state_lists, only: get_cum_ltd_rep_combination,cum_ltd_rep_combination,allocate_lists_omp_SQ
  use ham2vec, only: make_wk_loc_and_ele
  implicit none
  call read_ip
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  if(NODmax>0)then
    THS=cum_ltd_rep_combination(NOS,NODmax)
  else
    THS=1
  end if
  call allocate_lists_omp_SQ
  call make_wk_loc_and_ele(THS)
  write(*,'(a,i0)') 'ED2_CACHE_DRIVER_DIM=',THS
end program cubic_cache_driver
