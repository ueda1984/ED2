! Test-only native cache driver for full-column old/new matrix element proofs.
program cubic_cache_driver
  use iso_fortran_env, only: int64, real64
  use input_param, only: read_ip,NOS,NODmax,NODmin,local_NODmax,THS, &
    list_s,list_r,explist
  use state_lists, only: get_cum_ltd_rep_combination,cum_ltd_rep_combination,allocate_lists_omp_SQ
  use ham2vec, only: make_wk_loc_and_ele,ham_to_vec_wave_vector
  implicit none
  complex(real64), allocatable :: input_vector(:),output_vector(:)
  character(len=16) :: dump_setting
  integer :: status,u,ios,j
  integer(int64) :: re_bits,im_bits
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
  ! Test-only SHA256 oracle for every output component of a nonreal,
  ! deterministic trial vector. Default cache tests remain unchanged.
  dump_setting=''
  call get_environment_variable('ED2_CUBIC_MATVEC_PROOF_DUMP',dump_setting,status=status)
  if(status==0.and.trim(adjustl(dump_setting))=='1')then
    allocate(input_vector(THS),output_vector(THS))
    do j=1,THS
      input_vector(j)=cmplx(real(modulo(j*31,103)-51,real64)/103.0_real64, &
                           real(modulo(j*43,97)-48,real64)/97.0_real64,real64)
    end do
    call ham_to_vec_wave_vector(output_vector,input_vector,THS,NODmax,NODmin,list_s,list_r,explist)
    open(newunit=u,file='cubic_matvec_proof_dump.txt',status='replace',action='write',iostat=ios)
    if(ios/=0)error stop 'Failed to open matvec proof dump'
    do j=1,THS
      re_bits=transfer(real(output_vector(j),real64),0_int64)
      im_bits=transfer(aimag(output_vector(j)),0_int64)
      write(u,'(i0,1x,z16.16,1x,z16.16)',iostat=ios) j,re_bits,im_bits
      if(ios/=0)error stop 'Failed to write matvec proof dump'
    end do
    close(u,iostat=ios)
    if(ios/=0)error stop 'Failed to close matvec proof dump'
    write(*,'(a,i0)') 'ED2_CUBIC_MATVEC_PROOF_ROWS=',THS
  end if
end program cubic_cache_driver
