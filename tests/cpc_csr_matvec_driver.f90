! Independent test-only H*x driver built against the original ED2 modules.
! Does not depend on cubic symmetry optimizations or newer Lanczos features.
program cpc_csr_matvec_driver
  use iso_fortran_env, only: real64, int64
  use input_param, only: read_ip,NOS,NODmax,NODmin,local_NODmax,THS, &
    list_s,list_r,explist
  use state_lists, only: get_cum_ltd_rep_combination, &
    cum_ltd_rep_combination,allocate_lists_omp_SQ
  use ham2vec, only: make_wk_loc_and_ele,ham_to_vec_wave_vector, &
    make_full_hamiltonian,wk_ele,wk_loc,csr_coeff,csr_release_padded
  implicit none
  complex(real64), allocatable :: v(:),h_v(:),dense_h(:,:)
  real(real64) :: delta
  integer :: i,t,unit_no,io_status
  integer(int64) :: bits_re,bits_im

  call read_ip
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  if(NODmax>0)then
    THS=cum_ltd_rep_combination(NOS,NODmax)
  else
    THS=1
  end if
  if(THS<=0)error stop 'CSR test: invalid basis dimension'
  call allocate_lists_omp_SQ
  call make_wk_loc_and_ele(THS)
  if(allocated(wk_ele))then
    write(*,'(a,i0)') 'CPC_CSR_NATIVE_IMAG_COUNT=',count(aimag(wk_ele)/=0.0_real64)
    if(.not.allocated(wk_loc))error stop 'CSR test: partial padded cache retained'
    write(*,'(a)') 'CPC_CSR_PADDED_RETAINED_PASS'
  else
    if(.not.csr_release_padded)error stop 'CSR test: unexpected padded release'
    if(allocated(wk_loc))error stop 'CSR test: padded location retained'
    write(*,'(a,i0)') 'CPC_CSR_NATIVE_IMAG_COUNT=',count(aimag(csr_coeff)/=0.0_real64)
    write(*,'(a)') 'CPC_CSR_PADDED_RELEASE_PASS'
  end if

  allocate(v(THS),h_v(THS))
  if(THS<=216)then
    allocate(dense_h(THS,THS))
    dense_h=(0.0_real64,0.0_real64)
    call make_full_hamiltonian(THS,dense_h)
  end if
  open(newunit=unit_no,file='csr_matvec_bits.txt',status='replace', &
       action='write',iostat=io_status)
  if(io_status/=0)error stop 'CSR test: cannot open output'
  do t=1,2
    do i=1,THS
      if(t==1)then
        v(i)=cmplx(real(modulo(i*31,103)-51,real64)/103.0_real64, &
                   real(modulo(i*43,97)-48,real64)/97.0_real64,real64)
      else
        v(i)=cmplx(real(modulo(i*37,113)-56,real64)/113.0_real64, &
                   real(modulo(i*23,89)-44,real64)/89.0_real64,real64)
      end if
    end do
    call ham_to_vec_wave_vector(h_v,v,THS,NODmax,NODmin,list_s,list_r,explist)
    if(allocated(dense_h))then
      delta=maxval(abs(matmul(dense_h,v)-h_v))
      if(.not.(delta<=1.0e-10_real64))then
        write(*,*) 'CPC_CSR_DENSE_MISMATCH =',delta
        error stop 'CSR test: independent dense Hamiltonian mismatch'
      end if
      write(*,'(a,i0,a,es14.6)') &
        'CPC_CSR_DENSE_REFERENCE_PASS vector=',t,', max_delta=',delta
    end if
    do i=1,THS
      bits_re=transfer(real(h_v(i),real64),0_int64)
      bits_im=transfer(aimag(h_v(i)),0_int64)
      write(unit_no,'(i0,1x,i0,1x,z16.16,1x,z16.16)',iostat=io_status) &
        t,i,bits_re,bits_im
      if(io_status/=0)error stop 'CSR test: cannot write output'
    end do
  end do
  close(unit_no,iostat=io_status)
  if(io_status/=0)error stop 'CSR test: cannot close output'
  write(*,'(a,i0)') 'CPC_CSR_H_VECTOR_DIM=',THS
  write(*,'(a,i0)') 'CPC_CSR_H_VECTOR_COUNT=',2
end program cpc_csr_matvec_driver
