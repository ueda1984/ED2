! Test-only operator exporter; QS3ED2 and its algorithms are not instrumented.
program sparse_field_operator_check
  use input_param, only: read_ip, NOS, NODmax, NODmin, local_NODmax, THS, &
    list_s, list_r, explist, wk_dim
  use state_lists, only: get_cum_ltd_rep_combination, cum_ltd_rep_combination, &
    allocate_lists_omp_SQ, list_fly
  use ham2vec, only: make_full_hamiltonian, make_wk_loc_and_ele, ham_to_vec_wave_vector, wk_loc, wk_ele
  implicit none
  complex(8), allocatable :: h(:,:), x(:), y(:)
  integer, allocatable :: occupation(:)
  integer :: i, j, mode
  integer :: caches(3)=[0,3,-1]
  character(50) :: filename
  call read_ip
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  THS=cum_ltd_rep_combination(NOS,NODmax)
  call allocate_lists_omp_SQ
  allocate(h(THS,THS),x(THS),y(THS),occupation(NODmax))
  open(78,file='basis.dat',status='new')
  do i=1,THS
    occupation=list_fly(list_s(i),NODmax,NODmin,NOS)
    write(78,'(100i8)') i,occupation
  end do
  close(78)
  h=(0d0,0d0)
  call make_full_hamiltonian(THS,h)
  open(79,file='dense.dat',status='new')
  do j=1,THS
    do i=1,THS
      write(79,'(2i8,2es25.15)') i,j,real(h(i,j),8),aimag(h(i,j))
    end do
  end do
  close(79)
  do mode=1,size(caches)
    wk_dim=caches(mode)
    call make_wk_loc_and_ele(THS)
    write(filename,'(a,i0,a)') 'apply-',caches(mode),'.dat'
    open(80,file=trim(filename),status='new')
    do j=1,THS
      x=(0d0,0d0)
      x(j)=(1d0,0d0)
      call ham_to_vec_wave_vector(y,x,THS,NODmax,NODmin,list_s,list_r,explist)
      do i=1,THS
        write(80,'(2i8,2es25.15)') i,j,real(y(i),8),aimag(y(i))
      end do
    end do
    close(80)
    deallocate(wk_loc,wk_ele)
  end do
end program sparse_field_operator_check
