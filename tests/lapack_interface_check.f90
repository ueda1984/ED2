program lapack_interface_check
  use blas_lapack_interfaces
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  implicit none
  complex(real64) :: a(2,2), z(2,2), wq(1)
  complex(real64), allocatable :: work(:)
  real(real64) :: w(2), rq(1), d(2), e(1), tz(2,2), tw(10)
  real(real64) :: query
  integer :: support(4), iq(1), info, m, lw, lrw, liw, tiw(10), fail(2)
  real(real64), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  character(40) :: mode
  call get_command_argument(1, mode)
  query = 4.0_real64
  select case(trim(mode))
  case('info')
    call check_lapack_info('injected query', -18)
    stop 99
  case('nan')
    query = ieee_value(query,ieee_quiet_nan)
  case('inf')
    query = ieee_value(query,ieee_positive_inf)
  case('negative')
    query = -1.0_real64
  case('overflow')
    query = real(huge(lw),real64)+1.0_real64
  case('small')
    query = 3.0_real64
  case('integer')
    call zheevr_workspace_lengths(2, (4d0,0d0), 48d0, 19, lw, lrw, liw)
    stop 99
  case('dimension')
    call zheevr_workspace_lengths(huge(lw), (4d0,0d0), 48d0, 20, lw, lrw, liw)
    stop 99
  case('')
    ! A valid fractional recommendation must round up, never underallocate.
    call zheevr_workspace_lengths(2, (4.5d0,0d0), 48.5d0, 20, lw, lrw, liw)
    if(lw/=5 .or. lrw/=49 .or. liw/=20) stop 2
  case default
    stop 98
  end select
  if(len_trim(mode)>0)then
    call zheevr_workspace_lengths(2, cmplx(query,0d0,real64), 48d0, 20, lw, lrw, liw)
    stop 99
  end if

  ! Hermitian matrix [[2,i],[-i,2]], eigenvalues 1 and 3.
  a(1,:) = [(2d0,0d0),(0d0,1d0)]
  a(2,:) = [(0d0,-1d0),(2d0,0d0)]
  call zheevr('V','A','U',2,a,2,0d0,0d0,1,2,0d0,m,w,z,2,support,wq,-1,rq,-1,iq,-1,info)
  call check_lapack_info('test ZHEEVR query',info)
  call zheevr_workspace_lengths(2,wq(1),rq(1),iq(1),lw,lrw,liw)
  allocate(work(lw),rwork(lrw),iwork(liw))
  call zheevr('V','A','U',2,a,2,0d0,0d0,1,2,0d0,m,w,z,2,support,work,lw,rwork,lrw,iwork,liw,info)
  call check_lapack_info('test ZHEEVR solve',info)
  if(m/=2 .or. maxval(abs(w-[1d0,3d0]))>1d-12) stop 3
  a(1,:) = [(2d0,0d0),(0d0,1d0)]
  a(2,:) = [(0d0,-1d0),(2d0,0d0)]
  if(maxval(abs(matmul(a,z)-z*spread(w,1,2)))>1d-12) stop 4
  ! JOBZ=N and RANGE=V exercise correctly typed VL/VU and one-row Z.
  call zheevr('N','V','U',2,a,2,0d0,2d0,1,2,0d0,m,w,z,1,support,work,lw,rwork,lrw,iwork,liw,info)
  call check_lapack_info('test ZHEEVR values only',info)
  if(m/=1 .or. abs(w(1)-1d0)>1d-12) stop 5
  ! Tridiagonal matrix [[2,1],[1,2]], both full and index-selected spectra.
  d=2d0; e=1d0
  call dstevx('V','A',2,d,e,0d0,0d0,1,2,2d0*dlamch('S'),m,w,tz,2,tw,tiw,fail,info)
  call check_lapack_info('test DSTEVX full',info)
  if(m/=2 .or. maxval(abs(w-[1d0,3d0]))>1d-12) stop 6
  d=2d0; e=1d0
  call dstevx('V','I',2,d,e,0d0,0d0,2,2,2d0*dlamch('S'),m,w,tz,2,tw,tiw,fail,info)
  call check_lapack_info('test DSTEVX index',info)
  if(m/=1 .or. abs(w(1)-3d0)>1d-12) stop 7
  print *, 'PASS: ZHEEVR/DSTEVX spectra, residuals, query ranks and workspace guards'
end program lapack_interface_check
