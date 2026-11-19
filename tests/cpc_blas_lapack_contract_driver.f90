! Native regression for explicit BLAS/LAPACK LP64 interfaces in ORIGINAL ED2.
! Uses the existing production dense diagonalization, not a mock solver.
program cpc_blas_lapack_contract_driver
  use iso_fortran_env, only: real64
  use input_param, only: OUTDIR
  use eigen_solver, only: fulldia_routines_with_momentum
  use blas_lapack_interfaces, only: zdotc, dznrm2, zcopy, dlamch, &
    zheevr, zheevr_workspace_lengths, check_lapack_info
  implicit none
  complex(real64) :: h1(1,1), h2(2,2), h3(3,3), old3(3,3), x(3), y(3), cv(3)
  real(real64) :: ev1(1), ev2(2), ev3(3), err, dnorm, tol
  integer :: i, j
  OUTDIR='./'
  tol=2.0e-12_real64
  if(.not.(dlamch('E')>0.0_real64))error stop 'LAPACK machine epsilon unavailable'
  call check_small_workspace(0)
  call check_small_workspace(1)
  h1(1,1)=(-0.385_real64,0.0_real64)
  call fulldia_routines_with_momentum(1,h1,ev1,'V')
  if(abs(ev1(1)+0.385_real64)>tol)error stop 'ZHEEVR scalar eigenvalue mismatch'
  if(abs(abs(h1(1,1))-1.0_real64)>tol)error stop 'ZHEEVR scalar eigenvector not normalized'
  print '(a)', 'LAPACK_TYPED_ANALYTIC_1X1_PASS'
  h2=(0.0_real64,0.0_real64)
  h2(1,1)=(2.0_real64,0.0_real64)
  h2(2,2)=(4.0_real64,0.0_real64)
  h2(1,2)=(1.0_real64,1.0_real64)
  h2(2,1)=conjg(h2(1,2))
  call fulldia_routines_with_momentum(2,h2,ev2,'V')
  if(abs(ev2(1)-(3.0_real64-sqrt(3.0_real64)))>tol) &
    error stop 'Complex ZHEEVR analytic low eigenvalue mismatch'
  if(abs(ev2(2)-(3.0_real64+sqrt(3.0_real64)))>tol) &
    error stop 'Complex ZHEEVR analytic high eigenvalue mismatch'
  print '(a)', 'LAPACK_TYPED_ANALYTIC_2X2_PASS'

  h3=(0.0_real64,0.0_real64)
  h3(1,1)=(1.0_real64,0.0_real64)
  h3(2,2)=(3.0_real64,0.0_real64)
  h3(3,3)=(5.0_real64,0.0_real64)
  h3(1,2)=(0.25_real64,0.30_real64)
  h3(2,1)=conjg(h3(1,2))
  h3(2,3)=(-0.41_real64,0.18_real64)
  h3(3,2)=conjg(h3(2,3))
  h3(1,3)=(0.03_real64,-0.11_real64)
  h3(3,1)=conjg(h3(1,3))
  old3=h3
  call fulldia_routines_with_momentum(3,h3,ev3,'V')
  if(any(ev3(2:3)<ev3(1:2)))error stop 'ZHEEVR incorrect eigen ordering'
  do j=1,3
    cv=matmul(old3,h3(:,j))-ev3(j)*h3(:,j)
    err=maxval(abs(cv))
    if(.not.(err<tol))error stop 'ZHEEVR complex Hermitian residual failed'
    dnorm=sqrt(sum(abs(h3(:,j))**2))
    if(abs(dnorm-1.0_real64)>tol)error stop 'ZHEEVR eigenvector not normalized'
  end do
  do i=1,3
    do j=1,3
      err=abs(sum(conjg(h3(:,i))*h3(:,j)))
      if(i==j)err=abs(err-1.0_real64)
      if(.not.(err<tol))error stop 'ZHEEVR eigenvectors not orthogonal'
    end do
  end do
  print '(a)', 'LAPACK_TYPED_HERMITIAN_3X3_PASS'

  x=[(1.0_real64,2.0_real64),(3.0_real64,-1.0_real64),(-2.0_real64,1.0_real64)]
  y=[(2.0_real64,1.0_real64),(0.5_real64,2.0_real64),(1.0_real64,-3.0_real64)]
  err=abs(zdotc(3,x,1,y,1)-sum(conjg(x)*y))
  if(.not.(err<tol))error stop 'ZDOTC Fortran complex return ABI mismatch'
  err=abs(dznrm2(3,x,1)-sqrt(sum(abs(x)**2)))
  if(.not.(err<tol))error stop 'DZNRM2 argument or result kind mismatch'
  call zcopy(3,x,1,cv,1)
  if(any(cv/=x))error stop 'ZCOPY complex ABI mismatch'
  print '(a)', 'BLAS_TYPED_ZDOTC_ZCOPY_DZNRM2_PASS'
  print '(a)', 'BLAS_LAPACK_TYPED_INTERFACE_REGRESSION_PASS'
contains
  subroutine check_small_workspace(n)
    integer, intent(in) :: n
    integer :: m, info, lw, lrw, liw, iwq(1), isuppz(2)
    integer, allocatable :: iw(:)
    complex(real64) :: a(1,1), z(1,1), wq(1)
    complex(real64), allocatable :: work(:)
    real(real64) :: w(1), rwq(1)
    real(real64), allocatable :: rw(:)
    a(1,1)=(-0.385_real64,0.0_real64)
    ! Check the valid minimum query independently of backend block choices.
    call zheevr_workspace_lengths(n,(1.0_real64,0.0_real64),1.0_real64,1,lw,lrw,liw)
    if(lw/=1.or.lrw/=1.or.liw/=1)error stop 'Small ZHEEVR workspace minimum rejected'
    call zheevr('V','A','U',n,a,1,0.0_real64,0.0_real64,1,n,0.0_real64, &
      m,w,z,1,isuppz,wq,-1,rwq,-1,iwq,-1,info)
    call check_lapack_info('ZHEEVR small workspace query',info)
    call zheevr_workspace_lengths(n,wq(1),rwq(1),iwq(1),lw,lrw,liw)
    allocate(work(lw),rw(lrw),iw(liw))
    ! Execute with exactly the returned lengths, including one-element
    ! real/integer buffers when the backend reports the scalar quick return.
    call zheevr('V','A','U',n,a,1,0.0_real64,0.0_real64,1,n,0.0_real64, &
      m,w,z,1,isuppz,work,lw,rw,lrw,iw,liw,info)
    call check_lapack_info('ZHEEVR small solve',info)
    if(m/=n)error stop 'Small ZHEEVR returned wrong eigenvalue count'
    if(n==1)then
      if(abs(w(1)+0.385_real64)>tol)error stop 'Small ZHEEVR eigenvalue mismatch'
      if(abs(abs(z(1,1))-1.0_real64)>tol)error stop 'Small ZHEEVR vector not normalized'
    end if
    print '(a,i0,a,3(i0,1x))', 'LAPACK_TYPED_SMALL_WORKSPACE_PASS N=',n,' lengths=',lw,lrw,liw
  end subroutine check_small_workspace
end program cpc_blas_lapack_contract_driver
