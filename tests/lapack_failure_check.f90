! This executable overrides ZHEEVR only in the test link, not in QS3ED2.
module lapack_failure_control
  implicit none
  character(20) :: mode
end module lapack_failure_control

program lapack_failure_check
  use lapack_failure_control
  use eigen_solver, only: fulldia_routines_with_momentum
  implicit none
  complex(8) :: a(2,2)
  real(8) :: w(2)
  call get_command_argument(1,mode)
  a=(0d0,0d0)
  call fulldia_routines_with_momentum(2,a,w,'V')
  ! A successful exit here means the caller incorrectly ignored injected INFO.
end program lapack_failure_check

subroutine zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
                  m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
  use iso_fortran_env, only: real64, int32
  use lapack_failure_control
  implicit none
  character(1), intent(in) :: jobz, range, uplo
  integer(int32), intent(in) :: n, lda, il, iu, ldz, lwork, lrwork, liwork
  real(real64), intent(in) :: vl, vu, abstol
  complex(real64), intent(inout) :: a(lda,*)
  integer(int32), intent(out) :: m, isuppz(*), iwork(*), info
  real(real64), intent(out) :: w(*), rwork(*)
  complex(real64), intent(out) :: z(ldz,*), work(*)
  if(mode=='query' .and. lwork==-1)then
    ! Do not define query sizes: a failing query must not be consumed.
    info=-18
    return
  end if
  if(lwork==-1 .and. lrwork==-1 .and. liwork==-1)then
    work(1)=cmplx(max(1,2*n),0,real64)
    rwork(1)=real(max(1,24*n),real64)
    iwork(1)=max(1,10*n)
    info=0
  else
    info=2
  end if
end subroutine zheevr
