! Native Fortran LP64 interfaces for the BLAS/LAPACK subset used by ED2.
! No bind(C): in particular, ZDOTC uses the linked Fortran complex-result ABI.
module blas_lapack_interfaces
  use iso_fortran_env, only: real64, int32, int64, error_unit
  implicit none
  private
  public :: zheevr, dstevx, dlamch, dznrm2, zdotc
  public :: zgemm, zcopy, zaxpy, zscal, zdscal
  public :: check_lapack_info, zheevr_workspace_lengths

  interface
    subroutine zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
                      m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
      import :: real64, int32
      character(1), intent(in) :: jobz, range, uplo
      integer(int32), intent(in) :: n, lda, il, iu, ldz, lwork, lrwork, liwork
      real(real64), intent(in) :: vl, vu, abstol
      complex(real64), intent(inout) :: a(lda,*)
      integer(int32), intent(out) :: m, isuppz(*), iwork(*), info
      real(real64), intent(out) :: w(*), rwork(*)
      complex(real64), intent(out) :: z(ldz,*), work(*)
    end subroutine zheevr

    subroutine dstevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info)
      import :: real64, int32
      character(1), intent(in) :: jobz, range
      integer(int32), intent(in) :: n, il, iu, ldz
      real(real64), intent(in) :: vl, vu, abstol
      real(real64), intent(inout) :: d(*), e(*)
      integer(int32), intent(out) :: m, iwork(*), ifail(*), info
      real(real64), intent(out) :: w(*), z(ldz,*), work(*)
    end subroutine dstevx

    function dlamch(cmach) result(value)
      import :: real64
      character(1), intent(in) :: cmach
      real(real64) :: value
    end function dlamch

    function dznrm2(n, x, incx) result(value)
      import :: real64, int32
      integer(int32), intent(in) :: n, incx
      complex(real64), intent(in) :: x(*)
      real(real64) :: value
    end function dznrm2

    function zdotc(n, x, incx, y, incy) result(value)
      import :: real64, int32
      integer(int32), intent(in) :: n, incx, incy
      complex(real64), intent(in) :: x(*), y(*)
      complex(real64) :: value
    end function zdotc

    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: real64, int32
      character(1), intent(in) :: transa, transb
      integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
      complex(real64), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
      complex(real64), intent(inout) :: c(ldc,*)
    end subroutine zgemm

    subroutine zcopy(n, x, incx, y, incy)
      import :: real64, int32
      integer(int32), intent(in) :: n, incx, incy
      complex(real64), intent(in) :: x(*)
      complex(real64), intent(out) :: y(*)
    end subroutine zcopy

    subroutine zaxpy(n, alpha, x, incx, y, incy)
      import :: real64, int32
      integer(int32), intent(in) :: n, incx, incy
      complex(real64), intent(in) :: alpha, x(*)
      complex(real64), intent(inout) :: y(*)
    end subroutine zaxpy

    subroutine zscal(n, alpha, x, incx)
      import :: real64, int32
      integer(int32), intent(in) :: n, incx
      complex(real64), intent(in) :: alpha
      complex(real64), intent(inout) :: x(*)
    end subroutine zscal

    subroutine zdscal(n, alpha, x, incx)
      import :: real64, int32
      integer(int32), intent(in) :: n, incx
      real(real64), intent(in) :: alpha
      complex(real64), intent(inout) :: x(*)
    end subroutine zdscal
  end interface
contains
  subroutine check_lapack_info(routine, info)
    character(*), intent(in) :: routine
    integer(int32), intent(in) :: info
    if (info /= 0) then
      write(error_unit,'(a,a,i0)') trim(routine), ': LAPACK INFO = ', info
      error stop 1
    end if
  end subroutine check_lapack_info

  subroutine zheevr_workspace_lengths(n, work_query, rwork_query, iwork_query, lwork, lrwork, liwork)
    integer(int32), intent(in) :: n, iwork_query
    complex(real64), intent(in) :: work_query
    real(real64), intent(in) :: rwork_query
    integer(int32), intent(out) :: lwork, lrwork, liwork
    ! ZHEEVR requires max(1,2*N), max(1,24*N), max(1,10*N).
    ! Check before integer multiplication/conversion, including LP64 overflow.
    if (n < 0 .or. n > huge(n)/24) then
      write(error_unit,'(a)') 'ZHEEVR: dimension exceeds LP64 workspace limits'
      error stop 1
    end if
    lwork = checked_length(real(work_query,real64), max(1,2*n))
    lrwork = checked_length(rwork_query, max(1,24*n))
    if (iwork_query < max(1,10*n)) then
      write(error_unit,'(a)') 'ZHEEVR: invalid integer workspace query'
      error stop 1
    end if
    liwork = iwork_query
  end subroutine zheevr_workspace_lengths

  function checked_length(value, minimum) result(length)
    real(real64), intent(in) :: value
    integer(int32), intent(in) :: minimum
    integer(int32) :: length
    ! LAPACK double precision is IEEE binary64 on supported ED2 platforms.
    ! Use the exponent bits: -Ofast/-ffinite-math-only may otherwise fold an
    ! IEEE_IS_FINITE check to true, defeating validation of backend output.
    if (iand(transfer(value,0_int64),int(z'7ff0000000000000',int64)) == int(z'7ff0000000000000',int64)) then
      write(error_unit,'(a)') 'ZHEEVR: nonfinite workspace query'
      error stop 1
    end if
    if (value < real(minimum,real64) .or. value > real(huge(length),real64)) then
      write(error_unit,'(a)') 'ZHEEVR: invalid or overflowing workspace query'
      error stop 1
    end if
    ! Round up rather than underallocate if a backend returns a fractional value.
    length = ceiling(value,kind=int32)
  end function checked_length
end module blas_lapack_interfaces
