module stability_operator
  implicit none
  complex(8), allocatable :: h(:,:)
contains
  subroutine apply(dim,x,y)
    integer, intent(in) :: dim
    complex(8), intent(in) :: x(dim)
    complex(8), intent(out) :: y(dim)
    y=matmul(h,x)
  end subroutine apply
end module stability_operator

! Test seam to the exact numerical kernels used by ALG=1/2. H is supplied by
! Python independently; no production Hamiltonian/symmetry code is replaced.
program lanczos_stability_driver
  use stability_operator, only: h,apply
  use krylov_solvers, only: lanczos_ground, restarted_hermitian
  implicit none
  integer :: n,alg,wanted,keep,capacity,maxiter,miniter,interval,i,j,status,iterations
  real(8) :: tol,re,im
  logical :: use_start
  namelist /problem/ n,alg,wanted,keep,capacity,maxiter,miniter,interval,tol,use_start
  complex(8), allocatable :: start(:),vectors(:,:)
  real(8), allocatable :: values(:),residuals(:)
  character(120) :: reason
  read(*,problem)
  allocate(h(n,n),start(n),vectors(n,wanted),values(max(n,capacity)),residuals(wanted))
  do i=1,n
    do j=1,n
      read(*,*) re,im
      h(i,j)=cmplx(re,im,8)
    end do
  end do
  do i=1,n
    read(*,*) re,im
    start(i)=cmplx(re,im,8)
  end do
  if(alg==1)then
    if(use_start)then
      call lanczos_ground(apply,n,tol,maxiter,miniter,interval,values(1),vectors(:,1), &
                          iterations,status,reason,residuals(1),start)
    else
      call lanczos_ground(apply,n,tol,maxiter,miniter,interval,values(1),vectors(:,1), &
                          iterations,status,reason,residuals(1))
    end if
  else
    if(use_start)then
      call restarted_hermitian(apply,n,wanted,keep,capacity,tol,maxiter,values,vectors, &
                               iterations,status,reason,residuals,start)
    else
      call restarted_hermitian(apply,n,wanted,keep,capacity,tol,maxiter,values,vectors, &
                               iterations,status,reason,residuals)
    end if
  end if
  open(71,file='status.dat',status='new')
  write(71,'(i0)') status
  write(71,'(a)') trim(reason)
  write(71,'(i0)') iterations
  close(71)
  if(status/=0)stop 1
  open(72,file='eigenpairs.dat',status='new')
  do j=1,wanted
    write(72,'(i8,2es25.15)') j,values(j),residuals(j)
  end do
  close(72)
  open(73,file='vectors.dat',status='new')
  do j=1,wanted
    do i=1,n
      write(73,'(2i8,2es25.15)') j,i,real(vectors(i,j),8),aimag(vectors(i,j))
    end do
  end do
  close(73)
end program lanczos_stability_driver
