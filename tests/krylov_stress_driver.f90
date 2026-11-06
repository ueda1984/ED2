module stress_operator
  implicit none
  complex(8), allocatable :: h(:,:)
  integer :: h_calls=0
contains
  subroutine apply(dim,x,y)
    integer, intent(in) :: dim
    complex(8), intent(in) :: x(dim)
    complex(8), intent(out) :: y(dim)
    h_calls=h_calls+1
    y=matmul(h,x)
  end subroutine apply
end module stress_operator

! Test seam to the exact numerical kernels used by ALG=1/2. H is supplied by
! Python independently; no production Hamiltonian/symmetry code is replaced.
program krylov_stress_driver
  use stress_operator, only: h,apply,h_calls
  use omp_lib, only: omp_get_thread_num,omp_get_max_threads,omp_get_wtime
  use krylov_solvers, only: lanczos_ground, restarted_hermitian
  implicit none
  integer :: seed,seed_size,tid,k
  integer, allocatable :: seed_state(:),seed_records(:,:)
  real(8) :: started,elapsed
  integer :: n,alg,wanted,keep,capacity,maxiter,miniter,interval,i,j,status,iterations
  real(8) :: tol,re,im
  logical :: use_start
  namelist /problem/ n,alg,wanted,keep,capacity,maxiter,miniter,interval,tol,use_start,seed
  complex(8), allocatable :: start(:),vectors(:,:)
  real(8), allocatable :: values(:),residuals(:)
  character(120) :: reason
  read(*,problem)
  call random_seed(size=seed_size)
  allocate(seed_state(seed_size),seed_records(seed_size,omp_get_max_threads()))
  !$omp parallel private(tid,seed_state,k)
  tid=omp_get_thread_num()
  do k=1,seed_size
    seed_state(k)=modulo(seed+104729*tid+37*k,huge(1)-1)
  end do
  seed_records(:,tid+1)=seed_state
  call random_seed(put=seed_state)
  !$omp end parallel
  open(70,file="rng.dat",status="new")
  do k=1,size(seed_records,2)
    write(70,*) k-1,seed_records(:,k)
  end do
  close(70)
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
  started=omp_get_wtime()
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
  elapsed=omp_get_wtime()-started
  open(74,file='metrics.dat',status='new')
  write(74,*) h_calls,elapsed
  close(74)
  open(71,file='status.dat',status='new')
  write(71,'(i0)') status
  write(71,'(a)') trim(reason)
  write(71,'(i0)') iterations
  close(71)
  open(75,file='last-estimates.dat',status='new')
  do j=1,wanted
    write(75,'(2es25.15)') values(j),residuals(j)
  end do
  close(75)
  if(status==1 .and. alg==2)then
    open(76,file='last-vectors.dat',status='new')
    do j=1,wanted
      do i=1,n
        write(76,'(2es25.15)') real(vectors(i,j),8),aimag(vectors(i,j))
      end do
    end do
    close(76)
  end if
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
end program krylov_stress_driver
