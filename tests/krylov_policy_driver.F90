module policy_operator
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
end module policy_operator

! Test seam to the exact numerical kernels used by ALG=1/2. H is supplied by
! Python independently; no production Hamiltonian/symmetry code is replaced.
program krylov_policy_driver
  use iso_fortran_env, only: int64
  use policy_operator, only: h,apply,h_calls
  use omp_lib, only: omp_get_thread_num,omp_get_max_threads,omp_get_wtime
  use krylov_solvers, only: lanczos_ground, restarted_hermitian
#ifdef PHASE1F
  use krylov_observations, only: krylov_diagnostics,write_krylov_diagnostics
#endif
#ifdef PHASE1G
  use krylov_residual_policy, only: residual_policy
#endif
  implicit none
#ifdef PHASE1G
  type(residual_policy), allocatable :: policy
#endif
  logical :: policy_enabled=.false.,scale_provided=.false.
  real(8) :: residual_atol=0d0,residual_rtol=0d0,residual_scale=0d0,energy_origin=0d0
#ifdef PHASE1F
  type(krylov_diagnostics), allocatable :: observation
  integer, allocatable :: seed_option
#endif
  logical :: enable_diagnostics=.true.,use_seed_option=.false.
  integer :: seed,seed_size,tid,k,ios,peak_rss
  character(256) :: line
  integer, allocatable :: seed_state(:),seed_records(:,:)
  real(8) :: started,elapsed
  integer :: n,alg,wanted,keep,capacity,maxiter,miniter,interval,i,j,status,iterations
  real(8) :: tol,re,im
  logical :: use_start
  namelist /problem/ n,alg,wanted,keep,capacity,maxiter,miniter,interval,tol,use_start,seed,enable_diagnostics,use_seed_option, &
                       policy_enabled,scale_provided,residual_atol,residual_rtol,residual_scale,energy_origin
  complex(8), allocatable :: start(:),vectors(:,:)
  real(8), allocatable :: values(:),residuals(:)
  character(120) :: reason
  read(*,problem)
  call random_seed(size=seed_size)
  allocate(seed_state(seed_size),seed_records(seed_size,omp_get_max_threads()))
  !$omp parallel private(tid,seed_state,k)
  tid=omp_get_thread_num()
  do k=1,seed_size
    seed_state(k)=int(modulo(int(seed,int64)+104729_int64*tid+37_int64*k,int(huge(1)-1,int64)))
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
#ifdef PHASE1F
  if(enable_diagnostics)allocate(observation)
  if(use_seed_option)then
    allocate(seed_option);seed_option=seed
  end if
#endif
#ifdef PHASE1G
  if(policy_enabled)then
    allocate(policy)
    policy%atol=residual_atol;policy%rtol=residual_rtol
    policy%operator_scale=residual_scale;policy%origin=energy_origin;policy%scale_provided=scale_provided
  end if
#endif
  started=omp_get_wtime()
  if(alg==1)then
    if(use_start)then
      call lanczos_ground(apply,n,tol,maxiter,miniter,interval,values(1),vectors(:,1), &
                          iterations,status,reason,residuals(1),start &
#ifdef PHASE1F
                          ,diagnostics=observation,rng_seed=seed_option &
#endif
#ifdef PHASE1G
                          ,policy=policy &
#endif
                          )
    else
      call lanczos_ground(apply,n,tol,maxiter,miniter,interval,values(1),vectors(:,1), &
                          iterations,status,reason,residuals(1) &
#ifdef PHASE1F
                          ,diagnostics=observation,rng_seed=seed_option &
#endif
#ifdef PHASE1G
                          ,policy=policy &
#endif
                          )
    end if
  else
    if(use_start)then
      call restarted_hermitian(apply,n,wanted,keep,capacity,tol,maxiter,values,vectors, &
                               iterations,status,reason,residuals,start &
#ifdef PHASE1F
                               ,diagnostics=observation,rng_seed=seed_option &
#endif
#ifdef PHASE1G
                          ,policy=policy &
#endif
                          )
    else
      call restarted_hermitian(apply,n,wanted,keep,capacity,tol,maxiter,values,vectors, &
                               iterations,status,reason,residuals &
#ifdef PHASE1F
                               ,diagnostics=observation,rng_seed=seed_option &
#endif
#ifdef PHASE1G
                          ,policy=policy &
#endif
                          )
    end if
  end if
  elapsed=omp_get_wtime()-started
#ifdef PHASE1F
  if(allocated(observation))call write_krylov_diagnostics('diagnostics.dat',observation)
#endif
  !$omp parallel private(tid,seed_state)
  tid=omp_get_thread_num()
  call random_seed(get=seed_state)
  seed_records(:,tid+1)=seed_state
  !$omp end parallel
  open(70,file='rng-after.dat',status='new')
  do k=1,size(seed_records,2)
    write(70,*) k-1,seed_records(:,k)
  end do
  close(70)
  open(74,file='metrics.dat',status='new')
  peak_rss=-1
  open(77,file='/proc/self/status',status='old',iostat=ios)
  if(ios==0)then
    do
      read(77,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(index(line,'VmHWM:')==1)read(line(7:),*)peak_rss
    end do
    close(77)
  end if
  write(74,*) h_calls,elapsed,peak_rss
  close(74)
  open(71,file='status.dat',status='new')
  write(71,'(i0)') status
  write(71,'(a)') trim(reason)
  write(71,'(i0)') iterations
  close(71)
  open(75,file='last-estimates.dat',status='new')
  do j=1,wanted
    write(75,'(2es26.16e3)') values(j),residuals(j)
  end do
  close(75)
  if(status==1 .and. alg==2)then
    open(76,file='last-vectors.dat',status='new')
    do j=1,wanted
      do i=1,n
        write(76,'(2es26.16e3)') real(vectors(i,j),8),aimag(vectors(i,j))
      end do
    end do
    close(76)
  end if
  if(status/=0)stop 1
  open(72,file='eigenpairs.dat',status='new')
  do j=1,wanted
    write(72,'(i8,2es26.16e3)') j,values(j),residuals(j)
  end do
  close(72)
  open(73,file='vectors.dat',status='new')
  do j=1,wanted
    do i=1,n
      write(73,'(2i8,2es26.16e3)') j,i,real(vectors(i,j),8),aimag(vectors(i,j))
    end do
  end do
  close(73)
end program krylov_policy_driver
