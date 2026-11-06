! Test-only instrumentation for external benchmark copies. No solver arithmetic.
module phase1d_probe
  use omp_lib
  implicit none
  integer :: probe_h_calls=0
  real(8) :: started
contains
  subroutine probe_begin()
    integer :: n,i,tid,seed,rc
    integer, allocatable :: state(:),records(:,:)
    character(32) :: text
    call get_environment_variable('ED2_BENCH_SEED',text,status=rc)
    if(rc/=0) error stop 'benchmark seed missing'
    read(text,*) seed
    call random_seed(size=n)
    allocate(state(n),records(n,omp_get_max_threads()))
    !$omp parallel private(tid,state,i)
    tid=omp_get_thread_num()
    do i=1,n
      state(i)=modulo(seed+104729*tid+37*i,huge(1)-1)
    end do
    records(:,tid+1)=state
    call random_seed(put=state)
    !$omp end parallel
    open(91,file='benchmark-rng.dat',status='new')
    do i=1,size(records,2)
      write(91,*) i-1,records(:,i)
    end do
    close(91)
    probe_h_calls=0
    started=omp_get_wtime()
  end subroutine
  subroutine probe_end()
    integer :: ios,kib
    character(256) :: line
    kib=-1
    open(94,file='/proc/self/status',status='old',iostat=ios)
    if(ios==0)then
      do
        read(94,'(a)',iostat=ios) line
        if(ios/=0)exit
        if(index(line,'VmHWM:')==1)read(line(7:),*) kib
      end do
      close(94)
    end if
    open(92,file='benchmark-metrics.dat',status='new')
    write(92,*) probe_h_calls,omp_get_wtime()-started
    write(92,*) kib
    close(92)
  end subroutine
end module phase1d_probe
