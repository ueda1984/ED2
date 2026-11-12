! Internal telemetry/CLI controls; not a stable public library API.
module krylov_observations
  use iso_fortran_env, only: int64,compiler_version
  use blas_lapack_interfaces, only: dznrm2,zdotc
  use krylov_residual_policy, only: residual_policy
  !$ use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: krylov_diagnostics, begin_observation, finish_observation, observe_pair, observe_gram
  public :: read_krylov_controls, write_krylov_diagnostics,read_residual_controls
  public :: record_residual_policy,observe_requested_residual
  type :: krylov_diagnostics
    integer :: algorithm=0,dimension=0,wanted=0,status=-1,iterations=0,cycles=0
    integer :: requested_budget=0,effective_budget=0,miniter=0,interval=0,keep=0,capacity=0
    integer :: energy_consecutive_checks=0
    logical :: energy_gate_satisfied=.false.
    character(120) :: candidate_role='none observed'
    integer :: exact_closures=0,near_closures=0,random_vectors=0,max_threads=1,rng_active_threads=0
    integer(int64) :: h_applications=0_int64
    real(8) :: energy_tolerance=0d0,operator_scale=0d0,breakdown_threshold=0d0
    real(8) :: energy_change=0d0,last_beta=0d0,closure=0d0,orthogonality_error=0d0
    logical :: scale_available=.false.
    logical :: energy_change_available=.false.,orthogonality_available=.false.
    logical :: explicit_initial=.false.,seed_specified=.false.,roundoff_risk=.false.
    logical :: roundoff_limit_proven=.false.
    logical :: residual_policy_enabled=.false.,requested_scale_provided=.false.
    logical :: all_requested_residuals_met=.false.
    real(8) :: requested_atol=0d0,requested_rtol=0d0,requested_operator_scale=0d0,requested_origin=0d0
    real(8), allocatable :: requested_scale(:),requested_threshold(:)
    logical, allocatable :: requested_met(:),combined_certified(:)
    integer :: seed=0
    integer, allocatable :: rng_initial_state(:,:)
    character(120) :: reason='',termination='not_started'
    character(120) :: scale_provenance='maximum observed H action norm; estimate, not a spectral bound'
    character(120) :: inclusion='not_assessed; no independent lowest-state evidence'
    character(120) :: roundoff_assessment='not_assessed; no arithmetic-error enclosure'
    character(256) :: rng_implementation=''
    real(8), allocatable :: ritz(:),rayleigh(:),residual(:),norm(:),residual_threshold(:)
    logical, allocatable :: pair_available(:),residual_certified(:)
  end type
contains
  subroutine begin_observation(d,alg,n,wanted,budget,tol,explicit_initial,seed)
    type(krylov_diagnostics), intent(out) :: d
    integer, intent(in) :: alg,n,wanted,budget
    real(8), intent(in) :: tol
    logical, intent(in) :: explicit_initial
    integer, intent(in), optional :: seed
    d%algorithm=alg;d%dimension=n;d%wanted=wanted;d%requested_budget=budget
    d%energy_tolerance=tol;d%explicit_initial=explicit_initial
    d%rng_implementation=compiler_version()
    !$ d%max_threads=omp_get_max_threads()
    if(present(seed))then
      d%seed_specified=.true.;d%seed=seed
    end if
    ! Advisory only; does not change any convergence or failure decision.
    if(iand(transfer(tol,0_int64),int(z'7ff0000000000000',int64))/=int(z'7ff0000000000000',int64)) &
      d%roundoff_risk=tol>0d0 .and. tol<epsilon(1d0)
    if(d%roundoff_risk) d%roundoff_assessment='requested energy change below binary64 epsilon; heuristic risk, not proven cause'
  end subroutine

  subroutine finish_observation(d,status,reason,iterations)
    type(krylov_diagnostics), intent(inout) :: d
    integer, intent(in) :: status,iterations
    character(*), intent(in) :: reason
    d%status=status;d%reason=reason;d%iterations=iterations
    select case(status)
    case(-1)
      d%termination='invalid_input'
    case(1)
      d%termination='budget_exhausted'
    case(2)
      d%termination='numerical_failure'
    case(0)
      d%termination='legacy_policy_satisfied'
      if(d%residual_policy_enabled)d%termination='legacy_and_requested_residual_satisfied'
    end select
    if(index(reason,'happy breakdown')>0) d%termination='exact_subspace_closed'
    if(index(reason,'numerical near-breakdown')>0) d%termination='numerical_subspace_closed'
    if(status/=0 .and. index(reason,'breakdown')>0) d%termination='uncertified_numerical_breakdown'
    if(allocated(d%pair_available))then
      d%candidate_role='last observed diagnostic candidate; not a successful output'
      if(status==0)d%candidate_role='returned successful candidate under the legacy policy'
    end if
    if(d%residual_policy_enabled .and. allocated(d%requested_met))then
      d%all_requested_residuals_met=all(d%requested_met .and. d%pair_available)
      if(status==0)d%candidate_role='returned candidate under legacy plus requested residual policy'
    end if
    ! A residual certificate never supplies independent lowest-state evidence.
  end subroutine

  subroutine observe_pair(d,j,count,value,vector,action,residual,threshold)
    type(krylov_diagnostics), intent(inout) :: d
    integer, intent(in) :: j,count
    real(8), intent(in) :: value,residual,threshold
    complex(8), intent(in) :: vector(:),action(:)
    if(.not.allocated(d%ritz))then
      allocate(d%ritz(count),d%rayleigh(count),d%residual(count),d%norm(count), &
               d%residual_threshold(count),d%pair_available(count),d%residual_certified(count))
      d%ritz=0d0;d%rayleigh=0d0;d%residual=0d0;d%norm=0d0;d%residual_threshold=0d0
      d%pair_available=.false.;d%residual_certified=.false.
    end if
    d%ritz(j)=value;d%residual(j)=residual;d%residual_threshold(j)=threshold
    d%norm(j)=dznrm2(size(vector),vector,1)
    d%rayleigh(j)=real(zdotc(size(vector),vector,1,action,1),8)/(d%norm(j)**2)
    d%pair_available(j)=.true.;d%residual_certified(j)=residual<=threshold
    if(count==1)then
      d%orthogonality_error=abs(d%norm(j)**2-1d0);d%orthogonality_available=.true.
    end if
  end subroutine

  subroutine observe_gram(d,vectors)
    type(krylov_diagnostics), intent(inout) :: d
    complex(8), intent(in) :: vectors(:,:)
    integer :: i,j
    complex(8) :: g
    if(.not.allocated(d%pair_available))return
    if(.not.all(d%pair_available))return
    d%orthogonality_error=0d0
    do j=1,size(vectors,2)
      do i=1,size(vectors,2)
        g=zdotc(size(vectors,1),vectors(:,i),1,vectors(:,j),1)
        if(i==j)g=g-1d0
        d%orthogonality_error=max(d%orthogonality_error,abs(g))
      end do
    end do
    d%orthogonality_available=.true.
  end subroutine

  subroutine record_residual_policy(d,p)
    type(krylov_diagnostics), intent(inout) :: d
    type(residual_policy), intent(in) :: p
    d%residual_policy_enabled=.true.;d%requested_scale_provided=p%scale_provided
    d%requested_atol=p%atol;d%requested_rtol=p%rtol
    d%requested_operator_scale=p%operator_scale;d%requested_origin=p%origin
  end subroutine

  subroutine observe_requested_residual(d,j,count,scale,threshold,met)
    type(krylov_diagnostics), intent(inout) :: d
    integer, intent(in) :: j,count
    real(8), intent(in) :: scale,threshold
    logical, intent(in) :: met
    if(.not.allocated(d%requested_met))then
      allocate(d%requested_scale(count),d%requested_threshold(count),d%requested_met(count),d%combined_certified(count))
      d%requested_scale=0d0;d%requested_threshold=0d0;d%requested_met=.false.;d%combined_certified=.false.
    end if
    d%requested_scale(j)=scale;d%requested_threshold(j)=threshold;d%requested_met(j)=met
    d%combined_certified(j)=d%residual_certified(j) .and. met
  end subroutine

  subroutine read_control_real(text,value)
    character(*), intent(in) :: text
    real(8), intent(out) :: value
    character(1) :: extra
    integer :: rc
    real(8) :: parsed
    read(text,*,iostat=rc)value
    if(rc/=0)error stop 'invalid residual policy real value'
    read(text,*,iostat=rc)parsed,extra
    if(rc>=0)error stop 'extra data in residual policy real value'
  end subroutine

  subroutine read_residual_controls(p)
    type(residual_policy), allocatable, intent(out) :: p
    character(:), allocatable :: atol,rtol,scale,origin
    call environment_text('ED2_RESIDUAL_ATOL',atol)
    call environment_text('ED2_RESIDUAL_RTOL',rtol)
    call environment_text('ED2_RESIDUAL_SCALE',scale)
    call environment_text('ED2_ENERGY_ORIGIN',origin)
    if(len(atol)==0 .and. len(rtol)==0)then
      if(len(scale)>0 .or. len(origin)>0)error stop 'residual scale/origin requires an explicit residual tolerance'
      return
    end if
    allocate(p)
    if(len(atol)>0)call read_control_real(atol,p%atol)
    if(len(rtol)>0)call read_control_real(rtol,p%rtol)
    if(len(scale)>0)then
      call read_control_real(scale,p%operator_scale);p%scale_provided=.true.
    end if
    if(len(origin)>0)call read_control_real(origin,p%origin)
  end subroutine

  subroutine environment_text(name,text)
    character(*), intent(in) :: name
    character(:), allocatable, intent(out) :: text
    integer :: n,rc
    call get_environment_variable(name,length=n,status=rc)
    if(rc==1 .or. n==0)then
      text='';return
    end if
    if(rc/=0)error stop 'cannot inspect Krylov environment control'
    allocate(character(n)::text)
    call get_environment_variable(name,text,status=rc)
    if(rc/=0)error stop 'cannot read Krylov environment control'
  end subroutine

  subroutine read_krylov_controls(n,initial,seed,diagnostic_path)
    integer, intent(in) :: n
    complex(8), allocatable, intent(out) :: initial(:)
    integer, allocatable, intent(out) :: seed
    character(:), allocatable, intent(out) :: diagnostic_path
    character(:), allocatable :: text
    character(1) :: extra
    integer :: i,unit,rc,first
    real(8) :: re,im
    logical :: exists
    call environment_text('ED2_RNG_SEED',text)
    if(len(text)>0)then
      first=1
      if(text(1:1)=='-' .or. text(1:1)=='+')first=2
      if(first>len(text))error stop 'ED2_RNG_SEED must be an integer'
      do i=first,len(text)
        if(text(i:i)<'0' .or. text(i:i)>'9')error stop 'ED2_RNG_SEED must be an integer'
      end do
      allocate(seed);read(text,*,iostat=rc)seed
      if(rc/=0)error stop 'ED2_RNG_SEED is outside the default integer range'
    end if
    call environment_text('ED2_INITIAL_VECTOR',text)
    if(len(text)>0)then
      open(newunit=unit,file=text,status='old',action='read',iostat=rc)
      if(rc/=0)error stop 'cannot open ED2_INITIAL_VECTOR'
      allocate(initial(n))
      do i=1,n
        read(unit,*,iostat=rc) re,im
        if(rc/=0)error stop 'ED2_INITIAL_VECTOR needs n real/imaginary rows'
        initial(i)=cmplx(re,im,8)
      end do
      read(unit,*,iostat=rc)extra
      close(unit)
      if(rc>=0)error stop 'ED2_INITIAL_VECTOR has extra data'
    end if
    call environment_text('ED2_DIAGNOSTICS',diagnostic_path)
    if(len(diagnostic_path)>0)then
      inquire(file=diagnostic_path,exist=exists)
      if(exists)error stop 'ED2_DIAGNOSTICS file already exists'
    end if
  end subroutine

  subroutine write_krylov_diagnostics(path,d)
    character(*), intent(in) :: path
    type(krylov_diagnostics), intent(in) :: d
    integer :: unit,rc,j
    open(newunit=unit,file=path,status='new',action='write',iostat=rc)
    if(rc/=0)error stop 'cannot create ED2_DIAGNOSTICS file'
    write(unit,'(a)') 'schema=1 (internal telemetry, not stable API)'
    write(unit,'(a,i0)') 'algorithm=',d%algorithm
    write(unit,'(a,i0)') 'dimension=',d%dimension
    write(unit,'(a,i0)') 'wanted=',d%wanted
    write(unit,'(a,i0)') 'status=',d%status
    write(unit,'(a,a)') 'reason=',trim(d%reason)
    write(unit,'(a,a)') 'termination=',trim(d%termination)
    write(unit,'(a,a)') 'lowest_state_inclusion=',trim(d%inclusion)
    write(unit,'(a,i0)') 'iterations=',d%iterations
    write(unit,'(a,i0)') 'cycles=',d%cycles
    write(unit,'(a,i0)') 'h_applications=',d%h_applications
    write(unit,'(a,i0)') 'requested_budget=',d%requested_budget
    write(unit,'(a,i0)') 'effective_budget=',d%effective_budget
    write(unit,'(a,i0)') 'miniter=',d%miniter
    write(unit,'(a,i0)') 'interval=',d%interval
    write(unit,'(a,i0)') 'keep=',d%keep
    write(unit,'(a,i0)') 'capacity=',d%capacity
    write(unit,'(a,es26.16e3)') 'energy_tolerance=',d%energy_tolerance
    write(unit,'(a,l1)') 'scale_available=',d%scale_available
    write(unit,'(a,es26.16e3)') 'operator_scale=',d%operator_scale
    write(unit,'(a,a)') 'scale_provenance=',trim(d%scale_provenance)
    write(unit,'(a,es26.16e3)') 'breakdown_threshold=',d%breakdown_threshold
    write(unit,'(a,a)') 'candidate_role=',trim(d%candidate_role)
    write(unit,'(a,l1)') 'energy_gate_satisfied=',d%energy_gate_satisfied
    write(unit,'(a,i0)') 'energy_consecutive_checks=',d%energy_consecutive_checks
    write(unit,'(a,l1)') 'energy_change_available=',d%energy_change_available
    write(unit,'(a,es26.16e3)') 'last_energy_change=',d%energy_change
    write(unit,'(a,i0)') 'exact_closures=',d%exact_closures
    write(unit,'(a,i0)') 'near_closures=',d%near_closures
    write(unit,'(a,es26.16e3)') 'last_beta=',d%last_beta
    write(unit,'(a,es26.16e3)') 'closure=',d%closure
    write(unit,'(a,l1)') 'orthogonality_available=',d%orthogonality_available
    write(unit,'(a,es26.16e3)') 'max_abs_gram_minus_identity=',d%orthogonality_error
    write(unit,'(a,l1)') 'roundoff_risk=',d%roundoff_risk
    write(unit,'(a,l1)') 'roundoff_limit_proven=',d%roundoff_limit_proven
    write(unit,'(a,a)') 'roundoff_assessment=',trim(d%roundoff_assessment)
    write(unit,'(a,l1)') 'explicit_initial=',d%explicit_initial
    write(unit,'(a,l1)') 'seed_specified=',d%seed_specified
    write(unit,'(a,i0)') 'seed=',d%seed
    write(unit,'(a,i0)') 'random_vectors=',d%random_vectors
    write(unit,'(a,i0)') 'configured_max_threads=',d%max_threads
    write(unit,'(a,i0)') 'rng_active_threads=',d%rng_active_threads
    write(unit,'(a,a)') 'rng_implementation=',trim(d%rng_implementation)
    if(allocated(d%rng_initial_state))then
      do j=1,d%rng_active_threads
        write(unit,'(a,i0,a,*(1x,i0))') 'rng_initial_state_',j-1,'=',d%rng_initial_state(:,j)
      end do
    end if
    write(unit,'(a)') 'rng_mapping=seed+104729*tid+37*k+1000003*draw_index modulo (huge(integer)-1)'
    write(unit,'(a)') 'pair_columns=index available residual_certified Ritz Rayleigh residual norm residual_threshold'
    if(allocated(d%ritz))then
      do j=1,size(d%ritz)
        write(unit,'(i0,2(1x,l1),5(1x,es26.16e3))')j,d%pair_available(j),d%residual_certified(j), &
          d%ritz(j),d%rayleigh(j),d%residual(j),d%norm(j),d%residual_threshold(j)
      end do
    end if
    if(d%residual_policy_enabled)then
      write(unit,'(a,l1)') 'residual_policy_enabled=',d%residual_policy_enabled
      write(unit,'(a,es26.16e3)') 'requested_residual_atol=',d%requested_atol
      write(unit,'(a,es26.16e3)') 'requested_residual_rtol=',d%requested_rtol
      write(unit,'(a,l1)') 'requested_scale_provided=',d%requested_scale_provided
      write(unit,'(a,es26.16e3)') 'requested_operator_scale=',d%requested_operator_scale
      write(unit,'(a,es26.16e3)') 'requested_energy_origin=',d%requested_origin
      write(unit,'(a)') 'requested_scale_source=caller supplied H-origin*I scale; not independently verified'
      write(unit,'(a)') 'requested_scale_rule=rtol>0: max(supplied scale,abs(Ritz-origin)); rtol=0: unused (0)'
      write(unit,'(a,l1)') 'all_requested_residuals_met=',d%all_requested_residuals_met
      if(allocated(d%requested_met))then
        do j=1,size(d%requested_met)
          write(unit,'(a,i0,a,es26.16e3)') 'requested_residual_scale_',j,'=',d%requested_scale(j)
          write(unit,'(a,i0,a,es26.16e3)') 'requested_residual_threshold_',j,'=',d%requested_threshold(j)
          write(unit,'(a,i0,a,l1)') 'requested_residual_met_',j,'=',d%requested_met(j)
          write(unit,'(a,i0,a,l1)') 'combined_residual_certified_',j,'=',d%combined_certified(j)
        end do
      end if
    end if
    close(unit,iostat=rc)
    if(rc/=0)error stop 'cannot close ED2_DIAGNOSTICS file'
  end subroutine
end module
