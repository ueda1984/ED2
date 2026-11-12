! Internal numerical kernels and test seam, not a stable library API.
module krylov_solvers
  use iso_fortran_env, only: int64
  use blas_lapack_interfaces, only: dznrm2, zdotc, dstevx, zheevr, zgemm, dlamch, &
    check_lapack_info, zheevr_workspace_lengths
  use krylov_observations, only: krylov_diagnostics,begin_observation,finish_observation,observe_pair,observe_gram
  use krylov_residual_policy, only: residual_policy,valid_residual_policy,requested_residual_limit
  use krylov_observations, only: record_residual_policy,observe_requested_residual
  !$ use omp_lib, only: omp_get_thread_num,omp_get_num_threads
  implicit none
  private
  public :: lanczos_ground, restarted_hermitian, energy_change, residual_limit, krylov_diagnostics
  public :: finite_real, normalize_vector, residual_policy
  abstract interface
    subroutine operator_action(n,x,y)
      integer, intent(in) :: n
      complex(8), intent(in) :: x(n)
      complex(8), intent(out) :: y(n)
    end subroutine operator_action
  end interface
contains
  logical function finite_real(x)
    real(8), intent(in) :: x
    ! Keep the check effective under GNU -Ofast (IEEE binary64 platforms).
    finite_real = iand(transfer(x,0_int64),int(z'7ff0000000000000',int64)) /= int(z'7ff0000000000000',int64)
  end function finite_real

  logical function normalize_vector(x)
    complex(8), intent(inout) :: x(:)
    real(8) :: scale, norm
    integer :: i
    normalize_vector=.false.
    do i=1,size(x)
      if(.not.finite_real(real(x(i),8)) .or. .not.finite_real(aimag(x(i)))) return
    end do
    scale=max(maxval(abs(real(x,8))),maxval(abs(aimag(x))))
    if(scale<=0d0) return
    x=x/scale
    norm=dznrm2(size(x),x,1)
    if(.not.finite_real(norm) .or. norm<=0d0) return
    x=x/norm
    normalize_vector=.true.
  end function normalize_vector

  real(8) function energy_change(old,new,scale) result(change)
    real(8), intent(in) :: old,new,scale
    if(abs(new)>sqrt(epsilon(1d0))*max(scale,tiny(1d0)))then
      change=abs(old/new-1d0) ! The legacy relative-energy criterion away from zero.
    else
      change=abs(old-new)/max(scale,tiny(1d0))
    end if
  end function energy_change

  real(8) function residual_limit(tol,scale,value) result(limit)
    real(8), intent(in) :: tol,scale,value
    ! Energy errors are second order in eigenvector errors. This is an
    ! additional certificate, not a replacement for the energy-change gate.
    limit=max(256d0*epsilon(1d0),sqrt(tol))*max(scale,abs(value),tiny(1d0))
  end function residual_limit

  real(8) function near_limit(scale) result(limit)
    real(8), intent(in) :: scale
    limit=64d0*epsilon(1d0)*max(scale,tiny(1d0))
  end function near_limit

  subroutine random_start(x,seed,draw_index,diagnostics)
    complex(8), intent(out) :: x(:)
    integer, intent(in), optional :: seed
    integer, intent(inout) :: draw_index
    type(krylov_diagnostics), intent(inout), optional :: diagnostics
    integer :: i,k,tid,nseed
    integer, allocatable :: saved(:),state(:)
    real(8) :: v(2)
    if(present(seed))then
      call random_seed(size=nseed)
      if(present(diagnostics))then
        if(.not.allocated(diagnostics%rng_initial_state))then
          allocate(diagnostics%rng_initial_state(nseed,diagnostics%max_threads))
          diagnostics%rng_initial_state=0
        end if
      end if
      ! An explicitly seeded draw preserves each participating thread's ambient
      ! intrinsic RNG state. The default branch below is the original loop.
      !$omp parallel private(saved,state,tid,k,i,v)
      allocate(saved(nseed),state(nseed))
      tid=0
      !$ tid=omp_get_thread_num()
      call random_seed(get=saved)
      do k=1,nseed
        state(k)=int(modulo(int(seed,int64)+104729_int64*tid+37_int64*k+ &
                           1000003_int64*draw_index,int(huge(1)-1,int64)))
      end do
      if(present(diagnostics))then
        !$omp single
        diagnostics%rng_active_threads=1
        !$ diagnostics%rng_active_threads=omp_get_num_threads()
        !$omp end single
        if(draw_index==0)diagnostics%rng_initial_state(:,tid+1)=state
      end if
      call random_seed(put=state)
      !$omp do schedule(static)
      do i=1,size(x)
        call random_number(v)
        x(i)=cmplx(v(1)-0.5d0,v(2)-0.5d0,8)
      end do
      !$omp end do
      call random_seed(put=saved)
      deallocate(saved,state)
      !$omp end parallel
      draw_index=draw_index+1
      return
    end if
    !$omp parallel do private(i,v)
    do i=1,size(x)
      call random_number(v)
      x(i)=cmplx(v(1)-0.5d0,v(2)-0.5d0,8)
    end do
    draw_index=draw_index+1
  end subroutine random_start

  ! Optional trailing arguments preserve every existing call site.
  subroutine lanczos_ground(apply,n,tol,maxiter,miniter,interval,value,vector, &
                            iterations,status,reason,residual,initial,diagnostics,rng_seed,policy)
    procedure(operator_action) :: apply
    integer, intent(in) :: n,maxiter,miniter,interval
    real(8), intent(in) :: tol
    real(8), intent(out) :: value,residual
    complex(8), intent(out) :: vector(n)
    integer, intent(out) :: iterations,status
    character(*), intent(out) :: reason
    complex(8), intent(in), optional :: initial(n)
    type(krylov_diagnostics), intent(out), optional :: diagnostics
    integer, intent(in), optional :: rng_seed
    type(residual_policy), intent(in), optional :: policy
    if(present(diagnostics))call begin_observation(diagnostics,1,n,1,maxiter,tol,present(initial),rng_seed)
    if(present(policy))then
      if(present(diagnostics))call record_residual_policy(diagnostics,policy)
      if(.not.valid_residual_policy(policy,reason))then
        status=-1;iterations=0;value=0d0;residual=huge(1d0);vector=(0d0,0d0)
        if(present(diagnostics))call finish_observation(diagnostics,status,reason,iterations)
        return
      end if
    end if
    call lanczos_ground_core(apply,n,tol,maxiter,miniter,interval,value,vector, &
                             iterations,status,reason,residual,initial,diagnostics,rng_seed,policy)
    if(present(diagnostics))call finish_observation(diagnostics,status,reason,iterations)
  end subroutine

  subroutine restarted_hermitian(apply,n,wanted,keep,capacity,tol,maxcycles,values,vectors,iterations, &
                                status,reason,residuals,initial,diagnostics,rng_seed,policy)
    procedure(operator_action) :: apply
    integer, intent(in) :: n,wanted,keep,capacity,maxcycles
    real(8), intent(in) :: tol
    real(8), intent(out) :: values(:),residuals(wanted)
    complex(8), intent(out) :: vectors(n,wanted)
    integer, intent(out) :: iterations,status
    character(*), intent(out) :: reason
    complex(8), intent(in), optional :: initial(n)
    type(krylov_diagnostics), intent(out), optional :: diagnostics
    integer, intent(in), optional :: rng_seed
    type(residual_policy), intent(in), optional :: policy
    if(present(diagnostics))call begin_observation(diagnostics,2,n,wanted,maxcycles,tol,present(initial),rng_seed)
    if(present(policy))then
      if(present(diagnostics))call record_residual_policy(diagnostics,policy)
      if(.not.valid_residual_policy(policy,reason))then
        status=-1;iterations=0;values=0d0;residuals=huge(1d0);vectors=(0d0,0d0)
        if(present(diagnostics))call finish_observation(diagnostics,status,reason,iterations)
        return
      end if
    end if
    call restarted_hermitian_core(apply,n,wanted,keep,capacity,tol,maxcycles,values,vectors,iterations, &
                                  status,reason,residuals,initial,diagnostics,rng_seed,policy)
    if(present(diagnostics))then
      call observe_gram(diagnostics,vectors)
      call finish_observation(diagnostics,status,reason,iterations)
    end if
  end subroutine

  subroutine lanczos_ground_core(apply,n,tol,maxiter,miniter,interval,value,vector, &
                                 iterations,status,reason,residual,initial,diagnostics,rng_seed,policy)
    procedure(operator_action) :: apply
    integer, intent(in) :: n,maxiter,miniter,interval
    real(8), intent(in) :: tol
    real(8), intent(out) :: value,residual
    complex(8), intent(out) :: vector(n)
    integer, intent(out) :: iterations,status
    character(*), intent(out) :: reason
    complex(8), intent(in), optional :: initial(n)
    type(krylov_diagnostics), intent(inout), optional :: diagnostics
    integer, intent(in), optional :: rng_seed
    type(residual_policy), intent(in), optional :: policy
    integer :: draw_index
    complex(8), allocatable :: q(:,:),start(:),r(:),candidate(:),action(:)
    real(8), allocatable :: alpha(:),beta(:),d(:),e(:),w(:),z(:,:),work(:)
    integer, allocatable :: iwork(:),ifail(:)
    integer :: k,limit,info,m
    real(8) :: scale,hnorm,previous,change,bnorm
    logical :: breakdown,have_previous,energy_ok,qualified,requested_ok
    real(8) :: requested_scale,requested_threshold
    draw_index=0
    status=1;reason='maximum iterations without convergence';iterations=0
    value=0d0;residual=huge(1d0);vector=(0d0,0d0)
    if(n<1 .or. maxiter<1 .or. miniter<0 .or. interval<1 .or. .not.finite_real(tol))then
      status=-1;reason='invalid Lanczos parameters';return
    end if
    if(tol<=0d0)then
      status=-1;reason='invalid Lanczos tolerance';return
    end if
    limit=min(n,maxiter)
    if(present(diagnostics))then
      diagnostics%effective_budget=limit;diagnostics%miniter=miniter;diagnostics%interval=interval
      diagnostics%scale_provenance='max observed ||Hq|| and |alpha| in first recurrence; estimate, not a spectral bound'
    end if
    allocate(q(n,3),start(n),r(n),candidate(n),action(n),alpha(limit),beta(0:limit))
    if(present(initial))then
      start=initial
    else
      call random_start(start,rng_seed,draw_index,diagnostics)
    end if
    if(.not.normalize_vector(start))then
      status=-1;reason='invalid initial vector';return
    end if
    q=(0d0,0d0);q(:,1)=start;beta=0d0;scale=0d0;have_previous=.false.;previous=0d0
    do k=1,limit
      iterations=k
      if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
      call apply(n,q(:,1),action)
      hnorm=dznrm2(n,action,1)
      alpha(k)=real(zdotc(n,q(:,1),1,action,1),8)
      if(.not.finite_real(hnorm) .or. .not.finite_real(alpha(k)))then
        status=2;reason='nonfinite operator action';return
      end if
      scale=max(scale,hnorm,abs(alpha(k)))
      if(present(diagnostics))then
        diagnostics%operator_scale=scale;diagnostics%breakdown_threshold=near_limit(scale)
        diagnostics%scale_available=.true.
        diagnostics%random_vectors=draw_index
      end if
      r=action-alpha(k)*q(:,1)
      if(k>1) r=r-beta(k-1)*q(:,3)
      bnorm=dznrm2(n,r,1);beta(k)=bnorm
      if(.not.finite_real(bnorm))then
        status=2;reason='nonfinite Lanczos beta';return
      end if
      breakdown=bnorm<=near_limit(scale)
      if(present(diagnostics))then
        diagnostics%last_beta=bnorm
        if(breakdown .and. bnorm==0d0) diagnostics%exact_closures=diagnostics%exact_closures+1
        if(breakdown .and. bnorm>0d0) diagnostics%near_closures=diagnostics%near_closures+1
      end if
      if(breakdown .or. k==limit .or. (k>=miniter .and. mod(k,interval)==0))then
        allocate(d(k),e(max(1,k-1)),w(k),z(k,1),work(5*k),iwork(5*k),ifail(k))
        if(k==1)then
          w(1)=alpha(1);z(1,1)=1d0
        else
          d=alpha(1:k);e=beta(1:k-1)
          call dstevx('V','I',k,d,e,0d0,0d0,1,1,2d0*dlamch('S'),m,w,z,k,work,iwork,ifail,info)
          call check_lapack_info('DSTEVX (Lanczos)',info)
        end if
        value=w(1);change=1d0
        if(have_previous) change=energy_change(previous,value,scale)
        energy_ok=have_previous .and. change<tol
        if(present(diagnostics))then
          diagnostics%energy_change=change;diagnostics%energy_change_available=have_previous
          diagnostics%energy_gate_satisfied=energy_ok
        end if
        write(*,'("itr, ene, conv: ",i5,es25.15,es11.1)') k,value,change
        if(breakdown .or. k==limit .or. energy_ok)then
          call reconstruct(k,z(:,1),candidate)
          if(.not.normalize_vector(candidate))then
            status=2;reason='invalid reconstructed eigenvector';return
          end if
          if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
          call apply(n,candidate,action)
          residual=dznrm2(n,action-value*candidate,1)
          qualified=finite_real(residual)
          if(qualified)then
            if(present(diagnostics)) call observe_pair(diagnostics,1,1,value,candidate,action, &
                                                       residual,residual_limit(tol,scale,value))
          end if
          if(qualified) qualified=residual<=residual_limit(tol,scale,value)
          requested_ok=.true.
          if(present(policy))then
            if(.not.requested_residual_limit(policy,value,requested_scale,requested_threshold))then
              status=2;reason='nonfinite requested residual scale or threshold';return
            end if
            requested_ok=finite_real(residual)
            if(requested_ok)requested_ok=residual<=requested_threshold
            if(present(diagnostics) .and. finite_real(residual)) &
              call observe_requested_residual(diagnostics,1,1,requested_scale,requested_threshold,requested_ok)
            qualified=qualified .and. requested_ok
            if(k==limit .and. .not.requested_ok)reason='maximum iterations without requested residual convergence'
          end if
          if(qualified .and. (energy_ok .or. breakdown .or. k==n))then
            vector=candidate;status=0
            if(abs(value)<=sqrt(epsilon(1d0))*max(scale,tiny(1d0))) &
              write(*,'(a)') 'Near-zero energy: use the absolute residual; the relative variance ratio is ill-conditioned.'
            if(breakdown)then
              if(bnorm==0d0)then
                reason='happy breakdown (exact zero beta; residual certified)'
              else
                reason='numerical near-breakdown (residual certified)'
              end if
            else if(k==n)then
              reason='full Krylov dimension (residual certified)'
            else
              reason='energy convergence (residual certified)'
            end if
            return
          end if
          if(breakdown)then
            status=2;reason='numerical breakdown without residual convergence'
            if(present(policy))then
              if(.not.requested_ok)reason='breakdown without requested residual convergence'
            end if
            return
          end if
        end if
        previous=value;have_previous=.true.
        deallocate(d,e,w,z,work,iwork,ifail)
      end if
      if(k<limit)then
        ! Never divide after exact/near breakdown or after the final iteration.
        q(:,3)=q(:,1);q(:,1)=r/bnorm
      end if
    end do
  contains
    subroutine reconstruct(count,coeff,out)
      integer, intent(in) :: count
      real(8), intent(in) :: coeff(count)
      complex(8), intent(out) :: out(n)
      complex(8), allocatable :: basis(:,:),hv(:)
      integer :: step
      allocate(basis(n,3),hv(n))
      basis=(0d0,0d0);basis(:,1)=start;out=coeff(1)*start
      do step=1,count-1
        if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
        call apply(n,basis(:,1),hv)
        hv=hv-alpha(step)*basis(:,1)
        if(step>1) hv=hv-beta(step-1)*basis(:,3)
        basis(:,3)=basis(:,1)
        ! All replayed betas precede the stopping beta and were checked above.
        basis(:,1)=hv/beta(step)
        out=out+coeff(step+1)*basis(:,1)
      end do
    end subroutine reconstruct
  end subroutine lanczos_ground_core

  subroutine orthogonalize(x,q,m)
    complex(8), intent(inout) :: x(:)
    complex(8), intent(in) :: q(:,:)
    integer, intent(in) :: m
    integer :: pass,j
    do pass=1,2
      do j=1,m
        x=x-zdotc(size(x),q(:,j),1,x,1)*q(:,j)
      end do
    end do
  end subroutine orthogonalize

  logical function complement(q,m,x) result(found)
    complex(8), intent(in) :: q(:,:)
    integer, intent(in) :: m
    complex(8), intent(out) :: x(:)
    real(8), allocatable :: leverage(:)
    integer :: row,attempt
    ! Maximal remaining coordinate leverage avoids repeatedly selecting a
    ! coordinate already almost contained in the retained invariant space.
    allocate(leverage(size(x)))
    leverage=1d0
    if(m>0) leverage=leverage-sum(abs(q(:,1:m))**2,dim=2)
    found=.false.
    do attempt=1,size(x)
      row=maxloc(leverage,dim=1);leverage(row)=-huge(1d0)
      x=(0d0,0d0);x(row)=(1d0,0d0)
      call orthogonalize(x,q,m)
      if(dznrm2(size(x),x,1)>64d0*epsilon(1d0))then
        found=normalize_vector(x);return
      end if
    end do
  end function complement

  subroutine restarted_hermitian_core(apply,n,wanted,keep,capacity,tol,maxcycles,values,vectors,iterations, &
                                       status,reason,residuals,initial,diagnostics,rng_seed,policy)
    procedure(operator_action) :: apply
    integer, intent(in) :: n,wanted,keep,capacity,maxcycles
    real(8), intent(in) :: tol
    real(8), intent(out) :: values(:),residuals(wanted)
    complex(8), intent(out) :: vectors(n,wanted)
    integer, intent(out) :: iterations,status
    character(*), intent(out) :: reason
    complex(8), intent(in), optional :: initial(n)
    type(krylov_diagnostics), intent(inout), optional :: diagnostics
    integer, intent(in), optional :: rng_seed
    type(residual_policy), intent(in), optional :: policy
    integer :: draw_index
    complex(8), allocatable :: q(:,:),hq(:,:),t(:,:),z(:,:),ritz(:,:),x(:),work(:)
    real(8), allocatable :: eig(:),rwork(:),previous(:)
    integer, allocatable :: support(:),iwork(:)
    complex(8) :: wquery(1)
    real(8) :: rquery(1),scale,bnorm,change,closure,hnorm
    integer :: iquery(1),info,lw,lrw,liw,count,cap,retained,m,cursor,cycle,j,flag,seeds
    logical :: found,certified,have_previous,requested_ok,requested_all
    real(8) :: requested_scale,requested_threshold
    draw_index=0
    status=1;reason='maximum restart iterations without convergence';iterations=0
    values=0d0;vectors=(0d0,0d0);residuals=huge(1d0)
    if(n<1 .or. wanted<1 .or. wanted>n .or. keep<wanted .or. capacity<keep .or. maxcycles<1)then
      status=-1;reason='invalid thick-restart dimensions or iteration budget';return
    end if
    if(.not.finite_real(tol) .or. tol<=0d0)then
      status=-1;reason='invalid thick-restart tolerance';return
    end if
    cap=min(capacity,n);retained=min(keep,cap)
    if(present(diagnostics))then
      diagnostics%effective_budget=maxcycles;diagnostics%keep=retained;diagnostics%capacity=cap
    end if
    if(size(values)<cap .or. (retained==cap .and. cap<n))then
      status=-1;reason='no expansion space for thick restart';return
    end if
    allocate(q(n,cap),hq(n,cap),t(cap,cap),z(cap,cap),ritz(n,cap),x(n), &
             eig(cap),support(2*cap),previous(retained))
    m=0;scale=0d0;have_previous=.false.;flag=0
    ! A block of wanted independent starts is necessary to resolve multiplicity:
    ! a single-vector Krylov space has at most one direction per exact eigenvalue.
    seeds=min(wanted,cap)
    ! A supplied exact excited-state start may have no ground-space component.
    ! Keep it in addition to the generic block, not instead of one direction.
    if(present(initial)) seeds=min(wanted+1,cap)
    do j=1,seeds
      if(j==1 .and. present(initial))then
        x=initial
      else
        call random_start(x,rng_seed,draw_index,diagnostics)
      end if
      if(j==1)then
        if(.not.normalize_vector(x))then
          status=-1;reason='invalid initial vector';return
        end if
      end if
      call orthogonalize(x,q,m)
      found=dznrm2(n,x,1)>64d0*epsilon(1d0)
      if(found)then
        found=normalize_vector(x)
      else
        found=complement(q,m,x)
      end if
      if(.not.found)then
        status=2;reason='cannot initialize orthogonal restart block';return
      end if
      m=m+1;q(:,m)=x
      if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
      call apply(n,x,hq(:,m));iterations=iterations+1
      hnorm=dznrm2(n,hq(:,m),1)
      if(.not.finite_real(hnorm))then
        status=2;reason='nonfinite operator action';return
      end if
      scale=max(scale,hnorm)
      if(present(diagnostics))then
        diagnostics%operator_scale=scale;diagnostics%breakdown_threshold=near_limit(scale)
        diagnostics%scale_available=.true.
        diagnostics%random_vectors=draw_index
      end if
    end do
    do cycle=1,maxcycles
      if(present(diagnostics))then
        diagnostics%cycles=cycle
        if(allocated(diagnostics%pair_available))diagnostics%pair_available=.false.
        diagnostics%orthogonality_available=.false.
      end if
      cursor=1
      do while(m<cap)
        found=.false.
        do while(cursor<=m)
          x=hq(:,cursor);cursor=cursor+1
          call orthogonalize(x,q,m)
          bnorm=dznrm2(n,x,1)
          if(.not.finite_real(bnorm))then
            status=2;reason='nonfinite restart beta';return
          end if
          if(present(diagnostics))diagnostics%last_beta=bnorm
          if(bnorm>near_limit(scale))then
            found=normalize_vector(x);exit
          end if
          if(present(diagnostics))then
            if(bnorm==0d0)diagnostics%exact_closures=diagnostics%exact_closures+1
            if(bnorm>0d0)diagnostics%near_closures=diagnostics%near_closures+1
          end if
          if(bnorm==0d0)then
            write(*,'(a,i0)') 'Exact Krylov block closure; orthogonal continuation after column ',cursor-1
          else
            write(*,'(a,es12.4)') 'Numerical near-breakdown; orthogonal continuation, beta = ',bnorm
          end if
        end do
        if(.not.found) found=complement(q,m,x)
        if(.not.found) exit
        m=m+1;q(:,m)=x
        if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
        call apply(n,x,hq(:,m));iterations=iterations+1
        hnorm=dznrm2(n,hq(:,m),1)
        if(.not.finite_real(hnorm))then
          status=2;reason='nonfinite operator action';return
        end if
        scale=max(scale,hnorm)
        if(present(diagnostics))then
          diagnostics%operator_scale=scale;diagnostics%breakdown_threshold=near_limit(scale)
        diagnostics%scale_available=.true.
          diagnostics%random_vectors=draw_index
        end if
      end do
      ! Compute the Hermitian Rayleigh-Ritz projection from actual H actions.
      ! This also handles closed blocks and complex restart couplings correctly.
      call zgemm('C','N',m,m,n,(1d0,0d0),q,n,hq,n,(0d0,0d0),t,cap)
      do j=1,m
        t(j,j)=cmplx(real(t(j,j),8),0d0,8)
      end do
      if(m==1)then
        ! The projected scalar problem is exact and needs no LAPACK query.
        eig(1)=real(t(1,1),8);z(1,1)=(1d0,0d0)
      else
        call zheevr('V','A','U',m,t,cap,0d0,0d0,1,m,dlamch('S'),count,eig,z,cap, &
                     support,wquery,-1,rquery,-1,iquery,-1,info)
        call check_lapack_info('ZHEEVR (TRLan workspace query)',info)
        call zheevr_workspace_lengths(m,wquery(1),rquery(1),iquery(1),lw,lrw,liw)
        allocate(work(lw),rwork(lrw),iwork(liw))
        call zheevr('V','A','U',m,t,cap,0d0,0d0,1,m,dlamch('S'),count,eig,z,cap, &
                     support,work,lw,rwork,lrw,iwork,liw,info)
        call check_lapack_info('ZHEEVR (TRLan solve)',info)
        deallocate(work,rwork,iwork)
      end if
      if(m<wanted)then
        status=2;reason='insufficient independent Ritz vectors';return
      end if
      call zgemm('N','N',n,m,m,(1d0,0d0),q,n,z,cap,(0d0,0d0),ritz,n)
      values(1:m)=eig(1:m)
      vectors=ritz(:,1:wanted)
      certified=.true.;requested_all=.true.
      do j=1,wanted
        if(.not.normalize_vector(vectors(:,j)))then
          status=2;reason='invalid Ritz vector';return
        end if
        if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
        call apply(n,vectors(:,j),x);iterations=iterations+1
        residuals(j)=dznrm2(n,x-eig(j)*vectors(:,j),1)
        if(.not.finite_real(residuals(j)))then
          status=2;reason='nonfinite Ritz residual';return
        end if
        if(present(diagnostics))call observe_pair(diagnostics,j,wanted,eig(j),vectors(:,j),x, &
                                                 residuals(j),residual_limit(tol,scale,eig(j)))
        certified=certified .and. residuals(j)<=residual_limit(tol,scale,eig(j))
        if(present(policy))then
          if(.not.requested_residual_limit(policy,eig(j),requested_scale,requested_threshold))then
            status=2;reason='nonfinite requested residual scale or threshold';return
          end if
          requested_ok=residuals(j)<=requested_threshold
          if(present(diagnostics))call observe_requested_residual(diagnostics,j,wanted, &
                                               requested_scale,requested_threshold,requested_ok)
          requested_all=requested_all .and. requested_ok
        end if
      end do
      if(present(policy))certified=certified .and. requested_all
      change=1d0
      if(have_previous)then
        change=0d0
        do j=1,min(retained,m)
          change=max(change,energy_change(previous(j),eig(j),scale))
        end do
      end if
      if(present(diagnostics))then
        diagnostics%energy_change=change;diagnostics%energy_change_available=have_previous
      end if
      if(have_previous .and. change<tol)then
        flag=flag+1
      else
        flag=0
      end if
      if(present(diagnostics))then
        diagnostics%energy_consecutive_checks=flag;diagnostics%energy_gate_satisfied=flag>=2
      end if
      write(*,'("tritr, ene(NOK), max_conv, flag: ",i5,es25.15,es11.1,i5)') cycle,eig(min(retained,m)),change,flag
      ! Check closure of the entire computed subspace, not just one beta.
      closure=0d0
      do j=1,m
        x=hq(:,j)
        call orthogonalize(x,q,m)
        closure=max(closure,dznrm2(n,x,1))
      end do
      if(present(diagnostics)) diagnostics%closure=closure
      if(certified .and. (flag>=2 .or. m==n .or. closure<=near_limit(scale)))then
        status=0
        if(m==n)then
          reason='full invariant space (residuals certified)'
        else if(closure==0d0)then
          reason='happy breakdown (closed invariant subspace; residuals certified)'
        else if(closure<=near_limit(scale))then
          reason='numerical near-breakdown (residuals certified)'
        else
          reason='two energy-converged restarts (residuals certified)'
        end if
        return
      end if
      if(cycle==maxcycles)then
        if(present(policy))then
          if(.not.requested_all)reason='maximum restart iterations without requested residual convergence'
        end if
        return
      end if
      previous=eig(1:retained);have_previous=.true.
      m=retained
      do j=1,m
        x=ritz(:,j)
        call orthogonalize(x,q,j-1)
        if(.not.normalize_vector(x))then
          status=2;reason='loss of retained-subspace independence';return
        end if
        q(:,j)=x
        if(present(diagnostics)) diagnostics%h_applications=diagnostics%h_applications+1_int64
        call apply(n,x,hq(:,j));iterations=iterations+1
      end do
    end do
  end subroutine restarted_hermitian_core
end module krylov_solvers
