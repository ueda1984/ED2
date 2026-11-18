module ham2vec
  implicit none
  complex(8), parameter :: zu = (1.0d0,0.0d0), zi = (0.0d0,1.0d0), zz = (0.0d0,0.0d0)
  real(8), parameter :: eps = 1.0d-15
  integer, allocatable :: wk_loc(:,:)
  complex(8), allocatable :: wk_ele(:,:)
  !
contains

  subroutine make_wk_loc_and_ele(lv)
    use state_lists, only: j_flip_ni, list_fly, a_down_spin, &
      c_down_spin, a2_down_spin, c2_down_spin
    use input_param, only: NODmax,NODmin,NO_one,NO_two,p_one,p_two,list_s,list_r,explist, &
      Jint,hvec,NOS,local_NODmax,local_spin,wk_dim,MNTE
    use cubic_anchored_lookup, only: setup_cubic_anchored_lookup,resolve_cached_hop
    use omp_lib, only: omp_get_wtime
    implicit none
    integer, intent(in)::lv
    integer :: s, id, wk_dim1
    integer :: n0, n1, n2, i1, i2, dN
    integer :: i, j, ell(6), jd, max_jd
    integer, allocatable :: ni(:), st_list(:)
    complex(8) :: c, phase
    character(len=16) :: local_bond_setting
    integer :: local_bond_status
    logical :: use_local_bonds, cache_profile
    character(len=16) :: profile_setting
    integer :: profile_status
    real(8) :: lookup_begin
    !
    local_bond_setting=''
    call get_environment_variable('ED2_CUBIC_LOCAL_BOND_CACHE', &
      local_bond_setting,status=local_bond_status)
    use_local_bonds=local_bond_status==0.and. &
      trim(adjustl(local_bond_setting))=='1'
    profile_setting=''
    call get_environment_variable('ED2_CUBIC_CACHE_PROFILE',profile_setting,status=profile_status)
    if(profile_status==0)then
      if(trim(adjustl(profile_setting))/='0'.and.trim(adjustl(profile_setting))/='1') &
        error stop 'ED2_CUBIC_CACHE_PROFILE must be 0 or 1'
    end if
    cache_profile=profile_status==0.and.trim(adjustl(profile_setting))=='1'
    if(cache_profile.and..not.use_local_bonds) &
      error stop 'Cubic cache profile requires ED2_CUBIC_LOCAL_BOND_CACHE=1'
    if(wk_dim < 0)then
      wk_dim = lv
    else
      wk_dim = min(wk_dim,lv)
    end if
    !
    if(MNTE < 0)then
      wk_dim1= 2*NO_one+8*NO_two+1
    else
      wk_dim1= MNTE + 1
    end if
    !
    if(cache_profile)lookup_begin=omp_get_wtime()
    call setup_cubic_anchored_lookup(lv)
    if(cache_profile)write(*,'(a,es18.9)') &
      'Cubic cache profile: representative_key_setup_wall_s=',omp_get_wtime()-lookup_begin
    allocate(wk_loc(wk_dim1,wk_dim))
    allocate(wk_ele(wk_dim1,wk_dim))
    if(use_local_bonds)then
      call build_cubic_local_bond_cache(lv,wk_dim1)
      call dump_cubic_cache_if_requested()
      return
    end if
    !
    allocate(ni(NODmax),st_list(NODmax))
    dN = NODmax-NODmin
    !
    max_jd = 1
    !
    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,st_list,c,jd,phase) reduction(max:max_jd)
    do i = 1, wk_dim
      wk_loc(1,i)  = i
      wk_loc(2:,i) = 0
      wk_ele(:,i) = zz
      jd = 2
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0 )
      !
      do j = 1, NO_one
        ! j selects the field record; i1 selects the site acted on.
        i1 = p_one(j)
        n1 = count(st_list == i1)
        wk_ele(1,i) = wk_ele(1,i) + hvec(3,j) * (local_spin(i1)-n1)
        ! - term
        c = 0.5d0 * ( hvec(1,j) + zi * hvec(2,j) )
        if(n0 < dN .and. n1>0 .and. abs(c) > eps)then
          ni = a_down_spin(i1,st_list,NODmax)
                    call resolve_cached_hop(ni,lv,id,phase)
          if(id > 0)then
            if(jd > wk_dim1) stop "MNTE is too small!!"
            wk_loc(jd,i) = id
            wk_ele(jd,i) = c &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1) * list_r(id) / list_r(i) &
              * phase
            jd = jd + 1
          end if
        end if
        ! + term
        c = 0.5d0 * ( hvec(1,j) - zi * hvec(2,j) )
        if( n0>0 .and. n1<local_NODmax(i1) .and. abs(c) > eps )then
          ni = c_down_spin(i1,st_list,NODmax)
                    call resolve_cached_hop(ni,lv,id,phase)
          if(id > 0)then
            if(jd > wk_dim1) stop "MNTE is too small!!"
            wk_loc(jd,i) = id
            wk_ele(jd,i) = c &
              * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)) &
              * list_r(id) / list_r(i) &
              * phase
            jd = jd + 1
          end if
        end if
      end do
      !
      do j = 1, NO_two
        i1 = p_two(1,j)
        i2 = p_two(2,j)
        n1 = count(st_list == i1)
        n2 = count(st_list == i2)
        ! zz term
        wk_ele(1,i) = wk_ele(1,i) + Jint(3,j) * &
          (local_spin(i1)-n1)*(local_spin(i2)-n2)
        if(n1>0)then
          ! -z term
          c = 0.5d0 * ( zi * (Jint(4,j)+Jint(7,j)) - Jint(5,j) + Jint(8,j) )
          if(n0 < dN .and. abs(c) > eps)then
            ni = a_down_spin(i1,st_list,NODmax)
                        call resolve_cached_hop(ni,lv,id,phase)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1)*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * phase
              jd = jd + 1
            end if
          end if
          ! -- term
          c = 0.25d0*(Jint(1,j) - Jint(2,j)) + zi * 0.5d0 * Jint(9,j)
          if(n0+1<dN .and. n2>0 .and. abs(c) > eps)then
            ni = a2_down_spin(i1,i2,st_list,NODmax)
                        call resolve_cached_hop(ni,lv,id,phase)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(2.0d0*local_spin(i2)-n2+1)*n2) &
                * list_r(id) / list_r(i) &
                * phase
              jd = jd + 1
            end if
          end if
          ! -+ term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) - zi*0.5d0*Jint(6,j)
          if(n2<local_NODmax(i2) .and. abs(c) > eps)then
            ni = j_flip_ni(i1,i2,st_list,NODmax)
                        call resolve_cached_hop(ni,lv,id,phase)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(n2+1)*(2.0d0*local_spin(i2)-n2) ) &
                * list_r(id) / list_r(i) &
                * phase
              jd = jd + 1
            end if
          end if
        end if
        if(n2>0)then
          ! z- term
          c = 0.5d0 * ( -zi * (Jint(4,j) - Jint(7,j)) + Jint(5,j) + Jint(8,j) )
          if(n0<dN .and. abs(c) > eps)then
            ni = a_down_spin(i2,st_list,NODmax)
                        call resolve_cached_hop(ni,lv,id,phase)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2)*(local_spin(i1)-n1) &
                * list_r(id) / list_r(i) &
                * phase
              jd = jd + 1
            end if
          end if
          ! +- term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) + zi*0.5d0*Jint(6,j)
          if(n1<local_NODmax(i1) .and. abs(c) > eps)then
            ni = j_flip_ni(i2,i1,st_list,NODmax)
                        call resolve_cached_hop(ni,lv,id,phase)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2*(n1+1)*(2.0d0*local_spin(i1)-n1) ) &
                * list_r(id) / list_r(i) &
                * phase
              jd = jd + 1
            end if
          end if
        end if
        if( n0>0 )then
          ! +z term
          c = 0.5d0 * ( -zi * ( Jint(4,j) + Jint(7,j) ) - Jint(5,j) + Jint(8,j) )
          if( n1<local_NODmax(i1) .and. abs(c) > eps )then
            ni = c_down_spin(i1,st_list,NODmax)
                        call resolve_cached_hop(ni,lv,id,phase)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1))*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * phase
              jd = jd + 1
            end if
          end if
          if( n2<local_NODmax(i2) )then
            ! z+ term
            c = 0.5d0 * ( zi * ( Jint(4,j) - Jint(7,j) ) + Jint(5,j) + Jint(8,j) )
            if( abs(c) > eps )then
              ni = c_down_spin(i2,st_list,NODmax)
                            call resolve_cached_hop(ni,lv,id,phase)
              if(id > 0)then
                if(jd > wk_dim1) stop "MNTE is too small!!"
                wk_loc(jd,i) = id
                wk_ele(jd,i) = c &
                  * sqrt((n2+1)*(2.0d0*local_spin(i2)-n2))*(local_spin(i1)-n1) &
                  * list_r(id) / list_r(i) &
                  * phase
                jd = jd + 1
              end if
            end if
            ! ++ term
            c = 0.25d0 * ( Jint(1,j) - Jint(2,j) ) - zi * 0.5d0 * Jint(9,j)
            if( n1<local_NODmax(i1) .and. n0>1 .and. abs(c) > eps )then
              ni = c2_down_spin(i1,i2,st_list,NODmax)
                            call resolve_cached_hop(ni,lv,id,phase)
              if(id > 0)then
                if(jd > wk_dim1) stop "MNTE is too small!!"
                wk_loc(jd,i) = id
                wk_ele(jd,i) = c &
                  * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)*(n2+1)*(2.0d0*local_spin(i2)-n2)) &
                  * list_r(id) / list_r(i) &
                  * phase
                jd = jd + 1
              end if
            end if
          end if
        end if
        !
      end do
      !
      if(jd > max_jd) max_jd = jd
    end do

    write(*,*) "Current MNTE = ", wk_dim1 - 1
    write(*,*) "Optimal MNTE = ", max_jd - 1

    call dump_cubic_cache_if_requested()
    return
  end subroutine make_wk_loc_and_ele


  ! Exactly the fixed-Ndown=3 S=1/2 real XXZ fast path, off by default.
  ! Any generic ED2 case still runs make_wk_loc_and_ele above unchanged.
  !
  ! For this model all transverse operators exchange exactly one occupied
  ! (down) endpoint and one empty (up) endpoint. Each of 3 down spins has
  ! exactly 18 incident input bonds (6 nearest and 12 face diagonal).
  ! The 9000 interaction records can therefore be replaced, per basis state,
  ! by at most 54 sorted/uniqued *actual input bond record IDs*.
  !
  ! The diagonal follows exactly from spin half:
  !   (1/2-n_i)(1/2-n_j) = 1/4 - (n_i+n_j)/2 + n_i*n_j.
  ! The all-up contribution E_ref=Sum Jz/4 is computed once in original
  ! record order. Each bond with exactly one down endpoint adds -Jz/2;
  ! zero/two down endpoints yield no correction. This regrouping changes
  ! floating-point addition order ONLY for the diagonal (bounded numerically
  ! and tested against every original cache diagonal). Each off-diagonal
  ! cache ID, order and real/imaginary IEEE bit pattern is unchanged.
  subroutine build_cubic_local_bond_cache(lv,wk_dim1)
    use iso_fortran_env, only: real64,int64
    use omp_lib, only: omp_get_wtime,omp_get_thread_num,omp_get_max_threads
    use input_param, only: NOS,NODmin,NODmax,NO_one,NO_two,local_NODmax,local_spin, &
      p_two,Jint,list_s,list_r,wk_dim
    use state_lists, only: list_fly,insertion_sort,j_flip_ni
    use cubic_anchored_lookup, only: resolve_cached_hop
    implicit none
    integer, intent(in) :: lv,wk_dim1
    integer, parameter :: degree_expected=18
    integer, allocatable :: degree(:),neighbors(:,:)
    integer :: site,j,i,i1,i2,n1,n2,jd,id,k,nel,old_j,max_jd,nt
    integer :: spins(3),moved(3),incident(3*degree_expected)
    complex(8) :: c,phase
    real(real64) :: e_reference
    character(len=16) :: setting
    integer :: env_status
    integer :: n_sites,prof_status,thread_id,num_threads
    integer(int64), allocatable :: prof_counts(:,:)
    real(real64), allocatable :: prof_times(:,:)
    real(real64) :: t0,sampled_start,t_setup,t_columns,profile_start
    real(real64) :: sample_estimate
    character(len=16) :: prof_setting
    logical :: profile,do_sample
    prof_setting=''
    call get_environment_variable('ED2_CUBIC_CACHE_PROFILE',prof_setting,status=prof_status)
    profile=prof_status==0.and.trim(adjustl(prof_setting))=='1'
    if(profile)then
      num_threads=max(1,omp_get_max_threads())
      allocate(prof_counts(6,0:num_threads-1),prof_times(4,0:num_threads-1))
      prof_counts=0_int64
      prof_times=0.0_real64
      profile_start=omp_get_wtime()
    end if
    setting=''
    call get_environment_variable('ED2_CUBIC_ANCHORED_REP',setting,status=env_status)
    if(env_status/=0.or.trim(adjustl(setting))/='1') &
      error stop 'Local cubic cache requires validated ED2_CUBIC_ANCHORED_REP=1'
    setting=''
    call get_environment_variable('ED2_CUBIC_FAST_GROUP_VALIDATION',setting,status=env_status)
    if(env_status/=0.or.trim(adjustl(setting))/='1') &
      error stop 'Local cubic cache requires ED2_CUBIC_FAST_GROUP_VALIDATION=1'
    setting=''
    call get_environment_variable('ED2_CUBIC_ANCHORED_BASIS',setting,status=env_status)
    if(env_status/=0.or.trim(adjustl(setting))/='1') &
      error stop 'Local cubic cache requires ED2_CUBIC_ANCHORED_BASIS=1'
    if(NODmax/=3.or.NODmin/=3.or.NO_one/=0) &
      error stop 'Local cubic cache requires fixed three-magnon U(1) sector'
    if(NO_two/=9*NOS.or.NOS<64.or.wk_dim/=lv) &
      error stop 'Local cubic cache requires full cubic 9*N bonds and basis'
    if(any(local_NODmax/=1).or.any(local_spin/=0.5d0)) &
      error stop 'Local cubic cache requires spin half'
    ! Only real SxSx+SySy exchange and SzSz are supported. In particular,
    ! tests MUST reject transverse anisotropy, Dzyaloshinskii-Moriya terms
    ! and nonzero one-body fields rather than silently dropping terms.
    if(any(Jint(1,:)/=Jint(2,:)).or.any(Jint(4:9,:)/=0d0)) &
      error stop 'Local cubic cache requires exactly real XXZ bonds'
    allocate(degree(NOS),neighbors(degree_expected,NOS))
    degree=0
    neighbors=0
    e_reference=0d0
    do j=1,NO_two
      i1=p_two(1,j)
      i2=p_two(2,j)
      if(i1<1.or.i1>=i2.or.i2>NOS) &
        error stop 'Local cubic cache requires sorted distinct bond endpoints'
      if(degree(i1)>=degree_expected.or.degree(i2)>=degree_expected) &
        error stop 'Local cubic cache invalid bond degree'
      degree(i1)=degree(i1)+1
      degree(i2)=degree(i2)+1
      neighbors(degree(i1),i1)=j
      neighbors(degree(i2),i2)=j
      e_reference=e_reference+0.25d0*Jint(3,j)
    end do
    if(any(degree/=degree_expected))error stop 'Local cubic cache degree is not 18'
    if(wk_dim1<1)error stop 'Local cubic cache invalid capacity'
    if(profile)then
      t_setup=omp_get_wtime()-profile_start
      profile_start=omp_get_wtime()
    end if
    ! Work each projected state independently; the original input bond IDs
    ! are sorted before any matrix element, preserving legacy offdiag order.
    max_jd=1
    !$omp parallel do private(i,j,k,site,spins,moved,incident,nel,old_j,i1,i2,n1,n2,c,phase,jd,id, & 
    !$omp thread_id,t0,sampled_start,do_sample) &
    !$omp reduction(max:max_jd) schedule(static)
    do i=1,wk_dim
      if(profile)then
        thread_id=omp_get_thread_num()
        prof_counts(1,thread_id)=prof_counts(1,thread_id)+1_int64
        t0=omp_get_wtime()
      end if
      wk_loc(1,i)=i
      wk_loc(2:,i)=0
      wk_ele(:,i)=zz
      wk_ele(1,i)=cmplx(e_reference,0d0,kind=8)
      spins=list_fly(list_s(i),3,3,NOS)
      if(profile)then
        prof_times(1,thread_id)=prof_times(1,thread_id)+(omp_get_wtime()-t0)
        t0=omp_get_wtime()
      end if
      nel=0
      do k=1,3
        site=spins(k)
        incident(nel+1:nel+degree_expected)=neighbors(:,site)
        nel=nel+degree_expected
      end do
      call insertion_sort(incident,nel)
      if(profile)then
        prof_times(2,thread_id)=prof_times(2,thread_id)+(omp_get_wtime()-t0)
        t0=omp_get_wtime()
      end if
      old_j=0
      jd=2
      do k=1,nel
        j=incident(k)
        if(j==old_j)cycle
        if(profile)prof_counts(2,thread_id)=prof_counts(2,thread_id)+1_int64
        old_j=j
        i1=p_two(1,j)
        i2=p_two(2,j)
        n1=count(spins==i1)
        n2=count(spins==i2)
        if(n1==n2)cycle
        if(profile)prof_counts(3,thread_id)=prof_counts(3,thread_id)+1_int64
        wk_ele(1,i)=wk_ele(1,i)-0.5d0*Jint(3,j)
        if(n1==1)then
          c=0.25d0*(Jint(1,j)+Jint(2,j)) - zi*0.5d0*Jint(6,j)
          if(abs(c)<=eps)cycle
          moved=j_flip_ni(i1,i2,spins,3)
          if(profile)then
            prof_counts(4,thread_id)=prof_counts(4,thread_id)+1_int64
            do_sample=modulo(i+j,31)==0
            if(do_sample)sampled_start=omp_get_wtime()
          end if
          call resolve_cached_hop(moved,lv,id,phase)
          if(profile)then
            if(do_sample)then
              prof_times(4,thread_id)=prof_times(4,thread_id)+(omp_get_wtime()-sampled_start)
              prof_counts(5,thread_id)=prof_counts(5,thread_id)+1_int64
            end if
          end if
          if(id>0)then
            if(profile)prof_counts(6,thread_id)=prof_counts(6,thread_id)+1_int64
            if(jd>wk_dim1)error stop 'Local cubic cache MNTE too small'
            wk_loc(jd,i)=id
            wk_ele(jd,i)=c &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(n2+1)*(2.0d0*local_spin(i2)-n2)) &
              * list_r(id)/list_r(i)*phase
            jd=jd+1
          end if
        else
          c=0.25d0*(Jint(1,j)+Jint(2,j)) + zi*0.5d0*Jint(6,j)
          if(abs(c)<=eps)cycle
          moved=j_flip_ni(i2,i1,spins,3)
          if(profile)then
            prof_counts(4,thread_id)=prof_counts(4,thread_id)+1_int64
            do_sample=modulo(i+j,31)==0
            if(do_sample)sampled_start=omp_get_wtime()
          end if
          call resolve_cached_hop(moved,lv,id,phase)
          if(profile)then
            if(do_sample)then
              prof_times(4,thread_id)=prof_times(4,thread_id)+(omp_get_wtime()-sampled_start)
              prof_counts(5,thread_id)=prof_counts(5,thread_id)+1_int64
            end if
          end if
          if(id>0)then
            if(profile)prof_counts(6,thread_id)=prof_counts(6,thread_id)+1_int64
            if(jd>wk_dim1)error stop 'Local cubic cache MNTE too small'
            wk_loc(jd,i)=id
            wk_ele(jd,i)=c &
              * sqrt((2.0d0*local_spin(i2)-n2+1)*n2*(n1+1)*(2.0d0*local_spin(i1)-n1)) &
              * list_r(id)/list_r(i)*phase
            jd=jd+1
          end if
        end if
      end do
      if(profile)prof_times(3,thread_id)=prof_times(3,thread_id)+(omp_get_wtime()-t0)
      if(jd>max_jd)max_jd=jd
    end do
    !$omp end parallel do
    if(profile)then
      t_columns=omp_get_wtime()-profile_start
      sample_estimate=0.0_real64
      if(sum(prof_counts(5,:))>0_int64)sample_estimate= &
        sum(prof_times(4,:))*real(sum(prof_counts(4,:)),real64)/ &
        real(sum(prof_counts(5,:)),real64)
      write(*,'(a,i0,a,i0)') 'Cubic cache profile: columns=',sum(prof_counts(1,:)), &
        ', omp_max_threads=',num_threads
      write(*,'(a,es18.9)') 'Cubic cache profile: adjacency_setup_wall_s=',t_setup
      write(*,'(a,es18.9)') 'Cubic cache profile: column_parallel_wall_s=',t_columns
      write(*,'(a,es18.9)') 'Cubic cache profile: state_decode_thread_elapsed_sum_s=',sum(prof_times(1,:))
      write(*,'(a,es18.9)') 'Cubic cache profile: bond_union_sort_thread_elapsed_sum_s=',sum(prof_times(2,:))
      write(*,'(a,es18.9)') 'Cubic cache profile: bond_math_thread_elapsed_sum_s=',sum(prof_times(3,:))
      write(*,'(a,i0)') 'Cubic cache profile: unique_incident_bonds=',sum(prof_counts(2,:))
      write(*,'(a,i0)') 'Cubic cache profile: active_exchange_bonds=',sum(prof_counts(3,:))
      write(*,'(a,i0)') 'Cubic cache profile: resolver_calls=',sum(prof_counts(4,:))
      write(*,'(a,i0)') 'Cubic cache profile: resolver_sampled_calls=',sum(prof_counts(5,:))
      write(*,'(a,es18.9)') &
        'Cubic cache profile: resolver_sampled_thread_elapsed_sum_s=',sum(prof_times(4,:))
      write(*,'(a,es18.9)') &
        'Cubic cache profile: resolver_extrapolated_thread_elapsed_sum_s=',sample_estimate
      write(*,'(a,i0)') 'Cubic cache profile: stored_hops=',sum(prof_counts(6,:))
      write(*,'(a)') 'Cubic cache profile: resolver uses deterministic 1/31 bond-ID sampling;'
      write(*,'(a)') 'Cubic cache profile: thread elapsed sums are not parallel wall times;'
      write(*,'(a)') 'Cubic cache profile: probe overhead and sampling bias prohibit publication timing.'
    end if
    deallocate(degree,neighbors)
    write(*,'(a,i0,a)') 'Cubic cache builder: local-bond-sparse, degree=', &
      degree_expected,', sorted input-bond ordering retained'
    write(*,*) 'Current MNTE = ',wk_dim1-1
    write(*,*) 'Optimal MNTE = ',max_jd-1
  end subroutine build_cubic_local_bond_cache

  ! A test-only exact cache proof. Write all nonzero cache entries in the
  ! original per-column order, including the *bit patterns* of both complex
  ! components. The old and new offdiagonals must match byte-for-byte. The
  ! diagonal can differ by summation order but must agree within 1e-8.
  subroutine dump_cubic_cache_if_requested()
    use iso_fortran_env, only: int64,real64
    character(len=16) :: setting
    integer :: status,u,ios,i,j
    real(real64) :: re,im
    setting=''
    call get_environment_variable('ED2_CUBIC_CACHE_PROOF_DUMP',setting,status=status)
    if(status/=0.or.trim(adjustl(setting))/='1')return
    open(newunit=u,file='cubic_cache_proof_dump.txt',status='replace',action='write',iostat=ios)
    if(ios/=0)error stop 'Failed to open cubic cache proof dump'
    do i=1,size(wk_loc,2)
      do j=1,size(wk_loc,1)
        if(j>1.and.wk_loc(j,i)==0)cycle
        re=real(wk_ele(j,i),kind=real64)
        im=aimag(wk_ele(j,i))
        write(u,'(i0,1x,i0,1x,i0,1x,z16.16,1x,z16.16)',iostat=ios) &
          i,j,wk_loc(j,i),transfer(re,0_int64),transfer(im,0_int64)
        if(ios/=0)error stop 'Failed to write cubic cache proof dump'
      end do
    end do
    close(u,iostat=ios)
    if(ios/=0)error stop 'Failed to close cubic cache proof dump'
    write(*,'(a,i0)') 'Cubic cache proof dump: projected rows=',size(wk_loc,2)
  end subroutine dump_cubic_cache_if_requested

  subroutine make_full_hamiltonian(lv,Ham) 
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly, a_down_spin, &
      c_down_spin, a2_down_spin, c2_down_spin
    use input_param, only: NODmax,NODmin,NO_one,NO_two,p_one,p_two,list_s,list_r,explist, &
      Jint,hvec,NOS,local_NODmax,local_spin
    implicit none
    integer, intent(in)::lv
    complex(8),intent(inout)::Ham(lv,lv)
    integer::s,id
    integer :: n0, n1, n2, i1, i2, dN
    integer :: i, j, ell(6)
    integer, allocatable :: ni(:), st_list(:)
    complex(8) :: c
    allocate(ni(NODmax),st_list(NODmax))
    dN = NODmax-NODmin
    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,st_list,c)
    do i = 1, lv
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0 )
      !
      do j = 1, NO_one
        ! j selects the field record; i1 selects the site acted on.
        i1 = p_one(j)
        n1 = count(st_list == i1)
        Ham(i,i) = Ham(i,i) + hvec(3,j) * (local_spin(i1)-n1)
        ! - term
        c = 0.5d0 * ( hvec(1,j) + zi * hvec(2,j) )
        if(n1>0 .and. n0 < dN .and. abs(c) > eps)then
          ni = a_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            Ham(i,id) = Ham(i,id) + c &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1) &
              * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
          end if
        end if
        ! + term
        c = 0.5d0 * ( hvec(1,j) - zi * hvec(2,j) )
        if( n0 > 0 .and. n1<local_NODmax(i1) .and. abs(c) > eps )then
          ni = c_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            Ham(i,id) = Ham(i,id) + c &
              * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)) &
              * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
          end if
        end if
      end do
      !
      do j = 1, NO_two
        i1 = p_two(1,j)
        i2 = p_two(2,j)
        n1 = count(st_list == i1)
        n2 = count(st_list == i2)
        ! zz term
        Ham(i,i) = Ham(i,i) + Jint(3,j) * &
          (local_spin(i1)-n1)*(local_spin(i2)-n2)
        if(n1>0)then
          ! -z term
          c = 0.5d0 * ( zi * (Jint(4,j)+Jint(7,j)) - Jint(5,j) + Jint(8,j) )
          if(n0 < dN .and. abs(c) > eps)then
            ni = a_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              Ham(i,id) = Ham(i,id) + c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1)*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            end if
          end if
          ! -- term
          c = 0.25d0*(Jint(1,j) - Jint(2,j)) + zi * 0.5d0 * Jint(9,j)
          if(n0+1 < dN .and. n2>0 .and. abs(c) > eps)then
            ni = a2_down_spin(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              Ham(i,id) = Ham(i,id) + c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(2.0d0*local_spin(i2)-n2+1)*n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            end if
          end if
          ! -+ term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) - zi*0.5d0*Jint(6,j)
          if(n2<local_NODmax(i2) .and. abs(c) > eps)then
            ni = j_flip_ni(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              Ham(i,id) = Ham(i,id) + c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(n2+1)*(2.0d0*local_spin(i2)-n2) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            end if
          end if
        end if
        if(n2>0)then
          ! z- term
          c = 0.5d0 * ( -zi * (Jint(4,j) - Jint(7,j)) + Jint(5,j) + Jint(8,j) )
          if(n0 < dN .and. abs(c) > eps)then
            ni = a_down_spin(i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              Ham(i,id) = Ham(i,id) + c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2)*(local_spin(i1)-n1) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            end if
          end if
          ! +- term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) + zi*0.5d0*Jint(6,j)
          if(n1<local_NODmax(i1) .and. abs(c) > eps)then
            ni = j_flip_ni(i2,i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              Ham(i,id) = Ham(i,id) + c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2*(n1+1)*(2.0d0*local_spin(i1)-n1) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            end if
          end if
          !
        end if
        if( n0>0 )then
          ! +z term
          c = 0.5d0 * ( -zi * ( Jint(4,j) + Jint(7,j) ) - Jint(5,j) + Jint(8,j) )
          if( n1<local_NODmax(i1) .and. abs(c) > eps )then
            ni = c_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              Ham(i,id) = Ham(i,id) + c &
                * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1))*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            end if
          end if
          if( n2<local_NODmax(i2) )then
            ! z+ term
            c = 0.5d0 * ( zi * ( Jint(4,j) - Jint(7,j) ) + Jint(5,j) + Jint(8,j) )
            if( abs(c) > eps )then
              ni = c_down_spin(i2,st_list,NODmax)
              call representative_SQ(NODmax,s,ell,ni)
              call findstate(s,list_s,id,lv)
              if(id > 0)then
                Ham(i,id) = Ham(i,id) + c &
                  * sqrt((n2+1)*(2.0d0*local_spin(i2)-n2))*(local_spin(i1)-n1) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              end if
            end if
            ! ++ term
            c = 0.25d0 * ( Jint(1,j) - Jint(2,j) ) - zi * 0.5d0 * Jint(9,j)
            if( n1<local_NODmax(i1) .and. n0>1 .and. abs(c) > eps )then
              ni = c2_down_spin(i1,i2,st_list,NODmax)
              call representative_SQ(NODmax,s,ell,ni)
              call findstate(s,list_s,id,lv)
              if(id > 0)then
                Ham(i,id) = Ham(i,id) + c &
                  * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)*(n2+1)*(2.0d0*local_spin(i2)-n2)) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              end if
            end if
          end if
        end if
        !
      end do
    end do

    return
  end subroutine make_full_hamiltonian
  !
  subroutine ham_to_vec_wave_vector(v0,v1,lv,NODmax,NODmin,list_s,list_r,explist) 
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly,a_down_spin,c_down_spin,&
      a2_down_spin, c2_down_spin
    use input_param, only: NO_one,NO_two,p_one,p_two,Jint,hvec,NOS,local_NODmax,local_spin,&
      L1,L2,L3,L4,L5,L6,wk_dim
    integer, intent(in) :: lv
    complex(8), intent(in) :: v1(lv)
    complex(8), intent(out) :: v0(lv)
    integer, intent(in) :: NODmax,NODmin
    real(8), intent(in) :: list_r(1:lv)           
    integer, intent(in)::list_s(1:lv)          
    complex(8),intent(in)::explist(L1,L2,L3,L4,L5,L6)
    integer :: j, i, ell(6)
    integer, allocatable ::ni(:), st_list(:)
    integer :: s,id,wk_dim1
    integer :: n0, n1, n2, i1, i2, dN
    complex(8) :: c
    
    allocate(ni(NODmax),st_list(NODmax))
    dN = NODmax - NODmin

    wk_dim1= size(wk_ele,1)

    !$omp parallel do private(i,j)
    do i = 1, wk_dim
      v0(i) = 0.0d0
      do j = 1, wk_dim1
        if(wk_loc(j,i) == 0) exit
        v0(i) = v0(i) + wk_ele(j,i) * v1(wk_loc(j,i))
      end do
    end do
    !

    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,st_list,c)
    do i = wk_dim+1, lv
      v0(i) = (0.0d0, 0.0d0)
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0 )
      !
      do j = 1, NO_one
        ! j selects the field record; i1 selects the site acted on.
        i1 = p_one(j)
        n1 = count(st_list == i1)
        v0(i) = v0(i) + hvec(3,j) * (local_spin(i1)-n1) * v1(i)
        ! - term
        c = 0.5d0 * ( hvec(1,j) + zi * hvec(2,j) )
        if(n0 < dN .and. n1>0 .and. abs(c) > eps)then
          ni = a_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            v0(i) = v0(i) + c &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1) * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
          end if
        end if
        ! + term
        c = 0.5d0 * ( hvec(1,j) - zi * hvec(2,j) )
        if( n0>0 .and. n1<local_NODmax(i1) .and. abs(c) > eps )then
          ni = c_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            v0(i) = v0(i) + c &
              * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)) &
              * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
          end if
        end if
      end do
      !
      do j = 1, NO_two
        i1 = p_two(1,j)
        i2 = p_two(2,j)
        n1 = count(st_list == i1)
        n2 = count(st_list == i2)
        ! zz term
        v0(i) = v0(i) + Jint(3,j) * &
          (local_spin(i1)-n1)*(local_spin(i2)-n2) * v1(i)
        if(n1>0)then
          ! -z term
          c = 0.5d0 * ( zi * (Jint(4,j)+Jint(7,j)) - Jint(5,j) + Jint(8,j) )
          if(n0 < dN .and. abs(c) > eps)then
            ni = a_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              v0(i) = v0(i) + c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1)*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          ! -- term
          c = 0.25d0*(Jint(1,j) - Jint(2,j)) + zi * 0.5d0 * Jint(9,j)
          if(n0+1<dN .and. n2>0 .and. abs(c) > eps)then
            ni = a2_down_spin(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              v0(i) = v0(i) + c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(2.0d0*local_spin(i2)-n2+1)*n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          ! -+ term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) - zi*0.5d0*Jint(6,j)
          if(n2<local_NODmax(i2) .and. abs(c) > eps)then
            ni = j_flip_ni(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              v0(i) = v0(i) + c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(n2+1)*(2.0d0*local_spin(i2)-n2) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
        end if
        if(n2>0)then
          ! z- term
          c = 0.5d0 * ( -zi * (Jint(4,j) - Jint(7,j)) + Jint(5,j) + Jint(8,j) )
          if(n0<dN .and. abs(c) > eps)then
            ni = a_down_spin(i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              v0(i) = v0(i) + c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2)*(local_spin(i1)-n1) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          ! +- term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) + zi*0.5d0*Jint(6,j)
          if(n1<local_NODmax(i1) .and. abs(c) > eps)then
            ni = j_flip_ni(i2,i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              v0(i) = v0(i) + c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2*(n1+1)*(2.0d0*local_spin(i1)-n1) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          !
        end if
        if( n0>0 )then
          ! +z term
          c = 0.5d0 * ( -zi * ( Jint(4,j) + Jint(7,j) ) - Jint(5,j) + Jint(8,j) )
          if( n1<local_NODmax(i1) .and. abs(c) > eps )then
            ni = c_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              v0(i) = v0(i) + c &
                * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1))*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          if( n2<local_NODmax(i2) )then
            ! z+ term
            c = 0.5d0 * ( zi * ( Jint(4,j) - Jint(7,j) ) + Jint(5,j) + Jint(8,j) )
            if( abs(c) > eps )then
              ni = c_down_spin(i2,st_list,NODmax)
              call representative_SQ(NODmax,s,ell,ni)
              call findstate(s,list_s,id,lv)
              if(id > 0)then
                v0(i) = v0(i) + c &
                  * sqrt((n2+1)*(2.0d0*local_spin(i2)-n2))*(local_spin(i1)-n1) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
              end if
            end if
            ! ++ term
            c = 0.25d0 * ( Jint(1,j) - Jint(2,j) ) - zi * 0.5d0 * Jint(9,j)
            if( n1<local_NODmax(i1) .and. n0>1 .and. abs(c) > eps )then
              ni = c2_down_spin(i1,i2,st_list,NODmax)
              call representative_SQ(NODmax,s,ell,ni)
              call findstate(s,list_s,id,lv)
              if(id > 0)then
                v0(i) = v0(i) + c &
                  * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)*(n2+1)*(2.0d0*local_spin(i2)-n2)) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
              end if
            end if
          end if
        end if
        !
      end do
    end do

    return
  end subroutine ham_to_vec_wave_vector
  !
  subroutine calcu_lm(v1,lv,NODmax,NODmin,expe_lm,NOS,list_n,non)
    use state_lists, only: list_fly, a_down_spin, c_down_spin, representative_SQ, findstate
    use input_param, only: list_s, local_spin, explist, list_r, local_NODmax
    implicit none
    integer, intent(in) :: lv, NODmax, NODmin, NOS, non
    complex(8), intent(in) :: v1(lv)
    integer, intent(in) :: list_n(non)
    real(8), intent(inout) :: expe_lm(3,NOS)
    integer :: s, j, i, i1, id, ell(6)
    integer :: n0, n1, dN
    real(8) :: expval(3)
    integer, allocatable :: ni(:), st_list(:)
    complex(8) :: c
    allocate(ni(NODmax), st_list(NODmax))
    !
    expval = 0.0d0
    dN = NODmax-NODmin
    !
    !$omp parallel do private(n0,n1,i1,i,j,st_list,s,ni,id,ell,c) reduction(+:expval)
    do i = 1, lv
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0 )
      !
      do j = 1, non
        ! diagonal term
        i1 = list_n(j)
        n1 = count(st_list == i1)
        expval(3) = expval(3) + (local_spin(i1)-n1) * dble( conjg(v1(i))*v1(i) )
        ! - term
        if(n0<dN .and. n1>0)then
          ni = a_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            c = conjg(v1(i)) * 0.5d0 &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1) * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            expval(1) = expval(1) + dble(c)
            expval(2) = expval(2) - aimag(c)
          end if
        end if
        ! + term
        if( n0>0 .and. n1<local_NODmax(i1) )then
          ni = c_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            c = conjg(v1(i)) * 0.5d0 &
              * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)) &
              * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            expval(1) = expval(1) + dble(c)
            expval(2) = expval(2) + aimag(c)
          end if
        end if
        !
      end do
    end do
    !
    expval = expval / non
    !
    !$omp parallel do
    do i = 1, non
      expe_lm(1:3,list_n(i)) = expval(1:3)
    end do
    !
    return
  end subroutine calcu_lm
  !
  subroutine calcu_cf(v1,lv,NODmax,NODmin,expe_cf,p_two_cf)
    use state_lists, only: list_fly,a_down_spin,representative_SQ,findstate,a2_down_spin,&
      j_flip_ni,c_down_spin,c2_down_spin
    use input_param, only: local_spin,NOS,list_s,L1,L2,L3,explist,list_r,local_NODmax
    implicit none
    integer, intent(in) :: lv, NODmax, NODmin
    complex(8), intent(in) :: v1(lv)
    complex(8), intent(inout) :: expe_cf(3,3)
    integer, intent(in) :: p_two_cf(2)
    integer :: i1,i2,n0,n1,n2,s,i,id,ell(6),j,V,a,b,dN
    complex(8) :: expval(3,3)
    integer, allocatable :: ni(:), st_list(:), p_two(:,:)
    logical, allocatable :: fswap(:)
    V = L1*L2*L3
    allocate(ni(NODmax), st_list(NODmax),p_two(2,V),fswap(V))
    !
    expval = 0.0d0
    dN = NODmax-NODmin
    !
    call mk_p_two_and_fswap
    !
    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,st_list,a,b) reduction(+:expval)
    do i = 1, lv
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0)
      !
      do j = 1, V
        i1 = p_two(1,j)
        i2 = p_two(2,j)
        n1 = count(st_list == i1)
        n2 = count(st_list == i2)
        ! zz term
        expval(1,1) = expval(1,1) + conjg(v1(i)) * &
          (local_spin(i1)-n1)*(local_spin(i2)-n2) * v1(i)
        if(n1>0)then
          ! -z term
          if(n0<dN)then
            ni = a_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              call get_a_b(a,b,2,1,fswap(j))
              expval(a,b) = expval(a,b) + conjg(v1(i)) &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1)*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          ! -- term
          if(n0+1<dN .and. n2>0)then
            ni = a2_down_spin(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              expval(2,2) = expval(2,2) + conjg(v1(i)) &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(2.0d0*local_spin(i2)-n2+1)*n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          ! -+ term
          if(n2<local_NODmax(i2))then
            ni = j_flip_ni(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              call get_a_b(a,b,2,3,fswap(j))
              expval(a,b) = expval(a,b) + conjg(v1(i)) &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(n2+1)*(2.0d0*local_spin(i2)-n2) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
        end if
        if(n2>0)then
          ! z- term
          if(n0<dN)then
            ni = a_down_spin(i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              call get_a_b(a,b,1,2,fswap(j))
              expval(a,b) = expval(a,b) + conjg(v1(i)) &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2)*(local_spin(i1)-n1) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          ! +- term
          if(n1<local_NODmax(i1))then
            ni = j_flip_ni(i2,i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              call get_a_b(a,b,3,2,fswap(j))
              expval(a,b) = expval(a,b) + conjg(v1(i)) &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2*(n1+1)*(2.0d0*local_spin(i1)-n1) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          !
        end if
        if( n0>0 )then
          ! +z term
          if( n1<local_NODmax(i1) )then
            ni = c_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              call get_a_b(a,b,3,1,fswap(j))
              expval(a,b) = expval(a,b) + conjg(v1(i)) &
                * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1))*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
          end if
          if( n2<local_NODmax(i2) )then
            ! z+ term
            ni = c_down_spin(i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              call get_a_b(a,b,1,3,fswap(j))
              expval(a,b) = expval(a,b) + conjg(v1(i)) &
                * sqrt((n2+1)*(2.0d0*local_spin(i2)-n2))*(local_spin(i1)-n1) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
            end if
            ! ++ term
            if( n1<local_NODmax(i1) .and. n0>1 )then
              ni = c2_down_spin(i1,i2,st_list,NODmax)
              call representative_SQ(NODmax,s,ell,ni)
              call findstate(s,list_s,id,lv)
              if(id > 0)then
                expval(3,3) = expval(3,3) + conjg(v1(i)) &
                  * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)*(n2+1)*(2.0d0*local_spin(i2)-n2)) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6)) * v1(id)
              end if
            end if
          end if
        end if
        !
       end do
    end do
    !
    expe_cf = expval / V
    !
    return
    !
  contains
    !
    subroutine get_a_b(a,b,c,d,f)
      integer :: a, b, c, d
      logical :: f
      if(f)then
        a = d; b = c
      else
        a = c; b = d
      end if
    end subroutine get_a_b
    !
    subroutine mk_p_two_and_fswap
      use input_param, only: L1, L2, L3, shift_1, shift_2, shift_3
      implicit none
      integer :: i, j, k, t
      integer :: pp_two(2)
      t = 0
      pp_two(1:2) = p_two_cf(1:2)
      do k = 1, L3
        pp_two(:) = shift_3(pp_two(:))
        do j = 1, L2
          pp_two(:) = shift_2(pp_two(:))
          do i = 1, L1
            pp_two(:) = shift_1(pp_two(:))
            t = t + 1
            if(pp_two(1) < pp_two(2))then
              p_two(1:2,t) = pp_two(1:2)
              fswap(t) = .false.
            else
              p_two(1:2,t) = pp_two(2:1:-1)
              fswap(t) = .true.
            end if
          end do
        end do
      end do
      return
    end subroutine mk_p_two_and_fswap
    !
  end subroutine calcu_cf
  !
end module ham2vec
