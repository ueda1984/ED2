! Optional cubic Gamma/O_h, three-down-spin orbit-key lookup.
! This module is intentionally disjoint from the established basis builder.
! It returns the ORIGINAL representative index and the ORIGINAL character.
module cubic_anchored_lookup
  use iso_fortran_env, only: int64, real64
  implicit none
  private
  public :: setup_cubic_anchored_lookup, resolve_cached_hop
  integer, save :: side=0, nsites=0, table_size=0
  integer, allocatable, save :: points(:,:), characters(:)
  integer, allocatable, save :: indices(:), rep_characters(:)
  integer(int64), allocatable, save :: keys(:)
  logical, save :: initialized=.false., active=.false., verify=.false.
contains
  pure integer function site_index(x,y,z) result(s)
    integer, intent(in) :: x,y,z
    s=1+modulo(x,side)+side*(modulo(y,side)+side*modulo(z,side))
  end function site_index

  pure subroutine site_coordinates(s,x,y,z)
    integer, intent(in) :: s
    integer, intent(out) :: x,y,z
    x=mod(s-1,side)
    y=mod((s-1)/side,side)
    z=(s-1)/(side*side)
  end subroutine site_coordinates

  subroutine verify_cubic_generators()
    use input_param, only: shift_1,shift_2,shift_3,shift_4,shift_5,shift_6
    integer :: s,x,y,z
    do s=1,nsites
      call site_coordinates(s,x,y,z)
      if(shift_1(s)/=site_index(x+1,y,z))error stop 'anchored: invalid Tx'
      if(shift_2(s)/=site_index(x,y+1,z))error stop 'anchored: invalid Ty'
      if(shift_3(s)/=site_index(x,y,z+1))error stop 'anchored: invalid Tz'
      if(shift_4(s)/=site_index(-x,z,-y))error stop 'anchored: invalid I*C4x'
      if(shift_5(s)/=site_index(y,-x,z))error stop 'anchored: invalid C4z^-1'
      if(shift_6(s)/=site_index(z,x,y))error stop 'anchored: invalid C3[111]'
    end do
  end subroutine verify_cubic_generators

  subroutine make_point_actions()
    use input_param, only: shift_4,shift_5,shift_6,explist,L1,L2,L3
    integer :: a,b,c,g,s,k,t,previous
    complex(real64) :: chi
    allocate(points(nsites,48),characters(48))
    g=0
    do c=1,3
      do b=1,4
        do a=1,4
          g=g+1
          do s=1,nsites
            t=s
            do k=1,c
              t=shift_6(t)
            end do
            do k=1,b
              t=shift_5(t)
            end do
            do k=1,a
              t=shift_4(t)
            end do
            points(s,g)=t
          end do
          chi=explist(L1,L2,L3,a,b,c)
          if(abs(aimag(chi))>1.0e-10_real64)error stop 'anchored: nonreal character'
          if(abs(abs(real(chi,real64))-1.0_real64)>1.0e-10_real64) &
            error stop 'anchored: character is not +/-1'
          characters(g)=nint(real(chi,real64))
          do previous=1,g-1
            if(all(points(:,previous)==points(:,g)))error stop 'anchored: repeated point map'
          end do
        end do
      end do
    end do
    if(g/=48)error stop 'anchored: wrong point group size'
  end subroutine make_point_actions

  subroutine anchored_key_for_triple(spins,key,sign,conflict)
    use state_lists, only: insertion_sort
    integer, intent(in) :: spins(3)
    integer(int64), intent(out) :: key
    integer, intent(out) :: sign
    logical, intent(out) :: conflict
    integer :: g,j,a,point_xyz(3,3),candidate(3)
    integer :: tx,ty,tz
    integer(int64) :: candidate_key, wide
    wide=int(nsites,int64)
    key=huge(key)
    sign=0
    conflict=.false.
    do g=1,48
      do j=1,3
        call site_coordinates(points(spins(j),g), &
                              point_xyz(1,j),point_xyz(2,j),point_xyz(3,j))
      end do
      do a=1,3
        tx=point_xyz(1,a)
        ty=point_xyz(2,a)
        tz=point_xyz(3,a)
        do j=1,3
          candidate(j)=site_index(point_xyz(1,j)-tx,point_xyz(2,j)-ty, &
                                  point_xyz(3,j)-tz)-1
        end do
        call insertion_sort(candidate,3)
        candidate_key=(int(candidate(1),int64)*wide+int(candidate(2),int64))*wide+ &
                      int(candidate(3),int64)
        if(candidate_key<key)then
          key=candidate_key
          sign=characters(g)
          conflict=.false.
        else if(candidate_key==key.and.sign/=characters(g))then
          conflict=.true.
        end if
      end do
    end do
    if(sign==0)error stop 'anchored: empty key search'
  end subroutine anchored_key_for_triple

  subroutine setup_cubic_anchored_lookup(lv)
    use input_param, only: NOS,NODmin,NODmax,local_NODmax,local_spin, &
      L1,L2,L3,L4,L5,L6,M1,M2,M3,M4,M5,M6,list_s,NO_one,NO_two,Jint
    use state_lists, only: list_fly
    integer, intent(in) :: lv
    character(len=32) :: setting,check_setting
    integer :: status,status2,i,p,slot,k,sign
    integer :: triple(3)
    integer(int64) :: key
    logical :: conflict
    if(initialized)return
    initialized=.true.
    setting=''
    call get_environment_variable('ED2_CUBIC_ANCHORED_REP',setting,status=status)
    if(status/=0.or.trim(adjustl(setting))/='1')return
    check_setting=''
    call get_environment_variable('ED2_CUBIC_ANCHORED_VERIFY',check_setting,status=status2)
    verify=status2==0.and.trim(adjustl(check_setting))=='1'
    nsites=NOS
    side=nint(real(NOS,real64)**(1.0_real64/3.0_real64))
    if(side<4.or.side**3/=NOS)error stop 'anchored: not a cubic volume'
    if(NODmin/=3.or.NODmax/=3)error stop 'anchored: requires exactly 3 down spins'
    if(any(local_NODmax/=1).or.any(local_spin/=0.5_real64)) &
      error stop 'anchored: requires spin 1/2'
    if(L1/=side.or.L2/=side.or.L3/=side.or.L4/=4.or.L5/=4.or.L6/=3) &
      error stop 'anchored: incompatible symmetry periods'
    if(any([M1,M2,M3]/=0).or.M6/=0)error stop 'anchored: requires Gamma and real 1D O_h'
    if(.not.((M4==0.or.M4==2).and.(M5==0.or.M5==2))) &
      error stop 'anchored: unrecognized point group character'
    if(NO_one/=0.or.NO_two/=9*NOS)error stop 'anchored: incompatible interaction shape'
    if(any(abs(Jint(1,:)-Jint(2,:))>1.0e-12_real64).or. &
       any(abs(Jint(4:9,:))>1.0e-12_real64)) &
       error stop 'anchored: requires conserving real XXZ'
    if(lv/=size(list_s))error stop 'anchored: invalid projected basis dimension'
    call verify_cubic_generators()
    call make_point_actions()
    table_size=2*lv+1
    allocate(keys(table_size),indices(table_size),rep_characters(table_size))
    keys=-1_int64
    indices=0
    rep_characters=0
    do i=1,lv
      triple=list_fly(list_s(i),3,3,NOS)
      call anchored_key_for_triple(triple,key,sign,conflict)
      if(conflict)error stop 'anchored: invalid stabilizer character for basis representative'
      slot=1+int(modulo(key,int(table_size,int64)))
      do k=1,table_size
        if(indices(slot)==0)exit
        if(keys(slot)==key)error stop 'anchored: duplicate orbit key in projected basis'
        slot=mod(slot,table_size)+1
      end do
      if(indices(slot)/=0)error stop 'anchored: full lookup table'
      keys(slot)=key
      indices(slot)=i
      rep_characters(slot)=sign
    end do
    active=.true.
    write(*,'(a,i0,a,l1)') &
      'Representative lookup: cubic-anchored-orbits, basis=',lv, &
      ', O_h operations=48, per-hop legacy verification=',verify
  end subroutine setup_cubic_anchored_lookup

  subroutine resolve_cached_hop(spins,lv,id,phase)
    use input_param, only: list_s,explist
    use state_lists, only: representative_SQ, findstate
    integer, intent(inout) :: spins(:)
    integer, intent(in) :: lv
    integer, intent(out) :: id
    complex(real64), intent(out) :: phase
    integer :: s,ell(6),slot,k,sign,legacy_id,legacy_s,legacy_ell(6)
    integer :: old_spins(3)
    logical :: conflict
    integer(int64) :: key
    if(.not.active)then
      call representative_SQ(size(spins),s,ell,spins)
      call findstate(s,list_s,id,lv)
      phase=explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
      return
    end if
    if(size(spins)/=3)error stop 'anchored: non-three-spin cache hop'
    call anchored_key_for_triple(spins,key,sign,conflict)
    slot=1+int(modulo(key,int(table_size,int64)))
    id=-1
    phase=(1.0_real64,0.0_real64)
    do k=1,table_size
      if(indices(slot)==0)exit
      if(keys(slot)==key)then
        if(conflict)error stop 'anchored: inconsistent character in surviving orbit'
        id=indices(slot)
        phase=cmplx(real(sign*rep_characters(slot),real64),0.0_real64,kind=real64)
        exit
      end if
      slot=mod(slot,table_size)+1
    end do
    if(verify)then
      old_spins=spins
      call representative_SQ(3,legacy_s,legacy_ell,old_spins)
      call findstate(legacy_s,list_s,legacy_id,lv)
      if(legacy_id/=id)error stop 'anchored: reference representative ID disagrees'
      if(id>0)then
        if(abs(phase-explist(legacy_ell(1),legacy_ell(2),legacy_ell(3), &
                             legacy_ell(4),legacy_ell(5),legacy_ell(6)))>1.0e-10_real64) &
          error stop 'anchored: reference character disagrees'
      end if
    end if
  end subroutine resolve_cached_hop
end module cubic_anchored_lookup
