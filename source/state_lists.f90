module state_lists
  implicit none
  integer, allocatable :: cum_ltd_rep_combination(:,:)
  integer, allocatable :: ltd_rep_combination(:,:)
  ! Per-worker output avoids serial insertion and preserves deterministic sorting.
  type, private :: orbit_buffer
    integer :: count=0
    integer, allocatable :: states(:)
    real(8), allocatable :: norms(:)
  end type orbit_buffer
contains

  subroutine get_cum_ltd_rep_combination(NOS,NODmax,NODmin,m) 
    use iso_fortran_env, only: int64
    integer, intent(in) :: NOS, NODmax, NODmin, m(NOS)
    integer :: i, j, dN
    integer(int64) :: wide
    dN = NODmax-NODmin
    call get_ltd_rep_combination(NOS,NODmax,m)
    allocate(cum_ltd_rep_combination(0:NOS+1,0:NODmax+1))
    cum_ltd_rep_combination(:,:) = 0
    if(NODmin == 0)then
      j = 0
      cum_ltd_rep_combination(:,j) = 1
      do j = 1, NODmax
        do i = 0, NOS+1
          wide=int(cum_ltd_rep_combination(i,j-1),int64)+int(ltd_rep_combination(i,j),int64)
          if(wide>huge(1).or.wide<0)error stop "Hilbert count exceeds supported integer range"
          cum_ltd_rep_combination(i,j)=int(wide)
        end do
      end do
    else
      j = 0
      cum_ltd_rep_combination(:,j) = 1
      do j = 1, NODmax
        if(j <= dN)then
          do i = 0, NOS+1
            wide=int(cum_ltd_rep_combination(i,j-1),int64)+int(ltd_rep_combination(i,j),int64)
            if(wide>huge(1).or.wide<0)error stop "Hilbert count exceeds supported integer range"
            cum_ltd_rep_combination(i,j)=int(wide)
          end do
        else
          do i = 0, NOS+1
            wide=int(cum_ltd_rep_combination(i,j-1),int64)+int(ltd_rep_combination(i,j),int64) - &
              int(ltd_rep_combination(i,j-dN-1),int64)
            if(wide>huge(1).or.wide<0)error stop "Hilbert count exceeds supported integer range"
            cum_ltd_rep_combination(i,j)=int(wide)
          end do
        end if
      end do      
    end if
    !
  end subroutine get_cum_ltd_rep_combination

  subroutine get_ltd_rep_combination(NOS,NODmax,m) 
    use iso_fortran_env, only: int64
    integer, intent(in) :: NOS, NODmax, m(NOS)
    integer :: i, j
    integer(int64) :: wide
    integer, allocatable :: dp(:,:)
    allocate(ltd_rep_combination(0:NOS+1,0:NODmax+1))
    allocate(dp(0:NODmax+1,0:NOS+1))
    dp(:,:) = 0
    dp(0,:) = 1
    do i = 1, NOS
      do j = 1, NODmax
        if(j-1-m(i) >= 0)then
          wide=int(dp(j-1,i),int64)+int(dp(j,i-1),int64)-int(dp(j-1-m(i),i-1),int64)
        else
          wide=int(dp(j-1,i),int64)+int(dp(j,i-1),int64)
        end if
        if(wide>huge(1).or.wide<0)error stop "Hilbert count exceeds supported integer range"
        dp(j,i)=int(wide)
      end do
    end do
    ltd_rep_combination = transpose(dp)
  end subroutine get_ltd_rep_combination

  subroutine checkstate_SQ(s,ro,NODmax,NODmin,pk1,pk2,pk3,pk4,pk5,pk6) 
    use input_param, only: L1,L2,L3,L4,L5,L6,NOS,shift_1,shift_2,shift_3,shift_4,shift_5,shift_6
    integer, intent(in) :: s
    integer, intent(in) :: NODmax,NODmin
    real(8), intent(in) :: pk1,pk2,pk3,pk4,pk5,pk6
    real(8), intent(out) :: ro
    complex(8) :: coef
    integer :: i1,i2,i3,i4,i5,i6,l,count
    integer, allocatable :: n(:), np(:)
    real(8) :: p1, p2, p3, p4, p5
    allocate(n(NODmax),np(NODmax))
    n = list_fly(s,NODmax,NODmin,NOS)
    np = n
    !
    ro = 0.0d0
    count = 0
    coef = (0.0d0,0.0d0)
    !
    do i6 = 1, L6
      np(:) = shift_6(np(:))
      p5 = pk6 * dble(i6)
      do i5 = 1, L5
        np(:) = shift_5(np(:))
        p4 = pk5 * dble(i5) + p5
        do i4 = 1, L4
          np(:) = shift_4(np(:))
          p3 = pk4 * dble(i4) + p4
          do i3 = 1, L3
            np(:) = shift_3(np(:))
            p2 = pk3 * dble(i3) + p3
            do i2 = 1, L2
              np(:) = shift_2(np(:))
              p1 = pk2 * dble(i2) + p2
              do i1 = 1, L1
                np(:) = shift_1(np(:))
                call insertion_sort(np,NODmax)
                !
                do l = NODmax, 1, -1
                  if(np(l) < n(l))then
                    return
                  else if(np(l) > n(l))then
                    exit
                  end if
                end do
                if(l==0)then
                  coef = coef + exp( (0.0d0,-1.0d0) * ( pk1 * dble(i1) + p1) )
                  count = count + 1
                end if
                !
              end do
            end do
          end do
        end do
      end do
    end do
    !
    ro = (abs(coef)**2) *dble(L1*L2*L3*L4*L5*L6)/dble(count)
    !NOS/count='Ra', where, T^{Ra}|s>=|s>, and  ro='Na', see the eq. (118) in the Sandvik ref.
  end subroutine checkstate_SQ


  pure subroutine insertion_sort(a,NODmax)  
    integer, intent(in) :: NODmax
    integer, intent(inout) :: a(NODmax)
    integer :: temp, i, j
    do i = 2, NODmax
      j = i - 1
      temp = a(i)
      do while (a(j)>temp) 
        a(j+1) = a(j)
        j = j - 1
        if(j==0) exit
      end do
      a(j+1) = temp
    end do
  end subroutine insertion_sort

  integer function inv_list(ni,NODmax) 
    integer, intent(in) :: NODmax
    integer, intent(in) :: ni(NODmax)     
    integer :: i, j
    do i = 1, NODmax
      if(ni(i)>0) exit
    end do
    inv_list = 1
    if(i < NODmax+1) then
      do j = i, NODmax
        inv_list = inv_list + cum_ltd_rep_combination(ni(j)-1,j)
      end do
    end if
  end function inv_list

  recursive subroutine qsort_w_order(a, o, first, last) 
    integer :: a(*), o(*), first, last
    integer :: x, t8, i, j
    integer :: t
    x = a( (first+last) / 2 )
    i = first
    j = last
    do
      do while (a(i) < x)
        i=i+1
      end do
      do while (x < a(j))
        j=j-1
      end do
      if (i >= j) exit
      t8 = a(i);  a(i) = a(j);  a(j) = t8
      t  = o(i);  o(i) = o(j);  o(j) = t
      i=i+1
      j=j-1
    end do
    if (first < i - 1) call qsort_w_order(a, o, first, i - 1)
    if (j + 1 < last)  call qsort_w_order(a, o, j + 1, last)
  end subroutine qsort_w_order

  subroutine allocate_lists_omp_SQ 
    use input_param, only: NODmax,NODmin,rk1,rk2,rk3,rk4,rk5,rk6,THS,NOS,list_s,list_r,&
      L1,L2,L3,L4,L5,L6
    use omp_lib
    integer :: dim, a, i, j
    logical :: orbit_used
    real(8) :: r
    integer :: num = 1
    real(8), allocatable :: tmp_list_r(:)
    integer, allocatable :: order(:)
    real(8), allocatable :: tmp_list_r_pre(:)
    integer, allocatable :: list_s_pre(:)
    ! Optional experimental anchored cubic Gamma/O_h basis construction.
    ! Uses at most C(N-1,2) candidate translation-reduced shapes.
    call try_cubic_anchored_basis(orbit_used)
    if(orbit_used)return
    call try_orbit_basis(orbit_used)
    if(orbit_used)return
    num = omp_get_max_threads()
    write(*,*) "********************"
    write(*,*) "max_threads", num
    write(*,*) "********************"
    !
    dim = max(2*THS/(L1*L2*L3*L4*L5*L6),10000)
    allocate(list_s_pre(dim),tmp_list_r_pre(dim))
    write(*,'(" ### Count # of representative states. ")')
    a = 0
    !$omp parallel
    !$omp do private(j,i,r)
    do j = 1, num
      do i = j, THS, num
        call checkstate_SQ(i,r,NODmax,NODmin,rk1,rk2,rk3,rk4,rk5,rk6)
        if(r>1.0d-15) then
          !$omp critical
          a = a + 1
          list_s_pre(a) = i
          tmp_list_r_pre(a) = r
          !$omp end critical
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! !$omp parallel
    ! !$omp do private(j,i,r)
    ! do j = 1, num
    !   do i = j, THS, num
    !     call checkstate_SQ(i,r,NODmax,NODmin,rk1,rk2,rk3,rk4,rk5,rk6)
    !     if(r>1.0d-15) then
    !       !$omp atomic
    !       a = a + 1
    !     end if
    !   end do
    ! end do
    ! !$omp end do
    ! !$omp end parallel

    if(a==0) stop 'THS = 0.'

    !
    write(*,'(" ### Allocate work arrays for lists. ")')
    allocate(list_s(a),tmp_list_r(a))
    !
    write(*,'(" ### Store representative states. ")')
    ! a = 0
    ! !$omp parallel
    ! !$omp do private(j,i,r)
    ! do j = 1, num
    !   do i = j, THS, num
    !     call checkstate_SQ(i,r,NODmax,NODmin,rk1,rk2,rk3,rk4,rk5,rk6) 
    !     if(r>1.0d-15) then
    !       !$omp critical
    !       a = a + 1
    !       list_s(a) = i
    !       tmp_list_r(a) = r
    !       !$omp end critical
    !     end if
    !   end do
    ! end do
    ! !$omp end do
    ! !$omp end parallel

    !$omp parallel do
    do i = 1, a
       list_s(i) = list_s_pre(i)
       tmp_list_r(i) = tmp_list_r_pre(i)
    end do
    !
    deallocate(list_s_pre,tmp_list_r_pre)    
    !
    write(*,'(" ### Allocate arrays for lists. ")')
    THS = a
    allocate(order(THS),list_r(THS))
    !
    write(*,'(" ### Store reordered lists. ")')
    !$omp parallel do
    do i = 1, THS
      order(i) = i
    end do
    call qsort_w_order(list_s, order, 1, THS)
    !$omp parallel do
    do i = 1, THS
      list_r(i) = sqrt( tmp_list_r(order(i)) )
    end do
    !
    return
  end subroutine allocate_lists_omp_SQ

  function list_fly(t,NODmax,NODmin,NOS) result(ni) 
    integer, intent(in) :: t, NODmax, NODmin, NOS
    integer :: s
    integer :: i, j, b, j0
    integer :: ni(NODmax) 
    s = t
    j = NOS
    do i = NODmax, 1, -1
      call binary_search(s,cum_ltd_rep_combination(0,i), 0, j, b, j0)
      !write(*,*) "s, i, j, b, j0", s, i, j, b, j0
       j = j0
       if(j0 == 0)then
         ni(1:i) = 0
         return
       end if
       ni(i) = j0
       s = s - cum_ltd_rep_combination(j-1,i)
    end do
  end function list_fly
  !
  subroutine binary_search(s,list_s,ls,le,b,bmin)
    ! Select the FIRST cumulative count >= s; repeated counts occur with caps=0.
    integer, intent(out) :: b,bmin
    integer, intent(in) :: s,ls,le,list_s(ls:le)
    integer :: bmax
    bmin=ls;bmax=le+1
    do while(bmin<bmax)
      b=bmin+(bmax-bmin)/2
      if(list_s(b)<s)then
        bmin=b+1
      else
        bmax=b
      end if
    end do
    b=-1
    if(bmin<=le)then
      if(list_s(bmin)==s)b=bmin
    end if
  end subroutine binary_search

  function a_down_spin(i,n,NODmax) result(nd) 
    integer, intent(in) :: NODmax, i, n(NODmax)
    integer :: kr
    integer :: nd(NODmax)
    !
    do kr = NODmax, 1, -1
      if(i == n(kr)) exit
    end do
    nd = (/ 0, n(1:kr-1), n(kr+1:) /)
  end function a_down_spin

  function a2_down_spin(i,j,n,NODmax) result(nd) 
    integer, intent(in) :: NODmax, i, j, n(NODmax)
    integer :: kl,kr
    integer :: nd(NODmax)
    ! we assume i < j
    !
    do kr = NODmax, 1, -1
      if(j == n(kr)) exit
    end do
    do kl = kr-1, 1, -1
      if(i == n(kl)) exit
    end do
    nd = (/ 0, 0, n(1:kl-1), n(kl+1:kr-1), n(kr+1:) /)
  end function a2_down_spin

  function c_down_spin(i,n,NODmax) result(nd) 
    integer, intent(in) :: NODmax, i, n(NODmax)
    integer :: kr
    integer :: nd(NODmax)
    !    
    do kr = NODmax, 1, -1
      if(n(kr) <= i) exit
    end do
    nd = (/ n(2:kr), i, n(kr+1:) /)
  end function c_down_spin

  function c2_down_spin(i,j,n,NODmax) result(nd) 
    integer, intent(in) :: NODmax, i, j, n(NODmax)
    integer :: kl, kr
    integer :: nd(NODmax)
    ! we assume i < j
    !
    do kr = NODmax, 1, -1
      if(n(kr) <= j) exit
    end do
    do kl = kr, 1, -1
      if(n(kl) <= i) exit
    end do
    nd = (/ n(3:kl), i, n(kl+1:kr), j, n(kr+1:) /)
  end function c2_down_spin

  function j_flip_ni(i,j,n,NODmax) result(nd) 
    integer, intent(in) :: NODmax, i, j, n(NODmax)
    integer :: kr, kl
    integer :: nd(NODmax)
    nd = 1
    !
    ! i < j
    ! (....n_kl-1,n_kl,n_kl+1......n_kr-1,n_kr,n_kr+1,....)
    ! (....n_kl-1,n_kl+1......n_kr-1,n_kr,n_kl,n_kr+1,....)
    ! j < i
    ! (....n_kl-1,n_kl,n_kl+1......n_kr-1,n_kr,n_kr+1,....)
    ! (....n_kl-1,n_kr,n_kl,n_kl+1......n_kr-1,n_kr+1,....)
    !
    if(i < j)then
      do kl = 1, NODmax
        if(i == n(kl)) exit
      end do
      do kr = NODmax, 1, -1
        if(j >= n(kr)) exit
      end do
      !write(*,*) "kl,kr", kl, kr
      nd(kl:kr-1) = n(kl+1:kr)
      nd(kr) = j
    else if(j < i)then
      do kr = NODmax, 1, -1
        if(i == n(kr)) exit
      end do
      do kl = 1, NODmax
        if(j <= n(kl)) exit
      end do
      nd(kl) = j
      nd(kl+1:kr) = n(kl:kr-1)
    end if
    nd(:kl-1) = n(:kl-1)
    nd(kr+1:) = n(kr+1:)
  end function j_flip_ni

  subroutine findstate(s,list_s,b,m) 
    integer, intent(out) :: b
    integer, intent(in) :: s, m, list_s(m)
    integer :: bmin, bmax 
    if(m>0)then
      if(list_s(1)==1.and.list_s(m)==m)then
        b=-1
        if(s>=1.and.s<=m)b=s
        return
      end if
    end if
    bmin = 1; bmax = m
    do
      b = bmin + (bmax-bmin)/2
      if(s < list_s(b))then
        bmax = b-1
      else if(s > list_s(b))then
        bmin = b+1
      else
        return
      end if
      if(bmin>bmax)then
        b = -1; return
      end if
    end do
  end subroutine findstate

  subroutine representative_SQ(NODmax,t,ell,n) 
    use input_param, only: L1,L2,L3,L4,L5,L6,shift_1,shift_2,shift_3,shift_4,shift_5,shift_6
    integer,intent(in)::NODmax 
    integer,intent(inout)::n(NODmax) 
    integer,intent(out) :: t, ell(6)
    integer :: i1,i2,i3,i4,i5,i6,l
    integer, allocatable :: nd(:)
    if(all([L1,L2,L3,L4,L5,L6]==1))then
      t=inv_list(n,NODmax)
      ell=1
      return
    end if
    allocate(nd(NODmax))
    ell = (/L1,L2,L3,L4,L5,L6/)
    nd = n
    do i6 = 1, L6
      nd(:) = shift_6(nd(:))
      do i5 = 1, L5
        nd(:) = shift_5(nd(:))
        do i4 = 1, L4
          nd(:) = shift_4(nd(:))
          do i3 = 1, L3
            nd(:) = shift_3(nd(:))
            do i2 = 1, L2
              nd(:) = shift_2(nd(:))
              do i1 = 1, L1
                nd(:) = shift_1(nd(:))
                call insertion_sort(nd,NODmax)
                !
                do l = NODmax, 1, -1
                  if(nd(l)>n(l))then
                    exit
                  else if(nd(l)<n(l))then
                    ell(1)=i1
                    ell(2)=i2
                    ell(3)=i3
                    ell(4)=i4
                    ell(5)=i5
                    ell(6)=i6
                    n(1:l) = nd(1:l)
                    exit
                  end if
                end do
                !
              end do
            end do
          end do
        end do
      end do
    end do
    t = inv_list(n,NODmax)
  end subroutine representative_SQ
  !
  ! A product of cyclic generators need not be Abelian.  The existing orbit
  ! algorithm can also process noncommuting generators, but only when the
  ! supplied normal form enumerates each finite permutation exactly once and
  ! its prescribed phases are a genuine one-dimensional group character.
  ! Verify both conditions before allowing the fast path; otherwise the
  ! original full enumeration remains the backward-compatible fallback.
  subroutine validate_noncommuting_action(maps,periods,valid)
    use iso_fortran_env, only: int64
    use input_param, only: NOS,explist
    integer, intent(in) :: maps(0:NOS,6),periods(6)
    logical, intent(out) :: valid
    integer, allocatable :: elements(:,:),slots(:),composition(:)
    integer(int64), allocatable :: hashes1(:),hashes2(:)
    complex(8), allocatable :: chi(:)
    integer(int64) :: nwide,h1,h2,element_limit,site_image_limit
    integer :: ng,table_size,g,idx(6),q,b,site,loc,a,probe,j,large_status
    logical :: found
    character(len=16) :: permit_large
    complex(8) :: generator_chi
    valid=.false.
    nwide=1_int64
    do a=1,6
      if(periods(a)<1)return
      if(nwide>huge(nwide)/int(periods(a),int64))return
      nwide=nwide*int(periods(a),int64)
    end do
    ! The historical safe caps remain the default. A reviewed, explicitly
    ! opted-in large-group benchmark may use up to 64 million site images
    ! (~244 MiB for int32), with a hard bound on the total group order.
    ! Normal-form uniqueness, closure and phase checks are NEVER bypassed.
    element_limit=8192_int64
    site_image_limit=4000000_int64
    permit_large=''
    call get_environment_variable('ED2_ALLOW_LARGE_GROUP_VALIDATION', &
                                  permit_large,status=large_status)
    if(large_status==0 .and. trim(adjustl(permit_large))=='1')then
      element_limit=65536_int64
      site_image_limit=64000000_int64
    end if
    if(nwide>element_limit)return
    if(nwide>site_image_limit/int(NOS+1,int64))return
    ng=int(nwide)
    table_size=2*ng+1
    allocate(elements(0:NOS,ng),composition(0:NOS),slots(table_size))
    allocate(hashes1(ng),hashes2(ng),chi(ng))
    slots=0
    do g=1,ng
      ! L1 varies fastest, identical to checkstate_SQ and representative_SQ.
      q=g-1
      do a=1,6
        idx(a)=mod(q,periods(a))+1
        q=q/periods(a)
      end do
      do site=0,NOS
        loc=site
        do a=6,1,-1
          do b=1,idx(a)
            loc=maps(loc,a)
          end do
        end do
        elements(site,g)=loc
      end do
      chi(g)=explist(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
      call hash_permutation(elements(:,g),hashes1(g),hashes2(g))
      probe=1+int(modulo(hashes1(g),int(table_size,int64)))
      do
        j=slots(probe)
        if(j==0)then
          slots(probe)=g
          exit
        end if
        if(hashes1(g)==hashes1(j) .and. hashes2(g)==hashes2(j))then
          if(all(elements(:,g)==elements(:,j)))return
        end if
        probe=mod(probe,table_size)+1
      end do
    end do
    ! Last mixed-radix element has every exponent at its period: identity.
    do site=0,NOS
      if(elements(site,ng)/=site)return
    end do
    if(abs(chi(ng)-(1d0,0d0))>1d-10)return
    ! Left multiplication by each generator must remain in the enumerated
    ! set, and phase multiplication must respect the same group law. Because
    ! the set contains the identity and is finite, closure under generators
    ! certifies the finite group. Hash collisions are checked by exact maps.
    do a=1,6
      if(periods(a)==1)cycle
      idx=periods
      idx(a)=1
      generator_chi=explist(idx(1),idx(2),idx(3),idx(4),idx(5),idx(6))
      do g=1,ng
        composition=maps(elements(:,g),a)
        call hash_permutation(composition,h1,h2)
        probe=1+int(modulo(h1,int(table_size,int64)))
        found=.false.
        do
          j=slots(probe)
          if(j==0)exit
          if(h1==hashes1(j) .and. h2==hashes2(j))then
            if(all(composition==elements(:,j)))then
              found=.true.
              if(abs(chi(j)-generator_chi*chi(g))>1d-10)return
              exit
            end if
          end if
          probe=mod(probe,table_size)+1
        end do
        if(.not.found)return
      end do
    end do
    valid=.true.
  contains
    subroutine hash_permutation(perm,hash1,hash2)
      integer, intent(in) :: perm(0:)
      integer(int64), intent(out) :: hash1,hash2
      integer :: site_index
      hash1=17_int64
      hash2=97_int64
      do site_index=0,ubound(perm,1)
        hash1=modulo(hash1*65599_int64+int(perm(site_index),int64),1000000007_int64)
        hash2=modulo(hash2*65521_int64+int(perm(site_index),int64),1000000009_int64)
      end do
    end subroutine hash_permutation
  end subroutine validate_noncommuting_action

  ! Exact, small certificate for the *canonical* spin-1/2 cubic Gamma/O_h
  ! normal form, enabled only by ED2_CUBIC_FAST_GROUP_VALIDATION=1.
  !
  ! The six site permutations are checked on EVERY lattice site against
  ! (Tx, Ty, Tz, A, B, C), where
  !   A: (x,y,z) -> (-x,z,-y)   (improper order four)
  !   B: (x,y,z) -> (y,-x,z)   (proper order four)
  !   C: (x,y,z) -> (z,x,y)    (proper order three).
  !
  ! Since A/B/C fix the origin and are signed coordinate permutations,
  ! 48 distinct products A^a B^b C^c constitute all of O_h (48 elements).
  ! The translations act freely on the origin and comprise Z_L^3, giving
  ! a UNIQUE full affine normal form T_x^i T_y^j T_z^k A^a B^b C^c.
  ! At Gamma every translation has character +1. Explicitly checking
  ! chi(PQ)=chi(P)chi(Q) for all 48^2 point-element pairs certifies that
  ! the ED2 phases are an honest one-dimensional full-group character.
  !
  ! Unlike the generic validator, no 48*L^3 full-site maps are materialized.
  ! General noncommuting groups still use validate_noncommuting_action.
  pure function cubic_compose_signed(left,right) result(out)
    integer, intent(in) :: left(3),right(3)
    integer :: out(3),axis,i
    do i=1,3
      axis=abs(right(i))
      if(right(i)>0)then
        out(i)=left(axis)
      else
        out(i)=-left(axis)
      end if
    end do
  end function cubic_compose_signed

  subroutine validate_cubic_oh_gamma_certificate(maps,periods,valid)
    use iso_fortran_env, only: int64,real64
    use input_param, only: NOS,NODmin,NODmax,M1,M2,M3,M4,M5,M6, &
      local_NODmax,local_spin,explist
    integer, intent(in) :: maps(0:NOS,6),periods(6)
    logical, intent(out) :: valid
    integer, parameter :: ident(3)=[1,2,3]
    integer, parameter :: ag(3)=[-1,-3,2]
    integer, parameter :: bg(3)=[-2,1,3]
    integer, parameter :: cg(3)=[2,3,1]
    integer :: aside(3,0:4),bside(3,0:4),cside(3,0:3)
    integer :: point(3,48),product(3)
    complex(real64) :: character(48)
    integer :: side,site,coords(3),target(3),a,b,c,g,h,j,k,x,y,z
    valid=.false.
    side=nint(real(NOS,real64)**(1.0_real64/3.0_real64))
    if(side<4.or.side>10.or.mod(side,2)/=0)return
    if(int(side,int64)**3/=int(NOS,int64))return
    if(any(periods/=[side,side,side,4,4,3]))return
    if(NODmin/=3.or.NODmax/=3)return
    if(any(local_NODmax/=1))return
    if(any(abs(local_spin-0.5_real64)>1.0e-12_real64))return
    if(any([M1,M2,M3]/=0).or.M6/=0)return
    if(.not.((M4==0.or.M4==2).and.(M5==0.or.M5==2)))return
    if(any(maps(0,:)/=0))return

    ! Check the supplied *runtime* maps exactly, not a substitute model.
    ! All modulo operations are mathematically defined for negative signs.
    do site=1,NOS
      x=mod(site-1,side)
      y=mod((site-1)/side,side)
      z=(site-1)/(side*side)
      coords=[x,y,z]
      do a=1,6
        select case(a)
        case(1:3)
          target=coords
          target(a)=target(a)+1
        case(4)
          target=[-x,z,-y]
        case(5)
          target=[y,-x,z]
        case(6)
          target=[z,x,y]
        end select
        j=1+modulo(target(1),side)+side*(modulo(target(2),side)+ &
          side*modulo(target(3),side))
        if(maps(site,a)/=j)return
      end do
    end do

    aside(:,0)=ident
    bside(:,0)=ident
    cside(:,0)=ident
    do a=1,4
      aside(:,a)=cubic_compose_signed(ag,aside(:,a-1))
      bside(:,a)=cubic_compose_signed(bg,bside(:,a-1))
      if(a<=3)cside(:,a)=cubic_compose_signed(cg,cside(:,a-1))
    end do
    if(any(aside(:,4)/=ident).or.any(bside(:,4)/=ident))return
    if(any(cside(:,3)/=ident))return

    g=0
    do c=1,3
      do b=1,4
        do a=1,4
          g=g+1
          point(:,g)=cubic_compose_signed( &
            cubic_compose_signed(aside(:,a),bside(:,b)),cside(:,c))
          ! Gamma translation identity: indices 1:3 equal their periods.
          character(g)=explist(side,side,side,a,b,c)
          if(abs(aimag(character(g)))>1.0e-10_real64)return
          if(abs(abs(real(character(g),real64))-1.0_real64)>1.0e-10_real64)return
          do h=1,g-1
            if(all(point(:,h)==point(:,g)))return
          end do
        end do
      end do
    end do
    if(g/=48)return
    if(any(point(:,48)/=ident))return
    if(abs(character(48)-(1.0_real64,0.0_real64))>1.0e-10_real64)return

    ! A composition must land in the 48-element normal-form table with
    ! the product character. Repeated or missing elements fail closed.
    do g=1,48
      do h=1,48
        product=cubic_compose_signed(point(:,g),point(:,h))
        k=0
        do j=1,48
          if(all(product==point(:,j)))then
            k=j
            exit
          end if
        end do
        if(k==0)return
        if(abs(character(k)-character(g)*character(h))>1.0e-10_real64)return
      end do
    end do
    valid=.true.
  end subroutine validate_cubic_oh_gamma_certificate

  ! Opt-in experimental cubic 3-magnon basis builder.  The exact orbit key is
  ! one spin at site index zero plus an unordered pair of nonzero sites.
  ! Every affine orbit is enumerated once by marking its 48*3 anchored views.
  ! Its historical ED2 canonical representative is the minimum colex rank
  ! among 48 point transforms and <=3^3 translations.  This restriction is
  ! EXACT: at a lexically minimal translate, the minimum occupied x/y/z
  ! coordinate is zero in each axis separately (otherwise shifting the
  ! respective axis by -1 decreases all three rank indices).
  !
  ! The original stabilizer character sum is reconstructed in the exact
  ! historical normal-form order.  Gamma translations have phase one, so
  ! every stabilizing translation of a fixed point operation contributes
  ! exactly the same floating-point complex summand.
  subroutine try_cubic_anchored_basis(used)
    use iso_fortran_env, only: int64, real64
    use input_param, only: NOS,NODmax,NODmin,L1,L2,L3,L4,L5,L6, &
      rk1,rk2,rk3,rk4,rk5,rk6,shift_1,shift_2,shift_3, &
      shift_4,shift_5,shift_6, &
      list_s,list_r,THS
    use omp_lib, only: omp_get_wtime
    implicit none
    logical, intent(out) :: used
    integer :: maps(0:NOS,6),periods(6)
    integer, allocatable :: point_map(:,:), order(:)
    integer, allocatable :: sites_xyz(:,:), candidate_ranks(:)
    real(8), allocatable :: candidate_norms(:)
    logical, allocatable :: visited(:)
    logical :: certificate_valid,detail_timing,coord_cache_enabled
    character(len=16) :: enable,detail_setting,basis_dump_setting,legacy_setting,cert_setting, &
      coord_cache_setting
    integer :: status,detail_status,basis_dump_status,legacy_status,cert_status,coord_cache_status
    integer :: side,pairs,pair_a,pair_b,site,g,i,j,k,a,b,c,idx, &
      idx_a,idx_b,pi,pj,pk,best(3),shifted(3),point_sites(3), &
      point_xyz(3,3),tmp(3),source_xyz(3,3),target_xyz(3), &
      ntrans,stab_count,nrep,old_rank,basis_unit,basis_ios,delta_xyz(3,3,3)
    integer :: t
    real(8) :: p1,p2,p3,p4,p5,ro
    real(8) :: begin_s,validation_start,validation_end,scan_start,scan_end,finish_s
    complex(8) :: coef,phase
    used=.false.
    enable=''
    call get_environment_variable('ED2_CUBIC_ANCHORED_BASIS',enable,status=status)
    coord_cache_setting=''
    call get_environment_variable('ED2_CUBIC_BASIS_COORD_CACHE', &
      coord_cache_setting,status=coord_cache_status)
    if(coord_cache_status==0)then
      if(trim(adjustl(coord_cache_setting))/='0'.and. &
         trim(adjustl(coord_cache_setting))/='1') &
        error stop 'Basis coordinate cache flag must be 0 or 1'
    end if
    coord_cache_enabled=coord_cache_status==0.and. &
      trim(adjustl(coord_cache_setting))=='1'
    if(status/=0.or.trim(adjustl(enable))/='1')then
      if(coord_cache_enabled) &
        error stop 'Basis coordinate cache requires ED2_CUBIC_ANCHORED_BASIS=1'
      return
    end if
    legacy_setting=''
    call get_environment_variable('ED2_FORCE_LEGACY_ORBITS',legacy_setting,status=legacy_status)
    if(legacy_status==0.and.trim(adjustl(legacy_setting))=='1')then
      if(coord_cache_enabled) &
        error stop 'Basis coordinate cache is incompatible with legacy orbit mode'
      return
    end if
    cert_setting=''
    call get_environment_variable('ED2_CUBIC_FAST_GROUP_VALIDATION', &
      cert_setting,status=cert_status)
    if(cert_status/=0.or.trim(adjustl(cert_setting))/='1') &
      error stop 'Anchored cubic basis requires ED2_CUBIC_FAST_GROUP_VALIDATION=1'
    if(NODmin/=3.or.NODmax/=3) &
      error stop 'Anchored cubic basis requires exactly three down spins'
    detail_setting=''
    call get_environment_variable('ED2_BASIS_DETAIL_TIMING',detail_setting,status=detail_status)
    detail_timing=detail_status==0.and.trim(adjustl(detail_setting))=='1'
    if(detail_timing)begin_s=omp_get_wtime()
    periods=[L1,L2,L3,L4,L5,L6]
    ! A runtime copy of all six operators is used for exact certification.
    maps(:,1)=shift_1;maps(:,2)=shift_2;maps(:,3)=shift_3
    maps(:,4)=shift_4;maps(:,5)=shift_5;maps(:,6)=shift_6
    if(detail_timing)validation_start=omp_get_wtime()
    call validate_cubic_oh_gamma_certificate(maps,periods,certificate_valid)
    if(.not.certificate_valid)error stop 'Anchored cubic basis requires exact cubic Gamma/O_h certificate'
    write(*,'(a)') 'Symmetry action: noncommuting normal form and character validated'
    write(*,'(a)') 'Symmetry certificate: canonical cubic O_h semidirect product verified'
    if(detail_timing)validation_end=omp_get_wtime()
    side=L1
    pairs=(NOS-1)*(NOS-2)/2
    allocate(visited(pairs),point_map(NOS,48), &
      sites_xyz(3,NOS),candidate_ranks(pairs),candidate_norms(pairs))
    visited=.false.
    do i=1,NOS
      sites_xyz(1,i)=mod(i-1,side)
      sites_xyz(2,i)=mod((i-1)/side,side)
      sites_xyz(3,i)=(i-1)/(side*side)
    end do
    ! The exact ED2 point normal form uses a=1..4 fastest, then b,c.
    g=0
    do c=1,3
      do b=1,4
        do a=1,4
          g=g+1
          do site=1,NOS
            t=site
            do i=1,c
              t=shift_6(t)
            end do
            do i=1,b
              t=shift_5(t)
            end do
            do i=1,a
              t=shift_4(t)
            end do
            point_map(site,g)=t
          end do
        end do
      end do
    end do
    if(g/=48)error stop 'Anchored cubic basis point-group cardinality'
    if(detail_timing)scan_start=omp_get_wtime()
    nrep=0
    do pair_b=2,NOS-1
      do pair_a=1,pair_b-1
        idx=(pair_b-1)*(pair_b-2)/2+pair_a
        if(visited(idx))cycle
        best=[NOS+1,NOS+1,NOS+1]
        do g=1,48
          point_sites=[point_map(1,g),point_map(pair_a+1,g),point_map(pair_b+1,g)]
          point_xyz=sites_xyz(:,point_sites)
          ! Compute each periodic coordinate difference just once per point
          ! image, rather than repeating division/modulo in 3 anchor views
          ! and 27 canonical translations. Integer-only; rank and norm
          ! accumulation order remain exactly those of the legacy builder.
          if(coord_cache_enabled)then
            do i=1,3
              do j=1,3
                do k=1,3
                  delta_xyz(k,j,i)=modulo(point_xyz(k,j)-point_xyz(k,i),side)
                end do
              end do
            end do
          end if
          ! Mark the three translation anchors of each O_h image.
          do i=1,3
            if(coord_cache_enabled)then
              do j=1,3
                shifted(j)=delta_xyz(1,j,i)+ &
                  side*(delta_xyz(2,j,i)+side*delta_xyz(3,j,i))
              end do
            else
              do j=1,3
                shifted(j)=modulo(point_xyz(1,j)-point_xyz(1,i),side)+ &
                  side*(modulo(point_xyz(2,j)-point_xyz(2,i),side)+ &
                  side*modulo(point_xyz(3,j)-point_xyz(3,i),side))
              end do
            end if
            call insertion_sort(shifted,3)
            if(shifted(1)/=0)error stop 'Anchored cubic shape lacks origin'
            idx_a=shifted(2);idx_b=shifted(3)
            if(idx_a<1.or.idx_b>=NOS.or.idx_a>=idx_b) &
              error stop 'Anchored cubic shape malformed'
            idx=(idx_b-1)*(idx_b-2)/2+idx_a
            visited(idx)=.true.
          end do
          ! Exact global-rank minimization: one spin must be at coordinate
          ! zero along each axis at any optimal translation, giving <=27.
          do pi=1,3
            do pj=1,3
              do pk=1,3
                if(coord_cache_enabled)then
                  do j=1,3
                    tmp(j)=1+delta_xyz(1,j,pi)+ &
                      side*(delta_xyz(2,j,pj)+side*delta_xyz(3,j,pk))
                  end do
                else
                  do j=1,3
                    tmp(j)=1+modulo(point_xyz(1,j)-point_xyz(1,pi),side)+ &
                      side*(modulo(point_xyz(2,j)-point_xyz(2,pj),side)+ &
                      side*modulo(point_xyz(3,j)-point_xyz(3,pk),side))
                  end do
                end if
                call insertion_sort(tmp,3)
                do k=3,1,-1
                  if(tmp(k)<best(k))then
                    best=tmp
                    exit
                  else if(tmp(k)>best(k))then
                    exit
                  end if
                end do
              end do
            end do
          end do
        end do
        if(any(best>NOS))error stop 'Anchored cubic canonical state missing'
        old_rank=inv_list(best,3)
        ! Reproduce the original coef accumulation order i6/i5/i4/i3/i2/i1.
        do j=1,3
          source_xyz(:,j)=sites_xyz(:,best(j))
        end do
        coef=(0d0,0d0)
        stab_count=0
        do g=1,48
          a=mod(g-1,4)+1
          b=mod((g-1)/4,4)+1
          c=(g-1)/16+1
          point_sites=point_map(best,g)
          point_xyz=sites_xyz(:,point_sites)
          ntrans=0
          do i=1,3
            target_xyz=source_xyz(:,i)-point_xyz(:,1)
            do j=1,3
              shifted(j)=1+modulo(point_xyz(1,j)+target_xyz(1),side)+ &
                side*(modulo(point_xyz(2,j)+target_xyz(2),side)+ &
                side*modulo(point_xyz(3,j)+target_xyz(3),side))
            end do
            call insertion_sort(shifted,3)
            if(all(shifted==best))ntrans=ntrans+1
          end do
          if(ntrans==0)cycle
          ! The original p1 = rk2*i2 + (rk3*i3 + (rk4*i4 +
          ! (rk5*i5 + rk6*i6))) for each stabilizing translation.
          ! At Gamma rk1=rk2=rk3=0; using identity shift indices is exact.
          p5=rk6*dble(c)
          p4=rk5*dble(b)+p5
          p3=rk4*dble(a)+p4
          p2=rk3*dble(side)+p3
          p1=rk2*dble(side)+p2
          phase=exp((0d0,-1d0)*(rk1*dble(side)+p1))
          do i=1,ntrans
            coef=coef+phase
            stab_count=stab_count+1
          end do
        end do
        if(stab_count<1)error stop 'Anchored cubic stabilizer missing identity'
        ro=(abs(coef)**2)*dble(48*NOS)/dble(stab_count)
        if(ro>1d-15)then
          nrep=nrep+1
          candidate_ranks(nrep)=old_rank
          candidate_norms(nrep)=sqrt(ro)
        end if
      end do
    end do
    if(nrep==0)error stop 'Anchored cubic basis is empty'
    if(.not.all(visited))error stop 'Anchored cubic pair orbit coverage incomplete'
    if(detail_timing)scan_end=omp_get_wtime()
    allocate(list_s(nrep),list_r(nrep),order(nrep))
    list_s=candidate_ranks(:nrep)
    deallocate(visited,point_map,sites_xyz,candidate_ranks)
    do i=1,nrep
      order(i)=i
    end do
    call qsort_w_order(list_s,order,1,nrep)
    do i=1,nrep
      list_r(i)=candidate_norms(order(i))
    end do
    deallocate(order,candidate_norms)
    if(nrep>1)then
      if(any(list_s(2:)<=list_s(:nrep-1))) &
        error stop 'Anchored cubic basis nonunique historical ranks'
    end if
    THS=nrep
    basis_dump_setting=''
    call get_environment_variable('ED2_CUBIC_BASIS_PROOF_DUMP', &
      basis_dump_setting,status=basis_dump_status)
    if(basis_dump_status==0.and.trim(adjustl(basis_dump_setting))=='1')then
      open(newunit=basis_unit,file='cubic_basis_proof_dump.txt', &
        status='replace',action='write',iostat=basis_ios)
      if(basis_ios/=0)error stop 'Failed to open anchored cubic basis proof dump'
      do i=1,nrep
        write(basis_unit,'(i0,1x,z16.16)',iostat=basis_ios) &
          list_s(i),transfer(list_r(i),0_int64)
        if(basis_ios/=0)error stop 'Failed to write anchored cubic basis proof dump'
      end do
      close(basis_unit,iostat=basis_ios)
      if(basis_ios/=0)error stop 'Failed to close anchored cubic basis proof dump'
      write(*,'(a,i0)') 'Basis proof dump: exactly ranked states and bitwise norms=',nrep
    end if
    if(detail_timing)then
      finish_s=omp_get_wtime()
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE preflight_other wall_s=', &
        validation_start-begin_s+scan_start-validation_end
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE group_validation wall_s=', &
        validation_end-validation_start
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE orbit_scan wall_s=',scan_end-scan_start
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE finalize wall_s=',finish_s-scan_end
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE total wall_s=',finish_s-begin_s
    end if
    write(*,'(a,i0,a,i0)') 'Basis builder: anchored-cubic, orbit keys=',nrep, &
      ', translation-reduced candidates=',pairs
    if(coord_cache_enabled) &
      write(*,'(a)') 'Basis orbit coordinates: exact precomputed periodic deltas'
    used=.true.
  end subroutine try_cubic_anchored_basis

  ! Only a globally minimal representative may mark its entire orbit visited.
  ! Atomic bitset accesses are required even when different bits share a word.
  ! Canonicality is decided BEFORE publication, so a later-ranked worker can
  ! never suppress an unprocessed minimal representative. Output is sorted
  ! by the legacy rank after thread-local accumulation.
  subroutine try_orbit_basis(used)
    use iso_fortran_env, only: int64
    use input_param, only: NOS,NODmin,NODmax,THS,local_NODmax,local_spin, &
      L1,L2,L3,L4,L5,L6,shift_1,shift_2,shift_3,shift_4,shift_5,shift_6,list_s,list_r
    use omp_lib, only: omp_get_max_threads,omp_get_thread_num,omp_get_num_threads,omp_get_wtime
    logical, intent(out) :: used
    integer :: periods(6),i,j,a,b,cur,cycle_len,dim,nrep,nt,tid,actual_threads,capacity
    integer :: s,word,bit,g,offset
    integer(int64) :: group_order,nwords,seen_word,mask,total_count
    integer, allocatable :: maps(:,:),ns(:),np(:),ranks(:),order(:)
    integer(int64), allocatable :: visited(:)
    logical, allocatable :: seen(:)
    logical :: is_representative,commuting,valid_noncommuting
    integer :: env_status
    character(len=16) :: force_legacy
    real(8) :: ro
    ! Optional wall-clock breakdown. This is diagnostic only; physics is untouched.
    real(8) :: detail_begin,validation_begin,validation_sec,orbit_begin,orbit_finish,detail_finish
    logical :: detail_timing
    character(len=16) :: detail_setting,group_cert_setting,basis_dump_setting
    integer :: detail_status,group_cert_status,basis_dump_status,basis_unit,basis_i,basis_ios
    logical :: fast_cubic_certificate
    real(8), allocatable :: norms(:)
    type(orbit_buffer), allocatable :: buffers(:)
    used=.false.
    detail_setting=''
    call get_environment_variable('ED2_BASIS_DETAIL_TIMING',detail_setting,status=detail_status)
    detail_timing=detail_status==0.and.trim(adjustl(detail_setting))=='1'
    if(detail_timing)detail_begin=omp_get_wtime()
    group_cert_setting=''
    call get_environment_variable('ED2_CUBIC_FAST_GROUP_VALIDATION', &
      group_cert_setting,status=group_cert_status)
    fast_cubic_certificate=group_cert_status==0.and.trim(adjustl(group_cert_setting))=='1'
    validation_sec=0d0
    force_legacy=''
    call get_environment_variable('ED2_FORCE_LEGACY_ORBITS',force_legacy,status=env_status)
    if(env_status==0 .and. trim(adjustl(force_legacy))=='1')return
    periods=[L1,L2,L3,L4,L5,L6]
    if(any(periods<1).or.THS<1) return
    group_order=1_int64
    do i=1,6
      if(group_order>huge(1_int64)/periods(i))return
      group_order=group_order*periods(i)
    end do
    if(group_order==1_int64)then
      allocate(list_s(THS),list_r(THS))
      !$omp parallel do private(i)
      do i=1,THS
        list_s(i)=i
        list_r(i)=1d0
      end do
      write(*,'(a)') 'Basis builder: identity-direct'
      used=.true.
      return
    end if
    nt=omp_get_max_threads()
    if(group_order>int(huge(1),int64))return
    nwords=(int(THS,int64)+63_int64)/64_int64
    if(nwords>16777216_int64)return ! <=128 MiB shared bitmap, independent of nt
    if(group_order>16777216_int64/int(nt,int64))return ! <=64 MiB total rank workspaces
    allocate(maps(0:NOS,6),seen(NOS))
    maps(:,1)=shift_1;maps(:,2)=shift_2;maps(:,3)=shift_3
    maps(:,4)=shift_4;maps(:,5)=shift_5;maps(:,6)=shift_6
    if(any(maps(0,:)/=0))return
    commuting=.true.
    do a=1,6
      if(any(maps(1:,a)<1).or.any(maps(1:,a)>NOS))return
      seen=.false.
      do i=1,NOS
        j=maps(i,a)
        if(seen(j))return
        seen(j)=.true.
        if(local_NODmax(i)/=local_NODmax(j))return
        if(local_spin(i)/=local_spin(j))return
      end do
      seen=.false.
      do i=1,NOS
        if(seen(i))cycle
        cur=i;cycle_len=0
        do
          seen(cur)=.true.;cycle_len=cycle_len+1
          cur=maps(cur,a)
          if(cur==i)exit
        end do
        if(mod(periods(a),cycle_len)/=0)return
      end do
      do b=1,a-1
        do i=1,NOS
          if(maps(maps(i,a),b)/=maps(maps(i,b),a))commuting=.false.
        end do
      end do
    end do
    if(fast_cubic_certificate.and.commuting) &
      error stop 'Fast cubic group certificate requires noncommuting full O_h'
    if(.not.commuting)then
      if(detail_timing)validation_begin=omp_get_wtime()
      if(fast_cubic_certificate)then
        call validate_cubic_oh_gamma_certificate(maps,periods,valid_noncommuting)
        if(.not.valid_noncommuting) &
          error stop 'Fast cubic group certificate rejected noncanonical geometry or character'
      else
        call validate_noncommuting_action(maps,periods,valid_noncommuting)
      end if
      if(detail_timing)validation_sec=omp_get_wtime()-validation_begin
      if(.not.valid_noncommuting)return
      write(*,'(a)') 'Symmetry action: noncommuting normal form and character validated'
      if(fast_cubic_certificate) &
        write(*,'(a)') 'Symmetry certificate: canonical cubic O_h semidirect product verified'
    end if
    deallocate(seen,maps)
    dim=THS
    capacity=int(min(int(dim,int64),max(128_int64,2_int64*dim/group_order/int(nt,int64))))
    if(detail_timing)orbit_begin=omp_get_wtime()
    allocate(visited(nwords),buffers(nt))
    !$omp parallel do private(i)
    do i=1,int(nwords)
      visited(i)=0_int64
    end do
    actual_threads=1
    !$omp parallel default(shared) private(tid,ns,np,ranks,s,word,bit,seen_word,ro,is_representative,g,mask)
    tid=omp_get_thread_num()+1
    !$omp single
    actual_threads=omp_get_num_threads()
    !$omp end single
    allocate(ns(NODmax),np(NODmax),ranks(int(group_order)))
    allocate(buffers(tid)%states(capacity),buffers(tid)%norms(capacity))
    !$omp do schedule(dynamic,256)
    do s=1,dim
      word=(s-1)/64+1;bit=mod(s-1,64)
      !$omp atomic read
      seen_word=visited(word)
      if(btest(seen_word,bit))cycle
      call collect_canonical_orbit(s,ns,np,ranks,ro,is_representative)
      if(.not.is_representative)cycle
      do g=1,int(group_order)
        word=(ranks(g)-1)/64+1;bit=mod(ranks(g)-1,64)
        mask=ibset(0_int64,bit)
        !$omp atomic update
        visited(word)=ior(visited(word),mask)
      end do
      if(ro>1d-15)call append_orbit(buffers(tid),s,sqrt(ro),dim)
    end do
    !$omp end do
    deallocate(ns,np,ranks)
    !$omp end parallel
    if(detail_timing)orbit_finish=omp_get_wtime()
    total_count=0_int64
    do i=1,nt
      total_count=total_count+int(buffers(i)%count,int64)
    end do
    if(total_count>dim)error stop 'duplicate orbit representatives'
    nrep=int(total_count)
    if(nrep==0)stop 'THS = 0.'
    allocate(list_s(nrep),list_r(nrep),norms(nrep),order(nrep))
    offset=0
    do i=1,nt
      j=buffers(i)%count
      if(j>0)then
        list_s(offset+1:offset+j)=buffers(i)%states(:j)
        norms(offset+1:offset+j)=buffers(i)%norms(:j)
      end if
      offset=offset+j
    end do
    deallocate(buffers,visited)
    !$omp parallel do private(i)
    do i=1,nrep
      order(i)=i
    end do
    call qsort_w_order(list_s,order,1,nrep)
    !$omp parallel do private(i)
    do i=1,nrep
      list_r(i)=norms(order(i))
    end do
    if(nrep>1)then
      if(any(list_s(2:)<=list_s(:nrep-1)))error stop 'orbit representatives not strictly ordered'
    end if
    THS=nrep
    basis_dump_setting=''
    call get_environment_variable('ED2_CUBIC_BASIS_PROOF_DUMP', &
      basis_dump_setting,status=basis_dump_status)
    if(basis_dump_status==0.and.trim(adjustl(basis_dump_setting))=='1')then
      open(newunit=basis_unit,file='cubic_basis_proof_dump.txt', &
        status='replace',action='write',iostat=basis_ios)
      if(basis_ios/=0)error stop 'Failed to open cubic basis proof dump'
      do basis_i=1,nrep
        write(basis_unit,'(i0,1x,z16.16)',iostat=basis_ios) &
          list_s(basis_i),transfer(list_r(basis_i),0_int64)
        if(basis_ios/=0)error stop 'Failed to write cubic basis proof dump'
      end do
      close(basis_unit,iostat=basis_ios)
      if(basis_ios/=0)error stop 'Failed to close cubic basis proof dump'
      write(*,'(a,i0)') 'Basis proof dump: exactly ranked states and bitwise norms=',nrep
    end if
    if(detail_timing)then
      detail_finish=omp_get_wtime()
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE preflight_other wall_s=', &
        orbit_begin-detail_begin-validation_sec
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE group_validation wall_s=',validation_sec
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE orbit_scan wall_s=',orbit_finish-orbit_begin
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE finalize wall_s=',detail_finish-orbit_finish
      write(*,'(a,es23.15)') 'ED2_BASIS_PHASE total wall_s=',detail_finish-detail_begin
    end if
    write(*,'(a,i0)') 'Basis builder: atomic-orbits, threads=',actual_threads
    used=.true.
  end subroutine try_orbit_basis

  ! Thread-private enumeration with the same generator and phase-sum order as
  ! checkstate_SQ. Noncanonical candidates MUST NOT publish visited bits.
  ! A canonical orbit with zero stabilizer character may be marked: none of
  ! its members contributes a representative in this symmetry sector.
  subroutine collect_canonical_orbit(s,ns,np,ranks,ro,is_representative)
    use input_param, only: NOS,NODmax,NODmin,L1,L2,L3,L4,L5,L6, &
      rk1,rk2,rk3,rk4,rk5,rk6,shift_1,shift_2,shift_3,shift_4,shift_5,shift_6
    integer, intent(in) :: s
    integer, intent(out) :: ns(:),np(:),ranks(:)
    real(8), intent(out) :: ro
    logical, intent(out) :: is_representative
    integer :: i1,i2,i3,i4,i5,i6,l,count,g
    real(8) :: p1,p2,p3,p4,p5
    complex(8) :: coef
    is_representative=.false.;ro=0d0;count=0;g=0;coef=(0d0,0d0)
    ns=list_fly(s,NODmax,NODmin,NOS)
    if(inv_list(ns,NODmax)/=s)error stop 'invalid rank/unrank in orbit basis'
    np=ns
    do i6=1,L6
      np=shift_6(np);p5=rk6*dble(i6)
      do i5=1,L5
        np=shift_5(np);p4=rk5*dble(i5)+p5
        do i4=1,L4
          np=shift_4(np);p3=rk4*dble(i4)+p4
          do i3=1,L3
            np=shift_3(np);p2=rk3*dble(i3)+p3
            do i2=1,L2
              np=shift_2(np);p1=rk2*dble(i2)+p2
              do i1=1,L1
                np=shift_1(np)
                call insertion_sort(np,NODmax)
                do l=NODmax,1,-1
                  if(np(l)<ns(l))return
                  if(np(l)>ns(l))exit
                end do
                g=g+1;ranks(g)=inv_list(np,NODmax)
                if(l==0)then
                  coef=coef+exp((0d0,-1d0)*(rk1*dble(i1)+p1))
                  count=count+1
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    if(count<1.or.g/=size(ranks))error stop 'orbit missing identity or wrong order'
    ro=(abs(coef)**2)*dble(size(ranks))/dble(count)
    is_representative=.true.
  end subroutine collect_canonical_orbit

  subroutine append_orbit(buffer,state,norm,dim)
    use iso_fortran_env, only: int64
    type(orbit_buffer), intent(inout) :: buffer
    integer, intent(in) :: state,dim
    real(8), intent(in) :: norm
    integer :: capacity,n
    integer, allocatable :: grown_s(:)
    real(8), allocatable :: grown_r(:)
    n=buffer%count
    if(n==size(buffer%states))then
      capacity=int(min(int(dim,int64),2_int64*size(buffer%states)))
      if(capacity<=n)error stop 'orbit buffer exceeds raw basis'
      allocate(grown_s(capacity),grown_r(capacity))
      grown_s(:n)=buffer%states;grown_r(:n)=buffer%norms
      call move_alloc(grown_s,buffer%states);call move_alloc(grown_r,buffer%norms)
    end if
    buffer%count=n+1;buffer%states(n+1)=state;buffer%norms(n+1)=norm
  end subroutine append_orbit
end module state_lists
