module state_lists
  implicit none
  integer, allocatable :: cum_ltd_rep_combination(:,:)
  integer, allocatable :: ltd_rep_combination(:,:)
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
    real(8) :: r
    integer :: num = 1
    real(8), allocatable :: tmp_list_r(:)
    integer, allocatable :: order(:)
    real(8), allocatable :: tmp_list_r_pre(:)
    integer, allocatable :: list_s_pre(:)
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
    integer, intent(out) :: b, bmin
    integer, intent(in) :: s,ls,le
    integer, intent(in) :: list_s(ls:le)
    integer :: bmax
    bmin = ls; bmax = le
    do
       b = bmin + (bmax-bmin)/2
       if(s<list_s(b))then
          bmax = b-1
       else if(list_s(b)<s)then
          bmin = b+1
       else
          bmin = b
          return
       end if
       if(bmin>bmax)then
          b = -1
          return
       end if
    end do
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
end module state_lists

