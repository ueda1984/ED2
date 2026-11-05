module ham2vec
  implicit none
  complex(8), parameter :: zu = (1.0d0,0.0d0), zi = (0.0d0,1.0d0), zz = (0.0d0,0.0d0)
  real(8), parameter :: eps = 1.0d-15
  integer, allocatable :: wk_loc(:,:)
  complex(8), allocatable :: wk_ele(:,:)
  !
contains

  subroutine make_wk_loc_and_ele(lv)
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly, a_down_spin, &
      c_down_spin, a2_down_spin, c2_down_spin
    use input_param, only: NODmax,NODmin,NO_one,NO_two,p_one,p_two,list_s,list_r,explist, &
      Jint,hvec,NOS,local_NODmax,local_spin,wk_dim,MNTE
    implicit none
    integer, intent(in)::lv
    integer :: s, id, wk_dim1
    integer :: n0, n1, n2, i1, i2, dN
    integer :: i, j, ell(6), jd, max_jd
    integer, allocatable :: ni(:), st_list(:)
    complex(8) :: c
    !
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
    allocate(wk_loc(wk_dim1,wk_dim))
    allocate(wk_ele(wk_dim1,wk_dim))
    !
    allocate(ni(NODmax),st_list(NODmax))
    dN = NODmax-NODmin
    !
    max_jd = 1
    !
    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,st_list,c,jd) reduction(max:max_jd)
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
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            if(jd > wk_dim1) stop "MNTE is too small!!"
            wk_loc(jd,i) = id
            wk_ele(jd,i) = c &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1) * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            jd = jd + 1
          end if
        end if
        ! + term
        c = 0.5d0 * ( hvec(1,j) - zi * hvec(2,j) )
        if( n0>0 .and. n1<local_NODmax(i1) .and. abs(c) > eps )then
          ni = c_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            if(jd > wk_dim1) stop "MNTE is too small!!"
            wk_loc(jd,i) = id
            wk_ele(jd,i) = c &
              * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)) &
              * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
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
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1)*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              jd = jd + 1
            end if
          end if
          ! -- term
          c = 0.25d0*(Jint(1,j) - Jint(2,j)) + zi * 0.5d0 * Jint(9,j)
          if(n0+1<dN .and. n2>0 .and. abs(c) > eps)then
            ni = a2_down_spin(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(2.0d0*local_spin(i2)-n2+1)*n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              jd = jd + 1
            end if
          end if
          ! -+ term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) - zi*0.5d0*Jint(6,j)
          if(n2<local_NODmax(i2) .and. abs(c) > eps)then
            ni = j_flip_ni(i1,i2,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i1)-n1+1)*n1*(n2+1)*(2.0d0*local_spin(i2)-n2) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              jd = jd + 1
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
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2)*(local_spin(i1)-n1) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              jd = jd + 1
            end if
          end if
          ! +- term
          c = 0.25d0*(Jint(1,j)+Jint(2,j)) + zi*0.5d0*Jint(6,j)
          if(n1<local_NODmax(i1) .and. abs(c) > eps)then
            ni = j_flip_ni(i2,i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((2.0d0*local_spin(i2)-n2+1)*n2*(n1+1)*(2.0d0*local_spin(i1)-n1) ) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              jd = jd + 1
            end if
          end if
        end if
        if( n0>0 )then
          ! +z term
          c = 0.5d0 * ( -zi * ( Jint(4,j) + Jint(7,j) ) - Jint(5,j) + Jint(8,j) )
          if( n1<local_NODmax(i1) .and. abs(c) > eps )then
            ni = c_down_spin(i1,st_list,NODmax)
            call representative_SQ(NODmax,s,ell,ni)
            call findstate(s,list_s,id,lv)
            if(id > 0)then
              if(jd > wk_dim1) stop "MNTE is too small!!"
              wk_loc(jd,i) = id
              wk_ele(jd,i) = c &
                * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1))*(local_spin(i2)-n2) &
                * list_r(id) / list_r(i) &
                * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
              jd = jd + 1
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
                if(jd > wk_dim1) stop "MNTE is too small!!"
                wk_loc(jd,i) = id
                wk_ele(jd,i) = c &
                  * sqrt((n2+1)*(2.0d0*local_spin(i2)-n2))*(local_spin(i1)-n1) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
                jd = jd + 1
              end if
            end if
            ! ++ term
            c = 0.25d0 * ( Jint(1,j) - Jint(2,j) ) - zi * 0.5d0 * Jint(9,j)
            if( n1<local_NODmax(i1) .and. n0>1 .and. abs(c) > eps )then
              ni = c2_down_spin(i1,i2,st_list,NODmax)
              call representative_SQ(NODmax,s,ell,ni)
              call findstate(s,list_s,id,lv)
              if(id > 0)then
                if(jd > wk_dim1) stop "MNTE is too small!!"
                wk_loc(jd,i) = id
                wk_ele(jd,i) = c &
                  * sqrt((n1+1)*(2.0d0*local_spin(i1)-n1)*(n2+1)*(2.0d0*local_spin(i2)-n2)) &
                  * list_r(id) / list_r(i) &
                  * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
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

    return
  end subroutine make_wk_loc_and_ele

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
