module ham2vec
  implicit none
  complex(8), parameter :: zu = (1.0d0,0.0d0), zi = (0.0d0,1.0d0), zz = (0.0d0,0.0d0)
  real(8), parameter :: eps = 1.0d-15
  integer, allocatable :: wk_loc(:,:)
  complex(8), allocatable :: wk_ele(:,:)
  ! Optional compressed row-wise copy of fully cached COMPLEX Hamiltonian.
  ! Default is historical padded storage; no effect unless enabled.
  logical, save :: csr_complex_cache=.false.
  ! Separate optional release of historical padded cache after CSR packing.
  ! Keeping the arrays available is the default public module behavior.
  logical, save :: csr_release_padded=.false.
  ! Optional coefficient-certified upper bound on padded cache row width.
  logical, save :: cache_coeff_bound=.false.
  logical, save :: csr_direct_build=.false.
  integer, allocatable, save :: csr_rowptr(:),csr_target(:)
  complex(8), allocatable, save :: csr_coeff(:)
  !
contains

  subroutine configure_csr_complex_cache(lv)
    use iso_fortran_env, only: int64
    use input_param, only: wk_dim
    integer, intent(in) :: lv
    integer :: status,release_status,istat,i,j,nnz,p
    integer(int64) :: entries,extra_bytes,old_bytes
    character(len=16) :: setting,release_setting
    setting=''
    call get_environment_variable('ED2_CSR_COMPLEX_CACHE',setting,status=status)
    if(status==0)then
      if(trim(adjustl(setting))/='0'.and.trim(adjustl(setting))/='1') &
        error stop 'ED2_CSR_COMPLEX_CACHE must be 0 or 1'
    else if(status/=1)then
      error stop 'ED2_CSR_COMPLEX_CACHE flag too long or invalid'
    end if
    csr_complex_cache=status==0.and.trim(adjustl(setting))=='1'
    release_setting=''
    call get_environment_variable('ED2_CSR_RELEASE_PADDED_CACHE', &
      release_setting,status=release_status)
    if(release_status==0)then
      if(trim(adjustl(release_setting))/='0'.and. &
         trim(adjustl(release_setting))/='1') &
        error stop 'ED2_CSR_RELEASE_PADDED_CACHE must be 0 or 1'
    else if(release_status/=1)then
      error stop 'ED2_CSR_RELEASE_PADDED_CACHE flag too long or invalid'
    end if
    csr_release_padded=release_status==0.and.trim(adjustl(release_setting))=='1'
    if(csr_release_padded.and..not.csr_complex_cache) &
      error stop 'ED2_CSR_RELEASE_PADDED_CACHE requires ED2_CSR_COMPLEX_CACHE=1'
    if(allocated(csr_rowptr))deallocate(csr_rowptr)
    if(allocated(csr_target))deallocate(csr_target)
    if(allocated(csr_coeff))deallocate(csr_coeff)
    if(.not.csr_complex_cache)return
    if(wk_dim/=lv)error stop 'ED2_CSR_COMPLEX_CACHE requires fully cached matrix'

    ! Preserve all original non-sentinel terms, including explicit zeros.
    ! No sorting or floating-point reassociation is allowed.
    allocate(csr_rowptr(lv+1),stat=istat)
    if(istat/=0)error stop 'ED2_CSR_COMPLEX_CACHE row allocation failed'
    csr_rowptr(1)=1
    entries=0_int64
    do i=1,lv
      nnz=0
      do j=1,size(wk_loc,1)
        if(wk_loc(j,i)==0)exit
        if(wk_loc(j,i)<1.or.wk_loc(j,i)>lv) &
          error stop 'ED2_CSR_COMPLEX_CACHE invalid target'
        nnz=nnz+1
      end do
      entries=entries+int(nnz,int64)
      if(entries>=int(huge(i),int64)) &
        error stop 'ED2_CSR_COMPLEX_CACHE row index capacity exceeded'
      csr_rowptr(i+1)=int(entries)+1
    end do

    allocate(csr_target(int(entries)),csr_coeff(int(entries)),stat=istat)
    if(istat/=0)error stop 'ED2_CSR_COMPLEX_CACHE values allocation failed'
    do i=1,lv
      p=csr_rowptr(i)
      do j=1,size(wk_loc,1)
        if(wk_loc(j,i)==0)exit
        csr_target(p)=wk_loc(j,i)
        csr_coeff(p)=wk_ele(j,i)
        p=p+1
      end do
      if(p/=csr_rowptr(i+1))error stop 'ED2_CSR_COMPLEX_CACHE row count mismatch'
    end do
    extra_bytes=int(storage_size(csr_rowptr)/8,int64)*int(size(csr_rowptr),int64)+ &
      int(storage_size(csr_target)/8,int64)*int(size(csr_target),int64)+ &
      int(storage_size(csr_coeff)/8,int64)*int(size(csr_coeff),int64)
    write(*,'(a,i0,a,i0)') &
      'Hamiltonian action: csr-complex cache entries=',size(csr_coeff), &
      ', extra_bytes=',extra_bytes
    if(csr_release_padded)then
      old_bytes=int(storage_size(wk_loc)/8,int64)*int(size(wk_loc),int64)+ &
        int(storage_size(wk_ele)/8,int64)*int(size(wk_ele),int64)
      ! Peak construction memory still includes both representations.
      ! Deallocation reduces long-lived cached storage, not peak RSS.
      deallocate(wk_loc,wk_ele,stat=istat)
      if(istat/=0)error stop 'ED2_CSR_RELEASE_PADDED_CACHE deallocation failed'
      write(*,'(a,i0,a,i0)') &
        'Hamiltonian action: csr-padded-cache released_bytes=',old_bytes, &
        ', net_saved_bytes=',old_bytes-extra_bytes
    end if
  end subroutine configure_csr_complex_cache


  subroutine coefficient_certified_cache_width(width)
    use iso_fortran_env, only: int64
    use input_param, only: NO_one,NO_two,p_one,hvec,Jint
    integer, intent(out) :: width
    integer :: j,i1
    integer(int64) :: upper
    complex(8) :: c
    ! Each physical one-site and two-site operator can contribute at most
    ! ONE row term per applicable transition, irrespective of basis size.
    ! Exactly repeat the coefficient tests in make_wk_loc_and_ele.
    upper=1_int64 ! required diagonal record
    do j=1,NO_one
      i1=p_one(j)
      c=0.5d0*(hvec(1,i1)+zi*hvec(2,i1))
      if(abs(c)>eps) upper=upper+1_int64
      c=0.5d0*(hvec(1,i1)-zi*hvec(2,i1))
      if(abs(c)>eps) upper=upper+1_int64
    end do
    do j=1,NO_two
      c=0.5d0*(zi*(Jint(4,j)+Jint(7,j))-Jint(5,j)+Jint(8,j))
      if(abs(c)>eps) upper=upper+1_int64
      c=0.25d0*(Jint(1,j)-Jint(2,j))+zi*0.5d0*Jint(9,j)
      if(abs(c)>eps) upper=upper+1_int64
      c=0.25d0*(Jint(1,j)+Jint(2,j))-zi*0.5d0*Jint(6,j)
      if(abs(c)>eps) upper=upper+1_int64
      c=0.5d0*(-zi*(Jint(4,j)-Jint(7,j))+Jint(5,j)+Jint(8,j))
      if(abs(c)>eps) upper=upper+1_int64
      c=0.25d0*(Jint(1,j)+Jint(2,j))+zi*0.5d0*Jint(6,j)
      if(abs(c)>eps) upper=upper+1_int64
      c=0.5d0*(-zi*(Jint(4,j)+Jint(7,j))-Jint(5,j)+Jint(8,j))
      if(abs(c)>eps) upper=upper+1_int64
      c=0.5d0*(zi*(Jint(4,j)-Jint(7,j))+Jint(5,j)+Jint(8,j))
      if(abs(c)>eps) upper=upper+1_int64
      c=0.25d0*(Jint(1,j)-Jint(2,j))-zi*0.5d0*Jint(9,j)
      if(abs(c)>eps) upper=upper+1_int64
    end do
    if(upper>int(huge(width),int64)) &
      error stop 'ED2_CACHE_COEFF_BOUND upper bound exceeds index range'
    width=int(upper)
  end subroutine coefficient_certified_cache_width

  subroutine make_wk_loc_and_ele(lv)
    use state_lists, only: j_flip_ni, representative_SQ, findstate, list_fly, a_down_spin, &
      c_down_spin, a2_down_spin, c2_down_spin
    use input_param, only: NODmax,NODmin,NO_one,NO_two,p_one,p_two,list_s,list_r,explist, &
      Jint,hvec,NOS,local_NODmax,local_spin,wk_dim,MNTE
    use iso_fortran_env, only: int64
    implicit none
    integer, intent(in)::lv
    integer :: s, id, wk_dim1
    integer :: n0, n1, n2, i1, i2, dN
    integer :: i, j, ell(6), jd, max_jd,nominal_width,width_status
    integer :: direct_status,prereq_status,alloc_status
    integer :: pass,first,last,batch_rows,ib,batch_max,nnz,pos
    integer, parameter :: direct_batch_rows=64
    integer(int64) :: nominal_bytes,bounded_bytes,bytes_per_element
    integer(int64) :: entries,extra_bytes,stage_bytes,max_stage_bytes
    character(len=16) :: width_setting,direct_setting,prereq_setting
    integer, allocatable :: ni(:), st_list(:)
    complex(8) :: c
    !
    if(wk_dim < 0)then
      wk_dim = lv
    else
      wk_dim = min(wk_dim,lv)
    end if
    !
    width_setting=''
    call get_environment_variable('ED2_CACHE_COEFF_BOUND',width_setting,status=width_status)
    if(width_status==0)then
      if(trim(adjustl(width_setting))/='0'.and.trim(adjustl(width_setting))/='1') &
        error stop 'ED2_CACHE_COEFF_BOUND must be 0 or 1'
    else if(width_status/=1)then
      error stop 'ED2_CACHE_COEFF_BOUND flag too long or invalid'
    end if
    cache_coeff_bound=width_status==0.and.trim(adjustl(width_setting))=='1'
    if(MNTE < 0)then
      wk_dim1= 2*NO_one+8*NO_two+1
      if(cache_coeff_bound)then
        nominal_width=wk_dim1
        call coefficient_certified_cache_width(wk_dim1)
        if(wk_dim1<1.or.wk_dim1>nominal_width) &
          error stop 'ED2_CACHE_COEFF_BOUND invalid row width'
        bytes_per_element=int(storage_size(wk_dim1)/8,int64)+ &
          int(storage_size(c)/8,int64)
        nominal_bytes=bytes_per_element*int(nominal_width,int64)*int(wk_dim,int64)
        bounded_bytes=bytes_per_element*int(wk_dim1,int64)*int(wk_dim,int64)
        write(*,'(a,i0,a,i0,a,i0)') &
          'Hamiltonian cache: coefficient-bound width=',wk_dim1, &
          ', nominal_width=',nominal_width, &
          ', initial_bytes_saved=',nominal_bytes-bounded_bytes
      end if
    else
      wk_dim1= MNTE + 1
      if(cache_coeff_bound) &
        write(*,'(a)') 'Hamiltonian cache: explicit MNTE override preserved'
    end if
    !
    ! Explicit opt-in direct CSR build: keep only one small staging batch.
    ! The same row-generating Fortran code is called by both paths.
    direct_setting=''
    call get_environment_variable('ED2_CSR_DIRECT_BUILD',direct_setting,status=direct_status)
    if(direct_status==0)then
      if(trim(adjustl(direct_setting))/='0'.and.trim(adjustl(direct_setting))/='1') &
        error stop 'ED2_CSR_DIRECT_BUILD must be 0 or 1'
    else if(direct_status/=1)then
      error stop 'ED2_CSR_DIRECT_BUILD flag too long or invalid'
    end if
    csr_direct_build=direct_status==0.and.trim(adjustl(direct_setting))=='1'
    if(csr_direct_build)then
      prereq_setting=''
      call get_environment_variable('ED2_CSR_COMPLEX_CACHE',prereq_setting,status=prereq_status)
      if(prereq_status/=0.or.trim(adjustl(prereq_setting))/='1') &
        error stop 'ED2_CSR_DIRECT_BUILD requires ED2_CSR_COMPLEX_CACHE=1'
      prereq_setting=''
      call get_environment_variable('ED2_CSR_RELEASE_PADDED_CACHE',prereq_setting,status=prereq_status)
      if(prereq_status/=0.or.trim(adjustl(prereq_setting))/='1') &
        error stop 'ED2_CSR_DIRECT_BUILD requires ED2_CSR_RELEASE_PADDED_CACHE=1'
      if(wk_dim/=lv)error stop 'ED2_CSR_DIRECT_BUILD requires fully cached matrix'
      csr_complex_cache=.true.
      csr_release_padded=.true.
      if(allocated(csr_rowptr))deallocate(csr_rowptr)
      if(allocated(csr_target))deallocate(csr_target)
      if(allocated(csr_coeff))deallocate(csr_coeff)
    end if

    allocate(ni(NODmax),st_list(NODmax))
    dN = NODmax-NODmin
    max_jd=1
    if(csr_direct_build)then
      allocate(csr_rowptr(lv+1),stat=alloc_status)
      if(alloc_status/=0)error stop 'ED2_CSR_DIRECT_BUILD rowptr allocation failed'
      csr_rowptr(1)=1
      entries=0_int64
      max_stage_bytes=0_int64
      ! Pass 1 counts occupied records with no global padded matrix.
      ! Pass 2 allocates EXACT CSR capacity and generates the same rows.
      do pass=1,2
        do first=1,lv,direct_batch_rows
          last=min(lv,first+direct_batch_rows-1)
          batch_rows=last-first+1
          allocate(wk_loc(wk_dim1,batch_rows),wk_ele(wk_dim1,batch_rows),stat=alloc_status)
          if(alloc_status/=0)error stop 'ED2_CSR_DIRECT_BUILD batch allocation failed'
          stage_bytes=int(storage_size(wk_loc)/8,int64)*int(size(wk_loc),int64)+ &
            int(storage_size(wk_ele)/8,int64)*int(size(wk_ele),int64)
          max_stage_bytes=max(max_stage_bytes,stage_bytes)
          call fill_cache_batch(first,last,batch_max)
          max_jd=max(max_jd,batch_max)
          do i=first,last
            ib=i-first+1
            nnz=0
            do j=1,wk_dim1
              if(wk_loc(j,ib)==0)exit
              if(wk_loc(j,ib)<1.or.wk_loc(j,ib)>lv) &
                error stop 'ED2_CSR_DIRECT_BUILD invalid target index'
              nnz=nnz+1
            end do
            if(pass==1)then
              entries=entries+int(nnz,int64)
              if(entries>=int(huge(i),int64)) &
                error stop 'ED2_CSR_DIRECT_BUILD row offset capacity exceeded'
              csr_rowptr(i+1)=int(entries)+1
            else
              if(nnz/=csr_rowptr(i+1)-csr_rowptr(i)) &
                error stop 'ED2_CSR_DIRECT_BUILD two-pass row count mismatch'
              pos=csr_rowptr(i)
              do j=1,nnz
                csr_target(pos)=wk_loc(j,ib)
                csr_coeff(pos)=wk_ele(j,ib)
                pos=pos+1
              end do
              if(pos/=csr_rowptr(i+1)) &
                error stop 'ED2_CSR_DIRECT_BUILD row copy mismatch'
            end if
          end do
          deallocate(wk_loc,wk_ele)
        end do
        if(pass==1)then
          allocate(csr_target(int(entries)),csr_coeff(int(entries)),stat=alloc_status)
          if(alloc_status/=0)error stop 'ED2_CSR_DIRECT_BUILD CSR allocation failed'
        end if
      end do
      extra_bytes=int(storage_size(csr_rowptr)/8,int64)*int(size(csr_rowptr),int64)+ &
        int(storage_size(csr_target)/8,int64)*int(size(csr_target),int64)+ &
        int(storage_size(csr_coeff)/8,int64)*int(size(csr_coeff),int64)
      write(*,'(a,i0,a,i0)') &
        'Hamiltonian action: csr-complex cache entries=',size(csr_coeff), &
        ', extra_bytes=',extra_bytes
      write(*,'(a,i0,a,i0,a,i0)') &
        'Hamiltonian action: csr-direct-build passes=',2, &
        ', max_staging_bytes=',max_stage_bytes, &
        ', peak_cache_payload_upper_bytes=',extra_bytes+max_stage_bytes
    else
      allocate(wk_loc(wk_dim1,wk_dim))
      allocate(wk_ele(wk_dim1,wk_dim))
      call fill_cache_batch(1,wk_dim,max_jd)
    end if

    write(*,*) "Current MNTE = ", wk_dim1 - 1
    write(*,*) "Optimal MNTE = ", max_jd - 1

    if(.not.csr_direct_build)call configure_csr_complex_cache(lv)
    return
  contains
    subroutine fill_cache_batch(first,last,batch_max)
      integer, intent(in) :: first,last
      integer, intent(out) :: batch_max
      integer :: ib
      batch_max=1
    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,ib,st_list,c,jd) reduction(max:batch_max)
    do i = first, last
      ib = i - first + 1
      wk_loc(1,ib)  = i
      wk_loc(2:,ib) = 0
      wk_ele(:,ib) = zz
      jd = 2
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0 )
      !
      do j = 1, NO_one
        i1 = p_one(j)
        n1 = count(st_list == i1)
        wk_ele(1,ib) = wk_ele(1,ib) + hvec(3,i1) * (local_spin(i1)-n1)
        ! - term
        c = 0.5d0 * ( hvec(1,i1) + zi * hvec(2,i1) )
        if(n0 < dN .and. n1>0 .and. abs(c) > eps)then
          ni = a_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            if(jd > wk_dim1) stop "MNTE is too small!!"
            wk_loc(jd,ib) = id
            wk_ele(jd,ib) = c &
              * sqrt((2.0d0*local_spin(i1)-n1+1)*n1) * list_r(id) / list_r(i) &
              * explist(ell(1),ell(2),ell(3),ell(4),ell(5),ell(6))
            jd = jd + 1
          end if
        end if
        ! + term
        c = 0.5d0 * ( hvec(1,i1) - zi * hvec(2,i1) )
        if( n0>0 .and. n1<local_NODmax(i1) .and. abs(c) > eps )then
          ni = c_down_spin(i1,st_list,NODmax)
          call representative_SQ(NODmax,s,ell,ni)
          call findstate(s,list_s,id,lv)
          if(id > 0)then
            if(jd > wk_dim1) stop "MNTE is too small!!"
            wk_loc(jd,ib) = id
            wk_ele(jd,ib) = c &
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
        wk_ele(1,ib) = wk_ele(1,ib) + Jint(3,j) * &
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
              wk_loc(jd,ib) = id
              wk_ele(jd,ib) = c &
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
              wk_loc(jd,ib) = id
              wk_ele(jd,ib) = c &
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
              wk_loc(jd,ib) = id
              wk_ele(jd,ib) = c &
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
              wk_loc(jd,ib) = id
              wk_ele(jd,ib) = c &
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
              wk_loc(jd,ib) = id
              wk_ele(jd,ib) = c &
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
              wk_loc(jd,ib) = id
              wk_ele(jd,ib) = c &
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
                wk_loc(jd,ib) = id
                wk_ele(jd,ib) = c &
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
                wk_loc(jd,ib) = id
                wk_ele(jd,ib) = c &
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

    end subroutine fill_cache_batch
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
        i1 = p_one(j)
        n1 = count(st_list == i1)
        Ham(i,i) = Ham(i,i) + hvec(3,i1) * (local_spin(i1)-n1)
        ! - term
        c = 0.5d0 * ( hvec(1,i1) + zi * hvec(2,i1) )
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
        c = 0.5d0 * ( hvec(1,i1) - zi * hvec(2,i1) )
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
    
    ! Scratch state transitions are needed only for uncached rows.
    ! Fully cached H*x must not allocate these arrays per application.
    if(wk_dim<lv)allocate(ni(NODmax),st_list(NODmax))
    dN = NODmax - NODmin

    if(csr_complex_cache)then
      ! The padded array is intentionally absent in CSR release mode.
      wk_dim1=0
    else
      if(.not.allocated(wk_ele))error stop 'ED2 cached matrix unavailable'
      wk_dim1=size(wk_ele,1)
    end if

    if(csr_complex_cache)then
      !$omp parallel do private(i,j)
      do i = 1, wk_dim
        v0(i) = 0.0d0
        do j = csr_rowptr(i),csr_rowptr(i+1)-1
          v0(i) = v0(i) + csr_coeff(j) * v1(csr_target(j))
        end do
      end do
    else
      !$omp parallel do private(i,j)
      do i = 1, wk_dim
        v0(i) = 0.0d0
        do j = 1, wk_dim1
          if(wk_loc(j,i) == 0) exit
          v0(i) = v0(i) + wk_ele(j,i) * v1(wk_loc(j,i))
        end do
      end do
    end if
    !

    !$omp parallel do private(n0,n1,n2,i1,i2,s,ni,id,ell,i,j,st_list,c)
    do i = wk_dim+1, lv
      v0(i) = (0.0d0, 0.0d0)
      st_list = list_fly(list_s(i),NODmax,NODmin,NOS)
      n0 = count(st_list == 0 )
      !
      do j = 1, NO_one
        i1 = p_one(j)
        n1 = count(st_list == i1)
        v0(i) = v0(i) + hvec(3,i1) * (local_spin(i1)-n1) * v1(i)
        ! - term
        c = 0.5d0 * ( hvec(1,i1) + zi * hvec(2,i1) )
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
        c = 0.5d0 * ( hvec(1,i1) - zi * hvec(2,i1) )
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
