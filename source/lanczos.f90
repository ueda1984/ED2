module lanczos
  use input_param, only: real8_is_finite
  use blas_lapack_interfaces, only: dlamch, dznrm2, zdotc, dstevx, &
    zdscal, zcopy, zaxpy
  use ham2vec
  implicit none
  integer :: ite_lancz, maxitr, dim, minitr, itrint
  real(8) :: lnc_ene_conv
  namelist /input_lancz/ lnc_ene_conv, maxitr, minitr, itrint
contains
  subroutine read_lanczos_para 
    read(*,input_lancz)
    if(.not.real8_is_finite(lnc_ene_conv)) &
      error stop 'Lanczos tolerance lnc_ene_conv must be finite'
    write(*,input_lancz)
  end subroutine read_lanczos_para
  !
  subroutine lanczos_routines(ene,THS,psi) 
    use input_param, only: NODmax, NODmin, list_s, list_r, explist 
    real(8), allocatable, intent(out) :: ene(:)
    integer, intent(in) :: THS
    complex(8), allocatable, intent(out) :: psi(:,:)
    complex(8), allocatable :: xpsi(:,:)
    real(8) :: alpha(2), residual_norm
    dim = THS
    allocate(ene(2))
    ene = 99.0d0
    read(*,input_lancz)
    if(.not.real8_is_finite(lnc_ene_conv)) &
      error stop 'Lanczos tolerance lnc_ene_conv must be finite'
    write(*,input_lancz)
    write(*,'("************************************************************************************")')
    write(*,*) "Eigen solver: Lanczos"
    allocate(psi(dim,1))
    call lancz_eigen_val(ene(1),psi)
    write(*,'("************************************************************************************")')
    write(*,*) "total lanczos step:", ite_lancz
    write(*,'("************************************************************************************")')
    write(*,'("E_0 = ",es25.15)') ene(1)
    write(*,'("************************************************************************************")')
    !check_wave_vector
    allocate(xpsi(dim,2))
    call ham_to_vec_wave_vector(xpsi(:,1),psi(:,1),dim,NODmax,NODmin,list_s,list_r,explist) 
    call ham_to_vec_wave_vector(xpsi(:,2),xpsi(:,1),dim,NODmax,NODmin,list_s,list_r,explist) 
    alpha(1) = dble(zdotc(dim,psi(1,1),1,xpsi(1,1),1))
    alpha(2) = dble(zdotc(dim,psi(1,1),1,xpsi(1,2),1))
    write(*,'("************************************************************************************")')
    write(*,'("quality of eigenvector")')
    write(*,'("                   <GS| H |GS>  = ", es25.15)') alpha(1)
    if(abs(alpha(2)) > tiny(1.0d0))then
      write(*,'("|1 - (<GS|H|GS>)^2/<GS|H^2|GS>| = ", es25.15)') abs(1.0d0-(alpha(1)**2)/alpha(2))
    else
      write(*,*) "Relative H^2 diagnostic undefined at <H^2>=0; use the residual norm."
    end if
    xpsi(:,2) = xpsi(:,1) - ene(1)*psi(:,1)
    residual_norm = dznrm2(dim,xpsi(1,2),1)
    write(*,'("Residual norm ||H psi - E psi|| = ", es25.15)') residual_norm
    write(*,'("************************************************************************************")')
    return
  end subroutine lanczos_routines
  !
  subroutine lancz_eigen_val(ene,psi_io)
    use input_param, only: NODmax, NODmin, list_s, list_r, explist 
    real(8), intent(out) :: ene
    complex(8), intent(out) :: psi_io(dim)
    complex(8), allocatable :: psi(:,:)
    real(8), allocatable :: alpha(:), beta(:), z(:,:)
    integer :: i, m, info, ell, il = 1, iu = 1, mitr
    real(8) :: tmp(2), abstol, conv, hnorm, observed_scale, vl = 0.0d0, vu = 0.0d0
    logical :: breakdown, have_previous_ritz, small_centered_residual
    integer, allocatable :: iwork(:), ifail(:)
    real(8), allocatable :: work(:), w(:), d(:), e(:)
    abstol = 2.0d0*dlamch('S')
    if(maxitr<1.or.minitr<1.or.itrint<1.or.lnc_ene_conv<=0.0d0) &
      error stop 'Invalid Lanczos iteration counts or tolerance'
    if(dim<1)error stop 'Lanczos requires a nonempty Hilbert space'
    ite_lancz=0
    allocate( psi(dim,3), alpha(maxitr), beta(0:maxitr) )
    !$omp parallel do private(i,tmp)
    do i = 1, dim
      call random_number(tmp)
      psi(i,1) = cmplx(tmp(1),tmp(2),8)
      psi(i,1) = psi(i,1) - (0.5d0, 0.5d0)
    end do
    beta(0) = dznrm2( dim, psi(1,1), 1 )
    if(beta(0)<=tiny(1.0d0))error stop 'Lanczos starting vector has zero norm'
    call zdscal( dim, 1.0d0 / beta(0), psi(1,1), 1 )
    call zcopy(dim,psi(1,1),1,psi_io,1)
    tmp = 0.0d0
    have_previous_ritz = .false.
    mitr= min(maxitr,dim)
    do i = 1, mitr
      call ham_to_vec_wave_vector(psi(:,2),psi(:,1),dim,NODmax,NODmin,list_s,list_r,explist) 
      hnorm = dznrm2(dim,psi(1,2),1)
      alpha(i) = dble(zdotc(dim,psi(1,1),1,psi(1,2),1))
      if(i > 1)then
        !$omp parallel do
        do ell = 1, dim
          psi(ell,2) = psi(ell,2) - alpha(i)*psi(ell,1) - beta(i-1)*psi(ell,3)
        end do
      else
        call zaxpy( dim, (-1.0d0,0.0d0)*alpha(i), psi(1,1), 1, psi(1,2), 1 )
      end if
      beta(i) = dznrm2( dim, psi(1,2), 1 )
      ! Do not normalize a zero or numerically unresolved residual. A small
      ! value relative to ||H q|| alone cannot certify an invariant subspace:
      ! a large scalar energy offset can conceal resolved spectral structure.
      breakdown = beta(i) <= 32.0d0*epsilon(1.0d0)*max(hnorm,tiny(1.0d0))
      if(breakdown.and.i<dim.and.beta(i)>0.0d0)then
        small_centered_residual = .false.
        if(i>1)then
          ! After observing more than one direction, retain ordinary finite
          ! spectral closure only if the residual is small on a scale that
          ! does not grow under H -> H + C I. This is a local closure check,
          ! not a general accuracy certificate for the final eigenvector.
          observed_scale = max(maxval(abs(alpha(1:i)-alpha(1))),maxval(beta(1:i-1)))
          small_centered_residual = beta(i) <= &
            32.0d0*epsilon(1.0d0)*max(observed_scale,tiny(1.0d0))
        end if
        if(.not.small_centered_residual)then
          write(*,*) 'Remove an avoidable constant energy offset or use an independent solver.'
          write(*,*) 'ALG=3 is an option only when the dense matrix fits in memory.'
          error stop 'Lanczos numerical near-breakdown cannot certify an invariant subspace'
        end if
      end if
      call zcopy(dim,psi(1,1),1,psi(1,3),1)
      if(.not.breakdown.and.i<mitr)then
        !$omp parallel do
        do ell = 1, dim
          psi(ell,1) = psi(ell,2) / beta(i)
        end do
      end if
      if(i==mitr.or.breakdown.or.(i>1.and.i>=minitr.and.mod(i,itrint)==0))then
        allocate(w(i),d(i),e(i),iwork(5*i), ifail(i), work(5*i), z(i,1))
        if(i == 1)then
           if(.not.breakdown.and.dim>1) &
             error stop 'Lanczos iteration limit reached before convergence'
           write(*,*) "### one-dimensional invariant Krylov subspace"
           ene = alpha(1)
           ! psi_io already contains the complete normalized start vector.
           ite_lancz = 1
           conv = 0.0d0
           write(*,'("itr, ene, conv: ", i5, es25.15, es11.1)') i, ene, conv
           return
        end if
        d = alpha(1:i); e = beta(1:i)
        call dstevx('V', 'I', i, d, e, vl, vu, il, iu, abstol, &
          m, w, z, size(z,1), work, iwork, ifail, info)
        if(info /= 0)then
          write(*,*) "error in dstevx in lanczos.f90. info =", info
          error stop 'Lanczos projected eigensolver failed'
        end if
        ene = w(1)
        ! A first Ritz estimate has no predecessor and cannot certify convergence.
        if(have_previous_ritz)then
          conv = abs(tmp(1)-ene)/max(abs(tmp(1)),abs(ene),tiny(1.0d0))
        else
          conv = 1.0d0 ! explicit nonconvergence marker, not a measured change
        end if
        write(*,'("itr, ene, conv: ", i5, es25.15, es11.1)') i, ene, conv
        if(i>=mitr.and.mitr<dim.and..not.breakdown.and. &
          (.not.have_previous_ritz.or.conv>=lnc_ene_conv)) &
          error stop 'Lanczos iteration limit reached before convergence'
        if(i>=mitr.or.breakdown.or. &
          (have_previous_ritz.and.conv<lnc_ene_conv))then
          exit
        end if
        tmp(1) = ene
        have_previous_ritz = .true.
        deallocate(w, d, e, iwork, ifail, work, z)
      end if
    end do
    ite_lancz = i
    !
    ! calcu. eigenvector
    !
    call zcopy(dim,psi_io,1,psi(1,1),1)
    call zdscal(dim, z(1,1), psi_io, 1)
    tmp = 0.0d0
    do i = 1, ite_lancz-1
      call ham_to_vec_wave_vector(psi(:,2),psi(:,1),dim,NODmax,NODmin,list_s,list_r,explist) 
      if(i > 1)then
        !$omp parallel do
        do ell = 1, dim
          psi(ell,2) = psi(ell,2) - alpha(i)*psi(ell,1) - beta(i-1)*psi(ell,3)
        end do
      else
        !$omp parallel do
        do ell = 1, dim
          psi(ell,2) = psi(ell,2) - alpha(i)*psi(ell,1)
        end do
      end if
      call zcopy(dim,psi(1,1),1,psi(1,3),1)
      !$omp parallel do
      do ell = 1, dim
        psi(ell,1) = psi(ell,2) / beta(i)
      end do
      call zaxpy( dim, (1.0d0,0.0d0)*z(i+1,1), psi(1,1), 1, psi_io, 1 )
    end do
    beta(0) = dznrm2( dim, psi_io, 1 )
    call zdscal( dim, 1.0d0 / beta(0), psi_io, 1 )    
    return
  end subroutine lancz_eigen_val
  !
end module lanczos
