module eigen_solver
  use blas_lapack_interfaces, only: zheevr, dlamch, zcopy, check_lapack_info, zheevr_workspace_lengths
  use ham2vec
  implicit none
  private :: apply_hamiltonian
  integer :: NOE,NOK, NOM, maxitr, ite_lancz, i_vec_min, i_vec_max
  integer :: dim
  real(8) :: lnc_ene_conv
  namelist /input_TRLan/ NOE, NOK, NOM, maxitr, lnc_ene_conv, i_vec_min, i_vec_max
contains
  subroutine apply_hamiltonian(n,x,y)
    use input_param, only: NODmax,NODmin,list_s,list_r,explist
    integer, intent(in) :: n
    complex(8), intent(in) :: x(n)
    complex(8), intent(out) :: y(n)
    call ham_to_vec_wave_vector(y,x,n,NODmax,NODmin,list_s,list_r,explist)
  end subroutine apply_hamiltonian

  subroutine my_trlanczos_routines_complex(ene,psi,THS) 
    use input_param, only : OUTDIR
    real(8), allocatable, intent(out) :: ene(:)
    complex(8),allocatable, intent(out) :: psi(:,:)
    integer, intent(in) :: THS 
    integer :: i
    read(*,input_TRLan)
    write(*,input_TRLan)

    if(NOE > NOK) error stop 'NOE must not exceed NOK'
    if(NOM < 2*NOK) NOM = 2*NOK
    NOM=min(NOM,THS)
    NOK=min(NOK,THS)
    dim = THS
    if(NOE<1 .or. NOE>dim .or. NOK<NOE) error stop 'invalid TRLan eigenpair count'
    allocate( psi(dim,NOM+1), ene(NOM) )
    write(*,*) "Reduced THS=", THS 
    write(*,'("************************************************************************************")')
    write(*,*) "Eigen solver: my thick restart Lanczos"
    call lancz_eigen_val_complex(ene, psi) 
    write(*,'("***********************************************************************************")')
    write(*,*) "total lanczos step:", ite_lancz
    write(*,'("***********************************************************************************")')
    write(*,'("    i           E(i)")')
    do i = 1, NOE
      write(*,'(i5, 1es25.15)') i, ene(i)
    end do

    !open(10,file=trim(adjustl(OUTDIR))//'energy.dat',position='append')
    !do i = 1, NOE
    !  write(10,'(i5, 1es25.15)') i, ene(i)
    !end do
    !close(10)

    write(*,'("***********************************************************************************")')
    return
  end subroutine my_trlanczos_routines_complex

  subroutine read_input_TRLan 
    read(*,input_TRLan)
    write(*,input_TRLan)
  end subroutine read_input_TRLan

  subroutine Full_diag_routines(ene,psi,THS) 
    real(8),allocatable,intent(out)::ene(:)
    complex(8),allocatable,intent(out)::psi(:,:)
    integer,intent(in)::THS
    integer :: i
    dim=THS
    allocate(ene(dim),psi(dim,dim))
    !$omp parallel do
    do i = 1, dim
      psi(:,i) = (0.0d0,0.0d0)
      ene(i)=0.0d0
    end do
    !
    print *, "make ham"
    call make_full_hamiltonian(dim,psi)
    !
    print *, "diag starts"
    call fulldia_routines_with_momentum(dim,psi,ene,'V')

  end subroutine Full_diag_routines

  subroutine lancz_eigen_val_complex(ene, psi)
    use input_param, only: NODmax, NODmin, list_s, list_r, explist
    use krylov_solvers, only: restarted_hermitian
    real(8), intent(inout) :: ene(1:NOM)
    complex(8), intent(inout) :: psi(1:dim,1:NOM+1)
    integer :: status
    real(8) :: residuals(NOE)
    character(120) :: reason
    call restarted_hermitian(apply_hamiltonian,dim,NOE,NOK,NOM,lnc_ene_conv,maxitr,ene,psi(:,1:NOE), &
                              ite_lancz,status,reason,residuals)
    write(*,'(a,a)') 'TRLan status: ',trim(reason)
    write(*,'(a,100es25.15)') 'Eigenpair residual norms: ',residuals
    if(status/=0) error stop 'TRLan did not converge'
  end subroutine lancz_eigen_val_complex

  subroutine fulldia_routines_with_momentum(n,a,w,jobz) 
    use input_param, only : write_index_ene
    implicit none
    integer, intent(in) :: n
    real(8),intent(inout) :: w(n)
    complex(8),intent(inout) :: a(n,n)
    character(1), intent(in) :: jobz
    complex(8), allocatable :: z(:,:)
    !
    integer :: m, lda, il, iu, ldz, info
    integer, allocatable :: isuppz(:)
    real(8) :: vl, vu, abstol
    complex(8), allocatable :: work(:)
    integer :: lwork
    real(8), allocatable :: rwork(:)
    integer :: lrwork
    integer, allocatable :: iwork(:)
    integer :: liwork
    !
    complex(8) :: ztmp(1)
    real(8) :: dtmp(1)
    integer :: itmp(1)
    !
    if(jobz=='N')then
      allocate(z(1,1))
      ldz= 1
    else if(jobz=='V')then
      allocate(z(n,n))
      ldz= n
    else
      stop 888
    end if
    lda = n
    vl = -999.0d0
    vu =  999.0d0
    il = 1
    iu = n
    abstol = dlamch('S')
    allocate(isuppz(2*n))
    !     Use routine workspace query to get optimal workspace.
    lwork = -1
    lrwork = -1
    liwork = -1
    call zheevr(jobz, 'A', 'U', n, a, lda, vl, vu, il, iu, abstol, &
      m, w, z, ldz, isuppz, ztmp, lwork, dtmp, lrwork, itmp, liwork, info)

    call check_lapack_info('ZHEEVR (dense workspace query)', info)
    call zheevr_workspace_lengths(n, ztmp(1), dtmp(1), itmp(1), lwork, lrwork, liwork)
    allocate (work(lwork))
    allocate (rwork(lrwork))
    allocate (iwork(liwork))

    !     Solve the Hermitian eigenvalue problem
    call zheevr(jobz, 'A', 'U', n, a, lda, vl, vu, il, iu, abstol, &
      m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)

    call check_lapack_info('ZHEEVR (dense solve)', info)
    print *, "Eigenvalues by Full Diagonalization"
    call write_index_ene(w,n)
    if(jobz=='V') call zcopy(n*n,z,1,a,1)
  end subroutine fulldia_routines_with_momentum

end module eigen_solver
