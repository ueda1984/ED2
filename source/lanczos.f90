module lanczos
  use blas_lapack_interfaces, only: dznrm2, zdotc
  use ham2vec
  use krylov_observations, only: krylov_diagnostics,read_krylov_controls,write_krylov_diagnostics
  use krylov_residual_policy, only: residual_policy
  use krylov_observations, only: read_residual_controls
  implicit none
  private :: apply_hamiltonian
  integer :: ite_lancz, maxitr, dim, minitr, itrint
  real(8) :: lnc_ene_conv
  namelist /input_lancz/ lnc_ene_conv, maxitr, minitr, itrint
contains
  subroutine apply_hamiltonian(n,x,y)
    use input_param, only: NODmax,NODmin,list_s,list_r,explist
    integer, intent(in) :: n
    complex(8), intent(in) :: x(n)
    complex(8), intent(out) :: y(n)
    call ham_to_vec_wave_vector(y,x,n,NODmax,NODmin,list_s,list_r,explist)
  end subroutine apply_hamiltonian

  subroutine read_lanczos_para 
    read(*,input_lancz)
    write(*,input_lancz)
  end subroutine read_lanczos_para
  !
  subroutine lanczos_routines(ene,THS,psi) 
    use input_param, only: NODmax, NODmin, list_s, list_r, explist 
    real(8), allocatable, intent(out) :: ene(:)
    integer, intent(in) :: THS
    complex(8), allocatable, intent(out) :: psi(:,:)
    complex(8), allocatable :: xpsi(:,:)
    real(8) :: alpha(2)
    dim = THS
    allocate(ene(2))
    ene = 99.0d0
    read(*,input_lancz)
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
    allocate(xpsi(dim,1))
    call ham_to_vec_wave_vector(xpsi(:,1),psi(:,1),dim,NODmax,NODmin,list_s,list_r,explist) 
    alpha(1) = dble(zdotc(dim,psi(1,1),1,xpsi(1,1),1))
    alpha(2) = dznrm2(dim,xpsi(:,1),1)
    write(*,'("************************************************************************************")')
    write(*,'("quality of eigenvector")')
    write(*,'("                   <GS| H |GS>  = ", es25.15)') alpha(1)
    if(alpha(2)>0d0)then
      write(*,'("|1 - (<GS|H|GS>)^2/<GS|H^2|GS>| = ", es25.15)') abs(1d0-(alpha(1)/alpha(2))**2)
    else
      write(*,'("|1 - (<GS|H|GS>)^2/<GS|H^2|GS>| = ", es25.15)') 0d0
    end if
    write(*,'("************************************************************************************")')
    return
  end subroutine lanczos_routines
  !
  subroutine lancz_eigen_val(ene,psi_io)
    use input_param, only: NODmax, NODmin, list_s, list_r, explist
    use krylov_solvers, only: lanczos_ground
    real(8), intent(out) :: ene
    complex(8), intent(out) :: psi_io(dim)
    integer :: status
    real(8) :: residual
    character(120) :: reason
    type(krylov_diagnostics), allocatable :: observation
    type(residual_policy), allocatable :: policy
    complex(8), allocatable :: initial(:)
    integer, allocatable :: seed
    character(:), allocatable :: diagnostic_path
    call read_krylov_controls(dim,initial,seed,diagnostic_path)
    call read_residual_controls(policy)
    if(len(diagnostic_path)>0)allocate(observation)
    call lanczos_ground(apply_hamiltonian,dim,lnc_ene_conv,maxitr,minitr,itrint,ene,psi_io, &
                        ite_lancz,status,reason,residual,initial,observation,seed,policy)
    if(allocated(observation))call write_krylov_diagnostics(diagnostic_path,observation)
    write(*,'(a,a)') 'Lanczos status: ',trim(reason)
    write(*,'(a,es25.15)') 'Eigenpair residual norm: ',residual
    if(status/=0) error stop 'Lanczos did not converge'
  end subroutine lancz_eigen_val
end module lanczos
