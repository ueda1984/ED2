module input_param
  use iso_fortran_env, only: int64
  implicit none
  integer :: NOS       ! number of sites
  integer :: NODmax    ! Upper limit value for number of down spins
  integer :: NODmin    ! Lower limit value for number of down spins
  integer :: NOV       ! Number of lowest eigenvectors used for computing static physical quantities.
  integer :: THS       ! Total Hilbelt space of the target space 
  integer :: NO_one    ! number of one-body interactions
  integer :: NO_two    ! number of two-body interactions
  integer :: wk_dim    ! dimensions of the work array
  integer :: MNTE      ! Maximum number of transitional elements for a representative state
  integer :: ALG       ! 1:Conventional Lanczos,  2:Thick-Restarted Lanczos, 3:Full diagonalization
  integer, allocatable :: p_one(:)    ! locations of one-body interactions: p_one(1:NO_one)
  integer, allocatable :: p_two(:,:)  ! locations of two-body interactions: p_two(2,1:NO_two)
  real(8), allocatable :: Jint(:,:)   ! Jint(1:9 [Jx,Jy,Jz,Dx,Dy,Dz],1:NO_two)
  real(8), allocatable :: hvec(:,:)   ! magnetic_field(1:3[hx,hy,hz],1:NO_one)
  integer, allocatable :: local_NODmax(:)  ! local_NODmax(i) : NODmax for ith site 
  real(8), allocatable :: local_spin(:)    ! local_spin(i)   : spin for ith site 
  !
  integer :: cal_lm ! flag for expectation values of local magnetizations
  !
  integer :: re_wf, wr_wf
  !
  character(100) :: FILE_xyz_dm_gamma, FILE_hvec, FILE_p1, FILE_p12, FILE_NODmax, FILE_spin, FILE_wf
  character(100) :: FILE_S1="", FILE_S2="", FILE_S3="", FILE_S4="", FILE_S5="", FILE_S6=""
  character(100) :: OUTDIR
  !
  integer :: L1=1, L2=1, L3=1, L4=1, L5=1, L6=1
  integer :: M1=0, M2=0, M3=0, M4=0, M5=0, M6=0
  real(8) :: PI, rk1, rk2, rk3, rk4, rk5, rk6
  integer, allocatable :: list_s(:)
  real(8), allocatable :: list_r(:)
  integer, allocatable :: shift_1(:), shift_2(:), shift_3(:), shift_4(:), shift_5(:), shift_6(:)
  complex(8), allocatable :: explist(:,:,:,:,:,:)
  !
  integer :: cal_cf     ! flag for expectation values of two-body correlation functions
  integer :: NO_two_cf  ! number of two-body correlation functions we want to calculate
  integer, allocatable :: p_two_cf(:,:)  
                        ! (i,j) pairs of two-body correlation functions: p_two_cf(2,1:NO_two_cf)
  logical, allocatable :: swap_list(:)
  character(100) :: FILE_two_cf
  !
  namelist /input_parameters/ NOS,NODmax,NODmin,L1,L2,L3,L4,L5,L6,M1,M2,M3,M4,M5,M6,&
    NO_one,NO_two,wk_dim,MNTE,ALG, FILE_xyz_dm_gamma, FILE_hvec, &
    OUTDIR, FILE_S1, FILE_S2, FILE_S3, FILE_S4, FILE_S5, FILE_S6, FILE_NODmax, &
    FILE_spin, FILE_wf, re_wf, wr_wf, cal_lm, cal_cf, NO_two_cf, FILE_two_cf
  !
  real(8), allocatable :: ene(:)
  complex(8),allocatable :: psi(:,:) 
  !
contains
  !
  elemental logical function real8_is_finite(value)
    real(8), intent(in) :: value
    ! Supported ED2 platforms use IEEE binary64 for real(8). Check its
    ! exponent bits because -Ofast may assume floating-point inputs finite
    ! and optimize away an IEEE_IS_FINITE test.
    real8_is_finite = iand(transfer(value,0_int64),int(z'7ff0000000000000',int64)) /= &
      int(z'7ff0000000000000',int64)
  end function real8_is_finite
  !
  subroutine read_ip 
    integer :: i,j,k,l,m,n
    integer :: tmp
    !
    write(*,'(" ### Read input_parameters. ")')
    read(*,input_parameters)
    write(*,input_parameters)
    if(NOS<1)error stop 'NOS must be positive'
    if(ALG<1.or.ALG>3)error stop 'ALG must be 1, 2, or 3'
    if(NO_one<0.or.NO_two<0)error stop 'Interaction counts must be nonnegative'
    if(NODmin<0.or.NODmax<NODmin)error stop 'Require 0 <= NODmin <= NODmax'
    if(any([L1,L2,L3,L4,L5,L6]<1))error stop 'Symmetry periods must be positive'
    if(cal_cf==1)then
      if(NO_two_cf<0)error stop 'NO_two_cf must be nonnegative'
    end if
    !
    write(*,'(" ### Set wavevectors. ")')
    PI = acos(-1.0d0)
    rk1 = 2.0d0 * PI * dble(M1)/dble(L1) 
    rk2 = 2.0d0 * PI * dble(M2)/dble(L2) 
    rk3 = 2.0d0 * PI * dble(M3)/dble(L3)
    rk4 = 2.0d0 * PI * dble(M4)/dble(L4) 
    rk5 = 2.0d0 * PI * dble(M5)/dble(L5) 
    rk6 = 2.0d0 * PI * dble(M6)/dble(L6) 
    !
    write(*,'(" ### Set random_seed. ")')
    call random_seed
    !
    write(*,'(" ### Allocation arrays. ")')
    allocate(p_one(NO_one),p_two(2,NO_two),Jint(9,NO_two),hvec(3,NO_one),shift_1(0:NOS),&
      shift_2(0:NOS),shift_3(0:NOS),shift_4(0:NOS),shift_5(0:NOS),shift_6(0:NOS),&
      explist(L1,L2,L3,L4,L5,L6),local_NODmax(NOS),local_spin(NOS))
    !
    write(*,'(" ### Set phase factors. ")')
    forall(i=1:L1,j=1:L2,k=1:L3,l=1:L4,m=1:L5,n=1:L6) explist(i,j,k,l,m,n) = &
      exp((0.0d0,1.0d0)*(rk1*dble(i)+rk2*dble(j)+rk3*dble(k)+rk4*dble(l)+rk5*dble(m)+rk6*dble(n)))
    !
    write(*,'(" ### Read FILE_xyz_dm_gamma. ")')
    open(10, file=trim(adjustl(FILE_xyz_dm_gamma)),status='old')
    do i = 1, NO_two
      read(10,*) p_two(1,i), p_two(2,i), Jint(1:9,i)
      if(.not.all(real8_is_finite(Jint(:,i)))) &
        error stop 'FILE_xyz_dm_gamma contains a nonfinite coupling'
      if(any(p_two(:,i)<1).or.any(p_two(:,i)>NOS).or.p_two(1,i)==p_two(2,i)) &
        error stop 'Bond indices must be distinct sites in 1..NOS'
      if( p_two(1,i) > p_two(2,i) )then
        Jint(4:6,i) = -1.0d0 * Jint(4:6,i)
        tmp = p_two(1,i)
        p_two(1,i) = p_two(2,i)
        p_two(2,i) = tmp
      end if
    end do
    close(10)
    !
    write(*,'(" ### Read FILE_hvec. ")')
    open(10, file=trim(adjustl(FILE_hvec)),status='old')
    do i = 1, NO_one
      read(10,*) p_one(i), hvec(1:3,i)
      if(.not.all(real8_is_finite(hvec(:,i)))) &
        error stop 'FILE_hvec contains a nonfinite field'
      if(p_one(i)<1.or.p_one(i)>NOS)error stop 'Field index must be a site in 1..NOS'
    end do
    close(10)
    !
    write(*,'(" ### Set translational operations. ")')
    if(L1==1)then
      call set_identity_ope(shift_1,NOS)
    else
      open(10, file=trim(adjustl(FILE_S1)),status='old')
      shift_1(0) = 0
      do i = 1, NOS
        read(10,*) shift_1(i)
      end do
      close(10)
    end if
    !
    if(L2==1)then
      call set_identity_ope(shift_2,NOS)
    else
      open(10, file=trim(adjustl(FILE_S2)),status='old')
      shift_2(0) = 0
      do i = 1, NOS
        read(10,*) shift_2(i)
      end do
      close(10)
    end if
    !
    if(L3==1)then
      call set_identity_ope(shift_3,NOS)
    else
      open(10, file=trim(adjustl(FILE_S3)),status='old')
      shift_3(0) = 0
      do i = 1, NOS
        read(10,*) shift_3(i)
      end do
      close(10)
    end if
    !
    if(L4==1)then
      call set_identity_ope(shift_4,NOS)
    else
      open(10, file=trim(adjustl(FILE_S4)),status='old')
      shift_4(0) = 0
      do i = 1, NOS
        read(10,*) shift_4(i)
      end do
      close(10)
    end if
    !
    if(L5==1)then
      call set_identity_ope(shift_5,NOS)
    else
      open(10, file=trim(adjustl(FILE_S5)),status='old')
      shift_5(0) = 0
      do i = 1, NOS
        read(10,*) shift_5(i)
      end do
      close(10)
    end if
    !
    if(L6==1)then
      call set_identity_ope(shift_6,NOS)
    else
      open(10, file=trim(adjustl(FILE_S6)),status='old')
      shift_6(0) = 0
      do i = 1, NOS
        read(10,*) shift_6(i)
      end do
      close(10)
    end if
    !
    write(*,'(" ### Read FILE_NODmax. ")')
    open(10, file=trim(adjustl(FILE_NODmax)),status='old')
    do i = 1, NOS
      read(10,*) local_NODmax(i)
    end do
    close(10)
    !
    write(*,'(" ### Read FILE_spin. ")')
    open(10, file=trim(adjustl(FILE_spin)),status='old')
    do i = 1, NOS
      read(10,*) local_spin(i)
      if(.not.real8_is_finite(local_spin(i))) &
        error stop 'FILE_spin contains a nonfinite spin'
      ! Bound the value before doubling or converting it to default INTEGER.
      if(local_spin(i)<0.0d0.or.local_spin(i)>dble(huge(0))/2.0d0) &
        error stop 'Spin must be nonnegative with 2*S within the default INTEGER range'
      tmp = int(2.0d0*local_spin(i))
      if(dble(tmp)/=2.0d0*local_spin(i)) &
        error stop 'Spin must be an integer or half-integer'
      if(local_NODmax(i)<0.or.local_NODmax(i)>tmp) &
        error stop 'Local occupation cap must satisfy 0 <= cap <= 2*S'
    end do
    close(10)
    !
    if(cal_cf==1)then
      write(*,'(" ### Read FILE_two_cf. ")')
      allocate(p_two_cf(2,NO_two_cf), swap_list(NO_two_cf))
      open(10, file=trim(adjustl(FILE_two_cf)),status='old')
      do i = 1, NO_two_cf
        read(10,*) p_two_cf(1,i), p_two_cf(2,i)
        if(any(p_two_cf(:,i)<1).or.any(p_two_cf(:,i)>NOS).or.p_two_cf(1,i)==p_two_cf(2,i)) &
          error stop 'Correlation indices must be distinct sites in 1..NOS'
        if( p_two_cf(1,i) > p_two_cf(2,i) )then
          tmp = p_two_cf(1,i)
          p_two_cf(1,i) = p_two_cf(2,i)
          p_two_cf(2,i) = tmp
          swap_list(i) = .true.
        else
          swap_list(i) = .false.
        end if
      end do
      close(10)
    end if
    !
    return
  end subroutine read_ip

  subroutine set_identity_ope(shift_n,NOS)
    implicit none 
    integer :: NOS
    integer :: shift_n(0:NOS)
    integer :: i
    do i = 0, NOS
      shift_n(i) = i
    end do
    return
  end subroutine set_identity_ope

  subroutine write_wf(dim,i_vec_min,i_vec_max, NOE, ene) 
    integer,intent(in)::dim
    integer,intent(in)::i_vec_min,i_vec_max
    integer,intent(in)::NOE
    real(8),intent(in)::ene(1:NOE)
    character(200)::file_name
    character(500)::file_name1
    integer::i_vec,i
    if(wr_wf /= 1) return

    do i_vec=i_vec_min, i_vec_max
      write(file_name,'(A,A,I0,A)')trim(adjustl(FILE_wf)),'wf',i_vec,'.bin'
      write(*,*) " Write wave functions (",i_vec,"-th state) "
      print *, file_name
      open(10, file=trim(adjustl(file_name)), status='replace', action='write', form='unformatted')
      write(10) psi(1:dim,i_vec)
      close(10)
    end do

    write(file_name1,'(A300,A10)')trim(adjustl(FILE_wf)),'energy.dat'
    open(10, file=trim(adjustl(file_name1)), form='formatted')
    do i=1, NOE
      write(10,'(2ES24.15)') dble(i), ene(i) 
    end do
    close(10)

    write(*,'("************************************************************************************")')
    return
  end subroutine write_wf

  subroutine read_wf(dim,i_vec_min, i_vec_max,NOE) 
    integer,intent(in)::dim
    integer,intent(in)::i_vec_min, i_vec_max
    integer,intent(in)::NOE
    character(200)::file_name
    character(500)::file_name1
    integer::i_vec, i
    real(8)::dummy
    if(re_wf /= 1) return

    allocate(psi(1:dim,i_vec_min:i_vec_max))
    do i_vec=i_vec_min, i_vec_max
      write(file_name,'(A,A,I0,A)')trim(adjustl(FILE_wf)),'wf',i_vec,'.bin'
      write(*,*) " read wave functions (",i_vec,"-th state) "
      open(10, file=trim(adjustl(file_name)), status='old', action='read', form='unformatted')
      read(10) psi(1:dim,i_vec)
      close(10)
      if(.not.all(real8_is_finite(real(psi(:,i_vec),8))).or. &
        .not.all(real8_is_finite(aimag(psi(:,i_vec))))) &
        error stop 'Saved wavefunction contains a nonfinite coefficient'
      write(*,*) "first wavefunction coefficients"
      write(*,*) psi(1:min(10,dim),i_vec)
      write(*,*) " finish reading wave function!"
    end do

    write(file_name1,'(A300,A10)')trim(adjustl(FILE_wf)),'energy.dat'
    allocate(ene(1:NOE))
    ene(1:NOE)=0.0d0
    open(10, file=trim(adjustl(file_name1)),status='old')
    do i=1, NOE
      read(10,*) dummy, ene(i) 
      if(.not.real8_is_finite(ene(i))) &
        error stop 'Saved energy.dat contains a nonfinite energy'
      write(*,*) dummy, ene(i), "reading!"
    end do
    close(10)
    write(*,*) "finish reading energy!"

    write(*,'("************************************************************************************")')
    return
  end subroutine read_wf

  subroutine write_index_ene(ene,NOE) 
    integer, intent(in) :: NOE
    real(8), intent(in) :: ene(NOE)
    integer :: i
    open(111,file=trim(adjustl(OUTDIR))//'eigenvalues.dat',position='append')
    do i=1, NOE
      write(111,'(i10, es23.15)') i, ene(i)
    end do
    close(111)
  end subroutine write_index_ene

end module input_param
