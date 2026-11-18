program main
  use iso_fortran_env, only: int64
  !$ use omp_lib
  use state_lists, only: get_cum_ltd_rep_combination,cum_ltd_rep_combination,allocate_lists_omp_SQ
  use input_param, only: read_ip, NOS, NODmax, NODmin, local_NODmax, THS, ALG, re_wf, ene, psi, write_wf, &
    write_index_ene, read_wf, re_wf, wr_wf, cal_lm, cal_cf
  use eigen_solver, only: my_trlanczos_routines_complex, i_vec_min, i_vec_max, NOE, &
    read_input_TRLan, Full_diag_routines
  use ham2vec, only: make_wk_loc_and_ele
  use lanczos, only: read_lanczos_para, lanczos_routines
  use get_expectation_values, only: get_lm, get_cf
  implicit none
  ! Optional elapsed-wall-time telemetry; the default numerical path is unchanged.
  integer(int64) :: clock_origin, clock_previous, clock_rate, clock_now
  integer :: clock_status
  character(len=32) :: clock_switch
  logical :: clock_enabled
  !$ real(8) :: st, en
  !$ st = omp_get_wtime()
  clock_switch=''
  call get_environment_variable('ED2_PHASE_TIMING', clock_switch, status=clock_status)
  clock_enabled = clock_status==0 .and. trim(adjustl(clock_switch))=='1'
  if(clock_enabled)then
    call system_clock(count=clock_origin, count_rate=clock_rate)
    clock_previous=clock_origin
    if(clock_rate<=0_int64)error stop 'Invalid system_clock rate'
  end if
  write(*,'("******************************* Start QS^3 ED2 *************************************")')
  call read_ip
  call phase_mark('input')
  !
  write(*,'(" ### Store cum_ltd_rep_combination numbers. ")')
  call get_cum_ltd_rep_combination(NOS,NODmax,NODmin,local_NODmax)
  !
  write(*,'(" ### Set THS. ")')
  if(NODmax>0)then
    THS = cum_ltd_rep_combination(NOS,NODmax)
  else
    THS = 1
  end if
  write(*,*) "  THS   = ", THS
  !
  if(THS < 0) stop "THS is bigger than the upper limit of 8-bite integer!"
  call phase_mark('combinatorics')
  !
  write(*,'(" ### Allocate and Set arrays for state_list. ")')
  call allocate_lists_omp_SQ 
  write(*,*) "  THS(k)   = ", THS
  call phase_mark('basis')
  !
  If(ALG.eq.1)then !Conventional Lanczos
    !
    write(*,'(" ### Make wk_loc and wk_ele array. ")')
    call make_wk_loc_and_ele(THS)
    call phase_mark('operator_cache')
    !
    if(re_wf.ne.1)then
      write(*,'(" ### Start the Lanczos method. ")')
      call lanczos_routines(ene,THS,psi)
      call phase_mark('lanczos')
      !
      if(wr_wf.eq.1)then
        write(*,'(" ### Write wavevectors. ")')
        call write_wf(THS,1,1,1,ene)
      end if
      !
      call write_index_ene(ene,1)
      !
    else
      write(*,'(" ### Read lanczos_para. ")')
      call read_lanczos_para
      !
      write(*,'(" ### Read wavefunctions. ")')
      call read_wf(THS,1,1,1)
    end if
    !
  else if(ALG.eq.2 )then 
    !
    write(*,'(" ### Make wk_loc and wk_ele array. ")')
    call make_wk_loc_and_ele(THS)
    call phase_mark('operator_cache')
    !
    if(re_wf.ne.1)then
      write(*,'(" ### Start the thick-restarted Lanczos method. ")')
      call my_trlanczos_routines_complex(ene,psi,THS)
      call phase_mark('thick_restart_lanczos')
      !
      if(wr_wf.eq.1)then
        write(*,'(" ### Write wavevectors. ")')
        call write_wf(THS,i_vec_min, i_vec_max,NOE,ene)
      end if
      !
      call write_index_ene(ene,NOE)
      !
    else
      write(*,'(" ### Read input_TRLan. ")')
      call read_input_TRLan
      !
      write(*,'(" ### Read wavefunctions. ")')
      call read_wf(THS,i_vec_min,i_vec_max,NOE)
    end if
    !
  else if(ALG.eq.3)then
    write(*,'(" ### Start the full diagonalization. ")')
    call Full_diag_routines(ene,psi,THS)
    call phase_mark('dense_solve_and_operator')
    if(wr_wf.eq.1)then
      write(*,'(" ### Write wavevectors. ")')
      call write_wf(THS,1,THS,THS,ene)
    end if
  end if
  !
  if(cal_lm==1)then
    write(*,'(" ### Get local magnetizations. ")')
    call get_lm(psi,NODmax,NODmin,THS,1,NOS)
  end if
  !
  if(cal_cf==1)then
    write(*,'(" ### Get two-point correlation functions. ")')
    call get_cf(psi,NODmax,NODmin,THS,1)
  end if
  !
  call phase_mark('observables_and_output')
  if(clock_enabled)then
    call system_clock(count=clock_now)
    write(*,'(a,es18.10)') 'ED2_PHASE total wall_s=',real(clock_now-clock_origin,8)/real(clock_rate,8)
  end if
  !$ en = omp_get_wtime()
  !$ write(*,'("Time:",f10.3," [sec]")') en-st
  write(*,'("****************************************************************************************")')
  stop
contains
  subroutine phase_mark(name)
    character(len=*), intent(in) :: name
    integer(int64) :: next_tick
    if(.not.clock_enabled)return
    call system_clock(count=next_tick)
    write(*,'(a,a,a,es18.10)') 'ED2_PHASE ',trim(name),' wall_s=', &
      real(next_tick-clock_previous,8)/real(clock_rate,8)
    clock_previous=next_tick
  end subroutine phase_mark
end program main
