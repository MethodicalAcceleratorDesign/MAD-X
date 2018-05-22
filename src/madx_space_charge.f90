module spch_bbfi
  use name_lenfi
  use bbfi
  implicit none
  public
  logical :: lost_in_turn = .false., is_lost = .false.
  integer, save :: i_turn, N_macro_surv, N_for_I, N_spch, i_spch
  integer, parameter :: N_macro_max=16000
  double precision, save :: Ex_rms, Ey_rms, sigma_p, sigma_z
  double precision, save :: Ix_array(N_macro_max), Iy_array(N_macro_max)
  double precision, save :: dpi_array(N_macro_max), z_part_array(N_macro_max)
  double precision :: alpha, I_div_E_sum_max
  !  parameter(alpha=0.0, I_div_E_sum_max=7.0)
  double precision, save :: betx_bb(bbd_max), bety_bb(bbd_max), &
       alfx_bb(bbd_max), alfy_bb(bbd_max), &
       gamx_bb(bbd_max), gamy_bb(bbd_max), &
       dx_bb(bbd_max),   dy_bb(bbd_max)
  double precision,save :: rat_bb_n_ions=1d0
  double precision, save :: sigma_t=0.d0, mean_t=0.d0  ! calculate and transfer to BB
  character(len=name_len), save :: spch_bb_name(bbd_max)
end module spch_bbfi

module SpaceCharge

  use math_constfi

  implicit none
  !!private

  logical :: bb_sxy_update, virgin_state, emittance_update
  logical :: checkpnt_restart, exit_loss_turn

  logical :: sc_chrom_fix

  double precision, save :: ex_rms0=zero, ey_rms0=zero, sigma_p0=zero, sigma_z0=zero

  private :: table_input
  !private :: ixy_calcs
  !private :: ixy_fitting

contains

  subroutine SC_Init(first, run, dynap, turns, &
       betax_start, betay_start, &
       alfax_start, alfay_start, &
       gamax_start, gamay_start, &
       dx_start,    dpx_start,   &
       dy_start,    dpy_start)

    logical, intent(IN) :: run, dynap
    integer, intent(IN)  :: turns
    logical, intent(INOUT) :: first
    double precision, intent(IN)  :: betax_start, betay_start
    double precision, intent(IN)  :: alfax_start, alfay_start
    double precision, intent(IN)  :: gamax_start, gamay_start
    double precision, intent(IN)  :: dx_start,    dpx_start
    double precision, intent(IN)  :: dy_start,   dpy_start

    integer, external :: get_option
    double precision, external :: get_value

    exit_loss_turn = get_option('exit_loss_turn ') .ne. 0
    bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
    checkpnt_restart = get_value('run ', 'checkpnt_restart ') .ne. zero
    emittance_update = get_option('emittance_update ') .ne. 0
    virgin_state = get_value('run ', 'virgin_state ') .ne. zero

    !! ALSC: I think init should reset to zero these varibles
    ex_rms0 = zero;
    ey_rms0 = zero;
    sigma_p0 = zero;
    sigma_z0 = zero;

    if (run .and. bb_sxy_update) then
       open(90,file='checkpoint_restart.dat',form='unformatted',status='unknown')
    else if (dynap) then
       bb_sxy_update = .false.
       checkpnt_restart = .false.
    endif

    if (bb_sxy_update) then
       if (virgin_state) first=.true.
       call table_input( betax_start, betay_start, &
            alfax_start, alfay_start, &
            gamax_start, gamay_start, &
            dx_start,    dpx_start, &
            dy_start,    dpy_start)
       if (first) call make_bb6d_ixy(turns)
    endif

  end subroutine SC_Init

  subroutine SC_Update(betas, orbit, z, &
       betax_start, betay_start, &
       alfax_start, alfay_start, &
       gamax_start, gamay_start, &
       dx_start,    dpx_start,   &
       dy_start,    dpy_start)

    use spch_bbfi

    double precision :: betas, orbit(6), z(6,N_macro_surv)
    double precision :: betax_start, betay_start
    double precision :: alfax_start, alfay_start
    double precision :: gamax_start, gamay_start
    double precision :: dx_start,    dpx_start
    double precision :: dy_start,   dpy_start

    !frs on 04.06.2016 - fixing
    !a) bug concerning sigma_p
    !b) Filling data in file bb6d_ixy.txt even for "emittance_update = .false.",
    !   obviously without update!
    !c) Fixing checkpnt_restart for "emittance_update = .false." which
    !   worked for ".true." alright.
    ex_rms0=ex_rms
    ey_rms0=ey_rms
    sigma_z0=sigma_z
    sigma_p0=sigma_p
    if (bb_sxy_update .and. is_lost) then
       call ixy_calcs(betas, orbit, z,       &
            betax_start, betay_start, &
            alfax_start, alfay_start, &
            gamax_start, gamay_start, &
            dx_start,    dpx_start, &
            dy_start,    dpy_start)
       call ixy_fitting()
       is_lost = .false.
    endif

    if( .not. emittance_update) then
       ex_rms=ex_rms0
       ey_rms=ey_rms0
       sigma_z=sigma_z0
       sigma_p=sigma_p0
    endif
    !frs on 04.06.2016 - fixing
    !a) bug concerning sigma_p
    !b) Filling data in file bb6d_ixy.txt even for "emittance_update = .false.",
    !   obviously without update!
    !c) Fixing checkpnt_restart for "emittance_update = .false." which
    !   worked for ".true." alright.    
  end subroutine SC_Update

  subroutine table_input( betx_start, bety_start, &
       alfx_start, alfy_start, &
       gamx_start, gamy_start, &
       dx_start,    dpx_start, &
       dy_start,    dpy_start)
    use name_lenfi
    use bbfi
    use twtrrfi
    use time_varfi
    use spch_bbfi
    use math_constfi, only : zero, one
    implicit none
    double precision :: betx_start, bety_start
    double precision :: alfx_start, alfy_start
    double precision :: gamx_start, gamy_start
    double precision :: dx_start,   dpx_start
    double precision :: dy_start,   dpy_start

    integer :: i, ii, j, flag, range(2)
    double precision :: position
    character(len=name_len) :: name
    character(len=20) :: text

    integer :: double_from_table_row, string_from_table_row, advance_to_pos
    double precision, parameter :: cme10=1d-10

    MYFIELD = zero
    PHASE_TROMB = zero
    CAV_VOLT = zero
    TIME_VAR_M_IND = zero
    TIME_VAR_P_IND = zero
    TIME_VAR_C_IND = zero
    TIME_VAR_M_NT  = zero
    TIME_VAR_P_NT  = zero
    TIME_VAR_C_NT  = zero

    call table_range('spch_bb ', '#s/#e ', range)
    print *, 'Range for Table spch_bb : ', range(1), range(2)
    if (range(1).eq.0 .and. range(2).eq.0) print *," Info: Table spch_bb is empty "

    name=" "
    if (range(2).gt.bbd_max) then
       write(text, '(1p,i8)') bbd_max
       call fort_fail('TRRUN: Fatal: ', &
            'overrun of the number of BB elements in table spch_bb =' // text)
    endif
    N_spch = range(2)-range(1)
    if (N_spch.lt.1) &
         call fort_fail('TRRUN: Fatal: ', 'Table: spch_bb holds no BB elements')

    do i = range(1), range(2)
       j = advance_to_pos('spch_bb ', i)
       if (i .eq. range(1)) then
          flag = string_from_table_row('spch_bb ', 'name ', i, name); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 's ',i, position); if (flag.ne.0) goto 98
          if (name(:10).ne."RING$START" .and. position.ne.zero) then
             write(text, '(1p,i8)') i
             call fort_fail('TRRUN: Fatal: ', &
                  'Global TWISS not readable from table spch_bb'// text)
          endif
          flag = double_from_table_row('spch_bb ', 'betx ', i, betx_start); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'bety ', i, bety_start); if (flag.ne.0) goto 98
          if (abs(betx_start).lt.cme10 .or. abs(bety_start).lt.cme10) then
             write(text, '(1p,i8)') i
             call fort_fail('TRRUN: Fatal: ', &
                  'start beta values from TWISS table smaller than '// &
                  '1e-10, location: ' // text)
          endif
          flag = double_from_table_row('spch_bb ', 'alfx ', i, alfx_start); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'alfy ', i, alfy_start); if (flag.ne.0) goto 98
          gamx_start = (one+alfx_start*alfx_start)/betx_start
          gamy_start = (one+alfy_start*alfy_start)/bety_start
          flag = double_from_table_row('spch_bb ', 'dx ', i,  dx_start); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'dpx ', i, dpx_start); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'dy ', i,  dy_start); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'dpy ', i, dpy_start); if (flag.ne.0) goto 98
       else
          ii=i-1
          flag = string_from_table_row('spch_bb ', 'name ', i, spch_bb_name(ii)); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'betx ', i, betx_bb(ii)); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'bety ', i, bety_bb(ii)); if (flag.ne.0) goto 98
          if (abs(betx_bb(ii)).lt.cme10 .or. abs(bety_bb(ii)).lt.cme10) then
             write(text, '(1p,i8)') i
             call fort_fail('TRRUN: Fatal: ', &
                  'BB beta values from TWISS table smaller '// &
                  'than 1e-10 at location:  '//text)
          endif
          flag = double_from_table_row('spch_bb ', 'alfx ', i, alfx_bb(ii)); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'alfy ', i, alfy_bb(ii)); if (flag.ne.0) goto 98
          gamx_bb(ii) = (one+alfx_bb(ii)*alfx_bb(ii))/betx_bb(ii)
          gamy_bb(ii) = (one+alfy_bb(ii)*alfy_bb(ii))/bety_bb(ii)
          flag = double_from_table_row('spch_bb ', 'dx ', i,  dx_bb(ii)); if (flag.ne.0) goto 98
          flag = double_from_table_row('spch_bb ', 'dy ', i,  dy_bb(ii)); if (flag.ne.0) goto 98
       endif
    enddo
    goto 99

98  write(text, '(1p,i8)') i
    call fort_fail('TRRUN: Fatal: ', 'Table: spch_bb corrupted at row =' // text)

99  continue
    call table_range('time_var_mul ', '#s/#e ', range)
    print *, 'Range for Table time_var_mul: ', range(1), range(2)
    if (range(1).eq.0 .and. range(2).eq.0) then
       print*," Info: Table time_var_mul is empty "
       goto 102
    endif
    name=" "
    if (range(2).gt.n_time_var) then
       write(text, '(1p,i8)') n_time_var
       call fort_fail('TRRUN: Fatal: ', 'overrun of time varying mult arrays =' // text)
    endif
    do i=range(1),range(2)
       j = advance_to_pos('time_var_mul ', i)
       flag = double_from_table_row('time_var_mul ', 'number ',i, time_var_m_ind(i)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'turn ',i, time_var_m_nt(i)); if (flag.ne.0) goto 101
       flag = string_from_table_row('time_var_mul ', 'name ', i, time_var_m_ch(i)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn0 ',i, myfield(i,1,0)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks0 ',i, myfield(i,1,1)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn1 ',i, myfield(i,1,2)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks1 ',i, myfield(i,1,3)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn2 ',i, myfield(i,1,4)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks2 ',i, myfield(i,1,5)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn3 ',i, myfield(i,1,6)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks3 ',i, myfield(i,1,7)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn4 ',i, myfield(i,1,8)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks4 ',i, myfield(i,1,9)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn5 ',i, myfield(i,1,10)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks5 ',i, myfield(i,2,0)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn6 ',i, myfield(i,2,1)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks6 ',i, myfield(i,2,2)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn7 ',i, myfield(i,2,3)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks7 ',i, myfield(i,2,4)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn8 ',i, myfield(i,2,5)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks8 ',i, myfield(i,2,6)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn9 ',i, myfield(i,2,7)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks9 ',i, myfield(i,2,8)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dkn10 ',i, myfield(i,2,9)); if (flag.ne.0) goto 101
       flag = double_from_table_row('time_var_mul ', 'dks10 ',i, myfield(i,2,10)); if (flag.ne.0) goto 101
    enddo
    goto 102

101 write(text, '(1p,i8)') i
    call fort_fail('TRRUN: Fatal: ', 'Table: time_var_mul corrupted at row =' // text)

102 continue
    call table_range('time_var_pha ', '#s/#e ', range)
    print *, 'Range for Table time_var_pha: ', range(1), range(2)
    if (range(1).eq.0 .and. range(2).eq.0) then
       print*," Info: Table time_var_pha is empty "
       goto 104
    endif
    name=" "
    if (range(2).gt.n_time_var) then
       write(text, '(1p,i8)') n_time_var
       call fort_fail('TRRUN: Fatal: ', 'overrun of time varying pha arrays =' // text)
    endif
    do i=range(1),range(2)
       j = advance_to_pos('time_var_pha ', i)
       flag = double_from_table_row('time_var_pha ', 'number ',i, time_var_p_ind(i)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'turn ',i, time_var_p_nt(i)); if (flag.ne.0) goto 103
       flag = string_from_table_row('time_var_pha ', 'name ', i, time_var_p_ch(i)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm11 ',i, phase_tromb(i,1)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm12 ',i, phase_tromb(i,2)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm13 ',i, phase_tromb(i,3)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm14 ',i, phase_tromb(i,4)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm15 ',i, phase_tromb(i,5)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm16 ',i, phase_tromb(i,6)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm21 ',i, phase_tromb(i,7)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm22 ',i, phase_tromb(i,8)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm23 ',i, phase_tromb(i,9)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm24 ',i, phase_tromb(i,10)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm25 ',i, phase_tromb(i,11)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm26 ',i, phase_tromb(i,12)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm31 ',i, phase_tromb(i,13)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm32 ',i, phase_tromb(i,14)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm33 ',i, phase_tromb(i,15)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm34 ',i, phase_tromb(i,16)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm35 ',i, phase_tromb(i,17)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm36 ',i, phase_tromb(i,18)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm41 ',i, phase_tromb(i,19)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm42 ',i, phase_tromb(i,20)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm43 ',i, phase_tromb(i,21)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm44 ',i, phase_tromb(i,22)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm45 ',i, phase_tromb(i,23)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm46 ',i, phase_tromb(i,24)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm51 ',i, phase_tromb(i,25)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm52 ',i, phase_tromb(i,26)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm53 ',i, phase_tromb(i,27)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm54 ',i, phase_tromb(i,28)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm55 ',i, phase_tromb(i,29)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm56 ',i, phase_tromb(i,30)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm61 ',i, phase_tromb(i,31)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm62 ',i, phase_tromb(i,32)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm63 ',i, phase_tromb(i,33)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm64 ',i, phase_tromb(i,34)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm65 ',i, phase_tromb(i,35)); if (flag.ne.0) goto 103
       flag = double_from_table_row('time_var_pha ', 'rm66 ',i, phase_tromb(i,36)); if (flag.ne.0) goto 103
    enddo
    goto 104

103 write(text, '(1p,i8)') i
    call fort_fail('TRRUN: Fatal: ', 'Table: time_var_pha corrupted at row =' // text)

104 continue
    call table_range('time_var_cav ', '#s/#e ', range)
    print *, 'Range for Table time_var_cav: ', range(1), range(2)
    if (range(1).eq.0 .and. range(2).eq.0) then
       print*," Info: Table time_var_cav is empty "
       goto 106
    endif
    name=" "
    if (range(2).gt.n_time_var) then
       write(text, '(1p,i8)') n_time_var
       call fort_fail('TRRUN: Fatal: ', 'overrun of time varying cav arrays =' // text)
    endif
    do i=range(1),range(2)
       j = advance_to_pos('time_var_cav ', i)
       flag = double_from_table_row('time_var_cav ', 'number ',i, time_var_c_ind(i)); if (flag.ne.0) goto 105
       flag = double_from_table_row('time_var_cav ', 'turn ',i, time_var_c_nt(i)); if (flag.ne.0) goto 105
       flag = string_from_table_row('time_var_cav ', 'name ', i, time_var_c_ch(i)); if (flag.ne.0) goto 105
       flag = double_from_table_row('time_var_cav ', 'volt ',i, cav_volt(i)); if (flag.ne.0) goto 105
    enddo
    goto 106

105 write(text, '(1p,i8)') i
    call fort_fail('TRRUN: Fatal: ', 'Table: time_var_cav corrupted at row =' // text)

106 continue
  end subroutine table_input

  subroutine ixy_fitting()
    use name_lenfi
    use bbfi
    use spch_bbfi
    use math_constfi, only : zero, one
    implicit none

    integer :: i, iii, jjj
    integer :: i_for_I
    double precision :: Summ_dpi_square, Summ_z_part_square
    double precision :: Summ_x, Summ_y
    double precision :: Ix(N_macro_max), Iy(N_macro_max)
    double precision :: dpi(N_macro_max), z_part(N_macro_max)
    double precision :: Ix_sorted(N_macro_max), Iy_sorted(N_macro_max)
    double precision :: Ix_min,Iy_min, Ix_min_last,Iy_min_last
    ! I_x/Ex_rms+I_y/Ey_rms <= I_div_E_sum_max
    ! limit for particles taken for Ix, Iy evaluations
    double precision :: Ix_i, Iy_i, dpi_i, z_part_i, Ixy_rel_summ, N_for_I_dble
    double precision :: c_sumx, y_sumx, t_sumx, c_sumy, y_sumy, t_sumy, a_sum
    double precision :: J0x, J0y, Sum_Jx, Sum_Jy
    double precision, parameter :: ce10=1d10

    integer, external :: get_option
    double precision, external :: get_value

    !----------------------------------------------------------------------
    J0x=zero; J0y=zero; Sum_Jx=zero; Sum_Jy=zero



    i_for_I=0     ! I-evaluations
    I_div_E_sum_max = get_value('run ', 'i_div_e_sum_max ')
    DO i=1,N_macro_surv
       Ix_i=Ix_array(i)
       Iy_i=Iy_array(i)
       dpi_i=dpi_array(i)
       z_part_i=z_part_array(i)
       ! Lethal bug found by Valery!!!
       !     Ixy_rel_summ=Ix_i/Ex_rms+Iy_i/Ex_rms
       Ixy_rel_summ=Ix_i/Ex_rms+Iy_i/Ey_rms
       if (Ixy_rel_summ .LE. I_div_E_sum_max) then
          i_for_I=i_for_I+1
          Ix(i_for_I)=Ix_i
          Iy(i_for_I)=Iy_i
          dpi(i_for_I)=dpi_i
          z_part(i_for_I)=z_part_i
       endif
    ENDDO
    if (i_for_I.eq.0) then
       call fort_warn('trrun: ','the RMS emittances cannot be calculated: exit from IXY_FITTING');
       return
    endif
    N_for_I=i_for_I

    N_for_I_dble=dble(N_for_I)

    !     rms of the relative momentum
    Summ_dpi_square = zero
    sigma_p_loop: DO iii=1, N_for_I
       Summ_dpi_square = Summ_dpi_square + dpi(iii)**2
    ENDDO sigma_p_loop

    Summ_dpi_square = Summ_dpi_square / N_for_I_dble
    if (Summ_dpi_square .ge. zero) then
       sigma_p=sqrt(Summ_dpi_square)
    else
       call fort_fail('IXY_FITTING: Fatal: ','Summ_dpi_square<0')
    endif

    !     rms of the bunch length
    Summ_z_part_square = zero
    sigma_z_loop: DO iii=1, N_for_I
       Summ_z_part_square = Summ_z_part_square + z_part(iii)**2
    ENDDO sigma_z_loop
    Summ_z_part_square = Summ_z_part_square / N_for_I_dble
    if (Summ_z_part_square .ge. zero) then
       sigma_z = sqrt(Summ_z_part_square)
    else
       call fort_fail('IXY_FITTING: Fatal: ','Summ_z_part_square<0')
    endif

    !     SORTING Ix
    Ix_min_last = zero
    iii_loop: DO iii=1, N_for_I

       Ix_min=ce10
       !$OMP PARALLEL PRIVATE(jjj) SHARED(Ix_min)
       !$OMP DO REDUCTION(MIN:Ix_min)
       jjj_loop: DO jjj=1,N_for_I
          if (Ix(jjj) < Ix_min .AND. Ix(jjj) > Ix_min_last) then
             Ix_min=Ix(jjj)
          endif
       ENDDO jjj_loop
       !$OMP END DO
       !$OMP END PARALLEL
       Ix_sorted(iii) = Ix_min
       Ix_min_last = Ix_min
    ENDDO iii_loop

    Iy_min_last = zero
    iii_loop_y: DO iii=1, N_for_I

       Iy_min=ce10
       !$OMP PARALLEL PRIVATE(jjj) SHARED(Iy_min)
       !$OMP DO REDUCTION(MIN:Iy_min)
       jjj_loop_y: DO jjj=1,N_for_I
          if (Iy(jjj) < Iy_min .AND. Iy(jjj) > Iy_min_last) then
             Iy_min=Iy(jjj)
          endif
       ENDDO jjj_loop_y
       !$OMP END DO
       !$OMP END PARALLEL
       Iy_sorted(iii) = Iy_min
       Iy_min_last = Iy_min
    ENDDO iii_loop_y

    !New normalisation of the emittance calculation to exclude artificial collapses
    !R.Wasef
    J0x = zero
    J0y = zero
    jjj_loop_ray: DO jjj=1,N_for_I
       J0x = J0x + Ix_sorted(jjj)
       J0y = J0y + Iy_sorted(jjj)
    ENDDO jjj_loop_ray
    J0x = J0x / (N_for_I_dble*N_for_I_dble)
    J0y = J0y / (N_for_I_dble*N_for_I_dble)
    J0x = J0x*J0x
    J0y = J0y*J0y

    Sum_Jx = zero
    Sum_Jy = zero
    jjj_loop_ray2: DO jjj=1,N_for_I
       Sum_Jx = Sum_Jx + ( (Ix_sorted(jjj)*Ix_sorted(jjj)) /( J0x+(Ix_sorted(jjj)*Ix_sorted(jjj)) ) )
       Sum_Jy = Sum_Jy + ( (Iy_sorted(jjj)*Iy_sorted(jjj)) /( J0y+(Iy_sorted(jjj)*Iy_sorted(jjj)) ) )
    ENDDO jjj_loop_ray2

    !     Summ of step-function for Ex/Ey evaluation
    !     Kahan summation algorithm
    Summ_x = zero
    Summ_y = zero
    c_sumx = zero
    c_sumy = zero
    alpha = get_value('run ', 'alpha ')
!!!!!$OMP PARALLEL PRIVATE(iii,a_sum,c_sumx,c_sumy,y_sumx,y_sumy,t_sumx,t_sumy)
!!!!!$OMP DO REDUCTION(+:Summ_x,Summ_y)
    Summ_loop: DO iii=1,N_for_I    !!! the first particle is shifted
       a_sum = Log(one-(alpha+dble(iii-1))/N_for_I_dble)
       y_sumx = (Ix_sorted(iii)/(Ix_sorted(iii)*Ix_sorted(iii)+J0x))*a_sum/Sum_Jx-c_sumx
       y_sumy = (Iy_sorted(iii)/(Iy_sorted(iii)*Iy_sorted(iii)+J0y))*a_sum/Sum_Jy-c_sumy
       t_sumx = Summ_x + y_sumx
       t_sumy = Summ_y + y_sumy
       c_sumx = (t_sumx-Summ_x) - y_sumx
       c_sumy = (t_sumy-Summ_y) - y_sumy
       Summ_x = t_sumx
       Summ_y = t_sumy
    ENDDO Summ_loop
!!!!!$OMP END DO
!!!!!$OMP END PARALLEL
    Ey_rms = -one/Summ_y
    Ex_rms = -one/Summ_x

    return
  END subroutine ixy_fitting

  subroutine ixy_calcs(betas, orbit, z,          &
       betax_start, betay_start, &
       alfax_start, alfay_start, &
       gamax_start, gamay_start, &
       dx_start,    dpx_start,   &
       dy_start,    dpy_start)
    use name_lenfi
    use bbfi
    use spch_bbfi
    use math_constfi, only : two

    implicit none
    double precision :: betas, orbit(6), z(6,N_macro_surv)
    double precision :: betax_start, betay_start
    double precision :: alfax_start, alfay_start
    double precision :: gamax_start, gamay_start
    double precision :: dx_start,    dpx_start
    double precision :: dy_start,   dpy_start
    integer :: get_option

    integer :: i
    double precision :: dpi, xi, pxi, yi, pyi
    !-----------------------------------------------------------------------

    !      DPI=pt_000KK/beta_part_ini; z_part=t_000KK*beta_part_ini;
    !      XI=x_000KK-dx_start*DPI; PXI=px_000KK-dpx_start*DPI;
    !      YI=y_000KK-dy_start*DPI; PYI=py_000KK-dpy_start*DPI;
    !      Ix=(gamax_start*XI*XI+2*alfax_start*XI*PXI+betax_start*PXI*PXI)/2;
    !      Iy=(gamay_start*YI*YI+2*alfay_start*YI*PYI+betay_start*PYI*PYI)/2;

    sc_chrom_fix = get_option('sc_chrom_fix ') .ne. 0 !! frs add-on
    do i=1,N_macro_surv
       ! Exact formulation might be too computational time costly
       !        DPI=(sqrt((one+z(6,i)*betas)**2-gammas**(-2)))/betas-one
       if(sc_chrom_fix .eqv. .true.) then !! frs add-on
          DPI = z(6,i) - orbit(6)
       else
          DPI = (z(6,i) - orbit(6)) / betas !! orignal
       endif
       XI  =  z(1,i) - orbit(1) - dx_start  * DPI
       PXI =  z(2,i) - orbit(2) - dpx_start * DPI
       YI  =  z(3,i) - orbit(3) - dy_start  * DPI
       PYI =  z(4,i) - orbit(4) - dpy_start * DPI

       Ix_array(i) = (gamax_start*XI*XI + two*alfax_start*XI*PXI + betax_start*PXI*PXI) / two
       Iy_array(i) = (gamay_start*YI*YI + two*alfay_start*YI*PYI + betay_start*PYI*PYI) / two
       dpi_array(i) = DPI
       if(sc_chrom_fix .eqv. .true.) then !! frs add-on
          z_part_array(i)=z(5,i) - orbit(5)
       else
          z_part_array(i) = (z(5,i)-orbit(5))*betas !! orignal
       endif
    enddo

  end subroutine ixy_calcs

end module SpaceCharge
