module bbfi
  implicit none
  public
  integer, parameter :: bbd_max=100000
  integer :: bbd_loc(bbd_max)=0, bbd_cnt=0, bbd_flag=0, bbd_pos=0
  double precision :: bb_kick(2,bbd_max)=0.d0
  double precision, parameter :: explim=150.0d0   ! if x > explim, exp(-x) is outside machine limits.
end module bbfi

module spch_bbfi
  use name_lenfi
  use bbfi
  implicit none
  public
  integer N_ini,i_turn, N_macro_surv, N_for_I, N_macro_max, N_spch,i_spch,unit_chpt
  integer, parameter :: lu_max=1000
  parameter(N_macro_max=100000)
  double precision Ex_rms, Ey_rms, sigma_p, sigma_z
  double precision Ix_array(N_macro_max), Iy_array(N_macro_max),    &
       dpi_array(N_macro_max),                          &
       z_part_array(N_macro_max)
  double precision alpha, I_div_E_sum_max
  !  parameter(alpha=0.0, I_div_E_sum_max=7.0)
  double precision betx_bb(bbd_max), bety_bb(bbd_max), alfx_bb(bbd_max), &
       alfy_bb(bbd_max), gamx_bb(bbd_max), gamy_bb(bbd_max), dx_bb(bbd_max), &
       dy_bb(bbd_max), scsigx(bbd_max), scsigy(bbd_max), scsigz(bbd_max)
  double precision sc_intstr,sc_charge(bbd_max)
  double precision :: sc_map(bbd_max,6,6)=0.d0
  double precision sc_sect_map(bbd_max,6,6), trans_sc_sect_map(bbd_max,6,6)
  double precision sect_map(bbd_max,6,6), trans_sect_map(bbd_max,6,6)
  double precision rat_bb_n_ions,R_part
  double precision ::  sigma_t=0.d0, mean_t=0.d0  ! calculate and transfer to BB
  character*(name_len) spch_bb_name(bbd_max)
  save N_ini,i_turn,N_macro_surv,N_for_I,N_spch,i_spch,                                &
       Ex_rms,Ey_rms,sigma_p,sigma_z,                                           &
       Ix_array,Iy_array,dpi_array, z_part_array,                               &
       betx_bb,bety_bb,alfx_bb,alfy_bb,gamx_bb,gamy_bb,dx_bb,dy_bb,             &
       rat_bb_n_ions,sigma_t, mean_t,spch_bb_name,R_part
  data rat_bb_n_ions / 1d0 /
  integer I_NUMBER
  logical I_OPEN
end module spch_bbfi

module SpaceCharge

  use trackfi, only : max_part
  use math_constfi

  implicit none
  !!private

  logical :: bb_sxy_update, virgin_state, emittance_update
  logical :: checkpnt_restart, exit_loss_turn

  logical :: sc_chrom_fix

  double precision sgpr(6,6),npart,sc3d_emit(5)
  double precision :: ex_rms0 = zero, ey_rms0 = zero, sigma_p0 = zero, sigma_z0 = zero
  double precision :: N_ions_in_beam, Npart_gain, N_ions_ini, n_ions_macro, N_ions_for_bb
  double precision :: sigma_z_ini, z_factor, t_rms, pt_rms, z_keep(6,max_part)

  double precision, save :: betx_start=1d0, bety_start=1d0
  double precision, save :: alfx_start=0d0, alfy_start=0d0
  double precision, save :: gamx_start=0d0, gamy_start=0d0
  double precision, save :: dx_start=0d0,   dpx_start=0d0
  double precision, save :: dy_start=0d0,   dpy_start=0d0

  !VVK 20100321 -------------------------------------------------
  integer :: i_part                     ! local counter
  double precision  :: Summ_t_mean      ! local for mean value
  double precision  :: Summ_t_square    ! local for rms value
  !-------------------------------------------------------------------

  private :: table_input
  private :: ixy_calcs
  private :: ixy_fitting

contains

  subroutine SC_Init(first, run, dynap, turns)

    use trackfi
    use spch_bbfi
    use SCdat

    logical, intent(IN) :: run, dynap
    integer, intent(IN)  :: turns
    logical, intent(INOUT) :: first
    integer code, j, restart_sequ
    integer advance_node
    double precision node_value
    integer, external :: get_option
    double precision, external :: get_value
    integer, external :: get_file_unit

    !---- get options for space charge variables
    exit_loss_turn = get_option('exit_loss_turn ') .ne. 0
    bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
    sc_3d_kick = get_option('sc_3d_kick ') .ne. 0
    sc_3d_beamsize = get_option('sc_3d_beamsize ') .ne. 0
    if(sc_3d_beamsize) then
       i_spch = 0 !a special spch-update counter
       j = restart_sequ()
11     continue
       code = node_value('mad8_type ')
       if(code.eq.22) then
          i_spch = i_spch+1
          sc_intstr=2d0*0.5d0**1.5d0*arad * get_value('probe ', 'npart ') / get_value('probe ','gamma ')
          sc_charge(i_spch)=node_value('charge ')
          do j=1,6
             sc_map(i_spch,j,j)=1d0
          enddo
       endif
       if (advance_node().ne.0)  then
          j=j+1
          go to 11
       endif
       i_spch = 0 !a special spch-update counter
    endif
    sc_3d_periodic = get_option('sc_3d_periodic ') .ne. 0
    !  sc_3d_damp = get_option('sc_3d_damp ') .ne. 0 ! replaced by periodic approach
    !  sc_3d_damp_amp = get_value('run ','sc_3d_damp_amp ') ! replaced by periodic approach
    checkpnt_restart = get_value('run ', 'checkpnt_restart ') .ne. 0d0
    emittance_update = get_option('emittance_update ') .ne. 0
    if(sc_3d_kick) emittance_update=.FALSE.
    virgin_state = get_value('run ', 'virgin_state ') .ne. 0d0
    if(run) then
       ! 2015-Feb-23  16:20:19  ghislain: open file only when necessary
       if (bb_sxy_update) then
          inquire (file='checkpoint_restart.dat', OPENED=I_OPEN, NUMBER=I_NUMBER)
          if (I_OPEN) close(I_NUMBER)
          unit_chpt=get_file_unit(lu_max)
          if(checkpnt_restart) then
             open(unit_chpt,file='checkpoint_restart.dat',form='unformat&
                  &ted',status='old')
          else
             open(unit_chpt,file='checkpoint_restart.dat',form='unformatted')
          endif
       endif
    else
       bb_sxy_update = .false.
       checkpnt_restart = .false.
    endif
    if(bb_sxy_update) then
       if(virgin_state) first=.true.

       call table_input(                                 &
            betx_start, bety_start,                      &
            alfx_start, alfy_start,                      &
            gamx_start, gamy_start,                      &
            dx_start,    dpx_start,                      &
            dy_start,    dpy_start)

       if(sc_3d_kick.and.first) call mymap()

       if(first) call make_bb6d_ixy(turns)
    endif

  end subroutine SC_Init

  subroutine SC_Update(orbit, z)

    use spch_bbfi
    use trackfi
    use time_varfi

    double precision, intent(IN) :: orbit(6), z(6,N_macro_surv)

    !frs on 04.06.2016 - fixing
    !a) bug concerning sigma_p
    !b) Filling data in file bb6d_ixy.txt even for "emittance_update = !.false.",
    !   obviously without update!
    !c) Fixing checkpnt_restart for "emittance_update = !.false." which
    !   worked for ".true." alright.
    ex_rms0=ex_rms
    ey_rms0=ey_rms
    sigma_z0=sigma_z
    sigma_p0=sigma_p
    if (bb_sxy_update .and. is_lost) then
       call ixy_calcs(betas, orbit, z,                  &
            betx_start, bety_start,                      &
            alfx_start, alfy_start,                      &
            gamx_start, gamy_start,                      &
            dx_start,    dpx_start,                      &
            dy_start,    dpy_start)
       call ixy_fitting()
       is_lost = .false.
    endif
    if( .not. emittance_update) then
       ex_rms=ex_rms0
       ey_rms=ey_rms0
       rat_bb_n_ions=1d0
    endif
    sigma_z=sigma_z0
    sigma_p=sigma_p0
    !frs on 04.06.2016 - fixing
    !a) bug concerning sigma_p
    !b) Filling data in file bb6d_ixy.txt even for "emittance_update = !.false.",
    !   obviously without update!
    !c) Fixing checkpnt_restart for "emittance_update = !.false." which
    !   worked for ".true." alright.

  end subroutine SC_Update

  subroutine BB_Init(first)

    use trackfi
    use spch_bbfi
    use time_varfi

    logical, intent(INOUT) :: first

    double precision, external :: get_value

    if (bb_sxy_update) then
       trrun_nt = 0

       if (first) then
          time_var_m_cnt = 0 ; time_var_p_cnt = 0 ; time_var_c_cnt = 0
          ! <<N_macro_part_ini = N_macro_surv + N_macro_lost>>
          N_ions_in_beam = get_value('probe ', 'npart ') !BEAM->NPART
          if (N_ions_in_beam .lt. zero) call fort_fail('TRRUN: ','N_ions_in_beam .lt. zero')
          Npart_gain = get_value('run ', 'n_part_gain ')
          N_ions_ini = Npart_gain * N_ions_in_beam
          N_macro_surv = jmax    ! = number of START lines submitted
          n_ions_macro = N_ions_ini/N_macro_surv

          N_for_I = N_macro_surv ! at start (to be redefined in Ixy)
          if (N_macro_surv .gt. N_macro_max) &
               call fort_fail('TRRUN: ', 'Number N_macro_surv exceeds N_macro_max (array size)')

          ! 2015-Jul-03  18:07:00  ghislain: BUG or voluntary ?
          if (N_macro_surv .gt. N_macro_surv) &
               call fort_fail('TRRUN: ', 'Number START-lines exceeds the initial number of macroparticles N_macro_surv')

          t_rms = get_value('run ', 'sigma_z ') * beti
          pt_rms = get_value('run ', 'deltap_rms ')
          pt_rms = (sqrt((betas * (pt_rms + one))**2 + one/gammas**2) - one) * beti
          sigma_z_ini = t_rms !betas: BEAM->BETA
          sigma_z = sigma_z_ini !at start (to be redefined in Ixy)
          sigma_p = pt_rms       !default
          z_factor = one !at start sigma_z_ini/sigma_z
          Ex_rms = get_value('probe ', 'ex ') !BEAM->Ex
          Ey_rms = get_value('probe ', 'ey ') !BEAM->Ey
          if (checkpnt_restart.and.emittance_update) then
             Ex_rms = Ex_rms0
             Ey_rms = Ey_rms0
          endif

          ! write(8,'(4(g16.9,1x))') Ex_rms, Ey_rms,sigma_z,sigma_p

          first = .false.
       endif
    endif
  end subroutine BB_Init

  subroutine BB_Update(turn, orbit0, z)

    use spch_bbfi
    use trackfi
    use time_varfi
    use SCdat

    integer, intent(IN) :: turn

    double precision, intent(IN) :: orbit0(6), z(6,N_macro_surv)

    if(bb_sxy_update) then
       trrun_nt=turn
       time_var_m_lnt=0
       time_var_p_lnt=0
       time_var_c_lnt=0

       N_macro_surv=jmax
       i_spch = 0 !a special spch-update counter

       ex_rms0=ex_rms
       ey_rms0=ey_rms
       sigma_z0=sigma_z
       sigma_p0=sigma_p
       !sigma_p0 = sigma_p !CM, 3/11/14
       !fill, table=Ixy_unsorted; column=i_macro_part, Ix, Iy, dpi, z_part;
       !new on 3/31/14:
       !frs on 04.06.2016 - fixing
       !a) bug concerning sigma_p
       !b) Filling data in file bb6d_ixy.txt even for "emittance_update = !.false.",
       !   obviously without update!
       !c) Fixing checkpnt_restart for "emittance_update = !.false." which
       !   worked for ".true." alright.
       !        if (emittance_update) then

       if(sc_3d_kick) then
          call SigmaMatrixFit6Dv7(z,jmax,sgpr,sc3d_emit)
          if(sc_3d_beamsize) then
             call SigmaTransport2(sgpr)
          else
             call SigmaTransport(sgpr)
          endif
          call double_to_table_curr('bb6d_ixy ', 'turn ', dble(tot_turn+turn))
          call double_to_table_curr('bb6d_ixy ', 'n_macro_surv ', dble(N_ini))
          call double_to_table_curr('bb6d_ixy ', 'n_for_i ', dble(jmax))
          call double_to_table_curr('bb6d_ixy ', 'ex_rms ', sc3d_emit(1))
          call double_to_table_curr('bb6d_ixy ', 'ey_rms ', sc3d_emit(2))
          call double_to_table_curr('bb6d_ixy ', 'ez_rms ', sc3d_emit(3))
          call double_to_table_curr('bb6d_ixy ', 'sigma_ct ', sc3d_emit(4))
          call double_to_table_curr('bb6d_ixy ', 'sigma_pt ', sc3d_emit(5))
       else
          call ixy_calcs(betas, orbit0, z,                  &
               betx_start, bety_start,                      &
               alfx_start, alfy_start,                      &
               gamx_start, gamy_start,                      &
               dx_start,    dpx_start,                      &
               dy_start,    dpy_start)
          call ixy_fitting()
          call double_to_table_curr('bb6d_ixy ', 'turn ', dble(tot_turn+turn))
          call double_to_table_curr('bb6d_ixy ', 'n_macro_surv ', dble(n_macro_surv))
          call double_to_table_curr('bb6d_ixy ', 'n_for_i ', dble(n_for_i))
          call double_to_table_curr('bb6d_ixy ', 'ex_rms ', ex_rms)
          call double_to_table_curr('bb6d_ixy ', 'ey_rms ', ey_rms)
          call double_to_table_curr('bb6d_ixy ', 'ez_rms ', 0d0)
          call double_to_table_curr('bb6d_ixy ', 'sigma_ct ', sigma_z0)
          call double_to_table_curr('bb6d_ixy ', 'sigma_pt ', sigma_p0)
       endif
       call augment_count('bb6d_ixy ')

       if(sigma_p0.eq.0d0) sigma_p0=sigma_p
       !frs on 04.06.2016 - fixing
       !a) bug concerning sigma_p
       !b) Filling data in file bb6d_ixy.txt even for "emittance_update = !.false.",
       !   obviously without update!
       !c) Fixing checkpnt_restart for "emittance_update = !.false." which
       !   worked for ".true." alright.
       !new on 3/31/14:
       !        endif

       sigma_z=sigma_z0
       sigma_p=sigma_p0
       !frs on 04.06.2016 - fixing
       !a) bug concerning sigma_p
       !b) Filling data in file bb6d_ixy.txt even for "emittance_update = !.false.",
       !   obviously without update!
       !c) Fixing checkpnt_restart for "emittance_update = !.false." which
       !   worked for ".true." alright.
       !           sigma_p=sigma_p0
       if((sigma_z.GT.0d0).AND.(sigma_z_ini.GT.0d0)) then
          z_factor=sigma_z_ini/sigma_z
       else
          z_factor=1D0
       endif

       N_ions_for_bb=n_ions_macro*N_for_I*z_factor
       if(N_ions_in_beam.le.0d0) then
          rat_bb_n_ions=0d0
       else
          rat_bb_n_ions=N_ions_for_bb/N_ions_in_beam
       endif

       if(.not.emittance_update) then
          ex_rms=ex_rms0
          ey_rms=ey_rms0
          rat_bb_n_ions=1d0
       endif

       if(idnint(time_var_m_nt(time_var_m_cnt+1)).eq.tot_turn+turn) then
          time_var_m=.true.
       else
          time_var_m=.false.
       endif
       if(idnint(time_var_p_nt(time_var_p_cnt+1)).eq.tot_turn+turn) then
          time_var_p=.true.
       else
          time_var_p=.false.
       endif
       if(idnint(time_var_c_nt(time_var_c_cnt+1)).eq.tot_turn+turn) then
          time_var_c=.true.
       else
          time_var_c=.false.
       endif
    endif

    !ALSC     nlm = 0
    !ALSC     sum=0d0

    !frs on 04.06.2016 - fixing
    !a) bug concerning sigma_p
    !b) Filling data in file bb6d_ixy.txt even for "emittance_update = !.false.",
    !   obviously without update!
    !c) Fixing checkpnt_restart for "emittance_update = !.false." which
    !   worked for ".true." alright.
    if(bb_sxy_update) then
       if(emittance_update.or.(.not.emittance_update.and.mean_t.eq.0d0.and.sigma_t.eq.0d0)) then

          !     if((sigma_t.eq.0d0.or.emittance_update).and.bb_sxy_update) then
          !VVK 20100321 -------- Find RMS-value of t ----------------------
          ! if we do 1-turn tracking, orbit0(5)=0 always
          Summ_t_mean=0d0
          do i_part=1, jmax
             if (abs(z(5,i_part)) .GE. 0d0) then
                Summ_t_mean=Summ_t_mean+z(5,i_part)
             else
                Print *, 'NaN z(5,i) ? :', i_part, z(5,i_part)
             endif
          enddo
          mean_t=Summ_t_mean/dble(jmax)

          Summ_t_square=0d0
          do i_part=1, jmax
             if (abs(z(5,i_part)) .GE. 0d0) &
                  Summ_t_square=Summ_t_square+(z(5,i_part)-mean_t)**2
          enddo
          sigma_t=sqrt(Summ_t_square/dble(jmax))
          if(abs(sigma_t).eq.0d0) then
             sigma_t=t_max/2d0
             call fort_warn('TTRUN Frozen SC: sigma_t = zero: ','sigma_t set to L/track_harmon/betas/2')
          endif
       endif
       sigma_t=sigma_z
       !-----------------------------------------------------------------

    endif

  end subroutine BB_Update

  subroutine BB_Update2(turn, orbit0, z, part_id, last_turn)

    use spch_bbfi
    use trackfi
    use time_varfi
    use SCdat

    integer, intent(IN) :: turn
    integer :: i, j, part_id(*), last_turn(*)

    double precision, intent(IN) :: orbit0(6), z(6,N_macro_surv)

    if(bb_sxy_update) then
       tot_turn=tot_turn+turn
       !$OMP PARALLEL PRIVATE(i,j)
       !$OMP DO
       do i = 1, jmax
          part_id_keep(i)=part_id(i)
          last_turn_keep(i)=last_turn(i)
          do j=1,6
             z_keep(j,i)=z(j,i)
          enddo
       enddo
       !$OMP END DO
       !$OMP END PARALLEL
    endif

  end subroutine BB_Update2

  subroutine BB_Write(turn, orbit0, z)

    use spch_bbfi
    use trackfi
    use time_varfi
    use SCdat

    integer, intent(IN) :: turn
    integer :: i, j

    double precision, intent(IN) :: orbit0(6), z(6,N_macro_surv)

    if (bb_sxy_update) then
       rewind unit_chpt
       write(unit_chpt) jmax
       write(unit_chpt) Ex_rms
       write(unit_chpt) Ey_rms
       do i = 1, jmax
          do j=1,6
             write(unit_chpt) z(j,i)
          enddo
       enddo
       write(unit_chpt) sigma_t
       write(unit_chpt) mean_t
       write(unit_chpt) N_ini
    endif
  end subroutine BB_Write

  subroutine table_input(betx_start, bety_start, &
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
    double precision, intent(IN) :: betas, orbit(6), z(6,N_macro_surv)
    double precision, intent(IN) :: betax_start, betay_start
    double precision, intent(IN) :: alfax_start, alfay_start
    double precision, intent(IN) :: gamax_start, gamay_start
    double precision, intent(IN) :: dx_start,   dpx_start
    double precision, intent(IN) :: dy_start,   dpy_start
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

  subroutine SigmaMatrixFit6Dv7(z,jmax,sgnw,sc3d_emit)
    use SCdat
    use spch_bbfi
    implicit none
    ! Yuri Alexahin Oct 2017 Mathematica version
    !   SigmaMatrixFit6Dv7::usage = "SigmaMatrixFit6D[z,jmax,sgnw,sc3d_emit] makes
    !   Gaussian fit for jmax vectors z
    !   fcore::usage = "fraction of particles in the CORE provided by
    !   the fit to the calling program "(* fcore=1  for now *)
    ! Frank Schmidt Oct 2017 Fortran90 version
    ! Copyright Fermilab & CERN
    integer ip,k,ii,i,i6,j,jmax,nf1
    integer :: itr=0
    double precision s3,exf,etaf,etapr,sgnwc!frs 18.04.2019,trpr
    double precision :: aver(idim),dz(idim),apr(idim),s1(idim),avr(idim),t6(idim),&
         sgin(idim,idim),t66(idim,idim),sgpr(idim,idim),s2(idim,idim),sgnw(idim,idim),&
         sgnwa(idim,idim),reeig(idim),aieig(idim),am(idim,idim),sc3d_emit(5)
    double precision :: z(idim,*),dummy
    !frs 18.04.2019  double precision, external :: mytrace
    character*240 ch
    integer, external :: get_file_unit
    character(20) text
    !**********************************************************************************
    aver= sum(z(1:idim,1:jmax),dim=2)/dble(jmax)
    sgin=0d0
    do k=1,jmax
       t6=z(1:idim,k)-aver
       sgin=sgin+SPREAD(t6,1,idim)*TRANSPOSE(SPREAD(t6,1,idim))
    end do
    sgin=sgin/dble(jmax);!sgnw=sgin;goto 200;
    apr = aver; sgpr = sgin; etapr=1d0;
    do ii=1,isigmatfit
       call symsol(sgin,idim,eflag,work_1,work_2,work_3)
       !frs 18.04.2019     trpr = mytrace(sgpr)
       s1 = 0d0; s2 = 0d0; s3 = 0d0;
       do ip=1,jmax
          dz=z(1:idim,ip)-apr
          !        exf = Exp(-dot_product(dz,(MATMUL(sgin,dz)))/2d0); !old
          exf = Exp(-alfa*dot_product(dz,(MATMUL(sgin,dz)))); !frs 13.06.18
          s1 = s1 + z(1:idim,ip)*exf;
          s2 = s2 + SPREAD(dz,1,idim)*TRANSPOSE(SPREAD(dz,1,idim))*exf
          s3 = s3 + exf
       enddo
       !     avr = s1/s3; sgnw = s2/(s3 - etapr*dble(jmax)/16d0); !old
       avr = s1/s3; sgnw = s2/(s3 - etapr*dble(jmax)*efactor);etaf=1d0; !frs 13.06.18
       !     etaf = Min(1d0, 8d0*s3/dble(jmax));
       sgnwc=0d0
       !frs 18.04.2019
       do i6=1,idim
          if(sgpr(i6,i6).eq.0d0) then
             write(text, '(1p,i8)') i6
             call fort_fail('TRRUN: ','Fatal: sgpr zero component: ' // text)
          endif
          sgnwc=sgnwc+(sgnw(i6,i6)/sgpr(i6,i6)-1d0)**2
       enddo
       !frs 18.04.2019     If( Abs(mytrace(sgnw)/trpr - 1d0) + Abs(etaf - etapr) .lt. 2d-6) then
       If(sgnwc .lt. 4d-12) then!frs 18.04.2019
          goto 200
       else
          apr = (1d0 - dampf)*apr + dampf*avr;
          sgpr = (1d0 - dampf)*sgpr + dampf*sgnw;
          !        etapr = (1d0 - dampeta)*etapr + dampeta*etaf !frs 13.06.18
          sgin=sgpr
       endif
    enddo
200 continue
    fcore=etaf; print*,"etaf: ",etaf; print*,"Iteration: ",ii-1;
    sgnwa = matmul(sgnw,sr)
    reeig=0d0
    aieig=0d0
    am=0d0
    call ladeig(sgnwa,reeig,aieig,am)
    print*,"beta11:   ",  am(1,1)*am(1,1)+am(1,2)*am(1,2)
    print*,"beta12:   ",  am(1,3)*am(1,3)+am(1,4)*am(1,4)
    write(90,'(3(e25.17,1x))') sgnw(1,1),sgnw(3,3),sgnw(5,5)
    write(91,'(6(e25.17,1x))') ((sgnw(j,i), i=1,6), j=1,6)
    sc3d_emit(1)=aieig(1)
    sc3d_emit(2)=aieig(3)
    sc3d_emit(3)=aieig(5)
    sc3d_emit(4)=sgnw(5,5)
    sc3d_emit(5)=sgnw(6,6)
    if(sc_3d_periodic) then
       nf1=get_file_unit(lu_max)
       open(nf1,file="rmatrix.dat",STATUS="OLD",err=1000)
       goto 2000
1000   call fort_fail('TRRUN: Fatal: ',                                     &
            'File: rmatrix.dat non-existing')
2000   continue
       do i=1,47
          read(nf1,*) ch
       enddo
       sgnw=0d0
       read(nf1,*) ch,dummy,(sgnw(i,1:idim),i=1,idim)
       close(nf1)
       reeig=0d0
       aieig=0d0
       am=0d0
       call ladeig(sgnw,reeig,aieig,am)
       am=TRANSPOSE(am)
       sgnw=0d0
       do i=1,idim
          do j=1,idim
             do k=1,3
                sgnw(i,j)=sgnw(i,j)+sc3d_emit(k)*(am(2*k-1,i)*am(2*k-1,j)+am(2*k,i)*am(2*k,j))
             enddo
          enddo
       enddo
    endif
  end subroutine SigmaMatrixFit6Dv7
  subroutine SCkick3Dv4(track,ktrack,fk)
    use SCdat
    use spch_bbfi
    implicit none
    integer ktrack,itrack
    double precision, external :: eInt
    double precision, external :: lInt
    double precision, external :: bips
    double precision track(6,*),fk
    double precision h1,h2,h3,h4,hf
    integer, parameter :: nstt=100
    integer i,k,l,m,n,ip,mmz,mmzopt,m1,n1,it,kmax,get_option
    double precision a,sqa,zdig,xory,z1,z2,s2,sm,sn,dpdz1,dpdz11,dpdz2&
         &,dpdz21,sig1,sig2,zmax,r,u,v,uz1,uz2,aum,sqaum,rad
    double precision pot,pot1,pot2,pot3,efx,efy,efx1,efy1,efx2,efy2,ir&
         &,lInts,eInts,st0,weight,at,tst(0:nstt),tfc(0:nstt),sqat,fnx(0:nstt)&
         &,fny(0:nstt),fnp(0:nstt),exp2,dpdz12,dpdz22
    double precision x,y,ctm,x2,y2,t
    double precision :: sigx=0d0,sigy=0d0,sigct=0d0
    double precision :: iB(0:mmax,0:mmax+1),iBu(0:mmum,0:mmum+1)
    character(20) text
    !********************************************************************************
    kmax = get_option('sc_mult_ord ')
    if(kmax.gt.kmaxo) then
       write(text, '(1p,i4)') kmaxo
       call fort_fail('SCkick:', 'kmax > kmaxo =' // text)
    endif
    sigx=scsigx(i_spch); sigy=scsigy(i_spch); sigct=scsigz(i_spch);
    if(sigx.eq.0d0.or.sigy.eq.0d0.or.sigct.eq.0d0) call fort_fail('SCkick:', 'Some Sigma values = 0 ')
    if(sigy > sigx) then
       xory=0d0
    else
       xory=1d0
    endif
    if(xory==1d0) then
       sig1=sigx; sig2=sigy;
    else
       sig1=sigy; sig2=sigx;
    endif
    a=(sig2/sig1)**2-1d0; sqa=Sqrt(1d0+a);
    !(***** power expansion coefficients for small displacements *****)
    if(Abs(a)>5d-1) then
       h1=dble(mmax+1)
       h2=5d-1
       h3=dble(mmax+2)
       h4=-a
       !       call hygfx(mmax+1,1d0/2d0,mmax+2,-a,hf)
       call hygfx(h1,h2,h3,h4,hf)
       iB(mmax,0)=hf/dble(mmax+1)
    else
       iB(mmax,0)=bips(a, mmax);
    endif
    Do k=1,mmax
       m=mmax-k;iB(m,0)=(sqa-a*(dble(m)+3d0/2d0)*iB(m+1,0))/dble(m+1);
       Do l=0,mmax-k+1
          iB(m,l+1)=(1d0/sqa/(1d0+a)**l-(dble(m-l)+1d0/2d0)*iB(m,l))&
               &/(dble(l)+1d0/2d0)
       enddo
    enddo
    !(***** expansion coefficients for intermediate displacements *****)
    aum=a*um;sqaum=Sqrt(1d0+aum);
    iBu(mmum,0)=bips(aum,mmum);
    Do k=1,mmum
       m=mmum-k;
       iBu(m,0)=(sqaum-aum*(dble(m)+3d0/2d0)*iBu(m+1,0))/(dble(m)+1d0);
       Do l=0,mmum-k+1
          iBu(m,l+1)=(1d0/sqaum/(1d0+aum)**l-(dble(m-l)+1d0/2d0)*&
               &iBu(m,l))/(dble(l)+1d0/2d0)
       enddo
    enddo
!!! new stuff 05.18 for intermediate displacements !!!

    pot3=2d0*log(1d0+sqa)/(1d0+sqaum)+log(um)

    st0 = (1d0 - um)/dble(Nstt);
    Do it=0,Nstt
       tst(it) = um + st0*dble(it)
       at=a*tst(it);
       tfc(it)=tst(it)/(1d0 + at);
       if(it.eq.0.or.it.eq.Nstt) then
          weight=st0/3d0
       else
          weight=2d0*(dble(mod(it,2))+1d0)*st0/3d0
       endif
       sqat = Sqrt(1d0+at);
       fnx(it)=weight/sqat;
       fny(it)=weight/sqat**3;
       fnp(it)=weight/sqat/tst(it)
    enddo

!!! end of new stuff 05.18 for intermediate displacements !!!

    !(***** cycle over particles *****)
    !hrr
    !$OMP PARALLEL DO PRIVATE(itrack,x,y,ctm,zdig,mmzopt,mmz,z1,z2,m,Sm,Sn), &
    !$OMP& PRIVATE(n,pot,m1,dpdz1,n1,dpdz2,efx,efy,zmax,r,u,v,lInts,k,iR,it), &
    !$OMP& PRIVATE(eInts,uz1,uz2,dpdz11,dpdz21,rad,dpdz12,dpdz22,exp2,pot1,pot2)
    do itrack=1,ktrack
       x=track(1,itrack); y=track(3,itrack); ctm=track(5,itrack); !(* actually it is T=-ct *)
       zdig=Sqrt((x/sigx)**2+(y/sigy)**2);
       If(zdig <= dres) then
          !(* small displacements *)
          mmzopt=Ceiling(6.439952123143645d0 +1.9839807379159704d0*zdig+&
               &1.3522706712502284d0*zdig**2);
          mmz=Min(mmzopt, mmax);
          If(xory==1d0) then
             z1=(x/sig1)**2/2d0; z2=(y/sig1)**2/2d0
          else
             z1=(y/sig1)**2/2d0; z2=(x/sig1)**2/2d0;
          endif
          If(z1<epsz) z1=epsz;
          If(z2<epsz) z2=epsz;
          m=mmz;
          Sm=0d0;
          Do While(m>0)
             n=mmz-m; Sn=0d0;
             Do While(n>=0)
                Sn=-Sn*z2/dble(n+1)+iB(m+n-1,n)
                n = n-1;
             enddo
             Sm=-Sm*z1/dble(m+1)+Sn
             m = m-1
          enddo
          n=mmz;Sn=0d0;
          Do While(n>0)
             Sn=-Sn*z2/dble(n+1)+iB(n-1,n)
             n = n-1
          enddo
          pot=-Sm*z1-Sn*z2;
          m1=mmz-1; Sm=0d0;
          Do While(m1>=0)
             n=mmz-m1-1; Sn=0d0;
             Do While(n>=0)
                Sn=-Sn*z2/dble(n+1)+iB(m1+n,n)
                n = n-1
             enddo
             Sm=-Sm*z1/dble(m1+1)+Sn
             m1 = m1-1
          enddo
          dpdz1=-Sm;
          m=mmz-1; Sm=0d0;
          Do While(m>=0)
             n1=mmz-m-1; Sn=0d0;
             Do While(n1>=0)
                Sn=-Sn*z2/dble(n1+1)+iB(m+n1,n1+1)
                n1 = n1-1
             enddo
             Sm=-Sm*z1/(m+1)+Sn
             m = m-1
          enddo
          dpdz2=-Sm;
          efx=-(xory*dpdz1+(1d0-xory)*dpdz2)*x/sig1**2;
          efy=-(xory*dpdz2+(1d0-xory)*dpdz1)*y/sig1**2
       else
          zmax=Sqrt(x**2+y**2)/sig1;
          If(zmax>6) then
             !(* large displacements *)
             r=sigy/sigx;
             u=x/sigx;v=y/sigx;iR=1d0/(u**2+v**2);
             lInts=0d0
             do k=1,kmax
                lInts=lInts+lInt(k,r,u,v)*iR**(2*k)
             enddo
             pot=Log((1d0+r)**2*iR/2d0)-EulerGamma-lInts/Pi;
             eInts=0d0
             do k=1,kmax
                eInts=eInts+eInt(k,r,u,v)*iR**(2*k+1)
             enddo
             efx=(2d0*u*iR+eInts/Pi)/sigx;
             r=sigx/sigy;
             u=y/sigy;v=x/sigy;iR=1d0/(u**2+v**2);
             eInts=0d0
             do k=1,kmax
                eInts=eInts+eInt(k,r,u,v)*iR**(2*k+1)
             enddo
             efy=(2d0*u*iR+eInts/Pi)/sigy
          else
             !(* Intermediate displacements *)
             mmz=Min(3+NINT(1.65d0*zmax),mmum);
             If(xory==1d0) then
                z1=(x/sig1)**2/2d0;z2=(y/sig1)**2/2d0
             else
                z1=(y/sig1)**2/2d0;z2=(x/sig1)**2/2d0;
             endif
             If(z1<epsz) z1=epsz;If(z2<epsz) z2=epsz;
             uz1=um*z1; uz2=um*z2;
             m=mmz; Sm=0d0;
             Do While(m>0)
                n=mmz-m; Sn=0d0;
                Do While(n>=0)
                   Sn=-Sn*uz2/dble(n+1)+iBu(m+n-1,n);
                   n = n-1;
                enddo
                Sm=-Sm*uz1/dble(m+1)+Sn;
                m = m-1;
             enddo
             n=mmz;Sn=0d0;
             Do While(n>0)
                Sn=-Sn*uz2/dble(n+1)+iBu(n-1,n);
                n = n-1;
             enddo
             pot1=-Sm*uz1-Sn*uz2;
             m1=mmz-1; Sm=0d0;
             Do While(m1>=0)
                n=mmz-m1-1; Sn=0d0;
                Do While(n>=0)
                   Sn=-Sn*uz2/dble(n+1)+iBu(m1+n,n);
                   n = n-1;
                enddo
                Sm=-Sm*uz1/dble(m1+1)+Sn;
                m1 = m1-1;
             enddo
             dpdz11=-Sm*um;
             m=mmz-1; Sm=0d0;
             Do While(m>=0)
                n1=mmz-m-1; Sn=0d0;
                Do While(n1>=0)
                   Sn=-Sn*uz2/dble(n1+1)+iBu(m+n1,n1+1);
                   n1 = n1-1;
                enddo
                Sm=-Sm*uz1/(m+1)+Sn;
                m = m-1;
             enddo
             dpdz21=-Sm*um;
!!! new stuff 05.2018 for intermediate displacements !!!**)

             pot2 = 0d0; dpdz12 = 0d0; dpdz22 = 0d0;
             Do it=0,Nstt
                exp2=1d0/Exp(z1*tst(it)+z2*tfc(it));
                dpdz12 = dpdz12 - fnx(it)*exp2;
                dpdz22 = dpdz22 - fny(it)*exp2;
                pot2 = pot2 + fnp(it)*exp2
             enddo
             efx=-(xory*(dpdz11+dpdz12)+(1-xory)*(dpdz21+dpdz22))*&
                  x/sig1**2;
             efy=-(xory*(dpdz21+dpdz22)+(1-xory)*(dpdz11+dpdz12))*&
                  y/sig1**2;
             pot=pot1+pot2+pot3
          endif
       endif
       rad=fk*Exp(-(ctm/sigct)**2/2d0)
       track(2,itrack)=track(2,itrack)+rad*efx
       track(4,itrack)=track(4,itrack)+rad*efy
       track(6,itrack)=track(6,itrack)+rad*pot*ctm/sigct**2
    enddo
    !$OMP END PARALLEL DO
  end subroutine SCkick3Dv4
  subroutine mymap()
    use spch_bbfi
    implicit none
    integer i,j,k,m2,nf1
    double precision s,ree(6,6),ret(6,6),ret_trans(6,6),ret2(6,6),ret2_trans(6,6)
    character(240) ch
    character(24) name
    integer, external :: get_file_unit
    character(20) text
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nf1=get_file_unit(lu_max)
    open(nf1,file="my_sector.tfs",STATUS="OLD",err=1000)

    do i=1,8
       read(nf1,'(a)') ch
    enddo
    ret=0
    do i=1,6
       ret(i,i)=1d0
    enddo
    ret2=0
    do i=1,6
       ret2(i,i)=1d0
    enddo
    m2=0
    do
       read(nf1,*,end=100) name,s,((ree(j,k),k=1,6),j=1,6)
       if(name(1:7).eq."SPCH_BB") then
          m2=m2+1
          ret_trans=TRANSPOSE(ret)
          ret2_trans=TRANSPOSE(ret2)
          sc_sect_map(m2,1:6,1:6)=ret(1:6,1:6)
          trans_sc_sect_map(m2,1:6,1:6)=ret_trans(1:6,1:6)
          sect_map(m2,1:6,1:6)=ret2(1:6,1:6)
          trans_sect_map(m2,1:6,1:6)=ret2_trans(1:6,1:6)
          ret=ree
          ret2=0d0
          do i=1,6
             ret2(i,i)=1d0
          enddo
       else
          ret=matmul(ree,ret)
          ret2=matmul(ree,ret2)
       endif
    enddo
100 continue
    if(N_spch.ne.m2) then
       write(text, '(1p,a11,i6,1x,i6)') "N_spch m2: ",N_spch, m2
       call fort_fail('MYMAP: Fatal: ',                                  &
            'Wrong number of SC kicks: ' // text)
    endif
    goto 2000
1000 continue
    call fort_fail('TRRUN: Fatal: ',                                     &
         'File: my_sector.tfs non-existing')
2000 continue
  end subroutine mymap
  subroutine SigmaTransport(sgpr)
    use spch_bbfi
    implicit none
    integer i,j,k
    double precision sgpr(6,6), sgps(6,6),ret(6,6),ret_trans(6,6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N_spch
       ret(1:6,1:6)=sc_sect_map(i,1:6,1:6)
       ret_trans(1:6,1:6)=trans_sc_sect_map(i,1:6,1:6)
       sgps=0d0
       sgps=matmul(sgpr,ret_trans)
       sgpr=matmul(ret,sgps)
       scsigx(i)=sqrt(sgpr(1,1))
       scsigy(i)=sqrt(sgpr(3,3))
       scsigz(i)=sqrt(sgpr(5,5))
    enddo
  end subroutine SigmaTransport
  subroutine SigmaTransport2(sgpr)
    use spch_bbfi
    implicit none
    integer i,j,k
    double precision sgpr(6,6), sgps(6,6),ret(6,6),ret_trans(6,6),sigx,sigy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,N_spch
       ret(1:6,1:6)=sect_map(i,1:6,1:6)
       ret_trans(1:6,1:6)=trans_sect_map(i,1:6,1:6)
       sgps=0d0
       sgps=matmul(sgpr,ret_trans)
       sgpr=matmul(ret,sgps)
       sigx=sqrt(sgpr(1,1))
       sigy=sqrt(sgpr(3,3))
       sc_map(i,2,1) = sc_intstr*sc_charge(i)/(sigx*(sigx+sigy))
       sc_map(i,4,3) = sc_intstr*sc_charge(i)/(sigy*(sigx+sigy))
       ret(1:6,1:6)=sc_map(i,1:6,1:6)
       ret_trans=TRANSPOSE(ret)
       sgps=0d0
       sgps=matmul(sgpr,ret_trans)
       sgpr=matmul(ret,sgps)
       scsigx(i)=sqrt(sgpr(1,1))
       scsigy(i)=sqrt(sgpr(3,3))
       scsigz(i)=sqrt(sgpr(5,5))
    enddo
  end subroutine SigmaTransport2

end module SpaceCharge

subroutine hygfx ( a, b, c, x, hf )

  !*****************************************************************************80
  !
  !! HYGFX evaluates the hypergeometric function F(A,B,C,X).
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to
  !    incorporate this routine into a user program that the copyright
  !    is acknowledged.
  !
  !  Modified:
  !
  !    08 September 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) A, B, C, X, the arguments of the function.
  !    C must not be equal to a nonpositive integer.
  !    X < 1.
  !
  !    Output, real HF, the value of the function.
  !
  implicit none

  double precision a
  double precision a0
  double precision aa
  double precision b
  double precision bb
  double precision c
  double precision c0
  double precision c1
  double precision, parameter :: el = 0.5772156649015329D+00
  double precision eps
  double precision f0
  double precision f1
  double precision g0
  double precision g1
  double precision g2
  double precision g3
  double precision ga
  double precision gabc
  double precision gam
  double precision gb
  double precision gbm
  double precision gc
  double precision gca
  double precision gcab
  double precision gcb
  double precision gm
  double precision hf
  double precision hw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm
  double precision pa
  double precision pb
  double precision, parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision r0
  double precision r1
  double precision rm
  double precision rp
  double precision sm
  double precision sp
  double precision sp0
  double precision x
  double precision x1

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 .or. l1 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'HYGFX - Fatal error!'
     write ( *, '(a)' ) '  The hypergeometric series is divergent.'
     return
  end if

  if ( 0.95D+00 < x ) then
     eps = 1.0D-08
  else
     eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

     hf = 1.0D+00
     return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

     call gamma ( c, gc )
     call gamma ( c - a - b, gcab )
     call gamma ( c - a, gca )
     call gamma ( c - b, gcb )
     hf = gc * gcab /( gca *gcb )
     return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

     g0 = sqrt ( pi ) * 2.0D+00**( - a )
     call gamma ( c, g1 )
     call gamma ( 1.0D+00 + a / 2.0D+00 - b, g2 )
     call gamma ( 0.5D+00 + 0.5D+00 * a, g3 )
     hf = g0 * g1 / ( g2 * g3 )
     return

  else if ( l2 .or. l3 ) then

     if ( l2 ) then
        nm = int ( abs ( a ) )
     end if

     if ( l3 ) then
        nm = int ( abs ( b ) )
     end if

     hf = 1.0D+00
     r = 1.0D+00

     do k = 1, nm
        r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * x
        hf = hf + r
     end do

     return

  else if ( l4 .or. l5 ) then

     if ( l4 ) then
        nm = int ( abs ( c - a ) )
     end if

     if ( l5 ) then
        nm = int ( abs ( c - b ) )
     end if

     hf = 1.0D+00
     r  = 1.0D+00
     do k = 1, nm
        r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * x
        hf = hf + r
     end do
     hf = ( 1.0D+00 - x )**( c - a - b ) * hf
     return

  end if

  aa = a
  bb = b
  x1 = x
  !
  !  WARNING: ALTERATION OF INPUT ARGUMENTS A AND B, WHICH MIGHT BE CONSTANTS.
  !
  if ( x < 0.0D+00 ) then
     x = x / ( x - 1.0D+00 )
     if ( a < c .and. b < a .and. 0.0D+00 < b ) then
        a = bb
        b = aa
     end if
     b = c - b
  end if

  if ( 0.75D+00 <= x ) then

     gm = 0.0D+00

     if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

        m = int ( c - a - b )
        call gamma ( a, ga )
        call gamma ( b, gb )
        call gamma ( c, gc )
        call gamma ( a + m, gam )
        call gamma ( b + m, gbm )
        call psi ( a, pa )
        call psi ( b, pb )

        if ( m /= 0 ) then
           gm = 1.0D+00
        end if

        do j = 1, abs ( m ) - 1
           gm = gm * j
        end do

        rm = 1.0D+00
        do j = 1, abs ( m )
           rm = rm * j
        end do

        f0 = 1.0D+00
        r0 = 1.0D+00
        r1 = 1.0D+00
        sp0 = 0.0D+00
        sp = 0.0D+00

        if ( 0 <= m ) then

           c0 = gm * gc / ( gam * gbm )
           c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

           do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                   / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
           end do

           do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
                   + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / dble(k)
           end do

           f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
           hw = f1

           do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
                   + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                 sm = sm + ( 1.0D+00 - a ) &
                      / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
                      + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do

              rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

              r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
                   / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                 exit
              end if

              hw = f1

           end do

           hf = f0 * c0 + f1 * c1

        else if ( m < 0 ) then

           m = - m
           c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
           c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

           do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
                   / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
           end do

           do k = 1, m
              sp0 = sp0 + 1.0D+00 / dble ( k)
           end do

           f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )

           do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) &
                   / ( k * ( a + k - 1.0D+00 ) ) &
                   + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                 sm = sm + 1.0D+00 / dble ( j + k)
              end do

              rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                   / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
                 exit
              end if

              hw = f1

           end do

           hf = f0 * c0 + f1 * c1

        end if

     else

        call gamma ( a, ga )
        call gamma ( b, gb )
        call gamma ( c, gc )
        call gamma ( c - a, gca )
        call gamma ( c - b, gcb )
        call gamma ( c - a - b, gcab )
        call gamma ( a + b - c, gabc )
        c0 = gc * gcab / ( gca * gcb )
        c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
        hf = 0.0D0
        r0 = c0
        r1 = c1

        do k = 1, 250

           r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
                / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

           r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
                / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

           hf = hf + r0 + r1

           if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
              exit
           end if

           hw = hf

        end do

        hf = hf + c0 + c1

     end if

  else

     a0 = 1.0D+00

     if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

        a0 = ( 1.0D+00 - x )**( c - a - b )
        a = c - a
        b = c - b

     end if

     hf = 1.0D+00
     r = 1.0D+00

     do k = 1, 250

        r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
             / ( k * ( c + k - 1.0D+00 ) ) * x

        hf = hf + r

        if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
           exit
        end if

        hw = hf

     end do

     hf = a0 * hf

  end if

  if ( x1 < 0.0D+00 ) then
     x = x1
     c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
     hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'HYGFX - Warning!'
     write ( *, '(a)' ) '  A large number of iterations were needed.'
     write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end subroutine hygfx

subroutine gamma ( x, ga )

  !*****************************************************************************80
  !
  !! GAMMA evaluates the Gamma function.
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to
  !    incorporate this routine into a user program that the copyright
  !    is acknowledged.
  !
  !  Modified:
  !
  !    08 September 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument.
  !    X must not be 0, or any negative integer.
  !
  !    Output, real ( kind = 8 ) GA, the value of the Gamma function.
  !
  implicit none

  double precision, dimension ( 26 ) :: g = (/ &
       1.0D+00, &
       0.5772156649015329D+00, &
       -0.6558780715202538D+00, &
       -0.420026350340952D-01, &
       0.1665386113822915D+00, &
       -0.421977345555443D-01, &
       -0.96219715278770D-02, &
       0.72189432466630D-02, &
       -0.11651675918591D-02, &
       -0.2152416741149D-03, &
       0.1280502823882D-03, &
       -0.201348547807D-04, &
       -0.12504934821D-05, &
       0.11330272320D-05, &
       -0.2056338417D-06, &
       0.61160950D-08, &
       0.50020075D-08, &
       -0.11812746D-08, &
       0.1043427D-09, &
       0.77823D-11, &
       -0.36968D-11, &
       0.51D-12, &
       -0.206D-13, &
       -0.54D-14, &
       0.14D-14, &
       0.1D-15 /)
  double precision ga
  double precision gr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  double precision, parameter :: pi = 3.141592653589793D+00
  double precision r
  double precision x
  double precision z

  if ( x == aint ( x ) ) then

     if ( 0.0D+00 < x ) then
        ga = 1.0D+00
        m1 = int ( x ) - 1
        do k = 2, m1
           ga = ga * k
        end do
     else
        ga = 1.0D+300
     end if

  else

     if ( 1.0D+00 < abs ( x ) ) then
        z = abs ( x )
        m = int ( z )
        r = 1.0D+00
        do k = 1, m
           r = r * ( z - dble ( k))
        end do
        z = z - dble ( m)
     else
        z = x
     end if

     gr = g(26)
     do k = 25, 1, -1
        gr = gr * z + g(k)
     end do

     ga = 1.0D+00 / ( gr * z )

     if ( 1.0D+00 < abs ( x ) ) then
        ga = ga * r
        if ( x < 0.0D+00 ) then
           ga = - pi / ( x* ga * sin ( pi * x ) )
        end if
     end if

  end if

  return
end subroutine gamma
subroutine psi ( x, ps )

  !*****************************************************************************80
  !
  !! PSI computes the PSI function.
  !
  !  Licensing:
  !
  !    The original FORTRAN77 version of this routine is copyrighted by
  !    Shanjie Zhang and Jianming Jin.  However, they give permission to
  !    incorporate this routine into a user program that the copyright
  !    is acknowledged.
  !
  !  Modified:
  !
  !    08 September 2007
  !
  !  Author:
  !
  !    Original FORTRAN77 by Shanjie Zhang, Jianming Jin.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument.
  !
  !    Output, real ( kind = 8 ) PS, the value of the PSI function.
  !
  implicit none

  double precision, parameter :: a1 = -0.83333333333333333D-01
  double precision, parameter :: a2 =  0.83333333333333333D-02
  double precision, parameter :: a3 = -0.39682539682539683D-02
  double precision, parameter :: a4 =  0.41666666666666667D-02
  double precision, parameter :: a5 = -0.75757575757575758D-02
  double precision, parameter :: a6 =  0.21092796092796093D-01
  double precision, parameter :: a7 = -0.83333333333333333D-01
  double precision, parameter :: a8 =  0.4432598039215686D+00
  double precision, parameter :: el = 0.5772156649015329D+00
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  double precision, parameter :: pi = 3.141592653589793D+00
  double precision ps
  double precision s
  double precision x
  double precision x2
  double precision xa

  xa = abs ( x )
  s = 0.0D+00

  if ( x == aint ( x ) .and. x <= 0.0D+00 ) then

     ps = 1.0D+300
     return

  else if ( xa == aint ( xa ) ) then

     n = int ( xa )
     do k = 1, n - 1
        s = s + 1.0D+00 / dble ( k)
     end do

     ps = - el + s

  else if ( xa + 0.5D+00 == aint ( xa + 0.5D+00 ) ) then

     n = int ( xa - 0.5D+00 )

     do k = 1, n
        s = s + 1.0D+00 / dble ( 2 * k - 1)
     end do

     ps = - el + 2.0D+00 * s - 1.386294361119891D+00

  else

     if ( xa < 10.0D+00 ) then

        n = 10 - int ( xa )
        do k = 0, n - 1
           s = s + 1.0D+00 / ( xa + dble ( k) )
        end do

        xa = xa + dble ( n)

     end if

     x2 = 1.0D+00 / ( xa * xa )

     ps = log ( xa ) - 0.5D+00 / xa + x2 * ((((((( &
          a8   &
          * x2 + a7 ) &
          * x2 + a6 ) &
          * x2 + a5 ) &
          * x2 + a4 ) &
          * x2 + a3 ) &
          * x2 + a2 ) &
          * x2 + a1 )

     ps = ps - s

  end if

  if ( x < 0.0D+00 ) then
     ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D+00 / x
  end if

  return
end subroutine psi
