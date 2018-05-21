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

contains

  subroutine Init_SC()

    integer, external :: get_option
    double precision, external :: get_value

    exit_loss_turn = get_option('exit_loss_turn ') .ne. 0
    bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
    checkpnt_restart = get_value('run ', 'checkpnt_restart ') .ne. zero
    emittance_update = get_option('emittance_update ') .ne. 0
    virgin_state = get_value('run ', 'virgin_state ') .ne. zero

    !! ALSC: I think init should reset to zero these varibles
    ex_rms0=zero;
    ey_rms0=zero;
    sigma_p0=zero;
    sigma_z0=zero;
    
  end subroutine Init_SC

  subroutine SC_Update(betas, orbit, z,          &
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

end module SpaceCharge
