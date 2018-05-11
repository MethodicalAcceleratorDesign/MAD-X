module SpaceCharge

  implicit none
  private

  logical :: sc_chrom_fix

  exit_loss_turn = get_option('exit_loss_turn ') .ne. 0
  bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
  checkpnt_restart = get_value('run ', 'checkpnt_restart ') .ne. zero
  emittance_update = get_option('emittance_update ') .ne. 0
  virgin_state = get_value('run ', 'virgin_state ') .ne. zero

  
  public :: sc_chrom_fix


contains
  
  subroutine Init(c)
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
       call ixy_calcs(betas, orbit0, z,       &
            betx_start, bety_start, &
            alfx_start, alfy_start, &
            gamx_start, gamy_start, &
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

  end subroutine Init

end module SpaceCharge
