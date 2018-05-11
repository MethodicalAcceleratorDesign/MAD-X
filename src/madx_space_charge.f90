module SpaceCharge

  implicit none
  private

  logical :: sc_chrom_fix

  public :: sc_chrom_fix

contains
  
  subroutine Init_SC()
!!$    exit_loss_turn = get_option('exit_loss_turn ') .ne. 0
!!$    bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
!!$    checkpnt_restart = get_value('run ', 'checkpnt_restart ') .ne. zero
!!$    emittance_update = get_option('emittance_update ') .ne. 0
!!$    virgin_state = get_value('run ', 'virgin_state ') .ne. zero
  end subroutine Init_SC

end module SpaceCharge
