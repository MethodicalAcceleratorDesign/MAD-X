
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
  public ! hrr debug

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
! double precision, save :: dx_start=0d0,   dpx_start=0d0  !hrr debug
! double precision, save :: dy_start=0d0,   dpy_start=0d0  !hrr debug

  !VVK 20100321 -------------------------------------------------
  integer :: i_part                     ! local counter
  double precision  :: Summ_t_mean      ! local for mean value
  double precision  :: Summ_t_square    ! local for rms value
  !-------------------------------------------------------------------

! private :: table_input
! private :: ixy_calcs
! private :: ixy_fitting

end module SpaceCharge

module SpaceCharge2

  double precision, save :: dx_start=0d0,   dpx_start=0d0 !hrr debug
  double precision, save :: dy_start=0d0,   dpy_start=0d0 !hrr debug

end module SpaceCharge2



