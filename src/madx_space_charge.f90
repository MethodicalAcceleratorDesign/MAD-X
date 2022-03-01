
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

Module insert_madx_spch
  use name_lenfi
  implicit none
  Character (80) :: twiss_element_list_filename, &
                  title_madx_script, filename_ring_elements_seq
  Character (name_len) :: element_label(10000), element_keyword(10000)
  double precision ::           element_length(10000)
  Integer :: N_elements

  double precision :: ring_length

END Module insert_madx_spch
Subroutine Read_twiss_element_list

! The elenment list must be created with 
! the following MADX commands:
! SELECT, flag=twiss,clear;
! SELECT, flag=twiss, column=name, keyword, L;
! TWISS, file=deb_twiss_element_list.txt;
!

  USE insert_madx_spch
  use spch_bbfi
  implicit none

  Character (80) :: twiss_table_line
  Integer :: i_line, Label_length, i_element, &
       start_pos, end_pos, &
       keyword_start, keyword_end, &
       LENGTH_pos,nf1
  integer, external :: get_file_unit
  Character (30) :: string_temp
  Logical :: debug

  Namelist /input_filename/ twiss_element_list_filename

  Namelist /madx_job_title/ title_madx_script

  Namelist /ring_seq_filename/ filename_ring_elements_seq

  debug=.FALSE.
!debug=.TRUE.

  nf1=get_file_unit(lu_max)
  OPEN (nf1, File='madx_spch_input.madx', &
       STATUS='OLD', FORM='FORMATTED')

  READ(nf1,input_filename) 
  Print *, 'input_filename:'
  WRITE(*, input_filename)

  READ(nf1,madx_job_title) 
  Print *, 'madx_job_title:'
  WRITE(*, madx_job_title)

  READ(nf1,ring_seq_filename) 
  Print *, 'ring_seq_filename:'
  WRITE(*, ring_seq_filename)

  CLOSE(nf1)

  OPEN (nf1, File=TRIM(twiss_element_list_filename), &
       STATUS='OLD', FORM='FORMATTED')

  i_line=0
  i_element=0

  Read_lines: DO
  i_line=i_line+1
  read(nf1,FMT="(A80)", END=10) twiss_table_line

  if(debug) then
     Print *, TRIM(twiss_table_line)
     Print *, 'i_line=', i_line,' (2)=', twiss_table_line(2:2)
  endif !debug

      !  general twiss param
  if(twiss_table_line(1:1) .EQ. '@') then

	     ! LENGTH =============================================== 
     LENGTH_pos=INDEX(TRIM(twiss_table_line),'LENGTH')
     if (LENGTH_pos .EQ. 3) then
        read(twiss_table_line, *) string_temp, string_temp, &
             string_temp, ring_length
        print *, 'ring_length=', ring_length 
     endif ! LENGTH =========================================

  endif !  general twiss param


      ! elements parameters
  if(twiss_table_line(2:2) .EQ. '"') then

     start_pos=INDEX(TRIM(twiss_table_line),'$START')
     end_pos=INDEX(TRIM(twiss_table_line),'$END')

     if(debug) then
        Print *, 'start_pos=', start_pos, &
             '  end_pos=', end_pos
     endif !debug

     if(start_pos.EQ.0 .AND. end_pos.EQ.0) then
        i_element=i_element+1

	        !Label_length=LEN_TRIM( twiss_table_line)
        Label_length=INDEX(TRIM(twiss_table_line),'" ')

        if(debug) &
             Print *, 'Label_length=', Label_length

        element_label(i_element)= &
             twiss_table_line(3:Label_length-1)

			! relative position:
        keyword_start=INDEX(TRIM(twiss_table_line(Label_length+1:)),' "')
            ! absolute position:
        keyword_start=Label_length+1+keyword_start

            ! relative position:
        keyword_end=INDEX(TRIM(twiss_table_line(keyword_start+1:)),'"')
            ! absolute position:
        keyword_end=keyword_start+keyword_end

        element_keyword(i_element)= &
             twiss_table_line(keyword_start+1:keyword_end-1)

        read(twiss_table_line,*) string_temp, string_temp, &
			                         element_length(i_element)

        if(debug) then
           Print *, 'keyword_start, end=', keyword_start, keyword_end
           Print *, TRIM(twiss_table_line(keyword_start:keyword_end))
        endif ! debug

        Print *, 'element i=', i_element, &
             '  label=', TRIM(element_label(i_element)),  &
             '  keyword=', TRIM(element_keyword(i_element)), &
             '  l=', element_length(i_element)

     endif

  endif

ENDDO Read_lines


10 Continue
N_elements=i_element
Print *, 'Original ring has N_elements=', N_elements
CLOSE(nf1)

END Subroutine Read_twiss_element_list
Subroutine Write_Insert_SpCh
  USE insert_madx_spch
  use spch_bbfi
  implicit none

  integer i_element,nf1
  integer, external :: get_file_unit

  nf1=get_file_unit(lu_max)
  OPEN (nf1, File='spcharge_exec_insert_elements.madx', &
         STATUS='UNKNOWN', FORM='FORMATTED')

  DO i_element=1,N_elements
     Write (nf1, 10) TRIM(element_label(i_element)), &
	               TRIM(element_label(i_element))
10   Format (1x,'exec, Insert_SpCh (',(A),', ',(A), &
   '_count, S_position, L_div_max, sc_matrix, SC_count);')

  ENDDO
  Write(nf1,*) 'RETURN;'

  CLOSE(nf1)
END Subroutine Write_Insert_SpCh
SUBROUTINE Write_drifts_with_lengths
  USE insert_madx_spch
  use spch_bbfi
  implicit none

  integer i_element, drift_pos,nf1
  integer, external :: get_file_unit

  nf1=get_file_unit(lu_max)
  OPEN (nf1, File='spcharge_drift_elements.madx', &
       STATUS='UNKNOWN', FORM='FORMATTED')

  DO i_element=1,N_elements

     drift_pos=INDEX(element_keyword(i_element),'DRIFT')
     !Print *, 'drift_pos=', drift_pos
     if( drift_pos .EQ. 1) then

        Write (nf1, 10) TRIM(element_label(i_element)),   &
	               element_length(i_element)
     endif ! if( drift_pos .EQ. 1)

10   Format (1x,(A),': DRIFT, L=',f24.18,';')

  ENDDO
  Write(nf1,*) 'RETURN;'

  CLOSE(nf1)
END SUBROUTINE Write_drifts_with_lengths
SUBROUTINE Write_baseelem_codes
  USE insert_madx_spch
  use spch_bbfi
  implicit none

  integer i_element,nf1
  integer, external :: get_file_unit

  nf1=get_file_unit(lu_max)
  OPEN (nf1, File='spcharge_baseelem_codes.madx', &
       STATUS='UNKNOWN', FORM='FORMATTED')

  DO i_element=1,N_elements
     Write (nf1, 10) TRIM(element_keyword(i_element)), &
          TRIM(element_label(i_element)),   &
          TRIM(element_label(i_element)),   &
          TRIM(element_label(i_element))

10   Format (1x,'exec, set_code_',(A),'(',(A), &
          '); value, code_',(A),'; ',(A),'_count=1;')

  ENDDO
  Write(nf1,*) 'RETURN;'

  CLOSE(nf1)
END SUBROUTINE Write_baseelem_codes
Subroutine Write_SpCh_ring_seq
  USE insert_madx_spch
  use spch_bbfi
  implicit none

  integer i_element,nf1
  integer, external :: get_file_unit

  nf1=get_file_unit(lu_max)
  OPEN (nf1, File='spcharge_ring_seq.madx', &
       STATUS='UNKNOWN', FORM='FORMATTED')


  Write (nf1, 101) TRIM(title_madx_script)
101 format (1x, 'TITLE, "',(A),'";')

  Write(nf1,102); 102 format(/ &
  '! Space Charge =================================================='/)
  Write(nf1,103); 103 format( &
       'sc_matrix: MATRIX, L=0.0;       ! Class for linear space-charge')
  Write(nf1,104); 104 format( &
       'spch_bb: BEAMBEAM;              ! Class for bi-Gaussian space-charge')
  Write(nf1,105); 105 format( &
  '! ================================================================'/)

  Write(nf1,106) TRIM(filename_ring_elements_seq)
106 format(1x,'CALL, FILE=',(A),';'/ &
         1x,'CALL FILE=spcharge_drift_elements.madx;'/)

! drifts because they removed from sequence


!===========================================================
! codes

  Write(nf1,107)
107 format(1x,'CALL, FILE=spcharge_baseelem_codes.madx;'/)
  Write(nf1,*)

  Write(nf1,108); 108 format( &
       '! --- Counters initialization -----------------------------------')
  Write(nf1,109); 109 format( &
       'S_position=0;  ! Entry position for a current position')
  Write(nf1,110); 110 format( &
       'SC_count=1;    ! Counter for Space-charge'/)

! Ring length =====================================================

  Write(nf1,111) ring_length
111 format(1x,'L_ring=',f24.18,';'/)

!   Start sequence
  Write(nf1,112); 112 format( &
       ' ring: SEQUENCE,              L = L_ring;'/ &
       ' CALL FILE=spcharge_exec_insert_elements.madx;'/ &
       ' ENDSEQUENCE;'//' SAVE,SEQUENCE=ring,file=ring_seq.txt;'/ &
       ' CALL,file="ring_seq.txt";'/ &
       ' USE, PERIOD= ring;'/)

  Write(nf1,*) 'RETURN;'

  CLOSE(nf1)

END Subroutine Write_SpCh_ring_seq

Subroutine SC_Setup_Procedure
  CALL Read_twiss_element_list
  CALL Write_Insert_SpCh
  Call Write_baseelem_codes
  CALL Write_drifts_with_lengths
  CALL Write_SpCh_ring_seq
end Subroutine SC_Setup_Procedure

module SpaceCharge

  use trackfi, only : max_part
  use math_constfi

  implicit none
  !!private
  public ! hrr Sep 2021

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

  double precision, save :: dx_start=0d0,   dpx_start=0d0 !hrr Sep 2021
  double precision, save :: dy_start=0d0,   dpy_start=0d0 !hrr Sep 2021

end module SpaceCharge2



