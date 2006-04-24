MODULE madx_ptc_track_run_module
  ! This module serve as a COMMON block as in  F77
  ! It contains variables which exchange data between
  ! SUBROUTINE ptc_track_run and called from it
  ! external subroutines calculating particle interactions
  USE madx_ptc_module , ONLY: dp, lp, lnv, &
  !                          ! shorts for <double precision>, <logical>, 0D0 etc.
                   doublenum ! am temprorary double number for I/O with C-procedures  
implicit none
  SAVE
  PRIVATE

  PUBLIC :: ptc_track_run ! Subroutine inside the module

  !------------------------------------------------------------------!
  !  Variables from input files and probably corrected by this code: !
  !------------------------------------------------------------------!
  INTEGER, PUBLIC :: icase_ptc    ! Phase-space (4, 5 or 6) in input command
  !                               ! is NOT actul number of variables
  INTEGER, PUBLIC :: turns        ! The current turn and Total number of turns

  LOGICAL(lp), PUBLIC :: &
       closed_orbit, &          ! Switch to turn on closed orbit calculation
       element_by_element, &    ! element-by-element tracking (not over-one-turn one
       Radiation_PTC, &         ! Radiation is internally done by PTC
       Radiation_model1_FZ, &   ! Radiation according to model 1 (done by FZ)
       Radiation_Energy_Loss, & !
       Radiation_Quad, &
       Space_Charge, &
       beam_envelope, &
       NORM_OUT              ! IF true THEN output tables with normalizad coord,
  !                       ELSE (default) => output tables with standard coord.

  INTEGER :: Normal_Order_n0 ! The order of Normal Form ( Normal_Order_n0=1 => for Linear)

  integer :: ptc_ffile    ! periodicity of printing coordinates (every ffile-th turn)
  !                       ! integer for keyword 'ffile' in 'ptc_track' command
  logical(lp) :: ptc_onetable ! logical for keyword 'onetable'
  !                       ! in 'ptc_track' command

  !--------------------------------------------------!
  ! variables calculated or called inside this code  !
  !--------------------------------------------------!

  INTEGER, PUBLIC :: npara        ! icase_PTC is NOT actul number of variables
  !                               ! out from PTC subr. <init> called in subr.
  !                               ! Get_map_from_NormalForm
  INTEGER, PUBLIC :: i_th_turn  ! The current turn number

  Real (dp), public :: &
       Energy_rest_MeV, &    ! the rest and total energy of particles
       Energy_total_MeV, &   !
       deltap, &             ! The relative momentum deviation for off-momentum particles (5D case)
       dt                    ! Energy error, divided by the ref.momentum times the velocity of light

  integer, PUBLIC :: &
       jmax_numb_particl_at_i_th_turn, &  ! "jmax" may be reduced by particle loss
       j_tot_numb_starting_particles      ! "j_tot" saves initial number of particles
  !                                       ! meanining: <j_total_number_starting_particles>

  integer, PUBLIC,allocatable :: jmax_all_turns_numb_part(:) ! save jmax for all turns including init.

  INTEGER, allocatable :: particle_ID(:)     !(1:N_particle_max) ! => part_id (in TRRUN.F)
  INTEGER, allocatable :: part_ID_turns(:,:) ! save particle ID for every i-th turns (i_turn,j_part)

  integer :: tot_segm_one_table   ! for table printing (like in trrun)
  integer :: segment_one_table    !

  real(dp), PUBLIC ::  x_coord_co(1:6)  ! Closed orbit  at Start of ring
  !                                     ! => x0(1:6) in ptc_track; orbit0(1:6) in TRRUN.F

  real(dp), PRIVATE, allocatable :: x_co_at_all_observ(:,:)
  !                                     ! to keep closed orbit at all observations for
  !                                     ! <element_by_element> tracking

  real(dp), PRIVATE, allocatable :: x_all_incl_co_at0(:,:,:) ! save x at START of the ring
  !                                               ! (1:6,0:turns,1:jmax_numb_particl_at_i_th_turn)

  logical(lp), public :: ptc_track_debug=.False. !.TRUE.  .False.

  logical(lp) :: last_table_line_out ! = last_out => flag to avoid double entry of last line

  INTEGER, ALLOCATABLE :: elem_number_at_observ(:)   ! the sequent number of a ring element at
  !                                                  ! the observation point

  REAL(dp), allocatable :: sum_length_at_observ(:)   ! the sum length for the current observ.
  !                                                  ! point from the ring start

  CHARACTER(24), ALLOCATABLE  :: name_el_at_obsrv(:) ! the name of an element at observation point
  !                                                  ! contrary to <character (16)> in c-code

  !real(kind(1d0)) :: dble_num_C  ! to use as a temprorary double number for I/O with C-procedures 

CONTAINS

  SUBROUTINE ptc_track_run(max_obs)

    !USE MADX_PTC_MODULE ==================================================================!
    USE  madx_ptc_module, ONLY: universe, my_ring, default, index, c_                      !
    !                                                                                      !
    USE  madx_ptc_module, ONLY: &  ! "LAYOUT type (ring) => double linked list,            !
         FIBRE, &                  !  whose nodes (elements=magnets) of type FIBRE"        !
         NORMALFORM, & ! type for normalform                                               !
         REAL_8, &     ! type for map                                                      !
         damap, &      ! type for diff algebra                                             !
         RADIATION, STOCH_IN_REC, & ! type for radiation with quadrupoles in PTC           !
         BEAMENVELOPE, ENV_8        ! For beam envelope                                    !
    ! ======== functions ==================================================================!
    USE  madx_ptc_module, ONLY: &                                                          !
         print, find_orbit, track,UPDATE_STATES, my_state, &                               !
         PRODUCE_APERTURE_FLAG, ANALYSE_APERTURE_FLAG, &                                   !
         kill, daprint, alloc, Get_one, &                                                  !
         assignment(=), operator(+), operator(*), operator(.sub.), &                       !
         Find_Envelope, &                                                                  !
         ! Coord_MAD_to_PTC, Coord_PTC_to_MAD,  & => at the end of this module             !
         write_closed_orbit,Convert_dp_to_dt,mytime                                        !
    !======================================================================================!
    USE  madx_ptc_module, ONLY: &                                                          !
         c_1d_7,c_1D3,one,two, twopi, zero                                                 !
    !======================================================================================!

    IMPLICIT NONE

    integer, intent (IN) :: max_obs ! the maximum number of observation points >=1
    !                               ! one point at the end (beginning)+
    !                               ! points given in input file by the command
    !                               ! "ptc_observe,place=mark";

    include 'name_len.fi'   ! integer name_len;  parameter(name_len=24)

    include 'bb.fi'         ! integer bbd_loc,bbd_cnt,bbd_flag,bbd_pos,bbd_max;
    !                       !uses bbd_pos                parameter(bbd_max=200)
    !                       !real(kind(1d0)) bb_kick
    !                       ! common/bbi/bbd_loc(bbd_max),bbd_cnt,bbd_flag,bbd_pos
    !                       !common/bbr/bb_kick(2,bbd_max)


    integer ::flag_index_ptc_aperture, why_ptc_aperture(9)

    real(dp) :: deltap0

    ! INTEGER, parameter :: N_particle_max=1000, & ! boundary of particle arrays
    ! This parameter is not used anymore, instead allocatable arrays are used in
    ! in the int.subr <ptc_track_ini_conditions> to define number of tracks (particles)

    ! INTEGER, parameter :: N_elements_max=1000    ! boundary of element arrays
    ! is not used, it is known in PTC as <my_ring%n>

    real(dp),allocatable :: &
         x_coord_incl_co (:,:) ! (1:6,1:N_particle_max) ! z(6,npart) in TRRUN.F

    real (dp)::  current_x_coord_incl_co(1:6)
    ! buffer for coordinates of the curent particle VK

    real (dp):: x_input_0_by_ptc_track_command(1:6)
    ! buffer for input coordinates of the track command line VK

    integer :: mad8_code_current_node  ! = code <in trrun.F>
    ! code_of_the_current_node_as_mad8
    REAL(dp) :: el_element_length       ! = el

    Integer, allocatable  :: last_turn_of_lost_particle(:)! (1:N_particle_max)  != last_turn

    REAL(dp) , allocatable :: &
         last_position_of_lost_particle(:), & !(1:N_particle_max), & != last_pos
         last_orbit_of_lost_particle (:,:) !(1:6,1:N_particle_max)   != last_orbit

    !   /* C routines called from Fortran and C */
    real(kind(1d0)), external :: get_value, get_variable ! external c-functions
    INTEGER, external :: get_option, &   !  int get_option(char*);
         restart_sequ, & !  restart beamline and return number of beamline node
         advance_node    !  advance to the next node in expanded sequence
    !                    !  =0 (end of range), =1 (else)

    REAL(KIND(1d0)), external :: node_value !/*returns value for parameter par of current element */

    EXTERNAL :: comm_para ! subroutine needed for LF95

    !k    character*12 char_a
    ! variables added by V.Kapin                ! = in <trrun.F>
    integer :: ptc_switch                       ! = switch

    real(dp) :: maxaper(1:6) ! moved from int.subr.ST
    integer :: k_th_coord !, j_th_particle
    logical(lp) :: return_from_subr_ptc_track

    REAL(dp) :: sum_length  ! = sum  <in trrun.F>

    TYPE(fibre), POINTER :: current ! F90 pointer to the CURRENT beamline (fibre) element

    TYPE(normalform) :: Normal_Form_N  !  n =>  Normal_Form_N - local name
    TYPE(real_8)     :: Map_Y(6)       !  y => Map_Y - local name
    TYPE (damap)     :: A_t_map        ! x_stdt=A_t_map*x_norm
    TYPE (damap)     :: A_t_map_rev    ! x_norm=A_t_map_rev*x_norm
    !                                  ! A_t_map_rev=A_t_map^(-1)

    !integer,  dimension(4) :: iia     ! iia(2) is the number of variables in tracking (4,..,6)

    !hbu
    real(dp) :: spos_current_position ! s-position   =spos
    integer ::  nlm_current_element_number    ! line position =nlm
    character*(name_len) el_name
    !hbu
    character*4 vec_names(7)
    !hbu
    data vec_names  / 'x', 'px', 'y', 'py', 't', 'pt','s' / ! MADX
    !data vec_names / 'x', 'px', 'y', 'py', 'pt', 't','s' / ! PTC has a reverse order for pt and t

    !k    data char_a / ' ' /
    !-------------------------------------------------------------------------

    ! default value
    return_from_subr_ptc_track=.FALSE. ! used to RETURN from this subr.

    element_by_element =.FALSE. ! default is one-turn tracking

    last_table_line_out = .false.   ! flag to avoid double entry of last line

    NORM_OUT=.False.    ! IF true THEN output tables with normalizad coord,
    ! ELSE (default) output tables with standard coord.
    Normal_Order_n0=1   !
    x_coord_co(:)=zero
                     x_input_0_by_ptc_track_command(:)=zero
    CALL Values_from_ptc_track_command ! Int. proc. in this subr.(see below in the file)
    ! read command line from input file :
    ! Example: ptc_track,icase=5, coord=0.01,0.,0.01,0.,0.,0.;
    ! Keywords: "icase", "turns", "closed_orbit", "deltap", "coord"
    ! Values for: icase,   turns ,  closed_orbit,   deltap0, coord.
    ! NORM_OUT=.TRUE.  ! IF true THEN output tables with normalizad coord,
    !                  !  ELSE (default) output tables with standard coord.

    ALLOCATE (jmax_all_turns_numb_part(0:turns)) ! to save save jmax for all turns including init.

    CALL Check_ptc_universe_and_layout ! (before_tracking) =>  Int.Proc.(this subr.)
    if (return_from_subr_ptc_track) RETURN      ! (see below in this file)

    Call Call_my_state_and_update_states ! parameter "deltap" is defined now
    !                                    ! icase_ptc is changed to correct value

    Print *;  Print *,'  ================================================================'
    Print *, '  ptc_track: The current dimensionality of the problem is icase=', icase_ptc
    Print *,'  ================================================================'; Print *;

    warn_coordinate_system_changed: IF((.NOT. mytime) .AND.(icase_ptc.gt.4) ) THEN
      CALL FORT_WARN('time=false => coord. system: {-pathlength, delta_p} ', &
                     'the table headers mean:  PT -> delta_p, T -> pathlength')
    ENDIF warn_coordinate_system_changed

    ! initialize the closed orbit coordinates  at START of the ring
    x_coord_co(:)=zero
    if (ptc_track_debug) print *, " x_coord_co(:)=zero = ",x_coord_co

    ! Closed_orbit_at_START:
    IF(closed_orbit) CALL Find_Closed_Orbit   ! Calculates x_coord_co(1:6)

    Normal_forms: IF(closed_orbit) THEN  !-----------------------!
       !                                                         !
       !                                                         !
       call Get_map_from_NormalForm &                            !
            (ptc_track_debug, Normal_Order_n0, x_coord_co, &     !
            !INTENT:   IN,             IN,             IN        !
            !               1=> Normal_Order_n0  for Linear      !
            Map_Y, Normal_Form_N, A_t_map,   A_t_map_rev)        !
       !     OUT,     OUT,          OUT,       OUT               !
       !                                                         !
    END IF Normal_forms !----------------------------------------!

    ptc_switch = 1 ! set up the "switch"=1 as for "RUN" option
    !              ! In subr. trrun.F, it is input parameter equal to
    !              ! 1 (RUN); 2 (DYNAP fastune); 3 (DYNAP aperture)

    Call ffile_and_segm_for_switch !Int. proc. in this subr.(see below)

    CALL Prepare_Observation_points (max_obs, x_coord_co)
    ! getting parameters  at observation points and
    ! finding the closed orbits at observations for element_by_element tracking

    ! START TRACKING WITH PTC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    change_default: IF(Radiation_PTC) THEN
       DEFAULT=DEFAULT+RADIATION
       IF (Radiation_Quad) STOCH_IN_REC=.TRUE.
       !element_by_element=.FALSE. ! make PTC one-turn tracking
       print *, "################################################################"
       print *, "The PTC parameter DEFAULT before the tracking with the turn-loop"
       call print(default,6)
       CALL Find_Closed_Orbit ! Calculates x_coord_co(1:6)
    END IF change_default

    !============================================================================!
    ! getting inial particle coordinates ========================================!
    !    uses: call comm_para('coord ',nint,ndble,nchar,int_arr,x,char_a,char_l)
    !          returns the value for command parameter "name"

    call ptc_track_ini_conditions

    ! call trinicmd(switch,orbit0,eigen,jmax,z,turns,coords)
    !SUBR trinicmd(switch,orbit0,eigen,jend,z,turns,coords) - in this file
    !----------------------------------------------------------------------*
    ! Purpose: Define initial conditions for all particles to be tracked   *
    ! input:                                                               *
    !   switch (int)  1: run, 2: dynap fastune, 3: dynap aperture          *
    !   orbit0(6) - closed orbit                                           *
    !   x, px, y, py, t, deltap, fx, phix, fy, phiy, ft, phit              *
    !             - raw coordinates from start list                        *
    !   eigen     - Eigenvectors                                           *
    ! output:                                                              *
    !   jend      - number of particles to track                           *
    !   z(6,jend) - Transformed cartesian coordinates incl. c.o.           *
    !   coords      dp(6,0:turns,npart) (only switch > 1) particle coords. *
    !----------------------------------------------------------------------*


    Call Set_initial_particle_ID ! Int.subr. below in this subr.

    ALLOCATE ( last_turn_of_lost_particle(1:j_tot_numb_starting_particles), &
         last_position_of_lost_particle(1:j_tot_numb_starting_particles), &
         last_orbit_of_lost_particle (1:6,1:j_tot_numb_starting_particles))
    last_turn_of_lost_particle=zero
    last_position_of_lost_particle=zero
    last_orbit_of_lost_particle=zero


    Call Particle_Interactions_Ini ! Int.subr. - intialize data for int.subr. Particle_Interactions

    Call Init_info_for_tables       ! Int.subr. below in this subr.

    Call Start_Coord_to_TrackSum_Table
    ! The name of the table "tracksum".
    ! At debug stage it is printed out with the command <write> in input file,
    ! e.g. "write,table=tracksumm, file=track_summary.txt;"

    Call Enter_0st_turn_in_tables
    ! Subroutines "tt_putone" and "tt_puttab"  from the file "trrun.F" are used.
    ! (renamed as "tt_putone_coord" & "tt_puttab_coord"
    ! The name of the tables are "trackone" and "track.obs$$$$.p$$$$"
    ! At debug stage it is printed out with the command <write> in input file,
    ! e.g.  ??????

    !      Call Initialize_kinematics_and_orbit ! removed subr.

    IF (radiation_model1_FZ .OR. Space_Charge ) element_by_element=.TRUE.
    ! make element-by-element tracking (default=.false.),
    ! not one turn tracking
    ! PTC_radiation is not included, it is inherent to PTC, and is done with
    !                                tracking over the total ring

    c_%watch_user=.true.

    debug_print_1: IF (ptc_track_debug) then
       print *; print*,'Start loop over turns: '
       Print *, ' make tracking with element_by_element=', element_by_element
       print *,'The line <MY_RING> consists of <MY_RING%n>=', MY_RING%n, ' nodes'
       print *, "Number of observation points max_obs =", max_obs
    END IF debug_print_1

    loop_over_turns: do i_th_turn=1,turns ! ================loop over turn =======!
       debug_print_2: IF (ptc_track_debug) then                                   !
          print *; print*, 'start ',i_th_turn, '-th turn  '                       !
       end if debug_print_2 !                                                     !
       !                                                                          !
       track_over_total_ring: if ( .NOT.element_by_element) then !++++++++++++!   t
          !                                                                   !   u
          call One_turn_track_with_PTC                                        !   r
          !                                                                   !   n
       else ! +++++track ring element-by-element +++++++++++++++++++++++++++++!   s
          !                                                                   !   !
          call track_beam_elementwise_with_PTC                                !   !
          !                                                                   !   !
          !if (ptc_track_debug) pause ! STOP                                  !   !
          !                                                                   !   !
       end if track_over_total_ring ! ++++++++++++++++++++++++++++++++++++++++!   !
       !                                                                          !

       switch_1: if (ptc_switch .eq. 1) then ! switch=1 => for RUN command ===!   l
          !                                                                   !   o
          debug_print_3: if (ptc_track_debug) then                            !   o
             Print *; Print *, '   ptc_switch= ',ptc_switch                   !   p
             Print *, '  Call Write_tables_after_total_turn'                  !   !
          end if  debug_print_3                                               !   !
          !                                                                   !   o
          call Write_tables_after_total_turn                                  !   v
          !                                                                   !   e
       else ! === ENDIF (switch .ne. 1) => DYNAP (fastune,aperture)===========!   r
          !                                                                   !   !
          debug_print_4: if (ptc_track_debug) then !------------------!       !   !
             Print *; Print *, '   ptc_switch= ',ptc_switch           !       !   !
             Print *, '  Call Save_Coord_for_DYNAP_after_total_turn'  !       !   t
          end if debug_print_4 !--------------------------------------!       !   u
          !                                                                   !   r
          !  ?  Call Save_Coord_for_DYNAP_after_total_turn                    !   n
          !                                                                   !   !
       END IF switch_1!=== END if (switch) ==== RUN or DYNAP =================!   !
       !
       IF (jmax_numb_particl_at_i_th_turn.eq.0) EXIT Loop_over_turns              !
       ! all particles are lost                                                   !
       !
    end do Loop_over_turns !===========loop over turn ============================!


    Call Final_Coord_to_tables ! Complete all tables by one subroutine:
    !'trackone' filled before by by tt_puttab_coord, 'track.obs$$$$.p$$$$' by
    ! tt_putone_coord  and, the summary table 'tracksumm'


    Output_observ_with_PTC: IF(closed_orbit .AND. &
                               (.NOT.element_by_element).AND.(.NOT. Radiation_PTC)) THEN !-!
       debug_print_5: if (ptc_track_debug) then !----------------!                         !
          Print *, 'element_by_element=', element_by_element, &  !                         !
               ' Radiation_PTC=',Radiation_PTC                   !                         !
          Print *, ' Call Observation_with_PTC'                  !                         !
       end if debug_print_5 !------------------------------------!                         !
       Call Observation_with_PTC(max_obs, x_coord_co, Map_Y)                               !
    ENDIF Output_observ_with_PTC !---------------------------------------------------------!


    if (ptc_track_debug) Print *, 'Come to : <! Calculate beam envelope with PTC>'

    Beam_envelope_with_PTC: IF (beam_envelope) THEN !###############################!
       Radiat_PTC: IF ( Radiation_PTC) THEN !====================================!  !
          icase_6: IF (icase_ptc .EQ. 6 .AND.closed_orbit) THEN !--!             !  !
             Call  beam_enevelope_with_PTC                         !             !  !
          ELSE !---------------------------------------------------!             !  !
             Print *, ' Warning !: Option BEAM_ENVELOPE', &        !             !  !
                  ' require option ICASE=6', &                     !             !  !
                  ' and option CLOSED_ORBIT '                      !             !  !
             Print *, '    The code ignores BEAM_ENVELOPE option ' !             !  !
          ENDIF icase_6 !------------------------------------------!             !  !
       ELSE !====================================================================!  !
          Print *, ' Warning !!!: Option BEAM_ENVELOPE require option Radiation' !  !
          Print *, '              The code ignores BEAM_ENVELOPE option '        !  !
       ENDIF Radiat_PTC !========================================================!  !
    ENDIF Beam_envelope_with_PTC !##################################################!

    c_%watch_user=.false.

    Call DeAllocate_local_and_PTC_arrays ! Internal subroutine


    !=============================================================================
  CONTAINS ! Internal procedures for the HOST subr. "ptc_track"
    !                      (a new fortran-90/95 feature):
    !   1) Int.p. is entirely contained within the HOST.
    !   2) It is compiled together with the HOST;
    !   3) It can only be invoked from the host program unit;
    !
    !      No other procedure within the program can access it !!!
    !
    !   4) Int.p. may not CONTAIN futher subprograms;
    !   5) Int.p. automatically has a access to all the host entiteles;
    !   6) Int.p. is able to call other int.procedures of the same HOST ;
    !   7) the HOST does not have access to the local entities (declared
    !      inside int.p.) of any int.subroutine that it contains.
    !   8) The absense of IMPLICIT NONE within int.proc, because one in the
    !      host applies to the intern. proc. as well

    ! Internal procedures are used for this ptc_track_run in order to
    ! use modified subroutines imported from trrun.F with the same names

    !=============================================================================

    SUBROUTINE Check_ptc_universe_and_layout ! before_tracking => Int. proc. (f95)
      ! USE madx_ptc_module, ONLY: universe, index
      ! IMPLICIT NONE ! it already exists in the host

      ! It checks that the following commands has been performed
      ! before calling subr. ptc_track
      ! If not, then return_from_subr_ptc_track =>.TRUE.
      !                                        (accessible from the HOST subr.)

      ! 1) ptc_create_universe:   Needed to set-up PTC

      if (ptc_track_debug) then
         print *
         print *, '<subr. Check_ptc_universe_and_layout (before_tracking)>:'
         print *, 'PTC Universe universe=', universe
         print *, 'PTC Layout      index=', index
      end if
      if(universe.le.0) then
         call fort_warn('return from ptc_track: ',' no universe created')
         ! return
         return_from_subr_ptc_track=.TRUE.
      endif

      ! 2) ptc_create_layout:     creates PTC layout and fills it
      !                                          with current MAD-X sequence

      if(index.le.0) then
         call fort_warn('return from ptc_track: ',' no layout created')
         ! return
         return_from_subr_ptc_track=.TRUE.
      endif
    END  SUBROUTINE Check_ptc_universe_and_layout ! before_tracking =>  Int. proc. (f95)

    !=============================================================================

    SUBROUTINE Values_from_ptc_track_command ! Internal procedure (f95)
      USE madx_ptc_intstate_module, ONLY: getdebug  ! new debug control by PS (from 2006.03.20)
      ! IMPLICIT NONE => in host
      ! local variables
      !real(kind(1d0)) :: maxaper(1:6) move to HOST
      character*12 tol_a, char_a
      integer :: nint,ndble, nchar, int_arr(1),char_l
      data tol_a,char_a / 'maxaper ', ' ' /
      !defaults
      int_arr(:)=1; char_l=1
      ! Values of "ptc_track" command: (PTC thick lens tracking module)
      !                                                                      DEFAULT
      !  icase:         Phase-space (4, 5 or 6)                                (4)
      !  turns:         Number of turns                                        (1)
      !  closed_orbit:  Switch to turn on closed orbit calculation           (false)
      !  deltap:        The relative momentum deviation for
      !                 off-momentum particles (5D case)                ????   (4)
      !
      ! madxdict.h:408:"ptc_track: ptc_track none 0 0 "
      ! "ptc_track: ptc_track none 0 0 "
      ! "icase        = [i, 4], "
      ! "turns        = [i, 1], "
      ! "closed_orbit = [l, false, true], "
      ! "deltap       = [r, 0], "
      ! ! exluded: "coord        = [r, {0,0,0,0,0,0}]; "
      !
      ! added by VKapin
      !
      ! "element_by_element = [l, false, true], " => element_by_element
      ! "radiation = [l, false, true], "  => radiation_PTC
      ! "radiation_model1 = [l, false, true], " => radiation_model1_FZ
      ! "radiation_energy_loss = [l, false, true], "
      ! "radiation_quad = [l, false, true], "
      ! "space_charge = [l, false, true], "
      ! "beam_envelope = [l, false, true], "
      !
      ! "dump     = [l, false, true], " - c-procedure for input
      ! "onetable = [l, false, true], "
      ! "file     = [s, track, track], "
      ! "extension= [s, none, none], "
      ! "maxaper= [r, {0.1, 0.01, 0.1, 0.01, 1., 0.1}], "
      ! "ffile    = [i, 1], "
      ! "norm_out = [l, false, true], "
      ! "norm_no  = [i, 1], "
      ! "debug    = [l, false, true]; " - only for developer (not for User)

      
      icase_ptc    = get_value('ptc_track ','icase ')

      turns        = get_value('ptc_track ','turns ')

      closed_orbit = get_value('ptc_track ','closed_orbit ') .ne. 0

      deltap0      = get_value('ptc_track ','deltap ')

      element_by_element = get_value('ptc_track ', 'element_by_element ') .ne. 0

    IF(max_obs.gt.1.AND.(.NOT.CLOSED_ORBIT).AND.(.NOT. ELEMENT_BY_ELEMENT)) THEN
      Print *, ' '
      Print *, '===================================================================='
      Print *,'To perform tracking with observation points the option'
      Print *,'ELEMENT_BY_ELEMENT must be ON, if the option CLOSED_ORBIT is OFF' 
      call fort_warn(' ELEMENT_BY_ELEMENT',' has been switched ON by the code')
      Print *, '===================================================================='
      Print *, ' '
      ELEMENT_BY_ELEMENT=.TRUE.
    ENDIF

      Radiation_PTC = get_value('ptc_track ','radiation ') .ne. 0

      radiation_model1_FZ = get_value('ptc_track ','radiation_model1 ') .ne. 0
      Radiation_Energy_Loss = get_value('ptc_track ','radiation_energy_loss ') .ne. 0

      Radiation_Quad = get_value('ptc_track ','radiation_quad ') .ne. 0

      Space_Charge = get_value('ptc_track ','space_charge ') .ne. 0

      beam_envelope = get_value('ptc_track ','beam_envelope ') .ne. 0
      ! 'dump ' is done in c-code

      ptc_onetable = get_value('ptc_track ','onetable ') .ne. 0

      ! 'file ' and 'extension ' are done in c-code

      ptc_ffile    = get_value('ptc_track ', 'ffile ')

      norm_out     =  get_value('ptc_track ','norm_out ') .ne. 0

      norm_out_wrong: IF (norm_out .AND. (.NOT.closed_orbit) ) THEN
         Print *, '!!! Warning: The input options  norm_out=', norm_out, ' and closed_orbit=', &
              closed_orbit, ' are uncompatible !!!'
         norm_out=.FALSE.
         Print *, ' The input option is set to norm_out=', norm_out, ' and tracking is continued'
      END IF norm_out_wrong

      Normal_Order_n0  = get_value('ptc_track ','norm_no ')

      !ptc_track_debug  =  get_value('ptc_track ','debug ') .ne. 0  ! before 2006.03.20
      ptc_track_debug  =  getdebug () .ne. 0  ! ptc_track_debug=.T., if debuglevel.ge.1 ,i.e.,
      ! in the madx input: ptc_setswitch, debuglevel=1

      debug_print_1: if (ptc_track_debug) then
         print *,' '
         print*, '====================================================='
         print *, 'Start subr. <ptc_track(max_obs)> in a debugging regime '
         print *, 'Int.Subr <Values_from_ptc_track_command>:'
         print *, 'ptc_track_debug=', ptc_track_debug
      end if debug_print_1
      ! get_value =>
      !  file "madxd.h": getting values from INPUT-user file  !
      ! double get_value(char*, char*)                        !
      !/* C routines called from Fortran and C */             !

      !--- get vector of six coordinate maxapers (both RUN and DYNAP)
      call comm_para(tol_a, nint, ndble, nchar, int_arr, maxaper,char_a, char_l)
      !
      ! INPUT:  name     parameter name
      !        HERE: tol_a <= data tol_a,char_a / 'maxaper ', ' ' / (in this file)
      !
      !      EXAMPLES: RUN,maxaper=double array,turns=integer,ffile=integer;
      !                    maxaper upper limits for the six coordinates
      !                    default=(0.1, 0.01, 0.1, 0.01,1.0, 0.1)
      !              DYNAP,TURNS=real, FASTUNE=logical,LYAPUNOV=real,
      !                    MAXAPER:={..,..,..,..,..,..},ORBIT=logical;
      ! OUTPUT:
      !   (n_int # integers)                 HERE:  nint      USED later: no
      !   (n_double # doubl)                        ndbl                  no
      !   (n_string # string)                       nchar                 no
      !   (int_array := array for integer)          int_arr               no
      !   (double_array := array for double)        maxaper               Yes
      !   (strings :=
      !          one string for all, packed)        char_a                no
      !   (string_lengths:=
      !       length of each string in char)        char_l                no
      !
      ! comm_para
      ! madxd.h: /* C routines called from Fortran and C */
      !          void comm_para(char*, int*, int*, int*, int*, double*, char*, int*);
      ! madxn.c:
      ! void comm_para(char* name, int* n_int, int* n_double, int* n_string,
      !                int* int_array, double* double_array, char* strings,
      !        int* string_lengths)
      !  returns the value for command parameter "name" being either
      !  one or several integers (including logicals),
      !  one or several doubles,
      !  one or several strings (packed in one, with length array)
      !
      !    ATTENTION: no check on sufficient array sizes


      !    call comm_para('coord ',nint,ndble,nchar,int_arr,x,char_a,char_l)

      debug_print_2: if (ptc_track_debug) then
         print *; print *, "<subr. Values_from_ptc_track_command>:"
         print *, "icase_ptc=get_value('ptc_track ','icase ')= ", icase_ptc
         print *, "turns=get_value('ptc_track ','turns ')= ", turns
         print *, "closed_orbit = ", closed_orbit
         print *, &
              "deltap0 = get_value('ptc_track ','deltap ')= ",      deltap0
         print *, " Radiation_PTC =",  Radiation_PTC
         print *, " radiation_model1_FZ =",  radiation_model1_FZ
         print *, " Radiation_Energy_Loss =", Radiation_Energy_Loss
         print *, " Radiation_Quad =", Radiation_Quad
         print *, " Space_Charge =", Space_Charge
         print *, " maxaper =",  maxaper
         print *, " ptc_onetable =",  ptc_onetable
         print *, " ptc_ffile    =",  ptc_ffile
         print *, "  norm_out    =",   norm_out
         print *, "  Normal_Order_n0 =",   Normal_Order_n0
         print *, "   ptc_track_debug =",  ptc_track_debug
      end if debug_print_2

      ! getting particle coordinates with 'coord' keyword
      debug_print_3: if (ptc_track_debug) then
         print *; print*, '----------------------------------', &
              '----------------------------------'
         print *, "before call comm_para('coord ',", &
              "nint,ndble,nchar,int_arr,x,char_a,char_l)"
         print*, "x=", x_input_0_by_ptc_track_command, &
              "IS NOT DEFINED YET"; print *
      end if debug_print_3

      call comm_para('coord ',nint,ndble,nchar,int_arr, &
           x_input_0_by_ptc_track_command,char_a,char_l)

      !  returns the value for command parameter "name"
      debug_print_4: if (ptc_track_debug) then
         print *, "after call comm_para('coord ',", &
              "nint,ndble,nchar,int_arr,x,char_a,char_l)"
         print*, "nint =", nint, "not used later"
         print*, "ndble=", ndble, "not used later"
         print*, "nchar=", nchar, "not used later"
         print*, "int_arr=", int_arr, "not used later"
         print*, "x=", x_input_0_by_ptc_track_command, "MIGHT BE USED LATER"
         print*, "char_a=", char_a, "not used later"
         print*, "char_l=", char_l, "not used later"; print *
      end if debug_print_4
      ! The KeyWord 'coord ' from the "ptc_track"command is not used here
      ! any more...,  instead of it the command "ptc_start" is used

    END SUBROUTINE Values_from_ptc_track_command
    !=============================================================================

    SUBROUTINE Call_my_state_and_update_states
      ! USE  madx_ptc_module, ONLY:  my_state, UPDATE_STATES, default, print
      ! IMPLICIT NONE  => in host

      debug_print_1: if (ptc_track_debug) then
         print *,"before <call my_state(icase,deltap,deltap0)>", &
              "with parameters:  "
         print *, 'icase_ptc=',   icase_ptc, 'deltap=', deltap
         print *, 'deltap0=', deltap0 ; print *
         print*, 'my_state printing:'
      end if debug_print_1
      call my_state(icase_ptc,deltap,deltap0)
      !      IF (Radiation_PTC) DEFAULT=DEFAULT+RADIATION
      Print *, 'Radiation_PTC=', Radiation_PTC

      debug_print_2: if (ptc_track_debug) then
         print *, &
              "after call my_state(icase,deltap,deltap0)"
         print *, '=',icase_ptc, 'deltap=', deltap
         print *, 'deltap0=', deltap0
         print *; print*, '----------------------------------'
         print *, "before CALL UPDATE_STATES"
      end if debug_print_2

      CALL UPDATE_STATES
      debug_print_3: if (ptc_track_debug) then
         print *, "after CALL UPDATE_STATES"
         print *; print*, '----------------------------------'
         print *, "Prining by <call print(default,6)>:"
         call print(default,6)
         print *, "after call print(default,6)"
      end if debug_print_3

    END SUBROUTINE Call_my_state_and_update_states
    !=============================================================================

    !=============================================================================
    SUBROUTINE Find_Closed_Orbit
      ! USE madx_ptc_module, ONLY: dp, zero, find_orbit, my_ring,default
      ! IMPLICIT NONE => in the host
      INTEGER :: i_tmp ! the local counter in the DO-loop
      !====================================================================!
      !   initialize the closed orbit coordinates                          !
      ! x0(:)=zero                                                         !
      x_coord_co(:)=zero                                                   !
      if (ptc_track_debug) THEN !---------------------------!              !
         print *, " x_coord_co(:)=zero = "                  !              !
         CALL write_closed_orbit(icase_ptc,x_coord_co)      !              !
      end if !----------------------------------------------!              !
      !                                                                    !
      if(icase_ptc.eq.5) THEN !------------------------!                   !
         if(mytime) then !----------------------!      !                   !                   
            call Convert_dp_to_dt (deltap, dt)  !      !                   !
         else                                   !      !                   !
            dt=deltap                           !      !                   !
         endif !--------------------------------!      !                   !
         x_coord_co(5)=dt                              !                   !
      ENDIF !------------------------------------------!                   !
      !                                                                    !
      if (ptc_track_debug) then  !------------------------!                !
         print *, " if(icase.eq.5) ,x_coord_co(5)=deltap" !                !
         print *, "  ,x_coord_co(5),deltap=", &           !                !
              x_coord_co(5),deltap                        !                !
      end if !--------------------------------------------!                !
      !                                                                    !
      if(closed_orbit) then !------------------------------------!         !
         call find_orbit(my_ring,x_coord_co,1,default,c_1d_7)    !         !
         print*,"===== ptc_track ============================"   !         !
         CALL write_closed_orbit(icase_ptc,x_coord_co)           !         !
         print*,"============================================"   !         !
      endif  !---------------------------------------------------!         !
      !                                                                    !
      if (ptc_track_debug) print*,"After closed_orbit"; print *;           !
      !                                                                    !
      !END closed_orbit    which is logically (.not.ONEPASS)               !
      !====================================================================!

    END SUBROUTINE Find_Closed_Orbit
    !=============================================================================

    !=============================================================================
    SUBROUTINE DeAllocate_local_and_PTC_arrays
      ! USE madx_ptc_module, ONLY: kill
      ! IMPLICIT NONE => in the host

      if(closed_orbit) then
         call kill (Map_Y)
         call kill (Normal_Form_N)
         CALL kill (A_t_map)
         CALL kill (A_t_map_rev)
      endif

      Deallocate (particle_ID)  ! OK
      Deallocate (jmax_all_turns_numb_part)  ! OK
      Deallocate (x_coord_incl_co)  ! OK

      Deallocate (last_turn_of_lost_particle) !OK
      Deallocate (last_position_of_lost_particle) !OK
      Deallocate (last_orbit_of_lost_particle)   !OK

      IF (.NOT. element_by_element) THEN
         Deallocate ( x_all_incl_co_at0, part_ID_turns)
      ELSE
         Deallocate (x_co_at_all_observ)
      END IF

      Deallocate (elem_number_at_observ, sum_length_at_observ, name_el_at_obsrv)

    END SUBROUTINE DeAllocate_local_and_PTC_arrays
    !=============================================================================

    !=============================================================================

    !=============================================================================
    SUBROUTINE   ffile_and_segm_for_switch
      ! copy a fragment from trrun.f
      ! switch => input prameter to this subroutine trrun
      ! switch  (int) = 1: RUN,
      !                 2: DYNAP fastune,
      !                 3: dynap aperture
      ! IMPLICIT NONE => in the host

      if (ptc_switch .eq. 1)  then  ! like the "switch => =========================!
         !                                           input to the subroutine trrun !
         ! switch  == 1  => RUN                                                    !
         ! ffile = get_value('run ', 'ffile ') is commented                        !
         ptc_ffile = ptc_ffile                                                     !
         ! it is already read by Int.SUBR. <Values_from_ptc_track_command>         !
         !                                                                         !
         ! ffile - integer, local variable, the value is assigned from !           !
         ! "FFILE" option of the "RUN" command =>                      !           !
         ! periodicity of printing coordinates                         !           !
         !                                                                         !
      else   !  2: DYNAP fastune, or 3: dynap aperture ============================!
         ptc_ffile = 1                                                             !
         ! "FFILE" option => periodicity of printing coordinates                   !
      endif  !===== if (switch .eq. 1) ============================================!

      segment_one_table = 0        ! for one_table case
      tot_segm_one_table = turns / ptc_ffile + 1
      ! the output occurs every FFILE turns
      ! turns (integer) input parameter
      ! to this subroutine
      !           = number of turns to track
      ! tot_seg = interger part of [turns / ffile]+1
      !         => number of output segments with particle coordinates

      if (mod(turns, ptc_ffile) .ne. 0) &
           tot_segm_one_table = tot_segm_one_table + 1
      ! correction: the remainder of (turn/ffile), add number of segments

      if (ptc_track_debug) then
         print *; print *, "<subr. ffile_value_for_switch>:"
         print *, " ptc_switch =",  ptc_switch
         print *, " ptc_ffile    =",  ptc_ffile;
         print *, " segment_one_table    =",  segment_one_table;
         print *, " tot_segm_one_table    =",  tot_segm_one_table;
         print *;
      end if

    END SUBROUTINE   ffile_and_segm_for_switch
    !=============================================================================


    !=============================================================================
    SUBROUTINE Set_initial_particle_ID
      ! IMPLICIT NONE => in the host subr.

      integer :: j_th_particle ! local loop counter

      !  x(:)=x(:)+x0(:) ! old ptc_track by FS
      ! print*,"  Initial Coordinates: ", x
      if (ptc_track_debug) then
         print*; Print *, 'Initial particle coordinates of ', &
              jmax_numb_particl_at_i_th_turn, ' particles.'
         Print*, 'No.particle  No_coord,x_incl_co '

         Do j_th_particle=1,jmax_numb_particl_at_i_th_turn
            DO k_th_coord=1,6
               Print*,j_th_particle, k_th_coord, &
                    x_coord_incl_co(k_th_coord,j_th_particle)
            END DO
            CALL write_closed_orbit(icase_ptc,x_coord_co)
         end Do
         Print *, '================================================'
         Print *, ' '
      end if ! debug printing

      ! from trrun.F

      !--- jmax may be reduced by particle loss - keep number in j_tot
      !j_tot = jmax  ! j_tot => initial number of particles
      j_tot_numb_starting_particles =  &
           jmax_numb_particl_at_i_th_turn
      jmax_all_turns_numb_part(0)=j_tot_numb_starting_particles

      ALLOCATE (particle_ID(1:jmax_numb_particl_at_i_th_turn))
      particle_ID=0
      IF (.NOT. element_by_element) THEN
         ALLOCATE (part_ID_turns(0:turns,jmax_numb_particl_at_i_th_turn))
         part_ID_turns=0
      END IF

      !--- set particle id
      do j_th_particle=1, &     ! ==== loop over particles ================!
           jmax_numb_particl_at_i_th_turn                                  !
         !jmax_number_tracking_particles:=                                 !
         !=> the current number of particles                               !
         !                                                                 !
         particle_ID(j_th_particle) = j_th_particle                        !
         ! used inside and as OUTPUT parameter of TRRUN                    !
         IF (.NOT. element_by_element) &                                   !
              part_ID_turns(0,j_th_particle)=particle_ID(j_th_particle)    !
      end do !=============================================================!

    END SUBROUTINE Set_initial_particle_ID

    !=============================================================================

    SUBROUTINE Init_info_for_tables
      implicit none

      !hbu--- init info for tables initial s position is 0
      !hbu initial s position is 0
      spos_current_position = zero      ! = spos in <trrun.F>
      !hbu start of line, element 0
      nlm_current_element_number = 0 ! =  nlm in <trrun.F>
      !hbu
      el_name='start           '

    END SUBROUTINE Init_info_for_tables

    !=============================================================================

    SUBROUTINE Start_Coord_to_TrackSum_Table
      ! IMPLICIT NONE => in the host

      integer :: j_th_particle, k_th_coord ! counter
      ! real(dp) :: tmp_d ! temprorary dble vaiable
      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      REAL (dp) :: X_MAD(6), X_PTC(6)

      ! Summary Table
      !--- enter start coordinates in summary table
      do  j_th_particle = 1, & !=====loop over particles ==== ====================!
           j_tot_numb_starting_particles ! => initial number of particles         !
         !                                                                        !
         doublenum = j_th_particle                                                !
         call double_to_table('tracksumm ', 'number ', doublenum)                 !
         ! tmp_d = 1 <=  turn=1  in the original 2005 trrun.f                     !
         !                                                                        !
         doublenum = zero ! <=  turn=0  for starting particles                    !
         call double_to_table('tracksumm ', 'turn ', doublenum)                   !
         DO k_th_coord = 1, 6 !>>>>> loop over coord. components >>>>>>>>>>>>>!   !
            !tmp_d = z(k_th_coord,j_th_particle) - orbit0(k_th_coord)         !   !
            !z(1:6,1:j_tot) - coordinates                                     !   !
            !orbit0(1:6)    - (dble; 1:6) closed orbit                        !   !
            X_PTC(k_th_coord)=x_coord_incl_co(k_th_coord,j_th_particle) - &   !   !
                 x_coord_co(k_th_coord)                                       !   !
         ENDDO !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!   !
         !                                                                        !
         CALL Coord_PTC_to_MAD(X_PTC,X_MAD)                                       !
         !                                                                        !
         DO k_th_coord = 1, 6 !>>>>> loop over coord. components >>>>>>>>>>>>>!   !
            doublenum  = X_MAD(k_th_coord)                                    !   !
            !                                                                 !   !
           call double_to_table('tracksumm ',vec_names(k_th_coord),doublenum) !   !
            !madxn.c:1385: void                                               !   !
            !double_to_table(char* table,char* name,double* val)              !   !
            ! /* puts val at current position in column                       !   !
            !    with name "name". The table count is increased               !   !
            !    separately with "augment_count" */                           !   !
         enddo !>>>>> END loop over components >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!   !
         Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C, &                   !
              gamma0I,gambet) ! to get "energy" value                             !
                                                  doublenum=energy                !
         call double_to_table('tracksumm ', 'e ', doublenum)                      !
         !                                                                        !
         !hbu add s                                                               !  
                                    doublenum = spos_current_position             !
         call double_to_table('tracksumm ',vec_names(7),doublenum)                !
         ! madxd.h:12:#define double_to_table       double_to_table_              !
         ! madxd.h:150:void double_to_table(char*, char*, double*);               !
         ! ???????????????                                                        !
         call augment_count('tracksumm ')                                         !
         ! madxn.c:1094:void augment_count(char* table)                           !
         ! increase table occ. by 1, fill missing *                               !
      enddo ! ====loop over particle == do  i = 1,j_tot ==========================!

    END SUBROUTINE Start_Coord_to_TrackSum_Table

    !=============================================================================

    SUBROUTINE Enter_0st_turn_in_tables
      ! IMPLICIT NONE => in the host

      ! Local:
      integer :: j_th_particle ! loop counter
      if (ptc_track_debug) then
         print*; Print *, 'Start Subr.<Enter_0st_turn_in_tables>'
         print*, 'ptc_switch=', ptc_switch, '  ptc_ onetable=', ptc_onetable
         ! print *, 'eigen=', eigen, '<call track_pteigen(eigen)> Excluded !!!'
      endif ! debug printing

      !--- enter first turn, and possibly eigen in tables
      if (ptc_switch .eq. 1)  then !  switch=>RUN =====================================!
         if (ptc_onetable)  then !>>>>>>>>>>>> onetable=.TRUE.>>>>>>>>>>>>>>>>>>>>>>!  !
            !VK          call track_pteigen(eigen)                                  !  !
            ! madxn.c:6493:void track_pteigen(double* eigen)                        !  !
            !                                                                       !  !
            !hbu add s, node id and name                                            !  !
            !                                                                       !  !
            !VK          call tt_putone(jmax,    0, tot_segm, segment, part_id, &   !  !
            !VK     &      z, orbit0,spos, nlm, el_name) ! SUBROUTINE in this file  !  !
            !                                                                       !  !
            if (ptc_track_debug) then                                               !  !
               print*; Print *,'before <CALL tt_putone(...)> :'                     !  !
               Print *,'jmax_number_tracking_particles', &                          !  !
                    jmax_numb_particl_at_i_th_turn                                  !  !
               Print *,'tot_segm_one_tablet=', tot_segm_one_table                   !  !
               Print *,'segment_one_table=',   segment_one_table                    !  !
               Print *, 'particle_ID=', (particle_ID(j_th_particle), &              !  !
                    j_th_particle=1,jmax_numb_particl_at_i_th_turn )                !  !
               Print *, 'spos_current_position=', spos_current_position, &          !  !
                    ' nlm_current_element_number=', nlm_current_element_number, &   !  !
                    ' el_name= ', el_name                                           !  !
            end if !  debug printing                                                !  !
            !                                                                       !  !
            !VK Original:  call tt_putone(jmax,    0, tot_segm, segment, part_id, & !  !
            !VK         &      z, orbit0,spos, nlm, el_name)                        !  !
            !                                                                       !  !
            call tt_putone_coord(jmax_numb_particl_at_i_th_turn,    0, &            !  !
                 tot_segm_one_table, segment_one_table, particle_ID, &              !  !
                 x_coord_incl_co, x_coord_co,             &                         !  !
                 spos_current_position, nlm_current_element_number, el_name)        !  !
            !   SUBROUTINE in <trrun.F> file                                        !  !
            !   tt_putone(npart,turn, tot_segm, segment, part_id,                   !  !
            !z, orbit0,spos,ielem,el_name)                                          !  !
            !--- purpolse: enter all particle coordinates in one table  *           !  !
            !INPUT:                                                     *           !  !
            !    npart  (int)           number of particles             *           !  !
            !    turn   (int)           turn number                     *           !  !
            !    tot_segm (int)         total (target) number of entries*           !  !
            !    segment(int)           current segment count           *           !  !
            !    part_id (int array)    particle identifiers            *           !  !
            !    z (double (6,*))       particle orbits                 *           !  !
            !    orbit0 (double array)  reference orbit                 *           !  !
            !------------------------------------------------------------           !  !
            !                                                                       !  !
         else  !>>>>>(.not.onetable )>>>>>> onetable=.TRUE >>>>>>>>>>>>>>>>>>>>>>>>>!  !
            !                                                                       !  !
            do j_th_particle = 1, & !++++++loop over surviving particles +++++++!   !  !
                 jmax_numb_particl_at_i_th_turn                                 !   !  !
               !hbu                                                             !   !  !
               !call tt_puttab(part_id(i),   0,   1, z(1,i), orbit0,spos)       !   !  !
               call tt_puttab_coord(particle_ID(j_th_particle),   0,   1, &     !   !  !
                    x_coord_incl_co(1,j_th_particle), &                         !   !  !
                    x_coord_co, spos_current_position)                          !   !  !
               !SUBR.tt_puttab(     npart,turn,nobs,  orbit, orbit0,spos)       !   !  !
               !(in this file) - purpose: enter particle coordinates in table * !   !  !
               !    input:                                                    * !   !  !
               !    npart  (int)           particle number                    * !   !  !
               !    turn   (int)           turn number                        * !   !  !
               !    nobs   (int)           observation point number           * !   !  !
               !    orbit  (double array)  particle orbit                     * !   !  !
               !    orbit0 (double array)  reference orbit                    * !   !  !
               !--------------------------------------------------------------* !   !  !
            enddo !+++++END of loop over surviving particles +++++++++++++++++++!   !  !
         endif !>> END if (onetable) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!  !
      endif ! == switch =>RUN =========================================================!

    END SUBROUTINE Enter_0st_turn_in_tables

    !=============================================================================


    !=============================================================================

    SUBROUTINE One_turn_track_with_PTC  ! int.subroutine
      ! IMPLICIT NONE => in the host

      ! variables from the HOST subroutine (no need to be declared in int.subr !):

      ! real(dp), intent(INOUT) :: x_coord_incl_co(1:6,1:N_particle_max)
      ! integer,  intent(INOUT) :: jmax_number_tracking_particles
      ! number surviving particles
      ! DBLE,     intent(IN)    :: sum_length
      ! integer,  intent(IN)    :: i_th_turn    :: The number of the current turn
      ! INTEGER,  intent(INOUT) :: particle_ID (1:N_particle_max) ! numbers of surviving particles
      ! Integer,  intent(OUT)   :: last_turn_of_lost_particle (1:N_particle_max)
      ! DBLE,     intent(OUT)   :: last_position_of_lost_particle (1:N_particle_max)
      ! DBLE,     intent(OUT)   :: last_orbit_of_lost_particle (1:6,1:N_particle_max)
      !
      ! Usual tracking steps (pseudocode) in TRRUN (VK)
      !
      ! 1) Check for aperture and kill lost particles
      ! 2) Track all  particles through the element
      ! 3) add length - (not need for the track over turn)
      !
      ! In in this ptc_track there is internal subr.
      !  <PRODUCE_APERTURE_FLAG> and <ANALYSE_APERTURE_FLAG> providing
      !  an extra aperture control, so we add the following extra step:
      ! 4) <ANALYSE_APERTURE_FLAG> and kill lost particles
      !
      !
      ! Tracking over one turn (this internal subroutine) includes only steps 2 and 4,
      !                                           performing them inside particle loop.
      !
      ! !!! the logical structure of check-kill is similar to SUBR. TRCOLL (the file trrun.F.)
      !

      integer :: n_temp, j_last_particle_buffer,jmax_at_loop_start, j_particle

      if (ptc_track_debug) then ! debug printing --------------------------!
         print *; print *, 'Start SUBR.<One_turn_track_with_PTC>'
      end if

      n_temp=1 ! start particle loop with the first not-examined particle

      Exam_all_particles: DO ! ====<=======<========DO Exam_all_particles ====<===<=======<======!
         !(instead of <10 continue> in trrun.f)                                                  !
         !                                                                                       !
         Particle_loop: DO j_particle=n_temp, jmax_numb_particl_at_i_th_turn  !++++++++++++++!   !
            !                                                                                !   !
            jmax_at_loop_start = jmax_numb_particl_at_i_th_turn                              !   ^
            j_last_particle_buffer=j_particle ! remember index value after END DO            !   !
            !                                                                                !   !
            do k_th_coord=1,6 ! extract coords for the current particle -----!               +   ^
               current_x_coord_incl_co(k_th_coord)= &                        !               +   !
                    x_coord_incl_co(k_th_coord,j_particle)                   !               +   !
            end do !---------------------------------------------------------!               +   !
            !                                                                                !   !
            call track(my_ring,current_x_coord_incl_co,1,default)                            !   !
            ! The PTC subroutine " To TRACK the MY_RING for X coordinates                    +   !
            ! over one-turn in the DEFAULT state (citation, p. 25).                          +   !
            ! there is no any other an explicit description in KEK 2002-3 report             +   !
            !                                                                                !   ^
            do k_th_coord=1,6 ! save coordinates for the current particle ---!               +   !
               x_coord_incl_co(k_th_coord,j_particle)=  &                    !               +   !
                    current_x_coord_incl_co(k_th_coord)                      !               +   !
               !                                                             !               +   !
            end do !---------------------------------------------------------!               !   ^
            !                                                                                !   !
            if (ptc_track_debug) then ! debug printing ----------------------------!         !   !
               Print *,'DO j_particle=n_temp, jmax_numb_particl_at_i_th_turn:'     !         !   !
               Print *, 'DO ',j_particle,'=',n_temp,',', &                         !         !   !
                    jmax_numb_particl_at_i_th_turn                                 !         !   ^
            END if ! debug printing -----------------------------------------------!         !   !
            !                                                                                !   !
            if (ptc_track_debug) &                                                           !   !
                 print*,"before printing aperture flag!!!!!!!!",flag_index_ptc_aperture      !   !
            call PRODUCE_APERTURE_FLAG(flag_index_ptc_aperture)                              !   !
            if (ptc_track_debug) then                                                        !   !
               print*,"ready for printing aperture flag!!!!!!!",flag_index_ptc_aperture      !   !
               print*,"real aperture flag: ",c_%aperture_flag                                !   !
            end if                                                                           !   !
            ! Sa_extend_poly.f90:98:  SUBROUTINE PRODUCE_APERTURE_FLAG(I)                    !   !
            !                                                                                !   !
            if (ptc_track_debug) then ! debug printing ----------------------------!         !   ^
               print *, 'PTC: <PRODUCE_APERTURE_FLAG> => flag_index', &            !         +   !
                    flag_index_ptc_aperture                                        !         +   !
               !                                                                   !         +   !
               if(flag_index_ptc_aperture/=0) then !=== print diagnostics =======! !         !   !
                  call ANALYSE_APERTURE_FLAG &                                   ! !         !   !
                       (flag_index_ptc_aperture,why_ptc_aperture)  !             ! !         +   ^
                  !Sa_extend_poly.f90:79:  SUBROUTINE ANALYSE_APERTURE_FLAG(I,R) ! !         +   !
                  Write(6,*) "ptc_track unstable (tracking)-programs continues " ! !         +   !
                  Write(6,*) "why_ptc_aperture:", why_ptc_aperture               ! !         +   !
                  ! See produce aperture flag routine in sd_frame                ! !         +   !
                  ! goto 100 ! EXIT from the turns loop                          ! !         +   ^
               endif !=== print diagnostics =====================================! !         +   !
               !                                                                   !         +   !
            endif ! debug printing ------------------------------------------------!         !   !
            !                                                                                !   !
            if_ptc_track_unstable: IF (flag_index_ptc_aperture/=0) then ! =========!         +   ^
               ! => particle is lost !!(?)                                         !         +   !
               n_temp=j_last_particle_buffer                                       !         +   !
               !                                                                   !         +   !
               CALL kill_ptc_track &                                               !         +   !
                    (n_temp,i_th_turn,zero,jmax_numb_particl_at_i_th_turn, &       !         +   !
                    particle_ID, last_turn_of_lost_particle, &                     !         +   ^
                    last_position_of_lost_particle, last_orbit_of_lost_particle, & !         !   !
                    x_coord_incl_co)                                               !         +   !
               EXIT  Particle_loop !+++++>++++++++>+++++++++>+++++++++++>+++++++>++++++!     +   !
               !                                                                   !   V     +   ^
            END IF if_ptc_track_unstable ! === ptc_track unstable ===>=======>=====!   V     !   !
            !                                                                          V     !   !
         END DO  Particle_loop !+++++++++++++++++++++++++++++++++++++++++++++++++++++++!+++++!   !
         !                                                                             !         !
         !++++++++<+++++++++++<+++++++++++++<+++++++++++++++<+++++++++++<++++++++<+++++!         ^
         !                                                                                       !
         if (j_last_particle_buffer.GE.jmax_at_loop_start ) &                                    !
              EXIT Exam_all_particles !====>=====>======>===========>=========!                  !
         !                                                                    V                  !
      END DO Exam_all_particles ! ==>===DO Exam_all_particles ====>======>====V===>=======>======!
      !                                                                       V
      !<============<===============<===============<===============<=========V

      ! Save particles at exit of the current turn in order to do obs points using PTC maps

      DO j_particle=1, jmax_numb_particl_at_i_th_turn
         do k_th_coord=1,6
            x_all_incl_co_at0(k_th_coord,i_th_turn,j_particle)= &
                 x_coord_incl_co(k_th_coord,          j_particle)
         enddo
         part_ID_turns(i_th_turn,j_particle)=particle_ID(j_particle)
      ENDDO
      jmax_all_turns_numb_part(i_th_turn)=jmax_numb_particl_at_i_th_turn

    END SUBROUTINE One_turn_track_with_PTC

    !=============================================================================

    SUBROUTINE kill_ptc_track &
         (n_particle_kill,i_th_turn_kill,sum_accum_length_kill,j_max_kill, &
         part_ID_kill, last_turn_kill, last_position_of_kill,last_orbit_of_kill, &
         x_coord_incl_co_kill)
      ! IMPLICIT NONE => in the host

      ! the next variables are local
      Integer, intent (IN)    :: n_particle_kill, & ! number of particle to be kill
           i_th_turn_kill
      real(dp), &
           intent (IN)    :: sum_accum_length_kill ! input as 0

      Integer, intent (INOUT) :: j_max_kill, &   ! number surviving particles (here, j=j-1)
           part_ID_kill(*) ! particle number is reodered(1:N_particle_max)

      Integer, intent (OUT) :: last_turn_kill(*) ! last_turn_of_lost_particle (1:N_particle_max)

      real(dp), intent (OUT) :: last_position_of_kill(*), &
           last_orbit_of_kill(1:6,*)

      real(dp), intent (INOUT) :: x_coord_incl_co_kill(1:6,*)

      Integer :: j_th_kill, k_th_kill

      last_turn_kill(part_ID_kill(n_particle_kill))=i_th_turn_kill
      ! Save Number of this turn for n_particle_kill

      last_position_of_kill(part_ID_kill(n_particle_kill))=sum_accum_length_kill
      ! Save Position of n_particle_kill

      DO  k_th_kill=1,6
         last_orbit_of_kill(k_th_kill,part_ID_kill(n_particle_kill)) = &
              x_coord_incl_co_kill(k_th_kill,part_ID_kill(n_particle_kill))
      END DO

      ! Renumbering arrays
      Renumbering_loop: do j_th_kill = n_particle_kill+1, j_max_kill !=======!
         !                                                                   !
         part_ID_kill(j_th_kill-1) = part_ID_kill(j_th_kill)                 !
         !                                                                   !
         do k_th_kill = 1, 6 !--coords loop ----------------------!          !
            x_coord_incl_co_kill(k_th_kill,j_th_kill-1) = &       !          !
                 x_coord_incl_co_kill(k_th_kill,j_th_kill)        !          !
         enddo !--coords loop ------------------------------------!          !
         !                                                                   !
         save_ID_for_all_turns: IF (.NOT.element_by_element) THEN !== !      !
            part_ID_turns(i_th_turn,j_th_kill-1)= &                   !      !
                 part_ID_kill(j_th_kill-1)                            !      !
         END IF save_ID_for_all_turns !===============================!      !
         !                                                                   !
      enddo Renumbering_loop ! ==============================================!

      j_max_kill = j_max_kill - 1

    END SUBROUTINE kill_ptc_track
    !=============================================================================


    !=============================================================================
    SUBROUTINE track_beam_elementwise_with_PTC ! int.subroutine
      ! IMPLICIT NONE => in the host

      ! variables from the HOST subroutine (no need to be declared in int.subr !):

      ! real(dp), intent(INOUT) :: x_coord_incl_co(1:6,1:N_particle_max)
      ! integer,  intent(INOUT) :: jmax_number_tracking_particles
      ! number surviving particles
      ! DBLE,     intent(IN)    :: sum_length
      ! integer,  intent(IN)    :: i_th_turn    :: The number of the current turn
      ! INTEGER,  intent(INOUT) :: particle_ID (1:N_particle_max) ! numbers of surviving particles
      ! Integer,  intent(OUT)   :: last_turn_of_lost_particle (1:N_particle_max)
      ! DBLE,     intent(OUT)   :: last_position_of_lost_particle (1:N_particle_max)
      ! DBLE,     intent(OUT)   :: last_orbit_of_lost_particle (1:6,1:N_particle_max)
      !
      ! Usual tracking steps (pseudocode) in TRRUN (VK)
      !
      ! 1) Check for aperture and kill lost particles
      ! 2) Track all  particles through the element
      ! 3) add length - (not need for the track over turn)
      !
      ! In in this ptc_track there is internal subr.
      !  <PRODUCE_APERTURE_FLAG> and <ANALYSE_APERTURE_FLAG> providing
      !  an extra aperture control, so we add the following extra step:
      ! 4) <ANALYSE_APERTURE_FLAG> and kill lost particles
      !
      !
      ! This subroutine performs a element-by-element particle tracking over one turns
      ! using PTC (this internal subroutine) includes only steps 2,3 and 4,
      !                                           performing them inside particle loop.
      !
      ! !!! the logical structure of check-kill is similar to SUBR. TRCOLL (the file trrun.F.)
      !

      integer :: n_temp, j_last_particle_buffer,jmax_at_loop_start, i_current_elem, &
           iii_c_code, j_th_partic, j_part, number_observation_point
      character(16) :: name_curr_elem
      character(24) :: name_curr_elem_24
      LOGICAL(lp) ::  Entry_not_exit
      REAL(dp) :: length_curr_elem
      real (dp) :: x_coord_co_temp(1:6) ! buffer for the current values of CO

      x_coord_co_temp=zero

      if (ptc_track_debug) then ! debug printing --------------------------!
         print *; print *, 'Start SUBR.< track_beam_elementwise_with_PTC.>'
      end if


      current=>MY_RING%start ! F90 pointer to the CURRENT beamline element is set up
      ! to the starting point of the beamline

      sum_length=zero ! suml=zero the current length along the beamline

      iii_c_code=restart_sequ()  ! c-code restart the beamline sequence

      loop_over_elements: DO i_current_elem=1, MY_RING%n !************ loop over elements ************!
         !                                                                                            *
         ! at element ENTRY ================================================!                         *
         Entry_not_Exit=.TRUE. ! interaction at the entry of the curent     !                         *
         name_curr_elem=current%MAG%name                                    !                         *
         name_curr_elem_24=name_curr_elem                                   !                         *
         length_curr_elem=current%MAG%P%ld                                  !                         *
         if (ptc_track_debug) THEN !+++debug print++++++!                   !                         l
            Print *, &                                  !                   !                         o
                 'i_current_elem=', i_current_elem, &   !                   !                         o
                 ' Entry_not_exit=', Entry_not_exit, &  !                   !                         o
                 ' name_curr_elem=', name_curr_elem, &  !                   !                         p
                 ' sum_length=',sum_length, &           !                   !                         !
                 ' length_curr_elem=',length_curr_elem  !                   !                         !
            !Print *,'name=',current%MAG%name, &        !                   !                         e
            !        'l=',current%MAG%P%ld              !                   !                         l
         endif !+++end debug print++++++++++++++++++++++!                   !                         e
         Call Paricle_Interactions &                                        !                         m
              (i_current_elem, name_curr_elem, Entry_not_Exit,&             !                         e
              sum_length,length_curr_elem)                                  !                         n
         !                                                                  !                         t
         ! at element ENTRY ================================================!                         s
         !                                                                                            !
         n_temp=1 ! start particle loop with the first not-examined particle                          *
         !                                                                                            *
         Exam_all_particles: DO ! ===<=======<========DO Exam_all_particles ====<============<======! *
            !                                                                                       ! *
            Particle_loop: DO j_th_partic=n_temp, jmax_numb_particl_at_i_th_turn !+++++++++++++++!  ! !
               !                                                                                 !  ! !
               jmax_at_loop_start = jmax_numb_particl_at_i_th_turn                               !  ^ !
               j_last_particle_buffer=j_th_partic   ! remember index value after END DO          !  ! !
               !                                                                                 !  ! l
               do k_th_coord=1,6 ! extract coords for the current particle -----!                +  ^ o
                  current_x_coord_incl_co(k_th_coord)= &                        !                +  ! o
                       x_coord_incl_co(k_th_coord,j_th_partic)                  !                +  ! p
               end do !---------------------------------------------------------!                +  ! !
               !                                                                                 +  ^ !
               call track(my_ring,current_x_coord_incl_co, &                                     !  ! !
                    i_current_elem,i_current_elem+1,default)                                     !  ! o
               ! The PTC subroutine " To TRACK the MY_RING for X coordinates                     +  ! v
               ! throughout the element in the DEFAULT state (citation, p. 25).                  +  ! e
               !Print *, 'x=', current_x_coord_incl_co                                           +  ! r
               !                                                                                 !  ! !
               do k_th_coord=1,6 ! save coordinates for the current particle ---!                +  ! !
                  x_coord_incl_co(k_th_coord,j_th_partic)=  &                   !                +  ! !
                       current_x_coord_incl_co(k_th_coord)                      !                +  ! !
                  !                                                             !                +  ! !
               end do !---------------------------------------------------------!                !  ^ !
               !                                                                                 !  ! !
               !if (ptc_track_debug) &                                                           !  ! !
               !     print*,"before printing aperture flag!!!!!!!!",flag_index_ptc_aperture      !  ! !
               call PRODUCE_APERTURE_FLAG(flag_index_ptc_aperture)                               !  ! e
               ! Sa_extend_poly.f90:98:  SUBROUTINE PRODUCE_APERTURE_FLAG(I)                     !  ! l
               !if (ptc_track_debug) then                                                        !  ! e
               !     print*,"ready for printing aperture flag!!!!!!!",flag_index_ptc_aperture    !  ! m
               !     print*,"real aperture flag: ",c_%aperture_flag                              !  ! e
               !end if                                                                           !  ! n
               !                                                                                 !  ! t
               if (ptc_track_debug) then ! debug printing ----------------------------!          !  ^ s
                  !print *, 'PTC: <PRODUCE_APERTURE_FLAG> => flag_index', &           !          +  ! !
                  !                               flag_index_ptc_aperture             !          +  ! !
                  !                                                                   !          +  ! !
                  if(flag_index_ptc_aperture/=0) then !=== print diagnostics ======!  !          !  ! !
                     call ANALYSE_APERTURE_FLAG &                                  !  !          !  ! !
                          (flag_index_ptc_aperture,why_ptc_aperture)               !  !          +  ^ !
                     !Sa_extend_poly.f90:79:  SUBROUTINE ANALYSE_APERTURE_FLAG(I,R)!  !          +  ! !
                     Write(6,*) "ptc_track unstable (tracking)-programs continues "!  !          +  ! !
                     Write(6,*) "why_ptc_aperture:", why_ptc_aperture              !  !          +  ! !
                     ! See produce aperture flag routine in sd_frame               !  !          +  ! !
                     ! goto 100 ! EXIT from the turns loop                         !  !          +  ^ !
                  endif !=== print diagnostics ====================================!  !          +  ! !
                  !                                                                   !          +  ! !
               endif ! debug printing ------------------------------------------------!          !  ! !
               !                                                                                 !  ! !
               if_ptc_track_unstable: IF (flag_index_ptc_aperture/=0) then ! ========!           +  ^ !
                  ! => particle is lost !!(?)                                        !           +  ! !
                  n_temp=j_last_particle_buffer                                      !           +  ! l
                  !                                                                  !           +  ! o
                  CALL kill_ptc_track &                                              !           +  ! o
                       (n_temp,i_th_turn,sum_length, &                               !           +  ! p
                       jmax_numb_particl_at_i_th_turn, &                             !           +  ! !
                       particle_ID, last_turn_of_lost_particle, &                    !           +  ^ !
                       last_position_of_lost_particle, last_orbit_of_lost_particle, &!           !  ! o
                       x_coord_incl_co)                                              !           +  ! v
                  EXIT  Particle_loop !+++++>++++++++>+++++++++>+++++++++++>+++++++>++++++!      +  ! e
                  !                                                                  !    v      +  ^ r
               END IF if_ptc_track_unstable ! === ptc_track unstable ===>=======>====!    V      !  ! !
               !                                                                          v      !  ! !
            END DO  Particle_loop !+++++++++++++++++++++++++++++++++++++++++++++++++++++++v++++++!  ! e
            !                                                                             V         ! l
            !++++++++<+++++++++++<+++++++++++++<+++++++++++++++<+++++++++++<++++++++<+++++!         ^ e
            !                                                                                       ! m
            if (j_last_particle_buffer.GE.jmax_at_loop_start ) &                                    ! e
                 EXIT Exam_all_particles !====>=====>======!                                        ! n
            !                                                                                       ! t
         END DO Exam_all_particles ! ==>===DO Exam_all_particles ====>======>====V===>=======>======! s
         !                                                                                            !
         !<============<===============<===============<===============<=========V                    !
         !                                                                                            !
         sum_length=sum_length+current%MAG%P%ld                                                       !
         !                                                                                            !
         ! at element EXITY ================================================!                         *
         Entry_not_Exit=.FALSE. ! interaction at the entry of the curent    !                         *
         name_curr_elem=current%MAG%name                                    !                         *
         length_curr_elem=current%MAG%P%ld                                  !                         *
         if (ptc_track_debug) THEN !+++debug print+++++++!                  !                         *
            Print *, &                                   !                  !                         *
                 'i_current_elem=', i_current_elem, &    !                  !                         *
                 '  Entry_not_exit=', Entry_not_exit, &  !                  !                         *
                 '  name_curr_elem=', name_curr_elem, &  !                  !                         *
                 ' sum_length=',sum_length, &            !                  !                         !
                 '  length_curr_elem=',length_curr_elem  !                  !                         *
            !Print *,'name=',current%MAG%name, &         !                  !                         *
            !        'l=',current%MAG%P%ld               !                  !                         *
         endif !+++end debug print+++++++++++++++++++++++!                  !                         *
         Call Paricle_Interactions &                                        !                         *
              (i_current_elem, name_curr_elem, Entry_not_Exit,&             !                         *
              sum_length,length_curr_elem)                                  !                         *
         !                                                                  !                         *
         ! at element ENTRY ================================================!                         *
         !                                                                                            !
         number_observation_point=node_value('obs_point ')                                            !
         !                                                                                            !
         OBSERVATION_POINTS : IF (number_observation_point .GT. 1) THEN !##############!              !
            !                                                                          #              !
            if (mod(i_th_turn, ptc_ffile).EQ.0) then !:::: if(turn/ffile)=int ::::!    #              !
               ! mod() => the remainder of (turn/ffile),                          !    #              !
               ! i_th_turn := the current turn number (not total "turns")         !    #              !
               ! ffile := the periodicity of printing coordinate                  !    #              !
               !          ("FFILE"is the  option of the "RUN" command             !    #              !
               !                                                                  !    #              !
               if (i_th_turn .eq. turns)  last_table_line_out = .true.            !    #              *
               ! set the last turns = OK                                          !    #              *
               !                                                                  !    #              *
               extract_CO_at_observ_point: DO k_th_coord=1,6                      !    #              *
                  x_coord_co_temp(k_th_coord)= &                                  !    #              *
                       x_co_at_all_observ(k_th_coord,number_observation_point)    !    #              *
               ENDDO extract_CO_at_observ_point                                   !    #              *
                  if (ptc_track_debug) THEN !+++debug print+++++++!               !    #              *
                    Print *,'obs.No=',number_observation_point    !               !    #              *
                    Print *,'CO=', x_coord_co_temp                !               !    #              *
                    Print *,'x_coord_incl_co=', x_coord_incl_co   !               !    #              *                     
                  endif !+++++++++++++++++++++++++++++++++++++++++!               !    #              *
               !                                                                  !    #              *
               if (ptc_onetable) then !>>>>> onetable=.TRUE.>>>>>>>>>>>>>>>>!     !    #              *
                  !                                                         !     !    #              *
                  !hbu                                                      !     !    #              *
                  spos_current_position = sum_length ! sum_length           !     !    #              *
                  ! (dble) initial s position ???                           !     !    #              *
                  ! sum (dble) (in subr. ttmap) Accumulated length          !     !    #              *
                  !hbu spos added                                           !     !    #              *
                  !                                                         !     !    #              *
                  call tt_putone_coord(jmax_numb_particl_at_i_th_turn, &    !     !    #              *
                       i_th_turn, tot_segm_one_table,segment_one_table, &   !     !    #              *
                       particle_ID, x_coord_incl_co, &                      !     !    #              *
                       x_coord_co_temp, & ! x_coord_co, &                   !     !    #              *
                       spos_current_position, i_current_elem, &             !     !    #              *
                       name_curr_elem_24)                                   !     !    #              *
                  ! purpose: enter all particle coordinates in one table    !     !    #              *
                  !                                                         !     !    #              *
                  !                                                         !     !    #              *
               else !>>>>> onetable=.FALSE.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!     !    #              *
                  !                                                         !     !    #              *
                  !++loop over surviving particles +++!                     !     !    #              *
                  do j_part = 1, &                    !++===============!   !     !    #              *
                       jmax_numb_particl_at_i_th_turn                   !   !     !    #              *
                     !hbu                                               !   !     !    #              *
                     call tt_puttab_coord &                             !   !     !    #              *
                          (particle_ID(j_part),i_th_turn, &             !   !     !    #              *
                          number_observation_point, &                   !   !     !    #              *
                          x_coord_incl_co(1,j_part), &                  !   !     !    #              *
                          x_coord_co_temp, & ! x_coord_co, &            !   !     !    #              *
                          sum_length)                                   !   !     !    #              *
                     !SUBR.tt_puttab &                                  !   !     !    #              *
                     !        (npart, turn,nobs,  orbit, orbit0,spos)   !   !     !    #              *
                     !(in this file) - purpose: enter particle          !   !     !    #              *
                     !                          coordinates in table    !   !     !    #              *
                     !INPUT:                                            !   !     !    #              *
                     !  npart  (int)           particle number          !   !     !    #              *
                     !  turn   (int)           turn number              !   !     !    #              *
                     !  nobs   (int)           observation point number !   !     !    #              *
                     !  orbit  (double array)  particle orbit           !   !     !    #              *
                     !  orbit0 (double array)  reference orbit          !   !     !    #              *
                  enddo !+++++END of loop over surviving particles +++++!   !     !    #              *
                  !                                                         !     !    #              *
               endif ! >>>> ENDIF (onetable) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!     !    #              *
               !                                                                  !    #              *
            endif ! END :::: if(turn/ffile)=integer :^^:::::::::::::::::::::::::::!    #              *
            !                                                                          #              *
         ENDIF OBSERVATION_POINTS ! ###################################################!              *
         !                                                                                            *
         !                                                                                            *
         iii_c_code=advance_node()                                                                    !
         current=>current%next                                                                        !
         !                                                                                            !
      ENDDO loop_over_elements !************ loop over elements **************************************!

      jmax_all_turns_numb_part(i_th_turn)=jmax_numb_particl_at_i_th_turn

    END SUBROUTINE track_beam_elementwise_with_PTC
    !=============================================================================

    !=============================================================================
    SUBROUTINE Write_tables_after_total_turn
      ! segment of trrun.F
      ! IMPLICIT NONE => in the host

      INTEGER :: j_th_part ! local
      if (ptc_track_debug) then
         Print *, ' Start SUBR. <Write_tables_after_total_turn>'
         Print *, ' jmax=', jmax_numb_particl_at_i_th_turn
         Print *, ' spos_current_position=', spos_current_position
         Print *, ' nlm_current_element_number=', nlm_current_element_number
         Print *, ' el_name=', el_name
      end if

      if (mod(i_th_turn, ptc_ffile).EQ.0) then !:::: if(turn/ffile)=int ::::!     !       e
         ! mod() => the remainder of (turn/ffile),                          !     !       r
         ! turn  := the current turn number (not total "turns")             !     R       !
         ! ffile := the periodicity of printing coordinate                  !     U       !
         !          ("FFILE"is the  option of the "RUN" command             !     N       t
         !                                                                  !     !       u
         if (i_th_turn .eq. turns)  last_table_line_out = .true.            !     !       r
         ! set the last turns = OK                                          !     !       n
         !                                                                  !     !       s
         if (ptc_onetable) then !>>>>> onetable=.TRUE.>>>>>>>>>>>>>>>>>!    !     !       !
            !                                                          !    !     !       !
            spos_current_position = zero ! sum_length                  !    !     !       !
            ! (dble) initial s position ???                            !    !     !       !
            ! sum (dble) (in subr. ttmap) Accumulated length           !    !     !       !
            !hbu spos added                                            !    !     !       !
            !                                                          !    !     !       !
            call tt_putone_coord(jmax_numb_particl_at_i_th_turn, &     !    !     !       !
                 i_th_turn, tot_segm_one_table, segment_one_table, &   !    !     !       !
                 particle_ID, x_coord_incl_co, x_coord_co, &           !    !     !       !
                 spos_current_position,nlm_current_element_number, &   !    !     !       !
                 el_name)                                              !    !     R       !
            ! purpose: enter all particle coordinates in one table     !    !     U       l
            ! detailed comments see above                              !    !     N       o
            !                                                          !    !     !       o
         else !>>>>> onetable=.FALSE.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !     !       p
            !                                                          !    !     !       !
            !++loop over surviving particles +++!                      !    !     !       !
            do j_th_part = 1, &                 !++==============!     !    !     !       !
                 jmax_numb_particl_at_i_th_turn                  !     !    !     !       !
               !hbu                                              !     !    !     !       o
               call tt_puttab_coord &                            !     !    !     !       v
                    (particle_ID(j_th_part),i_th_turn,  1, &     !     !    !     !       e
                    x_coord_incl_co(1,j_th_part), &              !     !    !     !       r
                    x_coord_co, spos_current_position)           !     !    !     !       !
               !SUBR.tt_puttab &                                 !     !    !     !       !
               !        (npart, turn,nobs,  orbit, orbit0,spos)  !     !    !     !       !
               !(in this file) - purpose: enter particle         !     !    !     !       !
               !                          coordinates in table   !     !    !     !       t
               !INPUT:                                           !     !    !     !       u
               !  npart  (int)           particle number         !     !    !     !       r
               !  turn   (int)           turn number             !     !    !     !       n
               !  nobs   (int)           observation point number!     !    !     !       s
               !  orbit  (double array)  particle orbit          !     !    !     !       !
               !  orbit0 (double array)  reference orbit         !     !    !     !       !
            enddo !+++++END of loop over surviving particles ++++!     !    !     !       !
            !                                                          !    !     !       !
         endif ! >>>> ENDIF (onetable) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !     !       l
         !                                                                  !     !       o
      endif ! END :::: if(turn/ffile)=integer ::::::::::::::::::::::::::::::!     !       o
      !                                                                                   p
    END SUBROUTINE Write_tables_after_total_turn
    !=============================================================================

    !=============================================================================
    SUBROUTINE Final_Coord_to_tables ! Complete all tables by one subroutine:
      !'trackone' filled before by by tt_puttab_coord, 'track.obs$$$$.p$$$$' by
      ! tt_putone_coord  and, the summary table 'tracksumm'

      ! IMPLICIT NONE => in the host
      ! Local variables
      !real(dp) :: tmp_dble ! temprorary dble vaiable => global
      INTEGER :: j_part_tmp, turn_final, i_coord
      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      REAL (dp) :: X_MAD(6), X_PTC(6)
      sum_length =zero
      turn_final=min(turns,i_th_turn) ! the FORTRAN bad feature:
      ! if DO-loop ends, turn > turns
      debug_Final_Coord: if (ptc_track_debug) then
         Print *, ' Start SUBR. <<Final_Coord_to_tables>>:'
         Print *, ' jmax=', jmax_numb_particl_at_i_th_turn
         Print *, 'Total number of turns: turns=', turns
         Print *, 'The current turn number: turn=', i_th_turn
         Print *, 'The actual final turn number: turn_final=', turn_final
         Print *, 'The dimensionality of the task: icase_ptc=',icase_ptc
      ENDIF  debug_Final_Coord

      !do i         =1,jmax !++loop over surviving particles ++++!
      do j_part_tmp=1,jmax_numb_particl_at_i_th_turn             !
         !                                                       !
         !last_turn(part_id(i)) = min(turns, turn)               !
         last_turn_of_lost_particle(particle_ID(j_part_tmp))= &  !
              turn_final                                         !
         !last_pos(part_id(i)) = sum                             !
         last_position_of_lost_particle &                        !
              (particle_ID(j_part_tmp))= sum_length              !
         ! remember last turn and position of particles          !
         !                                                       !
         !do j = 1, 6 !>>> loop over coord. components >>>>>>>!  !
          do i_coord = 1, 6! icase_ptc                        !  !
            !last_orbit(j,part_id(i)) = z(j,i)                !  !
            last_orbit_of_lost_particle(i_coord,    &         !  !
                 particle_ID(j_part_tmp))=  &                 !  !
                 x_coord_incl_co(i_coord, j_part_tmp)         !  !
            ! last orbit => to finalize tables                !  !
         enddo ! END loop over coord. components >>>>>>>>>>>> !  !
         !                                                       !
      enddo ! ++ END loop over surviving particles ++++++++++++++!

      !turn = min(turn, turns)  => turn_final (see above)

      !--- enter last turn in tables if not done already
      !if (.not. last_out)  then ! ==== enter last turn in tables if not done =======!
      if (.not. last_table_line_out)  then                                           !
         ! flag to avoid double entry of last line                                   !
         !                                                                           !
         !if (switch .eq. 1)  then ! @@@@ switch=>RUN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@! !
         IF (ptc_switch.EQ. 1) THEN                                                ! !
            !                                                                      ! !
            !if (onetable)  then !>>>>>>> onetable=.TRUE.>>>>>>>>>>>>>>>>>>>>>>>!  ! !
            IF (ptc_onetable) THEN                    !>>> a single file >>>>>> !  ! !
               !                                                                !  ! !
               !hbu                                                             !  ! !
               !spos=sum ! (dble) initial s position ???                        !  ! !
               spos_current_position=sum_length                                 !  ! !
               ! sum (dble) (in subr. ttmap)   Accumulated length               !  ! !
               !                                                                !  ! !
               if (ptc_track_debug) then  ! ----- debug printing ---------- !   !  ! !
                  Print *, ' jmax=', jmax_numb_particl_at_i_th_turn         !   !  ! !
                  Print *, 'The final turn number: turn_final=', turn_final !   !  ! !
                  Print *, 'tot_segm_one_table=', tot_segm_one_table        !   !  ! !
                  Print *, 'segment_one_table=', segment_one_table          !   !  ! !
                  Print *, 'particle_ID=', particle_ID                      !   !  ! !
                  Print *, ' spos_current_position=', spos_current_position !   !  ! !
                  ! check later and compare with trrun                      !   !  ! !
                  Print *, ' nlm_current_element_number=', &                !   !  ! !
                       nlm_current_element_number; Print *, ' el_name=',&   !   !  ! !
                       el_name                                              !   !  ! !
               end if  ! ----- debug printing ----------------------------- !   !  ! !
               !                                                                !  ! !
               CALL tt_putone_coord(jmax_numb_particl_at_i_th_turn, &           !  ! !
                    turn_final, tot_segm_one_table, segment_one_table, &        !  ! !
                    particle_ID, x_coord_incl_co, x_coord_co, &                 !  ! !
                    spos_current_position, &                                    !  ! !
                    nlm_current_element_number, el_name)                        !  ! !
               !                                                                !  ! !
               !hbu get current node name                                       !  ! !
               ! call element_name(el_name,len(el_name))                        !  ! !
               !hbu spos added                                                  !  ! !
               ! call tt_putone(jmax, turn, tot_segm, segment, part_id,  &      !  ! !
               !                z, orbit0,spos,nlm,el_name)                     !  ! !
               ! SUBROUTINE in this file                                        !  ! !
               !      tt_putone(npart,turn, tot_segm, segment, part_id,         !  ! !
               !                z, orbit0,spos,ielem,el_name)                   !  ! !
               !--- purpose: enter all particle coordinates in one table   *    !  ! !
               !    input:                                                 *    !  ! !
               !    npart  (int)           number of particles             *    !  ! !
               !    turn   (int)           turn number                     *    !  ! !
               !    tot_segm (int)         total (target) number of entries*    !  ! !
               !    segment(int)           current segment count           *    !  ! !
               !    part_id (int array)    particle identifiers            *    !  ! !
               !    z (double (6,*))       particle orbits                 *    !  ! !
               !    orbit0 (double array)  reference orbit                 *    !  ! !
               !-----------------------------------------------------------*    !  ! !
               !                                                                !  ! !
            else ! >>>>>>>>>>>>>>>>>>> onetable=.TRUE. >>>>>>>>>>>>>>>>>>>>>>>>>!  ! !
               ! write all particle coordinates :            >>>>>>>>>>>>>>>>>  !  ! !
               ! one file per particle and observation point >>>>>>>>>>>>>>>>>  !  ! !
               !                                                                !  ! !
               spos_current_position=sum_length                                 !  ! !
               !                                                                !  ! !
               !do i = 1, jmax !#### loop over surv. particles ###############! !  ! !
               do j_part_tmp=1,jmax_numb_particl_at_i_th_turn                 ! !  ! !
                  !                                                           ! !  ! !
                  call tt_puttab_coord &                                      ! !  ! !
                       (particle_ID(j_part_tmp),turn_final,  1, &             ! !  ! !
                       x_coord_incl_co(1,j_part_tmp), &                       ! !  ! !
                       x_coord_co, spos_current_position)                     ! !  ! !
               enddo ! END loop over surv. particles #########################! !  ! !
               !                                                                !  ! !
            endif ! >>> END if (onetable) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!  ! !
            !                                                                      ! !
         endif ! END if (switch .eq. 1)  @@@@ switch=>RUN @@@@@@@@@@@@@@@@@@@@@@@@@! !
         !                                                                           !
      endif ! ====END:  enter last turn in tables if not done =======================!

      !    Summary table ===================!
      !--- enter last turn in summary table !---------------------------------------!
      !+      do  i          = 1,j_tot !#### loop over all started particles #######!
      do  j_part_tmp = 1,j_tot_numb_starting_particles  !###########################!
         !tmp_d = i ! convert INTEGER to DBLE                                       !
         doublenum  = j_part_tmp                                                    !
         call  double_to_table('tracksumm ', 'number ', doublenum)                  !
         !tmp_d = last_turn(i)                                                      !
         doublenum=last_turn_of_lost_particle(j_part_tmp)                           !
         call double_to_table('tracksumm ', 'turn ', doublenum)                     !
         ! call double_to_table('tracksumm ', 'turn ', tmp_d)                       !
         !                                                                          !
         !do j       = 1, 6 !>>>> loop over coord. components >>>>>>>>>>>>>>>>>!    !
         DO i_coord = 1, 6                                                     !    !
            !tmp_d =  last_orbit(j,i) - orbit0(j)                              !    !
            X_PTC(i_coord)=last_orbit_of_lost_particle(i_coord,j_part_tmp)- &  !    !
                 x_coord_co(i_coord)                                           !    !
         ENDDO !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !
         !                                                                          !
         CALL Coord_PTC_to_MAD(X_PTC,X_MAD)                                         !
         !                                                                          !
         DO i_coord = 1, 6 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !
            doublenum = X_MAD(i_coord)                                         !    !
            call double_to_table('tracksumm ', vec_names(i_coord), doublenum)  !    !
            !call double_to_table('tracksumm ', vec_names(j), tmp_d)           !    !
            !                                                                  !    !
         enddo ! END loop over coord. components >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !
         !hbu                                                                       !
         !  spos                  = last_pos(i)                                     !
         spos_current_position = last_position_of_lost_particle(j_part_tmp)         !
         !hbu                                                                       !
                                                 doublenum = spos_current_position  !
         call double_to_table('tracksumm ',vec_names(7),doublenum)                  !
         !                                                                          !
         ! to get "energy" value                                                    !
         Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)        !
         !                                                                          !
                                                   doublenum= energy                !
         call double_to_table('tracksumm ', 'e ',  doublenum)                       !
         !                                                                          !
         call augment_count('tracksumm ')                                           !
      enddo !#### loop over all started particles ##################################!

    END SUBROUTINE Final_Coord_to_tables
    !=============================================================================

    !==============================================================================
    SUBROUTINE Particle_Interactions_Ini
      !USE ptc_track_run_common, ONLY: &
      !            Energy_rest_MeV, Energy_total_MeV
      ! IMPLICIT NONE => in the host

      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet

      Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)
      ! real(dp) ,optional,INTENT(OUT)::MASS,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet

      Energy_rest_MeV=MASS_GeV*c_1D3
      Energy_total_MeV=ENERGY*c_1D3

    END SUBROUTINE Particle_Interactions_Ini
    !==============================================================================

    !==============================================================================
    SUBROUTINE Paricle_Interactions &
         (i_current_elem, name_curr_elem, Entry_not_Exit,&
         sum_length,length_curr_elem)
      ! This intenal subroutine of the host subr. <ptc_track_run>
      ! serves as service routine to call any particle interactions
      ! (radiations, space charge and etc.) in a between of the beamline elements)
      !USE madx_ptc_track_run_common, ONLY: &
      !    Energy_rest_MeV, Energy_total_MeV
      ! USE madx_ptc_module, ONLY: dp
      ! IMPLICIT NONE => in the host

      INTEGER, INTENT(IN)       :: i_current_elem
      CHARACTER(16), INTENT(IN) :: name_curr_elem !      current%MAG%    name
      ! TYPE:fibre   element Character(16)
      LOGICAL(lp), INTENT(IN)       ::  Entry_not_exit
      REAL(dp), INTENT(IN)      ::  sum_length, length_curr_elem

      Real (dp) :: B0_dipole, Quadr_k, TiltD_dipole, rad_curv_m, &
           SQRT_X2_Y2, x_comp, y_comp
      INTEGER :: i_elem_type, j_partic, ieave, iquasto


      real (dp), dimension(1:2) :: d_loss

      !Print *,MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      ! fv9.madX:0.938, 100    99.061

      Radiation_by_FZ_code: IF (radiation_model1_FZ) THEN
         !if (ptc_track_debug) Print *, 'Radiation_FZ=',Radiation_FZ, &
         !                             ' Radiation_Energy_Loss=', Radiation_Energy_Loss, &
         !                             ' Radiation_Quad=', Radiation_Quad

         ieave=0;   IF (Radiation_Energy_Loss) ieave=1;
         iquasto=0; IF (Radiation_Quad)      iquasto=1;

         B0_dipole = current%MAG%P%B0
         TiltD_dipole = current%MAG%P%TILTD
         IF (current%MAG%P%NMUL >= 2) THEN
            Quadr_K=current%MAG%BN(2)
         ELSE
            Quadr_k=zero
         ENDIF
         if (ptc_track_debug) Print *, 'B0_dipole=', B0_dipole, &
              'TiltD_dipole=', TiltD_dipole, &
              '  Quadr_k=',Quadr_k
         IF (B0_dipole.EQ.zero .AND.Quadr_k .EQ.zero ) i_elem_type=0
         IF (B0_dipole.NE.zero .AND.Quadr_k .EQ.zero ) i_elem_type=1
         IF (B0_dipole.EQ.zero .AND.Quadr_k .NE.zero ) i_elem_type=2
         IF (B0_dipole.NE.zero .AND.Quadr_k .NE.zero ) i_elem_type=3
         if (ptc_track_debug) Print *,'i_elem_type=',i_elem_type

         IF (i_elem_type .EQ. 0) RETURN

         IF ( (i_elem_type .EQ. 1) .OR. &
              ( (i_elem_type.EQ.3).AND.(.NOT.Radiation_Quad) ) ) THEN
            rad_curv_m = 1D0/B0_dipole
            Call photon (i_elem_type, rad_curv_m, length_curr_elem, &
                 Energy_total_MeV, Energy_rest_MeV, ieave,iquasto, d_loss(1), d_loss(2))
            DO j_partic=1, jmax_numb_particl_at_i_th_turn
               IF (Entry_not_Exit) THEN
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(1)
               ELSE
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(2)
               ENDIF
            ENDDO
         END IF

         IF (( i_elem_type .EQ. 2) .AND. Radiation_Quad) THEN
            if (ptc_track_debug) Print *,'jmax_numb_particl_at_i_th_turn=', &
                 jmax_numb_particl_at_i_th_turn
            DO j_partic=1, jmax_numb_particl_at_i_th_turn
               if (ptc_track_debug) Print *, 'j_partic=',j_partic
               SQRT_X2_Y2=SQRT(x_coord_incl_co(1,j_partic)*x_coord_incl_co(1,j_partic)+ &
                    x_coord_incl_co(3,j_partic)*x_coord_incl_co(3,j_partic))
               !IF (SQRT_X2_Y2 .EQ. zero) EXIT
               rad_curv_m = 1D0/Quadr_k/SQRT_X2_Y2

               Call photon (i_elem_type, rad_curv_m, length_curr_elem, &
                    Energy_total_MeV, Energy_rest_MeV, ieave,iquasto,d_loss(1), d_loss(2))
               IF (Entry_not_Exit) THEN
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(1)
               ELSE
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(2)
               ENDIF
            ENDDO
         END IF

         IF (( i_elem_type .EQ. 3) .AND. Radiation_Quad) THEN
            DO j_partic=1, jmax_numb_particl_at_i_th_turn

               x_comp=one/B0_dipole*cos(TiltD_dipole)+ &
                    one/Quadr_k/x_coord_incl_co(1,j_partic)

               y_comp=one/B0_dipole*sin(TiltD_dipole)+ &
                    one/Quadr_k/x_coord_incl_co(3,j_partic)

               rad_curv_m = SQRT(x_comp*x_comp+ y_comp*y_comp)

               Call photon (i_elem_type, rad_curv_m, length_curr_elem, &
                    Energy_total_MeV, Energy_rest_MeV,ieave,iquasto, d_loss(1), d_loss(2))
               IF (Entry_not_Exit) THEN
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(1)
               ELSE
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(2)
               ENDIF
            ENDDO
         END IF

         ! Call photon (i_Mag_type, rad_curv_m, dl_length, Energy_beam_MeV, Rest_mass_MeV,
         !                                                           ieave,iquasto,d1, d2)
         ! Radiation_Energy_Loss=TRUE => ieave=1; Radiation_Energy_Loss=FALSE => ieave=0;
         ! Radiation_Quad     =TRUE => iquasto=1; Radiation_Quad     =FALSE => iquasto=0
      END IF Radiation_by_FZ_code


      Space_Charge_Calculation: IF ( Space_Charge ) THEN
         ! Call Space_Charge
         if (ptc_track_debug) Print *, '  i_current_elem=',i_current_elem, &
              '  name_curr_elem=', name_curr_elem, &
              '  sum_length=', sum_length

      ENDIF Space_Charge_Calculation

      ! print * , ' SUBROUTINE Paricle_Interactions is under construction => RETURN'
      RETURN
    END SUBROUTINE Paricle_Interactions
    !==============================================================================

    !==============================================================================
    SUBROUTINE Prepare_Observation_points (max_obs, x_coord_co_at_START)
      ! getting parameters  at observation points and
      ! finding the closed orbits at observations for element_by_element tracking

      ! USE  madx_ptc_module, ONLY: dp, my_ring, track, default
      ! IMPLICIT NONE => in the host

      !type(fibre), POINTER :: current ! advance node in PTC - global in <PTC_TRACK_RUN>
      EXTERNAL             :: element_name

      integer, intent (IN) :: max_obs ! the maximum number of observation points >=1
      ! one point at the end (beginning) plus
      ! the points given in input file by the command
      ! "ptc_observe,place=mark";
      REAL (dp), INTENT( IN) :: x_coord_co_at_START(1:6) ! C.O. at START of ring
      real (dp) :: x_coord_co_temp(1:6) ! buffer for the current values of CO

      INTEGER :: iii_c_code, i_ring_element, number_obs, i_coord

      REAL(dp) :: sum_length_s  ! = sum  <in trrun.F> !!!! LOCAL in Observations

      REAL(dp) :: length_current_element_f90, length_current_element_c ! from two databases

      character (16) ::      name_16
      integer, parameter  :: name_length_16 = 16 ! 16

      debug_printing_1: if (ptc_track_debug) then
         Print *, 'Start subr. <Prepare_Observation_points> max_obs=', max_obs
         Print *, 'x_coord_co_at_START=', x_coord_co_at_START
      endif debug_printing_1
      ! Define the number of element at observation point

      allocate(elem_number_at_observ(1:max_obs)); elem_number_at_observ=1
      allocate(sum_length_at_observ(1:max_obs)); sum_length_at_observ(:)=zero
      allocate(name_el_at_obsrv(1:max_obs));
      el_by_el_1: IF (element_by_element) THEN
         allocate(x_co_at_all_observ(1:6,1:max_obs)); x_co_at_all_observ(:,:)=zero
      ENDIF el_by_el_1

      elem_number_at_observ(1)=1 ! The first point at begining (end) of ring

      Sum_length_S=zero

      number_obs=0

      x_coord_co_temp(:)=x_coord_co_at_START(:)

      iii_c_code=restart_sequ()  ! c-code restart the beamline sequence
      current=>MY_RING%start     ! F90 pointer to the CURRENT beamline element is set up
      ! to the starting point of the beamline

      element_number: DO i_ring_element=1, my_ring%n
         number_obs=node_value('obs_point ')
         IF (i_ring_element.EQ.1) number_obs=1 ! node_value gives 0 for 1st (?)

         length_current_element_f90=current%MAG%P%ld
         length_current_element_c=node_value('l ')
         Sum_length_S=Sum_length_S+length_current_element_f90

         find_CO_for_el_by_el: IF (element_by_element) THEN !===!
            IF (closed_orbit) THEN !-------------------------!  !
              Call track(my_ring,x_coord_co_temp, &          !  !
              i_ring_element, i_ring_element+1, default )    !  !
            ELSE                                             !  !
              x_coord_co_temp(:)=zero                        !  !
            ENDIF !------------------------------------------!  !
         ENDIF find_CO_for_el_by_el !===========================!

         IF(number_obs.GT.0) THEN
            elem_number_at_observ(number_obs)= i_ring_element
            sum_length_at_observ(number_obs) = Sum_length_S
            call element_name(name_16,name_length_16)
            name_el_at_obsrv(number_obs) = name_16

            save_CO_for_el_by_el: IF (element_by_element) THEN
               DO i_coord=1,icase_ptc
                  x_co_at_all_observ(i_coord,number_obs)=x_coord_co_temp(i_coord)
               ENDDO
            ENDIF save_CO_for_el_by_el

         ENDIF
         if (ptc_track_debug) then
            Print *, 'i_el ', i_ring_element,' num_obs=', number_obs
            Print *, ' l_c_code=',length_current_element_c,' l_f90=', &
                 length_current_element_f90, &
                 ' name_c=', name_16, ' &_f90=', current%MAG%name
            Print *, 'x_coord_co_temp=', x_coord_co_temp
                 
         endif


         iii_c_code=advance_node() ! c-code go to the next node
         current=>current%next     ! f90-code bring to the next mode
      ENDDO element_number 

      debug_print_2: if (ptc_track_debug) then
         Print *, 'elem_number_at_observ(i_obs)= ', elem_number_at_observ
         Print *, 'sum_length_at_observ(i_obs)= ', sum_length_at_observ
         print_CO_for_el_by_el: IF (element_by_element) THEN
            Print *, ' === CO at observations x_co_at_all_observ(i_coord,number_obs)='
            DO number_obs=1, max_obs
               Print *, 'number_obs=', number_obs, &
                    ' elem_numb=', elem_number_at_observ(number_obs), &
                    ' name=', name_el_at_obsrv(number_obs)
               Print *, (x_co_at_all_observ(i_coord,number_obs),i_coord=1,6)
            END DO
         ENDIF  print_CO_for_el_by_el
      endif debug_print_2

    END SUBROUTINE Prepare_Observation_points
    !==============================================================================


    !==============================================================================
    SUBROUTINE Observation_with_PTC(max_obs, x_coord_co_at_START, Map_Y_obs)
      !USE  madx_ptc_module, ONLY: my_ring, kill, default, REAL_8, alloc,  & !BERZ,
      !                            assignment(=), track, lnv, zero,  & !print, sub,
      !                            PRODUCE_APERTURE_FLAG, ANALYSE_APERTURE_FLAG, &
      !                            fibre, daprint, operator(*), operator(.sub.)
      ! IMPLICIT NONE => in the host

      integer, intent (IN) :: max_obs ! the maximum number of observation points >=1
      ! one point at the end (beginning) plus
      ! the points given in input file by the command
      ! "ptc_observe,place=mark";
      REAL (dp),     INTENT( IN) :: x_coord_co_at_START(1:6) ! => x0(1:6) in ptc_track at START
      TYPE (real_8), intent (INOUT) :: Map_Y_obs(6) !  y => Map_Y_obs - local name

      TYPE(damap) :: Map_damap

      INTEGER ::  iii_c_code, i_obs_point, i_coord, &
           i_from, i_till, i_dummy, i_Unit, i_turn_tmp, j_part_tmp, ielem
      !INTEGER :: nda, NormOrder, npara, mynd2
      INTEGER :: flag_index,why(9)

      INTEGER, ALLOCATABLE :: Temp_particle_ID(:)

      ! type(fibre), POINTER :: current ! advance node in PTC - global in <PTC_TRACK_RUN>

      REAL(dp) :: spos

      !TYPE (real_8) :: Map_Y_obs(6) !  y => Map_Y_obs - local name
      real(dp) :: X_co_temp (6) ! (36) !(6) !(lnv)

      real(dp), allocatable ::  X_co_observe(:,:) ! (i_coord, i_obs_point)
      integer, allocatable :: J(:) ! for extracting CO from the map
      REAL(dp),  allocatable :: Temp_X_incl_co_at_obs(:,:)
      Real(dp):: X_lnv_START(lnv), X_lnv_OBSRV(lnv)

      segment_one_table=0 ! for Printing to tables
      last_table_line_out = .false.

      spos=zero

      DO i_coord=1,6
         X_lnv_START(i_coord)=zero
         X_lnv_OBSRV(i_coord)=zero
      ENDDO


      ! Define maps between the beginning of the ring and the current observation point

      ! nda=0; NormOrder=1 npara=0;mynd2=0;
      ! call init(default,NormOrder,nda,BERZ,mynd2,npara) - already known
      if (ptc_track_debug) then
         write(31,*) 'Map_Y_obs IN: '; call daprint(Map_Y_obs,31);
      endif
      Map_Y_obs=npara
      Map_Y_obs=x_coord_co_at_START  ! Y=X

      if (ptc_track_debug) then
         Print *, ' x_coord_co_at_START=', x_coord_co_at_START
         write(32,*) 'Map_Y_obs=x_coord_co_at_START: '; call daprint(Map_Y_obs,32);
      endif

      Allocate (X_co_observe(1:6,1:max_obs), &
           Temp_particle_ID(1:j_tot_numb_starting_particles))
      X_co_observe(:,:)=zero; Temp_particle_ID=0

      assign_co_at_start_ring: DO i_coord=1,icase_ptc
         X_co_observe(i_coord,1)=x_coord_co_at_START(i_coord)
      END DO assign_co_at_start_ring

      CALL Alloc(Map_damap)

      iii_c_code=restart_sequ()  ! c-code restart the beamline sequence
      current=>MY_RING%start     ! F90 pointer to the CURRENT beamline element is set up

      obs_point_loop: DO i_obs_point=1, max_obs-1

         i_from=elem_number_at_observ(i_obs_point)
         i_till=elem_number_at_observ(i_obs_point+1)

         if (ptc_track_debug) then
            Print *, 'i_obs_point=', i_obs_point, ' name_f90=', current%MAG%name
            Print *, 'X_co_observe(:,i_obs_point) =', &
                 (X_co_observe(i_coord,i_obs_point), i_coord=1,6)
            Print *, 'Track from i_from=', i_from, 'i_till =',i_till
         end if

         call track(my_ring,Map_Y_obs,i_from,i_till,default)
         !call track(my_ring,Map_Y_obs,1,default)


         Map_damap=Map_Y_obs
         if (ptc_track_debug) then
            i_Unit=20+i_obs_point
            write(i_Unit,*) 'i_Unit=', i_Unit; Call daprint(Map_Y_obs,i_Unit);
         endif
         call PRODUCE_APERTURE_FLAG(flag_index)
         if(flag_index/=0) then
            call ANALYSE_APERTURE_FLAG(flag_index,why)
            Write(6,*) "Get_map_from_NormalForm is unstable (map production)-programs continues "
            Write(6,*) why ! See produce aperture flag routine in sd_frame
            CALL kill(Map_Y_obs) !CALL kill(y)
            ! c_%watch_user=.false.
            return
         endif

         ielem=elem_number_at_observ(i_obs_point+1)

         ! Call Extract_CO_from_Map (icase_PTC, Map_Y_obs, X_co_temp) - subr. removed
         allocate(J(icase_ptc)); J(:)=0
         X_co_temp(:)=zero
         DO i_coord=1, icase_ptc
            X_co_temp(i_coord) = (Map_Y_obs(i_coord)%T.sub.J) ! take line with all zero-order
            if (ptc_track_debug) Print *,'CO extracted: i_coord=', i_coord, &
                 ' X_co_temp(i_coord)=', X_co_temp(i_coord)
         ENDDO
         deallocate(J)
         Save_co_in_X_co_observe: DO i_coord=1, icase_ptc !icase_PTC
            X_co_observe(i_coord,i_obs_point+1)= X_co_temp(i_coord)
         END DO Save_co_in_X_co_observe

         ! transform the coordinates from the ring origin to observation point

         Local_loop_over_turns: DO i_turn_tmp=1,turns-1 ! don't provide (the turns+1)!
            if (ptc_track_debug) then
               Print *, 'Local_loop_over_turns i_turn_tmp= ', i_turn_tmp
               Print *, 'For the current jmax=', jmax_all_turns_numb_part(i_turn_tmp)
            endif
            ALLOCATE( Temp_X_incl_co_at_obs(1:6,1:jmax_all_turns_numb_part(i_turn_tmp)))
            Temp_X_incl_co_at_obs=zero

            Local_loop_for_particles: DO j_part_tmp=1,jmax_all_turns_numb_part(i_turn_tmp)
               if (ptc_track_debug) Print *, '   the current particle j_part_tmp= ', j_part_tmp



               !X_lnv_OBSRV=Map_Y_obs*X_lnv_START
               !             X_lnv_OBSRV=Map_Y_obs%T*X_lnv_START
               !Error: Operands of binary numeric operator '*' at (1) are TYPE(taylor)/real(kind(1d0))
               Read_START_coord: DO i_coord=1,icase_ptc
                  X_lnv_START(i_coord)= x_all_incl_co_at0(i_coord, i_turn_tmp, j_part_tmp) - &
                       x_coord_co_at_START(i_coord) ! X_out=M(x0,x)*x_in, where X=x0+x !
               ENDDO Read_START_coord
               ! Transform from START to observation point
               X_lnv_OBSRV=Map_damap*X_lnv_START
               !X_lnv_OBSRV=Map_Y_obs*X_lnv_START
               !             X_lnv_OBSRV=Map_Y_obs%T*X_lnv_START
               !Error: Operands of binary numeric operator '*' at (1) are TYPE(taylor)/real(kind(1d0))

               Loop_coord: DO i_coord=1,icase_ptc
                  !X_lnv_OBSRV(i_coord)=Map_Y_obs(i_coord)%T*X_lnv_START(i_coord)
                  !X_lnv_OBSRV(i_coord)=Map_Y_obs(i_coord)*X_lnv_START(i_coord)
                  !X_lnv_OBSRV(i_coord)=Map_damap*X_lnv_START(i_coord)
                  if (ptc_track_debug) Print *, 'i_coord=',i_coord, &
                       'X_lnv_START/OBSRV=', X_lnv_START(i_coord), X_lnv_OBSRV(i_coord)
                  Temp_X_incl_co_at_obs(i_coord,j_part_tmp)= X_lnv_OBSRV(i_coord)
               ENDDO Loop_coord

               Temp_particle_ID(j_part_tmp)=part_ID_turns(i_turn_tmp,j_part_tmp)

            ENDDO  Local_loop_for_particles

            ! Print coordinates at observation point
            Debug_Print_4: IF (PTC_track_debug) THEN
               Print *, 'Before writing to tables i_obs=', i_obs_point
               Print *, ' X_co_temp=',  X_co_temp
            ENDIF  Debug_Print_4

            if (mod(i_turn_tmp, ptc_ffile).EQ.0) then !:::: if(turn/ffile)=int :::!     !       e
               ! mod() => the remainder of (turn/ffile),                          !     R       r
               ! ffile := the periodicity of printing coordinate                  !     U       !
               !          ("FFILE"is the  option of the "RUN" command             !     N       t
               !                                                                  !     !       u
               if (i_turn_tmp .eq. turns)  last_table_line_out = .true.           !     !       r
               ! set the last turns = OK                                          !     !       n
               !                                                                  !     R       s
               if (ptc_onetable) then !>>>>> onetable=.TRUE.>>>>>>>>>>>>>>>>>!    !     U       !
                  !                                                          !    !     N       !
                  spos = sum_length_at_observ(i_obs_point+1)                 !    !     !       !
                  ! current_position = Accumulated length                    !    !     !       !
                  !                                                          !    !     !       !
                  call tt_putone_coord( &                                    !    !     !       !
                       jmax_all_turns_numb_part(i_turn_tmp), i_turn_tmp, &   !    !     !       !
                       tot_segm_one_table, segment_one_table, &              !    !     !       !
                       Temp_particle_ID, Temp_X_incl_co_at_obs, &            !    !     !       !
                       X_co_temp, spos, &                                    !    !     !       !
                       elem_number_at_observ(i_obs_point+1), &               !    !     !       !
                       name_el_at_obsrv(i_obs_point+1))                      !    !     R       !
                  ! purpose: enter all particle coordinates in one table     !    !     U       l
                  ! detailed comments see above                              !    !     N       o
                  !                                                          !    !     !       o
               else !>>>>> onetable=.FALSE.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !     !       p
                  !                                                          !    !     !       !
                  !++loop over surviving particles +++!                      !    !     !       !
                  do j_part_tmp = 1, &                !++=================!  !    !     !       !
                       jmax_all_turns_numb_part(i_turn_tmp)               !  !    !     !       !
                     !                                                    !  !    !     !       o
                     call tt_puttab_coord &                               !  !    !     !       v
                          (Temp_particle_ID(j_part_tmp),i_turn_tmp, &     !  !    !     !       e
                          i_obs_point+1,                             &    !  !    !     !       r
                          Temp_X_incl_co_at_obs(1,j_part_tmp),       &    !  !    !     !       !
                          X_co_temp, sum_length_at_observ(i_obs_point+1)) !  !    !     !       !
                     !SUBR.tt_puttab &                                    !  !    !     !       !
                     !        (npart, turn,nobs,  orbit, orbit0,spos)     !  !    !     !       !
                  enddo !+++++END of loop over surviving particles +++++++!  !    !     !       t
                  !                                                          !    !     !       u
               endif ! >>>> ENDIF (onetable) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !     !       r
               !                                                                  !     !       n
            endif ! END :::: if(turn/ffile)=integer ::::::::::::::::::::::::::::::!     !       s
            !                                                                                   !
            deALLOCATE( Temp_X_incl_co_at_obs)                                                  !
         ENDDO Local_loop_over_turns !==========================================================!

         Jump_to_next_observation_point: DO i_dummy=i_from,i_till-1
            iii_c_code=advance_node() ! c-code go to the next node
            current=>current%next     ! f90-code bring to the next mode
         END DO Jump_to_next_observation_point
         !
      END DO obs_point_loop !####################################################################!

      DEALLOCATE( X_co_observe,  Temp_particle_ID)

      CALL kill(Map_damap)
      !Call kill(Map_Y_obs)
    END SUBROUTINE Observation_with_PTC
    !==============================================================================

    !==============================================================================
    SUBROUTINE tt_putone_coord (npart,turn,tot_segm,segment,part_id,z,orbit0,&
         &spos,ielem,el_name)
      ! copied from TRRUN.F(tt_putone) and renamed to tt_putone_coord
      !hbu added spos, ielem, el_name

      ! USE  madx_ptc_module, ONLY: dp,lnv,zero, operator(*), assignment(=)
      ! IMPLICIT NONE => in the host
      !----------------------------------------------------------------------*
      !--- purpose: enter all particle coordinates in one table              *
      !    input:                                                            *
      !    npart  (int)           number of particles                        *
      !    turn   (int)           turn number                                *
      !    tot_segm (int)         total (target) number of entries           *
      !    segment(int)           current segment count                      *
      !    part_id (int array)    particle identifiers                       *
      !    z (double (6,*))       particle orbits                            *
      !    orbit0 (double array)  reference orbit                            *
      !----------------------------------------------------------------------*
      include 'name_len.fi'
      integer i,j,npart,turn,tot_segm,segment,part_id(*),length
      REAL(dp) :: z(6,*),orbit0(6),tmp,tt, ss
      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      !vvk
      REAL(dp) :: tmp_coord_array(lnv), tmp_norm_array(lnv), tmp_norm

      REAL (dp) :: X_MAD(6), X_PTC(6)

      !hbu was *36 allow longer info
      character*80 table_putone,comment
      !hbu
      integer :: ielem
      !hbu name of element
      character*(name_len) :: el_name
      !hbu
      REAL(dp) :: spos
      !hbu
      character*4 vec_names(7)
      !hbu
      data vec_names / 'x', 'px', 'y', 'py', 't', 'pt','s' / ! MAD order
      !data vec_names / 'x', 'px', 'y', 'py', 'pt', 't','s' / ! PTC has a reverse order for pt and t
      data table_putone / 'trackone' /

      !vk
      tmp_coord_array=zero; tmp_norm_array=zero;

      !hbu
      length = len(comment)
      segment = segment + 1
      !hbu
      !write(comment, '(''!segment'',4i8,1X,A)')   &
      write(comment, '(''#segment'',4i8,1X,A)')   &
           &segment,tot_segm,npart,ielem,el_name
      call comment_to_table(table_putone, comment, length)
      Call GET_ONE(MASS_GeV,energy,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet) ! to get "energy" value
      do i = 1, npart
                                                     doublenum=turn
         call double_to_table(table_putone, 'turn ', doublenum)
                                                  doublenum=energy
         call double_to_table(table_putone, 'e ', doublenum)
                                                       doublenum = part_id(i)
         call double_to_table(table_putone, 'number ', doublenum)

         do j = 1, 6
            tmp=zero
            !IF (j.LE.icase_ptc) tmp = z(j,i) - orbit0(j)
            tmp = z(j,i) - orbit0(j)
            tmp_coord_array(j)=tmp ! make array of coordinates
            X_PTC(j)=tmp
         end do
         CALL Coord_PTC_to_MAD(X_PTC,X_MAD) ! convert coordinates
         do j = 1, 6
            doublenum=X_MAD(j)
            IF( (.NOT.closed_orbit) .OR. (.NOT.NORM_OUT)) THEN
               call double_to_table(table_putone, vec_names(j), doublenum)
            END IF
         enddo

         ! make normalization, if necessary
         IF( closed_orbit .AND. NORM_OUT) THEN
            tmp_norm_array=A_t_map_rev*tmp_coord_array

            DO j = 1, 6
               X_PTC(j)=zero
               IF (j.LE.icase_ptc) X_PTC(j)=tmp_norm_array(j)
            ENDDO
            CALL Coord_PTC_to_MAD(X_PTC,X_MAD) ! convert coordinates
            DO j = 1, 6
               !tmp_norm=zero
               doublenum=X_MAD(j)
               call double_to_table(table_putone, vec_names(j), doublenum)
            END DO
         END IF

         !hbu spos
                                                        doublenum=spos
         call double_to_table(table_putone,vec_names(7),doublenum)
         call augment_count(table_putone)
      enddo
    END SUBROUTINE tt_putone_coord
    !==============================================================================

    !==============================================================================
    SUBROUTINE tt_puttab_coord (npart,turn,nobs,orbit,orbit0,spos)
      ! copied from TRRUN.F( tt_puttab) and renamed to tt_puttab_coord
      !
      ! USE  madx_ptc_module, ONLY: dp, lnv, zero, operator(*), assignment(=)
      ! IMPLICIT NONE => in the host
      !----------------------------------------------------------------------*
      !--- purpose: enter particle coordinates in table                      *
      !    input:                                                            *
      !    npart  (int)           particle number                            *
      !    turn   (int)           turn number                                *
      !    nobs   (int)           observation point number                   *
      !    orbit  (double array)  particle orbit                             *
      !    orbit0 (double array)  reference orbit                            *
      !----------------------------------------------------------------------*

      !vvk
      REAL(dp) :: tmp_coord_array(lnv), tmp_norm_array(lnv) !, tmp_norm

      integer  :: npart,turn,j,nobs
      REAL(dp) :: orbit(6),orbit0(6),tmp !,tt,tn

      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      REAL (dp) :: X_MAD(6), X_PTC(6)

      character*36 table_puttab
      !hbu
      REAL(dp) :: spos
      !hbu
      character*4 vec_names(7)
      !hbu
      data vec_names  / 'x', 'px', 'y', 'py', 't', 'pt','s' / ! MAD order
      !data vec_names / 'x', 'px', 'y', 'py', 'pt', 't','s' / ! PTC has a reverse order for pt and t
      data table_puttab / 'track.obs$$$$.p$$$$' /

      tmp_coord_array=zero; tmp_norm_array=zero

      !tt = turn; !tn = npart
      write(table_puttab(10:13), '(i4.4)') nobs    ! Write in the table head :
      write(table_puttab(16:19), '(i4.4)') npart   ! "@NAME ... "TRACK.OBS0001.P0005"

      Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet) ! to get "energy" value
                                               doublenum=energy
      call double_to_table(table_puttab, 'e ', doublenum)
                                                    doublenum=npart
      call double_to_table(table_puttab, 'number ', doublenum) ! number of the current particle
                                                  doublenum=turn
      call double_to_table(table_puttab, 'turn ', doublenum)   ! the number of the current turn
      do j = 1, 6
         tmp=zero
         !IF (j.LE.icase_ptc) tmp = orbit(j) - orbit0(j)
         tmp = orbit(j) - orbit0(j)
         tmp_coord_array(j)=tmp ! make array of coordinates
         X_PTC(j)=tmp
      end do

      IF( closed_orbit .AND. NORM_OUT) THEN
         tmp_norm_array=A_t_map_rev*tmp_coord_array
         DO j = 1, 6
            X_PTC(j)=zero
            IF (j.LE.icase_ptc) X_PTC(j)=tmp_norm_array(j)
         ENDDO

      ENDIF
      CALL Coord_PTC_to_MAD(X_PTC,X_MAD) ! convert coordinates (canonical or normal)
      do j = 1, 6
         IF (closed_orbit .AND. NORM_OUT) THEN
            !tmp_norm=zero
            doublenum=X_MAD(j)
            call double_to_table(table_puttab, vec_names(j), doublenum)
         ELSE
            !tmp=zero
            doublenum = X_MAD(j) ! orbit(j) - orbit0(j)
            call double_to_table(table_puttab, vec_names(j), doublenum)
         END IF
      enddo
      !hbu spos
                                                     doublenum=spos
      call double_to_table(table_puttab,vec_names(7),doublenum)
      call augment_count(table_puttab)
    END SUBROUTINE tt_puttab_coord
    !==============================================================================


    !==============================================================================
    SUBROUTINE beam_enevelope_with_PTC
      !use madx_ptc_module, ONLY: BEAMENVELOPE, ENV_8, REAL_8, &
      !                           Find_Envelope, my_ring, default, track, &
      !                           assignment(=), operator(+)
      ! IMPLICIT NONE => in the host

      TYPE (BEAMENVELOPE) :: ENV
      TYPE (ENV_8) :: YS(6)
      TYPE (REAL_8) :: A1(6), Y(6)
      ! REAL (dp)     :: x_coord_co(6)

      Print *, ' Subr. beam_enevelope_with_PTC: '
      Print*, 'x_coord_co=', x_coord_co
      Print *, ' before Call FIND_ENVELOPE(my_ring,YS,A1,x_coord_co, 1, default) '

      Y=x_coord_co ! new line at 2005018

      Call FIND_ENVELOPE(my_ring,YS,A1,x_coord_co, 1, default)
      Print *, ' after Call FIND_ENVELOPE(my_ring,YS,A1,x_coord_co, 1, default) '

      ENV=YS
      Print *, ' after ENV=YS '
      YS=Y
      Print *, ' after ENV=YS '
      YS=ENV%SIJ0
      Print *, ' after YS=ENV%SIJ0 '
      CALL TRACK (my_ring, YS, 1,  +default)
      Print *, ' after CALL TRACK (my_ring, YS, 1,  +default) '
      ENV%SIJ0=YS

      Print *, ' The equilibrium emittances are : ', ENV%EMITTANCE

    END SUBROUTINE beam_enevelope_with_PTC
    !==============================================================================


    !=============================================================================

    SUBROUTINE ptc_track_ini_conditions
      ! this subroutine is a truncated and commented copy of
      !k      subroutine trinicmd(switch,orbit0,eigen,jend,z,turns,coords)
      !
      !  First all code lines had been commented by !k, and then
      !  long fortran names in  a fashion of f90 has been introduced
      !  instead of hort names of f77.
      !
      ! USE madx_ptc_module, ONLY: dp, lnv, zero, twopi, &
      !                           operator(*), assignment(=)
      ! IMPLICIT NONE => in the host subr.

      !----------------------------------------------------------------------*
      ! Purpose:                                                             *
      !   Define initial conditions for all particles to be tracked          *
      ! input:                                                               *
      !   switch (int)  1: run, 2: dynap fastune, 3: dynap aperture          *
      !   orbit0(6) - closed orbit                                           *
      !   x, px, y, py, t, deltap, fx, phix, fy, phiy, ft, phit              *
      !             - raw coordinates from start list                        *
      !   eigen     - Eigenvectors                                           *
      ! output:                                                              *
      !   jend      - number of particles to track                           *
      !   z(6,jend) - Transformed cartesian coordinates incl. c.o.           *
      !   coords      dp(6,0:turns,npart) (only switch > 1) particle coords. *
      !----------------------------------------------------------------------*

      !k      include 'bb.fi': --------------------------------------!
      !k      integer bbd_loc,bbd_cnt,bbd_flag,bbd_pos,bbd_max       !
      !k      parameter(bbd_max=200)                                 !
      !k      real(kind(1d0)) bb_kick                                !
      !k      common/bbi/bbd_loc(bbd_max),bbd_cnt,bbd_flag,bbd_pos   !
      !k      common/bbr/bb_kick(2,bbd_max !-------------------------!
      !k
      !k      logical zgiv,zngiv
      !k      integer j,jend,k,kp,kq,next_start,itype(23),switch,turns
      INTEGER ::  j_particle_line_counter,kq,kp
      INTEGER ::  next_start ! int. function
      !k      real(kind(1d0)) phi,track(12),zstart(12),twopi,z(6,1000),zn(6),  &
      !k     &ex,ey,et,orbit0(6),eigen(6,6),x,px,y,py,t,deltae,fx,phix,fy,phiy, &
      !k     &ft,phit,get_value,get_variable,zero,deltax,coords(6,0:turns,*)
      !real(kind(1d0)) :: twoPi, Ex_horz_emi_m, Ey_vert_emi_m, Et_long_emi_m
      REAL(dp)  :: x_input,px_input,y_input,py_input,t_input,deltae_input, &
           fx_input,phix_input,fy_input,phiy_input,ft_input,phit_input, &
           phi_n

      ! real(kind(1d0)) :: get_value,get_variable ! dble c-functions in the HOST routine
      ! local temprorary variables
      INTEGER :: k_th_coord,  j_th_particle, & ! local
           itype_non_zero_flag(12) ! not itype(23)
      ! =1 for non-zero coordinates, =0 for zero
      real(dp) :: track_temp(12), z_start_sum(6)

      real(dp),allocatable :: x_coord_incl_co_temp(:,:) ! copy of array

      Real (dp) :: Z_norm_temp(lnv), Z_stdt_temp(lnv) ! Z_stdt=A_t_map*Z_norm

      REAL (dp) :: X_MAD(6), X_PTC(6)

      logical(lp) :: zgiv_exist, zngiv_exist ! existence of non-zero input for action-angle

      !k      parameter(zero=0d0)
      character*120 msg(2) ! text stings for messages
      !k
      debug_print_1: if (ptc_track_debug) then
         print *; print *, "<subr. ptc_track_ini_conditions>:"
         print *,'  Initialise orbit, emittances and eigenvectors etc.'
      end if debug_print_1
      !k      j = 0
      j_particle_line_counter=0 !k

      debug_print_2: if (ptc_track_debug) then
         print *, " j_particle_line_counter =", j_particle_line_counter
         print *, 'Get values and variables:  '
         print *, ' twopi  =',  twopi
         print *, ''
      end if debug_print_2

1     continue  ! <---- loop for particle reading -----<--------<-------<-------<-------!
      !                                                                                 !
      !k      j  =  next_start(x,px,y,py,t,deltae,fx,phix,fy,phiy,ft,phit)              !
      ! read command line from input file :                                             !
      ! ptc_start,x=0.001,px=4.492437849e-06;                                           l
      ! madxdict.h:415:"ptc_start: ptc_start none 0 0 "                                 o
      ! madxn.c:3384:                                                                   o
      !      int next_start(xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx)     p
      !/* returns the parameters of the next particle to track;                         !
      !        0 = none (there is no any more particles),                               !
      !        else = count the particle number */                                      f
      j_particle_line_counter  =  &                                                   ! o
           next_start(x_input,px_input,y_input,py_input,t_input,deltae_input, &       ! r
           fx_input,phix_input,fy_input,phiy_input,ft_input,phit_input)                 !
      debug_print_3: if (ptc_track_debug) then                                          !
         print *, "The next command line from input file is read."                    ! p
         print *, "The keyword <ptc_start> provides the following data:"              ! a
         print *, "j_particle_line_counter=",j_particle_line_counter                  ! r
         print *, "x,px,y,py,t,deltae,fx,phix,fy,phiy,ft,phit=", &                    ! t
              x_input,px_input,y_input,py_input,t_input,deltae_input, &               ! i
              fx_input,phix_input,fy_input,phiy_input,ft_input,phit_input             ! c
      end if debug_print_3                                                            ! l
      !                                                                               ! e
      !     the particle counter is inside <next_start> c-code                          !
      !                                                                                 !
      new_particle: if(j_particle_line_counter.ne.0) then ! if (new particle)=======!   !
         ! new particle is read                                                     !   r
         !                                                                          !   e
         if_switch_or_j: if(ptc_switch.lt.3.OR.j_particle_line_counter.eq.1) then   !   a
            ! if switch = 1,2 (RUN, DYNAP fastune) OR ! ++++++++++++++++++++++++!   !   d
            !        one particle for DYNAP aperture                            I   !   i
            !                                                                   F   n   n
            jmax_numb_particl_at_i_th_turn = j_particle_line_counter            !   e   g
            !                                                                   !   w   !
            !  if (jmax_numb_particl_at_i_th_turn .GE.N_particle_max) then      !   !   !
            !  Print*,'!!! ERROR number of particles >= maximum'; STOP; ENDIF   !   !   !
            !                                                                   !   !   !
            ALLOCATE &                                                          !   !   !
                 (x_coord_incl_co(1:6,1:jmax_numb_particl_at_i_th_turn))        !   !   !
            x_coord_incl_co=zero                                                !   !   !
            !                                                                   !   !   l
            track_temp(:)=zero                                                  !   !   o
            !---- Get start coordinates                                         s   !   o
            X_MAD(1)=x_input; X_MAD(2)=px_input                               ! w   p   p
            X_MAD(3)=y_input; X_MAD(4)=py_input                               ! i   a   !
            X_MAD(5)=t_input; X_MAD(6)=deltae_input                           ! t   r   !
            !                                                                 ! c   t   f
            IF(icase_ptc.eq.5 .AND. (.NOT.closed_orbit)) THEN !--!            ! h   i   o
               if(mytime) then !----------------------!          !              !                   
                  call Convert_dp_to_dt (deltap, dt)  !          !              !
               else                                   !          !              !
                  dt=deltap                           !          !              !
               endif !--------------------------------!          !              !
               X_MAD(6)=X_MAD(6)+dt                              !              !   l   !
            ENDIF !----------------------------------------------!              !   e   !
            !                                                                   !   !   p
            CALL Coord_MAD_to_PTC(X_MAD,X_PTC) ! convert coordinates          ! c   t   a
            DO k_th_coord=1,6                                                 ! h   i   r
               track_temp(k_th_coord)=X_PTC(k_th_coord)                         !   c   t
            END DO                                                              !   l   i
            track_temp(7) = fx_input                                            !   e   c
            track_temp(8) = phix_input                                          !   !   l
            track_temp(9) = fy_input                                          ! O   !   e
            track_temp(10) = phiy_input                                       ! n   !   !
            track_temp(11) = ft_input                                         ! e   !   !
            track_temp(12) = phit_input                                         !   !   l
            debug_print_4: if (ptc_track_debug) then !----!                     !   !   o
               Print *, 'icase_PTC=', icase_PTC           !                     !   !   o
               Print *, 'track_temp(1:12)=',track_temp    !                     !   !   p
            end if debug_print_4 !------------------------!                     !   !   !
            !                                                                   !   !   !
            watch_non_zero_coords: do k_th_coord = 1,12 ! $$$$$$$$$$$$$$$$$!    P   !   !
               ne_zero: if(abs(track_temp(k_th_coord)).ne.zero) then !-- !  !   a   !   !
                  itype_non_zero_flag(k_th_coord) = 1                    !  !   r   !   p
               else !----------------------------------------------------!  !   t   !   a
                  itype_non_zero_flag(k_th_coord) = 0                    !  !   i   !   r
               endif ne_zero !-------------------------------------------!  !   c   n   t
               !                                                            !   l   e   i
            enddo watch_non_zero_coords ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!   e   w   c
            !                                                                   !   !   l
            !---- Normalized coordinates.                                       !   !   e
            !k    Tranformation from the polar to cartesian coordinates         !   p   !
            !k    on the plane of the normalized coordinates with a negative    !   a   !
            !k    direction for rotation (Y moves to X)                         !   r   !
            !                                                                   !   t   !
            Z_norm_temp(:)=zero                                                 !   i   !
            polar_cart: do kq = 1, icase_ptc-1, 2 ! --------------!             !   c   !
               !do kq = 1, 5, 2                                   !             !   l   !
               !             !  kq=1,3,5 - coord                  !             !   e   !
               kp = kq + 1   ! kp=2,4,6 - momentum                !             !   !   !
               !k                                                 !             !   !   !
               !k      phi = twopi * track(kq+7)       ! angle    !             !   !   !
               phi_n = twopi *  track_temp(kq+7)                  !             !   !   !
               !                                                  !             !   !   !
               !k    zn(kq) =   track(kq+6) * cos(phi) ! coord    !             !   !   !
               Z_norm_temp(kq) = track_temp(kq+6)*cos(phi_n)      !             !   !   !
               !                                                  !             !   !   !
               !k    zn(kp) = - track(kq+6) * sin(phi) ! momentum !             !   !   !
               Z_norm_temp(kp) =-track_temp(kq+6)*sin(phi_n)      !             !   n   !
               !       the minus sign is due to coord.system      !             !   e   !
               !                    with a negative rotation      !             !   w   !
               !k                                                 !             !   !   !
            enddo polar_cart !------------------------------------!             !   !   !
            !k                                                                  !   p   !
            debug_print_5:if (ptc_track_debug) then !----------------------!    !   a   !
               Print *,'Z_norm_temp(k_th_coord)=', Z_norm_temp(:icase_ptc) !    !   r   !
            end if debug_print_5 !-----------------------------------------!    !   t   !
            !                                                                   !   i   !
            !-- Transform to unnormalized coordinates and refer to closed orbit.!   c   !
            Z_stdt_temp(:)=zero;                                                !   l   !
            !                                                                   !   e   !
            if(closed_orbit) Z_stdt_temp=A_t_map*Z_norm_temp(:icase_ptc)        !   !   !
            !                     ! ANY ORDER !!!                               !   !   !
            debug_print_6: if (ptc_track_debug) then ! --------!                !   !   !
               Print *, 'Z_stdt_temp=A_t_map*Z_norm_temp=', &  !                !   !   !
                    Z_stdt_temp(:icase_ptc)                    !                !   !   !
            end if debug_print_6 !-----------------------------!                !   !   !
            !                                                                   !   !   !
            zgiv_exist = .false. ! there is no non-zero                         !   !   !
            zngiv_exist = .false.! action-angle input                           !   !   !
            !k                                                                  !   !   !
            !k  See Ch7 in PhysGuide on Method by Edwards & Teng                !   !   !
            !k      MAD-8 User Ref. Manual, Ch13, p.13.4                        !   !   !
            !k    do k = 1, 6 !----------------------------------------------!  !   !   !
            !k      if (itype(k) .ne. 0) zgiv = .true.                       !  !   !   !
            !k      if (itype(k+6) .ne. 0) zngiv = .true.                    !  !   !   !
            !k      zstart(k) = track(k)                                   & !  !   !   !
            !k      + sqrt(ex) * (eigen(k,1) * zn(1) + eigen(k,2) * zn(2)) & !  !   !   !
            !k      + sqrt(ey) * (eigen(k,3) * zn(3) + eigen(k,4) * zn(4)) & !  !   !   !
            !k      + sqrt(et) * (eigen(k,5) * zn(5) + eigen(k,6) * zn(6))   !  !   !   !
            !k    enddo !----------------------------------------------------!  !   !   !
            !                                                                   !   !   !
            z_start_sum(:)=zero                                                 !   !   !
            SumDO: DO k_th_coord = 1,icase_ptc ! not 6 !----------------------! !   !   !
               IF (itype_non_zero_flag(k_th_coord).ne.0) zgiv_exist =.true.   ! !   !   !
               IF (itype_non_zero_flag(k_th_coord+6).ne.0) zngiv_exist =.true.! !   !   !
               z_start_sum(k_th_coord)=track_temp(k_th_coord)+ &              ! !   !   !
                    Z_stdt_temp(k_th_coord)                                   ! !   !   !
            ENDDO SumDO !-----------------------------------------------------! !   !   !
            !                                                                   !   !   !
            debug_print_7: if (ptc_track_debug) then !-----!                    !   !   !
               Print *, 'z_start_sum=',  z_start_sum       !                    !   !   !
            end if debug_print_7 !-------------------------!                    !   !   !
            !                                                                   !   !   !
            !k          if (switch .gt. 1)  then !------------!                 !   !   !
            !--- keep initial coordinates for dynap           !                 !   !   !
            !k            do k = 1, 6 !-----------------!     !                 !   !   !
            !k              coords(k,0,j) = zstart(k)   !     !                 !   !   !
            !k            enddo !-----------------------!     !                 !   !   !
            !k          endif !-------------------------------!                 !   !   !
            !                                                                   !   !   !
            !---- Warn user about possible data conflict.                       !   !   !
            Warn: if (zgiv_exist .and. zngiv_exist) then !------------!         !   !   !
               msg(1) = 'Absolute and normalized coordinates given,'  !         !   !   !
               msg(2) = 'Superposition used.'                         !         !   !   !
               call fort_warn('START-1: ', msg(1))                    !         !   !   !
               call fort_warn('START-2: ', msg(2))                    !         !   !   !
            endif Warn !----------------------------------------------!         !   !   !
            !                                                                   !   !   !
            ! fill a previous particles                                         !   !   !
            previous: IF (jmax_numb_particl_at_i_th_turn.GE.2) THEN !*********! !   !   !
               DO j_th_particle=1,jmax_numb_particl_at_i_th_turn-1 !=======!  ! !   !   !
                  DO k_th_coord = 1, icase_ptc ! not 6 !----------------!  !  ! !   !   !
                     x_coord_incl_co(k_th_coord, j_th_particle)= &      !  !  ! !   !   !
                          x_coord_incl_co_temp(k_th_coord,j_th_particle)!  !  ! !   !   !
                  END DO !----------------------------------------------!  !  ! !   !   !
               END DO !====================================================!  ! !   !   !
               DEALLOCATE (x_coord_incl_co_temp)                              ! !   !   !
            END IF previous !*************************************************! !   !   !
            !                                                                   !   !   !
            ! fill current particle                                             !   !   !
            current: do k_th_coord = 1, icase_ptc !6 !-------------!            !   !   !
               ! z(6,1:jend) - Transformed cartesian coordinates   !            !   !   !
               !               including closed orbit ( c.o.)      !            !   !   !
               !k            z(k,j) = orbit0(k) + zstart(k)        !            !   !   !
               x_coord_incl_co(k_th_coord,             &           !            !   !   !
                    jmax_numb_particl_at_i_th_turn)=   &           !            !   !   !
                    x_coord_co(k_th_coord)+z_start_sum(k_th_coord) !            !   !   !
            enddo current !----------------------------------------!            !   !   !
            !                                                                   !   !   !
            ! save already-read particles in the array                          !   !   !
            ALLOCATE &                                                          !   n   !
                 (x_coord_incl_co_temp(1:6,1:jmax_numb_particl_at_i_th_turn))   !   e   !
            x_coord_incl_co_temp=zero                                           !   w   !
            !                                                                   !   !   !
            save: DO j_th_particle=1,jmax_numb_particl_at_i_th_turn !==!        !   !   !
               DO k_th_coord = 1, icase_ptc ! not 6 !----------------! !        !   !   !
                  x_coord_incl_co_temp(k_th_coord, j_th_particle)= & ! !        !   !   !
                       x_coord_incl_co(k_th_coord,j_th_particle)     ! !        !   p   !
               END DO !----------------------------------------------! !        !   a   !
            END DO save !==============================================!        !   r   !
            DEALLOCATE (x_coord_incl_co)                                        !   t   !
            !                                                                   !   i   i
         END IF if_switch_or_j !++++++++++++++++++++++++++++++++++++++++++++++++!   c   n
         ! if switch = 1,2 (RUN, DYNAP fastune) OR                              !   l   g
         !                         one particle for DYNAP aperture              !   e   !
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!   !   !
         !                                                                          !   !
         goto 1  ! -->-- loop for particle reading ----->-------->------->----------!---!
         !                                                                          !
      endif new_particle ! ======END IF (new particle) =============================!
      ALLOCATE &
           (x_coord_incl_co(1:6,1:jmax_numb_particl_at_i_th_turn))
      x_coord_incl_co=zero
      x_coord_incl_co=x_coord_incl_co_temp

      Array_for_saving_all_particles: IF (.NOT.element_by_element) THEN !***********!
         ALLOCATE (x_all_incl_co_at0(1:6,0:turns,1:jmax_numb_particl_at_i_th_turn)) !
         x_all_incl_co_at0=zero                                                     !
         DO j_th_particle=1,jmax_numb_particl_at_i_th_turn !===========!            !
            DO k_th_coord = 1, icase_ptc   !------------------------!  !            !
               x_all_incl_co_at0(k_th_coord, 0, j_th_particle)=  &  !  !            !
                    x_coord_incl_co(k_th_coord,    j_th_particle)   !  !            !
            END DO !------------------------------------------------!  !            !
         END DO !======================================================!            !
      ENDIF Array_for_saving_all_particles !****************************************!

      DEALLOCATE (x_coord_incl_co_temp)

      !k      if (switch .eq. 3)  then
      !--- create second particle with x add-on
      !k        deltax = get_value('dynap ', 'lyapunov ')
      !k        jend = 2
      !k        z(1,jend) = z(1,1) + deltax
      !k        do k = 2, 6
      !k          z(k,jend) = z(k,1)
      !k        enddo
      !k      endif
      !k      end

      debug_printing: if (ptc_track_debug) then
         print *,' Subr. ptc_track_ini_conditions is finished'
         print *,' The number of particles is jmax_=', jmax_numb_particl_at_i_th_turn
         print *,' Array of coordinate is '
         part: DO j_th_particle=1,jmax_numb_particl_at_i_th_turn !---------------------!
            Print *,'Particle No.=', j_th_particle                                     !
            Print *, (x_coord_incl_co(k_th_coord,j_th_particle),k_th_coord=1,6)  !???  !
         END DO part !-----------------------------------------------------------------!
      end if debug_printing

    END SUBROUTINE ptc_track_ini_conditions

    !=============================================================================

    !=============================================================================

    SUBROUTINE Coord_MAD_to_PTC(X_MAD,X_PTC)
      ! Convert coordinates from MAD to PTC
      IMPLICIT NONE
      REAL(dp), INTENT(IN)  :: X_MAD(6)
      REAL(dp), INTENT(OUT) :: X_PTC(6)

      X_PTC=X_MAD ! for all elements

        if (ptc_track_debug) print *, &
         'Coord_MAD_to_PTC icase_ptc=', icase_ptc, ' mytimec=', mytime

      if (icase_ptc.eq.6) THEN
         X_PTC(5)=X_MAD(6); X_PTC(6)=X_MAD(5);
         if (mytime) X_PTC(6)=-X_PTC(6) ! reverse sign
      elseif(icase_ptc.eq.5) THEN
         X_PTC(5)=X_MAD(6); X_PTC(6)=zero
      elseif(icase_ptc.eq.4) THEN
         X_PTC(5)=zero; X_PTC(6)=zero      
      ENDIF
    END SUBROUTINE Coord_MAD_to_PTC

    SUBROUTINE Coord_PTC_to_MAD(X_PTC,X_MAD)
      ! Convert coordinates from PTC to MAD
      IMPLICIT NONE
      REAL(dp), INTENT(IN)  :: X_PTC(6)
      REAL(dp), INTENT(OUT) :: X_MAD(6)

      X_MAD=X_PTC ! for all elements

      if (ptc_track_debug) print *, &
          'Coord_PTC_to_MAD icase_ptc=', icase_ptc, ' mytime=', mytime

      IF (icase_ptc.eq.6) THEN
         X_MAD(5)=X_PTC(6); X_MAD(6)=X_PTC(5);
         if (mytime) X_MAD(5)=-X_MAD(5) ! reverse sign         
      elseif(icase_ptc.eq.5) THEN
         X_MAD(5)=X_PTC(6); X_MAD(6)=X_PTC(5)
         if (mytime) X_MAD(5)=-X_MAD(5) ! reverse sign
      elseif(icase_ptc.eq.4) THEN
         X_MAD(5)=zero; X_MAD(6)=zero      
      ENDIF
    END SUBROUTINE Coord_PTC_to_MAD
    !=============================================================================

  END subroutine ptc_track_run
  !==============================================================================
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  ! INDEPENTENT  SUBROUTINES of this module
  !==============================================================================
  !==============================================================================
  SUBROUTINE Get_map_from_NormalForm &
       (ptc_track_debug, Normal_Order_n0, x_coord_co, &
       Map_Y, Normal_Form_N, A_t_map, A_t_map_rev)

    USE  madx_ptc_module, ONLY: dp, real_8, normalform, damap,   &
         my_ring, default, BERZ, daprint, &
         assignment(=), operator(**), operator(.sub.), &
         track, init, alloc, kill, &
         PRODUCE_APERTURE_FLAG,  ANALYSE_APERTURE_FLAG
    ! USE ptc_results
    implicit none

    LOGICAL(lp),           INTENT(IN) :: ptc_track_debug
    INTEGER,           INTENT(IN) :: Normal_Order_n0   ! =1 for Linear
    REAL (dp),         INTENT(IN)  ::  x_coord_co(1:6) ! => x0(1:6) in ptc_track
    TYPE (real_8),     INTENT(OUT)  :: Map_Y(6)        !  y => Map_Y - local name
    TYPE (normalform), INTENT(OUT) :: Normal_Form_N    !  n =>  Normal_Form_N - local name
    TYPE (damap),      INTENT(OUT) :: A_t_map          ! x_stdt=A_t_map*x_norm
    TYPE (damap),      INTENT(OUT) :: A_t_map_rev      ! x_norm=A_t_map_rev*x_norm
    ! A_t_map_rev=A_t_map^(-1)

    integer :: mynd2,nda, flag_index,why(9) ! icase,
    ! integer :: npara ! Global in module

    REAL (dp) :: d_val(6) 
    INTEGER   :: i_vec, i_comp, ind(6)

    !------------------------------------------------------------------------------
    IF (ptc_track_debug)  print *, 'Start Subr.  Get_map_from_NormalForm '

    nda=0
    ! no = 1 =>  Normal_Order_n0 (INPUT)

    call init(default,Normal_Order_n0,nda,BERZ,mynd2,npara)
    if (ptc_track_debug) then
       Print *, 'aftercall init(default,Normal_Order_n0,nda,BERZ,mynd2,npara)'
       Print *, 'Normal_Order_n0,nda,mynd2,npara=', Normal_Order_n0,nda,mynd2,npara
    endif

    call alloc(Map_Y) ! call alloc(y)

    Map_Y=npara ! y=npara
    Map_Y=x_coord_co  ! Y=X
    IF (ptc_track_debug) call daprint(Map_Y,16);
    !    c_%watch_user=.true.

    call track(my_ring,Map_Y,1,default)
    IF (ptc_track_debug)call daprint(Map_Y,17);


    call PRODUCE_APERTURE_FLAG(flag_index)
    if(flag_index/=0) then
       call ANALYSE_APERTURE_FLAG(flag_index,why)
       Write(6,*) "Get_map_from_NormalForm is unstable (map production)-programs continues "
       Write(6,*) why ! See produce aperture flag routine in sd_frame
       CALL kill(Map_Y) !CALL kill(y)
       ! c_%watch_user=.false.
       return
    endif
    !    c_%watch_user=.false.

    call alloc(Normal_Form_N) ! call alloc(n)
    ! Normal_Form_N=0 !Error: Can't convert INTEGER(4) to TYPE(normalform)

    Normal_Form_N=Map_Y ! n=y

    call alloc(A_t_map)

    A_t_map=Normal_Form_N%A_t

    call alloc(A_t_map_rev)

    A_t_map_rev=Normal_Form_N%A_t**(-1)

    if (ptc_track_debug) then
       OPEN (UNIT=17,FILE='Normal_FORM_A.txt', STATUS='UNKNOWN')
       write(17,'(/a/)') 'Map_Y: '; call daprint(Map_Y,17);
       write(17,'(/a/)') 'Linear Normal_Form_N%A: '; call daprint(Normal_Form_N%A,17);     !
       write(17,'(/a/)') 'Linear Normal_Form_N%A_t: '; call daprint(Normal_Form_N%A_t,17); !
       write(17,'(/a/)') 'Linear A_t_map: '; call daprint(A_t_map,17);                     !
       write(17,'(/a/)') 'Linear A_t_map_rev: '; call daprint(A_t_map_rev,17);             !
       Close (17)
       ! CALL TEST_PTC_Normal(Normal_Form_N)
    end if

   !EigenVectors
    !allocate (d_val(1:icase_ptc))
    do i_vec=1,icase_ptc
      do i_comp=1,6 !icase_ptc
        ind(:)=0; ind(i_comp)=1
        d_val(i_comp)=Normal_Form_N%A_t%V(i_vec).sub.ind
      enddo
        if (ptc_track_debug) then
         WRITE(17,*) 'EigenVector V(',i_vec,')=', d_val  
        endif
      !enddo  
     end do
     !DEallocate (d_val)
  END subroutine Get_map_from_NormalForm
  !==============================================================================


END MODULE madx_ptc_track_run_module
!==============================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
