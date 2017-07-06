MODULE madx_ptc_track_run_module
  ! This module serve as a COMMON block as in  F77
  ! It contains variables which exchange data between
  ! SUBROUTINE ptc_track_run and called from it
  ! external subroutines calculating particle interactions
  use file_handler
  USE madx_ptc_module , ONLY: dp, lp, lnv, &
                                ! shorts for <double precision>, <logical>, 0D0 etc.
       doublenum ! am temprorary double number for I/O with C-procedures
  USE madx_ptc_intstate_module, ONLY: getdebug  ! new debug control by PS (from 2006.03.20)
  use name_lenfi
  use definition
  implicit none
  SAVE
  PRIVATE
  TYPE(INTERNAL_STATE) MYSTATE

  PUBLIC :: ptc_track_run ! Subroutine inside the module

  !------------------------------------------------------------------!
  !  Variables from input files and probably corrected by this code: !
  !------------------------------------------------------------------!
  INTEGER, PUBLIC :: icase_ptc    ! Phase-space (4, 5 or 6) in input command
  INTEGER, PUBLIC :: nvariables   !  actul number of variables
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
       recloss, &
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

  CHARACTER(name_len), ALLOCATABLE  :: name_el_at_obsrv(:) ! the name of an element at observation point
  !                                                  ! contrary to <character(16)> in c-code

  character(1000), private  :: whymsg

  !real(KIND(1d0)) :: dble_num_C  ! to use as a temprorary double number for I/O with C-procedures

CONTAINS

  SUBROUTINE ptc_track_run(max_obs)

    !USE MADX_PTC_MODULE ====================================================================!
    USE  madx_ptc_module, ONLY: universe, EXCEPTION, my_ring, default, index_mad, c_         !
    !                                                                                        !
    USE  madx_ptc_module, ONLY: &  ! "LAYOUT type (ring) => double linked list,              !
         FIBRE, &                  !  whose nodes (elements=magnets) of type FIBRE"          !
         NORMALFORM, & ! type for normalform                                                 !
         REAL_8, &     ! type for map                                                        !
         damap, &      ! type for diff algebra                                               !
         RADIATION ! STOCH_IN_REC, & ! type for radiation with quadrupoles in PTC            !
    ! ======== functions ==================================================================!
    USE  madx_ptc_module, ONLY: &                                                          !
         print, find_orbit,FIND_ORBIT_x, track,track_probe_x,UPDATE_STATES, my_state, &    !
         PRODUCE_APERTURE_FLAG, ANALYSE_APERTURE_FLAG, &                                   !
         kill, daprint, alloc, Get_one, &                                                  !
         assignment(=), operator(+), operator(*), operator(.sub.), &                       !
         ! Coord_MAD_to_PTC, Coord_PTC_to_MAD,  & => at the end of this module             !
         write_closed_orbit,Convert_dp_to_dt,mytime                                        !
    !======================================================================================!
    USE  madx_ptc_module, ONLY: &                                                          !
         c_1d_7,c_1D3,one, twopi, zero !, two                                              !
    !======================================================================================!

    USE Inf_NaN_Detection !VK20070328 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    use bbfi                ! integer bbd_loc,bbd_cnt,bbd_flag,bbd_pos,bbd_max;
    !                       !uses bbd_pos                parameter(bbd_max=200)
    !                       !real(dp) bb_kick
    !                       ! common/bbi/bbd_loc(bbd_max),bbd_cnt,bbd_flag,bbd_pos
    !                       !common/bbr/bb_kick(2,bbd_max)
    use name_lenfi   ! integer name_len;  parameter(name_len=24)

    IMPLICIT NONE

    integer, intent (IN) :: max_obs ! the maximum number of observation points >=1
    !                               ! one point at the end (beginning)+
    !                               ! points given in input file by the command
    !                               ! "ptc_observe,place=mark";


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

    !integer :: mad8_code_current_node  ! = code <in trrun.F>
    ! code_of_the_current_node_as_mad8
    ! REAL(dp) :: el_element_length       ! = el

    Integer, allocatable  :: last_turn_of_lost_particle(:)! (1:N_particle_max)  != last_turn

    REAL(dp) , allocatable :: &
         last_position_of_lost_particle(:), & !(1:N_particle_max), & != last_pos
         last_orbit_of_lost_particle (:,:) !(1:6,1:N_particle_max)   != last_orbit

    !   /* C routines called from Fortran and C */
    real(KIND(1d0)), external :: get_value, get_variable ! external c-functions
    INTEGER, external :: get_option, &   !  int get_option(char*);
         restart_sequ, & !  restart beamline and return number of beamline node
         advance_node    !  advance to the next node in expanded sequence
    !                    !  =0 (end of range), =1 (else)

    REAL(KIND(1d0)), external :: node_value !/*returns value for parameter par of current element */

    EXTERNAL :: comm_para ! subroutine needed for LF95

    !k    character(12) char_a
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
    REAL(dp) :: summ_ring_length      ! saved ring-length calculated before tracking
    ! at the last element in subr. Prepare_Observation_points
    integer ::  nlm_current_element_number    ! line position =nlm
    character(name_len) el_name
    !hbu
    character(4) vec_names(7)
    !hbu
    data vec_names  / 'x', 'px', 'y', 'py', 't', 'pt','s' / ! MADX
    !data vec_names / 'x', 'px', 'y', 'py', 'pt', 't','s' / ! PTC has a reverse order for pt and t
    logical(lp) rplot
    integer   mft ! debug output file
    !k    data char_a / ' ' /
    !-------------------------------------------------------------------------
    
    if (getdebug() > 2) then
      ptc_track_debug = .true.
      call kanalnummer(mft)
      open(unit=mft,file='ptc_track.log')
      
    endif

    rplot = get_value('ptc_track ','rootntuple ') .ne. 0
    if (rplot) then
       call newrplot()
    endif
    
    
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
    if (getdebug() > 1) then
      Print *;  Print *,'  ================================================================'
      Print *, '  ptc_track: The current dimensionality of the problem is icase=', icase_ptc
      Print *, '  ptc_track: nvariables=', nvariables
      Print *,'  ================================================================'; Print *;
    endif
    
    warn_coordinate_system_changed: IF((.NOT. mytime) .AND.(nvariables.gt.4) ) THEN
       CALL FORT_WARN('time=false => coord. system: {-pathlength, delta_p} ', &
            'the table headers mean:  PT -> delta_p, T -> pathlength')
    ENDIF warn_coordinate_system_changed

    ! initialize the closed orbit coordinates  at START of the ring
    x_coord_co(:)=zero
    if (ptc_track_debug) then
        print *, " x_coord_co(:)=zero = ",x_coord_co
    endif

    ! Closed_orbit_at_START:
    IF(closed_orbit) CALL Find_Closed_Orbit   ! Calculates x_coord_co(1:6)

    Normal_forms: IF(closed_orbit) THEN  !-----------------------!
       !                                                         !
       !                                                         !
       call Get_map_from_NormalForm &                            !
            (ptc_track_debug, Normal_Order_n0, x_coord_co, &     !
       ! INTENT:   IN,             IN,             IN            !
       !                1=> Normal_Order_n0  for Linear          !
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
!    change_default: IF(Radiation_PTC) THEN
!       MYSTATE=DEFAULT+RADIATION
!       IF (Radiation_Quad) STOCH_IN_REC=.TRUE.
!       !element_by_element=.FALSE. ! make PTC one-turn tracking
!       print *, "################################################################"
!       print *, "The PTC parameter DEFAULT before the tracking with the turn-loop"
!       call print(MYSTATE,6)
!       CALL Find_Closed_Orbit ! Calculates x_coord_co(1:6)
!    END IF change_default

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !!
    !!  WRITE DOWN
    !!

    Call Final_Coord_to_tables ! Complete all tables by one subroutine:
    !'trackone' filled before  by   tt_puttab_coord, 
    !'track.obs$$$$.p$$$$'     by   tt_putone_coord  
    !summary table             by  'tracksumm'
    

    Output_observ_with_PTC: IF(closed_orbit .AND. &
         (.NOT.element_by_element).AND.(.NOT. Radiation_PTC)) THEN !-!
       
       
       debug_print_5: if (ptc_track_debug) then !----------------!                         !
          Print *, 'element_by_element=', element_by_element, &  !                         !
               ' Radiation_PTC=',Radiation_PTC                   !                         !
          Print *, ' Call Observation_with_PTC'                  !                         !
       end if debug_print_5 !------------------------------------!                         !
       
       Call Observation_with_PTC(max_obs, x_coord_co, Map_Y)                               !

    ENDIF Output_observ_with_PTC !---------------------------------------------------------!


    if (ptc_track_debug) then
        Print *, 'Come to : <! Calculate beam envelope with PTC>'
    endif


    if (rplot) call rplotfinish()


    Beam_envelope_with_PTC: IF (beam_envelope) THEN 
        call fort_warn('ptc_track: ',' Calculation of Equilibrim emittance was moved to ptc_twiss')
    endif Beam_envelope_with_PTC

!     Beam_envelope_with_PTC: IF (beam_envelope) THEN !###############################!
!        Radiat_PTC: IF ( Radiation_PTC) THEN !====================================!  !
!           icase_6: IF (icase_ptc .EQ. 6 .AND.closed_orbit) THEN !--!             !  !
!              Call  beam_enevelope_with_PTC                         !             !  !
!           ELSE !---------------------------------------------------!             !  !
!              Print *, ' Warning !: Option BEAM_ENVELOPE', &        !             !  !
!                   ' require option ICASE=6', &                     !             !  !
!                   ' and option CLOSED_ORBIT '                      !             !  !
!              Print *, '    The code ignores BEAM_ENVELOPE option ' !             !  !
!           ENDIF icase_6 !------------------------------------------!             !  !
!        ELSE !====================================================================!  !
!           Print *, ' Warning !!!: Option BEAM_ENVELOPE require option Radiation' !  !
!           Print *, '              The code ignores BEAM_ENVELOPE option '        !  !
!        ENDIF Radiat_PTC !========================================================!  !
!     ENDIF Beam_envelope_with_PTC !##################################################!

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
      ! USE madx_ptc_module, ONLY: universe, index_mad
      implicit none

      ! It checks that the following commands has been performed
      ! before calling subr. ptc_track
      ! If not, then return_from_subr_ptc_track =>.TRUE.
      !                                        (accessible from the HOST subr.)

      ! 1) ptc_create_universe:   Needed to set-up PTC

      if (ptc_track_debug) then
         print *
         print *, '<subr. Check_ptc_universe_and_layout (before_tracking)>:'
         print *, 'PTC Universe universe=', universe
         print *, 'PTC Layout      index=', index_mad
      end if
      if(universe.le.0.or.EXCEPTION.ne.0) then
         call fort_warn('return from ptc_track: ',' no universe created')
         ! return
         return_from_subr_ptc_track=.TRUE.
      endif

      ! 2) ptc_create_layout:     creates PTC layout and fills it
      !                                          with current MAD-X sequence

      if(index_mad.le.0.or.EXCEPTION.ne.0) then
         call fort_warn('return from ptc_track: ',' no layout created')
         ! return
         return_from_subr_ptc_track=.TRUE.
      endif
    END  SUBROUTINE Check_ptc_universe_and_layout ! before_tracking =>  Int. proc. (f95)

    !=============================================================================

    SUBROUTINE Values_from_ptc_track_command ! Internal procedure (f95)
      USE madx_ptc_intstate_module, ONLY: getdebug  ! new debug control by PS (from 2006.03.20)
      implicit none
      ! local variables
      character(4) text
      character(12) tol_a, char_a
      integer :: nint,ndble, nchar, int_arr(1),char_l
      data tol_a,char_a / 'maxaper ', ' ' /
      real(dp) :: tmpr
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
      if(icase_ptc.eq.56) nvariables=5
      if(icase_ptc.ne.4.and.icase_ptc.ne.5.and.icase_ptc.ne.6.and.icase_ptc.ne.56) then
         write(text, '(i4)') icase_ptc
         call aafail('Values_from_ptc_track_command: ',' ICASE not 4, 5, 6 or 56 found: ' // text)
      endif
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

      IF ( Radiation_PTC .AND. radiation_model1_FZ ) THEN
         radiation_model1_FZ=.FALSE.
         Print *, ' '
         Print *, '===================================================================='
         Print *,'Both options RADIATION and RADIATION_MODEL1 are entered as TRUE !? '
         Print *,'RADIATION has a higher priority. RADIATION_MODEL1 set to be OFF.'
         call fort_warn(' RADIATION_MODEL1',' has been switched OFF by the code')
         Print *, '===================================================================='
         Print *, ' '        
      ENDIF

      Radiation_Energy_Loss = get_value('ptc_track ','radiation_energy_loss ') .ne. 0

      Radiation_Quad = get_value('ptc_track ','radiation_quad ') .ne. 0

      Space_Charge = get_value('ptc_track ','space_charge ') .ne. 0
      
      recloss      = get_value('ptc_track ','recloss ') .ne. 0

      IF (radiation_model1_FZ .OR. Space_Charge ) THEN 
         IF(.NOT.element_by_element) THEN
            element_by_element=.TRUE.
            Print *, ' '
            Print *, '===================================================================='
            call fort_warn(' ELEMENT_BY_ELEMENT',' has been switched ON by the code')
            Print *,'     Only element-by-element tracking can be performed with  '
            Print *,'     the options RADIATION_MODEL1 or SPACE_CHARGE.'
            Print *, '===================================================================='
            Print *, ' '        
         ENDIF
      ENDIF

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
      if (getdebug()>3) ptc_track_debug=.true.  ! ptc_track_debug=.T., if debuglevel.ge.3 ,i.e.,
      ! in the madx input: ptc_setswitch, debuglevel=4

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
      
      !swap 5 and 6
      tmpr = maxaper(5)
      maxaper(5) = maxaper(6)
      maxaper(6) = tmpr

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
         print *, " recloss =",  recloss
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
      USE  madx_ptc_module, ONLY:  radiation0
      implicit none
      character(4) text

      debug_print_1: if (ptc_track_debug) then
         print *,"before <call my_state(icase,deltap,deltap0)>", &
              "with parameters:  "
         print *, 'icase_ptc=',   icase_ptc, 'deltap=', deltap
         print *, 'deltap0=', deltap0 ; print *
         print*, 'my_state printing:'
      end if debug_print_1
      
      call my_state(icase_ptc,deltap,deltap0)
      
      nvariables = icase_ptc
      !icase 56: 
      if (nvariables .gt.6 ) nvariables = 6
      
      IF (Radiation_PTC) DEFAULT=DEFAULT+RADIATION0
      
      if(getdebug() > 1) then
        Print *, ' Radiation_PTC    =     ', Radiation_PTC
      endif
      
      debug_print_2: if (ptc_track_debug) then
         print *, &
              "after call my_state(icase,deltap,deltap0)"
         print *, 'icase_ptc=',icase_ptc, 'deltap=', deltap
         print *, 'deltap0=', deltap0
         print *; print*, '----------------------------------'
         print *, "before CALL UPDATE_STATES"
      end if debug_print_2
       
      CALL UPDATE_STATES
      MYSTATE=DEFAULT
      
      if (getdebug()>1) then
        print *, "after CALL UPDATE_STATES"
        print *; print*, '----------------------------------'
        print *, "Printing by <call print(default,6)>:"
        call print(MYSTATE,6)
        print *, "after call print(MYSTATE,6)"
      endif
      
      if ( MYSTATE%TOTALPATH .gt. 0 ) then 
        print *, "********************************************************************"
        print *, "MAXAPER check of T variable not to be done because TOTALPATH is used"
        print *, "********************************************************************"
      endif
      
      
    END SUBROUTINE Call_my_state_and_update_states
    !=============================================================================

    !=============================================================================
    SUBROUTINE Find_Closed_Orbit
      ! USE madx_ptc_module, ONLY: dp, zero, find_orbit, my_ring,default
      implicit none
      !
      !====================================================================!
      !   initialize the closed orbit coordinates                          !
      ! x0(:)=zero                                                         !
      x_coord_co(:)=zero                                                   !
      if (ptc_track_debug) THEN !---------------------------!              !
         print *, " x_coord_co(:)=zero = "                  !              !
         CALL write_closed_orbit(nvariables,x_coord_co)      !              !
      end if !----------------------------------------------!              !
      !                                                                    !
      if(nvariables.ge.5) THEN !------------------------!                   !
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
         CALL FIND_ORBIT_x(my_ring,x_coord_co,MYSTATE,c_1d_7,fibre1=1)     !
         !         call find_orbit(my_ring,x_coord_co,1,MYSTATE,c_1d_7)    !
         print*,"===== ptc_track ============================"   !         !
         CALL write_closed_orbit(nvariables,x_coord_co)           !         !
         print*,"============================================"   !         !
      endif  !---------------------------------------------------!         !
      !                                                                    !
      if (ptc_track_debug) then
          print*,"After closed_orbit"; print *;                            !
      endif                                                                !
      !                                                                    !
      !END closed_orbit    which is logically (.not.ONEPASS)               !
      !====================================================================!

    END SUBROUTINE Find_Closed_Orbit
    !=============================================================================

    !=============================================================================
    SUBROUTINE DeAllocate_local_and_PTC_arrays
      ! USE madx_ptc_module, ONLY: kill
      implicit none

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
      implicit none

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
      implicit none

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
            CALL write_closed_orbit(nvariables,x_coord_co)
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
      implicit none

      integer :: j_th_particle, k_th_coord ! counter
      ! real(dp) :: tmp_d ! temprorary dble vaiable
      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      REAL (dp) :: X_MAD(6), X_PTC(6)

      ! Summary Table
      !--- enter start coordinates in summary table
      do  j_th_particle = 1, & !=====loop over particles =====-====================!
           j_tot_numb_starting_particles ! => initial number of particles          !
         !                                                                         !
         doublenum = j_th_particle                                                 !
         call double_to_table_curr('tracksumm ', 'number ', doublenum)                  !
         ! tmp_d = 1 <=  turn=1  in the original 2005 trrun.f                      !
         !                                                                         !
         doublenum = zero ! <=  turn=0  for starting particles                     !
         call double_to_table_curr('tracksumm ', 'turn ', doublenum)                    !
         DO k_th_coord = 1, 6 !>>>>> loop over coord. components >>>>>>>>>>>>>>!   !
            !tmp_d = z(k_th_coord,j_th_particle) - orbit0(k_th_coord)          !   !
            !z(1:6,1:j_tot) - coordinates                                      !   !
            !orbit0(1:6)    - (dble; 1:6) closed orbit                         !   !
            X_PTC(k_th_coord)=x_coord_incl_co(k_th_coord,j_th_particle) - &    !   !
                 x_coord_co(k_th_coord)                                        !   !
         ENDDO !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!   !
         !                                                                         !
         CALL Coord_PTC_to_MAD(X_PTC,X_MAD)                                        !
         !                                                                         !
         DO k_th_coord = 1, 6 !>>>>> loop over coord. components >>>>>>>>>>>>>>!   !
            doublenum  = X_MAD(k_th_coord)                                     !   !
            !                                                                  !   !
            call double_to_table_curr('tracksumm ',vec_names(k_th_coord),doublenum) !   !
         enddo !>>>>> END loop over components >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!   !
         Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C, &                    !
              gamma0I,gambet) ! to get "energy" value                              !
         doublenum=energy                !
         call double_to_table_curr('tracksumm ', 'e ', doublenum)                       !
         !                                                                         !
         !hbu add s                                                                !
         doublenum = spos_current_position             !
         call double_to_table_curr('tracksumm ',vec_names(7),doublenum)                 !
         call augment_count('tracksumm ')                                          !
         ! madxn.c:1094:void augment_count(char* table)                            !
         ! increase table occ. by 1, fill missing *                                !
      enddo ! ====loop over particle == do  i = 1,j_tot ===========================!

    END SUBROUTINE Start_Coord_to_TrackSum_Table

    !=============================================================================

    SUBROUTINE Enter_0st_turn_in_tables
      implicit none

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

               call tt_puttab_coord(particle_ID(j_th_particle),   0,   1, &     !   !  !
                                    x_coord_incl_co(1,j_th_particle), &         !   !  !
                                    x_coord_co, spos_current_position)          !   !  !
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
      implicit none

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

      LOGICAL :: NaN_coord_after_track_VK=.False. !VK20070328 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      !NaN_coord_after_track_VK=.False. !

      if (ptc_track_debug) then ! debug printing --------------------------!
         print *; print *, 'Start SUBR.<One_turn_track_with_PTC>'
      end if

      n_temp=1 ! start particle loop with the first not-examined particle
      
      
      Exam_all_particles: DO ! ====<=======<========DO Exam_all_particles ====<===<=======<======!
         !(instead of <10 continue> in trrun.f)                                                  !
         !                                                                                       !
         Particle_loop: DO j_particle=n_temp, jmax_numb_particl_at_i_th_turn  !++++++++++++++!   !
            !                                                                                !   !
            NaN_coord_after_track_VK=.False. !VK20070709 XXXXXXXXXXXXXXXXXXXXX               !   !
            !                                                                                !   !
            jmax_at_loop_start = jmax_numb_particl_at_i_th_turn                              !   ^
            j_last_particle_buffer=j_particle ! remember index value after END DO            !   !
            !                                                                                !   !
            do k_th_coord=1,6 ! extract coords for the current particle -----!               +   ^
               current_x_coord_incl_co(k_th_coord)= &                        !               +   !
                    x_coord_incl_co(k_th_coord,j_particle)                   !               +   !
            end do !---------------------------------------------------------!               +   !
            !                                                                                !   !
            call track_probe_x(my_ring,current_x_coord_incl_co,MYSTATE,fibre1=1)             !   !
!            call track(my_ring,current_x_coord_incl_co,1,MYSTATE)                           !   !
            ! The PTC subroutine " To TRACK the MY_RING for X coordinates                    +   !
            ! over one-turn in the DEFAULT state (citation, p. 25).                          +   !
            ! there is no any other an explicit description in KEK 2002-3 report             +   !
            !                                                                                !   ^
            
            
            
            do k_th_coord=1,6 ! save coordinates for the current particle ---!               +   !
               if (ISNAN(current_x_coord_incl_co(k_th_coord))) then !VK20070328 XXXXXXXXX    +   !
                  ! BUG !? Aperture does not work, if lattice with spread multipoles XXXX    +   !
                  NaN_coord_after_track_VK=.TRUE.                                   !XXXX    +   !
                  x_coord_incl_co(k_th_coord,j_particle)=999                        !XXXX    !   !
               else  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    !   !
                  x_coord_incl_co(k_th_coord,j_particle)=  &                 !               +   !
                       current_x_coord_incl_co(k_th_coord)                   !               +   !
               endif !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    +   !
               !                                                             !               +   !
            end do !---------------------------------------------------------!               !   ^


            !                                                                                !   !
            if (ptc_track_debug) then ! debug printing ----------------------------!         !   !
               Print *,'DO j_particle=n_temp, jmax_numb_particl_at_i_th_turn:'     !         !   !
               Print *, 'DO ',j_particle,'=',n_temp,',', &                         !         !   !
                    jmax_numb_particl_at_i_th_turn                                 !         !   ^
               !                                                                   !         +   !
               do k_th_coord=1,6  !VK20070328 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  !         +   !
                  if (ISNAN(current_x_coord_incl_co(k_th_coord))) then         !X  !         +   !
                     Print *, 'NAN-value for coordinate number ', k_th_coord   !X  !         +   !
                  else                                                         !X  !         +   !
                     Print*, 'k_th_coord=', k_th_coord, &                      !X  !         +   !
                          'current_x_coord_incl_co=', &                        !X  !         !   ^
                          current_x_coord_incl_co(k_th_coord)                  !X  !         !   !
                  endif                                                        !X  !         !   !
               enddo !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!X  !         !   !
            END if ! debug printing -----------------------------------------------!         !   !
            !                                                                                !   !
            call PRODUCE_APERTURE_FLAG(flag_index_ptc_aperture)                              !   !
            !                                                                                !   !
!            if (NaN_coord_after_track_VK) flag_index_ptc_aperture=100 !VK20070328 XXXXXXXXX !   !
!            if (ptc_track_debug) then                                                 !XXXX !   !
!                 print *,'flag_index_ptc_aperture is set to', flag_index_ptc_aperture !XXXX !   !
!            endif                                                                           !   !
            !                                                                                !   !
            if(flag_index_ptc_aperture/=0) c_%watch_user=.false.                             !   !
            !                                                                                !   !
            if (ptc_track_debug) then                                                        !   !
               print*,"ready for printing aperture flag!!!!!!!",flag_index_ptc_aperture      !   !
               print*,"real aperture flag: ",c_%aperture_flag                                !   !
            ! Sa_extend_poly.f90:98:  SUBROUTINE PRODUCE_APERTURE_FLAG(I)                    !   !
            !                                                                                !   !
                            ! debug printing ----------------------------!         !   ^
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
            ! added 2012.09.24                                                               !   !
            if (NaN_coord_after_track_VK) flag_index_ptc_aperture=100 !VK20070328 XXXXXXXXX  !   !

            if (abs(current_x_coord_incl_co(1)).ge.maxaper(1).or. &                          !   !
                abs(current_x_coord_incl_co(2)).ge.maxaper(2).or. &                          !   !
                abs(current_x_coord_incl_co(3)).ge.maxaper(3).or. &                          !   !
                abs(current_x_coord_incl_co(4)).ge.maxaper(4)) flag_index_ptc_aperture = 51  !   !
                
            if (nvariables .gt. 4) then                            !   !
              if (abs(current_x_coord_incl_co(5)).ge.maxaper(5)) flag_index_ptc_aperture = 52!   !
            endif 

            if ( (nvariables .gt. 5) .and. (MYSTATE%TOTALPATH .eq. 0) ) then                  !   !
              if (abs(current_x_coord_incl_co(6)).ge.maxaper(6)) flag_index_ptc_aperture = 53!   !
            endif 

            if (ptc_track_debug) then                                                 !XXXX  !   !
                 print *,'flag_index_ptc_aperture is set to', flag_index_ptc_aperture !XXXX  !   !
            endif                                                                            !   !
            
            !                                                                                !   !
            if_ptc_track_unstable: IF (flag_index_ptc_aperture==0) then ! =========!         +   ^
              if (rplot) then
                call plottrack(j_particle, my_ring%n, i_th_turn, & 
                                     current_x_coord_incl_co(1), &
                                     current_x_coord_incl_co(2), &
                                     current_x_coord_incl_co(3), &
                                     current_x_coord_incl_co(4), &
                                     current_x_coord_incl_co(5), &
	                 sqrt(MY_RING%end%mag%p%p0c**2 + (MY_RING%end%mass)**2), &
	                 current_x_coord_incl_co(6))
               !   endif
              endif
            ELSE

               if (getdebug()>0) then
               
                 write(6,*) "Track ", particle_ID(j_particle) , " is lost at turn ", i_th_turn, &
		             " flag ", flag_index_ptc_aperture
               

                 if (getdebug() > 1) then

                   write(6,'(a,6(f9.6,1x))') "Track last coordinates ", current_x_coord_incl_co

                   if (flag_index_ptc_aperture < 50) then
                     write(6,*) "         PTC message :", messagelost
                   else
	 write(6,'(a,6(f9.6,1x))') "       lost on maxaper ", maxaper
                   endif

                 endif

               endif

               ! => particle is lost !!(?)                                         
               n_temp=j_last_particle_buffer
               !
                                                              
               CALL kill_ptc_track(n_temp,i_th_turn, zero ,jmax_numb_particl_at_i_th_turn,        &       
                                   particle_ID, last_turn_of_lost_particle,                     &                     
                                   last_position_of_lost_particle, last_orbit_of_lost_particle, & 
                                   x_coord_incl_co,MY_RING%end%mag%name,                        &
	               sqrt(MY_RING%end%mag%p%p0c**2 + (MY_RING%end%mass)**2))
               
               
               
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


    !=============================================================================
    SUBROUTINE track_beam_elementwise_with_PTC ! int.subroutine
      implicit none

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
      character(name_len) :: name_curr_elem
      LOGICAL(lp) ::  Entry_not_exit
      REAL(dp) :: length_curr_elem
      real (dp) :: x_coord_co_temp(1:6) ! buffer for the current values of CO

      LOGICAL :: NaN_coord_after_track_VK=.False. !VK20070328 XXXXXXXXXXXXXXXXXX
      !NaN_coord_after_track_VK=.False.

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
         length_curr_elem=current%MAG%P%ld                                  !                         *
         if (ptc_track_debug) THEN !+++debug print++++++!                   !                         l
            Print *, &                                  !                   !                         o
                 'i_current_elem=', i_current_elem, &   !                   !                         o
                 ' Entry_not_exit=', Entry_not_exit, &  !                   !                         o
                 ' name_curr_elem=', name_curr_elem, &  !                   !                         p
                 ' sum_length=',sum_length, &           !kill_ptc_track                   !                         !
                 ' length_curr_elem=',length_curr_elem  !                   !                         !
            !Print *,'name=',current%MAG%name, &        !                   !                         e
            !        'l=',current%MAG%P%ld              !                   !                         l
         endif !+++end debug print++++++++++++++++++++++!                   !                         e
         Call Particle_Interactions &                                       !                         m
              (i_current_elem, name_curr_elem, Entry_not_Exit,&             !                         e
              sum_length,length_curr_elem)                                  !                         n
         !                                                                  !                         t
         ! at element ENTRY ================================================!                         s
         !                                                                                            !
         
         ! length is needed here to report lost particles
         sum_length=sum_length+current%MAG%P%ld                                                       !
         
         n_temp=1 ! start particle loop with the first not-examined particle                          *
         !                                                                                            *
         Exam_all_particles: DO ! ===<=======<========DO Exam_all_particles ====<============<======! *
            !                                                                                       ! *
            Particle_loop: DO j_th_partic=n_temp, jmax_numb_particl_at_i_th_turn !+++++++++++++++!  ! !
               !                                                                                 !  ! !
               NaN_coord_after_track_VK=.False. !VK20070709 XXXXXXXXXXXXXXXXXXXXXX               !  ! !
               !                                                                                 !  ! !
               jmax_at_loop_start = jmax_numb_particl_at_i_th_turn                               !  ^ !
               j_last_particle_buffer=j_th_partic   ! remember index value after END DO          !  ! !
               !                                                                                 !  ! l
               do k_th_coord=1,6 ! extract coords for the current particle -----!                +  ^ o
                  current_x_coord_incl_co(k_th_coord)= &                        !                +  ! o
                       x_coord_incl_co(k_th_coord,j_th_partic)                  !                +  ! p
               end do !---------------------------------------------------------!                +  ! !
               !                                                                                 +  ^ !
               call track_probe_x(my_ring,current_x_coord_incl_co,MYSTATE, &                     !  ! !
                    fibre1=i_current_elem,fibre2=i_current_elem+1)                               !  ! !
               

               if (ptc_track_debug) then
                 write(mft,*) "t ", i_th_turn, " el ",i_current_elem+1," ",name_curr_elem," track ", j_th_partic, &
                              " : ", current_x_coord_incl_co
               endif  
	

               ! save coordinates for the current particle ---!
               do k_th_coord=1,6 
                  if (ISNAN(current_x_coord_incl_co(k_th_coord))) then 
                     ! BUG !? Aperture does not work, if lattice with spread multipoles XXXX     +  ! !
                     NaN_coord_after_track_VK=.TRUE.                           
                     x_coord_incl_co(k_th_coord,j_th_partic)=999               
                  else  
                     x_coord_incl_co(k_th_coord,j_th_partic)=  &     
                          current_x_coord_incl_co(k_th_coord)        
                  endif 
               end do 


               call PRODUCE_APERTURE_FLAG(flag_index_ptc_aperture)                               !  ! e

               if(flag_index_ptc_aperture/=0) c_%watch_user=.false. !VK20070709 XXXXXXXXXXXXXX   !  ! !
               !                                                                                 !  ! !
               if (ptc_track_debug) then ! debug printing ----------------------------!          !  ^ !
                  !print *, 'PTC: <PRODUCE_APERTURE_FLAG> => flag_index', &           !          +  ! !
                  !                               flag_index_ptc_aperture             !          +  ! !
                  !                                                                   !          +  ! !
                  if(flag_index_ptc_aperture/=0) then !=== print diagnostics ======!  !          !  ! !
                     call ANALYSE_APERTURE_FLAG &                                  !  !          !  ! !
                          (flag_index_ptc_aperture,why_ptc_aperture)               !  !          +  ^ !
                     Write(6,*) "ptc_track unstable (tracking)-programs continues "!  !          +  ! !
                     Write(6,*) "why_ptc_aperture:", why_ptc_aperture              !  !          +  ! !
                  endif !=== print diagnostics ====================================!  !          +  ! !
                  !                                                                   !          +  ! !
               endif ! debug printing ------------------------------------------------!          !  ! !


               if (NaN_coord_after_track_VK) flag_index_ptc_aperture=100 !VK20070328 XXXXXXXXX   !    !
               
               if (ptc_track_debug) then                                                 !XXXX   !    !
                  print *,'flag_index_ptc_aperture is set to', flag_index_ptc_aperture   !XXXX   !    !
               endif                                                                             !    !
               
               if (flag_index_ptc_aperture == 0) then
                 if (abs(current_x_coord_incl_co(1)).ge.maxaper(1) .or. &                          !   !
	 abs(current_x_coord_incl_co(2)).ge.maxaper(2) .or. &                          !   !
	 abs(current_x_coord_incl_co(3)).ge.maxaper(3) .or. &                          !   !
	 abs(current_x_coord_incl_co(4)).ge.maxaper(4)) flag_index_ptc_aperture = 51  !   !

                 if (nvariables .gt. 4) then                            !   !
                   if (abs(current_x_coord_incl_co(5)).ge.maxaper(5)) flag_index_ptc_aperture = 52!   !
                 endif 

                 if ( (nvariables .gt. 5) .and. (MYSTATE%TOTALPATH .eq. 0) ) then                  !   !
                   if (abs(current_x_coord_incl_co(6)).ge.maxaper(6)) flag_index_ptc_aperture = 53!   !
                 endif 
               endif 

               if_ptc_track_unstable: IF (flag_index_ptc_aperture==0) then ! ========!           +  ^ !
                  if (rplot) then
	  call plottrack(j_th_partic, i_current_elem+1, i_th_turn, & 
                                              current_x_coord_incl_co(1), &
                                              current_x_coord_incl_co(2), &
                                              current_x_coord_incl_co(3), &
                                              current_x_coord_incl_co(4), &
                                              current_x_coord_incl_co(5), &
		      sqrt(current%mag%p%p0c**2 + (current%mass)**2), &
		      current_x_coord_incl_co(6))
                  endif
               else
                  ! => particle is lost !!(?)                                        !           +  ! !
                  if (ptc_track_debug) then
                    write(mft,*) "Track ", particle_ID(j_th_partic) , " is lost at ", i_current_elem, " turn ", i_th_turn
                    write(mft,*) "maxaper : ", maxaper
                  endif

                  if (getdebug() > 0) then
                    write(6,*) "Track ", particle_ID(j_th_partic) , " is lost at ", i_current_elem, &
	                                   " turn ", i_th_turn, &
		               " flag ", flag_index_ptc_aperture
	if (getdebug() > 1) then
                      
	  write(6,'(a,6(f9.6,1x))') "Track last coordinates ", current_x_coord_incl_co
	  
	  if (flag_index_ptc_aperture < 50) then
                        write(6,*) "         PTC message :", messagelost
	  else
	    write(6,'(a,6(f9.6,1x))') "       lost on maxaper ", maxaper
	  endif
	
	endif

                  endif


                  
                  
                  n_temp=j_last_particle_buffer                                      !           +  ! l
                  !                                                                                 !           +  ! o
                  CALL kill_ptc_track(n_temp,i_th_turn, sum_length, &                               !           +  ! p
                                      jmax_numb_particl_at_i_th_turn, &                             !           +  ! !
                                      particle_ID, last_turn_of_lost_particle, &                    !           +  ^ !
                                      last_position_of_lost_particle, last_orbit_of_lost_particle, &!           !  ! o
                                      x_coord_incl_co, current%mag%name, &
	                  sqrt(current%mag%p%p0c**2 + (current%mass)**2))                                              !           +  ! v
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
         Call Particle_Interactions &                                       !                         *
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
               
               if (ptc_track_debug) THEN !+++debug print+++++++!                  !    #              *
                  Print *,'obs.No=',number_observation_point   !                  !    #              *
                  Print *,'CO=', x_coord_co_temp               !                  !    #              *
                  Print *,'x_coord_incl_co=', x_coord_incl_co  !                  !    #              *
               endif !+++++++++++++++++++++++++++++++++++++++++!                  !    #              *
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
                       name_curr_elem)                                      !     !    #              *
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
      implicit none

      INTEGER :: j_th_part ! local

      el_name='end             '
      nlm_current_element_number = MY_RING%n

      if (ptc_track_debug) then
         Print *, ' Start SUBR. <Write_tables_after_total_turn>'
         Print *, ' jmax=', jmax_numb_particl_at_i_th_turn
         Print *, ' spos_current_position=', spos_current_position
         Print *, ' nlm_current_element_number=', nlm_current_element_number
         Print *, ' el_name=', el_name
         Print *, ' summ_ring_length=', summ_ring_length
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
            !spos_current_position = zero ! sum_length                 !    !     !       !
            spos_current_position = summ_ring_length !VK 20070401      !    !     !       !
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

      implicit none
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
         Print *, 'The dimensionality of the task: nvariables=',nvariables
      ENDIF  debug_Final_Coord

        
!      particle_ID contains list of particles that survived
      do j_part_tmp=1,jmax_numb_particl_at_i_th_turn             
         last_turn_of_lost_particle(particle_ID(j_part_tmp))= turn_final		                 !
         last_position_of_lost_particle(particle_ID(j_part_tmp))= sum_length              
         ! remember last turn and position of particles          
          last_orbit_of_lost_particle(:, particle_ID(j_part_tmp))= x_coord_incl_co(:, j_part_tmp)
         !                                                       !
      enddo 

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
         call  double_to_table_curr('tracksumm ', 'number ', doublenum)                  !
         !tmp_d = last_turn(i)                                                      !
         doublenum=last_turn_of_lost_particle(j_part_tmp)                           !
         call double_to_table_curr('tracksumm ', 'turn ', doublenum)                     !
         ! call double_to_table_curr('tracksumm ', 'turn ', tmp_d)                       !
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
            call double_to_table_curr('tracksumm ', vec_names(i_coord), doublenum)  !    !
            !call double_to_table_curr('tracksumm ', vec_names(j), tmp_d)           !    !
            !                                                                  !    !
         enddo ! END loop over coord. components >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!    !
         !hbu                                                                       !
         !  spos                  = last_pos(i)                                     !
         spos_current_position = last_position_of_lost_particle(j_part_tmp)         !
         !hbu                                                                       !
         doublenum = spos_current_position  !
         call double_to_table_curr('tracksumm ',vec_names(7),doublenum)                  !
         !                                                                          !
         ! to get "energy" value                                                    !
         Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)        !
         !                                                                          !
         doublenum= energy                !
         call double_to_table_curr('tracksumm ', 'e ',  doublenum)                       !
         !                                                                          !
         call augment_count('tracksumm ')                                           !
      enddo !#### loop over all started particles ##################################!

    END SUBROUTINE Final_Coord_to_tables
    !=============================================================================

    !==============================================================================
    SUBROUTINE Particle_Interactions_Ini
      !USE ptc_track_run_common, ONLY: &
      !            Energy_rest_MeV, Energy_total_MeV
      implicit none

      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet

      Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)
      ! real(dp) ,optional,INTENT(OUT)::MASS,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet

      Energy_rest_MeV=MASS_GeV*c_1D3
      Energy_total_MeV=ENERGY*c_1D3

    END SUBROUTINE Particle_Interactions_Ini
    !==============================================================================

    !==============================================================================
    SUBROUTINE Particle_Interactions &
         (i_current_elem, name_curr_elem, Entry_not_Exit,&
         sum_length,length_curr_elem)
      ! This intenal subroutine of the host subr. <ptc_track_run>
      ! serves as service routine to call any particle interactions
      ! (radiations, space charge and etc.) in a between of the beamline elements)
      !USE madx_ptc_track_run_common, ONLY: &
      !    Energy_rest_MeV, Energy_total_MeV
      ! USE madx_ptc_module, ONLY: dp
      implicit none

      INTEGER, INTENT(IN)       :: i_current_elem
      CHARACTER(name_len), INTENT(IN) :: name_curr_elem !      current%MAG%    name
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
         if (ptc_track_debug) then
             Print *, 'B0_dipole=', B0_dipole, &
              'TiltD_dipole=', TiltD_dipole, &
              '  Quadr_k=',Quadr_k
         endif
         IF (B0_dipole.EQ.zero .AND.Quadr_k .EQ.zero ) i_elem_type=0
         IF (B0_dipole.NE.zero .AND.Quadr_k .EQ.zero ) i_elem_type=1
         IF (B0_dipole.EQ.zero .AND.Quadr_k .NE.zero ) i_elem_type=2
         IF (B0_dipole.NE.zero .AND.Quadr_k .NE.zero ) i_elem_type=3
         if (ptc_track_debug) then
             Print *,'i_elem_type=',i_elem_type
         endif

         IF (i_elem_type .EQ. 0) RETURN

         IF ( (i_elem_type .EQ. 1) .OR. &
              ( (i_elem_type.EQ.3).AND.(.NOT.Radiation_Quad) ) ) THEN
            rad_curv_m = 1D0/B0_dipole
            Call photon (i_elem_type, rad_curv_m, length_curr_elem, &
                                ! Energy_total_MeV, Energy_rest_MeV, ieave,iquasto, d_loss(1), d_loss(2)) !VK20070617
                 Energy_total_MeV,                  ieave,iquasto, d_loss(1), d_loss(2))
            DO j_partic=1, jmax_numb_particl_at_i_th_turn
               IF (Entry_not_Exit) THEN
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(1)
               ELSE
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(2)
               ENDIF
            ENDDO
         END IF

         IF (( i_elem_type .EQ. 2) .AND. Radiation_Quad) THEN
            if (ptc_track_debug) then
                Print *,'jmax_numb_particl_at_i_th_turn=', &
                 jmax_numb_particl_at_i_th_turn
            endif
            DO j_partic=1, jmax_numb_particl_at_i_th_turn
               if (ptc_track_debug) then
                   Print *, 'j_partic=',j_partic
               endif
               SQRT_X2_Y2=SQRT(x_coord_incl_co(1,j_partic)*x_coord_incl_co(1,j_partic)+ &
                    x_coord_incl_co(3,j_partic)*x_coord_incl_co(3,j_partic))
               !IF (SQRT_X2_Y2 .EQ. zero) EXIT
               rad_curv_m = 1D0/Quadr_k/SQRT_X2_Y2

               Call photon (i_elem_type, rad_curv_m, length_curr_elem, &
                                !Energy_total_MeV, Energy_rest_MeV, ieave,iquasto, d_loss(1), d_loss(2)) !VK20070617
                    Energy_total_MeV,                  ieave,iquasto, d_loss(1), d_loss(2))
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
                                !Energy_total_MeV, Energy_rest_MeV, ieave,iquasto, d_loss(1), d_loss(2)) !VK20070617
                    Energy_total_MeV,                  ieave,iquasto, d_loss(1), d_loss(2))
               IF (Entry_not_Exit) THEN
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(1)
               ELSE
                  x_coord_incl_co(5,j_partic)=x_coord_incl_co(5,j_partic)+d_loss(2)
               ENDIF
            ENDDO
         END IF

         ! Call photon (i_Mag_type, rad_curv_m, dl_length, Energy_beam_MeV,
         !                                                           ieave,iquasto,d1, d2)
         ! Radiation_Energy_Loss=TRUE => ieave=1; Radiation_Energy_Loss=FALSE => ieave=0;
         ! Radiation_Quad     =TRUE => iquasto=1; Radiation_Quad     =FALSE => iquasto=0
      END IF Radiation_by_FZ_code


      Space_Charge_Calculation: IF ( Space_Charge ) THEN
         ! Call Space_Charge
         if (ptc_track_debug) then
             Print *, '  i_current_elem=',i_current_elem, &
              '  name_curr_elem=', name_curr_elem, &
              '  sum_length=', sum_length
         endif

      ENDIF Space_Charge_Calculation

      ! print * , ' SUBROUTINE Particle_Interactions is under construction => RETURN'
      RETURN
    END SUBROUTINE Particle_Interactions
    !==============================================================================

    !==============================================================================
    SUBROUTINE Prepare_Observation_points (max_obs, x_coord_co_at_START)
      ! getting parameters  at observation points and
      ! finding the closed orbits at observations for element_by_element tracking

      ! USE  madx_ptc_module, ONLY: dp, my_ring, track, default
      implicit none

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

      character(name_len) ::      local_name

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

         find_CO_for_el_by_el: IF (element_by_element) THEN            !===!
            IF (closed_orbit) THEN !-----------------------------------!  !
               call track_probe_x(my_ring,x_coord_co_temp,MYSTATE, &   !  ! 
                    fibre1=i_ring_element,fibre2=i_ring_element+1)     !  ! 
!               Call track(my_ring,x_coord_co_temp, &                  !  !
!                    i_ring_element, i_ring_element+1, MYSTATE )       !  !
            ELSE                                                       !  !
               x_coord_co_temp(:)=zero                                 !  !
            ENDIF !----------------------------------------------------!  !
         ENDIF find_CO_for_el_by_el !=====================================!

         IF(number_obs.GT.0) THEN
            elem_number_at_observ(number_obs)= i_ring_element
            sum_length_at_observ(number_obs) = Sum_length_S
            call element_name(local_name,name_len)
            name_el_at_obsrv(number_obs) = local_name

            save_CO_for_el_by_el: IF (element_by_element) THEN
               DO i_coord=1,nvariables
                  x_co_at_all_observ(i_coord,number_obs)=x_coord_co_temp(i_coord)
               ENDDO
            ENDIF save_CO_for_el_by_el

         ENDIF
         if (ptc_track_debug) then
            Print *, 'i_el ', i_ring_element,' num_obs=', number_obs
            Print *, ' l_c_code=',length_current_element_c,' l_f90=', &
                 length_current_element_f90, &
                 ' name_c=', local_name, ' &_f90=', current%MAG%name
            Print *, 'x_coord_co_temp=', x_coord_co_temp

         endif

         if (i_ring_element .EQ. my_ring%n) summ_ring_length = Sum_length_S

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
      implicit none

      integer, intent (IN) :: max_obs ! the maximum number of observation points >=1
      ! one point at the end (beginning) plus
      ! the points given in input file by the command
      ! "ptc_observe,place=mark";
      REAL (dp),     INTENT( IN) :: x_coord_co_at_START(1:6) ! => x0(1:6) in ptc_track at START
      TYPE (real_8), intent (INOUT) :: Map_Y_obs(6) !  y => Map_Y_obs - local name
      TYPE(damap) :: Map_damap
      INTEGER ::  iii_c_code, i_obs_point, i_coord, &
                  i_from, i_till, i_dummy, i_turn_tmp, j_part_tmp, ielem
      INTEGER ::  flag_index,why(9)
      INTEGER, ALLOCATABLE :: Temp_particle_ID(:)
      integer  :: mf1,mf2,mf(max_obs)
      REAL(dp) :: spos
      real(dp) :: X_co_temp (6) ! (36) !(6) !(lnv)
      real(dp), allocatable ::  X_co_observe(:,:) ! (i_coord, i_obs_point)
      integer,  allocatable :: J(:) ! for extracting CO from the map
      REAL(dp),  allocatable :: Temp_X_incl_co_at_obs(:,:)
      Real(dp):: X_lnv_START(lnv), X_lnv_OBSRV(lnv)
      character(8) ch,ch1

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
         call kanalnummer(mf1)
         call kanalnummer(mf2)
         open(unit=mf1,file='map_y_obs_in.txt')
         open(unit=mf2,file='map_y_obs.txt')
         write(mf1,*) 'Map_Y_obs IN: '; call daprint(Map_Y_obs,mf1);
      endif
      Map_Y_obs=npara
      Map_Y_obs=x_coord_co_at_START  ! Y=X

      if (ptc_track_debug) then
         Print *, ' x_coord_co_at_START=', x_coord_co_at_START
         write(mf2,*) 'Map_Y_obs=x_coord_co_at_START: '; call daprint(Map_Y_obs,mf2);
      endif

      Allocate (X_co_observe(1:6,1:max_obs), Temp_particle_ID(1:j_tot_numb_starting_particles))
      X_co_observe(:,:)=zero; 
      Temp_particle_ID=0

      DO i_coord=1,nvariables
         X_co_observe(i_coord,1)=x_coord_co_at_START(i_coord)
      END DO 

      CALL Alloc(Map_damap)

      iii_c_code=restart_sequ()  ! c-code restart the beamline sequence
      current=>MY_RING%start     ! F90 pointer to the CURRENT beamline element is set up

      obs_point_loop: DO i_obs_point=1, max_obs-1

         if (ptc_track_debug) then
            call kanalnummer(mf(i_obs_point))
            ch=" "
            ch1=" "
            write(ch,'(i8)') i_obs_point
            ch1=adjustl(ch)
            open(unit=mf(i_obs_point),file="obs_point_"//ch1(:len_trim(ch1)))
         endif

         i_from=elem_number_at_observ(i_obs_point)
         i_till=elem_number_at_observ(i_obs_point+1)

         if (ptc_track_debug) then
            Print *, 'i_obs_point=', i_obs_point, ' name_f90=', current%MAG%name
            Print *, 'X_co_observe(:,i_obs_point) =', &
                 (X_co_observe(i_coord,i_obs_point), i_coord=1,6)
            Print *, 'Track from i_from=', i_from, 'i_till =',i_till
         end if

         call track_probe_x(my_ring,Map_Y_obs,MYSTATE,fibre1=i_from,fibre2=i_till)
         !call track(my_ring,Map_Y_obs,1,MYSTATE)


         Map_damap=Map_Y_obs
         

         if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
            write(whymsg,*) 'DA got unstable. ', &
	        ', PTC msg: ',messagelost(1:LEN_TRIM(messagelost))

            call fort_warn('ptc_track: ',whymsg(1:LEN_TRIM(whymsg))) 
            whymsg(LEN_TRIM(whymsg)+1:LEN_TRIM(whymsg)+1) = char(0)
            call seterrorflag(10,"ptc_trackline ",whymsg);

          endif

         
         
         if (ptc_track_debug) then
            write(mf(i_obs_point),*) 'i_Unit=', mf(i_obs_point); Call daprint(Map_Y_obs,mf(i_obs_point));
         endif
         call PRODUCE_APERTURE_FLAG(flag_index)
         if(flag_index/=0) then
            call ANALYSE_APERTURE_FLAG(flag_index,why)
            Write(6,*) "Get_map_from_NormalForm is unstable (map production)-programs continues "
            Write(6,*) why ! See produce aperture flag routine in sd_frame
            CALL kill(Map_Y_obs) !CALL kill(y)
            ! c_%watch_user=.false.
            if (ptc_track_debug) then
               close(mf1)
               close(mf2)
            endif
            close(mf(i_obs_point))
            return
         endif

         ielem=elem_number_at_observ(i_obs_point+1)

         ! Call Extract_CO_from_Map (icase_PTC, Map_Y_obs, X_co_temp) - subr. removed
         allocate(J(nvariables)); J(:)=0
         X_co_temp(:)=zero
         DO i_coord=1, nvariables
            X_co_temp(i_coord) = (Map_Y_obs(i_coord)%T.sub.J) ! take line with all zero-order
            if (ptc_track_debug) then
                Print *,'CO extracted: i_coord=', i_coord, &
                 ' X_co_temp(i_coord)=', X_co_temp(i_coord)
            endif
         ENDDO
         deallocate(J)
         Save_co_in_X_co_observe: DO i_coord=1, nvariables
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
               if (ptc_track_debug) then
                   Print *, '   the current particle j_part_tmp= ', j_part_tmp
               endif



               !X_lnv_OBSRV=Map_Y_obs*X_lnv_START
               !             X_lnv_OBSRV=Map_Y_obs%T*X_lnv_START
               !Error: Operands of binary numeric operator '*' at (1) are TYPE(taylor)/real(dp)
               Read_START_coord: DO i_coord=1,nvariables
                  X_lnv_START(i_coord)= x_all_incl_co_at0(i_coord, i_turn_tmp, j_part_tmp) - &
                       x_coord_co_at_START(i_coord) ! X_out=M(x0,x)*x_in, where X=x0+x !
               ENDDO Read_START_coord
               ! Transform from START to observation point
               X_lnv_OBSRV=Map_damap*X_lnv_START
               !X_lnv_OBSRV=Map_Y_obs*X_lnv_START
               !             X_lnv_OBSRV=Map_Y_obs%T*X_lnv_START
               !Error: Operands of binary numeric operator '*' at (1) are TYPE(taylor)/real(dp)

               Loop_coord: DO i_coord=1,nvariables
                  !X_lnv_OBSRV(i_coord)=Map_Y_obs(i_coord)%T*X_lnv_START(i_coord)
                  !X_lnv_OBSRV(i_coord)=Map_Y_obs(i_coord)*X_lnv_START(i_coord)
                  !X_lnv_OBSRV(i_coord)=Map_damap*X_lnv_START(i_coord)
                  if (ptc_track_debug) then
                      Print *, 'i_coord=',i_coord, &
                       'X_lnv_START/OBSRV=', X_lnv_START(i_coord), X_lnv_OBSRV(i_coord)
                  endif
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
         if (ptc_track_debug) close(mf(i_obs_point))
      END DO obs_point_loop !####################################################################!

      DEALLOCATE( X_co_observe,  Temp_particle_ID)

      CALL kill(Map_damap)
      if (ptc_track_debug) then
         close(mf1)
         close(mf2)
      endif
      !Call kill(Map_Y_obs)
    END SUBROUTINE Observation_with_PTC
    !==============================================================================

    !==============================================================================
    SUBROUTINE tt_putone_coord (npart,turn,tot_segm,segment,part_id,z,orbit0,&
         &spos,ielem,el_name)
      ! copied from TRRUN.F(tt_putone) and renamed to tt_putone_coord
      !hbu added spos, ielem, el_name

      ! USE  madx_ptc_module, ONLY: dp,lnv,zero, operator(*), assignment(=)
      implicit none
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
      integer i,j,npart,turn,tot_segm,segment,part_id(*),length
      REAL(dp) :: z(6,*),orbit0(6),tmp !, tt, ss
      real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet
      !vvk
      REAL(dp) :: tmp_coord_array(lnv), tmp_norm_array(lnv) ! tmp_norm

      REAL (dp) :: X_MAD(6), X_PTC(6)

      !hbu was *36 allow longer info
      character(80) table_putone
      character(41+name_len) comment
      !hbu
      integer :: ielem
      !hbu name of element
      character(name_len) :: el_name
      !hbu
      REAL(dp) :: spos
      !hbu
      character(4) vec_names(7)
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
      call comment_to_table_curr(table_putone, comment, length)
      Call GET_ONE(MASS_GeV,energy,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet) ! to get "energy" value
      do i = 1, npart
         doublenum=turn
         call double_to_table_curr(table_putone, 'turn ', doublenum)
         doublenum=energy
         call double_to_table_curr(table_putone, 'e ', doublenum)
         doublenum = part_id(i)
         call double_to_table_curr(table_putone, 'number ', doublenum)

         do j = 1, 6
            tmp=zero
            tmp = z(j,i) - orbit0(j)
            tmp_coord_array(j)=tmp ! make array of coordinates
            X_PTC(j)=tmp
         end do
         CALL Coord_PTC_to_MAD(X_PTC,X_MAD) ! convert coordinates
         do j = 1, 6
            doublenum=X_MAD(j)
            IF( (.NOT.closed_orbit) .OR. (.NOT.NORM_OUT)) THEN
               call double_to_table_curr(table_putone, vec_names(j), doublenum)
            END IF
         enddo

         ! make normalization, if necessary
         IF( closed_orbit .AND. NORM_OUT) THEN
            tmp_norm_array=A_t_map_rev*tmp_coord_array

            DO j = 1, 6
               X_PTC(j)=zero
               IF (j.LE.nvariables) X_PTC(j)=tmp_norm_array(j)
            ENDDO
            CALL Coord_PTC_to_MAD(X_PTC,X_MAD) ! convert coordinates
            DO j = 1, 6
               !tmp_norm=zero
               doublenum=X_MAD(j)
               call double_to_table_curr(table_putone, vec_names(j), doublenum)
            END DO
         END IF

         !hbu spos
         doublenum=spos
         call double_to_table_curr(table_putone,vec_names(7),doublenum)
         call augment_count(table_putone)
      enddo
    END SUBROUTINE tt_putone_coord
    !==============================================================================

    !==============================================================================
    SUBROUTINE tt_puttab_coord (npart,turn,nobs,orbit,orbit0,spos)
      ! save to table of this particular track (trackone == false)
      ! copied from TRRUN.F( tt_puttab) and renamed to tt_puttab_coord
      !
      ! USE  madx_ptc_module, ONLY: dp, lnv, zero, operator(*), assignment(=)
      implicit none
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

      character(36) table_puttab
      !hbu
      REAL(dp) :: spos
      !hbu
      character(4) vec_names(7)
      !hbu
      data vec_names  / 'x', 'px', 'y', 'py', 't', 'pt','s' / ! MAD order
      !data vec_names / 'x', 'px', 'y', 'py', 'pt', 't','s' / ! PTC has a reverse order for pt and t
      data table_puttab / 'track.obs$$$$.p$$$$' /

      !print*,"skowron: tt_puttab_coord"
      tmp_coord_array=zero; tmp_norm_array=zero

      !tt = turn; !tn = npart
      write(table_puttab(10:13), '(i4.4)') nobs    ! Write in the table head :
      write(table_puttab(16:19), '(i4.4)') npart   ! "@NAME ... "TRACK.OBS0001.P0005"

      Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet) ! to get "energy" value
      doublenum=energy
      call double_to_table_curr(table_puttab, 'e ', doublenum)
      doublenum=npart
      call double_to_table_curr(table_puttab, 'number ', doublenum) ! number of the current particle
      doublenum=turn
      call double_to_table_curr(table_puttab, 'turn ', doublenum)   ! the number of the current turn
      do j = 1, 6
         tmp=zero
         tmp = orbit(j) - orbit0(j)
         tmp_coord_array(j)=tmp ! make array of coordinates
         X_PTC(j)=tmp
      end do

      IF( closed_orbit .AND. NORM_OUT) THEN
         tmp_norm_array=A_t_map_rev*tmp_coord_array
         DO j = 1, 6
            X_PTC(j)=zero
            IF (j.LE.nvariables) X_PTC(j)=tmp_norm_array(j)
         ENDDO

      ENDIF
      CALL Coord_PTC_to_MAD(X_PTC,X_MAD) ! convert coordinates (canonical or normal)
      do j = 1, 6
         IF (closed_orbit .AND. NORM_OUT) THEN
            !tmp_norm=zero
            doublenum=X_MAD(j)
            call double_to_table_curr(table_puttab, vec_names(j), doublenum)
         ELSE
            !tmp=zero
            doublenum = X_MAD(j) ! orbit(j) - orbit0(j)
            call double_to_table_curr(table_puttab, vec_names(j), doublenum)
         END IF
      enddo
      !hbu spos
      doublenum=spos
      call double_to_table_curr(table_puttab,vec_names(7),doublenum)
      call augment_count(table_puttab)
    END SUBROUTINE tt_puttab_coord
    !==============================================================================
    

!Do not know how to implement it with new PTC, need to ask Etienne
    !==============================================================================
!     SUBROUTINE beam_enevelope_with_PTC
!       use ptc_spin
!       use madx_ptc_module, ONLY:ENVELOPE0,track_probe,track_probe_x,get_loss,&
!            FIND_ORBIT_x,init,print,alloc,kill
!       implicit none
!       TYPE (REAL_8) :: Y(6)
!       TYPE(INTERNAL_STATE) MYSTATE_ENV
!       type(damap) id
!       type(normalform) normal
!       TYPE (PROBE) XS0
!       TYPE (PROBE_8) XS
!       TYPE (DAMAPSPIN) M
!       TYPE (NORMAL_SPIN) nf
!       real(dp) x(6),energy,deltap
! 
!       CALL INIT(MYSTATE,1,0)
!       call alloc(id)
!       call alloc(y)
!       call alloc(normal)
!       !nf%stochastic=my_true
! 
!       x=zero
!       CALL FIND_ORBIT_x(my_ring,X,MYSTATE,1.0e-7_dp,fibre1=1)
!       WRITE(6,'(A)') " Closed orbit with Radiation "
!       WRITE(6,'(6(1x,E15.8))') x
!       call GET_loss(my_ring,energy,deltap)
!       write(6,'(a32,2(1x,E15.8))') "Energy loss: GEV and DeltaP/p0c ",energy,deltap
! 
!       id=1
!       y=x+id
! 
!       call TRACK_PROBE_X(my_ring,y,MYSTATE, FIBRE1=1)
! 
!       id=y
!       normal=id   ! Normal is a regular Normalform type
!       write(6,*) "Tunes "
!       write(6,*) normal%tune(1:3)
!       write(6,*) "Damping decrements"
!       write(6,*) normal%damping(1:3)
! 
!       call kill(id)
!       call kill(y)
!       call kill(normal)
! 
!       MYSTATE_ENV=MYSTATE+ENVELOPE0
!       CALL INIT(mystate_env,2,0)
!       call alloc(m)
!       call alloc(xs)
!       call alloc(nf)
!       
!       if (getdebug() > 1) then
!         PRINT*," "
!         WRITE(6,'(A)') " Print the state: MYSTATE_ENV "
!         PRINT*," "
!         call print(MYSTATE_ENV,6)
!       endif
!       
!       xs0=x
!       m=1        ! damapspin set to identity
!       xs=xs0+m   ! Probe_8 = closed orbit probe + Identity
! 
!       call track_probe(my_ring,xs,MYSTATE_ENV,fibre1=1)
! 
!       m=xs       ! damapspin = Probe_8 
!       nf=m       ! normal_spin = damapspin (Normalization including spin (if present) or radiation
!       ! envelope if present. (Spin without radiation)
!       write(6,*) ' Tunes : '
!       write(6,*) nf%n%tune(1:3)
!       write(6,*) ' The equilibrium emittances are : '
!       write(6,*) nf%emittance
!       write(6,*) ' nf%s_ij0(1,1),nf%s_ij0(5,5) : '
!       write(6,*) nf%s_ij0(1,1),nf%s_ij0(5,5)
! 
!       call kill(m)
!       call kill(xs)
!       call kill(nf)
! 
!     END SUBROUTINE beam_enevelope_with_PTC
!     !==============================================================================


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
      implicit none

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

      INTEGER ::  j_particle_line_counter,kq,kp
      INTEGER ::  next_start ! int. function
      REAL(KIND(1d0))  :: x_input,px_input,y_input,py_input,t_input,deltae_input
      REAL(KIND(1d0))  :: fx_input,phix_input,fy_input,phiy_input,ft_input,phit_input
      REAL(dp)  :: phi_n
      INTEGER :: k_th_coord,  j_th_particle ! local
      INTEGER :: itype_non_zero_flag(12) ! not itype(23)
      real(dp) :: track_temp(12), z_start_sum(6) ! =1 for non-zero coordinates, =0 for zero
      real(dp),allocatable :: x_coord_incl_co_temp(:,:) ! copy of array
      Real (dp) :: Z_norm_temp(lnv), Z_stdt_temp(lnv) ! Z_stdt=A_t_map*Z_norm
      REAL (dp) :: X_MAD(6), X_PTC(6)
      logical(lp) :: zgiv_exist, zngiv_exist ! existence of non-zero input for action-angle
      character(120) msg(2) ! text stings for messages

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

      j_particle_line_counter  =  &                                                   ! o
           next_start(x_input,px_input,y_input,py_input,t_input,deltae_input, &       ! r
                      fx_input,phix_input,fy_input,phiy_input,ft_input,phit_input)    !

      if (ptc_track_debug) then                                          !
         print *, "The next command line from input file is read."                    ! p
         print *, "The keyword <ptc_start> provides the following data:"              ! a
         print *, "j_particle_line_counter=",j_particle_line_counter                  ! r
         print *, "x,px,y,py,t,deltae,fx,phix,fy,phiy,ft,phit=", &                    ! t
              x_input,px_input,y_input,py_input,t_input,deltae_input, &               ! i
              fx_input,phix_input,fy_input,phiy_input,ft_input,phit_input             ! c
      endif                                                             ! l
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

            if (getdebug() > 2) then
               print*, 'X_MAD: ', X_MAD
            endif
            
            ! Very dodgy
            ! add momentum offset defined by 'deltap' to all tracks 
            IF(nvariables.gt.4 .AND. (.NOT.closed_orbit)) THEN !--!            ! h   i   o
               if(mytime) then !----------------------!          !              !
                  call Convert_dp_to_dt (deltap, dt)  !          !              !
               else                                   !          !              !
                  dt=deltap                           !          !              !
               endif !--------------------------------!          !              !
               X_MAD(6)=X_MAD(6)+dt                              !              !   l   !
            ENDIF !----------------------------------------------!              !   e   !
            
            !                                                                   !   !   p
            CALL Coord_MAD_to_PTC(X_MAD,X_PTC) ! convert coordinates          ! c   t   a
            
            if (getdebug() > 2) then
               print*, 'X_PTC: ', X_PTC
            endif
            
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
            polar_cart: do kq = 1, nvariables-1, 2 ! --------------!             !   c   !
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
               Print *,'Z_norm_temp(k_th_coord)=', Z_norm_temp(:nvariables) !    !   r   !
            end if debug_print_5 !-----------------------------------------!    !   t   !
            !                                                                   !   i   !
            !-- Transform to unnormalized coordinates and refer to closed orbit.!   c   !
            Z_stdt_temp(:)=zero;                                                !   l   !
            !                                                                   !   e   !
            if(closed_orbit) Z_stdt_temp=A_t_map*Z_norm_temp                    !   !   !
            !                     ! ANY ORDER !!!                               !   !   !
            debug_print_6: if (ptc_track_debug) then ! --------!                !   !   !
               Print *, 'Z_stdt_temp=A_t_map*Z_norm_temp=', &  !                !   !   !
                    Z_stdt_temp(:nvariables)                    !                !   !   !
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
            SumDO: DO k_th_coord = 1,6 ! not 6 !----------------------! !   !   !
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
            IF (jmax_numb_particl_at_i_th_turn.GE.2) THEN !*********! !   !   !

               x_coord_incl_co(:,1:jmax_numb_particl_at_i_th_turn-1) = & 
                    x_coord_incl_co_temp(:,1:jmax_numb_particl_at_i_th_turn-1)
               
               
               DEALLOCATE (x_coord_incl_co_temp)                              ! !   !   !
            END IF 
            !                                                                   !   !   !
            ! fill current particle                                             !   !   !
            x_coord_incl_co(:,jmax_numb_particl_at_i_th_turn)=x_coord_co(:)+z_start_sum(:) 
            !                                                                   !   !   !
            ! save already-read particles in the array                          !   !   !
            
            ALLOCATE(x_coord_incl_co_temp(1:6,1:jmax_numb_particl_at_i_th_turn))!   e   !
            x_coord_incl_co_temp=zero                                           !   w   !
            x_coord_incl_co_temp(:,:)= x_coord_incl_co(:,:)    
            DEALLOCATE (x_coord_incl_co)             !   t   !
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
            DO k_th_coord = 1, nvariables   !------------------------!  !            !
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
      ! only swap 5 with 6
      ! and than sign of 6
      IMPLICIT NONE
      REAL(dp), INTENT(IN)  :: X_MAD(6)
      REAL(dp), INTENT(OUT) :: X_PTC(6)

      X_PTC=X_MAD ! for all elements

      if (ptc_track_debug) then
          print *,'Coord_MAD_to_PTC icase_ptc=', icase_ptc, ' mytimec=', mytime
          
      endif

      !write(6,'(a10,1x,6(f9.6,1x))') 'X_MAD  ',X_MAD

      if (nvariables.gt.5) THEN ! 6 and 56
         
         X_PTC(5)=X_MAD(6); 
         X_PTC(6)=X_MAD(5);
         X_PTC(6)=-X_PTC(6) ! reverse sign
         
      elseif(nvariables.eq.5) THEN
      
         X_PTC(5)=X_MAD(6); 
         X_PTC(6)=zero
      
      elseif(nvariables.eq.4) THEN
         X_PTC(5)=zero; 
         X_PTC(6)=zero
      ENDIF

      !write(6,'(a10,1x,6(f9.6,1x))') 'X_PTC  ',X_PTC
      
    END SUBROUTINE Coord_MAD_to_PTC

    SUBROUTINE Coord_PTC_to_MAD(X_PTC,X_MAD)
      ! Convert coordinates from PTC to MAD
      IMPLICIT NONE
      REAL(dp), INTENT(IN)  :: X_PTC(6)
      REAL(dp), INTENT(OUT) :: X_MAD(6)

      X_MAD(1:4)=X_PTC(1:4) ! for transverse

      if (ptc_track_debug) then
          print *, &
           'Coord_PTC_to_MAD icase_ptc=', icase_ptc, ' mytime=', mytime
      endif

      IF (nvariables.gt.5) THEN
         X_MAD(5)=X_PTC(6); 
         X_MAD(6)=X_PTC(5);
         X_MAD(5)=-X_MAD(5) ! reverse sign
         
      elseif(nvariables.eq.5) THEN
      
         X_MAD(5)=zero; 
         X_MAD(6)=X_PTC(5)
      
      elseif(nvariables.eq.4) THEN
         X_MAD(5)=zero; 
         X_MAD(6)=zero
      ENDIF
      
    END SUBROUTINE Coord_PTC_to_MAD
    !=============================================================================

  END subroutine ptc_track_run
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
  !==============================================================================
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
         track, track_probe_x, init, alloc, kill, &
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
    INTEGER   :: i_vec, i_comp, ind(6), Number_of_eigenvectors
    integer   :: mf1,mf2,mf3
    !------------------------------------------------------------------------------
    IF (ptc_track_debug)  print *, 'Start Subr.  Get_map_from_NormalForm '

    nda=0
    ! no = 1 =>  Normal_Order_n0 (INPUT)

    call init(MYSTATE,Normal_Order_n0,nda,BERZ,mynd2,npara)
    if (ptc_track_debug) then
       Print *, 'aftercall init(MYSTATE,Normal_Order_n0,nda,BERZ,mynd2,npara)'
       Print *, 'Normal_Order_n0,nda,mynd2,npara=', Normal_Order_n0,nda,mynd2,npara
    endif

    call alloc(Map_Y) ! call alloc(y)

    Map_Y=npara ! y=npara
    Map_Y=x_coord_co  ! Y=X
    IF (ptc_track_debug) then
       call kanalnummer(mf1)
       call kanalnummer(mf2)
       open(unit=mf1,file='map_y_ini.txt')
       open(unit=mf2,file='map_y.txt')
       call daprint(Map_Y,mf1);
    endif
    !    c_%watch_user=.true.

    call track_probe_x(my_ring,Map_Y,MYSTATE,fibre1=1)
    IF (ptc_track_debug)call daprint(Map_Y,mf2);


    call PRODUCE_APERTURE_FLAG(flag_index)
    if(flag_index/=0) then
       call ANALYSE_APERTURE_FLAG(flag_index,why)
       Write(6,*) "Get_map_from_NormalForm is unstable (map production)-programs continues "
       Write(6,*) why ! See produce aperture flag routine in sd_frame
       CALL kill(Map_Y) !CALL kill(y)
       ! c_%watch_user=.false.
       IF (ptc_track_debug) then
          close(mf1)
          close(mf1)
       endif
       return
    endif
    !    c_%watch_user=.false.

    call alloc(Normal_Form_N) ! call alloc(n)
    ! Normal_Form_N=0 !Error: Can't convert INTEGER(4) to TYPE(normalform)

    Normal_Form_N=Map_Y ! n=y
    
    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then

       write(whymsg,*) 'DA got unstable when calculating the closed solution ', &
	   ', PTC msg: ',messagelost(1:LEN_TRIM(messagelost))

       call fort_warn('ptc_track: ',whymsg(1:LEN_TRIM(whymsg)))
       whymsg(LEN_TRIM(whymsg)+1:LEN_TRIM(whymsg)+1) = char(0)
       call seterrorflag(10,"ptc_track ",whymsg);
       call aafail("ptc_track","Fatal Error: Requested closed orbit but closed solution does not exist! program stops")
    endif

    call alloc(A_t_map)

    A_t_map=Normal_Form_N%A_t

    call alloc(A_t_map_rev)

    A_t_map_rev=Normal_Form_N%A_t**(-1)

    if (ptc_track_debug) then
       call kanalnummer(mf3)
       OPEN (UNIT=mf3,FILE='normal_form_a.txt', STATUS='UNKNOWN')
       write(mf3,'(/a/)') 'Map_Y: '; call daprint(Map_Y,mf3);
       write(mf3,'(/a/)') 'Linear Normal_Form_N%A: '; call daprint(Normal_Form_N%A,mf3);     !
       write(mf3,'(/a/)') 'Linear Normal_Form_N%A_t: '; call daprint(Normal_Form_N%A_t,mf3); !
       write(mf3,'(/a/)') 'Linear A_t_map: '; call daprint(A_t_map,mf3);                     !
       write(mf3,'(/a/)') 'Linear A_t_map_rev: '; call daprint(A_t_map_rev,mf3);             !
       Close (mf3)
       ! CALL TEST_PTC_Normal(Normal_Form_N)
    end if

    !EigenVectors

    Number_of_eigenvectors=nvariables
    if (icase_ptc .eq. 5)  Number_of_eigenvectors=4
    if (icase_ptc .eq. 56) Number_of_eigenvectors=4
    
    do i_vec=1,Number_of_eigenvectors 
       do i_comp=1,nvariables 
          ind(:)=0; ind(i_comp)=1
          d_val(i_comp)=Normal_Form_N%A_t%V(i_vec).sub.ind
       enddo
       if (ptc_track_debug) then
          WRITE(mf2,*) 'EigenVector V(',i_vec,')=', d_val
       endif
       !enddo
    end do
    !DEallocate (d_val)
    IF (ptc_track_debug) then
       close(mf1)
       close(mf2)
    endif
  END subroutine Get_map_from_NormalForm
  !==============================================================================

    SUBROUTINE kill_ptc_track (n_particle,       &
                               i_th_turn,        &
	           sum_accum_length, &
	           j_max,            &
                               part_ID,          &
	           last_turn,        &
	           last_position_of, &
	           last_orbit_of,    &
                               x_coord_incl_co, el_name,en )
      implicit none
      ! the next variables are local
      Integer, intent (IN)     :: n_particle  ! number of particle to be kill
      Integer, intent (IN)     :: i_th_turn
      real(dp),intent (IN)     :: sum_accum_length, en ! input as 0
      Integer, intent (INOUT)  :: j_max    ! number surviving particles (here, j=j-1)
      Integer, intent (INOUT)  :: part_ID(*) ! particle number is reodered(1:N_particle_max)
      Integer, intent (OUT)    :: last_turn(*) ! last_turn_of_lost_particle (1:N_particle_max)
      real(dp), intent (OUT)   :: last_position_of(*), last_orbit_of(1:6,*)
      real(dp), intent (INOUT) :: x_coord_incl_co(1:6,*)
      character(len=24), intent (IN) :: el_name ! PTC has longer names
      character(len=name_len) :: madx_name
      Integer :: j_th, k_th, np, namelen
      
      np = part_ID(n_particle)
      
      last_turn(np)=i_th_turn
      ! Save Number of this turn for n_particle

      last_position_of(np)=sum_accum_length
      ! Save Position of n_particle

      last_orbit_of(:,np) = x_coord_incl_co(:,np)

      if (recloss) then
        namelen = LEN_TRIM(el_name)
        if (namelen > name_len) namelen = name_len
        madx_name = el_name(1:namelen)
        call tp_ploss(np,i_th_turn, sum_accum_length, x_coord_incl_co(:,np), madx_name, en)
      endif 
      
      ! Renumbering arrays
      do j_th = n_particle+1, j_max 

         part_ID(j_th-1) = part_ID(j_th)              
         x_coord_incl_co(:,j_th-1) = x_coord_incl_co(:,j_th)  
         
         IF (.NOT.element_by_element) THEN 
            part_ID_turns(i_th_turn,j_th-1) = part_ID(j_th-1) 
         END IF
         
      enddo 

      j_max = j_max - 1

    END SUBROUTINE kill_ptc_track
    !=============================================================================


subroutine tp_ploss(npart,turn,spos,orbit,el_name, energy)
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !--- purpose: enter lost particle coordinates in table                 *
  !    input:                                                            *
  !    npart  (int)           particle number                            *
  !    turn   (int)           turn number                                *
  !    spos    (double)       s-coordinate when loss happens             *
  !    orbit  (double array)  particle orbit                             *
  !    orbit0 (double array)  reference orbit                            *
  !----------------------------------------------------------------------*
  integer :: npart, turn
  double precision :: spos, orbit(6)
  character(len=name_len) :: el_name

  integer :: j
  double precision :: tmp, tt, tn, energy
  character(len=120) :: table='trackloss'
  character(len=4) :: vec_names(6)
  data vec_names / 'x', 'px', 'y', 'py', 'pt', 't'/

  double precision, external :: get_value

  tn = npart
  tt = turn

  ! energy = get_value('probe ','energy ')

  ! the number of the current particle
  call double_to_table_curr(table, 'number ', tn)
  ! the number of the current turn
  call double_to_table_curr(table, 'turn ', tt)


  do j = 1, 6
     tmp = orbit(j)
     call double_to_table_curr(table, vec_names(j), tmp)
  enddo

  tmp = spos
  call double_to_table_curr(table,'s ',tmp)

  call double_to_table_curr(table, 'e ', energy)
  call string_to_table_curr(table, 'element ', el_name)

  call augment_count(table)
end subroutine tp_ploss



END MODULE madx_ptc_track_run_module
!==============================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

