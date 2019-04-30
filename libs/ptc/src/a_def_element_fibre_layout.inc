!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
! The concept of the fibre and the concept of the
! integration node, with their resulting link list types,
! the layout and the node_layout,
! are based on concepts first ellaborated  with J. Bengtsson.
! The node_layout is similar to the Lagrangian class that
! Bengtsson and Forest contemplated around 1990 for the
! C++ collaboration later known as CLASSIC.


TYPE girder
   INTEGER, pointer :: i  => null()
   INTEGER, pointer :: pos => null()
   INTEGER, pointer :: mul  => null()   ! number of magnets in the girder
   INTEGER, pointer :: added  => null()
   INTEGER, pointer :: discarded  => null()
   real(dp), pointer:: a(:) => null(),ent(:,:) => null()
   TYPE (girder_info),POINTER :: info  => null()
   TYPE (girder),POINTER :: PREVIOUS  => null()  ! Terminated link list
   TYPE (girder),POINTER :: NEXT  => null()
END TYPE girder

TYPE girder_info
   CHARACTER(nlp), POINTER ::  NAME  => null()   ! Identification
   real(dp), POINTER ::  d(:)  => null()
   type (element), pointer :: mag  => null()
   real(dp), pointer:: a(:) => null(),ent(:,:) => null()
   type (girder), pointer :: parent_girder => null()
   TYPE (girder_info),POINTER :: NEXT   => null()          ! Terminated link list
end TYPE girder_info

TYPE girder_list
   CHARACTER(120), POINTER ::  NAME  => null()! IDENTIFICATION
   INTEGER,  POINTER :: N   => null()   ! TOTAL ELEMENT IN THE CHAIN

   INTEGER, POINTER :: discarded  => null()  ! duplicate magnets
   INTEGER, POINTER :: LASTPOS  => null()  ! POSITION OF LAST VISITED
   TYPE (girder), POINTER :: LAST  => null()! LAST VISITED
   !
   TYPE (girder), POINTER :: END  => null()
   TYPE (girder), POINTER :: START  => null()
   type (fibre), pointer :: lastfibre  => null()

END TYPE girder_list



TYPE MUL_BLOCK
   ! stuff for setting multipole
   real(dp) AN(NMAX),BN(NMAX)
   INTEGER NMUL,NATURAL,ADD
END TYPE MUL_BLOCK

type work
   real(dp) beta0,energy,kinetic,p0c,brho,gamma0I,gambet
   real(dp) mass
   LOGICAL(lp) rescale
   integer power
end type work


TYPE INTERNAL_STATE
   INTEGER TOTALPATH   ! total time or path length is used
   LOGICAL(LP) TIME  ! Time is used instead of path length
   LOGICAL(LP) RADIATION ! Radiation is turned on
   LOGICAL(LP) NOCAVITY ! Cavity is turned into a drift
   LOGICAL(LP) FRINGE ! Fringe fields are turned on (mainly for quadrupoles)
   LOGICAL(LP) STOCHASTIC ! Random Stochastic kicks to x(5) 
   LOGICAL(LP) ENVELOPE ! Stochastic envelope terms tracked in probe_8 
   LOGICAL(LP) PARA_IN  ! If true, parameters in the map are included
   LOGICAL(LP) ONLY_4D  ! REAL_8 Taylor in (x,p_x,y,p_y)
   LOGICAL(LP) DELTA  ! REAL_8 Taylor in (x,p_x,y,p_y,delta)
   LOGICAL(LP) SPIN  ! Spin is tracked
   LOGICAL(LP) MODULATION !  One modulated family tracked by probe
   LOGICAL(LP) ONLY_2D  !  REAL_8 Taylor in (x,p_x)
   LOGICAL(LP) FULL_WAY  ! 
END TYPE INTERNAL_STATE



TYPE POL_BLOCK
   CHARACTER(nlp) NAME
   integer n_name
   CHARACTER(vp) VORNAME
   ! STUFF FOR SETTING MAGNET USING GLOBAL ARRAY TPSAFIT
   real(dp),DIMENSION(:), POINTER :: TPSAFIT => null()
   LOGICAL(lp), POINTER ::  SET_TPSAFIT  => null()
   LOGICAL(lp), POINTER ::  SET_ELEMENT => null()
   ! STUFF FOR PARAMETER DEPENDENCE
   INTEGER NPARA
   INTEGER IAN(NMAX),IBN(NMAX)
   real(dp) SAN(NMAX),SBN(NMAX)
   INTEGER IVOLT, IFREQ,IPHAS
   INTEGER IB_SOL
   real(dp) SVOLT, SFREQ,SPHAS
   real(dp) SB_SOL
   integer g,np,nb   !  group index  number of blocks
   ! User defined Functions
   TYPE(POL_BLOCK_sagan) sagan
END TYPE POL_BLOCK

TYPE POL_BLOCK_INICOND
   INTEGER       :: BETA(3)
   INTEGER       :: ALFA(3)
   INTEGER       :: DISPERSION(4)
END TYPE POL_BLOCK_INICOND


TYPE MADX_APERTURE
   integer,pointer :: pos  => null()
   INTEGER,pointer ::  KIND  => null()  ! 1,2,3,4
   REAL(DP),pointer :: R(:)  => null()
   REAL(DP),pointer :: X => null(),Y => null(),DX => null(),DY => null()
END TYPE MADX_APERTURE

TYPE S_APERTURE
   type(MADX_APERTURE), pointer:: APERTURE => null()
END TYPE S_APERTURE


TYPE MAGNET_CHART
   type(magnet_frame), pointer:: f  => null()
   type(MADX_APERTURE), pointer:: APERTURE  => null()
   type(S_APERTURE), pointer:: A(:)  => null()
   real(dp),pointer :: charge  => null() ! propagator
   integer,pointer :: dir  => null()   ! propagator
   real(dp), POINTER :: BETA0 => null(),GAMMA0I => null(),GAMBET => null(),MASS => null(),ag => null()
   real(dp), POINTER :: P0C => null()
   real(dp), POINTER :: LD => null(),B0 => null(),LC => null()         !
   real(dp), POINTER :: TILTD  => null()     ! INTERNAL FRAME
   real(dp),  DIMENSION(:), POINTER :: EDGE  => null()        ! INTERNAL FRAME

   !
   LOGICAL(lp), POINTER :: EXACT  => null()
   LOGICAL(lp), POINTER :: KILL_ENT_FRINGE => null(),KILL_EXI_FRINGE => null(), bend_fringe => null()
   LOGICAL(lp), POINTER :: KILL_ENT_SPIN => null(),KILL_EXI_SPIN => null() 
   integer, pointer :: permFRINGE => null(),highest_fringe => null()                 !
   !
   INTEGER, POINTER :: METHOD => null(),NST => null()                   ! METHOD OF INTEGRATION 2,4,OR 6 YOSHIDA
   INTEGER, POINTER :: NMUL => null()                         ! NUMBER OF MULTIPOLE

END TYPE MAGNET_CHART

TYPE tilting
   real(dp) tilt(0:nmax)
   LOGICAL(lp) natural                 ! for mad-like
END TYPE tilting

  TYPE time_energy
     real(dp),pointer :: time => null()
     real(dp),pointer :: energy => null()
     real(dp),pointer :: an(:) => null(),bn(:) => null()
     real(dp),pointer :: b_t => null()
  END TYPE time_energy
  
  TYPE ramping
     integer,pointer :: n => null()
     real(dp), pointer :: r => null(),t_max => null()    !, unit_time
     type(time_energy),pointer :: table(:) => null()
     character(255), pointer :: file => null()
  END TYPE ramping



!   Element
TYPE ELEMENT
   INTEGER, POINTER :: KIND => null()
   logical(lp), POINTER :: recut,even => null()
   LOGICAL(LP), POINTER :: PLOT => null()
   LOGICAL(LP), POINTER :: electric => null()
   character(120), pointer :: filef
   character(120), pointer ::  fileb
   ! common stuff to all element
   type(MAGNET_CHART), pointer :: P => null()
   CHARACTER(nlp), POINTER ::  NAME  => null()   ! Identification
   CHARACTER(vp), POINTER ::  VORNAME  => null()   ! Identification
   !
   !
   ! Length is common although certain things like Markers should not have a length
   ! Well let us say it is zero
   real(dp), POINTER ::  L  => null()                              ! Length of integration often same as LD
   !
   real(dp),   DIMENSION(:), POINTER:: AN => null(),BN => null()         !Multipole component
   real(dp), DIMENSION(:),   POINTER:: FINT => null(),HGAP => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(),H2 => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA => null(),VS => null()         ! use in quad fringe from sad
   !
   real(dp), POINTER :: VOLT => null(), FREQ => null(),PHAS => null(),DELTA_E => null(),LAG => null()       ! Cavity information
   real(dp), POINTER :: DC_ac => null(),A_ac => null(),theta_ac => null(),D_AC => null()      ! slow RF: amplitude and phase
   real(dp),   DIMENSION(:), POINTER:: D_AN => null(),D_BN => null(),D0_AN => null(),D0_BN => null()  !Multipole component
   logical(lp), POINTER :: THIN => null()
   integer(2), pointer :: slow_ac => null()
   real(dp), POINTER ::  B_SOL  => null()
   !  misalignements and rotation
   logical(lp), POINTER ::  MIS => null()
   !   real(dp),  DIMENSION(:), POINTER ::d => null(),r => null()                             !  Misalignements
   !storage  space
   !integer  twiss                                                            !
   ! TYPES OF MAGNETS
   TYPE(DRIFT1), POINTER :: D0  => null()               ! DRIFT
   TYPE(DKD2), POINTER :: K2    => null()            ! INTEGRATOR
   TYPE(KICKT3), POINTER :: K3  => null()              !  THIN KICK
   TYPE(CAV4), POINTER :: C4    => null()            ! CAVITY
   TYPE(SOL5), POINTER :: S5    => null()            ! solenoid
   TYPE(KTK), POINTER :: T6     => null()           ! INTEGRATOR   thick slow
   TYPE(TKTF), POINTER :: T7    => null()            ! INTEGRATOR   thick fast
   TYPE(NSMI), POINTER :: S8    => null()            ! NORMAL SMI
   TYPE(SSMI), POINTER :: S9    => null()            ! SKEW SMI
   TYPE(TEAPOT), POINTER :: TP10 => null()                ! sector teapot
   TYPE(MON), POINTER :: MON14   => null()            ! MONITOR OR INSTRUMENT
   TYPE(ESEPTUM), POINTER :: SEP15 => null()              ! MONITOR OR INSTRUMENT
   TYPE(STREX), POINTER :: K16   => null()             ! EXACT STRAIGHT INTEGRATOR
   TYPE(ENGE), POINTER :: ENGE17 => null()               ! SOLENOID SIXTRACK STYLE
   TYPE(RCOL), POINTER :: RCOL18 => null()             ! RCOLLIMATOR
   TYPE(ECOL), POINTER :: ECOL19 => null()             ! ECOLLIMATOR
   TYPE(CAV_TRAV), POINTER :: CAV21 => null()            ! CAVITY TRAVELLING WAVE
   TYPE(SAGAN), POINTER :: WI   => null()              ! USER DEFINED
   TYPE(PANCAKE), POINTER :: PA => null()          ! GENERAL B
   TYPE(ABELL), POINTER :: AB    => null()            ! TYPE ABELL 
   TYPE(HELICAL_DIPOLE), POINTER :: HE22  => null()         !  
   TYPE(SUPERDRIFT), POINTER :: SDR      => null()     ! GENERAL DRIFT
   TYPE(FIBRE), POINTER :: PARENT_FIBRE => null()
   TYPE(fibre_appearance), POINTER :: doko => null()
   type(element), pointer :: siamese => null()
   type(element), pointer :: girders => null()
   TYPE(AFFINE_FRAME), POINTER :: SIAMESE_FRAME => null()
   TYPE(AFFINE_FRAME), POINTER :: girder_FRAME => null()
   !   integer, POINTER ::girder_index => null()
   type(girder),pointer :: assembly => null()
   type(ramping), pointer :: ramp => null()
   type(tree_element), pointer :: forward(:)=> null(),backward(:)=> null()
   logical, pointer :: usef=> null(),useb=> null()
   integer, pointer :: skip_ptc_f=> null(),skip_ptc_b=> null()
   logical, pointer :: do1mapf=> null(),do1mapb=> null()
   real(dp), pointer :: ene   ! for tapering
END TYPE  ELEMENT


TYPE ELEMENTP
   INTEGER, POINTER :: KIND  => null()! WHAT IT IS
   logical(lp), POINTER :: KNOB  => null()! FALSE IF NO KNOB
   CHARACTER(nlp), POINTER ::  NAME  => null()   ! Identification
   CHARACTER(vp), POINTER ::  VORNAME  => null()   ! Identification
   LOGICAL(LP), POINTER :: electric => null()
   !
   !
   !
   TYPE(REAL_8), POINTER ::  L    => null() ! LENGTH OF INTEGRATION OFTEN SAME AS LD, CAN BE ZERO
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(),BN   => null()       !MULTIPOLE COMPONENT
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(),HGAP  => null()        !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(),H2 => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),    POINTER:: VA => null(),VS => null()         !valishev-like multipole
   !
   TYPE(REAL_8), POINTER :: VOLT => null(), FREQ => null(),PHAS => null() ! CAVITY INFORMATION
   real(dp), POINTER :: DELTA_E  => null()    ! CAVITY ENERGY GAIN
   TYPE(REAL_8), POINTER :: DC_ac => null(),A_ac => null(),theta_ac => null(),D_AC => null()      ! slow RF: amplitude and phase
   TYPE(REAL_8),   DIMENSION(:), POINTER:: D_AN => null(),D_BN => null(),D0_AN => null(),D0_BN => null()         !Multipole component
   !
   TYPE(REAL_8), POINTER :: B_SOL => null()
   logical(lp), POINTER :: THIN => null() 
   integer(2), pointer :: slow_ac => null()
   !  MISALIGNEMENTS AND ROTATION
   logical(lp), POINTER ::  MIS => null()
   !   real(dp),  DIMENSION(:), POINTER :: D,R

   TYPE(MAGNET_CHART), POINTER :: P => null()

   ! TYPES OF POLYMORPHIC MAGNETS
   TYPE(DRIFT1P), POINTER :: D0  => null()            ! DRIFT
   TYPE(DKD2P), POINTER :: K2    => null()            ! INTEGRATOR
   TYPE(KICKT3P), POINTER :: K3  => null()            ! THIN KICK
   TYPE(CAV4P), POINTER :: C4    => null()            ! DRIFT
   TYPE(SOL5P), POINTER :: S5    => null()            ! CAVITY
   TYPE(KTKP), POINTER :: T6     => null()            ! INTEGRATOR
   TYPE(TKTFP), POINTER :: T7    => null()            ! INTEGRATOR   THICK FAST
   TYPE(NSMIP), POINTER :: S8    => null()            ! NORMAL SMI
   TYPE(SSMIP), POINTER :: S9    => null()            ! SKEW SMI
   TYPE(TEAPOTP), POINTER :: TP10 => null()           ! SECTOR BEND WITH CYLINDRICAL GEOMETRY
   TYPE(MONP), POINTER :: MON14   => null()            ! MONITOR OR INSTRUMENT
   TYPE(ESEPTUMP), POINTER :: SEP15 => null()              ! MONITOR OR INSTRUMENT
   TYPE(STREXP), POINTER :: K16     => null()           ! EXACT STRAIGHT INTEGRATOR
   TYPE(ENGEP), POINTER :: ENGE17   => null()             ! SOLENOID SIXTRACK STYLE
   TYPE(RCOLP), POINTER :: RCOL18   => null()           ! RCOLLIMATOR
   TYPE(ECOLP), POINTER :: ECOL19   => null()           ! ECOLLIMATOR
   TYPE(CAV_TRAVP), POINTER :: CAV21 => null()            ! CAVITY TRAVELLING WAVE
   TYPE(HELICAL_DIPOLEP), POINTER :: HE22 => null()          ! GENERAL B
   TYPE(SUPERDRIFTP), POINTER :: SDR  => null()         ! GENERAL DRIFT
   TYPE(SAGANP), POINTER :: WI         => null()        ! USER DEFINED
   TYPE(PANCAKEP), POINTER :: PA       => null()    ! GENERAL B
      TYPE(ABELLP), POINTER :: AB    => null()            ! TYPE ABELL 
   TYPE(FIBRE), POINTER :: PARENT_FIBRE => null()
   !   TYPE(fibre_appearance), POINTER :: doko
   type(ramping), pointer :: ramp => null()
   type(tree_element), pointer :: forward(:)=> null(),backward(:)=> null()
   logical, pointer :: usef=> null(),useb=> null()
   integer, pointer :: skip_ptc_f=> null(),skip_ptc_b=> null()
   logical, pointer :: do1mapf=> null(),do1mapb=> null()
END TYPE  ELEMENTP

type fibre_appearance
   TYPE(FIBRE), POINTER :: PARENT_FIBRE => null()
   TYPE(fibre_appearance), POINTER :: next => null()
end type fibre_appearance

type info
   real(dp),pointer :: s  => null()
   real(dp),pointer ::  beta(:)  => null()
   real(dp),pointer ::  fix0(:)  => null()
   real(dp),pointer ::  fix(:)  => null()
   real(dp), pointer:: pos(:)  => null()
   real(dp), pointer:: m(:,:)  => null()
END type info


TYPE INTEGRATION_NODE
   INTEGER, POINTER :: pos_in_fibre => null(), CAS  => null()
   INTEGER, POINTER ::  pos => null(),lost  => null()
   real(dp), POINTER :: S(:)  => null()
   real(dp), POINTER :: ds_ac  => null()
   real(dp), POINTER :: ref(:)  => null()
   real(dp), pointer :: ent(:,:),a(:)  => null()
   real(dp), pointer :: exi(:,:),b(:)  => null()
   real(dp), POINTER :: delta_rad_in  => null()
   real(dp), POINTER :: delta_rad_out  => null()
   INTEGER, POINTER :: TEAPOT_LIKE  => null()
   TYPE (INTEGRATION_NODE), POINTER :: NEXT  => null()
   TYPE (INTEGRATION_NODE), POINTER :: PREVIOUS  => null()
   TYPE (NODE_LAYOUT), POINTER :: PARENT_NODE_LAYOUT  => null()
   TYPE(FIBRE), POINTER :: PARENT_FIBRE  => null()
   !   TYPE(EXTRA_WORK), POINTER :: WORK
   TYPE(BEAM_BEAM_NODE), POINTER :: BB  => null()
!   TYPE(tree_element), POINTER :: T
END TYPE INTEGRATION_NODE

TYPE FIBRE
   !  BELOW ARE THE DATA CARRIED BY THE NODE
   INTEGER,POINTER ::DIR => null()
   TYPE(PATCH),POINTER ::PATCH => null()
   TYPE(CHART),POINTER ::CHART => null()
   TYPE (ELEMENT), POINTER ::  MAG => null()
   TYPE (ELEMENTP),POINTER ::  MAGP => null()
   !  END OF DATA
   !  POINTER TO THE MAGNETS ON EACH SIDE OF THIS NODE
   TYPE (FIBRE),POINTER :: PREVIOUS => null()
   TYPE (FIBRE),POINTER :: NEXT => null()
   !  POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA
   TYPE (LAYOUT),POINTER :: PARENT_LAYOUT => null()
   type(info),pointer ::i => null()
   TYPE(INTEGRATION_NODE),POINTER :: T1 => null()  ! FIRST AND LAST INTEGRATION_NODE CHILDREN CORRESPOUNDING TO PATCHES
   TYPE(INTEGRATION_NODE),POINTER ::  T2 => null()  ! FIRST AND LAST INTEGRATION_NODE CHILDREN CORRESPOUNDING TO PATCHES
   TYPE(INTEGRATION_NODE),POINTER :: TM  => null()     ! MIDDLE INTEGRATION_NODE
   INTEGER,POINTER ::pos  => null()             ! POSITION IN LAYOUT
   ! NEW STUFF....
   real(dp), POINTER :: BETA0 => null(),GAMMA0I => null(),GAMBET => null(),MASS  => null()  !,P0C
   real(dp), POINTER :: CHARGE => null()
   real(dp), POINTER :: AG  => null()! spin g-2
   ! TO TIE LAYOUTS
   TYPE (FIBRE),POINTER :: P => null()
   TYPE (FIBRE),POINTER :: N => null()
   INTEGER,POINTER :: loc => null()

END TYPE FIBRE

TYPE LAYOUT
   CHARACTER(120), POINTER ::  NAME  => null()! IDENTIFICATION
   INTEGER, POINTER ::  INDEX  => null()! IDENTIFICATION, CHARGE SIGN
   REAL(DP),  POINTER :: HARMONIC_NUMBER  => null()!
   logical(lp),POINTER ::CLOSED  => null()
   INTEGER,  POINTER :: N   => null()   ! TOTAL ELEMENT IN THE CHAIN
   INTEGER,POINTER ::NTHIN  => null() ! NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)
   REAL(DP),  POINTER :: THIN  => null()   ! PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS
   !POINTERS OF LINK LAYOUT
   INTEGER, POINTER :: LASTPOS  => null()  ! POSITION OF LAST VISITED
   TYPE (FIBRE), POINTER :: LAST  => null() ! LAST VISITED
   !
   TYPE (FIBRE), POINTER :: END  => null()
   TYPE (FIBRE), POINTER :: START  => null()
   TYPE (FIBRE), POINTER :: START_GROUND  => null() ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
   TYPE (FIBRE), POINTER :: END_GROUND  => null()! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
   TYPE (LAYOUT), POINTER :: NEXT  => null()
   TYPE (LAYOUT), POINTER :: PREVIOUS  => null()
   type(NODE_LAYOUT), pointer :: T    => null()!  ASSOCIATED  CHILD THIN LENS LAYOUT
   TYPE (MAD_UNIVERSE), POINTER :: parent_universe  => null()
   TYPE(layout_array), POINTER :: DNA(:)  => null()
END TYPE LAYOUT

type layout_array
   type(layout), pointer :: L  => null()
   integer   :: counter
END type layout_array

type girder_siamese
   type(element), pointer :: mag  => null()
END type girder_siamese

TYPE MAD_UNIVERSE             ! THE MOTHER OF ALL STRUCTURES
   INTEGER,  POINTER :: N   => null()   ! NUMBERS OF LAYOUT
   INTEGER,  POINTER :: SHARED  => null()    ! NUMBERS OF REPOSITORY LAYOUTS
   TYPE (LAYOUT), POINTER :: END  => null()
   TYPE (LAYOUT), POINTER :: START  => null()
   !POINTERS OF UNIVERSE IN CONJUNCTION WITH THE POINTERS N AND P OF THE FIBRES
   INTEGER, POINTER :: NF  => null()  ! number of fibres in tied layouts
   INTEGER, POINTER :: LASTPOS  => null()  ! POSITION OF LAST VISITED
   TYPE (FIBRE), POINTER :: LAST  => null() ! LAST VISITED
END TYPE MAD_UNIVERSE



TYPE BEAM_LOCATION
   TYPE (INTEGRATION_NODE), POINTER :: NODE  => null()
END TYPE BEAM_LOCATION

TYPE NODE_LAYOUT
   CHARACTER(120), POINTER ::  NAME  => null()! IDENTIFICATION
   INTEGER, POINTER ::  INDEX  => null()! IDENTIFICATION
   logical(lp),POINTER ::CLOSED  => null()
   INTEGER,  POINTER :: N   => null()   ! TOTAL ELEMENT IN THE CHAIN
   !POINTERS OF LINK LAYOUT
   INTEGER, POINTER :: LASTPOS  => null()  ! POSITION OF LAST VISITED
   TYPE (INTEGRATION_NODE), POINTER :: LAST  => null()! LAST VISITED
   !
   TYPE (INTEGRATION_NODE), POINTER :: END  => null()
   TYPE (INTEGRATION_NODE), POINTER :: START  => null()
   TYPE (INTEGRATION_NODE), POINTER :: START_GROUND  => null()! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
   TYPE (INTEGRATION_NODE), POINTER :: END_GROUND  => null()! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
   TYPE (LAYOUT), POINTER :: PARENT_LAYOUT  => null()
   TYPE(ORBIT_LATTICE), POINTER :: ORBIT_LATTICE  => null()
END TYPE NODE_LAYOUT

TYPE ORBIT_NODE
   TYPE (INTEGRATION_NODE), POINTER :: NODE  => null()
   REAL(DP), POINTER :: LATTICE(:)  => null()
   integer, POINTER :: DPOS  => null()
   integer, POINTER :: ENTERING_TASK  => null()
   integer, POINTER :: PTC_TASK  => null()
   logical, POINTER :: cavity  => null()
END TYPE ORBIT_NODE

TYPE ORBIT_LATTICE
   TYPE(ORBIT_NODE), pointer :: ORBIT_NODES(:)  => null()
   INTEGER, pointer  ::  ORBIT_N_NODE   => null()! GET_N_NODE
   LOGICAL(lp), pointer  ::  ORBIT_USE_ORBIT_UNITS => null(),accel => null()
   INTEGER, pointer  ::  ORBIT_WARNING  => null()
   REAL(DP), pointer  ::  ORBIT_LMAX   => null()  ! GET_LMAX
   REAL(DP), pointer  ::  ORBIT_MAX_PATCH_TZ  => null()
   REAL(DP), pointer  ::  ORBIT_mass_in_amu  => null() ! GET_MASS_AMU
   REAL(DP), pointer  ::  ORBIT_gammat   => null()  ! GET_GAMMAT
   REAL(DP), pointer  ::  ORBIT_harmonic  => null() ! GET_HARMONIC
   REAL(DP), pointer  ::  ORBIT_L  => null()   ! GET_CIRCUMFERENCE
   real(dp),pointer ::  ORBIT_CHARGE  => null() ! GET_CHARGE
   REAL(DP), pointer  ::  ORBIT_OMEGA  => null() ! GET_omega
   REAL(DP), pointer  ::  ORBIT_P0C  => null()  ! GET_P0C
   REAL(DP), pointer  ::  ORBIT_BETA0  => null()! GET_BETA0
   REAL(DP), pointer  ::  orbit_kinetic  => null() ! GET_kinetic
   REAL(DP), pointer  ::  orbit_brho  => null()  ! GET_brho
   REAL(DP), pointer  ::  orbit_energy  => null()! GET_total_energy
   REAL(DP), pointer  ::  orbit_gamma  => null() ! GET_gamma
   !   REAL(DP), pointer  ::  orbit_dppfac ! GET_dppfac
   REAL(DP), pointer  ::  orbit_deltae  => null()! GET_deltae
   REAL(DP), pointer  ::  ORBIT_OMEGA_after => null()
   !   REAL(DP), pointer  ::  freqb,freqa,voltb,volta,phasa,phasb,xs6,dxs6
   TYPE(INTERNAL_STATE),pointer :: STATE => null()
   TYPE(INTEGRATION_NODE), pointer :: tp  => null()
   TYPE(layout), pointer :: parent_layout  => null()
   real(dp), pointer :: dt(:)  => null()
END TYPE ORBIT_LATTICE


!  BEAM BEAM

TYPE BEAM_BEAM_NODE
   REAL(DP), POINTER :: S  => null()
   REAL(DP), POINTER :: SX => null(),SY => null(),FK => null()
   REAL(DP), POINTER :: XM => null(),YM => null()  ! trivial transverse displacement: not needed
   REAL(DP), POINTER :: BBK(:)  => null() ! kick on closed orbit
   !    REAL(DP), pointer :: MID(:,:),O(:)
   REAL(DP), POINTER :: A(:)  => null()   ! patch angles
   REAL(DP), POINTER :: D(:)  => null()   ! patch translation
   INTEGER, POINTER :: A_X1 => null(), A_X2 => null() ! patch inverse propagator (not needed)
   LOGICAL(LP), POINTER :: PATCH  => null()! true= patch
END TYPE BEAM_BEAM_NODE

!  Additional work at node
TYPE EXTRA_WORK
   INTEGER, POINTER :: KIND  => null()
   TYPE(INTEGRATION_NODE), POINTER :: NODE  => null()
   TYPE(BEAM_BEAM_NODE), POINTER :: BB  => null()
   TYPE(MADX_APERTURE), POINTER :: A  => null()
   REAL(DP), POINTER :: ORB(:)  => null()
END TYPE EXTRA_WORK

!!!!  graphical arrays used in Su_hiromi only now
  real(sp), allocatable :: z_hiromi(:,:) ,z_hiromi0(:,:)
  integer n_hiromi,above_kicker,n_hiromi_lost
  logical(lp) :: fill_hiromi =my_false
  real(dp) :: r_hiromi=10.0,h_hiromi=10.0,kicker_hiromi=0.4d0,angle_hiromi(3)=0.d0
  real(dp) :: ent_hiromi(3,3),x0_hiromi(6),t0_hiromi,dt_hiromi
  type(internal_state) state_hiromi

!!! used by y_lin_beam only

TYPE E_BEAM
   integer N 
   type(probe), POINTER :: Z(:)  => null()
END TYPE E_BEAM

type(e_beam), pointer :: my_e_beam =>null()
