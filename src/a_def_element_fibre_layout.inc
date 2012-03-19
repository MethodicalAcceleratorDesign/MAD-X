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
   INTEGER, pointer :: i
   INTEGER, pointer :: pos
   INTEGER, pointer :: mul    ! number of magnets in the girder
   INTEGER, pointer :: added
   INTEGER, pointer :: discarded
   real(dp), pointer:: a(:),ent(:,:)
   TYPE (girder_info),POINTER :: info
   TYPE (girder),POINTER :: PREVIOUS   ! Terminated link list
   TYPE (girder),POINTER :: NEXT
END TYPE girder

TYPE girder_info
   CHARACTER(nlp), POINTER ::  NAME    ! Identification
   real(dp), POINTER ::  d(:)
   type (element), pointer :: mag
   real(dp), pointer:: a(:),ent(:,:)
   type (girder), pointer :: parent_girder
   TYPE (girder_info),POINTER :: NEXT            ! Terminated link list
end TYPE girder_info

TYPE girder_list
   CHARACTER(120), POINTER ::  NAME ! IDENTIFICATION
   INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN

   INTEGER, POINTER :: discarded   ! duplicate magnets
   INTEGER, POINTER :: LASTPOS   ! POSITION OF LAST VISITED
   TYPE (girder), POINTER :: LAST ! LAST VISITED
   !
   TYPE (girder), POINTER :: END
   TYPE (girder), POINTER :: START
   type (fibre), pointer :: lastfibre

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
   INTEGER TOTALPATH
   !    LOGICAL(LP) TIME,RADIATION,NOCAVITY,FRINGE,EXACTMIS,STOCHASTIC,ENVELOPE
   LOGICAL(LP) TIME,RADIATION,NOCAVITY,FRINGE,STOCHASTIC,ENVELOPE
   LOGICAL(LP) PARA_IN,ONLY_4D,DELTA
   LOGICAL(LP) SPIN,MODULATION
END TYPE INTERNAL_STATE



TYPE POL_BLOCK
   CHARACTER(nlp) NAME
   integer n_name
   CHARACTER(vp) VORNAME
   ! STUFF FOR SETTING MAGNET USING GLOBAL ARRAY TPSAFIT
   real(dp),DIMENSION(:), POINTER :: TPSAFIT
   LOGICAL(lp), POINTER ::  SET_TPSAFIT
   LOGICAL(lp), POINTER ::  SET_ELEMENT
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
   INTEGER,pointer ::  KIND   ! 1,2,3,4
   REAL(DP),pointer :: R(:)
   REAL(DP),pointer :: X,Y,DX,DY
END TYPE MADX_APERTURE

TYPE S_APERTURE
   type(MADX_APERTURE), pointer:: APERTURE
END TYPE S_APERTURE


TYPE MAGNET_CHART
   type(magnet_frame), pointer:: f
   type(MADX_APERTURE), pointer:: APERTURE
   type(S_APERTURE), pointer:: A(:)
   integer,pointer :: charge  ! propagator
   integer,pointer :: dir    ! propagator
   real(dp), POINTER :: BETA0,GAMMA0I,GAMBET,MASS,ag
   real(dp), POINTER :: P0C
   real(dp), POINTER :: LD,B0,LC         !
   real(dp), POINTER :: TILTD      ! INTERNAL FRAME
   real(dp),  DIMENSION(:), POINTER :: EDGE         ! INTERNAL FRAME

   !
   LOGICAL(lp), POINTER :: EXACT
   LOGICAL(lp), POINTER :: KILL_ENT_FRINGE,KILL_EXI_FRINGE, bend_fringe,permFRINGE                  !
   !
   INTEGER, POINTER :: METHOD,NST                   ! METHOD OF INTEGRATION 2,4,OR 6 YOSHIDA
   INTEGER, POINTER :: NMUL                         ! NUMBER OF MULTIPOLE

END TYPE MAGNET_CHART

TYPE tilting
   real(dp) tilt(0:nmax)
   LOGICAL(lp) natural                 ! for mad-like
END TYPE tilting

  TYPE time_energy
     real(dp),pointer :: time
     real(dp),pointer :: energy
     real(dp),pointer :: an(:),bn(:)
     real(dp),pointer :: b_t
  END TYPE time_energy
  
  TYPE ramping
     integer,pointer :: n
     real(dp), pointer :: r, unit_time,t_max    
     type(time_energy),pointer :: table(:)
     character(255), pointer :: file
  END TYPE ramping



!   Element
TYPE ELEMENT
   INTEGER, POINTER :: KIND
   logical(lp), POINTER :: recut,even
   LOGICAL(LP), POINTER :: PLOT
   LOGICAL(LP), POINTER :: electric

   ! common stuff to all element
   type(MAGNET_CHART), pointer :: P
   CHARACTER(nlp), POINTER ::  NAME    ! Identification
   CHARACTER(vp), POINTER ::  VORNAME    ! Identification
   !
   !
   ! Length is common although certain things like Markers should not have a length
   ! Well let us say it is zero
   real(dp), POINTER ::  L                               ! Length of integration often same as LD
   !
   real(dp),   DIMENSION(:), POINTER:: AN,BN         !Multipole component
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA,VS         !valishev-like multipole
   !
   real(dp), POINTER :: VOLT, FREQ,PHAS,DELTA_E,LAG       ! Cavity information
   real(dp), POINTER :: DC_ac,A_ac,theta_ac,D_AC      ! slow RF: amplitude and phase
   real(dp),   DIMENSION(:), POINTER:: D_AN,D_BN,D0_AN,D0_BN         !Multipole component
   logical(lp), POINTER :: THIN,slow_ac
   real(dp), POINTER ::  B_SOL                                          ! Solenoidal field
   !  misalignements and rotation
   logical(lp), POINTER ::  MIS
   !   real(dp),  DIMENSION(:), POINTER ::d,r                             !  Misalignements
   !storage  space
   !integer  twiss                                                            !
   ! TYPES OF MAGNETS
   TYPE(DRIFT1), POINTER :: D0               ! DRIFT
   TYPE(DKD2), POINTER :: K2               ! INTEGRATOR
   TYPE(KICKT3), POINTER :: K3               !  THIN KICK
   TYPE(CAV4), POINTER :: C4               ! CAVITY
   TYPE(SOL5), POINTER :: S5               ! solenoid
   TYPE(KTK), POINTER :: T6               ! INTEGRATOR   thick slow
   TYPE(TKTF), POINTER :: T7               ! INTEGRATOR   thick fast
   TYPE(NSMI), POINTER :: S8               ! NORMAL SMI
   TYPE(SSMI), POINTER :: S9               ! SKEW SMI
   TYPE(TEAPOT), POINTER :: TP10                ! sector teapot
   TYPE(MON), POINTER :: MON14              ! MONITOR OR INSTRUMENT
   TYPE(ESEPTUM), POINTER :: SEP15              ! MONITOR OR INSTRUMENT
   TYPE(STREX), POINTER :: K16               ! EXACT STRAIGHT INTEGRATOR
   TYPE(ENGE), POINTER :: ENGE17               ! SOLENOID SIXTRACK STYLE
   TYPE(RCOL), POINTER :: RCOL18             ! RCOLLIMATOR
   TYPE(ECOL), POINTER :: ECOL19             ! ECOLLIMATOR
   TYPE(CAV_TRAV), POINTER :: CAV21            ! CAVITY TRAVELLING WAVE
   TYPE(SAGAN), POINTER :: WI                ! USER DEFINED
   TYPE(PANCAKE), POINTER :: PA          ! GENERAL B
   TYPE(HELICAL_DIPOLE), POINTER :: HE22          ! GENERAL B
   TYPE(FIBRE), POINTER :: PARENT_FIBRE
   TYPE(fibre_appearance), POINTER :: doko
   type(element), pointer :: siamese
   type(element), pointer :: girders
   TYPE(AFFINE_FRAME), POINTER :: SIAMESE_FRAME
   TYPE(AFFINE_FRAME), POINTER :: girder_FRAME
   !   integer, POINTER ::girder_index
   type(girder),pointer :: assembly
   type(ramping), pointer :: ramp

END TYPE  ELEMENT


TYPE ELEMENTP
   INTEGER, POINTER :: KIND ! WHAT IT IS
   logical(lp), POINTER :: KNOB ! FALSE IF NO KNOB
   CHARACTER(nlp), POINTER ::  NAME    ! Identification
   CHARACTER(vp), POINTER ::  VORNAME    ! Identification
   LOGICAL(LP), POINTER :: electric
   !
   !
   !
   TYPE(REAL_8), POINTER ::  L    ! LENGTH OF INTEGRATION OFTEN SAME AS LD, CAN BE ZERO
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !MULTIPOLE COMPONENT
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),    POINTER:: VA,VS         !valishev-like multipole
   !
   TYPE(REAL_8), POINTER :: VOLT, FREQ,PHAS ! CAVITY INFORMATION
   real(dp), POINTER :: DELTA_E     ! CAVITY ENERGY GAIN
   TYPE(REAL_8), POINTER :: DC_ac,A_ac,theta_ac,D_AC      ! slow RF: amplitude and phase
   TYPE(REAL_8),   DIMENSION(:), POINTER:: D_AN,D_BN,D0_AN,D0_BN         !Multipole component
   !
   TYPE(REAL_8), POINTER :: B_SOL
   logical(lp), POINTER :: THIN,slow_ac

   !  MISALIGNEMENTS AND ROTATION
   logical(lp), POINTER ::  MIS
   !   real(dp),  DIMENSION(:), POINTER :: D,R

   TYPE(MAGNET_CHART), POINTER :: P

   ! TYPES OF POLYMORPHIC MAGNETS
   TYPE(DRIFT1P), POINTER :: D0             ! DRIFT
   TYPE(DKD2P), POINTER :: K2               ! INTEGRATOR
   TYPE(KICKT3P), POINTER :: K3             ! THIN KICK
   TYPE(CAV4P), POINTER :: C4               ! DRIFT
   TYPE(SOL5P), POINTER :: S5               ! CAVITY
   TYPE(KTKP), POINTER :: T6                ! INTEGRATOR
   TYPE(TKTFP), POINTER :: T7               ! INTEGRATOR   THICK FAST
   TYPE(NSMIP), POINTER :: S8               ! NORMAL SMI
   TYPE(SSMIP), POINTER :: S9               ! SKEW SMI
   TYPE(TEAPOTP), POINTER :: TP10           ! SECTOR BEND WITH CYLINDRICAL GEOMETRY
   TYPE(MONP), POINTER :: MON14              ! MONITOR OR INSTRUMENT
   TYPE(ESEPTUMP), POINTER :: SEP15              ! MONITOR OR INSTRUMENT
   TYPE(STREXP), POINTER :: K16               ! EXACT STRAIGHT INTEGRATOR
   TYPE(ENGEP), POINTER :: ENGE17               ! SOLENOID SIXTRACK STYLE
   TYPE(RCOLP), POINTER :: RCOL18             ! RCOLLIMATOR
   TYPE(ECOLP), POINTER :: ECOL19             ! ECOLLIMATOR
   TYPE(CAV_TRAVP), POINTER :: CAV21            ! CAVITY TRAVELLING WAVE
   TYPE(HELICAL_DIPOLEP), POINTER :: HE22          ! GENERAL B
   TYPE(SAGANP), POINTER :: WI                ! USER DEFINED
   TYPE(PANCAKEP), POINTER :: PA          ! GENERAL B
   TYPE(FIBRE), POINTER :: PARENT_FIBRE
   !   TYPE(fibre_appearance), POINTER :: doko
   type(ramping), pointer :: ramp
END TYPE  ELEMENTP

type fibre_appearance
   TYPE(FIBRE), POINTER :: PARENT_FIBRE
   TYPE(fibre_appearance), POINTER :: next
end type fibre_appearance

type info
   real(sp),pointer :: s
   real(sp),pointer ::  beta(:)
   real(sp),pointer ::  fix0(:)
   real(sp),pointer ::  fix(:)
   real(sp), pointer:: pos(:)
END type info

TYPE FIBRE
   !  BELOW ARE THE DATA CARRIED BY THE NODE
   INTEGER,POINTER ::DIR
   TYPE(PATCH),POINTER ::PATCH
   TYPE(CHART),POINTER ::CHART
   TYPE (ELEMENT), POINTER ::  MAG
   TYPE (ELEMENTP),POINTER ::  MAGP
   !  END OF DATA
   !  POINTER TO THE MAGNETS ON EACH SIDE OF THIS NODE
   TYPE (FIBRE),POINTER :: PREVIOUS
   TYPE (FIBRE),POINTER :: NEXT
   !  POINTING TO PARENT LAYOUT AND PARENT FIBRE DATA
   TYPE (LAYOUT),POINTER :: PARENT_LAYOUT
   type(info),pointer ::i
   TYPE(INTEGRATION_NODE),POINTER :: T1,T2   ! FIRST AND LAST INTEGRATION_NODE CHILDREN CORRESPOUNDING TO PATCHES
   TYPE(INTEGRATION_NODE),POINTER :: TM      ! MIDDLE INTEGRATION_NODE
   INTEGER,POINTER ::pos              ! POSITION IN LAYOUT
   ! NEW STUFF....
   real(dp), POINTER :: BETA0,GAMMA0I,GAMBET,MASS   !,P0C
   INTEGER, POINTER :: CHARGE
   real(dp), POINTER :: AG  ! spin g-2
   ! TO TIE LAYOUTS
   TYPE (FIBRE),POINTER :: P
   TYPE (FIBRE),POINTER :: N
   INTEGER,POINTER :: loc

END TYPE FIBRE

TYPE LAYOUT
   CHARACTER(120), POINTER ::  NAME ! IDENTIFICATION
   INTEGER, POINTER ::  INDEX,HARMONIC_NUMBER ! IDENTIFICATION, CHARGE SIGN
   logical(lp),POINTER ::CLOSED
   INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN
   INTEGER,POINTER ::NTHIN  ! NUMBER IF THIN LENSES IN COLLECTION  (FOR SPEED ESTIMATES)
   REAL(DP),  POINTER :: THIN    ! PARAMETER USED FOR AUTOMATIC CUTTING INTO THIN LENS
   !POINTERS OF LINK LAYOUT
   INTEGER, POINTER :: LASTPOS   ! POSITION OF LAST VISITED
   TYPE (FIBRE), POINTER :: LAST ! LAST VISITED
   !
   TYPE (FIBRE), POINTER :: END
   TYPE (FIBRE), POINTER :: START
   TYPE (FIBRE), POINTER :: START_GROUND ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
   TYPE (FIBRE), POINTER :: END_GROUND ! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
   TYPE (LAYOUT), POINTER :: NEXT
   TYPE (LAYOUT), POINTER :: PREVIOUS
   type(NODE_LAYOUT), pointer :: T   !  ASSOCIATED  CHILD THIN LENS LAYOUT
   TYPE (MAD_UNIVERSE), POINTER :: parent_universe
   TYPE(layout_array), POINTER :: DNA(:)
   !   TYPE(siamese_array), POINTER :: girder(:)
   !   TYPE(junction_array), POINTER :: con2(:)
END TYPE LAYOUT

type layout_array
   type(layout), pointer :: L
   integer   :: counter
END type layout_array

type girder_siamese
   type(element), pointer :: mag
END type girder_siamese

TYPE MAD_UNIVERSE             ! THE MOTHER OF ALL STRUCTURES
   INTEGER,  POINTER :: N     ! NUMBERS OF LAYOUT
   INTEGER,  POINTER :: SHARED     ! NUMBERS OF REPOSITORY LAYOUTS
   TYPE (LAYOUT), POINTER :: END
   TYPE (LAYOUT), POINTER :: START
   !POINTERS OF UNIVERSE IN CONJUNCTION WITH THE POINTERS N AND P OF THE FIBRES
   INTEGER, POINTER :: NF   ! number of fibres in tied layouts
   INTEGER, POINTER :: LASTPOS   ! POSITION OF LAST VISITED
   TYPE (FIBRE), POINTER :: LAST ! LAST VISITED
END TYPE MAD_UNIVERSE

TYPE BEAM
   !   TYPE(REAL_8), POINTER :: Y(:)
   REAL(DP), POINTER :: X(:,:)
   !   REAL(DP), POINTER :: SIGMA(:),DX(:),BBPAR,ORBIT(:)
   LOGICAL(LP), POINTER :: U(:)
   TYPE(BEAM_LOCATION), POINTER::POS(:)
   INTEGER, POINTER :: N,LOST
   REAL(DP), POINTER :: A(:),D(:)
   !   INTEGER, POINTER :: CHARGE
   !   LOGICAL(LP),POINTER :: TIME_INSTEAD_OF_S
   !   LOGICAL(LP),POINTER :: BEAM_BEAM,BBORBIT
END TYPE BEAM

TYPE INTEGRATION_NODE
   INTEGER, POINTER :: pos_in_fibre, CAS
   INTEGER, POINTER ::  pos,lost
   real(dp), POINTER :: S(:)
   real(dp), POINTER :: ds_ac
   real(dp), POINTER :: ref(:)
   real(dp), pointer :: ent(:,:),a(:)
   real(dp), pointer :: exi(:,:),b(:)
   real(dp), POINTER :: delta_rad_in
   real(dp), POINTER :: delta_rad_out
   INTEGER, POINTER :: TEAPOT_LIKE
   TYPE (INTEGRATION_NODE), POINTER :: NEXT
   TYPE (INTEGRATION_NODE), POINTER :: PREVIOUS
   TYPE (NODE_LAYOUT), POINTER :: PARENT_NODE_LAYOUT
   TYPE(FIBRE), POINTER :: PARENT_FIBRE
   !   TYPE(EXTRA_WORK), POINTER :: WORK
   TYPE(BEAM_BEAM_NODE), POINTER :: BB
   TYPE(tree_element), POINTER :: T
END TYPE INTEGRATION_NODE

TYPE BEAM_LOCATION
   TYPE (INTEGRATION_NODE), POINTER :: NODE
END TYPE BEAM_LOCATION

TYPE NODE_LAYOUT
   CHARACTER(120), POINTER ::  NAME ! IDENTIFICATION
   INTEGER, POINTER ::  INDEX ! IDENTIFICATION
   logical(lp),POINTER ::CLOSED
   INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN
   !POINTERS OF LINK LAYOUT
   INTEGER, POINTER :: LASTPOS   ! POSITION OF LAST VISITED
   TYPE (INTEGRATION_NODE), POINTER :: LAST ! LAST VISITED
   !
   TYPE (INTEGRATION_NODE), POINTER :: END
   TYPE (INTEGRATION_NODE), POINTER :: START
   TYPE (INTEGRATION_NODE), POINTER :: START_GROUND ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
   TYPE (INTEGRATION_NODE), POINTER :: END_GROUND ! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
   TYPE (LAYOUT), POINTER :: PARENT_LAYOUT
   TYPE(ORBIT_LATTICE), POINTER :: ORBIT_LATTICE
END TYPE NODE_LAYOUT

TYPE ORBIT_NODE
   TYPE (INTEGRATION_NODE), POINTER :: NODE
   REAL(DP), POINTER :: LATTICE(:)
   integer, POINTER :: DPOS
   integer, POINTER :: ENTERING_TASK
   integer, POINTER :: PTC_TASK
   logical, POINTER :: cavity
END TYPE ORBIT_NODE

TYPE ORBIT_LATTICE
   TYPE(ORBIT_NODE), pointer :: ORBIT_NODES(:)
   INTEGER, pointer  ::  ORBIT_N_NODE  ! GET_N_NODE
   LOGICAL(lp), pointer  ::  ORBIT_USE_ORBIT_UNITS,accel
   INTEGER, pointer  ::  ORBIT_WARNING
   REAL(DP), pointer  ::  ORBIT_LMAX    ! GET_LMAX
   REAL(DP), pointer  ::  ORBIT_MAX_PATCH_TZ
   REAL(DP), pointer  ::  ORBIT_mass_in_amu  ! GET_MASS_AMU
   REAL(DP), pointer  ::  ORBIT_gammat    ! GET_GAMMAT
   REAL(DP), pointer  ::  ORBIT_harmonic  ! GET_HARMONIC
   REAL(DP), pointer  ::  ORBIT_L    ! GET_CIRCUMFERENCE
   INTEGER,pointer ::  ORBIT_CHARGE  ! GET_CHARGE
   REAL(DP), pointer  ::  ORBIT_OMEGA  ! GET_omega
   REAL(DP), pointer  ::  ORBIT_P0C   ! GET_P0C
   REAL(DP), pointer  ::  ORBIT_BETA0 ! GET_BETA0
   REAL(DP), pointer  ::  orbit_kinetic  ! GET_kinetic
   REAL(DP), pointer  ::  orbit_brho   ! GET_brho
   REAL(DP), pointer  ::  orbit_energy ! GET_total_energy
   REAL(DP), pointer  ::  orbit_gamma ! GET_gamma
   !   REAL(DP), pointer  ::  orbit_dppfac ! GET_dppfac
   REAL(DP), pointer  ::  orbit_deltae ! GET_deltae
   REAL(DP), pointer  ::  ORBIT_OMEGA_after
   !   REAL(DP), pointer  ::  freqb,freqa,voltb,volta,phasa,phasb,xs6,dxs6
   TYPE(INTERNAL_STATE),pointer :: STATE
   TYPE(INTEGRATION_NODE), pointer :: tp
   real(dp), pointer :: dt(:)
END TYPE ORBIT_LATTICE


!  BEAM BEAM

TYPE BEAM_BEAM_NODE
   REAL(DP), POINTER :: S
   REAL(DP), POINTER :: SX,SY,FK
   REAL(DP), POINTER :: XM,YM   ! trivial transverse displacement: not needed
   REAL(DP), POINTER :: BBK(:)  ! kick on closed orbit
   !    REAL(DP), pointer :: MID(:,:),O(:)
   REAL(DP), POINTER :: A(:)    ! patch angles
   REAL(DP), POINTER :: D(:)    ! patch translation
   INTEGER, POINTER :: A_X1,A_X2 ! patch inverse propagator (not needed)
   LOGICAL(LP), POINTER :: PATCH ! true= patch
END TYPE BEAM_BEAM_NODE

!  Additional work at node
TYPE EXTRA_WORK
   INTEGER, POINTER :: KIND
   TYPE(INTEGRATION_NODE), POINTER :: NODE
   TYPE(BEAM_BEAM_NODE), POINTER :: BB
   TYPE(MADX_APERTURE), POINTER :: A
   REAL(DP), POINTER :: ORB(:)
END TYPE EXTRA_WORK

TYPE fibre_array
   type(fibre), pointer :: p
   integer, pointer :: pos
END TYPE fibre_array
