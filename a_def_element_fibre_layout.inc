TYPE INTERNAL_STATE
   LOGICAL(lp) TOTALPATH,TIME,RADIATION,NOCAVITY,FRINGE,EXACTMIS
   LOGICAL(lp) PARA_IN,ONLY_4D,DELTA  ! DA related
END TYPE INTERNAL_STATE

type work
   real(dp) beta0,energy,kinetic,p0c,brho,gamma0I,gambet
   real(dp) mass
   LOGICAL(lp) rescale
   integer power
end type work


TYPE POL_BLOCK
   CHARACTER(nlp) NAME
   integer n_name
   CHARACTER(vp) VORNAME
   ! STUFF FOR SETTING MAGNET USING GLOBAL ARRAY TPSAFIT
   real(dp),DIMENSION(:), POINTER :: TPSAFIT
   LOGICAL(lp), POINTER ::  SET_TPSAFIT
   ! STUFF FOR PARAMETER DEPENDENCE
   INTEGER NPARA
   INTEGER IAN(NMAX),IBN(NMAX)
   real(dp) SAN(NMAX),SBN(NMAX)
   INTEGER IVOLT, IFREQ,IPHAS
   INTEGER IB_SOL
   real(dp) SVOLT, SFREQ,SPHAS
   real(dp) SB_SOL
   ! User defined Functions
   TYPE(POL_BLOCK_sagan) sagan
   TYPE(POL_BLOCK1) USER1
   TYPE(POL_BLOCK2) USER2
END TYPE POL_BLOCK

TYPE MADX_APERTURE
   INTEGER,pointer ::  KIND   ! 1,2,3,4
   REAL(DP),pointer :: R(:)
   REAL(DP),pointer :: X,Y
END TYPE MADX_APERTURE

TYPE MAGNET_CHART
   type(magnet_frame), pointer:: f
   type(MADX_APERTURE), pointer:: APERTURE
   integer,pointer :: charge  ! propagator
   integer,pointer :: dir    ! propagator
   real(dp), POINTER :: LD,B0,LC         !
   real(dp), POINTER :: TILTD      ! INTERNAL FRAME
   real(dp), POINTER :: BETA0,GAMMA0I,GAMBET,P0C
   !frs     real(dp),  DIMENSION(:), POINTER :: EDGE(:)      ! INTERNAL FRAME
   real(dp),  DIMENSION(:), POINTER :: EDGE         ! INTERNAL FRAME

   !
   INTEGER, POINTER :: TOTALPATH                    !
   LOGICAL(lp), POINTER :: EXACT,RADIATION,NOCAVITY     !       STATE
   LOGICAL(lp), POINTER :: FRINGE,KILL_ENT_FRINGE,KILL_EXI_FRINGE,TIME, bend_fringe                  !
   !
   INTEGER, POINTER :: METHOD,NST                   ! METHOD OF INTEGRATION 2,4,OR 6 YOSHIDA
   INTEGER, POINTER :: NMUL                         ! NUMBER OF MULTIPOLE

END TYPE MAGNET_CHART


TYPE tilting
   real(dp) tilt(0:nmax)
   LOGICAL(lp) natural                 ! for mad-like
END TYPE tilting

!   Element
TYPE ELEMENT
   INTEGER, POINTER :: KIND
   ! common stuff to all element
   type(MAGNET_CHART), pointer :: P
   CHARACTER(nlp), POINTER ::  NAME    ! Identification
   CHARACTER(vp), POINTER ::  VORNAME    ! Identification
   !
   logical(lp), POINTER ::  PERMFRINGE
   !
   ! Length is common although certain things like Markers should not have a length
   ! Well let us say it is zero
   real(dp), POINTER ::  L                               ! Length of integration often same as LD
   !
   real(dp),   DIMENSION(:), POINTER:: AN,BN         !Multipole component
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2         !FRINGE FUDGE FOR MAD
   !
   real(dp), POINTER :: VOLT, FREQ,PHAS,DELTA_E       ! Cavity information
   real(dp), POINTER ::  B_SOL                                          ! Solenoidal field
   logical(lp), POINTER :: THIN
   !  misalignements and rotation
   logical(lp), POINTER ::  MIS,EXACTMIS
   real(dp),  DIMENSION(:), POINTER ::d,r                             !  Misalignements
   !storage  space
   !integer  twiss                                                            !
   ! TYPES OF MAGNETS
   TYPE(FITTED_MAGNET), POINTER :: BEND               ! Machida's magnet
   TYPE(DRIFT1), POINTER :: D0               ! DRIFT
   TYPE(DKD2), POINTER :: K2               ! INTEGRATOR
   TYPE(KICKT3), POINTER :: K3               !  THIN KICK
   TYPE(CAV4), POINTER :: C4               ! CAVITY
   TYPE(SOL5), POINTER :: S5               ! CAVITY
   TYPE(KTK), POINTER :: T6               ! INTEGRATOR   thick slow
   TYPE(TKTF), POINTER :: T7               ! INTEGRATOR   thick fast
   TYPE(NSMI), POINTER :: S8               ! NORMAL SMI
   TYPE(SSMI), POINTER :: S9               ! SKEW SMI
   TYPE(TEAPOT), POINTER :: TP10                ! sector teapot
   TYPE(MON), POINTER :: MON14              ! MONITOR OR INSTRUMENT
   TYPE(ESEPTUM), POINTER :: SEP15              ! MONITOR OR INSTRUMENT
   TYPE(STREX), POINTER :: K16               ! EXACT STRAIGHT INTEGRATOR
   TYPE(SOLT), POINTER :: S17               ! SOLENOID SIXTRACK STYLE
   TYPE(RCOL), POINTER :: RCOL18             ! RCOLLIMATOR
   TYPE(ECOL), POINTER :: ECOL19             ! ECOLLIMATOR
   TYPE(CAV_TRAV), POINTER :: CAV21            ! CAVITY TRAVELLING WAVE
   TYPE(USER1), POINTER :: U1                ! USER DEFINED
   TYPE(USER2), POINTER :: U2                ! USER DEFINED
   TYPE(SAGAN), POINTER :: WI                ! USER DEFINED
   TYPE(MARTIN), POINTER :: M22               ! TAYLOR MAP
   TYPE(MULTIP), POINTER :: MU             ! MULTIP ZGOUBI
   TYPE(layout), POINTER :: g23
   TYPE(FIBRE), POINTER :: PARENT_FIBRE
END TYPE  ELEMENT


TYPE ELEMENTP
   INTEGER, POINTER :: KIND ! WHAT IT IS
   logical(lp), POINTER :: KNOB ! FALSE IF NO KNOB
   CHARACTER(nlp), POINTER ::  NAME    ! Identification
   CHARACTER(vp), POINTER ::  VORNAME    ! Identification
   logical(lp), POINTER ::  PERMFRINGE
   !
   !
   !
   TYPE(REAL_8), POINTER ::  L    ! LENGTH OF INTEGRATION OFTEN SAME AS LD, CAN BE ZERO
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !MULTIPOLE COMPONENT
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2         !FRINGE FUDGE FOR MAD
   !
   TYPE(REAL_8), POINTER :: VOLT, FREQ,PHAS ! CAVITY INFORMATION
   real(dp), POINTER :: DELTA_E     ! CAVITY ENERGY GAIN
   !
   TYPE(REAL_8), POINTER :: B_SOL
   logical(lp), POINTER :: THIN

   !  MISALIGNEMENTS AND ROTATION
   logical(lp), POINTER ::  MIS,EXACTMIS
   real(dp),  DIMENSION(:), POINTER :: D,R

   TYPE(MAGNET_CHART), POINTER :: P

   ! TYPES OF POLYMORPHIC MAGNETS
   TYPE(FITTED_MAGNETP), POINTER :: BEND    ! MACHIDA'S FITTED MAGNET
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
   TYPE(SOLTP), POINTER :: S17               ! SOLENOID SIXTRACK STYLE
   TYPE(RCOLP), POINTER :: RCOL18             ! RCOLLIMATOR
   TYPE(ECOLP), POINTER :: ECOL19             ! ECOLLIMATOR
   TYPE(CAV_TRAVP), POINTER :: CAV21            ! CAVITY TRAVELLING WAVE
   TYPE(USER1P), POINTER :: U1                ! USER DEFINED
   TYPE(USER2P), POINTER :: U2                ! USER DEFINED
   TYPE(SAGANP), POINTER :: WI                ! USER DEFINED
   TYPE(MARTIN), POINTER :: M22               ! TAYLOR MAP
   TYPE(MULTIPP), POINTER :: MU               ! MULTIP ZGOUBI
   TYPE(layout), POINTER :: g23
   TYPE(FIBRE), POINTER :: PARENT_FIBRE
END TYPE  ELEMENTP


type info
   !   character(nlp),pointer :: name
   real(sp),pointer :: s
   real(sp),pointer ::  beta(:)
   real(sp),pointer ::  fix0(:)
   real(sp),pointer ::  fix(:)
   real(sp), pointer:: pos(:)
   !    type(madx_data), pointer :: d
END type info

TYPE FIBRE
   !  BELOW ARE THE DATA CARRIED BY THE NODE
   INTEGER,POINTER ::DIR
   !     REAL(DP),POINTER ::P0C,BETA0
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
   !   TYPE (FIBRE),POINTER ::  PARENT_PATCH
   !   TYPE (FIBRE),POINTER ::  PARENT_CHART
   !   TYPE (FIBRE),POINTER ::  PARENT_MAG
   type(info),pointer ::i
END TYPE FIBRE

TYPE LAYOUT
   CHARACTER(120), POINTER ::  NAME ! IDENTIFICATION
   INTEGER, POINTER ::  INDEX,CHARGE,HARMONIC_NUMBER ! IDENTIFICATION, CHARGE SIGN
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
   TYPE (MAD_UNIVERSE), POINTER :: parent_universe
END TYPE LAYOUT

TYPE MAD_UNIVERSE
   INTEGER,  POINTER :: N     ! NUMBERS OF LAYOUT
   INTEGER,  POINTER :: SHARED     ! NUMBERS OF REPOSITORY LAYOUTS
   TYPE (LAYOUT), POINTER :: END
   TYPE (LAYOUT), POINTER :: START
END TYPE MAD_UNIVERSE
