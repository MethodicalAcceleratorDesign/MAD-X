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
   TYPE (FIBRE),POINTER ::  PARENT_PATCH
   TYPE (FIBRE),POINTER ::  PARENT_CHART
   TYPE (FIBRE),POINTER ::  PARENT_MAG
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
END TYPE LAYOUT

