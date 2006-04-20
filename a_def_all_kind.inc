
! New stuff SixTrack90-like
TYPE DRIFT1
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp) , POINTER ::L                    !  Storage space for computed
END  TYPE DRIFT1

TYPE DRIFT1P
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L
END  TYPE DRIFT1P

TYPE DKD2
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE DKD2

TYPE DKD2P
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE DKD2P

TYPE KICKT3
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  DIMENSION(:), POINTER :: AN,BN                 !Multipole component
   real(dp),   POINTER:: thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
END  TYPE KICKT3

TYPE KICKT3P
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),   POINTER:: thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
END  TYPE KICKT3P

TYPE CAV4
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp) , POINTER ::L                    !  Storage space for computed
   real(dp), POINTER :: VOLT, FREQ,PHAS,DELTA_E
   integer, POINTER :: N_BESSEL                        ! Cavity information
END  TYPE CAV4

TYPE CAV4P
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8) , POINTER ::L
   TYPE(REAL_8), POINTER :: VOLT, FREQ,PHAS      ! Cavity Information
   real(dp), POINTER :: DELTA_E                        ! Cavity information
   integer, POINTER :: N_BESSEL                        ! Cavity information
END  TYPE CAV4P

TYPE CAV_TRAV
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp) , POINTER ::L                    !  Storage space for computed
   real(dp), POINTER :: VOLT, FREQ,PHAS,DELTA_E,DVDS
   real(dp), POINTER :: PSI,DPHAS
END  TYPE CAV_TRAV


TYPE CAV_TRAVP
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8) , POINTER ::L                    !  Storage space for computed
   TYPE(REAL_8), POINTER :: VOLT, FREQ,PHAS,DVDS
   TYPE(REAL_8), POINTER :: PSI,DPHAS
   real(dp), POINTER :: DELTA_E                        ! Cavity information
END  TYPE CAV_TRAVP



TYPE SOL5
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp), POINTER ::B_SOL
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE SOL5

TYPE SOL5P
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8), POINTER ::B_SOL
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE SOL5P

TYPE SOLT
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp), POINTER ::B_SOL
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MAT    !LINEAR MATRIX
   !frs     real(dp),  DIMENSION(:), POINTER :: LXY(:)
   real(dp),  DIMENSION(:), POINTER :: LXY
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE SOLT

TYPE SOLTP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8), POINTER ::B_SOL
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MAT
   !frs     TYPE(REAL_8), DIMENSION(:), POINTER :: lxy(:)
   TYPE(REAL_8), DIMENSION(:), POINTER :: lxy
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE SOLTP


TYPE KTK
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   !frs     real(dp),  DIMENSION(:), POINTER :: lx(:),ly(:)
   real(dp),  DIMENSION(:), POINTER :: lx,ly
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE KTK

TYPE KTKP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:), POINTER :: lx(:),ly(:)
   TYPE(REAL_8), DIMENSION(:), POINTER :: lx,ly
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE KTKP

TYPE TKTF
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   !frs     real(dp), DIMENSION(:), POINTER :: lx(:)
   real(dp), DIMENSION(:), POINTER :: lx
   real(dp), DIMENSION(:,:), POINTER :: RMATX,RMATY !LINEAR MATRIX
   !frs     real(dp), DIMENSION(:), POINTER ::   Rlx(:)
   real(dp), DIMENSION(:), POINTER ::   Rlx
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE TKTF

TYPE TKTFP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:), POINTER :: lx(:)
   TYPE(REAL_8), DIMENSION(:), POINTER :: lx
   TYPE(REAL_8), DIMENSION(:,:), POINTER :: RMATX,RMATY !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:),   POINTER :: Rlx(:)
   TYPE(REAL_8), DIMENSION(:),   POINTER :: Rlx
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE TKTFP

TYPE NSMI
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  DIMENSION(:), POINTER :: BN                 !Multipole component
END  TYPE NSMI

TYPE NSMIP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: BN         !Multipole component
END  TYPE NSMIP

TYPE SSMI
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  DIMENSION(:), POINTER :: AN          !Multipole component
END  TYPE SSMI

TYPE SSMIP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN         !Multipole component
END  TYPE SSMIP

TYPE TEAPOT
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:), POINTER :: bf_x,bf_y     ! B field polynomial
   !     INTEGER,POINTER :: SECTOR_NMUL
   logical(lp),   POINTER :: DRIFTKICK                   ! Split flag
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END TYPE TEAPOT

TYPE TEAPOTP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),  DIMENSION(:), POINTER :: bf_x,bf_y     ! B field polynomial
   !     INTEGER,POINTER :: SECTOR_NMUL
   logical(lp),   POINTER :: DRIFTKICK
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END TYPE TEAPOTP

TYPE MON
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp) , POINTER ::L   !  monitor and instrument of MAD
   real(dp), pointer :: x,y
END  TYPE MON

TYPE MONP
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L
   real(dp), pointer :: x,y
END  TYPE MONP

TYPE RCOL
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp) , POINTER :: L   !  monitor and instrument of MAD
   type(MADX_APERTURE), POINTER :: A
END  TYPE RCOL

TYPE RCOLP
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L
   type(MADX_APERTURE), POINTER :: A
END  TYPE RCOLP

TYPE ECOL
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp) , POINTER ::L   !  monitor and instrument of MAD
   type(MADX_APERTURE), POINTER :: A
END  TYPE ECOL

TYPE ECOLP
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L
   type(MADX_APERTURE), POINTER :: A
END  TYPE ECOLP

TYPE ESEPTUM               !
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp) , POINTER ::L
   real(dp) , POINTER ::VOLT   !  VOLTAGE IN MV/M
   real(dp) , POINTER ::PHAS    !  TILT TO MIMIC E_X AND E_Y OF MAD-X
END  TYPE ESEPTUM

TYPE ESEPTUMP               !
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L   !  monitor and instrument of MAD
   TYPE(REAL_8) , POINTER ::VOLT   !  VOLTAGE IN MV/M
   TYPE(REAL_8) , POINTER ::PHAS    !  TILT TO MIMIC E_X AND E_Y OF MAD-X
END  TYPE ESEPTUMP

TYPE STREX
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   logical(lp), POINTER :: DRIFTKICK,LIKEMAD
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE STREX

TYPE STREXP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   logical(lp),  POINTER :: DRIFTKICK,LIKEMAD
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
END  TYPE STREXP
! fitted magnet stuff

!  NO POLYMORPH NECESSARY: DESCRIBES A TAYLOR MAP AS IN COSY-INFINITY
TYPE MARTIN
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(tree_element), POINTER :: T
   TYPE(tree_element), POINTER :: T_RAD
   TYPE(tree_element), POINTER :: T_REV
   TYPE(tree_element), POINTER :: T_RAD_REV
   real(dp),pointer :: E(:,:)
   real(dp),pointer :: E_REV(:,:)
   logical(lp), POINTER :: deltamap
END  TYPE MARTIN

! GENERAL PANCAKE

TYPE PANCAKE
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  POINTER ::L
   type(tree_element),  POINTER :: B(:)
   type(tree_element),  POINTER :: Ax(:)
   type(tree_element),  POINTER :: Ay(:)
   REAL(DP), POINTER :: SCALE
   REAL(DP), POINTER :: D_IN(:),ANG_IN(:)
   REAL(DP), POINTER :: D_OUT(:),ANG_OUT(:)
END  TYPE PANCAKE

TYPE PANCAKEP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  POINTER ::L
   type(tree_element),  POINTER :: B(:)
   type(tree_element),  POINTER :: Ax(:)
   type(tree_element),  POINTER :: Ay(:)
   TYPE(REAL_8), POINTER :: SCALE
   REAL(DP), POINTER :: D_IN(:),ANG_IN(:)
   REAL(DP), POINTER :: D_OUT(:),ANG_OUT(:)
END  TYPE PANCAKEP
