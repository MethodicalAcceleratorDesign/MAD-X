!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

  TYPE temps_energie
     real(dp),pointer :: temps
     real(dp),pointer :: energie
     real(dp),pointer :: volt(:),phase(:),tc
  END TYPE temps_energie
  
  TYPE acceleration
     integer,pointer :: n,POS,nst 
     real(dp), pointer :: r,de(:),e_in(:)  !, unit_time
     type(work),pointer :: w1,w2      
     TYPE(FIBRE), pointer :: PREVIOUS    
     TYPE(FIBRE), pointer :: NEXT    
     type(temps_energie),pointer :: tableau(:)
     character(255), pointer :: fichier
  END TYPE acceleration


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
   real(dp),   POINTER:: VA,VS         !valishev-like multipole
   integer, pointer :: f
END  TYPE DKD2

TYPE DKD2P
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA,VS         !valishev-like multipole
   integer, pointer :: f
END  TYPE DKD2P

TYPE KICKT3
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  DIMENSION(:), POINTER :: AN,BN                 !Multipole component
   real(dp),   POINTER:: thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
   real(dp),   POINTER:: hf,vf  ! highly illegal additions by frs for fringe
   logical(lp),   POINTER:: patch ! highly illegal additions by frs
   real(dp), POINTER ::B_SOL
   real(dp),   POINTER::  ls
END  TYPE KICKT3

TYPE KICKT3P
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),   POINTER:: thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
   TYPE(REAL_8),   POINTER:: hf,vf  ! highly illegal additions by frs for fringe
   logical(lp),   POINTER:: patch ! highly illegal additions by frs
   TYPE(REAL_8), POINTER ::B_SOL
   real(dp),   POINTER::  ls
END  TYPE KICKT3P

TYPE CAV4
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp) , POINTER ::L                    !  Storage space for computed
   real(dp),  DIMENSION(:), POINTER :: AN,BN                 !Multipole component
   real(dp), POINTER :: VOLT, FREQ,PHAS,DELTA_E,phase0,t
   real(dp), POINTER :: F(:)
   real(dp), POINTER :: A,R
   real(dp), POINTER :: PH(:)
   INTEGER , POINTER :: NF
   integer, POINTER :: N_BESSEL,CAVITY_TOTALPATH     ! Cavity information
   logical, pointer :: always_on
   type(acceleration), pointer :: ACC
   
END  TYPE CAV4

TYPE CAV4P
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8) , POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8), POINTER :: VOLT, FREQ,PHAS      ! Cavity Information
   real(dp), POINTER :: DELTA_E,phase0,t                        ! Cavity information
   TYPE(REAL_8), POINTER :: F(:)
   TYPE(REAL_8), POINTER :: PH(:)
   TYPE(REAL_8), POINTER :: A,R
   INTEGER , POINTER :: NF
   integer, POINTER :: N_BESSEL,CAVITY_TOTALPATH  ! Cavity information
   logical, pointer :: always_on
   type(acceleration), pointer :: ACC
END  TYPE CAV4P

TYPE CAV_TRAV
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp) , POINTER ::L                    !  Storage space for computed
   real(dp), POINTER :: VOLT, FREQ,PHAS,DELTA_E,phase0
   real(dp), POINTER :: PSI,DPHAS,DVDS
   integer, POINTER :: CAVITY_TOTALPATH     ! Cavity information
END  TYPE CAV_TRAV


TYPE CAV_TRAVP
   logical(lp), POINTER :: THIN
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8) , POINTER ::L                    !  Storage space for computed
   TYPE(REAL_8), POINTER :: VOLT, FREQ,PHAS
   TYPE(REAL_8), POINTER :: PSI,DPHAS,DVDS
   real(dp), POINTER :: DELTA_E ,phase0                       ! Cavity information
   integer, POINTER :: CAVITY_TOTALPATH     ! Cavity information
END  TYPE CAV_TRAVP



TYPE SOL5
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp), POINTER ::B_SOL
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA,VS         !valishev-like multipole
END  TYPE SOL5

TYPE SOL5P
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8), POINTER ::B_SOL
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA,VS         !valishev-like multipole
END  TYPE SOL5P


TYPE KTK
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   !frs     real(dp),  DIMENSION(:), POINTER :: lx(:),ly(:)
   real(dp),  DIMENSION(:), POINTER :: lx,ly
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA,VS         !valishev-like multipole
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
   TYPE(REAL_8),   POINTER:: VA,VS         !valishev-like multipole
END  TYPE KTKP

TYPE TKTF
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   real(dp),  DIMENSION(:,:), POINTER :: MATX2,MATY2 !LINEAR MATRIX
   !frs     real(dp), DIMENSION(:), POINTER :: lx(:)
   real(dp), DIMENSION(:), POINTER :: lx
   real(dp), DIMENSION(:,:), POINTER :: RMATX,RMATY !LINEAR MATRIX
   !frs     real(dp), DIMENSION(:), POINTER ::   Rlx(:)
   real(dp), DIMENSION(:), POINTER ::   Rlx
   !     real(dp), DIMENSION(:), POINTER ::   dx(:)
   !     real(dp), DIMENSION(:), POINTER ::   dy(:)
   real(dp),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA,VS         !valishev-like multipole
   integer, pointer :: f
END  TYPE TKTF

TYPE TKTFP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX,MATY !LINEAR MATRIX
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX2,MATY2 !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:), POINTER :: lx(:)
   TYPE(REAL_8), DIMENSION(:), POINTER :: lx
   TYPE(REAL_8), DIMENSION(:,:), POINTER :: RMATX,RMATY !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:),   POINTER :: Rlx(:)
   TYPE(REAL_8), DIMENSION(:),   POINTER :: Rlx
   !     real(dp), DIMENSION(:), POINTER ::   dx(:)
   !     real(dp), DIMENSION(:), POINTER ::   dy(:)
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA,VS         !valishev-like multipole
   integer, pointer :: f
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
   integer, pointer :: f
   !!  electric info
   real(dp),  DIMENSION(:), POINTER :: AE,BE  ! an,bn for electric stuff
   real(dp),  DIMENSION(:,:,:), POINTER :: As,Bs
   real(dp),  POINTER :: e_x,e_y,PHI
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
   integer, pointer :: f
   !!  electric info
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AE,BE ! an,bn for electric stuff
   real(dp),  DIMENSION(:,:,:), POINTER :: As,Bs
   TYPE(REAL_8),  POINTER :: e_x,e_y,PHI
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
   real(dp) , POINTER :: L   !  collimator of MAD
  ! type(MADX_APERTURE), POINTER :: A
END  TYPE RCOL

TYPE RCOLP
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L
 !  type(MADX_APERTURE), POINTER :: A
END  TYPE RCOLP

TYPE ECOL
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp) , POINTER ::L   !  collimator of MAD
 !  type(MADX_APERTURE), POINTER :: A
END  TYPE ECOL

TYPE ECOLP
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8) , POINTER ::L
  ! type(MADX_APERTURE), POINTER :: A
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
   real(dp),   POINTER:: VA,VS         !valishev-like multipole
   integer, pointer :: f
END  TYPE STREX

TYPE STREXP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   logical(lp),  POINTER :: DRIFTKICK,LIKEMAD
   TYPE(REAL_8),   POINTER:: FINT,HGAP         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1,H2             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA,VS         !valishev-like multipole
   integer, pointer :: f
END  TYPE STREXP
! fitted magnet stuff

TYPE ENGE
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L,D
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp),  DIMENSION(:), POINTER :: A             !ENGE  COEFFICIENTS
   integer,  pointer :: nbessel                      !
   type(my_1D_taylor), POINTER :: F
END  TYPE ENGE

TYPE ENGEP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp), POINTER ::D
   real(dp),  DIMENSION(:), POINTER :: A             !ENGE  COEFFICIENTS
   integer,  pointer :: nbessel                      !
   type(my_1D_taylor), POINTER :: F
END  TYPE ENGEP



! GENERAL PANCAKE

TYPE PANCAKE
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp),  POINTER ::L
   type(tree_element),  POINTER :: B(:)
   !   type(tree_element),  POINTER :: Ax(:)
   !   type(tree_element),  POINTER :: Ay(:)
   REAL(DP), POINTER :: SCALE
   !   REAL(DP), POINTER :: D_IN(:),ANG_IN(:)
   !   REAL(DP), POINTER :: D_OUT(:),ANG_OUT(:)
END  TYPE PANCAKE

TYPE PANCAKEP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8),  POINTER ::L
   type(tree_element),  POINTER :: B(:)
   !   type(tree_element),  POINTER :: Ax(:)
   !   type(tree_element),  POINTER :: Ay(:)
   TYPE(REAL_8), POINTER :: SCALE
   !   REAL(DP), POINTER :: D_IN(:),ANG_IN(:)
   !   REAL(DP), POINTER :: D_OUT(:),ANG_OUT(:)
END  TYPE PANCAKEP


TYPE HELICAL_DIPOLE
   TYPE(MAGNET_CHART), POINTER:: P
   real(dp), POINTER ::L
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
   real(dp), POINTER :: FREQ,PHAS
   integer, POINTER :: N_BESSEL
END  TYPE HELICAL_DIPOLE

TYPE HELICAL_DIPOLEP
   TYPE(MAGNET_CHART), POINTER:: P
   TYPE(REAL_8), POINTER ::L
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN
   TYPE(REAL_8),  POINTER :: FREQ,PHAS
   integer, POINTER :: N_BESSEL
END  TYPE HELICAL_DIPOLEP
