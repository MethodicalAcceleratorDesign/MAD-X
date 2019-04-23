!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

  TYPE temps_energie
     real(dp),pointer :: temps=> null()
     real(dp),pointer :: energie=> null()
     real(dp),pointer :: volt(:)=> null(),phase(:)=> null(),tc=> null()
  END TYPE temps_energie
  
  TYPE acceleration
     integer,pointer :: n=> null(),POS=> null(),nst => null()
     real(dp), pointer :: r=> null(),de(:)=> null(),e_in(:)=> null()  !, unit_time
     type(work),pointer :: w1=> null() 
     type(work),pointer :: w2 => null()
     TYPE(FIBRE), pointer :: PREVIOUS => null()
     TYPE(FIBRE), pointer :: NEXT  => null()   
     type(temps_energie),pointer :: tableau(:) => null()
     character(255), pointer :: fichier => null()
  END TYPE acceleration


TYPE DRIFT1
   TYPE(MAGNET_CHART), POINTER :: P => null()
   real(dp) , POINTER ::L => null()                   !  Storage space for computed  
END  TYPE DRIFT1

TYPE DRIFT1P
   TYPE(MAGNET_CHART), POINTER :: P => null()
   TYPE(REAL_8) , POINTER ::L => null()
END  TYPE DRIFT1P

TYPE SUPERDRIFT
   TYPE(MAGNET_CHART), POINTER :: P => null()
   real(dp) , POINTER ::L => null()                    !  Storage space for computed
   REAL(DP),DIMENSION(:), POINTER:: D => null(),ang => null()     ! TRANSLATIONS  A_D(3) and Angles
   INTEGER, POINTER :: A_X1 => null() ,A_X2 => null()
END  TYPE SUPERDRIFT

TYPE SUPERDRIFTP
   TYPE(MAGNET_CHART), POINTER :: P => null()
   TYPE(REAL_8) , POINTER ::L => null()
   REAL(DP),DIMENSION(:), POINTER:: D => null(),ang => null()     ! TRANSLATIONS  A_D(3) and Angles
   INTEGER, POINTER :: A_X1 => null() ,A_X2 => null()
END  TYPE SUPERDRIFTP

TYPE DKD2
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(),BN => null()         !Multipole component
   real(dp), DIMENSION(:),   POINTER:: FINT => null(),HGAP => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(),H2 => null()             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA => null(),VS => null()         !valishev-like multipole
   integer, pointer :: f => null()
END  TYPE DKD2

TYPE DKD2P
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(),BN => null()         !Multipole component
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(),HGAP => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(),H2 => null()             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA => null(),VS => null()         !valishev-like multipole
   integer, pointer :: f => null()
END  TYPE DKD2P

TYPE KICKT3
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(),BN => null()                 !Multipole component
   real(dp),   POINTER:: thin_h_foc => null(),thin_v_foc => null(),thin_h_angle => null(),thin_v_angle => null()  ! highly illegal additions by frs
   real(dp),   POINTER:: hf => null(),vf => null()  ! highly illegal additions by frs for fringe
   logical(lp),   POINTER:: patch => null() ! highly illegal additions by frs
   real(dp), POINTER ::B_SOL => null()
   real(dp),   POINTER::  ls => null()
   real(dp),   POINTER:: dx => null(),dy => null(),pitch_x => null(),pitch_y => null()   ! SADISTIC
END  TYPE KICKT3

TYPE KICKT3P
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(),BN => null()         !Multipole component
   TYPE(REAL_8),   POINTER:: thin_h_foc => null(),thin_v_foc => null(),thin_h_angle => null(),thin_v_angle => null()  ! highly illegal additions by frs
   TYPE(REAL_8),   POINTER:: hf => null(),vf => null()  ! highly illegal additions by frs for fringe
   logical(lp),   POINTER:: patch => null() ! highly illegal additions by frs
   TYPE(REAL_8), POINTER ::B_SOL => null()
   real(dp),   POINTER::  ls => null()
   real(dp),   POINTER:: dx => null(),dy => null(),pitch_x => null(),pitch_y => null()   ! SADISTIC
END  TYPE KICKT3P

TYPE ABELL
   logical(lp), POINTER :: THIN => null()
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp) , POINTER ::L => null()                    !  Storage space for computed
   real(dp),  DIMENSION(:), POINTER :: AN => null(),BN => null()                 ! MULTIPLICATIVE FACTOR USUALLY 1
   real(dp), POINTER :: DZ(:) => null(), T(:) => null(), TE(:) => null() !, phi0(:) => null()
   complex(dp), POINTER :: B(:,:) => null(),E(:,:) => null()
   INTEGER , POINTER :: N,M => null()  
   REAL(DP), POINTER :: SCALE=> null(),angc=> null(),dc=> null(),hc=> null(),xc => null(),vc => null()
   logical,pointer :: xprime=> null()
END  TYPE ABELL

TYPE ABELLP
   logical(lp), POINTER :: THIN => null()
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()                    !  Storage space for computed
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(),BN => null()  
   real(dp), POINTER :: DZ(:) => null(),T(:) => null(), TE(:) => null()  !, phi0(:) => null()
   complex(dp), POINTER :: B(:,:) => null(),E(:,:) => null()
   INTEGER , POINTER :: N,M => null()  
   real(dp),  POINTER :: angc=> null(),dc=> null(),hc=> null(),xc => null(),vc => null()
   TYPE(REAL_8),  POINTER :: SCALE=> null()
   logical,pointer :: xprime=> null()
END  TYPE ABELLP

TYPE CAV4
   logical(lp), POINTER :: THIN => null()
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp) , POINTER ::L => null()                    !  Storage space for computed
   real(dp),  DIMENSION(:), POINTER :: AN => null(),BN => null()                 !Multipole component
   real(dp), POINTER :: VOLT => null(), FREQ => null(),PHAS => null(),DELTA_E => null(),phase0 => null(),t => null()
   real(dp), POINTER :: F(:) => null()
   real(dp), POINTER :: A => null(),R => null()
   real(dp), POINTER :: PH(:) => null()
   INTEGER , POINTER :: NF => null()
   integer, POINTER :: N_BESSEL => null(),CAVITY_TOTALPATH => null()     ! Cavity information
   logical, pointer :: always_on => null()
   type(acceleration), pointer :: ACC => null()
   
END  TYPE CAV4

TYPE CAV4P
   logical(lp), POINTER :: THIN => null()
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8) , POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(),BN => null()         !Multipole component
   TYPE(REAL_8), POINTER :: VOLT => null(), FREQ => null(),PHAS => null()      ! Cavity Information
   real(dp), POINTER :: DELTA_E => null(),phase0 => null(),t => null()                        ! Cavity information
   TYPE(REAL_8), POINTER :: F(:) => null()
   TYPE(REAL_8), POINTER :: PH(:) => null()
   TYPE(REAL_8), POINTER :: A => null(),R => null()
   INTEGER , POINTER :: NF => null()
   integer, POINTER :: N_BESSEL => null(),CAVITY_TOTALPATH => null()  ! Cavity information
   logical, pointer :: always_on => null()
   type(acceleration), pointer :: ACC => null()
END  TYPE CAV4P

TYPE CAV_TRAV
   logical(lp), POINTER :: THIN => null()
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp) , POINTER ::L => null()                    !  Storage space for computed
   real(dp), POINTER :: VOLT => null(), FREQ => null(),PHAS,DELTA_E => null(),phase0 => null()
   real(dp), POINTER :: PSI => null(),DPHAS => null(),DVDS => null()
   logical, pointer :: always_on => null()
   integer, POINTER :: CAVITY_TOTALPATH => null()     ! Cavity information
END  TYPE CAV_TRAV


TYPE CAV_TRAVP
   logical(lp), POINTER :: THIN => null()
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8) , POINTER ::L  => null()                   !  Storage space for computed
   TYPE(REAL_8), POINTER :: VOLT => null(), FREQ => null(),PHAS => null()
   TYPE(REAL_8), POINTER :: PSI => null(),DPHAS => null(),DVDS => null()
   real(dp), POINTER :: DELTA_E  => null(),phase0 => null()
   logical, pointer :: always_on => null()
   integer, POINTER :: CAVITY_TOTALPATH  => null()    ! Cavity information
END  TYPE CAV_TRAVP



TYPE SOL5
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null()
   real(dp), POINTER ::B_SOL => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(),BN => null()         !Multipole component
   real(dp), DIMENSION(:),    POINTER:: FINT => null(),HGAP => null()        !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(),H2 => null()             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA => null(),VS  => null()        ! sad f1,f2
   real(dp),   POINTER:: dx => null(),dy => null(),pitch_x => null(),pitch_y => null()   ! SADISTIC
END  TYPE SOL5

TYPE SOL5P
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8), POINTER ::B_SOL => null()
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   real(dp),   POINTER:: dx => null(),dy => null(),pitch_x => null(),pitch_y => null()   ! SADISTIC
END  TYPE SOL5P


TYPE KTK
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MATX => null(), MATY => null() !LINEAR MATRIX
   !frs     real(dp),  DIMENSION(:), POINTER :: lx(:) => null(), ly(:) => null()
   real(dp),  DIMENSION(:), POINTER :: lx => null(), ly => null()
   real(dp), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
END  TYPE KTK

TYPE KTKP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX => null(), MATY => null() !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:), POINTER :: lx(:) => null(), ly(:) => null()
   TYPE(REAL_8), DIMENSION(:), POINTER :: lx => null(), ly => null()
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
END  TYPE KTKP

TYPE TKTF
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   real(dp),  DIMENSION(:,:), POINTER :: MATX => null(), MATY => null() !LINEAR MATRIX
   real(dp),  DIMENSION(:,:), POINTER :: MATX2 => null(), MATY2 => null() !LINEAR MATRIX
   !frs     real(dp), DIMENSION(:), POINTER :: lx(:) => null()
   real(dp), DIMENSION(:), POINTER :: lx => null()
   real(dp), DIMENSION(:,:), POINTER :: RMATX => null(), RMATY => null() !LINEAR MATRIX
   !frs     real(dp), DIMENSION(:), POINTER ::   Rlx(:) => null()
   real(dp), DIMENSION(:), POINTER ::   Rlx => null()
   !     real(dp), DIMENSION(:), POINTER ::   dx(:) => null()
   !     real(dp), DIMENSION(:), POINTER ::   dy(:) => null()
   real(dp), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   integer, pointer :: f => null()
END  TYPE TKTF

TYPE TKTFP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX => null(), MATY => null() !LINEAR MATRIX
   TYPE(REAL_8),  DIMENSION(:,:), POINTER :: MATX2 => null(), MATY2 => null() !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:), POINTER :: lx(:) => null()
   TYPE(REAL_8), DIMENSION(:), POINTER :: lx => null()
   TYPE(REAL_8), DIMENSION(:,:), POINTER :: RMATX => null(), RMATY => null() !LINEAR MATRIX
   !frs     TYPE(REAL_8), DIMENSION(:),   POINTER :: Rlx(:) => null()
   TYPE(REAL_8), DIMENSION(:),   POINTER :: Rlx => null()
   !     real(dp), DIMENSION(:), POINTER ::   dx(:) => null()
   !     real(dp), DIMENSION(:), POINTER ::   dy(:) => null()
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   integer, pointer :: f => null()
END  TYPE TKTFP

TYPE NSMI
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp),  DIMENSION(:), POINTER :: BN => null()                 !Multipole component
END  TYPE NSMI

TYPE NSMIP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: BN => null()         !Multipole component
END  TYPE NSMIP

TYPE SSMI
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp),  DIMENSION(:), POINTER :: AN  => null()         !Multipole component
END  TYPE SSMI

TYPE SSMIP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null()         !Multipole component
END  TYPE SSMIP


TYPE TEAPOT
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp),  POINTER ::L => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   real(dp),  DIMENSION(:), POINTER :: bf_x => null(),bf_y => null()     ! B field polynomial
   !     INTEGER,POINTER :: SECTOR_NMUL => null()
   logical(lp),   POINTER :: DRIFTKICK => null()                ! Split flag
   real(dp), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   integer, pointer :: f => null()
   real(dp),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   !!  electric info
   real(dp),  DIMENSION(:), POINTER :: AE => null(), BE => null()  ! an,bn for electric stuff
   real(dp),  DIMENSION(:),  POINTER :: e_x => null(),e_y => null(),PHI => null(),vm => null()
   LOGICAL(LP), POINTER :: ELECTRIC => null()
END TYPE TEAPOT

TYPE TEAPOTP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8),  POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   TYPE(REAL_8),  DIMENSION(:), POINTER :: bf_x => null(),bf_y => null()     ! B field polynomial
   !     INTEGER,POINTER :: SECTOR_NMUL => null()
   logical(lp),   POINTER :: DRIFTKICK => null()
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   integer, pointer :: f => null()
   TYPE(REAL_8),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   !!  electric info
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AE => null(), BE => null()  ! an,bn for electric stuff
   TYPE(REAL_8),  DIMENSION(:),  POINTER :: e_x => null(),e_y => null(),PHI => null(),vm => null()
   LOGICAL(LP), POINTER :: ELECTRIC => null()
END TYPE TEAPOTP

TYPE MON
   TYPE(MAGNET_CHART), POINTER :: P => null()
   real(dp) , POINTER ::L => null()   !  monitor and instrument of MAD
   real(dp), pointer :: x => null(), y => null()
END  TYPE MON

TYPE MONP
   TYPE(MAGNET_CHART), POINTER :: P => null()
   TYPE(REAL_8) , POINTER ::L => null()
   real(dp), pointer :: x => null(), y => null()
END  TYPE MONP

TYPE RCOL
   TYPE(MAGNET_CHART), POINTER :: P => null()
   real(dp) , POINTER :: L => null()   !  collimator of MAD
  ! type(MADX_APERTURE), POINTER :: A => null()
END  TYPE RCOL

TYPE RCOLP
   TYPE(MAGNET_CHART), POINTER :: P => null()
   TYPE(REAL_8) , POINTER ::L => null()
 !  type(MADX_APERTURE), POINTER :: A => null()
END  TYPE RCOLP

TYPE ECOL
   TYPE(MAGNET_CHART), POINTER :: P => null()
   real(dp) , POINTER ::L => null()   !  collimator of MAD
 !  type(MADX_APERTURE), POINTER :: A => null()
END  TYPE ECOL

TYPE ECOLP
   TYPE(MAGNET_CHART), POINTER :: P => null()
   TYPE(REAL_8) , POINTER ::L => null()
  ! type(MADX_APERTURE), POINTER :: A => null()
END  TYPE ECOLP

TYPE ESEPTUM               !
   TYPE(MAGNET_CHART), POINTER :: P => null()
   real(dp) , POINTER ::L => null()
   real(dp) , POINTER ::VOLT => null()   !  VOLTAGE IN MV/M
   real(dp) , POINTER ::PHAS => null()    !  TILT TO MIMIC E_X AND E_Y OF MAD-X
END  TYPE ESEPTUM

TYPE ESEPTUMP               !
   TYPE(MAGNET_CHART), POINTER :: P => null()
   TYPE(REAL_8) , POINTER ::L  => null()  !  monitor and instrument of MAD
   TYPE(REAL_8) , POINTER ::VOLT => null()   !  VOLTAGE IN MV/M
   TYPE(REAL_8) , POINTER ::PHAS => null()   !  TILT TO MIMIC E_X AND E_Y OF MAD-X
END  TYPE ESEPTUMP

TYPE STREX
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   logical(lp), POINTER :: DRIFTKICK => null(), LIKEMAD => null()
   real(dp), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   real(dp),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   integer, pointer :: f => null()
END  TYPE STREX

TYPE STREXP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   logical(lp),  POINTER :: DRIFTKICK => null(), LIKEMAD => null()
   TYPE(REAL_8), DIMENSION(:),   POINTER:: FINT => null(), HGAP => null()         !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: H1 => null(), H2 => null()             !FRINGE FUDGE FOR MAD
   TYPE(REAL_8),   POINTER:: VA => null(), VS => null()         !valishev-like multipole
   integer, pointer :: f => null()
END  TYPE STREXP
! fitted magnet stuff

TYPE ENGE
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null(), D => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   real(dp),  DIMENSION(:), POINTER :: A => null()             !ENGE  COEFFICIENTS
   integer,  pointer :: nbessel => null()                      !
   type(my_1D_taylor), POINTER :: F => null()
END  TYPE ENGE

TYPE ENGEP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()         !Multipole component
   real(dp), POINTER ::D => null()
   real(dp),  DIMENSION(:), POINTER :: A => null()             !ENGE  COEFFICIENTS
   integer,  pointer :: nbessel => null()                      !
   type(my_1D_taylor), POINTER :: F => null()
END  TYPE ENGEP



! GENERAL PANCAKE

TYPE PANCAKE
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp),  POINTER ::L => null()
   type(tree_element),  POINTER :: B(:) => null()
   REAL(DP), POINTER :: SCALE=> null(),angc=> null(),dc=> null(),hc=> null(),xc => null(),vc => null()
   logical,pointer :: xprime=> null()
END  TYPE PANCAKE

TYPE PANCAKEP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8),  POINTER ::L => null()
   type(tree_element),  POINTER :: B(:) => null()
   TYPE(REAL_8), POINTER :: SCALE => null()
   real(dp),pointer :: angc=> null(),dc=> null(),hc=> null(),xc => null(),vc => null()
   logical,pointer :: xprime=> null()
END  TYPE PANCAKEP


TYPE HELICAL_DIPOLE
   TYPE(MAGNET_CHART), POINTER:: P => null()
   real(dp), POINTER ::L => null()
   real(dp),  DIMENSION(:), POINTER :: AN => null(), BN => null(),fake_shift=> null()         !Multipole component
   real(dp), POINTER :: FREQ => null(), PHAS => null()
   integer, POINTER :: N_BESSEL => null()
END  TYPE HELICAL_DIPOLE

TYPE HELICAL_DIPOLEP
   TYPE(MAGNET_CHART), POINTER:: P => null()
   TYPE(REAL_8), POINTER ::L => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN => null(), BN => null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: fake_shift=> null()   
   TYPE(REAL_8),  POINTER :: FREQ => null(), PHAS => null()
   integer, POINTER :: N_BESSEL => null()
END  TYPE HELICAL_DIPOLEP
