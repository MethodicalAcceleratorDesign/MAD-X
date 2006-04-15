!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90

module Mad_like
  USE S_TRACKING
  !USE file_handler
  IMPLICIT NONE
  public

  private QUADTILT, SOLTILT, EL_Q,EL_0,arbitrary_tilt
  private drft,mark,r_r !,rot
  PRIVATE SEXTTILT,OCTUTILT
  private HKICKTILT,VKICKTILT,GKICKTILT
  private GBTILT,SBTILT,pottilt,Set_mad_v
  PRIVATE RFCAVITYL,SMITILT,CHECKSMI,TWCAVITYL
  PRIVATE rectaETILT,recttilt
  PRIVATE B1,A1,A2,B2,A3,B3,A4,B4,A5,A6,A7,A8,A9,A10,B5,B6,B7,B8,B9,B10,BLTILT
  private fac
  private USER_1L,USER_2L,Taylor_maptilt
  PRIVATE MONIT,HMONIT,VMONIT,INSTRUMEN
  PRIVATE RCOLIT,ECOLIT
  ! linked
  private ADD_EE,EQUAL_L_L,add_Eb,add_BE,add_BB,MUL_B,mul_e,SUB_BB,makeitc,makeits
  private unary_subb
  PRIVATE GET_GAM
  logical(lp),PRIVATE ::  MADX= .FALSE.,MADX_MAGNET_ONLY=.FALSE.

  logical(lp),private::LIKEMAD =.false.,mad_list_killed =.true.,setmad = .false.,verbose=.FALSE.,&
       madkick=.false.,circular=.false.,makeit=.false.
  logical(lp)::DRIFT_KICK =.true.
  logical(lp),TARGET ::FIBRE_flip=.true.
  !  logical(lp) :: FIBRE_SURVEY=.true.
  INTEGER,TARGET ::FIBRE_DIR=1
  real(dp),PRIVATE::ENERGY,P0C,BRHO,KINETIC,gamma0I,gamBET,beta0,MC2
  !real(dp),PRIVATE::TOTAL_EPS
  character(80) file_fitted
  !  type(layout),save::mad_list
  type(layout),private::mad_list
  LOGICAL(LP) :: CURVED_ELEMENT=.FALSE.  !  TO SET UP BEND_FRINGE CORRECTLY FOR EXACT
  type(tree_element), PRIVATE :: mad_tree,mad_tree_rad
  type(tree_element),PRIVATE :: mad_tree_REV,mad_tree_rad_REV
  LOGICAL(LP) MAD_TREE_DELTAMAP
  logical(lp):: symplectic_print=.false.
  logical(lp):: symplectify=.false.
  integer :: symplectic_order = 0
  REAL(DP) :: symplectic_eps = -one
  REAL(DP)  MAD_TREE_LD , MAD_TREE_ANGLE
  INTEGER, PRIVATE, TARGET :: NPARA
  type(tree_element), private, allocatable :: t_e(:),t_ax(:),t_ay(:)
  TYPE EL_LIST
     real(dp) L,LD,LC,K(NMAX),KS(NMAX)
     real(dp) ang(3),t(3)
     real(dp) angi(3),ti(3)
     integer patchg
     real(dp) T1,T2,B0
     real(dp) volt,freq0,harmon,lag,DELTA_E,BSOL
     real(dp) tilt
     real(dp) FINT,hgap,h1,h2,X_COL,Y_COL
     real(dp) thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
     CHARACTER(120) file
     CHARACTER(120) file_rev
     CHARACTER(nlp) NAME
     CHARACTER(vp) VORNAME
     INTEGER KIND,nmul,nst,method
     LOGICAL(LP) APERTURE_ON
     INTEGER APERTURE_KIND
     REAL(DP) APERTURE_R(2),APERTURE_X,APERTURE_Y
     LOGICAL(LP) KILL_ENT_FRINGE,KILL_EXI_FRINGE,BEND_FRINGE,PERMFRINGE
     REAL(DP) DPHAS,PSI,dvds
     INTEGER N_BESSEL
     !     logical(lp) in,out
  END TYPE EL_LIST

  INTERFACE OPERATOR (+)
     !  linked
     MODULE PROCEDURE add_EE
     MODULE PROCEDURE add_Eb
     MODULE PROCEDURE add_BE
     MODULE PROCEDURE add_BB
  END INTERFACE



  INTERFACE OPERATOR (-)
     !  linked
     MODULE PROCEDURE SUB_BB
     MODULE PROCEDURE UNARY_SUBB
  END INTERFACE

  INTERFACE OPERATOR (*)
     !    linked
     MODULE PROCEDURE MUL_B
     MODULE PROCEDURE MUL_E
  END INTERFACE

  INTERFACE assignment (=)
     MODULE PROCEDURE EL_Q
     MODULE PROCEDURE EL_0
     !  linked
     MODULE PROCEDURE EQUAL_L_L
  end  INTERFACE

  INTERFACE OPERATOR (.ring.)
     MODULE PROCEDURE makeitc
  END INTERFACE

  INTERFACE OPERATOR (.line.)
     MODULE PROCEDURE makeits
  END INTERFACE



  INTERFACE operator (.is.)
     MODULE PROCEDURE r_r
  end  INTERFACE

  INTERFACE operator (.d.)
     MODULE PROCEDURE B1
  end  INTERFACE
  INTERFACE operator (.sd.)
     MODULE PROCEDURE a1
  end  INTERFACE
  INTERFACE operator (.Q.)
     MODULE PROCEDURE B2
  end  INTERFACE
  INTERFACE operator (.sQ.)
     MODULE PROCEDURE a2
  end  INTERFACE
  INTERFACE operator (.S.)
     MODULE PROCEDURE B3
  end  INTERFACE
  INTERFACE operator (.sS.)
     MODULE PROCEDURE a3
  end  INTERFACE
  INTERFACE operator (.O.)
     MODULE PROCEDURE B4
  end  INTERFACE
  INTERFACE operator (.sO.)
     MODULE PROCEDURE a4
  end  INTERFACE
  INTERFACE operator (.dE.)
     MODULE PROCEDURE B5
  end  INTERFACE
  INTERFACE operator (.sDe.)
     MODULE PROCEDURE a5
  end  INTERFACE
  INTERFACE operator (.Do.)
     MODULE PROCEDURE B6
  end  INTERFACE
  INTERFACE operator (.sDo.)
     MODULE PROCEDURE a6
  end  INTERFACE

  INTERFACE operator (.II.)
     MODULE PROCEDURE B1
  end  INTERFACE
  INTERFACE operator (.sII.)
     MODULE PROCEDURE a1
  end  INTERFACE
  INTERFACE operator (.IV.)
     MODULE PROCEDURE B2
  end  INTERFACE
  INTERFACE operator (.sIV.)
     MODULE PROCEDURE a2
  end  INTERFACE
  INTERFACE operator (.VI.)
     MODULE PROCEDURE B3
  end  INTERFACE
  INTERFACE operator (.sVI.)
     MODULE PROCEDURE a3
  end  INTERFACE
  INTERFACE operator (.VIII.)
     MODULE PROCEDURE B4
  end  INTERFACE
  INTERFACE operator (.sVIII.)
     MODULE PROCEDURE a4
  end  INTERFACE
  INTERFACE operator (.X.)
     MODULE PROCEDURE B5
  end  INTERFACE
  INTERFACE operator (.SX.)
     MODULE PROCEDURE a5
  end  INTERFACE
  INTERFACE operator (.XII.)
     MODULE PROCEDURE B6
  end  INTERFACE
  INTERFACE operator (.SXII.)
     MODULE PROCEDURE a6
  end  INTERFACE
  INTERFACE operator (.XIV.)
     MODULE PROCEDURE B7
  end  INTERFACE
  INTERFACE operator (.SXIV.)
     MODULE PROCEDURE a7
  end  INTERFACE
  INTERFACE operator (.XVI.)
     MODULE PROCEDURE B8
  end  INTERFACE
  INTERFACE operator (.SXVI.)
     MODULE PROCEDURE a8
  end  INTERFACE
  INTERFACE operator (.XVIII.)
     MODULE PROCEDURE B9
  end  INTERFACE
  INTERFACE operator (.SXVIII.)
     MODULE PROCEDURE a9
  end  INTERFACE
  INTERFACE operator (.XX.)
     MODULE PROCEDURE B10
  end  INTERFACE
  INTERFACE operator (.SXX.)
     MODULE PROCEDURE a10
  end  INTERFACE


  INTERFACE EL_Q_FOR_MADX
     MODULE PROCEDURE EL_Q
  end  INTERFACE

  INTERFACE OCTUPOLE
     MODULE PROCEDURE OCTUTILT
  end  INTERFACE

  INTERFACE SEXTUPOLE
     MODULE PROCEDURE SEXTTILT
  end  INTERFACE

  INTERFACE quadrupole
     MODULE PROCEDURE QUADTILT
  end  INTERFACE

  INTERFACE ZGOUBI_MULTIP
     MODULE PROCEDURE ZGOUBI_MULTIPTILT
  end  INTERFACE

  INTERFACE SOLENOID
     MODULE PROCEDURE SOLTILT
  end  INTERFACE

  INTERFACE SMI
     MODULE PROCEDURE SMITILT
  end  INTERFACE

  INTERFACE SINGLE_LENS
     MODULE PROCEDURE SMITILT
  end  INTERFACE

  INTERFACE multipole_block
     MODULE PROCEDURE BLTILT
  end  INTERFACE


  INTERFACE HKICKER
     MODULE PROCEDURE HKICKTILT
  end  INTERFACE

  INTERFACE VKICKER
     MODULE PROCEDURE VKICKTILT
  end  INTERFACE

  INTERFACE KICKER
     MODULE PROCEDURE GKICKTILT
  end  INTERFACE

  INTERFACE rbend
     !     MODULE PROCEDURE recttilt
     MODULE PROCEDURE rectaETILT
  end  INTERFACE

  INTERFACE sbend
     MODULE PROCEDURE sBtilt
  end  INTERFACE

  INTERFACE Gbend
     MODULE PROCEDURE GBtilt
  end  INTERFACE

  INTERFACE drift
     MODULE PROCEDURE drft
  end  INTERFACE

  INTERFACE marker
     MODULE PROCEDURE mark
  end  INTERFACE

  INTERFACE RCOLLIMATOR
     MODULE PROCEDURE RCOLIT
  end  INTERFACE
  INTERFACE ECOLLIMATOR
     MODULE PROCEDURE ECOLIT
  end  INTERFACE

  INTERFACE MONITOR
     MODULE PROCEDURE MONIT
  end  INTERFACE
  INTERFACE HMONITOR
     MODULE PROCEDURE HMONIT
  end  INTERFACE
  INTERFACE VMONITOR
     MODULE PROCEDURE VMONIT
  end  INTERFACE
  INTERFACE INSTRUMENT
     MODULE PROCEDURE INSTRUMEN
  end  INTERFACE

  INTERFACE RFCAVITY
     MODULE PROCEDURE RFCAVITYL
  end  INTERFACE

  INTERFACE TWCAVITY
     MODULE PROCEDURE TWCAVITYL
  end  INTERFACE

  INTERFACE ELSEPARATOR
     MODULE PROCEDURE ELSESTILT
  end  INTERFACE



  INTERFACE USER_1
     MODULE PROCEDURE USER_1L
  end  INTERFACE
  INTERFACE USER_2
     MODULE PROCEDURE USER_2L
  end  INTERFACE

  INTERFACE WIGGLER
     MODULE PROCEDURE WIGGLERL
  end  INTERFACE



  INTERFACE arbitrary
     MODULE PROCEDURE arbitrary_tilt
  end  INTERFACE

  !  Taylor map
  INTERFACE Taylor_map
     MODULE PROCEDURE  Taylor_maptilt
  end  INTERFACE



CONTAINS

  SUBROUTINE SET_MADX_(CONV,CONV1)
    IMPLICIT NONE
    logical(lp) CONV,CONV1
    MADX=CONV
    MADX_MAGNET_ONLY=CONV1
  END SUBROUTINE SET_MADX_



  FUNCTION r_r( S1, S2 )
    implicit none
    TYPE(TILTING) r_r
    TYPE(TILTING), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2


    r_r=S1
    R_R%TILT(0)=S2
    R_R%NATURAL=.FALSE.

  END FUNCTION r_r

  real(dp) function fac(n)    ! David Sagan
    implicit none
    integer n
    fac=one
    if(mad) then
       fac=madfac(iabs(n))
    endif

  end  function fac

  SUBROUTINE  CHECKSMI(S2,S1)
    implicit none
    type (EL_LIST),INTENT(IN):: S2
    INTEGER,INTENT(IN):: S1
    IF(S2%KIND==KIND8) THEN
       IF(S2%NMUL/=S1) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a24,1x,i4,a21,1x,i4)')  MYTYPE(KIND8),S2%NMUL,' DOES NOT ALLOW POLE ', 2*S1
          call write_e(KIND8)
       ENDIF
    ELSEIF(S2%KIND==KIND9) THEN
       IF(S2%NMUL/=-S1) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a24,1x,i4,a21,1x,i4)') MYTYPE(KIND9),S2%NMUL,' DOES NOT ALLOW POLE ',2*S1
          call write_e(KIND9)
       ENDIF
    ENDIF

  END SUBROUTINE CHECKSMI


  FUNCTION  A10(S2,S1)
    implicit none
    type (EL_LIST) A10
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-10)
    A10 =S2
    A10 %K(10)=A10%K(10)
    A10 %KS(10)=A10%KS(10)+S1/fac(10)
  END FUNCTION A10

  FUNCTION  B10(S2,S1)
    implicit none
    type (EL_LIST) B10
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,10)
    B10 =S2
    B10 %K(10)=B10 %K(10)+S1/fac(10)
    B10 %KS(10)=B10 %KS(10)
  END FUNCTION B10

  FUNCTION  A9(S2,S1)
    implicit none
    type (EL_LIST) A9
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-9)
    A9 =S2
    A9 %K(9)=A9%K(9)
    A9 %KS(9)=A9%KS(9)+S1/fac(9)
  END FUNCTION A9

  FUNCTION  B9(S2,S1)
    implicit none
    type (EL_LIST) B9
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,9)
    B9 =S2
    B9 %K(9)=B9 %K(9)+S1/fac(9)
    B9 %KS(9)=B9 %KS(9)
  END FUNCTION B9

  FUNCTION  A8(S2,S1)
    implicit none
    type (EL_LIST) A8
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-8)
    A8 =S2
    A8 %K(8)=A8%K(8)
    A8 %KS(8)=A8%KS(8)+S1/fac(8)
  END FUNCTION A8

  FUNCTION  B8(S2,S1)
    implicit none
    type (EL_LIST) B8
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,8)
    B8 =S2
    B8 %K(8)=B8 %K(8)+S1/fac(8)
    B8 %KS(8)=B8 %KS(8)
  END FUNCTION B8

  FUNCTION  A7(S2,S1)
    implicit none
    type (EL_LIST) A7
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-7)
    A7 =S2
    A7 %K(7)=A7%K(7)
    A7 %KS(7)=A7%KS(7)+S1/fac(7)
  END FUNCTION A7

  FUNCTION  B7(S2,S1)
    implicit none
    type (EL_LIST) B7
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,7)
    B7 =S2
    B7 %K(7)=B7 %K(7)+S1/fac(7)
    B7 %KS(7)=B7 %KS(7)
  END FUNCTION B7

  FUNCTION  A6(S2,S1)
    implicit none
    type (EL_LIST) A6
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-6)
    A6 =S2
    A6 %K(6)=A6%K(6)
    A6 %KS(6)=A6%KS(6)+S1/fac(6)
  END FUNCTION A6

  FUNCTION  B6(S2,S1)
    implicit none
    type (EL_LIST) B6
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,6)
    B6 =S2
    B6 %K(6)=B6 %K(6)+S1/fac(6)
    B6 %KS(6)=B6 %KS(6)
  END FUNCTION B6

  FUNCTION  A5(S2,S1)
    implicit none
    type (EL_LIST) A5
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-5)
    A5 =S2
    A5 %K(5)=A5%K(5)
    A5 %KS(5)=A5%KS(5)+S1/fac(5)
  END FUNCTION A5

  FUNCTION  B5(S2,S1)
    implicit none
    type (EL_LIST) B5
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,5)
    B5 =S2
    B5 %K(5)=B5 %K(5)+S1/fac(5)
    B5 %KS(5)=B5 %KS(5)
  END FUNCTION B5

  FUNCTION  A4(S2,S1)
    implicit none
    type (EL_LIST) A4
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-4)
    A4 =S2
    A4 %K(4)=A4%K(4)
    A4 %KS(4)=A4%KS(4)+S1/fac(4)
  END FUNCTION A4

  FUNCTION  B4(S2,S1)
    implicit none
    type (EL_LIST) B4
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,4)
    B4 =S2
    B4 %K(4)=B4 %K(4)+S1/fac(4)
    B4 %KS(4)=B4 %KS(4)
  END FUNCTION B4

  FUNCTION  A3(S2,S1)
    implicit none
    type (EL_LIST) A3
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-3)
    A3 =S2
    A3 %K(3)=A3%K(3)
    A3 %KS(3)=A3%KS(3)+S1/fac(3)
  END FUNCTION A3

  FUNCTION  B3(S2,S1)
    implicit none
    type (EL_LIST) B3
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,3)
    B3 =S2
    B3 %K(3)=B3 %K(3)+S1/fac(3)
    B3 %KS(3)=B3 %KS(3)
  END FUNCTION B3

  FUNCTION  A2(S2,S1)
    implicit none
    type (EL_LIST) A2
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,-2)
    A2 =S2
    A2 %K(2)=A2%K(2)
    A2 %KS(2)=A2%KS(2)+S1
  END FUNCTION A2

  FUNCTION  B2(S2,S1)
    implicit none
    type (EL_LIST) B2
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    CALL CHECKSMI(S2,2)
    B2 =S2
    B2 %K(2)=B2 %K(2)+S1
    B2 %KS(2)=B2 %KS(2)
  END FUNCTION B2

  FUNCTION  A1(S2,S1)
    implicit none
    type (EL_LIST) A1
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    real(dp) smad
    CALL CHECKSMI(S2,-1)
    smad=s1
    if(madkick) then
       if(s2%L/=0) smad=smad/s2%L
    endif
    A1 =S2
    A1 %K(1)=A1%K(1)
    A1 %KS(1)=A1%KS(1)+Smad
  END FUNCTION A1

  FUNCTION  B1(S2,S1)
    implicit none
    type (EL_LIST) B1
    type (EL_LIST),INTENT(IN):: S2
    real(dp),INTENT(IN):: S1
    real(dp) smad
    CALL CHECKSMI(S2,1)

    smad=s1
    if(madkick) then
       smad=-smad
       if(s2%L/=0) smad=smad/s2%L
    endif

    B1 =S2
    B1 %K(1)=B1 %K(1)+smad
    B1 %KS(1)=B1 %KS(1)
  END FUNCTION B1




  SUBROUTINE  EL_0(S2,S1)
    implicit none
    type (EL_LIST),INTENT(OUT):: S2
    INTEGER,INTENT(IN):: S1
    INTEGER I

    if(.not.setmad) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1) =  " Run the Set_mad routine first "
       call write_e(-1)
    endif

    IF(S1==0) THEN
       S2%ang=zero
       S2%t=zero
       S2%angi=zero
       S2%ti=zero
       S2%patchg=0
       S2%L=zero
       S2%LD=zero
       S2%LC=zero
       S2%TILT=zero
       DO I=1,NMAX
          S2%K(I)=zero;S2%KS(I)=zero
       ENDDO
       S2%T1=zero
       S2%T2=zero
       S2%B0=zero
       S2%BSOL=zero
       S2%volt=zero
       S2%freq0=zero
       S2%harmon=one
       S2%DELTA_E=zero
       S2%lag=zero
       S2%KIND=0
       S2%nmul=0
       S2%method=metd
       S2%nst=nstd
       s2%NAME=' '
       s2%VORNAME=' '
       s2%file=' '
       s2%file_rev=' '
       s2%FINT=half
       s2%hgap=zero
       s2%h1=zero
       s2%h2=zero
       s2%thin_h_foc=zero
       s2%thin_v_foc=zero
       s2%thin_h_angle=zero
       s2%thin_v_angle=zero
       s2%APERTURE_ON=.FALSE.
       s2%KILL_ENT_FRINGE=.FALSE.
       s2%KILL_EXI_FRINGE=.FALSE.
       s2%BEND_FRINGE=.FALSE.
       s2%PERMFRINGE=.FALSE.
       s2%DPHAS=ZERO
       s2%dvds=ZERO
       s2%PSI=ZERO
       s2%N_BESSEL=0

       s2%APERTURE_KIND=0
       S2%APERTURE_R=absolute_aperture
       S2%APERTURE_X=absolute_aperture
       S2%APERTURE_Y=absolute_aperture
    ENDIF
  END SUBROUTINE EL_0

  !  DEFINING ELEMEMTS

  FUNCTION  SMITILT(NAME,K1,N,T,LIST)
    implicit none
    type (EL_LIST) SMITILT
    type (EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    type (TILTING),optional, INTENT(IN):: T
    real(dp),optional, INTENT(IN):: K1
    INTEGER,optional,INTENT(IN):: N
    INTEGER NN,I
    LOGICAL(LP) SEARCH
    REAL(DP) K11
    NN=0
    K11=zero
    IF(PRESENT(N)) NN=N
    IF(PRESENT(K1)) K11=K1

    IF(PRESENT(LIST)) THEN   !
       SMITILT=LIST    !  SPECIAL SINCE SMI CAN ONLY BE A SINGLE POLE
       SMITILT%L=zero
       SMITILT%LD=zero
       SMITILT%LC=zero
       NN=1
       SEARCH=.TRUE.
       DO I=NMAX,1,-1
          IF(LIST%K(I)/=zero.AND.SEARCH) THEN
             SEARCH=.FALSE.
             K11=LIST%K(I)
             NN=I
          ENDIF
          IF(LIST%KS(I)/=zero.AND.SEARCH) THEN
             SEARCH=.FALSE.
             K11=LIST%KS(I)
             NN=-I
          ENDIF
       ENDDO

       IF(NN>=1.AND.NN<=10) THEN
          SMITILT%K(NN)=K11/fac(nN)
          SMITILT%KIND=kind8
          SMITILT%nmul=NN
       ELSEIF(NN<0.AND.NN>=-10) THEN
          SMITILT%KS(-NN)=K11/fac(nN)
          SMITILT%KIND=kind9
          SMITILT%nmul=-NN
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a21,1x,i4)') " FORBIDDEN 'SMITILT' ",NN
          call write_e(1221)
       ENDIF
       if(present(t)) SMITILT%tilt=t%tilt(0)

       IF(LEN(NAME)>nlp) THEN
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=name
          WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
          call write_i
          SMITILT%NAME=NAME(1:16)
       ELSE
          SMITILT%NAME=NAME
       ENDIF

    ELSE    !
       SMITILT=0
       SMITILT%L=zero
       SMITILT%LD=zero
       SMITILT%LC=zero
       IF(NN>=1.AND.NN<=10) THEN
          SMITILT%K(NN)=K11/fac(Nn)
          SMITILT%KIND=kind8
          SMITILT%nmul=NN
       ELSEIF(NN<0.AND.NN>=-10) THEN
          SMITILT%KS(-NN)=K11/fac(nN)
          SMITILT%KIND=kind9
          SMITILT%nmul=-NN
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a21,1x,i4)') " FORBIDDEN 'SMITILT' ",NN
          call write_e(1221)
       ENDIF
       if(present(t)) then
          IF(T%NATURAL) THEN
             SMITILT%tilt=t%tilt(iabs(Nn))
          ELSE
             SMITILT%tilt=t%tilt(0)
          ENDIF
       endif



       IF(LEN(NAME)>nlp) THEN
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=name
          WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
          call write_i
          SMITILT%NAME=NAME(1:16)
       ELSE
          SMITILT%NAME=NAME
       ENDIF

    ENDIF   !1
  END FUNCTION SMITILT

  FUNCTION  BLTILT(NAME,K,T,LIST)
    implicit none
    type (EL_LIST) BLTILT
    type (EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    type (TILTING),optional, INTENT(IN):: T
    TYPE(MUL_BLOCK),OPTIONAL, INTENT(IN):: K
    INTEGER I
    LOGICAL(LP) COUNT
    if(present(list)) then   !1
       BLTILT=list
       BLTILT%L=zero
       BLTILT%LD=zero
       BLTILT%LC=zero

       BLTILT%KIND=kind3
       BLTILT%nmul=LIST%NMUL
       COUNT=.TRUE.

       DO I=NMAX,1,-1
          BLTILT%K(I)=LIST%K(I)/fac(i)
          BLTILT%KS(I)=LIST%KS(I)/fac(i)
          IF(COUNT) THEN
             IF(BLTILT%K(I)/=zero.OR.BLTILT%KS(I)/=zero) THEN
                COUNT=.FALSE.
                BLTILT%nmul=I
             ENDIF
          ENDIF
       ENDDO

       if(present(t)) BLTILT%tilt=t%tilt(0)



       IF(LEN(NAME)>nlp) THEN
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=name
          WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
          call write_i
          BLTILT%NAME=NAME(1:16)
       ELSE
          BLTILT%NAME=NAME
       ENDIF

    else   !1
       BLTILT=0
       BLTILT%L=zero
       BLTILT%LD=zero
       BLTILT%LC=zero

       BLTILT%KIND=kind3
       BLTILT%nmul=K%NMUL
       DO I=1,K%NMUL
          BLTILT%K(I)=K%BN(I)/fac(i)
          BLTILT%KS(I)=K%AN(I)/fac(i)
       ENDDO

       if(present(t)) then
          IF(T%NATURAL) THEN
             BLTILT%tilt=t%tilt(K%NATURAL)
          ELSE
             BLTILT%tilt=t%tilt(0)
          ENDIF
       endif



       IF(LEN(NAME)>nlp) THEN
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=name
          WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
          call write_i
          BLTILT%NAME=NAME(1:16)
       ELSE
          BLTILT%NAME=NAME
       ENDIF
    endif    !1
  END FUNCTION BLTILT


  FUNCTION  HKICKTILT(NAME,L,kick,T)
    implicit none
    type (EL_LIST) HKICKTILT
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,OPTIONAL, INTENT(IN):: L,kick
    real(dp) L1,K11
    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(kick)) K11=kick
    madkick=.true.
    HKICKTILT=0
    HKICKTILT%L=L1
    HKICKTILT%LD=L1
    HKICKTILT%LC=L1
    IF(L1==zero) THEN
       HKICKTILT%K(1)=-K11        ! MAD convention K1>0 means px > 0
       HKICKTILT%KIND=MADKIND3N
       HKICKTILT%nmul=1
    ELSE
       HKICKTILT%K(1)=-K11/L1
       HKICKTILT%KIND=MADKIND2
       HKICKTILT%nmul=2
    ENDIF

    IF(PRESENT(T)) THEN
       IF(T%NATURAL) THEN
          HKICKTILT%tilt=t%tilt(1)
       ELSE
          HKICKTILT%tilt=t%tilt(0)
       ENDIF
    ENDIF

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       HKICKTILT%NAME=NAME(1:16)
    ELSE
       HKICKTILT%NAME=NAME
    ENDIF
  END FUNCTION HKICKTILT

  FUNCTION  VKICKTILT(NAME,L,kick,T)
    implicit none
    type (EL_LIST) VKICKTILT
    type (TILTING),OPTIONAL, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,OPTIONAL, INTENT(IN):: L,kick
    real(dp) L1,K11
    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(kick)) K11=kick

    madkick=.true.
    VKICKTILT=0
    VKICKTILT%L=L1
    VKICKTILT%LD=L1
    VKICKTILT%LC=L1
    IF(L1==zero) THEN
       VKICKTILT%KS(1)=K11        ! MAD convention K1>0 means px > 0
       VKICKTILT%KIND=MADKIND3S
       VKICKTILT%nmul=1
    ELSE
       VKICKTILT%KS(1)=K11/L1
       VKICKTILT%KIND=MADKIND2
       VKICKTILT%nmul=2
    ENDIF
    IF(PRESENT(T)) THEN
       IF(T%NATURAL) THEN
          VKICKTILT%tilt=t%tilt(1)
       ELSE
          VKICKTILT%tilt=t%tilt(0)
       ENDIF
    ENDIF

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       VKICKTILT%NAME=NAME(1:16)
    ELSE
       VKICKTILT%NAME=NAME
    ENDIF
  END FUNCTION VKICKTILT


  FUNCTION  GKICKTILT(NAME,L,hkick,vkick,T,LIST)
    implicit none
    type (EL_LIST) GKICKTILT
    type (EL_LIST), OPTIONAL,INTENT(IN):: LIST
    type (TILTING), OPTIONAL,INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,OPTIONAL, INTENT(IN):: L ,hkick ,vkick
    real(dp) L1,K11,K21
    L1=zero
    K11=zero
    K21=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(hkick)) K11=hkick
    IF(PRESENT(vkick)) K21=vkick
    madkick=.true.

    if(present(list)) then
       GKICKTILT=list
       l1=list%L
       K11=LIST%K(1)
       K21=LIST%KS(1)


    else
       GKICKTILT=0
    endif
    GKICKTILT%L=L1
    GKICKTILT%LD=L1
    GKICKTILT%LC=L1
    IF(L1==zero) THEN
       GKICKTILT%K(1)=-K11        ! MAD convention K1>0 means px > 0
       GKICKTILT%KS(1)=K21        ! MAD convention K1>0 means px > 0
       GKICKTILT%KIND=KIND3
       GKICKTILT%nmul=1
    ELSE
       GKICKTILT%K(1)=-K11/L1        ! MAD convention K1>0 means px > 0
       GKICKTILT%KS(1)=K21/L1        ! MAD convention K1>0 means px > 0
       GKICKTILT%KIND=MADKIND2
       GKICKTILT%nmul=2
    ENDIF
    IF(PRESENT(T)) THEN      !2002.11.09   BUG
       IF(T%NATURAL) THEN
          GKICKTILT%tilt=t%tilt(1)
       ELSE
          GKICKTILT%tilt=t%tilt(0)
       ENDIF
    ENDIF

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       GKICKTILT%NAME=NAME(1:16)
    ELSE
       GKICKTILT%NAME=NAME
    ENDIF
  END FUNCTION GKICKTILT


  FUNCTION  ZGOUBI_MULTIPTILT(NAME,L,K1,T,list)
    implicit none
    type (EL_LIST) ZGOUBI_MULTIPTILT
    type (EL_LIST),optional, INTENT(IN)::list
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,K1
    real(dp) L1,K11
    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(K1)) K11=K1
    if(present(list)) then
       ZGOUBI_MULTIPTILT=list
       l1=list%L
       K11=LIST%K(2)
    else
       ZGOUBI_MULTIPTILT=0
    endif
    ZGOUBI_MULTIPTILT%L=L1
    ZGOUBI_MULTIPTILT%LD=L1
    ZGOUBI_MULTIPTILT%LC=L1
    ZGOUBI_MULTIPTILT%K(2)=K11
    IF(L1==zero) THEN
       WRITE(6,*) " NOT PERMITTED IN ZGOUBI_MULTIPTILT "
       STOP 22
    ELSE
       ZGOUBI_MULTIPTILT%K(2)=K11
       ZGOUBI_MULTIPTILT%KIND=KINDMU
    ENDIF
    ZGOUBI_MULTIPTILT%nmul=2
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          WRITE(6,*) " NOT PERMITTED IN ZGOUBI_MULTIPTILT "
          STOP 23
       ELSE
          ZGOUBI_MULTIPTILT%tilt=t%tilt(0)
       ENDIF
    endif
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       ZGOUBI_MULTIPTILT%NAME=NAME(1:16)
    ELSE
       ZGOUBI_MULTIPTILT%NAME=NAME
    ENDIF
  END FUNCTION ZGOUBI_MULTIPTILT

  FUNCTION  QUADTILT(NAME,L,K1,T,list)
    implicit none
    type (EL_LIST) QUADTILT
    type (EL_LIST),optional, INTENT(IN)::list
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,K1
    real(dp) L1,K11
    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(K1)) K11=K1
    if(present(list)) then
       quadtilt=list
       l1=list%L
       K11=LIST%K(2)
    else
       QUADTILT=0
    endif
    QUADTILT%L=L1
    QUADTILT%LD=L1
    QUADTILT%LC=L1
    QUADTILT%K(2)=K11
    IF(L1==zero) THEN
       QUADTILT%K(2)=K11
       QUADTILT%KIND=MADKIND3N
    ELSE
       QUADTILT%K(2)=K11
       QUADTILT%KIND=MADKIND2
    ENDIF
    QUADTILT%nmul=2
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          QUADTILT%tilt=t%tilt(2)
       ELSE
          QUADTILT%tilt=t%tilt(0)
       ENDIF
    endif
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       QUADTILT%NAME=NAME(1:16)
    ELSE
       QUADTILT%NAME=NAME
    ENDIF
  END FUNCTION QUADTILT


  FUNCTION  SOLTILT(NAME,L,KS,K1,T,LIST)
    implicit none
    type (EL_LIST) SOLTILT
    type (EL_LIST),optional, INTENT(IN):: LIST
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,KS,K1
    real(dp) L1,K11,kq

    L1=zero
    K11=zero
    KQ=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(KS)) K11=KS
    IF(PRESENT(k1)) kq=K1

    if(present(list)) then
       SOLTILT=list
       l1=list%L
       K11=LIST%BSOL
       KQ=LIST%K(2)
    else
       SOLTILT=0
    endif
    SOLTILT%L=L1
    SOLTILT%LD=L1
    SOLTILT%LC=L1
    SOLTILT%BSOL=K11
    SOLTILT%nmul=2
    IF(L1==zero) THEN
       SOLTILT%KIND=KIND0
    ELSE
       SOLTILT%K(2)=KQ/FAC(2)    ! MAD FACTOR
       if(madkind2==kind2) then
          SOLTILT%KIND=KIND5
       else
          SOLTILT%KIND=KIND17
          !          SOLTILT%nmul=2
          !          SOLTILT%METHOD=2
       endif
    ENDIF
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          SOLTILT%tilt=zero   ! NO NATURAL TILT
       ELSE
          SOLTILT%tilt=t%tilt(0)
       ENDIF
    endif
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       SOLTILT%NAME=NAME(1:16)
    ELSE
       SOLTILT%NAME=NAME
    ENDIF
  END FUNCTION SOLTILT


  FUNCTION  SEXTTILT(NAME,L,K2,T,LIST)
    implicit none
    type (EL_LIST) SEXTTILT
    type (EL_LIST),optional, INTENT(IN)::list
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp),optional , INTENT(IN):: L,K2
    real(dp) L1,K11

    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(K2)) K11=K2
    if(present(list)) then
       SEXTTILT=list
       l1=list%L
       K11=LIST%K(3)
    else
       SEXTTILT=0
    endif
    SEXTTILT%L=L1
    SEXTTILT%LD=L1
    SEXTTILT%LC=L1
    IF(L1==zero) THEN
       SEXTTILT%K(3)=K11/FAC(3)    ! MAD FACTOR
       SEXTTILT%KIND=MADKIND3N
    ELSE
       SEXTTILT%K(3)=K11/FAC(3)         ! MAD FACTOR
       SEXTTILT%KIND=MADKIND2
    ENDIF
    SEXTTILT%nmul=3
    if(present(t)) then
       IF(T%NATURAL) THEN
          SEXTTILT%tilt=t%tilt(3)
       ELSE
          SEXTTILT%tilt=t%tilt(0)
       ENDIF
    endif

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       SEXTTILT%NAME=NAME(1:16)
    ELSE
       SEXTTILT%NAME=NAME
    ENDIF
  END FUNCTION SEXTTILT


  FUNCTION  OCTUTILT(NAME,L,K3,T,LIST)
    implicit none
    type (EL_LIST) OCTUTILT
    type (EL_LIST),optional, INTENT(IN)::list
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,K3
    real(dp) L1,K11
    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(K3)) K11=K3
    if(present(list)) then
       OCTUTILT=list
       l1=list%L
       K11=LIST%K(4)
    else
       OCTUTILT=0
    endif
    OCTUTILT%L=L1
    OCTUTILT%LD=L1
    OCTUTILT%LC=L1
    IF(L1==zero) THEN
       OCTUTILT%K(4)=K11/FAC(4)         ! MAD FACTOR
       OCTUTILT%KIND=MADKIND3N
    ELSE
       OCTUTILT%K(4)=K11/FAC(4)         ! MAD FACTOR
       OCTUTILT%KIND=MADKIND2
    ENDIF
    OCTUTILT%nmul=4
    if(present(t)) then
       IF(T%NATURAL) THEN
          OCTUTILT%tilt=t%tilt(4)
       ELSE
          OCTUTILT%tilt=t%tilt(0)
       ENDIF
    endif
    !  call rot(OCTUTILT%tilt,OCTUTILT%K,OCTUTILT%KS,OCTUTILT%C,OCTUTILT%S)
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       OCTUTILT%NAME=NAME(1:16)
    ELSE
       OCTUTILT%NAME=NAME
    ENDIF
  END FUNCTION OCTUTILT


  FUNCTION  SBTILT(NAME,L,ANGLE,E1,E2,T,LIST)
    implicit none
    type (EL_LIST) SBTILT
    type (EL_LIST),optional, INTENT(IN)::list
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,angle,E1,E2
    real(dp) L1,ANG1,E11,E22
    CURVED_ELEMENT=.TRUE.
    L1=zero
    ANG1=zero
    E11=zero
    E22=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(angle)) ANG1=angle

    IF(PRESENT(E1)) E11=E1
    IF(PRESENT(E2)) E22=E2




    if(present(list)) then
       SBTILT=list
       l1=list%L
       E11=LIST%T1
       E22=LIST%T2
       ANG1=LIST%B0
    else
       SBTILT=0
    endif

    if(present(t))then
       IF(EXACT_MODEL) THEN                 ! .and.madkind2==kind2
          SBTILT=POTTILT(NAME,L1,ANG1,E11,E22,T,LIST)
       ELSE
          SBTILT=GBEND(NAME,L1,ANG1,E11,E22,T,LIST)
       ENDIF
    else
       IF(EXACT_MODEL) THEN                 ! .and.madkind2==kind2
          SBTILT=POTTILT(NAME,L1,ANG1,E11,E22)
       ELSE
          SBTILT=GBEND(NAME,L1,ANG1,E11,E22)
       ENDIF
    endif

  END FUNCTION SBTILT


  FUNCTION  POTTILT(NAME,L,ANG,E1,E2,T,LIST)
    implicit none
    type (EL_LIST) POTTILT
    type (EL_LIST),optional, INTENT(IN)::list
    real(dp) ,optional, INTENT(IN):: E1,E2
    type (TILTING),optional, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp),optional , INTENT(IN):: L,ANG
    real(dp) E11,E22,L1,ANG1

    E11=zero
    E22=zero
    L1=zero
    ANG1=zero
    IF(PRESENT(E1)) E11=E1 ;
    IF(PRESENT(E2)) E22=E2 ;
    IF(PRESENT(ANG)) ANG1=ANG ;
    IF(PRESENT(L)) L1=L ;
    if(present(list)) then
       POTTILT=list
       l1=list%L
       ANG1=LIST%B0
       E11=LIST%T1
       E22=LIST%T2
    else
       POTTILT=0
    endif



    POTTILT%B0=ANG1/L1
    POTTILT%L=L1
    POTTILT%LD=L1
    POTTILT%T1=E11;
    POTTILT%T2=E22;

    IF(ANG/=zero) THEN
       POTTILT%LC=two*SIN(ANG/two)/POTTILT%B0
    ELSE
       POTTILT%LC=POTTILT%L
    ENDIF
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       POTTILT%NAME=NAME(1:16)
    ELSE
       POTTILT%NAME=NAME
    ENDIF

    if(present(t)) then
       IF(T%NATURAL) THEN
          POTTILT%tilt=t%tilt(1)
       ELSE
          POTTILT%tilt=t%tilt(0)
       ENDIF
    endif

    POTTILT%KIND=KIND10
    POTTILT%K(1)=POTTILT%B0+POTTILT%K(1)
    POTTILT%nmul=SECTOR_NMUL

  END FUNCTION POTTILT


  FUNCTION  GBTILT(NAME,L,ANGLE,e1,e2,T,LIST)
    implicit none
    type (EL_LIST) GBTILT
    type (EL_LIST),optional, INTENT(IN)::list
    type (TILTING), optional,INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,angle,e1,e2
    real(dp) L1,ANG1,t11,t21
    if(exact_model) then
       w_p=0
       w_p%nc=5
       w_p%fc='(4(1X,a72,/),(1X,a72))'
       w_p%c(1)= " *************************************************** "
       w_p%c(2)= " * In PTC, under the exact option                  * "
       w_p%c(3)= " * One must distinguish between RBEND and SBEND    * "
       w_p%c(4)= " * This is call is thus completely forbidden       * "
       w_p%c(5)= " *************************************************** "
       call write_e(101)
    endif
    L1=zero
    ANG1=zero
    t11=zero
    t21=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(angle)) ANG1=angle
    IF(PRESENT(e1)) t11=e1
    IF(PRESENT(e2)) t21=e2

    if(present(list)) then
       GBTILT=list
       l1=list%L
       ANG1=LIST%B0
       T11=LIST%T1
       T21=LIST%T2
    else
       GBTILT=0
    endif
    GBTILT%B0=ANG1/L1
    GBTILT%L=L1
    GBTILT%LD=L1
    IF(ANG1/=zero) THEN
       GBTILT%LC=two*SIN(ANG1/two)/GBTILT%B0
    ELSE
       GBTILT%LC=GBTILT%L
    ENDIF
    GBTILT%T1=T11 ; GBTILT%T2=T21;
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       GBTILT%NAME=NAME(1:16)
    ELSE
       GBTILT%NAME=NAME
    ENDIF
    GBTILT%K(1)=GBTILT%B0+GBTILT%K(1)   ! NEW IMPLEMENTATION FOR DIR=-1
    GBTILT%nmul=2

    GBTILT%KIND=MADKIND2
    if(present(t)) then
       IF(T%NATURAL) THEN
          GBTILT%tilt=t%tilt(1)
       ELSE
          GBTILT%tilt=t%tilt(0)
       ENDIF
    endif

  END FUNCTION GBTILT


  FUNCTION  RECTTILT(NAME,L,ANGLE,E1,E2,T)
    implicit none
    type (EL_LIST) RECTTILT
    type (TILTING),OPTIONAL, INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,angle,E1,E2
    real(dp) L1,LM,ANG1,E11,E22

    L1=zero
    ANG1=zero
    IF(PRESENT(L)) LM=L
    IF(PRESENT(angle)) ANG1=angle
    E11=zero
    E22=zero

    IF(PRESENT(E1)) E11=E1
    IF(PRESENT(E2)) E22=E2

    IF(MADLENGTH.or.ang1==zero) THEN
       L1=LM
    ELSE
       L1=two*LM*SIN(ANG1/two)/ANG1
    ENDIF

    RECTTILT=0
    RECTTILT%B0=two*SIN(ANG1/two)/L1
    !    IF(ANG1==zero) THEN
    !       RECTTILT%L=L1
    !       RECTTILT%LD=L1
    !       RECTTILT%LC=L1
    !    ELSE
    IF(EXACT_MODEL) THEN
       if(verbose) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/,1x,a72))'
          w_p%c(1)= NAME
          w_p%c(2)= " READ AS TRUE RECTANGULAR BEND "
          call write_i
       endif
       if(ang1==zero) then
          RECTTILT%LD=L1
       else
          RECTTILT%LD=ANG1/RECTTILT%B0
       endif
       RECTTILT%L=L1
       RECTTILT%LC=L1
       RECTTILT%K(1)=RECTTILT%B0+RECTTILT%K(1)
       if(LIKEMAD) then
          RECTTILT%T1=ANG1/two+E11    !one
          RECTTILT%T2=ANG1/two+E22    !zero
       else
          RECTTILT%T1=ANG1/two+E11    !one
          RECTTILT%T2=ANG1/two+E22    !zero

          !             RECTTILT%T1=one   !wrong???
          !             RECTTILT%T2=zero
       endif
       RECTTILT%nmul=2
    ELSE
       RECTTILT%LC=L1
       IF(ANG1==ZERO) THEN
          RECTTILT%L=L1
          RECTTILT%LD=L1
       ELSE
          RECTTILT%L=ANG1/RECTTILT%B0
          RECTTILT%LD=ANG1/RECTTILT%B0
       ENDIF
       RECTTILT%T1=ANG1/two+E11 ; RECTTILT%T2=ANG1/two+E22;
       RECTTILT%K(1)=RECTTILT%B0+RECTTILT%K(1) ! NEW IMPLEMENTATION FOR DIR=-1
       RECTTILT%nmul=2   ! 0 before
    ENDIF
    !    ENDIF
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       RECTTILT%NAME=NAME(1:16)
    ELSE
       RECTTILT%NAME=NAME
    ENDIF

    RECTTILT%KIND=MADKIND2
    IF(present(t)) THEN
       IF(T%NATURAL) THEN
          RECTTILT%tilt=t%tilt(1)
       ELSE
          RECTTILT%tilt=t%tilt(0)
       ENDIF
    endif
  END FUNCTION RECTTILT

  FUNCTION  RECTTILT_MADX(NAME,T,LIST)
    implicit none
    type (EL_LIST) LIST
    type (EL_LIST) RECTTILT_MADX
    type (TILTING)T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) L1,LM,ANG1,E11,E22

    IF(EXACT_MODEL) LIKEMAD=.true.

    !    LIKEMAD=.TRUE.                       ! etienne check it  out
    L1=zero

    LM=LIST%L
    ANG1=LIST%B0


    E11=LIST%T1
    E22=LIST%T2

    IF(MADLENGTH.OR.ANG1==ZERO) THEN
       L1=LM
    ELSE
       L1=two*LM*SIN(ANG1/two)/ANG1
    ENDIF

    RECTTILT_MADX=LIST
    RECTTILT_MADX%B0=two*SIN(ANG1/two)/L1
    !    IF(ANG1==zero) THEN
    !       RECTTILT_MADX%L=L1
    !       RECTTILT_MADX%LD=L1
    !       RECTTILT_MADX%LC=L1
    !    ELSE
    IF(EXACT_MODEL) THEN
       if(verbose) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/,1x,a72))'
          w_p%c(1)= NAME
          w_p%c(2)= " READ AS TRUE RECTANGULAR BEND "
          call write_i
       endif
       if(ang1==zero) then
          RECTTILT_MADX%LD=L1
       else
          RECTTILT_MADX%LD=ANG1/RECTTILT_MADX%B0
       endif
       RECTTILT_MADX%L=L1
       RECTTILT_MADX%LC=L1
       RECTTILT_MADX%K(1)=RECTTILT_MADX%B0+RECTTILT_MADX%K(1)
       if(LIKEMAD) then
          RECTTILT_MADX%T1=ANG1/two+E11    !one
          RECTTILT_MADX%T2=ANG1/two+E22    !zero
       else
          RECTTILT_MADX%T1=ANG1/two+E11    !one
          RECTTILT_MADX%T2=ANG1/two+E22    !zero
          !             RECTTILT_MADX%T1=one   !wrong???
          !             RECTTILT_MADX%T2=zero
       endif
       RECTTILT_MADX%nmul=2
    ELSE
       RECTTILT_MADX%LC=L1
       IF(ANG1==ZERO) THEN
          RECTTILT_MADX%L=L1
          RECTTILT_MADX%LD=L1
       ELSE
          RECTTILT_MADX%L=ANG1/RECTTILT_MADX%B0
          RECTTILT_MADX%LD=ANG1/RECTTILT_MADX%B0
       ENDIF
       RECTTILT_MADX%T1=ANG1/two+E11 ; RECTTILT_MADX%T2=ANG1/two+E22;
       RECTTILT_MADX%K(1)=RECTTILT_MADX%B0+RECTTILT_MADX%K(1) ! NEW IMPLEMENTATION FOR DIR=-1
       RECTTILT_MADX%nmul=2   ! 0 before
    ENDIF
    !    ENDIF
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       RECTTILT_MADX%NAME=NAME(1:16)
    ELSE
       RECTTILT_MADX%NAME=NAME
    ENDIF

    RECTTILT_MADX%KIND=MADKIND2
    RECTTILT_MADX%tilt=t%tilt(0)
  END FUNCTION RECTTILT_MADX



  FUNCTION  rectaETILT(NAME,L,ANGLE,E1,E2,T)
    implicit none
    type (EL_LIST) rectaETILT
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,ANGLE,E1,E2
    type (TILTING), optional,INTENT(IN):: T
    real(dp) ANGE,SPE
    real(dp) LM1,ANG1,ANGI1,e11,e22
    integer tempkind

    CURVED_ELEMENT=.TRUE.

    E11=zero
    E22=zero
    tempkind=madkind2
    IF(PRESENT(ANGLE)) THEN
       if(ANGLE==zero) then
          madkind2=kind2
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=name
          WRITE(w_p%c(2),'(a12,a16,a23)') ' ANGLE=0 IN ', NAME,' CHANGED TO DRIFT-KICK '
          call write_i

       endif
    ELSE
       madkind2=kind2
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a12,a16,a23)') ' ANGLE=0 IN ', NAME,' CHANGED TO DRIFT-KICK '
       call write_i
    ENDIF
    IF((PRESENT(E1).AND.PRESENT(E2)).OR.(.NOT.PRESENT(E1).AND.(.NOT.PRESENT(E2))) ) THEN !1
       if(present(e1).and.present(e2)) THEN
          IF(EXACT_MODEL) LIKEMAD=.true.
          E11=E1
          E22=E2
       endif

       IF(present(t)) then
          rectaETILT=RECTTILT(NAME,L,ANGLE,E11,E22,T)
       else
          rectaETILT=RECTTILT(NAME,L,ANGLE,E11,E22)
       endif
       return

    ELSE  !  1

       LM1=zero
       ANG1=zero
       IF(PRESENT(L)) LM1=L
       IF(PRESENT(angle)) ANG1=angle

       IF(PRESENT(E1)) ANGI1=e1
       IF(PRESENT(E2)) ANGI1=ANG1-e2

       rectaETILT=0
       ANGE=ANG1-ANGI1
       SPE=ANG1/two-ANGI1

       IF(MADLENGTH) THEN
          rectaETILT%L=LM1
          rectaETILT%LC=rectaETILT%L/COS(SPE)
          rectaETILT%B0=two*SIN(ANG1/two)/rectaETILT%LC
          if(ang1/=zero) then
             rectaETILT%LD=ANG1/rectaETILT%B0
          else
             rectaETILT%LD=rectaETILT%LC
          endif
       ELSE
          rectaETILT%LD=LM1
          rectaETILT%B0=ANG1/rectaETILT%LD
          if(ang1/=zero) then
             rectaETILT%LC=two*SIN(ANG1/two)/rectaETILT%B0
          else
             rectaETILT%LC=rectaETILT%LD
          endif
          rectaETILT%L=rectaETILT%LC*COS(SPE)
       ENDIF


       IF(EXACT_MODEL) THEN
          if(verbose) then
             w_p=0
             w_p%nc=2
             w_p%fc='((1X,a72,/,1x,a72))'
             w_p%c(1)= NAME
             w_p%c(2)= " READ AS TRUE RECTANGULAR BEND "
             call write_i
          endif
          rectaETILT%K(1)=rectaETILT%B0+rectaETILT%K(1) ! NEW IMPLEMENTATION FOR DIR=-1
          rectaETILT%nmul=2
          !         rectaETILT%T1=ANGI1/(ANG1/two)
          rectaETILT%T1=ANGI1
          rectaETILT%T2=ange

          !         rectaETILT%T2=rectaETILT%LC*SIN(SPE)
       ELSE
          rectaETILT%K(1)=rectaETILT%B0+rectaETILT%K(1)
          rectaETILT%L=rectaETILT%LD
          rectaETILT%T1=ANGI1 ; rectaETILT%T2=ANGE;
          rectaETILT%nmul=2   ! 0 before
       ENDIF

       IF(LEN(NAME)>nlp) THEN
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=name
          WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
          call write_i
          rectaETILT%NAME=NAME(1:16)
       ELSE
          rectaETILT%NAME=NAME
       ENDIF

       rectaETILT%KIND=MADKIND2
       if(present(t)) then
          IF(T%NATURAL) THEN
             rectaETILT%tilt=t%tilt(1)
          ELSE
             rectaETILT%tilt=t%tilt(0)
          ENDIF
       endif

    ENDIF !1
    madkind2=TEMPKIND

  END FUNCTION rectaETILT

  FUNCTION  TRUERECTILT_MADX(NAME,T,LIST)
    implicit none
    type (EL_LIST) TRUERECTILT_MADX
    type (EL_LIST),optional, INTENT(IN)::list
    CHARACTER(*), INTENT(IN):: NAME
    type (TILTING),INTENT(IN):: T
    real(dp) ANGE,SPE
    real(dp) LM1,ANG1,ANGI1,e1,e2

    CURVED_ELEMENT=.TRUE.

    E1=LIST%T1
    E2=LIST%T2


    LM1=LIST%L
    ANG1=LIST%B0

    IF(E1/=zero) THEN
       ANGI1=e1
    ELSEIF(E1/=zero) THEN
       ANGI1=ANG1-e2
    ELSE
       WRITE(6,*) " ERROR IN  TRUERECTILT_MADX INPUT "
       STOP 222
    ENDIF
    TRUERECTILT_MADX=LIST
    ANGE=ANG1-ANGI1
    SPE=ANG1/two-ANGI1

    IF(MADLENGTH) THEN
       TRUERECTILT_MADX%L=LM1
       TRUERECTILT_MADX%LC=TRUERECTILT_MADX%L/COS(SPE)
       TRUERECTILT_MADX%B0=two*SIN(ANG1/two)/TRUERECTILT_MADX%LC
       TRUERECTILT_MADX%LD=ANG1/TRUERECTILT_MADX%B0
    ELSE
       TRUERECTILT_MADX%LD=LM1
       TRUERECTILT_MADX%B0=ANG1/TRUERECTILT_MADX%LD
       TRUERECTILT_MADX%LC=two*SIN(ANG1/two)/TRUERECTILT_MADX%B0
       TRUERECTILT_MADX%L=TRUERECTILT_MADX%LC*COS(SPE)
    ENDIF


    IF(EXACT_MODEL) THEN
       if(verbose) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/,1x,a72))'
          w_p%c(1)= NAME
          w_p%c(2)= " READ AS TRUE RECTANGULAR BEND "
          call write_i
       endif
       TRUERECTILT_MADX%K(1)=TRUERECTILT_MADX%B0+TRUERECTILT_MADX%K(1) ! NEW IMPLEMENTATION FOR DIR=-1
       TRUERECTILT_MADX%nmul=2
       !         rectaETILT%T1=ANGI1/(ANG1/two)
       TRUERECTILT_MADX%T1=ANGI1
       TRUERECTILT_MADX%T2=ange

       !         rectaETILT%T2=rectaETILT%LC*SIN(SPE)
    ELSE
       TRUERECTILT_MADX%K(1)=TRUERECTILT_MADX%B0+TRUERECTILT_MADX%K(1)
       TRUERECTILT_MADX%L=TRUERECTILT_MADX%LD
       TRUERECTILT_MADX%T1=ANGI1 ; TRUERECTILT_MADX%T2=ANGE;
       TRUERECTILT_MADX%nmul=2   ! 0 before
    ENDIF

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       TRUERECTILT_MADX%NAME=NAME(1:16)
    ELSE
       TRUERECTILT_MADX%NAME=NAME
    ENDIF

    TRUERECTILT_MADX%KIND=MADKIND2
    TRUERECTILT_MADX%tilt=t%tilt(0)


  END FUNCTION TRUERECTILT_MADX





  FUNCTION  drft(NAME,L,LIST)
    implicit none
    type (EL_LIST) drft
    CHARACTER(*), INTENT(IN):: NAME
    TYPE(EL_LIST) ,optional, INTENT(IN):: LIST
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L

    if(present(list)) then
       drft=list
       l1=list%L
    else
       drft=0
    endif

    drft%L=L1
    drft%LD=L1
    drft%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       drft%NAME=NAME(1:16)
    ELSE
       drft%NAME=NAME
    ENDIF
    drft%KIND=KIND1

  END FUNCTION drft

  FUNCTION  RCOLIT(NAME,L,T,LIST)
    implicit none
    integer ipause, mypause
    type (EL_LIST) RCOLIT
    type (EL_LIST),OPTIONAL,INTENT(IN):: LIST
    type (TILTING),OPTIONAL,INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L

    if(present(list)) then
       RCOLIT=list
       l1=list%L
       WRITE(6,*) " WHAT ABOUT WRITING THE CODE USING X AND Y"
       ipause=mypause(0)
    else
       RCOLIT=0
    endif

    RCOLIT%L=L1
    RCOLIT%LD=L1
    RCOLIT%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       RCOLIT%NAME=NAME(1:16)
    ELSE
       RCOLIT%NAME=NAME
    ENDIF
    RCOLIT%KIND=KIND18
    if(present(t)) then
       RCOLIT%tilt=t%tilt(0)
    endif

  END FUNCTION RCOLIT

  FUNCTION  ECOLIT(NAME,L,T,LIST)
    implicit none
    integer ipause, mypause
    type (EL_LIST) ECOLIT
    type (EL_LIST),OPTIONAL,INTENT(IN):: LIST
    type (TILTING),OPTIONAL,INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L

    if(present(list)) then
       ECOLIT=list
       l1=list%L
       WRITE(6,*) " WHAT ABOUT WRITING THE CODE USING X AND Y"
       ipause=mypause(0)

    else
       ECOLIT=0
    endif

    ECOLIT%L=L1
    ECOLIT%LD=L1
    ECOLIT%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       ECOLIT%NAME=NAME(1:16)
    ELSE
       ECOLIT%NAME=NAME
    ENDIF
    ECOLIT%KIND=KIND19
    if(present(t)) then
       ECOLIT%tilt=t%tilt(0)
    endif

  END FUNCTION ECOLIT

  FUNCTION  MONIT(NAME,L,T,LIST)
    implicit none
    type (EL_LIST) MONIT
    type (EL_LIST),OPTIONAL,INTENT(IN):: LIST
    type (TILTING),OPTIONAL,INTENT(IN):: T
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L

    if(present(list)) then
       MONIT=list
       l1=list%L
    else
       MONIT=0
    endif

    MONIT%L=L1
    MONIT%LD=L1
    MONIT%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       MONIT%NAME=NAME(1:16)
    ELSE
       MONIT%NAME=NAME
    ENDIF
    MONIT%KIND=KIND11
    if(present(t)) then
       MONIT%tilt=t%tilt(0)
    endif

  END FUNCTION MONIT

  FUNCTION  hMONIT(NAME,L)
    implicit none
    type (EL_LIST) hMONIT
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L


    hMONIT=0
    hMONIT%L=L1
    hMONIT%LD=L1
    hMONIT%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       hMONIT%NAME=NAME(1:16)
    ELSE
       hMONIT%NAME=NAME
    ENDIF
    hMONIT%KIND=KIND12

  END FUNCTION hMONIT

  FUNCTION  VMONIT(NAME,L)
    implicit none
    type (EL_LIST) VMONIT
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L


    VMONIT=0
    VMONIT%L=L1
    VMONIT%LD=L1
    VMONIT%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       VMONIT%NAME=NAME(1:16)
    ELSE
       VMONIT%NAME=NAME
    ENDIF
    VMONIT%KIND=KIND13

  END FUNCTION VMONIT

  FUNCTION  INSTRUMEN(NAME,L)
    implicit none
    type (EL_LIST) INSTRUMEN
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L
    real(dp)  L1
    L1=zero
    IF(PRESENT(L)) L1=L


    INSTRUMEN=0
    INSTRUMEN%L=L1
    INSTRUMEN%LD=L1
    INSTRUMEN%LC=L1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       INSTRUMEN%NAME=NAME(1:16)
    ELSE
       INSTRUMEN%NAME=NAME
    ENDIF
    INSTRUMEN%KIND=KIND14

  END FUNCTION INSTRUMEN

  FUNCTION  mark(NAME)
    implicit none
    type (EL_LIST) mark
    CHARACTER(*), INTENT(IN):: NAME

    mark=0
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       mark%NAME=NAME(1:16)
    ELSE
       mark%NAME=NAME
    ENDIF

    mark%KIND=KIND0

  END FUNCTION mark

  FUNCTION  CHANGEREF(NAME,ANG,T,PATCHG)
    implicit none
    type (EL_LIST) CHANGEREF
    CHARACTER(*), INTENT(IN):: NAME
    REAL(DP) ANG(3),T(3)
    INTEGER PATCHG

    CHANGEREF=0
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       CHANGEREF%NAME=NAME(1:16)
    ELSE
       CHANGEREF%NAME=NAME
    ENDIF

    CHANGEREF%KIND=KIND0
    CHANGEREF%ANG=ANG
    CHANGEREF%T=T
    CHANGEREF%PATCHG=PATCHG

  END FUNCTION CHANGEREF

  subroutine  guirder(f,cell)
    implicit none
    type (fibre) f
    type (layout),target :: cell

    f%MAG%G23=>CELL
    f%MAGP%G23=>CELL
    f%MAG%KIND=KIND23
    f%MAGP%KIND=KIND23
    f%MAG%p%nst=1
    f%MAGP%p%nst=1
    f%chart%f%ent=1
    f%chart=0
    CALL SURVEY_no_patch(f)


  END  subroutine guirder

  FUNCTION  RFCAVITYL(NAME,L,VOLT,LAG,HARMON,REV_FREQ,DELTAE,LIST)
    implicit none
    type (EL_LIST) RFCAVITYL
    TYPE(EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,VOLT,LAG,REV_FREQ,DELTAE
    INTEGER,optional, INTENT(IN):: HARMON
    real(dp)  L1,VOLT1,LAG1,FREQ01
    INTEGER  HARMON1
    L1=zero
    VOLT1=zero
    LAG1=zero
    FREQ01=zero
    HARMON1=1
    IF(PRESENT(L)) L1=L
    IF(PRESENT(VOLT)) THEN
       VOLT1=VOLT
       IF(PRESENT(DELTAE)) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          w_p%c(1)= "Use either volt or deltae"
          call write_e(100)
       ENDIF
    elseIF(PRESENT(DELTAE)) THEN
       volt1=DELTAE*p0c
    endif
    IF(PRESENT(LAG)) LAG1=LAG
    IF(PRESENT(HARMON)) HARMON1=HARMON
    IF(PRESENT(REV_FREQ)) FREQ01=REV_FREQ

    if(present(list)) then
       RFCAVITYL=list
       l1=list%L
       VOLT1=LIST%VOLT
       LAG1=LIST%LAG
       FREQ01=LIST%FREQ0
       HARMON1=LIST%HARMON
       if(LIST%delta_e/=zero) then
          if(volt1==zero) then
             volt1=LIST%DELTA_E*p0c    ! DELTA_E used for two purposes, but OK
          else
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,a72))'
             w_p%c(1)= "Use either volt or deltae"
             call write_e(101)
          endif
       endif
    else
       RFCAVITYL=0
    endif

    RFCAVITYL%L=L1
    RFCAVITYL%LD=L1
    RFCAVITYL%LC=L1
    RFCAVITYL%KIND=KIND4
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       RFCAVITYL%NAME=NAME(1:16)
    ELSE
       RFCAVITYL%NAME=NAME
    ENDIF
    RFCAVITYL%VOLT=VOLT1
    RFCAVITYL%LAG=LAG1
    RFCAVITYL%HARMON=HARMON1
    RFCAVITYL%FREQ0=FREQ01
    !   RFCAVITYL%P0C=P0C
    RFCAVITYL%DELTA_E=zero

  END FUNCTION RFCAVITYL

  FUNCTION  TWCAVITYL(NAME,L,VOLT,LAG,HARMON,REV_FREQ,DELTAE,LIST)
    implicit none
    type (EL_LIST) TWCAVITYL
    TYPE(EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,VOLT,LAG,REV_FREQ,DELTAE
    INTEGER,optional, INTENT(IN):: HARMON
    real(dp)  L1,VOLT1,LAG1,FREQ01
    INTEGER  HARMON1
    L1=zero
    VOLT1=zero
    LAG1=zero
    FREQ01=zero
    HARMON1=1
    IF(PRESENT(L)) L1=L
    IF(PRESENT(VOLT)) THEN
       VOLT1=VOLT
       IF(PRESENT(DELTAE)) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          w_p%c(1)= "Use either volt or deltae"
          call write_e(100)
       ENDIF
    elseIF(PRESENT(DELTAE)) THEN
       volt1=DELTAE*p0c
    endif
    IF(PRESENT(LAG)) LAG1=LAG
    IF(PRESENT(HARMON)) HARMON1=HARMON
    IF(PRESENT(REV_FREQ)) FREQ01=REV_FREQ

    if(present(list)) then
       TWCAVITYL=list
       l1=list%L
       VOLT1=LIST%VOLT
       LAG1=LIST%LAG
       FREQ01=LIST%FREQ0
       HARMON1=LIST%HARMON
       if(LIST%delta_e/=zero) then
          if(volt1==zero) then
             volt1=LIST%DELTA_E*p0c    ! DELTA_E used for two purposes, but OK
          else
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,a72))'
             w_p%c(1)= "Use either volt or deltae"
             call write_e(101)
          endif
       endif
    else
       TWCAVITYL=0
    endif
    IF(L1==ZERO) THEN
       WRITE(6,*) " TWCAVITY MUST HAVE A LENGTH "
       STOP 555
    ENDIF

    TWCAVITYL%L=L1
    TWCAVITYL%LD=L1
    TWCAVITYL%LC=L1
    TWCAVITYL%KIND=KIND21
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       TWCAVITYL%NAME=NAME(1:16)
    ELSE
       TWCAVITYL%NAME=NAME
    ENDIF
    TWCAVITYL%VOLT=VOLT1
    TWCAVITYL%LAG=LAG1
    TWCAVITYL%HARMON=HARMON1
    TWCAVITYL%FREQ0=FREQ01
    !   RFCAVITYL%P0C=P0C
    TWCAVITYL%DELTA_E=zero

  END FUNCTION TWCAVITYL



  FUNCTION  ELSESTILT(NAME,L,E,T,LIST)
    implicit none
    type (TILTING),optional, INTENT(IN):: T
    type (EL_LIST),optional, INTENT(IN):: LIST
    type (EL_LIST) ELSESTILT
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L,E
    real(dp) L1,K11

    L1=zero
    K11=zero
    IF(PRESENT(L)) L1=L
    IF(PRESENT(E)) K11=E

    if(present(list)) then
       ELSESTILT=list
       l1=list%L
       K11=LIST%VOLT
    else
       ELSESTILT=0
    endif
    ELSESTILT%L=L1
    ELSESTILT%LD=L1
    ELSESTILT%LC=L1
    ELSESTILT%VOLT=K11
    ELSESTILT%KIND=KIND15

    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          ELSESTILT%tilt=t%tilt(1)
       ELSE
          ELSESTILT%tilt=t%tilt(0)
       ENDIF
    ENDIF

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       ELSESTILT%NAME=NAME(1:16)
    ELSE
       ELSESTILT%NAME=NAME
    ENDIF

  END FUNCTION ELSESTILT



  FUNCTION  USER_1L(NAME,L,T,list)
    implicit none
    type (EL_LIST) USER_1L
    type (TILTING),optional, INTENT(IN):: T
    type (EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L

    if(present(list)) then
       USER_1L=list
       USER_1L%L=list%L
    elseif(present(L)) then
       USER_1L=0
       USER_1L%L=L
       USER_1L%LD=L
       USER_1L%LC=L
    else
       write(6,*) " Error neither L nor list is present in USER_1L"
       stop 900
    endif
    USER_1L%KIND=KINDUSER1
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       USER_1L%NAME=NAME(1:16)
    ELSE
       USER_1L%NAME=NAME
    ENDIF
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          USER_1L%tilt=t%tilt(1)
       ELSE
          USER_1L%tilt=t%tilt(0)
       ENDIF
    ENDIF
  END FUNCTION USER_1L


  FUNCTION  USER_2L(NAME,L,T,list)
    implicit none
    type (EL_LIST) USER_2L
    type (TILTING),optional, INTENT(IN):: T
    type (EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L

    if(present(list)) then
       USER_2L=list
       USER_2L%L=list%L
    elseif(present(L)) then
       USER_2L=0
       USER_2L%L=L
       USER_2L%LD=L
       USER_2L%LC=L
    else
       write(6,*) " Error neither L nor list is present in USER_2L"
       stop 900
    endif
    USER_2L%KIND=KINDUSER2
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       USER_2L%NAME=NAME(1:16)
    ELSE
       USER_2L%NAME=NAME
    ENDIF
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          USER_2L%tilt=t%tilt(1)
       ELSE
          USER_2L%tilt=t%tilt(0)
       ENDIF
    ENDIF
  END FUNCTION USER_2L

  FUNCTION  WIGGLERL(NAME,L,T,list)
    implicit none
    type (EL_LIST) WIGGLERL
    type (TILTING),optional, INTENT(IN):: T
    type (EL_LIST),optional, INTENT(IN):: LIST
    CHARACTER(*), INTENT(IN):: NAME
    real(dp) ,optional, INTENT(IN):: L

    if(present(list)) then
       WIGGLERL=list
       WIGGLERL%L=list%L
    elseif(present(L)) then
       WIGGLERL=0
       WIGGLERL%L=L
       WIGGLERL%LD=L
       WIGGLERL%LC=L
    else
       write(6,*) " Error neither L nor list is present in WIGGLERL"
       stop 900
    endif
    WIGGLERL%KIND=KINDWIGGLER
    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       WIGGLERL%NAME=NAME(1:16)
    ELSE
       WIGGLERL%NAME=NAME
    ENDIF
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          WIGGLERL%tilt=t%tilt(1)
       ELSE
          WIGGLERL%tilt=t%tilt(0)
       ENDIF
    ENDIF
  END FUNCTION WIGGLERL



  SUBROUTINE  EL_Q(s22,S1)
    !changed
    implicit none
    type (fibre),target,INTENT(inOUT)::s22
    type (EL_LIST),INTENT(IN)::S1
    INTEGER I,flip
    logical(lp) DONE,THICKKICKTEMP
    type(element),pointer :: s2
    type(elementp), pointer :: s2p
    type(fibre), pointer::el
    integer ntot,ntot_rad,ntot_REV,ntot_rad_REV

    nullify(el);
    THICKKICKTEMP=.FALSE.

    nullify(s2); nullify(s2p);
    IF(MADX_MAGNET_ONLY) THEN
       S22%MAG=-1;     !  FIBRE AND  MUST ALREADY EXIST
       S22%MAGP=-1;    !  POINTER MUST STAY ALLOCATED OTHERWISE ALL HELL BREAKS LOOSE
    ELSE
       nullify(s22%mag); nullify(s22%magp);
       allocate(s22%mag);allocate(s22%magp);
       nullify(s22%CHART);nullify(s22%PATCH);
       allocate(s22%CHART);allocate(s22%PATCH);
       nullify(s22%dir);allocate(s22%dir);
       !     nullify(s22%P0C);allocate(s22%P0C);
       !    nullify(s22%BETA0);allocate(s22%BETA0);
       nullify(s22%PARENT_LAYOUT);
       !       nullify(s22%PARENT_PATCH);
       !       nullify(s22%PARENT_CHART);nullify(s22%PARENT_MAG);
       NULLIFY(S22%I)
       if(use_info) then
          allocate(s22%i);
          call alloc(s22%i)
       endif
    ENDIF

    IF(.NOT.MADX) then
       nullify(s22%next);
       nullify(s22%previous);
    endif
    ! CALL ALLOCATE_FIBRE(S22)
    ! CALL ALLOCATE_DATA_FIBRE(S22)  !ONLY ALLOWED ON POINTERS
    IF(.NOT.MADX_MAGNET_ONLY) THEN
       s22%dir=FIBRE_DIR    ! ALL THAT SHIT ALREADY EXISTS
       !     s22%P0C=P0C
       !     s22%BETA0=BETA0
       s22%CHART=0
       s22%PATCH=0
    ENDIF
    ! New stuff
    !Powering the CHART frame in MAG only
    !
    !
    flip=1
    if(FIBRE_flip) flip=FIBRE_dir
    s2=>s22%mag;
    s2p=>s22%magp;

    DONE=.FALSE.

    DO I=NMAX,1,-1
       IF(S1%K(I)/=zero.or.S1%KS(I)/=zero) THEN
          if(I>=S1%NMUL) THEN
             S2 = I
             DONE=.TRUE.
          ENDIF
          GOTO 100
       ENDIF
    ENDDO
100 CONTINUE

    IF(.NOT.DONE) S2 = S1%NMUL

    S2%P%B0=S1%B0
    !    if(S2%P%B0/=zero) S2%P%bend_fringe=.true.
    IF(CURVED_ELEMENT) THEN
       S2%P%bend_fringe=.true.
       CURVED_ELEMENT=.FALSE.
    ENDIF
    S2%KIND=S1%KIND; S2%P%METHOD=S1%METHOD ;        S2%P%NST=S1%NST ;
    S2%NAME=S1%NAME        ;S2%VORNAME=S1%VORNAME ;
    S2%L =S1%L ;S2%P%LD=S1%LD;S2%P%LC=S1%LC;

    S2%PERMFRINGE=S1%PERMFRINGE
    S2%P%KILL_EXI_FRINGE=S1%KILL_EXI_FRINGE
    S2%P%KILL_ENT_FRINGE=S1%KILL_ENT_FRINGE
    !    S2%P%BEND_FRINGE=S1%BEND_FRINGE    ! SET ON THE BASIS OF B0

    DO I=1,S2%P%NMUL
       S2 %BN(I)=flip*S1%K(I) ;S2 %AN(I)=flip*S1%KS(I)
    ENDDO
    S2%p%exact=EXACT_MODEL
    !    IF(S2%p%EXACT) THEN
    S2%P%EDGE(1)=(S1%T1)
    S2%P%EDGE(2)=(S1%T2)
    !     ENDIF
    ! S2%B0=S1%B0
    s2%P%tiltd=S1%tilt
    if(s1%kind==kind4) then
       ALLOCATE(S2%VOLT,S2%FREQ,S2%PHAS,S2%DELTA_E,S2%THIN)
       S2%volt=flip*S1%volt
       S2%freq=S1%freq0*S1%harmon
       S2%phas=-S1%lag
       !      S2%p0c=S1%p0c
       !frs
       S2%DELTA_E=S1%DELTA_E
       S2%THIN=.FALSE.
       IF(S2%L==zero) then
          S2%THIN=.TRUE.

       else
          S2%volt=S2%volt/S2%L
       endif
    endif

    if(s1%kind==kind21) then
       ALLOCATE(S2%VOLT,S2%FREQ,S2%PHAS,S2%LAG,S2%DELTA_E,S2%THIN)
       S2%volt=flip*S1%volt
       S2%freq=S1%freq0*S1%harmon
       S2%phas=-S1%lag
       !      S2%p0c=S1%p0c
       !frs
       S2%DELTA_E=S1%DELTA_E
       S2%THIN=.FALSE.
       !skowron 14.03.06
       S2%lag=s1%lag
       IF(S2%L==zero) then
          S2%THIN=.TRUE.
       else
          S2%volt=S2%volt/S2%L
       endif

    endif

    if(s1%kind==kind15) then
       ALLOCATE(S2%VOLT)
       S2%volt=S1%volt
       ALLOCATE(S2%phas)
       S2%phas=S1%lag
    endif

    if(s1%kind==kind5.or.s1%kind==kind17) then
       ALLOCATE(S2%B_SOL);
       S2%B_SOL=S1%BSOL
    endif


    CALL CONTEXT( S2%NAME )
    S2%P%BETA0=BETA0
    S2%P%gamma0I=gamma0I
    S2%P%gambet=gambet
    S2%P%p0c=p0c


    if(S2%KIND==KIND2.AND.EXACT_MODEL) then
       S2%KIND=KIND16
    endif

    if((S2%KIND==KIND6.or.S2%KIND==KIND7).AND.EXACT_MODEL.AND.S2%P%B0/=zero) then
       S2%KIND=KIND16
       THICKKICKTEMP=.TRUE.
    endif

    ntot=0; ntot_rad=0; ntot_REV=0 ; ntot_rad_REV=0;
    if(S2%KIND==KIND22) then
       IF(ASSOCIATED(mad_tree%CC)) ntot=mad_tree%n
       IF(ASSOCIATED(mad_tree_rad%CC)) ntot_rad=mad_tree_rad%n
       IF(ASSOCIATED(mad_tree_REV%CC)) ntot_REV=mad_tree_REV%n
       IF(ASSOCIATED(mad_tree_RAD_REV%CC)) ntot_rad_REV=mad_tree_RAD_REV%n
    endif

    !    CALL SETFAMILY(S2,ntot,ntot_rad,ntot_REV,ntot_rad_REV,6)
    if(s2%kind/=kindpa) then
       CALL SETFAMILY(S2,NTOT=ntot,ntot_rad=ntot_rad,NTOT_REV=ntot_REV,ntot_rad_REV=ntot_rad_REV,ND2=6)
    else
       CALL SETFAMILY(S2,t=T_E,T_ax=T_ax,T_ay=T_ay)
       S2%P%METHOD=4
       deallocate(T_E,t_ax,t_ay)
    endif
    IF(S2%KIND==KIND22) THEN
       S2%M22%DELTAMAP=MAD_TREE_DELTAMAP
       if(associated(s2%m22%t)) call copy_tree(mad_tree,s2%m22%t)
       call KILL(mad_tree)
       if(associated(s2%m22%t_rad)) call copy_tree(mad_tree_rad,s2%m22%t_rad)
       call KILL(mad_tree_rad)
       if(associated(s2%m22%t_rev)) call copy_tree(mad_tree_rev,s2%m22%t_rev)
       call KILL(mad_tree_rev)
       if(associated(s2%m22%t_rad_rev)) call copy_tree(mad_tree_rad_rev,s2%m22%t_rad_rev)
       call KILL(mad_tree_rad_rev)
    ENDIF

    IF(S2%KIND==KIND4) THEN
       S2%C4%N_BESSEL=S1%N_BESSEL
    ENDIF
    IF(S2%KIND==KIND21) THEN
       s2%CAV21%DPHAS=s1%DPHAS
       s2%CAV21%dvds=s1%dvds
       s2%CAV21%PSI=s1%PSI
    ENDIF

    if(LIKEMAD) then
       if(S2%KIND/=KIND16) then
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          w_p%c(1)= " Likemad is true and element is not STREX "
          call write_e(kind16)
       endif
       s2%k16%likemad=LIKEMAD
       S2%KIND=KIND20
       LIKEMAD=.false.
    endif


    if(S2%KIND==KIND10) then
       S2%TP10%DRIFTKICK=DRIFT_KICK
       IF(madkind2==kind6.or.madkind2==kind7)   S2%TP10%DRIFTKICK=.FALSE.   ! 2002.11.04
       IF(S2%p%b0==zero)   then
          S2%TP10%DRIFTKICK=.true.
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1x,a72))'
          w_p%c(1)=S2%name
          WRITE(w_p%c(2),'(a12,a16,a23)') ' ANGLE=0 IN ', S2%name,' CHANGED TO DRIFT-KICK '
          call write_i
       endif
    endif

    if(S2%KIND==KIND16.OR.S2%KIND==KIND20) then
       IF(S2%P%B0/=0.AND.(.NOT.DRIFT_KICK)) THEN
          S2%K16%DRIFTKICK=.FALSE.
       ELSE
          S2%K16%DRIFTKICK=.TRUE.
       ENDIF
       IF(THICKKICKTEMP)   S2%K16%DRIFTKICK=.FALSE.
    ENDIF

    IF(S2%KIND==KIND18) S2%RCOL18%A%KIND=2
    IF(S2%KIND==KIND19) S2%ECOL19%A%KIND=1

    IF(MADX) then
       s2%fint=s1%FINT
       s2%hgap=s1%hgap
       s2%h1=s1%h1
       s2%h2=s1%h2
       IF(S2%KIND==KIND3) THEN
          s2%K3%thin_h_foc=s1%thin_h_foc
          s2%K3%thin_v_foc=s1%thin_v_foc
          s2%K3%thin_h_angle=s1%thin_h_angle
          s2%K3%thin_v_angle=s1%thin_v_angle
       ENDIF
       if(s1%APERTURE_KIND/=0) then
          call alloc(s2%p%aperture)
          s2%p%aperture%kind = -s1%APERTURE_KIND
          if(s1%aperture_on) s2%p%aperture%kind =-s2%p%aperture%kind
          s2%p%aperture%r    = s1%APERTURE_R
          s2%p%aperture%x    = s1%APERTURE_X
          s2%p%aperture%y    = s1%APERTURE_y
       endif
    endif

    s2p=0
    call copy(s2,s2p)

    ! end of machida stuff here
    ! Default survey stuff here
    !         s22%CHART%A_XY=s2%P%tilTd      ! THAT SHIT SHOULD NOT BE CHANGED NORMALLY
    !         s22%CHART%L=s2%P%LC
    !         s22%CHART%ALPHA=s2%P%LD*s2%P%B0
    IF(.NOT.MADX_MAGNET_ONLY) THEN      ! THIS SHOULD BE DONE MANUALLY IF NECESSARY
       if(associated(s22%chart%f)) then
          s22%chart%f%ent=1
          !           s22%chart=1
          s22%chart=2
          CALL SURVEY_no_patch(S22)
       endif
    ENDIF


    if(s1%patchg/=0) then
       if(s1%patchg==4) then   ! zgoubi order using two patches
          s22%PATCH%B_ANG=s1%ang    !
          s22%PATCH%A_D=s1%t
          s22%PATCH%patch=3
       elseif(s1%patchg==2) then
          s22%PATCH%B_ANG=s1%ang    !
          s22%PATCH%B_D=s1%t
          s22%PATCH%patch=2
       elseif(s1%patchg==1) then
          s22%PATCH%A_ANG=s1%angi    !
          s22%PATCH%A_D=s1%ti
          s22%PATCH%patch=1
       elseif(s1%patchg==3) then
          s22%PATCH%A_ANG=s1%angi    !
          s22%PATCH%A_D=s1%ti
          s22%PATCH%B_ANG=s1%ang    !
          s22%PATCH%B_D=s1%t
          s22%PATCH%patch=3
       endif
    endif



    IF(.NOT.MADX) THEN
       el=>s22
       !    call APPEND_mad_like(mad_list,s22)
       call APPEND_mad_like(mad_list,el)
    ENDIF
    madkick=.false.

    if(s22%mag%kind==kind3) then
       s22%mag%p%nst=1
       s22%magp%p%nst=1
    endif
    if(s22%mag%L==ZERO) then
       s22%mag%p%nst=1
       s22%magp%p%nst=1
    endif


  END SUBROUTINE EL_Q





  SUBROUTINE  clean_up
    implicit none

    call kill(mad_list)

    mad_list_killed=.true.
  end SUBROUTINE  clean_up

  subroutine set_pointers
    implicit none
    call set_da_pointers
    c_%NP_pol => NP_pol
    c_%NPARA => NPARA
    c_%ROOT_CHECK => ROOT_CHECK
    c_%CHECK_STABLE => CHECK_STABLE
    c_%CHECK_MADX_APERTURE => CHECK_MADX_APERTURE
    c_%ROOT_CHECK => ROOT_CHECK
    c_%APERTURE_FLAG => APERTURE_FLAG
    c_%absolute_aperture => absolute_aperture
    c_%check_x_min => check_x_min
    c_%check_x_max => check_x_max
    c_%check_y_min => check_y_min
    c_%check_y_max => check_y_max
    c_%hyperbolic_aperture => hyperbolic_aperture
    c_%WATCH_USER => WATCH_USER
    c_%ALWAYS_EXACTMIS=> ALWAYS_EXACTMIS

    c_%x_prime => x_prime
    c_%NDPT_OTHER => NDPT_OTHER

    c_%CAVITY_TOTALPATH => CAVITY_TOTALPATH


    c_%other_program => other_program
    c_%NEW_METHOD => NEW_METHOD
    c_%MADTHICK => MADKIND2
    c_%MADTHIN_NORMAL => MADKIND3N
    c_%MADTHIN_SKEW => MADKIND3S
    c_%NSTD => NSTD
    c_%METD => METD
    c_%MADLENGTH => MADLENGTH
    c_%MAD => MAD
    c_%EXACT_MODEL => EXACT_MODEL
    c_%ALWAYS_EXACTMIS => ALWAYS_EXACTMIS
    c_%HIGHEST_FRINGE => HIGHEST_FRINGE
    c_%FIBRE_DIR => FIBRE_DIR
    c_%FIBRE_flip => FIBRE_flip
    c_%SECTOR_NMUL => SECTOR_NMUL
    c_%SECTOR_NMUL_MAX => SECTOR_NMUL_MAX
    c_%electron => electron
    c_%massfactor => muon
    c_%stoch_in_rec => stoch_in_rec
    c_%FEED_P0C => FEED_P0C
    c_%ALWAYS_EXACT_PATCHING => ALWAYS_EXACT_PATCHING
    c_%OLD_IMPLEMENTATION_OF_SIXTRACK => OLD_IMPLEMENTATION_OF_SIXTRACK
    c_%wedge_coeff => wedge_coeff
    c_%MAD8_WEDGE => MAD8_WEDGE
    c_%phase0 => phase0

  end subroutine set_pointers

  SUBROUTINE  Set_mad(Energy,kinetic,p0c,BRHO,BETa,noisy,method,step)
    implicit none
    real(dp) ,optional, INTENT(IN)::Energy,kinetic,BRHO,BETa,p0c
    integer, optional, INTENT(IN)::method,step
    logical(lp), optional, INTENT(IN)::noisy

    real(dp) Energy1,kinetic1,BRHO1,BETa1,p0c1
    logical(lp) verb
    integer met,ns
    logical(lp) all

    IF(MAD8_WEDGE) THEN
       WEDGE_COEFF(1)=ONE+ONE/FOUR
       WEDGE_COEFF(2)=TWO-HALF
    ELSE
       WEDGE_COEFF(1)=ONE
       WEDGE_COEFF(2)=ONE
    ENDIF

    call set_pointers

    CALL NULL_TREE(mad_tree)
    CALL NULL_TREE(mad_tree_rad)
    CALL NULL_TREE(mad_tree_REV)
    CALL NULL_TREE(mad_tree_rad_REV)


    ns=nstd
    met=METD
    verb=verbose
    Energy1=zero
    kinetic1=zero
    p0c1=zero
    BRHO1=zero
    BETa1=zero
    all=.true.
    if(present(Energy)) then
       Energy1=-Energy
    else
       all=.false.
    endif
    if(present(kinetic)) then
       kinetic1=-kinetic
    else
       all=.false.
    endif
    if(present(p0c)) then
       p0c1=-p0c
    else
       all=.false.
    endif
    if(present(BRHO)) then
       BRHO1=-BRHO
    else
       all=.false.
    endif
    if(present(BETa)) then
       BETa1=-BETa
    else
       all=.false.
    endif
    if(present(noisy)) then
       verb=noisy
    else
       all=.false.
    endif
    if(present(method)) then
       met=method
    else
       all=.false.
    endif
    if(present(step)) then
       ns=step
    else
       all=.false.
    endif
    if(all) then
       Energy1=-Energy1
       p0c1=-p0c1
       BRHO1=-BRHO1
       kinetic1=-kinetic1
       BETa1=-BETa1
    endif
    call Set_mad_v(Energy1,kinetic1,p0c1,BRHO1,BETa1,verb,met,ns)

  end SUBROUTINE  Set_mad

  SUBROUTINE  Set_madx(Energy,kinetic,p0c,BRHO,BETa,noisy,method,step)
    implicit none
    real(dp) ,optional, INTENT(IN)::Energy,kinetic,BRHO,BETa,p0c
    integer, optional, INTENT(IN)::method,step
    logical(lp), optional, INTENT(IN)::noisy

    real(dp) Energy1,kinetic1,BRHO1,BETa1,p0c1
    logical(lp) verb
    integer met,ns
    logical(lp) all

    call set_pointers


    ns=nstd
    met=METD
    verb=verbose
    Energy1=zero
    kinetic1=zero
    p0c1=zero
    BRHO1=zero
    BETa1=zero
    all=.true.
    if(present(Energy)) then
       Energy1=-Energy
    else
       all=.false.
    endif
    if(present(kinetic)) then
       kinetic1=-kinetic
    else
       all=.false.
    endif
    if(present(p0c)) then
       p0c1=-p0c
    else
       all=.false.
    endif
    if(present(BRHO)) then
       BRHO1=-BRHO
    else
       all=.false.
    endif
    if(present(BETa)) then
       BETa1=-BETa
    else
       all=.false.
    endif
    if(present(noisy)) then
       verb=noisy
    else
       all=.false.
    endif
    if(present(method)) then
       met=method
    else
       all=.false.
    endif
    if(present(step)) then
       ns=step
    else
       all=.false.
    endif
    if(all) then
       Energy1=-Energy1
       p0c1=-p0c1
       BRHO1=-BRHO1
       kinetic1=-kinetic1
       BETa1=-BETa1
    endif
    madx=.true.
    call Set_mad_v(Energy1,kinetic1,p0c1,BRHO1,BETa1,verb,met,ns)
    madx=.false.
  end SUBROUTINE  Set_madx





  SUBROUTINE GET_ENERGY(ENE,KIN,BRHOin,BET,P0CC)
    implicit none
    real(dp) ,INTENT(INOUT)::ENE,kin,BRHOin,BET,P0CC
    ENE=ENERGY
    KIN=KINETIC
    BRHOIN=BRHO
    BET=BETA0
    P0CC=P0C

  end SUBROUTINE  GET_ENERGY

  SUBROUTINE  GET_GAM(GAMI,GAMB)
    implicit none
    real(dp) ,INTENT(INOUT)::GAMI,GAMB
    GAMI=gamma0I
    GAMB=gambet

  end SUBROUTINE  GET_GAM

  SUBROUTINE  GET_ONE(MASS,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)
    implicit none
    real(dp) ,optional,INTENT(OUT)::ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet,MASS
    real(dp)  ENE,kin,BRHOin,BET,P0CC,GAMI,GAMB

    call GET_ENERGY(ENE,KIN,BRHOin,BET,P0CC)
    CALL GET_GAM(GAMI,GAMB)

    if(present(ENERGY))  ENERGY=ENE
    if(present(KINETIC)) KINETIC=kin
    if(present(BRHO)) BRHO=BRHOin
    if(present(BETA0)) BETA0=BET
    if(present(P0C)) P0C=P0CC
    if(present(gamma0I)) gamma0I=GAMI
    if(present(gambet)) gambet=GAMB
    if(present(MASS)) MASS=mc2

  end SUBROUTINE  GET_ONE

  SUBROUTINE  Set_mad_v(ENE,KIN,p0c1,BRHOin,BET,verb,met,ns)
    implicit none
    real(dp) ,INTENT(IN)::ENE,BRHOin,BET,p0c1
    real(dp) XMC2,cl,CU,ERG,beta0i,GAMMA,GAMMA2,CON ,KIN
    logical(lp) PROTON,verb
    integer met,ns


    METD=met
    nstd=ns
    if(mad_list_killed.and.(.not.madx)) then
       call set_up(mad_list)
       mad_list_killed=.false.
    endif
    setmad = .true.
    verbose=verb
    !total_EPS=c_1d_10

    ENERGY=ENE
    KINETIC=KIN
    beta0=BET
    brho=BRHOin
    p0c=p0c1

    PROTON=.NOT.ELECTRON
    cl=(clight/c_1d8)
    CU=c_55/c_24/SQRT(three)
    w_p=0
    w_p%nc=8
    w_p%fc='(7((1X,A72,/)),1X,A72)'
    if(electron) then
       XMC2=muon*pmae
       w_p%c(1)=" This is an electron "
    elseif(proton) then
       XMC2=pmap
       w_p%c(2)=" This is a proton! "
    endif
    if(energy<0) then
       energy=-energy
       erg=ENERGY
       p0c=SQRT(erg**2-xmc2**2)
    endif
    if(KINETIC<0) then
       KINETIC=-KINETIC
       erg=KINETIC+xmc2
       p0c=SQRT(erg**2-xmc2**2)
    endif
    if(brho<0) then
       brho=-brho
       p0c=SQRT(BRHO**2*(cl/ten)**2)
    endif
    if(beta0<0) then
       beta0=-beta0
       p0c=(one-beta0**2)
       if(p0c<=zero) then
          w_p=0
          w_p%nc=2
          w_p%fc='(((1X,A72,/)),1X,A72)'
          write(w_p%c(1),'(a9,1x,g20.14)') " Beta0 = ",beta0
          w_p%c(2) ="Beta0 is too close to 1 "
          call write_e(-567)
       endif
       p0c=xmc2*beta0/SQRT(p0c)
    endif
    if(p0c<0) p0c=-p0c
    erg=SQRT(p0c**2+XMC2**2)
    ENERGY=ERG
    KINETIC=ERG-xmc2
    beta0=SQRT(KINETIC**2+two*KINETIC*XMC2)/erg
    beta0i=one/beta0
    GAMMA=erg/XMC2
    write(W_P%C(2),'(A16,G20.14)') ' Kinetic Energy ',kinetic
    write(W_P%C(3),'(A7,G20.14)') ' gamma ',gamma
    write(W_P%C(4),'(A7,G20.14)')' beta0 ',BETa0
    CON=three*CU*CGAM*HBC/two*TWOPII/XMC2**3
    CRAD=CGAM*ERG**3*TWOPII
    CFLUC=CON*ERG**5
    GAMMA2=erg**2/XMC2**2
    BRHO=SQRT(ERG**2-XMC2**2)*ten/cl
    write(W_P%C(5),'(A7,G20.14)') ' p0c = ',p0c
    write(W_P%C(6),'(A9,G20.14)')' GAMMA = ',SQRT(GAMMA2)
    write(W_P%C(7),'(A8,G20.14)')' BRHO = ',brho
    write(W_P%C(8),'(A15,G20.14,1X,G20.14)')"CRAD AND CFLUC ", crad ,CFLUC
    CALL WRITE_I
    !END OF SET RADIATION STUFF  AND TIME OF FLIGHT SUFF

    gamma0I=XMC2*BETA0/P0C
    GAMBET=(XMC2/P0C)**2
    MC2=XMC2
  END SUBROUTINE Set_mad_v

  FUNCTION  Taylor_maptilt(NAME,file,file_rev,T)
    implicit none
    type (EL_LIST) Taylor_maptilt
    CHARACTER(*), INTENT(IN):: NAME
    CHARACTER(*),optional, INTENT(IN):: file,file_rev
    type (TILTING),optional, INTENT(IN):: T
    integer mf,no,n_map,nst
    real(dp) ld,ang
    type(damap) m,mr,id,id2
    INTEGER I,ndpt,time,timefac,K
    type(taylor) beta,gamma
    INTEGER  JS(6)
    Taylor_maptilt=0
    JS=0
    Taylor_maptilt%NST=1
    IF(PRESENT(FILE)) THEN
       call kanalnummer(mf)
       open(unit=mf,file=file)
       read(mf,*) n_map,no,ang,ld,MAD_TREE_DELTAMAP,nst   ! number of maps (1,2), no, ld=design length
       Taylor_maptilt%NST=nst
       IF(MAD_TREE_DELTAMAP) THEN
          read(mf,*) ndpt,time,timefac                  ! npdt, time=0,1 (1 uses time), timefac = +/- 1)
       ENDIF
       if(symplectic_order>0) no=symplectic_order
       if(symplectify) then
          if(symplectic_order>no) no=symplectic_order
       endif
       call init(NO,3,0,0)
       call alloc(m,mr,id,id2); call alloc(beta,gamma);
       IF(MAD_TREE_DELTAMAP) THEN   !
          !  ndpt=6 is changed
          id=1
          if(ndpt==6) then
             id%v(5)=(one.mono.'000001')
             id%v(6)=-(one.mono.'00001')
          endif
          id%v(5)=id%v(5)*timefac
          id%v(6)=id%v(6)*timefac



          ! $$$$$$$$$$$$$$$$$
          call  dainput_SPECIAL6(m,mf,time); m=id**(-1).o.m.o.id ;

          if(time==2) then   ! make it like time==1   COSY INFINITY MAP
             id2=1
             gamma=sqrt(one+ (ONE+TWO*id2%v(5)/beta0+id2%v(5)**2 )/gambet)
             id2%v(5)=gamma0I*(gamma-one)/(one-gamma0I)-one
             id2%v(5)=id2%v(5)-(id2%v(5).sub.0)  ! almost useless since -one above does the job
             id2%v(6)=beta0*id2%v(6)*one/(one+gamma0i)   ! perhaps mistake in cosy manual


             m=id2**(-1).o.m.o.id2 ;
             DO K=0,NO
                JS(5)=K
                M%V(6)=M%V(6)-((M%V(6).SUB.JS).MONO.JS)
             ENDDO
             id2%v(5)=(one.mono.'00001')
             BETA=SQRT(ONE+TWO*id2%v(5)/BETA0+id2%v(5)**2  )/(ONE/BETA0 +id2%v(5))
             m%v(6)=m%v(6)+ LD*(ONE/BETA-ONE/BETA0)
          endif


          if(time==1.or.time==2) then
             id2=1
             id2%v(5)=(TWO*id2%v(5)+id2%v(5)**2)/(SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )+ONE/BETA0)
             m=id2**(-1).o.m.o.id2 ;
             id2%v(5)=one.mono.'00001'
             BETA=(ONE+id2%v(5) )/SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )
             m%v(6)=m%v(6)-((m%v(6).sub.'000001').mono.'000001')
             m%v(6)=beta*(m%v(6)- LD*(ONE/BETA-ONE/BETA0) )
             m%v(6)=m%v(6)+(one.mono.'000001')
          endif

          if(symplectify)call symplectic(m,symplectic_eps,nst)

          call SET_TREE(mad_tree,M)

          if(n_map>=2) then
             call  dainput_SPECIAL6(mr,mf,time); mr=id**(-1).o.mr.o.id ;

             if(time==2) then   ! make it like time==1   COSY INFINITY MAP
                id2=1
                gamma=sqrt(one+ (ONE+TWO*id2%v(5)/beta0+id2%v(5)**2 )/gambet)
                id2%v(5)=gamma0I*(gamma-one)/(one-gamma0I)-one
                id2%v(5)=id2%v(5)-(id2%v(5).sub.0)  ! almost useless since -one above does the job
                id2%v(6)=beta0*id2%v(6)*one/(one+gamma0i)   ! perhaps mistake in cosy manual


                m=id2**(-1).o.m.o.id2 ;
                DO K=0,NO
                   JS(5)=K
                   MR%V(6)=MR%V(6)-((MR%V(6).SUB.JS).MONO.JS)
                ENDDO
                id2%v(5)=(one.mono.'00001')
                BETA=SQRT(ONE+TWO*id2%v(5)/BETA0+id2%v(5)**2  )/(ONE/BETA0 +id2%v(5))
                MR%v(6)=MR%v(6)+ LD*(ONE/BETA-ONE/BETA0)
             endif

             if(time==1) then
                id2=1
                id2%v(5)=(TWO*id2%v(5)+id2%v(5)**2)/(SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )+ONE/BETA0)
                mr=id2**(-1).o.mr.o.id2 ;
                id2%v(5)=one.mono.'00001'
                BETA=(ONE+id2%v(5) )/SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )
                mr%v(6)=mr%v(6)-((mr%v(6).sub.'000001').mono.'000001')
                mr%v(6)=beta*(mr%v(6)- LD*(ONE/BETA-ONE/BETA0) )
                mr%v(6)=mr%v(6)+(one.mono.'000001')
             endif

             call SET_TREE(mad_tree_rad,Mr)
          else
             call SET_TREE(mad_tree_rad,M)
          endif

       ELSE  ! MAD_TREE_DELTAMAP
          call  dainput_SPECIAL6(m,mf,time);

          if(symplectify)call symplectic(m,symplectic_eps,nst)

          call SET_TREE(mad_tree,M)
          if(n_map==2) then
             call  dainput_SPECIAL6(mr,mf,time);
             call SET_TREE(mad_tree_rad,Mr)
          else
             call SET_TREE(mad_tree_rad,M)
          endif
       ENDIF ! MAD_TREE_DELTAMAP
       close(mf)
       if(symplectic_print) then
          call kanalnummer(mf)
          open(unit=mf,file=file)
          write(mf,*) n_map,no,ang,ld,MAD_TREE_DELTAMAP,nst   ! number of maps (1,2), no, ld=design length
          Taylor_maptilt%NST=nst
          IF(MAD_TREE_DELTAMAP) THEN
             read(mf,*) ndpt,time,timefac                  ! npdt, time=0,1 (1 uses time), timefac = +/- 1)
          ENDIF
          call daprint(m,mf)

          close(mf)
       endif


       call kill(m,mr,id,id2); call kill(beta,gamma);


       MAD_TREE_LD=LD   ! put here for the logic of PRESENT(FILE_REV)==false
       MAD_TREE_ANGLE=ANG

    ELSE  ! PRESENT FILE
       LD=MAD_TREE_LD
       ANG=MAD_TREE_ANGLE
    ENDIF

    IF(PRESENT(FILE_REV)) THEN
       call kanalnummer(mf)
       open(unit=mf,file=FILE_REV)
       read(mf,*) n_map,no,ang,ld,MAD_TREE_DELTAMAP,nst   ! number of maps (1,2), no, ld=design length
       Taylor_maptilt%NST=nst
       IF(MAD_TREE_DELTAMAP) THEN
          read(mf,*) ndpt,time,timefac                  ! npdt, time=0,1 (1 uses time), timefac = +/- 1)
       ENDIF

       if(symplectic_order>0) no=symplectic_order
       if(symplectify) then
          if(symplectic_order>no) no=symplectic_order
       endif

       call init(NO,3,0,0)
       call alloc(m,mr,id,id2); call alloc(beta,gamma);
       IF(MAD_TREE_DELTAMAP) THEN   !
          !  ndpt=6 is changed
          id=1
          if(ndpt==6) then
             id%v(5)=(one.mono.'000001')
             id%v(6)=-(one.mono.'00001')
          endif
          id%v(5)=id%v(5)*timefac
          id%v(6)=id%v(6)*timefac


          ! $$$$$$$$$$$$$$$$$

          call  dainput_SPECIAL6(m,mf,time); m=id**(-1).o.m.o.id ;

          if(time==2) then   ! make it like time==1   COSY INFINITY MAP
             id2=1
             gamma=sqrt(one+ (ONE+TWO*id2%v(5)/beta0+id2%v(5)**2 )/gambet)
             id2%v(5)=gamma0I*(gamma-one)/(one-gamma0I)-one
             id2%v(5)=id2%v(5)-(id2%v(5).sub.0)  ! almost useless since -one above does the job
             id2%v(6)=beta0*id2%v(6)*one/(one+gamma0i)   ! perhaps mistake in cosy manual


             m=id2**(-1).o.m.o.id2 ;
             DO K=0,NO
                JS(5)=K
                M%V(6)=M%V(6)-((M%V(6).SUB.JS).MONO.JS)
             ENDDO
             id2%v(5)=(one.mono.'00001')
             BETA=SQRT(ONE+TWO*id2%v(5)/BETA0+id2%v(5)**2  )/(ONE/BETA0 +id2%v(5))
             m%v(6)=m%v(6)+ LD*(ONE/BETA-ONE/BETA0)
          endif

          if(time==1) then
             id2=1
             id2%v(5)=(TWO*id2%v(5)+id2%v(5)**2)/(SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )+ONE/BETA0)
             m=id2**(-1).o.m.o.id2 ;
             id2%v(5)=one.mono.'00001'
             BETA=(ONE+id2%v(5) )/SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )
             m%v(6)=m%v(6)-((m%v(6).sub.'000001').mono.'000001')
             m%v(6)=beta*(m%v(6)- LD*(ONE/BETA-ONE/BETA0) )
             m%v(6)=m%v(6)+(one.mono.'000001')
          endif


          if(symplectify)call symplectic(m,symplectic_eps,nst)

          call SET_TREE(mad_tree_rev,M)

          if(n_map>=2) then
             call  dainput_SPECIAL6(mr,mf,time); mr=id**(-1).o.mr.o.id ;

             if(time==2) then   ! make it like time==1   COSY INFINITY MAP
                id2=1
                gamma=sqrt(one+ (ONE+TWO*id2%v(5)/beta0+id2%v(5)**2 )/gambet)
                id2%v(5)=gamma0I*(gamma-one)/(one-gamma0I)-one
                id2%v(5)=id2%v(5)-(id2%v(5).sub.0)  ! almost useless since -one above does the job
                id2%v(6)=beta0*id2%v(6)*one/(one+gamma0i)   ! perhaps mistake in cosy manual


                m=id2**(-1).o.m.o.id2 ;
                DO K=0,NO
                   JS(5)=K
                   MR%V(6)=MR%V(6)-((MR%V(6).SUB.JS).MONO.JS)
                ENDDO
                id2%v(5)=(one.mono.'00001')
                BETA=SQRT(ONE+TWO*id2%v(5)/BETA0+id2%v(5)**2  )/(ONE/BETA0 +id2%v(5))
                MR%v(6)=MR%v(6)+ LD*(ONE/BETA-ONE/BETA0)
             endif

             if(time==1) then
                id2=1
                id2%v(5)=(TWO*id2%v(5)+id2%v(5)**2)/(SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )+ONE/BETA0)
                mr=id2**(-1).o.mr.o.id2 ;
                id2%v(5)=one.mono.'00001'
                BETA=(ONE+id2%v(5) )/SQRT(ONE/BETA0**2+TWO*id2%v(5)+id2%v(5)**2  )
                mr%v(6)=mr%v(6)-((mr%v(6).sub.'000001').mono.'000001')
                mr%v(6)=beta*(mr%v(6)- LD*(ONE/BETA-ONE/BETA0) )
                mr%v(6)=mr%v(6)+(one.mono.'000001')
             endif

             call SET_TREE(mad_tree_rad_rev,Mr)
          else
             call SET_TREE(mad_tree_rad_rev,M)
          endif

       ELSE  ! MAD_TREE_DELTAMAP
          call  dainput_SPECIAL6(m,mf,time);

          if(symplectify)call symplectic(m,symplectic_eps,nst)

          call SET_TREE(mad_tree_rev,M)
          if(n_map==2) then
             call  dainput_SPECIAL6(mr,mf,time);
             call SET_TREE(mad_tree_rad_rev,Mr)
          else
             call SET_TREE(mad_tree_rad_rev,M)
          endif
       ENDIF ! MAD_TREE_DELTAMAP
       close(mf)
       if(symplectic_print) then
          call kanalnummer(mf)
          open(unit=mf,file=file)
          write(mf,*) n_map,no,ang,ld,MAD_TREE_DELTAMAP,nst   ! number of maps (1,2), no, ld=design length
          IF(MAD_TREE_DELTAMAP) THEN
             read(mf,*) ndpt,time,timefac                  ! npdt, time=0,1 (1 uses time), timefac = +/- 1)
          ENDIF
          call daprint(m,mf)

          close(mf)
       endif
       call kill(m,mr,id,id2); call kill(beta,gamma);


    ELSE  ! PRESENT FILE
       LD=MAD_TREE_LD
       ANG=MAD_TREE_ANGLE
    ENDIF

    IF(ang/=zero) THEN
       Taylor_maptilt%LC=two*SIN(ANG/two)*ld/ang
       Taylor_maptilt%B0=ang/ld                     !COS(ANG/two)/R
       Taylor_maptilt%ld= ld
    ELSE
       Taylor_maptilt%LC=ld
       Taylor_maptilt%B0=zero                      !COS(ANG/two)/R
       Taylor_maptilt%ld= ld
    ENDIF

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:nlp)
       call write_i
       Taylor_maptilt%NAME=NAME(1:nlp)
    ELSE
       Taylor_maptilt%NAME=NAME
    ENDIF

    Taylor_maptilt%KIND=kind22
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          Taylor_maptilt%tilt=t%tilt(1)
       ELSE
          Taylor_maptilt%tilt=t%tilt(0)
       ENDIF
    ENDIF
  END FUNCTION Taylor_maptilt



  FUNCTION  arbitrary_tilt(NAME,file,T,no)
    implicit none
    type (EL_LIST) arbitrary_tilt
    CHARACTER(*), INTENT(IN):: NAME,file
    type (TILTING),optional, INTENT(IN):: T
    real(dp) L,ANGLE,HC
    integer mf,nst,I,ORDER
    integer, optional :: no
    LOGICAL(LP) REPEAT
    TYPE(TAYLOR) B(3),ax(2),ay(2)

    file_fitted=file
    arbitrary_tilt=0

    call kanalnummer(mf)
    open(unit=mf,file=file_fitted)
    read(mf,*) nst,L,ANGLE, ORDER,REPEAT
    if(present(no)) order=no
    CALL INIT(ORDER,2)
    CALL ALLOC(B)
    CALL ALLOC(ax)
    CALL ALLOC(ay)


    ALLOCATE(T_E(NST),T_ax(NST),T_ay(NST))

    DO I=1,NST
       IF(I==1.or.(.not.repeat)) THEN
          CALL READ(B(1),mf);CALL READ(B(2),mf);CALL READ(B(3),mf);
          CALL READ(Ax(1),mf);CALL READ(Ay(1),mf);CALL READ(Ax(2),mf);CALL READ(Ay(2),mf);
       ENDIF
       B(1)=B(1)/BRHO
       B(2)=B(2)/BRHO
       B(3)=B(3)/BRHO
       Ax(1)=Ax(1)/BRHO
       Ax(2)=Ax(2)/BRHO
       Ay(1)=Ay(1)/BRHO
       Ay(2)=Ay(2)/BRHO
       CALL SET_TREE_g(T_E(i),B)
       CALL SET_TREE_g(T_ax(i),ax)
       CALL SET_TREE_g(T_ay(i),ay)
    enddo
    call KILL(B)
    CALL KILL(ax)
    CALL KILL(ay)

    close(MF)



    HC=ANGLE/L


    !    IF(ANG/=zero.AND.R/=zero) THEN
    if(hc>zero) then
       arbitrary_tilt%LC=two*SIN(ANGLE/two)/hc
    else
       arbitrary_tilt%LC=L
    endif
    arbitrary_tilt%B0=hc                     !COS(ANG/two)/R
    arbitrary_tilt%LD=L
    arbitrary_tilt%L=arbitrary_tilt%LD

    IF(LEN(NAME)>nlp) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72,/),(1x,a72))'
       w_p%c(1)=name
       WRITE(w_p%c(2),'(a17,1x,a16)') ' IS TRUNCATED TO ', NAME(1:16)
       call write_i
       arbitrary_tilt%NAME=NAME(1:16)
    ELSE
       arbitrary_tilt%NAME=NAME
    ENDIF

    arbitrary_tilt%nst=NST
    arbitrary_tilt%KIND=KINDPA
    IF(PRESENT(t)) then
       IF(T%NATURAL) THEN
          arbitrary_tilt%tilt=t%tilt(1)
       ELSE
          arbitrary_tilt%tilt=t%tilt(0)
       ENDIF
    ENDIF
  END FUNCTION arbitrary_tilt
  ! linked


  SUBROUTINE  EQUAL_L(R,S1)
    implicit none
    type (layout),INTENT(inOUT)::R
    type (layout),INTENT(IN)::S1
    INTEGER I
    !    real(dp) gamma0I,gamBET
    TYPE (fibre), POINTER :: C   !,fitted
    !   logical(lp) firstfitted
    Nullify(C);    !Nullify(fitted);
    !   firstfitted=.true.
    CALL SET_UP(R)
    !    R%ENERGY=ENERGY
    !    R%KINETIC=KINETIC
    !    R%beta0=beta0
    !    R%brho=BRHO
    !    R%p0c=p0c
    !    gamma0I=SQRT(one-R%beta0**2)
    !    gambet =(gamma0I/R%beta0)**2

    !    R%CIRCUMFERENCE=zero
    c=>s1%start
    DO I=1,S1%N

       CALL APPEND( R, C )
       c=>c%next
    ENDDO


    if(use_info) then
       c=>R%start
       c%i%s=zero
       do i=1,R%n
          if(i<R%n.and.use_info) c%next%i%s=c%i%s+c%mag%p%ld

          c=>c%next
       enddo
    endif

  END SUBROUTINE EQUAL_L



  ! linked
  SUBROUTINE Set_Up_MAD( L )
    implicit none
    TYPE (layout) L
    NULLIFY(L%closed);  NULLIFY(L%lastpos);
    NULLIFY(L%NTHIN);NULLIFY(L%THIN);
    !    NULLIFY(L%ENERGY);NULLIFY(L%KINETIC);
    !    NULLIFY(L%P0C);NULLIFY(L%BRHO);NULLIFY(L%BETA0);
    NULLIFY(L%n);
    !    NULLIFY(L%circumference);
    allocate(l%n); l%n=0;
    allocate(l%closed); l%closed=.false.;
    allocate(l%lastpos); l%lastpos=0;

    NULLIFY( L % last )       ! layout is empty at first
    NULLIFY( L % end )       ! layout is empty at first
    NULLIFY( L % start )       ! layout is empty at first
    NULLIFY( L % start_ground )       ! layout is empty at first
    NULLIFY( L % PARENT_UNIVERSE )       ! layout is empty at first
  END SUBROUTINE Set_Up_MAD


  SUBROUTINE  EQUAL_L_L(R,S1)
    implicit none
    logical(lp) :: doneitt=.true.
    type (layout),INTENT(inOUT)::R
    type (layout),INTENT(IN)::S1
    INTEGER I
    TYPE (fibre), POINTER :: C

    if(makeit) then
       call equal_l(r,s1)
       r%closed=circular
       circular=.false.
       makeit=.false.
       CALL RING_L(R,doneitt)
       return
    endif

    Nullify(C);

    CALL SET_UP(R)

    c=>s1%start
    DO I=1,S1%N
       call APPEND_mad_like(R,C)
       C=>C%NEXT
    ENDDO

  END SUBROUTINE EQUAL_L_L

  FUNCTION add_EE( S1, S2 )
    implicit none
    TYPE (layout) add_EE
    TYPE (fibre), INTENT (IN) :: S1, S2

    call Set_Up_mad(add_ee)
    call APPEND_mad_like(add_ee,s1)
    call APPEND_mad_like(add_ee,s2)

  END FUNCTION add_EE

  FUNCTION add_EB( S1, S2 )
    implicit none
    TYPE (layout) add_EB
    TYPE (fibre), INTENT (IN) :: S1
    TYPE (layout), INTENT (IN) :: S2
    INTEGER I
    type(fibre), pointer ::c
    nullify(c)
    call Set_Up_MAD(add_EB)
    call APPEND_mad_like(add_EB,s1)

    c=>s2%start
    do i=1,s2%n
       call APPEND_mad_like(add_EB,c)
       c=>c%next
    enddo

  END FUNCTION add_EB

  FUNCTION add_BE( S2 , S1 )
    implicit none
    TYPE (layout) add_BE
    TYPE (fibre), INTENT (IN) :: S1
    TYPE (layout), INTENT (IN) :: S2
    INTEGER I
    type(fibre), pointer ::c
    nullify(c)
    call Set_Up_MAD(add_BE)

    c=>s2%start
    do i=1,s2%n
       call APPEND_mad_like(add_BE,c)
       c=>c%next
    enddo
    call APPEND_mad_like(add_BE,s1)

  END FUNCTION add_BE

  FUNCTION add_BB( S1 , S2 )
    implicit none
    TYPE (layout) add_BB
    TYPE (layout), INTENT (IN) :: S1
    TYPE (layout), INTENT (IN) :: S2
    INTEGER I
    type(fibre), pointer ::c
    nullify(c)
    call Set_Up_MAD(add_BB)

    c=>s1%start
    do i=1,s1%n
       call APPEND_mad_like(add_BB,c)
       c=>c%next
    enddo
    c=>s2%start
    do i=1,s2%n
       call APPEND_mad_like(add_BB,c)
       c=>c%next
    enddo

  END FUNCTION add_BB

  FUNCTION SUB_BB( S1 , S2 )
    implicit none
    TYPE (layout) SUB_BB
    TYPE (layout), INTENT (IN) :: S1
    TYPE (layout), INTENT (IN) :: S2
    INTEGER I
    type(fibre), pointer ::c
    nullify(c)
    call Set_Up_MAD(SUB_BB)

    c=>s1%start
    do i=1,s1%n
       call APPEND_mad_like(SUB_BB,c)
       c=>c%next
    enddo
    c=>s2%end
    do i=1,s2%n
       call APPEND_mad_like(SUB_BB,c)
       c=>c%previous
    enddo

  END FUNCTION SUB_BB




  FUNCTION MUL_B( S1, S2 )
    implicit none
    TYPE (layout) MUL_B
    integer, INTENT (IN) :: S1
    TYPE (layout), INTENT (IN) :: S2
    INTEGER I,j
    type(fibre), pointer ::c
    nullify(c)
    call Set_Up_MAD(MUL_B)
    if(s1>=0) then
       do j=1,s1
          c=>s2%start
          do i=1,s2%n
             call APPEND_mad_like(MUL_B,c)
             c=>c%next
          enddo
       enddo
    else
       do j=1,-s1
          c=>s2%end
          do i=1,s2%n
             call APPEND_mad_like(MUL_B,c)
             c=>c%previous
          enddo
       enddo
    endif

  END FUNCTION MUL_B

  FUNCTION MUL_E( S1, S2 )
    implicit none
    TYPE (layout) MUL_E
    integer, INTENT (IN) :: S1
    TYPE (fibre), INTENT (IN) :: S2
    INTEGER I
    call Set_Up_MAD(MUL_E)

    do I=1,IABS(s1)
       call APPEND_mad_like(MUL_E,S2)
    enddo

  END FUNCTION MUL_E


  FUNCTION UNARY_SUBB( S1 )
    implicit none
    TYPE (layout) UNARY_SUBB
    TYPE (layout), INTENT (IN) :: S1
    type(fibre), pointer ::c
    integer i
    nullify(c)
    call Set_Up_MAD(UNARY_SUBB)

    c=>s1%end
    do i=1,s1%n
       call APPEND_mad_like(UNARY_SUBB,c)
       c=>c%previous
    enddo

  END FUNCTION UNARY_SUBB

  FUNCTION makeitc( S1 )
    implicit none
    TYPE (layout) makeitc
    TYPE (layout), INTENT (IN) :: S1
    type(fibre), pointer ::c
    integer i
    nullify(c)
    call Set_Up_MAD(makeitc)

    makeit=.true.
    circular=.true.
    c=>s1%start
    do i=1,s1%n
       call APPEND_mad_like(makeitc,c)
       c=>c%next
    enddo

  END FUNCTION makeitc

  FUNCTION makeits( S1 )
    implicit none
    TYPE (layout) makeits
    TYPE (layout), INTENT (IN) :: S1
    type(fibre), pointer ::c
    integer i
    nullify(c)
    call Set_Up_MAD(makeits)

    makeit=.true.
    circular=.false.
    c=>s1%start
    do i=1,s1%n
       call APPEND_mad_like(makeits,c)
       c=>c%next
    enddo

  END FUNCTION makeits




end module Mad_like
