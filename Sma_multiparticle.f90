module ptc_multiparticle
  use S_TRACKING ,FRINGE_=>FRINGE__MULTI,FACE=>FACE_MULTI

  implicit none
  public
  private index
  PRIVATE fringer_STRAIGHT,FRINGE_CAV4,fringer_STREX,fringer_TEAPOT
  PRIVATE FRINGE_CAV_TRAV,ADJUST_TIME_CAV_TRAV_OUT,INTER_CAV_TRAV,ADJUST_WI
  PRIVATE INTER_DRIFT1,INTER_dkd2,INTER_KICKT3,INTER_CAV4,INTER_SOL5,INTER_KTK,INTER_TEAPOT &
       ,INTER_STREX,INTER_SOLT,INTER_NSMI,INTER_SSMI,INTER_MON_11_14,INTER_ESEPTUM,INTER_RCOL18 &
       ,INTER_ECOL19,INTER_WI
  PRIVATE TRACK_THIN_LAYOUT_1,TRACK_THIN_LAYOUT_12,TRACK_LAYOUT_12,TRACK_LAYOUT_S12
  INTEGER, PARAMETER :: CASE1=1,CASE2=2, CASE0=0, CASEP1=-1,CASEP2=-2
  INTEGER, PRIVATE :: TOTALPATH_FLAG
  integer :: index =0
  CHARACTER*27 CASE_NAME(-2:2)
  PRIVATE fuzzy_eq,fuzzy_neq




  INTERFACE TRACK_SLICE
     MODULE PROCEDURE INTER_DRIFT1
     MODULE PROCEDURE INTER_dkd2
     MODULE PROCEDURE INTER_KICKT3
     MODULE PROCEDURE INTER_CAV4
     MODULE PROCEDURE INTER_SOL5
     MODULE PROCEDURE INTER_KTK
     MODULE PROCEDURE INTER_TKTF
     MODULE PROCEDURE INTER_NSMI
     MODULE PROCEDURE INTER_SSMI
     MODULE PROCEDURE INTER_TEAPOT
     MODULE PROCEDURE INTER_STREX
     MODULE PROCEDURE INTER_SOLT
     MODULE PROCEDURE INTER_MON_11_14
     MODULE PROCEDURE INTER_ESEPTUM
     MODULE PROCEDURE INTER_CAV_TRAV
     MODULE PROCEDURE INTER_RCOL18
     MODULE PROCEDURE INTER_ECOL19
     MODULE PROCEDURE INTER_WI
  END INTERFACE

  INTERFACE TRACK_FRINGE
     MODULE PROCEDURE fringer_STRAIGHT   ! GOOD FOR DKD2,KTK, AND TKTF
     MODULE PROCEDURE FRINGE_CAV4
     MODULE PROCEDURE fringer_TEAPOT
     MODULE PROCEDURE fringer_strex
     MODULE PROCEDURE FRINGE_CAV_TRAV
  END INTERFACE

  INTERFACE TRACK_THIN_LAYOUT_S
     MODULE PROCEDURE TRACK_THIN_LAYOUT_12
     MODULE PROCEDURE TRACK_THIN_LAYOUT_1
  END INTERFACE


  INTERFACE TRACK_LAYOUT_USING_THIN_S
     MODULE PROCEDURE TRACK_LAYOUT_12
     MODULE PROCEDURE TRACK_LAYOUT_S12
  END INTERFACE

  INTERFACE TRACK
     MODULE PROCEDURE TRACK_LAYOUT_12
     MODULE PROCEDURE TRACK_LAYOUT_S12
     MODULE PROCEDURE TRACK_THIN_T
  END INTERFACE

  INTERFACE COPY
     MODULE PROCEDURE COPY_BEAM
  END INTERFACE

  INTERFACE OPERATOR (.feq.)
     MODULE PROCEDURE fuzzy_eq
  END INTERFACE
  INTERFACE OPERATOR (.fne.)
     MODULE PROCEDURE fuzzy_neq
  END INTERFACE


CONTAINS

  FUNCTION fuzzy_eq( S1, S2 )
    implicit none
    logical(lp) fuzzy_eq
    real(dp), INTENT (IN) :: S1,S2
    fuzzy_eq=.false.

    if(abs(s1-s2)<=c_%eps_pos) fuzzy_eq=.true.

  end FUNCTION fuzzy_eq

  FUNCTION fuzzy_neq( S1, S2 )
    implicit none
    logical(lp) fuzzy_neq
    real(dp), INTENT (IN) :: S1,S2
    fuzzy_neq=.false.

    if(abs(s1-s2)>c_%eps_pos) fuzzy_neq=.true.

  end FUNCTION fuzzy_neq



  SUBROUTINE move_to_s( L,s,current,i,ds ) ! Moves current to the i^th position
    implicit none
    TYPE (THIN_LENS), POINTER :: Current
    TYPE (THIN_layout) L
    real(dp) s,sp,ds
    integer i,k
    logical(lp) DOIT !,track_it

    !  track_it=.false.
    sp=mod(s,L%END%S(3))

    if(sp==zero.and.s/=zero) then
       current=>l%end
       i=l%n+1
       ds=zero
       !  track_it=.true.
       return
    endif

    if(sp==zero) then
       current=>l%start
       i=1
       ds=zero
       return
    endif


    nullify(current);
    Current => L%LAST

    k=L%LASTPOS
    I=K
    ds=zero

    IF(SP>CURRENT%S(3) ) then

       do i=k,l%n-1
          if(current%next%s(3)>=sp) exit
          current=>current%next
       enddo

       if(current%next%s(3)/=sp) ds=sp-current%s(3)
       if(current%next%s(3)==sp) THEN
          ds=ZERO
          CURRENT=>current%next
          I=I+1
       ENDIF

    elseif(SP<CURRENT%S(3)) then

       do i=k-1,1,-1
          current=>current%previous
          if(current%s(3)<=sp) exit
       enddo

       ds=sp-current%s(3)

    endif

    L%LASTPOS=I; L%LAST => Current;

    if(ds>zero) then
       if(CURRENT%S(4)-ds.feq.zero) then
          ds=zero
          current=>Current%next
          i=i+1
          L%LAST => Current;
       ELSEIF(ds.feq.zero) THEN
          DS=ZERO
       ENDIF
    endif


    DOIT=.TRUE.
    !DOIT=.FALSE.

    if(iabs(CURRENT%cas)==0.OR.iabs(CURRENT%cas)==1) then
       do while(DS==ZERO.AND.DOIT)    ! PUTS AT BEGINNING IF DS=ZERO
          CURRENT=>CURRENT%PREVIOUS
          IF(ASSOCIATED(CURRENT)) THEN

             IF((SP.FNE.CURRENT%S(3)).or.CURRENT%cas==-2) THEN
                CURRENT=>CURRENT%NEXT
                DOIT=.FALSE.
             ELSE
                I=I-1
             ENDIF

          ELSE
             CURRENT=>CURRENT%NEXT
             DOIT=.FALSE.
          ENDIF

       enddo
    elseif(iabs(CURRENT%cas)==2) then
       do while(DS==ZERO.AND.DOIT)    ! PUTS AT BEGINNING IF DS=ZERO
          CURRENT=>CURRENT%next
          IF(ASSOCIATED(CURRENT)) THEN

             IF((SP.FNE.CURRENT%S(3)).or.CURRENT%cas==1) THEN
                CURRENT=>CURRENT%previous
                DOIT=.FALSE.
             ELSE
                I=I+1
             ENDIF

          ELSE
             CURRENT=>CURRENT%previous
             DOIT=.FALSE.
          ENDIF

       enddo
    endif



    L%LASTPOS=I; L%LAST => Current;

    IF(I/=L%LAST%POS) THEN
       WRITE(6,*) " ERROR IN move_to_s ",I,L%LAST%POS
       STOP 999
    ENDIF

  END SUBROUTINE move_to_s



  ! tracking fringe areas

  SUBROUTINE fringer_STRAIGHT(EL,EL5,EL6,EL7,EL17,X,J)
    IMPLICIT NONE
    TYPE(DKD2),OPTIONAL,INTENT(IN):: EL
    TYPE(SOL5),OPTIONAL,INTENT(INOUT):: EL5
    TYPE(KTK),OPTIONAL,INTENT(INOUT):: EL6
    TYPE(TKTF),OPTIONAL,INTENT(INOUT):: EL7
    TYPE(SOLT),OPTIONAL,INTENT(INOUT):: EL17
    !      TYPE(BEAM), INTENT(INOUT) ::B
    integer,INTENT(IN):: J
    real(dp), INTENT(INOUT) :: X(6)
    integer  i
    ! J=1 front
    IF(PRESENT(EL)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL%P%DIR==1) THEN
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ELSE
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF
       else
          if(EL%P%DIR==1) THEN
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
          ENDIF
       ENDIF
       !         CALL X_IN_BEAM(X,B,I)

    ELSEIF(PRESENT(EL6)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL6%P%DIR==1) THEN
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,1,X)
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,1,X)
          ELSE
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,2,X)
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,2,X)
          ENDIF
       else
          if(EL6%P%DIR==1) THEN
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,2,X)
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,2,X)
          ELSE
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,1,X)
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,1,X)
          ENDIF
       ENDIF
       !         CALL X_IN_BEAM(X,B,I)

    ELSEIF(PRESENT(EL5)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL5%P%DIR==1) THEN
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,1,X)
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,1,X)
          ELSE
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,2,X)
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,2,X)
          ENDIF
       else
          if(EL5%P%DIR==1) THEN
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,2,X)
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,2,X)
          ELSE
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,1,X)
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,1,X)
          ENDIF
       ENDIF

    ELSEIF(PRESENT(EL7)) THEN

       if(J==1) then
          if(EL7%P%DIR==1) THEN
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,1,X)
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,1,X)
          ELSE
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,2,X)
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,2,X)
          ENDIF
       else
          if(EL7%P%DIR==1) THEN
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,2,X)
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,2,X)
          ELSE
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,1,X)
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,1,X)
          ENDIF
       ENDIF

    ELSEIF(PRESENT(EL17)) THEN

       if(J==1) then
          if(EL17%P%DIR==1) THEN
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,1,X)
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,1,X)
          ELSE
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,2,X)
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,2,X)
          ENDIF
       else
          if(EL17%P%DIR==1) THEN
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,2,X)
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,2,X)
          ELSE
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,1,X)
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,1,X)
          ENDIF
       ENDIF

    ENDIF
  END SUBROUTINE fringer_STRAIGHT

  SUBROUTINE FRINGE_CAV4(EL,X,J)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) ::  X(6)
    TYPE(CAV4),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i,JC
    REAL(DP) C1,S1,V,O

    JC=-2*J+3

    IF(EL%P%NOCAVITY) RETURN
    IF(.NOT.EL%P%FRINGE) RETURN
    IF(EL%THIN) RETURN
    IF(jC==1.AND.EL%P%KILL_ENT_FRINGE) RETURN
    IF(jC==-1.AND.EL%P%KILL_EXI_FRINGE) RETURN


    O=EL%freq*twopi/CLIGHT
    V=jC*EL%P%CHARGE*EL%volt*c_1d_3/EL%P%P0C



    C1=COS(O*(x(6))+EL%PHAS+phase0)
    S1=SIN(O*(x(6))+EL%PHAS+phase0)



    X(2)=X(2)+V*S1*X(1)
    X(4)=X(4)+V*S1*X(3)
    x(5)=x(5)-HALF*(X(1)**2+X(3)**2)*V*C1*O




  END SUBROUTINE FRINGE_CAV4

  SUBROUTINE ADJUST_TIME_CAV4(EL,X,J)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT)::  X(6)
    TYPE(CAV4),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i




    IF(J==1) THEN
       IF(EL%P%NOCAVITY) RETURN

       IF(EL%THIN) THEN
          CALL CAVITY(EL,X)
          RETURN
       ENDIF


    ELSE
       IF(EL%THIN) RETURN


       if(EL%P%TIME) then
          X(6)=X(6)-(CAVITY_TOTALPATH-EL%P%TOTALPATH)*EL%P%LD/EL%P%BETA0
       else
          X(6)=X(6)-(CAVITY_TOTALPATH-EL%P%TOTALPATH)*EL%P%LD
       endif

    ENDIF

  END SUBROUTINE ADJUST_TIME_CAV4



  SUBROUTINE fringer_TEAPOT(EL,X,J)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(TEAPOT),INTENT(IN):: EL
    integer,INTENT(IN):: J
    INTEGER I



    IF(EL%P%DIR==1) THEN

       IF(J==1) THEN

          IF(EL%P%EDGE(1)/=zero) THEN
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) then
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
                x(2)=x(2)+EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)-EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)+EL%P%EDGE(1)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)-EL%P%EDGE(1)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL WEDGE(-EL%P%EDGE(1),X,EL2=EL)
          ELSE
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ENDIF

       ELSE ! J=2

          IF(EL%P%EDGE(2)/=zero) THEN
             CALL WEDGE(-EL%P%EDGE(2),X,EL2=EL)
             IF(EL%P%FRINGE) then
                x(2)=x(2)+EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)-EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)+EL%P%EDGE(2)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)-EL%P%EDGE(2)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
          ENDIF
       ENDIF ! J=2
    ELSE


       IF(J==1) THEN

          IF(EL%P%EDGE(2)/=zero) THEN
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) then
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
                x(2)=x(2)-EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)+EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)-EL%P%EDGE(2)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)+EL%P%EDGE(2)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL WEDGE(-EL%P%EDGE(2),X,EL2=EL)
          ELSE
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF

       ELSE ! J=2

          IF(EL%P%EDGE(1)/=zero) THEN
             CALL WEDGE(-EL%P%EDGE(1),X,EL2=EL)
             IF(EL%P%FRINGE) then
                x(2)=x(2)-EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)+EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)-EL%P%EDGE(1)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)+EL%P%EDGE(1)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
          ENDIF

       ENDIF

    ENDIF



  END SUBROUTINE fringer_TEAPOT


  SUBROUTINE  fringer_STREX(EL,X,J)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    integer,INTENT(IN):: J
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(STREX),INTENT(IN):: EL
    real(dp) ANGH
    INTEGER I
    ! J=1 front



    IF(EL%P%DIR==1) THEN

       IF(J==1) THEN

          IF(EL%LIKEMAD) THEN

             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(1)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL  WEDGE(ANGH,X,EL1=EL)

          ELSE
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ENDIF

       ELSE  ! J==2


          IF(EL%LIKEMAD) THEN
             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(2)
             CALL  WEDGE(ANGH,X,EL1=EL)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
          ENDIF

       ENDIF ! J

    ELSE

       IF(J==1) THEN

          IF(EL%LIKEMAD) THEN

             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(2)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL  WEDGE(ANGH,X,EL1=EL)
          ELSE
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF


       ELSE ! J==2

          IF(EL%LIKEMAD) THEN
             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(1)
             CALL  WEDGE(ANGH,X,EL1=EL)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
          ENDIF

       ENDIF ! J


    ENDIF




  END SUBROUTINE fringer_STREX







  SUBROUTINE  FRINGE_CAV_TRAV(EL,X,J)
    IMPLICIT NONE
    !      TYPE(BEAM), INTENT(INOUT) ::B
    integer,INTENT(IN):: J
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAV),INTENT(INOUT):: EL
    INTEGER I
    ! J=1 front
    IF(J==1) THEN



       CALL FRINGECAV_TRAV(EL,EL%P%DIR,X)



    ELSE



       CALL FRINGECAV_TRAV(EL,-EL%P%DIR,X)



    ENDIF

  END SUBROUTINE FRINGE_CAV_TRAV

  SUBROUTINE ADJUST_TIME_CAV_TRAV_OUT(EL,X,J)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAV),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i

    IF(J==1) RETURN


    if(EL%P%TIME) then
       X(6)=X(6)-(1-EL%P%TOTALPATH)*EL%P%LD/EL%P%BETA0
    else
       X(6)=X(6)-(1-EL%P%TOTALPATH)*EL%P%LD
    endif



  END SUBROUTINE ADJUST_TIME_CAV_TRAV_OUT

  ! tracking one steps in the body


  SUBROUTINE INTER_DRIFT1(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(DRIFT1),INTENT(IN):: EL

    real(dp) DH,DD

    SELECT CASE(EL%P%METHOD)
    CASE(2,4,6)
       DH=EL%L/EL%P%NST
       DD=EL%P%LD/EL%P%NST




       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_DRIFT1

  SUBROUTINE INTER_dkd2 (EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(DKD2),INTENT(IN):: EL

    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST




       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,D,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(4)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST




       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK2,X)

       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(6)
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)

       DO J=4,2,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICK (EL,DK(J),X)
       ENDDO
       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK(1),X)

       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       DO J=2,4
          CALL KICK (EL,DK(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO




    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_dkd2

  SUBROUTINE INTER_KICKT3(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(KICKT3),INTENT(IN):: EL







    CALL TRACK(EL,X)





  END SUBROUTINE INTER_KICKT3

  SUBROUTINE INTER_NSMI(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(NSMI),INTENT(IN):: EL







    CALL TRACK(EL,X)




  END SUBROUTINE INTER_NSMI

  SUBROUTINE INTER_SSMI(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(SSMI),INTENT(IN):: EL

    INTEGER I






    CALL TRACK(EL,X)




  END SUBROUTINE INTER_SSMI

  SUBROUTINE INTER_CAV4(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  X(6)

    TYPE(CAV4),INTENT(INOUT):: EL
    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J


    TOTALPATH_FLAG=EL%P%TOTALPATH
    EL%P%TOTALPATH=CAVITY_TOTALPATH



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,D,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


    CASE(4)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK2,X)

       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(6)
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       DO J=4,2,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKCAV (EL,DK(J),X)
       ENDDO
       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK(1),X)

       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       DO J=2,4
          CALL KICKCAV(EL,DK(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO




    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    EL%P%TOTALPATH=TOTALPATH_FLAG


  END SUBROUTINE INTER_CAV4

  SUBROUTINE INTER_SOL5(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  X(6)
    TYPE(SOL5),INTENT(IN):: EL

    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2,D2H
    real(dp) dd1,dd2,DK(4),DF(4),DDF(4)
    INTEGER I,J




    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=(EL%P%LD)/two/EL%P%NST


       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,DH,X)
       CALL KICK_SOL(EL,DH,X)
       CALL KICKMUL(EL,D,X)

       CALL KICK_SOL(EL,DH,X)
       CALL SOL_ROT (EL,DH,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(4)

       D=EL%L/EL%P%NST
       D1=D*FD1
       D2=D*FD2
       DK1=D*FK1
       DK2=D*FK2
       D2H=DK2/two
       DD1=(EL%P%LD)/EL%P%NST*FD1
       DD2=(EL%P%LD)/EL%P%NST*FD2

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL KICKMUL(EL,DK1,X)

       CALL KICK_SOL(EL,D1,X)
       CALL SOL_ROT (EL,D1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D2H,X)
       CALL KICK_SOL(EL,D2H,X)
       CALL KICKMUL(EL,DK2,X)
       CALL KICK_SOL(EL,D2H,X)
       CALL SOL_ROT (EL,D2H,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL KICKMUL(EL,DK1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL SOL_ROT (EL,D1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(6)
       DO I =1,4
          DK(I)=EL%L*YOSK(I)/EL%P%NST
          DF(I)=DK(I)/two
          DDF(I)=EL%P%LD*YOSK(I)/two/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       DO J=4,1,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL KICKMUL(EL,DK(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO

       DO J=2,4
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL KICKMUL(EL,DK(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT



  END SUBROUTINE INTER_SOL5


  SUBROUTINE INTER_KTK(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(KTK),INTENT(INOUT):: EL

    INTEGER I
    real(dp) DK,DK2,DK6,DK4,DK5



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       IF(OLD_IMPLEMENTATION_OF_SIXTRACK) THEN
          DK2=EL%L/EL%P%NST
          DK=DK2/two


          CALL GETMAT(EL,X)
          CALL KICKKTK(EL,DK,X)  ! NEW
          CALL KICKPATH(EL,DK,X)
          CALL PUSHKTK(EL,X)

          CALL KICKPATH(EL,DK,X)
          if(EL%P%RADIATION) CALL GETMAT(EL,X)
          CALL KICKKTK(EL,DK,X)  ! NEW
       ELSE  ! OLD_IMPLEMENTATION_OF_SIXTRACK
          CALL GETMAT(EL,X)
          CALL PUSHKTK(EL,X)
          CALL KICKPATH(EL,DK,X)
          CALL KICKKTK(EL,DK2,X)  ! NEW
          CALL KICKPATH(EL,DK,X)
          CALL PUSHKTK(EL,X)
          if(EL%P%RADIATION) CALL GETMAT(EL,X)

       ENDIF ! OLD_IMPLEMENTATION_OF_SIXTRACK

    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two
       CALL GETMAT(EL,X)
       CALL KICKKTK(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKKTK(EL,DK6,X)
       CALL KICKPATH(EL,DK2,X)
       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKKTK(EL,DK,X) ! NEW
       if(EL%P%RADIATION) CALL GETMAT(EL,X)
    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two


       CALL GETMAT(EL,X)
       CALL KICKKTK(EL,DK,X)   ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHKTK(EL,X)

       CALL KICKKTK(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHKTK(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKKTK(EL,DK6,X)   ! SYMMETRY POINT
       CALL KICKPATH(EL,DK5,X)

       CALL PUSHKTK(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKKTK(EL,DK4,X)

       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKKTK(EL,DK,X)   ! NEW
       if(EL%P%RADIATION) CALL GETMAT(EL,X)

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_KTK

  SUBROUTINE INTER_TKTF(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(TKTF),INTENT(INOUT):: EL
    INTEGER I
    real(dp) DK,DK2,DK6,DK4,DK5



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DK2=EL%L/EL%P%NST
       DK=DK2/two

       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK2,X)

       CALL KICKPATH(EL,DK,X)
       CALL PUSHTKT7(EL,X)

    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two

       CALL KICKTKT7(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKTKT7(EL,DK6,X)

       CALL KICKPATH(EL,DK2,X)
       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK,X) ! NEW

    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two

       CALL KICKTKT7(EL,DK,X)  ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKTKT7(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKTKT7(EL,DK6,X)   ! SYMMETRY POINT

       CALL KICKPATH(EL,DK5,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKTKT7(EL,DK4,X)

       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK,X)  ! NEW

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_TKTF

  SUBROUTINE INTER_TEAPOT(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(TEAPOT),INTENT(IN):: EL
    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J


    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST

       CALL SSECH1(EL,DH,DD,X)
       CALL SKICK(EL,D,X)
       CALL SSECH1(EL,DH,DD,X)
    CASE(4)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST

       CALL SSECH1(EL,D1,DD1,X)
       CALL SKICK (EL,DK1,X)
       CALL SSECH1(EL,D2,DD2,X)
       CALL SKICK (EL,DK2,X)

       CALL SSECH1(EL,D2,DD2,X)
       CALL SKICK (EL,DK1,X)
       CALL SSECH1(EL,D1,DD1,X)



    CASE(6)
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       DO J=4,2,-1
          CALL SSECH1(EL,DF(J),DDF(J),X)
          CALL SKICK (EL,DK(J),X)
       ENDDO
       CALL SSECH1(EL,DF(1),DDF(1),X)
       CALL SKICK (EL,DK(1),X)

       CALL SSECH1(EL,DF(1),DDF(1),X)
       DO J=2,4
          CALL SKICK (EL,DK(J),X)
          CALL SSECH1(EL,DF(J),DDF(J),X)
       ENDDO



    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_TEAPOT

  SUBROUTINE INTER_STREX(EL,X)
    IMPLICIT NONE
    TYPE(STREX),INTENT(IN):: EL
    real(dp), INTENT(INOUT) ::  X(6)
    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J


    IF(EL%DRIFTKICK) THEN

       SELECT CASE(EL%P%METHOD)
       CASE(2)
          DH=EL%L/two/EL%P%NST
          D=EL%L/EL%P%NST
          DD=EL%P%LD/two/EL%P%NST




          CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,D,X)
          CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


       CASE(4)
          D1=EL%L*FD1/EL%P%NST
          D2=EL%L*FD2/EL%P%NST
          DD1=EL%P%LD*FD1/EL%P%NST
          DD2=EL%P%LD*FD2/EL%P%NST
          DK1=EL%L*FK1/EL%P%NST
          DK2=EL%L*FK2/EL%P%NST




          CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK1,X)
          CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK2,X)
          CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK1,X)
          CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



       CASE(6)
          DO I =1,4
             DF(I)=EL%L*YOSD(I)/EL%P%NST
             DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
             DK(I)=EL%L*YOSK(I)/EL%P%NST
          ENDDO




          DO J=4,2,-1
             CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
             CALL KICKEX (EL,DK(J),X)
          ENDDO
          CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK(1),X)
          CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          DO J=2,4
             CALL KICKEX (EL,DK(J),X)
             CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


          ENDDO


       CASE DEFAULT
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A72))'
          WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
          call write_e(357)
       END SELECT
    ELSE
       SELECT CASE(EL%P%METHOD)
       CASE(2)
          DH=EL%L/two/EL%P%NST
          D=EL%L/EL%P%NST
          DD=EL%P%LD/two/EL%P%NST




          CALL SPAR(EL,DH,DD,X)
          CALL KICKEX (EL,D,X)
          CALL SPAR(EL,DH,DD,X)


       CASE(4)
          D1=EL%L*FD1/EL%P%NST
          D2=EL%L*FD2/EL%P%NST
          DD1=EL%P%LD*FD1/EL%P%NST
          DD2=EL%P%LD*FD2/EL%P%NST
          DK1=EL%L*FK1/EL%P%NST
          DK2=EL%L*FK2/EL%P%NST




          CALL SPAR(EL,D1,DD1,X)
          CALL KICKEX (EL,DK1,X)
          CALL SPAR(EL,D2,DD2,X)
          CALL KICKEX (EL,DK2,X)
          CALL SPAR(EL,D2,DD2,X)
          CALL KICKEX (EL,DK1,X)
          CALL SPAR(EL,D1,DD1,X)



       CASE(6)
          DO I =1,4
             DF(I)=EL%L*YOSD(I)/EL%P%NST
             DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
             DK(I)=EL%L*YOSK(I)/EL%P%NST
          ENDDO




          DO J=4,2,-1
             CALL SPAR(EL,DF(J),DDF(J),X)
             CALL KICKEX (EL,DK(J),X)
          ENDDO
          CALL SPAR(EL,DF(1),DDF(1),X)
          CALL KICKEX (EL,DK(1),X)
          CALL SPAR(EL,DF(1),DDF(1),X)
          DO J=2,4
             CALL KICKEX (EL,DK(J),X)
             CALL SPAR(EL,DF(J),DDF(J),X)
          ENDDO





       CASE DEFAULT
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A72))'
          WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
          call write_e(357)
       END SELECT

    ENDIF



  END SUBROUTINE INTER_STREX

  SUBROUTINE INTER_SOLT(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(SOLT),INTENT(INOUT):: EL

    INTEGER I
    real(dp) DK,DK2,DK4,DK5,DK6,DKT



    DKT=EL%L/EL%P%NST

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DK2=EL%L/EL%P%NST
       DK=DK2/two





       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)  ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHSOL(EL,X)

       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DK2,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X)  ! NEW





    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two




       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKMUL(EL,DK6,X)
       CALL KICKPATH(EL,DK2,X)
       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X) ! NEW
       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DKT,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF



    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two





       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)   ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHSOL(EL,X)

       CALL KICKMUL(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHSOL(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKMUL(EL,DK6,X)   ! SYMMETRY POINT
       CALL KICKPATH(EL,DK5,X)

       CALL PUSHSOL(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKMUL(EL,DK4,X)

       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X)   ! NEW
       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DKT,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF







    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_SOLT

  SUBROUTINE INTER_MON_11_14(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(mon),INTENT(INOUT):: EL

    INTEGER I






    CALL MONTI(EL,X,I)




  END SUBROUTINE INTER_MON_11_14

  SUBROUTINE INTER_ESEPTUM(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(ESEPTUM),INTENT(INOUT):: EL

    INTEGER I






    CALL SEPTTRACK(EL,X,0)




  END SUBROUTINE INTER_ESEPTUM

  SUBROUTINE INTER_CAV_TRAV(EL,X,Z)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAV),INTENT(INOUT):: EL

    real(dp) , INTENT(IN) :: Z
    real(dp) DF(4),DK(4),DDF(4),D1
    INTEGER I,J
    REAL(DP) Z0
    INTEGER TOTALPATH

    Z0=Z

    TOTALPATH=EL%P%TOTALPATH
    EL%P%TOTALPATH=1
    SELECT CASE(EL%P%METHOD)
    CASE(2)
       D1=EL%L/EL%P%NST



       call rk2_cav(z0,d1,el,x)


    CASE(4)

       D1=EL%L/EL%P%NST



       call rk4_cav(z0,d1,el,x)



    CASE(6)

       D1=EL%L/EL%P%NST



       call rk6_cav(z0,d1,el,x)



    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    !    IF(EL%P%FRINGE)

    EL%P%TOTALPATH=TOTALPATH

  END SUBROUTINE INTER_CAV_TRAV

  SUBROUTINE INTER_RCOL18(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(RCOL),INTENT(INOUT):: EL
    INTEGER I

    CALL RCOLLIMATORI(EL,X,I)
  END SUBROUTINE INTER_RCOL18


  SUBROUTINE INTER_ECOL19(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(ECOL),INTENT(INOUT):: EL

    INTEGER I

    CALL ECOLLIMATORI(EL,X,I)

  END SUBROUTINE INTER_ECOL19


  SUBROUTINE INTER_WI(EL,X,Z)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::X(6)
    TYPE(sagan),INTENT(INOUT):: EL

    real(dp) , INTENT(IN) :: Z
    INTEGER I

    call INTR_SAGAN(EL,X,Z)

  END SUBROUTINE INTER_WI


  SUBROUTINE ADJUST_WI(EL,X,J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(sagan),INTENT(INOUT):: EL

    INTEGER, INTENT(IN) :: J
    INTEGER I

    IF(J==1) RETURN

    X(1)=X(1)-EL%INTERNAL(1)
    X(2)=X(2)-EL%INTERNAL(2)

  END SUBROUTINE ADJUST_WI





  ! MULTIPARTICLE AT THE FIBRE LEVEL

  ! front patch/misaglinments/tilt
  SUBROUTINE TRACK_FIBRE_FRONTR(C,X,K)
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    !    TYPE(BEAM),TARGET,INTENT(INOUT):: B
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    REAL(DP) ENT(3,3), A(3)
    real(dp) xp
    INTEGER I

    ! DIRECTIONAL VARIABLE
    !    C%MAG%P%DIR=>C%DIR
    !    if(present(charge)) then
    !       C%MAG%P%CHARGE=>CHARGE
    !    else
    !       charge1=1
    !       C%MAG%P%CHARGE=>CHARGE1
    !    endif
    !


    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

    ! PUSHING BEAM
    !



    IF(PATCHE/=0.AND.PATCHE/=2) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAG%P%P0C
             B0=>CN%MAG%P%BETA0

             X(2)=X(2)*P0/C%MAG%P%P0C
             X(4)=X(4)*P0/C%MAG%P%P0C
             IF(C%MAG%P%TIME)THEN
                X(5)=root(one+two*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-one !X(5) = DP/P0C_NEW
                X(5)=(two*X(5)+X(5)**2)/(root(one/C%MAG%P%BETA0**2+two*X(5)+X(5)**2)+one/C%MAG%P%BETA0)
             ELSE
                X(5)=(one+X(5))*P0/C%MAG%P%P0C-one
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG/=0.AND.PATCHG/=2) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%a_T
    ENDIF

    CALL DTILTD(C%DIR,C%MAG%P%TILTD,1,X)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       ou = K%EXACTMIS.or.C%MAG%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF




  END SUBROUTINE TRACK_FIBRE_FRONTR

  ! back patch/misaglinments/tilt

  SUBROUTINE TRACK_FIBRE_BACKR(C,X,K)
    implicit none
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    !    TYPE(BEAM),TARGET,INTENT(INOUT):: B
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    INTEGER I


    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAG%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%DIR,C%MAG%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%b_T
    ENDIF

    IF(PATCHG/=0.AND.PATCHG/=1) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely

    IF(PATCHE/=0.AND.PATCHE/=1) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAG%P%P0C
       B0=>CN%MAG%P%BETA0
       X(2)=X(2)*C%MAG%P%P0C/P0
       X(4)=X(4)*C%MAG%P%P0C/P0
       IF(C%MAG%P%TIME)THEN
          X(5)=root(one+two*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAG%P%P0C/P0-one !X(5) = DP/P0C_NEW
          X(5)=(two*X(5)+X(5)**2)/(root(one/B0**2+two*X(5)+X(5)**2)+one/B0)
       ELSE
          X(5)=(one+X(5))*C%MAG%P%P0C/P0-one
       ENDIF
    ENDIF




  END SUBROUTINE TRACK_FIBRE_BACKR

  ! thin lens tracking

  SUBROUTINE TRACK_THIN_SINGLE(T,X,K,CHARGE)
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(THIN_LENS), TARGET, INTENT(INOUT):: T
    REAL(DP),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    INTEGER, TARGET :: CHARGE
    INTEGER I


    T%PARENT_FIBRE%MAG=K
    T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR
    T%PARENT_FIBRE%MAG%P%CHARGE=>CHARGE

    SELECT CASE(T%CAS)
    CASE(CASEP1)
       CALL TRACK_FIBRE_FRONTR(T%PARENT_FIBRE,X,K)
    CASE(CASEP2)
       CALL TRACK_FIBRE_BACKR(T%PARENT_FIBRE,X,K)

    CASE(CASE1,CASE2)
       el=>T%PARENT_FIBRE%MAG
       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=X,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,X,1)
             CALL TRACK_FRINGE(EL%C4,X,J=1)
          ELSE
             CALL TRACK_FRINGE(EL%C4,X,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,X,2)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=X,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=X,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=X,J=T%CAS)
       case(KIND10)
          CALL TRACK_FRINGE(EL%TP10,X,T%CAS)
       case(KIND16,KIND20)
          CALL TRACK_FRINGE(EL%K16,X,T%CAS)
       case(KIND17)
          CALL TRACK_FRINGE(EL17=EL%S17,X=X,J=T%CAS)
       case(KIND21)
          CALL TRACK_FRINGE(EL%CAV21,X=X,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)
          CALL ADJUST_WI(EL%WI,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT

    CASE(CASE0)
       el=>T%PARENT_FIBRE%MAG
       SELECT CASE(EL%KIND)
       CASE(KIND0)
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X)
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X)
       case(KIND3)
          CALL TRACK_SLICE(EL%K3,X)
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X)
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X)
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X)
       case(KIND7)
          CALL TRACK_SLICE(EL%T7,X)
       case(KIND8)
          CALL TRACK_SLICE(EL%S8,X)
       case(KIND9)
          CALL TRACK_SLICE(EL%S9,X)
       case(KIND10)
          CALL TRACK_SLICE(EL%TP10,X)
       case(KIND11:KIND14)
          CALL TRACK_SLICE(EL%MON14,X)
       case(KIND15)
          CALL TRACK_SLICE(EL%SEP15,X)
       case(KIND16,KIND20)
          CALL TRACK_SLICE(EL%K16,X)
       case(KIND17)
          CALL TRACK_SLICE(EL%S17,X)
       case(KIND18)
          CALL TRACK_SLICE(EL%RCOL18,X)
       case(KIND19)
          CALL TRACK_SLICE(EL%ECOL19,X)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,T%S(2))
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,T%S(2))
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT


    END SELECT
    T%PARENT_FIBRE%MAG=DEFAULT

  END SUBROUTINE TRACK_THIN_SINGLE


  SUBROUTINE TRACK_THIN_SINGLE_FOR_time(B,I,DT,K)
    ! Tracks a single particle "I" of the beam for a time DT
    ! The particle is a location defined by the thin lens B%POS(I)%THINLENS
    ! and located B%X(I,7) metres in from of that thin lens
    implicit none
    TYPE(THIN_LENS), POINTER:: T
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP), INTENT(IN) :: DT
    REAL(DP) X(6),XT(6),DT0,YL,DT_BEFORE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    LOGICAL(LP) END_OF_LINE
    INTEGER I
    END_OF_LINE=.FALSE.

    IF(.NOT.K%TOTALPATH) STOP 451
    IF(B%U(i)) RETURN

    X=BEAM_IN_X(B,I)
    T=>B%POS(I)%THINLENS
    T%PARENT_FIBRE%MAG=K
    T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR
    T%PARENT_FIBRE%MAG%P%CHARGE=>B%CHARGE


    DT0=X(6)
    CALL DRIFT_BACK_TO_POSITION(T,B%X(I,7),X)
    B%X(I,7)=ZERO
    DT0=X(6)-DT0

    YL=zero

    DO WHILE(DT0<=DT)
       XT=X
       DT_BEFORE=DT0
       !         WRITE(6,*) " POS ",T%s(1),t%pos_in_fibre
       !         WRITE(6,*) " POS ",T%POS,T%CAS,T%PARENT_FIBRE%MAG%NAME
       CALL TRACK_THIN_SINGLE(T,X,K,B%CHARGE)
       DT0=DT0+(X(6)-XT(6))
       T=>T%NEXT
       IF(.NOT.ASSOCIATED(T%NEXT)) THEN
          END_OF_LINE=.TRUE.
          EXIT
       ENDIF
    ENDDO

    IF(.NOT.END_OF_LINE) THEN
       IF(DT0/=DT) THEN
          B%POS(I)%THINLENS=>T%PREVIOUS
          X=XT
          DT0=DT-DT_BEFORE
          !           WRITE(6,*) " DT0 ", DT0
          CALL DRIFT_TO_TIME(T,YL,DT0,X)
       ELSE
          B%POS(I)%THINLENS=>T
       ENDIF
    ELSE


       IF(DT0<DT) THEN
          B%POS(I)%THINLENS=>T%PREVIOUS
          X=XT
          DT0=DT-DT_BEFORE
          CALL DRIFT_TO_TIME(T,YL,DT0,X)
       ELSE
          B%POS(I)%THINLENS=>T
       ENDIF

    ENDIF



    CALL X_IN_BEAM(X,B,I,DL=YL)



    B%TIME_INSTEAD_OF_S=.TRUE.


  END SUBROUTINE TRACK_THIN_SINGLE_FOR_time



  SUBROUTINE TRACK_THIN_SINGLE_FOR_S(T,B,K)
    ! Tracks a full beam across a single thin lens
    ! Position of macroparticle put at next thin lens

    implicit none
    TYPE(THIN_LENS), TARGET, INTENT(INOUT):: T
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP) X(6)
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    INTEGER I

    DO I=1,B%N
       IF(B%U(i)) CYCLE
       X=BEAM_IN_X(B,I)

       CALL TRACK_THIN_SINGLE(T,X,K,B%CHARGE)
       CALL X_IN_BEAM(X,B,I,DL=ZERO)
       B%POS(I)%THINLENS=>T%NEXT

    ENDDO
    B%TIME_INSTEAD_OF_S=.FALSE.

  END SUBROUTINE TRACK_THIN_SINGLE_FOR_S

  SUBROUTINE TRACK_THIN_T(B,DT,K)
    ! Tracks to full beam for a time DT
    ! All the particles are at different locations
    ! Notice that the layout is hidden: this is consistant with time tracking
    ! Magnets are not ontological objects

    implicit none
    TYPE(BEAM),INTENT(INOUT):: B
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    real(dp),INTENT(IN):: DT

    INTEGER I
    DO I=1,B%N
       IF(B%U(i)) CYCLE
       call TRACK_THIN_SINGLE_FOR_time(B,I,DT,K)
    ENDDO
    B%TIME_INSTEAD_OF_S=.TRUE.

  END SUBROUTINE TRACK_THIN_T

  SUBROUTINE TRACK_LAYOUT_12( R,B,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2 )
    ! Tracks through the thin lens structure R%T of the layout R if it exists.
    ! Several posibilities:
    !1) Pos1 and Pos2 are given :  tracks from thin lens position pos1 to thin lens position pos2
    !2) Thin lens points T1 and T2 are given: Same as above pos1 and pos2 are derived from T1 and T2
    !3) P1 and P2 are fibres: Tracks from the first thin lens of P1 to the first of P2. Results should agree
    !   with plain PTC.
    !4) P1, IN_P1 and P2 , IN_P2 are give: movies to IN_P1 thin lens in P1 and tracks to IN_P2 position in p2
    !  For example, if the input is (P1,1) and (P2,1) the results is same as item #3 (same as plain PTC).

    ! In all the above cases one can elect to give only the first input. Then it tracks one turn around as in
    ! plain PTC.
    ! interfaced as TRACK_LAYOUT_USING_THIN_S
    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE(BEAM),INTENT(INOUT):: B
    INTEGER, OPTIONAL,INTENT(IN):: POS1,POS2,IN_P1,IN_P2
    TYPE(FIBRE), OPTIONAL,POINTER :: P1,P2
    TYPE(THIN_LENS),POINTER, OPTIONAL :: T1,T2
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER I1,I2
    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF
    I1=-2*R%T%N;I2=-2*R%T%N;
    IF(PRESENT(POS1)) I1=POS1
    IF(PRESENT(POS2)) I2=POS2
    IF(PRESENT(T1))   I1=T1%POS
    IF(PRESENT(T2))   I2=T2%POS
    IF(PRESENT(P1))   I1=P1%T1%POS
    IF(PRESENT(P2))   I2=P2%T1%POS
    IF(PRESENT(IN_P1))I1=I1-1+IN_P1
    IF(PRESENT(IN_P2))I2=I2-1+IN_P2
    IF(I1<=0) THEN
       WRITE(6,*) " I1 AND I2 =   ",I1,I2
       WRITE(6,*) " CHECK INPUT"
       RETURN
    ENDIF

    IF(I2>0) THEN
       IF(I2<I1) I2=I2+R%T%N
       CALL TRACK_THIN_LAYOUT_S( R%T,B,I1,I2,K )
    ELSE
       CALL TRACK_THIN_LAYOUT_S( R%T,B,I1,K )
    ENDIF
    B%TIME_INSTEAD_OF_S=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_12

  SUBROUTINE TRACK_LAYOUT_S12( R,B,K,S1,S2 )
    ! Tracks through the thin lens structure from position S1 to position S2 (defined as the S(3) variables
    ! of the thin lens.
    ! The final position is stored as in time tracking but is obviously the same for all the particles.
    ! interfaced as TRACK_LAYOUT_USING_THIN_S

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP),INTENT(IN) :: S1
    REAL(DP),OPTIONAL :: S2
    REAL(DP) DS1,DS2,SFINAL
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER I1,I2,NTURN,I
    TYPE(THIN_LENS), POINTER :: LAST
    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF
    DS2=ZERO
    LAST=>R%T%LAST
    CALL move_to_s( R%T,s1,LAST,i1,ds1 )
    NTURN=0

    I2=-2*R%T%N;
    IF(I1<=0) THEN
       WRITE(6,*) " I1  =   ",I1
       WRITE(6,*) " CHECK INPUT FOR S1"
       RETURN
    ENDIF

    IF(DS1/=ZERO) THEN
       CALL DRIFT_BEAM_BACK_TO_POSITION(LAST,DS1,B)
    ENDIF

    IF(PRESENT(S2)) THEN
       CALL move_to_s( R%T,s2,LAST,i2,ds2 )
       IF(R%CLOSED) THEN
          NTURN=INT((S2-DS2-S1+DS1)/R%T%END%S(3))
       ENDIF
       SFINAL=S2
    ELSE
       SFINAL=S1
    ENDIF

    DO I=1,NTURN
       CALL TRACK_THIN_LAYOUT_S( R%T,B,I1,K )
    ENDDO

    IF(I2>0) THEN
       if(i2<i1) then
          i2=r%t%n+i2
       endif
       CALL TRACK_THIN_LAYOUT_S( R%T,B,I1,I2,K )
    ELSE
       CALL TRACK_THIN_LAYOUT_S( R%T,B,I1,K )
    ENDIF

    IF(DS2/=ZERO) THEN
       DS2=-DS2
       CALL DRIFT_BEAM_BACK_TO_POSITION(LAST,DS2,B)
    ENDIF
    B%X(:,7:7)=SFINAL
    B%TIME_INSTEAD_OF_S=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_S12



  SUBROUTINE TRACK_THIN_LAYOUT_12( R,B,I1,I2,K ) ! T
    implicit none
    TYPE (THIN_LAYOUT), TARGET :: R
    TYPE (THIN_LENS), POINTER ::  C
    TYPE(BEAM),INTENT(INOUT):: B
    INTEGER, INTENT(IN):: I1,I2
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER I,J

    call move_to_THIN_LENS(r,c,I1)


    if(i2>i1) then
       J=I1

       DO  WHILE(J<I2.AND.ASSOCIATED(C))

          CALL TRACK_THIN_SINGLE_FOR_S(C,B,K)

          C=>C%NEXT
          J=J+1
       ENDDO
    ELSEIF(I1>I2) THEN
       WRITE(6,*) " BACKWARDS FORBIDDEN IN TRACK_THIN_LAYOUT"
       STOP 666
    ENDIF
    B%TIME_INSTEAD_OF_S=.FALSE.
  end SUBROUTINE TRACK_THIN_LAYOUT_12

  SUBROUTINE TRACK_THIN_LAYOUT_1( R,B,I1,K ) ! TRACKS LAYOUT WITHOUT ANY COLLECTIVE FRIVOLITES
    implicit none
    TYPE (THIN_LAYOUT), TARGET :: R
    TYPE(BEAM),INTENT(INOUT):: B
    INTEGER, INTENT(IN):: I1
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER II1,II2

    II1=I1
    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK_THIN_LAYOUT_s(R,B,II1,II2,k)

    B%TIME_INSTEAD_OF_S=.FALSE.

  end SUBROUTINE TRACK_THIN_LAYOUT_1

  !  STUFF ABOUT LIST AND STRUCTURES

  SUBROUTINE MAKE_THIN_LAYOUT( R) !
    ! Creates the thin layout and puts in R%T by calling MAKE_THIN_LAYOUT_2
    implicit none
    TYPE (LAYOUT), TARGET :: R

    call MAKE_THIN_LAYOUT_2( R,R%T )

  end SUBROUTINE MAKE_THIN_LAYOUT



  SUBROUTINE MAKE_THIN_LAYOUT_2( R,L ) !
    ! Creates a thin layout.
    ! At this point large patches would be problematic.

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE (THIN_LAYOUT), pointer :: L
    TYPE(FIBRE), POINTER :: P
    INTEGER I,J,k,TEAPOT_LIKE
    REAL(DP) S,DLD,DL,LI,SL
    LOGICAL(LP) CIRCULAR
    TYPE(THIN_LENS), POINTER :: T1,T2

    CASE_NAME(CASEP1)="THE ENTRANCE PATCH"
    CASE_NAME(CASEP2)="THE EXIT PATCH  "
    CASE_NAME(CASE1)= "THE ENTRANCE FRINGE"
    CASE_NAME(CASE2)="THE EXIT FRINGE"
    CASE_NAME(CASE0)="A STEP OF INTEGRATION BODY"

    if(associated(L)) then
       CALL kill_THIN_Layout(L)
       DEALLOCATE(L);
       NULLIFY(L);
    endif

    allocate(L)
    CALL Set_Up_THIN_LAYOUT( L )
    S=zero
    SL=ZERO  !  INTEGRATION LENGTH
    P=>R%START
    k=1
    DO I=1,R%N

       TEAPOT_LIKE=0
       IF(P%MAG%P%B0/=ZERO) TEAPOT_LIKE=1
       IF(P%MAG%KIND==KIND16.OR.P%MAG%KIND==KIND16)TEAPOT_LIKE=0

       IF(P%DIR==1) THEN
          LI=zero;
       ELSE
          LI=P%MAG%L;
       ENDIF
       DL=P%DIR*P%MAG%L/P%MAG%P%NST
       DLD=P%MAG%P%LD/P%MAG%P%NST
       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       T1=>L%END

       L%END%CAS=CASEP1
       L%END%pos_in_fibre=1
       L%END%pos=k;k=k+1;
       L%END%PARENT_THIN_LAYOUT=>L
       L%END%PARENT_FIBRE=>P

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       L%END%CAS=CASE1
       L%END%pos_in_fibre=2
       L%END%pos=k;k=k+1;
       L%END%PARENT_THIN_LAYOUT=>L
       L%END%PARENT_FIBRE=>P

       DO J=1,P%MAG%P%NST
          CALL APPEND_EMPTY_THIN( L )
          L%END%TEAPOT_LIKE=TEAPOT_LIKE
          L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=DL;
          L%END%CAS=CASE0
          L%END%pos_in_fibre=J+2
          L%END%pos=k;k=k+1;
          L%END%PARENT_THIN_LAYOUT=>L
          L%END%PARENT_FIBRE=>P
          S=S+DLD
          LI=LI+DL
          SL=SL+DL
       ENDDO

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       L%END%CAS=CASE2
       L%END%pos_in_fibre=P%MAG%P%NST+3
       L%END%pos=k;k=k+1;
       L%END%PARENT_FIBRE=>P

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       L%END%CAS=CASEP2
       L%END%pos_in_fibre=P%MAG%P%NST+4
       L%END%pos=k;k=k+1;
       L%END%PARENT_THIN_LAYOUT=>L
       L%END%PARENT_FIBRE=>P
       T2=>L%END

       P%T1=>T1
       P%T2=>T2

       P=>P%NEXT
    ENDDO
    L%N=k-1

    L%PARENT_LAYOUT=>R

    IF(R%CLOSED) THEN
       l%closed=.true.
       CIRCULAR=.TRUE.
       CALL RING_L_THIN(L,CIRCULAR)
    ENDIF

    call stat_THIN_LAYOUT(l)

  END SUBROUTINE MAKE_THIN_LAYOUT_2


  SUBROUTINE stat_THIN_LAYOUT( L )
    implicit none
    TYPE (THIN_LAYOUT), pointer :: L

    WRITE(6,*)  " PARENT LAYOUT NAME :", L%PARENT_LAYOUT%NAME(1:len_trim(L%PARENT_LAYOUT%NAME))
    WRITE(6,*) " NUMBER OF ORIGINAL LAYOUT ELEMENTS :", L%PARENT_LAYOUT%N
    WRITE(6,*) " NUMBER OF THIN OBJECTS :", L%N
    WRITE(6,*) " TOTAL IDEAL LENGTH OF STRUCTURE :", L%END%S(1)
    WRITE(6,*) " TOTAL INTEGERATION LENGTH OF STRUCTURE (mad8 style survey) :", L%END%S(3)

  end SUBROUTINE stat_THIN_LAYOUT


  SUBROUTINE DRIFT_TO_TIME(T,YL,DT,X)
    ! Drifts to a given time using either a regular or TEAPOT drift. (Cylindrical coordinate drift)
    ! The amount of the drift YL is computed to achieve a time DT.

    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(INOUT):: YL
    real(dp),INTENT(IN):: DT
    TYPE(THIN_LENS), pointer :: T
    TYPE(magnet_chart), pointer :: p
    real(dp) XN(6),PZ,PT
    real(dp)  A,b,R

    P=>T%PARENT_FIBRE%MAG%P

    IF(P%TIME) then
       B=P%BETA0
    ELSE
       B=one
    ENDIF


    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then


       PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
       R=one/P%B0
       A=(X(1)+R)*(one/b+x(5))/PZ
       A=ATAN(DT/(A+DT*X(2)/PZ) )
       YL=A/P%B0
       PT=one-X(2)*TAN(A)/PZ
       XN(1)=(X(1)+R)/COS(A)/PT-R
       XN(2)=X(2)*COS(A)+SIN(A)*PZ
       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    else


       PZ=ROOT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)

       YL=DT/(one/b+X(5))*PZ

       X(1)=X(1)+YL*X(2)/PZ
       X(3)=X(3)+YL*X(4)/PZ
       X(6)=X(6)+YL*(one/b+X(5))/PZ



    endif



  END SUBROUTINE DRIFT_TO_TIME

  SUBROUTINE DRIFT_BACK_TO_POSITION(T,YL,X)
    ! This is a regular drift
    ! It is used in time tracking to project back to the beginning of the thin lens
    ! and it is used in S tracking to drift in the middle of an step.

    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: YL
    TYPE(THIN_LENS), pointer :: T
    TYPE(magnet_chart), pointer :: p
    real(dp) XN(6),PZ,PT
    real(dp)  A,b,R

    P=>T%PARENT_FIBRE%MAG%P

    IF(P%TIME) then
       B=P%BETA0
    ELSE
       B=one
    ENDIF


    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then



       PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
       R=one/P%B0
       A=-YL*P%B0

       PT=one-X(2)*TAN(A)/PZ
       XN(1)=(X(1)+R)/COS(A)/PT-R
       XN(2)=X(2)*COS(A)+SIN(A)*PZ
       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
       !          WRITE(6,*) "XN(6)-X(6)-DT , DT",DT,XN(6)-X(6)-DT
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    else
       !       CALL DRIFT(YL,DL,P%beta0,1,P%EXACT,P%TIME,X)
       PZ=ROOT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)


       X(1)=X(1)-YL*X(2)/PZ
       X(3)=X(3)-YL*X(4)/PZ
       X(6)=X(6)-YL*(one/b+X(5))/PZ

    endif



  END SUBROUTINE DRIFT_BACK_TO_POSITION

  SUBROUTINE DRIFT_BEAM_BACK_TO_POSITION(T,YL,B)
    ! Same routine as above but applies to a beam
    IMPLICIT NONE
    real(dp)  X(6)
    TYPE(BEAM), INTENT(INOUT) :: B
    real(dp),INTENT(IN):: YL
    TYPE(THIN_LENS), pointer :: T
    INTEGER I

    DO I=1,B%N
       IF(B%U(i)) CYCLE
       X=BEAM_IN_X(B,I)

       CALL DRIFT_BACK_TO_POSITION(T,YL,X)

       B%POS(I)%THINLENS=>T

       CALL X_IN_BEAM(X,B,I,DL=ZERO)


    ENDDO




  END SUBROUTINE DRIFT_BEAM_BACK_TO_POSITION

  ! BEAM STUFF

  subroutine create_beam(B,N,CUT,SIG,T)
    USE gauss_dis
    implicit none
    INTEGER N,I,J
    REAL(DP) CUT,SIG,X
    TYPE(BEAM) B
    TYPE (THIN_LENS),optional,target::  T

    IF(.NOT.ASSOCIATED(B%N)) THEN
       CALL ALLOCATE_BEAM(B,N)
    ELSEIF(B%N/=N) THEN
       CALL KILL_BEAM(B)
       CALL ALLOCATE_BEAM(B,N)
    ENDIF

    DO I=1,N
       DO J=1,6
          CALL GRNF(X,cut)
          B%X(I,J)=X*SIG
       ENDDO
    ENDDO


    if(present(t)) then
       DO I=1,N
          if(associated(B%POS(I)%THINLENS))then
             B%POS(I)%THINLENS=>T
          endif
       ENDDO

    endif

  end    subroutine create_beam

  subroutine create_PANCAKE(B,N,CUT,SIG,T,A)
    USE gauss_dis
    implicit none
    INTEGER N,I,J
    REAL(DP) CUT,SIG(6),X,Y(LNV)
    TYPE(BEAM) B
    TYPE (THIN_LENS),optional,target::  T
    TYPE (DAMAP),OPTIONAL :: A

    IF(.NOT.ASSOCIATED(B%N)) THEN
       CALL ALLOCATE_BEAM(B,N)
    ELSEIF(B%N/=N) THEN
       CALL KILL_BEAM(B)
       CALL ALLOCATE_BEAM(B,N)
    ENDIF

    Y=ZERO
    DO I=1,N
       IF(.not.PRESENT(A)) THEN

          DO J=1,6
             CALL GRNF(X,cut)
             B%X(I,J)=X*SIG(J)
          ENDDO

       ELSE

          DO J=1,C_%ND2
             CALL GRNF(X,cut)
             Y(I)=X*SIG(I)
          ENDDO

          B%X(I,1:C_%ND2)=A*Y

          DO J=C_%ND2+1,6
             CALL GRNF(X,cut)
             B%X(I,J)=X*SIG(J)
          ENDDO


       ENDIF

       B%X(I,7)=ZERO
    ENDDO



    if(present(t)) then
       DO I=1,N
          if(associated(B%POS(I)%THINLENS))then
             B%POS(I)%THINLENS=>T
          endif
       ENDDO

    endif
  end    subroutine create_PANCAKE

  subroutine copy_beam(B1,B2)
    implicit none
    INTEGER I
    TYPE(BEAM), INTENT(INOUT) :: B1,B2

    IF(.NOT.ASSOCIATED(B2%N)) THEN
       CALL ALLOCATE_BEAM(B2,B1%N)
    ELSEIF(B1%N/=B2%N) THEN
       CALL KILL_BEAM(B2)
       CALL ALLOCATE_BEAM(B2,B1%N)
    ENDIF

    B2%X=B1%X
    B2%U=B1%U
    B2%N=B1%N
    B2%CHARGE=B1%CHARGE
    B2%LOST=B1%LOST
    DO I=1,B1%N
       if(associated(B1%POS(I)%THINLENS))then
          B2%POS(I)%THINLENS=>B1%POS(I)%THINLENS
       endif
    ENDDO

  END subroutine copy_beam

  subroutine PRINT_beam(B,MF,I)
    implicit none
    INTEGER K,MF,I1,I2
    INTEGER,OPTIONAL:: I
    TYPE(BEAM), INTENT(IN):: B
    TYPE(THIN_LENS),POINTER::T
    TYPE(FIBRE),POINTER::F

    I1=1
    I2=B%N

    IF(PRESENT(I)) THEN
       I1=I
       I2=I
    ENDIF
    IF(B%TIME_INSTEAD_OF_S) THEN
       WRITE(MF,*) "____________________________ TIME TRACKED BEAM __________________________________"
    ELSE
       WRITE(MF,*) "_________________ POSITION TRACKED BEAM (AS IN PTC PROPER)_______________________"
    ENDIF

    DO K=I1,I2
       IF(B%U(K)) THEN
          WRITE(MF,*) " PARTICLE # ",K, " IS LOST "
       ELSE
          T=>B%POS(K)%THINLENS
          F=>T%PARENT_FIBRE
          WRITE(MF,*) "_________________________________________________________________________"
          WRITE(MF,*) " PARTICLE # ",K, " IS LOCATED AT SLICE # ",T%POS," IN FIBRE  ",F%MAG%NAME
          WRITE(MF,*) " IN THE FIBRE POSITION  ",T%pos_in_fibre
          WRITE(MF,*) " IN ",CASE_NAME(T%CAS)
          IF(T%CAS==CASE0)WRITE(MF,*) " AT THE STEP NUMBER ",T%pos_in_fibre-2

          WRITE(MF,*) "........................................................................."
          IF(B%TIME_INSTEAD_OF_S) THEN
             WRITE(MF,*) " TIME AND POSITION AFTER THIN SLICE = ",B%X(K,6:7)
          ELSE
             WRITE(MF,*) " TIME AND POSITION  = ",B%X(K,6:7)
          ENDIF
          WRITE(MF,*) " X,Y = ",B%X(K,1),B%X(K,3)
          WRITE(MF,*) " PX,PY = ",B%X(K,2),B%X(K,4)
          WRITE(MF,*) " ENERGY VARIABLE = ",B%X(K,5)
       ENDIF
       WRITE(MF,*) "_________________________________________________________________________"
    ENDDO

  END subroutine PRINT_beam



  SUBROUTINE NULLIFY_BEAM(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    NULLIFY(B%N,B%LOST)
    NULLIFY(B%X)
    NULLIFY(B%U)
    NULLIFY(B%POS)
    NULLIFY(B%CHARGE)
    NULLIFY(B%TIME_INSTEAD_OF_S)
  END SUBROUTINE NULLIFY_BEAM

  SUBROUTINE NULLIFY_BEAMS(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B(:)
    INTEGER I
    DO I=1,SIZE(B)
       CALL NULLIFY_BEAM(B(i))
    ENDDO

  END SUBROUTINE NULLIFY_BEAMS

  SUBROUTINE ALLOCATE_BEAM(B,N)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    INTEGER , INTENT (IN) :: N
    INTEGER I
    ALLOCATE(B%N,B%LOST)
    B%N=N
    B%LOST=0
    ALLOCATE(B%X(N,7))
    ALLOCATE(B%U(N))
    ALLOCATE(B%POS(N))
    DO I=1,N
       NULLIFY(B%POS(i)%THINLENS)
    ENDDO
    ALLOCATE(B%CHARGE)
    ALLOCATE(B%TIME_INSTEAD_OF_S)

    B%X  = ZERO
    B%U  = .FALSE.
    B%CHARGE=1
    B%TIME_INSTEAD_OF_S=.FALSE.

  END SUBROUTINE ALLOCATE_BEAM

  SUBROUTINE KILL_BEAM(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    IF(ASSOCIATED(B%N)) DEALLOCATE(B%N,B%LOST,B%X,B%U,B%POS,B%CHARGE,B%TIME_INSTEAD_OF_S)
  END SUBROUTINE KILL_BEAM

  SUBROUTINE KILL_BEAMS(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B(:)
    INTEGER I
    DO I=1,SIZE(B)
       CALL KILL_BEAM(B(i))
    ENDDO
  END SUBROUTINE KILL_BEAMS


  FUNCTION BEAM_IN_X(B,I)
    IMPLICIT NONE
    REAL(DP) BEAM_IN_X(6)
    TYPE(BEAM), INTENT(INOUT) ::B
    INTEGER, INTENT(IN) :: I

    BEAM_IN_X=B%X(I,1:6)

  END  FUNCTION BEAM_IN_X

  SUBROUTINE X_IN_BEAM(X,B,I,DL)
    IMPLICIT NONE
    REAL(DP) X(6)
    REAL(DP),OPTIONAL:: DL
    TYPE(BEAM), INTENT(INOUT) ::B
    TYPE(THIN_LENS),POINTER :: T
    INTEGER, INTENT(IN) :: I

    B%X(I,1:6)=X(1:6)
    IF(PRESENT(DL)) B%X(I,7)=DL
    if(.not.CHECK_STABLE) then
       write(6,*) "unstable "
       CALL RESET_APERTURE_FLAG
       b%u(I)=.true.
       B%LOST=B%LOST+1
    endif

  END  SUBROUTINE X_IN_BEAM


end module ptc_multiparticle
