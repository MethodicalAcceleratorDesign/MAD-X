!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
MODULE S_TRACKING
  USE S_FAMILY

  IMPLICIT NONE
  logical(lp),TARGET :: ALWAYS_EXACT_PATCHING=.TRUE.
  type(fibre), pointer :: lost_fibre

  ! linked
  PRIVATE TRACK_LAYOUT_FLAG_R,TRACK_LAYOUT_FLAG_P,TRACK_LAYOUT_FLAG_S
  !  PRIVATE FIND_ORBIT_LAYOUT,FIND_ORBIT_M_LAYOUT,FIND_ENV_LAYOUT, FIND_ORBIT_LAYOUT_noda
  PRIVATE TRACK_LAYOUT_FLAG_R1,TRACK_LAYOUT_FLAG_P1,TRACK_LAYOUT_FLAG_S1
  PRIVATE MIS_FIBR,MIS_FIBP,MIS_FIBS
  PRIVATE TRACK_FIBRE_R,TRACK_FIBRE_P,TRACK_FIBRE_S
  private TRACK_LAYOUT_POINT_R1
  TYPE UPDATING
     logical(lp) UPDATE
  END TYPE UPDATING
  TYPE (UPDATING), PARAMETER ::  COMPUTE= UPDATING(.TRUE.)



  INTERFACE TRACK
     ! linked
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_R
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_P
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_S
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_R1
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_P1
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_S1
     MODULE PROCEDURE TRACK_FIBRE_R
     MODULE PROCEDURE TRACK_FIBRE_P
     MODULE PROCEDURE TRACK_FIBRE_S
     MODULE PROCEDURE TRACK_LAYOUT_POINT_R1  ! track one fibre
  END INTERFACE

  !  INTERFACE FIND_ORBIT
  !     ! LINKED
  !     ! no use of TPSA
  !     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
  !     ! returns linear matrix
  !     MODULE PROCEDURE FIND_ORBIT_M_LAYOUT
  !     !RETURN QUADRATIC YS
  !     MODULE PROCEDURE FIND_ENV_LAYOUT
  !
  !  END INTERFACE

  INTERFACE MIS_FIB
     MODULE PROCEDURE MIS_FIBR
     MODULE PROCEDURE MIS_FIBP
     MODULE PROCEDURE MIS_FIBS
  END INTERFACE


contains

  ! linked stuff

  SUBROUTINE TRACK_LAYOUT_POINT_R1(R,X,pIN,k) ! Tracks real(dp) from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE(FIBRE), POINTER:: pIN,P
    INTEGER II1,II2

    P=>R%START
    DO II1=1,R%N
       IF(ASSOCIATED(PIN,P)) EXIT
       P=>P%NEXT
    ENDDO

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k)
  END SUBROUTINE TRACK_LAYOUT_POINT_R1


  SUBROUTINE TRACK_LAYOUT_FLAG_R1(R,X,II1,k) ! Tracks real(dp) from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2


    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k)
  END SUBROUTINE TRACK_LAYOUT_FLAG_R1

  SUBROUTINE TRACK_LAYOUT_FLAG_P1(R,X,II1,k) ! Tracks polymorphs from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2

    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k)
  END SUBROUTINE TRACK_LAYOUT_FLAG_P1

  SUBROUTINE TRACK_LAYOUT_FLAG_S1(R,X,II1,k) ! Tracks envelope from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2


    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k)
  END SUBROUTINE TRACK_LAYOUT_FLAG_S1



  SUBROUTINE TRACK_LAYOUT_FLAG_R(R,X,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J
    TYPE (fibre), POINTER :: C


    CALL RESET_APERTURE_FLAG

    call move_to(r,c,MOD_N(I1,R%N))


    J=I1

    DO  WHILE(J<I2.AND.ASSOCIATED(C))

       CALL TRACK(C,X,K,R%CHARGE)

       C=>C%NEXT
       J=J+1
    ENDDO


  END SUBROUTINE TRACK_LAYOUT_FLAG_R



  SUBROUTINE TRACK_LAYOUT_FLAG_P(R,X,I1,I2,K) ! TRACKS POLYMORPHS FROM I1 TO I2 IN STATE K
    IMPLICIT NONE
    TYPE(LAYOUT),INTENT(INOUT):: R ;TYPE(REAL_8), INTENT(INOUT):: X(6);
    INTEGER, INTENT(IN):: I1,I2; TYPE(INTERNAL_STATE) K;
    INTEGER J;

    TYPE (FIBRE), POINTER :: C


    CALL RESET_APERTURE_FLAG

    CALL MOVE_TO(R,C,MOD_N(I1,R%N))


    J=I1

    DO  WHILE(J<I2.AND.ASSOCIATED(C))

       CALL TRACK(C,X,K,R%CHARGE)

       C=>C%NEXT
       J=J+1
    ENDDO


    ! PATCHES
  END SUBROUTINE TRACK_LAYOUT_FLAG_P

  SUBROUTINE TRACK_LAYOUT_FLAG_S(R,X,I1,I2,k) ! Tracks envelopes from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER I,J,M,N
    ! YS SPECIFIC STUFF
    TYPE(DAMAP) ID,XT,DISP
    TYPE(REAL_8) XR(6),X1,X3
    real(dp) V(6)
    TYPE (fibre), POINTER :: C



    CALL RESET_APERTURE_FLAG

    ! new stuff with kind=3
    IF(k%para_in ) knob=.true.
    ! end new stuff with kind=3


    call move_to(r,c,MOD_N(I1,R%N))


    J=I1

    DO  WHILE(J<I2.AND.ASSOCIATED(C))
       CALL TRACK(C,X,K,R%CHARGE)

       C=>C%NEXT
       J=J+1
    ENDDO
    ! new stuff with kind=3
    knob=.FALSE.
    ! end new stuff with kind=3

    ! Radiation

    CALL ALLOC(ID);CALL ALLOC(XT);CALL ALLOC(DISP);
    CALL ALLOC(X1,X3);CALL ALLOC(XR);
    xr=x
    v=xr
    id=0
    XT=XR
    DISP=XT*ID
    ID=1
    DISP=ID-DISP
    XT=DISP*XT

    xr=xt+v


    do i=1,6
       do j=1,6
          X(I)%SIGMAf(J)=zero
       enddo
    enddo

    do i=1,6
       do j=1,6
          DO M=1,6
             DO N=1,6
                X1=(xr(I)).par.ind_stoc(M)
                X3=(xr(j)).par.ind_stoc(n)
                X(I)%SIGMAf(J)=X(m)%E(n)*x1*x3+X(I)%SIGMAf(J)
                X(I)%SIGMAf(J)=X(m)%SIGMA0(n)*x1*x3+X(I)%SIGMAf(J)
             enddo
          enddo
       enddo
    enddo

    CALL KILL(XT);CALL KILL(DISP);
    CALL KILL(X1,X3);CALL KILL(XR);




  END SUBROUTINE TRACK_LAYOUT_FLAG_S

  SUBROUTINE TRACK_FIBRE_R(C,X,K,CHARGE)
    implicit none
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    INTEGER, target, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch,PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0

    IF(.NOT.CHECK_STABLE) return

    ! DIRECTIONAL VARIABLE
    C%MAG%P%DIR=>C%DIR
    C%MAG%P%CHARGE=>CHARGE
    !
    !    IF(.NOT.CHECK_STABLE) CHECK_STABLE=.TRUE.

    C%MAG=K
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=.FALSE. ; PATCHE=.FALSE. ;PATCHG=.FALSE.;
    ENDIF
    IF(PATCHE) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAG%P%P0C
             B0=>CN%MAG%P%BETA0

             X(2)=X(2)*P0/C%MAG%P%P0C
             X(4)=X(4)*P0/C%MAG%P%P0C
             IF(C%MAG%P%TIME)THEN
                X(5)=SQRT(one+two*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-one !X(5) = DP/P0C_NEW
                X(5)=(two*X(5)+X(5)**2)/(SQRT(one/C%MAG%P%BETA0**2+two*X(5)+X(5)**2)+one/C%MAG%P%BETA0)
             ELSE
                X(5)=(one+X(5))*P0/C%MAG%P%P0C-one
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

    ENDIF

    IF(PATCHG) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN
          X(1)=CN%DIR*C%DIR*X(1);X(2)=CN%DIR*C%DIR*X(2);
       ENDIF
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
    ENDIF
    IF(PATCHT) THEN
       X(6)=X(6)-C%PATCH%a_T
    ENDIF

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       ou = K%EXACTMIS.or.C%MAG%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF

    CALL TRACK(C%MAG,X)

    IF(C%MAG%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF

    IF(PATCHT) THEN
       X(6)=X(6)-C%PATCH%b_T
    ENDIF

    IF(PATCHG) THEN
       CN=>C%NEXT
       IF(ASSOCIATED(CN)) THEN
          X(1)=CN%DIR*C%DIR*X(1);X(2)=CN%DIR*C%DIR*X(2);
       ENDIF
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X,DONEITT)
       CALL TRANS(C%PATCH%B_D,X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
    ENDIF

    IF(PATCHE) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAG%P%P0C
       B0=>CN%MAG%P%BETA0
       X(2)=X(2)*C%MAG%P%P0C/P0
       X(4)=X(4)*C%MAG%P%P0C/P0
       IF(C%MAG%P%TIME)THEN
          X(5)=SQRT(one+two*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAG%P%P0C/P0-one !X(5) = DP/P0C_NEW
          X(5)=(two*X(5)+X(5)**2)/(SQRT(one/B0**2+two*X(5)+X(5)**2)+one/B0)
       ELSE
          X(5)=(one+X(5))*C%MAG%P%P0C/P0-one
       ENDIF
    ENDIF

    C%MAG=DEFAULT
    nullify(C%MAG%P%DIR)
    nullify(C%MAG%P%CHARGE)
    if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
       if(CHECK_MADX_APERTURE) c_%message="exceed absolute_aperture in TRACK_FIBRE_R"
       CHECK_STABLE=.false.
       lost_fibre=>c
    endif

  END SUBROUTINE TRACK_FIBRE_R

  SUBROUTINE TRACK_FIBRE_P(C,X,K,CHARGE)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    INTEGER, TARGET, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) OU,PATCH,PATCHT,PATCHG,PATCHE
    TYPE (FIBRE), POINTER :: CN
    REAL(DP), POINTER :: P0,B0

    IF(.NOT.CHECK_STABLE) return

    ! NEW STUFF WITH KIND=3: KNOB OF FPP IS SET TO TRUE IF NECESSARY
    IF(K%PARA_IN ) KNOB=.TRUE.
    ! END NEW STUFF WITH KIND=3

    ! DIRECTIONAL VARIABLE AND CHARGE IS PASSED TO THE ELEMENT
    C%MAGP%P%DIR=>C%DIR
    C%MAGP%P%CHARGE=>CHARGE
    !

    ! PASSING THE STATE K TO THE ELEMENT
    C%MAGP=K
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=.FALSE. ; PATCHE=.FALSE. ;PATCHG=.FALSE.;
    ENDIF
    ! ENERGY PATCH
    IF(PATCHE) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          IF(.NOT.CN%PATCH%ENERGY) THEN     ! NO NEED TO PATCH IF PATCHED BEFORE
             P0=>CN%MAGP%P%P0C
             B0=>CN%MAGP%P%BETA0

             X(2)=X(2)*P0/C%MAGP%P%P0C
             X(4)=X(4)*P0/C%MAGP%P%P0C
             IF(C%MAGP%P%TIME)THEN
                X(5)=SQRT(ONE+TWO*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAGP%P%P0C-ONE !X(5) = DP/P0C_NEW
                X(5)=(TWO*X(5)+X(5)**2)/(SQRT(ONE/C%MAGP%P%BETA0**2+TWO*X(5)+X(5)**2)+ONE/C%MAGP%P%BETA0)
             ELSE
                X(5)=(ONE+X(5))*P0/C%MAGP%P%P0C-ONE
             ENDIF
          ENDIF ! NO NEED TO PATCH
       ENDIF ! ASSOCIATED

    ENDIF


    ! POSITION PATCH
    IF(PATCHG) THEN
       PATCH=ALWAYS_EXACT_PATCHING.OR.C%MAGP%P%EXACT
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN
          X(1)=CN%DIR*C%DIR*X(1);X(2)=CN%DIR*C%DIR*X(2);
       ENDIF
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
    ENDIF
    ! TIME PATCH
    IF(PATCHT) THEN
       X(6)=X(6)-C%PATCH%A_T
    ENDIF

    ! MISALIGNMENTS AT THE ENTRANCE
    IF(C%MAGP%MIS) THEN
       OU = K%EXACTMIS.OR.C%MAGP%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF

    ! ************************************************************************
    !  THE ACTUAL MAGNET PROPAGATOR AS IT WOULD APPEAR IN A STANDARD CODE

    CALL TRACK(C%MAGP,X)
    !
    ! ************************************************************************

    ! MISALIGNMENTS AT THE EXIT
    IF(C%MAGP%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF

    !EXIT PATCH
    ! TIME PATCH
    IF(PATCHT) THEN
       X(6)=X(6)-C%PATCH%B_T
    ENDIF

    ! POSITION PATCH
    IF(PATCHG) THEN
       CN=>C%NEXT
       IF(ASSOCIATED(CN)) THEN
          X(1)=CN%DIR*C%DIR*X(1);X(2)=CN%DIR*C%DIR*X(2);
       ENDIF
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X,DONEITT)
       CALL TRANS(C%PATCH%B_D,X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
    ENDIF

    ! ENERGY PATCH
    IF(PATCHE) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAGP%P%P0C
       B0=>CN%MAGP%P%BETA0
       X(2)=X(2)*C%MAGP%P%P0C/P0
       X(4)=X(4)*C%MAGP%P%P0C/P0
       IF(C%MAGP%P%TIME)THEN
          X(5)=SQRT(ONE+TWO*X(5)/C%MAGP%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAGP%P%P0C/P0-ONE !X(5) = DP/P0C_NEW
          X(5)=(TWO*X(5)+X(5)**2)/(SQRT(ONE/B0**2+TWO*X(5)+X(5)**2)+ONE/B0)
       ELSE
          X(5)=(ONE+X(5))*C%MAGP%P%P0C/P0-ONE
       ENDIF
    ENDIF


    ! ELEMENT IS RESTAURED TO THE DEFAULT STATE
    C%MAGP=DEFAULT
    ! DIRECTIONAL VARIABLE AND CHARGE ARE ELIMINATED
    NULLIFY(C%MAGP%P%DIR)
    NULLIFY(C%MAGP%P%CHARGE)


    ! KNOB IS RETURNED TO THE PTC DEFAULT
    ! NEW STUFF WITH KIND=3
    KNOB=.FALSE.
    ! END NEW STUFF WITH KIND=3

    if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
       if(CHECK_MADX_APERTURE) c_%message="exceed absolute_aperture in TRACK_FIBRE_P"
       CHECK_STABLE=.false.
       lost_fibre=>c
    endif

  END SUBROUTINE TRACK_FIBRE_P


  SUBROUTINE TRACK_FIBRE_S(C,X,K,CHARGE,UPDATE)
    implicit none
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(REAL_8)  Y(6)
    INTEGER, target, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch,PATCHT,PATCHG,PATCHE
    TYPE(UPDATING), optional,intent(in):: UPDATE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    ! YS SPECIFIC STUFF
    TYPE(DAMAP) ID,XT,DISP
    TYPE(REAL_8) XR(6),X1,X3
    real(dp) V(6)
    INTEGER I,J,M,N

    IF(.NOT.CHECK_STABLE) return

    ! new stuff with kind=3
    IF(k%para_in ) knob=.true.
    ! end new stuff with kind=3

    ! DIRECTIONAL VARIABLE
    C%MAGP%P%DIR=>C%DIR
    C%MAGP%P%CHARGE=>CHARGE
    !

    C%MAGP=K
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=.FALSE. ; PATCHE=.FALSE. ;PATCHG=.FALSE.;
    ENDIF


    IF(PATCHE) THEN
       CALL ALLOC(Y)
       Y=X
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          IF(.NOT.CN%PATCH%ENERGY) THEN   ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAG%P%P0C
             B0=>CN%MAG%P%BETA0

             Y(2)=Y(2)*P0/C%MAGP%P%P0C
             Y(4)=Y(4)*P0/C%MAGP%P%P0C
             IF(C%MAGP%P%TIME)THEN
                Y(5)=SQRT(one+two*Y(5)/B0+Y(5)**2)  !Y(5) = 1+DP/P0C_OLD
                Y(5)=Y(5)*P0/C%MAGP%P%P0C-one !Y(5) = DP/P0C_NEW
                Y(5)=(two*Y(5)+Y(5)**2)/(SQRT(one/C%MAGP%P%BETA0**2+two*Y(5)+Y(5)**2)+one/C%MAGP%P%BETA0)
             ELSE
                Y(5)=(one+Y(5))*P0/C%MAGP%P%P0C-one
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

       X=Y
       CALL KILL(Y)
    ENDIF





    IF(PATCHG) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN
          X(1)%V=CN%DIR*C%DIR*X(1)%V;X(2)%V=CN%DIR*C%DIR*X(2)%V;
       ENDIF
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
    ENDIF
    IF(PATCHE) THEN
       CALL ALLOC(Y)
       Y=X
       Y(6)=Y(6)-C%PATCH%a_T
       X=Y
       CALL KILL(Y)
    ENDIF



    !      call TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAGP%MIS) THEN
       ou = K%EXACTMIS.or.C%MAGP%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF

    CALL TRACK(C%MAGP,X)


    IF(C%MAGP%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF

    IF(PATCHE) THEN
       CALL ALLOC(Y)
       Y=X
       Y(6)=Y(6)-C%PATCH%b_T
       X=Y
       CALL KILL(Y)
    ENDIF
    IF(PATCHG) THEN
       CN=>C%NEXT
       IF(ASSOCIATED(CN)) THEN
          X(1)%V=CN%DIR*C%DIR*X(1)%V;X(2)%V=CN%DIR*C%DIR*X(2)%V;
       ENDIF
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X,DONEITT)
       CALL TRANS(C%PATCH%B_D,X,C%MAGP%P%BETA0,PATCH,C%MAGP%P%TIME)
    ENDIF

    IF(PATCHE) THEN
       CALL ALLOC(Y)
       Y=X
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAGP%P%P0C
       B0=>CN%MAGP%P%BETA0
       Y(2)=Y(2)*C%MAGP%P%P0C/P0
       Y(4)=Y(4)*C%MAGP%P%P0C/P0
       IF(C%MAGP%P%TIME)THEN
          Y(5)=SQRT(one+two*Y(5)/C%MAGP%P%BETA0+Y(5)**2)  !Y(5) = 1+DP/P0C_OLD
          Y(5)=Y(5)*C%MAGP%P%P0C/P0-one !Y(5) = DP/P0C_NEW
          Y(5)=(two*Y(5)+Y(5)**2)/(SQRT(one/B0**2+two*Y(5)+Y(5)**2)+one/B0)
       ELSE
          Y(5)=(one+Y(5))*C%MAGP%P%P0C/P0-one
       ENDIF
       X=Y
       CALL KILL(Y)
    ENDIF


    nullify(C%MAGP%P%DIR)
    nullify(C%MAGP%P%CHARGE)
    C%MAGP=DEFAULT


    if(abs(x(1)%v)+abs(x(3)%v)>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
       if(CHECK_MADX_APERTURE) c_%message="exceed absolute_aperture in TRACK_FIBRE_S"
       CHECK_STABLE=.false.
       lost_fibre=>c
    endif

    ! new stuff with kind=3
    knob=.FALSE.
    ! end new stuff with kind=3
    if(present(update)) then
       ! Radiation

       CALL ALLOC(ID);CALL ALLOC(XT);CALL ALLOC(DISP);
       CALL ALLOC(X1,X3);CALL ALLOC(XR);
       xr=x
       v=xr
       id=0
       XT=XR
       DISP=XT*ID
       ID=1
       DISP=ID-DISP
       XT=DISP*XT

       xr=xt+v


       do i=1,6
          do j=1,6
             X(I)%SIGMAf(J)=zero
          enddo
       enddo

       do i=1,6
          do j=1,6
             DO M=1,6
                DO N=1,6
                   X1=(xr(I)).par.ind_stoc(M)
                   X3=(xr(j)).par.ind_stoc(n)
                   X(I)%SIGMAf(J)=X(m)%E(n)*x1*x3+X(I)%SIGMAf(J)
                   X(I)%SIGMAf(J)=X(m)%SIGMA0(n)*x1*x3+X(I)%SIGMAf(J)
                enddo
             enddo
          enddo
       enddo

       CALL KILL(XT);CALL KILL(DISP);
       CALL KILL(X1,X3);CALL KILL(XR);

    endif

  END SUBROUTINE TRACK_FIBRE_S





  !   Misalignment routines
  SUBROUTINE MIS_FIBR(C,X,OU,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: OU,ENTERING

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             CALL ROT_YZ(C%CHART%ANG_IN(1),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)   ! ROTATIONS
             CALL ROT_XZ(C%CHART%ANG_IN(2),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)
             CALL ROT_XY(C%CHART%ANG_IN(3),X,OU)
             CALL TRANS(C%CHART%D_IN,X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)         ! TRANSLATION
          ELSE
             CALL ROT_YZ(C%CHART%ANG_OUT(1),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)  ! ROTATIONS
             CALL ROT_XZ(C%CHART%ANG_OUT(2),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)
             CALL ROT_XY(C%CHART%ANG_OUT(3),X,OU)
             CALL TRANS(C%CHART%D_OUT,X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)        ! TRANSLATION
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
             CALL TRANS(C%CHART%D_OUT,X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)        ! TRANSLATION
             CALL ROT_XY(C%CHART%ANG_OUT(3),X,OU)
             CALL ROT_XZ(C%CHART%ANG_OUT(2),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)
             CALL ROT_YZ(C%CHART%ANG_OUT(1),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)  ! ROTATIONS
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
             CALL TRANS(C%CHART%D_IN,X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)         ! TRANSLATION
             CALL ROT_XY(C%CHART%ANG_IN(3),X,OU)
             CALL ROT_XZ(C%CHART%ANG_IN(2),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)
             CALL ROT_YZ(C%CHART%ANG_IN(1),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)   ! ROTATIONS
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_FIBR

  SUBROUTINE MIS_FIBP(C,X,OU,ENTERING)  ! Misaligns polymorphic fibres in PTC order for forward and backward fibres
    implicit none
    TYPE(FIBRE),INTENT(INOUT):: C
    type(REAL_8), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: OU,ENTERING

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN
          IF(ENTERING) THEN
             CALL ROT_YZ(C%CHART%ang_in(1),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                ! rotations
             CALL ROT_XZ(C%CHART%ang_in(2),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)
             CALL ROT_XY(C%CHART%ang_in(3),X,OU)
             CALL TRANS(C%CHART%d_in,X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                       !translation
          ELSE
             CALL ROT_YZ(C%CHART%ang_out(1),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                ! rotations
             CALL ROT_XZ(C%CHART%ang_out(2),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)
             CALL ROT_XY(C%CHART%ang_out(3),X,OU)
             CALL TRANS(C%CHART%d_out,X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                       !translation
          ENDIF
       ELSE
          IF(ENTERING) THEN
             C%CHART%d_out(1)=-C%CHART%d_out(1)
             C%CHART%d_out(2)=-C%CHART%d_out(2)
             C%CHART%ang_out(3)=-C%CHART%ang_out(3)
             CALL TRANS(C%CHART%d_out,X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                       !translation
             CALL ROT_XY(C%CHART%ang_out(3),X,OU)
             CALL ROT_XZ(C%CHART%ang_out(2),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)
             CALL ROT_YZ(C%CHART%ang_out(1),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                ! rotations
             C%CHART%d_out(1)=-C%CHART%d_out(1)
             C%CHART%d_out(2)=-C%CHART%d_out(2)
             C%CHART%ang_out(3)=-C%CHART%ang_out(3)
          ELSE
             C%CHART%d_in(1)=-C%CHART%d_in(1)
             C%CHART%d_in(2)=-C%CHART%d_in(2)
             C%CHART%ang_in(3)=-C%CHART%ang_in(3)
             CALL TRANS(C%CHART%d_in,X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                       !translation
             CALL ROT_XY(C%CHART%ang_in(3),X,OU)
             CALL ROT_XZ(C%CHART%ang_in(2),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)
             CALL ROT_YZ(C%CHART%ang_in(1),X,C%MAGP%P%BETA0,OU,C%MAGP%P%TIME)                ! rotations
             C%CHART%d_in(1)=-C%CHART%d_in(1)
             C%CHART%d_in(2)=-C%CHART%d_in(2)
             C%CHART%ang_in(3)=-C%CHART%ang_in(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_FIBP


  SUBROUTINE MIS_FIBS(C,Y,OU,ENTERING) ! Misaligns envelope fibres in PTC order for forward and backward fibres
    implicit none
    TYPE(FIBRE),INTENT(INOUT):: C
    type(ENV_8), INTENT(INOUT):: Y(6)
    type(REAL_8) X(6)
    logical(lp),INTENT(IN):: OU,ENTERING
    CALL ALLOC(X,6)
    X=Y
    CALL MIS_FIB(C,X,OU,ENTERING)
    Y=X
    CALL KILL(X,6)
  END SUBROUTINE MIS_FIBS

END MODULE S_TRACKING
