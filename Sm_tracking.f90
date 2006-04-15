!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90

MODULE S_TRACKING
  USE S_FAMILY

  IMPLICIT NONE
  public
  logical(lp),TARGET :: ALWAYS_EXACT_PATCHING=.TRUE.
  type(fibre), pointer :: lost_fibre

  ! linked
  PRIVATE TRACK_LAYOUT_FLAG_R,TRACK_LAYOUT_FLAG_P,TRACK_LAYOUT_FLAG_S
  !  PRIVATE FIND_ORBIT_LAYOUT,FIND_ORBIT_M_LAYOUT,FIND_ENV_LAYOUT, FIND_ORBIT_LAYOUT_noda
  PRIVATE TRACK_LAYOUT_FLAG_R1,TRACK_LAYOUT_FLAG_P1,TRACK_LAYOUT_FLAG_S1
  PRIVATE MIS_FIBR,MIS_FIBP,MIS_FIBS,PATCH_FIBR,PATCH_FIBP,PATCH_FIBS
  PRIVATE TRACK_FIBRE_R,TRACK_FIBRE_P,TRACK_FIBRE_S
  PRIVATE TRACK_LAYOUT_FLAG_R1f,TRACK_LAYOUT_FLAG_P1f,TRACK_LAYOUT_FLAG_S1f
  PRIVATE TRACK_LAYOUT_FLAG_Rf,TRACK_LAYOUT_FLAG_Pf,TRACK_LAYOUT_FLAG_Sf
  ! old Sj_elements
  PRIVATE TRACKR,TRACKP,TRACKS
  logical(lp),TARGET :: other_program=.false.
  logical(lp),TARGET :: x_prime=.false.
  integer j_global
  ! END old Sj_elements

  ! TYPE UPDATING
  !    logical(lp) UPDATE
  ! END TYPE UPDATING



  !  TYPE (UPDATING), PARAMETER ::  COMPUTE= UPDATING(.TRUE.)
  LOGICAL(LP) :: COMPUTE = .FALSE.

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
     ! old Sj_elements
     !  INTERFACE TRACK
     MODULE PROCEDURE TRACKR
     MODULE PROCEDURE TRACKP
     MODULE PROCEDURE TRACKS
     !  END INTERFACE
     ! END old Sj_elements
  END INTERFACE


  INTERFACE TRACK_FLAG
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_R1f
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_P1f
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_S1f
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_Rf
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_Pf
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_Sf
  END INTERFACE


  INTERFACE PATCH_FIB
     MODULE PROCEDURE PATCH_FIBR
     MODULE PROCEDURE PATCH_FIBP
     MODULE PROCEDURE PATCH_FIBS
  END INTERFACE

  INTERFACE MIS_FIB
     MODULE PROCEDURE MIS_FIBR
     MODULE PROCEDURE MIS_FIBP
     MODULE PROCEDURE MIS_FIBS
  END INTERFACE


contains
  ! old Sj_elements


  SUBROUTINE TRACKR(EL,X,MID,K)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(ELEMENT),INTENT(INOUT):: EL
    TYPE(WORM),OPTIONAL, INTENT(INOUT):: MID
    TYPE(INTERNAL_STATE) K

    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
    if(other_program) then
       call track_R(x)
       return
    endif
    SELECT CASE(EL%KIND)
    CASE(KIND0)
       IF(PRESENT(MID)) CALL XMID(MID,X,0)
    case(KIND1)
       CALL TRACK(EL%D0,X,MID)
    case(KIND2)
       CALL TRACK(EL%K2,X,MID)
    case(KIND3)
       CALL TRACK(EL%K3,X,MID)
    case(KIND4)
       CALL TRACK(EL%C4,X,MID)
    case(KIND5)
       CALL TRACK(EL%S5,X,MID)
    case(KIND6)
       CALL TRACK(EL%T6,X,MID)
    case(KIND7)
       CALL TRACK(EL%T7,X,MID)
    case(KIND8)
       CALL TRACK(EL%S8,X,MID)
    case(KIND9)
       CALL TRACK(EL%S9,X,MID)
    case(KIND10)
       CALL TRACK(EL%TP10,X,MID)
    CASE(KIND11:KIND14)
       call TRACK(EL%MON14,X,MID)
    CASE(KIND15)
       call TRACK(EL%SEP15,X,MID)
    CASE(KIND16,KIND20)
       call TRACK(EL%K16,X,MID)
    CASE(KIND17)
       call TRACK(EL%S17,X,MID)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X,MID)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X,MID)
    CASE(KIND21)
       call TRACK(EL%CAV21,X,MID)
    CASE(KIND22)
       call TRACK(EL%M22,X,MID)
    CASE(KIND23)
       IF(PRESENT(MID)) CALL XMID(MID,X,0)
       if(el%p%dir==1) then
          call TRACK(EL%G23,X,1,K)
       else
          call TRACK(EL%G23,X,EL%G23%N,0,K)
       endif
       IF(PRESENT(MID)) CALL XMID(MID,X,1)
    case(KINDUSER1)
       call TRACK(EL%U1,X,MID)
    case(KINDUSER2)
       call TRACK(EL%U2,X,MID)
    case(KINDWIGGLER)
       call TRACK(EL%WI,X,MID)
    case(KINDPA)
       call TRACK(EL%PA,X,MID)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKR"
       CALL WRITE_E(0)
    END SELECT
  END SUBROUTINE TRACKR

  SUBROUTINE TRACKP(EL,X,MID,K)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(ELEMENTP),INTENT(INOUT):: EL
    TYPE(WORM_8),OPTIONAL, INTENT(INOUT):: MID
    TYPE(INTERNAL_STATE) K

    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
    if(other_program) then
       call track_p(x)
       return
    endif
    SELECT CASE(EL%KIND)
    CASE(KIND0)
       IF(PRESENT(MID)) CALL XMID(MID,X,0)
    case(KIND1)
       CALL TRACK(EL%D0,X,MID)
    case(KIND2)
       CALL TRACK(EL%K2,X)
    case(KIND3)
       CALL TRACK(EL%K3,X,MID)
    case(KIND4)
       CALL TRACK(EL%C4,X,MID)
    case(KIND5)
       CALL TRACK(EL%S5,X,MID)
    case(KIND6)
       CALL TRACK(EL%T6,X,MID)
    case(KIND7)
       CALL TRACK(EL%T7,X,MID)
    case(KIND8)
       CALL TRACK(EL%S8,X,MID)
    case(KIND9)
       CALL TRACK(EL%S9,X,MID)
    case(KIND10)
       CALL TRACK(EL%TP10,X,MID)
    CASE(KIND11:KIND14)
       call TRACK(EL%MON14,X,MID)
    CASE(KIND15)
       call TRACK(EL%SEP15,X,MID)
    CASE(KIND16,KIND20)
       call TRACK(EL%K16,X,MID)
    CASE(KIND17)
       call TRACK(EL%S17,X,MID)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X,MID)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X,MID)
    CASE(KIND21)
       call TRACK(EL%CAV21,X,MID)
    CASE(KIND22)
       call TRACK(EL%M22,X,MID)
    CASE(KIND23)
       IF(PRESENT(MID)) CALL XMID(MID,X,0)
       if(el%p%dir==1) then
          call TRACK(EL%G23,X,1,K)
       else
          call TRACK(EL%G23,X,EL%G23%N,0,K)
       endif
       IF(PRESENT(MID)) CALL XMID(MID,X,1)
    case(KINDUSER1)
       call TRACK(EL%U1,X,MID)
    case(KINDUSER2)
       call TRACK(EL%U2,X,MID)
    case(KINDWIGGLER)
       call TRACK(EL%WI,X,MID)
    case(KINDPA)
       call TRACK(EL%PA,X,MID)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKP"
       CALL WRITE_E(0)
    END SELECT
  END SUBROUTINE TRACKP

  SUBROUTINE TRACKS(EL,X,MID,K)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: X(6)
    TYPE(ELEMENTP),INTENT(INOUT):: EL
    TYPE(INNER_ENV_8_DATA),OPTIONAL, INTENT(INOUT):: MID
    TYPE(INTERNAL_STATE) K

    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
    if(other_program) then
       write(6,*) " OTHER PROGRAM FORBIDDEN "
       STOP 350
    endif
    SELECT CASE(EL%KIND)
    CASE(KIND0)
       !       IF(PRESENT(MID)) CALL XMID(MID,X,0)
    case(KIND1)
       CALL TRACK(EL%D0,X)
    case(KIND2)
       CALL TRACK(EL%K2,X)
    case(KIND3)
       CALL TRACK(EL%K3,X)
    case(KIND4)
       CALL TRACK(EL%C4,X)
    case(KIND5)
       CALL TRACK(EL%S5,X)
    case(KIND6)
       CALL TRACK(EL%T6,X)
    case(KIND7)
       CALL TRACK(EL%T7,X)
    case(KIND8)
       CALL TRACK(EL%S8,X)
    case(KIND9)
       CALL TRACK(EL%S9,X)
    case(KIND10)
       CALL TRACK(EL%TP10,X)
    CASE(KIND11:KIND14)
       call TRACK(EL%MON14,X)
    CASE(KIND15)
       call TRACK(EL%SEP15,X)
    CASE(KIND16,KIND20)
       call TRACK(EL%K16,X)
    CASE(KIND17)
       call TRACK(EL%S17,X)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X)
    CASE(KIND21)
       call TRACK(EL%CAV21,X)
    CASE(KIND22)
       call TRACK(EL%M22,X)
    CASE(KIND23)
       if(el%p%dir==1) then
          call TRACK(EL%G23,X,1,K)
       else
          call TRACK(EL%G23,X,EL%G23%N,0,K)
       endif
    case(KINDUSER1)
       call TRACK(EL%U1,X)
    case(KINDUSER2)
       call TRACK(EL%U2,X)
    case(KINDWIGGLER)
       call TRACK(EL%WI,X)
    case(KINDFITTED)
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1)= "KINDFITTED not supported "
       CALL WRITE_E(0)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKS"
       CALL WRITE_E(0)
    END SELECT



  END SUBROUTINE TRACKS

  ! END old Sj_elements

  recursive integer function TRACK_LAYOUT_FLAG_R1f(R,X,II1,k,X_IN)
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1

    call track(R,X,II1,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_R1f)

  end  function TRACK_LAYOUT_FLAG_R1f

  recursive   integer function TRACK_LAYOUT_FLAG_P1f(R,X,II1,k,X_IN)
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1

    call track(R,X,II1,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_P1f)

  end  function TRACK_LAYOUT_FLAG_P1f

  recursive   integer function TRACK_LAYOUT_FLAG_S1f(R,X,II1,k,X_IN)
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INNER_ENV_8_DATA), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1


    call track(R,X,II1,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_S1f)

  end  function TRACK_LAYOUT_FLAG_S1f

  recursive   SUBROUTINE TRACK_LAYOUT_FLAG_R1(R,X,II1,k,X_IN) ! Tracks real(dp) from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2

    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k,X_IN)
    if(c_%watch_user) ALLOW_TRACKING=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_FLAG_R1

  recursive   SUBROUTINE TRACK_LAYOUT_FLAG_P1(R,X,II1,k,X_IN) ! Tracks polymorphs from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2

    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k,X_IN)
    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_LAYOUT_FLAG_P1

  recursive   SUBROUTINE TRACK_LAYOUT_FLAG_S1(R,X,II1,k,X_IN) ! Tracks envelope from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),INTENT(INOUT):: R
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INNER_ENV_8_DATA), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2


    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k,X_IN)
    if(c_%watch_user) ALLOW_TRACKING=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_FLAG_S1

  recursive   integer function TRACK_LAYOUT_FLAG_Rf(R,X,I1,I2,k,X_IN) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER, INTENT(IN):: I1,I2

    call track(R,X,I1,I2,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_Rf)

  end  function TRACK_LAYOUT_FLAG_Rf

  recursive   integer function TRACK_LAYOUT_FLAG_Pf(R,X,I1,I2,k,X_IN) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(LAYOUT),INTENT(INOUT):: R ;TYPE(REAL_8), INTENT(INOUT):: X(6);
    INTEGER, INTENT(IN):: I1,I2; TYPE(INTERNAL_STATE) K;
    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN

    call track(R,X,I1,I2,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_Pf)

  end  function TRACK_LAYOUT_FLAG_Pf

  recursive   integer function TRACK_LAYOUT_FLAG_Sf(R,X,I1,I2,k,X_IN) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INNER_ENV_8_DATA), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2

    call track(R,X,I1,I2,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_Sf)

  end  function TRACK_LAYOUT_FLAG_Sf


  recursive   SUBROUTINE TRACK_LAYOUT_FLAG_R(R,X,I1,I2,k,X_IN) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J
    TYPE (fibre), POINTER :: C


    CALL RESET_APERTURE_FLAG



    call move_to(r,c,MOD_N(I1,R%N))

    if(i2>i1) then
       J=I1

       DO  WHILE(J<I2.AND.ASSOCIATED(C))
          j_global=j
          CALL TRACK(C,X,K,R%CHARGE,X_IN)

          C=>C%NEXT
          J=J+1
       ENDDO
    else
       J=I1

       DO  WHILE(J>I2.AND.ASSOCIATED(C))
          j_global=j

          c%dir=-c%dir
          CALL TRACK(C,X,K,R%CHARGE,X_IN)
          c%dir=-c%dir

          C=>C%previous
          J=J-1
       ENDDO

    endif


    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_LAYOUT_FLAG_R



  recursive   SUBROUTINE TRACK_LAYOUT_FLAG_P(R,X,I1,I2,K,X_IN) ! TRACKS POLYMORPHS FROM I1 TO I2 IN STATE K
    IMPLICIT NONE
    TYPE(LAYOUT),INTENT(INOUT):: R ;TYPE(REAL_8), INTENT(INOUT):: X(6);
    INTEGER, INTENT(IN):: I1,I2; TYPE(INTERNAL_STATE) K;
    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER J;

    TYPE (FIBRE), POINTER :: C


    CALL RESET_APERTURE_FLAG

    call move_to(r,c,MOD_N(I1,R%N))

    if(i2>i1) then
       J=I1

       DO  WHILE(J<I2.AND.ASSOCIATED(C))
          j_global=j
          CALL TRACK(C,X,K,R%CHARGE,X_IN)

          C=>C%NEXT
          J=J+1
       ENDDO
    else
       J=I1

       DO  WHILE(J>I2.AND.ASSOCIATED(C))
          j_global=j

          c%dir=-c%dir
          CALL TRACK(C,X,K,R%CHARGE,X_IN)
          c%dir=-c%dir

          C=>C%previous
          J=J-1
       ENDDO

    endif

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

    ! PATCHES
  END SUBROUTINE TRACK_LAYOUT_FLAG_P

  recursive   SUBROUTINE TRACK_LAYOUT_FLAG_S(R,X,I1,I2,k,X_IN) ! Tracks envelopes from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INNER_ENV_8_DATA), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER I,J,M,N
    ! YS SPECIFIC STUFF
    TYPE(DAMAP) ID,XT,DISP
    TYPE(REAL_8) XR(6),X1,X3
    real(dp) V(6)
    TYPE (fibre), POINTER :: C



    CALL RESET_APERTURE_FLAG

    call move_to(r,c,MOD_N(I1,R%N))

    !    ! new stuff with kind=3
    !    IF(k%para_in ) knob=.true.
    !    ! end new stuff with kind=3


    if(i2>i1) then
       J=I1

       DO  WHILE(J<I2.AND.ASSOCIATED(C))
          j_global=j

          CALL TRACK(C,X,K,R%CHARGE,X_IN)

          C=>C%NEXT
          J=J+1
       ENDDO
    else
       J=I1

       DO  WHILE(J>I2.AND.ASSOCIATED(C))
          j_global=j

          c%dir=-c%dir
          CALL TRACK(C,X,K,R%CHARGE,X_IN)
          c%dir=-c%dir

          C=>C%previous
          J=J-1
       ENDDO

    endif
    !    ! new stuff with kind=3
    !    knob=.FALSE.
    !    ! end new stuff with kind=3

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

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_FLAG_S


  recursive   SUBROUTINE TRACK_FIBRE_R(C,X,K,CHARGE,X_IN)
    implicit none
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER,optional, target, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    REAL(DP) ENT(3,3), A(3)
    integer,target :: charge1
    real(dp) xp



    IF(.NOT.CHECK_STABLE) return

    if(c_%x_prime) then
       P0=>C%MAG%P%P0C
       B0=>C%MAG%P%BETA0
       IF(C%MAG%P%exact)THEN
          IF(C%MAG%P%TIME)THEN
             xp=x(2)/root(one+two*X(5)/B0+X(5)**2-x(2)**2-x(4)**2)
             x(4)=x(4)/root(one+two*X(5)/B0+X(5)**2-x(2)**2-x(4)**2)
             x(2)=xp
          else
             xp=x(2)/root((one+x(5))**2-x(2)**2-x(4)**2)
             x(4)=x(4)/root((one+x(5))**2-x(2)**2-x(4)**2)
             x(2)=xp
          endif
       else
          IF(C%MAG%P%TIME)THEN
             x(2)=x(2)/root(one+two*X(5)/B0+X(5)**2)
             x(4)=x(4)/root(one+two*X(5)/B0+X(5)**2)
          else
             x(2)=x(2)/(one+x(5))
             x(4)=x(4)/(one+x(5))
          endif
       endif
    endif


    IF(PRESENT(X_IN)) then
       X_IN%F=>c ; X_IN%E%F=>C; X_IN%NST=>X_IN%E%NST;
    endif

    ! DIRECTIONAL VARIABLE
    C%MAG%P%DIR=>C%DIR
    if(present(charge)) then
       C%MAG%P%CHARGE=>CHARGE
    else
       charge1=1
       C%MAG%P%CHARGE=>CHARGE1
    endif
    !
    !    IF(.NOT.CHECK_STABLE) CHECK_STABLE=.TRUE.
    C%MAG=K
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-6)
       X_IN%POS(1)=X_IN%nst
    endif

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
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-5)

    ! The chart frame of reference is located here implicitely
    IF(PATCHG/=0.AND.PATCHG/=2) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_TRUE)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-4)
    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%a_T
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-3)

    CALL DTILTD(C%DIR,C%MAG%P%TILTD,1,X)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-2)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       ou = K%EXACTMIS.or.C%MAG%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-1)
       X_IN%POS(2)=X_IN%nst
    endif

    CALL TRACK(C%MAG,X,X_IN,K)

    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(3)=X_IN%nst
    endif

    IF(C%MAG%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%DIR,C%MAG%P%TILTD,2,X)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%b_T
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHG/=0.AND.PATCHG/=1) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_FALSE)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

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

    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(4)=X_IN%nst
    endif

    IF(PRESENT(X_IN))  THEN
       IF(X_IN%E%DO_SURVEY) THEN
          CALL G_FRAME(X_IN%E,ENT,A,-7)
          CALL  SURVEY(C,ENT,A,E_IN=X_IN%E)
       ELSE
          CALL SURVEY_INNER_MAG(X_IN%E)
       ENDIF
    ENDIF

    if(c_%x_prime) then
       P0=>C%MAG%P%P0C
       B0=>C%MAG%P%BETA0
       IF(C%MAG%P%exact)THEN
          IF(C%MAG%P%TIME)THEN
             xp=root(one+two*X(5)/B0+X(5)**2)*x(2)/root(one+x(2)**2+x(4)**2)
             x(4)=root(one+two*X(5)/B0+X(5)**2)*x(4)/root(one+x(2)**2+x(4)**2)
             x(2)=xp
          else
             xp=(one+x(5))*x(2)/root(one+x(2)**2+x(4)**2)
             x(4)=(one+x(5))*x(4)/root(one+x(2)**2+x(4)**2)
             x(2)=xp
          endif
       else
          IF(C%MAG%P%TIME)THEN
             x(2)=root(one+two*X(5)/B0+X(5)**2)*x(2)
             x(4)=root(one+two*X(5)/B0+X(5)**2)*x(4)
          else
             x(2)=(one+x(5))*x(2)
             x(4)=(one+x(5))*x(4)
          endif
       endif
    endif


    C%MAG=DEFAULT
    nullify(C%MAG%P%DIR)
    nullify(C%MAG%P%CHARGE)
    if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
       if(CHECK_MADX_APERTURE) c_%message="exceed absolute_aperture in TRACK_FIBRE_R"
       CHECK_STABLE=.false.
    endif
    lost_fibre=>c
  END SUBROUTINE TRACK_FIBRE_R

  recursive   SUBROUTINE TRACK_FIBRE_P(C,X,K,CHARGE,X_IN)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER, optional,TARGET, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) OU,PATCH
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (FIBRE), POINTER :: CN
    REAL(DP), POINTER :: P0,B0
    REAL(DP) ENT(3,3), A(3)
    integer,target :: charge1
    TYPE(REAL_8) xp

    IF(.NOT.CHECK_STABLE) return

    if(c_%x_prime) then
       call alloc(xp)  ! deallocated below
       P0=>C%MAG%P%P0C
       B0=>C%MAG%P%BETA0
       IF(C%MAG%P%exact)THEN
          IF(C%MAG%P%TIME)THEN
             xp=x(2)/sqrt(one+two*X(5)/B0+X(5)**2-x(2)**2-x(4)**2)
             x(4)=x(4)/sqrt(one+two*X(5)/B0+X(5)**2-x(2)**2-x(4)**2)
             x(2)=xp
          else
             xp=x(2)/sqrt((one+x(5))**2-x(2)**2-x(4)**2)
             x(4)=x(4)/sqrt((one+x(5))**2-x(2)**2-x(4)**2)
             x(2)=xp
          endif
       else
          IF(C%MAG%P%TIME)THEN
             x(2)=x(2)/sqrt(one+two*X(5)/B0+X(5)**2)
             x(4)=x(4)/sqrt(one+two*X(5)/B0+X(5)**2)
          else
             x(2)=x(2)/(one+x(5))
             x(4)=x(4)/(one+x(5))
          endif
       endif
    endif


    IF(PRESENT(X_IN)) then
       X_IN%F=>c ; X_IN%E%F=>C; X_IN%NST=>X_IN%E%NST;
    endif

    ! NEW STUFF WITH KIND=3: KNOB OF FPP IS SET TO TRUE IF NECESSARY
    IF(K%PARA_IN ) KNOB=.TRUE.
    ! END NEW STUFF WITH KIND=3

    ! DIRECTIONAL VARIABLE AND CHARGE IS PASSED TO THE ELEMENT

    C%MAGP%P%DIR=>C%DIR
    if(present(charge)) then
       C%MAGP%P%CHARGE=>CHARGE
    else
       charge1=1
       C%MAGP%P%CHARGE=>CHARGE1
    endif
    !

    ! PASSING THE STATE K TO THE ELEMENT
    C%MAGP=K
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF
    ! ENERGY PATCH
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-6)
       X_IN%POS(1)=X_IN%nst
    endif
    IF(PATCHE/=0.AND.PATCHE/=2) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! NO NEED TO PATCH IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0) THEN     ! NO NEED TO PATCH IF PATCHED BEFORE
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
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-5)


    ! POSITION PATCH
    IF(PATCHG/=0.AND.PATCHG/=2) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_TRUE)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-4)
    ! TIME PATCH
    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%A_T
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-3)

    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,1,X)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-2)
    ! MISALIGNMENTS AT THE ENTRANCE
    IF(C%MAGP%MIS) THEN
       OU = K%EXACTMIS.OR.C%MAGP%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-1)
       X_IN%POS(2)=X_IN%nst
    endif
    ! ************************************************************************
    !  THE ACTUAL MAGNET PROPAGATOR AS IT WOULD APPEAR IN A STANDARD CODE

    CALL TRACK(C%MAGP,X,X_IN,K)
    !
    ! ************************************************************************
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(3)=X_IN%nst
    endif


    ! MISALIGNMENTS AT THE EXIT
    IF(C%MAGP%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,2,X)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    !EXIT PATCH
    ! TIME PATCH
    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%b_T
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    ! POSITION PATCH
    IF(PATCHG/=0.AND.PATCHG/=1) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_FALSE)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    ! ENERGY PATCH
    IF(PATCHE/=0.AND.PATCHE/=1) THEN
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

    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(4)=X_IN%nst
    endif

    IF(PRESENT(X_IN))  THEN
       IF(X_IN%E%DO_SURVEY) THEN
          CALL G_FRAME(X_IN%E,ENT,A,-7)
          CALL  SURVEY(C,ENT,A,E_IN=X_IN%E)
       ELSE
          CALL SURVEY_INNER_MAG(X_IN%E)
       ENDIF
    ENDIF

    if(c_%x_prime) then
       P0=>C%MAG%P%P0C
       B0=>C%MAG%P%BETA0
       IF(C%MAG%P%exact)THEN
          IF(C%MAG%P%TIME)THEN
             xp=sqrt(one+two*X(5)/B0+X(5)**2)*x(2)/sqrt(one+x(2)**2+x(4)**2)
             x(4)=sqrt(one+two*X(5)/B0+X(5)**2)*x(4)/sqrt(one+x(2)**2+x(4)**2)
             x(2)=xp
          else
             xp=(one+x(5))*x(2)/sqrt(one+x(2)**2+x(4)**2)
             x(4)=(one+x(5))*x(4)/sqrt(one+x(2)**2+x(4)**2)
             x(2)=xp
          endif
       else
          IF(C%MAG%P%TIME)THEN
             x(2)=sqrt(one+two*X(5)/B0+X(5)**2)*x(2)
             x(4)=sqrt(one+two*X(5)/B0+X(5)**2)*x(4)
          else
             x(2)=(one+x(5))*x(2)
             x(4)=(one+x(5))*x(4)
          endif
       endif
       call kill(xp)
    endif

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
    endif
    lost_fibre=>c

  END SUBROUTINE TRACK_FIBRE_P


  recursive   SUBROUTINE TRACK_FIBRE_S(C,X,K,CHARGE,X_IN)   !,UPDATE
    implicit none
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(ENV_8), INTENT(INOUT):: X(6)
    TYPE(INNER_ENV_8_DATA), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(REAL_8)  Y(6)
    INTEGER, optional,target, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    !    TYPE(UPDATING), optional,intent(in):: UPDATE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    ! YS SPECIFIC STUFF
    TYPE(DAMAP) ID,XT,DISP
    TYPE(REAL_8) XR(6),X1,X3
    real(dp) V(6)
    INTEGER I,J,M,N
    integer,target :: charge1
    TYPE(REAL_8) xp

    IF(.NOT.CHECK_STABLE) return

    if(c_%x_prime) then
       call alloc(xp)  ! deallocated below
       call alloc(y)
       P0=>C%MAG%P%P0C
       B0=>C%MAG%P%BETA0
       IF(C%MAG%P%exact)THEN
          IF(C%MAG%P%TIME)THEN
             xp=y(2)/sqrt(one+two*y(5)/B0+y(5)**2-y(2)**2-y(4)**2)
             y(4)=y(4)/sqrt(one+two*y(5)/B0+y(5)**2-y(2)**2-y(4)**2)
             y(2)=xp
          else
             xp=y(2)/sqrt((one+y(5))**2-y(2)**2-y(4)**2)
             y(4)=y(4)/sqrt((one+y(5))**2-y(2)**2-y(4)**2)
             y(2)=xp
          endif
       else
          IF(C%MAG%P%TIME)THEN
             y(2)=y(2)/sqrt(one+two*y(5)/B0+y(5)**2)
             y(4)=y(4)/sqrt(one+two*y(5)/B0+y(5)**2)
          else
             y(2)=y(2)/(one+y(5))
             y(4)=y(4)/(one+y(5))
          endif
       endif
       call kill(y)
    endif



    ! new stuff with kind=3
    IF(k%para_in ) knob=.true.
    ! end new stuff with kind=3

    ! DIRECTIONAL VARIABLE
    C%MAGP%P%DIR=>C%DIR
    if(present(charge)) then
       C%MAGP%P%CHARGE=>CHARGE
    else
       charge1=1
       C%MAGP%P%CHARGE=>CHARGE1
    endif
    !

    C%MAGP=K
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF


    !   IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-6)

    IF(PATCHE/=0.AND.PATCHE/=2) THEN
       CALL ALLOC(Y)
       Y=X
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN   ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0) THEN   ! No need to patch IF PATCHED BEFORE
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
    !       IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-5)






    IF(PATCHG/=0.AND.PATCHG/=2) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_TRUE)
    ENDIF
    !       IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-4)


    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       CALL ALLOC(Y)
       Y=X
       Y(6)=Y(6)+C%PATCH%A_T
       X=Y
       CALL KILL(Y)
    ENDIF

    !       IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-3)


    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,1,X)
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-2)

    IF(C%MAGP%MIS) THEN
       ou = K%EXACTMIS.or.C%MAGP%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-1)

    CALL TRACK(C%MAGP,X,X_IN,K)

    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(C%MAGP%MIS) THEN
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,2,X)
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(.NOT.K%TOTALPATH)) THEN
       CALL ALLOC(Y)
       Y=X
       Y(6)=Y(6)+C%PATCH%b_T
       X=Y
       CALL KILL(Y)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHG/=0.AND.PATCHG/=1) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_FALSE)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHE/=0.AND.PATCHE/=1) THEN
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


    if(c_%x_prime) then
       call alloc(y)
       P0=>C%MAG%P%P0C
       B0=>C%MAG%P%BETA0
       IF(C%MAG%P%exact)THEN
          IF(C%MAG%P%TIME)THEN
             xp=sqrt(one+two*y(5)/B0+y(5)**2)*y(2)/sqrt(one+y(2)**2+y(4)**2)
             y(4)=sqrt(one+two*y(5)/B0+y(5)**2)*y(4)/sqrt(one+y(2)**2+y(4)**2)
             y(2)=xp
          else
             xp=(one+y(5))*y(2)/sqrt(one+y(2)**2+y(4)**2)
             y(4)=(one+y(5))*y(4)/sqrt(one+y(2)**2+y(4)**2)
             y(2)=xp
          endif
       else
          IF(C%MAG%P%TIME)THEN
             y(2)=sqrt(one+two*y(5)/B0+y(5)**2)*y(2)
             y(4)=sqrt(one+two*y(5)/B0+y(5)**2)*y(4)
          else
             y(2)=(one+y(5))*y(2)
             y(4)=(one+y(5))*y(4)
          endif
       endif
       call kill(xp)
       call kill(y)
    endif



    C%MAGP=DEFAULT


    if(abs(x(1)%v)+abs(x(3)%v)>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
       if(CHECK_MADX_APERTURE) c_%message="exceed absolute_aperture in TRACK_FIBRE_S"
       CHECK_STABLE=.false.
    endif
    lost_fibre=>c

    ! new stuff with kind=3
    knob=.FALSE.
    ! end new stuff with kind=3
    if(COMPUTE) then
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
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

  END SUBROUTINE TRACK_FIBRE_S


  SUBROUTINE PATCH_FIBR(C,X,PATCH,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: PATCH,ENTERING

    IF(ENTERING) THEN
       X(3)=C%PATCH%A_X1*X(3);X(4)=C%PATCH%A_X1*X(4);
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       X(3)=C%PATCH%A_X2*X(3);X(4)=C%PATCH%A_X2*X(4);
    ELSE
       X(3)=C%PATCH%B_X1*X(3);X(4)=C%PATCH%B_X1*X(4);
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%B_D,X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       X(3)=C%PATCH%B_X2*X(3);X(4)=C%PATCH%B_X2*X(4);
    ENDIF


  END SUBROUTINE PATCH_FIBR


  SUBROUTINE PATCH_FIBP(C,X,PATCH,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: PATCH,ENTERING

    IF(ENTERING) THEN
       X(3)=C%PATCH%A_X1*X(3);X(4)=C%PATCH%A_X1*X(4);
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       X(3)=C%PATCH%A_X2*X(3);X(4)=C%PATCH%A_X2*X(4);
    ELSE
       X(3)=C%PATCH%B_X1*X(3);X(4)=C%PATCH%B_X1*X(4);
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X,PATCH)
       CALL TRANS(C%PATCH%B_D,X,C%MAG%P%BETA0,PATCH,C%MAG%P%TIME)
       X(3)=C%PATCH%B_X2*X(3);X(4)=C%PATCH%B_X2*X(4);
    ENDIF


  END SUBROUTINE PATCH_FIBP

  SUBROUTINE PATCH_FIBS(C,Y,PATCH,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(ENV_8), INTENT(INOUT):: Y(6)
    TYPE(REAL_8) X(6)
    logical(lp),INTENT(IN):: PATCH,ENTERING

    CALL ALLOC(X)
    X=Y
    CALL PATCH_FIB(C,X,PATCH,ENTERING)

    Y=X

    CALL KILL(X)

  END SUBROUTINE PATCH_FIBS





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

  SUBROUTINE TRACK_R(X)
    IMPLICIT NONE
    REAL(DP) X(6),x6,xp,yp,x5
    INTEGER icharef
    COMMON/ptc/ icharef


    if(j_global==1) return  ! skipping OBJECT OF ZGOUBI = TRACKING COMMAND INTERNAL TO ZGOUBI
    icharef=0

    x(1)=x(1)*c_100
    x(3)=x(3)*c_100
    x6=x(6)*c_100

    xp=x(2)/root((one+x(5))**2-x(2)**2-x(4)**2)
    yp=x(4)/root((one+x(5))**2-x(2)**2-x(4)**2)
    x(2)=atan(xp)*c_1d3
    x(4)=atan(yp/root(one+xp**2))*c_1d3

    x(6)=x(5)
    x(5)=x6

    !call track_z(x,j_global,j_global)

    x6=x(5)/c_100
    x(5)=x(6)
    x(6)=x6

    x(1)=x(1)/c_100
    x(3)=x(3)/c_100
    xp=tan(x(2)/c_1d3)
    yp=tan(x(4)/c_1d3)*root(one+xp**2)

    x(2)=(one+x(5))*xp/root(one+xp**2+yp**2)
    x(4)=(one+x(5))*yp/root(one+xp**2+yp**2)

    icharef=1

  END SUBROUTINE TRACK_R

  SUBROUTINE TRACK_P(X)
    IMPLICIT NONE
    TYPE(REAL_8) X(6)

    ! track_zp is a fortran external routine using numerical differentiation
    !call track_zp(x,j_global,j_global)
    WRITE(6,*) " NOT SUPPORTED "
    STOP 111
  END SUBROUTINE TRACK_P


END MODULE S_TRACKING
