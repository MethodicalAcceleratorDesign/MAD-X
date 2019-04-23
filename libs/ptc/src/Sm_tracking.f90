!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN


MODULE S_TRACKING
  USE S_FAMILY

  IMPLICIT NONE
  public
  logical(lp),TARGET :: ALWAYS_EXACT_PATCHING=.TRUE.
  !  type(fibre), pointer :: lost_fibre
  !  type(integration_node), pointer :: lost_node

  ! linked
  PRIVATE TRACK_LAYOUT_FLAG_R,TRACK_LAYOUT_FLAG_P
  !  PRIVATE FIND_ORBIT_LAYOUT,FIND_ORBIT_M_LAYOUT,FIND_ENV_LAYOUT, FIND_ORBIT_LAYOUT_noda
  PRIVATE TRACK_LAYOUT_FLAG_R1,TRACK_LAYOUT_FLAG_P1
  PRIVATE MIS_FIBR,MIS_FIBP,PATCH_FIBR,PATCH_FIBP
  PRIVATE TRACK_FIBRE_R,TRACK_FIBRE_P
  PRIVATE TRACK_LAYOUT_FLAG_R1f,TRACK_LAYOUT_FLAG_P1f
  PRIVATE TRACK_LAYOUT_FLAG_Rf,TRACK_LAYOUT_FLAG_Pf
  private TRACK_fibre_based_R,TRACK_fibre_based_P
  ! old Sj_elements
  ! END old Sj_elements

  ! TYPE UPDATING
  !    logical(lp) UPDATE
  ! END TYPE UPDATING



  !  TYPE (UPDATING), PARAMETER ::  COMPUTE= UPDATING(.TRUE.)
  LOGICAL :: COMPUTE = .FALSE.

  INTERFACE TRACK
     ! linked
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_R
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_P
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_R1
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_P1
     MODULE PROCEDURE TRACK_FIBRE_R
     MODULE PROCEDURE TRACK_FIBRE_P
     MODULE PROCEDURE TRACK_fibre_based_R
     MODULE PROCEDURE TRACK_fibre_based_P
     ! old Sj_elements
     ! END old Sj_elements
  END INTERFACE


  INTERFACE TRACK_FIBRE_SINGLE
     MODULE PROCEDURE TRACK_FIBRE_R
     MODULE PROCEDURE TRACK_FIBRE_P
  END INTERFACE

  INTERFACE TRACK_FLAG
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_R1f
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_P1f
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_Rf
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_Pf
  END INTERFACE


  INTERFACE PATCH_FIB
     MODULE PROCEDURE PATCH_FIBR
     MODULE PROCEDURE PATCH_FIBP
  END INTERFACE

  INTERFACE MIS_FIB
     MODULE PROCEDURE MIS_FIBR
     MODULE PROCEDURE MIS_FIBP
  END INTERFACE


contains
  ! old Sj_elements


  ! END old Sj_elements

  !  recursive
  integer function TRACK_LAYOUT_FLAG_R1f(R,X,II1,k,X_IN)
    implicit none
    TYPE(layout),target,INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1

    call track(R,X,II1,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_R1f)
    !    call RESET_APERTURE_FLAG(my_false)
  end  function TRACK_LAYOUT_FLAG_R1f

  !  recursive
  integer function TRACK_LAYOUT_FLAG_P1f(R,X,II1,k)
    implicit none
    TYPE(layout),target,INTENT(INOUT):: R
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    !    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1

    call track(R,X,II1,k)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_P1f)
    !    call RESET_APERTURE_FLAG(my_false)

  end  function TRACK_LAYOUT_FLAG_P1f

  !  recursive
  SUBROUTINE TRACK_LAYOUT_FLAG_R1(R,X,II1,k,X_IN) ! Tracks real(dp) from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),target,INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2

    !    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k,X_IN)
    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_FLAG_R1

  !  recursive
  SUBROUTINE TRACK_LAYOUT_FLAG_P1(R,X,II1,k) ! Tracks polymorphs from II1 to the end or back to II1 if closed
    implicit none
    TYPE(layout),target,INTENT(INOUT):: R
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    !    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: II1
    INTEGER II2

    !    CALL RESET_APERTURE_FLAG

    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK(R,X,II1,II2,k)
    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_LAYOUT_FLAG_P1

  !  recursive
  integer function TRACK_LAYOUT_FLAG_Rf(R,X,I1,I2,k,X_IN) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER, INTENT(IN):: I1,I2

    call track(R,X,I1,I2,k,X_IN)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_Rf)

  end  function TRACK_LAYOUT_FLAG_Rf

  !  recursive
  integer function TRACK_LAYOUT_FLAG_Pf(R,X,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(LAYOUT),target,INTENT(INOUT):: R ;
    TYPE(REAL_8), INTENT(INOUT):: X(6);
    INTEGER, INTENT(IN):: I1,I2; TYPE(INTERNAL_STATE) K;
    !    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN

    call track(R,X,I1,I2,k)
    call PRODUCE_APERTURE_FLAG(TRACK_LAYOUT_FLAG_Pf)

  end  function TRACK_LAYOUT_FLAG_Pf




  SUBROUTINE TRACK_fibre_based_R(X,k,fibre1,fibre2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE (fibre), POINTER :: fibre1
    TYPE (fibre), optional, POINTER :: fibre2
    TYPE (fibre), POINTER :: C,c1,c2,last

    c1=>fibre1
    if(present(fibre2)) then
       c2=>fibre2
       nullify(last)
    else
       if(fibre1%parent_layout%closed) then
          last=>fibre1%previous
          c2=>last
       else
          last=>fibre1%parent_layout%end
          c2=>fibre1%parent_layout%end
       endif
    endif


    c=>c1



    DO  WHILE(.not.ASSOCIATED(C,c2))

       CALL TRACK(C,X,K)

       if(.not.check_stable) then
         write(messagelost,*) "Error in tracking ",c%mag%name," ", messagelost(:len_trim(messagelost))
         exit
       endif  

       C=>C%NEXT
    ENDDO

    if(associated(last).and.check_stable) then
       CALL TRACK(last,X,K)
    endif

    C_%STABLE_DA=.true.



  END SUBROUTINE TRACK_fibre_based_R


  SUBROUTINE TRACK_fibre_based_p(X,k,fibre1,fibre2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    type(real_8), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE (fibre), POINTER :: fibre1
    TYPE (fibre), optional, POINTER :: fibre2
    TYPE (fibre), POINTER :: C,c1,c2,last

    c1=>fibre1
    if(present(fibre2)) then
       c2=>fibre2
       nullify(last)
    else
       if(fibre1%parent_layout%closed) then
          last=>fibre1%previous
          c2=>last
       else
          last=>fibre1%parent_layout%end
          c2=>fibre1%parent_layout%end
       endif
    endif


    c=>c1



    DO  WHILE(.not.ASSOCIATED(C,c2))

       CALL TRACK(C,X,K)

       if(.not.check_stable) then
         write(messagelost,*) "Error in tracking ",c%mag%name," ", messagelost(:len_trim(messagelost))
         exit
       endif  

       C=>C%NEXT

    ENDDO

    if(associated(last).and.check_stable) then
       CALL TRACK(last,X,K)
    endif

    C_%STABLE_DA=.true.


  END SUBROUTINE TRACK_fibre_based_p





  SUBROUTINE TRACK_LAYOUT_FLAG_R(R,X,I1,I2,k,X_IN) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (fibre), POINTER :: C


    ! CALL RESET_APERTURE_FLAG



    call move_to(r,c,I1)

    if(i2>=i1) then
       i22=i2
    else
       i22=r%n+i2
    endif

    !    if(i2>i1) then
    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))
       CALL TRACK(C,X,K,X_IN=X_IN)  !,C%CHARGE
       !       CALL TRACK(C,X,K,R%CHARGE,X_IN)

       if(.not.check_stable) then
         write(messagelost,*) "Error in tracking ",c%mag%name," ", messagelost(:len_trim(messagelost))
         exit
       endif  

       C=>C%NEXT
       J=J+1
    ENDDO

    C_%STABLE_DA=.true.

    !    else
    !       J=I1
    !
    !       DO  WHILE(J>I2.AND.ASSOCIATED(C))
    !          j_global=j
    !
    !          c%dir=-c%dir
    !          CALL TRACK(C,X,K,R%CHARGE,X_IN)
    !          c%dir=-c%dir
    !
    !          C=>C%previous
    !          J=J-1
    !       ENDDO
    !
    !    endif


    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_LAYOUT_FLAG_R



  !  recursive
  SUBROUTINE TRACK_LAYOUT_FLAG_P(R,X,I1,I2,K) ! TRACKS POLYMORPHS FROM I1 TO I2 IN STATE K
    IMPLICIT NONE
    TYPE(LAYOUT),target,INTENT(INOUT):: R ;TYPE(REAL_8), INTENT(INOUT):: X(6);
    INTEGER, INTENT(IN):: I1,I2; TYPE(INTERNAL_STATE) K;
    !    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    INTEGER J,I22

    TYPE (FIBRE), POINTER :: C


    !  CALL RESET_APERTURE_FLAG

    call move_to(r,c,I1)

    if(i2>=i1) then
       i22=i2
    else
       i22=r%n+i2
    endif

    !    if(i2>i1) then
    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))
       CALL TRACK(C,X,K)  !,C%CHARGE
       !       CALL TRACK(C,X,K,R%CHARGE)
       if(.not.check_stable) then
         write(messagelost,*) "Error in tracking ",c%mag%name," ", messagelost(:len_trim(messagelost))
         exit
       endif  

       C=>C%NEXT
       J=J+1
    ENDDO

    C_%STABLE_DA=.true.

    !    else
    !       J=I1

    !       DO  WHILE(J>I2.AND.ASSOCIATED(C))
    !          j_global=j

    !          c%dir=-c%dir
    !          CALL TRACK(C,X,K,R%CHARGE,X_IN)
    !          c%dir=-c%dir

    !          C=>C%previous
    !          J=J-1
    !       ENDDO

    !    endif

    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

    ! PATCHES
  END SUBROUTINE TRACK_LAYOUT_FLAG_P

  !  recursive
  !  SUBROUTINE TRACK_FIBRE_R(C,X,K,CHARGE,X_IN)
  SUBROUTINE TRACK_FIBRE_R(C,X,K,X_IN)
    implicit none
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    !    INTEGER,optional, target, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    REAL(DP) ENT(3,3), A(3)

    ! real(dp), POINTER :: BETA0,GAMMA0I,GAMBET,P0C,MASS0
    !INTEGER, POINTER :: CHARGE


    IF(.NOT.CHECK_STABLE) then
       CALL RESET_APERTURE_FLAG
    endif
    !    C%MAG%P%p0c=>c%p0c
    C%MAG%P%beta0=>c%beta0
    C%MAG%P%GAMMA0I=>c%GAMMA0I
    C%MAG%P%GAMBET=>c%GAMBET
    C%MAG%P%CHARGE=>c%CHARGE
    ! DIRECTIONAL VARIABLE
    C%MAG%P%DIR=>C%DIR



    IF(PRESENT(X_IN)) then
       X_IN%F=>c ; X_IN%E%F=>C; X_IN%NST=>X_IN%E%NST;
    endif

    !
    !    IF(.NOT.CHECK_STABLE) CHECK_STABLE=.TRUE.
    !FRONTAL PATCH
    !    IF(ASSOCIATED(C%PATCH)) THEN
    PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    !    ELSE
    !       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    !    ENDIF
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-6)
       X_IN%POS(1)=X_IN%nst
    endif

    IF(PATCHE/=0.AND.PATCHE/=2.AND.PATCHE/=5) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN).and.PATCHE/=4) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0.or.CN%PATCH%ENERGY==1.or.CN%PATCH%ENERGY==4) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAG%P%P0C
             B0=>CN%BETA0

             X(2)=X(2)*P0/C%MAG%P%P0C
             X(4)=X(4)*P0/C%MAG%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=root(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/C%MAG%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAG%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAG%P%P0C-1.0_dp
             ENDIF
          ENDIF ! No need to patch
       else  ! ASSOCIATED
             P0=>C%PATCH%P0b
             B0=>C%PATCH%B0b

             X(2)=X(2)*P0/C%MAG%P%P0C
             X(4)=X(4)*P0/C%MAG%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=root(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/C%MAG%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAG%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAG%P%P0C-1.0_dp
             ENDIF           
       ENDIF ! ASSOCIATED

    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-5)

    ! The chart frame of reference is located here implicitely
    IF(PATCHG==1.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-4)
    IF(PATCHT/=0.AND.PATCHT/=2.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%a_T  !/c%beta0
      else
       X(6)=X(6)-C%PATCH%a_L
      endif
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-3)

    CALL DTILTD(C%MAG%P%TILTD,1,X)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-2)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       ou = ALWAYS_EXACTMIS    !K%EXACTMIS.or.
       CALL MIS_FIB(C,X,k,OU,DONEITT)
    ENDIF
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-1)
       X_IN%POS(2)=X_IN%nst
    endif

    CALL TRACK(C%MAG,X,K,X_IN)
 

    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(3)=X_IN%nst
    endif

    IF(C%MAG%MIS) THEN
       CALL MIS_FIB(C,X,k,OU,DONEITF)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%MAG%P%TILTD,2,X)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%b_T   !/c%beta0
      else
       X(6)=X(6)-C%PATCH%b_L
      endif
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PATCHG==2.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_FALSE)
    ENDIF
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    ! The CHART frame of reference is located here implicitely

    IF(PATCHE/=0.AND.PATCHE/=1.AND.PATCHE/=4) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(ASSOCIATED(CN).and.PATCHE/=5) THEN ! ASSOCIATED       
!       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAG%P%P0C
       B0=>CN%BETA0
       X(2)=X(2)*C%MAG%P%P0C/P0
       X(4)=X(4)*C%MAG%P%P0C/P0
       IF(k%TIME.or.recirculator_cheat)THEN
          X(5)=root(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW
          X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp
       ENDIF
    else
             P0=>C%PATCH%P0b ! 8/31/2016
             B0=>C%PATCH%B0b ! 8/31/2016

             X(2)=X(2)*C%MAG%P%P0C/P0 ! 8/31/2016
             X(4)=X(4)*C%MAG%P%P0C/P0 ! 8/31/2016
       IF(k%TIME.or.recirculator_cheat)THEN ! 8/31/2016
          X(5)=root(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD ! 8/31/2016
          X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW ! 8/31/2016
          X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0) ! 8/31/2016
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp ! 8/31/2016
       ENDIF        
    ENDIF
    endif ! associated
    
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

    !    endif ! new 2010

    if(abs(x(1))+abs(x(3))>absolute_aperture.or.abs(x(6))>t_aperture) then   !.or.(.not.CHECK_MADX_APERTURE)) then
       messageLOST="exceed absolute_aperture in TRACK_FIBRE_R"
       xlost=x
       CHECK_STABLE=.false.
    endif
    if(.not.check_stable ) lost_fibre=>c

  END SUBROUTINE TRACK_FIBRE_R

  !  recursive
  !  SUBROUTINE TRACK_FIBRE_P(C,X,K,CHARGE)
  SUBROUTINE TRACK_FIBRE_P(C,X,K)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    !    TYPE(WORM_8), OPTIONAL,INTENT(INOUT):: X_IN
    !   INTEGER, optional,TARGET, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) OU,PATCH
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (FIBRE), POINTER :: CN
    REAL(DP), POINTER :: P0,B0

    IF(.NOT.CHECK_STABLE) then
       CALL RESET_APERTURE_FLAG
    endif
    !    C%MAGp%P%p0c=>c%p0c
    C%MAGp%P%beta0=>c%beta0
    C%MAGp%P%GAMMA0I=>c%GAMMA0I
    C%MAGp%P%GAMBET=>c%GAMBET
    C%MAGp%P%CHARGE=>c%CHARGE
    C%MAGP%P%DIR=>C%DIR
    !    if(present(charge)) then
    !       C%MAGP%P%CHARGE=>CHARGE
    !    endif

    ! NEW STUFF WITH KIND=3: KNOB OF FPP IS SET TO TRUE IF NECESSARY
    IF(K%PARA_IN ) KNOB=.TRUE.
    PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    IF(PATCHE/=0.AND.PATCHE/=2.AND.PATCHE/=5) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN).and.PATCHE/=4) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! NO NEED TO PATCH IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0.or.CN%PATCH%ENERGY==1.or.CN%PATCH%ENERGY==4) THEN     ! NO NEED TO PATCH IF PATCHED BEFORE
             P0=>CN%MAGP%P%P0C
             B0=>CN%BETA0

             X(2)=X(2)*P0/C%MAGP%P%P0C
             X(4)=X(4)*P0/C%MAGP%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=SQRT(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAGP%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(SQRT(1.0_dp/C%MAGP%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAGP%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAGP%P%P0C-1.0_dp
             ENDIF
          ENDIF ! NO NEED TO PATCH
       else ! ASSOCIATED
             P0=>C%PATCH%P0b
             B0=>C%PATCH%B0b

             X(2)=X(2)*P0/C%MAGP%P%P0C
             X(4)=X(4)*P0/C%MAGP%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=SQRT(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAGP%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(SQRT(1.0_dp/C%MAGP%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAGP%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAGP%P%P0C-1.0_dp
             ENDIF           
      ENDIF ! ASSOCIATED

    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-5)


    ! POSITION PATCH
    IF(PATCHG==1.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-4)
    ! TIME PATCH
    IF(PATCHT/=0.AND.PATCHT/=2.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%a_T    !/c%beta0
      else
       X(6)=X(6)-C%PATCH%a_L
      endif
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-3)

    CALL DTILTD(C%MAGP%P%TILTD,1,X)
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-2)
    ! MISALIGNMENTS AT THE ENTRANCE
    IF(C%MAGP%MIS) THEN
       OU =ALWAYS_EXACTMIS   ! K%EXACTMIS.OR.
       CALL MIS_FIB(C,X,k,OU,DONEITT)
    ENDIF

    CALL TRACK(C%MAGP,X,K)
    !    if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then ! new 2010
    !       if(CHECK_MADX_APERTURE) c_%message="exceed absolute_aperture in TRACK_FIBRE_P"
    !       CHECK_STABLE=.false.
    !    else ! new 2010



    ! MISALIGNMENTS AT THE EXIT
    IF(C%MAGP%MIS) THEN
       CALL MIS_FIB(C,X,k,OU,DONEITF)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    CALL DTILTD(C%MAGP%P%TILTD,2,X)
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    !EXIT PATCH
    ! TIME PATCH
    IF(PATCHT/=0.AND.PATCHT/=1.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%b_T   !/c%beta0
      else
       X(6)=X(6)-C%PATCH%b_L
      endif
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    ! POSITION PATCH
    IF(PATCHG==2.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_FALSE)
    ENDIF
    !    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    ! ENERGY PATCH
    IF(PATCHE/=0.AND.PATCHE/=1.AND.PATCHE/=4) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
!       IF(.NOT.ASSOCIATED(CN)) CN=>C
       IF(ASSOCIATED(CN).and.PATCHE/=5) THEN ! ASSOCIATED        
       P0=>CN%MAGP%P%P0C
       B0=>CN%BETA0
       X(2)=X(2)*C%MAGP%P%P0C/P0
       X(4)=X(4)*C%MAGP%P%P0C/P0
       IF(k%TIME.or.recirculator_cheat)THEN
          X(5)=SQRT(1.0_dp+2.0_dp*X(5)/C%MAGP%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAGP%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW
          X(5)=(2.0_dp*X(5)+X(5)**2)/(SQRT(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAGP%P%P0C/P0-1.0_dp
       ENDIF
    
       else ! ASSOCIATED
             P0=>C%PATCH%P0b ! 8/31/2016
             B0=>C%PATCH%B0b ! 8/31/2016

             X(2)=X(2)*C%MAG%P%P0C/P0 ! 8/31/2016
             X(4)=X(4)*C%MAG%P%P0C/P0 ! 8/31/2016
       IF(k%TIME.or.recirculator_cheat)THEN ! 8/31/2016
          X(5)=sqrt(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD ! 8/31/2016
          X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW ! 8/31/2016
          X(5)=(2.0_dp*X(5)+X(5)**2)/(sqrt(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0) ! 8/31/2016
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp ! 8/31/2016
       ENDIF        
      
      ENDIF ! ASSOCIATED  
ENDIF  
    !   endif ! new 2010


    ! KNOB IS RETURNED TO THE PTC DEFAULT
    ! NEW STUFF WITH KIND=3
    KNOB=ALWAYS_knobs
    ! END NEW STUFF WITH KIND=3

    ! new 2010
    if(abs(x(1))+abs(x(3))>absolute_aperture.or.abs(x(6))>t_aperture) then   !.or.(.not.CHECK_MADX_APERTURE)) then
       messageLOST="exceed absolute_aperture in TRACK_FIBRE_P"
       xlost=x
       CHECK_STABLE=.false.
    endif
    if(.not.check_stable ) lost_fibre=>c

  END SUBROUTINE TRACK_FIBRE_P



  SUBROUTINE PATCH_FIBR(C,X,k,PATCH,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: PATCH,ENTERING
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ENTERING) THEN
       X(3)=C%PATCH%A_X1*X(3);X(4)=C%PATCH%A_X1*X(4);
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAG%P%BETA0,PATCH,k%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAG%P%BETA0,PATCH,k%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X)  !,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAG%P%BETA0,PATCH,k%TIME)
       X(3)=C%PATCH%A_X2*X(3);X(4)=C%PATCH%A_X2*X(4);
    ELSE
       X(3)=C%PATCH%B_X1*X(3);X(4)=C%PATCH%B_X1*X(4);
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAG%P%BETA0,PATCH,k%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAG%P%BETA0,PATCH,k%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X)  !,PATCH)
       CALL TRANS(C%PATCH%B_D,X,C%MAG%P%BETA0,PATCH,k%TIME)
       X(3)=C%PATCH%B_X2*X(3);X(4)=C%PATCH%B_X2*X(4);
    ENDIF


  END SUBROUTINE PATCH_FIBR


  SUBROUTINE PATCH_FIBP(C,X,k,PATCH,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: PATCH,ENTERING
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ENTERING) THEN
       X(3)=C%PATCH%A_X1*X(3);X(4)=C%PATCH%A_X1*X(4);
       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAGP%P%BETA0,PATCH,k%TIME)
       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAGP%P%BETA0,PATCH,k%TIME)
       CALL ROT_XY(C%PATCH%A_ANG(3),X)  !,PATCH)
       CALL TRANS(C%PATCH%A_D,X,C%MAGP%P%BETA0,PATCH,k%TIME)
       X(3)=C%PATCH%A_X2*X(3);X(4)=C%PATCH%A_X2*X(4);
    ELSE
       X(3)=C%PATCH%B_X1*X(3);X(4)=C%PATCH%B_X1*X(4);
       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAGP%P%BETA0,PATCH,k%TIME)
       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAGP%P%BETA0,PATCH,k%TIME)
       CALL ROT_XY(C%PATCH%B_ANG(3),X)  !,PATCH)
       CALL TRANS(C%PATCH%B_D,X,C%MAGP%P%BETA0,PATCH,k%TIME)
       X(3)=C%PATCH%B_X2*X(3);X(4)=C%PATCH%B_X2*X(4);
    ENDIF


  END SUBROUTINE PATCH_FIBP

  !   Misalignment routines
  SUBROUTINE MIS_FIBR(C,X,k,OU,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: OU,ENTERING
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             CALL ROT_YZ(C%CHART%ANG_IN(1),X,C%MAG%P%BETA0,OU,k%TIME)   ! ROTATIONS
             CALL ROT_XZ(C%CHART%ANG_IN(2),X,C%MAG%P%BETA0,OU,k%TIME)
             CALL ROT_XY(C%CHART%ANG_IN(3),X)  !,OU)
             CALL TRANS(C%CHART%D_IN,X,C%MAG%P%BETA0,OU,k%TIME)         ! TRANSLATION
          ELSE
             CALL ROT_YZ(C%CHART%ANG_OUT(1),X,C%MAG%P%BETA0,OU,k%TIME)  ! ROTATIONS
             CALL ROT_XZ(C%CHART%ANG_OUT(2),X,C%MAG%P%BETA0,OU,k%TIME)
             CALL ROT_XY(C%CHART%ANG_OUT(3),X)  !,OU)
             CALL TRANS(C%CHART%D_OUT,X,C%MAG%P%BETA0,OU,k%TIME)        ! TRANSLATION
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
             CALL TRANS(C%CHART%D_OUT,X,C%MAG%P%BETA0,OU,k%TIME)        ! TRANSLATION
             CALL ROT_XY(C%CHART%ANG_OUT(3),X)  !,OU)
             CALL ROT_XZ(C%CHART%ANG_OUT(2),X,C%MAG%P%BETA0,OU,k%TIME)
             CALL ROT_YZ(C%CHART%ANG_OUT(1),X,C%MAG%P%BETA0,OU,k%TIME)  ! ROTATIONS
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
             CALL TRANS(C%CHART%D_IN,X,C%MAG%P%BETA0,OU,k%TIME)         ! TRANSLATION
             CALL ROT_XY(C%CHART%ANG_IN(3),X)  !,OU)
             CALL ROT_XZ(C%CHART%ANG_IN(2),X,C%MAG%P%BETA0,OU,k%TIME)
             CALL ROT_YZ(C%CHART%ANG_IN(1),X,C%MAG%P%BETA0,OU,k%TIME)   ! ROTATIONS
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_FIBR

  SUBROUTINE MIS_FIBP(C,X,k,OU,ENTERING)  ! Misaligns polymorphic fibres in PTC order for forward and backward fibres
    implicit none
    TYPE(FIBRE),INTENT(INOUT):: C
    type(REAL_8), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: OU,ENTERING
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN
          IF(ENTERING) THEN
             CALL ROT_YZ(C%CHART%ang_in(1),X,C%MAGP%P%BETA0,OU,k%TIME)                ! rotations
             CALL ROT_XZ(C%CHART%ang_in(2),X,C%MAGP%P%BETA0,OU,k%TIME)
             CALL ROT_XY(C%CHART%ang_in(3),X)  !,OU)
             CALL TRANS(C%CHART%d_in,X,C%MAGP%P%BETA0,OU,k%TIME)                       !translation
          ELSE
             CALL ROT_YZ(C%CHART%ang_out(1),X,C%MAGP%P%BETA0,OU,k%TIME)                ! rotations
             CALL ROT_XZ(C%CHART%ang_out(2),X,C%MAGP%P%BETA0,OU,k%TIME)
             CALL ROT_XY(C%CHART%ang_out(3),X)  !,OU)
             CALL TRANS(C%CHART%d_out,X,C%MAGP%P%BETA0,OU,k%TIME)                       !translation
          ENDIF
       ELSE
          IF(ENTERING) THEN
             C%CHART%d_out(1)=-C%CHART%d_out(1)
             C%CHART%d_out(2)=-C%CHART%d_out(2)
             C%CHART%ang_out(3)=-C%CHART%ang_out(3)
             CALL TRANS(C%CHART%d_out,X,C%MAGP%P%BETA0,OU,k%TIME)                       !translation
             CALL ROT_XY(C%CHART%ang_out(3),X)  !,OU)
             CALL ROT_XZ(C%CHART%ang_out(2),X,C%MAGP%P%BETA0,OU,k%TIME)
             CALL ROT_YZ(C%CHART%ang_out(1),X,C%MAGP%P%BETA0,OU,k%TIME)                ! rotations
             C%CHART%d_out(1)=-C%CHART%d_out(1)
             C%CHART%d_out(2)=-C%CHART%d_out(2)
             C%CHART%ang_out(3)=-C%CHART%ang_out(3)
          ELSE
             C%CHART%d_in(1)=-C%CHART%d_in(1)
             C%CHART%d_in(2)=-C%CHART%d_in(2)
             C%CHART%ang_in(3)=-C%CHART%ang_in(3)
             CALL TRANS(C%CHART%d_in,X,C%MAGP%P%BETA0,OU,k%TIME)                       !translation
             CALL ROT_XY(C%CHART%ang_in(3),X)  !,OU)
             CALL ROT_XZ(C%CHART%ang_in(2),X,C%MAGP%P%BETA0,OU,k%TIME)
             CALL ROT_YZ(C%CHART%ang_in(1),X,C%MAGP%P%BETA0,OU,k%TIME)                ! rotations
             C%CHART%d_in(1)=-C%CHART%d_in(1)
             C%CHART%d_in(2)=-C%CHART%d_in(2)
             C%CHART%ang_in(3)=-C%CHART%ang_in(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_FIBP



END MODULE S_TRACKING
