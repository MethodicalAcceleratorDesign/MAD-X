!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
MODULE S_ELEMENTS
  USE S_DEF_ELEMENT
  IMPLICIT NONE
  PRIVATE TRACKR,TRACKP,TRACKS

  INTERFACE TRACK
     MODULE PROCEDURE TRACKR
     MODULE PROCEDURE TRACKP
     MODULE PROCEDURE TRACKS
  END INTERFACE


contains

  SUBROUTINE TRACKR(EL,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(ELEMENT),INTENT(INOUT):: EL

    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
    SELECT CASE(EL%KIND)
    CASE(KIND0)
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
    CASE(KIND16)
       call TRACK(EL%K16,X)
    CASE(KIND17)
       call TRACK(EL%S17,X)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X)
    case(KINDFITTED)
       call track(el%bend,x)
    case(KINDUSER1)
       call TRACK(EL%U1,X)
    case(KINDUSER2)
       call TRACK(EL%U2,X)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKR"
       CALL WRITE_E(0)
    END SELECT
  END SUBROUTINE TRACKR

  SUBROUTINE TRACKP(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(ELEMENTP),INTENT(INOUT):: EL

    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
    SELECT CASE(EL%KIND)
    CASE(KIND0)
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
    CASE(KIND16)
       call TRACK(EL%K16,X)
    CASE(KIND17)
       call TRACK(EL%S17,X)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X)
    case(KINDFITTED)
       call tracK(el%bend,x)
    case(KINDUSER1)
       call TRACK(EL%U1,X)
    case(KINDUSER2)
       call TRACK(EL%U2,X)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKP"
       CALL WRITE_E(0)
    END SELECT
  END SUBROUTINE TRACKP

  SUBROUTINE TRACKS(EL,X)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: X(6)
    TYPE(ELEMENTP),INTENT(INOUT):: EL

    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
    SELECT CASE(EL%KIND)
    CASE(KIND0)
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
    CASE(KIND16)
       call TRACK(EL%K16,X)
    CASE(KIND17)
       call TRACK(EL%S17,X)
    CASE(KIND18)
       call TRACK(EL%RCOL18,X)
    CASE(KIND19)
       call TRACK(EL%ECOL19,X)
    case(KINDUSER1)
       call TRACK(EL%U1,X)
    case(KINDUSER2)
       call TRACK(EL%U2,X)
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







END  MODULE S_ELEMENTS
