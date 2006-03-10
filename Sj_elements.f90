!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90

MODULE S_ELEMENTS
  USE S_DEF_ELEMENT
  IMPLICIT NONE
  public
  !  PRIVATE TRACKR,TRACKP,TRACKS
  private SURVEY_mag

  INTERFACE TRACK
     !     MODULE PROCEDURE TRACKR
     !     MODULE PROCEDURE TRACKP
     !     MODULE PROCEDURE TRACKS
     MODULE PROCEDURE SURVEY_mag  ! tracks a chart for survey
  END INTERFACE



contains

  !  SUBROUTINE TRACKR(EL,X,MID)
  !    IMPLICIT NONE
  !    real(dp),INTENT(INOUT):: X(6)
  !    TYPE(ELEMENT),INTENT(INOUT):: EL
  !    TYPE(WORM),OPTIONAL, INTENT(INOUT):: MID
  !
  !    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
  !    SELECT CASE(EL%KIND)
  !    CASE(KIND0)
  !       IF(PRESENT(MID)) CALL XMID(MID,X,0)
  !    case(KIND1)
  !       CALL TRACK(EL%D0,X,MID)
  !    case(KIND2)
  !       CALL TRACK(EL%K2,X,MID)
  !    case(KIND3)
  !       CALL TRACK(EL%K3,X,MID)
  !    case(KIND4)
  !       CALL TRACK(EL%C4,X,MID)
  !    case(KIND5)
  !       CALL TRACK(EL%S5,X,MID)
  !    case(KIND6)
  !       CALL TRACK(EL%T6,X,MID)
  !    case(KIND7)
  !       CALL TRACK(EL%T7,X,MID)
  !    case(KIND8)
  !       CALL TRACK(EL%S8,X,MID)
  !    case(KIND9)
  !       CALL TRACK(EL%S9,X,MID)
  !    case(KIND10)
  !       CALL TRACK(EL%TP10,X,MID)
  !    CASE(KIND11:KIND14)
  !       call TRACK(EL%MON14,X,MID)
  !    CASE(KIND15)
  !       call TRACK(EL%SEP15,X,MID)
  !    CASE(KIND16,KIND20)
  !       call TRACK(EL%K16,X,MID)
  !    CASE(KIND17)
  !       call TRACK(EL%S17,X,MID)
  !    CASE(KIND18)
  !       call TRACK(EL%RCOL18,X,MID)
  !    CASE(KIND19)
  !       call TRACK(EL%ECOL19,X,MID)
  !    CASE(KIND21)
  !       call TRACK(EL%CAV21,X,MID)
  !    CASE(KIND22)
  !       call TRACK(EL%M22,X,MID)
  !    case(KINDFITTED)
  !       call track(el%bend,x,MID)
  !    case(KINDUSER1)
  !       call TRACK(EL%U1,X,MID)
  !    case(KINDUSER2)
  !       call TRACK(EL%U2,X,MID)
  !    case default
  !       w_p=0
  !       w_p%nc=1
  !       w_p%fc='(1((1X,a72)))'
  !       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKR"
  !       CALL WRITE_E(0)
  !    END SELECT
  !  END SUBROUTINE TRACKR
  !
  !  SUBROUTINE TRACKP(EL,X,MID)
  !    IMPLICIT NONE
  !    TYPE(REAL_8),INTENT(INOUT):: X(6)
  !    TYPE(ELEMENTP),INTENT(INOUT):: EL
  !    TYPE(WORM_8),OPTIONAL, INTENT(INOUT):: MID
  !
  !    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
  !    SELECT CASE(EL%KIND)
  !    CASE(KIND0)
  !       IF(PRESENT(MID)) CALL XMID(MID,X,0)
  !    case(KIND1)
  !       CALL TRACK(EL%D0,X,MID)
  !    case(KIND2)
  !       CALL TRACK(EL%K2,X)
  !    case(KIND3)
  !       CALL TRACK(EL%K3,X,MID)
  !    case(KIND4)
  !       CALL TRACK(EL%C4,X,MID)
  !    case(KIND5)
  !       CALL TRACK(EL%S5,X,MID)
  !    case(KIND6)
  !       CALL TRACK(EL%T6,X,MID)
  !    case(KIND7)
  !       CALL TRACK(EL%T7,X,MID)
  !    case(KIND8)
  !       CALL TRACK(EL%S8,X,MID)
  !    case(KIND9)
  !       CALL TRACK(EL%S9,X,MID)
  !    case(KIND10)
  !       CALL TRACK(EL%TP10,X,MID)
  !    CASE(KIND11:KIND14)
  !       call TRACK(EL%MON14,X,MID)
  !    CASE(KIND15)
  !       call TRACK(EL%SEP15,X,MID)
  !    CASE(KIND16,KIND20)
  !       call TRACK(EL%K16,X,MID)
  !    CASE(KIND17)
  !       call TRACK(EL%S17,X,MID)
  !    CASE(KIND18)
  !       call TRACK(EL%RCOL18,X,MID)
  !    CASE(KIND19)
  !       call TRACK(EL%ECOL19,X,MID)
  !    CASE(KIND21)
  !       call TRACK(EL%CAV21,X,MID)
  !    CASE(KIND22)
  !       call TRACK(EL%M22,X,MID)
  !    case(KINDFITTED)
  !       call tracK(el%bend,x,MID)
  !    case(KINDUSER1)
  !       call TRACK(EL%U1,X,MID)
  !    case(KINDUSER2)
  !       call TRACK(EL%U2,X,MID)
  !    case default
  !       w_p=0
  !       w_p%nc=1
  !       w_p%fc='(1((1X,a72)))'
  !       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKP"
  !       CALL WRITE_E(0)
  !    END SELECT
  !  END SUBROUTINE TRACKP
  !
  !  SUBROUTINE TRACKS(EL,X,MID)
  !    IMPLICIT NONE
  !    TYPE(ENV_8),INTENT(INOUT):: X(6)
  !    TYPE(ELEMENTP),INTENT(INOUT):: EL
  !    TYPE(INNER_ENV_8_DATA),OPTIONAL, INTENT(INOUT):: MID
  !
  !    if(associated(el%p%aperture)) call CHECK_APERTURE(EL%p%aperture,X)
  !    SELECT CASE(EL%KIND)
  !    CASE(KIND0)
  !       !       IF(PRESENT(MID)) CALL XMID(MID,X,0)
  !    case(KIND1)
  !       CALL TRACK(EL%D0,X)
  !    case(KIND2)
  !       CALL TRACK(EL%K2,X)
  !    case(KIND3)
  !       CALL TRACK(EL%K3,X)
  !    case(KIND4)
  !       CALL TRACK(EL%C4,X)
  !    case(KIND5)
  !       CALL TRACK(EL%S5,X)
  !    case(KIND6)
  !       CALL TRACK(EL%T6,X)
  !    case(KIND7)
  !       CALL TRACK(EL%T7,X)
  !    case(KIND8)
  !       CALL TRACK(EL%S8,X)
  !    case(KIND9)
  !       CALL TRACK(EL%S9,X)
  !    case(KIND10)
  !       CALL TRACK(EL%TP10,X)
  !    CASE(KIND11:KIND14)
  !       call TRACK(EL%MON14,X)
  !    CASE(KIND15)
  !       call TRACK(EL%SEP15,X)
  !    CASE(KIND16,KIND20)
  !       call TRACK(EL%K16,X)
  !    CASE(KIND17)
  !       call TRACK(EL%S17,X)
  !    CASE(KIND18)
  !       call TRACK(EL%RCOL18,X)
  !    CASE(KIND19)
  !       call TRACK(EL%ECOL19,X)
  !    CASE(KIND21)
  !       call TRACK(EL%CAV21,X)
  !    CASE(KIND22)
  !       call TRACK(EL%M22,X)
  !    case(KINDUSER1)
  !       call TRACK(EL%U1,X)
  !    case(KINDUSER2)
  !       call TRACK(EL%U2,X)
  !    case(KINDFITTED)
  !       w_p=0
  !       w_p%nc=1
  !       w_p%fc='(1((1X,a72)))'
  !       w_p%c(1)= "KINDFITTED not supported "
  !       CALL WRITE_E(0)
  !    case default
  !       w_p=0
  !       w_p%nc=1
  !       w_p%fc='(1((1X,a72)))'
  !       write(w_p%c(1),'(1x,i4,a21)') el%kind," not supported TRACKS"
  !       CALL WRITE_E(0)
  !    END SELECT
  !
  !
  !
  !  END SUBROUTINE TRACKS



  SUBROUTINE SURVEY_mag(C,el,dir,magnetframe,e_in) !  Tracks the chart through a magnet
    IMPLICIT NONE
    TYPE(CHART), TARGET ,optional, INTENT(INOUT):: C
    type(element),intent(inout) :: el
    TYPE(magnet_frame), OPTIONAL :: magnetframe
    INTEGER, intent(in):: dir
    TYPE(INNER_FRAME), OPTIONAL :: E_IN

    !  All PTC magnet have the same convention for the internal frame
    !  Show a user want to add a magnet with corckscrew survey
    !  his survey would have to be "caught" in this interface

    SELECT CASE(EL%KIND)
    case(kind0:KINDWIGGLER,kindmu)
       call SURVEY_chart(C,el%p,dir,magnetframe,E_IN)
    case default

       write(6,*) el%kind," not supported SURVEY_mag in S_ELEMENTS"

    END SELECT


  end SUBROUTINE SURVEY_mag



END  MODULE S_ELEMENTS
