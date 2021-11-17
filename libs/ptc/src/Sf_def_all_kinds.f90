!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_def_all_kinds
  use S_status
  implicit none
  public
 


 

 

 

 

contains

  !  RECURSIVE
  SUBROUTINE GET_LENGTH(R,L)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    REAL(DP), INTENT(OUT) :: L
    TYPE(FIBRE), POINTER:: P

    INTEGER I
    P=>R%START
    L=0.0_dp
    DO I=1,R%N
       IF(P%MAG%KIND/=KIND23) THEN
          L=L+P%MAG%P%LD
          !       ELSE
          !          CALL GET_LENGTH(P%MAG%G23,LG)
          !          L=L+LG
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE GET_LENGTH





end module S_def_all_kinds
