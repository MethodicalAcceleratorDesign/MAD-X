!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
module S_pol_user2
  use precision_constants
  implicit none
  private bLPOL2_0


  TYPE POL_BLOCK2
     INTEGER Iinternal
     real(dp) SInternal
  END TYPE POL_BLOCK2


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE bLPOL2_0
  END INTERFACE




contains




  SUBROUTINE  bLPOL2_0(S2,S1)
    implicit none
    type (POL_BLOCK2),INTENT(inOUT):: S2
    INTEGER,INTENT(IN):: S1

    S2%SInternal=one
    S2%Iinternal=0

  END SUBROUTINE bLPOL2_0


end module S_pol_user2
