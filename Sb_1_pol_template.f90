!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90
module S_pol_user1
  use precision_constants
  implicit none
  public
  private bLPOL1_0

  TYPE POL_BLOCK1
     INTEGER Iinternal
     real(dp) SInternal
  END TYPE POL_BLOCK1


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE bLPOL1_0
  END INTERFACE




CONTAINS

  SUBROUTINE  bLPOL1_0(S2,S1)
    implicit none
    type (POL_BLOCK1),INTENT(inOUT):: S2
    INTEGER,INTENT(IN):: S1
    S2%SInternal=one
    S2%Iinternal=0
  END SUBROUTINE bLPOL1_0


end module S_pol_user1
