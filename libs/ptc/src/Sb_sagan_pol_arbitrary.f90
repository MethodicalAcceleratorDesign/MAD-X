!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_pol_sagan
  !  use precision_constants
  use  definition
  implicit none
  public
  private bLPOL2_0

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE bLPOL2_0
  END INTERFACE




contains




  SUBROUTINE  bLPOL2_0(S2,S1)
    implicit none
    type (POL_BLOCK_sagan),INTENT(inOUT):: S2
    INTEGER,INTENT(IN):: S1

    S2%SInternal=1.0_dp
    S2%Iinternal=0
    S2%w%Sa=1.0_dp
    S2%w%ia=0

  END SUBROUTINE bLPOL2_0


end module S_pol_sagan
