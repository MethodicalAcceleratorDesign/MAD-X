!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_pol_sagan
  use precision_constants
  implicit none
  public
  private bLPOL2_0
  integer,private,parameter::n_max=10

  TYPE POL_sagan
     INTEGER ia(n_max)
     real(dp) Sa(n_max)
  END TYPE POL_sagan


  TYPE POL_BLOCK_sagan
     INTEGER Iinternal(2)
     real(dp) SInternal(2)
     type(POL_sagan) w
  END TYPE POL_BLOCK_sagan


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE bLPOL2_0
  END INTERFACE




contains




  SUBROUTINE  bLPOL2_0(S2,S1)
    implicit none
    type (POL_BLOCK_sagan),INTENT(inOUT):: S2
    INTEGER,INTENT(IN):: S1

    S2%SInternal=one
    S2%Iinternal=0
    S2%w%Sa=one
    S2%w%ia=0

  END SUBROUTINE bLPOL2_0


end module S_pol_sagan
