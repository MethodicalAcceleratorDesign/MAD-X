!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

TYPE POL_sagan
   INTEGER ia(n_max)
   real(dp) Sa(n_max)
END TYPE POL_sagan


TYPE POL_BLOCK_sagan
   INTEGER Iinternal(2)
   real(dp) SInternal(2)
   type(POL_sagan) w
END TYPE POL_BLOCK_sagan

TYPE undu_R
   real(dp), DIMENSION(:,:), pointer :: K
   real(dp), DIMENSION(:), pointer :: A
   real(dp), DIMENSION(:), pointer :: F
   real(dp), pointer :: offset
   integer, DIMENSION(:), pointer :: FORM
END TYPE undu_R

TYPE undu_p
   TYPE(REAL_8), DIMENSION(:,:), pointer :: K
   TYPE(REAL_8), DIMENSION(:), pointer :: A
   TYPE(REAL_8), DIMENSION(:), pointer :: F
   TYPE(REAL_8), pointer :: offset
   integer, DIMENSION(:), pointer :: FORM
END TYPE undu_p

TYPE SAGAN
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp), POINTER ::L      !  MUST ALWAYS BE THERE
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   real(dp),  DIMENSION(:),  POINTER ::INTERNAL       ! INTERNAL IS AN EXAMPLE
   TYPE(undu_R), POINTER :: W
END  TYPE SAGAN

TYPE SAGANP
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8), POINTER ::L       !  MUST ALWAYS BE THERE
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   TYPE(REAL_8), DIMENSION(:),POINTER ::INTERNAL        ! INTERNAL IS AN EXAMPLE
   TYPE(undu_p), POINTER :: W
END  TYPE SAGANP
