TYPE WIGGLER_R
   real(dp), DIMENSION(:,:), pointer :: K
   real(dp), DIMENSION(:), pointer :: A
   real(dp), DIMENSION(:), pointer :: F
   real(dp), pointer :: offset
   integer, DIMENSION(:), pointer :: FORM
END TYPE WIGGLER_R

TYPE WIGGLER_P
   TYPE(REAL_8), DIMENSION(:,:), pointer :: K
   TYPE(REAL_8), DIMENSION(:), pointer :: A
   TYPE(REAL_8), DIMENSION(:), pointer :: F
   TYPE(REAL_8), pointer :: offset
   integer, DIMENSION(:), pointer :: FORM
END TYPE WIGGLER_P

TYPE USER2
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp), POINTER ::L      !  MUST ALWAYS BE THERE
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   real(dp),  DIMENSION(:),  POINTER ::INTERNAL       ! INTERNAL IS AN EXAMPLE
   TYPE(WIGGLER_R), POINTER :: W
END  TYPE USER2

TYPE USER2P
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8), POINTER ::L       !  MUST ALWAYS BE THERE
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   TYPE(REAL_8), DIMENSION(:),POINTER ::INTERNAL        ! INTERNAL IS AN EXAMPLE
   TYPE(WIGGLER_P), POINTER :: W
END  TYPE USER2P
