! user1

TYPE USER1
   TYPE(MAGNET_CHART), POINTER :: P
   real(dp), POINTER ::L      !  MUST ALWAYS BE THERE
   real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   !           .........
   real(dp), POINTER ::INTERNAL           ! INTERNAL IS AN EXAMPLE
END  TYPE USER1

TYPE USER1P
   TYPE(MAGNET_CHART), POINTER :: P
   TYPE(REAL_8), POINTER ::L                             !  MUST ALWAYS BE THERE
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component (OPTIONAL but always defined)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   TYPE(REAL_8), POINTER ::INTERNAL        ! INTERNAL is an example of a variable specific to user1p
END  TYPE USER1P

