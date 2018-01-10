!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

TYPE POL_sagan
   INTEGER ia(n_max)
   real(dp) Sa(n_max)
END TYPE POL_sagan


TYPE POL_BLOCK_sagan
   INTEGER Iinternal(6)
   real(dp) SInternal(6)
   type(POL_sagan) w
END TYPE POL_BLOCK_sagan

TYPE undu_R
   real(dp), DIMENSION(:,:), pointer :: K=> null()
   real(dp), DIMENSION(:), pointer :: A=> null()
   real(dp), DIMENSION(:), pointer :: F=> null()
   real(dp), DIMENSION(:), pointer :: x0=> null()
   real(dp), DIMENSION(:), pointer :: y0=> null()
   real(dp), DIMENSION(:,:), pointer :: KE=> null()
   real(dp), DIMENSION(:), pointer :: AE=> null()
   real(dp), DIMENSION(:), pointer :: FE=> null()
   real(dp), DIMENSION(:), pointer :: x0E=> null()
   real(dp), DIMENSION(:), pointer :: y0E=> null()
   real(dp), DIMENSION(:), pointer :: ex=> null(),ey=> null()
   real(dp), pointer :: offset=> null()
      integer, pointer :: n=> null(),ne=> null()
   integer, DIMENSION(:), pointer :: FORM=> null()
      integer, DIMENSION(:), pointer :: FORME=> null()
END TYPE undu_R

TYPE undu_p
   TYPE(REAL_8), DIMENSION(:,:), pointer :: K=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: A=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: F=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: x0=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: y0=> null()
   TYPE(REAL_8), DIMENSION(:,:), pointer :: KE=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: AE=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: FE=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: x0E=> null()
   TYPE(REAL_8), DIMENSION(:), pointer :: y0E=> null()
   real(dp), DIMENSION(:), pointer :: ex=> null(),ey=> null()
   TYPE(REAL_8), pointer :: offset
         integer, pointer :: n=> null(),ne=> null()
   integer, DIMENSION(:), pointer :: FORM=> null()
      integer, DIMENSION(:), pointer :: FORME=> null()
END TYPE undu_p

TYPE SAGAN
   TYPE(MAGNET_CHART), POINTER :: P=> null()
   real(dp), POINTER ::L=> null()      !  MUST ALWAYS BE THERE
   integer, pointer :: n_min=> null()
   real(dp),  DIMENSION(:), POINTER :: AN=> null(),BN=> null()         !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   real(dp),  DIMENSION(:),  POINTER ::INTERNAL  => null()     ! INTERNAL IS AN EXAMPLE
   LOGICAL(lp), POINTER :: xprime  => null()
   TYPE(undu_R), POINTER :: W=> null()
END  TYPE SAGAN

TYPE SAGANP
   TYPE(MAGNET_CHART), POINTER :: P=> null()
   TYPE(REAL_8), POINTER ::L=> null()       !  MUST ALWAYS BE THERE
   integer, pointer :: n_min=> null()
   TYPE(REAL_8),  DIMENSION(:), POINTER :: AN=> null(),BN => null()        !Multipole component (OPTIONAL)
   !   ADD INTERNAL STUFF HERE AS POINTERS
   TYPE(REAL_8), DIMENSION(:),POINTER ::INTERNAL=> null()        ! INTERNAL IS AN EXAMPLE
   LOGICAL(lp), POINTER :: xprime  => null()
   TYPE(undu_p), POINTER :: W=> null()
END  TYPE SAGANP
