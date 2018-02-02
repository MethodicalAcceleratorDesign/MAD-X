!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

TYPE AFFINE_FRAME
   REAL(DP), POINTER,DIMENSION(:)  ::   ANGLE(:) => null()
   REAL(DP), POINTER,DIMENSION(:)  ::   D(:)     => null()
   REAL(DP), POINTER,DIMENSION(:)  ::   A        => null()
   REAL(DP), POINTER,DIMENSION(:,:)::   ENT      => null()
   REAL(DP), POINTER,DIMENSION(:)  ::   B        => null()
   REAL(DP), POINTER,DIMENSION(:,:)::   EXI      => null()
END TYPE AFFINE_FRAME


TYPE MAGNET_FRAME
   REAL(DP), POINTER,DIMENSION(:)  ::   A   => null()
   REAL(DP), POINTER,DIMENSION(:,:)::   ENT => null()
   REAL(DP), POINTER,DIMENSION(:)  ::   O   => null()
   REAL(DP), POINTER,DIMENSION(:,:)::   MID => null()
   REAL(DP), POINTER,DIMENSION(:)  ::   B   => null()
   REAL(DP), POINTER,DIMENSION(:,:)::   EXI => null() 
END TYPE MAGNET_FRAME

TYPE PATCH
   INTEGER(2), POINTER:: PATCH => null()     ! IF TRUE, SPACIAL PATCHES NEEDED
   INTEGER, POINTER :: A_X1 => null() ,A_X2 => null()    ! FOR ROTATION OF PI AT ENTRANCE = -1, DEFAULT = 1 ,
   INTEGER, POINTER :: B_X1 => null() ,B_X2 => null()    ! FOR ROTATION OF PI AT EXIT = -1    , DEFAULT = 1
   REAL(DP),DIMENSION(:), POINTER:: A_D => null() ,B_D => null()       !ENTRACE AND EXIT TRANSLATIONS  A_D(3)
   REAL(DP),DIMENSION(:), POINTER:: A_ANG => null() ,B_ANG => null()    !ENTRACE AND EXIT ROTATIONS    A_ANG(3)
   INTEGER(2), POINTER:: ENERGY => null()    ! IF TRUE, ENERGY PATCHES NEEDED
   INTEGER(2), POINTER:: TIME => null()      ! IF TRUE, TIME PATCHES NEEDED
   REAL(DP), POINTER:: A_T => null() ,B_T => null()      ! TIME SHIFT NEEDED SOMETIMES WHEN RELATIVE TIME IS USED
   REAL(DP), POINTER:: A_L => null() ,B_L => null() 
   REAL(DP), POINTER:: p0b => null() ,b0b => null()      ! for bmad patching
END TYPE PATCH

TYPE CHART
   type(magnet_frame), pointer :: f => null() 
   !  FIBRE MISALIGNMENTS
   real(dp),dimension(:),  POINTER::   D_IN => null() ,ANG_IN => null() 
   real(dp),dimension(:),  POINTER::   D_OUT => null() ,ANG_OUT => null() 
END TYPE CHART
