!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

TYPE INNER_FRAME
   INTEGER, POINTER ::   NST           => null()
   REAL(DP), POINTER :: ORIGIN(:,:)    => null()
   REAL(DP), POINTER :: FRAME(:,:,:)   => null()
   type(fibre), POINTER :: F           => null()
   REAL(DP),POINTER :: L(:)            => null()
   logical(LP),POINTER :: DO_SURVEY    => null()
END TYPE INNER_FRAME


TYPE worm
   INTEGER, POINTER ::   NST           => null()
   INTEGER, POINTER ::   POS(:)        => null()
   REAL(DP), POINTER :: RAY(:,:)       => null()
   TYPE(INNER_FRAME), POINTER :: E     => null()
   type(fibre), POINTER :: F           => null()
END TYPE worm	                   

!TYPE WORM_8
!   INTEGER, POINTER ::   NST          => null()
!   INTEGER, POINTER ::   POS(:)       => null()
!   TYPE(REAL_8), POINTER :: RAY(:,:)  => null()
!   TYPE(INNER_FRAME), POINTER :: E    => null()
!   type(fibre), POINTER :: F          => null()
!END TYPE WORM_8	                   

!TYPE INNER_ENV_8_DATA
!   INTEGER, POINTER ::   NST          => null()
!   INTEGER, POINTER ::   POS(:)       => null()
!   TYPE(ENV_8), POINTER :: RAY(:,:)   => null()
!   TYPE(INNER_FRAME), POINTER :: E    => null()
!   type(fibre), POINTER :: F          => null()
!END TYPE INNER_ENV_8_DATA             => null()
