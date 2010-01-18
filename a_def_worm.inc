!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

TYPE INNER_FRAME
   INTEGER, POINTER ::   NST
   REAL(DP), POINTER :: ORIGIN(:,:)
   REAL(DP), POINTER :: FRAME(:,:,:)
   type(fibre), POINTER :: F
   REAL(DP),POINTER :: L(:)
   logical(LP),POINTER :: DO_SURVEY
END TYPE INNER_FRAME


TYPE worm
   INTEGER, POINTER ::   NST
   INTEGER, POINTER ::   POS(:)
   REAL(DP), POINTER :: RAY(:,:)
   TYPE(INNER_FRAME), POINTER :: E
   type(fibre), POINTER :: F
END TYPE worm

!TYPE WORM_8
!   INTEGER, POINTER ::   NST
!   INTEGER, POINTER ::   POS(:)
!   TYPE(REAL_8), POINTER :: RAY(:,:)
!   TYPE(INNER_FRAME), POINTER :: E
!   type(fibre), POINTER :: F
!END TYPE WORM_8

!TYPE INNER_ENV_8_DATA
!   INTEGER, POINTER ::   NST
!   INTEGER, POINTER ::   POS(:)
!   TYPE(ENV_8), POINTER :: RAY(:,:)
!   TYPE(INNER_FRAME), POINTER :: E
!   type(fibre), POINTER :: F
!END TYPE INNER_ENV_8_DATA
