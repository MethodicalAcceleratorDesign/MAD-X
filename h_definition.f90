!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size

module definition
  use define_newda
  use scratch_size
  use DABNEW
  use lielib_berz, junk_no=>no,junk_nd=>nd,junk_nd2=>nd2,junk_ndpt=>ndpt,junk_nv=>nv
  use newda
  USE LIELIB_ETIENNE
  implicit none
  logical(lp) :: newread=.false. ,newprint =  .false. , first_time = .true.
  logical(lp) :: print77=.true. ,read77 =  .true.
  logical(lp) :: no_ndum_check = .false.
  logical(lp),TARGET :: setknob = .false., knob=.true.
  logical(lp),TARGET :: insane_PTC = .false.
  complex(dp) i_
  integer master
  integer,parameter::lnv=100
  !  scratch variables
  INTEGER iassdoluser(ndumt)
  integer DUMMY,temp
  integer DUMuser(ndumt,ndummax)
  integer iass0user(ndumt)
  integer,parameter::ndim2=2*ndim
  integer,parameter::mmmmmm1=1,mmmmmm2=2,mmmmmm3=3,mmmmmm4=4
  type (taylorlow) DUMluser(ndumt,ndummax)
  type (taylorlow) DUMMYl,templ             !,DUMl(ndum)
  private NDC,NDC2,NDT,IREF,itu,idpr,iflow,jtune,nres,ifilt
  private nplane,idsta,ista
  private xintex,dsta,sta,angle,rad,ps,rads,mx

  TYPE taylor
     ! integer is a pointer in old da-package of Berz
     INTEGER I
     type (taylorlow) j ! Taylorlow is the newda Taylor series
  END TYPE taylor

  TYPE UNIVERSAL_TAYLOR
     INTEGER, POINTER:: N,NV    !  Number of coeeficients and number of variables
     REAL(DP), POINTER,dimension(:)::C  ! Coefficients C(N)
     INTEGER, POINTER,dimension(:,:)::J ! Exponents of each coefficients J(N,NV)
  END TYPE UNIVERSAL_TAYLOR

  TYPE complextaylor
     type (taylor) r
     type (taylor) i
  END TYPE complextaylor

  ! this is a real polymorphic type

  TYPE REAL_8
     TYPE (TAYLOR) T   ! IF TAYLOR
     REAL(DP) R ! IF REAL
     INTEGER KIND  ! 0,1,2,3 (1=REAL,2=TAYLOR,3=TAYLOR KNOB, 0=SPECIAL)
     INTEGER I   ! USED FOR KNOBS AND KIND=0
     REAL(DP) S   ! SCALING FOR KNOBS AND KIND=0
     LOGICAL(lp) :: ALLOC  ! IF TAYLOR IS ALLOCATED IN DA-PACKAGE
  END TYPE REAL_8

  ! this is a real polymorphic type

  TYPE double_complex
     type (complextaylor) t
     complex(dp) r
     logical(lp) alloc
     integer kind
     integer i,j
     real(dp) s
  END TYPE double_complex

  type(taylor)        varf1,varf2
  type(complextaylor) varc1,varc2

  !Radiation
  TYPE ENV_8
     type (REAL_8) v
     type (REAL_8) e(ndim2)
     type (REAL_8) sigma0(ndim2)
     type (REAL_8) sigmaf(ndim2)
  END TYPE ENV_8

  type dummapping
     integer  d(ndim2)
     TYPE (TAYLORLOW) e(ndim2)
  end type dummapping




end module definition
