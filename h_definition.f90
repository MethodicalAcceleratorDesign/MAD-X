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
  integer iass0user(ndumt)
  integer,parameter::ndim2=2*ndim
  integer,parameter::mmmmmm1=1,mmmmmm2=2,mmmmmm3=3,mmmmmm4=4
  type (taylorlow) DUMMYl,templ             !,DUMl(ndum)
  private NDC,NDC2,NDT,IREF,itu,iflow,jtune,nres,ifilt  ,idpr
  private nplane,idsta,ista
  private xintex,dsta,sta,angle,rad,ps,rads,mx
  ! numerical differentiation by knobs
  logical(lp) :: knob_numerical=.false.
  real(dp) ::  knob_eps(lnv)=1.e-6_dp
  integer ::  knob_i =0

  !
  TYPE taylor
     INTEGER I !@1  integer I is a pointer in old da-package of Berz
     type (taylorlow) j !@1   Taylorlow is an experimental type not supported
  END TYPE taylor

  TYPE UNIVERSAL_TAYLOR
     INTEGER, POINTER:: N,NV    !  Number of coeeficients and number of variables
     REAL(DP), POINTER,dimension(:)::C  ! Coefficients C(N)
     INTEGER, POINTER,dimension(:,:)::J ! Exponents of each coefficients J(N,NV)
  END TYPE UNIVERSAL_TAYLOR

  TYPE complextaylor
     type (taylor) r  !@1 Real part
     type (taylor) i  !@1 Imaginary part
  END TYPE complextaylor

  ! this is a real polymorphic type

  TYPE REAL_8
     TYPE (TAYLOR) T   !@1  USED IF TAYLOR
     REAL(DP) R !@1    USED IF REAL
     INTEGER KIND  !@1  0,1,2,3 (1=REAL,2=TAYLOR,3=TAYLOR KNOB, 0=SPECIAL)
     INTEGER I   !@1   USED FOR KNOBS AND SPECIAL KIND=0
     REAL(DP) S   !@1   SCALING FOR KNOBS AND SPECIAL KIND=0
     LOGICAL(lp) :: ALLOC  !@1 IF TAYLOR IS ALLOCATED IN DA-PACKAGE
  END TYPE REAL_8

  ! this is a real polymorphic type

  TYPE double_complex
     type (complextaylor) t
     complex(dp) r
     logical(lp) alloc
     integer kind
     integer i,j
     complex(dp) s
  END TYPE double_complex

  type(taylor)        varf1,varf2
  type(complextaylor) varc1,varc2

  !Radiation
  TYPE ENV_8
     type (REAL_8) v
     type (REAL_8) e(ndim2)
     type (REAL_8) sigma0(ndim2)  !@2 added by yukiko nogiwa
     type (REAL_8) sigmaf(ndim2)
  END TYPE ENV_8

  !    scratch levels of DA using linked list

  type dascratch
     type(taylor), pointer :: t
     TYPE (dascratch),POINTER :: PREVIOUS
     TYPE (dascratch),POINTER :: NEXT
  end type dascratch

  TYPE dalevel
     INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN
     !
     logical(lp),POINTER ::CLOSED
     TYPE (dascratch), POINTER :: PRESENT
     TYPE (dascratch), POINTER :: END
     TYPE (dascratch), POINTER :: START
     TYPE (dascratch), POINTER :: START_GROUND ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
     TYPE (dascratch), POINTER :: END_GROUND ! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
  END TYPE dalevel

  TYPE DAMAP
     TYPE (TAYLOR) V(NDIM2)    ! Ndim2=6 but allocated to nd2=2,4,6
  END TYPE DAMAP

  TYPE vecfield
     type (taylor) v(ndim2)    !@1 <font face="Times New Roman">V<sub>i</sub>&#8706;<sub>i</sub></font> Operator
     integer ifac              !@1 Type of Factorization 0,1,-1 (One exponent, Dragt-Finn, Reversed Dragt-Finn)
  END TYPE vecfield

  TYPE pbfield
     type (taylor) h
     integer ifac
  END TYPE pbfield



  TYPE tree
     type (taylor) branch(ndim2)
  END TYPE tree

  !Radiation
  TYPE radtaylor
     type (taylor) v
     type (taylor) e(ndim2)
  END TYPE radtaylor

  TYPE DRAGTFINN
     real(dp)  constant(ndim2)
     type (damap) Linear
     type (vecfield) nonlinear
     type (pbfield)  pb
  END TYPE DRAGTFINN

  TYPE ONELIEEXPONENT
     !     real(dp)  CONSTANT(NDIM2),
     real(dp) EPS
     type (vecfield) VECTOR
     type (pbfield)  pb
  END TYPE ONELIEEXPONENT

  TYPE reversedragtfinn
     real(dp)  CONSTANT(NDIM2)
     type (damap) Linear
     type (vecfield) nonlinear
     type (pbfield)  pb
  END TYPE reversedragtfinn

  TYPE normalform
     type (damap) A_t
     type (damap) A1
     type (reversedragtfinn) A
     type (dragtfinn) NORMAL
     type (damap) DHDJ
     real(dp) TUNE(NDIM),DAMPING(NDIM)
     integer nord,jtune
     integer NRES,M(NDIM,NRESO),PLANE(NDIM)
     logical(lp) AUTO
  END TYPE normalform

  TYPE genfield
     type (taylor) h
     type (damap) m
     type (taylor) d(ndim,ndim)
     type (damap) linear
     type (damap) lineart
     type (damap) mt
     real(dp) constant(ndim2)
     integer ifac,imax
     logical(lp) linear_in
  END TYPE genfield



  TYPE pbresonance
     type (pbfield)  cos,sin
     integer ifac
  END TYPE pbresonance

  TYPE vecresonance
     type (vecfield)  cos,sin
     integer ifac
  END TYPE vecresonance

  TYPE taylorresonance
     type (taylor)  cos,sin
  END TYPE taylorresonance


  TYPE beamenvelope
     ! radiation normalization
     type (damap) transpose    ! Transpose of map which acts on polynomials
     type (taylor) bij         !  Represents the stochastic kick at the end of the turn  Env_f=M Env_f M^t + B
     TYPE (pbresonance) bijnr   !  Equilibrium beam sizes in resonance basis
     type (taylor) sij0  !  equilibrium beam sizes
     real(dp) emittance(3),tune(3),damping(3)
     logical(lp) AUTO,STOCHASTIC
     real(dp)  KICK(3)
     type (damap) STOCH
  END TYPE beamenvelope


  type  tree_element
     real(dp) ,  DIMENSION(:), POINTER :: CC
     integer,  DIMENSION(:), POINTER :: JL,JV
     INTEGER,POINTER :: N,ND2
  end  type tree_element


end module definition
