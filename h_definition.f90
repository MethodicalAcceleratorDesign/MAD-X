!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module definition
  !  use define_newda
  use precision_constants   ! added to replace use define_newda
  use scratch_size
  use DABNEW
  use lielib_berz, junk_no=>no,junk_nd=>nd,junk_nd2=>nd2,junk_ndpt=>ndpt,junk_nv=>nv
  !  use newda
  !  USE LIELIB_ETIENNE
  implicit none
  public
  logical(lp) :: newread=.false. ,newprint =  .false. , first_time = .true.
  logical(lp) :: print77=.true. ,read77 =  .true.
  logical(lp) :: no_ndum_check = .false.
  logical(lp),TARGET :: insane_PTC = .false.
  logical(lp),TARGET :: setknob = .false. !@1 Real part of knobs cannot set
  logical(lp),TARGET :: knob=.true. !@1 Knobs are effective
  integer, target :: npara_fpp !@1 position of last non-parameter tpsa variable
  !  complex(dp), parameter :: i_ = = cmplx(zero,one,kind=dp)
  !  complex(dp)  :: i_ = (zero,one)    ! cmplx(zero,one,kind=dp)
  complex(dp), parameter :: i_ = ( 0.0_dp,1.0_dp )    ! cmplx(zero,one,kind=dp)
  integer master
  integer,parameter::lnv=100
  !  scratch variables
  INTEGER iassdoluser(ndumt)
  integer DUMMY,temp
  integer iass0user(ndumt)
  integer,parameter::ndim2=2*ndim
  integer,parameter::mmmmmm1=1,mmmmmm2=2,mmmmmm3=3,mmmmmm4=4
  !  type (taylorlow) DUMMYl,templ             !,DUMl(ndum)
  !  private NDC,NDC2,NDT,IREF,itu,iflow,jtune,nres,ifilt  ,idpr
  !  private nplane,idsta,ista
  !  private xintex,dsta,sta,angle,rad,ps,rads,mx
  ! numerical differentiation by knobs
  logical(lp) :: knob_numerical=.false.
  real(dp) ::  knob_eps(lnv)=c_1d_6
  integer ::  knob_i =0

  !
  TYPE sub_taylor
     INTEGER j(lnv)
     INTEGER min,max
  END TYPE sub_taylor

  !!&1
  TYPE taylor
     INTEGER I !@1  integer I is a pointer in old da-package of Berz
     !     type (taylorlow) j !@1   Taylorlow is an experimental type not supported
  END TYPE taylor
  !@2  UNIVERSAL_TAYLOR is used by Sagan in BMAD Code at Cornell
  !@2  Also used by MAD-XP
  TYPE UNIVERSAL_TAYLOR
     INTEGER, POINTER:: N,NV    !  Number of coeeficients and number of variables
     REAL(DP), POINTER,dimension(:)::C  ! Coefficients C(N)
     INTEGER, POINTER,dimension(:,:)::J ! Exponents of each coefficients J(N,NV)
  END TYPE UNIVERSAL_TAYLOR
  !@3 ---------------------------------------------</br>
  TYPE complextaylor
     type (taylor) r  !@1 Real part
     type (taylor) i  !@1 Imaginary part
  END TYPE complextaylor

  !&2
  ! this is a real polymorphic type

  TYPE REAL_8
     TYPE (TAYLOR) T      !@1  USED IF TAYLOR
     REAL(DP) R           !@1    USED IF REAL
     !&2
     INTEGER KIND  !@1  0,1,2,3 (1=REAL,2=TAYLOR,3=TAYLOR KNOB, 0=SPECIAL)
     INTEGER I   !@1   USED FOR KNOBS AND SPECIAL KIND=0
     REAL(DP) S   !@1   SCALING FOR KNOBS AND SPECIAL KIND=0
     LOGICAL(lp) :: ALLOC  !@1 IF TAYLOR IS ALLOCATED IN DA-PACKAGE
     !&2
  END TYPE REAL_8
  !@3 ---------------------------------------------</br>
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
     type (REAL_8) sigma0(ndim2)
     type (REAL_8) sigmaf(ndim2)
  END TYPE ENV_8

  type spinor
     real(dp) x(3)
  end type spinor

  type spinor_8
     type(real_8) x(3)
     !   type(real_8) m(3,3)
  end type spinor_8

  type probe
     real(dp) x(6)
     type(spinor) s
  end type probe

  type probe_8
     type(real_8) x(6)
     type(spinor_8) s
     !  type (REAL_8) E_ij(ndim2,ndim2)
     real(dp) E_ij(ndim2,ndim2)
  end type probe_8

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

  !&1
  TYPE DAMAP
     TYPE (TAYLOR) V(ndim2)    ! Ndim2=6 but allocated to nd2=2,4,6 ! etienne_oct_2004
  END TYPE DAMAP
  !@3 ---------------------------------------------</br>

  TYPE GMAP
     TYPE (TAYLOR) V(lnv)    ! Ndim2=6 but allocated to nd2=2,4,6 ! etienne_oct_2004
     integer N
  END TYPE GMAP

  TYPE vecfield
     type (taylor) v(ndim2)          !@1 <font face="Times New Roman">V<sub>i</sub>&#8706;<sub>i</sub></font> Operator
     integer ifac                    !@1 Type of Factorization 0,1,-1 (One exponent, Dragt-Finn, Reversed Dragt-Finn)
  END TYPE vecfield
  !@3 ---------------------------------------------</br>
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

  !&5
  TYPE DRAGTFINN
     real(dp)  constant(ndim2)
     type (damap) Linear
     type (vecfield) nonlinear
     type (pbfield)  pb
  END TYPE DRAGTFINN
  !@3 ---------------------------------------------</br>
  TYPE reversedragtfinn
     real(dp)  CONSTANT(NDIM2)
     type (damap) Linear
     type (vecfield) nonlinear
     type (pbfield)  pb
  END TYPE reversedragtfinn
  !@3 ---------------------------------------------</br>
  TYPE ONELIEEXPONENT
     real(dp) EPS
     type (vecfield) VECTOR
     type (pbfield)  pb
  END TYPE ONELIEEXPONENT

  !&5

  TYPE normalform
     type (damap) A_t   ! Total A  :  A_t= A1 o A_rest
     type (damap) A1    ! Dispersion
     type (reversedragtfinn) A  ! Linear A_t and nonlinear A_t
     type (dragtfinn) NORMAL    ! Normal is the Normal Form R
     type (damap) DHDJ  ! Contains the tunes in convenient form: extracted from NORMAL (=R)
     !         .
     !         .
     !         .
     !         .
     !&3
     real(dp) TUNE(NDIM),DAMPING(NDIM)  ! linear tune and linear damping
     integer nord,jtune                 ! nord=1 A1 first order in parameters
     integer NRES,M(NDIM,NRESO),PLANE(NDIM)
     logical(lp) AUTO
     !&3
  END TYPE normalform

  TYPE genfield
     type (taylor) h
     type (damap) m
     type (taylor) d(ndim,ndim)
     type (damap) linear
     type (damap) lineart
     type (damap) mt
     real(dp) constant(ndim2),eps
     integer imax     !@1 imax=Maximum Number of Iteration (default=1000)
     integer ifac     !@1 ifac = the map is raised to the power 1/ifac and iterated ifac times (default=1)
     logical(lp) linear_in     !@1 Linear part is left in the map  (default=.false.)
     integer no_cut   !@1 Original map is not symplectic on and above no_cut
  END TYPE genfield



  TYPE pbresonance
     type (pbfield)  cos,sin
     integer ifac
  END TYPE pbresonance
  !@3 ---------------------------------------------</br>
  TYPE vecresonance
     type (vecfield)  cos,sin
     integer ifac
  END TYPE vecresonance
  !@3 ---------------------------------------------</br>
  TYPE taylorresonance
     type (taylor)  cos,sin
  END TYPE taylorresonance

  TYPE beamenvelope   !@2 A kind of Normal Form for Radiative Envelope
     ! radiation normalization
     type (damap) transpose    ! Transpose of map which acts on polynomials
     type (taylor) bij         !  Represents the stochastic kick at the end of the turn  Env_f=M Env_f M^t + B
     TYPE (pbresonance) bijnr   !  Equilibrium beam sizes in resonance basis
     type (taylor) sij0  !  equilibrium beam sizes
     real(dp) emittance(3),tune(3),damping(3)
     logical(lp) AUTO,STOCHASTIC
     real(dp)  KICK(3)   ! fake kicks for tracking stochastically
     type (damap) STOCH  ! Diagonalized of stochastic part of map for tracking stochastically
  END TYPE beamenvelope


  type  tree_element   !@1  USED FOR FAST TRACKING IN O_TREE_ELEMENT.F90
     real(dp) ,  DIMENSION(:), POINTER :: CC
     integer,  DIMENSION(:), POINTER :: JL,JV
     INTEGER,POINTER :: N,ND2,no
  end  type tree_element

  type daspin
     REAL(DP) X(6)
     type(damap) M
     type(real_8) s(3,3)
  end type daspin

  type normal_spin
     type(normalform) N
     type(real_8) NS(3,3)
     type(real_8) AS(3,3)
     real(dp) tune
     integer NRES,M(NDIM)

  end type normal_spin


end module definition
