!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module definition
  !  use define_newda
  !  use precision_constants   ! added to replace use define_newda
  !  use scratch_size
  !  use DABNEW
  !  use my_own_1D_TPSA
  use lielib_yang_berz, junk_no=>no,junk_nd=>nd,junk_nd2=>nd2,junk_ndpt=>ndpt,junk_nv=>nv
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
  !  integer,parameter::lnv=100
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
  real(dp) ::  knob_eps(lnv)=1e-6_dp
  integer ::  knob_i =0
  INTEGER,PARAMETER::NMAX=20
  integer,private,parameter::n_max=10   ! sagan stuff
  INTEGER, PARAMETER :: CASE1=1,CASE2=2, CASE0=0, CASEP1=-1,CASEP2=-2
  INTEGER, PARAMETER :: CASET=3,CASETF1=4,CASETF2=5
  INTEGER,PARAMETER  :: ISPIN0R=1,ISPIN1R=3
  logical(lp) :: doing_ac_modulation_in_ptc=.false.
  integer, target :: nb_ =0   ! global group index
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
     integer g,nb  !  group index, number in group
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
     integer g,nb  !  group index
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
     !  real(dp) G
  end type spinor

  type spinor_8
     type(real_8) x(3)
  end type spinor_8

  type res_spinor_8
     type(double_complex) x(3)
  end type res_spinor_8

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
  !&1
  !@3 ---------------------------------------------</br>

  TYPE GMAP
     TYPE (TAYLOR) V(lnv)    !
     integer N
  END TYPE GMAP

  !&4
  TYPE vecfield
     type (taylor) v(ndim2)          !@1 <font face="Times New Roman">V<sub>i</sub>&#8706;<sub>i</sub></font> Operator
     integer ifac                    !@1 Type of Factorization 0,1,-1 (One exponent, Dragt-Finn, Reversed Dragt-Finn)
  END TYPE vecfield
  !@3 ---------------------------------------------</br>
  TYPE pbfield
     type (taylor) h
     integer ifac
     integer nd_used
  END TYPE pbfield
  !&4


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
  !@3 ---------------------------------------------</br>

  TYPE normalform
     type (damap) A_t   ! Total A  :  A_t= A1 o A_rest
     type (damap) A1    !@1 Dispersion
     type (reversedragtfinn) A  !@1 Linear A_t and nonlinear A_t
     type (dragtfinn) NORMAL    !@1 Normal is the Normal Form R
     type (damap) DHDJ  !@1 Contains the tunes in convenient form: extracted from NORMAL (=R)
     real(dp) TUNE(NDIM),DAMPING(NDIM)  !@1 linear tune and linear damping
     integer nord,jtune                 !@1 nord=1 A1 first order in parameters
     integer NRES,M(NDIM,NRESO),PLANE(NDIM) !@1 NRES,M(NDIM,NRESO) -> resonances left in the map
     logical(lp) AUTO
  END TYPE normalform
  !@3 ---------------------------------------------</br>
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

  !@3 ---------------------------------------------</br>


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
     real(dp) s_ij0(6,6)  !  equilibrium beam sizes
     type (taylor) sij0  !  equilibrium beam sizes
     real(dp) emittance(3),tune(3),damping(3)
     logical(lp) AUTO,STOCHASTIC
     real(dp)  KICK(3)   ! fake kicks for tracking stochastically
     type (damap) STOCH  ! Diagonalized of stochastic part of map for tracking stochastically
  END TYPE beamenvelope


  type  tree_element   !@1  USED FOR FAST TRACKING IN O_TREE_ELEMENT.F90
     real(dp) ,  DIMENSION(:), POINTER :: CC
     real(dp) ,  DIMENSION(:), POINTER :: fix
     integer,  DIMENSION(:), POINTER :: JL,JV
     INTEGER,POINTER :: N,ND2,no
  end  type tree_element

  type spinmatrix
     type(real_8) s(3,3)
  end type spinmatrix

  type damapspin
     type(damap) M
     type(spinmatrix) s
     !     type(real_8) s(3,3)
     real(dp) e_ij(6,6)
  end type damapspin

  type normal_spin
     type(normalform) N   ! regular orbital normal form
     type(damapspin) a1   ! brings to fixed point
     type(damapspin) ar   ! normalises around the fixed point
     type(damapspin) as   ! pure spin map
     type(damapspin) a_t  ! !! (a_t%m,a_t%s) = (a1%m, I ) o (I ,as%s) o (ar%m,I)
!!!  extra spin info
     integer M(NDIM,NRESO),MS(NRESO),NRES  ! orbital and spin resonances to be left in the map
     type(real_8) n0(3)     ! n0 vector
     type(real_8) theta0    !  angle for the matrix around the orbit (analogous to linear tunes)
     real(dp) nu    !  spin tune
!!!Envelope radiation stuff
     real(dp) s_ij0(6,6)  !  equilibrium beam sizes
     complex(dp) s_ijr(6,6)  !  equilibrium beam sizes in resonance basis
     real(dp) emittance(3),tune(3),damping(3)   ! equilibrium emittances (partially well defined only for infinitesimal damping)
     logical(lp) AUTO,STOCHASTIC
     real(dp)  KICK(3)   ! fake kicks for tracking stochastically
     real(dp)  STOCH(6,6)  ! Diagonalized of stochastic part of map for tracking stochastically
     real(dp)  STOCH_inv(6,6)  ! Diagonalized of stochastic part of map for tracking stochastically
  end type normal_spin

  include "a_def_frame_patch_chart.inc"
  include "a_def_all_kind.inc"
  include "a_def_sagan.inc"
  include "a_def_element_fibre_layout.inc"

  type(fibre), pointer :: lost_fibre
  type(integration_node), pointer :: lost_node

  type rf_phasor
     real(dp) x(2)
     real(dp) om
     real(dp) t
  end type rf_phasor

  type rf_phasor_8
     type(real_8)  x(2)
     type(real_8) om
     reaL(DP) t
  end type rf_phasor_8

  type probe
     real(dp) x(6)
     type(spinor) s(ISPIN0R:ISPIN1R)
     type(rf_phasor) AC
     logical u
     type(integration_node),pointer :: lost_node
  end type probe

  type probe_8
     type(real_8) x(6)     ! Polymorphic orbital ray
     type(spinor_8) s(ISPIN0R:ISPIN1R)   ! Polymorphic spin
     type(rf_phasor_8) AC  ! Modulation
     real(dp) E_ij(6,6)   !  Envelope
     !   stuff for exception
     logical u
     type(integration_node),pointer :: lost_node
  end type probe_8

  type TEMPORAL_PROBE
     TYPE(probe)  XS
     TYPE(INTEGRATION_NODE), POINTER :: NODE
     real(DP)  DS,POS(6)
  END type TEMPORAL_PROBE

  type TEMPORAL_BEAM
     TYPE(TEMPORAL_PROBE), pointer :: TP(:)
     real(DP) a(3),ent(3,3),p0c,total_time
     integer n
     type(integration_node),pointer :: c   ! pointer close to a(3)
     type(internal_state)  state
  END type TEMPORAL_BEAM

contains

  SUBROUTINE RESET_APERTURE_FLAG(complete)
    IMPLICIT NONE
    logical(lp), optional :: complete
    logical(lp)  comp
    !    IF(c_%WATCH_USER) THEN
    !       IF(.NOT.ALLOW_TRACKING) THEN
    !          WRITE(6,*) "  EXECUTION OF THE CODE MUST BE INTERRUPTED AT YOUR REQUEST"
    !          WRITE(6,*) "  YOU DID NOT CHECK THE APERTURE STATUS"
    !          WRITE(6,*) "  USING A CALL TO PRODUCE_APERTURE_FLAG"
    !          WRITE(6,*) "  BEFORE CALLING A TRACKING FUNCTION"
    !          STOP 666
    !       ENDIF
    !       IF(.NOT.c_%check_stable) THEN
    !          WRITE(6,*) "  EXECUTION OF THE CODE MUST BE INTERRUPTED AT YOUR REQUEST"
    !          WRITE(6,*) " CODE MOST LIKELY DIED IN PURE DA/TPSA/LIE OPERATIONS  "
    !          STOP 667
    !       ENDIF
    !    ENDIF
    comp=.true.
    if(present(complete)) comp=complete
    c_%STABLE_DA =.TRUE.
    c_%CHECK_STABLE =.TRUE.
    c_%CHECK_MADX_APERTURE =.TRUE.
    c_%stable_da =.true.
    if(comp) then
       xlost=0.0_dp
       messagelost=" Aperture has been reset "
       nullify(lost_fibre)
       nullify(lost_node)

    endif
  END   SUBROUTINE RESET_APERTURE_FLAG

  SUBROUTINE PRODUCE_APERTURE_FLAG(I)
    IMPLICIT NONE
    INTEGER I
    I=0
    IF(.NOT.c_%CHECK_STABLE) THEN
       I=1
    ENDIF


    !    ALLOW_TRACKING=.TRUE.

  END   SUBROUTINE PRODUCE_APERTURE_FLAG

  ! moved here from sa_extend_poly.f90

  REAL(DP) FUNCTION  ROOT(X)  ! REPLACES SQRT(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       ROOT=1.0_dp
       return
    endif

    IF((X<0.0_dp).AND.c_%ROOT_CHECK) THEN
       ROOT=1.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Root undefined "
    ELSEIF(X>=0.0_dp) THEN
       ROOT=SQRT(X)
    ELSE      !  IF X IS NOT A NUMBER
       ROOT=1.0_dp
       c_%CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION ROOT

  REAL(DP) FUNCTION  ARCSIN(X)  ! REPLACES ASIN(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       ARCSIN=0.0_dp
       return
    endif
    IF((ABS(X)>1.0_dp).AND.c_%ROOT_CHECK) THEN
       ARCSIN=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Arcsin undefined "
    ELSEIF(ABS(X)<=1.0_dp) THEN
       ARCSIN=ASIN(X)
    ELSE      !  IF X IS NOT A NUMBER
       ARCSIN=0.0_dp
       c_%CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION ARCSIN

  REAL(DP) FUNCTION  ARCCOS(X)  ! REPLACES ACOS(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       ARCCOS=0.0_dp
       return
    endif
    IF((ABS(X)>1.0_dp).AND.c_%ROOT_CHECK) THEN
       ARCCOS=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Arccos undefined "
    ELSEIF(ABS(X)<=1.0_dp) THEN
       ARCCOS=ACOS(X)
    ELSE      !  IF X IS NOT A NUMBER
       ARCCOS=0.0_dp
       c_%CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION ARCCOS

  REAL(DP) FUNCTION  LOGE(X)  ! REPLACES ACOS(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       LOGE=0.0_dp
       return
    endif

    IF(X<=0.0_dp.AND.c_%ROOT_CHECK) THEN
       LOGE=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Log undefined "
    ELSE
       LOGE=LOG(X)
    ENDIF

  END FUNCTION LOGE



  REAL(DP) FUNCTION  COSEH(X) ! REPLACES COSH(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       COSEH=1.0_dp
       return
    endif

    IF((ABS(X)>c_%hyperbolic_aperture).AND.c_%ROOT_CHECK) THEN
       COSEH=1.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Coseh undefined "
    ELSEIF(ABS(X)<=c_%hyperbolic_aperture) THEN
       COSEH=COSH(X)
    ELSE      !  IF X IS NOT A NUMBER
       COSEH=1.0_dp
       c_%CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION COSEH

  REAL(DP) FUNCTION  SINEH(X) ! REPLACES SINH(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       SINEH=0.0_dp
       return
    endif

    IF((ABS(X)>c_%hyperbolic_aperture).AND.c_%ROOT_CHECK) THEN
       SINEH=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Sineh undefined "
    ELSEIF(ABS(X)<=c_%hyperbolic_aperture) THEN
       SINEH=SINH(X)
    ELSE      !  IF X IS NOT A NUMBER
       SINEH=0.0_dp
       c_%CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION SINEH

  REAL(DP) FUNCTION  arctan(X) ! REPLACES SINH(X)
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       arctan=0.0_dp
       return
    endif

    IF((ABS(X)>c_%hyperbolic_aperture).AND.c_%ROOT_CHECK) THEN
       arctan=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="Arctan undefined "
    ELSEIF(ABS(X)<=c_%hyperbolic_aperture) THEN
       arctan=atan(X)
    ELSE      !  IF X IS NOT A NUMBER
       arctan=0.0_dp
       c_%CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION arctan



end module definition
