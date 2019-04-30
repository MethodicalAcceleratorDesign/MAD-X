!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module definition
  !  use define_newda
  !  use precision_constants   ! added to replace use define_newda
  !  use scratch_size
  !  use DABNEW
  !  use my_own_1D_TPSA
  use lielib_yang_berz, junk_no=>no,junk_nd=>nd,junk_nd2=>nd2,junk_ndpt=>ndpt,junk_nv=>nv
  use c_dabnew
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
  integer master,c_master
  !  integer,parameter::lnv=100
  !  scratch variables
  INTEGER iassdoluser(ndumt)
  integer DUMMY  !,temp
  integer iass0user(ndumt)
  integer c_DUMMY !,c_temp
  integer c_iass0user(ndumt)
  INTEGER c_iassdoluser(ndumt)
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
  INTEGER,PARAMETER::NMAX=22
  integer,private,parameter::n_max=10   ! sagan stuff
  INTEGER, PARAMETER :: CASE1=1,CASE2=2, CASE0=0, CASEP1=-1,CASEP2=-2
  INTEGER, PARAMETER :: CASET=3,CASETF1=4,CASETF2=5
  INTEGER,PARAMETER  :: ISPIN0R=1,ISPIN1R=3
  logical(lp) :: doing_ac_modulation_in_ptc=.false.
  integer, target :: nb_ =0   ! global group index
  integer, parameter :: ndim2t=10   ! maximum complex size
  integer, parameter :: wiggler_suntao=24
  real(dp) :: global_e =0
  integer :: bmadparser = 0
  integer,parameter :: nacmax = 3
  logical :: tangent = .false.,force_rescale=.false.   ! force_rescale for vorname=HELICAL see fibre_work routine
  TYPE sub_taylor
     INTEGER j(lnv)
     INTEGER min,max
  END TYPE sub_taylor

  !!&1
  TYPE taylor
     INTEGER I !@1  integer I is a pointer in old da-package of Berz
     !     type (taylorlow) j !@1   Taylorlow is an experimental type not supported
  END TYPE taylor
  type(taylor) temp
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
 !    integer g,nb  !  group index, number in group
     !&2
  END TYPE REAL_8

 type  quaternion
  real(dp) x(0:3)
 end type  quaternion  

 type  complex_quaternion
  complex(dp) x(0:3)
 end type  complex_quaternion  

 type  quaternion_8 
  type(real_8) x(0:3)
END TYPE quaternion_8


  !@3 ---------------------------------------------</br>
  TYPE complex_8
     type (complextaylor) t
     complex(dp) r
     logical(lp) alloc
     integer kind
     integer i,j
     complex(dp) s
  !   integer g,nb  !  group index
  END TYPE complex_8

  type(taylor)        varf1,varf2
  type(complextaylor) varc1,varc2

  !Radiation
 ! TYPE ENV_8
 !    type (REAL_8) v
 !    type (REAL_8) e(ndim2)
 !    type (REAL_8) sigma0(ndim2)
 !    type (REAL_8) sigmaf(ndim2)
 ! END TYPE ENV_8

  type spinor
         real(dp) x(3)  ! x(3) = (s_x, s_y, s_z)   with  |s|=1   
  end type spinor

  type spinor_8
     type(real_8) x(3)  ! x(3) = (s_x, s_y, s_z)   with  |s|=1
  end type spinor_8

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
  !
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

!  TYPE beamenvelope   !@2 A kind of Normal Form for Radiative Envelope
     ! radiation normalization
!     type (damap) transpose    ! Transpose of map which acts on polynomials
!     type (taylor) bij         !  Represents the stochastic kick at the end of the turn  Env_f=M Env_f M^t + B
!     TYPE (pbresonance) bijnr   !  Equilibrium beam sizes in resonance basis
!     real(dp) s_ij0(6,6)  !  equilibrium beam sizes
!     type (taylor) sij0  !  equilibrium beam sizes
!     real(dp) emittance(3),tune(3),damping(3)
!     logical(lp) AUTO,STOCHASTIC
!     real(dp)  KICK(3)   ! fake kicks for tracking stochastically
!     type (damap) STOCH  ! Diagonalized of stochastic part of map for tracking stochastically
!  END TYPE beamenvelope

  !@3 ---------------------------------------------</br>
  type  tree_element   !@1  USED FOR FAST TRACKING IN O_TREE_ELEMENT.F90
  !   character(204) , pointer :: file
     real(dp) ,  DIMENSION(:), POINTER :: CC
     real(dp) ,  DIMENSION(:), POINTER :: fixr,fix,fix0
     integer,  DIMENSION(:), POINTER :: JL,JV
     INTEGER,POINTER :: N,NP,no
     real(dp), pointer :: e_ij(:,:)
     real(dp), pointer :: rad(:,:)
     real(dp), pointer :: ds,beta0,eps
     logical, pointer :: symptrack,usenonsymp,factored
  end  type tree_element
  !@3 ---------------------------------------------</br>
 
  !@3 ---------------------------------------------</br>
 
  !@3 ---------------------------------------------</br>

  !@3 ---------------------------------------------</br>
  include "a_def_frame_patch_chart.inc"
  include "a_def_sagan.inc"
  include "a_def_element_fibre_layout.inc"
  include "a_def_all_kind.inc"


  !@3 ---------------------------------------------</br>
  type(fibre), pointer :: lost_fibre=>null()
  type(integration_node), pointer :: lost_node=>null()

  type rf_phasor
     real(dp) x(2)
     real(dp) om
     real(dp) t
  end type rf_phasor
  !@3 ---------------------------------------------</br>
  type rf_phasor_8
     type(real_8)  x(2)  ! The two hands of the clock
     type(real_8) om     ! the omega of the modulation
     real(dp) t          ! the pseudo-time
  end type rf_phasor_8
  !@3 ---------------------------------------------</br>
  type probe
     real(dp) x(6)
     type(spinor) s(3)
     type(quaternion) q
     type(rf_phasor)  AC(nacmax)
     integer:: nac=0
     logical u,use_q
     type(integration_node),pointer :: last_node=>null()
      real(dp) e
  end type probe
  !@3 ---------------------------------------------</br>
type probe_8
   type(real_8) x(6)     ! Polymorphic orbital ray
   type(spinor_8) s(3)   ! Polymorphic spin s(1:3)
   type(quaternion_8) q
   type(rf_phasor_8)  ac(nacmax)  ! Modulation of magnet
   integer:: nac=0 !  number of modulated clocks <=nacmax
   real(dp) E_ij(6,6)   !  Envelope for stochastic radiation
   real(dp) x0(6) ! initial value of the ray for TPSA calculations with c_damap
   !   stuff for exception
   logical u,use_q
   type(integration_node),pointer :: last_node=>null()
  real(dp) e
end type probe_8
  !@3 ---------------------------------------------</br>
  type TEMPORAL_PROBE
     TYPE(probe)  XS   ! probe at r=0
     TYPE(INTEGRATION_NODE), POINTER :: NODE
     real(DP)  r, dt0   !  penetration ration, penetration time
     real(DP)   POS(6),T  ! (x,y,z,px,py,pz) at dt0 and total time
     real(DP)   IC(3)  ! (x,y,z,px,py,pz) at dt0
     type(spinor) s(3) ! spin vectors at dt0
  END type TEMPORAL_PROBE
  !@3 ---------------------------------------------</br>
  type TEMPORAL_BEAM
     TYPE(TEMPORAL_PROBE), pointer :: TP(:)
     real(DP) a(3),ent(3,3),p0c,total_time
     integer n
     type(integration_node),pointer :: c   ! pointer close to a(3)
     type(internal_state)  state
  END type TEMPORAL_BEAM
  !@3 ---------------------------------------------</br>
  TYPE C_taylor
     INTEGER I !@1  integer I is a pointer to the complexified Berz package
  END TYPE C_taylor
  type(c_taylor),pointer :: dz_c(:)=>null()
  type(real_8),pointer :: dz_8(:)=>null()
  type(taylor),pointer :: dz_t(:)=>null()
  !@3 ---------------------------------------------</br>
  type c_dascratch
     type(c_taylor), pointer :: t
     TYPE (c_dascratch),POINTER :: PREVIOUS
     TYPE (c_dascratch),POINTER :: NEXT
  end type c_dascratch
  !@3 ---------------------------------------------</br>
  TYPE c_dalevel
     INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN
     !
     logical(lp),POINTER ::CLOSED
     TYPE (c_dascratch), POINTER :: PRESENT
     TYPE (c_dascratch), POINTER :: END
     TYPE (c_dascratch), POINTER :: START
     TYPE (c_dascratch), POINTER :: START_GROUND ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
     TYPE (c_dascratch), POINTER :: END_GROUND ! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
  END TYPE c_dalevel
  !@3 ---------------------------------------------</br>
  type c_spinmatrix
     type(c_taylor) s(3,3)
  end type c_spinmatrix
  !@3 ---------------------------------------------</br>
  type c_spinor
     type(c_taylor) v(3)
  end type c_spinor
  !@3 ---------------------------------------------</br>

type c_yu_w
 type (c_taylor),pointer :: w(:,:)=> null() !@1 orbital part of the map 
 integer :: n=0 !@1 of non zero w
end type c_yu_w

 type  c_quaternion
  type(c_taylor) x(0:3)
END TYPE c_quaternion

type c_damap
 type (c_taylor) v(lnv) !@1 orbital part of the map 
 integer :: n=0 !@1 number of plane allocated
 type(c_spinmatrix) s !@1 spin matrix
 type(c_quaternion) q
 complex(dp) e_ij(6,6) !@1 stochastic fluctuation in radiation theory
 complex(dp) x0(lnv) 
 logical :: tpsa=.false.
end type c_damap

  !@3 ---------------------------------------------</br>
  TYPE c_vector_field  !@1
   !@1 n dimension used v(1:n) (nd2 by default) ; nrmax some big integer if eps<1  
   integer :: n=0,nrmax
   !@1 if eps=-integer  then |eps| # of Lie brackets are taken 
   !@ otherwise eps=eps_tpsalie=10^-9
   real(dp) eps
   type (c_taylor) v(lnv)  
 !  type(c_spinor) om 
   type(c_quaternion) q
  END TYPE c_vector_field
  !@3 ---------------------------------------------</br>
  TYPE c_vector_field_fourier  !@1 
      integer :: n=0
      type (c_vector_field), pointer :: f(:)  =>null()                       
  END TYPE c_vector_field_fourier
  !@3 ---------------------------------------------</br>
  TYPE c_factored_lie
      integer :: n = 0   
      integer :: dir= 0     
       type (c_vector_field), pointer :: f(:)=>null()                   
  END TYPE c_factored_lie
  !@3 ---------------------------------------------</br>
TYPE c_normal_form
 type(c_damap) a1   !@1 brings to fix point at least linear
 type(c_damap) a2   !@1 linear normal form 
 type(c_factored_lie) g   !@1 nonlinear part of a in phasors
 type(c_factored_lie) ker !@1  kernel i.e. normal form in phasors
 type(c_damap) a_t !@1 transformation a (m=a n a^-1) 
 type(c_damap) n   !@1 transformation n (m=a n a^-1)      
 type(c_damap) As  !@1  For Spin   (m = As a n a^-1 As^-1)  
 type(c_damap) Atot  !@1  For Spin   (m = Atot n Atot^-1)  
 integer NRES,M(NDIM2t/2,NRESO),ms(NRESO) !@1 stores resonances to be left in the map, including spin (ms)
 real(dp) tune(NDIM2t/2),damping(NDIM2t/2),spin_tune !@1 Stores simple information
 logical positive ! forces positive tunes (close to 1 if <0)
!!!Envelope radiation stuff to normalise radiation (Sand's like theory)
 complex(dp) s_ij0(6,6)  !@1  equilibrium beam sizes
 complex(dp) s_ijr(6,6)  !@1  equilibrium beam sizes in resonance basis
 real(dp) emittance(3)   !@1  Equilibrium emittances as defined by Chao (computed from s_ijr(2*i-1,2*i) i=1,2,3 )
END TYPE c_normal_form
  !@2 the routine c_canonize(at,a_cs,a0,a1,a2,phase) factors neatly the map "at"
  !@2 at= a_cs o rotation(phase) where  a_cs = a0 o a1 o a2 ; this gives the phase advance even nonlinear!
  !@3 ---------------------------------------------</br>
type(c_taylor) c_temp

 TYPE c_ray
  complex(dp) x(lnv)            !# orbital and/or magnet modulation clocks
  complex(dp) s1(3),s2(3),s3(3) !# 3 spin directions
  type(complex_quaternion) q    !# quaternion
  integer n                     !# of dimensions used in x(lnv)
  complex(dp) x0(lnv)           !# the initial orbit around which the map is computed
 end type c_ray


TYPE fibre_array
   type(fibre), pointer :: p  => null()
   integer, pointer :: pos  => null()
   real(dp),pointer :: v=> null() , vmax=> null(); 
   real(dp), pointer :: s(:)=> null()
   real(dp), pointer :: err=> null()
END TYPE fibre_array


TYPE node_array
   type(integration_node), pointer :: t  => null()
   integer, pointer :: pos  => null()
   real(dp),pointer :: v=> null() , vmax=> null(); 
   complex(dp), pointer :: s(:)=> null()
   real(dp), pointer :: err=> null()
   type(c_vector_field), pointer :: f => null()
   type(c_damap), pointer :: m => null()
END TYPE node_array


contains







 subroutine alloc_fibre_array(a,n,m)
 implicit none
 type(fibre_array), allocatable :: a(:)
 integer i,n,m

 allocate(a(n))

 do i=1,n
   allocate(a(i)%pos)
   allocate(a(i)%v,a(i)%vmax,a(i)%err,a(i)%s(m))
   a(i)%s=0.0_dp
   a(i)%v=0.0_dp;a(i)%vmax=1.d38;
   a(i)%pos=0
 enddo

 end  subroutine alloc_fibre_array

 subroutine kill_fibre_array(a)
 implicit none
 type(fibre_array), allocatable :: a(:)
 integer i

 do i=1,size(a)
   deallocate(a(i)%pos)
   deallocate(a(i)%v,a(i)%vmax,a(i)%err,a(i)%s)
 enddo

 end  subroutine kill_fibre_array

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
       messagelost="h_definition.f90 root : negative argument "
    ELSEIF(X>=0.0_dp) THEN
       ROOT=SQRT(X)
    ELSE      !  IF X IS NOT A NUMBER
       ROOT=1.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="h_definition.f90 root : NaN argument "
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
       messagelost="h_definition.f90 arcsin : abs(x)>1 "
    ELSEIF(ABS(X)<=1.0_dp) THEN
       ARCSIN=ASIN(X)
    ELSE      !  IF X IS NOT A NUMBER
       ARCSIN=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="h_definition.f90 arcsin : x is NaN "
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
       messagelost="h_definition.f90 arccos : abs(x)>1 "
    ELSEIF(ABS(X)<=1.0_dp) THEN
       ARCCOS=ACOS(X)
    ELSE      !  IF X IS NOT A NUMBER
       ARCCOS=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="h_definition.f90 arccos : x is NaN "
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
       messagelost="h_definition.f90 loge : negative argument "
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
       messagelost="h_definition.f90 coseh : abs(x)>hyperbolic_aperture "
    ELSEIF(ABS(X)<=c_%hyperbolic_aperture) THEN
       COSEH=COSH(X)
    ELSE      !  IF X IS NOT A NUMBER
       COSEH=1.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="h_definition.f90 coseh : x is NaN "
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
       messagelost="h_definition.f90 sineh : abs(x)>hyperbolic_aperture "
    ELSEIF(ABS(X)<=c_%hyperbolic_aperture) THEN
       SINEH=SINH(X)
    ELSE      !  IF X IS NOT A NUMBER
       SINEH=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="h_definition.f90 sineh : x is NaN "
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
       messagelost="h_definition.f90 arctan : abs(x)>hyperbolic_aperture "
    ELSEIF(ABS(X)<=c_%hyperbolic_aperture) THEN
       arctan=atan(X)
    ELSE      !  IF X IS NOT A NUMBER
       arctan=0.0_dp
       c_%CHECK_STABLE=.FALSE.
       messagelost="h_definition.f90 arctan : x is NaN "
    ENDIF

  END FUNCTION arctan



end module definition

