module duan_zhe_map
implicit none
  private    !  this private can be removed outside PTC
  public track_TREE_probe_complex_zhe,probe,tree_element, read_tree_zhe,kill_tree_zhe,print_tree_elements_zhe
  public print,zhe_ini,track_TREE_probe_complex_ptc, dp,INTERNAL_STATE !,EQUAL_PROBE_REAL6_zhe,
  public DEFAULT0,TOTALPATH0 ,TIME0,ONLY_4d0,DELTA0,SPIN0,MODULATION0,only_2d0   ,nrmax_zhe
  public RADIATION0, NOCAVITY0, FRINGE0 ,STOCHASTIC0,ENVELOPE0,gaussian_seed_zhe,nrmax_used_zhe,alloc_bunch,kill_bunch
  PUBLIC track_TREE_probe_complex_ji,track_TREE_probe_complex_ji_symp,TRACK_TREE_PROBE_COMPLEX_JI_VEC
 public file_zhe,number_zhe_maps,get_seed,set_seed
 character(255) ::    file_zhe="zhe"
  integer ::  number_zhe_maps = 1
  public use_ji
  logical :: use_ji =.false.
  private alex_ISEED
  integer :: alex_ISEED=1000

     public CHECK_STABLE_ZHE
     public c_verbose_zhe 
 ! public EQUALi_zhe,EQUALt_zhe
  public OPERATOR(+),operator(-), assignment(=)
  real(kind(1d0)) :: doublenum = 0d0
  integer,parameter::lp=4
  integer,parameter::dp=selected_real_kind(2*precision(1.e0))
  logical(lp),parameter:: my_true=.true.
  logical(lp),parameter:: my_false=.false.
   LOGICAL(lp),TARGET  :: CHECK_STABLE_ZHE=.TRUE.
   logical :: c_verbose_zhe = .false.
 complex(dp), parameter :: i_ = ( 0.0_dp,1.0_dp )    ! cmplx(zero,one,kind=dp)
  integer,parameter::lno=200,lnv=100
  public zhe_ISEED
  integer :: zhe_ISEED=1000,nrmax_used, ntot=100
public change_ntot
private subq,unarysubq,addq,unaryADDq,absq,absq2,mulq,divq,ranf
private EQUALq,EQUALqr,EQUALqi,powq,printq ,invq
real(dp),parameter::pi=3.141592653589793238462643383279502e0_dp
real(dp) :: cut_zhe=6.0_dp
 logical :: use_quaternion = .false.
 public bunch
TYPE bunch
 type(probe), pointer :: xs(:)
 integer n,r,reloaded
 logical, pointer :: stable(:)
 real(dp) , pointer :: turn(:)
end TYPE bunch

TYPE INTERNAL_STATE
   INTEGER TOTALPATH   ! total time or path length is used
   LOGICAL(LP) TIME  ! Time is used instead of path length
   LOGICAL(LP) RADIATION ! Radiation is turned on
   LOGICAL(LP) NOCAVITY ! Cavity is turned into a drift
   LOGICAL(LP) FRINGE ! Fringe fields are turned on (mainly for quadrupoles)
   LOGICAL(LP) STOCHASTIC ! Random Stochastic kicks to x(5) 
   LOGICAL(LP) ENVELOPE ! Stochastic envelope terms tracked in probe_8 
   LOGICAL(LP) PARA_IN  ! If true, parameters in the map are included
   LOGICAL(LP) ONLY_4D  ! REAL_8 Taylor in (x,p_x,y,p_y)
   LOGICAL(LP) DELTA  ! REAL_8 Taylor in (x,p_x,y,p_y,delta)
   LOGICAL(LP) SPIN  ! Spin is tracked
   LOGICAL(LP) MODULATION !  One modulated family tracked by probe
   LOGICAL(LP) ONLY_2D  !  REAL_8 Taylor in (x,p_x)
   LOGICAL(LP) FULL_WAY  ! 
END TYPE INTERNAL_STATE



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
 !    integer, pointer :: ng
  end  type tree_element

  type spinor
         real(dp) x(3)  ! x(3) = (s_x, s_y, s_z)   with  |s|=1   
  end type spinor

  !@3 ---------------------------------------------</br>
 type  quaternion
  real(dp) x(0:3)
 end type  quaternion 
  !@3 ---------------------------------------------</br>
  type probe
     real(dp) x(6)
     type(spinor) s(3)
     type(quaternion) q
     logical u,use_q
  !   type(integration_node),pointer :: last_node=>null()
      real(dp) e
  end type probe


   integer :: nrmax = 1000
 




 
  !PRIVATE DTILTR,DTILTP,DTILTS


  LOGICAL(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.

  TYPE(INTERNAL_STATE),PARAMETER::DEFAULT0=INTERNAL_STATE   (0,f,f,f,f,f,f,f,f,f,f,f,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::TOTALPATH0=INTERNAL_STATE (1,f,f,f,f,f,f,f,f,f,f,f,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::TIME0=INTERNAL_STATE      (0,t,f,f,f,f,f,f,f,f,f,f,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::RADIATION0=INTERNAL_STATE (0,f,t,f,f,f,f,f,f,f,f,f,f,t)
  TYPE(INTERNAL_STATE),PARAMETER::NOCAVITY0=INTERNAL_STATE  (0,f,f,t,f,f,f,f,f,f,f,f,f,t)
  TYPE(INTERNAL_STATE),PARAMETER::FRINGE0=INTERNAL_STATE    (0,f,f,f,t,f,f,f,f,f,f,f,f,t)
  TYPE(INTERNAL_STATE),PARAMETER::STOCHASTIC0=INTERNAL_STATE(0,f,f,f,f,t,f,f,f,f,f,f,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::ENVELOPE0=INTERNAL_STATE  (0,f,f,f,f,f,t,f,f,f,f,f,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::ONLY_4d0=INTERNAL_STATE   (0,f,f,t,f,f,f,f,t,f,f,f,f,t)
  TYPE(INTERNAL_STATE),PARAMETER::DELTA0=INTERNAL_STATE     (0,f,f,t,f,f,f,f,t,t,f,f,f,t)
  TYPE(INTERNAL_STATE),PARAMETER::SPIN0=INTERNAL_STATE      (0,f,f,f,f,f,f,f,f,f,t,f,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::MODULATION0=INTERNAL_STATE(0,f,f,f,f,f,f,f,f,f,f,t,f,F)
  TYPE(INTERNAL_STATE),PARAMETER::only_2d0   =INTERNAL_STATE(0,f,f,t,f,f,f,f,f,f,f,f,t,t)




  !  private s_init,S_init_berz,MAKE_STATES_0,MAKE_STATES_m,print_s,CONV
  private CONV,print_probe_zhe
 
private lubksb_nr,ludcmp_nr,matinv
   integer,private,parameter::nmax=400
  real(dp),private,parameter::tiny=1e-20_dp
 
 
  CHARACTER(24) MYTYPE(-100:100)
 

private orthonormaliser
  private track_TREE_probe_complexr
  integer :: size_tree=15,size_ji=6+6+9+9+1,size_ji_vec=13
  integer :: ind_spin(3,3),k1_spin(9),k2_spin(9),ind_ji(3,3)

 

 

  INTERFACE OPERATOR (.min.)
     MODULE PROCEDURE minu_zhe                       ! to define the minus of Schmidt
  END INTERFACE

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUALt_zhe
     MODULE PROCEDURE EQUALi_zhe
     MODULE PROCEDURE EQUAL_PROBE_REAL6_zhe
     MODULE PROCEDURE EQUAL_PROBE_REAL6_bunch
     MODULE PROCEDURE EQUALq
     MODULE PROCEDURE EQUALqi
     MODULE PROCEDURE EQUALqr
  end  INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_zhe
     MODULE PROCEDURE PARA_REMA_zhe
     MODULE PROCEDURE unaryADDq  
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE sub_zhe
     MODULE PROCEDURE unarySUBq
  END INTERFACE


   INTERFACE OPERATOR (*)
     MODULE PROCEDURE mulq
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE divq
  END INTERFACE
 
  INTERFACE print
     MODULE PROCEDURE print_s
     MODULE PROCEDURE print_probe_zhe
  END INTERFACE
 
  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POWq
  END INTERFACE

  INTERFACE abs
     MODULE PROCEDURE absq  
  END INTERFACE


  INTERFACE abs_square
     MODULE PROCEDURE absq2
  END INTERFACE


  INTERFACE PRINT
     MODULE PROCEDURE PRINTQ
  END INTERFACE


  INTERFACE track_TREE_probe_complex_ptc
     MODULE PROCEDURE track_TREE_probe_complexr
  END INTERFACE 

contains 


  SUBROUTINE ALLOC_TREE(T,N,np)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T
    INTEGER , INTENT(IN) :: N,np
    integer i
    !IF(N==0) RETURN


    ALLOCATE(T%CC(N),T%fix0(6),T%fix(6),T%fixr(6),T%JL(N),T%JV(N),T%N,T%ds,T%beta0,T%np,T%no, & 
  !  t%e_ij(c_%nd2,c_%nd2),T%rad(c_%nd2,c_%nd2),t%usenonsymp, t%symptrack, t%eps)  !,t%file)
     t%e_ij(6,6),T%rad(6,6),t%usenonsymp, t%symptrack,t%factored, t%eps)  !,t%file)
    t%cc=0
    t%jl=0
    t%jv=0
    T%N=N
    T%np=np
    T%no=0
    T%fix=0.0_dp
    T%fix0=0.0_dp
    T%fixr=0.0_dp
    T%e_ij=0.0_dp
    T%ds=0.0_dp
    T%beta0=0.0_dp
    T%rad=0.0_dp
!    T%file=' '
    do i=1,6
     T%rad(i,i)=1.0_dp
    enddo
    t%eps=1.d-5
    t%symptrack=.false.
    t%usenonsymp=.false.
     t%factored=.false.
 !    t%ng=1
  END SUBROUTINE ALLOC_TREE
 

  SUBROUTINE KILL_TREE(T)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T


     IF(ASSOCIATED(T%CC))DEALLOCATE(T%CC,T%fix0,T%fix,T%fixr,t%ds,t%beta0,T%JL,T%JV,T%N,T%NP, &
    T%No,t%e_ij,t%rad,t%eps,t%symptrack,t%usenonsymp,t%factored  )  !,t%file)


  END SUBROUTINE KILL_TREE

   subroutine print_probe_zhe(DS,MFF)
    implicit none
    TYPE(probe), INTENT(INOUT) :: DS
    INTEGER I,mf
    INTEGER, optional :: MFf
 
    mf=6
    if(present(mff)) mf=mff
    WRITE(MF,*) " ORBIT "
    do i=1,6
       write(mf,*) ' Variable ',i
       write(mf,'(6(1X,G20.13))') ds%x(i) 
    enddo
  if(ds%use_q) then
   call print(ds%q,mf)
else
    WRITE(MF,*) " SPIN X "
       write(mf,'(3(1X,G20.13))') ds%s(1)%x 
 
    WRITE(MF,*) " SPIN Y "
       write(mf,'(3(1X,G20.13))') ds%s(2)%x 
 
    WRITE(MF,*) " SPIN Z "
       write(mf,'(3(1X,G20.13))') ds%s(3)%x 

 endif

  END subroutine print_probe_zhe

  

  subroutine EQUAL_PROBE_REAL6_zhe (P,X)
    implicit none
    TYPE(PROBE), INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: X(6)
    INTEGER I
    P%u=my_false
    !       P%s(0)%x=zero
    !       P%s(0)%x(n0_normal)=one

    DO I=1,3
       P%s(i)%x=0.0_dp
       P%s(i)%x(i)=1.0_dp
    enddo
    P%X=X
    p%q%x=0.0_dp
    p%q%x(0)=1.0_dp
     p%use_q=use_quaternion
  END    subroutine EQUAL_PROBE_REAL6_zhe



  subroutine EQUAL_PROBE_REAL6_bunch(P,X)
    implicit none
    TYPE(bunch), INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: X(6)
    INTEGER I
    do i=1,p%n
     p%xs(i)=x
    enddo
  END    subroutine EQUAL_PROBE_REAL6_bunch

  subroutine alloc_bunch(P,n)
    implicit none
    integer n
    TYPE(bunch), INTENT(INOUT) :: P
    p%n=n
    p%r=n

    p%reloaded=0
     allocate(p%xs(n),p%stable(n),p%turn(n))
     p%xs(:)%use_q=use_quaternion
     p%stable =.true.
     p%turn=0
 
  END    subroutine alloc_bunch

   subroutine kill_bunch(P)
    implicit none
    integer n
    TYPE(bunch), INTENT(INOUT) :: P
    deallocate(p%xs,p%turn,p%stable)
  END    subroutine kill_bunch

  FUNCTION minu_zhe( S1,S2  )
    implicit none
    logical(lp) minu_zhe
    logical(lp), INTENT (IN) :: S1
    logical(lp), INTENT (IN) :: S2

    minu_zhe=.false.
    if(s1.and.(.not.s2)) minu_zhe=my_true

  END FUNCTION minu_zhe




  


  subroutine print_s(S,MF)
    implicit none
    type (INTERNAL_STATE) S
    INTEGER MF


    write(mf,*) "************ State Summary ****************"
 
 

 
 
    write(mf, '((1X,a20,1x,i4))' )  "      TOTALPATH   = ", S%TOTALPATH
    !    write(mf, '((1X,a20,1x,a5))' )  "      EXACTMIS    = ", CONV(S%EXACTMIS    )
    write(mf,'((1X,a20,1x,a5))' ) "      RADIATION   = ", CONV(S%RADIATION  )
    write(mf,'((1X,a20,1x,a5))' ) "      STOCHASTIC  = ", CONV(S%STOCHASTIC  )
    write(mf,'((1X,a20,1x,a5))' ) "      ENVELOPE    = ", CONV(S%ENVELOPE  )
    write(mf,'((1X,a20,1x,a5))' ) "      NOCAVITY    = ", CONV(S%NOCAVITY )
    write(mf,'((1X,a20,1x,a5))' ) "      TIME        = ", CONV(S%TIME )
    write(mf,'((1X,a20,1x,a5))' ) "      FRINGE      = ", CONV(S%FRINGE   )
    write(mf,'((1X,a20,1x,a5))' ) "      PARA_IN     = ", CONV(S%PARA_IN  )
    write(mf,'((1X,a20,1x,a5))' ) "      ONLY_2D     = ", CONV(S%ONLY_2D   )
    write(mf,'((1X,a20,1x,a5))' ) "      ONLY_4D     = ", CONV(S%ONLY_4D   )
    write(mf,'((1X,a20,1x,a5))' ) "      DELTA       = ", CONV(S%DELTA    )
    write(mf,'((1X,a20,1x,a5))' ) "      SPIN        = ", CONV(S%SPIN    )
    write(mf,'((1X,a20,1x,a5))' ) "      MODULATION  = ", CONV(S%MODULATION    )

    !    write(mf,'((1X,a20,1x,I4))' ) " SPIN DIMENSION   = ", S%SPIN_DIM
    !   ! call ! WRITE_I
  end subroutine print_s

  FUNCTION CONV(LOG)
    IMPLICIT NONE
    CHARACTER(5) CONV
    logical(lp) LOG
    CONV="FALSE"
    IF(LOG) CONV="TRUE "
  END FUNCTION CONV

  



  SUBROUTINE  EQUALt_zhe(S2,S1)
    implicit none
    type (INTERNAL_STATE),INTENT(OUT)::S2
    type (INTERNAL_STATE),INTENT(IN)::S1

    S2%TOTALPATH=   S1%TOTALPATH
    !   S2%EXACTMIS=       S1%EXACTMIS
    S2%RADIATION=     S1%RADIATION
    S2%NOCAVITY=    S1%NOCAVITY
    S2%TIME=        S1%TIME
    S2%FRINGE=           S1%FRINGE
    S2%stochastic=           S1%stochastic
    S2%ENVELOPE=           S1%ENVELOPE
    S2%PARA_IN=     S1%PARA_IN
    S2%ONLY_2D=      S1%ONLY_2D
    S2%ONLY_4D=      S1%ONLY_4D
    S2%DELTA=       S1%DELTA
    S2%SPIN=       S1%SPIN
    S2%MODULATION=       S1%MODULATION
    S2%FULL_WAY=       S1%FULL_WAY
    !    S2%spin_dim=       S1%spin_dim
  END SUBROUTINE EQUALt_zhe





  SUBROUTINE  EQUALi_zhe(S2,i)
    implicit none
    type (INTERNAL_STATE),INTENT(OUT)::S2
    integer, intent(in) :: i
    
    S2=default0
    select case(i) 
     case(0)
      S2=default0
    case(1)
         S2=TOTALPATH0
    case(2)
         S2=TIME0
    case(3)
         S2=RADIATION0
    case(4)
         S2=NOCAVITY0
    case(5)
         S2=FRINGE0
    case(6)
         S2=STOCHASTIC0
    case(7)
         S2=ENVELOPE0
    case(9)
         S2=ONLY_4d0
    case(10)
         S2=DELTA0
    case(11)
         S2=SPIN0
    case(12)
         S2=MODULATION0
    case(13)
         S2=only_2d0
    case default
      S2%TOTALPATH = -1
    end select
 
  END SUBROUTINE EQUALi_zhe



  FUNCTION add_zhe( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) add_zhe
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2

    if(s2%totalpath/=0.and.s2%totalpath/=1) then 
      add_zhe=s1
      return
    endif
    if(s1%totalpath/=0.and.s1%totalpath/=1) then 
      add_zhe=s1
      return
    endif
    add_zhe%TOTALPATH=0
    if((S1%TOTALPATH==1).OR.(S2%TOTALPATH==1)) add_zhe%TOTALPATH=1

    !    add_zhe%EXACT    =       S1%EXACT.OR.S2%EXACT
    !   add_zhe%EXACTMIS    =       (S1%EXACTMIS.OR.S2%EXACTMIS).and.(.not.sixtrack_compatible)
    add_zhe%RADIATION  =  S1%RADIATION.OR.S2%RADIATION
    add_zhe%NOCAVITY =  S1%NOCAVITY.OR.S2%NOCAVITY
    add_zhe%TIME     =  S1%TIME.OR.S2%TIME
    add_zhe%FRINGE   =       S1%FRINGE.OR.S2%FRINGE
    add_zhe%stochastic   =       S1%stochastic.OR.S2%stochastic
    add_zhe%ENVELOPE   =       S1%ENVELOPE.OR.S2%ENVELOPE
    add_zhe%ONLY_2D  =       S1%ONLY_2D.OR.S2%ONLY_2D
    add_zhe%ONLY_4D  =       S1%ONLY_4D.OR.S2%ONLY_4D
    add_zhe%DELTA  =       S1%DELTA.OR.S2%DELTA
    add_zhe%SPIN  =       S1%SPIN.OR.S2%SPIN
    add_zhe%MODULATION  =       S1%MODULATION.OR.S2%MODULATION
    add_zhe%PARA_IN  =       S1%PARA_IN.OR.S2%PARA_IN
    !    add_zhe%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(add_zhe%stochastic) THEN
       add_zhe%RADIATION=T
    ENDIF
  !  IF(add_zhe%ENVELOPE) THEN
  !     add_zhe%radiation=T
  !  ENDIF
    IF(add_zhe%stochastic) THEN
       add_zhe%radiation=T
    ENDIF
    IF(add_zhe%DELTA) THEN
        add_zhe%ONLY_4D=T
       add_zhe%NOCAVITY =  T
    ENDIF
    IF(add_zhe%ONLY_4D) THEN
       add_zhe%TOTALPATH=  0
       add_zhe%RADIATION  =  F
       add_zhe%NOCAVITY =  T
       add_zhe%stochastic   =  F
       add_zhe%ENVELOPE   =  F
    ENDIF
    IF(add_zhe%ONLY_2D) THEN
       add_zhe%TOTALPATH=  0
       add_zhe%RADIATION  =  F
       add_zhe%NOCAVITY =  T
       add_zhe%stochastic   =  F
       add_zhe%ENVELOPE   =  F
    ENDIF
    if(add_zhe%only_4d.and.add_zhe%only_2d) add_zhe%only_4d=my_false

    add_zhe%RADIATION  =  S1%RADIATION.OR.S2%RADIATION
    add_zhe%NOCAVITY =  S1%NOCAVITY.OR.S2%NOCAVITY
    add_zhe%TIME     =  S1%TIME.OR.S2%TIME
    add_zhe%FRINGE   =       S1%FRINGE.OR.S2%FRINGE
    add_zhe%stochastic   =       S1%stochastic.OR.S2%stochastic
    add_zhe%ENVELOPE   =       S1%ENVELOPE.OR.S2%ENVELOPE
    add_zhe%ONLY_2D  =       S1%ONLY_2D.OR.S2%ONLY_2D
    add_zhe%ONLY_4D  =       S1%ONLY_4D.OR.S2%ONLY_4D
    add_zhe%DELTA  =       S1%DELTA.OR.S2%DELTA
    add_zhe%SPIN  =       S1%SPIN.OR.S2%SPIN
    add_zhe%MODULATION  =       S1%MODULATION.OR.S2%MODULATION
    add_zhe%PARA_IN  =       S1%PARA_IN.OR.S2%PARA_IN

    add_zhe%FULL_WAY=add_zhe%RADIATION.OR.add_zhe%stochastic.OR.add_zhe%ENVELOPE.OR.add_zhe%SPIN.OR.add_zhe%MODULATION
  END FUNCTION add_zhe

  FUNCTION sub_zhe( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) sub_zhe
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2
    logical(lp) dum1,dum2,tt1,tt2

    if(s2%totalpath/=0.and.s2%totalpath/=1) then 
      sub_zhe=s1
      return
    endif
    if(s1%totalpath/=0.and.s1%totalpath/=1) then 
      sub_zhe=s1
      return
    endif

    tt1=s1%only_2d
    tt2=s2%only_2d

    sub_zhe%TOTALPATH=0
    dum1=S1%TOTALPATH==1
    dum2=S2%TOTALPATH==1
    if(dum1.min.dum2) sub_zhe%TOTALPATH=1

    !    sub_zhe%TOTALPATH=  S1%TOTALPATH.min.S2%TOTALPATH

    !   sub_zhe%EXACTMIS    =       (S1%EXACTMIS.min.S2%EXACTMIS).and.(.not.sixtrack_compatible)
    sub_zhe%RADIATION  =  S1%RADIATION.min.S2%RADIATION
    sub_zhe%NOCAVITY =  S1%NOCAVITY.min.S2%NOCAVITY
    sub_zhe%TIME     =  S1%TIME.min.S2%TIME
    sub_zhe%FRINGE   =       S1%FRINGE.min.S2%FRINGE
    sub_zhe%stochastic   =       S1%stochastic.min.S2%stochastic
    sub_zhe%ENVELOPE   =       S1%ENVELOPE.min.S2%ENVELOPE
    sub_zhe%ONLY_4D  =       S1%ONLY_4D.min.S2%ONLY_4D
    sub_zhe%ONLY_2D  =       S1%ONLY_2D.min.S2%ONLY_2D
    sub_zhe%DELTA  =       S1%DELTA.min.S2%DELTA
    sub_zhe%SPIN  =       S1%SPIN.min.S2%SPIN
    sub_zhe%MODULATION  = S1%MODULATION.min.S2%MODULATION
    sub_zhe%PARA_IN  =       (S1%PARA_IN.MIN.S2%PARA_IN)
    !    sub_zhe%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(sub_zhe%stochastic) THEN
       sub_zhe%RADIATION=T
    ENDIF
 !   IF(sub_zhe%ENVELOPE) THEN
 !      sub_zhe%RADIATION=T
 !   ENDIF
    IF(sub_zhe%DELTA) THEN
        if(tt1.or.tt2) then
          sub_zhe%ONLY_2D=T
        else
          sub_zhe%ONLY_4D=T
        endif
       sub_zhe%NOCAVITY =  T
    ENDIF

    IF(sub_zhe%ONLY_4D) THEN
       sub_zhe%TOTALPATH=  0
       sub_zhe%RADIATION  =  F
       sub_zhe%ENVELOPE  =  F
       sub_zhe%stochastic  =  F
       sub_zhe%NOCAVITY =  T
       sub_zhe%stochastic   =  F
    ENDIF

    IF(sub_zhe%ONLY_2D) THEN
       sub_zhe%TOTALPATH=  0
       sub_zhe%RADIATION  =  F
       sub_zhe%ENVELOPE  =  F
       sub_zhe%stochastic  =  F
       sub_zhe%NOCAVITY =  T
       sub_zhe%stochastic   =  F
    ENDIF
    sub_zhe%FULL_WAY=sub_zhe%RADIATION.OR.sub_zhe%stochastic.OR.sub_zhe%ENVELOPE.OR.sub_zhe%SPIN.OR.sub_zhe%MODULATION
  END FUNCTION sub_zhe

  FUNCTION PARA_REMA_zhe(S1)   ! UNARY +
    implicit none
    TYPE (INTERNAL_STATE) PARA_REMA_zhe
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1

    PARA_REMA_zhe         =    S1
    PARA_REMA_zhe%PARA_IN =     T

  END FUNCTION PARA_REMA_zhe

  

!!!!!!!!!!!!!!!!!!!!   tree tracking for PTC using stuff in 
  SUBROUTINE track_TREE_G_complex(T,XI)
 !   use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    REAL(DP), INTENT(INOUT) :: XI(:)
    REAL(DP) XT(lno),XF(lnv),XM(lno+1),XX
    INTEGER JC,I,IV

    XT=0.0_dp
    XF=0.0_dp
    XM=0.0_dp

    do i=1,T%np
       xt(i)=xi(i)
    enddo
    do i=1,T%np
       xf(i) = T%cc(i)
    enddo

    XM(1) = 1.0_dp
    JC=T%np

    do i=1,(T%N-T%np)/T%np
       !
       xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%np
          jc=jc+1
          xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo
    do i=1,size(xi)
       xI(i)=xF(i)
    enddo

  END SUBROUTINE track_TREE_G_complex



subroutine print_tree_element(t,mf)
implicit none
type(tree_element) t
 
integer i,mf
!   write(mf,'(a204)') t%file
write(mf,'(3(1X,i8))') t%N,t%NP,t%no ! ,t%ng
do i=1,t%n
 write(mf,'(1X,G20.13,1x,i8,1x,i8)')  t%cc(i),t%jl(i),t%jv(i)
enddo
write(mf,'(2(1X,L1))') t%symptrack,t%usenonsymp,t%factored
write(mf,'(18(1X,G20.13))') t%fix0,t%fix,t%fixr
do i=1,6
 write(mf,'(6(1X,G20.13))') t%e_ij(i,1:6)
enddo
do i=1,6
 write(mf,'(6(1X,G20.13))') t%rad(i,1:6)
enddo
 write(mf,'(3(1X,G20.13))') t%ds,t%beta0,t%eps

end subroutine print_tree_element

subroutine print_tree_elements_zhe(t,mf)
implicit none
type(tree_element) t(:)
 
integer i,mf

 do i=1,size(t)
  call print_tree_element(t(i),mf)
 enddo

end subroutine print_tree_elements_zhe
 
subroutine read_tree_element(t,mf)
implicit none
type(tree_element) t
 
integer i,mf
 
 ! read(mf,'(a204)') t%file
!read(mf,*) t%N,t%NP,t%no
do i=1,t%n
 read(mf,*)  t%cc(i),t%jl(i),t%jv(i) 
enddo
read(mf,*) t%symptrack,t%usenonsymp,t%factored
read(mf,'(18(1X,G20.13))') t%fix0,t%fix,t%fixr
do i=1,6
 read(mf,*) t%e_ij(i,1:6)
enddo
do i=1,6
 read(mf,*) t%rad(i,1:6)
enddo
 read(mf,*) t%ds,t%beta0,t%eps

end subroutine read_tree_element

subroutine read_tree_elements(t,mf)
implicit none
type(tree_element) t(:)
 
integer i,mf

 do i=1,size(t)
  call read_tree_element(t(i),mf)
 enddo

end subroutine read_tree_elements

  SUBROUTINE track_TREE_probe_complexr(T,xs,dofix0,dofix,sta,jump)
!    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T(:)
    logical, optional :: jump
    type(probe) xs
    real(dp) x(size_tree),x0(size_tree),s0(3,3),r(3,3),dx6,beta,q(3),p(3),qg(3),qf(3)
    real(dp) normb,norm 
    integer i,j,k,ier,is
    type(internal_state) sta
    logical dofix0,dofix,jumpnot
    type(quaternion)qu

   check_stable_zhe=.true.
       xs%u=.false.

    jumpnot=.true.
    if(present(jump)) jumpnot=.not.jump
    
 

   

    x=0.e0_dp
    x0=0.e0_dp
    do i=1,6
      x(i)=xs%x(i)
      x0(i)=xs%x(i)
    enddo
!      x0(1:6)=x(1:6)
      x(7:12)=x(1:6)
 if(jumpnot) then
     if(.not.sta%time) then
     dx6=x(6)
     x(5)=(2*x(5)+x(5)**2)/(sqrt(1.0_dp/t(1)%beta0**2+2.0_dp*x(5)+x(5)**2)+1.0_dp/t(1)%beta0)
     x(11)=x(5)
     x0(5)=x(5)
    endif

    if(dofix0) then
     do i=1,6
      x(i)=x(i)-t(1)%fix0(i)
      x0(i)=x0(i)-t(1)%fix0(i)
     enddo
      x(7:12)=x(1:6)
    endif

    if(sta%radiation) then
      x(1:6)=matmul(t(1)%rad,x(1:6))
      x0(1:6)=x(1:6)
!      x(1:6)=x(1:6)
      x(7:12)=x(1:6)
    endif
    

!!!

endif ! jumpnot

!!! symplectic here
if(t(3)%symptrack) then
    do i=1,3
     q(i)=x(2*i-1)
     p(i)=x(2*i)
    enddo
endif

 if(t(3)%usenonsymp.or..not.t(3)%symptrack) then
    call track_TREE_G_complex(T(1),X(1:6))
 else
    do i=1,3
     x(2*i-1)=0.d0   ! use non symplectic as approximation
    enddo
  endif
!!! symplectic here!! symplectic here
if(t(3)%symptrack) then
   do i=1,3
     qf(i)=x(2*i-1)   ! use non symplectic as approximation
    enddo
normb=1.d38
do is=1,nrmax
   do i=1,3
     x0(2*i)=p(i)
     x0(2*i-1)=qf(i)  
     qg(i)=0
    enddo
    call track_TREE_G_complex(T(3),X0(1:15))
 
    do i=1,3
    do j=1,3
     r(i,j)=x0(ind_spin(i,j))
    enddo
    enddo
    call matinv(r,r,3,3,ier)
    if(ier/=0) then
     write(6,*) "matinv failed in track_TREE_probe_complexr in zhe"
     xs%u=.true.
     check_stable_zhe=.false.
     stop
    endif
    do i=1,3
    do j=1,3
      qg(i)=r(i,j)*(q(j)-x0(2*j-1)) + qg(i)
    enddo
    enddo
    do i=1,3

     qf(i) = qf(i) + qg(i)
    enddo
   norm=abs(qg(1))+abs(qg(2))+abs(qg(3))


   if(norm>t(3)%eps) then
     normb=norm
   else
     if(normb<=norm) then 
       x(1)=qf(1)
       x(3)=qf(2)
       x(5)=qf(3)
       x(2)=x0(2)
       x(4)=x0(4)
       x(6)=x0(6)       

       x(1:6)=matmul(t(3)%rad,x(1:6))
       exit
     endif
     normb=norm
   endif


enddo  ! is 
 if(is>nrmax-10) then
   xs%u=.true.
  check_stable_zhe=.false.
 
  return
 endif
!!!    
 endif

if(jumpnot) then
    if(sta%spin) then  ! spin
    call track_TREE_G_complex(T(2),X(7:15))

     if(xs%use_q) then
       do k=0,3
         qu%x(k)=x(7+k)
       enddo 
 
       xs%q=qu*xs%q
       xs%q%x=xs%q%x/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)
     else

    s0=0.0e0_dp
 
    do i=1,3
    do j=1,3
     r(i,j)=x(ind_spin(i,j))
    enddo
    enddo

    call orthonormaliser(r)
    
    do k=1,3
     s0(k,1:3)=0.0e0_dp
     do i=1,3
     do j=1,3
!       s0(k,i)=x(ind_spin(i,j))*xs%s(k)%x(j)+s0(k,i)
        s0(k,i)=r(i,j)*xs%s(k)%x(j)+s0(k,i)
     enddo
    enddo
    enddo

    do k=1,3
     do j=1,3
       xs%s(k)%x(j)=s0(k,j)
     enddo
    enddo   
endif
    endif ! spin


    if(dofix) then
       if(sta%radiation) then
         do i=1,6
           x(i)=x(i)+t(1)%fixr(i)
         enddo
       else
         do i=1,6
           x(i)=x(i)+t(1)%fix(i)
         enddo
       endif
    endif


    if(.not.sta%time) then
     dx6=X(6)-dx6
      beta=sqrt(1.0_dp+2.0_dp*x(5)/t(1)%beta0+x(5)**2)/(1.0_dp/t(1)%BETA0 + x(5))
      x(6)=x(6)-dx6+beta*dx6 +  (beta/t(1)%beta0-1.0_dp)*t(1)%ds
      x(5)=(2.0_dp*x(5)/t(1)%beta0+x(5)**2)/(sqrt(1.0_dp+2.0_dp*x(5)/t(1)%beta0+x(5)**2)+1.0_dp)
       if(sta%totalpath==1) then
        x(6)=x(6)+t(1)%ds
       endif
    else
        if(sta%totalpath==1) then
        x(6)=x(6)+t(1)%ds/t(1)%beta0 
       endif     
    endif
endif ! jumpnot

    do i=1,6
      xs%x(i)=x(i)
    enddo

  end SUBROUTINE track_TREE_probe_complexr

  SUBROUTINE orthonormaliser(r)
   implicit none
   real(dp)  r(3,3),id(3,3),rt(3,3),eps,a,ab
   integer nmax,i,j,k
! Furmanizing the rotation 
    eps=1.d-8
    nmax=1000
    id=0
    do i=1,3
      id(i,i)=1.5e0_dp
    enddo
    ab=1.d8
    do i=1,nmax
     rt=matmul(r,transpose(r))
     r= matmul((id-0.5e0_dp*rt),r)

     a=-3.e0_dp
     do j=1,3
     do k=1,3
      a=a+abs(rt(j,k))
     enddo
     enddo
     a=abs(a)
     if(a<eps) then
      if(a>=ab) exit
      ab=a
     endif
    enddo
    if(i>nrmax-10) then
     write(6,*) i, a, "did not converge in orthonormaliser"
     read(5,*) i
      stop
    endif 
  end SUBROUTINE orthonormaliser



  SUBROUTINE track_TREE_G_complexr(T,XI)
 !   use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    REAL(DP), INTENT(INOUT) :: XI(:)
    REAL(DP) XT(lno),XF(lnv),XM(lno+1),XX
    INTEGER JC,I,IV

    XT=0.0_dp
    XF=0.0_dp
    XM=0.0_dp

    do i=1,T%np
       xt(i)=xi(i)
    enddo
    do i=1,T%np
       xf(i) = T%cc(i)
    enddo

    XM(1) = 1.0_dp
    JC=T%np

    do i=1,(T%N-T%np)/T%np
       !
       xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%np
          jc=jc+1
          xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo
    do i=1,size(xi)
       xI(i)=xF(i)
    enddo

  END SUBROUTINE track_TREE_G_complexr




 
  SUBROUTINE orthonormalise(r)
   implicit none
   real(dp)  r(3,3),id(3,3),rt(3,3),eps,a,ab
   integer nmax,i,j,k
! Furmanizing the rotation 
    eps=1.d-8
    nmax=1000
    id=0
    do i=1,3
      id(i,i)=1.5e0_dp
    enddo
    ab=1.d8
    do i=1,nmax
     rt=matmul(r,transpose(r))
     r= matmul((id-0.5e0_dp*rt),r)

     a=-3.e0_dp
     do j=1,3
     do k=1,3
      a=a+abs(rt(j,k))
     enddo
     enddo
     a=abs(a)
     if(a<eps) then
      if(a>=ab) exit
      ab=a
     endif
    enddo
    if(i>nrmax-10) then
     write(6,*) i, a, "did not converge in orthonormaliser"
      stop
    endif 
  end SUBROUTINE orthonormalise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   new zhe tracking   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE track_TREE_probe_complex_zhe(T,xs,spin,rad,stoch,linear,slim)
!    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT),target, INTENT(INout) :: T(3)
 
    type(probe) xs
    real(dp) x(size_tree),x0(size_tree),s0(3,3),r(3,3),dx6,beta,q(3),p(3),qg(3),qf(3)
    real(dp) normb,norm,x0_begin(size_tree),xr(6),normbb 
    integer i,j,k,ier,is
    logical, optional  :: spin,stoch,rad,linear,slim
    logical  spin0,stoch0,rad0,doit,as_is0,slim0
    integer no1
    type(quaternion) qu
    as_is0=t(1)%usenonsymp
    spin0=.true.
    stoch0=.true.
    rad0=.true.
    no1=t(1)%no
    slim0=.false.
    if(present(slim)) slim0=slim
    if(present(linear)) then
     if(linear) no1=1
    endif
    if(present(spin)) spin0=spin
    if(present(stoch)) stoch0=stoch
    if(present(rad)) rad0=rad

    if(as_is0) rad0=.true.
    doit=rad0.or.stoch0
    x=0
    x0=0
!    nrmax=1000
   check_stable_zhe=.true.
       xs%u=.false.
!!!! put stochastic kick in front per Sagan
 if(stoch0) then 

    do i=1,6
      x(i)=xs%x(i)-t(1)%fix0(i)
    enddo

    xr=0.0_dp
  do i=1,6
    xr(i)=GRNF_zhe()*t(2)%fix0(i)  
  enddo
    xr(1:6)=matmul(t(2)%rad,xr)

    x=x+xr(1:6)

    do i=1,6
      xs%x(i)=x(i)+t(1)%fix0(i)
    enddo
 endif
!!!!!!!!!!!!!!!!!!!
    x=0.e0_dp
    x0=0.e0_dp
x0_begin=0.0_dp
    do i=1,6
      x(i)=xs%x(i)
      x0(i)=xs%x(i)
      x0_begin(i)=xs%x(i)
    enddo
!      x0(1:6)=x(1:6)
   !   x(7:12)=x(1:6)  remove4/9/2018


if(doit) then

     do i=1,6
      x(i)=x(i)-t(1)%fix0(i)
      x0(i)=x0(i)-t(1)%fix0(i)
      x0_begin(i)=x0_begin(i)-t(1)%fix0(i)
     enddo
else
     do i=1,6
      x(i)=x(i)-t(3)%fix0(i)
      x0(i)=x0(i)-t(3)%fix0(i)
      x0_begin(i)=x0_begin(i)-t(3)%fix0(i)
     enddo
endif
      x(7:12)=x(1:6)
      x0_begin(7:12)= x0_begin(1:6)


!  if(rad0)   call track_TREE_G_complex(T(1),X(1:6))


      x0(1:6)=x(1:6)
      x(7:12)=x(1:6)
if(no1>1.and.(.not.as_is0)) then
    do i=1,3
     q(i)=x(2*i-1)
     p(i)=x(2*i)
    enddo

 !else
    do i=1,3
     x(2*i-1)=0.d0   ! use non symplectic as approximation
    enddo
 ! endif
!!! symplectic here!! symplectic here
! if(t(3)%symptrack) then
   do i=1,3
     qf(i)=x(2*i-1)   ! use non symplectic as approximation
    enddo
normb=1.d38
do is=1,nrmax
   do i=1,3
     x0(2*i)=p(i)
     x0(2*i-1)=qf(i)  
     qg(i)=0
    enddo

    call track_TREE_G_complex(T(3),X0(1:15))
 
    do i=1,3
    do j=1,3
     r(i,j)=x0(ind_spin(i,j))
    enddo
    enddo
    call matinv(r,r,3,3,ier)
    if(ier/=0) then
     write(6,*) "matinv failed in track_TREE_probe_complex_zhe"
       check_stable_zhe=.false.
       xs%u=.true.
      return
     stop
    endif
    do i=1,3
    do j=1,3
      qg(i)=r(i,j)*(q(j)-x0(2*j-1)) + qg(i)
    enddo
    enddo
    do i=1,3

     qf(i) = qf(i) + qg(i)
    enddo
   norm=abs(qg(1))+abs(qg(2))+abs(qg(3))
!write(6,*) is,normb,norm
   if(norm>t(3)%eps) then
      normbb=normb  ! saving for debugging
     normb=norm
   else
      normbb=abs(qf(1))+abs(qf(2))+abs(qf(3))
      
     if(normb<=norm) then 
       x(1)=qf(1)
       x(3)=qf(2)
       x(5)=qf(3)
       x(2)=x0(2)
       x(4)=x0(4)
       x(6)=x0(6)       

       x(1:6)=matmul(t(3)%rad,x(1:6))

       exit
     endif
     normb=norm
   endif

       nrmax_used=is

enddo  ! is 
 if(is>nrmax-10) then
   if(c_verbose_zhe) write(6,*) " Too many iterations ",normbb,norm,t(3)%eps
   xs%u=.true.
   check_stable_zhe=.false.
  return
 endif
elseif(.not.as_is0) then
       x(1:6)=matmul(t(3)%rad,x(1:6))
!!!    
 endif  ! no > 1



!if(jumpnot) then
    if(spin0) then  ! spin
    if(slim0) then
     if(xs%use_q) then
        ! actually using quaternion if SLIM not not used
          qf=0
          qf(1)=t(2)%e_ij(4,4)*xs%q%x(1)+t(2)%e_ij(4,6)*xs%q%x(3)
          qf(3)=t(2)%e_ij(6,4)*xs%q%x(1)+t(2)%e_ij(6,6)*xs%q%x(3)
          do i=1,3

            qf(1)=qf(1)+t(2)%e_ij(1,i)*x0_begin(i)
            qf(3)=qf(3)+t(2)%e_ij(3,i)*x0_begin(i)

          enddo
          qf(2)=xs%q%x(2)
          
         xs%q%x(1:3)=qf/sqrt(qf(1)**2+qf(2)**2+qf(3)**2)
      else
       write(6,*) "SLIM not permitted unless quaternion is the prime method "
          stop
      endif
    else
     if(xs%use_q) then
    call track_TREE_G_complex(T(2),x0_begin(7:15))

       do k=0,3
         qu%x(k)=x0_begin(7+k)
       enddo 
 
       xs%q=qu*xs%q
       xs%q%x=xs%q%x/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)
     else
    call track_TREE_G_complex(T(2),x0_begin(7:15))

    s0=0.0e0_dp
 
    do i=1,3
    do j=1,3
     r(i,j)=x0_begin(ind_spin(i,j))
    enddo
    enddo

    call orthonormaliser(r)
    
    do k=1,3
     s0(k,1:3)=0.0e0_dp
     do i=1,3
     do j=1,3
        s0(k,i)=r(i,j)*xs%s(k)%x(j)+s0(k,i)
     enddo
    enddo
    enddo

    do k=1,3
     do j=1,3
       xs%s(k)%x(j)=s0(k,j)
     enddo
    enddo 
     endif  
    endif ! slim
    endif ! spin
   
if(as_is0) then 
 if(no1>1) then
  call track_TREE_G_complex(T(1),X(1:6))
 else
       x(1:6)=matmul(t(3)%rad,x(1:6))
 endif
else
  if(rad0)   call track_TREE_G_complex(T(1),X(1:6))
endif

 norm=0
do i=1,6
 norm=norm+abs(x(i))
enddo

 if(norm>1) then
   if(c_verbose_zhe) write(6,*) " unstable " 
   xs%u=.true.
   check_stable_zhe=.false.
  return
 endif

if(doit) then

 
         do i=1,6
           x(i)=x(i)+t(1)%fix(i)
         enddo
else
 
         do i=1,6
           x(i)=x(i)+t(3)%fix(i)
         enddo
endif




    do i=1,6
      xs%x(i)=x(i)
    enddo

  end SUBROUTINE track_TREE_probe_complex_zhe

  SUBROUTINE track_TREE_probe_complex_ji_symp(T,xs)
!    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT),target, INTENT(INout) :: T 
 
    type(probe) xs
    real(dp) x(size_ji),x0(size_ji),s0(3,3),r(3,3),dx6,beta,q(3),p(3),qg(3),qf(3),efd(3),ef
    real(dp) normb,norm,xr(6),normbb
    integer i,j,k,ier,is
 
  

 

 

!    nrmax=1000
   check_stable_zhe=.true.
       xs%u=.false.

!!!!!!!!!!!!!!!!!!!
    x=0.e0_dp
    x0=0.e0_dp
    do i=1,6
      x(i)=xs%x(i)
      x0(i)=xs%x(i)
    enddo
 

     do i=1,6
      x(i)=x(i)-t%fix0(i)
      x0(i)=x0(i)-t%fix0(i)
     enddo

   !   x(7:12)=x(1:6)



      x0(1:6)=x(1:6)
  !    x(7:12)=x(1:6)
!if(t(1)%no>1) then
    do i=1,3
     q(i)=x(2*i-1)
     p(i)=x(2*i)
    enddo


!!! symplectic here!! symplectic here
! if(t(3)%symptrack) then
   do i=1,3
     qf(i)=x(2*i-1)   ! use non symplectic as approximation
    enddo
normb=1.d38
do is=1,nrmax
   do i=1,3
     x0(2*i)=p(i)
     x0(2*i-1)=qf(i)  
     qg(i)=0
    enddo
 
    call track_TREE_G_complex(T,X0(1:size_ji))
 
     ef=1.0_dp
    do i=1,3 
     ef=ef*exp(-t%fixr(i)*qf(i)**2)
    enddo
      
    do i=1,3
     efd(i)= -t%fixr(i)*2*qf(i)*ef
    enddo

    do i=1,3
    do j=1,3
     r(i,j)=ef*x0(ind_ji(i,j))+ efd(j)*x0(2*i-1)+x0(9+ind_ji(i,j))
    enddo
    enddo

    call matinv(r,r,3,3,ier)
    if(ier/=0) then
     write(6,*) "matinv failed in track_TREE_probe_complex_zhe"
       check_stable_zhe=.false.
       xs%u=.true.
      return
     stop
    endif
    do j=1,3
      x(2*j-1)=ef*x0(2*j-1)+x0(6+2*j-1)
    enddo
    do i=1,3
    do j=1,3
      qg(i)=r(i,j)*(q(j)-x(2*j-1)) + qg(i)
    enddo
    enddo
    do i=1,3
     qf(i) = qf(i) + qg(i)
    enddo
   norm=abs(qg(1))+abs(qg(2))+abs(qg(3))
!write(6,*) is,normb,norm
!write(6,*) norm>t%eps,norm,t%eps
   if(norm>t%eps) then
      normbb=normb  ! saving for debugging
     normb=norm
   else
     if(normb<=norm) then 
       x(1)=qf(1)
       x(3)=qf(2)
       x(5)=qf(3)
       do i=1,3
        x(2*i)=efd(i)*x0(size_ji) +  ef*x0(2*i) +x0(2*i+6)
       enddo
     

     !  x(1:6)=matmul(t(3)%rad,x(1:6))

       exit
     endif
     normb=norm
   endif

       nrmax_used=is

enddo  ! is 
 if(is>nrmax-10) then
   if(c_verbose_zhe) write(6,*) " Too many iterations ",normbb,norm,t%eps
   xs%u=.true.
   check_stable_zhe=.false.
  return
 endif
!!!    
 







         do i=1,6
           x(i)=x(i)+t%fix(i)
         enddo



    do i=1,6
      xs%x(i)=x(i)
    enddo

  end SUBROUTINE track_TREE_probe_complex_ji_symp

  SUBROUTINE change_ntot(n)
!    use da_arrays
    IMPLICIT NONE
    integer n
    ntot=n
 end   SUBROUTINE change_ntot

  SUBROUTINE track_TREE_probe_complex_ji_vec(T,xs)
!    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT),target, INTENT(INout) :: T(2) 
 
    type(probe) xs
    real(dp) x(6),h
    integer i
   
  
!    nrmax=1000
   check_stable_zhe=.true.
       xs%u=.false.

!!!!!!!!!!!!!!!!!!!
    x=0.e0_dp
    do i=1,6
      x(i)=xs%x(i)
    enddo
 

     do i=1,6
      x(i)=x(i)-t(1)%fix0(i)
     enddo

  !!!call track_TREE_G_complex(T(1),X0(1:size_ji_vec))

  h=1.0_dp/ntot
   
   do i=1,ntot
        call rk6_vec(x,t(1),h)
   enddo

   do i=1,ntot
        call rk6_vec(x,t(2),h)
   enddo

         do i=1,6
           x(i)=x(i)+t(1)%fix(i)
         enddo



    do i=1,6
      xs%x(i)=x(i)
    enddo

  end SUBROUTINE track_TREE_probe_complex_ji_vec

subroutine fev(x,f,t)
implicit none
TYPE(TREE_ELEMENT), INTENT(INout) :: t 
real(dp), INTENT(INout) ::  x(6)
real(dp)  f(6),e,de(6),x0(size_ji_vec)
integer i

!111 x=0
!write(6,*) "x(1:4)"
!read(5,*) x(1:4)
x0=0
x0(1:6)=x
 
  call track_TREE_G_complex(T,X0(1:size_ji_vec))
 
e=1
do i=1,6
 e=e*exp(-t%fixr(i)*x0(i)**2)
 de(i)=-2*t%fixr(i)*x0(i)
enddo
de=de*e

 
do i=1,3
f(2*i-1)=e*x0(2*i-1)+x0(2*i-1+6)+de(2*i-1)*x0(size_ji_vec)
f(2*i)= e*x0(2*i)+x0(2*i+6)-de(2*i)*x0(size_ji_vec)
enddo

!write(6,*) f(1:4)
!write(6,*) f(5:6)
!goto 111
end subroutine fev

  subroutine rk6_vec(y,gr,h)
    IMPLICIT none
 

    integer ne
    parameter (ne=6)
    real(dp), INTENT(INOUT)::  y(ne)
    real(dp)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne),e(ne),g(ne),o(ne),p(ne)
    TYPE(TREE_ELEMENT), INTENT(INout) :: gr 
    integer j
    real(dp), intent(inout) :: h

 
    call fev(y,f,gr)
 
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/9.0_dp
    enddo

    call fev(yt,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + (a(j) + 3.0_dp*b(j))/24.0_dp
    enddo

    call fev(yt,f,gr)
    do  j=1,ne
       c(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j)+(a(j)-3.0_dp*b(j)+4.0_dp*c(j))/6.0_dp
    enddo

    call fev(yt,f,gr)
    do  j=1,ne
       d(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j) + (-5.0_dp*a(j) + 27.0_dp*b(j) - 24.0_dp*c(j) + 6.0_dp*d(j))/8.0_dp
    enddo

    call fev(yt,f,gr)
    do  j=1,ne
       e(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j) + (221.0_dp*a(j) - 981.0_dp*b(j) + 867.0_dp*c(j)- 102.0_dp*d(j) + e(j))/9.0_dp
    enddo

    call fev(yt,f,gr)
    do   j=1,ne
       g(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j) = y(j)+(-183.0_dp*a(j)+678.0_dp*b(j)-472.0_dp*c(j)-66.0_dp*d(j)+80.0_dp*e(j) + 3.0_dp*g(j))/48.0_dp
    enddo

    call fev(yt,f,gr)
    do  j=1,ne
       o(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j) = y(j)+(716.0_dp*a(j)-2079.0_dp*b(j)+1002.0_dp*c(j)+834.0_dp*d(j)-454.0_dp*e(j)-9.0_dp*g(j)+72.0_dp*o(j))/82.0_dp
    enddo


    call fev(yt,f,gr)
    do  j=1,ne
       p(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+(41.0_dp*a(j)+216.0_dp*c(j)+27.0_dp*d(j)+272.0_dp*e(j)+27.0_dp*g(j)+216.0_dp*o(j)+41.0_dp*p(j))/840.0_dp
    enddo


    return
  end  subroutine rk6_vec



  SUBROUTINE track_TREE_probe_complex_ji(T,xs)
!    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT),target, INTENT(INout) :: T 
 
    type(probe) xs
    real(dp) x(12),ef
    integer i,j,k,ier,is,ns
    type(quaternion) qu

    ns=12



!    nrmax=1000
   check_stable_zhe=.true.
       xs%u=.false.


    x=0.0_dp
    do i=1,6
      x(i)=xs%x(i)
    enddo

     do i=1,6
      x(i)=x(i)-t%fix0(i)
     enddo
ef=1
do i=1,6
ef=ef*exp(t%fixr(i)*x(i)**2)
enddo  
  do i=1,6
  x(i)=x(i)*ef
  enddo  
 
 
   call track_TREE_G_complex(T,X(1:ns)) 


x(1:6)= x(1:6)+x(7:12)

         do i=1,6
           x(i)=x(i)+t%fix(i)
         enddo
 

    


    do i=1,6
      xs%x(i)=x(i)
    enddo

  end SUBROUTINE track_TREE_probe_complex_ji


  real(dp) FUNCTION RANF_zhe()
    implicit none
    integer ia,ic,iq,ir,ih,il,it
    DATA IA/16807/,IC/2147483647/,IQ/127773/,IR/2836/
    IH = zhe_ISEED/IQ
    IL = MOD(zhe_ISEED,IQ)
    IT = IA*IL-IR*IH
    IF(IT.GT.0) THEN
       zhe_ISEED = IT
    ELSE
       zhe_ISEED = IC+IT
    END IF
    RANF_zhe = zhe_ISEED/FLOAT(IC)

       !         t=sqrt(12.d0)*(RANF()-half)
       if(RANF_zhe>0.5_dp) then
          RANF_zhe=1.0_dp
       else
          RANF_zhe=-1.0_dp
       endif
    RETURN
  END FUNCTION RANF_zhe

  real(dp) FUNCTION RANF()
    implicit none
    integer ia,ic,iq,ir,ih,il,it
    DATA IA/16807/,IC/2147483647/,IQ/127773/,IR/2836/
    IH = zhe_ISEED/IQ
    IL = MOD(zhe_ISEED,IQ)
    IT = IA*IL-IR*IH
    IF(IT.GT.0) THEN
       zhe_ISEED = IT
    ELSE
       zhe_ISEED = IC+IT
    END IF
    RANF = zhe_ISEED/FLOAT(IC)
    RETURN
  END FUNCTION RANF

  SUBROUTINE gaussian_seed_zhe(seed,cut)
    implicit none
    integer seed
    real(dp) cut
    zhe_ISEED=seed
    cut_zhe=cut
  end SUBROUTINE gaussian_seed_zhe

  SUBROUTINE nrmax_zhe(nrmax_in)
    implicit none
    integer nrmax_in
    nrmax=nrmax_in
  end SUBROUTINE nrmax_zhe

  SUBROUTINE nrmax_used_zhe(nrmax_u)
    implicit none
    integer nrmax_u
    nrmax_u=nrmax_used
  end SUBROUTINE nrmax_used_zhe

   real(dp) function GRNF_zhe_gaussian()
    implicit none
    real(dp) r1,r2,x 


1   R1 = -LOG(1.0_dp-RANF())
    R2 = 2.0_dp*PI*RANF()
    R1 = SQRT(2.0_dp*R1)
    X  = R1*COS(R2)
    if(abs(x)>cut_zhe) goto 1
     GRNF_zhe_gaussian=x
     
 
    ! Y  = R1*SIN(R2)
    RETURN
  END function GRNF_zhe_gaussian

    subroutine get_seed(k)
    implicit none
    integer k
    k=alex_ISEED
   end  subroutine get_seed

    subroutine set_seed(k)
    implicit none
    integer k
    alex_ISEED=k
   end  subroutine set_seed
 

  
   real(dp) function GRNF_zhe()
    implicit none
    integer ia,ic,iq,ir,ih,il,it
    DATA IA/16807/,IC/2147483647/,IQ/127773/,IR/2836/
    IH = alex_ISEED/IQ
    IL = MOD(alex_ISEED,IQ)
    IT = IA*IL-IR*IH
    IF(IT.GT.0) THEN
       alex_ISEED = IT
    ELSE
       alex_ISEED = IC+IT
    END IF
    GRNF_zhe = alex_ISEED/FLOAT(IC)
    


       if(GRNF_zhe>0.5_dp) then
          GRNF_zhe=1.0_dp
       else
          GRNF_zhe=-1.0_dp
       endif
 
    ! Y  = R1*SIN(R2)
    RETURN
  END function GRNF_zhe
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of   new zhe tracking   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine matinv(a,ai,n,nmx,ier)

    implicit none

    integer i,ier,j,n,nmx
    integer,dimension(nmax)::indx
    real(dp) d
    real(dp),dimension(nmx,nmx)::a,ai
    real(dp),dimension(nmax,nmax)::aw
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif

    aw(1:n,1:n) = a(1:n,1:n)

    call ludcmp_nr(aw,n,nmax,indx,d,ier)
    if (ier .eq. 132) return

    ai(1:n,1:n) = 0.0_dp
    !    forall (i = 1:n) ai(i,i) = one
    do i=1,n
       ai(i,i) = 1.0_dp
    enddo

    do j=1,n
       call lubksb_nr(aw,n,nmax,indx,ai(1,j),nmx)
    enddo

  end subroutine matinv

  !
  subroutine ludcmp_nr(a,n,np,indx,d,ier)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
    !     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ier,imax,j,k,n,np
    integer,dimension(np)::indx
    real(dp) aamax,d,dum,sum
    real(dp),dimension(np,np)::a
    real(dp),dimension(nmax)::vv
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ier=0
    d=1.0_dp
    do i=1,n
       aamax=0.0_dp
       do j=1,n
          if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       enddo
       if(aamax.eq.0.0_dp) then
          ier=132
          return
       endif
       vv(i)=1.0_dp/aamax
    enddo
    do j=1,n
       if(j.gt.1) then
          do i=1,j-1
             sum=a(i,j)
             if(i.gt.1) then
                do k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                enddo
                a(i,j)=sum
             endif
          enddo
       endif
       aamax=0.0_dp
       do i=j,n
          sum=a(i,j)
          if (j.gt.1) then
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if(dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j.ne.imax) then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if(j.ne.n) then
          if(a(j,j).eq.0.0_dp) a(j,j)=tiny
          dum=1.0_dp/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    if(a(n,n).eq.0.0_dp) a(n,n)=tiny
    return
  end subroutine ludcmp_nr
  !
  subroutine lubksb_nr(a,n,np,indx,b,nmx)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
    !     INPUT A: NXN MATRIX IN lu FORM GIVEN BY ludcmp_nr
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,j,ll,n,nmx,np
    integer,dimension(np)::indx
    real(dp) sum
    real(dp),dimension(np,np)::a
    real(dp),dimension(nmx)::b
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ii = 0
    do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if(ii.ne.0) then
          do j=ii,i-1
             sum = sum-a(i,j)*b(j)
          enddo
       else if (sum.ne.0.0_dp) then
          ii = i
       endif
       b(i)=sum
    enddo
    do i=n,1,-1
       sum=b(i)
       if(i.lt.n) then
          do j=i+1,n
             sum = sum-a(i,j)*b(j)
          enddo
       endif

       b(i)=sum/a(i,i)

    enddo
    return
  end subroutine lubksb_nr
  !


!!!!!!!!!!!!!!!!!! old   special for lielib: others in sa_extend_poly

  SUBROUTINE KanalNummer(iff,file,old)
    implicit none
    integer, INTENT(OUT) :: iff
    character(*),optional :: file
    logical,optional :: old
    logical :: opened, exists
    integer :: i,ier

    DO i= 9999, 7, -1
       INQUIRE(UNIT= i, EXIST= exists, OPENED= opened)
       IF (exists .AND. (.NOT. opened)) GOTO 20
    ENDDO
    WRITE (UNIT= *, FMT= *) ' cannot find free unit within the range 7-9999..'
    CALL ReportOpenFiles
    STOP
20  CONTINUE
    iff= I

    if(present(file)) then
       if(present(old)) then
        if(old) then
        open(unit=iff,file=file,status='OLD',iostat=ier)
         if(ier/=0) then
          write(6,*) " file ",file," does not exist "
          stop 864
          endif
        else
         open(unit=iff,file=file,status='NEW',iostat=ier)
         if(ier/=0) then
          write(6,*) " file ",file,"  exists already "
          stop 865
         endif
        endif
       else
        open(unit=iff,file=file)
       endif
    endif
  END SUBROUTINE KanalNummer


  SUBROUTINE ReportOpenFiles
    implicit none

    logical :: opened, exists, named
    CHARACTER(LEN= 400) :: name
    integer :: i

    DO i= 9999, 7, -1
       INQUIRE(UNIT= i, EXIST= exists, OPENED= opened)
       IF (exists .AND. opened) THEN
          INQUIRE(UNIT= i, NAMED= named, NAME= name)
          write (*,7010) i, name(:LEN_TRIM(name))
       ENDIF
    ENDDO
7010 FORMAT(' iUnit:',I3,', name: "',A,'"')

  END SUBROUTINE ReportOpenFiles


subroutine zhe_ini(use_q)
implicit none
logical , optional ::use_q 
integer i,j,k
if(Present(use_q) )use_quaternion=use_q
    ind_spin(1,1)=1+6;ind_spin(1,2)=2+6;ind_spin(1,3)=3+6;
    ind_spin(2,1)=4+6;ind_spin(2,2)=5+6;ind_spin(2,3)=6+6;
    ind_spin(3,1)=7+6;ind_spin(3,2)=8+6;ind_spin(3,3)=9+6;    
    k1_spin(1)=1;k2_spin(1)=1;
    k1_spin(2)=1;k2_spin(2)=2;
    k1_spin(3)=1;k2_spin(3)=3;
    k1_spin(4)=2;k2_spin(4)=1;
    k1_spin(5)=2;k2_spin(5)=2;
    k1_spin(6)=2;k2_spin(6)=3;
    k1_spin(7)=3;k2_spin(7)=1;
    k1_spin(8)=3;k2_spin(8)=2;
    k1_spin(9)=3;k2_spin(9)=3;
     k=12
    do i=1,3
    do j=1,3
     k=k+1
     ind_ji(i,j)=k
    enddo
   enddo
end subroutine zhe_ini


subroutine read_tree_zhe(t,filename)
implicit none
TYPE(TREE_ELEMENT), INTENT(INOUT) :: T(:)
character(*) filename
integer inf,i,n,np,no
 call kanalnummer(inf,filename)
 
  do i=1,SIZE(T)
    read(inf,*) n,np,no
    CALL ALLOC_TREE(t(i),N,NP)
   ! t(i)%Ng=ng
    t(i)%N=n
    t(i)%NP=np
    t(i)%no=no
    call read_tree_element(t(i),inf)
  enddo
close(inf)
 
end subroutine read_tree_zhe

subroutine kill_tree_zhe(t)
implicit none
TYPE(TREE_ELEMENT), INTENT(INOUT) :: T(3)
integer i 

 
  do i=1,3
    CALL kill_tree(t(i))
  enddo

 
end subroutine kill_tree_zhe

!!! quaternion



  FUNCTION unaryADDq( S1 )
    implicit none
    TYPE (quaternion) unaryADDq
    TYPE (quaternion), INTENT (IN) :: S1
 
    unaryADDq=s1

  END FUNCTION unaryADDq

  FUNCTION unarySUBq( S1 )
    implicit none
    TYPE (quaternion) unarySUBq
    TYPE (quaternion), INTENT (IN) :: S1

         unarySUBq%x= -unarySUBq%x

  END FUNCTION unarySUBq


  FUNCTION invq( S1 )
    implicit none
    TYPE (quaternion) invq
    TYPE (quaternion), INTENT (IN) :: S1
    real(dp) norm
     integer i
     invq=s1
              do i=1,3
                invq%x(i)=-invq%x(i)
              enddo
                norm=abs_square(invq)
              do i=0,3
                invq%x(i)=invq%x(i)/norm
              enddo
      
  END FUNCTION invq


  FUNCTION absq( S1 )
    implicit none
    real(dp) absq
    TYPE (quaternion), INTENT (IN) :: S1
    integer i

 
   
     absq=sqrt(abs_square(s1))
  END FUNCTION absq


  FUNCTION absq2( S1 )
    implicit none
    real(dp) absq2
    TYPE (quaternion), INTENT (IN) :: S1
    integer i
 
           absq2=0
       do i=0,3
         absq2 = s1%x(i)**2+absq2
       enddo
  END FUNCTION absq2


  SUBROUTINE  EQUALq(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    type (quaternion),INTENT(IN)::S1
    integer i
     
    do i=0,3
    s2%x(i)=s1%x(i)
    enddo

  end SUBROUTINE  EQUALq

  SUBROUTINE  EQUALqr(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    integer i
 
    do i=0,3
    s2%x(i)=0
    enddo
    s2%x(0)=s1
  end SUBROUTINE  EQUALqr

  SUBROUTINE  EQUALqi(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
 
    do i=0,3
    s2%x(i)=0
    enddo
    s2%x(0)=s1
  end SUBROUTINE  EQUALqi


  FUNCTION POWq( S1, R2 )
    implicit none
    TYPE (quaternion) POWq,temp
    TYPE (quaternion), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster
      temp=1.0_dp

    R22=IABS(R2)
    DO I=1,R22
       temp=temp*s1
    ENDDO
    IF(R2.LT.0) THEN
       temp=invq(temp)
    ENDIF
     powq=temp
 
  END FUNCTION POWq


  FUNCTION addq( S1, S2 )
    implicit none
    TYPE (quaternion) addq
    TYPE (quaternion), INTENT (IN) :: S1, S2

 
       addq%x=s1%x+s2%x
  END FUNCTION addq


  FUNCTION subq( S1, S2 )
    implicit none
    TYPE (quaternion) subq
    TYPE (quaternion), INTENT (IN) :: S1, S2

 
          subq%x=s1%x+s2%x

  END FUNCTION subq

  FUNCTION mulq( S1, S2 )
    implicit none
    TYPE (quaternion) mulq
    TYPE (quaternion), INTENT (IN) :: S1, S2
    integer i
 
 
          mulq=0.0_dp

          mulq%x(0)=s1%x(0)*s2%x(0)-s1%x(1)*s2%x(1)-s1%x(2)*s2%x(2)-s1%x(3)*s2%x(3)

         mulq%x(1)=  s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
         mulq%x(2)=  s1%x(3)*s2%x(1)-s1%x(1)*s2%x(3)
         mulq%x(3)=  s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)

        do i=1,3
         mulq%x(i)= mulq%x(i) + s1%x(0)*s2%x(i)+ s1%x(i)*s2%x(0)
        enddo

  END FUNCTION mulq

  FUNCTION divq( S1, S2 )
    implicit none
    TYPE (quaternion) divq
    TYPE (quaternion), INTENT (IN) :: S1, S2

         
       divq=s1*invq(s2)

  END FUNCTION divq


  SUBROUTINE  printq(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (quaternion),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    INTEGER I,mfi
     mfi=6
     if(present(mfile)) mfi=mfile
      write(mfi,*) " real quaternion "
    DO I=0,3
      write(mfi,*) s1%x(i)
    ENDDO
  END SUBROUTINE printq

end module duan_zhe_map

