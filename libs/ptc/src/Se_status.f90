!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module S_status
  use s_frame
  USE S_extend_poly
!  use anbn
  ! use my_own_1D_TPSA
  !  USE S_pol_user1
  !  USE S_pol_user2
  Use S_pol_sagan
  implicit none
  public

  integer,private,parameter::ST=30
  !
  integer, parameter :: KIND0 = ST
  integer, parameter :: KIND1 = ST+1
  integer, parameter :: KIND2 = ST+2
  integer, parameter :: KIND3 = ST+3
  integer, parameter :: KIND4 = ST+4
  integer, parameter :: KIND5 = ST+5
  integer, parameter :: KIND6 = ST+6
  integer, parameter :: KIND7 = ST+7
  integer, parameter :: KIND8 = ST+8
  integer, parameter :: KIND9 = ST+9
  integer, parameter :: KIND10 =ST+10
  integer, parameter :: KIND11 =ST+11
  integer, parameter :: KIND12 =ST+12
  integer, parameter :: KIND13 =ST+13
  integer, parameter :: KIND14 =ST+14
  integer, parameter :: KIND15 =ST+15
  integer, parameter :: KIND16 =ST+16
  integer, parameter :: KIND17 =ST+17
  integer, parameter :: KIND18 =ST+18
  integer, parameter :: KIND19 =ST+19
  integer, parameter :: KIND20 =ST+20     !  MADLIKE wedges on RBEND
  integer, parameter :: KIND21 =ST+21     !  travelling wave cavity
  integer, parameter :: KIND22 =ST+22     !
  integer, parameter :: KIND23 =ST+23     !
  !  integer, parameter :: KINDFITTED = KIND23+1
  !  integer, parameter :: KINDUSER1 = KIND23+2
  !  integer, parameter :: KINDUSER2 = KIND23+3
  integer, parameter :: KINDwiggler = KIND23+2
  !  integer, parameter :: KINDmu      = KIND23+3
  integer, parameter :: KINDpa     = KIND23+3
  integer, parameter :: drift_kick_drift = kind2
  integer, parameter :: matrix_kick_matrix = kind7
  integer, parameter :: kick_sixtrack_kick = kind6
  character(6) ind_stoc(ndim2)
  !  Making PTC leaner is false
  !

  integer, target :: MADKIND2=KIND2
  integer, target :: MADKIND3N=KIND3
  integer, target :: MADKIND3S=KIND3
  integer, pointer :: MADTHICK
  integer, pointer :: MADTHIN_NORMAL
  integer, pointer :: MADTHIN_SKEW
  LOGICAL(lp), target :: MADLENGTH=.false.
  LOGICAL(lp), target :: MAD=.false.
  LOGICAL(lp), target :: EXACT_MODEL = .false.
  INTEGER, target:: NSTD,METD
  ! TYPE(B_CYL) SECTOR_B
  !TYPE(B_CYL),ALLOCATABLE ::  S_B(:)
  !  INTEGER, TARGET :: NDPT_OTHER = 0
  real(dp) CRAD,CFLUC
  !  real(dp) YOSK(0:4), YOSD(4)    ! FIRST 6TH ORDER OF YOSHIDA
  !  real(dp),PARAMETER::AAA=-0.25992104989487316476721060727823e0_dp  ! fourth order integrator
  !  real(dp),PARAMETER::FD1=half/(one+AAA),FD2=AAA*FD1,FK1=one/(one+AAA),FK2=(AAA-one)*FK1
  INTEGER , target :: CAVITY_TOTALPATH=1   !  default is fake

  LOGICAL(lp) :: firsttime_coef=.true.

  PRIVATE EQUALt,ADD,PARA_REMA,EQUALtilt
  !PRIVATE DTILTR,DTILTP,DTILTS
  PRIVATE DTILTR_EXTERNAL,DTILTP_EXTERNAL
  PRIVATE CHECK_APERTURE_R,CHECK_APERTURE_P !,CHECK_APERTURE_S
  LOGICAL(lp), target:: electron
  real(dp), target :: muon=1.0_dp
  LOGICAL(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.
  ! include "a_def_all_kind.inc"    ! sept 2007
  ! include "a_def_sagan.inc"
  ! include "a_def_element_fibre_layout.inc"
  !  !  include "a_def_user1.inc"
  !  !!  include "a_def_arbitrary.inc"
  !  !  include "a_def_user2.inc"

  TYPE(INTERNAL_STATE),PARAMETER::DEFAULT0=INTERNAL_STATE   (0,f,f,f,f,f,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::TOTALPATH0=INTERNAL_STATE (1,f,f,f,f,f,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::TIME0=INTERNAL_STATE      (0,t,f,f,f,f,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::RADIATION0=INTERNAL_STATE (0,f,t,f,f,f,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::NOCAVITY0=INTERNAL_STATE  (0,f,f,t,f,f,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::FRINGE0=INTERNAL_STATE    (0,f,f,f,t,f,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::STOCHASTIC0=INTERNAL_STATE(0,f,f,f,f,t,f,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::ENVELOPE0=INTERNAL_STATE  (0,f,f,f,f,f,t,f,f,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::ONLY_4d0=INTERNAL_STATE   (0,f,f,t,f,f,f,f,t,f,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::DELTA0=INTERNAL_STATE     (0,f,f,t,f,f,f,f,t,t,f,f)
  TYPE(INTERNAL_STATE),PARAMETER::SPIN0=INTERNAL_STATE      (0,f,f,f,f,f,f,f,f,f,t,f)
  TYPE(INTERNAL_STATE),PARAMETER::MODULATION0=INTERNAL_STATE(0,f,f,f,f,f,f,f,f,f,f,t)

  TYPE(INTERNAL_STATE), target ::  DEFAULT=DEFAULT0
  TYPE(INTERNAL_STATE), target ::  TOTALPATH=TOTALPATH0
  TYPE(INTERNAL_STATE), target ::  RADIATION=RADIATION0
  TYPE(INTERNAL_STATE), target ::  NOCAVITY=NOCAVITY0
  TYPE(INTERNAL_STATE), target ::  FRINGE=FRINGE0
  TYPE(INTERNAL_STATE), target ::  TIME=TIME0
  TYPE(INTERNAL_STATE), target ::  STOCHASTIC=STOCHASTIC0
  TYPE(INTERNAL_STATE), target ::  ENVELOPE=ENVELOPE0
  TYPE(INTERNAL_STATE), target ::  ONLY_4D=ONLY_4D0
  TYPE(INTERNAL_STATE), target ::  DELTA=DELTA0
  TYPE(INTERNAL_STATE), target ::  SPIN=SPIN0
  TYPE(INTERNAL_STATE), target ::  MODULATION=MODULATION0
   type(acceleration), pointer :: acc
   type(acceleration), pointer :: accFIRST
   type(fibre), pointer :: paccfirst
   type(fibre), pointer :: paccthen

  !  private s_init,S_init_berz,MAKE_STATES_0,MAKE_STATES_m,print_s,CONV
  private s_init,MAKE_STATES_0,MAKE_STATES_m,print_s,CONV
  LOGICAL(lp), target :: compute_stoch_kick = .true.
  private alloc_p,equal_p,dealloc_p,alloc_A,equal_A,dealloc_A
  PRIVATE KILL_S_APERTURE,ALLOC_S_APERTURE
  !,NULL_p
  PRIVATE B2PERPR,B2PERPP !,S_init_berz0
  type(tilting) tilt
  private minu
  real(dp) MADFAC(NMAX)
  CHARACTER(24) MYTYPE(-100:100)
  private check_S_APERTURE_r,check_S_APERTURE_out_r
  private check_S_APERTURE_p,check_S_APERTURE_out_p
  type(my_1D_taylor) val_del
  logical(lp) :: ramp=my_false
  logical(lp) :: accelerate=my_false, first_particle=my_false

  TYPE B_CYL
     integer firsttime
     integer, POINTER ::  nmul,n_mono
     integer, DIMENSION(:), POINTER   :: i,j
     real(dp), DIMENSION(:,:), POINTER   :: a_x,a_y,b_x,b_y
  END  TYPE B_CYL

  TYPE(B_CYL),ALLOCATABLE ::  S_B(:)

  INTERFACE OPERATOR (.min.)
     MODULE PROCEDURE minu                       ! to define the minus of Schmidt
  END INTERFACE

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUALt
     MODULE PROCEDURE EQUALtilt
     MODULE PROCEDURE equal_p
     MODULE PROCEDURE equal_A
  end  INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add
     MODULE PROCEDURE PARA_REMA
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE sub
  END INTERFACE

  INTERFACE MAKE_STATES
     MODULE PROCEDURE MAKE_STATES_m
     MODULE PROCEDURE MAKE_STATES_0
  END INTERFACE

  INTERFACE CHECK_APERTURE
     MODULE PROCEDURE CHECK_APERTURE_R
     MODULE PROCEDURE CHECK_APERTURE_P
     !     MODULE PROCEDURE CHECK_APERTURE_S
  END INTERFACE

  INTERFACE check_S_APERTURE
     MODULE PROCEDURE check_S_APERTURE_r
     MODULE PROCEDURE check_S_APERTURE_p
  END INTERFACE

  INTERFACE check_S_APERTURE_out
     MODULE PROCEDURE check_S_APERTURE_out_r
     MODULE PROCEDURE check_S_APERTURE_out_p
  END INTERFACE

  INTERFACE init
     MODULE PROCEDURE s_init
     !     MODULE PROCEDURE S_init_berz0
     !     MODULE PROCEDURE S_init_berz
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE print_s
  END INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE alloc_p
     MODULE PROCEDURE alloc_A
     MODULE PROCEDURE ALLOC_S_APERTURE
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE dealloc_p
     MODULE PROCEDURE dealloc_A
     MODULE PROCEDURE KILL_S_APERTURE
  END INTERFACE

  INTERFACE B2PERP
     MODULE PROCEDURE B2PERPR
     MODULE PROCEDURE B2PERPP
  END INTERFACE

  INTERFACE DTILTD
     MODULE PROCEDURE DTILTR_EXTERNAL
     MODULE PROCEDURE DTILTP_EXTERNAL       ! EXTERNAL
  END INTERFACE



CONTAINS

  real(dp) function cradf(p)
    implicit none
    type (MAGNET_CHART), pointer:: P
    cradf=crad*p%p0c**3
  end function cradf

  real(dp) function cflucf(p)
    implicit none
    type (MAGNET_CHART), pointer:: P
    cflucf=cfluc*p%p0c**5
  end function cflucf


  SUBROUTINE  NULL_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    nullify(P%KIND);nullify(P%R);nullify(P%X);nullify(P%Y);nullify(P%dX);nullify(P%dY);
  end subroutine NULL_A

  SUBROUTINE  alloc_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    nullify(p)
    allocate(p)
    CALL NULL_A(p)
    ALLOCATE(P%R(2));ALLOCATE(P%X);ALLOCATE(P%Y);ALLOCATE(P%KIND);
    P%KIND=0; P%R=0.0_dp;P%X=0.0_dp;P%Y=0.0_dp;
    ALLOCATE(P%DX);ALLOCATE(P%DY);
    P%DX=0.0_dp;P%DY=0.0_dp;
  end subroutine alloc_A

  SUBROUTINE  dealloc_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    if(associated(p%R)) then
       DEALLOCATE(P%R);DEALLOCATE(P%X);DEALLOCATE(P%Y);DEALLOCATE(P%KIND);
       DEALLOCATE(P%DX);DEALLOCATE(P%DY);
    endif
  end SUBROUTINE  dealloc_A



  SUBROUTINE  NULL_p(p)
    implicit none
    type (MAGNET_CHART), pointer:: P

    nullify(P%LD);nullify(P%B0);nullify(P%LC);
    nullify(P%TILTD);  nullify(P%dir);
    !    nullify(P%beta0);nullify(P%gamma0I);nullify(P%gambet);nullify(P%charge);
    nullify(P%P0C);
    nullify(P%EDGE)
    !    nullify(P%TOTALPATH)
    nullify(P%EXACT);  !nullify(P%RADIATION);nullify(P%NOCAVITY);
    nullify(P%permFRINGE);
    nullify(P%KILL_ENT_FRINGE);nullify(P%KILL_EXI_FRINGE);nullify(P%bend_fringe);  !nullify(P%TIME);
    nullify(P%METHOD);nullify(P%NST);
    nullify(P%NMUL);  !nullify(P%spin);
    nullify(P%F);
    nullify(P%APERTURE);
    nullify(P%A);

  end subroutine NULL_p




  SUBROUTINE  alloc_p(p)
    implicit none
    type (MAGNET_CHART), pointer:: P

    nullify(P);
    ALLOCATE(P);
    CALL NULL_P(P)

    ALLOCATE(P%LD);ALLOCATE(P%B0);ALLOCATE(P%LC);
    P%LD=0.0_dp;P%B0=0.0_dp;P%LC=0.0_dp;
    ALLOCATE(P%TILTD);P%TILTD=0.0_dp;
    ! ALLOCATE(P%beta0);
    ! ALLOCATE(P%MASS0);
    ! ALLOCATE(P%gamma0I);
    ! ALLOCATE(P%gambet);
    ALLOCATE(P%P0C);
    ! P%beta0 =one;
    ! P%MASS0 =one;
    !P%gamma0I=zero;P%gambet =zero;
    P%P0C =0.0_dp;
    ALLOCATE(P%EDGE(2));P%EDGE(1)=0.0_dp;P%EDGE(2)=0.0_dp;
    !    ALLOCATE(P%TOTALPATH); ! PART OF A STATE INITIALIZED BY EL=DEFAULT
    ALLOCATE(P%EXACT);  !ALLOCATE(P%RADIATION);ALLOCATE(P%NOCAVITY);
        ALLOCATE(P%permFRINGE);
    ALLOCATE(P%KILL_ENT_FRINGE);ALLOCATE(P%KILL_EXI_FRINGE);ALLOCATE(P%bend_fringe); !ALLOCATE(P%TIME);
    ALLOCATE(P%METHOD);ALLOCATE(P%NST);P%METHOD=2;P%NST=1;
    ALLOCATE(P%NMUL);P%NMUL=0;
    !    ALLOCATE(P%spin);
    !   ALLOCATE(P%TRACK);P%TRACK=.TRUE.;
    P%permFRINGE=.FALSE.
    P%KILL_ENT_FRINGE=.FALSE.
    P%KILL_EXI_FRINGE=.FALSE.
    P%bend_fringe=.false.
    call alloc(p%f)
    ! if(junk) ccc=ccc+1

  end subroutine alloc_p


  SUBROUTINE  dealloc_p(p)
    implicit none
    type (MAGNET_CHART), pointer:: P
    if(.not.associated(p)) return

    !    if(associated(P%dir)) then
    !    endif
    if(associated(p%LD)) DEALLOCATE(P%LD);
    if(associated(p%B0)) DEALLOCATE(P%B0);
    if(associated(p%LC)) DEALLOCATE(P%LC);
    if(associated(p%TILTD)) DEALLOCATE(P%TILTD);
    !    if(associated(p%beta0)) DEALLOCATE(P%beta0);
    !    if(associated(p%gamma0I)) DEALLOCATE(P%gamma0I);
    !    if(associated(p%gambet)) DEALLOCATE(P%gambet);
    if(associated(p%P0C)) DEALLOCATE(P%P0C);
    if(associated(p%f)) then
       call kill(p%f)
       DEALLOCATE(P%f);
    endif
    if(associated(p%APERTURE)) then
       CALL kill(p%APERTURE)
       DEALLOCATE(p%APERTURE);
    endif
    if(associated(p%A)) then
       CALL KILL(P%A)
    ENDIF

    if(associated(p%EDGE))DEALLOCATE(P%EDGE);
    !    DEALLOCATE(P%TOTALPATH);
    if(associated(p%EXACT))DEALLOCATE(P%EXACT);
    !DEALLOCATE(P%RADIATION);DEALLOCATE(P%NOCAVITY);

    if(associated(p%PERMFRINGE))DEALLOCATE(P%PERMFRINGE);
    if(associated(p%KILL_ENT_FRINGE))DEALLOCATE(P%KILL_ENT_FRINGE);
    if(associated(p%KILL_EXI_FRINGE))DEALLOCATE(P%KILL_EXI_FRINGE);
    if(associated(p%bend_fringe))DEALLOCATE(P%bend_fringe); !DEALLOCATE(P%TIME);
    if(associated(p%METHOD))DEALLOCATE(P%METHOD);
    !DEALLOCATE(P%spin);
    if(associated(p%NST))DEALLOCATE(P%NST);
    if(associated(p%NMUL))DEALLOCATE(P%NMUL)
    !    CALL NULL_P(P)
    DEALLOCATE(P)
    nullify(p);
    ! if(junk) ccc=ccc-1
  end subroutine dealloc_p


  SUBROUTINE  KILL_S_APERTURE(A)
    implicit none
    TYPE(S_APERTURE), POINTER :: A(:)
    INTEGER I
    DO I=1,SIZE(A)
       CALL kill(A(I)%APERTURE)
       DEALLOCATE(A(I)%APERTURE)
    ENDDO
    DEALLOCATE(A)
  END SUBROUTINE  KILL_S_APERTURE

  SUBROUTINE  ALLOC_S_APERTURE(A,N,APERTURE)
    implicit none
    TYPE(S_APERTURE), POINTER :: A(:)
    TYPE(MADX_APERTURE), OPTIONAL :: APERTURE
    INTEGER I,N
    ALLOCATE(A(N))
    DO I=1,SIZE(A)
       CALL ALLOC(A(I)%APERTURE)
       IF(PRESENT(APERTURE))A(I)%APERTURE=APERTURE
    ENDDO
  END SUBROUTINE  ALLOC_S_APERTURE

  SUBROUTINE  check_S_APERTURE_r(p,N,x)
    implicit none
    TYPE(magnet_chart), POINTER :: p
    real(dp),intent(inout):: x(6)
    INTEGER N
    if(p%dir==1) then
       call CHECK_APERTURE(p%a(n)%aperture,X)
    else
       call CHECK_APERTURE(p%a(p%nst+2-n)%aperture,X)
    endif

    !!   if(s_aperture_CHECK.and.associated(el%p%A)) call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2)
    !       if(associated(T%PARENT_FIBRE%MAG%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
  END SUBROUTINE  check_S_APERTURE_r

  SUBROUTINE  check_S_APERTURE_p(p,N,x)
    implicit none
    TYPE(magnet_chart), POINTER :: p
    type(real_8),intent(inout):: x(6)
    INTEGER N

    if(p%dir==1) then
       call CHECK_APERTURE(p%a(n)%aperture,X)
    else
       call CHECK_APERTURE(p%a(p%nst+2-n)%aperture,X)
    endif

    !!   if(s_aperture_CHECK.and.associated(el%p%A)) call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2)
    !       if(associated(T%PARENT_FIBRE%MAG%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
  END SUBROUTINE  check_S_APERTURE_p

  SUBROUTINE  check_S_APERTURE_out_r(p,N,x)
    implicit none
    TYPE(magnet_chart), POINTER :: p
    real(dp),intent(inout):: x(6)
    INTEGER N
    if(p%dir==1.and.n==p%nst) then
       call CHECK_APERTURE(p%a(n+1)%aperture,X)
    elseif(p%dir==-1.and.n==p%nst) then
       call CHECK_APERTURE(p%a(1)%aperture,X)
    endif
  END SUBROUTINE  check_S_APERTURE_out_r

  SUBROUTINE  check_S_APERTURE_out_p(p,N,x)
    implicit none
    TYPE(magnet_chart), POINTER :: p
    type(real_8),intent(inout):: x(6)
    INTEGER N
    if(p%dir==1.and.n==p%nst) then
       call CHECK_APERTURE(p%a(n+1)%aperture,X)
    elseif(p%dir==-1.and.n==p%nst) then
       call CHECK_APERTURE(p%a(1)%aperture,X)
    endif
  END SUBROUTINE  check_S_APERTURE_out_p

  SUBROUTINE  equal_A(elp,el)
    implicit none
    type (MADX_APERTURE),INTENT(inOUT)::elP
    type (MADX_APERTURE),INTENT(IN)::el
    ELP%KIND=EL%KIND
    ELP%R=EL%R
    ELP%X=EL%X
    ELP%Y=EL%Y
    ELP%DX=EL%DX
    ELP%DY=EL%DY
  END SUBROUTINE  equal_A



  SUBROUTINE  equal_p(elp,el)
    implicit none
    type (MAGNET_CHART),INTENT(inOUT)::elP
    type (MAGNET_CHART),INTENT(IN)::el

    !    elp%beta0 =el%beta0
    !    elp%gamma0I=el%gamma0I
    !    elp%gambet =el%gambet
    elp%P0C =el%P0C
    elp%EXACT=el%EXACT
    !   elp%RADIATION=el%RADIATION
    !   elp%TIME=el%TIME
    !   elp%NOCAVITY=el%NOCAVITY
    !   elp%spin=el%spin
       elp%permFRINGE=el%permFRINGE
    elp%KILL_ENT_FRINGE=el%KILL_ENT_FRINGE
    elp%KILL_EXI_FRINGE=el%KILL_EXI_FRINGE
    elp%bend_fringe=el%bend_fringe

    elp%LD=el%LD
    elp%LC=el%LC
    elp%TILTD=el%TILTD
    elp%B0=el%B0
    elp%EDGE(1)=el%EDGE(1)
    elp%EDGE(2)=el%EDGE(2)
    elp%METHOD=el%METHOD
    elp%NST=el%NST
    !   ELP%TOTALPATH=EL%TOTALPATH
    ELP%nmul=EL%nmul
    !    ELP%TRACK=EL%TRACK

    if(associated(el%f)) then
       elp%f=el%f
    endif
    if(associated(el%APERTURE)) then
       if(.not.associated(elp%APERTURE)) call alloc(elp%APERTURE)
       elp%APERTURE=el%APERTURE
    endif

  end SUBROUTINE  equal_p

  SUBROUTINE  CHECK_APERTURE_R(E,X)
    implicit none
    type (MADX_APERTURE),INTENT(IN)::E
    REAL(DP), INTENT(IN):: X(6)
    !    real(dp) xx,yy,dx,dy,ex,ey

    !  real(dp) :: xlost(6)=zero
    !  character(120) :: messagelost

    IF(CHECK_MADX_APERTURE.AND.APERTURE_FLAG) THEN
       SELECT CASE(E%KIND)
       CASE(1)  ! ellipse circles
          IF((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2/E%R(2)**2>1.0_dp) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             messagelost="Lost in real kind=1 elliptic Aperture"
          ENDIF
       CASE(2)  ! rectangle
          IF(ABS(X(1)-E%DX)>E%X.OR.ABS(X(3)-E%DY)>E%Y) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             messagelost="Lost in real kind=2 rectangular Aperture"
          ENDIF
       CASE(3)  ! RECTANGLE + ELLIPSE (CIRCLE)
          IF((ABS(X(1)-E%DX)>E%X).OR.(ABS(X(3)-E%DY)>E%Y).OR.  &
               ((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2**2/E%R(2)**2>1.0_dp)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             messagelost="Lost in real kind=3 rect-ellipse Aperture"
          ENDIF
       CASE(4) ! MARGUERITE
          IF(((X(1)-E%DX)**2/E%R(2)**2+(X(3)-E%DY)**2/E%R(1)**2>1.0_dp).OR.  &
               ((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2/E%R(2)**2>1.0_dp)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             messagelost="Lost in real kind=4 marguerite Aperture"
          ENDIF
       CASE(5) ! RACETRACK
          IF( (abs(x(1)-e%dx)) > (e%r(1)+e%x)                  &
               .or. abs(x(3)-e%dy) .gt. (e%y+e%r(1)) .or.                &
               ((((abs(x(1)-e%dx)-e%x)**2+                            &
               (abs(x(3)-e%dy)-e%y)**2) .gt. e%r(1)**2)                  &
               .and. (abs(x(1)-e%dx)) .gt. e%x                        &
               .and. abs(x(3)-e%dy) .gt. e%y)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             messagelost="Lost in real kind=5 racetrack Aperture"
          ENDIF

       CASE(6) ! PILES OF POINTS
          STOP 222
       CASE DEFAULT
          !   STOP 223
       END SELECT

    ENDIF


  END SUBROUTINE  CHECK_APERTURE_R

  SUBROUTINE  CHECK_APERTURE_P(E,X)
    implicit none
    type (MADX_APERTURE),INTENT(IN)::E
    TYPE(REAL_8), INTENT(IN):: X(6)
    REAL(DP) Y(6)

    Y=X
    CALL CHECK_APERTURE(E,Y)

  END SUBROUTINE  CHECK_APERTURE_P

  FUNCTION minu( S1,S2  )
    implicit none
    logical(lp) minu
    logical(lp), INTENT (IN) :: S1
    logical(lp), INTENT (IN) :: S2

    minu=.false.
    if(s1.and.(.not.s2)) minu=my_true

  END FUNCTION minu



  SUBROUTINE  EQUALTILT(S2,S1)
    implicit none
    INTEGER I
    type (TILTING),INTENT(OUT)::S2
    type (TILTING),INTENT(IN)::S1

    S2%NATURAL=S1%NATURAL
    DO I=0,NMAX
       S2%TILT(I)=S1%TILT(I)
    ENDDO

  END SUBROUTINE EQUALTILT


  SUBROUTINE MAKE_STATES_0(particle)
    USE   definition
    USE   da_arrays
    IMPLICIT NONE
    LOGICAL(lp) particle
    integer i,lda_old

    W_P=>W_I
    NULLIFY(ACC);       
    NULLIFY(ACCfirst);       
    NULLIFY(paccfirst);       
    NULLIFY(paccthen);       


    insane_PTC=.true.
    NEWFILE%MF=.FALSE.
    CLOSEFILE%MF=.FALSE.
    KNOB=.FALSE.
    SETKNOB=.TRUE.

    CALL MAKE_YOSHIDA
    MYTYPE(0)=" MARKER"
    MYTYPE(kind1)=" DRIFT"
    MYTYPE(kind2)=" DRIFT-KICK-DRIFT"
    MYTYPE(kind3)=" THIN ELEMENT"
    MYTYPE(kind4)=" RF CAVITY"
    MYTYPE(kind5)=" SOLENOID "
    MYTYPE(kind6)=" KICK-SixTrack-KICK "
    MYTYPE(kind7)=" MATRIX-KICK-MATRIX "
    MYTYPE(kind8)=" NORMAL SMI "
    MYTYPE(kind9)=" SKEW   SMI "
    MYTYPE(kind10)=" EXACT SECTOR "
    MYTYPE(kind11)=" MONITOR "
    MYTYPE(kind12)=" HORIZONTAL MONITOR "
    MYTYPE(kind13)=" VERTICAL MONITOR "
    MYTYPE(kind14)=" INSTRUMENT "
    MYTYPE(kind15)=" ELECTRIC SEPTUM "
    !               123456789012345678901234
    MYTYPE(kind16)=" TRUE PARAELLEL  BEND  "
    MYTYPE(kind20)=" STRAIGHT EXACT (BEND) "
    MYTYPE(kind17)=" SOLENOID SIXTRACK"
    !    MYTYPE(KINDFITTED)=" FITTED "
    !    MYTYPE(KINDUSER1)=" USER_1 "
    !    MYTYPE(KINDUSER2)=" USER_2 "
    ind_stoc(1)='100000'
    ind_stoc(2)='010000'
    ind_stoc(3)='001000'
    ind_stoc(4)='000100'
    ind_stoc(5)='000010'
    ind_stoc(6)='000001'
    MADTHICK=>MADKIND2
    MADTHIN_NORMAL=>MADKIND3N
    MADTHIN_SKEW=>MADKIND3S

    MADFAC=1.0_dp   ! to prevent overflow (David Sagan)
    DO I=2,NMAX
       MADFAC(I)=(I-1)*MADFAC(I-1)
    ENDDO

    !    w_p=0
    !    w_p%nc=1
    !    w_p%fc='(2((1X,A72,/)))'
    !    w_p%c(1) = " RBEND ARE READ DIFFERENTLY DEPENDING ON VARIABLE MADLENGTH"
    !    IF(MADLENGTH) THEN
    !       w_p%c(2) = " INPUT OF CARTESIAN LENGTH FOR RBEND  "
    !    ELSE
    !       w_p%c(2) = " INPUT OF POLAR LENGTH FOR RBEND  "
    !    ENDIF
    !    ! call ! WRITE_I

    NSTD=1
    METD=2
    electron=PARTICLE
    if(electron) then
       A_particle = A_ELECTRON
    else
       A_particle = A_PROTON
    endif
    !    w_p=0
    !    w_p%nc=1
    !    w_p%fc='(1((1X,A72)))'
    !    IF(ELECTRON) THEN
    !       w_p%c(1) = " THIS IS A ELECTRON (POSITRON ACTUALLY) "
    !    ELSE
    !       w_p%c(1) = " THIS IS A PROTON "
    !    ENDIF
    tilt%natural=.true.
    tilt%tilt(0)=0.0_dp
    do i=1,nmax
       tilt%tilt(i)=pih/i
    enddo
    !  SECTOR_B AND SECTOR_NMUL FOR TYPE TEAPOT

    change_sector=my_false

    IF(SECTOR_NMUL>0.and.firsttime_coef) THEN
     call alloc(e_muon_scale)
     e_muon_scale=1.0_dp
       !  verb=global_verbose
       !  global_verbose=.false.
       if(firsttime_coef.or.(.not.allocated(S_B))) then
       
         if(   SECTOR_NMUL==10.and. SECTOR_NMUL_MAX==10) then
             call set_s_b
         else 
         
         
         
          !          SECTOR_B%firsttime=0   !slightly unsafe
          ALLOCATE(S_B(SECTOR_NMUL_MAX))
          lda_old=lda_used
          lda_used=3000
       !        ALLOCATE(S_B0)
       !        S_B0%firsttime=0
       !        call nul_coef(S_B0)
       !        call make_coef(S_B0,SECTOR_NMUL_MAX)
       !        call get_bend_coeff(S_B0,SECTOR_NMUL_MAX)
          DO I=1,SECTOR_NMUL_MAX
             !             if(i==SECTOR_NMUL_MAX)     global_verbose=.true.
             S_B(I)%firsttime=0
             call nul_coef(S_B(I))
             call make_coef(S_B(I),I)
             call get_bend_coeff(S_B(I),I)
          ENDDO
          lda_used=lda_old
          call print_curv("Maxwellian_bend_for_ptc.txt")
       endif

       firsttime_coef=.FALSE.
    endif  !  calling preset routine
    else
         call kill(e_muon_scale)
         call alloc(e_muon_scale)
     e_muon_scale=1.0_dp
    ENDIF
    call clear_states
    !  global_verbose=verb

  END  SUBROUTINE MAKE_STATES_0
  
 SUBROUTINE print_curv(filename)
    IMPLICIT NONE
    INTEGER I,J,nmul,mf
    character(*) filename
    character(255) line
    call kanalnummer(mf,filename)
    
    
    nmul=SECTOR_NMUL
    write(mf,*) SECTOR_NMUL ,SECTOR_NMUL_MAX
         DO J=1,S_B(NMUL)%N_MONO
         write(line,*) "S_B(" ,NMUL, ")%I(",J,")=",S_B(NMUL)%I(J),";","S_B(" ,NMUL, ")%J(",J,")=",S_B(NMUL)%J(J),";"
         call context(line)
         write(mf,*) line(1:len_trim(line))
         enddo
    DO I=1,NMUL
       DO J=1,S_B(NMUL)%N_MONO
        if(S_B(NMUL)%A_X(I,J)/=0.0_dp) then
         write(line,*) "S_B(" , NMUL , ")%A_X(" ,I,",",J, ")=",S_B(NMUL)%A_X(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B(NMUL)%B_X(I,J)/=0.0_dp) then
         write(line,*) "S_B(" , NMUL , ")%B_X(" ,I,",",J, ")=",S_B(NMUL)%B_X(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B(NMUL)%A_y(I,J)/=0.0_dp) then
         write(line,*) "S_B(" , NMUL , ")%A_y(" ,I,",",J, ")=",S_B(NMUL)%A_y(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B(NMUL)%B_y(I,J)/=0.0_dp) then
         write(line,*) "S_B(" , NMUL , ")%B_y(" ,I,",",J, ")=",S_B(NMUL)%B_y(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
       ENDDO
    ENDDO
    
    close(mf)

  END SUBROUTINE print_curv

  subroutine make_set_coef(b,no)
    implicit none
    integer no
    integer i,j,k,a,m , ic

    type(B_CYL) b

    ic=1
    b%firsttime=-100
    allocate(b%nmul)
    allocate(b%n_mono)
    b%nmul=no
    b%n_mono=((no+2-ic)*(no+1-ic))/2
    allocate(b%i(b%n_mono),b%j(b%n_mono))
    allocate(b%a_x(no,b%n_mono),b%a_y(no,b%n_mono))
    allocate(b%b_x(no,b%n_mono),b%b_y(no,b%n_mono))
    b%i=0
    b%j=0
    b%a_x=0.0_dp
    b%b_x=0.0_dp
    b%a_y=0.0_dp
    b%b_y=0.0_dp
  end subroutine make_set_coef
  
 
  subroutine clear_states     !%nxyz
    implicit none
    DEFAULT=DEFAULT0
    TOTALPATH=TOTALPATH0
    RADIATION=RADIATION0
    NOCAVITY=NOCAVITY0
    ENVELOPE=ENVELOPE0
    FRINGE=FRINGE0
    TIME=TIME0
    ! EXACTMIS=EXACTMIS0
    !  exactmis%exactmis=exactmis%exactmis.and.(.not.sixtrack_compatible)
    ONLY_4D=ONLY_4d0
    DELTA=DELTA0
    SPIN=SPIN0
    MODULATION=MODULATION0
  end subroutine clear_states  !%nxyz

  subroutine print_s(S,MF)
    implicit none
    type (INTERNAL_STATE) S
    INTEGER MF


    write(mf,*) "************ State Summary ****************"
    write(mf,'((1X,a16,1x,i4,1x,a24))' ) "MADTHICK=>KIND =", MADKIND2,Mytype(MADKIND2)
    if(MADLENGTH) then
       write(mf,*) ' Rectangular Bend: input cartesian length '
    else
       write(mf,*) ' Rectangular Bend: input arc length (rho alpha) '
    endif
    write(mf,'((1X,a28,1x,i4))' ) ' Default integration method ',METD
    write(mf,'((1X,a28,1x,i4))' ) ' Default integration steps  ',NSTD
 !   if(CAVITY_TOTALPATH==1) write(mf,'((1X,a24))' ) ' Real Pill Box Cavities '
 !   if(CAVITY_TOTALPATH==0) write(mf,'((1X,a24))' ) ' Fake Pill Box Cavities '

    If(electron) then
       if(muon==1.0_dp)  then
          write(mf,*)"This is an electron (positron actually if charge=1) "
       else
          write(mf,'((1X,a21,1x,G21.14,1x,A24))' ) "This a particle with ",muon, "times the electron mass "
       endif
    else
       write(mf,*) "This is a proton "
    endif
    write(mf, '((1X,a20,1x,a5))' )  "      EXACT_MODEL = ", CONV(EXACT_MODEL    )
    write(mf, '((1X,a20,1x,i4))' )  "      TOTALPATH   = ", S%TOTALPATH
    !    write(mf, '((1X,a20,1x,a5))' )  "      EXACTMIS    = ", CONV(S%EXACTMIS    )
    write(mf,'((1X,a20,1x,a5))' ) "      RADIATION   = ", CONV(S%RADIATION  )
    write(mf,'((1X,a20,1x,a5))' ) "      STOCHASTIC  = ", CONV(S%STOCHASTIC  )
    write(mf,'((1X,a20,1x,a5))' ) "      ENVELOPE    = ", CONV(S%ENVELOPE  )
    write(mf,'((1X,a20,1x,a5))' ) "      NOCAVITY    = ", CONV(S%NOCAVITY )
    write(mf,'((1X,a20,1x,a5))' ) "      TIME        = ", CONV(S%TIME )
    write(mf,'((1X,a20,1x,a5))' ) "      FRINGE      = ", CONV(S%FRINGE   )
    write(mf,'((1X,a20,1x,a5))' ) "      PARA_IN     = ", CONV(S%PARA_IN  )
    write(mf,'((1X,a20,1x,a5))' ) "      ONLY_4D     = ", CONV(S%ONLY_4D   )
    write(mf,'((1X,a20,1x,a5))' ) "      DELTA       = ", CONV(S%DELTA    )
    write(mf,'((1X,a20,1x,a5))' ) "      SPIN        = ", CONV(S%SPIN    )
    write(mf,'((1X,a20,1x,a5))' ) "      MODULATION  = ", CONV(S%MODULATION    )
    write(mf, '((1X,a20,1x,a5))' )"      RAMPING     = "  , CONV(ramp    )
    write(mf, '((1X,a20,1x,a5))' )"      ACCELERATE  = "  , CONV(ACCELERATE    )
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

  SUBROUTINE MAKE_STATES_m(muonfactor)
    USE   definition
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    real(dp) muonfactor,MASSF
    CALL MAKE_YOSHIDA
    muon=muonfactor
    MASSF=muon*pmae
    call MAKE_STATES_0(doneitt)

    IF(ABS(MASSF-pmap)/PMAP<0.01E0_DP) THEN
       A_PARTICLE=A_PROTON
    ELSEIF(ABS(MASSF-pmae)/pmae<0.01E0_DP) THEN
       A_PARTICLE=A_ELECTRON
    ELSEIF(ABS(MASSF-pmaMUON)/pmaMUON<0.01E0_DP) THEN
       A_PARTICLE=A_MUON
    ENDIF

  END  SUBROUTINE MAKE_STATES_m

  SUBROUTINE update_STATES
    USE   definition
    IMPLICIT NONE

    TOTALPATH=  TOTALPATH+DEFAULT

    RADIATION=  RADIATION+DEFAULT

    NOCAVITY=  NOCAVITY+DEFAULT

    STOCHASTIC=  STOCHASTIC+DEFAULT

    ENVELOPE=  ENVELOPE+DEFAULT

    TIME=  TIME+DEFAULT

    FRINGE=  FRINGE+DEFAULT

    ONLY_4D= ONLY_4D+DEFAULT

    DELTA= DELTA+DEFAULT

    ! EXACTMIS= EXACTMIS+DEFAULT

    SPIN= SPIN+DEFAULT

    MODULATION= MODULATION+DEFAULT

  END  SUBROUTINE update_STATES


  SUBROUTINE  EQUALt(S2,S1)
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
    S2%ONLY_4D=      S1%ONLY_4D
    S2%DELTA=       S1%DELTA
    S2%SPIN=       S1%SPIN
    S2%MODULATION=       S1%MODULATION
    !    S2%spin_dim=       S1%spin_dim
  END SUBROUTINE EQUALt



  FUNCTION add( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) add
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2


    add%TOTALPATH=0
    if((S1%TOTALPATH==1).OR.(S2%TOTALPATH==1)) add%TOTALPATH=1

    !    add%EXACT    =       S1%EXACT.OR.S2%EXACT
    !   add%EXACTMIS    =       (S1%EXACTMIS.OR.S2%EXACTMIS).and.(.not.sixtrack_compatible)
    add%RADIATION  =  S1%RADIATION.OR.S2%RADIATION
    add%NOCAVITY =  S1%NOCAVITY.OR.S2%NOCAVITY
    add%TIME     =  S1%TIME.OR.S2%TIME
    add%FRINGE   =       S1%FRINGE.OR.S2%FRINGE
    add%stochastic   =       S1%stochastic.OR.S2%stochastic
    add%ENVELOPE   =       S1%ENVELOPE.OR.S2%ENVELOPE
    add%ONLY_4D  =       S1%ONLY_4D.OR.S2%ONLY_4D
    add%DELTA  =       S1%DELTA.OR.S2%DELTA
    add%SPIN  =       S1%SPIN.OR.S2%SPIN
    add%MODULATION  =       S1%MODULATION.OR.S2%MODULATION
    add%PARA_IN  =       S1%PARA_IN.OR.S2%PARA_IN.or.ALWAYS_knobs
    !    add%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(add%stochastic) THEN
       add%RADIATION=T
    ENDIF
  !  IF(add%ENVELOPE) THEN
  !     add%radiation=T
  !  ENDIF
    IF(add%stochastic) THEN
       add%radiation=T
    ENDIF
    IF(add%DELTA) THEN
       add%ONLY_4D=T
       add%NOCAVITY =  T
    ENDIF
    IF(add%ONLY_4D) THEN
       add%TOTALPATH=  0
       add%RADIATION  =  F
       add%NOCAVITY =  T
       add%stochastic   =  F
       add%ENVELOPE   =  F
    ENDIF
  END FUNCTION add

  FUNCTION sub( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) sub
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2
    logical(lp) dum1,dum2

    sub%TOTALPATH=0
    dum1=S1%TOTALPATH==1
    dum2=S2%TOTALPATH==1
    if(dum1.min.dum2) sub%TOTALPATH=1

    !    sub%TOTALPATH=  S1%TOTALPATH.min.S2%TOTALPATH

    !   sub%EXACTMIS    =       (S1%EXACTMIS.min.S2%EXACTMIS).and.(.not.sixtrack_compatible)
    sub%RADIATION  =  S1%RADIATION.min.S2%RADIATION
    sub%NOCAVITY =  S1%NOCAVITY.min.S2%NOCAVITY
    sub%TIME     =  S1%TIME.min.S2%TIME
    sub%FRINGE   =       S1%FRINGE.min.S2%FRINGE
    sub%stochastic   =       S1%stochastic.min.S2%stochastic
    sub%ENVELOPE   =       S1%ENVELOPE.min.S2%ENVELOPE
    sub%ONLY_4D  =       S1%ONLY_4D.min.S2%ONLY_4D
    sub%DELTA  =       S1%DELTA.min.S2%DELTA
    sub%SPIN  =       S1%SPIN.min.S2%SPIN
    sub%MODULATION  = S1%MODULATION.min.S2%MODULATION
    sub%PARA_IN  =       (S1%PARA_IN.MIN.S2%PARA_IN).or.ALWAYS_knobs
    !    sub%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(sub%stochastic) THEN
       sub%RADIATION=T
    ENDIF
 !   IF(sub%ENVELOPE) THEN
 !      sub%RADIATION=T
 !   ENDIF
    IF(sub%DELTA) THEN
       sub%ONLY_4D=T
       sub%NOCAVITY =  T
    ENDIF
    IF(sub%ONLY_4D) THEN
       sub%TOTALPATH=  0
       sub%RADIATION  =  F
       sub%ENVELOPE  =  F
       sub%stochastic  =  F
       sub%NOCAVITY =  T
       sub%stochastic   =  F
    ENDIF

  END FUNCTION sub

  FUNCTION PARA_REMA(S1)   ! UNARY +
    implicit none
    TYPE (INTERNAL_STATE) PARA_REMA
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1

    PARA_REMA         =    S1
    PARA_REMA%PARA_IN =     T

  END FUNCTION PARA_REMA



  subroutine S_init(STATE,NO1,NP1,pack,ND2,NPARA)
    !  subroutine S_init(STATE,NO1,NP1,PACKAGE,MAPINT,ND2,NPARA)
    implicit none
    TYPE (INTERNAL_STATE), INTENT(IN):: STATE
    LOGICAL(lp), optional, INTENT(IN):: pack
    INTEGER, INTENT(IN):: NO1,NP1
    INTEGER ND1,NDEL,NDPT1
    INTEGER,optional :: ND2,NPARA
    INTEGER  ND2l,NPARAl
    LOGICAL(lp) package
    call dd_p !valishev
    doing_ac_modulation_in_ptc=.false.
    package=my_true
    if(present(pack))     package=my_true


    NDEL=0
    NDPT1=0
    IF(STATE%NOCAVITY)  THEN
       IF(STATE%ONLY_4D) THEN
          IF(STATE%DELTA) THEN
             ND1=2
             NDEL=1
             !             MAPINT=5
          ELSE
             ND1=2
             NDEL=0
             !            MAPINT=4
          ENDIF
       ELSE
          ND1=3
          NDEL=0
          NDPT1=5  !+C_%NDPT_OTHER
          !         MAPINT=6
       ENDIF
    ELSE              ! CAVITY IN RING
       ND1=3
       NDPT1=0
       !       MAPINT=6
    ENDIF

    IF(STATE%modulation)  then
       doing_ac_modulation_in_ptc=.true.
       ND1=ND1+1
    endif

    !    write(6,*) NO1,ND1,NP1,NDEL,NDPT1
    !pause 678
    CALL INIT(NO1,ND1,NP1+NDEL,NDPT1,PACKAGE)

    ND2l=ND1*2
    NPARAl=ND2l+NDEL
    C_%NPARA=NPARAl
    C_%ND2=ND2l
    C_%npara_fpp=NPARAl
    C_%SPIN_POS=0
    C_%NSPIN=0

    if(present(nd2)) nd2=nd2l
    if(present(npara)) npara=nparal

  END  subroutine S_init


  subroutine init_default(STATE,NO1,NP1)
    implicit none
    TYPE (INTERNAL_STATE),OPTIONAL, INTENT(IN):: STATE
    INTEGER, OPTIONAL, INTENT(IN):: NO1,NP1
    INTEGER   ND2,NPARA,NO2,NP2
    TYPE (INTERNAL_STATE)   STATE2
    STATE2=DEFAULT
    NO2=1;NP2=0;
    IF(PRESENT(STATE))      STATE2=STATE
    IF(PRESENT(NO1))      NO2=NO1
    IF(PRESENT(NP1))      NP2=NP1
    call init(STATE2,NO2,NP2,my_true,ND2,NPARA)
    C_%NPARA=NPARA
    C_%ND2=ND2
    C_%npara_fpp=NPARA
  END  subroutine init_default

  SUBROUTINE B2PERPR(P,B,X,X5,B2)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: B2
    real(dp), INTENT(IN) :: X(6),B(3),X5
    TYPE(MAGNET_CHART), INTENT(IN) :: P
    real(dp) E(3)
    !---  GOOD FOR TIME=FALSE
    E(1)=P%DIR*X(2)/(1.0_dp+X5)
    E(2)=P%DIR*X(4)/(1.0_dp+X5)
    E(3)=P%DIR*ROOT((1.0_dp+X5)**2-X(2)**2-X(4)**2)/(1.0_dp+X5)

    B2=0.0_dp
    B2=(B(2)*E(3)-B(3)*E(2))**2+B2
    B2=(B(1)*E(2)-B(2)*E(1))**2+B2
    B2=(B(3)*E(1)-B(1)*E(3))**2+B2

    RETURN
  END       SUBROUTINE B2PERPR


  SUBROUTINE B2PERPP(P,B,X,X5,B2)
    IMPLICIT NONE
    TYPE(REAL_8) , INTENT(INOUT) :: B2
    TYPE(REAL_8) , INTENT(IN) :: X(6),B(3),X5
    TYPE(MAGNET_CHART), INTENT(IN) :: P
    TYPE(REAL_8)  E(3)
    CALL ALLOC(E,3)
    !---  GOOD FOR TIME=FALSE
    E(1)=P%DIR*X(2)/(1.0_dp+X5)
    E(2)=P%DIR*X(4)/(1.0_dp+X5)
    E(3)=P%DIR*SQRT((1.0_dp+X5)**2-X(2)**2-X(4)**2)/(1.0_dp+X5)


    B2=0.0_dp
    B2=(B(2)*E(3)-B(3)*E(2))**2+B2
    B2=(B(1)*E(2)-B(2)*E(1))**2+B2
    B2=(B(3)*E(1)-B(1)*E(3))**2+B2
    CALL KILL(E,3)
    RETURN
  END       SUBROUTINE B2PERPP

  SUBROUTINE DTILTR_EXTERNAL(DIR,TILTD,I,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    INTEGER,INTENT(IN):: I,DIR
    REAL(DP),INTENT(IN) :: TILTD
    real(dp) YS


    IF(TILTD==0.0_dp) RETURN
    IF(I==1) THEN
       ys=COS(DIR*TILTD)*x(1)+SIN(DIR*TILTD)*x(3)
       x(3)=COS(DIR*TILTD)*x(3)-SIN(DIR*TILTD)*x(1)
       x(1)=ys
       ys=COS(DIR*TILTD)*x(2)+SIN(DIR*TILTD)*x(4)
       x(4)=COS(DIR*TILTD)*x(4)-SIN(DIR*TILTD)*x(2)
       x(2)=ys
    ELSE
       ys=COS(DIR*TILTD)*x(1)-SIN(DIR*TILTD)*x(3)
       x(3)=COS(DIR*TILTD)*x(3)+SIN(DIR*TILTD)*x(1)
       x(1)=ys
       ys=COS(DIR*TILTD)*x(2)-SIN(DIR*TILTD)*x(4)
       x(4)=COS(DIR*TILTD)*x(4)+SIN(DIR*TILTD)*x(2)
       x(2)=ys
    ENDIF

  END SUBROUTINE DTILTR_EXTERNAL

  SUBROUTINE DTILTP_EXTERNAL(DIR,TILTD,I,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    INTEGER,INTENT(IN):: I,DIR
    REAL(DP),INTENT(IN) :: TILTD
    TYPE(REAL_8) YS

    IF(TILTD==0.0_dp) RETURN
    CALL ALLOC(YS)

    IF(I==1) THEN
       ys=COS(DIR*TILTD)*x(1)+SIN(DIR*TILTD)*x(3)
       x(3)=COS(DIR*TILTD)*x(3)-SIN(DIR*TILTD)*x(1)
       x(1)=ys
       ys=COS(DIR*TILTD)*x(2)+SIN(DIR*TILTD)*x(4)
       x(4)=COS(DIR*TILTD)*x(4)-SIN(DIR*TILTD)*x(2)
       x(2)=ys
    ELSE
       ys=COS(DIR*TILTD)*x(1)-SIN(DIR*TILTD)*x(3)
       x(3)=COS(DIR*TILTD)*x(3)+SIN(DIR*TILTD)*x(1)
       x(1)=ys
       ys=COS(DIR*TILTD)*x(2)-SIN(DIR*TILTD)*x(4)
       x(4)=COS(DIR*TILTD)*x(4)+SIN(DIR*TILTD)*x(2)
       x(2)=ys
    ENDIF
    CALL KILL(YS)

  END SUBROUTINE DTILTP_EXTERNAL

  SUBROUTINE dd_p   !(u,dd)    !valishev
    !use my_own_1D_TPSA

    ! Tracking subroutine for elliptical lens
    ! A.Valishev (valishev@fnal.gov) October 19, 2010
    ! Modified by E. Forest for PTC
    implicit none
    type(real_8) u,dd
    type(my_1D_taylor) del,logdel
    integer i,n

    call set_my_taylor_no(N_my_1D_taylor)
    del=0.0_dp
    del%a(1)=1.0_dp
    logdel=0.0_dp
    val_del=0.0_dp
    logdel=log(1.0_dp+del**2+del*sqrt(2.0_dp)*sqrt(1.0_dp+del**2/2.0_dp))

    do i=0,N_my_1D_taylor
       n=2*i+1
       if(n>N_my_1D_taylor) cycle
       val_del%a(i)=logdel%a(n)
    enddo
  end SUBROUTINE dd_p  !valishev

subroutine set_s_b
  implicit none
  integer i,s1,s2
 
  s1=SECTOR_NMUL_MAX
  s2=SECTOR_NMUL
  SECTOR_NMUL_MAX=10
  SECTOR_NMUL=10

  ALLOCATE(S_B(SECTOR_NMUL_MAX))
           DO I=1,SECTOR_NMUL_MAX
             call nul_coef(S_B(I))
             call make_set_coef(S_B(I),I)
             S_B(I)%firsttime=0
          ENDDO

 S_B(10)%I(1)=9;S_B(10)%J(1)=0;
 S_B(10)%I(2)=8;S_B(10)%J(2)=1;
 S_B(10)%I(3)=8;S_B(10)%J(3)=0;
 S_B(10)%I(4)=7;S_B(10)%J(4)=2;
 S_B(10)%I(5)=7;S_B(10)%J(5)=1;
 S_B(10)%I(6)=7;S_B(10)%J(6)=0;
 S_B(10)%I(7)=6;S_B(10)%J(7)=3;
 S_B(10)%I(8)=6;S_B(10)%J(8)=2;
 S_B(10)%I(9)=6;S_B(10)%J(9)=1;
 S_B(10)%I(10)=6;S_B(10)%J(10)=0;
 S_B(10)%I(11)=5;S_B(10)%J(11)=4;
 S_B(10)%I(12)=5;S_B(10)%J(12)=3;
 S_B(10)%I(13)=5;S_B(10)%J(13)=2;
 S_B(10)%I(14)=5;S_B(10)%J(14)=1;
 S_B(10)%I(15)=5;S_B(10)%J(15)=0;
 S_B(10)%I(16)=4;S_B(10)%J(16)=5;
 S_B(10)%I(17)=4;S_B(10)%J(17)=4;
 S_B(10)%I(18)=4;S_B(10)%J(18)=3;
 S_B(10)%I(19)=4;S_B(10)%J(19)=2;
 S_B(10)%I(20)=4;S_B(10)%J(20)=1;
 S_B(10)%I(21)=4;S_B(10)%J(21)=0;
 S_B(10)%I(22)=3;S_B(10)%J(22)=6;
 S_B(10)%I(23)=3;S_B(10)%J(23)=5;
 S_B(10)%I(24)=3;S_B(10)%J(24)=4;
 S_B(10)%I(25)=3;S_B(10)%J(25)=3;
 S_B(10)%I(26)=3;S_B(10)%J(26)=2;
 S_B(10)%I(27)=3;S_B(10)%J(27)=1;
 S_B(10)%I(28)=3;S_B(10)%J(28)=0;
 S_B(10)%I(29)=2;S_B(10)%J(29)=7;
 S_B(10)%I(30)=2;S_B(10)%J(30)=6;
 S_B(10)%I(31)=2;S_B(10)%J(31)=5;
 S_B(10)%I(32)=2;S_B(10)%J(32)=4;
 S_B(10)%I(33)=2;S_B(10)%J(33)=3;
 S_B(10)%I(34)=2;S_B(10)%J(34)=2;
 S_B(10)%I(35)=2;S_B(10)%J(35)=1;
 S_B(10)%I(36)=2;S_B(10)%J(36)=0;
 S_B(10)%I(37)=1;S_B(10)%J(37)=8;
 S_B(10)%I(38)=1;S_B(10)%J(38)=7;
 S_B(10)%I(39)=1;S_B(10)%J(39)=6;
 S_B(10)%I(40)=1;S_B(10)%J(40)=5;
 S_B(10)%I(41)=1;S_B(10)%J(41)=4;
 S_B(10)%I(42)=1;S_B(10)%J(42)=3;
 S_B(10)%I(43)=1;S_B(10)%J(43)=2;
 S_B(10)%I(44)=1;S_B(10)%J(44)=1;
 S_B(10)%I(45)=1;S_B(10)%J(45)=0;
 S_B(10)%I(46)=0;S_B(10)%J(46)=9;
 S_B(10)%I(47)=0;S_B(10)%J(47)=8;
 S_B(10)%I(48)=0;S_B(10)%J(48)=7;
 S_B(10)%I(49)=0;S_B(10)%J(49)=6;
 S_B(10)%I(50)=0;S_B(10)%J(50)=5;
 S_B(10)%I(51)=0;S_B(10)%J(51)=4;
 S_B(10)%I(52)=0;S_B(10)%J(52)=3;
 S_B(10)%I(53)=0;S_B(10)%J(53)=2;
 S_B(10)%I(54)=0;S_B(10)%J(54)=1;
 S_B(10)%I(55)=0;S_B(10)%J(55)=0;
 S_B(10)%A_Y(1,4)=-0.500000000000000E0_DP
 S_B(10)%A_X(1,7)=-1.16666666666667E0_DP
 S_B(10)%A_Y(1,8)=0.499999999999999E0_DP
 S_B(10)%A_Y(1,11)=2.62500000000000E0_DP
 S_B(10)%A_X(1,12)=0.999999999999999E0_DP
 S_B(10)%A_Y(1,13)=-0.500000000000000E0_DP
 S_B(10)%A_X(1,16)=2.62500000000000E0_DP
 S_B(10)%A_Y(1,17)=-1.87500000000000E0_DP
 S_B(10)%A_X(1,18)=-0.833333333333334E0_DP
 S_B(10)%A_Y(1,19)=0.500000000000000E0_DP
 S_B(10)%A_Y(1,22)=-2.18750000000000E0_DP
 S_B(10)%A_X(1,23)=-1.50000000000000E0_DP
 S_B(10)%A_Y(1,24)=1.25000000000000E0_DP
 S_B(10)%A_X(1,25)=0.666666666666667E0_DP
 S_B(10)%A_Y(1,26)=-0.500000000000000E0_DP
 S_B(10)%A_X(1,29)=-0.937500000000000E0_DP
 S_B(10)%A_Y(1,30)=0.937500000000000E0_DP
 S_B(10)%A_X(1,31)=0.750000000000000E0_DP
 S_B(10)%A_Y(1,32)=-0.750000000000000E0_DP
 S_B(10)%A_X(1,33)=-0.500000000000000E0_DP
 S_B(10)%A_Y(1,34)=0.500000000000000E0_DP
 S_B(10)%A_Y(1,37)=0.273437500000000E0_DP
 S_B(10)%A_X(1,38)=0.267857142857143E0_DP
 S_B(10)%A_Y(1,39)=-0.312500000000000E0_DP
 S_B(10)%A_X(1,40)=-0.300000000000000E0_DP
 S_B(10)%A_Y(1,41)=0.375000000000000E0_DP
 S_B(10)%A_X(1,42)=0.333333333333333E0_DP
 S_B(10)%A_Y(1,43)=-0.500000000000000E0_DP
 S_B(10)%B_X(1,45)=-1.00000000000000E0_DP
 S_B(10)%A_Y(1,45)=1.00000000000000E0_DP
 S_B(10)%A_X(1,46)=3.038194444444443E-002_DP
 S_B(10)%A_Y(1,47)=-3.906249999999999E-002_DP
 S_B(10)%A_X(1,48)=-4.464285714285716E-002_DP
 S_B(10)%A_Y(1,49)=6.250000000000000E-002_DP
 S_B(10)%A_X(1,50)=7.500000000000000E-002_DP
 S_B(10)%A_Y(1,51)=-0.125000000000000E0_DP
 S_B(10)%A_X(1,52)=-0.166666666666667E0_DP
 S_B(10)%A_Y(1,53)=0.500000000000000E0_DP
 S_B(10)%A_X(1,54)=1.00000000000000E0_DP
 S_B(10)%B_X(1,55)=-1.00000000000000E0_DP
 S_B(10)%A_Y(1,55)=1.00000000000000E0_DP
 S_B(10)%A_Y(2,4)=0.500000000000001E0_DP
 S_B(10)%A_X(2,7)=1.16666666666667E0_DP
 S_B(10)%B_Y(2,7)=0.166666666666667E0_DP
 S_B(10)%A_Y(2,8)=-0.499999999999998E0_DP
 S_B(10)%B_X(2,11)=0.250000000000001E0_DP
 S_B(10)%A_Y(2,11)=-2.62500000000000E0_DP
 S_B(10)%A_X(2,12)=-0.999999999999996E0_DP
 S_B(10)%B_Y(2,12)=-0.166666666666667E0_DP
 S_B(10)%A_Y(2,13)=0.500000000000001E0_DP
 S_B(10)%A_X(2,16)=-2.62500000000000E0_DP
 S_B(10)%B_Y(2,16)=-0.375000000000001E0_DP
 S_B(10)%B_X(2,17)=-0.208333333333333E0_DP
 S_B(10)%A_Y(2,17)=1.87500000000000E0_DP
 S_B(10)%A_X(2,18)=0.833333333333335E0_DP
 S_B(10)%B_Y(2,18)=0.166666666666666E0_DP
 S_B(10)%A_Y(2,19)=-0.500000000000001E0_DP
 S_B(10)%B_X(2,22)=-0.250000000000001E0_DP
 S_B(10)%A_Y(2,22)=2.18750000000000E0_DP
 S_B(10)%A_X(2,23)=1.50000000000000E0_DP
 S_B(10)%B_Y(2,23)=0.250000000000000E0_DP
 S_B(10)%B_X(2,24)=0.166666666666666E0_DP
 S_B(10)%A_Y(2,24)=-1.25000000000000E0_DP
 S_B(10)%A_X(2,25)=-0.666666666666667E0_DP
 S_B(10)%B_Y(2,25)=-0.166666666666667E0_DP
 S_B(10)%A_Y(2,26)=0.499999999999999E0_DP
 S_B(10)%A_X(2,29)=0.937500000000000E0_DP
 S_B(10)%B_Y(2,29)=0.133928571428572E0_DP
 S_B(10)%B_X(2,30)=0.125000000000000E0_DP
 S_B(10)%A_Y(2,30)=-0.937499999999998E0_DP
 S_B(10)%A_X(2,31)=-0.750000000000001E0_DP
 S_B(10)%B_Y(2,31)=-0.150000000000000E0_DP
 S_B(10)%B_X(2,32)=-0.125000000000001E0_DP
 S_B(10)%A_Y(2,32)=0.750000000000001E0_DP
 S_B(10)%A_X(2,33)=0.499999999999999E0_DP
 S_B(10)%B_Y(2,33)=0.166666666666667E0_DP
 S_B(10)%A_Y(2,34)=-0.500000000000000E0_DP
 S_B(10)%B_X(2,36)=-1.00000000000000E0_DP
 S_B(10)%A_Y(2,36)=1.00000000000000E0_DP
 S_B(10)%B_X(2,37)=3.348214285714293E-002_DP
 S_B(10)%A_Y(2,37)=-0.273437500000000E0_DP
 S_B(10)%A_X(2,38)=-0.267857142857142E0_DP
 S_B(10)%B_Y(2,38)=-4.464285714285710E-002_DP
 S_B(10)%B_X(2,39)=-4.999999999999985E-002_DP
 S_B(10)%A_Y(2,39)=0.312500000000000E0_DP
 S_B(10)%A_X(2,40)=0.300000000000000E0_DP
 S_B(10)%B_Y(2,40)=7.500000000000023E-002_DP
 S_B(10)%B_X(2,41)=8.333333333333333E-002_DP
 S_B(10)%A_Y(2,41)=-0.375000000000000E0_DP
 S_B(10)%A_X(2,42)=-0.333333333333333E0_DP
 S_B(10)%B_Y(2,42)=-0.166666666666667E0_DP
 S_B(10)%A_Y(2,43)=0.500000000000000E0_DP
 S_B(10)%A_X(2,44)=2.00000000000000E0_DP
 S_B(10)%B_Y(2,44)=1.00000000000000E0_DP
 S_B(10)%B_X(2,45)=-1.00000000000000E0_DP
 S_B(10)%A_Y(2,45)=1.00000000000000E0_DP
 S_B(10)%A_X(2,46)=-3.038194444444445E-002_DP
 S_B(10)%B_Y(2,46)=-4.340277777777785E-003_DP
 S_B(10)%B_X(2,47)=-5.580357142857138E-003_DP
 S_B(10)%A_Y(2,47)=3.906249999999994E-002_DP
 S_B(10)%A_X(2,48)=4.464285714285721E-002_DP
 S_B(10)%B_Y(2,48)=8.928571428571412E-003_DP
 S_B(10)%B_X(2,49)=1.250000000000004E-002_DP
 S_B(10)%A_Y(2,49)=-6.250000000000003E-002_DP
 S_B(10)%A_X(2,50)=-7.499999999999996E-002_DP
 S_B(10)%B_Y(2,50)=-2.500000000000000E-002_DP
 S_B(10)%B_X(2,51)=-4.166666666666665E-002_DP
 S_B(10)%A_Y(2,51)=0.125000000000000E0_DP
 S_B(10)%A_X(2,52)=0.166666666666667E0_DP
 S_B(10)%B_Y(2,52)=0.166666666666667E0_DP
 S_B(10)%B_X(2,53)=0.500000000000000E0_DP
 S_B(10)%A_Y(2,53)=-0.500000000000000E0_DP
 S_B(10)%A_X(2,54)=1.00000000000000E0_DP
 S_B(10)%B_Y(2,54)=1.00000000000000E0_DP
 S_B(10)%A_Y(3,4)=-0.499999999999996E0_DP
 S_B(10)%A_X(3,7)=-1.16666666666666E0_DP
 S_B(10)%B_Y(3,7)=-0.333333333333335E0_DP
 S_B(10)%A_Y(3,8)=0.500000000000001E0_DP
 S_B(10)%B_X(3,11)=-0.500000000000002E0_DP
 S_B(10)%A_Y(3,11)=2.74999999999999E0_DP
 S_B(10)%A_X(3,12)=1.00000000000000E0_DP
 S_B(10)%B_Y(3,12)=0.333333333333334E0_DP
 S_B(10)%A_Y(3,13)=-0.500000000000000E0_DP
 S_B(10)%A_X(3,16)=2.74999999999999E0_DP
 S_B(10)%B_Y(3,16)=0.750000000000003E0_DP
 S_B(10)%B_X(3,17)=0.416666666666667E0_DP
 S_B(10)%A_Y(3,17)=-2.00000000000000E0_DP
 S_B(10)%A_X(3,18)=-0.833333333333334E0_DP
 S_B(10)%B_Y(3,18)=-0.333333333333336E0_DP
 S_B(10)%A_Y(3,19)=0.499999999999999E0_DP
 S_B(10)%B_X(3,22)=0.500000000000002E0_DP
 S_B(10)%A_Y(3,22)=-2.31249999999999E0_DP
 S_B(10)%A_X(3,23)=-1.60000000000000E0_DP
 S_B(10)%B_Y(3,23)=-0.500000000000001E0_DP
 S_B(10)%B_X(3,24)=-0.333333333333336E0_DP
 S_B(10)%A_Y(3,24)=1.37500000000000E0_DP
 S_B(10)%A_X(3,25)=0.666666666666665E0_DP
 S_B(10)%B_Y(3,25)=0.333333333333333E0_DP
 S_B(10)%A_Y(3,26)=-0.500000000000000E0_DP
 S_B(10)%B_X(3,28)=-1.00000000000000E0_DP
 S_B(10)%A_Y(3,28)=1.00000000000000E0_DP
 S_B(10)%A_X(3,29)=-0.991071428571425E0_DP
 S_B(10)%B_Y(3,29)=-0.267857142857144E0_DP
 S_B(10)%B_X(3,30)=-0.250000000000000E0_DP
 S_B(10)%A_Y(3,30)=1.01250000000000E0_DP
 S_B(10)%A_X(3,31)=0.825000000000000E0_DP
 S_B(10)%B_Y(3,31)=0.300000000000001E0_DP
 S_B(10)%B_X(3,32)=0.249999999999999E0_DP
 S_B(10)%A_Y(3,32)=-0.874999999999999E0_DP
 S_B(10)%A_X(3,33)=-0.500000000000000E0_DP
 S_B(10)%B_Y(3,33)=-0.333333333333334E0_DP
 S_B(10)%A_Y(3,34)=0.500000000000000E0_DP
 S_B(10)%A_X(3,35)=3.00000000000000E0_DP
 S_B(10)%B_Y(3,35)=2.00000000000000E0_DP
 S_B(10)%B_X(3,36)=-1.00000000000000E0_DP
 S_B(10)%A_Y(3,36)=1.00000000000000E0_DP
 S_B(10)%B_X(3,37)=-6.696428571428595E-002_DP
 S_B(10)%A_Y(3,37)=0.290178571428571E0_DP
 S_B(10)%A_X(3,38)=0.289285714285715E0_DP
 S_B(10)%B_Y(3,38)=8.928571428571447E-002_DP
 S_B(10)%B_X(3,39)=0.100000000000000E0_DP
 S_B(10)%A_Y(3,39)=-0.350000000000000E0_DP
 S_B(10)%A_X(3,40)=-0.349999999999999E0_DP
 S_B(10)%B_Y(3,40)=-0.150000000000000E0_DP
 S_B(10)%B_X(3,41)=-0.166666666666667E0_DP
 S_B(10)%A_Y(3,41)=0.500000000000000E0_DP
 S_B(10)%A_X(3,42)=0.333333333333333E0_DP
 S_B(10)%B_Y(3,42)=0.333333333333333E0_DP
 S_B(10)%B_X(3,43)=2.00000000000000E0_DP
 S_B(10)%A_Y(3,43)=-2.00000000000000E0_DP
 S_B(10)%A_X(3,44)=2.00000000000000E0_DP
 S_B(10)%B_Y(3,44)=2.00000000000000E0_DP
 S_B(10)%A_X(3,46)=3.224206349206341E-002_DP
 S_B(10)%B_Y(3,46)=8.680555555555582E-003_DP
 S_B(10)%B_X(3,47)=1.116071428571431E-002_DP
 S_B(10)%A_Y(3,47)=-4.241071428571432E-002_DP
 S_B(10)%A_X(3,48)=-5.000000000000000E-002_DP
 S_B(10)%B_Y(3,48)=-1.785714285714292E-002_DP
 S_B(10)%B_X(3,49)=-2.499999999999997E-002_DP
 S_B(10)%A_Y(3,49)=7.499999999999990E-002_DP
 S_B(10)%A_X(3,50)=0.100000000000000E0_DP
 S_B(10)%B_Y(3,50)=5.000000000000010E-002_DP
 S_B(10)%B_X(3,51)=8.333333333333336E-002_DP
 S_B(10)%A_Y(3,51)=-0.250000000000000E0_DP
 S_B(10)%A_X(3,52)=-0.666666666666667E0_DP
 S_B(10)%B_Y(3,52)=-0.333333333333333E0_DP
 S_B(10)%B_X(3,53)=1.00000000000000E0_DP
 S_B(10)%A_Y(3,53)=-1.00000000000000E0_DP
 S_B(10)%A_Y(4,4)=0.499999999999999E0_DP
 S_B(10)%A_X(4,7)=1.16666666666667E0_DP
 S_B(10)%B_Y(4,7)=0.499999999999997E0_DP
 S_B(10)%A_Y(4,8)=-0.500000000000000E0_DP
 S_B(10)%B_X(4,11)=0.749999999999995E0_DP
 S_B(10)%A_Y(4,11)=-3.00000000000000E0_DP
 S_B(10)%A_X(4,12)=-1.00000000000000E0_DP
 S_B(10)%B_Y(4,12)=-0.499999999999999E0_DP
 S_B(10)%A_Y(4,13)=0.499999999999999E0_DP
 S_B(10)%A_X(4,16)=-3.00000000000000E0_DP
 S_B(10)%B_Y(4,16)=-1.19999999999999E0_DP
 S_B(10)%B_X(4,17)=-0.624999999999999E0_DP
 S_B(10)%A_Y(4,17)=2.25000000000000E0_DP
 S_B(10)%A_X(4,18)=0.833333333333331E0_DP
 S_B(10)%B_Y(4,18)=0.499999999999999E0_DP
 S_B(10)%A_Y(4,19)=-0.500000000000001E0_DP
 S_B(10)%B_X(4,21)=-1.00000000000000E0_DP
 S_B(10)%A_Y(4,21)=1.00000000000000E0_DP
 S_B(10)%B_X(4,22)=-0.799999999999996E0_DP
 S_B(10)%A_Y(4,22)=2.56250000000000E0_DP
 S_B(10)%A_X(4,23)=1.80000000000000E0_DP
 S_B(10)%B_Y(4,23)=0.824999999999999E0_DP
 S_B(10)%B_X(4,24)=0.499999999999999E0_DP
 S_B(10)%A_Y(4,24)=-1.62500000000000E0_DP
 S_B(10)%A_X(4,25)=-0.666666666666668E0_DP
 S_B(10)%B_Y(4,25)=-0.500000000000000E0_DP
 S_B(10)%A_Y(4,26)=0.500000000000000E0_DP
 S_B(10)%A_X(4,27)=4.00000000000000E0_DP
 S_B(10)%B_Y(4,27)=3.00000000000000E0_DP
 S_B(10)%B_X(4,28)=-1.00000000000000E0_DP
 S_B(10)%A_Y(4,28)=1.00000000000000E0_DP
 S_B(10)%A_X(4,29)=1.09821428571429E0_DP
 S_B(10)%B_Y(4,29)=0.433928571428570E0_DP
 S_B(10)%B_X(4,30)=0.412500000000000E0_DP
 S_B(10)%A_Y(4,30)=-1.16250000000000E0_DP
 S_B(10)%A_X(4,31)=-0.974999999999999E0_DP
 S_B(10)%B_Y(4,31)=-0.524999999999999E0_DP
 S_B(10)%B_X(4,32)=-0.375000000000000E0_DP
 S_B(10)%A_Y(4,32)=1.12500000000000E0_DP
 S_B(10)%A_X(4,33)=0.500000000000000E0_DP
 S_B(10)%B_Y(4,33)=0.499999999999998E0_DP
 S_B(10)%B_X(4,34)=4.50000000000000E0_DP
 S_B(10)%A_Y(4,34)=-4.50000000000000E0_DP
 S_B(10)%A_X(4,35)=3.00000000000000E0_DP
 S_B(10)%B_Y(4,35)=3.00000000000000E0_DP
 S_B(10)%B_X(4,37)=0.108482142857142E0_DP
 S_B(10)%A_Y(4,37)=-0.323660714285714E0_DP
 S_B(10)%A_X(4,38)=-0.332142857142857E0_DP
 S_B(10)%B_Y(4,38)=-0.150000000000000E0_DP
 S_B(10)%B_X(4,39)=-0.175000000000000E0_DP
 S_B(10)%A_Y(4,39)=0.425000000000000E0_DP
 S_B(10)%A_X(4,40)=0.450000000000000E0_DP
 S_B(10)%B_Y(4,40)=0.300000000000000E0_DP
 S_B(10)%B_X(4,41)=0.249999999999999E0_DP
 S_B(10)%A_Y(4,41)=-0.750000000000000E0_DP
 S_B(10)%A_X(4,42)=-3.00000000000000E0_DP
 S_B(10)%B_Y(4,42)=-2.00000000000000E0_DP
 S_B(10)%B_X(4,43)=3.00000000000000E0_DP
 S_B(10)%A_Y(4,43)=-3.00000000000000E0_DP
 S_B(10)%A_X(4,46)=-3.596230158730157E-002_DP
 S_B(10)%B_Y(4,46)=-1.413690476190471E-002_DP
 S_B(10)%B_X(4,47)=-1.874999999999998E-002_DP
 S_B(10)%A_Y(4,47)=4.910714285714286E-002_DP
 S_B(10)%A_X(4,48)=6.071428571428568E-002_DP
 S_B(10)%B_Y(4,48)=3.214285714285710E-002_DP
 S_B(10)%B_X(4,49)=4.999999999999998E-002_DP
 S_B(10)%A_Y(4,49)=-0.100000000000000E0_DP
 S_B(10)%A_X(4,50)=-0.150000000000000E0_DP
 S_B(10)%B_Y(4,50)=-0.150000000000000E0_DP
 S_B(10)%B_X(4,51)=-0.500000000000000E0_DP
 S_B(10)%A_Y(4,51)=0.500000000000000E0_DP
 S_B(10)%A_X(4,52)=-1.00000000000000E0_DP
 S_B(10)%B_Y(4,52)=-1.00000000000000E0_DP
 S_B(10)%A_Y(5,4)=-0.500000000000002E0_DP
 S_B(10)%A_X(5,7)=-1.16666666666667E0_DP
 S_B(10)%B_Y(5,7)=-0.666666666666668E0_DP
 S_B(10)%A_Y(5,8)=0.500000000000000E0_DP
 S_B(10)%B_X(5,11)=-1.00000000000000E0_DP
 S_B(10)%A_Y(5,11)=3.37500000000001E0_DP
 S_B(10)%A_X(5,12)=1.00000000000000E0_DP
 S_B(10)%B_Y(5,12)=0.666666666666664E0_DP
 S_B(10)%A_Y(5,13)=-0.500000000000000E0_DP
 S_B(10)%B_X(5,15)=-1.00000000000000E0_DP
 S_B(10)%A_Y(5,15)=1.00000000000000E0_DP
 S_B(10)%A_X(5,16)=3.37500000000001E0_DP
 S_B(10)%B_Y(5,16)=1.80000000000000E0_DP
 S_B(10)%B_X(5,17)=0.833333333333330E0_DP
 S_B(10)%A_Y(5,17)=-2.62500000000000E0_DP
 S_B(10)%A_X(5,18)=-0.833333333333334E0_DP
 S_B(10)%B_Y(5,18)=-0.666666666666667E0_DP
 S_B(10)%A_Y(5,19)=0.500000000000000E0_DP
 S_B(10)%A_X(5,20)=5.00000000000000E0_DP
 S_B(10)%B_Y(5,20)=4.00000000000000E0_DP
 S_B(10)%B_X(5,21)=-1.00000000000000E0_DP
 S_B(10)%A_Y(5,21)=1.00000000000000E0_DP
 S_B(10)%B_X(5,22)=1.20000000000000E0_DP
 S_B(10)%A_Y(5,22)=-3.00000000000001E0_DP
 S_B(10)%A_X(5,23)=-2.10000000000000E0_DP
 S_B(10)%B_Y(5,23)=-1.30000000000000E0_DP
 S_B(10)%B_X(5,24)=-0.666666666666667E0_DP
 S_B(10)%A_Y(5,24)=2.00000000000000E0_DP
 S_B(10)%A_X(5,25)=0.666666666666667E0_DP
 S_B(10)%B_Y(5,25)=0.666666666666667E0_DP
 S_B(10)%B_X(5,26)=8.00000000000000E0_DP
 S_B(10)%A_Y(5,26)=-8.00000000000000E0_DP
 S_B(10)%A_X(5,27)=4.00000000000000E0_DP
 S_B(10)%B_Y(5,27)=4.00000000000000E0_DP
 S_B(10)%A_X(5,29)=-1.28571428571429E0_DP
 S_B(10)%B_Y(5,29)=-0.664285714285714E0_DP
 S_B(10)%B_X(5,30)=-0.649999999999999E0_DP
 S_B(10)%A_Y(5,30)=1.45000000000000E0_DP
 S_B(10)%A_X(5,31)=1.20000000000000E0_DP
 S_B(10)%B_Y(5,31)=0.900000000000000E0_DP
 S_B(10)%B_X(5,32)=0.500000000000000E0_DP
 S_B(10)%A_Y(5,32)=-1.50000000000000E0_DP
 S_B(10)%A_X(5,33)=-8.00000000000000E0_DP
 S_B(10)%B_Y(5,33)=-6.00000000000000E0_DP
 S_B(10)%B_X(5,34)=6.00000000000000E0_DP
 S_B(10)%A_Y(5,34)=-6.00000000000000E0_DP
 S_B(10)%B_X(5,37)=-0.166071428571429E0_DP
 S_B(10)%A_Y(5,37)=0.383928571428572E0_DP
 S_B(10)%A_X(5,38)=0.414285714285714E0_DP
 S_B(10)%B_Y(5,38)=0.242857142857143E0_DP
 S_B(10)%B_X(5,39)=0.300000000000000E0_DP
 S_B(10)%A_Y(5,39)=-0.600000000000000E0_DP
 S_B(10)%A_X(5,40)=-0.600000000000000E0_DP
 S_B(10)%B_Y(5,40)=-0.600000000000000E0_DP
 S_B(10)%B_X(5,41)=-3.00000000000000E0_DP
 S_B(10)%A_Y(5,41)=3.00000000000000E0_DP
 S_B(10)%A_X(5,42)=-4.00000000000000E0_DP
 S_B(10)%B_Y(5,42)=-4.00000000000000E0_DP
 S_B(10)%A_X(5,46)=4.265873015873022E-002_DP
 S_B(10)%B_Y(5,46)=2.182539682539682E-002_DP
 S_B(10)%B_X(5,47)=3.035714285714282E-002_DP
 S_B(10)%A_Y(5,47)=-6.250000000000000E-002_DP
 S_B(10)%A_X(5,48)=-8.571428571428574E-002_DP
 S_B(10)%B_Y(5,48)=-5.714285714285715E-002_DP
 S_B(10)%B_X(5,49)=-0.100000000000000E0_DP
 S_B(10)%A_Y(5,49)=0.200000000000000E0_DP
 S_B(10)%A_X(5,50)=0.600000000000000E0_DP
 S_B(10)%B_Y(5,50)=0.400000000000000E0_DP
 S_B(10)%B_X(5,51)=-1.00000000000000E0_DP
 S_B(10)%A_Y(5,51)=1.00000000000000E0_DP
 S_B(10)%A_Y(6,4)=0.499999999999998E0_DP
 S_B(10)%A_X(6,7)=1.16666666666666E0_DP
 S_B(10)%B_Y(6,7)=0.833333333333330E0_DP
 S_B(10)%A_Y(6,8)=-0.500000000000004E0_DP
 S_B(10)%B_X(6,10)=-1.00000000000000E0_DP
 S_B(10)%A_Y(6,10)=1.00000000000000E0_DP
 S_B(10)%B_X(6,11)=1.24999999999999E0_DP
 S_B(10)%A_Y(6,11)=-3.87499999999999E0_DP
 S_B(10)%A_X(6,12)=-1.00000000000001E0_DP
 S_B(10)%B_Y(6,12)=-0.833333333333331E0_DP
 S_B(10)%A_Y(6,13)=0.499999999999999E0_DP
 S_B(10)%A_X(6,14)=6.00000000000000E0_DP
 S_B(10)%B_Y(6,14)=5.00000000000000E0_DP
 S_B(10)%B_X(6,15)=-1.00000000000000E0_DP
 S_B(10)%A_Y(6,15)=1.00000000000000E0_DP
 S_B(10)%A_X(6,16)=-3.87499999999999E0_DP
 S_B(10)%B_Y(6,16)=-2.62499999999999E0_DP
 S_B(10)%B_X(6,17)=-1.04166666666666E0_DP
 S_B(10)%A_Y(6,17)=3.12500000000001E0_DP
 S_B(10)%A_X(6,18)=0.833333333333331E0_DP
 S_B(10)%B_Y(6,18)=0.833333333333337E0_DP
 S_B(10)%B_X(6,19)=12.5000000000000E0_DP
 S_B(10)%A_Y(6,19)=-12.5000000000000E0_DP
 S_B(10)%A_X(6,20)=5.00000000000000E0_DP
 S_B(10)%B_Y(6,20)=5.00000000000000E0_DP
 S_B(10)%B_X(6,22)=-1.75000000000000E0_DP
 S_B(10)%A_Y(6,22)=3.75000000000000E0_DP
 S_B(10)%A_X(6,23)=2.50000000000001E0_DP
 S_B(10)%B_Y(6,23)=2.00000000000000E0_DP
 S_B(10)%B_X(6,24)=0.833333333333337E0_DP
 S_B(10)%A_Y(6,24)=-2.50000000000000E0_DP
 S_B(10)%A_X(6,25)=-16.6666666666667E0_DP
 S_B(10)%B_Y(6,25)=-13.3333333333333E0_DP
 S_B(10)%B_X(6,26)=10.0000000000000E0_DP
 S_B(10)%A_Y(6,26)=-10.0000000000000E0_DP
 S_B(10)%A_X(6,29)=1.60714285714286E0_DP
 S_B(10)%B_Y(6,29)=1.03571428571428E0_DP
 S_B(10)%B_X(6,30)=0.999999999999999E0_DP
 S_B(10)%A_Y(6,30)=-2.00000000000000E0_DP
 S_B(10)%A_X(6,31)=-1.50000000000000E0_DP
 S_B(10)%B_Y(6,31)=-1.50000000000000E0_DP
 S_B(10)%B_X(6,32)=-10.0000000000000E0_DP
 S_B(10)%A_Y(6,32)=10.0000000000000E0_DP
 S_B(10)%A_X(6,33)=-10.0000000000000E0_DP
 S_B(10)%B_Y(6,33)=-10.0000000000000E0_DP
 S_B(10)%B_X(6,37)=0.258928571428571E0_DP
 S_B(10)%A_Y(6,37)=-0.491071428571428E0_DP
 S_B(10)%A_X(6,38)=-0.571428571428572E0_DP
 S_B(10)%B_Y(6,38)=-0.428571428571428E0_DP
 S_B(10)%B_X(6,39)=-0.500000000000001E0_DP
 S_B(10)%A_Y(6,39)=1.00000000000000E0_DP
 S_B(10)%A_X(6,40)=4.00000000000000E0_DP
 S_B(10)%B_Y(6,40)=3.00000000000000E0_DP
 S_B(10)%B_X(6,41)=-5.00000000000000E0_DP
 S_B(10)%A_Y(6,41)=5.00000000000000E0_DP
 S_B(10)%A_X(6,46)=-5.456349206349204E-002_DP
 S_B(10)%B_Y(6,46)=-3.472222222222217E-002_DP
 S_B(10)%B_X(6,47)=-5.357142857142854E-002_DP
 S_B(10)%A_Y(6,47)=8.928571428571441E-002_DP
 S_B(10)%A_X(6,48)=0.142857142857143E0_DP
 S_B(10)%B_Y(6,48)=0.142857142857143E0_DP
 S_B(10)%B_X(6,49)=0.500000000000000E0_DP
 S_B(10)%A_Y(6,49)=-0.500000000000000E0_DP
 S_B(10)%A_X(6,50)=1.00000000000000E0_DP
 S_B(10)%B_Y(6,50)=1.00000000000000E0_DP
 S_B(10)%A_Y(7,4)=-0.499999999999998E0_DP
 S_B(10)%B_X(7,6)=-1.00000000000000E0_DP
 S_B(10)%A_Y(7,6)=1.00000000000000E0_DP
 S_B(10)%A_X(7,7)=-1.16666666666666E0_DP
 S_B(10)%B_Y(7,7)=-0.999999999999999E0_DP
 S_B(10)%A_Y(7,8)=0.499999999999999E0_DP
 S_B(10)%A_X(7,9)=7.00000000000000E0_DP
 S_B(10)%B_Y(7,9)=6.00000000000000E0_DP
 S_B(10)%B_X(7,10)=-1.00000000000000E0_DP
 S_B(10)%A_Y(7,10)=1.00000000000000E0_DP
 S_B(10)%B_X(7,11)=-1.50000000000000E0_DP
 S_B(10)%A_Y(7,11)=4.49999999999999E0_DP
 S_B(10)%A_X(7,12)=0.999999999999998E0_DP
 S_B(10)%B_Y(7,12)=1.00000000000000E0_DP
 S_B(10)%B_X(7,13)=18.0000000000000E0_DP
 S_B(10)%A_Y(7,13)=-18.0000000000000E0_DP
 S_B(10)%A_X(7,14)=6.00000000000000E0_DP
 S_B(10)%B_Y(7,14)=6.00000000000000E0_DP
 S_B(10)%A_X(7,16)=4.49999999999999E0_DP
 S_B(10)%B_Y(7,16)=3.75000000000000E0_DP
 S_B(10)%B_X(7,17)=1.25000000000000E0_DP
 S_B(10)%A_Y(7,17)=-3.75000000000000E0_DP
 S_B(10)%A_X(7,18)=-30.0000000000000E0_DP
 S_B(10)%B_Y(7,18)=-25.0000000000000E0_DP
 S_B(10)%B_X(7,19)=15.0000000000000E0_DP
 S_B(10)%A_Y(7,19)=-15.0000000000000E0_DP
 S_B(10)%B_X(7,22)=2.50000000000000E0_DP
 S_B(10)%A_Y(7,22)=-4.99999999999999E0_DP
 S_B(10)%A_X(7,23)=-3.00000000000000E0_DP
 S_B(10)%B_Y(7,23)=-3.00000000000000E0_DP
 S_B(10)%B_X(7,24)=-25.0000000000000E0_DP
 S_B(10)%A_Y(7,24)=25.0000000000000E0_DP
 S_B(10)%A_X(7,25)=-20.0000000000000E0_DP
 S_B(10)%B_Y(7,25)=-20.0000000000000E0_DP
 S_B(10)%A_X(7,29)=-2.14285714285714E0_DP
 S_B(10)%B_Y(7,29)=-1.71428571428571E0_DP
 S_B(10)%B_X(7,30)=-1.50000000000000E0_DP
 S_B(10)%A_Y(7,30)=3.00000000000000E0_DP
 S_B(10)%A_X(7,31)=15.0000000000000E0_DP
 S_B(10)%B_Y(7,31)=12.0000000000000E0_DP
 S_B(10)%B_X(7,32)=-15.0000000000000E0_DP
 S_B(10)%A_Y(7,32)=15.0000000000000E0_DP
 S_B(10)%B_X(7,37)=-0.428571428571428E0_DP
 S_B(10)%A_Y(7,37)=0.714285714285714E0_DP
 S_B(10)%A_X(7,38)=0.857142857142857E0_DP
 S_B(10)%B_Y(7,38)=0.857142857142857E0_DP
 S_B(10)%B_X(7,39)=4.00000000000000E0_DP
 S_B(10)%A_Y(7,39)=-4.00000000000000E0_DP
 S_B(10)%A_X(7,40)=6.00000000000000E0_DP
 S_B(10)%B_Y(7,40)=6.00000000000000E0_DP
 S_B(10)%A_X(7,46)=7.936507936507931E-002_DP
 S_B(10)%B_Y(7,46)=5.952380952380952E-002_DP
 S_B(10)%B_X(7,47)=0.107142857142857E0_DP
 S_B(10)%A_Y(7,47)=-0.178571428571429E0_DP
 S_B(10)%A_X(7,48)=-0.571428571428571E0_DP
 S_B(10)%B_Y(7,48)=-0.428571428571428E0_DP
 S_B(10)%B_X(7,49)=1.00000000000000E0_DP
 S_B(10)%A_Y(7,49)=-1.00000000000000E0_DP
 S_B(10)%B_X(8,3)=-1.00000000000000E0_DP
 S_B(10)%A_Y(8,3)=1.00000000000000E0_DP
 S_B(10)%A_Y(8,4)=0.500000000000004E0_DP
 S_B(10)%A_X(8,5)=8.00000000000000E0_DP
 S_B(10)%B_Y(8,5)=7.00000000000000E0_DP
 S_B(10)%B_X(8,6)=-1.00000000000000E0_DP
 S_B(10)%A_Y(8,6)=1.00000000000000E0_DP
 S_B(10)%A_X(8,7)=1.16666666666667E0_DP
 S_B(10)%B_Y(8,7)=1.16666666666667E0_DP
 S_B(10)%B_X(8,8)=24.5000000000000E0_DP
 S_B(10)%A_Y(8,8)=-24.5000000000000E0_DP
 S_B(10)%A_X(8,9)=7.00000000000000E0_DP
 S_B(10)%B_Y(8,9)=7.00000000000000E0_DP
 S_B(10)%B_X(8,11)=1.75000000000001E0_DP
 S_B(10)%A_Y(8,11)=-5.25000000000001E0_DP
 S_B(10)%A_X(8,12)=-49.0000000000000E0_DP
 S_B(10)%B_Y(8,12)=-42.0000000000000E0_DP
 S_B(10)%B_X(8,13)=21.0000000000000E0_DP
 S_B(10)%A_Y(8,13)=-21.0000000000000E0_DP
 S_B(10)%A_X(8,16)=-5.25000000000001E0_DP
 S_B(10)%B_Y(8,16)=-5.25000000000001E0_DP
 S_B(10)%B_X(8,17)=-52.5000000000000E0_DP
 S_B(10)%A_Y(8,17)=52.5000000000000E0_DP
 S_B(10)%A_X(8,18)=-35.0000000000000E0_DP
 S_B(10)%B_Y(8,18)=-35.0000000000000E0_DP
 S_B(10)%B_X(8,22)=-3.50000000000001E0_DP
 S_B(10)%A_Y(8,22)=7.00000000000001E0_DP
 S_B(10)%A_X(8,23)=42.0000000000000E0_DP
 S_B(10)%B_Y(8,23)=35.0000000000000E0_DP
 S_B(10)%B_X(8,24)=-35.0000000000000E0_DP
 S_B(10)%A_Y(8,24)=35.0000000000000E0_DP
 S_B(10)%A_X(8,29)=3.00000000000000E0_DP
 S_B(10)%B_Y(8,29)=3.00000000000000E0_DP
 S_B(10)%B_X(8,30)=17.5000000000000E0_DP
 S_B(10)%A_Y(8,30)=-17.5000000000000E0_DP
 S_B(10)%A_X(8,31)=21.0000000000000E0_DP
 S_B(10)%B_Y(8,31)=21.0000000000000E0_DP
 S_B(10)%B_X(8,37)=0.750000000000001E0_DP
 S_B(10)%A_Y(8,37)=-1.25000000000000E0_DP
 S_B(10)%A_X(8,38)=-5.00000000000000E0_DP
 S_B(10)%B_Y(8,38)=-4.00000000000000E0_DP
 S_B(10)%B_X(8,39)=7.00000000000000E0_DP
 S_B(10)%A_Y(8,39)=-7.00000000000000E0_DP
 S_B(10)%A_X(8,46)=-0.138888888888889E0_DP
 S_B(10)%B_Y(8,46)=-0.138888888888889E0_DP
 S_B(10)%B_X(8,47)=-0.500000000000000E0_DP
 S_B(10)%A_Y(8,47)=0.500000000000000E0_DP
 S_B(10)%A_X(8,48)=-1.00000000000000E0_DP
 S_B(10)%B_Y(8,48)=-1.00000000000000E0_DP
 S_B(10)%B_X(9,1)=-1.00000000000000E0_DP
 S_B(10)%A_Y(9,1)=1.00000000000000E0_DP
 S_B(10)%A_X(9,2)=9.00000000000000E0_DP
 S_B(10)%B_Y(9,2)=8.00000000000000E0_DP
 S_B(10)%B_X(9,3)=-1.00000000000000E0_DP
 S_B(10)%A_Y(9,3)=1.00000000000000E0_DP
 S_B(10)%B_X(9,4)=32.0000000000000E0_DP
 S_B(10)%A_Y(9,4)=-32.0000000000000E0_DP
 S_B(10)%A_X(9,5)=8.00000000000000E0_DP
 S_B(10)%B_Y(9,5)=8.00000000000000E0_DP
 S_B(10)%A_X(9,7)=-74.6666666666667E0_DP
 S_B(10)%B_Y(9,7)=-65.3333333333333E0_DP
 S_B(10)%B_X(9,8)=28.0000000000000E0_DP
 S_B(10)%A_Y(9,8)=-28.0000000000000E0_DP
 S_B(10)%B_X(9,11)=-98.0000000000000E0_DP
 S_B(10)%A_Y(9,11)=98.0000000000000E0_DP
 S_B(10)%A_X(9,12)=-56.0000000000000E0_DP
 S_B(10)%B_Y(9,12)=-56.0000000000000E0_DP
 S_B(10)%A_X(9,16)=98.0000000000000E0_DP
 S_B(10)%B_Y(9,16)=84.0000000000000E0_DP
 S_B(10)%B_X(9,17)=-70.0000000000000E0_DP
 S_B(10)%A_Y(9,17)=70.0000000000000E0_DP
 S_B(10)%B_X(9,22)=56.0000000000000E0_DP
 S_B(10)%A_Y(9,22)=-56.0000000000000E0_DP
 S_B(10)%A_X(9,23)=56.0000000000000E0_DP
 S_B(10)%B_Y(9,23)=56.0000000000000E0_DP
 S_B(10)%A_X(9,29)=-24.0000000000000E0_DP
 S_B(10)%B_Y(9,29)=-20.0000000000000E0_DP
 S_B(10)%B_X(9,30)=28.0000000000000E0_DP
 S_B(10)%A_Y(9,30)=-28.0000000000000E0_DP
 S_B(10)%B_X(9,37)=-5.00000000000000E0_DP
 S_B(10)%A_Y(9,37)=5.00000000000000E0_DP
 S_B(10)%A_X(9,38)=-8.00000000000000E0_DP
 S_B(10)%B_Y(9,38)=-8.00000000000000E0_DP
 S_B(10)%A_X(9,46)=0.555555555555555E0_DP
 S_B(10)%B_Y(9,46)=0.444444444444444E0_DP
 S_B(10)%B_X(9,47)=-1.00000000000000E0_DP
 S_B(10)%A_Y(9,47)=1.00000000000000E0_DP
 S_B(10)%B_X(10,1)=-1.00000000000000E0_DP
 S_B(10)%A_Y(10,1)=1.00000000000000E0_DP
 S_B(10)%A_X(10,2)=9.00000000000000E0_DP
 S_B(10)%B_Y(10,2)=9.00000000000000E0_DP
 S_B(10)%B_X(10,4)=36.0000000000000E0_DP
 S_B(10)%A_Y(10,4)=-36.0000000000000E0_DP
 S_B(10)%A_X(10,7)=-84.0000000000000E0_DP
 S_B(10)%B_Y(10,7)=-84.0000000000000E0_DP
 S_B(10)%B_X(10,11)=-126.000000000000E0_DP
 S_B(10)%A_Y(10,11)=126.000000000000E0_DP
 S_B(10)%A_X(10,16)=126.000000000000E0_DP
 S_B(10)%B_Y(10,16)=126.000000000000E0_DP
 S_B(10)%B_X(10,22)=84.0000000000000E0_DP
 S_B(10)%A_Y(10,22)=-84.0000000000000E0_DP
 S_B(10)%A_X(10,29)=-36.0000000000000E0_DP
 S_B(10)%B_Y(10,29)=-36.0000000000000E0_DP
 S_B(10)%B_X(10,37)=-9.00000000000000E0_DP
 S_B(10)%A_Y(10,37)=9.00000000000000E0_DP
 S_B(10)%A_X(10,46)=1.00000000000000E0_DP
 S_B(10)%B_Y(10,46)=1.00000000000000E0_DP

end subroutine set_s_b  
 


!!!!!!!!!!!!!!!  New routines for solving Maxwell's Equations !!!!!!!!!!!!!!!

 SUBROUTINE  get_bend_coeff(s_b0t,NO1)
    implicit none
    integer no,n,i,k,j(2),NO1,l
    type(taylor) x,y,h,df,ker,sol
    type(complextaylor) z 
    type(taylor) f,kick_x,kick_y
    type(damap) y0
    real(dp) h0,cker,cker0
     TYPE(B_CYL), intent(inout) :: s_b0t
  
    no=NO1


    call init(no,1,0,0)


    call alloc(x,y,kick_x,kick_y)
    call alloc(z)
    call alloc(f,h,df,ker,sol)
    call alloc(y0)


      x=1.d0.mono.1
      y=1.d0.mono.2
      y0=1
      y0%v(2)=0
     z=x-i_*y   
    
!!!  



    
    h0=1.d0
    h=(1.d0+h0*x)
    do k=1,no
    !  erect multipole
    f=dreal(-z**K/K)
    
    df=f

    do i=k,no-1 !k+1
     df=h0*(df.d.1)/h
     call  invert_laplace(df)
      sol=f+df
      sol=-(sol.d.1)/h
      j=0
      j(1)=i
      cker=(sol.sub.j)
      df=df-cker*dreal(-z**(I+1)/(I+1))
      f=f+df
     enddo

  !     write(16,*) " Erect ",k
  !     call clean_taylor(f,f,1.d-6)
  !     call print(f,16)
  !   y=((f.d.1).d.1)+((f.d.2).d.2)-h0*(f.d.1)/h
  !  y=y.cut.(no-1)
!       call print(y,6)
!       write(6,*) full_abs(y)
       kick_x=(f.d.1)
       kick_y=(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)
        s_b0t%b_x(k,l)=kick_x.sub.j
        s_b0t%b_y(k,l)=kick_y.sub.j
       enddo

   !    bx=bx*y0
   !    call clean_taylor(bx,bx,1.d-6)
   !    call print(bx,6)
enddo    


    do k=1,no
    !  Skew multipole 
    f=aimag(-z**K/K)
    df=f
    do i=k,no-1 !k+1
     df=h0*(df.d.1)/h
     call  invert_laplace(df)
      sol=f+df
      sol=(sol.d.2)/h
      j=0
      j(1)=i
      cker=(sol.sub.j)
      df=df-cker*aimag(-z**(I+1)/(I+1))
      f=f+df
     enddo
  !     write(16,*) " Skew ",k
  !     call clean_taylor(f,f,1.d-6)
  !     call print(f,16)
  !   y=((f.d.1).d.1)+((f.d.2).d.2)-h0*(f.d.1)/h
  !  y=y.cut.(no-1)
 !      call print(y,6)
  !     write(6,*) full_abs(y)
       kick_x=(f.d.1)
       kick_y=(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)
        s_b0t%a_x(k,l)=kick_x.sub.j
        s_b0t%a_y(k,l)=kick_y.sub.j
       enddo

     !  write(6,*) k
     !  write(6,*) S_B0%A_X(k,7)
     !  write(6,*) S_B0%A_Y(k,8)    !   bx=bx*y0
    !   call clean_taylor(bx,bx,1.d-6)
    !   call print(bx,6)
enddo    


    call kill(x,y,kick_x,kick_y)
    call kill(z)
    call kill(f,h,df,ker,sol)
    call kill(y0)

    end subroutine get_bend_coeff


    subroutine invert_laplace(df)
    implicit none
    type(taylor), intent(inout) :: df
    type(taylorresonance) dfr
     call alloc(dfr)

     dfr=df
       dfr%cos=((dfr%cos.i.1).i.2)/4.0_dp
       dfr%sin=((dfr%sin.i.1).i.2)/4.0_dp
      df=dfr

     call kill(dfr)
    end subroutine invert_laplace

  subroutine make_coef(b,no)
    implicit none
    integer no
    integer i,j,k,a,m , ic

    type(B_CYL) b

    ic=1
    b%firsttime=-100
    allocate(b%nmul)
    allocate(b%n_mono)
    b%nmul=no
    b%n_mono=((no+2-ic)*(no+1-ic))/2
    allocate(b%i(b%n_mono),b%j(b%n_mono))
    allocate(b%a_x(no,b%n_mono),b%a_y(no,b%n_mono))
    allocate(b%b_x(no,b%n_mono),b%b_y(no,b%n_mono))

    do i=1,no
       do j=1,b%n_mono
          b%a_x(i,j)=0.0_dp
          b%a_y(i,j)=0.0_dp
          b%b_x(i,j)=0.0_dp
          b%b_y(i,j)=0.0_dp
       enddo
    enddo


    k=0
    m=no-1
    do a=m,1,-1
       do j=m-a,1,-1
          k=k+1
          b%i(k)=a
          b%j(k)=j
       enddo
       k=k+1
       b%i(k)=a
       b%j(k)=0
    enddo
    do j=m,1,-1
       k=k+1
       b%i(k)=0
       b%j(k)=j
    enddo
    k=k+1
    b%i(k)=0
    b%j(k)=0

  end subroutine make_coef

  subroutine nul_coef(b)
    implicit none
    type(B_CYL) b
    if(b%firsttime/=-100) then
       nullify(b%nmul)
       nullify(b%n_mono)
       nullify(b%i)
       nullify(b%j)
       nullify(b%a_x)
       nullify(b%a_y)
       nullify(b%b_x)
       nullify(b%b_y)
       b%firsttime=-100
    else
       deallocate(b%nmul)
       deallocate(b%n_mono)
       deallocate(b%i)
       deallocate(b%j)
       deallocate(b%a_x)
       deallocate(b%a_y)
       deallocate(b%b_x)
       deallocate(b%b_y)

    endif


  end subroutine nul_coef


end module S_status
