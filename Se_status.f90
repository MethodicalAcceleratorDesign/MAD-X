!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module S_status
  use s_frame
  USE S_extend_poly
  use anbn
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
  TYPE(B_CYL),ALLOCATABLE ::  S_B(:)
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
  real(dp), target :: muon=one
  LOGICAL(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.
  ! include "a_def_all_kind.inc"    ! sept 2007
  ! include "a_def_sagan.inc"
  ! include "a_def_element_fibre_layout.inc"
  !  !  include "a_def_user1.inc"
  !  !!  include "a_def_arbitrary.inc"
  !  !  include "a_def_user2.inc"

  TYPE (INTERNAL_STATE), PARAMETER :: DEFAULT0 = INTERNAL_STATE &
       &(0,f,f,f,f,f,f,f,f,f,f,f,3)
  TYPE (INTERNAL_STATE), PARAMETER :: TOTALPATH0 = INTERNAL_STATE &
       &(1,f,f,f,f,f,f,f,f,f,f,f,3)
  TYPE (INTERNAL_STATE), PARAMETER :: TIME0 = INTERNAL_STATE &
       &(0,t,f,f,f,f,f,f,f,f,f,f,3)
  TYPE (INTERNAL_STATE), PARAMETER :: RADIATION0 = INTERNAL_STATE &
       &(0,f,t,f,f,f,f,f,f,f,f,f,3)
  TYPE (INTERNAL_STATE), PARAMETER :: NOCAVITY0 = INTERNAL_STATE &
       &(0,f,f,t,f,f,f,f,f,f,f,f,3)
  TYPE (INTERNAL_STATE), PARAMETER :: FRINGE0 = INTERNAL_STATE &
       &(0,f,f,f,t,f,f,f,f,f,f,f,3)
  TYPE (INTERNAL_STATE), PARAMETER :: EXACTMIS0 = INTERNAL_STATE &
       &(0,f,f,f,f,t,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: STOCHASTIC0 = INTERNAL_STATE &
       &(0,f,f,f,f,f,t,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: ONLY_4d0 = INTERNAL_STATE &
       &(0,f,f,t,f,f,f,f,t,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: DELTA0   = INTERNAL_STATE &
       &(0,f,f,t,f,f,f,f,t,t,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: SPIN0   = INTERNAL_STATE &
       &(0,f,f,f,f,f,f,f,f,f,t,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: MODULATION0   = INTERNAL_STATE &
       &(0,f,f,f,f,f,f,f,f,f,f,t,3)
  TYPE(INTERNAL_STATE), target ::  DEFAULT=DEFAULT0
  TYPE(INTERNAL_STATE), target ::  TOTALPATH=TOTALPATH0
  TYPE(INTERNAL_STATE), target ::  RADIATION=RADIATION0
  TYPE(INTERNAL_STATE), target ::  NOCAVITY=NOCAVITY0
  TYPE(INTERNAL_STATE), target ::  FRINGE=FRINGE0
  TYPE(INTERNAL_STATE), target ::  TIME=TIME0
  TYPE(INTERNAL_STATE), target ::  EXACTMIS=EXACTMIS0
  TYPE(INTERNAL_STATE), target ::  STOCHASTIC=STOCHASTIC0
  TYPE(INTERNAL_STATE), target ::  ONLY_4D=ONLY_4D0
  TYPE(INTERNAL_STATE), target ::  DELTA=DELTA0
  TYPE(INTERNAL_STATE), target ::  SPIN=SPIN0
  TYPE(INTERNAL_STATE), target ::  MODULATION=MODULATION0

  !  private s_init,S_init_berz,MAKE_STATES_0,MAKE_STATES_m,print_s,CONV
  private s_init,MAKE_STATES_0,MAKE_STATES_m,print_s,CONV
  LOGICAL(lp), target :: stoch_in_rec = .false.
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
    P%KIND=0; P%R=ZERO;P%X=ZERO;P%Y=ZERO;
    ALLOCATE(P%DX);ALLOCATE(P%DY);
    P%DX=ZERO;P%DY=ZERO;
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
    P%LD=zero;P%B0=zero;P%LC=zero;
    ALLOCATE(P%TILTD);P%TILTD=zero;
    ! ALLOCATE(P%beta0);
    ! ALLOCATE(P%MASS0);
    ! ALLOCATE(P%gamma0I);
    ! ALLOCATE(P%gambet);
    ALLOCATE(P%P0C);
    ! P%beta0 =one;
    ! P%MASS0 =one;
    !P%gamma0I=zero;P%gambet =zero;
    P%P0C =zero;
    ALLOCATE(P%EDGE(2));P%EDGE(1)=zero;P%EDGE(2)=zero;
    !    ALLOCATE(P%TOTALPATH); ! PART OF A STATE INITIALIZED BY EL=DEFAULT
    ALLOCATE(P%EXACT);  !ALLOCATE(P%RADIATION);ALLOCATE(P%NOCAVITY);
    !    ALLOCATE(P%permFRINGE);
    ALLOCATE(P%KILL_ENT_FRINGE);ALLOCATE(P%KILL_EXI_FRINGE);ALLOCATE(P%bend_fringe); !ALLOCATE(P%TIME);
    ALLOCATE(P%METHOD);ALLOCATE(P%NST);P%METHOD=2;P%NST=1;
    ALLOCATE(P%NMUL);P%NMUL=0;
    !    ALLOCATE(P%spin);
    !   ALLOCATE(P%TRACK);P%TRACK=.TRUE.;
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
    !   elp%permFRINGE=el%permFRINGE
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
          IF((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2/E%R(2)**2>ONE) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=zero
             xlost=x
             messagelost="Lost in real kind=1 elliptic Aperture"
          ENDIF
       CASE(2)  ! rectangle
          IF(ABS(X(1)-E%DX)>E%X.OR.ABS(X(3)-E%DY)>E%Y) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=zero
             xlost=x
             messagelost="Lost in real kind=2 rectangular Aperture"
          ENDIF
       CASE(3)  ! RECTANGLE + ELLIPSE (CIRCLE)
          IF((ABS(X(1)-E%DX)>E%X).OR.(ABS(X(3)-E%DY)>E%Y).OR.  &
               ((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2**2/E%R(2)**2>ONE)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=zero
             xlost=x
             messagelost="Lost in real kind=3 rect-ellipse Aperture"
          ENDIF
       CASE(4) ! MARGUERITE
          IF(((X(1)-E%DX)**2/E%R(2)**2+(X(3)-E%DY)**2/E%R(1)**2>ONE).OR.  &
               ((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2/E%R(2)**2>ONE)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=zero
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
             xlost=zero
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

    MADFAC=one   ! to prevent overflow (David Sagan)
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
    !    CALL WRITE_I

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
    tilt%tilt(0)=zero
    do i=1,nmax
       tilt%tilt(i)=pih/i
    enddo
    !  SECTOR_B AND SECTOR_NMUL FOR TYPE TEAPOT

    change_sector=my_false

    IF(SECTOR_NMUL>0.and.firsttime_coef) THEN
       !  verb=global_verbose
       !  global_verbose=.false.
       if(firsttime_coef.or.(.not.allocated(S_B))) then
          !          SECTOR_B%firsttime=0   !slightly unsafe
          ALLOCATE(S_B(SECTOR_NMUL_MAX))
          lda_old=lda_used
          lda_used=3000
          DO I=1,SECTOR_NMUL_MAX
             !             if(i==SECTOR_NMUL_MAX)     global_verbose=.true.
             S_B(I)%firsttime=0
             call nul_coef(S_B(I))
             call make_coef(S_B(I),I)
             call curvebend(S_B(I),I)
          ENDDO
          lda_used=lda_old
       endif

       ! integer firsttime
       !  integer, POINTER ::  nmul,n_mono
       ! integer, DIMENSION(:), POINTER   :: i,j
       ! real(dp), DIMENSION(:,:), POINTER   :: a_x,a_y,b_x,b_y

       !       call nul_coef(SECTOR_B)
       !       call make_coef(SECTOR_B,SECTOR_NMUL)
       !       call curvebend(SECTOR_B,SECTOR_NMUL)
       !       w_p=1
       !       w_p%nc=1
       !       w_p%fc='(1((1X,A34)))'
       !       w_p%fI='(2((1X,I4)),/)'
       !       W_P=(/SECTOR_NMUL,sector_b%n_mono/)
       !       w_p%c(1) = " Small Machine Sector Bend Order ="
       !       CALL WRITE_I
       firsttime_coef=.FALSE.
    ENDIF

    call clear_states
    !  global_verbose=verb

  END  SUBROUTINE MAKE_STATES_0

  subroutine clear_states     !%nxyz
    implicit none
    DEFAULT=DEFAULT0
    TOTALPATH=TOTALPATH0
    RADIATION=RADIATION0
    NOCAVITY=NOCAVITY0
    FRINGE=FRINGE0
    TIME=TIME0
    EXACTMIS=EXACTMIS0
    exactmis%exactmis=exactmis%exactmis.and.(.not.sixtrack_compatible)
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
    if(CAVITY_TOTALPATH==1) write(mf,'((1X,a24))' ) ' Real Pill Box Cavities '
    if(CAVITY_TOTALPATH==0) write(mf,'((1X,a24))' ) ' Fake Pill Box Cavities '

    If(electron) then
       if(muon==one)  then
          write(mf,*)"This is an electron (positron actually if charge=1) "
       else
          write(mf,'((1X,a21,1x,G20.14,1x,A24))' ) "This a particle with ",muon, "times the electron mass "
       endif
    else
       write(mf,*) "This is a proton "
    endif
    write(mf, '((1X,a20,1x,a5))' )  "      EXACT_MODEL = ", CONV(EXACT_MODEL    )
    write(mf, '((1X,a20,1x,i4))' )  "      TOTALPATH   = ", S%TOTALPATH
    write(mf, '((1X,a20,1x,a5))' )  "      EXACTMIS    = ", CONV(S%EXACTMIS    )
    write(mf,'((1X,a20,1x,a5))' ) "      RADIATION   = ", CONV(S%RADIATION  )
    write(mf,'((1X,a20,1x,a5))' ) "      NOCAVITY    = ", CONV(S%NOCAVITY )
    write(mf,'((1X,a20,1x,a5))' ) "      TIME        = ", CONV(S%TIME )
    write(mf,'((1X,a20,1x,a5))' ) "      FRINGE      = ", CONV(S%FRINGE   )
    write(mf,'((1X,a20,1x,a5))' ) "      PARA_IN     = ", CONV(S%PARA_IN  )
    write(mf,'((1X,a20,1x,a5))' ) "      ONLY_4D     = ", CONV(S%ONLY_4D   )
    write(mf,'((1X,a20,1x,a5))' ) "      DELTA       = ", CONV(S%DELTA    )
    write(mf,'((1X,a20,1x,a5))' ) "      SPIN        = ", CONV(S%SPIN    )
    write(mf,'((1X,a20,1x,a5))' ) "      MODULATION   = ", CONV(S%MODULATION    )
    write(mf,'((1X,a20,1x,I4))' ) " SPIN DIMENSION   = ", S%SPIN_DIM
    !   CALL WRITE_I
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

    TIME=  TIME+DEFAULT

    FRINGE=  FRINGE+DEFAULT

    ONLY_4D= ONLY_4D+DEFAULT

    DELTA= DELTA+DEFAULT

    EXACTMIS= EXACTMIS+DEFAULT

    SPIN= SPIN+DEFAULT

    MODULATION= MODULATION+DEFAULT

  END  SUBROUTINE update_STATES


  SUBROUTINE  EQUALt(S2,S1)
    implicit none
    type (INTERNAL_STATE),INTENT(OUT)::S2
    type (INTERNAL_STATE),INTENT(IN)::S1

    S2%TOTALPATH=   S1%TOTALPATH
    S2%EXACTMIS=       S1%EXACTMIS
    S2%RADIATION=     S1%RADIATION
    S2%NOCAVITY=    S1%NOCAVITY
    S2%TIME=        S1%TIME
    S2%FRINGE=           S1%FRINGE
    S2%stochastic=           S1%stochastic
    S2%PARA_IN=     S1%PARA_IN
    S2%ONLY_4D=      S1%ONLY_4D
    S2%DELTA=       S1%DELTA
    S2%SPIN=       S1%SPIN
    S2%MODULATION=       S1%MODULATION
    S2%spin_dim=       S1%spin_dim
  END SUBROUTINE EQUALt



  FUNCTION add( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) add
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2


    add%TOTALPATH=0
    if((S1%TOTALPATH==1).OR.(S2%TOTALPATH==1)) add%TOTALPATH=1

    !    add%EXACT    =       S1%EXACT.OR.S2%EXACT
    add%EXACTMIS    =       (S1%EXACTMIS.OR.S2%EXACTMIS).and.(.not.sixtrack_compatible)
    add%RADIATION  =  S1%RADIATION.OR.S2%RADIATION
    add%NOCAVITY =  S1%NOCAVITY.OR.S2%NOCAVITY
    add%TIME     =  S1%TIME.OR.S2%TIME
    add%FRINGE   =       S1%FRINGE.OR.S2%FRINGE
    add%stochastic   =       S1%stochastic.OR.S2%stochastic
    add%ONLY_4D  =       S1%ONLY_4D.OR.S2%ONLY_4D
    add%DELTA  =       S1%DELTA.OR.S2%DELTA
    add%SPIN  =       S1%SPIN.OR.S2%SPIN
    add%MODULATION  =       S1%MODULATION.OR.S2%MODULATION
    add%PARA_IN  =       S1%PARA_IN.OR.S2%PARA_IN.or.ALWAYS_knobs
    add%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(add%stochastic) THEN
       add%RADIATION=T
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

    sub%EXACTMIS    =       (S1%EXACTMIS.min.S2%EXACTMIS).and.(.not.sixtrack_compatible)
    sub%RADIATION  =  S1%RADIATION.min.S2%RADIATION
    sub%NOCAVITY =  S1%NOCAVITY.min.S2%NOCAVITY
    sub%TIME     =  S1%TIME.min.S2%TIME
    sub%FRINGE   =       S1%FRINGE.min.S2%FRINGE
    sub%stochastic   =       S1%stochastic.min.S2%stochastic
    sub%ONLY_4D  =       S1%ONLY_4D.min.S2%ONLY_4D
    sub%DELTA  =       S1%DELTA.min.S2%DELTA
    sub%SPIN  =       S1%SPIN.min.S2%SPIN
    sub%MODULATION  = S1%MODULATION.min.S2%MODULATION
    sub%PARA_IN  =       (S1%PARA_IN.MIN.S2%PARA_IN).or.ALWAYS_knobs
    sub%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(sub%stochastic) THEN
       sub%RADIATION=T
    ENDIF
    IF(sub%DELTA) THEN
       sub%ONLY_4D=T
       sub%NOCAVITY =  T
    ENDIF
    IF(sub%ONLY_4D) THEN
       sub%TOTALPATH=  0
       sub%RADIATION  =  F
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
    E(1)=P%DIR*X(2)/(one+X5)
    E(2)=P%DIR*X(4)/(one+X5)
    E(3)=P%DIR*ROOT((one+X5)**2-X(2)**2-X(4)**2)/(one+X5)

    B2=zero
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
    E(1)=P%DIR*X(2)/(one+X5)
    E(2)=P%DIR*X(4)/(one+X5)
    E(3)=P%DIR*SQRT((one+X5)**2-X(2)**2-X(4)**2)/(one+X5)


    B2=zero
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


    IF(TILTD==zero) RETURN
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

    IF(TILTD==zero) RETURN
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


end module S_status
