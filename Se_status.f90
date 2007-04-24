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
  INTEGER,TARGET :: SECTOR_NMUL_MAX=10
  TYPE(B_CYL),ALLOCATABLE ::  S_B(:)
  INTEGER, target :: SECTOR_NMUL = 4
  INTEGER, TARGET :: NDPT_OTHER = 0
  real(dp) CRAD,CFLUC
  !  real(dp) YOSK(0:4), YOSD(4)    ! FIRST 6TH ORDER OF YOSHIDA
  !  real(dp),PARAMETER::AAA=-0.25992104989487316476721060727823e0_dp  ! fourth order integrator
  !  real(dp),PARAMETER::FD1=half/(one+AAA),FD2=AAA*FD1,FK1=one/(one+AAA),FK2=(AAA-one)*FK1

  LOGICAL(lp) :: firsttime_coef=.true.

  PRIVATE EQUALt,ADD,PARA_REMA,EQUALtilt
  !PRIVATE DTILTR,DTILTP,DTILTS
  PRIVATE DTILTR_EXTERNAL,DTILTP_EXTERNAL,DTILTS_EXTERNAL
  PRIVATE CHECK_APERTURE_R,CHECK_APERTURE_P,CHECK_APERTURE_S
  LOGICAL(lp), target:: electron
  real(dp), target :: muon=one
  LOGICAL(lp),PRIVATE,PARAMETER::T=.TRUE.,F=.FALSE.
  INTEGER,PARAMETER::NMAX=20
  include "a_def_all_kind.inc"
  include "a_def_sagan.inc"
  !  include "a_def_user1.inc"
  !!  include "a_def_arbitrary.inc"
  !  include "a_def_user2.inc"
  include "a_def_element_fibre_layout.inc"
  TYPE(INTERNAL_STATE), target ::  DEFAULT
  TYPE(INTERNAL_STATE), target ::  TOTALPATH,RADIATION,NOCAVITY,FRINGE,TIME,EXACTMIS
  TYPE(INTERNAL_STATE), target ::  ONLY_4D,DELTA,SPIN,SPIN_ONLY

  TYPE (INTERNAL_STATE), PARAMETER :: DEFAULT0 = INTERNAL_STATE &
       &(f,f,f,f,f,f,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: TOTALPATH0 = INTERNAL_STATE &
       &(t,f,f,f,f,f,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: TIME0 = INTERNAL_STATE &
       &(f,t,f,f,f,f,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: RADIATION0 = INTERNAL_STATE &
       &(f,f,t,f,f,f,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: NOCAVITY0 = INTERNAL_STATE &
       &(f,f,f,t,f,f,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: FRINGE0 = INTERNAL_STATE &
       &(f,f,f,f,t,f,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: EXACTMIS0 = INTERNAL_STATE &
       &(f,f,f,f,f,t,f,f,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: ONLY_4d0 = INTERNAL_STATE &
       &(f,f,f,t,f,f,f,t,f,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: DELTA0   = INTERNAL_STATE &
       &(f,f,f,t,f,f,f,t,t,f,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: SPIN0   = INTERNAL_STATE &
       &(f,f,f,f,f,f,f,f,f,t,f,F,3)
  TYPE (INTERNAL_STATE), PARAMETER :: SPIN_ONLY0   = INTERNAL_STATE &
       &(f,f,f,f,f,f,f,f,f,f,t,F,3)
  private s_init,S_init_berz,MAKE_STATES_0,MAKE_STATES_m,print_s,CONV
  LOGICAL(lp), target :: stoch_in_rec = .false.
  private alloc_p,equal_p,dealloc_p,alloc_A,equal_A,dealloc_A !,NULL_p
  PRIVATE B2PERPR,B2PERPP !,S_init_berz0
  type(tilting) tilt
  private minu
  real(dp) MADFAC(NMAX)
  CHARACTER(24) MYTYPE(-100:100)

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
     MODULE PROCEDURE CHECK_APERTURE_S
  END INTERFACE

  INTERFACE init
     MODULE PROCEDURE s_init
     !     MODULE PROCEDURE S_init_berz0
     MODULE PROCEDURE S_init_berz
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE print_s
  END INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE alloc_p
     MODULE PROCEDURE alloc_A
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE dealloc_p
     MODULE PROCEDURE dealloc_A
  END INTERFACE

  INTERFACE B2PERP
     MODULE PROCEDURE B2PERPR
     MODULE PROCEDURE B2PERPP
  END INTERFACE

  !  INTERFACE DTILTD
  !     MODULE PROCEDURE DTILTR
  !     MODULE PROCEDURE DTILTP       ! DESIGN TILT
  !     MODULE PROCEDURE DTILTS
  !     MODULE PROCEDURE DTILTR_EXTERNAL
  !     MODULE PROCEDURE DTILTP_EXTERNAL       ! EXTERNAL
  !     MODULE PROCEDURE DTILTS_EXTERNAL
  !  END INTERFACE

  INTERFACE DTILTD
     MODULE PROCEDURE DTILTR_EXTERNAL
     MODULE PROCEDURE DTILTP_EXTERNAL       ! EXTERNAL
     MODULE PROCEDURE DTILTS_EXTERNAL
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

    nullify(P%KIND);nullify(P%R);nullify(P%X);nullify(P%Y);
  end subroutine NULL_A

  SUBROUTINE  alloc_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    nullify(p)
    allocate(p)
    CALL NULL_A(p)
    ALLOCATE(P%R(2));ALLOCATE(P%X);ALLOCATE(P%Y);ALLOCATE(P%KIND);
    P%KIND=0; P%R=ZERO;P%X=ZERO;P%Y=ZERO;

  end subroutine alloc_A

  SUBROUTINE  dealloc_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    if(associated(p%R)) then
       DEALLOCATE(P%R);DEALLOCATE(P%X);DEALLOCATE(P%Y);DEALLOCATE(P%KIND);
    endif
  end SUBROUTINE  dealloc_A



  SUBROUTINE  NULL_p(p)
    implicit none
    type (MAGNET_CHART), pointer:: P

    nullify(P%LD);nullify(P%B0);nullify(P%LC);
    nullify(P%TILTD);  nullify(P%dir);nullify(P%charge);
    nullify(P%beta0);nullify(P%gamma0I);nullify(P%gambet);nullify(P%P0C);
    nullify(P%EDGE)
    nullify(P%TOTALPATH)
    nullify(P%EXACT);nullify(P%RADIATION);nullify(P%RADIATION_NEW);nullify(P%NOCAVITY);
    nullify(P%FRINGE);nullify(P%KILL_ENT_FRINGE);nullify(P%KILL_EXI_FRINGE);nullify(P%bend_fringe);nullify(P%TIME);
    nullify(P%METHOD);nullify(P%NST);
    nullify(P%NMUL);nullify(P%spin);
    nullify(P%F);
    nullify(P%APERTURE);
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
    ALLOCATE(P%beta0);ALLOCATE(P%gamma0I);ALLOCATE(P%gambet);ALLOCATE(P%P0C);
    P%beta0 =one;P%gamma0I=zero;P%gambet =zero;P%P0C =zero;
    ALLOCATE(P%EDGE(2));P%EDGE(1)=zero;P%EDGE(2)=zero;
    ALLOCATE(P%TOTALPATH); ! PART OF A STATE INITIALIZED BY EL=DEFAULT
    ALLOCATE(P%EXACT);ALLOCATE(P%RADIATION);ALLOCATE(P%RADIATION_NEW);ALLOCATE(P%NOCAVITY);
    ALLOCATE(P%FRINGE);ALLOCATE(P%KILL_ENT_FRINGE);ALLOCATE(P%KILL_EXI_FRINGE);ALLOCATE(P%bend_fringe);ALLOCATE(P%TIME);
    ALLOCATE(P%METHOD);ALLOCATE(P%NST);P%METHOD=2;P%NST=1;
    ALLOCATE(P%NMUL);P%NMUL=0;
    ALLOCATE(P%spin);
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

    !    if(associated(P%dir)) then
    !    endif
    DEALLOCATE(P%LD);DEALLOCATE(P%B0);DEALLOCATE(P%LC);
    DEALLOCATE(P%TILTD);
    DEALLOCATE(P%beta0);DEALLOCATE(P%gamma0I);DEALLOCATE(P%gambet);DEALLOCATE(P%P0C);
    if(associated(p%f)) then
       call kill(p%f)
       DEALLOCATE(P%f);
    endif
    if(associated(p%APERTURE)) then
       CALL kill(p%APERTURE)
       DEALLOCATE(p%APERTURE);
    endif
    DEALLOCATE(P%EDGE);
    DEALLOCATE(P%TOTALPATH);
    DEALLOCATE(P%EXACT);DEALLOCATE(P%RADIATION);DEALLOCATE(P%RADIATION_NEW);DEALLOCATE(P%NOCAVITY);
    DEALLOCATE(P%FRINGE);DEALLOCATE(P%KILL_ENT_FRINGE);DEALLOCATE(P%KILL_EXI_FRINGE);DEALLOCATE(P%bend_fringe);DEALLOCATE(P%TIME);
    DEALLOCATE(P%METHOD);DEALLOCATE(P%spin);DEALLOCATE(P%NST);
    DEALLOCATE(P%NMUL)
    !    CALL NULL_P(P)
    DEALLOCATE(P)
    nullify(p);
    ! if(junk) ccc=ccc-1
  end subroutine dealloc_p

  SUBROUTINE  equal_A(elp,el)
    implicit none
    type (MADX_APERTURE),INTENT(inOUT)::elP
    type (MADX_APERTURE),INTENT(IN)::el
    ELP%KIND=EL%KIND
    ELP%R=EL%R
    ELP%X=EL%X
    ELP%Y=EL%Y
  END SUBROUTINE  equal_A



  SUBROUTINE  equal_p(elp,el)
    implicit none
    type (MAGNET_CHART),INTENT(inOUT)::elP
    type (MAGNET_CHART),INTENT(IN)::el

    elp%beta0 =el%beta0
    elp%gamma0I=el%gamma0I
    elp%gambet =el%gambet
    elp%P0C =el%P0C
    elp%EXACT=el%EXACT
    elp%RADIATION=el%RADIATION
    elp%RADIATION_NEW=el%RADIATION_NEW
    elp%TIME=el%TIME
    elp%NOCAVITY=el%NOCAVITY
    elp%spin=el%spin
    elp%FRINGE=el%FRINGE
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
    ELP%TOTALPATH=EL%TOTALPATH
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



    IF(CHECK_MADX_APERTURE.AND.APERTURE_FLAG) THEN

       SELECT CASE(E%KIND)
       CASE(1)  ! ellipse circles
          IF(X(1)**2/E%R(1)**2+X(3)**2/E%R(2)**2>ONE) THEN
             CHECK_STABLE=.FALSE.
             CHECK_MADX_APERTURE=.false.
          ENDIF
       CASE(2)  ! rectangle
          IF(ABS(X(1))>E%X.OR.ABS(X(3))>E%Y) THEN
             CHECK_STABLE=.FALSE.
             CHECK_MADX_APERTURE=.false.
          ENDIF
       CASE(3)  ! RECTANGLE + ELLIPSE (CIRCLE)
          IF((ABS(X(1))>E%X).OR.(ABS(X(3))>E%Y).OR.(X(1)**2/E%R(1)**2+X(3)**2/E%R(2)**2>ONE)) THEN
             CHECK_STABLE=.FALSE.
             CHECK_MADX_APERTURE=.false.
          ENDIF
       CASE(4) ! MARGUERITE
          IF((X(1)**2/E%R(2)**2+X(3)**2/E%R(1)**2>ONE).OR.  &
               (X(1)**2/E%R(1)**2+X(3)**2/E%R(2)**2>ONE)) THEN
             CHECK_STABLE=.FALSE.
             CHECK_MADX_APERTURE=.false.
          ENDIF
       CASE(5) ! PILES OF POINTS
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


  SUBROUTINE  CHECK_APERTURE_S(E,X)
    implicit none
    type (MADX_APERTURE),INTENT(IN)::E
    TYPE(ENV_8), INTENT(IN):: X(6)
    REAL(DP) Y(6)
    Y=X
    CALL CHECK_APERTURE(E,Y)

  END SUBROUTINE  CHECK_APERTURE_S

  FUNCTION minu( S1,S2  )
    implicit none
    logical(lp) minu
    logical(lp), INTENT (IN) :: S1
    logical(lp), INTENT (IN) :: S2

    minu=.false.
    if(s1.and.(.not.s2)) minu=.true.

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
    LOGICAL(lp) particle,verb
    integer i,MF,lda_old

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
    ONLY_4D=ONLY_4d0
    DELTA=DELTA0
    SPIN=SPIN0
    SPIN_ONLY=SPIN_ONLY0
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

    If(electron) then
       if(muon==one)  then
          write(mf,*)"This is an electron (positron actually) "
       else
          write(mf,'((1X,a21,1x,G20.14,1x,A24))' ) "This a particle with ",muon, "times the electron mass "
       endif
    else
       write(mf,*) "This is a proton "
    endif
    write(mf, '((1X,a20,1x,a5))' )  "      EXACT_MODEL = ", CONV(EXACT_MODEL    )
    write(mf, '((1X,a20,1x,a5))' )  "      TOTALPATH   = ", CONV(S%TOTALPATH)
    write(mf, '((1X,a20,1x,a5))' )  "      EXACTMIS    = ", CONV(S%EXACTMIS    )
    write(mf,'((1X,a20,1x,a5))' ) "      RADIATION   = ", CONV(S%RADIATION  )
    write(mf,'((1X,a20,1x,a5))' ) "      NOCAVITY    = ", CONV(S%NOCAVITY )
    write(mf,'((1X,a20,1x,a5))' ) "      TIME        = ", CONV(S%TIME )
    write(mf,'((1X,a20,1x,a5))' ) "      FRINGE      = ", CONV(S%FRINGE   )
    write(mf,'((1X,a20,1x,a5))' ) "      PARA_IN     = ", CONV(S%PARA_IN  )
    write(mf,'((1X,a20,1x,a5))' ) "      ONLY_4D     = ", CONV(S%ONLY_4D   )
    write(mf,'((1X,a20,1x,a5))' ) "      DELTA       = ", CONV(S%DELTA    )
    write(mf,'((1X,a20,1x,a5))' ) "      SPIN        = ", CONV(S%SPIN    )
    write(mf,'((1X,a20,1x,a5))' ) "      SPIN_ONLY   = ", CONV(S%SPIN_ONLY    )
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
    real(dp) muonfactor
    CALL MAKE_YOSHIDA
    muon=muonfactor
    call MAKE_STATES_0(doneitt)
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

    SPIN_ONLY= SPIN_ONLY+DEFAULT

  END  SUBROUTINE update_STATES


  SUBROUTINE  EQUALt(S2,S1)
    implicit none
    type (INTERNAL_STATE),INTENT(OUT)::S2
    type (INTERNAL_STATE),INTENT(IN)::S1

    S2%TOTALPATH=   S1%TOTALPATH
    S2%EXACTMIS=       S1%EXACTMIS
    S2%RADIATION=     S1%RADIATION
    S2%RADIATION_NEW=     S1%RADIATION_NEW
    S2%NOCAVITY=    S1%NOCAVITY
    S2%TIME=        S1%TIME
    S2%FRINGE=           S1%FRINGE
    S2%PARA_IN=     S1%PARA_IN
    S2%ONLY_4D=      S1%ONLY_4D
    S2%DELTA=       S1%DELTA
    S2%SPIN=       S1%SPIN
    S2%SPIN_ONLY=       S1%SPIN_ONLY
    S2%spin_dim=       S1%spin_dim
  END SUBROUTINE EQUALt



  FUNCTION add( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) add
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2


    add%TOTALPATH=  S1%TOTALPATH.OR.S2%TOTALPATH
    !    add%EXACT    =       S1%EXACT.OR.S2%EXACT
    add%EXACTMIS    =       S1%EXACTMIS.OR.S2%EXACTMIS
    add%RADIATION  =  S1%RADIATION.OR.S2%RADIATION
    add%NOCAVITY =  S1%NOCAVITY.OR.S2%NOCAVITY
    add%TIME     =  S1%TIME.OR.S2%TIME
    add%FRINGE   =       S1%FRINGE.OR.S2%FRINGE
    add%ONLY_4D  =       S1%ONLY_4D.OR.S2%ONLY_4D
    add%DELTA  =       S1%DELTA.OR.S2%DELTA
    add%SPIN  =       S1%SPIN.OR.S2%SPIN
    add%SPIN_ONLY  =       S1%SPIN_ONLY.OR.S2%SPIN_ONLY
    add%PARA_IN  =       S1%PARA_IN.OR.S2%PARA_IN.or.ALWAYS_knobs
    add%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(add%DELTA) THEN
       add%ONLY_4D=T
       add%NOCAVITY =  T
    ENDIF
    IF(add%ONLY_4D) THEN
       add%TOTALPATH=  F
       add%RADIATION  =  F
       add%NOCAVITY =  T
    ENDIF
  END FUNCTION add

  FUNCTION sub( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) sub
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2


    sub%TOTALPATH=  S1%TOTALPATH.min.S2%TOTALPATH
    sub%EXACTMIS    =       S1%EXACTMIS.min.S2%EXACTMIS
    sub%RADIATION  =  S1%RADIATION.min.S2%RADIATION
    sub%NOCAVITY =  S1%NOCAVITY.min.S2%NOCAVITY
    sub%TIME     =  S1%TIME.min.S2%TIME
    sub%FRINGE   =       S1%FRINGE.min.S2%FRINGE
    sub%ONLY_4D  =       S1%ONLY_4D.min.S2%ONLY_4D
    sub%DELTA  =       S1%DELTA.min.S2%DELTA
    sub%SPIN  =       S1%SPIN.min.S2%SPIN
    sub%SPIN_ONLY  = S1%SPIN_ONLY.min.S2%SPIN_ONLY
    sub%PARA_IN  =       (S1%PARA_IN.MIN.S2%PARA_IN).or.ALWAYS_knobs
    sub%SPIN_DIM  =       MAX(S1%SPIN_DIM,S2%SPIN_DIM)
    IF(sub%DELTA) THEN
       sub%ONLY_4D=T
       sub%NOCAVITY =  T
    ENDIF
    IF(sub%ONLY_4D) THEN
       sub%TOTALPATH=  F
       sub%RADIATION  =  F
       sub%NOCAVITY =  T
    ENDIF

  END FUNCTION sub

  FUNCTION PARA_REMA(S1)   ! UNARY +
    implicit none
    TYPE (INTERNAL_STATE) PARA_REMA
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1

    PARA_REMA         =    S1
    PARA_REMA%PARA_IN =     T

  END FUNCTION PARA_REMA



  subroutine S_init(STATE,NO1,NP1,PACKAGE,ND2,NPARA)
    !  subroutine S_init(STATE,NO1,NP1,PACKAGE,MAPINT,ND2,NPARA)
    implicit none
    TYPE (INTERNAL_STATE), INTENT(IN):: STATE
    LOGICAL(lp), INTENT(IN):: PACKAGE
    INTEGER, INTENT(IN):: NO1,NP1
    INTEGER ND1,NDEL,NDPT1
    INTEGER,optional :: ND2,NPARA
    INTEGER  ND2l,NPARAl,NSPIN1

    IF(STATE%SPIN_ONLY) THEN

       NDEL=0
       NDPT1=0
       ND1=0
       IF(STATE%NOCAVITY)  THEN
          IF(STATE%ONLY_4D) THEN
             IF(STATE%DELTA) THEN
                NDEL=1
                !             MAPINT=5
             ELSE
                NDEL=0
                !            MAPINT=4
             ENDIF
          ELSE
             NDEL=0
             !         MAPINT=6
          ENDIF
       ENDIF
       NSPIN1=STATE%SPIN_DIM
       CALL init_SPIN(NO1,ND1,NP1+NDEL,NSPIN1,NDPT1,PACKAGE)


       ND2l=ND1*2
       NPARAl=ND2l+NDEL+C_%NSPIN
       C_%NPARA=NPARAl
       C_%ND2=ND2l
       C_%npara_fpp=NPARAl
       C_%SPIN_POS=C_%NPARA-C_%NSPIN+1

       if(present(nd2)) nd2=nd2l
       if(present(npara)) npara=nparal
    ELSEIF(STATE%SPIN.AND.(.NOT.STATE%SPIN_ONLY)) THEN
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
             NDPT1=5+C_%NDPT_OTHER
             !         MAPINT=6
          ENDIF
       ELSE              ! CAVITY IN RING
          ND1=3
          NDPT1=0
          !       MAPINT=6
       ENDIF

       NSPIN1=STATE%SPIN_DIM
       CALL init_SPIN(NO1,ND1,NP1+NDEL,NSPIN1,NDPT1,PACKAGE)

       ND2l=ND1*2
       NPARAl=ND2l+NDEL+C_%NSPIN
       C_%NPARA=NPARAl
       C_%ND2=ND2l
       C_%npara_fpp=NPARAl
       C_%SPIN_POS=C_%NPARA-C_%NSPIN+1

       if(present(nd2)) nd2=nd2l
       if(present(npara)) npara=nparal
    ELSE
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
             NDPT1=5+C_%NDPT_OTHER
             !         MAPINT=6
          ENDIF
       ELSE              ! CAVITY IN RING
          ND1=3
          NDPT1=0
          !       MAPINT=6
       ENDIF

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
    ENDIF
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

  subroutine S_init_berz(STATE,NO1,NP1,ND2,NPARA)
    implicit none
    TYPE (INTERNAL_STATE), INTENT(IN):: STATE
    INTEGER, INTENT(IN):: NO1,NP1
    INTEGER, INTENT(OUT)::    ND2,NPARA

    call init(STATE,NO1,NP1,my_true,ND2,NPARA)
    C_%NPARA=NPARA
    C_%npara_fpp=NPARA
  END  subroutine S_init_berz


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

  SUBROUTINE DTILTS_EXTERNAL(DIR,TILTD,I,Y)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    INTEGER,INTENT(IN):: I,DIR
    REAL(DP),INTENT(IN) :: TILTD
    TYPE(REAL_8) X(6)

    IF(TILTD==zero) RETURN

    CALL ALLOC(X)
    X=Y
    CALL DTILTD(DIR,TILTD,I,X)
    Y=X
    CALL KILL(X)

  END SUBROUTINE DTILTS_EXTERNAL

end module S_status
