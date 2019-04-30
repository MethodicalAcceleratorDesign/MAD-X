!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module S_status
  use s_frame
  USE S_extend_poly
    use c_TPSA
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
  integer, parameter :: KINDhel = KIND22
  integer, parameter :: KINDwiggler = KIND23+2
  !  integer, parameter :: KINDmu      = KIND23+3
  integer, parameter :: KINDpa     = KIND23+3
  integer, parameter :: kindsuperdrift = KIND23+4
  integer, parameter :: kindabell= KIND23+5
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

  LOGICAL(lp) :: firsttime_coef=.true.,read_sector_info=my_true

  PRIVATE EQUALt,ADD,PARA_REMA,EQUALtilt,EQUALi
  !PRIVATE DTILTR,DTILTP,DTILTS
  PRIVATE DTILTR_EXTERNAL,DTILTP_EXTERNAL,orthonormaliser,orthonormalisep
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
  TYPE(INTERNAL_STATE), target ::  only_2d=only_2d0
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
  logical(lp) :: automatic_complex = my_true
  integer :: aperture_pos_default=0
  private track_TREE_G_complexr,track_TREE_G_complexp,track_TREE_probe_complexr,track_TREE_probe_complexp_new
  integer :: size_tree=15
  integer :: ind_spin(3,3),k1_spin(9),k2_spin(9)
  real(dp),TARGET ::INITIAL_CHARGE=1
  logical :: mcmillan=.false.
  real(dp) :: radfac=1   ! to fudge radiation (lower it)
  TYPE B_CYL
     integer firsttime
     integer, POINTER ::  nmul,n_mono   !,nmul_e,n_mono_e
     integer, DIMENSION(:), POINTER   :: i,j   !,ie,je
     real(dp), DIMENSION(:,:), POINTER   :: a_x,a_y,b_x,b_y,va,vb
  END  TYPE B_CYL

 ! TYPE(B_CYL),ALLOCATABLE ::  S_B(:) 
 ! TYPE(B_CYL) S_B,S_EB
  TYPE(B_CYL) S_E,S_B_from_V

  INTERFACE OPERATOR (.min.)
     MODULE PROCEDURE minu                       ! to define the minus of Schmidt
  END INTERFACE

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUALt
     MODULE PROCEDURE EQUALi
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

  INTERFACE orthonormalise
     MODULE PROCEDURE orthonormaliser
     MODULE PROCEDURE orthonormalisep
  END INTERFACE

  INTERFACE DTILTD
     MODULE PROCEDURE DTILTR_EXTERNAL
     MODULE PROCEDURE DTILTP_EXTERNAL       ! EXTERNAL
  END INTERFACE

  INTERFACE track_TREE_G_complex
     MODULE PROCEDURE track_TREE_G_complexr
     MODULE PROCEDURE track_TREE_G_complexp       ! EXTERNAL
  END INTERFACE

  INTERFACE track_TREE_probe_complex
     MODULE PROCEDURE track_TREE_probe_complexr
     MODULE PROCEDURE track_TREE_probe_complexp_new
  END INTERFACE 

CONTAINS

  real(dp) function cradf(p)
    implicit none
    type (MAGNET_CHART), pointer:: P
    cradf=radfac*crad*p%p0c**3
  end function cradf

  real(dp) function cflucf(p)
    implicit none
    type (MAGNET_CHART), pointer:: P
    cflucf=cfluc*p%p0c**5
  end function cflucf


  SUBROUTINE  NULL_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    nullify(P%KIND);nullify(P%R);nullify(P%X);nullify(P%Y);nullify(P%dX);nullify(P%dY);nullify(P%pos);
  end subroutine NULL_A

  SUBROUTINE  alloc_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    nullify(p)
    allocate(p)
    CALL NULL_A(p)
    ALLOCATE(P%R(2));ALLOCATE(P%X);ALLOCATE(P%Y);ALLOCATE(P%KIND);ALLOCATE(P%pos);
    P%KIND=0; P%R=0.0_dp;P%X=0.0_dp;P%Y=0.0_dp;P%pos=aperture_pos_default;
    ALLOCATE(P%DX);ALLOCATE(P%DY);
    P%DX=0.0_dp;P%DY=0.0_dp;
  end subroutine alloc_A

  SUBROUTINE  dealloc_A(p)
    implicit none
    type (MADX_APERTURE), pointer:: P

    if(associated(p%R)) then
       DEALLOCATE(P%R);DEALLOCATE(P%X);DEALLOCATE(P%Y);DEALLOCATE(P%KIND);
       DEALLOCATE(P%DX);DEALLOCATE(P%DY);DEALLOCATE(P%pos);
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
    nullify(P%permFRINGE,p%highest_fringe);
    nullify(P%KILL_ENT_FRINGE);nullify(P%KILL_EXI_FRINGE);nullify(P%bend_fringe);  !nullify(P%TIME);
    nullify(P%KILL_ENT_SPIN);nullify(P%KILL_EXI_SPIN);
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
        ALLOCATE(P%permFRINGE,p%highest_fringe);
    ALLOCATE(P%KILL_ENT_FRINGE);ALLOCATE(P%KILL_EXI_FRINGE);ALLOCATE(P%bend_fringe); !ALLOCATE(P%TIME);
    ALLOCATE(P%KILL_ENT_SPIN);ALLOCATE(P%KILL_EXI_SPIN);
    ALLOCATE(P%METHOD);ALLOCATE(P%NST);P%METHOD=2;P%NST=1;
    ALLOCATE(P%NMUL);P%NMUL=0;
    !    ALLOCATE(P%spin);
    !   ALLOCATE(P%TRACK);P%TRACK=.TRUE.;
    P%permFRINGE=0
    p%highest_fringe=highest_fringe
    P%KILL_ENT_FRINGE=.FALSE.
    P%KILL_EXI_FRINGE=.FALSE.
    P%KILL_ENT_SPIN=.FALSE.
    P%KILL_EXI_SPIN=.FALSE.
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
    if(associated(p%highest_fringe))DEALLOCATE(p%highest_fringe);
    if(associated(p%KILL_ENT_FRINGE))DEALLOCATE(P%KILL_ENT_FRINGE);
    if(associated(p%KILL_EXI_FRINGE))DEALLOCATE(P%KILL_EXI_FRINGE);
    if(associated(p%KILL_ENT_SPIN))DEALLOCATE(p%KILL_ENT_SPIN);
    if(associated(p%KILL_EXI_SPIN))DEALLOCATE(p%KILL_EXI_SPIN);
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
    ELP%pos=EL%pos
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
       elp%highest_fringe=el%highest_fringe

    elp%KILL_ENT_FRINGE=el%KILL_ENT_FRINGE
    elp%KILL_EXI_FRINGE=el%KILL_EXI_FRINGE
    elp%KILL_ENT_SPIN=el%KILL_ENT_SPIN
    elp%KILL_EXI_SPIN=el%KILL_EXI_SPIN
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
             !messagelost="Lost in real kind=1 elliptic Aperture"
             write(messagelost,*) "Se_status.f90 CHECK_APERTURE_R : Lost in real kind=1 elliptic aperture. ", &
                                   "Orbit: X=",X(1)," Y=",X(3)," Ap.: DX=",E%DX," DY=",E%DY," R1=",E%R(1)," R2=",E%R(2)
          ENDIF
       CASE(2)  ! rectangle
          IF(ABS(X(1)-E%DX)>E%X.OR.ABS(X(3)-E%DY)>E%Y) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             !messagelost="Lost in real kind=2 rectangular Aperture"
             write(messagelost,*) "Se_status.f90 CHECK_APERTURE_R : Lost in real kind=2 rectangular aperture. ", &
                                  "Orbit: X=",X(1)," Y=",X(3)," Ap.: DX=",E%DX," DY=",E%DY," X=",E%X," Y=",E%Y
          ENDIF
       CASE(3)  ! RECTANGLE + ELLIPSE (CIRCLE)
          IF((ABS(X(1)-E%DX)>E%X).OR.(ABS(X(3)-E%DY)>E%Y).OR.  &
               ((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2**2/E%R(2)**2>1.0_dp)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             !messagelost="Lost in real kind=3 rect-ellipse Aperture"
             write(messagelost,*) "Se_status.f90 CHECK_APERTURE_R : Lost in real kind=3 rect-ellipse aperture. ", &
                                  "Orbit: X=",X(1)," Y=",X(3)," Ap.: DX=",E%DX," DY=",E%DY," X=",E%X," Y=",E%Y," R=",E%R
          ENDIF
       CASE(4) ! MARGUERITE
          IF(((X(1)-E%DX)**2/E%R(2)**2+(X(3)-E%DY)**2/E%R(1)**2>1.0_dp).OR.  &
               ((X(1)-E%DX)**2/E%R(1)**2+(X(3)-E%DY)**2/E%R(2)**2>1.0_dp)) THEN
             CHECK_STABLE=.FALSE.
             STABLE_DA=.false.
             xlost=0.0_dp
             xlost=x
             !messagelost="Lost in real kind=4 marguerite Aperture"
             write(messagelost,*) "Se_status.f90 CHECK_APERTURE_R : Lost in real kind=4 marguerite Aperture. ", &
                                  "Orbit: X=",X(1)," Y=",X(3)," Ap.: DX=",E%DX," DY=",E%DY," X=",E%X," Y=",E%Y," R=",E%R
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
             !messagelost="Lost in real kind=5 racetrack Aperture"
             write(messagelost,*) "Se_status.f90 CHECK_APERTURE_R : Lost in real kind=5 racetrack Aperture. ", &
                                  "Orbit: X=",X(1)," Y=",X(3)," Ap.: DX=",E%DX," DY=",E%DY," X=",E%X," Y=",E%Y," R=",E%R
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
    logical :: change_first=.true.

!    W_P=>W_I
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
    MYTYPE(kindwiggler)=" Sagan Wiggler"

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

!    change_sector=my_false

    IF(SECTOR_NMUL>0.and.firsttime_coef) THEN
     call alloc(e_muon_scale)
     e_muon_scale=1.0_dp
     call alloc(a_spin_scale)
     a_spin_scale=1.0_dp
       !  verb=global_verbose
       !  global_verbose=.false.
       if(firsttime_coef) THEN ! .or.(.not.allocated(S_B))) then
       
         if(   SECTOR_NMUL==11.and.sector_nmul_max==22.and.read_sector_info) then

             if(mcmillan) then
             call set_s_e_mcmillan
             call set_s_b_mcmillan
              else
             call set_s_b
             call set_s_e
             endif
         else 
 if(change_first) then  
  if(lielib_print(11)==1) write(6,*) " recomputing with new SECTOR_NMUL and sector_nmul_max ",SECTOR_NMUL,SECTOR_NMUL_max
 change_first=.false.
 endif        

          lda_old=lda_used
          lda_used=3000

             i=SECTOR_NMUL_MAX
             !             if(i==SECTOR_NMUL_MAX)     global_verbose=.true.
!             S_B%firsttime=0
!             call nul_coef(S_B)
!             call make_coef(S_B,I,0)  !1)
!             call get_bend_coeff(S_B,I)
             s_e%firsttime = 0  ! Piotr 8.19.2014 
!             s_eB%firsttime = 0 
             call nul_coef(s_e)
             call make_coef(s_e,I,0)
!             call make_coef(s_eB,I,0)
             call make_coef(S_B_from_V,I,0)
      !       call make_coef_e(s_e,I)
!             call get_bend_coeff(S_EB,I,1.D0,MY_TRUE)
             call get_bend_electric_coeff(s_e,I)
             call get_bend_magnetic_potential(S_B_from_V,I,1.D0,MY_TRUE)

!          ENDDO
          lda_used=lda_old
          if (global_verbose ) then
            call print_curv("Maxwellian_bend_for_ptc.txt",S_B_from_V)
            call print_curv_elec("Maxwellian_bend_for_ptc_electric.txt",s_e)
           !call print_curv_elec("Maxwellian_bend_mag_from_pot.txt",S_B_from_V)
          endif
          
       endif

       firsttime_coef=.FALSE.
    endif  !  calling preset routine
    else
         call kill(e_muon_scale)
         call alloc(e_muon_scale)
         call kill(a_spin_scale)
         call alloc(a_spin_scale)
     a_spin_scale=1.0_dp
     e_muon_scale=1.0_dp
    ENDIF
    call clear_states
    !  global_verbose=verb

  END  SUBROUTINE MAKE_STATES_0
  
 SUBROUTINE print_curv(filename,s_b)
    IMPLICIT NONE
    INTEGER I,J,nmul,mf
    character(*) filename
    character(255) line
    type(b_cyl), target :: s_b 


    call kanalnummer(mf,filename)

    
    nmul=SECTOR_NMUL_MAX
    write(mf,*) SECTOR_NMUL ,SECTOR_NMUL_MAX
         DO J=1,S_B_from_V%N_MONO
         write(line,*) "S_B_from_V%I(",J,")=",S_B_from_V%I(J),";","S_B_from_V%J(",J,")=",S_B_from_V%J(J),";"
         call context(line)
         write(mf,*) line(1:len_trim(line))
         enddo
    DO I=1,NMUL
       DO J=1,S_B_from_V%N_MONO
        if(S_B_from_V%A_X(I,J)/=0.0_dp) then
         write(line,*) "S_B_from_V%A_X(" ,I,",",J, ")=",S_B_from_V%A_X(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B_from_V%B_X(I,J)/=0.0_dp) then
         write(line,*) "S_B_from_V%B_X(" ,I,",",J, ")=",S_B_from_V%B_X(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B_from_V%A_y(I,J)/=0.0_dp) then
         write(line,*) "S_B_from_V%A_y(" ,I,",",J, ")=",S_B_from_V%A_y(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B_from_V%B_y(I,J)/=0.0_dp) then
         write(line,*) "S_B_from_V%B_y(" ,I,",",J, ")=",S_B_from_V%B_y(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
       ENDDO
    ENDDO
    

    DO I=1,NMUL
       DO J=1,S_B_from_V%N_MONO
        if(S_B_from_V%VA(I,J)/=0.0_dp) then
         write(line,*) "S_B_from_V%VA(" ,I,",",J, ")=",S_B_from_V%VA(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_B_from_V%VB(I,J)/=0.0_dp) then
         write(line,*) "S_B_from_V%VB(" ,I,",",J, ")=",S_B_from_V%VB(I,J),"e0_dp"
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

 SUBROUTINE print_curv_elec(filename,S_E)
    IMPLICIT NONE
    INTEGER I,J,nmul,mf
    character(*) filename
    character(255) line
    type(b_cyl), target :: S_E 


    call kanalnummer(mf,filename)

    
    nmul=SECTOR_NMUL_MAX
    write(mf,*) SECTOR_NMUL  ,SECTOR_NMUL_MAX
         DO J=1,S_E%N_MONO
         write(line,*) "S_E%I(",J,")=",S_E%I(J),";","S_E%J(",J,")=",S_E%J(J),";"
         call context(line)
         write(mf,*) line(1:len_trim(line))
         enddo
    DO I=1,NMUL
       DO J=1,S_E%N_MONO
        if(S_E%A_X(I,J)/=0.0_dp) then
         write(line,*) "S_E%A_X(" ,I,",",J, ")=",S_E%A_X(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_E%B_X(I,J)/=0.0_dp) then
         write(line,*) "S_E%B_X(" ,I,",",J, ")=",S_E%B_X(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_E%A_y(I,J)/=0.0_dp) then
         write(line,*) "S_E%A_y(" ,I,",",J, ")=",S_E%A_y(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_E%B_y(I,J)/=0.0_dp) then
         write(line,*) "S_E%B_y(" ,I,",",J, ")=",S_E%B_y(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
       ENDDO
    ENDDO

    DO I=1,NMUL
       DO J=1,S_E%N_MONO
        if(S_E%VA(I,J)/=0.0_dp) then
         write(line,*) "S_E%VA(" ,I,",",J, ")=",S_E%VA(I,J),"e0_dp"
         call context(line)
         if(index(line,"E-")/=0) then
          line(len_trim(line)-4:len_trim(line))="_DP  "
          call context(line)
         endif
         write(mf,*) line(1:len_trim(line))
        endif
        if(S_E%VB(I,J)/=0.0_dp) then
         write(line,*) "S_E%VB(" ,I,",",J, ")=",S_E%VB(I,J),"e0_dp"
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

  END SUBROUTINE print_curv_elec

  subroutine make_set_coef(b,no,ic)
    implicit none
    integer no
    integer  ic

    type(B_CYL) b

 !   ic=1
    b%firsttime=-100
    allocate(b%nmul)
    allocate(b%n_mono)
    b%nmul=no
    b%n_mono=((no+2-ic)*(no+1-ic))/2
    allocate(b%i(b%n_mono),b%j(b%n_mono))
    allocate(b%a_x(no,b%n_mono),b%a_y(no,b%n_mono))
    allocate(b%b_x(no,b%n_mono),b%b_y(no,b%n_mono))
    allocate(b%va(no,b%n_mono),b%vb(no,b%n_mono))
    b%i=0
    b%j=0
    b%a_x=0.0_dp
    b%b_x=0.0_dp
    b%a_y=0.0_dp
    b%b_y=0.0_dp
    b%va=0.0_dp
    b%vb=0.0_dp
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
    ONLY_2D=ONLY_2d0
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
         if(abs(1836.1526740143d0-muon)<1.d-8)then
          write(mf,*) "This is a proton"
         else
          write(mf,'((1X,a21,1x,G21.14,1x,A24))' ) "This a particle with ",muon, "times the electron mass "
         endif
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
    write(mf,'((1X,a20,1x,a5))' ) "      ONLY_2D     = ", CONV(S%ONLY_2D   )
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

  SUBROUTINE MAKE_STATES_m(muonfactor,ag,ne)
    USE   definition
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    real(dp) muonfactor,MASSF
    real(dp), optional :: ag,ne
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
    elseif(present(ag)) then
     a_particle=ag
    else 
     write(6,*) "Cannot do spin : provide a=g-2. Now it is set to zero."
    ENDIF
     initial_charge=1
    if(present(ne)) then
     initial_charge=ne
    endif
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

    ONLY_2D= ONLY_2D+DEFAULT

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
    S2%ONLY_2D=      S1%ONLY_2D
    S2%ONLY_4D=      S1%ONLY_4D
    S2%DELTA=       S1%DELTA
    S2%SPIN=       S1%SPIN
    S2%MODULATION=       S1%MODULATION
    S2%FULL_WAY=       S1%FULL_WAY
    !    S2%spin_dim=       S1%spin_dim
  END SUBROUTINE EQUALt





  SUBROUTINE  EQUALi(S2,i)
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
 
  END SUBROUTINE EQUALi



  FUNCTION add( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) add
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2

    if(s2%totalpath/=0.and.s2%totalpath/=1) then 
      add=s1
      return
    endif
    if(s1%totalpath/=0.and.s1%totalpath/=1) then 
      add=s1
      return
    endif
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
    add%ONLY_2D  =       S1%ONLY_2D.OR.S2%ONLY_2D
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
    IF(add%ONLY_2D) THEN
       add%TOTALPATH=  0
       add%RADIATION  =  F
       add%NOCAVITY =  T
       add%stochastic   =  F
       add%ENVELOPE   =  F
    ENDIF

    add%RADIATION  =  S1%RADIATION.OR.S2%RADIATION
    add%NOCAVITY =  S1%NOCAVITY.OR.S2%NOCAVITY
    add%TIME     =  S1%TIME.OR.S2%TIME
    add%FRINGE   =       S1%FRINGE.OR.S2%FRINGE
    add%stochastic   =       S1%stochastic.OR.S2%stochastic
    add%ENVELOPE   =       S1%ENVELOPE.OR.S2%ENVELOPE
    add%ONLY_2D  =       S1%ONLY_2D.OR.S2%ONLY_2D
    add%ONLY_4D  =       S1%ONLY_4D.OR.S2%ONLY_4D
    add%DELTA  =       S1%DELTA.OR.S2%DELTA
    add%SPIN  =       S1%SPIN.OR.S2%SPIN
    add%MODULATION  =       S1%MODULATION.OR.S2%MODULATION
    add%PARA_IN  =       S1%PARA_IN.OR.S2%PARA_IN.or.ALWAYS_knobs
    if(add%only_4d.and.add%only_2d) add%only_4d=my_false

    ADD%FULL_WAY=ADD%RADIATION.OR.ADD%stochastic.OR.ADD%ENVELOPE.OR.ADD%SPIN.OR.ADD%MODULATION
  END FUNCTION add

  FUNCTION sub( S1, S2 )
    implicit none
    TYPE (INTERNAL_STATE) sub
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1, S2
    logical(lp) dum1,dum2,tt1,tt2

    if(s2%totalpath/=0.and.s2%totalpath/=1) then 
      sub=s1
      return
    endif
    if(s1%totalpath/=0.and.s1%totalpath/=1) then 
      sub=s1
      return
    endif

    tt1=s1%only_2d
    tt2=s2%only_2d

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
    sub%ONLY_2D  =       S1%ONLY_2D.min.S2%ONLY_2D
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
        if(tt1.or.tt2) then
          sub%ONLY_2D=T
        else
          sub%ONLY_4D=T
        endif
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

    IF(sub%ONLY_2D) THEN
       sub%TOTALPATH=  0
       sub%RADIATION  =  F
       sub%ENVELOPE  =  F
       sub%stochastic  =  F
       sub%NOCAVITY =  T
       sub%stochastic   =  F
    ENDIF
    sub%FULL_WAY=sub%RADIATION.OR.sub%stochastic.OR.sub%ENVELOPE.OR.sub%SPIN.OR.sub%MODULATION
  END FUNCTION sub

  FUNCTION PARA_REMA(S1)   ! UNARY +
    implicit none
    TYPE (INTERNAL_STATE) PARA_REMA
    TYPE (INTERNAL_STATE), INTENT (IN) :: S1

    PARA_REMA         =    S1
    PARA_REMA%PARA_IN =     T

  END FUNCTION PARA_REMA

  subroutine init_all(STATE,NO1,NP1,pack,ND2,NPARA,number_of_clocks)
    !  subroutine S_init(STATE,NO1,NP1,PACKAGE,MAPINT,ND2,NPARA)
    implicit none
    TYPE (INTERNAL_STATE), INTENT(IN):: STATE
    LOGICAL(lp), optional, INTENT(IN):: pack
    INTEGER, INTENT(IN):: NO1,NP1
    INTEGER,optional :: ND2,NPARA,number_of_clocks

    use_complex_in_ptc=my_true
    call S_init(STATE,NO1,NP1,pack,ND2,NPARA,number_of_clocks)
   end subroutine init_all

  subroutine S_init(STATE,NO1,NP1,pack,ND2,NPARA,number_of_clocks)
    !  subroutine S_init(STATE,NO1,NP1,PACKAGE,MAPINT,ND2,NPARA)
    implicit none
    TYPE (INTERNAL_STATE), INTENT(IN):: STATE
    LOGICAL(lp), optional, INTENT(IN):: pack
    INTEGER, INTENT(IN):: NO1,NP1
    INTEGER ND1,NDEL,NDPT1
    INTEGER,optional :: ND2,NPARA,number_of_clocks
    INTEGER  ND2l,NPARAl,n_acc,no1c,nv,i
    LOGICAL(lp) package
    n_rf=0
!    call dd_p !valishev
    doing_ac_modulation_in_ptc=.false.
    package=my_true
    if(present(pack))     package=my_true
    only2d=0
    n_acc=0
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
       ELSEif(STATE%ONLY_2d) then
          only2d=2
          IF(STATE%DELTA) THEN
             ND1=1
             NDEL=1
             !             MAPINT=5
          ELSE
             ND1=1
             NDEL=0
             !            MAPINT=4
          ENDIF
       ELSE
          ND1=3
          NDEL=0
          NDPT1=5 + ndpt_bmad
          !         MAPINT=6
       ENDIF
    ELSE              ! CAVITY IN RING
       ND1=3
       NDPT1=0
       !       MAPINT=6
    ENDIF
 n_acc=0
    IF(STATE%modulation)  then
       doing_ac_modulation_in_ptc=.true.
      n_acc=1
     !  ND1=ND1+1
     if(present(number_of_clocks)) n_acc=number_of_clocks 
        !1
    endif

   ! IF(STATE%spin.or.STATE%modulation.or.STATE%radiation.or.STATE%envelope)  then
       if(automatic_complex) use_complex_in_ptc=.true.
   ! endif
    !    write(6,*) NO1,ND1,NP1,NDEL,NDPT1
    !pause 678

    CALL INIT(NO1,ND1,NP1+NDEL+2*n_acc,NDPT1,PACKAGE)
    nv=2*nd1+NP1+NDEL+2*n_acc



    ND2l=ND1*2+2*n_acc
    NPARAl=ND2l+NDEL
    C_%NPARA=NPARAl
    C_%ND2=ND2l
    C_%npara_fpp=NPARAl
    C_%SPIN_POS=0
    C_%NSPIN=0

    if(present(nd2)) nd2=nd2l
    if(present(npara)) npara=nparal
! etienne

no1c=no1+complex_extra_order
ND1=ND1+n_acc
    if(use_complex_in_ptc) call c_init(NO1c,nd1,np1+ndel,ndpt1,n_acc,ptc=my_false)  ! PTC false because we will not use the real FPP for acc modulation
    n_rf=n_acc
 
  END  subroutine S_init

  subroutine kill_map_cp()
    implicit none

    if(associated(dz_8)) then
      call kill(dz_8)
      deallocate(dz_8)
      nullify(dz_8)
    endif
    
    if(associated(dz_t)) then
      call kill(dz_t)
      deallocate(dz_t)
      nullify(dz_t)
    endif    

    
    if(associated(dz_c)) then
      call kill(dz_c)
      deallocate(dz_c)
      nullify(dz_c)
    endif    

  end subroutine kill_map_cp


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

  SUBROUTINE DTILTR_EXTERNAL(TILTD,I,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    INTEGER,INTENT(IN):: I 
    REAL(DP),INTENT(IN) :: TILTD
    real(dp) YS


    IF(TILTD==0.0_dp) RETURN
    IF(I==1) THEN
       ys=COS(TILTD)*x(1)+SIN(TILTD)*x(3)
       x(3)=COS(TILTD)*x(3)-SIN(TILTD)*x(1)
       x(1)=ys
       ys=COS(TILTD)*x(2)+SIN(TILTD)*x(4)
       x(4)=COS(TILTD)*x(4)-SIN(TILTD)*x(2)
       x(2)=ys
    ELSE
       ys=COS(TILTD)*x(1)-SIN(TILTD)*x(3)
       x(3)=COS(TILTD)*x(3)+SIN(TILTD)*x(1)
       x(1)=ys
       ys=COS(TILTD)*x(2)-SIN(TILTD)*x(4)
       x(4)=COS(TILTD)*x(4)+SIN(TILTD)*x(2)
       x(2)=ys
    ENDIF
 

  END SUBROUTINE DTILTR_EXTERNAL

  SUBROUTINE DTILTP_EXTERNAL(TILTD,I,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    INTEGER,INTENT(IN):: I 
    REAL(DP),INTENT(IN) :: TILTD
    TYPE(REAL_8) YS

    IF(TILTD==0.0_dp) RETURN
    CALL ALLOC(YS)

    IF(I==1) THEN
       ys=COS(TILTD)*x(1)+SIN(TILTD)*x(3)
       x(3)=COS(TILTD)*x(3)-SIN(TILTD)*x(1)
       x(1)=ys
       ys=COS(TILTD)*x(2)+SIN(TILTD)*x(4)
       x(4)=COS(TILTD)*x(4)-SIN(TILTD)*x(2)
       x(2)=ys
    ELSE
       ys=COS(TILTD)*x(1)-SIN(TILTD)*x(3)
       x(3)=COS(TILTD)*x(3)+SIN(TILTD)*x(1)
       x(1)=ys
       ys=COS(TILTD)*x(2)-SIN(TILTD)*x(4)
       x(4)=COS(TILTD)*x(4)+SIN(TILTD)*x(2)
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
  integer i
 


  !ALLOCATE(S_B(SECTOR_NMUL_MAX))
i=SECTOR_NMUL_MAX
          ! DO I=1,SECTOR_NMUL_MAX

             S_B_from_V%firsttime = 1  ! Piotr 8.19.2014 
             call nul_coef(S_B_from_V)
             call make_set_coef(S_B_from_V,I,0)
             S_B_from_V%firsttime=0
        !  ENDDO
S_B_FROM_V%I(1)=22;S_B_FROM_V%J(1)=0;
 S_B_FROM_V%I(2)=21;S_B_FROM_V%J(2)=1;
 S_B_FROM_V%I(3)=21;S_B_FROM_V%J(3)=0;
 S_B_FROM_V%I(4)=20;S_B_FROM_V%J(4)=2;
 S_B_FROM_V%I(5)=20;S_B_FROM_V%J(5)=1;
 S_B_FROM_V%I(6)=20;S_B_FROM_V%J(6)=0;
 S_B_FROM_V%I(7)=19;S_B_FROM_V%J(7)=3;
 S_B_FROM_V%I(8)=19;S_B_FROM_V%J(8)=2;
 S_B_FROM_V%I(9)=19;S_B_FROM_V%J(9)=1;
 S_B_FROM_V%I(10)=19;S_B_FROM_V%J(10)=0;
 S_B_FROM_V%I(11)=18;S_B_FROM_V%J(11)=4;
 S_B_FROM_V%I(12)=18;S_B_FROM_V%J(12)=3;
 S_B_FROM_V%I(13)=18;S_B_FROM_V%J(13)=2;
 S_B_FROM_V%I(14)=18;S_B_FROM_V%J(14)=1;
 S_B_FROM_V%I(15)=18;S_B_FROM_V%J(15)=0;
 S_B_FROM_V%I(16)=17;S_B_FROM_V%J(16)=5;
 S_B_FROM_V%I(17)=17;S_B_FROM_V%J(17)=4;
 S_B_FROM_V%I(18)=17;S_B_FROM_V%J(18)=3;
 S_B_FROM_V%I(19)=17;S_B_FROM_V%J(19)=2;
 S_B_FROM_V%I(20)=17;S_B_FROM_V%J(20)=1;
 S_B_FROM_V%I(21)=17;S_B_FROM_V%J(21)=0;
 S_B_FROM_V%I(22)=16;S_B_FROM_V%J(22)=6;
 S_B_FROM_V%I(23)=16;S_B_FROM_V%J(23)=5;
 S_B_FROM_V%I(24)=16;S_B_FROM_V%J(24)=4;
 S_B_FROM_V%I(25)=16;S_B_FROM_V%J(25)=3;
 S_B_FROM_V%I(26)=16;S_B_FROM_V%J(26)=2;
 S_B_FROM_V%I(27)=16;S_B_FROM_V%J(27)=1;
 S_B_FROM_V%I(28)=16;S_B_FROM_V%J(28)=0;
 S_B_FROM_V%I(29)=15;S_B_FROM_V%J(29)=7;
 S_B_FROM_V%I(30)=15;S_B_FROM_V%J(30)=6;
 S_B_FROM_V%I(31)=15;S_B_FROM_V%J(31)=5;
 S_B_FROM_V%I(32)=15;S_B_FROM_V%J(32)=4;
 S_B_FROM_V%I(33)=15;S_B_FROM_V%J(33)=3;
 S_B_FROM_V%I(34)=15;S_B_FROM_V%J(34)=2;
 S_B_FROM_V%I(35)=15;S_B_FROM_V%J(35)=1;
 S_B_FROM_V%I(36)=15;S_B_FROM_V%J(36)=0;
 S_B_FROM_V%I(37)=14;S_B_FROM_V%J(37)=8;
 S_B_FROM_V%I(38)=14;S_B_FROM_V%J(38)=7;
 S_B_FROM_V%I(39)=14;S_B_FROM_V%J(39)=6;
 S_B_FROM_V%I(40)=14;S_B_FROM_V%J(40)=5;
 S_B_FROM_V%I(41)=14;S_B_FROM_V%J(41)=4;
 S_B_FROM_V%I(42)=14;S_B_FROM_V%J(42)=3;
 S_B_FROM_V%I(43)=14;S_B_FROM_V%J(43)=2;
 S_B_FROM_V%I(44)=14;S_B_FROM_V%J(44)=1;
 S_B_FROM_V%I(45)=14;S_B_FROM_V%J(45)=0;
 S_B_FROM_V%I(46)=13;S_B_FROM_V%J(46)=9;
 S_B_FROM_V%I(47)=13;S_B_FROM_V%J(47)=8;
 S_B_FROM_V%I(48)=13;S_B_FROM_V%J(48)=7;
 S_B_FROM_V%I(49)=13;S_B_FROM_V%J(49)=6;
 S_B_FROM_V%I(50)=13;S_B_FROM_V%J(50)=5;
 S_B_FROM_V%I(51)=13;S_B_FROM_V%J(51)=4;
 S_B_FROM_V%I(52)=13;S_B_FROM_V%J(52)=3;
 S_B_FROM_V%I(53)=13;S_B_FROM_V%J(53)=2;
 S_B_FROM_V%I(54)=13;S_B_FROM_V%J(54)=1;
 S_B_FROM_V%I(55)=13;S_B_FROM_V%J(55)=0;
 S_B_FROM_V%I(56)=12;S_B_FROM_V%J(56)=10;
 S_B_FROM_V%I(57)=12;S_B_FROM_V%J(57)=9;
 S_B_FROM_V%I(58)=12;S_B_FROM_V%J(58)=8;
 S_B_FROM_V%I(59)=12;S_B_FROM_V%J(59)=7;
 S_B_FROM_V%I(60)=12;S_B_FROM_V%J(60)=6;
 S_B_FROM_V%I(61)=12;S_B_FROM_V%J(61)=5;
 S_B_FROM_V%I(62)=12;S_B_FROM_V%J(62)=4;
 S_B_FROM_V%I(63)=12;S_B_FROM_V%J(63)=3;
 S_B_FROM_V%I(64)=12;S_B_FROM_V%J(64)=2;
 S_B_FROM_V%I(65)=12;S_B_FROM_V%J(65)=1;
 S_B_FROM_V%I(66)=12;S_B_FROM_V%J(66)=0;
 S_B_FROM_V%I(67)=11;S_B_FROM_V%J(67)=11;
 S_B_FROM_V%I(68)=11;S_B_FROM_V%J(68)=10;
 S_B_FROM_V%I(69)=11;S_B_FROM_V%J(69)=9;
 S_B_FROM_V%I(70)=11;S_B_FROM_V%J(70)=8;
 S_B_FROM_V%I(71)=11;S_B_FROM_V%J(71)=7;
 S_B_FROM_V%I(72)=11;S_B_FROM_V%J(72)=6;
 S_B_FROM_V%I(73)=11;S_B_FROM_V%J(73)=5;
 S_B_FROM_V%I(74)=11;S_B_FROM_V%J(74)=4;
 S_B_FROM_V%I(75)=11;S_B_FROM_V%J(75)=3;
 S_B_FROM_V%I(76)=11;S_B_FROM_V%J(76)=2;
 S_B_FROM_V%I(77)=11;S_B_FROM_V%J(77)=1;
 S_B_FROM_V%I(78)=11;S_B_FROM_V%J(78)=0;
 S_B_FROM_V%I(79)=10;S_B_FROM_V%J(79)=12;
 S_B_FROM_V%I(80)=10;S_B_FROM_V%J(80)=11;
 S_B_FROM_V%I(81)=10;S_B_FROM_V%J(81)=10;
 S_B_FROM_V%I(82)=10;S_B_FROM_V%J(82)=9;
 S_B_FROM_V%I(83)=10;S_B_FROM_V%J(83)=8;
 S_B_FROM_V%I(84)=10;S_B_FROM_V%J(84)=7;
 S_B_FROM_V%I(85)=10;S_B_FROM_V%J(85)=6;
 S_B_FROM_V%I(86)=10;S_B_FROM_V%J(86)=5;
 S_B_FROM_V%I(87)=10;S_B_FROM_V%J(87)=4;
 S_B_FROM_V%I(88)=10;S_B_FROM_V%J(88)=3;
 S_B_FROM_V%I(89)=10;S_B_FROM_V%J(89)=2;
 S_B_FROM_V%I(90)=10;S_B_FROM_V%J(90)=1;
 S_B_FROM_V%I(91)=10;S_B_FROM_V%J(91)=0;
 S_B_FROM_V%I(92)=9;S_B_FROM_V%J(92)=13;
 S_B_FROM_V%I(93)=9;S_B_FROM_V%J(93)=12;
 S_B_FROM_V%I(94)=9;S_B_FROM_V%J(94)=11;
 S_B_FROM_V%I(95)=9;S_B_FROM_V%J(95)=10;
 S_B_FROM_V%I(96)=9;S_B_FROM_V%J(96)=9;
 S_B_FROM_V%I(97)=9;S_B_FROM_V%J(97)=8;
 S_B_FROM_V%I(98)=9;S_B_FROM_V%J(98)=7;
 S_B_FROM_V%I(99)=9;S_B_FROM_V%J(99)=6;
 S_B_FROM_V%I(100)=9;S_B_FROM_V%J(100)=5;
 S_B_FROM_V%I(101)=9;S_B_FROM_V%J(101)=4;
 S_B_FROM_V%I(102)=9;S_B_FROM_V%J(102)=3;
 S_B_FROM_V%I(103)=9;S_B_FROM_V%J(103)=2;
 S_B_FROM_V%I(104)=9;S_B_FROM_V%J(104)=1;
 S_B_FROM_V%I(105)=9;S_B_FROM_V%J(105)=0;
 S_B_FROM_V%I(106)=8;S_B_FROM_V%J(106)=14;
 S_B_FROM_V%I(107)=8;S_B_FROM_V%J(107)=13;
 S_B_FROM_V%I(108)=8;S_B_FROM_V%J(108)=12;
 S_B_FROM_V%I(109)=8;S_B_FROM_V%J(109)=11;
 S_B_FROM_V%I(110)=8;S_B_FROM_V%J(110)=10;
 S_B_FROM_V%I(111)=8;S_B_FROM_V%J(111)=9;
 S_B_FROM_V%I(112)=8;S_B_FROM_V%J(112)=8;
 S_B_FROM_V%I(113)=8;S_B_FROM_V%J(113)=7;
 S_B_FROM_V%I(114)=8;S_B_FROM_V%J(114)=6;
 S_B_FROM_V%I(115)=8;S_B_FROM_V%J(115)=5;
 S_B_FROM_V%I(116)=8;S_B_FROM_V%J(116)=4;
 S_B_FROM_V%I(117)=8;S_B_FROM_V%J(117)=3;
 S_B_FROM_V%I(118)=8;S_B_FROM_V%J(118)=2;
 S_B_FROM_V%I(119)=8;S_B_FROM_V%J(119)=1;
 S_B_FROM_V%I(120)=8;S_B_FROM_V%J(120)=0;
 S_B_FROM_V%I(121)=7;S_B_FROM_V%J(121)=15;
 S_B_FROM_V%I(122)=7;S_B_FROM_V%J(122)=14;
 S_B_FROM_V%I(123)=7;S_B_FROM_V%J(123)=13;
 S_B_FROM_V%I(124)=7;S_B_FROM_V%J(124)=12;
 S_B_FROM_V%I(125)=7;S_B_FROM_V%J(125)=11;
 S_B_FROM_V%I(126)=7;S_B_FROM_V%J(126)=10;
 S_B_FROM_V%I(127)=7;S_B_FROM_V%J(127)=9;
 S_B_FROM_V%I(128)=7;S_B_FROM_V%J(128)=8;
 S_B_FROM_V%I(129)=7;S_B_FROM_V%J(129)=7;
 S_B_FROM_V%I(130)=7;S_B_FROM_V%J(130)=6;
 S_B_FROM_V%I(131)=7;S_B_FROM_V%J(131)=5;
 S_B_FROM_V%I(132)=7;S_B_FROM_V%J(132)=4;
 S_B_FROM_V%I(133)=7;S_B_FROM_V%J(133)=3;
 S_B_FROM_V%I(134)=7;S_B_FROM_V%J(134)=2;
 S_B_FROM_V%I(135)=7;S_B_FROM_V%J(135)=1;
 S_B_FROM_V%I(136)=7;S_B_FROM_V%J(136)=0;
 S_B_FROM_V%I(137)=6;S_B_FROM_V%J(137)=16;
 S_B_FROM_V%I(138)=6;S_B_FROM_V%J(138)=15;
 S_B_FROM_V%I(139)=6;S_B_FROM_V%J(139)=14;
 S_B_FROM_V%I(140)=6;S_B_FROM_V%J(140)=13;
 S_B_FROM_V%I(141)=6;S_B_FROM_V%J(141)=12;
 S_B_FROM_V%I(142)=6;S_B_FROM_V%J(142)=11;
 S_B_FROM_V%I(143)=6;S_B_FROM_V%J(143)=10;
 S_B_FROM_V%I(144)=6;S_B_FROM_V%J(144)=9;
 S_B_FROM_V%I(145)=6;S_B_FROM_V%J(145)=8;
 S_B_FROM_V%I(146)=6;S_B_FROM_V%J(146)=7;
 S_B_FROM_V%I(147)=6;S_B_FROM_V%J(147)=6;
 S_B_FROM_V%I(148)=6;S_B_FROM_V%J(148)=5;
 S_B_FROM_V%I(149)=6;S_B_FROM_V%J(149)=4;
 S_B_FROM_V%I(150)=6;S_B_FROM_V%J(150)=3;
 S_B_FROM_V%I(151)=6;S_B_FROM_V%J(151)=2;
 S_B_FROM_V%I(152)=6;S_B_FROM_V%J(152)=1;
 S_B_FROM_V%I(153)=6;S_B_FROM_V%J(153)=0;
 S_B_FROM_V%I(154)=5;S_B_FROM_V%J(154)=17;
 S_B_FROM_V%I(155)=5;S_B_FROM_V%J(155)=16;
 S_B_FROM_V%I(156)=5;S_B_FROM_V%J(156)=15;
 S_B_FROM_V%I(157)=5;S_B_FROM_V%J(157)=14;
 S_B_FROM_V%I(158)=5;S_B_FROM_V%J(158)=13;
 S_B_FROM_V%I(159)=5;S_B_FROM_V%J(159)=12;
 S_B_FROM_V%I(160)=5;S_B_FROM_V%J(160)=11;
 S_B_FROM_V%I(161)=5;S_B_FROM_V%J(161)=10;
 S_B_FROM_V%I(162)=5;S_B_FROM_V%J(162)=9;
 S_B_FROM_V%I(163)=5;S_B_FROM_V%J(163)=8;
 S_B_FROM_V%I(164)=5;S_B_FROM_V%J(164)=7;
 S_B_FROM_V%I(165)=5;S_B_FROM_V%J(165)=6;
 S_B_FROM_V%I(166)=5;S_B_FROM_V%J(166)=5;
 S_B_FROM_V%I(167)=5;S_B_FROM_V%J(167)=4;
 S_B_FROM_V%I(168)=5;S_B_FROM_V%J(168)=3;
 S_B_FROM_V%I(169)=5;S_B_FROM_V%J(169)=2;
 S_B_FROM_V%I(170)=5;S_B_FROM_V%J(170)=1;
 S_B_FROM_V%I(171)=5;S_B_FROM_V%J(171)=0;
 S_B_FROM_V%I(172)=4;S_B_FROM_V%J(172)=18;
 S_B_FROM_V%I(173)=4;S_B_FROM_V%J(173)=17;
 S_B_FROM_V%I(174)=4;S_B_FROM_V%J(174)=16;
 S_B_FROM_V%I(175)=4;S_B_FROM_V%J(175)=15;
 S_B_FROM_V%I(176)=4;S_B_FROM_V%J(176)=14;
 S_B_FROM_V%I(177)=4;S_B_FROM_V%J(177)=13;
 S_B_FROM_V%I(178)=4;S_B_FROM_V%J(178)=12;
 S_B_FROM_V%I(179)=4;S_B_FROM_V%J(179)=11;
 S_B_FROM_V%I(180)=4;S_B_FROM_V%J(180)=10;
 S_B_FROM_V%I(181)=4;S_B_FROM_V%J(181)=9;
 S_B_FROM_V%I(182)=4;S_B_FROM_V%J(182)=8;
 S_B_FROM_V%I(183)=4;S_B_FROM_V%J(183)=7;
 S_B_FROM_V%I(184)=4;S_B_FROM_V%J(184)=6;
 S_B_FROM_V%I(185)=4;S_B_FROM_V%J(185)=5;
 S_B_FROM_V%I(186)=4;S_B_FROM_V%J(186)=4;
 S_B_FROM_V%I(187)=4;S_B_FROM_V%J(187)=3;
 S_B_FROM_V%I(188)=4;S_B_FROM_V%J(188)=2;
 S_B_FROM_V%I(189)=4;S_B_FROM_V%J(189)=1;
 S_B_FROM_V%I(190)=4;S_B_FROM_V%J(190)=0;
 S_B_FROM_V%I(191)=3;S_B_FROM_V%J(191)=19;
 S_B_FROM_V%I(192)=3;S_B_FROM_V%J(192)=18;
 S_B_FROM_V%I(193)=3;S_B_FROM_V%J(193)=17;
 S_B_FROM_V%I(194)=3;S_B_FROM_V%J(194)=16;
 S_B_FROM_V%I(195)=3;S_B_FROM_V%J(195)=15;
 S_B_FROM_V%I(196)=3;S_B_FROM_V%J(196)=14;
 S_B_FROM_V%I(197)=3;S_B_FROM_V%J(197)=13;
 S_B_FROM_V%I(198)=3;S_B_FROM_V%J(198)=12;
 S_B_FROM_V%I(199)=3;S_B_FROM_V%J(199)=11;
 S_B_FROM_V%I(200)=3;S_B_FROM_V%J(200)=10;
 S_B_FROM_V%I(201)=3;S_B_FROM_V%J(201)=9;
 S_B_FROM_V%I(202)=3;S_B_FROM_V%J(202)=8;
 S_B_FROM_V%I(203)=3;S_B_FROM_V%J(203)=7;
 S_B_FROM_V%I(204)=3;S_B_FROM_V%J(204)=6;
 S_B_FROM_V%I(205)=3;S_B_FROM_V%J(205)=5;
 S_B_FROM_V%I(206)=3;S_B_FROM_V%J(206)=4;
 S_B_FROM_V%I(207)=3;S_B_FROM_V%J(207)=3;
 S_B_FROM_V%I(208)=3;S_B_FROM_V%J(208)=2;
 S_B_FROM_V%I(209)=3;S_B_FROM_V%J(209)=1;
 S_B_FROM_V%I(210)=3;S_B_FROM_V%J(210)=0;
 S_B_FROM_V%I(211)=2;S_B_FROM_V%J(211)=20;
 S_B_FROM_V%I(212)=2;S_B_FROM_V%J(212)=19;
 S_B_FROM_V%I(213)=2;S_B_FROM_V%J(213)=18;
 S_B_FROM_V%I(214)=2;S_B_FROM_V%J(214)=17;
 S_B_FROM_V%I(215)=2;S_B_FROM_V%J(215)=16;
 S_B_FROM_V%I(216)=2;S_B_FROM_V%J(216)=15;
 S_B_FROM_V%I(217)=2;S_B_FROM_V%J(217)=14;
 S_B_FROM_V%I(218)=2;S_B_FROM_V%J(218)=13;
 S_B_FROM_V%I(219)=2;S_B_FROM_V%J(219)=12;
 S_B_FROM_V%I(220)=2;S_B_FROM_V%J(220)=11;
 S_B_FROM_V%I(221)=2;S_B_FROM_V%J(221)=10;
 S_B_FROM_V%I(222)=2;S_B_FROM_V%J(222)=9;
 S_B_FROM_V%I(223)=2;S_B_FROM_V%J(223)=8;
 S_B_FROM_V%I(224)=2;S_B_FROM_V%J(224)=7;
 S_B_FROM_V%I(225)=2;S_B_FROM_V%J(225)=6;
 S_B_FROM_V%I(226)=2;S_B_FROM_V%J(226)=5;
 S_B_FROM_V%I(227)=2;S_B_FROM_V%J(227)=4;
 S_B_FROM_V%I(228)=2;S_B_FROM_V%J(228)=3;
 S_B_FROM_V%I(229)=2;S_B_FROM_V%J(229)=2;
 S_B_FROM_V%I(230)=2;S_B_FROM_V%J(230)=1;
 S_B_FROM_V%I(231)=2;S_B_FROM_V%J(231)=0;
 S_B_FROM_V%I(232)=1;S_B_FROM_V%J(232)=21;
 S_B_FROM_V%I(233)=1;S_B_FROM_V%J(233)=20;
 S_B_FROM_V%I(234)=1;S_B_FROM_V%J(234)=19;
 S_B_FROM_V%I(235)=1;S_B_FROM_V%J(235)=18;
 S_B_FROM_V%I(236)=1;S_B_FROM_V%J(236)=17;
 S_B_FROM_V%I(237)=1;S_B_FROM_V%J(237)=16;
 S_B_FROM_V%I(238)=1;S_B_FROM_V%J(238)=15;
 S_B_FROM_V%I(239)=1;S_B_FROM_V%J(239)=14;
 S_B_FROM_V%I(240)=1;S_B_FROM_V%J(240)=13;
 S_B_FROM_V%I(241)=1;S_B_FROM_V%J(241)=12;
 S_B_FROM_V%I(242)=1;S_B_FROM_V%J(242)=11;
 S_B_FROM_V%I(243)=1;S_B_FROM_V%J(243)=10;
 S_B_FROM_V%I(244)=1;S_B_FROM_V%J(244)=9;
 S_B_FROM_V%I(245)=1;S_B_FROM_V%J(245)=8;
 S_B_FROM_V%I(246)=1;S_B_FROM_V%J(246)=7;
 S_B_FROM_V%I(247)=1;S_B_FROM_V%J(247)=6;
 S_B_FROM_V%I(248)=1;S_B_FROM_V%J(248)=5;
 S_B_FROM_V%I(249)=1;S_B_FROM_V%J(249)=4;
 S_B_FROM_V%I(250)=1;S_B_FROM_V%J(250)=3;
 S_B_FROM_V%I(251)=1;S_B_FROM_V%J(251)=2;
 S_B_FROM_V%I(252)=1;S_B_FROM_V%J(252)=1;
 S_B_FROM_V%I(253)=1;S_B_FROM_V%J(253)=0;
 S_B_FROM_V%I(254)=0;S_B_FROM_V%J(254)=22;
 S_B_FROM_V%I(255)=0;S_B_FROM_V%J(255)=21;
 S_B_FROM_V%I(256)=0;S_B_FROM_V%J(256)=20;
 S_B_FROM_V%I(257)=0;S_B_FROM_V%J(257)=19;
 S_B_FROM_V%I(258)=0;S_B_FROM_V%J(258)=18;
 S_B_FROM_V%I(259)=0;S_B_FROM_V%J(259)=17;
 S_B_FROM_V%I(260)=0;S_B_FROM_V%J(260)=16;
 S_B_FROM_V%I(261)=0;S_B_FROM_V%J(261)=15;
 S_B_FROM_V%I(262)=0;S_B_FROM_V%J(262)=14;
 S_B_FROM_V%I(263)=0;S_B_FROM_V%J(263)=13;
 S_B_FROM_V%I(264)=0;S_B_FROM_V%J(264)=12;
 S_B_FROM_V%I(265)=0;S_B_FROM_V%J(265)=11;
 S_B_FROM_V%I(266)=0;S_B_FROM_V%J(266)=10;
 S_B_FROM_V%I(267)=0;S_B_FROM_V%J(267)=9;
 S_B_FROM_V%I(268)=0;S_B_FROM_V%J(268)=8;
 S_B_FROM_V%I(269)=0;S_B_FROM_V%J(269)=7;
 S_B_FROM_V%I(270)=0;S_B_FROM_V%J(270)=6;
 S_B_FROM_V%I(271)=0;S_B_FROM_V%J(271)=5;
 S_B_FROM_V%I(272)=0;S_B_FROM_V%J(272)=4;
 S_B_FROM_V%I(273)=0;S_B_FROM_V%J(273)=3;
 S_B_FROM_V%I(274)=0;S_B_FROM_V%J(274)=2;
 S_B_FROM_V%I(275)=0;S_B_FROM_V%J(275)=1;
 S_B_FROM_V%I(276)=0;S_B_FROM_V%J(276)=0;
 S_B_FROM_V%A_Y(1,104)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,118)=4.50000000000000E0_DP
 S_B_FROM_V%A_Y(1,119)=-0.999999999999999E0_DP
 S_B_FROM_V%A_Y(1,133)=-6.00000000000000E0_DP
 S_B_FROM_V%A_X(1,134)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,135)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,149)=-10.5000000000000E0_DP
 S_B_FROM_V%A_Y(1,150)=4.66666666666666E0_DP
 S_B_FROM_V%A_X(1,151)=3.50000000000000E0_DP
 S_B_FROM_V%A_Y(1,152)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,166)=9.45000000000000E0_DP
 S_B_FROM_V%A_X(1,167)=6.99999999999999E0_DP
 S_B_FROM_V%A_Y(1,168)=-3.50000000000000E0_DP
 S_B_FROM_V%A_X(1,169)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,170)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,184)=7.87500000000000E0_DP
 S_B_FROM_V%A_Y(1,185)=-5.24999999999999E0_DP
 S_B_FROM_V%A_X(1,186)=-4.37500000000000E0_DP
 S_B_FROM_V%A_Y(1,187)=2.50000000000000E0_DP
 S_B_FROM_V%A_X(1,188)=2.50000000000000E0_DP
 S_B_FROM_V%A_Y(1,189)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,203)=-3.75000000000000E0_DP
 S_B_FROM_V%A_X(1,204)=-3.49999999999999E0_DP
 S_B_FROM_V%A_Y(1,205)=2.62500000000000E0_DP
 S_B_FROM_V%A_X(1,206)=2.50000000000000E0_DP
 S_B_FROM_V%A_Y(1,207)=-1.66666666666667E0_DP
 S_B_FROM_V%A_X(1,208)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,209)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,223)=-1.40625000000000E0_DP
 S_B_FROM_V%A_Y(1,224)=1.25000000000000E0_DP
 S_B_FROM_V%A_X(1,225)=1.31250000000000E0_DP
 S_B_FROM_V%A_Y(1,226)=-1.12500000000000E0_DP
 S_B_FROM_V%A_X(1,227)=-1.25000000000000E0_DP
 S_B_FROM_V%A_Y(1,228)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,229)=1.50000000000000E0_DP
 S_B_FROM_V%A_Y(1,230)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,244)=0.273437500000000E0_DP
 S_B_FROM_V%A_X(1,245)=0.312499999999999E0_DP
 S_B_FROM_V%A_Y(1,246)=-0.312500000000000E0_DP
 S_B_FROM_V%A_X(1,247)=-0.375000000000000E0_DP
 S_B_FROM_V%A_Y(1,248)=0.375000000000000E0_DP
 S_B_FROM_V%A_X(1,249)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(1,250)=-0.500000000000000E0_DP
 S_B_FROM_V%A_X(1,251)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(1,252)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,266)=2.734375000000002E-002_DP
 S_B_FROM_V%A_Y(1,267)=-3.038194444444437E-002_DP
 S_B_FROM_V%A_X(1,268)=-3.906250000000002E-002_DP
 S_B_FROM_V%A_Y(1,269)=4.464285714285714E-002_DP
 S_B_FROM_V%A_X(1,270)=6.249999999999997E-002_DP
 S_B_FROM_V%A_Y(1,271)=-7.500000000000008E-002_DP
 S_B_FROM_V%A_X(1,272)=-0.125000000000000E0_DP
 S_B_FROM_V%A_Y(1,273)=0.166666666666667E0_DP
 S_B_FROM_V%A_X(1,274)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(1,275)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,276)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(1,276)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(2,104)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,118)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(2,118)=-0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,119)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(2,133)=-1.33333333333333E0_DP
 S_B_FROM_V%A_Y(2,133)=6.00000000000000E0_DP
 S_B_FROM_V%A_X(2,134)=4.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,134)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,135)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,149)=10.5000000000000E0_DP
 S_B_FROM_V%B_Y(2,149)=1.16666666666667E0_DP
 S_B_FROM_V%B_X(2,150)=1.16666666666667E0_DP
 S_B_FROM_V%A_Y(2,150)=-4.66666666666667E0_DP
 S_B_FROM_V%A_X(2,151)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(2,151)=-0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,152)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(2,166)=1.40000000000000E0_DP
 S_B_FROM_V%A_Y(2,166)=-9.44999999999999E0_DP
 S_B_FROM_V%A_X(2,167)=-7.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,167)=-0.875000000000000E0_DP
 S_B_FROM_V%B_X(2,168)=-0.999999999999999E0_DP
 S_B_FROM_V%A_Y(2,168)=3.50000000000000E0_DP
 S_B_FROM_V%A_X(2,169)=3.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,169)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,170)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,184)=-7.87499999999999E0_DP
 S_B_FROM_V%B_Y(2,184)=-0.875000000000000E0_DP
 S_B_FROM_V%B_X(2,185)=-0.875000000000000E0_DP
 S_B_FROM_V%A_Y(2,185)=5.25000000000000E0_DP
 S_B_FROM_V%A_X(2,186)=4.37500000000000E0_DP
 S_B_FROM_V%B_Y(2,186)=0.625000000000000E0_DP
 S_B_FROM_V%B_X(2,187)=0.833333333333333E0_DP
 S_B_FROM_V%A_Y(2,187)=-2.50000000000000E0_DP
 S_B_FROM_V%A_X(2,188)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(2,188)=-0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,189)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(2,203)=-0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,203)=3.75000000000000E0_DP
 S_B_FROM_V%A_X(2,204)=3.50000000000000E0_DP
 S_B_FROM_V%B_Y(2,204)=0.437500000000000E0_DP
 S_B_FROM_V%B_X(2,205)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,205)=-2.62500000000000E0_DP
 S_B_FROM_V%A_X(2,206)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(2,206)=-0.416666666666667E0_DP
 S_B_FROM_V%B_X(2,207)=-0.666666666666667E0_DP
 S_B_FROM_V%A_Y(2,207)=1.66666666666667E0_DP
 S_B_FROM_V%A_X(2,208)=2.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,208)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,209)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,223)=1.40625000000000E0_DP
 S_B_FROM_V%B_Y(2,223)=0.156250000000000E0_DP
 S_B_FROM_V%B_X(2,224)=0.187500000000000E0_DP
 S_B_FROM_V%A_Y(2,224)=-1.25000000000000E0_DP
 S_B_FROM_V%A_X(2,225)=-1.31250000000000E0_DP
 S_B_FROM_V%B_Y(2,225)=-0.187500000000000E0_DP
 S_B_FROM_V%B_X(2,226)=-0.250000000000000E0_DP
 S_B_FROM_V%A_Y(2,226)=1.12500000000000E0_DP
 S_B_FROM_V%A_X(2,227)=1.25000000000000E0_DP
 S_B_FROM_V%B_Y(2,227)=0.250000000000000E0_DP
 S_B_FROM_V%B_X(2,228)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,228)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,229)=-1.50000000000000E0_DP
 S_B_FROM_V%B_Y(2,229)=-0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,230)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(2,244)=3.472222222222226E-002_DP
 S_B_FROM_V%A_Y(2,244)=-0.273437500000000E0_DP
 S_B_FROM_V%A_X(2,245)=-0.312500000000000E0_DP
 S_B_FROM_V%B_Y(2,245)=-3.906249999999997E-002_DP
 S_B_FROM_V%B_X(2,246)=-5.357142857142853E-002_DP
 S_B_FROM_V%A_Y(2,246)=0.312500000000000E0_DP
 S_B_FROM_V%A_X(2,247)=0.375000000000000E0_DP
 S_B_FROM_V%B_Y(2,247)=6.250000000000000E-002_DP
 S_B_FROM_V%B_X(2,248)=0.100000000000000E0_DP
 S_B_FROM_V%A_Y(2,248)=-0.375000000000000E0_DP
 S_B_FROM_V%A_X(2,249)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,249)=-0.125000000000000E0_DP
 S_B_FROM_V%B_X(2,250)=-0.333333333333333E0_DP
 S_B_FROM_V%A_Y(2,250)=0.500000000000000E0_DP
 S_B_FROM_V%A_X(2,251)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,251)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(2,252)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,253)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,253)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,266)=-2.734374999999996E-002_DP
 S_B_FROM_V%B_Y(2,266)=-3.038194444444450E-003_DP
 S_B_FROM_V%B_X(2,267)=-4.340277777777775E-003_DP
 S_B_FROM_V%A_Y(2,267)=3.038194444444442E-002_DP
 S_B_FROM_V%A_X(2,268)=3.906250000000003E-002_DP
 S_B_FROM_V%B_Y(2,268)=5.580357142857133E-003_DP
 S_B_FROM_V%B_X(2,269)=8.928571428571428E-003_DP
 S_B_FROM_V%A_Y(2,269)=-4.464285714285711E-002_DP
 S_B_FROM_V%A_X(2,270)=-6.249999999999998E-002_DP
 S_B_FROM_V%B_Y(2,270)=-1.249999999999999E-002_DP
 S_B_FROM_V%B_X(2,271)=-2.500000000000000E-002_DP
 S_B_FROM_V%A_Y(2,271)=7.500000000000002E-002_DP
 S_B_FROM_V%A_X(2,272)=0.125000000000000E0_DP
 S_B_FROM_V%B_Y(2,272)=4.166666666666667E-002_DP
 S_B_FROM_V%B_X(2,273)=0.166666666666667E0_DP
 S_B_FROM_V%A_Y(2,273)=-0.166666666666667E0_DP
 S_B_FROM_V%A_X(2,274)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,274)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(2,275)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(2,275)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,104)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,118)=4.50000000000000E0_DP
 S_B_FROM_V%B_Y(3,118)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,119)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,133)=2.66666666666667E0_DP
 S_B_FROM_V%A_Y(3,133)=-6.50000000000000E0_DP
 S_B_FROM_V%A_X(3,134)=-4.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,134)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,135)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,149)=-11.3750000000000E0_DP
 S_B_FROM_V%B_Y(3,149)=-2.33333333333333E0_DP
 S_B_FROM_V%B_X(3,150)=-2.33333333333333E0_DP
 S_B_FROM_V%A_Y(3,150)=5.16666666666666E0_DP
 S_B_FROM_V%A_X(3,151)=3.50000000000000E0_DP
 S_B_FROM_V%B_Y(3,151)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,152)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,166)=-2.80000000000000E0_DP
 S_B_FROM_V%A_Y(3,166)=9.97499999999999E0_DP
 S_B_FROM_V%A_X(3,167)=7.75000000000000E0_DP
 S_B_FROM_V%B_Y(3,167)=1.75000000000000E0_DP
 S_B_FROM_V%B_X(3,168)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,168)=-4.00000000000000E0_DP
 S_B_FROM_V%A_X(3,169)=-3.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,169)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,170)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,184)=8.31250000000000E0_DP
 S_B_FROM_V%B_Y(3,184)=1.75000000000000E0_DP
 S_B_FROM_V%B_X(3,185)=1.75000000000000E0_DP
 S_B_FROM_V%A_Y(3,185)=-5.62499999999999E0_DP
 S_B_FROM_V%A_X(3,186)=-5.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,186)=-1.25000000000000E0_DP
 S_B_FROM_V%B_X(3,187)=-1.66666666666667E0_DP
 S_B_FROM_V%A_Y(3,187)=3.00000000000000E0_DP
 S_B_FROM_V%A_X(3,188)=2.50000000000000E0_DP
 S_B_FROM_V%B_Y(3,188)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,189)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,203)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,203)=-3.93750000000000E0_DP
 S_B_FROM_V%A_X(3,204)=-3.75000000000000E0_DP
 S_B_FROM_V%B_Y(3,204)=-0.875000000000000E0_DP
 S_B_FROM_V%B_X(3,205)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,205)=2.87500000000000E0_DP
 S_B_FROM_V%A_X(3,206)=3.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,206)=0.833333333333333E0_DP
 S_B_FROM_V%B_X(3,207)=1.33333333333333E0_DP
 S_B_FROM_V%A_Y(3,207)=-2.16666666666667E0_DP
 S_B_FROM_V%A_X(3,208)=-2.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,208)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,209)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,223)=-1.47656250000000E0_DP
 S_B_FROM_V%B_Y(3,223)=-0.312500000000000E0_DP
 S_B_FROM_V%B_X(3,224)=-0.375000000000000E0_DP
 S_B_FROM_V%A_Y(3,224)=1.33035714285714E0_DP
 S_B_FROM_V%A_X(3,225)=1.43750000000000E0_DP
 S_B_FROM_V%B_Y(3,225)=0.375000000000000E0_DP
 S_B_FROM_V%B_X(3,226)=0.500000000000000E0_DP
 S_B_FROM_V%A_Y(3,226)=-1.27500000000000E0_DP
 S_B_FROM_V%A_X(3,227)=-1.62500000000000E0_DP
 S_B_FROM_V%B_Y(3,227)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(3,228)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,228)=1.50000000000000E0_DP
 S_B_FROM_V%A_X(3,229)=1.50000000000000E0_DP
 S_B_FROM_V%B_Y(3,229)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,230)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,231)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,231)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,244)=-6.944444444444439E-002_DP
 S_B_FROM_V%A_Y(3,244)=0.286458333333333E0_DP
 S_B_FROM_V%A_X(3,245)=0.332589285714285E0_DP
 S_B_FROM_V%B_Y(3,245)=7.812500000000000E-002_DP
 S_B_FROM_V%B_X(3,246)=0.107142857142857E0_DP
 S_B_FROM_V%A_Y(3,246)=-0.339285714285715E0_DP
 S_B_FROM_V%A_X(3,247)=-0.425000000000000E0_DP
 S_B_FROM_V%B_Y(3,247)=-0.125000000000000E0_DP
 S_B_FROM_V%B_X(3,248)=-0.200000000000000E0_DP
 S_B_FROM_V%A_Y(3,248)=0.450000000000000E0_DP
 S_B_FROM_V%A_X(3,249)=0.750000000000000E0_DP
 S_B_FROM_V%B_Y(3,249)=0.250000000000000E0_DP
 S_B_FROM_V%B_X(3,250)=0.666666666666667E0_DP
 S_B_FROM_V%A_Y(3,250)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,251)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,251)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,252)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,252)=-2.00000000000000E0_DP
 S_B_FROM_V%A_X(3,266)=2.864583333333333E-002_DP
 S_B_FROM_V%B_Y(3,266)=6.076388888888886E-003_DP
 S_B_FROM_V%B_X(3,267)=8.680555555555556E-003_DP
 S_B_FROM_V%A_Y(3,267)=-3.224206349206343E-002_DP
 S_B_FROM_V%A_X(3,268)=-4.241071428571433E-002_DP
 S_B_FROM_V%B_Y(3,268)=-1.116071428571428E-002_DP
 S_B_FROM_V%B_X(3,269)=-1.785714285714286E-002_DP
 S_B_FROM_V%A_Y(3,269)=5.000000000000000E-002_DP
 S_B_FROM_V%A_X(3,270)=7.499999999999998E-002_DP
 S_B_FROM_V%B_Y(3,270)=2.500000000000000E-002_DP
 S_B_FROM_V%B_X(3,271)=5.000000000000000E-002_DP
 S_B_FROM_V%A_Y(3,271)=-0.100000000000000E0_DP
 S_B_FROM_V%A_X(3,272)=-0.250000000000000E0_DP
 S_B_FROM_V%B_Y(3,272)=-8.333333333333334E-002_DP
 S_B_FROM_V%B_X(3,273)=-0.333333333333333E0_DP
 S_B_FROM_V%A_Y(3,273)=0.666666666666667E0_DP
 S_B_FROM_V%A_X(3,274)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,274)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,104)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,118)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,118)=-1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,119)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(4,133)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,133)=7.50000000000000E0_DP
 S_B_FROM_V%A_X(4,134)=4.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,134)=1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,135)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,149)=13.1250000000000E0_DP
 S_B_FROM_V%B_Y(4,149)=3.87500000000000E0_DP
 S_B_FROM_V%B_X(4,150)=3.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,150)=-6.16666666666667E0_DP
 S_B_FROM_V%A_X(4,151)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,151)=-1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,152)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(4,166)=4.65000000000000E0_DP
 S_B_FROM_V%A_Y(4,166)=-11.0250000000000E0_DP
 S_B_FROM_V%A_X(4,167)=-9.25000000000001E0_DP
 S_B_FROM_V%B_Y(4,167)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(4,168)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,168)=5.00000000000000E0_DP
 S_B_FROM_V%A_X(4,169)=3.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,169)=1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,170)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,184)=-9.18749999999999E0_DP
 S_B_FROM_V%B_Y(4,184)=-2.81250000000000E0_DP
 S_B_FROM_V%B_X(4,185)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,185)=6.37500000000001E0_DP
 S_B_FROM_V%A_X(4,186)=6.25000000000000E0_DP
 S_B_FROM_V%B_Y(4,186)=2.25000000000000E0_DP
 S_B_FROM_V%B_X(4,187)=2.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,187)=-4.00000000000000E0_DP
 S_B_FROM_V%A_X(4,188)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,188)=-1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,189)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(4,203)=-1.60714285714286E0_DP
 S_B_FROM_V%A_Y(4,203)=4.31250000000000E0_DP
 S_B_FROM_V%A_X(4,204)=4.25000000000000E0_DP
 S_B_FROM_V%B_Y(4,204)=1.43750000000000E0_DP
 S_B_FROM_V%B_X(4,205)=1.80000000000000E0_DP
 S_B_FROM_V%A_Y(4,205)=-3.37500000000000E0_DP
 S_B_FROM_V%A_X(4,206)=-4.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,206)=-1.62500000000000E0_DP
 S_B_FROM_V%B_X(4,207)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,207)=3.16666666666667E0_DP
 S_B_FROM_V%A_X(4,208)=2.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,208)=1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,209)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,210)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,210)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,223)=1.61718750000000E0_DP
 S_B_FROM_V%B_Y(4,223)=0.498883928571429E0_DP
 S_B_FROM_V%B_X(4,224)=0.616071428571429E0_DP
 S_B_FROM_V%A_Y(4,224)=-1.49107142857143E0_DP
 S_B_FROM_V%A_X(4,225)=-1.68750000000000E0_DP
 S_B_FROM_V%B_Y(4,225)=-0.637500000000000E0_DP
 S_B_FROM_V%B_X(4,226)=-0.975000000000000E0_DP
 S_B_FROM_V%A_Y(4,226)=1.57500000000000E0_DP
 S_B_FROM_V%A_X(4,227)=2.37500000000000E0_DP
 S_B_FROM_V%B_Y(4,227)=1.12500000000000E0_DP
 S_B_FROM_V%B_X(4,228)=1.50000000000000E0_DP
 S_B_FROM_V%A_Y(4,228)=-2.50000000000000E0_DP
 S_B_FROM_V%A_X(4,229)=-1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,229)=-1.50000000000000E0_DP
 S_B_FROM_V%B_X(4,230)=3.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,230)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(4,244)=0.110863095238095E0_DP
 S_B_FROM_V%A_Y(4,244)=-0.312500000000000E0_DP
 S_B_FROM_V%A_X(4,245)=-0.372767857142858E0_DP
 S_B_FROM_V%B_Y(4,245)=-0.127232142857143E0_DP
 S_B_FROM_V%B_X(4,246)=-0.182142857142857E0_DP
 S_B_FROM_V%A_Y(4,246)=0.392857142857143E0_DP
 S_B_FROM_V%A_X(4,247)=0.525000000000000E0_DP
 S_B_FROM_V%B_Y(4,247)=0.225000000000000E0_DP
 S_B_FROM_V%B_X(4,248)=0.450000000000000E0_DP
 S_B_FROM_V%A_Y(4,248)=-0.600000000000000E0_DP
 S_B_FROM_V%A_X(4,249)=-1.25000000000000E0_DP
 S_B_FROM_V%B_Y(4,249)=-0.750000000000000E0_DP
 S_B_FROM_V%B_X(4,250)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,250)=2.00000000000000E0_DP
 S_B_FROM_V%A_X(4,251)=-3.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,251)=-3.00000000000000E0_DP
 S_B_FROM_V%A_X(4,266)=-3.124999999999997E-002_DP
 S_B_FROM_V%B_Y(4,266)=-9.672619047619055E-003_DP
 S_B_FROM_V%B_X(4,267)=-1.413690476190477E-002_DP
 S_B_FROM_V%A_Y(4,267)=3.596230158730165E-002_DP
 S_B_FROM_V%A_X(4,268)=4.910714285714287E-002_DP
 S_B_FROM_V%B_Y(4,268)=1.874999999999999E-002_DP
 S_B_FROM_V%B_X(4,269)=3.214285714285715E-002_DP
 S_B_FROM_V%A_Y(4,269)=-6.071428571428574E-002_DP
 S_B_FROM_V%A_X(4,270)=-0.100000000000000E0_DP
 S_B_FROM_V%B_Y(4,270)=-4.999999999999999E-002_DP
 S_B_FROM_V%B_X(4,271)=-0.150000000000000E0_DP
 S_B_FROM_V%A_Y(4,271)=0.150000000000000E0_DP
 S_B_FROM_V%A_X(4,272)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(4,272)=0.500000000000000E0_DP
 S_B_FROM_V%B_X(4,273)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,273)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,104)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,118)=4.50000000000000E0_DP
 S_B_FROM_V%B_Y(5,118)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,119)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,133)=5.33333333333333E0_DP
 S_B_FROM_V%A_Y(5,133)=-9.00000000000000E0_DP
 S_B_FROM_V%A_X(5,134)=-4.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,134)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,135)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,149)=-15.7500000000000E0_DP
 S_B_FROM_V%B_Y(5,149)=-6.16666666666667E0_DP
 S_B_FROM_V%B_X(5,150)=-4.66666666666667E0_DP
 S_B_FROM_V%A_Y(5,150)=7.66666666666667E0_DP
 S_B_FROM_V%A_X(5,151)=3.50000000000000E0_DP
 S_B_FROM_V%B_Y(5,151)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,152)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,166)=-7.40000000000000E0_DP
 S_B_FROM_V%A_Y(5,166)=12.9750000000000E0_DP
 S_B_FROM_V%A_X(5,167)=11.5000000000000E0_DP
 S_B_FROM_V%B_Y(5,167)=5.00000000000000E0_DP
 S_B_FROM_V%B_X(5,168)=4.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,168)=-6.50000000000000E0_DP
 S_B_FROM_V%A_X(5,169)=-3.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,169)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,170)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,184)=10.8125000000000E0_DP
 S_B_FROM_V%B_Y(5,184)=4.25000000000000E0_DP
 S_B_FROM_V%B_X(5,185)=5.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,185)=-7.87500000000000E0_DP
 S_B_FROM_V%A_X(5,186)=-8.12500000000000E0_DP
 S_B_FROM_V%B_Y(5,186)=-4.00000000000000E0_DP
 S_B_FROM_V%B_X(5,187)=-3.33333333333333E0_DP
 S_B_FROM_V%A_Y(5,187)=5.50000000000000E0_DP
 S_B_FROM_V%A_X(5,188)=2.50000000000000E0_DP
 S_B_FROM_V%B_Y(5,188)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,189)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,190)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,190)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,203)=2.42857142857143E0_DP
 S_B_FROM_V%A_Y(5,203)=-4.96428571428572E0_DP
 S_B_FROM_V%A_X(5,204)=-5.25000000000000E0_DP
 S_B_FROM_V%B_Y(5,204)=-2.25000000000000E0_DP
 S_B_FROM_V%B_X(5,205)=-3.20000000000000E0_DP
 S_B_FROM_V%A_Y(5,205)=4.50000000000000E0_DP
 S_B_FROM_V%A_X(5,206)=5.50000000000000E0_DP
 S_B_FROM_V%B_Y(5,206)=3.16666666666667E0_DP
 S_B_FROM_V%B_X(5,207)=2.66666666666667E0_DP
 S_B_FROM_V%A_Y(5,207)=-4.66666666666667E0_DP
 S_B_FROM_V%A_X(5,208)=-2.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,208)=-2.00000000000000E0_DP
 S_B_FROM_V%B_X(5,209)=4.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,209)=-4.00000000000000E0_DP
 S_B_FROM_V%A_X(5,223)=-1.86160714285714E0_DP
 S_B_FROM_V%B_Y(5,223)=-0.745535714285715E0_DP
 S_B_FROM_V%B_X(5,224)=-0.964285714285715E0_DP
 S_B_FROM_V%A_Y(5,224)=1.78571428571429E0_DP
 S_B_FROM_V%A_X(5,225)=2.25000000000000E0_DP
 S_B_FROM_V%B_Y(5,225)=1.05000000000000E0_DP
 S_B_FROM_V%B_X(5,226)=1.90000000000000E0_DP
 S_B_FROM_V%A_Y(5,226)=-2.40000000000000E0_DP
 S_B_FROM_V%A_X(5,227)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(5,227)=-2.50000000000000E0_DP
 S_B_FROM_V%B_X(5,228)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,228)=4.00000000000000E0_DP
 S_B_FROM_V%A_X(5,229)=-6.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,229)=-6.00000000000000E0_DP
 S_B_FROM_V%B_X(5,244)=-0.165674603174603E0_DP
 S_B_FROM_V%A_Y(5,244)=0.357142857142857E0_DP
 S_B_FROM_V%A_X(5,245)=0.446428571428572E0_DP
 S_B_FROM_V%B_Y(5,245)=0.196428571428571E0_DP
 S_B_FROM_V%B_X(5,246)=0.300000000000000E0_DP
 S_B_FROM_V%A_Y(5,246)=-0.500000000000000E0_DP
 S_B_FROM_V%A_X(5,247)=-0.800000000000000E0_DP
 S_B_FROM_V%B_Y(5,247)=-0.400000000000000E0_DP
 S_B_FROM_V%B_X(5,248)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,248)=1.20000000000000E0_DP
 S_B_FROM_V%A_X(5,249)=2.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,249)=2.00000000000000E0_DP
 S_B_FROM_V%B_X(5,250)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,250)=4.00000000000000E0_DP
 S_B_FROM_V%A_X(5,266)=3.571428571428573E-002_DP
 S_B_FROM_V%B_Y(5,266)=1.438492063492064E-002_DP
 S_B_FROM_V%B_X(5,267)=2.182539682539683E-002_DP
 S_B_FROM_V%A_Y(5,267)=-4.265873015873016E-002_DP
 S_B_FROM_V%A_X(5,268)=-6.249999999999996E-002_DP
 S_B_FROM_V%B_Y(5,268)=-3.035714285714283E-002_DP
 S_B_FROM_V%B_X(5,269)=-5.714285714285713E-002_DP
 S_B_FROM_V%A_Y(5,269)=8.571428571428572E-002_DP
 S_B_FROM_V%A_X(5,270)=0.200000000000000E0_DP
 S_B_FROM_V%B_Y(5,270)=0.100000000000000E0_DP
 S_B_FROM_V%B_X(5,271)=0.400000000000000E0_DP
 S_B_FROM_V%A_Y(5,271)=-0.600000000000000E0_DP
 S_B_FROM_V%A_X(5,272)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,272)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,104)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,118)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,118)=-2.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,119)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(6,133)=-6.66666666666667E0_DP
 S_B_FROM_V%A_Y(6,133)=11.0000000000000E0_DP
 S_B_FROM_V%A_X(6,134)=4.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,134)=2.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,135)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,149)=19.2500000000000E0_DP
 S_B_FROM_V%B_Y(6,149)=9.58333333333333E0_DP
 S_B_FROM_V%B_X(6,150)=5.83333333333333E0_DP
 S_B_FROM_V%A_Y(6,150)=-9.66666666666667E0_DP
 S_B_FROM_V%A_X(6,151)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,151)=-2.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,152)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(6,166)=11.5000000000000E0_DP
 S_B_FROM_V%A_Y(6,166)=-16.5750000000000E0_DP
 S_B_FROM_V%A_X(6,167)=-14.5000000000000E0_DP
 S_B_FROM_V%B_Y(6,167)=-8.12500000000000E0_DP
 S_B_FROM_V%B_X(6,168)=-5.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,168)=8.50000000000000E0_DP
 S_B_FROM_V%A_X(6,169)=3.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,169)=2.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,170)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,171)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,171)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,184)=-13.8125000000000E0_DP
 S_B_FROM_V%B_Y(6,184)=-6.56250000000000E0_DP
 S_B_FROM_V%B_X(6,185)=-8.12500000000000E0_DP
 S_B_FROM_V%A_Y(6,185)=10.8750000000000E0_DP
 S_B_FROM_V%A_X(6,186)=10.6250000000000E0_DP
 S_B_FROM_V%B_Y(6,186)=6.87500000000000E0_DP
 S_B_FROM_V%B_X(6,187)=4.16666666666667E0_DP
 S_B_FROM_V%A_Y(6,187)=-7.50000000000000E0_DP
 S_B_FROM_V%A_X(6,188)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,188)=-2.50000000000000E0_DP
 S_B_FROM_V%B_X(6,189)=5.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,189)=-5.00000000000000E0_DP
 S_B_FROM_V%B_X(6,203)=-3.75000000000000E0_DP
 S_B_FROM_V%A_Y(6,203)=6.07142857142857E0_DP
 S_B_FROM_V%A_X(6,204)=7.25000000000000E0_DP
 S_B_FROM_V%B_Y(6,204)=3.75000000000000E0_DP
 S_B_FROM_V%B_X(6,205)=5.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,205)=-7.00000000000000E0_DP
 S_B_FROM_V%A_X(6,206)=-7.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,206)=-5.83333333333333E0_DP
 S_B_FROM_V%B_X(6,207)=-3.33333333333333E0_DP
 S_B_FROM_V%A_Y(6,207)=6.66666666666667E0_DP
 S_B_FROM_V%A_X(6,208)=-10.0000000000000E0_DP
 S_B_FROM_V%B_Y(6,208)=-10.0000000000000E0_DP
 S_B_FROM_V%A_X(6,223)=2.27678571428571E0_DP
 S_B_FROM_V%B_Y(6,223)=1.11607142857143E0_DP
 S_B_FROM_V%B_X(6,224)=1.60714285714286E0_DP
 S_B_FROM_V%A_Y(6,224)=-2.32142857142857E0_DP
 S_B_FROM_V%A_X(6,225)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,225)=-2.00000000000000E0_DP
 S_B_FROM_V%B_X(6,226)=-3.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,226)=4.50000000000000E0_DP
 S_B_FROM_V%A_X(6,227)=5.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,227)=5.00000000000000E0_DP
 S_B_FROM_V%B_X(6,228)=-10.0000000000000E0_DP
 S_B_FROM_V%A_Y(6,228)=10.0000000000000E0_DP
 S_B_FROM_V%B_X(6,244)=0.248015873015873E0_DP
 S_B_FROM_V%A_Y(6,244)=-0.431547619047619E0_DP
 S_B_FROM_V%A_X(6,245)=-0.580357142857143E0_DP
 S_B_FROM_V%B_Y(6,245)=-0.312500000000000E0_DP
 S_B_FROM_V%B_X(6,246)=-0.571428571428571E0_DP
 S_B_FROM_V%A_Y(6,246)=0.714285714285714E0_DP
 S_B_FROM_V%A_X(6,247)=1.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,247)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(6,248)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,248)=-3.00000000000000E0_DP
 S_B_FROM_V%A_X(6,249)=5.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,249)=5.00000000000000E0_DP
 S_B_FROM_V%A_X(6,266)=-4.315476190476190E-002_DP
 S_B_FROM_V%B_Y(6,266)=-2.132936507936508E-002_DP
 S_B_FROM_V%B_X(6,267)=-3.472222222222222E-002_DP
 S_B_FROM_V%A_Y(6,267)=5.456349206349211E-002_DP
 S_B_FROM_V%A_X(6,268)=8.928571428571426E-002_DP
 S_B_FROM_V%B_Y(6,268)=5.357142857142858E-002_DP
 S_B_FROM_V%B_X(6,269)=0.142857142857143E0_DP
 S_B_FROM_V%A_Y(6,269)=-0.142857142857143E0_DP
 S_B_FROM_V%A_X(6,270)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(6,270)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(6,271)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,271)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,104)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(7,118)=4.50000000000000E0_DP
 S_B_FROM_V%B_Y(7,118)=3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,119)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(7,133)=8.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,133)=-13.5000000000000E0_DP
 S_B_FROM_V%A_X(7,134)=-4.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,134)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,135)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(7,149)=-23.6250000000000E0_DP
 S_B_FROM_V%B_Y(7,149)=-14.5000000000000E0_DP
 S_B_FROM_V%B_X(7,150)=-7.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,150)=12.1666666666667E0_DP
 S_B_FROM_V%A_X(7,151)=3.50000000000000E0_DP
 S_B_FROM_V%B_Y(7,151)=3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,152)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(7,153)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,153)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(7,166)=-17.4000000000000E0_DP
 S_B_FROM_V%A_Y(7,166)=22.9500000000000E0_DP
 S_B_FROM_V%A_X(7,167)=18.2500000000000E0_DP
 S_B_FROM_V%B_Y(7,167)=12.7500000000000E0_DP
 S_B_FROM_V%B_X(7,168)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,168)=-11.0000000000000E0_DP
 S_B_FROM_V%A_X(7,169)=-3.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,169)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(7,170)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,170)=-6.00000000000000E0_DP
 S_B_FROM_V%A_X(7,184)=19.1250000000000E0_DP
 S_B_FROM_V%B_Y(7,184)=10.8750000000000E0_DP
 S_B_FROM_V%B_X(7,185)=12.7500000000000E0_DP
 S_B_FROM_V%A_Y(7,185)=-16.5000000000000E0_DP
 S_B_FROM_V%A_X(7,186)=-13.7500000000000E0_DP
 S_B_FROM_V%B_Y(7,186)=-11.2500000000000E0_DP
 S_B_FROM_V%B_X(7,187)=-5.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,187)=10.0000000000000E0_DP
 S_B_FROM_V%A_X(7,188)=-15.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,188)=-15.0000000000000E0_DP
 S_B_FROM_V%B_X(7,203)=6.21428571428572E0_DP
 S_B_FROM_V%A_Y(7,203)=-8.21428571428572E0_DP
 S_B_FROM_V%A_X(7,204)=-11.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,204)=-7.00000000000000E0_DP
 S_B_FROM_V%B_X(7,205)=-9.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,205)=12.0000000000000E0_DP
 S_B_FROM_V%A_X(7,206)=10.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,206)=10.0000000000000E0_DP
 S_B_FROM_V%B_X(7,207)=-20.0000000000000E0_DP
 S_B_FROM_V%A_Y(7,207)=20.0000000000000E0_DP
 S_B_FROM_V%A_X(7,223)=-3.08035714285714E0_DP
 S_B_FROM_V%B_Y(7,223)=-1.74107142857143E0_DP
 S_B_FROM_V%B_X(7,224)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,224)=3.57142857142857E0_DP
 S_B_FROM_V%A_X(7,225)=6.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,225)=4.50000000000000E0_DP
 S_B_FROM_V%B_X(7,226)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,226)=-9.00000000000000E0_DP
 S_B_FROM_V%A_X(7,227)=15.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,227)=15.0000000000000E0_DP
 S_B_FROM_V%B_X(7,244)=-0.386904761904762E0_DP
 S_B_FROM_V%A_Y(7,244)=0.565476190476191E0_DP
 S_B_FROM_V%A_X(7,245)=0.892857142857142E0_DP
 S_B_FROM_V%B_Y(7,245)=0.535714285714286E0_DP
 S_B_FROM_V%B_X(7,246)=1.28571428571429E0_DP
 S_B_FROM_V%A_Y(7,246)=-1.42857142857143E0_DP
 S_B_FROM_V%A_X(7,247)=-3.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,247)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(7,248)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,248)=-6.00000000000000E0_DP
 S_B_FROM_V%A_X(7,266)=5.654761904761905E-002_DP
 S_B_FROM_V%B_Y(7,266)=3.273809523809525E-002_DP
 S_B_FROM_V%B_X(7,267)=5.952380952380951E-002_DP
 S_B_FROM_V%A_Y(7,267)=-7.936507936507933E-002_DP
 S_B_FROM_V%A_X(7,268)=-0.178571428571429E0_DP
 S_B_FROM_V%B_Y(7,268)=-0.107142857142857E0_DP
 S_B_FROM_V%B_X(7,269)=-0.428571428571429E0_DP
 S_B_FROM_V%A_Y(7,269)=0.571428571428571E0_DP
 S_B_FROM_V%A_X(7,270)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,270)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,104)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,118)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(8,118)=-3.50000000000000E0_DP
 S_B_FROM_V%A_Y(8,119)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(8,133)=-9.33333333333334E0_DP
 S_B_FROM_V%A_Y(8,133)=16.5000000000000E0_DP
 S_B_FROM_V%A_X(8,134)=4.00000000000000E0_DP
 S_B_FROM_V%B_Y(8,134)=3.50000000000000E0_DP
 S_B_FROM_V%A_Y(8,135)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,136)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(8,136)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,149)=28.8750000000000E0_DP
 S_B_FROM_V%B_Y(8,149)=21.2916666666667E0_DP
 S_B_FROM_V%B_X(8,150)=8.16666666666667E0_DP
 S_B_FROM_V%A_Y(8,150)=-15.1666666666667E0_DP
 S_B_FROM_V%A_X(8,151)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(8,151)=-3.50000000000000E0_DP
 S_B_FROM_V%B_X(8,152)=7.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,152)=-7.00000000000000E0_DP
 S_B_FROM_V%B_X(8,166)=25.5500000000000E0_DP
 S_B_FROM_V%A_Y(8,166)=-33.6000000000000E0_DP
 S_B_FROM_V%A_X(8,167)=-22.7500000000000E0_DP
 S_B_FROM_V%B_Y(8,167)=-19.2500000000000E0_DP
 S_B_FROM_V%B_X(8,168)=-7.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,168)=14.0000000000000E0_DP
 S_B_FROM_V%A_X(8,169)=-21.0000000000000E0_DP
 S_B_FROM_V%B_Y(8,169)=-21.0000000000000E0_DP
 S_B_FROM_V%A_X(8,184)=-28.0000000000000E0_DP
 S_B_FROM_V%B_Y(8,184)=-19.2500000000000E0_DP
 S_B_FROM_V%B_X(8,185)=-19.2500000000000E0_DP
 S_B_FROM_V%A_Y(8,185)=26.2500000000000E0_DP
 S_B_FROM_V%A_X(8,186)=17.5000000000000E0_DP
 S_B_FROM_V%B_Y(8,186)=17.5000000000000E0_DP
 S_B_FROM_V%B_X(8,187)=-35.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,187)=35.0000000000000E0_DP
 S_B_FROM_V%B_X(8,203)=-11.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,203)=13.0000000000000E0_DP
 S_B_FROM_V%A_X(8,204)=17.5000000000000E0_DP
 S_B_FROM_V%B_Y(8,204)=14.0000000000000E0_DP
 S_B_FROM_V%B_X(8,205)=14.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,205)=-21.0000000000000E0_DP
 S_B_FROM_V%A_X(8,206)=35.0000000000000E0_DP
 S_B_FROM_V%B_Y(8,206)=35.0000000000000E0_DP
 S_B_FROM_V%A_X(8,223)=4.87500000000000E0_DP
 S_B_FROM_V%B_Y(8,223)=3.12500000000000E0_DP
 S_B_FROM_V%B_X(8,224)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,224)=-7.00000000000000E0_DP
 S_B_FROM_V%A_X(8,225)=-10.5000000000000E0_DP
 S_B_FROM_V%B_Y(8,225)=-10.5000000000000E0_DP
 S_B_FROM_V%B_X(8,226)=21.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,226)=-21.0000000000000E0_DP
 S_B_FROM_V%B_X(8,244)=0.694444444444445E0_DP
 S_B_FROM_V%A_Y(8,244)=-0.833333333333333E0_DP
 S_B_FROM_V%A_X(8,245)=-1.75000000000000E0_DP
 S_B_FROM_V%B_Y(8,245)=-1.25000000000000E0_DP
 S_B_FROM_V%B_X(8,246)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,246)=4.00000000000000E0_DP
 S_B_FROM_V%A_X(8,247)=-7.00000000000000E0_DP
 S_B_FROM_V%B_Y(8,247)=-7.00000000000000E0_DP
 S_B_FROM_V%A_X(8,266)=-8.333333333333330E-002_DP
 S_B_FROM_V%B_Y(8,266)=-5.555555555555557E-002_DP
 S_B_FROM_V%B_X(8,267)=-0.138888888888889E0_DP
 S_B_FROM_V%A_Y(8,267)=0.138888888888889E0_DP
 S_B_FROM_V%A_X(8,268)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(8,268)=0.500000000000000E0_DP
 S_B_FROM_V%B_X(8,269)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,269)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,104)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(9,118)=4.50000000000000E0_DP
 S_B_FROM_V%B_Y(9,118)=4.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,119)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(9,120)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,120)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(9,133)=10.6666666666667E0_DP
 S_B_FROM_V%A_Y(9,133)=-20.0000000000000E0_DP
 S_B_FROM_V%A_X(9,134)=-4.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,134)=-4.00000000000000E0_DP
 S_B_FROM_V%B_X(9,135)=8.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,135)=-8.00000000000000E0_DP
 S_B_FROM_V%A_X(9,149)=-35.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,149)=-30.3333333333333E0_DP
 S_B_FROM_V%B_X(9,150)=-9.33333333333333E0_DP
 S_B_FROM_V%A_Y(9,150)=18.6666666666667E0_DP
 S_B_FROM_V%A_X(9,151)=-28.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,151)=-28.0000000000000E0_DP
 S_B_FROM_V%B_X(9,166)=-36.4000000000000E0_DP
 S_B_FROM_V%A_Y(9,166)=50.4000000000000E0_DP
 S_B_FROM_V%A_X(9,167)=28.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,167)=28.0000000000000E0_DP
 S_B_FROM_V%B_X(9,168)=-56.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,168)=56.0000000000000E0_DP
 S_B_FROM_V%A_X(9,184)=42.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,184)=35.0000000000000E0_DP
 S_B_FROM_V%B_X(9,185)=28.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,185)=-42.0000000000000E0_DP
 S_B_FROM_V%A_X(9,186)=70.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,186)=70.0000000000000E0_DP
 S_B_FROM_V%B_X(9,203)=20.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,203)=-24.0000000000000E0_DP
 S_B_FROM_V%A_X(9,204)=-28.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,204)=-28.0000000000000E0_DP
 S_B_FROM_V%B_X(9,205)=56.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,205)=-56.0000000000000E0_DP
 S_B_FROM_V%A_X(9,223)=-9.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,223)=-7.00000000000000E0_DP
 S_B_FROM_V%B_X(9,224)=-12.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,224)=16.0000000000000E0_DP
 S_B_FROM_V%A_X(9,225)=-28.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,225)=-28.0000000000000E0_DP
 S_B_FROM_V%B_X(9,244)=-1.55555555555556E0_DP
 S_B_FROM_V%A_Y(9,244)=1.66666666666667E0_DP
 S_B_FROM_V%A_X(9,245)=4.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,245)=4.00000000000000E0_DP
 S_B_FROM_V%B_X(9,246)=-8.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,246)=8.00000000000000E0_DP
 S_B_FROM_V%A_X(9,266)=0.166666666666667E0_DP
 S_B_FROM_V%B_Y(9,266)=0.111111111111111E0_DP
 S_B_FROM_V%B_X(9,267)=0.444444444444444E0_DP
 S_B_FROM_V%A_Y(9,267)=-0.555555555555556E0_DP
 S_B_FROM_V%A_X(9,268)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,268)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,104)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(10,105)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(10,105)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(10,118)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(10,118)=-4.50000000000000E0_DP
 S_B_FROM_V%B_X(10,119)=9.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,119)=-9.00000000000000E0_DP
 S_B_FROM_V%B_X(10,133)=-12.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,133)=24.0000000000000E0_DP
 S_B_FROM_V%A_X(10,134)=-36.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,134)=-36.0000000000000E0_DP
 S_B_FROM_V%A_X(10,149)=42.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,149)=42.0000000000000E0_DP
 S_B_FROM_V%B_X(10,150)=-84.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,150)=84.0000000000000E0_DP
 S_B_FROM_V%B_X(10,166)=50.4000000000000E0_DP
 S_B_FROM_V%A_Y(10,166)=-75.6000000000000E0_DP
 S_B_FROM_V%A_X(10,167)=126.000000000000E0_DP
 S_B_FROM_V%B_Y(10,167)=126.000000000000E0_DP
 S_B_FROM_V%A_X(10,184)=-63.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,184)=-63.0000000000000E0_DP
 S_B_FROM_V%B_X(10,185)=126.000000000000E0_DP
 S_B_FROM_V%A_Y(10,185)=-126.000000000000E0_DP
 S_B_FROM_V%B_X(10,203)=-36.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,203)=48.0000000000000E0_DP
 S_B_FROM_V%A_X(10,204)=-84.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,204)=-84.0000000000000E0_DP
 S_B_FROM_V%A_X(10,223)=18.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,223)=18.0000000000000E0_DP
 S_B_FROM_V%B_X(10,224)=-36.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,224)=36.0000000000000E0_DP
 S_B_FROM_V%B_X(10,244)=4.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,244)=-5.00000000000000E0_DP
 S_B_FROM_V%A_X(10,245)=9.00000000000000E0_DP
 S_B_FROM_V%B_Y(10,245)=9.00000000000000E0_DP
 S_B_FROM_V%A_X(10,266)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(10,266)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(10,267)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,267)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(11,91)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(11,91)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(11,104)=10.0000000000000E0_DP
 S_B_FROM_V%A_Y(11,104)=-10.0000000000000E0_DP
 S_B_FROM_V%A_X(11,118)=-45.0000000000000E0_DP
 S_B_FROM_V%B_Y(11,118)=-45.0000000000000E0_DP
 S_B_FROM_V%B_X(11,133)=-120.000000000000E0_DP
 S_B_FROM_V%A_Y(11,133)=120.000000000000E0_DP
 S_B_FROM_V%A_X(11,149)=210.000000000000E0_DP
 S_B_FROM_V%B_Y(11,149)=210.000000000000E0_DP
 S_B_FROM_V%B_X(11,166)=252.000000000000E0_DP
 S_B_FROM_V%A_Y(11,166)=-252.000000000000E0_DP
 S_B_FROM_V%A_X(11,184)=-210.000000000000E0_DP
 S_B_FROM_V%B_Y(11,184)=-210.000000000000E0_DP
 S_B_FROM_V%B_X(11,203)=-120.000000000000E0_DP
 S_B_FROM_V%A_Y(11,203)=120.000000000000E0_DP
 S_B_FROM_V%A_X(11,223)=45.0000000000000E0_DP
 S_B_FROM_V%B_Y(11,223)=45.0000000000000E0_DP
 S_B_FROM_V%B_X(11,244)=10.0000000000000E0_DP
 S_B_FROM_V%A_Y(11,244)=-10.0000000000000E0_DP
 S_B_FROM_V%A_X(11,266)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(11,266)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(12,78)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(12,78)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(12,90)=11.0000000000000E0_DP
 S_B_FROM_V%A_Y(12,90)=-11.0000000000000E0_DP
 S_B_FROM_V%A_X(12,103)=-55.0000000000000E0_DP
 S_B_FROM_V%B_Y(12,103)=-55.0000000000000E0_DP
 S_B_FROM_V%B_X(12,117)=-165.000000000000E0_DP
 S_B_FROM_V%A_Y(12,117)=165.000000000000E0_DP
 S_B_FROM_V%A_X(12,132)=330.000000000000E0_DP
 S_B_FROM_V%B_Y(12,132)=330.000000000000E0_DP
 S_B_FROM_V%B_X(12,148)=462.000000000000E0_DP
 S_B_FROM_V%A_Y(12,148)=-462.000000000000E0_DP
 S_B_FROM_V%A_X(12,165)=-462.000000000000E0_DP
 S_B_FROM_V%B_Y(12,165)=-462.000000000000E0_DP
 S_B_FROM_V%B_X(12,183)=-330.000000000000E0_DP
 S_B_FROM_V%A_Y(12,183)=330.000000000000E0_DP
 S_B_FROM_V%A_X(12,202)=165.000000000000E0_DP
 S_B_FROM_V%B_Y(12,202)=165.000000000000E0_DP
 S_B_FROM_V%B_X(12,222)=55.0000000000000E0_DP
 S_B_FROM_V%A_Y(12,222)=-55.0000000000000E0_DP
 S_B_FROM_V%A_X(12,243)=-11.0000000000000E0_DP
 S_B_FROM_V%B_Y(12,243)=-11.0000000000000E0_DP
 S_B_FROM_V%B_X(12,265)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(12,265)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(13,66)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(13,66)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(13,77)=12.0000000000000E0_DP
 S_B_FROM_V%A_Y(13,77)=-12.0000000000000E0_DP
 S_B_FROM_V%A_X(13,89)=-66.0000000000000E0_DP
 S_B_FROM_V%B_Y(13,89)=-66.0000000000000E0_DP
 S_B_FROM_V%B_X(13,102)=-220.000000000000E0_DP
 S_B_FROM_V%A_Y(13,102)=220.000000000000E0_DP
 S_B_FROM_V%A_X(13,116)=495.000000000000E0_DP
 S_B_FROM_V%B_Y(13,116)=495.000000000000E0_DP
 S_B_FROM_V%B_X(13,131)=792.000000000000E0_DP
 S_B_FROM_V%A_Y(13,131)=-792.000000000000E0_DP
 S_B_FROM_V%A_X(13,147)=-924.000000000000E0_DP
 S_B_FROM_V%B_Y(13,147)=-924.000000000000E0_DP
 S_B_FROM_V%B_X(13,164)=-792.000000000000E0_DP
 S_B_FROM_V%A_Y(13,164)=792.000000000000E0_DP
 S_B_FROM_V%A_X(13,182)=495.000000000000E0_DP
 S_B_FROM_V%B_Y(13,182)=495.000000000000E0_DP
 S_B_FROM_V%B_X(13,201)=220.000000000000E0_DP
 S_B_FROM_V%A_Y(13,201)=-220.000000000000E0_DP
 S_B_FROM_V%A_X(13,221)=-66.0000000000000E0_DP
 S_B_FROM_V%B_Y(13,221)=-66.0000000000000E0_DP
 S_B_FROM_V%B_X(13,242)=-12.0000000000000E0_DP
 S_B_FROM_V%A_Y(13,242)=12.0000000000000E0_DP
 S_B_FROM_V%A_X(13,264)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(13,264)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(14,55)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(14,55)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(14,65)=13.0000000000000E0_DP
 S_B_FROM_V%A_Y(14,65)=-13.0000000000000E0_DP
 S_B_FROM_V%A_X(14,76)=-78.0000000000000E0_DP
 S_B_FROM_V%B_Y(14,76)=-78.0000000000000E0_DP
 S_B_FROM_V%B_X(14,88)=-286.000000000000E0_DP
 S_B_FROM_V%A_Y(14,88)=286.000000000000E0_DP
 S_B_FROM_V%A_X(14,101)=715.000000000000E0_DP
 S_B_FROM_V%B_Y(14,101)=715.000000000000E0_DP
 S_B_FROM_V%B_X(14,115)=1287.00000000000E0_DP
 S_B_FROM_V%A_Y(14,115)=-1287.00000000000E0_DP
 S_B_FROM_V%A_X(14,130)=-1716.00000000000E0_DP
 S_B_FROM_V%B_Y(14,130)=-1716.00000000000E0_DP
 S_B_FROM_V%B_X(14,146)=-1716.00000000000E0_DP
 S_B_FROM_V%A_Y(14,146)=1716.00000000000E0_DP
 S_B_FROM_V%A_X(14,163)=1287.00000000000E0_DP
 S_B_FROM_V%B_Y(14,163)=1287.00000000000E0_DP
 S_B_FROM_V%B_X(14,181)=715.000000000000E0_DP
 S_B_FROM_V%A_Y(14,181)=-715.000000000000E0_DP
 S_B_FROM_V%A_X(14,200)=-286.000000000000E0_DP
 S_B_FROM_V%B_Y(14,200)=-286.000000000000E0_DP
 S_B_FROM_V%B_X(14,220)=-78.0000000000000E0_DP
 S_B_FROM_V%A_Y(14,220)=78.0000000000000E0_DP
 S_B_FROM_V%A_X(14,241)=13.0000000000000E0_DP
 S_B_FROM_V%B_Y(14,241)=13.0000000000000E0_DP
 S_B_FROM_V%B_X(14,263)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(14,263)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(15,45)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(15,45)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(15,54)=14.0000000000000E0_DP
 S_B_FROM_V%A_Y(15,54)=-14.0000000000000E0_DP
 S_B_FROM_V%A_X(15,64)=-91.0000000000000E0_DP
 S_B_FROM_V%B_Y(15,64)=-91.0000000000000E0_DP
 S_B_FROM_V%B_X(15,75)=-364.000000000000E0_DP
 S_B_FROM_V%A_Y(15,75)=364.000000000000E0_DP
 S_B_FROM_V%A_X(15,87)=1001.00000000000E0_DP
 S_B_FROM_V%B_Y(15,87)=1001.00000000000E0_DP
 S_B_FROM_V%B_X(15,100)=2002.00000000000E0_DP
 S_B_FROM_V%A_Y(15,100)=-2002.00000000000E0_DP
 S_B_FROM_V%A_X(15,114)=-3003.00000000000E0_DP
 S_B_FROM_V%B_Y(15,114)=-3003.00000000000E0_DP
 S_B_FROM_V%B_X(15,129)=-3432.00000000000E0_DP
 S_B_FROM_V%A_Y(15,129)=3432.00000000000E0_DP
 S_B_FROM_V%A_X(15,145)=3003.00000000000E0_DP
 S_B_FROM_V%B_Y(15,145)=3003.00000000000E0_DP
 S_B_FROM_V%B_X(15,162)=2002.00000000000E0_DP
 S_B_FROM_V%A_Y(15,162)=-2002.00000000000E0_DP
 S_B_FROM_V%A_X(15,180)=-1001.00000000000E0_DP
 S_B_FROM_V%B_Y(15,180)=-1001.00000000000E0_DP
 S_B_FROM_V%B_X(15,199)=-364.000000000000E0_DP
 S_B_FROM_V%A_Y(15,199)=364.000000000000E0_DP
 S_B_FROM_V%A_X(15,219)=91.0000000000000E0_DP
 S_B_FROM_V%B_Y(15,219)=91.0000000000000E0_DP
 S_B_FROM_V%B_X(15,240)=14.0000000000000E0_DP
 S_B_FROM_V%A_Y(15,240)=-14.0000000000000E0_DP
 S_B_FROM_V%A_X(15,262)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(15,262)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(16,36)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(16,36)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(16,44)=15.0000000000000E0_DP
 S_B_FROM_V%A_Y(16,44)=-15.0000000000000E0_DP
 S_B_FROM_V%A_X(16,53)=-105.000000000000E0_DP
 S_B_FROM_V%B_Y(16,53)=-105.000000000000E0_DP
 S_B_FROM_V%B_X(16,63)=-455.000000000000E0_DP
 S_B_FROM_V%A_Y(16,63)=455.000000000000E0_DP
 S_B_FROM_V%A_X(16,74)=1365.00000000000E0_DP
 S_B_FROM_V%B_Y(16,74)=1365.00000000000E0_DP
 S_B_FROM_V%B_X(16,86)=3003.00000000000E0_DP
 S_B_FROM_V%A_Y(16,86)=-3003.00000000000E0_DP
 S_B_FROM_V%A_X(16,99)=-5005.00000000000E0_DP
 S_B_FROM_V%B_Y(16,99)=-5005.00000000000E0_DP
 S_B_FROM_V%B_X(16,113)=-6435.00000000000E0_DP
 S_B_FROM_V%A_Y(16,113)=6435.00000000000E0_DP
 S_B_FROM_V%A_X(16,128)=6435.00000000000E0_DP
 S_B_FROM_V%B_Y(16,128)=6435.00000000000E0_DP
 S_B_FROM_V%B_X(16,144)=5005.00000000000E0_DP
 S_B_FROM_V%A_Y(16,144)=-5005.00000000000E0_DP
 S_B_FROM_V%A_X(16,161)=-3003.00000000000E0_DP
 S_B_FROM_V%B_Y(16,161)=-3003.00000000000E0_DP
 S_B_FROM_V%B_X(16,179)=-1365.00000000000E0_DP
 S_B_FROM_V%A_Y(16,179)=1365.00000000000E0_DP
 S_B_FROM_V%A_X(16,198)=455.000000000000E0_DP
 S_B_FROM_V%B_Y(16,198)=455.000000000000E0_DP
 S_B_FROM_V%B_X(16,218)=105.000000000000E0_DP
 S_B_FROM_V%A_Y(16,218)=-105.000000000000E0_DP
 S_B_FROM_V%A_X(16,239)=-15.0000000000000E0_DP
 S_B_FROM_V%B_Y(16,239)=-15.0000000000000E0_DP
 S_B_FROM_V%B_X(16,261)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(16,261)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(17,28)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(17,28)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(17,35)=16.0000000000000E0_DP
 S_B_FROM_V%A_Y(17,35)=-16.0000000000000E0_DP
 S_B_FROM_V%A_X(17,43)=-120.000000000000E0_DP
 S_B_FROM_V%B_Y(17,43)=-120.000000000000E0_DP
 S_B_FROM_V%B_X(17,52)=-560.000000000000E0_DP
 S_B_FROM_V%A_Y(17,52)=560.000000000000E0_DP
 S_B_FROM_V%A_X(17,62)=1820.00000000000E0_DP
 S_B_FROM_V%B_Y(17,62)=1820.00000000000E0_DP
 S_B_FROM_V%B_X(17,73)=4368.00000000000E0_DP
 S_B_FROM_V%A_Y(17,73)=-4368.00000000000E0_DP
 S_B_FROM_V%A_X(17,85)=-8008.00000000000E0_DP
 S_B_FROM_V%B_Y(17,85)=-8008.00000000000E0_DP
 S_B_FROM_V%B_X(17,98)=-11440.0000000000E0_DP
 S_B_FROM_V%A_Y(17,98)=11440.0000000000E0_DP
 S_B_FROM_V%A_X(17,112)=12870.0000000000E0_DP
 S_B_FROM_V%B_Y(17,112)=12870.0000000000E0_DP
 S_B_FROM_V%B_X(17,127)=11440.0000000000E0_DP
 S_B_FROM_V%A_Y(17,127)=-11440.0000000000E0_DP
 S_B_FROM_V%A_X(17,143)=-8008.00000000000E0_DP
 S_B_FROM_V%B_Y(17,143)=-8008.00000000000E0_DP
 S_B_FROM_V%B_X(17,160)=-4368.00000000000E0_DP
 S_B_FROM_V%A_Y(17,160)=4368.00000000000E0_DP
 S_B_FROM_V%A_X(17,178)=1820.00000000000E0_DP
 S_B_FROM_V%B_Y(17,178)=1820.00000000000E0_DP
 S_B_FROM_V%B_X(17,197)=560.000000000000E0_DP
 S_B_FROM_V%A_Y(17,197)=-560.000000000000E0_DP
 S_B_FROM_V%A_X(17,217)=-120.000000000000E0_DP
 S_B_FROM_V%B_Y(17,217)=-120.000000000000E0_DP
 S_B_FROM_V%B_X(17,238)=-16.0000000000000E0_DP
 S_B_FROM_V%A_Y(17,238)=16.0000000000000E0_DP
 S_B_FROM_V%A_X(17,260)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(17,260)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(18,21)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(18,21)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(18,27)=17.0000000000000E0_DP
 S_B_FROM_V%A_Y(18,27)=-17.0000000000000E0_DP
 S_B_FROM_V%A_X(18,34)=-136.000000000000E0_DP
 S_B_FROM_V%B_Y(18,34)=-136.000000000000E0_DP
 S_B_FROM_V%B_X(18,42)=-680.000000000000E0_DP
 S_B_FROM_V%A_Y(18,42)=680.000000000000E0_DP
 S_B_FROM_V%A_X(18,51)=2380.00000000000E0_DP
 S_B_FROM_V%B_Y(18,51)=2380.00000000000E0_DP
 S_B_FROM_V%B_X(18,61)=6188.00000000000E0_DP
 S_B_FROM_V%A_Y(18,61)=-6188.00000000000E0_DP
 S_B_FROM_V%A_X(18,72)=-12376.0000000000E0_DP
 S_B_FROM_V%B_Y(18,72)=-12376.0000000000E0_DP
 S_B_FROM_V%B_X(18,84)=-19448.0000000000E0_DP
 S_B_FROM_V%A_Y(18,84)=19448.0000000000E0_DP
 S_B_FROM_V%A_X(18,97)=24310.0000000000E0_DP
 S_B_FROM_V%B_Y(18,97)=24310.0000000000E0_DP
 S_B_FROM_V%B_X(18,111)=24310.0000000000E0_DP
 S_B_FROM_V%A_Y(18,111)=-24310.0000000000E0_DP
 S_B_FROM_V%A_X(18,126)=-19448.0000000000E0_DP
 S_B_FROM_V%B_Y(18,126)=-19448.0000000000E0_DP
 S_B_FROM_V%B_X(18,142)=-12376.0000000000E0_DP
 S_B_FROM_V%A_Y(18,142)=12376.0000000000E0_DP
 S_B_FROM_V%A_X(18,159)=6188.00000000000E0_DP
 S_B_FROM_V%B_Y(18,159)=6188.00000000000E0_DP
 S_B_FROM_V%B_X(18,177)=2380.00000000000E0_DP
 S_B_FROM_V%A_Y(18,177)=-2380.00000000000E0_DP
 S_B_FROM_V%A_X(18,196)=-680.000000000000E0_DP
 S_B_FROM_V%B_Y(18,196)=-680.000000000000E0_DP
 S_B_FROM_V%B_X(18,216)=-136.000000000000E0_DP
 S_B_FROM_V%A_Y(18,216)=136.000000000000E0_DP
 S_B_FROM_V%A_X(18,237)=17.0000000000000E0_DP
 S_B_FROM_V%B_Y(18,237)=17.0000000000000E0_DP
 S_B_FROM_V%B_X(18,259)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(18,259)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(19,15)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(19,15)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(19,20)=18.0000000000000E0_DP
 S_B_FROM_V%A_Y(19,20)=-18.0000000000000E0_DP
 S_B_FROM_V%A_X(19,26)=-153.000000000000E0_DP
 S_B_FROM_V%B_Y(19,26)=-153.000000000000E0_DP
 S_B_FROM_V%B_X(19,33)=-816.000000000000E0_DP
 S_B_FROM_V%A_Y(19,33)=816.000000000000E0_DP
 S_B_FROM_V%A_X(19,41)=3060.00000000000E0_DP
 S_B_FROM_V%B_Y(19,41)=3060.00000000000E0_DP
 S_B_FROM_V%B_X(19,50)=8568.00000000000E0_DP
 S_B_FROM_V%A_Y(19,50)=-8568.00000000000E0_DP
 S_B_FROM_V%A_X(19,60)=-18564.0000000000E0_DP
 S_B_FROM_V%B_Y(19,60)=-18564.0000000000E0_DP
 S_B_FROM_V%B_X(19,71)=-31824.0000000000E0_DP
 S_B_FROM_V%A_Y(19,71)=31824.0000000000E0_DP
 S_B_FROM_V%A_X(19,83)=43758.0000000000E0_DP
 S_B_FROM_V%B_Y(19,83)=43758.0000000000E0_DP
 S_B_FROM_V%B_X(19,96)=48620.0000000000E0_DP
 S_B_FROM_V%A_Y(19,96)=-48620.0000000000E0_DP
 S_B_FROM_V%A_X(19,110)=-43758.0000000000E0_DP
 S_B_FROM_V%B_Y(19,110)=-43758.0000000000E0_DP
 S_B_FROM_V%B_X(19,125)=-31824.0000000000E0_DP
 S_B_FROM_V%A_Y(19,125)=31824.0000000000E0_DP
 S_B_FROM_V%A_X(19,141)=18564.0000000000E0_DP
 S_B_FROM_V%B_Y(19,141)=18564.0000000000E0_DP
 S_B_FROM_V%B_X(19,158)=8568.00000000000E0_DP
 S_B_FROM_V%A_Y(19,158)=-8568.00000000000E0_DP
 S_B_FROM_V%A_X(19,176)=-3060.00000000000E0_DP
 S_B_FROM_V%B_Y(19,176)=-3060.00000000000E0_DP
 S_B_FROM_V%B_X(19,195)=-816.000000000000E0_DP
 S_B_FROM_V%A_Y(19,195)=816.000000000000E0_DP
 S_B_FROM_V%A_X(19,215)=153.000000000000E0_DP
 S_B_FROM_V%B_Y(19,215)=153.000000000000E0_DP
 S_B_FROM_V%B_X(19,236)=18.0000000000000E0_DP
 S_B_FROM_V%A_Y(19,236)=-18.0000000000000E0_DP
 S_B_FROM_V%A_X(19,258)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(19,258)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(20,10)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(20,10)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(20,14)=19.0000000000000E0_DP
 S_B_FROM_V%A_Y(20,14)=-19.0000000000000E0_DP
 S_B_FROM_V%A_X(20,19)=-171.000000000000E0_DP
 S_B_FROM_V%B_Y(20,19)=-171.000000000000E0_DP
 S_B_FROM_V%B_X(20,25)=-969.000000000000E0_DP
 S_B_FROM_V%A_Y(20,25)=969.000000000000E0_DP
 S_B_FROM_V%A_X(20,32)=3876.00000000000E0_DP
 S_B_FROM_V%B_Y(20,32)=3876.00000000000E0_DP
 S_B_FROM_V%B_X(20,40)=11628.0000000000E0_DP
 S_B_FROM_V%A_Y(20,40)=-11628.0000000000E0_DP
 S_B_FROM_V%A_X(20,49)=-27132.0000000000E0_DP
 S_B_FROM_V%B_Y(20,49)=-27132.0000000000E0_DP
 S_B_FROM_V%B_X(20,59)=-50388.0000000000E0_DP
 S_B_FROM_V%A_Y(20,59)=50388.0000000000E0_DP
 S_B_FROM_V%A_X(20,70)=75582.0000000000E0_DP
 S_B_FROM_V%B_Y(20,70)=75582.0000000000E0_DP
 S_B_FROM_V%B_X(20,82)=92378.0000000000E0_DP
 S_B_FROM_V%A_Y(20,82)=-92378.0000000000E0_DP
 S_B_FROM_V%A_X(20,95)=-92378.0000000000E0_DP
 S_B_FROM_V%B_Y(20,95)=-92378.0000000000E0_DP
 S_B_FROM_V%B_X(20,109)=-75582.0000000000E0_DP
 S_B_FROM_V%A_Y(20,109)=75582.0000000000E0_DP
 S_B_FROM_V%A_X(20,124)=50388.0000000000E0_DP
 S_B_FROM_V%B_Y(20,124)=50388.0000000000E0_DP
 S_B_FROM_V%B_X(20,140)=27132.0000000000E0_DP
 S_B_FROM_V%A_Y(20,140)=-27132.0000000000E0_DP
 S_B_FROM_V%A_X(20,157)=-11628.0000000000E0_DP
 S_B_FROM_V%B_Y(20,157)=-11628.0000000000E0_DP
 S_B_FROM_V%B_X(20,175)=-3876.00000000000E0_DP
 S_B_FROM_V%A_Y(20,175)=3876.00000000000E0_DP
 S_B_FROM_V%A_X(20,194)=969.000000000000E0_DP
 S_B_FROM_V%B_Y(20,194)=969.000000000000E0_DP
 S_B_FROM_V%B_X(20,214)=171.000000000000E0_DP
 S_B_FROM_V%A_Y(20,214)=-171.000000000000E0_DP
 S_B_FROM_V%A_X(20,235)=-19.0000000000000E0_DP
 S_B_FROM_V%B_Y(20,235)=-19.0000000000000E0_DP
 S_B_FROM_V%B_X(20,257)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(20,257)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(21,6)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(21,6)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(21,9)=20.0000000000000E0_DP
 S_B_FROM_V%A_Y(21,9)=-20.0000000000000E0_DP
 S_B_FROM_V%A_X(21,13)=-190.000000000000E0_DP
 S_B_FROM_V%B_Y(21,13)=-190.000000000000E0_DP
 S_B_FROM_V%B_X(21,18)=-1140.00000000000E0_DP
 S_B_FROM_V%A_Y(21,18)=1140.00000000000E0_DP
 S_B_FROM_V%A_X(21,24)=4845.00000000000E0_DP
 S_B_FROM_V%B_Y(21,24)=4845.00000000000E0_DP
 S_B_FROM_V%B_X(21,31)=15504.0000000000E0_DP
 S_B_FROM_V%A_Y(21,31)=-15504.0000000000E0_DP
 S_B_FROM_V%A_X(21,39)=-38760.0000000000E0_DP
 S_B_FROM_V%B_Y(21,39)=-38760.0000000000E0_DP
 S_B_FROM_V%B_X(21,48)=-77520.0000000000E0_DP
 S_B_FROM_V%A_Y(21,48)=77520.0000000000E0_DP
 S_B_FROM_V%A_X(21,58)=125970.000000000E0_DP
 S_B_FROM_V%B_Y(21,58)=125970.000000000E0_DP
 S_B_FROM_V%B_X(21,69)=167960.000000000E0_DP
 S_B_FROM_V%A_Y(21,69)=-167960.000000000E0_DP
 S_B_FROM_V%A_X(21,81)=-184756.000000000E0_DP
 S_B_FROM_V%B_Y(21,81)=-184756.000000000E0_DP
 S_B_FROM_V%B_X(21,94)=-167960.000000000E0_DP
 S_B_FROM_V%A_Y(21,94)=167960.000000000E0_DP
 S_B_FROM_V%A_X(21,108)=125970.000000000E0_DP
 S_B_FROM_V%B_Y(21,108)=125970.000000000E0_DP
 S_B_FROM_V%B_X(21,123)=77520.0000000000E0_DP
 S_B_FROM_V%A_Y(21,123)=-77520.0000000000E0_DP
 S_B_FROM_V%A_X(21,139)=-38760.0000000000E0_DP
 S_B_FROM_V%B_Y(21,139)=-38760.0000000000E0_DP
 S_B_FROM_V%B_X(21,156)=-15504.0000000000E0_DP
 S_B_FROM_V%A_Y(21,156)=15504.0000000000E0_DP
 S_B_FROM_V%A_X(21,174)=4845.00000000000E0_DP
 S_B_FROM_V%B_Y(21,174)=4845.00000000000E0_DP
 S_B_FROM_V%B_X(21,193)=1140.00000000000E0_DP
 S_B_FROM_V%A_Y(21,193)=-1140.00000000000E0_DP
 S_B_FROM_V%A_X(21,213)=-190.000000000000E0_DP
 S_B_FROM_V%B_Y(21,213)=-190.000000000000E0_DP
 S_B_FROM_V%B_X(21,234)=-20.0000000000000E0_DP
 S_B_FROM_V%A_Y(21,234)=20.0000000000000E0_DP
 S_B_FROM_V%A_X(21,256)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(21,256)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(22,3)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(22,3)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(22,5)=21.0000000000000E0_DP
 S_B_FROM_V%A_Y(22,5)=-21.0000000000000E0_DP
 S_B_FROM_V%A_X(22,8)=-210.000000000000E0_DP
 S_B_FROM_V%B_Y(22,8)=-210.000000000000E0_DP
 S_B_FROM_V%B_X(22,12)=-1330.00000000000E0_DP
 S_B_FROM_V%A_Y(22,12)=1330.00000000000E0_DP
 S_B_FROM_V%A_X(22,17)=5985.00000000000E0_DP
 S_B_FROM_V%B_Y(22,17)=5985.00000000000E0_DP
 S_B_FROM_V%B_X(22,23)=20349.0000000000E0_DP
 S_B_FROM_V%A_Y(22,23)=-20349.0000000000E0_DP
 S_B_FROM_V%A_X(22,30)=-54264.0000000000E0_DP
 S_B_FROM_V%B_Y(22,30)=-54264.0000000000E0_DP
 S_B_FROM_V%B_X(22,38)=-116280.000000000E0_DP
 S_B_FROM_V%A_Y(22,38)=116280.000000000E0_DP
 S_B_FROM_V%A_X(22,47)=203490.000000000E0_DP
 S_B_FROM_V%B_Y(22,47)=203490.000000000E0_DP
 S_B_FROM_V%B_X(22,57)=293930.000000000E0_DP
 S_B_FROM_V%A_Y(22,57)=-293930.000000000E0_DP
 S_B_FROM_V%A_X(22,68)=-352716.000000000E0_DP
 S_B_FROM_V%B_Y(22,68)=-352716.000000000E0_DP
 S_B_FROM_V%B_X(22,80)=-352716.000000000E0_DP
 S_B_FROM_V%A_Y(22,80)=352716.000000000E0_DP
 S_B_FROM_V%A_X(22,93)=293930.000000000E0_DP
 S_B_FROM_V%B_Y(22,93)=293930.000000000E0_DP
 S_B_FROM_V%B_X(22,107)=203490.000000000E0_DP
 S_B_FROM_V%A_Y(22,107)=-203490.000000000E0_DP
 S_B_FROM_V%A_X(22,122)=-116280.000000000E0_DP
 S_B_FROM_V%B_Y(22,122)=-116280.000000000E0_DP
 S_B_FROM_V%B_X(22,138)=-54264.0000000000E0_DP
 S_B_FROM_V%A_Y(22,138)=54264.0000000000E0_DP
 S_B_FROM_V%A_X(22,155)=20349.0000000000E0_DP
 S_B_FROM_V%B_Y(22,155)=20349.0000000000E0_DP
 S_B_FROM_V%B_X(22,173)=5985.00000000000E0_DP
 S_B_FROM_V%A_Y(22,173)=-5985.00000000000E0_DP
 S_B_FROM_V%A_X(22,192)=-1330.00000000000E0_DP
 S_B_FROM_V%B_Y(22,192)=-1330.00000000000E0_DP
 S_B_FROM_V%B_X(22,212)=-210.000000000000E0_DP
 S_B_FROM_V%A_Y(22,212)=210.000000000000E0_DP
 S_B_FROM_V%A_X(22,233)=21.0000000000000E0_DP
 S_B_FROM_V%B_Y(22,233)=21.0000000000000E0_DP
 S_B_FROM_V%B_X(22,255)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(22,255)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(1,103)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(1,118)=0.500000000000000E0_DP
 S_B_FROM_V%VA(1,132)=1.50000000000000E0_DP
 S_B_FROM_V%VA(1,134)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(1,149)=-1.16666666666667E0_DP
 S_B_FROM_V%VA(1,151)=0.500000000000000E0_DP
 S_B_FROM_V%VA(1,165)=-1.57500000000000E0_DP
 S_B_FROM_V%VA(1,167)=0.875000000000000E0_DP
 S_B_FROM_V%VA(1,169)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(1,184)=0.874999999999998E0_DP
 S_B_FROM_V%VA(1,186)=-0.625000000000000E0_DP
 S_B_FROM_V%VA(1,188)=0.500000000000000E0_DP
 S_B_FROM_V%VA(1,202)=0.468750000000000E0_DP
 S_B_FROM_V%VA(1,204)=-0.437500000000000E0_DP
 S_B_FROM_V%VA(1,206)=0.416666666666667E0_DP
 S_B_FROM_V%VA(1,208)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(1,223)=-0.156250000000000E0_DP
 S_B_FROM_V%VA(1,225)=0.187500000000000E0_DP
 S_B_FROM_V%VA(1,227)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(1,229)=0.500000000000000E0_DP
 S_B_FROM_V%VA(1,243)=-2.734375000000002E-002_DP
 S_B_FROM_V%VA(1,245)=3.906250000000002E-002_DP
 S_B_FROM_V%VA(1,247)=-6.249999999999997E-002_DP
 S_B_FROM_V%VA(1,249)=0.125000000000000E0_DP
 S_B_FROM_V%VA(1,251)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(1,253)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(1,266)=3.038194444444437E-003_DP
 S_B_FROM_V%VA(1,268)=-5.580357142857142E-003_DP
 S_B_FROM_V%VA(1,270)=1.250000000000001E-002_DP
 S_B_FROM_V%VA(1,272)=-4.166666666666664E-002_DP
 S_B_FROM_V%VA(1,274)=0.500000000000000E0_DP
 S_B_FROM_V%VB(1,275)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(2,103)=0.500000000000000E0_DP
 S_B_FROM_V%VB(2,117)=0.166666666666667E0_DP
 S_B_FROM_V%VA(2,118)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(2,132)=-1.50000000000000E0_DP
 S_B_FROM_V%VB(2,133)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(2,134)=0.500000000000000E0_DP
 S_B_FROM_V%VB(2,148)=-0.233333333333333E0_DP
 S_B_FROM_V%VA(2,149)=1.16666666666667E0_DP
 S_B_FROM_V%VB(2,150)=0.166666666666667E0_DP
 S_B_FROM_V%VA(2,151)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(2,165)=1.57500000000000E0_DP
 S_B_FROM_V%VB(2,166)=0.175000000000000E0_DP
 S_B_FROM_V%VA(2,167)=-0.875000000000000E0_DP
 S_B_FROM_V%VB(2,168)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(2,169)=0.500000000000000E0_DP
 S_B_FROM_V%VB(2,183)=0.125000000000000E0_DP
 S_B_FROM_V%VA(2,184)=-0.875000000000000E0_DP
 S_B_FROM_V%VB(2,185)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(2,186)=0.625000000000000E0_DP
 S_B_FROM_V%VB(2,187)=0.166666666666667E0_DP
 S_B_FROM_V%VA(2,188)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(2,202)=-0.468749999999999E0_DP
 S_B_FROM_V%VB(2,203)=-6.249999999999996E-002_DP
 S_B_FROM_V%VA(2,204)=0.437500000000000E0_DP
 S_B_FROM_V%VB(2,205)=8.333333333333333E-002_DP
 S_B_FROM_V%VA(2,206)=-0.416666666666667E0_DP
 S_B_FROM_V%VB(2,207)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(2,208)=0.500000000000000E0_DP
 S_B_FROM_V%VB(2,222)=-1.736111111111113E-002_DP
 S_B_FROM_V%VA(2,223)=0.156250000000000E0_DP
 S_B_FROM_V%VB(2,224)=2.678571428571427E-002_DP
 S_B_FROM_V%VA(2,225)=-0.187500000000000E0_DP
 S_B_FROM_V%VB(2,226)=-4.999999999999998E-002_DP
 S_B_FROM_V%VA(2,227)=0.250000000000000E0_DP
 S_B_FROM_V%VB(2,228)=0.166666666666667E0_DP
 S_B_FROM_V%VA(2,229)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(2,231)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(2,243)=2.734374999999996E-002_DP
 S_B_FROM_V%VB(2,244)=4.340277777777775E-003_DP
 S_B_FROM_V%VA(2,245)=-3.906250000000003E-002_DP
 S_B_FROM_V%VB(2,246)=-8.928571428571428E-003_DP
 S_B_FROM_V%VA(2,247)=6.249999999999998E-002_DP
 S_B_FROM_V%VB(2,248)=2.500000000000000E-002_DP
 S_B_FROM_V%VA(2,249)=-0.125000000000000E0_DP
 S_B_FROM_V%VB(2,250)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(2,251)=0.500000000000000E0_DP
 S_B_FROM_V%VB(2,252)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(2,265)=2.761994949494954E-004_DP
 S_B_FROM_V%VA(2,266)=-3.038194444444442E-003_DP
 S_B_FROM_V%VB(2,267)=-6.200396825396814E-004_DP
 S_B_FROM_V%VA(2,268)=5.580357142857139E-003_DP
 S_B_FROM_V%VB(2,269)=1.785714285714285E-003_DP
 S_B_FROM_V%VA(2,270)=-1.250000000000000E-002_DP
 S_B_FROM_V%VB(2,271)=-8.333333333333335E-003_DP
 S_B_FROM_V%VA(2,272)=4.166666666666666E-002_DP
 S_B_FROM_V%VB(2,273)=0.166666666666667E0_DP
 S_B_FROM_V%VA(2,274)=0.500000000000000E0_DP
 S_B_FROM_V%VA(3,103)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(3,117)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(3,118)=0.500000000000000E0_DP
 S_B_FROM_V%VA(3,132)=1.62500000000000E0_DP
 S_B_FROM_V%VB(3,133)=0.333333333333333E0_DP
 S_B_FROM_V%VA(3,134)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(3,148)=0.466666666666667E0_DP
 S_B_FROM_V%VA(3,149)=-1.29166666666667E0_DP
 S_B_FROM_V%VB(3,150)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(3,151)=0.500000000000000E0_DP
 S_B_FROM_V%VA(3,165)=-1.66250000000000E0_DP
 S_B_FROM_V%VB(3,166)=-0.350000000000000E0_DP
 S_B_FROM_V%VA(3,167)=1.00000000000000E0_DP
 S_B_FROM_V%VB(3,168)=0.333333333333333E0_DP
 S_B_FROM_V%VA(3,169)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(3,183)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(3,184)=0.937499999999999E0_DP
 S_B_FROM_V%VB(3,185)=0.250000000000000E0_DP
 S_B_FROM_V%VA(3,186)=-0.750000000000000E0_DP
 S_B_FROM_V%VB(3,187)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(3,188)=0.500000000000000E0_DP
 S_B_FROM_V%VA(3,202)=0.492187500000000E0_DP
 S_B_FROM_V%VB(3,203)=0.125000000000000E0_DP
 S_B_FROM_V%VA(3,204)=-0.479166666666667E0_DP
 S_B_FROM_V%VB(3,205)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(3,206)=0.541666666666667E0_DP
 S_B_FROM_V%VB(3,207)=0.333333333333333E0_DP
 S_B_FROM_V%VA(3,208)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(3,210)=-0.333333333333333E0_DP
 S_B_FROM_V%VB(3,222)=3.472222222222220E-002_DP
 S_B_FROM_V%VA(3,223)=-0.166294642857143E0_DP
 S_B_FROM_V%VB(3,224)=-5.357142857142856E-002_DP
 S_B_FROM_V%VA(3,225)=0.212500000000000E0_DP
 S_B_FROM_V%VB(3,226)=0.100000000000000E0_DP
 S_B_FROM_V%VA(3,227)=-0.375000000000000E0_DP
 S_B_FROM_V%VB(3,228)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(3,229)=0.500000000000000E0_DP
 S_B_FROM_V%VB(3,230)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(3,243)=-2.864583333333333E-002_DP
 S_B_FROM_V%VB(3,244)=-8.680555555555556E-003_DP
 S_B_FROM_V%VA(3,245)=4.241071428571433E-002_DP
 S_B_FROM_V%VB(3,246)=1.785714285714286E-002_DP
 S_B_FROM_V%VA(3,247)=-7.499999999999998E-002_DP
 S_B_FROM_V%VB(3,248)=-5.000000000000000E-002_DP
 S_B_FROM_V%VA(3,249)=0.250000000000000E0_DP
 S_B_FROM_V%VB(3,250)=0.333333333333333E0_DP
 S_B_FROM_V%VA(3,251)=1.00000000000000E0_DP
 S_B_FROM_V%VB(3,265)=-5.523989898989897E-004_DP
 S_B_FROM_V%VA(3,266)=3.224206349206343E-003_DP
 S_B_FROM_V%VB(3,267)=1.240079365079365E-003_DP
 S_B_FROM_V%VA(3,268)=-6.250000000000000E-003_DP
 S_B_FROM_V%VB(3,269)=-3.571428571428571E-003_DP
 S_B_FROM_V%VA(3,270)=1.666666666666667E-002_DP
 S_B_FROM_V%VB(3,271)=1.666666666666667E-002_DP
 S_B_FROM_V%VA(3,272)=-0.166666666666667E0_DP
 S_B_FROM_V%VB(3,273)=0.333333333333333E0_DP
 S_B_FROM_V%VA(4,103)=0.500000000000000E0_DP
 S_B_FROM_V%VB(4,117)=0.500000000000000E0_DP
 S_B_FROM_V%VA(4,118)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,132)=-1.87500000000000E0_DP
 S_B_FROM_V%VB(4,133)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,134)=0.500000000000000E0_DP
 S_B_FROM_V%VB(4,148)=-0.775000000000000E0_DP
 S_B_FROM_V%VA(4,149)=1.54166666666667E0_DP
 S_B_FROM_V%VB(4,150)=0.500000000000000E0_DP
 S_B_FROM_V%VA(4,151)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,165)=1.83750000000000E0_DP
 S_B_FROM_V%VB(4,166)=0.600000000000000E0_DP
 S_B_FROM_V%VA(4,167)=-1.25000000000000E0_DP
 S_B_FROM_V%VB(4,168)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,169)=0.500000000000000E0_DP
 S_B_FROM_V%VB(4,183)=0.401785714285715E0_DP
 S_B_FROM_V%VA(4,184)=-1.06250000000000E0_DP
 S_B_FROM_V%VB(4,185)=-0.450000000000000E0_DP
 S_B_FROM_V%VA(4,186)=1.00000000000000E0_DP
 S_B_FROM_V%VB(4,187)=0.500000000000000E0_DP
 S_B_FROM_V%VA(4,188)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,190)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(4,202)=-0.539062500000000E0_DP
 S_B_FROM_V%VB(4,203)=-0.205357142857143E0_DP
 S_B_FROM_V%VA(4,204)=0.562500000000000E0_DP
 S_B_FROM_V%VB(4,205)=0.325000000000000E0_DP
 S_B_FROM_V%VA(4,206)=-0.791666666666667E0_DP
 S_B_FROM_V%VB(4,207)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,208)=0.500000000000000E0_DP
 S_B_FROM_V%VB(4,209)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(4,222)=-5.543154761904764E-002_DP
 S_B_FROM_V%VA(4,223)=0.186383928571429E0_DP
 S_B_FROM_V%VB(4,224)=9.107142857142854E-002_DP
 S_B_FROM_V%VA(4,225)=-0.262500000000000E0_DP
 S_B_FROM_V%VB(4,226)=-0.225000000000000E0_DP
 S_B_FROM_V%VA(4,227)=0.625000000000000E0_DP
 S_B_FROM_V%VB(4,228)=0.500000000000000E0_DP
 S_B_FROM_V%VA(4,229)=1.50000000000000E0_DP
 S_B_FROM_V%VA(4,243)=3.124999999999997E-002_DP
 S_B_FROM_V%VB(4,244)=1.413690476190477E-002_DP
 S_B_FROM_V%VA(4,245)=-4.910714285714287E-002_DP
 S_B_FROM_V%VB(4,246)=-3.214285714285715E-002_DP
 S_B_FROM_V%VA(4,247)=0.100000000000000E0_DP
 S_B_FROM_V%VB(4,248)=0.150000000000000E0_DP
 S_B_FROM_V%VA(4,249)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(4,250)=1.00000000000000E0_DP
 S_B_FROM_V%VB(4,265)=8.793290043290049E-004_DP
 S_B_FROM_V%VA(4,266)=-3.596230158730165E-003_DP
 S_B_FROM_V%VB(4,267)=-2.083333333333332E-003_DP
 S_B_FROM_V%VA(4,268)=7.589285714285718E-003_DP
 S_B_FROM_V%VB(4,269)=7.142857142857141E-003_DP
 S_B_FROM_V%VA(4,270)=-2.500000000000000E-002_DP
 S_B_FROM_V%VB(4,271)=-0.100000000000000E0_DP
 S_B_FROM_V%VA(4,272)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(5,103)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(5,117)=-0.666666666666667E0_DP
 S_B_FROM_V%VA(5,118)=0.500000000000000E0_DP
 S_B_FROM_V%VA(5,132)=2.25000000000000E0_DP
 S_B_FROM_V%VB(5,133)=0.666666666666667E0_DP
 S_B_FROM_V%VA(5,134)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(5,148)=1.23333333333333E0_DP
 S_B_FROM_V%VA(5,149)=-1.91666666666667E0_DP
 S_B_FROM_V%VB(5,150)=-0.666666666666667E0_DP
 S_B_FROM_V%VA(5,151)=0.500000000000000E0_DP
 S_B_FROM_V%VA(5,165)=-2.16250000000000E0_DP
 S_B_FROM_V%VB(5,166)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(5,167)=1.62500000000000E0_DP
 S_B_FROM_V%VB(5,168)=0.666666666666667E0_DP
 S_B_FROM_V%VA(5,169)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(5,171)=-0.200000000000000E0_DP
 S_B_FROM_V%VB(5,183)=-0.607142857142857E0_DP
 S_B_FROM_V%VA(5,184)=1.31250000000000E0_DP
 S_B_FROM_V%VB(5,185)=0.800000000000000E0_DP
 S_B_FROM_V%VA(5,186)=-1.37500000000000E0_DP
 S_B_FROM_V%VB(5,187)=-0.666666666666667E0_DP
 S_B_FROM_V%VA(5,188)=0.500000000000000E0_DP
 S_B_FROM_V%VB(5,189)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(5,202)=0.620535714285715E0_DP
 S_B_FROM_V%VB(5,203)=0.321428571428572E0_DP
 S_B_FROM_V%VA(5,204)=-0.750000000000000E0_DP
 S_B_FROM_V%VB(5,205)=-0.633333333333333E0_DP
 S_B_FROM_V%VA(5,206)=1.16666666666667E0_DP
 S_B_FROM_V%VB(5,207)=0.666666666666667E0_DP
 S_B_FROM_V%VA(5,208)=2.00000000000000E0_DP
 S_B_FROM_V%VB(5,222)=8.283730158730163E-002_DP
 S_B_FROM_V%VA(5,223)=-0.223214285714286E0_DP
 S_B_FROM_V%VB(5,224)=-0.150000000000000E0_DP
 S_B_FROM_V%VA(5,225)=0.400000000000000E0_DP
 S_B_FROM_V%VB(5,226)=0.500000000000000E0_DP
 S_B_FROM_V%VA(5,227)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(5,228)=2.00000000000000E0_DP
 S_B_FROM_V%VA(5,243)=-3.571428571428573E-002_DP
 S_B_FROM_V%VB(5,244)=-2.182539682539683E-002_DP
 S_B_FROM_V%VA(5,245)=6.249999999999996E-002_DP
 S_B_FROM_V%VB(5,246)=5.714285714285713E-002_DP
 S_B_FROM_V%VA(5,247)=-0.200000000000000E0_DP
 S_B_FROM_V%VB(5,248)=-0.400000000000000E0_DP
 S_B_FROM_V%VA(5,249)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(5,265)=-1.307720057720058E-003_DP
 S_B_FROM_V%VA(5,266)=4.265873015873016E-003_DP
 S_B_FROM_V%VB(5,267)=3.373015873015870E-003_DP
 S_B_FROM_V%VA(5,268)=-1.071428571428571E-002_DP
 S_B_FROM_V%VB(5,269)=-1.428571428571428E-002_DP
 S_B_FROM_V%VA(5,270)=0.100000000000000E0_DP
 S_B_FROM_V%VB(5,271)=-0.200000000000000E0_DP
 S_B_FROM_V%VA(6,103)=0.500000000000000E0_DP
 S_B_FROM_V%VB(6,117)=0.833333333333333E0_DP
 S_B_FROM_V%VA(6,118)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(6,132)=-2.75000000000000E0_DP
 S_B_FROM_V%VB(6,133)=-0.833333333333333E0_DP
 S_B_FROM_V%VA(6,134)=0.500000000000000E0_DP
 S_B_FROM_V%VB(6,148)=-1.91666666666667E0_DP
 S_B_FROM_V%VA(6,149)=2.41666666666667E0_DP
 S_B_FROM_V%VB(6,150)=0.833333333333333E0_DP
 S_B_FROM_V%VA(6,151)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(6,153)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(6,165)=2.76250000000000E0_DP
 S_B_FROM_V%VB(6,166)=1.62500000000000E0_DP
 S_B_FROM_V%VA(6,167)=-2.12500000000000E0_DP
 S_B_FROM_V%VB(6,168)=-0.833333333333333E0_DP
 S_B_FROM_V%VA(6,169)=0.500000000000000E0_DP
 S_B_FROM_V%VB(6,170)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(6,183)=0.937500000000000E0_DP
 S_B_FROM_V%VA(6,184)=-1.81250000000000E0_DP
 S_B_FROM_V%VB(6,185)=-1.37500000000000E0_DP
 S_B_FROM_V%VA(6,186)=1.87500000000000E0_DP
 S_B_FROM_V%VB(6,187)=0.833333333333333E0_DP
 S_B_FROM_V%VA(6,188)=2.50000000000000E0_DP
 S_B_FROM_V%VA(6,202)=-0.758928571428572E0_DP
 S_B_FROM_V%VB(6,203)=-0.535714285714286E0_DP
 S_B_FROM_V%VA(6,204)=1.16666666666667E0_DP
 S_B_FROM_V%VB(6,205)=1.16666666666667E0_DP
 S_B_FROM_V%VA(6,206)=-1.66666666666667E0_DP
 S_B_FROM_V%VB(6,207)=3.33333333333333E0_DP
 S_B_FROM_V%VB(6,222)=-0.124007936507936E0_DP
 S_B_FROM_V%VA(6,223)=0.290178571428572E0_DP
 S_B_FROM_V%VB(6,224)=0.285714285714286E0_DP
 S_B_FROM_V%VA(6,225)=-0.750000000000000E0_DP
 S_B_FROM_V%VB(6,226)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(6,227)=-2.50000000000000E0_DP
 S_B_FROM_V%VA(6,243)=4.315476190476190E-002_DP
 S_B_FROM_V%VB(6,244)=3.472222222222222E-002_DP
 S_B_FROM_V%VA(6,245)=-8.928571428571426E-002_DP
 S_B_FROM_V%VB(6,246)=-0.142857142857143E0_DP
 S_B_FROM_V%VA(6,247)=0.500000000000000E0_DP
 S_B_FROM_V%VB(6,248)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(6,265)=1.939033189033189E-003_DP
 S_B_FROM_V%VA(6,266)=-5.456349206349210E-003_DP
 S_B_FROM_V%VB(6,267)=-5.952380952380953E-003_DP
 S_B_FROM_V%VA(6,268)=1.785714285714286E-002_DP
 S_B_FROM_V%VB(6,269)=7.142857142857142E-002_DP
 S_B_FROM_V%VA(6,270)=0.166666666666667E0_DP
 S_B_FROM_V%VA(7,103)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(7,117)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(7,118)=0.500000000000000E0_DP
 S_B_FROM_V%VA(7,132)=3.37500000000000E0_DP
 S_B_FROM_V%VB(7,133)=1.00000000000000E0_DP
 S_B_FROM_V%VA(7,134)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(7,136)=-0.142857142857143E0_DP
 S_B_FROM_V%VB(7,148)=2.90000000000000E0_DP
 S_B_FROM_V%VA(7,149)=-3.04166666666667E0_DP
 S_B_FROM_V%VB(7,150)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(7,151)=0.500000000000000E0_DP
 S_B_FROM_V%VB(7,152)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(7,165)=-3.82500000000000E0_DP
 S_B_FROM_V%VB(7,166)=-2.55000000000000E0_DP
 S_B_FROM_V%VA(7,167)=2.75000000000000E0_DP
 S_B_FROM_V%VB(7,168)=1.00000000000000E0_DP
 S_B_FROM_V%VA(7,169)=3.00000000000000E0_DP
 S_B_FROM_V%VB(7,183)=-1.55357142857143E0_DP
 S_B_FROM_V%VA(7,184)=2.75000000000000E0_DP
 S_B_FROM_V%VB(7,185)=2.25000000000000E0_DP
 S_B_FROM_V%VA(7,186)=-2.50000000000000E0_DP
 S_B_FROM_V%VB(7,187)=5.00000000000000E0_DP
 S_B_FROM_V%VA(7,202)=1.02678571428571E0_DP
 S_B_FROM_V%VB(7,203)=1.00000000000000E0_DP
 S_B_FROM_V%VA(7,204)=-2.00000000000000E0_DP
 S_B_FROM_V%VB(7,205)=-2.00000000000000E0_DP
 S_B_FROM_V%VA(7,206)=-5.00000000000000E0_DP
 S_B_FROM_V%VB(7,222)=0.193452380952381E0_DP
 S_B_FROM_V%VA(7,223)=-0.446428571428571E0_DP
 S_B_FROM_V%VB(7,224)=-0.642857142857143E0_DP
 S_B_FROM_V%VA(7,225)=1.50000000000000E0_DP
 S_B_FROM_V%VB(7,226)=-3.00000000000000E0_DP
 S_B_FROM_V%VA(7,243)=-5.654761904761905E-002_DP
 S_B_FROM_V%VB(7,244)=-5.952380952380951E-002_DP
 S_B_FROM_V%VA(7,245)=0.178571428571429E0_DP
 S_B_FROM_V%VB(7,246)=0.428571428571429E0_DP
 S_B_FROM_V%VA(7,247)=1.00000000000000E0_DP
 S_B_FROM_V%VB(7,265)=-2.976190476190478E-003_DP
 S_B_FROM_V%VA(7,266)=7.936507936507933E-003_DP
 S_B_FROM_V%VB(7,267)=1.190476190476191E-002_DP
 S_B_FROM_V%VA(7,268)=-7.142857142857142E-002_DP
 S_B_FROM_V%VB(7,269)=0.142857142857143E0_DP
 S_B_FROM_V%VA(8,103)=0.500000000000000E0_DP
 S_B_FROM_V%VB(8,117)=1.16666666666667E0_DP
 S_B_FROM_V%VA(8,118)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(8,120)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(8,132)=-4.12500000000000E0_DP
 S_B_FROM_V%VB(8,133)=-1.16666666666667E0_DP
 S_B_FROM_V%VA(8,134)=0.500000000000000E0_DP
 S_B_FROM_V%VB(8,135)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(8,148)=-4.25833333333333E0_DP
 S_B_FROM_V%VA(8,149)=3.79166666666667E0_DP
 S_B_FROM_V%VB(8,150)=1.16666666666667E0_DP
 S_B_FROM_V%VA(8,151)=3.50000000000000E0_DP
 S_B_FROM_V%VA(8,165)=5.60000000000000E0_DP
 S_B_FROM_V%VB(8,166)=3.85000000000000E0_DP
 S_B_FROM_V%VA(8,167)=-3.50000000000000E0_DP
 S_B_FROM_V%VB(8,168)=7.00000000000000E0_DP
 S_B_FROM_V%VB(8,183)=2.75000000000000E0_DP
 S_B_FROM_V%VA(8,184)=-4.37500000000000E0_DP
 S_B_FROM_V%VB(8,185)=-3.50000000000000E0_DP
 S_B_FROM_V%VA(8,186)=-8.75000000000000E0_DP
 S_B_FROM_V%VA(8,202)=-1.62500000000000E0_DP
 S_B_FROM_V%VB(8,203)=-2.00000000000000E0_DP
 S_B_FROM_V%VA(8,204)=3.50000000000000E0_DP
 S_B_FROM_V%VB(8,205)=-7.00000000000000E0_DP
 S_B_FROM_V%VB(8,222)=-0.347222222222222E0_DP
 S_B_FROM_V%VA(8,223)=0.875000000000000E0_DP
 S_B_FROM_V%VB(8,224)=1.50000000000000E0_DP
 S_B_FROM_V%VA(8,225)=3.50000000000000E0_DP
 S_B_FROM_V%VA(8,243)=8.333333333333330E-002_DP
 S_B_FROM_V%VB(8,244)=0.138888888888889E0_DP
 S_B_FROM_V%VA(8,245)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(8,246)=1.00000000000000E0_DP
 S_B_FROM_V%VB(8,265)=5.050505050505052E-003_DP
 S_B_FROM_V%VA(8,266)=-1.388888888888889E-002_DP
 S_B_FROM_V%VB(8,267)=-5.555555555555555E-002_DP
 S_B_FROM_V%VA(8,268)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(9,103)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(9,105)=-0.111111111111111E0_DP
 S_B_FROM_V%VB(9,117)=-1.33333333333333E0_DP
 S_B_FROM_V%VA(9,118)=0.500000000000000E0_DP
 S_B_FROM_V%VB(9,119)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(9,132)=5.00000000000000E0_DP
 S_B_FROM_V%VB(9,133)=1.33333333333333E0_DP
 S_B_FROM_V%VA(9,134)=4.00000000000000E0_DP
 S_B_FROM_V%VB(9,148)=6.06666666666667E0_DP
 S_B_FROM_V%VA(9,149)=-4.66666666666667E0_DP
 S_B_FROM_V%VB(9,150)=9.33333333333333E0_DP
 S_B_FROM_V%VA(9,165)=-8.40000000000000E0_DP
 S_B_FROM_V%VB(9,166)=-5.60000000000000E0_DP
 S_B_FROM_V%VA(9,167)=-14.0000000000000E0_DP
 S_B_FROM_V%VB(9,183)=-5.00000000000000E0_DP
 S_B_FROM_V%VA(9,184)=7.00000000000000E0_DP
 S_B_FROM_V%VB(9,185)=-14.0000000000000E0_DP
 S_B_FROM_V%VA(9,202)=3.00000000000000E0_DP
 S_B_FROM_V%VB(9,203)=4.00000000000000E0_DP
 S_B_FROM_V%VA(9,204)=9.33333333333333E0_DP
 S_B_FROM_V%VB(9,222)=0.777777777777778E0_DP
 S_B_FROM_V%VA(9,223)=-2.00000000000000E0_DP
 S_B_FROM_V%VB(9,224)=4.00000000000000E0_DP
 S_B_FROM_V%VA(9,243)=-0.166666666666667E0_DP
 S_B_FROM_V%VB(9,244)=-0.444444444444444E0_DP
 S_B_FROM_V%VA(9,245)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(9,265)=-1.010101010101010E-002_DP
 S_B_FROM_V%VA(9,266)=5.555555555555555E-002_DP
 S_B_FROM_V%VB(9,267)=-0.111111111111111E0_DP
 S_B_FROM_V%VA(10,91)=-0.100000000000000E0_DP
 S_B_FROM_V%VA(10,103)=0.500000000000000E0_DP
 S_B_FROM_V%VB(10,104)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(10,117)=1.50000000000000E0_DP
 S_B_FROM_V%VA(10,118)=4.50000000000000E0_DP
 S_B_FROM_V%VA(10,132)=-6.00000000000000E0_DP
 S_B_FROM_V%VB(10,133)=12.0000000000000E0_DP
 S_B_FROM_V%VB(10,148)=-8.40000000000000E0_DP
 S_B_FROM_V%VA(10,149)=-21.0000000000000E0_DP
 S_B_FROM_V%VA(10,165)=12.6000000000000E0_DP
 S_B_FROM_V%VB(10,166)=-25.2000000000000E0_DP
 S_B_FROM_V%VB(10,183)=9.00000000000000E0_DP
 S_B_FROM_V%VA(10,184)=21.0000000000000E0_DP
 S_B_FROM_V%VA(10,202)=-6.00000000000000E0_DP
 S_B_FROM_V%VB(10,203)=12.0000000000000E0_DP
 S_B_FROM_V%VB(10,222)=-2.00000000000000E0_DP
 S_B_FROM_V%VA(10,223)=-4.50000000000000E0_DP
 S_B_FROM_V%VA(10,243)=0.500000000000000E0_DP
 S_B_FROM_V%VB(10,244)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(10,265)=4.545454545454546E-002_DP
 S_B_FROM_V%VA(10,266)=0.100000000000000E0_DP
 S_B_FROM_V%VA(11,78)=-9.090909090909091E-002_DP
 S_B_FROM_V%VB(11,90)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(11,103)=5.00000000000000E0_DP
 S_B_FROM_V%VB(11,117)=15.0000000000000E0_DP
 S_B_FROM_V%VA(11,132)=-30.0000000000000E0_DP
 S_B_FROM_V%VB(11,148)=-42.0000000000000E0_DP
 S_B_FROM_V%VA(11,165)=42.0000000000000E0_DP
 S_B_FROM_V%VB(11,183)=30.0000000000000E0_DP
 S_B_FROM_V%VA(11,202)=-15.0000000000000E0_DP
 S_B_FROM_V%VB(11,222)=-5.00000000000000E0_DP
 S_B_FROM_V%VA(11,243)=1.00000000000000E0_DP
 S_B_FROM_V%VB(11,265)=9.090909090909091E-002_DP
 S_B_FROM_V%VA(12,66)=-8.333333333333333E-002_DP
 S_B_FROM_V%VB(12,77)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(12,89)=5.50000000000000E0_DP
 S_B_FROM_V%VB(12,102)=18.3333333333333E0_DP
 S_B_FROM_V%VA(12,116)=-41.2500000000000E0_DP
 S_B_FROM_V%VB(12,131)=-66.0000000000000E0_DP
 S_B_FROM_V%VA(12,147)=77.0000000000000E0_DP
 S_B_FROM_V%VB(12,164)=66.0000000000000E0_DP
 S_B_FROM_V%VA(12,182)=-41.2500000000000E0_DP
 S_B_FROM_V%VB(12,201)=-18.3333333333333E0_DP
 S_B_FROM_V%VA(12,221)=5.50000000000000E0_DP
 S_B_FROM_V%VB(12,242)=1.00000000000000E0_DP
 S_B_FROM_V%VA(12,264)=-8.333333333333333E-002_DP
 S_B_FROM_V%VA(13,55)=-7.692307692307693E-002_DP
 S_B_FROM_V%VB(13,65)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(13,76)=6.00000000000000E0_DP
 S_B_FROM_V%VB(13,88)=22.0000000000000E0_DP
 S_B_FROM_V%VA(13,101)=-55.0000000000000E0_DP
 S_B_FROM_V%VB(13,115)=-99.0000000000000E0_DP
 S_B_FROM_V%VA(13,130)=132.000000000000E0_DP
 S_B_FROM_V%VB(13,146)=132.000000000000E0_DP
 S_B_FROM_V%VA(13,163)=-99.0000000000000E0_DP
 S_B_FROM_V%VB(13,181)=-55.0000000000000E0_DP
 S_B_FROM_V%VA(13,200)=22.0000000000000E0_DP
 S_B_FROM_V%VB(13,220)=6.00000000000000E0_DP
 S_B_FROM_V%VA(13,241)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(13,263)=-7.692307692307693E-002_DP
 S_B_FROM_V%VA(14,45)=-7.142857142857142E-002_DP
 S_B_FROM_V%VB(14,54)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(14,64)=6.50000000000000E0_DP
 S_B_FROM_V%VB(14,75)=26.0000000000000E0_DP
 S_B_FROM_V%VA(14,87)=-71.5000000000000E0_DP
 S_B_FROM_V%VB(14,100)=-143.000000000000E0_DP
 S_B_FROM_V%VA(14,114)=214.500000000000E0_DP
 S_B_FROM_V%VB(14,129)=245.142857142857E0_DP
 S_B_FROM_V%VA(14,145)=-214.500000000000E0_DP
 S_B_FROM_V%VB(14,162)=-143.000000000000E0_DP
 S_B_FROM_V%VA(14,180)=71.5000000000000E0_DP
 S_B_FROM_V%VB(14,199)=26.0000000000000E0_DP
 S_B_FROM_V%VA(14,219)=-6.50000000000000E0_DP
 S_B_FROM_V%VB(14,240)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(14,262)=7.142857142857142E-002_DP
 S_B_FROM_V%VA(15,36)=-6.666666666666667E-002_DP
 S_B_FROM_V%VB(15,44)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(15,53)=7.00000000000000E0_DP
 S_B_FROM_V%VB(15,63)=30.3333333333333E0_DP
 S_B_FROM_V%VA(15,74)=-91.0000000000000E0_DP
 S_B_FROM_V%VB(15,86)=-200.200000000000E0_DP
 S_B_FROM_V%VA(15,99)=333.666666666667E0_DP
 S_B_FROM_V%VB(15,113)=429.000000000000E0_DP
 S_B_FROM_V%VA(15,128)=-429.000000000000E0_DP
 S_B_FROM_V%VB(15,144)=-333.666666666667E0_DP
 S_B_FROM_V%VA(15,161)=200.200000000000E0_DP
 S_B_FROM_V%VB(15,179)=91.0000000000000E0_DP
 S_B_FROM_V%VA(15,198)=-30.3333333333333E0_DP
 S_B_FROM_V%VB(15,218)=-7.00000000000000E0_DP
 S_B_FROM_V%VA(15,239)=1.00000000000000E0_DP
 S_B_FROM_V%VB(15,261)=6.666666666666667E-002_DP
 S_B_FROM_V%VA(16,28)=-6.250000000000000E-002_DP
 S_B_FROM_V%VB(16,35)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(16,43)=7.50000000000000E0_DP
 S_B_FROM_V%VB(16,52)=35.0000000000000E0_DP
 S_B_FROM_V%VA(16,62)=-113.750000000000E0_DP
 S_B_FROM_V%VB(16,73)=-273.000000000000E0_DP
 S_B_FROM_V%VA(16,85)=500.500000000000E0_DP
 S_B_FROM_V%VB(16,98)=715.000000000000E0_DP
 S_B_FROM_V%VA(16,112)=-804.375000000000E0_DP
 S_B_FROM_V%VB(16,127)=-715.000000000000E0_DP
 S_B_FROM_V%VA(16,143)=500.500000000000E0_DP
 S_B_FROM_V%VB(16,160)=273.000000000000E0_DP
 S_B_FROM_V%VA(16,178)=-113.750000000000E0_DP
 S_B_FROM_V%VB(16,197)=-35.0000000000000E0_DP
 S_B_FROM_V%VA(16,217)=7.50000000000000E0_DP
 S_B_FROM_V%VB(16,238)=1.00000000000000E0_DP
 S_B_FROM_V%VA(16,260)=-6.250000000000000E-002_DP
 S_B_FROM_V%VA(17,21)=-5.882352941176471E-002_DP
 S_B_FROM_V%VB(17,27)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(17,34)=8.00000000000000E0_DP
 S_B_FROM_V%VB(17,42)=40.0000000000000E0_DP
 S_B_FROM_V%VA(17,51)=-140.000000000000E0_DP
 S_B_FROM_V%VB(17,61)=-364.000000000000E0_DP
 S_B_FROM_V%VA(17,72)=728.000000000000E0_DP
 S_B_FROM_V%VB(17,84)=1144.00000000000E0_DP
 S_B_FROM_V%VA(17,97)=-1430.00000000000E0_DP
 S_B_FROM_V%VB(17,111)=-1430.00000000000E0_DP
 S_B_FROM_V%VA(17,126)=1144.00000000000E0_DP
 S_B_FROM_V%VB(17,142)=728.000000000000E0_DP
 S_B_FROM_V%VA(17,159)=-364.000000000000E0_DP
 S_B_FROM_V%VB(17,177)=-140.000000000000E0_DP
 S_B_FROM_V%VA(17,196)=40.0000000000000E0_DP
 S_B_FROM_V%VB(17,216)=8.00000000000000E0_DP
 S_B_FROM_V%VA(17,237)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(17,259)=-5.882352941176471E-002_DP
 S_B_FROM_V%VA(18,15)=-5.555555555555555E-002_DP
 S_B_FROM_V%VB(18,20)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(18,26)=8.50000000000000E0_DP
 S_B_FROM_V%VB(18,33)=45.3333333333333E0_DP
 S_B_FROM_V%VA(18,41)=-170.000000000000E0_DP
 S_B_FROM_V%VB(18,50)=-476.000000000000E0_DP
 S_B_FROM_V%VA(18,60)=1031.33333333333E0_DP
 S_B_FROM_V%VB(18,71)=1768.00000000000E0_DP
 S_B_FROM_V%VA(18,83)=-2431.00000000000E0_DP
 S_B_FROM_V%VB(18,96)=-2701.11111111111E0_DP
 S_B_FROM_V%VA(18,110)=2431.00000000000E0_DP
 S_B_FROM_V%VB(18,125)=1768.00000000000E0_DP
 S_B_FROM_V%VA(18,141)=-1031.33333333333E0_DP
 S_B_FROM_V%VB(18,158)=-476.000000000000E0_DP
 S_B_FROM_V%VA(18,176)=170.000000000000E0_DP
 S_B_FROM_V%VB(18,195)=45.3333333333333E0_DP
 S_B_FROM_V%VA(18,215)=-8.50000000000000E0_DP
 S_B_FROM_V%VB(18,236)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(18,258)=5.555555555555555E-002_DP
 S_B_FROM_V%VA(19,10)=-5.263157894736842E-002_DP
 S_B_FROM_V%VB(19,14)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(19,19)=9.00000000000000E0_DP
 S_B_FROM_V%VB(19,25)=51.0000000000000E0_DP
 S_B_FROM_V%VA(19,32)=-204.000000000000E0_DP
 S_B_FROM_V%VB(19,40)=-612.000000000000E0_DP
 S_B_FROM_V%VA(19,49)=1428.00000000000E0_DP
 S_B_FROM_V%VB(19,59)=2652.00000000000E0_DP
 S_B_FROM_V%VA(19,70)=-3978.00000000000E0_DP
 S_B_FROM_V%VB(19,82)=-4862.00000000000E0_DP
 S_B_FROM_V%VA(19,95)=4862.00000000000E0_DP
 S_B_FROM_V%VB(19,109)=3978.00000000000E0_DP
 S_B_FROM_V%VA(19,124)=-2652.00000000000E0_DP
 S_B_FROM_V%VB(19,140)=-1428.00000000000E0_DP
 S_B_FROM_V%VA(19,157)=612.000000000000E0_DP
 S_B_FROM_V%VB(19,175)=204.000000000000E0_DP
 S_B_FROM_V%VA(19,194)=-51.0000000000000E0_DP
 S_B_FROM_V%VB(19,214)=-9.00000000000000E0_DP
 S_B_FROM_V%VA(19,235)=1.00000000000000E0_DP
 S_B_FROM_V%VB(19,257)=5.263157894736842E-002_DP
 S_B_FROM_V%VA(20,6)=-5.000000000000000E-002_DP
 S_B_FROM_V%VB(20,9)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(20,13)=9.50000000000000E0_DP
 S_B_FROM_V%VB(20,18)=57.0000000000000E0_DP
 S_B_FROM_V%VA(20,24)=-242.250000000000E0_DP
 S_B_FROM_V%VB(20,31)=-775.200000000000E0_DP
 S_B_FROM_V%VA(20,39)=1938.00000000000E0_DP
 S_B_FROM_V%VB(20,48)=3876.00000000000E0_DP
 S_B_FROM_V%VA(20,58)=-6298.50000000000E0_DP
 S_B_FROM_V%VB(20,69)=-8398.00000000000E0_DP
 S_B_FROM_V%VA(20,81)=9237.80000000000E0_DP
 S_B_FROM_V%VB(20,94)=8398.00000000000E0_DP
 S_B_FROM_V%VA(20,108)=-6298.50000000000E0_DP
 S_B_FROM_V%VB(20,123)=-3876.00000000000E0_DP
 S_B_FROM_V%VA(20,139)=1938.00000000000E0_DP
 S_B_FROM_V%VB(20,156)=775.200000000000E0_DP
 S_B_FROM_V%VA(20,174)=-242.250000000000E0_DP
 S_B_FROM_V%VB(20,193)=-57.0000000000000E0_DP
 S_B_FROM_V%VA(20,213)=9.50000000000000E0_DP
 S_B_FROM_V%VB(20,234)=1.00000000000000E0_DP
 S_B_FROM_V%VA(20,256)=-5.000000000000000E-002_DP
 S_B_FROM_V%VA(21,3)=-4.761904761904762E-002_DP
 S_B_FROM_V%VB(21,5)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(21,8)=10.0000000000000E0_DP
 S_B_FROM_V%VB(21,12)=63.3333333333333E0_DP
 S_B_FROM_V%VA(21,17)=-285.000000000000E0_DP
 S_B_FROM_V%VB(21,23)=-969.000000000000E0_DP
 S_B_FROM_V%VA(21,30)=2584.00000000000E0_DP
 S_B_FROM_V%VB(21,38)=5537.14285714286E0_DP
 S_B_FROM_V%VA(21,47)=-9690.00000000000E0_DP
 S_B_FROM_V%VB(21,57)=-13996.6666666667E0_DP
 S_B_FROM_V%VA(21,68)=16796.0000000000E0_DP
 S_B_FROM_V%VB(21,80)=16796.0000000000E0_DP
 S_B_FROM_V%VA(21,93)=-13996.6666666667E0_DP
 S_B_FROM_V%VB(21,107)=-9690.00000000000E0_DP
 S_B_FROM_V%VA(21,122)=5537.14285714286E0_DP
 S_B_FROM_V%VB(21,138)=2584.00000000000E0_DP
 S_B_FROM_V%VA(21,155)=-969.000000000000E0_DP
 S_B_FROM_V%VB(21,173)=-285.000000000000E0_DP
 S_B_FROM_V%VA(21,192)=63.3333333333333E0_DP
 S_B_FROM_V%VB(21,212)=10.0000000000000E0_DP
 S_B_FROM_V%VA(21,233)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(21,255)=-4.761904761904762E-002_DP
 S_B_FROM_V%VA(22,1)=-4.545454545454546E-002_DP
 S_B_FROM_V%VB(22,2)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(22,4)=10.5000000000000E0_DP
 S_B_FROM_V%VB(22,7)=70.0000000000000E0_DP
 S_B_FROM_V%VA(22,11)=-332.500000000000E0_DP
 S_B_FROM_V%VB(22,16)=-1197.00000000000E0_DP
 S_B_FROM_V%VA(22,22)=3391.50000000000E0_DP
 S_B_FROM_V%VB(22,29)=7752.00000000000E0_DP
 S_B_FROM_V%VA(22,37)=-14535.0000000000E0_DP
 S_B_FROM_V%VB(22,46)=-22610.0000000000E0_DP
 S_B_FROM_V%VA(22,56)=29393.0000000000E0_DP
 S_B_FROM_V%VB(22,67)=32065.0909090909E0_DP
 S_B_FROM_V%VA(22,79)=-29393.0000000000E0_DP
 S_B_FROM_V%VB(22,92)=-22610.0000000000E0_DP
 S_B_FROM_V%VA(22,106)=14535.0000000000E0_DP
 S_B_FROM_V%VB(22,121)=7752.00000000000E0_DP
 S_B_FROM_V%VA(22,137)=-3391.50000000000E0_DP
 S_B_FROM_V%VB(22,154)=-1197.00000000000E0_DP
 S_B_FROM_V%VA(22,172)=332.500000000000E0_DP
 S_B_FROM_V%VB(22,191)=70.0000000000000E0_DP
 S_B_FROM_V%VA(22,211)=-10.5000000000000E0_DP
 S_B_FROM_V%VB(22,232)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(22,254)=4.545454545454546E-002_DP

end subroutine set_s_b  



subroutine set_s_e
  implicit none
  integer i
 

i=SECTOR_NMUL_MAX
!  ALLOCATE(S_e(SECTOR_NMUL_MAX))
!           DO I=1,SECTOR_NMUL_MAX
!           DO I=SECTOR_NMUL_MAX,SECTOR_NMUL_MAX   ! etienne 5/10/2015
             S_e%firsttime = 1  ! Piotr 8.19.2014 
             call nul_coef(S_e)
             call make_set_coef(S_e,I,0)
             S_e%firsttime=0
!          ENDDO

 S_E%I(1)=22;S_E%J(1)=0;
 S_E%I(2)=21;S_E%J(2)=1;
 S_E%I(3)=21;S_E%J(3)=0;
 S_E%I(4)=20;S_E%J(4)=2;
 S_E%I(5)=20;S_E%J(5)=1;
 S_E%I(6)=20;S_E%J(6)=0;
 S_E%I(7)=19;S_E%J(7)=3;
 S_E%I(8)=19;S_E%J(8)=2;
 S_E%I(9)=19;S_E%J(9)=1;
 S_E%I(10)=19;S_E%J(10)=0;
 S_E%I(11)=18;S_E%J(11)=4;
 S_E%I(12)=18;S_E%J(12)=3;
 S_E%I(13)=18;S_E%J(13)=2;
 S_E%I(14)=18;S_E%J(14)=1;
 S_E%I(15)=18;S_E%J(15)=0;
 S_E%I(16)=17;S_E%J(16)=5;
 S_E%I(17)=17;S_E%J(17)=4;
 S_E%I(18)=17;S_E%J(18)=3;
 S_E%I(19)=17;S_E%J(19)=2;
 S_E%I(20)=17;S_E%J(20)=1;
 S_E%I(21)=17;S_E%J(21)=0;
 S_E%I(22)=16;S_E%J(22)=6;
 S_E%I(23)=16;S_E%J(23)=5;
 S_E%I(24)=16;S_E%J(24)=4;
 S_E%I(25)=16;S_E%J(25)=3;
 S_E%I(26)=16;S_E%J(26)=2;
 S_E%I(27)=16;S_E%J(27)=1;
 S_E%I(28)=16;S_E%J(28)=0;
 S_E%I(29)=15;S_E%J(29)=7;
 S_E%I(30)=15;S_E%J(30)=6;
 S_E%I(31)=15;S_E%J(31)=5;
 S_E%I(32)=15;S_E%J(32)=4;
 S_E%I(33)=15;S_E%J(33)=3;
 S_E%I(34)=15;S_E%J(34)=2;
 S_E%I(35)=15;S_E%J(35)=1;
 S_E%I(36)=15;S_E%J(36)=0;
 S_E%I(37)=14;S_E%J(37)=8;
 S_E%I(38)=14;S_E%J(38)=7;
 S_E%I(39)=14;S_E%J(39)=6;
 S_E%I(40)=14;S_E%J(40)=5;
 S_E%I(41)=14;S_E%J(41)=4;
 S_E%I(42)=14;S_E%J(42)=3;
 S_E%I(43)=14;S_E%J(43)=2;
 S_E%I(44)=14;S_E%J(44)=1;
 S_E%I(45)=14;S_E%J(45)=0;
 S_E%I(46)=13;S_E%J(46)=9;
 S_E%I(47)=13;S_E%J(47)=8;
 S_E%I(48)=13;S_E%J(48)=7;
 S_E%I(49)=13;S_E%J(49)=6;
 S_E%I(50)=13;S_E%J(50)=5;
 S_E%I(51)=13;S_E%J(51)=4;
 S_E%I(52)=13;S_E%J(52)=3;
 S_E%I(53)=13;S_E%J(53)=2;
 S_E%I(54)=13;S_E%J(54)=1;
 S_E%I(55)=13;S_E%J(55)=0;
 S_E%I(56)=12;S_E%J(56)=10;
 S_E%I(57)=12;S_E%J(57)=9;
 S_E%I(58)=12;S_E%J(58)=8;
 S_E%I(59)=12;S_E%J(59)=7;
 S_E%I(60)=12;S_E%J(60)=6;
 S_E%I(61)=12;S_E%J(61)=5;
 S_E%I(62)=12;S_E%J(62)=4;
 S_E%I(63)=12;S_E%J(63)=3;
 S_E%I(64)=12;S_E%J(64)=2;
 S_E%I(65)=12;S_E%J(65)=1;
 S_E%I(66)=12;S_E%J(66)=0;
 S_E%I(67)=11;S_E%J(67)=11;
 S_E%I(68)=11;S_E%J(68)=10;
 S_E%I(69)=11;S_E%J(69)=9;
 S_E%I(70)=11;S_E%J(70)=8;
 S_E%I(71)=11;S_E%J(71)=7;
 S_E%I(72)=11;S_E%J(72)=6;
 S_E%I(73)=11;S_E%J(73)=5;
 S_E%I(74)=11;S_E%J(74)=4;
 S_E%I(75)=11;S_E%J(75)=3;
 S_E%I(76)=11;S_E%J(76)=2;
 S_E%I(77)=11;S_E%J(77)=1;
 S_E%I(78)=11;S_E%J(78)=0;
 S_E%I(79)=10;S_E%J(79)=12;
 S_E%I(80)=10;S_E%J(80)=11;
 S_E%I(81)=10;S_E%J(81)=10;
 S_E%I(82)=10;S_E%J(82)=9;
 S_E%I(83)=10;S_E%J(83)=8;
 S_E%I(84)=10;S_E%J(84)=7;
 S_E%I(85)=10;S_E%J(85)=6;
 S_E%I(86)=10;S_E%J(86)=5;
 S_E%I(87)=10;S_E%J(87)=4;
 S_E%I(88)=10;S_E%J(88)=3;
 S_E%I(89)=10;S_E%J(89)=2;
 S_E%I(90)=10;S_E%J(90)=1;
 S_E%I(91)=10;S_E%J(91)=0;
 S_E%I(92)=9;S_E%J(92)=13;
 S_E%I(93)=9;S_E%J(93)=12;
 S_E%I(94)=9;S_E%J(94)=11;
 S_E%I(95)=9;S_E%J(95)=10;
 S_E%I(96)=9;S_E%J(96)=9;
 S_E%I(97)=9;S_E%J(97)=8;
 S_E%I(98)=9;S_E%J(98)=7;
 S_E%I(99)=9;S_E%J(99)=6;
 S_E%I(100)=9;S_E%J(100)=5;
 S_E%I(101)=9;S_E%J(101)=4;
 S_E%I(102)=9;S_E%J(102)=3;
 S_E%I(103)=9;S_E%J(103)=2;
 S_E%I(104)=9;S_E%J(104)=1;
 S_E%I(105)=9;S_E%J(105)=0;
 S_E%I(106)=8;S_E%J(106)=14;
 S_E%I(107)=8;S_E%J(107)=13;
 S_E%I(108)=8;S_E%J(108)=12;
 S_E%I(109)=8;S_E%J(109)=11;
 S_E%I(110)=8;S_E%J(110)=10;
 S_E%I(111)=8;S_E%J(111)=9;
 S_E%I(112)=8;S_E%J(112)=8;
 S_E%I(113)=8;S_E%J(113)=7;
 S_E%I(114)=8;S_E%J(114)=6;
 S_E%I(115)=8;S_E%J(115)=5;
 S_E%I(116)=8;S_E%J(116)=4;
 S_E%I(117)=8;S_E%J(117)=3;
 S_E%I(118)=8;S_E%J(118)=2;
 S_E%I(119)=8;S_E%J(119)=1;
 S_E%I(120)=8;S_E%J(120)=0;
 S_E%I(121)=7;S_E%J(121)=15;
 S_E%I(122)=7;S_E%J(122)=14;
 S_E%I(123)=7;S_E%J(123)=13;
 S_E%I(124)=7;S_E%J(124)=12;
 S_E%I(125)=7;S_E%J(125)=11;
 S_E%I(126)=7;S_E%J(126)=10;
 S_E%I(127)=7;S_E%J(127)=9;
 S_E%I(128)=7;S_E%J(128)=8;
 S_E%I(129)=7;S_E%J(129)=7;
 S_E%I(130)=7;S_E%J(130)=6;
 S_E%I(131)=7;S_E%J(131)=5;
 S_E%I(132)=7;S_E%J(132)=4;
 S_E%I(133)=7;S_E%J(133)=3;
 S_E%I(134)=7;S_E%J(134)=2;
 S_E%I(135)=7;S_E%J(135)=1;
 S_E%I(136)=7;S_E%J(136)=0;
 S_E%I(137)=6;S_E%J(137)=16;
 S_E%I(138)=6;S_E%J(138)=15;
 S_E%I(139)=6;S_E%J(139)=14;
 S_E%I(140)=6;S_E%J(140)=13;
 S_E%I(141)=6;S_E%J(141)=12;
 S_E%I(142)=6;S_E%J(142)=11;
 S_E%I(143)=6;S_E%J(143)=10;
 S_E%I(144)=6;S_E%J(144)=9;
 S_E%I(145)=6;S_E%J(145)=8;
 S_E%I(146)=6;S_E%J(146)=7;
 S_E%I(147)=6;S_E%J(147)=6;
 S_E%I(148)=6;S_E%J(148)=5;
 S_E%I(149)=6;S_E%J(149)=4;
 S_E%I(150)=6;S_E%J(150)=3;
 S_E%I(151)=6;S_E%J(151)=2;
 S_E%I(152)=6;S_E%J(152)=1;
 S_E%I(153)=6;S_E%J(153)=0;
 S_E%I(154)=5;S_E%J(154)=17;
 S_E%I(155)=5;S_E%J(155)=16;
 S_E%I(156)=5;S_E%J(156)=15;
 S_E%I(157)=5;S_E%J(157)=14;
 S_E%I(158)=5;S_E%J(158)=13;
 S_E%I(159)=5;S_E%J(159)=12;
 S_E%I(160)=5;S_E%J(160)=11;
 S_E%I(161)=5;S_E%J(161)=10;
 S_E%I(162)=5;S_E%J(162)=9;
 S_E%I(163)=5;S_E%J(163)=8;
 S_E%I(164)=5;S_E%J(164)=7;
 S_E%I(165)=5;S_E%J(165)=6;
 S_E%I(166)=5;S_E%J(166)=5;
 S_E%I(167)=5;S_E%J(167)=4;
 S_E%I(168)=5;S_E%J(168)=3;
 S_E%I(169)=5;S_E%J(169)=2;
 S_E%I(170)=5;S_E%J(170)=1;
 S_E%I(171)=5;S_E%J(171)=0;
 S_E%I(172)=4;S_E%J(172)=18;
 S_E%I(173)=4;S_E%J(173)=17;
 S_E%I(174)=4;S_E%J(174)=16;
 S_E%I(175)=4;S_E%J(175)=15;
 S_E%I(176)=4;S_E%J(176)=14;
 S_E%I(177)=4;S_E%J(177)=13;
 S_E%I(178)=4;S_E%J(178)=12;
 S_E%I(179)=4;S_E%J(179)=11;
 S_E%I(180)=4;S_E%J(180)=10;
 S_E%I(181)=4;S_E%J(181)=9;
 S_E%I(182)=4;S_E%J(182)=8;
 S_E%I(183)=4;S_E%J(183)=7;
 S_E%I(184)=4;S_E%J(184)=6;
 S_E%I(185)=4;S_E%J(185)=5;
 S_E%I(186)=4;S_E%J(186)=4;
 S_E%I(187)=4;S_E%J(187)=3;
 S_E%I(188)=4;S_E%J(188)=2;
 S_E%I(189)=4;S_E%J(189)=1;
 S_E%I(190)=4;S_E%J(190)=0;
 S_E%I(191)=3;S_E%J(191)=19;
 S_E%I(192)=3;S_E%J(192)=18;
 S_E%I(193)=3;S_E%J(193)=17;
 S_E%I(194)=3;S_E%J(194)=16;
 S_E%I(195)=3;S_E%J(195)=15;
 S_E%I(196)=3;S_E%J(196)=14;
 S_E%I(197)=3;S_E%J(197)=13;
 S_E%I(198)=3;S_E%J(198)=12;
 S_E%I(199)=3;S_E%J(199)=11;
 S_E%I(200)=3;S_E%J(200)=10;
 S_E%I(201)=3;S_E%J(201)=9;
 S_E%I(202)=3;S_E%J(202)=8;
 S_E%I(203)=3;S_E%J(203)=7;
 S_E%I(204)=3;S_E%J(204)=6;
 S_E%I(205)=3;S_E%J(205)=5;
 S_E%I(206)=3;S_E%J(206)=4;
 S_E%I(207)=3;S_E%J(207)=3;
 S_E%I(208)=3;S_E%J(208)=2;
 S_E%I(209)=3;S_E%J(209)=1;
 S_E%I(210)=3;S_E%J(210)=0;
 S_E%I(211)=2;S_E%J(211)=20;
 S_E%I(212)=2;S_E%J(212)=19;
 S_E%I(213)=2;S_E%J(213)=18;
 S_E%I(214)=2;S_E%J(214)=17;
 S_E%I(215)=2;S_E%J(215)=16;
 S_E%I(216)=2;S_E%J(216)=15;
 S_E%I(217)=2;S_E%J(217)=14;
 S_E%I(218)=2;S_E%J(218)=13;
 S_E%I(219)=2;S_E%J(219)=12;
 S_E%I(220)=2;S_E%J(220)=11;
 S_E%I(221)=2;S_E%J(221)=10;
 S_E%I(222)=2;S_E%J(222)=9;
 S_E%I(223)=2;S_E%J(223)=8;
 S_E%I(224)=2;S_E%J(224)=7;
 S_E%I(225)=2;S_E%J(225)=6;
 S_E%I(226)=2;S_E%J(226)=5;
 S_E%I(227)=2;S_E%J(227)=4;
 S_E%I(228)=2;S_E%J(228)=3;
 S_E%I(229)=2;S_E%J(229)=2;
 S_E%I(230)=2;S_E%J(230)=1;
 S_E%I(231)=2;S_E%J(231)=0;
 S_E%I(232)=1;S_E%J(232)=21;
 S_E%I(233)=1;S_E%J(233)=20;
 S_E%I(234)=1;S_E%J(234)=19;
 S_E%I(235)=1;S_E%J(235)=18;
 S_E%I(236)=1;S_E%J(236)=17;
 S_E%I(237)=1;S_E%J(237)=16;
 S_E%I(238)=1;S_E%J(238)=15;
 S_E%I(239)=1;S_E%J(239)=14;
 S_E%I(240)=1;S_E%J(240)=13;
 S_E%I(241)=1;S_E%J(241)=12;
 S_E%I(242)=1;S_E%J(242)=11;
 S_E%I(243)=1;S_E%J(243)=10;
 S_E%I(244)=1;S_E%J(244)=9;
 S_E%I(245)=1;S_E%J(245)=8;
 S_E%I(246)=1;S_E%J(246)=7;
 S_E%I(247)=1;S_E%J(247)=6;
 S_E%I(248)=1;S_E%J(248)=5;
 S_E%I(249)=1;S_E%J(249)=4;
 S_E%I(250)=1;S_E%J(250)=3;
 S_E%I(251)=1;S_E%J(251)=2;
 S_E%I(252)=1;S_E%J(252)=1;
 S_E%I(253)=1;S_E%J(253)=0;
 S_E%I(254)=0;S_E%J(254)=22;
 S_E%I(255)=0;S_E%J(255)=21;
 S_E%I(256)=0;S_E%J(256)=20;
 S_E%I(257)=0;S_E%J(257)=19;
 S_E%I(258)=0;S_E%J(258)=18;
 S_E%I(259)=0;S_E%J(259)=17;
 S_E%I(260)=0;S_E%J(260)=16;
 S_E%I(261)=0;S_E%J(261)=15;
 S_E%I(262)=0;S_E%J(262)=14;
 S_E%I(263)=0;S_E%J(263)=13;
 S_E%I(264)=0;S_E%J(264)=12;
 S_E%I(265)=0;S_E%J(265)=11;
 S_E%I(266)=0;S_E%J(266)=10;
 S_E%I(267)=0;S_E%J(267)=9;
 S_E%I(268)=0;S_E%J(268)=8;
 S_E%I(269)=0;S_E%J(269)=7;
 S_E%I(270)=0;S_E%J(270)=6;
 S_E%I(271)=0;S_E%J(271)=5;
 S_E%I(272)=0;S_E%J(272)=4;
 S_E%I(273)=0;S_E%J(273)=3;
 S_E%I(274)=0;S_E%J(274)=2;
 S_E%I(275)=0;S_E%J(275)=1;
 S_E%I(276)=0;S_E%J(276)=0;
 S_E%B_Y(1,104)=1.00000000000000E0_DP
 S_E%B_X(1,118)=4.50000000000000E0_DP
 S_E%B_Y(1,119)=-0.999999999999999E0_DP
 S_E%B_Y(1,133)=-6.00000000000000E0_DP
 S_E%B_X(1,134)=-4.00000000000000E0_DP
 S_E%B_Y(1,135)=1.00000000000000E0_DP
 S_E%B_X(1,149)=-10.5000000000000E0_DP
 S_E%B_Y(1,150)=4.66666666666666E0_DP
 S_E%B_X(1,151)=3.50000000000000E0_DP
 S_E%B_Y(1,152)=-1.00000000000000E0_DP
 S_E%B_Y(1,166)=9.45000000000000E0_DP
 S_E%B_X(1,167)=6.99999999999999E0_DP
 S_E%B_Y(1,168)=-3.50000000000000E0_DP
 S_E%B_X(1,169)=-3.00000000000000E0_DP
 S_E%B_Y(1,170)=1.00000000000000E0_DP
 S_E%B_X(1,184)=7.87500000000000E0_DP
 S_E%B_Y(1,185)=-5.24999999999999E0_DP
 S_E%B_X(1,186)=-4.37500000000000E0_DP
 S_E%B_Y(1,187)=2.50000000000000E0_DP
 S_E%B_X(1,188)=2.50000000000000E0_DP
 S_E%B_Y(1,189)=-1.00000000000000E0_DP
 S_E%B_Y(1,203)=-3.75000000000000E0_DP
 S_E%B_X(1,204)=-3.49999999999999E0_DP
 S_E%B_Y(1,205)=2.62500000000000E0_DP
 S_E%B_X(1,206)=2.50000000000000E0_DP
 S_E%B_Y(1,207)=-1.66666666666667E0_DP
 S_E%B_X(1,208)=-2.00000000000000E0_DP
 S_E%B_Y(1,209)=1.00000000000000E0_DP
 S_E%B_X(1,223)=-1.40625000000000E0_DP
 S_E%B_Y(1,224)=1.25000000000000E0_DP
 S_E%B_X(1,225)=1.31250000000000E0_DP
 S_E%B_Y(1,226)=-1.12500000000000E0_DP
 S_E%B_X(1,227)=-1.25000000000000E0_DP
 S_E%B_Y(1,228)=1.00000000000000E0_DP
 S_E%B_X(1,229)=1.50000000000000E0_DP
 S_E%B_Y(1,230)=-1.00000000000000E0_DP
 S_E%B_Y(1,244)=0.273437500000000E0_DP
 S_E%B_X(1,245)=0.312499999999999E0_DP
 S_E%B_Y(1,246)=-0.312500000000000E0_DP
 S_E%B_X(1,247)=-0.375000000000000E0_DP
 S_E%B_Y(1,248)=0.375000000000000E0_DP
 S_E%B_X(1,249)=0.500000000000000E0_DP
 S_E%B_Y(1,250)=-0.500000000000000E0_DP
 S_E%B_X(1,251)=-1.00000000000000E0_DP
 S_E%B_Y(1,252)=1.00000000000000E0_DP
 S_E%B_X(1,266)=2.734375000000002E-002_DP
 S_E%B_Y(1,267)=-3.038194444444437E-002_DP
 S_E%B_X(1,268)=-3.906250000000002E-002_DP
 S_E%B_Y(1,269)=4.464285714285714E-002_DP
 S_E%B_X(1,270)=6.249999999999997E-002_DP
 S_E%B_Y(1,271)=-7.500000000000008E-002_DP
 S_E%B_X(1,272)=-0.125000000000000E0_DP
 S_E%B_Y(1,273)=0.166666666666667E0_DP
 S_E%B_X(1,274)=0.500000000000000E0_DP
 S_E%B_Y(1,275)=-1.00000000000000E0_DP
 S_E%B_X(1,276)=1.00000000000000E0_DP
 S_E%A_Y(1,276)=1.00000000000000E0_DP
 S_E%B_Y(2,104)=-1.00000000000000E0_DP
 S_E%B_X(2,118)=-4.50000000000000E0_DP
 S_E%A_Y(2,118)=-0.500000000000000E0_DP
 S_E%B_Y(2,119)=1.00000000000000E0_DP
 S_E%A_X(2,133)=-1.33333333333333E0_DP
 S_E%B_Y(2,133)=6.00000000000000E0_DP
 S_E%B_X(2,134)=4.00000000000000E0_DP
 S_E%A_Y(2,134)=0.500000000000000E0_DP
 S_E%B_Y(2,135)=-1.00000000000000E0_DP
 S_E%B_X(2,149)=10.5000000000000E0_DP
 S_E%A_Y(2,149)=1.16666666666667E0_DP
 S_E%A_X(2,150)=1.16666666666667E0_DP
 S_E%B_Y(2,150)=-4.66666666666667E0_DP
 S_E%B_X(2,151)=-3.50000000000000E0_DP
 S_E%A_Y(2,151)=-0.500000000000000E0_DP
 S_E%B_Y(2,152)=1.00000000000000E0_DP
 S_E%A_X(2,166)=1.40000000000000E0_DP
 S_E%B_Y(2,166)=-9.44999999999999E0_DP
 S_E%B_X(2,167)=-7.00000000000000E0_DP
 S_E%A_Y(2,167)=-0.875000000000000E0_DP
 S_E%A_X(2,168)=-0.999999999999999E0_DP
 S_E%B_Y(2,168)=3.50000000000000E0_DP
 S_E%B_X(2,169)=3.00000000000000E0_DP
 S_E%A_Y(2,169)=0.500000000000000E0_DP
 S_E%B_Y(2,170)=-1.00000000000000E0_DP
 S_E%B_X(2,184)=-7.87499999999999E0_DP
 S_E%A_Y(2,184)=-0.875000000000000E0_DP
 S_E%A_X(2,185)=-0.875000000000000E0_DP
 S_E%B_Y(2,185)=5.25000000000000E0_DP
 S_E%B_X(2,186)=4.37500000000000E0_DP
 S_E%A_Y(2,186)=0.625000000000000E0_DP
 S_E%A_X(2,187)=0.833333333333333E0_DP
 S_E%B_Y(2,187)=-2.50000000000000E0_DP
 S_E%B_X(2,188)=-2.50000000000000E0_DP
 S_E%A_Y(2,188)=-0.500000000000000E0_DP
 S_E%B_Y(2,189)=1.00000000000000E0_DP
 S_E%A_X(2,203)=-0.500000000000000E0_DP
 S_E%B_Y(2,203)=3.75000000000000E0_DP
 S_E%B_X(2,204)=3.50000000000000E0_DP
 S_E%A_Y(2,204)=0.437500000000000E0_DP
 S_E%A_X(2,205)=0.500000000000000E0_DP
 S_E%B_Y(2,205)=-2.62500000000000E0_DP
 S_E%B_X(2,206)=-2.50000000000000E0_DP
 S_E%A_Y(2,206)=-0.416666666666667E0_DP
 S_E%A_X(2,207)=-0.666666666666667E0_DP
 S_E%B_Y(2,207)=1.66666666666667E0_DP
 S_E%B_X(2,208)=2.00000000000000E0_DP
 S_E%A_Y(2,208)=0.500000000000000E0_DP
 S_E%B_Y(2,209)=-1.00000000000000E0_DP
 S_E%B_X(2,223)=1.40625000000000E0_DP
 S_E%A_Y(2,223)=0.156250000000000E0_DP
 S_E%A_X(2,224)=0.187500000000000E0_DP
 S_E%B_Y(2,224)=-1.25000000000000E0_DP
 S_E%B_X(2,225)=-1.31250000000000E0_DP
 S_E%A_Y(2,225)=-0.187500000000000E0_DP
 S_E%A_X(2,226)=-0.250000000000000E0_DP
 S_E%B_Y(2,226)=1.12500000000000E0_DP
 S_E%B_X(2,227)=1.25000000000000E0_DP
 S_E%A_Y(2,227)=0.250000000000000E0_DP
 S_E%A_X(2,228)=0.500000000000000E0_DP
 S_E%B_Y(2,228)=-1.00000000000000E0_DP
 S_E%B_X(2,229)=-1.50000000000000E0_DP
 S_E%A_Y(2,229)=-0.500000000000000E0_DP
 S_E%B_Y(2,230)=1.00000000000000E0_DP
 S_E%A_X(2,244)=3.472222222222226E-002_DP
 S_E%B_Y(2,244)=-0.273437500000000E0_DP
 S_E%B_X(2,245)=-0.312500000000000E0_DP
 S_E%A_Y(2,245)=-3.906249999999997E-002_DP
 S_E%A_X(2,246)=-5.357142857142853E-002_DP
 S_E%B_Y(2,246)=0.312500000000000E0_DP
 S_E%B_X(2,247)=0.375000000000000E0_DP
 S_E%A_Y(2,247)=6.250000000000000E-002_DP
 S_E%A_X(2,248)=0.100000000000000E0_DP
 S_E%B_Y(2,248)=-0.375000000000000E0_DP
 S_E%B_X(2,249)=-0.500000000000000E0_DP
 S_E%A_Y(2,249)=-0.125000000000000E0_DP
 S_E%A_X(2,250)=-0.333333333333333E0_DP
 S_E%B_Y(2,250)=0.500000000000000E0_DP
 S_E%B_X(2,251)=1.00000000000000E0_DP
 S_E%A_Y(2,251)=0.500000000000000E0_DP
 S_E%B_Y(2,252)=-1.00000000000000E0_DP
 S_E%B_X(2,253)=1.00000000000000E0_DP
 S_E%A_Y(2,253)=1.00000000000000E0_DP
 S_E%B_X(2,266)=-2.734374999999996E-002_DP
 S_E%A_Y(2,266)=-3.038194444444450E-003_DP
 S_E%A_X(2,267)=-4.340277777777775E-003_DP
 S_E%B_Y(2,267)=3.038194444444442E-002_DP
 S_E%B_X(2,268)=3.906250000000003E-002_DP
 S_E%A_Y(2,268)=5.580357142857133E-003_DP
 S_E%A_X(2,269)=8.928571428571428E-003_DP
 S_E%B_Y(2,269)=-4.464285714285711E-002_DP
 S_E%B_X(2,270)=-6.249999999999998E-002_DP
 S_E%A_Y(2,270)=-1.249999999999999E-002_DP
 S_E%A_X(2,271)=-2.500000000000000E-002_DP
 S_E%B_Y(2,271)=7.500000000000002E-002_DP
 S_E%B_X(2,272)=0.125000000000000E0_DP
 S_E%A_Y(2,272)=4.166666666666667E-002_DP
 S_E%A_X(2,273)=0.166666666666667E0_DP
 S_E%B_Y(2,273)=-0.166666666666667E0_DP
 S_E%B_X(2,274)=-0.500000000000000E0_DP
 S_E%A_Y(2,274)=-0.500000000000000E0_DP
 S_E%A_X(2,275)=1.00000000000000E0_DP
 S_E%B_Y(2,275)=-1.00000000000000E0_DP
 S_E%B_Y(3,104)=1.00000000000000E0_DP
 S_E%B_X(3,118)=4.50000000000000E0_DP
 S_E%A_Y(3,118)=1.00000000000000E0_DP
 S_E%B_Y(3,119)=-1.00000000000000E0_DP
 S_E%A_X(3,133)=2.66666666666667E0_DP
 S_E%B_Y(3,133)=-6.50000000000000E0_DP
 S_E%B_X(3,134)=-4.00000000000000E0_DP
 S_E%A_Y(3,134)=-1.00000000000000E0_DP
 S_E%B_Y(3,135)=1.00000000000000E0_DP
 S_E%B_X(3,149)=-11.3750000000000E0_DP
 S_E%A_Y(3,149)=-2.33333333333333E0_DP
 S_E%A_X(3,150)=-2.33333333333333E0_DP
 S_E%B_Y(3,150)=5.16666666666666E0_DP
 S_E%B_X(3,151)=3.50000000000000E0_DP
 S_E%A_Y(3,151)=1.00000000000000E0_DP
 S_E%B_Y(3,152)=-1.00000000000000E0_DP
 S_E%A_X(3,166)=-2.80000000000000E0_DP
 S_E%B_Y(3,166)=9.97499999999999E0_DP
 S_E%B_X(3,167)=7.75000000000000E0_DP
 S_E%A_Y(3,167)=1.75000000000000E0_DP
 S_E%A_X(3,168)=2.00000000000000E0_DP
 S_E%B_Y(3,168)=-4.00000000000000E0_DP
 S_E%B_X(3,169)=-3.00000000000000E0_DP
 S_E%A_Y(3,169)=-1.00000000000000E0_DP
 S_E%B_Y(3,170)=1.00000000000000E0_DP
 S_E%B_X(3,184)=8.31250000000000E0_DP
 S_E%A_Y(3,184)=1.75000000000000E0_DP
 S_E%A_X(3,185)=1.75000000000000E0_DP
 S_E%B_Y(3,185)=-5.62499999999999E0_DP
 S_E%B_X(3,186)=-5.00000000000000E0_DP
 S_E%A_Y(3,186)=-1.25000000000000E0_DP
 S_E%A_X(3,187)=-1.66666666666667E0_DP
 S_E%B_Y(3,187)=3.00000000000000E0_DP
 S_E%B_X(3,188)=2.50000000000000E0_DP
 S_E%A_Y(3,188)=1.00000000000000E0_DP
 S_E%B_Y(3,189)=-1.00000000000000E0_DP
 S_E%A_X(3,203)=1.00000000000000E0_DP
 S_E%B_Y(3,203)=-3.93750000000000E0_DP
 S_E%B_X(3,204)=-3.75000000000000E0_DP
 S_E%A_Y(3,204)=-0.875000000000000E0_DP
 S_E%A_X(3,205)=-1.00000000000000E0_DP
 S_E%B_Y(3,205)=2.87500000000000E0_DP
 S_E%B_X(3,206)=3.00000000000000E0_DP
 S_E%A_Y(3,206)=0.833333333333333E0_DP
 S_E%A_X(3,207)=1.33333333333333E0_DP
 S_E%B_Y(3,207)=-2.16666666666667E0_DP
 S_E%B_X(3,208)=-2.00000000000000E0_DP
 S_E%A_Y(3,208)=-1.00000000000000E0_DP
 S_E%B_Y(3,209)=1.00000000000000E0_DP
 S_E%B_X(3,223)=-1.47656250000000E0_DP
 S_E%A_Y(3,223)=-0.312500000000000E0_DP
 S_E%A_X(3,224)=-0.375000000000000E0_DP
 S_E%B_Y(3,224)=1.33035714285714E0_DP
 S_E%B_X(3,225)=1.43750000000000E0_DP
 S_E%A_Y(3,225)=0.375000000000000E0_DP
 S_E%A_X(3,226)=0.500000000000000E0_DP
 S_E%B_Y(3,226)=-1.27500000000000E0_DP
 S_E%B_X(3,227)=-1.62500000000000E0_DP
 S_E%A_Y(3,227)=-0.500000000000000E0_DP
 S_E%A_X(3,228)=-1.00000000000000E0_DP
 S_E%B_Y(3,228)=1.50000000000000E0_DP
 S_E%B_X(3,229)=1.50000000000000E0_DP
 S_E%A_Y(3,229)=1.00000000000000E0_DP
 S_E%B_Y(3,230)=-1.00000000000000E0_DP
 S_E%B_X(3,231)=1.00000000000000E0_DP
 S_E%A_Y(3,231)=1.00000000000000E0_DP
 S_E%A_X(3,244)=-6.944444444444439E-002_DP
 S_E%B_Y(3,244)=0.286458333333333E0_DP
 S_E%B_X(3,245)=0.332589285714285E0_DP
 S_E%A_Y(3,245)=7.812500000000000E-002_DP
 S_E%A_X(3,246)=0.107142857142857E0_DP
 S_E%B_Y(3,246)=-0.339285714285715E0_DP
 S_E%B_X(3,247)=-0.425000000000000E0_DP
 S_E%A_Y(3,247)=-0.125000000000000E0_DP
 S_E%A_X(3,248)=-0.200000000000000E0_DP
 S_E%B_Y(3,248)=0.450000000000000E0_DP
 S_E%B_X(3,249)=0.750000000000000E0_DP
 S_E%A_Y(3,249)=0.250000000000000E0_DP
 S_E%A_X(3,250)=0.666666666666667E0_DP
 S_E%B_Y(3,250)=-1.00000000000000E0_DP
 S_E%B_X(3,251)=-1.00000000000000E0_DP
 S_E%A_Y(3,251)=-1.00000000000000E0_DP
 S_E%A_X(3,252)=2.00000000000000E0_DP
 S_E%B_Y(3,252)=-2.00000000000000E0_DP
 S_E%B_X(3,266)=2.864583333333333E-002_DP
 S_E%A_Y(3,266)=6.076388888888886E-003_DP
 S_E%A_X(3,267)=8.680555555555556E-003_DP
 S_E%B_Y(3,267)=-3.224206349206343E-002_DP
 S_E%B_X(3,268)=-4.241071428571433E-002_DP
 S_E%A_Y(3,268)=-1.116071428571428E-002_DP
 S_E%A_X(3,269)=-1.785714285714286E-002_DP
 S_E%B_Y(3,269)=5.000000000000000E-002_DP
 S_E%B_X(3,270)=7.499999999999998E-002_DP
 S_E%A_Y(3,270)=2.500000000000000E-002_DP
 S_E%A_X(3,271)=5.000000000000000E-002_DP
 S_E%B_Y(3,271)=-0.100000000000000E0_DP
 S_E%B_X(3,272)=-0.250000000000000E0_DP
 S_E%A_Y(3,272)=-8.333333333333334E-002_DP
 S_E%A_X(3,273)=-0.333333333333333E0_DP
 S_E%B_Y(3,273)=0.666666666666667E0_DP
 S_E%B_X(3,274)=-1.00000000000000E0_DP
 S_E%A_Y(3,274)=-1.00000000000000E0_DP
 S_E%B_Y(4,104)=-1.00000000000000E0_DP
 S_E%B_X(4,118)=-4.50000000000000E0_DP
 S_E%A_Y(4,118)=-1.50000000000000E0_DP
 S_E%B_Y(4,119)=1.00000000000000E0_DP
 S_E%A_X(4,133)=-4.00000000000000E0_DP
 S_E%B_Y(4,133)=7.50000000000000E0_DP
 S_E%B_X(4,134)=4.00000000000000E0_DP
 S_E%A_Y(4,134)=1.50000000000000E0_DP
 S_E%B_Y(4,135)=-1.00000000000000E0_DP
 S_E%B_X(4,149)=13.1250000000000E0_DP
 S_E%A_Y(4,149)=3.87500000000000E0_DP
 S_E%A_X(4,150)=3.50000000000000E0_DP
 S_E%B_Y(4,150)=-6.16666666666667E0_DP
 S_E%B_X(4,151)=-3.50000000000000E0_DP
 S_E%A_Y(4,151)=-1.50000000000000E0_DP
 S_E%B_Y(4,152)=1.00000000000000E0_DP
 S_E%A_X(4,166)=4.65000000000000E0_DP
 S_E%B_Y(4,166)=-11.0250000000000E0_DP
 S_E%B_X(4,167)=-9.25000000000001E0_DP
 S_E%A_Y(4,167)=-3.00000000000000E0_DP
 S_E%A_X(4,168)=-3.00000000000000E0_DP
 S_E%B_Y(4,168)=5.00000000000000E0_DP
 S_E%B_X(4,169)=3.00000000000000E0_DP
 S_E%A_Y(4,169)=1.50000000000000E0_DP
 S_E%B_Y(4,170)=-1.00000000000000E0_DP
 S_E%B_X(4,184)=-9.18749999999999E0_DP
 S_E%A_Y(4,184)=-2.81250000000000E0_DP
 S_E%A_X(4,185)=-3.00000000000000E0_DP
 S_E%B_Y(4,185)=6.37500000000001E0_DP
 S_E%B_X(4,186)=6.25000000000000E0_DP
 S_E%A_Y(4,186)=2.25000000000000E0_DP
 S_E%A_X(4,187)=2.50000000000000E0_DP
 S_E%B_Y(4,187)=-4.00000000000000E0_DP
 S_E%B_X(4,188)=-2.50000000000000E0_DP
 S_E%A_Y(4,188)=-1.50000000000000E0_DP
 S_E%B_Y(4,189)=1.00000000000000E0_DP
 S_E%A_X(4,203)=-1.60714285714286E0_DP
 S_E%B_Y(4,203)=4.31250000000000E0_DP
 S_E%B_X(4,204)=4.25000000000000E0_DP
 S_E%A_Y(4,204)=1.43750000000000E0_DP
 S_E%A_X(4,205)=1.80000000000000E0_DP
 S_E%B_Y(4,205)=-3.37500000000000E0_DP
 S_E%B_X(4,206)=-4.00000000000000E0_DP
 S_E%A_Y(4,206)=-1.62500000000000E0_DP
 S_E%A_X(4,207)=-2.00000000000000E0_DP
 S_E%B_Y(4,207)=3.16666666666667E0_DP
 S_E%B_X(4,208)=2.00000000000000E0_DP
 S_E%A_Y(4,208)=1.50000000000000E0_DP
 S_E%B_Y(4,209)=-1.00000000000000E0_DP
 S_E%B_X(4,210)=1.00000000000000E0_DP
 S_E%A_Y(4,210)=1.00000000000000E0_DP
 S_E%B_X(4,223)=1.61718750000000E0_DP
 S_E%A_Y(4,223)=0.498883928571429E0_DP
 S_E%A_X(4,224)=0.616071428571429E0_DP
 S_E%B_Y(4,224)=-1.49107142857143E0_DP
 S_E%B_X(4,225)=-1.68750000000000E0_DP
 S_E%A_Y(4,225)=-0.637500000000000E0_DP
 S_E%A_X(4,226)=-0.975000000000000E0_DP
 S_E%B_Y(4,226)=1.57500000000000E0_DP
 S_E%B_X(4,227)=2.37500000000000E0_DP
 S_E%A_Y(4,227)=1.12500000000000E0_DP
 S_E%A_X(4,228)=1.50000000000000E0_DP
 S_E%B_Y(4,228)=-2.50000000000000E0_DP
 S_E%B_X(4,229)=-1.50000000000000E0_DP
 S_E%A_Y(4,229)=-1.50000000000000E0_DP
 S_E%A_X(4,230)=3.00000000000000E0_DP
 S_E%B_Y(4,230)=-3.00000000000000E0_DP
 S_E%A_X(4,244)=0.110863095238095E0_DP
 S_E%B_Y(4,244)=-0.312500000000000E0_DP
 S_E%B_X(4,245)=-0.372767857142858E0_DP
 S_E%A_Y(4,245)=-0.127232142857143E0_DP
 S_E%A_X(4,246)=-0.182142857142857E0_DP
 S_E%B_Y(4,246)=0.392857142857143E0_DP
 S_E%B_X(4,247)=0.525000000000000E0_DP
 S_E%A_Y(4,247)=0.225000000000000E0_DP
 S_E%A_X(4,248)=0.450000000000000E0_DP
 S_E%B_Y(4,248)=-0.600000000000000E0_DP
 S_E%B_X(4,249)=-1.25000000000000E0_DP
 S_E%A_Y(4,249)=-0.750000000000000E0_DP
 S_E%A_X(4,250)=-1.00000000000000E0_DP
 S_E%B_Y(4,250)=2.00000000000000E0_DP
 S_E%B_X(4,251)=-3.00000000000000E0_DP
 S_E%A_Y(4,251)=-3.00000000000000E0_DP
 S_E%B_X(4,266)=-3.124999999999997E-002_DP
 S_E%A_Y(4,266)=-9.672619047619055E-003_DP
 S_E%A_X(4,267)=-1.413690476190477E-002_DP
 S_E%B_Y(4,267)=3.596230158730165E-002_DP
 S_E%B_X(4,268)=4.910714285714287E-002_DP
 S_E%A_Y(4,268)=1.874999999999999E-002_DP
 S_E%A_X(4,269)=3.214285714285715E-002_DP
 S_E%B_Y(4,269)=-6.071428571428574E-002_DP
 S_E%B_X(4,270)=-0.100000000000000E0_DP
 S_E%A_Y(4,270)=-4.999999999999999E-002_DP
 S_E%A_X(4,271)=-0.150000000000000E0_DP
 S_E%B_Y(4,271)=0.150000000000000E0_DP
 S_E%B_X(4,272)=0.500000000000000E0_DP
 S_E%A_Y(4,272)=0.500000000000000E0_DP
 S_E%A_X(4,273)=-1.00000000000000E0_DP
 S_E%B_Y(4,273)=1.00000000000000E0_DP
 S_E%B_Y(5,104)=1.00000000000000E0_DP
 S_E%B_X(5,118)=4.50000000000000E0_DP
 S_E%A_Y(5,118)=2.00000000000000E0_DP
 S_E%B_Y(5,119)=-1.00000000000000E0_DP
 S_E%A_X(5,133)=5.33333333333333E0_DP
 S_E%B_Y(5,133)=-9.00000000000000E0_DP
 S_E%B_X(5,134)=-4.00000000000000E0_DP
 S_E%A_Y(5,134)=-2.00000000000000E0_DP
 S_E%B_Y(5,135)=1.00000000000000E0_DP
 S_E%B_X(5,149)=-15.7500000000000E0_DP
 S_E%A_Y(5,149)=-6.16666666666667E0_DP
 S_E%A_X(5,150)=-4.66666666666667E0_DP
 S_E%B_Y(5,150)=7.66666666666667E0_DP
 S_E%B_X(5,151)=3.50000000000000E0_DP
 S_E%A_Y(5,151)=2.00000000000000E0_DP
 S_E%B_Y(5,152)=-1.00000000000000E0_DP
 S_E%A_X(5,166)=-7.40000000000000E0_DP
 S_E%B_Y(5,166)=12.9750000000000E0_DP
 S_E%B_X(5,167)=11.5000000000000E0_DP
 S_E%A_Y(5,167)=5.00000000000000E0_DP
 S_E%A_X(5,168)=4.00000000000000E0_DP
 S_E%B_Y(5,168)=-6.50000000000000E0_DP
 S_E%B_X(5,169)=-3.00000000000000E0_DP
 S_E%A_Y(5,169)=-2.00000000000000E0_DP
 S_E%B_Y(5,170)=1.00000000000000E0_DP
 S_E%B_X(5,184)=10.8125000000000E0_DP
 S_E%A_Y(5,184)=4.25000000000000E0_DP
 S_E%A_X(5,185)=5.00000000000000E0_DP
 S_E%B_Y(5,185)=-7.87500000000000E0_DP
 S_E%B_X(5,186)=-8.12500000000000E0_DP
 S_E%A_Y(5,186)=-4.00000000000000E0_DP
 S_E%A_X(5,187)=-3.33333333333333E0_DP
 S_E%B_Y(5,187)=5.50000000000000E0_DP
 S_E%B_X(5,188)=2.50000000000000E0_DP
 S_E%A_Y(5,188)=2.00000000000000E0_DP
 S_E%B_Y(5,189)=-1.00000000000000E0_DP
 S_E%B_X(5,190)=1.00000000000000E0_DP
 S_E%A_Y(5,190)=1.00000000000000E0_DP
 S_E%A_X(5,203)=2.42857142857143E0_DP
 S_E%B_Y(5,203)=-4.96428571428572E0_DP
 S_E%B_X(5,204)=-5.25000000000000E0_DP
 S_E%A_Y(5,204)=-2.25000000000000E0_DP
 S_E%A_X(5,205)=-3.20000000000000E0_DP
 S_E%B_Y(5,205)=4.50000000000000E0_DP
 S_E%B_X(5,206)=5.50000000000000E0_DP
 S_E%A_Y(5,206)=3.16666666666667E0_DP
 S_E%A_X(5,207)=2.66666666666667E0_DP
 S_E%B_Y(5,207)=-4.66666666666667E0_DP
 S_E%B_X(5,208)=-2.00000000000000E0_DP
 S_E%A_Y(5,208)=-2.00000000000000E0_DP
 S_E%A_X(5,209)=4.00000000000000E0_DP
 S_E%B_Y(5,209)=-4.00000000000000E0_DP
 S_E%B_X(5,223)=-1.86160714285714E0_DP
 S_E%A_Y(5,223)=-0.745535714285715E0_DP
 S_E%A_X(5,224)=-0.964285714285715E0_DP
 S_E%B_Y(5,224)=1.78571428571429E0_DP
 S_E%B_X(5,225)=2.25000000000000E0_DP
 S_E%A_Y(5,225)=1.05000000000000E0_DP
 S_E%A_X(5,226)=1.90000000000000E0_DP
 S_E%B_Y(5,226)=-2.40000000000000E0_DP
 S_E%B_X(5,227)=-3.50000000000000E0_DP
 S_E%A_Y(5,227)=-2.50000000000000E0_DP
 S_E%A_X(5,228)=-2.00000000000000E0_DP
 S_E%B_Y(5,228)=4.00000000000000E0_DP
 S_E%B_X(5,229)=-6.00000000000000E0_DP
 S_E%A_Y(5,229)=-6.00000000000000E0_DP
 S_E%A_X(5,244)=-0.165674603174603E0_DP
 S_E%B_Y(5,244)=0.357142857142857E0_DP
 S_E%B_X(5,245)=0.446428571428572E0_DP
 S_E%A_Y(5,245)=0.196428571428571E0_DP
 S_E%A_X(5,246)=0.300000000000000E0_DP
 S_E%B_Y(5,246)=-0.500000000000000E0_DP
 S_E%B_X(5,247)=-0.800000000000000E0_DP
 S_E%A_Y(5,247)=-0.400000000000000E0_DP
 S_E%A_X(5,248)=-1.00000000000000E0_DP
 S_E%B_Y(5,248)=1.20000000000000E0_DP
 S_E%B_X(5,249)=2.00000000000000E0_DP
 S_E%A_Y(5,249)=2.00000000000000E0_DP
 S_E%A_X(5,250)=-4.00000000000000E0_DP
 S_E%B_Y(5,250)=4.00000000000000E0_DP
 S_E%B_X(5,266)=3.571428571428573E-002_DP
 S_E%A_Y(5,266)=1.438492063492064E-002_DP
 S_E%A_X(5,267)=2.182539682539683E-002_DP
 S_E%B_Y(5,267)=-4.265873015873016E-002_DP
 S_E%B_X(5,268)=-6.249999999999996E-002_DP
 S_E%A_Y(5,268)=-3.035714285714283E-002_DP
 S_E%A_X(5,269)=-5.714285714285713E-002_DP
 S_E%B_Y(5,269)=8.571428571428572E-002_DP
 S_E%B_X(5,270)=0.200000000000000E0_DP
 S_E%A_Y(5,270)=0.100000000000000E0_DP
 S_E%A_X(5,271)=0.400000000000000E0_DP
 S_E%B_Y(5,271)=-0.600000000000000E0_DP
 S_E%B_X(5,272)=1.00000000000000E0_DP
 S_E%A_Y(5,272)=1.00000000000000E0_DP
 S_E%B_Y(6,104)=-1.00000000000000E0_DP
 S_E%B_X(6,118)=-4.50000000000000E0_DP
 S_E%A_Y(6,118)=-2.50000000000000E0_DP
 S_E%B_Y(6,119)=1.00000000000000E0_DP
 S_E%A_X(6,133)=-6.66666666666667E0_DP
 S_E%B_Y(6,133)=11.0000000000000E0_DP
 S_E%B_X(6,134)=4.00000000000000E0_DP
 S_E%A_Y(6,134)=2.50000000000000E0_DP
 S_E%B_Y(6,135)=-1.00000000000000E0_DP
 S_E%B_X(6,149)=19.2500000000000E0_DP
 S_E%A_Y(6,149)=9.58333333333333E0_DP
 S_E%A_X(6,150)=5.83333333333333E0_DP
 S_E%B_Y(6,150)=-9.66666666666667E0_DP
 S_E%B_X(6,151)=-3.50000000000000E0_DP
 S_E%A_Y(6,151)=-2.50000000000000E0_DP
 S_E%B_Y(6,152)=1.00000000000000E0_DP
 S_E%A_X(6,166)=11.5000000000000E0_DP
 S_E%B_Y(6,166)=-16.5750000000000E0_DP
 S_E%B_X(6,167)=-14.5000000000000E0_DP
 S_E%A_Y(6,167)=-8.12500000000000E0_DP
 S_E%A_X(6,168)=-5.00000000000000E0_DP
 S_E%B_Y(6,168)=8.50000000000000E0_DP
 S_E%B_X(6,169)=3.00000000000000E0_DP
 S_E%A_Y(6,169)=2.50000000000000E0_DP
 S_E%B_Y(6,170)=-1.00000000000000E0_DP
 S_E%B_X(6,171)=1.00000000000000E0_DP
 S_E%A_Y(6,171)=1.00000000000000E0_DP
 S_E%B_X(6,184)=-13.8125000000000E0_DP
 S_E%A_Y(6,184)=-6.56250000000000E0_DP
 S_E%A_X(6,185)=-8.12500000000000E0_DP
 S_E%B_Y(6,185)=10.8750000000000E0_DP
 S_E%B_X(6,186)=10.6250000000000E0_DP
 S_E%A_Y(6,186)=6.87500000000000E0_DP
 S_E%A_X(6,187)=4.16666666666667E0_DP
 S_E%B_Y(6,187)=-7.50000000000000E0_DP
 S_E%B_X(6,188)=-2.50000000000000E0_DP
 S_E%A_Y(6,188)=-2.50000000000000E0_DP
 S_E%A_X(6,189)=5.00000000000000E0_DP
 S_E%B_Y(6,189)=-5.00000000000000E0_DP
 S_E%A_X(6,203)=-3.75000000000000E0_DP
 S_E%B_Y(6,203)=6.07142857142857E0_DP
 S_E%B_X(6,204)=7.25000000000000E0_DP
 S_E%A_Y(6,204)=3.75000000000000E0_DP
 S_E%A_X(6,205)=5.50000000000000E0_DP
 S_E%B_Y(6,205)=-7.00000000000000E0_DP
 S_E%B_X(6,206)=-7.50000000000000E0_DP
 S_E%A_Y(6,206)=-5.83333333333333E0_DP
 S_E%A_X(6,207)=-3.33333333333333E0_DP
 S_E%B_Y(6,207)=6.66666666666667E0_DP
 S_E%B_X(6,208)=-10.0000000000000E0_DP
 S_E%A_Y(6,208)=-10.0000000000000E0_DP
 S_E%B_X(6,223)=2.27678571428571E0_DP
 S_E%A_Y(6,223)=1.11607142857143E0_DP
 S_E%A_X(6,224)=1.60714285714286E0_DP
 S_E%B_Y(6,224)=-2.32142857142857E0_DP
 S_E%B_X(6,225)=-3.50000000000000E0_DP
 S_E%A_Y(6,225)=-2.00000000000000E0_DP
 S_E%A_X(6,226)=-3.50000000000000E0_DP
 S_E%B_Y(6,226)=4.50000000000000E0_DP
 S_E%B_X(6,227)=5.00000000000000E0_DP
 S_E%A_Y(6,227)=5.00000000000000E0_DP
 S_E%A_X(6,228)=-10.0000000000000E0_DP
 S_E%B_Y(6,228)=10.0000000000000E0_DP
 S_E%A_X(6,244)=0.248015873015873E0_DP
 S_E%B_Y(6,244)=-0.431547619047619E0_DP
 S_E%B_X(6,245)=-0.580357142857143E0_DP
 S_E%A_Y(6,245)=-0.312500000000000E0_DP
 S_E%A_X(6,246)=-0.571428571428571E0_DP
 S_E%B_Y(6,246)=0.714285714285714E0_DP
 S_E%B_X(6,247)=1.50000000000000E0_DP
 S_E%A_Y(6,247)=1.00000000000000E0_DP
 S_E%A_X(6,248)=2.00000000000000E0_DP
 S_E%B_Y(6,248)=-3.00000000000000E0_DP
 S_E%B_X(6,249)=5.00000000000000E0_DP
 S_E%A_Y(6,249)=5.00000000000000E0_DP
 S_E%B_X(6,266)=-4.315476190476190E-002_DP
 S_E%A_Y(6,266)=-2.132936507936508E-002_DP
 S_E%A_X(6,267)=-3.472222222222222E-002_DP
 S_E%B_Y(6,267)=5.456349206349211E-002_DP
 S_E%B_X(6,268)=8.928571428571426E-002_DP
 S_E%A_Y(6,268)=5.357142857142858E-002_DP
 S_E%A_X(6,269)=0.142857142857143E0_DP
 S_E%B_Y(6,269)=-0.142857142857143E0_DP
 S_E%B_X(6,270)=-0.500000000000000E0_DP
 S_E%A_Y(6,270)=-0.500000000000000E0_DP
 S_E%A_X(6,271)=1.00000000000000E0_DP
 S_E%B_Y(6,271)=-1.00000000000000E0_DP
 S_E%B_Y(7,104)=1.00000000000000E0_DP
 S_E%B_X(7,118)=4.50000000000000E0_DP
 S_E%A_Y(7,118)=3.00000000000000E0_DP
 S_E%B_Y(7,119)=-1.00000000000000E0_DP
 S_E%A_X(7,133)=8.00000000000000E0_DP
 S_E%B_Y(7,133)=-13.5000000000000E0_DP
 S_E%B_X(7,134)=-4.00000000000000E0_DP
 S_E%A_Y(7,134)=-3.00000000000000E0_DP
 S_E%B_Y(7,135)=1.00000000000000E0_DP
 S_E%B_X(7,149)=-23.6250000000000E0_DP
 S_E%A_Y(7,149)=-14.5000000000000E0_DP
 S_E%A_X(7,150)=-7.00000000000000E0_DP
 S_E%B_Y(7,150)=12.1666666666667E0_DP
 S_E%B_X(7,151)=3.50000000000000E0_DP
 S_E%A_Y(7,151)=3.00000000000000E0_DP
 S_E%B_Y(7,152)=-1.00000000000000E0_DP
 S_E%B_X(7,153)=1.00000000000000E0_DP
 S_E%A_Y(7,153)=1.00000000000000E0_DP
 S_E%A_X(7,166)=-17.4000000000000E0_DP
 S_E%B_Y(7,166)=22.9500000000000E0_DP
 S_E%B_X(7,167)=18.2500000000000E0_DP
 S_E%A_Y(7,167)=12.7500000000000E0_DP
 S_E%A_X(7,168)=6.00000000000000E0_DP
 S_E%B_Y(7,168)=-11.0000000000000E0_DP
 S_E%B_X(7,169)=-3.00000000000000E0_DP
 S_E%A_Y(7,169)=-3.00000000000000E0_DP
 S_E%A_X(7,170)=6.00000000000000E0_DP
 S_E%B_Y(7,170)=-6.00000000000000E0_DP
 S_E%B_X(7,184)=19.1250000000000E0_DP
 S_E%A_Y(7,184)=10.8750000000000E0_DP
 S_E%A_X(7,185)=12.7500000000000E0_DP
 S_E%B_Y(7,185)=-16.5000000000000E0_DP
 S_E%B_X(7,186)=-13.7500000000000E0_DP
 S_E%A_Y(7,186)=-11.2500000000000E0_DP
 S_E%A_X(7,187)=-5.00000000000000E0_DP
 S_E%B_Y(7,187)=10.0000000000000E0_DP
 S_E%B_X(7,188)=-15.0000000000000E0_DP
 S_E%A_Y(7,188)=-15.0000000000000E0_DP
 S_E%A_X(7,203)=6.21428571428572E0_DP
 S_E%B_Y(7,203)=-8.21428571428572E0_DP
 S_E%B_X(7,204)=-11.0000000000000E0_DP
 S_E%A_Y(7,204)=-7.00000000000000E0_DP
 S_E%A_X(7,205)=-9.00000000000000E0_DP
 S_E%B_Y(7,205)=12.0000000000000E0_DP
 S_E%B_X(7,206)=10.0000000000000E0_DP
 S_E%A_Y(7,206)=10.0000000000000E0_DP
 S_E%A_X(7,207)=-20.0000000000000E0_DP
 S_E%B_Y(7,207)=20.0000000000000E0_DP
 S_E%B_X(7,223)=-3.08035714285714E0_DP
 S_E%A_Y(7,223)=-1.74107142857143E0_DP
 S_E%A_X(7,224)=-3.00000000000000E0_DP
 S_E%B_Y(7,224)=3.57142857142857E0_DP
 S_E%B_X(7,225)=6.00000000000000E0_DP
 S_E%A_Y(7,225)=4.50000000000000E0_DP
 S_E%A_X(7,226)=6.00000000000000E0_DP
 S_E%B_Y(7,226)=-9.00000000000000E0_DP
 S_E%B_X(7,227)=15.0000000000000E0_DP
 S_E%A_Y(7,227)=15.0000000000000E0_DP
 S_E%A_X(7,244)=-0.386904761904762E0_DP
 S_E%B_Y(7,244)=0.565476190476191E0_DP
 S_E%B_X(7,245)=0.892857142857142E0_DP
 S_E%A_Y(7,245)=0.535714285714286E0_DP
 S_E%A_X(7,246)=1.28571428571429E0_DP
 S_E%B_Y(7,246)=-1.42857142857143E0_DP
 S_E%B_X(7,247)=-3.00000000000000E0_DP
 S_E%A_Y(7,247)=-3.00000000000000E0_DP
 S_E%A_X(7,248)=6.00000000000000E0_DP
 S_E%B_Y(7,248)=-6.00000000000000E0_DP
 S_E%B_X(7,266)=5.654761904761905E-002_DP
 S_E%A_Y(7,266)=3.273809523809525E-002_DP
 S_E%A_X(7,267)=5.952380952380951E-002_DP
 S_E%B_Y(7,267)=-7.936507936507933E-002_DP
 S_E%B_X(7,268)=-0.178571428571429E0_DP
 S_E%A_Y(7,268)=-0.107142857142857E0_DP
 S_E%A_X(7,269)=-0.428571428571429E0_DP
 S_E%B_Y(7,269)=0.571428571428571E0_DP
 S_E%B_X(7,270)=-1.00000000000000E0_DP
 S_E%A_Y(7,270)=-1.00000000000000E0_DP
 S_E%B_Y(8,104)=-1.00000000000000E0_DP
 S_E%B_X(8,118)=-4.50000000000000E0_DP
 S_E%A_Y(8,118)=-3.50000000000000E0_DP
 S_E%B_Y(8,119)=1.00000000000000E0_DP
 S_E%A_X(8,133)=-9.33333333333334E0_DP
 S_E%B_Y(8,133)=16.5000000000000E0_DP
 S_E%B_X(8,134)=4.00000000000000E0_DP
 S_E%A_Y(8,134)=3.50000000000000E0_DP
 S_E%B_Y(8,135)=-1.00000000000000E0_DP
 S_E%B_X(8,136)=1.00000000000000E0_DP
 S_E%A_Y(8,136)=1.00000000000000E0_DP
 S_E%B_X(8,149)=28.8750000000000E0_DP
 S_E%A_Y(8,149)=21.2916666666667E0_DP
 S_E%A_X(8,150)=8.16666666666667E0_DP
 S_E%B_Y(8,150)=-15.1666666666667E0_DP
 S_E%B_X(8,151)=-3.50000000000000E0_DP
 S_E%A_Y(8,151)=-3.50000000000000E0_DP
 S_E%A_X(8,152)=7.00000000000000E0_DP
 S_E%B_Y(8,152)=-7.00000000000000E0_DP
 S_E%A_X(8,166)=25.5500000000000E0_DP
 S_E%B_Y(8,166)=-33.6000000000000E0_DP
 S_E%B_X(8,167)=-22.7500000000000E0_DP
 S_E%A_Y(8,167)=-19.2500000000000E0_DP
 S_E%A_X(8,168)=-7.00000000000000E0_DP
 S_E%B_Y(8,168)=14.0000000000000E0_DP
 S_E%B_X(8,169)=-21.0000000000000E0_DP
 S_E%A_Y(8,169)=-21.0000000000000E0_DP
 S_E%B_X(8,184)=-28.0000000000000E0_DP
 S_E%A_Y(8,184)=-19.2500000000000E0_DP
 S_E%A_X(8,185)=-19.2500000000000E0_DP
 S_E%B_Y(8,185)=26.2500000000000E0_DP
 S_E%B_X(8,186)=17.5000000000000E0_DP
 S_E%A_Y(8,186)=17.5000000000000E0_DP
 S_E%A_X(8,187)=-35.0000000000000E0_DP
 S_E%B_Y(8,187)=35.0000000000000E0_DP
 S_E%A_X(8,203)=-11.0000000000000E0_DP
 S_E%B_Y(8,203)=13.0000000000000E0_DP
 S_E%B_X(8,204)=17.5000000000000E0_DP
 S_E%A_Y(8,204)=14.0000000000000E0_DP
 S_E%A_X(8,205)=14.0000000000000E0_DP
 S_E%B_Y(8,205)=-21.0000000000000E0_DP
 S_E%B_X(8,206)=35.0000000000000E0_DP
 S_E%A_Y(8,206)=35.0000000000000E0_DP
 S_E%B_X(8,223)=4.87500000000000E0_DP
 S_E%A_Y(8,223)=3.12500000000000E0_DP
 S_E%A_X(8,224)=6.00000000000000E0_DP
 S_E%B_Y(8,224)=-7.00000000000000E0_DP
 S_E%B_X(8,225)=-10.5000000000000E0_DP
 S_E%A_Y(8,225)=-10.5000000000000E0_DP
 S_E%A_X(8,226)=21.0000000000000E0_DP
 S_E%B_Y(8,226)=-21.0000000000000E0_DP
 S_E%A_X(8,244)=0.694444444444445E0_DP
 S_E%B_Y(8,244)=-0.833333333333333E0_DP
 S_E%B_X(8,245)=-1.75000000000000E0_DP
 S_E%A_Y(8,245)=-1.25000000000000E0_DP
 S_E%A_X(8,246)=-3.00000000000000E0_DP
 S_E%B_Y(8,246)=4.00000000000000E0_DP
 S_E%B_X(8,247)=-7.00000000000000E0_DP
 S_E%A_Y(8,247)=-7.00000000000000E0_DP
 S_E%B_X(8,266)=-8.333333333333330E-002_DP
 S_E%A_Y(8,266)=-5.555555555555557E-002_DP
 S_E%A_X(8,267)=-0.138888888888889E0_DP
 S_E%B_Y(8,267)=0.138888888888889E0_DP
 S_E%B_X(8,268)=0.500000000000000E0_DP
 S_E%A_Y(8,268)=0.500000000000000E0_DP
 S_E%A_X(8,269)=-1.00000000000000E0_DP
 S_E%B_Y(8,269)=1.00000000000000E0_DP
 S_E%B_Y(9,104)=1.00000000000000E0_DP
 S_E%B_X(9,118)=4.50000000000000E0_DP
 S_E%A_Y(9,118)=4.00000000000000E0_DP
 S_E%B_Y(9,119)=-1.00000000000000E0_DP
 S_E%B_X(9,120)=1.00000000000000E0_DP
 S_E%A_Y(9,120)=1.00000000000000E0_DP
 S_E%A_X(9,133)=10.6666666666667E0_DP
 S_E%B_Y(9,133)=-20.0000000000000E0_DP
 S_E%B_X(9,134)=-4.00000000000000E0_DP
 S_E%A_Y(9,134)=-4.00000000000000E0_DP
 S_E%A_X(9,135)=8.00000000000000E0_DP
 S_E%B_Y(9,135)=-8.00000000000000E0_DP
 S_E%B_X(9,149)=-35.0000000000000E0_DP
 S_E%A_Y(9,149)=-30.3333333333333E0_DP
 S_E%A_X(9,150)=-9.33333333333333E0_DP
 S_E%B_Y(9,150)=18.6666666666667E0_DP
 S_E%B_X(9,151)=-28.0000000000000E0_DP
 S_E%A_Y(9,151)=-28.0000000000000E0_DP
 S_E%A_X(9,166)=-36.4000000000000E0_DP
 S_E%B_Y(9,166)=50.4000000000000E0_DP
 S_E%B_X(9,167)=28.0000000000000E0_DP
 S_E%A_Y(9,167)=28.0000000000000E0_DP
 S_E%A_X(9,168)=-56.0000000000000E0_DP
 S_E%B_Y(9,168)=56.0000000000000E0_DP
 S_E%B_X(9,184)=42.0000000000000E0_DP
 S_E%A_Y(9,184)=35.0000000000000E0_DP
 S_E%A_X(9,185)=28.0000000000000E0_DP
 S_E%B_Y(9,185)=-42.0000000000000E0_DP
 S_E%B_X(9,186)=70.0000000000000E0_DP
 S_E%A_Y(9,186)=70.0000000000000E0_DP
 S_E%A_X(9,203)=20.0000000000000E0_DP
 S_E%B_Y(9,203)=-24.0000000000000E0_DP
 S_E%B_X(9,204)=-28.0000000000000E0_DP
 S_E%A_Y(9,204)=-28.0000000000000E0_DP
 S_E%A_X(9,205)=56.0000000000000E0_DP
 S_E%B_Y(9,205)=-56.0000000000000E0_DP
 S_E%B_X(9,223)=-9.00000000000000E0_DP
 S_E%A_Y(9,223)=-7.00000000000000E0_DP
 S_E%A_X(9,224)=-12.0000000000000E0_DP
 S_E%B_Y(9,224)=16.0000000000000E0_DP
 S_E%B_X(9,225)=-28.0000000000000E0_DP
 S_E%A_Y(9,225)=-28.0000000000000E0_DP
 S_E%A_X(9,244)=-1.55555555555556E0_DP
 S_E%B_Y(9,244)=1.66666666666667E0_DP
 S_E%B_X(9,245)=4.00000000000000E0_DP
 S_E%A_Y(9,245)=4.00000000000000E0_DP
 S_E%A_X(9,246)=-8.00000000000000E0_DP
 S_E%B_Y(9,246)=8.00000000000000E0_DP
 S_E%B_X(9,266)=0.166666666666667E0_DP
 S_E%A_Y(9,266)=0.111111111111111E0_DP
 S_E%A_X(9,267)=0.444444444444444E0_DP
 S_E%B_Y(9,267)=-0.555555555555556E0_DP
 S_E%B_X(9,268)=1.00000000000000E0_DP
 S_E%A_Y(9,268)=1.00000000000000E0_DP
 S_E%B_Y(10,104)=-1.00000000000000E0_DP
 S_E%B_X(10,105)=1.00000000000000E0_DP
 S_E%A_Y(10,105)=1.00000000000000E0_DP
 S_E%B_X(10,118)=-4.50000000000000E0_DP
 S_E%A_Y(10,118)=-4.50000000000000E0_DP
 S_E%A_X(10,119)=9.00000000000000E0_DP
 S_E%B_Y(10,119)=-9.00000000000000E0_DP
 S_E%A_X(10,133)=-12.0000000000000E0_DP
 S_E%B_Y(10,133)=24.0000000000000E0_DP
 S_E%B_X(10,134)=-36.0000000000000E0_DP
 S_E%A_Y(10,134)=-36.0000000000000E0_DP
 S_E%B_X(10,149)=42.0000000000000E0_DP
 S_E%A_Y(10,149)=42.0000000000000E0_DP
 S_E%A_X(10,150)=-84.0000000000000E0_DP
 S_E%B_Y(10,150)=84.0000000000000E0_DP
 S_E%A_X(10,166)=50.4000000000000E0_DP
 S_E%B_Y(10,166)=-75.6000000000000E0_DP
 S_E%B_X(10,167)=126.000000000000E0_DP
 S_E%A_Y(10,167)=126.000000000000E0_DP
 S_E%B_X(10,184)=-63.0000000000000E0_DP
 S_E%A_Y(10,184)=-63.0000000000000E0_DP
 S_E%A_X(10,185)=126.000000000000E0_DP
 S_E%B_Y(10,185)=-126.000000000000E0_DP
 S_E%A_X(10,203)=-36.0000000000000E0_DP
 S_E%B_Y(10,203)=48.0000000000000E0_DP
 S_E%B_X(10,204)=-84.0000000000000E0_DP
 S_E%A_Y(10,204)=-84.0000000000000E0_DP
 S_E%B_X(10,223)=18.0000000000000E0_DP
 S_E%A_Y(10,223)=18.0000000000000E0_DP
 S_E%A_X(10,224)=-36.0000000000000E0_DP
 S_E%B_Y(10,224)=36.0000000000000E0_DP
 S_E%A_X(10,244)=4.00000000000000E0_DP
 S_E%B_Y(10,244)=-5.00000000000000E0_DP
 S_E%B_X(10,245)=9.00000000000000E0_DP
 S_E%A_Y(10,245)=9.00000000000000E0_DP
 S_E%B_X(10,266)=-0.500000000000000E0_DP
 S_E%A_Y(10,266)=-0.500000000000000E0_DP
 S_E%A_X(10,267)=1.00000000000000E0_DP
 S_E%B_Y(10,267)=-1.00000000000000E0_DP
 S_E%B_X(11,91)=1.00000000000000E0_DP
 S_E%A_Y(11,91)=1.00000000000000E0_DP
 S_E%A_X(11,104)=10.0000000000000E0_DP
 S_E%B_Y(11,104)=-10.0000000000000E0_DP
 S_E%B_X(11,118)=-45.0000000000000E0_DP
 S_E%A_Y(11,118)=-45.0000000000000E0_DP
 S_E%A_X(11,133)=-120.000000000000E0_DP
 S_E%B_Y(11,133)=120.000000000000E0_DP
 S_E%B_X(11,149)=210.000000000000E0_DP
 S_E%A_Y(11,149)=210.000000000000E0_DP
 S_E%A_X(11,166)=252.000000000000E0_DP
 S_E%B_Y(11,166)=-252.000000000000E0_DP
 S_E%B_X(11,184)=-210.000000000000E0_DP
 S_E%A_Y(11,184)=-210.000000000000E0_DP
 S_E%A_X(11,203)=-120.000000000000E0_DP
 S_E%B_Y(11,203)=120.000000000000E0_DP
 S_E%B_X(11,223)=45.0000000000000E0_DP
 S_E%A_Y(11,223)=45.0000000000000E0_DP
 S_E%A_X(11,244)=10.0000000000000E0_DP
 S_E%B_Y(11,244)=-10.0000000000000E0_DP
 S_E%B_X(11,266)=-1.00000000000000E0_DP
 S_E%A_Y(11,266)=-1.00000000000000E0_DP
 S_E%B_X(12,78)=1.00000000000000E0_DP
 S_E%A_Y(12,78)=1.00000000000000E0_DP
 S_E%A_X(12,90)=11.0000000000000E0_DP
 S_E%B_Y(12,90)=-11.0000000000000E0_DP
 S_E%B_X(12,103)=-55.0000000000000E0_DP
 S_E%A_Y(12,103)=-55.0000000000000E0_DP
 S_E%A_X(12,117)=-165.000000000000E0_DP
 S_E%B_Y(12,117)=165.000000000000E0_DP
 S_E%B_X(12,132)=330.000000000000E0_DP
 S_E%A_Y(12,132)=330.000000000000E0_DP
 S_E%A_X(12,148)=462.000000000000E0_DP
 S_E%B_Y(12,148)=-462.000000000000E0_DP
 S_E%B_X(12,165)=-462.000000000000E0_DP
 S_E%A_Y(12,165)=-462.000000000000E0_DP
 S_E%A_X(12,183)=-330.000000000000E0_DP
 S_E%B_Y(12,183)=330.000000000000E0_DP
 S_E%B_X(12,202)=165.000000000000E0_DP
 S_E%A_Y(12,202)=165.000000000000E0_DP
 S_E%A_X(12,222)=55.0000000000000E0_DP
 S_E%B_Y(12,222)=-55.0000000000000E0_DP
 S_E%B_X(12,243)=-11.0000000000000E0_DP
 S_E%A_Y(12,243)=-11.0000000000000E0_DP
 S_E%A_X(12,265)=-1.00000000000000E0_DP
 S_E%B_Y(12,265)=1.00000000000000E0_DP
 S_E%B_X(13,66)=1.00000000000000E0_DP
 S_E%A_Y(13,66)=1.00000000000000E0_DP
 S_E%A_X(13,77)=12.0000000000000E0_DP
 S_E%B_Y(13,77)=-12.0000000000000E0_DP
 S_E%B_X(13,89)=-66.0000000000000E0_DP
 S_E%A_Y(13,89)=-66.0000000000000E0_DP
 S_E%A_X(13,102)=-220.000000000000E0_DP
 S_E%B_Y(13,102)=220.000000000000E0_DP
 S_E%B_X(13,116)=495.000000000000E0_DP
 S_E%A_Y(13,116)=495.000000000000E0_DP
 S_E%A_X(13,131)=792.000000000000E0_DP
 S_E%B_Y(13,131)=-792.000000000000E0_DP
 S_E%B_X(13,147)=-924.000000000000E0_DP
 S_E%A_Y(13,147)=-924.000000000000E0_DP
 S_E%A_X(13,164)=-792.000000000000E0_DP
 S_E%B_Y(13,164)=792.000000000000E0_DP
 S_E%B_X(13,182)=495.000000000000E0_DP
 S_E%A_Y(13,182)=495.000000000000E0_DP
 S_E%A_X(13,201)=220.000000000000E0_DP
 S_E%B_Y(13,201)=-220.000000000000E0_DP
 S_E%B_X(13,221)=-66.0000000000000E0_DP
 S_E%A_Y(13,221)=-66.0000000000000E0_DP
 S_E%A_X(13,242)=-12.0000000000000E0_DP
 S_E%B_Y(13,242)=12.0000000000000E0_DP
 S_E%B_X(13,264)=1.00000000000000E0_DP
 S_E%A_Y(13,264)=1.00000000000000E0_DP
 S_E%B_X(14,55)=1.00000000000000E0_DP
 S_E%A_Y(14,55)=1.00000000000000E0_DP
 S_E%A_X(14,65)=13.0000000000000E0_DP
 S_E%B_Y(14,65)=-13.0000000000000E0_DP
 S_E%B_X(14,76)=-78.0000000000000E0_DP
 S_E%A_Y(14,76)=-78.0000000000000E0_DP
 S_E%A_X(14,88)=-286.000000000000E0_DP
 S_E%B_Y(14,88)=286.000000000000E0_DP
 S_E%B_X(14,101)=715.000000000000E0_DP
 S_E%A_Y(14,101)=715.000000000000E0_DP
 S_E%A_X(14,115)=1287.00000000000E0_DP
 S_E%B_Y(14,115)=-1287.00000000000E0_DP
 S_E%B_X(14,130)=-1716.00000000000E0_DP
 S_E%A_Y(14,130)=-1716.00000000000E0_DP
 S_E%A_X(14,146)=-1716.00000000000E0_DP
 S_E%B_Y(14,146)=1716.00000000000E0_DP
 S_E%B_X(14,163)=1287.00000000000E0_DP
 S_E%A_Y(14,163)=1287.00000000000E0_DP
 S_E%A_X(14,181)=715.000000000000E0_DP
 S_E%B_Y(14,181)=-715.000000000000E0_DP
 S_E%B_X(14,200)=-286.000000000000E0_DP
 S_E%A_Y(14,200)=-286.000000000000E0_DP
 S_E%A_X(14,220)=-78.0000000000000E0_DP
 S_E%B_Y(14,220)=78.0000000000000E0_DP
 S_E%B_X(14,241)=13.0000000000000E0_DP
 S_E%A_Y(14,241)=13.0000000000000E0_DP
 S_E%A_X(14,263)=1.00000000000000E0_DP
 S_E%B_Y(14,263)=-1.00000000000000E0_DP
 S_E%B_X(15,45)=1.00000000000000E0_DP
 S_E%A_Y(15,45)=1.00000000000000E0_DP
 S_E%A_X(15,54)=14.0000000000000E0_DP
 S_E%B_Y(15,54)=-14.0000000000000E0_DP
 S_E%B_X(15,64)=-91.0000000000000E0_DP
 S_E%A_Y(15,64)=-91.0000000000000E0_DP
 S_E%A_X(15,75)=-364.000000000000E0_DP
 S_E%B_Y(15,75)=364.000000000000E0_DP
 S_E%B_X(15,87)=1001.00000000000E0_DP
 S_E%A_Y(15,87)=1001.00000000000E0_DP
 S_E%A_X(15,100)=2002.00000000000E0_DP
 S_E%B_Y(15,100)=-2002.00000000000E0_DP
 S_E%B_X(15,114)=-3003.00000000000E0_DP
 S_E%A_Y(15,114)=-3003.00000000000E0_DP
 S_E%A_X(15,129)=-3432.00000000000E0_DP
 S_E%B_Y(15,129)=3432.00000000000E0_DP
 S_E%B_X(15,145)=3003.00000000000E0_DP
 S_E%A_Y(15,145)=3003.00000000000E0_DP
 S_E%A_X(15,162)=2002.00000000000E0_DP
 S_E%B_Y(15,162)=-2002.00000000000E0_DP
 S_E%B_X(15,180)=-1001.00000000000E0_DP
 S_E%A_Y(15,180)=-1001.00000000000E0_DP
 S_E%A_X(15,199)=-364.000000000000E0_DP
 S_E%B_Y(15,199)=364.000000000000E0_DP
 S_E%B_X(15,219)=91.0000000000000E0_DP
 S_E%A_Y(15,219)=91.0000000000000E0_DP
 S_E%A_X(15,240)=14.0000000000000E0_DP
 S_E%B_Y(15,240)=-14.0000000000000E0_DP
 S_E%B_X(15,262)=-1.00000000000000E0_DP
 S_E%A_Y(15,262)=-1.00000000000000E0_DP
 S_E%B_X(16,36)=1.00000000000000E0_DP
 S_E%A_Y(16,36)=1.00000000000000E0_DP
 S_E%A_X(16,44)=15.0000000000000E0_DP
 S_E%B_Y(16,44)=-15.0000000000000E0_DP
 S_E%B_X(16,53)=-105.000000000000E0_DP
 S_E%A_Y(16,53)=-105.000000000000E0_DP
 S_E%A_X(16,63)=-455.000000000000E0_DP
 S_E%B_Y(16,63)=455.000000000000E0_DP
 S_E%B_X(16,74)=1365.00000000000E0_DP
 S_E%A_Y(16,74)=1365.00000000000E0_DP
 S_E%A_X(16,86)=3003.00000000000E0_DP
 S_E%B_Y(16,86)=-3003.00000000000E0_DP
 S_E%B_X(16,99)=-5005.00000000000E0_DP
 S_E%A_Y(16,99)=-5005.00000000000E0_DP
 S_E%A_X(16,113)=-6435.00000000000E0_DP
 S_E%B_Y(16,113)=6435.00000000000E0_DP
 S_E%B_X(16,128)=6435.00000000000E0_DP
 S_E%A_Y(16,128)=6435.00000000000E0_DP
 S_E%A_X(16,144)=5005.00000000000E0_DP
 S_E%B_Y(16,144)=-5005.00000000000E0_DP
 S_E%B_X(16,161)=-3003.00000000000E0_DP
 S_E%A_Y(16,161)=-3003.00000000000E0_DP
 S_E%A_X(16,179)=-1365.00000000000E0_DP
 S_E%B_Y(16,179)=1365.00000000000E0_DP
 S_E%B_X(16,198)=455.000000000000E0_DP
 S_E%A_Y(16,198)=455.000000000000E0_DP
 S_E%A_X(16,218)=105.000000000000E0_DP
 S_E%B_Y(16,218)=-105.000000000000E0_DP
 S_E%B_X(16,239)=-15.0000000000000E0_DP
 S_E%A_Y(16,239)=-15.0000000000000E0_DP
 S_E%A_X(16,261)=-1.00000000000000E0_DP
 S_E%B_Y(16,261)=1.00000000000000E0_DP
 S_E%B_X(17,28)=1.00000000000000E0_DP
 S_E%A_Y(17,28)=1.00000000000000E0_DP
 S_E%A_X(17,35)=16.0000000000000E0_DP
 S_E%B_Y(17,35)=-16.0000000000000E0_DP
 S_E%B_X(17,43)=-120.000000000000E0_DP
 S_E%A_Y(17,43)=-120.000000000000E0_DP
 S_E%A_X(17,52)=-560.000000000000E0_DP
 S_E%B_Y(17,52)=560.000000000000E0_DP
 S_E%B_X(17,62)=1820.00000000000E0_DP
 S_E%A_Y(17,62)=1820.00000000000E0_DP
 S_E%A_X(17,73)=4368.00000000000E0_DP
 S_E%B_Y(17,73)=-4368.00000000000E0_DP
 S_E%B_X(17,85)=-8008.00000000000E0_DP
 S_E%A_Y(17,85)=-8008.00000000000E0_DP
 S_E%A_X(17,98)=-11440.0000000000E0_DP
 S_E%B_Y(17,98)=11440.0000000000E0_DP
 S_E%B_X(17,112)=12870.0000000000E0_DP
 S_E%A_Y(17,112)=12870.0000000000E0_DP
 S_E%A_X(17,127)=11440.0000000000E0_DP
 S_E%B_Y(17,127)=-11440.0000000000E0_DP
 S_E%B_X(17,143)=-8008.00000000000E0_DP
 S_E%A_Y(17,143)=-8008.00000000000E0_DP
 S_E%A_X(17,160)=-4368.00000000000E0_DP
 S_E%B_Y(17,160)=4368.00000000000E0_DP
 S_E%B_X(17,178)=1820.00000000000E0_DP
 S_E%A_Y(17,178)=1820.00000000000E0_DP
 S_E%A_X(17,197)=560.000000000000E0_DP
 S_E%B_Y(17,197)=-560.000000000000E0_DP
 S_E%B_X(17,217)=-120.000000000000E0_DP
 S_E%A_Y(17,217)=-120.000000000000E0_DP
 S_E%A_X(17,238)=-16.0000000000000E0_DP
 S_E%B_Y(17,238)=16.0000000000000E0_DP
 S_E%B_X(17,260)=1.00000000000000E0_DP
 S_E%A_Y(17,260)=1.00000000000000E0_DP
 S_E%B_X(18,21)=1.00000000000000E0_DP
 S_E%A_Y(18,21)=1.00000000000000E0_DP
 S_E%A_X(18,27)=17.0000000000000E0_DP
 S_E%B_Y(18,27)=-17.0000000000000E0_DP
 S_E%B_X(18,34)=-136.000000000000E0_DP
 S_E%A_Y(18,34)=-136.000000000000E0_DP
 S_E%A_X(18,42)=-680.000000000000E0_DP
 S_E%B_Y(18,42)=680.000000000000E0_DP
 S_E%B_X(18,51)=2380.00000000000E0_DP
 S_E%A_Y(18,51)=2380.00000000000E0_DP
 S_E%A_X(18,61)=6188.00000000000E0_DP
 S_E%B_Y(18,61)=-6188.00000000000E0_DP
 S_E%B_X(18,72)=-12376.0000000000E0_DP
 S_E%A_Y(18,72)=-12376.0000000000E0_DP
 S_E%A_X(18,84)=-19448.0000000000E0_DP
 S_E%B_Y(18,84)=19448.0000000000E0_DP
 S_E%B_X(18,97)=24310.0000000000E0_DP
 S_E%A_Y(18,97)=24310.0000000000E0_DP
 S_E%A_X(18,111)=24310.0000000000E0_DP
 S_E%B_Y(18,111)=-24310.0000000000E0_DP
 S_E%B_X(18,126)=-19448.0000000000E0_DP
 S_E%A_Y(18,126)=-19448.0000000000E0_DP
 S_E%A_X(18,142)=-12376.0000000000E0_DP
 S_E%B_Y(18,142)=12376.0000000000E0_DP
 S_E%B_X(18,159)=6188.00000000000E0_DP
 S_E%A_Y(18,159)=6188.00000000000E0_DP
 S_E%A_X(18,177)=2380.00000000000E0_DP
 S_E%B_Y(18,177)=-2380.00000000000E0_DP
 S_E%B_X(18,196)=-680.000000000000E0_DP
 S_E%A_Y(18,196)=-680.000000000000E0_DP
 S_E%A_X(18,216)=-136.000000000000E0_DP
 S_E%B_Y(18,216)=136.000000000000E0_DP
 S_E%B_X(18,237)=17.0000000000000E0_DP
 S_E%A_Y(18,237)=17.0000000000000E0_DP
 S_E%A_X(18,259)=1.00000000000000E0_DP
 S_E%B_Y(18,259)=-1.00000000000000E0_DP
 S_E%B_X(19,15)=1.00000000000000E0_DP
 S_E%A_Y(19,15)=1.00000000000000E0_DP
 S_E%A_X(19,20)=18.0000000000000E0_DP
 S_E%B_Y(19,20)=-18.0000000000000E0_DP
 S_E%B_X(19,26)=-153.000000000000E0_DP
 S_E%A_Y(19,26)=-153.000000000000E0_DP
 S_E%A_X(19,33)=-816.000000000000E0_DP
 S_E%B_Y(19,33)=816.000000000000E0_DP
 S_E%B_X(19,41)=3060.00000000000E0_DP
 S_E%A_Y(19,41)=3060.00000000000E0_DP
 S_E%A_X(19,50)=8568.00000000000E0_DP
 S_E%B_Y(19,50)=-8568.00000000000E0_DP
 S_E%B_X(19,60)=-18564.0000000000E0_DP
 S_E%A_Y(19,60)=-18564.0000000000E0_DP
 S_E%A_X(19,71)=-31824.0000000000E0_DP
 S_E%B_Y(19,71)=31824.0000000000E0_DP
 S_E%B_X(19,83)=43758.0000000000E0_DP
 S_E%A_Y(19,83)=43758.0000000000E0_DP
 S_E%A_X(19,96)=48620.0000000000E0_DP
 S_E%B_Y(19,96)=-48620.0000000000E0_DP
 S_E%B_X(19,110)=-43758.0000000000E0_DP
 S_E%A_Y(19,110)=-43758.0000000000E0_DP
 S_E%A_X(19,125)=-31824.0000000000E0_DP
 S_E%B_Y(19,125)=31824.0000000000E0_DP
 S_E%B_X(19,141)=18564.0000000000E0_DP
 S_E%A_Y(19,141)=18564.0000000000E0_DP
 S_E%A_X(19,158)=8568.00000000000E0_DP
 S_E%B_Y(19,158)=-8568.00000000000E0_DP
 S_E%B_X(19,176)=-3060.00000000000E0_DP
 S_E%A_Y(19,176)=-3060.00000000000E0_DP
 S_E%A_X(19,195)=-816.000000000000E0_DP
 S_E%B_Y(19,195)=816.000000000000E0_DP
 S_E%B_X(19,215)=153.000000000000E0_DP
 S_E%A_Y(19,215)=153.000000000000E0_DP
 S_E%A_X(19,236)=18.0000000000000E0_DP
 S_E%B_Y(19,236)=-18.0000000000000E0_DP
 S_E%B_X(19,258)=-1.00000000000000E0_DP
 S_E%A_Y(19,258)=-1.00000000000000E0_DP
 S_E%B_X(20,10)=1.00000000000000E0_DP
 S_E%A_Y(20,10)=1.00000000000000E0_DP
 S_E%A_X(20,14)=19.0000000000000E0_DP
 S_E%B_Y(20,14)=-19.0000000000000E0_DP
 S_E%B_X(20,19)=-171.000000000000E0_DP
 S_E%A_Y(20,19)=-171.000000000000E0_DP
 S_E%A_X(20,25)=-969.000000000000E0_DP
 S_E%B_Y(20,25)=969.000000000000E0_DP
 S_E%B_X(20,32)=3876.00000000000E0_DP
 S_E%A_Y(20,32)=3876.00000000000E0_DP
 S_E%A_X(20,40)=11628.0000000000E0_DP
 S_E%B_Y(20,40)=-11628.0000000000E0_DP
 S_E%B_X(20,49)=-27132.0000000000E0_DP
 S_E%A_Y(20,49)=-27132.0000000000E0_DP
 S_E%A_X(20,59)=-50388.0000000000E0_DP
 S_E%B_Y(20,59)=50388.0000000000E0_DP
 S_E%B_X(20,70)=75582.0000000000E0_DP
 S_E%A_Y(20,70)=75582.0000000000E0_DP
 S_E%A_X(20,82)=92378.0000000000E0_DP
 S_E%B_Y(20,82)=-92378.0000000000E0_DP
 S_E%B_X(20,95)=-92378.0000000000E0_DP
 S_E%A_Y(20,95)=-92378.0000000000E0_DP
 S_E%A_X(20,109)=-75582.0000000000E0_DP
 S_E%B_Y(20,109)=75582.0000000000E0_DP
 S_E%B_X(20,124)=50388.0000000000E0_DP
 S_E%A_Y(20,124)=50388.0000000000E0_DP
 S_E%A_X(20,140)=27132.0000000000E0_DP
 S_E%B_Y(20,140)=-27132.0000000000E0_DP
 S_E%B_X(20,157)=-11628.0000000000E0_DP
 S_E%A_Y(20,157)=-11628.0000000000E0_DP
 S_E%A_X(20,175)=-3876.00000000000E0_DP
 S_E%B_Y(20,175)=3876.00000000000E0_DP
 S_E%B_X(20,194)=969.000000000000E0_DP
 S_E%A_Y(20,194)=969.000000000000E0_DP
 S_E%A_X(20,214)=171.000000000000E0_DP
 S_E%B_Y(20,214)=-171.000000000000E0_DP
 S_E%B_X(20,235)=-19.0000000000000E0_DP
 S_E%A_Y(20,235)=-19.0000000000000E0_DP
 S_E%A_X(20,257)=-1.00000000000000E0_DP
 S_E%B_Y(20,257)=1.00000000000000E0_DP
 S_E%B_X(21,6)=1.00000000000000E0_DP
 S_E%A_Y(21,6)=1.00000000000000E0_DP
 S_E%A_X(21,9)=20.0000000000000E0_DP
 S_E%B_Y(21,9)=-20.0000000000000E0_DP
 S_E%B_X(21,13)=-190.000000000000E0_DP
 S_E%A_Y(21,13)=-190.000000000000E0_DP
 S_E%A_X(21,18)=-1140.00000000000E0_DP
 S_E%B_Y(21,18)=1140.00000000000E0_DP
 S_E%B_X(21,24)=4845.00000000000E0_DP
 S_E%A_Y(21,24)=4845.00000000000E0_DP
 S_E%A_X(21,31)=15504.0000000000E0_DP
 S_E%B_Y(21,31)=-15504.0000000000E0_DP
 S_E%B_X(21,39)=-38760.0000000000E0_DP
 S_E%A_Y(21,39)=-38760.0000000000E0_DP
 S_E%A_X(21,48)=-77520.0000000000E0_DP
 S_E%B_Y(21,48)=77520.0000000000E0_DP
 S_E%B_X(21,58)=125970.000000000E0_DP
 S_E%A_Y(21,58)=125970.000000000E0_DP
 S_E%A_X(21,69)=167960.000000000E0_DP
 S_E%B_Y(21,69)=-167960.000000000E0_DP
 S_E%B_X(21,81)=-184756.000000000E0_DP
 S_E%A_Y(21,81)=-184756.000000000E0_DP
 S_E%A_X(21,94)=-167960.000000000E0_DP
 S_E%B_Y(21,94)=167960.000000000E0_DP
 S_E%B_X(21,108)=125970.000000000E0_DP
 S_E%A_Y(21,108)=125970.000000000E0_DP
 S_E%A_X(21,123)=77520.0000000000E0_DP
 S_E%B_Y(21,123)=-77520.0000000000E0_DP
 S_E%B_X(21,139)=-38760.0000000000E0_DP
 S_E%A_Y(21,139)=-38760.0000000000E0_DP
 S_E%A_X(21,156)=-15504.0000000000E0_DP
 S_E%B_Y(21,156)=15504.0000000000E0_DP
 S_E%B_X(21,174)=4845.00000000000E0_DP
 S_E%A_Y(21,174)=4845.00000000000E0_DP
 S_E%A_X(21,193)=1140.00000000000E0_DP
 S_E%B_Y(21,193)=-1140.00000000000E0_DP
 S_E%B_X(21,213)=-190.000000000000E0_DP
 S_E%A_Y(21,213)=-190.000000000000E0_DP
 S_E%A_X(21,234)=-20.0000000000000E0_DP
 S_E%B_Y(21,234)=20.0000000000000E0_DP
 S_E%B_X(21,256)=1.00000000000000E0_DP
 S_E%A_Y(21,256)=1.00000000000000E0_DP
 S_E%B_X(22,3)=1.00000000000000E0_DP
 S_E%A_Y(22,3)=1.00000000000000E0_DP
 S_E%A_X(22,5)=21.0000000000000E0_DP
 S_E%B_Y(22,5)=-21.0000000000000E0_DP
 S_E%B_X(22,8)=-210.000000000000E0_DP
 S_E%A_Y(22,8)=-210.000000000000E0_DP
 S_E%A_X(22,12)=-1330.00000000000E0_DP
 S_E%B_Y(22,12)=1330.00000000000E0_DP
 S_E%B_X(22,17)=5985.00000000000E0_DP
 S_E%A_Y(22,17)=5985.00000000000E0_DP
 S_E%A_X(22,23)=20349.0000000000E0_DP
 S_E%B_Y(22,23)=-20349.0000000000E0_DP
 S_E%B_X(22,30)=-54264.0000000000E0_DP
 S_E%A_Y(22,30)=-54264.0000000000E0_DP
 S_E%A_X(22,38)=-116280.000000000E0_DP
 S_E%B_Y(22,38)=116280.000000000E0_DP
 S_E%B_X(22,47)=203490.000000000E0_DP
 S_E%A_Y(22,47)=203490.000000000E0_DP
 S_E%A_X(22,57)=293930.000000000E0_DP
 S_E%B_Y(22,57)=-293930.000000000E0_DP
 S_E%B_X(22,68)=-352716.000000000E0_DP
 S_E%A_Y(22,68)=-352716.000000000E0_DP
 S_E%A_X(22,80)=-352716.000000000E0_DP
 S_E%B_Y(22,80)=352716.000000000E0_DP
 S_E%B_X(22,93)=293930.000000000E0_DP
 S_E%A_Y(22,93)=293930.000000000E0_DP
 S_E%A_X(22,107)=203490.000000000E0_DP
 S_E%B_Y(22,107)=-203490.000000000E0_DP
 S_E%B_X(22,122)=-116280.000000000E0_DP
 S_E%A_Y(22,122)=-116280.000000000E0_DP
 S_E%A_X(22,138)=-54264.0000000000E0_DP
 S_E%B_Y(22,138)=54264.0000000000E0_DP
 S_E%B_X(22,155)=20349.0000000000E0_DP
 S_E%A_Y(22,155)=20349.0000000000E0_DP
 S_E%A_X(22,173)=5985.00000000000E0_DP
 S_E%B_Y(22,173)=-5985.00000000000E0_DP
 S_E%B_X(22,192)=-1330.00000000000E0_DP
 S_E%A_Y(22,192)=-1330.00000000000E0_DP
 S_E%A_X(22,212)=-210.000000000000E0_DP
 S_E%B_Y(22,212)=210.000000000000E0_DP
 S_E%B_X(22,233)=21.0000000000000E0_DP
 S_E%A_Y(22,233)=21.0000000000000E0_DP
 S_E%A_X(22,255)=1.00000000000000E0_DP
 S_E%B_Y(22,255)=-1.00000000000000E0_DP
 S_E%VB(1,90)=-3.632077277826440E-018_DP
 S_E%VB(1,91)=-8.673617379884035E-019_DP
 S_E%VB(1,103)=-0.500000000000000E0_DP
 S_E%VB(1,104)=5.963111948670274E-019_DP
 S_E%VB(1,117)=1.517883041479706E-018_DP
 S_E%VB(1,118)=0.500000000000000E0_DP
 S_E%VB(1,119)=2.818925648462312E-018_DP
 S_E%VB(1,132)=1.50000000000000E0_DP
 S_E%VB(1,133)=6.505213034913027E-018_DP
 S_E%VB(1,134)=-0.500000000000000E0_DP
 S_E%VB(1,135)=-3.252606517456513E-018_DP
 S_E%VB(1,148)=-1.084202172485504E-019_DP
 S_E%VB(1,149)=-1.16666666666667E0_DP
 S_E%VB(1,150)=4.770489558936220E-018_DP
 S_E%VB(1,151)=0.500000000000000E0_DP
 S_E%VB(1,152)=-8.673617379884035E-019_DP
 S_E%VB(1,165)=-1.57500000000000E0_DP
 S_E%VB(1,166)=1.951563910473908E-018_DP
 S_E%VB(1,167)=0.875000000000000E0_DP
 S_E%VB(1,168)=-1.734723475976807E-018_DP
 S_E%VB(1,169)=-0.500000000000000E0_DP
 S_E%VB(1,170)=2.602085213965211E-018_DP
 S_E%VB(1,171)=6.938893903907228E-018_DP
 S_E%VB(1,183)=-9.324138683375338E-018_DP
 S_E%VB(1,184)=0.874999999999998E0_DP
 S_E%VB(1,185)=-9.324138683375338E-018_DP
 S_E%VB(1,186)=-0.625000000000000E0_DP
 S_E%VB(1,187)=4.336808689942018E-018_DP
 S_E%VB(1,188)=0.500000000000000E0_DP
 S_E%VB(1,202)=0.468750000000000E0_DP
 S_E%VB(1,203)=-4.770489558936220E-018_DP
 S_E%VB(1,204)=-0.437500000000000E0_DP
 S_E%VB(1,205)=-3.469446951953614E-018_DP
 S_E%VB(1,206)=0.416666666666667E0_DP
 S_E%VB(1,208)=-0.500000000000000E0_DP
 S_E%VB(1,222)=-1.301042606982605E-018_DP
 S_E%VB(1,223)=-0.156250000000000E0_DP
 S_E%VB(1,224)=-1.734723475976807E-018_DP
 S_E%VB(1,225)=0.187500000000000E0_DP
 S_E%VB(1,226)=-1.734723475976807E-018_DP
 S_E%VB(1,227)=-0.250000000000000E0_DP
 S_E%VB(1,228)=-3.469446951953614E-018_DP
 S_E%VB(1,229)=0.500000000000000E0_DP
 S_E%VB(1,243)=-2.734375000000002E-002_DP
 S_E%VB(1,244)=-2.764715539838036E-018_DP
 S_E%VB(1,245)=3.906250000000002E-002_DP
 S_E%VB(1,246)=-1.951563910473908E-018_DP
 S_E%VB(1,247)=-6.249999999999997E-002_DP
 S_E%VB(1,249)=0.125000000000000E0_DP
 S_E%VB(1,251)=-0.500000000000000E0_DP
 S_E%VB(1,253)=-1.00000000000000E0_DP
 S_E%VB(1,265)=-5.082197683525802E-020_DP
 S_E%VB(1,266)=3.038194444444437E-003_DP
 S_E%VB(1,267)=-2.710505431213761E-020_DP
 S_E%VB(1,268)=-5.580357142857142E-003_DP
 S_E%VB(1,269)=-5.421010862427522E-019_DP
 S_E%VB(1,270)=1.250000000000001E-002_DP
 S_E%VB(1,271)=8.673617379884035E-019_DP
 S_E%VB(1,272)=-4.166666666666664E-002_DP
 S_E%VB(1,274)=0.500000000000000E0_DP
 S_E%VA(1,275)=-1.00000000000000E0_DP
 S_E%VA(2,78)=-1.863472483959461E-020_DP
 S_E%VB(2,90)=1.382357769919018E-018_DP
 S_E%VA(2,91)=-1.118083490375676E-019_DP
 S_E%VA(2,103)=1.355252715606881E-018_DP
 S_E%VB(2,103)=0.500000000000000E0_DP
 S_E%VB(2,104)=-2.439454888092385E-018_DP
 S_E%VA(2,105)=1.490777987167569E-019_DP
 S_E%VB(2,105)=1.734723475976807E-018_DP
 S_E%VA(2,117)=0.166666666666667E0_DP
 S_E%VB(2,117)=8.673617379884035E-019_DP
 S_E%VA(2,118)=-9.757819552369540E-019_DP
 S_E%VB(2,118)=-0.500000000000000E0_DP
 S_E%VB(2,119)=8.944667923005412E-019_DP
 S_E%VA(2,120)=-3.252606517456513E-019_DP
 S_E%VA(2,132)=3.252606517456513E-019_DP
 S_E%VB(2,132)=-1.50000000000000E0_DP
 S_E%VA(2,133)=-0.166666666666667E0_DP
 S_E%VB(2,133)=6.396792817664476E-018_DP
 S_E%VA(2,134)=-2.168404344971009E-019_DP
 S_E%VB(2,134)=0.500000000000000E0_DP
 S_E%VB(2,135)=-8.131516293641283E-019_DP
 S_E%VA(2,136)=1.626303258728257E-019_DP
 S_E%VA(2,148)=-0.233333333333333E0_DP
 S_E%VB(2,148)=-2.764715539838036E-018_DP
 S_E%VA(2,149)=2.927345865710862E-018_DP
 S_E%VB(2,149)=1.16666666666667E0_DP
 S_E%VA(2,150)=0.166666666666667E0_DP
 S_E%VB(2,150)=-1.843143693225358E-018_DP
 S_E%VA(2,151)=2.602085213965211E-018_DP
 S_E%VB(2,151)=-0.500000000000000E0_DP
 S_E%VB(2,152)=6.505213034913027E-019_DP
 S_E%VA(2,153)=-1.084202172485504E-019_DP
 S_E%VA(2,165)=1.680513367352532E-018_DP
 S_E%VB(2,165)=1.57500000000000E0_DP
 S_E%VA(2,166)=0.175000000000000E0_DP
 S_E%VB(2,167)=-0.875000000000000E0_DP
 S_E%VA(2,168)=-0.166666666666667E0_DP
 S_E%VB(2,168)=-1.138412281109780E-018_DP
 S_E%VB(2,169)=0.500000000000000E0_DP
 S_E%VA(2,171)=8.673617379884035E-019_DP
 S_E%VA(2,183)=0.125000000000000E0_DP
 S_E%VB(2,183)=-4.119968255444917E-018_DP
 S_E%VA(2,184)=2.710505431213761E-018_DP
 S_E%VB(2,184)=-0.875000000000000E0_DP
 S_E%VA(2,185)=-0.125000000000000E0_DP
 S_E%VB(2,185)=2.547875105340935E-018_DP
 S_E%VA(2,186)=2.710505431213761E-020_DP
 S_E%VB(2,186)=0.625000000000000E0_DP
 S_E%VA(2,187)=0.166666666666667E0_DP
 S_E%VB(2,187)=-6.505213034913027E-019_DP
 S_E%VB(2,188)=-0.500000000000000E0_DP
 S_E%VA(2,202)=-5.421010862427522E-019_DP
 S_E%VB(2,202)=-0.468749999999999E0_DP
 S_E%VA(2,203)=-6.249999999999996E-002_DP
 S_E%VB(2,203)=8.348356728138384E-018_DP
 S_E%VA(2,204)=2.168404344971009E-019_DP
 S_E%VB(2,204)=0.437500000000000E0_DP
 S_E%VA(2,205)=8.333333333333333E-002_DP
 S_E%VB(2,205)=-1.138412281109780E-018_DP
 S_E%VA(2,206)=-1.734723475976807E-018_DP
 S_E%VB(2,206)=-0.416666666666667E0_DP
 S_E%VA(2,207)=-0.166666666666667E0_DP
 S_E%VB(2,208)=0.500000000000000E0_DP
 S_E%VA(2,222)=-1.736111111111113E-002_DP
 S_E%VB(2,222)=-3.144186300207963E-018_DP
 S_E%VB(2,223)=0.156250000000000E0_DP
 S_E%VA(2,224)=2.678571428571427E-002_DP
 S_E%VB(2,224)=-5.421010862427522E-020_DP
 S_E%VA(2,225)=-4.336808689942018E-019_DP
 S_E%VB(2,225)=-0.187500000000000E0_DP
 S_E%VA(2,226)=-4.999999999999998E-002_DP
 S_E%VB(2,226)=8.673617379884035E-019_DP
 S_E%VB(2,227)=0.250000000000000E0_DP
 S_E%VA(2,228)=0.166666666666667E0_DP
 S_E%VB(2,229)=-0.500000000000000E0_DP
 S_E%VB(2,231)=-0.500000000000000E0_DP
 S_E%VA(2,243)=2.303929616531697E-019_DP
 S_E%VB(2,243)=2.734374999999996E-002_DP
 S_E%VA(2,244)=4.340277777777775E-003_DP
 S_E%VB(2,244)=4.878909776184770E-019_DP
 S_E%VA(2,245)=-2.710505431213761E-019_DP
 S_E%VB(2,245)=-3.906250000000003E-002_DP
 S_E%VA(2,246)=-8.928571428571428E-003_DP
 S_E%VB(2,246)=9.215718466126788E-019_DP
 S_E%VA(2,247)=2.168404344971009E-019_DP
 S_E%VB(2,247)=6.249999999999998E-002_DP
 S_E%VA(2,248)=2.500000000000000E-002_DP
 S_E%VB(2,248)=-8.673617379884035E-019_DP
 S_E%VB(2,249)=-0.125000000000000E0_DP
 S_E%VA(2,250)=-0.166666666666667E0_DP
 S_E%VB(2,251)=0.500000000000000E0_DP
 S_E%VA(2,252)=-1.00000000000000E0_DP
 S_E%VA(2,265)=2.761994949494954E-004_DP
 S_E%VB(2,265)=-3.913292216314868E-019_DP
 S_E%VA(2,266)=-1.660184576618429E-019_DP
 S_E%VB(2,266)=-3.038194444444442E-003_DP
 S_E%VA(2,267)=-6.200396825396814E-004_DP
 S_E%VB(2,267)=-7.284483346386983E-020_DP
 S_E%VA(2,268)=5.963111948670274E-019_DP
 S_E%VB(2,268)=5.580357142857139E-003_DP
 S_E%VA(2,269)=1.785714285714285E-003_DP
 S_E%VB(2,269)=8.131516293641283E-020_DP
 S_E%VA(2,270)=-1.084202172485504E-019_DP
 S_E%VB(2,270)=-1.250000000000000E-002_DP
 S_E%VA(2,271)=-8.333333333333335E-003_DP
 S_E%VB(2,271)=-8.673617379884035E-019_DP
 S_E%VB(2,272)=4.166666666666666E-002_DP
 S_E%VA(2,273)=0.166666666666667E0_DP
 S_E%VB(2,274)=0.500000000000000E0_DP
 S_E%VA(3,78)=2.642742795433417E-019_DP
 S_E%VB(3,90)=3.827318372168556E-018_DP
 S_E%VA(3,91)=5.827586677109586E-019_DP
 S_E%VA(3,103)=-4.336808689942018E-019_DP
 S_E%VB(3,103)=-0.500000000000000E0_DP
 S_E%VB(3,104)=2.331034670843835E-018_DP
 S_E%VA(3,105)=2.710505431213761E-019_DP
 S_E%VA(3,117)=-0.333333333333333E0_DP
 S_E%VB(3,117)=2.727446090158847E-018_DP
 S_E%VA(3,118)=-1.734723475976807E-018_DP
 S_E%VB(3,118)=0.500000000000000E0_DP
 S_E%VB(3,119)=-9.757819552369540E-019_DP
 S_E%VA(3,120)=3.794707603699266E-019_DP
 S_E%VA(3,132)=-8.673617379884035E-019_DP
 S_E%VB(3,132)=1.62500000000000E0_DP
 S_E%VA(3,133)=0.333333333333333E0_DP
 S_E%VB(3,133)=6.830473686658678E-018_DP
 S_E%VB(3,134)=-0.500000000000000E0_DP
 S_E%VB(3,135)=8.239936510889834E-018_DP
 S_E%VA(3,136)=5.963111948670274E-019_DP
 S_E%VB(3,136)=-6.938893903907228E-018_DP
 S_E%VA(3,148)=0.466666666666667E0_DP
 S_E%VB(3,148)=-8.016319812814698E-018_DP
 S_E%VA(3,149)=-5.204170427930421E-018_DP
 S_E%VB(3,149)=-1.29166666666667E0_DP
 S_E%VA(3,150)=-0.333333333333333E0_DP
 S_E%VB(3,150)=-1.301042606982605E-018_DP
 S_E%VB(3,151)=0.500000000000000E0_DP
 S_E%VA(3,165)=4.336808689942018E-019_DP
 S_E%VB(3,165)=-1.66250000000000E0_DP
 S_E%VA(3,166)=-0.350000000000000E0_DP
 S_E%VA(3,167)=-1.084202172485504E-018_DP
 S_E%VB(3,167)=1.00000000000000E0_DP
 S_E%VA(3,168)=0.333333333333333E0_DP
 S_E%VB(3,168)=4.336808689942018E-019_DP
 S_E%VB(3,169)=-0.500000000000000E0_DP
 S_E%VB(3,170)=-3.469446951953614E-018_DP
 S_E%VA(3,171)=-1.734723475976807E-018_DP
 S_E%VB(3,171)=6.938893903907228E-018_DP
 S_E%VA(3,183)=-0.250000000000000E0_DP
 S_E%VB(3,183)=-2.927345865710862E-018_DP
 S_E%VA(3,184)=-6.938893903907228E-018_DP
 S_E%VB(3,184)=0.937499999999999E0_DP
 S_E%VA(3,185)=0.250000000000000E0_DP
 S_E%VB(3,185)=-1.517883041479706E-018_DP
 S_E%VB(3,186)=-0.750000000000000E0_DP
 S_E%VA(3,187)=-0.333333333333333E0_DP
 S_E%VB(3,187)=8.673617379884035E-019_DP
 S_E%VB(3,188)=0.500000000000000E0_DP
 S_E%VA(3,202)=-1.301042606982605E-018_DP
 S_E%VB(3,202)=0.492187500000000E0_DP
 S_E%VA(3,203)=0.125000000000000E0_DP
 S_E%VB(3,203)=1.084202172485504E-019_DP
 S_E%VA(3,204)=-4.336808689942018E-019_DP
 S_E%VB(3,204)=-0.479166666666667E0_DP
 S_E%VA(3,205)=-0.166666666666667E0_DP
 S_E%VB(3,205)=-4.770489558936220E-018_DP
 S_E%VA(3,206)=4.987329993433320E-018_DP
 S_E%VB(3,206)=0.541666666666667E0_DP
 S_E%VA(3,207)=0.333333333333333E0_DP
 S_E%VA(3,208)=8.673617379884035E-019_DP
 S_E%VB(3,208)=-0.500000000000000E0_DP
 S_E%VB(3,210)=-0.333333333333333E0_DP
 S_E%VA(3,222)=3.472222222222220E-002_DP
 S_E%VB(3,222)=-4.031876828930470E-019_DP
 S_E%VA(3,223)=-8.673617379884035E-019_DP
 S_E%VB(3,223)=-0.166294642857143E0_DP
 S_E%VA(3,224)=-5.357142857142856E-002_DP
 S_E%VB(3,224)=-2.168404344971009E-019_DP
 S_E%VA(3,225)=-8.673617379884035E-019_DP
 S_E%VB(3,225)=0.212500000000000E0_DP
 S_E%VA(3,226)=0.100000000000000E0_DP
 S_E%VB(3,226)=-4.336808689942018E-018_DP
 S_E%VA(3,227)=4.336808689942018E-019_DP
 S_E%VB(3,227)=-0.375000000000000E0_DP
 S_E%VA(3,228)=-0.333333333333333E0_DP
 S_E%VB(3,228)=-3.469446951953614E-018_DP
 S_E%VB(3,229)=0.500000000000000E0_DP
 S_E%VA(3,230)=-1.00000000000000E0_DP
 S_E%VA(3,243)=-1.626303258728257E-019_DP
 S_E%VB(3,243)=-2.864583333333333E-002_DP
 S_E%VA(3,244)=-8.680555555555556E-003_DP
 S_E%VB(3,244)=-3.794707603699266E-019_DP
 S_E%VA(3,245)=-3.252606517456513E-019_DP
 S_E%VB(3,245)=4.241071428571433E-002_DP
 S_E%VA(3,246)=1.785714285714286E-002_DP
 S_E%VB(3,246)=1.301042606982605E-018_DP
 S_E%VA(3,247)=2.168404344971009E-019_DP
 S_E%VB(3,247)=-7.499999999999998E-002_DP
 S_E%VA(3,248)=-5.000000000000000E-002_DP
 S_E%VB(3,248)=-5.204170427930421E-018_DP
 S_E%VA(3,249)=-8.673617379884035E-019_DP
 S_E%VB(3,249)=0.250000000000000E0_DP
 S_E%VA(3,250)=0.333333333333333E0_DP
 S_E%VB(3,251)=1.00000000000000E0_DP
 S_E%VA(3,265)=-5.523989898989897E-004_DP
 S_E%VB(3,265)=-3.412219288454746E-019_DP
 S_E%VA(3,266)=4.201283418381330E-019_DP
 S_E%VB(3,266)=3.224206349206343E-003_DP
 S_E%VA(3,267)=1.240079365079365E-003_DP
 S_E%VB(3,267)=-9.486769009248164E-020_DP
 S_E%VA(3,268)=-9.215718466126788E-019_DP
 S_E%VB(3,268)=-6.250000000000000E-003_DP
 S_E%VA(3,269)=-3.571428571428571E-003_DP
 S_E%VB(3,269)=4.336808689942018E-019_DP
 S_E%VB(3,270)=1.666666666666667E-002_DP
 S_E%VA(3,271)=1.666666666666667E-002_DP
 S_E%VB(3,271)=-8.673617379884035E-019_DP
 S_E%VB(3,272)=-0.166666666666667E0_DP
 S_E%VA(3,273)=0.333333333333333E0_DP
 S_E%VB(4,90)=2.005774019098183E-018_DP
 S_E%VA(4,91)=-2.574980159653073E-019_DP
 S_E%VA(4,103)=-1.517883041479706E-018_DP
 S_E%VB(4,103)=0.500000000000000E0_DP
 S_E%VB(4,104)=-2.493664996716660E-018_DP
 S_E%VA(4,105)=-3.320369153236857E-019_DP
 S_E%VB(4,105)=3.469446951953614E-018_DP
 S_E%VA(4,117)=0.500000000000000E0_DP
 S_E%VB(4,117)=-1.517883041479706E-018_DP
 S_E%VA(4,118)=4.770489558936220E-018_DP
 S_E%VB(4,118)=-0.500000000000000E0_DP
 S_E%VB(4,119)=1.734723475976807E-018_DP
 S_E%VA(4,120)=2.981555974335137E-019_DP
 S_E%VA(4,132)=-2.602085213965211E-018_DP
 S_E%VB(4,132)=-1.87500000000000E0_DP
 S_E%VA(4,133)=-0.500000000000000E0_DP
 S_E%VB(4,133)=-3.903127820947816E-018_DP
 S_E%VA(4,134)=-1.599198204416119E-018_DP
 S_E%VB(4,134)=0.500000000000000E0_DP
 S_E%VB(4,135)=1.734723475976807E-018_DP
 S_E%VA(4,136)=-9.757819552369540E-019_DP
 S_E%VA(4,148)=-0.775000000000000E0_DP
 S_E%VB(4,148)=-6.505213034913027E-018_DP
 S_E%VA(4,149)=3.903127820947816E-018_DP
 S_E%VB(4,149)=1.54166666666667E0_DP
 S_E%VA(4,150)=0.500000000000000E0_DP
 S_E%VA(4,151)=4.336808689942018E-019_DP
 S_E%VB(4,151)=-0.500000000000000E0_DP
 S_E%VB(4,153)=-6.938893903907228E-018_DP
 S_E%VA(4,165)=-2.818925648462312E-018_DP
 S_E%VB(4,165)=1.83750000000000E0_DP
 S_E%VA(4,166)=0.600000000000000E0_DP
 S_E%VA(4,167)=-9.486769009248164E-020_DP
 S_E%VB(4,167)=-1.25000000000000E0_DP
 S_E%VA(4,168)=-0.500000000000000E0_DP
 S_E%VA(4,169)=8.673617379884035E-019_DP
 S_E%VB(4,169)=0.500000000000000E0_DP
 S_E%VA(4,183)=0.401785714285715E0_DP
 S_E%VB(4,183)=4.336808689942018E-019_DP
 S_E%VA(4,184)=2.168404344971009E-018_DP
 S_E%VB(4,184)=-1.06250000000000E0_DP
 S_E%VA(4,185)=-0.450000000000000E0_DP
 S_E%VB(4,185)=8.673617379884035E-019_DP
 S_E%VB(4,186)=1.00000000000000E0_DP
 S_E%VA(4,187)=0.500000000000000E0_DP
 S_E%VB(4,187)=1.734723475976807E-018_DP
 S_E%VB(4,188)=-0.500000000000000E0_DP
 S_E%VB(4,190)=-0.250000000000000E0_DP
 S_E%VA(4,202)=3.035766082959412E-018_DP
 S_E%VB(4,202)=-0.539062500000000E0_DP
 S_E%VA(4,203)=-0.205357142857143E0_DP
 S_E%VB(4,203)=-2.168404344971009E-018_DP
 S_E%VA(4,204)=1.544988095791844E-018_DP
 S_E%VB(4,204)=0.562500000000000E0_DP
 S_E%VA(4,205)=0.325000000000000E0_DP
 S_E%VB(4,205)=3.469446951953614E-018_DP
 S_E%VA(4,206)=-4.336808689942018E-019_DP
 S_E%VB(4,206)=-0.791666666666667E0_DP
 S_E%VA(4,207)=-0.500000000000000E0_DP
 S_E%VB(4,208)=0.500000000000000E0_DP
 S_E%VA(4,209)=-1.00000000000000E0_DP
 S_E%VA(4,222)=-5.543154761904764E-002_DP
 S_E%VB(4,222)=6.505213034913027E-019_DP
 S_E%VA(4,223)=3.903127820947816E-018_DP
 S_E%VB(4,223)=0.186383928571429E0_DP
 S_E%VA(4,224)=9.107142857142854E-002_DP
 S_E%VB(4,224)=-3.469446951953614E-018_DP
 S_E%VA(4,225)=-4.336808689942018E-019_DP
 S_E%VB(4,225)=-0.262500000000000E0_DP
 S_E%VA(4,226)=-0.225000000000000E0_DP
 S_E%VB(4,227)=0.625000000000000E0_DP
 S_E%VA(4,228)=0.500000000000000E0_DP
 S_E%VB(4,229)=1.50000000000000E0_DP
 S_E%VB(4,243)=3.124999999999997E-002_DP
 S_E%VA(4,244)=1.413690476190477E-002_DP
 S_E%VB(4,244)=-1.084202172485504E-019_DP
 S_E%VA(4,245)=4.743384504624082E-020_DP
 S_E%VB(4,245)=-4.910714285714287E-002_DP
 S_E%VA(4,246)=-3.214285714285715E-002_DP
 S_E%VB(4,246)=-8.673617379884035E-019_DP
 S_E%VA(4,247)=4.336808689942018E-019_DP
 S_E%VB(4,247)=0.100000000000000E0_DP
 S_E%VA(4,248)=0.150000000000000E0_DP
 S_E%VB(4,249)=-0.500000000000000E0_DP
 S_E%VA(4,250)=1.00000000000000E0_DP
 S_E%VA(4,265)=8.793290043290049E-004_DP
 S_E%VB(4,265)=1.897353801849633E-019_DP
 S_E%VA(4,266)=-2.574980159653073E-019_DP
 S_E%VB(4,266)=-3.596230158730165E-003_DP
 S_E%VA(4,267)=-2.083333333333332E-003_DP
 S_E%VB(4,267)=5.421010862427522E-020_DP
 S_E%VA(4,268)=5.692061405548898E-019_DP
 S_E%VB(4,268)=7.589285714285718E-003_DP
 S_E%VA(4,269)=7.142857142857141E-003_DP
 S_E%VB(4,269)=-4.336808689942018E-019_DP
 S_E%VB(4,270)=-2.500000000000000E-002_DP
 S_E%VA(4,271)=-0.100000000000000E0_DP
 S_E%VB(4,272)=-0.250000000000000E0_DP
 S_E%VA(5,78)=-3.049318610115481E-020_DP
 S_E%VB(5,90)=2.168404344971009E-019_DP
 S_E%VA(5,91)=3.252606517456513E-019_DP
 S_E%VA(5,103)=-4.119968255444917E-018_DP
 S_E%VB(5,103)=-0.500000000000000E0_DP
 S_E%VB(5,104)=-4.336808689942018E-019_DP
 S_E%VA(5,105)=5.963111948670274E-019_DP
 S_E%VA(5,117)=-0.666666666666667E0_DP
 S_E%VB(5,117)=2.602085213965211E-018_DP
 S_E%VA(5,118)=-8.673617379884035E-019_DP
 S_E%VB(5,118)=0.500000000000000E0_DP
 S_E%VB(5,119)=-3.469446951953614E-018_DP
 S_E%VA(5,120)=-8.673617379884035E-019_DP
 S_E%VA(5,132)=3.794707603699266E-018_DP
 S_E%VB(5,132)=2.25000000000000E0_DP
 S_E%VA(5,133)=0.666666666666667E0_DP
 S_E%VB(5,134)=-0.500000000000000E0_DP
 S_E%VB(5,136)=-6.938893903907228E-018_DP
 S_E%VA(5,148)=1.23333333333333E0_DP
 S_E%VB(5,148)=1.734723475976807E-018_DP
 S_E%VA(5,149)=-6.938893903907228E-018_DP
 S_E%VB(5,149)=-1.91666666666667E0_DP
 S_E%VA(5,150)=-0.666666666666667E0_DP
 S_E%VB(5,151)=0.500000000000000E0_DP
 S_E%VA(5,165)=2.168404344971009E-018_DP
 S_E%VB(5,165)=-2.16250000000000E0_DP
 S_E%VA(5,166)=-1.00000000000000E0_DP
 S_E%VA(5,167)=8.673617379884035E-019_DP
 S_E%VB(5,167)=1.62500000000000E0_DP
 S_E%VA(5,168)=0.666666666666667E0_DP
 S_E%VB(5,169)=-0.500000000000000E0_DP
 S_E%VB(5,171)=-0.200000000000000E0_DP
 S_E%VA(5,183)=-0.607142857142857E0_DP
 S_E%VA(5,184)=-3.469446951953614E-018_DP
 S_E%VB(5,184)=1.31250000000000E0_DP
 S_E%VA(5,185)=0.800000000000000E0_DP
 S_E%VB(5,186)=-1.37500000000000E0_DP
 S_E%VA(5,187)=-0.666666666666667E0_DP
 S_E%VB(5,188)=0.500000000000000E0_DP
 S_E%VA(5,189)=-1.00000000000000E0_DP
 S_E%VA(5,202)=-2.059984127722458E-018_DP
 S_E%VB(5,202)=0.620535714285715E0_DP
 S_E%VA(5,203)=0.321428571428572E0_DP
 S_E%VB(5,203)=1.734723475976807E-018_DP
 S_E%VA(5,204)=1.734723475976807E-018_DP
 S_E%VB(5,204)=-0.750000000000000E0_DP
 S_E%VA(5,205)=-0.633333333333333E0_DP
 S_E%VB(5,206)=1.16666666666667E0_DP
 S_E%VA(5,207)=0.666666666666667E0_DP
 S_E%VB(5,208)=2.00000000000000E0_DP
 S_E%VA(5,222)=8.283730158730163E-002_DP
 S_E%VB(5,222)=8.673617379884035E-019_DP
 S_E%VA(5,223)=8.673617379884035E-019_DP
 S_E%VB(5,223)=-0.223214285714286E0_DP
 S_E%VA(5,224)=-0.150000000000000E0_DP
 S_E%VA(5,225)=-1.734723475976807E-018_DP
 S_E%VB(5,225)=0.400000000000000E0_DP
 S_E%VA(5,226)=0.500000000000000E0_DP
 S_E%VB(5,227)=-1.00000000000000E0_DP
 S_E%VA(5,228)=2.00000000000000E0_DP
 S_E%VA(5,243)=5.963111948670274E-019_DP
 S_E%VB(5,243)=-3.571428571428573E-002_DP
 S_E%VA(5,244)=-2.182539682539683E-002_DP
 S_E%VB(5,244)=-4.336808689942018E-019_DP
 S_E%VA(5,245)=-4.336808689942018E-019_DP
 S_E%VB(5,245)=6.249999999999996E-002_DP
 S_E%VA(5,246)=5.714285714285713E-002_DP
 S_E%VA(5,247)=-1.734723475976807E-018_DP
 S_E%VB(5,247)=-0.200000000000000E0_DP
 S_E%VA(5,248)=-0.400000000000000E0_DP
 S_E%VB(5,249)=-1.00000000000000E0_DP
 S_E%VA(5,265)=-1.307720057720058E-003_DP
 S_E%VB(5,265)=2.032879073410321E-019_DP
 S_E%VA(5,266)=2.168404344971009E-019_DP
 S_E%VB(5,266)=4.265873015873016E-003_DP
 S_E%VA(5,267)=3.373015873015870E-003_DP
 S_E%VB(5,267)=-1.084202172485504E-019_DP
 S_E%VA(5,268)=-4.336808689942018E-019_DP
 S_E%VB(5,268)=-1.071428571428571E-002_DP
 S_E%VA(5,269)=-1.428571428571428E-002_DP
 S_E%VB(5,269)=-8.673617379884035E-019_DP
 S_E%VB(5,270)=0.100000000000000E0_DP
 S_E%VA(5,271)=-0.200000000000000E0_DP
 S_E%VA(6,78)=-2.710505431213761E-020_DP
 S_E%VA(6,90)=-6.938893903907228E-018_DP
 S_E%VA(6,91)=2.168404344971009E-019_DP
 S_E%VA(6,103)=-4.336808689942018E-018_DP
 S_E%VB(6,103)=0.500000000000000E0_DP
 S_E%VB(6,104)=6.938893903907228E-018_DP
 S_E%VA(6,105)=1.734723475976807E-018_DP
 S_E%VA(6,117)=0.833333333333333E0_DP
 S_E%VB(6,117)=1.040834085586084E-017_DP
 S_E%VA(6,118)=-3.469446951953614E-018_DP
 S_E%VB(6,118)=-0.500000000000000E0_DP
 S_E%VB(6,119)=-3.469446951953614E-018_DP
 S_E%VA(6,132)=1.734723475976807E-018_DP
 S_E%VB(6,132)=-2.75000000000000E0_DP
 S_E%VA(6,133)=-0.833333333333333E0_DP
 S_E%VB(6,133)=1.040834085586084E-017_DP
 S_E%VB(6,134)=0.500000000000000E0_DP
 S_E%VB(6,135)=6.938893903907228E-018_DP
 S_E%VB(6,136)=6.938893903907228E-018_DP
 S_E%VA(6,148)=-1.91666666666667E0_DP
 S_E%VB(6,148)=1.561251128379126E-017_DP
 S_E%VB(6,149)=2.41666666666667E0_DP
 S_E%VA(6,150)=0.833333333333333E0_DP
 S_E%VB(6,151)=-0.500000000000000E0_DP
 S_E%VB(6,153)=-0.166666666666667E0_DP
 S_E%VA(6,165)=7.806255641895632E-018_DP
 S_E%VB(6,165)=2.76250000000000E0_DP
 S_E%VA(6,166)=1.62500000000000E0_DP
 S_E%VA(6,167)=3.469446951953614E-018_DP
 S_E%VB(6,167)=-2.12500000000000E0_DP
 S_E%VA(6,168)=-0.833333333333333E0_DP
 S_E%VB(6,168)=6.938893903907228E-018_DP
 S_E%VB(6,169)=0.500000000000000E0_DP
 S_E%VA(6,170)=-1.00000000000000E0_DP
 S_E%VA(6,183)=0.937500000000000E0_DP
 S_E%VB(6,183)=3.469446951953614E-018_DP
 S_E%VA(6,184)=3.469446951953614E-018_DP
 S_E%VB(6,184)=-1.81250000000000E0_DP
 S_E%VA(6,185)=-1.37500000000000E0_DP
 S_E%VB(6,186)=1.87500000000000E0_DP
 S_E%VA(6,187)=0.833333333333333E0_DP
 S_E%VB(6,188)=2.50000000000000E0_DP
 S_E%VA(6,202)=1.734723475976807E-018_DP
 S_E%VB(6,202)=-0.758928571428572E0_DP
 S_E%VA(6,203)=-0.535714285714286E0_DP
 S_E%VB(6,203)=3.469446951953614E-018_DP
 S_E%VA(6,204)=6.938893903907228E-018_DP
 S_E%VB(6,204)=1.16666666666667E0_DP
 S_E%VA(6,205)=1.16666666666667E0_DP
 S_E%VB(6,206)=-1.66666666666667E0_DP
 S_E%VA(6,207)=3.33333333333333E0_DP
 S_E%VA(6,222)=-0.124007936507936E0_DP
 S_E%VB(6,222)=3.469446951953614E-018_DP
 S_E%VB(6,223)=0.290178571428572E0_DP
 S_E%VA(6,224)=0.285714285714286E0_DP
 S_E%VB(6,224)=-6.938893903907228E-018_DP
 S_E%VB(6,225)=-0.750000000000000E0_DP
 S_E%VA(6,226)=-1.00000000000000E0_DP
 S_E%VB(6,227)=-2.50000000000000E0_DP
 S_E%VA(6,243)=2.168404344971009E-019_DP
 S_E%VB(6,243)=4.315476190476190E-002_DP
 S_E%VA(6,244)=3.472222222222222E-002_DP
 S_E%VB(6,245)=-8.928571428571426E-002_DP
 S_E%VA(6,246)=-0.142857142857143E0_DP
 S_E%VB(6,246)=-6.938893903907228E-018_DP
 S_E%VB(6,247)=0.500000000000000E0_DP
 S_E%VA(6,248)=-1.00000000000000E0_DP
 S_E%VA(6,265)=1.939033189033189E-003_DP
 S_E%VB(6,265)=-4.336808689942018E-019_DP
 S_E%VA(6,266)=-2.168404344971009E-019_DP
 S_E%VB(6,266)=-5.456349206349210E-003_DP
 S_E%VA(6,267)=-5.952380952380953E-003_DP
 S_E%VB(6,267)=2.168404344971009E-019_DP
 S_E%VA(6,268)=-4.336808689942018E-019_DP
 S_E%VB(6,268)=1.785714285714286E-002_DP
 S_E%VA(6,269)=7.142857142857142E-002_DP
 S_E%VB(6,270)=0.166666666666667E0_DP
 S_E%VA(7,78)=3.252606517456513E-019_DP
 S_E%VB(7,90)=1.734723475976807E-018_DP
 S_E%VB(7,103)=-0.500000000000000E0_DP
 S_E%VA(7,104)=1.387778780781446E-017_DP
 S_E%VB(7,105)=3.469446951953614E-018_DP
 S_E%VA(7,117)=-1.00000000000000E0_DP
 S_E%VB(7,118)=0.500000000000000E0_DP
 S_E%VB(7,119)=6.938893903907228E-018_DP
 S_E%VA(7,132)=-6.938893903907228E-018_DP
 S_E%VB(7,132)=3.37500000000000E0_DP
 S_E%VA(7,133)=1.00000000000000E0_DP
 S_E%VB(7,133)=2.775557561562891E-017_DP
 S_E%VB(7,134)=-0.500000000000000E0_DP
 S_E%VB(7,136)=-0.142857142857143E0_DP
 S_E%VA(7,148)=2.90000000000000E0_DP
 S_E%VB(7,148)=-6.938893903907228E-018_DP
 S_E%VB(7,149)=-3.04166666666667E0_DP
 S_E%VA(7,150)=-1.00000000000000E0_DP
 S_E%VB(7,151)=0.500000000000000E0_DP
 S_E%VA(7,152)=-1.00000000000000E0_DP
 S_E%VA(7,165)=-1.387778780781446E-017_DP
 S_E%VB(7,165)=-3.82500000000000E0_DP
 S_E%VA(7,166)=-2.55000000000000E0_DP
 S_E%VB(7,167)=2.75000000000000E0_DP
 S_E%VA(7,168)=1.00000000000000E0_DP
 S_E%VB(7,169)=3.00000000000000E0_DP
 S_E%VA(7,183)=-1.55357142857143E0_DP
 S_E%VB(7,183)=-6.938893903907228E-018_DP
 S_E%VB(7,184)=2.75000000000000E0_DP
 S_E%VA(7,185)=2.25000000000000E0_DP
 S_E%VB(7,186)=-2.50000000000000E0_DP
 S_E%VA(7,187)=5.00000000000000E0_DP
 S_E%VA(7,202)=3.469446951953614E-018_DP
 S_E%VB(7,202)=1.02678571428571E0_DP
 S_E%VA(7,203)=1.00000000000000E0_DP
 S_E%VB(7,204)=-2.00000000000000E0_DP
 S_E%VA(7,205)=-2.00000000000000E0_DP
 S_E%VB(7,206)=-5.00000000000000E0_DP
 S_E%VA(7,222)=0.193452380952381E0_DP
 S_E%VB(7,222)=-3.469446951953614E-018_DP
 S_E%VB(7,223)=-0.446428571428571E0_DP
 S_E%VA(7,224)=-0.642857142857143E0_DP
 S_E%VB(7,225)=1.50000000000000E0_DP
 S_E%VA(7,226)=-3.00000000000000E0_DP
 S_E%VA(7,243)=-4.336808689942018E-019_DP
 S_E%VB(7,243)=-5.654761904761905E-002_DP
 S_E%VA(7,244)=-5.952380952380951E-002_DP
 S_E%VB(7,244)=1.734723475976807E-018_DP
 S_E%VB(7,245)=0.178571428571429E0_DP
 S_E%VA(7,246)=0.428571428571429E0_DP
 S_E%VB(7,247)=1.00000000000000E0_DP
 S_E%VA(7,265)=-2.976190476190478E-003_DP
 S_E%VB(7,265)=5.421010862427522E-019_DP
 S_E%VB(7,266)=7.936507936507933E-003_DP
 S_E%VA(7,267)=1.190476190476191E-002_DP
 S_E%VB(7,268)=-7.142857142857142E-002_DP
 S_E%VA(7,269)=0.142857142857143E0_DP
 S_E%VB(8,90)=3.469446951953614E-018_DP
 S_E%VA(8,103)=-1.387778780781446E-017_DP
 S_E%VB(8,103)=0.500000000000000E0_DP
 S_E%VA(8,117)=1.16666666666667E0_DP
 S_E%VB(8,117)=-2.775557561562891E-017_DP
 S_E%VB(8,118)=-0.500000000000000E0_DP
 S_E%VB(8,120)=-0.125000000000000E0_DP
 S_E%VB(8,132)=-4.12500000000000E0_DP
 S_E%VA(8,133)=-1.16666666666667E0_DP
 S_E%VB(8,134)=0.500000000000000E0_DP
 S_E%VA(8,135)=-1.00000000000000E0_DP
 S_E%VA(8,148)=-4.25833333333333E0_DP
 S_E%VB(8,148)=-6.938893903907228E-018_DP
 S_E%VB(8,149)=3.79166666666667E0_DP
 S_E%VA(8,150)=1.16666666666667E0_DP
 S_E%VB(8,151)=3.50000000000000E0_DP
 S_E%VA(8,165)=-2.775557561562891E-017_DP
 S_E%VB(8,165)=5.60000000000000E0_DP
 S_E%VA(8,166)=3.85000000000000E0_DP
 S_E%VB(8,167)=-3.50000000000000E0_DP
 S_E%VA(8,168)=7.00000000000000E0_DP
 S_E%VA(8,183)=2.75000000000000E0_DP
 S_E%VB(8,183)=1.387778780781446E-017_DP
 S_E%VB(8,184)=-4.37500000000000E0_DP
 S_E%VA(8,185)=-3.50000000000000E0_DP
 S_E%VB(8,186)=-8.75000000000000E0_DP
 S_E%VB(8,202)=-1.62500000000000E0_DP
 S_E%VA(8,203)=-2.00000000000000E0_DP
 S_E%VB(8,204)=3.50000000000000E0_DP
 S_E%VA(8,205)=-7.00000000000000E0_DP
 S_E%VA(8,222)=-0.347222222222222E0_DP
 S_E%VB(8,223)=0.875000000000000E0_DP
 S_E%VA(8,224)=1.50000000000000E0_DP
 S_E%VB(8,225)=3.50000000000000E0_DP
 S_E%VB(8,243)=8.333333333333330E-002_DP
 S_E%VA(8,244)=0.138888888888889E0_DP
 S_E%VB(8,244)=3.469446951953614E-018_DP
 S_E%VB(8,245)=-0.500000000000000E0_DP
 S_E%VA(8,246)=1.00000000000000E0_DP
 S_E%VA(8,265)=5.050505050505052E-003_DP
 S_E%VA(8,266)=4.336808689942018E-019_DP
 S_E%VB(8,266)=-1.388888888888889E-002_DP
 S_E%VA(8,267)=-5.555555555555555E-002_DP
 S_E%VB(8,268)=-0.125000000000000E0_DP
 S_E%VB(9,90)=-6.938893903907228E-018_DP
 S_E%VB(9,103)=-0.500000000000000E0_DP
 S_E%VB(9,105)=-0.111111111111111E0_DP
 S_E%VA(9,117)=-1.33333333333333E0_DP
 S_E%VB(9,118)=0.500000000000000E0_DP
 S_E%VA(9,119)=-1.00000000000000E0_DP
 S_E%VB(9,132)=5.00000000000000E0_DP
 S_E%VA(9,133)=1.33333333333333E0_DP
 S_E%VB(9,134)=4.00000000000000E0_DP
 S_E%VA(9,148)=6.06666666666667E0_DP
 S_E%VB(9,149)=-4.66666666666667E0_DP
 S_E%VA(9,150)=9.33333333333333E0_DP
 S_E%VB(9,165)=-8.40000000000000E0_DP
 S_E%VA(9,166)=-5.60000000000000E0_DP
 S_E%VB(9,167)=-14.0000000000000E0_DP
 S_E%VA(9,183)=-5.00000000000000E0_DP
 S_E%VB(9,184)=7.00000000000000E0_DP
 S_E%VA(9,185)=-14.0000000000000E0_DP
 S_E%VA(9,202)=2.775557561562891E-017_DP
 S_E%VB(9,202)=3.00000000000000E0_DP
 S_E%VA(9,203)=4.00000000000000E0_DP
 S_E%VB(9,204)=9.33333333333333E0_DP
 S_E%VA(9,222)=0.777777777777778E0_DP
 S_E%VB(9,223)=-2.00000000000000E0_DP
 S_E%VA(9,224)=4.00000000000000E0_DP
 S_E%VB(9,243)=-0.166666666666667E0_DP
 S_E%VA(9,244)=-0.444444444444444E0_DP
 S_E%VB(9,245)=-1.00000000000000E0_DP
 S_E%VA(9,265)=-1.010101010101010E-002_DP
 S_E%VB(9,266)=5.555555555555555E-002_DP
 S_E%VA(9,267)=-0.111111111111111E0_DP
 S_E%VB(10,91)=-0.100000000000000E0_DP
 S_E%VB(10,103)=0.500000000000000E0_DP
 S_E%VA(10,104)=-1.00000000000000E0_DP
 S_E%VA(10,117)=1.50000000000000E0_DP
 S_E%VB(10,118)=4.50000000000000E0_DP
 S_E%VB(10,132)=-6.00000000000000E0_DP
 S_E%VA(10,133)=12.0000000000000E0_DP
 S_E%VA(10,148)=-8.40000000000000E0_DP
 S_E%VB(10,149)=-21.0000000000000E0_DP
 S_E%VB(10,165)=12.6000000000000E0_DP
 S_E%VA(10,166)=-25.2000000000000E0_DP
 S_E%VA(10,183)=9.00000000000000E0_DP
 S_E%VB(10,184)=21.0000000000000E0_DP
 S_E%VB(10,202)=-6.00000000000000E0_DP
 S_E%VA(10,203)=12.0000000000000E0_DP
 S_E%VA(10,222)=-2.00000000000000E0_DP
 S_E%VB(10,223)=-4.50000000000000E0_DP
 S_E%VB(10,243)=0.500000000000000E0_DP
 S_E%VA(10,244)=-1.00000000000000E0_DP
 S_E%VA(10,265)=4.545454545454546E-002_DP
 S_E%VB(10,266)=0.100000000000000E0_DP
 S_E%VB(11,78)=-9.090909090909091E-002_DP
 S_E%VA(11,90)=-1.00000000000000E0_DP
 S_E%VB(11,103)=5.00000000000000E0_DP
 S_E%VA(11,117)=15.0000000000000E0_DP
 S_E%VB(11,132)=-30.0000000000000E0_DP
 S_E%VA(11,148)=-42.0000000000000E0_DP
 S_E%VB(11,165)=42.0000000000000E0_DP
 S_E%VA(11,183)=30.0000000000000E0_DP
 S_E%VB(11,202)=-15.0000000000000E0_DP
 S_E%VA(11,222)=-5.00000000000000E0_DP
 S_E%VB(11,243)=1.00000000000000E0_DP
 S_E%VA(11,265)=9.090909090909091E-002_DP
 S_E%VB(12,66)=-8.333333333333333E-002_DP
 S_E%VA(12,77)=-1.00000000000000E0_DP
 S_E%VB(12,89)=5.50000000000000E0_DP
 S_E%VA(12,102)=18.3333333333333E0_DP
 S_E%VB(12,116)=-41.2500000000000E0_DP
 S_E%VA(12,131)=-66.0000000000000E0_DP
 S_E%VB(12,147)=77.0000000000000E0_DP
 S_E%VA(12,164)=66.0000000000000E0_DP
 S_E%VB(12,182)=-41.2500000000000E0_DP
 S_E%VA(12,201)=-18.3333333333333E0_DP
 S_E%VB(12,221)=5.50000000000000E0_DP
 S_E%VA(12,242)=1.00000000000000E0_DP
 S_E%VB(12,264)=-8.333333333333333E-002_DP
 S_E%VB(13,55)=-7.692307692307693E-002_DP
 S_E%VA(13,65)=-1.00000000000000E0_DP
 S_E%VB(13,76)=6.00000000000000E0_DP
 S_E%VA(13,88)=22.0000000000000E0_DP
 S_E%VB(13,101)=-55.0000000000000E0_DP
 S_E%VA(13,115)=-99.0000000000000E0_DP
 S_E%VB(13,130)=132.000000000000E0_DP
 S_E%VA(13,146)=132.000000000000E0_DP
 S_E%VB(13,163)=-99.0000000000000E0_DP
 S_E%VA(13,181)=-55.0000000000000E0_DP
 S_E%VB(13,200)=22.0000000000000E0_DP
 S_E%VA(13,220)=6.00000000000000E0_DP
 S_E%VB(13,241)=-1.00000000000000E0_DP
 S_E%VA(13,263)=-7.692307692307693E-002_DP
 S_E%VB(14,45)=-7.142857142857142E-002_DP
 S_E%VA(14,54)=-1.00000000000000E0_DP
 S_E%VB(14,64)=6.50000000000000E0_DP
 S_E%VA(14,75)=26.0000000000000E0_DP
 S_E%VB(14,87)=-71.5000000000000E0_DP
 S_E%VA(14,100)=-143.000000000000E0_DP
 S_E%VB(14,114)=214.500000000000E0_DP
 S_E%VA(14,129)=245.142857142857E0_DP
 S_E%VB(14,145)=-214.500000000000E0_DP
 S_E%VA(14,162)=-143.000000000000E0_DP
 S_E%VB(14,180)=71.5000000000000E0_DP
 S_E%VA(14,199)=26.0000000000000E0_DP
 S_E%VB(14,219)=-6.50000000000000E0_DP
 S_E%VA(14,240)=-1.00000000000000E0_DP
 S_E%VB(14,262)=7.142857142857142E-002_DP
 S_E%VB(15,36)=-6.666666666666667E-002_DP
 S_E%VA(15,44)=-1.00000000000000E0_DP
 S_E%VB(15,53)=7.00000000000000E0_DP
 S_E%VA(15,63)=30.3333333333333E0_DP
 S_E%VB(15,74)=-91.0000000000000E0_DP
 S_E%VA(15,86)=-200.200000000000E0_DP
 S_E%VB(15,99)=333.666666666667E0_DP
 S_E%VA(15,113)=429.000000000000E0_DP
 S_E%VB(15,128)=-429.000000000000E0_DP
 S_E%VA(15,144)=-333.666666666667E0_DP
 S_E%VB(15,161)=200.200000000000E0_DP
 S_E%VA(15,179)=91.0000000000000E0_DP
 S_E%VB(15,198)=-30.3333333333333E0_DP
 S_E%VA(15,218)=-7.00000000000000E0_DP
 S_E%VB(15,239)=1.00000000000000E0_DP
 S_E%VA(15,261)=6.666666666666667E-002_DP
 S_E%VB(16,28)=-6.250000000000000E-002_DP
 S_E%VA(16,35)=-1.00000000000000E0_DP
 S_E%VB(16,43)=7.50000000000000E0_DP
 S_E%VA(16,52)=35.0000000000000E0_DP
 S_E%VB(16,62)=-113.750000000000E0_DP
 S_E%VA(16,73)=-273.000000000000E0_DP
 S_E%VB(16,85)=500.500000000000E0_DP
 S_E%VA(16,98)=715.000000000000E0_DP
 S_E%VB(16,112)=-804.375000000000E0_DP
 S_E%VA(16,127)=-715.000000000000E0_DP
 S_E%VB(16,143)=500.500000000000E0_DP
 S_E%VA(16,160)=273.000000000000E0_DP
 S_E%VB(16,178)=-113.750000000000E0_DP
 S_E%VA(16,197)=-35.0000000000000E0_DP
 S_E%VB(16,217)=7.50000000000000E0_DP
 S_E%VA(16,238)=1.00000000000000E0_DP
 S_E%VB(16,260)=-6.250000000000000E-002_DP
 S_E%VB(17,21)=-5.882352941176471E-002_DP
 S_E%VA(17,27)=-1.00000000000000E0_DP
 S_E%VB(17,34)=8.00000000000000E0_DP
 S_E%VA(17,42)=40.0000000000000E0_DP
 S_E%VB(17,51)=-140.000000000000E0_DP
 S_E%VA(17,61)=-364.000000000000E0_DP
 S_E%VB(17,72)=728.000000000000E0_DP
 S_E%VA(17,84)=1144.00000000000E0_DP
 S_E%VB(17,97)=-1430.00000000000E0_DP
 S_E%VA(17,111)=-1430.00000000000E0_DP
 S_E%VB(17,126)=1144.00000000000E0_DP
 S_E%VA(17,142)=728.000000000000E0_DP
 S_E%VB(17,159)=-364.000000000000E0_DP
 S_E%VA(17,177)=-140.000000000000E0_DP
 S_E%VB(17,196)=40.0000000000000E0_DP
 S_E%VA(17,216)=8.00000000000000E0_DP
 S_E%VB(17,237)=-1.00000000000000E0_DP
 S_E%VA(17,259)=-5.882352941176471E-002_DP
 S_E%VB(18,15)=-5.555555555555555E-002_DP
 S_E%VA(18,20)=-1.00000000000000E0_DP
 S_E%VB(18,26)=8.50000000000000E0_DP
 S_E%VA(18,33)=45.3333333333333E0_DP
 S_E%VB(18,41)=-170.000000000000E0_DP
 S_E%VA(18,50)=-476.000000000000E0_DP
 S_E%VB(18,60)=1031.33333333333E0_DP
 S_E%VA(18,71)=1768.00000000000E0_DP
 S_E%VB(18,83)=-2431.00000000000E0_DP
 S_E%VA(18,96)=-2701.11111111111E0_DP
 S_E%VB(18,110)=2431.00000000000E0_DP
 S_E%VA(18,125)=1768.00000000000E0_DP
 S_E%VB(18,141)=-1031.33333333333E0_DP
 S_E%VA(18,158)=-476.000000000000E0_DP
 S_E%VB(18,176)=170.000000000000E0_DP
 S_E%VA(18,195)=45.3333333333333E0_DP
 S_E%VB(18,215)=-8.50000000000000E0_DP
 S_E%VA(18,236)=-1.00000000000000E0_DP
 S_E%VB(18,258)=5.555555555555555E-002_DP
 S_E%VB(19,10)=-5.263157894736842E-002_DP
 S_E%VA(19,14)=-1.00000000000000E0_DP
 S_E%VB(19,19)=9.00000000000000E0_DP
 S_E%VA(19,25)=51.0000000000000E0_DP
 S_E%VB(19,32)=-204.000000000000E0_DP
 S_E%VA(19,40)=-612.000000000000E0_DP
 S_E%VB(19,49)=1428.00000000000E0_DP
 S_E%VA(19,59)=2652.00000000000E0_DP
 S_E%VB(19,70)=-3978.00000000000E0_DP
 S_E%VA(19,82)=-4862.00000000000E0_DP
 S_E%VB(19,95)=4862.00000000000E0_DP
 S_E%VA(19,109)=3978.00000000000E0_DP
 S_E%VB(19,124)=-2652.00000000000E0_DP
 S_E%VA(19,140)=-1428.00000000000E0_DP
 S_E%VB(19,157)=612.000000000000E0_DP
 S_E%VA(19,175)=204.000000000000E0_DP
 S_E%VB(19,194)=-51.0000000000000E0_DP
 S_E%VA(19,214)=-9.00000000000000E0_DP
 S_E%VB(19,235)=1.00000000000000E0_DP
 S_E%VA(19,257)=5.263157894736842E-002_DP
 S_E%VB(20,6)=-5.000000000000000E-002_DP
 S_E%VA(20,9)=-1.00000000000000E0_DP
 S_E%VB(20,13)=9.50000000000000E0_DP
 S_E%VA(20,18)=57.0000000000000E0_DP
 S_E%VB(20,24)=-242.250000000000E0_DP
 S_E%VA(20,31)=-775.200000000000E0_DP
 S_E%VB(20,39)=1938.00000000000E0_DP
 S_E%VA(20,48)=3876.00000000000E0_DP
 S_E%VB(20,58)=-6298.50000000000E0_DP
 S_E%VA(20,69)=-8398.00000000000E0_DP
 S_E%VB(20,81)=9237.80000000000E0_DP
 S_E%VA(20,94)=8398.00000000000E0_DP
 S_E%VB(20,108)=-6298.50000000000E0_DP
 S_E%VA(20,123)=-3876.00000000000E0_DP
 S_E%VB(20,139)=1938.00000000000E0_DP
 S_E%VA(20,156)=775.200000000000E0_DP
 S_E%VB(20,174)=-242.250000000000E0_DP
 S_E%VA(20,193)=-57.0000000000000E0_DP
 S_E%VB(20,213)=9.50000000000000E0_DP
 S_E%VA(20,234)=1.00000000000000E0_DP
 S_E%VB(20,256)=-5.000000000000000E-002_DP
 S_E%VB(21,3)=-4.761904761904762E-002_DP
 S_E%VA(21,5)=-1.00000000000000E0_DP
 S_E%VB(21,8)=10.0000000000000E0_DP
 S_E%VA(21,12)=63.3333333333333E0_DP
 S_E%VB(21,17)=-285.000000000000E0_DP
 S_E%VA(21,23)=-969.000000000000E0_DP
 S_E%VB(21,30)=2584.00000000000E0_DP
 S_E%VA(21,38)=5537.14285714286E0_DP
 S_E%VB(21,47)=-9690.00000000000E0_DP
 S_E%VA(21,57)=-13996.6666666667E0_DP
 S_E%VB(21,68)=16796.0000000000E0_DP
 S_E%VA(21,80)=16796.0000000000E0_DP
 S_E%VB(21,93)=-13996.6666666667E0_DP
 S_E%VA(21,107)=-9690.00000000000E0_DP
 S_E%VB(21,122)=5537.14285714286E0_DP
 S_E%VA(21,138)=2584.00000000000E0_DP
 S_E%VB(21,155)=-969.000000000000E0_DP
 S_E%VA(21,173)=-285.000000000000E0_DP
 S_E%VB(21,192)=63.3333333333333E0_DP
 S_E%VA(21,212)=10.0000000000000E0_DP
 S_E%VB(21,233)=-1.00000000000000E0_DP
 S_E%VA(21,255)=-4.761904761904762E-002_DP
 S_E%VB(22,1)=-4.545454545454546E-002_DP
 S_E%VA(22,2)=-1.00000000000000E0_DP
 S_E%VB(22,4)=10.5000000000000E0_DP
 S_E%VA(22,7)=70.0000000000000E0_DP
 S_E%VB(22,11)=-332.500000000000E0_DP
 S_E%VA(22,16)=-1197.00000000000E0_DP
 S_E%VB(22,22)=3391.50000000000E0_DP
 S_E%VA(22,29)=7752.00000000000E0_DP
 S_E%VB(22,37)=-14535.0000000000E0_DP
 S_E%VA(22,46)=-22610.0000000000E0_DP
 S_E%VB(22,56)=29393.0000000000E0_DP
 S_E%VA(22,67)=32065.0909090909E0_DP
 S_E%VB(22,79)=-29393.0000000000E0_DP
 S_E%VA(22,92)=-22610.0000000000E0_DP
 S_E%VB(22,106)=14535.0000000000E0_DP
 S_E%VA(22,121)=7752.00000000000E0_DP
 S_E%VB(22,137)=-3391.50000000000E0_DP
 S_E%VA(22,154)=-1197.00000000000E0_DP
 S_E%VB(22,172)=332.500000000000E0_DP
 S_E%VA(22,191)=70.0000000000000E0_DP
 S_E%VB(22,211)=-10.5000000000000E0_DP
 S_E%VA(22,232)=-1.00000000000000E0_DP
 S_E%VB(22,254)=4.545454545454546E-002_DP

end subroutine set_s_e 

subroutine set_s_b_mcmillan
  implicit none
  integer i
 


  !ALLOCATE(S_B(SECTOR_NMUL_MAX))
i=SECTOR_NMUL_MAX
          ! DO I=1,SECTOR_NMUL_MAX

             S_B_from_V%firsttime = 1  ! Piotr 8.19.2014 
             call nul_coef(S_B_from_V)
             call make_set_coef(S_B_from_V,I,0)
             S_B_from_V%firsttime=0
        !  ENDDO
 
 S_B_FROM_V%I(1)=22;S_B_FROM_V%J(1)=0;
 S_B_FROM_V%I(2)=21;S_B_FROM_V%J(2)=1;
 S_B_FROM_V%I(3)=21;S_B_FROM_V%J(3)=0;
 S_B_FROM_V%I(4)=20;S_B_FROM_V%J(4)=2;
 S_B_FROM_V%I(5)=20;S_B_FROM_V%J(5)=1;
 S_B_FROM_V%I(6)=20;S_B_FROM_V%J(6)=0;
 S_B_FROM_V%I(7)=19;S_B_FROM_V%J(7)=3;
 S_B_FROM_V%I(8)=19;S_B_FROM_V%J(8)=2;
 S_B_FROM_V%I(9)=19;S_B_FROM_V%J(9)=1;
 S_B_FROM_V%I(10)=19;S_B_FROM_V%J(10)=0;
 S_B_FROM_V%I(11)=18;S_B_FROM_V%J(11)=4;
 S_B_FROM_V%I(12)=18;S_B_FROM_V%J(12)=3;
 S_B_FROM_V%I(13)=18;S_B_FROM_V%J(13)=2;
 S_B_FROM_V%I(14)=18;S_B_FROM_V%J(14)=1;
 S_B_FROM_V%I(15)=18;S_B_FROM_V%J(15)=0;
 S_B_FROM_V%I(16)=17;S_B_FROM_V%J(16)=5;
 S_B_FROM_V%I(17)=17;S_B_FROM_V%J(17)=4;
 S_B_FROM_V%I(18)=17;S_B_FROM_V%J(18)=3;
 S_B_FROM_V%I(19)=17;S_B_FROM_V%J(19)=2;
 S_B_FROM_V%I(20)=17;S_B_FROM_V%J(20)=1;
 S_B_FROM_V%I(21)=17;S_B_FROM_V%J(21)=0;
 S_B_FROM_V%I(22)=16;S_B_FROM_V%J(22)=6;
 S_B_FROM_V%I(23)=16;S_B_FROM_V%J(23)=5;
 S_B_FROM_V%I(24)=16;S_B_FROM_V%J(24)=4;
 S_B_FROM_V%I(25)=16;S_B_FROM_V%J(25)=3;
 S_B_FROM_V%I(26)=16;S_B_FROM_V%J(26)=2;
 S_B_FROM_V%I(27)=16;S_B_FROM_V%J(27)=1;
 S_B_FROM_V%I(28)=16;S_B_FROM_V%J(28)=0;
 S_B_FROM_V%I(29)=15;S_B_FROM_V%J(29)=7;
 S_B_FROM_V%I(30)=15;S_B_FROM_V%J(30)=6;
 S_B_FROM_V%I(31)=15;S_B_FROM_V%J(31)=5;
 S_B_FROM_V%I(32)=15;S_B_FROM_V%J(32)=4;
 S_B_FROM_V%I(33)=15;S_B_FROM_V%J(33)=3;
 S_B_FROM_V%I(34)=15;S_B_FROM_V%J(34)=2;
 S_B_FROM_V%I(35)=15;S_B_FROM_V%J(35)=1;
 S_B_FROM_V%I(36)=15;S_B_FROM_V%J(36)=0;
 S_B_FROM_V%I(37)=14;S_B_FROM_V%J(37)=8;
 S_B_FROM_V%I(38)=14;S_B_FROM_V%J(38)=7;
 S_B_FROM_V%I(39)=14;S_B_FROM_V%J(39)=6;
 S_B_FROM_V%I(40)=14;S_B_FROM_V%J(40)=5;
 S_B_FROM_V%I(41)=14;S_B_FROM_V%J(41)=4;
 S_B_FROM_V%I(42)=14;S_B_FROM_V%J(42)=3;
 S_B_FROM_V%I(43)=14;S_B_FROM_V%J(43)=2;
 S_B_FROM_V%I(44)=14;S_B_FROM_V%J(44)=1;
 S_B_FROM_V%I(45)=14;S_B_FROM_V%J(45)=0;
 S_B_FROM_V%I(46)=13;S_B_FROM_V%J(46)=9;
 S_B_FROM_V%I(47)=13;S_B_FROM_V%J(47)=8;
 S_B_FROM_V%I(48)=13;S_B_FROM_V%J(48)=7;
 S_B_FROM_V%I(49)=13;S_B_FROM_V%J(49)=6;
 S_B_FROM_V%I(50)=13;S_B_FROM_V%J(50)=5;
 S_B_FROM_V%I(51)=13;S_B_FROM_V%J(51)=4;
 S_B_FROM_V%I(52)=13;S_B_FROM_V%J(52)=3;
 S_B_FROM_V%I(53)=13;S_B_FROM_V%J(53)=2;
 S_B_FROM_V%I(54)=13;S_B_FROM_V%J(54)=1;
 S_B_FROM_V%I(55)=13;S_B_FROM_V%J(55)=0;
 S_B_FROM_V%I(56)=12;S_B_FROM_V%J(56)=10;
 S_B_FROM_V%I(57)=12;S_B_FROM_V%J(57)=9;
 S_B_FROM_V%I(58)=12;S_B_FROM_V%J(58)=8;
 S_B_FROM_V%I(59)=12;S_B_FROM_V%J(59)=7;
 S_B_FROM_V%I(60)=12;S_B_FROM_V%J(60)=6;
 S_B_FROM_V%I(61)=12;S_B_FROM_V%J(61)=5;
 S_B_FROM_V%I(62)=12;S_B_FROM_V%J(62)=4;
 S_B_FROM_V%I(63)=12;S_B_FROM_V%J(63)=3;
 S_B_FROM_V%I(64)=12;S_B_FROM_V%J(64)=2;
 S_B_FROM_V%I(65)=12;S_B_FROM_V%J(65)=1;
 S_B_FROM_V%I(66)=12;S_B_FROM_V%J(66)=0;
 S_B_FROM_V%I(67)=11;S_B_FROM_V%J(67)=11;
 S_B_FROM_V%I(68)=11;S_B_FROM_V%J(68)=10;
 S_B_FROM_V%I(69)=11;S_B_FROM_V%J(69)=9;
 S_B_FROM_V%I(70)=11;S_B_FROM_V%J(70)=8;
 S_B_FROM_V%I(71)=11;S_B_FROM_V%J(71)=7;
 S_B_FROM_V%I(72)=11;S_B_FROM_V%J(72)=6;
 S_B_FROM_V%I(73)=11;S_B_FROM_V%J(73)=5;
 S_B_FROM_V%I(74)=11;S_B_FROM_V%J(74)=4;
 S_B_FROM_V%I(75)=11;S_B_FROM_V%J(75)=3;
 S_B_FROM_V%I(76)=11;S_B_FROM_V%J(76)=2;
 S_B_FROM_V%I(77)=11;S_B_FROM_V%J(77)=1;
 S_B_FROM_V%I(78)=11;S_B_FROM_V%J(78)=0;
 S_B_FROM_V%I(79)=10;S_B_FROM_V%J(79)=12;
 S_B_FROM_V%I(80)=10;S_B_FROM_V%J(80)=11;
 S_B_FROM_V%I(81)=10;S_B_FROM_V%J(81)=10;
 S_B_FROM_V%I(82)=10;S_B_FROM_V%J(82)=9;
 S_B_FROM_V%I(83)=10;S_B_FROM_V%J(83)=8;
 S_B_FROM_V%I(84)=10;S_B_FROM_V%J(84)=7;
 S_B_FROM_V%I(85)=10;S_B_FROM_V%J(85)=6;
 S_B_FROM_V%I(86)=10;S_B_FROM_V%J(86)=5;
 S_B_FROM_V%I(87)=10;S_B_FROM_V%J(87)=4;
 S_B_FROM_V%I(88)=10;S_B_FROM_V%J(88)=3;
 S_B_FROM_V%I(89)=10;S_B_FROM_V%J(89)=2;
 S_B_FROM_V%I(90)=10;S_B_FROM_V%J(90)=1;
 S_B_FROM_V%I(91)=10;S_B_FROM_V%J(91)=0;
 S_B_FROM_V%I(92)=9;S_B_FROM_V%J(92)=13;
 S_B_FROM_V%I(93)=9;S_B_FROM_V%J(93)=12;
 S_B_FROM_V%I(94)=9;S_B_FROM_V%J(94)=11;
 S_B_FROM_V%I(95)=9;S_B_FROM_V%J(95)=10;
 S_B_FROM_V%I(96)=9;S_B_FROM_V%J(96)=9;
 S_B_FROM_V%I(97)=9;S_B_FROM_V%J(97)=8;
 S_B_FROM_V%I(98)=9;S_B_FROM_V%J(98)=7;
 S_B_FROM_V%I(99)=9;S_B_FROM_V%J(99)=6;
 S_B_FROM_V%I(100)=9;S_B_FROM_V%J(100)=5;
 S_B_FROM_V%I(101)=9;S_B_FROM_V%J(101)=4;
 S_B_FROM_V%I(102)=9;S_B_FROM_V%J(102)=3;
 S_B_FROM_V%I(103)=9;S_B_FROM_V%J(103)=2;
 S_B_FROM_V%I(104)=9;S_B_FROM_V%J(104)=1;
 S_B_FROM_V%I(105)=9;S_B_FROM_V%J(105)=0;
 S_B_FROM_V%I(106)=8;S_B_FROM_V%J(106)=14;
 S_B_FROM_V%I(107)=8;S_B_FROM_V%J(107)=13;
 S_B_FROM_V%I(108)=8;S_B_FROM_V%J(108)=12;
 S_B_FROM_V%I(109)=8;S_B_FROM_V%J(109)=11;
 S_B_FROM_V%I(110)=8;S_B_FROM_V%J(110)=10;
 S_B_FROM_V%I(111)=8;S_B_FROM_V%J(111)=9;
 S_B_FROM_V%I(112)=8;S_B_FROM_V%J(112)=8;
 S_B_FROM_V%I(113)=8;S_B_FROM_V%J(113)=7;
 S_B_FROM_V%I(114)=8;S_B_FROM_V%J(114)=6;
 S_B_FROM_V%I(115)=8;S_B_FROM_V%J(115)=5;
 S_B_FROM_V%I(116)=8;S_B_FROM_V%J(116)=4;
 S_B_FROM_V%I(117)=8;S_B_FROM_V%J(117)=3;
 S_B_FROM_V%I(118)=8;S_B_FROM_V%J(118)=2;
 S_B_FROM_V%I(119)=8;S_B_FROM_V%J(119)=1;
 S_B_FROM_V%I(120)=8;S_B_FROM_V%J(120)=0;
 S_B_FROM_V%I(121)=7;S_B_FROM_V%J(121)=15;
 S_B_FROM_V%I(122)=7;S_B_FROM_V%J(122)=14;
 S_B_FROM_V%I(123)=7;S_B_FROM_V%J(123)=13;
 S_B_FROM_V%I(124)=7;S_B_FROM_V%J(124)=12;
 S_B_FROM_V%I(125)=7;S_B_FROM_V%J(125)=11;
 S_B_FROM_V%I(126)=7;S_B_FROM_V%J(126)=10;
 S_B_FROM_V%I(127)=7;S_B_FROM_V%J(127)=9;
 S_B_FROM_V%I(128)=7;S_B_FROM_V%J(128)=8;
 S_B_FROM_V%I(129)=7;S_B_FROM_V%J(129)=7;
 S_B_FROM_V%I(130)=7;S_B_FROM_V%J(130)=6;
 S_B_FROM_V%I(131)=7;S_B_FROM_V%J(131)=5;
 S_B_FROM_V%I(132)=7;S_B_FROM_V%J(132)=4;
 S_B_FROM_V%I(133)=7;S_B_FROM_V%J(133)=3;
 S_B_FROM_V%I(134)=7;S_B_FROM_V%J(134)=2;
 S_B_FROM_V%I(135)=7;S_B_FROM_V%J(135)=1;
 S_B_FROM_V%I(136)=7;S_B_FROM_V%J(136)=0;
 S_B_FROM_V%I(137)=6;S_B_FROM_V%J(137)=16;
 S_B_FROM_V%I(138)=6;S_B_FROM_V%J(138)=15;
 S_B_FROM_V%I(139)=6;S_B_FROM_V%J(139)=14;
 S_B_FROM_V%I(140)=6;S_B_FROM_V%J(140)=13;
 S_B_FROM_V%I(141)=6;S_B_FROM_V%J(141)=12;
 S_B_FROM_V%I(142)=6;S_B_FROM_V%J(142)=11;
 S_B_FROM_V%I(143)=6;S_B_FROM_V%J(143)=10;
 S_B_FROM_V%I(144)=6;S_B_FROM_V%J(144)=9;
 S_B_FROM_V%I(145)=6;S_B_FROM_V%J(145)=8;
 S_B_FROM_V%I(146)=6;S_B_FROM_V%J(146)=7;
 S_B_FROM_V%I(147)=6;S_B_FROM_V%J(147)=6;
 S_B_FROM_V%I(148)=6;S_B_FROM_V%J(148)=5;
 S_B_FROM_V%I(149)=6;S_B_FROM_V%J(149)=4;
 S_B_FROM_V%I(150)=6;S_B_FROM_V%J(150)=3;
 S_B_FROM_V%I(151)=6;S_B_FROM_V%J(151)=2;
 S_B_FROM_V%I(152)=6;S_B_FROM_V%J(152)=1;
 S_B_FROM_V%I(153)=6;S_B_FROM_V%J(153)=0;
 S_B_FROM_V%I(154)=5;S_B_FROM_V%J(154)=17;
 S_B_FROM_V%I(155)=5;S_B_FROM_V%J(155)=16;
 S_B_FROM_V%I(156)=5;S_B_FROM_V%J(156)=15;
 S_B_FROM_V%I(157)=5;S_B_FROM_V%J(157)=14;
 S_B_FROM_V%I(158)=5;S_B_FROM_V%J(158)=13;
 S_B_FROM_V%I(159)=5;S_B_FROM_V%J(159)=12;
 S_B_FROM_V%I(160)=5;S_B_FROM_V%J(160)=11;
 S_B_FROM_V%I(161)=5;S_B_FROM_V%J(161)=10;
 S_B_FROM_V%I(162)=5;S_B_FROM_V%J(162)=9;
 S_B_FROM_V%I(163)=5;S_B_FROM_V%J(163)=8;
 S_B_FROM_V%I(164)=5;S_B_FROM_V%J(164)=7;
 S_B_FROM_V%I(165)=5;S_B_FROM_V%J(165)=6;
 S_B_FROM_V%I(166)=5;S_B_FROM_V%J(166)=5;
 S_B_FROM_V%I(167)=5;S_B_FROM_V%J(167)=4;
 S_B_FROM_V%I(168)=5;S_B_FROM_V%J(168)=3;
 S_B_FROM_V%I(169)=5;S_B_FROM_V%J(169)=2;
 S_B_FROM_V%I(170)=5;S_B_FROM_V%J(170)=1;
 S_B_FROM_V%I(171)=5;S_B_FROM_V%J(171)=0;
 S_B_FROM_V%I(172)=4;S_B_FROM_V%J(172)=18;
 S_B_FROM_V%I(173)=4;S_B_FROM_V%J(173)=17;
 S_B_FROM_V%I(174)=4;S_B_FROM_V%J(174)=16;
 S_B_FROM_V%I(175)=4;S_B_FROM_V%J(175)=15;
 S_B_FROM_V%I(176)=4;S_B_FROM_V%J(176)=14;
 S_B_FROM_V%I(177)=4;S_B_FROM_V%J(177)=13;
 S_B_FROM_V%I(178)=4;S_B_FROM_V%J(178)=12;
 S_B_FROM_V%I(179)=4;S_B_FROM_V%J(179)=11;
 S_B_FROM_V%I(180)=4;S_B_FROM_V%J(180)=10;
 S_B_FROM_V%I(181)=4;S_B_FROM_V%J(181)=9;
 S_B_FROM_V%I(182)=4;S_B_FROM_V%J(182)=8;
 S_B_FROM_V%I(183)=4;S_B_FROM_V%J(183)=7;
 S_B_FROM_V%I(184)=4;S_B_FROM_V%J(184)=6;
 S_B_FROM_V%I(185)=4;S_B_FROM_V%J(185)=5;
 S_B_FROM_V%I(186)=4;S_B_FROM_V%J(186)=4;
 S_B_FROM_V%I(187)=4;S_B_FROM_V%J(187)=3;
 S_B_FROM_V%I(188)=4;S_B_FROM_V%J(188)=2;
 S_B_FROM_V%I(189)=4;S_B_FROM_V%J(189)=1;
 S_B_FROM_V%I(190)=4;S_B_FROM_V%J(190)=0;
 S_B_FROM_V%I(191)=3;S_B_FROM_V%J(191)=19;
 S_B_FROM_V%I(192)=3;S_B_FROM_V%J(192)=18;
 S_B_FROM_V%I(193)=3;S_B_FROM_V%J(193)=17;
 S_B_FROM_V%I(194)=3;S_B_FROM_V%J(194)=16;
 S_B_FROM_V%I(195)=3;S_B_FROM_V%J(195)=15;
 S_B_FROM_V%I(196)=3;S_B_FROM_V%J(196)=14;
 S_B_FROM_V%I(197)=3;S_B_FROM_V%J(197)=13;
 S_B_FROM_V%I(198)=3;S_B_FROM_V%J(198)=12;
 S_B_FROM_V%I(199)=3;S_B_FROM_V%J(199)=11;
 S_B_FROM_V%I(200)=3;S_B_FROM_V%J(200)=10;
 S_B_FROM_V%I(201)=3;S_B_FROM_V%J(201)=9;
 S_B_FROM_V%I(202)=3;S_B_FROM_V%J(202)=8;
 S_B_FROM_V%I(203)=3;S_B_FROM_V%J(203)=7;
 S_B_FROM_V%I(204)=3;S_B_FROM_V%J(204)=6;
 S_B_FROM_V%I(205)=3;S_B_FROM_V%J(205)=5;
 S_B_FROM_V%I(206)=3;S_B_FROM_V%J(206)=4;
 S_B_FROM_V%I(207)=3;S_B_FROM_V%J(207)=3;
 S_B_FROM_V%I(208)=3;S_B_FROM_V%J(208)=2;
 S_B_FROM_V%I(209)=3;S_B_FROM_V%J(209)=1;
 S_B_FROM_V%I(210)=3;S_B_FROM_V%J(210)=0;
 S_B_FROM_V%I(211)=2;S_B_FROM_V%J(211)=20;
 S_B_FROM_V%I(212)=2;S_B_FROM_V%J(212)=19;
 S_B_FROM_V%I(213)=2;S_B_FROM_V%J(213)=18;
 S_B_FROM_V%I(214)=2;S_B_FROM_V%J(214)=17;
 S_B_FROM_V%I(215)=2;S_B_FROM_V%J(215)=16;
 S_B_FROM_V%I(216)=2;S_B_FROM_V%J(216)=15;
 S_B_FROM_V%I(217)=2;S_B_FROM_V%J(217)=14;
 S_B_FROM_V%I(218)=2;S_B_FROM_V%J(218)=13;
 S_B_FROM_V%I(219)=2;S_B_FROM_V%J(219)=12;
 S_B_FROM_V%I(220)=2;S_B_FROM_V%J(220)=11;
 S_B_FROM_V%I(221)=2;S_B_FROM_V%J(221)=10;
 S_B_FROM_V%I(222)=2;S_B_FROM_V%J(222)=9;
 S_B_FROM_V%I(223)=2;S_B_FROM_V%J(223)=8;
 S_B_FROM_V%I(224)=2;S_B_FROM_V%J(224)=7;
 S_B_FROM_V%I(225)=2;S_B_FROM_V%J(225)=6;
 S_B_FROM_V%I(226)=2;S_B_FROM_V%J(226)=5;
 S_B_FROM_V%I(227)=2;S_B_FROM_V%J(227)=4;
 S_B_FROM_V%I(228)=2;S_B_FROM_V%J(228)=3;
 S_B_FROM_V%I(229)=2;S_B_FROM_V%J(229)=2;
 S_B_FROM_V%I(230)=2;S_B_FROM_V%J(230)=1;
 S_B_FROM_V%I(231)=2;S_B_FROM_V%J(231)=0;
 S_B_FROM_V%I(232)=1;S_B_FROM_V%J(232)=21;
 S_B_FROM_V%I(233)=1;S_B_FROM_V%J(233)=20;
 S_B_FROM_V%I(234)=1;S_B_FROM_V%J(234)=19;
 S_B_FROM_V%I(235)=1;S_B_FROM_V%J(235)=18;
 S_B_FROM_V%I(236)=1;S_B_FROM_V%J(236)=17;
 S_B_FROM_V%I(237)=1;S_B_FROM_V%J(237)=16;
 S_B_FROM_V%I(238)=1;S_B_FROM_V%J(238)=15;
 S_B_FROM_V%I(239)=1;S_B_FROM_V%J(239)=14;
 S_B_FROM_V%I(240)=1;S_B_FROM_V%J(240)=13;
 S_B_FROM_V%I(241)=1;S_B_FROM_V%J(241)=12;
 S_B_FROM_V%I(242)=1;S_B_FROM_V%J(242)=11;
 S_B_FROM_V%I(243)=1;S_B_FROM_V%J(243)=10;
 S_B_FROM_V%I(244)=1;S_B_FROM_V%J(244)=9;
 S_B_FROM_V%I(245)=1;S_B_FROM_V%J(245)=8;
 S_B_FROM_V%I(246)=1;S_B_FROM_V%J(246)=7;
 S_B_FROM_V%I(247)=1;S_B_FROM_V%J(247)=6;
 S_B_FROM_V%I(248)=1;S_B_FROM_V%J(248)=5;
 S_B_FROM_V%I(249)=1;S_B_FROM_V%J(249)=4;
 S_B_FROM_V%I(250)=1;S_B_FROM_V%J(250)=3;
 S_B_FROM_V%I(251)=1;S_B_FROM_V%J(251)=2;
 S_B_FROM_V%I(252)=1;S_B_FROM_V%J(252)=1;
 S_B_FROM_V%I(253)=1;S_B_FROM_V%J(253)=0;
 S_B_FROM_V%I(254)=0;S_B_FROM_V%J(254)=22;
 S_B_FROM_V%I(255)=0;S_B_FROM_V%J(255)=21;
 S_B_FROM_V%I(256)=0;S_B_FROM_V%J(256)=20;
 S_B_FROM_V%I(257)=0;S_B_FROM_V%J(257)=19;
 S_B_FROM_V%I(258)=0;S_B_FROM_V%J(258)=18;
 S_B_FROM_V%I(259)=0;S_B_FROM_V%J(259)=17;
 S_B_FROM_V%I(260)=0;S_B_FROM_V%J(260)=16;
 S_B_FROM_V%I(261)=0;S_B_FROM_V%J(261)=15;
 S_B_FROM_V%I(262)=0;S_B_FROM_V%J(262)=14;
 S_B_FROM_V%I(263)=0;S_B_FROM_V%J(263)=13;
 S_B_FROM_V%I(264)=0;S_B_FROM_V%J(264)=12;
 S_B_FROM_V%I(265)=0;S_B_FROM_V%J(265)=11;
 S_B_FROM_V%I(266)=0;S_B_FROM_V%J(266)=10;
 S_B_FROM_V%I(267)=0;S_B_FROM_V%J(267)=9;
 S_B_FROM_V%I(268)=0;S_B_FROM_V%J(268)=8;
 S_B_FROM_V%I(269)=0;S_B_FROM_V%J(269)=7;
 S_B_FROM_V%I(270)=0;S_B_FROM_V%J(270)=6;
 S_B_FROM_V%I(271)=0;S_B_FROM_V%J(271)=5;
 S_B_FROM_V%I(272)=0;S_B_FROM_V%J(272)=4;
 S_B_FROM_V%I(273)=0;S_B_FROM_V%J(273)=3;
 S_B_FROM_V%I(274)=0;S_B_FROM_V%J(274)=2;
 S_B_FROM_V%I(275)=0;S_B_FROM_V%J(275)=1;
 S_B_FROM_V%I(276)=0;S_B_FROM_V%J(276)=0;
 S_B_FROM_V%A_X(1,91)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,105)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,120)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,136)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,153)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,171)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,190)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,210)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,231)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,253)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(1,276)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(1,276)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,91)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,91)=-0.100000000000000E0_DP
 S_B_FROM_V%B_X(2,104)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,105)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,105)=0.111111111111111E0_DP
 S_B_FROM_V%B_X(2,119)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,120)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,120)=-0.125000000000000E0_DP
 S_B_FROM_V%B_X(2,135)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,136)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,136)=0.142857142857143E0_DP
 S_B_FROM_V%B_X(2,152)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,153)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,153)=-0.166666666666667E0_DP
 S_B_FROM_V%B_X(2,170)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,171)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,171)=0.200000000000000E0_DP
 S_B_FROM_V%B_X(2,189)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,190)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,190)=-0.250000000000000E0_DP
 S_B_FROM_V%B_X(2,209)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,210)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,210)=0.333333333333333E0_DP
 S_B_FROM_V%B_X(2,230)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,231)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(2,231)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(2,252)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(2,253)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(2,253)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(2,275)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(2,275)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,91)=0.511111111111111E0_DP
 S_B_FROM_V%B_Y(3,91)=0.100000000000000E0_DP
 S_B_FROM_V%B_X(3,104)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,104)=-0.222222222222221E0_DP
 S_B_FROM_V%A_X(3,105)=-0.513888888888889E0_DP
 S_B_FROM_V%B_Y(3,105)=-0.111111111111111E0_DP
 S_B_FROM_V%A_X(3,118)=-0.999999999999996E0_DP
 S_B_FROM_V%B_X(3,119)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,119)=0.250000000000000E0_DP
 S_B_FROM_V%A_X(3,120)=0.517857142857143E0_DP
 S_B_FROM_V%B_Y(3,120)=0.125000000000000E0_DP
 S_B_FROM_V%A_X(3,134)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,135)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,135)=-0.285714285714286E0_DP
 S_B_FROM_V%A_X(3,136)=-0.523809523809524E0_DP
 S_B_FROM_V%B_Y(3,136)=-0.142857142857143E0_DP
 S_B_FROM_V%A_X(3,151)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,152)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,152)=0.333333333333333E0_DP
 S_B_FROM_V%A_X(3,153)=0.533333333333333E0_DP
 S_B_FROM_V%B_Y(3,153)=0.166666666666667E0_DP
 S_B_FROM_V%A_X(3,169)=0.999999999999999E0_DP
 S_B_FROM_V%B_X(3,170)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,170)=-0.399999999999999E0_DP
 S_B_FROM_V%A_X(3,171)=-0.550000000000000E0_DP
 S_B_FROM_V%B_Y(3,171)=-0.200000000000000E0_DP
 S_B_FROM_V%A_X(3,188)=-0.999999999999999E0_DP
 S_B_FROM_V%B_X(3,189)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,189)=0.500000000000000E0_DP
 S_B_FROM_V%A_X(3,190)=0.583333333333333E0_DP
 S_B_FROM_V%B_Y(3,190)=0.250000000000000E0_DP
 S_B_FROM_V%A_X(3,208)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,209)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,209)=-0.666666666666667E0_DP
 S_B_FROM_V%A_X(3,210)=-0.666666666666667E0_DP
 S_B_FROM_V%B_Y(3,210)=-0.333333333333333E0_DP
 S_B_FROM_V%A_X(3,229)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,230)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,230)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,231)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,231)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(3,251)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(3,252)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(3,252)=-2.00000000000000E0_DP
 S_B_FROM_V%A_X(3,274)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(3,274)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,91)=-0.391666666666667E0_DP
 S_B_FROM_V%B_Y(4,91)=-0.154166666666667E0_DP
 S_B_FROM_V%B_X(4,104)=-1.54166666666667E0_DP
 S_B_FROM_V%A_Y(4,104)=0.333333333333333E0_DP
 S_B_FROM_V%A_X(4,105)=0.395833333333333E0_DP
 S_B_FROM_V%B_Y(4,105)=0.172619047619048E0_DP
 S_B_FROM_V%A_X(4,118)=1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,118)=0.375000000000000E0_DP
 S_B_FROM_V%B_X(4,119)=1.55357142857143E0_DP
 S_B_FROM_V%A_Y(4,119)=-0.375000000000000E0_DP
 S_B_FROM_V%A_X(4,120)=-0.401785714285714E0_DP
 S_B_FROM_V%B_Y(4,120)=-0.196428571428571E0_DP
 S_B_FROM_V%B_X(4,133)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,134)=-1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,134)=-0.428571428571428E0_DP
 S_B_FROM_V%B_X(4,135)=-1.57142857142857E0_DP
 S_B_FROM_V%A_Y(4,135)=0.428571428571429E0_DP
 S_B_FROM_V%A_X(4,136)=0.410714285714286E0_DP
 S_B_FROM_V%B_Y(4,136)=0.228571428571429E0_DP
 S_B_FROM_V%B_X(4,150)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,151)=1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,151)=0.500000000000000E0_DP
 S_B_FROM_V%B_X(4,152)=1.60000000000000E0_DP
 S_B_FROM_V%A_Y(4,152)=-0.500000000000000E0_DP
 S_B_FROM_V%A_X(4,153)=-0.425000000000000E0_DP
 S_B_FROM_V%B_Y(4,153)=-0.275000000000000E0_DP
 S_B_FROM_V%B_X(4,168)=0.999999999999999E0_DP
 S_B_FROM_V%A_X(4,169)=-1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,169)=-0.600000000000000E0_DP
 S_B_FROM_V%B_X(4,170)=-1.65000000000000E0_DP
 S_B_FROM_V%A_Y(4,170)=0.600000000000000E0_DP
 S_B_FROM_V%A_X(4,171)=0.450000000000000E0_DP
 S_B_FROM_V%B_Y(4,171)=0.350000000000000E0_DP
 S_B_FROM_V%B_X(4,187)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,188)=1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,188)=0.750000000000000E0_DP
 S_B_FROM_V%B_X(4,189)=1.75000000000000E0_DP
 S_B_FROM_V%A_Y(4,189)=-0.750000000000000E0_DP
 S_B_FROM_V%A_X(4,190)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(4,190)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(4,207)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,208)=-1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,208)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(4,209)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,209)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,210)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,210)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(4,228)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,229)=1.50000000000000E0_DP
 S_B_FROM_V%B_Y(4,229)=1.50000000000000E0_DP
 S_B_FROM_V%B_X(4,230)=3.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,230)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(4,250)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(4,251)=-3.00000000000000E0_DP
 S_B_FROM_V%B_Y(4,251)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(4,273)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(4,273)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,91)=0.410119047619048E0_DP
 S_B_FROM_V%B_Y(5,91)=0.158333333333333E0_DP
 S_B_FROM_V%B_X(5,104)=1.58333333333333E0_DP
 S_B_FROM_V%A_Y(5,104)=-0.690476190476191E0_DP
 S_B_FROM_V%A_X(5,105)=-0.419642857142857E0_DP
 S_B_FROM_V%B_Y(5,105)=-0.178571428571429E0_DP
 S_B_FROM_V%A_X(5,118)=-3.10714285714286E0_DP
 S_B_FROM_V%B_Y(5,118)=-0.750000000000000E0_DP
 S_B_FROM_V%B_X(5,119)=-1.60714285714286E0_DP
 S_B_FROM_V%A_Y(5,119)=0.785714285714285E0_DP
 S_B_FROM_V%A_X(5,120)=0.433928571428571E0_DP
 S_B_FROM_V%B_Y(5,120)=0.205357142857143E0_DP
 S_B_FROM_V%B_X(5,133)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,133)=0.571428571428577E0_DP
 S_B_FROM_V%A_X(5,134)=3.14285714285714E0_DP
 S_B_FROM_V%B_Y(5,134)=0.857142857142856E0_DP
 S_B_FROM_V%B_X(5,135)=1.64285714285714E0_DP
 S_B_FROM_V%A_Y(5,135)=-0.914285714285714E0_DP
 S_B_FROM_V%A_X(5,136)=-0.457142857142857E0_DP
 S_B_FROM_V%B_Y(5,136)=-0.242857142857143E0_DP
 S_B_FROM_V%A_X(5,149)=1.00000000000001E0_DP
 S_B_FROM_V%B_X(5,150)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,150)=-0.666666666666657E0_DP
 S_B_FROM_V%A_X(5,151)=-3.20000000000000E0_DP
 S_B_FROM_V%B_Y(5,151)=-0.999999999999999E0_DP
 S_B_FROM_V%B_X(5,152)=-1.70000000000000E0_DP
 S_B_FROM_V%A_Y(5,152)=1.10000000000000E0_DP
 S_B_FROM_V%A_X(5,153)=0.500000000000000E0_DP
 S_B_FROM_V%B_Y(5,153)=0.300000000000000E0_DP
 S_B_FROM_V%A_X(5,167)=-0.999999999999986E0_DP
 S_B_FROM_V%B_X(5,168)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,168)=0.800000000000001E0_DP
 S_B_FROM_V%A_X(5,169)=3.30000000000000E0_DP
 S_B_FROM_V%B_Y(5,169)=1.20000000000000E0_DP
 S_B_FROM_V%B_X(5,170)=1.80000000000000E0_DP
 S_B_FROM_V%A_Y(5,170)=-1.40000000000000E0_DP
 S_B_FROM_V%A_X(5,171)=-0.600000000000000E0_DP
 S_B_FROM_V%B_Y(5,171)=-0.400000000000000E0_DP
 S_B_FROM_V%A_X(5,186)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,187)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,187)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,188)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(5,188)=-1.50000000000000E0_DP
 S_B_FROM_V%B_X(5,189)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,189)=2.00000000000000E0_DP
 S_B_FROM_V%A_X(5,190)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,190)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(5,206)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,207)=-2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,207)=1.33333333333333E0_DP
 S_B_FROM_V%A_X(5,208)=4.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,208)=2.00000000000000E0_DP
 S_B_FROM_V%B_X(5,209)=4.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,209)=-4.00000000000000E0_DP
 S_B_FROM_V%A_X(5,227)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,228)=2.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,228)=-2.00000000000000E0_DP
 S_B_FROM_V%A_X(5,229)=-6.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,229)=-6.00000000000000E0_DP
 S_B_FROM_V%A_X(5,249)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(5,250)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(5,250)=4.00000000000000E0_DP
 S_B_FROM_V%A_X(5,272)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(5,272)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,91)=-0.358630952380952E0_DP
 S_B_FROM_V%B_Y(6,91)=-0.209821428571429E0_DP
 S_B_FROM_V%B_X(6,104)=-2.09821428571429E0_DP
 S_B_FROM_V%A_Y(6,104)=0.892857142857143E0_DP
 S_B_FROM_V%A_X(6,105)=0.372023809523809E0_DP
 S_B_FROM_V%B_Y(6,105)=0.241071428571429E0_DP
 S_B_FROM_V%A_X(6,118)=4.01785714285714E0_DP
 S_B_FROM_V%B_Y(6,118)=1.96428571428571E0_DP
 S_B_FROM_V%B_X(6,119)=2.16964285714286E0_DP
 S_B_FROM_V%A_Y(6,119)=-1.02678571428571E0_DP
 S_B_FROM_V%A_X(6,120)=-0.392857142857143E0_DP
 S_B_FROM_V%B_Y(6,120)=-0.285714285714286E0_DP
 S_B_FROM_V%B_X(6,133)=5.23809523809524E0_DP
 S_B_FROM_V%A_Y(6,133)=-1.42857142857143E0_DP
 S_B_FROM_V%A_X(6,134)=-4.10714285714285E0_DP
 S_B_FROM_V%B_Y(6,134)=-2.28571428571428E0_DP
 S_B_FROM_V%B_X(6,135)=-2.28571428571429E0_DP
 S_B_FROM_V%A_Y(6,135)=1.21428571428571E0_DP
 S_B_FROM_V%A_X(6,136)=0.428571428571428E0_DP
 S_B_FROM_V%B_Y(6,136)=0.357142857142857E0_DP
 S_B_FROM_V%A_X(6,149)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,149)=-0.833333333333330E0_DP
 S_B_FROM_V%B_X(6,150)=-5.33333333333333E0_DP
 S_B_FROM_V%A_Y(6,150)=1.66666666666666E0_DP
 S_B_FROM_V%A_X(6,151)=4.25000000000000E0_DP
 S_B_FROM_V%B_Y(6,151)=2.75000000000000E0_DP
 S_B_FROM_V%B_X(6,152)=2.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,152)=-1.50000000000000E0_DP
 S_B_FROM_V%A_X(6,153)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(6,153)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(6,166)=-0.999999999999996E0_DP
 S_B_FROM_V%A_X(6,167)=2.49999999999999E0_DP
 S_B_FROM_V%B_Y(6,167)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(6,168)=5.50000000000000E0_DP
 S_B_FROM_V%A_Y(6,168)=-2.00000000000000E0_DP
 S_B_FROM_V%A_X(6,169)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,169)=-3.50000000000000E0_DP
 S_B_FROM_V%B_X(6,170)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,170)=2.00000000000000E0_DP
 S_B_FROM_V%A_X(6,171)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,171)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(6,185)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,186)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,186)=-1.25000000000000E0_DP
 S_B_FROM_V%B_X(6,187)=-5.83333333333333E0_DP
 S_B_FROM_V%A_Y(6,187)=2.50000000000000E0_DP
 S_B_FROM_V%A_X(6,188)=5.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,188)=5.00000000000000E0_DP
 S_B_FROM_V%B_X(6,189)=5.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,189)=-5.00000000000000E0_DP
 S_B_FROM_V%B_X(6,205)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,206)=2.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,206)=1.66666666666666E0_DP
 S_B_FROM_V%B_X(6,207)=6.66666666666667E0_DP
 S_B_FROM_V%A_Y(6,207)=-3.33333333333333E0_DP
 S_B_FROM_V%A_X(6,208)=-10.0000000000000E0_DP
 S_B_FROM_V%B_Y(6,208)=-10.0000000000000E0_DP
 S_B_FROM_V%B_X(6,226)=0.999999999999998E0_DP
 S_B_FROM_V%A_X(6,227)=-2.50000000000000E0_DP
 S_B_FROM_V%B_Y(6,227)=-2.50000000000000E0_DP
 S_B_FROM_V%B_X(6,228)=-10.0000000000000E0_DP
 S_B_FROM_V%A_Y(6,228)=10.0000000000000E0_DP
 S_B_FROM_V%B_X(6,248)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(6,249)=5.00000000000000E0_DP
 S_B_FROM_V%B_Y(6,249)=5.00000000000000E0_DP
 S_B_FROM_V%B_X(6,271)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(6,271)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(7,91)=0.389880952380952E0_DP
 S_B_FROM_V%B_Y(7,91)=0.223214285714286E0_DP
 S_B_FROM_V%B_X(7,104)=2.23214285714286E0_DP
 S_B_FROM_V%A_Y(7,104)=-1.44642857142857E0_DP
 S_B_FROM_V%A_X(7,105)=-0.416666666666667E0_DP
 S_B_FROM_V%B_Y(7,105)=-0.261904761904762E0_DP
 S_B_FROM_V%A_X(7,118)=-6.50892857142857E0_DP
 S_B_FROM_V%B_Y(7,118)=-3.08035714285714E0_DP
 S_B_FROM_V%B_X(7,119)=-2.35714285714286E0_DP
 S_B_FROM_V%A_Y(7,119)=1.71428571428571E0_DP
 S_B_FROM_V%A_X(7,120)=0.464285714285714E0_DP
 S_B_FROM_V%B_Y(7,120)=0.321428571428571E0_DP
 S_B_FROM_V%B_X(7,133)=-8.21428571428572E0_DP
 S_B_FROM_V%A_Y(7,133)=4.57142857142857E0_DP
 S_B_FROM_V%A_X(7,134)=6.85714285714285E0_DP
 S_B_FROM_V%B_Y(7,134)=3.64285714285714E0_DP
 S_B_FROM_V%B_X(7,135)=2.57142857142857E0_DP
 S_B_FROM_V%A_Y(7,135)=-2.14285714285714E0_DP
 S_B_FROM_V%A_X(7,136)=-0.571428571428571E0_DP
 S_B_FROM_V%B_Y(7,136)=-0.428571428571429E0_DP
 S_B_FROM_V%A_X(7,149)=8.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,149)=2.50000000000000E0_DP
 S_B_FROM_V%B_X(7,150)=8.50000000000000E0_DP
 S_B_FROM_V%A_Y(7,150)=-5.49999999999999E0_DP
 S_B_FROM_V%A_X(7,151)=-7.50000000000000E0_DP
 S_B_FROM_V%B_Y(7,151)=-4.50000000000000E0_DP
 S_B_FROM_V%B_X(7,152)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,152)=3.00000000000000E0_DP
 S_B_FROM_V%A_X(7,153)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,153)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(7,166)=3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,166)=-1.20000000000000E0_DP
 S_B_FROM_V%A_X(7,167)=-8.24999999999999E0_DP
 S_B_FROM_V%B_Y(7,167)=-3.00000000000000E0_DP
 S_B_FROM_V%B_X(7,168)=-9.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,168)=7.00000000000000E0_DP
 S_B_FROM_V%A_X(7,169)=9.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,169)=6.00000000000000E0_DP
 S_B_FROM_V%B_X(7,170)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,170)=-6.00000000000000E0_DP
 S_B_FROM_V%A_X(7,184)=-0.999999999999996E0_DP
 S_B_FROM_V%B_X(7,185)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,185)=1.49999999999999E0_DP
 S_B_FROM_V%A_X(7,186)=8.75000000000000E0_DP
 S_B_FROM_V%B_Y(7,186)=3.75000000000000E0_DP
 S_B_FROM_V%B_X(7,187)=10.0000000000000E0_DP
 S_B_FROM_V%A_Y(7,187)=-10.0000000000000E0_DP
 S_B_FROM_V%A_X(7,188)=-15.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,188)=-15.0000000000000E0_DP
 S_B_FROM_V%A_X(7,204)=0.999999999999993E0_DP
 S_B_FROM_V%B_X(7,205)=3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,205)=-2.00000000000000E0_DP
 S_B_FROM_V%A_X(7,206)=-10.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,206)=-5.00000000000000E0_DP
 S_B_FROM_V%B_X(7,207)=-20.0000000000000E0_DP
 S_B_FROM_V%A_Y(7,207)=20.0000000000000E0_DP
 S_B_FROM_V%A_X(7,225)=-1.00000000000000E0_DP
 S_B_FROM_V%B_X(7,226)=-3.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,226)=3.00000000000000E0_DP
 S_B_FROM_V%A_X(7,227)=15.0000000000000E0_DP
 S_B_FROM_V%B_Y(7,227)=15.0000000000000E0_DP
 S_B_FROM_V%A_X(7,247)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(7,248)=6.00000000000000E0_DP
 S_B_FROM_V%A_Y(7,248)=-6.00000000000000E0_DP
 S_B_FROM_V%A_X(7,270)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(7,270)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,91)=-0.375000000000000E0_DP
 S_B_FROM_V%B_Y(8,91)=-0.291666666666667E0_DP
 S_B_FROM_V%B_X(8,104)=-2.91666666666667E0_DP
 S_B_FROM_V%A_Y(8,104)=1.83333333333333E0_DP
 S_B_FROM_V%A_X(8,105)=0.416666666666667E0_DP
 S_B_FROM_V%B_Y(8,105)=0.361111111111111E0_DP
 S_B_FROM_V%A_X(8,118)=8.25000000000000E0_DP
 S_B_FROM_V%B_Y(8,118)=6.00000000000000E0_DP
 S_B_FROM_V%B_X(8,119)=3.25000000000000E0_DP
 S_B_FROM_V%A_Y(8,119)=-2.25000000000000E0_DP
 S_B_FROM_V%A_X(8,120)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(8,120)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(8,133)=16.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,133)=-8.50000000000000E0_DP
 S_B_FROM_V%A_X(8,134)=-9.00000000000000E0_DP
 S_B_FROM_V%B_Y(8,134)=-7.50000000000000E0_DP
 S_B_FROM_V%B_X(8,135)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,135)=3.00000000000000E0_DP
 S_B_FROM_V%A_X(8,136)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(8,136)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,149)=-14.8750000000000E0_DP
 S_B_FROM_V%B_Y(8,149)=-9.62500000000000E0_DP
 S_B_FROM_V%B_X(8,150)=-17.5000000000000E0_DP
 S_B_FROM_V%A_Y(8,150)=10.5000000000000E0_DP
 S_B_FROM_V%A_X(8,151)=10.5000000000000E0_DP
 S_B_FROM_V%B_Y(8,151)=10.5000000000000E0_DP
 S_B_FROM_V%B_X(8,152)=7.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,152)=-7.00000000000000E0_DP
 S_B_FROM_V%B_X(8,166)=-11.5500000000000E0_DP
 S_B_FROM_V%A_Y(8,166)=4.20000000000001E0_DP
 S_B_FROM_V%A_X(8,167)=15.7500000000000E0_DP
 S_B_FROM_V%B_Y(8,167)=12.2500000000000E0_DP
 S_B_FROM_V%B_X(8,168)=21.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,168)=-14.0000000000000E0_DP
 S_B_FROM_V%A_X(8,169)=-21.0000000000000E0_DP
 S_B_FROM_V%B_Y(8,169)=-21.0000000000000E0_DP
 S_B_FROM_V%A_X(8,184)=3.50000000000001E0_DP
 S_B_FROM_V%B_Y(8,184)=1.75000000000000E0_DP
 S_B_FROM_V%B_X(8,185)=12.2500000000000E0_DP
 S_B_FROM_V%A_Y(8,185)=-5.25000000000000E0_DP
 S_B_FROM_V%A_X(8,186)=-17.5000000000000E0_DP
 S_B_FROM_V%B_Y(8,186)=-17.5000000000000E0_DP
 S_B_FROM_V%B_X(8,187)=-35.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,187)=35.0000000000000E0_DP
 S_B_FROM_V%B_X(8,203)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,204)=-3.50000000000000E0_DP
 S_B_FROM_V%B_Y(8,204)=-2.33333333333333E0_DP
 S_B_FROM_V%B_X(8,205)=-14.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,205)=7.00000000000000E0_DP
 S_B_FROM_V%A_X(8,206)=35.0000000000000E0_DP
 S_B_FROM_V%B_Y(8,206)=35.0000000000000E0_DP
 S_B_FROM_V%B_X(8,224)=-0.999999999999999E0_DP
 S_B_FROM_V%A_X(8,225)=3.50000000000000E0_DP
 S_B_FROM_V%B_Y(8,225)=3.50000000000000E0_DP
 S_B_FROM_V%B_X(8,226)=21.0000000000000E0_DP
 S_B_FROM_V%A_Y(8,226)=-21.0000000000000E0_DP
 S_B_FROM_V%B_X(8,246)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(8,247)=-7.00000000000000E0_DP
 S_B_FROM_V%B_Y(8,247)=-7.00000000000000E0_DP
 S_B_FROM_V%B_X(8,269)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(8,269)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(9,91)=0.444444444444444E0_DP
 S_B_FROM_V%B_Y(9,91)=0.333333333333333E0_DP
 S_B_FROM_V%B_X(9,104)=3.33333333333333E0_DP
 S_B_FROM_V%A_Y(9,104)=-2.88888888888889E0_DP
 S_B_FROM_V%A_X(9,105)=-0.555555555555556E0_DP
 S_B_FROM_V%B_Y(9,105)=-0.444444444444444E0_DP
 S_B_FROM_V%A_X(9,118)=-13.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,118)=-9.00000000000000E0_DP
 S_B_FROM_V%B_X(9,119)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,119)=4.00000000000000E0_DP
 S_B_FROM_V%A_X(9,120)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,120)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(9,133)=-24.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,133)=20.0000000000000E0_DP
 S_B_FROM_V%A_X(9,134)=16.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,134)=12.0000000000000E0_DP
 S_B_FROM_V%B_X(9,135)=8.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,135)=-8.00000000000000E0_DP
 S_B_FROM_V%A_X(9,149)=35.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,149)=21.0000000000000E0_DP
 S_B_FROM_V%B_X(9,150)=28.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,150)=-28.0000000000000E0_DP
 S_B_FROM_V%A_X(9,151)=-28.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,151)=-28.0000000000000E0_DP
 S_B_FROM_V%B_X(9,166)=25.2000000000000E0_DP
 S_B_FROM_V%A_Y(9,166)=-19.6000000000000E0_DP
 S_B_FROM_V%A_X(9,167)=-42.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,167)=-28.0000000000000E0_DP
 S_B_FROM_V%B_X(9,168)=-56.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,168)=56.0000000000000E0_DP
 S_B_FROM_V%A_X(9,184)=-16.3333333333333E0_DP
 S_B_FROM_V%B_Y(9,184)=-7.00000000000000E0_DP
 S_B_FROM_V%B_X(9,185)=-28.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,185)=28.0000000000000E0_DP
 S_B_FROM_V%A_X(9,186)=70.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,186)=70.0000000000000E0_DP
 S_B_FROM_V%B_X(9,203)=-4.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,203)=2.66666666666666E0_DP
 S_B_FROM_V%A_X(9,204)=18.6666666666667E0_DP
 S_B_FROM_V%B_Y(9,204)=9.33333333333333E0_DP
 S_B_FROM_V%B_X(9,205)=56.0000000000000E0_DP
 S_B_FROM_V%A_Y(9,205)=-56.0000000000000E0_DP
 S_B_FROM_V%A_X(9,223)=0.999999999999999E0_DP
 S_B_FROM_V%B_X(9,224)=4.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,224)=-4.00000000000000E0_DP
 S_B_FROM_V%A_X(9,225)=-28.0000000000000E0_DP
 S_B_FROM_V%B_Y(9,225)=-28.0000000000000E0_DP
 S_B_FROM_V%A_X(9,245)=-0.999999999999999E0_DP
 S_B_FROM_V%B_X(9,246)=-8.00000000000000E0_DP
 S_B_FROM_V%A_Y(9,246)=8.00000000000000E0_DP
 S_B_FROM_V%A_X(9,268)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(9,268)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(10,91)=-0.500000000000000E0_DP
 S_B_FROM_V%B_Y(10,91)=-0.500000000000000E0_DP
 S_B_FROM_V%B_X(10,104)=-5.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,104)=4.00000000000000E0_DP
 S_B_FROM_V%A_X(10,105)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(10,105)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(10,118)=18.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,118)=18.0000000000000E0_DP
 S_B_FROM_V%B_X(10,119)=9.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,119)=-9.00000000000000E0_DP
 S_B_FROM_V%B_X(10,133)=48.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,133)=-36.0000000000000E0_DP
 S_B_FROM_V%A_X(10,134)=-36.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,134)=-36.0000000000000E0_DP
 S_B_FROM_V%A_X(10,149)=-63.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,149)=-63.0000000000000E0_DP
 S_B_FROM_V%B_X(10,150)=-84.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,150)=84.0000000000000E0_DP
 S_B_FROM_V%B_X(10,166)=-75.6000000000000E0_DP
 S_B_FROM_V%A_Y(10,166)=50.4000000000000E0_DP
 S_B_FROM_V%A_X(10,167)=126.000000000000E0_DP
 S_B_FROM_V%B_Y(10,167)=126.000000000000E0_DP
 S_B_FROM_V%A_X(10,184)=42.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,184)=42.0000000000000E0_DP
 S_B_FROM_V%B_X(10,185)=126.000000000000E0_DP
 S_B_FROM_V%A_Y(10,185)=-126.000000000000E0_DP
 S_B_FROM_V%B_X(10,203)=24.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,203)=-12.0000000000000E0_DP
 S_B_FROM_V%A_X(10,204)=-84.0000000000000E0_DP
 S_B_FROM_V%B_Y(10,204)=-84.0000000000000E0_DP
 S_B_FROM_V%A_X(10,223)=-4.50000000000000E0_DP
 S_B_FROM_V%B_Y(10,223)=-4.50000000000000E0_DP
 S_B_FROM_V%B_X(10,224)=-36.0000000000000E0_DP
 S_B_FROM_V%A_Y(10,224)=36.0000000000000E0_DP
 S_B_FROM_V%B_X(10,244)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(10,245)=9.00000000000000E0_DP
 S_B_FROM_V%B_Y(10,245)=9.00000000000000E0_DP
 S_B_FROM_V%B_X(10,267)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(10,267)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(11,91)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(11,91)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(11,104)=10.0000000000000E0_DP
 S_B_FROM_V%A_Y(11,104)=-10.0000000000000E0_DP
 S_B_FROM_V%A_X(11,118)=-45.0000000000000E0_DP
 S_B_FROM_V%B_Y(11,118)=-45.0000000000000E0_DP
 S_B_FROM_V%B_X(11,133)=-120.000000000000E0_DP
 S_B_FROM_V%A_Y(11,133)=120.000000000000E0_DP
 S_B_FROM_V%A_X(11,149)=210.000000000000E0_DP
 S_B_FROM_V%B_Y(11,149)=210.000000000000E0_DP
 S_B_FROM_V%B_X(11,166)=252.000000000000E0_DP
 S_B_FROM_V%A_Y(11,166)=-252.000000000000E0_DP
 S_B_FROM_V%A_X(11,184)=-210.000000000000E0_DP
 S_B_FROM_V%B_Y(11,184)=-210.000000000000E0_DP
 S_B_FROM_V%B_X(11,203)=-120.000000000000E0_DP
 S_B_FROM_V%A_Y(11,203)=120.000000000000E0_DP
 S_B_FROM_V%A_X(11,223)=45.0000000000000E0_DP
 S_B_FROM_V%B_Y(11,223)=45.0000000000000E0_DP
 S_B_FROM_V%B_X(11,244)=10.0000000000000E0_DP
 S_B_FROM_V%A_Y(11,244)=-10.0000000000000E0_DP
 S_B_FROM_V%A_X(11,266)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(11,266)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(12,78)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(12,78)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(12,90)=11.0000000000000E0_DP
 S_B_FROM_V%A_Y(12,90)=-11.0000000000000E0_DP
 S_B_FROM_V%A_X(12,103)=-55.0000000000000E0_DP
 S_B_FROM_V%B_Y(12,103)=-55.0000000000000E0_DP
 S_B_FROM_V%B_X(12,117)=-165.000000000000E0_DP
 S_B_FROM_V%A_Y(12,117)=165.000000000000E0_DP
 S_B_FROM_V%A_X(12,132)=330.000000000000E0_DP
 S_B_FROM_V%B_Y(12,132)=330.000000000000E0_DP
 S_B_FROM_V%B_X(12,148)=462.000000000000E0_DP
 S_B_FROM_V%A_Y(12,148)=-462.000000000000E0_DP
 S_B_FROM_V%A_X(12,165)=-462.000000000000E0_DP
 S_B_FROM_V%B_Y(12,165)=-462.000000000000E0_DP
 S_B_FROM_V%B_X(12,183)=-330.000000000000E0_DP
 S_B_FROM_V%A_Y(12,183)=330.000000000000E0_DP
 S_B_FROM_V%A_X(12,202)=165.000000000000E0_DP
 S_B_FROM_V%B_Y(12,202)=165.000000000000E0_DP
 S_B_FROM_V%B_X(12,222)=55.0000000000000E0_DP
 S_B_FROM_V%A_Y(12,222)=-55.0000000000000E0_DP
 S_B_FROM_V%A_X(12,243)=-11.0000000000000E0_DP
 S_B_FROM_V%B_Y(12,243)=-11.0000000000000E0_DP
 S_B_FROM_V%B_X(12,265)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(12,265)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(13,66)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(13,66)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(13,77)=12.0000000000000E0_DP
 S_B_FROM_V%A_Y(13,77)=-12.0000000000000E0_DP
 S_B_FROM_V%A_X(13,89)=-66.0000000000000E0_DP
 S_B_FROM_V%B_Y(13,89)=-66.0000000000000E0_DP
 S_B_FROM_V%B_X(13,102)=-220.000000000000E0_DP
 S_B_FROM_V%A_Y(13,102)=220.000000000000E0_DP
 S_B_FROM_V%A_X(13,116)=495.000000000000E0_DP
 S_B_FROM_V%B_Y(13,116)=495.000000000000E0_DP
 S_B_FROM_V%B_X(13,131)=792.000000000000E0_DP
 S_B_FROM_V%A_Y(13,131)=-792.000000000000E0_DP
 S_B_FROM_V%A_X(13,147)=-924.000000000000E0_DP
 S_B_FROM_V%B_Y(13,147)=-924.000000000000E0_DP
 S_B_FROM_V%B_X(13,164)=-792.000000000000E0_DP
 S_B_FROM_V%A_Y(13,164)=792.000000000000E0_DP
 S_B_FROM_V%A_X(13,182)=495.000000000000E0_DP
 S_B_FROM_V%B_Y(13,182)=495.000000000000E0_DP
 S_B_FROM_V%B_X(13,201)=220.000000000000E0_DP
 S_B_FROM_V%A_Y(13,201)=-220.000000000000E0_DP
 S_B_FROM_V%A_X(13,221)=-66.0000000000000E0_DP
 S_B_FROM_V%B_Y(13,221)=-66.0000000000000E0_DP
 S_B_FROM_V%B_X(13,242)=-12.0000000000000E0_DP
 S_B_FROM_V%A_Y(13,242)=12.0000000000000E0_DP
 S_B_FROM_V%A_X(13,264)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(13,264)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(14,55)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(14,55)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(14,65)=13.0000000000000E0_DP
 S_B_FROM_V%A_Y(14,65)=-13.0000000000000E0_DP
 S_B_FROM_V%A_X(14,76)=-78.0000000000000E0_DP
 S_B_FROM_V%B_Y(14,76)=-78.0000000000000E0_DP
 S_B_FROM_V%B_X(14,88)=-286.000000000000E0_DP
 S_B_FROM_V%A_Y(14,88)=286.000000000000E0_DP
 S_B_FROM_V%A_X(14,101)=715.000000000000E0_DP
 S_B_FROM_V%B_Y(14,101)=715.000000000000E0_DP
 S_B_FROM_V%B_X(14,115)=1287.00000000000E0_DP
 S_B_FROM_V%A_Y(14,115)=-1287.00000000000E0_DP
 S_B_FROM_V%A_X(14,130)=-1716.00000000000E0_DP
 S_B_FROM_V%B_Y(14,130)=-1716.00000000000E0_DP
 S_B_FROM_V%B_X(14,146)=-1716.00000000000E0_DP
 S_B_FROM_V%A_Y(14,146)=1716.00000000000E0_DP
 S_B_FROM_V%A_X(14,163)=1287.00000000000E0_DP
 S_B_FROM_V%B_Y(14,163)=1287.00000000000E0_DP
 S_B_FROM_V%B_X(14,181)=715.000000000000E0_DP
 S_B_FROM_V%A_Y(14,181)=-715.000000000000E0_DP
 S_B_FROM_V%A_X(14,200)=-286.000000000000E0_DP
 S_B_FROM_V%B_Y(14,200)=-286.000000000000E0_DP
 S_B_FROM_V%B_X(14,220)=-78.0000000000000E0_DP
 S_B_FROM_V%A_Y(14,220)=78.0000000000000E0_DP
 S_B_FROM_V%A_X(14,241)=13.0000000000000E0_DP
 S_B_FROM_V%B_Y(14,241)=13.0000000000000E0_DP
 S_B_FROM_V%B_X(14,263)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(14,263)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(15,45)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(15,45)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(15,54)=14.0000000000000E0_DP
 S_B_FROM_V%A_Y(15,54)=-14.0000000000000E0_DP
 S_B_FROM_V%A_X(15,64)=-91.0000000000000E0_DP
 S_B_FROM_V%B_Y(15,64)=-91.0000000000000E0_DP
 S_B_FROM_V%B_X(15,75)=-364.000000000000E0_DP
 S_B_FROM_V%A_Y(15,75)=364.000000000000E0_DP
 S_B_FROM_V%A_X(15,87)=1001.00000000000E0_DP
 S_B_FROM_V%B_Y(15,87)=1001.00000000000E0_DP
 S_B_FROM_V%B_X(15,100)=2002.00000000000E0_DP
 S_B_FROM_V%A_Y(15,100)=-2002.00000000000E0_DP
 S_B_FROM_V%A_X(15,114)=-3003.00000000000E0_DP
 S_B_FROM_V%B_Y(15,114)=-3003.00000000000E0_DP
 S_B_FROM_V%B_X(15,129)=-3432.00000000000E0_DP
 S_B_FROM_V%A_Y(15,129)=3432.00000000000E0_DP
 S_B_FROM_V%A_X(15,145)=3003.00000000000E0_DP
 S_B_FROM_V%B_Y(15,145)=3003.00000000000E0_DP
 S_B_FROM_V%B_X(15,162)=2002.00000000000E0_DP
 S_B_FROM_V%A_Y(15,162)=-2002.00000000000E0_DP
 S_B_FROM_V%A_X(15,180)=-1001.00000000000E0_DP
 S_B_FROM_V%B_Y(15,180)=-1001.00000000000E0_DP
 S_B_FROM_V%B_X(15,199)=-364.000000000000E0_DP
 S_B_FROM_V%A_Y(15,199)=364.000000000000E0_DP
 S_B_FROM_V%A_X(15,219)=91.0000000000000E0_DP
 S_B_FROM_V%B_Y(15,219)=91.0000000000000E0_DP
 S_B_FROM_V%B_X(15,240)=14.0000000000000E0_DP
 S_B_FROM_V%A_Y(15,240)=-14.0000000000000E0_DP
 S_B_FROM_V%A_X(15,262)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(15,262)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(16,36)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(16,36)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(16,44)=15.0000000000000E0_DP
 S_B_FROM_V%A_Y(16,44)=-15.0000000000000E0_DP
 S_B_FROM_V%A_X(16,53)=-105.000000000000E0_DP
 S_B_FROM_V%B_Y(16,53)=-105.000000000000E0_DP
 S_B_FROM_V%B_X(16,63)=-455.000000000000E0_DP
 S_B_FROM_V%A_Y(16,63)=455.000000000000E0_DP
 S_B_FROM_V%A_X(16,74)=1365.00000000000E0_DP
 S_B_FROM_V%B_Y(16,74)=1365.00000000000E0_DP
 S_B_FROM_V%B_X(16,86)=3003.00000000000E0_DP
 S_B_FROM_V%A_Y(16,86)=-3003.00000000000E0_DP
 S_B_FROM_V%A_X(16,99)=-5005.00000000000E0_DP
 S_B_FROM_V%B_Y(16,99)=-5005.00000000000E0_DP
 S_B_FROM_V%B_X(16,113)=-6435.00000000000E0_DP
 S_B_FROM_V%A_Y(16,113)=6435.00000000000E0_DP
 S_B_FROM_V%A_X(16,128)=6435.00000000000E0_DP
 S_B_FROM_V%B_Y(16,128)=6435.00000000000E0_DP
 S_B_FROM_V%B_X(16,144)=5005.00000000000E0_DP
 S_B_FROM_V%A_Y(16,144)=-5005.00000000000E0_DP
 S_B_FROM_V%A_X(16,161)=-3003.00000000000E0_DP
 S_B_FROM_V%B_Y(16,161)=-3003.00000000000E0_DP
 S_B_FROM_V%B_X(16,179)=-1365.00000000000E0_DP
 S_B_FROM_V%A_Y(16,179)=1365.00000000000E0_DP
 S_B_FROM_V%A_X(16,198)=455.000000000000E0_DP
 S_B_FROM_V%B_Y(16,198)=455.000000000000E0_DP
 S_B_FROM_V%B_X(16,218)=105.000000000000E0_DP
 S_B_FROM_V%A_Y(16,218)=-105.000000000000E0_DP
 S_B_FROM_V%A_X(16,239)=-15.0000000000000E0_DP
 S_B_FROM_V%B_Y(16,239)=-15.0000000000000E0_DP
 S_B_FROM_V%B_X(16,261)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(16,261)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(17,28)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(17,28)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(17,35)=16.0000000000000E0_DP
 S_B_FROM_V%A_Y(17,35)=-16.0000000000000E0_DP
 S_B_FROM_V%A_X(17,43)=-120.000000000000E0_DP
 S_B_FROM_V%B_Y(17,43)=-120.000000000000E0_DP
 S_B_FROM_V%B_X(17,52)=-560.000000000000E0_DP
 S_B_FROM_V%A_Y(17,52)=560.000000000000E0_DP
 S_B_FROM_V%A_X(17,62)=1820.00000000000E0_DP
 S_B_FROM_V%B_Y(17,62)=1820.00000000000E0_DP
 S_B_FROM_V%B_X(17,73)=4368.00000000000E0_DP
 S_B_FROM_V%A_Y(17,73)=-4368.00000000000E0_DP
 S_B_FROM_V%A_X(17,85)=-8008.00000000000E0_DP
 S_B_FROM_V%B_Y(17,85)=-8008.00000000000E0_DP
 S_B_FROM_V%B_X(17,98)=-11440.0000000000E0_DP
 S_B_FROM_V%A_Y(17,98)=11440.0000000000E0_DP
 S_B_FROM_V%A_X(17,112)=12870.0000000000E0_DP
 S_B_FROM_V%B_Y(17,112)=12870.0000000000E0_DP
 S_B_FROM_V%B_X(17,127)=11440.0000000000E0_DP
 S_B_FROM_V%A_Y(17,127)=-11440.0000000000E0_DP
 S_B_FROM_V%A_X(17,143)=-8008.00000000000E0_DP
 S_B_FROM_V%B_Y(17,143)=-8008.00000000000E0_DP
 S_B_FROM_V%B_X(17,160)=-4368.00000000000E0_DP
 S_B_FROM_V%A_Y(17,160)=4368.00000000000E0_DP
 S_B_FROM_V%A_X(17,178)=1820.00000000000E0_DP
 S_B_FROM_V%B_Y(17,178)=1820.00000000000E0_DP
 S_B_FROM_V%B_X(17,197)=560.000000000000E0_DP
 S_B_FROM_V%A_Y(17,197)=-560.000000000000E0_DP
 S_B_FROM_V%A_X(17,217)=-120.000000000000E0_DP
 S_B_FROM_V%B_Y(17,217)=-120.000000000000E0_DP
 S_B_FROM_V%B_X(17,238)=-16.0000000000000E0_DP
 S_B_FROM_V%A_Y(17,238)=16.0000000000000E0_DP
 S_B_FROM_V%A_X(17,260)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(17,260)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(18,21)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(18,21)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(18,27)=17.0000000000000E0_DP
 S_B_FROM_V%A_Y(18,27)=-17.0000000000000E0_DP
 S_B_FROM_V%A_X(18,34)=-136.000000000000E0_DP
 S_B_FROM_V%B_Y(18,34)=-136.000000000000E0_DP
 S_B_FROM_V%B_X(18,42)=-680.000000000000E0_DP
 S_B_FROM_V%A_Y(18,42)=680.000000000000E0_DP
 S_B_FROM_V%A_X(18,51)=2380.00000000000E0_DP
 S_B_FROM_V%B_Y(18,51)=2380.00000000000E0_DP
 S_B_FROM_V%B_X(18,61)=6188.00000000000E0_DP
 S_B_FROM_V%A_Y(18,61)=-6188.00000000000E0_DP
 S_B_FROM_V%A_X(18,72)=-12376.0000000000E0_DP
 S_B_FROM_V%B_Y(18,72)=-12376.0000000000E0_DP
 S_B_FROM_V%B_X(18,84)=-19448.0000000000E0_DP
 S_B_FROM_V%A_Y(18,84)=19448.0000000000E0_DP
 S_B_FROM_V%A_X(18,97)=24310.0000000000E0_DP
 S_B_FROM_V%B_Y(18,97)=24310.0000000000E0_DP
 S_B_FROM_V%B_X(18,111)=24310.0000000000E0_DP
 S_B_FROM_V%A_Y(18,111)=-24310.0000000000E0_DP
 S_B_FROM_V%A_X(18,126)=-19448.0000000000E0_DP
 S_B_FROM_V%B_Y(18,126)=-19448.0000000000E0_DP
 S_B_FROM_V%B_X(18,142)=-12376.0000000000E0_DP
 S_B_FROM_V%A_Y(18,142)=12376.0000000000E0_DP
 S_B_FROM_V%A_X(18,159)=6188.00000000000E0_DP
 S_B_FROM_V%B_Y(18,159)=6188.00000000000E0_DP
 S_B_FROM_V%B_X(18,177)=2380.00000000000E0_DP
 S_B_FROM_V%A_Y(18,177)=-2380.00000000000E0_DP
 S_B_FROM_V%A_X(18,196)=-680.000000000000E0_DP
 S_B_FROM_V%B_Y(18,196)=-680.000000000000E0_DP
 S_B_FROM_V%B_X(18,216)=-136.000000000000E0_DP
 S_B_FROM_V%A_Y(18,216)=136.000000000000E0_DP
 S_B_FROM_V%A_X(18,237)=17.0000000000000E0_DP
 S_B_FROM_V%B_Y(18,237)=17.0000000000000E0_DP
 S_B_FROM_V%B_X(18,259)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(18,259)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(19,15)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(19,15)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(19,20)=18.0000000000000E0_DP
 S_B_FROM_V%A_Y(19,20)=-18.0000000000000E0_DP
 S_B_FROM_V%A_X(19,26)=-153.000000000000E0_DP
 S_B_FROM_V%B_Y(19,26)=-153.000000000000E0_DP
 S_B_FROM_V%B_X(19,33)=-816.000000000000E0_DP
 S_B_FROM_V%A_Y(19,33)=816.000000000000E0_DP
 S_B_FROM_V%A_X(19,41)=3060.00000000000E0_DP
 S_B_FROM_V%B_Y(19,41)=3060.00000000000E0_DP
 S_B_FROM_V%B_X(19,50)=8568.00000000000E0_DP
 S_B_FROM_V%A_Y(19,50)=-8568.00000000000E0_DP
 S_B_FROM_V%A_X(19,60)=-18564.0000000000E0_DP
 S_B_FROM_V%B_Y(19,60)=-18564.0000000000E0_DP
 S_B_FROM_V%B_X(19,71)=-31824.0000000000E0_DP
 S_B_FROM_V%A_Y(19,71)=31824.0000000000E0_DP
 S_B_FROM_V%A_X(19,83)=43758.0000000000E0_DP
 S_B_FROM_V%B_Y(19,83)=43758.0000000000E0_DP
 S_B_FROM_V%B_X(19,96)=48620.0000000000E0_DP
 S_B_FROM_V%A_Y(19,96)=-48620.0000000000E0_DP
 S_B_FROM_V%A_X(19,110)=-43758.0000000000E0_DP
 S_B_FROM_V%B_Y(19,110)=-43758.0000000000E0_DP
 S_B_FROM_V%B_X(19,125)=-31824.0000000000E0_DP
 S_B_FROM_V%A_Y(19,125)=31824.0000000000E0_DP
 S_B_FROM_V%A_X(19,141)=18564.0000000000E0_DP
 S_B_FROM_V%B_Y(19,141)=18564.0000000000E0_DP
 S_B_FROM_V%B_X(19,158)=8568.00000000000E0_DP
 S_B_FROM_V%A_Y(19,158)=-8568.00000000000E0_DP
 S_B_FROM_V%A_X(19,176)=-3060.00000000000E0_DP
 S_B_FROM_V%B_Y(19,176)=-3060.00000000000E0_DP
 S_B_FROM_V%B_X(19,195)=-816.000000000000E0_DP
 S_B_FROM_V%A_Y(19,195)=816.000000000000E0_DP
 S_B_FROM_V%A_X(19,215)=153.000000000000E0_DP
 S_B_FROM_V%B_Y(19,215)=153.000000000000E0_DP
 S_B_FROM_V%B_X(19,236)=18.0000000000000E0_DP
 S_B_FROM_V%A_Y(19,236)=-18.0000000000000E0_DP
 S_B_FROM_V%A_X(19,258)=-1.00000000000000E0_DP
 S_B_FROM_V%B_Y(19,258)=-1.00000000000000E0_DP
 S_B_FROM_V%A_X(20,10)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(20,10)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(20,14)=19.0000000000000E0_DP
 S_B_FROM_V%A_Y(20,14)=-19.0000000000000E0_DP
 S_B_FROM_V%A_X(20,19)=-171.000000000000E0_DP
 S_B_FROM_V%B_Y(20,19)=-171.000000000000E0_DP
 S_B_FROM_V%B_X(20,25)=-969.000000000000E0_DP
 S_B_FROM_V%A_Y(20,25)=969.000000000000E0_DP
 S_B_FROM_V%A_X(20,32)=3876.00000000000E0_DP
 S_B_FROM_V%B_Y(20,32)=3876.00000000000E0_DP
 S_B_FROM_V%B_X(20,40)=11628.0000000000E0_DP
 S_B_FROM_V%A_Y(20,40)=-11628.0000000000E0_DP
 S_B_FROM_V%A_X(20,49)=-27132.0000000000E0_DP
 S_B_FROM_V%B_Y(20,49)=-27132.0000000000E0_DP
 S_B_FROM_V%B_X(20,59)=-50388.0000000000E0_DP
 S_B_FROM_V%A_Y(20,59)=50388.0000000000E0_DP
 S_B_FROM_V%A_X(20,70)=75582.0000000000E0_DP
 S_B_FROM_V%B_Y(20,70)=75582.0000000000E0_DP
 S_B_FROM_V%B_X(20,82)=92378.0000000000E0_DP
 S_B_FROM_V%A_Y(20,82)=-92378.0000000000E0_DP
 S_B_FROM_V%A_X(20,95)=-92378.0000000000E0_DP
 S_B_FROM_V%B_Y(20,95)=-92378.0000000000E0_DP
 S_B_FROM_V%B_X(20,109)=-75582.0000000000E0_DP
 S_B_FROM_V%A_Y(20,109)=75582.0000000000E0_DP
 S_B_FROM_V%A_X(20,124)=50388.0000000000E0_DP
 S_B_FROM_V%B_Y(20,124)=50388.0000000000E0_DP
 S_B_FROM_V%B_X(20,140)=27132.0000000000E0_DP
 S_B_FROM_V%A_Y(20,140)=-27132.0000000000E0_DP
 S_B_FROM_V%A_X(20,157)=-11628.0000000000E0_DP
 S_B_FROM_V%B_Y(20,157)=-11628.0000000000E0_DP
 S_B_FROM_V%B_X(20,175)=-3876.00000000000E0_DP
 S_B_FROM_V%A_Y(20,175)=3876.00000000000E0_DP
 S_B_FROM_V%A_X(20,194)=969.000000000000E0_DP
 S_B_FROM_V%B_Y(20,194)=969.000000000000E0_DP
 S_B_FROM_V%B_X(20,214)=171.000000000000E0_DP
 S_B_FROM_V%A_Y(20,214)=-171.000000000000E0_DP
 S_B_FROM_V%A_X(20,235)=-19.0000000000000E0_DP
 S_B_FROM_V%B_Y(20,235)=-19.0000000000000E0_DP
 S_B_FROM_V%B_X(20,257)=-1.00000000000000E0_DP
 S_B_FROM_V%A_Y(20,257)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(21,6)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(21,6)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(21,9)=20.0000000000000E0_DP
 S_B_FROM_V%A_Y(21,9)=-20.0000000000000E0_DP
 S_B_FROM_V%A_X(21,13)=-190.000000000000E0_DP
 S_B_FROM_V%B_Y(21,13)=-190.000000000000E0_DP
 S_B_FROM_V%B_X(21,18)=-1140.00000000000E0_DP
 S_B_FROM_V%A_Y(21,18)=1140.00000000000E0_DP
 S_B_FROM_V%A_X(21,24)=4845.00000000000E0_DP
 S_B_FROM_V%B_Y(21,24)=4845.00000000000E0_DP
 S_B_FROM_V%B_X(21,31)=15504.0000000000E0_DP
 S_B_FROM_V%A_Y(21,31)=-15504.0000000000E0_DP
 S_B_FROM_V%A_X(21,39)=-38760.0000000000E0_DP
 S_B_FROM_V%B_Y(21,39)=-38760.0000000000E0_DP
 S_B_FROM_V%B_X(21,48)=-77520.0000000000E0_DP
 S_B_FROM_V%A_Y(21,48)=77520.0000000000E0_DP
 S_B_FROM_V%A_X(21,58)=125970.000000000E0_DP
 S_B_FROM_V%B_Y(21,58)=125970.000000000E0_DP
 S_B_FROM_V%B_X(21,69)=167960.000000000E0_DP
 S_B_FROM_V%A_Y(21,69)=-167960.000000000E0_DP
 S_B_FROM_V%A_X(21,81)=-184756.000000000E0_DP
 S_B_FROM_V%B_Y(21,81)=-184756.000000000E0_DP
 S_B_FROM_V%B_X(21,94)=-167960.000000000E0_DP
 S_B_FROM_V%A_Y(21,94)=167960.000000000E0_DP
 S_B_FROM_V%A_X(21,108)=125970.000000000E0_DP
 S_B_FROM_V%B_Y(21,108)=125970.000000000E0_DP
 S_B_FROM_V%B_X(21,123)=77520.0000000000E0_DP
 S_B_FROM_V%A_Y(21,123)=-77520.0000000000E0_DP
 S_B_FROM_V%A_X(21,139)=-38760.0000000000E0_DP
 S_B_FROM_V%B_Y(21,139)=-38760.0000000000E0_DP
 S_B_FROM_V%B_X(21,156)=-15504.0000000000E0_DP
 S_B_FROM_V%A_Y(21,156)=15504.0000000000E0_DP
 S_B_FROM_V%A_X(21,174)=4845.00000000000E0_DP
 S_B_FROM_V%B_Y(21,174)=4845.00000000000E0_DP
 S_B_FROM_V%B_X(21,193)=1140.00000000000E0_DP
 S_B_FROM_V%A_Y(21,193)=-1140.00000000000E0_DP
 S_B_FROM_V%A_X(21,213)=-190.000000000000E0_DP
 S_B_FROM_V%B_Y(21,213)=-190.000000000000E0_DP
 S_B_FROM_V%B_X(21,234)=-20.0000000000000E0_DP
 S_B_FROM_V%A_Y(21,234)=20.0000000000000E0_DP
 S_B_FROM_V%A_X(21,256)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(21,256)=1.00000000000000E0_DP
 S_B_FROM_V%A_X(22,3)=1.00000000000000E0_DP
 S_B_FROM_V%B_Y(22,3)=1.00000000000000E0_DP
 S_B_FROM_V%B_X(22,5)=21.0000000000000E0_DP
 S_B_FROM_V%A_Y(22,5)=-21.0000000000000E0_DP
 S_B_FROM_V%A_X(22,8)=-210.000000000000E0_DP
 S_B_FROM_V%B_Y(22,8)=-210.000000000000E0_DP
 S_B_FROM_V%B_X(22,12)=-1330.00000000000E0_DP
 S_B_FROM_V%A_Y(22,12)=1330.00000000000E0_DP
 S_B_FROM_V%A_X(22,17)=5985.00000000000E0_DP
 S_B_FROM_V%B_Y(22,17)=5985.00000000000E0_DP
 S_B_FROM_V%B_X(22,23)=20349.0000000000E0_DP
 S_B_FROM_V%A_Y(22,23)=-20349.0000000000E0_DP
 S_B_FROM_V%A_X(22,30)=-54264.0000000000E0_DP
 S_B_FROM_V%B_Y(22,30)=-54264.0000000000E0_DP
 S_B_FROM_V%B_X(22,38)=-116280.000000000E0_DP
 S_B_FROM_V%A_Y(22,38)=116280.000000000E0_DP
 S_B_FROM_V%A_X(22,47)=203490.000000000E0_DP
 S_B_FROM_V%B_Y(22,47)=203490.000000000E0_DP
 S_B_FROM_V%B_X(22,57)=293930.000000000E0_DP
 S_B_FROM_V%A_Y(22,57)=-293930.000000000E0_DP
 S_B_FROM_V%A_X(22,68)=-352716.000000000E0_DP
 S_B_FROM_V%B_Y(22,68)=-352716.000000000E0_DP
 S_B_FROM_V%B_X(22,80)=-352716.000000000E0_DP
 S_B_FROM_V%A_Y(22,80)=352716.000000000E0_DP
 S_B_FROM_V%A_X(22,93)=293930.000000000E0_DP
 S_B_FROM_V%B_Y(22,93)=293930.000000000E0_DP
 S_B_FROM_V%B_X(22,107)=203490.000000000E0_DP
 S_B_FROM_V%A_Y(22,107)=-203490.000000000E0_DP
 S_B_FROM_V%A_X(22,122)=-116280.000000000E0_DP
 S_B_FROM_V%B_Y(22,122)=-116280.000000000E0_DP
 S_B_FROM_V%B_X(22,138)=-54264.0000000000E0_DP
 S_B_FROM_V%A_Y(22,138)=54264.0000000000E0_DP
 S_B_FROM_V%A_X(22,155)=20349.0000000000E0_DP
 S_B_FROM_V%B_Y(22,155)=20349.0000000000E0_DP
 S_B_FROM_V%B_X(22,173)=5985.00000000000E0_DP
 S_B_FROM_V%A_Y(22,173)=-5985.00000000000E0_DP
 S_B_FROM_V%A_X(22,192)=-1330.00000000000E0_DP
 S_B_FROM_V%B_Y(22,192)=-1330.00000000000E0_DP
 S_B_FROM_V%B_X(22,212)=-210.000000000000E0_DP
 S_B_FROM_V%A_Y(22,212)=210.000000000000E0_DP
 S_B_FROM_V%A_X(22,233)=21.0000000000000E0_DP
 S_B_FROM_V%B_Y(22,233)=21.0000000000000E0_DP
 S_B_FROM_V%B_X(22,255)=1.00000000000000E0_DP
 S_B_FROM_V%A_Y(22,255)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(1,78)=-9.090909090909090E-002_DP
 S_B_FROM_V%VA(1,91)=0.100000000000000E0_DP
 S_B_FROM_V%VA(1,105)=-0.111111111111111E0_DP
 S_B_FROM_V%VA(1,120)=0.125000000000000E0_DP
 S_B_FROM_V%VA(1,136)=-0.142857142857143E0_DP
 S_B_FROM_V%VA(1,153)=0.166666666666667E0_DP
 S_B_FROM_V%VA(1,171)=-0.200000000000000E0_DP
 S_B_FROM_V%VA(1,190)=0.250000000000000E0_DP
 S_B_FROM_V%VA(1,210)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(1,231)=0.500000000000000E0_DP
 S_B_FROM_V%VA(1,253)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(1,275)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(2,78)=4.545454545454545E-002_DP
 S_B_FROM_V%VB(2,90)=0.100000000000000E0_DP
 S_B_FROM_V%VA(2,91)=-5.000000000000000E-002_DP
 S_B_FROM_V%VB(2,104)=-0.111111111111111E0_DP
 S_B_FROM_V%VA(2,105)=5.555555555555555E-002_DP
 S_B_FROM_V%VB(2,119)=0.125000000000000E0_DP
 S_B_FROM_V%VA(2,120)=-6.249999999999999E-002_DP
 S_B_FROM_V%VB(2,135)=-0.142857142857143E0_DP
 S_B_FROM_V%VA(2,136)=7.142857142857141E-002_DP
 S_B_FROM_V%VB(2,152)=0.166666666666667E0_DP
 S_B_FROM_V%VA(2,153)=-8.333333333333334E-002_DP
 S_B_FROM_V%VB(2,170)=-0.200000000000000E0_DP
 S_B_FROM_V%VA(2,171)=0.100000000000000E0_DP
 S_B_FROM_V%VB(2,189)=0.250000000000000E0_DP
 S_B_FROM_V%VA(2,190)=-0.125000000000000E0_DP
 S_B_FROM_V%VB(2,209)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(2,210)=0.166666666666667E0_DP
 S_B_FROM_V%VB(2,230)=0.500000000000000E0_DP
 S_B_FROM_V%VA(2,231)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(2,252)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(2,274)=0.500000000000000E0_DP
 S_B_FROM_V%VA(3,78)=-4.646464646464647E-002_DP
 S_B_FROM_V%VB(3,90)=-0.100000000000000E0_DP
 S_B_FROM_V%VA(3,91)=5.138888888888889E-002_DP
 S_B_FROM_V%VA(3,103)=0.111111111111111E0_DP
 S_B_FROM_V%VB(3,104)=0.111111111111111E0_DP
 S_B_FROM_V%VA(3,105)=-5.753968253968253E-002_DP
 S_B_FROM_V%VA(3,118)=-0.125000000000000E0_DP
 S_B_FROM_V%VB(3,119)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(3,120)=6.547619047619047E-002_DP
 S_B_FROM_V%VA(3,134)=0.142857142857143E0_DP
 S_B_FROM_V%VB(3,135)=0.142857142857143E0_DP
 S_B_FROM_V%VA(3,136)=-7.619047619047617E-002_DP
 S_B_FROM_V%VA(3,151)=-0.166666666666667E0_DP
 S_B_FROM_V%VB(3,152)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(3,153)=9.166666666666667E-002_DP
 S_B_FROM_V%VA(3,169)=0.200000000000000E0_DP
 S_B_FROM_V%VB(3,170)=0.200000000000000E0_DP
 S_B_FROM_V%VA(3,171)=-0.116666666666667E0_DP
 S_B_FROM_V%VA(3,188)=-0.250000000000000E0_DP
 S_B_FROM_V%VB(3,189)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(3,190)=0.166666666666667E0_DP
 S_B_FROM_V%VA(3,208)=0.333333333333333E0_DP
 S_B_FROM_V%VB(3,209)=0.333333333333333E0_DP
 S_B_FROM_V%VA(3,210)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(3,229)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(3,230)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(3,251)=1.00000000000000E0_DP
 S_B_FROM_V%VB(3,273)=0.333333333333333E0_DP
 S_B_FROM_V%VA(4,78)=3.560606060606061E-002_DP
 S_B_FROM_V%VB(4,90)=0.154166666666667E0_DP
 S_B_FROM_V%VA(4,91)=-3.958333333333335E-002_DP
 S_B_FROM_V%VA(4,103)=-0.166666666666667E0_DP
 S_B_FROM_V%VB(4,104)=-0.172619047619048E0_DP
 S_B_FROM_V%VA(4,105)=4.464285714285714E-002_DP
 S_B_FROM_V%VB(4,117)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(4,118)=0.187500000000000E0_DP
 S_B_FROM_V%VB(4,119)=0.196428571428571E0_DP
 S_B_FROM_V%VA(4,120)=-5.133928571428571E-002_DP
 S_B_FROM_V%VB(4,133)=0.142857142857143E0_DP
 S_B_FROM_V%VA(4,134)=-0.214285714285714E0_DP
 S_B_FROM_V%VB(4,135)=-0.228571428571429E0_DP
 S_B_FROM_V%VA(4,136)=6.071428571428570E-002_DP
 S_B_FROM_V%VB(4,150)=-0.166666666666667E0_DP
 S_B_FROM_V%VA(4,151)=0.250000000000000E0_DP
 S_B_FROM_V%VB(4,152)=0.275000000000000E0_DP
 S_B_FROM_V%VA(4,153)=-7.500000000000000E-002_DP
 S_B_FROM_V%VB(4,168)=0.200000000000000E0_DP
 S_B_FROM_V%VA(4,169)=-0.300000000000000E0_DP
 S_B_FROM_V%VB(4,170)=-0.350000000000000E0_DP
 S_B_FROM_V%VA(4,171)=0.100000000000000E0_DP
 S_B_FROM_V%VB(4,187)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(4,188)=0.375000000000000E0_DP
 S_B_FROM_V%VB(4,189)=0.500000000000000E0_DP
 S_B_FROM_V%VA(4,190)=-0.250000000000000E0_DP
 S_B_FROM_V%VB(4,207)=0.333333333333333E0_DP
 S_B_FROM_V%VA(4,208)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(4,209)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(4,228)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(4,229)=1.50000000000000E0_DP
 S_B_FROM_V%VB(4,250)=1.00000000000000E0_DP
 S_B_FROM_V%VA(4,272)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(5,78)=-3.728354978354979E-002_DP
 S_B_FROM_V%VB(5,90)=-0.158333333333333E0_DP
 S_B_FROM_V%VA(5,91)=4.196428571428571E-002_DP
 S_B_FROM_V%VA(5,103)=0.345238095238095E0_DP
 S_B_FROM_V%VB(5,104)=0.178571428571429E0_DP
 S_B_FROM_V%VA(5,105)=-4.821428571428571E-002_DP
 S_B_FROM_V%VB(5,117)=0.250000000000000E0_DP
 S_B_FROM_V%VA(5,118)=-0.392857142857143E0_DP
 S_B_FROM_V%VB(5,119)=-0.205357142857143E0_DP
 S_B_FROM_V%VA(5,120)=5.714285714285714E-002_DP
 S_B_FROM_V%VA(5,132)=-0.142857142857144E0_DP
 S_B_FROM_V%VB(5,133)=-0.285714285714285E0_DP
 S_B_FROM_V%VA(5,134)=0.457142857142857E0_DP
 S_B_FROM_V%VB(5,135)=0.242857142857143E0_DP
 S_B_FROM_V%VA(5,136)=-7.142857142857142E-002_DP
 S_B_FROM_V%VA(5,149)=0.166666666666664E0_DP
 S_B_FROM_V%VB(5,150)=0.333333333333333E0_DP
 S_B_FROM_V%VA(5,151)=-0.550000000000000E0_DP
 S_B_FROM_V%VB(5,152)=-0.300000000000000E0_DP
 S_B_FROM_V%VA(5,153)=0.100000000000000E0_DP
 S_B_FROM_V%VA(5,167)=-0.200000000000000E0_DP
 S_B_FROM_V%VB(5,168)=-0.400000000000000E0_DP
 S_B_FROM_V%VA(5,169)=0.700000000000000E0_DP
 S_B_FROM_V%VB(5,170)=0.400000000000000E0_DP
 S_B_FROM_V%VA(5,171)=-0.200000000000000E0_DP
 S_B_FROM_V%VA(5,186)=0.250000000000000E0_DP
 S_B_FROM_V%VB(5,187)=0.500000000000000E0_DP
 S_B_FROM_V%VA(5,188)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(5,189)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(5,206)=-0.333333333333333E0_DP
 S_B_FROM_V%VB(5,207)=-0.666666666666667E0_DP
 S_B_FROM_V%VA(5,208)=2.00000000000000E0_DP
 S_B_FROM_V%VA(5,227)=0.500000000000000E0_DP
 S_B_FROM_V%VB(5,228)=2.00000000000000E0_DP
 S_B_FROM_V%VA(5,249)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(5,271)=-0.200000000000000E0_DP
 S_B_FROM_V%VA(6,78)=3.260281385281385E-002_DP
 S_B_FROM_V%VB(6,90)=0.209821428571429E0_DP
 S_B_FROM_V%VA(6,91)=-3.720238095238095E-002_DP
 S_B_FROM_V%VA(6,103)=-0.446428571428571E0_DP
 S_B_FROM_V%VB(6,104)=-0.241071428571429E0_DP
 S_B_FROM_V%VA(6,105)=4.365079365079364E-002_DP
 S_B_FROM_V%VB(6,117)=-0.654761904761905E0_DP
 S_B_FROM_V%VA(6,118)=0.513392857142857E0_DP
 S_B_FROM_V%VB(6,119)=0.285714285714286E0_DP
 S_B_FROM_V%VA(6,120)=-5.357142857142856E-002_DP
 S_B_FROM_V%VA(6,132)=0.357142857142858E0_DP
 S_B_FROM_V%VB(6,133)=0.761904761904761E0_DP
 S_B_FROM_V%VA(6,134)=-0.607142857142857E0_DP
 S_B_FROM_V%VB(6,135)=-0.357142857142857E0_DP
 S_B_FROM_V%VA(6,136)=7.142857142857141E-002_DP
 S_B_FROM_V%VB(6,148)=0.166666666666666E0_DP
 S_B_FROM_V%VA(6,149)=-0.416666666666664E0_DP
 S_B_FROM_V%VB(6,150)=-0.916666666666666E0_DP
 S_B_FROM_V%VA(6,151)=0.750000000000000E0_DP
 S_B_FROM_V%VB(6,152)=0.500000000000000E0_DP
 S_B_FROM_V%VA(6,153)=-0.166666666666667E0_DP
 S_B_FROM_V%VB(6,166)=-0.200000000000000E0_DP
 S_B_FROM_V%VA(6,167)=0.500000000000000E0_DP
 S_B_FROM_V%VB(6,168)=1.16666666666667E0_DP
 S_B_FROM_V%VA(6,169)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(6,170)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(6,185)=0.250000000000000E0_DP
 S_B_FROM_V%VA(6,186)=-0.624999999999999E0_DP
 S_B_FROM_V%VB(6,187)=-1.66666666666667E0_DP
 S_B_FROM_V%VA(6,188)=2.50000000000000E0_DP
 S_B_FROM_V%VB(6,205)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(6,206)=0.833333333333333E0_DP
 S_B_FROM_V%VB(6,207)=3.33333333333333E0_DP
 S_B_FROM_V%VB(6,226)=0.500000000000000E0_DP
 S_B_FROM_V%VA(6,227)=-2.50000000000000E0_DP
 S_B_FROM_V%VB(6,248)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(6,270)=0.166666666666667E0_DP
 S_B_FROM_V%VA(7,78)=-3.544372294372294E-002_DP
 S_B_FROM_V%VB(7,90)=-0.223214285714286E0_DP
 S_B_FROM_V%VA(7,91)=4.166666666666666E-002_DP
 S_B_FROM_V%VA(7,103)=0.723214285714286E0_DP
 S_B_FROM_V%VB(7,104)=0.261904761904762E0_DP
 S_B_FROM_V%VA(7,105)=-5.158730158730158E-002_DP
 S_B_FROM_V%VB(7,117)=1.02678571428571E0_DP
 S_B_FROM_V%VA(7,118)=-0.857142857142857E0_DP
 S_B_FROM_V%VB(7,119)=-0.321428571428571E0_DP
 S_B_FROM_V%VA(7,120)=7.142857142857142E-002_DP
 S_B_FROM_V%VA(7,132)=-1.14285714285714E0_DP
 S_B_FROM_V%VB(7,133)=-1.21428571428571E0_DP
 S_B_FROM_V%VA(7,134)=1.07142857142857E0_DP
 S_B_FROM_V%VB(7,135)=0.428571428571429E0_DP
 S_B_FROM_V%VA(7,136)=-0.142857142857143E0_DP
 S_B_FROM_V%VB(7,148)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(7,149)=1.37500000000000E0_DP
 S_B_FROM_V%VB(7,150)=1.50000000000000E0_DP
 S_B_FROM_V%VA(7,151)=-1.50000000000000E0_DP
 S_B_FROM_V%VB(7,152)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(7,165)=0.199999999999999E0_DP
 S_B_FROM_V%VB(7,166)=0.600000000000000E0_DP
 S_B_FROM_V%VA(7,167)=-1.75000000000000E0_DP
 S_B_FROM_V%VB(7,168)=-2.00000000000000E0_DP
 S_B_FROM_V%VA(7,169)=3.00000000000000E0_DP
 S_B_FROM_V%VA(7,184)=-0.249999999999998E0_DP
 S_B_FROM_V%VB(7,185)=-0.750000000000000E0_DP
 S_B_FROM_V%VA(7,186)=2.50000000000000E0_DP
 S_B_FROM_V%VB(7,187)=5.00000000000000E0_DP
 S_B_FROM_V%VA(7,204)=0.333333333333333E0_DP
 S_B_FROM_V%VB(7,205)=1.00000000000000E0_DP
 S_B_FROM_V%VA(7,206)=-5.00000000000000E0_DP
 S_B_FROM_V%VA(7,225)=-0.500000000000000E0_DP
 S_B_FROM_V%VB(7,226)=-3.00000000000000E0_DP
 S_B_FROM_V%VA(7,247)=1.00000000000000E0_DP
 S_B_FROM_V%VB(7,269)=0.142857142857143E0_DP
 S_B_FROM_V%VA(8,78)=3.409090909090909E-002_DP
 S_B_FROM_V%VB(8,90)=0.291666666666667E0_DP
 S_B_FROM_V%VA(8,91)=-4.166666666666666E-002_DP
 S_B_FROM_V%VA(8,103)=-0.916666666666667E0_DP
 S_B_FROM_V%VB(8,104)=-0.361111111111111E0_DP
 S_B_FROM_V%VA(8,105)=5.555555555555555E-002_DP
 S_B_FROM_V%VB(8,117)=-2.00000000000000E0_DP
 S_B_FROM_V%VA(8,118)=1.12500000000000E0_DP
 S_B_FROM_V%VB(8,119)=0.500000000000000E0_DP
 S_B_FROM_V%VA(8,120)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(8,132)=2.12500000000000E0_DP
 S_B_FROM_V%VB(8,133)=2.50000000000000E0_DP
 S_B_FROM_V%VA(8,134)=-1.50000000000000E0_DP
 S_B_FROM_V%VB(8,135)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(8,148)=1.92500000000000E0_DP
 S_B_FROM_V%VA(8,149)=-2.62500000000000E0_DP
 S_B_FROM_V%VB(8,150)=-3.50000000000000E0_DP
 S_B_FROM_V%VA(8,151)=3.50000000000000E0_DP
 S_B_FROM_V%VA(8,165)=-0.700000000000001E0_DP
 S_B_FROM_V%VB(8,166)=-2.45000000000000E0_DP
 S_B_FROM_V%VA(8,167)=3.50000000000000E0_DP
 S_B_FROM_V%VB(8,168)=7.00000000000000E0_DP
 S_B_FROM_V%VB(8,183)=-0.250000000000000E0_DP
 S_B_FROM_V%VA(8,184)=0.875000000000000E0_DP
 S_B_FROM_V%VB(8,185)=3.50000000000000E0_DP
 S_B_FROM_V%VA(8,186)=-8.75000000000000E0_DP
 S_B_FROM_V%VB(8,203)=0.333333333333333E0_DP
 S_B_FROM_V%VA(8,204)=-1.16666666666667E0_DP
 S_B_FROM_V%VB(8,205)=-7.00000000000000E0_DP
 S_B_FROM_V%VB(8,224)=-0.500000000000000E0_DP
 S_B_FROM_V%VA(8,225)=3.50000000000000E0_DP
 S_B_FROM_V%VB(8,246)=1.00000000000000E0_DP
 S_B_FROM_V%VA(8,268)=-0.125000000000000E0_DP
 S_B_FROM_V%VA(9,78)=-4.040404040404040E-002_DP
 S_B_FROM_V%VB(9,90)=-0.333333333333333E0_DP
 S_B_FROM_V%VA(9,91)=5.555555555555555E-002_DP
 S_B_FROM_V%VA(9,103)=1.44444444444444E0_DP
 S_B_FROM_V%VB(9,104)=0.444444444444444E0_DP
 S_B_FROM_V%VA(9,105)=-0.111111111111111E0_DP
 S_B_FROM_V%VB(9,117)=3.00000000000000E0_DP
 S_B_FROM_V%VA(9,118)=-2.00000000000000E0_DP
 S_B_FROM_V%VB(9,119)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(9,132)=-5.00000000000000E0_DP
 S_B_FROM_V%VB(9,133)=-4.00000000000000E0_DP
 S_B_FROM_V%VA(9,134)=4.00000000000000E0_DP
 S_B_FROM_V%VB(9,148)=-4.20000000000000E0_DP
 S_B_FROM_V%VA(9,149)=7.00000000000000E0_DP
 S_B_FROM_V%VB(9,150)=9.33333333333333E0_DP
 S_B_FROM_V%VA(9,165)=3.26666666666667E0_DP
 S_B_FROM_V%VB(9,166)=5.60000000000000E0_DP
 S_B_FROM_V%VA(9,167)=-14.0000000000000E0_DP
 S_B_FROM_V%VB(9,183)=1.00000000000000E0_DP
 S_B_FROM_V%VA(9,184)=-4.66666666666666E0_DP
 S_B_FROM_V%VB(9,185)=-14.0000000000000E0_DP
 S_B_FROM_V%VA(9,202)=-0.333333333333333E0_DP
 S_B_FROM_V%VB(9,203)=-1.33333333333333E0_DP
 S_B_FROM_V%VA(9,204)=9.33333333333333E0_DP
 S_B_FROM_V%VA(9,223)=0.500000000000000E0_DP
 S_B_FROM_V%VB(9,224)=4.00000000000000E0_DP
 S_B_FROM_V%VA(9,245)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(9,267)=-0.111111111111111E0_DP
 S_B_FROM_V%VA(10,78)=4.545454545454546E-002_DP
 S_B_FROM_V%VB(10,90)=0.500000000000000E0_DP
 S_B_FROM_V%VA(10,91)=-0.100000000000000E0_DP
 S_B_FROM_V%VA(10,103)=-2.00000000000000E0_DP
 S_B_FROM_V%VB(10,104)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(10,117)=-6.00000000000000E0_DP
 S_B_FROM_V%VA(10,118)=4.50000000000000E0_DP
 S_B_FROM_V%VA(10,132)=9.00000000000000E0_DP
 S_B_FROM_V%VB(10,133)=12.0000000000000E0_DP
 S_B_FROM_V%VB(10,148)=12.6000000000000E0_DP
 S_B_FROM_V%VA(10,149)=-21.0000000000000E0_DP
 S_B_FROM_V%VA(10,165)=-8.40000000000000E0_DP
 S_B_FROM_V%VB(10,166)=-25.2000000000000E0_DP
 S_B_FROM_V%VB(10,183)=-6.00000000000000E0_DP
 S_B_FROM_V%VA(10,184)=21.0000000000000E0_DP
 S_B_FROM_V%VA(10,202)=1.50000000000000E0_DP
 S_B_FROM_V%VB(10,203)=12.0000000000000E0_DP
 S_B_FROM_V%VB(10,222)=0.500000000000000E0_DP
 S_B_FROM_V%VA(10,223)=-4.50000000000000E0_DP
 S_B_FROM_V%VB(10,244)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(10,266)=0.100000000000000E0_DP
 S_B_FROM_V%VA(11,78)=-9.090909090909091E-002_DP
 S_B_FROM_V%VB(11,90)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(11,103)=5.00000000000000E0_DP
 S_B_FROM_V%VB(11,117)=15.0000000000000E0_DP
 S_B_FROM_V%VA(11,132)=-30.0000000000000E0_DP
 S_B_FROM_V%VB(11,148)=-42.0000000000000E0_DP
 S_B_FROM_V%VA(11,165)=42.0000000000000E0_DP
 S_B_FROM_V%VB(11,183)=30.0000000000000E0_DP
 S_B_FROM_V%VA(11,202)=-15.0000000000000E0_DP
 S_B_FROM_V%VB(11,222)=-5.00000000000000E0_DP
 S_B_FROM_V%VA(11,243)=1.00000000000000E0_DP
 S_B_FROM_V%VB(11,265)=9.090909090909091E-002_DP
 S_B_FROM_V%VA(12,66)=-8.333333333333333E-002_DP
 S_B_FROM_V%VB(12,77)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(12,89)=5.50000000000000E0_DP
 S_B_FROM_V%VB(12,102)=18.3333333333333E0_DP
 S_B_FROM_V%VA(12,116)=-41.2500000000000E0_DP
 S_B_FROM_V%VB(12,131)=-66.0000000000000E0_DP
 S_B_FROM_V%VA(12,147)=77.0000000000000E0_DP
 S_B_FROM_V%VB(12,164)=66.0000000000000E0_DP
 S_B_FROM_V%VA(12,182)=-41.2500000000000E0_DP
 S_B_FROM_V%VB(12,201)=-18.3333333333333E0_DP
 S_B_FROM_V%VA(12,221)=5.50000000000000E0_DP
 S_B_FROM_V%VB(12,242)=1.00000000000000E0_DP
 S_B_FROM_V%VA(12,264)=-8.333333333333333E-002_DP
 S_B_FROM_V%VA(13,55)=-7.692307692307693E-002_DP
 S_B_FROM_V%VB(13,65)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(13,76)=6.00000000000000E0_DP
 S_B_FROM_V%VB(13,88)=22.0000000000000E0_DP
 S_B_FROM_V%VA(13,101)=-55.0000000000000E0_DP
 S_B_FROM_V%VB(13,115)=-99.0000000000000E0_DP
 S_B_FROM_V%VA(13,130)=132.000000000000E0_DP
 S_B_FROM_V%VB(13,146)=132.000000000000E0_DP
 S_B_FROM_V%VA(13,163)=-99.0000000000000E0_DP
 S_B_FROM_V%VB(13,181)=-55.0000000000000E0_DP
 S_B_FROM_V%VA(13,200)=22.0000000000000E0_DP
 S_B_FROM_V%VB(13,220)=6.00000000000000E0_DP
 S_B_FROM_V%VA(13,241)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(13,263)=-7.692307692307693E-002_DP
 S_B_FROM_V%VA(14,45)=-7.142857142857142E-002_DP
 S_B_FROM_V%VB(14,54)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(14,64)=6.50000000000000E0_DP
 S_B_FROM_V%VB(14,75)=26.0000000000000E0_DP
 S_B_FROM_V%VA(14,87)=-71.5000000000000E0_DP
 S_B_FROM_V%VB(14,100)=-143.000000000000E0_DP
 S_B_FROM_V%VA(14,114)=214.500000000000E0_DP
 S_B_FROM_V%VB(14,129)=245.142857142857E0_DP
 S_B_FROM_V%VA(14,145)=-214.500000000000E0_DP
 S_B_FROM_V%VB(14,162)=-143.000000000000E0_DP
 S_B_FROM_V%VA(14,180)=71.5000000000000E0_DP
 S_B_FROM_V%VB(14,199)=26.0000000000000E0_DP
 S_B_FROM_V%VA(14,219)=-6.50000000000000E0_DP
 S_B_FROM_V%VB(14,240)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(14,262)=7.142857142857142E-002_DP
 S_B_FROM_V%VA(15,36)=-6.666666666666667E-002_DP
 S_B_FROM_V%VB(15,44)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(15,53)=7.00000000000000E0_DP
 S_B_FROM_V%VB(15,63)=30.3333333333333E0_DP
 S_B_FROM_V%VA(15,74)=-91.0000000000000E0_DP
 S_B_FROM_V%VB(15,86)=-200.200000000000E0_DP
 S_B_FROM_V%VA(15,99)=333.666666666667E0_DP
 S_B_FROM_V%VB(15,113)=429.000000000000E0_DP
 S_B_FROM_V%VA(15,128)=-429.000000000000E0_DP
 S_B_FROM_V%VB(15,144)=-333.666666666667E0_DP
 S_B_FROM_V%VA(15,161)=200.200000000000E0_DP
 S_B_FROM_V%VB(15,179)=91.0000000000000E0_DP
 S_B_FROM_V%VA(15,198)=-30.3333333333333E0_DP
 S_B_FROM_V%VB(15,218)=-7.00000000000000E0_DP
 S_B_FROM_V%VA(15,239)=1.00000000000000E0_DP
 S_B_FROM_V%VB(15,261)=6.666666666666667E-002_DP
 S_B_FROM_V%VA(16,28)=-6.250000000000000E-002_DP
 S_B_FROM_V%VB(16,35)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(16,43)=7.50000000000000E0_DP
 S_B_FROM_V%VB(16,52)=35.0000000000000E0_DP
 S_B_FROM_V%VA(16,62)=-113.750000000000E0_DP
 S_B_FROM_V%VB(16,73)=-273.000000000000E0_DP
 S_B_FROM_V%VA(16,85)=500.500000000000E0_DP
 S_B_FROM_V%VB(16,98)=715.000000000000E0_DP
 S_B_FROM_V%VA(16,112)=-804.375000000000E0_DP
 S_B_FROM_V%VB(16,127)=-715.000000000000E0_DP
 S_B_FROM_V%VA(16,143)=500.500000000000E0_DP
 S_B_FROM_V%VB(16,160)=273.000000000000E0_DP
 S_B_FROM_V%VA(16,178)=-113.750000000000E0_DP
 S_B_FROM_V%VB(16,197)=-35.0000000000000E0_DP
 S_B_FROM_V%VA(16,217)=7.50000000000000E0_DP
 S_B_FROM_V%VB(16,238)=1.00000000000000E0_DP
 S_B_FROM_V%VA(16,260)=-6.250000000000000E-002_DP
 S_B_FROM_V%VA(17,21)=-5.882352941176471E-002_DP
 S_B_FROM_V%VB(17,27)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(17,34)=8.00000000000000E0_DP
 S_B_FROM_V%VB(17,42)=40.0000000000000E0_DP
 S_B_FROM_V%VA(17,51)=-140.000000000000E0_DP
 S_B_FROM_V%VB(17,61)=-364.000000000000E0_DP
 S_B_FROM_V%VA(17,72)=728.000000000000E0_DP
 S_B_FROM_V%VB(17,84)=1144.00000000000E0_DP
 S_B_FROM_V%VA(17,97)=-1430.00000000000E0_DP
 S_B_FROM_V%VB(17,111)=-1430.00000000000E0_DP
 S_B_FROM_V%VA(17,126)=1144.00000000000E0_DP
 S_B_FROM_V%VB(17,142)=728.000000000000E0_DP
 S_B_FROM_V%VA(17,159)=-364.000000000000E0_DP
 S_B_FROM_V%VB(17,177)=-140.000000000000E0_DP
 S_B_FROM_V%VA(17,196)=40.0000000000000E0_DP
 S_B_FROM_V%VB(17,216)=8.00000000000000E0_DP
 S_B_FROM_V%VA(17,237)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(17,259)=-5.882352941176471E-002_DP
 S_B_FROM_V%VA(18,15)=-5.555555555555555E-002_DP
 S_B_FROM_V%VB(18,20)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(18,26)=8.50000000000000E0_DP
 S_B_FROM_V%VB(18,33)=45.3333333333333E0_DP
 S_B_FROM_V%VA(18,41)=-170.000000000000E0_DP
 S_B_FROM_V%VB(18,50)=-476.000000000000E0_DP
 S_B_FROM_V%VA(18,60)=1031.33333333333E0_DP
 S_B_FROM_V%VB(18,71)=1768.00000000000E0_DP
 S_B_FROM_V%VA(18,83)=-2431.00000000000E0_DP
 S_B_FROM_V%VB(18,96)=-2701.11111111111E0_DP
 S_B_FROM_V%VA(18,110)=2431.00000000000E0_DP
 S_B_FROM_V%VB(18,125)=1768.00000000000E0_DP
 S_B_FROM_V%VA(18,141)=-1031.33333333333E0_DP
 S_B_FROM_V%VB(18,158)=-476.000000000000E0_DP
 S_B_FROM_V%VA(18,176)=170.000000000000E0_DP
 S_B_FROM_V%VB(18,195)=45.3333333333333E0_DP
 S_B_FROM_V%VA(18,215)=-8.50000000000000E0_DP
 S_B_FROM_V%VB(18,236)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(18,258)=5.555555555555555E-002_DP
 S_B_FROM_V%VA(19,10)=-5.263157894736842E-002_DP
 S_B_FROM_V%VB(19,14)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(19,19)=9.00000000000000E0_DP
 S_B_FROM_V%VB(19,25)=51.0000000000000E0_DP
 S_B_FROM_V%VA(19,32)=-204.000000000000E0_DP
 S_B_FROM_V%VB(19,40)=-612.000000000000E0_DP
 S_B_FROM_V%VA(19,49)=1428.00000000000E0_DP
 S_B_FROM_V%VB(19,59)=2652.00000000000E0_DP
 S_B_FROM_V%VA(19,70)=-3978.00000000000E0_DP
 S_B_FROM_V%VB(19,82)=-4862.00000000000E0_DP
 S_B_FROM_V%VA(19,95)=4862.00000000000E0_DP
 S_B_FROM_V%VB(19,109)=3978.00000000000E0_DP
 S_B_FROM_V%VA(19,124)=-2652.00000000000E0_DP
 S_B_FROM_V%VB(19,140)=-1428.00000000000E0_DP
 S_B_FROM_V%VA(19,157)=612.000000000000E0_DP
 S_B_FROM_V%VB(19,175)=204.000000000000E0_DP
 S_B_FROM_V%VA(19,194)=-51.0000000000000E0_DP
 S_B_FROM_V%VB(19,214)=-9.00000000000000E0_DP
 S_B_FROM_V%VA(19,235)=1.00000000000000E0_DP
 S_B_FROM_V%VB(19,257)=5.263157894736842E-002_DP
 S_B_FROM_V%VA(20,6)=-5.000000000000000E-002_DP
 S_B_FROM_V%VB(20,9)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(20,13)=9.50000000000000E0_DP
 S_B_FROM_V%VB(20,18)=57.0000000000000E0_DP
 S_B_FROM_V%VA(20,24)=-242.250000000000E0_DP
 S_B_FROM_V%VB(20,31)=-775.200000000000E0_DP
 S_B_FROM_V%VA(20,39)=1938.00000000000E0_DP
 S_B_FROM_V%VB(20,48)=3876.00000000000E0_DP
 S_B_FROM_V%VA(20,58)=-6298.50000000000E0_DP
 S_B_FROM_V%VB(20,69)=-8398.00000000000E0_DP
 S_B_FROM_V%VA(20,81)=9237.80000000000E0_DP
 S_B_FROM_V%VB(20,94)=8398.00000000000E0_DP
 S_B_FROM_V%VA(20,108)=-6298.50000000000E0_DP
 S_B_FROM_V%VB(20,123)=-3876.00000000000E0_DP
 S_B_FROM_V%VA(20,139)=1938.00000000000E0_DP
 S_B_FROM_V%VB(20,156)=775.200000000000E0_DP
 S_B_FROM_V%VA(20,174)=-242.250000000000E0_DP
 S_B_FROM_V%VB(20,193)=-57.0000000000000E0_DP
 S_B_FROM_V%VA(20,213)=9.50000000000000E0_DP
 S_B_FROM_V%VB(20,234)=1.00000000000000E0_DP
 S_B_FROM_V%VA(20,256)=-5.000000000000000E-002_DP
 S_B_FROM_V%VA(21,3)=-4.761904761904762E-002_DP
 S_B_FROM_V%VB(21,5)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(21,8)=10.0000000000000E0_DP
 S_B_FROM_V%VB(21,12)=63.3333333333333E0_DP
 S_B_FROM_V%VA(21,17)=-285.000000000000E0_DP
 S_B_FROM_V%VB(21,23)=-969.000000000000E0_DP
 S_B_FROM_V%VA(21,30)=2584.00000000000E0_DP
 S_B_FROM_V%VB(21,38)=5537.14285714286E0_DP
 S_B_FROM_V%VA(21,47)=-9690.00000000000E0_DP
 S_B_FROM_V%VB(21,57)=-13996.6666666667E0_DP
 S_B_FROM_V%VA(21,68)=16796.0000000000E0_DP
 S_B_FROM_V%VB(21,80)=16796.0000000000E0_DP
 S_B_FROM_V%VA(21,93)=-13996.6666666667E0_DP
 S_B_FROM_V%VB(21,107)=-9690.00000000000E0_DP
 S_B_FROM_V%VA(21,122)=5537.14285714286E0_DP
 S_B_FROM_V%VB(21,138)=2584.00000000000E0_DP
 S_B_FROM_V%VA(21,155)=-969.000000000000E0_DP
 S_B_FROM_V%VB(21,173)=-285.000000000000E0_DP
 S_B_FROM_V%VA(21,192)=63.3333333333333E0_DP
 S_B_FROM_V%VB(21,212)=10.0000000000000E0_DP
 S_B_FROM_V%VA(21,233)=-1.00000000000000E0_DP
 S_B_FROM_V%VB(21,255)=-4.761904761904762E-002_DP
 S_B_FROM_V%VA(22,1)=-4.545454545454546E-002_DP
 S_B_FROM_V%VB(22,2)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(22,4)=10.5000000000000E0_DP
 S_B_FROM_V%VB(22,7)=70.0000000000000E0_DP
 S_B_FROM_V%VA(22,11)=-332.500000000000E0_DP
 S_B_FROM_V%VB(22,16)=-1197.00000000000E0_DP
 S_B_FROM_V%VA(22,22)=3391.50000000000E0_DP
 S_B_FROM_V%VB(22,29)=7752.00000000000E0_DP
 S_B_FROM_V%VA(22,37)=-14535.0000000000E0_DP
 S_B_FROM_V%VB(22,46)=-22610.0000000000E0_DP
 S_B_FROM_V%VA(22,56)=29393.0000000000E0_DP
 S_B_FROM_V%VB(22,67)=32065.0909090909E0_DP
 S_B_FROM_V%VA(22,79)=-29393.0000000000E0_DP
 S_B_FROM_V%VB(22,92)=-22610.0000000000E0_DP
 S_B_FROM_V%VA(22,106)=14535.0000000000E0_DP
 S_B_FROM_V%VB(22,121)=7752.00000000000E0_DP
 S_B_FROM_V%VA(22,137)=-3391.50000000000E0_DP
 S_B_FROM_V%VB(22,154)=-1197.00000000000E0_DP
 S_B_FROM_V%VA(22,172)=332.500000000000E0_DP
 S_B_FROM_V%VB(22,191)=70.0000000000000E0_DP
 S_B_FROM_V%VA(22,211)=-10.5000000000000E0_DP
 S_B_FROM_V%VB(22,232)=-1.00000000000000E0_DP
 S_B_FROM_V%VA(22,254)=4.545454545454546E-002_DP

end subroutine set_s_b_mcmillan

subroutine set_s_e_mcmillan
  implicit none
  integer i
 

i=SECTOR_NMUL_MAX
!  ALLOCATE(S_e(SECTOR_NMUL_MAX))
!           DO I=1,SECTOR_NMUL_MAX
!           DO I=SECTOR_NMUL_MAX,SECTOR_NMUL_MAX   ! etienne 5/10/2015
             S_e%firsttime = 1  ! Piotr 8.19.2014 
             call nul_coef(S_e)
             call make_set_coef(S_e,I,0)
             S_e%firsttime=0
!          ENDDO

 S_E%I(1)=22;S_E%J(1)=0;
 S_E%I(2)=21;S_E%J(2)=1;
 S_E%I(3)=21;S_E%J(3)=0;
 S_E%I(4)=20;S_E%J(4)=2;
 S_E%I(5)=20;S_E%J(5)=1;
 S_E%I(6)=20;S_E%J(6)=0;
 S_E%I(7)=19;S_E%J(7)=3;
 S_E%I(8)=19;S_E%J(8)=2;
 S_E%I(9)=19;S_E%J(9)=1;
 S_E%I(10)=19;S_E%J(10)=0;
 S_E%I(11)=18;S_E%J(11)=4;
 S_E%I(12)=18;S_E%J(12)=3;
 S_E%I(13)=18;S_E%J(13)=2;
 S_E%I(14)=18;S_E%J(14)=1;
 S_E%I(15)=18;S_E%J(15)=0;
 S_E%I(16)=17;S_E%J(16)=5;
 S_E%I(17)=17;S_E%J(17)=4;
 S_E%I(18)=17;S_E%J(18)=3;
 S_E%I(19)=17;S_E%J(19)=2;
 S_E%I(20)=17;S_E%J(20)=1;
 S_E%I(21)=17;S_E%J(21)=0;
 S_E%I(22)=16;S_E%J(22)=6;
 S_E%I(23)=16;S_E%J(23)=5;
 S_E%I(24)=16;S_E%J(24)=4;
 S_E%I(25)=16;S_E%J(25)=3;
 S_E%I(26)=16;S_E%J(26)=2;
 S_E%I(27)=16;S_E%J(27)=1;
 S_E%I(28)=16;S_E%J(28)=0;
 S_E%I(29)=15;S_E%J(29)=7;
 S_E%I(30)=15;S_E%J(30)=6;
 S_E%I(31)=15;S_E%J(31)=5;
 S_E%I(32)=15;S_E%J(32)=4;
 S_E%I(33)=15;S_E%J(33)=3;
 S_E%I(34)=15;S_E%J(34)=2;
 S_E%I(35)=15;S_E%J(35)=1;
 S_E%I(36)=15;S_E%J(36)=0;
 S_E%I(37)=14;S_E%J(37)=8;
 S_E%I(38)=14;S_E%J(38)=7;
 S_E%I(39)=14;S_E%J(39)=6;
 S_E%I(40)=14;S_E%J(40)=5;
 S_E%I(41)=14;S_E%J(41)=4;
 S_E%I(42)=14;S_E%J(42)=3;
 S_E%I(43)=14;S_E%J(43)=2;
 S_E%I(44)=14;S_E%J(44)=1;
 S_E%I(45)=14;S_E%J(45)=0;
 S_E%I(46)=13;S_E%J(46)=9;
 S_E%I(47)=13;S_E%J(47)=8;
 S_E%I(48)=13;S_E%J(48)=7;
 S_E%I(49)=13;S_E%J(49)=6;
 S_E%I(50)=13;S_E%J(50)=5;
 S_E%I(51)=13;S_E%J(51)=4;
 S_E%I(52)=13;S_E%J(52)=3;
 S_E%I(53)=13;S_E%J(53)=2;
 S_E%I(54)=13;S_E%J(54)=1;
 S_E%I(55)=13;S_E%J(55)=0;
 S_E%I(56)=12;S_E%J(56)=10;
 S_E%I(57)=12;S_E%J(57)=9;
 S_E%I(58)=12;S_E%J(58)=8;
 S_E%I(59)=12;S_E%J(59)=7;
 S_E%I(60)=12;S_E%J(60)=6;
 S_E%I(61)=12;S_E%J(61)=5;
 S_E%I(62)=12;S_E%J(62)=4;
 S_E%I(63)=12;S_E%J(63)=3;
 S_E%I(64)=12;S_E%J(64)=2;
 S_E%I(65)=12;S_E%J(65)=1;
 S_E%I(66)=12;S_E%J(66)=0;
 S_E%I(67)=11;S_E%J(67)=11;
 S_E%I(68)=11;S_E%J(68)=10;
 S_E%I(69)=11;S_E%J(69)=9;
 S_E%I(70)=11;S_E%J(70)=8;
 S_E%I(71)=11;S_E%J(71)=7;
 S_E%I(72)=11;S_E%J(72)=6;
 S_E%I(73)=11;S_E%J(73)=5;
 S_E%I(74)=11;S_E%J(74)=4;
 S_E%I(75)=11;S_E%J(75)=3;
 S_E%I(76)=11;S_E%J(76)=2;
 S_E%I(77)=11;S_E%J(77)=1;
 S_E%I(78)=11;S_E%J(78)=0;
 S_E%I(79)=10;S_E%J(79)=12;
 S_E%I(80)=10;S_E%J(80)=11;
 S_E%I(81)=10;S_E%J(81)=10;
 S_E%I(82)=10;S_E%J(82)=9;
 S_E%I(83)=10;S_E%J(83)=8;
 S_E%I(84)=10;S_E%J(84)=7;
 S_E%I(85)=10;S_E%J(85)=6;
 S_E%I(86)=10;S_E%J(86)=5;
 S_E%I(87)=10;S_E%J(87)=4;
 S_E%I(88)=10;S_E%J(88)=3;
 S_E%I(89)=10;S_E%J(89)=2;
 S_E%I(90)=10;S_E%J(90)=1;
 S_E%I(91)=10;S_E%J(91)=0;
 S_E%I(92)=9;S_E%J(92)=13;
 S_E%I(93)=9;S_E%J(93)=12;
 S_E%I(94)=9;S_E%J(94)=11;
 S_E%I(95)=9;S_E%J(95)=10;
 S_E%I(96)=9;S_E%J(96)=9;
 S_E%I(97)=9;S_E%J(97)=8;
 S_E%I(98)=9;S_E%J(98)=7;
 S_E%I(99)=9;S_E%J(99)=6;
 S_E%I(100)=9;S_E%J(100)=5;
 S_E%I(101)=9;S_E%J(101)=4;
 S_E%I(102)=9;S_E%J(102)=3;
 S_E%I(103)=9;S_E%J(103)=2;
 S_E%I(104)=9;S_E%J(104)=1;
 S_E%I(105)=9;S_E%J(105)=0;
 S_E%I(106)=8;S_E%J(106)=14;
 S_E%I(107)=8;S_E%J(107)=13;
 S_E%I(108)=8;S_E%J(108)=12;
 S_E%I(109)=8;S_E%J(109)=11;
 S_E%I(110)=8;S_E%J(110)=10;
 S_E%I(111)=8;S_E%J(111)=9;
 S_E%I(112)=8;S_E%J(112)=8;
 S_E%I(113)=8;S_E%J(113)=7;
 S_E%I(114)=8;S_E%J(114)=6;
 S_E%I(115)=8;S_E%J(115)=5;
 S_E%I(116)=8;S_E%J(116)=4;
 S_E%I(117)=8;S_E%J(117)=3;
 S_E%I(118)=8;S_E%J(118)=2;
 S_E%I(119)=8;S_E%J(119)=1;
 S_E%I(120)=8;S_E%J(120)=0;
 S_E%I(121)=7;S_E%J(121)=15;
 S_E%I(122)=7;S_E%J(122)=14;
 S_E%I(123)=7;S_E%J(123)=13;
 S_E%I(124)=7;S_E%J(124)=12;
 S_E%I(125)=7;S_E%J(125)=11;
 S_E%I(126)=7;S_E%J(126)=10;
 S_E%I(127)=7;S_E%J(127)=9;
 S_E%I(128)=7;S_E%J(128)=8;
 S_E%I(129)=7;S_E%J(129)=7;
 S_E%I(130)=7;S_E%J(130)=6;
 S_E%I(131)=7;S_E%J(131)=5;
 S_E%I(132)=7;S_E%J(132)=4;
 S_E%I(133)=7;S_E%J(133)=3;
 S_E%I(134)=7;S_E%J(134)=2;
 S_E%I(135)=7;S_E%J(135)=1;
 S_E%I(136)=7;S_E%J(136)=0;
 S_E%I(137)=6;S_E%J(137)=16;
 S_E%I(138)=6;S_E%J(138)=15;
 S_E%I(139)=6;S_E%J(139)=14;
 S_E%I(140)=6;S_E%J(140)=13;
 S_E%I(141)=6;S_E%J(141)=12;
 S_E%I(142)=6;S_E%J(142)=11;
 S_E%I(143)=6;S_E%J(143)=10;
 S_E%I(144)=6;S_E%J(144)=9;
 S_E%I(145)=6;S_E%J(145)=8;
 S_E%I(146)=6;S_E%J(146)=7;
 S_E%I(147)=6;S_E%J(147)=6;
 S_E%I(148)=6;S_E%J(148)=5;
 S_E%I(149)=6;S_E%J(149)=4;
 S_E%I(150)=6;S_E%J(150)=3;
 S_E%I(151)=6;S_E%J(151)=2;
 S_E%I(152)=6;S_E%J(152)=1;
 S_E%I(153)=6;S_E%J(153)=0;
 S_E%I(154)=5;S_E%J(154)=17;
 S_E%I(155)=5;S_E%J(155)=16;
 S_E%I(156)=5;S_E%J(156)=15;
 S_E%I(157)=5;S_E%J(157)=14;
 S_E%I(158)=5;S_E%J(158)=13;
 S_E%I(159)=5;S_E%J(159)=12;
 S_E%I(160)=5;S_E%J(160)=11;
 S_E%I(161)=5;S_E%J(161)=10;
 S_E%I(162)=5;S_E%J(162)=9;
 S_E%I(163)=5;S_E%J(163)=8;
 S_E%I(164)=5;S_E%J(164)=7;
 S_E%I(165)=5;S_E%J(165)=6;
 S_E%I(166)=5;S_E%J(166)=5;
 S_E%I(167)=5;S_E%J(167)=4;
 S_E%I(168)=5;S_E%J(168)=3;
 S_E%I(169)=5;S_E%J(169)=2;
 S_E%I(170)=5;S_E%J(170)=1;
 S_E%I(171)=5;S_E%J(171)=0;
 S_E%I(172)=4;S_E%J(172)=18;
 S_E%I(173)=4;S_E%J(173)=17;
 S_E%I(174)=4;S_E%J(174)=16;
 S_E%I(175)=4;S_E%J(175)=15;
 S_E%I(176)=4;S_E%J(176)=14;
 S_E%I(177)=4;S_E%J(177)=13;
 S_E%I(178)=4;S_E%J(178)=12;
 S_E%I(179)=4;S_E%J(179)=11;
 S_E%I(180)=4;S_E%J(180)=10;
 S_E%I(181)=4;S_E%J(181)=9;
 S_E%I(182)=4;S_E%J(182)=8;
 S_E%I(183)=4;S_E%J(183)=7;
 S_E%I(184)=4;S_E%J(184)=6;
 S_E%I(185)=4;S_E%J(185)=5;
 S_E%I(186)=4;S_E%J(186)=4;
 S_E%I(187)=4;S_E%J(187)=3;
 S_E%I(188)=4;S_E%J(188)=2;
 S_E%I(189)=4;S_E%J(189)=1;
 S_E%I(190)=4;S_E%J(190)=0;
 S_E%I(191)=3;S_E%J(191)=19;
 S_E%I(192)=3;S_E%J(192)=18;
 S_E%I(193)=3;S_E%J(193)=17;
 S_E%I(194)=3;S_E%J(194)=16;
 S_E%I(195)=3;S_E%J(195)=15;
 S_E%I(196)=3;S_E%J(196)=14;
 S_E%I(197)=3;S_E%J(197)=13;
 S_E%I(198)=3;S_E%J(198)=12;
 S_E%I(199)=3;S_E%J(199)=11;
 S_E%I(200)=3;S_E%J(200)=10;
 S_E%I(201)=3;S_E%J(201)=9;
 S_E%I(202)=3;S_E%J(202)=8;
 S_E%I(203)=3;S_E%J(203)=7;
 S_E%I(204)=3;S_E%J(204)=6;
 S_E%I(205)=3;S_E%J(205)=5;
 S_E%I(206)=3;S_E%J(206)=4;
 S_E%I(207)=3;S_E%J(207)=3;
 S_E%I(208)=3;S_E%J(208)=2;
 S_E%I(209)=3;S_E%J(209)=1;
 S_E%I(210)=3;S_E%J(210)=0;
 S_E%I(211)=2;S_E%J(211)=20;
 S_E%I(212)=2;S_E%J(212)=19;
 S_E%I(213)=2;S_E%J(213)=18;
 S_E%I(214)=2;S_E%J(214)=17;
 S_E%I(215)=2;S_E%J(215)=16;
 S_E%I(216)=2;S_E%J(216)=15;
 S_E%I(217)=2;S_E%J(217)=14;
 S_E%I(218)=2;S_E%J(218)=13;
 S_E%I(219)=2;S_E%J(219)=12;
 S_E%I(220)=2;S_E%J(220)=11;
 S_E%I(221)=2;S_E%J(221)=10;
 S_E%I(222)=2;S_E%J(222)=9;
 S_E%I(223)=2;S_E%J(223)=8;
 S_E%I(224)=2;S_E%J(224)=7;
 S_E%I(225)=2;S_E%J(225)=6;
 S_E%I(226)=2;S_E%J(226)=5;
 S_E%I(227)=2;S_E%J(227)=4;
 S_E%I(228)=2;S_E%J(228)=3;
 S_E%I(229)=2;S_E%J(229)=2;
 S_E%I(230)=2;S_E%J(230)=1;
 S_E%I(231)=2;S_E%J(231)=0;
 S_E%I(232)=1;S_E%J(232)=21;
 S_E%I(233)=1;S_E%J(233)=20;
 S_E%I(234)=1;S_E%J(234)=19;
 S_E%I(235)=1;S_E%J(235)=18;
 S_E%I(236)=1;S_E%J(236)=17;
 S_E%I(237)=1;S_E%J(237)=16;
 S_E%I(238)=1;S_E%J(238)=15;
 S_E%I(239)=1;S_E%J(239)=14;
 S_E%I(240)=1;S_E%J(240)=13;
 S_E%I(241)=1;S_E%J(241)=12;
 S_E%I(242)=1;S_E%J(242)=11;
 S_E%I(243)=1;S_E%J(243)=10;
 S_E%I(244)=1;S_E%J(244)=9;
 S_E%I(245)=1;S_E%J(245)=8;
 S_E%I(246)=1;S_E%J(246)=7;
 S_E%I(247)=1;S_E%J(247)=6;
 S_E%I(248)=1;S_E%J(248)=5;
 S_E%I(249)=1;S_E%J(249)=4;
 S_E%I(250)=1;S_E%J(250)=3;
 S_E%I(251)=1;S_E%J(251)=2;
 S_E%I(252)=1;S_E%J(252)=1;
 S_E%I(253)=1;S_E%J(253)=0;
 S_E%I(254)=0;S_E%J(254)=22;
 S_E%I(255)=0;S_E%J(255)=21;
 S_E%I(256)=0;S_E%J(256)=20;
 S_E%I(257)=0;S_E%J(257)=19;
 S_E%I(258)=0;S_E%J(258)=18;
 S_E%I(259)=0;S_E%J(259)=17;
 S_E%I(260)=0;S_E%J(260)=16;
 S_E%I(261)=0;S_E%J(261)=15;
 S_E%I(262)=0;S_E%J(262)=14;
 S_E%I(263)=0;S_E%J(263)=13;
 S_E%I(264)=0;S_E%J(264)=12;
 S_E%I(265)=0;S_E%J(265)=11;
 S_E%I(266)=0;S_E%J(266)=10;
 S_E%I(267)=0;S_E%J(267)=9;
 S_E%I(268)=0;S_E%J(268)=8;
 S_E%I(269)=0;S_E%J(269)=7;
 S_E%I(270)=0;S_E%J(270)=6;
 S_E%I(271)=0;S_E%J(271)=5;
 S_E%I(272)=0;S_E%J(272)=4;
 S_E%I(273)=0;S_E%J(273)=3;
 S_E%I(274)=0;S_E%J(274)=2;
 S_E%I(275)=0;S_E%J(275)=1;
 S_E%I(276)=0;S_E%J(276)=0;
 S_E%B_X(1,91)=1.00000000000000E0_DP
 S_E%B_X(1,105)=-1.00000000000000E0_DP
 S_E%B_X(1,120)=1.00000000000000E0_DP
 S_E%B_X(1,136)=-1.00000000000000E0_DP
 S_E%B_X(1,153)=1.00000000000000E0_DP
 S_E%B_X(1,171)=-1.00000000000000E0_DP
 S_E%B_X(1,190)=1.00000000000000E0_DP
 S_E%B_X(1,210)=-1.00000000000000E0_DP
 S_E%B_X(1,231)=1.00000000000000E0_DP
 S_E%B_X(1,253)=-1.00000000000000E0_DP
 S_E%B_X(1,276)=1.00000000000000E0_DP
 S_E%A_Y(1,276)=1.00000000000000E0_DP
 S_E%B_X(2,91)=-0.500000000000000E0_DP
 S_E%A_Y(2,91)=-0.100000000000000E0_DP
 S_E%A_X(2,104)=-1.00000000000000E0_DP
 S_E%B_X(2,105)=0.500000000000000E0_DP
 S_E%A_Y(2,105)=0.111111111111111E0_DP
 S_E%A_X(2,119)=1.00000000000000E0_DP
 S_E%B_X(2,120)=-0.500000000000000E0_DP
 S_E%A_Y(2,120)=-0.125000000000000E0_DP
 S_E%A_X(2,135)=-1.00000000000000E0_DP
 S_E%B_X(2,136)=0.500000000000000E0_DP
 S_E%A_Y(2,136)=0.142857142857143E0_DP
 S_E%A_X(2,152)=1.00000000000000E0_DP
 S_E%B_X(2,153)=-0.500000000000000E0_DP
 S_E%A_Y(2,153)=-0.166666666666667E0_DP
 S_E%A_X(2,170)=-1.00000000000000E0_DP
 S_E%B_X(2,171)=0.500000000000000E0_DP
 S_E%A_Y(2,171)=0.200000000000000E0_DP
 S_E%A_X(2,189)=1.00000000000000E0_DP
 S_E%B_X(2,190)=-0.500000000000000E0_DP
 S_E%A_Y(2,190)=-0.250000000000000E0_DP
 S_E%A_X(2,209)=-1.00000000000000E0_DP
 S_E%B_X(2,210)=0.500000000000000E0_DP
 S_E%A_Y(2,210)=0.333333333333333E0_DP
 S_E%A_X(2,230)=1.00000000000000E0_DP
 S_E%B_X(2,231)=-0.500000000000000E0_DP
 S_E%A_Y(2,231)=-0.500000000000000E0_DP
 S_E%A_X(2,252)=-1.00000000000000E0_DP
 S_E%B_X(2,253)=1.00000000000000E0_DP
 S_E%A_Y(2,253)=1.00000000000000E0_DP
 S_E%A_X(2,275)=1.00000000000000E0_DP
 S_E%B_Y(2,275)=-1.00000000000000E0_DP
 S_E%B_X(3,91)=0.511111111111111E0_DP
 S_E%A_Y(3,91)=0.100000000000000E0_DP
 S_E%A_X(3,104)=1.00000000000000E0_DP
 S_E%B_Y(3,104)=-0.222222222222221E0_DP
 S_E%B_X(3,105)=-0.513888888888889E0_DP
 S_E%A_Y(3,105)=-0.111111111111111E0_DP
 S_E%B_X(3,118)=-0.999999999999996E0_DP
 S_E%A_X(3,119)=-1.00000000000000E0_DP
 S_E%B_Y(3,119)=0.250000000000000E0_DP
 S_E%B_X(3,120)=0.517857142857143E0_DP
 S_E%A_Y(3,120)=0.125000000000000E0_DP
 S_E%B_X(3,134)=1.00000000000000E0_DP
 S_E%A_X(3,135)=1.00000000000000E0_DP
 S_E%B_Y(3,135)=-0.285714285714286E0_DP
 S_E%B_X(3,136)=-0.523809523809524E0_DP
 S_E%A_Y(3,136)=-0.142857142857143E0_DP
 S_E%B_X(3,151)=-1.00000000000000E0_DP
 S_E%A_X(3,152)=-1.00000000000000E0_DP
 S_E%B_Y(3,152)=0.333333333333333E0_DP
 S_E%B_X(3,153)=0.533333333333333E0_DP
 S_E%A_Y(3,153)=0.166666666666667E0_DP
 S_E%B_X(3,169)=0.999999999999999E0_DP
 S_E%A_X(3,170)=1.00000000000000E0_DP
 S_E%B_Y(3,170)=-0.399999999999999E0_DP
 S_E%B_X(3,171)=-0.550000000000000E0_DP
 S_E%A_Y(3,171)=-0.200000000000000E0_DP
 S_E%B_X(3,188)=-0.999999999999999E0_DP
 S_E%A_X(3,189)=-1.00000000000000E0_DP
 S_E%B_Y(3,189)=0.500000000000000E0_DP
 S_E%B_X(3,190)=0.583333333333333E0_DP
 S_E%A_Y(3,190)=0.250000000000000E0_DP
 S_E%B_X(3,208)=1.00000000000000E0_DP
 S_E%A_X(3,209)=1.00000000000000E0_DP
 S_E%B_Y(3,209)=-0.666666666666667E0_DP
 S_E%B_X(3,210)=-0.666666666666667E0_DP
 S_E%A_Y(3,210)=-0.333333333333333E0_DP
 S_E%B_X(3,229)=-1.00000000000000E0_DP
 S_E%A_X(3,230)=-1.00000000000000E0_DP
 S_E%B_Y(3,230)=1.00000000000000E0_DP
 S_E%B_X(3,231)=1.00000000000000E0_DP
 S_E%A_Y(3,231)=1.00000000000000E0_DP
 S_E%B_X(3,251)=1.00000000000000E0_DP
 S_E%A_X(3,252)=2.00000000000000E0_DP
 S_E%B_Y(3,252)=-2.00000000000000E0_DP
 S_E%B_X(3,274)=-1.00000000000000E0_DP
 S_E%A_Y(3,274)=-1.00000000000000E0_DP
 S_E%B_X(4,91)=-0.391666666666667E0_DP
 S_E%A_Y(4,91)=-0.154166666666667E0_DP
 S_E%A_X(4,104)=-1.54166666666667E0_DP
 S_E%B_Y(4,104)=0.333333333333333E0_DP
 S_E%B_X(4,105)=0.395833333333333E0_DP
 S_E%A_Y(4,105)=0.172619047619048E0_DP
 S_E%B_X(4,118)=1.50000000000000E0_DP
 S_E%A_Y(4,118)=0.375000000000000E0_DP
 S_E%A_X(4,119)=1.55357142857143E0_DP
 S_E%B_Y(4,119)=-0.375000000000000E0_DP
 S_E%B_X(4,120)=-0.401785714285714E0_DP
 S_E%A_Y(4,120)=-0.196428571428571E0_DP
 S_E%A_X(4,133)=1.00000000000000E0_DP
 S_E%B_X(4,134)=-1.50000000000000E0_DP
 S_E%A_Y(4,134)=-0.428571428571428E0_DP
 S_E%A_X(4,135)=-1.57142857142857E0_DP
 S_E%B_Y(4,135)=0.428571428571429E0_DP
 S_E%B_X(4,136)=0.410714285714286E0_DP
 S_E%A_Y(4,136)=0.228571428571429E0_DP
 S_E%A_X(4,150)=-1.00000000000000E0_DP
 S_E%B_X(4,151)=1.50000000000000E0_DP
 S_E%A_Y(4,151)=0.500000000000000E0_DP
 S_E%A_X(4,152)=1.60000000000000E0_DP
 S_E%B_Y(4,152)=-0.500000000000000E0_DP
 S_E%B_X(4,153)=-0.425000000000000E0_DP
 S_E%A_Y(4,153)=-0.275000000000000E0_DP
 S_E%A_X(4,168)=0.999999999999999E0_DP
 S_E%B_X(4,169)=-1.50000000000000E0_DP
 S_E%A_Y(4,169)=-0.600000000000000E0_DP
 S_E%A_X(4,170)=-1.65000000000000E0_DP
 S_E%B_Y(4,170)=0.600000000000000E0_DP
 S_E%B_X(4,171)=0.450000000000000E0_DP
 S_E%A_Y(4,171)=0.350000000000000E0_DP
 S_E%A_X(4,187)=-1.00000000000000E0_DP
 S_E%B_X(4,188)=1.50000000000000E0_DP
 S_E%A_Y(4,188)=0.750000000000000E0_DP
 S_E%A_X(4,189)=1.75000000000000E0_DP
 S_E%B_Y(4,189)=-0.750000000000000E0_DP
 S_E%B_X(4,190)=-0.500000000000000E0_DP
 S_E%A_Y(4,190)=-0.500000000000000E0_DP
 S_E%A_X(4,207)=1.00000000000000E0_DP
 S_E%B_X(4,208)=-1.50000000000000E0_DP
 S_E%A_Y(4,208)=-1.00000000000000E0_DP
 S_E%A_X(4,209)=-2.00000000000000E0_DP
 S_E%B_Y(4,209)=1.00000000000000E0_DP
 S_E%B_X(4,210)=1.00000000000000E0_DP
 S_E%A_Y(4,210)=1.00000000000000E0_DP
 S_E%A_X(4,228)=-1.00000000000000E0_DP
 S_E%B_X(4,229)=1.50000000000000E0_DP
 S_E%A_Y(4,229)=1.50000000000000E0_DP
 S_E%A_X(4,230)=3.00000000000000E0_DP
 S_E%B_Y(4,230)=-3.00000000000000E0_DP
 S_E%A_X(4,250)=1.00000000000000E0_DP
 S_E%B_X(4,251)=-3.00000000000000E0_DP
 S_E%A_Y(4,251)=-3.00000000000000E0_DP
 S_E%A_X(4,273)=-1.00000000000000E0_DP
 S_E%B_Y(4,273)=1.00000000000000E0_DP
 S_E%B_X(5,91)=0.410119047619048E0_DP
 S_E%A_Y(5,91)=0.158333333333333E0_DP
 S_E%A_X(5,104)=1.58333333333333E0_DP
 S_E%B_Y(5,104)=-0.690476190476191E0_DP
 S_E%B_X(5,105)=-0.419642857142857E0_DP
 S_E%A_Y(5,105)=-0.178571428571429E0_DP
 S_E%B_X(5,118)=-3.10714285714286E0_DP
 S_E%A_Y(5,118)=-0.750000000000000E0_DP
 S_E%A_X(5,119)=-1.60714285714286E0_DP
 S_E%B_Y(5,119)=0.785714285714285E0_DP
 S_E%B_X(5,120)=0.433928571428571E0_DP
 S_E%A_Y(5,120)=0.205357142857143E0_DP
 S_E%A_X(5,133)=-2.00000000000000E0_DP
 S_E%B_Y(5,133)=0.571428571428577E0_DP
 S_E%B_X(5,134)=3.14285714285714E0_DP
 S_E%A_Y(5,134)=0.857142857142856E0_DP
 S_E%A_X(5,135)=1.64285714285714E0_DP
 S_E%B_Y(5,135)=-0.914285714285714E0_DP
 S_E%B_X(5,136)=-0.457142857142857E0_DP
 S_E%A_Y(5,136)=-0.242857142857143E0_DP
 S_E%B_X(5,149)=1.00000000000001E0_DP
 S_E%A_X(5,150)=2.00000000000000E0_DP
 S_E%B_Y(5,150)=-0.666666666666657E0_DP
 S_E%B_X(5,151)=-3.20000000000000E0_DP
 S_E%A_Y(5,151)=-0.999999999999999E0_DP
 S_E%A_X(5,152)=-1.70000000000000E0_DP
 S_E%B_Y(5,152)=1.10000000000000E0_DP
 S_E%B_X(5,153)=0.500000000000000E0_DP
 S_E%A_Y(5,153)=0.300000000000000E0_DP
 S_E%B_X(5,167)=-0.999999999999986E0_DP
 S_E%A_X(5,168)=-2.00000000000000E0_DP
 S_E%B_Y(5,168)=0.800000000000001E0_DP
 S_E%B_X(5,169)=3.30000000000000E0_DP
 S_E%A_Y(5,169)=1.20000000000000E0_DP
 S_E%A_X(5,170)=1.80000000000000E0_DP
 S_E%B_Y(5,170)=-1.40000000000000E0_DP
 S_E%B_X(5,171)=-0.600000000000000E0_DP
 S_E%A_Y(5,171)=-0.400000000000000E0_DP
 S_E%B_X(5,186)=1.00000000000000E0_DP
 S_E%A_X(5,187)=2.00000000000000E0_DP
 S_E%B_Y(5,187)=-1.00000000000000E0_DP
 S_E%B_X(5,188)=-3.50000000000000E0_DP
 S_E%A_Y(5,188)=-1.50000000000000E0_DP
 S_E%A_X(5,189)=-2.00000000000000E0_DP
 S_E%B_Y(5,189)=2.00000000000000E0_DP
 S_E%B_X(5,190)=1.00000000000000E0_DP
 S_E%A_Y(5,190)=1.00000000000000E0_DP
 S_E%B_X(5,206)=-1.00000000000000E0_DP
 S_E%A_X(5,207)=-2.00000000000000E0_DP
 S_E%B_Y(5,207)=1.33333333333333E0_DP
 S_E%B_X(5,208)=4.00000000000000E0_DP
 S_E%A_Y(5,208)=2.00000000000000E0_DP
 S_E%A_X(5,209)=4.00000000000000E0_DP
 S_E%B_Y(5,209)=-4.00000000000000E0_DP
 S_E%B_X(5,227)=1.00000000000000E0_DP
 S_E%A_X(5,228)=2.00000000000000E0_DP
 S_E%B_Y(5,228)=-2.00000000000000E0_DP
 S_E%B_X(5,229)=-6.00000000000000E0_DP
 S_E%A_Y(5,229)=-6.00000000000000E0_DP
 S_E%B_X(5,249)=-1.00000000000000E0_DP
 S_E%A_X(5,250)=-4.00000000000000E0_DP
 S_E%B_Y(5,250)=4.00000000000000E0_DP
 S_E%B_X(5,272)=1.00000000000000E0_DP
 S_E%A_Y(5,272)=1.00000000000000E0_DP
 S_E%B_X(6,91)=-0.358630952380952E0_DP
 S_E%A_Y(6,91)=-0.209821428571429E0_DP
 S_E%A_X(6,104)=-2.09821428571429E0_DP
 S_E%B_Y(6,104)=0.892857142857143E0_DP
 S_E%B_X(6,105)=0.372023809523809E0_DP
 S_E%A_Y(6,105)=0.241071428571429E0_DP
 S_E%B_X(6,118)=4.01785714285714E0_DP
 S_E%A_Y(6,118)=1.96428571428571E0_DP
 S_E%A_X(6,119)=2.16964285714286E0_DP
 S_E%B_Y(6,119)=-1.02678571428571E0_DP
 S_E%B_X(6,120)=-0.392857142857143E0_DP
 S_E%A_Y(6,120)=-0.285714285714286E0_DP
 S_E%A_X(6,133)=5.23809523809524E0_DP
 S_E%B_Y(6,133)=-1.42857142857143E0_DP
 S_E%B_X(6,134)=-4.10714285714285E0_DP
 S_E%A_Y(6,134)=-2.28571428571428E0_DP
 S_E%A_X(6,135)=-2.28571428571429E0_DP
 S_E%B_Y(6,135)=1.21428571428571E0_DP
 S_E%B_X(6,136)=0.428571428571428E0_DP
 S_E%A_Y(6,136)=0.357142857142857E0_DP
 S_E%B_X(6,149)=-2.50000000000000E0_DP
 S_E%A_Y(6,149)=-0.833333333333330E0_DP
 S_E%A_X(6,150)=-5.33333333333333E0_DP
 S_E%B_Y(6,150)=1.66666666666666E0_DP
 S_E%B_X(6,151)=4.25000000000000E0_DP
 S_E%A_Y(6,151)=2.75000000000000E0_DP
 S_E%A_X(6,152)=2.50000000000000E0_DP
 S_E%B_Y(6,152)=-1.50000000000000E0_DP
 S_E%B_X(6,153)=-0.500000000000000E0_DP
 S_E%A_Y(6,153)=-0.500000000000000E0_DP
 S_E%A_X(6,166)=-0.999999999999996E0_DP
 S_E%B_X(6,167)=2.49999999999999E0_DP
 S_E%A_Y(6,167)=1.00000000000000E0_DP
 S_E%A_X(6,168)=5.50000000000000E0_DP
 S_E%B_Y(6,168)=-2.00000000000000E0_DP
 S_E%B_X(6,169)=-4.50000000000000E0_DP
 S_E%A_Y(6,169)=-3.50000000000000E0_DP
 S_E%A_X(6,170)=-3.00000000000000E0_DP
 S_E%B_Y(6,170)=2.00000000000000E0_DP
 S_E%B_X(6,171)=1.00000000000000E0_DP
 S_E%A_Y(6,171)=1.00000000000000E0_DP
 S_E%A_X(6,185)=1.00000000000000E0_DP
 S_E%B_X(6,186)=-2.50000000000000E0_DP
 S_E%A_Y(6,186)=-1.25000000000000E0_DP
 S_E%A_X(6,187)=-5.83333333333333E0_DP
 S_E%B_Y(6,187)=2.50000000000000E0_DP
 S_E%B_X(6,188)=5.00000000000000E0_DP
 S_E%A_Y(6,188)=5.00000000000000E0_DP
 S_E%A_X(6,189)=5.00000000000000E0_DP
 S_E%B_Y(6,189)=-5.00000000000000E0_DP
 S_E%A_X(6,205)=-1.00000000000000E0_DP
 S_E%B_X(6,206)=2.50000000000000E0_DP
 S_E%A_Y(6,206)=1.66666666666666E0_DP
 S_E%A_X(6,207)=6.66666666666667E0_DP
 S_E%B_Y(6,207)=-3.33333333333333E0_DP
 S_E%B_X(6,208)=-10.0000000000000E0_DP
 S_E%A_Y(6,208)=-10.0000000000000E0_DP
 S_E%A_X(6,226)=0.999999999999998E0_DP
 S_E%B_X(6,227)=-2.50000000000000E0_DP
 S_E%A_Y(6,227)=-2.50000000000000E0_DP
 S_E%A_X(6,228)=-10.0000000000000E0_DP
 S_E%B_Y(6,228)=10.0000000000000E0_DP
 S_E%A_X(6,248)=-1.00000000000000E0_DP
 S_E%B_X(6,249)=5.00000000000000E0_DP
 S_E%A_Y(6,249)=5.00000000000000E0_DP
 S_E%A_X(6,271)=1.00000000000000E0_DP
 S_E%B_Y(6,271)=-1.00000000000000E0_DP
 S_E%B_X(7,91)=0.389880952380952E0_DP
 S_E%A_Y(7,91)=0.223214285714286E0_DP
 S_E%A_X(7,104)=2.23214285714286E0_DP
 S_E%B_Y(7,104)=-1.44642857142857E0_DP
 S_E%B_X(7,105)=-0.416666666666667E0_DP
 S_E%A_Y(7,105)=-0.261904761904762E0_DP
 S_E%B_X(7,118)=-6.50892857142857E0_DP
 S_E%A_Y(7,118)=-3.08035714285714E0_DP
 S_E%A_X(7,119)=-2.35714285714286E0_DP
 S_E%B_Y(7,119)=1.71428571428571E0_DP
 S_E%B_X(7,120)=0.464285714285714E0_DP
 S_E%A_Y(7,120)=0.321428571428571E0_DP
 S_E%A_X(7,133)=-8.21428571428572E0_DP
 S_E%B_Y(7,133)=4.57142857142857E0_DP
 S_E%B_X(7,134)=6.85714285714285E0_DP
 S_E%A_Y(7,134)=3.64285714285714E0_DP
 S_E%A_X(7,135)=2.57142857142857E0_DP
 S_E%B_Y(7,135)=-2.14285714285714E0_DP
 S_E%B_X(7,136)=-0.571428571428571E0_DP
 S_E%A_Y(7,136)=-0.428571428571429E0_DP
 S_E%B_X(7,149)=8.00000000000000E0_DP
 S_E%A_Y(7,149)=2.50000000000000E0_DP
 S_E%A_X(7,150)=8.50000000000000E0_DP
 S_E%B_Y(7,150)=-5.49999999999999E0_DP
 S_E%B_X(7,151)=-7.50000000000000E0_DP
 S_E%A_Y(7,151)=-4.50000000000000E0_DP
 S_E%A_X(7,152)=-3.00000000000000E0_DP
 S_E%B_Y(7,152)=3.00000000000000E0_DP
 S_E%B_X(7,153)=1.00000000000000E0_DP
 S_E%A_Y(7,153)=1.00000000000000E0_DP
 S_E%A_X(7,166)=3.00000000000000E0_DP
 S_E%B_Y(7,166)=-1.20000000000000E0_DP
 S_E%B_X(7,167)=-8.24999999999999E0_DP
 S_E%A_Y(7,167)=-3.00000000000000E0_DP
 S_E%A_X(7,168)=-9.00000000000000E0_DP
 S_E%B_Y(7,168)=7.00000000000000E0_DP
 S_E%B_X(7,169)=9.00000000000000E0_DP
 S_E%A_Y(7,169)=6.00000000000000E0_DP
 S_E%A_X(7,170)=6.00000000000000E0_DP
 S_E%B_Y(7,170)=-6.00000000000000E0_DP
 S_E%B_X(7,184)=-0.999999999999996E0_DP
 S_E%A_X(7,185)=-3.00000000000000E0_DP
 S_E%B_Y(7,185)=1.49999999999999E0_DP
 S_E%B_X(7,186)=8.75000000000000E0_DP
 S_E%A_Y(7,186)=3.75000000000000E0_DP
 S_E%A_X(7,187)=10.0000000000000E0_DP
 S_E%B_Y(7,187)=-10.0000000000000E0_DP
 S_E%B_X(7,188)=-15.0000000000000E0_DP
 S_E%A_Y(7,188)=-15.0000000000000E0_DP
 S_E%B_X(7,204)=0.999999999999993E0_DP
 S_E%A_X(7,205)=3.00000000000000E0_DP
 S_E%B_Y(7,205)=-2.00000000000000E0_DP
 S_E%B_X(7,206)=-10.0000000000000E0_DP
 S_E%A_Y(7,206)=-5.00000000000000E0_DP
 S_E%A_X(7,207)=-20.0000000000000E0_DP
 S_E%B_Y(7,207)=20.0000000000000E0_DP
 S_E%B_X(7,225)=-1.00000000000000E0_DP
 S_E%A_X(7,226)=-3.00000000000000E0_DP
 S_E%B_Y(7,226)=3.00000000000000E0_DP
 S_E%B_X(7,227)=15.0000000000000E0_DP
 S_E%A_Y(7,227)=15.0000000000000E0_DP
 S_E%B_X(7,247)=1.00000000000000E0_DP
 S_E%A_X(7,248)=6.00000000000000E0_DP
 S_E%B_Y(7,248)=-6.00000000000000E0_DP
 S_E%B_X(7,270)=-1.00000000000000E0_DP
 S_E%A_Y(7,270)=-1.00000000000000E0_DP
 S_E%B_X(8,91)=-0.375000000000000E0_DP
 S_E%A_Y(8,91)=-0.291666666666667E0_DP
 S_E%A_X(8,104)=-2.91666666666667E0_DP
 S_E%B_Y(8,104)=1.83333333333333E0_DP
 S_E%B_X(8,105)=0.416666666666667E0_DP
 S_E%A_Y(8,105)=0.361111111111111E0_DP
 S_E%B_X(8,118)=8.25000000000000E0_DP
 S_E%A_Y(8,118)=6.00000000000000E0_DP
 S_E%A_X(8,119)=3.25000000000000E0_DP
 S_E%B_Y(8,119)=-2.25000000000000E0_DP
 S_E%B_X(8,120)=-0.500000000000000E0_DP
 S_E%A_Y(8,120)=-0.500000000000000E0_DP
 S_E%A_X(8,133)=16.0000000000000E0_DP
 S_E%B_Y(8,133)=-8.50000000000000E0_DP
 S_E%B_X(8,134)=-9.00000000000000E0_DP
 S_E%A_Y(8,134)=-7.50000000000000E0_DP
 S_E%A_X(8,135)=-4.00000000000000E0_DP
 S_E%B_Y(8,135)=3.00000000000000E0_DP
 S_E%B_X(8,136)=1.00000000000000E0_DP
 S_E%A_Y(8,136)=1.00000000000000E0_DP
 S_E%B_X(8,149)=-14.8750000000000E0_DP
 S_E%A_Y(8,149)=-9.62500000000000E0_DP
 S_E%A_X(8,150)=-17.5000000000000E0_DP
 S_E%B_Y(8,150)=10.5000000000000E0_DP
 S_E%B_X(8,151)=10.5000000000000E0_DP
 S_E%A_Y(8,151)=10.5000000000000E0_DP
 S_E%A_X(8,152)=7.00000000000000E0_DP
 S_E%B_Y(8,152)=-7.00000000000000E0_DP
 S_E%A_X(8,166)=-11.5500000000000E0_DP
 S_E%B_Y(8,166)=4.20000000000001E0_DP
 S_E%B_X(8,167)=15.7500000000000E0_DP
 S_E%A_Y(8,167)=12.2500000000000E0_DP
 S_E%A_X(8,168)=21.0000000000000E0_DP
 S_E%B_Y(8,168)=-14.0000000000000E0_DP
 S_E%B_X(8,169)=-21.0000000000000E0_DP
 S_E%A_Y(8,169)=-21.0000000000000E0_DP
 S_E%B_X(8,184)=3.50000000000001E0_DP
 S_E%A_Y(8,184)=1.75000000000000E0_DP
 S_E%A_X(8,185)=12.2500000000000E0_DP
 S_E%B_Y(8,185)=-5.25000000000000E0_DP
 S_E%B_X(8,186)=-17.5000000000000E0_DP
 S_E%A_Y(8,186)=-17.5000000000000E0_DP
 S_E%A_X(8,187)=-35.0000000000000E0_DP
 S_E%B_Y(8,187)=35.0000000000000E0_DP
 S_E%A_X(8,203)=1.00000000000000E0_DP
 S_E%B_X(8,204)=-3.50000000000000E0_DP
 S_E%A_Y(8,204)=-2.33333333333333E0_DP
 S_E%A_X(8,205)=-14.0000000000000E0_DP
 S_E%B_Y(8,205)=7.00000000000000E0_DP
 S_E%B_X(8,206)=35.0000000000000E0_DP
 S_E%A_Y(8,206)=35.0000000000000E0_DP
 S_E%A_X(8,224)=-0.999999999999999E0_DP
 S_E%B_X(8,225)=3.50000000000000E0_DP
 S_E%A_Y(8,225)=3.50000000000000E0_DP
 S_E%A_X(8,226)=21.0000000000000E0_DP
 S_E%B_Y(8,226)=-21.0000000000000E0_DP
 S_E%A_X(8,246)=1.00000000000000E0_DP
 S_E%B_X(8,247)=-7.00000000000000E0_DP
 S_E%A_Y(8,247)=-7.00000000000000E0_DP
 S_E%A_X(8,269)=-1.00000000000000E0_DP
 S_E%B_Y(8,269)=1.00000000000000E0_DP
 S_E%B_X(9,91)=0.444444444444444E0_DP
 S_E%A_Y(9,91)=0.333333333333333E0_DP
 S_E%A_X(9,104)=3.33333333333333E0_DP
 S_E%B_Y(9,104)=-2.88888888888889E0_DP
 S_E%B_X(9,105)=-0.555555555555556E0_DP
 S_E%A_Y(9,105)=-0.444444444444444E0_DP
 S_E%B_X(9,118)=-13.0000000000000E0_DP
 S_E%A_Y(9,118)=-9.00000000000000E0_DP
 S_E%A_X(9,119)=-4.00000000000000E0_DP
 S_E%B_Y(9,119)=4.00000000000000E0_DP
 S_E%B_X(9,120)=1.00000000000000E0_DP
 S_E%A_Y(9,120)=1.00000000000000E0_DP
 S_E%A_X(9,133)=-24.0000000000000E0_DP
 S_E%B_Y(9,133)=20.0000000000000E0_DP
 S_E%B_X(9,134)=16.0000000000000E0_DP
 S_E%A_Y(9,134)=12.0000000000000E0_DP
 S_E%A_X(9,135)=8.00000000000000E0_DP
 S_E%B_Y(9,135)=-8.00000000000000E0_DP
 S_E%B_X(9,149)=35.0000000000000E0_DP
 S_E%A_Y(9,149)=21.0000000000000E0_DP
 S_E%A_X(9,150)=28.0000000000000E0_DP
 S_E%B_Y(9,150)=-28.0000000000000E0_DP
 S_E%B_X(9,151)=-28.0000000000000E0_DP
 S_E%A_Y(9,151)=-28.0000000000000E0_DP
 S_E%A_X(9,166)=25.2000000000000E0_DP
 S_E%B_Y(9,166)=-19.6000000000000E0_DP
 S_E%B_X(9,167)=-42.0000000000000E0_DP
 S_E%A_Y(9,167)=-28.0000000000000E0_DP
 S_E%A_X(9,168)=-56.0000000000000E0_DP
 S_E%B_Y(9,168)=56.0000000000000E0_DP
 S_E%B_X(9,184)=-16.3333333333333E0_DP
 S_E%A_Y(9,184)=-7.00000000000000E0_DP
 S_E%A_X(9,185)=-28.0000000000000E0_DP
 S_E%B_Y(9,185)=28.0000000000000E0_DP
 S_E%B_X(9,186)=70.0000000000000E0_DP
 S_E%A_Y(9,186)=70.0000000000000E0_DP
 S_E%A_X(9,203)=-4.00000000000000E0_DP
 S_E%B_Y(9,203)=2.66666666666666E0_DP
 S_E%B_X(9,204)=18.6666666666667E0_DP
 S_E%A_Y(9,204)=9.33333333333333E0_DP
 S_E%A_X(9,205)=56.0000000000000E0_DP
 S_E%B_Y(9,205)=-56.0000000000000E0_DP
 S_E%B_X(9,223)=0.999999999999999E0_DP
 S_E%A_X(9,224)=4.00000000000000E0_DP
 S_E%B_Y(9,224)=-4.00000000000000E0_DP
 S_E%B_X(9,225)=-28.0000000000000E0_DP
 S_E%A_Y(9,225)=-28.0000000000000E0_DP
 S_E%B_X(9,245)=-0.999999999999999E0_DP
 S_E%A_X(9,246)=-8.00000000000000E0_DP
 S_E%B_Y(9,246)=8.00000000000000E0_DP
 S_E%B_X(9,268)=1.00000000000000E0_DP
 S_E%A_Y(9,268)=1.00000000000000E0_DP
 S_E%B_X(10,91)=-0.500000000000000E0_DP
 S_E%A_Y(10,91)=-0.500000000000000E0_DP
 S_E%A_X(10,104)=-5.00000000000000E0_DP
 S_E%B_Y(10,104)=4.00000000000000E0_DP
 S_E%B_X(10,105)=1.00000000000000E0_DP
 S_E%A_Y(10,105)=1.00000000000000E0_DP
 S_E%B_X(10,118)=18.0000000000000E0_DP
 S_E%A_Y(10,118)=18.0000000000000E0_DP
 S_E%A_X(10,119)=9.00000000000000E0_DP
 S_E%B_Y(10,119)=-9.00000000000000E0_DP
 S_E%A_X(10,133)=48.0000000000000E0_DP
 S_E%B_Y(10,133)=-36.0000000000000E0_DP
 S_E%B_X(10,134)=-36.0000000000000E0_DP
 S_E%A_Y(10,134)=-36.0000000000000E0_DP
 S_E%B_X(10,149)=-63.0000000000000E0_DP
 S_E%A_Y(10,149)=-63.0000000000000E0_DP
 S_E%A_X(10,150)=-84.0000000000000E0_DP
 S_E%B_Y(10,150)=84.0000000000000E0_DP
 S_E%A_X(10,166)=-75.6000000000000E0_DP
 S_E%B_Y(10,166)=50.4000000000000E0_DP
 S_E%B_X(10,167)=126.000000000000E0_DP
 S_E%A_Y(10,167)=126.000000000000E0_DP
 S_E%B_X(10,184)=42.0000000000000E0_DP
 S_E%A_Y(10,184)=42.0000000000000E0_DP
 S_E%A_X(10,185)=126.000000000000E0_DP
 S_E%B_Y(10,185)=-126.000000000000E0_DP
 S_E%A_X(10,203)=24.0000000000000E0_DP
 S_E%B_Y(10,203)=-12.0000000000000E0_DP
 S_E%B_X(10,204)=-84.0000000000000E0_DP
 S_E%A_Y(10,204)=-84.0000000000000E0_DP
 S_E%B_X(10,223)=-4.50000000000000E0_DP
 S_E%A_Y(10,223)=-4.50000000000000E0_DP
 S_E%A_X(10,224)=-36.0000000000000E0_DP
 S_E%B_Y(10,224)=36.0000000000000E0_DP
 S_E%A_X(10,244)=-1.00000000000000E0_DP
 S_E%B_X(10,245)=9.00000000000000E0_DP
 S_E%A_Y(10,245)=9.00000000000000E0_DP
 S_E%A_X(10,267)=1.00000000000000E0_DP
 S_E%B_Y(10,267)=-1.00000000000000E0_DP
 S_E%B_X(11,91)=1.00000000000000E0_DP
 S_E%A_Y(11,91)=1.00000000000000E0_DP
 S_E%A_X(11,104)=10.0000000000000E0_DP
 S_E%B_Y(11,104)=-10.0000000000000E0_DP
 S_E%B_X(11,118)=-45.0000000000000E0_DP
 S_E%A_Y(11,118)=-45.0000000000000E0_DP
 S_E%A_X(11,133)=-120.000000000000E0_DP
 S_E%B_Y(11,133)=120.000000000000E0_DP
 S_E%B_X(11,149)=210.000000000000E0_DP
 S_E%A_Y(11,149)=210.000000000000E0_DP
 S_E%A_X(11,166)=252.000000000000E0_DP
 S_E%B_Y(11,166)=-252.000000000000E0_DP
 S_E%B_X(11,184)=-210.000000000000E0_DP
 S_E%A_Y(11,184)=-210.000000000000E0_DP
 S_E%A_X(11,203)=-120.000000000000E0_DP
 S_E%B_Y(11,203)=120.000000000000E0_DP
 S_E%B_X(11,223)=45.0000000000000E0_DP
 S_E%A_Y(11,223)=45.0000000000000E0_DP
 S_E%A_X(11,244)=10.0000000000000E0_DP
 S_E%B_Y(11,244)=-10.0000000000000E0_DP
 S_E%B_X(11,266)=-1.00000000000000E0_DP
 S_E%A_Y(11,266)=-1.00000000000000E0_DP
 S_E%B_X(12,78)=1.00000000000000E0_DP
 S_E%A_Y(12,78)=1.00000000000000E0_DP
 S_E%A_X(12,90)=11.0000000000000E0_DP
 S_E%B_Y(12,90)=-11.0000000000000E0_DP
 S_E%B_X(12,103)=-55.0000000000000E0_DP
 S_E%A_Y(12,103)=-55.0000000000000E0_DP
 S_E%A_X(12,117)=-165.000000000000E0_DP
 S_E%B_Y(12,117)=165.000000000000E0_DP
 S_E%B_X(12,132)=330.000000000000E0_DP
 S_E%A_Y(12,132)=330.000000000000E0_DP
 S_E%A_X(12,148)=462.000000000000E0_DP
 S_E%B_Y(12,148)=-462.000000000000E0_DP
 S_E%B_X(12,165)=-462.000000000000E0_DP
 S_E%A_Y(12,165)=-462.000000000000E0_DP
 S_E%A_X(12,183)=-330.000000000000E0_DP
 S_E%B_Y(12,183)=330.000000000000E0_DP
 S_E%B_X(12,202)=165.000000000000E0_DP
 S_E%A_Y(12,202)=165.000000000000E0_DP
 S_E%A_X(12,222)=55.0000000000000E0_DP
 S_E%B_Y(12,222)=-55.0000000000000E0_DP
 S_E%B_X(12,243)=-11.0000000000000E0_DP
 S_E%A_Y(12,243)=-11.0000000000000E0_DP
 S_E%A_X(12,265)=-1.00000000000000E0_DP
 S_E%B_Y(12,265)=1.00000000000000E0_DP
 S_E%B_X(13,66)=1.00000000000000E0_DP
 S_E%A_Y(13,66)=1.00000000000000E0_DP
 S_E%A_X(13,77)=12.0000000000000E0_DP
 S_E%B_Y(13,77)=-12.0000000000000E0_DP
 S_E%B_X(13,89)=-66.0000000000000E0_DP
 S_E%A_Y(13,89)=-66.0000000000000E0_DP
 S_E%A_X(13,102)=-220.000000000000E0_DP
 S_E%B_Y(13,102)=220.000000000000E0_DP
 S_E%B_X(13,116)=495.000000000000E0_DP
 S_E%A_Y(13,116)=495.000000000000E0_DP
 S_E%A_X(13,131)=792.000000000000E0_DP
 S_E%B_Y(13,131)=-792.000000000000E0_DP
 S_E%B_X(13,147)=-924.000000000000E0_DP
 S_E%A_Y(13,147)=-924.000000000000E0_DP
 S_E%A_X(13,164)=-792.000000000000E0_DP
 S_E%B_Y(13,164)=792.000000000000E0_DP
 S_E%B_X(13,182)=495.000000000000E0_DP
 S_E%A_Y(13,182)=495.000000000000E0_DP
 S_E%A_X(13,201)=220.000000000000E0_DP
 S_E%B_Y(13,201)=-220.000000000000E0_DP
 S_E%B_X(13,221)=-66.0000000000000E0_DP
 S_E%A_Y(13,221)=-66.0000000000000E0_DP
 S_E%A_X(13,242)=-12.0000000000000E0_DP
 S_E%B_Y(13,242)=12.0000000000000E0_DP
 S_E%B_X(13,264)=1.00000000000000E0_DP
 S_E%A_Y(13,264)=1.00000000000000E0_DP
 S_E%B_X(14,55)=1.00000000000000E0_DP
 S_E%A_Y(14,55)=1.00000000000000E0_DP
 S_E%A_X(14,65)=13.0000000000000E0_DP
 S_E%B_Y(14,65)=-13.0000000000000E0_DP
 S_E%B_X(14,76)=-78.0000000000000E0_DP
 S_E%A_Y(14,76)=-78.0000000000000E0_DP
 S_E%A_X(14,88)=-286.000000000000E0_DP
 S_E%B_Y(14,88)=286.000000000000E0_DP
 S_E%B_X(14,101)=715.000000000000E0_DP
 S_E%A_Y(14,101)=715.000000000000E0_DP
 S_E%A_X(14,115)=1287.00000000000E0_DP
 S_E%B_Y(14,115)=-1287.00000000000E0_DP
 S_E%B_X(14,130)=-1716.00000000000E0_DP
 S_E%A_Y(14,130)=-1716.00000000000E0_DP
 S_E%A_X(14,146)=-1716.00000000000E0_DP
 S_E%B_Y(14,146)=1716.00000000000E0_DP
 S_E%B_X(14,163)=1287.00000000000E0_DP
 S_E%A_Y(14,163)=1287.00000000000E0_DP
 S_E%A_X(14,181)=715.000000000000E0_DP
 S_E%B_Y(14,181)=-715.000000000000E0_DP
 S_E%B_X(14,200)=-286.000000000000E0_DP
 S_E%A_Y(14,200)=-286.000000000000E0_DP
 S_E%A_X(14,220)=-78.0000000000000E0_DP
 S_E%B_Y(14,220)=78.0000000000000E0_DP
 S_E%B_X(14,241)=13.0000000000000E0_DP
 S_E%A_Y(14,241)=13.0000000000000E0_DP
 S_E%A_X(14,263)=1.00000000000000E0_DP
 S_E%B_Y(14,263)=-1.00000000000000E0_DP
 S_E%B_X(15,45)=1.00000000000000E0_DP
 S_E%A_Y(15,45)=1.00000000000000E0_DP
 S_E%A_X(15,54)=14.0000000000000E0_DP
 S_E%B_Y(15,54)=-14.0000000000000E0_DP
 S_E%B_X(15,64)=-91.0000000000000E0_DP
 S_E%A_Y(15,64)=-91.0000000000000E0_DP
 S_E%A_X(15,75)=-364.000000000000E0_DP
 S_E%B_Y(15,75)=364.000000000000E0_DP
 S_E%B_X(15,87)=1001.00000000000E0_DP
 S_E%A_Y(15,87)=1001.00000000000E0_DP
 S_E%A_X(15,100)=2002.00000000000E0_DP
 S_E%B_Y(15,100)=-2002.00000000000E0_DP
 S_E%B_X(15,114)=-3003.00000000000E0_DP
 S_E%A_Y(15,114)=-3003.00000000000E0_DP
 S_E%A_X(15,129)=-3432.00000000000E0_DP
 S_E%B_Y(15,129)=3432.00000000000E0_DP
 S_E%B_X(15,145)=3003.00000000000E0_DP
 S_E%A_Y(15,145)=3003.00000000000E0_DP
 S_E%A_X(15,162)=2002.00000000000E0_DP
 S_E%B_Y(15,162)=-2002.00000000000E0_DP
 S_E%B_X(15,180)=-1001.00000000000E0_DP
 S_E%A_Y(15,180)=-1001.00000000000E0_DP
 S_E%A_X(15,199)=-364.000000000000E0_DP
 S_E%B_Y(15,199)=364.000000000000E0_DP
 S_E%B_X(15,219)=91.0000000000000E0_DP
 S_E%A_Y(15,219)=91.0000000000000E0_DP
 S_E%A_X(15,240)=14.0000000000000E0_DP
 S_E%B_Y(15,240)=-14.0000000000000E0_DP
 S_E%B_X(15,262)=-1.00000000000000E0_DP
 S_E%A_Y(15,262)=-1.00000000000000E0_DP
 S_E%B_X(16,36)=1.00000000000000E0_DP
 S_E%A_Y(16,36)=1.00000000000000E0_DP
 S_E%A_X(16,44)=15.0000000000000E0_DP
 S_E%B_Y(16,44)=-15.0000000000000E0_DP
 S_E%B_X(16,53)=-105.000000000000E0_DP
 S_E%A_Y(16,53)=-105.000000000000E0_DP
 S_E%A_X(16,63)=-455.000000000000E0_DP
 S_E%B_Y(16,63)=455.000000000000E0_DP
 S_E%B_X(16,74)=1365.00000000000E0_DP
 S_E%A_Y(16,74)=1365.00000000000E0_DP
 S_E%A_X(16,86)=3003.00000000000E0_DP
 S_E%B_Y(16,86)=-3003.00000000000E0_DP
 S_E%B_X(16,99)=-5005.00000000000E0_DP
 S_E%A_Y(16,99)=-5005.00000000000E0_DP
 S_E%A_X(16,113)=-6435.00000000000E0_DP
 S_E%B_Y(16,113)=6435.00000000000E0_DP
 S_E%B_X(16,128)=6435.00000000000E0_DP
 S_E%A_Y(16,128)=6435.00000000000E0_DP
 S_E%A_X(16,144)=5005.00000000000E0_DP
 S_E%B_Y(16,144)=-5005.00000000000E0_DP
 S_E%B_X(16,161)=-3003.00000000000E0_DP
 S_E%A_Y(16,161)=-3003.00000000000E0_DP
 S_E%A_X(16,179)=-1365.00000000000E0_DP
 S_E%B_Y(16,179)=1365.00000000000E0_DP
 S_E%B_X(16,198)=455.000000000000E0_DP
 S_E%A_Y(16,198)=455.000000000000E0_DP
 S_E%A_X(16,218)=105.000000000000E0_DP
 S_E%B_Y(16,218)=-105.000000000000E0_DP
 S_E%B_X(16,239)=-15.0000000000000E0_DP
 S_E%A_Y(16,239)=-15.0000000000000E0_DP
 S_E%A_X(16,261)=-1.00000000000000E0_DP
 S_E%B_Y(16,261)=1.00000000000000E0_DP
 S_E%B_X(17,28)=1.00000000000000E0_DP
 S_E%A_Y(17,28)=1.00000000000000E0_DP
 S_E%A_X(17,35)=16.0000000000000E0_DP
 S_E%B_Y(17,35)=-16.0000000000000E0_DP
 S_E%B_X(17,43)=-120.000000000000E0_DP
 S_E%A_Y(17,43)=-120.000000000000E0_DP
 S_E%A_X(17,52)=-560.000000000000E0_DP
 S_E%B_Y(17,52)=560.000000000000E0_DP
 S_E%B_X(17,62)=1820.00000000000E0_DP
 S_E%A_Y(17,62)=1820.00000000000E0_DP
 S_E%A_X(17,73)=4368.00000000000E0_DP
 S_E%B_Y(17,73)=-4368.00000000000E0_DP
 S_E%B_X(17,85)=-8008.00000000000E0_DP
 S_E%A_Y(17,85)=-8008.00000000000E0_DP
 S_E%A_X(17,98)=-11440.0000000000E0_DP
 S_E%B_Y(17,98)=11440.0000000000E0_DP
 S_E%B_X(17,112)=12870.0000000000E0_DP
 S_E%A_Y(17,112)=12870.0000000000E0_DP
 S_E%A_X(17,127)=11440.0000000000E0_DP
 S_E%B_Y(17,127)=-11440.0000000000E0_DP
 S_E%B_X(17,143)=-8008.00000000000E0_DP
 S_E%A_Y(17,143)=-8008.00000000000E0_DP
 S_E%A_X(17,160)=-4368.00000000000E0_DP
 S_E%B_Y(17,160)=4368.00000000000E0_DP
 S_E%B_X(17,178)=1820.00000000000E0_DP
 S_E%A_Y(17,178)=1820.00000000000E0_DP
 S_E%A_X(17,197)=560.000000000000E0_DP
 S_E%B_Y(17,197)=-560.000000000000E0_DP
 S_E%B_X(17,217)=-120.000000000000E0_DP
 S_E%A_Y(17,217)=-120.000000000000E0_DP
 S_E%A_X(17,238)=-16.0000000000000E0_DP
 S_E%B_Y(17,238)=16.0000000000000E0_DP
 S_E%B_X(17,260)=1.00000000000000E0_DP
 S_E%A_Y(17,260)=1.00000000000000E0_DP
 S_E%B_X(18,21)=1.00000000000000E0_DP
 S_E%A_Y(18,21)=1.00000000000000E0_DP
 S_E%A_X(18,27)=17.0000000000000E0_DP
 S_E%B_Y(18,27)=-17.0000000000000E0_DP
 S_E%B_X(18,34)=-136.000000000000E0_DP
 S_E%A_Y(18,34)=-136.000000000000E0_DP
 S_E%A_X(18,42)=-680.000000000000E0_DP
 S_E%B_Y(18,42)=680.000000000000E0_DP
 S_E%B_X(18,51)=2380.00000000000E0_DP
 S_E%A_Y(18,51)=2380.00000000000E0_DP
 S_E%A_X(18,61)=6188.00000000000E0_DP
 S_E%B_Y(18,61)=-6188.00000000000E0_DP
 S_E%B_X(18,72)=-12376.0000000000E0_DP
 S_E%A_Y(18,72)=-12376.0000000000E0_DP
 S_E%A_X(18,84)=-19448.0000000000E0_DP
 S_E%B_Y(18,84)=19448.0000000000E0_DP
 S_E%B_X(18,97)=24310.0000000000E0_DP
 S_E%A_Y(18,97)=24310.0000000000E0_DP
 S_E%A_X(18,111)=24310.0000000000E0_DP
 S_E%B_Y(18,111)=-24310.0000000000E0_DP
 S_E%B_X(18,126)=-19448.0000000000E0_DP
 S_E%A_Y(18,126)=-19448.0000000000E0_DP
 S_E%A_X(18,142)=-12376.0000000000E0_DP
 S_E%B_Y(18,142)=12376.0000000000E0_DP
 S_E%B_X(18,159)=6188.00000000000E0_DP
 S_E%A_Y(18,159)=6188.00000000000E0_DP
 S_E%A_X(18,177)=2380.00000000000E0_DP
 S_E%B_Y(18,177)=-2380.00000000000E0_DP
 S_E%B_X(18,196)=-680.000000000000E0_DP
 S_E%A_Y(18,196)=-680.000000000000E0_DP
 S_E%A_X(18,216)=-136.000000000000E0_DP
 S_E%B_Y(18,216)=136.000000000000E0_DP
 S_E%B_X(18,237)=17.0000000000000E0_DP
 S_E%A_Y(18,237)=17.0000000000000E0_DP
 S_E%A_X(18,259)=1.00000000000000E0_DP
 S_E%B_Y(18,259)=-1.00000000000000E0_DP
 S_E%B_X(19,15)=1.00000000000000E0_DP
 S_E%A_Y(19,15)=1.00000000000000E0_DP
 S_E%A_X(19,20)=18.0000000000000E0_DP
 S_E%B_Y(19,20)=-18.0000000000000E0_DP
 S_E%B_X(19,26)=-153.000000000000E0_DP
 S_E%A_Y(19,26)=-153.000000000000E0_DP
 S_E%A_X(19,33)=-816.000000000000E0_DP
 S_E%B_Y(19,33)=816.000000000000E0_DP
 S_E%B_X(19,41)=3060.00000000000E0_DP
 S_E%A_Y(19,41)=3060.00000000000E0_DP
 S_E%A_X(19,50)=8568.00000000000E0_DP
 S_E%B_Y(19,50)=-8568.00000000000E0_DP
 S_E%B_X(19,60)=-18564.0000000000E0_DP
 S_E%A_Y(19,60)=-18564.0000000000E0_DP
 S_E%A_X(19,71)=-31824.0000000000E0_DP
 S_E%B_Y(19,71)=31824.0000000000E0_DP
 S_E%B_X(19,83)=43758.0000000000E0_DP
 S_E%A_Y(19,83)=43758.0000000000E0_DP
 S_E%A_X(19,96)=48620.0000000000E0_DP
 S_E%B_Y(19,96)=-48620.0000000000E0_DP
 S_E%B_X(19,110)=-43758.0000000000E0_DP
 S_E%A_Y(19,110)=-43758.0000000000E0_DP
 S_E%A_X(19,125)=-31824.0000000000E0_DP
 S_E%B_Y(19,125)=31824.0000000000E0_DP
 S_E%B_X(19,141)=18564.0000000000E0_DP
 S_E%A_Y(19,141)=18564.0000000000E0_DP
 S_E%A_X(19,158)=8568.00000000000E0_DP
 S_E%B_Y(19,158)=-8568.00000000000E0_DP
 S_E%B_X(19,176)=-3060.00000000000E0_DP
 S_E%A_Y(19,176)=-3060.00000000000E0_DP
 S_E%A_X(19,195)=-816.000000000000E0_DP
 S_E%B_Y(19,195)=816.000000000000E0_DP
 S_E%B_X(19,215)=153.000000000000E0_DP
 S_E%A_Y(19,215)=153.000000000000E0_DP
 S_E%A_X(19,236)=18.0000000000000E0_DP
 S_E%B_Y(19,236)=-18.0000000000000E0_DP
 S_E%B_X(19,258)=-1.00000000000000E0_DP
 S_E%A_Y(19,258)=-1.00000000000000E0_DP
 S_E%B_X(20,10)=1.00000000000000E0_DP
 S_E%A_Y(20,10)=1.00000000000000E0_DP
 S_E%A_X(20,14)=19.0000000000000E0_DP
 S_E%B_Y(20,14)=-19.0000000000000E0_DP
 S_E%B_X(20,19)=-171.000000000000E0_DP
 S_E%A_Y(20,19)=-171.000000000000E0_DP
 S_E%A_X(20,25)=-969.000000000000E0_DP
 S_E%B_Y(20,25)=969.000000000000E0_DP
 S_E%B_X(20,32)=3876.00000000000E0_DP
 S_E%A_Y(20,32)=3876.00000000000E0_DP
 S_E%A_X(20,40)=11628.0000000000E0_DP
 S_E%B_Y(20,40)=-11628.0000000000E0_DP
 S_E%B_X(20,49)=-27132.0000000000E0_DP
 S_E%A_Y(20,49)=-27132.0000000000E0_DP
 S_E%A_X(20,59)=-50388.0000000000E0_DP
 S_E%B_Y(20,59)=50388.0000000000E0_DP
 S_E%B_X(20,70)=75582.0000000000E0_DP
 S_E%A_Y(20,70)=75582.0000000000E0_DP
 S_E%A_X(20,82)=92378.0000000000E0_DP
 S_E%B_Y(20,82)=-92378.0000000000E0_DP
 S_E%B_X(20,95)=-92378.0000000000E0_DP
 S_E%A_Y(20,95)=-92378.0000000000E0_DP
 S_E%A_X(20,109)=-75582.0000000000E0_DP
 S_E%B_Y(20,109)=75582.0000000000E0_DP
 S_E%B_X(20,124)=50388.0000000000E0_DP
 S_E%A_Y(20,124)=50388.0000000000E0_DP
 S_E%A_X(20,140)=27132.0000000000E0_DP
 S_E%B_Y(20,140)=-27132.0000000000E0_DP
 S_E%B_X(20,157)=-11628.0000000000E0_DP
 S_E%A_Y(20,157)=-11628.0000000000E0_DP
 S_E%A_X(20,175)=-3876.00000000000E0_DP
 S_E%B_Y(20,175)=3876.00000000000E0_DP
 S_E%B_X(20,194)=969.000000000000E0_DP
 S_E%A_Y(20,194)=969.000000000000E0_DP
 S_E%A_X(20,214)=171.000000000000E0_DP
 S_E%B_Y(20,214)=-171.000000000000E0_DP
 S_E%B_X(20,235)=-19.0000000000000E0_DP
 S_E%A_Y(20,235)=-19.0000000000000E0_DP
 S_E%A_X(20,257)=-1.00000000000000E0_DP
 S_E%B_Y(20,257)=1.00000000000000E0_DP
 S_E%B_X(21,6)=1.00000000000000E0_DP
 S_E%A_Y(21,6)=1.00000000000000E0_DP
 S_E%A_X(21,9)=20.0000000000000E0_DP
 S_E%B_Y(21,9)=-20.0000000000000E0_DP
 S_E%B_X(21,13)=-190.000000000000E0_DP
 S_E%A_Y(21,13)=-190.000000000000E0_DP
 S_E%A_X(21,18)=-1140.00000000000E0_DP
 S_E%B_Y(21,18)=1140.00000000000E0_DP
 S_E%B_X(21,24)=4845.00000000000E0_DP
 S_E%A_Y(21,24)=4845.00000000000E0_DP
 S_E%A_X(21,31)=15504.0000000000E0_DP
 S_E%B_Y(21,31)=-15504.0000000000E0_DP
 S_E%B_X(21,39)=-38760.0000000000E0_DP
 S_E%A_Y(21,39)=-38760.0000000000E0_DP
 S_E%A_X(21,48)=-77520.0000000000E0_DP
 S_E%B_Y(21,48)=77520.0000000000E0_DP
 S_E%B_X(21,58)=125970.000000000E0_DP
 S_E%A_Y(21,58)=125970.000000000E0_DP
 S_E%A_X(21,69)=167960.000000000E0_DP
 S_E%B_Y(21,69)=-167960.000000000E0_DP
 S_E%B_X(21,81)=-184756.000000000E0_DP
 S_E%A_Y(21,81)=-184756.000000000E0_DP
 S_E%A_X(21,94)=-167960.000000000E0_DP
 S_E%B_Y(21,94)=167960.000000000E0_DP
 S_E%B_X(21,108)=125970.000000000E0_DP
 S_E%A_Y(21,108)=125970.000000000E0_DP
 S_E%A_X(21,123)=77520.0000000000E0_DP
 S_E%B_Y(21,123)=-77520.0000000000E0_DP
 S_E%B_X(21,139)=-38760.0000000000E0_DP
 S_E%A_Y(21,139)=-38760.0000000000E0_DP
 S_E%A_X(21,156)=-15504.0000000000E0_DP
 S_E%B_Y(21,156)=15504.0000000000E0_DP
 S_E%B_X(21,174)=4845.00000000000E0_DP
 S_E%A_Y(21,174)=4845.00000000000E0_DP
 S_E%A_X(21,193)=1140.00000000000E0_DP
 S_E%B_Y(21,193)=-1140.00000000000E0_DP
 S_E%B_X(21,213)=-190.000000000000E0_DP
 S_E%A_Y(21,213)=-190.000000000000E0_DP
 S_E%A_X(21,234)=-20.0000000000000E0_DP
 S_E%B_Y(21,234)=20.0000000000000E0_DP
 S_E%B_X(21,256)=1.00000000000000E0_DP
 S_E%A_Y(21,256)=1.00000000000000E0_DP
 S_E%B_X(22,3)=1.00000000000000E0_DP
 S_E%A_Y(22,3)=1.00000000000000E0_DP
 S_E%A_X(22,5)=21.0000000000000E0_DP
 S_E%B_Y(22,5)=-21.0000000000000E0_DP
 S_E%B_X(22,8)=-210.000000000000E0_DP
 S_E%A_Y(22,8)=-210.000000000000E0_DP
 S_E%A_X(22,12)=-1330.00000000000E0_DP
 S_E%B_Y(22,12)=1330.00000000000E0_DP
 S_E%B_X(22,17)=5985.00000000000E0_DP
 S_E%A_Y(22,17)=5985.00000000000E0_DP
 S_E%A_X(22,23)=20349.0000000000E0_DP
 S_E%B_Y(22,23)=-20349.0000000000E0_DP
 S_E%B_X(22,30)=-54264.0000000000E0_DP
 S_E%A_Y(22,30)=-54264.0000000000E0_DP
 S_E%A_X(22,38)=-116280.000000000E0_DP
 S_E%B_Y(22,38)=116280.000000000E0_DP
 S_E%B_X(22,47)=203490.000000000E0_DP
 S_E%A_Y(22,47)=203490.000000000E0_DP
 S_E%A_X(22,57)=293930.000000000E0_DP
 S_E%B_Y(22,57)=-293930.000000000E0_DP
 S_E%B_X(22,68)=-352716.000000000E0_DP
 S_E%A_Y(22,68)=-352716.000000000E0_DP
 S_E%A_X(22,80)=-352716.000000000E0_DP
 S_E%B_Y(22,80)=352716.000000000E0_DP
 S_E%B_X(22,93)=293930.000000000E0_DP
 S_E%A_Y(22,93)=293930.000000000E0_DP
 S_E%A_X(22,107)=203490.000000000E0_DP
 S_E%B_Y(22,107)=-203490.000000000E0_DP
 S_E%B_X(22,122)=-116280.000000000E0_DP
 S_E%A_Y(22,122)=-116280.000000000E0_DP
 S_E%A_X(22,138)=-54264.0000000000E0_DP
 S_E%B_Y(22,138)=54264.0000000000E0_DP
 S_E%B_X(22,155)=20349.0000000000E0_DP
 S_E%A_Y(22,155)=20349.0000000000E0_DP
 S_E%A_X(22,173)=5985.00000000000E0_DP
 S_E%B_Y(22,173)=-5985.00000000000E0_DP
 S_E%B_X(22,192)=-1330.00000000000E0_DP
 S_E%A_Y(22,192)=-1330.00000000000E0_DP
 S_E%A_X(22,212)=-210.000000000000E0_DP
 S_E%B_Y(22,212)=210.000000000000E0_DP
 S_E%B_X(22,233)=21.0000000000000E0_DP
 S_E%A_Y(22,233)=21.0000000000000E0_DP
 S_E%A_X(22,255)=1.00000000000000E0_DP
 S_E%B_Y(22,255)=-1.00000000000000E0_DP
 S_E%VB(1,78)=-9.090909090909090E-002_DP
 S_E%VB(1,91)=0.100000000000000E0_DP
 S_E%VB(1,105)=-0.111111111111111E0_DP
 S_E%VB(1,120)=0.125000000000000E0_DP
 S_E%VB(1,136)=-0.142857142857143E0_DP
 S_E%VB(1,153)=0.166666666666667E0_DP
 S_E%VB(1,171)=-0.200000000000000E0_DP
 S_E%VB(1,190)=0.250000000000000E0_DP
 S_E%VB(1,210)=-0.333333333333333E0_DP
 S_E%VB(1,231)=0.500000000000000E0_DP
 S_E%VB(1,253)=-1.00000000000000E0_DP
 S_E%VA(1,275)=-1.00000000000000E0_DP
 S_E%VB(2,78)=4.545454545454545E-002_DP
 S_E%VA(2,90)=0.100000000000000E0_DP
 S_E%VB(2,91)=-5.000000000000000E-002_DP
 S_E%VA(2,104)=-0.111111111111111E0_DP
 S_E%VB(2,105)=5.555555555555555E-002_DP
 S_E%VA(2,119)=0.125000000000000E0_DP
 S_E%VB(2,120)=-6.249999999999999E-002_DP
 S_E%VA(2,135)=-0.142857142857143E0_DP
 S_E%VB(2,136)=7.142857142857141E-002_DP
 S_E%VA(2,152)=0.166666666666667E0_DP
 S_E%VB(2,153)=-8.333333333333334E-002_DP
 S_E%VA(2,170)=-0.200000000000000E0_DP
 S_E%VB(2,171)=0.100000000000000E0_DP
 S_E%VA(2,189)=0.250000000000000E0_DP
 S_E%VB(2,190)=-0.125000000000000E0_DP
 S_E%VA(2,209)=-0.333333333333333E0_DP
 S_E%VB(2,210)=0.166666666666667E0_DP
 S_E%VA(2,230)=0.500000000000000E0_DP
 S_E%VB(2,231)=-0.500000000000000E0_DP
 S_E%VA(2,252)=-1.00000000000000E0_DP
 S_E%VB(2,274)=0.500000000000000E0_DP
 S_E%VB(3,78)=-4.646464646464647E-002_DP
 S_E%VA(3,90)=-0.100000000000000E0_DP
 S_E%VB(3,91)=5.138888888888889E-002_DP
 S_E%VB(3,103)=0.111111111111111E0_DP
 S_E%VA(3,104)=0.111111111111111E0_DP
 S_E%VB(3,105)=-5.753968253968253E-002_DP
 S_E%VB(3,118)=-0.125000000000000E0_DP
 S_E%VA(3,119)=-0.125000000000000E0_DP
 S_E%VB(3,120)=6.547619047619047E-002_DP
 S_E%VB(3,134)=0.142857142857143E0_DP
 S_E%VA(3,135)=0.142857142857143E0_DP
 S_E%VB(3,136)=-7.619047619047617E-002_DP
 S_E%VB(3,151)=-0.166666666666667E0_DP
 S_E%VA(3,152)=-0.166666666666667E0_DP
 S_E%VB(3,153)=9.166666666666667E-002_DP
 S_E%VB(3,169)=0.200000000000000E0_DP
 S_E%VA(3,170)=0.200000000000000E0_DP
 S_E%VB(3,171)=-0.116666666666667E0_DP
 S_E%VB(3,188)=-0.250000000000000E0_DP
 S_E%VA(3,189)=-0.250000000000000E0_DP
 S_E%VB(3,190)=0.166666666666667E0_DP
 S_E%VB(3,208)=0.333333333333333E0_DP
 S_E%VA(3,209)=0.333333333333333E0_DP
 S_E%VB(3,210)=-0.333333333333333E0_DP
 S_E%VB(3,229)=-0.500000000000000E0_DP
 S_E%VA(3,230)=-1.00000000000000E0_DP
 S_E%VB(3,251)=1.00000000000000E0_DP
 S_E%VA(3,273)=0.333333333333333E0_DP
 S_E%VB(4,78)=3.560606060606061E-002_DP
 S_E%VA(4,90)=0.154166666666667E0_DP
 S_E%VB(4,91)=-3.958333333333335E-002_DP
 S_E%VB(4,103)=-0.166666666666667E0_DP
 S_E%VA(4,104)=-0.172619047619048E0_DP
 S_E%VB(4,105)=4.464285714285714E-002_DP
 S_E%VA(4,117)=-0.125000000000000E0_DP
 S_E%VB(4,118)=0.187500000000000E0_DP
 S_E%VA(4,119)=0.196428571428571E0_DP
 S_E%VB(4,120)=-5.133928571428571E-002_DP
 S_E%VA(4,133)=0.142857142857143E0_DP
 S_E%VB(4,134)=-0.214285714285714E0_DP
 S_E%VA(4,135)=-0.228571428571429E0_DP
 S_E%VB(4,136)=6.071428571428570E-002_DP
 S_E%VA(4,150)=-0.166666666666667E0_DP
 S_E%VB(4,151)=0.250000000000000E0_DP
 S_E%VA(4,152)=0.275000000000000E0_DP
 S_E%VB(4,153)=-7.500000000000000E-002_DP
 S_E%VA(4,168)=0.200000000000000E0_DP
 S_E%VB(4,169)=-0.300000000000000E0_DP
 S_E%VA(4,170)=-0.350000000000000E0_DP
 S_E%VB(4,171)=0.100000000000000E0_DP
 S_E%VA(4,187)=-0.250000000000000E0_DP
 S_E%VB(4,188)=0.375000000000000E0_DP
 S_E%VA(4,189)=0.500000000000000E0_DP
 S_E%VB(4,190)=-0.250000000000000E0_DP
 S_E%VA(4,207)=0.333333333333333E0_DP
 S_E%VB(4,208)=-0.500000000000000E0_DP
 S_E%VA(4,209)=-1.00000000000000E0_DP
 S_E%VA(4,228)=-0.500000000000000E0_DP
 S_E%VB(4,229)=1.50000000000000E0_DP
 S_E%VA(4,250)=1.00000000000000E0_DP
 S_E%VB(4,272)=-0.250000000000000E0_DP
 S_E%VB(5,78)=-3.728354978354979E-002_DP
 S_E%VA(5,90)=-0.158333333333333E0_DP
 S_E%VB(5,91)=4.196428571428571E-002_DP
 S_E%VB(5,103)=0.345238095238095E0_DP
 S_E%VA(5,104)=0.178571428571429E0_DP
 S_E%VB(5,105)=-4.821428571428571E-002_DP
 S_E%VA(5,117)=0.250000000000000E0_DP
 S_E%VB(5,118)=-0.392857142857143E0_DP
 S_E%VA(5,119)=-0.205357142857143E0_DP
 S_E%VB(5,120)=5.714285714285714E-002_DP
 S_E%VB(5,132)=-0.142857142857144E0_DP
 S_E%VA(5,133)=-0.285714285714285E0_DP
 S_E%VB(5,134)=0.457142857142857E0_DP
 S_E%VA(5,135)=0.242857142857143E0_DP
 S_E%VB(5,136)=-7.142857142857142E-002_DP
 S_E%VB(5,149)=0.166666666666664E0_DP
 S_E%VA(5,150)=0.333333333333333E0_DP
 S_E%VB(5,151)=-0.550000000000000E0_DP
 S_E%VA(5,152)=-0.300000000000000E0_DP
 S_E%VB(5,153)=0.100000000000000E0_DP
 S_E%VB(5,167)=-0.200000000000000E0_DP
 S_E%VA(5,168)=-0.400000000000000E0_DP
 S_E%VB(5,169)=0.700000000000000E0_DP
 S_E%VA(5,170)=0.400000000000000E0_DP
 S_E%VB(5,171)=-0.200000000000000E0_DP
 S_E%VB(5,186)=0.250000000000000E0_DP
 S_E%VA(5,187)=0.500000000000000E0_DP
 S_E%VB(5,188)=-1.00000000000000E0_DP
 S_E%VA(5,189)=-1.00000000000000E0_DP
 S_E%VB(5,206)=-0.333333333333333E0_DP
 S_E%VA(5,207)=-0.666666666666667E0_DP
 S_E%VB(5,208)=2.00000000000000E0_DP
 S_E%VB(5,227)=0.500000000000000E0_DP
 S_E%VA(5,228)=2.00000000000000E0_DP
 S_E%VB(5,249)=-1.00000000000000E0_DP
 S_E%VA(5,271)=-0.200000000000000E0_DP
 S_E%VB(6,78)=3.260281385281385E-002_DP
 S_E%VA(6,90)=0.209821428571429E0_DP
 S_E%VB(6,91)=-3.720238095238095E-002_DP
 S_E%VB(6,103)=-0.446428571428571E0_DP
 S_E%VA(6,104)=-0.241071428571429E0_DP
 S_E%VB(6,105)=4.365079365079364E-002_DP
 S_E%VA(6,117)=-0.654761904761905E0_DP
 S_E%VB(6,118)=0.513392857142857E0_DP
 S_E%VA(6,119)=0.285714285714286E0_DP
 S_E%VB(6,120)=-5.357142857142856E-002_DP
 S_E%VB(6,132)=0.357142857142858E0_DP
 S_E%VA(6,133)=0.761904761904761E0_DP
 S_E%VB(6,134)=-0.607142857142857E0_DP
 S_E%VA(6,135)=-0.357142857142857E0_DP
 S_E%VB(6,136)=7.142857142857142E-002_DP
 S_E%VA(6,148)=0.166666666666666E0_DP
 S_E%VB(6,149)=-0.416666666666664E0_DP
 S_E%VA(6,150)=-0.916666666666666E0_DP
 S_E%VB(6,151)=0.750000000000000E0_DP
 S_E%VA(6,152)=0.500000000000000E0_DP
 S_E%VB(6,153)=-0.166666666666667E0_DP
 S_E%VA(6,166)=-0.200000000000000E0_DP
 S_E%VB(6,167)=0.500000000000000E0_DP
 S_E%VA(6,168)=1.16666666666667E0_DP
 S_E%VB(6,169)=-1.00000000000000E0_DP
 S_E%VA(6,170)=-1.00000000000000E0_DP
 S_E%VA(6,185)=0.250000000000000E0_DP
 S_E%VB(6,186)=-0.624999999999999E0_DP
 S_E%VA(6,187)=-1.66666666666667E0_DP
 S_E%VB(6,188)=2.50000000000000E0_DP
 S_E%VA(6,205)=-0.333333333333333E0_DP
 S_E%VB(6,206)=0.833333333333333E0_DP
 S_E%VA(6,207)=3.33333333333333E0_DP
 S_E%VA(6,226)=0.500000000000000E0_DP
 S_E%VB(6,227)=-2.50000000000000E0_DP
 S_E%VA(6,248)=-1.00000000000000E0_DP
 S_E%VB(6,270)=0.166666666666667E0_DP
 S_E%VB(7,78)=-3.544372294372294E-002_DP
 S_E%VA(7,90)=-0.223214285714286E0_DP
 S_E%VB(7,91)=4.166666666666666E-002_DP
 S_E%VB(7,103)=0.723214285714286E0_DP
 S_E%VA(7,104)=0.261904761904762E0_DP
 S_E%VB(7,105)=-5.158730158730158E-002_DP
 S_E%VA(7,117)=1.02678571428571E0_DP
 S_E%VB(7,118)=-0.857142857142857E0_DP
 S_E%VA(7,119)=-0.321428571428571E0_DP
 S_E%VB(7,120)=7.142857142857142E-002_DP
 S_E%VB(7,132)=-1.14285714285714E0_DP
 S_E%VA(7,133)=-1.21428571428571E0_DP
 S_E%VB(7,134)=1.07142857142857E0_DP
 S_E%VA(7,135)=0.428571428571429E0_DP
 S_E%VB(7,136)=-0.142857142857143E0_DP
 S_E%VA(7,148)=-0.500000000000000E0_DP
 S_E%VB(7,149)=1.37500000000000E0_DP
 S_E%VA(7,150)=1.50000000000000E0_DP
 S_E%VB(7,151)=-1.50000000000000E0_DP
 S_E%VA(7,152)=-1.00000000000000E0_DP
 S_E%VB(7,165)=0.199999999999999E0_DP
 S_E%VA(7,166)=0.600000000000000E0_DP
 S_E%VB(7,167)=-1.75000000000000E0_DP
 S_E%VA(7,168)=-2.00000000000000E0_DP
 S_E%VB(7,169)=3.00000000000000E0_DP
 S_E%VB(7,184)=-0.249999999999998E0_DP
 S_E%VA(7,185)=-0.750000000000000E0_DP
 S_E%VB(7,186)=2.50000000000000E0_DP
 S_E%VA(7,187)=5.00000000000000E0_DP
 S_E%VB(7,204)=0.333333333333333E0_DP
 S_E%VA(7,205)=1.00000000000000E0_DP
 S_E%VB(7,206)=-5.00000000000000E0_DP
 S_E%VB(7,225)=-0.500000000000000E0_DP
 S_E%VA(7,226)=-3.00000000000000E0_DP
 S_E%VB(7,247)=1.00000000000000E0_DP
 S_E%VA(7,269)=0.142857142857143E0_DP
 S_E%VB(8,78)=3.409090909090909E-002_DP
 S_E%VA(8,90)=0.291666666666667E0_DP
 S_E%VB(8,91)=-4.166666666666666E-002_DP
 S_E%VB(8,103)=-0.916666666666667E0_DP
 S_E%VA(8,104)=-0.361111111111111E0_DP
 S_E%VB(8,105)=5.555555555555555E-002_DP
 S_E%VA(8,117)=-2.00000000000000E0_DP
 S_E%VB(8,118)=1.12500000000000E0_DP
 S_E%VA(8,119)=0.500000000000000E0_DP
 S_E%VB(8,120)=-0.125000000000000E0_DP
 S_E%VB(8,132)=2.12500000000000E0_DP
 S_E%VA(8,133)=2.50000000000000E0_DP
 S_E%VB(8,134)=-1.50000000000000E0_DP
 S_E%VA(8,135)=-1.00000000000000E0_DP
 S_E%VA(8,148)=1.92500000000000E0_DP
 S_E%VB(8,149)=-2.62500000000000E0_DP
 S_E%VA(8,150)=-3.50000000000000E0_DP
 S_E%VB(8,151)=3.50000000000000E0_DP
 S_E%VB(8,165)=-0.700000000000001E0_DP
 S_E%VA(8,166)=-2.45000000000000E0_DP
 S_E%VB(8,167)=3.50000000000000E0_DP
 S_E%VA(8,168)=7.00000000000000E0_DP
 S_E%VA(8,183)=-0.250000000000000E0_DP
 S_E%VB(8,184)=0.875000000000000E0_DP
 S_E%VA(8,185)=3.50000000000000E0_DP
 S_E%VB(8,186)=-8.75000000000000E0_DP
 S_E%VA(8,203)=0.333333333333333E0_DP
 S_E%VB(8,204)=-1.16666666666667E0_DP
 S_E%VA(8,205)=-7.00000000000000E0_DP
 S_E%VA(8,224)=-0.500000000000000E0_DP
 S_E%VB(8,225)=3.50000000000000E0_DP
 S_E%VA(8,246)=1.00000000000000E0_DP
 S_E%VB(8,268)=-0.125000000000000E0_DP
 S_E%VB(9,78)=-4.040404040404040E-002_DP
 S_E%VA(9,90)=-0.333333333333333E0_DP
 S_E%VB(9,91)=5.555555555555555E-002_DP
 S_E%VB(9,103)=1.44444444444444E0_DP
 S_E%VA(9,104)=0.444444444444444E0_DP
 S_E%VB(9,105)=-0.111111111111111E0_DP
 S_E%VA(9,117)=3.00000000000000E0_DP
 S_E%VB(9,118)=-2.00000000000000E0_DP
 S_E%VA(9,119)=-1.00000000000000E0_DP
 S_E%VB(9,132)=-5.00000000000000E0_DP
 S_E%VA(9,133)=-4.00000000000000E0_DP
 S_E%VB(9,134)=4.00000000000000E0_DP
 S_E%VA(9,148)=-4.20000000000000E0_DP
 S_E%VB(9,149)=7.00000000000000E0_DP
 S_E%VA(9,150)=9.33333333333333E0_DP
 S_E%VB(9,165)=3.26666666666667E0_DP
 S_E%VA(9,166)=5.60000000000000E0_DP
 S_E%VB(9,167)=-14.0000000000000E0_DP
 S_E%VA(9,183)=1.00000000000000E0_DP
 S_E%VB(9,184)=-4.66666666666666E0_DP
 S_E%VA(9,185)=-14.0000000000000E0_DP
 S_E%VB(9,202)=-0.333333333333333E0_DP
 S_E%VA(9,203)=-1.33333333333333E0_DP
 S_E%VB(9,204)=9.33333333333333E0_DP
 S_E%VB(9,223)=0.500000000000000E0_DP
 S_E%VA(9,224)=4.00000000000000E0_DP
 S_E%VB(9,245)=-1.00000000000000E0_DP
 S_E%VA(9,267)=-0.111111111111111E0_DP
 S_E%VB(10,78)=4.545454545454546E-002_DP
 S_E%VA(10,90)=0.500000000000000E0_DP
 S_E%VB(10,91)=-0.100000000000000E0_DP
 S_E%VB(10,103)=-2.00000000000000E0_DP
 S_E%VA(10,104)=-1.00000000000000E0_DP
 S_E%VA(10,117)=-6.00000000000000E0_DP
 S_E%VB(10,118)=4.50000000000000E0_DP
 S_E%VB(10,132)=9.00000000000000E0_DP
 S_E%VA(10,133)=12.0000000000000E0_DP
 S_E%VA(10,148)=12.6000000000000E0_DP
 S_E%VB(10,149)=-21.0000000000000E0_DP
 S_E%VB(10,165)=-8.40000000000000E0_DP
 S_E%VA(10,166)=-25.2000000000000E0_DP
 S_E%VA(10,183)=-6.00000000000000E0_DP
 S_E%VB(10,184)=21.0000000000000E0_DP
 S_E%VB(10,202)=1.50000000000000E0_DP
 S_E%VA(10,203)=12.0000000000000E0_DP
 S_E%VA(10,222)=0.500000000000000E0_DP
 S_E%VB(10,223)=-4.50000000000000E0_DP
 S_E%VA(10,244)=-1.00000000000000E0_DP
 S_E%VB(10,266)=0.100000000000000E0_DP
 S_E%VB(11,78)=-9.090909090909091E-002_DP
 S_E%VA(11,90)=-1.00000000000000E0_DP
 S_E%VB(11,103)=5.00000000000000E0_DP
 S_E%VA(11,117)=15.0000000000000E0_DP
 S_E%VB(11,132)=-30.0000000000000E0_DP
 S_E%VA(11,148)=-42.0000000000000E0_DP
 S_E%VB(11,165)=42.0000000000000E0_DP
 S_E%VA(11,183)=30.0000000000000E0_DP
 S_E%VB(11,202)=-15.0000000000000E0_DP
 S_E%VA(11,222)=-5.00000000000000E0_DP
 S_E%VB(11,243)=1.00000000000000E0_DP
 S_E%VA(11,265)=9.090909090909091E-002_DP
 S_E%VB(12,66)=-8.333333333333333E-002_DP
 S_E%VA(12,77)=-1.00000000000000E0_DP
 S_E%VB(12,89)=5.50000000000000E0_DP
 S_E%VA(12,102)=18.3333333333333E0_DP
 S_E%VB(12,116)=-41.2500000000000E0_DP
 S_E%VA(12,131)=-66.0000000000000E0_DP
 S_E%VB(12,147)=77.0000000000000E0_DP
 S_E%VA(12,164)=66.0000000000000E0_DP
 S_E%VB(12,182)=-41.2500000000000E0_DP
 S_E%VA(12,201)=-18.3333333333333E0_DP
 S_E%VB(12,221)=5.50000000000000E0_DP
 S_E%VA(12,242)=1.00000000000000E0_DP
 S_E%VB(12,264)=-8.333333333333333E-002_DP
 S_E%VB(13,55)=-7.692307692307693E-002_DP
 S_E%VA(13,65)=-1.00000000000000E0_DP
 S_E%VB(13,76)=6.00000000000000E0_DP
 S_E%VA(13,88)=22.0000000000000E0_DP
 S_E%VB(13,101)=-55.0000000000000E0_DP
 S_E%VA(13,115)=-99.0000000000000E0_DP
 S_E%VB(13,130)=132.000000000000E0_DP
 S_E%VA(13,146)=132.000000000000E0_DP
 S_E%VB(13,163)=-99.0000000000000E0_DP
 S_E%VA(13,181)=-55.0000000000000E0_DP
 S_E%VB(13,200)=22.0000000000000E0_DP
 S_E%VA(13,220)=6.00000000000000E0_DP
 S_E%VB(13,241)=-1.00000000000000E0_DP
 S_E%VA(13,263)=-7.692307692307693E-002_DP
 S_E%VB(14,45)=-7.142857142857142E-002_DP
 S_E%VA(14,54)=-1.00000000000000E0_DP
 S_E%VB(14,64)=6.50000000000000E0_DP
 S_E%VA(14,75)=26.0000000000000E0_DP
 S_E%VB(14,87)=-71.5000000000000E0_DP
 S_E%VA(14,100)=-143.000000000000E0_DP
 S_E%VB(14,114)=214.500000000000E0_DP
 S_E%VA(14,129)=245.142857142857E0_DP
 S_E%VB(14,145)=-214.500000000000E0_DP
 S_E%VA(14,162)=-143.000000000000E0_DP
 S_E%VB(14,180)=71.5000000000000E0_DP
 S_E%VA(14,199)=26.0000000000000E0_DP
 S_E%VB(14,219)=-6.50000000000000E0_DP
 S_E%VA(14,240)=-1.00000000000000E0_DP
 S_E%VB(14,262)=7.142857142857142E-002_DP
 S_E%VB(15,36)=-6.666666666666667E-002_DP
 S_E%VA(15,44)=-1.00000000000000E0_DP
 S_E%VB(15,53)=7.00000000000000E0_DP
 S_E%VA(15,63)=30.3333333333333E0_DP
 S_E%VB(15,74)=-91.0000000000000E0_DP
 S_E%VA(15,86)=-200.200000000000E0_DP
 S_E%VB(15,99)=333.666666666667E0_DP
 S_E%VA(15,113)=429.000000000000E0_DP
 S_E%VB(15,128)=-429.000000000000E0_DP
 S_E%VA(15,144)=-333.666666666667E0_DP
 S_E%VB(15,161)=200.200000000000E0_DP
 S_E%VA(15,179)=91.0000000000000E0_DP
 S_E%VB(15,198)=-30.3333333333333E0_DP
 S_E%VA(15,218)=-7.00000000000000E0_DP
 S_E%VB(15,239)=1.00000000000000E0_DP
 S_E%VA(15,261)=6.666666666666667E-002_DP
 S_E%VB(16,28)=-6.250000000000000E-002_DP
 S_E%VA(16,35)=-1.00000000000000E0_DP
 S_E%VB(16,43)=7.50000000000000E0_DP
 S_E%VA(16,52)=35.0000000000000E0_DP
 S_E%VB(16,62)=-113.750000000000E0_DP
 S_E%VA(16,73)=-273.000000000000E0_DP
 S_E%VB(16,85)=500.500000000000E0_DP
 S_E%VA(16,98)=715.000000000000E0_DP
 S_E%VB(16,112)=-804.375000000000E0_DP
 S_E%VA(16,127)=-715.000000000000E0_DP
 S_E%VB(16,143)=500.500000000000E0_DP
 S_E%VA(16,160)=273.000000000000E0_DP
 S_E%VB(16,178)=-113.750000000000E0_DP
 S_E%VA(16,197)=-35.0000000000000E0_DP
 S_E%VB(16,217)=7.50000000000000E0_DP
 S_E%VA(16,238)=1.00000000000000E0_DP
 S_E%VB(16,260)=-6.250000000000000E-002_DP
 S_E%VB(17,21)=-5.882352941176471E-002_DP
 S_E%VA(17,27)=-1.00000000000000E0_DP
 S_E%VB(17,34)=8.00000000000000E0_DP
 S_E%VA(17,42)=40.0000000000000E0_DP
 S_E%VB(17,51)=-140.000000000000E0_DP
 S_E%VA(17,61)=-364.000000000000E0_DP
 S_E%VB(17,72)=728.000000000000E0_DP
 S_E%VA(17,84)=1144.00000000000E0_DP
 S_E%VB(17,97)=-1430.00000000000E0_DP
 S_E%VA(17,111)=-1430.00000000000E0_DP
 S_E%VB(17,126)=1144.00000000000E0_DP
 S_E%VA(17,142)=728.000000000000E0_DP
 S_E%VB(17,159)=-364.000000000000E0_DP
 S_E%VA(17,177)=-140.000000000000E0_DP
 S_E%VB(17,196)=40.0000000000000E0_DP
 S_E%VA(17,216)=8.00000000000000E0_DP
 S_E%VB(17,237)=-1.00000000000000E0_DP
 S_E%VA(17,259)=-5.882352941176471E-002_DP
 S_E%VB(18,15)=-5.555555555555555E-002_DP
 S_E%VA(18,20)=-1.00000000000000E0_DP
 S_E%VB(18,26)=8.50000000000000E0_DP
 S_E%VA(18,33)=45.3333333333333E0_DP
 S_E%VB(18,41)=-170.000000000000E0_DP
 S_E%VA(18,50)=-476.000000000000E0_DP
 S_E%VB(18,60)=1031.33333333333E0_DP
 S_E%VA(18,71)=1768.00000000000E0_DP
 S_E%VB(18,83)=-2431.00000000000E0_DP
 S_E%VA(18,96)=-2701.11111111111E0_DP
 S_E%VB(18,110)=2431.00000000000E0_DP
 S_E%VA(18,125)=1768.00000000000E0_DP
 S_E%VB(18,141)=-1031.33333333333E0_DP
 S_E%VA(18,158)=-476.000000000000E0_DP
 S_E%VB(18,176)=170.000000000000E0_DP
 S_E%VA(18,195)=45.3333333333333E0_DP
 S_E%VB(18,215)=-8.50000000000000E0_DP
 S_E%VA(18,236)=-1.00000000000000E0_DP
 S_E%VB(18,258)=5.555555555555555E-002_DP
 S_E%VB(19,10)=-5.263157894736842E-002_DP
 S_E%VA(19,14)=-1.00000000000000E0_DP
 S_E%VB(19,19)=9.00000000000000E0_DP
 S_E%VA(19,25)=51.0000000000000E0_DP
 S_E%VB(19,32)=-204.000000000000E0_DP
 S_E%VA(19,40)=-612.000000000000E0_DP
 S_E%VB(19,49)=1428.00000000000E0_DP
 S_E%VA(19,59)=2652.00000000000E0_DP
 S_E%VB(19,70)=-3978.00000000000E0_DP
 S_E%VA(19,82)=-4862.00000000000E0_DP
 S_E%VB(19,95)=4862.00000000000E0_DP
 S_E%VA(19,109)=3978.00000000000E0_DP
 S_E%VB(19,124)=-2652.00000000000E0_DP
 S_E%VA(19,140)=-1428.00000000000E0_DP
 S_E%VB(19,157)=612.000000000000E0_DP
 S_E%VA(19,175)=204.000000000000E0_DP
 S_E%VB(19,194)=-51.0000000000000E0_DP
 S_E%VA(19,214)=-9.00000000000000E0_DP
 S_E%VB(19,235)=1.00000000000000E0_DP
 S_E%VA(19,257)=5.263157894736842E-002_DP
 S_E%VB(20,6)=-5.000000000000000E-002_DP
 S_E%VA(20,9)=-1.00000000000000E0_DP
 S_E%VB(20,13)=9.50000000000000E0_DP
 S_E%VA(20,18)=57.0000000000000E0_DP
 S_E%VB(20,24)=-242.250000000000E0_DP
 S_E%VA(20,31)=-775.200000000000E0_DP
 S_E%VB(20,39)=1938.00000000000E0_DP
 S_E%VA(20,48)=3876.00000000000E0_DP
 S_E%VB(20,58)=-6298.50000000000E0_DP
 S_E%VA(20,69)=-8398.00000000000E0_DP
 S_E%VB(20,81)=9237.80000000000E0_DP
 S_E%VA(20,94)=8398.00000000000E0_DP
 S_E%VB(20,108)=-6298.50000000000E0_DP
 S_E%VA(20,123)=-3876.00000000000E0_DP
 S_E%VB(20,139)=1938.00000000000E0_DP
 S_E%VA(20,156)=775.200000000000E0_DP
 S_E%VB(20,174)=-242.250000000000E0_DP
 S_E%VA(20,193)=-57.0000000000000E0_DP
 S_E%VB(20,213)=9.50000000000000E0_DP
 S_E%VA(20,234)=1.00000000000000E0_DP
 S_E%VB(20,256)=-5.000000000000000E-002_DP
 S_E%VB(21,3)=-4.761904761904762E-002_DP
 S_E%VA(21,5)=-1.00000000000000E0_DP
 S_E%VB(21,8)=10.0000000000000E0_DP
 S_E%VA(21,12)=63.3333333333333E0_DP
 S_E%VB(21,17)=-285.000000000000E0_DP
 S_E%VA(21,23)=-969.000000000000E0_DP
 S_E%VB(21,30)=2584.00000000000E0_DP
 S_E%VA(21,38)=5537.14285714286E0_DP
 S_E%VB(21,47)=-9690.00000000000E0_DP
 S_E%VA(21,57)=-13996.6666666667E0_DP
 S_E%VB(21,68)=16796.0000000000E0_DP
 S_E%VA(21,80)=16796.0000000000E0_DP
 S_E%VB(21,93)=-13996.6666666667E0_DP
 S_E%VA(21,107)=-9690.00000000000E0_DP
 S_E%VB(21,122)=5537.14285714286E0_DP
 S_E%VA(21,138)=2584.00000000000E0_DP
 S_E%VB(21,155)=-969.000000000000E0_DP
 S_E%VA(21,173)=-285.000000000000E0_DP
 S_E%VB(21,192)=63.3333333333333E0_DP
 S_E%VA(21,212)=10.0000000000000E0_DP
 S_E%VB(21,233)=-1.00000000000000E0_DP
 S_E%VA(21,255)=-4.761904761904762E-002_DP
 S_E%VB(22,1)=-4.545454545454546E-002_DP
 S_E%VA(22,2)=-1.00000000000000E0_DP
 S_E%VB(22,4)=10.5000000000000E0_DP
 S_E%VA(22,7)=70.0000000000000E0_DP
 S_E%VB(22,11)=-332.500000000000E0_DP
 S_E%VA(22,16)=-1197.00000000000E0_DP
 S_E%VB(22,22)=3391.50000000000E0_DP
 S_E%VA(22,29)=7752.00000000000E0_DP
 S_E%VB(22,37)=-14535.0000000000E0_DP
 S_E%VA(22,46)=-22610.0000000000E0_DP
 S_E%VB(22,56)=29393.0000000000E0_DP
 S_E%VA(22,67)=32065.0909090909E0_DP
 S_E%VB(22,79)=-29393.0000000000E0_DP
 S_E%VA(22,92)=-22610.0000000000E0_DP
 S_E%VB(22,106)=14535.0000000000E0_DP
 S_E%VA(22,121)=7752.00000000000E0_DP
 S_E%VB(22,137)=-3391.50000000000E0_DP
 S_E%VA(22,154)=-1197.00000000000E0_DP
 S_E%VB(22,172)=332.500000000000E0_DP
 S_E%VA(22,191)=70.0000000000000E0_DP
 S_E%VB(22,211)=-10.5000000000000E0_DP
 S_E%VA(22,232)=-1.00000000000000E0_DP
 S_E%VB(22,254)=4.545454545454546E-002_DP

end  subroutine set_s_e_mcmillan
!!!!!!!!!!!!!!!  New routines for solving Maxwell's Equations !!!!!!!!!!!!!!!


 SUBROUTINE  get_bend_electric_coeff(s_b0t,NO1,h00,verb)
    implicit none
    integer no,i,k,j(2),NO1,l,mf,m,n
    type(taylor) x,y,h,df,ker,sol
    type(complextaylor) z 
    type(taylor) f,kick_x,kick_y
    type(damap) y0
    type(taylor), allocatable :: fs(:)
    real(dp) h0,cker,cker1
     TYPE(B_CYL), intent(inout) :: s_b0t
      logical(lp),optional :: verb    
      real(dp),optional :: h00    

    no=sector_nmul

      if(present(verb)) call kanalnummer(mf,"internal_sol_ele.txt")
    call init(no1,1,0,0)
if(mcmillan) then
 allocate(fs(no1))
call alloc(fs,no1)
endif

    call alloc(x,y,kick_x,kick_y)
    call alloc(z)
    call alloc(f,h,df,ker,sol)
    call alloc(y0)


      x=1.d0.mono.1
      y=1.d0.mono.2
      y0=1
      y0%v(2)=0
      z=x+i_*y      
!!!  


    h0=1.d0
    if(present(h00)) h0=h00

    h=(1.d0+h0*x)
    do k=1,no1
    !  erect multipole
    f=dreal(-z**K/K)
    
    df=f

    do i=k,no-1 !k+1
     df=-h0*(df.d.1)/h
     call  invert_laplace(df)
      sol=f+df
      sol=-(sol.d.1) 
      j=0
      j(1)=i
 
      cker=(sol.sub.j)
      df=df-cker*dreal(-z**(I+1)/(I+1))
      f=f+df
      f=f.cut.(no+1)
     enddo

    if(present(verb)) then

        write(mf,*) " Normal ",k
         call clean_taylor(f,f,1.d-10)
        call print(f,mf)
         y=((f.d.1).d.1)+((f.d.2).d.2)+h0*(f.d.1)/h
        y=y.cut.(no-1)
       ! call print(y,6)
        write(mf,*) full_abs(y)
        kick_x=-(f.d.1)   ! electric field
        kick_y=-(f.d.2)
        call print(kick_x,mf)
        call print(kick_Y,mf)
   endif

       kick_x=-(f.d.1)  ! electric field
       kick_y=-(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)

        s_b0t%b_x(k,l)=kick_x.sub.j
        s_b0t%b_y(k,l)=kick_y.sub.j
        s_b0t%vb(k,l)=f.sub.j
       enddo

if(mcmillan) fs(k)=f
enddo    



if(mcmillan) then

do k=1,no

 
f=fs(k)
 

 do m=k+1,no  !no,k+1,-1


 do n=m,0,-1
  j=0
  j(1)=m-n
  j(2)=n
  cker=fs(m).sub.j
  if(abs(cker)>1.d-10) exit
 enddo
 
 
  cker=f.sub.j
 
  cker1=fs(m).sub.j
  f=f-(cker/cker1)*fs(m)
 call clean_taylor(f,f,1.d-10)
 
 enddo
 
       kick_x=-(f.d.1)  ! electric field
       kick_y=-(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)

        s_b0t%b_x(k,l)=kick_x.sub.j
        s_b0t%b_y(k,l)=kick_y.sub.j
        s_b0t%vb(k,l)=f.sub.j
       enddo

enddo
        


endif

    do k=1,no1
    !  Skew multipole 
 
    f=aimag(-z**K/K)   ! check this
    df=f
    do i=k,no-1 !k+1
     df=-h0*(df.d.1)/h
     call  invert_laplace(df)
      sol=f+df
      sol=-(sol.d.2)
      j=0
      j(1)=i
      cker=(sol.sub.j)
      df=df-cker*aimag(-z**(I+1)/(I+1))
      f=f+df
      f=f.cut.(no+1)
     enddo
     if(present(verb)) then
        write(mf,*) " Skew ",k
        call clean_taylor(f,f,1.d-10)
        call print(f,mf)
         y=((f.d.1).d.1)+((f.d.2).d.2)+h0*(f.d.1)/h
        y=y.cut.(no-1)
      !  call print(y,6)
        write(mf,*) full_abs(y)
        kick_x=-(f.d.1)   ! electric field
        kick_y=-(f.d.2)
        call print(kick_x,mf)
        call print(kick_Y,mf)
       endif

       kick_x=-(f.d.1)   ! electric field
       kick_y=-(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)
        s_b0t%a_x(k,l)=kick_x.sub.j
        s_b0t%a_y(k,l)=kick_y.sub.j
        s_b0t%va(k,l)=f.sub.j
       enddo
    

if(mcmillan) fs(k)=f
enddo    



if(mcmillan) then

do k=1,no


f=fs(k)


 do m=k+1,no  !no,k+1,-1


 do n=m,0,-1
  j=0
  j(1)=m-n
  j(2)=n
  cker=fs(m).sub.j
  if(abs(cker)>1.d-10) exit
 enddo

  cker=f.sub.j

  cker1=fs(m).sub.j
  f=f-(cker/cker1)*fs(m)
 call clean_taylor(f,f,1.d-10)

 enddo

       kick_x=-(f.d.1)   ! electric field
       kick_y=-(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)
        s_b0t%a_x(k,l)=kick_x.sub.j
        s_b0t%a_y(k,l)=kick_y.sub.j
        s_b0t%va(k,l)=f.sub.j
       enddo

enddo
        


endif

    if(present(verb)) close(mf)

    call kill(x,y,kick_x,kick_y)
    call kill(z)
    call kill(f,h,df,ker,sol)
    call kill(y0)
if(mcmillan) then
 call kill(fs,no)
 deallocate(fs)
endif
    end subroutine get_bend_electric_coeff

 SUBROUTINE  get_bend_magnetic_potential(s_b0t,NO1,h00,verb)
! trying to get magnetic field from scalar potential
    implicit none
    integer no,i,k,j(2),NO1,l,mf,m,n
    type(taylor) x,y,h,df,ker,sol
    type(complextaylor) z 
    type(taylor) f,kick_x,kick_y
    type(damap) y0
    real(dp) h0,cker,cker1
    type(taylor), allocatable :: fs(:)
     TYPE(B_CYL), intent(inout) :: s_b0t
      logical(lp),optional :: verb    
      real(dp),optional :: h00    


      if(present(verb)) call kanalnummer(mf,"internal_mag_pot.txt")
    call init(no1,1,0,0)
if(mcmillan) then
 allocate(fs(no1))
call alloc(fs,no1)
endif

    no=sector_nmul

    call alloc(x,y,kick_x,kick_y)
    call alloc(z)
    call alloc(f,h,df,ker,sol)
    call alloc(y0)


      x=1.d0.mono.1
      y=1.d0.mono.2
      y0=1
      y0%v(2)=0
      z=x+i_*y      
!!!  


    h0=1.d0
    if(present(h00)) h0=h00

    h=(1.d0+h0*x)
    do k=1,no1
    !  erect multipole
    f=dreal(i_*z**K/K)
    
    df=f

    do i=k,no-1 !k+1
     df=-h0*(df.d.1)/h
     call  invert_laplace(df)
      sol=f+df
      sol=-(sol.d.2) 
      j=0
      j(1)=i
      cker=(sol.sub.j)
      df=df-cker*dreal(i_*z**(I+1)/(I+1))
      f=f+df
      f=f.cut.(no+1)
     enddo

    if(present(verb)) then

        write(mf,*) " Normal ",k
         call clean_taylor(f,f,1.d-10)
        call print(f,mf)
         y=((f.d.1).d.1)+((f.d.2).d.2)+h0*(f.d.1)/h
        y=y.cut.(no-1)
       ! call print(y,6)
        write(mf,*) full_abs(y)
       kick_x= h*(f.d.2)  ! magnetic kicks
       kick_y=-h*(f.d.1)
       kick_x=kick_x.cut.no
       kick_y=kick_y.cut.no
         call clean_taylor(kick_x,kick_x,1.d-10)
         call clean_taylor(kick_Y,kick_Y,1.d-10)
        call print(kick_x,mf)
        call print(kick_Y,mf)
       endif

       kick_x= -(f.d.1)  ! magnetic field
       kick_y= -(f.d.2)

       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)

        s_b0t%b_x(k,l)=kick_x.sub.j
        s_b0t%b_y(k,l)=kick_y.sub.j
        s_b0t%vb(k,l)=f.sub.j
       enddo
if(mcmillan) fs(k)=f
enddo    


if(mcmillan) then

do k=1,no


f=fs(k)

 do m=k+1,no  !no,k+1,-1


 do n=m,0,-1
  j=0
  j(1)=m-n
  j(2)=n
  cker=fs(m).sub.j
  if(abs(cker)>1.d-10) exit
 enddo
 
 
  cker=f.sub.j
 
  cker1=fs(m).sub.j
  f=f-(cker/cker1)*fs(m)
 call clean_taylor(f,f,1.d-10)

 enddo

       kick_x= -(f.d.1)  ! magnetic field
       kick_y= -(f.d.2)

       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)

        s_b0t%b_x(k,l)=kick_x.sub.j
        s_b0t%b_y(k,l)=kick_y.sub.j
        s_b0t%vb(k,l)=f.sub.j
       enddo

enddo
        


endif

    do k=1,no1
    !  Skew multipole 
 
    f=dreal(-z**K/K)   ! check this
    df=f
    do i=k,no-1 !k+1
     df=-h0*(df.d.1)/h
     call  invert_laplace(df)
      sol=f+df
      sol=-(sol.d.1)
      j=0
      j(1)=i
      cker=(sol.sub.j)
      df=df-cker*dreal(-z**(I+1)/(I+1))
      f=f+df
      f=f.cut.(no+1)
     enddo
     if(present(verb)) then
        write(mf,*) " Skew ",k
        call clean_taylor(f,f,1.d-10)
        call print(f,mf)
         y=((f.d.1).d.1)+((f.d.2).d.2)+h0*(f.d.1)/h
        y=y.cut.(no-1)
      !  call print(y,6)
        write(mf,*) full_abs(y)
       kick_x= h*(f.d.2)  ! magnetic kicks
       kick_y=-h*(f.d.1)
       kick_x=kick_x.cut.no
       kick_y=kick_y.cut.no
         call clean_taylor(kick_x,kick_x,1.d-10)
         call clean_taylor(kick_Y,kick_Y,1.d-10)
        call print(kick_x,mf)
        call print(kick_Y,mf)
       endif

       kick_x= -(f.d.1)  ! magnetic field
       kick_y= -(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)
        s_b0t%a_x(k,l)=kick_x.sub.j
        s_b0t%a_y(k,l)=kick_y.sub.j
        s_b0t%va(k,l)=f.sub.j
       enddo
if(mcmillan) fs(k)=f
enddo    


if(mcmillan) then

do k=1,no


f=fs(k)

 do m=k+1,no  !no,k+1,-1


 do n=m,0,-1
  j=0
  j(1)=m-n
  j(2)=n
  cker=fs(m).sub.j
  if(abs(cker)>1.d-10) exit
 enddo
 
  cker=f.sub.j
 
  cker1=fs(m).sub.j
  f=f-(cker/cker1)*fs(m)
 call clean_taylor(f,f,1.d-10)

 enddo

       kick_x= -(f.d.1)  ! magnetic field
       kick_y= -(f.d.2)
       call clean_taylor(kick_x,kick_x,1.d-6)
       call clean_taylor(kick_y,kick_y,1.d-6)
       do l=1,s_b0t%n_mono
        j(1)=s_b0t%i(l)
        j(2)=s_b0t%j(l)
        s_b0t%a_x(k,l)=kick_x.sub.j
        s_b0t%a_y(k,l)=kick_y.sub.j
        s_b0t%va(k,l)=f.sub.j
       enddo

enddo
        


endif

    if(present(verb)) close(mf)

    call kill(x,y,kick_x,kick_y)
    call kill(z)
    call kill(f,h,df,ker,sol)
    call kill(y0)
if(mcmillan) then
 call kill(fs,no)
 deallocate(fs)
endif
    end subroutine get_bend_magnetic_potential

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

 
  subroutine make_coef(b,no,ic)
    implicit none
    integer no
    integer i,j,k,a,m , ic

    type(B_CYL) b

 
    b%firsttime=-100
    allocate(b%nmul)
    allocate(b%n_mono)
    b%nmul=no
    b%n_mono=((no+2-ic)*(no+1-ic))/2
    allocate(b%i(b%n_mono),b%j(b%n_mono))
    allocate(b%a_x(no,b%n_mono),b%a_y(no,b%n_mono))
    allocate(b%b_x(no,b%n_mono),b%b_y(no,b%n_mono))    
    allocate(b%va(no,b%n_mono),b%vb(no,b%n_mono))

    do i=1,no
       do j=1,b%n_mono
          b%a_x(i,j)=0.0_dp
          b%a_y(i,j)=0.0_dp
          b%b_x(i,j)=0.0_dp
          b%b_y(i,j)=0.0_dp
          b%va(i,j)=0.0_dp
          b%vb(i,j)=0.0_dp
       enddo
    enddo


    k=0
    m=no-ic
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
       nullify(b%va)
       nullify(b%vb)
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
       deallocate(b%va)
       deallocate(b%vb)
    endif


  end subroutine nul_coef

!!!!!!!!!!!!!!!!!!!!   tree tracking for PTC using stuff in 
  SUBROUTINE SET_TREE_G_complex(T,Ma,factor)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T(:)
    TYPE(c_damap), INTENT(INOUT) :: Ma
    INTEGER N,NP,i,k,j,kq
    real(dp) norm,mat(6,6)
    TYPE(taylor), ALLOCATABLE :: M(:), MG(:)
    TYPE(damap) ms
    integer js(6)
    logical, optional :: factor
    logical fact
    fact=.false.
    if(present(factor)) fact=factor
    
!    np=ma%n+18
    if(ma%n/=6) then
     write(6,*) " you need a 6-d map in SET_TREE_G_complex for PTC "
     stop
    endif
    np=size_tree
! initialized in ptc ini
 !   ind_spin(1,1)=1+ma%n;ind_spin(1,2)=2+ma%n;ind_spin(1,3)=3+ma%n;
 !   ind_spin(2,1)=4+ma%n;ind_spin(2,2)=5+ma%n;ind_spin(2,3)=6+ma%n;
 !   ind_spin(3,1)=7+ma%n;ind_spin(3,2)=8+ma%n;ind_spin(3,3)=9+ma%n;    
 !   k1_spin(1)=1;k2_spin(1)=1;
 !   k1_spin(2)=1;k2_spin(2)=2;
 !   k1_spin(3)=1;k2_spin(3)=3;
 !   k1_spin(4)=2;k2_spin(4)=1;
 !   k1_spin(5)=2;k2_spin(5)=2;
 !   k1_spin(6)=2;k2_spin(6)=3;
 !   k1_spin(7)=3;k2_spin(7)=1;
 !   k1_spin(8)=3;k2_spin(8)=2;
 !   k1_spin(9)=3;k2_spin(9)=3;

   
    ALLOCATE(M(NP))
    CALL ALLOC(M,NP)
    ALLOCATE(Mg(NP))
    CALL ALLOC(mg,NP)
    do i=1,np
     m(i)=0.e0_dp
     mg(i)=0.e0_dp
    enddo
     do i=1,ma%n
      m(i)=ma%v(i)   ! orbital part
     enddo


 

if(use_quaternion) then
    call c_full_norm_quaternion(Ma%q,kq,norm)
    if(kq==-1) then
      do i=0,3
        m(ind_spin(1,1)+i)=ma%q%x(i)
      enddo
    elseif(kq/=-1) then
      m(ind_spin(1,1))=1.0_dp
      do i=ind_spin(1,1)+1,size_tree
        m(i)=0.0_dp
      enddo
    endif
else
    call c_full_norm_spin(Ma%s,k,norm)

    if(k==-1) then
      do i=1,3
      do j=1,3
        m(ind_spin(i,j))=ma%s%s(i,j)
      enddo
      enddo
    else
      do i=1,3
        m(ind_spin(i,i))=1.0e0_dp
      enddo
    endif
endif

 

  



      js=0
     js(1)=1;js(3)=1;js(5)=1; ! q_i(q_f,p_i) and p_f(q_f,p_i)
     call alloc(ms)
     if(fact) then
      ms=(ma.sub.1)**(-1)*ma
     else
       ms=ma
     endif
 

     ms=ms**js
!     do i=1,3
!      mg(i)=ms%v(2*i-1)   !  q_i(q_f,p_i)
!      mg(3+i)=ms%v(2*i)   !  p_f(q_f,p_i)
!     enddo
     do i=1,6
      mg(i)=ms%v(i) 
     enddo
    
     do i=1,3
     do j=1,3
       mg(ind_spin(i,j))=ms%v(2*i-1).d.(2*j-1)  !   Jacobian for Newton search
     enddo
     enddo
      call kill(ms)    
   

     call SET_TREE_g(T(1),m(1:6))

     call SET_TREE_g(T(2),m(7:15))

     call SET_TREE_g(T(3),mg(1:size_tree))

     if(fact) then
      t(3)%rad=ma
       t(3)%factored=.true.
     else
       do i=1,6
        t(3)%rad(i,i)=1.0_dp
       enddo
       t(3)%factored=.false.
     endif

       mat=ma**(-1)
       t(1)%e_ij=matmul(matmul(mat,ma%e_ij),transpose(mat))
 
     call kill(m); call kill(mg);
    deallocate(M);    deallocate(Mg);

  END SUBROUTINE SET_TREE_G_complex

subroutine print_tree_element(t,mf)
implicit none
type(tree_element) t
 
integer i,mf
!   write(mf,'(a204)') t%file
write(mf,'(3(1X,i8))') t%N,t%NP,t%no 
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

subroutine print_tree_elements(t,mf)
implicit none
type(tree_element) t(:)
 
integer i,mf

 do i=1,size(t)
  call print_tree_element(t(i),mf)
 enddo

end subroutine print_tree_elements
 
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

  SUBROUTINE track_TREE_probe_complexr(T,xs,dofix0,dofix,sta,jump,all_map)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T(:)
    logical, optional :: jump,all_map
    type(probe) xs
    real(dp) x(size_tree),x0(size_tree),s0(3,3),r(3,3),dx6,beta,q(3),p(3),qg(3),qf(3)
    real(dp) normb,norm 
    type(quaternion)qu
    integer i,j,k,ier,nrmax,is
    type(internal_state) sta
    logical dofix0,dofix,doit,jumpnot,allmap

    jumpnot=.true.
    if(present(jump)) jumpnot=.not.jump

 
    allmap=.true.
    if(present(all_map)) allmap=all_map

    nrmax=1000
    doit=.true.
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
     x(5)=(2*x(5)+x(5)**2)/(root(1.0_dp/t(1)%beta0**2+2.0_dp*x(5)+x(5)**2)+1.0_dp/t(1)%beta0)
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
     write(6,*) "matinv failed in track_TREE_probe_complexr (Se_status) "
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

!   if(norm>t(3)%eps.and.doit) then
!     if(normb<=norm) doit=.false.
!     normb=norm
!   else
!     if(normb<=norm) then 
!       x(1)=qf(1)
!       x(3)=qf(2)
!       x(5)=qf(3)
!       x(2)=x0(2)
!       x(4)=x0(4)
!       x(6)=x0(6)       

!       if(allmap) x(1:6)=matmul(t(3)%rad,x(1:6))
!       exit
!     endif
!     normb=norm
!   endif

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

       if(allmap) x(1:6)=matmul(t(3)%rad,x(1:6))
       exit
     endif
     normb=norm
   endif



enddo  ! is 
 if(is>nrmax-10) then
   xs%u=.true.
  check_stable=.false.
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

    call orthonormalise(r)
    
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
      beta=root(1.0_dp+2.0_dp*x(5)/t(1)%beta0+x(5)**2)/(1.0_dp/t(1)%BETA0 + x(5))
      x(6)=x(6)-dx6+beta*dx6 +  (beta/t(1)%beta0-1.0_dp)*t(1)%ds
      x(5)=(2.0_dp*x(5)/t(1)%beta0+x(5)**2)/(root(1.0_dp+2.0_dp*x(5)/t(1)%beta0+x(5)**2)+1.0_dp)
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
      stop
    endif 
  end SUBROUTINE orthonormaliser

SUBROUTINE track_TREE_probe_complexp_new(T,xs,dofix0,dofix,sta)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T(:)
    type(probe_8) xs
    type(probe) xs0
    type(real_8) x(size_tree),x0(size_tree),s0(3,3),r(3,3),dx6,beta,ds
    real(dp) m(6,6),xi(6),norm,z0(6)
    type(damap) dm,md,iq
    type(c_damap) m0,mt
    type(quaternion_8) qu
    integer i,j,k,o
    type(internal_state) sta
    logical dofix0,dofix
    integer, allocatable :: js(:)
 
    call alloc(x,size_tree)
    call alloc(x0,size_tree)
    call alloc(dx6,beta)
    do i=1,3
    do j=1,3
     call alloc(s0(i,j))
     call alloc(r(i,j))
    enddo
    enddo


  ! if(c_%nd2/=6) then
  !  write(6,*) " problem in track_TREE_probe_complexp for the moment if nd2/=6 "
  !  stop
  ! else
    call alloc(m0,mt)
    m0=xs

    do o=1,6
     z0(o)=xs%x(o)
     enddo
 !   call print(m0,6)
 !   write(6,*) z0
    call KILL(xs)
    call alloc(xs)
    mt=1
    xs0=0
    xs0=z0
    xs=mt+xs0
 !  endif


    if(sta%envelope.and.c_%no>1) then
    call alloc(dm)
     dm=xs
     m=(dm.sub.1)**(-1)
        xs%e_ij=xs%e_ij+matmul(matmul(m,t(1)%e_ij),transpose(m))
    call kill(dm)
     endif

    do i=1,6
      x(i)=xs%x(i)
      x0(i)=xs%x(i)
      x(i+6)=x(i)
    enddo

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
      x(i+6)=x(i)
     enddo
    endif

    if(sta%radiation) then
     do i=1,6
      x(i)=0.0_dp
     enddo

     do i=1,6
     do j=1,6
      x(i)=t(1)%rad(i,j)*x(j+6)+x(i)
     enddo
     enddo


      do i=1,6
       x0(i)=x(i)
       x(i+6)=x(i)
     enddo

    endif

 
    if(t(3)%symptrack) then
     xs0=0
     do i=1,6
      xs0%x(i)=x0(i)
      xi(i)=x0(i)
     enddo
 
      call  track_TREE_probe_complexr(T,xs0,.false.,.false.,sta,jump=.true.,all_map=.not.t(3)%factored)

!!! compute map  for speed up
     norm=0.d0
     do i=1,6
      norm=norm+abs(x(1).sub.'1')
     enddo
!     write(6,*) "norm = ",norm
     if(norm>0) then
       call alloc(dm,md,iq)
       allocate(js(c_%nd2))
              do i=1,3   !c_%nd
               xi(2*i-1)=xs0%x(2*i-1)
              enddo

              do i=1,c_%nd2
              x0(i)=xi(i)+(1.0_dp.mono.i)
              enddo

              if(c_%nd2==4.and.C_%NPARA==5) then
             !  x0(5)=xi(5)+(1.0_dp.mono.5)
               x0(6)=0.0_dp !xi(6)
              elseif(C_%NPARA==4) then
               x0(5)=xi(5)
               x0(6)=0.0_dp !xi(6)       
              endif

          call track_TREE_G_complex(T(3),X0(1:15))
       js=0
      do i=1,c_%nd2
          if(mod(i,2)==1) js(i)=1
          dm%v(i)=x0(i)-(x0(i).sub.'0') 
      enddo
        dm=dm**(js)
        do i=1,c_%nd2
          md%v(i)=x(i)-(x(i).sub.'0') 
        enddo 
          if(c_%nd2==4) then
            do i=1,c_%nd
              iq%v(2*i-1)=dm%v(2*i-1) 
              iq%v(2*i)=1.0_dp.mono.(2*i)
            enddo 
            x0(6)=x0(6)*iq  ! partial invertion undone
            x0(6)=x0(6)*md  ! previous line concatenated
          endif
          md=dm*md

        do i=1,c_%nd2
         x(i)=md%v(i)+xs0%x(i)
        enddo

      if(c_%nd2==4) then
       x(6)=x0(6)+x(6)
      endif

        call kill(dm,md,iq)
        deallocate(js)
     else
       do i=1,6  !c_%nd2
         x(i)=xs0%x(i)
        enddo
     endif
       do i=1,6
        x0(i)=0.0_dp
       enddo

       do i=1,6
        do j=1,6
       x0(i)=t(3)%rad(i,j)*x(j)+x0(i)  
       enddo
      enddo

        do i=1,6
         x(i)=x0(i)
        enddo

     else
       call track_TREE_G_complex(T(1),X(1:6))
     endif
 





 
    if(sta%spin) then  ! spin
    call track_TREE_G_complex(T(2),X(7:15))
      if(xs%use_q) then
call alloc(qu)
call alloc(ds)
       do k=0,3
         qu%x(k)=x(7+k)
       enddo 
 
       xs%q=qu*xs%q
        ds=1.0_dp/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)
            xs%q%x(1)=ds*xs%q%x(1)
            xs%q%x(2)=ds*xs%q%x(2)
            xs%q%x(3)=ds*xs%q%x(3)
            xs%q%x(0)=ds*xs%q%x(0)
call KILL(qu)
call KILL(ds)
   else

    do i=1,3
    do j=1,3
     r(i,j)=x(ind_spin(i,j))
    enddo
    enddo

    call orthonormalise(r)

    do k=1,3
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
endif ! spin
endif


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
    call kill(dx6)
    else
        if(sta%totalpath==1) then
        x(6)=x(6)+t(1)%ds/t(1)%beta0 
       endif     
    endif


    do i=1,6
      xs%x(i)=x(i)
      z0(i)=x(i)
    enddo

    mt=xs
    xs=mt*m0
    
    do i=1,6
     xs%x(i)=xs%x(i)-(xs%x(i).sub.'0')+z0(i)
    enddo
    


    call kill(m0,mt)
    call kill(dx6,beta)
    call kill(x0,size_tree)
    call kill(x,size_tree)  
    do i=1,3
    do j=1,3
     call kill(s0(i,j))
     call kill(r(i,j))
    enddo
    enddo

  end SUBROUTINE track_TREE_probe_complexp_new

  SUBROUTINE orthonormalisep(r)
   implicit none
   type(real_8)  r(3,3),id(3,3),rt(3,3)
    real(dp) eps,a,ab
   integer nmax,i,j,k
! Furmanizing the rotation  
    do i=1,3
    do j=1,3
     call alloc(id(i,j))
     call alloc(rt(i,j))
    enddo
    enddo
    eps=1.d-8
    nmax=1000
    do i=1,3
      id(i,i)=1.5e0_dp
    enddo
    ab=1.d8
    do i=1,nmax
    ! rt=matmul(r,transpose(r))
    ! r= matmul((id-0.5e0_dp*rt),r)

      call furman_rrt(r,r,rt)

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
     write(6,*) i, a, "did not converge in orthonormalisep"
     ! stop
    endif 
    do i=1,3
    do j=1,3
     call kill(id(i,j))
     call kill(rt(i,j))
    enddo
    enddo
  end SUBROUTINE orthonormalisep

  SUBROUTINE furman_rrt(r,s,rt)
   implicit none
   type(real_8)  r(3,3),s(3,3),rt(3,3),id(3,3),ik(3,3)
   integer i,j,k

    do i=1,3
    do j=1,3
     call alloc(id(i,j))
     call alloc(ik(i,j))
    enddo
    enddo

      do i=1,3
       ik(i,i)=1.5e0_dp
      do j=1,3
      do k=1,3
       id(i,k)=r(i,j)*r(k,j)+id(i,k)
      enddo
      enddo
      enddo

    do i=1,3
    do j=1,3
     rt(i,j)=id(i,j)
     id(i,j)=0.e0_dp
    enddo
    enddo

    do i=1,3
    do j=1,3
    do k=1,3
      id(i,k)=(ik(i,j)-0.5e0_dp*rt(i,j))*r(j,k)+ id(i,k)
    enddo
    enddo
    enddo



    do i=1,3
    do j=1,3
     s(i,j)=id(i,j)
     call kill(id(i,j))
     call kill(ik(i,j))
    enddo
    enddo

end   SUBROUTINE furman_rrt

  SUBROUTINE track_TREE_G_complexr(T,XI)
    use da_arrays
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

  SUBROUTINE track_TREE_G_complexp(T,XI)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    type(real_8), INTENT(INOUT) :: XI(:)
    type(real_8) XT(lno),XF(lnv),XM(lno+1),XX
    INTEGER JC,I,IV

    call alloc(xt)
    call alloc(xf)
    call alloc(xm)
    call alloc(xx)

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

    call kill(xt)
    call kill(xf)
    call kill(xm)
    call kill(xx)

  END SUBROUTINE track_TREE_G_complexp


end module S_status
