!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis

module fitted_MAG
  USE S_status
  implicit none
  PRIVATE read_d,get_adata2c_p
  PRIVATE POINTERS_user1R,POINTERS_user1P
  logical(lp) :: extrapolate =.true.
  INTEGER :: NX_0=0,NS_0=0,NY_0=0,fitted_order=1
  PRIVATE ALLOC_user1,KILL_user1,copy_el_elp,copy_elP_el,copy_el_el
  PRIVATE ZEROR_user1,ZEROP_user1,INTR,INTP,INTS   !POINTERS_D,
  PRIVATE fx,px_to_xpr0_r,px_to_xpp0_r,px_to_xp_r
  PRIVATE xp_to_px_r,xp_to_pxr0_r,xp_to_pxp0_r
  !  PRIVATE track_interpol_1h_r,track_interpol_1h_p
  ! Vector potential stuff
  PRIVATE get_adata2c_r,fxhp,px_to_xpr0_hr,px_to_xpr0_hp,xp_to_pxr0_hr   ,fxhr
  !frs real(dp),private::eps_fitted=c_1d_11
  private fxhpm,fxhrm, che, get_adatag,get_bdatag,get_adata,b_fun,get_bdev
  logical(lp)::point_at=.true.,first_fitted=.true.,COPY_FIT=.TRUE.
  logical(lp) , target ::  check_iteration=.true., check_interpolate_x = .true., check_interpolate_y = .true.
  logical(lp) :: midplane_enforcing=.true.,maxwell_enforcing=.true.
  REAL(DP) :: CHECKSIZE = 1.D30
  integer :: fitted_iteration =50
  integer :: er=0
  TYPE DATA_MAGNET
     INTEGER, POINTER :: NS,NX,NY
     real(dp), POINTER :: dtheta,dx,dy
     real(dp), POINTER,DIMENSION(:)::x,y,theta
     REAL(dp),POINTER,DIMENSION(:,:,:,:,:)::B
     real(dp),POINTER,DIMENSION(:,:,:,:):: A
     logical(lp),pointer :: read_in
  END TYPE DATA_MAGNET

  type(data_magnet), target :: only_data

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE ZEROr_user1
     MODULE PROCEDURE ZEROp_user1
  END INTERFACE

  INTERFACE read
     MODULE PROCEDURE read_d
  END INTERFACE

  INTERFACE get_adata
     MODULE PROCEDURE get_adata2c_r
     MODULE PROCEDURE get_adata2c_p
  END INTERFACE

  TYPE FITTED_MAGNET
     TYPE(MAGNET_CHART), POINTER:: P
     real(dp), POINTER ::L
     real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
     real(dp), POINTER ::SCALE,h
     logical(lp), POINTER ::symplectic
     TYPE(DATA_MAGNET), POINTER:: D
     real(dp),  DIMENSION(:,:), POINTER :: X
     real(dp), POINTER ::xmin,xmax,ymin,ymax
  END  TYPE FITTED_MAGNET

  TYPE FITTED_MAGNETP
     TYPE(MAGNET_CHART), POINTER:: P
     TYPE(REAL_8), POINTER ::L
     TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
     logical(lp), POINTER ::symplectic
     !  SPECIAL POINTERS
     TYPE(REAL_8), POINTER ::SCALE,h
     TYPE(DATA_MAGNET), POINTER:: D
     real(dp), POINTER ::xmin,xmax,ymin,ymax
  END  TYPE FITTED_MAGNETP


  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_user1
  END INTERFACE


  INTERFACE KILL
     MODULE PROCEDURE KILL_user1
  END INTERFACE

  INTERFACE POINTERS_FITTED_P
     MODULE PROCEDURE POINTERS_user1R
     MODULE PROCEDURE POINTERS_user1P
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_el_elp
     MODULE PROCEDURE copy_elp_el
     MODULE PROCEDURE copy_el_el
     MODULE PROCEDURE copy_D_D
  END INTERFACE


  INTERFACE TRACK
     MODULE PROCEDURE INTR
     MODULE PROCEDURE INTP
     MODULE PROCEDURE INTS
  END INTERFACE

  INTERFACE track_interpol_1
     !     MODULE PROCEDURE track_interpol_1h_r
     !     MODULE PROCEDURE track_interpol_1h_p
     MODULE PROCEDURE track_interpol_1h_r_new
     MODULE PROCEDURE track_interpol_1h_p_new
  END INTERFACE

  INTERFACE fx
     MODULE PROCEDURE fxhr
     MODULE PROCEDURE fxhp
     MODULE PROCEDURE fxhrm
     MODULE PROCEDURE fxhpm
  END INTERFACE

  INTERFACE px_to_xp_r
     MODULE PROCEDURE px_to_xpr0_r     !r
     MODULE PROCEDURE px_to_xpp0_r    !p
     MODULE PROCEDURE px_to_xpr0_hr
     MODULE PROCEDURE px_to_xpr0_hp
  end  INTERFACE

  INTERFACE xp_to_px_r
     MODULE PROCEDURE xp_to_pxr0_r     !r
     MODULE PROCEDURE xp_to_pxp0_r    !p
     MODULE PROCEDURE xp_to_pxr0_hr
     MODULE PROCEDURE xp_to_pxr0_hp
  end  INTERFACE

CONTAINS

  !!   OPERATIONS ON TYPE DATA_MAGNET
  SUBROUTINE copy_D_d(D,DD)
    IMPLICIT NONE
    TYPE(DATA_MAGNET),pointer ::DD
    TYPE(DATA_MAGNET), INTENT(in)::D

    if(point_at) then
       dd=>only_data
    else
       DD%NS=D%NS;DD%NX=D%NX;DD%NY=D%NY; dd%read_in=d%read_in

       DD%dtheta=D%dtheta;DD%dx=D%dx;DD%dy=D%dy;

       DD%theta=D%theta;DD%x=D%x;DD%y=D%Y;

       IF(COPY_FIT) THEN
          IF(ASSOCIATED(D%B)) DD%B=D%B
          IF(ASSOCIATED(D%A)) DD%A=D%A
       ENDIF
    endif
  END SUBROUTINE copy_D_d


  SUBROUTINE POINTERS_D(D)
    IMPLICIT NONE
    TYPE(DATA_MAGNET), pointer ::D

    if(point_at) then
       d=>only_data
    else

       NULLIFY(D%NS,D%NX,D%NY,d%read_in)
       ALLOCATE(D%NS,D%NX,D%NY,d%read_in)
       D%NS=NS_0;D%NX=NX_0;D%NY=NY_0; d%read_in=.false.
       NULLIFY(D%dtheta,D%dx,D%dy)
       ALLOCATE(D%dtheta,D%dx,D%dy)
       D%dtheta=zero;D%dx=zero;D%dy=zero;
       NULLIFY(D%theta,D%x,D%y)
       ALLOCATE(D%theta(2),D%x(2),D%y(2))
       D%theta=zero;D%x=zero;D%y=zero;
       NULLIFY(D%B)
       IF(COPY_FIT) ALLOCATE(D%B(3,0:3,NS_0,NX_0,NY_0))
       NULLIFY(D%a)
       IF(COPY_FIT) ALLOCATE(D%a(13,NS_0,NX_0,NY_0))
    endif
  END SUBROUTINE POINTERS_D


  SUBROUTINE read_d(D,brho,filen)
    IMPLICIT NONE
    TYPE(DATA_MAGNET), INTENT(inout)::D
    real(dp), INTENT(in)::brho
    character(*)  filen
    INTEGER mf
    integer i,j,k
    real(dp) x1,x2,x3

    mf=NEWFILE
    open(unit=mf,file=filen)

    read(mf,*) d%ns,D%theta(1),d%theta(2),D%nx,D%x(1),D%x(2),d%ny,D%y(1),D%y(2)
    w_p=0
    w_p%nc=2
    w_p%fc='(1(1X,A120))'
    write(w_p%c(1),'(3(1x,i4,2(1x,g16.10)))') d%ns,D%theta(1),D%theta(2),D%nx,D%x(1),D%x(2),d%ny,D%y(1),D%y(2)
    d%dtheta=(D%theta(2)-D%theta(1))/(D%ns-1)
    d%dx=(D%x(2)-D%x(1))/dble(D%nx-1)
    d%dy=(D%y(2)-D%y(1))/dble(D%ny-1)
    write(w_p%c(2),'(3(1x,g16.10))') d%dtheta,d%dx,d%dy
    call WRITE_I

    do i=1,d%ns
       read(mf,*) x1,x2,x3
       do j=1,d%nx
          do k=1,d%ny
             d%a(6,i,j,k)=zero
             read(mf,*) x1,x2,d%a(1,i,j,k),d%a(2,i,j,k),d%a(3,i,j,k)
             d%a(1,i,j,k)=d%a(1,i,j,k)/brho
             d%a(2,i,j,k)=d%a(2,i,j,k)/brho
             d%a(3,i,j,k)=d%a(3,i,j,k)/brho
          enddo
       enddo
    enddo
    close(mf)

    !posx=sss
    do i=1,d%ns
       do j=1,d%nx
          do k=1,d%ny
             call get_adatag(d,i,j,k)
          enddo
       enddo
    enddo
    do i=1,d%ns
       do j=1,d%nx
          do k=1,d%ny
             call get_bdatag(d,i,j,k)
          enddo
       enddo
    enddo

    do i=1,d%ns
       do j=1,d%nx
          do k=1,d%ny
             call get_bdev(d,i,j,k)
          enddo
       enddo
    enddo
    d%read_in=.true.
    write(6,*) er, " border errors ", d%ns*d%nx*d%ny*3*4, " calls"

  END SUBROUTINE read_d

  SUBROUTINE read_n(filen)
    IMPLICIT NONE
    character(*)  filen
    INTEGER mf
    real(dp) x1,x2,x3,x4,x5,x6

    mf=NEWFILE
    open(unit=mf,file=filen)
    read(mf,*) ns_0,x1,x2,nx_0,x3,x4,ny_0,x5,x6
    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A120))'
    write(w_p%c(1),'(3(1x,i4,2(1x,g16.10)))') ns_0,x1,x2,nx_0,x3,x4,ny_0,x5,x6
    call write_i
    close(mf)
  END SUBROUTINE read_n



  !!   OPERATIONS ON TYPE FITTED_MAGNET

  SUBROUTINE ZEROR_user1(EL,I)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%SCALE))  THEN
          deallocate(el%xmin,el%xmax,el%ymin,el%ymax)
          DEALLOCATE(EL%SCALE)
          DEALLOCATE(EL%h)
          DEALLOCATE(EL%symplectic)
          if(.not.point_at) then
             DEALLOCATE(EL%D%NS,EL%D%NX,EL%D%NY,el%d%read_in)
             DEALLOCATE(EL%D%dtheta,EL%D%dx,EL%D%dY)
             DEALLOCATE(EL%D%theta,EL%D%x ,EL%D%Y)
             DEALLOCATE(EL%X)
             !          IF(ASSOCIATED(EL%D%B)) DEALLOCATE(EL%D%B)
             IF(ASSOCIATED(EL%D%A)) DEALLOCATE(EL%D%A)
             DEALLOCATE(EL%D)
          endif
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%symplectic)
       NULLIFY(EL%SCALE)
       NULLIFY(el%xmin,el%xmax,el%ymin,el%ymax)
       NULLIFY(EL%h)
       NULLIFY(EL%D)
       NULLIFY(EL%X)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZEROR_user1

  SUBROUTINE ZEROP_user1(EL,I)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%SCALE))  THEN
          CALL KILL(EL)
          deallocate(el%xmin,el%xmax,el%ymin,el%ymax)
          DEALLOCATE(EL%SCALE)
          DEALLOCATE(EL%h)
          DEALLOCATE(EL%symplectic)
          if(.not.point_at) then
             DEALLOCATE(EL%D%NS,EL%D%NX,EL%D%NY,el%d%read_in)
             DEALLOCATE(EL%D%dtheta,EL%D%dx,EL%D%dY)
             DEALLOCATE(EL%D%theta,EL%D%x ,EL%D%Y)
             IF(ASSOCIATED(EL%D%A)) DEALLOCATE(EL%D%A)
             !          IF(ASSOCIATED(EL%D%B)) DEALLOCATE(EL%D%B)
             DEALLOCATE(EL%D)
          endif
       ENDIF
    elseif(i==0)       then
       ! nullifies pointers
       NULLIFY(EL%symplectic)
       NULLIFY(EL%SCALE)
       NULLIFY(el%xmin,el%xmax,el%ymin,el%ymax)
       NULLIFY(EL%h)
       NULLIFY(EL%D)
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZEROP_user1

  SUBROUTINE ALLOC_user1(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL
    CALL ALLOC(EL%SCALE)
    CALL ALLOC(EL%h)
    ! ALLOC INTERNAL POLYMORPHS IF ANY
  END SUBROUTINE ALLOC_user1

  SUBROUTINE KILL_user1(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL

    CALL KILL(EL%SCALE)
    CALL KILL(EL%h)
    ! KILL INTERNAL POLYMORPHS IF ANY

  END SUBROUTINE KILL_user1

  SUBROUTINE POINTERS_user1R(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(INOUT)::EL

    ALLOCATE(EL%SCALE)
    ALLOCATE(el%xmin,el%xmax,el%ymin,el%ymax)
    ALLOCATE(EL%h)
    ALLOCATE(EL%symplectic)
    ALLOCATE(EL%D)
    ALLOCATE(EL%X(2,NS_0))
    CALL POINTERS_D(EL%D)
    EL%h=0.0_dp

    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_user1R

  SUBROUTINE POINTERS_user1P(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL

    ALLOCATE(EL%SCALE)
    ALLOCATE(el%xmin,el%xmax,el%ymin,el%ymax)
    ALLOCATE(EL%h)
    ALLOCATE(EL%symplectic)
    ALLOCATE(EL%D)
    CALL POINTERS_D(EL%D)
    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_user1P

  SUBROUTINE reset_FITTED_MAGNETP(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL


    CALL resetpoly_R31(EL%SCALE)
    CALL resetpoly_R31(EL%h)
    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS

  END SUBROUTINE reset_FITTED_MAGNETP


  SUBROUTINE copy_el_elp(EL,ELP)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(in)::EL
    TYPE(FITTED_MAGNETP), INTENT(inout)::ELP

    ELP%xmin    =EL%xmin
    ELP%xmax    =EL%xmax
    ELP%ymin    =EL%ymin
    ELP%ymax    =EL%ymax
    ELP%h    =EL%h
    ELP%SYMPLECTIC    =EL%SYMPLECTIC
    CALL COPY(EL%D,ELP%D)

    !  COPY CODING HERE NO ALLOCATION OF POINTERS OR POLYMORPH NEEDED
    !  IF DONE CORRECTLY

  END SUBROUTINE copy_el_elp

  SUBROUTINE copy_elp_el(EL,ELP)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(in)::EL
    TYPE(FITTED_MAGNET), INTENT(inout)::ELP



    ELP%xmin    =EL%xmin
    ELP%xmax    =EL%xmax
    ELP%ymin    =EL%ymin
    ELP%ymax    =EL%ymax

    ELP%SCALE    =EL%SCALE
    ELP%h    =EL%h
    ELP%SYMPLECTIC    =EL%SYMPLECTIC
    CALL COPY(EL%D,ELP%D)
    !  COPY CODING HERE NO ALLOCATION OF POINTERS OR POLYMORPH NEEDED
    !  IF DONE CORRECTLY


  END SUBROUTINE copy_elp_el

  SUBROUTINE copy_el_el(EL,ELP)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(in)::EL
    TYPE(FITTED_MAGNET), INTENT(inout)::ELP

    ELP%xmin    =EL%xmin
    ELP%xmax    =EL%xmax
    ELP%ymin    =EL%ymin
    ELP%ymax    =EL%ymax
    ELP%SCALE    =EL%SCALE
    ELP%h    =EL%h
    ELP%SYMPLECTIC    =EL%SYMPLECTIC
    CALL COPY(EL%D,ELP%D)
    !  COPY CODING HERE NO ALLOCATION OF POINTERS


  END SUBROUTINE copy_el_el

  !  TRACKING STUFF

  SUBROUTINE INTR(EL,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(FITTED_MAGNET),INTENT(INOUT):: EL

    select case(el%p%method)
    case(1)
       call track_interpol_1(el,x)
    case(2)
       call track_interpol_1(el,x)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A120))'
       w_p%c(1)= "METHOD NOT SUPPORTED IN track_interpol"
       call write_e(101)
    END SELECT


  END SUBROUTINE INTR






  SUBROUTINE INTP(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(FITTED_MAGNETP),INTENT(INOUT):: EL


    select case(el%p%method)
    case(1)
       call track_interpol_1(el,x)
    case(2)
       call track_interpol_1(el,x)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A120))'
       w_p%c(1)= "METHOD NOT SUPPORTED IN track_interpol"
       call write_e(101)
    END SELECT




  END SUBROUTINE INTP

  SUBROUTINE INTS(EL,X)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: X(6)
    TYPE(FITTED_MAGNETP),INTENT(IN):: EL


    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A120))'
    w_p%c(1)=  "USER1P NOT DEFINED "
    call write_e(113)

  END SUBROUTINE INTS







  integer function mod_s(i,j)
    implicit none
    integer, intent(in) :: i,j
    logical(lp) super
    super=.true.
    if(super) then
       mod_s=mod(i,j)
       if(mod_s==0) mod_s=j
    else
       if(i>=j) mod_s=j
    endif
  end function  mod_s


  subroutine px_to_xpr0_r(x,hc,p)         ! convert in drift
    implicit none
    type(MAGNET_CHART), intent(in) :: p
    real(dp), intent(inout):: x(6)
    real(dp) pz,d
    real(dp) beta0
    real(dp), intent(in):: hc

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif


    d=root(one+two*x(5)/beta0+x(5)**2)
    pz=root(d**2-x(2)**2-x(4)**2)

    x(2)=(one+hc*x(1))*x(2)/pz
    x(4)=(one+hc*x(1))*x(4)/pz

  end subroutine px_to_xpr0_r

  subroutine px_to_xpp0_r(x,hc,p)         ! convert in drift
    implicit none
    type(MAGNET_CHART), intent(in) :: p

    type(real_8), intent(inout):: x(6)
    type(real_8) pz,d
    real(dp), intent(in):: hc
    real(dp) beta0

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif

    call alloc(pz)
    call alloc(d)

    d=SQRT(one+two*x(5)/beta0+x(5)**2)
    pz=SQRT(d**2-x(2)**2-x(4)**2)

    x(2)=(one+hc*x(1))*x(2)/pz
    x(4)=(one+hc*x(1))*x(4)/pz

    call kill(pz)
    call kill(d)

  end subroutine px_to_xpp0_r

  subroutine xp_to_pxr0_r(x,hc,p)         ! convert in drift
    implicit none

    real(dp), intent(inout):: x(6)
    type(MAGNET_CHART), intent(in) :: p
    real(dp) pz,d
    real(dp) beta0
    real(dp), intent(in):: hc

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif

    d=root(one+two*x(5)/beta0+x(5)**2)
    pz=root((one+hc*x(1))**2+x(2)**2+x(4)**2)

    x(2)=d*x(2)/pz
    x(4)=d*x(4)/pz

  end subroutine xp_to_pxr0_r

  subroutine xp_to_pxp0_r(x,hc,p)         ! convert in drift
    implicit none

    type(real_8), intent(inout):: x(6)
    type(MAGNET_CHART), intent(in) :: p
    type(real_8) pz,d
    real(dp) beta0
    real(dp), intent(in):: hc

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif

    call alloc(pz)
    call alloc(d)

    d=SQRT(one+two*x(5)/beta0+x(5)**2)
    pz=SQRT((one+hc*x(1))**2+x(2)**2+x(4)**2)

    x(2)=d*x(2)/pz
    x(4)=d*x(4)/pz

    call kill(pz)
    call kill(d)

  end subroutine xp_to_pxp0_r


  ! Hamiltonian shit


  subroutine fxhr(f,x,a,az,ax,r,hcurv,p)
    implicit none

    real(dp)  d(2),BETA0,GAMMA0I
    real(dp) ,intent(in) :: a(3),az(3),ax(3),r ,hcurv
    type(MAGNET_CHART), intent(inout) :: p
    real(dp) ,intent(inout) :: x(6)
    real(dp), intent(out):: f(6)

    if(p%time) then
       beta0=p%beta0;GAMMA0I=p%GAMMA0I;
    else
       beta0=one;GAMMA0I=zero;
    endif
    x(1)=x(1)-r


    d(1)=(one+hcurv*x(1))
    d(2)=(one+2*x(5)/beta0+x(5)**2-(x(2)-a(1))**2-(x(4)-a(2))**2)
    d(2)=root(d(2))
    f(1)=d(1)*(x(2)-a(1))/d(2)
    f(3)=d(1)*(x(4)-a(2))/d(2)

    f(2)=hcurv*d(2)+(f(1)*ax(1))+az(1)
    f(4)=(f(1)*ax(2))+az(2)
    f(5)=0
    f(6)=(one/beta0+x(5))*d(1)/d(2)


    x(1)=x(1)+r
  end subroutine fxhr




  subroutine fxhrm(f,m,x,a,az,ax,r,hcurv,p)
    implicit none

    real(dp)  d(2),BETA0,GAMMA0I
    real(dp) ,intent(in) :: a(3),az(3),ax(3),r ,hcurv
    type(MAGNET_CHART), intent(in) :: p
    real(dp) ,intent(inout) :: x(6)
    real(dp), intent(out):: f(6),m(2,2)

    if(p%time) then
       beta0=p%beta0;GAMMA0I=p%GAMMA0I;
    else
       beta0=one;GAMMA0I=zero;
    endif
    x(1)=x(1)-r


    d(1)=(one+hcurv*x(1))
    d(2)=root(one+2*x(5)/beta0+x(5)**2-(x(2)-a(1))**2-(x(4)-a(2))**2)

    f(1)=d(1)*(x(2)-a(1))/d(2)
    f(3)=d(1)*(x(4)-a(2))/d(2)
    f(2)=hcurv*d(2)+(f(1)*ax(1))+az(1)
    f(4)=(f(1)*ax(2))+az(2)
    f(5)=0
    f(6)=(one/beta0+x(5))*d(1)/d(2)

    m(1,1)=d(1)*(one+(x(2)-a(1))/d(2)**2)/d(2)
    m(2,1)=m(1,1)*ax(2)
    m(1,1)=m(1,1)*ax(1)-hcurv*(x(2)-a(1))/d(2)
    m(2,2)=d(1)*(x(2)-a(1))*x(4)/d(2)**3
    m(1,2)=-hcurv*ax(2)/d(2)+m(2,2)*ax(2)
    m(2,2)=m(2,2)*ax(2)

    x(1)=x(1)+r
  end subroutine fxhrm

  subroutine fxhp(f,x,a,az,ax,r,hcurv,p)
    implicit none

    real(dp) ,intent(in) :: r ,hcurv
    type(real_8) ,intent(in) :: a(3),az(3),ax(3)
    type(MAGNET_CHART), intent(in) :: p
    type(real_8),intent(inout) :: x(6)
    type(real_8), intent(out):: f(6)
    type(real_8)   d(2)
    real(dp)  BETA0,GAMMA0I

    call alloc(d,2)


    if(p%time) then
       beta0=p%beta0;GAMMA0I=p%GAMMA0I;
    else
       beta0=one;GAMMA0I=zero;
    endif
    x(1)=x(1)-r


    d(1)=(one+hcurv*x(1))
    d(2)=SQRT(one+2*x(5)/beta0+x(5)**2-(x(2)-a(1))**2-(x(4)-a(2))**2)

    f(1)=d(1)*(x(2)-a(1))/d(2)
    f(3)=d(1)*(x(4)-a(2))/d(2)
    f(2)=hcurv*d(2)+(f(1)*ax(1))+az(1)
    f(4)=(f(1)*ax(2))+az(2)
    f(5)=0
    f(6)=(one/beta0+x(5))*d(1)/d(2)


    x(1)=x(1)+r

    call kill(d,2)
  end subroutine fxhp

  subroutine fxhpm(f,m,x,a,az,ax,r,hcurv,p)
    implicit none

    real(dp) ,intent(in) :: r ,hcurv
    type(real_8) ,intent(in) :: a(3),az(3),ax(3)
    type(MAGNET_CHART), intent(in) :: p
    type(real_8),intent(inout) :: x(6)
    type(real_8), intent(out):: f(6)
    real(dp), intent(out):: m(2,2)
    type(real_8)   d(2)
    real(dp)  BETA0,GAMMA0I

    call alloc(d,2)


    if(p%time) then
       beta0=p%beta0;GAMMA0I=p%GAMMA0I;
    else
       beta0=one;GAMMA0I=zero;
    endif
    x(1)=x(1)-r


    d(1)=(one+hcurv*x(1))
    d(2)=SQRT(one+2*x(5)/beta0+x(5)**2-(x(2)-a(1))**2-(x(4)-a(2))**2)

    f(1)=d(1)*(x(2)-a(1))/d(2)
    f(3)=d(1)*(x(4)-a(2))/d(2)
    f(2)=hcurv*d(2)+(f(1)*ax(1))+az(1)
    f(4)=(f(1)*ax(2))+az(2)
    f(5)=0
    f(6)=(one/beta0+x(5))*d(1)/d(2)

    m(1,1)=d(1)*(one+(x(2)-a(1))/d(2)**2)/d(2)
    m(2,1)=m(1,1)*ax(2)
    m(1,1)=m(1,1)*ax(1)-hcurv*(x(2)-a(1))/d(2)
    m(2,2)=d(1)*(x(2)-a(1))*x(4)/d(2)**3
    m(1,2)=-hcurv*ax(2)/d(2)+m(2,2)*ax(2)
    m(2,2)=m(2,2)*ax(2)

    x(1)=x(1)+r

    call kill(d,2)
  end subroutine fxhpm

  subroutine track_interpol_1h_r_new(bend,x)
    implicit none
    type(fitted_magnet), intent(inOUT):: bend
    real(dp) , intent(inout)::x(6)
    real(dp) k1(6),xt(6),hc,rho,dz,DAL,ds
    real(dp) af(3),az(3),ax(3),m(2,2),mi(2,2),dif(2),t(2),det,bf(3),XI,ZETA
    integer  i,j,k,NF,KNF,it,ia
    real(dp) norm, norm0,rhob0
    logical(lp) doit
    logical(lp):: DONEITT
    real(dp) bet,eps,gam,alph,rg,theta_half
    IF(.NOT.CHECK_STABLE) return
    DONEITT=.true.
    theta_half=DEG_TO_RAD_*(bend%d%theta(2)-bend%d%theta(1))/2.0_dp
    !  x(1)=x(1)+one/el%p%b0-el%d%x(1)

    alph=bend%p%b0*bend%p%ld/2.0_dp

    gam=atan(sin(alph)/(cos(alph)+bend%p%b0*bend%h)    )
    eps=((bend%d%theta(2)-bend%d%theta(1))*DEG_TO_RAD_)/2.0_dp-gam
    eps=-eps
    bet=alph-gam
    rg=(cos(alph)+bend%p%b0*bend%h)/cos(gam)/bend%p%b0

    call ROT_XZ(bet,X,bend%p%BETA0,DONEITT,bend%p%TIME)

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)
    rhob0=one/bend%p%b0
    !  rotating
    nf=bend%p%nst
    x(1)=x(1)+rg
    call ROT_XZ(eps,X,bend%p%BETA0,DONEITT,bend%p%TIME)
    !     x(1)=x(1)-rg
    !    x(1)=x(1)+rho

    ds=(bend%d%dtheta*DEG_TO_RAD_)*rho
    DAL=bend%d%dtheta*DEG_TO_RAD_
    ia=1
    if(bend%symplectic) ia=0
    !tot=zero
    do i=1,bend%d%ns    -ia
       dz=ds/nf
       if((i==1.or.i==bend%d%ns).and.bend%symplectic) dz=dz/2.0_dp
       DO KNF=1,NF
          XI  =(X(1))*COS((I-1)*DAL-theta_half)-bend%h
          ZETA=(X(1))*SIN((I-1)*DAL-theta_half)
          BEND%X(1,I)=ATAN(ZETA/XI)
          BEND%X(2,I)=SQRT(XI**2+ZETA**2)-RHOB0

          !(pos,b,i,x,y,af,az,ax)

          call get_adata(bend%d,i,x(1),x(3),af,az,ax,bf,bend%scale)
          ! original formula for electrons (q<0)  now switch to proton


          ! fxhr(f,x,a,az,ax,r,hcurv,p)
          call fx(k1,x,af,az,ax,rho,hc,bend%p)
          if(.not.CHECK_STABLE) return

          if(bend%symplectic) then
             do j=1,3
                xt(2*j)=x(2*j)+dz*k1(2*j)         ! temporary
                xt(2*j-1)=x(2*j-1)
             enddo


             norm0=c_1d10
             doit=.true.

             do it=1,fitted_iteration
                call fx(k1,m,xt,af,az,ax,rho,hc,bend%p)
                if(.not.check_stable) return
                dif(1)=x(2)+dz*k1(2)-xt(2)
                dif(2)=x(4)+dz*k1(4)-xt(4)

                m=-dz*m;
                m(1,1)=one+m(1,1);m(2,2)=one+m(2,2);
                det=m(1,1)*m(2,2)-m(2,1)*m(1,2)
                m=m/det;
                mi(1,1)=m(2,2);mi(2,2)=m(1,1);mi(1,2)=-m(1,2);mi(2,1)=-m(2,1);
                t(1)=mi(1,1)*dif(1)+mi(1,2)*dif(2)
                t(2)=mi(2,1)*dif(1)+mi(2,2)*dif(2)
                xt(2)=xt(2)+t(1);    xt(4)=xt(4)+t(2);
                norm=abs(dif(1))+abs(dif(2))
                if(norm>eps_fitted.and.doit) then
                   norm0=norm
                else
                   doit=.false.
                   if(norm>=norm0) goto 100
                   norm0=norm
                endif
             enddo
             !             w_p=0
             !             w_p%nc=1
             !             w_p%fc='(1(1X,A72))'
             !             w_p%c(1)=  " Took more than 20 iterations "
             !             call write_e(100)
             check_iteration=.false.
             check_stable=.false.
             IF(.NOT.CHECK_STABLE) return
100          continue

             call fx(k1,xt,af,az,ax,rho,hc,bend%p)
             if(.not.check_stable) return

             x(5)=x(5)            ! temporary
             x(6)=x(6)+dz*k1(6)    ! temporary
             do j=1,2
                x(2*j-1)=x(2*j-1)+dz*k1(2*j-1)         ! temporary
                x(2*j) =xt(2*j)
             enddo
          else
             if(bend%p%method==2) then
                do j=1,6
                   xt(j)=x(j)+dz*k1(j)         ! temporary
                enddo
                do j=1,6
                   x(j)=x(j)+dz*k1(j)/2.d0         ! temporary
                enddo
                call get_adata(bend%d,i+1,xt(1),xt(3),af,az,ax,bf,bend%scale)

                call fx(k1,xt,af,az,ax,rho,hc,bend%p)
                do j=1,6
                   x(j)=x(j)+dz*k1(j)/2.d0         ! temporary
                enddo
             else

                do j=1,6
                   x(j)=x(j)+dz*k1(j)         ! temporary
                enddo

             endif
          endif

          if(che(x,bend)) then
             check_stable=.false.
             return
          endif
          !        tot=tot+bend%dz


       ENDDO
    enddo
    XI  =(X(1))*COS((bend%d%ns-1)*DAL-theta_half)-bend%h
    ZETA=(X(1))*SIN((bend%d%ns-1)*DAL-theta_half)
    BEND%X(1,bend%d%ns)=ATAN(ZETA/XI)
    BEND%X(2,bend%d%ns)=SQRT(XI**2+ZETA**2)-RHOB0

    !write(6,*) tot
    i=bend%d%ns-1

    !    x(1)=x(1)-rho


    !     x(1)=x(1)+rg
    call ROT_XZ(eps,X,bend%p%BETA0,DONEITT,bend%p%TIME)
    x(1)=x(1)-rg
    call ROT_XZ(bet,X,bend%p%BETA0,DONEITT,bend%p%TIME)

    if(bend%P%TIME) then
       X(6)=X(6)-(1-bend%P%TOTALPATH)*bend%P%LD/bend%P%BETA0
    else
       X(6)=X(6)-(1-bend%P%TOTALPATH)*bend%P%LD
    endif

  end subroutine track_interpol_1h_r_new

  subroutine track_interpol_1h_p_new(bend,x)
    implicit none
    type(fitted_magnetp), intent(inOUT):: bend
    type(real_8), intent(inout)::x(6)
    type(real_8) k1(6),af(3),az(3),ax(3),dif(2),t(2),xt(6),bf(3)
    real(dp) hc,rho,dz,ds
    real(dp) m(2,2),mi(2,2),det
    integer  i,j,k,NF,KNF,it,ia
    real(dp) norm, norm0
    logical(lp) doit
    logical(lp):: DONEITT
    real(dp) alph
    TYPE(REAL_8) bet,eps,gam,rg
    DONEITT=.true.
    IF(.NOT.CHECK_STABLE) return

    call alloc(bet,eps,gam,rg)
    call alloc(k1,6);call alloc(af,3);call alloc(az,3);call alloc(ax,3);call alloc(bf,3);
    call alloc(dif,2);call alloc(t,2);call alloc(xt,6);

    alph=bend%p%b0*bend%p%ld/2.0_dp

    gam=atan(sin(alph)/(cos(alph)+bend%p%b0*bend%h)    )
    eps=((bend%d%theta(2)-bend%d%theta(1))*DEG_TO_RAD_)/2.0_dp-gam
    eps=-eps
    bet=alph-gam
    rg=(cos(alph)+bend%p%b0*bend%h)/cos(gam)/bend%p%b0



    call ROT_XZ(bet,X,bend%p%BETA0,DONEITT,bend%p%TIME)

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)
    !    rhob0=one/bend%p%b0
    !  rotating
    nf=bend%p%nst
    x(1)=x(1)+rg
    call ROT_XZ(eps,X,bend%p%BETA0,DONEITT,bend%p%TIME)
    !     x(1)=x(1)-rg
    !    x(1)=x(1)+rho



    ds=(bend%d%dtheta*DEG_TO_RAD_)*rho

    ia=1
    if(bend%symplectic) ia=0

    !tot=zero
    do i=1,bend%d%ns   -ia
       dz=ds/nf
       if((i==1.or.i==bend%d%ns).and.bend%symplectic) dz=dz/2.0_dp
       DO KNF=1,NF

          call get_adata(bend%d,i,x(1),x(3),af,az,ax,bf,bend%scale)


          call fx(k1,x,af,az,ax,rho,hc,bend%p)

          if(bend%symplectic) then

             do j=1,3
                xt(2*j)=x(2*j)+dz*k1(2*j)         ! temporary
                xt(2*j-1)=x(2*j-1)
             enddo


             norm0=c_1d10
             doit=.true.

             do it=1,fitted_iteration
                call fx(k1,m,xt,af,az,ax,rho,hc,bend%p)
                dif(1)=x(2)+dz*k1(2)-xt(2)
                dif(2)=x(4)+dz*k1(4)-xt(4)

                m=-dz*m;
                m(1,1)=one+m(1,1);m(2,2)=one+m(2,2);
                det=m(1,1)*m(2,2)-m(2,1)*m(1,2)
                m=m/det;
                mi(1,1)=m(2,2);mi(2,2)=m(1,1);mi(1,2)=-m(1,2);mi(2,1)=-m(2,1);
                t(1)=mi(1,1)*dif(1)+mi(1,2)*dif(2)
                t(2)=mi(2,1)*dif(1)+mi(2,2)*dif(2)
                xt(2)=xt(2)+t(1);    xt(4)=xt(4)+t(2);
                norm=abs(dif(1))+abs(dif(2))
                if(norm>eps_fitted.and.doit) then
                   norm0=norm
                else
                   doit=.false.
                   if(norm>=norm0) goto 100
                   norm0=norm
                endif
             enddo
             !             w_p=0
             !             w_p%nc=1
             !             w_p%fc='(1(1X,A72))'
             !             w_p%c(1)=  " Took more than 20 iterations "
             !             call write_e(100)
             check_iteration=.false.
             check_stable=.false.
100          continue

             call fx(k1,xt,af,az,ax,rho,hc,bend%p)

             x(5)=x(5)            ! temporary
             x(6)=x(6)+dz*k1(6)    ! temporary
             do j=1,2
                x(2*j-1)=x(2*j-1)+dz*k1(2*j-1)         ! temporary
                x(2*j) =xt(2*j)
             enddo


          else
             if(bend%p%method==2) then
                do j=1,6
                   xt(j)=x(j)+dz*k1(j)         ! temporary
                enddo
                do j=1,6
                   x(j)=x(j)+dz*k1(j)/2.d0         ! temporary
                enddo
                call get_adata(bend%d,i+1,xt(1 ),xt(3),af,az,ax,bf,bend%scale)

                call fx(k1,xt,af,az,ax,rho,hc,bend%p)
                do j=1,6
                   x(j)=x(j)+dz*k1(j)/2.d0         ! temporary
                enddo
             else
                do j=1,6
                   x(j)=x(j)+dz*k1(j)         ! temporary
                enddo

             endif
          endif


          !        tot=tot+bend%dz

       ENDDO
    enddo

    !write(6,*) tot
    i=bend%d%ns-1

    !    x(1)=x(1)-rho

    !     x(1)=x(1)+rg
    call ROT_XZ(eps,X,bend%p%BETA0,DONEITT,bend%p%TIME)
    x(1)=x(1)-rg
    call ROT_XZ(bet,X,bend%p%BETA0,DONEITT,bend%p%TIME)


    if(bend%P%TIME) then
       X(6)=X(6)-(1-bend%P%TOTALPATH)*bend%P%LD/bend%P%BETA0
    else
       X(6)=X(6)-(1-bend%P%TOTALPATH)*bend%P%LD
    endif


    call kill(bet,eps,gam,rg)

    call kill(dif,2);call kill(t,2);call kill(xt,6);
    call kill(k1,6);call kill(af,3);call kill(az,3);call kill(ax,3);call kill(bf,3);

  end subroutine track_interpol_1h_p_new

  subroutine get_adatag(b,i,posx,posy)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i,posx,posy
    integer j,ix,iy
    real(dp) db(3,3),x0,y0,rho
    real(dp) rm,rp,r0,b0(3)
    !  second order
    real(dp) db2(3,3)

    if(b%read_in) then
       write(6,*) " This is forbidden since already the data was read in for the fitted magnet"
       stop 222
    endif
    db(:,:)=0.0_dp
    db2(:,:)=0.0_dp
    rho=b%X(1)
    !  posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    !  posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    !  interpolate
    ix=(posx)
    iy=(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
    !write(6,*) "R: x0,y0",x0,y0
    rm=(x0-b%dx)/rho
    rp=(x0+b%dx)/rho
    !    rm2=(x0-b%dx)/rho
    !    rp2=(x0+b%dx)/rho
    r0=x0/rho
    !  first order
    if(fitted_order==1) then
       if(posx<2.or.posx>b%nx-1) then
          if(extrapolate) then
             ix=2
             if(posx>2) ix=b%nx-2
             x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A72))'
             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " ix,b%nx ",ix,b%nx
             call write_e(2001)
          endif
       endif

       if(posy<2.or.posy>b%ny-1) then
          if(extrapolate) then
             iy=2
             if(posy>2) iy=b%ny-2
             y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A72))'
             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " iy,b%ny ",iy,b%ny
             call write_e(2002)
          endif
       endif

       do j=1,2
          b0(j)=r0*b%a(j,i,ix,iy)
       enddo
       j=3
       b0(j)=b%a(j,i,ix,iy)


       do j=1,2
          db(j,1)=(rp*b%a(j,i,ix+1,iy)-rm*b%a(j,i,ix-1,iy))/b%dx/2
          db(j,2)=r0*(b%a(j,i,ix,iy+1)-b%a(j,i,ix,iy-1))/b%dy/2
       enddo
       j=3
       db(j,1)=(b%a(j,i,ix+1,iy)-b%a(j,i,ix-1,iy))/b%dx/2
       db(j,2)=(b%a(j,i,ix,iy+1)-b%a(j,i,ix,iy-1))/b%dy/2

    endif
    if(fitted_order==2) then
       ! second order
       if(posx<2.or.posx>b%nx-1) then
          if(extrapolate) then
             ix=2
             if(posx>2) ix=b%nx-2
             x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A120))'
             write(w_p%c(1),'(a9,2(1x,i4))') " ix,b%nx ",ix,b%nx
             call write_e(2003)
          endif
       endif

       if(posy<2.or.posy>b%ny-1) then
          if(extrapolate) then
             iy=2
             if(posy>2) iy=b%ny-2
             y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A120))'
             write(w_p%c(1),'(a9,2(1x,i4))') " iy,b%ny ",iy,b%ny
             call write_e(2004)
          endif
       endif

       do j=1,2
          b0(j)=r0*b%a(j,i,ix,iy)
       enddo
       j=3
       b0(j)=b%a(j,i,ix,iy)


       do j=1,2
          db(j,1)=(rp*b%a(j,i,ix+1,iy)-rm*b%a(j,i,ix-1,iy))/b%dx/2
          db(j,2)=r0*(b%a(j,i,ix,iy+1)-b%a(j,i,ix,iy-1))/b%dy/2
       enddo
       j=3
       db(j,1)=(b%a(j,i,ix+1,iy)-b%a(j,i,ix-1,iy))/b%dx/2
       db(j,2)=(b%a(j,i,ix,iy+1)-b%a(j,i,ix,iy-1))/b%dy/2


       do j=1,2
          db2(j,1)=(rp*b%a(j,i,ix+1,iy)+rm*b%a(j,i,ix-1,iy)-r0*2*b%a(j,i,ix,iy))/b%dx**2/2
          db2(j,2)=r0*(b%a(j,i,ix,iy+1)+b%a(j,i,ix,iy-1)-2*b%a(j,i,ix,iy))/b%dy**2/2
       enddo
       j=3
       db2(j,1)=(b%a(j,i,ix+1,iy)+b%a(j,i,ix-1,iy)-2*b%a(j,i,ix,iy))/b%dx**2/2
       db2(j,2)=(b%a(j,i,ix,iy+1)+b%a(j,i,ix,iy-1)-2*b%a(j,i,ix,iy))/b%dy**2/2


       !       do j=1,2
       !          db2(j,3)=(rp*b%a(j,i,ix+1,iy+1)+rm*b%a(j,i,ix-1,iy-1)-r0*2*b%a(j,i,ix,iy))/2  ! this is asymetric
       !          db2(j,3)=(db2(j,3)-db2(j,1)*b%dx**2-db2(j,2)*b%dy**2)/b%dx/b%dy
       !       enddo


       do j=1,2
          db2(j,3)=rp*(b%a(j,i,ix+1,iy+1)-b%a(j,i,ix+1,iy-1))/2+rm*(b%a(j,i,ix-1,iy-1)-b%a(j,i,ix-1,iy+1))/2  ! this is asymetric
          db2(j,3)=(db2(j,3))/b%dx/b%dy/2
       enddo





       j=3
       !          db2(j,3)=(b%a(j,i,ix+1,iy+1)+b%a(j,i,ix-1,iy-1)-2*b%a(j,i,ix,iy))/2
       !          db2(j,3)=(db2(j,3)-db2(j,1)*b%dx**2-db2(j,2)*b%dy**2)/b%dx/b%dy
       db2(j,3)=(b%a(j,i,ix+1,iy+1)-b%a(j,i,ix+1,iy-1))/2+(b%a(j,i,ix-1,iy-1)-b%a(j,i,ix-1,iy+1))/2  ! this is asymetric
       db2(j,3)=(db2(j,3))/b%dx/b%dy/2

       if(posy<3.or.posy>b%ny-2) then
          if(extrapolate) then
             iy=3
             if(posy>3) iy=b%ny-3
             y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A72))'
             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " iy,b%ny ",iy,b%ny
             call write_e(2005)
          endif
       endif
       !       j=3
       !        write(6,*) "here",db(j,2)
       !       db(j,2)=( 8*(b%a(j,i,ix,iy+1)-b%a(j,i,ix,iy-1))-(b%a(j,i,ix,iy+2)-b%a(j,i,ix,iy-2)) )/b%dy/12
       !        write(6,*) db(j,2)
    endif



    ! Maxwell's equation
    ! insuring curl=0
    if(maxwell_enforcing) then
       rm=db(2,1)-db(1,2)-b0(2)/x0
    else
       rm=0.d0
    endif
    db(1,2)=db(1,2)+rm/2.d0
    db(2,1)=db(2,1)-rm/2.d0
    db(2,3)=db(3,2)*r0
    db(1,3)=db(3,1)*r0

    if(maxwell_enforcing) then
       rm=db2(2,3)-2*db2(1,2)-db(2,2)/x0
    else
       rm=0.d0
    endif
    db2(1,2)=db2(1,2)+rm/4.d0
    db2(2,3)=db2(2,3)-rm/2.d0

    if(maxwell_enforcing) then
       rm=2*db2(2,1)-db2(1,3)-db(1,2)/x0
    else
       rm=0.d0
    endif
    db2(2,1)=db2(2,1)-rm/4.d0
    db2(1,3)=db2(1,3)+rm/2.d0

    !insuring div =0
    db(3,3)=-db(1,1)-db(2,2)

    ! vector potential  ( (1+hx)bx ,(1+hx)by , bz )

    !   b%a(1,i,ix,iy)=b0(1)
    !   b%a(2,i,ix,iy)=b0(2)
    !   b%a(3,i,ix,iy)=b0(3)

    b%a(4,i,ix,iy)=db(1,1)
    b%a(5,i,ix,iy)=db(1,2)
    b%a(6,i,ix,iy)=db(3,1)
    b%a(7,i,ix,iy)=db(3,2)


    b%a(8,i,ix,iy)=db2(3,1)
    b%a(9,i,ix,iy)=db2(3,2)
    b%a(10,i,ix,iy)=db2(1,1)
    b%a(11,i,ix,iy)=db2(1,2)
    b%a(12,i,ix,iy)=db2(1,3)
    b%a(13,i,ix,iy)=db2(3,3)

    DO J=4,7
       IF(ABS(b%a(J,i,ix,iy))>CHECKSIZE) THEN
          WRITE(6,*) J,I,IX,IY
          WRITE(6,*) b%a(J,i,ix,iy)
          STOP 999
       ENDIF
    ENDDO


  end subroutine get_adatag

  subroutine get_bdatag(b,i,posx,posy)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i,posx,posy
    integer j,ix,iy,co(3),pos(3),k
    real(dp) r0,b0(3),x0,y0,rho

    if(b%read_in) then
       write(6,*) " This is forbidden since already the data was read in for the fitted magnet"
       stop 222
    endif
    rho=b%X(1)
    !  posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    !  posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    !  interpolate
    ix=(posx)
    iy=(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)

    r0=x0/rho

    do j=1,2
       b0(j)=r0*b%a(j,i,ix,iy)
    enddo
    j=3
    b0(j)=b%a(j,i,ix,iy)

    b%b(1,0,i,ix,iy)= b%a(1,i,ix,iy)
    b%b(2,0,i,ix,iy)= b%a(2,i,ix,iy)
    b%b(3,0,i,ix,iy)= b%a(3,i,ix,iy)



    ! vector potential  ( (1+hx)bx ,(1+hx)by , bz )

    b%a(1,i,ix,iy)=b0(1)
    b%a(2,i,ix,iy)=b0(2)
    b%a(3,i,ix,iy)=b0(3)

    pos(1)=ix
    pos(2)=iy
    pos(3)=i
    co=0
    do j=1,3
       do k=1,3
          co(1)=k
          b%b(j,k,i,ix,iy)=dev_f(b,pos,j,co,1.d0)
       enddo
    enddo



  end subroutine get_bdatag



  subroutine get_bdev(b,i,posx,posy)

    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i,posx,posy
    integer j,ix,iy,co(3),pos(3),k

    ix=(posx)
    iy=(posy)
    pos(1)=ix
    pos(2)=iy
    pos(3)=i
    co=0
    do j=1,3
       do k=1,3
          co(1)=k
          b%b(j,k,i,ix,iy)=dev_f(b,pos,j,co,1.d0)
       enddo
    enddo

  end subroutine get_bdev


  subroutine get_adata2c_r(b,i,x,y,af,az,ax,bf,scale)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i
    real(dp), intent(in) :: x,y,scale
    real(dp), intent(out) :: af(3),az(3),ax(3),bf(3)
    integer ix,iy,k
    real(dp) posx,posy,x0,y0,rho,dx,dy

    rho=b%X(1)
    posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one
    !
    if(posy<(b%ny-1)/2+1.and.midplane_enforcing) posy=posy+one

    !  interpolate
    ix=int(posx)
    iy=int(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
    !    write(6,*) "R: x0,y0",x0,y0
    !    write(6,*) "posx, posy",posx, posy


    !  first order

    if(fitted_order/=1.and.fitted_order/=2) stop 889
    if(posx<2.or.posx>b%nx-1) then
       if(extrapolate) then
          ix=2
          if(posx>2) ix=b%nx-2
          x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
       else
          !             w_p=0
          !             w_p%nc=1
          !             w_p%fc='(1(1X,A72))'
          !             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " ix,b%nx ",ix,b%nx
          !             call write_e(2234)
          write(6,*) " check_interpolate_x ",x
          check_interpolate_x=.false.
          check_stable=.false.
       endif
    endif

    if(posy<2.or.posy>b%ny-1) then
       if(extrapolate) then
          iy=2
          if(posy>2) iy=b%ny-2
          y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
       else
          !             w_p=0
          !             w_p%nc=1
          !             w_p%fc='(1(1X,A72))'
          !             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " iy,b%ny ",iy,b%ny
          !             call write_e(2235)
          write(6,*) " check_interpolate_y ",y
          check_interpolate_y=.false.
          check_stable=.false.
       endif
    endif

    !az(n) grad az
    !af(n)  a
    !ax(n)  grad ax

    !    write(6,*) size(b%a,1),size(b%a,2),size(b%a,3),size(b%a,4)
    !    write(6,*) i,ix,iy
    !     pause 333
    az=zero;af=zero;ax=zero;
    dx=x-x0
    dy=y-y0

    az(1)=-b%a(2,i,ix,iy)+b%a(4,i,ix,iy)*dy-b%a(5,i,ix,iy)*dx
    az(2)= b%a(1,i,ix,iy)+b%a(4,i,ix,iy)*dx+b%a(5,i,ix,iy)*dy
    af(1)=-b%a(3,i,ix,iy)*dy-b%a(6,i,ix,iy)*dx*dy-b%a(7,i,ix,iy)*dy**2/two
    ax(1)=-b%a(6,i,ix,iy)*dy
    ax(2)=-b%a(3,i,ix,iy)-b%a(6,i,ix,iy)*dx-b%a(7,i,ix,iy)*dy



    af(1)=af(1)-b%a(8,i,ix,iy)*dx**2*dy-b%a(9,i,ix,iy)*dy**3/3.d0    -b%a(13,i,ix,iy)*dx*dy**2/2.d0
    ax(1)=ax(1)-2.d0*b%a(8,i,ix,iy)*dx*dy                            -b%a(13,i,ix,iy)*dy**2/2.d0
    ax(2)=ax(2)-b%a(8,i,ix,iy)*dx**2-b%a(9,i,ix,iy)*dy**2   -b%a(13,i,ix,iy)*dx*dy

    az(1)=az(1)+b%a(10,i,ix,iy)*2.d0*dx*dy-b%a(12,i,ix,iy)*(dx**2-dy**2)/2.d0-b%a(5,i,ix,iy)*dx**2/2.0_dp/x0
    az(2)=az(2)+b%a(10,i,ix,iy)*dx**2+b%a(11,i,ix,iy)*dy**2+b%a(12,i,ix,iy)*dx*dy



    bf=0.d0
    do k=1,3
       af(k)=-scale*af(k)
       az(k)=-scale*az(k)
       ax(k)=-scale*ax(k)
       bf(k)=-scale*bf(k)
    enddo

  end subroutine get_adata2c_r

  subroutine get_adata2c_p(b,i,xx,yy,af,az,ax,bf,scale)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i
    type(real_8), intent(in) :: xx,yy,scale
    type(real_8), intent(out) :: af(3),az(3),ax(3),bf(3)
    integer ix,iy,k
    real(dp) posx,posy,x0,y0,x,y,rho
    type(real_8) dx,dy
    !  second order
    call alloc(dx,dy)

    rho=b%X(1)
    x=xx
    y=yy
    posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    if(posy<(b%ny-1)/2+1.and.midplane_enforcing) posy=posy+one

    !  interpolate
    ix=int(posx)
    iy=int(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
    !write(6,*) "R: x0,y0",x0,y0

    !  first order

    if(fitted_order/=1.and.fitted_order/=2) stop 890
    if(posx<2.or.posx>b%nx-1) then
       if(extrapolate) then
          ix=2
          if(posx>2) ix=b%nx-2
          x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
       else
          !             w_p=0
          !             w_p%nc=1
          !             w_p%fc='(1(1X,A72))'
          !             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " ix,b%nx ",ix,b%nx
          !             call write_e(2234)
          check_interpolate_x=.false.
          check_stable=.false.
       endif
    endif

    if(posy<2.or.posy>b%ny-1) then
       if(extrapolate) then
          iy=2
          if(posy>2) iy=b%ny-2
          y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
       else
          !             w_p=0
          !             w_p%nc=1
          !             w_p%fc='(1(1X,A72))'
          !             write(w_p%c(1),'(A9,1X,I4,1X,I4)') " iy,b%ny ",iy,b%ny
          !             call write_e(2234)
          check_interpolate_y=.false.
          check_stable=.false.
       endif
    endif


    !az(n) grad az
    !af(n)  a
    !ax(n)  grad ax
    af(1)=zero; af(2)=zero; af(3)=zero;
    az(1)=zero; az(2)=zero; az(3)=zero;
    ax(1)=zero; ax(2)=zero; ax(3)=zero;
    dx=xx-x0
    dy=yy-y0

    az(1)=-b%a(2,i,ix,iy)+b%a(4,i,ix,iy)*dy-b%a(5,i,ix,iy)*dx
    az(2)= b%a(1,i,ix,iy)+b%a(4,i,ix,iy)*dx+b%a(5,i,ix,iy)*dy
    af(1)=-b%a(3,i,ix,iy)*dy-b%a(6,i,ix,iy)*dx*dy-b%a(7,i,ix,iy)*dy**2/two
    ax(1)=-b%a(6,i,ix,iy)*dy
    ax(2)=-b%a(3,i,ix,iy)-b%a(6,i,ix,iy)*dx-b%a(7,i,ix,iy)*dy



    af(1)=af(1)-b%a(8,i,ix,iy)*dx**2*dy-b%a(9,i,ix,iy)*dy**3/3.d0    -b%a(13,i,ix,iy)*dx*dy**2/2.d0
    ax(1)=ax(1)-2.d0*b%a(8,i,ix,iy)*dx*dy                            -b%a(13,i,ix,iy)*dy**2/2.d0
    ax(2)=ax(2)-b%a(8,i,ix,iy)*dx**2-b%a(9,i,ix,iy)*dy**2   -b%a(13,i,ix,iy)*dx*dy

    az(1)=az(1)+b%a(10,i,ix,iy)*2.d0*dx*dy-b%a(12,i,ix,iy)*(dx**2-dy**2)/2.d0-b%a(5,i,ix,iy)*dx**2/2.0_dp/x0
    az(2)=az(2)+b%a(10,i,ix,iy)*dx**2+b%a(11,i,ix,iy)*dy**2+b%a(12,i,ix,iy)*dx*dy



    bf(1)=0.d0
    bf(2)=0.d0
    bf(3)=0.d0

    do k=1,3
       af(k)=-scale*af(k)
       az(k)=-scale*az(k)
       ax(k)=-scale*ax(k)
       bf(k)=-scale*bf(k)
    enddo


    call kill(dx,dy)

  end subroutine get_adata2c_p

  subroutine px_to_xpr0_hr(x,rho,hc,af,p)         ! convert in drift
    implicit none
    type(MAGNET_CHART), intent(in) :: p
    real(dp), intent(inout):: x(6)
    real(dp), intent(in):: af(3)
    real(dp), intent(in):: hc,rho
    real(dp) pz,d
    real(dp) beta0

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif


    d=root(one+two*x(5)/beta0+x(5)**2)
    pz=root(d**2-(x(2)-af(1))**2-(x(4)-af(2))**2)

    x(2)=(one+hc*(x(1)-rho))*(x(2)-af(1))/pz
    x(4)=(one+hc*(x(1)-rho))*(x(4)-af(2))/pz

  end subroutine px_to_xpr0_hr

  subroutine px_to_xpr0_hp(x,rho,hc,af,p)         ! convert in drift
    implicit none
    type(MAGNET_CHART), intent(in) :: p
    type(real_8), intent(inout):: x(6)
    type(real_8), intent(in):: af(3)
    real(dp), intent(in):: hc,rho
    type(real_8) pz,d
    real(dp) beta0

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif
    call alloc(pz,d)

    d=SQRT(one+two*x(5)/beta0+x(5)**2)
    pz=SQRT(d**2-(x(2)-af(1))**2-(x(4)-af(2))**2)

    x(2)=(one+hc*(x(1)-rho))*(x(2)-af(1))/pz
    x(4)=(one+hc*(x(1)-rho))*(x(4)-af(2))/pz
    call kill(pz,d)

  end subroutine px_to_xpr0_hp

  subroutine xp_to_pxr0_hr(x,rho,hc,af,p)         ! convert in drift
    implicit none

    real(dp), intent(inout):: x(6)
    real(dp), intent(in):: af(3)
    type(MAGNET_CHART), intent(in) :: p
    real(dp) pz,d
    real(dp) beta0
    real(dp), intent(in):: hc,rho

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif

    d=root(one+two*x(5)/beta0+x(5)**2)
    pz=root((one+hc*(x(1)-rho))**2+x(2)**2+x(4)**2)

    x(2)=d*x(2)/pz+af(1)
    x(4)=d*x(4)/pz+af(2)

  end subroutine xp_to_pxr0_hr

  subroutine xp_to_pxr0_hp(x,rho,hc,af,p)         ! convert in drift
    implicit none

    type(real_8), intent(inout):: x(6)
    type(real_8), intent(in):: af(3)
    type(MAGNET_CHART), intent(in) :: p
    type(real_8) pz,d
    real(dp) beta0
    real(dp), intent(in):: hc,rho

    if(p%time) then
       beta0=p%beta0 ;
    else
       beta0=one;
    endif
    call alloc(pz,d)
    d=SQRT(one+two*x(5)/beta0+x(5)**2)
    pz=SQRT((one+hc*(x(1)-rho))**2+x(2)**2+x(4)**2)

    x(2)=d*x(2)/pz+af(1)
    x(4)=d*x(4)/pz+af(2)
    call kill(pz,d)
  end subroutine xp_to_pxr0_hp

  logical(lp) function che(x,b)
    implicit none

    real(dp), intent(in):: x(6)
    type(FITTED_MAGNET), intent(in) :: b
    real(dp) :: s=c_1_2

    che=.false.

    if(x(1)<b%xmin) che=.true.
    if(x(1)>b%xmax) che=.true.
    if(x(3)<b%ymin) che=.true.
    if(x(3)>b%ymax) che=.true.

  end function che


  function dev_f(d,pos,ind,i,brho)
    use precision_constants
    implicit none
    real(dp) dev_f,brho,dh,ds
    integer  i(:),ii
    integer c,j,m(3),ind,pos(3)
    type(DATA_MAGNET) d


    dev_f=0.0_dp


    if(size(i)/=3) stop 888
    c=0
    do j=1,size(i)
       if(i(j)/=0) c=c+1
    enddo
    !write(6,*) " c = ",c
    select case(c)
    case(0)
       m=0
       dev_f=b_fun(m,d,pos,ind)
    case(1)
       !   m=0
       !   dev_f=-b_fun(m,d,pos,ind)
       !   m(i(1))=1
       !   write(6,*) i(1),dev_f*brho
       !   dev_f=dev_f+b_fun(m,d,pos,ind)
       !   write(6,*) m,dev_f*brho,b_fun(m,d,pos,ind)*brho
       m=0
       dev_f=-11.0_dp/6.0_dp*b_fun(m,d,pos,ind)
       m(i(1))=1
       dev_f=dev_f+3*b_fun(m,d,pos,ind)
       m(i(1))=2
       dev_f=dev_f-3*b_fun(m,d,pos,ind)/2.0_dp
       m(i(1))=3
       dev_f=dev_f+b_fun(m,d,pos,ind)/3.0_dp
       m=0
       m(i(1))=1+m(i(1))
    case(2)
       m=0
       m(i(1))=1+m(i(1))
       m(i(2))=1+m(i(2))
       dev_f=4*b_fun(m,d,pos,ind)
       m=0
       m(i(2))=1+m(i(2))
       dev_f=dev_f-3*b_fun(m,d,pos,ind)
       m=0
       m(i(1))=2+m(i(1))
       m(i(2))=1+m(i(2))
       dev_f=dev_f-b_fun(m,d,pos,ind)
       m=0
       m(i(1))=1+m(i(1))
       dev_f=dev_f-3*b_fun(m,d,pos,ind)
       m=0
       dev_f=dev_f+9.0_dp/4.0_dp*b_fun(m,d,pos,ind)
       m=0
       m(i(1))=2+m(i(1))
       dev_f=dev_f+3.0_dp*b_fun(m,d,pos,ind)/4.0_dp
       m=0
       m(i(1))=1+m(i(1))
       m(i(2))=2+m(i(2))
       dev_f=dev_f-b_fun(m,d,pos,ind)
       m=0
       m(i(2))=2+m(i(2))
       dev_f=dev_f+3.0_dp*b_fun(m,d,pos,ind)/4.0_dp
       m=0
       m(i(1))=2+m(i(1))
       m(i(2))=2+m(i(2))
       dev_f=dev_f+b_fun(m,d,pos,ind)/4.0_dp
       m=0
       m(i(1))=1+m(i(1))
       m(i(2))=1+m(i(2))
    case(3)
       m=0
       m(i(1))=1+m(i(1))
       m(i(2))=1+m(i(2))
       m(i(3))=1+m(i(3))
       dev_f=b_fun(m,d,pos,ind)

       m=0
       m(i(2))=1+m(i(2))
       m(i(3))=1+m(i(3))
       dev_f=dev_f-b_fun(m,d,pos,ind)
       m=0
       m(i(1))=1+m(i(1))
       m(i(3))=1+m(i(3))
       dev_f=dev_f-b_fun(m,d,pos,ind)
       m=0
       m(i(1))=1+m(i(1))
       m(i(2))=1+m(i(2))
       dev_f=dev_f-b_fun(m,d,pos,ind)

       m=0
       m(i(1))=1+m(i(1))
       dev_f=dev_f+b_fun(m,d,pos,ind)
       m=0
       m(i(2))=1+m(i(2))
       dev_f=dev_f+b_fun(m,d,pos,ind)
       m=0
       m(i(3))=1+m(i(3))
       dev_f=dev_f+b_fun(m,d,pos,ind)

       m=0
       dev_f=dev_f-b_fun(m,d,pos,ind)

       m=0
       m(i(1))=1+m(i(1))
       m(i(2))=1+m(i(2))
       m(i(3))=1+m(i(3))

    end select

    dev_f=dev_f*brho

    ds=d%x(1)*d%dtheta*DEG_TO_RAD_

    dh=1.0_dp

    !write(6,*) m

    do ii=1,m(1)
       dh=dh*d%dx
    enddo
    do ii=1,m(2)
       dh=dh*d%dy
    enddo
    do ii=1,m(3)
       dh=dh*ds
    enddo
    dev_f=dev_f/dh
    !write(6,*) "dh = ",dh
  end function dev_f


  function b_fun(x,d,pos,ind)
    use precision_constants
    implicit none
    type(DATA_MAGNET) d
    real(dp) b_fun
    integer  x(3),error,pos(3),ind,t(3),i

    error=0



    t=pos+x
    if(t(1)<0) then
       t(1)=1
       error=-1
    endif
    if(t(2)<0) then
       t(2)=1
       error=-2
    endif
    if(t(3)<0) then
       t(3)=1
       error=-3
    endif
    if(t(1)>d%nx) then
       t(1)=d%nx
       error=1
    endif
    if(t(2)>d%ny) then
       t(2)=d%ny
       error=2
    endif
    if(t(3)>d%ns) then
       t(3)=d%ns
       error=3
    endif
    if(error/=0) then
       er=er+1
    endif
    !write(6,*) "*********************************"
    !write(6,*) d%x(1),d%dtheta,d%dx,d%dy,ds,dh
    !write(6,*) "*********************************"

    b_fun=D%B(ind,0,t(3),t(1),t(2))  ! times brho temp


    !write(6,*) D%B(ind,t(3),t(1),t(2))



  end function b_fun


end module fitted_MAG
