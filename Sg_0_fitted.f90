!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis

module fitted_MAG
  USE S_status
  implicit none
  PRIVATE read_d,fit,get_bdata2c_r,get_bdata2c_p,get_adata2c_p   !  ,get_bdata
  PRIVATE POINTERS_user1R,POINTERS_user1P
  logical(lp) :: extrapolate =.true.
  INTEGER :: NX_0=0,NS_0=0,NY_0=0
  PRIVATE ALLOC_user1,KILL_user1,copy_el_elp,copy_elP_el,copy_el_el
  PRIVATE POINTERS_D,ZEROR_user1,ZEROP_user1,INTR,INTP,INTS
  PRIVATE fxr,fxp,fx,px_to_xpr0_r,px_to_xpp0_r,px_to_xp_r
  PRIVATE xp_to_px_r,xp_to_pxr0_r,xp_to_pxp0_r
  PRIVATE track_interpol_2,track_interpol_2_r,track_interpol_2_p
  PRIVATE track_interpol_4,track_interpol_4_r,track_interpol_4_P
  PRIVATE track_interpol_1h_r,track_interpol_1h_p
  ! Vector potential stuff
  PRIVATE get_adata2c_r,fxhr,fxhp,px_to_xpr0_hr,px_to_xpr0_hp,xp_to_pxr0_hr
  !frs real(dp),private::eps_fitted=c_1d_11
  private get_adatag,fxhpm,fxhrm, che
  logical(lp)::point_at=.true.,COPY_FIT=.TRUE.
  TYPE DATA_MAGNET
     INTEGER, POINTER :: NS,NX,NY,ORDER
     real(dp), POINTER :: dtheta,dx,dy
     real(dp), POINTER,DIMENSION(:)::x,y,theta
     REAL(sp),POINTER,DIMENSION(:,:,:,:)::B
     real(dp),POINTER,DIMENSION(:,:,:,:):: A
  END TYPE DATA_MAGNET


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE ZEROr_user1
     MODULE PROCEDURE ZEROp_user1
  END INTERFACE

  INTERFACE read
     MODULE PROCEDURE read_d
  END INTERFACE

  INTERFACE get_bdata
     MODULE PROCEDURE get_bdata2c_r
     MODULE PROCEDURE get_bdata2c_p
     MODULE PROCEDURE get_adata2c_r
     MODULE PROCEDURE get_adata2c_p
  END INTERFACE

  TYPE FITTED_MAGNET
     TYPE(MAGNET_CHART), POINTER:: P
     real(dp), POINTER ::L
     real(dp),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
     real(dp), POINTER ::SCALE
     logical(lp), POINTER ::symplectic
     TYPE(DATA_MAGNET), POINTER:: D
     real(dp),  DIMENSION(:,:), POINTER :: X
  END  TYPE FITTED_MAGNET

  TYPE FITTED_MAGNETP
     TYPE(MAGNET_CHART), POINTER:: P
     TYPE(REAL_8), POINTER ::L
     TYPE(REAL_8),  DIMENSION(:), POINTER :: AN,BN         !Multipole component
     logical(lp), POINTER ::symplectic
     !  SPECIAL POINTERS
     TYPE(REAL_8), POINTER ::SCALE
     TYPE(DATA_MAGNET), POINTER:: D
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
     MODULE PROCEDURE track_interpol_1h_r
     MODULE PROCEDURE track_interpol_1h_p
  END INTERFACE
  INTERFACE track_interpol_2
     MODULE PROCEDURE track_interpol_2_r
     MODULE PROCEDURE track_interpol_2_p
  END INTERFACE

  INTERFACE track_interpol_4
     MODULE PROCEDURE track_interpol_4_r
     MODULE PROCEDURE track_interpol_4_P
  END INTERFACE

  INTERFACE fx
     MODULE PROCEDURE fxr
     MODULE PROCEDURE fxp
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
    TYPE(DATA_MAGNET), INTENT(inout)::DD
    TYPE(DATA_MAGNET), INTENT(in)::D

    DD%NS=D%NS;DD%NX=D%NX;DD%NY=D%NY; DD%ORDER=D%ORDER;

    DD%dtheta=D%dtheta;DD%dx=D%dx;DD%dy=D%dy;

    DD%theta=D%theta;DD%x=D%x;DD%y=D%Y;

    IF(COPY_FIT) THEN
       IF(ASSOCIATED(D%B)) DD%B=D%B
       IF(ASSOCIATED(D%A)) DD%A=D%A
    ENDIF
  END SUBROUTINE copy_D_d


  SUBROUTINE POINTERS_D(D)
    IMPLICIT NONE
    TYPE(DATA_MAGNET), INTENT(inout)::D

    NULLIFY(D%NS,D%NX,D%NY,D%ORDER)
    ALLOCATE(D%NS,D%NX,D%NY,D%ORDER)
    D%NS=NS_0;D%NX=NX_0;D%NY=NY_0; D%ORDER=1;
    NULLIFY(D%dtheta,D%dx,D%dy)
    ALLOCATE(D%dtheta,D%dx,D%dy)
    D%dtheta=zero;D%dx=zero;D%dy=zero;
    NULLIFY(D%theta,D%x,D%y)
    ALLOCATE(D%theta(2),D%x(2),D%y(2))
    D%theta=zero;D%x=zero;D%y=zero;
    NULLIFY(D%B)
    IF(COPY_FIT) ALLOCATE(D%B(3,NS_0,NX_0,NY_0))
    NULLIFY(D%a)
    IF(COPY_FIT) ALLOCATE(D%a(9,NS_0,NX_0,NY_0))

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
    d%dx=(D%x(2)-D%x(1))/(D%nx-1)
    d%dy=(D%y(2)-D%y(1))/(D%ny-1)
    write(w_p%c(2),'(3(1x,g16.10))') d%dtheta,d%dx,d%dy
    call WRITE_I

    do i=1,d%ns
       read(mf,*) x1,x2,x3
       do j=1,d%nx
          do k=1,d%ny
             d%a(6,i,j,k)=zero
             read(mf,*) x1,x2,d%b(1,i,j,k),d%b(2,i,j,k),d%b(3,i,j,k)
             d%b(1,i,j,k)=d%b(1,i,j,k)/brho
             d%b(2,i,j,k)=d%b(2,i,j,k)/brho
             d%b(3,i,j,k)=d%b(3,i,j,k)/brho
          enddo
       enddo
    enddo
    close(mf)

    !posx=sss
    do i=1,d%ns
       do j=1,d%nx
          do k=1,d%ny
             call get_adatag(d,i,j,k)
             x1=-d%a(7,i,j,k)*d%dy-d%a(9,i,j,k)*d%dy**2/two
             if(k<d%ny) d%a(6,i,j,k+1)=d%a(6,i,j,k)+x1
          enddo
       enddo
    enddo

    !    CALL VECTOR_D(D)

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



  subroutine get_bdata2c_r(order,b,i,x,y,bf)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i,order
    real(dp), intent(in) :: x,y
    real(dp), intent(out) :: bf(3)
    integer j,k,l,ix,iy
    real(dp) posx,posy,db(3,2),x0,y0
    !  second order
    real(dp) db2(3,3),dt(0:4,0:4,3),d(0:4,0:4,3),f(-2:2,-2:2,3)
    real(dp) xm2,xm1,x1,x2,ym2,ym1,y1,y2



    posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    !  interpolate
    ix=int(posx)
    iy=int(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
    !write(6,*) "R: x0,y0",x0,y0

    if(order==1) then
       !  first order

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
             call write_e(2234)
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
             call write_e(2235)
          endif
       endif


       do j=1,3
          !db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix-1,iy))/b%dx/2
          !db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2
          db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix,iy))/b%dx
          db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy))/b%dy
       enddo

       do j=1,3
          bf(j)=b%b(j,i,ix,iy)+ db(j,1)*(x-x0)+db(j,2)*(y-y0)
       enddo
    elseif(order==2) then
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
             call write_e(2236)
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
             call write_e(2237)
          endif
       endif


       do j=1,3
          db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix-1,iy))/b%dx/2
          db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2
       enddo
       do j=1,3
          db2(j,1)=(b%b(j,i,ix+1,iy)+b%b(j,i,ix-1,iy)-2*b%b(j,i,ix,iy))/b%dx**2/2
          db2(j,2)=(b%b(j,i,ix,iy+1)+b%b(j,i,ix,iy-1)-2*b%b(j,i,ix,iy))/b%dy**2/2
       enddo
       do j=1,3
          db2(j,3)=(b%b(j,i,ix+1,iy+1)+b%b(j,i,ix-1,iy-1)-2*b%b(j,i,ix,iy))/2
          db2(j,3)=(db2(j,3)-db2(j,1)*b%dx**2-db2(j,2)*b%dy**2)/b%dx/b%dy
       enddo
       do j=1,3
          bf(j)=b%b(j,i,ix,iy)+ db(j,1)*(x-x0)+db(j,2)*(y-y0)
          bf(j)=bf(j)+db2(j,1)*(x-x0)**2+db2(j,2)*(y-y0)**2+db2(j,3)*(x-x0)*(y-y0)
       enddo


    elseIF(ORDER==4) THEN


       if(posx<3.or.posx>b%nx-2) then
          if(extrapolate) then
             ix=3
             if(posx>3) ix=b%nx-2
             x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A120))'
             write(w_p%c(1),'(a9,2(1x,i4))') " ix,b%nx ",ix,b%nx
             call write_e(2238)
          endif
       endif

       if(posy<3.or.posy>b%ny-2) then
          if(extrapolate) then
             iy=3
             if(posy>3) ix=b%ny-2
             y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A120))'
             write(w_p%c(1),'(a9,2(1x,i4))') " iy,b%ny ",iy,b%ny
             call write_e(2239)
          endif
       endif

       do k=0,4
          do l=0,4
             do j=1,3
                d(k,l,j)=zero
             enddo
          enddo
       enddo

       do j=1,3
          do k=-2,2
             do l=-2,2
                f(k,l,j)=b%b(j,i,ix+k,iy+l)
             enddo
          enddo
       enddo

       xm2=-2*b%dx
       xm1=-b%dx
       x1 =b%dx
       x2 =2*b%dx
       ym2=-2*b%dy
       ym1=-b%dy
       y1 =b%dy
       y2 =2*b%dy

       ! first order  !
       do j=1,3
          d(0,0,j)=f(0,0,j)
       enddo

       do j=1,3
          d(1,0,j)=(eight*(f(1,0,j)-f(-1,0,j))-(f(2,0,j)-f(-2,0,j)))/twelve/b%dx
          d(2,0,j)= (c_16*(f(1,0,j)+f(-1,0,j))-(f(2,0,j)+f(-2,0,j))-c_30*f(0,0,j))/c_24/b%dx**2
          d(3,0,j)= ((f(2,0,j)-f(-2,0,j))-two*(f(1,0,j)-f(-1,0,j)))/twelve/b%dx**3
          d(4,0,j)= ((f(2,0,j)+f(-2,0,j))-four*(f(1,0,j)+f(-1,0,j))+six*f(0,0,j))/c_24/b%dx**4
          d(0,1,j)= (eight*(f(0,1,j)-f(0,-1,j))-(f(0,2,j)-f(0,-2,j)))/twelve/b%dy
          d(0,2,j)= (c_16*(f(0,1,j)+f(0,-1,j))-(f(0,2,j)+f(0,-2,j))-c_30*f(0,0,j))/c_24/b%dy**2
          d(0,3,j)= ((f(0,2,j)-f(0,-2,j))-two*(f(0,1,j)-f(0,-1,j)))/twelve/b%dy**3
          d(0,4,j)= ((f(0,2,j)+f(0,-2,j))-four*(f(0,1,j)+f(0,-1,j))+six*f(0,0,j))/c_24/b%dy**4
       enddo

       !goto  1000

       do j=1,3
          dt(1,2,j)=( ( (f(1,2,j)+f(1,-2,j)) -   (f(2,1,j)+f(2,-1,j)))  + fit(d,x2,y1,j)+fit(d,x2,ym1,j)-(fit(d,x1,y2,j)+&
               fit(d,x1,ym2,j)) )/4/b%dx/b%dy**2
          dt(2,1,j)=( ( (f(2,1,j)+f(-2,1,j)) -   (f(1,2,j)+f(-1,2,j))) + fit(d,x1,y2,j)+fit(d,xm1,y2,j)-(fit(d,x2,y1,j)+&
               fit(d,xm2,y1,j)) )/4/b%dy/b%dx**2
          dt(2,2,j)=( ( (f(1,2,j)+f(1,-2,j)) - 2*(f(2,1,j)+f(2,-1,j))) + 2*(fit(d,x2,y1,j)+fit(d,x2,ym1,j))-(fit(d,x1,y2,j)+&
               fit(d,x1,ym2,j)) )
          dt(2,2,j)=dt(2,2,j)+( ( (f(2,1,j)+f(-2,1,j)) - 2*(f(1,2,j)+f(-1,2,j))) + 2*(fit(d,x1,y2,j)+fit(d,xm1,y2,j))-&
               (fit(d,x2,y1,j)+fit(d,xm2,y1,j)) )
          dt(2,2,j)=-dt(2,2,j)/eight/b%dx**2/b%dy**2/two
       enddo

       do j=1,3
          d(1,2,j)=dt(1,2,j)
          d(2,1,j)=dt(2,1,j)
          d(2,2,j)=dt(2,2,j)
       enddo

       do j=1,3
          dt(1,1,j)=(six*(f(1,2,j)+f(2,1,j))-five*(f(2,2,j)-four*f(1,1,j))-six*(fit(d,x2,y1,j)+fit(d,x1,y2,j))+&
               five*(fit(d,x2,y2,j) -four*fit(d,x1,y1,j)))/c_24/b%dx/b%dy
          d(1,1,j)=dt(1,1,j)
       enddo


       do j=1,3
          dt(1,3,j)=        -(four*f(1,-2,j)-f(2,-1,j)-four*fit(d,x1,ym2,j)+fit(d,x2,ym1,j)   )
          dt(1,3,j)=(dt(1,3,j)+(f(-2,1,j)-four*f(-1,2,j)-fit(d,xm2,y1,j)+four*fit(d,xm1,y2,j)   ))/60/b%dx/b%dy**3
          dt(3,1,j)=        -(-f(1,-2,j)+four*f(2,-1,j)+fit(d,x1,ym2,j)-four*fit(d,x2,ym1,j)   )
          dt(3,1,j)=(dt(3,1,j)+(-four*f(-2,1,j)+f(-1,2,j)+four*fit(d,xm2,y1,j)-fit(d,xm1,y2,j)    ))/60/b%dx**3/b%dy
          d(1,3,j)=dt(1,3,j)
          d(3,1,j)=dt(3,1,j)
       enddo

1000   continue

       do j=1,3
          bf(j)=zero
          do l=0,4
             do k=0,4
                bf(j)=bf(j)+d(l,k,j)*(x-x0)**l*(y-y0)**k
             enddo
          enddo
       enddo


    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A120))'
       write(w_p%c(1),'(a30,(1x,i4))') "FORBIDDEN INTERPOLATION ORDER ",order
       call write_e(2240)
    endif



  end subroutine get_bdata2c_r

  subroutine get_bdata2c_p(order,b,i,xx,yy,bf)
    implicit none
    type(DATA_MAGNET) b
    type(real_8), intent(in):: xx,yy
    type(real_8), intent(out):: bf(3)
    integer, intent(in) :: i,order

    integer j,k ,ix,iy,l
    real(dp) posx,posy,db(3,2),x,y,x0,y0
    !  second order
    real(dp) db2(3,3),dt(0:4,0:4,3),d(0:4,0:4,3),f(-2:2,-2:2,3)
    real(dp) xm2,xm1,x1,x2,ym2,ym1,y1,y2




    x=xx
    y=yy

    posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    !  interpolate
    ix=int(posx)
    iy=int(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)

    if(order==1) then
       !  first order


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
             call write_e(2234)
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
             call write_e(2235)
          endif
       endif

       do j=1,3
          !         db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix-1,iy))/b%dx/2
          !         db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2
          db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix,iy))/b%dx
          db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy))/b%dy
       enddo

       do j=1,3
          bf(j)=b%b(j,i,ix,iy)+ db(j,1)*(xx-x0)+db(j,2)*(yy-y0)
       enddo
    elseif(order==2) then
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
             call write_e(2236)
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
             call write_e(2237)
          endif
       endif



       do j=1,3
          db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix-1,iy))/b%dx/2
          db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2
       enddo
       do j=1,3
          db2(j,1)=(b%b(j,i,ix+1,iy)+b%b(j,i,ix-1,iy)-2*b%b(j,i,ix,iy))/b%dx**2/2
          db2(j,2)=(b%b(j,i,ix,iy+1)+b%b(j,i,ix,iy-1)-2*b%b(j,i,ix,iy))/b%dy**2/2
       enddo
       do j=1,3
          db2(j,3)=(b%b(j,i,ix+1,iy+1)+b%b(j,i,ix-1,iy-1)-2*b%b(j,i,ix,iy))/2
          db2(j,3)=(db2(j,3)-db2(j,1)*b%dx**2-db2(j,2)*b%dy**2)/b%dx/b%dy
       enddo
       do j=1,3
          bf(j)=b%b(j,i,ix,iy)+ db(j,1)*(xx-x0)+db(j,2)*(yy-y0)
          bf(j)=bf(j)+db2(j,1)*(xx-x0)**2+db2(j,2)*(yy-y0)**2+db2(j,3)*(xx-x0)*(yy-y0)
       enddo



    elseif(order==4) then


       if(posx<3.or.posx>b%nx-2) then
          if(extrapolate) then
             ix=3
             if(posx>3) ix=b%nx-2
             x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A120))'
             write(w_p%c(1),'(a9,2(1x,i4))') " ix,b%nx ",ix,b%nx
             call write_e(2237)
          endif
       endif

       if(posy<3.or.posy>b%ny-2) then
          if(extrapolate) then
             iy=3
             if(posy>3) ix=b%ny-2
             y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
          else
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A120))'
             write(w_p%c(1),'(a9,2(1x,i4))') " iy,b%ny ",iy,b%ny
             call write_e(2237)
          endif
       endif

       do k=0,4
          do l=0,4
             do j=1,3
                d(k,l,j)=zero
             enddo
          enddo
       enddo

       do j=1,3
          do k=-2,2
             do l=-2,2
                f(k,l,j)=b%b(j,i,ix+k,iy+l)
             enddo
          enddo
       enddo

       xm2=-2*b%dx
       xm1=-b%dx
       x1 =b%dx
       x2 =2*b%dx
       ym2=-2*b%dy
       ym1=-b%dy
       y1 =b%dy
       y2 =2*b%dy

       ! first order  !
       do j=1,3
          d(0,0,j)=f(0,0,j)
       enddo

       do j=1,3
          d(1,0,j)=(eight*(f(1,0,j)-f(-1,0,j))-(f(2,0,j)-f(-2,0,j)))/twelve/b%dx
          d(2,0,j)= (c_16*(f(1,0,j)+f(-1,0,j))-(f(2,0,j)+f(-2,0,j))-c_30*f(0,0,j))/c_24/b%dx**2
          d(3,0,j)= ((f(2,0,j)-f(-2,0,j))-two*(f(1,0,j)-f(-1,0,j)))/twelve/b%dx**3
          d(4,0,j)= ((f(2,0,j)+f(-2,0,j))-four*(f(1,0,j)+f(-1,0,j))+six*f(0,0,j))/c_24/b%dx**4
          d(0,1,j)= (eight*(f(0,1,j)-f(0,-1,j))-(f(0,2,j)-f(0,-2,j)))/twelve/b%dy
          d(0,2,j)= (c_16*(f(0,1,j)+f(0,-1,j))-(f(0,2,j)+f(0,-2,j))-c_30*f(0,0,j))/c_24/b%dy**2
          d(0,3,j)= ((f(0,2,j)-f(0,-2,j))-two*(f(0,1,j)-f(0,-1,j)))/twelve/b%dy**3
          d(0,4,j)= ((f(0,2,j)+f(0,-2,j))-four*(f(0,1,j)+f(0,-1,j))+six*f(0,0,j))/c_24/b%dy**4
       enddo


       do j=1,3
          dt(1,2,j)=( ( (f(1,2,j)+f(1,-2,j)) -   (f(2,1,j)+f(2,-1,j)))  + fit(d,x2,y1,j)+fit(d,x2,ym1,j)-(fit(d,x1,y2,j)+&
               fit(d,x1,ym2,j)) )/4/b%dx/b%dy**2
          dt(2,1,j)=( ( (f(2,1,j)+f(-2,1,j)) -   (f(1,2,j)+f(-1,2,j))) + fit(d,x1,y2,j)+fit(d,xm1,y2,j)-(fit(d,x2,y1,j)+&
               fit(d,xm2,y1,j)) )/4/b%dy/b%dx**2
          dt(2,2,j)=( ( (f(1,2,j)+f(1,-2,j)) - 2*(f(2,1,j)+f(2,-1,j))) + 2*(fit(d,x2,y1,j)+fit(d,x2,ym1,j))-(fit(d,x1,y2,j)+&
               fit(d,x1,ym2,j)) )
          dt(2,2,j)=dt(2,2,j)+( ( (f(2,1,j)+f(-2,1,j)) - 2*(f(1,2,j)+f(-1,2,j))) + 2*(fit(d,x1,y2,j)+fit(d,xm1,y2,j))-&
               (fit(d,x2,y1,j)+fit(d,xm2,y1,j)) )
          dt(2,2,j)=-dt(2,2,j)/eight/b%dx**2/b%dy**2/two
       enddo

       do j=1,3
          d(1,2,j)=dt(1,2,j)
          d(2,1,j)=dt(2,1,j)
          d(2,2,j)=dt(2,2,j)
       enddo

       do j=1,3
          dt(1,1,j)=(six*(f(1,2,j)+f(2,1,j))-five*(f(2,2,j)-four*f(1,1,j))-six*(fit(d,x2,y1,j)+fit(d,x1,y2,j))+&
               five*(fit(d,x2,y2,j) -four*fit(d,x1,y1,j)))/c_24/b%dx/b%dy
          d(1,1,j)=dt(1,1,j)
          !write(16,*) j, d(1,1,j),d(2,2,j),d(2,0,j)
       enddo
       !goto 1000


       do j=1,3
          dt(1,3,j)=        -(four*f(1,-2,j)-f(2,-1,j)-four*fit(d,x1,ym2,j)+fit(d,x2,ym1,j)   )
          dt(1,3,j)=(dt(1,3,j)+(f(-2,1,j)-four*f(-1,2,j)-fit(d,xm2,y1,j)+four*fit(d,xm1,y2,j)   ))/60/b%dx/b%dy**3
          dt(3,1,j)=        -(-f(1,-2,j)+four*f(2,-1,j)+fit(d,x1,ym2,j)-four*fit(d,x2,ym1,j)   )
          dt(3,1,j)=(dt(3,1,j)+(-four*f(-2,1,j)+f(-1,2,j)+four*fit(d,xm2,y1,j)-fit(d,xm1,y2,j)    ))/60/b%dx**3/b%dy
          d(1,3,j)=dt(1,3,j)
          d(3,1,j)=dt(3,1,j)
       enddo

1000   continue

       do j=1,3
          bf(j)=zero
          do l=0,4
             do k=0,4
                bf(j)=bf(j)+d(l,k,j)*(xx-x0)**l*(yy-y0)**k
             enddo
          enddo
       enddo

    ELSE

       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A120))'
       write(w_p%c(1),'(a30,(1x,i4))') "FORBIDDEN INTERPOLATION ORDER ",order
       call write_e(2240)
    endif


  end subroutine get_bdata2c_p


  real(dp) function fit(d,x,y,k)
    implicit none
    real(dp) x,y,d(0:4,0:4,3)
    integer i,j,k


    fit=zero
    do i=0,4
       do j=0,4
          fit=fit+d(i,j,k)*x**i*y**j
       enddo
    enddo

  end function  fit


  !!   OPERATIONS ON TYPE FITTED_MAGNET

  SUBROUTINE ZEROR_user1(EL,I)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%SCALE))  THEN
          DEALLOCATE(EL%SCALE)
          DEALLOCATE(EL%symplectic)
          DEALLOCATE(EL%D%NS,EL%D%NX,EL%D%NY)
          DEALLOCATE(EL%D%dtheta,EL%D%dx,EL%D%dY)
          DEALLOCATE(EL%D%theta,EL%D%x ,EL%D%Y)
          DEALLOCATE(EL%X)
          DEALLOCATE(EL%d%order)
          IF(ASSOCIATED(EL%D%B)) DEALLOCATE(EL%D%B)
          IF(ASSOCIATED(EL%D%A)) DEALLOCATE(EL%D%A)
          DEALLOCATE(EL%D)
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%symplectic)
       NULLIFY(EL%SCALE)
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
          DEALLOCATE(EL%SCALE)
          DEALLOCATE(EL%symplectic)
          DEALLOCATE(EL%D%NS,EL%D%NX,EL%D%NY)
          DEALLOCATE(EL%D%dtheta,EL%D%dx,EL%D%dY)
          DEALLOCATE(EL%D%theta,EL%D%x ,EL%D%Y)
          DEALLOCATE(EL%D%order)
          IF(ASSOCIATED(EL%D%A)) DEALLOCATE(EL%D%A)
          IF(ASSOCIATED(EL%D%B)) DEALLOCATE(EL%D%B)
          DEALLOCATE(EL%D)
       ENDIF
    elseif(i==0)       then
       ! nullifies pointers
       NULLIFY(EL%symplectic)
       NULLIFY(EL%SCALE)
       NULLIFY(EL%D)
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZEROP_user1

  SUBROUTINE ALLOC_user1(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL
    CALL ALLOC(EL%SCALE)
    ! ALLOC INTERNAL POLYMORPHS IF ANY
  END SUBROUTINE ALLOC_user1

  SUBROUTINE KILL_user1(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL

    CALL KILL(EL%SCALE)
    ! KILL INTERNAL POLYMORPHS IF ANY

  END SUBROUTINE KILL_user1

  SUBROUTINE POINTERS_user1R(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(INOUT)::EL

    ALLOCATE(EL%SCALE)
    ALLOCATE(EL%symplectic)
    ALLOCATE(EL%D)
    ALLOCATE(EL%X(2,NS_0))
    CALL POINTERS_D(EL%D)


    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_user1R

  SUBROUTINE POINTERS_user1P(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL

    ALLOCATE(EL%SCALE)
    ALLOCATE(EL%symplectic)
    ALLOCATE(EL%D)
    CALL POINTERS_D(EL%D)
    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_user1P

  SUBROUTINE reset_FITTED_MAGNETP(EL)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(INOUT)::EL


    CALL resetpoly_R31(EL%SCALE)
    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS

  END SUBROUTINE reset_FITTED_MAGNETP


  SUBROUTINE copy_el_elp(EL,ELP)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(in)::EL
    TYPE(FITTED_MAGNETP), INTENT(inout)::ELP

    ELP%SCALE    =EL%SCALE
    ELP%SYMPLECTIC    =EL%SYMPLECTIC
    CALL COPY(EL%D,ELP%D)

    !  COPY CODING HERE NO ALLOCATION OF POINTERS OR POLYMORPH NEEDED
    !  IF DONE CORRECTLY

  END SUBROUTINE copy_el_elp

  SUBROUTINE copy_elp_el(EL,ELP)
    IMPLICIT NONE
    TYPE(FITTED_MAGNETP), INTENT(in)::EL
    TYPE(FITTED_MAGNET), INTENT(inout)::ELP




    ELP%SCALE    =EL%SCALE
    ELP%SYMPLECTIC    =EL%SYMPLECTIC
    CALL COPY(EL%D,ELP%D)
    !  COPY CODING HERE NO ALLOCATION OF POINTERS OR POLYMORPH NEEDED
    !  IF DONE CORRECTLY


  END SUBROUTINE copy_elp_el

  SUBROUTINE copy_el_el(EL,ELP)
    IMPLICIT NONE
    TYPE(FITTED_MAGNET), INTENT(in)::EL
    TYPE(FITTED_MAGNET), INTENT(inout)::ELP

    ELP%SCALE    =EL%SCALE
    ELP%SYMPLECTIC    =EL%SYMPLECTIC
    CALL COPY(EL%D,ELP%D)
    !  COPY CODING HERE NO ALLOCATION OF POINTERS


  END SUBROUTINE copy_el_el

  !  TRACKING STUFF

  SUBROUTINE INTR(EL,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(FITTED_MAGNET),INTENT(INOUT):: EL

    x(1)=x(1)+one/el%p%b0-el%d%x(1)


    select case(el%p%method)
    case(1)
       if(el%d%order/=1 ) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only order 1 possible "
          call write_e(500)
       endif
       call track_interpol_1(el,x)
    case(2)
       if(el%symplectic)then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only Method 1  possible "
          call write_e(600)
       endif
       if(el%d%order/=1.and.el%d%order/=2.and.el%d%order/=4 ) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only order 1,2,4 possible "
          call write_e(501)
       endif
       call track_interpol_2(el,x)
    case(4)
       if(el%symplectic)then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only Method 1  possible "
          call write_e(601)
       endif
       if(el%d%order/=1.and.el%d%order/=2.and.el%d%order/=4 ) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only order 1,2,4 possible "
          call write_e(502)
       endif
       call track_interpol_4(el,x)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A120))'
       w_p%c(1)= "METHOD NOT SUPPORTED IN track_interpol"
       call write_e(101)
    END SELECT

    x(1)=x(1)-one/el%p%b0+el%d%x(1)

  END SUBROUTINE INTR






  SUBROUTINE INTP(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(FITTED_MAGNETP),INTENT(INOUT):: EL

    x(1)=x(1)+one/el%p%b0-el%d%x(1)

    select case(el%p%method)
    case(1)
       if(el%d%order/=1 ) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only order 1 possible "
          call write_e(500)
       endif
       call track_interpol_1(el,x)
    case(2)
       if(el%symplectic)then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only Method 1  possible "
          call write_e(600)
       endif
       if(el%d%order/=1.and.el%d%order/=2.and.el%d%order/=4 ) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only order 1,2,4 possible "
          call write_e(501)
       endif
       call track_interpol_2(el,x)
    case(4)
       if(el%symplectic)then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only Method 1  possible "
          call write_e(601)
       endif
       if(el%d%order/=1.and.el%d%order/=2.and.el%d%order/=4 ) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A120))'
          w_p%c(1)= " Only order 1,2,4 possible "
          call write_e(502)
       endif
       call track_interpol_4(el,x)
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A120))'
       w_p%c(1)= "METHOD NOT SUPPORTED IN track_interpol"
       call write_e(101)
    END SELECT


    x(1)=x(1)-one/el%p%b0+el%d%x(1)

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



  subroutine fxr(f,x,b,r,hcurv,p)
    implicit none

    real(dp)  d(3),c(6),BETA0,GAMMA0I
    real(dp) ,intent(in) :: b(3),r ,hcurv
    type(MAGNET_CHART), intent(in) :: p
    real(dp) ,intent(inout) :: x(6)
    real(dp), intent(out):: f(6)

    if(p%time) then
       beta0=p%beta0;GAMMA0I=p%GAMMA0I;
    else
       beta0=one;GAMMA0I=zero;
    endif
    x(1)=x(1)-r
    d(1)=root(x(2)**2+x(4)**2+(one+hcurv*x(1))**2)
    d(2)=(d(1)**3)/root(one+2*x(5)/beta0+x(5)**2)
    d(3)=one+hcurv*x(1)

    c(1)=d(1)**2-x(2)**2
    c(2)=-x(2)*x(4)
    c(3)= x(2)*x(4)
    c(4)=-d(1)**2+x(4)**2
    c(5)=d(2)*(x(4)*b(3)-d(3)*b(2)) +hcurv*d(3)*(d(1)**2+x(2)**2)
    c(6)=d(2)*(x(2)*b(3)-d(3)*b(1)) -hcurv*d(3)*c(3)

    d(3)=c(1)*c(4)-c(2)*c(3)
    f(1)=x(2)
    f(2)=(c(4)*c(5)-c(2)*c(6))/d(3)
    f(3)=x(4)
    f(4)=(c(1)*c(6)-c(3)*c(5))/d(3)
    f(5)=0
    d(2)=one+two*x(5)/beta0+x(5)**2
    !   d(2)=SQRT((one+d(2)*gambet)/d(2)/gambet)
    !   f(6)=d(2)*d(1)

    d(2)=gamma0I/beta0/d(2)
    f(6)=root((1+d(2)**2))*d(1)  ! (time)-prime = dt/dz
    x(1)=x(1)+r
  end subroutine fxr

  subroutine fxp(f,x,b,r,hcurv,p)
    implicit none

    type(real_8)  d(3),c(6)
    type(real_8) ,intent(inout) :: x(6)
    type(real_8) ,intent(in) :: b(3)
    real(dp) ,intent(in) ::  hcurv,r
    real(dp)   BETA0,GAMMA0I
    type(real_8), intent(out):: f(6)
    type(MAGNET_CHART), intent(in) :: p

    call alloc(d,3)
    call alloc(c,6)

    if(p%time) then
       beta0=p%beta0;GAMMA0I=p%GAMMA0I;
    else
       beta0=one;GAMMA0I=zero;
    endif

    x(1)=x(1)-r
    d(1)=SQRT(x(2)**2+x(4)**2+(one+hcurv*x(1))**2)
    d(2)=(d(1)**3)/SQRT(one+2*x(5)/beta0+x(5)**2)
    d(3)=one+hcurv*x(1)

    c(1)=d(1)**2-x(2)**2
    c(2)=-x(2)*x(4)
    c(3)= x(2)*x(4)
    c(4)=-d(1)**2+x(4)**2
    c(5)=d(2)*(x(4)*b(3)-d(3)*b(2)) +hcurv*d(3)*(d(1)**2+x(2)**2)
    c(6)=d(2)*(x(2)*b(3)-d(3)*b(1)) -hcurv*d(3)*c(3)

    d(3)=c(1)*c(4)-c(2)*c(3)
    f(1)=x(2)
    f(2)=(c(4)*c(5)-c(2)*c(6))/d(3)
    f(3)=x(4)
    f(4)=(c(1)*c(6)-c(3)*c(5))/d(3)
    f(5)=0
    d(2)=one+two*x(5)/beta0+x(5)**2
    !   d(2)=SQRT((one+d(2)*gambet)/d(2)/gambet)
    !   f(6)=d(2)*d(1)

    d(2)=gamma0I/beta0/d(2)
    f(6)=SQRT((1+d(2)**2))*d(1)  ! (time)-prime = dt/dz
    x(1)=x(1)+r
    call kill(d,3)
    call kill(c,6)
  end subroutine fxp


  subroutine track_interpol_2_r(bend,x)
    implicit none
    type(fitted_magnet), intent(inOUT):: bend
    real(dp) , intent(inout)::x(6)
    real(dp) b(3),b1(3),f(6),k1(6),x1(6),hc,rho,dz,DAL,rhob0
    integer  i,j,k
    integer knf,nf

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)
    rhob0=one/bend%p%b0
    !  rotating
    nf=bend%p%nst

    !  going into non-canonical coordinates
    call px_to_xp_r(x,hc,bend%p)

    x(1)=x(1)+rho


    dz=(bend%d%dtheta*DEG_TO_RAD_)*rho/nf
    DAL=bend%d%dtheta*DEG_TO_RAD_

    !tot=zero
    do i=1,bend%d%ns-1
       BEND%X(1,I)=(X(1))*COS((I-1)*DAL)-rhob0
       BEND%X(2,I)=(X(1))*SIN((I-1)*DAL)
       do knf=1,nf
          call get_bdata(bend%D%order,bend%d,i,x(1),x(3),b)
          ! original formula for electrons (q<0)  now switch to proton

          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo

          call get_bdata(bend%D%order,bend%d,i+1,x(1),x(3),b1)
          ! original formula for electrons (q<0)  now switch to proton
          do k=1,3
             b1(k)=-bend%scale*b1(k)
          enddo

          do k=1,3
             b(k)=b(k)+(b1(k)-b(k))*REAL(knf-1,kind=DP)/nf
          enddo
          call fx(k1,x,b,rho,hc,bend%p)

          do j=1,6
             x1(j)=x(j)+dz*k1(j)         ! temporary
          enddo

          call get_bdata(bend%d%order,bend%d,i,x1(1),x1(3),b)

          ! original formula for electrons (q<0)  now switch to proton
          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo

          call get_bdata(bend%d%order,bend%d,i+1,x1(1),x1(3),b1)

          ! original formula for electrons (q<0)  now switch to proton
          do k=1,3
             b1(k)=-bend%scale*b1(k)
          enddo
          do k=1,3
             b(k)=b(k)+(b1(k)-b(k))*REAL(knf,kind=DP)/nf
          enddo
          !
          call fx(f,x1,b,rho,hc,bend%p)


          do j=1,6
             x(j)=x(j)+dz*(f(j)+k1(j))/two
          enddo
          !        tot=tot+bend%dz
       enddo !nf


    enddo
    BEND%X(1,bend%d%ns)=(X(1))*COS((bend%d%ns-1)*DAL)-rhob0
    BEND%X(2,bend%d%ns)=(X(1))*SIN((bend%d%ns-1)*DAL)

    !write(6,*) tot
    x(1)=x(1)-rho
    call xp_to_px_r(x,hc,bend%p)



  end subroutine track_interpol_2_r

  subroutine track_interpol_2_p(bend,x)
    implicit none
    type(fitted_magnetp), intent(inOUT):: bend
    type(real_8), intent(inout)::x(6)
    real(dp) hc,rho,dz,DAL
    integer  i,j,k
    integer knf,nf
    type(real_8) b(3),b1(3),f(6),k1(6),x1(6)

    call alloc(b,3);call alloc(b1,3);call alloc(f,6);
    call alloc(k1,6);call alloc(x1,6);

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)


    !  rotating
    nf=bend%p%nst

    !  going into non-canonical coordinates
    call px_to_xp_r(x,hc,bend%p)

    x(1)=x(1)+rho


    dz=(bend%d%dtheta*DEG_TO_RAD_)*rho/nf
    DAL=bend%d%dtheta*DEG_TO_RAD_



    !tot=zero
    do i=1,bend%d%ns-1
       do knf=1,nf
          call get_bdata(bend%D%order,bend%d,i,x(1),x(3),b)
          ! original formula for electrons (q<0)  now switch to proton

          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo

          call get_bdata(bend%D%order,bend%d,i+1,x(1),x(3),b1)
          ! original formula for electrons (q<0)  now switch to proton
          do k=1,3
             b1(k)=-bend%scale*b1(k)
          enddo

          do k=1,3
             b(k)=b(k)+(b1(k)-b(k))*REAL(knf-1,kind=DP)/nf
          enddo
          call fx(k1,x,b,rho,hc,bend%p)

          do j=1,6
             x1(j)=x(j)+dz*k1(j)         ! temporary
          enddo

          call get_bdata(bend%d%order,bend%d,i,x1(1),x1(3),b)

          ! original formula for electrons (q<0)  now switch to proton
          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo

          call get_bdata(bend%d%order,bend%d,i+1,x1(1),x1(3),b1)

          ! original formula for electrons (q<0)  now switch to proton
          do k=1,3
             b1(k)=-bend%scale*b1(k)
          enddo
          do k=1,3
             b(k)=b(k)+(b1(k)-b(k))*REAL(knf,kind=DP)/nf
          enddo
          !
          call fx(f,x1,b,rho,hc,bend%p)


          do j=1,6
             x(j)=x(j)+dz*(f(j)+k1(j))/two
          enddo
          !        tot=tot+bend%dz
       enddo !nf


    enddo

    x(1)=x(1)-rho
    call xp_to_px_r(x,hc,bend%p)


    call kill(b,3);call kill(b1,3);call kill(f,6);
    call kill(k1,6);call kill(x1,6);

  end subroutine track_interpol_2_p


  subroutine track_interpol_4_r(bend,x)
    implicit none
    type(fitted_magnet), intent(inout):: bend
    real(dp) , intent(inout)::x(6)
    real(dp)        b(3),f(6),K1(6),K2(6),K3(6),K4(6),X1(6),X2(6),X3(6)
    real(dp) hc,rho,dz,dal,rhob0
    integer  i,j,k
    integer knf,nf

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)
    rhob0=one/bend%p%b0


    !  rotating
    nf=bend%p%nst

    !  going into non-canonical coordinates
    call px_to_xp_r(x,hc,bend%p)

    x(1)=x(1)+rho


    dz=(bend%d%dtheta*DEG_TO_RAD_)*rho/nf
    DAL=bend%d%dtheta*DEG_TO_RAD_

    do i=1,bend%d%ns-2,2
       BEND%X(1,I)=(X(1))*COS((I-1)*DAL)-rhob0
       BEND%X(2,I)=(X(1))*SIN((I-1)*DAL)
       do knf=1,nf
          call get_bdata(bend%d%order,bend%d,i,x(1),x(3),b)
          ! original formula for electrons (q<0)  now switch to proton

          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          !
          call fx(f,x,b,rho,hc,bend%p)
          do j=1,6
             k1(j)=two*dz*f(j)                ! temporary
             x1(j)=x(j)+dz*f(j)
          enddo
          call get_bdata(bend%d%order,bend%d,mod_s(I+1,bend%d%ns),x1(1),x1(3),b)
          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          call fx(f,x1,b,rho,hc,bend%p)

          do j=1,6
             k2(j)=two*dz*f(j)
             x2(j)=x(j)+dz*f(j)
          enddo

          call get_bdata(bend%d%order,bend%d,mod_s(I+1,bend%d%ns),x2(1),x2(3),b)
          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          call fx(f,x2,b,rho,hc,bend%p)

          do j=1,6
             k3(j)=two*dz*f(j)
             x3(j)=x(j)+k3(j)
          enddo

          call get_bdata(bend%d%order,bend%d,mod_s(I+2,bend%d%ns),x3(1),x3(3),b)

          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          call fx(f,x3,b,rho,hc,bend%p)

          do j=1,6
             k4(j)=two*dz*f(j)
             x(j)=x(j)+(k1(j)+two*k2(j)+two*k3(j)+K4(J))/six
          enddo
       enddo !nf
    enddo
    BEND%X(1,bend%d%ns)=(X(1))*COS((bend%d%ns-1)*DAL)-rhob0
    BEND%X(2,bend%d%ns)=(X(1))*SIN((bend%d%ns-1)*DAL)

    x(1)=x(1)-rho
    call xp_to_px_r(x,hc,bend%p)

    do i=1,bend%d%ns-2,2
       BEND%X(1,I+1)=half*(BEND%X(1,I)+BEND%X(1,I+2))
       BEND%X(2,I+1)=half*(BEND%X(2,I)+BEND%X(2,I+2))
    ENDDO


  end subroutine track_interpol_4_r

  subroutine track_interpol_4_P(bend,x)
    implicit none
    type(fitted_magnetP), intent(inout):: bend
    TYPE(REAL_8)  , intent(inout)::x(6)
    TYPE(REAL_8)  b(3),f(6),K1(6),K2(6),K3(6),K4(6),X1(6),X2(6),X3(6)
    real(dp) hc,rho,dz,dal
    integer  i,j,k
    integer knf,nf
    call alloc(b,3);call alloc(f,6);
    call alloc(k1,6);call alloc(K2,6);
    call alloc(k3,6);call alloc(K4,6);
    call alloc(X1,6);call alloc(X2,6);call alloc(X3,6);

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)


    !  rotating
    nf=bend%p%nst

    !  going into non-canonical coordinates
    call px_to_xp_r(x,hc,bend%p)

    x(1)=x(1)+rho


    dz=(bend%d%dtheta*DEG_TO_RAD_)*rho/nf
    DAL=bend%d%dtheta*DEG_TO_RAD_

    do i=1,bend%d%ns-2,2
       do knf=1,nf
          call get_bdata(bend%d%order,bend%d,i,x(1),x(3),b)
          ! original formula for electrons (q<0)  now switch to proton

          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          !
          call fx(f,x,b,rho,hc,bend%p)
          do j=1,6
             k1(j)=two*dz*f(j)                ! temporary
             x1(j)=x(j)+dz*f(j)
          enddo
          call get_bdata(bend%d%order,bend%d,mod_s(I+1,bend%d%ns),x1(1),x1(3),b)
          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          call fx(f,x1,b,rho,hc,bend%p)

          do j=1,6
             k2(j)=two*dz*f(j)
             x2(j)=x(j)+dz*f(j)
          enddo

          call get_bdata(bend%d%order,bend%d,mod_s(I+1,bend%d%ns),x2(1),x2(3),b)
          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          call fx(f,x2,b,rho,hc,bend%p)

          do j=1,6
             k3(j)=two*dz*f(j)
             x3(j)=x(j)+k3(j)
          enddo

          call get_bdata(bend%d%order,bend%d,mod_s(I+2,bend%d%ns),x3(1),x3(3),b)

          do k=1,3
             b(k)=-bend%scale*b(k)
          enddo
          call fx(f,x3,b,rho,hc,bend%p)

          do j=1,6
             k4(j)=two*dz*f(j)
             x(j)=x(j)+(k1(j)+two*k2(j)+two*k3(j)+K4(J))/six
          enddo
       enddo !nf
    enddo

    x(1)=x(1)-rho
    call xp_to_px_r(x,hc,bend%p)

    call KILL(b,3);call KILL(f,6);
    call KILL(k1,6);call KILL(K2,6);
    call KILL(k3,6);call KILL(K4,6);
    call KILL(X1,6);call KILL(X2,6);call KILL(X3,6);

  end subroutine track_interpol_4_P


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
    f(6)=(one/beta0)*d(1)/d(2)


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
    f(6)=(one/beta0)*d(1)/d(2)

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
    f(6)=(one/beta0)*d(1)/d(2)


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
    f(6)=(one/beta0)*d(1)/d(2)

    m(1,1)=d(1)*(one+(x(2)-a(1))/d(2)**2)/d(2)
    m(2,1)=m(1,1)*ax(2)
    m(1,1)=m(1,1)*ax(1)-hcurv*(x(2)-a(1))/d(2)
    m(2,2)=d(1)*(x(2)-a(1))*x(4)/d(2)**3
    m(1,2)=-hcurv*ax(2)/d(2)+m(2,2)*ax(2)
    m(2,2)=m(2,2)*ax(2)

    x(1)=x(1)+r

    call kill(d,2)
  end subroutine fxhpm

  subroutine track_interpol_1h_r(bend,x)
    implicit none
    type(fitted_magnet), intent(inOUT):: bend
    real(dp) , intent(inout)::x(6)
    real(dp) k1(6),xt(6),hc,rho,dz,DAL
    real(dp) af(3),az(3),ay(3),m(2,2),mi(2,2),dif(2),t(2),det
    integer  i,j,k,NF,KNF,it
    real(dp) norm, norm0,rhob0
    logical(lp) doit

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)
    rhob0=one/bend%p%b0
    !  rotating
    nf=bend%p%nst
    x(1)=x(1)+rho

    !  going into non-canonical coordinates
    call get_bdata(bend%d,1,x(1),x(3),af,az,ay)
    do k=1,3
       af(k)=-bend%scale*af(k)
       az(k)=-bend%scale*az(k)
       ay(k)=-bend%scale*ay(k)
    enddo
    call px_to_xp_r(x,rho,hc,af,bend%p)



    dz=(bend%d%dtheta*DEG_TO_RAD_)*rho/nf
    DAL=bend%d%dtheta*DEG_TO_RAD_

    !tot=zero
    do i=1,bend%d%ns-1
       DO KNF=1,NF
          BEND%X(1,I)=(X(1))*COS((I-1)*DAL)-rhob0
          BEND%X(2,I)=(X(1))*SIN((I-1)*DAL)
          call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
          do k=1,3
             af(k)=-bend%scale*af(k)
             az(k)=-bend%scale*az(k)
             ay(k)=-bend%scale*ay(k)
          enddo
          call xp_to_px_r(x,rho,hc,af,bend%p)

          !(pos,b,i,x,y,af,az,ay)

          call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
          ! original formula for electrons (q<0)  now switch to proton

          do k=1,3
             af(k)=-bend%scale*af(k)
             az(k)=-bend%scale*az(k)
             ay(k)=-bend%scale*ay(k)
          enddo


          ! fxhr(f,x,a,az,ay,r,hcurv,p)
          call fx(k1,x,af,az,ay,rho,hc,bend%p)
          if(.not.CHECK_STABLE) return

          if(bend%symplectic) then
             do j=1,3
                xt(2*j)=x(2*j)+dz*k1(2*j)         ! temporary
                xt(2*j-1)=x(2*j-1)
             enddo


             norm0=c_1d10
             doit=.true.

             do it=1,20
                call fx(k1,m,xt,af,az,ay,rho,hc,bend%p)
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
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A72))'
             w_p%c(1)=  " Took more than 20 iterations "
             call write_e(100)
100          continue

             call fx(k1,xt,af,az,ay,rho,hc,bend%p)
             if(.not.check_stable) return

             do j=1,3
                x(2*j-1)=x(2*j-1)+dz*k1(2*j-1)         ! temporary
                x(2*j) =xt(2*j)
             enddo
          else
             do j=1,6
                x(j)=x(j)+dz*k1(j)         ! temporary
             enddo
          endif
          if(che(x,bend)) then
             check_stable=.false.
             return
          endif
          !        tot=tot+bend%dz

          call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
          do k=1,3
             af(k)=-bend%scale*af(k)
             az(k)=-bend%scale*az(k)
             ay(k)=-bend%scale*ay(k)
          enddo
          call px_to_xp_r(x,rho,hc,af,bend%p)
       ENDDO
    enddo
    BEND%X(1,bend%d%ns)=(X(1))*COS((bend%d%ns-1)*DAL)-rhob0
    BEND%X(2,bend%d%ns)=(X(1))*SIN((bend%d%ns-1)*DAL)

    !write(6,*) tot
    i=bend%d%ns-1
    call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
    do k=1,3
       af(k)=-bend%scale*af(k)
       az(k)=-bend%scale*az(k)
       ay(k)=-bend%scale*ay(k)
    enddo
    call xp_to_px_r(x,rho,hc,af,bend%p)


    x(1)=x(1)-rho


  end subroutine track_interpol_1h_r

  subroutine track_interpol_1h_p(bend,x)
    implicit none
    type(fitted_magnetp), intent(inOUT):: bend
    type(real_8), intent(inout)::x(6)
    type(real_8) k1(6),af(3),az(3),ay(3),dif(2),t(2),xt(6)
    real(dp) hc,rho,dz
    real(dp) m(2,2),mi(2,2),det
    integer  i,j,k,NF,KNF,it
    real(dp) norm, norm0
    logical(lp) doit


    call alloc(k1,6);call alloc(af,3);call alloc(az,3);call alloc(ay,3);
    call alloc(dif,2);call alloc(t,2);call alloc(xt,6);

    hc=one/bend%D%X(1)
    rho=bend%D%X(1)
    !  rotating
    nf=bend%p%nst
    x(1)=x(1)+rho

    !  going into non-canonical coordinates
    call get_bdata(bend%d,1,x(1),x(3),af,az,ay)
    do k=1,3
       af(k)=-bend%scale*af(k)
       az(k)=-bend%scale*az(k)
       ay(k)=-bend%scale*ay(k)
    enddo
    call px_to_xp_r(x,rho,hc,af,bend%p)



    dz=(bend%d%dtheta*DEG_TO_RAD_)*rho/nf


    !tot=zero
    do i=1,bend%d%ns-1
       DO KNF=1,NF
          call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
          do k=1,3
             af(k)=-bend%scale*af(k)
             az(k)=-bend%scale*az(k)
             ay(k)=-bend%scale*ay(k)
          enddo
          call xp_to_px_r(x,rho,hc,af,bend%p)

          !(pos,b,i,x,y,af,az,ay)

          call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
          ! original formula for electrons (q<0)  now switch to proton

          do k=1,3
             af(k)=-bend%scale*af(k)
             az(k)=-bend%scale*az(k)
             ay(k)=-bend%scale*ay(k)
          enddo


          ! fxhr(f,x,a,az,ay,r,hcurv,p)
          !          call fx(k1,x,af,az,ay,rho,hc,bend%p)

          !           do j=1,6
          !             x(j)=x(j)+dz*k1(j)         ! temporary
          !          enddo
          call fx(k1,x,af,az,ay,rho,hc,bend%p)

          if(bend%symplectic) then

             do j=1,3
                xt(2*j)=x(2*j)+dz*k1(2*j)         ! temporary
                xt(2*j-1)=x(2*j-1)
             enddo


             norm0=c_1d10
             doit=.true.

             do it=1,20
                call fx(k1,m,xt,af,az,ay,rho,hc,bend%p)
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
             w_p=0
             w_p%nc=1
             w_p%fc='(1(1X,A72))'
             w_p%c(1)=  " Took more than 20 iterations "
             call write_e(100)
100          continue

             call fx(k1,xt,af,az,ay,rho,hc,bend%p)

             do j=1,3
                x(2*j-1)=x(2*j-1)+dz*k1(2*j-1)         ! temporary
                x(2*j) =xt(2*j)
             enddo

          else
             do j=1,6
                x(j)=x(j)+dz*k1(j)         ! temporary
             enddo
          endif


          !        tot=tot+bend%dz

          call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
          do k=1,3
             af(k)=-bend%scale*af(k)
             az(k)=-bend%scale*az(k)
             ay(k)=-bend%scale*ay(k)
          enddo
          call px_to_xp_r(x,rho,hc,af,bend%p)
       ENDDO
    enddo

    !write(6,*) tot
    i=bend%d%ns-1
    call get_bdata(bend%d,i,x(1),x(3),af,az,ay)
    do k=1,3
       af(k)=-bend%scale*af(k)
       az(k)=-bend%scale*az(k)
       ay(k)=-bend%scale*ay(k)
    enddo
    call xp_to_px_r(x,rho,hc,af,bend%p)


    x(1)=x(1)-rho


    call kill(dif,2);call kill(t,2);call kill(xt,6);
    call kill(k1,6);call kill(af,3);call kill(az,3);call kill(ay,3);

  end subroutine track_interpol_1h_p

  subroutine get_adatag(b,i,posx,posy)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i,posx,posy
    integer j,ix,iy
    real(dp) db(3,2),x0,y0,rho
    real(dp) rm,rp,r0,b0(3)
    real(dp) db2(3,3)

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
    r0=x0/rho

    !  first order
    if(iabs(b%order)==1) then
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
             call write_e(2234)
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
             call write_e(2235)
          endif
       endif

       do j=1,2
          !          db(j,1)=(rp*b%b(j,i,ix+1,iy)-rm*b%b(j,i,ix-1,iy))/b%dx/2
          !          db(j,2)=(r0*b%b(j,i,ix,iy+1)-r0*b%b(j,i,ix,iy-1))/b%dy/2
          db(j,1)=(rp*b%b(j,i,ix+1,iy)-r0*b%b(j,i,ix,iy))/b%dx
          db(j,2)=(r0*b%b(j,i,ix,iy+1)-r0*b%b(j,i,ix,iy))/b%dy
       enddo
       j=3
       !          db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix-1,iy))/b%dx/2
       !          db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2
       db(j,1)=(rp*b%b(j,i,ix+1,iy)-r0*b%b(j,i,ix,iy))/b%dx
       db(j,2)=(r0*b%b(j,i,ix,iy+1)-r0*b%b(j,i,ix,iy))/b%dy

       do j=1,2
          b0(j)=r0*b%b(j,i,ix,iy)
       enddo
       j=3
       b0(j)=b%b(j,i,ix,iy)
    elseif(iabs(b%order)==2) then
       ! second order
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
             call write_e(2234)
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
             call write_e(2235)
          endif
       endif


       do j=1,2
          db(j,1)=(rp*b%b(j,i,ix+1,iy)-rm*b%b(j,i,ix-1,iy))/b%dx/2
          db(j,2)=r0*(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2
       enddo
       j=3
       db(j,1)=(b%b(j,i,ix+1,iy)-b%b(j,i,ix-1,iy))/b%dx/2
       db(j,2)=(b%b(j,i,ix,iy+1)-b%b(j,i,ix,iy-1))/b%dy/2

       do j=1,2
          db2(j,1)=(rp*b%b(j,i,ix+1,iy)+rm*b%b(j,i,ix-1,iy)-2*r0*b%b(j,i,ix,iy))/b%dx**2/2
          db2(j,2)=r0*(b%b(j,i,ix,iy+1)+b%b(j,i,ix,iy-1)-2*b%b(j,i,ix,iy))/b%dy**2/2
       enddo
       j=3
       db2(j,1)=(b%b(j,i,ix+1,iy)+b%b(j,i,ix-1,iy)-2*b%b(j,i,ix,iy))/b%dx**2/2
       db2(j,2)=(b%b(j,i,ix,iy+1)+b%b(j,i,ix,iy-1)-2*b%b(j,i,ix,iy))/b%dy**2/2
       do j=1,2
          db2(j,3)=(rp*b%b(j,i,ix+1,iy+1)+rm*b%b(j,i,ix-1,iy-1)-2*r0*b%b(j,i,ix,iy))/2
          db2(j,3)=(db2(j,3)-db2(j,1)*b%dx**2-db2(j,2)*b%dy**2)/b%dx/b%dy
       enddo
       j=3
       db2(j,3)=(b%b(j,i,ix+1,iy+1)+b%b(j,i,ix-1,iy-1)-2*b%b(j,i,ix,iy))/2
       db2(j,3)=(db2(j,3)-db2(j,1)*b%dx**2-db2(j,2)*b%dy**2)/b%dx/b%dy

       do j=1,2
          b0(j)=r0*b%b(j,i,ix,iy)
       enddo
       j=3
       b0(j)=b%b(j,i,ix,iy)

    endif


    ! vector potential

    b%a(1,i,ix,iy)=b0(2)
    b%a(2,i,ix,iy)=db(2,1)
    b%a(3,i,ix,iy)=db(1,1)
    b%a(4,i,ix,iy)=b0(1)
    b%a(5,i,ix,iy)=db(1,2)
    !       b%a(6,i,ix,iy)=zero
    b%a(7,i,ix,iy)=b0(3)
    b%a(8,i,ix,iy)=db(3,1)
    b%a(9,i,ix,iy)=db(3,2)


  end subroutine get_adatag

  subroutine get_adata2c_r(b,i,x,y,af,az,ax)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i
    real(dp), intent(in) :: x,y
    real(dp), intent(out) :: af(3),az(3),ax(3)
    integer ix,iy
    real(dp) posx,posy,x0,y0,rho,dx,dy
    !  second order

    rho=b%X(1)
    posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    !  interpolate
    ix=int(posx)
    iy=int(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
    !write(6,*) "R: x0,y0",x0,y0

    !  first order

    if(iabs(b%order)==1) then
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
             call write_e(2234)
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
             call write_e(2235)
          endif
       endif
    elseif(iabs(b%order)==2) then
       ! second order
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
             call write_e(2234)
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
             call write_e(2234)
          endif
       endif

    endif
    az=zero;af=zero;ax=zero;
    dx=x-x0
    dy=y-y0
    az(1)=b%a(3,i,ix,iy)*dy
    az(1)=az(1)-b%a(1,i,ix,iy)-b%a(2,i,ix,iy)*dx
    az(2)=b%a(4,i,ix,iy)+b%a(3,i,ix,iy)*dx+b%a(5,i,ix,iy)*dy
    af(1)=-b%a(7,i,ix,iy)*dy-b%a(8,i,ix,iy)*dx*dy-b%a(9,i,ix,iy)*dy**2/two
    ax(2)=-b%a(7,i,ix,iy)-b%a(8,i,ix,iy)*dx-b%a(9,i,ix,iy)*dy



  end subroutine get_adata2c_r

  subroutine get_adata2c_p(b,i,xx,yy,af,az,ax)
    implicit none
    type(DATA_MAGNET) b
    integer, intent(in) :: i
    type(real_8), intent(in) :: xx,yy
    type(real_8), intent(out) :: af(3),az(3),ax(3)
    integer ix,iy
    real(dp) posx,posy,x0,y0,x,y,rho
    type(real_8) bf(3),dx,dy
    !  second order
    call alloc(bf,3)
    call alloc(dx,dy)

    rho=b%X(1)
    x=xx
    y=yy
    posx= (x-b%x(1))/(b%x(2)-b%x(1))*(b%nx-1)+one
    posy= (y-b%y(1))/(b%y(2)-b%y(1))*(b%ny-1)+one

    !  interpolate
    ix=int(posx)
    iy=int(posy)
    x0=REAL((ix-1),kind=DP)/REAL((b%nx-1),kind=DP)*(b%x(2)-b%x(1))+b%x(1)
    y0=REAL((iy-1),kind=DP)/REAL((b%ny-1),kind=DP)*(b%y(2)-b%y(1))+b%y(1)
    !write(6,*) "R: x0,y0",x0,y0

    !  first order

    if(iabs(b%order)==1) then
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
             call write_e(2234)
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
             call write_e(2234)
          endif
       endif
    elseif(iabs(b%order)==2) then
       ! second order
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
             call write_e(2234)
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
             call write_e(2234)
          endif
       endif

    endif



    ! vector potential
    af(1)=zero; af(2)=zero; af(3)=zero;
    az(1)=zero; az(2)=zero; az(3)=zero;
    ax(1)=zero; ax(2)=zero; ax(3)=zero;


    dx=xx-x0
    dy=yy-y0
    az(1)=b%a(3,i,ix,iy)*dy
    az(1)=az(1)-b%a(1,i,ix,iy)-b%a(2,i,ix,iy)*dx
    az(2)=b%a(4,i,ix,iy)+b%a(3,i,ix,iy)*dx+b%a(5,i,ix,iy)*dy
    af(1)=-b%a(7,i,ix,iy)*dy-b%a(8,i,ix,iy)*dx*dy-b%a(9,i,ix,iy)*dy**2/two
    ax(2)=-b%a(7,i,ix,iy)-b%a(8,i,ix,iy)*dx-b%a(9,i,ix,iy)*dy



    call kill(dx,dy)
    call kill(bf,3)

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

    if(x(1)<b%d%x(1)/s) che=.true.
    if(x(1)>s*b%d%x(2)) che=.true.
    if(x(3)<s*b%d%y(1)) che=.true.
    if(x(3)>s*b%d%y(2)) che=.true.

  end function che

end module fitted_MAG
