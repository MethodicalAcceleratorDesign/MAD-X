!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90
module fitted_MAG_1
  USE S_def_all_kinds
  implicit none
  public
  integer, parameter :: dpm =dp
  integer mff
  private fxr
  logical(lp) :: use_beamlet=.false.
  TYPE magnet_b
     INTEGER, POINTER :: nb(:)
     real(dp), POINTER :: blim(:,:),hb,x_c
     real(dpm), POINTER :: b(:,:,:,:)
     real(dp), pointer :: grid(:,:,:)
     real(dp),pointer :: s(:)
     real(dp),pointer :: ds,dsh,scale
     real(dp), pointer :: dx,dy
  END TYPE magnet_b

  INTERFACE get_b
     MODULE PROCEDURE get_b_int_r
     !     MODULE PROCEDURE get_b_int_p
  END INTERFACE

  INTERFACE fxb
     MODULE PROCEDURE fxr
     !     MODULE PROCEDURE fxp
  END INTERFACE




contains

  subroutine pointer_magnet_b(mag_b,hb,x_c,nb,blim)
    implicit none
    type(magnet_b), intent(inout) :: mag_b
    integer, intent(in):: nb(3)
    real(dpm), intent(in):: hb,blim(3,2),x_c

    allocate(mag_b%nb(3),mag_b%blim(3,2),mag_b%hb,mag_b%x_c)

    mag_b%blim=blim
    mag_b%nb=nb
    mag_b%hb=hb
    mag_b%x_c=x_c

    allocate(mag_b%b(3,mag_b%nb(1),mag_b%nb(2),mag_b%nb(3)),mag_b%s(nb(3)))

    allocate(mag_b%ds,mag_b%dsh,mag_b%scale,mag_b%grid(mag_b%nb(1),mag_b%nb(2),2))
    allocate(mag_b%dx,mag_b%dy)

  end subroutine pointer_magnet_b

  subroutine get_b_int_r(b,i,x,y,bf)
    implicit none
    type(magnet_b) b
    integer, intent(in) :: i
    real(dp), intent(in) :: x,y
    real(dp), intent(out) :: bf(3)
    integer ix,iy,k,j,l
    real(dp) posx,posy,x0,y0,dx,dy,dz,f(3,-2:2,-2:2),sb

    posx= (x-b%blim(1,1))/(b%blim(1,2)-b%blim(1,1))*(b%nb(1)-1)+one
    posy= (y-b%blim(2,1))/(b%blim(2,2)-b%blim(2,1))*(b%nb(2)-1)+one

    if(.not.use_beamlet) then



       ix=int(posx)
       iy=int(posy)
       x0=REAL((ix-1),kind=DP)/REAL((b%nb(1)-1),kind=DP)*(b%blim(1,2)-b%blim(1,1))+b%blim(1,1)
       y0=REAL((iy-1),kind=DP)/REAL((b%nb(2)-1),kind=DP)*(b%blim(2,2)-b%blim(2,1))+b%blim(2,1)

       !  first order


       if(posx<3.or.posx>b%nb(1)-2) then
          ix=3
          if(posx>3) ix=b%nb(1)-3
          x0=REAL((ix-1),kind=DP)/REAL((b%nb(1)-1),kind=DP)*(b%blim(1,2)-b%blim(1,1))+b%blim(1,1)
       endif

       if(posy<3.or.posy>b%nb(2)-2) then
          iy=3
          if(posy>3) iy=b%nb(2)-3
          y0=REAL((iy-1),kind=DP)/REAL((b%nb(2)-1),kind=DP)*(b%blim(2,2)-b%blim(2,1))+b%blim(2,1)
       endif

       do k=1,3
          do j=-2,2
             do l=-2,2
                !        write(mff,*) k,ix+j,iy+l,i

                f(k,j,l)=b%b(k,ix+j,iy+l,i)
             enddo
          enddo
       enddo

       !az(n) grad az
       !af(n)  a
       !ax(n)  grad ax

       !    write(6,*) size(b%a,1),size(b%a,2),size(b%a,3),size(b%a,4)
       !    write(6,*) i,ix,iy
       !   pause 333
       bf=zero;
       dx=(x-x0)/(b%blim(1,2)-b%blim(1,1))*(b%nb(1)-1)
       dy=(y-y0)/(b%blim(2,2)-b%blim(2,1))*(b%nb(2)-1)
       dz=zero
       do k=1,3
          ! bf(k)=dx*(dx-one)/two*f(K,-1,0)+dy*(dy-one)/two*f(K,0,-1)+(ONE+DX*DY-DX**2-DY**2)*F(K,0,0)
          ! bf(k)=  bf(k)+DX*(DX-TWO*DY+one)/TWO*F(K,1,0)+DY*(DY-TWO*DX+one)/TWO*F(K,0,1)+DX*DY*F(K,1,1)

          ! bf(k)=(one-dx-dy)*f(k,0,0)+dx*f(k,1,0)+dy*f(k,0,1)
          bf(k)=(one-dx**2-dy**2)*f(k,0,0)+f(k,1,0)*(4*dx+3*dx**2-dx**3)/6+f(k,-1,0)*(-4*dx+3*dx**2+dx**3)/6
          bf(k)= bf(k) +f(k,2,0)*(dx**3-dx)/12-f(k,-2,0)*(dx**3-dx)/12
          bf(k)=bf(k) +f(k,0,1)*(4*dy+3*dy**2-dy**3)/6+f(k,0,-1)*(-4*dy+3*dy**2+dy**3)/6
          bf(k)= bf(k) +f(k,0,2)*(dy**3-dy)/12-f(k,0,-2)*(dy**3-dy)/12

          ! makes less symplectic  xy term only h**3 accurate
          bf(k)=bf(k)+dx*dy*(2*f(k,0,0)+f(k,1,1)+f(k,-1,-1)-f(k,0,1)-f(k,1,0)-f(k,0,-1)-f(k,-1,0))/two
          ! correction of f_xx and f_yy  ! f_xxxx and f_yyyy
          bf(k)=bf(k)+(dx**2-dx**4)*(-6*f(k,0,0)+4*f(k,1,0)+4*f(k,-1,0)-f(k,2,0)-f(k,-2,0))/c_24
          bf(k)=bf(k)+(dy**2-dy**4)*(-6*f(k,0,0)+4*f(k,0,1)+4*f(k,0,-1)-f(k,0,2)-f(k,0,-2))/c_24


       enddo
    else   ! beamlets
       sb=c_1_8



       bf=zero
       do ix=1,b%nb(1)
          do iy=1,b%nb(2)

             dz=exp(-(x-b%grid(ix,iy,1))**2/b%dx**2/sb**2-(y-b%grid(ix,iy,2))**2/b%dy**2/sb**2)/pi/sb**2
             bf(1)=real(b%b(1,ix,iy,i),kind=dp)*dz + bf(1)
             bf(2)=real(b%b(2,ix,iy,i),kind=dp)*dz + bf(2)
             bf(3)=real(b%b(3,ix,iy,i),kind=dp)*dz + bf(3)
          enddo
       enddo
    endif

    do k=1,3
       bf(k)=-b%scale*bf(k)
    enddo

  end subroutine get_b_int_r

  subroutine read_magnet_b(mag_b,filename)
    implicit none
    integer mf
    integer nb(3)
    real(dp) hb,blim(3,2),x,y,x_c
    character*(*) filename
    character*255 line
    type(magnet_b) mag_b
    integer i,j,k

    call kanalnummer(mf)

    open(unit=mf, file=filename)
    read(mf,'(a255)') line

    write(6,*) line
    read(line,*) x_c, hb,nb(3),blim(3,1),blim(3,2),nb(1),blim(1,1),blim(1,2),nb(2),blim(2,1),blim(2,2)
    write(6,*)x_c, hb,nb(3),blim(3,1),blim(3,2),nb(1),blim(1,1),blim(1,2),nb(2),blim(2,1),blim(2,2)

    if(mod(nb(3),2)/=1) then
       write(6,*) " Need a odd number of slices in fitted map"
       stop 1000
    endif

    call  pointer_magnet_b(mag_b,hb,x_c,nb,blim)

    do k=1,mag_b%nb(3)
       read(mf,*)   mag_b%s(k)
       do i=1,mag_b%nb(1)
          do j=1,mag_b%nb(2)
             read(mf,*) x,y,mag_b%b(1,i,j,k),mag_b%b(2,i,j,k),mag_b%b(3,i,j,k)
             if(k==1) then
                mag_b%grid(i,j,1)=x
                mag_b%grid(i,j,2)=y
             endif
          enddo
       enddo
    enddo

    close(mf)

    mag_b%ds = (blim(3,2)-blim(3,1))/(mag_b%nb(3)-1)
    mag_b%dx = (blim(1,2)-blim(1,1))/(mag_b%nb(1)-1)
    mag_b%dy = (blim(2,2)-blim(2,1))/(mag_b%nb(2)-1)
    mag_b%dsh=mag_b%ds /two
    mag_b%scale = one
  end subroutine read_magnet_b


  subroutine fxr(f,x,b,hcurv,BETA0_in,GAMMA0I_in,time)
    implicit none

    real(dp)  d(3),c(6),BETA0,GAMMA0I     ! variable r is a distance not necessarily hcurv**-1
    real(dp) ,intent(in) :: b(3),hcurv
    real(dp), intent(in) :: BETA0_in,GAMMA0I_in
    real(dp) ,intent(inout) :: x(6)
    real(dp), intent(out):: f(6)
    logical(lp), intent(in) ::  time

    if(time) then
       beta0=BETA0_in;GAMMA0I=GAMMA0I_in;
    else
       beta0=one;GAMMA0I=zero;
    endif

    ! using absolute coordinate
    !    d(1)=root(x(2)**2+x(4)**2+(one+hcurv*x(1))**2)
    if(hcurv==zero) then
       d(1)=root(x(2)**2+x(4)**2+one)
    else
       d(1)=root(x(2)**2+x(4)**2+(hcurv*x(1))**2)
    endif
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


    d(2)=gamma0I/beta0/d(2)
    f(6)=root((1+d(2)**2))*d(1)  ! (time)-prime = dt/dz
  end subroutine fxr

  !

  subroutine track_4r(bend,x,BETA0_in,GAMMA0I_in,time)
    implicit none
    type(magnet_b), intent(inOUT):: bend
    real(dp) , intent(inout)::x(6)
    real(dp) , intent(in)::BETA0_in,GAMMA0I_in
    logical(lp), intent(in) :: time
    real(dp) k1(6),k2(6),k3(6),k4(6),xt(6),hc,rho,dz,DAL,ds
    real(dp) bf(3),XI,ZETA,rhob0
    integer  i,j,ns
    real(dp) theta_half

    IF(.NOT.CHECK_STABLE) return


    ds=bend%ds
    ns=bend%nb(3)
    hc=bend%hb
    x(1)=x(1)+bend%x_c
    do i=1,ns-1,2
       dz=two*ds

       call get_b(bend,i,x(1),x(3),bf)
       ! original formula for electrons (q<0)  now switch to proton

       call fxb(k1,x,bf,hc,BETA0_in,GAMMA0I_in,time)
       do j=1,6
          xt(j)=x(j)+dz*k1(j)/two         ! temporary
       enddo

       call get_b(bend,i+1,xt(1),xt(3),bf)
       call fxb(k2,xt,bf,hc,BETA0_in,GAMMA0I_in,time)
       do j=1,6
          xt(j)=x(j)+dz*k2(j)/two         ! temporary
       enddo

       call get_b(bend,i+1,xt(1),xt(3),bf)
       call fxb(k3,xt,bf,hc,BETA0_in,GAMMA0I_in,time)
       do j=1,6
          xt(j)=x(j)+dz*k3(j)            ! temporary
       enddo

       call get_b(bend,i+2,xt(1),xt(3),bf)
       call fxb(k4,xt,bf,hc,BETA0_in,GAMMA0I_in,time)
       write(17,'(i4,4(1x,E15.8))') i+2,bf(2),xt(1)
       do j=1,6
          x(j)=x(j)+dz*(k1(j)+two*k2(j)+two*k3(j)+k4(j))/six     ! temporary
       enddo

       if(.not.CHECK_STABLE) return


       !          if(che(x,bend)) then
       !             check_stable=.false.
       !             return
       !          endif
       write(mff,*) bend%s(i),x(1)

    ENDDO

    x(1)=x(1)-bend%x_c




  end subroutine track_4r



end module fitted_MAG_1
