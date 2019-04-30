module beam_beam_ptc
  !  use orbit_ptc
  use S_TRACKING
  ! use madx_ptc_module
  implicit none
  public
  private BBKICKP, BBKICKR,PATCH_BBR,PATCH_BBP
  private ccperrfP, ccperrfr ,ccperrf
  !  private TRACK_NODE_LAYOUT_FLAG_R,TRACK_NODE_LAYOUT_FLAG_P
  private imax
  integer ::imax=1000

  INTERFACE ccperrf
     MODULE PROCEDURE ccperrfP
     MODULE PROCEDURE ccperrfr
  END INTERFACE


  INTERFACE BBKICK
     MODULE PROCEDURE BBKICKR
     MODULE PROCEDURE BBKICKP
  END INTERFACE
  INTERFACE PATCH_BB
     MODULE PROCEDURE PATCH_BBR
     MODULE PROCEDURE PATCH_BBP
  END INTERFACE

  logical(lp), target :: do_beam_beam= my_false

contains

  SUBROUTINE PATCH_BBR(B,X,k,BETA0,exact,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: B
    real(dp), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: exact,ENTERING
    REAL(DP),INTENT(IN):: BETA0
    REAL(DP) A(3),D(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ENTERING) THEN
       X(3)=B%A_X1*X(3);X(4)=B%A_X1*X(4);
       CALL ROT_YZ(B%A(1),X,BETA0,exact,k%TIME)
       CALL ROT_XZ(B%A(2),X,BETA0,exact,k%TIME)
       CALL ROT_XY(B%A(3),X)  !,exact)
       CALL TRANS(B%D,X,BETA0,exact,k%TIME)
       X(3)=B%A_X2*X(3);X(4)=B%A_X2*X(4);
    ELSE
       A=-B%A
       D=-B%D
       X(3)=B%A_X2*X(3);X(4)=B%A_X2*X(4);
       CALL TRANS(D,X,BETA0,exact,k%TIME)
       CALL ROT_XY(A(3),X)  !,exact)
       CALL ROT_XZ(A(2),X,BETA0,exact,k%TIME)
       CALL ROT_YZ(A(1),X,BETA0,exact,k%TIME)
       X(3)=B%A_X1*X(3);X(4)=B%A_X1*X(4);
    ENDIF


  END SUBROUTINE PATCH_BBR

  SUBROUTINE PATCH_BBP(B,X,k,BETA0,exact,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: B
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: exact,ENTERING
    REAL(DP),INTENT(IN):: BETA0
    REAL(DP) A(3),D(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ENTERING) THEN
       X(3)=B%A_X1*X(3);X(4)=B%A_X1*X(4);
       CALL ROT_YZ(B%A(1),X,BETA0,exact,k%TIME)
       CALL ROT_XZ(B%A(2),X,BETA0,exact,k%TIME)
       CALL ROT_XY(B%A(3),X)  !,exact)
       CALL TRANS(B%D,X,BETA0,exact,k%TIME)
       X(3)=B%A_X2*X(3);X(4)=B%A_X2*X(4);
    ELSE
       A=-B%A
       D=-B%D
       X(3)=B%A_X2*X(3);X(4)=B%A_X2*X(4);
       CALL TRANS(D,X,BETA0,exact,k%TIME)
       CALL ROT_XY(A(3),X)  !,exact)
       CALL ROT_XZ(A(2),X,BETA0,exact,k%TIME)
       CALL ROT_YZ(A(1),X,BETA0,exact,k%TIME)
       X(3)=B%A_X1*X(3);X(4)=B%A_X1*X(4);
    ENDIF


  END SUBROUTINE PATCH_BBP


  subroutine BBKICKR(BB,X)

    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   track a set of particle through a beam-beam interaction region.    *
    !   see mad physicist's manual for the formulas used.                  *
    !input:                                                                *
    ! input/output:                                                        *
    !   b%x(:,6)(double)  track coordinates: (x, px, y, py,  pt,t).        *
    !   b%n    (integer) number of tracks.                                 *
    !----------------------------------------------------------------------*
    real(dp) sx2,sy2,xs,ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,&
         cbx,cby,ten3m,explim,sx,sy,xm,ym
    TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
    REAL(DP), INTENT(INOUT) :: X(6)
    parameter(ten3m=1.0e-3_dp,explim=150.0_dp)
    SX=BB%SX
    SY=BB%SY
    XM=BB%XM
    YM=BB%YM
    FK=BB%FK
    !    write(6,*) "bb%FK = " ,bb%FK

    if (fk == 0.0_dp)  return
    sx2 = sx*sx
    sy2 = sy*sy
    !---- limit formulae for sigma(x) = sigma(y).
    if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then !ten3m = 1.0d-3
       xs = x(1) - xm
       ys = x(3) - ym
       rho2 = xs * xs + ys * ys
       tk = rho2 / (sx2 + sy2)
       if (tk .gt. explim) then
          phix = xs * fk / rho2
          phiy = ys * fk / rho2
       else if (rho2 .ne. 0.0_dp) then
          phix = xs * fk / rho2 * (1.0_dp - exp(-tk) )
          phiy = ys * fk / rho2 * (1.0_dp - exp(-tk) )
       else
          phix = 0.0_dp
          phiy = 0.0_dp
       endif
       phix = phix - bb%bbk(1) ! subtract horizontal bb kick
       phiy = phiy - bb%bbk(2) ! subtract vertical co
       x(2) = x(2) + phix
       x(4) = x(4) + phiy

       !---- case sigma(x) > sigma(y).
    else if (sx2 > sy2) then
       r2 = 2.0_dp * (sx2 - sy2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       xs = x(1) - xm
       ys = x(3) - ym
       xr = abs(xs) / r
       yr = abs(ys) / r
       call ccperrf(xr, yr, crx, cry)
       tk = (xs * xs / sx2 + ys * ys / sy2) / 2.0_dp
       if (tk .gt. explim) then
          phix = rk * cry
          phiy = rk * crx
       else
          xb = (sy / sx) * xr
          yb = (sx / sy) * yr
          call ccperrf(xb, yb, cbx, cby)
          phix = rk * (cry - exp(-tk) * cby)
          phiy = rk * (crx - exp(-tk) * cbx)
       endif
       x(2) = x(2) + phix * sign(1.0_dp,xs)
       x(4) = x(4) + phiy * sign(1.0_dp,ys)
       x(2) = x(2) - BB%bbk(1)
       x(4) = x(4) - BB%bbk(2)

       !---- case sigma(x) < sigma(y).
    else
       r2 = 2.0_dp * (sy2 - sx2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       !       do itrack = 1, b%n
       !         IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       xr = abs(xs) / r
       yr = abs(ys) / r
       call ccperrf(yr, xr, cry, crx)
       tk = (xs * xs / sx2 + ys * ys / sy2) / 2.0_dp
       if (tk .gt. explim) then
          phix = rk * cry
          phiy = rk * crx
       else
          xb  = (sy / sx) * xr
          yb  = (sx / sy) * yr
          call ccperrf(yb, xb, cby, cbx)
          phix = rk * (cry - exp(-tk) * cby)
          phiy = rk * (crx - exp(-tk) * cbx)
       endif

       x(2) = x(2) + phix * sign(1.0_dp,xs)
       x(4) = x(4) + phiy * sign(1.0_dp,ys)

       x(2) = x(2) - BB%bbk(1)
       x(4) = x(4) - BB%bbk(2)

    endif
    !    IF(DZ/=ZERO) CALL DRIFT_BEAM_BACK_TO_POSITION(th,DZ,B)

  end subroutine BBKICKR

  subroutine ccperrfr(xx, yy, wx, wy)
    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   modification of wwerf, double precision complex error function,    *
    !   written at cern by K. Koelbig.                                     *
    ! input:                                                               *
    !   xx, yy    (double)    real + imag argument                         *
    ! output:                                                              *
    !   wx, wy    (double)    real + imag function result                  *
    !----------------------------------------------------------------------*
    integer n,nc,nu
    real(dp), INTENT(INOUT):: xx,yy,wx,wy
    real(dp) x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux,  &
         rx(33),ry(33),cc,xlim,ylim,fac1,fac2,fac3
    parameter(cc=1.12837916709551_dp,        &
         xlim=5.33_dp,ylim=4.29_dp,fac1=3.2_dp,fac2=23.0_dp,fac3=21.0_dp)

    x = abs(xx)
    y = abs(yy)

    if (y .lt. ylim  .and.  x .lt. xlim) then
       q  = (1.0_dp - y / ylim) * sqrt(1.0_dp - (x/xlim)**2)
       h  = 1.0_dp / (fac1 * q)
       nc = 7 + int(fac2*q)
       xl = h**(1 - nc)
       xh = y + 0.5_dp/h
       yh = x
       nu = 10 + int(fac3*q)
       rx(nu+1) = 0.0_dp
       ry(nu+1) = 0.0_dp

       do n = nu, 1, -1
          tx = xh + n * rx(n+1)
          ty = yh - n * ry(n+1)
          tn = tx*tx + ty*ty
          rx(n) = 0.5_dp * tx / tn
          ry(n) = 0.5_dp * ty / tn
       enddo

       sx = 0.0_dp
       sy = 0.0_dp

       do n = nc, 1, -1
          saux = sx + xl
          sx = rx(n) * saux - ry(n) * sy
          sy = rx(n) * sy + ry(n) * saux
          xl = h * xl
       enddo

       wx = cc * sx
       wy = cc * sy
    else
       xh = y
       yh = x
       rx(1) = 0.0_dp
       ry(1) = 0.0_dp

       do n = 9, 1, -1
          tx = xh + n * rx(1)
          ty = yh - n * ry(1)
          tn = tx*tx + ty*ty
          rx(1) = 0.5_dp * tx / tn
          ry(1) = 0.5_dp * ty / tn
       enddo

       wx = cc * rx(1)
       wy = cc * ry(1)
    endif

    !      if(y .eq. zero) wx = exp(-x**2)
    if(yy .lt. 0.0_dp) then
       wx =   2.0_dp * exp(y*y-x*x) * cos(2.0_dp*x*y) - wx
       wy = - 2.0_dp * exp(y*y-x*x) * sin(2.0_dp*x*y) - wy
       if(xx .gt. 0.0_dp) wy = -wy
    else
       if(xx .lt. 0.0_dp) wy = -wy
    endif

  end SUBROUTINE ccperrfr



  subroutine BBKICKP(BB,X)

    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   track a set of particle through a beam-beam interaction region.    *
    !   see mad physicist's manual for the formulas used.                  *
    !input:                                                                *
    ! input/output:                                                        *
    !   b%x(:,6)(double)  track coordinates: (x, px, y, py,  pt,t).      *
    !   b%n    (integer) number of tracks.                              *
    !----------------------------------------------------------------------*
    integer it
    TYPE(REAL_8) xs,ys,tk,phix,phiy,xb,yb,crx,cry
    TYPE(REAL_8) xr,yr,cbx,cby,rho2
    REAL(DP) sx2,sy2,sx,sy,xm,ym,fk,ten3m,explim,xn1,xn2,xs1,xs2,arglim,rk
    REAL(DP) r,r2,n
    TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    parameter(ten3m=1.0e-3_dp,arglim=1.0e-2_dp,explim=150.0_dp)

    if (BB%fk == 0.0_dp)  return

    CALL ALLOC(xr,yr,cbx,cby,rho2)
    CALL ALLOC(xs,ys,tk,phix,phiy)
    CALL ALLOC(xb,yb,crx,cry)


    SX=BB%SX
    SY=BB%SY
    XM=BB%XM
    YM=BB%YM
    FK=BB%FK

    sx2 = sx*sx
    sy2 = sy*sy
    !---- limit formulae for sigma(x) = sigma(y).
    xn1=abs(sx2 - sy2)
    xn2= ten3m * (sx2 + sy2)
    xs1=sx2
    xs2= sy2
    if (xn1 .le.xn2 ) then !ten3m = 1.0d-3
       xs = x(1) - xm
       ys = x(3) - ym
       rho2 = xs * xs + ys * ys
       tk = rho2 / (sx2 + sy2)
       if (tk .gt. explim) then
          phix = xs * fk / rho2
          phiy = ys * fk / rho2
       else if (tk > arglim) then
          phix = xs * fk / rho2 * (1.0_dp - exp(-tk) )
          phiy = ys * fk / rho2 * (1.0_dp - exp(-tk) )
       else

          xr=1.0_dp
          yr=1.0_dp

          n=mybig
          do it=1,imax
             xr=-xr*tk/(it+1)
             yr=yr+xr
             if(it>10)n=full_abs(xr)
             if(n<=puny) exit
          enddo
          if(it>imax-2) then
             write(6,*) it,n
             write(6,*) " Stopped in Beam-Beam "
          endif
          phix = xs * fk / (2.0_dp * sx2) * YR ! fudge
          phiY = Ys * fk / (2.0_dp * sx2) * YR ! fudge
       endif

       phix = phix - bb%bbk(1)
       phiy = phiy - bb%bbk(2)

       x(2) = x(2) + phix
       x(4) = x(4) + phiy


       !---- case sigma(x) > sigma(y).
    else


       xs = x(1) - xm
       ys = x(3) - ym
       tk = (xs * xs / sx2 + ys * ys / sy2) / 2.0_dp

       if(xs1 > xs2) then
          r2 = 2.0_dp * (sx2 - sy2)
          r  = sqrt(r2)
          rk = fk * sqrt(pi) / r
          xr = xs / r   !
          yr = ys / r   !
          
          if( (xr.sub.'0') < 0) then
            xr = -xr
          endif
          if( (yr.sub.'0') < 0) then
            yr = -yr
          endif
          
          call ccperrf(xr, yr, crx, cry)
          if (tk .gt. explim) then
             phix = rk * cry
             phiy = rk * crx
          else
             xb = (sy / sx) * xr
             yb = (sx / sy) * yr
             call ccperrf(xb, yb, cbx, cby)
             phix = rk * (cry - exp(-tk) * cby)
             phiy = rk * (crx - exp(-tk) * cbx)
          endif
          x(2) = x(2) + phix
          x(4) = x(4) + phiy
          !          if (.NOT.bborbit)  then
          x(2) = x(2) - BB%bbk(1)
          x(4) = x(4) - BB%bbk(2)
          !          endif
          !       enddo

          !---- case sigma(x) < sigma(y).
       else
          r2 = 2.0_dp * (sy2 - sx2)
          r  = sqrt(r2)
          rk = fk * sqrt(pi) / r

          !       xs = x(1) - xm
          !       ys = x(3) - ym
          xr = xs / r !abs
          yr = ys / r !abs
          if( (xr.sub.'0') < 0) then
            xr = -xr
          endif
          if( (yr.sub.'0') < 0) then
            yr = -yr
          endif

          call ccperrf(yr, xr, cry, crx)
          !       tk = (xs * xs / sx2 + ys * ys / sy2) / two
          if (tk .gt. explim) then
             phix = rk * cry
             phiy = rk * crx
          else
             xb  = (sy / sx) * xr
             yb  = (sx / sy) * yr
             call ccperrf(yb, xb, cby, cbx)
             phix = rk * (cry - exp(-tk) * cby)
             phiy = rk * (crx - exp(-tk) * cbx)
          endif
          !          track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
          !          track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
          x(2) = x(2) + phix
          x(4) = x(4) + phiy
          !          if (.NOT.bborbit)  then
          !--- subtract closed orbit kick
          !            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
          !            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
          x(2) = x(2) - BB%bbk(1)
          x(4) = x(4) - BB%bbk(2)
          !          endif
          !       enddo
       endif
    endif

    CALL KILL(xr,yr,cbx,cby,rho2)
    CALL KILL(xs,ys,tk,phix,phiy)
    CALL KILL(xb,yb,crx,cry)

  end subroutine BBKICKP

  subroutine ccperrfP(xx, yy, wx, wy)
    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   modification of wwerf, double precision complex error function,    *
    !   written at cern by K. Koelbig.                                     *
    ! input:                                                               *
    !   xx, yy    (double)    real + imag argument                         *
    ! output:                                                              *
    !   wx, wy    (double)    real + imag function result                  *
    !----------------------------------------------------------------------*
    TYPE(REAL_8),INTENT(INOUT):: xx,yy,wx,wy
    TYPE(complex_8) z,zt,w
    complex(dp) z0,w0,w1,wt0
    real(dp) xx0, yy0, wx0, wy0
    integer i
    call alloc( z,zt,w)

    z=xx+i_*yy
    z0=z
    z=z-z0

    xx0=real(z0)
    yy0=aimag(z0)
    call ccperrf(xx0, yy0, wx0, wy0)

    w0=wx0+i_*wy0

    w1=-2.0_dp*z0*w0+2.0_dp*i_/sqrt(pi)

    w=w0+w1*z

    zt=z

    do i=2,c_%no

       zt=z*zt
       wt0=  -2.0_dp*(w0+z0*w1)/i

       w=w+wt0*zt

       w0=w1
       w1=wt0

    enddo

    wx=real(w)
    wy=aimag(w)


    call kill( z,zt,w)

  end SUBROUTINE ccperrfP
end module  beam_beam_ptc
