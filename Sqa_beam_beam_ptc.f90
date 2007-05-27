module beam_beam_ptc
  use orbit_ptc
  ! use madx_ptc_module
  implicit none
  public
  private BBKICKP, BBKICKR !,BBKICK
  private ccperrfP, ccperrfr,ccperrf
  private TRACK_NODE_LAYOUT_FLAG_R,TRACK_NODE_LAYOUT_FLAG_P

  INTERFACE ccperrf
     MODULE PROCEDURE ccperrfP
     MODULE PROCEDURE ccperrfr
  END INTERFACE

  INTERFACE BBKICK
     MODULE PROCEDURE BBKICKR
     MODULE PROCEDURE BBKICKP
  END INTERFACE

  INTERFACE track_x_bb
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_R
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_P
  END INTERFACE



contains

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
    logical(lp) bborbit
    integer itrack
    real(dp) sx2,sy2,xs,ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,&
         bbpar,cbx,cby,ten3m,explim,sx,sy,xm,ym,DZ
    TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
    REAL(DP), INTENT(INOUT) :: X(6)
    parameter(ten3m=1.0e-3_dp,explim=150.0_dp)
    SX=BB%SX
    SY=BB%SY
    XM=BB%XM
    YM=BB%YM
    FK=BB%FK
    write(6,*) "bb%FK = " ,bb%FK
    !    CALL BEAM_BEAM(B,TH%BT,sx,sy,xm,ym,bbpar,BBORBIT,MY_TRUE)
    !---- initialize.
    !      bborbit = get_option('bborbit ') .ne. 0
    !      pi=get_variable('pi ')   ! value of pi
    !      sx = node_value('sigx ') ! Sigma-x
    !      sy = node_value('sigy ') ! Sigma-y
    !      xm = node_value('xma ')  ! x-offset between beams
    !      ym = node_value('yma ')  ! y-offset between beams
    !    fk = two * bbpar !parvec(5) * parvec(6) / parvec(7)
    !parvec(5)=get_value('probe ', 'arad ')
    !arad = 1.e-16 * charge * charge * get_variable("qelect")
    !qelect=1.602176462e-19
    !parvec(6)=node_value('charge ') * get_value('probe ', 'npart ')
    !parvec(7)=get_value('probe ','gamma ')
    !Yes, it is all in a_scratch_size
    !
    !parvec(5) = c_1d_16 * charge * charge * qelect
    ! The "npart" is the number of particles in the other beam and has to be  passed
    !parvec(6) = charge * get_value('probe ', 'npart ')
    !parvec(7) = gamma


    if (fk == zero)  return
    !    DZ=-TH%BT%DX(3)
    !    IF(DZ/=ZERO) CALL DRIFT_BEAM_BACK_TO_POSITION(th,DZ,B)
    !    DZ=-DZ
    !    if (.not. bborbit)  then
    !What counts is that if (.not. bborbit) then subtract closed orbit THIS IS THE DEFAULT
    !closed is tricky and really should be treated selfconsistently, I.e. out of reach here
    !--- find position of closed orbit bb_kick
    !        do ipos = 1, bbd_cnt
    !          if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
    !        enddo
    !        ipos = 0
    !    1   continue
    !    endif
    sx2 = sx*sx
    sy2 = sy*sy
    !---- limit formulae for sigma(x) = sigma(y).
    if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then !ten3m = 1.0d-3
       !       do itrack = 1, b%n
       !          IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       rho2 = xs * xs + ys * ys
       tk = rho2 / (two * sx2)
       if (tk .gt. explim) then
          phix = xs * fk / rho2
          phiy = ys * fk / rho2
       else if (rho2 .ne. zero) then
          phix = xs * fk / rho2 * (one - exp(-tk) )
          phiy = ys * fk / rho2 * (one - exp(-tk) )
       else
          phix = zero
          phiy = zero
       endif
       !          if (.NOT.bborbit)  then
       !--- subtract closed orbit kick
       phix = phix - bb%bbk(1) ! subtract horizontal bb kick
       phiy = phiy - bb%bbk(2) ! subtract vertical co
       !          endif
       x(2) = x(2) + phix
       x(4) = x(4) + phiy
       !       enddo

       !---- case sigma(x) > sigma(y).
    else if (sx2 > sy2) then
       r2 = two * (sx2 - sy2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       !       do itrack = 1, b%n
       !          IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       xr = abs(xs) / r
       yr = abs(ys) / r
       call ccperrf(xr, yr, crx, cry)
       tk = (xs * xs / sx2 + ys * ys / sy2) / two
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
       !          track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
       !          track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
       x(2) = x(2) + phix * sign(one,xs)
       x(4) = x(4) + phiy * sign(one,ys)
       !          if (.NOT.bborbit)  then
       !--- subtract closed orbit kick
       !            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
       !            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
       x(2) = x(2) - BB%bbk(1)
       x(4) = x(4) - BB%bbk(2)
       !          endif
       !       enddo

       !---- case sigma(x) < sigma(y).
    else
       r2 = two * (sy2 - sx2)
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
       tk = (xs * xs / sx2 + ys * ys / sy2) / two
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
       x(2) = x(2) + phix * sign(one,xs)
       x(4) = x(4) + phiy * sign(one,ys)
       !          if (.NOT.bborbit)  then
       !--- subtract closed orbit kick
       !            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
       !            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
       x(2) = x(2) - BB%bbk(1)
       x(4) = x(4) - BB%bbk(2)
       !          endif
       !       enddo
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
    real(dp) xx,yy,wx,wy,x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux,&
         &rx(33),ry(33),cc,xlim,ylim,fac1,fac2,fac3
    parameter(cc=1.12837916709551_dp,        &
         xlim=5.33_dp,ylim=4.29_dp,fac1=3.2_dp,fac2=23.0_dp,fac3=21.0_dp)

    x = abs(xx)
    y = abs(yy)

    if (y .lt. ylim  .and.  x .lt. xlim) then
       q  = (one - y / ylim) * sqrt(one - (x/xlim)**2)
       h  = one / (fac1 * q)
       nc = 7 + int(fac2*q)
       xl = h**(1 - nc)
       xh = y + half/h
       yh = x
       nu = 10 + int(fac3*q)
       rx(nu+1) = zero
       ry(nu+1) = zero

       do n = nu, 1, -1
          tx = xh + n * rx(n+1)
          ty = yh - n * ry(n+1)
          tn = tx*tx + ty*ty
          rx(n) = half * tx / tn
          ry(n) = half * ty / tn
       enddo

       sx = zero
       sy = zero

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
       rx(1) = zero
       ry(1) = zero

       do n = 9, 1, -1
          tx = xh + n * rx(1)
          ty = yh - n * ry(1)
          tn = tx*tx + ty*ty
          rx(1) = half * tx / tn
          ry(1) = half * ty / tn
       enddo

       wx = cc * rx(1)
       wy = cc * ry(1)
    endif

    !      if(y .eq. zero) wx = exp(-x**2)
    if(yy .lt. zero) then
       wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx
       wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy
       if(xx .gt. zero) wy = -wy
    else
       if(xx .lt. zero) wy = -wy
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
    logical(lp) bborbit
    integer itrack
    TYPE(REAL_8) sx2,sy2,xs,ys,rho2,tk,phix,phiy,rk,xb,yb,crx,cry
    TYPE(REAL_8) xr,yr,r,r2,bbpar,cbx,cby,DZ
    REAL(DP) sx,sy,xm,ym,fk,ten3m,explim
    TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    parameter(ten3m=1.0e-3_dp,explim=150.0_dp)

    if (BB%fk == zero)  return

    CALL ALLOC(xr,yr,r,r2,bbpar,cbx,cby,DZ)
    CALL ALLOC(sx2,sy2,xs,ys,rho2,tk,phix,phiy)
    CALL ALLOC(rk,xb,yb,crx,cry)


    SX=BB%SX
    SY=BB%SY
    XM=BB%XM
    YM=BB%YM
    FK=BB%FK

    sx2 = sx*sx
    sy2 = sy*sy
    !---- limit formulae for sigma(x) = sigma(y).
    if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then !ten3m = 1.0d-3
       !       do itrack = 1, b%n
       !          IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       rho2 = xs * xs + ys * ys
       tk = rho2 / (two * sx2)
       if (tk .gt. explim) then
          phix = xs * fk / rho2
          phiy = ys * fk / rho2
       else if (rho2 .ne. zero) then
          phix = xs * fk / rho2 * (one - exp(-tk) )
          phiy = ys * fk / rho2 * (one - exp(-tk) )
       else
          phix = zero
          phiy = zero
       endif
       !          if (.NOT.bborbit)  then
       !--- subtract closed orbit kick
       phix = phix - bb%bbk(1) ! subtract horizontal bb kick
       phiy = phiy - bb%bbk(2) ! subtract vertical co
       !          endif
       x(2) = x(2) + phix
       x(4) = x(4) + phiy
       !       enddo

       !---- case sigma(x) > sigma(y).
    else if (sx2 > sy2) then
       r2 = two * (sx2 - sy2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       !       do itrack = 1, b%n
       !          IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       xr = (xs) / r   !abs
       yr = (ys) / r   !abs
       call ccperrf(xr, yr, crx, cry)
       tk = (xs * xs / sx2 + ys * ys / sy2) / two
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
       !          track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
       !          track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
       IF(XS>=ZERO) THEN
          x(2) = x(2) + phix * XS
       ELSE
          x(2) = x(2) - phix * XS
       ENDIF
       IF(YS>=ZERO) THEN
          x(4) = x(4) + phiy * YS
       ELSE
          x(4) = x(4) - phiy * YS
       ENDIF
       !          if (.NOT.bborbit)  then
       !--- subtract closed orbit kick
       !            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
       !            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
       x(2) = x(2) - BB%bbk(1)
       x(4) = x(4) - BB%bbk(2)
       !          endif
       !       enddo

       !---- case sigma(x) < sigma(y).
    else
       r2 = two * (sy2 - sx2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       !       do itrack = 1, b%n
       !         IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       xr = (xs) / r !abs
       yr = (ys) / r !abs
       call ccperrf(yr, xr, cry, crx)
       tk = (xs * xs / sx2 + ys * ys / sy2) / two
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
       IF(XS>=ZERO) THEN
          x(2) = x(2) + phix * XS
       ELSE
          x(2) = x(2) - phix * XS
       ENDIF
       IF(YS>=ZERO) THEN
          x(4) = x(4) + phiy * YS
       ELSE
          x(4) = x(4) - phiy * YS
       ENDIF
       !          if (.NOT.bborbit)  then
       !--- subtract closed orbit kick
       !            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
       !            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
       x(2) = x(2) - BB%bbk(1)
       x(4) = x(4) - BB%bbk(2)
       !          endif
       !       enddo
    endif
    !    IF(DZ/=ZERO) CALL DRIFT_BEAM_BACK_TO_POSITION(th,DZ,B)
    CALL KILL(xr,yr,r,r2,bbpar,cbx,cby,DZ)
    CALL KILL(sx2,sy2,xs,ys,rho2,tk,phix,phiy)
    CALL KILL(rk,xb,yb,crx,cry)

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
    integer n,nc,nu
    TYPE(REAL_8) xx,yy,wx,wy,x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux,&
         &rx(33),ry(33)
    REAL(DP) xlim,ylim,fac1,fac2,fac3,cc,qr
    parameter(cc=1.12837916709551_dp,        &
         xlim=5.33_dp,ylim=4.29_dp,fac1=3.2_dp,fac2=23.0_dp,fac3=21.0_dp)
    CALL ALLOC(xx,yy,wx,wy,x,y,q,h,xl)
    CALL ALLOC(xh,yh,tx,ty,tn,sx,sy,saux)
    CALL ALLOC(RX,33)
    CALL ALLOC(RY,33)
    !    x = abs(xx)
    !    y = abs(yy)

    IF( X>=ZERO) THEN
       X=XX
    ELSE
       X=-XX
    ENDIF
    IF( Y>=ZERO) THEN
       Y=YY
    ELSE
       Y=-YY
    ENDIF
    if (y .lt. ylim  .and.  x .lt. xlim) then
       q  = (one - y / ylim) * sqrt(one - (x/xlim)**2)
       h  = one / (fac1 * q)
       qr=q
       nc = 7 + int(fac2*qr)
       xl = h**(1 - nc)
       xh = y + half/h
       yh = x
       nu = 10 + int(fac3*qr)
       rx(nu+1) = zero
       ry(nu+1) = zero

       do n = nu, 1, -1
          tx = xh + n * rx(n+1)
          ty = yh - n * ry(n+1)
          tn = tx*tx + ty*ty
          rx(n) = half * tx / tn
          ry(n) = half * ty / tn
       enddo

       sx = zero
       sy = zero

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
       rx(1) = zero
       ry(1) = zero

       do n = 9, 1, -1
          tx = xh + n * rx(1)
          ty = yh - n * ry(1)
          tn = tx*tx + ty*ty
          rx(1) = half * tx / tn
          ry(1) = half * ty / tn
       enddo

       wx = cc * rx(1)
       wy = cc * ry(1)
    endif

    !      if(y .eq. zero) wx = exp(-x**2)
    if(yy .lt. zero) then
       wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx
       wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy
       if(xx .gt. zero) wy = -wy
    else
       if(xx .lt. zero) wy = -wy
    endif

    CALL KILL(xx,yy,wx,wy,x,y,q,h,xl)
    CALL KILL(xh,yh,tx,ty,tn,sx,sy,saux)
    CALL KILL(RX,33)
    CALL KILL(RY,33)

  end SUBROUTINE ccperrfP

!!!!!!!!!!!!!!!  track single node !!!!!!!!!!!!!!!

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_R(R,X,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C


    CALL RESET_APERTURE_FLAG

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       if(associated(c%bb)) call BBKICK(c%BB,X)

       CALL TRACK_NODE_SINGLE(C,X,K,R%CHARGE)

       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_R

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_P(R,X,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C


    CALL RESET_APERTURE_FLAG

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       if(associated(c%bb)) call BBKICK(c%BB,X)

       CALL TRACK_NODE_SINGLE(C,X,K,R%CHARGE)

       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_P

  SUBROUTINE TRACK_FLAG_Beam_bb(R,b,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP) X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22,i
    TYPE (INTEGRATION_NODE), POINTER :: C


    CALL RESET_APERTURE_FLAG

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       DO I=1,B%N
          IF(B%U(i)) CYCLE
          X=BEAM_IN_X(B,I)!
          if(associated(c%bb)) call BBKICK(c%BB,X)

          CALL TRACK_NODE_SINGLE(C,X,K,R%CHARGE)
          CALL X_IN_BEAM(B,X,I,DL=ZERO,T=c%NEXT)

       ENDDO



       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_FLAG_Beam_bb



end module  beam_beam_ptc
