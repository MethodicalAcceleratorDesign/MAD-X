subroutine ttbb(track,ktrack,parvec)

  implicit none

  !----------------------------------------------------------------------*
  ! purpose:                                                             *
  !   track a set of particle through a beam-beam interaction region.    *
  !   see mad physicist's manual for the formulas used.                  *
  !input:                                                                *
  ! input/output:                                                        *
  !   track(6,*)(double)  track coordinates: (x, px, y, py, t, pt).      *
  !   ktrack    (integer) number of tracks.                              *
  !----------------------------------------------------------------------*
  integer ktrack,beamshape,b_dir_int
  double precision track(6,*),parvec(*),fk,dp
  double precision gamma0,beta0,beta_dp,ptot,b_dir
  double precision q,q_prime,get_value,node_value,get_variable
  double precision zero,one,two
  logical first
  save first
  data first / .true. /
  parameter(zero=0.0d0,one=1.0d0,two=2.0d0)
  !     if x > explim, exp(-x) is outside machine limits.

  !---- Calculate momentum deviation and according changes
  !     of the relativistic factor beta0
  dp  = get_variable('track_deltap ')
  q = get_value('probe ','charge ')
  q_prime = node_value('charge ')
  gamma0 = parvec(7)
  beta0 = sqrt(one-one/gamma0**2)
  ptot = beta0*gamma0*(one+dp)
  beta_dp = ptot / sqrt(one + ptot**2)
  b_dir_int = node_value('bbdir ')
  b_dir=dble(b_dir_int)
  b_dir = b_dir/sqrt(b_dir*b_dir + 1.0d-32)
  !---- pre-factor, if zero, anything else does not need to be calculated
  fk = two*parvec(5)*parvec(6)/parvec(7)/beta0/(one+dp)/q*          &
       (one-beta0*beta_dp*b_dir)/(beta_dp+0.5*(b_dir-one)*b_dir*beta0)
  !
  if (fk .eq. zero)  return
  !---- choose beamshape: 1-Gaussian, 2-flattop=trapezoidal, 3-hollow-parabolic
  beamshape = node_value('bbshape ')
  if(beamshape.lt.1.or.beamshape.gt.3) then
     beamshape=1
     if(first) then
        first = .false.
        call aawarn('TTBB: ',                                         &
             'beamshape out of range, set to default=1')
     endif
  endif
  if(beamshape.eq.1) call ttbb_gauss(track,ktrack,fk)
  if(beamshape.eq.2) call ttbb_flattop(track,ktrack,fk)
  if(beamshape.eq.3) call ttbb_hollowparabolic(track,ktrack,fk)

end subroutine ttbb
subroutine ttbb_gauss(track,ktrack,fk)

  use bbfi
  implicit none

  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having a Gaussian perpendicular shape                       *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  logical bborbit
  integer ktrack,itrack,ipos,get_option
  double precision track(6,*),pi,sx,sy,xm,ym,sx2,sy2,xs,            &
       ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,cbx,cby,      &
       get_variable,node_value,zero,one,two,three,ten3m,explim
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
       explim=150d0)
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  pi=get_variable('pi ')
  sx = node_value('sigx ')
  sy = node_value('sigy ')
  xm = node_value('xma ')
  ym = node_value('yma ')
  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif
  sx2 = sx*sx
  sy2 = sy*sy
  !---- limit formulae for sigma(x) = sigma(y).
  if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
     do itrack = 1, ktrack
        xs = track(1,itrack) - xm
        ys = track(3,itrack) - ym
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
        if (ipos .ne. 0)  then
           !--- subtract closed orbit kick
           phix = phix - bb_kick(1,ipos)
           phiy = phiy - bb_kick(2,ipos)
        endif
        track(2,itrack) = track(2,itrack) + phix
        track(4,itrack) = track(4,itrack) + phiy
     enddo

     !---- case sigma(x) > sigma(y).
  else if (sx2 .gt. sy2) then
     r2 = two * (sx2 - sy2)
     r  = sqrt(r2)
     rk = fk * sqrt(pi) / r
     do itrack = 1, ktrack
        xs = track(1,itrack) - xm
        ys = track(3,itrack) - ym
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
        track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
        track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
        if (ipos .ne. 0)  then
           !--- subtract closed orbit kick
           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
        endif
     enddo

     !---- case sigma(x) < sigma(y).
  else
     r2 = two * (sy2 - sx2)
     r  = sqrt(r2)
     rk = fk * sqrt(pi) / r
     do itrack = 1, ktrack
        xs = track(1,itrack) - xm
        ys = track(3,itrack) - ym
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
        track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
        track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
        if (ipos .ne. 0)  then
           !--- subtract closed orbit kick
           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
        endif
     enddo
  endif
end subroutine ttbb_gauss
subroutine ttbb_flattop(track,ktrack,fk)

  use bbfi
  implicit none

  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having an trapezoidal and, so flat top radial profile       *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  logical bborbit,first
  integer ktrack,itrack,ipos,get_option
  double precision track(6,*),pi,r0x,r0y,wi,wx,wy,xm,ym,            &
       r0x2,r0y2,xs,ys,rho,rho2,fk,phir,phix,phiy,get_variable,          &
       node_value,zero,one,two,three,ten3m,explim,norm,r1,zz
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
       explim=150d0)
  save first
  data first / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  pi=get_variable('pi ')
  ! mean radii of the is given via variables sigx and sigy
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  wi = node_value('width ')
  xm = node_value('xma ')
  ym = node_value('yma ')
  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif
  r0x2 = r0x*r0x
  r0y2 = r0y*r0y
  wx = r0x*wi
  wy = r0y*wi
  !---- limit formulae for mean radius(x) = mean radius(y),
  !-----      preliminary the only case considered.
  !
  if (abs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
     zz = 0.5*(r0x + r0y)
     r0x=zz
     r0y=zz
     r0x2=r0x*r0x
     r0y2=r0y*r0y
     if(first) then
        first=.false.
        call aawarn('TTBB_FLATTOP: ','beam is assumed to be circular')
     endif
  endif
  norm = (12.0*r0x**2 + wx**2)/24.0
  r1=r0x-wx/2.0
  do itrack = 1, ktrack
     xs = track(1,itrack) - xm
     ys = track(3,itrack) - ym
     rho2 = xs * xs + ys * ys
     rho  = sqrt(rho2)
     if(rho.le.r1) then
        phir = 0.5/norm
        phix = phir*xs
        phiy = phir*ys
     else if(rho.gt.r1.and.rho.lt.r1+wx) then
        phir = ((r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0 +            &
             wx**2/48.0)/rho2 + 0.25 + 0.5*r0x/wx -                            &
             rho/3.0/wx)/norm
        phix = phir*xs
        phiy = phir*ys
     else if(rho.ge.r1+wx) then
        phir = 1.0/rho2
        phix = xs*phir
        phiy = ys*phir
     endif
     track(2,itrack) = track(2,itrack)+phix*fk
     track(4,itrack) = track(4,itrack)+phiy*fk
  end do

end subroutine ttbb_flattop
subroutine ttbb_hollowparabolic(track,ktrack,fk)

  use bbfi
  implicit none

  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having a hollow-parabolic perpendicular shape               *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  logical bborbit,first
  integer ktrack,itrack,ipos,get_option
  double precision track(6,*),pi,r0x,r0y,wi,wx,wy,xm,ym,            &
       r0x2,r0y2,xs,ys,rho,rho2,fk,phir,phix,phiy,get_variable,          &
       node_value,zero,one,two,three,ten3m,explim,zz
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
       explim=150d0)
  save first
  data first / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  pi=get_variable('pi ')
  ! mean radii of the is given via variables sigx and sigy
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  wi = node_value('width ')
  ! width is given as FWHM of parabolic density profile, but formulas were
  ! derived with half width at the bottom of the parabolic density profile
  wi = wi/sqrt(2.0)
  xm = node_value('xma ')
  ym = node_value('yma ')
  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif
  r0x2 = r0x*r0x
  r0y2 = r0y*r0y
  wx  = wi*r0x
  wy  = wi*r0y
  !---- limit formulae for mean radius(x) = mean radius(y),
  !-----      preliminary the only case considered.
  !
  if (abs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
     zz = 0.5*(r0x + r0y)
     r0x=zz
     r0y=zz
     r0x2=r0x*r0x
     r0y2=r0y*r0y
     if(first) then
        first=.false.
        call aawarn('TTBB_HOLLOWPARABOLIC: ',                         &
             'beam is assumed to be circular')
     endif
  endif
  do itrack = 1, ktrack
     xs = track(1,itrack) - xm
     ys = track(3,itrack) - ym
     rho2 = xs * xs + ys * ys
     rho  = sqrt(rho2)
     if(rho.le.r0x-wx) then
        phix = zero
        phiy = zero
     else if(rho.gt.r0x-wx.and.rho.lt.r0x+wx) then
        phir=0.75/wx/r0x/rho2*(r0x**4/12.0/wx**2 - r0x**2/2.0 +       &
             2.0*r0x*wx/3.0 - wx**2/4.0 + rho2/2.0*(1.0 -                      &
             r0x**2/wx**2) + rho**3/3.0*2.0*r0x/wx**2 -                        &
             rho**4/4.0/wx**2)
        phix = phir*xs
        phiy = phir*ys
     else
        phir = 1.0/rho2
        phix = xs*phir
        phiy = ys*phir
     endif
     track(2,itrack) = track(2,itrack)+phix*fk
     track(4,itrack) = track(4,itrack)+phiy*fk
  end do

end subroutine ttbb_hollowparabolic
!!$subroutine ttbb_old(track,ktrack,parvec)
!!$
!!$  use bbfi
!!$  implicit none
!!$
!!$  !----------------------------------------------------------------------*
!!$  ! purpose:                                                             *
!!$  !   track a set of particle through a beam-beam interaction region.    *
!!$  !   see mad physicist's manual for the formulas used.                  *
!!$  !input:                                                                *
!!$  ! input/output:                                                        *
!!$  !   track(6,*)(double)  track coordinates: (x, px, y, py, t, pt).      *
!!$  !   ktrack    (integer) number of tracks.                              *
!!$  !----------------------------------------------------------------------*
!!$  logical bborbit
!!$  integer ktrack,itrack,ipos,get_option
!!$  double precision track(6,*),parvec(*),pi,sx,sy,xm,ym,sx2,sy2,xs,  &
!!$       ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,cbx,cby,      &
!!$       get_variable,node_value,zero,one,two,three,ten3m,explim
!!$  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
!!$       explim=150d0)
!!$  !     if x > explim, exp(-x) is outside machine limits.
!!$
!!$  !---- initialize.
!!$  bborbit = get_option('bborbit ') .ne. 0
!!$  pi=get_variable('pi ')
!!$  sx = node_value('sigx ')
!!$  sy = node_value('sigy ')
!!$  xm = node_value('xma ')
!!$  ym = node_value('yma ')
!!$  fk = two * parvec(5) * parvec(6) / parvec(7)
!!$  if (fk .eq. zero)  return
!!$  ipos = 0
!!$  if (.not. bborbit)  then
!!$     !--- find position of closed orbit bb_kick
!!$     do ipos = 1, bbd_cnt
!!$        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
!!$     enddo
!!$     ipos = 0
!!$1    continue
!!$  endif
!!$  sx2 = sx*sx
!!$  sy2 = sy*sy
!!$  !---- limit formulae for sigma(x) = sigma(y).
!!$  if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
!!$     do itrack = 1, ktrack
!!$        xs = track(1,itrack) - xm
!!$        ys = track(3,itrack) - ym
!!$        rho2 = xs * xs + ys * ys
!!$        tk = rho2 / (two * sx2)
!!$        if (tk .gt. explim) then
!!$           phix = xs * fk / rho2
!!$           phiy = ys * fk / rho2
!!$        else if (rho2 .ne. zero) then
!!$           phix = xs * fk / rho2 * (one - exp(-tk) )
!!$           phiy = ys * fk / rho2 * (one - exp(-tk) )
!!$        else
!!$           phix = zero
!!$           phiy = zero
!!$        endif
!!$        if (ipos .ne. 0)  then
!!$           !--- subtract closed orbit kick
!!$           phix = phix - bb_kick(1,ipos)
!!$           phiy = phiy - bb_kick(2,ipos)
!!$        endif
!!$        track(2,itrack) = track(2,itrack) + phix
!!$        track(4,itrack) = track(4,itrack) + phiy
!!$     enddo
!!$
!!$     !---- case sigma(x) > sigma(y).
!!$  else if (sx2 .gt. sy2) then
!!$     r2 = two * (sx2 - sy2)
!!$     r  = sqrt(r2)
!!$     rk = fk * sqrt(pi) / r
!!$     do itrack = 1, ktrack
!!$        xs = track(1,itrack) - xm
!!$        ys = track(3,itrack) - ym
!!$        xr = abs(xs) / r
!!$        yr = abs(ys) / r
!!$        call ccperrf(xr, yr, crx, cry)
!!$        tk = (xs * xs / sx2 + ys * ys / sy2) / two
!!$        if (tk .gt. explim) then
!!$           phix = rk * cry
!!$           phiy = rk * crx
!!$        else
!!$           xb = (sy / sx) * xr
!!$           yb = (sx / sy) * yr
!!$           call ccperrf(xb, yb, cbx, cby)
!!$           phix = rk * (cry - exp(-tk) * cby)
!!$           phiy = rk * (crx - exp(-tk) * cbx)
!!$        endif
!!$        track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
!!$        track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
!!$        if (ipos .ne. 0)  then
!!$           !--- subtract closed orbit kick
!!$           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
!!$           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
!!$        endif
!!$     enddo
!!$
!!$     !---- case sigma(x) < sigma(y).
!!$  else
!!$     r2 = two * (sy2 - sx2)
!!$     r  = sqrt(r2)
!!$     rk = fk * sqrt(pi) / r
!!$     do itrack = 1, ktrack
!!$        xs = track(1,itrack) - xm
!!$        ys = track(3,itrack) - ym
!!$        xr = abs(xs) / r
!!$        yr = abs(ys) / r
!!$        call ccperrf(yr, xr, cry, crx)
!!$        tk = (xs * xs / sx2 + ys * ys / sy2) / two
!!$        if (tk .gt. explim) then
!!$           phix = rk * cry
!!$           phiy = rk * crx
!!$        else
!!$           xb  = (sy / sx) * xr
!!$           yb  = (sx / sy) * yr
!!$           call ccperrf(yb, xb, cby, cbx)
!!$           phix = rk * (cry - exp(-tk) * cby)
!!$           phiy = rk * (crx - exp(-tk) * cbx)
!!$        endif
!!$        track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
!!$        track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
!!$        if (ipos .ne. 0)  then
!!$           !--- subtract closed orbit kick
!!$           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
!!$           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
!!$        endif
!!$     enddo
!!$  endif
!!$end subroutine ttbb_old
