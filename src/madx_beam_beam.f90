module BeamBeam

contains

  SUBROUTINE tmbb(fsec,ftrk,orbit,fmap,re,te)
    use trackfi, only : fsecarb
    use twissbeamfi, only : gamma, arad, charge, npart
    use math_constfi, only : zero, one, two
    implicit none
    !----------------------------------------------------------------------*
    !     purpose:                                                         *
    !     transport map for beam-beam element.                             *
    !     input:                                                           *
    !     fsec      (logical) must be .true. for this purpose              *
    !     ftrk      (logical) must be true for this purpose                *
    !     input/output:                                                    *
    !     orbit(6)  (double)  closed orbit (only kick is added             *
    !     since BB thin)                                                   *
    !     output:                                                          *
    !     fmap      (logical) if map has been calculated correctly         *
    !     re(6,6)   (double)  transfer matrix.                             *
    !     te(6,6,6) (double)  second-order terms.                          *
    !----------------------------------------------------------------------*
    logical, intent(IN)  :: fsec, ftrk
    logical, intent(OUT) :: fmap
    double precision, intent(IN OUT) :: orbit(6)
    double precision, intent(OUT) :: re(6,6), te(6,6,6)
    integer :: beamshape, b_dir_int
    logical, save :: first=.true.
    logical :: bb_ultra_relati
    double precision :: parvec(26), fk, q, q_prime, dp
    double precision :: gamma0, beta0, beta_dp, ptot, b_dir
    integer, external :: get_option
    double precision, external :: node_value, get_variable
    !---  standard 4D
    q = charge
    q_prime = node_value('charge ')
    parvec(5) = arad
    parvec(6) = q_prime * npart
    parvec(7) = gamma
    !---- Calculate momentum deviation and according changes
    !     of the relativistic factor beta0
    dp = get_variable('track_deltap ')
    gamma0 = parvec(7)
    beta0 = sqrt(one-one/gamma0**2)
    ptot = beta0*gamma0*(one+dp)
    beta_dp = ptot / sqrt(one + ptot**2)
    b_dir_int = node_value('bbdir ')
    b_dir = dble(b_dir_int)
    b_dir = b_dir/sqrt(b_dir*b_dir + 1.0d-32)
    bb_ultra_relati=get_option('bb_ultra_relati ').ne.0
    if (bb_ultra_relati) then
       fk = two * parvec(5) * parvec(6) / parvec(7)
    else
       fk = two*parvec(5)*parvec(6)/parvec(7)/beta0/(one+dp)/q*          &
            (one-beta0*beta_dp*b_dir)/(beta_dp+0.5*(b_dir-one)*b_dir*beta0)
    endif
    !---- choose beamshape: 1-Gaussian (default), 2-flattop=trapezoidal, 3-hollow-parabolic
    beamshape = node_value('bbshape ')
    select case (beamshape)
    case (1)
       call tmbb_gauss(fsec,ftrk,orbit,fmap,re,te,fk)
    case (2)
       call tmbb_flattop(fsec,ftrk,orbit,fmap,re,te,fk)
    case (3)
       call tmbb_hollowparabolic(fsec,ftrk,orbit,fmap,re,te,fk)
    case default
       if (first) then
          first = .false.
          call fort_warn('TMBB: ', 'beamshape out of range, set to default=1')
       endif
       beamshape=1
       call tmbb_gauss(fsec,ftrk,orbit,fmap,re,te,fk)
    end select
  end SUBROUTINE tmbb
  SUBROUTINE tmbb_gauss(fsec,ftrk,orbit,fmap,re,te,fk)
    use bbfi
    use twisslfi
    use spch_bbfi
    use fasterror
    use math_constfi, only : zero, one, two, three, ten3m, pi
    implicit none
    logical :: fsec, ftrk, fmap
    double precision :: orbit(6), re(6,6), te(6,6,6)
    logical :: bborbit, bb_sxy_update, long_coup_on
    integer :: mylen
    double precision :: sx,sy, xm, ym, sx2, sy2, xs, ys, rho2, fk, tk, exk
    double precision :: phix, phiy, rho4, phixx, phixy, phiyy, rho6, rk, exkc
    double precision :: xb, yb, phixxx, phixxy, phixyy, phiyyy, crx, cry, xr
    double precision :: yr, r, r2, cbx, cby
    character(len=20) :: text
    character(len=name_len) :: name
    integer, external ::  get_option
    double precision, external :: node_value
    double precision, external :: get_value
    !---- initialize.
    bborbit = get_option('bborbit ') .ne. 0
    if (bbd_flag.ne.0 .and. .not.bborbit)  then
       if (bbd_cnt .eq. bbd_max)  then
          call fort_warn('TMBB_GAUSS: ','maximum bb number reached')
       else
          bbd_cnt = bbd_cnt + 1
          bbd_loc(bbd_cnt) = bbd_pos
          bb_kick(1,bbd_cnt) = zero
          bb_kick(2,bbd_cnt) = zero
       endif
    endif
    fmap = .true.
    bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
    long_coup_on  = get_option('long_coup_off ') .eq. 0
    !frs on 06.06.2016
    !  safeguard TWISS from failing due to undefined SC elements
    if (bb_sxy_update .and. long_coup_on .and. N_spch .gt. 0) then
       fk = fk * rat_bb_n_ions !Ratio_for_bb_N_ions
       name=' '
       call element_name(name,len(name))
       mylen = len_trim(name)
       i_spch = i_spch + 1
       if (i_spch .gt. N_spch) then
          write(text, '(1p,i8)') i_spch
          call fort_fail('TMBB: ', 'Table with too few BB elements: '//text)
       endif
       if (spch_bb_name(i_spch)(:mylen) .ne. name(:mylen)) &
            call fort_fail('TMBB: ', 'wrong element name in Table: spch_bb')
       sx = sqrt(betx_bb(i_spch)*Ex_rms + (dx_bb(i_spch)*sigma_p)**2)
       sy = sqrt(bety_bb(i_spch)*Ey_rms + (dy_bb(i_spch)*sigma_p)**2)
    else
       sx = node_value('sigx ')
       sy = node_value('sigy ')
    endif
    if (sx < 1e-16 .or. sy < 1e-16) then
       re = zero;
       re(1,1) = 1;
       re(2,2) = 1;
       re(3,3) = 1;
       re(4,4) = 1;
       return;
    endif
    xm = node_value('xma ')
    ym = node_value('yma ')
    if (fk .ne. zero)  then
       !---- if tracking is desired ...
       if (ftrk) then
          sx2 = sx * sx
          sy2 = sy * sy
          xs  = orbit(1) - xm
          ys  = orbit(3) - ym
          !---- limit formulas for sigma(x) = sigma(y).
          if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
             rho2 = xs * xs + ys * ys
             !---- limit case for xs = ys = 0.
             if (rho2 .eq. zero) then
                re(2,1) = fk / (two * sx2)
                re(4,3) = fk / (two * sx2)
                !---- general case.
             else
                tk = rho2 / (two * sx2)
                if (tk .gt. explim) then
                   ! if x > explim, exp(-x) is outside machine limits.
                   exk  = zero
                   exkc = one
                   phix = xs * fk / rho2
                   phiy = ys * fk / rho2
                else
                   exk  = exp(-tk)
                   exkc = one - exk
                   phix = xs * fk / rho2 * exkc
                   phiy = ys * fk / rho2 * exkc
                endif
                !---- orbit kick - only applied if option bborbit (HG 5/12/01),
                !     else stored
                if (bborbit) then
                   orbit(2) = orbit(2) + phix
                   orbit(4) = orbit(4) + phiy
                elseif (bbd_flag .ne. 0)  then
                   bb_kick(1,bbd_cnt) = phix
                   bb_kick(2,bbd_cnt) = phiy
                endif
                !---- first-order effects.
                rho4 = rho2 * rho2
                phixx = fk * (- exkc * (xs*xs - ys*ys) / rho4 + exk * xs*xs / (rho2 * sx2))
                phixy = fk * (- exkc * two * xs * ys   / rho4 + exk * xs*ys / (rho2 * sx2))
                phiyy = fk * (+ exkc * (xs*xs - ys*ys) / rho4 + exk * ys*ys / (rho2 * sx2))
                re(2,1) = phixx
                re(2,3) = phixy
                re(4,1) = phixy
                re(4,3) = phiyy
                !---- second-order effects.
                if (fsec) then
                   rho6 = rho4 * rho2
                   phixxx = fk*xs * (+ exkc * (xs*xs - three*ys*ys) / rho6 &
                        - exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2) &
                        - exk * xs*xs / (two * rho2 * sx2**2))
                   phixxy = fk*ys * (+ exkc * (three*xs*xs - ys*ys) / rho6 &
                        - exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2) &
                        - exk * xs*xs / (two * rho2 * sx2**2))
                   phixyy = fk*xs * (- exkc * (xs*xs - three*ys*ys) / rho6 &
                        + exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2) &
                        - exk * ys*ys / (two * rho2 * sx2**2))
                   phiyyy = fk*ys * (- exkc * (three*xs*xs - ys*ys) / rho6 &
                        + exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2) &
                        - exk * ys*ys / (two * rho2 * sx2**2))
                   te(2,1,1) = phixxx
                   te(2,1,3) = phixxy
                   te(2,3,1) = phixxy
                   te(4,1,1) = phixxy
                   te(2,3,3) = phixyy
                   te(4,1,3) = phixyy
                   te(4,3,1) = phixyy
                   te(4,3,3) = phiyyy
                endif
             endif
             !---- case sigma(x) > sigma(y).
          else
             r2 = two * (sx2 - sy2)
             if (sx2 .gt. sy2) then
                r  = sqrt(r2)
                rk = fk * sqrt(pi) / r
                xr = abs(xs) / r
                yr = abs(ys) / r
                if (fasterror_on) then
                   call wzsub(xr, yr, crx, cry)
                else
                   call ccperrf(xr, yr, crx, cry)
                endif
                tk = (xs * xs / sx2 + ys * ys / sy2) / two
                if (tk .gt. explim) then
                   ! if x > explim, exp(-x) is outside machine limits.
                   exk = zero
                   cbx = zero
                   cby = zero
                else
                   exk = exp(-tk)
                   xb  = (sy / sx) * xr
                   yb  = (sx / sy) * yr
                   if (fasterror_on) then
                      call wzsub(xb, yb, cbx, cby)
                   else
                      call ccperrf(xb, yb, cbx, cby)
                   endif
                endif
                !---- case sigma(x) < sigma(y).
             else
                r  = sqrt(-r2)
                rk = fk * sqrt(pi) / r
                xr = abs(xs) / r
                yr = abs(ys) / r
                if (fasterror_on) then
                   call wzsub(yr, xr, cry, crx)
                else
                   call ccperrf(yr, xr, cry, crx)
                endif
                tk = (xs * xs / sx2 + ys * ys / sy2) / two
                if (tk .gt. explim) then
                   ! if x > explim, exp(-x) is outside machine limits.
                   exk = zero
                   cbx = zero
                   cby = zero
                else
                   exk = exp(-tk)
                   xb  = (sy / sx) * xr
                   yb  = (sx / sy) * yr
                   if (fasterror_on) then
                      call wzsub(yb, xb, cby, cbx)
                   else
                      call ccperrf(yb, xb, cby, cbx)
                   endif
                endif
             endif
             phix = rk * (cry - exk * cby) * sign(one, xs)
             phiy = rk * (crx - exk * cbx) * sign(one, ys)
             !---- orbit kick - only applied if option bborbit (HG 5/12/01),
             !     else stored
             if (bborbit) then
                orbit(2) = orbit(2) + phix
                orbit(4) = orbit(4) + phiy
             elseif (bbd_flag .ne. 0)  then
                bb_kick(1,bbd_cnt) = phix
                bb_kick(2,bbd_cnt) = phiy
             endif
             !---- first-order effects.
             phixx = (two / r2) * (- (xs * phix + ys * phiy) + fk * (one - (sy / sx) * exk))
             phixy = (two / r2) * (- (xs * phiy - ys * phix))
             phiyy = (two / r2) * (+ (xs * phix + ys * phiy) - fk * (one - (sx / sy) * exk))
             re(2,1) = phixx
             re(2,3) = phixy
             re(4,1) = phixy
             re(4,3) = phiyy
             !---- second-order effects.
             if (fsec) then
                phixxx = (- phix - (xs * phixx + ys * phixy) + fk * xs * sy * exk / sx**3) / r2
                phixxy = (- phiy - (xs * phixy - ys * phixx)) / r2
                phixyy = (+ phix - (xs * phiyy - ys * phixy)) / r2
                phiyyy = (+ phiy + (xs * phixy + ys * phiyy) - fk * ys * sx * exk / sy**3) / r2
                te(2,1,1) = phixxx
                te(2,1,3) = phixxy
                te(2,3,1) = phixxy
                te(4,1,1) = phixxy
                te(2,3,3) = phixyy
                te(4,1,3) = phixyy
                te(4,3,1) = phixyy
                te(4,3,3) = phiyyy
             endif
          endif
          !---- no tracking desired.
       else
          re(2,1) = fk / (sx * (sx + sy))
          re(4,3) = fk / (sy * (sx + sy))
       endif
    endif
    !  print *, 'bborbit=', bborbit, 'bbd_flag=', bbd_flag, ', fk=', fk, &
    !    ', bbd_pos=', bbd_pos, ', bbd_cnt=', bbd_cnt, &
    !    ', bb_kick_x=', bb_kick(1,bbd_cnt), ', bb_kick_y=', bb_kick(2,bbd_cnt)
  end SUBROUTINE tmbb_gauss
  SUBROUTINE tmbb_flattop(fsec,ftrk,orbit,fmap,re,te,fk)
    use bbfi
    use twisslfi
    use math_constfi, only : zero, one, two, three, ten3m, pi
    implicit none
    logical :: fsec, ftrk, fmap
    double precision :: fk
    double precision :: orbit(6), re(6,6), te(6,6,6)
    logical :: bborbit
    logical, save :: firstflag=.true.
    double precision :: r0x, r0y, xm, ym, r0x2, r0y2, xs, ys, rho2
    double precision :: phix, phiy, phixx, phixy, phiyy, phixxx, phixxy, phixyy, phiyyy
    double precision :: wi, rho, wx, wy, norm, phir, phirr, phirrr, zz
    integer, external ::  get_option
    double precision, external :: node_value
    !---- initialize.
    phix = zero
    phiy = zero
    bborbit = get_option('bborbit ') .ne. 0
    if (bbd_flag.ne.0 .and. .not.bborbit)  then
       if (bbd_cnt .eq. bbd_max)  then
          call fort_warn('TMBB_FLATTOP: ','maximum bb number reached')
       else
          bbd_cnt = bbd_cnt + 1
          bbd_loc(bbd_cnt) = bbd_pos
          bb_kick(1,bbd_cnt) = zero
          bb_kick(2,bbd_cnt) = zero
       endif
    endif
    fmap = .true.
    r0x = node_value('sigx ')
    r0y = node_value('sigy ')
    if (r0x < 1e-16 .or. r0y < 1e-16) then
       re = zero;
       re(1,1) = 1;
       re(2,2) = 1;
       re(3,3) = 1;
       re(4,4) = 1;
       return;
    endif
    wi = node_value('width ')
    xm = node_value('xma ')
    ym = node_value('yma ')
    if (fk .ne. zero)  then
       r0x2 = r0x * r0x
       r0y2 = r0y * r0y
       wx  = r0x * wi
       wy  = r0y * wi
       xs  = orbit(1) - xm
       ys  = orbit(3) - ym
       !---- limit formulas for sigma(x) = sigma(y).
       !-----preliminary the only case considered.
       if (dabs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
          zz = 0.5*(r0x + r0y)
          r0x = zz
          r0y = zz
          r0x2 = r0x*r0x
          r0y2 = r0y*r0y
          if (firstflag) then
             firstflag=.false.
             call fort_warn('TMBB_FLATTOP: ', 'beam is assumed to be circular')
          endif
       endif
       norm = (12.0*r0x**2 + wx**2)/24.0
       !---- if tracking is desired ...
       if (ftrk) then
          rho2 = xs * xs + ys * ys
          rho  = sqrt(rho2)
          if (rho.le.r0x-wx/2.0) then
             phir = 0.5/norm
             phix = phir*xs
             phiy = phir*ys
             phixx = 0.5/norm
             phixy = zero
             phiyy = 0.5/norm
             re(2,1) = phixx*fk
             re(2,3) = phixy*fk
             re(4,1) = phixy*fk
             re(4,3) = phiyy*fk
             if (fsec) then
                te(2,1,1) = zero
                te(2,1,3) = zero
                te(2,3,1) = zero
                te(4,1,1) = zero
                te(2,3,3) = zero
                te(4,1,3) = zero
                te(4,3,1) = zero
                te(4,3,3) = zero
             endif
          else if (rho.gt.r0x-wx/2.0 .and. rho.lt.r0x+wx/2.0) then
             phir = ((r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0 + &
                  wx**2/48.0)/rho2 + 0.25 + 0.5*r0x/wx - rho/3.0/wx)/norm
             phix = phir*xs
             phiy = phir*ys
             phirr = - (2.0*(r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0 +   &
                  wx**2/48.0)/rho2**2 + 1.0/rho/3.0/wx)/norm
             phixx = phir + xs*xs*phirr
             phixy = xs*ys*phirr
             phiyy = phir + ys*ys*phirr
             re(2,1) = phixx*fk
             re(2,3) = phixy*fk
             re(4,1) = phixy*fk
             re(4,3) = phiyy*fk
             !     2nd order effects
             if (fsec) then
                phirrr =  (8.0*(r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0   &
                     + wx**2/48.0)/rho2**3 + 1.0/3.0/wx/rho**3)/norm
                phixxx = 3.0*xs*phirr + xs**3*phirrr
                phixxy = ys*phirr + xs*xs*ys*phirrr
                phixyy = xs*phirr + xs*ys*ys*phirrr
                phiyyy = 3.0*ys*phirr + ys**3*phirrr
                te(2,1,1) = phixxx*fk
                te(2,1,3) = phixxy*fk
                te(2,3,1) = phixxy*fk
                te(4,1,1) = phixxy*fk
                te(2,3,3) = phixyy*fk
                te(4,1,3) = phixyy*fk
                te(4,3,1) = phixyy*fk
                te(4,3,3) = phiyyy*fk
             endif
          else if (rho .ge. r0x+wx/2.0) then
             phir = 1.0/rho2
             phix = xs*phir
             phiy = ys*phir
             phirr = -2.0/rho2**2
             phixx = phir+xs*xs*phirr
             phixy = xs*ys*phirr
             phiyy = phir+ys*ys*phirr
             re(2,1) = phixx*fk
             re(2,3) = phixy*fk
             re(4,1) = phixy*fk
             re(4,3) = phiyy*fk
             !     2nd order effects
             if (fsec) then
                phirrr = 8.0/rho2**3
                phixxx = 3.0*xs*phirr + xs**3*phirrr
                phixxy = ys*phirr + xs*xs*ys*phirrr
                phixyy = xs*phirr + xs*ys*ys*phirrr
                phiyyy = 3.0*ys*phirr + ys**3*phirrr
                te(2,1,1) = phixxx*fk
                te(2,1,3) = phixxy*fk
                te(2,3,1) = phixxy*fk
                te(4,1,1) = phixxy*fk
                te(2,3,3) = phixyy*fk
                te(4,1,3) = phixyy*fk
                te(4,3,1) = phixyy*fk
                te(4,3,3) = phiyyy*fk
             end if
          endif
          !---- orbit kick - only applied if option bborbit (HG 5/12/01),
          !     else stored
          if (bborbit) then
             orbit(2) = orbit(2) + phix*fk
             orbit(4) = orbit(4) + phiy*fk
          elseif (bbd_flag .ne. 0)  then
             bb_kick(1,bbd_cnt) = phix*fk
             bb_kick(2,bbd_cnt) = phiy*fk
          endif
       else
          !---- no tracking desired.
          phixx=0.5/norm
          phiyy=0.5/norm
          re(2,1) = phixx*fk
          re(4,3) = phiyy*fk
       endif
    endif
  end SUBROUTINE tmbb_flattop
  SUBROUTINE tmbb_hollowparabolic(fsec,ftrk,orbit,fmap,re,te,fk)
    use bbfi
    use twisslfi
    use math_constfi, only : zero, one, two, three, ten3m, pi
    implicit none
    logical :: fsec, ftrk, fmap
    double precision :: fk
    double precision :: orbit(6), re(6,6), te(6,6,6)
    logical :: bborbit
    logical, save :: firstflag=.true.
    double precision :: r0x, r0y, xm, ym, r0x2, r0y2, xs, ys, rho2
    double precision :: phix, phiy, phixx, phixy, phiyy, phixxx, phixxy, phixyy, phiyyy
    double precision :: wi, rho, wx, wy, phir, phirr, phirrr, zz
    integer, external :: get_option
    double precision, external :: node_value
    !---- initialize.
    phix = zero
    phiy = zero
    bborbit = get_option('bborbit ') .ne. 0
    if (bbd_flag.ne.0 .and. .not.bborbit)  then
       if (bbd_cnt .eq. bbd_max)  then
          call fort_warn('TMBB_HOLLOWPARABOLIC: ', 'maximum bb number reached')
       else
          bbd_cnt = bbd_cnt + 1
          bbd_loc(bbd_cnt) = bbd_pos
          bb_kick(1,bbd_cnt) = zero
          bb_kick(2,bbd_cnt) = zero
       endif
    endif
    fmap = .true.
    r0x = node_value('sigx ')
    r0y = node_value('sigy ')
    if (r0x < 1e-16 .or. r0y < 1e-16) then
       re = zero;
       re(1,1) = 1;
       re(2,2) = 1;
       re(3,3) = 1;
       re(4,4) = 1;
       return;
    endif
    !     width is given as FWHM of the parabolic density profile,
    !     but formulas were derived with half width at the bottom of this
    !     density profile
    wi = node_value('width ')
    wi = wi/sqrt(2.0)
    xm = node_value('xma ')
    ym = node_value('yma ')
    if (fk .ne. zero)  then
       !---- if tracking is desired ...
       if (ftrk) then
          r0x2 = r0x * r0x
          r0y2 = r0y * r0y
          wx  = r0x * wi
          wy  = r0y * wi
          xs  = orbit(1) - xm
          ys  = orbit(3) - ym
          !---- limit formulas for sigma(x) = sigma(y).
          !-----preliminary the only case considered.
          if (abs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
             zz = 0.5*(r0x + r0y)
             r0x=zz
             r0y=zz
             r0x2=r0x*r0x
             r0y2=r0y*r0y
             if (firstflag) then
                firstflag=.false.
                call fort_warn('TMBB_HOLLOWPARABOLIC: ', 'beam is assumed to be circular')
             endif
          endif
          rho2 = xs * xs + ys * ys
          rho  = sqrt(rho2)
          if (rho.le.r0x-wx) then
             re(2,1) = zero
             re(4,3) = zero
             re(2,3) = zero
             re(4,1) = zero
             if (fsec) then
                te(2,1,1) = zero
                te(2,1,3) = zero
                te(2,3,1) = zero
                te(4,1,1) = zero
                te(2,3,3) = zero
                te(4,1,3) = zero
                te(4,3,1) = zero
                te(4,3,3) = zero
             end if
          else if (rho.gt.r0x-wx .and. rho.lt.r0x+wx) then
             phir=0.75/wx/r0x/rho2*(r0x**4/12.0/wx**2 - r0x**2/2.0 + &
                  2.0*r0x*wx/3.0 - wx**2/4.0 + rho2/2.0*(1.0 -       &
                  r0x**2/wx**2) + rho**3/3.0*2.0*r0x/wx**2 -         &
                  rho**4/4.0/wx**2)
             phix = phir*xs
             phiy = phir*ys
             phirr = 0.75/wx/r0x*(-2.0/rho**4*(r0x**4/12.0/wx**2 -   &
                  r0x**2/2.0+2.0*r0x*wx/3.0-wx**2/4.0) +             &
                  2.0*r0x/wx**2/3.0/rho-0.5/wx**2)
             phixx = phir+xs*xs*phirr
             phixy = xs*ys*phirr
             phiyy = phir+ys*ys*phirr
             re(2,1) = phixx*fk
             re(2,3) = phixy*fk
             re(4,1) = phixy*fk
             re(4,3) = phiyy*fk
             ! 2nd order effects
             if (fsec) then
                phirrr = 0.75/wx/r0x*(8.0/rho2**3*(r0x**4/12.0/           &
                     wx**2-r0x**2/2.0+2.0*r0x*wx/3.0-wx**2/4.0) -         &
                     2.0*r0x/3.0/wx**2/rho**3)
                phixxx = 3.0*xs*phirr + xs**3*phirrr
                phixxy = ys*phirr + xs*xs*ys*phirrr
                phixyy = xs*phirr + xs*ys*ys*phirrr
                phiyyy = 3.0*ys*phirr + ys**3*phirrr
                te(2,1,1) = phixxx*fk
                te(2,1,3) = phixxy*fk
                te(2,3,1) = phixxy*fk
                te(4,1,1) = phixxy*fk
                te(2,3,3) = phixyy*fk
                te(4,1,3) = phixyy*fk
                te(4,3,1) = phixyy*fk
                te(4,3,3) = phiyyy*fk
             end if
          else if (rho .ge. r0x+wx) then
             phir = 1.0/rho2
             phix = xs*phir
             phiy = ys*phir
             phirr = -2.0/rho2**2
             phixx = phir+xs*xs*phirr
             phixy = xs*ys*phirr
             phiyy = phir+ys*ys*phirr
             re(2,1) = phixx*fk
             re(2,3) = phixy*fk
             re(4,1) = phixy*fk
             re(4,3) = phiyy*fk
             ! 2nd order effects
             if (fsec) then
                phirrr = 8.0/rho2**3
                phixxx = 3.0*xs*phirr + xs**3*phirrr
                phixxy = ys*phirr + xs*xs*ys*phirrr
                phixyy = xs*phirr + xs*ys*ys*phirrr
                phiyyy = 3.0*ys*phirr + ys**3*phirrr
                te(2,1,1) = phixxx*fk
                te(2,1,3) = phixxy*fk
                te(2,3,1) = phixxy*fk
                te(4,1,1) = phixxy*fk
                te(2,3,3) = phixyy*fk
                te(4,1,3) = phixyy*fk
                te(4,3,1) = phixyy*fk
                te(4,3,3) = phiyyy*fk
             end if
          endif
          !---- orbit kick - only applied if option bborbit (HG 5/12/01),
          !     else stored
          if (bborbit) then
             orbit(2) = orbit(2) + phix*fk
             orbit(4) = orbit(4) + phiy*fk
          elseif (bbd_flag .ne. 0)  then
             bb_kick(1,bbd_cnt) = phix*fk
             bb_kick(2,bbd_cnt) = phiy*fk
          endif
       else
          !---- no tracking desired.
          re(2,1) = zero
          re(4,3) = zero
       endif
    endif
  end SUBROUTINE tmbb_hollowparabolic

end module BeamBeam
