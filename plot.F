      subroutine pearwe (ibk, actwin)
!----------------------------------------------------------------------*
! Purpose:
!   Plot Arnold web
!
!--- Input
!   IBK       (integer) array containing:
!                       number of superperiods
!                       number of constraints  N
!                       N constraints:
!                       minimum
!                       maximum
!                       step
!                       Length L
!                       expression in polish notation of length L,
!                       coded as 1+, 2-, 3*, 4/, 1 KX, 2 KY, 3 KS
!                       + MQADD + 4, e.g. 100006 = KY
!   ACTWIN       (real) active user window (WC)
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer    ibk(*)
      real       actwin(4)
      integer klwid, ipm, ipp, ip, j
      integer ifp, kx, ky, kax, kay, ksmax, k1, k2, k3, iselct, ierr
      integer k(3), l(2), lstyl(0:3)
      real s, qs, v1, v2, vm, vp, sclwid
      real ql(2), qu(2), x(2), y(2)
      data lstyl / 1, 2, 4, 3 /

!--- set flag for first error message (following suppressed)
      ifp = 0
      if (horname(:2) .eq. 'QX' .or. horname(:2) .eq. 'qx')  then
        kx = 1
      else
        kx = 2
      endif
      ky  = 3 - kx
      kax = 2 * kx -1
      kay = 2 * ky -1
      ql(kx) = actwin(kax)
      qu(kx) = actwin(kax+1)
      ql(ky) = actwin(kay)
      qu(ky) = actwin(kay+1)
      qs = qsval
      s   = ibk(1)
      if (qs .eq. 0.)  then
        ksmax = 0
      else
        ksmax = mksmax
      endif
      do k1 = -ntmax, ntmax
        k(kx) = k1
        l(kx) = ntmax - abs(k1)
        do k2 = -l(kx), l(kx)
          if (k1 .ne. 0 .or. k2 .ne. 0)  then
            k(ky) = k2
            l(ky) = min(l(kx) - abs(k2), ksmax)
            do k3 = -l(ky), l(ky)
              k(3) = k3
              call peqcon(kx, ky, k, ibk, iselct, ierr)
              if (iselct .ne. 0 .or. ierr .ne. 0)  then
                call jsln(lstyl(min(abs(k3),3)))
                if (ierr .ne. 0 .and. ifp .eq. 0)  then
                  ifp = 1
!                  call pewarn ('PLARWE', 1,
!     +            'Illegal constraint --- all constraints ignored.')
                endif
!--- set line width scale factor
                klwid = max(1, 6 - abs(k1) - abs(k2))
                sclwid = klwid
                call jslwsc(sclwid)
                if (k1*k2 .lt. 0)  then
                  v1 = k(kx) * ql(kx) + k(ky) * qu(ky) + k3 * qs
                  v2 = k(kx) * qu(kx) + k(ky) * ql(ky) + k3 * qs
                else
                  v1 = k(kx) * qu(kx) + k(ky) * qu(ky) + k3 * qs
                  v2 = k(kx) * ql(kx) + k(ky) * ql(ky) + k3 * qs
                endif
                vm = min(v1, v2)
                vp = max(v1, v2)
                ipm = vm / s + 0.499 * (1. + sign(1., vm))
                ipp = vp / s - 0.499 * (1. - sign(1., vm))
                do ip = max(0, ipm), ipp
                  if (k(2) .eq. 0)  then
!--- vertical line
                    y(1) = ql(2)
                    y(2) = qu(2)
                    x(1) = (ip * s - k3 * qs) / k(1)
                    x(2) = x(1)
                  else
                    x(1) = ql(1)
                    x(2) = qu(1)
                    y(1) = (ip * s - k(1) * ql(1) - k3 * qs) / k(2)
                    y(2) = (ip * s - k(1) * qu(1) - k3 * qs) / k(2)
                    if (k(1) .ne. 0)  then
                      do j = 1, 2
                        if (y(j) .lt. ql(2))  then
                          x(j) = (ip * s - k(2) * ql(2) - k3 * qs) / k  &
     &                    (1)
                          y(j) = ql(2)
                        elseif (y(j) .gt. qu(2))  then
                          x(j) = (ip * s - k(2) * qu(2) - k3 * qs) / k  &
     &                    (1)
                          y(j) = qu(2)
                        endif
                      enddo
                    endif
                  endif
                  call gvpl(2,x,y)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!--- reset line width scale factor
      call jslwsc(1.)
      end
      subroutine pecat1(rb, ra, rd)
!----------------------------------------------------------------------*
! purpose:
!   concatenate two transport maps
! input:
!   rb(6,6) second map in beam line order.
!   ra(6,6) first map in beam line order.
! output:
!   rd(6,6) result map.
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      double precision rb(6,6), ra(6,6), rd(6,6)
      integer k, j

      do k = 1, 6
        do j = 1, 6
          rd(j,k) = rb(j,1) * ra(1,k) + rb(j,2) * ra(2,k)               &
     &    + rb(j,3) * ra(3,k) + rb(j,4) * ra(4,k)                       &
     &    + rb(j,5) * ra(5,k) + rb(j,6) * ra(6,k)
        enddo
      enddo
      end
      subroutine pecurv (ncc, spname, annh, usex, sych, ippar,          &
     &np, xval, yval, window, actwin, ierr)
!----------------------------------------------------------------------*
! Purpose:
!   Plot one curve
! Input:
!   NCC      (integer)  current curve count (1,2, etc.)
!   SPNAME    (char)     curve annotation string
!   ANNH     (real)     character height
!   USEX     (real)     user character height expansion
!   SYCH     (real)     symbol character height
!   IPPAR    (integer)  array containing the plot parameters
!   NP       (integer)  no. of points to plot
!   XVAL     (real)     x values
!   YVAL     (real)     y values
!   WINDOW   (real)     array containing the window to use
!   ACTWIN   (real)     active (inside frame) window
! Output:
!   IERR     (integer)  0 if OK, else GXPLOT error
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      real annh, sych, xval(*), yval(*), window(*)
      real actwin(*), gxcubv, usex
      real xpos, ypos, xwpos, ywpos, xf, wsclx, wscly
      real xaux, yaux
      real rsave(20), act(4), xpl(2), ypl(2), winnor(4)
      real xreal(maxpnt), yreal(maxpnt), yy1d(maxpnt)
      real yy2d(maxpnt)
      real xmd
      integer ippar(*)
      integer isymb, icolr, ilb, ispli, kft, klt, kact,                 &
     &kf, kl, npt, j, iecub
      integer              ncc, np, ierr, istyl, ipbar, isave(20)
      character         spname*(*), symloc*1
      save                 act
      data winnor /0., 1., 0., 1./

!--- save GKS settings
      call gxsave (isave, rsave, ierr)
      call gxswnd (window)
      wsclx = 1. / (window(2) - window(1))
      wscly = 1. / (window(4) - window(3))
      xmd = 1.e-8 * (window(2) - window(1))**2
      if (ncc .eq. 1)  then
!--- first curve in frame - reset label position array, get act.
!    window in NDC
        act(1) = (actwin(1) - window(1)) * wsclx
        act(2) = (actwin(2) - window(1)) * wsclx
        act(3) = (actwin(3) - window(3)) * wscly
        act(4) = (actwin(4) - window(3)) * wscly
      endif
      istyl = ippar(1)
      ipbar = ippar(3)
      isymb = ippar(4)
      icolr = ippar(5)
      if(icolr .eq. 100)  icolr = mod(ncc-1,4) + 1
      icolr = max(1, min(icolr, 7))
      ilb = -1
      if (istyl .ne. 0)  then
!--- polyline requested
        if (np .lt. 2) goto 999
        if(istyl .eq. 100)  istyl = mod(ncc-1,4) + 1
        ispli = ippar(2)
!--- get first and last blank in annotation
        call gxpnbl (spname, kft, klt)
        if (kft .ne. 0) then
!--- annotation exists
          ilb    = 0
        endif
!--- set line style
        call jsln (max (1, min (4, istyl)))
!--- set line colour
        call jsplci(icolr)
        kact = 1
   10   continue
!--- get first and last point inside
        call peiact(kact, np, xval, yval, actwin, kf, kl)
!--- quit if no points inside
        if (kf .eq. 0) goto 40
        kf = max(1, kf - 1)
        kl = min(np, kl + 1)
        npt      = 1
        xreal(1) = xval(kf)
        yreal(1) = yval(kf)
        do j = kf + 1, kl
!--- avoid identical points
          if ((xreal(npt) - xval(j))**2 +                               &
     &    (yreal(npt) - yval(j))**2 .gt. xmd) then
            npt        = npt + 1
            xreal(npt) = xval(j)
            yreal(npt) = yval(j)
          endif
          if ((j .eq. kl .and. npt .ge. 2) .or. npt .eq. maxpnt) then
!--- plot - get first curve annotation position
            if (ilb .eq. 0)                                             &
     &      call pegacn(ncc, window, act, xreal, yreal, npt, usex,      &
     &      xwpos, xpos, ypos, ilb)
            if (splinef .eq. 0 .or. npt .eq. 2 .or. ispli .ne. 0) then
!--- no spline
              call gxpl (npt, xreal, yreal, actwin)
              if (ilb .gt. 0) then
!--- get y pos. on curve for label
                ywpos = yreal(ilb - 1) + (yreal(ilb) - yreal(ilb-1))    &
     &          * (xwpos  - xreal(ilb - 1))                             &
     &          / (xreal(ilb) - xreal(ilb - 1))
                ilb = -2
              endif
            else
!--- spline
              call gxplt1 (npt, xreal, yreal, actwin)
              if (ilb .gt. 0) then
!--- get y pos. on curve for label
!                call gxcubi (npt, xreal, yreal, yy1d, yy2d, iecub)
!                ywpos = gxcubv (xwpos, npt, xreal, yreal, yy1d, yy2d)
                ilb   = -2
              endif
            endif
            xreal(1) = xreal(npt)
            yreal(1) = yreal(npt)
            npt = 1
          endif
        enddo
        if (kl .lt. np)  then
          kact = kl + 1
          goto 10
        endif
      else
!--- no polyline
        if (np .eq. 0) goto 999
      endif
!--- plot symbols or bars if requested
      if (ipbar .ne. 0)  then
        call jsln (1)
!--- set line colour
        call jsplci(icolr)
        do j = 1, np
          xpl(1) = xval(j)
          xpl(2) = xval(j)
          ypl(1) = yval(j)
          ypl(2) = actwin(3)
          call gvpl (2, xpl, ypl)
        enddo
      endif
   40 continue
      if (isymb .ne. 0)  then
        if (isymb .le. 5)  then
          call jsmk (isymb)
          call gxpmsw (np, xval, yval, actwin)
        elseif (isymb .eq. 100)  then
          if (istyl .ne. 0)  then
!--- use current curve count
            write (symloc, '(I1)')  mod (ncc, 10)
          endif
!--- set marker colour
          call jspmci(icolr)
!--- plot one character symbol
!    switch to normalized window
          call gxswnd (winnor)
!--- set character height
          call jschh (sych)
!--- text alignment
          call jstxal (2, 3)
!--- text expansion factor - mind distorted viewports
          call gxqrvp (xf)
          call jschxp (xf)
          do j = 1, np
            if (isymb .eq. 100 .and. istyl .eq. 0)  then
!--- use current point number
              write (symloc, '(I1)')  mod (j, 10)
            endif
            xaux = wsclx * (xval(j) - window(1))
            yaux = wscly * (yval(j) - window(3))
            if (xaux .gt. act(1) .and. xaux .lt. act(2)                 &
     &      .and. yaux .gt. act(3) .and. yaux .lt. act(4))              &
     &      call gxtx (xaux, yaux, symloc)
          enddo
        endif
      endif
      if (ilb .eq. -2)  then
!--- plot annotation
!    switch to normalized window
        call gxswnd (winnor)
!--- set character height
        call jschh (annh)
!--- text alignment
        call jstxal (2, 5)
!--- text expansion factor - mind distorted viewports
        call gxqrvp (xf)
        call jschxp (xf)
!--- set marker colour
        call jstxci(icolr)
!--- plot annotation string
        call gxtx (xpos, ypos, spname(kft:klt))
!--- connect to curve
        xpl(1) = xpos
        xpl(2) = xpos
        ypl(1) = ypos
        ypl(2) = (ywpos - window(3)) * wscly
        if (ypl(2) .gt. ypl(1))  ypl(1) = ypl(1) + .02
!--- set dotted line
        call jsln (3)
!--- set line colour
        call jsplci(icolr)
!--- plot line
        call gxpl (2, xpl, ypl, act)
      endif
!--- restore
      call gxrest (isave, rsave)
  999 continue
      end
      subroutine peelma(pos, itp, temp, dstp, am)
!----------------------------------------------------------------------*
! purpose:
!   return 4x6 element matrix
! input:
!   pos         (int)   position: 1 entry, 2 body, 3 exit
!   itp         (int)   type:
!                       1 straight sbend,
!                       7 vert. sbend,
!                       21 straight rbend,
!                       27 vert. rbend,
!                       2 foc. quad,
!                       3 defoc. quad
!                         else drift
!   temp        (d.p.)  element parameters
!       1    l      [m]
!       2    rhoinv [1/m]
!       3    e1
!       4    e2
!       5    h1
!       6    h2
!       7    tilt
!       8    ks
!       9    hgap [m]
!      10    fint [Tm]
!      11    angle = K_0*l
!      12    lrad
!      13    k0 or k0*l (l=0)
!      14    k0s or k0s*l
!      15    k1 or k1*l (l=0)
!      16    k1s or k1s*l
!      17    k2 or k2*l
!      18    k2s  etc.
!   dstp        (d.p.)  sublength
! output:
!   am          (d.p.)  matrix
!----------------------------------------------------------------------*
      implicit none
      include 'twissl.fi'
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      double precision pi, zero, eps, one, two, twopi
      parameter         (pi = 3.1415926535898d0)
      parameter         (zero = 0.d0, eps = 1.d-5)
      parameter         (one = 1.d0, two = 2.d0, twopi = two * pi)
      integer pos, itp
      double precision temp(*), dstp, am(6,6)

      double precision e2, rk2, cangle, sangle, tilt, d, f 
      double precision ek(6), te(6,6,6), orbit(6), am2(6,6)
      double precision node_value
      integer i, j, k
      logical drift, bend, quad, fmap,ftrk,cecptk,cebttk

!--- Initialise

      ftrk = .false.
      do i = 1, 6
        do j = 1, 6
          am(j,i) = zero
          am2(j,i) = zero
        enddo
        am(i,i) = one
        am2(i,i) = one
        orbit(i) = zero
      enddo

!---  save values of centre_cptk & centre_bttk and set them to .false.

      cecptk = centre_cptk
      cebttk = centre_bttk
      centre_cptk = .false.
      centre_bttk = .false.

!--- End of initialisation

!--- treat quads with non-zero tilt, and bends with non-zero tilt
!    different from pi/2 as drifts

      k = mod(itp, 20)
      bend = k .eq. 1 .or. k .eq. 7
      quad = itp .eq. 2 .or. itp .eq. 3
      drift = itp .eq. 0 .or. .not. (bend .or. quad)
      if (drift)  then
        am(1,2) = dstp
        am(3,4) = dstp
      else
         if (bend) then
            if (pos .eq. 1) then
              e2 = node_value('e2 ')
              call store_node_value('e2 ', zero)
            endif  
            if (pos .eq. 3) call store_node_value('e2 ', e2)
            call tmbend(.false.,orbit,fmap,dstp,ek,am,te)
         else

!--- quad -------------------------------------------

            rk2  = temp(15)
            tilt = temp(7)
            call tmquad(.false.,.false.,tilt,orbit,fmap,dstp,ek,am,te)
            call tmfoc(dstp, rk2, cangle, sangle, d, f)
            am2(1,1) = cangle
            am2(2,2) = cangle
            am2(1,2) = sangle
            am2(2,1) = -rk2 * sangle
            rk2 = -rk2
            call tmfoc(dstp, rk2, cangle, sangle, d, f)
            am2(3,3) = cangle
            am2(4,4) = cangle
            am2(3,4) = sangle
            am2(4,3) = -rk2 * sangle
         endif   
      endif
   
!---  restore values of centre_cptk & centre_bttk 

      centre_cptk = cecptk
      centre_bttk = cebttk

      end
 
      subroutine pefill(ierr)

#ifdef _WIN32
      include 'win32calls.fi'
#endif


!----------------------------------------------------------------------*
! Purpose:
!   fill plot arrays with coordinate values, set up machine plot
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer mtype, m_adble
      parameter (mtype = 50, m_adble = 20)
      double precision zero, half
      parameter (zero = 0.d0, half = 0.5d0)
      integer ierr,i,j,k,l,ndble,new1,new2,currtyp,nint,crow,pltyp
      integer double_from_table, advance_to_pos,                        &
     &restart_sequ
      integer ilist(mtype), p(mxplot), proc_n(2, mxcurv)
      integer lastnb, table_length
      double precision fact, currpos, currleng, node_value,             &
     &d_val, d_val1, adble(m_adble)
      real tval, step
      logical machp, rselect
!--- codes see in peschm
      data ilist /                                                      &
     &0, 21, 1, 0, 2, 10, 12, 8, 0, 9,                                  &
     &6, 0, 0, 0, 0,  0,  4, 4, 4, 0,                                   &
     &0, 0, 0, 0, 0,  0, 14, 0, 0, 0,                                   &
     &20 * 0 /

      do i=1,mxplot
        p(i)=0
      enddo

      ierr = 0
      k = double_from_table(tabname, horname, 1, d_val)
      if (k .lt. 0)  then
        if (k .eq. -1)  then
          print *, 'Warning: table ', tabname, ' not found'
        elseif (k .eq. -2)  then
          print *, 'Warning: hor. variable ', horname,                  &
     &    ' not in table ', tabname
        else
          print *, 'Warning: table ', tabname, ' is empty'
        endif
        ierr = 1
        goto 999
      endif
!--- save process flags, proc_flag may be modified in peintp
      do i = 1, mxcurv
        do j = 1, 2
          proc_n(j,i) = proc_flag(j,i)
        enddo
      enddo
      machp = itbv .ne. 0
      rselect = itbv .ne. 0 .and. hrange(2) .gt. hrange(1)
      do l = 1, nivvar
        k = double_from_table(tabname, sname(l), 1, d_val)
        if (k .lt. 0)  then
          print *, 'Warning: vertical variable: ',                      &
     &    sname(l)(:lastnb(sname(l))), ' not in table ',                &
     &    tabname
          ierr = 1
          goto 999
        endif
      enddo
      if (rselect)  then
!------- adjust element range to horizontal range
        new1 = nrrang(1)
        new2 = nrrang(2)
        crow = nrrang(1)
        do j = nrrang(1), nrrang(2)
          k = double_from_table(tabname, horname, j, d_val)
          tval = d_val
          if (tval .lt. hrange(1)) new1 = j
          if (tval .lt. hrange(2)) new2 = j
        enddo
        nrrang(1) = new1
        if (nrrang(2) .gt. new2+2) nrrang(2) = new2 + 2
      endif
      if (itbv .eq. 0)  then
        nrrang(1) = 1
        nrrang(2) = table_length(tabname)
      endif
      if (nrrang(1) .eq. 0) nrrang(1) = 1
!--- get interpolation interval size
      if (machp)  then
        k = double_from_table(tabname, horname, nrrang(1), d_val)
        k = double_from_table(tabname, horname, nrrang(2), d_val1)
        step = (d_val1 - d_val) / (maxpnt / 2)
      endif
      do j = 1, mxcurv
        nqval(j) = 0
      enddo
      nelmach = 0
      pos_flag = 3
      fact = pos_flag - 1
      j = restart_sequ()
      do j = nrrang(1), nrrang(2)
        if (itbv .eq. 1 .and. advance_to_pos(tabname, j) .eq. 0)        &
     &  goto 10
        k = double_from_table(tabname, horname, j, currpos)
        if (itbv .eq. 1)  then
          currtyp = node_value('mad8_type ')
          if (currtyp .le. mtype) then
            pltyp = ilist(currtyp)
          else
            pltyp = 0
          endif
!--- get element parameters
          ndble = m_adble
          call pelfill(tabname, currtyp, j, ndble, adble)
          currleng = adble(1)
!--- sbend, rbend
          if (mod(pltyp,20) .eq. 1 .and. adble(7) .ne. zero)            &
     &    pltyp = pltyp + 6
!--- quad
          if (pltyp .eq. 2 .and. min(adble(15), adble(16)) .lt. zero)   &
     &    pltyp = 3
!--- sext
          if (pltyp .eq. 10 .and. min(adble(17), adble(18)) .lt. zero)  &
     &    pltyp = 11
!--- oct
          if (pltyp .eq. 12 .and. min(adble(19), adble(20)) .lt. zero)  &
     &    pltyp = 13
          if (pltyp .gt. 0 .and. machp)  then
            nelmach = nelmach + 1
            estart(nelmach) = currpos - fact * half * currleng
            eend(nelmach)   = estart(nelmach) + currleng
            ieltyp(nelmach) = pltyp
          endif
          if (machp .and. j .gt. nrrang(1) .and. currleng .gt. zero     &
     &    .and. ipparm(2,1) .gt. 0 .and. step .gt. zero)  then
            nint = currleng / step
            if (nint .lt. 2) nint = 2
            call peintp(crow, nint, pltyp, proc_n, adble, ierr)
            if (ierr .eq. 1)  then
              ierr = 0
              print *, 'Warning: plot buffer full, plot truncated'
              goto 100
            elseif (ierr .ne. 0)  then
              goto 999
            endif
          endif
        endif
        do l = 1, nivvar
          if (nqval(l) .eq. maxseql)  then
            print *, 'Warning: plot buffer full, plot truncated'
            goto 100
          elseif (nqval(l) .eq. 0)  then
            nqval(l) = nqval(l) + 1
            qhval(nqval(l),l) = currpos
            k = double_from_table(tabname, sname(l), j, d_val)
            k = p(l)
            qvval(nqval(l),l) = d_val
            if (proc_flag(1,l) .eq. 1) then
              qvval(nqval(l),l) = sqrt(abs(qvval(nqval(l),l)))
            endif
          elseif (itbv .eq. 0 .or. currpos - qhval(nqval(l),l)          &
     &    .gt. 0.1d0 * step) then
            nqval(l) = nqval(l) + 1
            qhval(nqval(l),l) = currpos
            k = double_from_table(tabname, sname(l), j, d_val)
            k = p(l)
            qvval(nqval(l),l) = d_val
            if (proc_flag(1,l) .eq. 1) then
              qvval(nqval(l),l) = sqrt(abs(qvval(nqval(l),l)))
            endif
          endif
        enddo
        crow = j
      enddo
   10 continue
  100 continue
      fpmach = machp .and. nelmach .gt. 0 .and. noline .eq. 0
  999 continue
      end
      subroutine pegacn(ncc, window, act, xreal, yreal, np, usex,       &
     &xwpos, xpos, ypos, ilb)
!----------------------------------------------------------------------*
! Purpose:
!   Find suitable position for the curve annotation
! Input:
!   NCC      (integer)  current curve count (1,2, etc.)
!   WINDOW   (real)     array containing the window to use
!   ACT      (real)     window in NDC
!   XREAL    (real)     x values of curve
!   YREAL    (real)     y values of curve
!   NP       (integer)  no. of points to plot
!   USEX     (real)     user character height expansion
! Output:
!   XWPOS    (real)     x position of label in world coords.
!   XPOS     (real)     x pos. of label in NDC
!   YPOS     (real)     y pos. of label in NDC
!   ILB      (integer)  number of point behind label, or 0 if no
!                       label possible
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer mposx, mposy, mpost
      parameter (mposx = 8, mposy = 3, mpost = mposx * mposy)

      real window(*), act(*), xreal(*), yreal(*)
      real xwpos, xpos, ypos, usex
      real ywpos, xmax, xmin, xdiff, ydiff, d, t, eps
      real xdiag(2,2), ydiag(2,2)
      real xadd, yadd
      integer ncc, np, ilb, i, iapos, iposx, iposy, j, iy
      integer iucomp, kapos(mposx, mposy)
      save kapos

!--- reset position array if first curve in frame
      if (ncc .eq. 1)  then
        do i = 1, mposx
          do j = 1, mposy
            kapos(i,j) = 0
          enddo
        enddo
      endif
      xdiff = window(2) - window(1)
      ydiff = window(4) - window(3)
      eps = 1.e-6 * max(xdiff, ydiff)
      xmax  = xreal(1)
      xmin  = xmax
      do i = 2, np
        xmin = min(xmin, xreal(i))
        xmax = max(xmax, xreal(i))
      enddo
   20 continue
!--- find first unoccupied position
      iapos  = iucomp(0, kapos, mpost)
      if (iapos .eq. 0)  then
        ilb = 0
      else
        iposx  = mod (iapos-1, mposx) + 1
        iposy  = (iapos-1) / mposx + 1
        kapos(iposx,iposy) = -1
!--- annot. pos. in NDC
        xpos = act(1) +                                                 &
     &  0.125 * usex * (iposx - .5) * (act(2) - act(1))
        ypos = act(4) -                                                 &
     &  usex * (0.05 * (act(4) - act(3)) + 0.03 * (iposy - 1))
!---- annot. position in world coord.
        xwpos = window(1) + xpos * xdiff
!--- get next if outside x values of curve
        if (xwpos .le. xmin .or. xwpos .gt. xmax) goto 20
        ywpos = window(3) + ypos * ydiff
!--- get endpoint of both diagonals of box
        xadd = 0.0625 * xdiff
        yadd = 0.03  * ydiff
        xdiag(1,1) = xwpos - xadd
        xdiag(2,1) = xwpos + xadd
        xdiag(1,2) = xwpos - xadd
        xdiag(2,2) = xwpos + xadd
        ydiag(1,1) = ywpos
        ydiag(2,1) = ywpos + yadd
        ydiag(1,2) = ywpos + yadd
        ydiag(2,2) = ywpos
!--- make sure no part of curve cuts these lines (curve approx. by
!    straight line segments)
        do i = 2, np
          if (xwpos .gt. xreal(i-1) .and. xwpos .le. xreal(i)) ilb = i
          do j = 1, 2
            d = (xdiag(2,j) - xdiag(1,j)) * (yreal(i-1) - yreal(i)) -   &
     &      (ydiag(2,j) - ydiag(1,j)) * (xreal(i-1) - xreal(i))
            if (abs(d) .lt. eps) goto 30
            t = (xreal(i-1) - xdiag(1,j)) * (yreal(i-1) - yreal(i)) -   &
     &      (yreal(i-1) - ydiag(1,j)) * (xreal(i-1) - xreal(i))
            t = t / d
            if (t .lt. 0. .or. t .gt. 1.) goto 30
            t = (xdiag(2,j) - xdiag(1,j)) * (yreal(i-1) - ydiag(1,j)) - &
     &      (ydiag(2,j) - ydiag(1,j)) * (xreal(i-1) - xdiag(1,j))
            t = t / d
            if (t .ge. 0. .and. t .le. 1.) goto 20
   30       continue
          enddo
        enddo
      endif
      if (ilb .gt. 0)  then
        do iy = 1, mposy
          kapos(iposx, iy) = 1
        enddo
      endif
      do i = 1, mposx
        do j = 1, mposy
          kapos(i,j) = max(0, kapos(i,j))
        enddo
      enddo
      end
      subroutine pegaxn (nax, vax, sax, ns)
!----------------------------------------------------------------------*
! Purpose:
!   Returns compound vertical axis annotation
!
!--- Input
!   NAX       (integer) no. of vert. var. names in VAX
!   VAX          (char) vert. var. names
!---Output
!   SAX          (char) remaining (possibly truncated) names
!   NS        (integer) no. of names in SAX
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      character * 16 vax(*), sax(*), scut, saloc
      integer nax, ns, i, k, k1, k2, j, k1f, k2f
      ns = 0
      if (nax .le. 0)  then
        sax(1) = ' '
      else
        do  i = 1, nax
          saloc = vax(i)
          call gxpnbl(saloc, k1, k2)
          if (k2 .gt. 1 .and. index('XY', saloc(k2:k2)) .ne. 0)  then
            scut = saloc(:k2-1)
            do j = 1, ns
              if (scut .eq. sax(j))  goto 10
            enddo
            do j = i + 1, nax
              call gxpnbl(vax(j), k1f, k2f)
              if (k2 .eq. k2f)  then
                if (index('XY', vax(j)(k2:k2)) .ne. 0)  then
                  if (saloc(:k2-1) .eq. vax(j)(:k2-1))  then
                    saloc = scut
                    do k = 1, ns
                      if (saloc .eq. sax(k))  goto 10
                    enddo
                  endif
                endif
              endif
            enddo
          endif
          ns      = ns + 1
          sax(ns) = saloc
   10     continue
        enddo
      endif
      end
      subroutine  pegetn (iflag, svar, it, ipflg, sovar, sname)
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Finds variable, dependent variables, axis and curve annotations    *
! Input:                                                               *
!   IFLAG    (integer)  0 for dependent variables and process flag,    *
!                       1 for axis, 2 for curve, 3 for trunc. name,    *
!                       4 to print the axis names on IQLOG             *
!   SVAR        (char)  variable to be looked up.                      *
!   IT          (int)   table number (see PLGTBS).                     *
! Output:                                                              *
!   IPFLG(1) (integer)  process flag: 0 as is, 1 take root, else call  *
!                       function PLPVAL                                *
!   IPFLG(2) (integer)  interpol. flag: 0 spline, else call            *
!                       function PEINTP                                *
!   SOVAR       (char)  array of (up to MXDEP) dependent variables     *
!   SNAME       (char)  requested annotation                           *
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer iflag, it, ipflg(2)
      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/

      character         svar*(mcnam), sovar(*)*(mcnam),                 &
     &sname * (*)
      character         svlabl(mnvar)*(mxlabl)
      character         svanno(mnvar)*(mxlabl)
      character         svname(mnvar)*(mcnam)
!--- strings:
!   SVLABL   plot prescriptions for variables on axis labels
!   SVANNO   plot prescriptions for variables in annotations
!   SVNAME   names of variables known to the program
      integer i, iref, k1, k2, k1f, k2f, j
      integer              iproc(mnvar,3), intpo(mnvar)
      integer              ivdep(mnvar,mxdep,3)

      data (svname(j), j = 1, 32) /                                     &
     &'s', 'size', 'deltap',                                            &
     &'qs', 'x', 'y', 'xsize', 'ysize',                                 &
     &'dt', 'xn', 'yn', 'pxn', 'pyn',                                   &
     &'gammatr', 'xrms', 'yrms',                                        &
     &'xmax', 'ymax', 'bxmax',                                          &
     &'bymax', 'dxmax', 'dymax',                                        &
     &'tn', 't', 'turns', 'particle', 'alfa',                           &
     &'ptn', 'wt', 'phit',                                              &
     &'rbxmax',                                                         &
     &'rbymax' /
      data (svname(j), j = 33, mnvar) /                                 &
     &'betx', 'rbetx',                                                  &
     &'alfx', 'mux', 'dx',                                              &
     &'dpx', 'qx', 'px', 'wx',                                          &
     &'phix', 'dmux',                                                   &
     &'ddx', 'ddpx', 'iwx',                                             &
     &'xix',                                                            &
     &'bety', 'rbety',                                                  &
     &'alfy', 'muy', 'dy',                                              &
     &'dpy', 'qy', 'py', 'wy',                                          &
     &'phiy', 'dmuy',                                                   &
     &'ddy', 'ddpy', 'iwy',                                             &
     &'xiy', 'xns', 'pxns', 'wxs',                                      &
     &'yns', 'pyns', 'wys',                                             &
     &'energy', 'spintune',                                             &
     &'poltotal', 'poldiffx', 'poldiffy', 'poldiffs' /

      data (svlabl(j), j = 1, 32) /                                     &
     &'s (m)', 'n<G>s<G> (mm)', '<G>d<G><?>E<?>/p<?>0<?>c',             &
     &'Q<?>s<?>', 'x (m)', 'y (m)', 'n<G>s<G> (mm)', 'n<G>s<G> (mm)',   &
     &'ct (m)', 'x<?>n<?>', 'y<?>n<?>', 'p<?>xn<?>', 'p<?>yn<?>',       &
     &'<G>g<G><?>tr<?>', 'X<?>rms<?> (m)', 'Y<?>rms<?> (m)',            &
     &'X<?>max<?> (m)', 'Y<?>max<?> (m)', '<G>b<G><?>x_max<?> (m)',     &
     &'<G>b<G><?>y_max<?> (m)', 'D<?>x_max<?> (m)', 'D<?>y_max<?> (m)', &
     &'t<?>n<?>', 'ct (m)', 'turns', 'particle', '<G>a<G>',             &
     &'p<?>t_n<?>', 'W<?>t<?>', '<G>F<G><?>t<?> (rad/2<G>p<G>)',        &
     &'<G>b<G><?>x_max<?><!>1/2<!> (m<!>1/2<!>)',                       &
     &'<G>b<G><?>y_max<?><!>1/2<!> (m<!>1/2<!>)' /
      data (svlabl(j), j = 33, mnvar) /                                 &
     &'<G>b<G><?>x<?> (m)', '<G>b<G><?>x<?><!>1/2<!> (m<!>1/2<!>)',     &
     &'<G>a<G><?>x<?>', '<G>m<G><?>x<?> (rad/2<G>p<G>)', 'D<?>x<?> (m)',&
     &'D<?>px<?>', 'Q<?>x<?>', 'p<?>x<?>/p<?>0<?>', 'W<?>x<?>',         &
     &'<G>F<G><?>x<?> (rad/2<G>p<G>)', 'd<G>m<G><?>x<?>/d<G>d<G>',      &
     &'dD<?>x<?>/d<G>d<D> (m)', 'dD<?>px<?>/d<G>d<G>', 'W<?>x<?> (m)',  &
     &'XI<?>x<?>',                                                      &
     &'<G>b<G><?>y<?> (m)', '<G>b<G><?>y<?><!>1/2<!> (m<!>1/2<!>)',     &
     &'<G>a<G><?>y<?>', '<G>m<G><?>y<?> (rad/2<G>p<G>)', 'D<?>y<?> (m)',&
     &'D<?>py<?>', 'Q<?>y<?>', 'p<?>y<?>/p<?>0<?>', 'W<?>y<?>',         &
     &'<G>F<G><?>y<?> (rad/2<G>p<G>)', 'd<G>m<G><?>y<?>/d<G>d<G>',      &
     &'dD<?>y<?>/d<G>d<D> (m)', 'dD<?>py<?>/d<G>d<G>', 'W<?>y<?> (m)',  &
     &'XI<?>y<?>', 'x<?>ns<?>', 'p<?>x_ns<?>', 'W<?>xs<?>',             &
     &'y<?>ns<?>', 'p<?>y_ns<?>', 'W<?>ys<?>',                          &
     &'E[GeV]', 'spintune',                                             &
     &'polarization','polarization','polarization','polarization' /

      data (svanno(j), j = 1, 32) /                                     &
     &'s', 'n<G>s<G>', '<G>d<G>',                                       &
     &'Q<?>s<?>', 'x', 'y', 'n<G>s<G><?>x<?>', 'n<G>s<G><?>y<?>',       &
     &'ct', 'x<?>n<?>', 'y<?>n<?>', 'p<?>xn<?>', 'p<?>yn<?>',           &
     &'<G>g<G><?>tr<?>', 'X<?>rms<?>', 'Y<?>rms<?>',                    &
     &'X<?>max<?>', 'Y<?>max<?>', '<G>b<G><?>x_max<?>',                 &
     &'<G>b<G><?>y_max<?>', 'D<?>x_max<?>', 'D<?>y_max<?>',             &
     &'t<?>n<?>', 't', 'turns', 'particle', '<G>a<G>',                  &
     &'p<?>t_n<?>', 'W<?>t<?>', '<G>F<G><?>t<?>',                       &
     &'<G>b<G><?>x_max<?><!>1/2<!>',                                    &
     &'<G>b<G><?>y_max<?><!>1/2<!>' /
      data (svanno(j), j = 33, mnvar) /                                 &
     &'<G>b<G><?>x<?>', '<G>b<G><?>x<?><!>1/2<!>',                      &
     &'<G>a<G><?>x<?>', '<G>m<G><?>x<?>', 'D<?>x<?>',                   &
     &'D<?>px<?>', 'Q<?>x<?>', 'p<?>x<?>', 'W<?>x<?>',                  &
     &'<G>F<G><?>x<?>', '<G>m<G><?>x<?>''',                             &
     &'D<?>x<?>''', 'D<?>px<?>''', 'W<?>x<?>',                          &
     &'XI<?>x<?>',                                                      &
     &'<G>b<G><?>y<?>', '<G>b<G><?>y<?><!>1/2<!>',                      &
     &'<G>a<G><?>y<?>', '<G>m<G><?>y<?>', 'D<?>y<?>',                   &
     &'D<?>py<?>', 'Q<?>y<?>', 'p<?>y<?>', 'W<?>y<?>',                  &
     &'<G>F<G><?>y<?>', '<G>m<G><?>y<?>''',                             &
     &'D<?>y<?>''', 'D<?>py<?>''', 'W<?>y<?>',                          &
     &'XI<?>y<?>', 'x<?>ns<?>', 'p<?>x_ns<?>', 'W<?>xs<?>',             &
     &'y<?>ns<?>', 'p<?>y_ns<?>', 'W<?>ys<?>',                          &
     &' ', ' ',                                                         &
     &'p<?>tot<?>', 'p<?>diff_x<?>', 'p<?>diff_y<?>', 'p<?>diff_s<?>'/

      data (iproc(j,1), j = 1, 32) /                                    &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 2, 3, 4, 5,                                                    &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &6, 0, 0, 0, 0,                                                    &
     &7, 8, 9,                                                          &
     &1,                                                                &
     &1 /

      data (iproc(j,1), j = 33, mnvar) /                                &
     &0, 1,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0,                                                       &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &0, 1,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0,                                                       &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 14, 15, 16,                                                    &
     &17, 18, 19,                                                       &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      data (iproc(j,2), j = 1, 32) /                                    &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 2, 3, 4, 5,                                                    &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &6, 0, 0, 0, 0,                                                    &
     &7, 8, 9,                                                          &
     &1,                                                                &
     &1 /

      data (iproc(j,2), j = 33, mnvar) /                                &
     &0, 1,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 10,                                                      &
     &11, 0,                                                            &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &0, 1,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 12,                                                      &
     &13, 0,                                                            &
     &0, 0, 0,                                                          &
     &0, 14, 15, 16,                                                    &
     &17, 18, 19,                                                       &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      data (iproc(j,3), j = 1, 32) /                                    &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 2, 3, 4, 5,                                                    &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &6, 0, 0, 0, 0,                                                    &
     &7, 8, 9,                                                          &
     &1,                                                                &
     &1 /

      data (iproc(j,3), j = 33, mnvar) /                                &
     &0, 1,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 10,                                                      &
     &11, 0,                                                            &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &0, 1,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 12,                                                      &
     &13, 0,                                                            &
     &0, 0, 0,                                                          &
     &0, 14, 15, 16,                                                    &
     &17, 18, 19,                                                       &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      data (intpo(j), j = 1, 32) / 32 * 0 /
!--- in INTPO, n+100 means: take SQRT of var. n
      data (intpo(j), j = 33, mnvar) /                                  &
     &1, 101,                                                           &
     &2, 3, 4,                                                          &
     &5, 0, 0, 0,                                                       &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &6, 106,                                                           &
     &7, 8, 9,                                                          &
     &10, 0, 0, 0,                                                      &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0,                                                       &
     &0, 0, 0,                                                          &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      data (ivdep(j,1,1), j = 1, 32) /                                  &
     &1, 2, 3,                                                          &
     &4, 5, 6, 7, 8,                                                    &
     &9, 5, 6, 5, 6,                                                    &
     &14, 15, 16,                                                       &
     &17, 18, 19,                                                       &
     &20, 21, 22,                                                       &
     &24, 24, 25, 26, 27,                                               &
     &3, 3, 3,                                                          &
     &19,                                                               &
     &20 /
      data (ivdep(j,2,1), j = 1, 32) /                                  &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 0, 0, 40, 55,                                                  &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 24, 24,                                                        &
     &0,                                                                &
     &0 /
      data (ivdep(j,1,1), j = 33, mnvar) /                              &
     &33, 33,                                                           &
     &35, 36, 37,                                                       &
     &38, 39, 40, 41,                                                   &
     &42, 43,                                                           &
     &44, 45, 46,                                                       &
     &47,                                                               &
     &48, 48,                                                           &
     &50, 51, 52,                                                       &
     &53, 54, 55, 56,                                                   &
     &57, 58,                                                           &
     &59, 60, 61,                                                       &
     &62, 5, 5, 5,                                                      &
     &6, 6, 6,                                                          &
     &69, 70,                                                           &
     &71, 72, 73, 74 /
      data (ivdep(j,2,1), j = 33, mnvar) /                              &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0,                                                       &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0,                                                       &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 40, 40,                                                     &
     &0, 55, 55,                                                        &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      data (ivdep(j,1,2), j = 1, 32) /                                  &
     &1, 2, 3,                                                          &
     &4, 5, 6, 7, 8,                                                    &
     &9, 5, 6, 5, 6,                                                    &
     &14, 15, 16,                                                       &
     &17, 18, 19,                                                       &
     &20, 21, 22,                                                       &
     &24, 24, 25, 26, 27,                                               &
     &3, 3, 3,                                                          &
     &19,                                                               &
     &20 /
      data (ivdep(j,2,2), j = 1, 32) /                                  &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 0, 0, 40, 55,                                                  &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 24, 24,                                                        &
     &0,                                                                &
     &0 /
      data (ivdep(j,1,2), j = 33, mnvar) /                              &
     &33, 33,                                                           &
     &35, 36, 37,                                                       &
     &38, 39, 40, 5,                                                    &
     &5, 43,                                                            &
     &44, 45, 46,                                                       &
     &47,                                                               &
     &48, 48,                                                           &
     &50, 51, 52,                                                       &
     &53, 54, 55, 6,                                                    &
     &6, 58,                                                            &
     &59, 60, 61,                                                       &
     &62, 5, 5, 5,                                                      &
     &6, 6, 6,                                                          &
     &69, 70,                                                           &
     &71, 72, 73, 74 /
      data (ivdep(j,2,2), j = 33, mnvar) /                              &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 40,                                                      &
     &40, 0,                                                            &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 55,                                                      &
     &55, 0,                                                            &
     &0, 0, 0,                                                          &
     &0, 0, 40, 40,                                                     &
     &0, 55, 55,                                                        &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      data (ivdep(j,1,3), j = 1, 32) /                                  &
     &1, 2, 3,                                                          &
     &4, 5, 6, 7, 8,                                                    &
     &9, 5, 6, 5, 6,                                                    &
     &14, 15, 16,                                                       &
     &17, 18, 19,                                                       &
     &20, 21, 22,                                                       &
     &24, 24, 25, 26, 27,                                               &
     &3, 3, 3,                                                          &
     &19,                                                               &
     &20 /
      data (ivdep(j,2,3), j = 1, 32) /                                  &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 0, 0, 40, 55,                                                  &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0,                                                          &
     &0, 0, 0, 0, 0,                                                    &
     &0, 24, 24,                                                        &
     &0,                                                                &
     &0 /
      data (ivdep(j,1,3), j = 33, mnvar) /                              &
     &33, 33,                                                           &
     &35, 36, 37,                                                       &
     &38, 39, 40, 5,                                                    &
     &5, 43,                                                            &
     &44, 45, 46,                                                       &
     &47,                                                               &
     &48, 48,                                                           &
     &50, 51, 52,                                                       &
     &53, 54, 55, 6,                                                    &
     &6, 58,                                                            &
     &59, 60, 61,                                                       &
     &62, 5, 5, 5,                                                      &
     &6, 6, 6,                                                          &
     &69, 70,                                                           &
     &71, 72, 73, 74 /
      data (ivdep(j,2,3), j = 33, mnvar) /                              &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 40,                                                      &
     &40, 0,                                                            &
     &0, 0, 0,                                                          &
     &0,                                                                &
     &0, 0,                                                             &
     &0, 0, 0,                                                          &
     &0, 0, 0, 55,                                                      &
     &55, 0,                                                            &
     &0, 0, 0,                                                          &
     &0, 0, 40, 40,                                                     &
     &0, 55, 55,                                                        &
     &0, 0,                                                             &
     &0, 0, 0, 0 /

      if (it .le. 0 .or. it .gt. 3)  then
        sovar(1) = svar
        sovar(2) = ' '
        sname    = svar
        ipflg(1) = 0
        ipflg(2) = 0
        goto 999
      endif
      sovar(1) = ' '
      sname = svar
!--- search in list of known variables
      do  iref = 1, mnvar
        if (svar .eq. svname(iref))  goto 9
      enddo
      call pupnbl(svar, k1, k2)
      do  iref = 1, mnvar
        call pupnbl(svname(iref), k1f, k2f)
        if (k2 + 1 .eq. k2f)  then
          if (index('xy', svname(iref)(k2f:k2f)) .ne. 0)  then
            if (svar(:k2) .eq. svname(iref)(:k2))  goto 9
          endif
        endif
      enddo
      goto 999
    9 continue
      if (iflag .eq. 0)  then
        sname = svname(iref)
        ipflg(1) = iproc(iref,it)
        ipflg(2) = intpo(iref)
        do  j = 1, mxdep
          if (ivdep(iref,j,it) .eq. 0)  then
            sovar(j) = ' '
          else
            sovar(j) = svname(ivdep(iref,j,it))
          endif
        enddo
      elseif (iflag .eq. 1) then
        sname = svlabl(iref)
        if (svar .ne. svname(iref))  then
!--- incomplete match
!    replace x or y in name by blank
          call pupnbl(sname, k1, k2)
          do  i = 2, k2
            if (index('XYxy', sname(i:i)) .ne. 0)  then
              sname(i:i) = ' '
            endif
          enddo
        endif
      elseif (iflag .eq. 2) then
        sname = svanno(iref)
      elseif (iflag .eq. 3) then
        if (svar .eq. svname(iref))  then
          sname = svname(iref)
        else
          sname = svname(iref)(:k2)
        endif
      else
        sname = svar
      endif
  999 end
      subroutine peiact(kact, np, x, y, ac, kf, kl)
!----------------------------------------------------------------------*
! Purpose:
!   Return first and last point of curve inside active window
! Input:
!   KACT        (int)   starting point for check
!   NP          (int)   number of points in XVAL, YVAL
!   X           (real)  x values
!   Y           (real)  y values
!   AC          (real)  active window in WC
! Output:
!   KF          (int)   first point inside, or 0
!   KL          (int)   last  point inside, or 0
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      real x(*), y(*), ac(4)
      real toleps, xtol, ytol
      parameter (toleps = 1.e-5)
      integer i, kact, np, kf, kl
      xtol = toleps * (ac(2) - ac(1))
      ytol = toleps * (ac(4) - ac(3))
      kf = 0
      kl = 0
      do i = kact, np
        if(x(i) + xtol .lt. ac(1)) goto 10
        if(x(i) - xtol .gt. ac(2)) goto 10
        if(y(i) + ytol .lt. ac(3)) goto 10
        if(y(i) - ytol .gt. ac(4)) goto 10
        kf = i
        goto 20
   10   continue
      enddo
!--- no point inside
      goto 999
   20 continue
      do i = kf, np
        if(x(i) + xtol .lt. ac(1)) goto 40
        if(x(i) - xtol .gt. ac(2)) goto 40
        if(y(i) + ytol .lt. ac(3)) goto 40
        if(y(i) - ytol .gt. ac(4)) goto 40
      enddo
   40 kl = i - 1
  999 continue
      end
      subroutine peintp(crow, nint, type, proc, telpar, ierr)

#ifdef _WIN32
      include 'win32calls.fi'
#endif


!----------------------------------------------------------------------*
! purpose:
!   interpolate variables plotted against s
! input:
!   crow        (int)   table row number at start of element
!   nint        (int)   number of interpolation intervals
!   type        (int)   (local) element type
!   proc        (int)   original process flags
!   step        (d.p.)  max. dist. between two successive hor. values
!   telpar      (d.p.)  temporary element parameters
! output:
!   ierr        (int)   0 if ok, else > 0
!   the results are stored in qhval and qvval
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      double precision pi, zero, eps, one, two, twopi
      parameter         (pi = 3.1415926535898d0)
      parameter         (zero = 0.d0, eps = 1.d-5)
      parameter         (one = 1.d0, two = 2.d0, twopi = two * pi)
      integer nint, crow, mat, type, ierr, proc(2,*)
      integer double_from_table
      double precision telpar(*)
      double precision get_value
      double precision tw0(mintpl), tw1(mintpl), dmu, step, am(6,6)
      double precision ex, ey, bx0, ax0, ux0, dx0, dpx0,                &
     &by0, ay0, uy0, dy0, dpy0, x0, px0, y0, py0, xn0, pxn0, yn0, pyn0
      double precision bx1, ax1, ux1, dx1, dpx1,                        &
     &by1, ay1, uy1, dy1, dpy1, x1, px1, y1, py1, xn1, pxn1, yn1, pyn1
      equivalence       (bx0, tw0(1)), (ax0, tw0(2)), (ux0, tw0(3)),    &
     &(dx0, tw0(4)), (dpx0, tw0(5)),                                    &
     &(by0, tw0(6)), (ay0, tw0(7)), (uy0, tw0(8)),                      &
     &(dy0, tw0(9)), (dpy0, tw0(10)), (x0,tw0(11)),                     &
     &(px0, tw0(12)), (y0, tw0(13)), (py0, tw0(14)),                    &
     &(xn0, tw0(15)), (pxn0, tw0(16)),                                  &
     &(yn0, tw0(17)), (pyn0, tw0(18))
      equivalence       (bx1, tw1(1)), (ax1, tw1(2)), (ux1, tw1(3)),    &
     &(dx1, tw1(4)), (dpx1, tw1(5)),                                    &
     &(by1, tw1(6)), (ay1, tw1(7)), (uy1, tw1(8)),                      &
     &(dy1, tw1(9)), (dpy1, tw1(10)), (x1,tw1(11)),                     &
     &(px1, tw1(12)), (y1, tw1(13)), (py1, tw1(14)),                    &
     &(xn1,tw1(15)), (pxn1, tw1(16)),                                   &
     &(yn1, tw1(17)), (pyn1, tw1(18))
      double precision s0, dstp, gamx, gamy
      integer i, j, k, ipc
      integer jrp,jcp
      logical elmflg
      double precision bxx0,byy0

      ierr = 0
      if (telpar(1) .eq. zero)  goto 999
      step = telpar(1) / nint
!---  set flag for correct interpolation
      elmflg = .false.
      do i = 1, nivvar
        elmflg = elmflg .or. proc_flag(2,i) .gt. 0
      enddo
      if (elmflg)  then
        k = double_from_table(tabname, 'x ', crow, x0)
        k = double_from_table(tabname, 'px ', crow, px0)
        k = double_from_table(tabname, 'betx ', crow, bx0)
        k = double_from_table(tabname, 'alfx ', crow, ax0)
        k = double_from_table(tabname, 'mux ', crow, ux0)
        k = double_from_table(tabname, 'dx ', crow, dx0)
        k = double_from_table(tabname, 'dpx ', crow, dpx0)
        k = double_from_table(tabname, 'y ', crow, y0)
        k = double_from_table(tabname, 'py ', crow, py0)
        k = double_from_table(tabname, 'bety ', crow, by0)
        k = double_from_table(tabname, 'alfy ', crow, ay0)
        k = double_from_table(tabname, 'muy ', crow, uy0)
        k = double_from_table(tabname, 'dy ', crow, dy0)
        k = double_from_table(tabname, 'dpy ', crow, dpy0)
        k = double_from_table(tabname, 's ', crow, s0)
        ex = get_value('beam ','ex ')
        ey = get_value('beam ','ey ')

!--- xn, pxn, yn, pyn
        if (ex * bx0 .eq. zero)  then
          xn0 = zero
        else
          xn0 = x0 / sqrt(ex*abs(bx0))
        endif
        if (ey * by0 .eq. zero)  then
          yn0 = zero
        else
          yn0 = y0 / sqrt(ey*abs(by0))
        endif
        if (bx0 .ne. zero)  then
          gamx = (one + ax0**2) / bx0
        else
          gamx = zero
        endif
        if (by0 .ne. zero)  then
          gamy = (one + ay0**2) / by0
        else
          gamy = zero
        endif
!--- get intermediate s values, and interpolate twiss parameters

        do i = 1, nint
          if (i .eq. 1)  then
            mat = 1
          elseif (i .eq. nint)  then
            mat = 3
          else
            mat = 2
          endif
          dstp = i * step
!--- element matrices 6x6: entry, body, exit
          call peelma(mat, type, telpar, dstp, am(1,1))

!--- interpolate twiss parameters
!    beta_x, beta_y
          bx1 = -two * am(1,1) * am(1,2) * ax0                          &
     &    + am(1,1)**2 * bx0 + am(1,2)**2 * gamx
          by1 = -two * am(3,3) * am(3,4) * ay0                          &
     &    + am(3,3)**2 * by0 + am(3,4)**2 * gamy

!--- alfa_x, alfa_y
          ax1 = (am(1,1) * am(2,2) + am(1,2) * am(2,1))                 &
     &    * ax0 - am(1,1) * am(2,1) * bx0                               &
     &    - am(1,2) * am(2,2) * gamx
          ay1 = (am(3,3) * am(4,4) + am(3,4) * am(4,3))                 &
     &    * ay0 - am(3,3) * am(4,3) * by0                               &
     &    - am(3,4) * am(4,4) * gamy
!--- mu_x, mu_y
          dmu = atan2(am(1,2), bx0 * am(1,1)                            &
     &    - ax0 * am(1,2))
          if (dmu .lt. zero)  dmu = dmu + twopi
          ux1 = ux0 + dmu / twopi
          dmu = atan2(am(3,4), by0 * am(3,3)                            &
     &    - ay0 * am(3,4))
          if (dmu .lt. zero)  dmu = dmu + twopi
          uy1 = uy0 + dmu / twopi
!--- d-x, d-y
          dx1 = am(1,1) * dx0 + am(1,2) * dpx0 + am(1,6)
          dy1 = am(3,3) * dy0 + am(3,4) * dpy0 + am(3,6)
!--- d'-x, d'-y
          dpx1 = am(2,1) * dx0 + am(2,2) * dpx0 + am(2,6)
          dpy1 = am(4,3) * dy0 + am(4,4) * dpy0 + am(4,6)
!--- x, px, y, py
          x1 = am(1,1) * x0 + am(1,2) * px0                             &
     &    + am(1,6) * currdp
          y1 = am(3,3) * y0 + am(3,4) * py0                             &
     &    + am(3,6) * currdp
          px1 = am(2,1) * x0 + am(2,2) * px0                            &
     &    + am(2,6) * currdp
          py1 = am(4,3) * y0 + am(4,4) * py0                            &
     &    + am(4,6) * currdp
!--- xn, pxn, yn, pyn
          if (ex * bx1.eq. zero)  then
            xn1 = zero
          else
            xn1 = x1 / sqrt(ex*abs(bx1))
          endif
          if (ey * by1 .eq. zero)  then
            yn1 = zero
          else
            yn1 = y1 / sqrt(ey*abs(by1))
          endif
          pxn1 = px1 * gamx
          pyn1 = py1 * gamy
!--- loop over variables, interpolate those with codes
          do j = 1, nivvar
            ipc = mod(proc(2,j), 100)
            if (ipc .gt. 0)  then
              if (nqval(j) .eq. maxseql)  then
                ierr = 1
                goto 999
              endif
              ipparm(2,j) = 1
              nqval(j) = nqval(j) + 1
              qhval(nqval(j),j) = s0 + dstp
              if (proc(1,j) .gt. 0)  then
                qvval(nqval(j), j) = sqrt(abs(tw1(ipc)))
              else
                qvval(nqval(j), j) = tw1(ipc)
              endif
            else
              ipparm(2,j) = 0
            endif
          enddo
        enddo
      endif
  999 end

      subroutine pelfill(tabname, type, row, n, a)

#ifdef _WIN32
      include 'win32calls.fi'
#endif


!----------------------------------------------------------------------*
! Purpose:
!   Returns element parameters needed for plotting
!--- Input:
!   tabname char  table name
!   type  int     MAD-8 element type
!   row   int     table row number
!   n     int     dimension of a
!--- I/O:
!    a     double  element parameters (as far as needed):
!  word / e_type = 1            = 2                               = 3
! C    F
! 0    1    l      [m]          l [m]                             l [m]
! 1    2    rhoinv [1/m]        volt [MV]                          deltap/p
! 2    3    e1                  ex [MV/m]                          kick
! 3    4    e2                  ey [MV/m]                           .
! 4    5    h1                  freq [MHz]                          .
! 5    6    h2                  lag [2 Pi]                          .
! 6    7    tilt                tilt                                .
! 7    8    ks                  betrf                              kick
! 8    9    hgap [m]            pg {MW]                 C:8-43 F:9-44: rm
! 9   10    fint [Tm]           shunt [MOhm/m]          C:44-259 F:45-260: tm
!10   11    angle = K_0*l       tfill [micro sec]
!11   12    lrad                harmon
!12   13    k0 or k0*l (l=0)    xsize (coll.) or xma (beam-beam) or x (mon.)
!13   14    k0s or k0s*l        ysize (coll.) or yma (beam_beam) or y (mon.)
!14   15    k1 or k1*l (l=0)    sigx
!15   16    k1s or k1s*l        sigy
!16   17    k2 or k2*l          fractional charge
!17   18    k2s  etc.           npart (# particles in opposite beam)
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer type,n
      double precision mi, mx, tmp1, tmp2, a(n)
      double precision node_value
      integer j, row, double_from_table
      character *(*) tabname
      do j = 1, n
        a(j) = 0.d0
      enddo
      tmp1 = 0.d0
      tmp2 = 0.d0
      a(1) = node_value('l ')
      if (a(1) .gt. 0.d0 .and. type .gt. 1 .and. type .lt. 8)  then
        if (n .ge. 11) then
!frs   use angle&tilt instead of k0l&k0sl since the later are no longer
!frs   calculated from the former 25.04.2003
!          j = double_from_table(tabname, 'k0l ' , row, tmp1)
!          j = double_from_table(tabname, 'k0sl ' , row, tmp2)
!          if (tmp1 .ne. 0.d0)  then
!            a(11) = tmp1
!          else
!            a(11) = tmp2
!          endif

          a(11) = node_value('angle ')
          a(7) = node_value('tilt ')
          a(3) = node_value('e1 ')
          a(4) = node_value('e2 ')
!          if (abs(tmp2) .ne. 0.d0) then
!            a(11) = cos(tmp2) * tmp1
!          else
!            a(11) = tmp1
!          endif
        endif
        if (n .ge. 15)                                                  &
     &  j = double_from_table(tabname, 'k1l ' , row, a(15))
        if (n .ge. 16)                                                  &
     &  j = double_from_table(tabname, 'k1sl ', row, a(16))
        if (n .ge. 17)                                                  &
     &  j = double_from_table(tabname, 'k2l ' , row, a(17))
        if (n .ge. 18)                                                  &
     &  j = double_from_table(tabname, 'k2sl ', row, a(18))
        if (n .ge. 19)                                                  &
     &  j = double_from_table(tabname, 'k3l ' , row, a(19))
        if (n .ge. 20)                                                  &
     &  j = double_from_table(tabname, 'k3sl ', row, a(20))
        do j = 15, n
          a(j) = a(j) / a(1)
        enddo
      endif
      mi = min(tmp2, a(16), a(18), a(20))
      mx = max(tmp2, a(16), a(18), a(20))
      if (mi .ne. mx) a(7) = 1.d0
      end
      subroutine pemima
!----------------------------------------------------------------------*
! Purpose:
!   Constrain axis reference, find minima and maxima of coordinates,
!   construct axis labels
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer i, j, k, iv, idum(2), ns, k1, k2, i1, i2, it(4)
      character * (mtitl)  s
      character * (mxlabl) slab
      character * (mcnam) sdum(mxcurv), saxis(mxcurv)

      do i = 1, 4
        it(i) = 0
      enddo
      numax = 0
      do  j = 1, nivvar
        do i = 1, numax
          if (it(i) .eq. naxref(j))  goto 10
        enddo
        if (numax .eq. 4)  then
          naxref(j) = it(4)
        else
          numax = numax + 1
          it(numax) = naxref(j)
        endif
   10   continue
      enddo
      do i = 1, 4
        do j = 1, numax - 1
          if (it(j) .gt. it(j+1))  then
            k = it(j)
            it(j) = it(j+1)
            it(j+1) = k
          endif
        enddo
      enddo
      do j = 1, nivvar
        do i = 1, numax
          if (naxref(j) .eq. it(i))  then
            naxref(j) = i
            goto 50
          endif
        enddo
   50   continue
      enddo
      do j = 1, nivvar
        k = naxref(j)
        do i = 1, nqval(j)
          hmima(1) = min(hmima(1), qhval(i,j))
          hmima(2) = max(hmima(2), qhval(i,j))
          vmima(1,k) = min(vmima(1,k), qvval(i,j))
          vmima(2,k) = max(vmima(2,k), qvval(i,j))
        enddo
      enddo
!--- get axis annotation
      do j = 1, nivvar
        k = naxref(j)
        nvvar(k) = nvvar(k) + 1
        vaxis(nvvar(k),k) = slabl(j)
      enddo
      do iv = 1, 4
        if (nvvar(iv) .gt. 0)  then
          if (nvvar(iv) .eq. 1)  then
            call pegetn (1, vaxis(1,iv), itbv, idum, sdum, slab)
            ns = 1
          else
            call pegaxn (nvvar(iv), vaxis(1,iv), saxis, ns)
            call pegetn (1, saxis(1), itbv, idum, sdum, slab)
          endif
          call gxpnbl (slab, k1, k2)
          s  = '<#>' // slab
          k2 = k2 + 3
          do i = 2, ns
            call pegetn (1, saxis(i), itbv, idum, sdum, slab)
            call gxpnbl (slab, i1, i2)
            if (index(s(:k2),slab(:i2)) .eq. 0)  then
              s(k2 + 1:) = ', ' // slab(:i2)
              k2 = k2 + i2 + 2
            endif
          enddo
          axlabel(iv) = s
        endif
      enddo
      end
      subroutine peplot

#ifdef _WIN32
      include 'win32calls.fi'
#endif


!----------------------------------------------------------------------*
! Purpose:
!   Plot all types of graphs from MAD.
!   Uses GXPLOT with underlying X-Windows (PostScript)
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer mlsize, mtsize, masize
      parameter            (mlsize = 13,mtsize = 13,masize = 20)
!--- character sizes:
!   MLSIZE    label character height
!   MTSIZE    text   - " -
!   MASIZE    annotation - " -
      character svar*(mcnam)
      character *(mxlabl) slocn, slname
      character stemp*(mtitl), stext*300, sform*20
      character sdum(mxdep)*(mcnam)
      real prmach, symch, tmpval
      double precision plot_option
      integer lastnb, iaxseq(4)
      integer idum, k1dum, k2dum, k3dum, i, npar, ivvar, nvax, ivax,    &
     &ierr, vdum(2)
!--- strings:
!   SVAR     buffer for variable names etc.
!   SLOCN    local name buffer (without leading "_")
!   STEMP    temporary buffer for titles
!   STEXT    buffer for labels etc.
!   SFORM    format buffer
!--- reals:
!   PRMACH fraction of viewport taken by machine plot
!   SYMCH  preset symbol character height

      data prmach /0.1/, symch /0.01/
      data iaxseq / 1, 4, 2, 3 /

!--- reset axis and curve defaults
      call gxsdef ('AXIS', 0)
      call gxsdef ('CURVE', 0)
!--- set "new line" character (change default = '/')
      call gxsvar ('SDEFNL', idum, fdum, '%')
!--- set top of viewport - leave space to plot machine if required
      if (fpmach)  then
        yvtop = 1. - prmach
      else
        yvtop = 1.
      endif
!--- set line width scale factor
      tmpval = plot_option('lwidth ')
      if (tmpval .eq. 0.) tmpval = 1.
      call jslwsc (tmpval)
!--- loop over frames
!--- set viewport
      vpt(1) = 0.
      vpt(2) = 1.
      vpt(3) = 0.
      vpt(4) = yvtop
      call gxsvpt (vpt)
!--- find variable name in list
      svar = horname
      call pegetn (1, svar, itbv, vdum, sdum, slname)
      slocn = ' '
      call gxpnbl(slname, k1dum, k2dum)
      k3dum = 0
      do idum = k1dum, k2dum
        if (slname(idum:idum) .ne. '_') then
          k3dum = k3dum + 1
          slocn(k3dum:k3dum) = slname(idum:idum)
        endif
      enddo
!--- prepare horizontal axis
      do i = 1, 4, 3
        call gxqaxs ('X', i, npar, ipar, range, stext, sform)
!--- set character sizes for labels and text including user requests
        ipar(7) = max (mlsize * qlscl + .01, 1.1)
        ipar(13) = max( mtsize * qtscl + .01, 1.1)
!--- text left adjusted
        ipar(10) = 1
!--- font
        ipar(11) = plot_option('font ')
        if (ipar(11) .eq. 0) ipar(11) = 1
!--- axis ref. number
        ipar(21) = 1
!--- range centre etc.
        if (hrange(1) .lt. hrange(2)) then
!--- use range as is
          ipar(23) = 1
          range(1) = hrange(1)
          range(2) = hrange(2)
!--- set min. and max. for horizontal axis
          xax(1) = hrange(1)
          xax(2) = hrange(2)
        else
          xax(1) = hmima(1)
          xax(2) = hmima(2)
        endif
        if (i .eq. 1) then
!--- bottom title
          stext = '<#>' // slocn(:lastnb(slocn))
        else
!--- suppress labels on upper axis
          ipar(3) = 0
!--- ticks below axis
          ipar(4) = 1
!--- top title
          stext = toptitle
        endif
!--- set axis parameters
        call gxsaxs ('X', i, npar, ipar, range, stext, sform)
      enddo
      do nvax = 1, numax
!--- set curve parameters for frame call
        ivax = iaxseq(nvax)
        call gxqcrv (nvax, npar, ipar, ssymb)
        ipar(2) = ivax
        call gxscrv (nvax, npar, ipar, ' ')
        call gxqaxs ('Y', ivax, npar, ipar, range, stext, sform)
!--- set character sizes for labels and text including user requests
        ipar(7) = max (mlsize * qlscl + .01, 1.1)
        ipar(13) = max (mtsize * qtscl + .01, 1.1)
!--- right adjusted label
        ipar(10) = 3
!--- font
        ipar(11) = plot_option('font ')
        if (ipar(11) .eq. 0) ipar(11) = 1
!--- range centre etc.
        if (vrange(1,nvax) .lt. vrange(2,nvax)) then
!--- use range as is
          ipar(23) = 1
          range(1) = vrange(1,nvax)
          range(2) = vrange(2,nvax)
!--- store y values for frame scaling
          yax(2 * nvax - 1) = vrange(1,nvax)
          yax(2 * nvax) = vrange(2,nvax)
        else
!--- store y values for frame scaling
          yax(2 * nvax - 1) = vmima(1,nvax)
          yax(2 * nvax) = vmima(2,nvax)
        endif
!--- get axis annotation
        slocn = axlabel(nvax)
        stemp = ' '
        call gxpnbl(slocn, k1dum, k2dum)
        k3dum = 0
        do idum = k1dum, k2dum
          if (slocn(idum:idum) .ne. '_') then
            k3dum = k3dum + 1
            stemp(k3dum:k3dum) = slocn(idum:idum)
          endif
        enddo
        if (nvax .eq. 1) then
          stext = '%' // stemp
        else
          stext = stemp
        endif
        call gxsaxs ('Y', ivax, npar, ipar, range, stext, sform)
        nptval(nvax) = 2
        ipxval(nvax) = 1
        ipyval(nvax) = 2 * nvax - 1
        icvref(nvax) = nvax
      enddo
!--- if only one y axis, plot right axis with ticks only
      if (numax .eq. 1) then
        ivax = 4
        call gxqaxs ('Y', ivax, npar, ipar, range, stext, sform)
        ipar(3) = 0
        ipar(4) = 1
        ipar(21) = 1
        call gxsaxs ('Y', ivax, npar, ipar, range, stext, sform)
      endif
!--- plot frame, keep windows for curves + clipping
      call gxfrm1 (numax, nptval, ipxval, ipyval, icvref, xax, yax,     &
     &window, actwin, ierr)
      if (ierr .ne. 0) goto 120
!--- now loop over vertical variables for real curve plotting
      do ivvar = 1, nivvar
        nvax = naxref(ivvar)
        ivax = iaxseq(nvax)
!--- find variable name in list for annotation
        svar = slabl(ivvar)
        call pegetn (2, svar, itbv, vdum, sdum, slname)
        slocn = ' '
        call gxpnbl(slname, k1dum, k2dum)
        k3dum = 0
        do idum = k1dum, k2dum
          if (slname(idum:idum) .ne. '_') then
            k3dum = k3dum + 1
            slocn(k3dum:k3dum) = slname(idum:idum)
          endif
        enddo
!--- character height including user request
        chh = 0.001 * masize * qascl
!--- call curve plot routine with simple arrays and flags
        call pecurv (ivvar, slocn, chh, qascl,                          &
     &  symch * qsscl, ipparm(1,ivvar), nqval(ivvar), qhval(1,ivvar),   &
     &  qvval(1,ivvar), window(1,nvax), actwin(1,nvax), ierr)
        if (ierr .ne. 0) goto 150
      enddo
      if (fpmach)  then
        vpt(1) = 0.
        vpt(2) = 1.
        vpt(3) = yvtop
        vpt(4) = 1.
        call gxsvpt (vpt)
        window(3,1) = -1.
        window(4,1) = 1.
        call gxswnd (window)
        call peschm (nelmach, ieltyp, xax, estart, eend, actwin)
      endif
!--- plot Arnold web if requested
      if (ngrid .ne. 0 .and. lgrid .eq. 0)  then
        call pearwe(igrid, actwin)
      endif
      goto 999

  120 continue
!--- curve for vert. var. missing
  150 continue
  999 continue
      end
      subroutine peqcon (kx, ky, kxys, ibk, iselct, ierr)
!----------------------------------------------------------------------*
! Purpose:
!   Check user constraints on Arnold web
!
!--- Input
!   KX        (integer) position of K_x in KXYS
!   KY        (integer) position of K_y in KXYS
!   KXYS      (integer) array containing K_x at KX, K_y at KY, KS at 3
!   IBK       (integer) array containing:
!                       number of superperiods
!                       number of constraints  N
!                       N constraints:
!                       minimum
!                       maximum
!                       step
!                       Length L
!                       expression in polish notation of length L,
!                       coded as 1+, 2-, 3*, 4/, 1 KX, 2 KY, 3 KS
!                       + MQADD + 4, e.g. 100006 = KY
!--- Output
!   ISELCT    (integer) 1 if combination of KX, KY, KS accepted,
!                       0 if not accepted
!   IERR                0 if OK, 1 if expression illegal, 2 if / 0
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer ibk(*), kxys(3)
      integer kx, ky, iselct, ierr, ipt, mqadd4, ic, low, lup, lsp, nw, &
     &nstack, i, kc, k, kop
      integer kloc(3), istack(100)
      ierr = 0
      kloc(1) = kx
      kloc(2) = ky
      kloc(3) = 3
      iselct = 1
      ipt = 2
      mqadd4 = mqadd + 4
      do ic = 1, ibk(2)
        low = ibk(ipt+1)
        lup = ibk(ipt+2)
        lsp = max(1, ibk(ipt+3))
        nw  = ibk(ipt+4)
        ipt = ipt + 4
!--- calculate expression given in inverse Polish
        nstack = 0
        do i = 1, nw
          kc = ibk(ipt+i)
          if (kc .le. mqadd)  then
!--- simple integer
            nstack = nstack + 1
            istack(nstack) = kc
          elseif (kc .gt. mqadd4)  then
!--- K_x, K_y, or K_s
            k = kc - mqadd4
            nstack = nstack + 1
            istack(nstack) = kxys(kloc(k))
          else
!--- operator
            if (nstack .lt. 2)  then
              ierr = 1
              goto 999
            endif
            kop = kc - mqadd
            if (kop .eq. 1)  then
              istack(nstack-1) =                                        &
     &        istack(nstack-1) + istack(nstack)
            elseif (kop .eq. 2)  then
              istack(nstack-1) =                                        &
     &        istack(nstack-1) - istack(nstack)
            elseif (kop .eq. 3)  then
              istack(nstack-1) =                                        &
     &        istack(nstack-1) * istack(nstack)
            elseif (kop .eq. 4)  then
              if (istack(nstack) .eq. 0)  then
                ierr = 2
                goto 999
              endif
              istack(nstack-1) =                                        &
     &        istack(nstack-1) / istack(nstack)
            endif
            nstack = nstack - 1
          endif
        enddo
        if (nstack .ne. 1)  then
          ierr = 1
          goto 999
        endif
        do i = low, lup, lsp
          if (istack(1) .eq. i) goto 30
        enddo
!--- test failed
        iselct = 0
        goto 999
   30   ipt = ipt + nw
      enddo
  999 end
      subroutine  peschm (nel, ityp, hr, es, ee, actwin)
!----------------------------------------------------------------------*
! Purpose:
!   Plot schema
! Input:
!   nel      (integer)  no. of elements
!   ityp     (integer)  array with element types:
!                       0: drift                                       *
!                       1: sbend, zero tilt                            *
!                       2: focussing quad                              *
!                       3: defocussing quad                            *
!                       4: monitor                                     *
!                       5: collimator                                  *
!                       6: electrostatic separator                     *
!                       7: sbend, non-zero tilt                        *
!                       8: multipole                                   *
!                       9: RF cavity                                   *
!                       10: positive sext                              *
!                       11: negative sext                              *
!                       12: positive oct                               *
!                       13: negative oct                               *
!                       14: lcavity                                    *
!                       21: rbend, zero tilt                           *
!                       27: rbend, non-zero tilt                       *
!   hr          (real)  horizontal range (lower and upper)
!   es          (real)  array with element start position
!   ee          (real)  array with element end position
!   actwin      (real)  active window for curve plot (array of 4)
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer ityp(*)
      integer nel, i, it, j, mobj, msize
      parameter (mobj = 14, msize = 88)
      integer npst(mobj), npnd(mobj), npsl(msize)
      real ell, hr(2), es(*), ee(*), actwin(4)
      real shapex(msize), shapey(msize)
      real txp(2), typ(2), typz(2)

      data npst   / 1,  6, 11, 16, 21,                                  &
     &33, 43, 48,                                                       &
     &50,                                                               &
     &64, 69, 74, 79, 84 /
      data npnd   / 5, 10, 15, 20, 32,                                  &
     &42, 47, 49,                                                       &
     &63,                                                               &
     &68, 73, 78, 83, 88 /
      data npsl   /5 * 1, 5 * 1, 5 * 1, 5 * 3, 5 * 1, 0, 4 * 1, 0, 1,   &
     &1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 5 * 1, 2 * 1,                       &
     &6 * 1, 0, 5 * 1, 0, 1,                                            &
     &5 * 1, 5 * 1, 5 * 1, 5 * 1, 5 * 1 /
      data typz   / 2 * 0. /
      data shapex /0., 1., 1., 0., 0.,                                  &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0., 0., 1., 1., 0., 0., 0., 1.,                   &
     &0., 1., 0.5, 0.5, 0., 1., 0.5, 0.5, 0., 1.,                       &
     &0., 1., 1., 0., 0.,                                               &
     &0., 0.,                                                           &
     &0., 0.25, 0.25, 0.75, 0.75, 1.,                                   &
     &0., 0.25, 0.25, 0.75, 0.75, 1., 0., 1.,                           &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0.,                                               &
     &0., 1., 1., 0., 0. /
      data shapey /0.6, 0.6, -0.6, -0.6, 0.6,                           &
     &0., 0., 0.8, 0.8, 0.,                                             &
     &0., 0., -0.8, -0.8, 0.,                                           &
     &0.6, 0.6, -0.6, -0.6, 0.6,                                        &
     &0.8, 0.8, 0.4, 0.4, 0.8, -0.8, -0.8, -0.4, -0.4, -0.8, 0., 0.,    &
     &0.4, 0.4, 0.8, 0.4, -0.4, -0.4, -0.8, -0.4, 0., 0.,               &
     &0.5, 0.5, -0.5, -0.5, 0.5,                                        &
     &0.5, -0.5,                                                        &
     &0.2, 0.2, 0.8, 0.8, 0.2, 0.2,                                     &
     &-0.2, -0.2, -0.8, -0.8, -0.2, -0.2, 0., 0.,                       &
     &0., 0., 0.5, 0.5, 0.,                                             &
     &0., 0., -0.5, -0.5, 0.,                                           &
     &0., 0., 0.25, 0.25, 0.,                                           &
     &0., 0., -0.25, -0.25, 0.,                                         &
     &0.2, 0.2, -0.2, -0.2, 0.2 /

!--- set line style to solid
      call jsln(1)
      do i = 1, nel
        it = mod(ityp(i), 20)
        ell = ee(i) - es(i)
        if (i .eq. 1) then
          if(es(1) .gt. hr(1))  then
            txp(1) = hr(1)
            txp(2) = es(1)
            call gvpl (2, txp, typz)
          endif
        else
          if (ee(i-1) .lt. es(i))  then
            txp(1) = ee(i-1)
            txp(2) = es(i)
            call gvpl (2, txp, typz)
          endif
        endif
        if (es(i) .gt. actwin(2)) goto 50
        if (ee(i) .ge. actwin(1)) then
          txp(1) = es(i) + shapex(npst(it)) * ell
          typ(1) = shapey(npst(it))
          do  j = npst(it)+1, npnd(it)
            txp(2) = es(i) + shapex(j) * ell
            typ(2) = shapey(j)
            if (npsl(j) .gt. 0)  then
              call jsln(npsl(j))
              call gvpl(2, txp, typ)
            endif
            txp(1) = txp(2)
            typ(1) = typ(2)
          enddo
        endif
      enddo
   50 continue
      call jsln(1)
      if (ee(nel) .lt. hr(2))  then
        txp(1) = ee(nel)
        txp(2) = hr(2)
        call gvpl (2, txp, typz)
      endif
      end
      subroutine pesopt(ierr)

#ifdef _WIN32
      include 'win32calls.fi'
#endif

!----------------------------------------------------------------------*
! Purpose:
!   Stores plot options and values, checks
!
! Output:  ierr  (int)     =0: OK, >0: error
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,noversi,                   &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
!   ITYP     (integer)  array with element types:
!                       0: drift
!                       1: bend, zero tilt
!                       2: focussing quad
!                       3: defocussing quad
!                       4: monitor
!                       5: collimator
!                       6: electrostatic separator
!                       7: bend, non-zero tilt
!                       8: multipole
      integer ierr, table_org
      integer i, j, k
      character * (mcnam) sdum(mxcurv)
      integer nint, ndble, int_arr(100), char_l(100)
      double precision d_arr(100)
      double precision plot_option
      character * 400 char_a, version
      character * 10 tmp_a
      ierr = 0
      nivaxs = 0
      nivvar = 0
      nelmach = 0
      fpmach = .false.
      ndeltap = 1
      toptitle = ' '
      hmima(1) = 1.e20
      hmima(2) = -hmima(1)
      do i = 1, 4
        vmima(1,i) = 1.e20
        vmima(2,i) = -1.e20
        nvvar(i) = 0
      enddo
      char_a = ' '
      tmp_a = 'notitle '
      notitle = 0
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (nint .gt. 0) notitle = int_arr(1)

!--- get noversi
      char_a = ' '
      tmp_a = 'noversion '
      noversi = 0
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (nint .gt. 0) noversi = int_arr(1)

!--- any table - for hor = s plot machine
      char_a = ' '
      tabname = ' '
      tmp_a = 'table '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (k .gt. 0)  tabname = char_a

      char_a = ' '
      tmp_a = 'haxis '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (k .eq. 0)  then
        print *, 'no horizontal variable'
        ierr = 1
        goto 999
      else
        horname = char_a
      endif
      if (horname .eq. 's' .and. table_org(tabname) .eq. 0)  then
        itbv = 1
      else
        itbv = 0
      endif

!--- Prepare title
      if (notitle .eq. 0)  then
        char_a = ' '
        tmp_a = 'title '
        call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,           &
     &  char_a, char_l)
        if (k .eq. 0) then
          call get_title(char_a, k)
        else
          k = char_l(1)
        endif
        if (noversi .eq. 0) then
          call get_version(version, j)
          if (k .gt. 0)  then
            toptitle = char_a(:k) // '<#>' // version(:j)
          else
            toptitle = '<#>' // version(:j)
          endif
        else
          if (k .gt. 0) then
            toptitle = char_a(:k)
          else
            toptitle = ' '
          endif
        endif
      endif

      qascl = plot_option('ascale ')
      qlscl = plot_option('lscale ')
      qsscl = plot_option('sscale ')
      qtscl = plot_option('rscale ')

      char_a = ' '
      tmp_a = 'range '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      call table_range(tabname, char_a, nrrang)
      if (nrrang(1) .eq. 0 .and. nrrang(2) .eq. 0)  then
        print *, 'unknown table or illegal range, skipped'
        ierr = 1
        goto 999
      endif

      char_a = ' '
      tmp_a = 'noline '
      noline = 0
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (nint .gt. 0) noline = int_arr(1)

      char_a = ' '
      tmp_a = 'qs '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      qsval = d_arr(1)

      hrange(1) = 0
      tmp_a = 'hmin '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (ndble .gt. 0) hrange(1) = d_arr(1)

      hrange(2) = 0
      tmp_a = 'hmax '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (ndble .gt. 0) hrange(2) = d_arr(1)

      tmp_a = 'vmin '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      do i = 1, 4
        vrange(1,i) = 0
      enddo
      do i = 1, ndble
        vrange(1,i) = d_arr(i)
      enddo

      tmp_a = 'vmax '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      do i = 1, 4
        vrange(2,i) = 0
      enddo
      do i = 1, ndble
        vrange(2,i) = d_arr(i)
      enddo

      do i = 1, 5
        ipparm(i,1) = 0
      enddo
      char_a = ' '
      tmp_a = 'style '
      call comm_para(tmp_a, nint, ndble, k, ipparm(1,1), d_arr,         &
     &char_a, char_l)
      tmp_a = 'spline '
      call comm_para(tmp_a, nint, ndble, k, ipparm(2,1), d_arr,         &
     &char_a, char_l)
      tmp_a = 'bars '
      call comm_para(tmp_a, nint, ndble, k, ipparm(3,1), d_arr,         &
     &char_a, char_l)
      tmp_a = 'symbol '
      call comm_para(tmp_a, nint, ndble, k, ipparm(4,1), d_arr,         &
     &char_a, char_l)
      tmp_a = 'colour '
      call comm_para(tmp_a, nint, ndble, k, ipparm(5,1), d_arr,         &
     &char_a, char_l)
      do i = 2, mxcurv
        ipparm(1,i) = ipparm(1,1)
        ipparm(2,i) = ipparm(2,1)
        ipparm(3,i) = ipparm(3,1)
        ipparm(4,i) = ipparm(4,1)
        ipparm(5,i) = ipparm(5,1)
      enddo
      splinef = ipparm(2,1)

      char_a = ' '
      tmp_a = 'vaxis '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (k .gt. 0)  then
        nivaxs = 1
        nivvar = min(k, mxcurv)
        call pesplit(k, char_a, char_l, slabl)
        do j = 1, nivvar
          naxref(j) = 1
        enddo
      else
        char_a = ' '
        tmp_a = 'vaxis1 '
        call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,           &
     &  char_a, char_l)
        if (k .gt. 0)  then
          if (nivvar+k .gt. mxcurv) goto 100
          nivaxs = nivaxs + 1
          call pesplit(k, char_a, char_l, slabl(nivvar+1))
          do j = 1, k
            nivvar = nivvar + 1
            naxref(nivvar) = 1
          enddo
        endif
        char_a = ' '
        tmp_a = 'vaxis2 '
        call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,           &
     &  char_a, char_l)
        if (k .gt. 0)  then
          if (nivvar+k .gt. mxcurv) goto 100
          nivaxs = nivaxs + 1
          call pesplit(k, char_a, char_l, slabl(nivvar+1))
          do j = 1, k
            nivvar = nivvar + 1
            naxref(nivvar) = 2
          enddo
        endif
        char_a = ' '
        tmp_a = 'vaxis3 '
        call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,           &
     &  char_a, char_l)
        if (k .gt. 0)  then
          if (nivvar+k .gt. mxcurv) goto 100
          nivaxs = nivaxs + 1
          call pesplit(k, char_a, char_l, slabl(nivvar+1))
          do j = 1, k
            nivvar = nivvar + 1
            naxref(nivvar) = 3
          enddo
        endif
        char_a = ' '
        tmp_a = 'vaxis4 '
        call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,           &
     &  char_a, char_l)
        if (k .gt. 0)  then
          if (nivvar+k .gt. mxcurv) goto 100
          nivaxs = nivaxs + 1
          call pesplit(k, char_a, char_l, slabl(nivvar+1))
          do j = 1, k
            nivvar = nivvar + 1
            naxref(nivvar) = 4
          enddo
        endif
      endif
      if (nivvar .eq. 0)  then
        print *, 'Warning: no vertical plot variables, plot skipped'
        ierr = 1
        goto 999
      endif
      goto 110
  100 continue
      print *, 'Warning: # vertical variables cut at = ', nivvar
  110 continue
      do j = 1, nivvar
        call pegetn (0, slabl(j), itbv, proc_flag(1,j), sname(j),       &
     &  sdum(1))
        if (slabl(j) .eq. 'rbetx')  then
          sname(j) = 'betx'
          proc_flag(1,j) = 1
        else if (slabl(j) .eq. 'rbety')  then
          sname(j) = 'bety'
          proc_flag(1,j) = 1
        else
          sname(j) = slabl(j)
          proc_flag(1,j) = 0
        endif
      enddo
  999 end
      subroutine pesplit(n_str, char_a, char_l, char_buff)
      implicit none
      integer n_str, char_l(*)
      character*(*) char_a, char_buff(*)
      integer i, k, l
      k = 0
      do i = 1, n_str
        l = char_l(i)
        char_buff(i) = char_a(k+1:k+l)
        k = k+l
      enddo
      end
      subroutine plginit

#ifdef _WIN32
      include 'win32calls.fi'
#endif

!----------------------------------------------------------------------*
! Purpose:
!   Overall initialization
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      logical intrac
      integer ipseps, iset, nint, ndble, k, int_arr(100), char_l(100)
      double precision plot_option
      double precision d_arr(100)
      real tmpval
      character * 40 char_a
      character * 8 tmp_a
      data iset / 0 /
      call gxtint
      call gxsvar ('INUNIT', 5, 0., ' ')
      call gxsvar ('IOUNIT', 6, 0., ' ')
      char_a = ' '
      tmp_a = 'file '
      call comm_para(tmp_a, nint, ndble, k, int_arr, d_arr,             &
     &char_a, char_l)
      if (k .gt. 0) then
        plfnam = char_a(:char_l(1))
      else
        plfnam = 'madx'
      endif
      ipseps = plot_option('post ')
      if (ipseps .eq. 0 .and. .not. intrac())  then
        ipseps = 2
      endif
      if (iset .eq. 0 .and. ipseps .ne. 0) then
        iset = 1
        call gxsvar ('SMETNM', 0, 0., plfnam)
        call gxsvar('IPSEPS', ipseps, 0., ' ')
      endif
      if (intrac())  then
!--- set wait time to 1 sec.
        call gxsvar ('WTTIME', 0, 1., ' ')
        call gxasku
      endif
!--- reduce window size (only X11)
      call gxsvar('NYPIX', 670, 0., ' ')
!--- set bounding box (only X11)
      tmpval=plot_option('xsize ')
      call gxsvar('XMETAF', 0, tmpval, ' ')
      tmpval=plot_option('ysize ')
      call gxsvar('YMETAF', 0, tmpval, ' ')
!--- inhibit initial X-Window (only X11)
      call gxsvar('ITSEOP', 1, 0., ' ')
      call gxinit
      call gxclos
      end
      subroutine plotit(initfl)
!----------------------------------------------------------------------*
! Purpose:
!   Plots on screen and/or file
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax,       &
     &mplred, mplout
      parameter (mtitl  = 128, mxlabl = 80)
      parameter (mnvar = 74, mxdep = 2)
      parameter (mqadd = 100000)
      parameter (mntmax = 20, mksmax = 10)
      parameter (mplred = 46, mplout = 47)

      common /peotcl/ fpmach
      save   /peotcl/
      logical fpmach
      common /peotcr/ yvtop, fdum, chh,                                 &
     &vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)
      save   /peotcr/

      real yvtop, fdum, chh
      real vpt, window, actwin, range, xax, yax

      common /peotci/ ivnarw, ipar(50), nptval(4), ipxval(4),           &
     &ipyval(4), icvref(4)
      integer ivnarw, ipar, nptval, ipxval, ipyval, icvref
      save   /peotci/
      integer maxseql, mtwcol, mpparm,                                  &
     &mxcurv, mopt, mfile, marg, maxarg,                                &
     &mxdp, mxplot
      parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
     &mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
     &mxdp = 25, mxplot = 100)
      integer mintpl
      parameter (mintpl = 18)
      common / peaddi / itbv, ntmax, mared, nivvar, nivaxs,             &
     &nelmach, ngrid, lgrid,        ndeltap, ndeltas, numax,            &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag,            &
     &n_names, splinef, noline, notitle,                                &
     &nqval(mxcurv), nvvar(4), nrrang(2), proc_flag(2, mxcurv),         &
     &lvcol(maxseql), ipparm(mpparm,mxcurv), naxref(mxcurv),            &
     &ieltyp(maxseql), igrid(maxseql), intopt(mxdp,mopt),               &
     &numopt(mopt), countopt(mopt), pstart(mxplot), pend(mxplot),       &
     &occnt(maxseql), twtpos(mintpl)
      integer itbv, ntmax, mared, nivvar, nivaxs, nelmach,              &
     &nqval, ngrid, lgrid,        ndeltap, ndeltas, numax,              &
     &nplot, currplot, suml, nrow, ncolumn, nvcol, pos_flag, lvcol,     &
     &n_names, splinef, nvvar, nrrang, proc_flag, ipparm, naxref,       &
     &ieltyp, igrid, intopt, noline, notitle,                           &
     &numopt, countopt, pstart, pend,                                   &
     &occnt, twtpos
      common / peaddr / qascl, qlscl, qsscl, qtscl, qsval, currdp,      &
     &hrange(2), vrange(2,4), hmima(2), vmima(2,4),                     &
     &qhval(maxseql,mxcurv), qvval(maxseql,mxcurv),                     &
     &estart(maxseql), eend(maxseql)
      real qascl, qlscl, qsscl, qtscl, qsval, currdp,                   &
     &hrange, vrange, hmima, vmima,                                     &
     &qhval, qvval,                                                     &
     &estart, eend
      common / peaddc / horname, tabname, toptitle,                     &
     &inname, outname, dbname, sequ_name, tabtype, plfnam, plpnam,      &
     &arglist(maxarg), elnames(maxseql),                                &
     &axlabel(4), sname(mxcurv), slabl(mxcurv),                         &
     &vaxis(mxcurv,4), charopt(mxcurv,mopt),                            &
     &defopt(mopt), c_names(mtwcol), ssymb
      character * (mfile) plfnam, plpnam
      character * (mcnam) horname, tabname, sname, slabl,               &
     &vaxis, defopt, elnames, sequ_name, c_names, tabtype
      character * (marg) arglist
      character * (mfile) inname, outname, dbname
      character * (mxlabl) axlabel
      character * (mtitl) toptitle, charopt
      character * 1 ssymb
      save /peaddi/, /peaddr/, /peaddc/
      integer initfl
      if (initfl .eq. 0)  then
!--- overall initalization
        call plginit
        plpnam = plfnam
      endif
      if (plpnam .ne. plfnam)  then
        call gxsvar ('SMETNM', 0, 0., plfnam)
!--- close current .ps file if any
        call gxterm
        plpnam = plfnam
        call gxinit
      endif
      call gxopen
      call peplot
      call gxwait
      call gxclrw
      call gxclos
      end
      subroutine pupnbl(string,ifirst,ilast)
!***********************************************************************
!
!   Purpose: returns position of first and last non-blank in STRING
!
!--- Input
!   string     character string
!--- Output
!   ifirst     first non-blank in string, or 0 if only blanks
!   ilast      last non-blank
!
!   Author: H. Grote / CERN                        date: June 16, 1987
!                                             last mod: Sept. 13, 2001
!
!***********************************************************************
      character *(*)  string
      ifirst=0
      ilast=0
      do i=1,len(string)
        if(string(i:i).ne.' ') then
          ifirst=i
          goto 20
        endif
      enddo
      goto 999
   20 continue
      do i=len(string),1,-1
        if(string(i:i).ne.' ') then
          ilast=i
          goto 999
        endif
      enddo
  999 end
      integer function iucomp(comp, arr, n)
!----------------------------------------------------------------------*
! Purpose:
!   Find first occurrence of integer in integer array
!
!---Input:
!   comp        integer being looked up
!   arr         integer array being searched
!   n           length of arr
!
!  returns 0 if not found, else position in arr
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      integer comp, arr(*), n
      integer j
      iucomp = 0
      do j = 1, n
        if (comp .eq. arr(j))  then
          iucomp = j
          return
        endif
      enddo
      end
