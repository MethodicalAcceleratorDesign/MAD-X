subroutine plotit(initfl)
  use plotfi
  use plot_bfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Plots on screen and/or file
  !
  ! calls the routines gxsvar, gxterm, gxinit, gxopen, gxwait, gxclrw and
  ! gxclos in file gxx11.F.
  ! calls  the routines plginit and peplot in this file.
  !
  ! it is called by the routine exec_plot in file madxn.c after pemima.
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  integer :: initfl

  !--- type definitions of local variables
  character(len=mfile),save :: plpnam
  integer, save :: plot_No

  !--- Routine body

  if (initfl .eq. 0)  then
     !--- overall initialization
     plot_No = 0
     call plginit
     plpnam = plfnam
  endif

  plot_No = plot_No + 1
  print *,"plot number = ",plot_No

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
end subroutine plotit

subroutine pecurv (ncc, spname, annh, usex, sych, ippar, np, xval, yval, window, actwin, ierr)
  use plotfi
  implicit none
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
  !
  ! calls the routines peiact and pegacn in this file.
  ! calls the routines gxsave, gxswnd, gxpnbl, jsln, jsplci,  gxpl,
  !                    gxplt1, gvpl, jsmk, gxpmsw, jspmci, jschh, jstxal,
  !                    gxqrvp, jschxp, gxtx, jstxci, gxrest
  !                    defined in file gxx11.F.
  !
  ! it is called by the routine peplot in this file.
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  integer  :: ncc, np, ierr, ippar(*)
  character(len=*) :: spname
  real :: annh, sych, usex, xval(*), yval(*), window(*), actwin(*)

  !--- type definitions of local variables
  real :: xpos, ypos, xwpos, ywpos, xf, wsclx, wscly, xaux, yaux, xmd
  real :: rsave(20), xpl(2), ypl(2), winnor(4)
  real, save :: act(4)
  real :: xreal(maxpnt), yreal(maxpnt)
  integer :: isymb, icolr, ilb, ispli, kft, klt, kact, kf, kl, npt, j, iecub
  integer :: istyl, ipbar, isave(20)
  character(len=1) :: symloc
  logical :: zero_suppr, marker_plot

  double precision, external :: get_value

  !--- Initialisation of local variables
  data winnor /0., 1., 0., 1./

  iecub = 0
  zero_suppr = get_value('plot ','zero_suppr ') .ne. 0.d0
  marker_plot = get_value('plot ','marker_plot ') .ne. 0.d0

  !--- Output initialisation
  ierr = 0

  !--- Routine body

  !--- save GKS settings
  call gxsave (isave, rsave, ierr)
  call gxswnd (window)
  wsclx = 1. / (window(2) - window(1))
  wscly = 1. / (window(4) - window(3))
  xmd = 1.e-8 * (window(2) - window(1))**2

  if (ncc .eq. 1)  then
     !--- first curve in frame - reset label position array, get act. window in NDC
     act(1) = (actwin(1) - window(1)) * wsclx
     act(2) = (actwin(2) - window(1)) * wsclx
     act(3) = (actwin(3) - window(3)) * wscly
     act(4) = (actwin(4) - window(3)) * wscly
  endif

  istyl = ippar(1)
  ipbar = ippar(3)
  isymb = ippar(4)
  icolr = ippar(5)

  if (icolr .eq. 100)  icolr = mod(ncc-1,4) + 1
  icolr = max(1, min(icolr, 7))
  ilb = -1

  if (istyl .ne. 0)  then
     !--- polyline requested
     if (np .lt. 2) return

     if (istyl .eq. 100)  istyl = mod(ncc-1,4) + 1
     ispli = ippar(2)

     !--- get first and last blank in annotation
     call gxpnbl (spname, kft, klt)
     if (kft .ne. 0) ilb = 0 !--- annotation exists

     !--- set line style
     call jsln (max (1, min (4, istyl)))

     !--- set line colour
     call jsplci(icolr)
     kact = 1

10   continue
     !--- get first and last point inside
     call peiact(kact, np, xval, yval, actwin, kf, kl)
     if (kf .eq. 0) goto 40      !--- break if no points inside

     kf = max(1,  kf - 1)
     kl = min(np, kl + 1)
     npt      = 1
     xreal(1) = xval(kf)
     yreal(1) = yval(kf)
     do j = kf + 1, kl
        !--- avoid identical points
        if ( (marker_plot .and. .not.zero_suppr) .or.  &
             (xreal(npt) - xval(j))**2 + (yreal(npt) - yval(j))**2 .gt. xmd   .and. &
             (yval(j).ne. 0 .or. .not.zero_suppr) ) then ! hbu optionally skip 0 points
           npt        = npt + 1
           xreal(npt) = xval(j)
           yreal(npt) = yval(j)
        endif

        if ((j.eq.kl .and. npt.ge.2) .or. npt.eq.maxpnt) then
           !--- plot - get first curve annotation position
           if (ilb .eq. 0) &
                call pegacn(ncc, window, act, xreal, yreal, npt, usex, xwpos, xpos, ypos, ilb)

           if (interf.eq.0 .or. npt.eq.2 .or. ispli.ne.0) then
              !--- no spline
              call gxpl (npt, xreal, yreal, actwin)
              if (ilb .gt. 0) then
                 !--- get y pos. on curve for label
                 ywpos = yreal(ilb-1) + &
                      (yreal(ilb) - yreal(ilb-1)) * (xwpos - xreal(ilb-1)) &
                       / (xreal(ilb) - xreal(ilb-1))
                 ilb = -2
              endif
           else
              !--- spline
              call gxplt1 (npt, xreal, yreal, actwin)
              if (ilb .gt. 0) ilb   = -2
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
     if (np .eq. 0) return
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

40 continue ! jump from polyline and no point inside

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
           if (xaux.gt.act(1) .and. xaux.lt.act(2) .and. &
               yaux.gt.act(3) .and. yaux.lt.act(4)) &
               call gxtx (xaux, yaux, symloc)
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
     if (ypl(2) .gt. ypl(1))  ypl(1) = ypl(1) + .02d0
     !--- set dotted line
     call jsln (3)
     !--- set line colour
     call jsplci(icolr)
     !--- plot line
     call gxpl (2, xpl, ypl, act)
  endif

  !--- restore
  call gxrest (isave, rsave)

end subroutine pecurv

subroutine pefill(ierr)
  use plotfi
  use plot_bfi
  use math_constfi, only : zero, half, ten
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   fill plot arrays with coordinate values, set up machine plot
  !
  ! Output:  ierr  (int)     =0: OK, >0: error
  !
  ! calls the functions double_from_table_row, advance_to_pos
  ! table_length and restart_sequ defined in file madxn.c.
  ! calls the function lastnb defined in file util.F.
  ! calls the routines peintp in this file.
  !
  ! it is called by the routine exec_plot in file madxn.c after pesopt.
  !----------------------------------------------------------------------*
  integer, intent(OUT) :: ierr

  integer :: i, j, k, l, nint, new1, new2, currtyp, crow, pltyp, pos_flag
  integer :: ilist(mtype), p(mxplot), proc_n(2, mxcurv)
  double precision :: fact, d_val, d_val1
  double precision :: currpos, currleng, currtilt
  double precision :: currk1l, currk1sl, currk2l, currk2sl, currk3l, currk3sl
  real :: tval, step, mystep
  logical :: machp, rselect, marker_plot, range_plot
  character (len=120) :: warnstr

  !--- definitions of function primitives
  integer, external ::  double_from_table_row, restart_sequ, advance_to_pos
  integer, external ::  lastnb, table_length, advance_node, get_option
  double precision, external ::  node_value, get_value

  !--- codes see in peschm
  ! plot types (pltyp) as a function of mad code
  data ilist / 0, 21, 1, 0, 2, 10, 12, 8, 0, 9,  &
               6,  0, 0, 0, 0,  0,  4, 4, 4, 0,  &
               0,  0, 0, 0, 0,  0, 14, 0, 0, 0,  &
               0,  0, 0, 0, 0,  0,  9, 0, 0, 0,  &
               10 * 0 /

  !--- Initialize return code
  ierr = 0

  !--- Initialize marker_plot logical
  marker_plot = get_value('plot ','marker_plot ') .ne. zero
  range_plot  = get_value('plot ','range_plot ') .ne. zero

  !--- Initialisation of variables in common peaddi
  nelmach = 0
  NQVAL(1:mxcurv) = 0
  IELTYP(1:maxseql) = 0

  !--- Initialisation of variables in common peaddr
  ESTART(1:maxseql) = 0.0
  EEND(1:maxseql)   = 0.0
  QHVAL(1:maxseql,1:mxcurv) = 0.0
  QVVAL(1:maxseql,1:mxcurv) = 0.0

  !--- Initialisation of variables in common peotcl
  dpp_flag = .false.

  !--- Initialisation of local variables
  currtilt = zero
  currk1l = zero;  currk1sl = zero
  currk2l = zero;  currk2sl = zero
  currk3l = zero;  currk3sl = zero

  new1 = 0;  new2 = 0
  currtyp = 0
  nint = 0
  step = zero
  crow = 0
  pltyp = 0
  pos_flag = 3
  machp = itbv .ne. 0
  P(1:mxplot) = 0

  !--- save process flags, proc_flag may be modified in peintp
  PROC_N(1:2,1:mxcurv) = PROC_FLAG(1:2,1:mxcurv)

  !--- Routine body

  !--- No interpolation if centre option is set
  if (get_option('centre ') .ne. 0) then
     if (interf .eq. 1) &
        call fort_warn('PLOT: ','Interpolation is not compatible with the Twiss centre option')
     interf = 0
     pos_flag = 2
  endif

  k = double_from_table_row(tabname, horname, 1, d_val)
  if (k .lt. 0)  then
     select case (k)
       case (-1)
          write (warnstr,'(3a)') 'table ',tabname,' not found'
       case (-2)
          write (warnstr,'(4a)') 'hor. variable ',horname,' not in table ',tabname
       case default
          write (warnstr,'(3a)') 'table ',tabname,' is empty'
       end select
       call fort_warn('PLOT: ',warnstr)
       ierr = 1
       return
  endif

  !rdemaria 1/6/2005: if the horizontal name is deltap don't plot momentum offset.
  if (horname.eq.'dpp' .or. horname.eq.'deltap') dpp_flag = .true.

  rselect = machp .and. hrange(2).gt.hrange(1)

  do l = 1, nivvar
     k = double_from_table_row(tabname, sname(l), 1, d_val)
     if (k .lt. 0)  then
        write (warnstr,'(4a)') 'vertical variable: ',sname(l)(:lastnb(sname(l))),' not in table ',tabname
        call fort_warn('PLOT: ',warnstr)
        ierr = 1
        return
     endif
     if (sname(l) .eq. 'dpp') dpp_flag = .true.
  enddo

  if (rselect .or. range_plot)  then
     !-------adjust element range to horizontal range
     new1 = nrrang(1)
     new2 = nrrang(2)
     crow = nrrang(1)
     do j = nrrang(1), nrrang(2)
        k = double_from_table_row(tabname, horname, j, d_val)
        tval = d_val
        if (tval .lt. hrange(1)) new1 = j
        if (tval .lt. hrange(2)) new2 = j
     enddo
     nrrang(1) = new1
     if (nrrang(2) .gt. new2+2) nrrang(2) = new2 + 2
  endif

  if (itbv.eq.0 .and. .not.range_plot)  then
     nrrang(1) = 1
     nrrang(2) = table_length(tabname)
  endif

  if (nrrang(1) .eq. 0) nrrang(1) = 1

  !--- get interpolation interval size
  if (machp)  then
     k = double_from_table_row(tabname, horname, nrrang(1), d_val)
     k = double_from_table_row(tabname, horname, nrrang(2), d_val1)
     step = (d_val1 - d_val) / (maxpnt / 2)
  endif
  fact = pos_flag - 1
  j = restart_sequ()

  do j = nrrang(1), nrrang(2)
     crow = j
     if (itbv.eq.1 .and. advance_to_pos(tabname, j).eq.0) goto 100

     k = double_from_table_row(tabname, horname, j, currpos)

     if (itbv .eq. 1)  then
        currtyp = node_value('mad8_type ')
        if (currtyp .eq. 39) currtyp = 15
        if (currtyp .eq. 38) currtyp = 24
        if (currtyp .le. mtype) pltyp = ilist(currtyp)

        !--- get element parameters & build up pltyp (to be used by the routine peschm)
        currleng = node_value('l ')
        ! for codes 1 to 8 only and elements with length
        if (currleng.gt.0.d0 .and. currtyp.gt.1 .and. currtyp.lt.8) then
           currtilt = node_value('tilt ')
           k = double_from_table_row(tabname, 'k1l ' , j, currk1l);  currk1l  = currk1l /currleng
           k = double_from_table_row(tabname, 'k1sl ', j, currk1sl); currk1sl = currk1sl/currleng
           k = double_from_table_row(tabname, 'k2l ' , j, currk2l);  currk2l  = currk2l /currleng
           k = double_from_table_row(tabname, 'k2sl ', j, currk2sl); currk2sl = currk2sl/currleng
           k = double_from_table_row(tabname, 'k3l ' , j, currk3l);  currk3l  = currk3l /currleng
           k = double_from_table_row(tabname, 'k3sl ', j, currk3sl); currk3sl = currk3sl/currleng
        endif

        !--- sbend, rbend
        if (mod(pltyp,20) .eq. 1 .and. currtilt .ne. zero)        pltyp = pltyp + 6
        !--- quad (focusing pltyp=2 or defocusing pltyp=3 ?)
        if (pltyp .eq. 2 .and. min(currk1l, currk1sl) .lt. zero)  pltyp = 3
        !--- sext (focusing pltyp=10 or defocusing pltyp=11 ?)
        if (pltyp .eq. 10 .and. min(currk2l, currk2sl) .lt. zero) pltyp = 11
        !--- oct (focusing pltyp=12 or defocusing pltyp=13 ?)
        if (pltyp .eq. 12 .and. min(currk3l, currk3sl) .lt. zero) pltyp = 13

        !--- Compute the start & end position
        if (machp)  then
           nelmach = nelmach + 1
           estart(nelmach) = currpos - fact * half * currleng
           eend(nelmach)   = estart(nelmach) + currleng
           ieltyp(nelmach) = pltyp
        endif

        !--- Interpolation if required
        if (machp .and. j.gt.nrrang(1) .and. currleng.gt.zero     &
             .and. interf.gt.0 .and. step.gt.zero .and. .not.ptc_flag) then

           nint = currleng / step ! implicit type promotion
           if (nint .lt. 2) nint = 2
           call peintp(crow, nint, proc_n, currleng, ierr)

           if (ierr .eq. 1)  then
              ierr = 0
              call fort_warn('PLOT: ','plot buffer full, plot truncated')
              goto 100
           elseif (ierr .ne. 0)  then
              return
           endif
        endif
     endif

     mystep = step / ten
     do l = 1, nivvar
        if (nqval(l) .eq. maxseql)  then
           call fort_warn('PLOT: ','plot buffer full, plot truncated')
           goto 100
        elseif (nqval(l) .eq. 0)  then
           nqval(l) = nqval(l) + 1
           qhval(nqval(l),l) = currpos
           k = double_from_table_row(tabname, sname(l), j, d_val)
           k = p(l)
           qvval(nqval(l),l) = d_val
           if (proc_flag(1,l) .eq. 1) then ! case of rbetx or rbety
              qvval(nqval(l),l) = sqrt(abs(qvval(nqval(l),l)))
           endif
        elseif (itbv .eq. 0 .or. &
                currpos - qhval(nqval(l),l) .gt. mystep   .or. &
                (marker_plot .and. currtyp .eq. 25)) then
           nqval(l) = nqval(l) + 1
           qhval(nqval(l),l) = currpos
           k = double_from_table_row(tabname, sname(l), j, d_val)
           k = p(l)
           qvval(nqval(l),l) = d_val
           if (proc_flag(1,l) .eq. 1) then ! case of rbetx or rebty
              qvval(nqval(l),l) = sqrt(abs(qvval(nqval(l),l)))
           endif
        endif
     enddo

     k = advance_node()

     if (itbv.eq.1 .and. k.eq.0) goto 100
  enddo

100 continue
  fpmach = machp .and. nelmach.gt.0 .and. noline.eq.0

end subroutine pefill

subroutine pegacn(ncc, window, act, xreal, yreal, np, usex, xwpos, xpos, ypos, ilb)
  use plotfi
  implicit none
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
  !
  ! calls the utility routine iucomp in this file.
  !
  ! it is called by the routine pecurv in this file.
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  integer :: ncc, np, ilb
  real :: window(*), act(*), xreal(*), yreal(*)
  real :: usex, xwpos, xpos, ypos

  !--- type definitions of local variables
  integer :: i, j, iy, iapos, iposx, iposy, kapos(mposx, mposy)
  real :: ywpos, xmax, xmin, xdiff, ydiff, d, t, eps, xadd, yadd
  real :: xdiag(2,2), ydiag(2,2)

  integer, external :: iucomp

  !--- Initialisation of local variables
  save kapos

  !--- Output initialisation
  ilb = 0

  !--- Routine body

  !--- reset position array if first curve in frame
  if (ncc .eq. 1) KAPOS = 0

  xdiff = window(2) - window(1)
  ydiff = window(4) - window(3)
  eps = 1e-6 * max(xdiff, ydiff)
  xmax = maxval(XREAL(1:np))
  xmin = minval(XREAL(1:np))

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
     xpos = act(1) + 0.125 * usex * (iposx - 0.5) * (act(2) - act(1))
     ypos = act(4) - usex * (0.05 * (act(4) - act(3)) + 0.03 * (iposy - 1))

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
               (ydiag(2,j) - ydiag(1,j)) * (xreal(i-1) - xreal(i))
           if (abs(d) .lt. eps) goto 30
           t = (xreal(i-1) - xdiag(1,j)) * (yreal(i-1) - yreal(i)) -   &
               (yreal(i-1) - ydiag(1,j)) * (xreal(i-1) - xreal(i))
           t = t / d
           if (t .lt. 0. .or. t .gt. 1.) goto 30
           t = (xdiag(2,j) - xdiag(1,j)) * (yreal(i-1) - ydiag(1,j)) - &
               (ydiag(2,j) - ydiag(1,j)) * (xreal(i-1) - xdiag(1,j))
           t = t / d
           if (t .ge. 0. .and. t .le. 1.) goto 20
30         continue
        enddo
     enddo
  endif

  if (ilb .gt. 0) KAPOS(iposx, 1:mposy) = 1

  do i = 1, mposx
     do j = 1, mposy
        kapos(i,j) = max(0, kapos(i,j))
     enddo
  enddo
end subroutine pegacn

subroutine pegaxn (nax, vax, sax, ns)
  use plotfi
  implicit none
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
  ! calls the routines gxpnbl in file gxx11.F.
  !
  ! it is called by the routine pemima in this file.
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  character(len=mcnam) :: vax(*), sax(*)
  integer :: nax, ns

  !--- type definitions of local variables
  character(len=mcnam) :: scut, saloc
  integer :: i, k, k1, k2, j, k1f, k2f

  !--- Initialisation of local variables
  ns = 0
  sax(1) = ' '

  !--- Routine body

  if (nax .le. 0) return

  do  i = 1, nax
     saloc = vax(i)
     call gxpnbl(saloc, k1, k2)

     ! try to compound similar axes annotation into one, eg betx and bety into bet
     if (k2 .gt. 1 .and. index('XY', saloc(k2:k2)) .ne. 0)  then
        scut = saloc(:k2-1)
        do j = 1, ns
           if (scut .eq. sax(j))  goto 10
        enddo
        do j = i + 1, nax
           call gxpnbl(vax(j), k1f, k2f)
           if (k2.eq.k2f .and. index('XY',vax(j)(k2:k2)).ne.0 .and. saloc(:k2-1).eq.vax(j)(:k2-1)) then
              saloc = scut
              do k = 1, ns
                 if (saloc .eq. sax(k))  goto 10
              enddo
           endif
        enddo
     endif

     ns      = ns + 1
     sax(ns) = saloc

10   continue
  enddo

end subroutine pegaxn

subroutine pegetn (iflag, svar, it, ipflg, sovar, reqann)
  use plotfi
  use plot_bfi
  implicit none
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
  !   reqann      (char)  requested annotation                           *
  !                                                                      *
  ! calls the utility routine pupnbl in this file.                       *
  !                                                                      *
  ! it is called by the routines pemima, peplot and pesopt in this file. *
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  integer :: iflag, it, ipflg(2)
  character(len=mcnam) :: svar, sovar(*)
  character(len=*) :: reqann

  !--- type definitions of local variables
  character(len=mxlabl) :: svlabl(mnvar), svanno(mnvar)
  character(len=mcnam)  :: svname(mnvar)
  !--- strings:
  !   SVLABL   plot prescriptions for variables on axis labels
  !   SVANNO   plot prescriptions for variables in annotations
  !   SVNAME   names of variables known to the program
  integer :: i, iref, k1, k2, k1f, k2f, j
  integer :: iproc(mnvar,3), intpo(mnvar), ivdep(mnvar,mxdep,3)

  !--- Initialisation of local variables
  data (svname(j), j = 1, 32) /                                     &
       's', 'size', 'deltap',                                            &
       'qs', 'x', 'y', 'xsize', 'ysize',                                 &
       'dt', 'xn', 'yn', 'pxn', 'pyn',                                   &
       'gammatr', 'xrms', 'yrms',                                        &
       'xmax', 'ymax', 'bxmax',                                          &
       'bymax', 'dxmax', 'dymax',                                        &
       'tn', 't', 'turns', 'particle', 'alfa',                           &
       'ptn', 'wt', 'phit',                                              &
       'rbxmax',                                                         &
       'rbymax' /

  data (svname(j), j = 33, mnvar) /                                 &
       'betx', 'rbetx',                                                  &
       'alfx', 'mux', 'dx',                                              &
       'dpx', 'qx', 'px', 'wx',                                          &
       'phix', 'dmux',                                                   &
       'ddx', 'ddpx', 'iwx',                                             &
       'xix',                                                            &
       'bety', 'rbety',                                                  &
       'alfy', 'muy', 'dy',                                              &
       'dpy', 'qy', 'py', 'wy',                                          &
       'phiy', 'dmuy',                                                   &
       'ddy', 'ddpy', 'iwy',                                             &
       'xiy', 'xns', 'pxns', 'wxs',                                      &
       'yns', 'pyns', 'wys',                                             &
       'energy', 'spintune',                                             &
       'poltotal', 'poldiffx', 'poldiffy', 'poldiffs' /

  data (svlabl(j), j = 1, 32) /                                     &
       's (m)', 'n<G>s<G> (mm)', '<G>d<G><?>E<?>/p<?>0<?>c',             &
       'Q<?>s<?>', 'x (m)', 'y (m)', 'n<G>s<G> (mm)', 'n<G>s<G> (mm)',   &
       'ct (m)', 'x<?>n<?>', 'y<?>n<?>', 'p<?>xn<?>', 'p<?>yn<?>',       &
       '<G>g<G><?>tr<?>', 'X<?>rms<?> (m)', 'Y<?>rms<?> (m)',            &
       'X<?>max<?> (m)', 'Y<?>max<?> (m)', '<G>b<G><?>x_max<?> (m)',     &
       '<G>b<G><?>y_max<?> (m)', 'D<?>x_max<?> (m)', 'D<?>y_max<?> (m)', &
       't<?>n<?>', 'ct (m)', 'turns', 'particle', '<G>a<G>',             &
       'p<?>t_n<?>', 'W<?>t<?>', '<G>F<G><?>t<?> (rad/2<G>p<G>)',        &
       '<G>b<G><?>x_max<?><!>1/2<!> (m<!>1/2<!>)',                       &
       '<G>b<G><?>y_max<?><!>1/2<!> (m<!>1/2<!>)' /

  data (svlabl(j), j = 33, mnvar) /                                 &
       '<G>b<G><?>x<?> (m)', '<G>b<G><?>x<?><!>1/2<!> (m<!>1/2<!>)',     &
       '<G>a<G><?>x<?>', '<G>m<G><?>x<?> (rad/2<G>p<G>)', 'D<?>x<?> (m)',&
       'D<?>px<?>', 'Q<?>x<?>', 'p<?>x<?>/p<?>0<?>', 'W<?>x<?>',         &
       '<G>F<G><?>x<?> (rad/2<G>p<G>)', 'd<G>m<G><?>x<?>/d<G>d<G>',      &
       'dD<?>x<?>/d<G>d<D> (m)', 'dD<?>px<?>/d<G>d<G>', 'W<?>x<?> (m)',  &
       'XI<?>x<?>',                                                      &
       '<G>b<G><?>y<?> (m)', '<G>b<G><?>y<?><!>1/2<!> (m<!>1/2<!>)',     &
       '<G>a<G><?>y<?>', '<G>m<G><?>y<?> (rad/2<G>p<G>)', 'D<?>y<?> (m)',&
       'D<?>py<?>', 'Q<?>y<?>', 'p<?>y<?>/p<?>0<?>', 'W<?>y<?>',         &
       '<G>F<G><?>y<?> (rad/2<G>p<G>)', 'd<G>m<G><?>y<?>/d<G>d<G>',      &
       'dD<?>y<?>/d<G>d<D> (m)', 'dD<?>py<?>/d<G>d<G>', 'W<?>y<?> (m)',  &
       'XI<?>y<?>', 'x<?>ns<?>', 'p<?>x_ns<?>', 'W<?>xs<?>',             &
       'y<?>ns<?>', 'p<?>y_ns<?>', 'W<?>ys<?>',                          &
       'E[GeV]', 'spintune',                                             &
       'polarization','polarization','polarization','polarization' /

  data (svanno(j), j = 1, 32) /                                     &
       's', 'n<G>s<G>', '<G>d<G>',                                       &
       'Q<?>s<?>', 'x', 'y', 'n<G>s<G><?>x<?>', 'n<G>s<G><?>y<?>',       &
       'ct', 'x<?>n<?>', 'y<?>n<?>', 'p<?>xn<?>', 'p<?>yn<?>',           &
       '<G>g<G><?>tr<?>', 'X<?>rms<?>', 'Y<?>rms<?>',                    &
       'X<?>max<?>', 'Y<?>max<?>', '<G>b<G><?>x_max<?>',                 &
       '<G>b<G><?>y_max<?>', 'D<?>x_max<?>', 'D<?>y_max<?>',             &
       't<?>n<?>', 't', 'turns', 'particle', '<G>a<G>',                  &
       'p<?>t_n<?>', 'W<?>t<?>', '<G>F<G><?>t<?>',                       &
       '<G>b<G><?>x_max<?><!>1/2<!>',                                    &
       '<G>b<G><?>y_max<?><!>1/2<!>' /

  data (svanno(j), j = 33, mnvar) /                                 &
       '<G>b<G><?>x<?>', '<G>b<G><?>x<?><!>1/2<!>',                      &
       '<G>a<G><?>x<?>', '<G>m<G><?>x<?>', 'D<?>x<?>',                   &
       'D<?>px<?>', 'Q<?>x<?>', 'p<?>x<?>', 'W<?>x<?>',                  &
       '<G>F<G><?>x<?>', '<G>m<G><?>x<?>''',                             &
       'D<?>x<?>''', 'D<?>px<?>''', 'W<?>x<?>',                          &
       'XI<?>x<?>',                                                      &
       '<G>b<G><?>y<?>', '<G>b<G><?>y<?><!>1/2<!>',                      &
       '<G>a<G><?>y<?>', '<G>m<G><?>y<?>', 'D<?>y<?>',                   &
       'D<?>py<?>', 'Q<?>y<?>', 'p<?>y<?>', 'W<?>y<?>',                  &
       '<G>F<G><?>y<?>', '<G>m<G><?>y<?>''',                             &
       'D<?>y<?>''', 'D<?>py<?>''', 'W<?>y<?>',                          &
       'XI<?>y<?>', 'x<?>ns<?>', 'p<?>x_ns<?>', 'W<?>xs<?>',             &
       'y<?>ns<?>', 'p<?>y_ns<?>', 'W<?>ys<?>',                          &
       'E', ' ',                                                         &
       'p<?>tot<?>', 'p<?>diff_x<?>', 'p<?>diff_y<?>', 'p<?>diff_s<?>'/

  data (iproc(j,1), j = 1, 32) /                                    &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 2, 3, 4, 5,                                                    &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       6, 0, 0, 0, 0,                                                    &
       7, 8, 9,                                                          &
       1,                                                                &
       1 /

  data (iproc(j,1), j = 33, mnvar) /                                &
       0, 1,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 0,                                                       &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0,                                                                &
       0, 1,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 0,                                                       &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 14, 15, 16,                                                    &
       17, 18, 19,                                                       &
       0, 0,                                                             &
       0, 0, 0, 0 /

  data (iproc(j,2), j = 1, 32) /                                    &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 2, 3, 4, 5,                                                    &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       6, 0, 0, 0, 0,                                                    &
       7, 8, 9,                                                          &
       1,                                                                &
       1 /

  data (iproc(j,2), j = 33, mnvar) /                                &
       0, 1,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 10,                                                      &
       11, 0,                                                            &
       0, 0, 0,                                                          &
       0,                                                                &
       0, 1,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 12,                                                      &
       13, 0,                                                            &
       0, 0, 0,                                                          &
       0, 14, 15, 16,                                                    &
       17, 18, 19,                                                       &
       0, 0,                                                             &
       0, 0, 0, 0 /

  data (iproc(j,3), j = 1, 32) /                                    &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 2, 3, 4, 5,                                                    &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       6, 0, 0, 0, 0,                                                    &
       7, 8, 9,                                                          &
       1,                                                                &
       1 /

  data (iproc(j,3), j = 33, mnvar) /                                &
       0, 1,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 10,                                                      &
       11, 0,                                                            &
       0, 0, 0,                                                          &
       0,                                                                &
       0, 1,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 12,                                                      &
       13, 0,                                                            &
       0, 0, 0,                                                          &
       0, 14, 15, 16,                                                    &
       17, 18, 19,                                                       &
       0, 0,                                                             &
       0, 0, 0, 0 /

  !--- in INTPO, n+100 means: take SQRT of var. n
  ! 2014-May-05  16:04:47  ghislain: added proper interpolation for x (code 11),
  ! px (12), y (13), py(14)
  data (intpo(j), j = 1, 32) / &
       0, 0, 0, &
       0, 11, 13, 0, 0, &
       0, 0, 0, 0, 0, &
       0, 0, 0, &
       0, 0, 0, &
       0, 0, 0, &
       0, 0, 0, 0, 0, &
       0, 0, 0, &
       0, &
       0/

  data (intpo(j), j = 33, mnvar) /                                  &
       1, 101,                                                           &
       2, 3, 4,                                                          &
       5, 0, 12, 0,                                                       &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0,                                                                &
       6, 106,                                                           &
       7, 8, 9,                                                          &
       10, 0, 14, 0,                                                      &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 0,                                                       &
       0, 0, 0,                                                          &
       0, 0,                                                             &
       0, 0, 0, 0 /

  data (ivdep(j,1,1), j = 1, 32) /                                  &
       1, 2, 3,                                                          &
       4, 5, 6, 7, 8,                                                    &
       9, 5, 6, 5, 6,                                                    &
       14, 15, 16,                                                       &
       17, 18, 19,                                                       &
       20, 21, 22,                                                       &
       24, 24, 25, 26, 27,                                               &
       3, 3, 3,                                                          &
       19,                                                               &
       20 /

  data (ivdep(j,2,1), j = 1, 32) /                                  &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 0, 0, 40, 55,                                                  &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 24, 24,                                                        &
       0,                                                                &
       0 /

  data (ivdep(j,1,1), j = 33, mnvar) /                              &
       33, 33,                                                           &
       35, 36, 37,                                                       &
       38, 39, 40, 41,                                                   &
       42, 43,                                                           &
       44, 45, 46,                                                       &
       47,                                                               &
       48, 48,                                                           &
       50, 51, 52,                                                       &
       53, 54, 55, 56,                                                   &
       57, 58,                                                           &
       59, 60, 61,                                                       &
       62, 5, 5, 5,                                                      &
       6, 6, 6,                                                          &
       69, 70,                                                           &
       71, 72, 73, 74 /

  data (ivdep(j,2,1), j = 33, mnvar) /                              &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 0,                                                       &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0,                                                                &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 0,                                                       &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 40, 40,                                                     &
       0, 55, 55,                                                        &
       0, 0,                                                             &
       0, 0, 0, 0 /

  data (ivdep(j,1,2), j = 1, 32) /                                  &
       1, 2, 3,                                                          &
       4, 5, 6, 7, 8,                                                    &
       9, 5, 6, 5, 6,                                                    &
       14, 15, 16,                                                       &
       17, 18, 19,                                                       &
       20, 21, 22,                                                       &
       24, 24, 25, 26, 27,                                               &
       3, 3, 3,                                                          &
       19,                                                               &
       20 /

  data (ivdep(j,2,2), j = 1, 32) /                                  &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 0, 0, 40, 55,                                                  &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 24, 24,                                                        &
       0,                                                                &
       0 /

  data (ivdep(j,1,2), j = 33, mnvar) /                              &
       33, 33,                                                           &
       35, 36, 37,                                                       &
       38, 39, 40, 5,                                                    &
       5, 43,                                                            &
       44, 45, 46,                                                       &
       47,                                                               &
       48, 48,                                                           &
       50, 51, 52,                                                       &
       53, 54, 55, 6,                                                    &
       6, 58,                                                            &
       59, 60, 61,                                                       &
       62, 5, 5, 5,                                                      &
       6, 6, 6,                                                          &
       69, 70,                                                           &
       71, 72, 73, 74 /

  data (ivdep(j,2,2), j = 33, mnvar) /                              &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 40,                                                      &
       40, 0,                                                            &
       0, 0, 0,                                                          &
       0,                                                                &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 55,                                                      &
       55, 0,                                                            &
       0, 0, 0,                                                          &
       0, 0, 40, 40,                                                     &
       0, 55, 55,                                                        &
       0, 0,                                                             &
       0, 0, 0, 0 /

  data (ivdep(j,1,3), j = 1, 32) /                                  &
       1, 2, 3,                                                          &
       4, 5, 6, 7, 8,                                                    &
       9, 5, 6, 5, 6,                                                    &
       14, 15, 16,                                                       &
       17, 18, 19,                                                       &
       20, 21, 22,                                                       &
       24, 24, 25, 26, 27,                                               &
       3, 3, 3,                                                          &
       19,                                                               &
       20 /

  data (ivdep(j,2,3), j = 1, 32) /                                  &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 0, 0, 40, 55,                                                  &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0,                                                          &
       0, 0, 0, 0, 0,                                                    &
       0, 24, 24,                                                        &
       0,                                                                &
       0 /

  data (ivdep(j,1,3), j = 33, mnvar) /                              &
       33, 33,                                                           &
       35, 36, 37,                                                       &
       38, 39, 40, 5,                                                    &
       5, 43,                                                            &
       44, 45, 46,                                                       &
       47,                                                               &
       48, 48,                                                           &
       50, 51, 52,                                                       &
       53, 54, 55, 6,                                                    &
       6, 58,                                                            &
       59, 60, 61,                                                       &
       62, 5, 5, 5,                                                      &
       6, 6, 6,                                                          &
       69, 70,                                                           &
       71, 72, 73, 74 /

 data (ivdep(j,2,3), j = 33, mnvar) /                              &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 40,                                                      &
       40, 0,                                                            &
       0, 0, 0,                                                          &
       0,                                                                &
       0, 0,                                                             &
       0, 0, 0,                                                          &
       0, 0, 0, 55,                                                      &
       55, 0,                                                            &
       0, 0, 0,                                                          &
       0, 0, 40, 40,                                                     &
       0, 55, 55,                                                        &
       0, 0,                                                             &
       0, 0, 0, 0 /

  !--- Routine body

  if (it .le. 0 .or. it .gt. 3)  then
     sovar(1) = svar
     sovar(2) = ' '
     reqann   = svar
     ipflg(1) = 0
     ipflg(2) = 0
     return
  endif

  sovar(1) = ' '
  reqann = svar

  !--- search in list of known variables
  do  iref = 1, mnvar
     if (svar .eq. svname(iref))  goto 9
  enddo

  call pupnbl(svar, k1, k2)
  do  iref = 1, mnvar
     call pupnbl(svname(iref), k1f, k2f)
     if ( k2+1 .eq. k2f .and. svar(:k2) .eq. svname(iref)(:k2) .and. &
          index('xy', svname(iref)(k2f:k2f)) .ne. 0 ) goto 9
  enddo

  return

9 continue ! found variable name as svname(iref)

  select case (iflag)
    case (0)
       reqann = svname(iref)
       ipflg(1) = iproc(iref,it)
       ipflg(2) = intpo(iref)
       do  j = 1, mxdep
        sovar(j) = ' '
        if (ivdep(iref,j,it) .ne. 0) &
             sovar(j) = svname(ivdep(iref,j,it))
       enddo

    case (1)
       reqann = svlabl(iref)
       if (svar .ne. svname(iref))  then
          !--- incomplete match
          !    replace x or y in name by blank
          call pupnbl(reqann, k1, k2)
          do  i = 2, k2
             if (index('XYxy', reqann(i:i)) .ne. 0) reqann(i:i) = ' '
          enddo
       endif

    case(2)
       reqann = svanno(iref)

    case (3)
       if (svar .eq. svname(iref))  then
          reqann = svname(iref)
       else
          reqann = svname(iref)(:k2)
       endif

    case default
       reqann = svar

  end select

end subroutine pegetn

subroutine peiact(kact, np, x, y, ac, kf, kl)
  implicit none
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

  !--- type definition of the routine arguments
  integer :: kact, np, kf, kl
  real :: x(*), y(*), ac(4)

  !--- type definitions of local variables
  integer :: i
  real :: xtol, ytol
  real, parameter :: toleps=1.e-5

  !--- Routine body

  xtol = toleps * (ac(2) - ac(1))
  ytol = toleps * (ac(4) - ac(3))
  kf = 0
  kl = 0
  do i = kact, np
     if (x(i) + xtol .lt. ac(1)) goto 10
     if (x(i) - xtol .gt. ac(2)) goto 10
     if (y(i) + ytol .lt. ac(3)) goto 10
     if (y(i) - ytol .gt. ac(4)) goto 10
     kf = i
     goto 20
10   continue
  enddo

  !--- no point inside
  return

20 continue
  do i = kf, np
     if (x(i) + xtol .lt. ac(1)) goto 40
     if (x(i) - xtol .gt. ac(2)) goto 40
     if (y(i) + ytol .lt. ac(3)) goto 40
     if (y(i) - ytol .gt. ac(4)) goto 40
  enddo

40 kl = i - 1

end subroutine peiact

subroutine peintp(crow, nint, proc, length, ierr)
  use plotfi
  use plot_bfi
  use math_constfi, only : zero, one
  implicit none
  !----------------------------------------------------------------------*
  !     purpose:
  !     interpolate variables plotted against s
  !     input:
  !     crow        (int)   table row number at start of element
  !     nint        (int)   number of interpolation intervals
  !     type        (int)   (local) element type
  !     proc        (int)   original process flags
  !     step        (d.p.)  max. dist. between two successive hor. values
  !     length      (d.p.)  element length
  !     output:
  !     ierr        (int)   0 if ok, else > 0
  !     the results are stored in qhval and qvval
  !
  !     calls the functions double_from_table_row and get_value
  !     defined in file madxn.c.
  !     calls the function peelma in this file.
  !
  !     it is called by the routine pefill in this file.
  !----------------------------------------------------------------------*

  !---  type definition of the routine arguments
  integer :: crow, nint, proc(2,*), ierr
  double precision :: length

  !---  type definitions of local variables
  integer :: i, j, k, ipc
  double precision :: tw1(mintpl)
  double precision :: ex, ey
  double precision :: xn1, pxn1, yn1, pyn1
  double precision :: s_elem, s_incr, s, gamx, gamy

  !---  definitions of function primitives
  integer, external :: double_from_table_row, interpolate_node, reset_interpolation
  integer, external :: embedded_twiss
  double precision, external :: get_value

  !---  Output initialisation
  ierr = 0

  !---  Routine body

  if (crow .eq. 1) then
     k = double_from_table_row(tabname, 'betx ', 1, tw1(1))
     k = double_from_table_row(tabname, 'alfx ', 1, tw1(2))
     k = double_from_table_row(tabname, 'mux ',  1, tw1(3))
     k = double_from_table_row(tabname, 'dx ',   1, tw1(4))
     k = double_from_table_row(tabname, 'dpx ',  1, tw1(5))

     k = double_from_table_row(tabname, 'bety ', 1, tw1(6))
     k = double_from_table_row(tabname, 'alfy ', 1, tw1(7))
     k = double_from_table_row(tabname, 'muy ',  1, tw1(8))
     k = double_from_table_row(tabname, 'dy ',   1, tw1(9))
     k = double_from_table_row(tabname, 'dpy ',  1, tw1(10))

     k = double_from_table_row(tabname, 'x ',    1, tw1(11))
     k = double_from_table_row(tabname, 'px ',   1, tw1(12))

     k = double_from_table_row(tabname, 'y ',    1, tw1(13))
     k = double_from_table_row(tabname, 'py ',   1, tw1(14))

     k = double_from_table_row(tabname, 's ', 1, s_incr)
     s = zero
     ex = get_value('beam ','ex ')
     ey = get_value('beam ','ey ')

     !---  xn, pxn, yn, pyn
     xn1 = zero
     if (ex*tw1(1) .ne. zero) xn1 = tw1(11) / sqrt(ex*tw1(1))
     yn1 = zero
     if (ey*tw1(6) .ne. zero) yn1 = tw1(13) / sqrt(ey*tw1(6))

     ! 2015-Aug-09  10:56:43  ghislain: added initialisation of gamx gamy
     gamx = zero
     if (tw1(1).ne.zero) gamx = (one + tw1(2)**2) / tw1(1)
     gamy = zero
     if (tw1(6).ne.zero) gamy = (one + tw1(7)**2) / tw1(6)

     pxn1 = tw1(12) * gamx
     pyn1 = tw1(14) * gamy
     tw1(15) = xn1
     tw1(16) = pxn1
     tw1(17) = yn1
     tw1(18) = pyn1

     !---  loop over variables, interpolate those with codes
     do j = 1, nivvar
        ipc = mod(proc(2,j), 100)
        if (ipc .gt. 0)  then
           if (nqval(j) .eq. maxseql)  then
              ierr = 1
              return
           endif

           ipparm(2,j) = 1
           nqval(j) = nqval(j) + 1
           qhval(nqval(j),j) = s
           if (proc(1,j) .gt. 0)  then
              qvval(nqval(j), j) = sqrt(abs(tw1(ipc)))
           else
              qvval(nqval(j), j) = tw1(ipc)
           endif
        else
           ipparm(2,j) = 0
        endif
     enddo

    return
  endif ! crow .eq. 1

  if (length .eq. zero)  return

  ! now crow was not 1 and length is non-zero

  k = double_from_table_row(tabname, horname, crow - 1, s_elem)
  !---  set flag for correct interpolation
  !---  get intermediate s values, and interpolate twiss parameters

  k = interpolate_node(nint)
  k = embedded_twiss()
  do i = 1, nint
     k = double_from_table_row('embedded_twiss_table ', 'betx ', i, tw1(1))
     k = double_from_table_row('embedded_twiss_table ', 'alfx ', i, tw1(2))
     k = double_from_table_row('embedded_twiss_table ', 'mux ',  i, tw1(3))
     k = double_from_table_row('embedded_twiss_table ', 'dx ',   i, tw1(4))
     k = double_from_table_row('embedded_twiss_table ', 'dpx ',  i, tw1(5))

     k = double_from_table_row('embedded_twiss_table ', 'bety ', i, tw1(6))
     k = double_from_table_row('embedded_twiss_table ', 'alfy ', i, tw1(7))
     k = double_from_table_row('embedded_twiss_table ', 'muy ',  i, tw1(8))
     k = double_from_table_row('embedded_twiss_table ', 'dy ',   i, tw1(9))
     k = double_from_table_row('embedded_twiss_table ', 'dpy ',  i, tw1(10))

     k = double_from_table_row('embedded_twiss_table ', 'x ',    i, tw1(11))
     k = double_from_table_row('embedded_twiss_table ', 'px ',   i, tw1(12))

     k = double_from_table_row('embedded_twiss_table ', 'y ',    i, tw1(13))
     k = double_from_table_row('embedded_twiss_table ', 'py ',   i, tw1(14))

     k = double_from_table_row('embedded_twiss_table ', 's ', i, s_incr)
     s = s_elem + s_incr
     ex = get_value('beam ','ex ')
     ey = get_value('beam ','ey ')

     !---  xn, pxn, yn, pyn
     xn1 = zero
     if (ex*tw1(1) .ne. zero) xn1 = tw1(11) / sqrt(ex*tw1(1))
     yn1 = zero
     if (ey*tw1(6) .ne. zero) yn1 = tw1(13) / sqrt(ey*tw1(6))

     gamx = zero
     if (tw1(1) .ne. zero) gamx = (one + tw1(2)**2) / tw1(1)
     gamy = zero
     if (tw1(6) .ne. zero) gamy = (one + tw1(7)**2) / tw1(6)

     pxn1 = tw1(12) * gamx
     pyn1 = tw1(14) * gamy
     tw1(15) = xn1
     tw1(16) = pxn1
     tw1(17) = yn1
     tw1(18) = pyn1

     !---  loop over variables, interpolate those with codes
     do j = 1, nivvar
        ipc = mod(proc(2,j), 100)
        if (ipc .gt. 0)  then
           if (nqval(j) .eq. maxseql)  then
              ierr = 1
              return
           endif

           ipparm(2,j) = 1
           nqval(j) = nqval(j) + 1
           qhval(nqval(j),j) = s
           if (proc(1,j) .gt. 0)  then
              qvval(nqval(j), j) = sqrt(abs(tw1(ipc)))
           else
              qvval(nqval(j), j) = tw1(ipc)
           endif
        else
           ipparm(2,j) = 0
        endif
     enddo

  enddo ! loop over nint
  k = reset_interpolation()

end subroutine peintp

subroutine pemima
  use plotfi
  use plot_bfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Constrain axis reference, find minima and maxima of coordinates,   *
  !   construct axis labels                                              *
  ! calls the routine gxpnbl in file gxx11.F.                            *
  ! calls the routines pegetn and pegaxn in this file.                   *
  !                                                                      *
  ! it is called by the routine exec_plot in file madxn.c after pefill.  *
  !----------------------------------------------------------------------*

  !--- type definitions of local variables
  integer :: i, j, k, iv, idum(2), ns, k1, k2, i1, i2, it(4)
  character(len=mtitl) :: s
  character(len=mxlabl) :: slab
  character(len=mcnam) :: sdum(mxcurv), saxis(mxcurv), vaxis(mxcurv,4)

  !--- Initialisation of variables in common peaddi
  numax = 0

  !--- Initialisation of local variables
  IT = 0

  !--- Routine body
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
        if (it(j) .gt. it(j+1)) call swapint(it(j), it(j+1))
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

  do j = 1, nivvar ! loop over all variables to be plotted
     k = naxref(j)
     do i = 1, nqval(j)
        ! find overall hmin and hmax
        hmima(1) = min(hmima(1), qhval(i,j))
        hmima(2) = max(hmima(2), qhval(i,j))
        ! find vmin and vmax for each axis
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

        if (nvvar(iv) .eq. 1)  then ! only one variable to be plotted
           call pegetn (1, vaxis(1,iv), itbv, idum, sdum, slab)
           ns = 1
        else ! several variables to be plotted
           call pegaxn (nvvar(iv), vaxis(1,iv), saxis, ns)
           call pegetn (1, saxis(1), itbv, idum, sdum, slab)
        endif

        call gxpnbl (slab, k1, k2) ! extract indices for first and last non-blank characters in string
        s  = '<#>' // slab
        k2 = k2 + 3

        do i = 2, ns ! extract the labels for other variables
           call pegetn (1, saxis(i), itbv, idum, sdum, slab)
           call gxpnbl (slab, i1, i2)
           if (k2+i2+2 .gt. mtitl) then
              call fort_warn('PLOT: ','Array index larger than array limit; label is truncated')
              goto 100 ! stop building label
           endif
           if (index(s(:k2),slab(:i2)) .eq. 0)  then
              s(k2+1:) = ', ' // slab(:i2)
              k2 = k2 + i2 + 2
           endif
        enddo

100     axlabel(iv) = s
     endif
  enddo
end subroutine pemima

subroutine peplot
  use plotfi
  use plot_bfi
  use math_constfi, only : zero, one, ten3m
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Plot all types of graphs from MAD.
  !   Uses GXPLOT with underlying X-Windows (PostScript)
  !
  ! calls the function plot_option and the routine comm_para
  ! defined in file madxn.c.
  ! calls the routines pegetn, pecurv and peschm in this file.
  ! calls the routines gxsdef, gxsvar, jslwsc, gxsvpt, gxpnbl, gxqaxs,
  !                    gxsaxs, iaxseq, gxqcrv, gxscrv, axlabel, gxfrm1,
  !                    gxswnd  defined in file gxx11.F.
  !
  ! it is called by the routine plotit in this file.
  !----------------------------------------------------------------------*
  !--- type definitions of local variables
  !--- strings:
  !   SVAR     buffer for variable names etc.
  !   SLOCN    local name buffer (without leading "_")
  !   STEMP    temporary buffer for titles
  !   STEXT    buffer for labels etc.
  !   SFORM    format buffer
  !--- reals:
  !   PRMACH fraction of viewport taken by machine plot
  !   SYMCH  preset symbol character height

  integer :: ipar(50), nptval(4), ipxval(4), ipyval(4), icvref(4)
  integer :: lastnb, iaxseq(4)
  integer :: idum, k1dum, k2dum, k3dum, vdum(2)
  integer :: i, j, k, npar, ivvar, nvax, ivax, ierr
  double precision :: deltap
  real :: prmach, symch, tmpval, yvtop, fdum, chh
  real :: vpt(4), window(4,4), actwin(4,4), range(2), xax(2), yax(8)

  character(len=mcnam)  :: svar
  character(len=mxlabl) :: slocn, slname
  character(len=mtitl)  :: stemp
  character(len=300)    :: stext
  character(len=20)     :: sform
  character(len=mcnam)  :: sdum(mxdep)
  character(len=1)      :: ssymb
  character(len=80)     :: ch

  !--- definitions of function primitives
  double precision, external :: plot_option
  integer, external :: double_from_table_row, double_from_table_header
  integer, external :: table_column_exists, table_header_exists

  !--- Initialisation of local variables
  data prmach /0.1/, symch /0.01/
  data iaxseq / 1, 4, 2, 3 /

  ssymb = ' '
  fdum = 0.e0
  chh = 0.e0
  VPT = 0.e0
  NPTVAL = 0
  IPXVAL = 0
  IPYVAL = 0
  ICVREF = 0
  WINDOW = 0.e0
  ACTWIN = 0.e0
  RANGE = 0.e0
  XAX = 0.e0
  YAX = 0.e0
  IPAR = 0

  !--- Routine body

  !--- Acquire deltap
  deltap = zero
  if (.not.ptc_flag) then
     if (tabname .eq. "summ") then
        if (table_column_exists('summ ', 'deltap ').ne.0) &
             k = double_from_table_row('summ ', 'deltap ', 1, deltap)
     else if (table_header_exists(tabname, 'deltap ').ne.0) then
        k = double_from_table_header(tabname, 'deltap ', deltap)
     endif
  else if (table_header_exists(tabname, 'deltap ').ne.0) then
     k = double_from_table_header(tabname, 'deltap ', deltap)
  endif

  write (ch,'(f9.4)') deltap * 100.d0
  if (dpp_flag) then
     slocn = ' <#>'
  else
     slocn = 'Momentum offset = '//ch(:7)//' %'//'&'//'<#>'
  endif

  !--- reset axis and curve defaults
  call gxsdef ('AXIS', 0)
  call gxsdef ('CURVE', 0)

  !--- set "new line" character (change default = '/')
  call gxsvar ('SDEFNL', idum, fdum, '&')

  !--- set top of viewport - leave space to plot machine if required
  yvtop = 1.
  if (fpmach) yvtop = 1. - prmach

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
  slname = "_"
  call pegetn (1, svar, itbv, vdum, sdum, slname)
  call gxpnbl(slname, k1dum, k2dum)
  k3dum = 31
  if (dpp_flag) k3dum = 4
  slocn = ' '
  if (k1dum.gt.0 .and. k2dum.gt.0) then
     do idum = k1dum, k2dum
        if (slname(idum:idum) .ne. '_') then
           k3dum = k3dum + 1
           slocn(k3dum:k3dum) = slname(idum:idum)
        endif
     enddo
  endif

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
        stext = slocn(:lastnb(slocn))
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

     !--- range centre etc
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
     stext = stemp
     if (nvax .eq. 1) stext = '&' // stemp
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
  call gxfrm1 (numax, nptval, ipxval, ipyval, icvref, xax, yax, window, actwin, ierr)
  if (ierr .ne. 0) return

  !--- now loop over vertical variables for real curve plotting
  do ivvar = 1, nivvar
     nvax = naxref(ivvar)
     ivax = iaxseq(nvax)

     !--- find variable name in list for annotation
     svar = slabl(ivvar)
     slname = "_"
     call pegetn (2, svar, itbv, vdum, sdum, slname)
     call gxpnbl(slname, k1dum, k2dum)
     k3dum = 0
     slocn = ' '
     if (k1dum.gt.0 .and. k2dum.gt.0) then
        do idum = k1dum, k2dum
           if (slname(idum:idum) .ne. '_') then
              k3dum = k3dum + 1
              slocn(k3dum:k3dum) = slname(idum:idum)
           endif
        enddo
     endif

     !--- character height including user request
     chh = 0.001 * masize * qascl

     !--- call curve plot routine with simple arrays and flags
     call pecurv (ivvar, slocn, chh, qascl,                          &
          symch * qsscl, ipparm(1,ivvar), nqval(ivvar), qhval(1,ivvar),     &
          qvval(1,ivvar), window(1,nvax), actwin(1,nvax), ierr)

     if (ierr .ne. 0) return
  enddo

  if (fpmach)  then
     vpt(1) = zero
     vpt(2) = one
     vpt(3) = yvtop
     vpt(4) = one
     call gxsvpt (vpt)
     window(3,1) = -one
     window(4,1) = one
     call gxswnd (window)
     call peschm (nelmach, ieltyp, xax, estart, eend, actwin)
  endif

end subroutine peplot

subroutine peschm (nel, ityp, hr, es, ee, actwin)
  use plotfi
  implicit none
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
  !                       37: crab cavity                                *
  !   hr          (real)  horizontal range (lower and upper)
  !   es          (real)  array with element start position
  !   ee          (real)  array with element end position
  !   actwin      (real)  active window for curve plot (array of 4)
  !
  ! calls the routines jsln, gvpl defined in file gx11.F
  ! it is called by the routine peplot in this file.
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  integer :: nel, ityp(*)
  real :: hr(2), es(*), ee(*), actwin(4)

  !--- type definitions of local variables
  integer :: i, j, it, j_nodrift, im1
  integer :: npst(mobj), npnd(mobj), npsl(msize), i_nodrift(maxseql)
  real :: ell, shapex(msize), shapey(msize)
  real :: txp(2), typ(2), typz(2)

  !--- Initialisation of local variables
  data npst   / 1,  6, 11, 16, 21,                                  &
       33, 43, 48,                                                       &
       50,                                                               &
       64, 69, 74, 79, 84 /
  data npnd   / 5, 10, 15, 20, 32,                                  &
       42, 47, 49,                                                       &
       63,                                                               &
       68, 73, 78, 83, 88 /
  data npsl   /5 * 1, 5 * 1, 5 * 1, 5 * 3, 5 * 1, 0, 4 * 1, 0, 1,   &
       1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 5 * 1, 2 * 1,                       &
       6 * 1, 0, 5 * 1, 0, 1,                                            &
       5 * 1, 5 * 1, 5 * 1, 5 * 1, 5 * 1 /
  data typz   / 2 * 0. /

  data shapex /0., 1., 1., 0., 0.,                                  &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0., 0., 1., 1., 0., 0., 0., 1.,                   &
       0., 1., 0.5, 0.5, 0., 1., 0.5, 0.5, 0., 1.,                       &
       0., 1., 1., 0., 0.,                                               &
       0., 0.,                                                           &
       0., 0.25, 0.25, 0.75, 0.75, 1.,                                   &
       0., 0.25, 0.25, 0.75, 0.75, 1., 0., 1.,                           &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0.,                                               &
       0., 1., 1., 0., 0. /

  data shapey /0.6, 0.6, -0.6, -0.6, 0.6,                           &
       0., 0., 0.8, 0.8, 0.,                                             &
       0., 0., -0.8, -0.8, 0.,                                           &
       0.6, 0.6, -0.6, -0.6, 0.6,                                        &
       0.8, 0.8, 0.4, 0.4, 0.8, -0.8, -0.8, -0.4, -0.4, -0.8, 0., 0.,    &
       0.4, 0.4, 0.8, 0.4, -0.4, -0.4, -0.8, -0.4, 0., 0.,               &
       0.5, 0.5, -0.5, -0.5, 0.5,                                        &
       0.5, -0.5,                                                        &
       0.2, 0.2, 0.8, 0.8, 0.2, 0.2,                                     &
       -0.2, -0.2, -0.8, -0.8, -0.2, -0.2, 0., 0.,                       &
       0., 0., 0.5, 0.5, 0.,                                             &
       0., 0., -0.5, -0.5, 0.,                                           &
       0., 0., 0.25, 0.25, 0.,                                           &
       0., 0., -0.25, -0.25, 0.,                                         &
       0.2, 0.2, -0.2, -0.2, 0.2 /

  !--- Routine body

  j_nodrift = 0
  im1 = 0

  !--- set line style to solid
  do i = 1, nel
     call jsln(1)
     it = mod(ityp(i), 20)
     if (it .eq. 0) goto 10

     j_nodrift = j_nodrift + 1
     i_nodrift(j_nodrift) = i
     if (j_nodrift .gt. 1) im1 = i_nodrift(j_nodrift-1)
     ell = ee(i) - es(i)
     if (j_nodrift .eq. 1) then
        if (es(i) .gt. hr(1))  then
           txp(1) = hr(1)
           txp(2) = es(i)
           call gvpl (2, txp, typz)
        endif
     else
        if (ee(im1) .lt. es(i))  then
           txp(1) = ee(im1)
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
10   continue
  enddo

50 continue

  call jsln(1)

  if (j_nodrift .lt. 1) then
   print*, "plot.f90 peschm j_nodrift must be > 0 and it is ", j_nodrift
   stop;
  endif

  j = i_nodrift(j_nodrift)
  if (ee(j) .lt. hr(2))  then
     txp(1) = ee(j)
     txp(2) = hr(2)
     call gvpl (2, txp, typz)
  endif
end subroutine peschm

subroutine pesopt(ierr)
  use plotfi
  use plot_bfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Stores plot options and values, checks
  !
  ! Output:  ierr  (int)     =0: OK, >0: error
  !
  ! calls the function plot_option and the routines comm_para,
  ! get_title, get_version and table_range defined in file madxn.c.
  ! calls the utility routine pesplit and the routine pegetn in this file.
  !
  ! it is called by the routine exec_plot in file madxn.c
  !----------------------------------------------------------------------*
  integer :: ierr

  integer :: i, j, k, nint, ndble
  integer :: notitle, noversi, nivaxs !, inter_setplot not used
  integer :: int_arr(szcompar), char_l(szcompar)
  integer :: plot_style(szcompar), plot_symbol(szcompar)
  double precision :: d_arr(szcompar)
  character(len=mcnam) :: sdum(1)
  character(len=szchara) :: char_a, version
  character(len=8) :: vaxisi
  character (len=120) :: warnstr

  double precision, external :: plot_option

  ierr = 0

  !--- Initialisation of variables in common peaddi
  itbv = 0
  nivvar = 0
  interf = 0
  noline = 0
  NVVAR = 0
  PROC_FLAG = 0
  IPPARM = 0
  NAXREF = 0
  NRRANG = 0

  !--- Initialisation of variables in common peaddr
  qascl = 0.0
  qlscl = 0.0
  qsscl = 0.0
  qtscl = 0.0
  HRANGE = 0.0
  VRANGE = 0.0
  hmima(1)   =  1.e20
  hmima(2)   = -1.e20
  VMIMA(1,:) =  1.e20
  VMIMA(2,:) = -1.e20

  !--- Initialisation of variables in common peaddc

  horname = ' '
  tabname = ' '
  toptitle = ' '
  SNAME = ' '

  !--- Initialisation of variables in common peotcl
  fpmach = .false.
  ptc_flag = .false.

  !--- Initialisation of local variables
  nivaxs = 0
  notitle = 0
  noversi = 0

  INT_ARR = 0
  CHAR_L = 0
  D_ARR = 0.0d0

  char_a = ' '
  sdum(1) = ' '

  !--- Routine body

  !--- get notitle
  call comm_para('notitle ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (nint .gt. 0) notitle = int_arr(1)

  !--- get noversion
  call comm_para('noversion ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (nint .gt. 0) noversi = int_arr(1)

  !--- ptc flag setting
  call comm_para('ptc ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (nint .gt. 0 .and. int_arr(1) .eq. 1) ptc_flag = .true.

  !--- if ptc flag is on look for the ptc_table
  if (ptc_flag) then
     call comm_para('ptc_table ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
     if (k .gt. 0) tabname = char_a
  else
     !--- else normal twiss treatment : any table
     call comm_para('table ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
     if (k .gt. 0) tabname = char_a
  endif

  !--- Horizontal variable  - for hor = s plot machine
  char_a = ' '
  call comm_para( 'haxis ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (k .eq. 0)  then
     call fort_warn('PLOT: ','no horizontal variable selected')
     ierr = 1
     return
  else
     horname = char_a
  endif

  if (horname .eq. 's')  itbv = 1

  !--- Prepare title
  if (notitle .eq. 0)  then
     char_a = ' '
     call comm_para('title ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
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
        if (k .gt. 0) toptitle = char_a(:k)
     endif
  endif

  qascl = plot_option('ascale ')
  qlscl = plot_option('lscale ')
  qsscl = plot_option('sscale ')
  qtscl = plot_option('rscale ')

  char_a = ' '
  call comm_para('range ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  call table_range(tabname, char_a, nrrang)
  if (nrrang(1) .eq. 0 .and. nrrang(2) .eq. 0)  then
     call fort_warn('PLOT: ','illegal RANGE')
     ierr = 1
     return
  endif

  char_a = ' '
  call comm_para('noline ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (nint .gt. 0) noline = int_arr(1)

  call comm_para('hmin ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (ndble .gt. 0) hrange(1) = d_arr(1)

  call comm_para('hmax ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (ndble .gt. 0) hrange(2) = d_arr(1)

  call comm_para('vmin ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  VRANGE(1,1:ndble) = D_ARR(1:ndble)

  call comm_para('vmax ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  VRANGE(2,1:ndble) = D_ARR(1:ndble)

  !--- Check that STYLE & SYMBOL are both non zero
  call comm_para('style ', nint, ndble, k, plot_style, d_arr, char_a, char_l)
  call comm_para('symbol ', nint, ndble, k, plot_symbol, d_arr, char_a, char_l)
  if (plot_style(1).eq.0 .and. plot_symbol(1).eq.0) then
     call fort_warn('PLOT: ','style & symbol attributes would make plot invisible. STYLE reset to 1.')
     plot_style(1) = 1
  endif
  ipparm(1,1) = plot_style(1)
  ipparm(4,1) = plot_symbol(1)

  char_a = ' '
  call comm_para('bars ', nint, ndble, k, ipparm(3,1), d_arr, char_a, char_l)

  call comm_para('colour ', nint, ndble, k, ipparm(5,1), d_arr, char_a, char_l)

  !--- if ptc_flag is on, no interpolation and check only ptc-related attributes
  if (ptc_flag .and. itbv .eq. 0) then
   call fort_warn('PLOT: ','no interpolation available with PTC.')
   return
  endif
  !--- Spline is obsolete
  call comm_para('spline ', nint,ndble,k,int_arr,d_arr, char_a,char_l)
  if (int_arr(1) .eq. 1) call fort_warn('PLOT: ', &
       'SPLINE attribute is obsolete, no action taken; use INTERPOLATE attribute instead.')

  !--- Interpolate: priority is given to the SETOPT option.
  !  If False, the option of the PLOT command will be considered.
  !  If True, the option of the PLOT command is ignored.
  ipparm(2,1) = plot_option('interpolate ')
  if (ipparm(2,1) .eq. 0) call comm_para('interpolate ', nint, ndble, k, ipparm(2,1), d_arr,char_a, char_l)

  !--- Interpolation is not possible for ptc twiss variables; reset to false
  if (ptc_flag .and. ipparm(2,1).eq.1) then
     ipparm(2,1) = 0
     call fort_warn('PLOT: ','Interpolation is not compatible with PTC_TWISS tables; INTERPOLATE ignored')
  endif

  !--- another interpolation flag is used in other subroutines
  interf = ipparm(2,1)

  !--- Continue fetching variables to be plotted
  !--- First test the variables on vaxis
  char_a = ' '
  call comm_para('vaxis ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  if (k .gt. 0)  then
     nivaxs = 1
     !--- 2014-May-05  16:34:59  ghislain: in order to avoid creating spurious ps files
     ! with basename equal to the first ignored variable (eg dy.ps), it is better to discard the whole plot
     ! and let the user fix the problem (for now at least...)
     if (k .gt. mxcurv) then
        write (warnstr,'(i2,a,i2,a)') k," variables on VAXIS would overflow max number (",&
             mxcurv,"). all variables ignored."
        call fort_warn('PLOT: ',warnstr)
        goto 110
     endif
     nivvar = min(k, mxcurv)
     call pesplit(k, char_a, char_l, slabl)
     NAXREF(1:nivvar) = 1
     !--- Do not even look further for variables on vaxis_i.
  else
     !--- Look for variables on vaxis_i
     do i = 1, 4
        write (vaxisi,'("vaxis",i1," ")') i
        char_a = ' '
        call comm_para(vaxisi, nint, ndble, k, int_arr, d_arr, char_a, char_l)
        if (k .gt. 0)  then
           !-- we need to test for overflow in number of variables before parsing the labels...
           if (nivvar+k .gt. mxcurv) then
              write (warnstr,'(i2,a,i1,a,i2,a)') k," variables on VAXIS",i,&
                   " would overflow max number(",mxcurv,"). all variables ignored."
              call fort_warn('PLOT: ',warnstr)
              goto 110
           endif
           nivaxs = nivaxs + 1
           call pesplit(k, char_a, char_l, slabl(nivvar+1))
           do j = 1, k
              nivvar = nivvar + 1
              naxref(nivvar) = i
           enddo
        endif
     enddo
  endif

  IPPARM(1,2:mxcurv) = ipparm(1,1)
  IPPARM(2,2:mxcurv) = ipparm(2,1)
  IPPARM(3,2:mxcurv) = ipparm(3,1)
  IPPARM(4,2:mxcurv) = ipparm(4,1)
  IPPARM(5,2:mxcurv) = ipparm(5,1)

110 continue

  if (nivvar .eq. 0) then
     !--- nothing to be plotted, probably because of too many variables on a single axis
     call fort_warn('PLOT: ','no vertical plot variables, plot skipped')
     ierr=1
     return
  endif

  do j = 1, nivvar
     call pegetn (0, slabl(j), itbv, proc_flag(1,j), sname(j), sdum(1))
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

end subroutine pesopt

subroutine pesplit(n_str, char_a, char_l, char_buff)
  implicit none
  !----------------------------------------------------------------------*
  !
  !   Utility routine
  !   Purpose: splits a string in several sub-strings
  !
  !--- Input
  !   n_str      number of sub-strings
  !   char_l     number of characters of each sub-string
  !   char_a     character string
  !--- Output
  !   cahr_buf   sub_strings
  !
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  integer :: n_str, char_l(*)
  character(len=*) :: char_a, char_buff(*)

  !--- type definitions of local variables
  integer :: i, k, l

  !--- Routine body
  k=0
  do i = 1, n_str
     l = char_l(i)
     char_buff(i) = char_a(k+1:k+l)
     k = k+l
  enddo

end subroutine pesplit

subroutine plginit
  use plotfi
  use plot_bfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Overall initialization
  !
  ! calls the function plot_option and the routine comm_para
  ! defined in file madxn.c.
  ! calls the function intrac defined in file madxu.c.
  ! calls the routines gxtint, gxsvar, gxasku, gxinit, gxclos
  ! defined in file gxx11.F
  !
  ! it is called by the routine plotit in this file.
  !----------------------------------------------------------------------*

  !--- type definitions of local variables
  integer :: ipseps, iset=0, nint, ndble, k, int_arr(100), char_l(100)
  double precision :: d_arr(100)
  real :: tmpval
  character(len=40) :: char_a

  !--- definitions of function primitives
  double precision, external :: plot_option
  logical, external :: intrac

  !--- Routine body
  call gxtint
  call gxsvar ('INUNIT', 5, 0., ' ')
  call gxsvar ('IOUNIT', 6, 0., ' ')
  char_a = ' '
  call comm_para('file ', nint, ndble, k, int_arr, d_arr, char_a, char_l)
  plfnam = 'madx'
  if (k .gt. 0)  plfnam = char_a(:char_l(1))

  ipseps = plot_option('post ')
  if (ipseps .eq. 0 .and. .not.intrac()) ipseps = 2

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
  tmpval = plot_option('xsize ')
  call gxsvar('XMETAF', 0, tmpval, ' ')
  tmpval = plot_option('ysize ')
  call gxsvar('YMETAF', 0, tmpval, ' ')

  !--- inhibit initial X-Window (only X11)
  call gxsvar('ITSEOP', 1, 0., ' ')
  call gxinit
  call gxclos
end subroutine plginit

subroutine pupnbl(string,ifirst,ilast)
  implicit none
  !----------------------------------------------------------------------*
  !
  !   Utility routine
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
  !----------------------------------------------------------------------*

  !--- type definition of the routine arguments
  character(len=*)  :: string
  integer, intent(OUT) :: ifirst, ilast

  !--- type definitions of local variables
  integer :: i

  !--- Output initialisation
  ifirst = 0
  ilast  = 0

  !--- Routine body

  do i=1,len(string)
     if (string(i:i).ne.' ') then
        ifirst=i
        goto 20
     endif
  enddo

  return

20 continue
  do i=len(string),1,-1
     if (string(i:i).ne.' ') then
        ilast=i
        return ! normal exit when ifirst and ilast are found
     endif
  enddo

end subroutine pupnbl
  !***********************************************************************

integer function iucomp(comp, arr, n)
  implicit none
  !----------------------------------------------------------------------*
  ! Utility function
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

  !--- type definition of the routine arguments
  integer :: comp, arr(*), n

  !--- type definitions of local variables
  integer :: j

  !--- Output initialisation
  iucomp = 0

  !--- Routine body

  do j = 1, n
     if (comp .eq. arr(j))  then
        iucomp = j
        return
     endif
  enddo
end function iucomp
