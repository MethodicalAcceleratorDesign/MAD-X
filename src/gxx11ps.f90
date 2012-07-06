subroutine gvfa(np, x, y)
  use gxx11_common
  implicit none
  integer i,icol,ierr,np
  real fx,fy,xs,ys
  !***********************************************************************
  !
  !   Purpose: Fill area plot with viewport emulation for HIGZ
  !
  !--- Input
  !    NP, x, y: as for GFA
  !   Author: H. Grote / CERN                        date: Nov. 30, 1993
  !                                              last mod: Nov. 30, 1993
  !***********************************************************************
  real x(*), y(*)

  real w(4), v(4)

  !--- set proper colour index
  call jqlctp(i)
  if (i .ne. 2)  then
     call jslctp(2)
     call jqplci(ierr, icol)
     call gxscol(icol)
  endif
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  !--- get current window
  call jqnt(1, ierr, w, v)
  !--- transform
  xs = w(2) - w(1)
  fx = vpfacx / xs
  ys = w(4) - w(3)
  fy = vpfacy / ys
  do  i = 1, np
     xvp(i) = w(1) + xs * (vploc(1) + fx * (x(i) - w(1)))
     yvp(i) = w(3) + ys * (vploc(3) + fy * (y(i) - w(3)))
  enddo
  call gfa(np, xvp, yvp)
  !--- set flag for clear permission
  iclear = 1
end subroutine gvfa
subroutine gvpl(np, x, y)
  use gxx11_common
  implicit none
  integer i,icol,ierr,iloop,ilow,n,np,nup
  real fx,fy,xs,ys
  !***********************************************************************
  !
  !   Purpose: Plot polyline and emulate viewports for HIGZ or X11
  !
  !--- Input
  !   np, x, y: as for GPL
  !   Author: H. Grote / CERN                        date: Nov. 18, 1992
  !                                           last mod: May 13, 1993
  !***********************************************************************
  real x(*), y(*)

  real w(4), v(4)
  !--- set proper colour index
  call jqlctp(i)
  if (i .ne. 2)  then
     call jslctp(2)
     call jqplci(ierr, icol)
     call gxscol(icol)
  endif
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  !--- get current window
  call jqnt(1, ierr, w, v)
  !--- transform
  xs = w(2) - w(1)
  fx = vpfacx / xs
  ys = w(4) - w(3)
  fy = vpfacy / ys
  do  iloop=1, np, madim2
     nup = min(np, iloop + madim2 - 1)
     ilow = max(1, iloop - 1)
     n = 0
     do  i = ilow, nup
        n = n + 1
        xvp(n) = w(1) + xs * (vploc(1) + fx * (x(i) - w(1)))
        yvp(n) = w(3) + ys * (vploc(3) + fy * (y(i) - w(3)))
     enddo
     call gpl(n, xvp, yvp)
  enddo
  !--- set flag for clear permission
  iclear = 1
end subroutine gvpl
subroutine gvpm(np, x, y)
  use gxx11_common
  implicit none
  integer i,icol,ierr,iloop,n,np,nup
  real fx,fy,xs,ys
  !***********************************************************************
  !
  !   Purpose: Plot marker symbol and emulate viewports for HIGZ or X11
  !
  !--- Input
  !   np, x, y: as for GPM
  !   Author: H. Grote / CERN                        date: Nov. 18, 1992
  !                                           last mod: May 13, 1993
  !***********************************************************************
  real x(*), y(*)

  real w(4), v(4)

  !--- set proper colour index
  call jqlctp(i)
  if (i .ne. 2)  then
     call jslctp(2)
     call jqpmci(ierr, icol)
     call gxscol(icol)
  endif
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  !--- get current window
  call jqnt(1, ierr, w, v)
  !--- transform
  xs = w(2) - w(1)
  fx = vpfacx / xs
  ys = w(4) - w(3)
  fy = vpfacy / ys
  do  iloop=1, np, madim2
     nup = min(np, iloop + madim2 - 1)
     n = 0
     do  i = iloop, nup
        n = n + 1
        xvp(n) = w(1) + xs * (vploc(1) + fx * (x(i) - w(1)))
        yvp(n) = w(3) + ys * (vploc(3) + fy * (y(i) - w(3)))
     enddo
     call gxwpm(n, xvp, yvp)
  enddo
  !--- set flag for clear permission
  iclear = 1
end subroutine gvpm
subroutine gvtx(x, y, sss)
  use gxx11_common
  implicit none
  integer i,icol,ierr
  real chh,chux,chuy,fx,fy,hfac,x,xs,y,ys
  !***********************************************************************
  !
  !   Purpose: Plot text and emulate viewports for HIGZ or X11
  !
  !--- Input
  !   x, y, s: as for GTX
  !   Author: H. Grote / CERN                        date: Nov. 18, 1992
  !                                           last mod: May 13, 1993
  !***********************************************************************
  character(*) sss

  real w(4), v(4)

  !--- set proper colour index
  call jqlctp(i)
  if (i .ne. 2)  then
     call jslctp(2)
     call jqtxci(ierr, icol)
     call gxscol(icol)
  endif
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  !--- get current window
  call jqnt(1, ierr, w, v)
  !--- get current character height and text orientation
  call jqchh(ierr, chh)
  call jqchup(ierr, chux, chuy)
  !--- transform
  xs = w(2) - w(1)
  fx = vpfacx / xs
  ys = w(4) - w(3)
  fy = vpfacy / ys
  if (chux .eq. 0.)  then
     hfac = fy
  else
     hfac = fx
  endif
  call jschh(hfac * chh)
  xvp(1) = w(1) + xs * (vploc(1) + fx * (x - w(1)))
  yvp(1) = w(3) + ys * (vploc(3) + fy * (y - w(3)))
  call gtx(xvp(1), yvp(1), sss)
  call jschh(chh)
  !--- set flag for clear permission
  iclear = 1
end subroutine gvtx
subroutine gxarng(nopt,rmini,rmaxi,rmin,rmax,nint)
  implicit none
  integer nint,nopt
  real rmax,rmaxi,rmin,rmini
  !***********************************************************************
  !
  !   Purpose: calculates axis ranges
  !
  !--- Input
  !   nopt         =0: normal
  !              =1: start or terminate axis at 0. if possible
  !              =2: centre axis around 0.
  !   rmini        minimum (x or y) value to consider
  !   rmaxi        maximum (x or y) value to consider
  !--- Output
  !   rmin         lower end of axis
  !   rmax         upper end of axis
  !   nint         no. of intervals as returned by GXSCAL
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Aug. 8, 1988
  !
  !***********************************************************************
  nint=10
  if(rmini.ge.rmaxi) then
     rmin=rmini
     rmax=rmini+1.
  else
     if(nopt.eq.0) then
        call gxscal(rmini,rmaxi,rmin,rmax,nint)
     elseif(nopt.eq.1) then
        !--- start or terminate at 0.
        if(rmini.gt.0.) then
           call gxscal(0.,rmaxi,rmin,rmax,nint)
        elseif(rmaxi.lt.0.) then
           call gxscal(0.,-rmini,rmin,rmax,nint)
           rmin=-rmax
           rmax=0.
        else
           call gxscal(rmini,rmaxi,rmin,rmax,nint)
        endif
     else
        call gxscal(0.,max(abs(rmini),abs(rmaxi)),rmin,rmax,nint)
        rmin=-rmax
     endif
  endif
end subroutine gxarng
subroutine gxasku
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: asks user interactively for plot options
  !
  !   must be called before GXINIT if at all
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************
  logical intrac

  call gxundf
  if (intrac())  then
     call gxask1
     call gxask2
  endif
end subroutine gxasku
subroutine gxask1
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: asks user interactively plot window
  !
  !   called by GXASKU
  !
  !   Author: H. Grote / CERN                        date: May 13, 1993
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************

  !--- Input and Output unit definition
  if(lnunit.ne.lundef)  call gxsvar('INUNIT',miunit,0.,' ')
  if(lounit.ne.lundef)  call gxsvar('IOUNIT',mounit,0.,' ')
  itermt=0
  interm=0
  lnterm=lundef
  ltermt=lundef
  goto 999
30 continue
  write(iounit,*) ' Error on Input, stop.'
  stop
10000 format(//' GX (X11 based) plot package initialization'/)
10010 format(/' Do you want to plot on your terminal ? (<CR> = yes>:'/)
999 end subroutine gxask1
subroutine gxask2
  use gxx11_common
  implicit none
  integer iamx,iamy,ierr,ifirst,ilast
  real xax,yax
  !***********************************************************************
  !
  !   Purpose: asks user interactively for Postscript file
  !
  !   called by GXASKU
  !
  !   Author: H. Grote / CERN                        date: May 13, 1993
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************
  character sline*80
  character(60) gxform, sform
  logical affirm

  if (lpseps .ne. lundef)  then
     write(iounit,10040)
     call gxrdtx(inunit,sline,ierr)
     if(ierr.ne.0) goto 30
     call gxpnbl(sline, ifirst, ilast)
     if(ifirst .eq. 0)  then
        ipseps = 0
     elseif (index('0123456789', sline(:ifirst)) .eq. 0)  then
        ipseps = 0
     else
        sform = '(I$$)'
        write (sform(3:4), '(I2.2)')  ilast
        read (sline, sform) ipseps
        if (ipseps .ge. 1 .and. ipseps .le. 2)  then
           call jswks(1)
           inmeta = mtmeta
           if(lmetop .ne. lundef)  then
              if (lmetnm .ne. lundef)  then
                 smetnm = 'gxx11'
                 call gxpnbl(smetnm, ifirst, ilast)
                 write(iounit,10100) smetnm(ifirst:ilast),               &
                      smetnm(ifirst:ilast)
                 call gxrdtx(inunit,sline,ierr)
                 if(sline(1:1).ne.' ')                                   &
                      call gxsvar('SMETNM',0,0.,sline)
                 lmetnm = lundef
              endif
           endif
        else
           ipseps = 0
        endif
     endif
  endif
  lpseps = lundef
  lnmeta=lundef
  if (ipseps .eq. 0)  then
     inmeta = 0
  else
     inmeta = mtmeta
     !--- paper size (only if not set already)
     if(lmetax.ne.lundef.or.lmetay.ne.lundef)  then
        iamx = mxsize
        iamy = mysize
        xax = iamx
        yax = iamy
        write(iounit,10050) iamx, iamy
        call gxrdtx(inunit,sline,ierr)
        if(ierr.ne.0) goto 30
        if(affirm(sline(1:1))) then
           write(iounit,10070)
           call gxrdtx(inunit,sline,ierr)
           if(ierr.ne.0) goto 30
           sform=gxform(sline)
           if(index(sform,'I').ne.0)  then
              read(sline,sform)  iamx
              xax = iamx
           else
              read(sline,sform) xax
           endif
           write(iounit,10080)
           call gxrdtx(inunit,sline,ierr)
           if(ierr.ne.0) goto 30
           sform=gxform(sline)
           if(index(sform,'I').ne.0)  then
              read(sline,sform)  iamy
              yax = iamy
           else
              read(sline,sform) yax
           endif
        endif
        call gxsvar('XMETAF', 0, xax, ' ')
        call gxsvar('YMETAF', 0, yax, ' ')
     endif
  endif
  goto 999
30 continue
  write(iounit,*) ' Error on Input, stop.'
  stop
10040 format(/' Do you want to write a .ps file <1>, .eps files <2>,',  &
       ' or none <CR>:'/)
10050 format                                                            &
       (/' specify bounding box size (default:',i3,                      &
       '(x) by', i3, '(y) cm)?'/' (<CR>=no):'/)
10070 format(/' enter bounding box x size in cm:'/)
10080 format(/' enter bounding box y size in cm:'/)
10100 format(/' enter postscript or eps file name (leading part)' /     &
       ' (<CR> gives "',a,'.ps" resp. "',a,'nn.eps"):'/)
999 end subroutine gxask2
subroutine gxaxis(type,axlow,axup,axpos,ipos,fmt,textin,sepchr,iparm,ierr)
  implicit none
  integer i,ia,ialow,iaup,iaxort,ie,ierr,ietick,ifircl,ifirst,ifont,&
       ifs,ilabl,ilast,ilbort,ils,impfl,in,intrep,intv,ipos,irf,iscloc,  &
       islbl,islpc,isp,ispchl,isradc,itext,itick,ival,k1,k2,l,l1,l2,     &
       naxal,nchct,nlines
  real a1,a1b,a2,alp,amxx,axlow,axpos,axup,chhigh,chwdth,cthigh,    &
       cuhigh,diff,diffe,diffn,fcw,fwc,hxf,pfact,ptick,sgspac,sk,space,  &
       sphlin,spmlog,spwlin,spwlog,tickl,wbused,wsused
  !***********************************************************************
  !
  !   Purpose: plots an axis with tick marks, numbers, and title
  !
  !-- Input
  !   type        'X' for an x-axis, 'Y' for a y-axis
  !   axlow       lower end of axis in current world coords.
  !   axup        upper end of axis in current world coords.
  !   axpos       position of axis in the other coordinate
  !   ipos        =0: AXPOS value given in normalized dev. coord. [0.,1.]
  !            =1: AXPOS value given in current world coords.
  !   fmt         (floating point) format for axis labels (=numbers)
  !            including the brackets, e.g. '(F6.3)'. If blank,
  !            a reasonable default is used.
  !   textin      axis title (trailing blanks will be suppressed)
  !   sepchr      will start a new line when encountered in TEXT
  !   iparm       axis parameters:
  !   1        0 = linear scaling
  !         1 = logarithmic scale
  !+++ if the user scale makes no sense, scaling is switched to automatic
  !
  !   2        if = 0, no tick marks. If < 0, the number of intervals
  !         will be chosen automatically.
  !         if linear scaling and > 0, no. of major tick mark intervals
  !         (labels are only written at major tick marks).
  !         if log. scaling and <> 0, major tick marks at the powers of
  !         ten in the scale, i.e. at all integer values.
  !
  !   3        0 no labels (scale numbers), 1 hor. labels, 2 vertical
  !
  !   4        odd for tick mark below (x-axis) or at left (y axis),
  !         even for above resp. at right.
  !         if = 0, no ticks
  !
  !   5        as 4, but for labels (=scale numbers)
  !
  !   6        as 4, but for the axis text (title)
  !         the text is written horizontally for x-, vertically for y-axes
  !
  !   7        character height in normalized pixels
  !
  !   8        tick mark length in normalized pixels
  !         a normalized pixel is defined as follows: imagine your
  !         default (square) screen (device) area devided into
  !         1000 x 1000 pixels, i.e. one pixel is 0.001 x 0.001 in NDC
  !
  !   9        Linear scaling: no. of extra intervals with half-size ticks
  !         between main ticks
  !         log. scaling: if > 0 : flag that extra ticks are to be
  !         plotted at the positions log(2), log(3),..., log(9)
  !
  !   10       =1 : adjust axis titles at left, =2 centre, =3 adjust
  !         them at right. If the string '<#>' is found inside a line,
  !         the text to the left of it will always be left adjusted,
  !         the text to the right of it right adjusted.
  !
  !   11       text font for axis labels, default = 1 (see GKS for details)
  !
  !   12       if 0 (default) no minor labels at minor ticks
  !         for log. scale, if > 0 yes
  !
  !   13       if > 0 and not > 1000: character height for axis text, else
  !         parameter 7 is used
  !
  !-- Output
  !   ierr        =0: everything OK
  !            =-1: AXLOW.GE.AXUP
  !            else the corresponding GKS error
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Feb. 3, 1993
  !
  !***********************************************************************
  character(*)  textin,fmt,type,sepchr
  integer iparm(*)
  character   stext*240,text*240,sltext*40,                         &
       fmtloc*60,stsep*1,slog(9)*1
  integer isave(20),ihoral(2,2,2),iveral(2,2,2)
  real atext(4),tick(4),alabl(4),cnt(4),rsave(20),tetick(4),        &
       alogv(10)
  logical xaxis,linscl, labflg
  save alogv, ifircl, ihoral, iveral, slog
  !--- horizontal and vertical alignment as function of orientation (I),
  !   position above or below etc. (J) and axis (x or y) (K), e.g.
  !   iveral(1,2,1) vertical alignment for horizontal labels above an x axis
  data ihoral/2,3,2,1,3,2,1,2/
  data iveral/1,3,5,3,3,5,3,1/
  data slog/'1','2','3','4','5','6','7','8','9'/
  data ifircl/0/
  stext = ' '
  text = ' '
  sltext = ' '
  fmtloc = ' '
  stsep = ' '
  do i = 1, 20
     isave(i) = 0
  enddo
  do i = 1, 9
     slog(i) = ' '
  enddo
  do i=1,4
     atext(i)=0.
  enddo
  do i=1,20
     rsave(i)=0.
  enddo
  space=0.
  if(ifircl.eq.0)  then
     !--- set logarithms
     do i=1,10
        alogv(i)=log10(float(i))
     enddo
     ifircl=1
  endif
  !
  !--- check for reasonable axis range
  !
  if(axlow.ge.axup)  then
     ierr=-1
     goto 999
  endif
  !
  !--- get current user settings and keep them
  call gxsave(isave,rsave,ierr)
  if(ierr.ne.0) goto 999
  !
  !--- get Input parameters
  !
  stsep=sepchr
  xaxis=type(:1).eq.'X'
  if(xaxis)  then
     iaxort=1
  else
     iaxort=2
  endif
  linscl=iparm(1).eq.0
  intv=iparm(2)
  if(intv.lt.0)  then
     !--- choose number of intervals automatically if linear
     if(linscl)  then
        call gxdint(axlow,axup,intv)
     else
        intv=1
     endif
  endif
  labflg = .false.
  fmtloc=fmt
  call gxival(fmtloc,ival)
  if(ival.eq.0 .and. linscl)  then
     !--- use reasonable default as format
     call gxdfmt(axlow,axup,intv,ival,iscloc,fmtloc)
  else
     iscloc = 0
  endif
  ilbort=iparm(3)
  if(ilbort.gt.0)  ilbort=mod(ilbort-1,2)+1
  itick=iparm(4)
  ietick=iparm(9)
  if(itick.gt.0)  itick=mod(itick-1,2)+1
  ilabl=iparm(5)
  if(ilabl.gt.0)  ilabl=mod(ilabl-1,2)+1
  itext=iparm(6)
  if(itext.gt.0)  itext=mod(itext-1,2)+1
  cuhigh=.001*iparm(7)
  if(iparm(13).gt.0.and.iparm(13).le.1000)  then
     cthigh=.001*iparm(13)
  else
     cthigh=cuhigh
  endif
  tickl =.001*iparm(8)
  if(itick.eq.0)  tickl=0.
  naxal=iparm(10)
  ifont=iparm(11)
  islpc=iparm(12)
  if(xaxis)  then
     k1=1
     k2=3
  else
     k1=3
     k2=1
     !--- apply expansion factor to tick marks if y axis
     call gxqrvp(hxf)
     tickl=hxf*tickl
  endif
  !
  !--- transform into normalized window
  !
  if(isave(1).ne.0)  then
     !--- rsave(1...4) contains the window
     !   fwc converts a length along the axis from world to NDC
     !   fcw does the inverse of FWC
     fcw=(rsave(k1+1)-rsave(k1))
     fwc=1./fcw
     cnt(k1)=fwc*(axlow-rsave(k1))
     cnt(k1+1)=fwc*(axup-rsave(k1))
     if(ipos.gt.0)  then
        cnt(k2)=(axpos-rsave(k2))/(rsave(k2+1)-rsave(k2))
     else
        cnt(k2)=axpos
     endif
     call jswn(isave(1),0.,1.,0.,1.)
     !     call jselnt(isave(1))
  else
     fcw=1.
     fwc=1.
     cnt(k1)=axlow
     cnt(k1+1)=axup
     cnt(k2)=axpos
  endif
  cnt(k2+1)=cnt(k2)
  !
  !--- set line style
  !
  call gxspmt
  !
  !   set font and precision
  !
  call jstxfp(ifont,2)
  !
  !--- plot a line for the axis
  !
  call gvpl(2,cnt(1),cnt(3))
  !
  !--- plot tick marks, labels, and title
  !
  tick(k2)=cnt(k2)
  tick(k2+1)=cnt(k2+1)
  tetick(k2)=cnt(k2)
  tetick(k2+1)=cnt(k2+1)
  alabl(k2)=cnt(k2)
  alabl(k2+1)=cnt(k2+1)
  if(intv.gt.0)  then
     if(itick.eq.2) then
        !--- plot tick marks above x-axis, or to the right of y-axis
        tick(k2+1)=cnt(k2)+tickl
        tetick(k2+1)=cnt(k2)+.5*tickl
     else
        tick(k2+1)=cnt(k2)-tickl
        tetick(k2+1)=cnt(k2)-.5*tickl
     endif
     if(linscl)  then
        diff=(axup-axlow)/intv
        diffn=(cnt(k1+1)-cnt(k1))/intv
        diffe=diffn/max(ietick,1)
        if(tickl.gt.0.) then
           do i=0,intv
              !--- tick marks
              tick(k1)=cnt(k1)+diffn*i
              tick(k1+1)=tick(k1)
              call gvpl(2,tick(1),tick(3))
              !--- extra ticks
              if(i.lt.intv) then
                 do ie=1,ietick-1
                    tetick(k1)=tick(k1)+ie*diffe
                    tetick(k1+1)=tetick(k1)
                    call gvpl(2,tetick(1),tetick(3))
                 enddo
              endif
           enddo
        endif
     elseif(tickl.gt.0.) then
        !--- log scale
        ialow=axlow
        if(axlow.lt.float(ialow))  ialow=ialow-1
        ialow=sign(min(abs(ialow),99),ialow)
        iaup=axup
        if(axup.lt.float(iaup))  iaup=iaup-1
        iaup=sign(min(abs(iaup),99),iaup)
        ia=ialow
        !--- start loop
40      continue
        do i=1,9
           alp=ia+alogv(i)
           if(alp.gt.axup) goto 60
           if(alp.ge.axlow) then
              if(i.eq.1.or.ietick.ne.0)  then
                 tick(k1)=(alp-rsave(k1))*fwc
                 tick(k1+1)=tick(k1)
                 call gvpl(2,tick(1),tick(3))
              endif
           endif
        enddo
        ia=ia+1
        goto 40
        !--- end loop
60      continue
     endif
     !--- labels
     if(ilbort.eq.0.or.ilabl.eq.0.or.cuhigh.eq.0.) then
        alabl(k2)=tick(k2)
        alabl(k2+1)=alabl(k2)
     else
        labflg = .true.
        !--- set correct character height (viewports !), get width
        call gxschf(1,ilbort,cuhigh,chhigh,chwdth)
        space=.5*chhigh
        if(iaxort.eq.ilbort) then
           sphlin=chhigh
           spwlin=chwdth*ival
           spwlog=2.5*chwdth
           spmlog=chwdth
        else
           spwlin=chhigh
           sphlin=chwdth*ival
           spwlog=chhigh
           spmlog=chhigh
        endif
        if(ilabl.eq.2) then
           !--- plot labels above x-axis, or to the right of y-axis
           alabl(k2)=max(tick(k2),tick(k2+1))+space
           alabl(k2+1)=alabl(k2)+sphlin
        else
           alabl(k2)=min(tick(k2),tick(k2+1))-space
           alabl(k2+1)=alabl(k2)-sphlin
        endif
        !--- set text alignment
        call jstxal(ihoral(ilbort,ilabl,iaxort),                      &
             iveral(ilbort,ilabl,iaxort))
        if(linscl)  then
           !--- linear scale
           intrep=spwlin/diffn+.99999
           amxx=max(abs(axlow),abs(axup))
           pfact = 10.**(-iscloc)
           do i=0,intv,intrep
              !--- centered figures
              ptick=axlow+diff*i
              alabl(k1)=cnt(k1)+diffn*i
              if(abs(ptick)/amxx.lt.1.e-5)  then
                 sltext='0.0'
              else
                 write(sltext,fmtloc) pfact * ptick
              endif
              call gxpnbl(sltext,ifirst,ilast)
              call gxstx(alabl(1),alabl(3),sltext(ifirst:ilast))
           enddo
        else
           !---  log scale
           !   impfl = 0 if power of ten to be plotted with first minor label,
           !   else > 0
           impfl=iaup-ialow
           !--- islbl = 0 if no secondary label to be plotted else 1
           if(ietick.eq.0.or.islpc.eq.0)  then
              islbl=0
           else
              islbl=1
           endif
           !--- wbused = half space used by major label in WC
           wbused=spwlog*fcw
           !--- wsused = half space used by minor label in WC
           wsused=2.*spmlog*fcw
           ia=ialow
           !--- interval which is free for major labels
           a1b=ialow-1000.
           !--- start loop
80         continue
           !--- interval which is free for minor labels
           a1=ia+wbused
           a2=ia+1.-wbused
           do i=1,9
              !--- label position in WC
              alp=ia+alogv(i)
              if(alp.gt.axup) goto 100
              if(alp.ge.axlow)  then
                 !--- label position in NDC
                 alabl(k1)=(alp-rsave(k1))*fwc
                 if(i.eq.1.and.alp.ge.a1b+wbused) then
                    !--- major label (at integers)
                    call gxppow(alabl,ia)
                    a1b=alp+wbused
                 else
                    if(islbl.ne.0.and.alp.ge.a1+wsused .and.alp.le.       &
                         a2-wsused) then
                       !--- minor label (at 2, 3, ..., 9)
                       call gxstx(alabl(1),alabl(3),slog(i))
                       a1=a1+wsused
                    endif
                    if(impfl.eq.0) then
                       impfl=1
                       !--- plot power of ten between first and second minor label
                       alabl(k1)=alabl(k1)                                 &
                            +.5*(min(axup,alp+alogv(i+1))-alp)*fwc
                       call gxppow(alabl,ia)
                       a1=a1+spwlog
                    endif
                 endif
              endif
           enddo
           ia=ia+1
           goto 80
100        continue
        endif
     endif
  endif
  !--- axis title
  call gxpnbl(textin, ifirst, ilast)
  if (iscloc .ne. 0 .and. labflg)  then
     ispchl = max(1, index(textin(:ilast), stsep))
     isradc = index(textin(ispchl:ilast), '<#>')
     if (isradc .eq. 0)  then
        sltext = '<#>[*10**($$$)]'
     else
        sltext = '   [*10**($$$)]'
     endif
     write(sltext(11:13), '(i3)') iscloc
     if (ifirst .eq. 0)  then
        text = sltext
     else
        text = textin(ifirst:ilast) // sltext
     endif
  else
     text = textin
  endif
  call gxpnbl(text,ifirst,ilast)
  if(ifirst.ne.0.and.itext.ne.0.and.cuhigh.gt.0.)  then
     !--- set correct character height (viewports !), get width
     call gxschf(1,iaxort,cthigh,chhigh,chwdth)
     if(xaxis)  then
        sgspac=-1.75 * chhigh
     else
        sgspac=1.5 * chhigh
     endif
     if(naxal.lt.1.or.naxal.gt.3)  naxal=2
     if(naxal.eq.1)  then
        atext(k1)=cnt(k1)
     elseif(naxal.eq.2)  then
        atext(k1)=.5*(cnt(k1)+cnt(k1+1))
     else
        atext(k1)=cnt(k1+1)
     endif
     if(itext.eq.2)  then
        atext(k2)=max(tick(k2),tick(k2+1),alabl(k2),alabl(k2+1))      &
             +space
     else
        atext(k2)=min(tick(k2),tick(k2+1),alabl(k2),alabl(k2+1))      &
             - 2.5 * space
     endif
     !--- get number of separation characters
     call gxchct(text(ifirst:ilast),stsep,nchct)
     nlines=nchct+1
     if(xaxis)  then
        if(itext.eq.2)  atext(k2)=atext(k2)+1.5*nlines*chhigh
     else
        if(itext.eq.1)  atext(k2)=atext(k2)-1.5*nlines*chhigh
     endif
     !--- write line by line
     irf=1
110  continue
     isp=index(text(irf:ilast),stsep)
     if(isp.eq.0)  then
        isp=ilast
     else
        isp=irf+isp-2
     endif
     if(isp.ge.irf)  then
        stext=text(irf:isp)
        l=isp+1-irf
        call gxpnbl(stext(:l),ifs,ils)
        in=index(stext(:l),'<#>')
        if(in.eq.0)  then
           call jstxal(naxal,1)
           if(ifs.gt.0)  call gxtx(atext(1),atext(3),stext(:ils))
        else
           !--- split line into left and right adjusted part
           l1=in-1
           l2=in+3
           if(l1.gt.0)  then
              call jstxal(1,1)
              sk=atext(k1)
              atext(k1)=cnt(k1)
              call gxtx(atext(1),atext(3),stext(:l1))
              atext(k1)=sk
           endif
           if(l2.le.ils)  then
              call jstxal(3,1)
              sk=atext(k1)
              atext(k1)=cnt(k1+1)
              call gxtx(atext(1),atext(3),stext(l2:ils))
              atext(k1)=sk
           endif
        endif
     endif
     irf=isp+2
     if(irf.le.ilast)  then
        atext(k2)=atext(k2)+sgspac
        goto 110
     endif
  endif
  !
  !--- restore user settings
  !
  call gxrest(isave,rsave)
999 end subroutine gxaxis
subroutine gxchct(stext,sch,n)
  implicit none
  integer ilast,irf,isp,n
  !***********************************************************************
  !
  !   Purpose: counts number of given characters in a string
  !
  !--- Input
  !   stext      string
  !   sch        special character
  !--- Output
  !   n          number of occurences
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************
  character stext*(*),sch*1
  n=0
  irf=1
  ilast=len(stext)
10 continue
  isp=index(stext(irf:),sch)
  if(isp.gt.0) then
     n=n+1
     irf=irf+isp
     if(irf.le.ilast) goto 10
  endif
end subroutine gxchct
subroutine gxclos
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: close terminal workstation
  !
  !   Author: H. Grote / CERN                        date: Feb. 26, 1988
  !                                           last mod: Feb. 26, 1988
  !
  !***********************************************************************
  call gxundf
  if(lacttm.eq.lundef)  then
     if(interm.gt.0)  then
        call wdawk(interm)
        call wclwk(interm)
        lacttm=0
     endif
  endif
end subroutine gxclos
subroutine gxclrw
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: clears open workstations, sets new picture name
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: March 2, 1988
  !
  !***********************************************************************

  call gxundf
  if(iclear .ne. 0) then
     call gclrwk(0, 1)
     iclear = 0
     if (iepsop .eq. 2)  then
        call gxopps(0, 0)
        call wclwk(inmeta)
        close(imetps)
        iepsop = -iepsop
     endif
  endif
end subroutine gxclrw
subroutine gxcubi(npoint,xx,yy,yy1d,yy2d,ierror)
  implicit none
  integer i,ierror,n,npoint
  real xx,yy,yy1d,yy2d
  !***********************************************************************
  !
  !   new (internal) double precision version 29.1.88
  !
  !   calculates a third order natural spline through function values
  !   yy(i), i=1,...,NPOINT, at knots XX(I) (XX(I+1) > XX(I))
  !
  !---  Input
  !   npoint  number of knots and dimension of XX, YY, YY1D, YY2D
  !          minimum value = 3
  !   xx(i)   x values,  XX(I) < XX(I+1)  for all I
  !   yy(i)   function values
  !
  !---  Output
  !   yy1d(i) first derivative of third order pol. in interval I,
  !          at point XX(I)
  !   yy2d(i) second derivative
  !   ierror  0 if everything OK, else number of the first x value found
  !          that is smaller or equal to the previous one, or -1 if NPOINT < 3
  !
  !++++++ warning: all first and second derivatives set to zero if the
  !    condition XX(I+1) > XX(I) is not fulfilled for I = [1,n-1] (IERROR > 0)
  !
  !--- remark: very near to routine SPLIN3 in CERN library, E209
  !
  !   Author hG                         13.11.86   last mod. 29.1.88
  !
  !***********************************************************************
  dimension xx(*),yy(*),yy1d(*),yy2d(*)
  double precision zero,half,one,three,third,dfac,dx1,dx2,dy1,dy2,  &
       dd,dyx1,dyx2,divdif,alf,bet
  save zero,half,one,three
  data zero,half,one,three/0.d0,0.5d0,1.d0,3.d0/
  third=one/three
  ierror=0
  n=npoint
  yy2d(1)=0.
  yy2d(n)=0.
  yy1d(1)=0.
  !
  !--- method: see long write-up of E209. Basically, the second
  !   derivatives are found first from the solution of N-2 equations.
  !   the first and last second order derivative are set to zero (hence
  !   natural spline). The equations form a three-diagonal matrix.
  !   in a first pass, all but the latest unknown are
  !   eliminated, in the second pass all are then calculated by going
  !   backwards.
  !--- yy1d serves temporarily as intermediate storage for the factors
  !   in the first pass of this process.
  !
  if(n.eq.3)  then
     !--- only three points - direct solution
     dx1=xx(2)-xx(1)
     if(dx1.le.zero)  then
        ierror=2
        goto 40
     endif
     dx2=xx(3)-xx(2)
     if(dx2.le.zero)  then
        ierror=3
        goto 40
     endif
     dy1=yy(2)-yy(1)
     dy2=yy(3)-yy(2)
     dd=one/(dx1+dx2)
     dyx1=dy1/dx1
     dyx2=dy2/dx2
     yy2d(2)=three*dd*(dyx2-dyx1)
     yy1d(1)=dyx1-dx1*yy2d(2)*half*third
     yy1d(3)=dyx2+dx2*yy2d(2)*half*third
     yy1d(2)=yy1d(1)+half*dx1*yy2d(2)
  elseif(npoint.gt.3)  then
     dx2=xx(2)-xx(1)
     if(dx2.le.zero)  then
        ierror=2
        goto 40
     endif
     dyx2=(yy(2)-yy(1))/dx2
     do i=2,n-1
        dx1=dx2
        dx2=xx(i+1)-xx(i)
        if(dx2.le.zero)  then
           ierror=i+1
           goto 40
        endif
        dyx1=dyx2
        dyx2=(yy(i+1)-yy(i))/dx2
        dd=one/(dx1+dx2)
        divdif=dd*(dyx2-dyx1)
        alf=half*dd*dx1
        bet=half-alf
        !
        !--- the following IF is only necessary for splines other than natural
        !
        if(i.eq.2)  then
           divdif=divdif-third*alf*yy2d(1)
        elseif(i.eq.n-1)  then
           divdif=divdif-third*bet*yy2d(n)
        endif
        dfac=one/(one+alf*yy1d(i-1))
        yy1d(i)=-dfac*bet
        yy2d(i)=dfac*(three*divdif-alf*yy2d(i-1))
     enddo
     !
     !--- now the last unknown derivative, YY2D(N-1), has been calculated.
     !   the others follow from going up the system.
     do i=n-2,1,-1
        dd=yy1d(i)
        yy2d(i)=dd*yy2d(i+1)+yy2d(i)
     enddo
     !
     !--- now the first derivatives from a direct equation (not the one
     !   given in the E-209 writeup - it can be simplified)
     !
     do i=1,n-1
        dx2=xx(i+1)-xx(i)
        dyx2=(yy(i+1)-yy(i))/dx2
        yy1d(i)=dyx2-dx2*third*(yy2d(i)+half*yy2d(i+1))
     enddo
     dx2=xx(n)-xx(n-1)
     dyx2=(yy(n)-yy(n-1))/dx2
     yy1d(n)=dyx2+dx2*third*(yy2d(n)+half*yy2d(n-1))
  else
     !
     !--- n < 3: error exit as well
     ierror=-1
     goto 40
  endif
  goto 999
40 continue
  !
  !--- error condition: all first and second derivatives to zero
  !
  do i=1,n
     yy1d(i)=0.
     yy2d(i)=0.
  enddo
999 end subroutine gxcubi
function gxcubv(x,npoint,xx,yy,yy1d,yy2d)
  implicit none
  integer i,npoint
  real gxcubv,x,xx,yy,yy1d,yy2d
  !***********************************************************************
  !
  !   new (internal) double precision version
  !
  !   calculates the value of a third order spline at X. The routine
  !   gxcubi must be called beforehand.
  !
  !---  Input
  !   x       abscissa value. For X outside [XX(1),XX(npoint)], a linear
  !          extrapolation is performed.
  !   npoint  number of knots and dimension of XX, YY, YY1D, YY2D
  !          minimum value = 3
  !   xx(i)   x values,  XX(I) < XX(I+1)  for all I
  !   yy(i)   function values
  !   yy1d(i) first derivative of third order pol. in interval I,
  !          at point XX(I), from GXCUBI
  !   yy2d(i) second derivative, from GXCUBI
  !
  !   Author hG                         13.11.86   last mod. 29.1.88
  !
  !***********************************************************************
  dimension xx(*),yy(*),yy1d(*),yy2d(*)
  double precision half,dx,h2,h3,h4,h5,h6
  save half
  data half/0.5d0/
10 continue
  if(x.le.xx(1))  then
     dx=x-xx(1)
     gxcubv=yy(1)+yy1d(1)*dx
  elseif(x.ge.xx(npoint))  then
     dx=x-xx(npoint)
     gxcubv=yy(npoint)+yy1d(npoint)*dx
  else
     do i=1,npoint-1
        if(x.lt.xx(i+1)) goto 30
     enddo
30   continue
     dx=x-xx(i)
     h4=yy2d(i)
     h6=(yy2d(i+1)-h4)/(xx(i+1)-xx(i))
     h2=half*h4
     h5=half*h6
     h3=h5/3.d0
     gxcubv=((h3*dx+h2)*dx+yy1d(i))*dx+yy(i)
     !   first derivative in X  = (H5*DX+H4)*DX+YY1D(I)
     !   second    "      "  "  = H6*DX+H4
  endif
end function gxcubv
subroutine gxcrv1(nset,nptval,ipxval,ipyval,icvref,xval,yval,window,actwin,ierr)
  use gxx11_common
  implicit none
  integer ibar,ic,ierr,isplin,isym,j,kset,line,nset,mark
  real dum1,dum2,fsx,fsy,xs,ys
  !***********************************************************************
  !
  !   Purpose: plots curves into an existing frame, clips
  !
  !--- Input
  !   nset       number of curves (=ordered sets of (x,y) pairs) to plot
  !   nptval(i)  number of points ((x,y) pairs) in set I
  !   ipxval(i)  first x value of set I in array XVAL
  !   ipyval(i)  first y value of set I in array YVAL
  !   icvref(i)  number of the parameter set to be used for curve I. This value
  !            will be forced into [1,MAXSET].
  !            The x and y axis reference numbers of set I will be taken
  !            from this parameter set. All x and y axes with the
  !            corresponding reference numbers will be (scaled if automatic)
  !            and plotted together with set I.
  !            If no x resp. y axis exists with the reference number
  !            from the parameter set, the curve will be plotted with
  !            automatic scaling, but without x resp. y axis.
  !   xval       array containing the x values for all sets
  !   yval       array containing the y values for all sets
  !   window(j,I) GKS window (J=1...4) to be used with curve I. These values
  !            can be obtained from routines GXFRAM or GXFRM1
  !   actwin(j,I) active window (J=1...4) to clip curve I. These values
  !            can be obtained from routine GXFRM1
  !--- Output
  !   ierr       0 if everything OK, else GKS error, or
  !            1 : GXINIT not called (initialization)
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: Dec. 9, 1988
  !
  !***********************************************************************

  integer nptval(*),ipxval(*),ipyval(*),icvref(*)
  real xval(*),yval(*),window(4,*),actwin(4,*)
  real wn(4),ac(4),xx(2),yy(2),rsave(20)
  integer isave(20)
  character sss*1
  !
  !--- get current user settings and keep them
  call gxsave(isave,rsave,ierr)
  if(ierr.ne.0)  goto 999
  !--- set reasonable defaults for plot style
  call gxspmt
  !--- loop over curves
  do ic=1,nset
     !--- get curve parameter set ref.
     kset=max(1,min(maxset,icvref(ic)))
     !--- get curve plot parameters
     line=icvpar(4,kset)
     ibar=icvpar(7,kset)
     isplin=icvpar(5,kset)
     mark=min(5,icvpar(8,kset))
     if(mark.ne.0) then
        isym=0
     else
        isym=icvpar(9,kset)
     endif
     !--- set window
     do j=1,4
        ac(j)=actwin(j,kset)
        wn(j)=window(j,kset)
     enddo
     call jswn(inormt,wn(1),wn(2),wn(3),wn(4))
     !     call jselnt(inormt)
     !
     !--- plot curves
     !
     !   color index and line thickness
     if(line.ne.0.or.ibar.ne.0) then
        call jsplci(icvpar(6,kset))
        call jslwsc(float(icvpar(3,kset)))
     endif
     if(line.ne.0) then
        !   polyline style
        call jsln(line)
        if(isplin.eq.0)  then
           !   plot polyline
           call gxpl(nptval(ic),xval(ipxval(ic)),yval(ipyval(ic)),ac)
        else
           !   smooth with a third order spline
           call gxplt1(nptval(ic),xval(ipxval(ic)),yval(ipyval(ic)),   &
                ac)
        endif
     endif
     if(mark.ne.0) then
        !   set marker type
        call jsmk(mark)
        !   plot marker at point positions
        call gxpm(nptval(ic),xval(ipxval(ic)),yval(ipyval(ic)),ac)
     endif
     if(ibar.ne.0) then
        !   vertical bars to lower x axis position (whether x axis plotted or not)
        yy(1)=axwndy(1,kset)
        call jsln(1)
        do j=0,nptval(ic)-1
           xx(1)=xval(ipxval(ic)+j)
           xx(2)=xx(1)
           yy(2)=yval(ipyval(ic)+j)
           call gxpl(2,xx,yy,ac)
        enddo
     endif
     if(isym.ne.0) then
        !--- center character on point
        call jstxal(2,3)
        !--- set character height
        call gxschf(1,1,0.001*icvpar(10,kset),dum1,dum2)
        !--- get plot character
        sss=splotc(kset:kset)
        !--- set ndc because of character sizes  (curves with different scales)
        call jswn(inormt,0.,1.,0.,1.)
        !        call jselnt(inormt)
        fsx=1./(wn(2)-wn(1))
        fsy=1./(wn(4)-wn(3))
        do j=0,nptval(ic)-1
           xs=fsx*(xval(ipxval(ic)+j)-wn(1))
           ys=fsy*(yval(ipyval(ic)+j)-wn(3))
           call gxtx1(xs,ys,sss,ac)
        enddo
     endif
  enddo
  !--- restore previous settings
  call gxrest(isave,rsave)
999 end subroutine gxcrv1
subroutine gxcurv(nset,nptval,ipxval,ipyval,icvref,xval,yval,window,ierr)
  use gxx11_common
  implicit none
  integer ibar,ic,ierr,isplin,isym,j,kset,line,nset,mark
  real dum1,dum2,fsx,fsy,xs,ys
  !***********************************************************************
  !
  !   Purpose: plots curves into an existing frame
  !
  !--- Input
  !   nset       number of curves (=ordered sets of (x,y) pairs) to plot
  !   nptval(i)  number of points ((x,y) pairs) in set I
  !   ipxval(i)  first x value of set I in array XVAL
  !   ipyval(i)  first y value of set I in array YVAL
  !   icvref(i)  number of the parameter set to be used for curve I. This value
  !            will be forced into [1,MAXSET].
  !            The x and y axis reference numbers of set I will be taken
  !            from this parameter set. All x and y axes with the
  !            corresponding reference numbers will be (scaled if automatic)
  !            and plotted together with set I.
  !            If no x resp. y axis exists with the reference number
  !            from the parameter set, the curve will be plotted with
  !            automatic scaling, but without x resp. y axis.
  !   xval       array containing the x values for all sets
  !   yval       array containing the y values for all sets
  !   window(j,I) GKS window (J=1...4) to be used with curve I. These values
  !            can be obtained from routine GXFRAM
  !--- Output
  !   ierr       0 if everything OK, else GKS error, or
  !            1 : GXINIT not called (initialization)
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: March 7, 1988
  !
  !***********************************************************************

  integer nptval(*),ipxval(*),ipyval(*),icvref(*)
  real xval(*),yval(*),window(4,*)
  real wn(4),xx(2),yy(2),rsave(20)
  integer isave(20)
  character sss*1
  !
  !--- get current user settings and keep them
  call gxsave(isave,rsave,ierr)
  if(ierr.ne.0)  goto 999
  !--- set reasonable defaults for plot style
  call gxspmt
  !--- loop over curves
  do ic=1,nset
     !--- get curve parameter set ref.
     kset=max(1,min(maxset,icvref(ic)))
     !--- get curve plot parameters
     line=icvpar(4,kset)
     ibar=icvpar(7,kset)
     isplin=icvpar(5,kset)
     mark=min(5,icvpar(8,kset))
     if(mark.ne.0) then
        isym=0
     else
        isym=icvpar(9,kset)
     endif
     !--- set window
     do j=1,4
        wn(j)=window(j,kset)
     enddo
     call jswn(inormt,wn(1),wn(2),wn(3),wn(4))
     !     call jselnt(inormt)
     !
     !--- plot curves
     !
     !   color index and line width
     if(line.ne.0.or.ibar.ne.0) then
        call jsplci(icvpar(6,kset))
        call jslwsc(float(icvpar(3,kset)))
     endif
     if(line.ne.0) then
        !   polyline style
        call jsln(line)
        if(isplin.eq.0)  then
           !   plot polyline
           call gvpl(nptval(ic),xval(ipxval(ic)),yval(ipyval(ic)))
        else
           !   smooth with a third order spline
           call gxplts(nptval(ic),xval(ipxval(ic)),yval(ipyval(ic)))
        endif
     endif
     if(mark.ne.0) then
        !   set marker type
        call jsmk(mark)
        !   plot marker at point positions
        call gvpm(nptval(ic),xval(ipxval(ic)),yval(ipyval(ic)))
     endif
     if(ibar.ne.0) then
        !   vertical bars to lower x axis position (whether x axis plotted or not)
        yy(1)=axwndy(1,kset)
        call jsln(1)
        do j=0,nptval(ic)-1
           xx(1)=xval(ipxval(ic)+j)
           xx(2)=xx(1)
           yy(2)=yval(ipyval(ic)+j)
           call gvpl(2,xx,yy)
        enddo
     endif
     if(isym.ne.0) then
        !--- center character on point
        call jstxal(2,3)
        !--- set character height
        call gxschf(1,1,0.001*icvpar(10,kset),dum1,dum2)
        !--- get plot character
        sss=splotc(kset:kset)
        !--- set ndc because of character sizes  (curves with different scales)
        call jswn(inormt,0.,1.,0.,1.)
        !        call jselnt(inormt)
        fsx=1./(wn(2)-wn(1))
        fsy=1./(wn(4)-wn(3))
        do j=0,nptval(ic)-1
           xs=fsx*(xval(ipxval(ic)+j)-wn(1))
           ys=fsy*(yval(ipyval(ic)+j)-wn(3))
           call gxstx(xs,ys,sss)
        enddo
     endif
  enddo
  !--- restore previous settings
  call gxrest(isave,rsave)
999 end subroutine gxcurv
subroutine gxdfmt(axlow,axup,intv,ival,iscal,fmt)
  implicit none
  integer i,i1,i2,ii,intv,iscal,ival,j,mform
  real axl,axlow,axu,axup,fact,step,up,x,y
  !***********************************************************************
  !
  !   Purpose: calculates reasonable format for axis labels
  !
  !--- Input
  !   axlow, aXUP  axis range
  !   intv         no. of intervals, or 0 if not known
  !--- Output
  !   ival         length of format (e.g. 8 for F8.2)
  !   iscal        power of 10 extracted from axlow and axup
  !   fmt          format in correct form, e.g. '(F8.2)'
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Feb. 25, 1991
  !
  !***********************************************************************
  character fmt *(*)
  parameter (mform=8)
  character form(mform)*8
  integer iv(mform),ic(mform)
  save form, iv, ic, up
  data form/'(G10.4)','(F6.1)','(F6.2)','(F5.2)','(F6.3)', '(F7.4)',&
       '(F8.5)','(F9.6)'/
  data iv/10,6,6,5,6,7,8,9/
  data ic/4,1,2,2,3,4,5,6/
  data up/999./
  x=max(abs(axlow),abs(axup))
  if (x .eq. 0.)  then
     iscal = 0
     i = 1
     goto 30
  endif
  i1 = log10(x)
  if (i1 .gt. 3 .or. i1 .le. -3)  then
     iscal = 3 * (i1 / 3)
  else
     iscal = 0
  endif
  fact = 10.**(-iscal)
  axl = fact * axlow
  axu = fact * axup
  y=axu-axl
  if (intv .gt. 0)  then
     !--- get all digits of step if possible
     step = y / intv
     do  i1 = 0, 4
        if (step .ge. 0.99                                            &
             .and. step - int(step + 0.5) .lt. 0.01)  goto 2
        step = 10. * step
     enddo
2    step = y / intv + abs(axl) - int(abs(axl))
     do  i2 = 0, 4
        if (step - int(step + 0.0001) .lt. 0.01)  goto 4
        step = 10. * step
     enddo
4    i = max(i1, i2)
     ii = abs(log10(x)) + 1.001
     if (axl .lt. 0.)                                                &
          ii = max( max(log10(x), log10(-axl) + 1.) + 1.001, 2.001)
     if (i + ii .ge. 9)  then
        i = 1
        goto 30
     else
        ival = i + ii + 1
        fmt = ' '
        write(fmt(:6), '(''(F'',I1,''.'',I1,'')'')')  ival, i
     endif
     goto 999
  else
     do i=1,mform
        if(x.ge.up) goto 20
        x=10.*x
     enddo
     i=1
     goto 30
20   continue
     ii=i
     do  j=ii,mform
        if(y.ge.10.**(1-ic(i)))  goto 30
        i=i+1
     enddo
     i=1
  endif
30 continue
  fmt=form(i)
  ival=iv(i)
999 end subroutine gxdfmt
subroutine gxdfvm(sin,sout,nml)
  implicit none
  integer i1,i2,i3,i4,jb1,jb2,nml
  !***********************************************************************
  !
  !   Purpose: returns the VM filename (fn ft fm)
  !
  !--- Input
  !   sin        ruser Input - either fn, or fn ft, or fn ft fm
  !--- Output
  !   sout       complete fn ft fm filename
  !   nml        last character of file name in SIN
  !
  !   Author: H. Grote / CERN                        date: April 7, 1988
  !                                           last mod: April 7, 1988
  !
  !***********************************************************************
  character(*)  sin,sout
  character(20) sloc
  call gxpnbl(sin,i1,i2)
  if(i1.eq.0)  then
     !--- user Input is totally blank
     sloc='GXMETA   METAFILE A'
     nml=1
     goto 500
  endif
  jb1=index(sin(i1:i2),' ')
  if(jb1.eq.0)  then
     !--- user Input is one piece
     sloc=sin(i1:i2)
     sloc(9:)=' METAFILE A'
     nml=i2
     goto 500
  endif
  jb1=jb1+i1-1
  call gxpnbl(sin(jb1:),i3,i4)
  i3=i3+jb1-1
  jb2=index(sin(i3:i2),' ')
  nml=jb1-1
  if(jb2.eq.0)  then
     !--- user Input two pieces
     sloc=sin(i1:i2)//' A'
  else
     !--- user Input three pieces
     sloc=sin(i1:i2)
  endif
500 sout=sloc
end subroutine gxdfvm
subroutine gxdint(axlow,axup,intv)
  implicit none
  integer i,intv,mrange
  real axlow,axup,d,dn,x
  !***********************************************************************
  !
  !   Purpose: calculates reasonable number of axis intervals
  !
  !--- Input
  !   axlow, aXUP  axis range
  !--- Output
  !   intv         number of intervals
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************
  parameter (mrange=10)
  integer iv(mrange)
  real rangl(mrange)
  save iv, rangl
  data iv/10,6,8,10,10,6,8,10,6,8/
  data rangl/1.,1.2,1.6,2.,2.5,3.,4.,5.,6.,8./
  d=axup-axlow
  if(d.gt.0.)  then
     x=100.+log10(d)
     dn=10.**(x-int(x + 0.0001))
     do i=1,mrange
        if(abs(dn-rangl(i))/rangl(i).lt.1.e-3)  then
           intv=iv(i)
           goto 999
        endif
     enddo
  endif
  intv=10
999 end subroutine gxdint
subroutine gxeopn(string,number)
  use gxx11_common
  implicit none
  integer number
  !***********************************************************************
  !
  !   Purpose: transfers unit number to common block for files opened
  !         externally
  !
  !--- Input
  !   string  (character) option :
  !         'mETA' for metafile, 'ERROR' for error file
  !   number  unit number
  !
  !   Author: H. Grote / CERN                        date: Dec. 21, 1987
  !
  !***********************************************************************

  character string *(*),sloc *4
  call gxundf
  sloc=string
  if(sloc.eq.'META')  then
     lmetop=lundef
     call gxsvar('IMETUN',number,0.,' ')
  elseif(sloc.eq.'ERRO')  then
     lerrop=lundef
     call gxsvar('IERRUN',number,0.,' ')
  endif
end subroutine gxeopn
subroutine gxfchr(imode, ch, ifont, width, np, ipen, x, y, ierr)
  implicit none
  integer i,ierr,ifont,imode,ip,ipos,isel,istr,iwid,j,k,kbit,kword, &
       lx,ly,np
  real width
  !***********************************************************************
  !
  !   Purpose: returns the polygon for a character
  !
  !--- Input:
  !   imode    =0: give character widths only, else give all
  !   ch       character
  !   ifont    font (only 1 or -13)
  !--- Output:
  !   width    character width
  !   np       # points in polygon
  !   ipen     0 for pen up, 1 for pen down
  !   x        x coordinates
  !   y        y coordinates
  !   ierr     =0: OK, =1: wrong font, =2: character not found
  !
  !***********************************************************************
  character(1) ch
  real x(*), y(*)
  integer ipen(*)
  integer nchinf(2), ichinf(95,2), ichcod(652,2)
  character(100) chstr(2)
  save chstr, nchinf, ichinf, ichcod
  data nchinf / 95, 91 /
  data chstr /                                                      &
       ' !"#$%&''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[/]^_`abcdefghijklmnopqrstuvwxyz{|}~',&
       ' !"#$%&''()*+,-./0123456789:;<=>?@ABCDEFGHIKLMNOPRSTUVWXYZ[/]^_`abcdefghiklmnoprstuvwxyz{|}~' /
  data (ichinf(j,1), j = 1,  95) /  541065217, 556276737, 574633992,&
       592715798, 609247262, 627078198, 644384851, 656940149, 674769020, &
       691546246, 708843664, 728240278, 740827290, 761792674, 774378660, &
       794298537, 810566827, 827330748, 844118208, 860896462, 877663453, &
       894452962, 911236339, 927994122, 944796942, 961568043, 975710530, &
       992490828,1012927833,1030230364,1046482272,1061703011,1080862070, &
       1095244191,1112823205,1129597370,1146370508,1162617306,1179130338,&
       1196709352,1213733373,1226836483,1245719045,1263802895,1279529493,&
       1298145817,1314396705,1331189287,1347695164,1364745800,1381251679,&
       1397772909,1413485185,1431841413,1447563919,1465918099,1481642651,&
       1497896607,1515199140,1530403498,1549273776,1563957938,1582829240,&
       1600654011,1613242045,1632387780,1649165012,1665678052,1682719474,&
       1699235586,1714429715,1733056282,1749821231,1763712824,1781017407,&
       1799625545,1814039375,1836598097,1850484577,1867269994,1884046203,&
       1900823435,1916018587,1933854626,1949310899,1967925178,1983910851,&
       2002265031,2017727439,2034246611,2051283931,2067290081,2082474998,&
       2100844536,2120232973/
  data (ichinf(j,2), j = 1,  95) /  541065217, 556276737, 574633992,&
       592715798, 609247262, 627078198, 644122707, 656940149, 674769020, &
       691546246, 708843664, 727978134, 740827290, 761530530, 774378660, &
       794298537, 810566827, 827330748, 844118208, 860896462, 877663453, &
       894452962, 911236339, 927994122, 944796942, 961568043, 975710530, &
       992490828,1012927833,1029968220,1046482272,1061703011,1080862070, &
       1095244191,1112823205,1129320890,1145575870,1162617284,1179667916,&
       1195643359,1213733347,1226836457,1263802859,1279791601,1298145781,&
       1314396669,1330659843,1347951123,1381249561,1397235237,1413485100,&
       1431852592,1447577157,1465401942,1481120358,1498963564,1515199107,&
       1530403465,1549273743,1563957905,1582829207,1600391834,1613242012,&
       1632918179,1649177273,1665675989,1682464481,1698711288,1717064457,&
       1733050141,1750091564,1764500285,1799897925,1816141653,1834240860,&
       1850225518,1868322682,1884826511,1917340569,1934641067,1951407036,&
       1967409090,1984970707,2002537442,2017487860,2035574798,2050774056,&
       2067289149,2082474066,2100843604,2120233065, 4 * 0/
  data (ichcod(j,1), j = 1, 117) /   43336327,  42091009,1115702017,&
       1115816725,1117012498,1108361871,1125139089, 102057364,1158825232,&
       1167017488,1167132057,1111951513,1162281484,1225523590,1216742425,&
       1145308697,1178863762,1200899605,1142243988,1100104080,1108230797,&
       1133266570,1200179208,1216759939,1199654400,1140867713,1099106965,&
       1098908693,1158890769,1150239630,1116619152,1100104212,1125467157,&
       1158956691,1209223572,1251281031,1199982340,1191331840,1224755713,&
       1250118277,1233602695, 193743757,1259227790,1242384779,1216759683,&
       1182877056,1132479105,1107444100,1099317768,1116292621,1183729424,&
       1192380052,1167410324,1142047760,1150109066,1208174849,1241533184,&
       1266764674,  51724948,1108492816,1116685072,1116800409,1150763924,&
       1116750347,1107772034,1136805061,1170669977,1117209492,1150305547,&
       1158104194,1136804549,1103561749,1141440914,1183581842,1099695762,&
       1182794249,1258881793,1115701761,1115833089,1128350403,1111753225,&
       1258881666,1107378816,1124156034, 169427271,  76890900,1108427148,&
       1099514372,1124156544,1166034689,1208240265,1217153041,1192510869,&
       1150616337,1142048149,1166017040,1108427411,1125401621,1175799572,&
       1200834577,1208960909,1183465856,1216348821,1209353485,1183664012,&
       1208699016,1216759811,1191265664,1140867713,1107444100, 110444935,&
       1225197205,1182795669,1117078028,1116554254,1166952205,1208699016/
  data (ichcod(j,1),j= 118, 234) / 1216759811,1191265664,1140867713,&
       1107444100, 135415700,1175799061,1133789841,1108099591,1115898753,&
       1157645696,1191266307,1216759943,1208633100,1166886157,1133265546,&
       1107757205,1132462485,1217725461,1117012498,1108361870,1133331852,&
       1191921673,1216825476,1208108929,1174422528,1115767298,1099186567,&
       1107903243,1150043789,1200506896,1209157524,1175798805, 135153547,&
       1183401224,1149780745,1108033934,1099907602,1125401749,1159022228,&
       1200769038,1208567684,1182876928,1140867713,1107493518,1108165260,&
       1124942478,  42091009,1115702017,1115816590,1108165260,1124942478,&
       50414208,1107378818,1124156225,1120092740, 168968713,1241514508,  &
       1259078150,1258684946,1242120704,  26231185,1108558484,1133856149,&
       1184122643,1200703375,1192052364,1149912199,  75645953,1149256961,&
       1149372685,1217349520,1175471375,1150174219,1141392518,1166362373,&
       1216760072, 152062214,1233472133,1267158026,1275874191,1259424275,&
       1226065813,1175798932,1133724305,1108296076,1099514374,1115964290,&
       1149322752,1199589633,1241664131,  76890240,  76892288,  34031367,&
       34947584,  34948757,1209288851,1225869583,1217218572,1183515147,  &
       1183533066,1216956679,1225017474,1208043136,1107298576,1217546132,&
       1184187541,1133789842,1108361613,1099448837,1115898753,1149257344,&
       1199655043,1225064981,1107296789,1167410964,1209157776,1225607432/
  data (ichcod(j,1),j= 235, 351) / 1216694275,1191265664,1107296789,&
       1107296789,1217724939,1175126528,1216348693,1107296789,1217724939,&
       1175128336,1217546132,1184187541,1133789842,1108361613,1099448837,&
       1115898753,1149257344,1199655043,1225083144, 109594888,  34947584,&
       152389888,  34294027,  34947584, 102057477,1166165249,1140867840, &
       1107378562,1090863367,  34947584, 152388103,  76302592,  34947584,&
       33572864,  34947584,  34948608, 169166336, 169167360,  34947584,  &
       34949376, 152389888,  76891028,1116881424,1099776392,1107640963,  &
       1132545152,1182812033,1216563461,1233668493,1225803922,1200899733,&
       1150616085,1107296789,1184188436,1217612049,1225672844,1208698506,&
       1107952789,1133789842,1108361613,1099448837,1115898753,1149257344,&
       1199655043,1225083272,1233996048,1217546132,1184187541, 100944194,&
       34947584,  34948757,1209288851,1225869583,1217218572,1183531531,  &
       93014272, 143804308,1175798805,1117012370,1099973134,1116554124,  &
       1183467401,1208502406,1216563073,1174422528,1115767171,  68502528,&
       9783189,  34947590,1115898753,1157645824,1199655043,1225148693,   &
       9782400, 144000128,  18170752, 102056832, 102058112, 185944192,   &
       26560640, 143999360,   9782411,1149241493,1149962389,1098908053,  &
       1217724800,1216348697,1111949849,1167655495,1170669849,1246168345,&
       1162281369,1159266759,1162281351,1166821767,  37964611,  34554512/
  data (ichcod(j,1),j= 352, 468) / 1125270292,1117078036,1116931982,&
       1199572875,1183663502,1141785357,1108033928,1099317763,1124156416,&
       1166034561,1199768085,1107296779,1124942862,1166952077,1200310280,&
       1208371075,1182877056,1140867841,1107494795,1183663502,1141785357,&
       1108033928,1099317763,1124156416,1166034561,1199769493,1199572875,&
       1183663502,1141785357,1108033928,1099317763,1124156416,1166034561,&
       1199767944,1200113546,1191986829,1166951438,1124942347,1099448710,&
       1107510017,1140868480,1182877571,  85279765,1125401233,1115685134,&
       1150158734,1203914565,1187399111,1145520966, 126568077,1166951438,&
       1124942347,1099448710,1107510017,1140868480,1182877571,  34947584,&
       34227085,1150174734,1192052618,1199571349,1108624021,1108754837,  &
       34488832,  43336468,1133855510,1117061902,1128481478,1103577287,  &
       34947584, 118374916,  67651456,  34947584,  34488832,  34227085,  &
       1150174734,1192052618,1199572874,1225607694,1267616909,1292520704,&
       34488832,  34227085,1150174734,1192052618,1199571982,1124942347,  &
       1099448710,1107510017,1140868480,1182877571,1208371208,1200309901,&
       1166951438,  34488903,  34292493,1141785998,1183664011,1208502278,&
       1199785601,1166033920,1124155907, 126764999, 126568077,1166951438,&
       1124942347,1099448710,1107510017,1140868480,1182877571,  34488832,&
       34095755,1133331598,1175324427,1183663374,1133396493,1099645449/
  data (ichcod(j,1),j= 469, 585) / 1124615559,1183205124,1191396993,&
       1157645184,1107378563,  43336324,1124156416,1157628174,1150157326,&
       1107575425,1132479744,1174488964, 126764928,  17712128, 118375424,&
       26100608,  93209472,  93210496, 160319360,  26101504, 118374784,  &
       17712128, 118375424,1128546886,1095188679, 118374784,  26101518,  &
       25184000,  77153176,1125597845,1116947217,1133528078,1141654282,  &
       1107903240,1141261316,1132610305,1119961795,1128612806,1153892889,&
       1111949977,1134052375,1150633107,1141982096,1125008140,1141523721,&
       1141392134,1124352898,1140933825,1153647685,1137066695,  25575816,&
       1108034316,1141654795,1191659526,1225148935,1250577036,  76891028,&
       1116881424,1099776392,1107640963,1132545152,1182812033,1216563461,&
       1233668493,1225803922,1200899733,1150615746,1251410569,1191772427,&
       1267400966,1267073296,1267728644,1266942219,1267400978,1267859723,&
       1108165390,1166952077,1200310144, 126437128,1174881160,1116225926,&
       1090797827,1098990208,1174423297,1199702549,1107378816,1132463633,&
       1175602580,1150632853,1117012498,1108361870,1133266187,1175012999,&
       1183139331,1157694221,1108033929,1099383301,1124287618,1166034498,&
       1178879430,1153909703,1120289348,1111688341,1082131605,1216348679,&
       1191641368,1100431896,1100562712, 119031703,1209550745,1192756373,&
       1133789842,1108361613,1099448837,1115898753,1149257344,1199655043/
  data (ichcod(j,1),j= 586, 652) / 1225083272,1233996048,1217546132,&
       1184187541,  35144343,1125663385,1108871192,1217874200,1218005016,&
       34947590,1115898753,1157645824,1199655043,1225148693,  43533079,  &
       1134052121,1117259672,1209485464,1209616280, 126764928, 126568077,&
       1166951438,1124942347,1099448710,1107510017,1140868480,1182877571,&
       34947732,1125466774,1108674069,1184122645,1184253461,  68043533,  &
       1108033928,1099317763,1124156416,1166034561,1199785990,1208502155,&
       1183663502,1141768725,1117012757,1117143573, 110446356,1200965398,&
       1184170510,1107575425,1132479744,1174488964, 126764928,  43336468,&
       1133855510,1117062677,1184122645,1184253461,  60113556,1108427264,&
       60114069,1167345170,1175471502,1150092173,1166886540,1200244743,  &
       1208371075,1182877056,1149256449,1115881472/
  data (ichcod(j,2),j=   1, 117) /   43336327,  42091009,1115702017,&
       1115816725,1117012498,1108361871,1125139089, 102057364,1158825232,&
       1167017488,1167132057,1111951513,1162281484,1225523590,1216742425,&
       1145308697,1178863762,1200899605,1142243988,1100104080,1108230797,&
       1133266570,1200179208,1216759939,1199654400,1140867713,1099106965,&
       1098908693,1158890769,1150239630,1116619152,1100104212,1125467157,&
       1158956691,1209223572,1251281031,1199982340,1191331840,1224755713,&
       1250118277,1233602695, 193743757,1259227790,1242384779,1216759683,&
       1182877056,1132479105,1107444100,1099317768,1116292621,1183729424,&
       1192380052,1167410324,1142047760,1150109066,1208174849,1241533184,&
       1266764674,  51724948,1108492816,1116685072,1116800409,1150763924,&
       1116750347,1107772034,1136805061,1170669977,1117209492,1150305547,&
       1158104194,1136804549,1103561749,1141440914,1183581842,1099695762,&
       1182794249,1258881793,1115701761,1115833089,1128350403,1111753225,&
       1258881666,1107378816,1124156034, 169427271,  76890900,1108427148,&
       1099514372,1124156544,1166034689,1208240265,1217153041,1192510869,&
       1150616337,1142048149,1166017040,1108427411,1125401621,1175799572,&
       1200834577,1208960909,1183465856,1216348821,1209353485,1183664012,&
       1208699016,1216759811,1191265664,1140867713,1107444100, 110444935,&
       1225197205,1182795669,1117078028,1116554254,1166952205,1208699016/
  data (ichcod(j,2),j= 118, 234) / 1216759811,1191265664,1140867713,&
       1107444100, 135415700,1175799061,1133789841,1108099591,1115898753,&
       1157645696,1191266307,1216759943,1208633100,1166886157,1133265546,&
       1107757205,1132462485,1217725461,1117012498,1108361870,1133331852,&
       1191921673,1216825476,1208108929,1174422528,1115767298,1099186567,&
       1107903243,1150043789,1200506896,1209157524,1175798805, 135153547,&
       1183401224,1149780745,1108033934,1099907602,1125401749,1159022228,&
       1200769038,1208567684,1182876928,1140867713,1107493518,1108165260,&
       1124942478,  42091009,1115702017,1115816590,1108165260,1124942478,&
       50414208,1107378818,1124156225,1120092740, 168968713,1241514508,  &
       1259078150,1258684946,1242120704,  26231185,1108558484,1133856149,&
       1184122643,1200703375,1192052364,1149912199,  75645953,1149256961,&
       1149372685,1217349520,1175471375,1150174219,1141392518,1166362373,&
       1216760072, 152062214,1233472133,1267158026,1275874191,1259424275,&
       1226065813,1175798932,1133724305,1108296076,1099514374,1115964290,&
       1149322752,1199589633,1241664131,  76890240,  76892288,  34031367,&
       34947584,  34948757,1209288851,1225869583,1217218572,1183515147,  &
       1183533066,1216956679,1225017474,1208043136,1107296661,1216348544,&
       1217725589,1082131605,1216348288,1216348693,1107296789,1217724939,&
       1175126528,1216349461,1157628944,1116684814,1099710857,1107772038/
  data (ichcod(j,2),j= 235, 351) / 1141196293,1199982599,1216956556,&
       1208895375,1175471120,  34947584,  34949141,  34947584, 152389888,&
       34294027,  34947584,  34947584, 152388103,  76302592,  76890240,  &
       76892288,  34947584,  34948608, 169166336, 169167360,  34947584,  &
       34949376, 152389888,  25183104,1107771787,1099907602,1125401749,  &
       1167410964,1209157775,1217087495,1182812288,  34947584, 152389888,&
       34949397,  34947584,  34948757,1209288851,1225869582,1217153035,  &
       1183465994,  18171019,1090519317,1209336064,1207960597,1140850837,&
       1200948373,1133789842,1108361613,1099448837,1115898753,1149257344,&
       1199655043,1225083272,1233996048,1217546132,1184187541,  17842450,&
       1100235285,1125467028,1142047886,1149241360,1209157524,1192576533,&
       1167344914,1150158229,1166016783,1099907598,1116357384,1132938502,&
       1174816647,1208502410,1225673103,1242497301,1209336587,1175126272,&
       1207960725,1133789842,1108361613,1099448837,1115898753,1149257344,&
       1199655043,1225083272,1233996048,1217546132,1184187541,  67847947,&
       143999360,  26560661,  25184384,  35209799,  35210649,  38225351, &
       18434631,  85542215,  26821913,  29836615,  25642380,1233584707,  &
       1262682639,1116750610,1125401237,1108624019,  76432269,1116422665,&
       1099317635,1107378944,1140868353,1183074183,1217087758,  76432782,&
       1175275147,1199785985,1216366848, 102057236,1142047502,1116422663/
  data (ichcod(j,2),j= 352, 468) / 1098989895, 102057749,1209223184,&
       1200506637,1175209100,  76301707,1183401735,1191462530,1174488320,&
       1140867841,1115832837,  17711630,1124877893,1195853895, 143542284,&
       1191789122,1095057607,  93209614,1124942347,1099448709,1107444353,&
       1132479616,1166100099,1191593737,1183597966,1150305298,1142178965,&
       1167410836,1200752268,1175274766,1133396621,1116422921,1149764744,&
       1116160389,1099121153,1124091008,1166100099,  67977996,1107968391,&
       1099186690,1115767680,1157645953,1208174854,1233734028,1217283982,&
       1183597960,1149453127,   9126285,1116619534,1141720204,1158235397,&
       1149241486,1208698761,1149256644,1128726666,1091322382,1125008269,&
       1133200135,1107297031,1141589261,1175340814,1208764425,1199851079,&
       51266055,1099121025,1107313408,1140999300,  51265792, 135088014,  &
       1192117773,1141457672,1116209800,1132938246,1157711232,1174423169,&
       9781653,1117012755,1191183374,1090519950,1086784266,1116029570,   &
       1132479616,1166100099,1200031886,1200047875,1191266176,1216366978,&
       1241776526,1125008008,1107509632, 135153547,1191790086,1149453057,&
       1098908686,1124942346,1099383172,1107378816,1132479617,1166280200,&
       1166296577,1182812032,1216432516,1241991690,1233996046,  76432000,&
       118376328,1208174720,  17515021,1133398542,  34095621,1115833089, &
       1140868352,1174488835,1199982473,1191986829,1166951566,1133331083/
  data (ichcod(j,2),j= 469, 570) / 1107836999, 151929870,1124942347,&
       1099448709,1107444353,1132479616,1166100099,1191593737,1183598093,&
       1158546830,1140850955,1108165518,1225655310,1124942347,1099448709,&
       1107444353,1132479616,1166100099,1191593737,1183598093,1158562830,&
       9060620,1108230926,1133331339,1116029570,1132479616,1174488835,   &
       1208436875,1217267733,1145503882,1091322382,1125008269,1133200134,&
       1124287361,1149257088,1191266307,1225148939,1250821397,1142178707,&
       1133659153,1167083280,  93340687,1125008012,1116357512,1158104583,&
       84361990,1107640707,1098990273,1153647940,1162232903,1128726666,  &
       1091322382,1125008269,1133200134,1124287361,1140868352,1174488836,&
       1199982601,1217284241,1209288469,1175799187,1167148558,1191921673,&
       1233585429,1142178707,1133659153,1167083280, 118506766,1133265417,&
       1099317636,1107444480,1153582404,1162233031,1137132357,  77153176,&
       1125597845,1116947217,1133528078,1141654282,1107903240,1141261316,&
       1132610305,1119961795,1128612806,1153892889,1111949977,1134052375,&
       1150633107,1141982096,1125008140,1141523721,1141392134,1124352898,&
       1140933825,1153647685,1137066695,  25575816,1108034316,1141654795,&
       1191659526,1225148935,1250577036/
  data (ichcod(j,2),j= 571,652) / 82 * 0 /

  chstr(1)(61:61) = '\\'
  chstr(2)(59:59) = '\\'
  if (ifont .eq. 1)  then
     isel = 1
  elseif (ifont .eq. -13)  then
     isel = 2
  else
     ierr = 1
     goto 999
  endif
  ipos = index(chstr(isel)(:nchinf(isel)), ch)
  if (ipos .eq. 0)  then
     ierr = 2
  else
     ierr = 0
     iwid = 0
     np = 0
     istr = 0
     call cbyt(ichinf(ipos,isel), 19, iwid, 1, 6)
     width = 0.01 * iwid
     if (imode .ne. 0)  then
        call cbyt(ichinf(ipos,isel), 11, np, 1, 8)
        call cbyt(ichinf(ipos,isel),  1, istr, 1, 10)
        do i = 1, np
           k     = istr + i
           kword = k / 2
           kbit  = 17 - 16 * (k - 2 * kword)
           ipen(i) = 0
           call cbyt(ichcod(kword,isel), kbit, ipen(i), 1, 16)
           if (ipen(i) .ge. 16384)  then
              ip = 1
              ipen(i) = ipen(i) - 16384
           else
              ip = 0
           endif
           lx   = ipen(i) / 128
           x(i) = 0.01 * lx
           ly   = ipen(i) - 128 * lx
           if (ly .ge. 64)  then
              y(i) = 0.01 * (64 - ly)
           else
              y(i) = 0.01 * ly
           endif
           ipen(i) = ip
        enddo
     endif
  endif
999 end subroutine gxfchr
character(60) function gxform(string)
  implicit none
  integer i,ipt,kmant,l,n
  !
  !   creates the correct format for Input variables  contained  in  a
  !   character  variable  in  free format (blank characters acting as
  !   separators !). It accepts I,F,E,D,L, and A format variables, the
  !   latter  without  any  quotes,   being   just   those   character
  !   combinations  which  cannot  be  attributed  to any of the other
  !   formats.
  !
  !   gxform overcomes the short-coming of FORTRAN77 to not allow  free
  !   format  READ  statements  from  internal  files,  i.e. character
  !   variables.  If STRING is a character variable, then
  !
  !   read(strING,*)  N,A,etc.
  !
  !   is  not  legal  ANSI  FORTRAN  (although   supported   by   some
  !   compilers).  In this case one can use
  !
  !   read(strING,GXFORM(STRING))  N,A,etc.
  !
  !   which is legal. GXFORM has to be  declared  CHARACTER*60  in  the
  !   calling routine.
  !
  !   Input
  !   string    character type Input line
  !   Output
  !   gxform    FORMAT, e.g. (I4,A12,L4,F8.3,I3,E12.4,D24.8)
  !
  !   restrictions: the maximum length of the complete  format  is  60
  !   characters.
  !
  !   Author    HG      4.3.86   last mod.: 9.6.86
  !
  character string*(*),stemp*1,sfchar*1,form*80
  logical count,realfl,expfl,number
  form='(A1)'
  ipt=1
  n=1
  kmant=0
  realfl=.false.
  expfl=.false.
  count=.false.
  number=.false.
  sfchar='I'
  do i=1,len(string)
     if(ipt.ge.60) goto 20
     stemp=string(i:i)
     if(stemp.eq.' ') then
        if(count) then
           l=i-n
           ipt=ipt+1
           form(ipt:ipt)=sfchar
           if(l.ge.10) then
              write(form(ipt+1:),'(I2)') l
              ipt=ipt+2
           else
              write(form(ipt+1:),'(I1)') l
              ipt=ipt+1
           endif
           if(realfl.or.expfl)  then
              write(form(ipt+1:),'(I2)')  min(9,kmant)
              form(ipt+1:ipt+1)='.'
              ipt=ipt+2
           endif
           ipt=ipt+1
           form(ipt:ipt)=','
           n=i
           expfl=.false.
           realfl=.false.
           count=.false.
           number=.false.
           sfchar='I'
           kmant=0
        endif
     else
        if(.not.count) then
           !--- first character of a new variable
           if(stemp.eq.'.') then
              !--- could be floating, or logical
              if(index(string(i:),'.T.').eq.1.or.                       &
                   index(string(i:),'.F.').eq.1.or.                                  &
                   index(string(i:),'.TRUE.').eq.1.or.                               &
                   index(string(i:),'.FALSE.').eq.1) then

                 sfchar='L'
              elseif(index('0123456789',string(i+1:i+1)).ne.0) then
                 number=.true.
              else
                 sfchar='A'
              endif
           elseif(index('+-0123456789',stemp).ne.0) then
              number=.true.
           else
              sfchar='A'
           endif
        endif
        count=.true.
        if(number) then
           if(stemp.eq.'E'.or.stemp.eq.'D') then
              expfl=.true.
              sfchar=stemp
           endif
           if(realfl.and..not.expfl) then
              kmant=kmant+1
              sfchar='F'
           endif
           realfl=realfl.or.stemp.eq.'.'
           if(realfl.and..not.expfl)  sfchar='F'
        endif
     endif
  enddo
  if(ipt.ge.4)  form(ipt:ipt)=')'
20 gxform=form
end function gxform
subroutine gxfram(ncurv,nptval,ipxval,ipyval,icvref,xval,yval,window,ierr)
  use gxx11_common
  implicit none
  integer i,iaxr,iayr,ierr,j,jc,kset,ncurv
  real axpos,d,fx,fy
  !***********************************************************************
  !
  !   Purpose: plots one frame with several axes, returns GKS windows
  !
  !--- Input
  !   ncurv      number of curves (=ordered sets of (x,y) pairs) to plot
  !   nptval(i)  number of points ((x,y) pairs) in set I
  !   ipxval(i)  first x value of set I in array XVAL
  !   ipyval(i)  first y value of set I in array YVAL
  !   icvref(i)  number of the parameter set to be used for curve I. This value
  !            will be forced into [1,MAXSET].
  !            The x and y axis reference numbers of set I will be taken
  !            from this parameter set. All x and y axes with the
  !            corresponding reference numbers will be (scaled if automatic)
  !            and plotted together with set I.
  !            If no x resp. y axis exists with the reference number
  !            from the parameter set, the curve will be plotted with
  !            automatic scaling, but without x resp. y axis.
  !   xval       array containing the x values for all sets
  !   yval       array containing the y values for all sets
  !--- Output
  !   window(j,I) GKS window (J=1...4) to be used with curve I. These values
  !            are used by routine GXCURV
  !   ierr       0 if everything OK, else GKS error, or
  !            1 : GXINIT not called (initialization)
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************

  integer nptval(*),ipxval(*),ipyval(*),icvref(*)
  real xval(*),yval(*),window(4,*)
  integer ixax(mxaxs),iyax(myaxs),ixaref(maxset),iyaref(maxset),    &
       ilpar(30)
  real wn(4)
  !
  do i=1,30
     ilpar(i)=0
  enddo
  ierr=0
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  !--- set reasonable defaults for plot style
  call gxspmt
  !--- set axis flags to "not plotted"
  do i=1,mxaxs
     ixax(i)=0
  enddo
  do i=1,myaxs
     iyax(i)=0
  enddo
  !--- set axis curve references to "not set"
  do i=1,maxset
     ixaref(i)=0
     iyaref(i)=0
  enddo
  !--- get x and y axis reference numbers
  do i=1,ncurv
     kset=max(1,min(maxset,icvref(i)))
     ixaref(kset)=icvpar(1,kset)
     iyaref(kset)=icvpar(2,kset)
  enddo
  !--- get all window values for x
  call gxprwn(1,ncurv,icvref,nptval,ipxval,ipxval,xval,xval,        &
       mxaxs,ixaref,ixapar,rangex,axwndx)
  !--- get all window values for y
  call gxprwn(2,ncurv,icvref,nptval,ipxval,ipyval,xval,yval,        &
       myaxs,iyaref,iyapar,rangey,axwndy)
  !--- get the window in NDC into which the plot has to fit, depending on
  !   axis positions, labels, tick marks, etc.
  call gxmarg(ixaref,iyaref,axwndx,axwndy,actwnd)
  !--- identical x and y ratios if requested
  if(isqflg.gt.0)  then
     d=min(actwnd(2)-actwnd(1),actwnd(4)-actwnd(3))
     actwnd(2)=actwnd(1)+d
     actwnd(4)=actwnd(3)+d
  endif
  !--- set window factors
  fx=1./(actwnd(2)-actwnd(1))
  fy=1./(actwnd(4)-actwnd(3))
  !--- loop over curve sets, plot axes
  do kset=1,maxset
     iaxr=ixaref(kset)
     if(iaxr.eq.0) goto 100
     !--- get window according to margin
     wn(1)=(actwnd(2)*axwndx(1,kset)-actwnd(1)*axwndx(2,kset))*fx
     wn(2)=((1.-actwnd(1))*axwndx(2,kset)- (1.-actwnd(2))*axwndx     &
          (1,kset))*fx
     wn(3)=(actwnd(4)*axwndy(1,kset)-actwnd(3)*axwndy(2,kset))*fy
     wn(4)=((1.-actwnd(3))*axwndy(2,kset)- (1.-actwnd(4))*axwndy     &
          (1,kset))*fy
     call jswn(inormt,wn(1),wn(2),wn(3),wn(4))
     !     call jselnt(inormt)
     !--- keep
     do j=1,4
        window(j,kset)=wn(j)
     enddo
     !--- plot x axes
     do i=1,mxaxs
        if(ixapar(21,i).eq.iaxr) then
           if(ixax(i).eq.0) then
              ixax(i)=1
              !--- x axis no. 1 and 2 at bottom, 3 and 4 at top of frame
              if(i.le.2) then
                 axpos=axwndy(1,kset)
              else
                 axpos=axwndy(2,kset)
              endif
              !--- set parameters, get interval number if scaling automatic
              do jc=1,mpaxs
                 ilpar(jc)=ixapar(jc,i)
              enddo
              ilpar(2)=ixapar(19,i)
              call gxaxis('X',axwndx(1,kset), axwndx(2,kset),axpos,     &
                   1,sxform(i),sxtext(i),sdefnl,ilpar,ierr)
              if(ierr.ne.0) goto 999
           endif
        endif
     enddo
     !--- plot y axes
     iayr=iyaref(kset)
     do i=1,myaxs
        if(iyapar(21,i).eq.iayr) then
           if(iyax(i).eq.0) then
              iyax(i)=1
              !--- y axis 1 at left,annotation at left, 2 at left, ann. at right,
              !   3 at right, ann. at left, 4 at right, ann. at right
              if(i.le.2) then
                 axpos=axwndx(1,kset)
              else
                 axpos=axwndx(2,kset)
              endif
              !--- set parameters, get interval number if scaling automatic
              do jc=1,mpaxs
                 ilpar(jc)=iyapar(jc,i)
              enddo
              ilpar(2)=iyapar(19,i)
              call gxaxis('Y',axwndy(1,kset),axwndy(2,kset), axpos,     &
                   1,syform(i), sytext(i),sdefnl,ilpar,ierr)

              if(ierr.ne.0) goto 999
           endif
        endif
     enddo
100  continue
  enddo
999 end subroutine gxfram
subroutine gxfrm1(ncurv,nptval,ipxval,ipyval,icvref,xval,yval,window,actwin,ierr)
  use gxx11_common
  implicit none
  integer i,iaxr,iayr,ierr,j,jc,kset,ncurv
  !***********************************************************************
  !
  !   Purpose: plots one frame with several axes, returns GKS and active
  !         windows.
  !
  !--- Input
  !   ncurv      number of curves (=ordered sets of (x,y) pairs) to plot
  !   nptval(i)  number of points ((x,y) pairs) in set I
  !   ipxval(i)  first x value of set I in array XVAL
  !   ipyval(i)  first y value of set I in array YVAL
  !   icvref(i)  number of the parameter set to be used for curve I. This value
  !            will be forced into [1,MAXSET].
  !            The x and y axis reference numbers of set I will be taken
  !            from this parameter set. All x and y axes with the
  !            corresponding reference numbers will be (scaled if automatic)
  !            and plotted together with set I.
  !            If no x resp. y axis exists with the reference number
  !            from the parameter set, the curve will be plotted with
  !            automatic scaling, but without x resp. y axis.
  !   xval       array containing the x values for all sets
  !   yval       array containing the y values for all sets
  !--- Output
  !   window(j,I) GKS window (J=1...4) to be used with curve I. These values
  !            are used by routines GXCURV and GXCRV1.
  !   actwin(j,I) active window (J=1...4) used to clip curve I. These values
  !            are used by routine GXCRV1.
  !   ierr       0 if everything OK, else GKS error, or
  !            1 : GXINIT not called (initialization)
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************

  integer nptval(*),ipxval(*),ipyval(*),icvref(*)
  real xval(*),yval(*),window(4,*),actwin(4,*)
  integer ixax(mxaxs),iyax(myaxs),ixaref(maxset),iyaref(maxset),    &
       ilpar(30)
  real wn(4),axpos,d,fx,fy
  !
  do i=1,30
     ilpar(i)=0
  enddo
  ierr=0
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  !--- set reasonable defaults for plot style
  call gxspmt
  !--- set axis flags to "not plotted"
  do i=1,mxaxs
     ixax(i)=0
  enddo
  do i=1,myaxs
     iyax(i)=0
  enddo
  !--- set axis curve references to "not set"
  do i=1,maxset
     ixaref(i)=0
     iyaref(i)=0
  enddo
  !--- get x and y axis reference numbers
  do i=1,ncurv
     kset=max(1,min(maxset,icvref(i)))
     ixaref(kset)=icvpar(1,kset)
     iyaref(kset)=icvpar(2,kset)
  enddo
  !--- get all window values for x
  call gxprwn(1,ncurv,icvref,nptval,ipxval,ipxval,xval,xval,        &
       mxaxs,ixaref,ixapar,rangex,axwndx)
  !--- get all window values for y
  call gxprwn(2,ncurv,icvref,nptval,ipxval,ipyval,xval,yval,        &
       myaxs,iyaref,iyapar,rangey,axwndy)
  !--- get the window in NDC into which the plot has to fit, depending on
  !   axis positions, labels, tick marks, etc.
  call gxmarg(ixaref,iyaref,axwndx,axwndy,actwnd)
  !--- identical x and y ratios if requested
  if(isqflg.gt.0)  then
     d=min(actwnd(2)-actwnd(1),actwnd(4)-actwnd(3))
     actwnd(2)=actwnd(1)+d
     actwnd(4)=actwnd(3)+d
  endif
  !--- set window factors
  fx=1./(actwnd(2)-actwnd(1))
  fy=1./(actwnd(4)-actwnd(3))
  !--- loop over curve sets, plot axes
  do kset=1,maxset
     iaxr=ixaref(kset)
     if(iaxr.eq.0) goto 100
     !--- get window according to margin
     wn(1)=(actwnd(2)*axwndx(1,kset)-actwnd(1)*axwndx(2,kset))*fx
     wn(2)=((1.-actwnd(1))*axwndx(2,kset)- (1.-actwnd(2))*axwndx     &
          (1,kset))*fx
     wn(3)=(actwnd(4)*axwndy(1,kset)-actwnd(3)*axwndy(2,kset))*fy
     wn(4)=((1.-actwnd(3))*axwndy(2,kset)- (1.-actwnd(4))*axwndy     &
          (1,kset))*fy
     call jswn(inormt,wn(1),wn(2),wn(3),wn(4))
     !     call jselnt(inormt)
     !--- keep
     do j=1,4
        window(j,kset)=wn(j)
     enddo
     !--- active window in user coordiantes
     actwin(1,kset)=wn(1)+(wn(2)-wn(1))*actwnd(1)
     actwin(2,kset)=wn(1)+(wn(2)-wn(1))*actwnd(2)
     actwin(3,kset)=wn(3)+(wn(4)-wn(3))*actwnd(3)
     actwin(4,kset)=wn(3)+(wn(4)-wn(3))*actwnd(4)
     !--- plot x axes
     do i=1,mxaxs
        if(ixapar(21,i).eq.iaxr) then
           if(ixax(i).eq.0) then
              ixax(i)=1
              !--- x axis no. 1 and 2 at bottom, 3 and 4 at top of frame
              if(i.le.2) then
                 axpos=axwndy(1,kset)
              else
                 axpos=axwndy(2,kset)
              endif
              !--- set parameters, get interval number if scaling automatic
              do jc=1,mpaxs
                 ilpar(jc)=ixapar(jc,i)
              enddo
              ilpar(2)=ixapar(19,i)
              call gxaxis('X',axwndx(1,kset), axwndx(2,kset),axpos,     &
                   1,sxform(i),sxtext(i),sdefnl,ilpar,ierr)
              if(ierr.ne.0) goto 999
           endif
        endif
     enddo
     !--- plot y axes
     iayr=iyaref(kset)
     do i=1,myaxs
        if(iyapar(21,i).eq.iayr) then
           if(iyax(i).eq.0) then
              iyax(i)=1
              !--- y axis 1 at left,annotation at left, 2 at left, ann. at right,
              !   3 at right, ann. at left, 4 at right, ann. at right
              if(i.le.2) then
                 axpos=axwndx(1,kset)
              else
                 axpos=axwndx(2,kset)
              endif
              !--- set parameters, get interval number if scaling automatic
              do jc=1,mpaxs
                 ilpar(jc)=iyapar(jc,i)
              enddo
              ilpar(2)=iyapar(19,i)
              call gxaxis('Y',axwndy(1,kset),axwndy(2,kset), axpos,     &
                   1,syform(i), sytext(i),sdefnl,ilpar,ierr)

              if(ierr.ne.0) goto 999
           endif
        endif
     enddo
100  continue
  enddo
999 end subroutine gxfrm1
subroutine gxinit
  use gxx11_common
  implicit none
  integer idummy,ierr
  !***********************************************************************
  !
  !   Purpose: initializes GKS PLOT package
  !
  !   the default is to open the plot package for metafile writing only.
  !   the corresponding parameters (unit number, file name, status, paper
  !   width and length) can be set by calls to gxsvar beforehand.
  !
  !   for interactive usage (plot Output on screen) it is mandatory to
  !   call gxaSKU before calling GXINIT. Parameters can be modified after
  !   the gxasKU call (but before the GXINIT call) by calling gxsvar.
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  if (ltotin .ne. lundef)  then
     print '(/'' GXPLOT-X11 '',F5.2,'' initialized''/)',versio
  endif
  !--- reset open flag for .eps files
  iepsop = 0
  iepscf = 0
  iclear = 0
  ipage = 0
  call wopks(ierrun, idummy)
  call gxundf
  !--- set default options
  call gxsdef('OPTINIT',0)
  if(interm.ne.0)  then
     itermt = 1
  else
     itermt = 0
  endif
  if (ltseop .ne. lundef .or. itseop .eq. 0)  then
     if (itermt .ne. 0)  call wopwk(interm, mconid, itermt)
  endif
  wxfact=1.
  wyfact=1.
  wfact=1.
  if(inmeta .ne. 0)  then
     !--- orientation (portrait or landscape)
     if (xmetaf .gt. ymetaf) then
        ipstyp = 115
     else
        ipstyp = 114
     endif
     iorips = ipstyp - 113
     ibbox(1) = mlbb1
     ibbox(2) = mlbb2
     ibbox(3) = mubb1
     ibbox(4) = mubb2
     if (xmetaf .gt. 0.)  then
        if (iorips .eq. 1)  then
           ibbox(3) = mlbb1 + mwid1 * xmetaf / mysize + 0.5
        else
           ibbox(4) = mlbb2 + mwid2 * xmetaf / mxsize + 0.5
        endif
     endif
     if (ymetaf .gt. 0.)  then
        if (iorips .eq. 1)  then
           ibbox(4) = mlbb2 + mwid2 * ymetaf / mxsize + 0.5
        else
           ibbox(3) = mlbb1 + mwid1 * ymetaf / mysize + 0.5
        endif
     endif
     if (interm .eq. 0)  then
        imetun = -abs(imetun)
     endif
     if (lpseps .ne. lundef)  then
        ipseps = 1
        lpseps = lundef
     endif
     if (ipseps .eq. 1)  then
        if(lmetop.ne.lundef)  call gxsfop('PSFILE','UNKNOWN',ierr)
        call gxopps(imetun, ipstyp)
     elseif (ipseps .eq. 2)  then
        ipstyp = 113
        iepsop = -2
     endif
  endif
  !--- activate workstations
  if (ltseop .ne. lundef .or. itseop .eq. 0)  then
     if(interm.gt.0) then
        call wacwk(interm)
        lacttm=lundef
     endif
  endif
  !--- set default window and viewport, aspect source flags,
  !   norm. transf. number
  call gxsdef('DEVICE',0)
  !--- axis default values
  call gxsdef('AXIS',0)
  !--- curve defaults
  call gxsdef('CURVE',0)
  !--- set flag that GXINIT has been called
  ltotin=lundef
10000 format(//' GKS error number =',i4,' returned for terminal',       &
       ' device =',i8,'  STOP')
end subroutine gxinit
subroutine gxival(string,ivalex)
  implicit none
  integer i,ifnd,ivalex,n
  !***********************************************************************
  !
  !   Purpose: extracts a positive integer from a string
  !
  !--- Input
  !   string    arbitrary character string
  !--- Output
  !   ivalex    first integer found in string (terminated by any
  !           non-numeric character, including blank)
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Sept. 8, 1987
  !
  !***********************************************************************
  character(*) string
  character snumer*10
  save snumer
  data snumer/'0123456789'/
  ivalex=0
  ifnd=0
  do i=1,len(string)
     n=index(snumer,string(i:i))
     if(n.eq.0) then
        if(ifnd.ne.0) goto 999
     else
        ifnd=1
        ivalex=10*ivalex+n-1
     endif
  enddo
999 end subroutine gxival
subroutine gxmarg(ixref,iyref,wnx,wny,active)
  use gxx11_common
  implicit none
  integer i,iax,iaxr,ifirst,ilast,ind,intv,iscal,ival,j,k,kset,nax, &
       nchct,nint,nref
  real add,ahi,alo,fact,gap,gapt,hgap,hgapt,hwid,hwidt,txf,vgap,    &
       vgapt,vwid,vwidt,xf
  !***********************************************************************
  !
  !   Purpose: calculates window margins from axes specifications
  !
  !--- Input
  !   ixref      x axis reference numbers of curve sets
  !   iyref      y axis reference numbers of curve sets
  !   wnx(2,i)   lower and upper x value, curve set I
  !   wny(2,i)   lower and upper y value, curve set I
  !--- Output
  !   active     (1...4) = window in NDC that can be used for curves
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: March 3, 1988
  !
  !***********************************************************************

  integer ixref(*),iyref(*)
  real wnx(2,*),wny(2,*)
  real active(4),bmin(2)
  integer iref(mxaxs+myaxs),iapar(mpaxs)
  real border(4)
  logical flag
  character fmt*20,text*300,fmtloc*60
  save fact
  !--- fact includes one character height plus the gap of half that height
  data fact/.0015/
  !--- get viewport ratio
  call gxqrvp(xf)
  do i=1,4
     border(i)=0.
  enddo
  do i=1,2
     bmin(i)=0.
  enddo
  do iax=1,2
     nref=0
     do kset=1,maxset
        if(iax.eq.1) then
           !--- x axis
           iaxr=ixref(kset)
        else
           !--- y axis
           iaxr=iyref(kset)
        endif
        if(iaxr.eq.0) goto 70
        do j=1,nref
           if(iaxr.eq.iref(j)) goto 70
        enddo
        nref=nref+1
        iref(nref)=iaxr
        if(iax.eq.1) then
           nax=mxaxs
        else
           nax=myaxs
        endif
        do i=1,nax
           add=0.
           if(iax.eq.1) then
              !--- x axis
              if(iaxr.ne.ixapar(21,i)) goto 60
              ind=(4-i)/2
              k=4-ind
              fmt=sxform(i)
              text=sxtext(i)
              !--- tick mark expansion
              txf=1.
              do j=1,mpaxs
                 iapar(j)=ixapar(j,i)
              enddo
           else
              !--- y axis
              if(iaxr.ne.iyapar(21,i)) goto 60
              ind=(4-i)/2
              k=2-ind
              fmt=syform(i)
              text=sytext(i)
              !--- tick mark expansion
              txf=xf
              do j=1,mpaxs
                 iapar(j)=iyapar(j,i)
              enddo
           endif
           gap=fact*iapar(7)
           !--- use separate character height for axis text if specified
           if(iapar(13).gt.0.and.iapar(13).le.1000)  then
              gapt=fact*iapar(13)
           else
              gapt=gap
           endif
           !--- get character height and width for hor. and vert. text
           call gxschf(0,1,gap,hgap,hwid)
           call gxschf(0,2,gap,vgap,vwid)
           call gxschf(0,1,gapt,hgapt,hwidt)
           call gxschf(0,2,gapt,vgapt,vwidt)
           !--- tick marks
           if(iapar(2).ne.0.and.iapar(4).ne.0) then
              if(mod(iapar(4),2).eq.ind) add=add+txf*fact*iapar(8)
           endif
           !--- labels
           if(iapar(3).ne.0.and.iapar(5).ne.0) then
              flag=mod(iapar(5),2).eq.ind
              if(iapar(1).eq.0)  then
                 !--- linear scale - use format if given by user, else calculate
                 call gxival(fmt,ival)
                 if(ival.eq.0) then
                    if(iax.eq.1)  then
                       call gxscal(wnx(1,kset),wnx(2,kset),alo,ahi,        &
                            nint)
                    else
                       call gxscal(wny(1,kset),wny(2,kset),alo,ahi,        &
                            nint)
                    endif
                    intv = iapar(2)
                    if (intv .le. 0)  call gxdint(alo,ahi,intv)
                    call gxdfmt(alo,ahi,intv,ival,iscal,fmtloc)
                 endif
              else
                 !--- log. scale - use powers of 10
                 ival=5
              endif
              if(iax.eq.iapar(3)) then
                 !--- labels parallel to axis
                 if(iax.eq.1) then
                    if(flag) add=add+hgap
                    !--- keep minimum border for labels on perpendicular axes
                    bmin(iax)=max(bmin(iax),.5*ival*hwid)
                 else
                    if(flag) add=add+vgap
                    bmin(iax)=max(bmin(iax),.5*ival*vwid)
                 endif
              else
                 !--- labels perpendicular to axis
                 if(iax.eq.1) then
                    if(flag) add=add+ival*vwid
                    !--- keep minimum border for labels on perpendicular axes
                    bmin(iax)=max(bmin(iax),.5*vgap)
                 else
                    if(flag) add=add+ival*hwid
                    !--- keep minimum border for labels on perpendicular axes
                    bmin(iax)=max(bmin(iax),.5*hgap)
                 endif
              endif
           endif
           !--- text - always parallel to axis
           call gxpnbl(text,ifirst,ilast)
           if(iapar(6).ne.0.and.mod(iapar(6),2).eq.ind                 &
                .and.ifirst.ne.0)  then
              !--- add space for one line and extra space for each line separator
              call gxchct(text(ifirst:ilast),sdefnl,nchct)
              if(iax.eq.1) then
                 add=add+(nchct+1)*hgapt
                 if(i.ge.3) add=add+.5*hgapt
              else
                 add=add+(nchct+1)*vgapt
              endif
           endif
           !--- take largest margin
           border(k)=max(border(k),add)
60         continue
        enddo
70      continue
     enddo
  enddo
  do i=1,3,2
     active(i)=max(border(i),bmin((i+1)/2))
     active(i+1)=1.-max(border(i+1),bmin((i+1)/2))
     !--- protect against too large borders
     if(active(i).ge.active(i+1))  then
        active(i)=.4
        active(i+1)=.6
     endif
  enddo
999 end subroutine gxmarg
subroutine gxopen
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: open  terminal workstation
  !
  !   Author: H. Grote / CERN                        date: Feb. 26, 1988
  !                                           last mod: Feb. 26, 1988
  !
  !***********************************************************************

  call gxundf
  if(lacttm.ne.lundef)  then
     if(interm.gt.0)  then
        call wopwk(interm, mconid, itermt)
        call wacwk(interm)
        lacttm=lundef
     endif
     call gxclrw
  endif
end subroutine gxopen
subroutine gxopps(iun, ityp)
  use gxx11_common
  implicit none
  integer imun,ityp,iun
  !***********************************************************************
  !
  !   Purpose: open or close .ps or .eps Output unit
  !
  !--- Input
  !   iun       +- utput unit number, if = 0: close
  !   ityp      type of Output: 113 = eps,
  !           else ps with 114 = portrait, 115 = landscape
  !
  !   Author: H. Grote / CERN                        date: Apr.  6, 1995
  !                                           last mod: Apr. 27, 1995
  !
  !***********************************************************************

  imun = abs(iun)
  call gxwpep(imun, ityp)
end subroutine gxopps
subroutine gxplot(ncurv,nptval,ipxval,ipyval,icvref,xval,yval,ierr)
  use gxx11_common
  implicit none
  integer ierr,ncurv
  !***********************************************************************
  !
  !   Purpose: plots one frame with several curves and axes
  !
  !--- Input
  !   ncurv      number of curves (=ordered sets of (x,y) pairs) to plot
  !   nptval(i)  number of points ((x,y) pairs) in set I
  !   ipxval(i)  first x value of set I in array XVAL
  !   ipyval(i)  first y value of set I in array YVAL
  !   icvref(i)  number of the parameter set to be used for curve I. This value
  !            will be forced into [1,MAXSET].
  !            The x and y axis reference numbers of set I will be taken
  !            from this parameter set. All x and y axes with the
  !            corresponding reference numbers will be (scaled if automatic)
  !            and plotted together with set I.
  !            If no x resp. y axis exists with the reference number
  !            from the parameter set, the curve will be plotted with
  !            automatic scaling, but without x resp. y axis.
  !   xval       array containing the x values for all sets
  !   yval       array containing the y values for all sets
  !--- Output
  !   ierr       0 if everything OK, else GKS error, or
  !            1 : GXINIT not called (initialization)
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  integer nptval(*),ipxval(*),ipyval(*),icvref(*)
  real xval(*),yval(*)
  !
  call gxundf
  !--- exit if not initialized
  if(ltotin.ne.lundef)  then
     ierr=1
     goto 999
  endif
  !--- clear work station(s) if requested
  if(iclflg.gt.0)  call gxclrw
  !--- plot frame
  call gxfram(ncurv,nptval,ipxval,ipyval,icvref,xval,yval,          &
       cvwnwd,ierr)
  !--- plot curves
  call gxcurv(ncurv,nptval,ipxval,ipyval,icvref,xval,yval,          &
       cvwnwd,ierr)
  !
  !--- wait for <CR> if interactive
  !
  call gxwait
999 end subroutine gxplot
subroutine gxplts(np,xp1,yp1)
  use gxx11_common
  implicit none
  integer i,ierror,j,k,nextra,np
  real d,gxcubv,screen,selem,sg,sl,step,xmax,xmin
  !***********************************************************************
  !
  !   Purpose: plots a smoothed polyline (3rd order cubic splines)
  !
  !--- Input
  !   np         number of points
  !   xp1        x values
  !   yp1        y values
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************
  real xp1(*),yp1(*)

  logical curlfl
  save screen, selem
  !--- screen is the dimension of a reasonable screen, SELEM the length
  !   of a curve piece such that it looks smooth
  data screen/30./, selem/.2/
  !
  if(np.le.2.or.np.gt.madim1) goto 70
  xmin=xp1(1)
  xmax=xp1(1)
  curlfl=.false.
  sg=sign(1.,xp1(2)-xp1(1))
  do i=2,np
     curlfl=curlfl.or.sg*xp1(i).le.sg*xp1(i-1)
     xmin=min(xmin,xp1(i))
     xmax=max(xmax,xp1(i))
  enddo
  if(xmax.eq.xmin) goto 999
  if(curlfl)  then
     !--- y is not a unique function of x - spline x and y independetly
     !   as function of s
     s(1)=0.
     do i=2,np
        s(i)=s(i-1)+sqrt((xp1(i)-xp1(i-1))**2+(yp1(i)-yp1(i-1))**2)
     enddo
     !--- step at which extra points should occur
     step=s(np)*selem/screen
     call gxcubi(np,s,xp1,yy1d(1,1),yy2d(1,1),ierror)
     if(ierror.ne.0) goto 70
     call gxcubi(np,s,yp1,yy1d(1,2),yy2d(1,2),ierror)
     if(ierror.ne.0) goto 70
     k=1
     p(1,1)=xp1(1)
     p(1,2)=yp1(1)
     do i=2,np
        !--- number of extra points to be plotted
        nextra=(s(i)-s(i-1))/step
        d=(s(i)-s(i-1))/(nextra+1)
        do j=1,nextra
           k=k+1
           sl=s(i-1)+j*d
           p(k,1)=gxcubv(sl,np,s,xp1,yy1d(1,1),yy2d(1,1))
           p(k,2)=gxcubv(sl,np,s,yp1,yy1d(1,2),yy2d(1,2))
        enddo
        k=k+1
        p(k,1)=xp1(i)
        p(k,2)=yp1(i)
     enddo
  else
     !--- step at which extra points should occur
     step=(xmax-xmin)*selem/screen
     call gxcubi(np,xp1,yp1,yy1d,yy2d,ierror)
     if(ierror.ne.0) goto 70
     k=1
     p(1,1)=xp1(1)
     p(1,2)=yp1(1)
     do i=2,np
        !--- number of extra points to be plotted
        nextra=abs(xp1(i)-xp1(i-1))/step
        d=(xp1(i)-xp1(i-1))/(nextra+1)
        do j=1,nextra
           k=k+1
           p(k,1)=xp1(i-1)+j*d
           p(k,2)=gxcubv(p(k,1),np,xp1,yp1,yy1d,yy2d)
        enddo
        k=k+1
        p(k,1)=xp1(i)
        p(k,2)=yp1(i)
     enddo
  endif
  call gvpl(k,p(1,1),p(1,2))
  goto 999
70 continue
  !
  !--- error condition - not enough, too many, or identical points
  !
  call gvpl(np,xp1,yp1)
999 end subroutine gxplts
subroutine gxplt1(np,xp1,yp1,ac)
  use gxx11_common
  implicit none
  integer i,ierror,j,k,nextra,np
  real d,gxcubv,screen,selem,sg,sl,step,xmax,xmin
  !***********************************************************************
  !
  !   Purpose: plots a smoothed polyline (3rd order cubic splines), clips
  !
  !--- Input
  !   np         number of points
  !   xp1         x values
  !   yp1         y values
  !   ac         active window for clipping ---> routine GXPL
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************
  real xp1(*),yp1(*),ac(4)

  logical curlfl
  save screen, selem
  !--- screen is the dimension of a reasonable screen, SELEM the length
  !   of a curve piece such that it looks smooth
  data screen/30./, selem/.2/
  !
  if(np.le.2.or.np.gt.madim1) goto 70
  xmin=xp1(1)
  xmax=xp1(1)
  curlfl=.false.
  sg=sign(1.,xp1(2)-xp1(1))
  do i=2,np
     curlfl=curlfl.or.sg*xp1(i).le.sg*xp1(i-1)
     xmin=min(xmin,xp1(i))
     xmax=max(xmax,xp1(i))
  enddo
  if(xmax.eq.xmin) goto 999
  if(curlfl)  then
     !--- y is not a unique function of x - spline x and y independetly
     !   as function of s
     s(1)=0.
     do i=2,np
        s(i)=s(i-1)+sqrt((xp1(i)-xp1(i-1))**2+(yp1(i)-yp1(i-1))**2)
     enddo
     !--- step at which extra points should occur
     step=s(np)*selem/screen
     call gxcubi(np,s,xp1,yy1d(1,1),yy2d(1,1),ierror)
     if(ierror.ne.0) goto 70
     call gxcubi(np,s,yp1,yy1d(1,2),yy2d(1,2),ierror)
     if(ierror.ne.0) goto 70
     k=1
     p(1,1)=xp1(1)
     p(1,2)=yp1(1)
     do i=2,np
        !--- number of extra points to be plotted
        nextra=(s(i)-s(i-1))/step
        d=(s(i)-s(i-1))/(nextra+1)
        do j=1,nextra
           k=k+1
           sl=s(i-1)+j*d
           p(k,1)=gxcubv(sl,np,s,xp1,yy1d(1,1),yy2d(1,1))
           p(k,2)=gxcubv(sl,np,s,yp1,yy1d(1,2),yy2d(1,2))
        enddo
        k=k+1
        p(k,1)=xp1(i)
        p(k,2)=yp1(i)
     enddo
  else
     !--- step at which extra points should occur
     step=(xmax-xmin)*selem/screen
     call gxcubi(np,xp1,yp1,yy1d,yy2d,ierror)
     if(ierror.ne.0) goto 70
     k = 1
     p(1,1)=xp1(1)
     p(1,2)=yp1(1)
     do i=2,np
        !--- number of extra points to be plotted
        nextra=abs(xp1(i)-xp1(i-1))/step
        d=(xp1(i)-xp1(i-1))/(nextra+1)
        do j=1,nextra
           k=k+1
           p(k,1)=xp1(i-1)+j*d
           p(k,2)=gxcubv(p(k,1),np,xp1,yp1,yy1d,yy2d)
        enddo
        if (k .eq. madim1)  then
           !--- flush buffer
           call gxpl(k, p(1,1), p(1,2), ac)
           p(1,1) = p(k,1)
           p(1,2) = p(k,2)
           k = 1
        endif
        k=k+1
        p(k,1)=xp1(i)
        p(k,2)=yp1(i)
     enddo
  endif
  if (k .gt. 1)  call gxpl(k,p(1,1),p(1,2),ac)
  goto 999
70 continue
  !
  !--- error condition - not enough, too many, or identical points
  !
  call gxpl(np,xp1,yp1,ac)
999 end subroutine gxplt1
subroutine gxpl(n,x,y,ac)
  use gxx11_common
  implicit none
  integer i,ilow,j,k,n
  real xtol,ytol
  !***********************************************************************
  !
  !   Purpose: plots a polyline, clips at active window.
  !
  !--- Input:
  !   n          no. of points
  !   x          x positions
  !   y          y positions
  !   ac         active window
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: Dec. 9, 1988
  !
  !***********************************************************************
  real x(*),y(*),ac(4)

  xtol = toleps * (ac(2) - ac(1))
  ytol = toleps * (ac(4) - ac(3))
  if (n .gt. 1)  then
     ilow = 1
10   continue
     do  i = ilow, n
        if(x(i) + xtol .lt. ac(1))  goto 30
        if(x(i) - xtol .gt. ac(2))  goto 30
        if(y(i) + ytol .lt. ac(3))  goto 30
        if(y(i) - ytol .gt. ac(4))  goto 30
     enddo
30   continue
     if (i - 1 .gt. ilow)  then
        call gvpl(i - ilow, x(ilow), y(ilow))
     endif
     do  j = max(i, 2), n
        call gxplxx(x(j-1),y(j-1), ac, xp, yp, k)
        if (k .eq. 2)  then
           call gvpl(2, xp, yp)
        endif
        if(x(j) + xtol .lt. ac(1))  goto 40
        if(x(j) - xtol .gt. ac(2))  goto 40
        if(y(j) + ytol .lt. ac(3))  goto 40
        if(y(j) - ytol .gt. ac(4))  goto 40
        ilow = j
        goto 10
40      continue
     enddo
  endif
end subroutine gxpl
subroutine gxplxx(xin, yin, ac, xout,yout, kp)
  implicit none
  integer i,j,k,kp
  real t,xc,xtol,yc,ytol
  !***********************************************************************
  !
  !   Purpose: returns two points inside or on border of active window
  !
  !--- Input:
  !   xin        x positions of Input points
  !   yin        y positions of Input points
  !   ac         active window
  !--- Output
  !   xout       x positions of Output points
  !   yout       y positions of Output points
  !   kp         no. of points to plot (2 if OK, else less)
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: Dec. 3, 1992
  !
  !***********************************************************************
  real toleps
  parameter (toleps = 1.e-5)
  real ac(4), xin(2), yin(2), xout(2), yout(2)
  real x(2), y(2), xr(4), yr(4)

  xtol = toleps * (ac(2) - ac(1))
  ytol = toleps * (ac(4) - ac(3))
  kp = 0
  do  i = 1, 2
     if(xin(i) + xtol .lt. ac(1))  goto 10
     if(xin(i) - xtol .gt. ac(2))  goto 10
     if(yin(i) + ytol .lt. ac(3))  goto 10
     if(yin(i) - ytol .gt. ac(4))  goto 10
     kp = kp + 1
     j = i
     xout(kp) = xin(i)
     yout(kp) = yin(i)
10   continue
  enddo
  if (kp .lt. 2)  then
     if (kp .eq. 0)  then
        !--- both outside
        do  i = 1, 2
           x(i) = xin(i)
           y(i) = yin(i)
        enddo
     else
        x(1) = xin(j)
        y(1) = yin(j)
        x(2) = xin(3-j)
        y(2) = yin(3-j)
     endif
     !--- treat four cases = four sides of ac: low, up, left, right
     k = 0
     if (abs(y(2) - y(1)) .gt. ytol)  then
        t = (ac(3) - y(1)) / (y(2) - y(1))
        if (t .ge. 0. .and. t .lt. 1.)  then
           xc = x(1) + t * (x(2) - x(1))
           if (xc + xtol .ge. ac(1) .and. xc - xtol .le. ac(2))  then
              k = k + 1
              xr(k) = xc
              yr(k) = ac(3)
           endif
        endif
        t = (ac(4) - y(1)) / (y(2) - y(1))
        if (t .ge. 0. .and. t .lt. 1.)  then
           xc = x(1) + t * (x(2) - x(1))
           if (xc + xtol .ge. ac(1) .and. xc - xtol .le. ac(2))  then
              k = k + 1
              xr(k) = xc
              yr(k) = ac(4)
           endif
        endif
     endif
     if (abs(x(2) - x(1)) .gt. xtol)  then
        t = (ac(1) - x(1)) / (x(2) - x(1))
        if (t .ge. 0. .and. t .lt. 1.)  then
           yc = y(1) + t * (y(2) - y(1))
           if (yc + ytol .ge. ac(3) .and. yc - ytol .le. ac(4))  then
              k = k + 1
              yr(k) = yc
              xr(k) = ac(1)
           endif
        endif
        t = (ac(2) - x(1)) / (x(2) - x(1))
        if (t .ge. 0. .and. t .lt. 1.)  then
           yc = y(1) + t * (y(2) - y(1))
           if (yc + ytol .ge. ac(3) .and. yc - ytol .le. ac(4))  then
              k = k + 1
              yr(k) = yc
              xr(k) = ac(2)
           endif
        endif
     endif
     if (kp .eq. 0)  then
        if( k .eq. 2)  then
           do  i = 1, 2
              xout(i) = xr(i)
              yout(i) = yr(i)
           enddo
           kp = 2
        endif
     elseif (k .eq. 1)  then
        xout(2) = xr(1)
        yout(2) = yr(1)
        kp = 2
     elseif (k .gt. 1)  then
        if (abs(xr(1) - xout(1)) .lt. abs(xr(2) - xout(1)))  then
           xout(2) = xr(2)
           yout(2) = yr(2)
           kp = 2
        endif
     endif
  endif
end subroutine gxplxx
subroutine gxpm(n,x,y,ac)
  use gxx11_common
  implicit none
  integer i,iloop,k,n,nup
  real xerr,yerr
  !***********************************************************************
  !
  !   Purpose: plots a marker symbol if inside active window.
  !
  !--- Input:
  !   n          no. of marker symbols
  !   x          x positions
  !   y          y positions
  !   ac         active window
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: Dec. 9, 1988
  !
  !***********************************************************************
  real x(*),y(*),ac(4)

  xerr=1.e-3*(ac(2)-ac(1))
  yerr=1.e-3*(ac(4)-ac(3))
  do  iloop=1,n,madim2
     nup=min(n,iloop+madim2-1)
     k=0
     do  i=iloop,nup
        if(x(i).lt.ac(1)-xerr)  goto 20
        if(x(i).gt.ac(2)+xerr)  goto 20
        if(y(i).lt.ac(3)-yerr)  goto 20
        if(y(i).gt.ac(4)+yerr)  goto 20
        k=k+1
        xp(k)=x(i)
        yp(k)=y(i)
20      continue
     enddo
     if(k.gt.0) call gvpm(k,xp,yp)
  enddo
end subroutine gxpm
subroutine gxpmsw(n,x,y,ac)
  implicit none
  integer n
  !***********************************************************************
  !
  !   Purpose: plots a software (!) marker symbol if inside active window.
  !         this is necessary where the scaling of hardware symbols
  !         does not work (e.g. Apollo with GTS-GRAL).
  !         only symbols 1 to 5 (.+*ox) are supported.
  !
  !--- Input:
  !   n          no. of marker symbols
  !   x          x positions
  !   y          y positions
  !   ac         active window
  !
  !   Author: H. Grote / CERN                          date: May 11, 1989
  !                                               last mod: July 10, 1995
  !
  !***********************************************************************
  real x(*),y(*),ac(4)
  call gxpm(n,x,y,ac)
end subroutine gxpmsw
subroutine gxpnbl(string,ifirst,ilast)
  implicit none
  integer i,ifirst,ilast
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
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************
  character(*)  string
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
999 end subroutine gxpnbl
subroutine gxppow(alabl,ipower)
  implicit none
  integer i1,i2,ipower
  !***********************************************************************
  !
  !   Purpose: plots a power of ten as label
  !
  !--- Input
  !   alabl      text coordinates
  !   ipower     power to plot
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 8, 1988
  !
  !***********************************************************************
  real alabl(4)
  character sdumm*10
  if(ipower.eq.0)  then
     call gxstx(alabl(1),alabl(3),' 1  ')
  elseif(ipower.eq.1)  then
     call gxstx(alabl(1),alabl(3),'10  ')
  else
     sdumm=' '
     write(sdumm,'(I10)')  ipower
     call gxpnbl(sdumm,i1,i2)
     call gxtx(alabl(1),alabl(3),'10<!>'//sdumm(i1:i2)//'<!>')
  endif
end subroutine gxppow
subroutine gxprwn(ixy,ncrv,icvref,nptval,ipcval,ipval,cval,val,nax,iaref,iapar,range,wn)
  use gxx11_common
  implicit none
  integer i,iadd,iauto,iaxr,ic,icurv,ifl,ip,ipc,iscalf,ixy,j,k,kset,&
       nax,ncrv,nextop,npint,npt,nref,nx
  real rmax,rmaxi,rmaxt,rmin,rmini,rmint,tolo
  !***********************************************************************
  !
  !   Purpose: calculates all window values for either x or y
  !
  !--- Input
  !   ixy        flag: 1 = x, 2 = y (only valid x points for y range)
  !   ncrv       number of curves (=ordered sets of (x,y) pairs) to plot
  !   icvref(i)  curve set for curve I
  !   nptval(i)  number of points ((x,y) pairs) in set I
  !   ipcval(i)  first x value of set I in array VAL (for range check)
  !   ipval(i)   first x or y value of set I in array VAL
  !   cval       array containing the x values for all sets (range check)
  !   val        array containing the x or y values for all sets
  !   nax        max. no. of x or y axes
  !   iaref(i)   axis reference number of set I. All x or y axes with
  !            this reference number will be scaled if automatic.
  !   iapar(j,K) axis parameters of x or y axis K
  !   range(j,K) lower and upper limit of axis K
  !--- Output
  !   wn(j,i)    lower and upper window values for set I
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Nov. 18, 1994
  !
  !***********************************************************************

  integer nptval(*), ipcval(*), ipval(*), iaref(*), iapar(mpaxs,*), &
       icvref(*)
  real cval(*), val(*), range(2,*), wn(2,*)
  real rng(2)
  integer iref(mxaxs+myaxs)
  !
  !--- nref counts the number of axis reference numbers
  nref=0
  do icurv=1,ncrv
     kset=max(1,min(maxset,icvref(icurv)))
     iaxr=iaref(kset)
     !--- check whether there are axes with this ref. number
     do j=1,nref
        if(iaxr.eq.iref(j)) goto 80
     enddo
     !--- not yet in list
     nref=nref+1
     iref(nref)=iaxr
     !--- loop over related axes for scaling etc.
     !--- nextop gives the highest option for zero adjustment
     nextop=0
     !--- iscalf gives the highest scaling type request
     iscalf=0
     !--- iauto is 0 for automatic, 1 for hand scaling
     iauto=0
     do i=1,nax
        if(iapar(21,i).eq.iaxr) then
           nextop=max(nextop,iapar(22,i))
           iscalf=max(iscalf,iapar(1,i))
           iauto =max(iauto ,iapar(23,i))
        endif
     enddo
     if(iauto.ne.0)  then
        !--- hand scaling requested - get extrema of all ranges given
        ifl=0
        do i=1,nax
           if(iapar(21,i).eq.iaxr) then
              if(range(1,i).lt.range(2,i)) then
                 if(ifl.eq.0) then
                    ifl=1
                    rng(1)=range(1,i)
                    rng(2)=range(2,i)
                 endif
                 rng(1)=min(rng(1),range(1,i))
                 rng(2)=max(rng(2),range(2,i))
              endif
           endif
        enddo
        if(ifl.eq.0) then
           !--- no valid range found - use automatic scaling
           iauto=0
        endif
     endif
     if(iauto.eq.0) then
        !--- find minima and maxima of coordinates
        ifl=0
        do ic=icurv,ncrv
           k=max(1,min(maxset,icvref(ic)))
           if(iaxr.eq.iaref(k)) then
              !--- find min. and max. x or y values
              npt=nptval(ic)
              ip=ipval(ic)
              ipc = ipcval(ic)
              iadd = 0
              if (ixy .eq. 2)  then
                 !--- take only y values paired with valid x values
                 tolo = toleps * (axwndx(2,kset) - axwndx(1,kset))
                 do  j = 0, npt - 1
                    if (cval(ipc+j) .ge. axwndx(1, kset) - tolo           &
                         .and. cval(ipc+j) .le. axwndx(2, kset) + tolo)                    &
                         then
                       iadd = j
                       goto 32
                    endif
                 enddo
32               continue
              endif
              if(ifl.eq.0) then
                 rmini=val(ip+iadd)
                 rmaxi=rmini
                 ifl=1
              endif
              do j=iadd,npt-1
                 if (ixy .eq. 1)  then
                    rmini=min(rmini,val(ip+j))
                    rmaxi=max(rmaxi,val(ip+j))
                 elseif (cval(ipc+j) .ge. axwndx(1, kset)                &
                      .and. cval(ipc+j) .le. axwndx(2, kset))  then
                    rmini=min(rmini,val(ip+j))
                    rmaxi=max(rmaxi,val(ip+j))
                 endif
              enddo
           endif
        enddo
     endif
     !   nx counts the axes belonging to IAXR
     nx=0
     do i=1,nax
        !--- keep interval number as given by user
        iapar(19,i)=iapar(2,i)
        if(iapar(21,i).eq.iaxr) then
           nx=nx+1
           if(iauto.eq.0) then
              !--- automatic scaling of this axis requested
              call gxarng(nextop,rmini,rmaxi,rmin,rmax,npint)
              if(iapar(2,i).lt.0)  iapar(19,i)=npint
           else
              !--- hand scaling
              rmin=rng(1)
              rmax=rng(2)
              if(nextop.eq.1) then
                 !--- start or end axis at 0. if possible
                 if(rmin.gt.0.) then
                    rmin=0.
                 elseif(rmax.lt.0.) then
                    rmax=0.
                 endif
              elseif(nextop.eq.2) then
                 !--- make axis symmetric around 0.
                 rmax=max(abs(rmin),abs(rmax))
                 rmin=-rmax
              endif
           endif
           !--- keep overall min. and max.
           if(nx.eq.1) then
              rmint=rmin
              rmaxt=rmax
           else
              rmint=min(rmint,rmin)
              rmaxt=max(rmaxt,rmax)
           endif
        endif
     enddo
     if(nx.eq.0) then
        !--- no axis found for this ref. number - automatic scale
        call gxarng(0,rmini,rmaxi,rmint,rmaxt,npint)
     endif
     !--- set windows
     do ic=icurv,ncrv
        k=max(1,min(maxset,icvref(ic)))
        if(iaref(k).eq.iaxr) then
           wn(1,k)=rmint
           wn(2,k)=rmaxt
        endif
     enddo
80   continue
  enddo
999 end subroutine gxprwn
subroutine gxqaxs(type,naxis,npar,ipar,range,stext,sform)
  use gxx11_common
  implicit none
  integer i,naxis,npar
  !***********************************************************************
  !
  !   Purpose: returns axis parameters
  !
  !--- Input
  !   type     'X' for an x-axis, 'Y' for a y-axis
  !   naxis    axis number
  !--- Output
  !   npar     no. of axis parameters in IPAR
  !          or = 0 if NAXIS and/or TYPE are wrong, in which case the other
  !          Output parameters will not be set
  !   ipar     parameter list
  !   range(1) lower axis limit
  !      (2) upper axis limit
  !   stext    axis text
  !   sform    axis label format, e.g. '(F6.2)'
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  character(*)  type,stext,sform
  integer ipar(*)
  real range(2)
  npar=0
  if(type.eq.'X')  then
     if(naxis.gt.0.and.naxis.le.mxaxs) then
        npar=mpaxs
        stext=sxtext(naxis)
        sform=sxform(naxis)
        do i=1,mpaxs
           ipar(i)=ixapar(i,naxis)
        enddo
        do i=1,2
           range(i)=rangex(i,naxis)
        enddo
     endif
  elseif(type.eq.'Y')  then
     if(naxis.gt.0.and.naxis.le.myaxs) then
        npar=mpaxs
        stext=sytext(naxis)
        sform=syform(naxis)
        do i=1,mpaxs
           ipar(i)=iyapar(i,naxis)
        enddo
        do i=1,2
           range(i)=rangey(i,naxis)
        enddo
     endif
  endif
end subroutine gxqaxs
subroutine gxqcrv(nset,npar,ipar,symb)
  use gxx11_common
  implicit none
  integer i,npar,nset
  !***********************************************************************
  !
  !   Purpose: inquire curve set parameters
  !
  !--- Input
  !   nset    curve set number
  !--- Output
  !   npar     number of curve set parameters returned in IPAR
  !   ipar     parameter list
  !   symb     plot symbol
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  character(1) symb
  integer ipar(*)
  npar=0
  if(nset.gt.0.and.nset.le.maxset)  then
     npar=mpcurv
     do i=1,mpcurv
        ipar(i)=icvpar(i,nset)
     enddo
     symb=splotc(nset:nset)
  endif
end subroutine gxqcrv
subroutine gxqrvp(xf)
  implicit none
  integer ict,ierr
  real xf
  !***********************************************************************
  !
  !   Purpose: inquire view port ratio (y to x extension)
  !
  !--- Output
  !   xf       expansion factor
  !
  !   Author: H. Grote / CERN                        date: March 2, 1988
  !                                           last mod: March 2, 1988
  !
  !***********************************************************************
  real w(4),v(4)
  xf=1.
  !--- get current norm. transf. number
  call jqcntn(ierr,ict)
  if(ierr.ne.0) goto 999
  !--- get current window and viewport
  call jqnt(ict,ierr,w,v)
  if(ierr.ne.0) goto 999
  if(v(2).gt.v(1).and.v(4).gt.v(3))  then
     xf=(v(4)-v(3))/(v(2)-v(1))
  endif
999 end subroutine gxqrvp
subroutine gxqvar(name,intv,realv,charv)
  use gxx11_common
  implicit none
  integer intv
  real realv
  !***********************************************************************
  !
  !   Purpose: returns values of certain variables in common GXCOMM
  !
  !--- Input:
  !   name     name of the variable (character):
  !   = itermt   terminal workstation type (default = MTERMT)
  !   = interm   terminal workstation number (default  = MTTERM if
  !            GXASKU called, 0 otherwise for batch)
  !            if = 0, no graphics display on terminal
  !   = inmeta   metafile workstation number  (default = MTMETA)
  !            if = 0, no metafile written
  !   = ierrun   GKS error file unit number (default = MERRUN)
  !   = imetun   metafile unit  (default = METAUN)
  !   = inunit   terminal or default READ unit (default = 5)
  !   = iounit   terminal or default PRINT unit (default = 6)
  !   = isfflg   =0 (default) for square, 1 for full screen area
  !   = isqflg   =0 (default) for independent window optimization in x and y,
  !             =1 for an identical window range in x and y.
  !         this means that if:
  !                            ISFFLG=0, ISQFLG=1
  !                            and the viewport has not been tampered with
  !                            and the x and y scales are identical
  !         then
  !            (on a plotter) a circle will be plotted as a circle (!)
  !                            if GXPLOT is called
  !   = iwtflg   if = 0 (default), no action.
  !         if = 1 (set by GXASKU if interactive), GXPLOT will wait for some
  !         Input from the keyboard (e.g. <CR>) before returning so that you
  !         can look at the picture. The waiting routine GXWAIT can be called
  !         separately.
  !   = iclflg   =0 : no action; = 1 (default): causes a "clear workstations"
  !         at the end of GXPLOT. This is simply done by
  !         if(INTERM.GT.0)  CALL GCLRWK(INTERM,0)
  !         if(INMETA.GT.0)  CALL GCLRWK(INMETA,0)
  !         in case you want to do it separately.
  !   = inormt   normalization transformation number (default=MNORMT)
  !   = ipseps   .ps (1), .eps (2), else no Output
  !   = idinit   treat first GXINIT call as dummy if not zero
  !   = nxpix    x size of window in pixels (X11)
  !   = nypix    y size of window in pixels (X11)
  !   = xmetaf   paper length in cm for metafile plotting
  !   = ymetaf   paper width in cm for metafile plotting
  !            if either XMETAF or YMETAF = 0. (default), then the
  !            default square will be plotted
  !   = serrnm   GKS error file name (default GXFERR)
  !   = smetnm   Metafile name (default GXMETA)
  !   = sdefnl   new line start default in axis titles
  !
  !--- Output:
  !   intv     integer value if the variable is INTEGER
  !   realv    real value if the variable is REAL
  !   charv    if the variable is CHARACTER
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: May 12, 1993
  !
  !***********************************************************************

  character(*) name,charv
  character(6) code
  code=name
  if    (code.eq.'ITERMT')  then
     intv=itermt
  elseif(code.eq.'INTERM')  then
     intv=interm
  elseif(code.eq.'INMETA')  then
     intv=inmeta
  elseif(code.eq.'IERRUN')  then
     intv=ierrun
  elseif(code.eq.'IMETUN')  then
     intv=imetun
  elseif(code.eq.'INUNIT')  then
     intv=inunit
  elseif(code.eq.'IOUNIT')  then
     intv=iounit
  elseif(code.eq.'ISFFLG')  then
     intv=isfflg
  elseif(code.eq.'ISQFLG')  then
     intv=isqflg
  elseif(code.eq.'IWTFLG')  then
     intv=iwtflg
  elseif(code.eq.'ICLFLG')  then
     intv=iclflg
  elseif(code.eq.'INORMT')  then
     intv=inormt
  elseif(code.eq.'IPSEPS')  then
     intv=ipseps
  elseif(code.eq.'IDINIT')  then
     intv=idinit
  elseif(code.eq.'NXPIX')  then
     intv = nxpix
  elseif(code.eq.'NYPIX')  then
     intv = nypix
  elseif(code.eq.'XMETAF')  then
     realv=xmetaf
  elseif(code.eq.'YMETAF')  then
     realv=ymetaf
  elseif(code.eq.'SERRNM')  then
     charv=serrnm
  elseif(code.eq.'SMETNM')  then
     charv=smetnm
  elseif(code.eq.'SDEFNL')  then
     charv=sdefnl
  endif
end subroutine gxqvar
subroutine gxqwac(wact)
  use gxx11_common
  implicit none
  integer i
  !***********************************************************************
  !
  !   Purpose: returns current active user area (inside frame) in NDC
  !
  !--- Output
  !   wact     active window in NDC
  !
  !   Author: H. Grote / CERN                        date: Dec 17, 1987
  !                                           last mod: Dec 17, 1987
  !
  !***********************************************************************

  real wact(4)
  do i=1,4
     wact(i)=actwnd(i)
  enddo
end subroutine gxqwac
subroutine gxrdtx(lrunit,sline,ierr)
  implicit none
  integer ierr,lrunit
  !***********************************************************************
  !
  !   reads a character string from LRUNIT, recovers from
  !   empty carriage return under VM (string = blank).
  !
  !--- Input
  !   lrunit       Input unit
  !--- Output
  !   sline        string read, or blank if <CR>
  !   ierr         = 0 if no error, else = 1
  !
  !   Author hG  11.2.86
  !
  !***********************************************************************
  character sline*(*)
  ierr=0
  sline=' '
  read(lrunit,'(A)',err=30,end=30) sline
  goto 999
30 ierr=1
999 end subroutine gxrdtx
subroutine gxrest(isave,rsave)
  implicit none
  !***********************************************************************
  !
  !   Purpose: restores GKS settings
  !
  !--- Input
  !   isave     integer list of saved values  (GXINIT def. in brackets)
  !       1   norm. transf. number  (1)
  !       2   line style             (1)
  !       3,4 hor. and vert. text alignment  (0,0)
  !       5,6 font and precision  (1,0)
  !         7 text colour index   (1)
  !         8 marker colour index (1)
  !         9 polyline colour index (1)
  !        10 marker type           (3)
  !        11 text path             (0)
  !        12 fill area interior style (0)
  !        13 fill area style index  (if interior style = 2)
  !   rsave     floating list of saved values
  !       1-4 window               (0.,1.,0.,1.)
  !       5-8 viewport             (0.,1.,0.,WFACT)
  !       9   character height   (0.01)
  !     10,11 character up vector  (0.,1.)
  !        12 line width scale factor  (1.)
  !        13 marker scale factor      (1.)
  !        14 character spacing factor (0.)
  !        15 character expansion factor (1.)
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: March 7, 1988
  !
  !***********************************************************************
  integer isave(*)
  real    rsave(*)
  !
  if(isave(1).ne.0)  then
     !     call jselnt(isave(1))
     call jswn(isave(1),rsave(1),rsave(2),rsave(3),rsave(4))
  endif
  call jsln(isave(2))
  call jstxal(isave(3),isave(4))
  call jstxfp(isave(5),isave(6))
  call jstxci(isave(7))
  call jspmci(isave(8))
  call jsplci(isave(9))
  call jsmk(isave(10))
  !  call jstxp(isave(11))
  call jsfais(isave(12))
  !  call jsfasi(isave(13))
  !
  call jschh(rsave(9))
  call jschup(rsave(10),rsave(11))
  call jslwsc(rsave(12))
  call jsmksc(rsave(13))
  !  call jschsp(rsave(14))
  call jschxp(rsave(15))
999 end subroutine gxrest
subroutine gxsave(isave,rsave,ierr)
  implicit none
  integer ierr
  !***********************************************************************
  !
  !   Purpose: saves current GKS settings
  !
  !--- Output
  !   isave     integer list of saved values
  !       1   norm. transf. number
  !       2   line style
  !       3,4 hor. and vert. text alignment
  !       5,6 font and precision
  !         7 text colour index
  !         8 marker colour index
  !         9 polyline colour index
  !        10 marker type
  !        11 text path
  !        12 fill area interior style
  !        13 fill area style index  (if interior style = 2)
  !   rsave     floating list of saved values
  !       1-4 window               (0.,1.,0.,1.)
  !       5-8 viewport             (0.,1.,0.,WFACT)
  !       9   character height   (0.01)
  !     10,11 character up vector  (0.,1.)
  !        12 line width scale factor  (1.)
  !        13 marker scale factor      (1.)
  !        14 character spacing factor (0.)
  !        15 character expansion factor (1.)
  !   ierr      0 if OK, or GKS error number
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: March 7, 1988
  !
  !***********************************************************************
  integer isave(*)
  real    rsave(*)
  !
  call jqcntn(ierr,isave(1))
  if(ierr.ne.0) goto 999
  call jqln(ierr,isave(2))
  if(ierr.ne.0) goto 999
  call jqtxal(ierr,isave(3),isave(4))
  if(ierr.ne.0) goto 999
  call jqtxfp(ierr,isave(5),isave(6))
  if(ierr.ne.0) goto 999
  call jqtxci(ierr,isave(7))
  if(ierr.ne.0) goto 999
  call jqpmci(ierr,isave(8))
  if(ierr.ne.0) goto 999
  call jqplci(ierr,isave(9))
  if(ierr.ne.0) goto 999
  call jqmk(ierr,isave(10))
  if(ierr.ne.0) goto 999
  !  call jqtxp(ierr,isave(11))
  if(ierr.ne.0) goto 999
  call jqfais(ierr,isave(12))
  if(ierr.ne.0) goto 999
  !  call jqfasi(ierr,isave(13))
  if(ierr.ne.0) goto 999
  !
  call jqnt(isave(1),ierr,rsave(1),rsave(5))
  if(ierr.ne.0) goto 999
  call jqchh(ierr,rsave(9))
  if(ierr.ne.0) goto 999
  call jqchup(ierr,rsave(10),rsave(11))
  if(ierr.ne.0) goto 999
  call jqlwsc(ierr,rsave(12))
  if(ierr.ne.0) goto 999
  call jqmksc(ierr,rsave(13))
  if(ierr.ne.0) goto 999
  !  call jqchsp(ierr,rsave(14))
  if(ierr.ne.0) goto 999
  call jqchxp(ierr,rsave(15))
  if(ierr.ne.0) goto 999
999 end subroutine gxsave
subroutine gxsaxs(type,naxis,npar,ipar,range,stext,sform)
  use gxx11_common
  implicit none
  integer i,naxis,npar
  !***********************************************************************
  !
  !   Purpose: set axis parameters
  !
  !--- Input
  !   type     'X' for an x-axis, 'Y' for a y-axis
  !   naxis    axis number
  !   npar     parameters 1 to NPAR will be taken from IPAR
  !   ipar     parameter list
  !   range(1) lower axis limit
  !      (2) upper axis limit
  !   stext    axis text - only plotted if first character not blank
  !   sform    axis label format, e.g. '(F6.2)'
  !
  !----------- REMARK.
  !          no action if TYPE and/or NAXIS are wrong.
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  character(*)  type,stext,sform
  integer ipar(*)
  real range(2)
  if(type.eq.'X')  then
     if(naxis.gt.0.and.naxis.le.mxaxs) then
        sxtext(naxis)=stext
        sxform(naxis)=sform
        do i=1,min(npar,mpaxs)
           ixapar(i,naxis)=ipar(i)
        enddo
        do i=1,2
           rangex(i,naxis)=range(i)
        enddo
     endif
  elseif(type.eq.'Y')  then
     if(naxis.gt.0.and.naxis.le.myaxs) then
        sytext(naxis)=stext
        syform(naxis)=sform
        do i=1,min(npar,mpaxs)
           iyapar(i,naxis)=ipar(i)
        enddo
        do i=1,2
           rangey(i,naxis)=range(i)
        enddo
     endif
  endif
end subroutine gxsaxs
subroutine gxscal(xmin,xmax,xlo,xhi,nint)
  implicit none
  integer i,icase,idistl,mrange,n1,n2,nint,niv
  real x1,x2,xhi,xlo,xmax,xmin
  !***********************************************************************
  !
  !   Purpose: for a given arbitrary interval, finds the nearest interval
  !         in round numbers and a good subdivision
  !
  !--- Input:
  !   xmin          lower corner of interval
  !   xmax          upper corner of interval
  !
  !--- Output:
  !   xlo           lower corner of (possibly) extended interval, round
  !   xhi           upper corner (dito)
  !   nint          proposed number of intervals (zero will be precisely
  !               on the limit of a subinterval if included in range)
  !
  !   Author: H. Grote / CERN                        date: Sept. 7, 1987
  !   last modification:                                   Jan. 12, 1988
  !
  !***********************************************************************
  parameter (mrange=10)
  integer iv(mrange+1)
  double precision rangl(mrange+1)
  double precision xhigh,xlow,ztol1,err1,dist,distl,distn,          &
       distr,sdist
  save iv, rangl, err1
  data iv/10,12,8,10,10,12,8,10,12,8,10/
  data rangl/1.d0,1.2d0,1.6d0,2.d0,2.5d0,3.d0,4.d0,5.d0,            &
       6.d0,8.d0,10.d0/
  data ztol1,err1/5.d-4,1.d-8/
  !
  !--- prepare Input data: consider only the case of at least one positive
  !   value. a boundary too near to zero will be set to zero.
  !
  x1=min(xmin,xmax)
  x2=max(xmin,xmax)
  if(x1.eq.x2)  then
     if(x1.lt.0.)  then
        x2=0.
     elseif(x1.gt.0.)  then
        x1=0.
     else
        x1=0.
        x2=1.
     endif
  endif
  if(x2.ge.0.)  then
     icase=1
     xlow=x1
     xhigh=x2
  else
     icase=2
     xlow=-x2
     xhigh=-x1
  endif
  if(abs(xlow).lt.xhigh)  then
     if(abs(xlow)/xhigh.lt.2.d0*ztol1)  xlow=0.d0
  else
     if(xhigh/abs(xlow).lt.2.d0*ztol1)  xhigh=0.d0
  endif
  if(icase.eq.1.and.xhigh.eq.0.d0)  then
     icase=2
     xhigh=-xlow
     xlow=0.d0
  endif
  !
  !--- choose a reasonable (round) distance
  !
  distl=log10(xhigh-xlow)
  idistl=distl
  if(distl.lt.0.d0)  idistl=idistl-1
  distr=10.d0**idistl
  distn=(xhigh-xlow)/distr
  do i=1,mrange
     if(distn.le.rangl(i)+err1) goto 20
  enddo
20 continue
  dist=rangl(i)*distr
  if(xlow.eq.0.d0)  then
     !
     !--- first case: one of the boundary values is zero
     !
     nint=iv(i)
     if(icase.eq.1)  then
        xlo=0.
        xhi=dist
     else
        xlo=-dist
        xhi=0.
     endif
  elseif(xlow.lt.0.d0)  then
     !
     !--- second case: zero is included in interval and should therefore come
     !   to coincide with a subinterval boundary
     !
     niv=iv(i)
     sdist=dist/niv
     n1=abs(xlow)/sdist+1.d0-ztol1
     n2=xhigh/sdist+1.d0-ztol1
     nint=n1+n2
     xlo=-n1*sdist
     xhi=n2*sdist
  else
     !
     !--- both boundaries are below or above zero
     !
     niv=iv(i)
     sdist=dist/niv
     n1=xlow/sdist+ztol1
     n2=xhigh/sdist+1.d0-ztol1
     nint=n2-n1
     xlow=n1*sdist
     xhigh=n2*sdist
     if(icase.eq.1)  then
        xlo=xlow
        xhi=xhigh
     else
        xlo=-xhigh
        xhi=-xlow
     endif
  endif
end subroutine gxscal
subroutine gxschf(imode,iort,ch,chret,chwid)
  implicit none
  integer ierr,iffo,imode,iort,ippr
  real ch,chhxf,chret,chwid,hxf,xf
  !***********************************************************************
  !
  !   Purpose: set character height and orientation with correct scales
  !
  !--- Input
  !   imode    if 0, inquire only. If > 0, set height and exp. factor
  !   iort     text orientation: 1 horizontal, 2 vertical
  !   ch       character height for a standard viewport
  !--- Output
  !   chret    character height actually set
  !   chwid    character width
  !
  !   Author: H. Grote / CERN                          date: March 2, 1988
  !                                                last mod: May 16,  1995
  !
  !***********************************************************************
  !--- set height expansion factor if font
  call jqtxfp(ierr, iffo, ippr)
  if (imode .eq. 0 .and. iffo .lt. 0)  then
     chhxf = 1.5
  else
     chhxf = 1.
  endif
  call gxqrvp(hxf)
  if (iort .eq. 1)  then
     chret = ch * chhxf
     xf = hxf
     if(imode .gt. 0) call jschup(0., 1.)
  else
     chret = hxf * ch
     xf = 1. / hxf
     if(imode .gt. 0) call jschup(-1., 0.)
  endif
  chwid = .9 * chret * xf
  if (imode .gt. 0)  then
     call jschh(chret)
     call jschxp(xf)
  endif
999 end subroutine gxschf
subroutine gxscol(icol)
  use gxx11_common
  implicit none
  integer icol
  !***********************************************************************
  !
  !   Purpose: set foreground colour
  !
  !--- Input
  !   icol     requested colour
  !
  !   Author: H. Grote / CERN                        date: Apr. 7, 1995
  !                                           last mod: Apr. 7, 1995
  !
  !***********************************************************************

  !      character chst * 16
  icucol = max(1, mod(icol - 1, mcolor) + 1)
end subroutine gxscol
subroutine gxscrv(nset,npar,ipar,symb)
  use gxx11_common
  implicit none
  integer i,npar,nset
  !***********************************************************************
  !
  !   Purpose: set curve set parameters
  !
  !--- Input
  !   nset     curve set number
  !   npar     parameters 1 to NPAR will be taken from IPAR
  !   ipar     parameter list
  !   symb     plot symbol
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  character(1) symb
  integer ipar(*)
  if(nset.gt.0.and.nset.le.maxset)  then
     do i=1,min(npar,mpcurv)
        icvpar(i,nset)=ipar(i)
     enddo
     splotc(nset:nset)=symb
  endif
end subroutine gxscrv
subroutine gxsdef(sitem,item)
  use gxx11_common
  implicit none
  integer i,i1,i2,item,j,k
  !***********************************************************************
  !
  !   Purpose: sets undefined variables or restores default values
  !
  !--- Input
  !   sitem    selects the set of values to be defined:
  !          'OPTION' for GKS options other than curves and axes: defaults
  !          'OPTINI' as for OPTION, but only undefined values are set
  !          'CURVE'  for curve set ITEM: set defaults
  !          'AXIS'   for both x and y axis no. ITEM: set defaults
  !          'XAXIS'  for x axis no. ITEM
  !          'YAXIS'  for y axis no. ITEM
  !          'DEVICE' for viewport: set defaults
  !   item     curve set or axis number. If = 0 , all curve sets or axes
  !   iflag    if = 0, define only undefined. if = 1, reset to defaults
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  character(*)  sitem
  character titem*50,sxtdef*60,sytdef*60,sxfdef*20,syfdef*20
  character spldef*(maxset)
  integer ixadef(mpaxs), iyadef(mpaxs), icvdef(mpcurv)
  real vp(4), vpdef(4), rxdef(2), rydef(2)
  real rgb(3,mcolor)
  character col(mcolor) * 16
  logical defaul,xaxis,yaxis
  save spldef, sxtdef, sxfdef, vp, ixadef, iyadef, icvdef,          &
       rxdef, rydef, rgb, col
  data rgb /3*0., 1.,0.,0., 0.,1.,0., 1.,0.,1., 0.,1.,1., 1.,1.,0./
  data col / 'black', 'red', 'green', 'blue', 'cyan',               &
       'magenta' /
  data spldef/'********************'/
  data sxtdef/' '/, sytdef/' '/
  data sxfdef/' '/, syfdef/' '/
  data vp / 0., 1., 0., 1. /
  data  vpdef/ 0.05, 0.95 ,0.05, 0.95 /
  data ixadef/0,-1,1,-99,-99,-99,10,10,2,2,1,9*0,                   &
       -99,0,0/
  data iyadef/0,-1,1,-99,-99,-99,10,10,2,2,1,9*0,                   &
       -99,0,0/
  data icvdef/1,1,1,1,0,1,0,0,0,10/
  data rxdef/0.,0./, rydef/0.,0./

  call gxundf
  titem=sitem
  !------------------------------------------------
  !
  !    options
  !
  !------------------------------------------------
  if(titem(:6).eq.'OPTION'.or.titem(:6).eq.'OPTINI')  then
     defaul=titem(:6).eq.'OPTION'
     spsnam = ' '
     if(lnunit.ne.lundef.or.defaul) then
        lnunit=lundef
        inunit=miunit
     endif
     if(lounit.ne.lundef.or.defaul) then
        lounit=lundef
        iounit=mounit
     endif
     if(lnormt.ne.lundef.or.defaul) then
        lnormt=lundef
        inormt=mnormt
     endif
     if(lmetun.ne.lundef.or.defaul) then
        lmetun=lundef
        imetun=metaun
     endif
     if(lerrun.ne.lundef.or.defaul) then
        lerrun=lundef
        ierrun=merrun
     endif
     if(lmetax.ne.lundef.or.defaul) then
        lmetax=lundef
        xmetaf = mxsize
     endif
     if(lmetay.ne.lundef.or.defaul) then
        lmetay=lundef
        ymetaf = mysize
     endif
     if(ltermt.ne.lundef.or.defaul) then
        ltermt=lundef
        itermt=mtermt
     endif
     if(lnterm.ne.lundef.or.defaul) then
        lnterm=lundef
        interm=0
     endif
     if(lsfflg.ne.lundef.or.defaul) then
        lsfflg=lundef
        isfflg=0
     endif
     if(lsqflg.ne.lundef.or.defaul) then
        lsqflg=lundef
        isqflg=0
     endif
     if(lwtflg.ne.lundef.or.defaul) then
        lwtflg=lundef
        iwtflg=0
     endif
     if(lclflg.ne.lundef.or.defaul) then
        lclflg=lundef
        iclflg=1
     endif
     if(lnmeta.ne.lundef.or.defaul) then
        lnmeta=lundef
        inmeta=mtmeta
     endif
     if(lpseps.ne.lundef.or.defaul) then
        lpseps = lundef
        ipseps = 1
     endif
     if(ldinit.ne.lundef.or.defaul) then
        ldinit = lundef
        idinit = 0
     endif
     if(lxpix.ne.lundef.or.defaul) then
        lxpix = lundef
        nxpix = mxpix
     endif
     if(lypix.ne.lundef.or.defaul) then
        lypix = lundef
        nypix = mypix
     endif
     if(lerrnm.ne.lundef.or.defaul) then
        lerrnm=lundef
        serrnm='GXFERR'
     endif
     if(lmetnm.ne.lundef.or.defaul) then
        lmetnm=lundef
        smetnm='GXMETA'
     endif
     if(ldefnl.ne.lundef.or.defaul) then
        lmetnm=lundef
        sdefnl='//'
     endif
     if(lmpict.ne.lundef.or.defaul) then
        lmpict=lundef
     endif
     if(lttime.ne.lundef.or.defaul) then
        lttime=lundef
        !--- wait time for Higz in sec.
        wttime=0.5
     endif
     icucol = 1
     do  k = 1, mcolor
        do  i = 1, 3
           rgbcol(i,k) = rgb(i,k)
        enddo
        colour(k) = col(k)
     enddo
     !------------------------------------------------
     !
     !    device
     !
     !------------------------------------------------
  elseif(titem(:6).eq.'DEVICE') then
     !--- set normalization transformation
     !     call jselnt(inormt)
     do  i = 1, 4
        vptdef(i) = vpdef(i)
     enddo
     call gxsvpt(vp)
     !------------------------------------------------
     !
     !    curve sETS
     !
     !------------------------------------------------
  elseif(titem(:5).eq.'CURVE') then
     if(item.eq.0) then
        i1=1
        i2=maxset
     else
        i1=item
        i2=item
     endif
     do i=i1,i2
        do j=1,mpcurv
           icvpar(j,i)=icvdef(j)
        enddo
     enddo
     splotc(i1:i2)=spldef(i1:i2)
     !------------------------------------------------
     !
     !    axes
     !
     !------------------------------------------------
  else
     xaxis=titem(:3).eq.'XAX'.or.titem(:4).eq.'AXIS'
     yaxis=titem(:3).eq.'YAX'.or.titem(:4).eq.'AXIS'
     if(xaxis)  then
        if(item.eq.0) then
           i1=1
           i2=mxaxs
        else
           i1=item
           i2=item
        endif
        do i=i1,i2
           rangex(1,i)=rxdef(1)
           rangex(2,i)=rxdef(2)
           sxtext(i)=sxtdef
           sxform(i)=sxfdef
           do j=1,mpaxs
              if(ixadef(j).eq.-99) then
                 ixapar(j,i)=i
              else
                 ixapar(j,i)=ixadef(j)
              endif
           enddo
        enddo
     endif
     if(yaxis)  then
        if(item.eq.0) then
           i1=1
           i2=myaxs
        else
           i1=item
           i2=item
        endif
        do i=i1,i2
           rangey(1,i)=rydef(1)
           rangey(2,i)=rydef(2)
           sytext(i)=sytdef
           syform(i)=syfdef
           do j=1,mpaxs
              if(iyadef(j).eq.-99) then
                 iyapar(j,i)=i
              else
                 iyapar(j,i)=iyadef(j)
              endif
           enddo
        enddo
     endif
  endif
end subroutine gxsdef
subroutine gxsfop(fnparm,status,ierr)
  use gxx11_common
  implicit none
  integer ierr,ifirst,lgth
  !***********************************************************************
  !
  !   Purpose: opens the metafile or error file unit
  !
  !--- Input
  !   fnparm    selection parameter: starting 'MET' for metafile,
  !           'ERR' for error file, 'PSF' for Postscript, 'EPS' for
  !           encapsulated Postscript
  !   status    if ne ' ' then the unit will be opened with this,
  !           else 'UNKNOWN' will be used as default
  !--- Output
  !   ierr      =0 if everything OK
  !           =1 if FNPARM invalid
  !           =2 if the corresponding unit is not defined
  !
  !   Author: H. Grote / CERN                        date: Sept. 9, 1987
  !                                           last mod: May 12, 1993
  !
  !***********************************************************************

  character(*) fnparm,status
  character smlocn * 80, sub*3, stat*7
  call gxundf
  sub=fnparm
  ierr=0
  !--- choose positive unit number (see higz)
  imetps = abs(imetun)
  if(status(:1).eq.' ')  then
     stat='UNKNOWN'
  else
     stat=status
  endif
  if(sub.eq.'MET')  then
     if(lmetun.eq.lundef) then
        if(lmetnm.eq.lundef) then
           lgth=len(smetnm)
           open(unit=imetps,file=smetnm(:lgth),status=stat)
        else
           open(unit=imetps,status=stat)
        endif
        lmetop=lundef
     else
        ierr=2
     endif
  elseif(sub .eq. 'EPS')  then
     if(lmetun .eq. lundef) then
        if(lmetnm .eq. lundef) then
           iepscf = iepscf + 1
           smlocn = smetnm
           call gxpnbl(smlocn, ifirst, lgth)
           write(smlocn(lgth+1:lgth+2), '(I2.2)')  iepscf
           spsnam = smlocn(:lgth+2) // '.eps'
           open(unit = imetps, file = spsnam, status = stat)
        else
           open(unit = imetps, status = stat)
        endif
        iepsop = 2
        lmetop = lundef
     else
        ierr=2
     endif
  elseif (sub .eq. 'PSF')  then
     if(lmetun.eq.lundef) then
        if(lmetnm .eq. lundef) then
           call gxpnbl(smetnm, ifirst, lgth)
           spsnam = smetnm(:lgth) // '.ps'
           open(unit = imetps, file = spsnam, status = stat)
        else
           open(unit = imetps, status = stat)
        endif
        iepsop = 1
        lmetop = lundef
     else
        ierr=2
     endif
  elseif(sub.eq.'ERR')  then
     if(lerrun.eq.lundef) then
        if(lerrnm.eq.lundef) then
           lgth=len(serrnm)
           open(unit=ierrun,file=serrnm(:lgth),status=stat)
        else
           open(unit=ierrun,status=stat)
        endif
        lerrop=lundef
     else
        ierr=2
     endif
  else
     ierr=1
  endif
end subroutine gxsfop
subroutine gxspmt
  implicit none
  !***********************************************************************
  !
  !   Purpose: sets defaults for line style, marker, text, and colour
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************
  call jschup(0.,1.)
  call jsln(1)
  call jstxal(0,0)
  call jstxfp(1,2)
  call jstxci(1)
  call jspmci(1)
  call jsplci(1)
  call jsmksc(.5)
end subroutine gxspmt
subroutine gxsvar(name,intv,realv,charv)
  use gxx11_common
  implicit none
  integer intv
  real realv
  !***********************************************************************
  !
  !   Purpose: sets selected variables for plotting (in common GXCOMM)
  !
  !--- Input:
  !   name     name of the variable (character):
  !   = itermt   terminal workstation type (default = MTERMT)
  !   = interm   terminal workstation number (default  = MTTERM if
  !            GXASKU called, 0 otherwise for batch)
  !            if = 0, no graphics display on terminal
  !   = inmeta   metafile workstation number  (default = MTMETA)
  !            if = 0, no metafile written
  !   = iczebr   call (with any value) sets to LUNDEF (no more MZEBRA
  !             call)
  !   = wttime   Higz wait time in sec. before plotting into a window
  !   = ierrun   GKS error file unit number (default = MERRUN)
  !   = imetun   metafile unit  (default = METAUN)
  !   = inunit   terminal or default READ unit (default = 5)
  !   = iounit   terminal or default PRINT unit (default = 6)
  !   = isfflg   =0 (default) for square, 1 for full screen area
  !   = isqflg   =0 (default) for independent window optimization
  !             in x and y,
  !             =1 for an identical window range in x and y.
  !         this means that if:
  !                            ISFFLG=0, ISQFLG=1
  !                            and the viewport has not been tampered
  !                            with and the x and y scales are identical
  !         then
  !            (on a plotter) a circle will be plotted as a circle (!)
  !                            if GXPLOT is called
  !   = iwtflg   if = 0 (default), no action.
  !         if = 1 (set by GXASKU if interactive), GXPLOT will wait
  !         for some Input from the keyboard (e.g. <CR>) before
  !         returning so that you can look at the picture. The
  !         waiting routine GXWAIT can be called separately.
  !   = iclflg   =0 : no action; = 1 (default): causes a
  !         "clear workstations"
  !         at the end of GXPLOT. This is simply done by
  !         if(INTERM.GT.0)  CALL GCLRWK(INTERM,0)
  !         if(INMETA.GT.0)  CALL GCLRWK(INMETA,0)
  !         in case you want to do it separately.
  !   = inormt   normalization transformation number (default=MNORMT)
  !   = ipseps   .ps (1), .eps (2), else no Output
  !   = idinit   treat first GXINIT call as dummy if not zero
  !   = nxpix    x size of window in pixels (X11)
  !   = nypix    y size of window in pixels (X11)
  !   = xmetaf   paper length in cm for metafile plotting
  !   = ymetaf   paper width in cm for metafile plotting
  !            if either XMETAF or YMETAF = 0. (default), then the
  !            default square will be plotted
  !   = serrnm   GKS error file name (default GXFERR)
  !   = smetnm   Metafile name (default GXMETA)
  !   = sdefnl   new line start default in axis titles
  !
  !   intv     integer value if the variable is INTEGER
  !   realv    real value if the variable is REAL
  !   charv    if the variable is CHARACTER
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: May 24, 1993
  !
  !***********************************************************************

  character(*) name,charv
  character(6) code
  call gxundf
  code=name
  if    (code.eq.'ITERMT')  then
     itermt=intv
     ltermt=lundef
  elseif(code.eq.'INTERM')  then
     interm=intv
     lnterm=lundef
  elseif(code.eq.'INMETA')  then
     inmeta=intv
     lnmeta=lundef
  elseif(code.eq.'ICZEBR')  then
     iczebr=lundef
  elseif(code.eq.'WTTIME')  then
     wttime = realv
     lttime = lundef
  elseif(code.eq.'IERRUN')  then
     ierrun=intv
     lerrun=lundef
  elseif(code.eq.'IMETUN')  then
     imetun=intv
     lmetun=lundef
  elseif(code.eq.'INUNIT')  then
     inunit=intv
     lnunit=lundef
  elseif(code.eq.'ITSEOP')  then
     itseop=intv
     ltseop=lundef
  elseif(code.eq.'IOUNIT')  then
     iounit=intv
     lounit=lundef
  elseif(code.eq.'ISFFLG')  then
     isfflg=intv
     lsfflg=lundef
  elseif(code.eq.'ISQFLG')  then
     isqflg=intv
     lsqflg=lundef
  elseif(code.eq.'IWTFLG')  then
     iwtflg=intv
     lwtflg=lundef
  elseif(code.eq.'ICLFLG')  then
     iclflg=intv
     lclflg=lundef
  elseif(code.eq.'INORMT')  then
     inormt=intv
     lnormt=lundef
  elseif(code.eq.'IPSEPS')  then
     ipseps=intv
     lpseps=lundef
  elseif(code.eq.'IDINIT')  then
     idinit=intv
     ldinit=lundef
  elseif(code.eq.'NXPIX')  then
     nxpix = intv
     lxpix = lundef
  elseif(code.eq.'NYPIX')  then
     nypix = intv
     lypix = lundef
  elseif(code.eq.'XMETAF')  then
     xmetaf=realv
     lmetax=lundef
  elseif(code.eq.'YMETAF')  then
     ymetaf=realv
     lmetay=lundef
  elseif(code.eq.'SERRNM')  then
     serrnm=charv
     lerrnm=lundef
  elseif(code.eq.'SMETNM')  then
     smetnm=charv
     lmetnm=lundef
  elseif(code.eq.'SDEFNL')  then
     sdefnl=charv
     ldefnl=lundef
  endif
end subroutine gxsvar
subroutine gxstep
  use gxx11_common
  implicit none
  integer ierr
  !***********************************************************************
  !
  !   Purpose: opens .eps file
  !
  !   Author: H. Grote / CERN                        date: May 12, 1993
  !                                           last mod: May 12, 1993
  !
  !***********************************************************************

  if (iepsop .eq. -1)  then
     call gxsfop('PSFILE','UNKNOWN',ierr)
  else
     call gxsfop('EPSFILE','UNKNOWN',ierr)
  endif
  call gxopps(imetun, ipstyp)
end subroutine gxstep
subroutine gxstx(xpch, ypch, ch)
  use gxx11_common
  implicit none
  integer i,ich,ie,ierr,ifont,ifttmp,ihort,inttmp,iprec,ipstmp,     &
       ivert,k,lch,np
  real chsize,chux,chuy,cosa,enorm,exfact,sina,wid,width,xcoord,    &
       xmult,xpch,xshift,ycoord,ymult,ypch,yshift
  !***********************************************************************
  !
  !   Purpose: writes a software character string for fonts 1, -13
  !
  !--- Input:
  !   xpch     x position
  !   ypch     y position
  !   ch       string
  !
  !   Author: H. Grote / CERN                        date: Dec. 14, 1990
  !                                           last mod: May 13, 1993
  !**********************************************************************
  character(*) ch

  real xp1(200), yp1(200), xpl(200), ypl(200), rsave(20)
  integer isave(20), ipen(200)
  real yfact(5)
  data yfact / 1.185, 1.,0.5, 0., -0.315 /
  wid=0.
  do i=1,200
     xp1(i)=0.
     yp1(i)=0.
     xpl(i)=0.
     ypl(i)=0.
  enddo
  do i=1,20
     rsave(i)=0.
  enddo
  !--- keep (e)ps flag
  ipstmp = ipseps
  !--- open .eps file if requested
  if (iepsop .lt. 0) call gxstep
  call jqtxfp(ie, ifont, iprec)
  if (ie .ne. 0) goto 999
  !--- start mod 7.6.96
  !    always write .ps or .eps file with a font, sw-characters only
  !    on screen
  if (ifont .ne. 1 .and. ifont .ne. -13)  then
     call gvtx(xpch, ypch, ch)
     goto 999
  else
     !--- switch terminal off, write file
     inttmp = interm
     interm = 0
     if (ifont .eq. 1)  then
        ifttmp = -1
     else
        ifttmp = -12
     endif
     call jstxfp(ifttmp, 2)
     call gvtx(xpch, ypch, ch)
     call jstxfp(ifont, 2)
     interm = inttmp
     !--- switch output file off
     ipseps = 0
     !--- end mod 7.6.96
  endif
  call gxsave(isave, rsave, ierr)
  call jsln(1)
  lch = len(ch)
  ihort = isave(3)
  ivert = isave(4)
  if (ihort .eq. 0) ihort = 1
  if (ivert .eq. 0) ivert = 4
  chsize = rsave(9)
  exfact = rsave(15)
  chux = rsave(10)
  chuy = rsave(11)
  enorm = 1. / sqrt(chux**2 + chuy**2)
  sina = -chux * enorm
  cosa = chuy * enorm
  ymult = chsize / 0.22
  xmult = exfact * ymult
  width = 0.
  ierr = 0
  do ich = 1, lch
     call gxfchr(0, ch(ich:ich), ifont, wid, np, ipen, xp1, yp1, ie)
     ierr = ierr + ie
     width = width + wid
  enddo
  if (ierr .eq. 0)  then
     xshift = 0.5 * (1 - ihort) * width * xmult
     yshift = - chsize * yfact(ivert)
     do ich = 1, lch
        call gxfchr(1, ch(ich:ich), ifont, wid, np, ipen, xp1, yp1,     &
             ierr)
        k = 0
        do i = 1, np
           if (ipen(i) .eq. 0)  then
              !--- pen up
              if (k .gt. 1)  call gvpl(k, xpl, ypl)
              k = 1
           else
              k = k + 1
           endif
           xcoord = xmult * xp1(i) + xshift
           ycoord = ymult * yp1(i) + yshift
           xpl(k) = xpch + cosa * xcoord - sina * ycoord
           ypl(k) = ypch + cosa * ycoord + sina * xcoord
        enddo
        if (k .gt. 1)  call gvpl(k, xpl, ypl)
        xshift = xshift + wid * xmult
     enddo
  else
     call gvtx(xpch, ypch, ch)
  endif
  ipseps = ipstmp
  call gxrest(isave, rsave)
999 end subroutine gxstx
subroutine gxsvpt(vp)
  use gxx11_common
  implicit none
  integer i
  real fdx,fdy
  !***********************************************************************
  !
  !   Purpose: sets workstation viewport. This needs a separate routine
  !         because of the screen size ratio (see below).
  !
  !--- Input
  !   vp       xlow, xup, ylow, yup of viewport.
  !          the factor of the screen ratio WFACT is applied to the
  !          gSVP call directly. The biggest possible screen size
  !          (square or full) is therefore always obtained with
  !          vP = (0.,1.,0.,1.) .
  !          in addition, if a metafile viewport other than the default
  !          is specified, this overrides the screen values.
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Nov 5, 1987
  !
  !***********************************************************************

  real vp(4)
  fdx = vptdef(2) - vptdef(1)
  fdy = vptdef(4) - vptdef(3)
  do  i = 1, 2
     vploc(i) = vptdef(1) + fdx * vp(i)
     vploc(i+2) = vptdef(3) + fdy * vp(i+2)
  enddo
  vpfacx = vploc(2) - vploc(1)
  vpfacy = vploc(4) - vploc(3)
  call jsvp(1, vploc(1), vploc(2), vploc(3), vploc(4))
end subroutine gxsvpt
subroutine gxswnd(window)
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: sets the window for INORMT
  !
  !--- Input
  !   window    vector of length four (xlow, xup, ylow, yup)
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: June 16, 1987
  !
  !***********************************************************************

  real window(4)
  call jswn(inormt,window(1),window(2),window(3),window(4))
end subroutine gxswnd
subroutine gxterm
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: terminates GKS PLOT package
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************

  call gxundf
  if(ltotin.eq.lundef)  then
     if(interm.gt.0.and.lacttm.eq.lundef)  then
        call wdawk(interm)
        call wclwk(interm)
     endif
     if(inmeta.gt.0)  then
        if (iepsop .gt. 0)  then
           call gxopps(0, 0)
        endif
     endif
     call wclks
     if (iepsop .eq. 1 .or. iepsop .eq. 2)  then
        close(imetps)
        iepsop = - iepsop
     endif
     lmetop = 0
  endif
end subroutine gxterm
subroutine gxtint
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   set defaults for a few variables - first routine to call
  !
  !***********************************************************************

  integer i, j
  do i = 1, maxset
     do j = 1, 2
        axwndx(j,i) = 0
        axwndy(j,i) = 0
     enddo
     do j = 1, 4
        cvwnwd(j,i) = 0
     enddo
  enddo
  do i = 1, mxaxs
     do j = 1, 2
        rangex(j,i) = 0
        rangey(j,i) = 0
     enddo
  enddo
  smetnm = ' '
  serrnm = ' '
  do i = 1, mxaxs
     sxtext(i) = ' '
     sxform(i) = ' '
  enddo
  do i = 1, myaxs
     sytext(i) = ' '
     syform(i) = ' '
  enddo
  splotc = ' '
  stortx = ' '
  sdefnl = ' '
end subroutine gxtint
subroutine gxtx(xt,yt,text)
  implicit none
  integer i,ierr,ifdef,ifloc,ifont,igts,ihor,ihoru,in,ip,ipass,ipk, &
       iprec,ivert,ivertu,k,kfpos,l,last,ln,lnk,mchar,mfont,ms
  real ax,axup,ay,ayup,cang,chf,chh,chux,chuy,chxp,crf,ctf,cwf,cxl, &
       cyl,f,falign,sang,sq,t,x,xlift,xt,y,ylift,yt
  !***********************************************************************
  !
  !   Purpose: plots mixture of Roman and Greek text, superscripts
  !         and subscripts. Arguments exactly as GTX.
  !
  !--- Input:
  !   xt       text x position (as for GTX)
  !   yt       text y position (as for GTX)
  !   text     text (as for GTX)
  !          strings included in <G>...<G> appear as Greek,
  !          strings included in <!>...<!> appear as superscripts,
  !          strings included in <?>...<?> as subscripts.
  !          example:  TEXT='<G>a<!>b<!><G>'
  !          gives "alpha to the power of beta"
  !
  !   Author: H. Grote / CERN                        date: June 7, 1988
  !                                           last mod: May 13, 1993
  !
  !***********************************************************************
  parameter (ms=3,mchar=95,mfont=2)
  character text*(*),search(ms)*3,schar*(mchar),stemp*1
  integer ls(ms),is(ms)
  logical toggle(ms)
  real chgtsw(mchar,mfont)
  save crf, is, ls, chf, schar
  !--- set roman as default font
  !   data ifdEF/1/
  data crf/.5/
  data is/0,0,-13/
  data ls/3,3,3/
  data search/'<!>','<?>','<G>'/
  data chf/1.5/
  !--- list all possible keybord characters
  data schar/                                                       &
       ' 1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!@#$%~&*()_+-={}[]:"|;''$><,.?,./'/
  !--- character widths for Roman and Greek fonts.
  data chgtsw/.7273,.9091,.9091,.9091,.9091,.9091,.9091,.9091,.9091,&
       .9091,.9091,.8182,.9545,.9545,.9545,.8636,.8182,.9545,1.,.3636,   &
       .7273,.9545,.7727,1.0909,1.,1.,.9545,1.,.9545,.9091,              &
       .7273,1.,.8182,1.0909,.9091,.8182,.9091,.8636,.8636,.8182,        &
       .8636,.8182,.5455,.8636,.8636,.3636,.4545,.7727,.3636,1.3636,     &
       .8636,.8636,.8636,.8636,.5909,.7727,.5455,.8636,.7273,1.,         &
       .7727,.7273,.7727,.4545,1.2273,.9545,.9091,1.0909,1.,1.1818,      &
       .7273,.6364,.6364,1.1818,1.1818,1.1818,1.1818,.6364,.6364,.6364,  &
       .6364,.4545,.7273,.3636,.4545,.4545,1.,1.0909,1.0909,.4545,       &
       .4545,.8182,.4545,.4545,1.,                                       &
       .7273,.9091,.9091,.9091,.9091,.9091,.9091,.9091,.9091,.9091,      &
       .9091,.8182,.9545,.9091,.8182,.8636,.9091,.7727,1.,.3636,         &
       .6818,.9545,.8182,1.0909,1.,.9091,1.,.6818,.9545,.8182,           &
       .7273,1.,.8182,1.,.8182,1.,.9091,.9545,.8636,.8182,               &
       .8182,.7273,1.,.8636,.9091,.5000,.6818,.8182,.7273,.9545,         &
       .8182,1.0455,1.,.6818,.8182,.9091,.9091,.7727,.9091,1.0455,       &
       .7273,.9545,.6818,.4545,1.2273,.9545,.9091,1.0909,1.,1.1364,      &
       .7273,.6364,.6364,1.1364,1.1364,1.1364,1.1364,.6364,.6364,.6364,  &
       .6364,.4545,.7273,.3636,.4545,.4545,1.,1.0909,1.0909,.4545,       &
       .4545,.8182,.4545,.4545,1./
  stemp = '\\'
  schar(87:87) = stemp
  lnk=999
  ipk=0
  last=len(text)
  if(last.eq.0)  goto 999
  do i=1,ms
     in=index(text,search(i)(:ls(i)))
     if(in.ne.0.and.in.lt.lnk) then
        lnk=in
        ipk=i
     endif
  enddo
  !--- call normal text routine if no special characters
  if(ipk.eq.0)  then
     call gxstx(xt,yt,text)
     goto 999
  endif
  !--- start processing
  call jqtxfp(ierr,ifont,iprec)
  if(ierr.ne.0) goto 999
  ifdef=ifont
  call jqtxal(ierr,ihoru,ivertu)
  if(ierr.ne.0) goto 999
  call jqchh(ierr,chh)
  if(ierr.ne.0) goto 999
  call jqchxp(ierr,chxp)
  if(ierr.ne.0) goto 999
  call jqchup(ierr,chux,chuy)
  if(ierr.ne.0) goto 999
  sq=sqrt(chux**2+chuy**2)
  sang=-chux/sq
  cang=chuy/sq
  if(ihoru.eq.0) then
     ihor=1
  else
     ihor=ihoru
  endif
  if(ivertu.eq.0) then
     ivert=5
  else
     ivert=ivertu
  endif
  call jstxal(1,5)
  !--- total character spacing in y direction
  cyl=chf*chh
  !--- alignment factor from vertical alignment definition
  falign=.5*(ivert-5)
  !--- x and y displacement in case of superscript
  axup=-.5*cyl*sang
  ayup=.5*cyl*cang
  !--- set initial starting point
  x=xt
  y=yt
  do ipass=1,2
     !--- reset to defaults
     ip=ipk
     ln=lnk
     k=0
     do i=1,ms
        toggle(i)=.false.
     enddo
     xlift=0.
     ylift=0.
     ctf=1.
     kfpos=1
     !
     !--- loop over text
     !
30   if(k.ge.last) goto 70
     if(ip.eq.0)  then
        l=last
     else
        l=ln-1
     endif
     do i=1,l-k
        stemp=text(k+i:k+i)
        igts=index(schar,stemp)
        if(igts.eq.0) then
           cwf=1.
        else
           cwf=chgtsw(igts,kfpos)
        endif
        cxl=cwf*chxp*chh
        t=ctf*cxl
        ax=falign*axup
        ay=falign*ayup
        if(ipass.eq.2) then
           call gxstx(x+xlift,y+ylift,stemp)
        endif
        x=x+t*cang
        y=y+t*sang
     enddo
     k=l
     if(ip.ne.0)  then
        k=k+ls(ip)
        toggle(ip)=.not.toggle(ip)
        if(ip.eq.1)  then
           if(toggle(1)) then
              xlift=axup
              ylift=ayup
           else
              xlift=0.
              ylift=0.
           endif
        endif
        if(ip.le.2) then
           if(toggle(1).or.toggle(2))  then
              call jschh(crf*chh)
              ctf=crf
           else
              call jschh(chh)
              ctf=1.
           endif
        endif
        do i=3,ms
           if(i.eq.ip)  then
              if(toggle(i)) then
                 ifloc=is(i)
                 kfpos=i-1
                 if (i .eq. 3 .and. ifdef .ne. 1) ifloc = -12
              else
                 ifloc=ifdef
                 kfpos=1
              endif
              call jstxfp(ifloc,2)
           endif
        enddo
     endif
     ip=0
     ln=999
     if(k.lt.last)  then
        do i=1,ms
           in=index(text(k+1:last),search(i)(:ls(i)))
           if(in.ne.0.and.in.lt.ln)  then
              ln=in
              ip=i
           endif
        enddo
     endif
     if(ip.ne.0)  ln=ln+k
     goto 30
70   continue
     f=.5*(ihor-1)
     x=xt+f*(xt-x)+ax
     y=yt+f*(yt-y)+ay
  enddo
90 continue
  !--- restore defaults
  call jstxal(ihoru,ivertu)
  call jstxfp(ifont,iprec)
  call jschh(chh)
999 end subroutine gxtx
subroutine gxtx1(x,y,s,ac)
  implicit none
  real x,y
  !***********************************************************************
  !
  !   Purpose: plots a text if reference point inside active window.
  !
  !--- Input:
  !   x          x position
  !   y          y position
  !   s          text
  !   ac         active window
  !
  !   Author: H. Grote / CERN                        date: Dec. 9, 1988
  !                                           last mod: Dec. 9, 1988
  !
  !***********************************************************************
  real ac(4)
  character(*)  s
  if(x.ge.ac(1).and.x.le.ac(2).and.y.ge.ac(3).and.y.le.ac(4))       &
       call gxstx(x,y,s)
end subroutine gxtx1
subroutine gxundf
  use gxx11_common
  implicit none
  integer ifirst
  !***********************************************************************
  !
  !   Purpose: sets an integer for testing undefined variables
  !
  !   Author: H. Grote / CERN                        date: April 7, 1988
  !                                           last mod: April 7, 1988
  !
  !***********************************************************************

  save ifirst
  data ifirst/0/
  if(ifirst.eq.0)  then
     ifirst=1
     if(lundef.eq.654321)  then
        lundef=654320
     else
        lundef=654321
     endif
  endif
end subroutine gxundf
subroutine gxwait
  use gxx11_common
  implicit none
  !***********************************************************************
  !
  !   Purpose: waits for Input from keyboard  (e.g. <CR>) if interactive
  !   allows emergency stop when entering STOP
  !
  !   Author: H. Grote / CERN                        date: June 16, 1987
  !                                           last mod: Feb. 26, 1988
  !
  !***********************************************************************

  call gxundf
end subroutine gxwait
subroutine gxwclr(iun, ityp)
  use gxx11_common
  implicit none
  integer ityp,iun
  !***********************************************************************
  !
  !   Purpose: write .ps page end (clear)
  !
  !--- Input
  !   iun       not used
  !   ityp      not used
  !
  !   Author: H. Grote / CERN                          date: Apr. 27, 1995
  !                                                last mod: Apr. 27, 1995
  !
  !***********************************************************************

  write(iutlps, '(a)')  'gs showpage gr'
  ! flush in a portable way
  endfile(iutlps)
  backspace(iutlps)
  ipage = ipage + 1
  stortx = '%%Page: number'
  write(stortx(15:20), '(i6)')  ipage
  istotx = 20
end subroutine gxwclr
subroutine gxwpep(iun, ityp)
  use gxx11_common
  implicit none
  integer i,iday,ihour,iii,imonth,isec,ittp,ityp,iun,iyear,l,minute
  real fsc
  !***********************************************************************
  !
  !   Purpose: write .ps or .eps file prologue or epilogue
  !
  !--- Input
  !   iun       Output unit number, if = 0: epilogue, > 0: prologue
  !   ityp      type of Output: 113 = eps,
  !           else ps with 114 = portrait, 115 = landscape
  !
  !   Author: H. Grote / CERN                        date: Apr. 27, 1995
  !                                           last mod: Apr. 27, 1995
  !
  !***********************************************************************

  character(mlpro) pspro(mpspro), eppro(meppro), psdict(mdict),  &
       psepi(mpsep), epepi(mepep)
  character orient(2) * 12, head(mhead) * 60
  save ittp, orient, head, eppro, pspro, epepi, psepi, psdict
  data head / '$GX11psBegin',                                       &
       '0 setlinecap 0 setlinejoin', '<scale> ', ' ' /
  data orient / 'Portrait', 'Landscape' /
  data pspro / '%!PS-Adobe-2.0', '%%Title: ',                       &
       '%%Creator: gx11 version nnnn',                                   &
       '%%CreationDate: dd/mm/yy hh:mm',                                 &
       '%%Orientation: Landscape', '%%BoundingBox: nnn nnn nnn nnn',     &
       '%%Pages: (atend)', '%%EndComments' /
  data psepi / '%%Trailer', 'end', '%%EOF' /
  data eppro /'%!PS-Adobe-2.0 EPSF-2.0', '%%Title: filename',       &
       '%%Creator: gx11 version nnnn',                                   &
       '%%CreationDate: dd/mm/yy hh:mm',                                 &
       '%%Orientation: Landscape', '%%BoundingBox: nnn nnn nnn nnn',     &
       '%%Pages: 0', '%%EndComments' /
  data epepi / '$GX11psEnd showpage', '%%EOF' /
  data psdict / '/$GX11psDict 200 dict def', '$GX11psDict begin',   &
       '$GX11psDict /mtrx matrix put /l {lineto} bind def',              &
       '/m {moveto} bind def /s {stroke} bind def',                      &
       '/n {newpath} bind def /gs {gsave} bind def',                     &
       '/gr {grestore} bind def /clp {closepath} bind def',              &
       '/t {translate} bind def /sd {setdash} bind def',                 &
       '/fft {findfont} bind def /col-1 {} def /r {rotate} bind def',    &
       '/sf {scalefont setfont} bind def /sw {stringwidth} bind def',    &
       '/stwn { /fs exch def /fn exch def /text exch def fn findfont ',  &
       'fs sf text sw pop xs add /xs exch def} bind def',                &
       '/black {0 0 0 setrgbcolor} bind def',                            &
       '/blue {0 0 1 setrgbcolor} bind def',                             &
       '/green {0 1 0 setrgbcolor} bind def',                            &
       '/cyan {0 1 1 setrgbcolor} bind def',                             &
       '/red {1 0 0 setrgbcolor} bind def',                              &
       '/magenta {1 0 1 setrgbcolor} bind def',                          &
       '/yellow {1 1 0 setrgbcolor} bind def',                           &
       '/white {1 1 1 setrgbcolor} bind def',                            &
       '        end',                                                    &
       '/$GX11psBegin',                                                  &
       '     {$GX11psDict begin /$GX11psEnteredState save def} def',     &
       '/$GX11psEnd {$GX11psEnteredState restore end} def',              &
       '%%EndProlog' /

  if (iun .gt. 0)  then
     !--- prologue
     fsc = 1. / msfact
     do  i = 1, mhead
        pshead(i) = head(i)
     enddo
     !        call datime(idate, itime)
     !        call datesp(idate, iyear, imonth, iday)
     !        call timesp(itime, ihour, minute)
     call mydtime(iyear, imonth, iday, ihour, minute, isec)
     if (ityp .eq. 113)  then
        !--- eps
        eppro(2)(10:) = spsnam
        write(eppro(3)(25:29), '(f5.2)')  versio
        write(eppro(4)(17:24), '(i2.2,''/'',i2.2,''/'',i2.2)')        &
             iyear, imonth, iday
        write(eppro(4)(26:33), '(i2.2,'':'',i2.2)')                   &
             ihour, minute
        eppro(5)(16:24) = orient(iorips)
        write(eppro(6)(15:30), '(4i4)')  ibbox
     else
        pspro(2)(10:) = spsnam
        write(pspro(3)(25:29), '(f5.2)')  versio
        write(pspro(4)(17:24), '(i2.2,''/'',i2.2,''/'',i2.2)')        &
             iyear, imonth, iday
        write(pspro(4)(26:33), '(i2.2,'':'',i2.2)')                   &
             ihour, minute
        pspro(5)(16:24) = orient(iorips)
        write(pspro(6)(15:30), '(4i4)')  ibbox
     endif
     iutlps = iun
     ittp = ityp
     !--- write prologue
     if (ittp .eq. 113)  then
        do  i = 1, meppro
           call gxpnbl(eppro(i), iii, l)
           write(iutlps, '(a)')  eppro(i)(:l)
        enddo
     else
        do  i = 1, mpspro
           call gxpnbl(pspro(i), iii, l)
           write(iutlps, '(a)')  pspro(i)(:l)
        enddo
     endif
     do  i = 1, mdict
        call gxpnbl(psdict(i), iii, l)
        write(iutlps, '(a)')  psdict(i)(:l)
     enddo
     write(pshead(3)(:20), '(2f7.3, '' scale'')')  fsc, fsc
     do  i = 1, mhead
        call gxpnbl(pshead(i), iii, l)
        write(iutlps, '(a)')  pshead(i)(:l)
     enddo
     if (ittp .ne. 113)  then
        ipage = ipage + 1
        write(iutlps, '(a,i5)') '%%Page: number', ipage
     endif
  else
     if (ittp .eq. 113)  then
        !--- write epilogue
        do  i = 1, mepep
           call gxpnbl(epepi(i), iii, l)
           write(iutlps, '(a)')  epepi(i)(:l)
        enddo
     else
        if(iclear .ne. 0) then
           call gclrwk(0, 1)
           iclear = 0
        endif
        do  i = 1, mpsep
           call gxpnbl(psepi(i), iii, l)
           write(iutlps, '(a)')  psepi(i)(:l)
        enddo
     endif
  endif
end subroutine gxwpep
subroutine gxwpl(np, xp1, yp1, ifill)
  use gxx11_common
  implicit none
  integer i,i1,i2,ierr,ifill,iforl,iii,k,kadd,l,np
  real f1,f2,r,v1,v2
  !***********************************************************************
  !
  !   Purpose: write polyline into .ps or .eps file
  !
  !--- Input
  !   np      number of points
  !   xp1      x coordinates
  !   yp1      y coordinates
  !   ifill   fill area request: 0 = no, 1 = yes
  !
  !   Author: H. Grote / CERN                        date: Apr. 27, 1995
  !                                              last mod: May  23, 1995
  !
  !***********************************************************************

  real xp1(*), yp1(*)
  character eloc * 40, sline * (mline), style(4) * 24,              &
       formt1 * 24, formt2 * 24
  data style / '[] 0 sd', '[20] 0 sd', '[2 10] 0 sd',               &
       '[20 10 2 10] 0 sd' /
  data eloc / 'col-1 s' /
  data formt1 /'('' n'', 2ix, '' m'')'/
  data formt2 /'(2ix, '' l'')'/

  if (istotx .gt. 0)  then
     write(iutlps, '(a)')  stortx(:istotx)
     istotx = 0
  endif
  if (iorips .eq. 1)  then
     v1 = xp1(1) - rx11pr(1)
     v2 = yp1(1) - rx11pr(3)
     f1 = msfact * (ibbox(3) - ibbox(1)) / (rx11pr(2) - rx11pr(1))
     f2 = msfact * (ibbox(4) - ibbox(2)) / (rx11pr(4) - rx11pr(3))
  else
     v1 = xp1(1) - rx11pr(1)
     v2 = yp1(1) - rx11pr(3)
     f1 = msfact * (ibbox(4) - ibbox(2)) / (rx11pr(2) - rx11pr(1))
     f2 = msfact * (ibbox(3) - ibbox(1)) / (rx11pr(4) - rx11pr(3))
  endif
  sline = colour(icucol)
  call gxpnbl(sline, iii, k)
  if (iorips .eq. 1)  then
     i1 = f1 * v1 + msfact * ibbox(2) + 0.5
     i2 = f2 * v2 + msfact * ibbox(1) + 0.5
  else
     i2 = f1 * v1 + msfact * ibbox(2) + 0.5
     i1 = f2 * v2 + msfact * ibbox(1) + 0.5
     i1 = msfact * ibbox(3) - i1
  endif
  call jqlwsc(ierr, r)
  call jqln(ierr, i)
  i = max(1, min(i,4))
  write(iutlps, '(f6.3, a, a)')  0.75 * r,                          &
       ' setlinewidth ', style(i)
  iforl =                                                           &
       max(log10(float(abs(i1)+1))+3., log10(float(abs(i2)+1))+3.)
  write(formt1(10:10), '(i1)') iforl
  kadd = 2 * (iforl + 2)
  write(sline(k+1:k+kadd), formt1)  i1, i2
  k = k + kadd
  do  i = 2, np
     v1 = xp1(i) - rx11pr(1)
     v2 = yp1(i) - rx11pr(3)

     if (k + 16 .gt. mline)  then
        call gxpnbl(sline, iii, l)
        write(iutlps, '(a)')  sline(:l)
        k = 0
        sline = ' '
     endif
     if (iorips .eq. 1)  then
        i1 = f1 * v1 + msfact * ibbox(2) + 0.5
        i2 = f2 * v2 + msfact * ibbox(1) + 0.5
     else
        i2 = f1 * v1 + msfact * ibbox(2) + 0.5
        i1 = f2 * v2 + msfact * ibbox(1) + 0.5
        i1 = msfact * ibbox(3) - i1
     endif
     iforl =                                                         &
          max(log10(float(abs(i1)+1))+3., log10(float(abs(i2)+1))+3.)
     write(formt2(4:4), '(i1)') iforl
     kadd = 2 * (iforl + 1)
     write(sline(k+1:k+kadd), formt2)  i1, i2
     k = k + kadd
  enddo
  if (k .gt. 0) then
     call gxpnbl(sline, iii, l)
     write(iutlps, '(a)')  sline(:l)
  endif
  if (ifill .eq. 0)  then
     write(iutlps, '(a)')  eloc
  else
     write(iutlps, '(a)')  'fill ' // eloc
  endif
end subroutine gxwpl
subroutine gxwpm(np, xp1, yp1)
  use gxx11_common
  implicit none
  integer i,ifill,ip,is,j,n,np
  real sq21,xf,xmf,xms,xmsm,xmsq,xmsqm,xpsf,xs,xsm,xsq,xsqm,yf,ymf, &
       yms,ymsm,ymsq,ymsqm,ys,ysm,ysq,ysqm
  !***********************************************************************
  !
  !   Purpose: plot marker symbol on display and/or PostScript output
  !
  !--- Input
  !   np      number of marker symbols
  !   xp1      x coordinates
  !   yp1      y coordinates
  !
  !   Author: H. Grote / CERN                        date: July 6, 1995
  !                                              last mod: July 6, 1995
  !
  !***********************************************************************

  parameter (xs = 0.005, ys = 0.005, sq21 = 0.4142135,              &
       xms = -xs, yms = -ys,                                             &
       xsm = 0.1 * xs, ysm = 0.1 * ys,                                   &
       xmsm = -xsm, ymsm = -ysm,                                         &
       xsq = xs * sq21, ysq = ys * sq21,                                 &
       xmsq = -xsq, ymsq = -ysq,                                         &
       xsqm = xsm * sq21, ysqm = ysm * sq21,                             &
       xmsqm = -xsqm, ymsqm = -ysqm )
  real xp1(*), yp1(*)
  real xsym(9,4,5), ysym(9,4,5), xlps(9), yloc(9), xlwd(9)
  integer nsym(4,5)
  character(1) dum
  data xsym /                                                       &
       xmsm, xmsm, xmsqm, xsqm, xsm, xsm, xsqm, xmsqm, xmsm, 27 * 0.,    &
       xms, xs, 7*0., 0., 0., 7*0., 18 * 0.,                             &
       xms, xs, 7*0., 0., 0., 7*0., xms, xs, 7*0., xms, xs, 7*0.,        &
       xms, xms, xmsq, xsq, xs, xs, xsq, xmsq, xms, 27 * 0.,             &
       xms, xs, 7*0., xms, xs, 7*0., 18 * 0. /
  data ysym /                                                       &
       ymsqm, ysqm, ysm, ysm, ysqm, ymsqm, ymsm, ymsm, ymsqm, 27 * 0.,   &
       0., 0., 7*0., ys, yms, 7*0., 18 * 0.,                             &
       0., 0., 7*0., ys, yms, 7*0., yms, ys, 7*0., ys, yms, 7*0.,        &
       ymsq, ysq, ys, ys, ysq, ymsq, yms, yms, ymsq, 27 * 0.,            &
       yms, ys, 7*0., ys, yms, 7*0, 18 * 0. /
  data nsym / 9, 0, 0, 0,                                           &
       2, 2, 0, 0,                                                       &
       2, 2, 2, 2,                                                       &
       9, 0, 0, 0,                                                       &
       2, 2, 0, 0 /

  call gxqvar('XMETAF', i, xmf, dum)
  call gxqvar('YMETAF', i, ymf, dum)
  if (xmf .ne. 0.)  then
     xpsf = ymf / xmf
  else
     xpsf = 1.
  endif
  xf = rx11pr(8) * (rx11pr(2) - rx11pr(1))
  yf = rx11pr(8) * (rx11pr(4) - rx11pr(3))
  is = mod(ix11pr(4) - 1, 5) + 1
  if (is .eq. 1)  then
     ifill = 1
  else
     ifill = 0
  endif
  do  ip = 1, np
     do  i = 1, 4
        n = nsym(i,is)
        if (n .gt. 0)  then
           do  j = 1, n
              xlwd(j) = xp1(ip) + xf * xsym(j,i,is)
              yloc(j) = yp1(ip) + yf * ysym(j,i,is)
              xlps(j) = xp1(ip) + xpsf * xf * xsym(j,i,is)
           enddo
           if (ipseps .ne. 0)  call gxwpl(n, xlps, yloc, ifill)
        endif
     enddo
  enddo
end subroutine gxwpm
subroutine gxwtx(xp1, yp1, txin)
  use gxx11_common
  implicit none
  integer i,i1,i11,i2,i22,ia,iang1,iangle,ie,ifont,ifos,ihl,iii,    &
       iprec,ivl,k,lf,lint,lt,lt1,lt2,lt3,mltx,mltx2
  real chh,f1,f2,v1,v2,x,xp1,xup,y,yp1,yup
  !***********************************************************************
  !
  !   Purpose: write text with predefined font .ps or .eps file
  !
  !--- Input
  !   xp1      x coordinate
  !   yp1      y coordinate
  !   txin    text
  !
  !   Author: H. Grote / CERN                        date: May 10, 1995
  !                                              last mod: May 10, 1995
  !
  !***********************************************************************

  parameter (mltx = 120, mltx2 = 2 * mltx)
  character txin * (*)
  character sline * (mline), txlc * (mltx2)
  character(24) tloc1, tloc2, tloc3, sfont(mtfont)
  integer ifosiz(mtfont)
  data sfont / '/Times-Italic', '/Times-Bold', '/Times-BoldItalic', &
       '/Helvetica', '/Helvetica-Oblique', '/Helvetica-Bold',            &
       '/Helvetica-BoldOblique', '/Courier', '/Courier-Oblique',         &
       '/Courier-Bold',  '/Courier-BoldOblique', '/Symbol' /
  data ifosiz / 1030, 1000, 1025, 930, 930,                         &
       930, 930, 1205, 1205, 1170, 1165, 1005 /

  if (istotx .gt. 0)  then
     write(iutlps, '(a)')  stortx(:istotx)
     istotx = 0
  endif
  !--- copy text to local, treat ()
  txlc = ' '
  call gxpnbl(txin, iii, lint)
  lt = 0
  do  i = 1, min(mltx, lint)
     if(txin(i:i) .eq. '(' .or. txin(i:i) .eq. ')')  then
        lt = lt + 1
        txlc(lt:lt) = '\\'
     endif
     lt = lt + 1
     txlc(lt:lt) = txin(i:i)
  enddo
  call jqtxfp(ie, ifont, iprec)
  ifont = max(1, min(mtfont, abs(ifont)))
  call jqtxal(ie, ihl, ivl)
  call jqchh(ie, chh)
  call jqchup(ie, xup, yup)
  iang1 = atan2(-xup, yup) * 45. / atan(1.) + 0.5
  iangle = iang1 + (iorips - 1) * 90 + 0.5
  ifos = msfact * ifosiz(ifont) * chh / (rx11pr(4) - rx11pr(3))     &
       + 0.5
  if (ihl .eq. 2)  then
     tloc1 = ' xs 2 div'
     lt1 = 9
  elseif ( ihl .eq. 3)  then
     tloc1 = ' xs'
     lt1 = 3
  endif
  tloc2 = ' neg 0'
  lt2 = 6
  tloc3 = ' t 0 0 m'
  lt3 = 8
  call gxpnbl(sfont(ifont), iii, lf)
  x = xp1
  y = yp1 - 0.25 * (5 - ivl) * chh
  if (iorips .eq. 1)  then
     v1 = x - rx11pr(1)
     v2 = y - rx11pr(3)
     f1 = msfact * (ibbox(3) - ibbox(1)) / (rx11pr(2) - rx11pr(1))
     f2 = msfact * (ibbox(4) - ibbox(2)) / (rx11pr(4) - rx11pr(3))
  else
     v1 = x - rx11pr(1)
     v2 = y - rx11pr(3)
     f1 = msfact * (ibbox(4) - ibbox(2)) / (rx11pr(2) - rx11pr(1))
     f2 = msfact * (ibbox(3) - ibbox(1)) / (rx11pr(4) - rx11pr(3))
  endif
  !--- horizontal alignment - uses font width from font def.
  if (ihl .gt. 1)  then
     write(iutlps, '(a)')  '/xs 0 def'
     write(iutlps, '(a, a, a)')  '(', txlc(:lt), ')'
     write(iutlps, '(a, i6, a)') sfont(ifont)(:lf), ifos, ' stwn'
  endif
  sline = colour(icucol)
  call gxpnbl(sline, iii, k)
  if (iorips .eq. 1)  then
     i1 = f1 * v1 + msfact * ibbox(2) + 0.5
     i2 = f2 * v2 + msfact * ibbox(1) + 0.5
  else
     i2 = f1 * v1 + msfact * ibbox(2) + 0.5
     i1 = f2 * v2 + msfact * ibbox(1) + 0.5
     i1 = msfact * ibbox(3) - i1
  endif
  write(sline(k+1:k+10), '(2i5)')  i1, i2
  k = k + 10
  sline(k+1:k+3) = ' t '
  k = k + 3
  write(sline(k+1:k+5), '(i5)')  iangle
  k = k + 5
  sline(k+1:k+2) = ' r'
  k = k + 2
  if (ihl .gt. 1)  then
     sline(k+1:) =                                                   &
          tloc1(:lt1) // tloc2(:lt2) // tloc3(:lt3)
     k = k + lt1 + lt2 + lt3
  endif
  sline(k+1:) = ' 0 0 m'
  write(iutlps, '(a)')  sline(:k)
  write(iutlps, '(a, a, i6, a)')  sfont(ifont)(:lf), ' fft ',       &
       ifos, ' sf 0 0 m'
  write(iutlps, '(a, a, a)')  '(', txlc(:lt), ')'
  sline = 'show'
  k = 4
  tloc2 = ' 0'
  lt2 = 2
  if (ihl .gt. 1)  then
     sline(k+1:) =                                                   &
          tloc1(:lt1) // tloc2(:lt2) // tloc3(:lt3)
     k = k + lt1 + lt2 + lt3
  endif
  ia = -iangle
  write(sline(k+1:k+6), '(i6)')  ia
  k = k + 6
  sline(k+1:k+2) = ' r'
  k = k + 2
  i11 = -i1
  i22 = -i2
  write(sline(k+1:k+12), '(2i6)')  i11, i22
  k = k + 12
  sline(k+1:k+2) = ' t'
  k = k + 2
  write(iutlps, '(a)')  sline(:k)
end subroutine gxwtx
subroutine jqmk(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1
  ierr = 0
  i1 = ivals(1)
end subroutine jqmk
subroutine jqfais(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1
  ierr = 0
  i1 = ivals(2)
end subroutine jqfais
subroutine jqtxal(ierr, i1, i2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr, i1,i2
  ierr = 0
  i1 = ivals(3)
  i2 = ivals(4)
end subroutine jqtxal
subroutine jqtxfp(ierr, i1, i2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr, i1,i2
  ierr = 0
  i1 = ivals(5)
  i2 = ivals(6)
end subroutine jqtxfp
subroutine jqpmci(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr, i1
  ierr = 0
  i1 = ivals(7)
end subroutine jqpmci
subroutine jswks(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(8) = i1
end subroutine jswks
subroutine jqwks(i1, ierr, i2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr, i1,i2
  ierr = 0
  i2 = ivals(8)
end subroutine jqwks
subroutine jqchup(ierr, r1, r2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr
  real r1,r2
  ierr = 0
  r1 = rvals(1)
  r2 = rvals(2)
end subroutine jqchup
subroutine jqchh(ierr, r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr
  real r1
  ierr = 0
  r1 = rvals(3)
end subroutine jqchh
subroutine jqtxci(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1
  ierr = 0
  i1 = ivals(9)
end subroutine jqtxci
subroutine jqnt(i1, ierr, ar1, ar2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1,i
  real ar1(4), ar2(4)
  ierr = 0
  do  i = 1, 4
     ar1(i) = rvals(i+3)
     ar2(i) = rvals(i+7)
  enddo
end subroutine jqnt
subroutine jqmksc(ierr, r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr
  real r1
  ierr = 0
  r1 = rvals(14)
end subroutine jqmksc
subroutine jsvp(i1, r1, r2, r3, r4)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  real r1,r2,r3,r4
  rvals(8) = r1
  rvals(9) = r2
  rvals(10) = r3
  rvals(11) = r4
end subroutine jsvp
subroutine jqplci(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1
  ierr = 0
  i1 = ivals(11)
end subroutine jqplci
subroutine jschxp(r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  real r1
  rvals(12) = r1
end subroutine jschxp
subroutine jqchxp(ierr, r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr
  real r1
  ierr = 0
  r1 = rvals(12)
end subroutine jqchxp
subroutine jqlwsc(ierr, r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr
  real r1
  ierr = 0
  r1 = rvals(13)
end subroutine jqlwsc
subroutine jqln(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1
  ierr = 0
  i1 = ivals(12)
end subroutine jqln
subroutine jqcntn(ierr, i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer ierr,i1
  ierr = 0
  i1 = ivals(13)
end subroutine jqcntn
subroutine gclrwk(i1, i2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1,i2
  if (ipseps .eq. 1)  call gxwclr(i1, i2)
end subroutine gclrwk
subroutine gtx(r1, r2, string)
  use gxx11_common
  use gxx11_aux
  implicit none
  real r1,r2
  character(*) string
  if (ipseps .ne. 0)  call gxwtx(r1, r2, string)
end subroutine gtx
subroutine gfa(i1, ar1, ar2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  real ar1(*),ar2(*)
  if (ipseps .ne. 0)  call gxwpl(i1, ar1, ar2, ivals(2))
end subroutine gfa
subroutine gpl(i1, ar1, ar2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  real ar1(*),ar2(*)
  if (ipseps .ne. 0)  call gxwpl(i1, ar1, ar2, 0)
end subroutine gpl
subroutine jsmk(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(1) = i1
  ix11pr(4) = i1
end subroutine jsmk
subroutine jsfais(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(2) = i1
end subroutine jsfais
subroutine jstxal(i1, i2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1,i2
  ivals(3) = i1
  ivals(4) = i2
  ix11pr(1) = i1
  ix11pr(2) = i2
end subroutine jstxal
subroutine jstxfp(i1, i2)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1,i2
  ivals(5) = i1
  ivals(6) = i2
end subroutine jstxfp
subroutine jspmci(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(7) = i1
  ivals(14) = 0
end subroutine jspmci
subroutine jschup(r1, r2)
  use gxx11_common
  use gxx11_aux
  implicit none
  real r1,r2
  rvals(1) = r1
  rvals(2) = r2
  rx11pr(5) = r1
  rx11pr(6) = r2
end subroutine jschup
subroutine jschh(r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  real r1
  rvals(3) = r1
  rx11pr(7) = r1
end subroutine jschh
subroutine jstxci(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(9) = i1
  ivals(14) = 0
end subroutine jstxci
subroutine jsmksc(r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  real r1
  rvals(14) = r1
  rx11pr(8) = r1
end subroutine jsmksc
subroutine jswn(i1, r1, r2, r3, r4)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  real r1,r2,r3,r4
  rvals(4) = r1
  rvals(5) = r2
  rvals(6) = r3
  rvals(7) = r4
  rx11pr(1) = r1
  rx11pr(2) = r2
  rx11pr(3) = r3
  rx11pr(4) = r4
  if (r2 .gt. r1)  then
     fxpix = nxpix / (r2 - r1)
  else
     fxpix = 1.
  endif
  if (r4 .gt. r3)  then
     fypix = nypix / (r4 - r3)
  else
     fypix = 1.
  endif
end subroutine jswn
subroutine jsplci(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(11) = i1
  ivals(14) = 0
end subroutine jsplci
subroutine jslwsc(r1)
  use gxx11_common
  use gxx11_aux
  implicit none
  real r1
  rvals(13) = r1
end subroutine jslwsc
subroutine jsln(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1,iz
  ivals(12) = i1
  iz = i1 - 1
end subroutine jsln
subroutine jslctp(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  ivals(14) = i1
end subroutine jslctp
subroutine jqlctp(i1)
  use gxx11_common
  use gxx11_aux
  implicit none
  integer i1
  i1 = ivals(14)
end subroutine jqlctp
subroutine wacwk(iw)
  use gxx11_common
  implicit none
  integer iw
  !***********************************************************************
  !
  !   Purpose: Activate workstation
  !
  !--- Input
  !   iw       workstation number
  !   Author: H. Grote / CERN                        date: Jan. 25, 1994
  !                                           last mod: Jan. 25, 1994
  !***********************************************************************

  if (iw .gt. 0 .and. iw .le. mx11tf)  then
     ix11tf(iw) = 1
  endif
end subroutine wacwk
subroutine wclks
  implicit none
  !***********************************************************************
  !
  !   Purpose: Close X11 package
  !
  !   Author: H. Grote / CERN                        date: Jan. 25, 1994
  !                                           last mod: Jan. 25, 1994
  !***********************************************************************
end subroutine wclks
subroutine wclwk(iw)
  use gxx11_common
  implicit none
  integer iw
  !***********************************************************************
  !
  !   Purpose: Close workstation
  !
  !--- Input
  !   iw       workstation number
  !   Author: H. Grote / CERN                        date: Jan. 25, 1994
  !                                           last mod: Jan. 25, 1994
  !***********************************************************************

end subroutine wclwk
subroutine wdawk(iw)
  use gxx11_common
  implicit none
  integer iw
  !***********************************************************************
  !
  !   Purpose: Deactivate workstation
  !
  !--- Input
  !   iw       workstation number
  !   Author: H. Grote / CERN                        date: Jan. 25, 1994
  !                                           last mod: Jan. 25, 1994
  !***********************************************************************

  if (iw .gt. 0 .and. iw .le. mx11tf)  then
     ix11tf(iw) = 0
  endif
end subroutine wdawk
subroutine wopks(ieu, idum)
  use gxx11_common
  implicit none
  integer i,idum,ieu
  !***********************************************************************
  !
  !   Purpose: Open X11 package
  !
  !--- Input
  !   ieu      error file unit (not used)
  !   idum     dummy
  !   Author: H. Grote / CERN                        date: Jan. 25, 1994
  !                                           last mod: Jan. 25, 1994
  !***********************************************************************

  !--- preset workstations to inactive and closed
  do  i = 1, mx11tf
     ix11tf(i) = 0
     ix11op(i) = 0
  enddo
end subroutine wopks
logical function affirm(sus)
  implicit none
  character(1) sus
  affirm = sus.eq.'y'.or.sus.eq.'Y'.or.sus.eq.'o'.or.sus.eq.'O'
end function affirm
subroutine wopwk(iw, icont, it)
  use gxx11_common
  implicit none
  integer icont,it,iw
  !      integer ix,iy
  !      real r
  !***********************************************************************
  !
  !   Purpose: Open workstation
  !
  !--- Input
  !   iw       workstation number
  !   icont    dummy
  !   it       dummy
  !   Author: H. Grote / CERN                        date: Jan. 25, 1994
  !                                           last mod: Jan. 25, 1994
  !***********************************************************************

end subroutine wopwk
