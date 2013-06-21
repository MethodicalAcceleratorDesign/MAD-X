    subroutine setup(resp,a,im,ic,nm,nc)
! ****************************************************
!                                                    *
!    Set DOUBLE PRECISION (A)                        *
!    response matrix in FORTRAN storage format       *
!    RESP comes from "C" routine                     *
!                                                    *
!     Author: WFH  05.02.02                          *
!                                                    *
! ****************************************************
      implicit none

      integer im,ic,nm,nc    ! here we are
      double precision resp,a(nm, nc)

!     write(*,*) 'in setup: ',resp, im+1, ic+1
      a(im+1,ic+1) = resp

      return
    end subroutine setup

    subroutine setupi(resp,a,im,ic,nm,nc)
! ****************************************************
!                                                    *
!    Set INTEGER (A)                                 *
!    matrix in FORTRAN storage format                *
!    RESP comes from "C" routine                     *
!                                                    *
!     Author: WFH  05.02.02                          *
!                                                    *
! ****************************************************
      implicit none

      integer im,ic,nm,nc
      integer resp,a(nm, nc)

!     write(*,*) 'in setupi: ',resp, im+1, ic+1
      a(im+1,ic+1) = resp
!     write(*,*) 'done in setupi'

      return
    end subroutine setupi

    subroutine micit(a,conm,xin,cin,res,nx,rms,im,ic,iter,ny,ax,cinx, &
     &xinx,resx,rho,ptop,rmss,xrms,xptp,xiter,ifail)
! ****************************************************
!                                                    *
!    Driving routine for MICADO correction           *
!                                                    *
!     Author: WFH  05.02.02                          *
!                                                    *
! ****************************************************

      implicit none

      integer im,ic,iter,i,j,nx(ic),ny(ic)
      real rms,ax(im,ic),cinx(ic),xinx(im),resx(im),rho(3*ic),ptop(ic), &
     &rmss(ic),xrms(ic),xptp(ic),xiter(ic),rzero
      parameter(rzero=0e0)
      double precision a(im,ic),xin(im),cin(ic),res(im)
      character*16 conm(ic)
      integer      n
      integer      ifail

      do  j = 1,ic
         call f_ctof(n, conm(j), 16)
      enddo

      do  i = 1,im
         do  j = 1,ic
            ax(i,j) = a(i,j)
!           write(*,*) i,j,ax(i,j),a(i,j)
!           ny(j) = j-1
            ny(j) = j
            cinx(j) = rzero
         enddo
      enddo
      
      do  i = 1,im
         xinx(i) = xin(i)
         resx(i) = rzero
      enddo

      write(*,*) ' '
      write(*,*) 'start MICADO correction with ',iter,' correctors'
      write(*,*) ' '

      call micado(ax,conm,xinx,resx,cinx,ny,rms,im,ic,iter,rho,ptop,    &
     &rmss,xrms,xptp,xiter,ifail)

      do  i = 1,ic
         cin(i) = cinx(i)
         nx(ny(i)) = i
      enddo

      do  i = 1,im
         res(i) = resx(i)
      enddo

      return
    end subroutine micit

    subroutine haveit(a,xin,cin,res,nx,im,ic,cb,xmeas,xres,y,z,xd)
! ****************************************************
!                                                    *
!    Driving routine for LSQ correction              *
!                                                    *
!     Author: WFH  05.02.02                          *
!                                                    *
! ****************************************************
      implicit none

      integer im,ic,nx(ic),i,j
      double precision a(im,ic),xin(im),cin(ic),res(im),cb(ic),         &
     &xmeas(im),xres(im),y(ic,im),z(ic,ic),xd(ic),zero
      parameter(zero=0d0)

!      do  i = 1,im
!        do  j = 1,ic
!          write(*,*) i,j,a(i,j)
!        enddo
!      enddo
      
      do  i = 1,im
         res(i) = zero
      enddo

!     write(*,*) '==> ',res
      write(*,*) ' '
      write(*,*) 'start LEAST SQUARES correction with all correctors'
      write(*,*) ' '

      call solsql(im,ic,a,xin,res,cin,cb,xmeas,xres,y,z,xd)

! 6001 format(1X,'Corrector: ',I4,'   strength: ',F12.8)

      write(*,*) ' '
      write(*,*) 'end LEAST SQUARES correction with all correctors'
      write(*,*) ' '

      do  i = 1,ic
!         write(*,*) i,cin(i)
          nx(i) = i
      enddo

      return
    end subroutine haveit
      
    subroutine svddec_m(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,    &
         &ws,wvec,sortw,                                              &
         &sngcut, sngval,                                             &
         &im,ic,iflag,sing,dbg)
! ****************************************************
!                                                    *
!    Performs SVD and analysis for matrix with more  *
!    monitors than correctors.                       *
!                                                    *
!     Author: WFH  12.09.02                          *
!                                                    *
! ****************************************************
      implicit none
      integer im,ic,i,j,jj,ii
      integer iflag,sing(2,ic)
      double precision a(im,ic)
      double precision svdmat(im,ic)
      double precision umat(im,ic), utmat(ic,im)
      double precision vmat(im,ic), vtmat(ic,im)
      double precision wmat(im,ic), wtmat(ic,im)
      double precision wvec(ic)
      double precision rat, zero, sngval, sngcut
      double precision ws(ic)
      integer amater, svdmx, svdnx, svdnm
      integer sortw(ic)
      integer nsing
      logical matu, matv
      integer dbg
      parameter(zero = 0d0)
      parameter(nsing = 5)

      if(dbg.eq.2) then
         write(*,*) 'SVD parameters: '
         write(*,*) 'SNGCUT:         ',sngcut
         write(*,*) 'SNGVAL:         ',sngval
      endif

      MATU = .TRUE.
      MATV = .TRUE.
      iflag = 0

      SVDNM = max(ic,im)
      svdmx = im
      svdnx = ic

      do  i = 1,im
        do  j = 1,ic
           svdmat(i,j) = a(i,j)
        enddo
      enddo

      if(dbg.eq.2) then
          write(*,*) 'A0:'
          do  j = 1,im
             write(*,6003) (svdmat(j,i),i=1,ic)
          enddo
      endif

      call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,                 &
     &matv,vmat,amater,ws)

 6001 format(1X,'Corrector: ',I4,'   sing: ',F12.4)
 6002 format('VMAT: ',I4,I4,5X,F12.6,2X,F12.6)
 6003 format(16(2X,F7.2))
 6004 format(16(2X,F7.2))

      if(amater.ne.0) then
        write(*,*) 'end SVD with error code: ',amater
      endif

      if(dbg.eq.1) then
         do  i = 1,ic
!            write(*,*) i,wvec(i)
            write(*,6001) i,wvec(i)
            wmat(i,i) = wvec(i)
         enddo
      endif

      call rvord(wvec,sortw,ws,ic)

      if(dbg.eq.1) then
         do  i = 1,ic
            write(*,*) i,sortw(i),wvec(sortw(i))
         enddo
      endif

      do  ii = 1,min(nsing,ic)
         i = sortw(ii)
 
         if(dbg.eq.1) then
            write(61,*) wvec(i)
         endif

         if(abs(wvec(i)).lt.sngval) then
            if(dbg.eq.1) then
               do  j = 1,ic
                  write(*,6002) i,j,vmat(j,i)
               enddo
            endif
            do  j = 1,ic-1
               do  jj = j+1,ic
                  if(abs(vmat(j,i)).gt.1.0E-4) then
                     rat = abs(vmat(j,i)) + abs(vmat(jj,i))
!                    rat = abs(vmat(j,i) - vmat(jj,i))
!                    rat = abs(vmat(j,i) + vmat(jj,i))
                     rat = rat/abs(abs(vmat(j,i)) - abs(vmat(jj,i)))
                     if(rat.gt.sngcut) then
                        if(dbg.eq.1) then
                           write(*,*) 'dependent pair: ',i,j,jj,rat
                           write(65,*) 'dependent pair: ',i,j,jj,rat
                        endif
                        if(iflag.lt.(ic*ic*ic)) then
                           iflag = iflag + 1
                           sing(1,iflag) =  j - 1
                           sing(2,iflag) = jj - 1
                        endif
                     endif
                  endif
               enddo
            enddo
         endif
      enddo

      if(dbg.eq.1) then
         do  j=1,iflag
            write(66,*) j,sing(1,j),sing(2,j)
         enddo
      endif
      
      return
    end subroutine svddec_m

    subroutine svddec_c(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,    &
         &ws,wvec,sortw,                                              &
         &sngcut, sngval,                                             &
         &im,ic,iflag,sing,dbg)
! ****************************************************
!                                                    *
!    Performs SVD and analysis for matrix with more  *
!    correctors than monitors.                       *
!                                                    *
!     Author: WFH  12.09.02                          *
!                                                    *
! ****************************************************
      implicit none
      integer im,ic,i,j,jj,ii
      integer iflag,sing(2,ic)
      double precision a(im,ic)
      double precision svdmat(ic,im)
      double precision umat(ic,im), utmat(im,ic)
      double precision vmat(ic,im), vtmat(im,ic)
      double precision wmat(ic,im), wtmat(im,ic)
      double precision wvec(ic)
      double precision rat, zero, sngcut, sngval
      double precision ws(ic)
      integer sortw(ic)
      integer amater, svdmx, svdnx, svdnm
      integer nsing
      logical matu, matv
      integer dbg
      parameter(zero = 0d0)
      parameter(nsing = 5)

      if(dbg.eq.2) then
        write(*,*) 'SVD parameters: '
        write(*,*) 'SNGCUT:         ',sngcut
        write(*,*) 'SNGVAL:         ',sngval
      endif

      MATU = .TRUE.
      MATV = .TRUE.
      iflag = 0

      SVDNM = max(ic,im)
      svdmx = ic
      svdnx = im

      do  i = 1,im
        do  j = 1,ic
            svdmat(j,i) = a(i,j)
!           svdmat(i,j) = a(i,j)
        enddo
      enddo

 !8373 continue

      if(dbg.eq.2) then
         write(*,*) 'A0:'
         do  j = 1,ic
            write(*,6003) (svdmat(j,i),i=1,im)
         enddo
      endif

      call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,                 &
     &matv,vmat,amater,ws)

 6001 format(1X,'Corrector: ',I4,'   sing: ',F12.4)
 6002 format('UMAT: ',I4,I4,5X,F12.6,2X,F12.6)
 6003 format(16(2X,F7.2))
 6004 format(16(2X,F7.2))

      if(amater.ne.0) then
        write(*,*) 'end SVD with error code: ',amater
      endif

      if(dbg.eq.1) then
         do  i = 1,im
!            write(*,*) i,wvec(I)
            write(*,6001) i,wvec(I)
            wmat(i,i) = wvec(i)
         enddo
      endif

      call rvord(wvec,sortw,ws,im)

      if(dbg.eq.1) then
         do  i = 1,im
            write(*,*) i,sortw(i),wvec(sortw(i))
         enddo
      endif

      do  ii = 1,min(nsing,im)
         i = sortw(ii)

         if(dbg.eq.1) then
            write(61,*) wvec(i)
         endif

         if(abs(wvec(i)).lt.sngval) then
            if(dbg.eq.1) then
               do  j = 1,ic
                  write(*,6002) i,j,umat(j,i)
               enddo
            endif
            do  j = 1,ic-1
               do  jj = j+1,ic
                  if(abs(umat(j,i)).gt.1.0E-4) then
!                    rat = abs(umat(j,i) - umat(jj,i))
                     rat = abs(umat(j,i)) +  abs(umat(jj,i))
                     rat = rat/abs(abs(umat(j,i)) - abs(umat(jj,i)))

                   if(dbg.eq.1) then
                     write(62,*) wvec(i),rat
                   endif

                   if(rat.gt.sngcut) then
                      if(dbg.eq.1) then
                         write(*,*) 'dependent pair: ',j,jj,rat
                         write(65,*) 'dependent pair: ',j,jj,rat
                      endif

! Ghislain : wanring compare this line with equivalent  if(iflag.lt.(ic*ic*ic)) then
                      if(iflag.lt.ic) then
                         iflag = iflag + 1
                         sing(1,iflag) =  j - 1
                         sing(2,iflag) = jj - 1
                      endif
                   endif
                endif
             enddo
          enddo
       endif
    enddo

    return
  end subroutine svddec_c

  subroutine svdcorr_m(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,   &
       &xin,xc,xout,xa,xb,xpred,ws,wvec,                            &
       &sortw,nx,im,ic,iflag,dbg)
! ******************************************************
!                                                      *
!    Performs SVD and correction for matrix with more  *
!    monitors than correctors.                         *
!                                                      *
!     Author: WFH  12.09.02                            *
!                                                      *
! ******************************************************
      implicit none
      integer im,ic,nx(ic),i,j
      integer iflag
      double precision a(im,ic)
      double precision svdmat(im,ic)
      double precision umat(im,ic), utmat(ic,im)
      double precision vmat(im,ic), vtmat(ic,im)
      double precision wmat(im,ic), wtmat(ic,im)
      double precision xin(im),xout(im),xc(ic)
      double precision xa(im),xb(im),xpred(im),ws(ic),wvec(ic)
      integer amater, svdmx, svdnx, svdnm
      integer sortw(ic)
      logical matu, matv
      integer dbg

      MATU = .TRUE.
      MATV = .TRUE.
      iflag = 0
      write(*,*) ' '
      write(*,*) 'start SVD correction using ',ic,' correctors'
      write(*,*) ' '

      SVDNM = max(ic,im)
      svdmx = im
      svdnx = ic

      do  i = 1,ic
         nx(i) = i
      enddo

      do  i = 1,im
        do  j = 1,ic
            svdmat(i,j) = a(i,j)
        enddo
      enddo

      if(dbg.eq.1) then
         write(*,*) 'A0:'
         do  j = 1,im
            write(*,6003) (svdmat(j,i),i=1,ic)
         enddo
      endif

      call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,                 &
     &matv,vmat,amater,ws)

 6001 format(1X,'Corrector: ',I4,'   sing: ',F12.4)
 6002 format('VMAT: ',I4,I4,5X,F12.6,2X,F12.6)
 6003 format(16(2X,F7.2))
 6004 format(16(2X,F7.2))

      if(amater.ne.0) then
         write(*,*) 'end SVD with error code: ',amater
      endif

      do  i = 1,ic
          if(dbg.eq.1) then
!             write(*,*) i,wvec(i)
             write(*,6001) i,wvec(i)
          endif
          wmat(i,i) = wvec(i)
          if(abs(wvec(i)).gt.1.0001) then
                 wtmat(i,i) = 1/wvec(i)
          else
                 wtmat(i,i) = 0.0
          endif
      enddo

      if(dbg.eq.1) then
         call rvord(wvec,sortw,ws,ic)
         do  i = 1,ic
            write(*,*) i,sortw(i),wvec(sortw(i))
         enddo
      endif

      do  i = 1,ic
          do  j = 1,im
              vtmat(i,j) = vmat(j,i)
              utmat(i,j) = umat(j,i)
          enddo
      enddo

      if(dbg.eq.1) then
         write(*,*) 'A1:'
         do  j = 1,im
            write(*,6003) (svdmat(j,i),i=1,ic)
         enddo
         write(*,*) ' '
         
         write(*,*) 'Va:'
         do  j = 1,ic
            write(*,6004) (vmat(j,i),i=1,ic)
         enddo
         write(*,*) 'Vt:'
         do  j = 1,ic
            write(*,6004) (vtmat(j,i),i=1,ic)
         enddo
         write(*,*) 'W:'
         do  j = 1,im
            write(*,6004) (wmat(j,i),i=1,ic)
         enddo
         write(*,*) 'Wt:'
         do  j = 1,ic
            write(*,6004) (wtmat(j,i),i=1,im)
         enddo
         write(*,*) 'U:'
         do  j = 1,im
            write(*,6004) (umat(j,i),i=1,im)
         enddo
         write(*,*) 'Ut:'
         do  j = 1,im
            write(*,6004) (utmat(j,i),i=1,im)
         enddo
      endif

      call dmmpy(im,im,utmat(1,1),utmat(1,2),utmat(2,1),                &
           &xin(1),xin(2),xa(1),xa(2))
      call dmmpy(ic,im,wtmat(1,1),wtmat(1,2),wtmat(2,1),                &
           &xa(1),xa(2),xb(1),xb(2))
      call dmmpy(ic,ic,vmat(1,1),vmat(1,2),vmat(2,1),                   &
           &xb(1),xb(2),xc(1),xc(2))
      call dmmpy(im,ic,svdmat(1,1),svdmat(1,2),svdmat(2,1),             &
           &xc(1),xc(2),xpred(1),xpred(2))

      if(dbg.eq.1) then
         write(*,*) xc
         write(*,*) xpred
      endif

      do  i = 1,im
         xout(i) = xin(i) - xpred(i)
      enddo
      do  i = 1,ic
         xc(i) = -xc(i)
      enddo

      return
    end subroutine svdcorr_m

    subroutine svdcorr_c(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,   &
         &xin,xc,xout,xa,xb,xpred,ws,wvec,                            &
         &sortw,nx,im,ic,iflag,dbg)
! ******************************************************
!                                                      *
!    Performs SVD and correction for matrix with more  *
!    correctors than monitors.                         *
!                                                      *
!     Author: WFH  12.09.02                            *
!                                                      *
! ******************************************************
      implicit none
      integer im,ic,nx(ic),i,j
      integer iflag
      double precision a(im,ic)
      double precision svdmat(ic,im)
      double precision umat(ic,im), utmat(im,ic)
      double precision vmat(ic,im), vtmat(im,ic)
      double precision wmat(ic,im), wtmat(im,ic)
      double precision xin(im),xout(im),xc(ic)
      double precision xa(im),xpred(im),ws(ic),wvec(ic)
      integer sortw(ic)
      integer amater, svdmx, svdnx, svdnm
      logical matu, matv
      integer dbg
      double precision xb(2000)

      MATU = .TRUE.
      MATV = .TRUE.
      iflag = 0
      write(*,*) ' '
      write(*,*) 'start SVD correction using ',ic,' correctors'
      write(*,*) ' '

      SVDNM = max(ic,im)
      svdmx = ic
      svdnx = im

      do  i = 1,im
        do  j = 1,ic
            svdmat(j,i) = a(i,j)
        enddo
      enddo
      do  i = 1,ic
        nx(i) = i
      enddo

 8373 continue

      if(dbg.eq.1) then
          write(*,*) 'A0:'
          do  j = 1,ic
            write(*,6003) (svdmat(j,i),i=1,im)
          enddo
      endif

      call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,                 &
           &matv,vmat,amater,ws)

 6001 format(1X,'Corrector: ',I4,'   sing: ',F12.4)
 6002 format('UMAT: ',I4,I4,5X,F12.6,2X,F12.6)
 6003 format(16(2X,F7.2))
 6004 format(16(2X,F7.2))

      if(amater.ne.0) then
         write(*,*) 'end SVD with error code: ',amater
      endif

      do  i = 1,im
         if(dbg.eq.1) then
            write(*,*) i,wvec(I)
            write(*,6001) i,wvec(I)
         endif
         wmat(i,i) = wvec(i)
         if(abs(wvec(i)).gt.1.0001) then
            wtmat(i,i) = 1/wvec(i)
         else
            wtmat(i,i) = 0.0
          endif
       enddo

       if(dbg.eq.1) then
          call rvord(wvec,sortw,ws,im)
          do  i = 1,im
             write(*,*) i,sortw(i),wvec(sortw(i))
          enddo
       endif

       do  i = 1,im
          do  j = 1,ic
             vtmat(i,j) = vmat(j,i)
             utmat(i,j) = umat(j,i)
          enddo
      enddo

      if(dbg.eq.1) then
         write(*,*) 'A1:'
         do  j = 1,ic
            write(*,6003) (svdmat(j,i),i=1,im)
         enddo
         write(*,*) ' '

         write(*,*) 'Va:'
         do  j = 1,im
            write(*,6004) (vmat(j,i),i=1,im)
         enddo
         write(*,*) 'Vt:'
         do  j = 1,im
            write(*,6004) (vtmat(j,i),i=1,im)
         enddo
         write(*,*) 'W:'
         do  j = 1,ic
            write(*,6004) (wmat(j,i),i=1,im)
         enddo
         write(*,*) 'Wt:'
         do  j = 1,im
            write(*,6004) (wtmat(j,i),i=1,ic)
         enddo
         write(*,*) 'U:'
         do  j = 1,ic
            write(*,6004) (umat(j,i),i=1,ic)
         enddo
         write(*,*) 'Ut:'
         do  j = 1,ic
            write(*,6004) (utmat(j,i),i=1,ic)
         enddo
      endif
      
      call dmmpy(im,im,vtmat(1,1),vtmat(1,2),vtmat(2,1),                &
           &xin(1),xin(2),xa(1),xa(2))
      call dmmpy(im,im,wtmat(1,1),wtmat(1,2),wtmat(2,1),                &
           &xa(1),xa(2),xb(1),xb(2))
      call dmmpy(ic,im,umat(1,1),umat(1,2),umat(2,1),                   &
           &xb(1),xb(2),xc(1),xc(2))
      call dmmpy(im,ic,a(1,1),a(1,2),a(2,1),                            &
           &xc(1),xc(2),xpred(1),xpred(2))

      if(dbg.eq.1) then
         write(*,*) xc
         write(*,*) xpred
      endif

      do  i = 1,im
         xout(i) = xin(i) - xpred(i)
      enddo
      do  i = 1,ic
         xc(i) = -xc(i)
      enddo

      return
    end subroutine svdcorr_c


    subroutine solsql(m,n,xad,orb0,orbr,xinc,cb,xmeas,xres,y,z,xd)
      implicit none
!*********************************************************************
!     Subroutine SOLSQL to solve least sq. problem for orbit         *
!     after matrix has been reconditioned                            *
!                                                                    *
!     Authors:     WFH                 Date:  21.03.1995             *
!                                                                    *
!*********************************************************************
      integer m,n,i,j,ifail
      double precision xad(m,n),orb0(m),orbr(m),xinc(n),cb(n),xmeas(m), &
           &xres(m),y(n,m),z(n,n),xd(n),zero
      parameter(zero=0d0)

! --- copy from original matrix, sizes are not equal !
      do i = 1,m
         xmeas(i) = orb0(i)
         xres(i)  = zero
      enddo
!      do  i = 1,m
!        do  j = 1,n
!          write(*,*) i,j,xad(i,j)
!        enddo
!      enddo

      call lstsql(m,n,xad,xmeas,cb,xres,ifail,y,z,xd)
      write(*,*) 'IFAIL from lstsql: ',ifail

      do i=1,n
         xinc(i) = cb(i)
!        write(*,*) i,cb(i)
      enddo
      do j = 1,m
         orbr(j) = xres(j) + xmeas(j)
      enddo

 6001 format(3x,i4,2x,f10.5)
 6002 format(3x,i4,2x,f10.5,2x,f10.5,2x,f10.5)

      return
    end subroutine solsql

    subroutine lstsql(m,n,x,d,cb,xpred,ifail,y,z,xd)
      implicit none
!*********************************************************************
!     Subroutine LSTSQR to make a least sqared minimization          *
!                                                                    *
!     Authors:     WFH                 Date:  21.03.1995             *
!                                                                    *
!*********************************************************************
      integer m,n,ifail,i,j,w(200000)
      double precision x(m,n),d(m),cb(n),xpred(m),y(n,m),z(n,n),xd(n),  &
           &dw(100000)
      equivalence (w(1),dw(1))

      do  i = 1,m
         do  j = 1,n
            y(j,i) = x(i,j)
         enddo
      enddo

      call dmmlt(n,m,n,y,y(1,2),y(2,1),x,x(1,2),x(2,1),z,               &
           &z(1,2),z(2,1),dw)

      call dinv(n,z,n,w,ifail)

      call dmmpy(n,m,y(1,1),y(1,2),y(2,1),d(1),d(2),xd(1),xd(2))

      call dmmpy(n,n,z(1,1),z(1,2),z(2,1),xd(1),xd(2),cb(1),cb(2))

      do i=1,n
         cb(i)=-cb(i)
!        write(*,*) '=> ',i,cb(i)
      enddo

      call dmmpy(m,n,x(1,1),x(1,2),x(2,1),cb(1),cb(2),xpred(1),xpred(2))

      return
    end subroutine lstsql

    subroutine micado(a,conm,b,orbr,xinc,nx,rms,m,n,iter,rho,ptop,    &
         &rmss,xrms,xptp,xiter,ifail)
      implicit none
!*********************************************************************
!     Subroutine MICADO to run MICADO minimisation                   *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     interface between ORBCOR and HTLS routine
!     ARRAY and loop dimensions are transmitted as subroutine arguments
!     A(M,N)    : response matrix correctors ==>> monitors
!     B(M)         : orbit to be corrected
!     ORBR(M)      : residual orbit
!     XINC(N)      : strength of correctors
!     NX(N)        : sequence number of correctors
      integer m,n,nx(n),prtlev,iter
      integer ifail
      real a(m,n),b(m),orbr(m),xinc(n),rms,rho(3*n),ptop(n),rmss(n),    &
           &xrms(n),xptp(n),xiter(n)
      character*16 conm(n)

      prtlev = 3

      call htls(a,conm,b,m,n,xinc,nx,orbr,rms,prtlev,iter,rho,ptop,rmss,&
           &xrms,xptp,xiter,ifail)

! --- energy shift caused by corrector strength changes
!     (inhibited)
!     if(iplane.eq.2) go to 85
!     do  75 ij1=1,iter
!  75 dp=dp-apc(nx(ij1))*xinc(ij1)*mcfl
!     write(61,*)' DP NEW CORR=',dp

      return
    end subroutine micado

    subroutine htls(a,conm,b,m,n,x,ipiv,r,rms,prtlev,iter,rho,ptop,   &
         &rmss,xrms,xptp,xiter,ifail)
      implicit none
!*********************************************************************
!     Subroutine HTLS to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     dimension of array RHO should be 3*N
!     M  = NMTOT nr available monitors
!     N  = NCTOT nr available independent correctors
      logical interm
      integer m,n,ipiv(n),prtlev,iter,ij1,k2,k,i,kpiv,k3,j,ip,j1,kk,ki, &
           &iii,kkk
      real a(m,n),b(m),x(n),r(m),rms,rho(3*n),ptop(n),rmss(n),xrms(n),  &
           &xptp(n),xiter(n),xxcal,ptp,g,h,sig,beta,piv,pivt,rm,pt,rzero, &
           &reps7
      parameter(rzero=0e0,reps7=1e-7)
      character*4 units
      character*16 conm(n)
      integer      ifail

      interm = .true.
      ifail = 0
      units = 'mrad'
      ptp = rzero

      call calrms(b,m,rm,pt)

      if(rm.le.rms) then
         write(*,*) '++++++ WARNING: RMS already smaller than desired '
         write(*,*) '++++++ WARNING: no correction is done            '
         rms = rm
         iter = 0
         ifail = -2
         return
      endif

! --- calculate first pivot
!==========================

      do ij1=1,3*n
         rho(ij1)=rzero
      enddo

      k2=n + 1
      piv=rzero
      
      do k=1,n
         ipiv(k)=k
         h=rzero
         g=rzero
         do i=1,m
            h=h+a(i,k)*a(i,k)
            g=g+a(i,k)*b(i)
         enddo
         rho(k)=h
         rho(k2) = g
         if(h.ne.rzero) then
            pivt = g*g/h
         else
            pivt = rzero
         endif
         if(pivt-piv.gt.rzero) then
            piv = pivt
            kpiv=k
         endif
         k2 = k2 + 1
      enddo

! --- boucle pour chaque iteration

      do k=1,iter

!X    WRITE(*,*) ' Start iteration ',K
         if (kpiv.eq.k) go to 8
!X    write(*,*) 'change row, KPIV, K: ',KPIV, K

! --- on echange les K et KPIV si KPIV plus grand que K
        h=rho(k)
        rho(k)=rho(kpiv)
        rho(kpiv)=h
        k2=n+k
        k3=n+kpiv
        g = rho(k2)
        rho(k2) = rho(k3)
        rho(k3) = g
        do i=1,m
           h=a(i,k)
           a(i,k)=a(i,kpiv)
           a(i,kpiv)=h
        enddo

! --- calcul de beta,sigma et uk dans htul
    8   continue
!X    write(*,*) 'call HTUL'
        call htul(a,m,n,k,sig,beta)
!X    write(*,*) 'back from HTUL'

! --- on garde SIGMA dans RHO(N+K)
        j=n+k
        rho(j)=-sig
        ip=ipiv(kpiv)
        ipiv(kpiv)=ipiv(k)
        ipiv(k)=ip
        if(k.eq.n) go to 13

! --- transformation de A dans HTAL
!X    write(*,*) 'call HTAL'
        call htal(a,m,n,k,beta)
!X    write(*,*) 'back from HTAL'

! --- transformation de B dans HTBL
   13   continue
!X    write(*,*) 'call HTBL'
        call htbl(a,b,m,n,k,beta)
!X    write(*,*) 'back from HTBL'

! --- recherche du pivot (K+1)
!=============================

        rho(k)=sqrt(piv)
        if(k.eq.n) go to 11
        piv=rzero
        kpiv = k + 1
        j1 = kpiv
        k2=n + j1
!x    write(*,*) 'loop 18, ',j1,n
        do j=j1,n
           h=rho(j)-(a(k,j))*(a(k,j))
!X    write(*,*) 'K,J: ',K,J
!X    write(*,*) RHO(J),A(K,J)
!X    write(*,*) 'H: ',H

           if(abs(h).lt.reps7) then
              write(*,*) 'Correction process aborted'
              write(*,*) 'during ',k,'th iteration'
              write(*,*) 'Last r.m.s.: ',rmss(k-1)
              write(*,*) 'Last p-t-p.: ',ptop(k-1)
              write(*,*) 'Division by zero expected'
              write(*,*) 'Probably two kickers too close: ',h
              write(*,*) 'SUSPECTED KICKER: ',J,'  '
              write(*,*) conm(ipiv(j))
              ifail = -1
              return
!           stop 777
           endif

           rho(j)=h
           g=rho(k2)-(a(k,j))*(b(k))
           rho(k2) = g
           if(h.ne.rzero) then
              pivt = g*g/h
           else
              pivt = rzero
           endif
!X    write(*,*) 'compare for pivot K,J,PIVT,PIV: ',K,J,PIVT,PIV
           if(pivt.lt.piv)go to 18
!X    write(*,*) 'comparison succeeded, set pivot'
           kpiv=j
           piv=pivt
18         continue
           k2 = k2 + 1
        enddo

! --- calcul des x
11      x(k)=b(k)/rho(n+k)
        if(k.eq.1)go to 27
        do i=2,k
           kk=k-i+1
           x(kk)=b(kk)
           ki=kk+1
           do j=ki,k
              x(kk)=x(kk)-a(kk,j)*x(j)
           enddo
           x(kk)=x(kk)/rho(n+kk)
        enddo
!       write(*,*) 'after 15'
27      continue

! --- save residual orbit and inverse sign of corrections (convention!)
        do iii= 1,m
           r(iii) = b(iii)
        enddo
        do iii= 1,k
           x(iii) =-x(iii)
        enddo

! --- calcul du vecteur residuel dans htrl
!=========================================

!     transform orbit R back to "normal space"
!     write(*,*) 'transform back to normal space'
        call htrl(a,r,m,n,k,rho)
        call calrms(r,m,rmss(k),ptop(k))
!       WRITE(61,'(I10,2F15.2)')K,RMSS(K),PTOP(K)
        if(k.lt.n) then
           xiter(k+1) = k
           xrms(k+1)  = rmss(k)
           xptp(k+1)  = ptop(k)
        endif
!       write(*,*) 'orbit back in normal space ',prtlev,k

! --- write intermediate results to 61 files
        if(k.eq.1)then
           if (prtlev .ge. 2) then
              write(61,52)
              write(61,54)units
              write(61,'(4x,"0",42x,f12.8,f15.8)')rm,pt
           endif
52         FORMAT(/' ***********    start MICADO    ***********'/)
54         FORMAT(' iter',5X,'corrector',13X,A4,6X,'mrad',               &
                &5X,"  rms",10X," ptop",/)
        endif

        if (prtlev .ge. 2) then
           write(61,'(/,1x,72("-"),/,                                     &
                &1x,i4,3x,38(" "),1x,f12.8,f15.8)')k,rmss(k),ptop(k)
        endif

        do kkk = 1,k
           xxcal=x(kkk)
           if (prtlev .ge. 2) then
              write(61,'(I3,1X,A16,9x,f8.4,2x,f8.4)')                      &
                   &kkk,conm(ipiv(kkk)),x(kkk),xxcal
           endif
        enddo

        if (interm) then
           if (prtlev .ge. 2) then
              write(61,58)k
              write(61,'(1x,8f9.3)')(r(kkk),kkk=1,m)
           endif
58         format(/,' residual orbit after iteration ',i2,':')
        endif

        if(k.eq.iter) then
           if (prtlev .ge. 2) then
              write(61,53)
          endif
        endif
 53     format(/' ***********    end   MICADO    ***********'/)


        if(ptop(k).le.ptp)go to 202
        if(rmss(k).le.rms)go to 202
     enddo
     return
     
! --- correction is already good enough:
!=======================================

202  ptp=ptop(k)
     rms=rmss(k)
     iter=k
     return
   end subroutine htls


   subroutine htal(a,m,n,k,beta)
     implicit none
!*********************************************************************
!     Subroutine HTAL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     Householder transform of matrix A
      integer m,n,nc,k,j,k1
      real beta,a(m,n),h,rzero
      parameter(rzero=0e0)

      nc=n-k

      do j=1,nc
         h=rzero
         
         do k1=k,m
            h=h+a(k1,k)*a(k1,k+j)
         enddo
         
         h=beta*h
         do k1=k,m
            a(k1,k+j)=a(k1,k+j)-a(k1,k)*h
         enddo
      enddo
      
    end subroutine htal


    subroutine htbl(a,b,m,n,k,beta)
      implicit none
!*********************************************************************
!     Subroutine HTBL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     Householder transform of vector B
      integer m,n,k,k1
      real beta,a(m,n),b(m),h,rzero
      parameter(rzero=0e0)

      h=rzero

      do k1=k,m
         h=h+a(k1,k)*b(k1)
      enddo
      
      h=beta*h
      
      do k1=k,m
         b(k1)=b(k1)-a(k1,k)*h
      enddo

    end subroutine htbl


    subroutine htrl(a,b,m,n,k,rho)
      implicit none
!*********************************************************************
!     Subroutine HTRL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     calculate residual orbit vector
      integer m,n,k,kk,kl,i,lv,kn
      real beta,a(m,n),b(m),rho(3*n),rzero
      parameter(rzero=0e0)

      do i= 1,k,1
         b(i)=rzero
      enddo

      do kk=1,k
         lv=m-k+kk
         kn=n+k-kk+1
         kl=k-kk+1
         
         beta=-1d0/(rho(kn)*a(kl,kl))
         call htbl(a,b,m,n,kl,beta)
      enddo

    end subroutine htrl


    subroutine htul(a,m,n,k,sig,beta)
      implicit none
!*********************************************************************
!     Subroutine HTUL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     calculate vector U
      integer m,n,k,i
      real beta,sig,a(m,n),h,rzero
      parameter(rzero=0e0)

      sig=rzero

      do i=k,m
         sig=sig+a(i,k)* a(i,k)
      enddo
      
      sig=sqrt(sig)
!     on choisit le signe correct pour sig:
      h=a(k,k)
      if (h.lt.rzero) sig=-sig
      beta=h + sig
      a(k,k)=beta
      beta=1d0/(sig*beta)
    end subroutine htul


    subroutine lequ2(a,ifail,b)
      implicit none
!*********************************************************************
!     Subroutine LEQU2 to solve X * A = b                            *
!     for 2 dimensions                                               *
!                                                                    *
!     Authors:     WFH                 Date:  26.02.1992             *
!                                                                    *
!*********************************************************************
      integer ifail
      real a(2,2),b(2),x(2),det,reps6
      parameter(reps6=1e-6)

      det = a(2,2)*a(1,1) - a(1,2)*a(2,1)
      if(abs(det).lt.reps6) then
         ifail = -1
         return
      endif

      x(1) = (b(1)*a(2,2) - b(2)*a(1,2))/det
      x(2) = (b(2) - a(2,1)*x(1))/a(2,2)
      
      b(1) = x(1)
      b(2) = x(2)
      ifail = 0
      return
    end subroutine lequ2


    subroutine calrms(r,m,rms,ptp)
      implicit none
!*********************************************************************
!     Subroutine CALRMS to calculate rms                             *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!*********************************************************************
!     calculates rms and p.to.p value of R(1) .... R(M)
      integer m,i,imax,imin,maxmin
      real r(m),xave,xrms,rms,ptp,ave,rzero
      parameter(rzero=0e0)

      xave = rzero
      xrms = rzero
      
      do i=1,m
         xave = xave + r(i)
         xrms = xrms + (r(i)*r(i))
      enddo
      
      ave = xave / float(m)
      rms = xrms / float(m)
      
      imax=maxmin(r(1),m,1)
      imin=maxmin(r(1),m,0)
      ptp=r(imax)-r(imin)
      rms=sqrt(rms)
      return
    end subroutine calrms


    function maxmin (a,n,m)
      implicit none
!*********************************************************************
!     Subroutine MAXMIN to find maximum and minimum of a list        *
!                                                                    *
!     Authors:     WFH                 Date:  21.08.2000             *
!                                                                    *
!*********************************************************************
!     if M=0, MAXMIN=lowest index of minimum element in A
!     if M=1, MAXMIN=lowest index of maximun element in A
!     if N<1, MAXMIN=1
      integer m,n,maxmin,i
      real a(n),curent

      maxmin=1
      if (n.lt.1) return
      curent=a(1)
      do i=2,n
         if ((m.eq.0).and.(a(i).ge.curent)) go to 10
         if ((m.eq.1).and.(a(i).le.curent)) go to 10
         curent=a(i)
         maxmin=i
10       continue
      enddo
      return
    end function maxmin

    subroutine dinv(n,a,idim,r,ifail)
      implicit none
!******************************************************************
!
!     REPLACES A BY ITS INVERSE.
!
!     (PARAMETERS AS FOR DEQINV.)
!
!     CALLS ... DFACT, DFINV, F010PR, ABEND.
!
!******************************************************************
      integer n,idim,ifail,jfail
      integer r(n),t1,t2,t3
      double precision a(idim,n),det,temp,s,c11,c12,c13,c21,c22,c23,c31,&
           &c32,c33,zero
      parameter(zero=0d0)

!
!  TEST FOR PARAMETER ERRORS.
!
      if((n.lt.1).or.(n.gt.idim)) go to 7
!
!  TEST FOR N.LE.3.
!
      if(n.gt.3) go to 6
      ifail=0
      if(n.lt.3) go to 4
!
!  N=3 CASE.
!
!     COMPUTE COFACTORS.
      c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
      c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
      c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
      c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      t1=abs(sngl(a(1,1)))
      t2=abs(sngl(a(2,1)))
      t3=abs(sngl(a(3,1)))
!
!     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      if(t1.ge.t2) go to 1
      if(t3.ge.t2) go to 2
!        (PIVOT IS A21)
      temp=a(2,1)
      det=c13*c32-c12*c33
      go to 3

1     if(t3.ge.t1) go to 2
!     (PIVOT IS A11)
      temp=a(1,1)
      det=c22*c33-c23*c32
      go to 3

!     (PIVOT IS A31)
2     temp=a(3,1)
      det=c23*c12-c22*c13
!
!     SET ELEMENTS OF INVERSE IN A.
3     if(det.eq.zero) go to 8
      s=temp/det
      a(1,1)=s*c11
      a(1,2)=s*c21
      a(1,3)=s*c31
      a(2,1)=s*c12
      a(2,2)=s*c22
      a(2,3)=s*c32
      a(3,1)=s*c13
      a(3,2)=s*c23
      a(3,3)=s*c33
      return
!
4     if(n.lt.2) go to 5
!
!  N=2 CASE BY CRAMERS RULE.
!
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(det.eq.zero) go to 8
      s=1d0/det
      c11   =s*a(2,2)
      a(1,2)=-s*a(1,2)
      a(2,1)=-s*a(2,1)
      a(2,2)=s*a(1,1)
      a(1,1)=c11
      return
!
!  N=1 CASE.
!
5     if(a(1,1).eq.zero) go to 8
      a(1,1)=1d0/a(1,1)
      return
!
!  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
!
6     call dfact(n,a,idim,r,ifail,det,jfail)
      if(ifail.ne.0) return
      call dfinv(n,a,idim,r)
      return
!
!  ERROR EXITS.
!
7     ifail=+1
      return
!
8     ifail=-1
      return
!
    end subroutine dinv

    subroutine dfact(n,a,idim,ir,ifail,det,jfail)
      implicit none
      integer idim,j,k,n,ifail,nxch,jp1,i,l,jm1,jfail,ir(*),normal,     &
           &imposs,jrange,jover,junder
      parameter(normal=0,imposs=-1,jrange=0,jover=1,junder=-1)
      real p,q,t,g1,g2
      parameter(g1=1e-19,g2=1e19)
      double precision a(idim,*),det,tf,s11,s12,zero,one
      parameter(zero=0d0,one=1d0)
      character*6 hname
      data hname /  ' DFACT'  /

      if(idim .ge. n  .and.  n .gt. 0)  goto 110
      call tmprnt(hname,n,idim,0)
      return

110   ifail  =  normal
      jfail  =  jrange
      nxch   =  0
      det    =  one
      do    j  =  1, n
120      k  =  j
         p  =  abs(sngl(a(j,j)))
         if(j .eq. n)  goto 122
         jp1  =  j+1
         do    i  =  jp1, n
            q  =  abs(sngl(a(i,j)))
            if(q .le. p)  goto 121
            k  =  i
            p  =  q
121         continue
         enddo
         if(k .ne. j)  goto 123
122      if(p .gt. zero)  goto 130
         det    =  zero
         ifail  =  imposs
         jfail  =  jrange
         return
123      do    l  =  1, n
            tf      =  a(j,l)
            a(j,l)  =  a(k,l)
            a(k,l)  =  tf
         enddo
         nxch      =  nxch + 1
         ir(nxch)  =  j*2**12 + k
130      det     =  det * a(j,j)
         a(j,j)  =  one / a(j,j)
         t  =  abs(sngl(det))
         if(t .lt. g1)  then
            det    =  zero
            if(jfail .eq. jrange)  jfail  =  junder
         elseif(t .gt. g2)  then
            det    =  one
            if(jfail .eq. jrange)  jfail  =  jover
         endif
         if(j .eq. n)  goto 144
         jm1  =  j-1
         jp1  =  j+1
         do   k  =  jp1, n
            s11  =  -a(j,k)
            s12  =  -a(k,j+1)
            if(j .eq. 1)  goto 142
            do  i  =  1, jm1
               s11  =  a(i,k)*a(j,i)+s11
               s12  =  a(i,j+1)*a(k,i)+s12
            enddo
142         a(j,k)    =  -s11 * a(j,j)
            a(k,j+1)  =  -(a(j,j+1)*a(k,j)+s12)
         enddo
144      continue
      enddo
150   if(mod(nxch,2) .ne. 0)  det  =  -det
      if(jfail .ne. jrange)   det  =  zero
      ir(n)  =  nxch
      return
    end subroutine dfact

    subroutine dfeqn(n,a,idim,ir,k,b)
      implicit none
      integer idim,i,j,k,n,m,nxch,ij,l,im1,nm1,nmi,nmjp1,ir(*)
      double precision a(idim,*),b(idim,*),te,s21, s22
      character*6 hname
      data hname /  ' DFEQN'  /

      if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0)  goto 210

      call tmprnt(hname,n,idim,k)
      return

210   nxch  =  ir(n)
      if(nxch .eq. 0)  goto 220
      do    m  =  1, nxch
         ij  =  ir(m)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do   l  =  1, k
            te      =  b(i,l)
            b(i,l)  =  b(j,l)
            b(j,l)  =  te
         enddo
      enddo
220   do    l  =  1, k
         b(1,l)  =  a(1,1)*b(1,l)
      enddo
      if(n .eq. 1)  goto 299
      do    l  =  1, k
         do   i  =  2, n
            im1  =  i-1
            s21  =  - b(i,l)
            do   j  =  1, im1
               s21  =  a(i,j)*b(j,l)+s21
            enddo
            b(i,l)  =  - a(i,i)*s21
         enddo
         nm1  =  n-1
         do   i  =  1, nm1
            nmi  =  n-i
            s22  =  - b(nmi,l)
            do   j  =  1, i
               nmjp1  =  n - j+1
               s22    =  a(nmi,nmjp1)*b(nmjp1,l)+s22
            enddo
            b(nmi,l)  =  - s22
         enddo
      enddo
299   continue
      return
    end subroutine dfeqn

    subroutine dfinv(n,a,idim,ir)
      implicit none
      integer idim,i,j,k,n,m,nxch,ij,nm1,nmi,im2,ir(*)
      double precision a(idim,*),ti,s31,s32,s33,s34,zero
      parameter(zero=0d0)
      character*6 hname
      data hname /  ' DFINV'  /

      if(idim .ge. n  .and.  n .gt. 0)  goto 310

      call tmprnt(hname,n,idim,0)
      return

310   if(n .eq. 1)  return
      a(2,1)  =  -a(2,2) * a(1,1) * a(2,1)
      a(1,2)  =  -a(1,2)
      if(n .eq. 2)  goto 330
      do    i  =  3, n
         im2  =  i-2
         do j  =  1, im2
            s31  =  zero
            s32  =  a(j,i)
            do  k  =  j, im2
               s31  =  a(k,j)*a(i,k)+s31
               s32  =  a(j,k+1)*a(k+1,i)+s32
            enddo
            a(i,j)  =  -a(i,i) * (a(i-1,j)*a(i,i-1)+s31)
            a(j,i)  =  -s32
         enddo
         a(i,i-1)  =  -a(i,i) * a(i-1,i-1) * a(i,i-1)
         a(i-1,i)  =  -a(i-1,i)
      enddo
330   nm1  =  n-1
      do   i  =  1, nm1
         nmi  =  n-i
         do   j  =  1, i
            s33  =  a(i,j)
            do   k  =  1, nmi
               s33  =  a(i+k,j)*a(i,i+k)+s33
            enddo
            a(i,j)  =  s33
         enddo
         do   j  =  1, nmi
            s34  =  zero
            do   k  =  j, nmi
               s34  =  a(i+k,i+j)*a(i,i+k)+s34
            enddo
            a(i,i+j)  =  s34
         enddo
      enddo
      nxch  =  ir(n)
      if(nxch .eq. 0)  return
      do m  =  1, nxch
         k   =  nxch - m+1
         ij  =  ir(k)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do  k  =  1, n
            ti      =  a(k,i)
            a(k,i)  =  a(k,j)
            a(k,j)  =  ti
         enddo
      enddo
      return
    end subroutine dfinv

    subroutine tmprnt(name,n,idim,k)
      implicit none
      logical mflag,rflag
      integer idim,k,n,ifmt
      character*6 name

      mflag=.true.
      rflag=.true.
      if(name(3:6) .eq. 'FEQN') ifmt=1001
      if(name(3:6) .ne. 'FEQN') ifmt=1002
      if(mflag) then
         if(name(3:6) .eq. 'feqn') then
            if(ifmt.eq.1001) write(*,1001) name, n, idim, k
            if(ifmt.eq.1002) write(*,1002) name, n, idim, k
         else
            if(ifmt.eq.1001) write(*,1001) name, n, idim
            if(ifmt.eq.1002) write(*,1002) name, n, idim
         endif
      endif
      if(.not. rflag) call abend
      return

1001  FORMAT(7X," PARAMETER ERROR IN SUBROUTINE ", A6,                  &
           &" ... (N.LT.1 OR IDIM.LT.N).",                              &
           &5X,"N =", I4, 5X,"IDIM =", I4,".")
1002  FORMAT(7X," PARAMETER ERROR IN SUBROUTINE ", A6,                  &
           &" ... (N.LT.1 OR IDIM.LT.N OR K.LT.1).",                    &
           &5X,"N =", I4, 5X,"IDIM =", I4, 5X,"K =", I4,".")
    end subroutine tmprnt

    subroutine abend
      implicit none
      write(*,*) 'Abnormal end ...'
      stop 888
    end subroutine abend

    subroutine dmmpy(m,n,x,x12,x21,y,y2,z,z2)
      implicit none
      integer m,n,ix,jx,jy,iz,lxi1,lzi,i,lxij,lyj,j,locf
      double precision x(*),x12(*),x21(*),y(*),y2(*),z(*),z2(*),sum,zero
      parameter(zero=0d0)

      if(m .le. 0  .or.  n .le. 0)  return
      ix  =  (locf(x21) - locf(x)) / 2
      jx  =  (locf(x12) - locf(x)) / 2
      jy  =  (locf(y2) - locf(y)) / 2
      iz  =  (locf(z2)  - locf(z)) / 2
      lxi1  =  1
      lzi   =  1
      do     i  =  1, m
         lxij  =  lxi1
         lyj   =  1
         sum   =  zero
         do  j  =  1, n
            sum  =  x(lxij)*y(lyj)+sum
            lxij =  lxij + jx
            lyj  =  lyj + jy
         enddo
         z(lzi)  =  sum
         lxi1    =  lxi1 + ix
         lzi     =  lzi + iz
      enddo
      return
    end subroutine dmmpy

    subroutine dmmlt(m,n,k,x,x12,x21,y,y12,y21,z,z12,z21,t)
      implicit none
      integer m,n,k,ix,jx,jy,iz,lz,ly1l,lz1l,l,lxi1,lzil,i,lxij,lyjl,   &
           &j,locf,lzii,lxk1,lzik,kdash,lxkj,ltl,lxil,lti,lyil,lxii,ltk,lxik, &
           &lxki,locx,locy,ly,lzki,locz
      double precision x(*),x12(*),x21(*),y(*),y12(*),y21(*),z(*),      &
           &z12(*),z21(*),t(*),s11,s21,s22,s31,s41,s51,s52,zero
      parameter(zero=0d0)

      if(min0(m,n,k) .le. 0)  return
      locx  =  locf(x(1))
      locy  =  locf(y(1))
      locz  =  locf(z(1))
      ix  =  (locf(x21(1)) - locx) / 2
      jx  =  (locf(x12(1)) - locx) / 2
      jy  =  (locf(y21(1)) - locy) / 2
      ly  =  (locf(y12(1)) - locy) / 2
      iz  =  (locf(z21(1)) - locz) / 2
      lz  =  (locf(z12(1)) - locz) / 2
      if(locz .eq. locx)  goto 30
      if(locz .eq. locy)  goto 40
      if(locx .eq. locy)  goto 20
10    ly1l  =  1
      lz1l  =  1
      do     l  =  1, k
         lxi1  =  1
         lzil  =  lz1l
         do  i  =  1, m
            s11   =  zero
            lxij  =  lxi1
            lyjl  =  ly1l
            do  j  =  1, n
               s11   =  x(lxij)*y(lyjl)+s11
               lxij  =  lxij + jx
               lyjl  =  lyjl + jy
            enddo
            z(lzil)  =  s11
            lxi1     =  lxi1 + ix
            lzil     =  lzil + iz
         enddo
         ly1l  =  ly1l + ly
         lz1l  =  lz1l + lz
      enddo
      return

20    if(m .ne. k  .or.  ix .ne. ly  .or.  jx .ne. jy)  goto 10
      lxi1  =  1
      lzii  =  1
      do     i  =  1, m
         s21   =  zero
         lxij  =  lxi1
         do  j  =  1, n
            s21   =  x(lxij)*x(lxij)+s21
            lxij  =  lxij + jx
         enddo
         z(lzii)  =  s21
         if(i .eq. m)  goto 24
         lxk1  =  lxi1 + ix
         lzik  =  lzii + lz
         lzki  =  lzii + iz
         do  kdash  =  i+1, m
            s22   =  zero
            lxij  =  lxi1
            lxkj  =  lxk1
            do  j  =  1, n
               s22   =  x(lxij)*x(lxkj)+s22
               lxij  =  lxij + jx
               lxkj  =  lxkj + jx
            enddo
            z(lzik)  =  s22
            z(lzki)  =  z(lzik)
            lxk1  =  lxk1 + ix
            lzik  =  lzik + lz
            lzki  =  lzki + iz
         enddo
         lxi1  =  lxi1 + ix
         lzii  =  lzii + iz + lz
24       continue
      enddo
      return

30    if(locx .eq. locy)  goto 50
      lxi1  =  1
      do     i  =  1, m
         ly1l  =  1
         ltl   =  1
         do  l  =  1, k
            s31   =  zero
            lxij  =  lxi1
            lyjl  =  ly1l
            do  j  =  1, n
               s31   =  x(lxij)*y(lyjl)+s31
               lxij  =  lxij + jx
               lyjl  =  lyjl + jy
            enddo
            t(ltl)  =  s31
            ly1l    =  ly1l + ly
            ltl     =  ltl + 1
         enddo
         lxil  =  lxi1
         ltl   =  1
         do  l  =  1, k
            x(lxil)  =  t(ltl)
            lxil     =  lxil + jx
            ltl      =  ltl + 1
         enddo
         lxi1  =  lxi1 + ix
      enddo
      return

40    ly1l  =  1
      do     l  =  1, k
         lxi1  =  1
         lti   =  1
         do  i  =  1, m
            s41   =  zero
            lxij  =  lxi1
            lyjl  =  ly1l
            do  j  =  1, n
               s41   =  x(lxij)*y(lyjl)+s41
               lxij  =  lxij + jx
               lyjl  =  lyjl + jy
            enddo
            t(lti)  =  s41
            lxi1    =  lxi1 + ix
            lti     =  lti + 1
         enddo
         lyil  =  ly1l
         lti   =  1
         do  i  =  1, m
            y(lyil)  =  t(lti)
            lyil     =  lyil + jy
            lti      =  lti + 1
         enddo
         ly1l  =  ly1l + ly
      enddo
      return

50    lxi1  =  1
      lxii  =  1
      do     i  =  1, m
         s51   =  zero
         lxij  =  lxi1
         do  j  =  1, n
            s51   =  x(lxij)*x(lxij)+s51
            lxij  =  lxij + jx
         enddo
         t(1)  =  s51
         if(i .eq. m)  goto 54
         lxk1  =  lxi1 + ix
         ltk  =  2
         do  kdash  =  i+1, m
            s52   =  zero
            lxij  =  lxi1
            lxkj  =  lxk1
            do  j  =  1, n
               s52   =  x(lxij)*x(lxkj)+s52
               lxij  =  lxij + jx
               lxkj  =  lxkj + jx
            enddo
            t(ltk)  =  s52
            lxk1    =  lxk1 + ix
            ltk     =  ltk + 1
         enddo
54       lxik  =  lxii
         ltk   =  1
        do  kdash  =  i, m
           x(lxik)  =  t(ltk)
           lxik     =  lxik + jx
           ltk      =  ltk + 1
        enddo
        lxi1     =  lxi1 + ix
        lxii     =  lxii + ix + jx
     enddo
     if(m .eq. 1)  return
     lxii  =  1
     do     i  =  1, m-1
        lxik  =  lxii + jx
        lxki  =  lxii + ix
        do  kdash  =  i+1, m
           x(lxki)  =  x(lxik)
           lxik     =  lxik + jx
           lxki     =  lxki + ix
        enddo
        lxii  =  lxii + ix + jx
     enddo
     return
   end subroutine dmmlt


      SUBROUTINE SVD(NM,M,N,A,W,MATU,U,MATV,V,IERR,RV1)
      implicit  none
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      double precision c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag
      logical matu,matv
! ------------------------------------------------------------------
!     this subroutine is a translation of the algol procedure svd,
!     NUM. MATH. 14, 403-420(1970) by golub and reinsch.
!     handbook for auto. comp., vol 2 -linear algebra, 134-151(1971).
!
!     this subroutine determines the singular value decomposition
!          t
!     a=usv  of a real m by n rectangular matrix.  householder
!     bidiagonalization and a variant of the qr algorithm are used.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.
!
!        m is the number of rows of a (and u).
!
!        n is the number of columns of a (and u) and the order of v.
!
!        a contains the rectangular input matrix to be decomposed.
!
!        matu should be set to .true. if the u matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!        matv should be set to .true. if the v matrix in the
!          decomposition is desired, and to .false. otherwise.
!
!     on output
!
!        a is unaltered (unless overwritten by u or v).
!
!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.
!
!        u contains the matrix u (orthogonal column vectors) of the
!          decomposition if matu has been set to .true.  otherwise
!          u is used as a temporary array.  u may coincide with a.
!          if an error exit is made, the columns of u corresponding
!          to indices of correct singular values should be correct.
!
!        v contains the matrix v (orthogonal) of the decomposition if
!          matv has been set to .true.  otherwise v is not referenced.
!          v may also coincide with a if u is not needed.  if an error
!          exit is made, the columns of v corresponding to indices of
!          correct singular values should be correct.
!
!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
! ------------------------------------------------------------------
!
      ierr = 0
!
      do i = 1, m
        do j = 1, n
          u(i,j) = a(i,j)
        enddo
      enddo
!     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
!
      do i = 1, n
        l = i + 1
        rv1(i) = scale * g
        g = 0.0d0
        s = 0.0d0
        scale = 0.0d0
        if (i .gt. m) go to 210
!
        do k = i, m
          scale = scale + dabs(u(k,i))
        enddo
!
        if (scale .eq. 0.0d0) go to 210
!
        do k = i, m
          u(k,i) = u(k,i) / scale
          s = s + u(k,i)**2
        enddo
!
        f = u(i,i)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,i) = f - g
        if (i .eq. n) go to 190
!
        do j = l, n
          s = 0.0d0
!
          do k = i, m
            s = s + u(k,i) * u(k,j)
          enddo
!
          f = s / h
!
          do k = i, m
            u(k,j) = u(k,j) + f * u(k,i)
          enddo
        enddo
!
  190   do k = i, m
          u(k,i) = scale * u(k,i)
        enddo
!
  210   w(i) = scale * g
        g = 0.0d0
        s = 0.0d0
        scale = 0.0d0
        if (i .gt. m .or. i .eq. n) go to 290
!
        do k = l, n
          scale = scale + dabs(u(i,k))
        enddo
!
        if (scale .eq. 0.0d0) go to 290
!
        do k = l, n
          u(i,k) = u(i,k) / scale
          s = s + u(i,k)**2
        enddo
!
        f = u(i,l)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,l) = f - g
!
        do k = l, n
          rv1(k) = u(i,k) / h
        enddo
!
        if (i .eq. m) go to 270
!
        do j = l, m
          s = 0.0d0
!
          do k = l, n
            s = s + u(j,k) * u(i,k)
          enddo
!
          do k = l, n
            u(j,k) = u(j,k) + s * rv1(k)
          enddo
        enddo
!
  270   do k = l, n
          u(i,k) = scale * u(i,k)
        enddo
!
  290   x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
      enddo
!     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
!     .......... for i=n step -1 until 1 do -- ..........
      do ii = 1, n
        i = n + 1 - ii
        if (i .eq. n) go to 390
        if (g .eq. 0.0d0) go to 360
!
        do j = l, n
!     .......... double division avoids possible underflow ..........
          v(j,i) = (u(i,j) / u(i,l)) / g
        enddo
!
        do j = l, n
          s = 0.0d0
!
          do k = l, n
            s = s + u(i,k) * v(k,j)
          enddo
!
          do k = l, n
            v(k,j) = v(k,j) + s * v(k,i)
          enddo
        enddo
!
  360   do j = l, n
          v(i,j) = 0.0d0
          v(j,i) = 0.0d0
        enddo
!
  390   v(i,i) = 1.0d0
        g = rv1(i)
        l = i
      enddo
!     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
!     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
!
      do ii = 1, mn
        i = mn + 1 - ii
        l = i + 1
        g = w(i)
        if (i .eq. n) go to 430
!
        do j = l, n
          u(i,j) = 0.0d0
        enddo
!
  430   if (g .eq. 0.0d0) go to 475
        if (i .eq. mn) go to 460
!
        do j = l, n
          s = 0.0d0
!
          do k = l, m
            s = s + u(k,i) * u(k,j)
          enddo
!     .......... double division avoids possible underflow ..........
          f = (s / u(i,i)) / g
!
          do k = i, m
            u(k,j) = u(k,j) + f * u(k,i)
          enddo
        enddo
!
  460   do j = i, m
          u(j,i) = u(j,i) / g
        enddo
!
        go to 490
!
  475   do j = i, m
          u(j,i) = 0.0d0
        enddo
!
  490   u(i,i) = u(i,i) + 1.0d0
      enddo
!     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
!     .......... for k=n step -1 until 1 do -- ..........
      do kk = 1, n
        k1 = n - kk
        k = k1 + 1
        its = 0
!     .......... test for splitting.
!                for l=k step -1 until 1 do -- ..........
  520   do ll = 1, k
          l1 = k - ll
          l = l1 + 1
          tst2 = tst1 + dabs(rv1(l))
          if (tst2 .eq. tst1) go to 565
!     .......... rv1(1) is always zero, so there is no exit
!                through the bottom of the loop ..........
          tst2 = tst1 + dabs(w(l1))
          if (tst2 .eq. tst1) go to 540
        enddo
!     .......... cancellation of rv1(l) if l greater than 1 ..........
  540   c = 0.0d0
        s = 1.0d0
!
        do i = l, k
          f = s * rv1(i)
          rv1(i) = c * rv1(i)
          tst2 = tst1 + dabs(f)
          if (tst2 .eq. tst1) go to 565
          g = w(i)
          h = pythag(f,g)
          w(i) = h
          c = g / h
          s = -f / h
          if (.not. matu) go to 560
!
          do j = 1, m
            y = u(j,l1)
            z = u(j,i)
            u(j,l1) = y * c + z * s
            u(j,i) = -y * s + z * c
          enddo
!
  560     continue
        enddo
!     .......... test for convergence ..........
  565   z = w(k)
        if (l .eq. k) go to 650
!     .......... shift from bottom 2 by 2 minor ..........
        if (its .eq. 30) go to 1000
        its = its + 1
        x = w(l)
        y = w(k1)
        g = rv1(k1)
        h = rv1(k)
        f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
        g = pythag(f,1.0d0)
        f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
!     .......... next qr transformation ..........
        c = 1.0d0
        s = 1.0d0
!
        do i1 = l, k1
          i = i1 + 1
          g = rv1(i)
          y = w(i)
          h = s * g
          g = c * g
          z = pythag(f,h)
          rv1(i1) = z
          c = f / z
          s = h / z
          f = x * c + g * s
          g = -x * s + g * c
          h = y * s
          y = y * c
          if (.not. matv) go to 575
!
          do j = 1, n
            x = v(j,i1)
            z = v(j,i)
            v(j,i1) = x * c + z * s
            v(j,i) = -x * s + z * c
          enddo
!
  575     z = pythag(f,h)
          w(i1) = z
!     .......... rotation can be arbitrary if z is zero ..........
          if (z .eq. 0.0d0) go to 580
          c = f / z
          s = h / z
  580     f = c * g + s * y
          x = -s * g + c * y
          if (.not. matu) go to 600
!
          do j = 1, m
            y = u(j,i1)
            z = u(j,i)
            u(j,i1) = y * c + z * s
            u(j,i) = -y * s + z * c
          enddo
!
  600     continue
        enddo
!
        rv1(l) = 0.0d0
        rv1(k) = f
        w(k) = x
        go to 520
!     .......... convergence ..........
  650   if (z .ge. 0.0d0) go to 700
!     .......... w(k) is made non-negative ..........
        w(k) = -z
        if (.not. matv) go to 700
!
        do j = 1, n
          v(j,k) = -v(j,k)
        enddo
!
  700   continue
      enddo
!
      go to 1001
!     .......... set error -- no convergence to a
!                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
    end SUBROUTINE SVD


    DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      implicit  none
      double precision a,b
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
10    continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10

20    pythag = p

      return
    end FUNCTION PYTHAG

    subroutine rvord(inv,outv,ws,n)
      implicit   none
      integer    n
      double precision   inv(n), ws(n)
      integer  i,j,jmax
      integer  outv(n)

      do  i = 1,n
         ws(i) = inv(i)
      enddo

      do  j = 1,n
         jmax = 1
         do  i = 1,n
            if(ws(i).gt.ws(jmax)) then
               jmax = i
            endif
         enddo
         outv(n-j+1) = jmax
         ws(jmax) = 0.0
      enddo
      
      return
    end subroutine rvord

    subroutine primat(a,nc,nm)
      implicit   none
      integer i, j
      integer nm,nc
      integer a(nc, nm)

      do i = 1,nc
         write(*,*) (a(i,j),j=1,nm)
      enddo
      
      return
    end subroutine primat

    subroutine prdmat(a,nc,nm)
      implicit   none
      integer i, j
      integer nm,nc
      double precision a(nc, nm)

      do i = 1,nc
         write(*,*) (a(i,j),j=1,nm)
      enddo
      
      return
    end subroutine prdmat
