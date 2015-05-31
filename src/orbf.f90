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

    integer im,ic,nm,nc 
    double precision resp,a(nm, nc)

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

    a(im+1,ic+1) = resp
 
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
    ! RMS = value of tolerance for correction

    implicit none

    integer :: im, ic, iter, i, j
    integer :: nx(ic),ny(ic)
    real ax(im,ic), cinx(ic), xinx(im), resx(im), rho(3*ic), ptop(ic)
    real rmss(ic), xrms(ic), xptp(ic), xiter(ic)

    real rms, rzero, pt, rm
    parameter(rzero=0e0)

    double precision a(im,ic),xin(im),cin(ic),res(im)
    character(16) :: conm(ic)

    integer      n
    integer      ifail
    real         calrms

    ! translate corrector names to fortran strings
    do  j = 1,ic
       call f_ctof(n, conm(j), 16)
    enddo

    AX  = A
    CINX = rzero

    NY(1:ic) = (/ (i, i = 1, ic) /) ! NY(i) = i 

    XINX = XIN
    RESX = rzero

    write(*,*) ' '
    write(*,*) 'start MICADO correction with ',iter,' correctors'
    write(*,*) ' '

    rm = calrms(xinx,im)

    if(rm.le.rms) then
       write(*,*) '++++++ WARNING: RMS already smaller than desired '
       write(*,*) '++++++ WARNING: no correction is done            '
       rms = rm
       iter = 0
       ifail = -2
    else
       open(61,file='fort.61')     
       call htls(ax,conm,xinx,im,ic,cinx,ny,resx,rms,3,iter,rho,ptop,rmss,&
            &xrms,xptp,xiter,ifail)
       close(61)
    endif

    CIN = CINX
    RES = RESX
    NX(NY(1:ic)) = (/ (i, i = 1, ic) /) ! NX(NY(i)) = i 
    
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

    integer im, ic, nx(ic), i, j, ifail
    double precision a(im,ic), xin(im), cin(ic), res(im), cb(ic)
    double precision xmeas(im), xres(im), y(ic,im), z(ic,ic), xd(ic)

    double precision zero
    parameter(zero=0d0)

    integer w(ic)
      
    write(*,*) ' '
    write(*,*) 'start LEAST SQUARES correction with all correctors'
    write(*,*) ' '

    RES = zero

    XMEAS = XIN
    XRES = zero
    Y = TRANSPOSE(A) 
    Z = MATMUL(Y, A) 

    call dinv(ic,Z,ic,w,ifail)
    if ( ifail .ne. 0 ) then
       write(*,*) 'IFAIL from dinv: ',ifail
    endif
    
    XD = MATMUL(Y, XMEAS) 
    CB = MATMUL(Z, XD)
    CB = -CB 
    XRES = MATMUL(A, CB)
 
    CIN = CB
    RES = XRES + XMEAS
    NX(1:ic) = (/ (i, i = 1, ic) /) ! NX(i) = i 

    write(*,*) ' '
    write(*,*) 'end LEAST SQUARES correction with all correctors'
    write(*,*) ' '

    return
end subroutine haveit
      
subroutine svddec_m(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,    &
    &ws,wvec,sortw,sngcut,sngval,im,ic,iflag,sing,dbg)
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

    ! 2013-Dec-19  09:46:02  ghislain: explicit opening of fort.61 for Windows
    if(dbg.gt.0) then
       open(61,file='fort.61')
       
       write(*,*) 'SVD parameters: '
       write(*,*) 'SNGCUT:         ',sngcut
       write(*,*) 'SNGVAL:         ',sngval
    endif

    matu = .TRUE.
    matv = .TRUE.
    iflag = 0

    svdnm = max(ic,im)
    svdmx = im
    svdnx = ic

    SVDMAT = A

    if(dbg.gt.0) then
       write(*,*) 'A0:'
       do  j = 1,im
          write(*,'(16(2X,F7.2))') (svdmat(j,i),i=1,ic)
       enddo
    endif

    call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,matv,vmat,amater,ws)

    if(amater.ne.0) then
        write(*,*) 'end SVD with error code: ',amater
    endif

    if(dbg.gt.0) then
       do  i = 1,ic
          write(*,'(1X,"Corrector: ",I4,"   sing: ",F12.4)') i,wvec(i)
          wmat(i,i) = wvec(i) ! 2015-May-30  19:05:43  ghislain: strange!!!
       enddo
    endif

    call rvord(wvec,sortw,ws,ic)

    if(dbg.gt.0) then
       do  i = 1,ic
          write(*,*) i,sortw(i),wvec(sortw(i))
       enddo
    endif

    do  ii = 1,min(nsing,ic)
       i = sortw(ii)
       
       if(dbg.gt.0) then
          write(61,*) wvec(i)
       endif

       if(abs(wvec(i)).lt.sngval) then
          
          if(dbg.gt.0) then
              do  j = 1,ic
                 write(*,'(A7,2I4,5X,2(F12.6,2X))') 'VMAT: ',i,j,vmat(j,i)
                 ! ghislain: strange we print the transpose here!!
              enddo
           endif

           do  j = 1,ic-1
              do  jj = j+1,ic

                 if(abs(vmat(j,i)).gt.1.0E-4) then
                    rat = abs(vmat(j,i)) + abs(vmat(jj,i))
                    ! rat = abs(vmat(j,i) - vmat(jj,i))
                    ! rat = abs(vmat(j,i) + vmat(jj,i))
                    rat = rat/abs(abs(vmat(j,i)) - abs(vmat(jj,i)))

                    if(dbg.eq.1) then
                       write(62,*) wvec(i),rat
                    endif
        
                    if(rat.gt.sngcut) then

                       if(dbg.gt.0) then
                          write(*,*)  'dependent pair: ',i,j,jj,rat
                          write(65,*) 'dependent pair: ',i,j,jj,rat
                       endif
                       
                       ! Ghislain : was  "if(iflag.lt.(ic*ic*ic)) then"
                       ! triggered a bug on MICADO with ncond=1
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

     if(dbg.gt.0) then
        do  j=1,iflag
           write(66,*) j,sing(1,j),sing(2,j)
        enddo
        ! 2013-Dec-19  09:46:02  ghislain: explicit closing of fort.61 for Windows
        close(61)
     endif
    
    return
end subroutine svddec_m

subroutine svddec_c(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,    &
    &ws,wvec,sortw,sngcut,sngval,im,ic,iflag,sing,dbg)
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

    ! 2013-Dec-19  09:46:02  ghislain: explicit opening of fort.61 for Windows.
    if(dbg.gt.0) then
       open(61,file='fort.61')

        write(*,*) 'SVD parameters: '
        write(*,*) 'SNGCUT:         ',sngcut
        write(*,*) 'SNGVAL:         ',sngval
    endif

    matu = .TRUE.
    matv = .TRUE.
    iflag = 0

    svdnm = max(ic,im)
    svdmx = ic
    svdnx = im

    SVDMAT = transpose(A)
    
    if(dbg.gt.0) then
       write(*,*) 'A0:'
       do  j = 1,ic
          write(*,'(16(2X,F7.2))') (svdmat(j,i),i=1,im)
       enddo
    endif

    call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,matv,vmat,amater,ws)

    if(amater.ne.0) then
        write(*,*) 'end SVD with error code: ',amater
    endif

    if(dbg.gt.0) then
       do  i = 1,im
          write(*,'(1X,"Corrector: ",I4,"   sing: ",F12.4)') i,wvec(i)
          wmat(i,i) = wvec(i) ! 2015-May-30  19:06:33  ghislain: strange !!!
       enddo
    endif

    call rvord(wvec,sortw,ws,im)

    if(dbg.gt.0) then
        do  i = 1,im
            write(*,*) i,sortw(i),wvec(sortw(i))
        enddo
    endif

    do  ii = 1,min(nsing,im)
        i = sortw(ii)

        if(dbg.gt.0) then
            write(61,*) wvec(i)
        endif

        if(abs(wvec(i)).lt.sngval) then
           
           if(dbg.gt.0) then
              do  j = 1,ic
                 write(*,'(A7,2I4,5X,2(F12.6,2X))') 'UMAT: ',i,j,umat(j,i)
                 ! ghislain: why do we print transpose of UMAT silently ? 
              enddo
           endif

           do  j = 1,ic-1
              do  jj = j+1,ic

                 if(abs(umat(j,i)).gt.1.0E-4) then
                    ! rat = abs(umat(j,i) - umat(jj,i))
                    rat = abs(umat(j,i)) +  abs(umat(jj,i))
                    rat = rat/abs(abs(umat(j,i)) - abs(umat(jj,i)))
                    
                    if(dbg.gt.0) then
                       write(62,*) wvec(i),rat
                    endif
                    
                    if(rat.gt.sngcut) then
                       if(dbg.gt.0) then
                          write(*,*)  'dependent pair: ',j,jj,rat
                          write(65,*) 'dependent pair: ',j,jj,rat
                       endif
                        
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

    if(dbg.gt.0) then
       do  j=1,iflag
          write(66,*) j,sing(1,j),sing(2,j)
       enddo
       
       ! 2013-Dec-19  09:46:02  ghislain: explicit closing of fort.61 for Windows
       close(61)
    endif

    return
end subroutine svddec_c

subroutine svdcorr_m(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,   &
    &xin,xc,xout,xa,xb,xpred,ws,wvec,sortw,nx,im,ic,iflag,dbg)
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
    integer sortw(ic)
    integer amater, svdmx, svdnx, svdnm
    logical matu, matv
    integer dbg

    matu = .TRUE.
    matv = .TRUE.
    iflag = 0

    write(*,*) ' '
    write(*,*) 'start SVD correction using ',ic,' correctors'
    write(*,*) ' '

    svdnm = max(ic,im)
    svdmx = im
    svdnx = ic

    do  i = 1,ic
        nx(i) = i
    enddo

    SVDMAT = A

    if(dbg.gt.0) then
       write(*,*) 'A0:'
       do  j = 1,im
          write(*,'(16(2X,F7.2))') (svdmat(j,i),i=1,ic)
       enddo
    endif

    call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,matv,vmat,amater,ws)

    if(amater.ne.0) then
        write(*,*) 'end SVD with error code: ',amater
    endif

    do  i = 1,ic
       if(dbg.gt.0) write(*,'(1X,"Corrector: ",I4,"   sing: ",F12.4)') i,wvec(i)

       wmat(i,i) = wvec(i)

       if(abs(wvec(i)).gt.1.0001) then
          wtmat(i,i) = 1/wvec(i)
       else
          wtmat(i,i) = 0.0
       endif
    enddo

    if(dbg.gt.0) then
       call rvord(wvec,sortw,ws,ic)
       do  i = 1,ic
          write(*,*) i,sortw(i),wvec(sortw(i))
       enddo
    endif

    VTMAT = transpose(VMAT)
    UTMAT = transpose(UMAT)

    if(dbg.gt.0) then
       write(*,*) 'A1:'
       do  j = 1,im
          write(*,'(16(2X,F7.2))') (svdmat(j,i),i=1,ic)
       enddo
       write(*,*) ' '       
       write(*,*) 'Va:'
       do  j = 1,ic
          write(*,'(16(2X,F7.2))') (vmat(j,i),i=1,ic)
       enddo
       write(*,*) 'Vt:'
       do  j = 1,ic
          write(*,'(16(2X,F7.2))') (vtmat(j,i),i=1,ic)
       enddo
       write(*,*) 'W:'
       do  j = 1,im
          write(*,'(16(2X,F7.2))') (wmat(j,i),i=1,ic)
       enddo
       write(*,*) 'Wt:'
       do  j = 1,ic
          write(*,'(16(2X,F7.2))') (wtmat(j,i),i=1,im)
       enddo
       write(*,*) 'U:'
       do  j = 1,im
          write(*,'(16(2X,F7.2))') (umat(j,i),i=1,im)
       enddo
       write(*,*) 'Ut:'
       do  j = 1,im
          write(*,'(16(2X,F7.2))') (utmat(j,i),i=1,im)
       enddo
    endif

    XA(:im) = matmul(UTMAT(:im,:im), XIN(:im))
    XB(:ic) = matmul(WTMAT(:ic,:im), XA(:im))
    XC(:ic) = matmul(VMAT(:ic,:ic), XB(:ic))
    XPRED(:im) = matmul(SVDMAT(:im,:ic), XC(:ic))

    if(dbg.gt.0) then
       write(*,*) "correctors: ", xc
       write(*,*) "monitors:   ", xpred
    endif

    XOUT = XIN - XPRED    
    XC = - XC

    return
end subroutine svdcorr_m

subroutine svdcorr_c(a,svdmat,umat,vmat,wmat,utmat,vtmat,wtmat,   &
    &xin,xc,xout,xa,xb,xpred,ws,wvec,sortw,nx,im,ic,iflag,dbg)
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
    double precision xa(im),xb(im),xpred(im),ws(ic),wvec(ic)
    integer sortw(ic)
    integer amater, svdmx, svdnx, svdnm
    logical matu, matv
    integer dbg
    
    matu = .TRUE.
    matv = .TRUE.
    iflag = 0

    write(*,*) ' '
    write(*,*) 'start SVD correction using ',ic,' correctors'
    write(*,*) ' '

    svdnm = max(ic,im)
    svdmx = ic
    svdnx = im

    do  i = 1,ic
        nx(i) = i
    enddo

    SVDMAT = transpose(A)

    if(dbg.gt.0) then
       write(*,*) 'A0:'
       do  j = 1,ic
          write(*,'(16(2X,F7.2))') (svdmat(j,i),i=1,im)
       enddo
    endif

     call svd(svdnm,svdmx,svdnx,svdmat,wvec,matu,umat,matv,vmat,amater,ws)

     if(amater.ne.0) then
         write(*,*) 'end SVD with error code: ',amater
     endif

     do  i = 1,im
        if(dbg.gt.0) write(*,'(1X,"Corrector: ",I4,"   sing: ",F12.4)') i,wvec(i)

        wmat(i,i) = wvec(i)

        if(abs(wvec(i)).gt.1.0001) then
           wtmat(i,i) = 1/wvec(i)
        else
           wtmat(i,i) = 0.0
        endif
     enddo

     if(dbg.gt.0) then
         call rvord(wvec,sortw,ws,im)
         do  i = 1,im
             write(*,*) i,sortw(i),wvec(sortw(i))
         enddo
     endif

     VTMAT = transpose(VMAT)
     UTMAT = transpose(UMAT)

     if(dbg.gt.0) then
        write(*,*) 'A1:'
        do  j = 1,ic
           write(*,'(16(2X,F7.2))') (svdmat(j,i),i=1,im)
        enddo
        write(*,*) ' '        
        write(*,*) 'Va:'
        do  j = 1,im
           write(*,'(16(2X,F7.2))') (vmat(j,i),i=1,im)
        enddo
        write(*,*) 'Vt:'
        do  j = 1,im
           write(*,'(16(2X,F7.2))') (vtmat(j,i),i=1,im)
        enddo
        write(*,*) 'W:'
        do  j = 1,ic
           write(*,'(16(2X,F7.2))') (wmat(j,i),i=1,im)
        enddo
        write(*,*) 'Wt:'
        do  j = 1,im
           write(*,'(16(2X,F7.2))') (wtmat(j,i),i=1,ic)
        enddo
        write(*,*) 'U:'
        do  j = 1,ic
           write(*,'(16(2X,F7.2))') (umat(j,i),i=1,ic)
        enddo
        write(*,*) 'Ut:'
        do  j = 1,ic
           write(*,'(16(2X,F7.2))') (utmat(j,i),i=1,ic)
        enddo
     endif
      
     XA(:im) = matmul(VTMAT(:im,:im), XIN(:im))
     XB(:im) = matmul(WTMAT(:im,:im), XA(:im))
     XC(:ic) = matmul(UMAT(:ic,:im),XB(:im))
     XPRED(:im) = matmul(A(:im,:ic), XC(:ic))

     if(dbg.gt.0) then
        write(*,*) "correctors: ", xc
        write(*,*) "monitors:   ", xpred
     endif

     XOUT = XIN - XPRED
     XC = -XC

     return
 end subroutine svdcorr_c

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
     integer m,n,ipiv(n),prtlev,iter,ij1,k2,k,i,kpiv,k3,j,ip,j1,kk,ki, &
         &iii,kkk
     real a(m,n),b(m),x(n),r(m),rms,rho(3*n),ptop(n),rmss(n),xrms(n),  &
         &xptp(n),xiter(n),xxcal,ptp,g,h,sig,beta,piv,pivt,rm,pt

     real rzero, reps7
     parameter(rzero=0e0,reps7=1e-7)

     character(4) :: units
     character(16) :: conm(n)
     integer      ifail
     real calrms

     integer k1, kn, kl, lv

     ifail = 0
     units = 'mrad'
     ptp = rzero

     rm = calrms(b,m)
     pt = MAXVAL(B)-MINVAL(B)

     ! --- calculate first pivot
     !==========================

     RHO = rzero

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
        
        if (kpiv.eq.k) go to 8

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
8       continue
        call htul(a,m,n,k,sig,beta)
        
        ! --- on garde SIGMA dans RHO(N+K)
        j=n+k
        rho(j)=-sig
        ip=ipiv(kpiv)
        ipiv(kpiv)=ipiv(k)
        ipiv(k)=ip
        if(k.eq.n) go to 13

        ! --- transformation de A dans HTAL
        call htal(a,m,n,k,beta)
        
        ! --- transformation de B dans HTBL
13      continue
        call htbl(a,b,m,n,k,beta)

        ! --- recherche du pivot (K+1)
        !=============================

        rho(k)=sqrt(piv)
        if(k.eq.n) go to 11

        piv=rzero
        kpiv = k + 1
        j1 = kpiv
        k2=n + j1

        do j=j1,n
           h=rho(j)-(a(k,j))*(a(k,j))

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
           endif

           rho(j)=h
           g=rho(k2)-(a(k,j))*(b(k))
           rho(k2) = g
           if(h.ne.rzero) then
              pivt = g*g/h
           else
              pivt = rzero
           endif

           if(pivt.lt.piv)go to 18

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

        rmss(k) = calrms(r,m)
        ptop(k)=MAXVAL(R)-MINVAL(R)

        if(k.lt.n) then
           xiter(k+1) = k
           xrms(k+1)  = rmss(k)
           xptp(k+1)  = ptop(k)
        endif
        
        ! --- write intermediate results to 61 files
        if(k.eq.1) then
           if (prtlev .ge. 2) then
              write(61,52)
52            FORMAT(/' ***********    start MICADO    ***********'/)
              write(61,54) units
54            FORMAT(' iter',5X,'corrector',13X,A4,6X,'mrad', 5X,"  rms",10X," ptop",/)
              write(61,'(4x,"0",42x,f12.8,f15.8)')rm,pt
           endif
        endif

        if (prtlev .ge. 2) then
           write(61,'(/,1x,72("-"),/,1x,i4,3x,38(" "),1x,f12.8,f15.8)') k,rmss(k),ptop(k)
        endif

        do kkk = 1,k
           xxcal=x(kkk)
           if (prtlev .ge. 2) then
              write(61,'(I3,1X,A16,9x,F8.4,2X,F8.4)') kkk,conm(ipiv(kkk)),x(kkk),xxcal
           endif
        enddo
        
        if (prtlev .ge. 2) then
           write(61,58) k
58         format(/,' residual orbit after iteration ',i4,':')
           write(61,'(1x,8f9.3)')(r(kkk),kkk=1,m)
        endif
        
        if(k.eq.iter) then
           if (prtlev .ge. 2) then
              write(61,53)
53            format(/' ***********    end   MICADO    ***********'/)
           endif
        endif

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
           sig=sig+a(i,k)*a(i,k)
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

  function calrms(r,m)
       implicit none
       !*********************************************************************
       !     Function CALRMS to calculate rms of double precision values    *
       !     in array R of size m                                           *
       !                                                                    *
       !     Authors:     GJR adapted from subroutine calrms by WFH         * 
       !     Date:        2014-Feb-28                                       *
       !                                                                    *
       !*********************************************************************

       real calrms
       integer m,i
       real r(m), rzero
       parameter(rzero=0e0)

       calrms = rzero
      
       do i=1,m
           calrms = calrms + (r(i)*r(i))
       enddo
      
       calrms = sqrt(calrms / float(m))       
       return
   end function calrms

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

1      if(t3.ge.t1) go to 2
       !     (PIVOT IS A11)
       temp=a(1,1)
       det=c22*c33-c23*c32
       go to 3

       !     (PIVOT IS A31)
2      temp=a(3,1)
       det=c23*c12-c22*c13
       !
       !     SET ELEMENTS OF INVERSE IN A.
3      if(det.eq.zero) go to 8
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
4      if(n.lt.2) go to 5
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
5      if(a(1,1).eq.zero) go to 8
       a(1,1)=1d0/a(1,1)
       return
       !
       !  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
       !
6      call dfact(n,a,idim,r,ifail,det,jfail)
       if(ifail.ne.0) return
       call dfinv(n,a,idim,r)
       return
       !
       !  ERROR EXITS.
       !
7      ifail=+1
       return
       !
8      ifail=-1
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

       if(idim .lt. n  .or.  n .lt. 1)  then
1001      FORMAT(7X," PARAMETER ERROR IN SUBROUTINE ", A6, &
               &" ... (N.LT.1 OR IDIM.LT.N).", &
               &5X,"N =", I4, 5X,"IDIM =", I4,".")
          write(*,1001) "DFACT", n, idim          
          return
       endif

110    ifail  =  normal
       jfail  =  jrange
       nxch   =  0
       det    =  one
       do    j  =  1, n
120        k  =  j
           p  =  abs(sngl(a(j,j)))
           if(j .eq. n)  goto 122
           jp1  =  j+1
           do    i  =  jp1, n
               q  =  abs(sngl(a(i,j)))
               if(q .le. p)  goto 121
               k  =  i
               p  =  q
121        continue
           enddo
           if(k .ne. j)  goto 123
122        if(p .gt. zero)  goto 130
           det    =  zero
           ifail  =  imposs
           jfail  =  jrange
           return
123        do    l  =  1, n
               tf      =  a(j,l)
               a(j,l)  =  a(k,l)
               a(k,l)  =  tf
           enddo
           nxch      =  nxch + 1
           ir(nxch)  =  j*2**12 + k
130        det     =  det * a(j,j)
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
142            a(j,k)    =  -s11 * a(j,j)
               a(k,j+1)  =  -(a(j,j+1)*a(k,j)+s12)
           enddo
144    continue
       enddo
150    if(mod(nxch,2) .ne. 0)  det  =  -det
       if(jfail .ne. jrange)   det  =  zero
       ir(n)  =  nxch
       return
   end subroutine dfact

   subroutine dfeqn(n,a,idim,ir,k,b)
     implicit none
     integer idim,i,j,k,n,m,nxch,ij,l,im1,nm1,nmi,nmjp1,ir(*)
     double precision a(idim,*),b(idim,*),te,s21, s22
     
     if (idim .lt. n  .or.  n .lt. 1  .or.  k .lt. 1)  then
1002    FORMAT(7X," PARAMETER ERROR IN SUBROUTINE ", A6, &
             &" ... (N.LT.1 OR IDIM.LT.N).", &
             &5X,"N =", I4, 5X,"IDIM =", I4, 5X,"K =", I4,".")
        write(*,1002) "DFEQN", n, idim, k
        return
     endif
     
     nxch  =  ir(n)
     
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
     
220  do    l  =  1, k
        b(1,l)  =  a(1,1)*b(1,l)
     enddo
     
     if(n .eq. 1)  return
     
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
     
     return
end subroutine dfeqn

subroutine dfinv(n,a,idim,ir)
    implicit none
    integer idim,i,j,k,n,m,nxch,ij,nm1,nmi,im2,ir(*)
    double precision a(idim,*),ti,s31,s32,s33,s34,zero
    parameter(zero=0d0)
    
    if (idim < n  .or.  n < 1 )  then   
1001   FORMAT(7X," PARAMETER ERROR IN SUBROUTINE ", A6, &
            &" ... (N.LT.1 OR IDIM.LT.N).", &
            &5X,"N =", I4, 5X,"IDIM =", I4,".")
       write(*,1001) "DFINV", n, idim
       return
    endif

    if(n .eq. 1)  return
    
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
330 nm1  =  n-1
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

subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
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

        if (scale .eq. 0.0d0) go to 210

        do k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
        enddo

        f = u(i,i)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,i) = f - g
        if (i .eq. n) go to 190

        do j = l, n
            s = 0.0d0

            do k = i, m
                s = s + u(k,i) * u(k,j)
            enddo

            f = s / h

            do k = i, m
                u(k,j) = u(k,j) + f * u(k,i)
            enddo
        enddo

        190   do k = i, m
            u(k,i) = scale * u(k,i)
        enddo

210     w(i) = scale * g
        g = 0.0d0
        s = 0.0d0
        scale = 0.0d0
        if (i .gt. m .or. i .eq. n) go to 290

        do k = l, n
            scale = scale + dabs(u(i,k))
        enddo

        if (scale .eq. 0.0d0) go to 290

        do k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
        enddo

        f = u(i,l)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,l) = f - g

        do k = l, n
            rv1(k) = u(i,k) / h
        enddo

        if (i .eq. m) go to 270

        do j = l, m
            s = 0.0d0

            do k = l, n
                s = s + u(j,k) * u(i,k)
            enddo

            do k = l, n
                u(j,k) = u(j,k) + s * rv1(k)
            enddo
        enddo

        270   do k = l, n
            u(i,k) = scale * u(i,k)
        enddo

290     x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
    enddo
    !     .......... accumulation of right-hand transformations ..........
    if (.not. matv) go to 410
    !     .......... for i=n step -1 until 1 do -- ..........
    do ii = 1, n
        i = n + 1 - ii
        if (i .eq. n) go to 390
        if (g .eq. 0.0d0) go to 360

        do j = l, n
            !     .......... double division avoids possible underflow ..........
            v(j,i) = (u(i,j) / u(i,l)) / g
        enddo

        do j = l, n
            s = 0.0d0

            do k = l, n
                s = s + u(i,k) * v(k,j)
            enddo

            do k = l, n
                v(k,j) = v(k,j) + s * v(k,i)
            enddo
        enddo

360     do j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
        enddo

390     v(i,i) = 1.0d0
        g = rv1(i)
        l = i
    enddo
    !     .......... accumulation of left-hand transformations ..........
410 if (.not. matu) go to 510
    !     ..........for i=min(m,n) step -1 until 1 do -- ..........
    mn = n
    if (m .lt. n) mn = m

    do ii = 1, mn
        i = mn + 1 - ii
        l = i + 1
        g = w(i)
        if (i .eq. n) go to 430

        do j = l, n
            u(i,j) = 0.0d0
        enddo

430     if (g .eq. 0.0d0) go to 475
        if (i .eq. mn) go to 460

        do j = l, n
            s = 0.0d0

            do k = l, m
                s = s + u(k,i) * u(k,j)
            enddo
            !     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g

            do k = i, m
                u(k,j) = u(k,j) + f * u(k,i)
            enddo
        enddo

460     do j = i, m
            u(j,i) = u(j,i) / g
        enddo

        go to 490

475     do j = i, m
            u(j,i) = 0.0d0
        enddo

490     u(i,i) = u(i,i) + 1.0d0
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
520     do ll = 1, k
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
540     c = 0.0d0
        s = 1.0d0

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

            do j = 1, m
                y = u(j,l1)
                z = u(j,i)
                u(j,l1) = y * c + z * s
                u(j,i) = -y * s + z * c
            enddo

560     continue
        enddo

        !     .......... test for convergence ..........
565     z = w(k)
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

            do j = 1, n
                x = v(j,i1)
                z = v(j,i)
                v(j,i1) = x * c + z * s
                v(j,i) = -x * s + z * c
            enddo

575         z = pythag(f,h)
            w(i1) = z
            !     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
580         f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600

            do j = 1, m
                y = u(j,i1)
                z = u(j,i)
                u(j,i1) = y * c + z * s
                u(j,i) = -y * s + z * c
            enddo

600     continue
        enddo

        rv1(l) = 0.0d0
        rv1(k) = f
        w(k) = x
        go to 520
        !     .......... convergence ..........
650     if (z .ge. 0.0d0) go to 700
        !     .......... w(k) is made non-negative ..........
        w(k) = -z
        if (.not. matv) go to 700

        do j = 1, n
            v(j,k) = -v(j,k)
        enddo

700 continue
    enddo

    go to 1001
    !     .......... set error -- no convergence to a
    !                singular value after 30 iterations ..........
1000 ierr = k
1001 return
end subroutine svd


double precision function pythag(a,b)
    implicit  none
    double precision a,b
    !
    !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
    !
    double precision p,r,s,t,u
    p = dmax1(dabs(a),dabs(b))
    if (p .eq. 0.0d0) go to 20

    r = (dmin1(dabs(a),dabs(b))/p)**2

10  continue
    t = 4.0d0 + r
    if (t .eq. 4.0d0) go to 20
    s = r/t
    u = 1.0d0 + 2.0d0*s
    p = u*p
    r = (s/u)**2 * r
    go to 10

20  pythag = p

   return
   end function pythag

   subroutine rvord(inv,outv,ws,n)
     ! 2015-Apr-28  15:45:16  ghislain: analysis
     ! subroutine to sort the indexes of elements stored in input vector INV 
     ! in reverse order in output vector OUTV.
     ! WS is a workspace vector, N is the dimension of the vectors.
     ! Supposition that INV contains positive numbers only!
       implicit   none
       integer    n
       double precision   inv(n), ws(n)
       integer  i,j,jmax
       integer  outv(n)

       WS = INV

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

subroutine abend
    implicit none
    write(*,*) 'Abnormal end ...'
    stop 888
end subroutine abend
