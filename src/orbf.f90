subroutine micit(a, conm, xin, cin, res, nx, rms, im, ic, iter, ny, ax, cinx, &
                 xinx, resx, rho, ptop, rmss, xrms, xptp, xiter, ifail)
    use math_constfi, only : zero
    use io_units, only : orbfout
    implicit none
    ! ****************************************************
    !                                                    *
    !    Driving routine for MICADO correction           *
    !                                                    *
    !     Author: WFH  05.02.02                          *
    !                                                    *
    ! ****************************************************
    ! RMS = value of tolerance for correction
    integer :: im, ic, iter, ifail
    double precision :: a(im,ic), xin(im), cin(ic), res(im)
    character(len=16) :: conm(ic)
    integer :: nx(ic)
    double precision :: ax(im,ic), cinx(ic), xinx(im), resx(im), rho(3*ic), ptop(ic)
    double precision :: rmss(ic), xrms(ic), xptp(ic), xiter(ic)

    integer :: n, i, j
    integer :: ny(ic)
    double precision :: rms, pt, rm
 
    write (*,'(/A, I5, A/)') 'start MICADO correction with ',iter,' correctors'

    ! translate corrector names to fortran strings
    do  j = 1,ic
       call f_ctof(n, conm(j), 16)
    enddo

    AX(:im,:ic) = A(:im,:ic)
    CINX(:ic) = zero

    NY(1:ic) = (/ (i, i = 1, ic) /) ! NY(i) = i 

    XINX(:im) = XIN(:im)
    RESX(:im) = zero

    rm = sqrt(dot_product(XINX,XINX)/float(im))
    if(rm.le.rms) then
       write(*,*) '++++++ WARNING: RMS already smaller than desired '
       write(*,*) '++++++ WARNING: no correction is done            '
       rms = rm
       iter = 0
       ifail = -2
    else
       open(unit=orbfout,file='fort.61')
       call htls(ax, conm, xinx, im, ic, cinx, ny, resx, rms, 3, & 
                 iter, rho, ptop, rmss, xrms, xptp, xiter, ifail)
       close(unit=orbfout)
    endif

    CIN(:ic) = CINX(:ic)
    RES(:im) = RESX(:im)
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
    use math_constfi, only : zero
    implicit none

    integer :: im, ic
    integer :: nx(ic)
    double precision :: a(im,ic), xin(im), cin(ic), res(im), cb(ic)
    double precision :: xmeas(im), xres(im), y(ic,im), z(ic,ic), xd(ic)

    ! Note XD is no longer used

    integer :: i, j, ifail
    integer :: w(ic)
          
    write (*,'(/A/)') 'start LEAST SQUARES correction with all correctors'

    RES(:im) = zero

    XMEAS(:im) = XIN(:im)
    XRES(:im) = zero
 
    Y(:ic,:im) = TRANSPOSE(A(:im,:ic)) 
    Z(:ic,:ic) = MATMUL(Y(:ic,:im), A(:im,:ic)) 

    call dinv(ic,Z,ic,w,ifail)
    if ( ifail .ne. 0 ) & !write(*,*) 'IFAIL from dinv: ',ifail 
         call fort_warn('LSQ in HAVEIT: ','DINV returned failure code '//char(ifail+ichar('0')))

    ! pseudoinverse of A is (At * A)^-1 * At
    ! and solution for correctors is
    !CB = - ( Z * Y ) * XMEAS = - (At A)-1 * At * XMEAS
    
    CB = -matmul( matmul(Z,Y) , XMEAS)
    !XD(:ic) = MATMUL(Y(:ic,:im), XMEAS(:im)) 
    !CB(:ic) = MATMUL(Z(:ic,:ic), XD(:ic))
    !CB(:ic) = -CB(:ic) 
    
    ! and resulting orbit change is A*CB
    XRES(:im) = MATMUL(A(:im,:ic), CB(:ic))
 
    CIN(:ic) = CB(:ic) 
    RES(:im) = XRES(:im) + XMEAS(:im) ! residual orbit
    NX(1:ic) = (/ (i, i = 1, ic) /) ! NX(i) = i 

    write (*,'(/A/)') 'end LEAST SQUARES correction with all correctors'

    return
end subroutine haveit

subroutine svddec(a, svdmat, umat, vmat, ws, wvec, sortw, &
                  sngcut, sngval, im, ic, iflag, sing, dbg)
    use math_constfi, only : zero
    implicit none
    ! ****************************************************
    !                                                    *
    !    Performs SVD analysis for any response matrix   *      
    !                                                    *
    !     Author: GJR, 2015-07-01                        *
    !             based on work by WFH  12.09.02         *
    !                                                    *
    ! ****************************************************
    integer :: im, ic
    double precision :: a(im,ic), svdmat(im,ic)
    double precision :: umat(im,ic)
    double precision :: vmat(ic,ic)
    double precision :: wmat(ic,ic)
    double precision :: ws(ic), wvec(ic)
    double precision :: sngval, sngcut
    integer :: sortw(ic), iflag, sing(2,ic)
    integer :: dbg

    integer :: i, j, jj, ii, errflag
    integer, parameter :: nsing=5
    double precision :: rat

    SVDMAT(:im,:ic) = A(:im,:ic)

    call prepsvd(im, ic, svdmat, wvec, umat, vmat, errflag, ws)
    if (errflag .ne. 0) write(*,*) 'end SVD with error code: ',errflag

    call rvord(wvec,sortw,ws,ic)

    iflag = 0
    do  ii = 1, min(nsing,ic)
       i = sortw(ii)
       if ( abs(wvec(i)) .lt. sngval) then
           do  j = 1, ic-1
              do  jj = j+1, ic

                 if( abs(vmat(j,i)) .gt. 1.0d-4) then
                    rat = abs(vmat(j,i)) + abs(vmat(jj,i))
                    rat = rat/abs(abs(vmat(j,i)) - abs(vmat(jj,i)))

                    if (rat .gt. sngcut) then
                       if (iflag .lt. ic) then
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

end subroutine svddec

subroutine svdcorr(a, svdmat, umat, vmat, wmat, utmat, vtmat, wtmat, &
     xin, xc, xout, xpred, ws, wvec, sortw, nx, im, ic, iflag, dbg)
    use math_constfi, only : zero
    implicit none
    ! ******************************************************
    !                                                      *
    !  Performs SVD and correction for any response matrix *      
    !                                                      *
    !     Author: GJR, 2015-07-01                          *
    !             based on work by WFH  12.09.02           *
    !                                                      *
    ! ******************************************************
    integer :: im, ic 
    double precision :: a(im,ic), svdmat(im,ic)
    double precision :: umat(im,ic), utmat(ic,im) 
    double precision :: vmat(ic,ic), vtmat(ic,ic)     
    double precision :: wmat(ic,ic), wtmat(ic,ic)
    double precision :: xin(im), xout(im), xpred(im)
    double precision :: xc(ic), ws(ic), wvec(ic)
    integer :: sortw(ic), nx(ic), iflag, dbg 

    integer :: i, j
    double precision :: ainv(ic,im)

    write (*,'(/A,I5,A/)') 'start SVD correction using ',ic,' correctors'

    SVDMAT = A

    call prepsvd(im, ic, svdmat, wvec, umat, vmat, iflag, ws)

    WMAT = zero ; WTMAT = zero
    do  i = 1,ic
       wmat(i,i) = wvec(i) 
       if(abs(wvec(i)).gt.1.0001) wtmat(i,i) = 1/wvec(i) 
    enddo

    VTMAT = transpose(VMAT)
    UTMAT = transpose(UMAT)

    AINV = matmul(VMAT, matmul(WTMAT,UTMAT))

    XC = matmul(AINV,XIN)
    XPRED = matmul(SVDMAT, XC)

    XOUT = XIN - XPRED
    XC = - XC

    ! NX(i) = i in FORTRAN, used only in C such that nx(0) = 1
    NX(1:ic) = (/ (i, i = 1, ic) /) 

    if(dbg.gt.0) then
       write(*,*) 'SVDMAT:'; call prdmat(svdmat,im,ic); write(*,*) ' '       
       write(*,*) 'VMAT:';   call prdmat(vmat,ic,ic); write(*,*) ' '       
       write(*,*) 'VTMAT:';  call prdmat(vtmat,ic,ic); write(*,*) ' ' 
       write(*,*) 'WMAT:';   call prdmat(wmat,ic,ic); write(*,*) ' ' 
       write(*,*) 'WTMAT:';  call prdmat(wtmat,ic,ic); write(*,*) ' ' 
       write(*,*) 'UMAT:';   call prdmat(umat,im,ic); write(*,*) ' ' 
       write(*,*) 'UTMAT:';  call prdmat(utmat,ic,im); write(*,*) ' ' 
       do i = 1, ic
          write(*,'(1X,"Corrector: ",I4,"   sing: ",F12.4)') i, wvec(i)
       enddo
       call rvord(wvec,sortw,ws,ic)
       write (*,*) 'Sorted vector of singular values: i sortw(i) wvec(sortw(i)) '
       do  i = 1, ic
          write(*,*) i, sortw(i), wvec(sortw(i))
       enddo
       write(*,*) ' ' ; write(*,*) "correctors: ", xc
       write(*,*) ' ' ; write(*,*) "monitors:   ", xpred
    endif

end subroutine svdcorr

subroutine prepsvd(im, ic, svdmat, wvec, umat, vmat, iflag, ws)
    use math_constfi, only : zero
    implicit none
    ! ******************************************************
    !                                                      *
    !   Prepare SVD and correction for any response matrix * 
    !   interface routine to the SVD subroutine that takes *
    !   extended matrices such that the number of lines is *
    !   always the largest of im and ic
    !                                                      *
    !     Author: GJR, 2015-07-01                          *
    !                                                      *
    ! ******************************************************
    integer :: im, ic, iflag 
    double precision :: svdmat(im,ic), umat(im,ic), vmat(ic,ic)
    double precision :: ws(ic), wvec(ic)

    double precision :: vmat2(im,ic) 
    double precision :: umat2(ic,ic), svdmat2(ic,ic)

    if (im .ge. ic) then
       VMAT2 = zero ; UMAT = zero
       call svd(im, im, ic, svdmat, wvec, .true., umat, .true., vmat2, iflag, ws)
       VMAT(:ic,:ic) = VMAT2(:ic,:ic)

    elseif (im .lt. ic) then
       SVDMAT2 = zero ;
       SVDMAT2(:im,:ic) = SVDMAT(:im,:ic)
       UMAT2 = zero ;  VMAT = zero
       call svd(ic, im, ic, svdmat2, wvec, .true., umat2, .true., vmat, iflag, ws)
       UMAT(:im,:ic) = UMAT2(:im,:ic)
       SVDMAT(:im,:ic) = SVDMAT2(:im,:ic)
    endif

    return
end subroutine prepsvd

subroutine htls(a, conm, b, m, n, x, ipiv, r, rms, prtlev, iter, rho, ptop, &
                rmss, xrms, xptp, xiter, ifail)
     use math_constfi, only : zero, one
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

     integer, intent(IN)  :: m, n, prtlev
     integer, intent(OUT) :: iter, ifail
     double precision :: a(m,n), b(m), x(n), r(m)
     character(len=16) :: conm(n)     
     integer :: ipiv(n)
     double precision :: rms 
     double precision :: rho(3*n), ptop(n), rmss(n), xrms(n), xptp(n), xiter(n)

     integer :: i, j, j1, j2, k, k2, k3, ij1, kpiv, kk, ki, iii, kkk
     integer :: k1, kn, kl
     double precision :: ptp, g, h, h2, h3, sig, beta, piv, pivt, rm, pt, temp

     double precision, parameter :: reps7=1.d-7
     character(len=4) :: units='mrad'

     ifail = 0
     ptp = zero

     rm = sqrt(dot_product(B(:m),B(:m))/float(m))
     pt = MAXVAL(B(:m))-MINVAL(B(:m))

     ! --- calculate first pivot
     !==========================

     RHO(:3*n) = zero

     k2 = n + 1
     piv = zero
      
     do k = 1, n ! initialisation loop over correctors
        h = dot_product(A(1:m,k),A(1:m,k)) ; rho(k)  = h
        g = dot_product(A(1:m,k),B(1:m))   ; rho(k2) = g
        if (h .ne. zero) then
           pivt = g*g/h
        else
           pivt = zero
        endif
        if(pivt-piv .gt. zero) then
           piv = pivt
           kpiv = k
        endif
        k2 = k2 + 1
     enddo

     ! --- boucle pour chaque iteration
     do k = 1, iter

        if (kpiv .ne. k) then
           ! --- on echange les K et KPIV si KPIV plus grand que K
           call swapreal(rho(k),rho(kpiv))
           k2 = n + k;   k3 = n + kpiv
           call swapreal(rho(k2), rho(k3))
           do i = 1, m ! swap a(i,k) and a(i,kpiv) for all monitors
              call swapreal(a(i,k),a(i,kpiv))
           enddo
        endif

        ! --- calcul de beta,sigma et uk
        sig = sqrt(dot_product(A(k:m,k),A(k:m,k)))
        h2 = a(k,k)
        sig = sign(sig,h2)
        beta = h2 + sig
        a(k,k) = beta
        beta = one/(sig*beta)
        
        ! --- on garde SIGMA dans RHO(N+K)
        j = n + k
        rho(j) = -sig

        !--- swap ipiv(k) and ipiv(kpiv)
        call swapint(ipiv(kpiv),ipiv(k))

        if (k .ne. n) then
           ! --- transformation de A
           do i = 1, n-k
              h = beta * dot_product(A(k:m,k),A(k:m,k+i))         
              A(k:m,k+i) = A(k:m,k+i) - A(k:m,k)*h
           enddo
        endif
        ! --- transformation de B dans HTBL
        h3 = beta * dot_product(A(k:m,k),B(k:m))
        B(k:m) = B(k:m) - A(k:m,k)*h3 
       
        ! --- recherche du pivot (K+1)
        !=============================
        rho(k) = sqrt(piv)

        if (k .eq. n) go to 11

        piv = zero
        kpiv = k + 1
        j1 = kpiv
        k2 = n + j1

        do j = j1, n
           h = rho(j) - a(k,j)*a(k,j)

           if(abs(h) .lt. reps7) then
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

           rho(j) = h
           g = rho(k2) - a(k,j)*b(k)
           rho(k2) = g
           if (h .ne. zero) then
              pivt = g*g/h
           else
              pivt = zero
           endif

           if (pivt .ge. piv) then
              kpiv=j
              piv=pivt
           endif

           k2 = k2 + 1
        enddo

        ! --- calcul des x
11      x(k) = b(k)/rho(n+k)

        if (k .ne. 1) then
           do i = 2,k
              kk = k - i + 1
              x(kk) = b(kk)
              ki = kk + 1
              do j = ki, k
                 x(kk) = x(kk) - a(kk,j)*x(j)
              enddo
              x(kk) = x(kk) / rho(n+kk)
           enddo
        endif
        
        ! --- save residual orbit and inverse sign of corrections (convention!)
        R(:m) = B(:m)
        X(:k) = -X(:k)
        
        ! --- calcul du vecteur residuel
        !===============================
        !     transform orbit R back to "normal space"
        R(1:k) = zero
        do i = 1, k
           kl = k - i + 1
           kn = k - i + 1 + n           
           beta = -one / (rho(kn)*a(kl,kl))           
           h = beta * dot_product(A(kl:m,kl),R(kl:m))           
           R(kl:m) = R(kl:m) - A(kl:m,kl)*h            
        enddo

        rmss(k) = sqrt(dot_product(R,R)/float(m))    
        ptop(k) = MAXVAL(R) - MINVAL(R)

        if (k .lt. n) then
           xiter(k+1) = k
           xrms(k+1)  = rmss(k)
           xptp(k+1)  = ptop(k)
        endif

        ! --- write intermediate results to fort.61 file
        if (prtlev .ge. 2) then
           if (k .eq. 1) then
              write(61,'(/" ***********    start MICADO    ***********"/)')
              write(61,'(" iter",5X,"corrector",13X,A4,6X,"mrad",7X,"rms",11X,"ptop",/)') units
              write(61,'(4x,"0",42x,f12.8,f15.8)') rm, pt
           endif

           write(61,'(/,1x,72("-"),/)')
           write(61,'(1x,i4,3x,38(" "),1x,f12.8,f15.8)') k,rmss(k),ptop(k)

           write(61,'(I3,1X,A16,9X,F18.14)') (i, conm(ipiv(i)), x(i), i=1,k)
        
           write(61,'(/," residual orbit after iteration ",i4,":")') k
           write(61,'(1x,8f9.3)') (r(i),i=1,m)
        
           if (k .eq. iter) &
                write(61,'(/" ***********    end   MICADO    ***********"/)')

        endif

        if (ptop(k).le.ptp .or. rmss(k).le.rms ) then
           ! --- correction is already good enough:
           !=======================================
           ptp = ptop(k)
           rms = rmss(k)
           iter = k
           return
        endif

     enddo

     return     
end subroutine htls

subroutine dinv(n,a,idim,r,ifail)
       use math_constfi, only : zero, one
       implicit none
       !******************************************************************
       !
       !     REPLACES A BY ITS INVERSE.
       !
       !     (PARAMETERS AS FOR DEQINV.)
       !
       !     CALLS ... DFACT, DFINV, F010PR.
       !
       !******************************************************************
       integer :: n, idim, ifail
       double precision :: a(idim,n)
       integer :: r(n)

       integer :: jfail, t1, t2, t3
       double precision :: det, temp, s
       double precision :: c11, c12, c13, c21, c22, c23, c31, c32, c33
       
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
       c11 = a(2,2)*a(3,3) - a(2,3)*a(3,2)
       c12 = a(2,3)*a(3,1) - a(2,1)*a(3,3)
       c13 = a(2,1)*a(3,2) - a(2,2)*a(3,1)
       c21 = a(3,2)*a(1,3) - a(3,3)*a(1,2)
       c22 = a(3,3)*a(1,1) - a(3,1)*a(1,3)
       c23 = a(3,1)*a(1,2) - a(3,2)*a(1,1)
       c31 = a(1,2)*a(2,3) - a(1,3)*a(2,2)
       c32 = a(1,3)*a(2,1) - a(1,1)*a(2,3)
       c33 = a(1,1)*a(2,2) - a(1,2)*a(2,1)
       t1 = abs(int(a(1,1)))
       t2 = abs(int(a(2,1)))
       t3 = abs(int(a(3,1)))
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
       s=one/det
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
       a(1,1)=one/a(1,1)
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

end subroutine dinv

subroutine dfact(n,a,idim,ir,ifail,det,jfail)
       use math_constfi, only : zero, one
       implicit none

       integer :: n, idim, ir(*), ifail, jfail 
       double precision :: a(idim,*), det

       integer ::  i, j, k, l, nxch, jp1, jm1
       double precision :: tf, s11, s12, p, q, t

       integer, parameter :: normal=0, imposs=-1, jrange=0, jover=1, junder=-1
       double precision, parameter :: g1=1e-19, g2=1e19       

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
120       k  =  j
          p  =  abs(a(j,j))
          if (j .eq. n)  goto 122
           jp1  =  j+1
           do    i  =  jp1, n
               q  =  abs(a(i,j))
               if (q .le. p)  goto 121
               k  =  i
               p  =  q
121        continue
           enddo
           if (k .ne. j)  goto 123
122        if (p .gt. zero)  goto 130
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
           t  =  abs(det)
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
     integer :: n, idim, ir(*), k
     double precision :: a(idim,*), b(idim,*)

     integer :: i, j, m, nxch, ij, l, im1, nm1, nmi, nmjp1
     double precision :: te, s21, s22
     
     if (idim .lt. n  .or.  n .lt. 1  .or.  k .lt. 1)  then
1002    FORMAT(7X," PARAMETER ERROR IN SUBROUTINE ", A6, &
             &" ... (N.LT.1 OR IDIM.LT.N).", &
             &5X,"N =", I4, 5X,"IDIM =", I4, 5X,"K =", I4,".")
        write(*,1002) "DFEQN", n, idim, k
        return
     endif
     
     nxch  =  ir(n)
     
     if (nxch .eq. 0)  goto 220
     
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
     
     if (n .eq. 1)  return
     
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
    use math_constfi, only : zero
    implicit none

    integer :: n, idim, ir(*)
    double precision :: a(idim,*)
 
    integer :: i, j, k, m, nxch, ij, nm1, nmi, im2
    double precision :: ti, s31, s32, s33, s34
        
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
    integer :: nm, m, n
    double precision :: a(nm,n), w(n), u(nm,n), v(nm,n), rv1(n)
    logical :: matu, matv

    integer :: i, j, k, l, ii, i1, kk, k1, ll, l1, mn, its, ierr
    double precision :: c, f, g, h, s, x, y, z, tst1, tst2, scale
    double precision :: pythag

    ierr = 0

    U(1:m,1:n) = A(1:m,1:n)
    !     .......... householder reduction to bidiagonal form ..........
    g = 0.0d0
    scale = 0.0d0
    x = 0.0d0

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

        U(i:m,i) = U(i:m,i) / scale
        s = dot_product(U(i:m,i),U(i:m,i))        

        f = u(i,i)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,i) = f - g
        if (i .eq. n) go to 190

        do j = l, n
           s = dot_product(U(i:m,i),U(i:m,j))           
           f = s / h          
           U(i:m,j) = U(i:m,j) + f*U(i:m,i)
        enddo

190     continue
        U(i:m,i) = scale * U(i:m,i)

210     w(i) = scale * g
        g = 0.0d0
        s = 0.0d0
        scale = 0.0d0
        if (i .gt. m .or. i .eq. n) go to 290

        do k = l, n
            scale = scale + dabs(u(i,k))
        enddo

        if (scale .eq. 0.0d0) go to 290

        U(i,l:n) = U(i,l:n) / scale
        s = dot_product(U(i,l:n),U(i,l:n))

        f = u(i,l)
        g = -dsign(dsqrt(s),f)
        h = f * g - s
        u(i,l) = f - g

        RV1(l:n) = U(i,l:n) / h

        if (i .eq. m) go to 270

        do j = l, m
           s = dot_product(U(j,l:n),U(i,l:n))
           U(j,l:n) = U(j,l:n) + s * RV1(l:n)
        enddo

270     continue
        U(i,l:n) = scale * U(i,l:n)

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
           s = dot_product(U(i,l:n), V(l:n,j))
           V(l:n,j) = V(l:n,j) + s*V(l:n,i)
        enddo

360     continue
        V(i,l:n) = 0.0d0
        V(l:n,i) = 0.d0

390     continue
        v(i,i) = 1.0d0
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

        U(i,l:n) = 0.0d0

430     if (g .eq. 0.0d0) go to 475
        if (i .eq. mn) go to 460

        do j = l, n
           s = dot_product(U(l:m,i),U(l:m,j))
           !     .......... double division avoids possible underflow ..........
           f = (s / u(i,i)) / g

           U(i:m,j) = U(i:m,j) + f*U(i:m,i)
        enddo

460     continue
        U(i:m,i) = U(i:m,i) / g 

        go to 490

475     continue
        U(i:m,i) = 0.0d0

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

        V(1:n,k) = -V(1:n,k)

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
    !
    !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
    !
    double precision :: a, b

    double precision :: p, r, s, t, u

    p = max(dabs(a),dabs(b))
    if (p .eq. 0.0d0) go to 20

    r = (min(dabs(a),dabs(b))/p)**2

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
     implicit   none
     ! 2015-Apr-28  15:45:16  ghislain: analysis
     ! subroutine to sort the indexes of elements stored in input vector INV 
     ! in reverse order in output vector OUTV.
     ! WS is a workspace vector, N is the dimension of the vectors.
     ! Supposition that INV contains positive numbers only!
     integer  :: n
     double precision  ::  inv(n), ws(n)
     integer  :: outv(n)
     
     integer  :: i, j, jmax
     
     WS(:n) = INV(:n)
     
     do j = 1,n
        jmax = 1
        do i = 1,n
           if (ws(i).gt.ws(jmax)) jmax = i
        enddo
        outv(n-j+1) = jmax
        ws(jmax) = 0.0
     enddo
     
     return
end subroutine rvord

subroutine setup(resp,a,ir,ic,nr,nc)
    implicit none
    ! ****************************************************
    !                                                    *
    !    Set DOUBLE PRECISION (A)                        *
    !    response matrix in FORTRAN storage format       *
    !    RESP comes from "C" routine                     *
    !                                                    *
    !     Author: WFH  05.02.02                          *
    !                                                    *
    ! ****************************************************
    integer :: ir, ic, nr, nc 
    double precision :: resp, a(nr, nc)

    a(ir+1,ic+1) = resp

    return
end subroutine setup

subroutine setupi(resp,a,ir,ic,nr,nc)
    implicit none
    ! ****************************************************
    !                                                    *
    !    Set INTEGER (A)                                 *
    !    matrix in FORTRAN storage format                *
    !    RESP comes from "C" routine                     *
    !                                                    *
    !     Author: WFH  05.02.02                          *
    !                                                    *
    ! ****************************************************
    integer :: ir, ic, nr, nc
    integer :: resp, a(nr, nc)

    a(ir+1,ic+1) = resp
 
    return
end subroutine setupi

subroutine primat(a,nc,nr) 
     implicit   none
     integer :: nr, nc
     integer :: a(nc, nr)     
     integer :: i, j
     
     do i = 1, nc
        write(*,*) (a(i,j), j=1, nr)
     enddo
     
     return
end subroutine primat
   
subroutine prdmat(a,nc,nr)
     implicit   none
     integer :: nr, nc
     double precision :: a(nc, nr)
     
     integer :: i, j
     
     do i = 1, nc
        write(*,*) (a(i,j), j=1, nr)
     enddo
     
     return
end subroutine prdmat

SUBROUTINE swapint(a, b)
  INTEGER, INTENT(IN OUT) :: a, b
  INTEGER :: temp
  temp = a ; a = b ; b = temp
END SUBROUTINE swapint

SUBROUTINE swapreal(a, b)
  DOUBLE PRECISION, INTENT(IN OUT) :: a, b
  DOUBLE PRECISION :: temp
  temp = a ; a = b ; b = temp
END SUBROUTINE swapreal

