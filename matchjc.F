      subroutine mtjac(ncon,nvar,strategy,cool,balance,random,          &
     &nrep,bisec,match_mode,                                            &
     &tol,calls,call_lim,                                               &
     &vect,dvect,fun_vec,                                               &
     &w_ifjac,w_iwa4,fval,xstart,xold)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   JACOBIAN command.                                                  *
! Attributes:                                                          *
!   ncon      (int)     # constraints                                  *
!   nvar      (int)     # variables                                    *
!   strategy  (int)     # strategy   1 normal                          *
!                                    2 print jacobian                  *
!                                    3 cancel variables                *
!   balance   (real)    # balance cooling factor, <0 use opt values    *
!   cool      (real)    # cooling factor                               *
!   random    (real)    # random  factor                               *
!   bisec     (int)     # bisec iteration number                       *
!   match_mode(int)     # mode use_macro=2                             *
!   nrep      (int)     # number of repetition                         *
!   tol       (real)    Final tolerance for match.                     *
!   calls     (int)     current call count                             *
!   call_lim  (int)     current call limit                             *
!   vect      (real)    variable values                                *
!   dvect     (real)    variable steps                                 *
!   fun_vect  (real)    function values                                *
!   all other working spaces for jacobian                              *
!----------------------------------------------------------------------*
      integer calls,call_lim,ncon,nvar
      integer strategy,nrep,i
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision tol,vect(*),dvect(*),fun_vec(*)
      double precision w_ifjac(*),w_iwa4(*)
      double precision fval(*),xold(*),xstart(*)
      double precision random,cool,balance
      integer bisec,match_mode
      include 'match.fi'
      external mtfcn

! mtfcn store the parameter using mtputi(x) and compute the function
! using mtcond

      icovar = 0
      ilevel = 0
! check if the variables are in the contraint and reset them if necessary
      call mtgeti(vect,dvect)
! call the main routine
      write(*,*) "JACOBIAN Strategy =", strategy

!  Repeat the program ntimes
!      write(*,*) nrep
      do i=1,nrep
        if (strategy .ge. 1) then
          calls=0
          call jacob(mtfcn,ncon,nvar,strategy,calls,call_lim,           &
     &vect,fun_vec,tol,                                                 &
     &w_ifjac,w_iwa4,fval,                                              &
     &xstart,xold,cool,balance,random,bisec,match_mode)
        endif
      enddo
      end


      subroutine jacob(fcn,m,n,strategy,calls,call_lim,                 &
     &x,fvec,epsfcn,                                                    &
     &fjac,wa4,fval,                                                    &
     &xstart,xold,cool,balance,random,bisec,match_mode)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   main JACOBIAN routine.                                             *
! Attributes:                                                          *
!   fcn       (real)    # function                                     *
!   m         (int)     # constraints                                  *
!   n         (int)     # variables                                    *
!   strategy  (int)     # strategy   1 normal,2 print jacobian,3 smart *
!   tol       (real)    Final tolerance for match.                     *
!   calls     (int)     current call count                             *
!   call_lim  (int)     current call limit                             *
!   x         (real)    variable values                                *
!   dvect     (real)    variable steps                                 *
!   fvec     (real)    function values                                 *
!                      mtfjac2 fjac,wa4                                *
!   cool      (real)    # cooling factor                               *
!   balance   (real)    # balance cooling factor                       *
!   random    (real)    # random  factor                               *
!   bisec     (int)     # bisec iteration number                       *
!   match_mode(int)     # mode use_macro=2                             *
!----------------------------------------------------------------------*

      integer izero
      integer i,iflag,info,j,level,m,n,calls,call_lim
      integer ireset,strategy,bisec,match_mode
      double precision fval(m),fmin_old,epsfcn
      double precision fnorm,ftol,gtol
      double precision dxnorm,xnorm,dx(n),fmin_start,fmin_old2
      double precision vmod
      double precision xtol,x(n),xstart(n),xold(n),fvec(m)
      double precision fjac(m,n),wa4(m),zero,one,two
      double precision epsil,epsmch,cool,balance,random
      parameter(zero=0d0,one=1d0,two=2d0, epsil=1d-9,epsmch=1d-16)
      include 'match.fi'
      external fcn, mtcond

!      double precision xdiff(n)
      double precision effjac(M,N),effsol(M+N),effrhs(M+N)
      double precision work(2*(M+N)),DNRM2
      integer effcon, effvar,coninfo(M),varinfo(N)

      ireset = 0
      izero = 0
      info = 0
      level = 0
      ftol = epsfcn
      gtol = epsil
      xtol = epsil

!---- Store the starting values in xstart for mtslope
      do j = 1, n
        xstart(j) = x(j)
      enddo

!---- Apply a cooling factor: reduce distance from solution to a point
!---- defined by balance and the limit or opt values
      call mtcool(x,cool,balance)

!---- Apply a random factor:
      call mtrandom(x,random)


!---- Check if the limit is within the constraint and
!---- reset the illegal values
      call mtlimit(x,ireset)

!---- Compute matching functions in fvec (penalty values)
!      call mtcond(izero, m, fvec, iflag)
      call FCN(M,N,X,fvec,IFLAG)
      calls=calls+1
      if (iflag .ne. 0) then
        call aawarn('JACOBIAN', ' stopped, possibly unstable')
        info = - 1
        go to 300
      endif

!---- Compute the norm of the function values (penalty values)
      fnorm = vmod(m, fvec)
      fmin = fnorm**2
      write(*,831) fmin
      fmin_start=fmin
!      fmin_old = fmin
!      edm = fmin
!     the first calls usually have a worse penalty due perhaps
!      to dimensional problems
!      fmin = 1E24


!---- Check the input parameters for errors.
      if (n .le. 0                                                      &
     &.or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero       &
     &.or. call_lim .le. 0                      ) then
        call aawarn('JACOBIAN', ' error in the input parameters')
        go to 300
      endif

!---- Quit, when initial value is already OK
!---- Do not aprly for calculating jacobian
      if (strategy.ne.2) then
        if (fmin .le. ftol) then
          call aawarn('JACOBIAN', ' penalty function already ok')
          go to 300
        endif
      endif
      if (ilevel .ge. 1) call mtprnt('old', n, x)

!---- Quit when call==1 (call set to 0)
!      if (calls .eq. 1) then
!        info=5
!        go to 300
!      endif

!      write(*,*) "fvec=",fvec
!      write(*,*) fmin,ftol
!      write(*,*) "x=",x

!---- Start loop
  20  continue
!---- Reset ireset
      ireset=0

!---- Calculate the jacobian
      call fdjac2(fcn,m,n,x,fvec,fjac,m,iflag,xtol,wa4)
      if (strategy.eq.2) then
!---- Print the jacobian and exit
        call jacob_print(m,n,fjac,match_mode)
        goto 300
      endif

!---- Reset solution vector
      do i=1,N+M
        effsol(i)=0
      enddo
!---- Reset varinfo
      do i=1,N
        varinfo(i)=0
      enddo

!---- All the variables are affective
      effvar=N

!---- Cancel the zero lines in the jacobian corresponding to
!---- an inequality that is not effective or variable not effective
!---- Reset constraint counter
      effcon=0
      do i=1,M
!---- Assume bad constraint: coninfo=1
        coninfo(i)=1
!---- Compute the norm of a row
!---- DNRM2 (N, DX, INCX)
!---- 2-dim array are stored by colums
        if (DNRM2(N,fjac(i,1),M).ge.1D-16) then
!      if (DNRM2(N,fjac(i,1),M).ge.0) then
!---- Good constraint: coninfo=0
          coninfo(i)=0
!---- Increase constraint counter
          effcon=effcon+1
          do j=1,N
!---- Copy RHS in the solution vector
            effsol(effcon)=fvec(i)
!---- Save a copy of RHS
            effrhs(effcon)=fvec(i)
!---- Update the effective jacobian
            effjac(effcon,j)=fjac(i,j)
          enddo
!        write(*,*) i," effective constraint"
        endif
      end do

 33   continue
!---- Check dimension of the system
!      if(effvar .lt. effcon) then
!        write(*,*) "constraints:",effcon,"  variables:",effvar
!      else
!        write(*,*) "Warning: underdetermined system"
!        write(*,*) "constraints:",effcon,"  variables:",effvar
!      endif

!---- Solve the least square problem and put solution in effsol vector
!---- Debug
!      write(*,*) "!!! pre solve routine"
!      write(*,*) "!!!effcon,effvar"
!      write(*,*) "!!!",effcon,effvar
!      write(*,*) "!!!i,j,effjac(i,j),effsol(j)"
!      do i=1,m
!      do j=1,n
!      write(*,*) "!!!",i,j,effjac(i,j),effsol(j)
!      enddo
!      enddo
!---- Debug
!---- CALL DGELS(TRANSA, M, N, NRHS, DA, LDA, DB, LDB, DWORK,
!---- LDWORK,INFO)
      write(*,*) "Solve system with ",effcon,"con,",effvar,"var"
      call DGELS ('N',effcon,effvar,1,effjac,M,effsol,N+M,              &
     &WORK,2*(N+M),INFO)
!---- Debug
!      write(*,*) "!!! solve routine"
!      write(*,*) "!!!",info,effsol(1)
!      do j=1,n
!      write(*,*) "!!!",j,effsol(j)
!      enddo
!---- Debug
      if (info.lt.0) then
        call aawarn('JACOBIAN', ' system solving routine failure')
        print *, '++++++++++ JACOBIAN ended: DGELS failure'
        print *, '++++++++++ JACOBIAN ended: info = ', info
        goto 300
      endif
      if (effsol(1).ne.effsol(1)) then
        print *, '++++++++++ JACOBIAN ended: NaN in system solving'
        goto 300
      endif


!---- Update the starting point
      effvar=0
!      write(*,*) "!!!,i,x(i),effvar,effsol(effvar)"
      do i=1,N
!----   Save in xold
        xold(i)=x(i)
        if (varinfo(i).eq.0) then
          effvar=effvar+1
          x(i)=x(i)-effsol(effvar)
!          write(*,*) "!!!",i,x(i),effvar,effsol(effvar)
        endif
      enddo

      call FCN(M,N,x,fvec,IFLAG)
      fnorm = vmod(m, fvec)
      fmin_old=fmin
      fmin = fnorm**2

      xnorm=DNRM2(N, x, 1)
      dxnorm=DNRM2(effvar, effsol, 1)/xnorm
      write(*,*) 'Step length', dxnorm
!      write(*,*) "tar=",fmin,fmin_old,ftol
!      write(*,830) calls,fmin,effcon,N
!      do i=1,N
!        xdiff(i)=xstart(i)-x(i)
!      enddo
!      write(*,830) calls,fmin,DNRM2(N, xdiff, 1)/DNRM2(N, xstart, 1)

      ! first fmin_old equal to fmin_start
      ! Bisection search
      fmin_old2=1E20
      j=0
36    continue
      if(fmin.ge.fmin_old .and. j.lt. bisec ) then
        ! go back to average solution
        do i=1,N
          x(i)=(x(i)+xold(i))*.5
          dx(i)=(x(i)-xold(i))*.5
        enddo
        j=j+1
        call FCN(M,N,x,fvec,IFLAG)
!        calls=calls+1
        fnorm = vmod(m, fvec)
        fmin = fnorm**2
        xnorm=DNRM2(N, x, 1)
        dxnorm=DNRM2(N, dx, 1)/xnorm
!        write(*,*) dxnorm,xnorm,fmin,fmin_old2
        if (fmin.gt.fmin_old2) then
!          do i=1,N
!            x(i)=2*x(i)-xold(i)
!          enddo
!          call FCN(M,N,x,fvec,IFLAG)
!!        calls=calls+1
!          fnorm = vmod(m, fvec)
!          fmin = fnorm**2
!          xnorm=DNRM2(N, x, 1)
!          write(*,*) dxnorm,xnorm,fmin,fmin_old2
          goto 37
        endif
        fmin_old2=fmin
        goto 36
      endif
37    continue

      if (j.gt.0) then
        write(*,*) 'Bisec iteration', j
      endif



!     Check for slope and limits and set the results in varinfo and
!     the number of effective variables in effvar
      if(strategy.eq.3) then
        call mtvarinfo(x,xstart,varinfo,effvar)
      endif

!      write(*,*) "!!! post info routne"
!      write(*,*) "!!!effcon,effvar,N"
!      write(*,*) "!!!",effcon,effvar,M,N

!     If needed recalculate the solution excluding some variables
      if ((effvar.lt.N).and.(effvar.ge.effcon)) then
        write(*,*) "Reset system to ",effcon,"con,",effvar,"var"
!      write(*,*) "effcon,effvar,i,j,effrhs(effcon),fjac(i,j)"

!----   Recalculate  effvar
        effvar=0
        do j=1,N
!        write(*,*) j,varinfo(j)
          if (varinfo(j).eq.0) then
            effvar=effvar+1
          endif
        enddo

!----   Restore the starting point
        do i=1,N
          x(i)=xold(i)
        enddo

!----   Reset RHS and jacobian
        do i=1,M
          do j=1,N
            effjac(i,j)=0
          enddo
        enddo
        do i=1,N+M
          effsol(i)=0
        enddo

!----   Rewrite the effective jacobian
        effcon=0
        do i=1,M
          if (coninfo(i).eq.0) then
            effcon=effcon+1
            effvar=0
!            write(*,*) i,effcon
            do j=1,N
              !reset the arrays
!            write(*,*) i,effcon,j
              if (varinfo(j).eq.0) then
                effvar=effvar+1
!                write(*,*) i,effcon,j,effvar
                effsol(effcon)=effrhs(effcon)
                effjac(effcon,effvar)=fjac(i,j)
!                write(*,*) "!!!",effcon,effvar,i,j,effrhs(effcon),fjac(i,j)
              endif
            enddo
          endif
        enddo

        if (ireset.eq.20) then
          write(*,*) "Too loops in system resizing, set strategy=1"
          strategy=1
!----     Exit var cancel loop
          goto 34
        else
          ireset=ireset+1
        endif

        if (effvar.lt.effcon) then
!----     Impossible to get  better
          write(*,*) "Too var to exclude, set strategy=1"
          strategy=1
!----     Exit var cancel loop
          goto 34
        endif

!----   Solve the system again
        goto 33
      endif

  34  continue

!      if(strategy.ne.3) then
!---- Check if the solution respect the slope and
!---- reset the illegal values
      call mtslope(x,xstart)

!---- Check if the limit is within the constraint and
!---- reset the illegal values
      call mtlimit(x,ireset)
!      write(*,*) "effsol=",effsol

!      endif

!---- calculate the target function and set new values
      call FCN(M,N,x,fvec,IFLAG)
      calls=calls+1
!      write(*,*) "fvec=",fvec

      fnorm = vmod(m, fvec)
      fmin_old=fmin
      fmin = fnorm**2

      xnorm=DNRM2(N, x, 1)
      dxnorm=DNRM2(effvar, effsol, 1)/xnorm
!      write(*,*) "tar=",fmin,fmin_old,ftol
!      write(*,830) calls,fmin,effcon,N
      write(*,830) calls,dxnorm,fmin
!      do i=1,N
!        xdiff(i)=xstart(i)-x(i)
!      enddo
!      write(*,830) calls,fmin,DNRM2(N, xdiff, 1)/DNRM2(N, xstart, 1)

!      ! Bisection search
!35    continue
!      if(fmin.ge.fmin_old .and. dxnorm.gt. bisec ) then
!        ! go back to average solution
!        do i=1,N
!          x(i)=(x(i)+xold(i))*.5
!          dx(i)=(x(i)-xold(i))*.5
!        enddo
!        call FCN(M,N,x,fvec,IFLAG)
!!        calls=calls+1
!        fnorm = vmod(m, fvec)
!        fmin = fnorm**2
!        xnorm=DNRM2(N, x, 1)
!        dxnorm=DNRM2(N, dx, 1)/xnorm
!        write(*,829) dxnorm,fmin
!        fmin_old=fmin
!!        write(*,*) "iterating",oldxnorm,xnorm,fmin_old,fmin
!        goto 35
!      endif
      

      ! check if the target explode
      if((fmin .gt. 1D20)) then
!       write(*,*) fmin, 1/xtol, (fmin .gt. 1D30)
        ! go back to the preovius solution
        do i=1,N
          x(i)=xold(i)
        enddo
        ! set values to the previous solution
        call FCN(M,N,X,fvec,IFLAG)
        print *, '++++++++++ JACOBIAN ended: infinity penalty function'
        goto 300
      endif

      ! check if the function converge
      if(fmin .le. ftol) then
        print *, '++++++++++ JACOBIAN ended: converged successfully'
        goto 300
      endif
      ! Iteration do not converge
      ! For the first calls it can happens due to not fisical penalty
      if(fmin.ge.fmin_old .and. fmin.lt.1D-18) then
        ! go back to the preovius solution
        do i=1,N
          x(i)=x(i)+effsol(i)
        enddo
        ! set values to the previous solution
        call FCN(M,N,X,fvec,IFLAG)
        print *, '++++++++++ JACOBIAN ended: low penalty oscillations'
        goto 300
      endif
      if(calls.ge.call_lim) then
        print *, '++++++++++ JACOBIAN ended: call limit'
        goto 300
      endif

!     restart loop
      goto 20

!      open(25,file="jacobian.dat")
!      write(25,*) "# ", m,n
!      do i=1,m
!      do j=1,n
!      write(25,999) i,j,fjac(i,j)
!      enddo
!      enddo
!      close(25)

!      call jacob_print(m,n,fjac,match_mode)
!  999 format(i3,i3,e16.8)
  300 continue

!---- Store the final distance
      do i = 1, n
        dx(i)=(x(i)-xstart(i))
      enddo
      
      xnorm=DNRM2(N, x, 1)
      dxnorm=sqrt(DNRM2(N, dx, 1)/xnorm)
      write(*,*) 'final relative difference norm:',dxnorm

  831 format('Initial Penalty Function = ',e16.8,//)
  830 format('call: ',I5, ' Dx = ', e16.8,                              &
     &'  Penalty function = ',e16.8)
      end

      subroutine jacob_print2(m,n,fjac,match_mode)

      implicit none

      integer m,n,ivar,nvar,match_mode
      double precision fjac(m,n)
      integer i,j

      do i=1,m
      do j=1,n
        write(*,*) i,j,fjac(i,j)
      enddo
      enddo

      end

      subroutine jacob_print(m,n,fjac,match_mode)
!      subroutine jacob_print(ncon,fsum,fvect)

      implicit none

      integer m,n,ivar,nvar,match_mode
      double precision fjac(m,n)
      logical local
      integer ncon,next_constraint,next_global,i,j,pos,type,range(2),   &
     &flag,get_option,restart_sequ,advance_to_pos,name_l,char_from_table
!      integer double_from_table
      parameter(name_l=24)
      double precision value,c_min,c_max,weight
      include 'name_len.fi'
      character*(name_len) namevar,name,node_name
      integer next_vary,slope
      double precision step,opt


      if(match_mode.eq.1) then

      write(*,997) "Node Name    ","Constraint      ",                  &
     &"Variable        ","Derivative    "
      write(*,996) "----------------------------------------------------&
     &----------------"
!      write(25,*) "# ", m,n
      ncon=1
      local=get_option('match_local ') .ne. 0
      call table_range('twiss ','#s/#e ',range)
      if(local) then
        j=restart_sequ()
        do pos=range(1),range(2)
          j=advance_to_pos('twiss ',pos)
 20       continue
          i=next_constraint(name,name_l,type,value,c_min,c_max,weight)
          if(i.ne.0)  then
            flag=char_from_table('twiss ','name ',pos,node_name)
            do nvar=1,n
 22           ivar=next_vary(namevar,name_l,c_min,c_max,step,slope,opt)
              if (ivar.eq.0) then
                goto 22
              endif
              write (*,998) node_name,name,namevar,fjac(ncon,ivar)
!            write(25,999) ncon,ivar,fjac(ncon,ivar)
            enddo
            ncon=ncon+1
            goto 20
          endif
        enddo
      endif
 30   continue
      i=next_global(name,name_l,type,value,c_min,c_max,weight)
!      write(*,*) i,name
      if(i.ne.0)  then
        pos=1
!          flag=double_from_table('summ ',name,pos,val)
        do nvar=1,n
 32       ivar=next_vary(namevar,name_l,c_min,c_max,step,slope,opt)
          if (ivar.eq.0) then
            goto 32
          endif
!            write(*,*) ncon,ivar,fjac(ncon,ivar)
          write (*,998) 'Global',name,namevar,fjac(ncon,ivar)
        enddo
        ncon=ncon+1
        goto 30
      endif

!      do nvar=1,n
! 23   ivar=next_vary(namevar,name_l,c_min,c_max,step,slope,opt)
!      if (ivar.eq.0) then
!        goto 23
!      endif
!      write(25,*) namevar
!      enddo
!      close(25)
      endif

  996 format(a)
  997 format(3(a16,1x),a16)
  998 format(3(a16,1x),e16.8)
  999 format(i3,i3,e16.8)
      end



      subroutine mtcool(vect,cool,balance)

      implicit none


      integer j,next_vary,name_l,slope
      parameter(name_l=24)
      double precision vect(*),c_min,c_max,step,opt
      double precision val,cool,balance
      include 'name_len.fi'
      character*(name_len) name

 1    continue
      j = next_vary(name,name_l,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        if (balance.ge.0) then
          val = (1-cool)*vect(j)+cool*((1-balance)*c_max+balance*c_min)
        else
          val =(1-cool)*vect(j)+cool*opt
        endif
        vect(j) = val
        goto 1
      endif
      end

      subroutine mtrandom(vect,random)

      implicit none


      integer j,next_vary,name_l,slope
      parameter(name_l=24)
      double precision vect(*),c_min,c_max,step,opt
      double precision val,random,frndm
      include 'name_len.fi'
      character*(name_len) name

 1    continue
      j = next_vary(name,name_l,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        val = (1+ random *( frndm() - 0.5) ) * vect(j)
        vect(j) = val
        goto 1
      endif
      end

      subroutine mtslope(x,xstart)

      implicit none


      integer j,next_vary,name_l,slope
      parameter(name_l=24)
      double precision x(*),xstart(*),c_min,c_max,step,opt
      double precision diff
      include 'name_len.fi'
      character*(name_len) name

 1    continue
      j = next_vary(name,name_l,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        if (slope.ne.0) then
          diff=x(j)-xstart(j)
          if(slope*diff.lt.0) then
            write(*,831) "reset parameter:",name,                       &
     &"from",x(j),"to",xstart(j)
            x(j)=xstart(j)
          endif
        endif
        goto 1
      endif
  831 format(a16,1x,a24,a4,e16.8,a4,e16.8)
      end


      subroutine mtvarinfo(x,xstart,varinfo,effvar)

      implicit none


      integer j,next_vary,name_l,slope,varinfo(*),effvar
      parameter(name_l=24)
      double precision x(*),xstart(*),c_min,c_max,step,opt
      double precision diff,val,oldval
      include 'name_len.fi'
      character*(name_len) name

      effvar=0
 1    continue
      j = next_vary(name,name_l,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
!        varinfo(j)=0
        effvar=effvar+1
        val=x(j)
        oldval=xstart(j)
        if (slope.ne.0) then
          diff=val-oldval
          if(slope*diff.lt.0) then
            write(*,*) "exlude parameter:",name,"bad slope"
            varinfo(j)=1
            effvar=effvar-1
          endif
        endif
        if (val.lt.c_min) then
          write(*,*) "exlude parameter:",name,"hit minimum"
          varinfo(j)=1
          effvar=effvar-1
        endif
        if (val.gt.c_max) then
          write(*,*) "exlude parameter:",name,"hit maximum"
          varinfo(j)=1
          effvar=effvar-1
        endif
        goto 1
      endif
  831 format(a16,1x,a24,a4,e16.8,a4,e16.8)
      end

!  TODO
!  Find target function in collect
!  Find variable constraint in vary

!  Loop for divide the funtion in steps


