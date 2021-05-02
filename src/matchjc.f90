      subroutine mtjac(ncon,nvar,strategy,cool,balance,random,          &
     &nrep,bisec,cond,match_mode,                                       &
     &tol,calls,call_lim,                                               &
     &vect,dvect,fun_vec,                                               &
     &w_ifjac,w_iwa4,fval,xstart,xold)

      use matchfi
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
!   cond      (real)    # condition number for rank                    *
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
      double precision random,cool,balance,cond
      integer bisec,match_mode
      external mtfcn

      fval(1)=0

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
!          calls=0
          call jacob(mtfcn,ncon,nvar,strategy,calls,call_lim,           &
     &vect,fun_vec,tol,                                                 &
     &w_ifjac,w_iwa4,                                                   &
     &xstart,xold,cool,balance,random,bisec,cond,match_mode)
        endif
      enddo
      end


      subroutine jacob(fcn,m,n,strategy,calls,call_lim,                 &
     &x,fvec,epsfcn,                                                    &
     &fjac,wa4,                                                         &
     &xstart,xold,cool,balance,random,bisec,cond,match_mode)

      use matchfi
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
!   cond      (real)    # condition number for rank                    *
!   match_mode(int)     # mode use_macro=2                             *
!----------------------------------------------------------------------*

      integer izero
      integer i,iflag,info,j,level,m,n,calls,call_lim,ecalls
      integer ireset,strategy,bisec,match_mode
      double precision fmin_old,epsfcn
      double precision ftol,gtol
      double precision dxnorm,xnorm,dx(n),fmin_start,fmin_old2
      double precision vdot
      double precision xtol,x(n),xstart(n),xold(n),fvec(m)
      double precision xopt(n),xbest(n),fminbest,condnum
      double precision fjac(m,n),wa4(m),zero,one,two
      double precision epsil,epsmch,cool,balance,random
      parameter(zero=0d0,one=1d0,two=2d0, epsil=1d-9,epsmch=1d-16)
      external fcn, mtcond

!      double precision xdiff(n)
      double precision effjac(M,N),effsol(M+N),effrhs(M+N)
      double precision DNRM2
      double precision WORK(1000*(N+M)),SV(N+M),COND
      integer IWORK(30*(N+M)),RANK
      integer effcon, effvar,coninfo(M),varinfo(N)
      integer debug
      integer, external :: get_option

      debug = get_option('debug ')
      ecalls = 0
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
      call mtcool(x,cool,balance,xopt)

!---- Apply a random factor:
      call mtrandom(x,random)

!---- Check if the limit is within the constraint and
!---- reset the illegal values
      call mtlimit(x,ireset)

!---- Compute matching functions in fvec (penalty values)
      call FCN(M,N,X,fvec,IFLAG)
      calls=calls+1
      if (iflag .ne. 0) then
        call fort_warn('JACOBIAN', ' stopped, possibly unstable')
        info = - 1
        go to 300
      endif

!---- Compute the norm of the function values (penalty values)
      fmin = vdot(m,fvec,fvec)
      write(*,'("Initial Penalty Function = ",e24.16,//)') fmin
      fmin_start=fmin

!---- store the first best set
      fminbest=fmin
      do j = 1, n
        xbest(j) = x(j)
      enddo

!---- Check the input parameters for errors.
      if (n .lt. 0                                                      &
     &.or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero       &
     &.or. call_lim .le. 0                      ) then
        call fort_warn('JACOBIAN', ' error in the input parameters')
        go to 300
      endif

!---- Quit, when initial value is already OK
!---- Do not apply for calculating jacobian
      if (strategy.ne.2) then
        if (fmin .le. ftol) then
          call fort_warn('JACOBIAN', ' penalty function already ok')
          go to 300
        endif
      endif
      if (ilevel .ge. 1) call mtprnt('old', n, x)


!---- Start main loop
  20  continue

      if (debug .ne. 0) then
        write(*,*) "Current variables values"
        do i=1,N
          write(*,'(i5,e20.12)') i, x(i)
        enddo
      endif

!---- Reset ireset
      ireset=0

!---- Calculate the jacobian
      call fdjac2(fcn,m,n,x,fvec,fjac,m,iflag,xtol,wa4)
      ecalls=ecalls+n
      if (iflag .ne. 0) then
        call fort_warn('JACOBIAN', ' stopped, possibly unstable')
        info = - 1
        go to 300
      endif

      if (strategy.eq.2) then
!        call DGESDD('A',M,N,fjac,M,SV,U,M,VT,N,                         &
!     & WORK, 1000*(N+M), IWORK, INFO )
!---- Print the jacobian on the match2 variables and exit
        call jacob_print(n,match_mode)
        goto 300
      endif

      if (debug .ne. 0) then
        write(*,*) "Jacobian: "
        do i=1,M
          write(*,'(100e20.12)') (fjac(i,j),j=1,N)
        enddo
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
!---- 2-dim array are stored by colums
        if (DNRM2(N,fjac(i,1),M).ge.1D-16) then
!---- Good constraint: coninfo=0
          coninfo(i)=0
!---- Increase constraint counter
          effcon=effcon+1
!---- Copy RHS in the solution vector
          effsol(effcon)=fvec(i)
          do j=1,N
!---- Update the effective jacobian
            effjac(effcon,j)=fjac(i,j)
          enddo
!---- Save a copy of RHS
          effrhs(effcon)=effsol(effcon)
        endif
      end do

!-- BUG, effjac and effsol can be empty in case jacobian is e.g. full of zero...

 33   continue


!---- Solve the least square problem and put solution in effsol vector
!---- Debug
!      write(*,*) "!!! pre solve routine"
!      write(*,*) "!!!effcon,effvar"
!      write(*,*) "!!!",effcon,effvar
!      write(*,*) "!!!i,j,effjac(i,j),effsol(j)"
!      do i=1,effcon
!      do j=1,effvar
!      write(*,*) "!!!",i,j,effjac(i,j)
!      enddo
!      enddo
!---- Debug
!---- CALL DGELS(TRANSA, M, N, NRHS, DA, LDA, DB, LDB, DWORK,
!---- LDWORK,INFO)
      write(*,*) "Solve system with ",effcon,"con,",effvar,"var"

      if (debug .ne. 0) then
        write(*,*) "Effective Jacobian: "
        do i=1,effcon
          write(*,'(100e20.12)') (effjac(i,j),j=1,effvar)
        enddo
      endif

!      call DGELS ('N',effcon,effvar,1,effjac,M,effsol,N+M,              &
!     &WORK,2*(N+M),INFO)
!      DGELSD( M, N, NRHS, A, LDA, B, LDB, S,  RCOND,
!      RANK,  WORK, LWORK, IWORK, INFO )
      call DGELSD (effcon,effvar,1,effjac,M,effsol,N+M,SV,COND,         &
     &RANK,WORK,1000*(N+M),IWORK,INFO)
      condnum=SV(1)/SV(min(effcon,effvar))
!      call DGELSS (effcon,effvar,1,effjac,M,effsol,N+M,SV,RCOND,        &
!     &RANK,WORK,1000*(N+M),INFO)
      write(*,*) "Rank  ",RANK,                                         &
     &"  Condition number ",condnum

!---- Debug
!      write(*,*) "!!! solve routine"
!      write(*,*) "!!!",info,effsol(1)
!      do j=1,n
!      write(*,*) "!!!",j,effsol(j)
!      enddo
!---- Debug
      if (info.lt.0) then
        call fort_warn('JACOBIAN', ' system solving routine failure')
        print *, '++++++++++ JACOBIAN ended: DGELSD failure'
        print *, '++++++++++ JACOBIAN ended: info = ', info
        goto 300
      endif
      if (effsol(1).ne.effsol(1)) then
        print *, '++++++++++ JACOBIAN ended: NaN in system solving'
        goto 300
      endif


!---- Update the starting point
      effvar=0
      do i=1,N
!----   Save in xold
        xold(i)=x(i)
        if (varinfo(i).eq.0) then
          effvar=effvar+1
          x(i)=x(i)-effsol(effvar)
        endif
      enddo


!     Check for slope and limits and set the results in varinfo and
!     the number of effective variables in effvar
      if(strategy.eq.3) then
        call mtvarinfo(x,xstart,varinfo,effvar)
      endif

!     If needed recalculate the solution excluding some variables
      if ((effvar.lt.N).and.(effvar.ge.effcon)) then
        write(*,*) "Reset system to ",effcon,"con,",effvar,"var"

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
            do j=1,N
              if (varinfo(j).eq.0) then
                effvar=effvar+1
                effsol(effcon)=effrhs(effcon)
                effjac(effcon,effvar)=fjac(i,j)
              endif
            enddo
          endif
        enddo

        if (ireset.eq.20) then
          write(*,*) "Too many loops in system resizing, set strategy=1"
          strategy=1
!----     Exit var cancel loop
          goto 34
        else
          ireset=ireset+1
        endif

        if ((effvar.lt.effcon+1).or.(condnum.gt.1E10)) then
!----     Impossible to get  better
          write(*,*) "Too many variables to exclude, set strategy=1"
          strategy=1
!----     Exit var cancel loop
          goto 34
        endif

!----   Solve the system again
        goto 33
      endif

  34  continue

!---- Check if the solution respect the slope and
!---- reset the illegal values
      call mtslope(x,xstart)

!---- Check if the limit is within the constraint and
!---- reset the illegal values
      call mtlimit(x,ireset)

!----- Calculate the penalty function
      call FCN(M,N,x,fvec,IFLAG)
      ecalls=ecalls+1
      if (iflag .ne. 0) then
        call fort_warn('JACOBIAN', ' stopped, possibly unstable')
        info = - 1
        go to 300
      endif

      fmin_old=fmin
      fmin = vdot(m, fvec,fvec)

      xnorm=DNRM2(N, x, 1)
      dxnorm=DNRM2(effvar, effsol, 1)/xnorm
      write(*,*) 'Step length ', dxnorm

      ! Bisection search
      fmin_old2=1E20
      j=0
36    continue
      if(fmin.ge.fmin_old .and. j.lt. bisec ) then
        ! go back to average solution
        do i=1,N
          x(i)=(x(i)+xold(i))*.5
          dx(i)=(x(i)-xold(i))*.5                   ! LD: never used
        enddo
        j=j+1
        call FCN(M,N,x,fvec,IFLAG)
        ecalls=ecalls+1
        if (iflag .ne. 0) then
          call fort_warn('JACOBIAN', ' stopped, possibly unstable')
          info = - 1
          go to 300
        endif
        fmin = vdot(m,fvec,fvec)
        dxnorm=sqrt(vdot(N, dx, dx)/vdot(N, x, x))  ! LD: never used
        if (fmin.gt.fmin_old2) then
          goto 37
        endif
        fmin_old2=fmin
        goto 36
      endif
37    continue

      if (j.gt.0) then
        write(*,*) 'Bisec iteration ', j
      endif


!---- calculate the target function and set new values
      call FCN(M,N,x,fvec,IFLAG)
      calls=calls+1
      if (iflag .ne. 0) then
        call fort_warn('JACOBIAN', ' stopped, possibly unstable')
        info = - 1
        go to 300
      endif

      fmin = vdot(m,fvec,fvec)
      xnorm=DNRM2(N, x, 1)
      dxnorm=DNRM2(effvar, effsol, 1)/xnorm
      write(*,"('call: ',I5,' Dx = ', e16.8,                            &
     &'  Penalty function =',e24.16)")  calls,dxnorm,fmin

!      Store the best point:
      if ((fmin .lt. fminbest)) then
        do i=1,N
          xbest(i)=x(i)
        enddo
        fminbest=fmin
      endif

      ! check if the target explode
      if((fmin .gt. 1D20)) then
        print *, '++++++++++ JACOBIAN ended: infinity penalty function'
        goto 300
      endif

      ! check if the function converge
      if(fmin .le. ftol) then
        print *, '++++++++++ JACOBIAN ended: converged successfully'
        goto 300
      endif
      ! Iteration do not converge
      ! For the first calls it can happen due to not physical penalty
      if((calls.ge.3).and.(abs(1-fmin/fmin_old).lt.1D-15) ) then
!        print *, 'dbg',fmin,fmin_old,abs(1-fmin/fmin_old)
        print *, '++++++++++ JACOBIAN ended: minimum found'
        goto 300
      endif
      if(calls.ge.call_lim) then
        print *, '++++++++++ JACOBIAN ended: call limit'
        goto 300
      endif

!     restart main loop
      goto 20

!      open(25,file="jacobian.dat")
!      write(25,*) "# ", m,n
!      do i=1,m
!      do j=1,n
!      write(25,999) i,j,fjac(i,j)
!      enddo
!      enddo
!      close(25)

!      call jacob_print(n,match_mode)
!  999 format(i3,i3,e16.8)
  300 continue

!     go back to the best solution
      do i=1,N
        x(i)=xbest(i)
      enddo
      ! set values to the previous solution
      call FCN(M,N,X,fvec,IFLAG)
      ecalls=ecalls+1
      if (iflag .ne. 0) then
        call fort_warn('JACOBIAN', ' stopped, possibly unstable')
        info = - 1
      endif

!---- Store the final distance
      do i = 1, n
        dx(i)=(x(i)-xstart(i))
      enddo

      if (debug .ne. 0) then
        write(*,*) 'Effective number of calls:', calls+ecalls
      endif

!      xnorm=DNRM2(N, x, 1)
      dxnorm=sqrt(DNRM2(N, dx, 1))
      write(*,*) 'Final difference norm:',dxnorm

      end

      subroutine jacob_print(n,match_mode)

      use name_lenfi
      implicit none

      integer n,ivar,nvar,match_mode
      logical local
      integer ncon,next_constraint,next_global,i,j,pos,type,range(2),   &
     &flag,get_option,restart_sequ,advance_to_pos,string_from_table_row
!      integer double_from_table_row
      double precision value,c_min,c_max,weight,val
      character*(name_len) namevar,name,node_name
      integer next_vary,slope
      double precision step,opt
      integer oldpos,nnode,mtputconsname,void


      if(match_mode.eq.1) then

        ncon=1
        nnode=0
        local=get_option('match_local ') .ne. 0
        call table_range('twiss ','#s/#e ',range)
        if(local) then
          j=restart_sequ()
          oldpos=range(1)
          do pos=range(1),range(2)
            j=advance_to_pos('twiss ',pos)
 20         continue
            i=next_constraint(name,name_len,type,value,c_min,c_max,weight, &
     &                        pos,val,node_name,name_len)
            if(i.ne.0)  then
              if (pos.ne.oldpos) then
                nnode=nnode+1
                ncon=1
                oldpos=pos
              endif
              do nvar=1,n
 22             ivar=next_vary(namevar,name_len,c_min,c_max,step,slope,opt)
                if (ivar.eq.0) then
                  goto 22
                endif
                void=mtputconsname(node_name,nnode,name,ncon)
              enddo
              ncon=ncon+1
              goto 20
            endif
          enddo
        endif
        nnode=nnode+1
 30     continue
        i=next_global(name,name_len,type,value,c_min,c_max,weight)
        if(i.ne.0)  then
          pos=1
          do nvar=1,n
 32         ivar=next_vary(namevar,name_len,c_min,c_max,step,slope,opt)
            if (ivar.eq.0) then
              goto 32
            endif
            void=mtputconsname('Global              ',nnode,name,ncon)
          enddo
          ncon=ncon+1
          goto 30
        endif

      endif

  996 format(a)
  997 format(3(a16,1x),a16)
  998 format(3(a16,1x),e16.8)
  999 format(i3,i3,e16.8)
      end

      subroutine mtmove(vect,varinfo,dir,balance)
      use name_lenfi
      implicit none
      integer j,next_vary,slope
      integer varinfo,effvar
      double precision vect(*),c_min,c_max,step,opt
      double precision val,dir,balance
      character*(name_len) name

      effvar=0
 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        if (varinfo(j).eq.0) then
          effvar=effvar+1
          if (opt.ge.0) then
            val =vect(effvar)+dir*opt
          else
            val =vect(effvar)+dir*((1-balance)*c_max+balance*c_min)
          endif
          vect(effvar) = val
          goto 1
        endif
      endif
      end


      subroutine mtcool(vect,cool,balance,xopt)
      use name_lenfi
      implicit none
      integer j,next_vary,slope
      double precision vect(*),c_min,c_max,step,opt
      double precision val,xopt(*),cool,balance
      character*(name_len) name

      if (cool.gt.0) then
        write(*,'(4a16)') "name","oldvalue","opt value","new value"
      endif
 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        if (opt.gt.0) then
          xopt(j)=opt
        else
          xopt(j)=(1-balance)*c_max+balance*c_min
        endif
        val=(1-cool)*vect(j)+cool*xopt(j)
        if (cool.gt.0) then
          write(*,'(a15,3e16.5)') name,vect(j),xopt(j),val
        endif
        vect(j) = val
        goto 1
      endif
      end

      subroutine mtrandom(vect,random)
      use name_lenfi
      implicit none
      integer j,next_vary,slope
      double precision vect(*),c_min,c_max,step,opt
      double precision val,random,frndm
      character*(name_len) name

 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        val = (1+ random *( frndm() - 0.5) ) * vect(j)
        vect(j) = val
        goto 1
      endif
      end

      subroutine mtslope(x,xstart)

      use name_lenfi
      implicit none


      integer j,next_vary,slope
      double precision x(*),xstart(*),c_min,c_max,step,opt
      double precision diff
      character*(name_len) name

 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
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

      use name_lenfi
      implicit none


      integer j,next_vary,slope,varinfo(*),effvar
      double precision x(*),xstart(*),c_min,c_max,step,opt
      double precision diff,val,oldval
      character*(name_len) name

      effvar=0
 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
!        varinfo(j)=0
        effvar=effvar+1
        val=x(j)
        oldval=xstart(j)
        if (slope.ne.0) then
          diff=val-oldval
          if(slope*diff.lt.0) then
            write(*,*) "exclude parameter:",name,"bad slope"
            varinfo(j)=1
            effvar=effvar-1
          endif
        endif
        if (val.lt.c_min) then
          write(*,*) "exclude parameter:",name,"hit minimum"
          varinfo(j)=1
          effvar=effvar-1
        endif
        if (val.gt.c_max) then
          write(*,*) "exclude parameter:",name,"hit maximum"
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



      subroutine mtsvd(M,N,fjac,SV,U,VT)
      implicit none
      integer M,N,IWORK(30*(N+M)),INFO
      double precision fjac(m,n),SV(N+M),U(M,M),VT(N,N),WORK(1000*(N+M))
      call DGESDD('A',M,N,fjac,M,SV,U,M,VT,N,WORK,1000*(N+M),IWORK,INFO)
      end
