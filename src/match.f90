      subroutine mtgetc(vect,dvect)
!----------------------------------------------------------------------*
! This is exactly equivalent to mtgeti but this version is used for    *
! calling from C.                                                      *
!                                                                      *
! This seems to be the only way deal with the problem that originally  *
! the subroutine mtgeti was called both from C and from Fortran.       *
! Because madX is setup with Linux-specific name-mangling, the         *
! adjustments made with the !DEC_dollar ATTRIBUTES do not seem able to *
! cope with such a case.                                               *
! Subroutine added at 10:38:56 on 8 Apr 2003 by JMJ                    *
!----------------------------------------------------------------------*
      implicit none
      double precision vect(*),dvect(*)
      call mtgeti(vect,dvect)
      end

      subroutine mtgeti(vect,dvect)
      use name_lenfi
      implicit none
! This subroutine should only be called from Fortran.  There is an
! equivalent mtgetc for the occasions when it needs to be called from C.
! Modified at 10:38:56 on 8 Apr 2003 by JMJ


      logical psum
      integer j,next_vary,get_option,slope
      double precision get_variable,vect(*),dvect(*),c_min,c_max,step,  &
     &dval,val,s_fact,valold,eps,eps2,stplim, vmax, vmin, opt
      parameter(s_fact=5d-1)
      parameter(eps = 1.0d-10,eps2 = 1.0d-1,stplim = 2.0d-1)
      parameter(vmax=1.e+20,vmin=-1.e+20)
      character*(name_len) name

      psum=get_option('match_summary ') .ne. 0
 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        val = get_variable(name)
        if (val .ge. c_max) then
          valold = val
          dval = min(step, (val - c_max)*s_fact)
!          write(*,*) step, c_min, val, s_fact
          val = c_max - 2*dval
          write(*,*) "reset parameter:",name,"from",valold,"to",val
        elseif (val .le. c_min) then
          valold = val
          dval = min(step, (c_min - val)*s_fact)
!          write(*,*) step, c_min, val, s_fact
          val = c_min + 2*dval
          write(*,831) "reset parameter:",name,"from",valold,"to",val
        else
          dval = step
        endif
        if(psum) write(*,830) name,val,c_min,c_max
        vect(j) = val
        dvect(j) = dval
        goto 1
      endif
!  830 format(a23,1x,1p,'=',e16.8,';!',2x,e16.8,2x,e16.8)
  830 format(a24,1x,1p,e16.8,3x,e16.8,3x,e16.8)
  831 format(a16,1x,a24,a4,e16.8,a4,e16.8)
      end



      subroutine mtlimit(vect,ireset)

      use name_lenfi
      implicit none


      integer j,next_vary,ireset,slope
      double precision vect(*),c_min,c_max,step,                        &
     &dval,val,s_fact,valold,eps,eps2,stplim, vmax, vmin, opt
      parameter(s_fact=5d-1)
      parameter(eps = 1.0d-10,eps2 = 1.0d-1,stplim = 2.0d-1)
      parameter(vmax=1.e+20,vmin=-1.e+20)
      character*(name_len) name

 1    continue
      j = next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if (j .ne. 0)  then
        val = vect(j)
        if (val .ge. c_max) then
          valold = val
          dval = min(step, (val - c_max)*s_fact)
          val = c_max - 2*dval
          write(*,831) "reset parameter:",name,"from",valold,"to",val
          ireset = ireset + 1
        elseif (val .le. c_min) then
          valold = val
          dval = min(step, (c_min - val)*s_fact)
          val = c_min + 2*dval
          write(*,831) "reset parameter:",name,"from",valold,"to",val
          ireset = ireset + 1
        endif
        vect(j) = val
        goto 1
      endif
  831 format(a16,1x,a24,a4,e16.8,a4,e16.8)
      end


      subroutine collect(ncon,fsum,fvect)

      use name_lenfi
      use twiss0fi
      use twisscfi
      implicit none


      logical fprt,local,psum, slow_match
      integer ncon,next_constraint,next_global,i,j,pos,type,                &
     &flag,get_option,restart_sequ,advance_to_pos,double_from_table_row,    &
     &string_from_table_row
      double precision fsum,fvect(*),val,valhg,c_min,c_max,weight,f_val
      character*(name_len) name, node_name
      integer n_pos, next_constr_namepos, advance_node
      local=get_option('match_local ') .ne. 0
      fprt=get_option('match_print ') .ne. 0
      psum=get_option('match_summary ') .ne. 0
      slow_match = get_option('slow_match ') .ne. 0
      if (local) then
        j=restart_sequ()
        pos=1
        do while (j .gt. 0)
          if (slow_match) j=advance_to_pos('twiss ',pos) ! (expensive) NOP?
          do while (next_constraint(                                    &
     &                  name,name_len,type, valhg,c_min,c_max,weight,   &
     &                  slow_match,pos,val,node_name,name_len).ne.0)
            if (type.eq.1) then
              f_val =weight*dim(c_min,val)
              if(fprt) write(*,880) name,weight,val,c_min,f_val**2
            elseif(type.eq.2) then
              f_val=weight*dim(val,c_max)
              if(fprt) write(*,890) name,weight,val,c_max,f_val**2
            elseif(type.eq.3) then
              f_val=weight*dim(c_min,val)+weight*dim(val,c_max)
              if(fprt) write(*,840) name,weight,val,c_min,c_max,f_val**2
            elseif(type.eq.4) then
              f_val=weight*(val-valhg)
              if(fprt) write(*,840) name,weight,val,valhg,valhg,f_val**2
            endif
            ncon=ncon+1
            fvect(ncon)=f_val
            fsum=fsum+f_val**2
            if(psum .and. type.eq.4)                                    &
     &write(*,830) node_name,name,type,valhg,val,f_val**2
            if(psum .and. type.eq.2)                                    &
     &write(*,830) node_name,name,type,c_max,val,f_val**2
            if(psum .and. type.eq.1)                                    &
     &write(*,830) node_name,name,type,c_min,val,f_val**2
            if(psum .and. type.eq.3)                                    &
     &write(*,832) node_name,name,type,c_min,c_max,val,f_val**2
          end do
          j=advance_node()
          pos=pos+1
        end do
      endif
 30   continue
      i=next_global(name,name_len,type,valhg,c_min,c_max,weight)
      if(i.ne.0)  then
        pos=1
        flag=double_from_table_row('summ ',name,pos,val)
        if(type.eq.1) then
          f_val=weight*dim(c_min,val)
          if(fprt) write(*,880) name,weight,val,c_min,f_val**2
        elseif(type.eq.2) then
          f_val=weight*dim(val,c_max)
          if(fprt) write(*,890) name,weight,val,c_max,f_val**2
        elseif(type.eq.3) then
          f_val=weight*dim(c_min,val)+ weight*dim(val,c_max)
          if(fprt) write(*,840) name,weight,val,c_min,c_max,f_val**2
        elseif(type.eq.4) then
          f_val=weight*(val-valhg)
          if(fprt) write(*,840) name,weight,val,valhg,valhg,f_val**2
        endif
        ncon=ncon+1
        fvect(ncon)=f_val
        fsum=fsum+f_val**2

        if(psum)                                                        &
     &write(*,830) "Global constraint:      ",name,type,valhg,val,      &
     &f_val**2

        goto 30
      endif
!      if (psum) write(*,831) fsum
  830 format(a24,3x,a6,5x,i3,3x,1p,e16.8,3x,e16.8,3x,e16.8)
  832 format(a24,3x,a6,5x,i3,3x,1p,e16.8,3x,e16.8,3x,e16.8,3x,e16.8)
  831 format(//,'Final Penalty Function = ',e16.8,//)
  840 format(a10,3x,1p,5e16.6)
  880 format(a10,3x,1p,3e16.6,16x,e16.6)
  890 format(a10,3x,1p,2e16.6,16x,2e16.6)
      end


      subroutine mtfcn(nf,nx,x,fval,iflag)
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Compute matching functions.                                        *
! Input:                                                               *
!   NF        (integer) Number of functions to be computed.            *
!   NX        (integer) Number of input parameters.                    *
!   X(NX)     (real)    Input parameters.                              *
! Output:                                                              *
!   FVAL(NF)  (real)    Matching functions computed.                   *
!   IFLAG     (integer) Stability flag.                                *
!----------------------------------------------------------------------*
      integer iflag,nf,nx,izero,ione
      double precision fval(nf),x(nx)

      izero = 0
      ione = 1
!---- Store parameter values in data structure.
      call mtputi(x)

!---- Compute matching functions.

      call mtcond(izero, nf, fval, iflag)

 9999 end
      subroutine mtputi(vect)

      use name_lenfi
      implicit none


      integer j,next_vary,slope
      double precision vect(*),c_min,c_max,step,s_fact,opt
      parameter(s_fact=5d-1)
      character*(name_len) name

 1    continue
      j=next_vary(name,name_len,c_min,c_max,step,slope,opt)
      if(j.ne.0)  then
        call set_variable(name,vect(j))
        goto 1
      endif
      end
      subroutine mtlmdf(ncon,nvar,tol,calls,call_lim,vect,dvect,fun_vec,&
     &diag,w_ifjac,w_ipvt,w_qtf,w_iwa1,w_iwa2,w_iwa3,w_iwa4,xold)

      use matchfi
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   LMDIF command.                                                     *
! Attributes:                                                          *
!   ncon      (int)     # constraints                                  *
!   nvar      (int)     # variables                                    *
!   tol       (real)    Final tolerance for match.                     *
!   calls     (int)     current call count                             *
!   call_lim  (int)     current call limit                             *
!   vect      (real)    variable values                                *
!   dvect     (real)    variable steps                                 *
!   fun_vect  (real)    function values                                *
!   all other working spaces for lmdif                                 *
!----------------------------------------------------------------------*
      integer calls,call_lim,ncon,nvar,i,ipvt(nvar)
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision tol,vect(*),dvect(*),fun_vec(*),diag(*),         &
     &w_ifjac(*),w_qtf(*),w_iwa1(*),w_iwa2(*),w_iwa3(*),w_iwa4(*),one,  &
     &xold(*),w_ipvt(*)
      parameter(one=1d0)
      external mtfcn

      icovar = 0
      ilevel = 0
      ipvt(:) = 0
      call mtgeti(vect,dvect)
      call lmdif(mtfcn,ncon,nvar,calls,call_lim,vect,fun_vec,tol,diag,  &
     &one,w_ifjac,ncon,ipvt,w_qtf,w_iwa1,w_iwa2,w_iwa3,w_iwa4,xold)
      do i=1,nvar
         w_ipvt(i)=ipvt(i)
      enddo
 9999 end

      subroutine lmdif(fcn,m,n,calls,call_lim,x,fvec,epsfcn,diag,factor,&
     &fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4,xold)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   The purpose of LMDIF is to minimize the sum of the squares of      *
!   M nonlinear functions in N variables by a modification of          *
!   the Levenberg-Marquardt algorithm. The user must provide a         *
!   subroutine which calculates the functions. The Jacobian is         *
!   then calculated by a forward-difference approximation.             *
!                                                                      *
!       FCN is the name of the user-supplied subroutine which          *
!         calculates the functions. FCN must be declared               *
!         in an external statement in the user calling                 *
!         program, and should be written as follows:                   *
!                                                                      *
!         SUBROUTINE FCN(M,N,X,FVEC,IFLAG)                             *
!         DIMENSION X(N),FVEC(M)                                       *
!         CALCULATE THE FUNCTIONS AT X AND                             *
!         RETURN THIS VECTOR IN FVEC.                                  *
!         RETURN                                                       *
!         END                                                          *
!                                                                      *
!         The value of IFLAG should be set to zero, unless there       *
!         is an error in evaluation of the function.                   *
!                                                                      *
!       M is a positive integer input variable set to the number       *
!         of functions.                                                *
!                                                                      *
!       N is a positive integer input variable set to the number       *
!         of variables. N must not exceed M.                           *
!                                                                      *
!       X is an array of length N. On input X must contain             *
!         an initial estimate of the solution vector. On output X      *
!         contains the final estimate of the solution vector.          *
!                                                                      *
!       FVEC is an output array of length M which contains             *
!         the functions evaluated at the output X.                     *
!                                                                      *
!       EPSFCN is an input variable used in determining a suitable     *
!         step length for the forward-difference approximation. This   *
!         approximation assumes that the relative errors in the        *
!         functions are of the order of EPSFCN. If EPSFCN is less      *
!         than the machine precision, it is assumed that the relative  *
!         errors in the functions are of the order of the machine      *
!         precision.                                                   *
!                                                                      *
!       DIAG is an array of length N. If MODE = 1 (see                 *
!         below), DIAG is internally set. If MODE = 2, DIAG            *
!         must contain positive entries that serve as                  *
!         multiplicative scale factors for the variables.              *
!                                                                      *
!       FACTOR is a positive input variable used in determining the    *
!         initial step bound. This bound is set to the product of      *
!         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else  *
!         to FACTOR itself. In most cases FACTOR should lie in the     *
!         interval (.1,100.). 100. Is a generally recommended value.   *
!                                                                      *
!       FJAC is an output M by N array. The upper N by N submatrix     *
!         of FJAC contains an upper triangular matrix R with           *
!         diagonal elements of nonincreasing magnitude such that       *
!                                                                      *
!                T     T           T                                   *
!               P *(JAC *JAC)*P = R *R,                                *
!                                                                      *
!         where P is a permutation matrix and JAC is the final         *
!         calculated Jacobian. column J of P is column IPVT(J)         *
!         (see below) of the identity matrix. The lower trapezoidal    *
!         part of FJAC contains information generated during           *
!         the computation of R.                                        *
!                                                                      *
!       LDFJAC is a positive integer input variable not less than M    *
!         which specifies the leading dimension of the array FJAC.     *
!                                                                      *
!       IPVT is an integer output array of length N. IPVT              *
!         defines a permutation matrix P such that JAC*P = Q*r,        *
!         where JAC is the final calculated Jacobian, Q is             *
!         orthogonal (not stored), and R is upper triangular           *
!         with diagonal elements of nonincreasing magnitude.           *
!         column J of P is column IPVT(J) of the identity matrix.      *
!                                                                      *
!       QTF is an output array of length N which contains              *
!         the first N elements of the vector (Q transpose)*FVEC.       *
!                                                                      *
!       WA1, WA2, and WA3 are work arrays of length N.                 *
!                                                                      *
!       WA4 is a work array of length M.                               *
! Source:                                                              *
!   Argonne National Laboratory. MINPACK Project. March 1980.          *
!   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.             *
!----------------------------------------------------------------------*
      integer izero
      integer i,iflag,info,iter,j,l,ldfjac,level,m,n,calls,call_lim,    &
     &ipvt(n), ireset
      double precision fmin_old,actred,delta,dirder,epsfcn,factor,fnorm,&
     &fnorm1,ftol,gnorm,gtol,par,pnorm,prered,ratio,sum,temp,temp1,     &
     &temp2,vmod,xnorm,xtol,x(n),xold(n),fvec(m),diag(n),fjac(ldfjac,n),&
     &qtf(n),wa1(n),wa2(n),wa3(n),wa4(m),zero,one,two,p1,p25,p5,p75,p90,&
     &p0001,epsil,epsmch
      parameter(zero=0d0,one=1d0,two=2d0,p1=0.1d0,p5=0.5d0,p25=0.25d0,  &
     &p75=0.75d0,p90=0.9d0,p0001=0.0001d0,epsil=1d-8,epsmch=1d-16)
      external fcn, mtcond

      ireset = 0
      izero = 0
      info = 0
      level = 0
      ftol = epsfcn
      gtol = epsil
      xtol = epsil
!---- Compute matching functions.
      call mtcond(izero, m, fvec, iflag)
!      call fcn(m,n,x,fvec,iflag)
      calls = calls + 1
      if (iflag .ne. 0) then
        call fort_warn('LMDIF', ' stopped, possibly unstable')
        info = - 1
        go to 300
      endif
      fnorm = vmod(m, fvec)
      fmin = fnorm**2
      fmin_old = fmin
      edm = fmin
      write(*,831) fmin

!---- Check the input parameters for errors.
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m                     &
     &.or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero       &
     &.or. call_lim .le. 0 .or. factor .le. zero) go to 300

!---- Quit, when initial value is already OK.
      if (fmin .le. ftol) then
        info = 4
        do j = 1, n
          xold(j) = x(j)
        enddo
        go to 300
      endif
      if (ilevel .ge. 1) call mtprnt('old', n, x)

!---- Initialize Levenberg-Marquardt parameter and iteration count
      par = zero
      iter = 1

!---- Beginning of the outer loop.
   30 continue
!      write(*,*) 'outer ',calls,fmin,fmin_old

!---- Calculate the Jacobian matrix.
      call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,xtol,wa4)
      calls = calls + n
      if (iflag .ne. 0) then
        info = - 1
        go to 300
      endif

!---- Compute the QR factorization of the Jacobian.
      call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)

!---- On the first iteration scale according to the norms
!     of the columns of the initial Jacobian.
!     Calculate the norm of the scaled X
!     and initialize the step bound delta.
      if (iter .eq. 1) then
        do j = 1, n
          xold(j) = x(j)
          diag(j) = wa2(j)
          if (wa2(j) .eq. zero) diag(j) = one
          wa3(j) = diag(j)*x(j)
        enddo
        xnorm = vmod(n, wa3)
        delta = factor*xnorm
        if (delta .eq. zero) delta = factor
      endif

!---- Form (Q transpose)*FVEC and store the first N components in QTF.
      do i = 1, m
        wa4(i) = fvec(i)
      enddo
      do j = 1, n
        if (fjac(j,j) .ne. zero) then
          sum = zero
          do i = j, m
            sum = sum + fjac(i,j)*wa4(i)
          enddo
          temp = -sum/fjac(j,j)
          do i = j, m
            wa4(i) = wa4(i) + fjac(i,j)*temp
          enddo
        endif
        fjac(j,j) = wa1(j)
        qtf(j) = wa4(j)
      enddo

!---- Compute the norm of the scaled gradient.
      gnorm = zero
      if (fnorm .ne. zero) then
        do j = 1, n
          l = ipvt(j)
          if (wa2(l) .ne. zero) then
            sum = zero
            do i = 1, j
              sum = sum + fjac(i,j)*(qtf(i)/fnorm)
            enddo
            gnorm = max(gnorm,abs(sum/wa2(l)))
          endif
        enddo
      endif

!---- Test for convergence of the gradient norm.
      if (gnorm .le. gtol) info = 4
      if (info .ne. 0) go to 300

!---- Rescale if necessary.
      do j = 1, n
        diag(j) = max(diag(j),wa2(j))
      enddo

!---- Beginning of the inner loop.
  200 continue
!      write(*,*) 'inner ',calls,fmin,fmin_old
!      if (fmin.lt.p90*fmin_old) write(*,830) calls,fmin
      fmin_old = fmin

!---- Determine the Levenberg-Marquardt parameter.
      call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,         &
     &wa3,wa4)

!---- Store the direction P and X + P. Calculate the norm of P.
      do j = 1, n
        wa1(j) = -wa1(j)
        wa2(j) = x(j) + wa1(j)
        wa3(j) = diag(j)*wa1(j)
      enddo
      pnorm = vmod(n, wa3)

!---- On the first iteration, adjust the initial step bound.
      if (iter .eq. 1) delta = min(delta,pnorm)

!---- Evaluate the function at X + P and calculate its norm.
! 23.5.2002: Oliver Bruening: inserted a routine that checks for
!                       minimum and maximum values of the parameters:
      call mtlimit(wa2,ireset)
      call fcn(m,n,wa2,wa4,iflag)
      calls = calls + 1
      if (iflag .ne. 0) then
        fnorm1 = two * fnorm
      else
        fnorm1 = vmod(m, wa4)
      endif

!---- Compute the scaled actual reduction.
      actred = -one
      if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2

!---- Compute the scaled predicted reduction and
!     the scaled directional derivative.
      do j = 1, n
        wa3(j) = zero
        l = ipvt(j)
        temp = wa1(l)
        do i = 1, j
          wa3(i) = wa3(i) + fjac(i,j)*temp
        enddo
      enddo
      temp1 = vmod(n, wa3)/fnorm
      temp2 = (sqrt(par)*pnorm)/fnorm
      prered = temp1**2 + temp2**2/p5
      dirder = -(temp1**2 + temp2**2)

!---- Compute the ratio of the actual to the predicted reduction.
      ratio = zero
      if (prered .ne. zero) ratio = actred/prered

!---- Update the step bound.
      if (ratio .le. p25) then
        if (actred .ge. zero) temp = p5
        if (actred .lt. zero)                                           &
     &temp = p5*dirder/(dirder + p5*actred)
        if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
        delta = temp*min(delta,pnorm/p1)
        par = par/temp
      else if (par .eq. zero .or. ratio .ge. p75) then
        delta = pnorm/p5
        par = p5*par
      endif

!---- Test for successful iteration.
      if (ratio .ge. p0001) then

!---- Successful iteration. Update X, FVEC, and their norms.
        do j = 1, n
          x(j) = wa2(j)
          xold(j) = wa2(j)
          wa2(j) = diag(j)*x(j)
        enddo
        do i = 1, m
          fvec(i) = wa4(i)
        enddo
        xnorm = vmod(n, wa2)
        fnorm = fnorm1
        iter = iter + 1

!---- If requested, print iterates.
        fmin = fnorm**2
        edm = gnorm * fmin
        level = 3
        if (mod(iter,10) .eq. 0) level = 2
        if (ilevel .ge. level) call mtprnt('inter', n, x)
        write(*,830) calls,fmin
      endif

!---- Tests for convergence.
      if (abs(actred) .le. ftol .and. prered .le. ftol                  &
     &.and. p5*ratio .le. one) info = 1
      if (delta .le. xtol*xnorm) info = 3
      if (abs(actred) .le. ftol .and. prered .le. ftol                  &
     &.and. p5*ratio .le. one .and. info .eq. 2) info = 2
      if (fmin .le. ftol) info = 4
      if (info .ne. 0) go to 300

!---- Tests for termination and stringent tolerances.
      if (calls .ge. call_lim) info = 5
      if (abs(actred) .le. epsmch .and. prered .le. epsmch              &
     &.and. p5*ratio .le. one) info = 6
      if (delta .le. epsmch*xnorm) info = 7
      if (gnorm .le. epsmch) info = 8
      if (ireset .gt. 20) info = 10
      if (info .ne. 0) go to 300

!---- End of the inner loop. Repeat if iteration unsuccessful.
      if (ratio .lt. p0001) go to 200
      if(info.eq.10) goto 300

!---- End of the outer loop.
      go to 30

!---- Termination, either normal or user imposed.
  300 continue
      if (info .lt. 0) then
        print *, '++++++++++ LMDIF ended: unstable'
      else if (info .eq. 0) then
        print *, '++++++++++ LMDIF ended: error'
      else if (info .eq. 4) then
        print *, '++++++++++ LMDIF ended: converged successfully'
      else if (info .eq. 3) then
        print *, '++++++++++ LMDIF ended: converged without success'
      else if (info .lt. 3) then
        print *, '++++++++++ LMDIF ended: converged: info = ', info
      else if (info .eq. 5) then
        print *, '++++++++++ LMDIF ended: call limit'
      else if (info .eq. 10) then
        print *, '++++++++++'
        print *, '++++++++++ LMDIF ended: variables too close to limit'
        print *, '++++++++++'
      else
        print *, '++++++++++ LMDIF ended: accuracy limit'
      endif

      call fcn(m,n,xold,fvec,iflag)
      fnorm = vmod(m, fvec)
      fmin = fnorm**2
!      call mtputi(xold)
      write(*,830) calls,fmin

      if (ilevel .ge. 1) call mtprnt('last',n, x)

  831 format('Initial Penalty Function = ',e16.8,//)
  830 format('call:',I8,3x,'Penalty function = ',e16.8)
      end
      subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   This subroutine computes a forward-difference approximation        *
!   to the M by N Jacobian matrix associated with a specified          *
!   problem of M functions in N variables.                             *
! Input:                                                               *
!       FCN is the name of the user-supplied subroutine which          *
!         calculates the functions. FCN must be declared               *
!         in an external statement in the user calling                 *
!         program, and should be written as follows:                   *
!                                                                      *
!         SUBROUTINE FCN(M,N,X,FVEC,IFLAG)                             *
!         DIMENSION X(N),FVEC(M)                                       *
!         CALCULATE THE FUNCTIONS AT X AND                             *
!         RETURN THIS VECTOR IN FVEC.                                  *
!         RETURN                                                       *
!         END                                                          *
!                                                                      *
!         The value of IFLAG should be set to zero, unless there       *
!         is an error in evaluation of the function.                   *
!                                                                      *
!       M is a positive integer input variable set to the number       *
!         of functions.                                                *
!                                                                      *
!       N is a positive integer input variable set to the number       *
!         of variables. N must not exceed M.                           *
!                                                                      *
!       X is an input array of length N.                               *
!                                                                      *
!       FVEC is an input array of length M which must contain the      *
!         functions evaluated at X.                                    *
!                                                                      *
!       FJAC is an output M by N array which contains the              *
!         approximation to the Jacobian matrix evaluated at X.         *
!                                                                      *
!       LDFJAC is a positive integer input variable not less than M    *
!         which specifies the leading dimension of the array FJAC.     *
!                                                                      *
!       IFLAG is an integer variable which tells the calling program   *
!         wether the approximation is valid.                           *
!                                                                      *
!       EPSFCN is an input variable used in determining a suitable     *
!         step length for the forward-difference approximation. This   *
!         approximation assumes that the relative errors in the        *
!         functions are of the order of EPSFCN. If EPSFCN is less      *
!         than the machine precision, it is assumed that the relative  *
!         errors in the functions are of the order of the machine      *
!         precision.                                                   *
!                                                                      *
!       WA is a work array of length M.                                *
! Source:                                                              *
!   Argonne National Laboratory. MINPACK Project. March 1980.          *
!   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.             *
!----------------------------------------------------------------------*
      integer i,iflag,j,ldfjac,m,n
      double precision eps,epsfcn,fjac(ldfjac,n),fvec(m),h,temp,wa(m),  &
     &x(n),zero,epsmch
      parameter(zero=0d0,epsmch=1d-16)
      external fcn

      eps = sqrt(max(epsfcn,epsmch))
      iflag = 0

      do j = 1, n
        temp = x(j)
        h = eps*abs(temp)
        if (h .eq. zero) h = eps
        x(j) = temp + h
        call fcn(m,n,x,wa,iflag)
        x(j) = temp
        if (iflag .ne. 0) go to 30
        do i = 1, m
          fjac(i,j) = (wa(i) - fvec(i))/h
        enddo
      enddo
   30 continue

      end
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   This subroutine uses Householder transformations with column       *
!   pivoting (optional) to compute a QR factorization of the           *
!   M by N matrix a. That is, QRFAC determines an orthogonal           *
!   matrix Q, a permutation matrix P, and an upper trapezoidal         *
!   matrix R with diagonal elements of nonincreasing magnitude,        *
!   such that A*P = Q*R. The Householder transformation for            *
!   column K, K = 1,2,...,min(M,N), is of the form                     *
!                                                                      *
!                           T                                          *
!           I - (1/U(K))*U*U                                           *
!                                                                      *
!   where U has zeros in the first K-1 positions. The form of          *
!   this transformation and the method of pivoting first               *
!   appeared in the corresponding LINPACK subroutine.                  *
!                                                                      *
!       M is a positive integer input variable set to the number       *
!         of rows of A.                                                *
!                                                                      *
!       N is a positive integer input variable set to the number       *
!         of columns of A.                                             *
!                                                                      *
!       A is an M by N array. On input A contains the matrix for       *
!         which the QR factorization is to be computed. On output      *
!         The strict upper trapezoidal part of A contains the strict   *
!         upper trapezoidal part of R, and the lower trapezoidal       *
!         part of A contains a factored form of Q (the non-trivial     *
!         elements of the U vectors described above).                  *
!                                                                      *
!       LDA is a positive integer input variable not less than M       *
!         which specifies the leading dimension of the array A.        *
!                                                                      *
!       PIVOT is a logical input variable. If PIVOT is set true,       *
!         then column pivoting is enforced. If PIVOT is set false,     *
!         then no column pivoting is done.                             *
!                                                                      *
!       IPVT is an integer output array of length LIPVT. Ipvt          *
!         defines the permutation matrix P such that a*p = Q*r.        *
!         column J of P is column IPVT(J) of the identity matrix.      *
!         if PIVOT is false, IPVT is not referenced.                   *
!                                                                      *
!       LIPVT is a positive integer input variable. If PIVOT is false, *
!         then LIPVT may be as small as 1. If PIVOT is true, then      *
!         LIPVT must be at least N.                                    *
!                                                                      *
!       RDIAG is an output array of length N which contains the        *
!         diagonal elements of R.                                      *
!                                                                      *
!       ACNORM is an output array of length N which contains the       *
!         norms of the corresponding columns of the input matrix A.    *
!         If this information is not needed, then ACNORM can coincide  *
!         with RDIAG.                                                  *
!                                                                      *
!       WA is a work array of length N. If PIVOT is false, then WA     *
!         can coincide with RDIAG.                                     *
! Source:                                                              *
!   Argonne National Laboratory. MINPACK Project. March 1980.          *
!   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.             *
!----------------------------------------------------------------------*
      logical pivot
      integer i,j,k,kmax,lda,lipvt,m,minmn,n,ipvt(lipvt)
      double precision ajnorm,sum,temp,vmod,a(lda,n),rdiag(n),acnorm(n),&
     &wa(n),zero,one,p05,epsmch
      parameter(zero=0d0,one=1d0,p05=0.05d0,epsmch=1d-16)

!---- Compute the initial column norms and initialize several arrays.
      do j = 1, n
        acnorm(j) = vmod(m, a(1,j))
        rdiag(j) = acnorm(j)
        wa(j) = rdiag(j)
        if (pivot) ipvt(j) = j
      enddo
!---- Reduce A to R with Householder transformations.
      minmn = min(m,n)
      do j = 1, minmn
        if (pivot) then

!---- Bring the column of largest norm into the pivot position.
          kmax = j
          do k = j, n
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k
          enddo
          if (kmax .ne. j) then
            do i = 1, m
              temp = a(i,j)
              a(i,j) = a(i,kmax)
              a(i,kmax) = temp
            enddo
            rdiag(kmax) = rdiag(j)
            wa(kmax) = wa(j)
            k = ipvt(j)
            ipvt(j) = ipvt(kmax)
            ipvt(kmax) = k
          endif
        endif

!---- Compute the Householder transformation to reduce the
!     J-th column of A to a multiple of the J-th unit vector.
        ajnorm = vmod(m-j+1, a(j,j))
        if (ajnorm .ne. zero) then
          if (a(j,j) .lt. zero) ajnorm = -ajnorm
          do i = j, m
            a(i,j) = a(i,j)/ajnorm
          enddo
          a(j,j) = a(j,j) + one

!---- Apply the transformation to the remaining columns
!     and update the norms.
          do k = j + 1, n
            sum = zero
            do i = j, m
              sum = sum + a(i,j)*a(i,k)
            enddo
            temp = sum/a(j,j)
            do i = j, m
              a(i,k) = a(i,k) - temp*a(i,j)
            enddo
            if (pivot .and. rdiag(k) .ne. zero) then
              temp = a(j,k)/rdiag(k)
              rdiag(k) = rdiag(k)*sqrt(max(zero,one-temp**2))
              if (p05*(rdiag(k)/wa(k))**2 .le. epsmch) then
                rdiag(k) = vmod(m-j, a(j+1,k))
                wa(k) = rdiag(k)
              endif
            endif
          enddo
        endif
        rdiag(j) = -ajnorm
      enddo

      end
      subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)


      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Given an M by N matrix A, an N by N nonsingular diagonal           *
!   matrix D, an M-vector B, and a positive number DELTA,              *
!   the problem is to determine a value for the parameter              *
!   PAR such that if X solves the system                               *
!                                                                      *
!           A*X = B,     SQRT(PAR)*D*X = 0,                            *
!                                                                      *
!   in the least squares sense, and DXNORM is the Euclidean            *
!   norm of D*X, then either PAR is zero and                           *
!                                                                      *
!           (DXNORM-DELTA) .LE. 0.1*DELTA,                             *
!                                                                      *
!   or PAR is positive and                                             *
!                                                                      *
!           ABS(DXNORM-DELTA) .LE. 0.1*DELTA.                          *
!                                                                      *
!   This subroutine completes the solution of the problem              *
!   if it is provided with the necessary information from the          *
!   QR factorization, with column pivoting, of A. That is, if          *
!   A*P = Q*R, where P is a permutation matrix, Q has orthogonal       *
!   columns, and R is an upper triangular matrix with diagonal         *
!   elements of nonincreasing magnitude, then LMPAR expects            *
!   the full upper triangle of R, the permutation matrix P,            *
!   and the first N components of (Q transpose)*B. On output           *
!   LMPAR also provides an upper triangular matrix S such that         *
!                                                                      *
!            T   T                   T                                 *
!           P *(A *A + PAR*D*D)*P = S *S.                              *
!                                                                      *
!   S is employed within LMPAR and may be of separate interest.        *
!                                                                      *
!   Only a few iterations are generally needed for convergence         *
!   of the algorithm. If, however, the limit of 10 iterations          *
!   is reached, then the output PAR will contain the best              *
!   value obtained so far.                                             *
!                                                                      *
!       N is a positive integer input variable set to the order of R.  *
!                                                                      *
!       R is an N by N array. On input the full upper triangle         *
!         must contain the full upper triangle of the matrix R.        *
!         on output the full upper triangle is unaltered, and the      *
!         strict lower triangle contains the strict upper triangle     *
!         (transposed) of the upper triangular matrix S.               *
!                                                                      *
!       LDR is a positive integer input variable not less than N       *
!         which specifies the leading dimension of the array R.        *
!                                                                      *
!       IPVT is an integer input array of length N which defines the   *
!         permutation matrix P such that A*P = Q*R. column J of P      *
!         is column IPVT(J) of the identity matrix.                    *
!                                                                      *
!       DIAG is an input array of length N which must contain the      *
!         diagonal elements of the matrix D.                           *
!                                                                      *
!       QTB is an input array of length N which must contain the first *
!         N elements of the vector (Q transpose)*B.                    *
!                                                                      *
!       DELTA is a positive input variable which specifies an upper    *
!         bound on the Euclidean norm of D*X.                          *
!                                                                      *
!       PAR is a nonnegative variable. On input PAR contains an        *
!         initial estimate of the Levenberg-Marquardt parameter.       *
!         on output PAR contains the final estimate.                   *
!                                                                      *
!       X is an output array of length N which contains the least      *
!         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0,   *
!         for the output PAR.                                          *
!                                                                      *
!       SDIAG is an output array of length N which contains the        *
!         diagonal elements of the upper triangular matrix S.          *
!                                                                      *
!       WA1 and WA2 are work arrays of length N.                       *
! Source:                                                              *
!   Argonne National Laboratory. MINPACK Project. March 1980.          *
!   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.             *
!----------------------------------------------------------------------*
      integer i,iter,j,k,l,ldr,n,nsing,ipvt(n)
      double precision delta,diag(n),dxnorm,fp,gnorm,par,parc,parl,paru,&
     &qtb(n),r(ldr,n),sdiag(n),sum,temp,vmod,wa1(n),wa2(n),x(n),zero,   &
     &fltmin,p1,p001
      parameter(zero=0d0,fltmin=1d-35,p1=0.1d0,p001=0.001d0)

!---- Compute and store in X the Gauss-Newton direction. If the
!     Jacobian is rank-deficient, obtain a least squares solution.
      nsing = n
      do j = 1, n
        wa1(j) = qtb(j)
        if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1
        if (nsing .lt. n) wa1(j) = zero
      enddo
      do k = 1, nsing
        j = nsing - k + 1
        wa1(j) = wa1(j)/r(j,j)
        temp = wa1(j)
        do i = 1, j - 1
          wa1(i) = wa1(i) - r(i,j)*temp
        enddo
      enddo
      do j = 1, n
        l = ipvt(j)
        x(l) = wa1(j)
      enddo

!---- Initialize the iteration counter.
!     Evaluate the function at the origin, and test
!     for acceptance of the Gauss-Newton direction.
      iter = 0
      do j = 1, n
        wa2(j) = diag(j)*x(j)
      enddo
      dxnorm = vmod(n, wa2)
      fp = dxnorm - delta
      if (fp .le. p1*delta) go to 220

!---- If the Jacobian is not rank deficient, the Newton
!     step provides a lower bound, PARL, for the zero of
!     the function. Otherwise set this bound to zero.
      parl = zero
      if (nsing .ge. n) then
        do j = 1, n
          l = ipvt(j)
          wa1(j) = diag(l)*(wa2(l)/dxnorm)
        enddo
        do j = 1, n
          sum = zero
          do i = 1, j - 1
            sum = sum + r(i,j)*wa1(i)
          enddo
          wa1(j) = (wa1(j) - sum)/r(j,j)
        enddo
        temp = vmod(n, wa1)
        parl = ((fp/delta)/temp)/temp
      endif

!---- Calculate an upper bound, PARU, for the zero of the function.
      do j = 1, n
        sum = zero
        do i = 1, j
          sum = sum + r(i,j)*qtb(i)
        enddo
        l = ipvt(j)
        wa1(j) = sum/diag(l)
      enddo
      gnorm = vmod(n, wa1)
      paru = gnorm/delta
      if (paru .eq. zero) paru = fltmin/min(delta,p1)

!---- If the input PAR lies outside of the interval (PARL,PARU),
!     set PAR to the closer endpoint.
      par = max(par,parl)
      par = min(par,paru)
      if (par .eq. zero) par = gnorm/dxnorm

!---- Beginning of an iteration.
  150 continue
      iter = iter + 1

!---- Evaluate the function at the current value of PAR.
      if (par .eq. zero) par = max(fltmin,p001*paru)
      temp = sqrt(par)
      do j = 1, n
        wa1(j) = temp*diag(j)
      enddo
      call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
      do j = 1, n
        wa2(j) = diag(j)*x(j)
      enddo
      dxnorm = vmod(n, wa2)
      temp = fp
      fp = dxnorm - delta

!---- If the function is small enough, accept the current value
!     of PAR. also test for the exceptional cases where PARL
!     is zero or the number of iterations has reached 10.
      if (abs(fp) .le. p1*delta                                         &
     &.or. parl .eq. zero .and. fp .le. temp                            &
     &.and. temp .lt. zero .or. iter .eq. 10) go to 220

!---- Compute the Newton correction.
      do j = 1, n
        l = ipvt(j)
        wa1(j) = diag(l)*(wa2(l)/dxnorm)
      enddo
      do j = 1, n
        wa1(j) = wa1(j)/sdiag(j)
        temp = wa1(j)
        do i = j + 1, n
          wa1(i) = wa1(i) - r(i,j)*temp
        enddo
      enddo
      temp = vmod(n, wa1)
      parc = ((fp/delta)/temp)/temp

!---- Depending on the sign of the function, update PARL or PARU.
      if (fp .gt. zero) parl = max(parl,par)
      if (fp .lt. zero) paru = min(paru,par)

!---- Compute an improved estimate for PAR.
      par = max(parl,par+parc)

!---- End of an iteration.
      go to 150
  220 continue

!---- Termination.
      if (iter .eq. 0) par = zero

      end
      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Given an M by N matrix A, an N by N diagonal matrix D,             *
!   and an M-vector B, the problem is to determine an X which          *
!   solves the system                                                  *
!                                                                      *
!           A*X = B,     D*X = 0,                                      *
!                                                                      *
!   in the least squares sense.                                        *
!                                                                      *
!   This subroutine completes the solution of the problem              *
!   if it is provided with the necessary information from the          *
!   QR factorization, with column pivoting, of A. That is, if          *
!   A*P = Q*R, where P is a permutation matrix, Q has orthogonal       *
!   columns, and R is an upper triangular matrix with diagonal         *
!   elements of nonincreasing magnitude, then QRSOLV expects           *
!   the full upper triangle of R, the permutation matrix P,            *
!   and the first N components of (Q transpose)*B. The system          *
!   A*X = B, D*X = 0, is then equivalent to                            *
!                                                                      *
!                  T      T                                            *
!           R*Z = Q *B,  P *D*P*Z = 0,                                 *
!                                                                      *
!   where X = P*Z. If this system does not have full rank,             *
!   then a least squares solution is obtained. On output QRSOLV        *
!   also provides an upper triangular matrix S such that               *
!                                                                      *
!            T   T               T                                     *
!           P *(A *A + D*D)*P = S *S.                                  *
!                                                                      *
!     S is computed within QRSOLV and may be of separate interest.     *
!                                                                      *
!       N is a positive integer input variable set to the order of R.  *
!                                                                      *
!       R is an N by N array. On input the full upper triangle         *
!         must contain the full upper triangle of the matrix R.        *
!         On output the full upper triangle is unaltered, and the      *
!         strict lower triangle contains the strict upper triangle     *
!         (transposed) of the upper triangular matrix S.               *
!                                                                      *
!       LDR is a positive integer input variable not less than N       *
!         which specifies the leading dimension of the array R.        *
!                                                                      *
!       IPVT is an integer input array of length N which defines the   *
!         permutation matrix P such that A*P = Q*R. Column J of P      *
!         is column IPVT(J) of the identity matrix.                    *
!                                                                      *
!       DIAG is an input array of length N which must contain the      *
!         diagonal elements of the matrix D.                           *
!                                                                      *
!       QTB is an input array of length N which must contain the first *
!         N elements of the vector (Q transpose)*B.                    *
!                                                                      *
!       X is an output array of length N which contains the least      *
!         squares solution of the system A*X = B, D*X = 0.             *
!                                                                      *
!       SDIAG is an output array of length N which contains the        *
!         diagonal elements of the upper triangular matrix S.          *
!                                                                      *
!       WA is a work array of length N.                                *
! Source:                                                              *
!   Argonne National Laboratory. MINPACK Project. March 1980.          *
!   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.             *
!----------------------------------------------------------------------*
      integer i,j,k,l,ldr,n,nsing,ipvt(n)
      double precision cos,cotan,diag(n),qtb(n),qtbpj,r(ldr,n),sdiag(n),&
     &sin,sum,tan,temp,wa(n),x(n),zero,p5,p25
      parameter(zero=0d0,p5=0.5d0,p25=0.25d0)

!---- Copy R and (Q transpose)*B to preserve input and initialize S.
!     In particular, save the diagonal elements of R in X.
      do j = 1, n
        do i = j, n
          r(i,j) = r(j,i)
        enddo
        x(j) = r(j,j)
        wa(j) = qtb(j)
      enddo

!---- Eliminate the diagonal matrix D using a Givens rotation.
      do j = 1, n

!---- Prepare the row of D to be eliminated, locating the
!     diagonal element using P from the QR factorization.
        l = ipvt(j)
        if (diag(l) .ne. zero) then
          do k = j, n
            sdiag(k) = zero
          enddo
          sdiag(j) = diag(l)

!---- The transformations to eliminate the row of D
!     modify only a single element of (Q transpose)*B
!     beyond the first N, which is initially zero.
          qtbpj = zero
          do k = j, n

!---- Determine a Givens rotation which eliminates the
!     appropriate element in the current row of D.
            if (sdiag(k) .ne. zero) then
              if (abs(r(k,k)) .lt. abs(sdiag(k))) then
                cotan = r(k,k)/sdiag(k)
                sin = p5/sqrt(p25+p25*cotan**2)
                cos = sin*cotan
              else
                tan = sdiag(k)/r(k,k)
                cos = p5/sqrt(p25+p25*tan**2)
                sin = cos*tan
              endif

!---- Compute the modified diagonal element of R and
!     the modified element of ((Q transpose)*b,0).
              r(k,k) = cos*r(k,k) + sin*sdiag(k)
              temp = cos*wa(k) + sin*qtbpj
              qtbpj = -sin*wa(k) + cos*qtbpj
              wa(k) = temp

!---- Accumulate the tranformation in the row of S.
              do i = k + 1, n
                temp = cos*r(i,k) + sin*sdiag(i)
                sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
                r(i,k) = temp
              enddo
            endif
          enddo
        endif

!---- Store the diagonal element of S and restore
!     the corresponding diagonal element of R.
        sdiag(j) = r(j,j)
        r(j,j) = x(j)
      enddo

!---- Solve the triangular system for z. If the system is
!     singular, then obtain a least squares solution.
      nsing = n
      do j = 1, n
        if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
        if (nsing .lt. n) wa(j) = zero
      enddo
      do j = nsing, 1, - 1
        sum = zero
        do i = j + 1, nsing
          sum = sum + r(i,j)*wa(i)
        enddo
        wa(j) = (wa(j) - sum)/sdiag(j)
      enddo

!---- Permute the components of Z back to components of X.
      do j = 1, n
        l = ipvt(j)
        x(l) = wa(j)
      enddo

      end
      subroutine mtprnt(text,n,x)
      implicit none
      integer n
      double precision x(n)
      character*(*) text

      print *, text, ' variable values: ', x
      end
      subroutine mtmigr(ncon,nvar,strategy,tol,calls,call_lim,vect,     &
     &dvect,fun_vect,w_iwa1,w_iwa2,w_iwa3,w_iwa4,w_iwa5,w_iwa6,w_iwa7,  &
     &w_iwa8)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   MIGRAD command.                                                    *
!   ncon      (int)     # constraints                                  *
!   nvar      (int)     # variables                                    *
!   strategy  (int)     strategy selection (see minuit manual).        *
!   tol       (real)    Final tolerance for match.                     *
!   calls     (int)     current call count                             *
!   call_lim  (int)     current call limit                             *
!   vect      (real)    variable values                                *
!   dvect     (real)    variable steps                                 *
!   fun_vect  (real)    function values                                *
!   all other working spaces for mtmig1                                *
!----------------------------------------------------------------------*
      integer calls,call_lim,ncon,nvar,strategy
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision tol,vect(*),dvect(*),fun_vect(*),w_iwa1(*),      &
     &w_iwa2(*),w_iwa3(*),w_iwa4(*),w_iwa5(*),w_iwa6(*),w_iwa7(*),      &
     &w_iwa8(*)
      external mtfcn

      icovar = 0
      ilevel = 0
!---- Too many variable parameters?
      if (nvar .gt. ncon) &
        call fort_warn('MTMIGR', &
          'More variables than constraints seen. MIGRAD may not converge to optimal solution.')

!---- Call minimization routine.
      call mtgeti(vect, dvect)
      call mtmig1(mtfcn, ncon, nvar,  strategy, tol, calls, call_lim,   &
     &vect, dvect, fun_vect, w_iwa1, w_iwa2, w_iwa3, w_iwa4, w_iwa5,    &
     &w_iwa6, w_iwa7, w_iwa8)

 9999 end
      subroutine mtmig1(fcn,nf,nx,strategy,tol,calls,call_lim,x,dx,fvec,&
     &covar, wa,work_1,work_2,work_3,work_4,work_5,work_6)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Minimization by MIGRAD method by Davidon/Fletcher/Powell.          *
!   (Computer Journal 13, 317 (1970).                                  *
! Input:                                                               *
!   FCN       (subr)    Returns value of penalty function.             *
!   NF        (integer) Number of functions.                           *
!   NX        (integer) Number of parameters.                          *
!   X(NX)     (real)    Parameter values. On output, best estimate.    *
!   DX(NX)    (real)    Parameter errors. On output, error estimate.   *
! Output:                                                              *
!   FVEC(NF)  (real)    Vector of function values in best point.       *
! Working arrays:                                                      *
!   COVAR(NX,NX)        Covariance matrix.                             *
!   WA(NX,7)            Working vectors.                               *
!----------------------------------------------------------------------*
      logical eflag
      integer i,iflag,improv,iter,j,level,nf,npsdf,nrstrt,nx,strategy,  &
     &calls,call_lim,mgrd,mg2,mvg,mflnu,mgsave,mxsave
      parameter(mgrd=1,mg2=2,mvg=3,mflnu=5,mgsave=6,mxsave=7)
! ilevel: print level
      double precision covar(nx,nx),d,delgam,dgi,dx(nx),fvec(nf),gdel,  &
     &gssq,gvg,sum,vdot,vgi,wa(nx,7),x(nx),tol,work_1(nx),work_2(nx),   &
     &work_3(nx),work_4(*),work_5(*),work_6(*),zero,one,two,half,epsmch,&
     &eps1,eps2
      parameter(zero=0d0,one=1d0,two=2d0,half=0.5d0,epsmch=1d-16,       &
     &eps1=1d-3,eps2=1d-4)
      external fcn

!---- Initialize penalty function.
      call fcn(nf, nx, x, fvec, iflag)
      calls = calls + 1
      if (iflag .ne. 0) then
        call fort_warn('MTMIG1','Matching stopped -- start point seems to be unstable')
        go to 500
      endif
      fmin = vdot(nf, fvec, fvec)
      edm = fmin

!---- Start MIGRAD algorithm.
      nrstrt = 0
      npsdf = 0

!---- Come here to restart algorithm.
  100 continue
      if (strategy .eq. 2  .or.  strategy.gt.2 .and. icovar.lt.2) then
        call mthess(fcn, nf, nx, calls, covar, x,                       &
     &wa(1,mgrd), wa(1,mg2), fvec, wa(1,mvg), work_1, work_2, work_3,   &
     &work_4, work_5,work_6)
        npsdf = 0
      else
        call mtderi(fcn, nf, nx, calls, x, wa(1,mgrd), wa(1,mg2), fvec)
        if (icovar .lt. 2) then
          do i = 1, nx
            do j = 1, nx
              covar(i,j) = zero
            enddo
            if (wa(i,mg2) .eq. zero) wa(i,mg2) = one
            covar(i,i) = one / wa(i,mg2)
          enddo
        endif
      endif

!---- Initialize for first iteration.
      improv = 0
      edm = zero
      do i = 1, nx
        sum = zero
        do j = 1, nx
          sum = sum + covar(i,j) * wa(j,mgrd)
        enddo
        edm = edm + sum * wa(i,mgrd)
      enddo
      edm = min(half * edm, fmin)

!---- Print after initialization.
      if (ilevel .ge. 1) then
        call mtprnt('init', nx, x)
      endif
      iter = 0

!==== Start main iteration loop: Check for call limit.
  200 if (calls .lt. call_lim) then

!---- Find step size according to Newton's method.
        gdel = zero
        gssq = zero
        do i = 1, nx
          sum = zero
          wa(i,mgsave) = wa(i,mgrd)
          gssq = gssq + wa(i,mgrd)**2
          do j = 1, nx
            sum = sum + covar(i,j) * wa(j,mgrd)
          enddo
          dx(i) = - sum
          gdel = gdel + dx(i) * wa(i,mgrd)
        enddo

!---- First derivatives all zero?
        if (gssq .eq. zero) go to 400

!---- If GDEL .GE. 0 matrix is not positive definite.
        if (gdel .ge. zero) then
          if (npsdf .eq. 0) then
            call symsol(covar, nx, eflag, work_1, work_2, work_3)
            call mtpsdf(covar,nx,work_4,work_5,work_6)
            call symsol(covar, nx, eflag, work_1, work_2, work_3)
            npsdf = 1
            go to 200
          else
            nrstrt = nrstrt + 1
            if (nrstrt .gt. strategy) go to 500
            go to 100
          endif
        endif

!---- Search for minimum along predicted line.
        call mtline(fcn, nf, nx, calls, x, dx, fvec,                    &
     &wa(1,mxsave), iflag)

!---- No improvement found.
        if (iflag .ne. 0) then
          if (edm .lt. eps1 * tol) go to 400
          print *, 'accuracy limit'
          if (edm .lt. two * epsmch * fmin) go to 500
          if (strategy .eq. 0  .and.  nrstrt .eq. 0) then
            strategy = 1
            nrstrt = 1
            go to 100
          endif
          print *, 'failed'
          go to 500
        endif

!---- Find gradient in new point.
        call mtderi(fcn, nf, nx, calls, x, wa(1,mgrd), wa(1,mg2), fvec)
        npsdf = 0

!---- Estimated distance to minimum.
  300   continue
        edm = zero
        gvg = zero
        delgam = zero
        do i = 1, nx
          vgi = zero
          sum = zero
          do j = 1, nx
            vgi = vgi + covar(i,j) * (wa(j,mgrd) - wa(j,mgsave))
            sum = sum + covar(i,j) * wa(j,mgrd)
          enddo
          wa(i,mvg) = vgi
          dgi = wa(i,mgrd) - wa(i,mgsave)
          gvg = gvg + vgi * dgi
          delgam = delgam + dx(i) * dgi
          edm = edm + sum * wa(i,mgrd)
        enddo
        edm = min(half * edm, fmin)

!---- Test for convergence and print-out.
        if (edm .ge. zero  .and.  edm .lt. eps2 * tol) go to 400
        iter = iter + 1
        level = 3
        if (mod(iter,10) .eq. 0) level = 2
        if (ilevel .ge. level) then
          call mtprnt('progress', nx, x)
        endif
        write(*,830) calls,fmin

!---- Force positive definiteness.
        if (edm .lt. zero  .or. gvg .le. zero) then
          icovar = 0
          if (npsdf .eq. 1) go to 500
          call symsol(covar, nx, eflag, work_1, work_2, work_3)
          call mtpsdf(covar,nx,work_4,work_5,work_6)
          call symsol(covar, nx, eflag, work_1, work_2, work_3)
          npsdf = 1
          go to 300
        endif

!---- Update covariance matrix.
        do i = 1, nx
          do j = 1, nx
            d = dx(i) * dx(j) / delgam - wa(i,mvg) * wa(j,mvg) / gvg
            covar(i,j) = covar(i,j) + d
          enddo
        enddo
        if (delgam .gt. gvg) then
          do i = 1, nx
            wa(i,mflnu) = dx(i) / delgam - wa(i,mvg) / gvg
          enddo
          do i = 1, nx
            do j = 1, nx
              d = gvg * wa(i,mflnu) * wa(j,mflnu)
              covar(i,j) = covar(i,j) + d + d
            enddo
          enddo
        endif
        improv = improv + 1
        if (improv .ge. nx) icovar = 3
        go to 200
      endif

!---- Call limit reached.
      print *, 'call limit'
      go to 500

!==== End of main iteration loop; Check covariance matrix.
  400 continue
      if (strategy .ge. 2 .or. (strategy.eq.1 .and. icovar.lt.3)) then
        print *, 'verify'
        call mthess(fcn, nf, nx, calls, covar, x,                       &
     &wa(1,mgrd), wa(1,mg2), fvec, wa(1,mvg), work_1, work_2, work_3,   &
     &work_4, work_5,work_6)
        npsdf = 0
        if (edm .gt. eps1 * tol) go to 100
      endif
      print *, 'converged'

!---- Common exit point; final print-out.
  500 continue
      call mtputi(x)
      if (ilevel .ge. 1) call mtprnt('final', nx, x)

      write(*,830) calls,fmin

  830 format('call:',I8,3x,'Penalty function = ',e16.8)
      end
      subroutine mthess(fcn,nf,nx,calls,covar,x,grd,g2,fvec,wa,work_1,  &
     &work_2,work_3,work_4,work_5,work_6)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Build covariance matrix.                                           *
! input:                                                               *
!   fcn       (subr)    returns value of penalty function.             *
!   nf        (integer) number of functions.                           *
!   nx        (integer) number of parameters.                          *
!   calls     (integer) number of calls (increased)                    *
!   x(nx)     (real)    parameter values. on output, best estimate.    *
! output:                                                              *
!   covar(nx,nx)        covariance matrix.                             *
!   grd(nx)   (real)    gradient of penalty function                   *
!                       w.r.t. internal parameter values.              *
!   g2(nx)    (real)    second derivatives of penalty function         *
!                       w.r.t. internal parameter values.              *
! working array:                                                       *
!   fvec(nf)  (real)    function values.                               *
!   wa(nx,2)  (real)    working vectors.                               *
!----------------------------------------------------------------------*
      logical eflag
      integer i,icycle,iflag,j,nf,nx,calls
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision eps,f1,f2,fij,vdot,xs1,xs2,xsave,xstep,          &
     &covar(nx,nx),x(nx),grd(nx),g2(nx),wa(nx,2),fvec(nf),work_1(nx),   &
     &work_2(nx),work_3(nx),work_4(*),work_5(*),work_6(*),zero,one,two, &
     &half,epsmch
      parameter(zero=0d0,one=1d0,two=2d0,half=0.5d0,epsmch=1d-16)
      external fcn

      eps = sqrt(epsmch)

      do i = 1, nx
        xsave = x(i)
        xstep = eps * max(abs(xsave), one)
        do icycle = 1, 10
          x(i) = xsave + xstep
          call fcn(nf, nx, x, fvec, iflag)
          calls = calls + 1
          if (iflag .eq. 0) then
            f2 = vdot(nf, fvec, fvec)
            x(i) = xsave - xstep
            call fcn(nf, nx, x, fvec, iflag)
            calls = calls + 1
            if (iflag .eq. 0) then
              f1 = vdot(nf, fvec, fvec)
              go to 40
            endif
          endif
          xstep = half * xstep
        enddo
        f1 = fmin
        f2 = fmin
   40   continue
        grd(i) = (f2 - f1) / (two * xstep)
        g2(i) = (f2 - two * fmin + f1) / xstep**2
        if (g2(i) .eq. zero) g2(i) = one
        x(i) = xsave
        covar(i,i) = g2(i)
        wa(i,1) = f2
        wa(i,2) = xstep
      enddo

!---- Off-diagonal elements.
      do i = 1, nx - 1
        xs1 = x(i)
        x(i) = xs1 + wa(i,2)
        do j = i + 1, nx
          xs2 = x(j)
          x(j) = xs2 + wa(j,2)
          call fcn(nf, nx, x, fvec, iflag)
          calls = calls + 1
          if (iflag .eq. 0) then
            fij = vdot(nf, fvec, fvec)
            covar(i,j) = (fij+fmin-wa(i,1)-wa(j,1)) / (wa(i,2)*wa(j,2))
            covar(j,i) = covar(i,j)
          else
            covar(i,j) = zero
            covar(j,i) = zero
          endif
          x(j) = xs2
        enddo
        x(i) = xs1
      enddo

!---- Restore original point.
      call mtputi(x)

!---- Ensure positive definiteness and invert.
      call mtpsdf(covar,nx,work_4,work_5,work_6)
      call symsol(covar, nx, eflag, work_1, work_2, work_3)

      end
      subroutine mtderi(fcn,nf,nx,calls,x,grd,g2,fvec)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Find first derivatives of penalty function.                        *
! Input:                                                               *
!   FCN       (subr)    Returns value of penalty function.             *
!   NF        (integer) Number of functions.                           *
!   NX        (integer) Number of parameters.                          *
!   X(NX)     (real)    Parameter values. On output, best estimate.    *
! Output:                                                              *
!   GRD(*)    (real)    Gradient of penalty function                   *
!                       w.r.t. internal parameter values.              *
!   G2(*)     (real)    Second derivatives of penalty function         *
!                       w.r.t. internal parameter values.              *
! Working array:                                                       *
!   FVEC(NF)  (real)    Function values.                               *
!----------------------------------------------------------------------*
      integer i,icycle,iflag,nf,nx,calls
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision eps,f1,f2,fvec(nf),g2(nx),grd(nx),vdot,x(nx),    &
     &xsave,xstep,zero,one,two,half,epsmch
      parameter(zero=0d0,one=1d0,two=2d0,half=0.5d0,epsmch=1d-16)
      external fcn

      eps = sqrt(epsmch)

      do i = 1, nx
        xsave = x(i)
        xstep = eps * abs(xsave)
        if (xstep .eq. zero) xstep = eps
        do icycle = 1, 10
          x(i) = xsave + xstep
          call fcn(nf, nx, x, fvec, iflag)
          calls = calls + 1
          if (iflag .eq. 0) then
            f2 = vdot(nf, fvec, fvec)
            x(i) = xsave - xstep
            call fcn(nf, nx, x, fvec, iflag)
            calls = calls + 1
            if (iflag .eq. 0) then
              f1 = vdot(nf, fvec, fvec)
              go to 60
            endif
          endif
          xstep = half * xstep
        enddo
        f2 = fmin
        f1 = fmin
   60   continue
        grd(i) = (f2 - f1) / (two * xstep)
        g2(i) = (f2 - two * fmin + f1) / xstep**2
        if (g2(i) .eq. zero) g2(i) = one
        x(i) = xsave
      enddo

      call mtputi(x)

      end
      subroutine mtpsdf(covar,nx,work_4,work_5,work_6)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Force covariance matrix to be positive definite.                   *
! Updated:                                                             *
!   COVAR(*,*)        Covariance matrix.                               *
!----------------------------------------------------------------------*
      integer i,j,nval,nx
      double precision add,pmax,pmin,covar(nx,nx),work_4(nx,nx),        &
     &work_5(nx),work_6(nx),one,eps,epsmch
      parameter(one=1d0,eps=1d-3,epsmch=1d-16)

!---- Copy matrix and find eigenvalues.
      do i = 1, nx
        do j = 1, nx
          work_4(j,i) = covar(j,i)
        enddo
      enddo
      call symeig(work_4,nx,nx,work_5,nval,work_6)

!---- Enforce positive definiteness.
      pmin = work_5(1)
      pmax = work_5(1)
      do i = 1, nx
        if (work_5(i) .lt. pmin) pmin = work_5(i)
        if (work_5(i) .gt. pmax) pmax = work_5(i)
      enddo
      pmax = max(abs(pmax), one)
      if (pmin .le. epsmch * pmax) then
        add = eps * pmax - pmin
        do i = 1, nx
          covar(i,i) = covar(i,i) + add
        enddo
        print *, 'not posdef'
      endif

      end
      subroutine mtline(fcn,nf,nx,calls,x,dx,fvec,xsave,iflag)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Search for minimum along predicted direction.                      *
! Input:                                                               *
!   FCN       (subr)    Returns value of penalty function.             *
!   NF        (integer) Number of functions.                           *
!   NX        (integer) Number of parameters.                          *
!   X(NX)     (real)    Parameter values. On output, best estimate.    *
!   DX(NX)    (real)    Initial direction.                             *
!   FVEC(NF)  (real)    Function values.                               *
!   IFLAG     (integer) Error flag.                                    *
! Working array:                                                       *
!   XSAVE(NX) (real)    Save area for initial point.                   *
!----------------------------------------------------------------------*
      integer i,iflag,ipt,nf,npts,nvmax,nx,calls,maxpt
      parameter(maxpt=12)
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision c1,c2,den,dx(nx),f3,fval(3),fvec(nf),fvmin,      &
     &overal,ratio,s13,s21,s32,slam,slamax,slamin,sval(3),svmin,tol9,   &
     &undral,vdot,x(nx),xsave(nx),zero,one,two,half,tol8,alpha,slambg,  &
     &p100,p1000,epsmch
      parameter(zero=0d0,one=1d0,two=2d0,half=0.5d0,tol8=5d-2,alpha=2d0,&
     &slambg=5d0,p100=1d2,p1000=1d3,epsmch=1d-16)
      external fcn

!---- Initialize.
      overal = p1000
      undral = - p100
      sval(1) = zero
      fval(1) = fmin
      svmin = zero
      fvmin = fmin
      npts = 0

      slamin = zero
      do i = 1, nx
        xsave(i) = x(i)
        if (dx(i) .ne. zero) then
          ratio = abs(x(i) / dx(i))
          if (slamin .eq. zero  .or.  ratio .lt. slamin) slamin = ratio
        endif
      enddo
      if (slamin .eq. zero) slamin = epsmch
      slamin = slamin * epsmch
      slamax = slambg

!---- Compute function for move by DX.
      slam = one
   20 continue
      sval(2) = slam
      do i = 1, nx
        x(i) = xsave(i) + slam * dx(i)
      enddo
      call fcn(nf, nx, x, fvec, iflag)
      calls = calls + 1
      npts = npts + 1

!---- If machine becomes unstable, cut step.
      if (iflag .ne. 0) then
        slam = half * slam
        if (slam .gt. slamin) go to 20
        go to 400
      endif
      fval(2) = vdot(nf, fvec, fvec)
      if (fval(2) .lt. fvmin) then
        svmin = sval(2)
        fvmin = fval(2)
      endif
      if (slam .lt. one) go to 400

!---- Compute function for move by 1/2 DX.
      slam = half * slam
      sval(3) = slam
      do i = 1, nx
        x(i) = xsave(i) + slam * dx(i)
      enddo
      call fcn(nf, nx, x, fvec, iflag)
      calls = calls + 1
      npts = npts + 1
      if (iflag .ne. 0) go to 400
      fval(3) = vdot(nf, fvec, fvec)
      if (fval(3) .lt. fvmin) then
        svmin = sval(3)
        fvmin = fval(3)
      endif

!---- Begin iteration.
  200 continue
      slamax = max(slamax, alpha * abs(svmin))

!---- Quadratic interpolation using three points.
      s21 = sval(2) - sval(1)
      s32 = sval(3) - sval(2)
      s13 = sval(1) - sval(3)
      den = s21 * s32 * s13
      c2 = (s32 * fval(1) + s13 * fval(2) + s21 * fval(3)) / den
      c1 = ((sval(3) + sval(2)) * s32 * fval(1) +                       &
     &(sval(1) + sval(3)) * s13 * fval(2) +                             &
     &(sval(2) + sval(1)) * s21 * fval(3)) / den
      if (c2 .ge. zero) then
        slam = svmin + sign(slamax, c1 - two * c2 * svmin)
      else
        slam = c1 / (two * c2)
        if (slam .gt. svmin + slamax) slam = svmin + slamax
        if (slam .le. svmin - slamax) slam = svmin - slamax
      endif
      if (slam .gt. zero) then
        if (slam .gt. overal) slam = overal
      else
        if (slam .lt. undral) slam = undral
      endif

!---- See if new point coincides with a previous one.
  300 continue
      tol9 = tol8 * max(one, slam)
      do ipt = 1, 3
        if (abs(slam - sval(ipt)) .lt. tol9) go to 400
      enddo

!---- Compute function for interpolated point.
      do i = 1, nx
        x(i) = xsave(i) + slam * dx(i)
      enddo
      call fcn(nf, nx, x, fvec, iflag)
      calls = calls + 1
      npts = npts + 1
      if (iflag .ne. 0) go to 400
      f3 = vdot(nf, fvec, fvec)

!---- Find worst point of previous three.
      nvmax = 1
      if (fval(2) .gt. fval(nvmax)) nvmax = 2
      if (fval(3) .gt. fval(nvmax)) nvmax = 3

!---- If no improvement, cut interval.
      if (f3 .ge. fval(nvmax)) then
        if (npts .ge. maxpt) go to 400
        if (slam .gt. svmin) overal = min(overal, slam - tol8)
        if (slam .le. svmin) undral = max(undral, slam + tol8)
        slam = half * (slam + svmin)
        go to 300
      endif

!---- Accept new point; replace previous worst point.
      sval(nvmax) = slam
      fval(nvmax) = f3
      if (f3 .lt. fvmin) then
        svmin = slam
        fvmin = f3
      else
        if (slam .gt. svmin) overal = min(overal, slam - tol8)
        if (slam .lt. svmin) undral = max(undral, slam + tol8)
      endif
      if (npts .lt. maxpt) go to 200

!---- Common exit point: Return best point and step used.
  400 continue
      fmin = fvmin
      do i = 1, nx
        dx(i) = svmin * dx(i)
        x(i) = xsave(i) + dx(i)
      enddo
      call mtputi(x)

!---- Return Failure indication.
      iflag = 0
      if (svmin .eq. zero) iflag = 2

      end
      subroutine mtsimp(ncon,nvar,tol,calls,call_lim,vect,dvect,        &
     &fun_vect,w_iwa1,w_iwa2,w_iwa3)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Control routine for simplex minimization; SIMPLEX command.         *
! Attributes:                                                          *
!----------------------------------------------------------------------*
      integer calls,call_lim,ncon,nvar
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision tol,vect(*),dvect(*),fun_vect(*),w_iwa1(*),      &
     &w_iwa2(*),w_iwa3(*)
      external mtfcn

      icovar = 0
      ilevel = 0

!---- Too many variable parameters?
      if (nvar .gt. ncon) &
           call fort_warn('MTSIMP', &
           'More variables than constraints seen. SIMPLEX may not converge to optimal solution.')

!---- Call minimization routine.
      call mtgeti(vect, dvect)
      call mtsim1(mtfcn, ncon, nvar, calls, call_lim, tol,              &
     &vect, dvect, fun_vect, w_iwa1, w_iwa2, w_iwa3)

 9999 end

      subroutine mtsim1(fcn,nf,nx,calls,call_lim,tol,x,dx,fvec,psim,    &
     &fsim,wa)

      use matchfi
      implicit none


!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Minimization using the SIMPLEX method by Nelder and Mead.          *
!   (Computer Journal 7, 308 (1965).                                   *
! Input:                                                               *
!   FCN       (subr)    Returns value of function to be minimized.     *
!   NF        (integer) Number of functions.                           *
!   NX        (integer) Number of parameters.                          *
!   X(NX)     (real)    Parameter values. On output, best estimate.    *
!   DX(NX)    (real)    Parameter errors. On output, error estimate.   *
! Output:                                                              *
!   FVEC(NF)  (real)    Vector of function values in best point.       *
! Working arrays:                                                      *
!   PSIM(NX,0:NX)       Coordinates of simplex vertices.               *
!   FSIM(0:NX)          Function values in simplex vertices.           *
!   WA(NX,4)            Working vectors.                               *
!----------------------------------------------------------------------*
      integer i,idir,iflag,j,jh,jhold,jl,k,level,ncycl,nf,nrstrt,ns,nx, &
     &calls,call_lim,mbar,mstar,mstst,mrho, izero
      parameter(mbar=1,mstar=2,mstst=3,mrho=4)
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision f,f1,f2,fbar,fbest,frho,fstar,fstst,pb,pbest,    &
     &pmax,pmin,rho,step,vdot,tol,x(nx),dx(nx),fvec(nf),psim(nx,0:nx),  &
     &fsim(0:nx),wa(nx,4),alpha,beta,gamma,rhomin,rhomax,rho1,rho2,zero,&
     &two,three,half,p01,eps,epsmch
      parameter(alpha=1d0,beta=0.5d0,gamma=2d0,rhomin=4d0,rhomax=8d0,   &
     &rho1=1d0+alpha,rho2=rho1+alpha*gamma,zero=0d0,two=2d0,three=3d0,  &
     &half=5d-1,p01=1d-1,eps=1d-8,epsmch=1d-16)
      external fcn

      izero = 0
!---- Initialize penalty function.
      call mtcond(izero, nf, fvec, iflag)
!      call fcn(nf, nx, x, fvec, iflag)
      calls = calls + 1
      if (iflag .ne. 0) then
        call fort_warn('MTSIMP', ' stopped, possibly unstable')
        go to 400
      endif
      fmin = vdot(nf, fvec, fvec)
      edm = fmin
      nrstrt = 0

!---- Choose the initial simplex using single-parameter searches.
!     Keep initial point in PBAR.
  100 continue
      if (ilevel .ge. 1) call mtprnt('init', nx, x)
      fbar = fmin
      do i = 1, nx
        wa(i,mbar) = x(i)
        pbest = x(i)
        fbest = fmin
        step  = dx(i)

!---- Find proper initial direction and step.
        do idir = 1, 12
          x(i) = pbest + step
          call fcn(nf, nx, x, fvec, iflag)
          calls = calls + 1
          if (iflag .eq. 0) then
            f = vdot(nf, fvec, fvec)
            if (f .le. fbest) go to 120
          endif
          if (mod(idir,2) .eq. 0) step = p01 * step
          step = - step
        enddo
        go to 160

!---- Improvement found; attempt increasing steps.
  120   continue
        do ns = 1, 3
          pbest = x(i)
          fbest = f
          step = step * three
          x(i) = x(i) + step
          call fcn(nf, nx, x, fvec, iflag)
          calls = calls + 1
          if (iflag .ne. 0) go to 140
          f = vdot(nf, fvec, fvec)
          if (f .gt. fbest) go to 140
        enddo
        go to 160

!---- Backtrack to best point.
  140   continue
        x(i) = pbest
        f = fbest

!---- Store local minimum in i'th direction.
  160   continue
        fsim(i) = f
        do k = 1, nx
          psim(k,i) = x(k)
        enddo
      enddo

!---- Store initial point as 0'th vertex.
      jh = 0
      call mtrazz(nx, fbar, wa(1,mbar), fsim, psim, jh, jl)

!---- Extract best point.
      do i = 1, nx
        x(i) = psim(i,jl)
      enddo

!---- Print-out after setting up simplex.
      if (ilevel .ge. 2) then
        call mtprnt('progress',nx, x)
      endif
      ncycl = 0

!==== Start main loop.
  200 continue
      if (edm .lt. tol) then
        print *, 'converged'
      else if (calls .gt. call_lim) then
        print *, 'call limit'
      else

!---- Calculate PBAR and P*.
        do i = 1, nx
          pb = psim(i,0)
          do j = 1, nx
            pb = pb + psim(i,j)
          enddo
          wa(i,mbar) = (pb - psim(i,jh)) / float(nx)
          wa(i,mstar) = wa(i,mbar) + alpha * (wa(i,mbar) - psim(i,jh))
        enddo
        call fcn(nf, nx, wa(1,mstar), fvec, iflag)
        calls = calls + 1
        if (iflag .ne. 0) then
          fstar = two * fsim(jh)
        else
          fstar = vdot(nf, fvec, fvec)
        endif

!---- Point P* is better than point PSIM(*,JL).
        jhold = jh
        if (fstar .lt. fsim(jl)) then

!---- Try expanded point P**.
          do i = 1, nx
            wa(i,mstst) = wa(i,mbar) + gamma * (wa(i,mstar)-wa(i,mbar))
          enddo
          call fcn(nf, nx, wa(1,mstst), fvec, iflag)
          calls = calls + 1
          if (iflag .ne. 0) then
            fstst = two * fsim(jh)
            rho = zero
          else
            fstst = vdot(nf, fvec, fvec)

!---- Fit a parabola through FSIM(JH), F*, F**; minimum = RHO.
            f1 = (fstar - fsim(jh)) * rho2
            f2 = (fstst - fsim(jh)) * rho1
            rho = half * (rho2 * f1 - rho1 * f2) / (f1 - f2)
          endif

!---- Minimum inexistent ot too close to PBAR;
!     Use P** if it gives improvement; otherwise use P*.
          if (rho .lt. rhomin) then
            if (fstst .lt. fsim(jl)) then
              call mtrazz(nx, fstst, wa(1,mstst), fsim, psim,           &
     &jh, jl)
            else
              call mtrazz(nx, fstar, wa(1,mstar), fsim, psim,           &
     &jh, jl)
            endif

!---- Usable minimum found.
          else
            if (rho .gt. rhomax) rho = rhomax
            do i = 1, nx
              wa(i,mrho) = psim(i,jh) + rho * (wa(i,mbar) - psim(i,jh))
            enddo
            call fcn(nf, nx, wa(1,mrho), fvec, iflag)
            calls = calls + 1
            if (iflag .ne. 0) then
              frho = two * fsim(jh)
            else
              frho = vdot(nf, fvec, fvec)
            endif

!---- Select farthest point which gives decent improvement.
            if (frho .lt. fsim(jl) .and. frho .lt. fstst) then
              call mtrazz(nx, frho, wa(1,mrho), fsim, psim, jh, jl)
            else if (fstst .lt. fsim(jl)) then
              call mtrazz(nx, fstst, wa(1,mstst), fsim, psim,           &
     &jh, jl)
            else
              call mtrazz(nx, fstar, wa(1,mstar), fsim, psim,           &
     &jh, jl)
            endif
          endif

!---- F* is higher than FSIM(JL).
        else
          if (fstar .lt. fsim(jh)) then
            call mtrazz(nx, fstar, wa(1,mstar), fsim, psim, jh, jl)
          endif

!---- If F* is still highest value, try contraction,
!     giving point P** = PWRK(*,3).
          if (jhold .eq. jh) then
            do i = 1, nx
              wa(i,mstst) = wa(i,mbar)                                  &
     &+ beta * (psim(i,jh) - wa(i,mbar))
            enddo
            call fcn(nf, nx, wa(1,mstst), fvec, iflag)
            calls = calls + 1
            if (iflag .ne. 0) then
              fstst = two * fsim(jh)
            else
              fstst = vdot(nf, fvec, fvec)
            endif

!---- Restart algorithm, if F** is higher; otherwise use it.
            if (fstst .gt. fsim(jh)) then
              print *, 'failed'
              if (nrstrt .ne. 0) go to 300
              nrstrt = 1
              do j = 1, nx
                x(j) = psim(j,jl)
              enddo
              go to 100
            endif
            call mtrazz(nx, fstst, wa(1,mstst), fsim, psim, jh, jl)
          endif
        endif

!---- New minimum found.
        if (jl .eq. jhold) then
          nrstrt = 0
          ncycl = ncycl + 1
          level = 3
          if (mod(ncycl,10) .eq. 0) level = 2
          if (ilevel .ge. level) then
            call mtprnt('progress',nx, x)
          endif
!          if(calls.gt.(calls_old+dcalls)) then
          write(*,830) calls,fmin
!             calls_old=calls
!          endif
        endif
        go to 200
      endif

!==== End main loop: Try central point of simplex.
  300 continue
      do i = 1, nx
        pmin = psim(i,0)
        pmax = psim(i,0)
        pb = psim(i,0)
        do j = 1, nx
          pb = pb + psim(i,j)
        enddo
        wa(i,mbar) = (pb - psim(i,jh)) / float(nx)
        dx(i) = pmax - pmin
      enddo
      call fcn(nf, nx, wa(1,mbar), fvec, iflag)
      calls = calls + 1
      if (iflag .eq. 0) then
        fbar = vdot(nf, fvec, fvec)
        if (fbar .lt. fsim(jl)) then
          call mtrazz(nx, fbar, wa(1,mbar), fsim, psim, jh, jl)
        endif
      endif

!---- Recompute step sizes and extract best point.
      do i = 1, nx
        pmin = psim(i,0)
        pmax = psim(i,0)
        do j = 1, nx
          pmin = min(psim(i,j),pmin)
          pmax = max(psim(i,j),pmax)
        enddo
        dx(i) = pmax - pmin
        x(i) = psim(i,jl)
      enddo
      fmin = fsim(jl)
      call mtputi(x)

!---- Check for necessity to restart after final change.
      if (calls + 3 * nx .lt. call_lim  .and.  nrstrt .eq. 0  .and.     &
     &fmin .gt. two * (tol + epsmch)) then
        nrstrt = 1
        go to 100
      endif

!---- Final print-out.
  400 continue
      if (ilevel .ge. 1) then
        call mtprnt('final',nx, x)
      endif

      write(*,830) calls,fmin

  830 format('call:',I8,3x,'Penalty function = ',e16.8)
      end
      subroutine mtrazz(nvrr,fnew,pnew,fsim,psim,jh,jl)

      use matchfi
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Replace vertex in simplex whose function is highest.               *
!   Return indices of new best and worst points.                       *
! Input:                                                               *
!   NVAR      (integer) Number of parameters.                          *
!   FNEW      (real)    Function value for new vertex.                 *
!   PNEW(*)   (real)    Coordinates of new vertex.                     *
! Updated:                                                             *
!   FSIM(*)   (real)    Function values in (NVAR + 1) vertices.        *
!   PSIM(*,*) (real)    Coordinates of (NVAR + 1) vertices.            *
!   JH        (integer) Index of highest function value.               *
!   JL        (integer) Index of lowest function value.                *
!----------------------------------------------------------------------*
      integer i,jh,jl,nvrr
! icovar: functionality still unclear  HG 28.2.02
! ilevel: print level
      double precision fnew,fsim(0:nvrr),pnew(nvrr),psim(nvrr,0:nvrr),  &
     &ten
      parameter(ten=1d1)

!---- Replace vertex with highest function value.
      do i = 1, nvrr
        psim(i,jh) = pnew(i)
      enddo
      fsim(jh) = fnew

!---- Find indices of lowest and highest function value.
      jl = 0
      jh = 0
      do i = 1, nvrr
        if (fsim(i) .lt. fsim(jl)) jl = i
        if (fsim(i) .gt. fsim(jh)) jh = i
      enddo

!---- Get best value and estimated distance to minimum.
      fmin = fsim(jl)
      edm = min(ten * (fsim(jh) - fmin), fmin)

      end
