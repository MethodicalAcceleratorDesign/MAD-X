      subroutine fort_info(t1, t2)
      character *(*) t1, t2
      integer get_option
      if (get_option('info ') .ne. 0 .and. get_option('warn ') .ne. 0)  &
     &print '(a,1x,a,1x,a)', '++++++ info:', t1, t2
      end
      subroutine fort_warn(t1, t2)
      character *(*) t1, t2
      integer get_option
      if (get_option('warn ') .ne. 0)                                   &
     &print '(a,1x,a,1x,a)', '++++++ warning:', t1, t2
      end
      subroutine getclor(orbit0, rt, tt, error)
!----------------------------------------------------------------------*
! Purpose:
!   Get periodic closed orbit (e.g. at start of Twiss), 
!   first + second order one-turn map
! Input:
!   orbit0(6)   (real)  initial guess
! Output:
!   rt(6,6)     (real)  one-turn matrix
!   tt(6,6,6)   (real)  one-turn second-order map
!   error       (int)   error flag (0: OK, else != 0)
!----------------------------------------------------------------------*
      implicit none
      include 'twiss0.fi'
      double precision orbit0(6), rt(6,6), tt(6,6,6)
      double precision opt(fundim)
      integer error
      call m66one(rt)
      call dzero(opt,fundim)
      call tmclor(orbit0, .true., .true., opt, rt, tt, error)
      end
      subroutine m66add(term1,term2,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Add two matrices.
! Input:
!   TERM1(6,6)  (real)  First term.
!   TERM2(6,6)  (real)  Second term.
! Output:
!   TARGET(6,6) (real)  Sum: TARGET = TERM1 + TERM2.
!----------------------------------------------------------------------*
      integer i,j
      double precision target(6,6),term1(6,6),term2(6,6)

      do i = 1, 6
        do j = 1, 6
          target(i,j) = term1(i,j) + term2(i,j)
        enddo
      enddo

      end
      subroutine m66byv(amat,avec,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Multiply matrix times vector.
! Input:
!   AMAT(6,6)   (real)  Input matrix.
!   AVEC(6)     (real)  Input vector.
! Output:
!   TARGET(6)   (real)  Output vector: TARGET = AMAT * AVEC.
!----------------------------------------------------------------------*
      integer i,j
      double precision amat(6,6),avec(6),target(6),temp(6)

      call dzero(temp,6)
      do i = 1, 6
        do j = 1, 6
          temp(i) = temp(i) + amat(i,j) * avec(j)
        enddo
      enddo

      call dcopy(temp,target,6)

      end
      subroutine m66cpy(source,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Copy matrix.
! Input:
!   SOURCE(6,6) (real)  Input matrix.
! Output:
!   TARGET(6,6) (real)  Output matrix: TARGET = SOURCE.
!----------------------------------------------------------------------*
      integer i,j
      double precision source(6,6),target(6,6)

      do i = 1, 6
        do j = 1, 6
          target(i,j) = source(i,j)
        enddo
      enddo

      end
      subroutine m66div(anum,aden,target,eflag)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   "Divide" matrices, i. e. postmultiply with inverse of denominator.
! Input:
!   ANUM(6,6)   (real)  "Numerator" matrix.
!   ADEN(6,6)   (real)  "Denominator" matrix.
! Output:
!   TARGET(6,6) (real)  "Quotient" matrix: TARGET = ANUM * ADEN**(-1).
!   EFLAG       (logical) Error flag.
!----------------------------------------------------------------------*
      logical eflag
      integer i,irank,j
      double precision aden(6,6),anum(6,6),augmat(6,12),target(6,6)

!---- Copy input to local array.
      do i = 1, 6
        do j = 1, 6
          augmat(i,j)   = aden(i,j)
          augmat(i,j+6) = anum(i,j)
        enddo
      enddo

!---- Solve resulting system.
      call solver(augmat,6,6,irank)
      if (irank .lt. 6) then
        eflag = .true.

!---- Copy result.
      else
        eflag = .false.
        do i = 1, 6
          do j = 1, 6
            target(i,j) = augmat(i,j+6)
          enddo
        enddo
      endif

      end
      subroutine m66exp(source,target,eflag)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   "Exponentiate" matrix.
!   Original author:    Liam Healy.
! Input:
!   SOURCE(6,6) (real)  Input matrix.
! Output:
!   TARGET(6,6) (real)  Output matrix: TARGET = exp(SOURCE).
!   EFLAG     (logical) Error flag.
!----------------------------------------------------------------------*
      logical eflag
      integer i,j
      double precision b(6,6),c(6,6),source(6,6),target(6,6),one,two,   &
     &twelve
      parameter(one=1d0,two=2d0,twelve=12d0)

      call m66mpy(source,source,b)
      call m66mpy(source,b,c)
      do j = 1, 6
        do i = 1, 6
          b(i,j) = (source(i,j) - c(i,j) / twelve) / two
          c(i,j) = - b(i,j)
        enddo
        b(j,j) = b(j,j) + one
        c(j,j) = c(j,j) + one
      enddo
      call m66div(b,c,target,eflag)

      end
      subroutine m66inv(source,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Invert symplectic matrix.
! Input:
!   SOURCE(6,6) (real)  Input matrix.
! Output:
!   TARGET(6,6) (real)  Output matrix: TARGET = tr(J) * tr(SOURCE) * J.
!----------------------------------------------------------------------*
      integer i
      double precision source(6,6),target(6,6),temp(6,6)

!---- TEMP = transpose(SOURCE) * J.
      do i = 1, 6
        temp(i,1) = - source(2,i)
        temp(i,2) = + source(1,i)
        temp(i,3) = - source(4,i)
        temp(i,4) = + source(3,i)
        temp(i,5) = - source(6,i)
        temp(i,6) = + source(5,i)
      enddo

!---- TARGET = transpose(J) * TEMP.
      do i = 1, 6
        target(1,i) = - temp(2,i)
        target(2,i) = + temp(1,i)
        target(3,i) = - temp(4,i)
        target(4,i) = + temp(3,i)
        target(5,i) = - temp(6,i)
        target(6,i) = + temp(5,i)
      enddo

      end
      subroutine m66mak(f2,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Compute matrix TARGET corresponding to Lie polynomial F2.
!   Original author:    Liam Healy.
! Input:
!   F2          (poly)  Polynomial of order 2.
! Output:
!   TARGET(6,6) (real)  Output matrix: TARGET * v = - [J,v].
!----------------------------------------------------------------------*
      double precision f2(*),target(6,6),two
      parameter(two=2d0)

      target(1,1) = - f2(8)
      target(1,2) = - two * f2(13)
      target(1,3) = - f2(14)
      target(1,4) = - f2(15)
      target(1,5) = - f2(16)
      target(1,6) = - f2(17)
      target(2,1) = two * f2(7)
      target(2,2) = f2(8)
      target(2,3) = f2(9)
      target(2,4) = f2(10)
      target(2,5) = f2(11)
      target(2,6) = f2(12)
      target(3,1) = - f2(10)
      target(3,2) = - f2(15)
      target(3,3) = - f2(19)
      target(3,4) = - two * f2(22)
      target(3,5) = - f2(23)
      target(3,6) = - f2(24)
      target(4,1) = f2(9)
      target(4,2) = f2(14)
      target(4,3) = two * f2(18)
      target(4,4) = f2(19)
      target(4,5) = f2(20)
      target(4,6) = f2(21)
      target(5,1) = - f2(12)
      target(5,2) = - f2(17)
      target(5,3) = - f2(21)
      target(5,4) = - f2(24)
      target(5,5) = - f2(26)
      target(5,6) = - two * f2(27)
      target(6,1) = f2(11)
      target(6,2) = f2(16)
      target(6,3) = f2(20)
      target(6,4) = f2(23)
      target(6,5) = two * f2(25)
      target(6,6) = f2(26)

      end
      subroutine m66mpy(fact1,fact2,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Multiply two matrices.
!   TARGET may coincide with one of the factors.
! Input:
!   FACT1(6,6)  (real)  First factor.
!   FACT2(6,6)  (real)  Second factor.
! Output:
!   TARGET(6,6) (real)  Product matrix: TARGET = FACT1 * FACT2.
!----------------------------------------------------------------------*
      integer i,j,k
      double precision fact1(6,6),fact2(6,6),target(6,6),temp(6,6)

      call dzero(temp,36)
      do k = 1, 6
        do j = 1, 6
          do i = 1, 6
            temp(i,k) = temp(i,k) + fact1(i,j) * fact2(j,k)
          enddo
        enddo
      enddo
      call dcopy(temp,target,36)

      end
      subroutine m66mtr(fact1,fact2,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Multiply a matrix with the transpose of another matrix.
!   TARGET must not coincide with either factor.
! Input:
!   FACT1(6,6)  (real)  First factor.
!   FACT2(6,6)  (real)  Second factor (will be transposed).
! Output:
!   TARGET(6,6) (real)  Product: TARGET = FACT1 * tr(FACT2).
!----------------------------------------------------------------------*
      integer i,j,k
      double precision fact1(6,6),fact2(6,6),target(6,6)

      call dzero(target,36)
      do j = 1, 6
        do k = 1, 6
          do i = 1, 6
            target(i,j) = target(i,j) + fact1(i,k) * fact2(j,k)
          enddo
        enddo
      enddo

      end
      subroutine m66nrm(fm,res)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Computes the norm of a matrix.
!   Reference:          L. Collatz,
!                       Functional Analysis & Numerical Mathematics.
!   Source:             MARYLIE, Version 3.0.
! Input:
!   FM(6,6)     (real)  Input matrix.
! Output:
!   RES         (real)  Norm of FM: RES = max abs column sum.
!----------------------------------------------------------------------*
      integer i,j
      double precision fm(6,6),res,sum,zero
      parameter(zero=0d0)

      res = zero
      do j = 1, 6
        sum = zero
        do i = 1, 6
          sum = sum + abs(fm(i,j))
        enddo
        res = max(res,sum)
      enddo

      end
      subroutine m66one(target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Set matrix to unity.
! Output:
!   TARGET(6,6) (real)  Unit matrix: TARGET = I.
!----------------------------------------------------------------------*
      integer j
      double precision target(6,6),one
      parameter(one=1d0)

      call dzero(target,36)
      do j = 1, 6
        target(j,j) = one
      enddo

      end
      subroutine m66ref(source,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Reflect symplectic first order transform.
! Input:
!   SOURCE(6,6) (real)  Input matrix.
! Output:
!   TARGET(6,6) (real)  Reflected matrix.
!----------------------------------------------------------------------*
      integer i
      double precision source(6,6),target(6,6),temp(6,6)

!---- TEMP = transpose(SOURCE) * J * signs.
      do i = 1, 6
        temp(i,1) =   source(2,i)
        temp(i,2) =   source(1,i)
        temp(i,3) =   source(4,i)
        temp(i,4) =   source(3,i)
        temp(i,5) = - source(6,i)
        temp(i,6) = - source(5,i)
      enddo

!---- TARGET = signs * transpose(J) * TEMP.
      do i = 1, 6
        target(1,i) =   temp(2,i)
        target(2,i) =   temp(1,i)
        target(3,i) =   temp(4,i)
        target(4,i) =   temp(3,i)
        target(5,i) = - temp(6,i)
        target(6,i) = - temp(5,i)
      enddo

      end
      subroutine m66scl(scalar,source,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Multiply matrix by scalar.
! Input:
!   SCALAR      (real)  Scale factor.
!   SOURCE(6,6) (real)  Input matrix.
! Output:
!   TARGET(6,6) (real)  Scaled matrix: TARGET = SCALAR * SOURCE.
!----------------------------------------------------------------------*
      integer i,j
      double precision scalar,source(6,6),target(6,6)

      do i = 1, 6
        do j = 1, 6
          target(i,j) = scalar * source(i,j)
        enddo
      enddo

      end
      logical function m66sta(amat)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Check effect of a matrix on momentum.
! Input:
!   AMAT(6,6)   (real)  Input matrix.
! Result:
!   .TRUE.              For static case     (constant p).
!   .FALSE.             For dynamic case    (variable p).
!----------------------------------------------------------------------*
      integer j
      double precision amat(6,6),tol,one
      parameter(one=1d0,tol=1d-12)

      m66sta = abs(amat(6,6) - one) .le. tol
      do j = 1, 5
        m66sta = m66sta .and. abs(amat(6,j)) .le. tol
      enddo

      end
      subroutine m66sub(term1,term2,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Subtract two matrices.
! Input:
!   TERM1(6,6)  (real)  Minuend matrix.
!   TERM2(6,6)  (real)  Subtrahend matrix.
! Output:
!   TARGET(6,6) (real)  Difference matrix: TARGET = TERM1 - TERM2.
!----------------------------------------------------------------------*
      integer i,j
      double precision target(6,6),term1(6,6),term2(6,6)

      do j = 1, 6
        do i = 1, 6
          target(i,j) = term1(i,j) - term2(i,j)
        enddo
      enddo

      end
      subroutine m66trm(fact1,fact2,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Multiply the transpose of a matrix with another matrix.
!   TARGET must not coincide with either factor.
! Input:
!   FACT1(6,6)  (real)  First factor (will be transposed).
!   FACT2(6,6)  (real)  Second factor.
! Output:
!   TARGET(6,6) (real)  Product: TARGET = tr(FACT1) * FACT2.
!----------------------------------------------------------------------*
      integer i,j,k
      double precision fact1(6,6),fact2(6,6),target(6,6)

      call dzero(target,36)
      do j = 1, 6
        do k = 1, 6
          do i = 1, 6
            target(i,j) = target(i,j) + fact1(k,i) * fact2(k,j)
          enddo
        enddo
      enddo

      end
      subroutine m66tp(source,target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Transpose a matrix.
!   TARGET and SOURCE may overlap.
! Input:
!   SOURCE(6,6) (real)  Input matrix.
! Output:
!   TARGET(6,6) (real)  Transposed matrix: TARGET = tr(SOURCE).
!----------------------------------------------------------------------*
      integer i,j
      double precision source(6,6),target(6,6),temp(6,6)

      do i = 1, 6
        do j = 1, 6
          temp(j,i) = source(i,j)
        enddo
      enddo
      call m66cpy(temp,target)

      end
      subroutine m66zro(target)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Clear a matrix to zero.
! Output:
!   TARGET(6,6) (real)  Zero matrix: TARGET = 0.
!----------------------------------------------------------------------*
      double precision target(6,6)

      call dzero(target,36)

      end
      subroutine solver(augmat,ndim,mdim,irank)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Solve the linear equation  A * X = B.
! Input:
!   AUGMAT(n,n+m)       A(n,n), augmented by B(n,m).
!   NDIM, MDIM          n, m.
! Output:
!   AUGMAT(n,n+m)       Identity(n,n), augmented by X(n,m).
!   IRANK               Rank of A.
!----------------------------------------------------------------------*
      integer ic,ip,ir,irank,it,mdim,nc,ndim,nr
      double precision augmat(ndim,ndim+mdim),h,pivot,zero
      parameter(zero=0d0)

      nr = ndim
      nc = ndim + mdim
      irank = 0
      do it = 1, nr
        pivot = zero
        ip = 0
        do ir = it, nr
          if (abs(augmat(ir,it)) .ge. abs(pivot)) then
            pivot = augmat(ir,it)
            ip = ir
          endif
        enddo

        if (pivot .eq. zero) go to 9999
        irank = it

        do ic = 1, nc
          augmat(ip,ic) = augmat(ip,ic) / pivot
        enddo

        if (ip .ne. it) then
          do ic = 1, nc
            h = augmat(ip,ic)
            augmat(ip,ic) = augmat(it,ic)
            augmat(it,ic) = h
          enddo
        endif

        do ir = 1, nr
          if (ir .ne. it) then
            h = augmat(ir,it)
            do ic = 1, nc
              augmat(ir,ic) = augmat(ir,ic) - h * augmat(it,ic)
            enddo
          endif
        enddo
      enddo

      irank = ndim

 9999 end
      subroutine symsol(a,n,eflag,work_1,work_2,work_3)
      implicit none
!----------------------------------------------------------------------*
! Purpose:
!   Invert symmetric matrix.
! Input:
!   A(*,*)    (real)    Matrix to be inverted.
!   N         (integer) Actual size of A.
! Output:
!   A(*,*)    (real)    Inverted matrix.
!   EFLAG     (logical) Error flag.
!----------------------------------------------------------------------*
      logical eflag
      integer i,j,k,n
      double precision a(n,n),si,work_1(n),work_2(n),work_3(n),zero,one
      parameter(zero=0d0,one=1d0)

!---- Scale upper triangle.
      eflag = .true.
      do i = 1, n
        si = a(i,i)
        if (si .le. zero) go to 100
        work_1(i) = one / sqrt(si)
      enddo
      do i = 1, n
        do j = i, n
          a(i,j) = a(i,j) * work_1(i) * work_1(j)
        enddo
      enddo

!---- Invert upper triangle.
      do i = 1, n
        if (a(i,i) .eq. zero) go to 100
        work_2(i) = one
        work_3(i) = one / a(i,i)
        a(i,i) = zero
        do j = 1, n
          if (j .lt. i) then
            work_2(j) = a(j,i)
            work_3(j) = work_2(j) * work_3(i)
            a(j,i) = zero
          else if (j .gt. i) then
            work_2(j) = a(i,j)
            work_3(j) = - work_2(j) * work_3(i)
            a(i,j) = zero
          endif
        enddo
        do j = 1, n
          do k = j, n
            a(j,k) = a(j,k) + work_2(j) * work_3(k)
          enddo
        enddo
      enddo

!---- Rescale upper triangle and symmetrize.
      do i = 1, n
        do j = i, n
          a(i,j) = a(i,j) * work_1(i) * work_1(j)
          a(j,i) = a(i,j)
        enddo
      enddo
      eflag = .false.

  100 continue

      end
      subroutine symeig(a,nd,n,eigen,nval,work)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Eigenvalues of a real symmetric matrix in ascending order.         *
! Input:                                                               *
!   A(ND,ND)  (real)    Symmetric input matrix; destroyed by call.     *
!   N         (integer) Rank of matrix.                                *
! Output:                                                              *
!   EIGEN(*)  (real)    Eigenvalues of A in descending order.          *
!   NVAL      (integer) Number of eigenvalues found.                   *
!----------------------------------------------------------------------*
      integer i,it,j,k,l,m,n,nd,nval,itmax
      parameter(itmax=15)
      double precision b,c,f,g,h,p,r,s,work(nd),a(nd,nd),eigen(nd),zero,&
     &one,two,four,big,eps
      parameter(zero=0d0,one=1d0,two=2d0,four=4d0,big=1d10,eps=1d-20)

!---- Matrix is 1 * 1.
      nval = n
      if (n .le. 0) go to 300
      if (n .eq. 1) then
        eigen(1) = a(1,1)
        go to 300
      endif

!---- Matrix is 2 * 2.
      if (n .eq. 2) then
        f = a(1,1) + a(2,2)
        g = sqrt((a(1,1) - a(2,2))**2 + four * a(2,1)**2)
        eigen(1) = (f - g) / two
        eigen(2) = (f + g) / two
        go to 300
      endif

!---- N is at least 3, reduce to tridiagonal form.
      do i = n, 3, -1
        g = zero
        do k = 1, i-2
          g = g + a(i,k)**2
        enddo
        eigen(i) = a(i,i)
        if (g .eq. zero) then
          work(i) = a(i,i-1)
        else
          h = g + a(i,i-1)**2
          work(i) = sign(sqrt(h),a(i,i-1))
          h = h + a(i,i-1) * work(i)
          a(i,i-1) = a(i,i-1) + work(i)
          f = zero
          do j = 1, i-1
            g = zero
            do k = 1, i-1
              if (k .le. j) then
                g = g + a(j,k) * a(i,k)
              else
                g = g + a(k,j) * a(i,k)
              endif
            enddo
            work(j) = g / h
            f = f + work(j) * a(i,j)
          enddo
          do j = 1, i-1
            work(j) = work(j) - (f / (h + h)) * a(i,j)
            do k = 1, j
              a(j,k) = a(j,k) - a(i,j) * work(k) - work(j) * a(i,k)
            enddo
          enddo
        endif
      enddo
      work(2) = a(2,1)
      work(1) = zero
      eigen(2) = a(2,2)
      eigen(1) = a(1,1)

!---- Iterate on tridiagonal matrix.
      do i = 2, n
        work(i-1) = work(i)
      enddo

      work(n) = zero
      f = zero
      b = zero
      do l = 1, n
        b = max(eps*(abs(eigen(l))+abs(work(l))),b)
        do m = l, n
          if (abs(work(m)) .le. b) go to 130
        enddo
        m = n
  130   if (m .ne. l) then
          do it = 1, itmax
            p = (eigen(l+1) - eigen(l)) / (two * work(l))
            if (abs(p) .gt. big) then
              r = abs(p)
            else
              r = sqrt(p*p+one)
            endif
            h = eigen(l) - work(l) / (p + sign(r,p))
            do i = l, n
              eigen(i) = eigen(i) - h
            enddo
            f = f + h
            p = eigen(m)
            c = one
            s = zero
            do i = m-1, l, -1
              g = c * work(i)
              h = c * p
              r = sqrt(work(i)**2+p**2)
              work(i+1) = s * r
              s = work(i) / r
              c = p / r
              p = c * eigen(i) - s * g
              eigen(i+1) = h + s * (c * g + s * eigen(i))
            enddo
            work(l) = s * p
            eigen(l) = c * p
            if (abs(work(l)) .le. b) go to 170
          enddo
          nval = l - 1
          go to 300
        endif
  170   p = eigen(l) + f
        do i = l, 2, -1
          if (p .ge. eigen(i-1)) go to 190
          eigen(i) = eigen(i-1)
        enddo
        i = 1
  190   eigen(i) = p
      enddo
  300 continue

      end
      subroutine dcopy(in,out,n)
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Copy arrays.                                                       *
! Input:                                                               *
!   in  (double)    array to be copied.                                *
!   n   (integer)   array length.                                      *
! Output:                                                              *
!   out (double)    target array.                                      *
!----------------------------------------------------------------------*
      implicit none
      integer n, i
      double precision in(*), out(*)

      do i = 1, n
        out(i) = in(i)
      enddo

      end
      subroutine dzero(vector,n)
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Copy an array.                                                     *
! Input:                                                               *
!   n      (integer) array length.                                     *
! Input/output:                                                        *
!   vector (double)  array to be zeroed.                               *
!----------------------------------------------------------------------*
      implicit none
      integer n, i
      double precision vector(*),zero
      parameter(zero=0d0)

      do i = 1, n
        vector(i) = zero
      enddo

      end
      subroutine aawarn(rout,text)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Print warning message.                                             *
! Input:                                                               *
!   ROUT      (char)    Calling routine name.                          *
!   TEXT      (char)    Message.                                       *
!----------------------------------------------------------------------*
      character*(*) rout,text

      print *,  '++++++ warning: ',rout,text

      end
      subroutine aafail(rout,text)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Print fatal error message.                                         *
! Input:                                                               *
!   ROUT      (char)    Calling routine name.                          *
!   TEXT      (char)    Message.                                       *
!----------------------------------------------------------------------*
      character*(*) rout,text

      print *,  '+-+-+- fatal: ',rout,text
      stop

      end
      double precision function proxim(x,y)
!----------------------------------------------------------------------*
!   Proximity function of x and y.                                     *
!   If angle is larger than pi between vector x and y, 2pi is added to *
!   to this angle                                                      *
!----------------------------------------------------------------------*
      implicit none
      double precision x,y,twopi,get_variable

      twopi=get_variable('twopi ')
      proxim = x+twopi*anint((y-x)/twopi)

      end
      character * 48 function charconv(tint)
!----------------------------------------------------------------------*
! purpose:                                                             *
!   converts integer array to string (based on ascii)                  *
! input:                                                               *
!   tint  (int array)  1 = length, rest = string                       *
!----------------------------------------------------------------------*
      implicit none
      integer tint(*)
      integer i, j, m, n
      parameter (m = 128)
      character *(m) letter
      data letter /                                                     &
     &'                                !"#$%&''()*+,-./0123456789:;<=>?@&
     &ABCDEFGHIJKLMNOPQRSTUVWXYZ[ ]^_`abcdefghijklmnopqrstuvwxyz{|}~'/
      charconv = ' '
      n = tint(1)
      do i = 1, n
        j = tint(i+1)
        if (j .lt. m)  charconv(i:i) = letter(j:j)
      enddo
      end
      subroutine laseig(fm,reeig,aieig,am)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Return eigenvalues and eigenvectors of a 4x4 matrix.               *
! Input:                                                               *
!   FM(6,6)   (real)    Matrix to be transformed.                      *
! Output:                                                              *
!   REEIG(6)  (real)    Real parts of eigenvalues.                     *
!   AIEIG(6)  (real)    Imaginary parts of eigenvalues.                *
!   AM(6,6)   (real)    Transforming matrix, contains eigenvectors.    *
!----------------------------------------------------------------------*
      integer i,ihi,ilo,info,ipind,iqind,j,k,mdim,nn,kpnt(6)
      double precision fm(6,6),reeig(6),aieig(6),am(6,6),aival(6),big,c,&
     &d(6),dx,dy,pb,reval(6),s,tm(6,6),zero,one
      parameter(zero=0d0,one=1d0,ilo=1,ihi=4,mdim=6,nn=4)

!---- Compute eigenvalues and vectors.
      call m66cpy(fm,tm)
      call m66one(am)
      call orthes(mdim,nn,ilo,ihi,tm,d)
      call ortran(mdim,nn,ilo,ihi,tm,d,am)
      call hqr2(mdim,nn,ilo,ihi,tm,reval,aival,am,info)
      if (info .ne. 0) then
        write (6, 910) ((fm(i,k), k = 1, 6), i = 1, 6)
  910   format('Unable to find eigenvalues for matrix:'/                &
     &  (6f12.6))
        call aafail('LASEIG',' Unable to find eigenvalues for matrix')
        go to 999
      endif
!---- Normalize the eigenvectors.
      do k = 1, 5, 2
        pb = zero
        do ipind = 2, 6, 2
          iqind = ipind - 1
          pb = pb + am(iqind,k) * am(ipind,k+1)                         &
     &    - am(ipind,k) * am(iqind,k+1)
        enddo
        s = sqrt(abs(pb))
        if (pb .lt. zero) then
          aival(k) = - aival(k)
          aival(k+1) = - aival(k+1)
        endif
        do i = 1, 6
          am(i,k)   = am(i,k) / s
          am(i,k+1) = am(i,k+1) * (s / pb)
        enddo
      enddo
!---- Sort these eigenvectors.
      call m66cpy(am,tm)
!---- Find the eigenvectors with the largest vertical component.
      big = zero
      kpnt(3) = 1
      do i = 1, 3, 2
        c = tm(3,i)**2 + tm(3,i+1)**2 + tm(4,i)**2 + tm(4,i+1)**2
        if (c .gt. big) then
          big = c
          kpnt(3) = i
        endif
      enddo
!---- Find the remaining vector.
      do i = 1, 3, 2
        if (i .ne. kpnt(3)) kpnt(1) = i
      enddo
!---- Reorder vectors.
      do i = 1, 3, 2
        k = kpnt(i)
        reeig(i) = reval(k)
        aieig(i) = aival(k)
        reeig(i+1) = reval(k+1)
        aieig(i+1) = aival(k+1)
        do j = 1, 6
          am(j,i) = tm(j,k)
          am(j,i+1) = tm(j,k+1)
        enddo
      enddo
      reeig(5) = one
      aieig(5) = zero
      reeig(6) = one
      aieig(6) = zero
!---- Rephase the result.
      call m66one(tm)
      dx = sqrt(am(1,1)**2 + am(1,2)**2)
      tm(1,1) = am(1,1) / dx
      tm(2,1) = am(1,2) / dx
      tm(1,2) = - tm(2,1)
      tm(2,2) = tm(1,1)
      dy = sqrt(am(3,3)**2 + am(3,4)**2)
      tm(3,3) = am(3,3) / dy
      tm(4,3) = am(3,4) / dy
      tm(3,4) = - tm(4,3)
      tm(4,4) = tm(3,3)
      call m66mpy(am,tm,am)
  999 end
      subroutine ladeig(fm,reeig,aieig,am)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Return eigenvalues and eigenvectors of a 6x6 matrix.               *
! Input:                                                               *
!   FM(6,6)   (real)    Matrix to be transformed.                      *
! Output:                                                              *
!   REEIG(6)  (real)    Real parts of eigenvalues.                     *
!   AIEIG(6)  (real)    Imaginary parts of eigenvalues.                *
!   AM(6,6)   (real)    Transforming matrix, contains eigenvectors.    *
!----------------------------------------------------------------------*
      integer i,ihi,ilo,info,j,k,mdim,nn,kpnt(6)
      double precision fm(6,6),reeig(6),aieig(6),am(6,6),aival(6),big,c,&
     &d(6),dt,dx,dy,pb,reval(6),s,tm(6,6),zero
      parameter(zero=0d0,ilo=1,ihi=6,mdim=6,nn=6)

!---- Compute eigenvalues and eigenvectors.
      call m66cpy(fm,tm)
      call orthes(mdim,nn,ilo,ihi,tm,d)
      call ortran(mdim,nn,ilo,ihi,tm,d,am)
      call hqr2(mdim,nn,ilo,ihi,tm,reval,aival,am,info)
      if (info .ne. 0) then
        write (6, 910) ((fm(i,k), k = 1, 6), i = 1, 6)
  910   format('Unable to find eigenvalues for matrix:'/                &
     &  (6f12.6))
        call aafail('LADEIG',' Unable to find eigenvalues for matrix')
        go to 9999
      endif
!---- Normalize the eigenvectors.
      do k = 1, 5, 2
        pb = zero
        do i = 1, 5, 2
          pb = pb + am(i,k) * am(i+1,k+1)                               &
     &    - am(i+1,k) * am(i,k+1)
        enddo
        s = sqrt(abs(pb))
        if (pb .lt. zero) then
          aival(k) = - aival(k)
          aival(k+1) = - aival(k+1)
        endif
        do i = 1, 6
          am(i,k)   = am(i,k) / s
          am(i,k+1) = am(i,k+1) * (s / pb)
        enddo
      enddo
!---- Copy vectors to temporary array.
      call m66cpy(am,tm)
!---- Find the vector with the largest vertical component.
      big = zero
      kpnt(3) = 1
      do i = 1, 5, 2
        c = tm(3,i)**2 + tm(3,i+1)**2 + tm(4,i)**2 + tm(4,i+1)**2
        if (c .gt. big) then
          big = c
          kpnt(3) = i
        endif
      enddo
!---- Find  the vector with the largest horizontal component.
      kpnt(1) = 1
      big = zero
      do i = 1, 5, 2
        if (i .ne. kpnt(3)) then
          c = tm(1,i)**2 + tm(1,i+1)**2 + tm(2,i)**2 + tm(2,i+1)**2
          if (c .gt. big) then
            big = c
            kpnt(1) = i
          endif
        endif
      enddo
!---- Find the remaining vector.
      do i = 1, 5, 2
        if (i .ne. kpnt(3)  .and.  i .ne. kpnt(1)) kpnt(5) = i
      enddo
!---- Reorder vectors.
      do i = 1, 5, 2
        k = kpnt(i)
        reeig(i) = reval(k)
        aieig(i) = aival(k)
        reeig(i+1) = reval(k+1)
        aieig(i+1) = aival(k+1)
        do j = 1, 6
          am(j,i) = tm(j,k)
          am(j,i+1) = tm(j,k+1)
        enddo
      enddo
!---- Rephase the result.
      call m66one(tm)
      dx = sqrt(am(1,1)**2 + am(1,2)**2)
      tm(1,1) = am(1,1) / dx
      tm(2,1) = am(1,2) / dx
      tm(1,2) = - tm(2,1)
      tm(2,2) = tm(1,1)
      dy = sqrt(am(3,3)**2 + am(3,4)**2)
      tm(3,3) = am(3,3) / dy
      tm(4,3) = am(3,4) / dy
      tm(3,4) = - tm(4,3)
      tm(4,4) = tm(3,3)
      dt = sqrt(am(5,5)**2 + am(5,6)**2)
      tm(5,5) = am(5,5) / dt
      tm(6,5) = am(5,6) / dt
      tm(5,6) = - tm(6,5)
      tm(6,6) = tm(5,5)
      call m66mpy(am,tm,am)
 9999 end
      subroutine orthes(ndim,n,ilow,iupp,a,d)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Converts an unsymmetric real matrix, A, to upper Hessenberg form   *
!   applying successive orthogonal transformations.                    *
!                                                                      *
!   Translation of the ALGOL procedure ORTHES in:                      *
!   Handbook Series Linear Algebra,                                    *
!   Num. Math. 12, 349-368 (1968) by R. S. Martin and J. H. Wilkinson. *
! Input:                                                               *
!   N         (integer) Order of the matrix A.                         *
!   ILOW,IUPP (integer) Determine a submatrix, set by BALANC.          *
!                       May be set to 1 and N respectively.            *
!   A(NDIM,N) (real)    Input matrix.                                  *
! Output:                                                              *
!   A(NDIM,N) (real)    The matrix A, converted to upper Hessenberg.   *
!                       The lower triangle contains information        *
!                       about the orthogonal transformations.          *
!   D(N)      (real)    Further information.                           *
!----------------------------------------------------------------------*
      integer i,ilow,iupp,j,m,n,ndim
      double precision a(ndim,n),d(n),f,g,h,scale,zero
      parameter(zero=0d0)

      do m = ilow + 1, iupp - 1
        h = zero
        d(m) = zero
!---- Find scale factor.
        scale = zero
        do i = m, iupp
          scale = scale + abs(a(i,m-1))
        enddo
        if (scale .ne. zero) then
          do i = iupp, m, - 1
            d(i) = a(i,m-1) / scale
            h = h + d(i) * d(i)
          enddo
          g = sign(sqrt(h),d(m))
          h = h + d(m) * g
          d(m) = d(m) + g
!---- Form (I - (u*uT) / h) * A.
          do j = m, n
            f = zero
            do i = iupp, m, - 1
              f = f + d(i) * a(i,j)
            enddo
            f = f / h
            do i = m, iupp
              a(i,j) = a(i,j) - f * d(i)
            enddo
          enddo
!---- Form (I - (u*uT) / h) * A * (I - (u*uT) / h).
          do i = 1, iupp
            f = zero
            do j = iupp, m, - 1
              f = f + d(j) * a(i,j)
            enddo
            f = f / h
            do j = m, iupp
              a(i,j) = a(i,j) - f * d(j)
            enddo
          enddo
          d(m) = scale * d(m)
          a(m,m-1) = - scale * g
        endif
      enddo
      end
      subroutine ortran(ndim,n,ilow,iupp,h,d,v)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Accumulate the orthogonal similarity transformation used by        *
!   ORTHES to reduce a general real matrix A to upper Hessenberg form. *
!                                                                      *
!   Translation of the ALGOL procedure ORTRANS in:                     *
!   Handbook Series Linear Algebra,                                    *
!   Num. Math. 16, 181-204 (1970) by G. Peters and J. H. Wilkinson.    *
! Input:                                                               *
!   N         (integer) Order of the matrices A and V.                 *
!   ILOW,IUPP (integer) Determine a sub-matrix set by BALANC.          *
!                       May be set to 1 and N respectively.            *
!   H(NDIM,N) (real)    The matrix resulting from running ORTHES.      *
!   D(N)      (real)    Further information about the transformation.  *
! Output:                                                              *
!   V(NDIM,N) (real)    The accumulated transformation.                *
!   D(N)      (real)    Destroyed.                                     *
!----------------------------------------------------------------------*
      integer i,ilow,iupp,j,k,m,n,ndim
      double precision d(n),h(ndim,n),v(ndim,n),x,y,zero,one
      parameter(zero=0d0,one=1d0)

!---- Initialize V to identity matrix.
      do i = 1, n
        do j = 1, n
          v(i,j) = zero
        enddo
        v(i,i) = one
      enddo
!---- Accumulate transformations.
      do k = iupp - 2, ilow, - 1
        m = k + 1
        y = h(m,k)
        if (y .ne. zero) then
          y = y * d(m)
          do i = k + 2, iupp
            d(i) = h(i,k)
          enddo
          do j = m, iupp
            x = zero
            do i = m, iupp
              x = x + d(i) * v(i,j)
            enddo
            x = x / y
            do i = m, iupp
              v(i,j) = v(i,j) + x * d(i)
            enddo
          enddo
        endif
      enddo
      end
      subroutine hqr2(ndim,n,ilow,iupp,h,wr,wi,vecs,ierr)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Finds eigenvalues and eigenvectors of an unsymmetric real matrix,  *
!   A which has been reduced to upper Hessenberg form, H, by the       *
!   subroutine ORTHES. The orthogonal transformations must be placed   *
!   in the array VECS by subroutine ORTRAN.                            *
!                                                                      *
!   Translation of the ALGOL procedure HQR2 in:                        *
!   Handbook Series Linear Algebra,                                    *
!   Num. Math. 16, 181 - 204 (1970) by G. Peters and J. H. Wilkinson.  *
! Input:                                                               *
!   N         (integer) Order of the Hessenberg matrix H.              *
!   ILOW,IUPP (integer)                                                *
!   H(NDIM,N) (real)    The Hessenberg matrix produced by ORTHES.      *
!   VECS(NDIM,N) (real) A square matrix of order N containing the      *
!                       similarity transformation from A to H          *
! Output:                                                              *
!   H(NDIM,N) (real)    Modified.                                      *
!   WR(N)     (real)    Real parts of eigenvalues of H (or A).         *
!   WI(N)     (real)    Imaginary parts of eigenvalues of H (or A).    *
!   VECS(NDIM,N) (real) The unnormalized eigenvectors of A.            *
!                       Complex vectors are stored as pairs of reals.  *
!----------------------------------------------------------------------*
      integer i,ien,ierr,ilow,its,iupp,j,k,l,m,n,na,ndim
      double precision den,h(ndim,n),hnorm,p,q,r,ra,s,sa,t,temp,tempi,  &
     &tempr,vecs(ndim,n),vi,vr,w,wi(n),wr(n),x,y,z,epsmch,zero,one,two, &
     &triqua,fac1
      parameter(epsmch=1d-16,zero=0d0,one=1d0,two=2d0,triqua=.75d0,     &
     &fac1=.4375d0)

      ierr = 0
!---- Store isolated roots.
      do i = 1, n
        if (i .lt. ilow  .or.  i .gt. iupp) then
          wr(i) = h(i,i)
          wi(i) = zero
        endif
      enddo
      ien = iupp
      t = zero
!---- Next eigenvalue.
   60 if (ien .ge. ilow) then
        its = 0
        na = ien - 1
!---- Next iteration; look for single small sub-diagonal element.
   70   continue
        do l = ien, ilow + 1, -1
          if (abs(h(l,l-1)) .le.                                        &
     &    epsmch * (abs(h(l-1,l-1)) + abs(h(l,l)))) go to 100
        enddo
        l = ilow
  100   continue
        x = h(ien,ien)
        if (l .eq. ien) go to 270
        y = h(na,na)
        w = h(ien,na) * h(na,ien)
        if (l .eq. na) go to 280
        if (its .eq. 30) then
          ierr = ien
          go to 9999
        endif
!---- Form exceptional shift.
        if (its .eq. 10  .or.  its .eq. 20) then
          t = t + x
          do i = ilow, ien
            h(i,i) = h(i,i) - x
          enddo
          s = abs(h(ien,na)) + abs(h(na,ien-2))
          x = triqua * s
          y = x
          w = - fac1 * s * s
        endif
        its = its + 1
!---- Look for two consecutive small sub-diagonal elements.
        do m = ien - 2, l, - 1
          z = h(m,m)
          r = x - z
          s = y - z
          p = (r * s - w) / h(m+1,m) + h(m,m+1)
          q = h(m+1,m+1) - z - r - s
          r = h(m+2,m+1)
          s = abs(p) + abs(q) + abs(r)
          p = p / s
          q = q / s
          r = r / s
          if (m .eq. l) go to 150
          if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. epsmch * abs(p)    &
     &    * (abs(h(m-1,m-1)) + abs(z) + abs(h(m+1,m+1)))) go to 150
        enddo
  150   continue
        h(m+2,m) = zero
        do i = m + 3, ien
          h(i,i-2) = zero
          h(i,i-3) = zero
        enddo
!---- Double QR step involving rows L to IEN and columns M to IEN.
        do k = m, na
          if (k .ne. m) then
            p = h(k,k-1)
            q = h(k+1,k-1)
            if (k .ne. na) then
              r = h(k+2,k-1)
            else
              r = zero
            endif
            x = abs(p) + abs(q) + abs(r)
            if (x .eq. zero) go to 260
            p = p / x
            q = q / x
            r = r / x
          endif
          s = sign(sqrt(p**2+q**2+r**2),p)
          if (k .ne. m) then
            h(k,k-1) = - s * x
          else if (l .ne. m) then
            h(k,k-1) = - h(k,k-1)
          endif
          p = p + s
          x = p / s
          y = q / s
          z = r / s
          q = q / p
          r = r / p
!---- Row modification.
          do j = k, n
            p = h(k,j) + q * h(k+1,j)
            if (k .ne. na) then
              p = p + r * h(k+2,j)
              h(k+2,j) = h(k+2,j) - p * z
            endif
            h(k+1,j) = h(k+1,j) - p * y
            h(k,j) = h(k,j) - p * x
          enddo
!---- Column modification.
          j = min(ien,k+3)
          do i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            if (k .ne. na) then
              p = p + z * h(i,k+2)
              h(i,k+2) = h(i,k+2) - p * r
            endif
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k) = h(i,k) - p
          enddo
!---- Accumulate transformations.
          do i = ilow, iupp
            p = x * vecs(i,k) + y * vecs(i,k+1)
            if (k .ne. na) then
              p = p + z * vecs(i,k+2)
              vecs(i,k+2) = vecs(i,k+2) - p * r
            endif
            vecs(i,k+1) = vecs(i,k+1) - p * q
            vecs(i,k) = vecs(i,k) - p
          enddo
  260     continue
        enddo
!---- Go to next iteration.
        go to 70
!==== One real root found.
  270   h(ien,ien) = x + t
        wr(ien) = h(ien,ien)
        wi(ien) = zero
        ien = na
        go to 60
!==== Two roots (real pair or complex conjugate) found.
  280   p = (y - x) / two
        q = p**2 + w
        z = sqrt(abs(q))
        x = x + t
        h(ien,ien) = x
        h(na,na) = y + t
!---- Real pair.
        if (q .gt. zero) then
          z = p + sign(z,p)
          wr(na) = x + z
          wr(ien) = x - w / z
          wi(na) = zero
          wi(ien) = zero
          x = h(ien,na)
          r = sqrt(x**2+z**2)
          p = x / r
          q = z / r
!---- Row modification.
          do j = na, n
            z = h(na,j)
            h(na,j) = q * z + p * h(ien,j)
            h(ien,j) = q * h(ien,j) - p * z
          enddo
!---- Column modification.
          do i = 1, ien
            z = h(i,na)
            h(i,na) = q * z + p * h(i,ien)
            h(i,ien) = q * h(i,ien) - p * z
          enddo
!---- Accumulate transformations.
          do i = ilow, iupp
            z = vecs(i,na)
            vecs(i,na) = q * z + p * vecs(i,ien)
            vecs(i,ien) = q * vecs(i,ien) - p * z
          enddo
!---- Complex pair.
        else
          wr(na) = x + p
          wr(ien) = x + p
          wi(na) = z
          wi(ien) = -z
        endif
!----- Go to next root.
        ien = ien - 2
        go to 60
      endif
!==== Compute matrix norm.
      hnorm = zero
      k = 1
      do i = 1, n
        do j = k, n
          hnorm = hnorm + abs(h(i,j))
        enddo
        k = i
      enddo
!==== Back substitution.
      do ien = n, 1, -1
        p = wr(ien)
        q = wi(ien)
        na = ien - 1
!---- Real vector.
        if (q .eq. zero) then
          m = ien
          h(ien,ien) = one
          do i = na, 1, -1
            w = h(i,i) - p
            r = h(i,ien)
            do j = m, na
              r = r + h(i,j) * h(j,ien)
            enddo
            if (wi(i) .lt. zero) then
              z = w
              s = r
            else
              m = i
              if (wi(i) .eq. zero) then
                temp = w
                if (w .eq. zero) temp = epsmch * hnorm
                h(i,ien) = - r / temp
              else
                x = h(i,i+1)
                y = h(i+1,i)
                q = (wr(i) - p)**2 + wi(i)**2
                t = (x * s - z * r) / q
                h(i,ien) = t
                if (abs(x) .gt. abs(z)) then
                  h(i+1,ien) = - (r + w * t) / x
                else
                  h(i+1,ien) = - (s + y * t) / z
                endif
              endif
            endif
          enddo
!---- Complex vector associated with lamda = P - i * Q.
        else if (q .lt. zero) then
          m = na
          if (abs(h(ien,na)) .gt. abs(h(na,ien))) then
            h(na,na) = - (h(ien,ien) - p) / h(ien,na)
            h(na,ien) = - q / h(ien,na)
          else
            den = (h(na,na) - p)**2 + q**2
            h(na,na) = - h(na,ien) * (h(na,na) - p) / den
            h(na,ien) = h(na,ien) * q / den
          endif
          h(ien,na) = one
          h(ien,ien) = zero
          do i = ien - 2, 1, - 1
            w = h(i,i) - p
            ra = h(i,ien)
            sa = zero
            do j = m, na
              ra = ra + h(i,j) * h(j,na)
              sa = sa + h(i,j) * h(j,ien)
            enddo
            if (wi(i) .lt. zero) then
              z = w
              r = ra
              s = sa
            else
              m = i
              if (wi(i) .eq. zero) then
                den = w**2 + q**2
                h(i,na) = - (ra * w + sa * q) / den
                h(i,ien) = (ra * q - sa * w) / den
              else
                x = h(i,i+1)
                y = h(i+1,i)
                vr = (wr(i) - p)**2 + wi(i)**2 - q**2
                vi = two * (wr(i) - p) * q
                if (vr .eq. zero  .and.  vi .eq. zero) then
                  vr = epsmch * hnorm                                   &
     &            * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z))
                endif
                tempr = x * r - z * ra + q * sa
                tempi = x * s - z * sa - q * ra
                den = vr**2 + vi**2
                h(i,na) = (tempr * vr + tempi * vi) / den
                h(i,ien) = (tempi * vr - tempr * vi) / den
                if (abs(x) .gt. abs(z) + abs(q)) then
                  h(i+1,na) = (- ra - w * h(i,na) + q * h(i,ien)) / x
                  h(i+1,ien) = (- sa - w * h(i,ien) - q * h(i,na)) / x
                else
                  tempr = - r - y * h(i,na)
                  tempi = - s - y * h(i,ien)
                  den = z**2 + q**2
                  h(i+1,na) = (tempr * z + tempi * q) / den
                  h(i+1,ien) = (tempi * z - tempr * q) / den
                endif
              endif
            endif
          enddo
        endif
      enddo
!==== Vectors of isolated roots.
      do i = 1, n
        if (i .lt. ilow  .or.  i .gt. iupp) then
          do j = i, n
            vecs(i,j) = h(i,j)
          enddo
        endif
      enddo
!==== Multiply by transformation matrix to give eigenvectors of the
!     original full matrix.
      do j = n, ilow, - 1
        m = min(j,iupp)
        if (wi(j) .lt. zero) then
          l = j - 1
          do i = ilow, iupp
            y = zero
            z = zero
            do k = ilow, m
              y = y + vecs(i,k) * h(k,l)
              z = z + vecs(i,k) * h(k,j)
            enddo
            vecs(i,l) = y
            vecs(i,j) = z
          enddo
        else if (wi(j) .eq. zero) then
          do i = ilow, iupp
            z = zero
            do k = ilow, m
              z = z + vecs(i,k) * h(k,j)
            enddo
            vecs(i,j) = z
          enddo
        endif
      enddo
 9999 end

      subroutine suelem(el, ve, we,tilt)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Compute Displacement and rotation for one element.                 *
! Output:                                                              *
!   EL        (real)    Element length along design orbit.             *
!   VE(3)     (real)    Displacement of exit w.r.t. entry.             *
!   WE(3,3)   (real)    Rotation of exit w.r.t. entry.                 *
! Reference pointer used:                                              *
!   LCELM     /REFER/   Current element bank.                          *
! Local links:                                                         *
!   LSEQ                Beam lines sequence for a lump.                *
!----------------------------------------------------------------------*
! Modified: 28-DEC-1998, T. Raubenheimer (SLAC)                        *
!   Added LCAVITY element at ISP 27                                    *
!----------------------------------------------------------------------*
      include 'twtrr.fi'
      integer code,nn
      double precision angle,cospsi,costhe,ds,dx,sinpsi,sinthe,tilt,    &
     &ve(3),we(3,3),node_value,el,normal(0:maxmul),skew(0:maxmul)       &
     &,zero,one
      parameter(zero=0d0,one=1d0)
!---- Branch on subprocess code.
      tilt = zero
      angle = zero
      code = node_value('mad8_type ')
      go to ( 10,  20,  20,  40,  50,  60,  70,  80,  90, 100,          &
     &110, 120, 130, 140, 150, 160, 170, 180, 190, 200,                 &
     &210, 220, 230, 240, 250,  20, 270, 280, 290, 300,                 &
     &310, 310, 310, 310, 310, 310, 310, 310, 310, 310), code

!---- Drift space.
   10 continue
!---- Arbitrary matrix.
   40 continue
!---- Quadrupole.
   50 continue

!---- Sextupole.
   60 continue

!---- Octupole.
   70 continue

!---- Solenoid.
   90 continue

!---- RF cavity.
  100 continue

!---- Electrostatic separator.
  110 continue

!---- Kickers.
  140 continue
  150 continue
  160 continue

!---- Monitors.
  170 continue
  180 continue
  190 continue

!---- Apertures.
  200 continue
  210 continue

!---- Marker.
  250 continue

!---- Beam-beam.
  220 continue

!---- lcavity
  270 continue
!---- Reserved.
  280 continue
  290 continue
  300 continue

!---- Beam instrument.
  240 continue

!---- Lump.
  230 continue

!---- User-defined elements.
  310 continue

!****** end of straight elements ***************
      ve(1) = 0
      ve(2) = 0
      ve(3) = el
      we(1,1) = one
      we(2,1) = 0
      we(3,1) = 0
      we(1,2) = 0
      we(2,2) = one
      we(3,2) = 0
      we(1,3) = 0
      we(2,3) = 0
      we(3,3) = one
      go to 500

!---- multipoles , introduced  17.09.02 / AV
   80 continue
      call dzero(normal,maxmul+1)
      call dzero(skew,maxmul+1)
      call get_node_vector('knl ',nn,normal)
!      call get_node_vector('ksl ',ns,skew)
!      print *,"mult ",code,"  angle",normal(0),"  skew ",ns
!-----  dipole_bv introduced to suppress SU in MADX input (AV  7.10.02)
      angle = normal(0)*node_value('dipole_bv ')
      if (abs(angle) .lt. 1d-13) then
        tilt = zero
      else
        tilt =  node_value('tilt ')
      endif
! As el=0, there is no dx and no ds
      dx = zero
      ds = zero
      go to 490

!---- Any kind of  bend. 
   20 continue
!--------------  dipole_bv introduced to suppress SU (AV  7.10.02)
      angle = node_value('angle ')*node_value('dipole_bv ')
!      print *,"SUELEM dipole : angle =",angle
      if (abs(angle) .lt. 1d-13) then
        dx = zero
        ds = el
        tilt = zero
      else
        tilt =  node_value('tilt ')
!      print *,"SUELEM dipole : tilt =",tilt," length= ",el
        dx = el * (cos(angle)-one)/angle
        ds = el * sin(angle)/angle
      endif
!      print *,"SUELEM dipole : tilt =",tilt," length= ",
!     &el," angv = ",angv," bv =",node_value('dipole_bv ')
      go to 490

!---- Rotation around S-axis.
  120 continue
!        call ucopy(q(lcelm+meangr), angle, mwflt)
      we(1,1) = cos(angle)
      we(2,1) = sin(angle)
      we(1,2) = - we(2,1)
      we(2,2) = + we(1,1)
      go to 500

!---- Rotation around Y-axis.
  130 continue
      angle=0
      we(1,1) = cos(angle)
      we(3,1) = sin(angle)
      we(1,3) = - we(3,1)
      we(3,3) = + we(1,1)
      go to 500

!---- Common for bends and multipoles: Displacement and rotation matrix.
  490 continue
      cospsi = cos(tilt)
      sinpsi = sin(tilt)
      costhe = cos(angle)
      sinthe = sin(angle)
      ve(1) = dx * cospsi
      ve(2) = dx * sinpsi
      ve(3) = ds
      we(1,1) = costhe * cospsi*cospsi + sinpsi*sinpsi
      we(2,1) = (costhe - one) * cospsi * sinpsi
      we(3,1) = sinthe * cospsi
      we(1,2) = we(2,1)
      we(2,2) = costhe * sinpsi*sinpsi + cospsi*cospsi
      we(3,2) =  sinthe * sinpsi
      we(1,3) = - we(3,1)
      we(2,3) = - we(3,2)
      we(3,3) = costhe
  500 continue
      end
!-----------------  end of suelem subroutine --------------------------

!**********************************************************************
      subroutine sumtrx(the, phi, psi, w)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Given three survey angles, compute rotation matrix.                *
! Input:                                                               *
!   THE       (real)    Azimuthal angle.                               *
!   PHI       (real)    Elevation angle.                               *
!   PSI       (real)    Roll angle.                                    *
! Output:                                                              *
!   W(3,3)    (real)    Rotation matrix.                               *
!----------------------------------------------------------------------*
      double precision cosphi,cospsi,costhe,phi,psi,sinphi,sinpsi,      &
     &sinthe,the,w(3,3)

      costhe = cos(the)
      sinthe = sin(the)
      cosphi = cos(phi)
      sinphi = sin(phi)
      cospsi = cos(psi)
      sinpsi = sin(psi)
      w(1,1) = + costhe * cospsi - sinthe * sinphi * sinpsi
      w(1,2) = - costhe * sinpsi - sinthe * sinphi * cospsi
      w(1,3) =                     sinthe * cosphi
      w(2,1) =                              cosphi * sinpsi
      w(2,2) =                              cosphi * cospsi
      w(2,3) =                              sinphi
      w(3,1) = - sinthe * cospsi - costhe * sinphi * sinpsi
      w(3,2) = + sinthe * sinpsi - costhe * sinphi * cospsi
      w(3,3) =                     costhe * cosphi

      end
!-----------------  end of sumtrx subroutine --------------------------

!**********************************************************************
      subroutine sutran(w, v, we)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Transform rotation W and displacement V from entrance to exit.     *
! Input:                                                               *
!   W(3,3)    (real)    Rotation matrix w.r.t. input system.           *
!   V(3)      (real)    Displacement w.r.t. input system.              *
!   WE(3,3)   (real)    Rotation matrix due to element.                *
! Output:                                                              *
!   W(3,3)    (real)    Rotation matrix w.r.t. output system.          *
!   V(3)      (real)    Displacement w.r.t. output system.             *
!----------------------------------------------------------------------*
      integer i,k
      double precision v(3),vt(3),w(3,3),we(3,3),wt(3,3)

!---- VT := transpose(WE) * V;
!     WT := transpose(WE) * W;
      do i = 1, 3
        vt(i) = we(1,i)*v(1) + we(2,i)*v(2) + we(3,i)*v(3)
        do k = 1, 3
          wt(i,k) = we(1,i)*w(1,k) + we(2,i)*w(2,k) + we(3,i)*w(3,k)
        enddo
      enddo

!---- V := VT       [= transpose(WE) * V];
!     W := WT * WE  [= transpose(WE) * W * WE];
      do i = 1, 3
        v(i) = vt(i)
        do k = 1, 3
          w(i,k) = wt(i,1)*we(1,k) + wt(i,2)*we(2,k) + wt(i,3)*we(3,k)
        enddo
      enddo
      end
!-----------------  end of sutran subroutine --------------------------
      integer function lastnb(t)
!----------------------------------------------------------------------*
! Purpose:
!   Find last non-blank in string
!
!----------------------------------------------------------------------*
      implicit none
      integer mcnam, maxpnt
      parameter (mcnam = 16, maxpnt = 500)

      character *(*) t
      integer i
      do i = len(t), 1, -1
        if (t(i:i) .ne. ' ') goto 20
      enddo
      i = 1
   20 lastnb = i
      end
