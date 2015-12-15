subroutine photon(iele,rad,dlength,energy,ieave,iquasto,d1,d2)
  !SUBROUTINE photon (i_elem_type, rad_curv_m, length_curr_elem, &
  !                        Energy_beam_MeV, mass_rest_MeV, d1, d2)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     subroutine to determine the energy losses at entrance
  !     and exit of magnet (dipole or quadrupole)
  !
  !
  !     input:
  !            iele - Magnet Type, iele=2 means a Quadrupole and
  !                   its Radiation can be switched off with iquasto flag
  !                   set in the Initialisation
  !            rad - radius of curvature
  !                  (rho for bend
  !                   1/(abs(k1)*sqrt(x**2+y**2)) for quadrupole
  !
  !            dlength - length of element
  !            energy - beam energy in MeV
  !     output:
  !            d1 - momentum loss at entrance (normalized to beam energy)
  !            d2 - momentum loss at exit (normalized to beam energy)
  !
  !                                                           FZ/28.02.2001
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  USE precision_constants, ONLY:CLIGHT,dhbar,pmae,dp,sp,qelect,c_1d6,c_1d_10,c_1d_3,c_0_3079,c_137,one,two,three,five
  implicit none
  include 'photoni.inc'

  INTEGER, INTENT(IN) :: iele, ieave, iquasto
  REAL(dp), INTENT(IN) :: dlength,energy
  REAL(dp), INTENT(OUT) :: d1,d2
  REAL(dp), INTENT(INOUT) :: rad

  REAL(dp) :: gamma,theta, dnpho, dnphohalf,ec,eave,ep1,ep2, ra2, &
       ran,xl,en,gamfac

  REAL(dp), parameter :: el=qelect, vl=CLIGHT
  real(sp) amu
  Double Precision :: ran2 ! function
  !      integer i,j,nmax,n1,n2,npmax,iseed,idumy,ieave,iquasto,iele
  !      double precision pi,vl,dhbar,el,gamma,energy,xilog,sollog,b,sol,  &
  !     &xil,sll,ak,bk,theta,dlength,rad,dnpho,dnphohalf,ec,eave,ep1,ep2,  &
  !     &d1,d2,ra2,ran,ran2,xl,en,gamfac
  !      parameter (nmax=100)

  !      common / lookup / xil(nmax), sll(nmax)
  !      common / lookup2 / ak(nmax), bk(nmax)
  !      common / rann / iseed, idumy
  !      common / trick / pi, vl, dhbar, el, gamma
  !      common / eave1 / ieave,iquasto

  !      real amu, amax
  !      integer ierr,nphoton
  INTEGER :: ierr, nphoton, n1, n2, npmax, i,j, idumy
  external rnpset, rnpssn
  ierr=0
  rad=abs(rad)
  !      print*,"RRRRRRRRRRRRRRRRRad in photon",iele,rad,dlength,energy

  ! for electrons
  gamma=energy/pmae*c_1d_3

  if (gamma.le.1) then
     write(*,*) ' error in subroutine photon - no initialization '
     return
  endif

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     Stop Quadrupole Radiation
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (iquasto.eq.0.and.iele.eq.2) return
  ! it is check already before CALL this subr.

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     compute mean number of photons emitted
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

989 continue

  !     trajectory bending angle
  theta = dlength/rad

  dnpho = five/two/sqrt(three)*gamma/c_137*theta
  !
  dnphohalf = dnpho/two
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     compute critical energy and average (half) energy loss
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     critical photon energy as fraction of beam energy

  ec = three/two*vl*gamma**3/rad*dhbar/(energy*el*c_1d6)
  !
  eave = dnphohalf*c_0_3079*ec
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     compute number of photons n1 and n2 actually emitted
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  amu = dnphohalf
  !      write(*,*) ' call rnpssn 1 ',amu,ierr
  call rnpssn(amu,nphoton,ierr)
  !      write(*,*) ' call rnpssn 2 ',amu,ierr,nphoton

  n1 = nphoton
  call rnpssn(amu,nphoton,ierr)
  !      write(*,*) ' after calling rnpssn ',amu,ierr,nphoton
  n2 = nphoton

  npmax = n1 + n2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     toss randum number to determine photon energy
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ep1 = 0.
  ep2 = 0.

  do i = 1, npmax

     ra2 = ran2(idumy)
     if (ra2.lt.c_1d_10) goto 989

     ran = log(ra2)
     do j = 1, nmax
        if (ran.lt.sll(j)) then
           if(j.eq.1) then
              xl = (ran-bk(1))/ak(1)
              goto 990
           endif
           xl = (ran-bk(j))/ak(j)
           goto 990
        endif
     enddo
     xl = (ran-bk(nmax))/ak(nmax)
990  continue

     !
     !     note: xl is logarithm of photon energy
     !     divided by critical energy
     !

     if (xl.lt.(-46.).or.(xl.gt.(2.3))) goto 989

     !     photon energy as fraction of beam energy
     en = -exp(xl)*ec

     if (i.le.n1) then
        ep1 = ep1 + en
     else
        ep2 = ep2 + en
     endif

  enddo

  !     compute total relative photon energy (average loss is
  !     added for compensation)
  !
  gamfac = gamma*gamma/(gamma*gamma-one)
  if( ieave.eq.1 ) then
     d1 = (ep1 + eave) * gamfac
     d2 = (ep2 + eave) * gamfac
  else
     d1 = ep1 * gamfac
     d2 = ep2 * gamfac
  endif

  return
end subroutine photon


double precision function   ran2(idum)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !   The Function Returns a Uniform Random Number between 0.0 and 1.0.
  !
  !   The Function uses a 'Subtractive Method'!!!
  !
  !   Set 'idum' to a negative value to initialize or reinitialize the Sequence.
  !
  !   (Numerical Recipies ran3(idumy), pg. 273.)
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  USE precision_constants, ONLY:one
  implicit none
  integer    idum
  integer    mbig, mseed, mz
  double precision     fac
  parameter(mbig=1000000000,mseed=161803398,mz=0,fac=one/mbig)
  integer    i, iff, ii, inext, inextp, k
  integer    mj, mk, ma(55)
  save       iff, inext, inextp, ma
  data       iff /0/

  !     Initialization:
  if(idum.lt.0.or.iff.eq.0) then
     iff = 1
     mj = mseed-iabs(idum)
     mj = mod(mj,mbig)
     mk = 1
     do i=1,54
        ii = mod(2*i,55)
        ma(ii) = mk
        mk = mj - mk
        if(mk.lt.mz) mk = mk + mbig
        mj = ma(ii)
     enddo

     do k=1,4
        do i=1,55
           ma(i) = ma(i) - ma(1+mod(i+30,55))
           if(ma(i).lt.mz) ma(i) = ma(i) + mbig
        enddo
     enddo
     inext = 0
     inextp = 31
     idum = 1
  endif

  inext = inext + 1
  if(inext.eq.56) inext = 1
  inextp = inextp + 1
  if(inextp.eq.56) inextp = 1
  mj = ma(inext) - ma(inextp)
  if(mj.lt.mz) mj = mj + mbig
  ma(inext) = mj
  ran2 = mj * fac

  return
END  function   ran2

! Following routines extracted from CERNLIB and used to be in poisson.f90

      SUBROUTINE RNPSSN(AMU,N,IERR)
      implicit none
      integer ierr,n,j
      real amu,AMAX,AMU0,EMU,AMXA,amx,r,p
      EXTERNAL RANLUX

      SAVE EMU,AMU0,AMAX

      DATA AMU0 /-12345.67/, AMAX /88/
      PARAMETER (AMXA = 88)
      dimension r(1)

      IERR=0
      IF(AMU .LE. 0) THEN
        IERR=1
        J=0
      ELSEIF(AMU .GT. AMAX) THEN
        CALL RNORMX(R,1,RANLUX)
        J=R(1)*SQRT(AMU)+AMU+0.5
      ELSE
        IF(AMU .NE. AMU0) THEN
          AMU0=AMU
          EMU=EXP(-AMU)
        ENDIF
        P=1
        J=-1
    1   J=J+1
        CALL RANLUX(R,1)
        P=P*R(1)
        IF(P .GT. EMU) GO TO 1
      ENDIF
      N=J
      RETURN

      ENTRY RNPSET(AMX)
      AMAX=MIN(AMX,AMXA)
      WRITE(6,'(/7X,''+++++ CERN V136 RNPSSN :  SWITCH TO '',           &
     &''NORMAL APPROXIMATION FOR      AMU > '',F7.2/)') AMAX
      RETURN
      END

      SUBROUTINE RNORMX(DEVIAS,NDEV,ROUTIN)
      implicit none
      integer ndev,IDEV
      real r1,r2,a,b,q,s,t,u,v,x,y,DEVIAS,DEVIAT
!        Generator of a vector of independent Gaussian-distributed
!        (pseudo-)random numbers, of mean zero and variance one,
!        making use of a uniform pseudo-random generator (RANMAR).
!        The algorithm for converting uniform numbers to Gaussian
!        is that of "Ratio of Uniforms with Quadratic Bounds."  The
!        method is in principle exact (apart from rounding errors),
!        and is based on the variant published by Joseph Leva in
!        ACM TOMS vol. 18(1992), page 449 for the method and 454 for
!        the Fortran algorithm (ACM No. 712).
!        It requires at least 2 and on average 2.74 uniform deviates
!        per Gaussian (normal) deviate.
!   WARNING -- The uniform generator should not produce exact zeroes,
!   since the pair (0.0, 0.5) provokes a floating point exception.
      SAVE  S, T, A, B, R1, R2
      DIMENSION U(2), DEVIAS(*)
      EXTERNAL ROUTIN
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2/ 0.27597, 0.27846/
!         generate pair of uniform deviates
      DO IDEV = 1, NDEV
   50   CALL ROUTIN(U,2)
        V = 1.7156 * (U(2) - 0.5)
        X = U(1) - S
        Y = ABS(V) - T
        Q = X**2 + Y*(A*Y - B*X)
!           accept P if inside inner ellipse
        IF (Q .LT. R1)  GO TO 100
!           reject P if outside outer ellipse
        IF (Q .GT. R2)  GO TO 50
!           reject P if outside acceptance region
        IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 50
!           ratio of P's coordinates is normal deviate
  100   DEVIAT = V/U(1)
        DEVIAS(IDEV) = DEVIAT
      enddo
      RETURN
      END

      SUBROUTINE RANLUX(RVEC,LENV)
      implicit none
      integer LENV,ISDEXT,ISEEDS,MAXLEV,NDSKIP,NEXT,ITWO24,J24,I24,     &
     &INSEED,MKOUNT,IN24,NSKIP,LXDFLT,JSDFLT,JSEED,LP,I,K,ICONS,IVEC,   &
     &ISK,IGIGA,ISD,K2,K1,INOUT,LOUT,INS,LUX,ILX,IOUTER,INNER,IZIP,     &
     &IZIP2,KOUNT
      real RVEC,UNI,TWOP12,TWOM12,TWOM24,CARRY,SEEDS
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
!                               default
!  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
!orresponds to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
        NOTYET = .FALSE.
        JSEED = JSDFLT
        INSEED = JSEED
        WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
        LUXLEV = LXDFLT
        NSKIP = NDSKIP(LUXLEV)
        LP = NSKIP + 24
        IN24 = 0
        KOUNT = 0
        MKOUNT = 0
        WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',     &
     &LUXLEV,'      p =',LP
        TWOM24 = 1.
        DO I= 1, 24
          TWOM24 = TWOM24 * 0.5
          K = JSEED/53668
          JSEED = 40014*(JSEED-K*53668) -K*12211
          IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
          ISEEDS(I) = MOD(JSEED,ITWO24)
        enddo
        TWOM12 = TWOM24 * 4096.
        DO I= 1,24
          SEEDS(I) = REAL(ISEEDS(I))*TWOM24
          NEXT(I) = I-1
        enddo
        NEXT(1) = 24
        I24 = 24
        J24 = 10
        CARRY = 0.
        IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
      DO IVEC= 1, LENV
        UNI = SEEDS(J24) - SEEDS(I24) - CARRY
        IF (UNI .LT. 0.)  THEN
          UNI = UNI + 1.0
          CARRY = TWOM24
        ELSE
          CARRY = 0.
        ENDIF
        SEEDS(I24) = UNI
        I24 = NEXT(I24)
        J24 = NEXT(J24)
        RVEC(IVEC) = UNI
!  small numbers (with less than 12 "significant" bits) are "padded".
        IF (UNI .LT. TWOM12)  THEN
          RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
!        and zero is forbidden in case someone takes a logarithm
          IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
        ENDIF
!        Skipping to luxury.  As proposed by Martin Luscher.
        IN24 = IN24 + 1
        IF (IN24 .EQ. 24)  THEN
          IN24 = 0
          KOUNT = KOUNT + NSKIP
          DO ISK= 1, NSKIP
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
              UNI = UNI + 1.0
              CARRY = TWOM24
            ELSE
              CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
          enddo
        ENDIF
      enddo
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
        MKOUNT = MKOUNT + 1
        KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
!
!           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
      NOTYET = .FALSE.
      TWOM24 = 1.
      DO I= 1, 24
        NEXT(I) = I-1
        TWOM24 = TWOM24 * 0.5
      enddo
      NEXT(1) = 24
      TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO I= 1, 24
        SEEDS(I) = REAL(ISDEXT(I))*TWOM24
      enddo
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
      IF (LUXLEV .LE. MAXLEV) THEN
        NSKIP = NDSKIP(LUXLEV)
        WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',   &
     &LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
        INSEED = -1
        RETURN
!
!                    Entry to ouput seeds as integers
        ENTRY RLUXUT(ISDEXT)
        DO I= 1, 24
          ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
        enddo
        ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
        IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
        RETURN
!
!                    Entry to output the "convenient" restart point
        ENTRY RLUXAT(LOUT,INOUT,K1,K2)
        LOUT = LUXLEV
        INOUT = INSEED
        K1 = KOUNT
        K2 = MKOUNT
        RETURN
!
!                    Entry to initialize from one or three integers
        ENTRY RLUXGO(LUX,INS,K1,K2)
        IF (LUX .LT. 0) THEN
          LUXLEV = LXDFLT
        ELSE IF (LUX .LE. MAXLEV) THEN
          LUXLEV = LUX
        ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
          LUXLEV = MAXLEV
          WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
        ELSE
          LUXLEV = LUX
          DO ILX= 0, MAXLEV
            IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
          enddo
        ENDIF
        IF (LUXLEV .LE. MAXLEV)  THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',&
     &LUXLEV,'     P=', NSKIP+24
        ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
        ENDIF
        IN24 = 0
        IF (INS .LT. 0)  WRITE (6,'(A)')                                &
     &' Illegal initialization by RLUXGO, negative input seed'
        IF (INS .GT. 0)  THEN
          JSEED = INS
          WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEED',&
     &'S',                                                              &
     &JSEED, K1,K2
        ELSE
          JSEED = JSDFLT
          WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SE',&
     &'ED'
        ENDIF
        INSEED = JSEED
        NOTYET = .FALSE.
        TWOM24 = 1.
        DO I= 1, 24
          TWOM24 = TWOM24 * 0.5
          K = JSEED/53668
          JSEED = 40014*(JSEED-K*53668) -K*12211
          IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
          ISEEDS(I) = MOD(JSEED,ITWO24)
        enddo
        TWOM12 = TWOM24 * 4096.
        DO I= 1,24
          SEEDS(I) = REAL(ISEEDS(I))*TWOM24
          NEXT(I) = I-1
        enddo
        NEXT(1) = 24
        I24 = 24
        J24 = 10
        CARRY = 0.
        IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
        KOUNT = K1
        MKOUNT = K2
        IF (K1+K2 .NE. 0)  THEN
          DO IOUTER= 1, K2+1
            INNER = IGIGA
            IF (IOUTER .EQ. K2+1)  INNER = K1
            DO ISK= 1, INNER
              UNI = SEEDS(J24) - SEEDS(I24) - CARRY
              IF (UNI .LT. 0.)  THEN
                UNI = UNI + 1.0
                CARRY = TWOM24
              ELSE
                CARRY = 0.
              ENDIF
              SEEDS(I24) = UNI
              I24 = NEXT(I24)
              J24 = NEXT(J24)
            enddo
          enddo
!         Get the right value of IN24 by direct calculation
          IN24 = MOD(KOUNT, NSKIP+24)
          IF (MKOUNT .GT. 0)  THEN
            IZIP = MOD(IGIGA, NSKIP+24)
            IZIP2 = MKOUNT*IZIP + IN24
            IN24 = MOD(IZIP2, NSKIP+24)
          ENDIF
!       Now IN24 had better be between zero and 23 inclusive
          IF (IN24 .GT. 23) THEN
            WRITE (6,'(A/A,3I11,A,I5)')                                 &
     &'  Error in RESTARTING with RLUXGO:','  The values', INS,         &
     &K1, K2, ' cannot occur at luxury level', LUXLEV
            IN24 = 0
          ENDIF
        ENDIF
        RETURN
        END
