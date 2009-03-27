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
