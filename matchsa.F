      subroutine mtsa(ncon,nvar,tol,calls,call_lim,vect,fun_vect,iseed  &
     &,iprint,lb,nacp,ub,xopt,c,vm,xp)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   SA command.                                                        *
! Attributes :                                                         *
!   ncon      (int)     # constraints                                  *
!   nvar      (int)     # variables                                    *
!   iseed     (int)     initial seed                                   *
!   iprint    (int)     print flag  (0 = minimum info printed)         *
!   tol       (real)    Final tolerance for match.                     *
!   calls     (int)     current call count                             *
!   call_lim  (int)     current call limit                             *
!   vect      (real)    variable values                                *
!   fun_vect  (real)    function values                                *
!   all other           working spaces for SA                          *
!----------------------------------------------------------------------*
      integer calls,call_lim,ncon,nvar,nacp(*),neps,n,m,ns,nt,          &
     &iprint,iseed,iseed2,nacc,nobds,ier,i,iflag
      parameter(neps=4)
      double precision tol,t,vect(*),fun_vect(*),fopt,vmod
      double precision rt,lb(*),ub(*),xopt(*),c(*),                     &
     &vm(*),xp(*),fstar(neps)
      integer j,next_vary,get_option,name_l
      double precision get_variable,val,c_min,c_max,step
      logical logmax,psum
      common/forsa/m
      parameter(name_l=24)
!     include 'name_len.fi'
!     character*(name_len) name !uncomment when name_len.fi present
      character*(name_l) name   !clear
      n=nvar
      m=ncon
!     include 'match.fi'         !not used here i guess
      psum=get_option('match_summary ') .ne. 0
 1    continue
      j = next_vary(name,name_l,c_min,c_max,step)
      if (j .ne. 0)  then
        val = get_variable(name)
        lb(j)=c_min
        ub(j)=c_max
        c(j) =2.d0
        vm(j)=1.d0
        if(psum) write(*,830) name,val,c_min,c_max
        vect(j) = val
        goto 1
      endif
  830 format(a24,1x,1p,e16.8,3x,e16.8,3x,e16.8)
! Choose SA parameters next:
! The value of rt and ns as suggested by Corana
      rt=0.85d0
      ns = 20
! The value of nt suggested by Corana is max(100, 5*nvar)
      nt = max(100, 5*nvar)
!      iprint=1
!      iseed=1
      iseed2 = iseed+1
      logmax=.FALSE.
! The initial temeprature T is best to be an user-defined parameter
! set it to 10 here for brevery
      t=10.
      write(*,8800) n, t, rt, tol, ns, nt, neps, call_lim, iprint,      &
     &iseed, iseed2
8800  format(/,' simulated annealing ',/,                               &
     &/,' number of parameters: ',i3,                                   &
     &/,' initial temp: ',g8.2, ' rt: ',g8.2,' eps:',g8.2,              &
     &/,' ns: ',i3, '   nt: ',i4, '   neps: ',i2,                       &
     &/,' calls: ',i10, '   iprint: ',i1, '   iseed: ',i4,              &
     &'   iseed2: ',i4)
8001  format(/,' optimal function value: ',g20.13                       &
     &/,' number of function evaluations:     ',i10,                    &
     &/,' number of accepted evaluations:     ',i10,                    &
     &/,' number of out of bound evaluations: ',i10,                    &
     &/,' final temp: ', g20.13,'  ier: ', i3)
      call prtvec(vect,n,'starting values')
      call prtvec(vm,n,'initial step length')
      call prtvec(lb,n,'lower bound')
      call prtvec(ub,n,'upper bound')
      !call prtvec(c,n,'c vector')
      call sa(nvar,vect,logmax,rt,tol,ns,nt,neps,call_lim,lb,ub,        &
     &c,iprint,iseed,iseed2,t,vm,xopt,fopt,nacc,calls,nobds,ier,        &
     &fstar,xp,nacp)
      call prtvec(xopt,n,'solution')
      !call prtvec(vm,n,'final step length')
      write(*,8001) fopt, calls, nacc, nobds, t, ier
      do i=1,nvar
        vect(i)=xopt(i)
      enddo

      call mtfcn(ncon,nvar,vect,fun_vect,iflag)

      fopt = vmod(m,fun_vect)

      print *," fopt again = ", fopt

 9999 end

      subroutine fcnsa(n,x,f)
      implicit none
      integer iflag,nf,n
      double precision f,x(n),fvec(1000),vmod
      common/forsa/nf
      iflag=0
      call mtfcn(nf,n,x,fvec,iflag)
      f = vmod(nf, fvec)
!     f=0
!     do i=1,nf
!     f=f+fvec(i)**2
!     enddo
!      f = dsqrt(f)
      RETURN
      END

!* =====================================================================
! GAMS : Source of SIMANN in OPT from NETLIB
! ======================================================================
! NIST Guide to Available Math Software.
! Source for module SIMANN from package OPT.
! Retrieved from NETLIB on Wed Oct  2 21:25:58 1996.
! ======================================================================
! Logic slightly modified by Dobrin Kaltchev in Nov 1998.
! ABSTRACT:
!   Simulated annealing is a global optimization method that distinguishes
!   between different local optima. Starting from an initial point, the
!   algorithm takes a step and the function is evaluated. When minimizing a
!   function, any downhill step is accepted and the process repeats from this
!   new point. An uphill step may be accepted. Thus, it can escape from local
!   optima. This uphill decision is made by the Metropolis criteria. As the
!   optimization process proceeds, the length of the steps decline and the
!   algorithm closes in on the global optimum. Since the algorithm makes very
!   few assumptions regarding the function to be optimized, it is quite
!   robust with respect to non-quadratic surfaces. The degree of robustness
!   can be adjusted by the user. In fact, simulated annealing can be used as
!   a local optimizer for difficult functions.
!
!   This implementation of simulated annealing was used in "Global Optimization
!   of Statistical Functions with Simulated Annealing," Goffe, Ferrier and
!   Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp.
!   65-100. Briefly, we found it competitive, if not superior, to multiple
!   restarts of conventional optimization routines for difficult optimization
!   problems.
!
!   For more information on this routine, contact its author:
!   Bill Goffe, bgoffe@whale.st.usm.edu
!
!      PROGRAM SIMANN
!  This file is an example of the Corana et al. simulated annealing
!  algorithm for multimodal and robust optimization as implemented
!  and modified by Goffe, Ferrier and Rogers. Counting the above line
!  ABSTRACT as 1, the routine itself (SA), with its supplementary
!  routines, is on lines 232-990. A multimodal example from Judge et al.
!  (FCN) is on lines 150-231. The rest of this file (lines 1-149) is a
!  driver routine with values appropriate for the Judge example. Thus, this
!  example is ready to run.
!
!  To understand the algorithm, the documentation for SA on lines 236-
!  484 should be read along with the parts of the paper that describe
!  simulated annealing. Then the following lines will then aid the user
!  in becomming proficient with this implementation of simulated
!  annealing.
!
!  Learning to use SA:
!      Use the sample function from Judge with the following suggestions
!  to get a feel for how SA works. When you've done this, you should be
!  ready to use it on most any function with a fair amount of expertise.
!    1. Run the program as is to make sure it runs okay. Take a look at
!       the intermediate output and see how it optimizes as temperature
!       (T) falls. Notice how the optimal point is reached and how
!       falling T reduces VM.
!    2. Look through the documentation to SA so the following makes a
!       bit of sense. In line with the paper, it shouldn't be that hard
!       to figure out. The core of the algorithm is described on pp. 68-70
!       and on pp. 94-95. Also see Corana et al. pp. 264-9.
!    3. To see how it selects points and makes decisions about uphill
!       and downhill moves, set IPRINT = 3 (very detailed intermediate
!       output) and MAXEVL = 100 (only 100 function evaluations to limit
!       output).
!    4. To see the importance of different temperatures, try starting
!       with a very low one (say T = 10E-5). You'll see (i) it never
!       escapes from the local optima (in annealing terminology, it
!       quenches) & (ii) the step length (VM) will be quite small. This
!       is a key part of the algorithm: as temperature (T) falls, step
!       length does too. In a minor point here, note how VM is quickly
!       reset from its initial value. Thus, the input VM is not very
!       important. This is all the more reason to examine VM once the
!       algorithm is underway.
!    5. To see the effect of different parameters and their effect on
!       the speed of the algorithm, try RT = .95 & RT = .1. Notice the
!       vastly different speed for optimization. Also try NT = 20. Note
!       that this sample function is quite easy to optimize, so it will
!       tolerate big changes in these parameters. RT and NT are the
!       parameters one should adjust to modify the runtime of the
!       algorithm and its robustness.
!    6. Try constraining the algorithm with either LB or UB.

      SUBROUTINE SA(N,X,MAX,RT,EPS,NS,NT,NEPS,MAXEVL,LB,UB,C,IPRINT,    &
     &ISEED1,ISEED2,T,VM,XOPT,FOPT,NACC,NFCNEV,NOBDS,IER,               &
     &FSTAR,XP,NACP)
      implicit none

!  Version: 3.2
!  Date: 1/22/94.
!  Differences compared to Version 2.0:
!     1. If a trial is out of bounds, a point is randomly selected
!        from LB(i) to UB(i). Unlike in version 2.0, this trial is
!        evaluated and is counted in acceptances and rejections.
!        All corresponding documentation was changed as well.
!  Differences compared to Version 3.0:
!     1. If VM(i) > (UB(i) - LB(i)), VM is set to UB(i) - LB(i).
!        The idea is that if T is high relative to LB & UB, most
!        points will be accepted, causing VM to rise. But, in this
!        situation, VM has little meaning; particularly if VM is
!        larger than the acceptable region. Setting VM to this size
!        still allows all parts of the allowable region to be selected.
!  Differences compared to Version 3.1:
!     1. Test made to see if the initial temperature is positive.
!     2. WRITE statements prettied up.
!     3. References to paper updated.
!
!  Synopsis:
!  This routine implements the continuous simulated annealing global
!  optimization algorithm described in Corana et al.'s article
!  "Minimizing Multimodal Functions of Continuous Variables with the
!  "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
!  no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
!  Software.
!
!  A very quick (perhaps too quick) overview of SA:
!     SA tries to find the global optimum of an N dimensional function.
!  It moves both up and downhill and as the optimization process
!  proceeds, it focuses on the most promising area.
!     To start, it randomly chooses a trial point within the step length
!  VM (a vector of length N) of the user selected starting point. The
!  function is evaluated at this trial point and its value is compared
!  to its value at the initial point.
!     In a maximization problem, all uphill moves are accepted and the
!  algorithm continues from that trial point. Downhill moves may be
!  accepted; the decision is made by the Metropolis criteria. It uses T
!  (temperature) and the size of the downhill move in a probabilistic
!  manner. The smaller T and the size of the downhill move are, the more
!  likely that move will be accepted. If the trial is accepted, the
!  algorithm moves on from that point. If it is rejected, another point
!  is chosen instead for a trial evaluation.
!     Each element of VM periodically adjusted so that half of all
!  function evaluations in that direction are accepted.
!     A fall in T is imposed upon the system with the RT variable by
!  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
!  downhill moves are less likely to be accepted and the percentage of
!  rejections rise. Given the scheme for the selection for VM, VM falls.
!  Thus, as T declines, VM falls and SA focuses upon the most promising
!  area for optimization.
!
!  The importance of the parameter T:
!     The parameter T is crucial in using SA successfully. It influences
!  VM, the step length over which the algorithm searches for optima. For
!  a small intial T, the step length may be too small; thus not enough
!  of the function might be evaluated to find the global optima. The user
!  should carefully examine VM in the intermediate output (set IPRINT =
!  1) to make sure that VM is appropriate. The relationship between the
!  initial temperature and the resulting step length is function
!  dependent.
!     To determine the starting temperature that is consistent with
!  optimizing a function, it is worthwhile to run a trial run first. Set
!  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
!  rises as well. Then select the T that produces a large enough VM.
!
!  For modifications to the algorithm and many details on its use,
!  (particularly for econometric applications) see Goffe, Ferrier
!  and Rogers, "Global Optimization of Statistical Functions with
!  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
!  Jan./Feb. 1994, pp. 65-100.
!  For more information, contact
!              Bill Goffe
!              Department of Economics and International Business
!              University of Southern Mississippi
!              Hattiesburg, MS  39506-5072
!              (601) 266-4484 (office)
!              (601) 266-4920 (fax)
!              bgoffe@whale.st.usm.edu (Internet)
!
!  As far as possible, the parameters here have the same name as in
!  the description of the algorithm on pp. 266-8 of Corana et al.
!
!  In this description, SP is single precision, DP is double precision,
!  INT is integer, L is logical and (N) denotes an array of length n.
!  Thus, DP(N) denotes a double precision array of length n.
!
!  Input Parameters:
!    Note: The suggested values generally come from Corana et al. To
!          drastically reduce runtime, see Goffe et al., pp. 90-1 for
!          suggestions on choosing the appropriate RT and NT.
!    N - Number of variables in the function to be optimized. (INT)
!    X - The starting values for the variables of the function to be
!        optimized. (DP(N))
!    MAX - Denotes whether the function should be maximized or
!          minimized. A true value denotes maximization while a false
!          value denotes minimization. Intermediate output (see IPRINT)
!          takes this into account. (L)
!    RT - The temperature reduction factor. The value suggested by
!         Corana et al. is .85. See Goffe et al. for more advice. (DP)
!    EPS - Error tolerance for termination. If the final function
!          values from the last neps temperatures differ from the
!          corresponding value at the current temperature by less than
!          EPS and the final function value at the current temperature
!          differs from the current optimal function value by less than
!          EPS, execution terminates and IER = 0 is returned. (EP)
!    NS - Number of cycles. After NS*N function evaluations, each
!         element of VM is adjusted so that approximately half of
!         all function evaluations are accepted. The suggested value
!         is 20. (INT)
!    NT - Number of iterations before temperature reduction. After
!         NT*NS*N function evaluations, temperature (T) is changed
!         by the factor RT. Value suggested by Corana et al. is
!         MAX(100, 5*N). See Goffe et al. for further advice. (INT)
!    NEPS - Number of final function values used to decide upon termi-
!           nation. See EPS. Suggested value is 4. (INT)
!    MAXEVL - The maximum number of function evaluations. If it is
!             exceeded, IER = 1. (INT)
!    LB - The lower bound for the allowable solution variables. (DP(N))
!    UB - The upper bound for the allowable solution variables. (DP(N))
!         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
!         I = 1, N, a point is from inside is randomly selected. This
!         This focuses the algorithm on the region inside UB and LB.
!         Unless the user wishes to concentrate the search to a par-
!         ticular region, UB and LB should be set to very large positive
!         and negative values, respectively. Note that the starting
!         vector X should be inside this region. Also note that LB and
!         UB are fixed in position, while VM is centered on the last
!         accepted trial set of variables that optimizes the function.
!    C - Vector that controls the step length adjustment. The suggested
!        value for all elements is 2.0. (DP(N))
!    IPRINT - controls printing inside SA. (INT)
!             Values: 0 - Nothing printed.
!                     1 - Function value for the starting value and
!                         summary results before each temperature
!                         reduction. This includes the optimal
!                         function value found so far, the total
!                         number of moves (broken up into uphill,
!                         downhill, accepted and rejected), the
!                         number of out of bounds trials, the
!                         number of new optima found at this
!                         temperature, the current optimal X and
!                         the step length VM. Note that there are
!                         N*NS*NT function evalutations before each
!                         temperature reduction. Finally, notice is
!                         is also given upon achieveing the termination
!                         criteria.
!                     2 - Each new step length (VM), the current optimal
!                         X (XOPT) and the current trial X (X). This
!                         gives the user some idea about how far X
!                         strays from XOPT as well as how VM is adapting
!                         to the function.
!                     3 - Each function evaluation, its acceptance or
!                         rejection and new optima. For many problems,
!                         this option will likely require a small tree
!                         if hard copy is used. This option is best
!                         used to learn about the algorithm. A small
!                         value for MAXEVL is thus recommended when
!                         using IPRINT = 3.
!             Suggested value: 1
!             Note: For a given value of IPRINT, the lower valued
!                   options (other than 0) are utilized.
!    ISEED1 - The first seed for the random number generator RANMAR.
!             0 .LE. ISEED1 .LE. 31328. (INT)
!    ISEED2 - The second seed for the random number generator RANMAR.
!             0 .LE. ISEED2 .LE. 30081. Different values for ISEED1
!             and ISEED2 will lead to an entirely different sequence
!             of trial points and decisions on downhill moves (when
!             maximizing). See Goffe et al. on how this can be used
!             to test the results of SA. (INT)
!
!  Input/Output Parameters:
!    T - On input, the initial temperature. See Goffe et al. for advice.
!        On output, the final temperature. (DP)
!    VM - The step length vector. On input it should encompass the
!         region of interest given the starting value X. For point
!         X(I), the next trial point is selected is from X(I) - VM(I)
!         to  X(I) + VM(I). Since VM is adjusted so that about half
!         of all points are accepted, the input value is not very
!         important (i.e. is the value is off, SA adjusts VM to the
!         correct value). (DP(N))
!
!  Output Parameters:
!    XOPT - The variables that optimize the function. (DP(N))
!    FOPT - The optimal value of the function. (DP)
!    NACC - The number of accepted function evaluations. (INT)
!    NFCNEV - The total number of function evaluations. In a minor
!             point, note that the first evaluation is not used in the
!             core of the algorithm; it simply initializes the
!             algorithm. (INT).
!    NOBDS - The total number of trial function evaluations that
!            would have been out of bounds of LB and UB. Note that
!            a trial point is randomly selected between LB and UB.
!            (INT)
!    IER - The error return number. (INT)
!          Values: 0 - Normal return; termination criteria achieved.
!                  1 - Number of function evaluations (NFCNEV) is
!                      greater than the maximum number (MAXEVL).
!                  2 - The starting value (X) is not inside the
!                      bounds (LB and UB).
!                  3 - The initial temperature is not positive.
!                  99 - Should not be seen; only used internally.
!
!  Work arrays that must be dimensioned in the calling routine:
!       RWK1 (DP(NEPS))  (FSTAR in SA)
!       RWK2 (DP(N))     (XP    "  " )
!       IWK  (INT(N))    (NACP  "  " )
!
!  Required Functions (included):
!    EXPREP - Replaces the function EXP to avoid under- and overflows.
!             It may have to be modified for non IBM-type main-
!             frames. (DP)
!    RMARIN - Initializes the random number generator RANMAR.
!    RANMAR - The actual random number generator. Note that
!             RMARIN must run first (SA does this). It produces uniform
!             random numbers on [0,1]. These routines are from
!             Usenet's comp.lang.fortran. For a reference, see
!             "Toward a Universal Random Number Generator"
!             by George Marsaglia and Arif Zaman, Florida State
!             University Report: FSU-SCRI-87-50 (1987).
!             It was later modified by F. James and published in
!             "A Review of Pseudo-random Number Generators." For
!             further information, contact stuart@ads.com. These
!             routines are designed to be portable on any machine
!             with a 24-bit or more mantissa. I have found it produces
!             identical results on a IBM 3081 and a Cray Y-MP.
!
!  Required Subroutines (included):
!    PRTVEC - Prints vectors.
!    PRT1 ... PRT10 - Prints intermediate output.
!    FCNSA - Function to be optimized. The form is
!            SUBROUTINE FCNSA(N,X,F)
!            INTEGER N
!            DOUBLE PRECISION  X(N), F
!            ...
!            function code with F = F(X)
!            ...
!            RETURN
!            END
!          Note: This is the same form used in the multivariable
!          minimization algorithms in the IMSL edition 10 library.
!
!  Machine Specific Features:
!    1. EXPREP may have to be modified if used on non-IBM type main-
!       frames. Watch for under- and overflows in EXPREP.
!    2. Some FORMAT statements use G25.18; this may be excessive for
!       some machines.
!    3. RMARIN and RANMAR are designed to be protable; they should not
!       cause any problems.

!  Type all external variables.
      DOUBLE PRECISION  X(*), LB(*), UB(*), C(*), VM(*), FSTAR(*),      &
     &XOPT(*), XP(*), T, EPS, RT, FOPT
      INTEGER  NACP(*), N, NS, NT, NEPS, NACC, MAXEVL, IPRINT,          &
     &NOBDS, IER, NFCNEV, ISEED1, ISEED2
      LOGICAL  MAX
!  Type all internal variables.
      DOUBLE PRECISION  F, FP, P, PP, RATIO
      INTEGER  NUP, NDOWN, NREJ, NNEW, LNOBDS, H, I, J, M
      LOGICAL  QUIT

!  Type all functions.
      DOUBLE PRECISION  EXPREP
      REAL  RANMAR

!  Initialize the random number generator RANMAR.
      CALL RMARIN(ISEED1,ISEED2)
!  Set initial values.
      NACC = 0
      NOBDS = 0
      NFCNEV = 0
      IER = 99

      DO I = 1, N
        XOPT(I) = X(I)
        NACP(I) = 0
      enddo

      DO I = 1, NEPS
        FSTAR(I) = 1.0D+20
      enddo

!  If the initial temperature is not positive, notify the user and
!  return to the calling routine.
      IF (T .LE. 0.0) THEN
        WRITE(*,'(/,''  THE INITIAL TEMPERATURE IS NOT POSITIVE.'')')
!         '',
!     1             /,''  RESET THE VARIABLE T. ''/)')
        IER = 3
        RETURN
      END IF

!  If the initial value is out of bounds, notify the user and return
!  to the calling routine.
      DO I = 1, N
        IF ((X(I) .GT. UB(I)) .OR. (X(I) .LT. LB(I))) THEN
          CALL PRT1
          IER = 2
          RETURN
        END IF

      enddo

!  Evaluate the function with input X and return value as F.

      CALL FCNSA(N,X,F)

!  If the function is to be minimized, switch the sign of the function.
!  Note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.
      IF(.NOT. MAX) F = -F
      NFCNEV = NFCNEV + 1
      FOPT = F
      FSTAR(1) = F
      IF(IPRINT .GE. 1) CALL PRT2(MAX,N,X,F)

!  Start the main loop. Note that it terminates if (i) the algorithm
!  succesfully optimizes the function or (ii) there are too many
!  function evaluations (more than MAXEVL).
100   NUP = 0
      NREJ = 0
      NNEW = 0
      NDOWN = 0
      LNOBDS = 0

      DO M = 1,NT
        DO J = 1, NS
          DO H = 1, N

!  Generate XP, the trial value of X. Note use of VM to choose XP.
            DO I = 1, N
              IF (I .EQ. H) THEN
999             XP(I) = X(I) + (RANMAR()*2.- 1.) * VM(I)
              ELSE
                XP(I) = X(I)
              END IF

!  If XP is out of bounds, select a point in bounds for the trial.

              IF((XP(I) .LT. LB(I)) .OR. (XP(I) .GT. UB(I))) THEN
                XP(I) = LB(I) + (UB(I) - LB(I))*RANMAR()
                LNOBDS = LNOBDS + 1
                NOBDS = NOBDS + 1
                IF(IPRINT .GE. 3) CALL PRT3(MAX,N,XP,X,FP,F)
              END IF

            enddo

            IF((XP(h) .LT. LB(h)) .OR. (XP(h) .GT. UB(h))) THEN
              XP(h) = LB(h) + (UB(h) - LB(h))*RANMAR()
              LNOBDS = LNOBDS + 1
              NOBDS = NOBDS + 1
              IF(IPRINT .GE. 3) CALL PRT3(MAX,N,XP,X,FP,F)
            END IF


!  Evaluate the function with the trial point XP and return as FP.
            CALL FCNSA(N,XP,FP)
            IF(.NOT. MAX) FP = -FP
            NFCNEV = NFCNEV + 1
            IF(IPRINT .GE. 3) CALL PRT4(MAX,N,XP,X,FP,F)

!  If too many function evaluations occur, terminate the algorithm.
            IF(NFCNEV .GE. MAXEVL) THEN
              CALL PRT5
              IF (.NOT. MAX) FOPT = -FOPT
              IER = 1
              RETURN
            END IF

!  Accept the new point if the function value increases.
            IF(FP .Gt. F) THEN
              IF(IPRINT .GE. 3) THEN
                WRITE(*,'(''  POINT ACCEPTED'')')
              END IF
              DO I = 1, N
                X(I) = XP(I)
              enddo
              F = FP
              NACC = NACC + 1
              NACP(H) = NACP(H) + 1
              NUP = NUP + 1

!  If greater than any other point, record as new optimum.
              IF (FP .GT. FOPT) THEN
                IF(IPRINT .GE. 3) THEN
                  WRITE(*,'(''  NEW OPTIMUM'')')
                END IF
                DO I = 1, N
                  XOPT(I) = XP(I)
                enddo
                FOPT = FP

                NNEW = NNEW + 1
              END IF

!  If the point is lower, use the Metropolis criteria to decide on
!  acceptance or rejection.
            ELSE
              P = EXPREP((FP - F)/T)
              PP = RANMAR()
              IF (PP .LT. P) THEN
                IF(IPRINT .GE. 3) CALL PRT6(MAX)
                DO I = 1, N
                  X(I) = XP(I)
                enddo
                F = FP
                NACC = NACC + 1
                NACP(H) = NACP(H) + 1
                NDOWN = NDOWN + 1
              ELSE
                NREJ = NREJ + 1
                IF(IPRINT .GE. 3) CALL PRT7(MAX)
              END IF
            END IF

          enddo
        enddo

!  Adjust VM so that approximately half of all evaluations are accepted.
        DO I = 1, N
          RATIO = DBLE(NACP(I)) /DBLE(NS)
          IF (RATIO .GT. .6) THEN
            VM(I) = VM(I)*(1. + C(I)*(RATIO - .6)/.4)
          ELSE IF (RATIO .LT. .4) THEN
            VM(I) = VM(I)/(1. + C(I)*((.4 - RATIO)/.4))
          END IF
          IF (VM(I) .GT. (UB(I)-LB(I))) THEN
            VM(I) = UB(I) - LB(I)
          END IF
        enddo

        IF(IPRINT .GE. 2) THEN
          CALL PRT8(N,VM,XOPT,X)
        END IF

        DO I = 1, N
          NACP(I) = 0
        enddo

      enddo

      IF(IPRINT .GE. 1) THEN
        CALL PRT9(MAX,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
      END IF

!  Check termination criteria.
      QUIT = .FALSE.
      FSTAR(1) = F
      IF ((FOPT - FSTAR(1)) .LE. EPS) QUIT = .TRUE.
      DO I = 1, NEPS
        IF (ABS(F - FSTAR(I)) .GT. EPS) QUIT = .FALSE.
      enddo

!  Terminate SA if appropriate.
      IF (QUIT) THEN
        DO I = 1, N
          X(I) = XOPT(I)
        enddo
        IER = 0
        IF (.NOT. MAX) FOPT = -FOPT
        IF(IPRINT .GE. 1) CALL PRT10
        RETURN
      END IF

!  If termination criteria is not met, prepare for another loop.
      T = RT*T
      DO I = NEPS, 2, -1
        FSTAR(I) = FSTAR(I-1)
      enddo
      F = FOPT
      DO I = 1, N
        X(I) = XOPT(I)
      enddo

!  Loop again.
      GO TO 100

      END

      FUNCTION  EXPREP(RDUM)
      implicit none
!  This function replaces exp to avoid under- and overflows and is
!  designed for IBM 370 type machines. It may be necessary to modify
!  it for other machines. Note that the maximum and minimum values of
!  EXPREP are such that they has no effect on the algorithm.

      DOUBLE PRECISION  RDUM, EXPREP

      IF (RDUM .GT. 174.) THEN
        EXPREP = 3.69D+75
      ELSE IF (RDUM .LT. -180.) THEN
        EXPREP = 0.0
      ELSE
        EXPREP = EXP(RDUM)
      END IF

      RETURN
      END

      subroutine RMARIN(IJ,KL)
      implicit none
!  This subroutine and the next function generate random numbers. See
!  the comments for SA for more information. The only changes from the
!  orginal code is that (1) the test to make sure that RMARIN runs first
!  was taken out since SA assures that this is done (this test didn't
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was
!  taken out since ivec isn't used. With these exceptions, all following
!  lines are original.

! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
      real s,t
      real U(97), C, CD, CM
      integer IJ,KL,i,j,k,l,ii,jj,m
      integer I97, J97
      common /raset1/ U, C, CD, CM, I97, J97
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.                          &
     &KL .lt. 0  .or.  KL .gt. 30081 ) then
        print '(A)',' The first random number seed must have a value'
        print '(A)',' between 0 and 31328'
        print '(A)',' The second seed must have a value between 0 and'
        print '(A)',' 30081'
        stop
      endif
      i = mod(IJ/177, 177) + 2
      j = mod(IJ    , 177) + 2
      k = mod(KL/169, 178) + 1
      l = mod(KL,     169)
      do ii = 1, 97
        s = 0.0
        t = 0.5
        do jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          if (mod(l*m, 64) .ge. 32) then
            s = s + t
          endif
          t = 0.5 * t
        enddo
        U(ii) = s
      enddo
      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0
      I97 = 97
      J97 = 33
      return
      end

      real function ranmar()
      implicit none
      real U(97), C, CD, CM, uni
      integer I97, J97
      common /raset1/ U, C, CD, CM, I97, J97
      uni = U(I97) - U(J97)
      if( uni .lt. 0.0 ) uni = uni + 1.0
      U(I97) = uni
      I97 = I97 - 1
      if(I97 .eq. 0) I97 = 97
      J97 = J97 - 1
      if(J97 .eq. 0) J97 = 97
      C = C - CD
      if( C .lt. 0.0 ) C = C + CM
      uni = uni - C
      if( uni .lt. 0.0 ) uni = uni + 1.0
      RANMAR = uni
      return
      END

      SUBROUTINE PRT1
      implicit none
!  This subroutine prints intermediate output, as does PRT2 through
!  PRT10. Note that if SA is minimizing the function, the sign of the
!  function value and the directions (up/down) are reversed in all
!  output to correspond with the actual function optimization. This
!  correction is because SA was written to maximize functions and
!  it minimizes by maximizing the negative a function.

      WRITE(*,'(/,''THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS '')')
      WRITE(*,'(/,''(LB AND UB). EXECUTION TERMINATED WITHOUT ANY'')')
      WRITE(*,'(/,''OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  '')')
      WRITE(*,'(/,''LB(I) .LT. X(I) .LT. UB(I), I = 1, N. ''/)')

      RETURN
      END

      SUBROUTINE PRT2(MAX,N,X,F)
      implicit none

      DOUBLE PRECISION  X(*), F
      INTEGER  N
      LOGICAL  MAX

      WRITE(*,'(''  '')')
      CALL PRTVEC(X,N,'INITIAL X')
      IF (MAX) THEN
        WRITE(*,'(''  INITIAL F: '',/, G25.18)') F
      ELSE
        WRITE(*,'(''  INITIAL F: '',/, G25.18)') -F
      END IF

      RETURN
      END

      SUBROUTINE PRT3(MAX,N,XP,X,FP,F)
      implicit none

      DOUBLE PRECISION  XP(*), X(*), FP, F
      INTEGER  N
      LOGICAL  MAX

      WRITE(*,'(''  '')')
      CALL PRTVEC(X,N,'CURRENT X')
      IF (MAX) THEN
        WRITE(*,'(''  CURRENT F: '',G25.18)') F
      ELSE
        WRITE(*,'(''  CURRENT F: '',G25.18)') -F
      END IF
      CALL PRTVEC(XP,N,'TRIAL X')
      WRITE(*,'(''  POINT REJECTED SINCE OUT OF BOUNDS'')')

      RETURN
      END

      SUBROUTINE PRT4(MAX,N,XP,X,FP,F)
      implicit none

      DOUBLE PRECISION  XP(*), X(*), FP, F
      INTEGER  N
      LOGICAL  MAX

      WRITE(*,'(''  '')')
      CALL PRTVEC(X,N,'CURRENT X')
      IF (MAX) THEN
        WRITE(*,'(''  CURRENT F: '',G25.18)') F
        CALL PRTVEC(XP,N,'TRIAL X')
        WRITE(*,'(''  RESULTING F: '',G25.18)') FP
      ELSE
        WRITE(*,'(''  CURRENT F: '',G25.18)') -F
        CALL PRTVEC(XP,N,'TRIAL X')
        WRITE(*,'(''  RESULTING F: '',G25.18)') -FP
      END IF

      RETURN
      END

      SUBROUTINE PRT5
      implicit none

      WRITE(*, '(/,''TOO MANY FUNCTION EVALUATIONS; CONSIDER '',/)')
      WRITE(*, '(/,''INCREASING MAXEVL OR EPS, OR DECREASING '',/)')
      WRITE(*, '(/,''NT OR RT. THESE RESULTS ARE LIKELY TO BE '',/)')
      WRITE(*, '(/,''POOR.'',/)')

      RETURN
      END

      SUBROUTINE PRT6(MAX)
      implicit none

      LOGICAL  MAX

      IF (MAX) THEN
        WRITE(*,'(''  THOUGH LOWER, POINT ACCEPTED'')')
      ELSE
        WRITE(*,'(''  THOUGH HIGHER, POINT ACCEPTED'')')
      END IF

      RETURN
      END

      SUBROUTINE PRT7(MAX)
      implicit none

      LOGICAL  MAX

      IF (MAX) THEN
        WRITE(*,'(''  LOWER POINT REJECTED'')')
      ELSE
        WRITE(*,'(''  HIGHER POINT REJECTED'')')
      END IF

      RETURN
      END

      SUBROUTINE PRT8(N,VM,XOPT,X)
      implicit none

      DOUBLE PRECISION  VM(*), XOPT(*), X(*)
      INTEGER  N

      WRITE(*,'(/''INTERMEDIATE RESULTS AFTER'',/)')
      WRITE(*,'(/'' STEP LENGTH ADJUSTMENT'',/)')
      CALL PRTVEC(VM,N,'NEW STEP LENGTH (VM)')
      CALL PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
      CALL PRTVEC(X,N,'CURRENT X')
      WRITE(*,'('' '')')

      RETURN
      END

      SUBROUTINE PRT9(MAX,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
      implicit none

      DOUBLE PRECISION  XOPT(*), VM(*), T, FOPT
      INTEGER  N, NUP, NDOWN, NREJ, LNOBDS, NNEW, TOTMOV
      LOGICAL  MAX

      TOTMOV = NUP + NDOWN + NREJ

      WRITE(*,'(/,''INTERMEDIATE RESULTS BEFORE'',/)')
      WRITE(*,'(/,  ''  NEXT TEMPERATURE REDUCTION'',/)')
      WRITE(*,'(''  CURRENT TEMPERATURE:            '',G12.5)') T
      IF (MAX) THEN
        WRITE(*,'(''  MAX FUNCTION VALUE SO FAR:  '',G25.18)') FOPT
        WRITE(*,'(''  TOTAL MOVES:                '',I8)') TOTMOV
        WRITE(*,'(''     UPHILL:                  '',I8)') NUP
        WRITE(*,'(''     ACCEPTED DOWNHILL:       '',I8)') NDOWN
        WRITE(*,'(''     REJECTED DOWNHILL:       '',I8)') NREJ
        WRITE(*,'(''  OUT OF BOUNDS TRIALS:       '',I8)') LNOBDS
        WRITE(*,'(''  NEW MAXIMA THIS TEMPERATURE:'',I8)') NNEW
      ELSE
        WRITE(*,'(''  MIN FUNCTION VALUE SO FAR:  '',D25.18)') -FOPT
        WRITE(*,'(''  TOTAL MOVES:                '',I8)') TOTMOV
        WRITE(*,'(''     DOWNHILL:                '',I8)')  NUP
        WRITE(*,'(''     ACCEPTED UPHILL:         '',I8)')  NDOWN
        WRITE(*,'(''     REJECTED UPHILL:         '',I8)')  NREJ
        WRITE(*,'(''  TRIALS OUT OF BOUNDS:       '',I8)')  LNOBDS
        WRITE(*,'(''  NEW MINIMA THIS TEMPERATURE:'',I8)')  NNEW
      END IF
      CALL PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
      CALL PRTVEC(VM,N,'STEP LENGTH (VM)')
      WRITE(*,'('' '')')

      RETURN
      END

      SUBROUTINE PRT10
      implicit none

      WRITE(*,'(/,''  SA ACHIEVED TERMINATION CRITERIA. IER = 0. '',/)')

      RETURN
      END

      SUBROUTINE PRTVEC(VECTOR,NCOLS,NAME)
      implicit none
!  This subroutine prints the double precision vector named VECTOR.
!  Elements 1 thru NCOLS will be printed. NAME is a character variable
!  that describes VECTOR. Note that if NAME is given in the call to
!  PRTVEC, it must be enclosed in quotes. If there are more than 10
!  elements in VECTOR, 10 elements will be printed on each line.

      integer i,j,ll,lines
      INTEGER NCOLS
      DOUBLE PRECISION VECTOR(NCOLS)
      CHARACTER *(*) NAME

      WRITE(*,1001) NAME

      IF (NCOLS .GT. 10) THEN
        LINES = INT(NCOLS/10.)

        DO I = 1, LINES
          LL = 10*(I - 1)
          WRITE(*,1000) (VECTOR(J),J = 1+LL, 10+LL)
        enddo

        WRITE(*,1000) (VECTOR(J),J = 11+LL, NCOLS)
      ELSE
        WRITE(*,1000) (VECTOR(J),J = 1, NCOLS)
      END IF

 1000 FORMAT( 10(G12.5,1X))
 1001 FORMAT(/,25X,A)

      RETURN
      END

