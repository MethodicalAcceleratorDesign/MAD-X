SUBROUTINE twiss(rt,disp0,tab_name,sector_tab_name)
  use twiss0fi
  use twissafi
  use twisslfi
  use twisscfi
  use twissotmfi
  use twissbeamfi
  use trackfi, only : fsecarb
  use fasterror
  use matrices, only : EYE
  use math_constfi, only : zero, one, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TWISS command: Track linear lattice parameters.                  *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6), disp0(6), cp_thrd=1d-12 ! coupling limit
  integer :: tab_name(*)
  integer :: sector_tab_name(*) ! holds sectormap data

  integer :: i, j, ithr_on
  integer :: chrom, eflag, chrom_warn
  double precision :: orbit0(6), orbit(6), tt(6,6,6), ddisp0(6), r0mat(2,2)
  double precision :: s0mat(6,6) ! initial sigma matrix
  character(len=48) :: charconv
  character(len=150) :: warnstr
  logical :: fast_error_func

  double precision, external :: get_value
  integer, external :: get_option

  integer, parameter :: izero=0, ione=1
  !---- Initialization
  table_name = charconv(tab_name)
  sectorTableName = charconv(sector_tab_name)
  chrom=0
  chrom_warn=1
  eflag=0
  ithr_on=0
  fsecarb=.false.

  ORBIT0 = zero
  call get_node_vector('orbit0 ', 6, orbit0)
  RT = EYE
  RW = EYE
  TT = zero
  call get_disp0(disp0)
  DDISP0 = zero
  R0MAT = zero
  S0MAT = zero
  OPT_FUN0 = zero
  OPT_FUN = zero
  DISP = zero
  DDISP = zero
  RMAT = zero
  SIGMAT = zero
  betx=zero; alfx=zero; amux=zero; bxmax=zero; dxmax=zero; xcomax=zero
  sigxco=zero; sigdx=zero; cosmux=zero; wx=zero; phix=zero; dmux=zero
  qx=zero; sinmux=zero; xix=zero

  bety=zero; alfy=zero; amuy=zero; bymax=zero; dymax=zero; ycomax=zero
  sigyco=zero; sigdy=zero; cosmuy=zero; wy=zero; phiy=zero; dmuy=zero
  qy=zero; sinmuy=zero; xiy=zero

  gammacp=one
  nmode_flip = 0
  mode_flip  =.false.

  synch_1=zero;  synch_2=zero;  synch_3=zero;  synch_4=zero;  synch_5=zero
  synch_6=zero;  synch_8=zero

  suml=zero; circ=zero; eta=zero; alfa=zero; gamtr=zero; wgt=zero

  !---- Track chromatic functions
  chrom = get_option('twiss_chrom ')

  !---- flag if called from match process
  !---- get match flag for storing variables in nodes
  match_is_on = get_option('match_is_on ') .ne. 0
  if (match_is_on) chrom_warn = get_option('chrom_match ')

  !---- flags for writing cumulative or lumped matrices
  rmatrix = get_value('twiss ','rmatrix ').ne.zero
  sectormap = get_option('twiss_sector ').ne.zero

  !---- Get circumference
  circ   = get_value('probe ','circ ')
  if (circ .eq. zero) call fort_fail('TWISS: ', 'Zero length sequence.')

  !---- Get beam parameters
  radiate = get_value('probe ','radiate ') .ne. zero
  energy  = get_value('probe ','energy ')
  deltap  = get_value('probe ','deltap ')
  beta    = get_value('probe ','beta ')
  gamma   = get_value('probe ','gamma ')
  pc      = get_value('probe ','pc ')
  arad    = get_value('probe ','arad ')
  dtbyds  = get_value('probe ','dtbyds ')
  charge  = get_value('probe ','charge ')
  npart   = get_value('probe ','npart ')

  !---- Set fast_error_func flag to use faster error function
  !---- including tables. Thanks to late G. Erskine
  fast_error_func = get_option('fast_error_func ') .ne. 0
  if (fast_error_func .and. .not.fasterror_on) then
     call wzset
     fasterror_on = .true.
  endif

  if (get_option('twiss_inval ') .ne. 0) then
     !---- Initial values from command attributes.
     call twinifun(opt_fun0, rt)
     if (get_option('twiss_print ') .ne. 0)  then
        print *, ' '
        print '(''open line - error with deltap: '',1p,e14.6)', deltap
        print '(''initial orbit vector: '', 1p,6e14.6)', orbit0
     endif
     call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,0,0,ithr_on); if (eflag.ne.0) go to 900
     if (get_option('twiss_print ') .ne. 0) &
        print '(''final orbit vector:   '', 1p,6e14.6)', orbit
     call tmsigma(s0mat) !-- from initial conditions in opt_fun0
  else
     !---- Initial values from periodic solution.
     call tmclor(orbit0,.true.,.true.,opt_fun0,rt,tt,eflag);          if (eflag.ne.0) go to 900
!    LD: useless, identical to last call to tmfrst in tmclor...
!    call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,0,0,ithr_on); if (eflag.ne.0) go to 900
     call twcpin(rt,disp0,r0mat,eflag);                               if (eflag.ne.0) go to 900
     !----IrinaTecker: For Sagan-Rubin and Talman versions of coupling
!    call twcpin_sagan(rt,disp0,r0mat,eflag);                         if (eflag.ne.0) go to 900
!    call twcpin_talman(rt,disp0,r0mat,eflag);                        if (eflag.ne.0) go to 900
     !---- Initialize opt_fun0
     call twinifun(opt_fun0,rt)
     call tmsigma_emit(rt,s0mat)
  endif

  if (sectormap)  then
     SORB = ORBIT0
     SRMAT = EYE
     STMAT = zero
  endif

  ! save sigma matrix
  do i= 1,6
  do j= 1,6
    opt_fun0(74 + (i-1)*6 + j) = s0mat(i,j)
  enddo
  enddo
  sigmat = s0mat

  !---- Build table of lattice functions, coupled.
  call twcpgo(rt,orbit0)


  !---- List chromatic functions.
  if (chrom .ne. 0) then
    if ( any(rt(1:2,3:4) .gt. cp_thrd)  .or.  any(rt(3:4,1:2) .gt. cp_thrd) ) then
       if (chrom_warn .eq. 1) then
          write (warnstr, '(a)') 'Calculation of Wx, Wy etc. could be inaccurate due to coupling!'
          call fort_warn('TWISS: ', warnstr)
       endif
    endif
     call twbtin(rt,tt)
     call twchgo
  endif

  !---- Print summary
  call double_to_table_curr('summ ','nflips ' , dble(nmode_flip))

  if (MOD(nmode_flip, 2).ne.0) then
     write (warnstr, '(a,i5)') 'Total number of modes flips is not even! Nflips = ', nmode_flip
     print*, 'Total number of modes flips is not even! Nflips = ', nmode_flip
  endif

  if (get_option('twiss_summ ') .ne. 0) call tw_summ(rt,tt)

  if (get_option('keeporbit ') .ne. 0)  then
     call store_node_vector('orbit0 ', 6, opt_fun0(9))
  endif

  call set_option('twiss_success ', ione)
  return

900 call set_option('twiss_success ', izero)

end SUBROUTINE twiss

SUBROUTINE tmrefe(rt)
  use twissbeamfi
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Transfer matrix w.r.t. ideal orbit for one period.               *
  !     Ignores cavities, radiation, and imperfections.                  *
  !     entry point for mad_twiss, mad_emit and mad_beam                 *
  !     Output:                                                          *
  !     rt(6,6) (double) transfer matrix.                                *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6)

  integer :: eflag, ithr_on
  double precision :: orbit0(6), orbit(6), tt(6,6,6)

  double precision, external :: get_value

  !---- Get beam parameters
  radiate  = get_value('probe ','radiate ') .ne. zero
  energy = get_value('probe ','energy ')
  deltap = get_value('probe ','deltap ')
  beta   = get_value('probe ','beta ')
  gamma  = get_value('probe ','gamma ')
  pc     = get_value('probe ','pc ')
  arad   = get_value('probe ','arad ')
  dtbyds = get_value('probe ','dtbyds ')
  charge = get_value('probe ','charge ')
  npart  = get_value('probe ','npart ')

  ithr_on = 0
  ORBIT0 = zero ; ORBIT = zero ; TT = zero

  !---- Get transfer matrix.
  call tmfrst(orbit0,orbit,.false.,.false.,rt,tt,eflag,0,0,ithr_on)
end SUBROUTINE tmrefe

SUBROUTINE tmrefo(kobs,orbit0,orbit,rt)
  use twiss0fi
  use twissbeamfi
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Transfer matrix w.r.t. ideal orbit for one period.               *
  !     entry point for mad_dynap and mad_track                          *
  !     Input:                                                           *
  !     kobs    if > 0, track until node with this obs. point number     *
  !     Output:                                                          *
  !     orbit0(6) (double) closed orbit at start=end                     *
  !     orbit(6)  (double) closed orbit at obs. point kobs, or at end    *
  !     rt(6,6) (double) transfer matrix.                                *
  !----------------------------------------------------------------------*
  integer :: kobs
  double precision :: orbit0(6), orbit(6), rt(6,6)

  integer :: eflag, ithr_on
  double precision :: opt_fun0(fundim), tt(6,6,6)

  double precision, external :: get_value

  integer, parameter :: izero=0, ione=1

  !---- Get beam parameters
  radiate  = get_value('probe ','radiate ') .ne. zero
  energy = get_value('probe ','energy ')
  deltap = get_value('probe ','deltap ')
  beta   = get_value('probe ','beta ')
  gamma  = get_value('probe ','gamma ')
  pc     = get_value('probe ','pc ')
  arad   = get_value('probe ','arad ')
  dtbyds = get_value('probe ','dtbyds ')
  charge = get_value('probe ','charge ')
  npart  = get_value('probe ','npart ')

  ithr_on = izero
  ORBIT0 = zero
  !---- Get closed orbit and coupled transfer matrix.
  call tmclor(orbit0,.true.,.true.,opt_fun0,rt,tt,eflag)
  call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,kobs,0,ithr_on) ! useless?
end SUBROUTINE tmrefo

SUBROUTINE twinifun(opt_fun0,rt)
  use twiss0fi
  use twisslfi
  use twissbeamfi, only : energy
  use math_constfi, only : zero, twopi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial twiss parameters put to opt_fun0.                        *
  !     Input/output:                                                    *
  !     opt_fun0(fundim) (double) initial optical values:                *
  !     betx0,alfx0,amux0,bety0,alfy0,amuy0, etc.                        *
  !----------------------------------------------------------------------*
  double precision :: opt_fun0(*), rt(6,6)

  integer :: i, j
  double precision :: betx, alfx, mux, bety, alfy, muy, dx, dpx, dy, dpy
  double precision :: x, px, y, py, t, pt
  double precision :: wx, phix, dmux, wy, phiy, dmuy, ddx, ddpx, ddy, ddpy
  double precision :: r(2,2)

  double precision, external :: get_value

  !---- Initialize
  betx = get_value('twiss ','betx ')
  bety = get_value('twiss ','bety ')

  if (betx.gt.zero) opt_fun0(3) = betx
  if (bety.gt.zero) opt_fun0(6) = bety

  if (opt_fun0(3).le.zero .or. opt_fun0(6).le.zero) &
       call fort_fail('TWINIFUN: ', 'BETX and BETY must be both larger than zero.')

  alfx=  get_value('twiss ','alfx ')
  mux=   get_value('twiss ','mux ')
  alfy=  get_value('twiss ','alfy ')
  muy=   get_value('twiss ','muy ')
  x=     get_value('twiss ','x ')
  px=    get_value('twiss ','px ')
  y=     get_value('twiss ','y ')
  py=    get_value('twiss ','py ')
  t=     get_value('twiss ','t ')
  pt=    get_value('twiss ','pt ')
  dx=    get_value('twiss ','dx ')
  dpx=   get_value('twiss ','dpx ')
  dy=    get_value('twiss ','dy ')
  dpy=   get_value('twiss ','dpy ')
  wx=    get_value('twiss ','wx ')
  phix=  get_value('twiss ','phix ')
  dmux=  get_value('twiss ','dmux ')
  wy=    get_value('twiss ','wy ')
  phiy=  get_value('twiss ','phiy ')
  dmuy=  get_value('twiss ','dmuy ')
  ddx=   get_value('twiss ','ddx ')
  ddpx=  get_value('twiss ','ddpx ')
  ddy=   get_value('twiss ','ddy ')
  ddpy=  get_value('twiss ','ddpy ')
  r(1,1) = get_value('twiss ','r11 ')
  r(1,2) = get_value('twiss ','r12 ')
  r(2,1) = get_value('twiss ','r21 ')
  r(2,2) = get_value('twiss ','r22 ')

  if (alfx  .ne.zero) opt_fun0(4 ) = alfx
  if (mux   .ne.zero) opt_fun0(5 ) = mux * twopi
  if (alfy  .ne.zero) opt_fun0(7 ) = alfy
  if (muy   .ne.zero) opt_fun0(8 ) = muy * twopi
  if (x     .ne.zero) opt_fun0(9 ) = x
  if (px    .ne.zero) opt_fun0(10) = px
  if (y     .ne.zero) opt_fun0(11) = y
  if (py    .ne.zero) opt_fun0(12) = py
  if (t     .ne.zero) opt_fun0(13) = t
  if (pt    .ne.zero) opt_fun0(14) = pt
  if (dx    .ne.zero) opt_fun0(15) = dx
  if (dpx   .ne.zero) opt_fun0(16) = dpx
  if (dy    .ne.zero) opt_fun0(17) = dy
  if (dpy   .ne.zero) opt_fun0(18) = dpy
  if (wx    .ne.zero) opt_fun0(19) = wx
  if (phix  .ne.zero) opt_fun0(20) = phix * twopi
  if (dmux  .ne.zero) opt_fun0(21) = dmux * twopi
  if (wy    .ne.zero) opt_fun0(22) = wy
  if (phiy  .ne.zero) opt_fun0(23) = phiy * twopi
  if (dmuy  .ne.zero) opt_fun0(24) = dmuy * twopi
  if (ddx   .ne.zero) opt_fun0(25) = ddx
  if (ddpx  .ne.zero) opt_fun0(26) = ddpx
  if (ddy   .ne.zero) opt_fun0(27) = ddy
  if (ddpy  .ne.zero) opt_fun0(28) = ddpy
  if (r(1,1).ne.zero) opt_fun0(29) = r(1,1)
  if (r(1,2).ne.zero) opt_fun0(30) = r(1,2)
  if (r(2,1).ne.zero) opt_fun0(31) = r(2,1)
  if (r(2,2).ne.zero) opt_fun0(32) = r(2,2)
  if (energy.ne.zero) opt_fun0(33) = energy
  if (rmatrix) then
     do i= 1,6
        do j= 1,6
           opt_fun0(33 + (i-1)*6 + j) = rt(i,j)
        enddo
     enddo
  endif

end SUBROUTINE twinifun

SUBROUTINE twprep(save,case,opt_fun,position,ii)
  use twissafi
  use twisslfi
  use twissotmfi
  use math_constfi, only : zero, twopi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Finalize twiss parameters, if save flag fill                     *
  !     table with twiss parameters.                                     *
  !     Input:                                                           *
  !     case        (integer) =1 fill from twcpgo; =2 fill from twchgo   *
  !     position    (double)  end position of element                    *
  !     ii          (integer) index of interpolated segment
  !     Input/output:                                                    *
  !     opt_fun(fundim) (double) optical values:                         *
  !     betx,alfx,amux,bety,alfy,amuy, etc.                              *
  !----------------------------------------------------------------------*
  integer :: save, case, ii
  double precision :: opt_fun(*), position

  integer :: i
  double precision :: opt5, opt8, opt20, opt21, opt23, opt24

  if (case .eq. 1) then
     !--- fill with data from twcpgo (Twiss Couple)
     opt_fun(2) = position
     opt5 = opt_fun(5) ; opt_fun(5) = opt_fun(5) / twopi
     opt8 = opt_fun(8) ; opt_fun(8) = opt_fun(8) / twopi
     if (save .ne. 0) call twfill(case,opt_fun,position)
     if (match_is_on) call copy_twiss_data(opt_fun, 0, 110, ii)
     opt_fun(5) = opt5
     opt_fun(8) = opt8

  elseif (case.eq.2) then
     !--- fill with data from twchgo (Twiss Chrom)
     opt20 = opt_fun(20) ; opt_fun(20) = opt_fun(20) / twopi
     opt21 = opt_fun(21) ; opt_fun(21) = opt_fun(21) / twopi
     opt23 = opt_fun(23) ; opt_fun(23) = opt_fun(23) / twopi
     opt24 = opt_fun(24) ; opt_fun(24) = opt_fun(24) / twopi
     if (save .ne. 0) call twfill(case,opt_fun,position)
     if (match_is_on) call copy_twiss_data(opt_fun, 18, 10, ii)! fill 10 values starting with wx
     opt_fun(20) = opt20
     opt_fun(21) = opt21
     opt_fun(23) = opt23
     opt_fun(24) = opt24
  endif
end SUBROUTINE twprep

SUBROUTINE twfill(case,opt_fun,position)
  use twissafi
  use twisslfi
  use twissotmfi
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Fill twiss table with twiss parameters.                          *
  !     Input:                                                           *
  !     case        (integer) =1 fill from twcpgo; =2 fill from twchgo   *
  !     position    (double)  end position of element                    *
  !     Input/output:                                                    *
  !     opt_fun(fundim) (double) optical values:                         *
  !     betx,alfx,amux,bety,alfy,amuy, etc.                              *
  !----------------------------------------------------------------------*
  integer :: case
  double precision :: opt_fun(*), position

  double precision, external :: get_value

  ripken = get_value('twiss ','ripken ') .ne. zero

  if (case .eq. 1) then
     call vector_to_table_curr(table_name, 's '    , opt_fun( 2), 17) ! fill 17 values starting with s
     call vector_to_table_curr(table_name, 'r11 '  , opt_fun(29), 5 ) ! fill 5  values starting with r11
     call vector_to_table_curr(table_name, 'sig11 ', opt_fun(75), 36) ! fill 36 values starting with sigmat11
     call vector_to_table_curr(table_name, 'kmax ' , opt_fun(70), 5 ) ! fill 5  values starting with kmax
     if (rmatrix) call vector_to_table_curr(table_name, 're11 ', opt_fun(34), 36) ! fill Rmatrix
     if (ripken)  call twfill_ripken(opt_fun)
  elseif (case .eq. 2) then
     call vector_to_table_curr(table_name, 'wx ', opt_fun(19), 10) ! fill 10 values starting with wx
  endif

  !---- Augment table twiss
  call augment_count(table_name)

end SUBROUTINE twfill

SUBROUTINE twfill_ripken(opt_fun)
  use math_constfi, only : zero, one, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Fill twiss table with Ripken-Mais twiss parameters.              *
  !     beta11, beta12, beta21, beta22,                                  *
  !     alfa11, alfa12, alfa21, alfa22,                                  *
  !     gama11, gama12, gama21, gama22,                                  *
  !----------------------------------------------------------------------*
  double precision :: opt_fun(*)

  double precision :: kappa, betx, bety, alfx, alfy, gamx, gamy, r11, r12, r21, r22
  double precision :: beta11, beta12, beta21, beta22, &
                      alfa11, alfa12, alfa21, alfa22, &
                      gama11, gama12, gama21, gama22

  betx = opt_fun(3);   bety = opt_fun(6);   alfx = opt_fun(4);   alfy = opt_fun(7)
  r11 =  opt_fun(29);  r12 =  opt_fun(30);  r21 =  opt_fun(31);  r22 =  opt_fun(32)

  kappa = one/(one + (r11*r22 - r12*r21))

  gamx = (one + alfx**2) / betx;  gamy = (one + alfy**2) / bety

  beta11 = kappa * betx
  beta22 = kappa * bety
  beta12 = kappa * ( r22**2*bety + two*r12*r22*alfy + r12**2*gamy )
  beta21 = kappa * ( r11**2*betx - two*r12*r11*alfx + r12**2*gamx )

  alfa11 =  kappa * alfx
  alfa22 =  kappa * alfy
  alfa12 =  kappa * ( r21*r22*bety + (r12*r21 + r11*r22)*alfy + r11*r12*gamy )
  alfa21 = -kappa * ( r21*r11*betx - (r12*r21 + r11*r22)*alfx + r12*r22*gamx )

  gama11 = kappa * gamx
  gama22 = kappa * gamy
  gama12 = zero;  if (beta12 .ne. zero) gama12 = ((one-kappa)**2 + alfa12**2) / beta12
  gama21 = zero;  if (beta21 .ne. zero) gama21 = ((one-kappa)**2 + alfa21**2) / beta21

  call double_to_table_curr('twiss ','beta11 ' ,beta11)
  call double_to_table_curr('twiss ','beta12 ' ,beta12)
  call double_to_table_curr('twiss ','beta21 ' ,beta21)
  call double_to_table_curr('twiss ','beta22 ' ,beta22)
  call double_to_table_curr('twiss ','alfa11 ' ,alfa11)
  call double_to_table_curr('twiss ','alfa12 ' ,alfa12)
  call double_to_table_curr('twiss ','alfa21 ' ,alfa21)
  call double_to_table_curr('twiss ','alfa22 ' ,alfa22)
  call double_to_table_curr('twiss ','gama11 ' ,gama11)
  call double_to_table_curr('twiss ','gama12 ' ,gama12)
  call double_to_table_curr('twiss ','gama21 ' ,gama21)
  call double_to_table_curr('twiss ','gama22 ' ,gama22)

end SUBROUTINE twfill_ripken

SUBROUTINE tmclor(guess,fsec,ftrk,opt_fun0,rt,tt,eflag)
  use twissbeamfi, only : deltap
  use matrices, only : EYE
  use math_constfi, only : zero, one
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Find closed orbit for a beam line sequence.                      *
  !     Called from emit as well
  !     Input:                                                           *
  !     guess(6)     (double)  first guess for orbit start               *
  !     fsec         (logical) if true, return second order terms.       *
  !     ftrk         (logical) if true, track orbit.                     *
  !     Output:                                                          *
  !     guess(6)     (double)  final orbit                               *
  !     opt_fun0(fundim) (double)  initial optical values:               *
  !     betx0,alfx0,amux0,bety0,alfy0,amuy0, etc.                        *
  !     rt(6,6)      (double)  transfer matrix.                          *
  !     tt(6,6,6)    (double)  second order terms.                       *
  !     eflag        (integer) error flag (0: OK, else != 0)             *
  !----------------------------------------------------------------------*
  double precision :: guess(6), opt_fun0(*), rt(6,6), tt(6,6,6)
  logical :: fsec, ftrk
  integer :: eflag

  logical :: pflag
  integer :: i, k, irank, itra, thr_on, ithr_on, save_opt
  double precision :: cotol, err
  double precision :: orbit0(6), orbit(6), a(6,7), b(4,5)

  integer, external :: get_option
  double precision, external :: get_variable, get_value
  logical, external :: m66sta
  integer, parameter :: itmax=20

  deltap = get_value('probe ','deltap ')

  !---- Initialize.
  ithr_on = 0
  thr_on = get_option('threader ')
  pflag = get_option('twiss_print ') .ne. 0
  cotol = get_variable('twiss_tol ')
  eflag = 0

  !---- Initialize guess.
  ORBIT0 = GUESS

  !---- Iteration for closed orbit.
  iterate: do itra = 1, itmax

     !---- Track orbit and transfer matrix.
     call tmfrst(orbit0,orbit,fsec,ftrk,rt,tt,eflag,0,0,thr_on)

     if (eflag.ne.0) return
     ! turn off threader immediately after first iteration, ie after first turn
     thr_on = 0

     if (.not.m66sta(rt)) then
        !---- Solve for dynamic case.
        err = zero
        A(:,:6) = RT(:,:6) - EYE
        A(1:6,7) = ORBIT(1:6) - ORBIT0(1:6)
        err = maxval(abs(A(1:6,7)))

        call solver(a,6,1,irank)
        if (irank.lt.6) then
           print *, 'Singular matrix occurred during closed orbit search.'
           eflag = 1
           return
        endif
        ORBIT0(1:6) = ORBIT0(1:6) - A(1:6,7)

     else
        !---- Solve for static case.
        err = zero
        B(:4,:4) = RT(:4,:4) - EYE(:4,:4)
        B(1:4,5) = ORBIT(1:4) - ORBIT0(1:4)
        err = maxval(abs(B(1:4,5)))

        call solver(b,4,1,irank)
        if (irank.lt.4) then
           print *, 'Singular matrix occurred during closed orbit search.'
           eflag = 1
           return
        endif
        ORBIT0(1:4) = ORBIT0(1:4) - B(1:4,5)

     endif

     !---- Message and convergence test.
     if (pflag)  then
        print *, ' '
        print '(''iteration: '',i3,'' error: '',1p,e14.6,'' deltap: '',1p,e14.6)',itra,err,deltap
        print '(''orbit: '', 1p,6e14.6)', orbit0
     endif

     if (err.lt.cotol) then
        save_opt = get_option('keeporbit ')
        call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,0,save_opt,ithr_on)
        OPT_FUN0(9:14) = ORBIT0(1:6)
        GUESS = ORBIT0
        return ! normal exit
     endif

  enddo iterate

  !---- No convergence.
  print '(''Closed orbit did not converge in '', i3, '' iterations'')', itmax
  OPT_FUN0(9:14) = zero
  return

end SUBROUTINE tmclor

SUBROUTINE tmfrst(orbit0,orbit,fsec,ftrk,rt,tt,eflag,kobs,save,thr_on)
  use bbfi
  use twiss0fi
  use twissbeamfi, only : beta, gamma, arad, charge, npart, pc, energy
  use name_lenfi
  use twisscfi
  use spch_bbfi
  use matrices, only : EYE, symp_thrd  ! , symp_thrd_orbit
  use math_constfi, only : zero
  use code_constfi
  use twtapering
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Transfer matrix w.r.t. actual orbit for one (half) superperiod.  *
  !     Misalignment and field errors are considered.                    *
  !     Input:                                                           *
  !     orbit0(6) (double)  first guess for orbit start                  *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Output:                                                          *
  !     orbit(6)  (double)  orbit after one turn starting from orbit0    *
  !     rt(6,6)   (double)  transfer matrix.                             *
  !     tt(6,6,6) (double)  second order terms.                          *
  !     eflag     (integer) error flag.                                  *
  !     kobs      (integer) if > 0, stop at node with this obs. point #  *
  !     save      (integer) if > 0, save orbit at all BPMs               *
  !     thr_on    (integer) if > 0, use threader                         *
  !----------------------------------------------------------------------*
  double precision :: orbit0(6), orbit(6)
  logical :: fsec, ftrk
  double precision :: rt(6,6), tt(6,6,6)
  integer :: eflag, kobs, save, thr_on, i

  logical :: fmap, istaper
  character(len=28) :: tmptxt1, tmptxt2, tmptxt3
  character(len=150) :: warnstr
  character(len=name_len) :: c_name(2), p_name, el_name
  character(len=2) :: ptxt(2)=(/'x-','y-'/)
  integer :: j, code, n_align, nobs, node, old, poc_cnt, debug
  integer :: kpro, corr_pick(2), enable, coc_cnt(2), lastnb, rep_cnt(2)
  double precision :: orbit2(6), ek(6), re(6,6), te(6,6,6), orbitori(6)
  double precision :: al_errors(align_max), dptemp
  double precision :: el, cick, err, nrm, nrm0, anglet, newk0, tempk
  double precision :: parvec(26),  vector(10), reforb(6)
  double precision :: restsum(2), restorb(6,2), restm(6,6,2), restt(6,6,6,2)
  double precision :: cmatr(6,6,2), pmatr(6,6), dorb(6)

  integer, external :: restart_sequ, advance_node, node_al_errors, get_vector, get_option
  double precision, external :: node_value, get_value
  double precision, parameter :: orb_limit=1d1
  integer, parameter :: max_rep=100

  111 continue
  debug = get_option('debug ')

  !---- Initialize
  !---- corr_pick stores for both projection the last pickup used by the
  !     threader in order to avoid corrections in front of it when
  !     restarting.

  CORR_PICK = 0
  REP_CNT = 0
  RESTSUM = zero
  nrm0    = zero
  if (thr_on .gt. 0)  then
     VECTOR = zero
     j = get_vector('threader ', 'vector ', vector)
     if (j .lt. 3) thr_on = 0
  endif
  istaper = get_value('twiss ','tapering ').ne.zero
  !istaper = .true.
  TT = zero
  RT = EYE

  eflag = 0
  suml = zero

  ORBIT = ORBIT0

  parvec(5) = arad
  parvec(6) = charge * npart
  parvec(7) = gamma
  bbd_cnt = 0
  bbd_flag = 1
  i_spch = 0 !!!!!!!!!!

  !---  start
  node = restart_sequ()

  !---  loop over nodes
10 continue

  bbd_pos = node !--- for space charge

  code = node_value('mad8_type ')

!  if (code .eq. code_tkicker)     code = code_kicker ! TKICKER treated as KICKER
  if (code .eq. code_placeholder) code = code_instrument ! PLACEHOLDER treated as INSTRUMENT

  if (thr_on .gt. 0)  then
     !---  threader is on - keep position,orbit,matrix,and tensor for restart
     select case (code)
     case (code_hkicker)
            restsum(1) = suml
            RESTORB(:,1) = ORBIT(:)
            RESTM(:,:,1) = RT(:,:)
            RESTT(:,:,:,1) = TT(:,:,:)

     case (code_kicker)
            restsum(1) = suml
            RESTORB(:,1) = ORBIT(:)
            RESTM(:,:,1) = RT(:,:)
            RESTT(:,:,:,1) = TT(:,:,:)

            restsum(2) = suml
            RESTORB(:,2) = ORBIT(:)
            RESTM(:,:,2) = RT(:,:)
            RESTT(:,:,:,2) = TT(:,:,:)

     case (code_vkicker)
            restsum(2) = suml
            RESTORB(:,2) = ORBIT(:)
            RESTM(:,:,2) = RT(:,:)
            RESTT(:,:,:,2) = TT(:,:,:)

     end select
  endif

  !---- Element length.
  el = node_value('l ')

  nobs = node_value('obs_point ')

  n_align = node_al_errors(al_errors)
  if (n_align .ne. 0)  then
    ORBIT2 = ORBIT
    call tmali1(orbit2,al_errors,beta,gamma,orbit,re)
    RT = matmul(RE,RT)
    !--- IT check if RT is symplectic
    if (thr_on .eq. 0 .and. debug .ne. 0) then
      call m66symp(rt,nrm)
      if ((nrm-nrm0) .gt. symp_thrd) then
        call element_name(el_name,len(el_name))
        write (warnstr,'(a,e13.6,a,a)') "Symplectic deviation: ", nrm, " in element ", el_name
        call fort_warn('THREADER-1: ', warnstr)
      endif
      nrm0 = nrm
    endif
  endif

  !---- Element matrix
  if (istaper) then
    orbitori = orbit

    select case (code)

    case (code_rbend, code_sbend)
      anglet = node_value('angle ')
      do i=1,3
        call tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te,.false.,el)
        dptemp = (orbit(6)+orbitori(6))/(2*beta)
        newk0 = (1+dptemp)*anglet/el
        call store_node_value('k0 ',newk0 )
        orbit = orbitori
  	  enddo
      
    case (code_quadrupole)
      call tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te,.false.,el)
      dptemp = (orbit(6)+orbitori(6))/(2*beta)
      tempk = node_value('k1 ')
      newk0 = tempk/(1-dptemp)-tempk
      call store_node_value('k1tap ',newk0 )

      tempk = node_value('k1s ')
      newk0 = tempk/(1-dptemp)-tempk
      call store_node_value('k1stap ',newk0 )

    case (code_sextupole)
      call tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te,.false.,el)
      dptemp = (orbit(6)+orbitori(6))/(2*beta)
      tempk = node_value('k2 ')
      newk0 = tempk/(1-dptemp) - tempk
      call store_node_value('k2tap ',newk0 )  

      tempk = node_value('k2s ')
      newk0 = tempk/(1-dptemp) - tempk
      call store_node_value('k2stap ',newk0 )

    end select
    orbit = orbitori
  endif

  call tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te,.false.,el)

  !--- if element has a map, concatenate
  if (fmap) then
    call tmcat(.true.,re,te,rt,tt,rt,tt)
    suml = suml + el
    !--- IT check if RT is symplectic
    if (thr_on .eq. zero .and. debug .ne. 0) then
      call m66symp(rt,nrm)
      if ((nrm-nrm0) .gt. symp_thrd) then
        call element_name(el_name,len(el_name))
        write (warnstr,'(a,e13.6,a,a)') "Symplectic deviation: ", nrm, " in element ", el_name
        call fort_warn('THREADER-M: ', warnstr)
      endif
      nrm0 = nrm
    endif
  endif

  if (n_align .ne. 0)  then
    ORBIT2 = ORBIT
    call tmali2(el,orbit2,al_errors,beta,gamma,orbit,re)
    RT = matmul(RE,RT)
    !--- IT check if RT is symplectic
    if (thr_on .eq. 0 .and. debug .ne. 0) then
      call m66symp(rt,nrm)
      if ((nrm-nrm0) .gt. symp_thrd) then
        call element_name(el_name,len(el_name))
        write (warnstr,'(a,e13.6,a,a)') "Symplectic deviation: ", nrm, " in element ", el_name
        call fort_warn('THREADER-2: ', warnstr)
      endif
      nrm0 = nrm
    endif
  endif

  if (kobs.gt.0 .and. kobs.eq.nobs) return

  if (thr_on .gt. 0)  then !---  threader is on;  keep matrix
     select case (code)
     case (code_hkicker)
            CMATR(1:6,1:6,1) = RT(1:6,1:6)
            j = name_len
            call element_name(c_name(1),j)
            coc_cnt(1) = node_value('occ_cnt ')
     case (code_kicker)
            CMATR(1:6,1:6,1) = RT(1:6,1:6)
            CMATR(1:6,1:6,2) = RT(1:6,1:6)
            j = name_len
            call element_name(c_name(1),j)
            coc_cnt(1) = node_value('occ_cnt ')
            call element_name(c_name(2),j)
            coc_cnt(2) = node_value('occ_cnt ')
     case (code_vkicker)
            CMATR(1:6,1:6,2) = RT(1:6,1:6)
            j = name_len
            call element_name(c_name(2),j)
            coc_cnt(2) = node_value('occ_cnt ')
     end select
  endif

  if (code .ge. code_monitor - 1 .and. code .le. code_monitor + 1)  then !---  monitors (code 17 to 19)
     enable = node_value('enable ')
     if (save .gt. 0 .and. enable .gt. 0) &
          call store_node_vector('orbit_ref ', 6, orbit)

     if (thr_on .gt. 0 .and. enable .gt. 0) then
        !--- threader is on; test for overflow w.r.t. stored orbit
        REFORB = zero
        call get_node_vector('orbit_ref ', 6, reforb)

        do kpro = 1, 2 ! projection plane: x or y
           if (      kpro .eq. 1 .and. code .le. code_monitor &        !--- plane is x and monitor or hmonitor
                .or. kpro .eq. 2 .and. code .ge. code_monitor)  then   !--- plane is y and monitor or vmonitor

              j = 2*kpro-1
              dorb(j) = orbit(j) - reforb(j) !--- orbit distortion

              if (abs(dorb(j)) .gt. vector(kpro) .and. node .ge. corr_pick(kpro))  then

                 if (debug .ne. 0) &
                      print *,' node = ',node,' kpro = ',kpro,' code = ',code,' dorb = ',dorb(j)

                 !---  reset count if new pickup
                 if (node .gt. corr_pick(kpro)) rep_cnt(kpro) = 0

                 !---  check for max. repetition
                 rep_cnt(kpro) = rep_cnt(kpro) + 1

                 if (rep_cnt(kpro) .gt. max_rep)  then
                    write(warnstr,'(a,i6,a)') 'pickup skipped after ',max_rep,' correction attempts'
                    call fort_warn('THREADER: ',warnstr)
                    goto 20
                 endif

                 !---  keep matrix
                 PMATR = RT
                 old = node
                 j = name_len
                 call element_name(p_name,j)
                 poc_cnt = node_value('occ_cnt ')
                 call tmthrd(kpro,dorb,cmatr(1,1,kpro),pmatr,vector,node,cick,err)

                 if (err .eq. 1)  then
                    write (warnstr, '(a,i6)') 'no corrector before pickup at node ', old
                    call fort_warn('THREADER: ',warnstr)

                 elseif (err .eq. 2)  then
                    call fort_warn('THREADER: ','kicker is at start of sequence')

                 elseif (err .eq. 0)  then
                    corr_pick(kpro) = old
                    write(tmptxt1, '(1p,6g14.6)') cick
                    tmptxt2 = c_name(kpro)
                    write(tmptxt2(lastnb(c_name(kpro))+1:), '(''['',i2,'']'')') coc_cnt(kpro)
                    tmptxt3 = p_name
                    write(tmptxt3(lastnb(p_name)+1:), '(''['',i2,'']'')') poc_cnt
                    call fort_info('-threader-','pickup: ' //tmptxt3(:lastnb(tmptxt3)) &
                         // '  kicker: '//tmptxt2(:lastnb(tmptxt2))//' total ' &
                         //ptxt(kpro)//'kick:'//tmptxt1)

                    !---  restore restart values
                    suml = restsum(kpro)
                    ORBIT = RESTORB(:,kpro)
                    RT = RESTM(:,:,kpro)
                    TT = RESTT(:,:,:,kpro)
                    goto 20 ! break the kpro loop
                 endif

              endif
           endif
        enddo
     endif
  endif

20 continue

  !---- Test for overflow.
  if (maxval(abs(ORBIT)) .ge. orb_limit) then
     eflag=1
     return
  endif

  if (advance_node().ne.0) then
     node=node+1
     goto 10 ! loop over nodes
  endif


  if(orbit(6) .gt. 1e-10 .and. istaper) then
    endpt=endpt+orbit(6)
    orderrun = orderrun+1
    goto 111
  endif
  bbd_flag=0

end SUBROUTINE tmfrst

SUBROUTINE tmthrd(kpro,dorb,cmatr,pmatr,thrvec,node,cick,error)
  use matrices
  use code_constfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Correct orbit position at bpm in first pass (threader)           *
  !     Input:                                                           *
  !     kpro       (integer) projection (1 = x, 2 = y)                   *
  !     dorb       (d.p.)    difference orbit at pickup                  *
  !     cmatr      (d.p.)    transfer matrix up to preceding corrector   *
  !     pmatr      (d.p.)    transfer matrix up to current pickup        *
  !     thrvec     (d.p.)    threader parameter vector (see madxdict.h)  *
  !     Input/Output:                                                    *
  !     node       (integer) current node                                *
  !     Output:                                                          *
  !     cick         (d.p.)    kick value applied                        *
  !     error        (integer) 0: OK, else aborted                       *
  !----------------------------------------------------------------------*
  integer :: kpro, node
  double precision ::  dorb(6), cmatr(6,6), pmatr(6,6), thrvec(3)
  double precision :: cick, error

  integer :: i, j, itp, lc, lc1, lc2, npick, code, ncorr
  double precision :: atemp(6,6)

  integer, external :: advance_node, retreat_node
  double precision, external :: node_value
  double precision, parameter :: tol1min=1.d-2

  error = 0

  !---  keep position of current pickup
  npick = node

  itp = code_kicker + (-1)**kpro ! 14 for x plane or 16 for y plane

  lc  = 2 * (kpro - 1)
  lc1 = lc + 1
  lc2 = lc + 2

  !---  look for one preceding corrector of correct (MAD-8) type
10 continue

  if (retreat_node() .eq. 0)  then !--- we are at start already
     do i = node+1, npick
        j = advance_node()
     enddo
     node = npick
     error = 1
     return
  endif

  node = node - 1
  code = node_value('mad8_type ')

!  if (code .eq. code_tkicker)     code = code_kicker
  if (code .eq. code_placeholder) code = code_instrument

  !-- if wrong corrector or not a corrector, loop
  if (code .ne. itp .and. code .ne. code_kicker) goto 10

  !--- corrector found for the projection plane
  ncorr = node

  !---  transport matrix from kicker to pickup
  ATEMP = matmul(JMATT, matmul(transpose(CMATR),JMAT)) ! invert symplectic matrix
  ATEMP = matmul(PMATR,ATEMP)

  !--- if kicker is not efficient enough - loop
  if (abs(atemp(lc1,lc2)) .lt. tol1min) goto 10

  !---  now we got one good corrector - get kick with attenuation factor
  cick = -thrvec(3) * dorb(2*kpro-1) / atemp(lc1,lc2)

  !---  add kick to kicker strengths
  if (kpro .eq. 1)  then
     cick = node_value('other_bv ') * cick + node_value('chkick ')
     call store_node_value('chkick ', cick)
  else if (kpro .eq. 2) then
     cick = node_value('other_bv ') * cick + node_value('cvkick ')
     call store_node_value('cvkick ', cick)
  endif

  !---  set node for restart in front of kicker
  if (retreat_node() .eq. 0)  then !--- we are at start already
     do i = node+1, npick
        j = advance_node()
     enddo
     node = npick
     error = 2
     return
  endif
  node = node -1

end SUBROUTINE tmthrd

SUBROUTINE twcpin(rt,disp0,r0mat,eflag)
  use twiss0fi
  use twisscfi
  use twissbeamfi, only : deltap
  use matrices, only : EYE, SMAT, symp_thrd
  use math_constfi, only : zero, one, two, quarter
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values for linear coupling parameters.                   *
  !     Entry point for mad_emit                                         *
  !     Input:                                                           *
  !     rt(6,6)      (double)  one turn transfer matrix.                 *
  !     Output:                                                          *
  !     disp0        (double)  initial dispersion vector                 *
  !     r0mat(2,2)   (double)  coupling matrix                           *
  !     eflag        (integer) error flag.                               *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: rt(6,6)
  double precision, intent(OUT) :: disp0(6), r0mat(2,2)
  integer, intent(OUT) :: eflag

  logical :: stabx, staby
  double precision :: a(2,2), b(2,2), c(2,2), d(2,2), ra(6,6), r_eig(6,6)
  double precision :: e(2,2), aux(2,2), f(2,2), bbar(2,2)
  double precision :: arg, den, det, r_det, dtr, sinmu2
  double precision :: betx0, alfx0, amux0
  double precision :: bety0, alfy0, amuy0
  character(len=150) :: msg

  integer, external :: get_option
  double precision, external :: get_value
  double precision, parameter :: eps=1d-8, diff_cos = 1d-5

  double precision  :: em(6,6),  cosmu1_eig, cosmu2_eig, nrm
  double precision  :: reval(6), aival(6) ! re and im parts
  logical, external :: m66sta
  character(len=180):: warnstr

  !--- initialize deltap because twcpin can be called directly from mad_emit
  deltap = get_value('probe ','deltap ')

  !---- Initialization
  eflag   = 0
  gammacp = one
  RA = RT
  R_EIG = eye
  betx0=zero; alfx0=zero; amux0=zero
  bety0=zero; alfy0=zero; amuy0=zero
  cosmu1_eig=zero; cosmu2_eig=zero
  stabx=.false.; staby=.false.

  call m66symp(rt,nrm)
  if (nrm .gt. symp_thrd) then
     write (warnstr,'(a, e13.6, a)') "One-turn map R symplectic deviation: ", &
                                     nrm, " (symplectifying R)"
     call fort_warn('TWCPIN: ', warnstr)
  endif

  if (nrm .gt. zero) then
    call tmsymp(RA)
  endif

  ! RA(1:4,1:4) = ( A , B
  !                 C , D)
  A = RA(1:2,1:2) ; B = RA(1:2,3:4)
  C = RA(3:4,1:2) ; D = RA(3:4,3:4)

  !---- Initial dispersion.
  call twdisp_ini(RA, disp0) !-- LD: 2016.04.18

  !---- Matrix C + B(bar) and its determinant (for R_A)
  BBAR = matmul(matmul(-SMAT,transpose(B)),SMAT)
  aux = C + BBAR
  det = aux(1,1) * aux(2,2) - aux(1,2) * aux(2,1)

  !---- Coupling matrix.
  !--- (trA - trD) / 2
  dtr = (A(1,1) + A(2,2) - D(1,1) - D(2,2)) / two

  !--- det|C+BBAR| + [(trA - trD) / 2]**2:
  arg = det + dtr**2

  if (arg .lt. zero) then
     !---- Unstable due to coupling.
     stabx = .false. ; staby = .false.

  else
     if (arg .eq. zero) then
        R0MAT(:2,:2) = EYE(:2,:2)
     else
        den     = - (dtr + sign(sqrt(arg),dtr))
        R0MAT   = AUX / den
     endif

     r_det     = r0mat(1,1) * r0mat(2,2) - r0mat(1,2) * r0mat(2,1)
     gammacp   = 1 / sqrt( 1 + r_det )


     ! !---- Decouple: Find diagonal blocks.
     ! a(1,1) = rt(1,1) - rt(1,3)*r0mat(1,1) - rt(1,4)*r0mat(2,1)
     ! a(1,2) = rt(1,2) - rt(1,3)*r0mat(1,2) - rt(1,4)*r0mat(2,2)
     ! a(2,1) = rt(2,1) - rt(2,3)*r0mat(1,1) - rt(2,4)*r0mat(2,1)
     ! a(2,2) = rt(2,2) - rt(2,3)*r0mat(1,2) - rt(2,4)*r0mat(2,2)

     ! d(1,1) = rt(3,3) + r0mat(1,1)*rt(1,3) + r0mat(1,2)*rt(2,3)
     ! d(1,2) = rt(3,4) + r0mat(1,1)*rt(1,4) + r0mat(1,2)*rt(2,4)
     ! d(2,1) = rt(4,3) + r0mat(2,1)*rt(1,3) + r0mat(2,2)*rt(2,3)
     ! d(2,2) = rt(4,4) + r0mat(2,1)*rt(1,4) + r0mat(2,2)*rt(2,4)

     !--The same as above a and d, just rewriten in matrix formalism (e = a, f = d)
     !-- In d callucaltion to have the same values as below in f one has to add parenthesis around R0MAT*B
     ! i.e d(1,1) = rt(3,3) +  ( r0mat(1,1)*rt(1,3) + r0mat(1,2)*rt(2,3) )
     e = A - matmul( B,     R0MAT )
     f = D + matmul( R0MAT, B     )

     !---- First mode.
     cosmux = (e(1,1) + e(2,2)) / two
     stabx = abs(cosmux).lt.one
     if (stabx) then
        sinmu2 = - e(1,2)*e(2,1) - quarter*(e(1,1) - e(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmux = sign(sqrt(sinmu2), e(1,2))
        betx0 = e(1,2) / sinmux
        alfx0 = (e(1,1) - e(2,2)) / (two * sinmux)
     else
        betx0 = zero
        alfx0 = zero
     endif

     !---- Second mode.
     cosmuy = (f(1,1) + f(2,2)) / two
     staby = abs(cosmuy).lt.one
     if (staby) then
        sinmu2 = - f(1,2)*f(2,1) - quarter*(f(1,1) - f(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmuy = sign(sqrt(sinmu2), f(1,2))
        bety0 = f(1,2) / sinmuy
        alfy0 = (f(1,1) - f(2,2)) / (two * sinmuy)
     else
        bety0 = zero
        alfy0 = zero
     endif

  endif

  ! Create decoupled matrix R_EIG =[E 0; 0 F ] for eigenvalues calculation
   R_EIG(1:2, 1:2) = E
   R_EIG(3:4, 3:4) = F

  !---- Find eigenvectors at initial position.
  reval = zero
  aival = zero


  if (get_option('debug ') .ne. 0) then
    call laseig(r_eig, reval, aival, em)
    cosmu1_eig = ( reval(1)+ aival(1) + reval(2) + aival(2) )/ 2
    cosmu2_eig = ( reval(3)+ aival(3) + reval(4) + aival(4) )/ 2

    write (warnstr,'(a,e13.6,a,e13.6)') "cosmux =  ", cosmux, ", cosmuy =", cosmuy
    call fort_warn('TWCPIN: ', warnstr)
    write (warnstr,'(a,e13.6,a,e13.6)')  "cosmu1_eig =", cosmu1_eig,  ", cosmu2_eig =", cosmu2_eig
    call fort_warn('TWCPIN: ', warnstr)

    if (.not. ((abs(cosmux - cosmu1_eig) .lt. diff_cos .and. abs(cosmuy - cosmu2_eig) .lt. diff_cos) .or.  &
          (abs(cosmuy - cosmu1_eig) .lt. diff_cos .and. abs(cosmux - cosmu2_eig) .lt. diff_cos))) then
        write (warnstr,'(a)') "Difference in the calculation of cosmux/cosmuy based of R_EIG eigen values!  "
        call fort_warn('TWCPIN: ', warnstr)
        write (warnstr,'(a,e13.6, a, e13.6)') "cosmux-cosmu1_eig =", cosmux-cosmu1_eig, "cosmux-cosmu2_eig =", cosmux-cosmu2_eig
        call fort_warn('TWCPIN: ', warnstr)
        write (warnstr,'(a,e13.6, a, e13.6)') "cosmuy-cosmu1_eig =", cosmuy-cosmu1_eig, "cosmuy-cosmu2_eig =", cosmuy-cosmu2_eig
        call fort_warn('TWCPIN: ', warnstr)
    endif
  endif

  ! call twcpin_print(rt,r0mat)

  !---- Give message, if unstable.
  if (.not.stabx .or. .not.staby) then
     if (.not.stabx .and. .not.staby) then
        write (msg,'(3(a,f12.6))') "Both modes are unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     elseif (.not.stabx) then
        write (msg,'(3(a,f12.6))') "Mode 1 is unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     elseif (.not.staby) then
        write (msg,'(3(a,f12.6))') "Mode 2 is unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     endif
     call fort_warn('TWCPIN: ',msg)
     eflag = 1
  endif

  opt_fun0(3)=betx0
  opt_fun0(4)=alfx0
  opt_fun0(5)=amux0
  opt_fun0(6)=bety0
  opt_fun0(7)=alfy0
  opt_fun0(8)=amuy0

  OPT_FUN0(15:18) = DISP0(1:4)

  opt_fun0(29)=r0mat(1,1)
  opt_fun0(30)=r0mat(1,2)
  opt_fun0(31)=r0mat(2,1)
  opt_fun0(32)=r0mat(2,2)

end SUBROUTINE twcpin

SUBROUTINE twcpin_talman(rt,disp0,r0mat,eflag)
  use twiss0fi
  use twisscfi
  use twissbeamfi, only : deltap
  use matrices, only : EYE, SMAT
  use math_constfi, only : zero, one, two, quarter
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values for linear coupling parameters.                   *
  !     Irina Tecker: Talman's approach                                  *
  !     Entry point for mad_emit                                         *
  !     Input:                                                           *
  !     rt(6,6)      (double)  one turn transfer matrix.                 *
  !     Output:                                                          *
  !     disp0        (double)  initial dispersion vector                 *
  !     r0mat(2,2)   (double)  coupling matrix                           *
  !     eflag        (integer) error flag.                               *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: rt(6,6)
  double precision, intent(OUT) :: disp0(6), r0mat(2,2)
  integer, intent(OUT) :: eflag

  logical :: stabx, staby
  double precision :: a(2,2),  b(2,2),  c(2,2),  d(2,2)
  double precision :: bbar(2,2), r0mat_bar(2,2)
  double precision :: e(2,2), aux(2,2), f(2,2)
  double precision :: arg, den, det, dtr, sinmu2, r_det
  double precision :: betx0, alfx0, amux0
  double precision :: bety0, alfy0, amuy0
  character(len=150) :: msg

  double precision :: det_e, det_f

  integer, external :: get_option
  double precision, external :: get_value
  double precision, parameter :: eps=1d-8

  !--- initialize deltap because twcpin can be called directly from mad_emit
  deltap = get_value('probe ','deltap ')

  !---- Initialization
  stabx=.false.; staby=.false.
  eflag=0
  gammacp=one
  betx0=zero; alfx0=zero; amux0=zero
  bety0=zero; alfy0=zero; amuy0=zero
  det_e=zero; det_f=zero

  !---- Initial dispersion.
  call twdisp_ini(rt,disp0) !-- LD: 2016.04.18

  ! RT(1:4,1:4) = ( A , B
  !                 C , D)
  A = RT(1:2,1:2) ; B = RT(1:2,3:4)
  C = RT(3:4,1:2) ; D = RT(3:4,3:4)

  !---- Matrix C + B(bar) and its determinant (for R_A)
  bbar = matmul(matmul(-SMAT,transpose(b)),SMAT)
  aux(:2,:2) = C(:2,:2) + BBAR(:2,:2)
  det = aux(1,1) * aux(2,2) - aux(1,2) * aux(2,1)


  !---- Coupling matrix.
  !--- (trA - trD) / 2
  dtr = (A(1,1) + A(2,2) - D(1,1) - D(2,2)) / two
  !--- det|C+BBAR| + [(trA - trD) / 2]**2:
  arg = det + dtr**2


  if (arg .lt. zero) then
     !---- Unstable due to coupling.
     stabx = .false. ; staby = .false.

  else
     if (arg .eq. zero) then
        R0MAT(:2,:2) = EYE(:2,:2)
        r_det     = r0mat(1,1) * r0mat(2,2) - r0mat(1,2) * r0mat(2,1)
        gammacp   = 1 / sqrt( 1 + r_det )

     else
        ! --- - [(trA - trD) / 2 + sign(trA - trD)*sqrt(det|C+BBAR| + [(trA - trD) / 2]**2)]
        den     = ( dtr + sign(sqrt(arg),dtr) ) ! --- (lamda_A - trd) p.13, Talman to have -R_A = (C+BBAR)/ - (lamda_A - trd)
        gammacp = sqrt(abs(den)/abs(-2*sign(sqrt(arg),dtr)))! -- numerical coefitient in front of R0MAT (p.16 ,Talman)

        ! V = [ gammacp* I, RD ; -RD_BAR, gammacp*I]
        !R0MAT = gammacp*AUX / den  ! ---  RD
        !R0MAT_BAR = - matmul(matmul(-SMAT,transpose(R0MAT)),SMAT) ! ---- RA = - RD_BAR
        R0MAT = - AUX / den  ! --- - RA
     endif
     R0MAT_BAR =  matmul(matmul(-SMAT,transpose(R0MAT)),SMAT) ! --- RD = - RA_BAR

     !---- Decouple: Find diagonal blocks.
     !A : g2(A + BRA - RDC - RFRA)
     !D = g2(-RERD - RAB + CRD + D ) .
     !     e =   gammacp**2 * A    + gammacp * matmul(B,R0MAT_BAR) - gammacp * matmul(R0MAT,C) - matmul(matmul(R0MAT,D),R0MAT_BAR)
     !     f =  -matmul(matmul(R0MAT_BAR,A),R0MAT) - gammacp * matmul(R0MAT_BAR,B) + gammacp * matmul(C,R0MAT) + gammacp**2 * D
     e =   gammacp**2 * ( A    -  matmul(B,R0MAT) - matmul(R0MAT_BAR,C) + matmul(matmul(R0MAT_BAR,D),R0MAT))
     f =   gammacp**2 * (matmul(matmul(R0MAT,A),R0MAT_BAR) +  matmul(R0MAT,B) + matmul(C,R0MAT_BAR) + D)

     det_e = e(1,1) * e(2,2) - e(1,2) * e(2,1)
     det_f = f(1,1) * f(2,2) - f(1,2) * f(2,1)

     !---- First mode.
     cosmux = (e(1,1) + e(2,2)) / two
     stabx = abs(cosmux).lt.one
     if (stabx) then
        sinmu2 = - e(1,2)*e(2,1) - quarter*(e(1,1) - e(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmux = sign(sqrt(sinmu2), e(1,2))
        betx0 = e(1,2) / sinmux
        alfx0 = (e(1,1) - e(2,2)) / (two * sinmux)
     else
        betx0 = zero
        alfx0 = zero
     endif

     !---- Second mode.
     cosmuy = (f(1,1) + f(2,2)) / two
     staby = abs(cosmuy).lt.one
     if (staby) then
        sinmu2 = - f(1,2)*f(2,1) - quarter*(f(1,1) - f(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmuy = sign(sqrt(sinmu2), f(1,2))
        bety0 = f(1,2) / sinmuy
        alfy0 = (f(1,1) - f(2,2)) / (two * sinmuy)
     else
        bety0 = zero
        alfy0 = zero
     endif

  endif

  !---- Give message, if unstable.
  if (.not.stabx .or. .not.staby) then
     if (.not.stabx .and. .not.staby) then
        write (msg,'(3(a,f12.6))') "Both modes are unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     elseif (.not.stabx) then
        write (msg,'(3(a,f12.6))') "Mode 1 is unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     elseif (.not.staby) then
        write (msg,'(3(a,f12.6))') "Mode 2 is unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     endif
     call fort_warn('TWCPIN: ',msg)
     eflag = 1
  endif

  opt_fun0(3)=betx0
  opt_fun0(4)=alfx0
  opt_fun0(5)=amux0
  opt_fun0(6)=bety0
  opt_fun0(7)=alfy0
  opt_fun0(8)=amuy0

  OPT_FUN0(15:18) = DISP0(1:4)

  opt_fun0(29)=r0mat(1,1)
  opt_fun0(30)=r0mat(1,2)
  opt_fun0(31)=r0mat(2,1)
  opt_fun0(32)=r0mat(2,2)

end SUBROUTINE twcpin_talman

SUBROUTINE twcpin_sagan(rt,disp0,r0mat, eflag)
  use twiss0fi
  use twisscfi
  use twissbeamfi, only : deltap
  use matrices, only : SMAT, SMATT, EYE
  use math_constfi, only : zero, one, two, four, half, quarter
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values for linear coupling parameters.                   *
  !     version from Sagan formalism by Irina Tecker                     *
  !     Input:                                                           *
  !     rt(6,6)      (double)  one turn transfer matrix.                 *
  !     Output:                                                          *
  !     disp0        (double)  initial dispersion vector                 *
  !     r0mat(2,2)   (double)  coupling matrix                           *
  !     gammacp      (double)  coefficient for the coupling matrix       *
  !     eflag        (integer) error flag.                               *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: rt(6,6)
  double precision, intent(OUT) :: disp0(6), r0mat(2,2)
  integer, intent(OUT) :: eflag
  integer :: i
  logical :: stabx, staby
  double precision :: e(2,2), aux(2,2), f(2,2), r0mat_bar(2,2)
  double precision :: a(2,2), b(2,2), c(2,2), d(2,2), cbar(2,2)
  double precision :: v(4, 4), u (4,4), vbar(4,4)
  double precision :: arg, den, det, dtr, sinmu2, det_e, det_f
  double precision :: betx0, alfx0, amux0
  double precision :: bety0, alfy0, amuy0
  character(len=150) :: msg

  integer, external :: get_option
  double precision, external :: get_value
  double precision, parameter :: eps=1d-8

  !---- Initialization
  stabx=.false.; staby=.false.
  eflag = 0
  gammacp = one
  betx0=zero; alfx0=zero; amux0=zero
  bety0=zero; alfy0=zero; amuy0=zero
  deltap   = get_value('probe ','deltap ')

  !---- Initial dispersion.
  call twdisp_ini(rt,disp0)  !-- LD: 2016.04.18

  ! RT(1:4,1:4) = ( A , B
  !                 C , D)
  A = RT(1:2,1:2) ; B = RT(1:2,3:4)
  C = RT(3:4,1:2) ; D = RT(3:4,3:4)

  !---- Matrix H = B + C(bar) = m + n(bar) (Sagan) and its determinant.
  cbar = matmul(matmul(-SMAT,transpose(c)),SMAT)
  aux(:2,:2) = B(:2,:2) + CBAR(:2,:2)
  !-- det(H) = det (B + C(bar))
  det = aux(1,1) * aux(2,2) - aux(1,2) * aux(2,1)

  !---- Coupling matrix.
  !---- trace(M-N)(Sagan) = trace(A-D)
  dtr = A(1,1) + A(2,2) - D(1,1) - D(2,2)
  !--- det|C+BBAR| + [(trA - trD) / 2]**2:
  arg = dtr**2 + four * det

  !---- If arg <= 0, the motion is unstable, we are in a stop band
  if (arg .le. zero) then

     !---- Unstable due to coupling.
     stabx = .false. ; staby = .false.

  else
     gammacp = sqrt(half + half*sqrt(dtr**2/arg))
     den     =  gammacp * sqrt(arg)
     R0MAT   = - sign(AUX, dtr) / den !-- r0mat = C (Sagan)

     !--R0MAT_BAR = SRMAT^{T}S^{T} - symplectic conjugate of R0MAT
     R0MAT_BAR = matmul(-SMAT, matmul(transpose(R0MAT),SMAT))


     ! --Compute matrix V ("rotation", [ gammacp*I R0MAT ; -R0MAT_BAR gammacp*I]
     V(1:4, 1:4)= zero
     do i = 1, 4
        V(i,i) = gammacp
     enddo

     V(1:2,3:4) =  R0MAT
     V(3:4,1:2) = -R0MAT_BAR

     ! --Compute matrix inv(V)
     VBAR(1:4, 1:4)= zero
     do i = 1, 4
          VBAR(i,i) = gammacp
     enddo
     VBAR(1:2,3:4) = -R0MAT
     VBAR(3:4,1:2) = R0MAT_BAR

     ! --Compute uncoupled block-diagonal U matrix [E 0, F 0]

     U = matmul (matmul (VBAR, RT(1:4,1:4)), V)

     !---- Find diagonal blocks.

     E = U(1:2, 1:2)
     F = U(3:4, 3:4)

     !---- First mode.
     cosmux = (e(1,1) + e(2,2)) / two
     stabx = abs(cosmux).lt.one
     if (stabx) then
        sinmu2 = - e(1,2)*e(2,1) - quarter*(e(1,1) - e(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmux = sign(sqrt(sinmu2), e(1,2))
        betx0 = e(1,2) / sinmux
        alfx0 = (e(1,1) - e(2,2)) / (two * sinmux)
     else
        betx0 = zero
        alfx0 = zero
     endif

     !---- Second mode.
     cosmuy = (f(1,1) + f(2,2)) / two
     staby = abs(cosmuy).lt.one
     if (staby) then
        sinmu2 = - f(1,2)*f(2,1) - quarter*(f(1,1) - f(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmuy = sign(sqrt(sinmu2), f(1,2))
        bety0 = f(1,2) / sinmuy
        alfy0 = (f(1,1) - f(2,2)) / (two * sinmuy)
     else
        bety0 = zero
        alfy0 = zero
     endif

  endif

  !---- Give message, if unstable.
  if (.not.stabx .or. .not.staby) then
     if (.not.stabx .and. .not.staby) then
        write (msg,'(3(a,f12.6))') "Both modes are unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     elseif (.not.stabx) then
        write (msg,'(3(a,f12.6))') "Mode 1 is unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     elseif (.not.staby) then
        write (msg,'(3(a,f12.6))') "Mode 2 is unstable for delta(p)/p = ",deltap, &
             ": cosmux = ",cosmux,", cosmuy = ",cosmuy
     endif
     call fort_warn('TWCPIN_SAGAN: ',msg)
     eflag = 1
  endif

  opt_fun0(3)=betx0
  opt_fun0(4)=alfx0
  opt_fun0(5)=amux0
  opt_fun0(6)=bety0
  opt_fun0(7)=alfy0
  opt_fun0(8)=amuy0

  OPT_FUN0(15:18) = DISP0(1:4)

  opt_fun0(29)=r0mat(1,1)
  opt_fun0(30)=r0mat(1,2)
  opt_fun0(31)=r0mat(2,1)
  opt_fun0(32)=r0mat(2,2)

end SUBROUTINE twcpin_sagan

SUBROUTINE twdisp(rt,vect,disp)
  use matrices, only : EYE
  use math_constfi, only : zero, one
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values for dispersion or its first derivative by delta.  *
  !     Only the first four components of the vectors are set.           *
  !     Input:                                                           *
  !     rt(6,6)   (double)    one turn transfer matrix.                  *
  !     vect(6)   (double)    right-hand side:                           *
  !     column 6 of rt for dispersion,                                   *
  !     auxiliary vector for derivative of dipersion.                    *
  !     Output:                                                          *
  !     disp(6)   (double)    dispersion vector.                         *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6), vect(6), disp(6)
  integer, external :: get_option

  integer :: i, j, irank
  double precision :: a(4,5)

  A(:4,:4) = RT(:4,:4) - EYE(:4,:4)
  A(1:4,5) = -VECT(1:4)

  call solver(a,4,1,irank)

  if (irank.ge.4) then
     DISP(1:4) = A(1:4,5)
  else
     if (get_option('info ') .ne. 0) then
       print *, 'TWDISP: Unable to compute intial dispersion --- dispersion set to zero.'
     endif
     DISP(1:4) = zero
  endif

end SUBROUTINE twdisp

SUBROUTINE twdisp_ini(rt,disp0)
  use math_constfi, only : zero, one
  use twiss0fi
  use twisscfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values for dispersion.                                   *
  !     Input:                                                           *
  !     rt(6,6)   (double)    one turn transfer matrix.                  *
  !     Output:                                                          *
  !     disp(6)   (double)    dispersion vector.                         *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6), disp0(6)
  integer, external :: get_option

  !---- Initial dispersion.
  if (get_option('twiss_inval ') .eq. 0) then
     call twdisp(rt,rt(:,6),disp0)
  else
     DISP0(1:4) = OPT_FUN(15:18)
  endif
  disp0(5) = zero
  disp0(6) = one

end SUBROUTINE twdisp_ini

SUBROUTINE twcpgo(rt,orbit0)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twiss_elpfi
  use twissotmfi
  use twissbeamfi, only : radiate, beta, gamma
  use spch_bbfi
  use name_lenfi
  use matrices, only: EYE
  use math_constfi, only : zero, one, two
  use code_constfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track Twiss parameters with optional output, version with        *
  !     coupling.                                                        *
  !     Input:                                                           *
  !     rt(6,6)   (double)  one turn transfer matrix.                    *
  ! 2014-Jun-11  11:10:37  ghislain: added:                              *
  !     orbit0(6) (double)  initial orbit vector                         *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6), orbit0(6)

  logical :: fmap, cplxy, cplxt, sector_sel
  integer :: i, i1, i2, iecnt, code, save, n_align, elpar_vl
  double precision :: ek(6), re(6,6), rwi(6,6), rc(6,6), te(6,6,6)
  double precision :: orbit00(6), ek00(6), re00(6,6), te00(6,6,6), disp00(6)
  double precision :: rw0(6,6), rmat0(2,2), sigmat00(6,6), pos0
  double precision :: srmat0(6,6)
  double precision :: alfx0, betx0, amux0
  double precision :: alfy0, bety0, amuy0
  double precision :: orbit(6), orbit2(6)
  double precision :: bvk, sumloc, sd, el, dl, currpos
  double precision :: al_errors(align_max)
  ! character(len=name_len) el_name
  character(len=130) :: msg

  integer, external :: el_par_vector, advance_node, restart_sequ, get_option, node_al_errors
  integer, external :: start_interp_node, fetch_interp_node
  double precision, external :: node_value, get_value
  double precision, parameter :: eps=1d-16

  !---- Initialization
  sumloc=zero
  amux=zero
  amuy=zero
  currpos=zero
  centre = get_value('twiss ','centre ').ne.zero
  RWI = EYE
  RC = EYE
  cplxy=.false.; cplxt=.false.

  !---- Store one-turn-map
  ROTM = RT

  !---- Create internal table for lattice functions if requested
  save = get_option('twiss_save ')
  OPT_FUN(:fundim) = OPT_FUN0(:fundim)

  !---- Initial values for lattice functions.
  betx    = opt_fun(3)
  alfx    = opt_fun(4)
  amux    = opt_fun(5)
  bety    = opt_fun(6)
  alfy    = opt_fun(7)
  amuy    = opt_fun(8)

  ORBIT     = OPT_FUN(9:14)
  DISP(1:4) = OPT_FUN(15:18)
  disp(5) = zero
  disp(6) = one

  rmat(1,1) = opt_fun(29)
  rmat(1,2) = opt_fun(30)
  rmat(2,1) = opt_fun(31)
  rmat(2,2) = opt_fun(32)

  ! IT sigma matrix(1:6, 1:6) = opt_fun(75:110)
  do i1=1,6
  do i2=1,6
    sigmat(i1,i2) = opt_fun(74 + (i1-1)*6 + i2)
  enddo
  enddo

  !--- 2014-May-30  15:19:52  ghislain: if initial values are provided, copy the initial orbit
  !                 because opt_fun contains only the values given on the twiss command itself
  !                 but does not know the values coming from COGUESS or USEORBIT
  if (get_option('twiss_inval ') .ne. 0) ORBIT = ORBIT0

  !---- Maximum and r.m.s. values.
  bxmax = betx
  dxmax = disp(1)
  bymax = bety
  dymax = disp(3)
  xcomax= zero
  ycomax= zero
  sigxco= zero
  sigyco= zero
  sigdx = zero
  sigdy = zero

  !---- Loop over positions.
  iecnt = 0
  i = restart_sequ()
  i_spch = 0

  i = 1
  do while (i .ne. 0)
    el = node_value('l ')
    if (start_interp_node(i) .ne. 0) then
      do while (fetch_interp_node(i, dl) .ne. 0)
        call backup_optics()
        call track_one_element(dl, .true., .true.)
        call restore_optics()
      end do
      call track_one_element(el, .false., .false.)
    else
      i = 1
      call track_one_element(el, .not. centre, .true.)
    endif
    i = advance_node()
  end do

  call compute_summary()

contains

subroutine track_one_element(el, fexit, contrib_rms)
  double precision, intent(in) :: el
  logical :: fexit
  logical :: contrib_rms

  sector_sel = node_value('sel_sector ') .ne. zero .and. sectormap
  code = node_value('mad8_type ')
!  if (code .eq. code_tkicker)     code = code_kicker
  if (code .eq. code_placeholder) code = code_instrument
  bvk = node_value('other_bv ')
  elpar_vl = el_par_vector(g_polarity, g_elpar)
  ele_body = el .gt. eps

  !--- 2013-Nov-14  10:34:00  ghislain: add acquisition of name of element here.
  !call element_name(el_name,len(el_name))

  opt_fun(70) = g_elpar(g_kmax)
  opt_fun(71) = g_elpar(g_kmin)
  opt_fun(72) = g_elpar(g_calib)
  opt_fun(73) = g_elpar(g_polarity)

  n_align = node_al_errors(al_errors)
  if (n_align .ne. 0)  then
     !print*, "coupl1: Element = "
     ele_body = .false.
     orbit2 = orbit
     call tmali1(orbit2,al_errors,beta,gamma,orbit,re)
     call twcptk(re,orbit)
     if (sectormap) SRMAT = matmul(RE,SRMAT)
  endif

  if (centre) then
    call backup_optics()

    call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te,.true.,el/two)

    if (fmap) call twcptk(re,orbit)

    call save_opt_fun()
    call twprep(save,1,opt_fun,currpos+el/two,i)

    call restore_optics()
  endif

  ! now do exact calculation with full length:
  call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te,.false.,el)

  if (fmap) then
     call twcptk(re,orbit)
     !print*, "couplfmap: Element = ", el_name
     if (sectormap) call tmcat(.true.,re,te,srmat,stmat,srmat,stmat)
  endif

  if (n_align .ne. 0)  then
     !print*, "coupl2: Element = ", el_name
     ele_body = .false.
     orbit2 = orbit
     call tmali2(el,orbit2,al_errors,beta,gamma,orbit,re)
     call twcptk(re,orbit)
     if (sectormap) SRMAT = matmul(RE,SRMAT)
  endif

  sumloc = sumloc + el
  if (sector_sel) call twwmap(sumloc, orbit)

  currpos=currpos+el

  ! Avoid adding the same value at the exit twice end when interpolate is on:
  if (contrib_rms) then
    iecnt=iecnt+1

    !--- 2013-Nov-14  14:18:07  ghislain: should only calculate the RMS values for active elements,
    !                           skipping markers(25), beambeam(22), changeref(35), translation(36),
    !                           srotation(12), yrotation(13) etc...
    !                           But leave drifts(1) in for now. They are the most numerous elements in machines
    sigxco = sigxco + orbit(1)**2
    sigyco = sigyco + orbit(3)**2
    sigdx  = sigdx  + disp(1)**2
    sigdy  = sigdy  + disp(3)**2
  endif

  call save_opt_fun()
  if (fexit) then
     call twprep(save,1,opt_fun,currpos,i)
  endif
end subroutine track_one_element

subroutine compute_summary
  wgt    = max(iecnt, 1)
  sigxco = sqrt(sigxco / wgt)
  sigyco = sqrt(sigyco / wgt)
  sigdx  = sqrt(sigdx / wgt)
  sigdy  = sqrt(sigdy / wgt)
  ! IT. remove cosmux/cosmuy calculation here as it's not used anywhere and only messes up check for coupled periodic lattices
  ! cosmux = (rt(1,1) + rt(2,2)) / two
  ! cosmuy = (rt(3,3) + rt(4,4)) / two

  !---- Warning messages.
  if (cplxt .or. radiate) &
       call fort_warn('TWCPGO: ','TWISS uses the RF system or synchrotron radiation only '// &
                       'to find the closed orbit, for optical calculations it ignores both.')
end subroutine compute_summary

subroutine backup_optics()
  orbit00 = orbit ; ek00 = ek ; re00 = re ; te00 = te
  pos0=currpos
  betx0=betx; alfx0=alfx; amux0=amux
  bety0=bety; alfy0=alfy; amuy0=amuy
  rmat0 = rmat ; disp00 = disp
  if (sectormap) srmat0 = srmat
  if (rmatrix) rw0 = rw
  sigmat00 = sigmat
end subroutine backup_optics

subroutine restore_optics()
  orbit = orbit00 ; ek = ek00 ; re = re00 ; te = te00
  currpos = pos0
  betx=betx0; alfx=alfx0; amux=amux0
  bety=bety0; alfy=alfy0; amuy=amuy0
  rmat = rmat0 ; disp = disp00
  if (sectormap) srmat = srmat0
  if (rmatrix) rw = rw0
  sigmat = sigmat00
end subroutine restore_optics

subroutine save_opt_fun()
  integer :: i1, i2

  !--- save maxima and name of elements where they occur
  bxmax  = max(bxmax,  betx)
  bymax  = max(bymax,  bety)
  dxmax  = max(dxmax,  abs(disp(1)))
  dymax  = max(dymax,  abs(disp(3)))
  xcomax = max(xcomax, abs(orbit(1)))
  ycomax = max(ycomax, abs(orbit(3)))

     opt_fun(3) = betx
     opt_fun(4) = alfx
     opt_fun(5) = amux
     opt_fun(6) = bety
     opt_fun(7) = alfy
     opt_fun(8) = amuy
     opt_fun(9:14) = ORBIT
     opt_fun(15:18) = disp(1:4)

     opt_fun(29) = rmat(1,1)
     opt_fun(30) = rmat(1,2)
     opt_fun(31) = rmat(2,1)
     opt_fun(32) = rmat(2,2)

     ! !IT sigma matrix(1:6, 1:6) = opt_fun(75:110)
     do i1=1,6
        do i2=1,6
           opt_fun(74 + (i1-1)*6 + i2) = sigmat(i1,i2)
        enddo
     enddo

  sd = rt(5,6) + dot_product(RT(5,1:4),DISP(1:4))
  eta = - sd * beta**2 / circ
  alfa = one / gamma**2 + eta
  opt_fun(74) = alfa
end subroutine save_opt_fun

end SUBROUTINE twcpgo

SUBROUTINE twcptk(re,orbit)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  use matrices, only : JMAT, JMATT, SMAT, SMATT
  use math_constfi, only : zero, one, two, twopi
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track coupled lattice functions.                                 *
  !     Irina Tecker : new version rewritten in matrix formalism         *
  !                   + added mode flip                                  *
  !                   + no "guard" for negative betx                     *
  !     Input:                                                           *
  !     re(6,6)  (double)   transfer matrix of element.                  *
  !     orbit(6) (double)   closed orbit                                 *
  !----------------------------------------------------------------------*
  double precision :: re(6,6), orbit(6)

  integer :: i, i1, i2, j, irank
  double precision :: rwi(6,6), rc(6,6), dt(6)
  double precision :: a(2,2), b(2,2), c(2,2), d(2,2), ra(4,8)
  double precision :: rmat0(2,2), e(2,2), f(2,2), cd(2,2), tmp(2,2)
  double precision :: rmat_bar(2,2), ebar(2,2), fbar(2,2), hmat(2,2)
  double precision :: edet, fdet, tempa, tempb, det
  double precision :: alfx_ini, betx_ini, amux_ini
  double precision :: alfy_ini, bety_ini, amuy_ini
  character(len=name_len) :: name
  character(len=200)      :: warnstr

  integer, external :: get_option
  double precision, parameter :: eps=1d-36

  logical :: mode_flip_ele
  logical :: cp_error

  double precision :: trace_e, trace_f
  integer eflag
  integer, parameter :: izero = 0
  !---- Initialization

  alfx_ini=zero; betx_ini=zero; amux_ini=zero
  alfy_ini=zero; bety_ini=zero; amuy_ini=zero
  mode_flip_ele = mode_flip
  cp_error=.false.
  eflag = 0

  call element_name(name,len(name))

  !---- Dispersion.
  DT = matmul(RE, DISP)

  ! RE(1:4,1:4) = ( A , B
  !                 C , D)
  A = RE(1:2,1:2) ; B = RE(1:2,3:4)
  C = RE(3:4,1:2) ; D = RE(3:4,3:4)

  !---- Check RE rank
  RA(1:4,1:4) = RE(1:4,1:4)
  RA(1:4,5:8) = zero
  RA(1,5)     = one
  RA(2,6)     = one
  RA(3,7)     = one
  RA(4,8)     = one

  call solver(ra,4,4,irank)
  if (irank.lt.4) then
     write (warnstr, '(a)') 'Coupling failed: Rank deficient matrix RE, irank < 4 '
     call fort_warn('TWCPTK: ', warnstr)
     return
  endif

  if   (all( B == zero  .and.  C == zero)  ) then

     E(1:2,1:2) = A(1:2,1:2)
     F(1:2,1:2) = D(1:2,1:2)

     ! symplectic conjugate of E = S*E^T*S^T
     EBAR = matmul(SMAT, matmul(transpose(E),SMATT))
     edet = e(1,1) * e(2,2) - e(1,2) * e(2,1)
     CD   = - matmul(F, RMAT)
     RMAT = - matmul(CD, EBAR) / edet

  else


     ! R is not simplectic but this RMAT_BAr is correct
     RMAT_BAR  = matmul(SMAT, matmul(transpose(RMAT),SMATT))
     TMP       = gammacp*(A -  matmul(B,RMAT))
     det       = tmp(1,1) * tmp(2,2) - tmp(1,2) * tmp(2,1)


     if (det .gt. 0.9 .or. (det .gt. 0.1 .and. .not. mode_flip) ) then  !--- normal mode
        E       = A - matmul( B,RMAT )                           ! former a
        edet    = e(1,1) * e(2,2) - e(1,2) * e(2,1)              ! former adet
        EBAR    = matmul( SMAT, matmul(transpose( E ), SMATT ))  ! symplectic conjugate of E = S*E^T*S^T
        CD      = C - matmul( D, RMAT )        ! former b
        F       = D + matmul( C, RMAT_BAR )    ! former c
        RMAT    = - matmul( CD, EBAR ) /  edet
        gammacp = sqrt( det )

     else  !---flip mode

        mode_flip_ele = .not. mode_flip
        nmode_flip = nmode_flip + 1

        TMP = gammacp*(matmul(A,RMAT_BAR) + B)
        det = tmp(1,1) * tmp(2,2) - tmp(1,2) * tmp(2,1)

        if (det .lt. zero) then
           write (warnstr, '(a,a)') 'Negative determinant in ', name
           call fort_warn('TWCPTK: ', warnstr)
           !print*, "Negative determinant! Rounding error?!?!? det = ", det
        endif

        F       = matmul (A, RMAT_BAR) + B
        fdet    = f(1,1) * f(2,2) - f(1,2) * f(2,1)
        FBAR    = matmul(SMAT, matmul(transpose(F),SMATT)) ! symplectic conjugate of F = S*F^T*S^T
        E       = C -  matmul(D, RMAT)
        CD      = matmul(C, RMAT_BAR) + D
        RMAT    = - matmul(CD, FBAR)/ abs( fdet )
        gammacp = sqrt (abs(det))

        write (warnstr, '(a, a, a, i3)') 'Mode flip in the element ', name, ', nflips up to now = ', nmode_flip
        call fort_warn('TWCPTK: ', warnstr)

        trace_e = E(1,1)+E(2,2)
        trace_f = F(1,1)+F(2,2)
     endif
  endif

  amux_ini = amux
  amuy_ini = amuy

  if(mode_flip) then
     call twcptk_twiss(f, e, cp_error)
  else
     call twcptk_twiss(e, f, cp_error)
  endif

  if (cp_error) then
     ! print *, '+++ det of block diagonal matrix is zero in ', name
     write (warnstr, '(a, a, a)') 'Det of block diagonal matrix is zero in ', name, ', twiss parameter might be unphysical! '
     call fort_warn('TWCPTK: ', warnstr)
  endif

  ! When we are comming out of a flipped mode, the phase is often off by a factor of twopi.
  ! Unfortunately there is no definitive way to calcuate what the "right" way to handle this is but "-twopi" is better than nothing ((c), Sagan)
  ! if(mode_flip .and. .not. mode_flip_ele ) then

  !   if ((amux-amux_ini) .ge. twopi ) then
  !     amux = amux - twopi
  !   elseif ((amux-amux_ini) .lt. zero) then
  !     write (warnstr,'(a,e13.6,a,a)') "Negative phase advance in x-plane (mode flip) of ", amux-amux_ini, " in element ", name
  !     call fort_warn('TWCPTK: ', warnstr)
  !     amux = amux + twopi
  !   endif

  !   if ((amuy-amuy_ini) .ge. twopi ) then
  !     amuy = amuy - twopi
  !   elseif ((amuy-amuy_ini) .lt. zero) then
  !     write (warnstr,'(a,e13.6,a,a)') "Negative phase advance in y-plane (mode flip) of ", amuy-amuy_ini, " in element ", name
  !     call fort_warn('TWCPTK: ', warnstr)
  !     amuy = amuy + twopi
  !   endif

  ! endif

  mode_flip = mode_flip_ele

  if (betx.lt.eps .or. bety.lt.eps) then
     write (warnstr, '(a, a, a)') 'Negative beta in ', name, ' Twiss parameter might be unphysical!'
     call fort_warn('TWCPTK: ', warnstr)
    ! print *, '+++ negative beta: name', name, 'betx=', betx, ', bety=', bety
  endif

  sigmat = matmul(RE, matmul(sigmat,transpose(RE)))

     ! !---- Auxiliary matrices.
     ! !LD:  a = re_x - re_xy*rmat
     ! a(1,1) = re(1,1) - (re(1,3) * rmat(1,1) + re(1,4) * rmat(2,1))
     ! a(1,2) = re(1,2) - (re(1,3) * rmat(1,2) + re(1,4) * rmat(2,2))
     ! a(2,1) = re(2,1) - (re(2,3) * rmat(1,1) + re(2,4) * rmat(2,1))
     ! a(2,2) = re(2,2) - (re(2,3) * rmat(1,2) + re(2,4) * rmat(2,2))

     ! b(1,1) = re(3,1) - (re(3,3) * rmat(1,1) + re(3,4) * rmat(2,1))
     ! b(1,2) = re(3,2) - (re(3,3) * rmat(1,2) + re(3,4) * rmat(2,2))
     ! b(2,1) = re(4,1) - (re(4,3) * rmat(1,1) + re(4,4) * rmat(2,1))
     ! b(2,2) = re(4,2) - (re(4,3) * rmat(1,2) + re(4,4) * rmat(2,2))

     ! c(1,1) = re(3,3) + (re(3,1) * rmat(2,2) - re(3,2) * rmat(2,1))
     ! c(1,2) = re(3,4) - (re(3,1) * rmat(1,2) - re(3,2) * rmat(1,1))
     ! c(2,1) = re(4,3) + (re(4,1) * rmat(2,2) - re(4,2) * rmat(2,1))
     ! c(2,2) = re(4,4) - (re(4,1) * rmat(1,2) - re(4,2) * rmat(1,1))


     ! !---- Track R matrix.
     ! adet = a(1,1) * a(2,2) - a(1,2) * a(2,1)
     ! if (abs(adet) .gt. eps) then
     !    rmat(1,1) = - (b(1,1) * a(2,2) - b(1,2) * a(2,1)) / adet
     !    rmat(1,2) =   (b(1,1) * a(1,2) - b(1,2) * a(1,1)) / adet
     !    rmat(2,1) = - (b(2,1) * a(2,2) - b(2,2) * a(2,1)) / adet
     !    rmat(2,2) =   (b(2,1) * a(1,2) - b(2,2) * a(1,1)) / adet

     !    !---- Mode 1.
     !    tempb = a(1,1) * betx - a(1,2) * alfx
     !    tempa = a(2,1) * betx - a(2,2) * alfx
     !    alfx = - (tempa * tempb + a(1,2) * a(2,2)) / (adet * betx)
     !    betx =   (tempb * tempb + a(1,2) * a(1,2)) / (adet * betx)
     !    if (abs(a(1,2)) .gt. eps) amux = amux + atan2(a(1,2),tempb)

     !    !---- Mode 2.
     !    tempb = c(1,1) * bety - c(1,2) * alfy
     !    tempa = c(2,1) * bety - c(2,2) * alfy
     !    alfy = - (tempa * tempb + c(1,2) * c(2,2)) / (adet * bety)
     !    bety =   (tempb * tempb + c(1,2) * c(1,2)) / (adet * bety)
     !    if (abs(c(1,2)) .gt. eps) amuy = amuy + atan2(c(1,2),tempb)
     ! else
     !   !LD: 09.2015
     !   call element_name(name,len(name))
     !   print *, 'coupling too strong for element ', name, ' (adet= ', adet, ')'
     !   call fort_warn('TWCPTK: ','twiss parameter might be unphysical (split element)')
     ! endif

  !---- Cumulative R matrix and one-turn map at element location.
  if (rmatrix) then
     RW = matmul(RE,RW)
     if (get_option('twiss_inval ') .ne. 0) then
        RC = RW
     else
        RWI = matmul(JMATT, matmul(transpose(RW),JMAT)) ! invert symplectic matrix
        RC = matmul(RW,matmul(ROTM,RWI))
     endif
  endif

     DISP(1:4) = DT(1:4)
     disp(5) = zero
     disp(6) = one

  if (rmatrix) then
    do i1=1,6
    do i2=1,6
      opt_fun(33 + (i1-1)*6 + i2) = rc(i1,i2)
    enddo
    enddo
  endif

end SUBROUTINE twcptk
SUBROUTINE twcptk_twiss(matx, maty, error)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  use math_constfi, only : zero, twopi
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Calculation of twiss parameters by Irina Tecker                  *
  !     Irina Tecker : new version rewritten in matrix formalism         *
  !                   + added mode flip                                  *
  !                   + no "guard" for negative betx                     *
  !     Input:                                                           *
  !     matx(2,2)  (double)   X-plane matrix of block-diagonal           *
  !     maty(2,2)  (double)   Y-plane matrix of block-diagonal           *
  !----------------------------------------------------------------------*

  double precision :: matx(2,2), maty(2,2)
  double precision :: matx11, matx12, matx21, matx22
  double precision :: maty11, maty12, maty21, maty22
  double precision :: alfx_ini, betx_ini, tempa
  double precision :: alfy_ini, bety_ini, tempb
  double precision :: detx, dety
  logical          :: error
  double precision, parameter :: eps=1d-36
  character(len=name_len) :: name
  character(len=180)      :: warnstr
  alfx_ini=zero; betx_ini=zero
  alfy_ini=zero; bety_ini=zero

  error = .true.

  detx = matx(1,1) * matx(2,2) - matx(1,2) * matx(2,1)
  dety = maty(1,1) * maty(2,2) - maty(1,2) * maty(2,1)

  call element_name(name,len(name))

  if (detx == 0) return
  if (dety == 0) return
  betx_ini = betx ; alfx_ini = alfx ;
  bety_ini = bety ; alfy_ini = alfy ;


  matx11 = matx(1,1)
  matx12 = matx(1,2)
  matx21 = matx(2,1)
  matx22 = matx(2,2)

  maty11 = maty(1,1)
  maty12 = maty(1,2)
  maty21 = maty(2,1)
  maty22 = maty(2,2)

  !---- Mode 1.
  tempb = matx11 * betx_ini - matx12 * alfx_ini
  tempa = matx21 * betx_ini - matx22 * alfx_ini
  alfx = - (tempa * tempb + matx12 * matx22) / (detx*betx_ini)
  betx =   (tempb * tempb + matx12 * matx12) / (detx*betx_ini)
  !if (abs(matx12).gt.eps) amux = amux + atan2(matx12,tempb)
  if (abs(matx12).gt.eps) then
     amux = amux + atan2(matx12,tempb)

     if (atan2(matx12,tempb) .lt. zero)then
        if (ele_body .and. abs(atan2(matx12,tempb)) > 0.1) then
           write (warnstr,'(a,e13.6,a,a)') "Negative phase advance in x-plane ", &
                atan2(matx12,tempb), " in the element ", name
           call fort_warn('TWCPTK_TWISS: ', warnstr)
           ! amux = amux + twopi
           ! print*, " matx12 =", matx12, " tempb = ", tempb , "matx11 * betx_ini =",  matx11 * betx_ini, &
           !      "  matx12 * alfx_ini = ",  matx12 * alfx_ini
           ! print*, " betx_ ini = " , betx_ini, " alfx_ini = ", alfx_ini
        endif
     endif
  endif




  !---- Mode 2.
  tempb = maty11 * bety_ini - maty12 * alfy_ini
  tempa = maty21 * bety_ini - maty22 * alfy_ini
  alfy = - (tempa * tempb + maty12 * maty22) / (detx*bety_ini)
  bety =   (tempb * tempb + maty12 * maty12) / (detx*bety_ini)
!  if (abs(maty12).gt.eps) amuy = amuy + atan2(maty12,tempb)
  if (abs(maty12).gt.eps) then
     amuy = amuy + atan2(maty12,tempb)

     if (atan2(maty12,tempb) .lt. zero) then
        if (ele_body .and. abs(atan2(maty12,tempb)) > 0.1 ) then
           write (warnstr,'(a,e13.6,a,a)') "Negative phase advance in y-plane ", &
                atan2(maty12,tempb), " in the element ", name
           call fort_warn('TWCPTK_TWISS: ', warnstr)
           ! print*, " maty12 =", maty12, " tempb = ", tempb , "maty11 * bety_ini =",  maty11 * bety_ini, &
           !      "  maty12 * alfy_ini = ",  maty12 * alfy_ini
           ! print*, " bety_ ini = " , bety_ini, " alfy_ini = ", alfy_ini

        endif
     endif
  endif

  error = .false.
END SUBROUTINE twcptk_twiss

SUBROUTINE tmsigma(s0mat)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  use math_constfi, only : zero, twopi, one, two
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Calculation of sigma (beam) by Irina Tecker                      *
  !     Couling is also included into the calculation following          *
  !     Lebedev's approach                                               *
  !     Output:                                                          *
  !     sigma(6,6)  - initial sigma beam matrix                          *
  !                 = [(beta, -alpha) (-alpha, gamma)]                   *
  !                   where beta*gamma-alpha**2 =1                       *
  !----------------------------------------------------------------------*
  double precision, intent (OUT) :: s0mat(6, 6)
  double precision, external :: get_value
  double precision :: e1, e2
  double precision :: betx0, alfx0, bety0, alfy0
  double precision :: beta1x, beta2x, beta1y, beta2y, &
                      alfa1x, alfa2x, alfa1y, alfa2y

  double precision ::gamx0, kappa, gamy0
  double precision ::r11, r12, r21, r22, u, v1, v2, r11new, sumrelement
  double precision, parameter :: eps=1d-36
  betx0 = opt_fun0(3)
  bety0 = opt_fun0(6)
  alfx0 = opt_fun0(4)
  alfy0 = opt_fun0(7)
  r11 =  opt_fun0(29);  r12 =  opt_fun0(30);  r21 =  opt_fun0(31);  r22 =  opt_fun0(32)

  e1 = get_value('probe ','ex ')!BEAM->Ex
  e2 = get_value('probe ','ey ')!BEAM->Ey

  ! We have to check if there is any coupling. If not the coupled formula fails
  ! since we devide by 0. If you go to the limit the formula works but since 0/0 is undefined it doesn't work.

  sumrelement = abs(r11) + abs(r12) + abs(r21) + abs(r22)
  !Uncoupled case
  if ( sumrelement < eps ) then
      s0mat(1, 1) =  e1*betx0
      s0mat(2, 2) =  e1*(1 + alfx0**2)/betx0
      s0mat(1, 2) =  -e1*alfx0
      s0mat(2, 1) = s0mat(1, 2)

      s0mat(3, 3) =  e2*bety0
      s0mat(4, 4) =  e2*(1 + alfy0**2)/bety0
      s0mat(4, 3) =  -e2*alfy0
      s0mat(3, 4) = s0mat(4, 3)
  ! Coupled case
  else

      kappa = one/(one + (r11*r22 - r12*r21))
      u = one - kappa
      gamx0 = (one + alfx0**2) / betx0;  gamy0 = (one + alfy0**2) / bety0


      beta1x = kappa * betx0
      beta2y = kappa * bety0
      beta2x = kappa * ( r22**2*bety0 + two*r12*r22*alfy0 + r12**2*gamy0 )
      beta1y = kappa * ( r11**2*betx0 - two*r12*r11*alfx0 + r12**2*gamx0 )
      alfa1x =  kappa * alfx0
      alfa2y =  kappa * alfy0
      alfa2x =  kappa * ( r21*r22*bety0 + (r12*r21 + r11*r22)*alfy0 + r11*r12*gamy0 )
      alfa1y = -kappa * ( r21*r11*betx0 - (r12*r21 + r11*r22)*alfx0 + r12*r22*gamx0 )

      ! There are two choises of eigenvectors and this is to chose the correspondent one.
      v2 = asin( r12 * (1-u) / sqrt(beta1x*beta1y) )
      r11new = sqrt( beta2y / beta2x ) * ( alfa2x*sin(v2) + u*cos(v2) ) / (1-u)

      if ( abs(r11-r11new) < abs(r11+r11new) ) then
        v2 = v2 + twopi/2
        v1 = asin ( sqrt(beta2x* beta2y)* sin(v2) / sqrt(beta1x*beta1y) ) + twopi/2
      else
        v2 = v2
        v1 = asin ( sqrt(beta2x* beta2y)* sin(v2) / sqrt(beta1x*beta1y) )
      end if

      ! This is the element-by-element of the matrix multiplication V*Sigma(uncoupled)*V^T
      s0mat(1, 1) =  e1 * beta1x + e2 * beta2x
      s0mat(2, 2) =  e1 * ((one-u)**2+alfa1x**2) / beta1x + e2 * (u**2+alfa2x**2) / beta2x
      s0mat(3, 3) =  e1 * beta1y + e2 * beta2y
      s0mat(4, 4) =  e1 * (u**2+alfa1y**2)/beta1y+e2 * ((1-u)**2+alfa2y**2) / beta2y

      s0mat(1, 2) =  -e1*alfa1x - e2*alfa2x
      s0mat(2, 1) = s0mat(1, 2)

      s0mat(1, 3) = e1 * sqrt(beta1x*beta1y) * cos(v1)-e2 * sqrt(beta2x*beta2y)*cos(v2)
      s0mat(3, 1) = s0mat(1, 3)

      s0mat(3, 4) = -e1 * alfa1y - e2 * alfa2y
      s0mat(4, 3) = s0mat(3, 4)

      s0mat(1, 4) =  e1 * sqrt(beta1x / beta1y) * ( u*sin(v1) -alfa1y*cos(v1) ) &
      - e2 * sqrt(beta2x / beta2y)*( (1-u)*sin(v2) -alfa2y*cos(v2) )
      s0mat(4, 1) = s0mat(1, 4)

      s0mat(2, 3) = -e1*sqrt(beta1y / beta1x)*( (1-u)*sin(v1)+ alfa1x*cos(v1) ) &
      -e2*sqrt(beta2y / beta2x) * ( (u)*sin(v2) -alfa2x*cos(v2) )
      s0mat(3, 2) = s0mat(2, 3)

      s0mat(2, 4) = e1 * ( (alfa1y*(1-u)-alfa2x*u)*sin(v1)+(u*(1-u)+alfa1x*alfa1y) * cos(v1) ) / sqrt(beta2x*beta1y) &
      - e2 * ( (alfa2x*(1-u) - alfa2y*u) * sin(v2) + (u*(1-u)+alfa2x*alfa2y)*cos(v2) ) / sqrt(beta2x*beta2y)
      s0mat(4, 2) = s0mat(2, 4)

    endif

END SUBROUTINE tmsigma

SUBROUTINE tmsigma_emit(rt, s0mat)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  use math_constfi, only : zero, twopi
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Convert eigenvectors to internal sigma matrix form.                *
  ! Input:                                                               *
  !   EM(6,6)   (real)    Eigenvector matrix.                            *
  !   EX        (real)    Horizontal emittance.                          *
  !   EY        (real)    Vertical emittance.                            *
  !   ET        (real)    Longitudinal emittance.                        *
  ! Output:                                                              *
  !   SIGMA(6,6)(real)    Beam matrix in internal form.                  *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: rt(6,6)
  double precision :: ex, ey, et, em(6,6), s0mat(6,6), tmp_e(36)
  double precision :: reval(6), aival(6) ! re and im parts
  double precision, external :: get_value
  logical, external :: m66sta
  external :: print_eigenvectors
  logical:: saveig

  integer :: j, k

  !---- Find eigenvectors at initial position.
  if (m66sta(rt)) then
     call laseig(rt, reval, aival, em)
  else
     call ladeig(rt, reval, aival, em)
  endif

  ex = get_value('probe ','ex ')!BEAM->Ex
  ey = get_value('probe ','ey ')!BEAM->Ey
  et = get_value('probe ','et ')!BEAM->Ez

  saveig = get_value('twiss ','eigenvector ').ne.zero
  if(saveig) then
    tmp_e = RESHAPE(em, shape(tmp_e))
    call print_eigenvectors(tmp_e)
  end if

  do j = 1, 6
    do k = 1, 6
      s0mat(j,k) = ex * (em(j,1)*em(k,1) + em(j,2)*em(k,2)) + &
                   ey * (em(j,3)*em(k,3) + em(j,4)*em(k,4))
      !---- Solve for dynamic case
      if (.not.m66sta(rt)) then
        s0mat(j,k) = s0mat(j,k) + et * (em(j,5)*em(k,5) + em(j,6)*em(k,6))
      endif
    enddo
  enddo
END SUBROUTINE tmsigma_emit

SUBROUTINE twcptk_twiss_new(matx, maty, error)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  use math_constfi, only : zero, twopi
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Calculation of twiss parameters by Irina Tecker                  *
  !     Irina Tecker : new version rewritten in matrix formalism         *
  !                   + added mode flip                                  *
  !                   + no "guard" for negative betx                     *
  !     Input:                                                           *
  !     matx(2,2)  (double)   X-plane matrix of block-diagonal           *
  !     maty(2,2)  (double)   Y-plane matrix of block-diagonal           *
  !----------------------------------------------------------------------*

  double precision :: matx(2,2), maty(2,2)
  double precision :: matx11, matx12, matx21, matx22
  double precision :: maty11, maty12, maty21, maty22
  double precision :: alfx_ini, betx_ini, tempa
  double precision :: alfy_ini, bety_ini, tempb
  double precision :: detx, dety
  logical          :: error
  double precision, parameter :: eps=1d-36
  character(len=name_len) :: name

  !---Initialization
  alfx_ini=zero; betx_ini=zero
  alfy_ini=zero; bety_ini=zero
  error = .true.

  detx = matx(1,1) * matx(2,2) - matx(1,2) * matx(2,1)
  dety = maty(1,1) * maty(2,2) - maty(1,2) * maty(2,1)

  call element_name(name,len(name))

  if (detx == 0) return
  if (dety == 0) return
  betx_ini = betx ; alfx_ini = alfx ;
  bety_ini = bety ; alfy_ini = alfy ;


  matx11 = matx(1,1)
  matx12 = matx(1,2)
  matx21 = matx(2,1)
  matx22 = matx(2,2)

  maty11 = maty(1,1)
  maty12 = maty(1,2)
  maty21 = maty(2,1)
  maty22 = maty(2,2)

  !---- Mode 1.
  tempb = matx11 * betx_ini - matx12 * alfx_ini
  tempa = matx21 * betx_ini - matx22 * alfx_ini
  alfx = - (tempa * tempb + matx12 * matx22) / (betx_ini)
  betx =   (tempb * tempb + matx12 * matx12) / (betx_ini)
  if (abs(matx12).gt.eps) amux = amux + atan2(matx12,tempb)

  !---- Mode 2.
  tempb = maty11 * bety_ini - maty12 * alfy_ini
  tempa = maty21 * bety_ini - maty22 * alfy_ini
  alfy = - (tempa * tempb + maty12 * maty22) / (bety_ini)
  bety =   (tempb * tempb + maty12 * maty12) / (bety_ini)
  if (abs(maty12).gt.eps) amuy = amuy + atan2(maty12,tempb)

  error = .false.
end SUBROUTINE twcptk_twiss_new

SUBROUTINE twcptk_sagan(re,orbit) ! new, RD matrix, talman, sagan
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  use matrices, only : JMAT, JMATT, SMAT, SMATT
  use math_constfi, only : zero, one, two, twopi
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track coupled lattice functions.                                 *
  !     Irina Tecker : new version rewritten in matrix formalism         *
  !                   + Talman/Sagan                                     *
  !                   + added mode flip                                  *
  !                   + no "guard" for negative betx                     *
  !     Input:                                                           *
  !     re(6,6)  (double)   transfer matrix of element.                  *
  !     orbit(6) (double)   closed orbit                                 *
  !----------------------------------------------------------------------*
  double precision :: re(6,6), orbit(6)

  integer :: i, i1, i2, j, irank
  double precision :: rwi(6,6), rc(6,6), dt(6)
  double precision :: a(2,2), b(2,2), c(2,2), d(2,2),ra(4,8)
  double precision :: rmat0(2,2), e(2,2), f(2,2), cd(2,2), tmp(2,2)
  double precision :: rmat_bar(2,2), ebar(2,2), fbar(2,2), dbar(2,2)
  double precision :: edet, fdet, tempa, tempb, det, det44, gamma_2
  double precision :: alfx0, betx0, amux0, alfx_ini, betx_ini
  double precision :: alfy0, bety0, amuy0, alfy_ini, bety_ini

  character(len=name_len) :: name
  character(len=150)      :: warnstr

  integer, external :: get_option
  double precision, parameter :: eps=1d-36
  logical :: mode_flip_ele =.false.
  logical :: cp_error=.false.

  !---Initialization
  alfx0=zero; betx0=zero; amux0=zero; alfx_ini=zero; betx_ini=zero
  alfy0=zero; bety0=zero; amuy0=zero; alfy_ini=zero; bety_ini=zero

  call element_name(name,len(name))
  !---- Dispersion.
  DT = matmul(RE, DISP)

  ! RE(1:4,1:4) = ( A , B
  !                 C , D)
  A = RE(1:2,1:2) ; B = RE(1:2,3:4)
  C = RE(3:4,1:2) ; D = RE(3:4,3:4)

  !---- Check RE rank
  RA(1:4,1:4) = RE(1:4,1:4)
  RA(1:4,5:8)  = zero
  RA(1,5)     = one
  RA(2,6)     = one
  RA(3,7)     = one
  RA(4,8)     = one

  call solver(ra,4,4,irank)
  if (irank.lt.4) then
     write (warnstr, '(a)') 'Coupling failed: Rank deficient matrix RE, irank < 4 '
     call fort_warn('TWCPTK: ', warnstr)
     return
  endif

  if   (all( B == zero  .and.  C == zero)  ) then

     E(:2,:2) = A(:2,:2)
     F(:2,:2) = D(:2,:2)
     DBAR = matmul(SMAT, matmul(transpose(D),SMATT)) ! symplectic conjugate of D = S*D^T*S^T
     RMAT = matmul(matmul(E, RMAT), DBAR)

  else

     RMAT_BAR  = matmul(SMAT, matmul(transpose(RMAT),SMATT)) ! symplectic conjugate of RMAT
     TMP       = gammacp*A -  matmul(B,RMAT_BAR)
     det       = (tmp(1,1) * tmp(2,2) - tmp(1,2) * tmp(2,1))

     if (det .gt. 0.9 .or. (det .gt. 0.1 .and. .not. mode_flip) ) then  !--- normal mode

        gamma_2 = sqrt (det)
        E = TMP / gamma_2
        F  = ( gammacp*D + matmul(C, RMAT)) / gamma_2

        CD   = matmul(A, RMAT) + gammacp*B
        FBAR = matmul(SMAT, matmul(transpose(F),SMATT)) ! symplectic conjugate of F = S*F^T*S^T
        RMAT = matmul(CD, FBAR)

        gammacp = gamma_2

     else  !---flip mode

        mode_flip_ele = .not. mode_flip
        nmode_flip = nmode_flip + 1

        TMP = matmul(A,RMAT) + gammacp*B
        det = (tmp(1,1) * tmp(2,2) - tmp(1,2) * tmp(2,1))

        if (det .lt. zero) then
           write (warnstr, '(a, a, a, f12.6)') 'Negative determinant in ', name, ', det = ', det
           call fort_warn('TWCPTK: ', warnstr)
           !print*, "Negative deterinant! Rounding error?!?!?"
        endif

        gamma_2 = sqrt (abs(det))
        E  = ( gammacp*C - matmul(D, RMAT_BAR)) / gamma_2
        F = TMP / gamma_2

        CD = gammacp*A - matmul(B, RMAT_BAR)
        EBAR = matmul(SMAT, matmul(transpose(E),SMATT)) ! symplectic conjugate of E = S*E^T*S^T
        RMAT = matmul(CD, EBAR)

        gammacp = gamma_2

        write (warnstr, '(a, a, a, i3)') 'Mode flip in the element ', name, ', nflips up to now = ', nmode_flip
        call fort_warn('TWCPTK: ', warnstr)
        !print *, '+++ mode flip in the element ', name, ", nflips up to now = ", nmode_flip
     endif
  endif


  if(mode_flip) then
     call twcptk_twiss_new(f, e, cp_error)
  else
     call twcptk_twiss_new(e, f, cp_error)
  endif

  if (cp_error) then
     ! print *, '+++ det of block diagonal matrix is zero in ', name
     write (warnstr, '(a, a, a)') 'Det of block diagonal matrix is zero in ', name, ', twiss parameter might be unphysical! '
     call fort_warn('TWCPTK: ', warnstr)
  endif

  ! When we are comming out of a flipped mode, the phase is often off by a factor of twopi. As twopi doesn't have any effect on physics we can substract twopi.
  ! Unfortunately there is no definitive way to calcuate what the "right" way to handle this is but "-twopi" is better than nothing (c, Sagan)
  if(mode_flip .and. .not. mode_flip_ele) then

     amux = amux - twopi
     amuy = amuy - twopi

  endif

  mode_flip = mode_flip_ele

  EDET = E(1,1) * E(2,2) - E(1,2) * E(2,1)           ! former adet
  FDET = F(1,1) * F(2,2) - F(1,2) * F(2,1)

  if ( EDET .gt. two .or. FDET .gt. two .or. EDET .lt. zero .or. FDET .lt. zero) then
     ! print *, '+++ det of block diagonal matrix is > 2 in ', name
     write (warnstr, '(a, a)') 'Det of block diagonal matrix is > 2 or  < 0 in ', name
     call fort_warn('TWCPTK: ', warnstr)
  endif

  if (betx.lt.eps .or. bety.lt.eps) then
     write (warnstr, '(a, a, a)') 'Negative beta in ', name, ' Twiss parameter might be unphysical!'
     call fort_warn('TWCPTK: ', warnstr)
     print *, '+++ negative beta: name', name, 'betx=', betx, ', bety=', bety
  endif

  ! !---- Auxiliary matrices.
  ! !LD:  a = re_x - re_xy*rmat
  ! a(1,1) = re(1,1) - (re(1,3) * rmat(1,1) + re(1,4) * rmat(2,1))
  ! a(1,2) = re(1,2) - (re(1,3) * rmat(1,2) + re(1,4) * rmat(2,2))
  ! a(2,1) = re(2,1) - (re(2,3) * rmat(1,1) + re(2,4) * rmat(2,1))
  ! a(2,2) = re(2,2) - (re(2,3) * rmat(1,2) + re(2,4) * rmat(2,2))

  ! b(1,1) = re(3,1) - (re(3,3) * rmat(1,1) + re(3,4) * rmat(2,1))
  ! b(1,2) = re(3,2) - (re(3,3) * rmat(1,2) + re(3,4) * rmat(2,2))
  ! b(2,1) = re(4,1) - (re(4,3) * rmat(1,1) + re(4,4) * rmat(2,1))
  ! b(2,2) = re(4,2) - (re(4,3) * rmat(1,2) + re(4,4) * rmat(2,2))

  ! c(1,1) = re(3,3) + (re(3,1) * rmat(2,2) - re(3,2) * rmat(2,1))
  ! c(1,2) = re(3,4) - (re(3,1) * rmat(1,2) - re(3,2) * rmat(1,1))
  ! c(2,1) = re(4,3) + (re(4,1) * rmat(2,2) - re(4,2) * rmat(2,1))
  ! c(2,2) = re(4,4) - (re(4,1) * rmat(1,2) - re(4,2) * rmat(1,1))


  ! !---- Track R matrix.
  ! adet = a(1,1) * a(2,2) - a(1,2) * a(2,1)
  ! if (abs(adet) .gt. eps) then
  !    rmat(1,1) = - (b(1,1) * a(2,2) - b(1,2) * a(2,1)) / adet
  !    rmat(1,2) =   (b(1,1) * a(1,2) - b(1,2) * a(1,1)) / adet
  !    rmat(2,1) = - (b(2,1) * a(2,2) - b(2,2) * a(2,1)) / adet
  !    rmat(2,2) =   (b(2,1) * a(1,2) - b(2,2) * a(1,1)) / adet

  !    !---- Mode 1.
  !    tempb = a(1,1) * betx - a(1,2) * alfx
  !    tempa = a(2,1) * betx - a(2,2) * alfx
  !    alfx = - (tempa * tempb + a(1,2) * a(2,2)) / (adet * betx)
  !    betx =   (tempb * tempb + a(1,2) * a(1,2)) / (adet * betx)
  !    if (abs(a(1,2)) .gt. eps) amux = amux + atan2(a(1,2),tempb)

  !    !---- Mode 2.
  !    tempb = c(1,1) * bety - c(1,2) * alfy
  !    tempa = c(2,1) * bety - c(2,2) * alfy
  !    alfy = - (tempa * tempb + c(1,2) * c(2,2)) / (adet * bety)
  !    bety =   (tempb * tempb + c(1,2) * c(1,2)) / (adet * bety)
  !    if (abs(c(1,2)) .gt. eps) amuy = amuy + atan2(c(1,2),tempb)
  ! else
  !   !LD: 09.2015
  !   call element_name(name,len(name))
  !   print *, 'coupling too strong for element ', name, ' (adet= ', adet, ')'
  !   call fort_warn('TWCPTK: ','twiss parameter might be unphysical (split element)')
  ! endif



  !---- Cumulative R matrix and one-turn map at element location.
  if (rmatrix) then
     RW = matmul(RE,RW)
     if (get_option('twiss_inval ') .ne. 0) then
        RC = RW
     else
        RWI = matmul(JMATT, matmul(transpose(RW),JMAT)) ! invert symplectic matrix
        RC = matmul(RW,matmul(ROTM,RWI))
     endif
  endif

  if (rmatrix) then
     do i1=1,6
        do i2=1,6
           opt_fun(33 + (i1-1)*6 + i2) = rc(i1,i2)
        enddo
     enddo
  endif

end SUBROUTINE twcptk_sagan


SUBROUTINE twbtin(rt,tt)
  use twiss0fi
  use twisscfi
  use math_constfi, only : zero, one, two, quarter
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values of Chromatic functions.                           *
  !     Input:                                                           *
  !     rt(6,6)   (double)  transfer matrix.                             *
  !     tt(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6), tt(6,6,6)

  logical :: stabx, staby
  integer :: i, j, k
  double precision :: disp0(6), ddisp0(6), rtp(6,6), aux(6)
  double precision :: sinmu2, bx, ax, by, ay, temp

  integer, external :: get_option
  double precision, parameter :: eps=1d-8

  !---- Initialization
  betx   = opt_fun0(3 )
  alfx   = opt_fun0(4 )
  amux   = opt_fun0(5 )
  bety   = opt_fun0(6 )
  alfy   = opt_fun0(7 )
  amuy   = opt_fun0(8 )
  wx     = opt_fun0(19)
  phix   = opt_fun0(20)
  dmux   = opt_fun0(21)
  wy     = opt_fun0(22)
  phiy   = opt_fun0(23)
  dmuy   = opt_fun0(24)

  !---- Initial value flag.
  if (get_option('twiss_inval ') .ne. 0) then
     DISP(:4)  = OPT_FUN0(15:18)
     DDISP(:4) = OPT_FUN0(25:28)
     disp(5)  = zero ; disp(6)  = one
     ddisp(5) = zero ; ddisp(6) = zero
     return
  endif

  !---- Initial dispersion.
  call twdisp(rt,rt(:,6),disp0)
  disp0(5) = zero ;  disp0(6) = one

  !---- Derivative of transfer matrix w.r.t. delta(p)/p.
  AUX = zero
  do i = 1, 6
     do k = 1, 6
        temp = dot_product(TT(i,:,k),DISP0)
        aux(i) = aux(i) + temp * disp0(k)
        rtp(i,k) = two * temp
     enddo
  enddo

  !---- Derivative of dispersion.
  call twdisp(rt,aux,ddisp0)

  ddisp0(5) = zero
  ddisp0(6) = zero
  DISP  = DISP0
  DDISP = DDISP0

  !---- Horizontal motion.
  cosmux = (rt(1,1) + rt(2,2)) / two
  stabx = abs(cosmux) .lt. one

  if (stabx) then
     sinmu2 = - rt(1,2)*rt(2,1) - quarter*(rt(1,1) - rt(2,2))**2
     if (sinmu2.lt.0) sinmu2 = eps
     sinmux = sign(sqrt(sinmu2), rt(1,2))
     betx = rt(1,2) / sinmux
     alfx = (rt(1,1) - rt(2,2)) / (two * sinmux)
     bx = rtp(1,2) / rt(1,2) + (rtp(1,1) + rtp(2,2)) * cosmux / (two * sinmu2)
     ax = (rtp(1,1) - rtp(2,2)) / (two * sinmux) - alfx * rtp(1,2) / rt(1,2)
     wx = sqrt(bx**2 + ax**2)
     if (wx.gt.eps) phix = atan2(ax,bx)
  endif

  !---- Vertical motion.
  cosmuy = (rt(3,3) + rt(4,4)) / two
  staby = abs(cosmuy) .lt. one

  if (staby) then
     sinmu2 = - rt(3,4)*rt(4,3) - quarter*(rt(3,3) - rt(4,4))**2
     if (sinmu2.lt.0) sinmu2 = eps
     sinmuy = sign(sqrt(sinmu2), rt(3,4))
     bety = rt(3,4) / sinmuy
     alfy = (rt(3,3) - rt(4,4)) / (two * sinmuy)
     by = rtp(3,4) / rt(3,4) + (rtp(3,3) + rtp(4,4)) * cosmuy / (two * sinmu2)
     ay = (rtp(3,3) - rtp(4,4)) / (two * sinmuy) - alfy * rtp(3,4) / rt(3,4)
     wy = sqrt(by**2 + ay**2)
     if (wy.gt.eps) phiy = atan2(ay,by)
  endif

  !---- Fill optics function array
  opt_fun0(19) = wx
  opt_fun0(20) = phix
  opt_fun0(22) = wy
  opt_fun0(23) = phiy
  opt_fun0(25) = ddisp(1)
  opt_fun0(26) = ddisp(2)
  opt_fun0(27) = ddisp(3)
  opt_fun0(28) = ddisp(4)

end SUBROUTINE twbtin

SUBROUTINE twchgo
  use twiss0fi
  use twisslfi
  use twissafi
  use twisscfi
  use twissbeamfi, only : radiate, deltap, beta, gamma
  use spch_bbfi, only : i_spch
  use math_constfi, only : zero, one, two
  use code_constfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track Chromatic functions.                                       *
  !----------------------------------------------------------------------*
  logical :: fmap, cplxy, cplxt
  integer :: i, code, save, n_align
  double precision :: orbit(6), orbit2(6), ek(6), re(6,6), te(6,6,6)
  double precision :: orbit00(6), ek00(6), re00(6,6), te00(6,6,6), disp00(6), ddisp00(6)
  double precision :: rmat0(2,2), sigmat00(6,6)
  double precision :: al_errors(align_max), el, dl
  character(len=130) :: msg
  double precision :: betx0, alfx0, amux0, wx0, dmux0, phix0
  double precision :: bety0, alfy0, amuy0, wy0, dmuy0, phiy0
  integer, external :: restart_sequ, advance_node, get_option, node_al_errors
  integer, external :: start_interp_node, fetch_interp_node
  double precision, external :: node_value, get_value

  !---- If save requested reset table
  save = get_option('twiss_save ')
  if (save .ne. 0) call reset_count(table_name)

  centre = get_value('twiss ','centre ').ne.zero

  !---- Initial values for lattice functions.
  amux = zero
  amuy = zero
  ORBIT = OPT_FUN0(9:14)
  DISP(1:4) = OPT_FUN0(15:18)
  disp(5) = zero
  disp(6) = one
  TE      = zero
  cplxy = .false.; cplxt = .false.

  !---- Initial values for chromatic functions.
  opt_fun(19) = wx
  opt_fun(20) = phix
  opt_fun(21) = dmux
  opt_fun(22) = wy
  opt_fun(23) = phiy
  opt_fun(24) = dmuy

  OPT_FUN(25:28) = DDISP(1:4)

  !---- and synchrotron radiation integrals
  synch_1 = zero; synch_2 = zero; synch_3 = zero; synch_4 = zero; synch_5 = zero
  synch_6=zero;  synch_8=zero

  !---- Loop over positions.
  i = restart_sequ()
  i_spch=0

  i = 1
  do while (i .ne. 0)
    el = node_value('l ')
    if (start_interp_node(i) .ne. 0) then
      do while (fetch_interp_node(i, dl) .ne. 0)
        call backup_optics()
        call track_one_element(dl, .true.)
        call restore_optics()
      end do
      call track_one_element(el, .false.)
    else
      i = 1
      call track_one_element(el, .not. centre)
    endif
    i = advance_node()
  end do

contains

subroutine track_one_element(el, fexit)
  double precision, intent(in) :: el
  logical :: fexit

  code = node_value('mad8_type ')
!  if (code .eq. code_tkicker)     code = code_kicker ! TKICKER is a KICKER
  if (code .eq. code_placeholder) code = code_instrument ! PLACEHOLDER is an INSTRUMENT

  !---- Physical element.
  n_align = node_al_errors(al_errors)
  if (n_align .ne. 0)  then
     ORBIT2 = ORBIT
     call tmali1(orbit2,al_errors,beta,gamma,orbit,re)
     call tw_synch_int()
     call twbttk(re,te)
  endif

  if (centre) then
     call backup_optics()

     call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te,.true.,el/two)

     if (fmap) call twbttk(re,te)

     call save_opt_fun()
     call twprep(save,2,opt_fun,zero,i)

     call restore_optics()
  endif

  call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te,.false.,el)

  if (fmap) then
      call tw_synch_int()
      call twbttk(re,te)
  endif

  if (n_align.ne.0)  then
     ORBIT2 = ORBIT
     call tmali2(el,orbit2,al_errors,beta,gamma,orbit,re)
     call tw_synch_int()
     call twbttk(re,te)
  endif

  call save_opt_fun()
  if (.not.centre) then
     call twprep(save,2,opt_fun,zero,i)
  endif
end subroutine track_one_element

subroutine check_summary()
  !---- Warning, if system is coupled.
  if (cplxy) then
     write (msg,'(a,f12.6,a)') "TWISS found transverse coupling for delta(p)/p =",deltap, &
          "chromatic functions may be wrong."
     call fort_warn('TWCHGO: ',msg)
  endif
  if (cplxt .or. radiate) then
     call fort_warn('TWCHGO: ','TWISS uses the RF system or synchrotron radiation '// &
                 'only to find the closed orbit, for optical calculations it '// &
                 'ignores both.')
  endif
end subroutine check_summary

subroutine backup_optics()
     ORBIT00 = ORBIT ; EK00 = EK ; RE00 = RE ; TE00 = TE
     betx0=betx; alfx0=alfx; amux0=amux; wx0=wx; dmux0=dmux; phix0=phix
     bety0=bety; alfy0=alfy; amuy0=amuy; wy0=wy; dmuy0=dmuy; phiy0=phiy
     RMAT0 = RMAT ; disp00 = disp ; ddisp00 = ddisp
     sigmat00 = sigmat
end subroutine backup_optics

subroutine restore_optics()
     betx=betx0; alfx=alfx0; amux=amux0; wx=wx0; dmux=dmux0; phix=phix0
     bety=bety0; alfy=alfy0; amuy=amuy0; wy=wy0; dmuy=dmuy0; phiy=phiy0
     RMAT = RMAT0 ; disp = disp00 ; ddisp = ddisp00
     ORBIT = ORBIT00 ; EK = EK00 ; RE = RE00 ; TE = TE00
     sigmat = sigmat00
end subroutine restore_optics

subroutine save_opt_fun()
     opt_fun(19) = wx
     opt_fun(20) = phix
     opt_fun(21) = dmux
     opt_fun(22) = wy
     opt_fun(23) = phiy
     opt_fun(24) = dmuy

     OPT_FUN(25:28) = ddisp(1:4)
end subroutine

end SUBROUTINE twchgo

SUBROUTINE tw_synch_int()
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissbeamfi, only : beta
  use math_constfi, only : zero, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Compute and add synchrotron radiation integral                   *
  !----------------------------------------------------------------------*
  double precision :: an, e1, e2, sk1, rhoinv, blen
  double precision :: syncint(8)

  double precision, external :: node_value

  !---- Initialisation
  blen = node_value('blen ')
  rhoinv = node_value('rhoinv ')
  sk1 = node_value('k1 ') + node_value('k1tap ')
  e1 = node_value('e1 ')
  e2 = node_value('e2 ')
  an = node_value('angle ')
  if (node_value('mad8_type ') .eq. 2) then ! RBEND
     e1 = e1 + an / two
     e2 = e2 + an / two
  endif

  !---- Synchrotron radiation integrals through bending magnets.
  
     ! Note that calcsyncint expects dx and dpx as derivatives wrt deltap.
     ! since MAD take disp(1) and disp(2) as derivatives wrt pt, they must be
     ! multiplied by beta before the call to calcsyncint.
     syncint = zero
     call calcsyncint(rhoinv,blen,sk1,e1,e2,betx,alfx,disp(1)*beta,disp(2)*beta,syncint)
     synch_1 = synch_1 + syncint(1)
     synch_2 = synch_2 + syncint(2)
     synch_3 = synch_3 + syncint(3)
     synch_4 = synch_4 + syncint(4)
     synch_5 = synch_5 + syncint(5)
     synch_6 = synch_6 + syncint(6)
     synch_8 = synch_8 + syncint(8)

end subroutine tw_synch_int


SUBROUTINE twbttk(re,te)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissbeamfi, only : beta
  use math_constfi, only : zero, one, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track lattice functions, including chromatic effects.            *
  !     Input:                                                           *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  double precision :: re(6,6), te(6,6,6)

  integer :: i,k
  double precision :: auxp(6), rep(6,6), fre(6,6), frep(6,6)
  double precision :: ax1, ax2, ay1, ay2, bx1, bx2, by1, by2
  double precision :: t2, ta, tb, temp, tg
  double precision :: detl, f

  double precision, external :: proxim
  double precision, parameter :: eps=1d-8

  !---- Initialisation
  AUXP = zero
  do i = 1, 6
     do k = 1, 6
        temp = dot_product(TE(i,:,k),DISP(:))
        auxp(i) = auxp(i) + temp*disp(k) + re(i,k)*ddisp(k)
        rep(i,k) = two*temp
     enddo
  enddo

  DISP = matmul(re, disp)
  DDISP = AUXP

  !---- Tor: modified to cancel energy change
  disp(6) = one

  !---- Tor/MDW: scale by square root of the determinant of the
  !     longitudinal 2x2 part of the R-matrix
  detl = re(5,5)*re(6,6) - re(5,6)*re(6,5)
  f = one / sqrt(detl)
  FRE  = f * RE
  FREP = f * REP

  !---- Track horizontal functions including energy scaling.
  tb = fre(1,1)*betx - fre(1,2)*alfx
  ta = fre(2,1)*betx - fre(2,2)*alfx
  t2 = tb**2 + fre(1,2)**2
  tg = fre(1,1)*alfx - fre(1,2)*(one + alfx**2) / betx

  !---- Linear functions.
  alfx = - (tb*ta + fre(1,2)*fre(2,2)) / betx
  betx = t2 / betx
  if (fre(1,2).ne.zero .or. tb.ne.zero) amux = amux + atan2(fre(1,2),tb)
  bx1 = wx*cos(phix)
  ax1 = wx*sin(phix)
  bx2 = ((tb**2 - fre(1,2)**2)*bx1                                  &
       - two*tb*fre(1,2)*ax1) / t2                                  &
       + two*(tb*frep(1,1) - tg*frep(1,2)) / betx
  ax2 = ((tb**2 - fre(1,2)**2)*ax1                                  &
       + two*tb*fre(1,2)*bx1) / t2                                  &
       - (tb*(frep(1,1)*alfx + frep(2,1)*betx)                      &
       - tg*(frep(1,2)*alfx + frep(2,2)*betx)                       &
       + fre(1,1)*frep(1,2) - fre(1,2)*frep(1,1)) / betx
  wx = sqrt(ax2**2 + bx2**2)
  if (wx.gt.eps) phix = proxim(atan2(ax2, bx2), phix)
  dmux = dmux + fre(1,2)*(fre(1,2)*ax1 - tb*bx1) / t2               &
       + (fre(1,1)*frep(1,2) - fre(1,2)*frep(1,1)) / betx

  !---- Track vertical functions including energy scaling.
  tb = fre(3,3)*bety - fre(3,4)*alfy
  ta = fre(4,3)*bety - fre(4,4)*alfy
  t2 = tb**2 + fre(3,4)**2
  tg = fre(3,3)*alfy - fre(3,4)*(one + alfy**2) / bety

  !---- Linear functions.
  alfy = - (tb*ta + fre(3,4)*fre(4,4)) / bety
  bety = t2 / bety
  if (fre(3,4).ne.zero .or. tb.ne.zero) amuy = amuy + atan2(fre(3,4),tb)

  by1 = wy*cos(phiy)
  ay1 = wy*sin(phiy)
  by2 = ((tb**2 - fre(3,4)**2)*by1                                  &
       - two*tb*fre(3,4)*ay1) / t2                                  &
       + two*(tb*frep(3,3) - tg*frep(3,4)) / bety
  ay2 = ((tb**2 - fre(3,4)**2)*ay1                                  &
       + two*tb*fre(3,4)*by1) / t2                                  &
       - (tb*(frep(3,3)*alfy + frep(4,3)*bety)                      &
       - tg*(frep(3,4)*alfy + frep(4,4)*bety)                       &
       + fre(3,3)*frep(3,4) - fre(3,4)*frep(3,3)) / bety
  wy = sqrt(ay2**2 + by2**2)
  if (wy.gt.eps) phiy = proxim(atan2(ay2, by2), phiy)
  dmuy = dmuy + fre(3,4)*(fre(3,4)*ay1 - tb*by1) / t2               &
       + (fre(3,3)*frep(3,4) - fre(3,4)*frep(3,3)) / bety

end SUBROUTINE twbttk

SUBROUTINE tw_summ(rt,tt)
  use twiss0fi
  use twisscfi
  use twissbeamfi, only : deltap, beta, gamma
  use math_constfi, only : zero, one, two, twopi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Compute summary data for TWISS and OPTICS commands.              *
  !     Input:                                                           *
  !     rt(6,6)   (double)  one turn transfer matrix.                    *
  !     tt(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6), tt(6,6,6)

  integer :: i
  double precision :: sd, detl, f, tb, t2
  double precision :: disp0(6), frt(6,6), frtp(6,6), rtp(6,6)
  double precision :: bx0, ax0, by0, ay0, sx, sy, orbit5
  double precision, parameter :: eps=1d-16, diff_cos=5d-5
  character(len=150) :: warnstr

  integer, external :: get_option

  !---- Initialization chromatic part
  RTP = zero
  FRT = zero
  FRTP = zero

  DISP0(1:4) = OPT_FUN0(15:18)

  wx      = opt_fun(19)
  phix    = opt_fun(20)
  dmux    = opt_fun(21)
  wy      = opt_fun(22)
  phiy    = opt_fun(23)
  dmuy    = opt_fun(24)

  DDISP(1:4) = OPT_FUN(25:28)

  !---- Summary data for non-periodic case.
  if (get_option('twiss_inval ') .ne. 0) then
     detl = rt(5,5) * rt(6,6) - rt(5,6) * rt(6,5)
     f = one / sqrt(detl)
     FRT  = f * RT
     FRTP = f * RTP
     tb = frt(1,1) * betx - frt(1,2) * alfx
     t2 = tb**2 + frt(1,2)**2
     bx0 = wx * cos(phix)
     ax0 = wx * sin(phix)
     xix = dmux + frt(1,2) * (frt(1,2) * ax0 - tb * bx0) / t2        &
          + (frt(1,1) * frtp(1,2) - frt(1,2) * frtp(1,1)) / betx
     xix = xix / twopi
     tb = frt(3,3) * bety - frt(3,4) * alfy
     t2 = tb**2 + frt(3,4)**2
     by0 = wy * cos(phiy)
     ay0 = wy * sin(phiy)
     xiy = dmuy + frt(3,4) * (frt(3,4) * ay0 - tb * by0) / t2        &
          + (frt(3,3) * frtp(3,4) - frt(3,4) * frtp(3,3)) / bety
     xiy = xiy / twopi
     alfa  = zero
     gamtr = zero
     cosmux = zero
     cosmuy = zero

     !---- Summary data for periodic case.
  else
     sd = rt(5,6)
     sx = tt(1,1,6) + tt(2,2,6)
     sy = tt(3,3,6) + tt(4,4,6)

     do i = 1, 4
        sd = sd + rt(5,i) * disp(i)
        sx = sx + (tt(1,1,i) + tt(2,2,i)) * disp0(i)
        sy = sy + (tt(3,3,i) + tt(4,4,i)) * disp0(i)
     enddo

     xix = - sx / (twopi * sinmux)
     xiy = - sy / (twopi * sinmuy)
     eta = - sd * beta**2 / suml

     alfa = one / gamma**2 + eta

     if (abs(alfa) .lt. eps) then
        alfa  = zero
        gamtr = zero
     else
        gamtr = sign(one,alfa) * sqrt( one / abs(alfa))
     endif

     if (get_option('info  ') .ne. 0) then
        if (abs(cosmux - cos(amux)) .gt. diff_cos) then
           write (warnstr,'(a,e13.6)') "Difference in the calculation of cosmux: cosmux - cos(amux) =  ", cosmux - cos(amux)
           call fort_warn('TW_SUMM: ', warnstr)
           write (warnstr,'(a,e13.6,a,e13.6)') "cosmux  =  ", cosmux, ", cos(amux) = ", cos(amux)
           call fort_warn('TW_SUMM: ', warnstr)
        endif
        if ( abs(cosmuy - cos(amuy)) .gt. diff_cos) then
           write (warnstr,'(a,e13.6)') "Difference in the calculation of cosmuy: cosmuy - cos(amuy) =  ", cosmuy - cos(amuy)
           call fort_warn('TW_SUMM: ', warnstr)
           write (warnstr,'(a,e13.6,a,e13.6)') "cosmuy  =  ", cosmuy, ", cos(amuy) = ", cos(amuy)
           call fort_warn('TW_SUMM: ', warnstr)
        endif
     endif

  endif

  !---- Initialization transverse
  !---  fix length problem - HG 14.4.08 ! ghislain : ???
  suml    = opt_fun (2)
  betx    = opt_fun0(3)
  alfx    = opt_fun0(4)

  amux    = opt_fun (5)

  bety    = opt_fun0(6)
  alfy    = opt_fun0(7)

  amuy    = opt_fun (8)

  qx = amux / twopi
  qy = amuy / twopi

  !---- Adjust values
  orbit5 = -opt_fun0(13)

  !---- check rmat
  ! if (abs(opt_fun0(29)).gt.1d-8 .or. abs(opt_fun0(30)).gt.1d-8 .or. &
  !     abs(opt_fun0(31)).gt.1d-8 .or. abs(opt_fun0(32)).gt.1d-8) then
  !     print *, 'rmat=', opt_fun0(29:32)
  !     call fort_warn('Chromaticity calculation wrong due to coupling, ',&
  !                    'use chrom option or manual calculation')
  ! endif

  !---- Fill summary table
  call double_to_table_curr('summ ','length ' ,suml)
  call double_to_table_curr('summ ','orbit5 ' ,orbit5)
  call double_to_table_curr('summ ','alfa '   ,alfa)
  call double_to_table_curr('summ ','gammatr ',gamtr)
  call double_to_table_curr('summ ','q1 '     ,qx)
  call double_to_table_curr('summ ','dq1 '    ,xix)
  call double_to_table_curr('summ ','betxmax ',bxmax)
  call double_to_table_curr('summ ','dxmax '  ,dxmax)
  call double_to_table_curr('summ ','dxrms '  ,sigdx)
  call double_to_table_curr('summ ','xcomax ' ,xcomax)
  call double_to_table_curr('summ ','xcorms ' ,sigxco)
  call double_to_table_curr('summ ','q2 '     ,qy)
  call double_to_table_curr('summ ','dq2 '    ,xiy)
  call double_to_table_curr('summ ','betymax ',bymax)
  call double_to_table_curr('summ ','dymax '  ,dymax)
  call double_to_table_curr('summ ','dyrms '  ,sigdy)
  call double_to_table_curr('summ ','ycomax ' ,ycomax)
  call double_to_table_curr('summ ','ycorms ' ,sigyco)
  call double_to_table_curr('summ ','deltap ' ,deltap)
  call double_to_table_curr('summ ','synch_1 ' ,synch_1)
  call double_to_table_curr('summ ','synch_2 ' ,synch_2)
  call double_to_table_curr('summ ','synch_3 ' ,synch_3)
  call double_to_table_curr('summ ','synch_4 ' ,synch_4)
  call double_to_table_curr('summ ','synch_5 ' ,synch_5)
  call double_to_table_curr('summ ','synch_6 ' ,synch_6)
  call double_to_table_curr('summ ','synch_8 ' ,synch_8)

end SUBROUTINE tw_summ


SUBROUTINE tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te,fcentre,dl)
  use twtrrfi
  use name_lenfi
  use time_varfi
  use matrices, only: EYE
  use math_constfi, only : zero
  use code_constfi
  use BeamBeam
  implicit none
  !----------------------------------------------------------------------*
  !     purpose:                                                         *
  !     transport map for a complete element.                            *
  !     optionally, follow orbit.                                        *
  !     input:                                                           *
  !     code                element type code                            *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     dl        (double)  slice length.                                *
  !     input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     bval                basic element parameters                     *
  !     fval                element force parameters                     *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  integer :: code
  logical :: fsec, ftrk, fmap, fcentre
  double precision :: dl
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: plot_tilt, el
  double precision :: node_value

  integer, external :: get_option

  !---- Initialization
  EK = zero
  RE = EYE
  TE = zero
  plot_tilt = zero
  fmap = .false.
  time_var_p = .false.

  el = node_value('l ')

  !---- Select element type.
  select case (code)

     case (code_drift, code_hmonitor:code_rcollimator, code_instrument, code_twcavity, &
        code_slmonitor:code_imonitor, code_placeholder, code_collimator)
        !---- Drift space, monitors and derivatives, collimators, instrument
        call tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)

     case (code_rbend, code_sbend)
        call tmbend(ftrk,fcentre,orbit,fmap,el,dl,ek,re,te, code)

     case (code_matrix)
        call tmarb(fsec,ftrk,orbit,fmap,ek,re,te)

     case (code_quadrupole)
        call tmquad(fsec,ftrk,fcentre,plot_tilt,orbit,fmap,el,dl,ek,re,te)

     case (code_sextupole)
        call tmsext(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)

     case (code_octupole)
        call tmoct(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)

     case (code_multipole)
        if(get_option('thin_cf ').ne.zero .and. node_value('lrad ') .gt. zero) then
          !call tmmult_cf_short(fsec,ftrk,orbit,fmap,re,te)
          call tmmult_cf(fsec,ftrk,orbit,fmap,re,te)
        else
          call tmmult(fsec,ftrk,orbit,fmap,re,te)
        endif
     case (code_solenoid)
        call tmsol(fsec,ftrk,orbit,fmap,dl,ek,re,te)

     case (code_rfcavity)
        call tmrf(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)

     case (code_elseparator)
        call tmsep(fsec,ftrk,fcentre,orbit,fmap,dl,ek,re,te)

     case (code_srotation)
        call tmsrot(ftrk,orbit,fmap,ek,re,te)

     case (code_yrotation)
        call tmyrot(ftrk,orbit,fmap,ek,re,te)

     case (code_xrotation)
        call tmxrot(ftrk,orbit,fmap,ek,re,te)

     case (code_hkicker, code_vkicker, code_kicker, code_tkicker)
        call tmcorr(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)

     case (code_beambeam)
        !---- (Particles/bunch taken for the opposite beam).
        call tmbb(fsec,ftrk,orbit,fmap,re,te)

     case (code_marker)

        ! nothing on purpose!

     case (code_wire)
        ! nothing for now...

     case (code_dipedge)
        call tmdpdg(ftrk,orbit,fmap,ek,re,te)

     case (code_translation)
        call tmtrans(fsec,ftrk,orbit,fmap,ek,re,te)

      case(code_changeref)
        call fort_warn('TWISS: ','Changeref is nto implemented for MAD-X twiss.')

     case (code_crabcavity)
        call tmcrab(fsec,ftrk,orbit,fmap,dl,ek,re,te)

     case (code_hacdipole, code_vacdipole)
        ! nothing in MAD-X only used for conversion to sixtrack

     case (code_nllens)
        call tmnll(fsec,ftrk,orbit,fmap,ek,re,te)

     case (code_rfmultipole)
        call tmrfmult(fsec,ftrk,orbit,fmap,ek,re,te)
     case (code_changerefp0)
        call tmchp0(ftrk,orbit,fmap,ek,re, te)
     case default !--- anything else:
        ! nil (23, 28, 34)

     end select

end SUBROUTINE tmmap

SUBROUTINE tmbend(ftrk,fcentre,orbit,fmap,el,dl,ek,re,te,code)
  use twtrrfi
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : radiate, deltap, gamma, arad
  use matrices
  use math_constfi, only : zero, one, two, three
  use code_constfi
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for sector bending magnets                         *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     el        (double)  element length.                              *
  !     dl        (double)  slice length.                                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: ftrk, fmap, fcentre
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6), el, dl

  logical :: cplxy
  logical :: kill_ent_fringe, kill_exi_fringe
  integer :: elpar_vl
  integer :: nd, n_ferr, code
  double precision :: f_errors(0:maxferr)
  double precision :: rw(6,6), tw(6,6,6), ek0(6)
  double precision :: x, y
  double precision :: an, sk0, sk1, sk2, sks, tilt, e1, e2, h, h1, h2, hgap, fint, fintx, rhoinv, blen, bvk
  double precision :: dh, corr, ct, st, hx, hy, rfac, pt, h_k

  integer, external :: el_par_vector, node_fd_errors
  double precision, external :: node_value, get_value
  character(len=name_len) :: name
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !---- Initialize.
  EK0 = zero
  RW = EYE
  TW = zero
  ct=0.d0; st=0.d0

  !code = node_value('mad8_type ')
  kill_ent_fringe = node_value('kill_ent_fringe ') .ne. 0d0
  kill_exi_fringe = node_value('kill_exi_fringe ') .ne. 0d0 .or. fcentre



  !---- Test for non-zero length.
  fmap = el .ne. zero
  if (.not. fmap) return

     !-- get element parameters
     elpar_vl = el_par_vector(b_k3s, g_elpar)
     bvk = node_value('other_bv ')
     an = bvk * g_elpar(b_angle)
     tilt = g_elpar(b_tilt)
     e1 = g_elpar(b_e1)
     e2 = g_elpar(b_e2)

     if (code .eq. code_rbend) then
        e1 = e1 + an / two
        e2 = e2 + an / two
     endif

     !---  bvk also applied further down
     sk0 = g_elpar(b_k0)
     sk1 = g_elpar(b_k1)
     sk2 = g_elpar(b_k2)
     h1 = g_elpar(b_h1)
     h2 = g_elpar(b_h2)
     hgap = g_elpar(b_hgap)
     fint = g_elpar(b_fint)
     fintx = g_elpar(b_fintx)
     sks = g_elpar(b_k1s)
     h = an / el
     h_k = h * bvk
     !---- Apply field errors and change coefficients using DELTAP.
     F_ERRORS = zero
     n_ferr = node_fd_errors(f_errors)
     if (sk0 .ne. 0) then 
      f_errors(0) = f_errors(0) + sk0*el - g_elpar(b_angle)
      h_k = sk0 * bvk
    endif



     
!!     if (sk0*el .ne. g_elpar(b_angle)) then
!!        call element_name(name,len(name))
!!        print *, name, ': k0l ~= angle, delta= ', sk0*el - g_elpar(b_angle), g_elpar(b_angle)
!!     endif

     dh = (- h * deltap + bvk * f_errors(0) / el) / (one + deltap) ! dipole term
     sk1 = bvk * (sk1 + f_errors(2) / el) / (one + deltap) ! quad term
     sk2 = bvk * (sk2 + f_errors(4) / el) / (one + deltap) ! sext term
     sks = bvk * (sks + f_errors(3) / el) / (one + deltap) ! skew quad term

     !---  calculate body slice from start (no exit fringe field):
     if (dl .lt. el .and. .not. fcentre) then
       el = dl
       kill_exi_fringe = .true.
     endif

     !---- Half radiation effects at entrance.
     if (ftrk .and. radiate) then
        ct = cos(tilt)
        st = sin(tilt)
        x =   orbit(1) * ct + orbit(3) * st
        y = - orbit(1) * st + orbit(3) * ct
        hx = h + dh + sk1*(x - h*y**2/two) + sks*y + sk2*(x**2 - y**2)/two
        hy = sks * x - sk1*y - sk2*x*y
        rfac = (arad * gamma**3 * el / three) * &
             (hx**2 + hy**2) * (one + h*x) * (one - tan(e1)*x)
        pt = orbit(6)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif

     !---- Body of the dipole.
     !---- Get map for body section
     call tmsect(.true.,dl,h,dh,sk1,sk2,ek,re,te)

     !---- Get map for entrance fringe field and concatenate
     if (.not.kill_ent_fringe) then
        corr = (h_k + h_k) * hgap * fint
        call tmfrng(.true.,h_k,sk1,e1,h1,one,corr,rw,tw)
        call tmcat1(.true.,ek,re,te,ek0,rw,tw,ek,re,te)
     endif

   !---- Get map for exit fringe fields and concatenate
     if (.not.kill_exi_fringe) then
        if (fintx .lt. 0) fintx = fint
        corr = (h_k + h_k) * hgap * fintx
        call tmfrng(.true.,h_k,sk1,e2,h2,-one,corr,rw,tw)
        call tmcat1(.true.,ek0,rw,tw,ek,re,te,ek,re,te)
     endif

     !---- Apply tilt.
     if (tilt .ne. zero) then
        call tmtilt(.true.,tilt,ek,re,te)
        cplxy = .true.
     endif

     !---- Track orbit.
     if (ftrk) then
        call tmtrak(ek,re,te,orbit,orbit)
     endif

     if (fcentre) return

     !---- Half radiation effects at exit.
     if (ftrk .and. radiate) then
        x =   orbit(1) * ct + orbit(3) * st
        y = - orbit(1) * st + orbit(3) * ct
        hx = h + dh + sk1*(x - h*y**2/two) + sks*y + sk2*(x**2 - y**2)/two
        hy = sks * x - sk1*y - sk2*x*y
        rfac = (arad * gamma**3 * el / three) * &
             (hx**2 + hy**2) * (one + h*x) * (one - tan(e2)*x)
        pt = orbit(6)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif

end SUBROUTINE tmbend

SUBROUTINE tmsect(fsec,el,h,dh,sk1,sk2,ek,re,te)
  use twissbeamfi, only : beta, gamma, dtbyds
  use matrices, only: EYE
  use math_constfi, only : zero, one, two, three, four, six, nine, twelve, fifteen
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for a sector dipole without fringe fields.         *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     el        (double)  element length.                              *
  !     h         (double)  reference curvature of magnet.               *
  !     dh        (double)  dipole field error.                          *
  !     sk1       (double)  quadrupole strength.                         *
  !     sk2       (double)  sextupole strengh.                           *
  !     Output:                                                          *
  !     ek(6)     (double)  kick due to dipole.                          *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  logical, intent(IN) :: fsec
  double precision :: el, h, dh, sk1, sk2
  double precision :: ek(6), re(6,6), te(6,6,6)

  double precision :: bi, bi2, bi2gi2
  double precision :: cm, cp, cx, cy, cyy, dd, difsq, dm, dp, dx, dyy
  double precision :: fm, fp, fx, fyy, gx, h2, hx, sm, sp, sumsq, sx, sy, syy
  double precision :: t1, t116, t126, t166, t2, t216, t226, t266
  double precision :: t336, t346, t436, t446, t5, t516, t526, t566
  double precision :: xk, xkl, xklsq, xksq, xs6
  double precision :: y0, y1, y2, y2klsq, y2ksq, yk, ykl, yklsq, yksq, ys2
  double precision :: zc, zd, zf, zs

  double precision, parameter :: twty=20d0, twty2=22d0, twty4=24d0, thty=30d0
  double precision, parameter :: foty2=42d0, fvty6=56d0, svty2=72d0, httwty=120d0
  double precision, parameter :: c1=one, c2=one/two, c3=one/twty4, c4=one/720d0
  double precision, parameter :: s1=one, s2=one/six, s3=one/httwty, s4=one/5040d0
  double precision, parameter :: cg0=one/twty, cg1=5d0/840d0, cg2=21d0/60480d0
  double precision, parameter :: ch0=one/fvty6, ch1=14d0/4032d0, ch2=147d0/443520d0

  !---- Initialize.
  EK = zero
  RE = EYE
  if (fsec) TE = zero

  bi = one / beta
  bi2 = bi * bi
  bi2gi2 = one / (beta * gamma) ** 2

  !---- Horizontal.
  xksq = h**2 + sk1
  xk = sqrt(abs(xksq))
  xkl = xk * el
  xklsq = xksq * el**2

  if (abs(xklsq) .lt. 1e-2) then
     cx = (c1 - xklsq * (c2 - xklsq*c3))
     sx = (s1 - xklsq * (s2 - xklsq*s3)) * el
     dx = (c2 - xklsq * (c3 - xklsq*c4)) * el**2
     fx = (s2 - xklsq * (s3 - xklsq*s4)) * el**3
     gx = (cg0 - xklsq * (cg1 - xklsq*cg2)) * el**5
     hx = (ch0 - xklsq * (ch1 - xklsq*ch2)) * el**7
  else
     if (xklsq .gt. zero) then
        cx = cos(xkl)
        sx = sin(xkl) / xk
     else
        cx = cosh(xkl)
        sx = sinh(xkl) / xk
     endif
     dx = (one - cx) / xksq
     fx = (el  - sx) / xksq
     gx = (three*el - sx*(four-cx)) / (two*xksq**2)
     hx = (fifteen*el - sx*(twty2-nine*cx+two*cx**2)) / (six*xksq**3)
  endif

  re(1,1) = cx
  re(1,2) = sx
  re(1,6) = h * dx * bi
  re(2,1) = - xksq * sx
  re(2,2) = cx
  re(2,6) = h * sx * bi
  re(5,2) = - re(1,6)
  re(5,1) = - re(2,6)
  re(5,6) = el*bi2gi2 - h**2*fx*bi2

  ek(1) = - dh*dx
  ek(2) = - dh*sx
  ek(5) =   h*dh*fx*bi + el*dtbyds

  !---- Vertical.
  yksq = - sk1
  yk = sqrt(abs(yksq))
  ykl = yk*el
  yklsq = yksq*el**2

  if (abs(yklsq) .lt. 1e-2) then
     cy = (c1 - yklsq * (c2 - yklsq*c3))
     sy = (s1 - yklsq * (s2 - yklsq*s3)) * el
  else if (yklsq .gt. zero) then
     cy = cos(ykl)
     sy = sin(ykl) / yk
  else
     cy = cosh(ykl)
     sy = sinh(ykl) / yk
  endif

  re(3,3) = cy
  re(3,4) = sy
  re(4,3) = - yksq * sy
  re(4,4) = cy

  ek(3)   = zero
  ek(4)   = zero

  !---- Second-order terms.
  if (fsec) then
     !---- Pure horizontal terms.
     xs6 = (sk2 + two*h*sk1) / six
     ys2 = (sk2 +     h*sk1) / two
     h2 = h / two

     t116 = xs6 * (three*sx*fx - dx**2) - h * sx**2
     t126 = xs6 * (sx*dx**2 - two*cx*gx) - h * sx * dx
     t166 = xs6 * (dx**3 - two*sx*gx) - h2 * dx**2
     t216 = xs6 * (three*cx*fx + sx*dx)
     t226 = xs6 * (three*sx*fx + dx**2)
     t266 = xs6 * (sx*dx**2 - two*cx*gx)
     t516 = h * xs6 * (three*dx*fx - four*gx) + (sk1/two) * (fx + sx*dx)
     t526 = h * xs6 * (dx**3 - two*sx*gx) + (sk1/two) * dx**2
     t566 = h * xs6 * (three*hx - two*dx*gx) + (sk1/two) * gx - fx

     t1 = (sk1/two) * (dx**2 - sx*fx) - dx
     t2 = (sk1/two) * (el*dx - fx)
     t5 = fx - sk1 * (gx - fx*dx / two)

     te(1,1,1) = - xs6 * (sx**2 + dx) - h2*xksq*sx**2
     te(1,1,2) = (- xs6*dx + h2*cx) * sx
     te(1,2,2) = (- xs6*dx + h2*cx) * dx
     te(1,1,6) = (- h2*t116 + (sk1/four)*el*sx) * bi
     te(1,2,6) = (- h2*t126 + (sk1/four) * (el*dx - fx) - sx/two) *bi
     te(1,6,6) = (- h**2*t166 + h*t1) * bi2 - h2 * dx * bi2gi2
     te(2,1,1) = - xs6 * (one + two*cx) * sx
     te(2,1,2) = - xs6 * (one + two*cx) * dx
     te(2,2,2) = - (two*xs6*dx + h2) * sx
     te(2,1,6) = (- h2*t216 - (sk1/four) * (sx - el*cx)) * bi
     te(2,2,6) = (- h2*t226 + (sk1/four) * el * sx) * bi
     te(2,6,6) = (- h**2*t266 + h*t2) * bi2 - h2 * sx * bi2gi2
     te(5,1,1) = (h2*xs6 * (sx*dx + three*fx) - (sk1/four) * (el - cx*sx)) * bi
     te(5,1,2) = (h2*xs6*dx**2 + (sk1/four)*sx**2) * bi
     te(5,2,2) = (h*xs6*gx - sk1 * (fx - sx*dx) / four - sx/two) * bi
     te(5,1,6) = h2 * ((t516 - sk1 * (el*dx - fx) / two) * bi2 + sx * bi2gi2)
     te(5,2,6) = h2 * ((t526 - sk1 * (dx**2 - sx*fx) / two) * bi2 + dx * bi2gi2)
     te(5,6,6) = (h**2 * (t566 + t5) * bi2 + (three/two) * (h**2*fx - el) * bi2gi2) * bi

     !---- Mixed terms.
     y2ksq = four * yksq
     call tmfoc(el,y2ksq,cyy,syy,dyy,fyy)
     y2klsq = y2ksq * el**2

     if (max(abs(y2klsq),abs(xklsq)) .le. 1e-2) then
        y0 = one
        y1 = xklsq + y2klsq
        y2 = xklsq**2 + xklsq*y2klsq + y2klsq**2
        zc = (y0 - (y1 - y2 / thty) / twelve) * el**2 /   two
        zs = (y0 - (y1 - y2 / foty2) / twty) * el**3 /   six
        zd = (y0 - (y1 - y2 / fvty6) / thty) * el**4 /  twty4
        zf = (y0 - (y1 - y2 / svty2) / foty2) * el**5 / httwty
     else if (xksq .le. zero  .or.  yksq .le. zero) then
        dd = xksq - y2ksq
        zc = (cyy - cx) / dd
        zs = (syy - sx) / dd
        zd = (dyy - dx) / dd
        zf = (fyy - fx) / dd
     else
        sumsq = (xk/two + yk) ** 2
        difsq = (xk/two - yk) ** 2
        call tmfoc(el,sumsq,cp,sp,dp,fp)
        call tmfoc(el,difsq,cm,sm,dm,fm)
        zc = sp * sm / two
        zs = (sp*cm - cp*sm) / (four*xk*yk)
        if (xksq .gt. y2ksq) then
           zd = (dyy - zc) / xksq
           zf = (fyy - zs) / xksq
        else
           zd = (dx - zc) / y2ksq
           zf = (fx - zs) / y2ksq
        endif
     endif

     t336 = sk2 * (cy*zd - two*sk1*sy*zf) + h * sk1 * fx * sy
     t346 = sk2 * (sy*zd - two*cy*zf) + h * fx * cy
     t436 = two * ys2 * fx * cy - sk2 * sk1 * (sy*zd - two*cy*zf)
     t446 = two * ys2 * fx * sy - sk2 * (cy*zd - two*sk1*sy*zf)

     te(1,3,3) = + sk2*sk1*zd + ys2*dx
     te(1,3,4) = + sk2*zs/two
     te(1,4,4) = + sk2*zd - h2*dx
     te(2,3,3) = + sk2*sk1*zs + ys2*sx
     te(2,3,4) = + sk2*zc/two
     te(2,4,4) = + sk2*zs - h2*sx
     te(3,1,3) = + sk2*(cy*zc/two - sk1*sy*zs) + h2*sk1*sx*sy
     te(3,1,4) = + sk2*(sy*zc/two - cy*zs) + h2*sx*cy
     te(3,2,3) = + sk2*(cy*zs/two - sk1*sy*zd) + h2*sk1*dx*sy
     te(3,2,4) = + sk2*(sy*zs/two - cy*zd) + h2*dx*cy
     te(3,3,6) = (h2*t336 - sk1*el*sy/four) * bi
     te(3,4,6) = (h2*t346 - (sy + el*cy) / four) * bi
     te(4,1,3) = sk2*sk1*(cy*zs - sy*zc/two) + ys2*sx*cy
     te(4,1,4) = sk2*(sk1*sy*zs - cy*zc/two) + ys2*sx*sy
     te(4,2,3) = sk2*sk1*(cy*zd - sy*zs/two) + ys2*dx*cy
     te(4,2,4) = sk2*(sk1*sy*zd - cy*zs/two) + ys2*dx*sy
     te(4,3,6) = (h2*t436 + sk1 * (sy - el*cy) / four) * bi
     te(4,4,6) = (h2*t446 - sk1*el*sy/four) * bi
     te(5,3,3) = (- h*sk2*sk1*zf - h*ys2*fx + sk1*(el-cy*sy)/four)*bi
     te(5,3,4) = (- h*sk2*zd/two - sk1*sy**2/four) * bi
     te(5,4,4) = (- h*sk2*zf + h*h2*fx - (el + sy*cy)/four) * bi
     call tmsymm(te)

     !---- Effect of dipole error.
     if (dh .ne. zero) then
        re(1,1) = re(1,1) + dh * t116
        re(1,2) = re(1,2) + dh * t126
        re(1,6) = re(1,6) + dh * (two*h*t166 - t1) * bi
        re(2,1) = re(2,1) + dh * (t216 - h*sx)
        re(2,2) = re(2,2) + dh * t226
        re(2,6) = re(2,6) + dh * (two*h*t266 - t2) * bi
        re(5,1) = re(5,1) - dh * t516 * bi
        re(5,2) = re(5,2) - dh * (t526 - dx) * bi
        re(5,6) = re(5,6) - dh * h * ((two*t566 + t5) * bi2 + fx * bi2gi2)
        re(3,3) = re(3,3) - dh * t336
        re(3,4) = re(3,4) - dh * t346
        re(4,3) = re(4,3) - dh * t436
        re(4,4) = re(4,4) - dh * t446

        ek(1) = ek(1) - dh**2 * t166
        ek(2) = ek(2) - dh**2 * t266
        ek(5) = ek(5) + dh**2 * t566 * bi
     endif
  endif

end SUBROUTINE tmsect

SUBROUTINE tmfrng(fsec,h,sk1,edge,he,sig,corr,re,te)
  use matrices, only : EYE
  use math_constfi, only : zero, one, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for fringe field of a dipole.                      *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     h         (double)  curvature of magnet body.                    *
  !     sk1       (double)  quadrupole strength in magnet body.          *
  !     edge      (double)  edge focussing angle.                        *
  !     he        (double)  curvature of pole face.                      *
  !     sig       (double)  sign: +1 for entry, -1 for exit.             *
  !     corr      (double)  correction factor according to slac 75.      *
  !     Output:                                                          *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec
  double precision :: h, sk1, edge, he, sig, corr
  double precision :: re(6,6), te(6,6,6)

  double precision :: hh, psip, secedg, tanedg

  !---- Linear terms.
  RE = EYE
  tanedg = tan(edge)
  secedg = one / cos(edge)
  psip = edge - corr * secedg * (one + sin(edge)**2)
  re(2,1) = + h * tanedg
  re(4,3) = - h * tan(psip)

  !---- Second-order terms.
  if (fsec) then
     TE = zero
     hh = sig * (h/two)
     te(1,1,1) = - hh * tanedg**2
     te(1,3,3) = + hh * secedg**2
     te(2,1,1) = (h/two) * he * secedg**3 + sk1 * tanedg
     te(2,1,2) = - te(1,1,1)
     te(2,3,3) = hh * h * tanedg**3 - te(2,1,1)
     te(2,3,4) = + te(1,1,1)
     te(3,1,3) = - te(1,1,1)
     te(4,1,3) = - te(2,1,1)
     te(4,1,4) = + te(1,1,1)
     te(4,2,3) = - te(1,3,3)
     if (sig .gt. zero) then
        te(2,3,3) = te(2,3,3) + (h*secedg)**2 * tanedg/two
     else
        te(2,1,1) = te(2,1,1) - (h*tanedg)**2 * tanedg/two
        te(4,1,3) = te(4,1,3) + (h*secedg)**2 * tanedg/two
     endif
     call tmsymm(te)
  endif

end SUBROUTINE tmfrng

SUBROUTINE tmtilt(fsec,tilt,ek,r,t)
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Apply TILT to a TRANSPORT map.                                   *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     tilt      (double)  roll angle. The inverse is used to rotate    *
  !     the matrix and tensor                                            *
  !     ek(6)     (double)  element kick, unrotated.                     *
  !     r(6,6)    (double)  transfer matrix, unrotated.                  *
  !     t(6,6,6)  (double)  second order terms, unrotated.               *
  !     Output:                                                          *
  !     ek(6)     (double)  element kick, rotated.                       *
  !     r(6,6)    (double)  transfer matrix, rotated.                    *
  !     t(6,6,6)  (double)  second order terms, rotated.                 *
  !----------------------------------------------------------------------*
  logical :: fsec
  double precision :: tilt
  double precision, intent(IN OUT) :: ek(6), r(6,6), t(6,6,6)

  integer :: i, j, k
  double precision :: c, s, r1j, r2j, ri1, ri2, xx
  double precision :: t1jk, t2jk, ti1k, ti2k, tij1, tij2

  !---- don't bother tilting by 0 degrees..
  if (tilt .eq. zero) return

  c = cos(tilt)
  s = sin(tilt)

  !---- Rotate at entrance.
  do i = 1, 6
     ri1 = r(i,1)
     r(i,1) = ri1 * c - r(i,3) * s
     r(i,3) = ri1 * s + r(i,3) * c
     ri2 = r(i,2)
     r(i,2) = ri2 * c - r(i,4) * s
     r(i,4) = ri2 * s + r(i,4) * c
     if (fsec) then
        do k = 1, 6
           ti1k = t(i,1,k)
           t(i,1,k) = ti1k * c - t(i,3,k) * s
           t(i,3,k) = ti1k * s + t(i,3,k) * c
           ti2k = t(i,2,k)
           t(i,2,k) = ti2k * c - t(i,4,k) * s
           t(i,4,k) = ti2k * s + t(i,4,k) * c
        enddo
        do j = 1, 6
           tij1 = t(i,j,1)
           t(i,j,1) = tij1 * c - t(i,j,3) * s
           t(i,j,3) = tij1 * s + t(i,j,3) * c
           tij2 = t(i,j,2)
           t(i,j,2) = tij2 * c - t(i,j,4) * s
           t(i,j,4) = tij2 * s + t(i,j,4) * c
        enddo
     endif
  enddo

  !---- Rotate kick.
  xx = ek(1)
  ek(1) = xx * c - ek(3) * s
  ek(3) = xx * s + ek(3) * c
  xx = ek(2)
  ek(2) = xx * c - ek(4) * s
  ek(4) = xx * s + ek(4) * c

  !---- Rotate at exit.
  do j = 1, 6
     r1j = r(1,j)
     r(1,j) = c * r1j - s * r(3,j)
     r(3,j) = s * r1j + c * r(3,j)
     r2j = r(2,j)
     r(2,j) = c * r2j - s * r(4,j)
     r(4,j) = s * r2j + c * r(4,j)
     if (fsec) then
        do k = 1, 6
           t1jk = t(1,j,k)
           t(1,j,k) = c * t1jk - s * t(3,j,k)
           t(3,j,k) = s * t1jk + c * t(3,j,k)
           t2jk = t(2,j,k)
           t(2,j,k) = c * t2jk - s * t(4,j,k)
           t(4,j,k) = s * t2jk + c * t(4,j,k)
        enddo
     endif
  enddo

end SUBROUTINE tmtilt

SUBROUTINE tmcorr(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)
  use twtrrfi
  use math_constfi, only : zero, one, two, three, half
  use twissbeamfi, only : radiate, deltap, gamma, arad
  use code_constfi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for orbit correctors.                              *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     el        (double)  element length.                              *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical, intent(IN) :: fsec, ftrk, fcentre
  logical, intent(OUT) :: fmap
  double precision, intent(IN OUT) :: orbit(6), el, dl
  double precision, intent(OUT) :: ek(6), re(6,6), te(6,6,6)

  logical :: cplxy
  integer :: i, code, n_ferr
  double precision :: f_errors(0:maxferr)
  double precision :: rfac, pt, tilt, bvk
  double precision :: xkick, ykick, dpx, dpy, xau, div

  integer, external :: node_fd_errors
  double precision, external :: node_value, get_value
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !--- Initialization
  rfac=0.d0

  if ( .not. ftrk) then
     !---- No orbit track desired, use drift map.
     call tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)

  else
     !---- Tracking desired, use corrector map.
     !---- Initialize.
     div = el ; if (el .eq. zero) div = one
     bvk = node_value('other_bv ')
     tilt = -node_value('tilt ')

     F_ERRORS = zero
     n_ferr = node_fd_errors(f_errors)

     !---- Original setting.
     code = node_value('mad8_type ')
     select case (code)

       case (code_hkicker)
          xkick=bvk*(node_value('kick ')+node_value('chkick '))
          ykick=zero

       case (code_kicker, code_tkicker)
          xkick=bvk*(node_value('hkick ')+node_value('chkick '))
          ykick=bvk*(node_value('vkick ')+node_value('cvkick '))

       case (code_vkicker)
          xkick=zero
          ykick=bvk*(node_value('kick ')+node_value('cvkick '))

       case default
          xkick=zero
          ykick=zero

     end select

     xkick=xkick+bvk*(f_errors(0)/div);
     ykick=ykick+bvk*(f_errors(1)/div);

     xau = xkick
     xkick = xkick*cos(tilt)+ykick*sin(tilt)
     ykick =  -xau*sin(tilt)+ykick*cos(tilt)

     !---- Sum up total kicks.
     dpx = xkick / (one + deltap)
     dpy = ykick / (one + deltap)
     if (dpy .ne. zero) cplxy = .true.

     !---- Half kick at entrance.
     orbit(2) = orbit(2) + half * dpx
     orbit(4) = orbit(4) + half * dpy

     !---- Half radiation effects at entrance.
     if (radiate  .and.  el.ne.zero) then
        rfac = arad * gamma**3 * (dpx**2 + dpy**2) / (three * el)
        pt = orbit(6)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif

     !---- Drift to end.
     if (el .ne. zero) then
       call tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)
       if (fcentre) return
     endif

     !---- Half radiation effects at exit.
     if (radiate  .and.  el.ne.zero) then
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif

     !---- Half kick at exit.
     orbit(2) = orbit(2) + half * dpx
     orbit(4) = orbit(4) + half * dpy

     fmap = .true.
  endif

end SUBROUTINE tmcorr

INTEGER FUNCTION Factorial(n)
  implicit none
  integer, intent(in) :: n
  integer :: i, Ans

  Ans = 1
  do i = 1, n
    Ans = Ans * i
  enddo
  Factorial = Ans
END FUNCTION Factorial

SUBROUTINE tmmult_cf_short(fsec, ftrk, orbit, fmap, re, te)
  use twtrrfi, only : maxmul, maxferr
  use twissbeamfi, only : deltap, beta
  use math_constfi, only : zero, one, two, three
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Computes thin-lens kick through combined-function magnet.        *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !     Detailed description:                                            *
  !     See Phys. Rev. AccelBeams 19.054002 by M.Titze                   *
  !                                                                      *
  !     Note:                                                            *
  !     This 'short' version yields the correct chromaticity of a        *
  !     sliced sector CFM only if TWISS is used without the 'chrom'      *
  !     option. For a correct treatment of the chrom option,             * 
  !     tmmult_cf must be used.                                     *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  integer :: nord, k, j, nn, ns, bvk, iord, n_ferr
  integer, external :: Factorial
  double precision :: dpx, dpy, tilt, kx, ky, elrad, bp1, etahat, h0
  double precision :: an, angle, dtmp
  double precision :: normal(0:maxmul), skew(0:maxmul), f_errors(0:maxferr)
  double precision :: orbit(6), re(6,6), te(6,6,6), tilt2
  double complex :: kappa, barkappa, sum0, del_p_g, pkick, dxdpg, dydpg, &
                    dxx, dxy, dyy, rp, rm
  double complex :: lambda(0:maxmul)
  double complex :: g(0:maxmul, 0:maxmul)

  double precision, external :: node_value
  integer, external :: node_fd_errors
  fmap = .true.

  ! Read magnetic field components & fill lambda's according to field
  ! components relative to given plane
  normal = zero ; call get_node_vector('knl ', nn, normal)
  skew   = zero ; call get_node_vector('ksl ', ns, skew)
  nord = max(nn, ns)
  tilt = node_value('tilt ')
  elrad = node_value('lrad ')

  F_ERRORS(0:maxferr) = zero
  n_ferr = node_fd_errors(f_errors)
  bvk = node_value('other_bv ')
  tilt2 = 0 ! A parameter describing the relative tilt between
            ! the dipole component and the higher-order components of the CFM

  !####SETTING UP THE MULTIPOLES
  an = node_value('angle ')
  f_errors(0) = f_errors(0) + normal(0) ! The zero-component of the B-field

  !Below here should not be commented output
  !---- Other components and errors.
  nord = 0
  ! that loop should start at one since nominal dipole strength already taken into account above
  !needs to be here though
  do iord = 0, max(nn, ns, n_ferr/2-1)
  !   get the maximum effective order; loop runs over maximum of user given values
     if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero .or. &
          normal(iord).ne.zero .or. skew(iord).ne.zero) nord = iord+1 !  why  +1
  enddo

  do iord = 1, nord
     f_errors(2*iord)   = (normal(iord) + f_errors(2*iord))
     f_errors(2*iord+1) = (skew(iord)   + f_errors(2*iord+1))
     if (tilt .ne. zero) then
        if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero) then
           angle = atan2(f_errors(2*iord+1), f_errors(2*iord)) / (iord+1) - tilt
        else
           angle = -tilt
        endif
        angle = (iord+1) * angle
        dtmp = sqrt(f_errors(2*iord)**2 + f_errors(2*iord+1)**2)
        f_errors(2*iord)   = dtmp * cos(angle)
        f_errors(2*iord+1) = dtmp * sin(angle)
     endif
     f_errors(2*iord)   = bvk * f_errors(2*iord)
     f_errors(2*iord+1) = bvk * f_errors(2*iord + 1)
  enddo
  !Done with all the setting up...

  ! The "normal" components are considered here as the expansion coefficients of
  ! B_y wrt. the reference plane, while the "skew" components are considered as the
  ! corresponding expansion coefficients of B_x, see documentation. This can
  ! be modified in the future, in particular to use tilted components,
  ! but bare in mind that the bending curvature (the
  ! lambda(0) terms) should be unchanged.
  !
  ! The above means precisely, that we currently implemented the following scheme:
  !
  ! B_y |_{\varphi = tilt} + i B_x |_{\varphi = tilt} = \sum_{k = 0}^nord \lambda_k r^k
  !
  ! with complex coefficients \lambda_k and in which
  !
  ! Im[\lambda_k] = 1/k! \partial^k B_x / \partial_r^k |_{\varphi = tilt} ,
  ! Re[\lambda_k] = 1/k! \partial^k B_y / \partial_r^k |_{\varphi = tilt} .
  !
  ! play the role as the k'th skew- and normal field component.

  do k = 0, nord
     ! The factor (one + deltap) below is taken from the original MAD-X routine.
     lambda(k) = (f_errors(2*k) + (0, 1)*f_errors(2*k+1))/elrad/Factorial(k)/(one + deltap)
  enddo

  ! Set the curvature of the MAD-X reference trajectory; Due to MAD-X standards,
  ! this trajectory should be independent on deltap
  if (an .eq. 0) then
     ! The special case an == 0 is treated differently for the purpose of backwards compatibility
     kx = real(lambda(0))*(one + deltap)
     ky = -aimag(lambda(0))*(one + deltap)
  else
     kx = an/elrad*cos(tilt)
     ky = an/elrad*sin(tilt)
  endif

  !  N.B. B_y |_{\varphi = tilt, r = 0} = kx
  !       B_x |_{\varphi = tilt, r = 0} = -ky, see Eqs. (18) in
  ! Phys. Rev. AccelBeams 19.054002

  kappa = kx + (0, 1)*ky
  barkappa = conjg(kappa)

  ! Now fill up the g_{ij}'s for j = 0, ..., i and i = 0, ..., nord + 1.
  g(0, 0) = (0, 0)
  g(1, 0) = -lambda(0)
  g(1, 1) = conjg(g(1, 0))

  do k = 1, nord
     do j = 0, k - 1
        ! Eq. (6), in Ref. above
        g(k + 1, j + 1) = (barkappa*g(k, j + 1)*(j + one)*(j - k + three/two) +  &
             kappa*g(k, j)*(k - j)*(one/two - j))/(k - j)/(j + one)
     enddo
     ! Eq. (8) in Ref. above
     sum0 = 0
     do j = 1, k
       sum0 = sum0 - (k + 1 - j)*g(k + 1, j)*exp(-two*(0, 1)*j*tilt2)
     enddo
     g(k + 1, 0) = ( sum0 - two**k*exp(-(0, 1)*k*tilt2)*( lambda(k) &
                    + one/two*(barkappa*exp((0, 1)*tilt2) + kappa*exp(-(0, 1)*tilt2)) &
                    *lambda(k - 1) ) )/(k + one)
     g(k + 1, k + 1) = conjg(g(k + 1, 0))
  enddo

  etahat = sqrt(two*orbit(6)/beta + orbit(6)**2 + one) - one ! etahat = deltap of individual particle
  h0 = sqrt((one + etahat)**2 - orbit(2)**2 - orbit(4)**2)

  if (ftrk) then
     rp = (orbit(1) + (0, 1)*orbit(3))/two
     rm = conjg(rp)
     ! Compute \partial_+ G using Eq. (7) in Ref. above
     del_p_g = 0
     do k = 1, nord
        sum0 = 0
        do j = 0, k - 1
           sum0 = sum0 + (k - j)*g(k, j)*rp**(k - 1 - j)*rm**j
        enddo
        del_p_g = del_p_g + sum0
     enddo
     ! Now compute kick (Eqs. (38) in Ref. above)
     pkick = elrad*(barkappa*h0 + del_p_g)
     dpx = real(pkick)
     dpy = - aimag(pkick)
     orbit(1) = orbit(1) + elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)/h0
     orbit(2) = orbit(2) + dpx
     orbit(3) = orbit(3) + elrad*(kx*orbit(1) + ky*orbit(3))*orbit(4)/h0
     orbit(4) = orbit(4) + dpy
     ! N.B. orbit(5) = \sigma/beta and orbit(6) = beta*p_\sigma
     orbit(5) = orbit(5) - elrad*(kx*orbit(1) + ky*orbit(3)) &
                *(one + beta*orbit(6))/beta/h0
  endif
  ! First-order terms by derivation of Eqs. (39) in Ref. above, at zero
  ! re(6,6) is assumed to be a unit matrix as input
  ! Eqs. (39) are required to obtain agreement to the thick sectormap.

  if (nord .ge. 1) then
     ! The next two expressions emerge by the first
     ! derivative of \partial_+ G wrt. x and y, see documentation.
     ! the factors (one + deltap) means that we take the derivative with respect to a given
     ! energy-offset, so that if deltap != 0, the derivative is evaluated at a larger radius.
     dxdpg = elrad/two*(two*g(2, 0) + g(2, 1))
     dydpg = elrad/two*(0, 1)*(two*g(2, 0) - g(2, 1))
     re(2, 1) = real(dxdpg)
     re(2, 3) = real(dydpg)
     re(2, 6) = elrad*kx/beta*(one + beta*orbit(6))/h0
     re(4, 1) = - aimag(dxdpg)
     re(4, 3) = - aimag(dydpg)
     re(4, 6) = elrad*ky/beta*(one + beta*orbit(6))/h0
     re(5, 1) = - elrad*kx/beta*(one + beta*orbit(6))/h0
     re(5, 3) = - elrad*ky/beta*(one + beta*orbit(6))/h0
  endif
  ! Second-order terms by derivation of Eqs. (39) in Ref. above, at zero
  ! te(6,6,6) is assumed to be a zero tensor as input
  if ((fsec) .and. (nord .ge. 2)) then
     ! The next three expressions emerge by the second
     ! derivative of \partial_+ G wrt. x and y at zero, see documentation.
     ! The additional factor two in the end accounts for the fact
     ! that the te(6,6,6) tensor and the 2nd derivative are related
     ! by the usual Taylor-expansion factorials.
     dxx = elrad/two*(three*g(3, 0) + two*g(3, 1) + g(3, 2))/two
     dxy = elrad/two*(0, 1)*(three*g(3, 0) - g(3, 2))/two
     dyy = -elrad/two*(three*g(3, 0) - two*g(3, 1) + g(3, 2))/two
     bp1 = one - one/beta**2
     te(1, 1, 2) = elrad*kx*(one/h0 + orbit(2)**2/h0**3)/two
     te(1, 2, 1) = te(1, 1, 2)
     te(1, 2, 3) = elrad*ky*(one/h0 + orbit(2)**2/h0**3)/two
     te(1, 3, 2) = te(1, 2, 3)
     te(2, 1, 1) = real(dxx)   ! cf
     te(2, 1, 3) = real(dxy)   ! cf
     te(2, 3, 1) = te(2, 1, 3) ! cf
     te(2, 3, 3) = real(dyy)   ! cf
     te(2, 2, 2) = - te(1, 1, 2)
     te(2, 4, 4) = - te(1, 1, 2)
     te(2, 6, 6) = te(1, 1, 2)*bp1
     te(3, 1, 4) = te(1, 1, 2)
     te(3, 3, 4) = te(1, 2, 3)
     te(3, 4, 1) = te(1, 1, 2)
     te(3, 4, 3) = te(1, 2, 3)
     te(4, 1, 1) = - aimag(dxx)  ! cf
     te(4, 1, 3) = - aimag(dxy)  ! cf
     te(4, 3, 1) = te(4, 1, 3)   ! cf
     te(4, 3, 3) = - aimag(dyy)  ! cf
     te(4, 2, 2) = - te(1, 2, 3)
     te(4, 4, 4) = - te(1, 2, 3)
     te(4, 6, 6) = te(1, 2, 3)*bp1
     te(5, 1, 6) = - te(1, 1, 2)*bp1
     te(5, 3, 6) = - te(1, 2, 3)*bp1
     te(5, 6, 1) = - te(1, 1, 2)*bp1
     te(5, 6, 3) = - te(1, 2, 3)*bp1
  endif
end SUBROUTINE tmmult_cf_short

SUBROUTINE tmmult_cf(fsec, ftrk, orbit, fmap, re, te)
  use twtrrfi, only : maxmul, maxferr
  use twissbeamfi, only : deltap, beta
  use math_constfi, only : zero, one, two, three
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Computes thin-lens kick through combined-function magnet.        *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !     Detailed description:                                            *
  !     See Phys. Rev. AccelBeams 19.054002 by M.Titze                   *
  !     Implementation of full map (Eq. (2.38))                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  integer :: nord, k, j, nn, ns, bvk, iord, n_ferr
  integer, external :: Factorial
  double precision :: dpx, dpy, tilt, kx, ky, elrad, bp1, etahat, h0
  double precision :: an, angle, dtmp
  double precision :: normal(0:maxmul), skew(0:maxmul), f_errors(0:maxferr)
  double precision :: orbit(6), re(6,6), te(6,6,6), tilt2
  double complex :: kappa, barkappa, sum0, del_p_g, pkick, dxdpg, dydpg, &
                    dxx, dxy, dyy, rp, rm
  double complex :: lambda(0:maxmul)
  double complex :: g(0:maxmul, 0:maxmul)

  double complex :: sigma(0:maxmul, 0:maxmul), &
                    dxsigma(0:maxmul, 0:maxmul), &
                    dysigma(0:maxmul, 0:maxmul), &
                    dxxsigma(0:maxmul, 0:maxmul), &
                    dxysigma(0:maxmul, 0:maxmul), &
                    dyysigma(0:maxmul, 0:maxmul)

  double precision, external :: node_value
  integer, external :: node_fd_errors
  fmap = .true.

  ! Read magnetic field components & fill lambda's according to field
  ! components relative to given plane
  normal = zero ; call get_node_vector('knl ', nn, normal)
  skew   = zero ; call get_node_vector('ksl ', ns, skew)
  nord = max(nn, ns)
  tilt = node_value('tilt ')
  elrad = node_value('lrad ')

  F_ERRORS(0:maxferr) = zero
  n_ferr = node_fd_errors(f_errors)
  bvk = node_value('other_bv ')
  tilt2 = 0 ! A parameter describing the relative tilt between
            ! the dipole component and the higher-order components of the CFM

  !####SETTING UP THE MULTIPOLES
  an = node_value('angle ')
  f_errors(0) = f_errors(0) + normal(0) ! The zero-component of the B-field

  !Below here should not be commented output
  !---- Other components and errors.
  nord = 0
  ! that loop should start at one since nominal dipole strength already taken into account above
  !needs to be here though
  do iord = 0, max(nn, ns, n_ferr/2-1)
  !   get the maximum effective order; loop runs over maximum of user given values
     if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero .or. &
          normal(iord).ne.zero .or. skew(iord).ne.zero) nord = iord+1 !  why  +1
  enddo

  do iord = 1, nord
     f_errors(2*iord)   = (normal(iord) + f_errors(2*iord))
     f_errors(2*iord+1) = (skew(iord)   + f_errors(2*iord+1))
     if (tilt .ne. zero) then
        if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero) then
           angle = atan2(f_errors(2*iord+1), f_errors(2*iord)) / (iord+1) - tilt
        else
           angle = -tilt
        endif
        angle = (iord+1) * angle
        dtmp = sqrt(f_errors(2*iord)**2 + f_errors(2*iord+1)**2)
        f_errors(2*iord)   = dtmp * cos(angle)
        f_errors(2*iord+1) = dtmp * sin(angle)
     endif
     f_errors(2*iord)   = bvk * f_errors(2*iord)
     f_errors(2*iord+1) = bvk * f_errors(2*iord + 1)
  enddo
  !Done with all the setting up...

  ! The "normal" components are considered here as the expansion coefficients of
  ! B_y wrt. the reference plane, while the "skew" components are considered as the
  ! corresponding expansion coefficients of B_x, see documentation. This can
  ! be modified in the future, in particular to use tilted components,
  ! but bare in mind that the bending curvature (the
  ! lambda(0) terms) should be unchanged.
  !
  ! The above means precisely, that we currently implemented the following scheme:
  !
  ! B_y |_{\varphi = tilt} + i B_x |_{\varphi = tilt} = \sum_{k = 0}^nord \lambda_k r^k
  !
  ! with complex coefficients \lambda_k and in which
  !
  ! Im[\lambda_k] = 1/k! \partial^k B_x / \partial_r^k |_{\varphi = tilt} ,
  ! Re[\lambda_k] = 1/k! \partial^k B_y / \partial_r^k |_{\varphi = tilt} .
  !
  ! play the role as the k'th skew- and normal field component.

  do k = 0, nord
     ! The factor (one + deltap) below is taken from the original MAD-X routine.
     lambda(k) = (f_errors(2*k) + (0, 1)*f_errors(2*k+1))/elrad/Factorial(k)/(one + deltap)
  enddo

  ! Set the curvature of the MAD-X reference trajectory; Due to MAD-X standards,
  ! this trajectory should be independent on deltap
  if (an .eq. 0) then
     ! The special case an == 0 is treated differently for the purpose of backwards compatibility
     kx = real(lambda(0))*(one + deltap)
     ky = -aimag(lambda(0))*(one + deltap)
  else
     kx = an/elrad*cos(tilt)
     ky = an/elrad*sin(tilt)
  endif

  !  N.B. B_y |_{\varphi = tilt, r = 0} = kx
  !       B_x |_{\varphi = tilt, r = 0} = -ky, see Eqs. (18) in
  ! Phys. Rev. AccelBeams 19.054002

  kappa = kx + (0, 1)*ky
  barkappa = conjg(kappa)

  ! Now fill up the g_{ij}'s for j = 0, ..., i and i = 0, ..., nord + 1.
  g(0, 0) = (0, 0)
  g(1, 0) = -lambda(0)
  g(1, 1) = conjg(g(1, 0))

  do k = 1, nord
     do j = 0, k - 1
        ! Eq. (6), in Ref. above
        g(k + 1, j + 1) = (barkappa*g(k, j + 1)*(j + one)*(j - k + three/two) +  &
             kappa*g(k, j)*(k - j)*(one/two - j))/(k - j)/(j + one)
     enddo
     ! Eq. (8) in Ref. above
     sum0 = 0
     do j = 1, k
       sum0 = sum0 - (k + 1 - j)*g(k + 1, j)*exp(-two*(0, 1)*j*tilt2)
     enddo
     g(k + 1, 0) = ( sum0 - two**k*exp(-(0, 1)*k*tilt2)*( lambda(k) &
                    + one/two*(barkappa*exp((0, 1)*tilt2) + kappa*exp(-(0, 1)*tilt2)) &
                    *lambda(k - 1) ) )/(k + one)
     g(k + 1, k + 1) = conjg(g(k + 1, 0))
  enddo

  etahat = sqrt(two*orbit(6)/beta + orbit(6)**2 + one) - one ! etahat = deltap of individual particle
  h0 = sqrt((one + etahat)**2 - orbit(2)**2 - orbit(4)**2)

  if (ftrk) then
     rp = (orbit(1) + (0, 1)*orbit(3))/two
     rm = conjg(rp)
     ! Compute \partial_+ G using Eq. (7) in Ref. above
     del_p_g = 0
     do k = 1, nord
        sum0 = 0
        do j = 0, k - 1
           sum0 = sum0 + (k - j)*g(k, j)*rp**(k - 1 - j)*rm**j
        enddo
        del_p_g = del_p_g + sum0
     enddo
     ! Now compute kick (Eqs. (38) in Ref. above)
     pkick = elrad*(barkappa*h0 + del_p_g)
     dpx = real(pkick)
     dpy = - aimag(pkick)
     orbit(1) = orbit(1) + elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)/h0
     orbit(2) = orbit(2) + dpx
     orbit(3) = orbit(3) + elrad*(kx*orbit(1) + ky*orbit(3))*orbit(4)/h0
     orbit(4) = orbit(4) + dpy
     ! N.B. orbit(5) = \sigma/beta and orbit(6) = beta*p_\sigma
     orbit(5) = orbit(5) - elrad*(kx*orbit(1) + ky*orbit(3)) &
                *(one + beta*orbit(6))/beta/h0
  endif
  ! First-order terms by derivation of Eqs. (39) in Ref. above, at zero
  ! re(6,6) is assumed to be a unit matrix as input
  ! Eqs. (39) are required to obtain agreement to the thick sectormap.

  if (nord .ge. 1) then
     ! The next two expressions emerge by the first
     ! derivative of \partial_+ G wrt. x and y, see documentation.
     ! the factors (one + deltap) means that we take the derivative with respect to a given
     ! energy-offset, so that if deltap != 0, the derivative is evaluated at a larger radius.

     do k = 0, nord
         do j = 0, k
             ! N.B.: j = k + 1 leads to sigma(k, k + 1) = 0, so the inner loop will go only up to k
             sigma(k, j) = elrad*(k + one - j)*g(k + 1, j)
         enddo
     enddo

     do k = 0, nord - 1
         do j = 0, k
             dxsigma(k, j) = ((k + one - j)*sigma(k + 1, j) + (j + one)*sigma(k + 1, j + 1))/two
             dysigma(k, j) = (0, 1)*((k + one - j)*sigma(k + 1, j) - (j + one)*sigma(k + 1, j + 1))/two
         enddo
     enddo

     dxdpg = 0
     dydpg = 0
     do k = 0, nord - 1
         do j = 0, k
             dxdpg = dxdpg + dxsigma(k, j)*rp**(k - j)*rm**j
             dydpg = dydpg + dysigma(k, j)*rp**(k - j)*rm**j
         enddo
     enddo

     ! if evaluated at zero, one can instead use:
     !dxdpg = elrad/two*(two*g(2, 0) + g(2, 1))
     !dydpg = elrad/two*(0, 1)*(two*g(2, 0) - g(2, 1))

     re(1, 1) = one + elrad*kx*orbit(2)/h0
     re(1, 2) = elrad*(kx*orbit(1) + ky*orbit(3))*(one/h0 + orbit(2)**2/h0**3)
     re(1, 3) = elrad*ky*orbit(2)/h0
     re(1, 4) = elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)*orbit(4)/h0**3
     re(1, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)/h0**3/beta*(one + beta*orbit(6))

     re(2, 1) = real(dxdpg)
     re(2, 2) = one - elrad*kx*orbit(2)/h0 
     re(2, 3) = real(dydpg)
     re(2, 4) = -elrad*kx*orbit(4)/h0 
     re(2, 6) = elrad*kx/beta*(one + beta*orbit(6))/h0

     re(3, 1) = elrad*kx*orbit(4)/h0
     re(3, 2) = elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)*orbit(4)/h0**3
     re(3, 3) = one + elrad*ky*orbit(4)/h0
     re(3, 4) = elrad*(kx*orbit(1) + ky*orbit(3))*(one/h0 + orbit(4)**2/h0**3)
     re(3, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*orbit(4)/h0**3/beta*(one + beta*orbit(6))

     re(4, 1) = - aimag(dxdpg)
     re(4, 2) = elrad*ky*orbit(2)/h0
     re(4, 3) = - aimag(dydpg)
     re(4, 4) = one - elrad*ky*orbit(4)/h0
     re(4, 6) = elrad*ky/beta*(one + beta*orbit(6))/h0

     re(5, 1) = -elrad*kx/beta*(one + beta*orbit(6))/h0
     re(5, 2) = -elrad*orbit(2)*(kx*orbit(1) + ky*orbit(3))*(one + beta*orbit(6))/beta/h0**3
     re(5, 3) = -elrad*ky/beta*(one + beta*orbit(6))/h0
     re(5, 4) = -elrad*orbit(4)*(kx*orbit(1) + ky*orbit(3))*(one + beta*orbit(6))/beta/h0**3
     re(5, 5) = one
     re(5, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*(beta**2/h0 - (one + beta*orbit(6))**2/h0**3)/beta**2

     re(6, 6) = one
  endif
  ! Second-order terms by derivation of Eqs. (39) in Ref. above, at zero
  ! te(6,6,6) is assumed to be a zero tensor as input
  if ((fsec) .and. (nord .ge. 2)) then
     ! The next three expressions emerge by the second
     ! derivative of \partial_+ G wrt. x and y at zero, see documentation.
     ! The additional factor two in the end accounts for the fact
     ! that the te(6,6,6) tensor and the 2nd derivative are related
     ! by the usual Taylor-expansion factorials.
     do k = 0, nord - 2
         do j = 0, k
             dxxsigma(k, j) = ((k + one - j)*((k + two - j)*sigma(k + 2, j) + two*(j + one)*sigma(k + 2, j + 1)) + &
                              (j + one)*(j + two)*sigma(k + 2, j + 2))/4
             dxysigma(k, j) = (0, 1)*((k + one - j)*(k + two - j)*sigma(k + 2, j) &
                              - (j + one)*(j + two)*sigma(k + 2, j + 2))/4
             dyysigma(k, j) = -((k + one - j)*((k + two - j)*sigma(k + 2, j) - two*(j + one)*sigma(k + 2, j + 1)) + &
                              (j + one)*(j + two)*sigma(k + 2, j + 2))/4
         enddo
     enddo

     dxx = 0
     dxy = 0
     dyy = 0
     do k = 0, nord - 2
         do j = 0, k
             dxx = dxx + dxxsigma(k, j)*rp**(k - j)*rm**j
             dxy = dxy + dxysigma(k, j)*rp**(k - j)*rm**j
             dyy = dyy + dyysigma(k, j)*rp**(k - j)*rm**j
         enddo
     enddo
     dxx = dxx/two
     dxy = dxy/two
     dyy = dyy/two

     ! if evaluated at zero, one can use instead:
     !dxx = elrad/two*(three*g(3, 0) + two*g(3, 1) + g(3, 2))/two
     !dxy = elrad/two*(0, 1)*(three*g(3, 0) - g(3, 2))/two
     !dyy = -elrad/two*(three*g(3, 0) - two*g(3, 1) + g(3, 2))/two

     te(1, 1, 2) = elrad*kx*(one/h0 + orbit(2)**2/h0**3)/two
     te(1, 1, 4) = elrad*kx*orbit(2)*orbit(4)/h0**3/two
     te(1, 1, 6) = -elrad*kx*orbit(2)/h0**3*(one + beta*orbit(6))/beta/two

     te(1, 2, 1) = te(1, 1, 2)
     te(1, 2, 2) = elrad*(kx*orbit(1) + ky*orbit(3))*3*orbit(2)/h0**3*(one + orbit(2)**2/h0**2)/two
     te(1, 2, 3) = elrad*ky*(one/h0 + orbit(2)**2/h0**3)/two
     te(1, 2, 4) = elrad*(kx*orbit(1) + ky*orbit(3))*(orbit(4)/h0**3 + 3*orbit(2)**2*orbit(4)/h0**5)/two
     te(1, 2, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*(one/h0**2 + 3*orbit(2)**2/h0**4)/beta/h0*(one + beta*orbit(6))/two

     te(1, 3, 2) = te(1, 2, 3)
     te(1, 3, 4) = elrad*ky*orbit(2)*orbit(4)/h0**3/two
     te(1, 3, 6) = elrad*ky*orbit(2)/beta/h0*(one + beta*orbit(6))/two

     te(1, 4, 1) = te(1, 1, 4)
     te(1, 4, 2) = te(1, 2, 4)
     te(1, 4, 3) = te(1, 3, 4)
     te(1, 4, 4) = elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)*(one/h0**3 + 3*orbit(4)**2/h0**5)/two
     te(1, 4, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)*orbit(4)*3/h0**5/beta*(one + beta*orbit(6))/two

     te(1, 6, 1) = te(1, 1, 6)
     te(1, 6, 2) = te(1, 2, 6)
     te(1, 6, 3) = te(1, 3, 6)
     te(1, 6, 4) = te(1, 4, 6)
     te(1, 6, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)/beta/h0**3*(beta - 3/h0/beta**2*(one + orbit(6))**2)/two

     te(2, 1, 1) = real(dxx)   ! cf
     te(2, 1, 3) = real(dxy)   ! cf

     te(2, 2, 2) = - te(1, 1, 2) ! checked
     te(2, 2, 4) = -elrad*kx*orbit(2)*orbit(4)/h0**3/two
     te(2, 2, 6) = elrad*kx*orbit(2)/h0**3/beta*(one + beta*orbit(6))/two

     te(2, 3, 1) = te(2, 1, 3) ! cf
     te(2, 3, 3) = real(dyy)   ! cf

     te(2, 4, 2) = te(2, 2, 4)
     te(2, 4, 4) = -elrad*kx*(one/h0 + orbit(4)**2/h0**3)/two ! - te(1, 1, 2) old
     te(2, 4, 6) = elrad*kx*orbit(4)/h0**3/beta*(one + beta*orbit(6))/two
     
     te(2, 6, 2) = te(2, 2, 6)
     te(2, 6, 4) = te(2, 4, 6)
     te(2, 6, 6) = elrad*kx*(one/h0 - (one + beta*orbit(6))**2/beta**2/h0**3)/two  ! = te(1, 1, 2)*bp1 old

     te(3, 1, 2) = elrad*kx*orbit(2)*orbit(4)/h0**3/two
     te(3, 1, 4) = elrad*kx*(one/h0 + orbit(4)**2/h0**3)/two
     te(3, 1, 6) = -elrad*kx*orbit(4)/h0**3/beta*(one + beta*orbit(6))/two

     te(3, 2, 1) = te(3, 1, 2)

     te(3, 2, 2) = elrad*(kx*orbit(1) + ky*orbit(3))*orbit(4)/h0**3*(one + 3*orbit(2)**2/h0**2)/two
     te(3, 2, 3) = elrad*ky*orbit(2)*orbit(4)/h0**3/two
     te(3, 2, 4) = elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)/h0**3*(one + 3*orbit(4)**2/h0**2)/two
     te(3, 2, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*orbit(2)*orbit(4)/h0**5/beta*(one + beta*orbit(6))/two

     te(3, 3, 2) = te(3, 2, 3)
     te(3, 3, 4) = te(1, 2, 3) ! checked
     te(3, 3, 6) = -elrad*ky*orbit(4)/h0**3/beta*(one + beta*orbit(6))/two

     te(3, 4, 1) = te(3, 1, 4)
     te(3, 4, 2) = te(3, 2, 4)
     te(3, 4, 3) = te(3, 3, 4)

     te(3, 4, 4) = elrad*(kx*orbit(1) + ky*orbit(3))*3*orbit(4)/h0**3*(one + orbit(4)**2/h0**2)/two
     te(3, 4, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))/beta/h0**3*(one + beta*orbit(6))*(one + 3*orbit(4)**2/h0**2)/two

     te(3, 6, 1) = te(3, 1, 6)
     te(3, 6, 2) = te(3, 2, 6)
     te(3, 6, 3) = te(3, 3, 6)
     te(3, 6, 4) = te(3, 4, 6)
     te(3, 6, 6) = -elrad*(kx*orbit(1) + ky*orbit(3))*orbit(4)/beta/h0**3*(beta - 3/beta/h0**2*(one + beta*orbit(6)))/two

     te(4, 1, 1) = - aimag(dxx)  ! cf
     te(4, 1, 3) = - aimag(dxy)  ! cf

     te(4, 2, 2) = elrad*ky/h0*(one + orbit(2)**2/h0**2)/two
     te(4, 2, 4) = elrad*ky*orbit(2)*orbit(4)/h0**3/two
     te(4, 2, 6) = -elrad*ky*orbit(2)/beta/h0**3*(one + beta*orbit(6))/two

     te(4, 3, 1) = te(4, 1, 3)   ! cf
     te(4, 3, 3) = - aimag(dyy)  ! cf

     te(4, 4, 2) = te(4, 2, 4)
     te(4, 4, 4) = -elrad*ky/h0*(one + orbit(4)**2/h0**2)/two
     te(4, 4, 6) = elrad*ky*orbit(4)/beta/h0**3*(one + beta*orbit(6))/two

     te(4, 6, 2) = te(4, 2, 6)
     te(4, 6, 4) = te(4, 4, 6)
     te(4, 6, 6) = elrad*ky/beta/h0*(beta - one/beta/h0**2*(one + beta*orbit(6)))/two

     te(5, 1, 2) = -elrad*kx/beta*(one + beta*orbit(6))*orbit(2)/h0**3/two
     te(5, 1, 4) = -elrad*kx/beta*(one + beta*orbit(6))*orbit(4)/h0**3/two
     te(5, 1, 6) = -elrad*kx/beta/h0*(beta - one/beta/h0**2*(one + beta*orbit(6)))/two

     te(5, 2, 1) = te(5, 1, 2)
     te(5, 2, 2) = -elrad*(kx*orbit(1) + ky*orbit(3))*(one + beta*orbit(6))/beta/h0**3*(one + 3*orbit(2)**2/h0**2)/two
     te(5, 2, 3) = -elrad*orbit(2)*ky*(one + beta*orbit(6))/beta/h0**3/two
     te(5, 2, 4) = -elrad*orbit(2)*(kx*orbit(1) + ky*orbit(3))*(one + beta*orbit(6))/beta*3/h0**5*orbit(4)/two
     te(5, 2, 6) = -elrad*orbit(2)*(kx*orbit(1) + ky*orbit(3))/beta/h0**3*(beta - 3/beta/h0**2*(one + beta*orbit(6))**2)/two

     ! 5, 3 = ky/kx * 5, 1
     te(5, 3, 2) = te(5, 2, 3)
     te(5, 3, 4) = -elrad*ky/beta*(one + beta*orbit(6))*orbit(4)/h0**3/two
     te(5, 3, 6) = -elrad*ky/beta/h0*(beta - one/beta/h0**2*(one + beta*orbit(6)))/two

     te(5, 4, 1) = te(5, 1, 4)
     te(5, 4, 2) = te(5, 2, 4)
     te(5, 4, 3) = te(5, 3, 4)
     te(5, 4, 4) = -elrad*(kx*orbit(1) + ky*orbit(3))*(one + beta*orbit(6))/beta/h0**3*(one + 3*orbit(4)**2/h0**2)/two
     te(5, 4, 6) = -elrad*orbit(4)*(kx*orbit(1) + ky*orbit(3))/beta/h0**3*(beta - 3/beta/h0**2*(one + beta*orbit(6))**2)/two

     te(5, 6, 1) = te(5, 1, 6)
     te(5, 6, 2) = te(5, 2, 6)
     te(5, 6, 3) = te(5, 3, 6)
     te(5, 6, 4) = te(5, 4, 6)
     te(5, 6, 6) = elrad*3*(kx*orbit(1) + ky*orbit(3))/beta**2/h0**3*(one + beta*orbit(6))*&
                   (beta - one/beta/h0**2*(one + beta*orbit(6))**2)/two
  endif
end SUBROUTINE tmmult_cf

SUBROUTINE tmmult(fsec,ftrk,orbit,fmap,re,te)
  use twtrrfi
  use twisslfi
  use twissbeamfi, only : radiate, deltap, beta, gamma, arad
  use math_constfi, only : zero, one, two, three
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for thin multipole.                                *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  double precision :: orbit(6), re(6,6), te(6,6,6)

  integer :: n_ferr, nord, iord, j, nd, nn, ns
  double precision :: f_errors(0:maxferr)
  double precision :: normal(0:maxmul), skew(0:maxmul)
  double precision :: bi, pt, rfac, bvk, elrad, tilt, angle, an, anr, ani
  double precision :: x, y, dbr, dbi, dipr, dipi, dr, di, drt, dpx, dpy, dpxr, dpyr, dtmp

  integer, external :: get_option, node_fd_errors
  double precision, external :: node_value, get_value
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !---- Initialize
  rfac = zero
  F_ERRORS(0:maxferr) = zero
  n_ferr = node_fd_errors(f_errors)
  bvk = node_value('other_bv ')

  !---- Multipole length for radiation.
  elrad = node_value('lrad ')

  bi = one / beta

  fmap = .true.

  !---- Multipole components.
  NORMAL = zero ; call get_node_vector('knl ',nn,normal)
  SKEW   = zero ; call get_node_vector('ksl ',ns,skew)
  tilt = node_value('tilt ')

  nd = 2 * max(nn, ns, n_ferr/2-1)

  !---- Angle (bvk applied later)
  an = node_value('angle ')
  if (an .ne. 0) then 
    anr = an
    f_errors(0) = f_errors(0) + normal(0) - an
  endif
  !---- Dipole error.
  dbr = f_errors(0) / (one + deltap)
  dbi = f_errors(1) / (one + deltap)

  !---- Nominal dipole strength.
  dipr = normal(0) / (one + deltap)
  dipi = skew(0)   / (one + deltap)

  if (tilt .ne. zero)  then
     if (dipi.ne.zero .or. dipr.ne.zero) then
        angle = atan2(dipi, dipr) - tilt
     else
        angle = -tilt
     endif
     dtmp = sqrt(dipi**2 + dipr**2)
     dipr = dtmp * cos(angle)
     dipi = dtmp * sin(angle)
     dtmp = sqrt(dbi**2 + dbr**2)
     dbr = dtmp * cos(angle)
     dbi = dtmp * sin(angle)
     anr = an * cos(angle)
     ani = an * sin(angle)
     anr   = bvk * anr
     ani   = bvk * ani
  endif


  dbr  = bvk * dbr
  dbi  = bvk * dbi
  dipr = bvk * dipr
  dipi = bvk * dipi
  !---- Other components and errors.
  nord = 0
  ! that loop should start at one since nominal dipole strength already taken into account above
  do iord = 0, max(nn, ns, n_ferr/2-1)
     ! get the maximum effective order; loop runs over maximum of user given values
     if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero .or. &
          normal(iord).ne.zero .or. skew(iord).ne.zero) nord = iord
  enddo

  do iord = 1, nord
     f_errors(2*iord)   = (normal(iord) + f_errors(2*iord))   / (one + deltap)
     f_errors(2*iord+1) = (skew(iord)   + f_errors(2*iord+1)) / (one + deltap)
     if (tilt .ne. zero) then
        if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero) then
           angle = atan2(f_errors(2*iord+1), f_errors(2*iord)) / (iord+1) - tilt
        else
           angle = -tilt
        endif
        angle = (iord+1) * angle
        dtmp = sqrt(f_errors(2*iord)**2 + f_errors(2*iord+1)**2)
        f_errors(2*iord)   = dtmp * cos(angle)
        f_errors(2*iord+1) = dtmp * sin(angle)
     endif
     f_errors(2*iord)   = bvk * f_errors(2*iord)
     f_errors(2*iord+1) = bvk * f_errors(2*iord+1)
  enddo

  if (ftrk) then !---- Track orbit.
     x = orbit(1)
     y = orbit(3)

     !---- Multipole kick.
     dr = zero
     di = zero

     do iord = nord, 1, -1
        drt = (dr * x - di * y) / (iord+1) + f_errors(2*iord)
        di  = (dr * y + di * x) / (iord+1) + f_errors(2*iord+1)
        dr  = drt
     enddo
     dpx = dbr + (dr * x - di * y)
     dpy = dbi + (di * x + dr * y)


     !---- Radiation effects at entrance.
     if (radiate  .and.  elrad.ne.zero) then
        dpxr = dpx + dipr
        dpyr = dpy + dipi
        rfac = arad * gamma**3 * (dpxr**2+dpyr**2) / (three*elrad)
        pt = orbit(6)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif

     !---- Track orbit.
     orbit(2) = orbit(2) - dpx + dipr * (deltap + bi*orbit(6))
     orbit(4) = orbit(4) + dpy - dipi * (deltap + bi*orbit(6))
     orbit(5) = orbit(5) - (dipr*x + dipi*y) * bi

     !---- Add the missing focussing component of thin dipoles for co
     if (elrad.gt.zero .and. get_option('thin_foc ').eq.1) then
        if (an .ne. 0) then
          orbit(2) = orbit(2) - anr*dipr/elrad * x ! 
          orbit(4) = orbit(4) - ani*dipi/elrad * y
        else
          orbit(2) = orbit(2) - (one+deltap)*dipr*dipr/elrad * x ! 
          orbit(4) = orbit(4) - (one+deltap)*dipi*dipi/elrad * y
        endif
     endif
     !---- Radiation effects at exit.
     if (radiate  .and.  elrad.ne.zero) then
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif

  else !---- Orbit not wanted.
     x = zero
     y = zero
     nord = min(nord, 2)
  endif

  !---- First-order terms (use X,Y from orbit tracking).
  if (nord .ge. 1) then
     dr = zero
     di = zero
     do iord = nord, 1, -1
        drt = (dr * x - di * y) / (iord) + f_errors(2*iord)
        di  = (dr * y + di * x) / (iord) + f_errors(2*iord+1)
        dr  = drt
     enddo
     re(2,1) = - dr
     re(2,3) = + di
     re(4,1) = + di
     re(4,3) = + dr
  endif

  !---- Add the missing focussing component of thin dipoles
  !---- The (1+deltap) is from that the term is h*k0 (so one geometrical and one is bending strength)
  if (elrad.gt.zero.and.get_option('thin_foc ').eq.1) then
    if (an .ne. 0) then
      re(2,1) = re(2,1) - anr*dipr/elrad
      re(4,3) = re(4,3) - ani*dipi/elrad
    else
      re(2,1) = re(2,1) - (one+deltap)*dipr*dipr/elrad
      re(4,3) = re(4,3) - (one+deltap)*dipi*dipi/elrad
    endif

  endif
  re(2,6) = + dipr * bi
  re(4,6) = - dipi * bi
  re(5,1) = - re(2,6)
  re(5,3) = - re(4,6)

  !---- Second-order terms (use X,Y from orbit tracking).

  if (fsec) then
     if (nord .ge. 2) then
        dr = zero
        di = zero
        do iord = nord, 2, -1
           drt = (dr * x - di * y) / (iord-1) + f_errors(2*iord)
           di  = (dr * y + di * x) / (iord-1) + f_errors(2*iord+1)
           dr  = drt
        enddo
        dr = dr / two
        di = di / two
        te(2,1,1) = - dr
        te(2,1,3) = + di
        te(2,3,1) = + di
        te(2,3,3) = + dr
        te(4,1,1) = + di
        te(4,1,3) = + dr
        te(4,3,1) = + dr
        te(4,3,3) = - di
     endif
  endif
end SUBROUTINE tmmult

SUBROUTINE tmoct(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)
  use twtrrfi
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : radiate, deltap, gamma, arad
  use matrices, only : EYE
  use math_constfi, only : zero, one, two, three, four, six
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for octupole element.                              *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is NOT rotated before and after this routine   *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical, intent(IN) :: fsec, ftrk, fcentre
  logical, intent(OUT) :: fmap
  double precision, intent(IN) :: el, dl
  double precision, intent(IN OUT) :: orbit(6)
  double precision, intent(OUT) :: ek(6), re(6,6), te(6,6,6)

  logical :: cplxy
  integer :: i, n_ferr, elpar_vl
  double precision :: f_errors(0:maxferr)
  double precision :: rw(6,6), tw(6,6,6)
  double precision :: sk3, sk3l, sk3s, octr, octi, posr, posi, cr, ci
  double precision :: rfac, pt, bvk, tilt4

  integer, external :: node_fd_errors, el_par_vector
  double precision, external :: node_value, get_value
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !---Initializasion
  rfac=0.d0

  if ( .not. ftrk) then
     !---- No orbit track requested, use drift map and return
     call tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)
     return
  endif

  !---- Initialize.
  fmap = el .ne. zero
  if (.not. fmap) return

  !-- get element parameters
  bvk = node_value('other_bv ')
  elpar_vl = el_par_vector(o_k3s, g_elpar)

  !---- Field error.
  F_ERRORS = zero
  n_ferr = node_fd_errors(f_errors)

  !---- Set up octupole strength.
  sk3  = bvk * ( g_elpar(o_k3)  + f_errors(6)/el)
  sk3s = bvk * ( g_elpar(o_k3s) + f_errors(7)/el)
  tilt4 = -four * node_value('tilt ')

  if (sk3s.ne.zero) then
     tilt4 = atan2(sk3s, sk3) + tilt4
     sk3 = sqrt(sk3**2 + sk3s**2)
  endif
  sk3l = (el * sk3) / (one + deltap)

  !---- Normal and skew components of octupole.
  octr = sk3l * cos(tilt4)
  octi = sk3l * sin(tilt4)

  !---- Half kick at entrance.
  posr = orbit(1) * (orbit(1)**2 - three*orbit(3)**2) / six
  posi = orbit(3) * (three*orbit(1)**2 - orbit(3)**2) / six
  cr = octr * posr - octi * posi
  ci = octr * posi + octi * posr
  orbit(2) = orbit(2) - cr / two
  orbit(4) = orbit(4) + ci / two

  !---- Half radiation effects at entrance.
  if (radiate) then
     rfac = arad * gamma**3 * (cr**2 + ci**2) / (three * el)
     pt = orbit(6)
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

  !---- First-order terms w.r.t. orbit.
  RW = EYE
  posr = (orbit(1)**2 - orbit(3)**2) / four
  posi = orbit(1) * orbit(3) / two
  cr = octr * posr - octi * posi
  ci = octr * posi + octi * posr

  rw(2,1) = - cr
  rw(2,3) = + ci
  rw(4,1) = + ci
  rw(4,3) = + cr

  if (ci .ne. zero) cplxy = .true.

  TW = zero
  !---- Second-order terms w.r.t. orbit.
  if (fsec) then
     cr = (octr * orbit(1) - octi * orbit(3)) / four
     ci = (octr * orbit(3) + octi * orbit(1)) / four

     tw(2,1,1) = - cr
     tw(2,1,3) = + ci
     tw(2,3,1) = + ci
     tw(2,3,3) = + cr

     tw(4,1,1) = + ci
     tw(4,1,3) = + cr
     tw(4,3,1) = + cr
     tw(4,3,3) = - ci
  endif

  !---- Concatenate with drift map.
  call tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)
  call tmcat(fsec,re,te,rw,tw,re,te)
  if (fcentre) return

  !---- Half kick at exit.
  posr = orbit(1) * (orbit(1)**2 - three*orbit(3)**2) / six
  posi = orbit(3) * (three*orbit(1)**2 - orbit(3)**2) / six
  cr = octr * posr - octi * posi
  ci = octr * posi + octi * posr
  orbit(2) = orbit(2) - cr / two
  orbit(4) = orbit(4) + ci / two

  !---- Half radiation effects.
  if (radiate) then
     rfac = arad * gamma**3 * (cr**2 + ci**2) / (three * el)
     pt = orbit(6)
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

  !---- First-order terms w.r.t. orbit.
  RW = EYE
  posr = (orbit(1)**2 - orbit(3)**2) / four
  posi = orbit(1) * orbit(3) / two
  cr = octr * posr - octi * posi
  ci = octr * posi + octi * posr

  rw(2,1) = - cr
  rw(2,3) = + ci
  rw(4,1) = + ci
  rw(4,3) = + cr

  if (ci .ne. zero) cplxy = .true.

  !---- Second-order terms w.r.t. orbit.
  TW = zero
  if (fsec) then
     cr = (octr * orbit(1) - octi * orbit(3)) / four
     ci = (octr * orbit(3) + octi * orbit(1)) / four

     tw(2,1,1) = - cr
     tw(2,1,3) = + ci
     tw(2,3,1) = + ci
     tw(2,3,3) = + cr

     tw(4,1,1) = + ci
     tw(4,1,3) = + cr
     tw(4,3,1) = + cr
     tw(4,3,3) = - ci
  endif
  call tmcat(fsec,rw,tw,re,te,re,te)

end SUBROUTINE tmoct

SUBROUTINE tmarb(fsec,ftrk,orbit,fmap,ek,re,te)
  use trackfi, only : fsecarb
  use twtrrfi
  use name_lenfi
  use time_varfi
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for arbitrary element.                             *
  !     Input:                                                           *
  !     FSEC      (logical) If true, return second order terms.          *
  !     FTRK      (logical) If true, track orbit.                        *
  !     Input/output:                                                    *
  !     ORBIT(6)  (real)    Closed orbit.                                *
  !     Output:                                                          *
  !     FMAP      (logical) If true, element has a map.                  *
  !     EK(6)     (real)    Kick due to element.                         *
  !     RE(6,6)   (real)    Transfer matrix.                             *
  !     TE(6,6,6) (real)    Second-order terms.                          *
  !----------------------------------------------------------------------*
  logical, intent(IN) :: fsec, ftrk
  logical, intent(OUT) :: fmap
  double precision, intent(IN OUT) :: orbit(6)
  double precision, intent(OUT) :: ek(6), re(6,6), te(6,6,6)

  logical :: time_var
  integer :: i1, i2, i3, mylen
  double precision :: temp
  character(len=name_len) :: name

  double precision, external :: node_value

  !---- Element Kick
  ek(1)=node_value('kick1 ')
  ek(2)=node_value('kick2 ')
  ek(3)=node_value('kick3 ')
  ek(4)=node_value('kick4 ')
  ek(5)=node_value('kick5 ')
  ek(6)=node_value('kick6 ')

  !---- Matrix
  re(1,1)=node_value('rm11 ')
  re(1,2)=node_value('rm12 ')
  re(1,3)=node_value('rm13 ')
  re(1,4)=node_value('rm14 ')
  re(1,5)=node_value('rm15 ')
  re(1,6)=node_value('rm16 ')
  re(2,1)=node_value('rm21 ')
  re(2,2)=node_value('rm22 ')
  re(2,3)=node_value('rm23 ')
  re(2,4)=node_value('rm24 ')
  re(2,5)=node_value('rm25 ')
  re(2,6)=node_value('rm26 ')
  re(3,1)=node_value('rm31 ')
  re(3,2)=node_value('rm32 ')
  re(3,3)=node_value('rm33 ')
  re(3,4)=node_value('rm34 ')
  re(3,5)=node_value('rm35 ')
  re(3,6)=node_value('rm36 ')
  re(4,1)=node_value('rm41 ')
  re(4,2)=node_value('rm42 ')
  re(4,3)=node_value('rm43 ')
  re(4,4)=node_value('rm44 ')
  re(4,5)=node_value('rm45 ')
  re(4,6)=node_value('rm46 ')
  re(5,1)=node_value('rm51 ')
  re(5,2)=node_value('rm52 ')
  re(5,3)=node_value('rm53 ')
  re(5,4)=node_value('rm54 ')
  re(5,5)=node_value('rm55 ')
  re(5,6)=node_value('rm56 ')
  re(6,1)=node_value('rm61 ')
  re(6,2)=node_value('rm62 ')
  re(6,3)=node_value('rm63 ')
  re(6,4)=node_value('rm64 ')
  re(6,5)=node_value('rm65 ')
  re(6,6)=node_value('rm66 ')

  time_var = node_value('time_var ') .ne. zero

  timevar:    if (time_var.and.time_var_p) then
     time_var_p_cnt = time_var_p_cnt+1
     time_var_p_lnt = time_var_p_lnt+1

     if (idnint(time_var_p_ind(time_var_p_cnt)) .ne. time_var_p_lnt)    &
          call fort_fail('TMARB: ', 'wrong index in Table: time_var_pha')
     call element_name(name,len(name))
     mylen=len_trim(name)
     if (time_var_p_ch(time_var_p_cnt)(:mylen) .ne. name(:mylen))       &
          call fort_fail('TMARB: ', 'wrong element name in Table: time_var_pha')

!---- Matrix
     re(1,1)=phase_tromb(time_var_p_cnt,1)
     re(1,2)=phase_tromb(time_var_p_cnt,2)
     re(1,3)=phase_tromb(time_var_p_cnt,3)
     re(1,4)=phase_tromb(time_var_p_cnt,4)
     re(1,5)=phase_tromb(time_var_p_cnt,5)
     re(1,6)=phase_tromb(time_var_p_cnt,6)
     re(2,1)=phase_tromb(time_var_p_cnt,7)
     re(2,2)=phase_tromb(time_var_p_cnt,8)
     re(2,3)=phase_tromb(time_var_p_cnt,9)
     re(2,4)=phase_tromb(time_var_p_cnt,10)
     re(2,5)=phase_tromb(time_var_p_cnt,11)
     re(2,6)=phase_tromb(time_var_p_cnt,12)
     re(3,1)=phase_tromb(time_var_p_cnt,13)
     re(3,2)=phase_tromb(time_var_p_cnt,14)
     re(3,3)=phase_tromb(time_var_p_cnt,15)
     re(3,4)=phase_tromb(time_var_p_cnt,16)
     re(3,5)=phase_tromb(time_var_p_cnt,17)
     re(3,6)=phase_tromb(time_var_p_cnt,18)
     re(4,1)=phase_tromb(time_var_p_cnt,19)
     re(4,2)=phase_tromb(time_var_p_cnt,20)
     re(4,3)=phase_tromb(time_var_p_cnt,21)
     re(4,4)=phase_tromb(time_var_p_cnt,22)
     re(4,5)=phase_tromb(time_var_p_cnt,23)
     re(4,6)=phase_tromb(time_var_p_cnt,24)
     re(5,1)=phase_tromb(time_var_p_cnt,25)
     re(5,2)=phase_tromb(time_var_p_cnt,26)
     re(5,3)=phase_tromb(time_var_p_cnt,27)
     re(5,4)=phase_tromb(time_var_p_cnt,28)
     re(5,5)=phase_tromb(time_var_p_cnt,29)
     re(5,6)=phase_tromb(time_var_p_cnt,30)
     re(6,1)=phase_tromb(time_var_p_cnt,31)
     re(6,2)=phase_tromb(time_var_p_cnt,32)
     re(6,3)=phase_tromb(time_var_p_cnt,33)
     re(6,4)=phase_tromb(time_var_p_cnt,34)
     re(6,5)=phase_tromb(time_var_p_cnt,35)
     re(6,6)=phase_tromb(time_var_p_cnt,36)

!---- Matrix
     call store_node_value('rm11 ',re(1,1))
     call store_node_value('rm12 ',re(1,2))
     call store_node_value('rm13 ',re(1,3))
     call store_node_value('rm14 ',re(1,4))
     call store_node_value('rm15 ',re(1,5))
     call store_node_value('rm16 ',re(1,6))
     call store_node_value('rm21 ',re(2,1))
     call store_node_value('rm22 ',re(2,2))
     call store_node_value('rm23 ',re(2,3))
     call store_node_value('rm24 ',re(2,4))
     call store_node_value('rm25 ',re(2,5))
     call store_node_value('rm26 ',re(2,6))
     call store_node_value('rm31 ',re(3,1))
     call store_node_value('rm32 ',re(3,2))
     call store_node_value('rm33 ',re(3,3))
     call store_node_value('rm34 ',re(3,4))
     call store_node_value('rm35 ',re(3,5))
     call store_node_value('rm36 ',re(3,6))
     call store_node_value('rm41 ',re(4,1))
     call store_node_value('rm42 ',re(4,2))
     call store_node_value('rm43 ',re(4,3))
     call store_node_value('rm44 ',re(4,4))
     call store_node_value('rm45 ',re(4,5))
     call store_node_value('rm46 ',re(4,6))
     call store_node_value('rm51 ',re(5,1))
     call store_node_value('rm52 ',re(5,2))
     call store_node_value('rm53 ',re(5,3))
     call store_node_value('rm54 ',re(5,4))
     call store_node_value('rm55 ',re(5,5))
     call store_node_value('rm56 ',re(5,6))
     call store_node_value('rm61 ',re(6,1))
     call store_node_value('rm62 ',re(6,2))
     call store_node_value('rm63 ',re(6,3))
     call store_node_value('rm64 ',re(6,4))
     call store_node_value('rm65 ',re(6,5))
     call store_node_value('rm66 ',re(6,6))
  endif timevar

  !---- Second Order Terms
  if (fsec) then
     te(1,1,1)=node_value('tm111 ')
     te(1,1,2)=node_value('tm112 ')
     te(1,1,3)=node_value('tm113 ')
     te(1,1,4)=node_value('tm114 ')
     te(1,1,5)=node_value('tm115 ')
     te(1,1,6)=node_value('tm116 ')
     te(1,2,1)=node_value('tm121 ')
     te(1,2,2)=node_value('tm122 ')
     te(1,2,3)=node_value('tm123 ')
     te(1,2,4)=node_value('tm124 ')
     te(1,2,5)=node_value('tm125 ')
     te(1,2,6)=node_value('tm126 ')
     te(1,3,1)=node_value('tm131 ')
     te(1,3,2)=node_value('tm132 ')
     te(1,3,3)=node_value('tm133 ')
     te(1,3,4)=node_value('tm134 ')
     te(1,3,5)=node_value('tm135 ')
     te(1,3,6)=node_value('tm136 ')
     te(1,4,1)=node_value('tm141 ')
     te(1,4,2)=node_value('tm142 ')
     te(1,4,3)=node_value('tm143 ')
     te(1,4,4)=node_value('tm144 ')
     te(1,4,5)=node_value('tm145 ')
     te(1,4,6)=node_value('tm146 ')
     te(1,5,1)=node_value('tm151 ')
     te(1,5,2)=node_value('tm152 ')
     te(1,5,3)=node_value('tm153 ')
     te(1,5,4)=node_value('tm154 ')
     te(1,5,5)=node_value('tm155 ')
     te(1,5,6)=node_value('tm156 ')
     te(1,6,1)=node_value('tm161 ')
     te(1,6,2)=node_value('tm162 ')
     te(1,6,3)=node_value('tm163 ')
     te(1,6,4)=node_value('tm164 ')
     te(1,6,5)=node_value('tm165 ')
     te(1,6,6)=node_value('tm166 ')
     te(2,1,1)=node_value('tm211 ')
     te(2,1,2)=node_value('tm212 ')
     te(2,1,3)=node_value('tm213 ')
     te(2,1,4)=node_value('tm214 ')
     te(2,1,5)=node_value('tm215 ')
     te(2,1,6)=node_value('tm216 ')
     te(2,2,1)=node_value('tm221 ')
     te(2,2,2)=node_value('tm222 ')
     te(2,2,3)=node_value('tm223 ')
     te(2,2,4)=node_value('tm224 ')
     te(2,2,5)=node_value('tm225 ')
     te(2,2,6)=node_value('tm226 ')
     te(2,3,1)=node_value('tm231 ')
     te(2,3,2)=node_value('tm232 ')
     te(2,3,3)=node_value('tm233 ')
     te(2,3,4)=node_value('tm234 ')
     te(2,3,5)=node_value('tm235 ')
     te(2,3,6)=node_value('tm236 ')
     te(2,4,1)=node_value('tm241 ')
     te(2,4,2)=node_value('tm242 ')
     te(2,4,3)=node_value('tm243 ')
     te(2,4,4)=node_value('tm244 ')
     te(2,4,5)=node_value('tm245 ')
     te(2,4,6)=node_value('tm246 ')
     te(2,5,1)=node_value('tm251 ')
     te(2,5,2)=node_value('tm252 ')
     te(2,5,3)=node_value('tm253 ')
     te(2,5,4)=node_value('tm254 ')
     te(2,5,5)=node_value('tm255 ')
     te(2,5,6)=node_value('tm256 ')
     te(2,6,1)=node_value('tm261 ')
     te(2,6,2)=node_value('tm262 ')
     te(2,6,3)=node_value('tm263 ')
     te(2,6,4)=node_value('tm264 ')
     te(2,6,5)=node_value('tm265 ')
     te(2,6,6)=node_value('tm266 ')
     te(3,1,1)=node_value('tm311 ')
     te(3,1,2)=node_value('tm312 ')
     te(3,1,3)=node_value('tm313 ')
     te(3,1,4)=node_value('tm314 ')
     te(3,1,5)=node_value('tm315 ')
     te(3,1,6)=node_value('tm316 ')
     te(3,2,1)=node_value('tm321 ')
     te(3,2,2)=node_value('tm322 ')
     te(3,2,3)=node_value('tm323 ')
     te(3,2,4)=node_value('tm324 ')
     te(3,2,5)=node_value('tm325 ')
     te(3,2,6)=node_value('tm326 ')
     te(3,3,1)=node_value('tm331 ')
     te(3,3,2)=node_value('tm332 ')
     te(3,3,3)=node_value('tm333 ')
     te(3,3,4)=node_value('tm334 ')
     te(3,3,5)=node_value('tm335 ')
     te(3,3,6)=node_value('tm336 ')
     te(3,4,1)=node_value('tm341 ')
     te(3,4,2)=node_value('tm342 ')
     te(3,4,3)=node_value('tm343 ')
     te(3,4,4)=node_value('tm344 ')
     te(3,4,5)=node_value('tm345 ')
     te(3,4,6)=node_value('tm346 ')
     te(3,5,1)=node_value('tm351 ')
     te(3,5,2)=node_value('tm352 ')
     te(3,5,3)=node_value('tm353 ')
     te(3,5,4)=node_value('tm354 ')
     te(3,5,5)=node_value('tm355 ')
     te(3,5,6)=node_value('tm356 ')
     te(3,6,1)=node_value('tm361 ')
     te(3,6,2)=node_value('tm362 ')
     te(3,6,3)=node_value('tm363 ')
     te(3,6,4)=node_value('tm364 ')
     te(3,6,5)=node_value('tm365 ')
     te(3,6,6)=node_value('tm366 ')
     te(4,1,1)=node_value('tm411 ')
     te(4,1,2)=node_value('tm412 ')
     te(4,1,3)=node_value('tm413 ')
     te(4,1,4)=node_value('tm414 ')
     te(4,1,5)=node_value('tm415 ')
     te(4,1,6)=node_value('tm416 ')
     te(4,2,1)=node_value('tm421 ')
     te(4,2,2)=node_value('tm422 ')
     te(4,2,3)=node_value('tm423 ')
     te(4,2,4)=node_value('tm424 ')
     te(4,2,5)=node_value('tm425 ')
     te(4,2,6)=node_value('tm426 ')
     te(4,3,1)=node_value('tm431 ')
     te(4,3,2)=node_value('tm432 ')
     te(4,3,3)=node_value('tm433 ')
     te(4,3,4)=node_value('tm434 ')
     te(4,3,5)=node_value('tm435 ')
     te(4,3,6)=node_value('tm436 ')
     te(4,4,1)=node_value('tm441 ')
     te(4,4,2)=node_value('tm442 ')
     te(4,4,3)=node_value('tm443 ')
     te(4,4,4)=node_value('tm444 ')
     te(4,4,5)=node_value('tm445 ')
     te(4,4,6)=node_value('tm446 ')
     te(4,5,1)=node_value('tm451 ')
     te(4,5,2)=node_value('tm452 ')
     te(4,5,3)=node_value('tm453 ')
     te(4,5,4)=node_value('tm454 ')
     te(4,5,5)=node_value('tm455 ')
     te(4,5,6)=node_value('tm456 ')
     te(4,6,1)=node_value('tm461 ')
     te(4,6,2)=node_value('tm462 ')
     te(4,6,3)=node_value('tm463 ')
     te(4,6,4)=node_value('tm464 ')
     te(4,6,5)=node_value('tm465 ')
     te(4,6,6)=node_value('tm466 ')
     te(5,1,1)=node_value('tm511 ')
     te(5,1,2)=node_value('tm512 ')
     te(5,1,3)=node_value('tm513 ')
     te(5,1,4)=node_value('tm514 ')
     te(5,1,5)=node_value('tm515 ')
     te(5,1,6)=node_value('tm516 ')
     te(5,2,1)=node_value('tm521 ')
     te(5,2,2)=node_value('tm522 ')
     te(5,2,3)=node_value('tm523 ')
     te(5,2,4)=node_value('tm524 ')
     te(5,2,5)=node_value('tm525 ')
     te(5,2,6)=node_value('tm526 ')
     te(5,3,1)=node_value('tm531 ')
     te(5,3,2)=node_value('tm532 ')
     te(5,3,3)=node_value('tm533 ')
     te(5,3,4)=node_value('tm534 ')
     te(5,3,5)=node_value('tm535 ')
     te(5,3,6)=node_value('tm536 ')
     te(5,4,1)=node_value('tm541 ')
     te(5,4,2)=node_value('tm542 ')
     te(5,4,3)=node_value('tm543 ')
     te(5,4,4)=node_value('tm544 ')
     te(5,4,5)=node_value('tm545 ')
     te(5,4,6)=node_value('tm546 ')
     te(5,5,1)=node_value('tm551 ')
     te(5,5,2)=node_value('tm552 ')
     te(5,5,3)=node_value('tm553 ')
     te(5,5,4)=node_value('tm554 ')
     te(5,5,5)=node_value('tm555 ')
     te(5,5,6)=node_value('tm556 ')
     te(5,6,1)=node_value('tm561 ')
     te(5,6,2)=node_value('tm562 ')
     te(5,6,3)=node_value('tm563 ')
     te(5,6,4)=node_value('tm564 ')
     te(5,6,5)=node_value('tm565 ')
     te(5,6,6)=node_value('tm566 ')
     te(6,1,1)=node_value('tm611 ')
     te(6,1,2)=node_value('tm612 ')
     te(6,1,3)=node_value('tm613 ')
     te(6,1,4)=node_value('tm614 ')
     te(6,1,5)=node_value('tm615 ')
     te(6,1,6)=node_value('tm616 ')
     te(6,2,1)=node_value('tm621 ')
     te(6,2,2)=node_value('tm622 ')
     te(6,2,3)=node_value('tm623 ')
     te(6,2,4)=node_value('tm624 ')
     te(6,2,5)=node_value('tm625 ')
     te(6,2,6)=node_value('tm626 ')
     te(6,3,1)=node_value('tm631 ')
     te(6,3,2)=node_value('tm632 ')
     te(6,3,3)=node_value('tm633 ')
     te(6,3,4)=node_value('tm634 ')
     te(6,3,5)=node_value('tm635 ')
     te(6,3,6)=node_value('tm636 ')
     te(6,4,1)=node_value('tm641 ')
     te(6,4,2)=node_value('tm642 ')
     te(6,4,3)=node_value('tm643 ')
     te(6,4,4)=node_value('tm644 ')
     te(6,4,5)=node_value('tm645 ')
     te(6,4,6)=node_value('tm646 ')
     te(6,5,1)=node_value('tm651 ')
     te(6,5,2)=node_value('tm652 ')
     te(6,5,3)=node_value('tm653 ')
     te(6,5,4)=node_value('tm654 ')
     te(6,5,5)=node_value('tm655 ')
     te(6,5,6)=node_value('tm656 ')
     te(6,6,1)=node_value('tm661 ')
     te(6,6,2)=node_value('tm662 ')
     te(6,6,3)=node_value('tm663 ')
     te(6,6,4)=node_value('tm664 ')
     te(6,6,5)=node_value('tm665 ')
     te(6,6,6)=node_value('tm666 ')

     if (maxval(abs(TE)) .gt. zero) fsecarb = .true.

  endif

  !---- Track orbit.
  fmap = .true.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmarb

SUBROUTINE tmquad(fsec,ftrk,fcentre,plot_tilt,orbit,fmap,el,dl,ek,re,te)
  use twtrrfi
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : radiate, deltap, gamma, arad
  use math_constfi, only : zero, one, two, three
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for quadrupole element.                            *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     plot_tilt (double)  external tilt needed for plot                *
  !     el        (double)  element length.                              *
  !     dl        (double)  slice length.                                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is rotated (in tmmap) before and after this    *
  !     routine                                                          *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap, fcentre
  double precision :: plot_tilt, el, dl
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  logical :: cplxy
  integer :: i, j, n_ferr, elpar_vl
  double precision :: ct, st, tmp
  double precision :: f_errors(0:maxferr)
  double precision :: tilt, sk1, pt, sk1s, bvk, rfac

  integer, external :: node_fd_errors
  integer, external :: el_par_vector
  double precision, external :: node_value, get_value
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !---- Initialize.
  st = zero
  ct = zero
  cplxy = .false.
  fmap = el .ne. zero
  if (.not. fmap) return

  !---- Field error.
  F_ERRORS = zero
  n_ferr = node_fd_errors(f_errors)

  !-- element paramters
  elpar_vl = el_par_vector(q_k1st, g_elpar)
  bvk = node_value('other_bv ')
  sk1  = bvk * ( g_elpar(q_k1)  + g_elpar(q_k1t)  + f_errors(2)/el)
  sk1s = bvk * ( g_elpar(q_k1s) + g_elpar(q_k1st) + f_errors(3)/el)
  tilt = g_elpar(q_tilt)
  if (sk1s .ne. zero) then
     tilt = -atan2(sk1s, sk1)/two + tilt
     sk1 = sqrt(sk1**2 + sk1s**2)
  endif
  if (tilt .ne. zero)  then
     !---  rotate orbit before entry
     st = sin(tilt)
     ct = cos(tilt)
     tmp = orbit(1)
     orbit(1) = ct * tmp + st * orbit(3)
     orbit(3) = ct * orbit(3) - st * tmp
     tmp = orbit(2)
     orbit(2) = ct * tmp + st * orbit(4)
     orbit(4) = ct * orbit(4) - st * tmp
  endif

  cplxy = cplxy .or. tilt .ne. zero
  tilt = tilt + plot_tilt

  sk1 = sk1 / (one + deltap)

  !---- Half radiation effect at entry.
  if (radiate .and. ftrk) then
     rfac = (arad * gamma**3 * sk1**2 * el / three) * (orbit(1)**2 + orbit(3)**2)
     pt = orbit(6)
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

  call qdbody(fsec,ftrk,tilt,sk1,orbit,dl,ek,re,te)
  if (fcentre) return

  !---- Half radiation effect at exit.
  if (radiate .and. ftrk) then
     rfac = (arad * gamma**3 * sk1**2 * el / three) * (orbit(1)**2 + orbit(3)**2)
     pt = orbit(6)
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

  if (tilt .ne. zero)  then
     !---  rotate orbit at exit
     tmp = orbit(1)
     orbit(1) = ct * tmp - st * orbit(3)
     orbit(3) = ct * orbit(3) + st * tmp
     tmp = orbit(2)
     orbit(2) = ct * tmp - st * orbit(4)
     orbit(4) = ct * orbit(4) + st * tmp
  endif

end SUBROUTINE tmquad

SUBROUTINE qdbody(fsec,ftrk,tilt,sk1,orbit,el,ek,re,te)
  use twissbeamfi, only : beta, gamma, dtbyds
  use math_constfi, only : zero, one, two, four, six, ten3m
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for quadrupole element.                            *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     tilt      (double)  map tilt angle                               *
  !     sk1       (double)  processed quad kick                          *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk
  double precision :: tilt, sk1, el
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: qk, qkl, qkl2
  double precision :: cx, sx, cy, sy, biby4

  !---- Set up c's and s's.
  qk = sqrt(abs(sk1))
  qkl = qk * el
  if (abs(qkl) .lt. ten3m) then
     qkl2 = sk1 * el**2
     cx = (one - qkl2 / two)
     sx = (one - qkl2 / six) * el
     cy = (one + qkl2 / two)
     sy = (one + qkl2 / six) * el
  else if (sk1 .gt. zero) then
     cx = cos(qkl)
     sx = sin(qkl) / qk
     cy = cosh(qkl)
     sy = sinh(qkl) / qk
  else
     cx = cosh(qkl)
     sx = sinh(qkl) / qk
     cy = cos(qkl)
     sy = sin(qkl) / qk
  endif

  !---- First-order terms.
  re(1,1) = cx
  re(1,2) = sx
  re(2,1) = - sk1 * sx
  re(2,2) = cx
  re(3,3) = cy
  re(3,4) = sy
  re(4,3) = + sk1 * sy
  re(4,4) = cy
  re(5,6) = el/(beta*gamma)**2

  ek(5) = el*dtbyds

  !---- Second-order terms.
  if (fsec) then
     biby4 = one / (four * beta)

     te(1,1,6) = + sk1 * el * sx * biby4
     te(1,6,1) = te(1,1,6)
     te(2,2,6) = te(1,1,6)
     te(2,6,2) = te(1,1,6)
     te(1,2,6) = - (sx + el*cx) * biby4
     te(1,6,2) = te(1,2,6)
     te(2,1,6) = - sk1 * (sx - el*cx) * biby4
     te(2,6,1) = te(2,1,6)

     te(3,3,6) = - sk1 * el * sy * biby4
     te(3,6,3) = te(3,3,6)
     te(4,4,6) = te(3,3,6)
     te(4,6,4) = te(3,3,6)
     te(3,4,6) = - (sy + el*cy) * biby4
     te(3,6,4) = te(3,4,6)
     te(4,3,6) = + sk1 * (sy - el*cy) * biby4
     te(4,6,3) = te(4,3,6)

     te(5,1,1) = - sk1 * (el - sx*cx) * biby4
     te(5,1,2) = + sk1 * sx**2 * biby4
     te(5,2,1) = te(5,1,2)
     te(5,2,2) = - (el + sx*cx) * biby4
     te(5,3,3) = + sk1 * (el - sy*cy) * biby4
     te(5,3,4) = - sk1 * sy**2 * biby4
     te(5,4,3) = te(5,3,4)
     te(5,4,4) = - (el + sy*cy) * biby4
     te(5,6,6) = (- six * re(5,6)) * biby4
  endif

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)
  !---- Apply tilt.
  if (tilt .ne. zero) call tmtilt(fsec,tilt,ek,re,te)

end SUBROUTINE qdbody

SUBROUTINE tmsep(fsec,ftrk,fcentre,orbit,fmap,dl,ek,re,te)
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : deltap, pc, charge
  use math_constfi, only : zero, one, two, ten3m
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for electrostatic separator.                       *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     dl        (double)  slice length.                                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is rotated (in tmmap) before and after this    *
  !     routine                                                          *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap, fcentre
  double precision :: dl
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  logical :: cplxy
  integer :: elpar_vl
  double precision :: tilt, ekick, efield, exfld, eyfld
  double precision :: ct, st, tmp

  integer, external :: el_par_vector
  double precision, external :: node_value

  !---- Initialize.
  st = zero
  ct = zero
  cplxy = .false.

  fmap = dl .ne. zero
  if (.not. fmap) return

  if (ftrk) then
     !-- get element parameters
     elpar_vl = el_par_vector(e_ey, g_elpar)
     !---- Strength and tilt.
     exfld = g_elpar(e_ey) !--This is a correct. Needs to be like this because of how the tilt is defined.
     eyfld = g_elpar(e_ex) !--This is a correct. Needs to be like this because of how the tilt is defined.
     tilt = g_elpar(e_tilt)
     if (eyfld.ne.zero) then
        tilt = -atan2(eyfld, exfld) + tilt
     endif

     if (tilt .ne. zero)  then
        !---  rotate orbit before entry
        st = sin(tilt)
        ct = cos(tilt)
        tmp = orbit(1)
        orbit(1) = ct * tmp + st * orbit(3)
        orbit(3) = ct * orbit(3) - st * tmp
        tmp = orbit(2)
        orbit(2) = ct * tmp + st * orbit(4)
        orbit(4) = ct * orbit(4) - st * tmp
     endif

     efield = sqrt(exfld**2 + eyfld**2)

  else
     exfld = zero
     eyfld = zero
     efield = zero
     tilt = zero
  endif

  cplxy = cplxy .or. tilt .ne. zero

  ekick  = efield * ten3m * charge / (pc * (one + deltap))

  call spbody(fsec,ftrk,tilt,ekick,orbit,dl,ek,re,te)
  if (fcentre) return

  if (tilt .ne. zero)  then
     !---  rotate orbit at exit
     tmp = orbit(1)
     orbit(1) = ct * tmp - st * orbit(3)
     orbit(3) = ct * orbit(3) + st * tmp
     tmp = orbit(2)
     orbit(2) = ct * tmp - st * orbit(4)
     orbit(4) = ct * orbit(4) + st * tmp
  endif

end SUBROUTINE tmsep

SUBROUTINE spbody(fsec,ftrk,tilt,ekick,orbit,el,ek,re,te)
  use twissbeamfi, only : beta, gamma, dtbyds
  use math_constfi, only : zero, one, two, three
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for electrostatic separator.                       *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     tilt      (double)  map tilt angle                               *
  !     ekick     (double)  processed electrostatic kick                 *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk
  double precision :: el, tilt, ekick
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: ekl, ch, sh, sy, dy, fact

  double precision, parameter :: by2=1d0/2d0, by6=1d0/6d0, by24=1d0/24d0, eps=1d-4

  !---- Prepare linear transformation parameters.
  !     DY = (COSH(K*L) - 1) / K.
  ekl = ekick * el
  if (abs(ekl) .gt. eps) then
     ch = cosh(ekl)
     sh = sinh(ekl)
     sy = sh / ekick
     dy = (ch - one) / ekick**2
  else
     ch = (one + by2  * ekl**2)
     sy = (one + by6  * ekl**2) * el
     sh = sy * ekick
     dy = (by2 + by24 * ekl**2) * el**2
  endif

  !---- Kicks.
  ek(3) = dy * (ekick / beta)
  ek(4) = sy * (ekick / beta)
  ek(5) = el * dtbyds

  !---- First-order terms.
  re(1,2) = el
  re(3,3) = ch - ekl * sh / beta**2
  re(3,4) = sy
  re(3,6) = (dy - el * sy / beta**2) * ekick
  re(4,3) = (sh - ekl * ch / beta**2) * ekick
  re(4,4) = ch
  re(4,6) = (sh - ekl * ch / beta**2)
  re(5,3) = - re(4,6)
  re(5,4) = - dy * ekick
  re(5,6) = - (sy - el * ch / beta**2)

  !---- Second-order terms.
  if (fsec) then
     fact = el / (two * beta)
     te(1,2,3) = - fact * ekick
     te(1,2,6) = - fact
     fact = el * (three*sh/gamma**2 + ekl*ch) / (two*beta**3)
     te(3,3,3) = fact * ekick**2
     te(3,3,6) = fact * ekick
     te(3,6,6) = fact
     fact = el * (three*ch/gamma**2 + ekl*sh) / (two*beta**3)
     te(4,3,3) = fact * ekick**3
     te(4,3,6) = fact * ekick**2
     te(4,6,6) = fact * ekick
     te(5,3,3) = - fact * ekick**2
     te(5,3,6) = - fact * ekick
     te(5,6,6) = - fact
     fact = el * sh / (two * beta)
     te(3,2,2) = fact
     te(3,4,4) = fact
     te(4,3,4) = - fact * ekick**2
     te(4,4,6) = - fact * ekick
     te(5,3,4) = fact * ekick
     te(5,4,6) = fact
     fact = el * ch / (two * beta)
     te(3,3,4) = - fact * ekick
     te(3,4,6) = - fact
     te(4,2,2) = fact * ekick
     te(4,4,4) = fact * ekick
     te(5,2,2) = - fact
     te(5,4,4) = - fact
     call tmsymm(te)
  endif

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)
  !---- Apply tilt.
  if (tilt .ne. zero) call tmtilt(fsec,tilt,ek,re,te)

end SUBROUTINE spbody

SUBROUTINE tmsext(fsec,ftrk,fcentre,orbit,fmap,el,dl,ek,re,te)
  use twtrrfi
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : radiate, deltap, gamma, arad
  use math_constfi, only : zero, one, two, three, twelve
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for sextupole element.                             *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     el        (double)  element length.                              *
  !     dl        (double)  slice length.                                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is rotated (in tmmap) before and after this    *
  !     routine                                                          *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap, fcentre
  double precision :: el, dl
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  logical :: cplxy
  integer :: i, j, n_ferr, elpar_vl
  double precision :: ct, st, tmp
  double precision :: f_errors(0:maxferr)
  double precision :: tilt, sk2, pt, sk2s, bvk, rfac

  integer, external :: el_par_vector, node_fd_errors
  double precision, external :: node_value, get_value
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')


  !---- Initialize.
  st = zero
  ct = zero
  cplxy = .false.
  fmap = el .ne. zero
  if (.not. fmap) return

  !---- Field error.
  F_ERRORS = zero
  n_ferr = node_fd_errors(f_errors)

  !-- get element parameters
  elpar_vl = el_par_vector(s_k2st, g_elpar)
  bvk = node_value('other_bv ')
  sk2  = bvk * ( g_elpar(s_k2)  + g_elpar(s_k2t)  +  f_errors(4)/el )
  sk2s = bvk * ( g_elpar(s_k2s) + g_elpar(s_k2st) +  f_errors(5)/el )
  tilt = node_value('tilt ')
  if (sk2s .ne. zero) then
     tilt = -atan2(sk2s, sk2)/three + tilt
     sk2 = sqrt(sk2**2 + sk2s**2)
  endif

  if (tilt .ne. zero)  then
     !---  rotate orbit before entry
     st = sin(tilt)
     ct = cos(tilt)
     tmp = orbit(1)
     orbit(1) = ct * tmp + st * orbit(3)
     orbit(3) = ct * orbit(3) - st * tmp
     tmp = orbit(2)
     orbit(2) = ct * tmp + st * orbit(4)
     orbit(4) = ct * orbit(4) - st * tmp
  endif
  cplxy = cplxy .or. (tilt .ne. zero)

  sk2 = sk2 / (one + deltap)

  !---- Half radiation effects at entrance.
  if (ftrk .and. radiate) then
     rfac = arad * gamma**3 * sk2**2 * el * (orbit(1)**2 + orbit(3)**2)**2 / twelve
     pt = orbit(6)
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

  call sxbody(fsec,ftrk,tilt,sk2,orbit,dl,ek,re,te)
  if (fcentre) return

  !---- Half radiation effects at exit.
  if (ftrk) then
     if (radiate) then
        rfac = arad * gamma**3 * sk2**2 * el * (orbit(1)**2 + orbit(3)**2)**2 / twelve
        pt = orbit(6)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        orbit(2) = orbit(2) * f_damp_t;
        orbit(4) = orbit(4) * f_damp_t;
        orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
     endif
  endif

  if (tilt .ne. zero)  then
     !---  rotate orbit at exit
     tmp = orbit(1)
     orbit(1) = ct * tmp - st * orbit(3)
     orbit(3) = ct * orbit(3) + st * tmp
     tmp = orbit(2)
     orbit(2) = ct * tmp - st * orbit(4)
     orbit(4) = ct * orbit(4) + st * tmp
  endif

end SUBROUTINE tmsext

SUBROUTINE sxbody(fsec,ftrk,tilt,sk2,orbit,el,ek,re,te)
  use twissbeamfi, only : beta, gamma, dtbyds
  use math_constfi, only : zero, two, three, four
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for sextupole element.                             *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     tilt      (double)  map tilt angle                               *
  !     sk2       (double)  processed sextupole kick                     *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk
  double precision :: tilt, sk2, el
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: skl, s1, s2, s3, s4

  !---- First-order terms.
  re(1,2) = el
  re(3,4) = el
  re(5,6) = el/(beta*gamma)**2
  ek(5) = el*dtbyds

  !---- Second-order terms.
  if (fsec) then
     skl = sk2 * el
     if (skl .ne. zero) then
        s1 = skl / two
        s2 = s1 * el / two
        s3 = s2 * el / three
        s4 = s3 * el / four
        te(1,1,1) = - s2
        te(1,1,2) = - s3
        te(1,2,2) = - two * s4
        te(1,3,3) = + s2
        te(1,3,4) = + s3
        te(1,4,4) = + two * s4
        te(2,1,1) = - s1
        te(2,1,2) = - s2
        te(2,2,2) = - two * s3
        te(2,3,3) = + s1
        te(2,3,4) = + s2
        te(2,4,4) = + two * s3
        te(3,1,3) = + s2
        te(3,1,4) = + s3
        te(3,2,3) = + s3
        te(3,2,4) = + two * s4
        te(4,1,3) = + s1
        te(4,1,4) = + s2
        te(4,2,3) = + s2
        te(4,2,4) = + two * s3
     endif
     te(1,2,6) = - el / (two * beta)
     te(3,4,6) = te(1,2,6)
     te(5,2,2) = te(1,2,6)
     te(5,4,4) = te(1,2,6)
     te(5,6,6) = - three * re(5,6) / (two * beta)
     call tmsymm(te)
  endif

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)
  !---- Apply tilt.
  if (tilt .ne. zero) call tmtilt(fsec,tilt,ek,re,te)

end SUBROUTINE sxbody

SUBROUTINE tmsol(fsec,ftrk,orbit,fmap,dl,ek,re,te)
  use twisslfi
  use math_constfi, only : zero, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for solenoid element, with centre option.          *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     dl        (double)  slice length.                                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  double precision :: dl
  double precision :: orbit(6), ek(6),re(6,6),te(6,6,6)

  if (dl .eq. zero) then
     call tmsol_th(ftrk,orbit,fmap,ek,re,te)
     return
  endif

  call tmsol0(fsec,ftrk,orbit,fmap,dl,ek,re,te)

end SUBROUTINE tmsol

SUBROUTINE tmsol0(fsec,ftrk,orbit,fmap,el,ek,re,te)
  use twissbeamfi, only : radiate, deltap, beta, gamma, dtbyds, arad
  use math_constfi, only : zero, one, two, three, six
    use matrices, only : EYE
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for solenoid element, no centre option.            *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  double precision :: el
  double precision :: orbit(6), ek(6), re(6,6), ek_t1(6), ek_t2(6), re_t1(6,6)
  double precision :: re_t2(6,6), te(6,6,6), te_t1(6,6,6), ek2(6)
  double precision :: ek_s(6), re_s(6,6), te_s(6,6,6)
  logical :: cplxy
  double precision :: sks, sk, skl, bvk, pxbeta, beta0, startrot
  double precision :: co, si, sibk, temp, xtilt,xtilt_rad, dl

  double precision, external :: node_value, get_value
  double precision, parameter :: ten5m=1d-5
  double precision :: rfac, kx, ky
  double precision :: pt, bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  beta0   = get_value('probe ','beta ')

  !---- Initialize.
  fmap = el .ne. zero
  if (.not. fmap) return
  EK = zero
  RE = EYE
  ek_s = zero
  re_s = EYE
  !---- Strength.
  sks = node_value('ks ')
  xtilt_rad = node_value('xtilt ')
  startrot =node_value('rot_start ')

  re_t1 = EYE
  re_t2 = EYE
  ek_t1 = zero
  ek_t2 = zero
  te_s = zero

  if (sks .ne. zero) cplxy = .true.

  !---- BV flag
  bvk = node_value('other_bv ')
  sks = sks * bvk

  !---- Set up C's and S's.
  sk = sks / two / (one + deltap)
  skl = sk * el
  co = cos(skl)
  si = sin(skl)
  if (abs(skl) .lt. ten5m) then
     sibk = (one - skl**2/six) * el
  else
     sibk = si/sk
  endif

  !---- Half radiation effect at entry.
  if (radiate .and. ftrk) then
     kx = ((sk**2)*orbit(1)-sk*orbit(4))*el;
     ky = ((sk**2)*orbit(3)+sk*orbit(2))*el;
     rfac = (arad * gamma**3 / three) * (kx**2 + ky**2) / el;
     pt = orbit(6);
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

  !---- First-order terms.
  re_s(1,1) = co**2
  re_s(2,2) = re_s(1,1)
  re_s(3,3) = re_s(1,1)
  re_s(4,4) = re_s(1,1)

  re_s(1,2) = co * sibk
  re_s(3,4) = re_s(1,2)
  re_s(1,3) = co * si
  re_s(2,4) = re_s(1,3)
  re_s(3,1) = - re_s(1,3)
  re_s(4,2) = re_s(3,1)
  re_s(2,1) = sk * re_s(3,1)
  re_s(4,3) = re_s(2,1)
  re_s(1,4) = si * sibk
  re_s(3,2) = - re_s(1,4)
  re_s(4,1) = sk * si**2
  re_s(2,3) = - re_s(4,1)
  re_s(5,6) = el/(beta*gamma)**2

  ek_s(5) = el*dtbyds

  !---- Second-order terms.
  if (fsec) then
     temp = el * co * si / beta
     te_s(1,4,6) = - temp
     te_s(3,2,6) =   temp
     te_s(1,1,6) =   temp * sk
     te_s(2,2,6) =   temp * sk
     te_s(3,3,6) =   temp * sk
     te_s(4,4,6) =   temp * sk
     te_s(2,3,6) =   temp * sk**2
     te_s(4,1,6) = - temp * sk**2

     temp = el * (co**2 - si**2) / (two * beta)
     te_s(1,2,6) = - temp
     te_s(3,4,6) = - temp
     te_s(1,3,6) = - temp * sk
     te_s(2,4,6) = - temp * sk
     te_s(3,1,6) =   temp * sk
     te_s(4,2,6) =   temp * sk
     te_s(2,1,6) =   temp * sk**2
     te_s(4,3,6) =   temp * sk**2

     temp = el / (two * beta)
     te_s(5,2,2) = - temp
     te_s(5,4,4) = - temp
     te_s(5,1,4) =   temp * sk
     te_s(5,2,3) = - temp * sk
     te_s(5,1,1) = - temp * sk**2
     te_s(5,3,3) = - temp * sk**2
     te_s(5,6,6) = - three * re_s(5,6) / (two * beta)
     call tmsymm(te_s)
  endif


  !---- Track orbit.

  if (ftrk) then

    if(abs(xtilt_rad) > ten5m) then
      te_t1 = zero
      xtilt = -sin(xtilt_rad)

      pxbeta = xtilt*startrot/beta
      ek_t1(1) =  startrot*xtilt
      ek_t1(2) =  xtilt
      ek_t1(5) = -0.5d0*pxbeta*xtilt
      re_t1(1,6) = -pxbeta
      re_t1(5,2) = -pxbeta
      call tmtrak(ek_t1,re_t1,te_t1,orbit,orbit)
      call tmcat(.true.,re_t1,te_t1,re,te,re,te)


      call tmtrak(ek_s,re_s,te_s,orbit,orbit) ! Calls the normal solenoid
      call tmcat(.true.,re_s,te_s,re,te,re,te)

      !To tilt it back
      xtilt=-xtilt
      pxbeta = xtilt*(el+startrot)/beta
      ek_t2(1) =  (el+startrot)*xtilt
      ek_t2(2) =  xtilt
      ek_t2(5) = -0.5d0*pxbeta*xtilt
      re_t2(1,6) = -pxbeta
      re_t2(5,2) = -pxbeta

      call tmtrak(ek_t2,re_t2,te_t1 ,orbit,orbit)
      call tmcat(.true.,re_t2,te_t1,re,te,re,te)

    else
      ek=ek_s
      re=re_s
      te=te_s
      call tmtrak(ek,re,te,orbit,orbit)
    endif
  else
    ek=ek_s
    re=re_s
    te=te_s
  endif


  !---- Half radiation effect at exit.
  if (radiate .and. ftrk) then
     kx = ((sk**2)*orbit(1)-sk*orbit(4))*el;
     ky = ((sk**2)*orbit(3)+sk*orbit(2))*el;
     rfac = (arad * gamma**3 / three) * (kx**2 + ky**2) / el;
     pt = orbit(6);
     bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
     f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
     orbit(2) = orbit(2) * f_damp_t;
     orbit(4) = orbit(4) * f_damp_t;
     orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
  endif

end SUBROUTINE tmsol0

SUBROUTINE tmtrans(fsec,ftrk,orbit,fmap,ek,re,te)
  use twisslfi
  use twissbeamfi, only : beta
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for translation.                         *
  !     Treated in a purely linear way.                                  *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: ftrk, fmap,fsec
  double precision :: orbit(6);

  double precision :: x, y, z
  double precision :: node_value, ek(6), re(6,6), te(6,6,6)


 !---- Get translation parameters
 x    = node_value('dx ')
 y    = node_value('dy ')
 z    = node_value('ds ')

 ek(1) = ek(1) - x
 ek(3) = ek(3) - y
 ek(5) = ek(5) - z/beta

  !---- Track orbit.
 if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmtrans

SUBROUTINE tmsrot(ftrk,orbit,fmap,ek,re,te)
  use twisslfi
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for rotation about S-axis.                         *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: ftrk, fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  logical :: cplxy
  double precision :: psi, ct, st

  double precision, external :: node_value

  !---- Initialize.
  psi = node_value('angle ')
  fmap = psi .ne. zero
  if (.not. fmap) return

  !---- First-order terms.
  cplxy = .true.
  ct = cos(psi)
  st = sin(psi)

  re(1,1) = ct
  re(1,3) = st
  re(3,1) = -st
  re(3,3) = ct
  re(2,2) = ct
  re(2,4) = st
  re(4,2) = -st
  re(4,4) = ct


  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmsrot
SUBROUTINE tmxrot(ftrk,orbit,fmap,ek,re,te)
  use twisslfi
  use twissbeamfi, only : beta
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for rotation about X-axis.                         *
  !     Treated in a purely linear way.                                  *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: ftrk, fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: angle, ca, sa, ta
  double precision :: node_value

  !---- Initialize.
  angle = node_value('angle ')
  if (angle .eq. 0) return

  angle = angle * node_value('other_bv ')

  !---- Kick.
  ca = cos(angle)
  sa = sin(angle)
  ta = tan(angle)

  ek(4) = sa

  !---- Transfer matrix.
  re(3,3) = 1/ca
  re(4,4) =   ca
  re(4,6) =   sa/beta
  re(5,3) =  -ta/beta

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmxrot

SUBROUTINE tmyrot(ftrk,orbit,fmap,ek,re,te)
  use twisslfi
  use twissbeamfi, only : beta
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for rotation about Y-axis.                         *
  !     Treated in a purely linear way.                                  *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: ftrk, fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: angle, ca, sa, ta
  double precision :: node_value

  !---- Initialize.
  angle = node_value('angle ')
  if (angle .eq. 0) return

  angle = angle * node_value('other_bv ')

  !---- Kick.
  ca = cos(angle)
  sa = sin(angle)
  ta = tan(angle)

  ek(2) = sa

  !---- Transfer matrix.
  re(1,1) = 1/ca
  re(2,2) =   ca
  re(2,6) =   sa/beta
  re(5,1) =  -ta/beta

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmyrot

SUBROUTINE tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)
  use twissbeamfi, only : beta, gamma, dtbyds
  use matrices, only : EYE
  use math_constfi, only : zero, two, three
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for drift space, no centre option.                 *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     dl        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  double precision :: dl
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  !---- Initialize.
  EK = zero
  RE = EYE
  if (fsec) TE = zero
  fmap = dl .ne. zero

  !---- First-order terms.

  re(1,2) = dl
  re(3,4) = dl
  re(5,6) = dl/(beta*gamma)**2

  ek(5) = dl*dtbyds

  !---- Second-order terms.
  if (fsec) then
     te(1,2,6) = - dl / (two * beta)
     te(1,6,2) = te(1,2,6)
     te(3,4,6) = te(1,2,6)
     te(3,6,4) = te(3,4,6)
     te(5,2,2) = te(1,2,6)
     te(5,4,4) = te(5,2,2)
     te(5,6,6) = te(1,2,6) * three / (beta * gamma) ** 2
  endif
  
  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmdrf

SUBROUTINE tmchp0(ftrk,orbit,fmap,ek,re, te)
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : deltap, pc, gamma, energy, beta
  use matrices, only : EYE
  use math_constfi, only : zero, one
  use phys_constfi, only : clight
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for change of reference Energy                     *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap, fcentre, fringe
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)


  double precision :: energy1, pc1, gamma1, pt, mass, beta1i
  double precision, external :: node_value, get_value
  integer, external :: el_par_vector

  mass     = get_value('probe ','mass ')
  re = EYE
  ek = zero
  te = zero
  pt = orbit(6) 

  !---- Transfer map.
  fmap = .true.
  
  energy1 = pt*pc+energy
  pc1 = sqrt(energy1**2-mass**2)
  gamma1 = energy1/mass
  beta1i = energy1/pc1

  re(2,2) = pc/pc1
  re(4,4) = pc/pc1
  re(6,6) = pc/pc1

  ek(6) = beta1i*(gamma/gamma1-one)

  if(ftrk) call tmtrak(ek,re,te,orbit,orbit)


end SUBROUTINE tmchp0


SUBROUTINE tmrf(fsec,ftrk,fcentre,orbit,fmap,el,ds,ek,re,te)
  use twisslfi
  use twtapering
  use twiss_elpfi
  use twissbeamfi, only : deltap, pc, radiate
  use matrices, only : EYE
  use math_constfi, only : zero, one, two, half, ten6p, ten3m, pi, twopi
  use phys_constfi, only : clight
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for RF cavity.                                     *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     el        (double)  element length.                              *
  !     ds        (double)  slice length.                                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*

  logical :: fsec, ftrk, fmap, fcentre, istaper, fringe
  double precision :: el, ds
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)
  double precision :: ek_f(6), re_f(6,6), te_f(6,6,6)
  double precision :: ek_tmp(6), re_tmp(6,6), te_tmp(6,6,6)
  double precision :: ek_ch(6), re_ch(6,6), te_ch(6,6,6)

  integer :: elpar_vl, ncav
  double precision :: rfv, rff, rfl, dl, omega, vrf, phirf, bvk
  double precision :: ek0(6), rw(6,6), tw(6,6,6)
  double precision :: c0, c1, c2, tmpphase

  double precision, external :: node_value, get_value
  integer, external :: el_par_vector, get_ncavities

  !-- get element parameters
  elpar_vl = el_par_vector(r_lagt, g_elpar)

  !---- Fetch voltage.
  rfv = g_elpar(r_volt)

  !---- Cavity not excited, use drift map.
  if (rfv .eq. zero) then
    call tmdrf(fsec,ftrk,orbit,fmap,ds,ek,re,te)
    return
  endif

  !---- Cavity is excited, use full map.

  !---- Initialize.
  EK0 = zero
  RW = EYE
  TW = zero

  ek_f = zero
  re_f = EYE
  te_f = zero

  !---- BV flag
  rff = g_elpar(r_freq)
  rfl = g_elpar(r_lag) +  g_elpar(r_lagt)
  bvk = node_value('other_bv ')
  

  !-- LD: 20.6.2014 (bvk=-1: not -V -> V but lag -> pi-lag)
  if (bvk .eq. -one) then
    rfl = 0.5-rfl
  endif

  !---- Set up.
  omega = rff * ten6p * twopi / clight
  vrf   = rfv * ten3m / (pc * (one + deltap))
  phirf = rfl * twopi - omega * orbit(5)

  istaper = get_value('twiss ','tapering ').ne.zero
  
  if(istaper .and. ftrk) then
    ncav = get_ncavities()   
    phirf = asin((sin(pi*half)*vrf - endpt/ncav)/vrf)
    tmpphase = (phirf+omega*orbit(5))/twopi-g_elpar(r_lag)
    call store_node_value('lagtap ', tmpphase)
  endif
  
  c0 =   vrf * sin(phirf)
  c1 = - vrf * cos(phirf) * omega
  c2 = - vrf * sin(phirf) * omega**2 * half

  !---- Transfer map.
  fmap = .true.

  !---- Sandwich cavity between two drifts.
  if (el .ne. zero) then
    fringe = node_value('fringe ') .gt. zero
    ! TODO: generalize for ds!=0.5
    dl = el / two
    
    if (fringe) then
      call tmrffringe(fsec,ftrk,orbit, fmap, el, one, ek, re_f, te_f)
      call tmdrf(fsec,ftrk,orbit,fmap,dl,ek0,rw,tw)
      call tmcat(fsec,rw,tw,re_f,te_f,rw,tw)
    else
      call tmdrf(fsec,ftrk,orbit,fmap,dl,ek0,rw,tw)
    endif


  if (ftrk) then
    orbit(6) = orbit(6) + c0
    ek(6) = c0
    re(6,5) = c1
    if (fsec) te(6,5,5) = c2
  else
    ek(6) = c0 - c1 * orbit(5) + c2 * orbit(5)**2
    re(6,5) = c1 - two * c2 * orbit(5)
    if (fsec) te(6,5,5) = c2
  endif

    call tmcat(fsec,re,te,rw,tw,re,te)

    
     ! call tmchp0(ftrk,orbit,fmap,ek_ch,re_ch, te_ch)
      !call tmcat(fsec,re_ch,te_ch,re,te,re,te)
    

    if (fcentre) return

    call tmdrf(fsec,ftrk,orbit,fmap,dl,ek0,rw,tw)
    call tmcat(fsec,rw,tw,re,te,re,te)
    if (fringe) then
      call tmrffringe(fsec,ftrk,orbit, fmap, el, -one, ek, re_f, te_f) 
      call tmcat(fsec,re_f,te_f,re,te,re,te)
    endif

  else
    if (ftrk) then
      orbit(6) = orbit(6) + c0
      ek(6) = c0
      re(6,5) = c1
      if (fsec) te(6,5,5) = c2
    else
      ek(6) = c0 - c1 * orbit(5) + c2 * orbit(5)**2
      re(6,5) = c1 - two * c2 * orbit(5)
      if (fsec) te(6,5,5) = c2
   endif
  endif

end SUBROUTINE tmrf

SUBROUTINE tmrffringe(fsec,ftrk,orbit, fmap, el, jc, ek, re, te) 
  use twisslfi
  use twiss_elpfi
  use twissbeamfi, only : deltap, pc, beta
  use matrices, only : EYE
  use math_constfi, only : zero, one, two, half, ten6p, ten3m, pi, twopi
  use phys_constfi, only : clight
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for RF cavity fringe.                              *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     fcentre   (logical) legacy centre behaviour (no exit effects).   *
  !     el        (double)  element length.                              *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  double precision :: el,  V, dpxy, dptxy, s1, c1, tcorr
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  integer :: elpar_vl
  double precision :: rff, rfl, dl, omega, vrf, phirf, bvk, jc
  double precision :: ek0(6), rw(6,6), tw(6,6,6)
  double precision :: rfv, x, px, y, py, t, pt

  double precision, external :: node_value
  integer, external :: el_par_vector
  ek = zero
  te = zero
  re = EYE

  !-- get element parameters
  elpar_vl = el_par_vector(r_freq, g_elpar)

  !---- Fetch voltage.
  rfv = g_elpar(r_volt)
  rff = g_elpar(r_freq)
  rfl = g_elpar(r_lag)

  omega = rff * (ten6p * twopi / clight)
  ! vrf   = rfv * ten3m / (pc * (one + deltas))
  vrf   = rfv * ten3m
  phirf = rfl * twopi
  t  = orbit(5);

  tcorr = jc*el/(2*beta)
  V = jc*vrf/(pc*el* (one + deltap))
  s1 = sin(phirf - omega*(t+tcorr))
  c1 = cos(phirf - omega*(t+tcorr))

  dpxy = -V*s1*half

  re(2,1)   = dpxy
  re(4,3)   = dpxy
  te(6,1,1) = 0.25d0*V*c1*omega
  te(6,3,3) = 0.25d0*V*c1*omega

  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmrffringe

SUBROUTINE tmcat(fsec,rb,tb,ra,ta,rd,td)
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Concatenate two TRANSPORT maps.                                  *
  !     This routine is time-critical and is carefully optimized.        *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     rb(6,6), tb(6,6,6)  second map in beam line order.               *
  !     ra(6,6), ta(6,6,6)  first map in beam line order.                *
  !     Output:                                                          *
  !     rd(6,6), td(6,6,6)  result map.                                  *
  !----------------------------------------------------------------------*
  logical :: fsec
  double precision, intent(IN)  :: ra(6,6), ta(6,6,6)
  double precision, intent(IN)  :: rb(6,6), tb(36,6)
  double precision, intent(OUT) :: rd(6,6), td(6,6,6)

  double precision :: rw(6,6), tw(6,6,6), ts(36,6)
  integer :: i1, i2, i3

  !---- Initialization
  TW = zero

  !---- Transfer matrix.
  do i2 = 1, 6
     do i1 = 1, 6
        rw(i1,i2) = rb(i1,1) * ra(1,i2) + rb(i1,2) * ra(2,i2)         &
             + rb(i1,3) * ra(3,i2) + rb(i1,4) * ra(4,i2)              &
             + rb(i1,5) * ra(5,i2) + rb(i1,6) * ra(6,i2)
     enddo
  enddo

  !---- Second order terms.
  if (fsec) then
     do i3 = 1, 6
        do i1 = 1, 36
           ts(i1,i3) = tb(i1,1) * ra(1,i3) + tb(i1,2) * ra(2,i3)       &
                + tb(i1,3) * ra(3,i3) + tb(i1,4) * ra(4,i3)            &
                + tb(i1,5) * ra(5,i3) + tb(i1,6) * ra(6,i3)
        enddo
     enddo
     do i2 = 1, 6
        do i3 = i2, 6
           do i1 = 1, 6
              tw(i1,i2,i3) =                                            &
                   rb(i1,1) * ta(1,i2,i3) + rb(i1,2) * ta(2,i2,i3)      &
                   + rb(i1,3) * ta(3,i2,i3) + rb(i1,4) * ta(4,i2,i3)    &
                   + rb(i1,5) * ta(5,i2,i3) + rb(i1,6) * ta(6,i2,i3)    &
                   + ts(i1,   i2) * ra(1,i3) + ts(i1+ 6,i2) * ra(2,i3)  &
                   + ts(i1+12,i2) * ra(3,i3) + ts(i1+18,i2) * ra(4,i3)  &
                   + ts(i1+24,i2) * ra(5,i3) + ts(i1+30,i2) * ra(6,i3)
              tw(i1,i3,i2) = tw(i1,i2,i3)
           enddo
        enddo
     enddo
  endif

  !---- Copy result.
  RD = RW
  TD = TW

end SUBROUTINE tmcat

SUBROUTINE tmcat1(fsec,eb,rb,tb,ea,ra,ta,ed,rd,td)
  use math_constfi, only : two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Concatenate two TRANSPORT maps including zero-order terms.       *
  !     This routine is time-critical and is carefully optimized.        *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     eb(6), rb(6,6), tb(6,6,6)  second map in beam line order.        *
  !     ea(6), ra(6,6), ta(6,6,6)  first map in beam line order.         *
  !     Output:                                                          *
  !     ed(6), rd(6,6), td(6,6,6)  result map.                           *
  !----------------------------------------------------------------------*
  logical :: fsec
  double precision, intent(IN)  :: ea(6), ra(6,6), ta(6,6,6)
  double precision, intent(IN)  :: eb(6), rb(6,6), tb(36,6)
  double precision, intent(OUT) :: ed(6), rd(6,6), td(6,6,6)

  integer :: i, ij, j, k
  double precision :: ew(6), rw(6,6), tw(6,6,6), es(6,6), ts(36,6)

  if ( .not. fsec) then
     !---- First order only
     do k = 1, 6
        !---- Zero-order terms.
        ew(k) = eb(k) + rb(k,1) * ea(1) + rb(k,2) * ea(2) &
                      + rb(k,3) * ea(3) + rb(k,4) * ea(4) &
                      + rb(k,5) * ea(5) + rb(k,6) * ea(6)

        !---- First-order terms.
        do j = 1, 6
           rw(j,k) = rb(j,1) * ra(1,k) + rb(j,2) * ra(2,k) &
                   + rb(j,3) * ra(3,k) + rb(j,4) * ra(4,k) &
                   + rb(j,5) * ra(5,k) + rb(j,6) * ra(6,k)
        enddo
     enddo

  else
     !---- Second order terms.

     !---- Auxiliary terms.
     do k = 1, 6

        !---- Sum over S of TB(I,S,K) * EA(S).
        do i = 1, 6
           es(i,k) = tb(i   ,k) * ea(1) + tb(i+ 6,k) * ea(2) &
                   + tb(i+12,k) * ea(3) + tb(i+18,k) * ea(4) &
                   + tb(i+24,k) * ea(5) + tb(i+30,k) * ea(6)
        enddo

        !---- Sum over S of TB(I,J,S) * RA(S,K).
        do ij = 1, 36
           ts(ij,k) = tb(ij,1) * ra(1,k) + tb(ij,2) * ra(2,k) &
                    + tb(ij,3) * ra(3,k) + tb(ij,4) * ra(4,k) &
                    + tb(ij,5) * ra(5,k) + tb(ij,6) * ra(6,k)
        enddo
     enddo

     !---- Final values.
     do k = 1, 6

        !---- Zero-order terms.
        ew(k) = eb(k) + (rb(k,1) + es(k,1)) * ea(1) &
                      + (rb(k,2) + es(k,2)) * ea(2) &
                      + (rb(k,3) + es(k,3)) * ea(3) &
                      + (rb(k,4) + es(k,4)) * ea(4) &
                      + (rb(k,5) + es(k,5)) * ea(5) &
                      + (rb(k,6) + es(k,6)) * ea(6)

        !---- First-order terms.
        do j = 1, 6
           rw(j,k) = (rb(j,1) + two * es(j,1)) * ra(1,k) &
                   + (rb(j,2) + two * es(j,2)) * ra(2,k) &
                   + (rb(j,3) + two * es(j,3)) * ra(3,k) &
                   + (rb(j,4) + two * es(j,4)) * ra(4,k) &
                   + (rb(j,5) + two * es(j,5)) * ra(5,k) &
                   + (rb(j,6) + two * es(j,6)) * ra(6,k)
        enddo

        !---- Second-order terms.
        do j = k, 6
           do i = 1, 6
              tw(i,j,k) = + (rb(i,1)+two*es(i,1))*ta(1,j,k) + ts(i   ,j)*ra(1,k) &
                          + (rb(i,2)+two*es(i,2))*ta(2,j,k) + ts(i+ 6,j)*ra(2,k) &
                          + (rb(i,3)+two*es(i,3))*ta(3,j,k) + ts(i+12,j)*ra(3,k) &
                          + (rb(i,4)+two*es(i,4))*ta(4,j,k) + ts(i+18,j)*ra(4,k) &
                          + (rb(i,5)+two*es(i,5))*ta(5,j,k) + ts(i+24,j)*ra(5,k) &
                          + (rb(i,6)+two*es(i,6))*ta(6,j,k) + ts(i+30,j)*ra(6,k)
              tw(i,k,j) = tw(i,j,k)
           enddo
        enddo
     enddo

     !---- Copy second-order terms.
     TD = TW

  endif

  !---- Copy zero- and first-order terms.
  ED = EW
  RD = RW

end SUBROUTINE tmcat1

SUBROUTINE tmtrak(ek,re,te,orb1,orb2)
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track orbit and change reference for RE matrix.                  *
  !     Input:                                                           *
  !     ek(6)     (double)  kick on orbit.                               *
  !     re(6,6)   (double)  transfer matrix before update.               *
  !     te(6,6,6) (double)  second order terms.                          *
  !     orb1(6)   (double)  orbit before element.                        *
  !     Output:                                                          *
  !     orb2(6)   (double)  orbit after element.                         *
  !     re(6,6)   (double)  transfer matrix after update.                *
  !----------------------------------------------------------------------*
  double precision, intent(IN)     :: ek(6), te(6,6,6), orb1(6)
  double precision, intent(IN OUT) :: re(6,6)
  double precision, intent(OUT)    :: orb2(6)

  integer :: i, k, l
  double precision :: sum1, sum2, temp(6)

  integer, external :: get_option

  do i = 1, 6
     sum2 = ek(i)
     do k = 1, 6
        sum1 = zero
        do l = 1, 6
           sum1 = sum1 + te(i,k,l) * orb1(l)
        enddo
        sum2 = sum2 + (re(i,k) + sum1) * orb1(k)
        re(i,k) = re(i,k) + sum1 + sum1
     enddo
     temp(i) = sum2
  enddo

  ORB2(1:6) = TEMP(1:6)

  !---- Symplectify transfer matrix.
  if (get_option('sympl ') .ne. 0) call tmsymp(re)

end SUBROUTINE tmtrak

SUBROUTINE tmsymp(r)
  use matrices
  use math_constfi, only : zero, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Symplectify a 6 by 6 matrix R.                                   *
  !     Algorithm described in the doctoral thesis by Liam Healey.       *
  !     Input:                                                           *
  !     r(6,6)    (double)  matrix to be symplectified.                  *
  !     Output:                                                          *
  !     r(6,6)    (double)  the symplectified matrix.                    *
  !----------------------------------------------------------------------*
  double precision, intent(IN OUT) :: r(6,6)

  logical :: eflag, error
  integer :: i, j
  double precision :: a(6,6), b(6,6), v(6,6), nrm

  A = -R + EYE
  B =  R + EYE

  call m66div(A,B,V,eflag) ! V = A * B-1
  if (eflag) goto 100

  ! A = Jt*Vt*J ; ie A = V**(-1)
  A = matmul(JMATT, matmul(transpose(V),JMAT)) ! invert symplectic matrix

  A = ( A - V ) / two
  B = -A
  A = A + EYE
  B = B + EYE

  call m66div(A,B,V,eflag) ! V = A * B-1
  if (eflag) goto 100

  R = V
  return

100 continue
  call m66symp(r,nrm)
  if (nrm .gt. zero) then
     print *," Singular matrix occurred during symplectification of R (left unchanged)."
     print *," The column norm of R'*J*R-J is ",nrm
  endif

end SUBROUTINE tmsymp

SUBROUTINE tmsymm(t)
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Symmetrize second-order array T.                                 *
  !     Input:                                                           *
  !     t(6,6,6)  (double)  array to be symmetrized.                     *
  !     Output:                                                          *
  !     t(6,6,6)  (double)  symmetrized array.                           *
  !----------------------------------------------------------------------*
  double precision :: t(6,6,6)

  integer :: i, j, k

  do k = 1, 5
     do j = k+1, 6
        do i = 1, 6
           t(i,j,k) = t(i,k,j)
        enddo
     enddo
  enddo

end SUBROUTINE tmsymm

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
  double precision :: anum(6,6), aden(6,6), target(6,6)
  logical :: eflag

  integer :: i, j, irank
  double precision :: augmat(6,12)

  !---- Copy input to local array.
  AUGMAT(:,1:6) = ADEN
  AUGMAT(:,7:12) = ANUM

  !---- Solve resulting system.
  call solver(AUGMAT,6,6,irank)
  if (irank .lt. 6) then
     eflag = .true.
  else
     !---- Copy result.
     eflag = .false.
     TARGET = AUGMAT(:,7:12)
  endif

end subroutine m66div

subroutine m66symp(r,nrm) ! can be moved to twiss.f90
  use matrices, only : JMAT
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Check if a 6 by 6 matrix R is symplectic.
  ! Input:
  !   r(6,6)    (double)  Matrix R to check
  ! Output:
  !   nrm       (double)  The column norm of (R'*J*R)-J
  !----------------------------------------------------------------------*
  double precision :: R(6,6), nrm

  double precision :: T(6,6)
  double precision :: sum
  double precision, parameter :: zero=0d0
  integer :: i, j

  ! calculate norm( (Rt J R) - J)

  T = matmul( matmul(transpose(R), JMAT), R) - JMAT

  nrm = zero
  do j = 1, 6
     sum = zero
     do i = 1, 6
        sum = sum + abs(T(i,j))
     enddo
     nrm = max(nrm,sum)
  enddo

end subroutine m66symp

SUBROUTINE tmali1(orb1, errors, beta, gamma, orb2, rm)
  use twiss0fi, only : align_max
  use matrices, only : EYE
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for orbit displacement at entry of an element.     *
  !       called from EMIT, TRRUN and TWISS modules.                     *
  !     Input:                                                           *
  !     orb1(6)   (real)    Orbit before misalignment.                   *
  !     errors(align_max) (real)    alignment errors                     *
  !     beta (real)         beam beta relativistic factor                *
  !     gamma (real)        beam gamma relativistic factor               *
  !     Output:                                                          *
  !     orb2(6)   (real)    Orbit after misalignment.                    *
  !     rm(6,6)   (real)    First order transfer matrix w.r.t. orbit.    *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: orb1(6), errors(align_max)
  double precision, intent(IN)  :: beta, gamma
  double precision, intent(OUT) :: orb2(6),  rm(6,6)

  double precision :: orbt(6), w(3,3)
  double precision :: ds, dx, dy, the, phi, psi, s2

  !---- Build rotation matrix and compute additional drift length.
  dx  = errors(1)
  dy  = errors(2)
  ds  = errors(3)
  the = errors(5)
  phi = errors(4)
  psi = errors(6)
  call sumtrx(the, phi, psi, w)
  s2 = (w(1,3) * dx + w(2,3) * dy + w(3,3) * ds) / w(3,3)

  !---- F2 terms (transfer matrix).
  RM = EYE
  rm(2,2) = w(1,1)
  rm(2,4) = w(2,1)
  rm(2,6) = w(3,1) / beta
  rm(4,2) = w(1,2)
  rm(4,4) = w(2,2)
  rm(4,6) = w(3,2) / beta

  rm(1,1) =   w(2,2) / w(3,3)
  rm(1,2) = rm(1,1) * s2
  rm(1,3) = - w(1,2) / w(3,3)
  rm(1,4) = rm(1,3) * s2
  rm(3,1) = - w(2,1) / w(3,3)
  rm(3,2) = rm(3,1) * s2
  rm(3,3) =   w(1,1) / w(3,3)
  rm(3,4) = rm(3,3) * s2
  rm(5,1) = w(1,3) / (w(3,3) * beta)
  rm(5,2) = rm(5,1) * s2
  rm(5,3) = w(2,3) / (w(3,3) * beta)
  rm(5,4) = rm(5,3) * s2
  rm(5,6) = - s2 / (beta * gamma)**2

  !---- Track orbit.
  ORBT = matmul(RM,ORB1)
  orb2(1) = orbt(1) - (w(2,2) * dx - w(1,2) * dy) / w(3,3)
  orb2(2) = orbt(2) + w(3,1)
  orb2(3) = orbt(3) - (w(1,1) * dy - w(2,1) * dx) / w(3,3)
  orb2(4) = orbt(4) + w(3,2)
  orb2(5) = orbt(5) - s2 / beta
  orb2(6) = orbt(6)


end SUBROUTINE tmali1

SUBROUTINE tmali2(el, orb1, errors, beta, gamma, orb2, rm)
  use twiss0fi, only : align_max
  use matrices, only : EYE
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for orbit displacement at exit of an element.      *
  !       called from EMIT, TRRUN and TWISS modules.                     *
  !     Input:                                                           *
  !     el (real)           Element length                               *
  !     orb1(6)   (real)    Orbit before misalignment.                   *
  !     errors(align_max) (real)    alignment errors                     *
  !     beta (real)         beam beta relativistic factor                *
  !     gamma (real)        beam gamma relativistic factor               *
  !     Output:                                                          *
  !     orb2(6)   (real)    Orbit after misalignment.                    *
  !     rm(6,6)   (real)    First order transfer matrix w.r.t. orbit.    *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: el, orb1(6), errors(align_max)
  double precision, intent(IN)  :: beta, gamma
  double precision, intent(OUT) :: orb2(6), rm(6,6)

  double precision :: orbt(6), v(3), ve(3), w(3,3), we(3,3)
  double precision :: ds, dx, dy, the, phi, psi, s2, tilt

  !---- Misalignment rotation matrix w.r.t. entrance system.
  dx  = errors(1)
  dy  = errors(2)
  ds  = errors(3)
  the = errors(5)
  phi = errors(4)
  psi = errors(6)
  tilt = zero
  call sumtrx(the, phi, psi, w)

  !---- VE and WE represent the change of reference.
  call suelem(el, ve, we, tilt)

  !---- Misalignment displacements at exit w.r.t. entrance system.
  v(1) = dx + w(1,1)*ve(1)+w(1,2)*ve(2)+w(1,3)*ve(3)-ve(1)
  v(2) = dy + w(2,1)*ve(1)+w(2,2)*ve(2)+w(2,3)*ve(3)-ve(2)
  v(3) = ds + w(3,1)*ve(1)+w(3,2)*ve(2)+w(3,3)*ve(3)-ve(3)

  !---- Convert all references to exit, build additional drift.
  !call sutran(w, v, we)
  V = matmul(transpose(WE),V)
  W = matmul(matmul(transpose(WE),W),WE)
  s2 = - (w(1,3) * v(1) + w(2,3) * v(2) + w(3,3) * v(3)) / w(3,3)

  !---- Transfer matrix.
  RM = EYE
  rm(1,1) = w(1,1)
  rm(3,1) = w(2,1)
  rm(5,1) = w(3,1) / beta
  rm(1,3) = w(1,2)
  rm(3,3) = w(2,2)
  rm(5,3) = w(3,2) / beta

  rm(2,2) =   w(2,2) / w(3,3)
  rm(1,2) = rm(2,2) * s2
  rm(4,2) = - w(1,2) / w(3,3)
  rm(3,2) = rm(4,2) * s2
  rm(2,4) = - w(2,1) / w(3,3)
  rm(1,4) = rm(2,4) * s2
  rm(4,4) =   w(1,1) / w(3,3)
  rm(3,4) = rm(4,4) * s2
  rm(2,6) = w(1,3) / (w(3,3) * beta)
  rm(1,6) = rm(2,6) * s2
  rm(4,6) = w(2,3) / (w(3,3) * beta)
  rm(3,6) = rm(4,6) * s2
  rm(5,6) = - s2 / (beta * gamma)**2

  !---- Track orbit.
  orbt(1) = orb1(1) + (w(2,2) * v(1) - w(1,2) * v(2)) / w(3,3)
  orbt(2) = orb1(2) - w(3,1)
  orbt(3) = orb1(3) + (w(1,1) * v(2) - w(2,1) * v(1)) / w(3,3)
  orbt(4) = orb1(4) - w(3,2)
  orbt(5) = orb1(5) - s2 / beta
  orbt(6) = orb1(6)
  ORB2 = matmul(RM,ORBT)

end SUBROUTINE tmali2

SUBROUTINE ccperrf(xx, yy, wx, wy)
  use math_constfi, only : zero, one, two, half
  implicit none
  !----------------------------------------------------------------------*
  !     purpose:                                                         *
  !     modification of wwerf, double precision complex error function,  *
  !     written at cern by K. Koelbig.                                   *
  !     input:                                                           *
  !     xx, yy    (double)    real + imag argument                       *
  !     output:                                                          *
  !     wx, wy    (double)    real + imag function result                *
  !----------------------------------------------------------------------*
  double precision :: xx, yy, wx, wy

  integer :: n, nc, nu
  double precision :: x, y, q, h, xl, xh, yh, tx, ty, tn, sx, sy, saux
  double precision :: rx(33), ry(33)

  double precision, parameter :: xlim=5.33d0, ylim=4.29d0, fac1=3.2d0, fac2=23d0, fac3=21d0
  double precision, parameter :: cc=1.12837916709551d0

  x = abs(xx)
  y = abs(yy)

  if (y .lt. ylim  .and.  x .lt. xlim) then
     q  = (one - y / ylim) * sqrt(one - (x/xlim)**2)
     h  = one / (fac1 * q)
     nc = 7 + int(fac2*q)
     xl = h**(1 - nc)
     xh = y + half/h
     yh = x
     nu = 10 + int(fac3*q)
     rx(nu+1) = zero
     ry(nu+1) = zero

     do n = nu, 1, -1
        tx = xh + n * rx(n+1)
        ty = yh - n * ry(n+1)
        tn = tx*tx + ty*ty
        rx(n) = half * tx / tn
        ry(n) = half * ty / tn
     enddo

     sx = zero
     sy = zero

     do n = nc, 1, -1
        saux = sx + xl
        sx = rx(n) * saux - ry(n) * sy
        sy = rx(n) * sy + ry(n) * saux
        xl = h * xl
     enddo

     wx = cc * sx
     wy = cc * sy
  else
     xh = y
     yh = x
     rx(1) = zero
     ry(1) = zero

     do n = 9, 1, -1
        tx = xh + n * rx(1)
        ty = yh - n * ry(1)
        tn = tx*tx + ty*ty
        rx(1) = half * tx / tn
        ry(1) = half * ty / tn
     enddo

     wx = cc * rx(1)
     wy = cc * ry(1)
  endif

  !     if (y .eq. zero) wx = exp(-x**2)
  if (yy .lt. zero) then
     wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx
     wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy
     if (xx .gt. zero) wy = -wy
  else
     if (xx .lt. zero) wy = -wy
  endif

end SUBROUTINE ccperrf

double precision function gauss_erf(x)
  use math_constfi, only : zero, one, two, half
  implicit none
  !---- returns the value of the Gauss error integral:
  !     1/sqrt(2*pi) Int[-inf, x] exp(-x**2/2) dx
  double precision, intent(IN) :: x

  double precision :: xx, re, im

  xx = x / sqrt(two)
  call ccperrf(zero,xx,re,im)
  gauss_erf = one - half * exp(-xx**2) * re

end function gauss_erf

SUBROUTINE twwmap(pos, orbit)
  use twissotmfi
  use twissafi
  use matrices, only : EYE
  use math_constfi, only : zero, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Save concatenated sectormap (kick, rmatrix, tmatrix)             *
  !     Input:                                                           *
  !     pos   (double)  position                                         *
  !     orbit (double)  current orbit                                    *
  !     Further input in twissotm.fi                                     *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: pos, orbit(6)

  integer :: i, k, l
  double precision :: sum1, sum2, ek(6)
  double precision, external :: get_value
  logical :: accmap, sectorpure

  sectorpure = get_value('twiss ','sectorpure ') .ne. zero

  do i = 1, 6
     sum2 = orbit(i)
     do k = 1, 6
        sum1 = zero
        do l = 1, 6
           sum1 = sum1 + stmat(i,k,l) * sorb(l)
        enddo
        sum2 = sum2 - (srmat(i,k) - sum1) * sorb(k)
        if(sectorpure) srmat(i,k) = srmat(i,k) - two * sum1
     enddo
     ek(i) = sum2
  enddo
  SORB = ORBIT

  !---  jluc: note that twiss can be called through madxn.c's (1)
  !     pro_twiss or (2) pro_embedded_twiss
  !     In the first case, sectorTableName is set by the user or
  !     by default. In the second case,
  !     the parameter is meaningless and the following sector data
  !     storage should be skipped.
  !     To this end, we test whether the passed sectorTableName is
  !     equal to 'dummy'
  if ( sectorTableName .ne. 'dummy') then
     !---- Output.
     call sector_out(sectorTableName, pos, ek, srmat, stmat)
  endif
  !     else should issue a warning albeit the
  !     test must be true
  !     write (99, '(g20.6)') pos
  !     write (99, '(6e16.8)') ek
  !     write (99, '(6e16.8)') srmat
  !     write (99, '(6e16.8)') stmat

  !---- re-initialize map.
  accmap = get_value('twiss ','sectoracc ') .ne. zero
  if (accmap) then
  !  print *,"accumulate sector maps"
  else
  !  print *,"clearing sector maps"
    SRMAT = EYE
    STMAT = zero
  endif

end SUBROUTINE twwmap

SUBROUTINE tmdpdg(ftrk,orbit,fmap,ek,re,te)
  use matrices, only : EYE
  use math_constfi, only : zero, one
  implicit none
  !----------------------------------------------------------------------*
  !     dipedge element                                                  *
  !     Purpose: computes the effect of the dipole edge                  *
  !     Input:                                                           *
  !     ftrk         (logical) if true, track orbit.                     *
  !     Output:                                                          *
  !     fmap         (logical) if true, element has a map.               *
  !     re(6,6)      (double)  transfer matrix.                          *
  !----------------------------------------------------------------------*
  logical, intent(IN)  :: ftrk
  logical, intent(OUT) :: fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: bvk, e1, h, hgap, fint, tilt, corr
  double precision :: ek0(6), rw(6,6), tw(6,6,6)

  double precision, external :: node_value

  !-- LD: modified to include corrections from Frank.
  EK0 = zero
  RW = EYE
  TW = zero

  bvk  = node_value('other_bv ')
  e1   = bvk * node_value('e1 ')
  h    = bvk * node_value('h ')
  hgap = node_value('hgap ')
  fint = node_value('fint ')
  tilt = node_value('tilt ')
  corr = (h + h) * hgap * fint

  fmap = h .ne. zero .and. ( e1 .ne. zero .or. corr .ne. zero )
  if (.not. fmap) return

  !---- Fringe fields effects computed from the TWISS routine tmfrng
  !     tmfrng returns the matrix elements rw(used) and tw(unused)
  !     No radiation effects as it is a pure thin lens with no lrad
  call tmfrng(.false.,h,zero,e1,zero,zero,corr,rw,tw)
  call tmcat1(.true.,ek,re,te,ek0,rw,tw,ek,re,te)

  !---- Apply tilt.
  if (tilt .ne. zero) then
     call tmtilt(.false.,tilt,ek,re,te)
  endif

  if (ftrk) then
     call tmtrak(ek,re,te,orbit,orbit)
  endif

end SUBROUTINE tmdpdg

SUBROUTINE tmsol_th(ftrk,orbit,fmap,ek,re,te)
  use twissbeamfi, only : deltap, radiate, gamma, arad
  use math_constfi, only : zero, one, two, three
  implicit none
  !     Stolen from trrun.F courtesy Alex Koschick
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for solenoid element, no centre option.            *
  !     Input:                                                           *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: ftrk, fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: bvk, sk, skl, sks, sksl, cosTh, sinTh, Q0, Q

  double precision, external :: node_value, get_value
  double precision :: elrad, rfac, kx, ky
  double precision :: pt, bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  fmap = .true.

  !---- Get solenoid parameters
  sksl    = node_value('ksi ')
  sks     = node_value('ks ')
  elrad  = node_value('lrad ')

  !---- BV flag
  bvk = node_value('other_bv ')
  sks = sks * bvk
  sksl = sksl * bvk

  !---- Set up strengths
  sk    = sks / two
  skl   = sksl / two
  Q0 = skl/(one+deltap)
  cosTh = cos(Q0)
  sinTh = sin(Q0)
  Q = -sk*Q0

  !---- Half radiation effect at entry.
  if (radiate .and. ftrk) then
    if(elrad .eq. zero) then
      call fort_warn('TWCPGO: ','Radiation effects ignored for solenoid '// &
                         'with l=0, lrad=0 and radiate=true')
    else
       kx = ((sk**2)*orbit(1)-sk*orbit(4))*elrad;
       ky = ((sk**2)*orbit(3)+sk*orbit(2))*elrad;
       rfac = (arad * gamma**3 / three) * (kx**2 + ky**2) / elrad;
       pt = orbit(6);
       bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
       f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
       orbit(2) = orbit(2) * f_damp_t;
       orbit(4) = orbit(4) * f_damp_t;
       orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
    endif
  endif

  !---- First-order terms.
  re(1,1) = cosTh
  re(2,2) = re(1,1)
  re(3,3) = re(1,1)
  re(4,4) = re(1,1)
  re(1,3) = sinTh
  re(2,4) = re(1,3)
  re(3,1) = -re(1,3)
  re(4,2) = -re(1,3)
  re(2,1) = Q*re(1,1)
  re(4,3) = Q*re(1,1)
  re(4,1) = -Q*re(1,3)
  re(2,3) = Q*re(1,3)
  re(1,2) = zero
  re(1,4) = zero
  re(3,2) = zero
  re(3,4) = zero

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

  !---- Half radiation effect at exit.
  if (radiate .and. ftrk) then
      if(elrad .eq. zero) then
      call fort_warn('TWCPGO: ','Radiation effects ignored for solenoid '// &
                         'with l=0, lrad=0 and radiate=true')
    else
       kx = ((sk**2)*orbit(1)-sk*orbit(4))*elrad;
       ky = ((sk**2)*orbit(3)+sk*orbit(2))*elrad;
       rfac = (arad * gamma**3 / three) * (kx**2 + ky**2) / elrad;
       pt = orbit(6);
       bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
       f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
       orbit(2) = orbit(2) * f_damp_t;
       orbit(4) = orbit(4) * f_damp_t;
       orbit(6) = orbit(6) * (one - rfac) - rfac / bet0;
    endif
  endif

end SUBROUTINE tmsol_th

SUBROUTINE tmnll(fsec,ftrk,orbit,fmap,ek,re,te)
  use math_constfi, only : zero, one, two, pi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for elliptical lens A.Valishev Oct.19,2010         *
  !     Input:                                                           *
  !     FSEC      (logical) If true, return second order terms.          *
  !     FTRK      (logical) If true, track orbit.                        *
  !     Input/output:                                                    *
  !     ORBIT(6)  (real)    Closed orbit.                                *
  !     Output:                                                          *
  !     FMAP      (logical) If true, element has a map.                  *
  !     EK(6)     (real)    Kick due to element.                         *
  !     RE(6,6)   (real)    Transfer matrix.                             *
  !     TE(6,6,6) (real)    Second-order terms.                          *
  !----------------------------------------------------------------------*
  logical :: fsec, ftrk, fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  double precision :: knll, cnll, k1, dd, u, v, dUu, dUv, dux, duy, dvx, dvy, x, y

  double precision, external :: node_value

  !---- Matrix
  fmap = .true.
  cnll=node_value('cnll ')
  knll=node_value('knll ')/cnll
  k1 = knll/cnll
  re(1,1) = one
  re(2,1) =-two*k1
  re(2,2) = one
  re(3,3) = one
  re(4,3) = two*k1
  re(4,4) = one
  re(5,5) = one
  re(6,6) = one

  !---- Second Order Terms
  !  if (fsec) then
  ! no second order terms
  !  endif

  if (.not.ftrk) return

  !---- Track orbit.
  x = orbit(1)/cnll
  y = orbit(3)/cnll

  u = 0.5*sqrt((x-1)**2+y**2)+0.5*sqrt((x+1)**2+y**2)
  v = 0.5*sqrt((x+1)**2+y**2)-0.5*sqrt((x-1)**2+y**2)

  if (u .eq. 1d0) then
     dd=0
  else
     dd=u**2*log(u+sqrt(u*u-1))/sqrt(u**2-1)
  endif

  dUu = (u+log(u+sqrt(u*u-1))*sqrt(u**2-1)+dd)/(u**2-v**2) &
       -2*u*(u*log(u+sqrt(u*u-1))*sqrt(u**2-1) &
       +v*(acos(v)-0.5*pi)*sqrt(1-v**2)) /(u**2-v**2)**2

  dUv = 2*v*(u*log(u+sqrt(u*u-1))*sqrt(u**2-1) &
       +v*(acos(v)-0.5*pi)*sqrt(1-v**2)) /(u**2-v**2)**2 &
       -(v-(acos(v)-0.5*pi)*sqrt(1-v**2)+v**2*(acos(v)-0.5*pi)/sqrt(1-v**2)) &
       /(u**2-v**2)

  dux = 0.5*(x-1)/sqrt((x-1)**2+y**2) +0.5*(x+1)/sqrt((x+1)**2+y**2)
  duy = 0.5*y/sqrt((x-1)**2+y**2) +0.5*y/sqrt((x+1)**2+y**2)
  dvx = 0.5*(x+1)/sqrt((x+1)**2+y**2) -0.5*(x-1)/sqrt((x-1)**2+y**2)
  dvy = 0.5*y/sqrt((x+1)**2+y**2) -0.5*y/sqrt((x-1)**2+y**2)

  orbit(2) = orbit(2) + knll*(dUu*dux+dUv*dvx)
  orbit(4) = orbit(4) + knll*(dUu*duy+dUv*dvy)

end SUBROUTINE tmnll

SUBROUTINE tmrfmult(fsec,ftrk,orbit,fmap,ek,re,te)
  use twtrrfi
  use twisslfi
  use twissbeamfi, only : radiate, deltap, beta, gamma, pc, arad
  use math_constfi, only : zero, one, two, three, ten3m, ten6p, pi, twopi
  use phys_constfi, only : clight
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for thin rf-multipole.                             *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  ! Strategy to implement BV-flag
  ! 0) apply bv-flag to voltage and multipole strengths (the inverse map has V = -V; K?N|S = -K?N|S)
  ! 1) track orbit(6) : just as in track: P o inverse(M) * P
  ! 2) ek, re, te : create the vector / matrix / tensor elements for
  !    the inverse map (see 0), then apply transformation P to each element

  logical :: fsec, ftrk, fmap
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  integer :: nord, i, j, nn, ns, n_ferr, dummyi, ii, jj, kk
  double precision :: elrad, rfac, bvk, tilt, angle
  double precision :: f_errors(0:maxferr)
  double precision :: normal(0:maxmul), skew(0:maxmul)
  double precision :: cangle, sangle, dtmp

  double precision :: krf, vrf
  double precision :: x, y, z, px, py, pt, dpx, dpy, dpt
  double precision :: freq, volt, lag, harmon
  double precision :: field_cos(2,0:maxmul), field_sin(2,0:maxmul)
  double precision :: pnl(0:maxmul), psl(0:maxmul)
  double precision :: P(6)
  double complex :: Cm2, Sm2, Cm1, Sm1, Cp0, Sp0, Cp1, Sp1

  integer, external :: node_fd_errors
  double precision, external :: node_value, get_value
  double complex, parameter :: icomp=(0d0,1d0) ! imaginary
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !---- Zero the arrays
  NORMAL = zero
  SKEW = zero
  PNL = zero
  PSL = zero
  F_ERRORS = zero

  !---- Read-in the parameters
  freq   = node_value('freq ');
  volt   = node_value('volt ');
  lag    = node_value('lag ');
  harmon = node_value('harmon ');
  bvk    = node_value('other_bv ')
  elrad  = node_value('lrad ')
  tilt   = node_value('tilt ')

  n_ferr = node_fd_errors(f_errors);
  call get_node_vector('knl ', nn, normal)
  call get_node_vector('ksl ', ns, skew)
  nord = max(nn, ns, n_ferr/2-1);

  call get_node_vector('pnl ', dummyi, pnl); ! NOTE !!!!! THIS DOES NOT MAKE USE OF NODE->PNL and NODE->PSL
  call get_node_vector('psl ', dummyi, psl);

  rfac = zero
  fmap = .true.

  !---- Set-up some parameters
  krf = 2 * pi * freq * ten6p/clight;
  vrf = bvk * volt * ten3m / ( pc * (one+deltap) );

  !---- Particle's coordinates
  if (ftrk) then
    ! if bvk = -1 apply the transformation P: (-1, 1, 1, -1, -1, 1) * X
    x  = orbit(1) * bvk;
    px = orbit(2);
    y  = orbit(3);
    py = orbit(4) * bvk;
    z  = orbit(5) * bvk;
    pt = orbit(6);
  else
    x  = zero;
    px = zero;
    y  = zero;
    py = zero;
    z  = zero;
    pt = zero;
  endif

  !---- Vector with strengths + field errors
  do i = 0, nord;
     field_cos(1,i) = bvk * (normal(i) * cos( pnl(i)*twopi - krf * z ) + f_errors(2*i))   / (one + deltap);
     field_sin(1,i) = bvk * (normal(i) * sin( pnl(i)*twopi - krf * z ) )                  / (one + deltap);
     field_cos(2,i) = bvk * (skew(i)   * cos( psl(i)*twopi - krf * z ) + f_errors(2*i+1)) / (one + deltap);
     field_sin(2,i) = bvk * (skew(i)   * sin( psl(i)*twopi - krf * z ) )                  / (one + deltap);

     if (tilt.ne.zero)  then
        angle = (i+1) * (-tilt);
        cangle = cos(angle);
        sangle = sin(angle);
        dtmp           = field_cos(1,i) * cangle - field_cos(2,i) * sangle;
        field_cos(2,i) = field_cos(1,i) * sangle + field_cos(2,i) * cangle;
        field_cos(1,i) = dtmp;
        dtmp           = field_sin(1,i) * cangle - field_sin(2,i) * sangle;
        field_sin(2,i) = field_sin(1,i) * sangle + field_sin(2,i) * cangle;
        field_sin(1,i) = dtmp;
     endif
  enddo

  !---- Prepare to calculate the kick and the matrix elements
  Cm2 = 0d0;
  Sm2 = 0d0;
  Cm1 = 0d0;
  Sm1 = 0d0;
  Cp0 = 0d0;
  Sp0 = 0d0;
  Cp1 = 0d0;
  Sp1 = 0d0;

  do i = nord, 0, -1
    if (i .ge. 2) then
      Cm2 = Cm2 * (x+icomp*y) / (i-1) + field_cos(1,i)+icomp*field_cos(2,i);
      Sm2 = Sm2 * (x+icomp*y) / (i-1) + field_sin(1,i)+icomp*field_sin(2,i);
    endif

    if (i .ge. 1) then
      Cm1 = Cm1 * (x+icomp*y) / i   + field_cos(1,i)+icomp*field_cos(2,i);
      Sm1 = Sm1 * (x+icomp*y) / i   + field_sin(1,i)+icomp*field_sin(2,i);
    endif

    Cp0 = Cp0 * (x+icomp*y) / (i+1)   + field_cos(1,i)+icomp*field_cos(2,i);
    Sp0 = Sp0 * (x+icomp*y) / (i+1)   + field_sin(1,i)+icomp*field_sin(2,i);
    Cp1 = Cp1 * (x+icomp*y) / (i+2)   + field_cos(1,i)+icomp*field_cos(2,i);
    Sp1 = Sp1 * (x+icomp*y) / (i+2)   + field_sin(1,i)+icomp*field_sin(2,i);
  enddo

  Sp1 = Sp1 * (x+icomp*y);
  Cp1 = Cp1 * (x+icomp*y);

  !---- Track orbit.
  if (ftrk) then

     !---- The kick
     dpx = -REAL(Cp0);
     dpy = AIMAG(Cp0);
     dpt =  vrf * sin(lag * twopi - krf * z) - krf * REAL(Sp1);

     !---- Radiation effects at entrance.
     if (radiate  .and.  elrad .ne. zero) then
        rfac = arad * gamma**3 * (dpx**2+dpy**2) / (three*elrad)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        px = px * f_damp_t;
        py = py * f_damp_t;
        pt = pt * (one - rfac) - rfac / bet0;
     endif

     !---- Apply the kick
     px = px + dpx
     py = py + dpy
     pt = pt + dpt

     !---- Radiation effects at exit.
     if (radiate  .and.  elrad .ne. zero) then
        px = px * f_damp_t;
        py = py * f_damp_t;
        pt = pt * (one - rfac) - rfac / bet0;
     endif

    ! apply the transformation P: (-1, 1, 1, -1, -1, 1) * X
    orbit(1) = x  * bvk;
    orbit(2) = px;
    orbit(3) = y;
    orbit(4) = py * bvk;
    orbit(5) = z  * bvk;
    orbit(6) = pt;

  endif

  !---- Element Kick
  !ek(2) = -REAL(Cp0);
  !ek(4) = AIMAG(Cp0);
  !ek(6) =  vrf * sin(lag * twopi - krf * z) - krf * REAL(Sp1);

  !---- First-order terms
  re(2,1) = -REAL(Cm1);
  re(2,3) =  AIMAG(Cm1);
  
  re(2,5) = -krf * REAL(Sp0);
  re(4,1) =  re(2,3);
  re(4,3) = -re(2,1);
  re(4,5) =  krf * AIMAG(Sp0);
  re(6,1) =  re(2,5);
  re(6,3) =  re(4,5);
  re(6,5) = -krf * vrf * cos(lag * twopi - krf * z) + krf * krf * REAL(Cp1);

  !---- Second-order terms (use X,Y from orbit tracking).
  if (fsec) then
     te(2,1,1) = 0.5 * (-REAL(Cm2));
     te(2,1,3) = 0.5 * ( AIMAG(Cm2));
     te(2,1,5) = 0.5 * (-krf * REAL(Sm1));
     te(2,3,1) = te(2,1,3);
     te(2,3,3) = 0.5 * ( REAL(Cm2));
     te(2,3,5) = 0.5 *   krf * AIMAG(Sm1);
     te(2,5,1) = te(2,1,5);
     te(2,5,3) = te(2,3,5);
     te(2,5,5) = 0.5 * krf * krf * REAL(Cp0);
     te(4,1,1) =  te(2,1,3);
     te(4,1,3) = -te(2,1,1);
     te(4,1,5) =  te(2,3,5);
     te(4,3,1) =  te(4,1,3);
     te(4,3,3) = -te(2,1,3);
     te(4,3,5) = -te(2,1,5);
     te(4,5,1) =  te(4,1,5);
     te(4,5,3) =  te(4,3,5);
     te(4,5,5) =  0.5 * (-krf * krf * AIMAG(Cp0));
     te(6,1,1) =  te(2,1,5);
     te(6,1,3) =  te(2,3,5);
     te(6,1,5) =  te(2,5,5);
     te(6,3,1) =  te(6,1,3);
     te(6,3,3) = -te(2,1,5);
     te(6,3,5) =  te(4,5,5);
     te(6,5,1) =  te(6,1,5);
     te(6,5,3) =  te(6,3,5);
     te(6,5,5) =  0.5 * (-krf * krf * vrf * sin(lag * twopi - krf * z) + krf * krf * krf * REAL(Sp1));
  endif

  ! Apply P tranformation to each index
  if (bvk .eq. -one) then
    P(1) = -one;
    P(2) =  one;
    P(3) =  one;
    P(4) = -one;
    P(5) = -one;
    P(6) =  one;
    do ii=1,6
      do jj=1,6
        do kk=1,6
          te(ii,jj,kk) = te(ii,jj,kk) * P(ii) * P(jj) * P(kk);
        enddo
        re(ii,jj) = re(ii,jj) * P(ii) * P(jj);
      enddo
      ek(ii) = ek(ii) * P(ii);
    enddo
  endif

end SUBROUTINE tmrfmult

SUBROUTINE calcsyncint(rhoinv,blen,k1,e1,e2,betxi,alfxi,dxi,dpxi,I)
  use math_constfi, only : zero, one, two
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Calculate synchrotron radiation integrals contribution of        *
  !     single element with parameters passed as input.                  *
  !     Method is implemented from SLAC-Pub-1193 where integration is    *
  !     done explicitly and includes effect of poleface rotations        *
  !                                                                      *
  !     Input:                                                           *
  !     rhoinv (double) inverse radius of curvature                      *
  !     blen (double) length of element                                  *
  !     k1 (double) gradient of element                                  *
  !     e1, e2 (double) pole face rotations at entrance and exit         *
  !     betxi, alfxi, dxi, dpxi (double) twiss parameters in x plane     *
  !                                                                      *
  ! ATTENTION:                                                           *
  ! because of the choice of canonical variables in MAD-X, (PT instead   *
  ! of DELTAP) the internal dispersion and dispersion derivatives must   *
  ! be multiplied by beta to match the functions quoted in litterature.  *
  ! This subroutine expects the dispersion as quoted in litterature and  *
  ! the calling routine must provide the dispersion and derivatives      *
  ! accordingly!                                                         *
  !                                                                      *
  !     Output:                                                          *
  !     I(5) (double) contributions of element to the 5 synchrotron      *
  !          radiation integrals.                                        *
  !                                                                      *
  !                                                                      *
  ! internal calculation is done with complex numbers to cater for both  *
  ! globally focusing and defocusing cases.                              *
  !                                                                      *
  ! Note: relies on implicit type promotion between real and complex for *
  !       k2.  dispaverage and curlyhaverage are explicitely converted   *
  !       to doubles.                                                    *
  !                                                                      *
  !                                                                      *
  !     Author: Ghislain Roy - June 2014                                 *
  !----------------------------------------------------------------------*

  double precision, intent(IN) :: rhoinv, blen, k1, e1, e2, betxi, alfxi, dxi, dpxi
  double precision, intent(OUT) :: I(8)

  ! local variables
  double precision :: dx2, gamx, dispaverage, curlyhaverage, lq
  double precision :: betx, alfx, dx, dpx, u0x, u1x, u2x
  double precision :: gammai, betxaverage, k1n
  double complex :: k2, k, kl

  integer, external :: get_option
  double precision, external :: node_value

  betx = betxi
  dx = dxi

  ! effect of poleface rotation
  alfx = alfxi - betxi*rhoinv*tan(e1)
  dpx = dpxi + dxi*rhoinv*tan(e1)

  gamx = (1+alfx**2)/betx

  ! global gradient combining weak focusing and dipole gradient
  ! k2 can be positive or negative and k can be real or imaginary
  k2 = rhoinv*rhoinv + 2*k1
  k = sqrt(k2)
  kl = k*blen

  ! propagation of dispersion at exit
  dx2 = real(dx*cos(kl) + dpx*sin(kl)/k + rhoinv*(1-cos(kl))/(k*k))

  dispaverage = real(dx * sin(kl)/kl &
             + dpx * (1 - cos(kl))/(k*kl) &
             + rhoinv * (kl - sin(kl))/(k2*kl))

  curlyhaverage = real( gamx*dx*dx + 2*alfx*dx*dpx + betx*dpx*dpx &
                + 2*rhoinv*blen*( -(gamx*dx + alfx*dpx)*(kl-sin(kl))/(kl*kl*k) &
                                 + (alfx*dx + betx*dpx)*(1-cos(kl))/(kl*kl)) &
                + blen*blen*rhoinv*rhoinv*( &
                       gamx*(3*kl - 4*sin(kl) + sin(kl)*cos(kl))/(2*k2*kl**3) &
                     - alfx*(1-cos(kl))**2/(k*kl**3) &
                     + betx*(kl-cos(kl)*sin(kl))/(2*kl**3)))
  if (rhoinv .ne. 0.d0) then
  I(1) = dispaverage * rhoinv * blen
  I(2) = rhoinv*rhoinv * blen
  I(3) = abs(rhoinv)**3 * blen
  I(4) = dispaverage*rhoinv*(rhoinv**2 + 2*k1) * blen &
           - rhoinv*rhoinv*(dx*tan(e1) + dx2*tan(e2))
  I(5) = curlyhaverage * abs(rhoinv)**3 * blen
  end if

  if(k1 .ne. 0.d0) then
    lq = node_value('l ')
    gammai = (one+alfxi*alfxi)/betxi;
    dx2 = real(dx*cos(kl) + dpx*sin(kl)/k)
    dispaverage = real(dx * sin(kl)/kl &
             + dpxi * (1 - cos(kl))/(k*kl))
  if(k1 .ge. zero) then
    k1n = k1
    u0x = (one + sin(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/two
    u1x = sin(sqrt(k1n)*lq)**two/(k1n*lq)
    u2x = (one - sin(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/(two*k1n)
    dx2 = cos(sqrt(k1n)*lq)*dxi + (one/sqrt(k1n))*sin(sqrt(k1n)*lq)*dpxi
    dispaverage = (dxi+dx2)/two
  else
    k1n = -k1
    u0x = (one + sinh(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/two
    u1x = sinh(sqrt(k1n)*lq)**two/(k1n*lq)
    u2x = -(one - sinh(two*sqrt(k1n)*lq)/(two*sqrt(k1n)*lq))/(two*k1n)
    dx2 = cosh(sqrt(k1n)*lq)*dxi + (one/sqrt(k1n))*sinh(sqrt(k1n)*lq)*dpxi
    dispaverage = (dxi+dx2)/two
  endif

    betxaverage = betxi*u0x - alfxi*u1x  + gammai*u2x

    I(6) = (k1n**2)*betxaverage*lq
    I(8) = (k1n**2)*dispaverage**2*lq

  endif

  if (get_option('debug ') .ne. 0) then
     print *, ' '
     print *, 'Input:  rhoinv = ', rhoinv, 'k1 = ', k1, 'e1 =', e1, 'e2 = ', e2, 'blen = ', blen
     print *, '        betxi = ', betxi, 'alfxi = ', alfxi, 'dxi = ', dxi, 'dpxi = ', dpxi
     print *, ' --> '
     print *, '        k2 = ', k2, '  k = ', k, 'k*l = ', kl
     print *, '        alfx = ', alfx, 'dpx = ', dpx, 'gamx = ', gamx, 'dx2 = ', dx2
     print *, '        dispaverage = ', dispaverage, 'curlyhaverage = ', curlyhaverage
     print *, 'Contributions to Radiation Integrals:', I(1), I(2), I(3), I(4), I(5)
     print *, ' '
  endif

END SUBROUTINE calcsyncint

subroutine tmfoc(el,sk1,c,s,d,f)
  use math_constfi, only : zero, one, two, six, twelve
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Compute linear focussing functions.                                *
  ! Input:                                                               *
  !   el        (double)  element length.                                *
  !   sk1       (double)  quadrupole strength.                           *
  ! Output:                                                              *
  !   c         (double)  cosine-like function.             c(k,l)       *
  !   s         (double)  sine-like function.               s(k,l)       *
  !   d         (double)  dispersion function.              d(k,l)       *
  !   f         (double)  integral of dispersion function.  f(k,l)       *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: el, sk1
  double precision, intent(OUT) :: c, s, d, f

  double precision :: qk, qkl, qkl2

  double precision, parameter :: twty=20d0, thty=30d0, foty2=42d0

  qk = sqrt(abs(sk1))
  qkl = qk * el
  qkl2 = sk1 * el**2

  if (abs(qkl2) .le. 1e-2) then
     c = (one - qkl2 * (one - qkl2 / twelve) /  two)
     s = (one - qkl2 * (one - qkl2 / twty) /  six) * el
     d = (one - qkl2 * (one - qkl2 / thty) / twelve) * el**2 / two
     f = (one - qkl2 * (one - qkl2 / foty2) / twty) * el**3 / six
  else
     if (qkl2 .gt. zero) then
        c = cos(qkl)
        s = sin(qkl) / qk
     else
        c = cosh(qkl)
        s = sinh(qkl) / qk
     endif

     d = (one - c) / sk1
     f = (el  - s) / sk1
  endif

end subroutine tmfoc

SUBROUTINE tmcrab(fsec,ftrk,orbit,fmap,el,ek,re,te)
  use twtrrfi
  use twisslfi
  use twissbeamfi, only : radiate, deltap, pc, beta, gamma, arad
  use math_constfi, only : zero, one, two, three, half, quarter, ten3m, ten3p, ten6p, pi, twopi
  use phys_constfi, only : clight
  use matrices, only : EYE
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for thin crab cavity.                             *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  ! Strategy to implement BV-flag
  ! 0) apply bv-flag to voltage and multipole strengths (the inverse map has V = -V; K?N|S = -K?N|S)
  ! 1) track orbit(6) : just as in track: P o inverse(M) * P
  ! 2) ek, re, te : create the vector / matrix / tensor elements for
  !    the inverse map (see 0), then apply transformation P to each element

  logical :: fsec, ftrk, fmap
  double precision :: el
  double precision :: orbit(6), ek(6), re(6,6), te(6,6,6)

  integer :: j, ii, jj, kk, dummyi, n_ferr
  double precision :: elrad, rfac, bvk, tilt, cangle, sangle, dtmp
  double precision :: f_errors(0:maxferr)
  double precision :: ed(6), rd(6,6), td(6,6,6)

  double precision :: krf
  double precision :: x, y, z, px, py, pt, dpx, dpy, dpt
  double precision :: freq, rfv, rfl, harmon
  double precision :: field_cos(2), field_sin(2)
  double precision :: kn0l, pn0
  double precision :: P(6)
  double complex :: Cp0, Sp0, Cp1, Sp1

  integer, external :: node_fd_errors
  double precision, external :: node_value, get_value
  double complex, parameter :: icomp=(0d0,1d0) ! imaginary
  double precision :: bet0, bet_sqr, f_damp_t

  bet0  =  get_value('beam ','beta ')

  !---- Zero the arrays
  F_ERRORS = zero
  !FIELD = zero
  TE = zero

  !---- drift matrix
  ED = zero
  RD = EYE
  TD = zero

  call tmdrf(fsec,ftrk,orbit,fmap,el/two,ed,rd,td);

  !---- Read-in the parameters
  harmon = node_value('harmon ');
  bvk = node_value('other_bv ')
  elrad = node_value('lrad ')
  tilt = node_value('tilt ')

  rfv  = node_value('volt ')
  freq = node_value('freq ')
  rfl  = node_value('lag ')

  kn0l = rfv / pc / ten3p; ! MeV / 1d3 / GeV = rad
  pn0  = quarter + rfl; ! 2pi/4 + rfl

  n_ferr = node_fd_errors(f_errors);

  rfac = zero
  fmap = .true.

  !---- Set-up some parameters
  krf = twopi * freq * ten6p/clight;

  !if (n_ferr .gt. 0) then
  !   call dcopy(f_errors, field, min(2, n_ferr))
  !endif

  !---- Particle's coordinates
  if (ftrk) then
    ! apply the transformation P: diag(-1, 1, 1, -1, -1, 1) * X
    x  = orbit(1) * bvk;
    px = orbit(2);
    y  = orbit(3);
    py = orbit(4) * bvk;
    z  = orbit(5) * bvk;
    pt = orbit(6);
  else
    x  = zero;
    px = zero;
    y  = zero;
    py = zero;
    z  = zero;
    pt = zero;
  endif

  !---- Vector with strengths + field errors
  field_cos(1) = bvk * (kn0l * cos(pn0 * twopi - krf * z) + f_errors(0)) / (one + deltap);
  field_sin(1) = bvk * (kn0l * sin(pn0 * twopi - krf * z))               / (one + deltap);
  field_cos(2) = zero;
  field_sin(2) = zero;
  if (tilt .ne. zero)  then
     cangle = cos(-tilt);
     sangle = sin(-tilt);

     dtmp         = field_cos(1) * cangle;
     field_cos(2) = field_cos(1) * sangle;
     field_cos(1) = dtmp;

     dtmp         = field_sin(1) * cangle - field_sin(2) * sangle;
     field_sin(2) = field_sin(1) * sangle + field_sin(2) * cangle;
     field_sin(1) = dtmp;
  endif

  !---- Prepare to calculate the kick and the matrix elements
  Cp0 = field_cos(1) + icomp*field_cos(2);
  Sp0 = field_sin(1) + icomp*field_sin(2);
  Cp1 = Cp0 * (x+icomp*y);
  Sp1 = Sp0 * (x+icomp*y);

  !---- Track orbit.
  if (ftrk) then

     !---- The kick
     dpx = -REAL(Cp0);
     dpy = AIMAG(Cp0);
     dpt = - krf * REAL(Sp1);

     !---- Radiation effects at entrance.
     if (radiate  .and.  elrad .ne. zero) then
        rfac = arad * gamma**3 * (dpx**2+dpy**2) / (three*elrad)
        bet_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / bet_sqr);
        px = px * f_damp_t;
        py = py * f_damp_t;
        pt = pt * (one - rfac) - rfac / bet0;
     endif

     !---- Apply the kick
     px = px + dpx
     py = py + dpy
     pt = pt + dpt

     !---- Radiation effects at exit.
     if (radiate  .and.  elrad .ne. zero) then
        px = px * f_damp_t;
        py = py * f_damp_t;
        pt = pt * (one - rfac) - rfac / bet0;
     endif

    ! apply the transformation P: diag(-1, 1, 1, -1, -1, 1) * X
    orbit(1) = x  * bvk;
    orbit(2) = px;
    orbit(3) = y;
    orbit(4) = py * bvk;
    orbit(5) = z  * bvk;
    orbit(6) = pt;

  endif

  !---- Element Kick
  ek(2) = -REAL(Cp0);
  ek(4) = AIMAG(Cp0);
  ek(6) = -krf * REAL(Sp1);

  !---- First-order terms
  re(2,5) = -krf * REAL(Sp0);
  re(4,5) =  krf * AIMAG(Sp0);
  re(6,1) =  re(2,5);
  re(6,3) =  re(4,5);
  re(6,5) =  krf * krf * REAL(Cp1);

  !---- Second-order terms
  if (fsec) then
     te(2,5,5) =  ( krf * krf * REAL(Cp0)) / two;
     te(4,5,5) =  (-krf * krf * AIMAG(Cp0)) / two;
     te(6,5,5) =  (-krf + krf * krf * krf * REAL(Sp1)) / two;
  endif

  ! Apply P tranformation to each index
  if (bvk .eq. -one) then
    P(1) = -one;
    P(2) =  one;
    P(3) =  one;
    P(4) = -one;
    P(5) = -one;
    P(6) =  one;
    do ii=1,6
      do jj=1,6
        do kk=1,6
          te(ii,jj,kk) = te(ii,jj,kk) * P(ii) * P(jj) * P(kk);
        enddo
        re(ii,jj) = re(ii,jj) * P(ii) * P(jj);
      enddo
      ek(ii) = ek(ii) * P(ii);
    enddo
  endif

  ! Add half a drift space before and after the Crab kick
  call tmcat1(fsec,ed,rd,td,ek,re,te,ek,re,te);
  call tmdrf(fsec,ftrk,orbit,fmap,el/two,ed,rd,td);
  call tmcat1(fsec,ek,re,te,ed,rd,td,ek,re,te);

end SUBROUTINE tmcrab
SUBROUTINE twcpin_print(rt,r0mat )
  use twiss0fi
  use twisscfi
  use twissbeamfi, only : deltap
  use matrices, only : EYE, SMAT
  use math_constfi, only : zero, one, two, quarter
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Print initial values for linear coupling parameters.                   *
  !     Entry point for mad_emit                                         *
  !     Input:                                                           *
  !     rt(6,6)      (double)  one turn transfer matrix.                 *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: rt(6,6), r0mat(2,2)
  double precision :: e(2,2),  f(2,2), edet, fdet, r0mat_bar(2,2)
  double precision :: e1(2,2),  f1(2,2),  e1det ,  f1det
  double precision :: e11(2,2), f11(2,2), e11det , f11det
  double precision :: et(2,2),  ft(2,2),  etdet ,  ftdet
  double precision :: em(2,2),  fm(2,2),  emdet ,  fmdet
  double precision :: e_12(2,2), e_21(2,2)
  double precision :: f_12(2,2), f_21(2,2)
  double precision :: e1_12(2,2), e1_21(2,2)
  double precision :: f1_12(2,2), f1_21(2,2)
  double precision :: e11_12(2,2), e11_21(2,2)
  double precision :: f11_12(2,2), f11_21(2,2)
  double precision :: et_12(2,2), et_21(2,2)
  double precision :: ft_12(2,2), ft_21(2,2)
  double precision :: em_12(2,2), em_21(2,2)
  double precision :: fm_12(2,2), fm_21(2,2)
  double precision :: a(2,2), b(2,2), c(2,2), d(2,2)
  double precision :: ra(4,4), ss(4,4)
  double precision :: u(4,4), v(4,4), vbar(4,4), vu (4,4), tmp(2,2), r_Det
  double precision :: gamma

  edet   = zero; fdet   = zero
  e1det  = zero; f1det  = zero
  e11det = zero; f11det = zero
  etdet  = zero; ftdet  = zero
  emdet  = zero; fmdet  = zero
  gamma  = one


  open (unit = 2, file = "afterclean_twcpin.out")
  write(2,*) "After clean fort "
  !-- simplecticity
  RA  = RT(1:4, 1:4)
  ss  = zero
  ss(1:2, 1:2)  = SMAT
  ss(3:4, 3:4)  = SMAT


  ! RT(1:4,1:4) = ( A , B
  !                 C , D)
  A  = RT(1:2,1:2) ; B  = RT(1:2,3:4)
  C  = RT(3:4,1:2) ; D  = RT(3:4,3:4)

  R0MAT_BAR =  matmul(matmul(-SMAT,transpose(R0MAT)),SMAT)
  r_det = r0mat(1,1) * r0mat(2,2) - r0mat(1,2) * r0mat(2,1)
  gamma = 1 /sqrt(1 + r_det)

  write(2,*) "check symplecticity Rt(1:4,1:4) " , matmul(transpose(RA),matmul(SS,RA)) - SS
  write(2,*) "R0MAT   = ", R0MAT
  write(2,*) "gammacp = ", gamma

  e   = zero
  f   = zero
  e1  = zero
  f1  = zero
  e11 = zero
  f11 = zero
  et  = zero
  et  = zero
  em  = zero
  em  = zero

  !--MAD-X Legacy
  e = (A - matmul(B,R0MAT))
  f = (D + matmul(R0MAT, B))
  edet = e(1,1) * e(2,2) - e(1,2) * e(2,1)
  fdet = f(1,1) * f(2,2) - f(1,2) * f(2,1)

  ! --Compute matrix V ("rotation", gamma*[ I R0MAT_BAR ; -R0MAT I]
  V(1:4, 1:4) = zero
  V(1:2, 1:2) = EYE(1:2, 1:2)
  V(3:4, 3:4) = EYE(3:4, 3:4)
  V(1:2, 3:4) =  R0MAT_BAR
  V(3:4, 1:2) = -R0MAT

  ! --Compute matrix inv(V)= VBAR
  VBAR(1:4, 1:4) = zero
  VBAR(1:2, 1:2) = EYE(1:2, 1:2)
  VBAR(3:4, 3:4) = EYE(3:4, 3:4)
  VBAR(1:2, 3:4) = -R0MAT_BAR
  VBAR(3:4, 1:2) = R0MAT

  ! --Compute VU = MV
  VU = matmul(RT(1:4,1:4),V)
  write(2,*) "MADX   E=VU(1:2,1:2)        " , VU(1:2, 1:2)
  write(2,*) "MADX   F=VU(3:4,3:4)        " , VU(3:4, 3:4)
  write(2,*) "MADX -RE=VU(3:4,1:2)        " , VU(3:4, 1:2)
  write(2,*) "MADX   E=VU(1:2,1:2)/(-R)   " , -matmul(R0MAT_BAR,VU(3:4, 1:2))/r_det
  write(2,*) "det V = ", v(1,1) * v(2,2) - v(1,2) * v(2,1)

     ! --Compute uncoupled block-diagonal U = g**2*[E 0, F 0]
     U = zero
     U = matmul(VBAR,VU)/(1 + r_det)

     !---- Find diagonal blocks.
     EM = U(1:2, 1:2)
     FM = U(3:4, 3:4)

     write(2,*) "MADX E=U(1:2,1:2) " , U(1:2, 1:2)
     write(2,*) "MADX F=U(3:4,3:4) " , U(1:2, 1:2)
     write(2,*) "MADX 0=U(1:2,3:4) " , U(1:2, 3:4)
     write(2,*) "MADX 0=U(3:4,1:2) " , U(3:4, 1:2)

     tmp = matmul(D,R0MAT) - C
     e1 = matmul(R0MAT_BAR, TMP)/r_det ! another way to express E =

     e11 = VU(1:2, 1:2)
     f11 = VU(3:4, 3:4)

     ! U = VBAR*M*V (Talman)
     !A = gamma2*(A        + B*RA - RD*C - RD*D*RA)(talman) = gamma2*(A        -B*R - RBAR*C + RBAR*D*R) (madx)
     !D = gamma2*(-RA*A*RD - RA*B + C*RD + D )     (talman) = gamma2*(R*A*RBAR +R*B + C*RBAR + D)        (madx)
     et = gamma**2*(A - matmul(B,R0MAT) - matmul(R0MAT_BAR,C) + matmul(R0MAT_BAR,matmul(D,R0MAT)))
     ft = gamma**2*(D + matmul(R0MAT,B) + matmul(C, R0MAT_BAR) + matmul(R0MAT,matmul(A,R0MAT_BAR)))

     write(2,*) "Talman 0=U(1:2,3:4) ", U(1:2, 3:4)
     write(2,*) "Talman 0=U(3:4,1:2) ", U(3:4,1:2)

     write(2,*) "E legacy (E=A-BR)         " , e
     write(2,*) "E (E1=RBAR(DR - C)/||R||) " , e1
     write(2,*) "E (E11 from VU=MV)        " , e11
     write(2,*) "E (Et)(talman)            " , et
     write(2,*) "E (Em)(E from E =VBARMV)  " , em

end SUBROUTINE twcpin_print
SUBROUTINE twcptk_print(re,r0mat, e, f)
  use twiss0fi
  use twisscfi
  use twissbeamfi, only : deltap
  use matrices, only : EYE, SMAT
  use math_constfi, only : zero, one, two, quarter
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Print initial values for linear coupling parameters.                   *
  !     Entry point for mad_emit                                         *
  !     Input:                                                           *
  !     rt(6,6)      (double)  one turn transfer matrix.                 *
  !     Output:                                                          *
  !     disp0        (double)  initial dispersion vector                 *
  !     r0mat(2,2)   (double)  coupling matrix                           *
  !     eflag        (integer) error flag.                               *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: re(6,6), r0mat(2,2), e(2,2), f(2,2)
  double precision ::  edet, fdet, r0mat_bar(2,2)
  double precision :: a(2,2), b(2,2), c(2,2), d(2,2)
  double precision :: ra(4,4), ss(4,4)

  edet = zero; fdet = zero

  open (unit = 3, file = "afterclean_twcptk.out")
  write(3,*) "After clean fort tk "
  !-- simplecticity
  RA  = RE(1:4, 1:4)
  ss  = zero
  ss(1:2, 1:2)  = SMAT
  ss(3:4, 3:4)  = SMAT

  R0MAT_BAR =  matmul(matmul(-SMAT,transpose(R0MAT)),SMAT)

  write(3,*) "check symplecticity RE(1:4,1:4) " , matmul(transpose(RA),matmul(SS,RA)) - SS
  write(3,*) "R0MAT = ", R0MAT

  edet = e(1,1) * e(2,2) - e(1,2) * e(2,1)
  fdet = f(1,1) * f(2,2) - f(1,2) * f(2,1)

  write(3,*) "EDET tk =       " , edet
  write(3,*) "FDET tk =       " , fdet


end SUBROUTINE twcptk_print

