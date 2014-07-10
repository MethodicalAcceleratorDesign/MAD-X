!     FS 15.03.2004: correcting warning printout
SUBROUTINE twiss(rt,disp0,tab_name,sector_tab_name)

  use twiss0fi
  use twissafi
  use twisslfi
  use twisscfi
  use twissotmfi
  use trackfi
  use fasterror
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TWISS command: Track linear lattice parameters.                  *
  !----------------------------------------------------------------------*
  logical fast_error_func
  integer i,ithr_on
  integer tab_name(*),chrom,eflag,get_option,izero,ione
  double precision rt(6,6),disp0(6),orbit0(6),orbit(6),tt(6,6,6), &
       ddisp0(6),r0mat(2,2),zero,one,two,get_value
  integer sector_tab_name(*) ! holds sectormap data
  character(48) charconv

  parameter (zero=0d0, one=1d0, two=2d0)

  data izero, ione / 0, 1 /

  !---- Initialization
  table_name = charconv(tab_name)
  sectorTableName = charconv(sector_tab_name)
  chrom=0
  eflag=0
  ithr_on=0
  fsecarb=.false.
  i = 6

  call dzero(orbit0,6)
  call get_node_vector('orbit0 ', i, orbit0)
  call m66one(rt)
  call m66one(rw)
  call dzero(tt,216)
  call get_disp0(disp0)
  call dzero(ddisp0,6)
  call dzero(r0mat,4)
  call dzero(opt_fun0,fundim)
  call dzero(opt_fun,fundim)
  call dzero(disp,6)
  call dzero(ddisp,6)
  call dzero(rmat,4)

  betx=zero
  alfx=zero
  amux=zero
  bety=zero
  alfy=zero
  amuy=zero
  bxmax=zero
  dxmax=zero
  bymax=zero
  dymax=zero
  xcomax=zero
  ycomax=zero
  sigxco=zero
  sigyco=zero
  sigdx=zero
  sigdy=zero
  wgt=zero
  cosmux=zero
  cosmuy=zero
  wx=zero
  phix=zero
  dmux=zero
  wy=zero
  phiy=zero
  dmuy=zero
  synch_1=zero
  synch_2=zero
  synch_3=zero
  synch_4=zero
  synch_5=zero
  suml=zero
  circ=zero
  eta=zero
  alfa=zero
  gamtr=zero
  qx=zero
  qy=zero
  sinmux=zero
  sinmuy=zero
  xix=zero
  xiy =zero

  !---- Track chromatic functions
  chrom=get_option('twiss_chrom ')
  
  !---- flag if called from match process
  !---- get match flag for storing variables in nodes
  match_is_on = get_option('match_is_on ') .ne. 0

  !---- flags for writing cumulative or lumped matrices
  rmatrix=get_value('twiss ','rmatrix ').ne.zero
  sectormap=get_option('twiss_sector ').ne.zero

  !---- Get circumference
  circ=get_value('probe ','circ ')
  if (circ.eq.zero) call aafail('TWISS: ', 'Zero length sequence.')

  !---- Set fast_error_func flag to use faster error function
  !---- including tables. Thanks to late G. Erskine
  fast_error_func = get_option('fast_error_func ') .ne. 0
  if(fast_error_func.and..not.fasterror_on) then
     call wzset
     fasterror_on = .true.
  endif

  if (get_option('twiss_inval ') .ne. 0) then   
     !---- Initial values from command attributes.
     call twinifun(opt_fun0, rt)
     if (get_option('twiss_print ') .ne. 0)  then
        print *, ' '
        print '(''open line - error with deltap: '',1p,e14.6)', get_value('probe ','deltap ')
        print '(''initial orbit vector: '', 1p,6e14.6)', orbit0
     endif
     call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,0,0,ithr_on)
     if(eflag.ne.0) go to 900
     if (get_option('twiss_print ') .ne. 0)  then
        print '(''final orbit vector:   '', 1p,6e14.6)', orbit
     endif
  else 
     !---- Initial values from periodic solution.
     call tmclor(orbit0,.true.,.true.,opt_fun0,rt,tt,eflag)
     if(eflag.ne.0) go to 900
     call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,0,0,ithr_on)
     if(eflag.ne.0) go to 900
     call twcpin(rt,disp0,r0mat,eflag)
     if(eflag.ne.0) go to 900
     !---- Initialize opt_fun0
     call twinifun(opt_fun0,rt)
  endif

  if (sectormap)  then
     call dcopy(orbit0, sorb, 6)
     call m66one(srmat)
     call dzero(stmat, 216)
  endif

  !---- Build table of lattice functions, coupled.
  call twcpgo(rt,orbit0)

  !---- List chromatic functions.
  if (chrom.ne.0) then
     call twbtin(rt,tt)
     call twchgo
  endif

  !---- Print summary
  if(get_option('twiss_summ ') .ne. 0) call tw_summ(rt,tt)

  if (get_option('keeporbit ') .ne. 0)  then
     i = 6
     call store_node_vector('orbit0 ', i, opt_fun0(9))
  endif

  call set_option('twiss_success ', ione)
  return

900 call set_option('twiss_success ', izero)

 end SUBROUTINE twiss



SUBROUTINE tmrefe(rt)

  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Transfer matrix w.r.t. ideal orbit for one period.               *
  !     Ignores cavities, radiation, and imperfections.                  *
  !     Output:                                                          *
  !     rt(6,6) (double) transfer matrix.                                *
  !----------------------------------------------------------------------*
  integer eflag,ithr_on
  double precision orbit0(6),orbit(6),rt(6,6),tt(6,6,6)

  ithr_on=0
  call dzero(orbit0,6)
  !---- Get transfer matrix.
  call tmfrst(orbit0,orbit,.false.,.false.,rt,tt,eflag,0,0,ithr_on)
end SUBROUTINE tmrefe
SUBROUTINE tmrefo(kobs,orbit0,orbit,rt)


  use twiss0fi
  implicit none
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Transfer matrix w.r.t. ideal orbit for one period.               *
  !     Input:                                                           *
  !     kobs    if > 0, track until node with this obs. point number     *
  !     Output:                                                          *
  !     orbit0(6) (double) closed orbit at start=end                     *
  !     orbit(6)  (double) closed orbit at obs. point kobs, or at end    *
  !     rt(6,6) (double) transfer matrix.                                *
  !----------------------------------------------------------------------*
  integer eflag,kobs,izero,ione,ithr_on
  double precision opt_fun0(fundim),orbit0(6),orbit(6),rt(6,6),     &
       tt(6,6,6)
  data izero, ione / 0, 1 /

  ithr_on=0
  call dzero(orbit0,6)
  !---- Get closed orbit and coupled transfer matrix.
  call tmclor(orbit0,.true.,.true.,opt_fun0,rt,tt,eflag)
  call set_option('bbd_flag ', ione)
  call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,kobs,0,ithr_on)
  call set_option('bbd_flag ', izero)
end SUBROUTINE tmrefo
SUBROUTINE twinifun(opt_fun0,rt)

  use twiss0fi
  use twisslfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial twiss parameters put to opt_fun0.                        *
  !     Input/output:                                                    *
  !     opt_fun0(fundim) (double) initial optical values:                *
  !     betx0,alfx0,amux0,bety0,alfy0,amuy0, etc.                        *
  !----------------------------------------------------------------------*
  integer i1,i2
  double precision opt_fun0(*),rt(6,6),betx,alfx,mux,bety,alfy,muy, &
       x,px,y,py,t,pt,dx,dpx,dy,dpy,wx,phix,dmux,wy,phiy,dmuy,ddx,ddpx,  &
       ddy,ddpy,r(2,2),energy,get_value,zero,twopi,get_variable
  parameter(zero=0d0)

  !---- Initialize
  twopi=get_variable('twopi ')
  betx=get_value('twiss ','betx ')
  bety=get_value('twiss ','bety ')

  if(betx.gt.zero) opt_fun0(3)=betx
  if(bety.gt.zero) opt_fun0(6)=bety

  if(opt_fun0(3).le.zero .or. opt_fun0(6).le.zero) &
       call aafail('TWINIFUN: ', 'BETX and BETY must be both larger than zero.')

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
  r(1,1)=get_value('twiss ','r11 ')
  r(1,2)=get_value('twiss ','r12 ')
  r(2,1)=get_value('twiss ','r21 ')
  r(2,2)=get_value('twiss ','r22 ')
  energy=get_value('probe ','energy ')

  if(alfx  .ne.zero) opt_fun0(4 )=alfx
  if(mux   .ne.zero) opt_fun0(5 )=mux * twopi
  if(alfy  .ne.zero) opt_fun0(7 )=alfy
  if(muy   .ne.zero) opt_fun0(8 )=muy * twopi
  if(x     .ne.zero) opt_fun0(9 )=x
  if(px    .ne.zero) opt_fun0(10)=px
  if(y     .ne.zero) opt_fun0(11)=y
  if(py    .ne.zero) opt_fun0(12)=py
  if(t     .ne.zero) opt_fun0(13)=t
  if(pt    .ne.zero) opt_fun0(14)=pt
  if(dx    .ne.zero) opt_fun0(15)=dx
  if(dpx   .ne.zero) opt_fun0(16)=dpx
  if(dy    .ne.zero) opt_fun0(17)=dy
  if(dpy   .ne.zero) opt_fun0(18)=dpy
  if(wx    .ne.zero) opt_fun0(19)=wx
  if(phix  .ne.zero) opt_fun0(20)=phix * twopi
  if(dmux  .ne.zero) opt_fun0(21)=dmux * twopi
  if(wy    .ne.zero) opt_fun0(22)=wy
  if(phiy  .ne.zero) opt_fun0(23)=phiy * twopi
  if(dmuy  .ne.zero) opt_fun0(24)=dmuy * twopi
  if(ddx   .ne.zero) opt_fun0(25)=ddx
  if(ddpx  .ne.zero) opt_fun0(26)=ddpx
  if(ddy   .ne.zero) opt_fun0(27)=ddy
  if(ddpy  .ne.zero) opt_fun0(28)=ddpy
  if(r(1,1).ne.zero) opt_fun0(29)=r(1,1)
  if(r(1,2).ne.zero) opt_fun0(30)=r(1,2)
  if(r(2,1).ne.zero) opt_fun0(31)=r(2,1)
  if(r(2,2).ne.zero) opt_fun0(32)=r(2,2)
  if(energy.ne.zero) opt_fun0(33)=energy
  if(rmatrix) then
     do i1=1,6
        do i2=1,6
           opt_fun0(33+(i1-1)*6+i2)=rt(i1,i2)
        enddo
     enddo
  endif

end SUBROUTINE twinifun
SUBROUTINE twprep(save,case,opt_fun,position)

  use twissafi
  use twisslfi
  use twissotmfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Finalize twiss parameters, if save flag fill                     *
  !     table with twiss parameters.                                     *
  !     Input:                                                           *
  !     case        (integer) =1 fill from twcpgo; =2 fill from twchgo   *
  !     position    (double)  end position of element                    *
  !     Input/output:                                                    *
  !     opt_fun(fundim) (double) optical values:                         *
  !     betx,alfx,amux,bety,alfy,amuy, etc.                              *
  !----------------------------------------------------------------------*

  ! 2014-May-15  15:40:41  ghislain: suppressed flag
  !     that was always set to 1 in calling functions and never modified.
  !     flag        (integer) fill flag: 0 no, !=0 yes                   *


  integer save,case,i
  double precision opt_fun(*),position,twopi,opt5,opt8,opt20,opt21, &
       opt23,opt24,get_variable,zero
  parameter(zero=0d0)

  !---- Initialize
  twopi=get_variable('twopi ')

  if(case.eq.1) then
     !--- fill with data from twcpgo (Twiss Couple)
     opt_fun(2)=position
     opt5 = opt_fun(5)
     opt_fun(5)= opt_fun(5) / twopi
     opt8 = opt_fun(8)
     opt_fun(8)= opt_fun(8) / twopi
     if(save.ne.0) call twfill(case,opt_fun,position)
     if (match_is_on)  call copy_twiss_data(opt_fun)
     opt_fun(5)= opt5
     opt_fun(8)= opt8
  elseif(case.eq.2) then
     !--- fill with data from twchgo (Twiss Chrom)
     opt20 = opt_fun(20)
     opt_fun(20)= opt_fun(20) / twopi
     opt21 = opt_fun(21)
     opt_fun(21)= opt_fun(21) / twopi
     opt23 = opt_fun(23)
     opt_fun(23)= opt_fun(23) / twopi
     opt24 = opt_fun(24)
     opt_fun(24)= opt_fun(24) / twopi
     if(save.ne.0) call twfill(case,opt_fun,position)
     if (match_is_on)  call copy_twiss_data(opt_fun)
     opt_fun(20)= opt20
     opt_fun(21)= opt21
     opt_fun(23)= opt23
     opt_fun(24)= opt24
  endif
end SUBROUTINE twprep
SUBROUTINE twfill(case,opt_fun,position)

  use twissafi
  use twisslfi
  use twissotmfi
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

  ! 2014-May-15  15:40:41  ghislain: suppressed flag
  !     that was always set to 1 in calling functions and never modified.
  !     flag        (integer) fill flag: 0 no, !=0 yes                   *

  integer case,i
  double precision opt_fun(*),position,twopi,opt5,opt8,opt20,opt21,opt23,opt24,zero,get_value
  parameter(zero=0d0)

  ripken=get_value('twiss ','ripken ').ne.zero

  if(case.eq.1) then
     i = 17
     call vector_to_table_curr(table_name, 's ', opt_fun(2), i)
     i = 5
     call vector_to_table_curr(table_name, 'r11 ', opt_fun(29), i)
     i = 5
     call vector_to_table_curr(table_name, 'kmax ', opt_fun(70), i)
     if(rmatrix) then
        i = 36
        call vector_to_table_curr(table_name, 're11 ', opt_fun(34), i)
     endif
     if(ripken) call twfill_ripken(opt_fun)
  elseif(case.eq.2) then
     i = 10
     call vector_to_table_curr(table_name, 'wx ', opt_fun(19), i)
  endif

  !---- Augment table twiss
  call augment_count(table_name)

end SUBROUTINE twfill
SUBROUTINE twfill_ripken(opt_fun)

  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Fill twiss table with Ripken-Mais twiss parameters.              *
  !     beta11, beta12, beta21, beta22,                                  *
  !     alfa11, alfa12, alfa21, alfa22,                                  *
  !     gama11, gama12, gama21, gama22,                                  *
  !----------------------------------------------------------------------*

  double precision opt_fun(*)
  double precision kappa, betx, bety, alfx, alfy, gamx, gamy, r11, r12, r21, r22 
  double precision beta11, beta12, beta21, beta22, & 
                   alfa11, alfa12, alfa21, alfa22, &
                   gama11, gama12, gama21, gama22

  betx = opt_fun(3)
  bety = opt_fun(6)
  alfx = opt_fun(4) 
  alfy = opt_fun(7) 
  r11 =  opt_fun(29)
  r12 =  opt_fun(30)
  r21 =  opt_fun(31)
  r22 =  opt_fun(32)

  kappa = 1./(1. + (r11*r22-r12*r21))

  gamx = (1.+alfx**2)/betx
  gamy = (1.+alfy**2)/bety

  beta11 = kappa * betx
  beta22 = kappa * bety
  beta12 = kappa * ( r22**2*bety + 2*r12*r22*alfy + r12**2*gamy )
  beta21 = kappa * ( r11**2*betx - 2*r12*r11*alfx + r12**2*gamx )
  
  alfa11 =  kappa * alfx
  alfa22 =  kappa * alfy
  alfa12 =  kappa * ( r21*r22*bety + (r12*r21 + r11*r22)*alfy + r11*r12*gamy )
  alfa21 = -kappa * ( r21*r11*betx - (r12*r21 + r11*r22)*alfx + r12*r22*gamx ) 

  gama11 = kappa * gamx
  gama22 = kappa * gamy
  gama12 = 0.
  if (beta12 .ne. 0.) gama12 = ((1. - kappa)**2 + alfa12**2) / beta12
  gama21 = 0.
  if (beta21 .ne. 0.) gama21 = ((1. - kappa)**2 + alfa21**2) / beta21
   
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

  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Find closed orbit for a beam line sequence.                      *
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
  logical fsec,ftrk,m66sta,pflag
  integer eflag,i,k,irank,itra,itmax,get_option,save_opt,thr_on,ithr_on
  parameter(itmax=20)
  double precision guess(6),opt_fun0(*),rt(6,6),tt(6,6,6),cotol,err,&
       orbit0(6),orbit(6),a(6,7),b(4,5),as(3,4),bs(2,3),deltap,get_value,&
       zero,one,get_variable
  parameter(zero=0d0,one=1d0)
  equivalence(a(1,1),b(1,1),as(1,1),bs(1,1))

  !---- Initialize.
  ithr_on=0
  thr_on = get_option('threader ')
  pflag = get_option('twiss_print ') .ne. 0
  deltap = get_value('probe ','deltap ')
  cotol = get_variable('twiss_tol ')
  eflag = 0
  !---- Initialize guess.
  call dcopy(guess,orbit0,6)

  !---- Iteration for closed orbit.
  iterate: do itra = 1, itmax

     !---- Track orbit and transfer matrix.
     call tmfrst(orbit0,orbit,fsec,ftrk,rt,tt,eflag,0,0,thr_on)
     if (eflag.ne.0)  return

     ! 2014-Mar-28  20:41:21  ghislain: turn off the threader immediately after first iteration, ie first turn
     thr_on = 0

     if (.not.m66sta(rt)) then 
        !---- Solve for dynamic case.
        err = zero
        call dcopy(rt,a,36)
        do i = 1, 6
           a(i,i) = a(i,i) - one
           a(i,7) = orbit(i) - orbit0(i)
           err = max(abs(a(i,7)), err)
        enddo
        call solver(a,6,1,irank)
        if (irank.lt.6) then
           print *, 'Singular matrix occurred during closed orbit search.'
           eflag = 1
           return
        endif
        do i = 1, 6
           orbit0(i) = orbit0(i) - a(i,7)
        enddo
        
     else                      
        !---- Solve for static case.
        err = zero
        do i = 1, 4
           do k = 1, 4
              b(i,k) = rt(i,k)
           enddo
           b(i,i) = b(i,i) - one
           b(i,5) = orbit(i) - orbit0(i)
           err = max(abs(b(i,5)), err)
        enddo
        call solver(b,4,1,irank)
        if (irank.lt.4) then
           print *, 'Singular matrix occurred during closed orbit search.'
           eflag = 1
           return
        endif
        do i = 1, 4
           orbit0(i) = orbit0(i) - b(i,5)
        enddo
     endif

     !---- Message and convergence test.

     if (pflag)  then
        print *, ' '
        print '(''iteration: '',i3,'' error: '',1p,e14.6,'' deltap: '',1p,e14.6)',itra,err,deltap
        print '(''orbit: '', 1p,6e14.6)', orbit0
     endif

     if (err.lt.cotol) then
        save_opt=get_option('keeporbit ')
        call tmfrst(orbit0,orbit,.true.,.true.,rt,tt,eflag,0,save_opt,ithr_on)
        opt_fun0(9 )=orbit0(1)
        opt_fun0(10)=orbit0(2)
        opt_fun0(11)=orbit0(3)
        opt_fun0(12)=orbit0(4)
        opt_fun0(13)=orbit0(5)
        opt_fun0(14)=orbit0(6)
        guess(1)=orbit0(1)
        guess(2)=orbit0(2)
        guess(3)=orbit0(3)
        guess(4)=orbit0(4)
        guess(5)=orbit0(5)
        guess(6)=orbit0(6)
        return ! normal exit
     endif

  enddo iterate

  !---- No convergence.
  print '(''Closed orbit did not converge in '', i3, '' iterations'')', itmax
  opt_fun0(9 )=zero
  opt_fun0(10)=zero
  opt_fun0(11)=zero
  opt_fun0(12)=zero
  opt_fun0(13)=zero
  opt_fun0(14)=zero
  return

end SUBROUTINE tmclor

SUBROUTINE tmfrst(orbit0,orbit,fsec,ftrk,rt,tt,eflag,kobs,save,thr_on)
  use bbfi
  use twiss0fi
  use name_lenfi
  use twisscfi
  use spch_bbfi
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
  logical fsec,ftrk,fmap
  character(28) tmptxt1, tmptxt2, tmptxt3
  character(2) ptxt(2)
  character(name_len) c_name(2), p_name
  integer eflag,j,code,restart_sequ,advance_node,node_al_errors,    &
       n_align,kobs,nobs,node,save,thr_on,ccode,pcode,get_vector,err,old,&
       kpro,corr_pick(2),enable,coc_cnt(2),lastnb,rep_cnt(2),max_rep,    &
       poc_cnt,get_option,debug
  double precision orbit0(6),orbit(6),orbit2(6),rt(6,6),tt(6,6,6),el,ek(6),   &
       re(6,6),te(6,6,6),al_errors(align_max),betas,gammas,node_value,   &
       get_value,parvec(26),orb_limit,zero,vector(10),reforb(6),         &
       restsum(2),restorb(6,2),restm(6,6,2),restt(6,6,6,2),cmatr(6,6,2), &
       pmatr(6,6),dorb(6),cick
  parameter(orb_limit=1d1,zero=0d0,ccode=15,pcode=18)
  parameter(max_rep=100)
  data ptxt / 'x-','y-'/

  debug = get_option('debug ')

  !---- Initialize
  !---- corr_pick stores for both projection the last pickup used by the
  !     threader in order to avoid corrections in front of it when
  !     restarting.

  do j = 1, 2
     corr_pick(j) = 0
     rep_cnt(j) = 0
     restsum(j) = zero
  enddo

  if (thr_on .gt. 0)  then
     call dzero(vector,10)
     j = get_vector('threader ', 'vector ', vector)
     if (j .lt. 3) thr_on = 0
  endif
 
  betas = get_value('probe ','beta ')
  gammas= get_value('probe ','gamma ')

  call dzero(tt,216)
  call m66one(rt)

  eflag = 0
  suml = zero

  call dcopy(orbit0,orbit,6)
  
  parvec(5)=get_value('probe ', 'arad ')
  parvec(6)=get_value('probe ', 'charge ') * get_value('probe ', 'npart ')
  parvec(7)=get_value('probe ', 'gamma ')
  bbd_cnt=0
  bbd_flag=1
  i_spch=0

  !---  start
  node = restart_sequ()

  !---  loop over nodes
10 continue

  bbd_pos=node !--- for space charge 

  code = node_value('mad8_type ')
  
  if(code.eq.39) code=15 ! tkicker treated as kicker
  if(code.eq.38) code=24 ! placeholder treated as instrument
  
  if (code .ge. ccode-1 .and. code .le. ccode+1)  then
     !---  kicker (code 14 to 16: hkicker, kicker, vkicker)
     if (thr_on .gt. 0)  then
        !---  threader is on - keep position,orbit,matrix,and tensor for restart
        if (code .le. ccode) then !--- kicker or hkicker
           restsum(1) = suml
           call dcopy(orbit,restorb(1,1),6)
           call dcopy(rt,restm(1,1,1),36)
           call dcopy(tt,restt(1,1,1,1),216)
        endif
        if (code .ge. ccode) then !--- kicker or vkicker
           restsum(2) = suml
           call dcopy(orbit,restorb(1,2),6)
           call dcopy(rt,restm(1,1,2),36)
           call dcopy(tt,restt(1,1,1,2),216)
        endif
     endif
  endif

  el = node_value('l ')
  nobs = node_value('obs_point ')
  
  n_align = node_al_errors(al_errors)
  if (n_align.ne.0)  then
     call dcopy(orbit,orbit2,6)
     call tmali1(orbit2,al_errors,betas,gammas,orbit,re)
     call m66mpy(re,rt,rt)
  endif
  
  !---- Element matrix and length.
  call tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te)
  
  !--- if element has a map
  if (fmap) then 
     !---- call m66mpy(re,rt,rt)
     call tmcat(.true.,re,te,rt,tt,rt,tt)
     suml = suml + el
  endif

  if (n_align.ne.0)  then
     call dcopy(orbit,orbit2,6)
     call tmali2(el,orbit2,al_errors,betas,gammas,orbit,re)
     call m66mpy(re,rt,rt)
  endif
  
  if (kobs.gt.0.and.kobs.eq.nobs) return
  
  if (code .ge. ccode-1 .and. code .le. ccode+1)  then !---  kickers (code 14 to 16)
     
     if (thr_on .gt. 0)  then !---  threader is on;  keep matrix
        if (code .le. ccode) then !--- kicker (15) or hkicker (14)
           call dcopy(rt,cmatr(1,1,1),36)
           j = name_len
           call element_name(c_name(1),j)
           coc_cnt(1) = node_value('occ_cnt ')
        endif
        if (code .ge. ccode) then !--- kicker (15) or vkicker (16)
           call dcopy(rt,cmatr(1,1,2),36)
           j = name_len
           call element_name(c_name(2),j)
           coc_cnt(2) = node_value('occ_cnt ')
        endif
     endif
     
  elseif (code .ge. pcode-1 .and. code .le. pcode+1)  then !---  monitors (code 18 to 19)
     
     enable = node_value('enable ')
     if (save .gt. 0 .and. enable .gt. 0) then
        j = 6
        call store_node_vector('orbit_ref ', j, orbit)
     endif

     if (thr_on .gt. 0 .and. enable .gt. 0) then 
        !--- threader is on; test for overflow w.r.t. stored orbit
        call dzero(reforb,6)
        j = 6
        call get_node_vector('orbit_ref ', j, reforb)
        
        do kpro = 1, 2 ! projection plane: x or y
           if (      kpro .eq. 1 .and. code .le. pcode &        !--- plane is x and monitor or hmonitor 
                .or. kpro .eq. 2 .and. code .ge. pcode)  then   !--- plane is y and monitor or vmonitor
              
              j = 2*kpro-1
              dorb(j) = orbit(j) - reforb(j) !--- orbit distortion
              
              if (abs(dorb(j)) .gt. vector(kpro) .and. node .ge. corr_pick(kpro))  then

                 if (debug .ne. 0) then
                    print *,' node = ',node,' kpro = ',kpro,' code = ',code,' dorb = ',dorb(j)
                 endif                    

                 !---  reset count if new pickup
                 if (node .gt. corr_pick(kpro)) rep_cnt(kpro) = 0

                 !---  check for max. repetition
                 rep_cnt(kpro) = rep_cnt(kpro) + 1                 
                 if (rep_cnt(kpro) .gt. max_rep)  then
                    write(tmptxt1, '(i6)') max_rep
                    call fort_warn('threader: pickup skipped after',      &
                         tmptxt1(:6)//' correction attempts')
                    goto 20
                 endif
                 
                 !---  keep matrix
                 call dcopy(rt,pmatr,36)
                 old = node
                 j = name_len
                 call element_name(p_name,j)
                 poc_cnt = node_value('occ_cnt ')
                 call tmthrd(kpro,dorb,cmatr(1,1,kpro),pmatr,vector,node,cick,err)
                 
                 if (err .eq. 1)  then
                    write(tmptxt1, '(i6)') old
                    call fort_warn( 'threader: no corrector before pickup at node ', tmptxt1)
                    ! no return ? 
                 elseif (err .eq. 2)  then
                    call fort_warn('threader: kicker is at start', 'of sequence')
                    ! no return ? 
                 elseif (err .eq. 0)  then
                    corr_pick(kpro) = old
                    write(tmptxt1, '(1p,6g14.6)') cick
                    tmptxt2 = c_name(kpro)
                    write(tmptxt2(lastnb(c_name(kpro))+1:),               &
                         '(''['',i2,'']'')') coc_cnt(kpro)
                    tmptxt3 = p_name
                    write(tmptxt3(lastnb(p_name)+1:),                     &
                         '(''['',i2,'']'')') poc_cnt
                    call fort_info(                                       &
                         '-threader- pickup: ' //tmptxt3(:lastnb(tmptxt3))                 &
                         // '  kicker: '//tmptxt2(:lastnb(tmptxt2))//' total '             &
                         //ptxt(kpro)//'kick:', tmptxt1)
                    !---  restore restart values
                    suml = restsum(kpro)
                    call dcopy(restorb(1,kpro),orbit,6)
                    call dcopy(restm(1,1,kpro),rt,36)
                    call dcopy(restt(1,1,1,kpro),tt,216)
                    goto 20 ! bad! this breaks the loop and kpro=2 might not be treated. 
                 endif

              endif
           endif
        enddo
     endif
  endif
  
20 continue
  
  !---- Test for overflow.
  do j = 1, 6
     if (abs(orbit(j)).ge.orb_limit) then
        eflag = 1
        return
     endif
  enddo

  if (advance_node().ne.0) then
     node=node+1
     goto 10 ! loop over nodes
  endif
  
  bbd_flag=0
end SUBROUTINE tmfrst

SUBROUTINE tmthrd(kpro,dorb,cmatr,pmatr,thrvec,node,cick,error)

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
  integer error,node,kpro
  integer i,j,npick,advance_node,ccode
  integer both,code,itp,lc,lc1,lc2,retreat_node,ncorr
  double precision thrvec(3),dorb(6),cmatr(6,6),pmatr(6,6)
  double precision cick,tol1min,node_value,atemp(6,6)
  parameter (tol1min=1.d-2 , ccode=15)

  error = 0
  
  !---  keep position of current pickup
  npick = node

  both = ccode            ! 15
  itp = both + (-1)**kpro ! 14 or 16
 
  lc = 2 * (kpro - 1)
  lc1 = lc + 1
  lc2 = lc + 2

  !---  look for one preceding corrector of correct (MAD-8) type
10 continue

  if (retreat_node() .eq. 0)  then !--- we are at start already
     do i = node+1,npick
        j = advance_node()
     enddo
     node = npick
     error = 1
     return
  endif

  node = node - 1
  code = node_value('mad8_type ')

  if(code.eq.39) code=15
  if(code.eq.38) code=24

  !-- if wrong corrector or not a corrector, loop
  if (code .ne. itp .and. code .ne. both) goto 10 

  !--- corrector found for the projection plane
  ncorr = node

  !---  transport matrix from kicker to pickup
  call m66inv(cmatr,atemp)
  call m66mpy(pmatr,atemp,atemp)

  !--- if kicker is not efficient enough - loop
  if (abs(atemp(lc1,lc2)) .lt. tol1min) goto 10 

  !---  now we got one good corrector - get kick with attenuation factor
  cick = -thrvec(3) * dorb(2*kpro-1) / atemp(lc1,lc2)

  !---  add kick to kicker strengths
  if (kpro .eq. 1)  then
     cick = node_value('other_bv ') * cick + node_value('chkick ')
     call store_node_value('chkick ', cick)
  else
     cick = node_value('other_bv ') * cick + node_value('cvkick ')
     call store_node_value('cvkick ', cick)
  endif

  !---  set node for restart in front of kicker
  if (retreat_node() .eq. 0)  then !--- we are at start already
     do i = node+1,npick
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
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values for linear coupling parameters.                   *
  !     Input:                                                           *
  !     rt(6,6)      (double)  one turn transfer matrix.                 *
  !     Output:                                                          *
  !     disp0        (double)  initial dispersion vector                 *
  !     r0mat(2,2)   (double)  coupling matrix                           *
  !     eflag        (integer) error flag.                               *
  !----------------------------------------------------------------------*
  integer eflag,stabx,staby,get_option
  double precision rt(6,6),disp0(6),r0mat(2,2),a(2,2),arg,aux(2,2), &
       d(2,2),den,det,dtr,sinmu2,betx0,alfx0,amux0,bety0,alfy0,amuy0,    &
       deltap,get_value,eps,zero,one,two,fourth
  parameter(eps=1d-8,zero=0d0,one=1d0,two=2d0,fourth=0.25d0)
  character(120) msg

  !---- Initialization
  betx0=zero
  bety0=zero
  alfx0=zero
  alfy0=zero
  deltap = get_value('probe ','deltap ')

  !---- Initial dispersion.
  if (get_option('twiss_inval ') .eq. 0) then
     call twdisp(rt,rt(1,6),disp0)
  else
     call dcopy(opt_fun0(15), disp0, 4)
  endif
  disp0(5) = zero
  disp0(6) = one

  !---- Matrix C + B(bar) and its determinant.
  aux(1,1) = rt(3,1) + rt(2,4)
  aux(1,2) = rt(3,2) - rt(1,4)
  aux(2,1) = rt(4,1) - rt(2,3)
  aux(2,2) = rt(4,2) + rt(1,3)
  det = aux(1,1) * aux(2,2) - aux(1,2) * aux(2,1)

  !---- Coupling matrix.
  dtr = (rt(1,1) + rt(2,2) - rt(3,3) - rt(4,4)) / two
  arg = det + dtr**2
  if (arg.ge.zero) then
     if (arg .eq. zero) then
        r0mat(1,1) = one
        r0mat(2,2) = one
        r0mat(1,2) = zero
        r0mat(2,1) = zero
     else
        den = - (dtr + sign(sqrt(arg),dtr))
        r0mat(1,1) = aux(1,1) / den
        r0mat(1,2) = aux(1,2) / den
        r0mat(2,1) = aux(2,1) / den
        r0mat(2,2) = aux(2,2) / den
     endif

     !---- Decouple: Find diagonal blocks.
     a(1,1) = rt(1,1) - rt(1,3)*r0mat(1,1) - rt(1,4)*r0mat(2,1)
     a(1,2) = rt(1,2) - rt(1,3)*r0mat(1,2) - rt(1,4)*r0mat(2,2)
     a(2,1) = rt(2,1) - rt(2,3)*r0mat(1,1) - rt(2,4)*r0mat(2,1)
     a(2,2) = rt(2,2) - rt(2,3)*r0mat(1,2) - rt(2,4)*r0mat(2,2)
     d(1,1) = rt(3,3) + r0mat(1,1)*rt(1,3) + r0mat(1,2)*rt(2,3)
     d(1,2) = rt(3,4) + r0mat(1,1)*rt(1,4) + r0mat(1,2)*rt(2,4)
     d(2,1) = rt(4,3) + r0mat(2,1)*rt(1,3) + r0mat(2,2)*rt(2,3)
     d(2,2) = rt(4,4) + r0mat(2,1)*rt(1,4) + r0mat(2,2)*rt(2,4)

     !---- First mode.
     cosmux = (a(1,1) + a(2,2)) / two
     stabx=0
     if(abs(cosmux).lt.one) stabx=1
     if (stabx.ne.0) then
        sinmu2 = - a(1,2)*a(2,1) - fourth*(a(1,1) - a(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmux = sign(sqrt(sinmu2), a(1,2))
        betx0 = a(1,2) / sinmux
        alfx0 = (a(1,1) - a(2,2)) / (two * sinmux)
     else
        betx0 = zero
        alfx0 = zero
     endif

     !---- Second mode.
     cosmuy = (d(1,1) + d(2,2)) / two
     staby=0
     if(abs(cosmuy).lt.one) staby=1
     if (staby.ne.0) then
        sinmu2 = - d(1,2)*d(2,1) - fourth*(d(1,1) - d(2,2))**2
        if (sinmu2.lt.zero) sinmu2 = eps
        sinmuy = sign(sqrt(sinmu2), d(1,2))
        bety0 = d(1,2) / sinmuy
        alfy0 = (d(1,1) - d(2,2)) / (two * sinmuy)
     else
        bety0 = zero
        alfy0 = zero
     endif

     !---- Unstable due to coupling.
  else
     stabx = 0
     staby = 0
  endif

  !---- Initial phase angles.
  amux0 = zero
  amuy0 = zero

  !---- Give message, if unstable.
  eflag = 0
  if (stabx+staby.lt.2) then
     eflag = 1
     if (staby.ne.0) then
        write (msg, 910) 1,deltap,cosmux,cosmuy
        call aawarn('TWCPIN: ',msg)
     else if (stabx.ne.0) then
        write (msg, 910) 2,deltap,cosmux,cosmuy
        call aawarn('TWCPIN: ',msg)
     else
        write (msg, 920) deltap,cosmux,cosmuy
        call aawarn('TWCPIN: ',msg)
     endif
  endif
  opt_fun0(3)=betx0
  opt_fun0(4)=alfx0
  opt_fun0(5)=amux0
  opt_fun0(6)=bety0
  opt_fun0(7)=alfy0
  opt_fun0(8)=amuy0
  opt_fun0(15)=disp0(1)
  opt_fun0(16)=disp0(2)
  opt_fun0(17)=disp0(3)
  opt_fun0(18)=disp0(4)
  opt_fun0(29)=r0mat(1,1)
  opt_fun0(30)=r0mat(1,2)
  opt_fun0(31)=r0mat(2,1)
  opt_fun0(32)=r0mat(2,2)

910 format('Mode ',i1,' is unstable for delta(p)/p =',f12.6,          &
       ', cosmux = ',f12.6,', cosmuy = ',f12.6)
920 format('Both modes are unstable for delta(p)/p = ',f12.6,         &
       ', cosmux = ',f12.6,', cosmuy = ',f12.6)
end SUBROUTINE twcpin
SUBROUTINE twdisp(rt,vect,disp)

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
  integer i,j,irank
  double precision rt(6,6),vect(6),disp(6),a(4,5),zero,one
  parameter(zero=0d0,one=1d0)

  do i = 1, 4
     do j = 1, 4
        a(i,j) = rt(i,j)
     enddo
     a(i,i) = a(i,i) - one
     a(i,5) = - vect(i)
  enddo

  call solver(a,4,1,irank)

  if (irank.ge.4) then
     do i = 1, 4
        disp(i) = a(i,5)
     enddo
  else
     call aawarn('TWDISP: ',                                         &
          'Unable to compute dispersion --- dispersion set to zero.')
     call dzero(disp,4)
  endif

end SUBROUTINE twdisp
SUBROUTINE twcpgo(rt,orbit0)

  use twiss0fi
  use twisslfi
  use twisscfi
  use twiss_elpfi
  use twissotmfi
  use spch_bbfi

  use name_lenfi
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
  logical fmap,cplxy,cplxt,dorad,sector_sel,mycentre_cptk
  integer i,iecnt,code,save,advance_node,restart_sequ,get_option,   &
       node_al_errors,n_align
  double precision rt(6,6),ek(6),re(6,6),rwi(6,6),rc(6,6),  &
       te(6,6,6),el,orbit(6),orbit2(6),betas,gammas,                 &
       al_errors(align_max),bvk,sumloc,pos0,node_value,get_value,sd,zero,&
       one,two
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0)
  character(130) msg

  !--- 2013-Nov-14  10:36:35  ghislain: adding name of element for output of name
  character*(name_len) el_name, bxmax_name, bymax_name, dxmax_name, dymax_name, xcomax_name, ycomax_name

  double precision orbit0(6)
  
  !---- Initialization
  pos0=zero
  amux=zero
  amuy=zero
  currpos=zero
  sumloc = zero
  dorad = get_value('probe ','radiate ').ne.zero
  centre = get_value('twiss ','centre ').ne.zero
  call m66one(rwi)
  call m66one(rc)

  !---- Store one-turn-map
  call dcopy(rt,rotm,36)

  !---- Create internal table for lattice functions if requested
  save=get_option('twiss_save ')
  call dcopy(opt_fun0,opt_fun,fundim)

  !---- Initial values for lattice functions.
  betx    =opt_fun(3 )
  alfx    =opt_fun(4 )
  amux    =opt_fun(5 )
  bety    =opt_fun(6 )
  alfy    =opt_fun(7 )
  amuy    =opt_fun(8 )
  orbit(1)=opt_fun(9 )
  orbit(2)=opt_fun(10)
  orbit(3)=opt_fun(11)
  orbit(4)=opt_fun(12)
  orbit(5)=opt_fun(13)
  orbit(6)=opt_fun(14)
  disp(1) =opt_fun(15)
  disp(2) =opt_fun(16)
  disp(3) =opt_fun(17)
  disp(4) =opt_fun(18)
  disp(5) =zero
  disp(6) =one
  rmat(1,1)=opt_fun(29)
  rmat(1,2)=opt_fun(30)
  rmat(2,1)=opt_fun(31)
  rmat(2,2)=opt_fun(32)

  !--- 2014-May-30  15:19:52  ghislain: if initial values are provided, copy the initial orbit
  !                 because opt_fun contains only the values given on the twiss command itself
  !                 but does not know the values coming from COGUESS or USEORBIT
  if (get_option('twiss_inval ') .ne. 0) orbit=orbit0

  !---- Maximum and r.m.s. values.
  bxmax =betx
  dxmax =disp(1)
  bymax =bety
  dymax =disp(3)
  xcomax=zero
  ycomax=zero
  sigxco=zero
  sigyco=zero
  sigdx =zero
  sigdy =zero

  !--- 2013-Nov-14  10:42:16  ghislain: need to initialize names of elements where maxima occur
  bxmax_name = 'nil '
  bymax_name = 'nil ' 
  dxmax_name = 'nil '
  dymax_name = 'nil ' 
  xcomax_name = 'nil '
  ycomax_name = 'nil '

  !---- Loop over positions.
  cplxy = .false.
  cplxt = .false.
  iecnt=0
  betas = get_value('probe ','beta ')
  gammas= get_value('probe ','gamma ')
  centre_cptk=.false.
  i = restart_sequ()
  i_spch=0

10 continue

  sector_sel = node_value('sel_sector ') .ne. zero .and. sectormap
  code = node_value('mad8_type ')
  if(code.eq.39) code=15
  if(code.eq.38) code=24
  bvk = node_value('other_bv ')
  elpar_vl = el_par_vector(g_polarity, g_elpar)
  el = node_value('l ')

  !--- 2013-Nov-14  10:34:00  ghislain: add acquistion of name of element here.  
  call element_name(el_name,len(el_name))

  opt_fun(70) = g_elpar(g_kmax)
  opt_fun(71) = g_elpar(g_kmin)
  opt_fun(72) = g_elpar(g_calib)
  opt_fun(73) = g_elpar(g_polarity)

  n_align = node_al_errors(al_errors)
  if (n_align.ne.0)  then
     call dcopy(orbit,orbit2,6)
     call tmali1(orbit2,al_errors,betas,gammas,orbit,re)
     mycentre_cptk = centre_cptk
     centre_cptk = .false.
     call twcptk(re,orbit)
     centre_cptk = mycentre_cptk
     if (sectormap) call m66mpy(re,srmat,srmat)
  endif

  if(centre) centre_cptk=.true.

  call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te)

  if(centre) then
     pos0=currpos
     currpos=currpos+el/two
     sd = rt(5,6)
     do i = 1, 4
        sd = sd + rt(5,i) * disp(i)
     enddo
     eta = - sd * betas**2 / circ
     alfa = one / gammas**2 + eta
     opt_fun(74)=alfa
     call twprep(save,1,opt_fun,currpos)
  endif

  centre_cptk=.false.

  if (fmap) then
     call twcptk(re,orbit)
     if (sectormap) call tmcat(.true.,re,te,srmat,stmat,srmat,stmat)
  endif

  if (n_align.ne.0)  then
     call dcopy(orbit,orbit2,6)
     call tmali2(el,orbit2,al_errors,betas,gammas,orbit,re)
     mycentre_cptk = centre_cptk
     centre_cptk = .false.
     call twcptk(re,orbit)
     centre_cptk = mycentre_cptk
     if (sectormap) call m66mpy(re,srmat,srmat)
  endif

  sumloc = sumloc + el
  if (sector_sel) call twwmap(sumloc, orbit)
  sd = rt(5,6)
  do i = 1, 4
     sd = sd + rt(5,i) * disp(i)
  enddo
  eta = - sd * betas**2 / circ
  alfa = one / gammas**2 + eta
  opt_fun(74)=alfa
  if(centre) then
     bxmax =max(opt_fun(3)      ,bxmax)
     bymax =max(opt_fun(6)      ,bymax)
     dxmax =max(abs(opt_fun(15)),dxmax)
     dymax =max(abs(opt_fun(17)),dymax)
     xcomax=max(abs(opt_fun(9 )),xcomax)
     ycomax=max(abs(opt_fun(11)),ycomax)
  else
     opt_fun(9 )=orbit(1)
     opt_fun(10)=orbit(2)
     opt_fun(11)=orbit(3)
     opt_fun(12)=orbit(4)
     opt_fun(13)=orbit(5)
     opt_fun(14)=orbit(6)
     currpos=currpos+el
  endif
  iecnt=iecnt+1

  !--- 2013-Nov-14  10:41:42  ghislain: change to save name of element where max is happening
  ! bxmax =max(betx         ,bxmax)
  ! bymax =max(bety         ,bymax)
  ! dxmax =max(abs(disp(1)) ,dxmax)
  ! dymax =max(abs(disp(3)) ,dymax)
  ! xcomax=max(abs(orbit(1)),xcomax)
  ! ycomax=max(abs(orbit(3)),ycomax)
  !
  if (betx .gt. bxmax) then
     bxmax = betx
     bxmax_name = el_name
  endif
  if (bety .gt. bymax) then
     bymax = bety
     bymax_name = el_name
  endif
  if (abs(disp(1)) .gt. dxmax) then
     dxmax = abs(disp(1))
     dxmax_name = el_name
  endif
  if (abs(disp(3)) .gt. dymax) then
     dymax = abs(disp(3))
     dymax_name = el_name
  endif
  if (abs(orbit(1)) .gt. xcomax) then
     xcomax = abs(orbit(1))
     xcomax_name = el_name
  endif
  if (abs(orbit(3)) .gt. ycomax) then
     ycomax = abs(orbit(3))
     ycomax_name = el_name
  endif
  !--- 2013-Nov-14  10:41:42  ghislain: end of change

  !--- 2013-Nov-14  14:18:07  ghislain: should only calculate the RMS values for active elements, 
  !                           skipping markers(25), beambeam(22), changeref(35), translation(36),
  !                           srotation(12), yrotation(13) etc...
  !                           But leave drifts(1) in for now. They are the most numerous elements in machines
  ! if ( & ! code .ne. 1 .and. & 
  !      code .ne. 22 .and. & 
  !      code .ne. 25 .and. & 
  !      code .ne. 35 .and. code .ne. 36 .and. &
  !      code .ne. 12 .and. code .ne. 13) then
     sigxco=sigxco+orbit(1)**2
     sigyco=sigyco+orbit(3)**2
     sigdx =sigdx + disp(1)**2
     sigdy =sigdy + disp(3)**2
  !endif

  if(.not.centre) call twprep(save,1,opt_fun,currpos)
  if(centre) then
     currpos=pos0+el
     opt_fun(2 )=currpos
     opt_fun(3 )=betx
     opt_fun(4 )=alfx
     opt_fun(5 )=amux
     opt_fun(6 )=bety
     opt_fun(7 )=alfy
     opt_fun(8 )=amuy
     opt_fun(9 )=orbit(1)
     opt_fun(10)=orbit(2)
     opt_fun(11)=orbit(3)
     opt_fun(12)=orbit(4)
     opt_fun(13)=orbit(5)
     opt_fun(14)=orbit(6)
     opt_fun(15)=disp(1)
     opt_fun(16)=disp(2)
     opt_fun(17)=disp(3)
     opt_fun(18)=disp(4)
     opt_fun(29)=rmat(1,1)
     opt_fun(30)=rmat(1,2)
     opt_fun(31)=rmat(2,1)
     opt_fun(32)=rmat(2,2)
  endif

  if(advance_node().ne.0) goto 10

  !---- Compute summary.
  wgt = max(iecnt, 1)
  sigxco=sqrt(sigxco / wgt)
  sigyco=sqrt(sigyco / wgt)
  sigdx =sqrt(sigdx / wgt)
  sigdy =sqrt(sigdy / wgt)
  cosmux=(rt(1,1) + rt(2,2)) / two
  cosmuy=(rt(3,3) + rt(4,4)) / two

  !---- Warning messages.
  if (cplxt.or.dorad) then
     write (msg, 910)
     call aawarn('TWCPGO: ',msg)
  endif

910 format('TWISS uses the RF system or synchrotron radiation ',      &
       'only to find the closed orbit, for optical calculations it ',    &
       'ignores both.')

  ! write (6,*) bxmax, bxmax_name
  ! write (6,*) bymax, bymax_name
  ! write (6,*) dxmax, dxmax_name
  ! write (6,*) dymax, dymax_name
  ! write (6,*) xcomax, xcomax_name
  ! write (6,*) ycomax, ycomax_name


end SUBROUTINE twcpgo
SUBROUTINE twcptk(re,orbit)
  use twiss0fi
  use twisslfi
  use twisscfi
  use twissotmfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track coupled lattice functions.                                 *
  !     Input:                                                           *
  !     re(6,6)  (double)   transfer matrix of element.                  *
  !     orbit(6) (double)   closed orbit                                 *
  !----------------------------------------------------------------------*
  integer i,i1,i2,j,get_option
  double precision re(6,6),orbit(6),rw0(6,6),rwi(6,6),rc(6,6),      &
       rmat0(2,2),a(2,2),adet,b(2,2),c(2,2),dt(6),tempa,tempb,alfx0,     &
       alfy0,betx0,bety0,amux0,amuy0,zero,one,eps
  parameter(zero=0d0,one=1d0,eps=1d-36)

  !initialize
  bety0=zero
  betx0=zero
  amux0=zero
  amuy0=zero
  alfy0=zero
  alfx0=zero

  !---- Dispersion.
  call dzero(dt,6)
  do i = 1, 6
     do j = 1, 6
        dt(i) = dt(i) + re(i,j) * disp(j)
     enddo
  enddo
  if(.not.centre.or.centre_cptk) then
     opt_fun(15)=dt(1)
     opt_fun(16)=dt(2)
     opt_fun(17)=dt(3)
     opt_fun(18)=dt(4)
  endif
  if(centre_cptk) then
     alfx0=alfx
     alfy0=alfy
     betx0=betx
     bety0=bety
     amux0=amux
     amuy0=amuy
     call dcopy(rmat,rmat0,4)
     if(rmatrix) call dcopy(rw,rw0,36)
  else
     call dcopy(dt,disp,6)
     disp(5) = zero
     disp(6) = one
  endif

  !---- Auxiliary matrices.
  a(1,1) = re(1,1) - (re(1,3) * rmat(1,1) + re(1,4) * rmat(2,1))
  a(1,2) = re(1,2) - (re(1,3) * rmat(1,2) + re(1,4) * rmat(2,2))
  a(2,1) = re(2,1) - (re(2,3) * rmat(1,1) + re(2,4) * rmat(2,1))
  a(2,2) = re(2,2) - (re(2,3) * rmat(1,2) + re(2,4) * rmat(2,2))
  b(1,1) = re(3,1) - (re(3,3) * rmat(1,1) + re(3,4) * rmat(2,1))
  b(1,2) = re(3,2) - (re(3,3) * rmat(1,2) + re(3,4) * rmat(2,2))
  b(2,1) = re(4,1) - (re(4,3) * rmat(1,1) + re(4,4) * rmat(2,1))
  b(2,2) = re(4,2) - (re(4,3) * rmat(1,2) + re(4,4) * rmat(2,2))
  c(1,1) = re(3,3) + (re(3,1) * rmat(2,2) - re(3,2) * rmat(2,1))
  c(1,2) = re(3,4) - (re(3,1) * rmat(1,2) - re(3,2) * rmat(1,1))
  c(2,1) = re(4,3) + (re(4,1) * rmat(2,2) - re(4,2) * rmat(2,1))
  c(2,2) = re(4,4) - (re(4,1) * rmat(1,2) - re(4,2) * rmat(1,1))

  !---- Track R matrix.
  adet = a(1,1) * a(2,2) - a(1,2) * a(2,1)
  if (abs(adet).gt.eps) then
     rmat(1,1) = - (b(1,1) * a(2,2) - b(1,2) * a(2,1)) / adet
     rmat(1,2) =   (b(1,1) * a(1,2) - b(1,2) * a(1,1)) / adet
     rmat(2,1) = - (b(2,1) * a(2,2) - b(2,2) * a(2,1)) / adet
     rmat(2,2) =   (b(2,1) * a(1,2) - b(2,2) * a(1,1)) / adet
    
     !---- Mode 1.
     tempb = a(1,1) * betx - a(1,2) * alfx
     tempa = a(2,1) * betx - a(2,2) * alfx
     alfx = - (tempa * tempb + a(1,2) * a(2,2)) / (adet * betx)
     betx =   (tempb * tempb + a(1,2) * a(1,2)) / (adet * betx)
     if(abs(a(1,2)).gt.eps) amux=amux+atan2(a(1,2),tempb)
     
     !---- Mode 2.
     tempb = c(1,1) * bety - c(1,2) * alfy
     tempa = c(2,1) * bety - c(2,2) * alfy
     alfy = - (tempa * tempb + c(1,2) * c(2,2)) / (adet * bety)
     bety =   (tempb * tempb + c(1,2) * c(1,2)) / (adet * bety)
     if(abs(c(1,2)).gt.eps) amuy=amuy+atan2(c(1,2),tempb)
  endif
  
  !---- Cummulative R matrix and one-turn map at element location.
  if(rmatrix) then
     call m66mpy(re,rw,rw)
     if (get_option('twiss_inval ') .ne. 0) then
        call dcopy(rw,rc,36)
     else
        call m66inv(rw,rwi)
        call m66mpy(rotm,rwi,rc)
        call m66mpy(rw,rc,rc)
     endif
  endif
  
  if(.not.centre.or.centre_cptk) then
     opt_fun(3 )=betx
     opt_fun(4 )=alfx
     opt_fun(5 )=amux
     opt_fun(6 )=bety
     opt_fun(7 )=alfy
     opt_fun(8 )=amuy
     opt_fun(29)=rmat(1,1)
     opt_fun(30)=rmat(1,2)
     opt_fun(31)=rmat(2,1)
     opt_fun(32)=rmat(2,2)
  endif
  if(rmatrix) then
     do i1=1,6
        do i2=1,6
           opt_fun(33+(i1-1)*6+i2)=rc(i1,i2)
        enddo
     enddo
  endif
  if(centre_cptk) then
     opt_fun(9 )=orbit(1)
     opt_fun(10)=orbit(2)
     opt_fun(11)=orbit(3)
     opt_fun(12)=orbit(4)
     opt_fun(13)=orbit(5)
     opt_fun(14)=orbit(6)
     alfx=alfx0
     alfy=alfy0
     betx=betx0
     bety=bety0
     amux=amux0
     amuy=amuy0
     call dcopy(rmat0,rmat,4)
     if(rmatrix) call dcopy(rw0,rw,36)
  endif

end SUBROUTINE twcptk
SUBROUTINE twbtin(rt,tt)

  use twiss0fi
  use twisscfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Initial values of Chromatic functions.                           *
  !     Input:                                                           *
  !     rt(6,6)   (double)  transfer matrix.                             *
  !     tt(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  logical stabx,staby
  integer i,k,j,get_option
  double precision rt(6,6),tt(6,6,6),disp0(6),ddisp0(6),rtp(6,6),   &
       sinmu2,bx,ax,by,ay,eps,temp,aux(6),zero,one,two,fourth,           &
       get_variable,twopi
  parameter(eps=1d-8,zero=0d0,one=1d0,two=2d0,fourth=0.25d0)

  !---- Initialization
  twopi=get_variable('twopi ')
  betx   =opt_fun0(3 )
  alfx   =opt_fun0(4 )
  amux   =opt_fun0(5 )
  bety   =opt_fun0(6 )
  alfy   =opt_fun0(7 )
  amuy   =opt_fun0(8 )
  wx     =opt_fun0(19)
  phix   =opt_fun0(20)
  dmux   =opt_fun0(21)
  wy     =opt_fun0(22)
  phiy   =opt_fun0(23)
  dmuy   =opt_fun0(24)

  !---- Initial value flag.
  if (get_option('twiss_inval ') .ne. 0) then
     call dcopy(opt_fun0(15),disp,4)
     call dcopy(opt_fun0(25),ddisp,4)
     disp(5) = zero
     disp(6) = one
     ddisp(5) = zero
     ddisp(6) = zero
     return
  endif

  !---- Initial dispersion.
  call twdisp(rt,rt(1,6),disp0)
  disp0(5) = zero
  disp0(6) = one

  !---- Derivative of transfer matrix w.r.t. delta(p)/p.
  call dzero(aux,6)
  do i = 1, 6
     do k = 1, 6
        temp = zero
        do j = 1, 6
           temp = temp + tt(i,j,k) * disp0(j)
        enddo
        aux(i) = aux(i) + temp * disp0(k)
        rtp(i,k) = two * temp
     enddo
  enddo

  !---- Derivative of dispersion.
  call twdisp(rt,aux,ddisp0)

  ddisp0(5) = zero
  ddisp0(6) = zero
  call dcopy(disp0,disp,6)
  call dcopy(ddisp0,ddisp,6)

  !---- Horizontal motion.
  cosmux = (rt(1,1) + rt(2,2)) / two
  stabx = abs(cosmux) .lt. one
  if (stabx) then
     sinmu2 = - rt(1,2)*rt(2,1) - fourth*(rt(1,1) - rt(2,2))**2
     if (sinmu2.lt.0) sinmu2 = eps
     sinmux = sign(sqrt(sinmu2), rt(1,2))
     betx = rt(1,2) / sinmux
     alfx = (rt(1,1) - rt(2,2)) / (two * sinmux)
     bx = rtp(1,2) / rt(1,2) +                                       &
          (rtp(1,1) + rtp(2,2)) * cosmux / (two * sinmu2)
     ax = (rtp(1,1) - rtp(2,2)) / (two * sinmux) -                   &
          alfx * rtp(1,2) / rt(1,2)
     wx = sqrt(bx**2 + ax**2)
     if (wx.gt.eps) phix = atan2(ax,bx)
  endif

  !---- Vertical motion.
  cosmuy = (rt(3,3) + rt(4,4)) / two
  staby = abs(cosmuy) .lt. one
  if (staby) then
     sinmu2 = - rt(3,4)*rt(4,3) - fourth*(rt(3,3) - rt(4,4))**2
     if (sinmu2.lt.0) sinmu2 = eps
     sinmuy = sign(sqrt(sinmu2), rt(3,4))
     bety = rt(3,4) / sinmuy
     alfy = (rt(3,3) - rt(4,4)) / (two * sinmuy)
     by = rtp(3,4) / rt(3,4) +                                       &
          (rtp(3,3) + rtp(4,4)) * cosmuy / (two * sinmu2)
     ay = (rtp(3,3) - rtp(4,4)) / (two * sinmuy) -                   &
          alfy * rtp(3,4) / rt(3,4)
     wy = sqrt(by**2 + ay**2)
     if (wy.gt.eps) phiy = atan2(ay,by)
  endif

  !---- Fill optics function array
  opt_fun0(19)=wx
  opt_fun0(20)=phix
  opt_fun0(22)=wy
  opt_fun0(23)=phiy
  opt_fun0(25)=ddisp(1)
  opt_fun0(26)=ddisp(2)
  opt_fun0(27)=ddisp(3)
  opt_fun0(28)=ddisp(4)

end SUBROUTINE twbtin
SUBROUTINE twchgo

  use twiss0fi
  use twisslfi
  use twissafi
  use twisscfi
  use spch_bbfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track Chromatic functions.                                       *
  !----------------------------------------------------------------------*
  logical fmap,cplxy,cplxt,dorad,mycentre_bttk
  integer i,code,save,restart_sequ,advance_node,get_option,n_align, &
       node_al_errors
  double precision orbit(6),orbit2(6),ek(6),re(6,6),        &
       te(6,6,6),al_errors(align_max),deltap,el,betas,gammas,  &
       node_value,get_value,pos0,zero,one,two
  parameter(zero=0d0,one=1d0,two=2d0)
  character(130) msg

  !---- If save requested reset table
  save=get_option('twiss_save ')
  if(save.ne.0) call reset_count(table_name)
  deltap = get_value('probe ','deltap ')
  dorad = get_value('probe ','radiate ').ne.zero
  centre = get_value('twiss ','centre ').ne.zero

  !---- Initial values for lattice functions.
  pos0=zero
  amux=zero
  amuy=zero
  orbit(1) =opt_fun0(9 )
  orbit(2) =opt_fun0(10)
  orbit(3) =opt_fun0(11)
  orbit(4) =opt_fun0(12)
  orbit(5) =opt_fun0(13)
  orbit(6) =opt_fun0(14)
  disp(1)  =opt_fun0(15)
  disp(2)  =opt_fun0(16)
  disp(3)  =opt_fun0(17)
  disp(4)  =opt_fun0(18)
  disp(5)  =zero
  disp(6)  =one
  call dzero(te,216)

  !---- Initial values for chromatic functions.
  opt_fun(19)=wx
  opt_fun(20)=phix
  opt_fun(21)=dmux
  opt_fun(22)=wy
  opt_fun(23)=phiy
  opt_fun(24)=dmuy
  opt_fun(25)=ddisp(1)
  opt_fun(26)=ddisp(2)
  opt_fun(27)=ddisp(3)
  opt_fun(28)=ddisp(4)
  synch_1 = zero
  synch_2 = zero
  synch_3 = zero
  synch_4 = zero
  synch_5 = zero

  !---- Loop over positions.
  cplxy = .false.
  cplxt = .false.

  betas = get_value('probe ','beta ')
  gammas= get_value('probe ','gamma ')
  centre_bttk=.false.
  i = restart_sequ()
  if(centre) currpos = zero
  i_spch=0

10 continue
  el = node_value('l ')
  code = node_value('mad8_type ')
  if(code.eq.39) code=15
  if(code.eq.38) code=24

  !---- Physical element.
  n_align = node_al_errors(al_errors)
  if (n_align.ne.0)  then
     call dcopy(orbit,orbit2,6)
     call tmali1(orbit2,al_errors,betas,gammas,orbit,re)
     mycentre_bttk = centre_bttk
     centre_bttk = .false.
     call twbttk(re,te)
     centre_bttk = mycentre_bttk
  endif
  if(centre) centre_bttk=.true.
  call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te)
  if(centre) then
     pos0=currpos
     currpos=currpos+el/two
     call twprep(save,2,opt_fun,currpos)
  endif
  centre_bttk=.false.
  if (fmap) then
     call twbttk(re,te)
  endif
  if (n_align.ne.0)  then
     call dcopy(orbit,orbit2,6)
     call tmali2(el,orbit2,al_errors,betas,gammas,orbit,re)
     mycentre_bttk = centre_bttk
     centre_bttk = .false.
     call twbttk(re,te)
     centre_bttk = mycentre_bttk
  endif
  if(.not.centre) call twprep(save,2,opt_fun,zero)
  if(centre) then
     currpos=pos0+el
     opt_fun(2 )=currpos
     opt_fun(3 )=betx
     opt_fun(4 )=alfx
     opt_fun(5 )=amux
     opt_fun(6 )=bety
     opt_fun(7 )=alfy
     opt_fun(8 )=amuy
     opt_fun(9 )=orbit(1)
     opt_fun(10)=orbit(2)
     opt_fun(11)=orbit(3)
     opt_fun(12)=orbit(4)
     opt_fun(13)=orbit(5)
     opt_fun(14)=orbit(6)
     opt_fun(15)=disp(1)
     opt_fun(16)=disp(2)
     opt_fun(17)=disp(3)
     opt_fun(18)=disp(4)
     opt_fun(19)=wx
     opt_fun(20)=phix
     opt_fun(21)=dmux
     opt_fun(22)=wy
     opt_fun(23)=phiy
     opt_fun(24)=dmuy
     opt_fun(25)=ddisp(1)
     opt_fun(26)=ddisp(2)
     opt_fun(27)=ddisp(3)
     opt_fun(28)=ddisp(4)
     opt_fun(29)=rmat(1,1)
     opt_fun(30)=rmat(1,2)
     opt_fun(31)=rmat(2,1)
     opt_fun(32)=rmat(2,2)
  endif
  if (advance_node().ne.0)  goto 10

  !---- Warning, if system is coupled.
  if (cplxy) then
     write (msg, 910) deltap
     call aawarn('TWCHGO: ',msg)
  endif
  if (cplxt.or.dorad) then
     write (msg, 920)
     call aawarn('TWCHGO: ',msg)
  endif

910 format('TWISS found transverse coupling for delta(p)/p =',f12.6,  &
       'chromatic functions may be wrong.')
920 format('TWISS uses the RF system or synchrotron radiation ',      &
       'only to find the closed orbit, for optical calculations it ',    &
       'ignores both.')

end SUBROUTINE twchgo
SUBROUTINE twbttk(re,te)

  use twiss0fi
  use twisslfi
  use twisscfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Track lattice functions, including chromatic effects.            *
  !     Input:                                                           *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  integer save,i,j,k,get_option
  double precision re(6,6),te(6,6,6),aux(6),auxp(6),ax1,ax2,ay1,ay2,&
       bx1,bx2,by1,by2,proxim,rep(6,6),t2,ta,tb,temp,tg,fre(6,6),        &
       frep(6,6),curlyh,detl,f,rhoinv,blen,alfx0,alfy0,betx0,bety0,amux0,&
       amuy0,wx0,wy0,dmux0,dmuy0,rmat0(2,2),phix0,phiy0,node_value,eps,  &
       zero,one,two
  parameter(eps=1d-8,zero=0d0,one=1d0,two=2d0)

  double precision an, e1, e2, sk1
  double precision syncint(5), beta, get_value
  
  !initialize
  wx0=zero
  wy0=zero
  phiy0=zero
  phix0=zero
  dmuy0=zero
  dmux0=zero
  betx0=zero
  bety0=zero
  amuy0=zero
  amux0=zero
  alfx0=zero
  alfy0=zero
  !---- Create internal table for lattice functions if requested
  save=get_option('twiss_save ')

  !---- Tor: needed for synchrotron integrals
  blen=node_value('blen ')
  rhoinv=node_value('rhoinv ')
  ! 2014-May-15  16:50:36  ghislain: added 
  sk1 = node_value('k1 ')
  e1 = node_value('e1 ')
  e2 = node_value('e2 ')
  
  an = node_value('angle ')
  if(node_value('mad8_type ').eq.2) then
     e1 = e1 + an / two
     e2 = e2 + an / two
  endif

  !---- Synchrotron radiation integrals through bending magnets.    
  if(.not.centre_bttk .and. rhoinv .ne. 0.d0) then
     ! Note that calcsyncint expects dx and dpx as derivatives wrt deltap.
     ! since MAD take disp(1) and disp(2) as derivatives wrt pt, they must be 
     ! multiplied by beta before the call to calcsyncint.
     beta = get_value('probe ','beta ')
     call calcsyncint(rhoinv,blen,sk1,e1,e2,betx,alfx,disp(1)*beta,disp(2)*beta,syncint)
     synch_1 = synch_1 + syncint(1)
     synch_2 = synch_2 + syncint(2)
     synch_3 = synch_3 + syncint(3)
     synch_4 = synch_4 + syncint(4)
     synch_5 = synch_5 + syncint(5)
  endif

  call dzero(aux,6)
  call dzero(auxp,6)
  do i = 1, 6
     do k = 1, 6
        temp = zero
        do j = 1, 6
           temp = temp + te(i,j,k)*disp(j)
        enddo
        aux(i) = aux(i) + re(i,k)*disp(k)
        auxp(i) = auxp(i) + temp*disp(k) + re(i,k)*ddisp(k)
        rep(i,k) = two*temp
     enddo
  enddo
  if(.not.centre_bttk) then
     call dcopy(aux,disp,6)
     call dcopy(auxp,ddisp,6)
  else
     alfx0=alfx
     alfy0=alfy
     betx0=betx
     bety0=bety
     amux0=amux
     amuy0=amuy
     wx0=wx
     wy0=wy
     dmux0=dmux
     dmuy0=dmuy
     phix0=phix
     phiy0=phiy
     do i=1,2
        do j=1,2
           rmat0(i,j)=rmat(i,j)
        enddo
     enddo
  endif

  !---- Tor: modified to cancel energy change
  disp(6) = one

  !---- Tor/MDW: scale by square root of the determinant of the
  !     longitudinal 2x2 part of the R-matrix
  detl = re(5,5)*re(6,6) - re(5,6)*re(6,5)
  f = one / sqrt(detl)
  call m66scl(f,re,fre)
  call m66scl(f,rep,frep)

  !---- Track horizontal functions including energy scaling.
  tb = fre(1,1)*betx - fre(1,2)*alfx
  ta = fre(2,1)*betx - fre(2,2)*alfx
  t2 = tb**2 + fre(1,2)**2
  tg = fre(1,1)*alfx - fre(1,2)*(one + alfx**2) / betx

  !---- Linear functions.
  alfx = - (tb*ta + fre(1,2)*fre(2,2)) / betx
  betx = t2 / betx
  if(fre(1,2).ne.zero.or.tb.ne.zero) amux=amux+atan2(fre(1,2),tb)
  bx1 = wx*cos(phix)
  ax1 = wx*sin(phix)
  bx2 = ((tb**2 - fre(1,2)**2)*bx1                                  &
       - two*tb*fre(1,2)*ax1) / t2                                       &
       + two*(tb*frep(1,1) - tg*frep(1,2)) / betx
  ax2 = ((tb**2 - fre(1,2)**2)*ax1                                  &
       + two*tb*fre(1,2)*bx1) / t2                                       &
       - (tb*(frep(1,1)*alfx + frep(2,1)*betx)                           &
       - tg*(frep(1,2)*alfx + frep(2,2)*betx)                            &
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
  if(fre(3,4).ne.zero.or.tb.ne.zero) amuy=amuy+atan2(fre(3,4),tb)

  by1 = wy*cos(phiy)
  ay1 = wy*sin(phiy)
  by2 = ((tb**2 - fre(3,4)**2)*by1                                  &
       - two*tb*fre(3,4)*ay1) / t2                                       &
       + two*(tb*frep(3,3) - tg*frep(3,4)) / bety
  ay2 = ((tb**2 - fre(3,4)**2)*ay1                                  &
       + two*tb*fre(3,4)*by1) / t2                                       &
       - (tb*(frep(3,3)*alfy + frep(4,3)*bety)                           &
       - tg*(frep(3,4)*alfy + frep(4,4)*bety)                            &
       + fre(3,3)*frep(3,4) - fre(3,4)*frep(3,3)) / bety
  wy = sqrt(ay2**2 + by2**2)
  if (wy.gt.eps) phiy = proxim(atan2(ay2, by2), phiy)
  dmuy = dmuy + fre(3,4)*(fre(3,4)*ay1 - tb*by1) / t2               &
       + (fre(3,3)*frep(3,4) - fre(3,4)*frep(3,3)) / bety

  !---- Fill optics function array
  if(.not.centre.or.centre_bttk) then
     opt_fun(19)=wx
     opt_fun(20)=phix
     opt_fun(21)=dmux
     opt_fun(22)=wy
     opt_fun(23)=phiy
     opt_fun(24)=dmuy
     opt_fun(25)=auxp(1)
     opt_fun(26)=auxp(2)
     opt_fun(27)=auxp(3)
     opt_fun(28)=auxp(4)
  endif
  if(centre_bttk) then
     alfx=alfx0
     alfy=alfy0
     betx=betx0
     bety=bety0
     amux=amux0
     amuy=amuy0
     wx=wx0
     wy=wy0
     dmux=dmux0
     dmuy=dmuy0
     phix=phix0
     phiy=phiy0
     do i=1,2
        do j=1,2
           rmat(i,j)=rmat0(i,j)
        enddo
     enddo
  endif

end SUBROUTINE twbttk


SUBROUTINE tw_summ(rt,tt)

  use twiss0fi
  use twisscfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Compute summary data for TWISS and OPTICS commands.              *
  !     Input:                                                           *
  !     rt(6,6)   (double)  one turn transfer matrix.                    *
  !     tt(6,6,6) (double)  second order terms.                          *
  !----------------------------------------------------------------------*
  integer i,get_option
  double precision rt(6,6),tt(6,6,6),deltap,sd,betas,gammas,        &
       disp0(6),detl,f,tb,frt(6,6),frtp(6,6),rtp(6,6),t2,bx0,ax0,by0,ay0,&
       sx,sy,orbit5,twopi,get_variable,get_value,zero,one,two
  parameter(zero=0d0,one=1d0,two=2d0)

  !---- Initialization chromatic part
  twopi=get_variable('twopi ')
  call dzero(rtp,36)
  call dzero(frt,36)
  call dzero(frtp,36)
  deltap = get_value('probe ','deltap ')
  betas = get_value('probe ','beta ')
  gammas= get_value('probe ','gamma ')
  disp0(1)=opt_fun0(15)
  disp0(2)=opt_fun0(16)
  disp0(3)=opt_fun0(17)
  disp0(4)=opt_fun0(18)
  wx      =opt_fun(19)
  phix    =opt_fun(20)
  dmux    =opt_fun(21)
  wy      =opt_fun(22)
  phiy    =opt_fun(23)
  dmuy    =opt_fun(24)
  ddisp(1)=opt_fun(25)
  ddisp(2)=opt_fun(26)
  ddisp(3)=opt_fun(27)
  ddisp(4)=opt_fun(28)

  !---- Summary data for non-periodic case.
  if(get_option('twiss_inval ') .ne. 0) then
     detl = rt(5,5) * rt(6,6) - rt(5,6) * rt(6,5)
     f = one / sqrt(detl)
     call m66scl(f,rt,frt)
     call m66scl(f,rtp,frtp)
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
     alfa  =zero
     gamtr =zero
     cosmux=zero
     cosmuy=zero

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
     eta = - sd * betas**2 / suml
     alfa = one / gammas**2 + eta
     if (alfa.gt.zero) then
        gamtr = sqrt(one / alfa)
     else if (alfa .eq. zero) then
        gamtr = zero
     else
        gamtr = - sqrt(- one / alfa)
     endif
  endif

  !---- Initialization transverse part
  !---  fix length problem - HG 14.4.08
  suml    = currpos
  !     suml    = currpos-get_value('sequence ','range_start ')
  !---  end of fix length problem - HG 14.4.08
  betx    = opt_fun0(3)
  alfx    = opt_fun0(4)
  amux    = opt_fun(5)
  bety    = opt_fun0(6)
  alfy    = opt_fun0(7)
  amuy    = opt_fun(8)
  qx = amux / twopi
  qy = amuy / twopi

  !---- Adjust values
  orbit5 = -opt_fun0(13)
  xcomax = xcomax
  sigxco = sigxco
  ycomax = ycomax
  sigyco = sigyco

  !     if(opt_fun0(29).ne.zero.or.opt_fun0(30).ne.zero.or.             &
  !      opt_fun0(31).ne.zero.or.opt_fun0(32).ne.zero) then
  !     call fort_warn('Chromaticity calculation wrong due to coupling',&
  !      'use manual derivative or madxp')
  !     frs         xix=zero ! The punishment was too harsh!!!!
  !     xiy=zero
  !     endif

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

end SUBROUTINE tw_summ


SUBROUTINE tmmap(code,fsec,ftrk,orbit,fmap,ek,re,te)
  use twtrrfi
  use name_lenfi
  use time_varfi
  implicit none

  !----------------------------------------------------------------------*
  !     purpose:                                                         *
  !     transport map for a complete element.                            *
  !     optionally, follow orbit.                                        *
  !     input:                                                           *
  !     code                element type code                            *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
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
  integer code
  logical fsec, ftrk, fmap
  double precision node_value,el,orbit(6),ek(6),re(6,6),te(6,6,6)
  double precision plot_tilt,zero
  parameter (zero = 0.d0)
  !---- Initialization
  time_var_p=.false.     
  call dzero(ek,6)
  call m66one(re)
  call dzero(te,216)
  plot_tilt=zero
  fmap=.false.
  el = node_value('l ')
  !---- Select element type.
  go to ( 10,  20,  30,  40,  50,  60,  70,  80,  90, 100,      &
       110, 120, 130, 140, 150, 160, 170, 180, 190, 200,        &
       210, 220, 230, 240, 250, 260,  10, 280, 290, 310,        &
       310, 310, 300, 310, 310, 310, 310, 310, 310, 310,	&
       310, 420, 430), code  ! Enable non-linear thin lens and RF-Multipole  
  !     310, 310, 310), ! Disable non-linear thin lens and RF-Multipole
   
  !---- Drift space, monitor, collimator, or beam instrument.
10 continue
170 continue
180 continue
190 continue
200 continue
210 continue
240 continue
  call tmdrf(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Bending magnet.
20 continue
30 continue
  call tmbend(ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Arbitrary matrix.
40 continue
  call tmarb(fsec,ftrk,orbit,fmap,ek,re,te)
  go to 500

  !---- Quadrupole.
50 continue
  call tmquad(fsec,ftrk,plot_tilt,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Sextupole.
60 continue
  call tmsext(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Octupole.
70 continue
  call tmoct(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Multipole.
80 continue
  call tmmult(fsec,ftrk,orbit,fmap,re,te)
  go to 500

  !---- Solenoid.
90 continue
  call tmsol(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- RF cavity.
100 continue
  call tmrf(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Electrostatic separator.
110 continue
  call tmsep(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Rotation around s-axis.
120 continue
  call tmsrot(ftrk,orbit,fmap,ek,re,te)
  go to 500

  !---- Rotation around y-axis.
130 continue
  call tmyrot(ftrk,orbit,fmap,ek,re,te)
  go to 500

  !---- Correctors.
140 continue
150 continue
160 continue
  call tmcorr(fsec,ftrk,orbit,fmap,el,ek,re,te)
  go to 500

  !---- Beam-beam.
220 continue
  !     (Particles/bunch taken for the opposite beam).
  call tmbb(fsec,ftrk,orbit,fmap,re,te)
  go to 500

  !---- Lump.
230 continue
  go to 500

  !---- Marker.
250 continue
  go to 500

  !---- General bend (dipole, quadrupole, and skew quadrupole).
260 continue
  go to 500

  !---- LCAV cavity.
  continue

  !---- Reserved.
280 continue
290 continue
  go to 500
300 call tmdpdg(ftrk,orbit,fmap,ek,re,te)
  go to 500

  !---- non-linear thin lens
420 continue
  call tmnll(fsec,ftrk,orbit,fmap,ek,re,te)
  go to 500

  !---- RF-Multipole.
430 continue
  call tmrfmult(fsec,ftrk,orbit,fmap,ek,re,te)
  go to 500

  !---- User-defined elements.
310 continue

  !---- End of element calculation;
500 continue
end SUBROUTINE tmmap

SUBROUTINE tmbend(ftrk,orbit,fmap,el,ek,re,te)
  use twtrrfi
  use twisslfi
  use twiss_elpfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for sector bending magnets                         *
  !     Input:                                                           *
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
  logical ftrk,fmap,cplxy,dorad
  integer nd,n_ferr,node_fd_errors,code
  double precision orbit(6),f_errors(0:maxferr),ek(6),re(6,6),      &
       te(6,6,6),rw(6,6),tw(6,6,6),x,y,deltap,field(2,0:maxmul),fintx,   &
       el,tilt,e1,e2,sk1,sk2,h1,h2,hgap,fint,sks,an,h,dh,corr,ek0(6),ct, &
       st,hx,hy,rfac,arad,gamma,pt,rhoinv,blen,node_value,get_value,bvk, &
       el0,orbit0(6),zero,one,two,three
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0)

  !---- Initialize.
  ct=zero
  st=zero
  code = node_value('mad8_type ')
  if(code.eq.39) code=15
  if(code.eq.38) code=24
  deltap=zero
  call dzero(ek0,6)
  call m66one(rw)
  call dzero(tw,216)

  !---- Test for non-zero length.

  fmap = el .ne. zero
  if (fmap) then
     call dzero(f_errors,maxferr+1)
     n_ferr = node_fd_errors(f_errors)
     !-- get element parameters
     elpar_vl = el_par_vector(b_k3s, g_elpar)
     bvk = node_value('other_bv ')
     arad = get_value('probe ','arad ')
     deltap = get_value('probe ','deltap ')
     gamma = get_value('probe ','gamma ')
     dorad = get_value('probe ','radiate ') .ne. zero
     an = bvk * g_elpar(b_angle)
     tilt = g_elpar(b_tilt)
     e1 = g_elpar(b_e1)
     e2 = g_elpar(b_e2)

     if(code.eq.2) then
        e1 = e1 + an / two
        e2 = e2 + an / two
     endif

     !---  bvk also applied further down

     sk1 = g_elpar(b_k1)
     sk2 = g_elpar(b_k2)
     h1 = g_elpar(b_h1)
     h2 = g_elpar(b_h2)
     hgap = g_elpar(b_hgap)
     fint = g_elpar(b_fint)
     fintx = g_elpar(b_fintx)
     sks = g_elpar(b_k1s)
     h = an / el

      !---- Apply field errors and change coefficients using DELTAP.
     if (n_ferr .gt. 0) then
        nd = n_ferr
        call dzero(field,nd)
        call dcopy(f_errors,field,n_ferr)
        dh = (- h * deltap + bvk * field(1,0) / el) / (one + deltap)
        sk1 = (sk1 + field(1,1) / el) / (one + deltap)
        sk2 = (sk2 + field(1,2) / el) / (one + deltap)
        sks = (sks + field(2,1) / el) / (one + deltap)
     else
        dh = - h * deltap / (one + deltap)
        sk1 = sk1 / (one + deltap)
        sk2 = sk2 / (one + deltap)
        sks = sks / (one + deltap)
     endif
     sk1 = bvk * sk1
     sk2 = bvk * sk2
     sks = bvk * sks

     !---- Half radiation effects at entrance.
     if (ftrk .and. dorad) then
        ct = cos(tilt)
        st = sin(tilt)
        x =   orbit(1) * ct + orbit(3) * st
        y = - orbit(1) * st + orbit(3) * ct
        hx = h + dh + sk1*(x - h*y**2/two) + sks*y +                  &
             sk2*(x**2 - y**2)/two
        hy = sks * x - sk1*y - sk2*x*y
        rfac = (arad * gamma**3 * el / three)                         &
             * (hx**2 + hy**2) * (one + h*x) * (one - tan(e1)*x)
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- Body of the dipole.
     !---- centre option
     if(centre_cptk.or.centre_bttk) then
        call dcopy(orbit,orbit00,6)
        call dcopy(ek,ek00,6)
        call dcopy(re,re00,36)
        call dcopy(te,te00,216)
        el0=el/two
        call tmsect(.true.,el0,h,dh,sk1,sk2,ek,re,te)
        !---- Fringe fields.
        corr = (h + h) * hgap * fint
        call tmfrng(.true.,h,sk1,e1,h1,one,corr,rw,tw)
        call tmcat1(.true.,ek,re,te,ek0,rw,tw,ek,re,te)
        !---- Apply tilt.
        if (tilt .ne. zero) then
           call tmtilt(.true.,tilt,ek,re,te)
           cplxy = .true.
        endif
        !---- Track orbit.
        call dcopy(orbit,orbit0,6)
        if (ftrk) call tmtrak(ek,re,te,orbit0,orbit0)
        if(centre_cptk) call twcptk(re,orbit0)
        if(centre_bttk) call twbttk(re,te)
        call dcopy(orbit00,orbit,6)
        call dcopy(ek00,ek,6)
        call dcopy(re00,re,36)
        call dcopy(te00,te,216)
     endif
     !---- End

     call tmsect(.true.,el,h,dh,sk1,sk2,ek,re,te)

     !---- Fringe fields.
     corr = (h + h) * hgap * fint
     call tmfrng(.true.,h,sk1,e1,h1,one,corr,rw,tw)
     call tmcat1(.true.,ek,re,te,ek0,rw,tw,ek,re,te)
     !---- Tor: use FINTX if set
     if (fintx .ge. 0) then
        corr = (h + h) * hgap * fintx
     else
        corr = (h + h) * hgap * fint
     endif
     call tmfrng(.true.,h,sk1,e2,h2,-one,corr,rw,tw)
     call tmcat1(.true.,ek0,rw,tw,ek,re,te,ek,re,te)

     !---- Apply tilt.
     if (tilt .ne. zero) then
        call tmtilt(.true.,tilt,ek,re,te)
        cplxy = .true.
     endif

     !---- Track orbit.
     if (ftrk) then
        call tmtrak(ek,re,te,orbit,orbit)

        !---- Half radiation effects at exit.
        if (ftrk .and. dorad) then
           x =   orbit(1) * ct + orbit(3) * st
           y = - orbit(1) * st + orbit(3) * ct
           hx = h + dh + sk1*(x - h*y**2/two) + sks*y +                &
                sk2*(x**2 - y**2)/two
           hy = sks * x - sk1*y - sk2*x*y
           rfac = (arad * gamma**3 * el / three)                       &
                * (hx**2 + hy**2) * (one + h*x) * (one - tan(e2)*x)
           pt = orbit(6)
           orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
           orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
           orbit(6) = orbit(6) - rfac * (one + pt) ** 2
        endif
     endif
  endif

  !---- Tor: set parameters for sychrotron integral calculations
!  rhoinv = h
!  blen = el

end SUBROUTINE tmbend
SUBROUTINE tmsect(fsec,el,h,dh,sk1,sk2,ek,re,te)

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
  logical fsec
  double precision beta,gamma,dtbyds,bi,bi2,bi2gi2,cm,cp,cx,cy,cyy, &
       dd,dh,difsq,dm,dp,dx,dyy,ek(6),el,fm,fp,fx,fyy,gx,h,h2,hx,re(6,6),&
       sk1,sk2,sm,sp,sumsq,sx,sy,syy,t1,t116,t126,t166,t2,t216,t226,t266,&
       t336,t346,t436,t446,t5,t516,t526,t566,te(6,6,6),xk,xkl,xklsq,xksq,&
       xs6,y0,y1,y2,y2klsq,y2ksq,yk,ykl,yklsq,yksq,ys2,zc,zd,zf,zs,      &
       get_value,zero,one,two,three,four,six,nine,twelve,fifteen,twty,   &
       twty2,twty4,thty,foty2,fvty6,svty2,httwty,c1,c2,c3,c4,s1,s2,s3,s4,&
       cg0,cg1,cg2,ch0,ch1,ch2
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,four=4d0,six=6d0,    &
       nine=9d0,twelve=12d0,fifteen=15d0,twty=20d0,twty2=22d0,twty4=24d0,&
       thty=30d0,foty2=42d0,fvty6=56d0,svty2=72d0,httwty=120d0,c1=one,   &
       c2=one/two,c3=one/twty4,c4=one/720d0,s1=one,s2=one/six,           &
       s3=one/httwty,s4=one/5040d0,cg0=one/twty,cg1=5d0/840d0,           &
       cg2=21d0/60480d0,ch0=one/fvty6,ch1=14d0/4032d0,ch2=147d0/443520d0)

  !---- Initialize.
  call dzero(ek, 6)
  call m66one(re)
  if (fsec) call dzero(te, 216)
  beta = get_value('probe ','beta ')
  gamma = get_value('probe ','gamma ')
  dtbyds = get_value('probe ','dtbyds ')
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
     t516 = h * xs6 * (three*dx*fx - four*gx) +                      &
          (sk1/two) * (fx + sx*dx)
     t526 = h * xs6 * (dx**3 - two*sx*gx) + (sk1/two) * dx**2
     t566 = h * xs6 * (three*hx - two*dx*gx) +                       &
          (sk1/two) * gx - fx
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
     te(5,1,1) = (h2*xs6 * (sx*dx + three*fx) -                      &
          (sk1/four) * (el - cx*sx)) * bi
     te(5,1,2) = (h2*xs6*dx**2 + (sk1/four)*sx**2) * bi
     te(5,2,2) = (h*xs6*gx - sk1 * (fx - sx*dx) / four - sx/two) * bi
     te(5,1,6) = h2 * ((t516 - sk1 * (el*dx - fx) / two) * bi2 +     &
          sx * bi2gi2)
     te(5,2,6) = h2 * ((t526 - sk1 * (dx**2 - sx*fx) / two) * bi2 +  &
          dx * bi2gi2)
     te(5,6,6) = (h**2 * (t566 + t5) * bi2 +                         &
          (three/two) * (h**2*fx - el) * bi2gi2) * bi

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
        re(5,6) = re(5,6) -                                           &
             dh * h * ((two*t566 + t5) * bi2 + fx * bi2gi2)
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
  integer i,k,l
  double precision t(6,6,6)

  do k = 1, 5
     do l = k+1, 6
        do i = 1, 6
           t(i,l,k) = t(i,k,l)
        enddo
     enddo
  enddo

end SUBROUTINE tmsymm
SUBROUTINE tmfrng(fsec,h,sk1,edge,he,sig,corr,re,te)

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
  logical fsec
  double precision corr,edge,h,he,hh,psip,re(6,6),secedg,sig,sk1,   &
       tanedg,te(6,6,6),zero,one
  parameter(zero=0d0,one=1d0)

  call m66one(re)
  if (fsec) call dzero(te, 216)
  !---- Linear terms.
  tanedg = tan(edge)
  secedg = one / cos(edge)
  psip = edge - corr * secedg * (one + sin(edge)**2)
  re(2,1) = + h * tanedg
  re(4,3) = - h * tan(psip)

  !---- Second-order terms.
  if (fsec) then
     hh = sig * (h/2)
     te(1,1,1) = - hh * tanedg**2
     te(1,3,3) = + hh * secedg**2
     te(2,1,1) = (h/2) * he * secedg**3 + sk1 * tanedg
     te(2,1,2) = - te(1,1,1)
     te(2,3,3) = hh * h * tanedg**3 - te(2,1,1)
     te(2,3,4) = + te(1,1,1)
     te(3,1,3) = - te(1,1,1)
     te(4,1,3) = - te(2,1,1)
     te(4,1,4) = + te(1,1,1)
     te(4,2,3) = - te(1,3,3)
     if (sig .gt. zero) then
        te(2,3,3) = te(2,3,3) + (h*secedg)**2 * tanedg/2
     else
        te(2,1,1) = te(2,1,1) - (h*tanedg)**2 * tanedg/2
        te(4,1,3) = te(4,1,3) + (h*secedg)**2 * tanedg/2
     endif
     call tmsymm(te)
  endif

end SUBROUTINE tmfrng
SUBROUTINE tmcat1(fsec,eb,rb,tb,ea,ra,ta,ed,rd,td)

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
  logical fsec
  integer i,ij,j,k
  double precision ea(6),eb(6),ed(6),ew(6),es(6,6),ra(6,6),rb(6,6), &
       rd(6,6),rw(6,6),ta(6,6,6),td(6,6,6),tw(6,6,6),tb(36,6),ts(36,6),  &
       two
  parameter(two=2d0)

  !---- Second order terms.
  if (fsec) then

     !---- Auxiliary terms.
     do k = 1, 6

        !---- Sum over S of TB(I,S,K) * EA(S).
        do i = 1, 6
           es(i,k) = tb(i   ,k) * ea(1) + tb(i+ 6,k) * ea(2)           &
                + tb(i+12,k) * ea(3) + tb(i+18,k) * ea(4)                         &
                + tb(i+24,k) * ea(5) + tb(i+30,k) * ea(6)
        enddo

        !---- Sum over S of TB(I,J,S) * RA(S,K).
        do ij = 1, 36
           ts(ij,k) = tb(ij,1) * ra(1,k) + tb(ij,2) * ra(2,k)          &
                + tb(ij,3) * ra(3,k) + tb(ij,4) * ra(4,k)                         &
                + tb(ij,5) * ra(5,k) + tb(ij,6) * ra(6,k)
        enddo
     enddo

     !---- Final values.
     do k = 1, 6

        !---- Zero-order terms.
        ew(k) = eb(k) + (rb(k,1) + es(k,1)) * ea(1)                   &
             + (rb(k,2) + es(k,2)) * ea(2)                                     &
             + (rb(k,3) + es(k,3)) * ea(3)                                     &
             + (rb(k,4) + es(k,4)) * ea(4)                                     &
             + (rb(k,5) + es(k,5)) * ea(5)                                     &
             + (rb(k,6) + es(k,6)) * ea(6)

        !---- First-order terms.
        do j = 1, 6
           rw(j,k) = (rb(j,1) + two * es(j,1)) * ra(1,k)               &
                + (rb(j,2) + two * es(j,2)) * ra(2,k)                             &
                + (rb(j,3) + two * es(j,3)) * ra(3,k)                             &
                + (rb(j,4) + two * es(j,4)) * ra(4,k)                             &
                + (rb(j,5) + two * es(j,5)) * ra(5,k)                             &
                + (rb(j,6) + two * es(j,6)) * ra(6,k)
        enddo

        !---- Second-order terms.
        do j = k, 6
           do i = 1, 6
              tw(i,j,k) =                                               &
                   + (rb(i,1)+two*es(i,1))*ta(1,j,k) + ts(i   ,j)*ra(1,k)            &
                   + (rb(i,2)+two*es(i,2))*ta(2,j,k) + ts(i+ 6,j)*ra(2,k)            &
                   + (rb(i,3)+two*es(i,3))*ta(3,j,k) + ts(i+12,j)*ra(3,k)            &
                   + (rb(i,4)+two*es(i,4))*ta(4,j,k) + ts(i+18,j)*ra(4,k)            &
                   + (rb(i,5)+two*es(i,5))*ta(5,j,k) + ts(i+24,j)*ra(5,k)            &
                   + (rb(i,6)+two*es(i,6))*ta(6,j,k) + ts(i+30,j)*ra(6,k)
              tw(i,k,j) = tw(i,j,k)
           enddo
        enddo
     enddo

     !---- Copy second-order terms.
     call dcopy(tw,td,216)

     !---- Second-order not desired.
  else
     do k = 1, 6

        !---- Zero-order terms.
        ew(k) = eb(k) + rb(k,1) * ea(1) + rb(k,2) * ea(2)             &
             + rb(k,3) * ea(3) + rb(k,4) * ea(4)                               &
             + rb(k,5) * ea(5) + rb(k,6) * ea(6)

        !---- First-order terms.
        do j = 1, 6
           rw(j,k) = rb(j,1) * ra(1,k) + rb(j,2) * ra(2,k)             &
                + rb(j,3) * ra(3,k) + rb(j,4) * ra(4,k)                           &
                + rb(j,5) * ra(5,k) + rb(j,6) * ra(6,k)
        enddo
     enddo
  endif

  !---- Copy zero- and first-order terms.
  call dcopy(ew,ed,6)
  call dcopy(rw,rd,36)

end SUBROUTINE tmcat1
SUBROUTINE tmtilt(fsec,tilt,ek,r,t)

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
  logical fsec
  integer i,j,k
  double precision c,ek(6),r(6,6),r1j,r2j,ri1,ri2,s,t(6,6,6),t1jk,  &
       t2jk,ti1k,ti2k,tij1,tij2,tilt,xx

  !YIL13: make sure we don't bother tilting by 0 degrees..
  if (tilt.eq.0) return

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
SUBROUTINE tmcorr(fsec,ftrk,orbit,fmap,el,ek,re,te)

  use twtrrfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for orbit correctors.                              *
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
  logical fsec,ftrk,fmap,cplxy,dorad
  integer i, n_ferr,code,node_fd_errors
  double precision orbit(6),f_errors(0:maxferr),ek(6),re(6,6),      &
       te(6,6,6),deltap,gamma,arad,el,rfac,pt,xkick,ykick,dpx,dpy,       &
       node_value,get_value,bvk,field(2),tilt,xau,div,zero,one,three,half
  parameter(zero=0d0,one=1d0,three=3d0,half=5d-1)

  !---- Initialize.
  rfac=zero
  call dzero(f_errors,maxferr+1)
  n_ferr = node_fd_errors(f_errors)
  if (el .eq. zero)  then
     div = one
  else
     div = el
  endif
  bvk = node_value('other_bv ')
  deltap = get_value('probe ','deltap ')
  arad = get_value('probe ','arad ')
  arad = get_value('probe ','arad ')
  gamma = get_value('probe ','gamma ')
  dorad = get_value('probe ','radiate ') .ne. zero
  tilt = -node_value('tilt ')

  !---- Tracking desired, use corrector map.
  if (ftrk) then
     do i = 1, 2
        field(i) = zero
     enddo
     if (n_ferr .gt. 0) call dcopy(f_errors, field, min(2, n_ferr))
     !---- Original setting.
     code = node_value('mad8_type ')

     if(code.eq.39) code=15
     if(code.eq.38) code=24

     if(code.eq.14) then
        xkick=bvk*(node_value('kick ')+node_value('chkick ')+field(1)/div)
        ykick=zero
     else if(code.eq.15) then
        xkick=bvk*(node_value('hkick ')+node_value('chkick ')+field(1)/div)
        ykick=bvk*(node_value('vkick ')+node_value('cvkick ')+field(2)/div)
     else if(code.eq.16) then
        xkick=zero
        ykick=bvk*(node_value('kick ')+node_value('cvkick ')+field(2)/div)
     else
        xkick=zero
        ykick=zero
     endif
     xau=xkick
     xkick= xkick*cos(tilt)+ykick*sin(tilt)
     ykick=-xau*sin(tilt)+ykick*cos(tilt)
     !---- Sum up total kicks.
     dpx = xkick / (one + deltap)
     dpy = ykick / (one + deltap)
     if (dpy .ne. zero) cplxy = .true.

     !---- Half kick at entrance.
     orbit(2) = orbit(2) + half * dpx
     orbit(4) = orbit(4) + half * dpy

     !---- Half radiation effects at entrance.
     if (dorad  .and.  el .ne. zero) then
        rfac = arad * gamma**3 * (dpx**2 + dpy**2) / (three * el)
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- Drift to end.
     if (el .ne. zero) then
        call tmdrf(fsec,ftrk,orbit,fmap,el,ek,re,te)
     endif
     fmap = .true.

     !---- Half radiation effects at exit.
     if (dorad  .and.  el .ne. zero) then
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- Half kick at exit.
     orbit(2) = orbit(2) + half * dpx
     orbit(4) = orbit(4) + half * dpy

     !---- No orbit track desired, use drift map.
  else
     call tmdrf(fsec,ftrk,orbit,fmap,el,ek,re,te)
  endif

end SUBROUTINE tmcorr
SUBROUTINE tmmult(fsec,ftrk,orbit,fmap,re,te)

  use twtrrfi
  use twisslfi
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
  logical fsec,ftrk,fmap,dorad
  integer n_ferr,nord,iord,j,nd,nn,ns,node_fd_errors,get_option
  double precision orbit(6),f_errors(0:maxferr),re(6,6),te(6,6,6),x,&
       y,dbr,dbi,dipr,dipi,dr,di,drt,dpx,dpy,elrad,beta,bi,deltap,       &
       vals(2,0:maxmul),field(2,0:maxmul),normal(0:maxmul),orbit0(6),    &
       skew(0:maxmul),node_value,get_value,pt,arad,gammas,rfac,bvk,dpxr, &
       dpyr,zero,one,two,three, tilt, angle, dtmp
  double precision orbit00(6),re00(6,6),te00(6,6,6)
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0)

  !---- Initialize
  rfac=zero
  call dzero(f_errors,maxferr+1)
  n_ferr = node_fd_errors(f_errors)
  bvk = node_value('other_bv ')
  !---- Multipole length for radiation.
  elrad = node_value('lrad ')
  dorad = get_value('probe ','radiate ') .ne. zero
  arad = get_value('probe ','arad ')
  gammas= get_value('probe ','gamma ')
  deltap = get_value('probe ', 'deltap ')
  beta = get_value('probe ','beta ')
  fmap = .true.
  bi = one / beta

  !---- Multipole components.
  call dzero(normal,maxmul+1)
  call dzero(skew,maxmul+1)
  call get_node_vector('knl ',nn,normal)
  call get_node_vector('ksl ',ns,skew)
  tilt = node_value('tilt ')
  call dzero(vals,2*(maxmul+1))
  do iord = 0, nn
     vals(1,iord) = normal(iord)
  enddo
  do iord = 0, ns
     vals(2,iord) = skew(iord)
  enddo
  !---- Field error vals.
  call dzero(field,2*(maxmul+1))
  if (n_ferr .gt. 0) then
     call dcopy(f_errors,field,n_ferr)
  endif
  nd = 2 * max(nn, ns, n_ferr/2-1)

  !---- Dipole error.
  dbr = field(1,0) / (one + deltap)
  dbi = field(2,0) / (one + deltap)

  !---- Nominal dipole strength.
  dipr = vals(1,0) / (one + deltap)
  dipi = vals(2,0) / (one + deltap)

  if (tilt .ne. zero)  then
     if(dipi.ne.zero.or.dipr.ne.zero) then
        angle = atan2(dipi, dipr) - tilt
     else
        angle = -tilt
     endif
     dtmp = sqrt(dipi**2+dipr**2)
     dipr = dtmp * cos(angle)
     dipi = dtmp * sin(angle)
     dtmp = sqrt(dbi**2+dbr**2)
     dbr = dtmp * cos(angle)
     dbi = dtmp * sin(angle)
  endif

  dbr = bvk * dbr
  dbi = bvk * dbi
  dipr = bvk * dipr
  dipi = bvk * dipi

  !---- Other components and errors.
  nord = 0
  do iord = 1, nd/2
     do j = 1, 2
        if (field(j,iord) .ne. zero)  nord = iord
        if (vals(j,iord) .ne. zero)  nord = iord
     enddo
  enddo
  do iord = 1, nord
     do j = 1, 2
        field(j,iord) = (vals(j,iord) + field(j,iord))                &
             / (one + deltap)
     enddo
     if (tilt .ne. zero)  then
        if(field(2,iord).ne.zero.or.field(1,iord).ne.zero) then
           angle = atan2(field(2,iord), field(1,iord))/(iord+1) - tilt
        else
           angle = -tilt
        endif
        dtmp = sqrt(field(1,iord)**2+field(2,iord)**2)
        angle = (iord+1) * angle
        field(1,iord) = dtmp * cos(angle)
        field(2,iord) = dtmp * sin(angle)
     endif
     do j = 1, 2
        field(j,iord) = bvk * field(j,iord)
     enddo
  enddo
  !---- Track orbit.
  if (ftrk) then
     x = orbit(1)
     y = orbit(3)

     !---- Multipole kick.
     dr = zero
     di = zero
     do iord = nord, 1, -1
        drt = (dr * x - di * y) / (iord+1) + field(1,iord)
        di  = (dr * y + di * x) / (iord+1) + field(2,iord)
        dr  = drt
     enddo
     dpx = dbr + (dr * x - di * y)
     dpy = dbi + (di * x + dr * y)


     !---- Radiation effects at entrance.
     if (dorad  .and.  elrad .ne. zero) then
        dpxr = dpx + dipr
        dpyr = dpy + dipi
        rfac = arad * gammas**3 * (dpxr**2+dpyr**2) / (three*elrad)
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- Track orbit.
     orbit(2) = orbit(2) - dpx + dipr * (deltap + bi*orbit(6))
     orbit(4) = orbit(4) + dpy - dipi * (deltap + bi*orbit(6))
     orbit(5) = orbit(5) - (dipr*x + dipi*y) * bi
     !---- Add the missing focussing component of thin dipoles for co
     if(elrad.gt.zero.and.get_option('thin_foc ').eq.1) then
        orbit(2) = orbit(2)-dipr*dipr/elrad*x
        orbit(4) = orbit(4)-dipi*dipi/elrad*y
     endif

     !---- Radiation effects at exit.
     if (dorad  .and.  elrad .ne. zero) then
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- Orbit not wanted.
  else
     x = zero
     y = zero
     nord = min(nord, 2)
  endif

  !---- First-order terms (use X,Y from orbit tracking).
  if (nord .ge. 1) then
     dr = zero
     di = zero
     do iord = nord, 1, -1
        drt = (dr * x - di * y) / (iord) + field(1,iord)
        di  = (dr * y + di * x) / (iord) + field(2,iord)
        dr  = drt
     enddo
     re(2,1) = - dr
     re(2,3) = + di
     re(4,1) = + di
     re(4,3) = + dr
  endif
  !---- Add the missing focussing component of thin dipoles
  if(elrad.gt.zero.and.get_option('thin_foc ').eq.1) then
     re(2,1)=re(2,1)-dipr*dipr/elrad
     re(4,3)=re(4,3)-dipi*dipi/elrad
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
           drt = (dr * x - di * y) / (iord-1) + field(1,iord)
           di  = (dr * y + di * x) / (iord-1) + field(2,iord)
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
  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     if(centre_cptk) then
        call dcopy(orbit,orbit0,6)
        call twcptk(re,orbit0)
     endif
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif

end SUBROUTINE tmmult
SUBROUTINE tmoct(fsec,ftrk,orbit,fmap,el,ek,re,te)

  use twtrrfi
  use twisslfi
  use twiss_elpfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for octupole element.                              *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
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
  logical fsec,ftrk,fmap,cplxy,dorad
  integer i,j,n_ferr,node_fd_errors
  double precision orbit(6),f_errors(0:maxferr),ek(6),re(6,6),      &
       te(6,6,6),rw(6,6),tw(6,6,6),deltap,el,sk3,sk3l,rfac,arad,gamma,pt,&
       octr,octi,posr,posi,cr,ci,tilt4,node_value,get_value,sk3s,bvk,    &
       field(2,0:3),el0,orbit0(6),zero,one,two,three,four,six
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,four=4d0,six=6d0)

  !---- Initialize.
  deltap = zero
  fmap = el .ne. zero
  if (.not. fmap) return
  call dzero(f_errors,maxferr+1)
  n_ferr = node_fd_errors(f_errors)
  bvk = node_value('other_bv ')
  !-- get element parameters
  elpar_vl = el_par_vector(o_k3s, g_elpar)
  !---- Set up half octupole strength.
  if (ftrk) then
     !---- Field error.
     do i = 0,3
        do j = 1, 2
           field(j,i) = zero
        enddo
     enddo
     if (n_ferr .gt. 0) call dcopy(f_errors, field, min(8, n_ferr))
     arad = get_value('probe ','arad ')
     gamma = get_value('probe ','gamma ')
     deltap = get_value('probe ','deltap ')
     dorad = get_value('probe ','radiate ') .ne. zero
     sk3 = bvk * g_elpar(o_k3)
     sk3s = bvk * g_elpar(o_k3s)
     sk3 = sk3 + bvk * field(1,3)/el
     sk3s = sk3s + bvk * field(2,3)/el
     tilt4 = -4*node_value('tilt ')
     if(sk3s.ne.zero) then
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
     if (dorad) then
        rfac = arad * gamma**3 * (cr**2 + ci**2) / (three * el)
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- First-order terms w.r.t. orbit.
     call m66one(rw)
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
     call dzero(tw,216)
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
     !---- centre option
     if(centre_cptk.or.centre_bttk) then
        call dcopy(orbit,orbit00,6)
        call dcopy(ek,ek00,6)
        call dcopy(re,re00,36)
        call dcopy(te,te00,216)
        el0=el/two
        call dcopy(orbit,orbit0,6)
        call tmdrf0(fsec,ftrk,orbit0,fmap,el0,ek,re,te)
        call tmcat(fsec,re,te,rw,tw,re,te)
        if(centre_cptk) call twcptk(re,orbit0)
        if(centre_bttk) call twbttk(re,te)
        call dcopy(orbit00,orbit,6)
        call dcopy(ek00,ek,6)
        call dcopy(re00,re,36)
        call dcopy(te00,te,216)
     endif
     call tmdrf0(fsec,ftrk,orbit,fmap,el,ek,re,te)
     call tmcat(fsec,re,te,rw,tw,re,te)

     !---- Half kick at exit.
     posr = orbit(1) * (orbit(1)**2 - three*orbit(3)**2) / six
     posi = orbit(3) * (three*orbit(1)**2 - orbit(3)**2) / six
     cr = octr * posr - octi * posi
     ci = octr * posi + octi * posr
     orbit(2) = orbit(2) - cr / two
     orbit(4) = orbit(4) + ci / two

     !---- Half radiation effects.
     if (dorad) then
        rfac = arad * gamma**3 * (cr**2 + ci**2) / (three * el)
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
     endif

     !---- First-order terms w.r.t. orbit.
     call m66one(rw)
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
     call dzero(tw,216)
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

     !---- No orbit track requested, use drift map.
  else
     call tmdrf(fsec,ftrk,orbit,fmap,el,ek,re,te)
  endif

end SUBROUTINE tmoct
SUBROUTINE tmarb(fsec,ftrk,orbit,fmap,ek,re,te)

  use trackfi
  use twtrrfi
  use name_lenfi
  use time_varfi  
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
  logical fsec,ftrk,fmap,time_var
  integer i1,i2,i3,mylen
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),node_value,     &
       test,zero
  character(name_len) name
  parameter(zero=0d0)

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
  if(time_var.and.time_var_p) then
     time_var_p_cnt=time_var_p_cnt+1
     time_var_p_lnt=time_var_p_lnt+1
     if(idnint(time_var_p_ind(time_var_p_cnt)).ne.time_var_p_lnt)    &
          call aafail('TMARB: ', 'wrong index in Table: time_var_pha')
     call element_name(name,len(name))
     mylen=len_trim(name)
     if(time_var_p_ch(time_var_p_cnt)(:mylen).ne.name(:mylen))       & 
          call aafail('TMARB: ', 'wrong element name in Table: '//   &
          'time_var_pha')
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
  endif

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
     test=zero
     do i1=1,6
        do i2=1,6
           do i3=1,6
              test=test+abs(te(i1,i2,i3))
           enddo
        enddo
     enddo
     if(test.ne.zero) fsecarb=.true.
  endif

  !---- Track orbit.
  fmap = .true.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmarb
SUBROUTINE tmquad(fsec,ftrk,plot_tilt,orbit,fmap,el,ek,re,te)

  use twtrrfi
  use twisslfi
  use twiss_elpfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for quadrupole element.                            *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     plot_tilt (double)  external tilt needed for plot                *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is rotated (in tmmap) before and after this    *
  !     routine                                                          *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical fsec,ftrk,fmap,cplxy,dorad
  integer i,j,n_ferr,node_fd_errors
  double precision ct, st, tmp
  double precision orbit(6),orbit0(6),f_errors(0:maxferr),ek(6),    &
       re(6,6),te(6,6,6),deltap,el,el0,tilt,sk1,rfac,arad,gamma,pt,sk1s, &
       bvk,field(2,0:1),node_value,get_value,plot_tilt,zero,one,two,three
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0)

  !---- Initialize.
  st=zero
  ct=zero
  cplxy=.false.
  fmap = el .ne. zero
  if (.not. fmap) return
  !---- Field error.
  call dzero(f_errors,maxferr+1)
  n_ferr = node_fd_errors(f_errors)
  do i = 0, 1
     do j = 1, 2
        field(j,i) = zero
     enddo
  enddo
  if (n_ferr .gt. 0) call dcopy(f_errors, field, min(4,n_ferr))
  !-- element paramters
  elpar_vl = el_par_vector(q_k1s, g_elpar)
  bvk = node_value('other_bv ')
  sk1 = bvk * g_elpar(q_k1)
  sk1s = bvk * g_elpar(q_k1s)
  sk1 = sk1 + bvk * field(1,1)/el
  sk1s = sk1s + bvk * field(2,1)/el
  tilt = g_elpar(q_tilt)
  if(sk1s.ne.zero) then
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
  tilt=tilt+plot_tilt
  arad = get_value('probe ','arad ')
  gamma = get_value('probe ','gamma ')
  deltap = get_value('probe ','deltap ')
  dorad = get_value('probe ','radiate ') .ne. zero
  sk1 = sk1 / (one + deltap)

  !---- Half radiation effect at entry.
  if (dorad.and.ftrk) then
     rfac = (arad * gamma**3 * sk1**2 * el / three) * (orbit(1)**2 + &
          orbit(3)**2)
     pt = orbit(6)
     orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
     orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
     orbit(6) = orbit(6) - rfac * (one + pt) ** 2
  endif

  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     el0=el/two
     call dcopy(orbit,orbit0,6)
     call qdbody(fsec,ftrk,tilt,sk1,orbit0,el0,ek,re,te)
     if(centre_cptk) call twcptk(re,orbit0)
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif
  call qdbody(fsec,ftrk,tilt,sk1,orbit,el,ek,re,te)

  !---- Half radiation effect at exit.
  if (dorad.and.ftrk) then
     rfac = (arad * gamma**3 * sk1**2 * el / three) * (orbit(1)**2 + &
          orbit(3)**2)
     pt = orbit(6)
     orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
     orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
     orbit(6) = orbit(6) - rfac * (one + pt) ** 2
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
  logical fsec,ftrk
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),el,tilt,sk1,    &
       gamma,beta,dtbyds,qk,qkl,qkl2,cx,sx,cy,sy,biby4,get_value,zero,   &
       one,two,four,six,ten3m
  parameter(zero=0d0,one=1d0,two=2d0,four=4d0,six=6d0,ten3m=1d-3)
  !---- Initialize.
  beta = get_value('probe ','beta ')
  gamma = get_value('probe ','gamma ')
  dtbyds = get_value('probe ','dtbyds ')

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
SUBROUTINE tmsep(fsec,ftrk,orbit,fmap,el,ek,re,te)

  use twisslfi
  use twiss_elpfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for electrostatic separator.                       *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is rotated (in tmmap) before and after this    *
  !     routine                                                          *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical fsec,ftrk,fmap,cplxy
  double precision ct, st, tmp
  double precision orbit(6),orbit0(6),ek(6),re(6,6),te(6,6,6),      &
       deltap,el,el0,tilt,ekick,charge,pc,efield,exfld,eyfld,node_value, &
       get_value,zero,one,two,ten3m
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0,ten3m=1d-3)

  !---- Initialize.
  st=zero
  ct=zero
  cplxy=.false.
  charge = get_value('probe ','charge ')
  pc = get_value('probe ','pc ')
  deltap = get_value('probe ','deltap ')

  fmap = el .ne. zero
  if (.not. fmap) return
  if (ftrk) then
     !-- get element parameters
     elpar_vl = el_par_vector(e_ey, g_elpar)
     !---- Strength and tilt.
     exfld = g_elpar(e_ex)
     eyfld = g_elpar(e_ey)
     tilt = g_elpar(e_tilt)
     if(eyfld.ne.zero) then
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
  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     el0=el/two
     call dcopy(orbit,orbit0,6)
     call spbody(fsec,ftrk,tilt,ekick,orbit0,el0,ek,re,te)
     if(centre_cptk) call twcptk(re,orbit0)
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif
  call spbody(fsec,ftrk,tilt,ekick,orbit,el,ek,re,te)
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
  logical fsec,ftrk
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),el,tilt,ekl,ch, &
       sh,sy,ekick,beta,gamma,dtbyds,dy,fact,get_value,zero,one,two,     &
       three,by2,by24,by6,eps
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,by2=1d0/2d0,         &
       by6=1d0/6d0,by24=1d0/24d0,eps=1d-4)

  !---- Initialize.
  dtbyds = get_value('probe ','dtbyds ')
  gamma = get_value('probe ','gamma ')
  beta = get_value('probe ','beta ')

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
SUBROUTINE tmsext(fsec,ftrk,orbit,fmap,el,ek,re,te)

  use twtrrfi
  use twisslfi
  use twiss_elpfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for sextupole element.                             *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Remark: the orbit is rotated (in tmmap) before and after this    *
  !     routine                                                          *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     el        (double)  element length.                              *
  !     ek(6)     (double)  kick due to element.                         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  logical fsec,ftrk,fmap,cplxy,dorad
  integer i,j,n_ferr,node_fd_errors
  double precision ct, st, tmp
  double precision orbit(6),orbit0(6),f_errors(0:maxferr),ek(6),    &
       re(6,6),te(6,6,6),deltap,el,el0,tilt,sk2,rfac,arad,gamma,pt,sk2s, &
       bvk,field(2,0:2),node_value,get_value,zero,one,two,three,twelve
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,twelve=12d0)

  !---- Initialize.
  st=zero
  ct=zero
  cplxy=.false.
  fmap = el .ne. zero
  if (.not. fmap) return
  !---- Field error.
  call dzero(f_errors,maxferr+1)
  n_ferr = node_fd_errors(f_errors)
  do i = 0, 2
     do j = 1,2
        field(j,i) = zero
     enddo
  enddo
  if (n_ferr .gt. 0) call dcopy(f_errors, field, min(6,n_ferr))
  !-- get element parameters
  elpar_vl = el_par_vector(s_k2s, g_elpar)
  bvk = node_value('other_bv ')
  sk2 = bvk * g_elpar(s_k2)
  sk2s = bvk * g_elpar(s_k2s)
  sk2 = sk2 + bvk * field(1,2)/el
  sk2s = sk2s + bvk * field(2,2)/el
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
  cplxy = cplxy .or. tilt .ne. zero
  arad = get_value('probe ','arad ')
  gamma = get_value('probe ','gamma ')
  deltap = get_value('probe ','deltap ')
  dorad = get_value('probe ','radiate ') .ne. zero
  sk2 = sk2 / (one + deltap)

  !---- Half radiation effects at entrance.
  if (ftrk .and. dorad) then
     rfac = arad * gamma**3 * sk2**2 * el * (orbit(1)**2 + orbit(3)  &
          **2)**2 / twelve
     pt = orbit(6)
     orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
     orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
     orbit(6) = orbit(6) - rfac * (one + pt) ** 2
  endif

  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     el0=el/two
     call dcopy(orbit,orbit0,6)
     call sxbody(fsec,ftrk,tilt,sk2,orbit0,el0,ek,re,te)
     if(centre_cptk) call twcptk(re,orbit0)
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif
  call sxbody(fsec,ftrk,tilt,sk2,orbit,el,ek,re,te)

  !---- Half radiation effects at exit.
  if (ftrk) then
     if (dorad) then
        rfac = arad * gamma**3 * sk2**2 * el * (orbit(1)**2 + orbit   &
             (3)**2)**2 / twelve
        pt = orbit(6)
        orbit(2) = orbit(2) - rfac * (one + pt) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + pt) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + pt) ** 2
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
  logical fsec,ftrk
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),el,tilt,sk2,skl,&
       gamma,dtbyds,beta,s1,s2,s3,s4,get_value,zero,two,three,four
  parameter(zero=0d0,two=2d0,three=3d0,four=4d0)

  !---- Initialize.
  beta = get_value('probe ','beta ')
  gamma = get_value('probe ','gamma ')
  dtbyds = get_value('probe ','dtbyds ')

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
SUBROUTINE tmsol(fsec,ftrk,orbit,fmap,el,ek,re,te)

  use twisslfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for solenoid element, with centre option.          *
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
  logical fsec,ftrk,fmap
  double precision orbit(6),orbit0(6),ek(6),re(6,6),te(6,6,6),      &
       el,el0,zero,two
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  parameter(zero=0d0,two=2d0)

  if(el.eq.zero) then
     call tmsol_th(ftrk,orbit,fmap,ek,re,te)
  else
     !---- centre option
     if(centre_cptk.or.centre_bttk) then
        call dcopy(orbit,orbit00,6)
        call dcopy(ek,ek00,6)
        call dcopy(re,re00,36)
        call dcopy(te,te00,216)
        el0=el/two
        call dcopy(orbit,orbit0,6)
        call tmsol0(fsec,ftrk,orbit0,fmap,el0,ek,re,te)
        if(centre_cptk) call twcptk(re,orbit0)
        if(centre_bttk) call twbttk(re,te)
        call dcopy(orbit00,orbit,6)
        call dcopy(ek00,ek,6)
        call dcopy(re00,re,36)
        call dcopy(te00,te,216)
     endif
     call tmsol0(fsec,ftrk,orbit,fmap,el,ek,re,te)
  endif

end SUBROUTINE tmsol
SUBROUTINE tmsol0(fsec,ftrk,orbit,fmap,el,ek,re,te)

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
  logical fsec,ftrk,fmap,cplxy
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),deltap,el,sks,  &
       sk,gamma,skl,beta,co,si,sibk,temp,dtbyds,node_value,get_value,bvk,&
       zero,one,two,three,six,ten5m
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,six=6d0,ten5m=1d-5)

  !---- Initialize.
  fmap = el .ne. zero
  if (.not. fmap) return
  dtbyds = get_value('probe ','dtbyds ')
  gamma = get_value('probe ','gamma ')
  beta = get_value('probe ','beta ')
  deltap = get_value('probe ','deltap ')

  !---- Strength.
  sks = node_value('ks ')
  if (sks .ne. zero) then
     cplxy = .true.
  endif

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

  !---- First-order terms.
  re(1,1) = co**2
  re(2,2) = re(1,1)
  re(3,3) = re(1,1)
  re(4,4) = re(1,1)
  re(1,2) = co * sibk
  re(3,4) = re(1,2)
  re(1,3) = co * si
  re(2,4) = re(1,3)
  re(3,1) = - re(1,3)
  re(4,2) = re(3,1)
  re(2,1) = sk * re(3,1)
  re(4,3) = re(2,1)
  re(1,4) = si * sibk
  re(3,2) = - re(1,4)
  re(4,1) = sk * si**2
  re(2,3) = - re(4,1)
  re(5,6) = el/(beta*gamma)**2
  ek(5) = el*dtbyds

  !---- Second-order terms.
  if (fsec) then
     temp = el * co * si / beta
     te(1,4,6) = - temp
     te(3,2,6) =   temp
     te(1,1,6) =   temp * sk
     te(2,2,6) =   temp * sk
     te(3,3,6) =   temp * sk
     te(4,4,6) =   temp * sk
     te(2,3,6) =   temp * sk**2
     te(4,1,6) = - temp * sk**2

     temp = el * (co**2 - si**2) / (two * beta)
     te(1,2,6) = - temp
     te(3,4,6) = - temp
     te(1,3,6) = - temp * sk
     te(2,4,6) = - temp * sk
     te(3,1,6) =   temp * sk
     te(4,2,6) =   temp * sk
     te(2,1,6) =   temp * sk**2
     te(4,3,6) =   temp * sk**2

     temp = el / (two * beta)
     te(5,2,2) = - temp
     te(5,4,4) = - temp
     te(5,1,4) =   temp * sk
     te(5,2,3) = - temp * sk
     te(5,1,1) = - temp * sk**2
     te(5,3,3) = - temp * sk**2
     te(5,6,6) = - three * re(5,6) / (two * beta)
     call tmsymm(te)
  endif

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

end SUBROUTINE tmsol0
SUBROUTINE tmsrot(ftrk,orbit,fmap,ek,re,te)

  use twisslfi
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
  logical ftrk,fmap,cplxy
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),theta,ct,st,    &
       orbit0(6),node_value,zero
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  parameter(zero=0d0)

  !---- Initialize.
  theta = node_value('angle ')
  fmap = theta .ne. zero
  if (.not. fmap) return

  !---- First-order terms.
  cplxy = .true.
  ct = cos(theta)
  st = sin(theta)
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
  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     if(centre_cptk) then
        call dcopy(orbit,orbit0,6)
        call twcptk(re,orbit0)
     endif
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif
end SUBROUTINE tmsrot
SUBROUTINE tmyrot(ftrk,orbit,fmap,ek,re,te)

  use twisslfi
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
  logical ftrk,fmap
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),phi,cosphi,     &
       orbit0(6),sinphi,tanphi,beta,node_value,get_value,zero,one
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  parameter(zero=0d0,one=1d0)

  !---- Initialize.
  phi = node_value('angle ')
  fmap = phi .ne. zero
  if (.not. fmap) return
  beta = get_value('probe ','beta ')

  !---- Kick.
  cosphi = cos(phi)
  sinphi = sin(phi)
  tanphi = sinphi / cosphi
  ek(2) = - sinphi

  !---- Transfer matrix.
  re(1,1) = one / cosphi
  re(2,2) = cosphi
  re(2,6) = - sinphi / beta
  re(5,1) = tanphi / beta

  !---- Track orbit.
  if (ftrk) call tmtrak(ek,re,te,orbit,orbit)

  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     if(centre_cptk) then
        call dcopy(orbit,orbit0,6)
        call twcptk(re,orbit0)
     endif
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif

end SUBROUTINE tmyrot
SUBROUTINE tmdrf0(fsec,ftrk,orbit,fmap,dl,ek,re,te)

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
  logical fsec,ftrk,fmap
  double precision dl,beta,gamma,dtbyds,orbit(6),ek(6),re(6,6),     &
       te(6,6,6),get_value,zero,two,three
  parameter(zero=0d0,two=2d0,three=3d0)

  !---- Initialize.
  call dzero(ek, 6)
  call m66one(re)
  if (fsec) call dzero(te, 216)
  fmap = dl .ne. zero
  dtbyds = get_value('probe ', 'dtbyds ')
  gamma = get_value('probe ', 'gamma ')
  beta = get_value('probe ', 'beta ')

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

end SUBROUTINE tmdrf0
SUBROUTINE tmdrf(fsec,ftrk,orbit,fmap,dl,ek,re,te)

  use twisslfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for drift space, with centre option.               *
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
  logical fsec,ftrk,fmap
  double precision dl,dl0,orbit(6),orbit0(6),ek(6),re(6,6),         &
       te(6,6,6),two
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  parameter(two=2d0)

  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     dl0=dl/two
     call dcopy(orbit,orbit0,6)
     call tmdrf0(fsec,ftrk,orbit0,fmap,dl0,ek,re,te)
     if(centre_cptk) call twcptk(re,orbit0)
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif
  call tmdrf0(fsec,ftrk,orbit,fmap,dl,ek,re,te)

end SUBROUTINE tmdrf
SUBROUTINE tmrf(fsec,ftrk,orbit,fmap,el,ek,re,te)

  use twisslfi
  use twiss_elpfi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for RF cavity.                                     *
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
  logical fsec,ftrk,fmap
  double precision orbit(6),orbit0(6),ek(6),re(6,6),te(6,6,6),      &
       rw(6,6),tw(6,6,6),el,rfv,rff,rfl,dl,omega,vrf,phirf,pc,deltap,c0, &
       c1,c2,ek0(6),ten6p,clight,pi,twopi,zero,one,two,half,ten3m,bvk
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  double precision, external :: get_variable,get_value,node_value
  integer, external :: el_par_vector
  integer elpar_vl
  parameter(zero=0d0,one=1d0,two=2d0,half=5d-1,ten6p=1d6,ten3m=1d-3)

  !-- get element parameters
  elpar_vl = el_par_vector(r_freq, g_elpar)

  !---- Fetch voltage.
  rfv = g_elpar(r_volt)

  !---- Cavity not excited, use drift map.
  if (rfv .eq. zero) then
    call tmdrf(fsec,ftrk,orbit,fmap,el,ek,re,te)
    return
  endif

  !---- Cavity is excited, use full map.

  !---- Initialize.
  call dzero(ek0,6)
  call m66one(rw)
  call dzero(tw,216)

  !---- BV flag
  clight=get_variable('clight ')
  pi=get_variable('pi ')
  twopi=two*pi

  rff = g_elpar(r_freq)
  rfl = g_elpar(r_lag)
  bvk = node_value('other_bv ')
  deltap = get_value('probe ','deltap ')
  pc = get_value('probe ','pc ')

  !-- LD: 20.6.2014 (bvk=-1: not -V -> V but lag -> pi-lag)
  if (bvk .eq. -one) then
    rfl = 0.5-rfl
  endif

  !---- Set up.
  omega = rff * ten6p * twopi / clight
  vrf   = rfv * ten3m / (pc * (one + deltap))
  phirf = rfl * twopi - omega * orbit(5)
  c0 =   vrf * sin(phirf)
  c1 = - vrf * cos(phirf) * omega
  c2 = - vrf * sin(phirf) * omega**2 * half

  !---- Transfer map.
  fmap = .true.
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

  !---- Sandwich cavity between two drifts.
  if (el .ne. zero) then
    dl = el / two
    call tmdrf0(fsec,ftrk,orbit,fmap,dl,ek0,rw,tw)
    call tmcat(fsec,re,te,rw,tw,re,te)
    if(centre_cptk.or.centre_bttk) then
      call dcopy(orbit,orbit00,6)
      call dcopy(ek,ek00,6)
      call dcopy(re,re00,36)
      call dcopy(te,te00,216)
      if(centre_cptk) then
        call dcopy(orbit,orbit0,6)
        call twcptk(re,orbit0)
      endif
      if(centre_bttk) call twbttk(re,te)
      call dcopy(orbit00,orbit,6)
      call dcopy(ek00,ek,6)
      call dcopy(re00,re,36)
      call dcopy(te00,te,216)
    endif
    call tmdrf0(fsec,ftrk,orbit,fmap,dl,ek0,rw,tw)
    call tmcat(fsec,rw,tw,re,te,re,te)
  endif
end SUBROUTINE tmrf

SUBROUTINE tmcat(fsec,rb,tb,ra,ta,rd,td)

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
  logical fsec
  integer i1,i2,i3
  double precision ra(6,6),rb(6,6),rd(6,6),rw(6,6),ta(6,6,6),       &
       tb(36,6),td(6,6,6),ts(36,6),tw(6,6,6)

  !---- Initialization
  call dzero(tw,216)

  !---- Transfer matrix.
  do i2 = 1, 6
     do i1 = 1, 6
        rw(i1,i2) = rb(i1,1) * ra(1,i2) + rb(i1,2) * ra(2,i2)         &
             + rb(i1,3) * ra(3,i2) + rb(i1,4) * ra(4,i2)                       &
             + rb(i1,5) * ra(5,i2) + rb(i1,6) * ra(6,i2)
     enddo
  enddo

  !---- Second order terms.
  if (fsec) then
     do i3 = 1, 6
        do i1 = 1, 36
           ts(i1,i3) = tb(i1,1) * ra(1,i3) + tb(i1,2) * ra(2,i3)       &
                + tb(i1,3) * ra(3,i3) + tb(i1,4) * ra(4,i3)                       &
                + tb(i1,5) * ra(5,i3) + tb(i1,6) * ra(6,i3)
        enddo
     enddo
     do i2 = 1, 6
        do i3 = i2, 6
           do i1 = 1, 6
              tw(i1,i2,i3) =                                            &
                   rb(i1,1) * ta(1,i2,i3) + rb(i1,2) * ta(2,i2,i3)                   &
                   + rb(i1,3) * ta(3,i2,i3) + rb(i1,4) * ta(4,i2,i3)                 &
                   + rb(i1,5) * ta(5,i2,i3) + rb(i1,6) * ta(6,i2,i3)                 &
                   + ts(i1,   i2) * ra(1,i3) + ts(i1+ 6,i2) * ra(2,i3)               &
                   + ts(i1+12,i2) * ra(3,i3) + ts(i1+18,i2) * ra(4,i3)               &
                   + ts(i1+24,i2) * ra(5,i3) + ts(i1+30,i2) * ra(6,i3)
              tw(i1,i3,i2) = tw(i1,i2,i3)
           enddo
        enddo
     enddo
  endif

  !---- Copy result.
  call dcopy(rw,rd,36)
  call dcopy(tw,td,216)

end SUBROUTINE tmcat
SUBROUTINE tmtrak(ek,re,te,orb1,orb2)

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
  integer i,k,l,get_option
  double precision sum1,sum2,ek(6),re(6,6),te(6,6,6),orb1(6),       &
       orb2(6),temp(6),zero
  parameter(zero=0d0)

  !---- initialize.
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
  call dcopy(temp,orb2,6)

  !---- Symplectify transfer matrix.
  if(get_option('sympl ').ne.0) call tmsymp(re)

end SUBROUTINE tmtrak
SUBROUTINE tmsymp(r)

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
  logical eflag
  integer i,j
  double precision a(6,6),b(6,6),r(6,6),v(6,6),zero,one,two,nrm
  parameter(zero=0d0,one=1d0,two=2d0)

  do i = 1, 6
     do j = 1, 6
        a(i,j) = - r(i,j)
        b(i,j) = + r(i,j)
     enddo
     a(i,i) = a(i,i) + one
     b(i,i) = b(i,i) + one
  enddo
  call m66div(a,b,v,eflag)
  if (eflag) goto 100
  call m66inv(v,a)
  do i = 1, 6
     do j = 1, 6
        a(i,j) = (a(i,j) - v(i,j)) / two
        b(i,j) = - a(i,j)
     enddo
     b(i,i) = b(i,i) + one
     a(i,i) = a(i,i) + one
  enddo
  call m66div(a,b,v,eflag)
  if (eflag) goto 100
  call dcopy(v,r,36)
  goto 999

100 continue
  call m66symp(r,nrm)
  if (nrm .gt. zero) then
    print *," Singular matrix occurred during simplectification of R (left unchanged)."
    print *," The column norm of R'*J*R-J is ",nrm
  endif

999 end SUBROUTINE tmsymp
SUBROUTINE tmali1(orb1, errors, betas, gammas, orb2, rm)

  use twiss0fi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for orbit displacement at entry of an element.     *
  !     Input:                                                           *
  !     orb1(6)   (real)    Orbit before misalignment.                   *
  !     errors(align_max) (real)    alignment errors                     *
  !     betas     (real)    current beam beta                            *
  !     gammas    (real)    current beam gamma                           *
  !     Output:                                                          *
  !     orb2(6)   (real)    Orbit after misalignment.                    *
  !     rm(6,6)   (real)    First order transfer matrix w.r.t. orbit.    *
  !----------------------------------------------------------------------*
  double precision ds,dx,dy,phi,psi,rm(6,6),s2,the,betas, gammas,   &
       w(3,3),orb1(6),orb2(6),orbt(6),errors(align_max)

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
  call m66one(rm)
  rm(2,2) = w(1,1)
  rm(2,4) = w(2,1)
  rm(2,6) = w(3,1) / betas
  rm(4,2) = w(1,2)
  rm(4,4) = w(2,2)
  rm(4,6) = w(3,2) / betas

  rm(1,1) =   w(2,2) / w(3,3)
  rm(1,2) = rm(1,1) * s2
  rm(1,3) = - w(1,2) / w(3,3)
  rm(1,4) = rm(1,3) * s2
  rm(3,1) = - w(2,1) / w(3,3)
  rm(3,2) = rm(3,1) * s2
  rm(3,3) =   w(1,1) / w(3,3)
  rm(3,4) = rm(3,3) * s2
  rm(5,1) = w(1,3) / (w(3,3) * betas)
  rm(5,2) = rm(5,1) * s2
  rm(5,3) = w(2,3) / (w(3,3) * betas)
  rm(5,4) = rm(5,3) * s2
  rm(5,6) = - s2 / (betas * gammas)**2

  !---- Track orbit.
  call m66byv(rm, orb1, orbt)
  orb2(1) = orbt(1) - (w(2,2) * dx - w(1,2) * dy) / w(3,3)
  orb2(2) = orbt(2) + w(3,1)
  orb2(3) = orbt(3) - (w(1,1) * dy - w(2,1) * dx) / w(3,3)
  orb2(4) = orbt(4) + w(3,2)
  orb2(5) = orbt(5) - s2 / betas
  orb2(6) = orbt(6)

end SUBROUTINE tmali1
SUBROUTINE tmali2(el, orb1, errors, betas, gammas, orb2, rm)

  use twiss0fi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     TRANSPORT map for orbit displacement at exit of an element.      *
  !     Input:                                                           *
  !     orb1(6)   (real)    Orbit before misalignment.                   *
  !     errors(align_max) (real)    alignment errors                     *
  !     betas     (real)    current beam beta                            *
  !     gammas    (real)    current beam gamma                           *
  !     Output:                                                          *
  !     orb2(6)   (real)    Orbit after misalignment.                    *
  !     rm(6,6)   (real)    First order transfer matrix w.r.t. orbit.    *
  !----------------------------------------------------------------------*
  double precision ds,dx,dy,el,phi,psi,s2,the,tilt,orb1(6),orb2(6), &
       rm(6,6),v(3),ve(3),w(3,3),we(3,3),orbt(6),errors(align_max),betas,&
       gammas,zero
  parameter(zero=0d0)

  !---- Misalignment rotation matrix w.r.t. entrance system.
  dx  = errors(1)
  dy  = errors(2)
  ds  = errors(3)
  the = errors(5)
  phi = errors(4)
  psi = errors(6)
  tilt=zero
  call sumtrx(the, phi, psi, w)
  !---- VE and WE represent the change of reference.
  call suelem(el, ve, we, tilt)
  !---- Misalignment displacements at exit w.r.t. entrance system.
  v(1) = dx + w(1,1)*ve(1)+w(1,2)*ve(2)+w(1,3)*ve(3)-ve(1)
  v(2) = dy + w(2,1)*ve(1)+w(2,2)*ve(2)+w(2,3)*ve(3)-ve(2)
  v(3) = ds + w(3,1)*ve(1)+w(3,2)*ve(2)+w(3,3)*ve(3)-ve(3)

  !---- Convert all references to exit, build additional drift.
  call sutran(w, v, we)
  s2 = - (w(1,3) * v(1) + w(2,3) * v(2) + w(3,3) * v(3)) / w(3,3)

  !---- Transfer matrix.
  call m66one(rm)
  rm(1,1) = w(1,1)
  rm(3,1) = w(2,1)
  rm(5,1) = w(3,1) / betas
  rm(1,3) = w(1,2)
  rm(3,3) = w(2,2)
  rm(5,3) = w(3,2) / betas

  rm(2,2) =   w(2,2) / w(3,3)
  rm(1,2) = rm(2,2) * s2
  rm(4,2) = - w(1,2) / w(3,3)
  rm(3,2) = rm(4,2) * s2
  rm(2,4) = - w(2,1) / w(3,3)
  rm(1,4) = rm(2,4) * s2
  rm(4,4) =   w(1,1) / w(3,3)
  rm(3,4) = rm(4,4) * s2
  rm(2,6) = w(1,3) / (w(3,3) * betas)
  rm(1,6) = rm(2,6) * s2
  rm(4,6) = w(2,3) / (w(3,3) * betas)
  rm(3,6) = rm(4,6) * s2
  rm(5,6) = - s2 / (betas * gammas)**2

  !---- Track orbit.
  orbt(1) = orb1(1) + (w(2,2) * v(1) - w(1,2) * v(2)) / w(3,3)
  orbt(2) = orb1(2) - w(3,1)
  orbt(3) = orb1(3) + (w(1,1) * v(2) - w(2,1) * v(1)) / w(3,3)
  orbt(4) = orb1(4) - w(3,2)
  orbt(5) = orb1(5) - s2 / betas
  orbt(6) = orb1(6)
  call m66byv(rm, orbt, orb2)

end SUBROUTINE tmali2
SUBROUTINE tmbb(fsec,ftrk,orbit,fmap,re,te)

  use trackfi
  implicit none

  !----------------------------------------------------------------------*
  !     purpose:                                                         *
  !     transport map for beam-beam element.                             *
  !     input:                                                           *
  !     fsec      (logical) must be .true. for this purpose              *
  !     ftrk      (logical) must be true for this purpose                *
  !     input/output:                                                    *
  !     orbit(6)  (double)  closed orbit (only kick is added             *
  !     since BB thin)                                                   *
  !     output:                                                          *
  !     fmap      (logical) if map has been calculated correctly         *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !----------------------------------------------------------------------*
  integer beamshape, b_dir_int,get_option
  logical fsec,ftrk,fmap,first,bb_ultra_relati
  double precision orbit(6),re(6,6),te(6,6,6),node_value,get_value
  double precision fk,parvec(26),q,q_prime,dp,get_variable
  double precision gamma0,beta0,beta_dp,ptot,b_dir
  double precision zero,one,two
  parameter(zero=0.0d0,one=1.0d0,two=2.0d0)
  save first
  data first / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---  standard 4D
  q=get_value('probe ','charge ')
  q_prime=node_value('charge ')
  parvec(5)=get_value('probe ', 'arad ')
  parvec(6)=node_value('charge ') * get_value('probe ', 'npart ')
  parvec(7)=get_value('probe ','gamma ')
  !---- Calculate momentum deviation and according changes
  !     of the relativistic factor beta0
  dp  = get_variable('track_deltap ')
  gamma0 = parvec(7)
  beta0 = sqrt(one-one/gamma0**2)
  ptot = beta0*gamma0*(one+dp)
  beta_dp = ptot / sqrt(one + ptot**2)
  b_dir_int = node_value('bbdir ')
  b_dir = dble(b_dir_int)
  b_dir = b_dir/sqrt(b_dir*b_dir + 1.0d-32)
  !
  bb_ultra_relati=get_option('bb_ultra_relati ').ne.0
  if(bb_ultra_relati) then
     fk = two * parvec(5) * parvec(6) / parvec(7)
  else
     fk = two*parvec(5)*parvec(6)/parvec(7)/beta0/(one+dp)/q*          &
       (one-beta0*beta_dp*b_dir)/(beta_dp+0.5*(b_dir-one)*b_dir*beta0)
  endif
  !---- chose beamshape
  beamshape = node_value('bbshape ')
  if(beamshape.lt.1.or.beamshape.gt.3) then
     beamshape=1
     if(first) then
        first = .false.
        call aawarn('TMBB: ',                                         &
             'beamshape out of range, set to default=1')
     endif
  endif
  if(beamshape.eq.1) call tmbb_gauss(fsec,ftrk,orbit,fmap,re,te,fk)
  if(beamshape.eq.2) call tmbb_flattop(fsec,ftrk,orbit,fmap,re,te,  &
       fk)
  if(beamshape.eq.3) call tmbb_hollowparabolic(fsec,ftrk,orbit,fmap,&
       re,te,fk)
end SUBROUTINE tmbb
SUBROUTINE tmbb_gauss(fsec,ftrk,orbit,fmap,re,te,fk)

  use bbfi
  use twisslfi
  use spch_bbfi
  use fasterror
  implicit none

  logical fsec,ftrk,fmap,bborbit,bb_sxy_update
  integer get_option,mylen
  double precision orbit(6),re(6,6),te(6,6,6),pi,sx,sy,  &
       xm,ym,sx2,sy2,xs,ys,rho2,fk,tk,exk,phix,phiy,rho4,phixx,phixy,    &
       phiyy,rho6,rk,exkc,xb,yb,phixxx,phixxy,phixyy,phiyyy,crx,cry,xr,  &
       yr,r,r2,cbx,cby,get_variable,node_value,zero,one,two,   &
       orbit0(6),three,ten3m,explim
  double precision orbit00(6),re00(6,6),te00(6,6,6)
  character*20 text
  character*(name_len) name
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
       explim=150d0)
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  if (bbd_flag .ne. 0 .and. .not. bborbit)  then
     if (bbd_cnt .eq. bbd_max)  then
        call aawarn('TMBB_GAUSS: ','maximum bb number reached')
     else
        bbd_cnt = bbd_cnt + 1
        bbd_loc(bbd_cnt) = bbd_pos
        bb_kick(1,bbd_cnt) = zero
        bb_kick(2,bbd_cnt) = zero
     endif
  endif
  pi=get_variable('pi ')
  fmap = .true.
  bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
  if(bb_sxy_update) then
     fk = fk * rat_bb_n_ions !Ratio_for_bb_N_ions
     name=' '
     call element_name(name,len(name))
     mylen=len_trim(name)
     i_spch=i_spch+1
     
     if(i_spch.gt.N_spch) then
        write(text, '(1p,i8)') i_spch
        call aafail('TMBB: ', 'Table with too few BB elements: '//  &
             text)
     endif
     if(spch_bb_name(i_spch)(:mylen).ne.name(:mylen)) then
        call aafail('TMBB: ', 'wrong element name in Table: '//     &
             'spch_bb')
     endif

     sx=sqrt(betx_bb(i_spch)*Ex_rms+(dx_bb(i_spch)*sigma_p)**2)
     sy=sqrt(bety_bb(i_spch)*Ey_rms+(dy_bb(i_spch)*sigma_p)**2)
  else
     sx = node_value('sigx ')
     sy = node_value('sigy ')
  endif
  xm = node_value('xma ')
  ym = node_value('yma ')
  if (fk .ne. zero)  then
     !---- if tracking is desired ...
     if (ftrk) then
        sx2 = sx * sx
        sy2 = sy * sy
        xs  = orbit(1) - xm
        ys  = orbit(3) - ym

        !---- limit formulas for sigma(x) = sigma(y).
        if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
           rho2 = xs * xs + ys * ys

           !---- limit case for xs = ys = 0.
           if (rho2 .eq. zero) then
              re(2,1) = fk / (two * sx2)
              re(4,3) = fk / (two * sx2)

              !---- general case.
           else
              tk = rho2 / (two * sx2)
              if (tk .gt. explim) then
                 exk  = zero
                 exkc = one
                 phix = xs * fk / rho2
                 phiy = ys * fk / rho2
              else
                 exk  = exp(-tk)
                 exkc = one - exk
                 phix = xs * fk / rho2 * exkc
                 phiy = ys * fk / rho2 * exkc
              endif

              !---- orbit kick - only applied if option bborbit (HG 5/12/01),
              !     else stored
              if (bborbit) then
                 orbit(2) = orbit(2) + phix
                 orbit(4) = orbit(4) + phiy
              elseif (bbd_flag .ne. 0)  then
                 bb_kick(1,bbd_cnt) = phix
                 bb_kick(2,bbd_cnt) = phiy
              endif
              !---- first-order effects.
              rho4 = rho2 * rho2
              phixx = fk * (- exkc * (xs*xs - ys*ys) / rho4             &
                   + exk * xs*xs / (rho2 * sx2))
              phixy = fk * (- exkc * two * xs * ys / rho4               &
                   + exk * xs*ys / (rho2 * sx2))
              phiyy = fk * (+ exkc * (xs*xs - ys*ys) / rho4             &
                   + exk * ys*ys / (rho2 * sx2))
              re(2,1) = phixx
              re(2,3) = phixy
              re(4,1) = phixy
              re(4,3) = phiyy

              !---- second-order effects.
              if (fsec) then
                 rho6 = rho4 * rho2
                 phixxx = fk*xs * (+ exkc * (xs*xs - three*ys*ys) / rho6 &
                      - exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2)                &
                      - exk * xs*xs / (two * rho2 * sx2**2))
                 phixxy = fk*ys * (+ exkc * (three*xs*xs - ys*ys) / rho6 &
                      - exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2)                &
                      - exk * xs*xs / (two * rho2 * sx2**2))
                 phixyy = fk*xs * (- exkc * (xs*xs - three*ys*ys) / rho6 &
                      + exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2)                &
                      - exk * ys*ys / (two * rho2 * sx2**2))
                 phiyyy = fk*ys * (- exkc * (three*xs*xs - ys*ys) / rho6 &
                      + exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2)                &
                      - exk * ys*ys / (two * rho2 * sx2**2))
                 te(2,1,1) = phixxx
                 te(2,1,3) = phixxy
                 te(2,3,1) = phixxy
                 te(4,1,1) = phixxy
                 te(2,3,3) = phixyy
                 te(4,1,3) = phixyy
                 te(4,3,1) = phixyy
                 te(4,3,3) = phiyyy
              endif
           endif

           !---- case sigma(x) > sigma(y).
        else
           r2 = two * (sx2 - sy2)
           if (sx2 .gt. sy2) then
              r  = sqrt(r2)
              rk = fk * sqrt(pi) / r
              xr = abs(xs) / r
              yr = abs(ys) / r

              if(fasterror_on) then
                 call wzsub(xr, yr, crx, cry)
              else
                 call ccperrf(xr, yr, crx, cry)
              endif
           
              tk = (xs * xs / sx2 + ys * ys / sy2) / two
              if (tk .gt. explim) then
                 exk = zero
                 cbx = zero
                 cby = zero
              else
                 exk = exp(-tk)
                 xb  = (sy / sx) * xr
                 yb  = (sx / sy) * yr

                 if(fasterror_on) then
                    call wzsub(xb, yb, cbx, cby)
                 else
                    call ccperrf(xb, yb, cbx, cby)
                 endif

              endif

              !---- case sigma(x) < sigma(y).
           else
              r  = sqrt(-r2)
              rk = fk * sqrt(pi) / r
              xr = abs(xs) / r
              yr = abs(ys) / r

              if(fasterror_on) then
                 call wzsub(yr, xr, cry, crx)
              else
                 call ccperrf(yr, xr, cry, crx)
              endif
           
              tk = (xs * xs / sx2 + ys * ys / sy2) / two
              if (tk .gt. explim) then
                 exk = zero
                 cbx = zero
                 cby = zero
              else
                 exk = exp(-tk)
                 xb  = (sy / sx) * xr
                 yb  = (sx / sy) * yr

                 if(fasterror_on) then
                    call wzsub(yb, xb, cby, cbx)
                 else
                    call ccperrf(yb, xb, cby, cbx)
                 endif

              endif
           endif

           phix = rk * (cry - exk * cby) * sign(one, xs)
           phiy = rk * (crx - exk * cbx) * sign(one, ys)
           !---- orbit kick - only applied if option bborbit (HG 5/12/01),
           !     else stored
           if (bborbit) then
              orbit(2) = orbit(2) + phix
              orbit(4) = orbit(4) + phiy
           elseif (bbd_flag .ne. 0)  then
              bb_kick(1,bbd_cnt) = phix
              bb_kick(2,bbd_cnt) = phiy
           endif

           !---- first-order effects.
           phixx = (two / r2) * (- (xs * phix + ys * phiy)             &
                + fk * (one - (sy / sx) * exk))
           phixy = (two / r2) * (- (xs * phiy - ys * phix))
           phiyy = (two / r2) * (+ (xs * phix + ys * phiy)             &
                - fk * (one - (sx / sy) * exk))
           re(2,1) = phixx
           re(2,3) = phixy
           re(4,1) = phixy
           re(4,3) = phiyy

           !---- second-order effects.
           if (fsec) then
              phixxx = (- phix - (xs * phixx + ys * phixy)              &
                   + fk * xs * sy * exk / sx**3) / r2
              phixxy = (- phiy - (xs * phixy - ys * phixx)) / r2
              phixyy = (+ phix - (xs * phiyy - ys * phixy)) / r2
              phiyyy = (+ phiy + (xs * phixy + ys * phiyy)              &
                   - fk * ys * sx * exk / sy**3) / r2
              te(2,1,1) = phixxx
              te(2,1,3) = phixxy
              te(2,3,1) = phixxy
              te(4,1,1) = phixxy
              te(2,3,3) = phixyy
              te(4,1,3) = phixyy
              te(4,3,1) = phixyy
              te(4,3,3) = phiyyy
           endif
        endif

        !---- no tracking desired.
     else
        re(2,1) = fk / (sx * (sx + sy))
        re(4,3) = fk / (sy * (sx + sy))
     endif
     !---- centre option
     if(centre_cptk.or.centre_bttk) then
        call dcopy(orbit,orbit00,6)
        call dcopy(re,re00,36)
        call dcopy(te,te00,216)
        if(centre_cptk) then
           call dcopy(orbit,orbit0,6)
           call twcptk(re,orbit0)
        endif
        if(centre_bttk) call twbttk(re,te)
        call dcopy(orbit00,orbit,6)
        call dcopy(re00,re,36)
        call dcopy(te00,te,216)
     endif
  endif

end SUBROUTINE tmbb_gauss
SUBROUTINE tmbb_flattop(fsec,ftrk,orbit,fmap,re,te,fk)

  use bbfi
  use twisslfi
  implicit none

  logical fsec,ftrk,fmap,bborbit,firstflag
  integer get_option
  double precision orbit(6),re(6,6),te(6,6,6),pi,r0x,r0y,&
       xm,ym,r0x2,r0y2,xs,ys,rho2,fk,phix,phiy,phixx,phixy,phiyy,phixxx, &
       phixxy,phixyy,phiyyy,get_variable,node_value,zero,one,  &
       two,three,ten3m,explim,wi,rho,wx,wy,norm,phir,phirr,phirrr,zz
  double precision orbit0(6),orbit00(6),re00(6,6),te00(6,6,6)
  parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,ten3m=1.0d-3,&
       explim=150.0d0)
  save firstflag
  data firstflag / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  phix=zero
  phiy=zero
  bborbit = get_option('bborbit ') .ne. 0
  if (bbd_flag .ne. 0 .and. .not. bborbit)  then
     if (bbd_cnt .eq. bbd_max)  then
        call aawarn('TMBB_FLATTOP: ','maximum bb number reached')
     else
        bbd_cnt = bbd_cnt + 1
        bbd_loc(bbd_cnt) = bbd_pos
        bb_kick(1,bbd_cnt) = zero
        bb_kick(2,bbd_cnt) = zero
     endif
  endif
  pi=get_variable('pi ')
  fmap = .true.
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  wi = node_value('width ')
  xm = node_value('xma ')
  ym = node_value('yma ')
  if (fk .ne. zero)  then
     r0x2 = r0x * r0x
     r0y2 = r0y * r0y
     wx  = r0x * wi
     wy  = r0y * wi
     xs  = orbit(1) - xm
     ys  = orbit(3) - ym
     !---- limit formulas for sigma(x) = sigma(y).
     !-----preliminary the only case considered.
     if (dabs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
        zz = 0.5*(r0x + r0y)
        r0x=zz
        r0y=zz
        r0x2=r0x*r0x
        r0y2=r0y*r0y
        if(firstflag) then
           firstflag=.false.
           call aawarn('TMBB_FLATTOP: ',                               &
                'beam is assumed to be circular')
        endif
     endif
     norm = (12.0*r0x**2+wx**2)/24.0
     !---- if tracking is desired ...
     if (ftrk) then
        rho2 = xs * xs + ys * ys
        rho  = sqrt(rho2)
        if(rho.le.r0x-wx/2.0) then
           phir=0.5/norm
           phix=phir*xs
           phiy=phir*ys
           phixx=0.5/norm
           phixy=zero
           phiyy=0.5/norm
           re(2,1) = phixx*fk
           re(2,3) = phixy*fk
           re(4,1) = phixy*fk
           re(4,3) = phiyy*fk
           if(fsec) then
              te(2,1,1) = zero
              te(2,1,3) = zero
              te(2,3,1) = zero
              te(4,1,1) = zero
              te(2,3,3) = zero
              te(4,1,3) = zero
              te(4,3,1) = zero
              te(4,3,3) = zero
           endif
        else if(rho.gt.r0x-wx/2.0.and.rho.lt.r0x+wx/2.0) then
           phir = ((r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0 +          &
                wx**2/48.0)/rho2 + 0.25 + 0.5*r0x/wx -                            &
                rho/3.0/wx)/norm
           phix = phir*xs
           phiy = phir*ys
           phirr = - (2.0*(r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0 +   &
                wx**2/48.0)/rho2**2 + 1.0/rho/3.0/wx)/norm
           phixx = phir + xs*xs*phirr
           phixy = xs*ys*phirr
           phiyy = phir + ys*ys*phirr
           re(2,1) = phixx*fk
           re(2,3) = phixy*fk
           re(4,1) = phixy*fk
           re(4,3) = phiyy*fk
           !     2nd order effects
           if(fsec) then
              phirrr =  (8.0*(r0x**2/4.0 - r0x**3/6.0/wx - r0x*wx/8.0   &
                   + wx**2/48.0)/rho2**3 + 1.0/3.0/wx/rho**3)/norm
              phixxx = 3.0*xs*phirr + xs**3*phirrr
              phixxy = ys*phirr + xs*xs*ys*phirrr
              phixyy = xs*phirr + xs*ys*ys*phirrr
              phiyyy = 3.0*ys*phirr + ys**3*phirrr
              te(2,1,1) = phixxx*fk
              te(2,1,3) = phixxy*fk
              te(2,3,1) = phixxy*fk
              te(4,1,1) = phixxy*fk
              te(2,3,3) = phixyy*fk
              te(4,1,3) = phixyy*fk
              te(4,3,1) = phixyy*fk
              te(4,3,3) = phiyyy*fk
           endif
        else if(rho.ge.r0x+wx/2.0) then
           phir = 1.0/rho2
           phix = xs*phir
           phiy = ys*phir
           phirr = -2.0/rho2**2
           phixx = phir+xs*xs*phirr
           phixy = xs*ys*phirr
           phiyy = phir+ys*ys*phirr
           re(2,1) = phixx*fk
           re(2,3) = phixy*fk
           re(4,1) = phixy*fk
           re(4,3) = phiyy*fk
           !     2nd order effects
           if(fsec) then
              phirrr = 8.0/rho2**3
              phixxx = 3.0*xs*phirr + xs**3*phirrr
              phixxy = ys*phirr + xs*xs*ys*phirrr
              phixyy = xs*phirr + xs*ys*ys*phirrr
              phiyyy = 3.0*ys*phirr + ys**3*phirrr
              te(2,1,1) = phixxx*fk
              te(2,1,3) = phixxy*fk
              te(2,3,1) = phixxy*fk
              te(4,1,1) = phixxy*fk
              te(2,3,3) = phixyy*fk
              te(4,1,3) = phixyy*fk
              te(4,3,1) = phixyy*fk
              te(4,3,3) = phiyyy*fk
           end if
        endif
        !---- orbit kick - only applied if option bborbit (HG 5/12/01),
        !     else stored
        if (bborbit) then
           orbit(2) = orbit(2) + phix*fk
           orbit(4) = orbit(4) + phiy*fk
        elseif (bbd_flag .ne. 0)  then
           bb_kick(1,bbd_cnt) = phix*fk
           bb_kick(2,bbd_cnt) = phiy*fk
        endif
        !---- no tracking desired.
     else
        phixx=0.5/norm
        phiyy=0.5/norm
        re(2,1) = phixx*fk
        re(4,3) = phiyy*fk
     endif
     !---- centre option
     if(centre_cptk.or.centre_bttk) then
        call dcopy(orbit,orbit00,6)
        call dcopy(re,re00,36)
        call dcopy(te,te00,216)
        if(centre_cptk) then
           call dcopy(orbit,orbit0,6)
           call twcptk(re,orbit0)
        endif
        if(centre_bttk) call twbttk(re,te)
        call dcopy(orbit00,orbit,6)
        call dcopy(re00,re,36)
        call dcopy(te00,te,216)
     endif
  endif
  !
end SUBROUTINE tmbb_flattop
SUBROUTINE tmbb_hollowparabolic(fsec,ftrk,orbit,fmap,re,te,fk)

  use bbfi
  use twisslfi
  implicit none

  logical fsec,ftrk,fmap,bborbit,firstflag
  integer get_option
  double precision orbit(6),re(6,6),te(6,6,6),pi,r0x,r0y,&
       xm,ym,r0x2,r0y2,xs,ys,rho2,fk,phix,phiy,phixx,phixy,phiyy,phixxx, &
       phixxy,phixyy,phiyyy,get_variable,node_value,zero,one,  &
       two,orbit0(6),three,ten3m,explim,wi,rho,wx,wy,phir,phirr,phirrr,zz
  double precision orbit00(6),re00(6,6),te00(6,6,6)
  parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,ten3m=1.0d-3,&
       explim=150.0d0)
  save firstflag
  data firstflag / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  phix=zero
  phiy=zero
  bborbit = get_option('bborbit ') .ne. 0
  if (bbd_flag .ne. 0 .and. .not. bborbit)  then
     if (bbd_cnt .eq. bbd_max)  then
        call aawarn('TMBB_HOLLOWPARABOLIC: ',                         &
             'maximum bb number reached')
     else
        bbd_cnt = bbd_cnt + 1
        bbd_loc(bbd_cnt) = bbd_pos
        bb_kick(1,bbd_cnt) = zero
        bb_kick(2,bbd_cnt) = zero
     endif
  endif
  pi=get_variable('pi ')
  fmap = .true.
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  !     width is given as FWHM of the parabolic density profile,
  !     but formulas were derived with half width at the bottom of this
  !     density profile
  wi = node_value('width ')
  wi = wi/sqrt(2.0)
  xm = node_value('xma ')
  ym = node_value('yma ')
  if (fk .ne. zero)  then
     !---- if tracking is desired ...
     if (ftrk) then
        r0x2 = r0x * r0x
        r0y2 = r0y * r0y
        wx  = r0x * wi
        wy  = r0y * wi
        xs  = orbit(1) - xm
        ys  = orbit(3) - ym
        !---- limit formulas for sigma(x) = sigma(y).
        !-----preliminary the only case considered.
        if (abs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
           zz = 0.5*(r0x + r0y)
           r0x=zz
           r0y=zz
           r0x2=r0x*r0x
           r0y2=r0y*r0y
           if(firstflag) then
              firstflag=.false.
              call aawarn('TMBB_HOLLOWPARABOLIC: ',                     &
                   'beam is assumed to be circular')
           endif
        endif
        rho2 = xs * xs + ys * ys
        rho  = sqrt(rho2)
        if(rho.le.r0x-wx) then
           re(2,1) = zero
           re(4,3) = zero
           re(2,3) = zero
           re(4,1) = zero
           if(fsec) then
              te(2,1,1) = zero
              te(2,1,3) = zero
              te(2,3,1) = zero
              te(4,1,1) = zero
              te(2,3,3) = zero
              te(4,1,3) = zero
              te(4,3,1) = zero
              te(4,3,3) = zero
           end if
        else if(rho.gt.r0x-wx.and.rho.lt.r0x+wx) then
           phir=0.75/wx/r0x/rho2*(r0x**4/12.0/wx**2 - r0x**2/2.0+      &
                2.0*r0x*wx/3.0 - wx**2/4.0 + rho2/2.0*(1.0 -                      &
                r0x**2/wx**2) + rho**3/3.0*2.0*r0x/wx**2 -                        &
                rho**4/4.0/wx**2)
           phix = phir*xs
           phiy = phir*ys
           phirr = 0.75/wx/r0x*(-2.0/rho**4*(r0x**4/12.0/wx**2-        &
                r0x**2/2.0+2.0*r0x*wx/3.0-wx**2/4.0)+                             &
                2.0*r0x/wx**2/3.0/rho-0.5/wx**2)
           phixx = phir+xs*xs*phirr
           phixy = xs*ys*phirr
           phiyy = phir+ys*ys*phirr
           re(2,1) = phixx*fk
           re(2,3) = phixy*fk
           re(4,1) = phixy*fk
           re(4,3) = phiyy*fk
           ! 2nd order effects
           if(fsec) then
              phirrr = 0.75/wx/r0x*(8.0/rho2**3*(r0x**4/12.0/           &
                   wx**2-r0x**2/2.0+2.0*r0x*wx/3.0-wx**2/4.0) -                      &
                   2.0*r0x/3.0/wx**2/rho**3)
              phixxx = 3.0*xs*phirr + xs**3*phirrr
              phixxy = ys*phirr + xs*xs*ys*phirrr
              phixyy = xs*phirr + xs*ys*ys*phirrr
              phiyyy = 3.0*ys*phirr + ys**3*phirrr
              te(2,1,1) = phixxx*fk
              te(2,1,3) = phixxy*fk
              te(2,3,1) = phixxy*fk
              te(4,1,1) = phixxy*fk
              te(2,3,3) = phixyy*fk
              te(4,1,3) = phixyy*fk
              te(4,3,1) = phixyy*fk
              te(4,3,3) = phiyyy*fk
           end if
        else if(rho.ge.r0x+wx) then
           phir = 1.0/rho2
           phix = xs*phir
           phiy = ys*phir
           phirr = -2.0/rho2**2
           phixx = phir+xs*xs*phirr
           phixy = xs*ys*phirr
           phiyy = phir+ys*ys*phirr
           re(2,1) = phixx*fk
           re(2,3) = phixy*fk
           re(4,1) = phixy*fk
           re(4,3) = phiyy*fk
           !     2nd order effects
           if(fsec) then
              phirrr = 8.0/rho2**3
              phixxx = 3.0*xs*phirr + xs**3*phirrr
              phixxy = ys*phirr + xs*xs*ys*phirrr
              phixyy = xs*phirr + xs*ys*ys*phirrr
              phiyyy = 3.0*ys*phirr + ys**3*phirrr
              te(2,1,1) = phixxx*fk
              te(2,1,3) = phixxy*fk
              te(2,3,1) = phixxy*fk
              te(4,1,1) = phixxy*fk
              te(2,3,3) = phixyy*fk
              te(4,1,3) = phixyy*fk
              te(4,3,1) = phixyy*fk
              te(4,3,3) = phiyyy*fk
           end if
        endif
        !---- orbit kick - only applied if option bborbit (HG 5/12/01),
        !     else stored
        if (bborbit) then
           orbit(2) = orbit(2) + phix*fk
           orbit(4) = orbit(4) + phiy*fk
        elseif (bbd_flag .ne. 0)  then
           bb_kick(1,bbd_cnt) = phix*fk
           bb_kick(2,bbd_cnt) = phiy*fk
        endif
        !---- no tracking desired.
     else
        re(2,1) = zero
        re(4,3) = zero
     endif
     !---- centre option
     if(centre_cptk.or.centre_bttk) then
        call dcopy(orbit,orbit00,6)
        call dcopy(re,re00,36)
        call dcopy(te,te00,216)
        if(centre_cptk) then
           call dcopy(orbit,orbit0,6)
           call twcptk(re,orbit0)
        endif
        if(centre_bttk) call twbttk(re,te)
        call dcopy(orbit00,orbit,6)
        call dcopy(re00,re,36)
        call dcopy(te00,te,216)
     endif
  endif
  !
end SUBROUTINE tmbb_hollowparabolic
!!$SUBROUTINE tmbb_old(fsec,ftrk,orbit,fmap,re,te)
!!$
!!$  use bbfi
!!$  use twisslfi
!!$  implicit none
!!$
!!$  !----------------------------------------------------------------------*
!!$  !     purpose:                                                         *
!!$  !     transport map for beam-beam element.                             *
!!$  !     input:                                                           *
!!$  !     fsec      (logical) must be .true. for this purpose              *
!!$  !     ftrk      (logical) must be true for this purpose                *
!!$  !     input/output:                                                    *
!!$  !     orbit(6)  (double)  closed orbit (only kick is added             *
!!$  !     since BB thin)                                                   *
!!$  !     output:                                                          *
!!$  !     fmap      (logical) if map has been calculated correctly         *
!!$  !     re(6,6)   (double)  transfer matrix.                             *
!!$  !     te(6,6,6) (double)  second-order terms.                          *
!!$  !----------------------------------------------------------------------*
!!$  logical fsec,ftrk,fmap,bborbit
!!$  integer get_option
!!$  double precision parvec(26),orbit(6),re(6,6),te(6,6,6),pi,sx,sy,  &
!!$       xm,ym,sx2,sy2,xs,ys,rho2,fk,tk,exk,phix,phiy,rho4,phixx,phixy,    &
!!$       phiyy,rho6,rk,exkc,xb,yb,phixxx,phixxy,phixyy,phiyyy,crx,cry,xr,  &
!!$       yr,r,r2,cbx,cby,get_variable,get_value,node_value,zero,one,two,   &
!!$       orbit0(6),three,ten3m,explim
!!$  double precision orbit00(6),re00(6,6),te00(6,6,6)
!!$  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
!!$       explim=150d0)
!!$  !     if x > explim, exp(-x) is outside machine limits.
!!$
!!$  !---- initialize.
!!$  bborbit = get_option('bborbit ') .ne. 0
!!$  if (bbd_flag .ne. 0 .and. .not. bborbit)  then
!!$     if (bbd_cnt .eq. bbd_max)  then
!!$        call aawarn('TMBB: ','maximum bb number reached')
!!$     else
!!$        bbd_cnt = bbd_cnt + 1
!!$        bbd_loc(bbd_cnt) = bbd_pos
!!$        bb_kick(1,bbd_cnt) = zero
!!$        bb_kick(2,bbd_cnt) = zero
!!$     endif
!!$  endif
!!$  pi=get_variable('pi ')
!!$  fmap = .true.
!!$  sx = node_value('sigx ')
!!$  sy = node_value('sigy ')
!!$  xm = node_value('xma ')
!!$  ym = node_value('yma ')
!!$  !---  standard 4D
!!$  parvec(5)=get_value('probe ', 'arad ')
!!$  parvec(6)=node_value('charge ') * get_value('probe ', 'npart ')
!!$  parvec(7)=get_value('probe ','gamma ')
!!$  fk = two * parvec(5) * parvec(6) / parvec(7)
!!$  if (fk .ne. zero)  then
!!$     !---- if tracking is desired ...
!!$     if (ftrk) then
!!$        sx2 = sx * sx
!!$        sy2 = sy * sy
!!$        xs  = orbit(1) - xm
!!$        ys  = orbit(3) - ym
!!$
!!$        !---- limit formulas for sigma(x) = sigma(y).
!!$        if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
!!$           rho2 = xs * xs + ys * ys
!!$
!!$           !---- limit case for xs = ys = 0.
!!$           if (rho2 .eq. zero) then
!!$              re(2,1) = fk / (two * sx2)
!!$              re(4,3) = fk / (two * sx2)
!!$
!!$              !---- general case.
!!$           else
!!$              tk = rho2 / (two * sx2)
!!$              if (tk .gt. explim) then
!!$                 exk  = zero
!!$                 exkc = one
!!$                 phix = xs * fk / rho2
!!$                 phiy = ys * fk / rho2
!!$              else
!!$                 exk  = exp(-tk)
!!$                 exkc = one - exk
!!$                 phix = xs * fk / rho2 * exkc
!!$                 phiy = ys * fk / rho2 * exkc
!!$              endif
!!$
!!$              !---- orbit kick - only applied if option bborbit (HG 5/12/01),
!!$              !     else stored
!!$              if (bborbit) then
!!$                 orbit(2) = orbit(2) + phix
!!$                 orbit(4) = orbit(4) + phiy
!!$              elseif (bbd_flag .ne. 0)  then
!!$                 bb_kick(1,bbd_cnt) = phix
!!$                 bb_kick(2,bbd_cnt) = phiy
!!$              endif
!!$              !---- first-order effects.
!!$              rho4 = rho2 * rho2
!!$              phixx = fk * (- exkc * (xs*xs - ys*ys) / rho4             &
!!$                   + exk * xs*xs / (rho2 * sx2))
!!$              phixy = fk * (- exkc * two * xs * ys / rho4               &
!!$                   + exk * xs*ys / (rho2 * sx2))
!!$              phiyy = fk * (+ exkc * (xs*xs - ys*ys) / rho4             &
!!$                   + exk * ys*ys / (rho2 * sx2))
!!$              re(2,1) = phixx
!!$              re(2,3) = phixy
!!$              re(4,1) = phixy
!!$              re(4,3) = phiyy
!!$
!!$              !---- second-order effects.
!!$              if (fsec) then
!!$                 rho6 = rho4 * rho2
!!$                 phixxx = fk*xs * (+ exkc * (xs*xs - three*ys*ys) / rho6 &
!!$                      - exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2)                &
!!$                      - exk * xs*xs / (two * rho2 * sx2**2))
!!$                 phixxy = fk*ys * (+ exkc * (three*xs*xs - ys*ys) / rho6 &
!!$                      - exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2)                &
!!$                      - exk * xs*xs / (two * rho2 * sx2**2))
!!$                 phixyy = fk*xs * (- exkc * (xs*xs - three*ys*ys) / rho6 &
!!$                      + exk * (xs*xs - three*ys*ys) / (two * rho4 * sx2)                &
!!$                      - exk * ys*ys / (two * rho2 * sx2**2))
!!$                 phiyyy = fk*ys * (- exkc * (three*xs*xs - ys*ys) / rho6 &
!!$                      + exk * (three*xs*xs - ys*ys) / (two * rho4 * sx2)                &
!!$                      - exk * ys*ys / (two * rho2 * sx2**2))
!!$                 te(2,1,1) = phixxx
!!$                 te(2,1,3) = phixxy
!!$                 te(2,3,1) = phixxy
!!$                 te(4,1,1) = phixxy
!!$                 te(2,3,3) = phixyy
!!$                 te(4,1,3) = phixyy
!!$                 te(4,3,1) = phixyy
!!$                 te(4,3,3) = phiyyy
!!$              endif
!!$           endif
!!$
!!$           !---- case sigma(x) > sigma(y).
!!$        else
!!$           r2 = two * (sx2 - sy2)
!!$           if (sx2 .gt. sy2) then
!!$              r  = sqrt(r2)
!!$              rk = fk * sqrt(pi) / r
!!$              xr = abs(xs) / r
!!$              yr = abs(ys) / r
!!$              call ccperrf(xr, yr, crx, cry)
!!$              tk = (xs * xs / sx2 + ys * ys / sy2) / two
!!$              if (tk .gt. explim) then
!!$                 exk = zero
!!$                 cbx = zero
!!$                 cby = zero
!!$              else
!!$                 exk = exp(-tk)
!!$                 xb  = (sy / sx) * xr
!!$                 yb  = (sx / sy) * yr
!!$                 call ccperrf(xb, yb, cbx, cby)
!!$              endif
!!$
!!$              !---- case sigma(x) < sigma(y).
!!$           else
!!$              r  = sqrt(-r2)
!!$              rk = fk * sqrt(pi) / r
!!$              xr = abs(xs) / r
!!$              yr = abs(ys) / r
!!$              call ccperrf(yr, xr, cry, crx)
!!$              tk = (xs * xs / sx2 + ys * ys / sy2) / two
!!$              if (tk .gt. explim) then
!!$                 exk = zero
!!$                 cbx = zero
!!$                 cby = zero
!!$              else
!!$                 exk = exp(-tk)
!!$                 xb  = (sy / sx) * xr
!!$                 yb  = (sx / sy) * yr
!!$                 call ccperrf(yb, xb, cby, cbx)
!!$              endif
!!$           endif
!!$
!!$           phix = rk * (cry - exk * cby) * sign(one, xs)
!!$           phiy = rk * (crx - exk * cbx) * sign(one, ys)
!!$           !---- orbit kick - only applied if option bborbit (HG 5/12/01),
!!$           !     else stored
!!$           if (bborbit) then
!!$              orbit(2) = orbit(2) + phix
!!$              orbit(4) = orbit(4) + phiy
!!$           elseif (bbd_flag .ne. 0)  then
!!$              bb_kick(1,bbd_cnt) = phix
!!$              bb_kick(2,bbd_cnt) = phiy
!!$           endif
!!$
!!$           !---- first-order effects.
!!$           phixx = (two / r2) * (- (xs * phix + ys * phiy)             &
!!$                + fk * (one - (sy / sx) * exk))
!!$           phixy = (two / r2) * (- (xs * phiy - ys * phix))
!!$           phiyy = (two / r2) * (+ (xs * phix + ys * phiy)             &
!!$                - fk * (one - (sx / sy) * exk))
!!$           re(2,1) = phixx
!!$           re(2,3) = phixy
!!$           re(4,1) = phixy
!!$           re(4,3) = phiyy
!!$
!!$           !---- second-order effects.
!!$           if (fsec) then
!!$              phixxx = (- phix - (xs * phixx + ys * phixy)              &
!!$                   + fk * xs * sy * exk / sx**3) / r2
!!$              phixxy = (- phiy - (xs * phixy - ys * phixx)) / r2
!!$              phixyy = (+ phix - (xs * phiyy - ys * phixy)) / r2
!!$              phiyyy = (+ phiy + (xs * phixy + ys * phiyy)              &
!!$                   - fk * ys * sx * exk / sy**3) / r2
!!$              te(2,1,1) = phixxx
!!$              te(2,1,3) = phixxy
!!$              te(2,3,1) = phixxy
!!$              te(4,1,1) = phixxy
!!$              te(2,3,3) = phixyy
!!$              te(4,1,3) = phixyy
!!$              te(4,3,1) = phixyy
!!$              te(4,3,3) = phiyyy
!!$           endif
!!$        endif
!!$
!!$        !---- no tracking desired.
!!$     else
!!$        re(2,1) = fk / (sx * (sx + sy))
!!$        re(4,3) = fk / (sy * (sx + sy))
!!$     endif
!!$     !---- centre option
!!$     if(centre_cptk.or.centre_bttk) then
!!$        call dcopy(orbit,orbit00,6)
!!$        call dcopy(re,re00,36)
!!$        call dcopy(te,te00,216)
!!$        if(centre_cptk) then
!!$           call dcopy(orbit,orbit0,6)
!!$           call twcptk(re,orbit0)
!!$        endif
!!$        if(centre_bttk) call twbttk(re,te)
!!$        call dcopy(orbit00,orbit,6)
!!$        call dcopy(re00,re,36)
!!$        call dcopy(te00,te,216)
!!$     endif
!!$  endif
!!$
!!$end SUBROUTINE tmbb_old
SUBROUTINE ccperrf(xx, yy, wx, wy)

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
  integer n,nc,nu
  double precision xx,yy,wx,wy,x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux,&
       rx(33),ry(33),cc,zero,one,two,half,xlim,ylim,fac1,fac2,fac3
  parameter(cc=1.12837916709551d0,zero=0d0,one=1d0,two=2d0,         &
       half=5d-1,xlim=5.33d0,ylim=4.29d0,fac1=3.2d0,fac2=23d0,fac3=21d0)

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

  !     if(y .eq. zero) wx = exp(-x**2)
  if(yy .lt. zero) then
     wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx
     wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy
     if(xx .gt. zero) wy = -wy
  else
     if(xx .lt. zero) wy = -wy
  endif

end SUBROUTINE ccperrf
double precision function gauss_erf(x)

  implicit none

  !---- returns the value of the Gauss error integral:
  !     1/sqrt(2*pi) Int[-inf, x] exp(-x**2/2) dx
  double precision x,xx,re,im,zero,one,two,half
  parameter(zero=0d0,one=1d0,two=2d0,half=5d-1)

  xx = x / sqrt(two)
  call ccperrf(zero,xx,re,im)
  gauss_erf = one - half * exp(-xx**2) * re
end function gauss_erf
SUBROUTINE twwmap(pos, orbit)

  use twissotmfi
  use twissafi
  implicit none

  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Save concatenated sectormap (kick, rmatrix, tmatrix)             *
  !     Input:                                                           *
  !     pos   (double)  position                                         *
  !     orbit (double)  current orbit                                    *
  !     Further input in twissotm.fi                                     *
  !----------------------------------------------------------------------*
  integer i, k, l
  double precision sum1, sum2, pos, orbit(6)
  double precision ek(6)
  double precision zero,two
  parameter(zero=0d0,two=2d0)

  !---- Track ORBIT0 using zero kick.
  do i = 1, 6
     sum2 = orbit(i)
     do k = 1, 6
        sum1 = zero
        do l = 1, 6
           sum1 = sum1 + stmat(i,k,l) * sorb(l)
        enddo
        sum2 = sum2 - (srmat(i,k) - sum1) * sorb(k)
        !     srmat(i,k) = srmat(i,k) - two * sum1
     enddo
     ek(i) = sum2
  enddo
  call dcopy(orbit, sorb, 6)

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
  call m66one(srmat)
  call dzero(stmat, 216)

end SUBROUTINE twwmap
SUBROUTINE tmdpdg(ftrk,orbit,fmap,ek,re,te)

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
  logical ftrk,fmap
  double precision fint,tilt,e1,h,hgap,corr,ek0(6),rw(6,6),tw(6,6,6),orbit(6),  &
       ek(6),re(6,6),te(6,6,6),bvk,node_value,zero,one
  parameter(zero=0d0,one=1d0)

  !-- LD: modified to include corrections from Frank.
  call dzero(ek0,6)
  call m66one(rw)
  call dzero(tw,216)

  e1 = node_value('e1 ')
  h = node_value('h ')
  bvk = node_value('other_bv ')
  h = bvk * h
  hgap = node_value('hgap ')
  fint = node_value('fint ')
  tilt = node_value('tilt ')
  corr = (h + h) * hgap * fint
  fmap = h .ne. zero .and. ( e1 .ne. zero .or. corr .ne. zero )
  if (.not. fmap) return
  !     print *,"------------------------------------------ "
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
!  if (ftrk) then
!     orbit(2) = orbit(2) + rw(2,1) * orbit(1)
!     orbit(4) = orbit(4) + rw(4,3) * orbit(3)
!  endif
  return
end SUBROUTINE tmdpdg
SUBROUTINE tmsol_th(ftrk,orbit,fmap,ek,re,te)

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
  logical ftrk,fmap
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),deltap,bvk,     &
       get_value,node_value,sk,skl,sks,sksl,cosTh,sinTh,Q0,Q,zero,one,two
  parameter(zero=0d0,one=1d0,two=2d0)

  fmap=.true.

  !
  !---- Initialize.
  deltap  = get_value('probe ','deltap ')
  !
  !---- Get solenoid parameters
  sksl    = node_value('ksi ')
  sks     = node_value('ks ')

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

end SUBROUTINE tmsol_th

SUBROUTINE tmnll(fsec,ftrk,orbit,fmap,ek,re,te)

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
  logical fsec,ftrk,fmap
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),node_value,get_variable,&
       pi,zero,one,two,knll,cnll,k1, &
       dd, u, v, dUu, dUv, dux, duy, dvx, dvy, x, y
  parameter(zero=0d0,one=1d0,two=2d0)

  !---- Matrix
  fmap = .true.
  cnll=node_value('cnll ')
  knll=node_value('knll ')/cnll
  k1=knll/cnll
  re(1,1)=one
  re(2,1)=-two*k1
  re(2,2)=one
  re(3,3)=one
  re(4,3)=two*k1
  re(4,4)=one
  re(5,5)=one
  re(6,6)=one

  !---- Second Order Terms
!  if (fsec) then
!  endif

  !---- Track orbit.
  if (ftrk) then
    pi=get_variable('pi ')
    x = orbit(1)/cnll
    y = orbit(3)/cnll
    u=0.5*sqrt((x-1)**2+y**2)+0.5*sqrt((x+1)**2+y**2)
    v=0.5*sqrt((x+1)**2+y**2)-0.5*sqrt((x-1)**2+y**2)
    if (u.eq.1d0) then 
       dd=0 
    else 
       dd=u**2*log(u+sqrt(u*u-1))/sqrt(u**2-1)
    endif
    dUu=(u+log(u+sqrt(u*u-1))*sqrt(u**2-1)+dd)/(u**2-v**2) &
     -2*u*(u*log(u+sqrt(u*u-1))*sqrt(u**2-1) &
     +v*(acos(v)-0.5*pi)*sqrt(1-v**2)) /(u**2-v**2)**2
    dUv=2*v*(u*log(u+sqrt(u*u-1))*sqrt(u**2-1) &
         +v*(acos(v)-0.5*pi)*sqrt(1-v**2)) /(u**2-v**2)**2 &
         -(v-(acos(v)-0.5*pi)*sqrt(1-v**2)+v**2*(acos(v)-0.5*pi)/sqrt(1-v**2))&
         /(u**2-v**2)
    dux=0.5*(x-1)/sqrt((x-1)**2+y**2) +0.5*(x+1)/sqrt((x+1)**2+y**2)
    duy=0.5*y/sqrt((x-1)**2+y**2) +0.5*y/sqrt((x+1)**2+y**2)
    dvx=0.5*(x+1)/sqrt((x+1)**2+y**2) -0.5*(x-1)/sqrt((x-1)**2+y**2)
    dvy=0.5*y/sqrt((x+1)**2+y**2) -0.5*y/sqrt((x-1)**2+y**2)

    orbit(2)=orbit(2)+knll*(dUu*dux+dUv*dvx)
    orbit(4)=orbit(4)+knll*(dUu*duy+dUv*dvy)
  endif

end SUBROUTINE tmnll

SUBROUTINE tmrfmult(fsec,ftrk,orbit,fmap,ek,re,te)

  use twtrrfi
  use twisslfi
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
  logical fsec,ftrk,fmap,dorad
  integer nord,iord,j,nn,ns,node_fd_errors
  double precision orbit(6),ek(6),re(6,6),te(6,6,6),elrad,beta,deltap,  &
       field(2,0:maxmul),normal(0:maxmul),orbit0(6),                    &
       skew(0:maxmul),arad,gammas,rfac,bvk,                             &
       zero,one,two,three, tilt, angle, cangle, sangle, dtmp
  double precision orbit00(6),ek00(6),re00(6,6),te00(6,6,6)
  double precision f_errors(0:maxferr)
  !double precisio vals(2,0:maxmul)
  double precision, external :: node_value,get_value,get_variable
  integer n_ferr
  
  !--- AL: RF-multipole
  integer dummyi
  double precision pc, krf, vrf
  double precision pi, twopi, clight, ten3m
  double precision x, y, z, dpx, dpy, dpt
  double precision freq, volt, lag, harmon
  double precision field_cos(2,0:maxmul)
  double precision field_sin(2,0:maxmul)
  double precision pnl(0:maxmul), psl(0:maxmul)
  complex*16 ii, Cm2, Sm2, Cm1, Sm1, Cp0, Sp0, Cp1, Sp1
    
  parameter ( zero=0d0, one=1d0, two=2d0, three=3d0, ten3m=1d-3)
  parameter ( ii=(0d0,1d0) )

  !---- Zero the arrays
  call dzero(normal,maxmul+1)
  call dzero(skew,maxmul+1)
  call dzero(pnl,maxmul+1)
  call dzero(psl,maxmul+1)
  call dzero(f_errors,maxferr+1)
  call dzero(field,2*(maxmul+1))
  
  !---- Read-in the parameters
  clight=get_variable('clight ')
  pi=get_variable('pi ')
  twopi=pi*two

  freq = node_value('freq ');
  volt = node_value('volt ');
  lag = node_value('lag ');
  harmon = node_value('harmon ');
  bvk = node_value('other_bv ')
  elrad = node_value('lrad ')
  deltap = get_value('probe ', 'deltap ')
  dorad = get_value('probe ','radiate ') .ne. zero
  arad = get_value('probe ','arad ')
  gammas = get_value('probe ','gamma ')
  beta = get_value('probe ','beta ')
  pc = get_value('probe ','pc ')
  tilt = node_value('tilt ')
  n_ferr = node_fd_errors(f_errors);
  call get_node_vector('knl ', nn, normal)
  call get_node_vector('ksl ', ns, skew)
  call get_node_vector('pnl ', dummyi, pnl); ! NOTE !!!!! THIS DOES NOT MAKE USE OF NODE->PNL and NODE->PSL
  call get_node_vector('psl ', dummyi, psl);
 
  rfac = zero
  fmap = .true.
  
  !---- Set-up some parameters
  !-- LD: 20.6.2014 (bvk=-1: not -V -> V but lag -> pi-lag)
  !-- AL: 30.6.2014 (bvk=-1: z -> -z)
  if (bvk .eq. -one) then
    orbit(5) = -orbit(5);
  endif

  krf = 2*pi*freq*1d6/clight;
  vrf = volt*ten3m/(pc*(one+deltap));
  
  if (n_ferr.gt.0) then
     call dcopy(f_errors,field,n_ferr)
  endif
  nord = max(nn, ns, n_ferr/2-1);
  
  !---- Particle's coordinates
  if (ftrk) then
    x = orbit(1);
    y = orbit(3);
    z = orbit(5);
  else
    x = zero;
    y = zero;
    z = zero;
  endif
  
  !---- Vector with strengths + field errors
  do iord = 0, nord;
     field_cos(1,iord) = bvk * (normal(iord) * cos(pnl(iord) * twopi - krf * z) + field(1,iord)) / (one + deltap);
     field_sin(1,iord) = bvk * (normal(iord) * sin(pnl(iord) * twopi - krf * z))                 / (one + deltap);
     field_cos(2,iord) = bvk * (skew(iord)   * cos(psl(iord) * twopi - krf * z) + field(2,iord)) / (one + deltap);
     field_sin(2,iord) = bvk * (skew(iord)   * sin(psl(iord) * twopi - krf * z))                 / (one + deltap);
     if (tilt.ne.zero)  then
        angle = (iord+1) * (-tilt);
        cangle = cos(angle);
        sangle = sin(angle);
        dtmp              = field_cos(1,iord) * cangle - field_cos(2,iord) * sangle;
        field_cos(2,iord) = field_cos(1,iord) * sangle + field_cos(2,iord) * cangle;
        field_cos(1,iord) = dtmp;
        dtmp              = field_sin(1,iord) * cangle - field_sin(2,iord) * sangle;
        field_sin(2,iord) = field_sin(1,iord) * sangle + field_sin(2,iord) * cangle;
        field_sin(1,iord) = dtmp;
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
  do iord = nord, 0, -1
    if (iord.ge.2) then
      Cm2 = Cm2 * (x+ii*y) / (iord-1) + field_cos(1,iord)+ii*field_cos(2,iord);
      Sm2 = Sm2 * (x+ii*y) / (iord-1) + field_sin(1,iord)+ii*field_sin(2,iord);
    endif
    if (iord.ge.1) then
      Cm1 = Cm1 * (x+ii*y) / (iord)   + field_cos(1,iord)+ii*field_cos(2,iord);
      Sm1 = Sm1 * (x+ii*y) / (iord)   + field_sin(1,iord)+ii*field_sin(2,iord);
    endif
    Cp0 = Cp0 * (x+ii*y) / (iord+1)   + field_cos(1,iord)+ii*field_cos(2,iord);
    Sp0 = Sp0 * (x+ii*y) / (iord+1)   + field_sin(1,iord)+ii*field_sin(2,iord);
    Cp1 = Cp1 * (x+ii*y) / (iord+2)   + field_cos(1,iord)+ii*field_cos(2,iord);
    Sp1 = Sp1 * (x+ii*y) / (iord+2)   + field_sin(1,iord)+ii*field_sin(2,iord);
  enddo
  Sp1 = Sp1 * (x+ii*y);
  Cp1 = Cp1 * (x+ii*y);
  
  !---- Track orbit.
  if (ftrk) then
     !---- The kick
     dpx = -REAL(Cp0);
     dpy = AIMAG(Cp0);
     dpt =  vrf * sin(lag * twopi - krf * z) - krf * REAL(Sp1);
     
     !---- Radiation effects at entrance.
     if (dorad  .and.  elrad .ne. zero) then
        rfac = arad * gammas**3 * (dpx**2+dpy**2) / (three*elrad)
        orbit(2) = orbit(2) - rfac * (one + orbit(6)) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + orbit(6)) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + orbit(6)) ** 2
     endif
     
     !---- Apply the kick
     orbit(2) = orbit(2) + dpx
     orbit(4) = orbit(4) + dpy
     orbit(6) = orbit(6) + dpt
  
     !---- Radiation effects at exit.
     if (dorad  .and.  elrad .ne. zero) then
        orbit(2) = orbit(2) - rfac * (one + orbit(6)) * orbit(2)
        orbit(4) = orbit(4) - rfac * (one + orbit(6)) * orbit(4)
        orbit(6) = orbit(6) - rfac * (one + orbit(6)) ** 2
     endif
  endif
  
  !---- Element Kick
  ek(2) = -REAL(Cp0);
  ek(4) = AIMAG(Cp0);
  ek(6) =  vrf * sin(lag * twopi - krf * z) - krf * REAL(Sp1);

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

  !-- AL: 30.6.2014 (bvk=-1: z -> -z)
  if (bvk .eq. -one) then
    orbit(5) = -orbit(5);
  endif
  
  !---- centre option
  if(centre_cptk.or.centre_bttk) then
     call dcopy(orbit,orbit00,6)
     call dcopy(ek,ek00,6)
     call dcopy(re,re00,36)
     call dcopy(te,te00,216)
     if(centre_cptk) then
        call dcopy(orbit,orbit0,6)
        call twcptk(re,orbit0)
     endif
     if(centre_bttk) call twbttk(re,te)
     call dcopy(orbit00,orbit,6)
     call dcopy(ek00,ek,6)
     call dcopy(re00,re,36)
     call dcopy(te00,te,216)
  endif

end SUBROUTINE tmrfmult

SUBROUTINE calcsyncint(rhoinv,blen,k1,e1,e2,betxi,alfxi,dxi,dpxi,I)
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

  double precision :: rhoinv, blen, k1, e1, e2, betxi, alfxi, dxi, dpxi
  double precision :: I(5)

  ! local variables
  double precision :: dx2, gamx, dispaverage, curlyhaverage
  double precision :: betx, alfx, dx, dpx

  complex*16 :: k2, k, kl

  integer :: get_option
 
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

  I(1) = dispaverage * rhoinv * blen
  I(2) = rhoinv*rhoinv * blen
  I(3) = abs(rhoinv)**3 * blen
  I(4) = dispaverage*rhoinv*(rhoinv**2 + 2*k1) * blen & 
           - rhoinv*rhoinv*(dx*tan(e1) + dx2*tan(e2))
  I(5) = curlyhaverage * abs(rhoinv)**3 * blen
 
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

  return
END SUBROUTINE calcsyncint
