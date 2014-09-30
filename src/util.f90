module Inf_NaN_Detection

  !!     Inf_NaN_Detection module
  !!     Copyright(c) 2003, Lahey Computer Systems, Inc.
  !!     Copies of this source code, or standalone compiled files
  !!     derived from this source may not be sold without permission
  !!     from Lahey Computers Systems. All or part of this module may be
  !!     freely incorporated into executable programs which are offered
  !!     for sale. Otherwise, distribution of all or part of this file is
  !!     permitted, provided this copyright notice and header are included.

  !!     This module exposes four elemental functions:
  !!
  !!     isnan(x)    - test for a "not a number" value
  !!
  !!     isinf(x)    - test for either a positive or negative "infinite" value
  !!
  !!     isposinf(x) - test for a positive "infinite" value
  !!
  !!     isneginf(x) - test for a negative "infinite" value
  !!
  !!     Each function accepts a single or double precision real argument, and
  !!     returns a true or false value to indicate the presence of the value
  !!     being tested for. If the argument is array valued, the function returns
  !!     a conformable logical array, suitable for use with the ANY function, or
  !!     as a logical mask.
  !!
  !!     Each function operates by transferring the bit pattern from a real
  !!     variable to an integer container. Unless testing for + or - infinity,
  !!     the sign bit is cleared to zero. The value is exclusive ORed with
  !!     the value being tested for. The integer result of the IEOR function is
  !!     converted to a logical result by comparing it to zero.
  !!
  implicit none

  private


  public :: isnan, isinf, isposinf, isneginf, sp, dp


  ! Order set-up
  integer,parameter::sp=kind(1.e0)
  integer,parameter::dp=selected_real_kind(2*precision(1.e0))

  ! Kind numbers for single and double precision integer containers
  integer, parameter :: Single = selected_int_kind(precision(1.e0))
  integer, parameter :: Double = selected_int_kind(2*precision(1.e0))

  !VK20070611: The below lines are not accepted by NAG-compiler with <Makefile_nag>
  ! Single precision IEEE values
  integer(Single), parameter :: sNaN    = Z"7FC00000"
  integer(Single), parameter :: sPosInf = Z"7F800000"
  integer(Single), parameter :: sNegInf = Z"FF800000"

  ! Double precision IEEE values
  integer(Double), parameter :: dNaN    = Z"7FF8000000000000"
  integer(Double), parameter :: dPosInf = Z"7FF0000000000000"
  integer(Double), parameter :: dNegInf = Z"FFF0000000000000"

  ! Locatation of single and double precision sign bit (Intel)
  ! Subtract one because bit numbering starts at zero
  integer, parameter :: SPSB = bit_size(sNaN) - 1
  integer, parameter :: DPSB = bit_size(dNaN) - 1

  interface isnan
     module procedure sisnan
     module procedure disnan
  end interface isnan

  interface isinf
     module procedure sisinf
     module procedure disinf
  end interface isinf

  interface isposinf
     module procedure sisposinf
     module procedure disposinf
  end interface isposinf

  interface isneginf
     module procedure sisneginf
     module procedure disneginf
  end interface isneginf

contains

  ! Single precision test for NaN
  elemental function sisnan(x) result(res)
    implicit none
    real(sp), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sNan),SPSB), sNaN) == 0
  end function sisnan

  ! Double precision test for NaN
!DEC$ ATTRIBUTES FORCEINLINE :: disnan 
  elemental function disnan(d) result(res)
    implicit none
    real(dp), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dNaN),DPSB), dNaN) == 0
  end function disnan

  ! Single precision test for Inf
  elemental function sisinf(x) result(res)
    implicit none
    real(sp), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sPosInf),SPSB), sPosInf) == 0
  end function sisinf

  ! Double precision test for Inf
  elemental function disinf(d) result(res)
    implicit none
    real(dp), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dPosInf),DPSB), dPosInf) == 0
  end function disinf

  ! Single precision test for +Inf
  elemental function sisposinf(x) result(res)
    implicit none
    real(sp), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sPosInf), sPosInf) == 0
  end function sisposinf

  ! Double precision test for +Inf
  elemental function disposinf(d) result(res)
    implicit none
    real(dp), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dPosInf), dPosInf) == 0
  end function disposinf

  ! Single precision test for -Inf
  elemental function sisneginf(x) result(res)
    implicit none
    real(sp), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sNegInf), sNegInf) == 0
  end function sisneginf

  ! Double precision test for -Inf
  elemental function disneginf(d) result(res)
    implicit none
    real(dp), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dNegInf), dNegInf) == 0
  end function disneginf

end module Inf_NaN_Detection

module bbfi
  implicit none
  public
  integer bbd_max
  parameter(bbd_max=100000)
  integer :: bbd_loc(bbd_max)=0,bbd_cnt=0,bbd_flag=0,bbd_pos=0
  double precision :: bb_kick(2,bbd_max)=0.d0
end module bbfi
module deltrafi
  implicit none
  public
  logical :: dorad=.false.,dodamp=.false.,dorand=.false.,fastune=.false.
  double precision :: deltax=0.d0
end module deltrafi
module dyntabfi
  implicit none
  public
  double precision :: dynapfrac=0.d0,dktrturns=0.d0,xend=0.d0,pxend=0.d0,&
       yend=0.d0,pyend=0.d0,tend=0.d0,ptend=0.d0,smear=0.d0,yapunov=0.d0
end module dyntabfi
module wmaxmin0fi
  implicit none
  public
  double precision :: wxmax=0.d0,wxmin=0.d0,wymax=0.d0,wymin=0.d0,&
       wxymax=0.d0,wxymin=0.d0
end module wmaxmin0fi
module tunesfi
  implicit none
  public
  double precision :: x0=0.d0,y0=0.d0,tunx=0.d0,tuny=0.d0,dtune=0.d0
end module tunesfi
module twiss0fi
  implicit none
  public
  integer align_max,fundim
  parameter(align_max=14,fundim = 74)
end module twiss0fi
module twissafi
  implicit none
  public
  character(48) :: table_name=' ',sectorTableName=' '
  logical :: match_is_on=.false.
end module twissafi
module twisslfi
  implicit none
  public
  logical ::  centre=.false.,centre_cptk=.false.,centre_bttk=.false.,first,&
       rmatrix=.false.,sectormap=.false.,ripken=.false.
end module twisslfi
module twisscfi
  use twiss0fi
  implicit none
  public
  double precision :: opt_fun0(fundim)=0.d0,opt_fun(fundim)=0.d0,disp(6)=0.d0,&
       ddisp(6)=0.d0,rmat(2,2)=0.d0,betx=0.d0,alfx=0.d0,amux=0.d0,bety=0.d0,&
       alfy=0.d0,amuy=0.d0,bxmax=0.d0,dxmax=0.d0,bymax=0.d0,dymax=0.d0,&
       xcomax=0.d0,ycomax=0.d0,sigxco=0.d0,sigyco=0.d0,sigdx=0.d0,sigdy=0.d0,&
       wgt=0.d0,cosmux=0.d0,cosmuy=0.d0,wx=0.d0,phix=0.d0,dmux=0.d0,wy=0.d0,&
       phiy=0.d0,dmuy=0.d0,synch_1=0.d0,synch_2=0.d0,synch_3=0.d0,synch_4=0.d0,&
       synch_5=0.d0,suml=0.d0,circ=0.d0,eta=0.d0,alfa=0.d0,gamtr=0.d0,qx=0.d0,&
       qy=0.d0,sinmux=0.d0,sinmuy=0.d0,xix=0.d0,xiy=0.d0,currpos=0.d0
end module twisscfi
module twissotmfi
  implicit none
  public
  double precision :: rotm(6,6)=0.d0,rw(6,6)=0.d0,skick(6)=0.d0,sorb(6)=0.d0,&
       srmat(6,6)=0.d0,stmat(6,6,6)=0.d0
end module twissotmfi
module max_iterate
  implicit none
  public
  integer MAXITER
  parameter(MAXITER=150)
end module max_iterate
module twiss_elpfi
  implicit none
  public
  !---fixed positions for element parameters-----------------------------*
  double precision :: g_elpar(50)=0.d0
  !-general--------------------------------------------------------------*
  integer    g_el, g_kmax, g_kmin, g_calib, g_polarity
  parameter (g_el = 2, g_kmax = 3, g_kmin = 4, g_calib = 5, g_polarity = 6)
  !-bend-----------------------------------------------------------------*
  integer    b_angle , b_tilt , b_k0 , b_k0s ,                      &
       b_k1 , b_k1s , b_e1 , b_e2 , b_k2 ,                    &
       b_k2s , b_h1 , b_h2 , b_hgap ,                         &
       b_fint , b_fintx , b_k3 , b_k3s
  parameter (b_angle = 7, b_tilt = 8, b_k0 = 9, b_k0s = 10,         &
       b_k1 = 11, b_k1s = 12, b_e1 = 13 , b_e2 = 14, b_k2 =15,&
       b_k2s = 16, b_h1 = 17, b_h2 = 18, b_hgap = 19,         &
       b_fint = 20, b_fintx = 21, b_k3 = 22, b_k3s = 23)
  !-quad-----------------------------------------------------------------*
  integer    q_tilt, q_k1 , q_k1s
  parameter (q_tilt = 7, q_k1 = 8, q_k1s = 9)
  !-sext-----------------------------------------------------------------*
  integer    s_tilt, s_k2 , s_k2s
  parameter (s_tilt = 7, s_k2 = 8, s_k2s = 9)
  !-oct------------------------------------------------------------------*
  integer    o_tilt, o_k3 , o_k3s
  parameter (o_tilt = 7, o_k3 = 8, o_k3s = 9)
  !-mult-----------------------------------------------------------------*
  integer    m_tilt, m_lrad
  parameter (m_tilt = 7, m_lrad = 8)
  !-sol------------------------------------------------------------------*
  integer    so_lrad, so_ks, so_ksi
  parameter (so_lrad = 7, so_ks = 8, so_ksi = 9)
  !-rfc------------------------------------------------------------------*
  integer    r_volt, r_lag , r_freq
  parameter (r_volt = 7, r_lag = 8, r_freq = 9)
  !-elsep----------------------------------------------------------------*
  integer    e_tilt, e_ex , e_ey
  parameter (e_tilt = 7, e_ex = 8, e_ey = 9)
  !-hkick----------------------------------------------------------------*
  integer    h_tilt, h_lrad , h_kick, h_hkick, h_chkick
  parameter (h_tilt = 7, h_lrad = 8, h_kick = 9, h_hkick = 10,       &
       h_chkick = 11)
  !-vkick----------------------------------------------------------------*
  integer    v_tilt, v_lrad , v_kick, v_vkick, v_cvkick
  parameter (v_tilt = 7, v_lrad = 8, v_kick = 9, v_vkick = 10,       &
       v_cvkick = 11)
  !-kick-----------------------------------------------------------------*
  integer    k_tilt, k_lrad , k_hkick, k_vkick, k_chkick, k_cvkick
  parameter (k_tilt = 7, k_lrad = 8, k_hkick = 9, k_vkick = 10,      &
       k_chkick = 11, k_cvkick = 12)
end module twiss_elpfi
module emitfi
  implicit none
  public
  double precision :: qx=0.d0,qy=0.d0,qs=0.d0,cg=0.d0,sum(3)=0.d0,sumu0=0.d0
  save qx, qy, qs, cg,sum,sumu0
end module emitfi
module twtrrfi
  implicit none
  public
  !---- maxmul is the maximum multipole order both in twiss and trrun
  integer maxmul,maxferr,maxnaper
  parameter(maxmul=20,maxferr=50,maxnaper=100)
end module twtrrfi
module ibsdbfi
  implicit none
  public
  integer :: bunch=0
  double precision :: circ=0.d0,clight=0.d0,arad=0.d0,freq0=0.d0,alpha=0.d0,&
       amass=0.d0,charge=0.d0,en0=0.d0,gammas=0.d0,gamma=0.d0,ex=0.d0,ey=0.d0,&
       et=0.d0,sigt=0.d0,sige=0.d0,betas=0.d0,beta=0.d0,parnum=0.d0,&
       currnt=0.d0,sigx=0.d0,sigy=0.d0,alfa=0.d0
end module ibsdbfi
module matchfi
  implicit none
  public
  integer :: icovar=0,ilevel=0
  double precision :: edm=0.d0,fmin=0.d0
end module matchfi
module name_lenfi
  implicit none
  public
  integer name_len
  parameter(name_len=48)
end module name_lenfi
module physconsfi
  implicit none
  public
  double precision :: amu0=0.d0,elamda=0.d0,emass=0.d0,eps0=0.d0,erad=0.d0,&
       hbar=0.d0,plamda=0.d0,pmass=0.d0,prad=0.d0,qelect=0.d0,mumass=0.d0
end module physconsfi
module touschekfi
  implicit none
  public
  integer :: bunch=0
  double precision :: circ=0.d0,clight=0.d0,arad=0.d0,freq0=0.d0,amass=0.d0,&
       charge=0.d0,en0=0.d0,gammas=0.d0,gamma=0.d0,ex=0.d0,ey=0.d0,et=0.d0,&
       sigt=0.d0,sige=0.d0,betas=0.d0,beta=0.d0,parnum=0.d0,currnt=0.d0,&
       alfa=0.d0,um1=0.d0,deltap=0.d0,fb1=0.d0,fb2=0.d0
end module touschekfi
module trackfi
  implicit none
  public
  double precision :: arad=0.d0,betas=0.d0,beti=0.d0,gammas=0.d0,dtbyds=0.d0,&
       deltas=0.d0,bet0=0.d0,bet0i=0.d0,t_max=1.d20,pt_max=1.d20
  logical :: dodamp=.false.,dorad=.false.,dorand=.false.,fsecarb=.false.
  save arad,betas,beti,gammas,dtbyds,bet0,bet0i
end module trackfi
module time_varfi
  use twtrrfi
  use name_lenfi
  implicit none
  public
  logical time_var_m,time_var_p,time_var_c
  integer n_time_var
  parameter (n_time_var = 10000)
  integer time_var_m_cnt,time_var_p_cnt,time_var_c_cnt,                  &
       time_var_m_lnt,time_var_p_lnt,time_var_c_lnt,                     &
       trrun_nt
  double precision myfield(n_time_var,2,0:maxmul),                       &
       phase_tromb(n_time_var,36),cav_volt(n_time_var),                  &
       time_var_m_ind(n_time_var),time_var_p_ind(n_time_var),            &
       time_var_c_ind(n_time_var),                                       &
       time_var_m_nt(n_time_var),time_var_p_nt(n_time_var),              &
       time_var_c_nt(n_time_var)
  character*(name_len) time_var_m_ch(n_time_var),                        &
       time_var_p_ch(n_time_var),time_var_c_ch(n_time_var)
  save time_var_m_cnt,time_var_p_cnt,time_var_c_cnt,                     &
       time_var_m_lnt,time_var_p_lnt,time_var_c_lnt,trrun_nt,            &
       myfield,phase_tromb,cav_volt,time_var_m_ind,                      &
       time_var_p_ind,time_var_c_ind,time_var_m_nt,time_var_p_nt,        &
       time_var_c_nt,time_var_m_ch,time_var_p_ch,time_var_c_ch
end module time_varfi
module spch_bbfi
  use name_lenfi
  use bbfi
  implicit none
  public
  logical :: lost_in_turn = .false., is_lost = .false.
  integer i_turn, N_macro_surv, N_for_I, N_macro_max, N_spch, i_spch 
  parameter(N_macro_max=16000)
  double precision Ex_rms, Ey_rms, sigma_p, sigma_z
  double precision Ix_array(N_macro_max), Iy_array(N_macro_max),    &
       dpi_array(N_macro_max),                          &
       z_part_array(N_macro_max)
  double precision alpha, I_div_E_sum_max
!  parameter(alpha=0.0, I_div_E_sum_max=7.0)
  double precision betx_bb(bbd_max), bety_bb(bbd_max),              &
       alfx_bb(bbd_max), alfy_bb(bbd_max),              &
       gamx_bb(bbd_max), gamy_bb(bbd_max),              &
       dx_bb(bbd_max),   dy_bb(bbd_max)
  double precision 	rat_bb_n_ions
  double precision   sigma_t, mean_t  ! calculate and transfer to BB
  character*(name_len) spch_bb_name(bbd_max)
  save i_turn,N_macro_surv,N_for_I,N_spch,i_spch,                                &
        Ex_rms,Ey_rms,sigma_p,sigma_z,                                           &
        Ix_array,Iy_array,dpi_array, z_part_array,                               &
        betx_bb,bety_bb,alfx_bb,alfy_bb,gamx_bb,gamy_bb,dx_bb,dy_bb,             &
        rat_bb_n_ions,sigma_t, mean_t,spch_bb_name
  data rat_bb_n_ions / 1d0 /
end module spch_bbfi
module plotfi
  implicit none
  public
  !--- m_adble is the number of different types of elements
  integer mtype, m_adble
  parameter (mtype = 50, m_adble = 20)

  !--- mcnam adjusted to NAME_L
  integer mcnam, maxpnt
  parameter (mcnam = 48, maxpnt = 500)

  !--- szcompar is the size of the arrays returned by the routine comm_para
  integer szcompar
  parameter (szcompar = 100)

  !--- szchara is the size of the character strings char_a & version in
  !--- the routine pesopt
  integer szchara
  parameter (szchara = 400)

  !--- character sizes:
  !   MLSIZE    label character height
  !   MTSIZE    text   - " -
  !   MASIZE    annotation - " -
  integer mlsize, mtsize, masize
  parameter  (mlsize = 13,mtsize = 13,masize = 20)

  !--- parameters used in the routine pegacn in file plot.F

  integer mposx, mposy, mpost
  parameter (mposx = 8, mposy = 3, mpost = mposx * mposy)

  !--- parameters used in the routine peschm in file plot.F

  integer mobj, msize
  parameter (mobj = 14, msize = 88)

  integer mtitl, mxlabl, mnvar, mxdep, mqadd, mntmax, mksmax, mplred, mplout

  parameter (mtitl  = 128, mxlabl = 160)
  parameter (mnvar = 74, mxdep = 2)
  parameter (mqadd = 100000)
  parameter (mntmax = 20, mksmax = 10)
  parameter (mplred = 46, mplout = 47)

  integer maxseql, mtwcol, mpparm, mxcurv, mopt, mfile, marg, maxarg, mxdp, mxplot

  parameter (maxseql = 20000, mtwcol = 46, mpparm = 10,             &
       mxcurv = 10, mopt = 60, mfile = 120, marg = 60, maxarg = 1000,    &
       mxdp = 25, mxplot = 100)

  integer mintpl
  parameter (mintpl = 18)

  !--- Definition of common / peaddi /
  !--- itbv set in routine pesopt, used in routines pemima, peplot
  !---      and pesopt.
  !--- nivvar set in routine pesopt, used in routines pefill, peintp,
  !---        pemima, peplot and pesopt.
  !--- nelmach set in routine pefill, used in routines pefill, peplot.
  !--- numax set in routine pemima, used in routines pemima, peplot.
  !--- interf set in routine pesopt, used in routines pecurv, pefill.
  !--- noline set in routine pesopt, used in routines pefill, pesopt.
  !--- nqval set in routine pefill and peintp, used in routines pefill,
  !---        peintp, pemima and peplot.
  !--- nvvar set in routine pesopt and pemima, used in routine pemima.
  !--- nrrang set in routine pefill, used in routine pesopt and pefill.
  !--- proc_flag set in routine pesopt, used in routine pefill, peintp
  !---           and pesopt.
  !--- ipparm set in routine peintp and pesopt, used in routine peplot
  !---           and pesopt.
  !--- naxref set in routine pemima and pesopt, used in routine pemima
  !---           and pesopt.
  !--- ieltyp set in routine pefill, used in routine psplot.

  integer :: itbv=0,nivvar=0,nelmach=0,numax=0,interf=0,noline=0, &
       nqval(mxcurv)=0,nvvar(4)=0,nrrang(2)=0,                         &
       proc_flag(2,mxcurv)=0,ipparm(mpparm,mxcurv)=0,                  &
       naxref(mxcurv)=0,ieltyp(maxseql)=0

  !--- Definition of common / peaddr /
  !--- qascl set in routine pesopt, used in routine peplot.
  !--- qlscl set in routine pesopt, used in routine peplot.
  !--- qsscl set in routine pesopt, used in routine peplot.
  !--- qtscl set in routine pesopt, used in routine peplot.
  !--- hrange set in routine pesopt, used in routines pefill, peplot.
  !--- vrange set in routine pesopt, used in routine peplot.
  !--- hmima set in routines pesopt and pemima,
  !---       used in routines peplot, pesopt and pemima.
  !--- vmima set in routines pesopt and pemima,
  !---       used in routines peplot and pemima.
  !--- qhval set in routines pefill and peintp,
  !---       used in routines peplot, pefill and pemima.
  !--- qvval set in routines pefill and peintp,
  !---       used in routines peplot, pefill and pemima.
  !--- estart set in routine pefill, and peintp,
  !---        used in routines peplot and pefill.
  !--- eend set in routine pefill, used in routine peplot.

  real :: qascl=0.,qlscl=0.,qsscl=0.,qtscl=0.,                    &
       hrange(2)=0.,vrange(2,4)=0.,hmima(2)=0.,vmima(2,4)=0.,          &
       qhval(maxseql,mxcurv)=0.,qvval(maxseql,mxcurv)=0.,              &
       estart(maxseql)=0.,eend(maxseql)=0.

  !--- Definition of common / peaddc /
  !--- horname set in routine pesopt,
  !---         used in routines pefill, peplot and pesopt.
  !--- tabname set in routine pesopt,
  !---         used in routines pefill, peintp, pelfill and pesopt.
  !--- toptitle set in routine pesopt, used in routine peplot.
  !--- plfnam set in routine plginit, used in routines plotit and plginit.
  !--- axlabel set in routine pemima, used in routine peplot.
  !--- sname set in routine pesopt,
  !---         used in routines pefill and pesopt.
  !--- slabl set in routine pesplit,
  !---         used in routines peplot, pemima and pesopt.

  character(mfile) :: plfnam=' '
  character(mcnam) :: horname=' ',tabname=' ',sname(mxcurv)=' ',slabl(mxcurv)=' '
  character(mxlabl) :: axlabel(4)=' '
  character(mtitl) :: toptitle=' '

  save itbv,nivvar,nelmach,numax,interf,noline,nqval,nvvar,nrrang,proc_flag,&
       ipparm,naxref,ieltyp,qascl,qlscl,qsscl,qtscl,hrange,vrange,hmima,&
       vmima,qhval,qvval,estart,eend,horname,tabname,toptitle,plfnam,axlabel,&
       sname,slabl
end module plotfi
module plot_bfi
  implicit none
  public
  !--- Definition of common / peotcl /
  !--- fpmach set in routines pesopt and pefill, used in routine peplot
  !--- ddp_flag set in routine pefill, used in routine peplot
  !--- ptc_flag set in routines pesopt, used in routine pefill

  logical :: fpmach=.false.,dpp_flag=.false.,ptc_flag=.false.
  save fpmach,dpp_flag,ptc_flag
end module plot_bfi
module plot_cfi
  implicit none
  public
  !--- Definition of common / e2save /
  !--- e2s initialised in routine pefill, used in routine peelma

  double precision :: e2s=0.d0
end module plot_cfi
module plot_mathfi
  implicit none
  public
  !--- Definitions of mathematical constants

  double precision pi, zero, eps, one, two, twopi, half
  parameter         (pi = 3.1415926535898d0)
  parameter         (zero = 0.d0, half = 0.5d0, eps = 1.d-5)
  parameter         (one = 1.d0, two = 2.d0, twopi = two * pi)
end module plot_mathfi
module resindexfi
  implicit none
  public
  integer mnres,mymorder
  parameter (mnres=1000,mymorder=20)
end module resindexfi
module gxx11_common
  implicit none
  public
  !
  integer madim1,madim2,maxset,mconid,merrun,metaun,miunit,mmetat,&
       normt,mounit,mpaxs,mpcurv,mtermt,mtick,mtmeta,mtterm,mxaxs,mxpix,&
       mxsize,myaxs,mypix,mysize,mnormt,mx11pr,mx11tf,mxxpix,mxypix,&
       mcolor,mpspro,meppro,mdict,mlpro,&
       mpsep,mepep,mhead,mline,msfact,mlbb1,mlbb2,mubb1,mubb2,mtfont,&
       mwid1,mwid2
  real toleps,versio
  parameter (mxaxs = 4, myaxs = 4, mpaxs = 23, mpcurv = 10,&
       maxset = 20, mtterm = 1, mmetat = 4,&
       mtermt = 101, mtmeta = 2, mconid = 7, mtick = 10, metaun = 11,&
       mxpix = 1000, mypix = 1000, mxsize = 27, mysize = 19,&
       madim1 = 500, toleps = 1.e-5,&
       merrun = 10, miunit = 5, mounit = 6, versio = 1.50,&
       mx11pr = 10, mx11tf = 10, mxxpix = 1200, mxypix = 1000,&
       mcolor = 6, mpspro = 8, meppro = 8, mdict = 24, mlpro = 68,&
       mpsep = 3, mepep = 2, mhead = 4, mline = 72, msfact = 4,&
       mlbb1 = 17, mlbb2 = 30, mubb1 = 573, mubb2 = 790, mtfont = 12,&
       mwid1 = mubb1 - mlbb1, mwid2 = mubb2 - mlbb2 )
  parameter (mnormt = 2, madim2 = 100)
  !
  integer :: &
       itermt=0, interm=0, inmeta=0, ierrun=0, imetun=0, inunit=0, iounit=0, ipage=0,&
       isfflg=0, isqflg=0, iwtflg=0, iclflg=0, inormt=0, ipseps=0, iepsop=0, itseop=0,&
       iepscf=0, imetps=0, ipctct=0, iczebr=0, idinit=0, ipstyp=0, iclear=0, istotx=0,&
       lmpict=0, ltermt=0, lnterm=0, lnmeta=0, lerrun=0, lmetun=0, lnunit=0, lounit=0,&
       lsfflg=0, lsqflg=0, lwtflg=0, lclflg=0, lnormt=0, lmetax=0, lmetay=0, lmetnm=0,&
       lerrnm=0, ldefnl=0, lerrop=0, lmetop=0, ltotin=0, lacttm=0, lpseps=0, lundef=0,&
       lttime=0, ldinit=0, ltseop=0,&
       ixapar(mpaxs,mxaxs)=0, iyapar(mpaxs,myaxs)=0, icvpar(mpcurv ,maxset)=0
  !
  integer :: &
       nxpix=0, nypix=0, lxpix=0, lypix=0, icucol=0, iorips=0,      &
       iutlps=0, ibbox(4)=0, ix11pr(mx11pr)=0, ix11tf(mx11tf)=0, ix11op(mx11tf)=0  !
  !
  real :: &
       fxpix=0., fypix=0., rx11pr(mx11pr)=0., rgbcol(3,mcolor)=0.
  !
  real :: &
       xmetaf=0., ymetaf=0., xsterm=0., ysterm=0., wfact=0., wttime=0., wxfact=0., wyfact=0.,&
       vpfacx=0., vpfacy=0.,&
       vptdef(4)=0., vploc(4)=0., actwnd(4)=0., rangex(2,mxaxs)=0., rangey(2 ,myaxs)=0.,&
       cvwnwd(4,maxset)=0., axwndx(2,maxset), axwndy(2,maxset)=0.
  !
  character :: &
       smetnm*256=" ", serrnm*256=" ", sxtext(mxaxs)*300=" ", sytext(myaxs)*300=" ",&
       sxform(mxaxs)*20=" ", syform(myaxs)*20=" ", splotc*(maxset)=" ", stortx * 20=" ",&
       sdefnl*1=" ",spsnam * 256=" ", colour(mcolor) * 16=" ", pshead(mhead) * 60=" "
  !
  real :: xp(madim2+1)=0.,xvp(madim2+1)=0,yp(madim2+1)=0.,yvp(madim2+1)=0
  !
  real :: p(madim1,2)=0.,s(madim1)=0.,yy1d(madim1,2)=0.,yy2d(madim1,2)=0.
end module gxx11_common
module gxx11_aux
  implicit none
  public
  !
  character(100) strloc
  integer, dimension(14) :: ivals=(/ 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 1, 0 /)
  !   ivals(1)       marker type
  !     (2)       fill area interior style
  !     (3)       horizontal text alignment
  !     (4)       vertical text alignment
  !     (5)       text font
  !     (6)       text precision
  !     (7)       marker colour index
  !     (8)       metafile status (0 closed, 1 open)
  !     (9)       text colour index
  !    (10)       free
  !    (11)       polyline colour index
  !    (12)       polyline style
  !    (13)       current normalisation transformation number
  !    (14)       last call type: 0 undef., 1 line, 2 text, 3 marker

  real, dimension(14) :: rvals=(/ 0., 1., 0.01, 0., 1., 0., 1., 0., 1., 0., 1., 1., 1., 1. /)
  !   rvals(1-2)  chup vector
  !     3         character height
  !     4-7       window
  !     8-11      viewport
  !     12        character expansion factor
  !     13        line width scale factor
  !     14        marker scale factor
  save    ivals, rvals
end module gxx11_aux
module fasterror
  implicit none
  logical :: fasterror_on = .false.
  integer idim,nx,ny,kstep
  double precision hrecip,wtimag,wtreal
  parameter ( nx = 490, ny = 470 )
  parameter ( idim = (nx+2)*(ny+2) )
  public
  common /wzcom1/ hrecip, kstep
  common /wzcom2/ wtreal(idim), wtimag(idim)
end module fasterror
subroutine fort_info(t1, t2)
  implicit none


  character(*) t1, t2
  integer get_option
  if (get_option('info ') .ne. 0 .and. get_option('warn ') .ne. 0)  &
       print '(a,1x,a,1x,a)', '++++++ info:', t1, t2
end subroutine fort_info
subroutine fort_warn(t1, t2)
  implicit none



  character(*) t1, t2
  integer get_option
  if (get_option('warn ') .ne. 0) then
     print '(a,1x,a,1x,a)', '++++++ warning:', t1, t2
     call augmentfwarn()
  endif
end subroutine fort_warn
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
  use twiss0fi
  implicit none

  double precision orbit0(6), rt(6,6), tt(6,6,6)
  double precision opt(fundim)
  integer error
  call m66one(rt)
  call dzero(opt,fundim)
  call tmclor(orbit0, .true., .true., opt, rt, tt, error)
end subroutine getclor
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

end subroutine m66add
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

end subroutine m66byv
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

end subroutine m66cpy
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

end subroutine m66div
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
  double precision b(6,6),c(6,6),source(6,6),target(6,6),one,two,twelve
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

end subroutine m66exp
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

end subroutine m66inv
subroutine m66symp(r,nrm)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Check if a 6 by 6 matrix R is symplectic.
  ! Input:
  !   r(6,6)    (double)  Matrix R to check
  ! Output:
  !   nrm       (double)  The column norm of R'*J*R-J
  !----------------------------------------------------------------------*
  double precision R(6,6),J(6,6),T(6,6),nrm,z,o,n
  parameter(z=0d0,o=1d0,n=-1d0)
  J = reshape((/ z, o, z, z, z, z, &
               & n, z, z, z, z, z, &
               & z, z, z, o, z, z, &
               & z, z, n, z, z, z, &
               & z, z, z, z, z, o, &
               & z, z, z, z, n, z /), shape(J))
  call m66trm(R,J,T)
  call m66mpy(T,R,T)
  call m66sub(T,J,T)
  call m66nrm(T,nrm)
end subroutine m66symp
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

end subroutine m66mak
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

end subroutine m66mpy
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

end subroutine m66mtr
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

end subroutine m66nrm
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

end subroutine m66one
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

end subroutine m66ref
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

end subroutine m66scl
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

end function m66sta
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

end subroutine m66sub
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

end subroutine m66trm
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

end subroutine m66tp
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

end subroutine m66zro
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

9999 end subroutine solver
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

end subroutine symsol
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
       one,two,four,big,eps
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
130  if (m .ne. l) then
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
170  p = eigen(l) + f
     do i = l, 2, -1
        if (p .ge. eigen(i-1)) go to 190
        eigen(i) = eigen(i-1)
     enddo
     i = 1
190  eigen(i) = p
  enddo
300 continue

end subroutine symeig
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

end subroutine dcopy
subroutine dzero(vector,n)
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Zero an array.                                                     *
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

end subroutine dzero
subroutine aawarn(rout,text)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Print warning message.                                             *
  ! Input:                                                               *
  !   ROUT      (char)    Calling routine name.                          *
  !   TEXT      (char)    Message.                                       *
  !----------------------------------------------------------------------*
  character(*) rout,text

  print *,  '++++++ warning: ',rout,text
  call augmentfwarn()

end subroutine aawarn
subroutine aafail(rout,text)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Print fatal error message.                                         *
  ! Input:                                                               *
  !   ROUT      (char)    Calling routine name.                          *
  !   TEXT      (char)    Message.                                       *
  !----------------------------------------------------------------------*
  character(*) rout,text
  print *,' '
  print *,  '+-+-+- fatal: ',rout,text
  print *,' '
  print *,' '
  stop  1

end subroutine aafail
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

end function proxim
character(48) function charconv(tint)
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
  character(m) letter
  data letter / &
       '                                !"#$%&''()*+,-./0123456789:;<=>?@&
       &ABCDEFGHIJKLMNOPQRSTUVWXYZ[ ]^_`abcdefghijklmnopqrstuvwxyz{|}~'/
  charconv = ' '
  n = tint(1)
  do i = 1, n
     j = tint(i+1)
     if (j .lt. m)  charconv(i:i) = letter(j:j)
  enddo
end function charconv
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
       d(6),dx,dy,pb,reval(6),s,tm(6,6),zero,one
  parameter(zero=0d0,one=1d0,ilo=1,ihi=4,mdim=6,nn=4)

  !---- Compute eigenvalues and vectors.
  call m66cpy(fm,tm)
  call m66one(am)
  call orthes(mdim,nn,ilo,ihi,tm,d)
  call ortran(mdim,nn,ilo,ihi,tm,d,am)
  call hqr2(mdim,nn,ilo,ihi,tm,reval,aival,am,info)
  if (info .ne. 0) then
     write (6, 910) ((fm(i,k), k = 1, 6), i = 1, 6)
910  format('Unable to find eigenvalues for matrix:'/(6f12.6))
     call aafail('LASEIG',' Unable to find eigenvalues for matrix')
     go to 999
  endif
  !---- Normalize the eigenvectors.
  do k = 1, 5, 2
     pb = zero
     do ipind = 2, 6, 2
        iqind = ipind - 1
        pb = pb + am(iqind,k) * am(ipind,k+1) - am(ipind,k) * am(iqind,k+1)
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
999 end subroutine laseig
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
       d(6),dt,dx,dy,pb,reval(6),s,tm(6,6),zero
  parameter(zero=0d0,ilo=1,ihi=6,mdim=6,nn=6)

  !---- Compute eigenvalues and eigenvectors.
  call m66cpy(fm,tm)
  call orthes(mdim,nn,ilo,ihi,tm,d)
  call ortran(mdim,nn,ilo,ihi,tm,d,am)
  call hqr2(mdim,nn,ilo,ihi,tm,reval,aival,am,info)
  if (info .ne. 0) then
     write (6, 910) ((fm(i,k), k = 1, 6), i = 1, 6)
910  format('Unable to find eigenvalues for matrix:'/(6f12.6))
     call aafail('LADEIG',' Unable to find eigenvalues for matrix')
     go to 9999
  endif
  !---- Normalize the eigenvectors.
  do k = 1, 5, 2
     pb = zero
     do i = 1, 5, 2
        pb = pb + am(i,k) * am(i+1,k+1) - am(i+1,k) * am(i,k+1)
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
9999 end subroutine ladeig
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
end subroutine orthes
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
end subroutine ortran
subroutine hqr2(ndim,n,ilow,iupp,h,wr,wi,vecs,ierr)
  use max_iterate
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
       tempr,vecs(ndim,n),vi,vr,w,wi(n),wr(n),x,y,z,epsmch,zero,one,two, &
       triqua,fac1
  parameter(epsmch=1d-16,zero=0d0,one=1d0,two=2d0,triqua=.75d0,fac1=.4375d0)

  !Initialize
  z=zero
  s=zero
  p=zero
  q=zero
  r=zero

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
        if (abs(h(l,l-1)) .le. epsmch * (abs(h(l-1,l-1)) + abs(h(l,l)))) go to 100
     enddo
     l = ilow
100  continue
     x = h(ien,ien)
     if (l .eq. ien) go to 270
     y = h(na,na)
     w = h(ien,na) * h(na,ien)
     if (l .eq. na) go to 280
     if (its .eq. MAXITER) then
        write(6,*) "Maximum Iteration exceeded in HQR2, increase MAXITER: ",MAXITER
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
             * (abs(h(m-1,m-1)) + abs(z) + abs(h(m+1,m+1)))) go to 150
     enddo
150  continue
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
270  h(ien,ien) = x + t
     wr(ien) = h(ien,ien)
     wi(ien) = zero
     ien = na
     go to 60
     !==== Two roots (real pair or complex conjugate) found.
280  p = (y - x) / two
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
                         * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z))
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
9999 end subroutine hqr2

integer function lastnb(t)
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Find last non-blank in string
  !
  !----------------------------------------------------------------------*
  implicit none

  character(*) t
  integer i
  do i = len(t), 1, -1
     if (t(i:i) .ne. ' ') goto 20
  enddo
  i = 1
20 lastnb = i
end function lastnb
subroutine tmfoc(el,sk1,c,s,d,f)
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
  double precision c,d,el,f,qk,qkl,qkl2,s,sk1,zero,one,two,six,     &
       twelve,twty,thty,foty2
  parameter(zero=0d0,one=1d0,two=2d0,six=6d0,twelve=12d0,twty=20d0, &
       thty=30d0,foty2=42d0)

  !---- Initialize.
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

subroutine f77flush(i,option)
  implicit none
  integer i,ios
  real a
  logical ostat, fexist,option
  character(20) faccess,fform
  character(255) fname
  character(1) c
  inquire(err=5,iostat=ios,unit=i,opened=ostat,exist=fexist)
  if (.not.ostat.or..not.fexist) return
  inquire(err=6,iostat=ios,unit=i,access=faccess,form=fform,name=fname)
  close (unit=i,err=7,iostat=ios)
  !     write (*,*) 'Re-opening ',i,' ',faccess,fform,fname
  open(err=8,iostat=ios,unit=i,access=faccess,form=fform,file=fname,status='old')
  if (option) then
     if (fform.eq.'FORMATTED') then
3       read (i,100,err=9,iostat=ios,end=4) c
        go to 3
     else
2       read (i,err=10,iostat=ios,end=1) a
        go to 2
     endif
4    backspace i
1    continue
  endif
  return
100 format (a1)
5 write (*,*) ' F77FLUSH 1st INQUIRE FAILED with IOSTAT ',ios,' on UNIT ',i
  stop
6 write (*,*) ' F77FLUSH 2nd INQUIRE FAILED with IOSTAT ', ios,' on UNIT ',i
  stop
7 write (*,*) ' F77FLUSH CLOSE FAILED with IOSTAT ',ios,' on UNIT ',i
  stop
8 write (*,*) ' F77FLUSH RE-OPEN FAILED with IOSTAT ',ios,' on UNIT ',i
  stop
9 write (*,*) ' F77FLUSH FORMATTED READ FAILED with IOSTAT ',ios,' on UNIT ',i
  stop
10 write (*,*) ' F77FLUSH UNFORMATTED READ FAILED with IOSTAT ',ios,' on UNIT ',i
  stop
end subroutine f77flush


subroutine seterrorflag(errorcode,from,descr)
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Puts global error flag in c code.                                *
  !     Input:                                                           *
  !     errorcode                                                        *
  !     from - name of a routine where the error occured                 *
  !     descr - description of the error that has occured                *
  !     Input/output:                                                    *
  !----------------------------------------------------------------------*
  implicit none
  integer :: errorcode
  character(*) :: from
  character(*) :: descr
  integer  n,m

  n = LEN(from)
  m = LEN(descr)
  call seterrorflagfort(errorcode,from,n,descr,m)

end subroutine seterrorflag
