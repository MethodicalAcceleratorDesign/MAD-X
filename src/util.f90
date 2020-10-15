module math_constfi ! 2015-Aug-06 Ghislain
  implicit none
  public
  !--- Definition of useful mathematical constants
  double precision, parameter :: zero=0d0, one=1d0, two=2d0, three=3d0, four=4d0, five=5d0
  double precision, parameter :: six=6d0, seven=7d0, eight=8d0, nine=9d0, ten=10d0
  double precision, parameter :: eleven=11d0, twelve=12d0, thirteen=13d0, fourteen=14d0, fifteen=15d0
  double precision, parameter :: sixteen=16d0, seventeen=17d0, eighteen=18d0, nineteen=19d0, twenty=20d0
  double precision, parameter :: half=0.5d0, quarter=0.25d0
  double precision, parameter :: ten3m=1d-3, ten6m=1d-6, ten9m=1d-9
  double precision, parameter :: ten3p=1d3,  ten6p=1d6,  ten9p=1d9
  double precision, parameter :: pi = 3.141592653589793238462643383279502884197169399375105820974944d0
  double precision, parameter :: twopi = two*pi
  double precision, parameter :: degrad = 180d0/pi, raddeg = pi/180d0
  double precision, parameter :: e = 2.718281828459045235360287471352662497757247093699959574966967d0
end module math_constfi

module phys_constfi
  use math_constfi, only : pi
  implicit none
  public
  !--- Definition of physical constants
  ! sources  :
  ! https://physics.nist.gov/cuu/pdf/wall_2018.pdf
  double precision, parameter :: clight = 299792458d0        ! Speed of light in vacuum [m/s]
  double precision, parameter :: qelect = 1.602176634d-19   ! Elementary charge [A*s]

  double precision, parameter :: hbar   = 6.582119569e-25    ! Reduced Plack's constant [GeV*s]

  double precision, parameter :: amu0   = 4d-7 * pi          ! Permeability of vacuum [V*s/A*m]
  ! Rest Mass [GeV]
  double precision, parameter ::  emass  = 0.51099895000d-3   ! GeV CODATA 2018

  double precision, parameter ::  pmass  = 0.93827208816     ! GeV CODATA 2018

  double precision, parameter ::  nmass  = 0.93956542052      ! GeV CODATA 2018
  double precision, parameter ::  mumass = 0.1056583755     ! GeV CODATA 2018

  ! Classical radius [m]
  double precision, parameter :: erad   = 2.8179403262d-15
  double precision, parameter :: prad   = erad*emass/pmass
end module phys_constfi

module matrices
  ! useful matrices: Identity EYE(6,6) and Symplectic JMAT(6.6) and SMAT(2,2)
  implicit none
  double precision, parameter :: EYE(6,6)=reshape((/1d0,0d0,0d0,0d0,0d0,0d0,&
                                                    0d0,1d0,0d0,0d0,0d0,0d0,&
                                                    0d0,0d0,1d0,0d0,0d0,0d0,&
                                                    0d0,0d0,0d0,1d0,0d0,0d0,&
                                                    0d0,0d0,0d0,0d0,1d0,0d0,&
                                                    0d0,0d0,0d0,0d0,0d0,1d0 /), shape(EYE))
  double precision, parameter :: JMAT(6,6)=reshape((/0d0,1d0,0d0,0d0,0d0,0d0, &
                                                    -1d0,0d0,0d0,0d0,0d0,0d0, &
                                                     0d0,0d0,0d0,1d0,0d0,0d0, &
                                                     0d0,0d0,-1d0,0d0,0d0,0d0, &
                                                     0d0,0d0,0d0,0d0,0d0,1d0, &
                                                     0d0,0d0,0d0,0d0,-1d0,0d0 /), shape(JMAT))
  double precision, parameter :: JMATINV(6,6) = -JMAT
  double precision, parameter :: JMATT(6,6)   = -JMAT
  double precision, parameter :: SMAT(2,2)=reshape((/0d0,1d0, &
                                                    -1d0,0d0 /), shape(SMAT))
  double precision, parameter :: SMATINV(2,2) = -SMAT
  double precision, parameter :: SMATT(2,2)   = -SMAT

  double precision, parameter :: symp_thrd = 1d-12 !, symp_thrd_orbit = 1d-10 ! threshold during closed orbit search
end module matrices

module code_constfi
  implicit none
  public
  !--- Definition of mad8 codes in a more readable format
  integer, parameter :: code_drift = 1
  integer, parameter :: code_rbend = 2
  integer, parameter :: code_sbend = 3
  integer, parameter :: code_matrix = 4
  integer, parameter :: code_quadrupole = 5
  integer, parameter :: code_sextupole = 6
  integer, parameter :: code_octupole = 7
  integer, parameter :: code_multipole = 8
  integer, parameter :: code_solenoid = 9
  integer, parameter :: code_rfcavity = 10
  integer, parameter :: code_elseparator = 11
  integer, parameter :: code_srotation = 12
  integer, parameter :: code_yrotation = 13
  integer, parameter :: code_hkicker = 14
  integer, parameter :: code_kicker = 15
  integer, parameter :: code_vkicker = 16
  integer, parameter :: code_hmonitor = 17
  integer, parameter :: code_monitor = 18
  integer, parameter :: code_vmonitor = 19
  integer, parameter :: code_ecollimator = 20
  integer, parameter :: code_rcollimator = 21
  integer, parameter :: code_beambeam = 22
  ! code 23 is free
  integer, parameter :: code_instrument = 24
  integer, parameter :: code_marker = 25
  ! code 26 is free
  integer, parameter :: code_twcavity = 27
  ! code 28 is free
  integer, parameter :: code_wire = 29
  integer, parameter :: code_slmonitor = 30
  integer, parameter :: code_blmonitor = 31
  integer, parameter :: code_imonitor = 32
  integer, parameter :: code_dipedge = 33
  integer, parameter :: code_xrotation = 34
  integer, parameter :: code_changeref = 35
  integer, parameter :: code_translation = 36
  integer, parameter :: code_crabcavity = 37
  integer, parameter :: code_placeholder = 38
  integer, parameter :: code_tkicker = 39
  integer, parameter :: code_hacdipole = 40
  integer, parameter :: code_vacdipole = 41
  integer, parameter :: code_nllens = 42
  integer, parameter :: code_rfmultipole = 43
  integer, parameter :: code_collimator = 44
  integer, parameter :: code_changerefp0 = 45
  integer, parameter :: code_sixmarker = 46
end module code_constfi

module aperture_enums
  implicit none
  public
 
  integer, parameter :: ap_circle = 0
  integer, parameter :: ap_ellipse = 1
  integer, parameter :: ap_rectangle = 2
  integer, parameter :: ap_lhcscreen = 3
  integer, parameter :: ap_rectcircle = 4
  integer, parameter :: ap_rectellipse = 5
  integer, parameter :: ap_racetrack = 6
  integer, parameter :: ap_octagon = 7
  integer, parameter :: ap_custom = 8
  integer, parameter :: ap_notset = 9
  integer, parameter :: ap_custom_inter = 10
end module aperture_enums
module track_enums
  implicit none
  public
  integer, parameter :: enum_other_bv = 1
  integer, parameter :: enum_lrad = 2
  integer, parameter :: enum_noise = 3
  integer, parameter :: enum_angle = 4
  integer, parameter :: enum_time_var = 5
  integer, parameter :: enum_sinkick = 6
  integer, parameter :: enum_kick = 7
  integer, parameter :: enum_chkick = 8
  integer, parameter :: enum_cvkick = 9
  integer, parameter :: enum_hkick = 10 
  integer, parameter :: enum_vkick = 11  
  integer, parameter :: total_enums = 11
end module track_enums
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
  integer, parameter :: sp=kind(1.e0)
  integer, parameter :: dp=selected_real_kind(2*precision(1.e0))

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

module deltrafi
  implicit none
  public
  logical :: radiate=.false., damp=.false., quantum=.false., fastune=.false.
  double precision :: deltax=0.d0
end module deltrafi

module dyntabfi
  implicit none
  public
  double precision :: dynapfrac=0.d0, dktrturns=0.d0
  double precision :: xend=0.d0, pxend=0.d0, yend=0.d0, pyend=0.d0, tend=0.d0, ptend=0.d0
  !double precision :: smear=0.d0, yapunov=0.d0
  double precision :: smear=0.d0, lyapunov=0.d0
end module dyntabfi

module wmaxmin0fi
  implicit none
  public
  double precision :: wxmax=0.d0,  wymax=0.d0,  wxymax=0.d0
  double precision :: wxmin=1.d20, wymin=1.d20, wxymin=1.d20
end module wmaxmin0fi

module tunesfi
  implicit none
  public
  double precision :: x0=0.d0, y0=0.d0, tunx=0.d0, tuny=0.d0, dtune=0.d0
end module tunesfi

module twiss0fi
  implicit none
  public
  !integer, parameter :: align_max=14, fundim=74
  !IT increase fundim to 110 to fit in sigma matrices 6x6
  integer, parameter :: align_max=14, fundim=110

end module twiss0fi

module twissafi
  implicit none
  public
  character(len=48) :: table_name=' ', sectorTableName=' '
  logical :: match_is_on=.false.
end module twissafi

module twisslfi
  implicit none
  public
  logical :: centre=.false., first
  logical :: rmatrix=.false., sectormap=.false., ripken=.false.
  logical :: mode_flip=.false.
  logical :: ele_body=.false.


end module twisslfi

module twisscfi
  use twiss0fi
  implicit none
  public
  double precision :: opt_fun0(fundim)=0.d0, opt_fun(fundim)=0.d0
  double precision :: disp(6)=0.d0, ddisp(6)=0.d0
  double precision :: rmat(2,2)=0.d0
  double precision :: sigmat(6,6)=0.d0
  double precision :: betx=0.d0, alfx=0.d0, amux=0.d0, cosmux=0.d0, sinmux=0.d0, qx=0.d0
  double precision :: bety=0.d0, alfy=0.d0, amuy=0.d0, cosmuy=0.d0, sinmuy=0.d0, qy=0.d0
  double precision :: bxmax=0.d0, dxmax=0.d0, bymax=0.d0, dymax=0.d0
  double precision :: xcomax=0.d0, ycomax=0.d0, sigxco=0.d0, sigyco=0.d0
  double precision :: sigdx=0.d0, sigdy=0.d0
  double precision :: wgt=0.d0, suml=0.d0, circ=0.d0, eta=0.d0, alfa=0.d0, gamtr=0.d0
  double precision :: wx=0.d0, phix=0.d0, dmux=0.d0, xix=0.d0, wy=0.d0, phiy=0.d0, dmuy=0.d0, xiy=0.d0
  double precision :: synch_1=0.d0, synch_2=0.d0, synch_3=0.d0, synch_4=0.d0, synch_5=0.d0
  double precision :: synch_6=0.d0, synch_7=0.d0, synch_8=0.d0
  double precision :: gammacp=1.d0
  integer :: nmode_flip=0
end module twisscfi

module twissotmfi
  implicit none
  public
  double precision :: rotm(6,6)=0.d0, rw(6,6)=0.d0, skick(6)=0.d0, sorb(6)=0.d0
  double precision :: srmat(6,6)=0.d0, stmat(6,6,6)=0.d0
end module twissotmfi

module twissbeamfi
  use math_constfi, only : zero
  implicit none
  public
  logical, save :: radiate=.false.
  double precision, save :: energy=zero, deltap=zero, beta=zero, gamma=zero, pc=zero
  double precision, save :: arad=zero, dtbyds=zero, charge=zero, npart=zero
end module twissbeamfi

module max_iterate
  implicit none
  public
  integer, parameter :: maxiter=250
end module max_iterate

module twiss_elpfi
  implicit none
  public
  !---fixed positions for element parameters
  double precision :: g_elpar(50)=0.d0
  !-general
  integer, parameter :: g_el=2, g_kmax=3, g_kmin=4, g_calib=5, g_polarity=6
  !-bend
  integer, parameter :: b_angle=7, b_tilt=8, b_k0=9, b_k0s=10
  integer, parameter :: b_k1=11, b_k1s=12, b_e1=13 , b_e2=14, b_k2=15
  integer, parameter :: b_k2s=16, b_h1=17, b_h2=18, b_hgap=19
  integer, parameter :: b_fint=20, b_fintx=21, b_k3=22, b_k3s=23
  !-quad
  integer, parameter :: q_tilt=7, q_k1=8, q_k1s=9, q_k1t=10, q_k1st=11
  !-sext
  integer, parameter :: s_tilt=7, s_k2=8, s_k2s=9, s_k2t=10, s_k2st=11
  !-oct
  integer, parameter :: o_tilt=7, o_k3=8, o_k3s=9
  !-mult
  integer, parameter :: m_tilt=7, m_lrad=8
  !-sol
  integer, parameter :: so_lrad=7, so_ks=8, so_ksi=9
  !-rfc
  integer, parameter :: r_volt=7, r_lag=8, r_freq=9, r_lagt=10
  !-elsep
  integer, parameter :: e_tilt=7, e_ex=8, e_ey=9
  !-hkick
  integer, parameter :: h_tilt=7, h_lrad=8, h_kick=9, h_hkick=10, h_chkick=11
  !-vkick
  integer, parameter :: v_tilt=7, v_lrad=8, v_kick=9, v_vkick=10, v_cvkick=11
  !-kick
  integer, parameter :: k_tilt=7, k_lrad=8, k_hkick=9, k_vkick=10, k_chkick=11, k_cvkick=12
end module twiss_elpfi

module emitfi
  implicit none
  public
  double precision, save :: qx=0.d0, qy=0.d0, qs=0.d0, cg=0.d0, sum(3)=0.d0, sumu0=0.d0
end module emitfi

module twtrrfi
  implicit none
  public
  !---- maxmul is the maximum multipole order both in twiss and trrun
  integer, parameter :: maxmul=20, maxferr=50, maxnaper=100
end module twtrrfi

module twtapering
  implicit none
  public
  !---- for tapering
  integer :: totrfcav =0, orderrun=0
  double precision ::  endpt=0
end module twtapering

module ibsdbfi
  implicit none
  public
  integer :: bunch=0
  double precision :: circ=0.d0, arad=0.d0, freq0=0.d0, alpha=0.d0
  double precision :: amass=0.d0, charge=0.d0, en0=0.d0, gammas=0.d0, gamma=0.d0
  double precision :: ex=0.d0, ey=0.d0, et=0.d0, sigt=0.d0, sige=0.d0, sigx=0.d0, sigy=0.d0
  double precision :: betas=0.d0, beta=0.d0
  double precision :: parnum=0.d0, currnt=0.d0, alfa=0.d0
end module ibsdbfi

module matchfi
  implicit none
  public
  integer :: icovar=0, ilevel=0
  double precision :: edm=0.d0, fmin=0.d0
end module matchfi

module name_lenfi
  implicit none
  public
  integer, parameter :: name_len=48
end module name_lenfi

module warncolim
  implicit none
  public
  integer :: warnede = 0
  integer :: warnedr = 0
end module warncolim

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
  double precision :: circ=0.d0, arad=0.d0, freq0=0.d0, amass=0.d0
  double precision :: charge=0.d0, en0=0.d0, gammas=0.d0, gamma=0.d0, deltap=0.d0
  double precision :: ex=0.d0, ey=0.d0, et=0.d0, sigt=0.d0, sige=0.d0
  double precision :: betas=0.d0, beta=0.d0, parnum=0.d0, currnt=0.d0, alfa=0.d0
  double precision :: um1=0.d0, fb1=0.d0, fb2=0.d0
end module touschekfi

module trackfi
  implicit none
  public
  double precision, save :: arad=0.d0, betas=0.d0, beti=0.d0, gammas=0.d0, dtbyds=0.d0
  double precision, save :: bet0=0.d0, bet0i=0.d0
  double precision :: deltas=0.d0, t_max=1.d20, pt_max=1.d20
  logical :: radiate=.false., damp=.false., quantum=.false., fsecarb=.false., thin_cf=.false.

  !! AL
  integer, save :: tot_turn = 0;
  integer, save :: jmax, segment
  logical :: lost_in_turn = .false., is_lost = .false.

  integer, parameter :: max_part=20000
  integer, save :: part_id_keep(max_part), last_turn_keep(max_part)
end module trackfi

module time_varfi
  use twtrrfi
  use name_lenfi
  implicit none
  public
  logical time_var_m, time_var_p, time_var_c
  integer, parameter :: n_time_var = 10000
  integer, save :: time_var_m_cnt, time_var_p_cnt, time_var_c_cnt
  integer, save :: time_var_m_lnt, time_var_p_lnt, time_var_c_lnt, trrun_nt
  double precision, save :: myfield(n_time_var,2,0:maxmul), phase_tromb(n_time_var,36), cav_volt(n_time_var)
  double precision, save :: time_var_m_ind(n_time_var), time_var_p_ind(n_time_var), time_var_c_ind(n_time_var)
  double precision, save :: time_var_m_nt(n_time_var), time_var_p_nt(n_time_var),  time_var_c_nt(n_time_var)
  character(len=name_len), save :: time_var_m_ch(n_time_var), time_var_p_ch(n_time_var), time_var_c_ch(n_time_var)
end module time_varfi

module plotfi
  implicit none
  public
  !--- m_adble is the number of different types of elements
  integer, parameter :: mtype = 50, m_adble = 20

  !--- mcnam adjusted to NAME_L
  integer, parameter :: mcnam = 48, maxpnt = 500

  !--- szcompar is the size of the arrays returned by the routine comm_para
  integer, parameter :: szcompar = 100

  !--- szchara is the size of the character strings char_a & version in
  !--- the routine pesopt
  integer, parameter :: szchara = 400

  !--- character sizes:
  !   MLSIZE    label character height
  !   MTSIZE    text   - " -
  !   MASIZE    annotation - " -
  integer, parameter :: mlsize = 13 ,mtsize = 13, masize = 20

  !--- parameters used in the routine pegacn in file plot.F

  integer, parameter :: mposx = 8, mposy = 3, mpost = mposx * mposy

  !--- parameters used in the routine peschm in file plot.F

  integer, parameter :: mobj = 14, msize = 88
  integer, parameter :: mtitl  = 128, mxlabl = 300 !hbu mxlabl = 160 increased to 300 as required by plot.f90
  integer, parameter :: mnvar = 74, mxdep = 2
  integer, parameter :: mqadd = 100000
  integer, parameter :: mntmax = 20, mksmax = 10
  integer, parameter :: mplred = 46, mplout = 47

  integer, parameter :: maxseql = 150000, mtwcol = 46, mpparm = 10
  integer, parameter :: mxcurv = 10, mopt = 60, mfile = 120, marg = 60
  integer, parameter :: maxarg = 1000, mxdp = 25, mxplot = 100

  integer, parameter :: mintpl = 18

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
  integer, save :: itbv=0, nivvar=0, nelmach=0, numax=0, interf=0, noline=0
  integer, save :: nqval(mxcurv)=0, nvvar(4)=0, nrrang(2)=0
  integer, save :: proc_flag(2,mxcurv)=0, ipparm(mpparm,mxcurv)=0
  integer, save :: naxref(mxcurv)=0, ieltyp(maxseql)=0

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
  real, save :: qascl=0., qlscl=0., qsscl=0., qtscl=0.
  real, save :: hrange(2)=0., vrange(2,4)=0., hmima(2)=0., vmima(2,4)=0.
  real, save :: qhval(maxseql,mxcurv)=0., qvval(maxseql,mxcurv)=0.
  real, save :: estart(maxseql)=0., eend(maxseql)=0.

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
  character(len=mfile), save  :: plfnam=' '
  character(len=mcnam), save  :: horname=' ', tabname=' ', sname(mxcurv)=' ', slabl(mxcurv)=' '
  character(len=mxlabl), save :: axlabel(4)=' '
  character(len=mtitl), save  :: toptitle=' '
end module plotfi

module plot_bfi
  implicit none
  public
  !--- Definition of common / peotcl /
  !--- fpmach set in routines pesopt and pefill, used in routine peplot
  !--- ddp_flag set in routine pefill, used in routine peplot
  !--- ptc_flag set in routines pesopt, used in routine pefill
  logical, save :: fpmach=.false., dpp_flag=.false., ptc_flag=.false.
end module plot_bfi

module resindexfi
  implicit none
  public
  integer, parameter :: mnres=1000, mymorder=20
end module resindexfi

module gxx11_common
  implicit none
  public
  real, parameter :: toleps=1.e-5, versio=1.50

  integer :: normt
  integer, parameter :: madim1=500, madim2=100, maxset=20, mconid=7, merrun=10
  integer, parameter :: metaun=11, miunit=5, mmetat=4, mounit=6, mpaxs=23, mpcurv=10
  integer, parameter :: mtermt=101, mtick=10, mtmeta=2, mtterm=1, mtfont=12
  integer, parameter :: mxaxs=4, mxpix=1000, mxsize=27, myaxs=4, mypix=1000, mysize=19
  integer, parameter :: mnormt=2, mx11pr=10, mx11tf=10, mxxpix=1200, mxypix=1000
  integer, parameter :: mcolor=6, mpspro=8, meppro=8, mdict=24, mlpro=68, mpsep=3, mepep=2
  integer, parameter :: mhead=4, mline=72, msfact=4, mlbb1=17, mlbb2=30, mubb1=573, mubb2=790
  integer, parameter :: mwid1 = mubb1 - mlbb1, mwid2 = mubb2 - mlbb2


  integer :: itermt=0, interm=0, inmeta=0, ierrun=0, imetun=0, inunit=0, iounit=0, ipage=0
  integer :: isfflg=0, isqflg=0, iwtflg=0, iclflg=0, inormt=0, ipseps=0, iepsop=0, itseop=0
  integer :: iepscf=0, imetps=0, ipctct=0, iczebr=0, idinit=0, ipstyp=0, iclear=0, istotx=0
  integer :: lmpict=0, ltermt=0, lnterm=0, lnmeta=0, lerrun=0, lmetun=0, lnunit=0, lounit=0
  integer :: lsfflg=0, lsqflg=0, lwtflg=0, lclflg=0, lnormt=0, lmetax=0, lmetay=0, lmetnm=0
  integer :: lerrnm=0, ldefnl=0, lerrop=0, lmetop=0, ltotin=0, lacttm=0, lpseps=0, lundef=0
  integer :: lttime=0, ldinit=0, ltseop=0

  integer :: ixapar(mpaxs,mxaxs)=0, iyapar(mpaxs,myaxs)=0, icvpar(mpcurv ,maxset)=0

  integer :: nxpix=0, nypix=0, lxpix=0, lypix=0, icucol=0, iorips=0
  integer :: iutlps=0, ibbox(4)=0, ix11pr(mx11pr)=0, ix11tf(mx11tf)=0, ix11op(mx11tf)=0

  real :: fxpix=0., fypix=0., rx11pr(mx11pr)=0., rgbcol(3,mcolor)=0.

  real :: xmetaf=0., ymetaf=0., xsterm=0., ysterm=0., wfact=0., wttime=0., wxfact=0., wyfact=0.
  real :: vpfacx=0., vpfacy=0., vptdef(4)=0., vploc(4)=0., actwnd(4)=0.
  real :: rangex(2,mxaxs)=0., rangey(2 ,myaxs)=0.
  real :: cvwnwd(4,maxset)=0., axwndx(2,maxset), axwndy(2,maxset)=0.

  real :: xp(madim2+1)=0., xvp(madim2+1)=0, yp(madim2+1)=0., yvp(madim2+1)=0

  real :: p(madim1,2)=0., s(madim1)=0., yy1d(madim1,2)=0., yy2d(madim1,2)=0.

  character(len=256) :: smetnm=" ", serrnm=" ", spsnam=" "
  character(len=300) :: sxtext(mxaxs)=" ", sytext(myaxs)=" "
  character(len=20)  :: sxform(mxaxs)=" ", syform(myaxs)=" ", stortx=" "
  character(len=1)   :: sdefnl=" "
  character(len=16)  :: colour(mcolor)=" "
  character(len=60)  :: pshead(mhead)=" "
  character(len=maxset) :: splotc=" "
end module gxx11_common

module gxx11_aux
  implicit none
  public
  !
  character(len=100) :: strloc
  integer, save, dimension(14) :: ivals=(/ 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 1, 0 /)
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

  real, save, dimension(14) :: rvals=(/ 0., 1., 0.01, 0., 1., 0., 1., 0., 1., 0., 1., 1., 1., 1. /)
  !   rvals(1-2)  chup vector
  !     3         character height
  !     4-7       window
  !     8-11      viewport
  !     12        character expansion factor
  !     13        line width scale factor
  !     14        marker scale factor
end module gxx11_aux

module fasterror
  implicit none
  logical :: fasterror_on = .false.
  integer :: kstep
  integer, parameter :: nx=490, ny=470, idim=(nx+2)*(ny+2)
  double precision :: hrecip, wtimag(idim), wtreal(idim)
end module fasterror

module SCdat
  ! Yuri Alexahin Oct 2017 Mathematica version
  ! Frank Schmidt Oct 2017 Fortran90 version
  ! Copyright Fermilab & CERN
  implicit none
  public
  logical eflag,sc_3d_kick,sc_3d_periodic,sc_3d_beamsize
  integer, parameter :: idim=6,isigmatfit=50000,nptot=500000
  integer, parameter :: Nintegrate=1000000
  integer, parameter :: mmax=77,mmum=15,kmaxo=20;
  double precision, parameter :: Pi = 4d0 * atan(1d0), EulerGamma =&
       & 0.57721566490153286060651209008240243104215933593992d0
  double precision, parameter :: dampf=0.85d0, dampeta=0.25d0, epsz=1d-15 !(* damping factor *)
  double precision, parameter :: alfa=0.25d0 !(* not optimized yet *)
  double precision, parameter :: efactor=2d0*alfa/(1d0 + 2d0*alfa)**4
  double precision, parameter :: dmax=6d0
  double precision, parameter :: um=.1d0;
  double precision, parameter :: dres=6.5d0;
  double precision work_1(idim),work_2(idim),work_3(idim)
  double precision fcore
  double precision, dimension(6,6) :: sr = reshape(&
        (/0d0,-1d0,0d0,0d0,0d0,0d0, 1d0,0d0, 0d0,0d0, 0d0,0d0,&
        0d0,0d0,0d0,-1d0,0d0,0d0, 0d0,0d0,1d0,0d0, 0d0,0d0,&
        0d0,0d0,0d0,0d0,0d0,-1d0, 0d0,0d0, 0d0,0d0,1d0,0d0/), shape(sr))
  double precision, dimension(kmaxo) :: lIntdfact=(/ 3.14159265358979d0,&
       4.71238898038469d0,15.7079632679490d0,82.4668071567321d0,&
       593.761011528471d0,5442.80927234432d0,60648.4461775510d0,&
       796010.856080356d0,12028608.4918809d0,205689205.211164d0,&
       3926793917.66768d0,82789905097.4935d0,1910536271480.62d0,&
       47899873663549.8d0,1.296489913826749d15,3.767923812058989d16,&
       1.170272807510086d18,3.868401780380561d19,1.355976624070239d21,&
       5.023893392180235d22 /)
  double precision, dimension(kmaxo) :: eIntdfact=(/ 6.28318530717959d0,&
       18.8495559215388d0,94.2477796076938d0,659.734457253857d0,&
       5937.61011528471d0,65313.7112681318d0,849078.246485713d0,&
       12736173.6972857d0,216514952.853857d0,4113784104.22328d0,&
       86389466188.6889d0,1986957722339.84d0,49673943058496.1d0,&
       1.341196462579395d15,3.889469741480246d16,1.205735619858876d18,&
       3.978927545534292d19,1.392624640937002d21,5.152711171466908d22,&
       2.009557356872094d24 /)
  double precision, DIMENSION(2*kmaxo, kmaxo+1) :: scbinom=reshape( (/ &
                                !1
       1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !2
       -1.d0,3.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !3
       -1.d0,6.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !4
       1.d0,-10.d0,5.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !5
       1.d0,-15.d0,15.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !6
       -1.d0,21.d0,-35.d0,7.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !7
       -1.d0,28.d0,-70.d0,28.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !8
       1.d0,-36.d0,126.d0,-84.d0,9.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !9
       1.d0,-45.d0,210.d0,-210.d0,45.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !10
       -1.d0,55.d0,-330.d0,462.d0,-165.d0,11.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !11
       -1.d0,66.d0,-495.d0,924.d0,-495.d0,66.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !12
       1.d0,-78.d0,715.d0,-1716.d0,1287.d0,-286.d0,13.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !13
       1.d0,-91.d0,1001.d0,-3003.d0,3003.d0,-1001.d0,91.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !14
       -1.d0,105.d0,-1365.d0,5005.d0,-6435.d0,3003.d0,-455.d0,15.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !15
       -1.d0,120.d0,-1820.d0,8008.d0,-12870.d0,8008.d0,-1820.d0,120.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !16
       1.d0,-136.d0,2380.d0,-12376.d0,24310.d0,-19448.d0,6188.d0,-680.d0,17.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,0.d0,&
                                !17
       1.d0,-153.d0,3060.d0,-18564.d0,43758.d0,-43758.d0,18564.d0,-3060.d0,153.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !18
       -1.d0,171.d0,-3876.d0,27132.d0,-75582.d0,92378.d0,-50388.d0,11628.d0,-969.d0,19.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
       0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !19
       -1.d0,190.d0,-4845.d0,38760.d0,-125970.d0,184756.d0,-125970.d0,38760.d0,-4845.d0,&
       190.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !20
       1.d0,-210.d0,5985.d0,-54264.d0,203490.d0,-352716.d0,293930.d0,-116280.d0,20349.d0,&
       -1330.d0,21.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !21
       1.d0,-231.d0,7315.d0,-74613.d0,319770.d0,-646646.d0,646646.d0,-319770.d0,74613.d0,&
       -7315.d0,231.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !22
       -1.d0,253.d0,-8855.d0,100947.d0,-490314.d0,1144066.d0,-1352078.d0,817190.d0,-245157.d0,&
       33649.d0,-1771.d0,23.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !23
       -1.d0,276.d0,-10626.d0,134596.d0,-735471.d0,1961256.d0,-2704156.d0,1961256.d0,-735471.d0,&
       134596.d0,-10626.d0,276.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !24
       1.d0,-300.d0,12650.d0,-177100.d0,1081575.d0,-3268760.d0,5200300.d0,-4457400.d0,2042975.d0,&
       -480700.d0,53130.d0,-2300.d0,25.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !25
       1.d0,-325.d0,14950.d0,-230230.d0,1562275.d0,-5311735.d0,9657700.d0,-9657700.d0,5311735.d0,&
       -1562275.d0,230230.d0,-14950.d0,325.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !26
       -1.d0,351.d0,-17550.d0,296010.d0,-2220075.d0,8436285.d0,-17383860.d0,20058300.d0,-13037895.d0,&
       4686825.d0,-888030.d0,80730.d0,-2925.d0,27.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !27
       -1.d0,378.d0,-20475.d0,376740.d0,-3108105.d0,13123110.d0,-30421755.d0,40116600.d0,-30421755.d0,&
       13123110.d0,-3108105.d0,376740.d0,-20475.d0,378.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !28
       1.d0,-406.d0,23751.d0,-475020.d0,4292145.d0,-20030010.d0,51895935.d0,-77558760.d0,67863915.d0,&
       -34597290.d0,10015005.d0,-1560780.d0,118755.d0,-3654.d0,29.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !29
       1.d0,-435.d0,27405.d0,-593775.d0,5852925.d0,-30045015.d0,86493225.d0,-145422675.d0,145422675.d0,&
       -86493225.d0,30045015.d0,-5852925.d0,593775.d0,-27405.d0,435.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !30
       -1.d0,465.d0,-31465.d0,736281.d0,-7888725.d0,44352165.d0,-141120525.d0,265182525.d0,-300540195.d0,&
       206253075.d0,-84672315.d0,20160075.d0,-2629575.d0,169911.d0,-4495.d0,31.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                                !31
       -1.d0,496.d0,-35960.d0,906192.d0,-10518300.d0,64512240.d0,-225792840.d0,471435600.d0,-601080390.d0,&
       471435600.d0,-225792840.d0,64512240.d0,-10518300.d0,906192.d0,-35960.d0,496.d0,-1.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !32
       1.d0,-528.d0,40920.d0,-1107568.d0,13884156.d0,-92561040.d0,354817320.d0,-818809200.d0,1166803110.d0,&
       -1037158320.d0,573166440.d0,-193536720.d0,38567100.d0,-4272048.d0,237336.d0,-5456.d0,33.d0,0.d0,&
       0.d0,0.d0,0.d0,&
                                !33
       1.d0,-561.d0,46376.d0,-1344904.d0,18156204.d0,-131128140.d0,548354040.d0,-1391975640.d0,2203961430.d0,&
       -2203961430.d0,1391975640.d0,-548354040.d0,131128140.d0,-18156204.d0,1344904.d0,-46376.d0,561.d0,-1.d0,&
       0.d0,0.d0,0.d0,&
                                !34
       -1.d0,595.d0,-52360.d0,1623160.d0,-23535820.d0,183579396.d0,-834451800.d0,2319959400.d0,-4059928950.d0,&
       4537567650.d0,-3247943160.d0,1476337800.d0,-417225900.d0,70607460.d0,-6724520.d0,324632.d0,-6545.d0,35.d0,&
       0.d0,0.d0,0.d0,&
                                !35
       -1.d0,630.d0,-58905.d0,1947792.d0,-30260340.d0,254186856.d0,-1251677700.d0,3796297200.d0,-7307872110.d0,&
       9075135300.d0,-7307872110.d0,3796297200.d0,-1251677700.d0,254186856.d0,-30260340.d0,1947792.d0,&
       -58905.d0,630.d0,-1.d0,0.d0,0.d0,&
                                !36
       1.d0,-666.d0,66045.d0,-2324784.d0,38608020.d0,-348330136.d0,1852482996.d0,-6107086800.d0,12875774670.d0,&
       -17672631900.d0,15905368710.d0,-9364199760.d0,3562467300.d0,-854992152.d0,124403620.d0,-10295472.d0,&
       435897.d0,-7770.d0,37.d0,0.d0,0.d0,&
                                !37
       1.d0,-703.d0,73815.d0,-2760681.d0,48903492.d0,-472733756.d0,2707475148.d0,-9669554100.d0,22239974430.d0,&
       -33578000610.d0,33578000610.d0,-22239974430.d0,9669554100.d0,-2707475148.d0,472733756.d0,-48903492.d0,&
       2760681.d0,-73815.d0,703.d0,-1.d0,0.d0,&
                                !38
       -1.d0,741.d0,-82251.d0,3262623.d0,-61523748.d0,635745396.d0,-3910797436.d0,15084504396.d0,-37711260990.d0,&
       62359143990.d0,-68923264410.d0,51021117810.d0,-25140840660.d0,8122425444.d0,-1676056044.d0,211915132.d0,&
       -15380937.d0,575757.d0,-9139.d0,39.d0,0.d0,&
                                !39
       -1.d0,780.d0,-91390.d0,3838380.d0,-76904685.d0,847660528.d0,-5586853480.d0,23206929840.d0,-62852101650.d0,&
       113380261800.d0,-137846528820.d0,113380261800.d0,-62852101650.d0,23206929840.d0,-5586853480.d0,847660528.d0,&
       -76904685.d0,3838380.d0,-91390.d0,780.d0,-1.d0,&
                                !40
       1.d0,-820.d0,101270.d0,-4496388.d0,95548245.d0,-1121099408.d0,7898654920.d0,-35240152720.d0,103077446706.d0,&
       -202112640600.d0,269128937220.d0,-244662670200.d0,151584480450.d0,-63432274896.d0,17620076360.d0,-3159461968.d0,&
       350343565.d0,-22481940d0,749398.d0,-10660.d0,41.d0 /), &
       shape(scbinom), order=(/2,1/) )
  save sc_3d_kick,sc_3d_periodic,sc_3d_beamsize
end module SCdat

! SUBROUTINES


subroutine fort_info(t1, t2)
  implicit none
  character(*) :: t1, t2
  integer, external :: get_option

  if (get_option('info ') .ne. 0 .and. get_option('warn ') .ne. 0)  &
       print '(a,1x,a,1x,a)', '++++++ info:', t1, t2
end subroutine fort_info

subroutine fort_warn(t1, t2)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Print warning message.                                             *
  ! Input:                                                               *
  !   T1      (char)    usually calling routine name or feature          *
  !   T2      (char)    Message.                                         *
  !----------------------------------------------------------------------*
  character(*) :: t1, t2
  integer, external :: get_option

  if (get_option('warn ') .ne. 0) then
     print '(a,1x,a,1x,a)', '++++++ warning:', t1, t2
     call augmentfwarn()
  endif
end subroutine fort_warn

subroutine fort_fail(t1,t2)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Print fatal error message.                                         *
  ! Input:                                                               *
  !   T1      (char)    usually calling routine name or feature          *
  !   T2      (char)    Message                                          *
  !----------------------------------------------------------------------*
  character(*) :: t1, t2
  integer, external :: get_option

  print *,' '
  print *,  '+-+-+- fatal: ',t1,t2
  print *,' '

  if (get_option('no_fatal_stop ') .eq. 0) stop  1
end subroutine fort_fail

subroutine aafail(rout,text)
  implicit none
  character(*) :: rout, text
  call fort_fail(rout,text)
end subroutine aafail

logical function m66sta(amat)
  use math_constfi , only : one
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
  double precision, intent(IN) :: amat(6,6)

  integer :: j
  double precision, parameter :: tol=1d-12

  m66sta = abs(amat(6,6) - one) .le. tol
  do j = 1, 5
     m66sta = m66sta .and. abs(amat(6,j)) .le. tol
  enddo

end function m66sta

subroutine dcopy(in,out,n)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Copy arrays.                                                       *
  ! Input:                                                               *
  !   in  (double)    array to be copied.                                *
  !   n   (integer)   array length.                                      *
  ! Output:                                                              *
  !   out (double)    target array.                                      *
  !----------------------------------------------------------------------*
  integer, intent(IN) :: n
  double precision, intent(IN)  :: in(*)
  double precision, intent(OUT) :: out(*)

  OUT(1:n) = IN(1:n)

end subroutine dcopy

subroutine solver(augmat,ndim,mdim,irank)
  use math_constfi, only : zero
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
  integer, intent(IN)  :: ndim, mdim
  integer, intent(OUT) :: irank
  double precision, intent(IN OUT) :: augmat(ndim,ndim+mdim)

  double precision :: h, pivot
  integer :: ic, ip, ir, it, nc, nr

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

     if (pivot .eq. zero) return
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

end subroutine solver

subroutine symsol(a,n,eflag,work_1,work_2,work_3)
  use math_constfi, only : zero, one
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
  integer, intent(IN)  :: n
  double precision, intent(IN OUT) :: a(n,n)
  logical, intent(OUT) :: eflag

  integer :: i, j, k
  double precision :: si, work_1(n), work_2(n), work_3(n)

  !---- Scale upper triangle.
  eflag = .true.
  do i = 1, n
     si = a(i,i)
     if (si .le. zero) return
     work_1(i) = one / sqrt(si)
  enddo
  do i = 1, n
     do j = i, n
        a(i,j) = a(i,j) * work_1(i) * work_1(j)
     enddo
  enddo

  !---- Invert upper triangle.
  do i = 1, n
     if (a(i,i) .eq. zero) return
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

end subroutine symsol

subroutine symeig(a,nd,n,eigen,nval,work)
  use math_constfi, only : zero, one, two, four
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
  integer, intent(IN)  :: nd, n
  integer, intent(OUT) :: nval
  double precision, intent(IN OUT)  :: a(nd,nd)
  double precision, intent(OUT) :: eigen(nd), work(nd)

  integer :: i, it, j, k, l, m
  double precision :: b, c, f, g, h, p, r, s

  double precision, parameter :: big=1d10, eps=1d-20
  integer, parameter :: itmax=15

  nval = n
  if (n .le. 0) return

  !---- Matrix is 1 * 1.
  if (n .eq. 1) then
     eigen(1) = a(1,1)
     return
  endif

  !---- Matrix is 2 * 2.
  if (n .eq. 2) then
     f = a(1,1) + a(2,2)
     g = sqrt((a(1,1) - a(2,2))**2 + four * a(2,1)**2)
     eigen(1) = (f - g) / two
     eigen(2) = (f + g) / two
     return
  endif

  !---- n is at least 3, reduce to tridiagonal form.
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
        return
     endif

170  p = eigen(l) + f

     do i = l, 2, -1
        if (p .ge. eigen(i-1)) go to 190
        eigen(i) = eigen(i-1)
     enddo
     i = 1
190  eigen(i) = p

  enddo

end subroutine symeig

double precision function proxim(x,y)
  use math_constfi, only : twopi
  !----------------------------------------------------------------------*
  !   Proximity function of x and y.                                     *
  !   If angle is larger than pi between vector x and y, 2pi is added to *
  !   to this angle                                                      *
  !----------------------------------------------------------------------*
  implicit none
  double precision :: x, y

  proxim = x + twopi*anint((y-x)/twopi)

end function proxim

character(len=48) function charconv(tint)
  implicit none
  !----------------------------------------------------------------------*
  ! purpose:                                                             *
  !   converts integer array to string (based on ascii)                  *
  ! input:                                                               *
  !   tint  (int array)  1 = length, rest = string                       *
  !----------------------------------------------------------------------*
  integer :: tint(*)

  integer :: i, j, n
  integer, parameter :: m = 128
  character(len=m) :: letter
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
  use matrices, only : EYE
  use math_constfi, only : zero, one
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
  double precision, intent(IN)  :: fm(6,6)
  double precision, intent(OUT) :: reeig(6), aieig(6), am(6,6)

  integer :: i, info, ipind, iqind, j, k, kpnt(6)
  double precision :: aival(6), d(6), reval(6), tm(6,6)
  double precision :: big, c, dx, dy, pb, s

  integer, parameter :: ilo=1, ihi=4, mdim=6, nn=4

  !---- Compute eigenvalues and vectors.
  TM = FM
  AM = EYE
  call orthes(mdim,nn,ilo,ihi,tm,d)
  call ortran(mdim,nn,ilo,ihi,tm,d,am)
  call hqr2(mdim,nn,ilo,ihi,tm,reval,aival,am,info)

  if (info .ne. 0) then
     write (6, 910) ((fm(i,k), k = 1, 6), i = 1, 6)
910  format('Unable to find eigenvalues for matrix:'/(6f12.6))
     call fort_fail('LASEIG',' Unable to find eigenvalues for matrix')
     return
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
  TM = AM
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
  TM = EYE
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

  AM = matmul(AM,TM)

end subroutine laseig

subroutine ladeig(fm,reeig,aieig,am)
  use matrices, only : EYE
  use math_constfi, only : zero
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
  double precision, intent(IN)  :: fm(6,6)
  double precision, intent(OUT) :: reeig(6), aieig(6), am(6,6)

  integer i, info, j, k, kpnt(6)
  double precision :: aival(6), d(6), reval(6), tm(6,6)
  double precision :: big, c, dt, dx, dy, pb, s

  integer, parameter :: ilo=1, ihi=6, mdim=6, nn=6

  !---- Compute eigenvalues and eigenvectors.
  TM = FM
  call orthes(mdim,nn,ilo,ihi,tm,d)
  call ortran(mdim,nn,ilo,ihi,tm,d,am)
  call hqr2(mdim,nn,ilo,ihi,tm,reval,aival,am,info)
  if (info .ne. 0) then
     write (6, 910) ((fm(i,k), k = 1, 6), i = 1, 6)
910  format('Unable to find eigenvalues for matrix:'/(6f12.6))
     call fort_fail('LADEIG',' Unable to find eigenvalues for matrix')
     return
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
  TM = AM
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
  TM = EYE
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
  AM = matmul(AM,TM)

end subroutine ladeig

subroutine orthes(ndim,n,ilow,iupp,a,d)
  use math_constfi, only : zero
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
  integer, intent (IN) :: ndim, n, ilow, iupp
  double precision, intent (IN OUT) :: a(ndim,n)
  double precision, intent (OUT) :: d(n)

  integer :: i, j, m
  double precision :: f, g, h, scale

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
  use math_constfi, only : zero, one
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
  integer, intent(IN) :: ndim, n, ilow, iupp
  double precision, intent(IN)  :: h(ndim,n)
  double precision, intent(OUT) :: d(n), v(ndim,n)

  integer :: i, j, k, m
  double precision :: x, y

  !---- Initialize V to identity matrix.
  V(:n,:n) = zero
  do i = 1, n
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
  use math_constfi, only : zero, one, two
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
  integer, intent(IN)  :: ndim, n, ilow, iupp
  integer, intent(OUT) :: ierr
  double precision, intent(OUT) :: h(ndim,n), wi(n), wr(n), vecs(ndim,n)

  integer :: i, j, k, l, m, na, ien, its
  double precision :: den, hnorm, p, q, r, ra, s, sa, t, temp, tempi
  double precision :: tempr, vi, vr, w, x, y, z

  double precision, parameter :: epsmch=1d-16, triqua=.75d0, fac1=.4375d0

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
        return
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

end subroutine hqr2

integer function lastnb(t)
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Find last non-blank in string
  !----------------------------------------------------------------------*
  implicit none

  character(*) :: t
  integer :: i

  do i = len(t), 1, -1
     if (t(i:i) .ne. ' ') goto 20
  enddo
  i = 1

20 lastnb = i
end function lastnb

subroutine f77flush(i, option)
  implicit none
  integer :: i, ios
  logical :: ostat, fexist, option
  character(len=20) :: faccess, fform
  character(len=255) :: fname
  character(len=1) :: c
  real :: a ! dp ?

  inquire( err=5, iostat=ios, unit=i, opened=ostat, exist=fexist)
  if (.not.ostat .or. .not.fexist) return

  inquire(err=6, iostat=ios, unit=i, access=faccess, form=fform, name=fname)

  close (unit=i, err=7, iostat=ios)

  open(err=8, iostat=ios, unit=i, access=faccess, form=fform, file=fname, status='old')

  if (option) then
     if (fform .eq. 'FORMATTED') then
3       read (i, 100, err=9, iostat=ios, end=4) c
        go to 3
     else
2       read (i, err=10, iostat=ios, end=1) a
        go to 2
     endif
4    backspace i
1    continue
  endif

  return

100 format (a1)

5  write (*,*) ' F77FLUSH 1st INQUIRE FAILED with IOSTAT ',ios,' on UNIT ',i ;  stop
6  write (*,*) ' F77FLUSH 2nd INQUIRE FAILED with IOSTAT ', ios,' on UNIT ',i ;  stop
7  write (*,*) ' F77FLUSH CLOSE FAILED with IOSTAT ',ios,' on UNIT ',i ; stop
8  write (*,*) ' F77FLUSH RE-OPEN FAILED with IOSTAT ',ios,' on UNIT ',i ; stop
9  write (*,*) ' F77FLUSH FORMATTED READ FAILED with IOSTAT ',ios,' on UNIT ',i ; stop
10 write (*,*) ' F77FLUSH UNFORMATTED READ FAILED with IOSTAT ',ios,' on UNIT ',i ;  stop

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
  character(*) :: from, descr

  integer :: n,m

  n = LEN(from)
  m = LEN(descr)
  call seterrorflagfort(errorcode,from,n,descr,m)

end subroutine seterrorflag

integer function get_file_unit (lu_max)
  !
  !   get_file_unit returns a unit number that is not in use
  integer lu_max,  lu, m, iostat
  logical   opened
  !
  m = lu_max  ;  if (m < 1) m = 97
  do lu = m,1,-1
     inquire (unit=lu, opened=opened, iostat=iostat)
     if (iostat.ne.0) cycle
     if (.not.opened) exit
  end do
  !
  get_file_unit = lu
  return
end function get_file_unit

double precision function lInt(n,r,u,v)
  use SCdat
  implicit none
  integer n,i,ii,i1,i2,j
  double precision r,u,v,result
  !********************************************************************************
  result=0.d0
  SELECT CASE (n)
  CASE (1:20)
     i=2*n
     ii=i-1
     i1=i
     i2=0
     result=result+scbinom(ii,1)*u**i
     do j=2,n
        i1=i1-2
        i2=i2+2
        result=result+scbinom(ii,j)*u**i1*v**i2
     enddo
     result=result+scbinom(ii,n+1)*v**i
     lInt=lIntdfact(n)*(-1d0+r**2)**n*result
  CASE DEFAULT
     call aafail('lInt: Fatal: ',                                     &
          'Out of range in function lInt: Program stops')
     stop
  END SELECT
end function lInt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function eInt(n,r,u,v)
  use SCdat
  implicit none
  integer n,i,i1,i2,j
  double precision r,u,v,result
  !********************************************************************************
  result=0.d0
  SELECT CASE (n)
  CASE (1:20)
     i=2*n
     i1=i
     i2=0
     result=result+scbinom(i,1)*u**i
     do j=2,n
        i1=i1-2
        i2=i2+2
        result=result+scbinom(i,j)*u**i1*v**i2
     enddo
     result=result+scbinom(i,n+1)*v**i
     eInt=eIntdfact(n)*(-1d0+r**2)**n*u*result
  CASE DEFAULT
     call aafail('eInt: Fatal: ',                                     &
          'Out of range in function eInt: Program stops')
     stop
  END SELECT
end function eInt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function bips(a, mmax1)

  implicit none
  integer mmax1
  double precision a,r,u,v
  !********************************************************************************
  bips=(8388608d0/dble(1+mmax1)+a*(-(4194304d0/dble(2+mmax1))+&
       &a*(3145728d0/dble(3+mmax1)+a*(-(2621440d0/dble(4+mmax1))+&
       &a*(2293760d0/dble(5+mmax1)+a*(-(2064384d0/dble(6+mmax1))+&
       &a*(1892352d0/dble(7+mmax1)+a*(-(1757184d0/dble(8+mmax1))+&
       &a*(1647360d0/dble(9+mmax1)+17d0*a*(-(91520d0/dble(10+mmax1))+19d0*&
       &a*(4576d0/dble(11+mmax1)+7d0*a*(-(624d0/dble(12+mmax1))+(598d0*&
       &a)/dble(13+mmax1)-(575d0*a**2)/dble(14+mmax1)))))))))))))/8388608d0;
end function bips
