  TYPE fibrelist
     real(dp) GAMMA0I_GAMBET_MASS_AG(4)  !GAMMA0I,GAMBET,MASS ,AG  BETA0 is computed
     real(dp)  CHARGE
     integer  DIR !DIR,CHARGE
     integer patch
  END TYPE fibrelist
type(fibrelist) fib0
namelist / fibrename / fib0

 TYPE patchlist
     integer  A_X1,A_X2,B_X1,B_X2
     real(dp) A_D(3),B_D(3),A_ANG(3),B_ANG(3), A_L,B_L , A_T,B_T 
     INTEGER(2)  ENERGY,TIME,GEOMETRY
  END TYPE patchlist
type(patchlist) patch0
namelist / patchname / patch0

 TYPE CHARTlist
     real(dp) D_IN(3),D_OUT(3),ANG_IN(3),ANG_OUT(3) 
  END TYPE CHARTlist
type(CHARTlist) CHART0
namelist / CHARTname / CHART0

  TYPE MAGNET_CHARTLIST
    real(dp) LC_LD_B0_P0(4)   ! LC LD B0 P0C
    real(dp) TILTD_EDGE(3)   ! TILTD EDGE
    LOGICAL(Lp) KIN_KEX_BENDFRINGE_EXACT(4)  ! KILL_ENT_FRINGE, KILL_EXI_FRINGE, bend_fringe,EXACT
	INTEGER METHOD_NST_NMUL_permfringe_highest(5)   ! METHOD,NST,NMUL,permfringr, highest_fringe
    LOGICAL(Lp) kill_spin(2) 
  END TYPE MAGNET_CHARTLIST

type(MAGNET_CHARTLIST) MAGL0
namelist / MAGLname / MAGL0

  integer, private, parameter :: nmul_max=NMAX
  TYPE ele_list
    INTEGER KIND
!    character(nlp+3) name_vorname(2)  
    character(nlp) name_vorname(2)  
	real(dp) L,B_SOL
	real(dp) an(nmul_max)  
	real(dp) bn(nmul_max)  
    real(dp) VOLT_FREQ_PHAS(3)   ! alg not used 
	LOGICAL(LP) THIN 
    real(dp) fint_hgap_h1_h2_va_vs(8)
	logical(lp) recut_even_electric_MIS(4) 
	integer(2) slow_ac 
     logical(lp) usebf_do1bf(4)
	 integer skipptcbf(2)
     character(120) filef
     character(120) fileb
  END TYPE ele_list
type(ele_list) ELE0
namelist / ELEname / ELE0
     
  integer, private, parameter :: harm_max=20  
  TYPE cav_list
   integer N_BESSEL,NF,CAVITY_TOTALPATH
   real(dp) phase0,t
   logical always_on
   real(dp) F(harm_max)
   real(dp) PH(harm_max)   
   real(dp) A,R
  END TYPE cav_list
type(cav_list) cav0
namelist / CAVname / cav0
 
  TYPE hel_list
   integer N_BESSEL
   real(dp) fake_shift(6)
  END TYPE hel_list
type(hel_list) hel0
namelist / HELname / hel0

  TYPE thin3_list
   real(dp) thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle,hf,vf,ls
   logical(lp) patch
   real(dp) dx_dy_pitchx_pitchy(4)
  END TYPE thin3_list
type(thin3_list) thin30
namelist / thin30name / thin30


  TYPE sol5_list
   real(dp) dx_dy_pitchx_pitchy(4)
  END TYPE sol5_list
type(sol5_list) sol50
namelist / sol50name / sol50

  TYPE tp10_list
   logical(lp) DRIFTKICK
	real(dp) ae(nmul_max)  
	real(dp) be(nmul_max)  
  END TYPE tp10_list
type(tp10_list) tp100
namelist / tp100name / tp100

  TYPE ab_list
    integer n_m(2)
	complex(dp) b(nmul_max,nmul_max) 
	complex(dp) E(nmul_max,nmul_max)  
	real(dp) dz_t_te(3*nmul_max+3)  
	real(dp) scale_angc_xc_vc_dc_hc(6)
  END TYPE ab_list
type(ab_list) ab0
namelist / ab0name / ab0
 
  TYPE k16_list
   logical(lp) DRIFTKICK,LIKEMAD
  END TYPE k16_list
type(k16_list) k160
namelist / k160name / k160

  TYPE ap_list
   integer KIND
   real(dp) r(2)
   real(dp) x,y,dx,dy
   logical on
  END TYPE ap_list
 type(ap_list) aplist
 namelist / aperturename / aplist

  TYPE tcav_list   
   real(dp) PSI_DPHAS_DVDS(3)
   logical always_on
  END TYPE tcav_list
type(tcav_list) tcav0
namelist / tCAVname / tcav0

  TYPE siam_list   
   real(dp) angle(3)
   real(dp) d(3)
  END TYPE siam_list
type(siam_list) siam0
namelist / siamname / siam0

  TYPE track_list   
   character(120) name
   integer range(3),direction
  END TYPE track_list
 type(track_list) track0
 namelist / trackname / track0

  integer, private, parameter :: mode_max=200  

  TYPE wig_list
   real(dp) offset
   real(dp) internal(6)
   real(dp) ex(wiggler_suntao),ey(wiggler_suntao)
   integer n
   real(dp) a(mode_max)
   real(dp) f(mode_max)   
   integer form(mode_max)
   real(dp) k(3,mode_max)
   integer ne
   real(dp) ae(mode_max)
   real(dp) fe(mode_max)   
   integer forme(mode_max)
   real(dp) ke(3,mode_max)
  END TYPE wig_list
type(wig_list) wig0
namelist / wigname / wig0

