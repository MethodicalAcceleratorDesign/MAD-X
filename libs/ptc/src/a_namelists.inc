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
     real(dp) A_D(3),B_D(3),A_ANG(3),B_ANG(3), A_T,B_T 
     INTEGER(2)  ENERGY,TIME
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
    LOGICAL(Lp) KIN_KEX_BENDFRINGE_permFRINGE_EXACT(5)  ! KILL_ENT_FRINGE, KILL_EXI_FRINGE, bend_fringe,permFRINGE,EXACT
	INTEGER METHOD_NST_NMUL(3)   ! METHOD,NST,NMUL
  END TYPE MAGNET_CHARTLIST
type(MAGNET_CHARTLIST) MAGL0
namelist / MAGLname / MAGL0

  integer, private, parameter :: nmul_max=20  
  TYPE ele_list
    INTEGER KIND
    character(nlp) name_vorname(2)  
	real(dp) L,B_SOL
	real(dp) an(nmul_max)  
	real(dp) bn(nmul_max)  
    real(dp) VOLT_FREQ_PHAS(3)   ! alg not used 
	LOGICAL(LP) THIN 
    real(dp) fint_hgap_h1_h2(4)
	logical(lp) slowac_recut_even_electric_MIS(5)  
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

  TYPE thin3_list
   real(dp) thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle,hf,vf,ls
   logical(lp) patch
  END TYPE thin3_list
type(thin3_list) thin30
namelist / thin30name / thin30

  TYPE tp10_list
   logical(lp) DRIFTKICK
  END TYPE tp10_list
type(tp10_list) tp100
namelist / tp100name / tp100

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