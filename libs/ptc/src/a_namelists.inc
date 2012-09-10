  TYPE fibrelist
     real(dp) GAMMA0I_GAMBET_MASS_AG(4)  !GAMMA0I,GAMBET,MASS ,AG  BETA0 is computed
     integer  DIR_CHARGE(2)  !DIR,CHARGE
!	 integer pos  ! layout position
!	 integer loc  ! tied universed
  END TYPE fibrelist
type(fibrelist) fib0
namelist / fibrename / fib0

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
    real(dp) VOLT_FREQ_PHAS_LAG(4) 
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
namelist / tp100_list / tp100

  TYPE k16_list
   logical(lp) DRIFTKICK,LIKEMAD
  END TYPE k16_list
type(k16_list) k160
namelist / k160_list / k160

