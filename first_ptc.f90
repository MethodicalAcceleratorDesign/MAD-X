subroutine ttwm()
  !--- PTC entry point
  USE madx_keywords
  implicit none
  logical(lp) ok
  integer iii,kindi,icav,blah
  real(dp) beta0,get_value
  type(layout) LHC
  iii=0
  print77=.false.
  read77 =.false.
  print*,"Now PTC"
  kindi=get_value('ptc ','kindi ')
  print*,kindi
  if(kindi==2) kindi=kind2
  if(kindi==6) kindi=kind6
  if(kindi==7) kindi=kind7

  call ptc_input(lhc,icav)
  call ptc_twiss(lhc,icav)
  call kill(lhc)
  call kill_tpsa
  call nul_coef(SECTOR_B)
end subroutine ttwm
