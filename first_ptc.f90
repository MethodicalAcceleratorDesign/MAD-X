subroutine ttwm()
  !--- PTC entry point
  USE madx_keywords
  implicit none
  integer iii,icav
  integer :: EXCEPTION=0
  real(dp) beta0
  real(kind(1d0)) get_value
  type(layout) LHC
  iii=0
  print77=.false.
  read77 =.false.
  print*,"Now PTC"

  call ptc_input(lhc,icav,EXCEPTION)
  if(EXCEPTION.eq.1) then
     call fort_warn('wrong magnet type KINDI which must be: ','1, 2, 3')
     return
  endif
  if(get_value('ptc ','ptc_twiss ') .ne. 0) call ptc_twiss(lhc,icav)
  if(get_value('ptc ','normal ') .ne. 0) call ptc_normal(lhc,icav)
  call kill(lhc)
  call kill_tpsa
  call nul_coef(SECTOR_B)
end subroutine ttwm
