subroutine ttwm()
  !--- Fortran main program to call c modules and ptc Fortran90 module
  USE madx_keywords
  implicit none
  logical(lp) ok
  integer iii,kindi,icav,blah
  real(dp) beta0
  type(layout) LHC!,lhc1
  iii=0
  print77=.false.
  read77 =.false.
  print*,"Now PTC"
  OK=.TRUE.
  DO WHILE(OK)
     WRITE(6,*)  " Give magnet kind type 2, 6 or 7 "
     READ(98,*)   kindi
     if(kindi==2) kindi=kind2
     if(kindi==6) kindi=kind6
     if(kindi==7) kindi=kind7
     if(kindi==kind2.or.kindi==kind6.or.kindi==kind7) ok=.false.
  enddo
  !  do iii=1,1000

  call ptc_input(lhc,icav)
  !     call kill(lhc)
  !     rewind 98
  !     READ(98,*) blah
  !  enddo
  !  lhc1=0
  !  call copy(lhc,lhc1)
  !  call ptc_twiss(lhc1,icav)
  !   do iii=1,1000
!  do iii=1,1000
!     rewind 98
!     READ(98,*)   blah
!     READ(98,*)   blah
!     READ(98,*)   blah
!     READ(98,*)   blah
     call ptc_twiss(lhc,icav)
!  enddo
  rewind 98
  READ(98,*)   blah
  READ(98,*)   blah
  READ(98,*)   blah
  rewind 16
  rewind 20
  !   enddo
  call kill(lhc)
  call kill_tpsa
  call nul_coef(SECTOR_B)
  !  call kill(lhc1)
  !  print*,"kill(lhc1)"
  !  isleep=usleep(10)
  !  do isleep=1,10000000
  !     blah=sqrt(float(isleep/isleep)+float(isleep/isleep))
  !  enddo
end subroutine ttwm
