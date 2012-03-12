TYPE fibrelist
   real(dp) t(4)  !GAMMA0I,GAMBET,MASS ,AG  BETA0 is computed
   integer  i(2)  !DIR,CHARGE
   !      integer pos  ! layout position
   !      integer loc  ! tied universed
END TYPE fibrelist
type(fibrelist) fib0
namelist / fibrename / fib0




!TYPE MAGNET_CHART
!   type(magnet_frame), pointer:: f
!   type(MADX_APERTURE), pointer:: APERTURE
!   type(S_APERTURE), pointer:: A(:)s
!   integer,pointer :: charge  ! propagator
!   integer,pointer :: dir    ! propagator
!   real(dp), POINTER :: BETA0,GAMMA0I,GAMBET,MASS,ag
!   real(dp), POINTER :: P0C
!   real(dp), POINTER :: LD,B0,LC         !
!   real(dp), POINTER :: TILTD      ! INTERNAL FRAME
!   real(dp),  DIMENSION(:), POINTER :: EDGE         ! INTERNAL FRAME!

!
!   LOGICAL(lp), POINTER :: EXACT
!   LOGICAL(lp), POINTER :: KILL_ENT_FRINGE,KILL_EXI_FRINGE, bend_fringe,permFRINGE                  !
!
!   INTEGER, POINTER :: METHOD,NST                   ! METHOD OF INTEGRATION 2,4,OR 6 YOSHIDA
!   INTEGER, POINTER :: NMUL                         ! NUMBER OF MULTIPOLE

!END TYPE MAGNET_CHART
