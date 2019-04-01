subroutine touschek
  !--- 2013-Nov-26  15:41:21  ghislain:
  !    rewrite of the original touschek routine (F.Zimmermann, C. Milardi) with different logic:
  !    instead of calculating an exact inverse lifetime based on the average of
  !    beam parameters (beta,alfa,dispersion,derivative of dispersion) across the element,
  !    we either (w/o centre flag in twiss) calculate the exact inverse lifetime
  !      at entrance and exit of element, based on exact beam parameters at both locations ,
  !      and then only take the average inverse lifetime between entrance and exit,
  !    or (w centre flag in twiss) calculate the exact inverse lifetime
  !      at the center of the element, based on exact beam parameters at the center.
  !
  !    The integration of the inverse lifetime along the elememt is then performed
  !      considering a constant inverse lifetime along the element (see above)
  !      and the length of the element
  use name_lenfi
  use touschekfi
  use math_constfi, only : zero, one, two, four, eight, half, pi
  use phys_constfi, only : clight
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   TOUSCHEK SCATTERING, TOUSCHEK Command                              *
  !   These routines implement the formalism from Piwinski               *
  !   (DESY 98-179 & A. Chao/M. Tigner, Handbook of Acc. Physics)        *
  ! Attribute:                                                           *
  !   TABLE     (name)    Name of Twiss table.                           *
  !----------------------------------------------------------------------*
  integer :: i, j, flag, iflag, range(2), table_output, lp, centre
  double precision :: ccost, fact, rr, beta2, gamma2, tolerance, pi2
  double precision :: uloss, km, um, bx, by, ax, ay, dx, dpx, dy, dpy, l, s
  double precision :: sigx2, sigy2, ddx2, ddy2, sigh2
  double precision :: litousch, litouschp, tlitouschek, litouschw, tltouschek
  character(len=name_len) :: name, sequ_name

  integer, external :: get_option, double_from_table_row, restart_sequ
  integer, external :: string_from_table_row, advance_to_pos, get_string
  double precision, external ::  get_value, DGAUSS, ftousch

  ! ************* Get the parameters for the common blocks *************
  ! *************         /machin/ and /beamdb/            *************

  lp = get_string('beam ', 'particle ', sequ_name)

  charge   = get_value('probe ', 'charge ')
  gammas   = get_value('probe ', 'gamma ')
  gamma    = get_value('probe ', 'gamma ')
  en0      = get_value('probe ', 'energy ')
  amass    = get_value('probe ', 'mass ')
  ex       = get_value('probe ', 'ex ')
  ey       = get_value('probe ', 'ey ')
  et       = get_value('probe ', 'et ')
  sigt     = get_value('probe ', 'sigt ')
  sige     = get_value('probe ', 'sige ')
  parnum   = get_value('probe ', 'npart ')
  circ     = get_value('probe ', 'circ ')
  currnt   = get_value('probe ', 'bcurrent ')
  betas    = get_value('probe ', 'beta ')
  beta     = get_value('probe ', 'beta ')
  arad     = get_value('probe ', 'arad ')
  alfa     = get_value('probe ', 'alfa ')
  freq0    = get_value('probe ', 'freq0 ')
  bunch    = get_value('probe ', 'kbunch ')
  deltap   = get_value('probe ','deltap ')

  print *, ''
  print *, 'TOUSCHEK MODULE PARAMETERS'
  print *, ''
  print *, 'particle ', sequ_name(:lp)
  print *, 'charge   ', charge
  print *, 'gammas   ', gammas
  print *, 'gamma    ', gamma
  print *, 'energy   ', en0
  print *, 'mass     ', amass
  print *, 'Ex       ', ex
  print *, 'Ey       ', ey
  print *, 'Et       ', et
  print *, 'sigt     ', sigt
  print *, 'sige     ', sige
  print *, 'parnum   ', parnum
  print *, 'circ     ', circ
  print *, 'currnt   ', currnt
  print *, 'betas    ', betas
  print *, 'beta     ', beta
  print *, 'clight   ', clight
  print *, 'arad     ', arad
  print *, 'alfa     ', alfa
  print *, 'freq0    ', freq0
  print *, 'kbunch   ', bunch
  print *, 'deltap   ', deltap
  print *, ''

  beta2  = beta*beta
  gamma2 = gamma*gamma
  ccost  = arad*arad * clight * parnum / (eight*sqrt(pi) * gamma2*gamma2 * beta2)

  tolerance    = get_value('touschek ', 'tolerance ')
  table_output = get_option('touschek_table ')
  centre       = get_option('centre ')

  if (centre .ne. 0) then
     print *, ' TOUSCHEK will use optical functions at center of elements.'
  else
     print *, ' TOUSCHEK will use optical functions at entrance and end of element'
  endif

  !--- Look for RF cavities and setup RF-bucket boundaries
  call cavtouschek(um,uloss,iflag)
  um1 = um
  !--- setup boundaries for numerical integration
  km = ATAN(sqrt(um1))
  pi2 = pi/two

  if (um1.eq.0) then
     call fort_warn('TOUSCHEK ', '  rf voltage = 0, rest skipped ')
     return
  endif

  if (iflag .eq.1) then
     call fort_warn('TOUSCHEK ', ' uloss = 0 missing chrom in twiss ')
     return
  endif

  !--- Start new Twiss Table reading
  call table_range('twiss ', '#s/#e ', range)
  if (get_option('debug ') .ne. 0) then
     print *, 'Range for Table ', range(1), range(2)
  endif

  j = restart_sequ()

  !--- Initialization (p indices stands for previous element)
  !    A zero value for previous inverse lifetime here does not matter because
  !    first element in sequence is always a marker of zero length whose
  !    contribution will be zero by integration over the length.
  tlitouschek = zero
  litousch  = zero
  litouschp = zero
  litouschw = zero

  !--- Start loop over elements in range
  do i = range(1), range(2)
     j = advance_to_pos('twiss ', i)

     if (string_from_table_row('twiss ', 'name ', i, name) .ne. 0) goto 102
     if (double_from_table_row('twiss ', 's ',    i, s)   .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'betx ', i, bx)  .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'bety ', i, by)  .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'alfx ', i, ax)  .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'alfy ', i, ay)  .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'dx ',   i, dx)  .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'dpx ',  i, dpx) .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'dy ',   i, dy)  .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'dpy ',  i, dpy) .ne. 0) goto 102
     if (double_from_table_row('twiss ', 'l ',    i, l)   .ne. 0) goto 102

     !--- Calculate factors and arguments to numerical integrator
     sigx2 = ex*bx
     sigy2 = ey*by

     ddx2 = (dpx*bx+dx*ax)**2
     ddy2 = (dpy*by+dy*ay)**2

     sigh2 = one / ((one/sige**2)+((dx**2+ddx2)/sigx2)+((dy**2+ddy2)/sigy2))

     fact = sqrt(sigh2)/(sigt*sige*ex*ey)

     !--- Setup the parameters B1 and B2.
     !    They will be used inside the ftousch routine that will be integrated numerically

     fb1 = ( (sigx2-sigh2*ddx2) / ex**2 + (sigy2-sigh2*ddy2) / ey**2 )  &
            / (two*beta2*gamma2)

     ! LD: use first version of eq. 34, more stable for B_2^2 (improvement from G. Roy)
     fb2 = sqrt( &
          ( ( (sigx2-sigh2*ddx2) / ex**2 - (sigy2-sigh2*ddy2) / ey**2 )**2 &
            + (four * sigh2**2 * ddx2 * ddy2) / (ex**2 * ey**2) ) &
            / (four * beta2**2 * gamma2**2) )

     if ( get_option('debug ') .ne. 0 ) &
          write (10,*) s, sigx2, sigy2, ddx2, ddy2, sigh2, fact, fb1, fb2

     !--- Calculate Inverse Touschek Lifetime at current location by numerical integration
     litousch = ccost * fact * DGAUSS(ftousch,km,pi2,tolerance)

     if (centre .ne. 0) then
        !--- Calculate contribution of current element by taking Loss Rate
        !    calculated at centre of element, and normalising by element length
        !    over total length of beamline
        litouschw = litousch * l / circ
        if ( get_option('debug ') .ne. 0 ) &
             write (7,*)  s, litousch, zero, l, litouschw * circ, fb1, fb2
     else
        !--- Calculate contribution of current element by averaging Loss Rate
        !    between exit of current element and entrance of current element,
        !    or end of previous element, and normalising by element length
        !    over total length of beamline
        litouschw = half*(litousch + litouschp) * l / circ
        if ( get_option('debug ') .ne. 0 ) &
             write (8,*)  s, litousch, litouschp, l, litouschw * circ, fb1, fb2
     endif

     !---- Accumulate contributions.
     tlitouschek = tlitouschek + litouschw

     ! *************** Fill "touschek_table"  *********************

     if(table_output.ne.0) then
        call string_to_table_curr('touschek ','name ', name )
        call double_to_table_curr('touschek ','s ', s)               ! position of calculation (end or centre)
        call double_to_table_curr('touschek ','tli ', litousch)      ! Loss rate at current location
        call double_to_table_curr('touschek ','tliw ', litouschw)    ! Weighted Loss Rate contribution of element
        call double_to_table_curr('touschek ','tlitot ', tlitouschek)! Accumulated Loss Rate to present element
        call augment_count('touschek ')
     endif

     !--- Prepare for next element (for calculation at end of element)
     litouschp = litousch

  enddo

  tltouschek = 1d0/tlitouschek

  print *, ' '
  print *, 'Energy radiated per turn  ', uloss,       '[MeV]'
  print *, 'Touschek Inverse Lifetime ', tlitouschek, '[seconds-1]'
  print *, 'Touschek Lifetime         ', tltouschek,  '[seconds]   '
  print *, '                          ', tltouschek/3600.,'[hours]'

  RETURN

102 continue
  call fort_warn('TOUSCHEK ', 'table value not found, rest skipped ')
  return

end subroutine touschek

subroutine cavtouschek (um,uloss,iflag)
  use name_lenfi
  use touschekfi
  use math_constfi, only : zero, one, two, three, half, ten3p, ten6p, ten3m, pi, twopi
  use phys_constfi, only : clight
  use code_constfi
  implicit none

  double precision, intent(OUT) :: um, uloss
  integer, intent(OUT) :: iflag

  integer :: i, lg, code, flagerr
  double precision :: el, rfv, rff, rfl, harmonl, phirf, c0, vrf, pc, omega, orbit5
  double precision :: eta, qover, fq, vrfsum, harmonlm, umt, synch_2
  character(len=name_len) :: sequ_name, el_name

  integer, external :: get_string, restart_sequ, advance_node, double_from_table_row
  double precision, external :: get_value, node_value

  !---- Initialize.
  qover = zero
  vrfsum = zero
  harmonlm = ten6p
  umt = zero
  iflag = 0
  um = zero

  flagerr = double_from_table_row('summ ','synch_2 ',1,synch_2)
  if (synch_2 .eq. zero) then
      uloss = zero ;     iflag = 1
  else
     uloss = two/three * arad * en0**4 * beta**3 * synch_2 * ten3p / amass**3
  endif

  lg = get_string('sequence ', 'name ', sequ_name)
  i = restart_sequ()


  ! treat only cavities
10 if (node_value('mad8_type ') .eq. code_rfcavity) then
     lg = get_string('element ', 'name ', el_name)
     el  = node_value('l ')
     rfv = node_value('volt ')
     rff = node_value('freq ')
     rfl = node_value('lag ')

     if (rff.eq.zero .or. rfv.eq.zero) goto 11 ! skip

     harmonl = ten6p * rff * circ/clight

     pc = get_value('probe ','pc ')

     omega = rff * ten6p * twopi / clight
     vrf   = rfv * ten3m / (pc * (one + deltap))

     flagerr = double_from_table_row('twiss ','t ',1,orbit5)
     phirf = rfl * twopi - omega * orbit5
     c0    = vrf * charge
     if (cos(phirf) .lt. 0) vrf = -vrf
     eta = alfa - one / gammas**2
     if (uloss .ne. zero) then
        qover = qover +  rfv/uloss
        vrfsum = vrfsum +  rfv/harmonl !Charge is not inccluded sine it is not used in twiss or track for the rf-cavity 
        harmonlm = min(harmonl, harmonlm)
     else
        umt = umt + (two * c0) / (harmonl * eta * pi)
     endif

  endif

11 if (advance_node().ne.0)  goto 10

  if (uloss.ne.zero) then
      ! C. Steier, et al., "Measuring and optimizing the momentum aperture in a particle accelerator"
      ! https://journals.aps.org/pre/pdf/10.1103/PhysRevE.65.056506 
     !fq = (two*vrfsum)*(sqrt(one - one/qover**2)) - two*uloss*acos(one/qover)
     fq = two * (sqrt(one - one/qover**2) * vrfsum * harmonlm - uloss * acos(one/qover))
     um = ten3m / (harmonlm * eta * pi) * fq / (pc * (one + deltap))
  else
     um = umt
  endif

  um = abs(um) * beta**2

END subroutine cavtouschek

double precision function ftousch(k)
  use touschekfi
  use math_constfi, only : zero, one, two
  implicit none

  double precision, intent(IN) :: k

  integer :: iflag
  double precision :: z, aftoush
  double precision :: ZR, ZI, BJOR, BJOI, BJIR, BJII, BYOR, BYOI, BYIR, BYII

  z = TAN(k)**2

  ZI = fb2 * z
  ZR = zero

  call CJYDBB(ZR, ZI, BJOR, BJOI, BJIR, BJII, BYOR, BYOI, BYIR, BYII, IFLAG)

  aftoush = two * sqrt(one+z) * &
       ( ((two*z+one)**2 * (z/um1/(one+z) - one)/z) + z - &
       sqrt(z*um1*(one+z)) - (two + one/(two*z))*log(z/um1/(one+z)) )

  if (iflag .eq. 0) then
     ftousch = aftoush * exp(-fb1*z) * BJOR
  else
     ftousch = aftoush * BJOR * (exp(-(fb1-fb2)*z) + exp(-(fb1+fb2)*z))/two
  end if

  return
end function ftousch

! *********************************************************************
!
! $Id: touschek.f90,v 1.4 2010-03-13 02:25:54 frs Exp $
!
! $Log: not supported by cvs2svn $
! Revision 1.3  2009/07/01 17:50:07  frs
! Introducing elements:
! placeholder 38 equivalent with instrument 24
! tkicker 39 equivalent with kicker 15
!
! Revision 1.2  2009/04/06 23:30:07  frs
! Fortran Clean-up: indenting, remove potentially uninitialized variable and
! also remove unused variables
!
! Revision 1.1  2009/03/27 09:37:53  frs
! New files needed for MAD-X Version 4
!
! Revision 1.8  2009/03/21 00:03:43  frs
! Clean-up
!
! Revision 1.7  2008/03/07 10:09:53  frankz
! add ion charge in rf bucket height calculation of Touschek module
!
! Revision 1.6  2005/08/24 13:42:19  frs
! simple clean-up
!
! Revision 1.5  2005/02/05 00:54:43  frs
! Fortran clean-up with mymod, f2c, lf95, f95
!
! Revision 1.4  2005/01/31 09:13:53  frs
! Well tested touschek module
!
! Revision 1.1.1.1  1996/02/15 17:48:17  mclareni
! Kernlib
!
!
!#include "kernnum/pilot.h"
DOUBLE PRECISION FUNCTION DGAUSS(F,A,B,EPS)
  implicit none
  DOUBLE PRECISION :: F, A, B, EPS
  EXTERNAL F

  integer :: i, LGFILE
  DOUBLE PRECISION :: W(12), X(12), AA, BB, C1, C2, U, S8, S16, CONST
  LOGICAL :: MFLAG, RFLAG
  !
  !     ******************************************************************
  !
  !     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.
  !
  !     DGAUSS IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
  !     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
  !     EPS.
  !
  !     ******************************************************************
  !
  !
  DATA W / 0.1012285362903762591525313543D0,                        &
       0.2223810344533744705443559944D0,                                 &
       0.3137066458778872873379622020D0,                                 &
       0.3626837833783619829651504493D0,                                 &
       0.2715245941175409485178057246D-1,                                &
       0.6225352393864789286284383699D-1,                                &
       0.9515851168249278480992510760D-1,                                &
       0.1246289712555338720524762822D0,                                 &
       0.1495959888165767320815017305D0,                                 &
       0.1691565193950025381893120790D0,                                 &
       0.1826034150449235888667636680D0,                                 &
       0.1894506104550684962853967232D0/
  !
  DATA X / 0.9602898564975362316835608686D0,                        &
       0.7966664774136267395915539365D0,                                 &
       0.5255324099163289858177390492D0,                                 &
       0.1834346424956498049394761424D0,                                 &
       0.9894009349916499325961541735D0,                                 &
       0.9445750230732325760779884155D0,                                 &
       0.8656312023878317438804678977D0,                                 &
       0.7554044083550030338951011948D0,                                 &
       0.6178762444026437484466717640D0,                                 &
       0.4580167776572273863424194430D0,                                 &
       0.2816035507792589132304605015D0,                                 &
       0.9501250983763744018531933543D-1/
  !
  !     ******************************************************************
  !
  !  START.
  DGAUSS=0.0D0
  IF(B.EQ.A) RETURN
  CONST=0.005D0/(B-A)
  BB=A
  !
  !  COMPUTATIONAL LOOP.
1 AA=BB
  BB=B
2 C1=0.5D0*(BB+AA)
  C2=0.5D0*(BB-AA)
  S8=0.0D0
  DO I=1,4
     U=C2*X(I)
     S8=S8+W(I)*(F(C1+U)+F(C1-U))
  enddo
  S8=C2*S8
  S16=0.0D0
  DO I=5,12
     U=C2*X(I)
     S16=S16+W(I)*(F(C1+U)+F(C1-U))
  enddo
  S16=C2*S16
  IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5
  BB=C1
  IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2
  DGAUSS=0.0D0
  CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG)
  IF(MFLAG) THEN
     IF(LGFILE.EQ.0) THEN
        WRITE(*,6)
     ELSE
        WRITE(LGFILE,6)
     ENDIF
  ENDIF
  IF(.NOT. RFLAG) CALL ABEND
  RETURN
5 DGAUSS=DGAUSS+S16
  IF(BB.NE.B) GO TO 1
  RETURN
  !
6 FORMAT( 4X, 'FUNCTION DGAUSS ... TOO HIGH ACCURACY REQUIRED')
END FUNCTION DGAUSS

!********************************************************************
!
! $Id: touschek.f90,v 1.4 2010-03-13 02:25:54 frs Exp $
!
! $Log: not supported by cvs2svn $
! Revision 1.3  2009/07/01 17:50:07  frs
! Introducing elements:
! placeholder 38 equivalent with instrument 24
! tkicker 39 equivalent with kicker 15
!
! Revision 1.2  2009/04/06 23:30:07  frs
! Fortran Clean-up: indenting, remove potentially uninitialized variable and
! also remove unused variables
!
! Revision 1.1  2009/03/27 09:37:53  frs
! New files needed for MAD-X Version 4
!
! Revision 1.8  2009/03/21 00:03:43  frs
! Clean-up
!
! Revision 1.7  2008/03/07 10:09:53  frankz
! add ion charge in rf bucket height calculation of Touschek module
!
! Revision 1.6  2005/08/24 13:42:19  frs
! simple clean-up
!
! Revision 1.5  2005/02/05 00:54:43  frs
! Fortran clean-up with mymod, f2c, lf95, f95
!
! Revision 1.4  2005/01/31 09:13:53  frs
! Well tested touschek module
!
! Revision 1.1.1.1  1996/02/15 17:48:35  mclareni
! Kernlib
!
!
SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
  IMPLICIT NONE

  INTEGER KOUNTE,L,I,LIMITM,LIMITR,LOG,LOGF,LGFILE
  PARAMETER(KOUNTE  =  27)
  CHARACTER(6)        ERCODE,   CODE(KOUNTE)
  LOGICAL             MFLAG,    RFLAG
  INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)
  DATA      LOGF      /  0  /
  DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 255, 255 /
  DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 255, 255 /
  DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 255, 255 /
  DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 255, 255 /
  DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 255, 255 /
  DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 255, 255 /
  DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 255, 255 /
  DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 255, 255 /
  DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 255, 255 /
  DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 255, 255 /
  DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 255, 255 /
  DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 255, 255 /
  DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 255, 255 /
  DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 255, 255 /
  DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 255, 255 /
  DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 255, 255 /
  DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 255, 255 /
  DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 255, 255 /
  DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 255, 255 /
  DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 255, 255 /
  DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 255, 255 /
  DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 255,   0 /
  DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 255,   0 /
  DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 255,   0 /
  DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 255,   0 /
  DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 255, 255 /
  DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 255, 255 /
  LOGF  =  LGFILE

  L  =  0
  IF(ERCODE .NE. ' ')  THEN
     DO  L = 1, 6
        IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
     enddo
12   CONTINUE
  ENDIF
  DO     I  =  1, KOUNTE
     IF(L .EQ. 0)  GOTO 13
     IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
13   IF(LIMITM.GE.0) KNTM(I)  =  LIMITM
     IF(LIMITR.GE.0) KNTR(I)  =  LIMITR
14   CONTINUE
  enddo
  RETURN

  ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
  LOG  =  LOGF
  DO     I  =  1, KOUNTE
     IF(ERCODE .EQ. CODE(I))  GOTO 21
  enddo
  WRITE(*,1000)  ERCODE
  CALL ABEND
  RETURN
21 RFLAG  =  KNTR(I) .GE. 1
  IF(RFLAG  .AND.  (KNTR(I) .LT. 255))  KNTR(I)  =  KNTR(I) - 1
  MFLAG  =  KNTM(I) .GE. 1
  IF(MFLAG  .AND.  (KNTM(I) .LT. 255))  KNTM(I)  =  KNTM(I) - 1
  IF(.NOT. RFLAG)  THEN
     IF(LOGF .LT. 1)  THEN
        WRITE(*,1001)  CODE(I)
     ELSE
        WRITE(LOGF,1001)  CODE(I)
     ENDIF
  ENDIF
  IF(MFLAG .AND. RFLAG)  THEN
     IF(LOGF .LT. 1)  THEN
        WRITE(*,1002)  CODE(I)
     ELSE
        WRITE(LOGF,1002)  CODE(I)
     ENDIF
  ENDIF
  RETURN
1000 FORMAT(' KERNLIB LIBRARY ERROR. ' /                               &
       ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',                    &
       ' ERROR MONITOR. RUN ABORTED.')
1001 FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',           &
       'CONDITION ',A6)
1002 FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
END SUBROUTINE KERSET


SUBROUTINE CJYDBB(ZR,ZI,BJOR,BJOI,BJIR,BJII,BYOR,BYOI,BYIR,BYII,IFLAG)
!_________________________________________________
!
!  given the arguments:  (Re=ZR, Im=ZI),  ZI can be .GT. 150
!  returns the modified Bessel functions:
!  J0 (Re=BJOR,Im= BJOI), J1(Re=BJIR, Im=BJII)
!  Y0 (Re=BYOR,Im= BYOI), Y1(Re=BYIR, Im=BYII)
!  if ZI .GT. 150 then iflag = 1 and   BJOR = BJOR/cosh(ZI)
!_________________________________________________

  IMPLICIT NONE

  INTEGER :: iflag, k, j, l, n, m

  DOUBLE PRECISION ZR,ZI,BJOR,BJOI,BJIR,BJII,                       &
       BYOR,BYOI,BYIR,BYII,zmag,angz,cang,                               &
       cangz,sang,sangz,aabs,cterm0,cterm2,                              &
       ang,u,ui,angu,ur,szi,zzi,hcosx,hsinx,                             &
       zimag,sxir,cxir,sxor,cxor,ti,t,zir,zor,                           &
       z2mag,angt,psi,fk,rk,cterm,cterm1,                                &
       bjiit,byort,byoit,byirt,byiit,pi2,gam,tyr,                        &
       tyi,sum1,sum0,cmu0,cmu1,tk2,s,smag,                               &
       spor,spir,spii,sqor,sqoi,sqir,sqii,                               &
       por,poi,spoi,pii,pi4,htanx,byiti,byitr,                           &
       bjiti,bjitr,byoti,byotr,bjoti,bjotr,                              &
       sii,sir,cii,cir,soi,sor,coi,cor,z1,qii,                           &
       bjirt,rl,angs,tr,pir,qor,qoi,qir,sink

  DOUBLE PRECISION, DIMENSION(18) :: CJOR, CJOI, CJIR, CJII
  DOUBLE PRECISION, DIMENSION(18) :: CYOR, CYOI, CYIR, CYII

  DOUBLE PRECISION, DIMENSION(12) :: CPO, CPI, CPOR, CPOI, CPIR, CPII
  DOUBLE PRECISION, DIMENSION(12) :: CQO, CQI, CQOR, CQOI, CQIR, CQII

  DOUBLE PRECISION, PARAMETER :: zero=0.d0

  !---- Initialize.
  CJOR = zero
  CJOI = zero
  CJIR = zero
  CJII = zero

  CYOR = zero
  CYOI = zero
  CYIR = zero
  CYII = zero

  CPO  = zero
  CPI  = zero
  CPOR = zero
  CPOI = zero
  CPIR = zero
  CPII = zero

  CQO  = zero
  CQI  = zero
  CQOR = zero
  CQOI = zero
  CQIR = zero
  CQII = zero

  iflag = 0

  Z2MAG=ZR*ZR+ZI*ZI
  ZMAG=SQRT(Z2MAG)
  ANGZ=ATAN2(ZI,ZR)
  IF(ZMAG - .1D+02.gt.0) goto 21
  T=Z2MAG/.4D+01
  ANGT=ANGZ*.2D+01
  PSI=0.0D+00
  FK=.1D+01
  DO K=1,18
     RK=DBLE(FLOAT(K))
     ANG=ANGT*RK
     SINK=(-.1D+01)**K
     FK=FK*RK
     SANG=SIN(ANG)
     CANG=COS(ANG)
     CTERM=SINK/(FK*FK)
     CJOR(K)=CTERM*CANG
     CJOI(K)=CTERM*SANG
     CTERM1=CTERM/(RK+.1D+01)
     CJIR(K)=CTERM1*CANG
     CJII(K)=CTERM1*SANG
     PSI=PSI+.1D+01/RK
     CTERM0=CTERM*PSI
     CYOR(K)=CTERM0*CANG
     CYOI(K)=CTERM0*SANG
     CTERM2=CTERM1*(PSI+.5D+00/(RK+.1D+01))
     CYIR(K)=CTERM2*CANG
     CYII(K)=CTERM2*SANG
  enddo
  BJOR=(((((((((((((((((CJOR(18)*T+CJOR(17))*T+                     &
       CJOR(16))*T+CJOR(15))*T+CJOR(14))*T+CJOR(13))*T                   &
       +CJOR(12))*T+CJOR(11))*T+CJOR(10))*T+CJOR(9))*T                   &
       +CJOR(8))*T+CJOR(7))*T+CJOR(6))*T+CJOR(5))*T+                     &
       CJOR(4))*T+CJOR(3))*T+CJOR(2))*T+CJOR(1))*T+                      &
       .1D+01
  BJOI=(((((((((((((((((CJOI(18)*T+CJOI(17))*T+                     &
       CJOI(16))*T+CJOI(15))*T+CJOI(14))*T+CJOI(13))*T                   &
       +CJOI(12))*T+CJOI(11))*T+CJOI(10))*T+CJOI(9))*T                   &
       +CJOI(8))*T+CJOI(7))*T+CJOI(6))*T+CJOI(5))*T+                     &
       CJOI(4))*T+CJOI(3))*T+CJOI(2))*T+CJOI(1))*T
  BJIRT=(((((((((((((((((CJIR(18)*T+CJIR(17))*T+                    &
       CJIR(16))*T+CJIR(15))*T+CJIR(14))*T+CJIR(13))*T                   &
       +CJIR(12))*T+CJIR(11))*T+CJIR(10))*T+CJIR(9))*T                   &
       +CJIR(8))*T+CJIR(7))*T+CJIR(6))*T+CJIR(5))*T+                     &
       CJIR(4))*T+CJIR(3))*T+CJIR(2))*T+CJIR(1))*T+                      &
       .1D+01
  BJIIT=(((((((((((((((((CJII(18)*T+CJII(17))*T+CJII(16))           &
       *T+CJII(15))*T+CJII(14))*T+CJII(13))*T+CJII(12))*T                &
       +CJII(11))*T+CJII(10))*T+CJII(9))*T+CJII(8))*T+                   &
       CJII(7))*T+CJII(6))*T+CJII(5))*T+CJII(4))*T+                      &
       CJII(3))*T+CJII(2))*T+CJII(1))*T
  BYORT=(((((((((((((((((CYOR(18)*T+CYOR(17))*T+CYOR(16))           &
       *T+CYOR(15))*T+CYOR(14))*T+CYOR(13))*T+CYOR(12))*T+               &
       CYOR(11))*T+CYOR(10))*T+CYOR(9))*T+CYOR(8))*T+                    &
       CYOR(7))*T+CYOR(6))*T+CYOR(5))*T+CYOR(4))*T+CYOR(3))              &
       *T+CYOR(2))*T+CYOR(1))*T
  BYOIT=(((((((((((((((((CYOI(18)*T+CYOI(17))*T+CYOI(16))*T         &
       +CYOI(15))*T+CYOI(14))*T+CYOI(13))*T+CYOI(12))*T+                 &
       CYOI(11))*T+CYOI(10))*T+CYOI(9))*T+CYOI(8))*T+CYOI(7)             &
       )*T+CYOI(6))*T+CYOI(5))*T+CYOI(4))*T+CYOI(3))*                    &
       T+CYOI(2))*T+CYOI(1))*T
  BYIRT=(((((((((((((((((CYIR(18)*T+CYIR(17))*T+CYIR(16))*T+        &
       CYIR(15))*T+CYIR(14))*T+CYIR(13))*T+CYIR(12))*T+                  &
       CYIR(11))*T+CYIR(10))*T+CYIR(9))*T+CYIR(8))*T+CYIR(7)             &
       )*T+CYIR(6))*T+CYIR(5))*T+CYIR(4))*T+CYIR(3))*T+                  &
       CYIR(2))*T+CYIR(1))*T+.5D+00
  BYIIT=(((((((((((((((((CYII(18)*T+CYII(17))*T+CYII(16))*T         &
       +CYII(15))*T+CYII(14))*T+CYII(13))*T+CYII(12))*T+                 &
       CYII(11))*T+CYII(10))*T+CYII(9))*T+CYII(8))*T+CYII(7)             &
       )*T+CYII(6))*T+CYII(5))*T+CYII(4))*T+CYII(3))*T+                  &
       CYII(2))*T+CYII(1))*T
  CANGZ=COS(ANGZ)
  SANGZ=SIN(ANGZ)
  BJIR=(ZR*BJIRT-ZI*BJIIT)/.2D+01
  BJII=(ZI*BJIRT+ZR*BJIIT)/.2D+01
  PI2=.2D+01*ATAN(.1D+01)
  GAM=.577215664901533D+00
  TYR=LOG(ZMAG/.2D+01)+GAM
  TYI=ANGZ
  BYOR=(TYR*BJOR-TYI*BJOI-BYORT)/PI2
  BYOI=(TYI*BJOR+TYR*BJOI-BYOIT)/PI2
  BYIR=(TYR*BJIR-TYI*BJII-CANGZ/ZMAG-(BYIRT*ZR-                     &
       BYIIT*ZI)/.2D+01)/PI2
  BYII=(TYI*BJIR+TYR*BJII+SANGZ/ZMAG-(BYIIT*ZR+BYIRT*ZI)/.2D+01)/   &
       &PI2
  RETURN
21 FK=.1D+01
  if (ZR .ne.0) then
     AABS=ABS(ZI/ZR)
     ANG=ATAN(AABS)
  else
     ang = ASIN(1d0)
  end if
  IF(ZI.lt.0) goto 50
  ANGZ=ANG
  GOTO 52
50 ANGZ=-ANG
52 CONTINUE
  SUM0=.1D+01
  SUM1=.1D+01
  CMU0=0.0
  CMU1=.4D+01
  L=1
  DO K=1,24
     L=-L
     RK=DBLE(FLOAT(K))
     FK=FK
     TK2=4.*RK-4.+1./RK
     SUM0=SUM0*(CMU0-TK2)
     SUM1=SUM1*(CMU1/RK-TK2)
     IF(L.gt.0) goto 23
     N=(K+1)/2
     CQO(N)=SUM0/FK
     CQI(N)=SUM1/FK
     GOTO 30
23   M=K/2
     CPO(M)=SUM0/FK
     CPI(M)=SUM1/FK
30   CONTINUE
  enddo
  RL=-.1D+01
  DO J=1,12
     RL=-RL
     CQO(J)=CQO(J)*RL
     CQI(J)=CQI(J)*RL
     CPO(J)=CPO(J)*RL
     CPI(J)=CPI(J)*RL
  enddo
  ANGS=-ANGZ*.2D+01
  SMAG=.1D+01/(.64D+02*Z2MAG)
  S=SMAG
  DO L=1,12
     RL=DBLE(FLOAT(L))
     ANG=ANGS*RL
     CANG=COS(ANG)
     SANG=SIN(ANG)
     CPOR(L)=CPO(L)*CANG
     CPOI(L)=CPO(L)*SANG
     CPIR(L)=CPI(L)*CANG
     CPII(L)=CPI(L)*SANG
     CQOR(L)=CQO(L)*CANG
     CQOI(L)=CQO(L)*SANG
     CQIR(L)=CQI(L)*CANG
     CQII(L)=CQI(L)*SANG
  enddo
  SPOR=(((((((((((CPOR(12)*S+CPOR(11))*S+CPOR(10))*S+               &
       CPOR(9))*S+CPOR(8))*S+CPOR(7))*S+CPOR(6))*S+CPOR(5))              &
       *S+CPOR(4))*S+CPOR(3))*S+CPOR(2))*S+CPOR(1))*S
  SPOI=(((((((((((CPOI(12)*S+CPOI(11))*S+CPOI(10))*S+               &
       CPOI(9))*S+CPOI(8))*S+CPOI(7))*S+CPOI(6))*S+CPOI(5))              &
       *S+CPOI(4))*S+CPOI(3))*S+CPOI(2))*S+CPOI(1))*S
  SPIR=(((((((((((CPIR(12)*S+CPIR(11))*S+CPIR(10))*S+               &
       CPIR(9))*S+CPIR(8))*S+CPIR(7))*S+CPIR(6))*S+CPIR(5))              &
       *S+CPIR(4))*S+CPIR(3))*S+CPIR(2))*S+CPIR(1))*S
  SPII=(((((((((((CPII(12)*S+CPII(11))*S+CPII(10))*S+               &
       CPII(9))*S+CPII(8))*S+CPII(7))*S+CPII(6))*S+CPII(5))              &
       *S+CPII(4))*S+CPII(3))*S+CPII(2))*S+CPII(1))*S
  SQOR=(((((((((((CQOR(12)*S+CQOR(11))*S+CQOR(10))*S+               &
       CQOR(9))*S+CQOR(8))*S+CQOR(7))*S+CQOR(6))*S+CQOR(5))              &
       *S+CQOR(4))*S+CQOR(3))*S+CQOR(2))*S+CQOR(1))*S
  SQOI=(((((((((((CQOI(12)*S+CQOI(11))*S+CQOI(10))*S+               &
       CQOI(9))*S+CQOI(8))*S+CQOI(7))*S+CQOI(6))*S+CQOI(5))              &
       *S+CQOI(4))*S+CQOI(3))*S+CQOI(2))*S+CQOI(1))*S
  SQIR=(((((((((((CQIR(12)*S+CQIR(11))*S+CQIR(10))*S+               &
       CQIR(9))*S+CQIR(8))*S+CQIR(7))*S+CQIR(6))*S+CQIR(5))              &
       *S+CQIR(4))*S+CQIR(3))*S+CQIR(2))*S+CQIR(1))*S
  SQII=(((((((((((CQII(12)*S+CQII(11))*S+CQII(10))*S+               &
       CQII(9))*S+CQII(8))*S+CQII(7))*S+CQII(6))*S+CQII(5))              &
       *S+CQII(4))*S+CQII(3))*S+CQII(2))*S+CQII(1))*S
  T=.8D+01*ZMAG
  ANGT=ANGZ
  TR=T*COS(ANGT)
  TI=T*SIN(ANGT)
  POR=.1D+01-SPOR
  POI=-SPOI
  PIR=.1D+01-SPIR
  PII=-SPII
  QOR=TR*SQOR-TI*SQOI
  QOI=TI*SQOR+TR*SQOI
  QIR=TR*SQIR-TI*SQII
  QII=TI*SQIR+TR*SQII
  PI4=ATAN(.1D+01)
  PI2=.2D+01*PI4
  Z1=ABS(ZR)
  ZOR=Z1-PI4
  ZIR=Z1-(PI2+PI4)
  ZZI=ZI
  CXOR=COS(ZOR)
  SXOR=SIN(ZOR)
  CXIR=COS(ZIR)
  SXIR=SIN(ZIR)
  ZIMAG=ABS(ZZI)
  !***
  IF(ZIMAG .GT. .15D+03) then
     !         WRITE(*,100)
     !         WRITE(*,*) iflag,zimag
     iflag = 1
  end if
100 FORMAT(1X,'ZIMAG .GT. .15D+03 (SUBROUTINE CJYDBB)')
  !***
  IF(ZIMAG - .15D+03.gt.0) goto 61
  HSINX=(EXP(ZZI)-EXP(-ZZI))/.2D+01
  HCOSX=(EXP(ZZI)+EXP(-ZZI))/.2D+01
  COR=CXOR*HCOSX
  COI=-SXOR*HSINX
  SOR=SXOR*HCOSX
  SOI=CXOR*HSINX
  CIR=CXIR*HCOSX
  CII=-SXIR*HSINX
  SIR=SXIR*HCOSX
  SII=CXIR*HSINX
  GOTO 62
61 SZI=SNGL(ZZI)
  HTANX=DBLE(TANH(SZI))
  COR=CXOR
  COI=-SXOR*HTANX
  SOR=SXOR
  SOI=CXOR*HTANX
  CIR=CXIR
  CII=-SXIR*HTANX
  SIR=SXIR
  SII=CXIR*HTANX
62 CONTINUE
  BJOTR=(POR*COR-POI*COI)-(QOR*SOR-QOI*SOI)
  BJOTI=(POI*COR+POR*COI)-(QOI*SOR+QOR*SOI)
  BYOTR=(POR*SOR-POI*SOI)+(QOR*COR-QOI*COI)
  BYOTI=(POI*SOR+POR*SOI)+(QOI*COR+QOR*COI)
  BJITR=(PIR*CIR-PII*CII)-(QIR*SIR-QII*SII)
  BJITI=(PII*CIR+PIR*CII)-(QII*SIR+QIR*SII)
  BYITR=(PIR*SIR-PII*SII)+(QIR*CIR-QII*CII)
  BYITI=(PII*SIR+PIR*SII)+(QII*CIR+QIR*CII)
  U=SQRT(.1D+01/(PI2*ZMAG))
  ANGU=-ANGZ/.2D+01
  UR=U*COS(ANGU)
  UI=U*SIN(ANGU)
  BJOR=(UR*BJOTR-UI*BJOTI)
  BJOI=(UI*BJOTR+UR*BJOTI)
  BYOR=(UR*BYOTR-UI*BYOTI)
  BYOI=(UI*BYOTR+UR*BYOTI)
  BJIR=(UR*BJITR-UI*BJITI)
  BJII=(UI*BJITR+UR*BJITI)
  BYIR=(UR*BYITR-UI*BYITI)
  BYII=(UI*BYITR+UR*BYITI)
  IF(ZR.lt.0) goto 40
  RETURN
40 BYOR=BYOR+.2D+01*BJOI
  BYOI=-BYOI+.2D+01*BJOR
  BYIR=-BYIR-.2D+01*BJII
  BYII=BYII-.2D+01*BJIR
  BJOI=-BJOI
  BJIR=-BJIR

  RETURN
END SUBROUTINE CJYDBB

subroutine abend
    implicit none
    write(*,*) 'Abnormal end ...'
    stop 888
end subroutine abend
