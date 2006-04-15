!===================================================================================
subroutine photon(iele,rad,dlength,energy,dmass,ieave,iquasto,d1,d2)
  !SUBROUTINE photon (i_elem_type, rad_curv_m, length_curr_elem, &
  !                        Energy_beam_MeV, mass_rest_MeV, d1, d2)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     subroutine to determine the energy losses at entrance
  !     and exit of magnet (dipole or quadrupole)
  !
  !
  !     input:
  !            iele - Magnet Type, iele=2 means a Quadrupole and
  !                   its Radiation can be switched off with iquasto flag
  !                   set in the Initialisation
  !            rad - radius of curvature
  !                  (rho for bend
  !                   1/(abs(k1)*sqrt(x**2+y**2)) for quadrupole
  !
  !            dlength - length of element
  !            energy - beam energy in MeV
  !            dmass - particle mass in MeV
  !     output:
  !            d1 - momentum loss at entrance (normalized to beam energy)
  !            d2 - momentum loss at exit (normalized to beam energy)
  !
  !                                                           FZ/28.02.2001
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  USE precision_constants, ONLY:CLIGHT,dhbar,pmae,dp,sp,qelect,c_1d6,c_1d_10,c_1d_3,c_0_3079,c_137,one,two,three,five
  implicit none
  include 'photoni.inc'

  INTEGER, INTENT(IN) :: iele, ieave, iquasto
  REAL(dp), INTENT(IN) :: dlength,energy,dmass
  REAL(dp), INTENT(OUT) :: d1,d2
  REAL(dp), INTENT(INOUT) :: rad

  REAL(dp) :: gamma,theta, dnpho, dnphohalf,ec,eave,ep1,ep2, ra2, &
       ran,xl,en,gamfac

  REAL(dp), parameter :: el=qelect, vl=CLIGHT
  real(sp) amu
  Double Precision :: ran2 ! function
  !      integer i,j,nmax,n1,n2,npmax,iseed,idumy,ieave,iquasto,iele
  !      double precision pi,vl,dhbar,el,gamma,energy,xilog,sollog,b,sol,  &
  !     &xil,sll,ak,bk,theta,dlength,rad,dnpho,dnphohalf,ec,eave,ep1,ep2,  &
  !     &dmass,d1,d2,ra2,ran,ran2,xl,en,gamfac
  !      parameter (nmax=100)

  !      common / lookup / xil(nmax), sll(nmax)
  !      common / lookup2 / ak(nmax), bk(nmax)
  !      common / rann / iseed, idumy
  !      common / trick / pi, vl, dhbar, el, gamma
  !      common / eave1 / ieave,iquasto

  !      real amu, amax
  !      integer ierr,nphoton
  INTEGER :: ierr, nphoton, n1, n2, npmax, i,j, idumy
  external rnpset, rnpssn
  ierr=0
  rad=abs(rad)
  !      print*,"RRRRRRRRRRRRRRRRRad in photon",iele,rad,dlength,energy,dmass

  ! for electrons
  gamma=energy/pmae*c_1d_3

  if (gamma.le.1) then
     write(*,*) ' error in subroutine photon - no initialization '
     return
  endif

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     Stop Quadrupole Radiation
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (iquasto.eq.0.and.iele.eq.2) return
  ! it is check already before CALL this subr.

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     compute mean number of photons emitted
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

989 continue

  !     trajectory bending angle
  theta = dlength/rad

  dnpho = five/two/sqrt(three)*gamma/c_137*theta
  !
  dnphohalf = dnpho/two
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     compute critical energy and average (half) energy loss
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     critical photon energy as fraction of beam energy

  ec = three/two*vl*gamma**3/rad*dhbar/(energy*el*c_1d6)
  !
  eave = dnphohalf*c_0_3079*ec
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     compute number of photons n1 and n2 actually emitted
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  amu = dnphohalf
  !      write(*,*) ' call rnpssn 1 ',amu,ierr
  call rnpssn(amu,nphoton,ierr)
  !      write(*,*) ' call rnpssn 2 ',amu,ierr,nphoton

  n1 = nphoton
  call rnpssn(amu,nphoton,ierr)
  !      write(*,*) ' after calling rnpssn ',amu,ierr,nphoton
  n2 = nphoton

  npmax = n1 + n2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     toss randum number to determine photon energy
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ep1 = 0.
  ep2 = 0.

  do i = 1, npmax

     ra2 = ran2(idumy)
     if (ra2.lt.c_1d_10) goto 989

     ran = log(ra2)
     do j = 1, nmax
        if (ran.lt.sll(j)) then
           if(j.eq.1) then
              xl = (ran-bk(1))/ak(1)
              goto 990
           endif
           xl = (ran-bk(j))/ak(j)
           goto 990
        endif
     enddo
     xl = (ran-bk(nmax))/ak(nmax)
990  continue

     !
     !     note: xl is logarithm of photon energy
     !     divided by critical energy
     !

     if (xl.lt.(-46.).or.(xl.gt.(2.3))) goto 989

     !     photon energy as fraction of beam energy
     en = -exp(xl)*ec

     if (i.le.n1) then
        ep1 = ep1 + en
     else
        ep2 = ep2 + en
     endif

  enddo

  !     compute total relative photon energy (average loss is
  !     added for compensation)
  !
  gamfac = gamma*gamma/(gamma*gamma-one)
  if( ieave.eq.1 ) then
     d1 = (ep1 + eave) * gamfac
     d2 = (ep2 + eave) * gamfac
  else
     d1 = ep1 * gamfac
     d2 = ep2 * gamfac
  endif

  return
end subroutine photon
!===================================================================================


!===================================================================================
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function   ran2(idum)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !   The Function Returns a Uniform Random Number between 0.0 and 1.0.
  !
  !   The Function uses a 'Subtractive Method'!!!
  !
  !   Set 'idum' to a negative value to initialize or reinitialize the Sequence.
  !
  !   (Numerical Recipies ran3(idumy), pg. 273.)
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  USE precision_constants, ONLY:CLIGHT,dhbar,pmae,dp,sp,qelect,c_1d6,c_1d_10,c_1d_3,c_0_3079,c_137,one,two,three,five
  implicit none
  integer    idum
  integer    mbig, mseed, mz
  double precision     fac
  parameter(mbig=1000000000,mseed=161803398,mz=0,fac=one/mbig)
  integer    i, iff, ii, inext, inextp, k
  integer    mj, mk, ma(55)
  save       iff, inext, inextp, ma
  data       iff /0/

  !     Initialization:
  if(idum.lt.0.or.iff.eq.0) then
     iff = 1
     mj = mseed-iabs(idum)
     mj = mod(mj,mbig)
     mk = 1
     do i=1,54
        ii = mod(2*i,55)
        ma(ii) = mk
        mk = mj - mk
        if(mk.lt.mz) mk = mk + mbig
        mj = ma(ii)
11      continue
     enddo

     do k=1,4
        do i=1,55
           ma(i) = ma(i) - ma(1+mod(i+30,55))
           if(ma(i).lt.mz) ma(i) = ma(i) + mbig
12         continue
        enddo
13      continue
     enddo
     inext = 0
     inextp = 31
     idum = 1
  endif

  inext = inext + 1
  if(inext.eq.56) inext = 1
  inextp = inextp + 1
  if(inextp.eq.56) inextp = 1
  mj = ma(inext) - ma(inextp)
  if(mj.lt.mz) mj = mj + mbig
  ma(inext) = mj
  ran2 = mj * fac

  return
END  function   ran2
!========================================================================================
