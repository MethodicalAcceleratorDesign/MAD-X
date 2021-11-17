subroutine trdynrun (eigen,coords,turns,npart,distvect,zn,onelog,turnnumber,dq)
  use io_units, only : dynapout
  use deltrafi
  use wmaxmin0fi
  use dyntabfi
  use tunesfi
  use math_constfi, only : zero, one, two, pi, twopi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !           PLEASE DOCUMENT THE CODE!!!
  !----------------------------------------------------------------------*
  integer :: turns, npart
  double precision :: eigen(6,6), coords(6,0:turns,*), distvect(turns)
  double precision :: zn(turns,6), onelog(turns), turnnumber(turns), dq(2*turns)

  integer i, j, k, ix, iy, initt, ktrturns, nturnhalf
  double precision :: znt(6), track(6)
  double precision :: fitlyap, templyap
  double precision :: tunx1, tunx2, tuny1, tuny2, tuneabt, tuneabt2
  double precision :: dphi(3),dphi1,dphi2

  double precision, external :: get_value

  write(*,*) ' entered dynap '

  fastune = get_value('dynap ','fastune ') .ne. 0
  deltax  = get_value('dynap ','lyapunov ')

  dktrturns = turns

  !---- Initialize max and min betatron invariants.
  wxmax = zero ;  wymax = zero ;  wxymax = zero
  ! min should be initialized to huge value, not zero
  wxmin = 1.d20 ; wymin = 1.d20 ; wxymin = 1.d20

  !---- compute minimum and maximum amplitudes and normalized coordinates
  do k = 1, npart, 2 ! select only main particles, not their evaluation partner
     x0 = coords(1,0,k)
     y0 = coords(3,0,k)
     do i = 1, turns
        TRACK = COORDS(:,i,k)
        call wmaxmin(track,eigen,znt)
        ZN(i,:) = ZNT
     end do

     !---- compute 'smear'
     smear = two * (wxymax - wxymin) / (wxymax + wxymin)

     !---- Fast tune calculation by interpolated FFT.
     if (fastune) then

        ktrturns = turns
        ix = 1 ; iy = 3 
        initt = 0

        if (ktrturns .le. 64) then
           tunx = tuneabt(zn, ix, initt, ktrturns, turns, dq)
           tuny = tuneabt(zn, iy, initt, ktrturns, turns, dq)
        else
           tunx = tuneabt2(zn, ix, initt, ktrturns, turns, dq)
           tuny = tuneabt2(zn, iy, initt, ktrturns, turns, dq)
        endif

        !---- Fast tune variation over half the number of turns.
        nturnhalf = ktrturns / 2
        if (nturnhalf .le. 64) then
           initt = 0
           tunx1 = tuneabt(zn, ix, initt, nturnhalf, turns, dq)
           tuny1 = tuneabt(zn, iy, initt, nturnhalf, turns, dq)
           initt = nturnhalf
           tunx2 = tuneabt(zn, ix, initt, nturnhalf, turns, dq)
           tuny2 = tuneabt(zn, iy, initt, nturnhalf, turns, dq)
        else
           initt = 0
           tunx1 = tuneabt2(zn, ix, initt, nturnhalf, turns, dq)
           tuny1 = tuneabt2(zn, iy, initt, nturnhalf, turns, dq)
           initt = nturnhalf
           tunx2 = tuneabt2(zn, ix, initt, nturnhalf, turns, dq)
           tuny2 = tuneabt2(zn, iy, initt, nturnhalf, turns, dq)
        endif

        if (abs(tunx1-tunx) .gt. 0.4) tunx1 = 1-tunx1
        if (abs(tunx2-tunx) .gt. 0.4) tunx2 = 1-tunx2
        if (abs(tuny1-tuny) .gt. 0.4) tuny1 = 1-tuny1
        if (abs(tuny2-tuny) .gt. 0.4) tuny2 = 1-tuny2

        dtune = sqrt((tunx2 - tunx1)**2 + (tuny2 - tuny1)**2)

        call dynaptunefill
     endif

     !---- Lyapunov exponent calculation

     open(unit=dynapout,file="lyapunov.data")

     do i = 1, turns
        TRACK = COORDS(:,i,k+1) 
        call wmaxmin(track,eigen,znt) 
        do j = 1, 3
           dphi1 = zero ; dphi2 = zero
           if (znt(2*j).ne.zero .or. znt(2*j-1).ne.zero) & 
                dphi1 = atan2(znt(2*j),znt(2*j-1))
           if (zn(i,2*j).ne.zero .or. zn(i,2*j-1).ne.zero) & 
                dphi2 = atan2(zn(i,2*j),zn(i,2*j-1))           

           dphi(j) =  dphi1 - dphi2
           if (dphi(j) .gt.  pi) dphi(j) = dphi(j) - twopi
           if (dphi(j) .lt. -pi) dphi(j) = dphi(j) + twopi
        end do
        distvect(i) = sqrt(dot_product(dphi, dphi))
        write(dynapout,*) i, distvect(i)
     end do

     lyapunov = fitlyap(distvect,onelog,turnnumber,turns,deltax)
     call dynapfill()

  enddo

  close(unit=dynapout, status='keep')
  write(*,*) ' end dynap '

end subroutine trdynrun

subroutine dynapfill()
  use dyntabfi
  use wmaxmin0fi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   writes the dynap output in the table "dynap"                       *
  ! Output:                                                              *
  !   dynapfrac (real)    : fractional dynamic aperture                  *
  !   dktrturns (real)    : number of turns after tracking               *
  !   xend      (real)    : final x position w.r.t. closed orbit         *
  !   pxend     (real)    : final x momentum w.r.t. closed orbit         *
  !   yend      (real)    : final y position w.r.t. closed orbit         *
  !   pyend     (real)    : final y momentum w.r.t. closed orbit         *
  !   tend      (real)    : final longit. position w.r.t. closed orbit   *
  !   ptend     (real)    : final longit. momentum  w.r.t. closed orbit  *
  !   wxmin     (real)    : minimum x betatron invariant during tracking *
  !   wxmax     (real)    : maximum x betatron invariant during tracking *
  !   wymin     (real)    : minimum y betatron invariant during tracking *
  !   wymax     (real)    : maximum y betatron invariant during tracking *
  !   wxymin    (real)    : minimum of (wx + wy) during tracking         *
  !   wxymax    (real)    : maximum of (wx + wy) during tracking         *
  !   smear     (real)    : 2.0 * (wxymax - wxymin) / (wxymax + wxymin)  *
  !   lyapunov  (real)    : interpolated Lyapunov exponent               * ghislain
  !----------------------------------------------------------------------*

  call double_to_table_curr('dynap ', 'dktrturns ', dktrturns)
  call double_to_table_curr('dynap ', 'xend ',      xend )
  call double_to_table_curr('dynap ', 'pxend ',     pxend )
  call double_to_table_curr('dynap ', 'yend ',      yend )
  call double_to_table_curr('dynap ', 'pyend ',     pyend )
  !call double_to_table_curr('dynap ', 'tend ',     ptend ) ! ??? 
  call double_to_table_curr('dynap ', 'tend ',      tend )
  call double_to_table_curr('dynap ', 'ptend ',     ptend )
  call double_to_table_curr('dynap ', 'wxmin ',     wxmin )
  call double_to_table_curr('dynap ', 'wxmax ',     wxmax )
  call double_to_table_curr('dynap ', 'wymin ',     wymin )
  call double_to_table_curr('dynap ', 'wymax ',     wymax )
  call double_to_table_curr('dynap ', 'wxymin ',    wxymin )
  call double_to_table_curr('dynap ', 'wxymax ',    wxymax )
  call double_to_table_curr('dynap ', 'smear ',     smear )
  call double_to_table_curr('dynap ', 'lyapunov ',  lyapunov )

  call augment_count('dynap ')
end subroutine dynapfill

subroutine dynaptunefill
  use tunesfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   writes the dynap tunes in the table "dynaptune"                    *
  ! Output:                                                              *
  !   tunx      (real)    : x tune from FFT                              *
  !   tuny      (real)    : y tune from FFT                              *
  !   dtune     (real)    : tune error                                   *
  !----------------------------------------------------------------------*
  double precision :: tx_tmp, ty_tmp

  tx_tmp = tunx ; ty_tmp = tuny
  if (tx_tmp .gt. 0.5d0)  tx_tmp = 1.d0 - tx_tmp
  if (ty_tmp .gt. 0.5d0)  ty_tmp = 1.d0 - ty_tmp
  call double_to_table_curr('dynaptune ', 'x '    , x0)
  call double_to_table_curr('dynaptune ', 'y '    , y0)
  call double_to_table_curr('dynaptune ', 'tunx ' , tx_tmp)
  call double_to_table_curr('dynaptune ', 'tuny ' , ty_tmp)
  call double_to_table_curr('dynaptune ', 'dtune ', dtune )

  call augment_count('dynaptune ')
end subroutine dynaptunefill

subroutine fft(data, nn, isign)
  use math_constfi, only : zero, one, two, half, twopi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes the FFT                                                   *
  ! Author: numerical receipes,  pg. 395                                 *
  !   DATA:     is a real array with the signal on input                 *
  !             with the Fourier transform on output.                    *
  !   N:        is the number of data: must be a power of 2              *
  !   ISIGN=1:  direct Fourier transform                                 *
  !   ISIGN=-1: inverse Fourier transform                                *
  !----------------------------------------------------------------------*
  !     frankz: only removed a set of parameter statements - keep rest
  !             unchanged
  !
  !---- Double precision version.
  integer :: nn,isign
  double precision :: data(*)

  integer :: i, j, m, n, istep, mmax
  double precision :: tempi, tempr, theta, wi, wpi, wpr, wr, wtemp

  !---- Rearrange the data points.
  n = 2 * nn
  j = 1
  do i = 1, n, 2
     if(j .gt. i) then
        tempr = data(j)
        tempi = data(j+1)
        data(j) = data(i)
        data(j+1) = data(i+1)
        data(i) = tempr
        data(i+1) = tempi
     endif
     m = n / 2
1    if (m .ge. 2  .and.  j .gt. m) then
        j = j - m
        m = m / 2
        go to 1
     endif
     j = j + m
  enddo
  mmax = 2
2 if (n .gt. mmax) then
     istep = 2 * mmax
     theta = twopi / (isign * mmax)
     wpr = - two * sin(half * theta)**2
     wpi = sin(theta)
     wr = one
     wi = zero
     do m = 1, mmax, 2
        do i = m, n, istep
           j = i + mmax
           tempr = wr * data(j)   - wi * data(j+1)
           tempi = wr * data(j+1) + wi * data(j)
           data(j)   = data(i)   - tempr
           data(j+1) = data(i+1) - tempi
           data(i)   = data(i)   + tempr
           data(i+1) = data(i+1) + tempi
        enddo
        wtemp = wr
        wr = wr * wpr - wi    * wpi + wr
        wi = wi * wpr + wtemp * wpi + wi
     enddo
     mmax = istep
     go to 2
  endif

end subroutine fft

double precision function tuneabt(zn, ixy, initt, maxn, turns, dq)
  use math_constfi, only : zero, one, two, pi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes the tune using formula (18) of CERN SL/95-84 (AP).        *
  !   No filter.                                                         *
  ! (best suited for MAXN <= 64 TURNS)                                   *
  ! zn(*,4) contains the 4 transverse coordinates for different turns    *
  ! maxn is the number of turns                                          *
  ! Authors:                                                             *
  !   R. Bartolini - CERN and Bologna University,                        *
  !   E. TODESCO   - INFN and CERN.                                      *
  !----------------------------------------------------------------------*
  integer :: ixy, initt, maxn, turns
  double precision :: zn(turns,6), dq(2*turns)

  integer :: mft, nft, nftmax, npoint, mf
  double precision :: ftmax, temp, cf1, cf2, cf3, arg, assk

  !---- Use first NPOINT points.
  mft = int(log(float(maxn)) / log(two))
  npoint = 2**mft

  !---- Copy data to dq
  do mf = 1, npoint
     dq(2*mf-1) = zn(mf+initt,ixy)
     dq(2*mf)   = zn(mf+initt,ixy+1)
  enddo
  call fft(dq, npoint, -1)

  !---- Search for maximum of Fourier spectrum.
  ftmax = zero
  nftmax = 0
  do nft = 1, npoint
     temp = sqrt(dq(2*nft-1)**2 + dq(2*nft)**2)
     if (temp .gt. ftmax) then
        ftmax  = temp
        nftmax = nft
     endif
  enddo

  !---- Improve estimate by interpolation.
  cf1 = sqrt(dq(2*nftmax-3)**2 + dq(2*nftmax-2)**2)
  cf2 = sqrt(dq(2*nftmax-1)**2 + dq(2*nftmax)**2)
  cf3 = sqrt(dq(2*nftmax+1)**2 + dq(2*nftmax+2)**2)
  if (cf3 .gt. cf1) then
     arg  = sin(pi / npoint) / (cf2 / cf3 + cos(pi / npoint))
     assk = float(nftmax) + npoint / pi * atan(arg)
  else
     arg  = sin(pi / npoint) / (cf1 / cf2 + cos(pi / npoint))
     assk = float(nftmax-1) + npoint / pi * atan(arg)
  endif
  tuneabt = one - (assk - one) / float(npoint)

end function tuneabt

double precision function tuneabt2(zn, ixy, initt, maxn, turns,dq)
  use math_constfi, only : zero, one, two, pi, twopi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes the tune using the interpolated FFT with Hanning filter.  *
  !   See CERN SL/95-84 formula (25).                                    *
  !   (best suited for MAXN > 64 turns)                                  *
  ! X, XP are the coordinates of the orbit,                              *
  ! MAXN  is the length of the orbit.                                    *
  ! Authors:                                                             *
  !   R. Bartolini - CERN and Bologna University,                        *
  !   E. TODESCO   - INFN and CERN.                                      *
  !----------------------------------------------------------------------*
  integer :: ixy, initt, maxn, turns
  double precision :: zn(turns,6), dq(2*turns)

  integer :: mft, nft, mf, nftmax, npoint, nn
  double precision :: cf1, cf2, cf3, scra1, scra2, scra3, scra4
  double precision :: ftmax, temp, step, assk, co, si, p1, p2

  !---- Use first NPOINT points.
  mft = int(log(float(maxn)) / log(two))
  npoint = 2**mft

  !---- Copy data to local storage using Hanning filter.
  step = pi / npoint
  do mf = 1, npoint
     temp = sin(mf * step)**2
     dq(2*mf-1) = temp * zn(mf+initt,ixy)
     dq(2*mf)   = temp * zn(mf+initt,ixy+1)
  enddo
  call fft(dq, npoint, -1)

  !---- Search for maximum of Fourier spectrum.
  ftmax = zero
  nftmax = 0
  do nft = 1, npoint
     temp = sqrt(dq(2*nft-1)**2 + dq(2*nft)**2)
     if (temp .gt. ftmax) then
        ftmax  = temp
        nftmax = nft
     endif
  enddo
  cf1 = sqrt(dq(2*nftmax-3)**2 + dq(2*nftmax-2)**2)
  cf2 = sqrt(dq(2*nftmax-1)**2 + dq(2*nftmax)**2)
  cf3 = sqrt(dq(2*nftmax+1)**2 + dq(2*nftmax+2)**2)
  if (cf3 .gt. cf1) then
     p1 = cf2
     p2 = cf3
     nn = nftmax
  else
     p1 = cf1
     p2 = cf2
     nn = nftmax-1
  endif

  !---- Interpolation.
  co = cos(twopi / npoint)
  si = sin(twopi / npoint)
  scra1 = co**2 * (p1 + p2)**2 - 2*p1*p2*(2*co**2 - co - one)
  scra2 = (p1 + p2*co) * (p1 - p2)
  scra3 = p1**2 + p2**2 + 2*p1*p2*co
  scra4 = (-scra2 + p2*sqrt(scra1)) / scra3
  assk = nn + npoint / twopi * asin(si * scra4)
  tuneabt2 = one - (assk - one) / float(npoint)

end function tuneabt2

double precision function fitlyap(distvect, onelog, turnnumber, nturn, deltax)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Computes interpolated Lyapunov exponent.
  !   DISTVECT (normalized) distance between two companion particles
  !   NTURN is the number of turns.
  !----------------------------------------------------------------------*
  integer :: nturn
  double precision :: distvect(*), onelog(*), turnnumber(*), deltax

  integer :: i
  double precision :: deltalog(6)
  double precision :: turnlog(nturn), test(nturn)
  double precision :: fitlyap2, fitlyap3

  double precision, external :: slopexy
  double precision, parameter :: zero=0d0, one=1d0, dlmax=1d-5

  TURNNUMBER(:nturn) = (/ (dble(i), i = 1, nturn) /)
  TURNLOG(:nturn) = (/ (log(dble(i)), i = 1, nturn) /)
  ONELOG(:nturn) = zero
  do i = 1, nturn
     if (distvect(i) .ne. zero) onelog(i) = log(distvect(i))
  enddo

  !---- Loglog fit over 3 subsequent periods of i = NTURN/4 turns
  !     starting at turn i + 1, i.e. at the second fourth
  i = int(nturn / 4)

  !---- DELTALOG = loglog slope for the last three periods
  deltalog(1) = slopexy(turnnumber(  i + 1), onelog(  i + 1), i)
  deltalog(2) = slopexy(turnnumber(2*i + 1), onelog(2*i + 1), i)
  deltalog(3) = slopexy(turnnumber(3*i + 1), onelog(3*i + 1), i)


  ! 1/turnnumber ln (distvec(i)/distvec(0))

  fitlyap = zero
  !---- if loglog slope is significantly different from zero anywhere
  if ( maxval(DELTALOG(1:3)) .ge. dlmax ) & 
       fitlyap = maxval(DELTALOG(1:3))
  
  !--- DELTALOG = loglog slope for the last three periods
  deltalog(4) = slopexy(turnlog(  i + 1), onelog(  i + 1), i)
  deltalog(5) = slopexy(turnlog(2*i + 1), onelog(2*i + 1), i)
  deltalog(6) = slopexy(turnlog(3*i + 1), onelog(3*i + 1), i)

  fitlyap2 = zero
  !---- if loglog slope is significantly different from zero anywhere
  if ( maxval(DELTALOG(4:6)) + 1.d0 .ge. dlmax ) & 
       fitlyap2 = maxval(DELTALOG(4:6))
  
!  write(69,*) 'deltalogs: ', deltalog(1:6), 'fitlyaps: ', fitlyap, fitlyap2, & 
!       ' nturn and i:', nturn, i
end function fitlyap

subroutine wmaxmin(track,eigen,znt)
  use wmaxmin0fi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes maximum and minimum betatron invariants during tracking.  *
  ! Input:                                                               *
  !   TRACK(6,*)(real)    Track coordinates: (X, PX, Y, PY, T, PT).      *
  !----------------------------------------------------------------------*
  double precision, intent(IN)  :: track(6), eigen(6,6)
  double precision, intent(OUT) :: znt(6)

  integer :: i
  double precision :: wx, wy, wxy

  !---- Convert to normalized values.
  do i = 1, 3 
     znt(2*i-1) = eigen(2,2*i) * track(1) - eigen(1,2*i) * track(2) &
                + eigen(4,2*i) * track(3) - eigen(3,2*i) * track(4) &
                + eigen(6,2*i) * track(5) - eigen(5,2*i) * track(6)
     znt(2*i)   = eigen(1,2*i-1) * track(2) - eigen(2,2*i-1) * track(1) &
                + eigen(3,2*i-1) * track(4) - eigen(4,2*i-1) * track(3) &
                + eigen(5,2*i-1) * track(6) - eigen(6,2*i-1) * track(5)
  enddo

  !---- Convert to amplitudes (and phases: not computed).
  wx = znt(1)*znt(1) + znt(2)*znt(2)
  wy = znt(3)*znt(3) + znt(4)*znt(4)
  wxy = wx + wy

  !---- Compare to and redefine WMIN and WMAX in TRDYNAP.
  wxmax  = max(wx,  wxmax)  ;  wxmin  = min(wx,  wxmin)
  wymax  = max(wy,  wymax)  ;  wymin  = min(wy,  wymin)
  wxymax = max(wxy, wxymax) ;  wxymin = min(wxy, wxymin)

end subroutine wmaxmin

double precision function slopexy(vectorx, vectory, nturn)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes slope a for linear fit of Y = a X + b.                    *
  ! VECTORX  is the array of abscissas X,                                *
  ! VECTORY  is the array of ordinates Y to interpolate,                 *
  ! NTURN    is their dimension.                                         *
  !----------------------------------------------------------------------*
  integer, intent(IN) :: nturn
  double precision, intent(IN) :: vectorx(nturn), vectory(nturn)

  double precision :: x2mean, xmean, xymean, ymean

  xmean = sum(vectorx) / nturn
  ymean = sum(vectory) / nturn
  x2mean = dot_product(vectorx, vectorx) / nturn
  xymean = dot_product(vectorx, vectory) / nturn

  slopexy = (xymean - xmean * ymean) / (x2mean - xmean**2)

end function slopexy

