!----------------------------------------------------------------------*
subroutine dynap(eigen, coords, turns, npart, distvect, zn,dq,    &
     onelog, turnnumber)
  use deltrafi
  implicit none



  integer turns,npart
  double precision eigen(6,6),coords(6,0:turns,*),get_value
  double precision zn(turns,6), distvect(turns),dq(2*turns),        &
       onelog(turns), turnnumber(turns)

  write(*,*) ' entered dynap '

  !      print *, 'eigenvectors'
  !      do i = 1, 6
  !        print '(1p,6e12.4)', (eigen(i,j), j=1,6)
  !      enddo
  !      do k = 1, 2
  !        print *, 'particle: ', k
  !        do i = 1, turns
  !          print '(1p,6e12.4)', (coords(j,i,k), j = 1, 6)
  !        enddo
  !      enddo

  fastune = get_value('dynap ','fastune ') .ne. 0
  deltax = get_value('dynap ','lyapunov ')

  call trdynrun (eigen,coords,turns,npart,distvect,zn,onelog,       &
       turnnumber,dq)
  !      call dynapfill()

  write(*,*) ' end dynap '

end subroutine dynap

!**********************************************************************
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
  !   yapunov   (real)    : interpolated Lyapunov exponent               *
  !----------------------------------------------------------------------*

  call double_to_table_curr('dynap ', 'dynapfrac ', dynapfrac)
  call double_to_table_curr('dynap ', 'dktrturns ', dktrturns)
  call double_to_table_curr('dynap ', 'xend ', xend )
  call double_to_table_curr('dynap ', 'pxend ', pxend )
  call double_to_table_curr('dynap ', 'yend ', yend )
  call double_to_table_curr('dynap ', 'pyend ',pyend )
  call double_to_table_curr('dynap ', 'tend ', ptend )
  call double_to_table_curr('dynap ', 'wxmin ', wxmin )
  call double_to_table_curr('dynap ', 'wxmax ', wxmax )
  call double_to_table_curr('dynap ', 'wymin ', wymin )
  call double_to_table_curr('dynap ', 'wymax ', wymax )
  call double_to_table_curr('dynap ', 'wxymin ', wxymin )
  call double_to_table_curr('dynap ', 'wxymax ', wxymax )
  call double_to_table_curr('dynap ', 'smear ', smear )
  call double_to_table_curr('dynap ', 'yapunov ', yapunov )

  write(*,*) ' yapunov = ', yapunov

  call augment_count('dynap ')
end subroutine dynapfill
!-----------------  end of dynapfill subroutine --------------------------
!**********************************************************************
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
  double precision tx_tmp, ty_tmp
  tx_tmp = tunx
  if (tunx .gt. 0.5d0)  tx_tmp = 1.d0 - tunx
  ty_tmp = tuny
  if (tuny .gt. 0.5d0)  ty_tmp = 1.d0 - tuny
  call double_to_table_curr('dynaptune ', 'x '    , x0)
  call double_to_table_curr('dynaptune ', 'y '    , y0)
  call double_to_table_curr('dynaptune ', 'tunx ', tx_tmp)
  call double_to_table_curr('dynaptune ', 'tuny ', ty_tmp)
  call double_to_table_curr('dynaptune ', 'dtune ', dtune )

  write(*,*) ' tunes = ', tunx, tuny, dtune

  call augment_count('dynaptune ')
end subroutine dynaptunefill
!-----------------  end of dynaptunefill subroutine --------------------------


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine trdynrun (eigen,coords,turns,npart,distvect,zn,onelog, &
     turnnumber,dq)

  use deltrafi
  use wmaxmin0fi
  use dyntabfi
  use tunesfi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !
  !----------------------------------------------------------------------*
  integer k,ix,iy,initt,turns,ktrturns,nturnhalf,i,j,npart
  double precision eigen(6,6),coords(6,0:turns,*),distvect(turns),  &
       zn(turns,6),znt(6),track(6),fitlyap,     &
       templyap,tunx1,tunx2,tuny1,tuny2,tuneabt,tuneabt2,zero,one,two,   &
       onelog(turns), turnnumber(turns),dq(2*turns),dphi(3),dphi1,dphi2, &
       get_variable,pi,twopi
  parameter(zero=0d0,one=1d0,two=2d0)

  pi=get_variable('pi ')
  twopi=get_variable('twopi ')
  dktrturns = turns

  !---- Initialize max and min betatron invariants.
  wxmax = zero
  wymax = zero
  wxymax = zero
  wxmin = zero
  wymin = zero
  wxymin = zero

  !--- distance for second particle with x add-on
  !        deltax = get_value('dynap ', 'lyapunov ')
  !        jend = 2
  !        z(1,jend) = z(1,1) + deltax
  !        do k = 2, 6
  !          z(k,jend) = z(k,1)
  !        enddo

  !---- compute minimum and maximum amplitudes
  !     and normalized coordinates
  !     for all particles
  do k = 1, npart, 2
     x0 = coords(1,0,k)
     y0 = coords(3,0,k)
     do i = 1, turns
        do j = 1, 6
           track(j) = coords(j,i,k)
           !          write(*,*) track(j)
        end do
        call wmaxmin(track,eigen,znt)
        do j = 1, 6
           zn(i,j)=znt(j)
           !          write(*,*) zn(i,j)
        end do
     end do

     !---- compute 'smear'
     smear = two * (wxymax - wxymin) / (wxymax + wxymin)

     !      write(*,*)' smear = ', smear

     !---- Fast tune calculation by interpolated FFT.
     if (fastune) then

        ktrturns = turns

        ix = 1
        iy =3
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
        if (abs(tunx1-tunx).gt.0.4) tunx1 = 1-tunx1
        if (abs(tunx2-tunx).gt.0.4) tunx2 = 1-tunx2
        if (abs(tuny1-tuny).gt.0.4) tuny1 = 1-tuny1
        if (abs(tuny2-tuny).gt.0.4) tuny2 = 1-tuny2
        dtune = sqrt((tunx2 - tunx1)**2 + (tuny2 - tuny1)**2)
        write(*,*) tunx1,tunx2, ' ',tuny1,tuny2
        call dynaptunefill
     endif

     !---- Lyapunov exponent calculation

     open(50,file="lyapunov.data")

     do i = 1, turns
        do j = 1, 6
           track(j) = coords(j,i,k+1)
           !          write(*,*) track(j)
        end do
        call wmaxmin(track,eigen,znt)
        templyap = zero
        do j = 1, 3
           if(znt(2*j).eq.zero.and.znt(2*j-1).eq.zero) then
              dphi1=zero
           else
              dphi1=atan2(znt(2*j),znt(2*j-1))
           endif
           if(zn(i,2*j).eq.zero.and.zn(i,2*j-1).eq.zero) then
              dphi2=zero
           else
              dphi2=atan2(zn(i,2*j),zn(i,2*j-1))
           endif
           dphi(j) =  dphi1-dphi2
           if (dphi(j).gt.pi) dphi(j)=dphi(j)-twopi
           if (dphi(j).lt.-pi) dphi(j)=dphi(j)+twopi
           templyap = templyap + (dphi(j))**2
        end do
        distvect(i) = sqrt(templyap)
        write(50,*) i,distvect(i)
        !            write(*,*) distvect(i)
     end do
     !        write(*,*) ' templyap =', templyap
     yapunov = fitlyap(distvect,onelog,turnnumber,turns)
     call dynapfill()

  enddo
end subroutine trdynrun

!-----------------  end of trdynrun subroutine --------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fft(data, nn, isign)

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
  integer i,isign,istep,j,m,mmax,n,nn
  double precision data(*),tempi,tempr,theta,wi,wpi,wpr,wr,wtemp,   &
       get_variable,zero,half,one,two,twopi
  parameter(zero=0d0,half=0.5d0,one=1d0,two=2d0)

  !---- Initialize
  twopi=get_variable('twopi ')

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

!-----------------  end of fft subroutine --------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function tuneabt(zn, ixy, initt, maxn, turns, dq)

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
  integer mft,nft,ixy,initt,maxn,nftmax,npoint,mf,turns
  !     choose dimensions for now
  double precision zn(turns,6),dq(2*turns),ftmax,temp,cf1,cf2,cf3,  &
       arg,assk,pi,get_variable,zero,one,two
  parameter(zero=0d0,one=1d0,two=2d0)

  !---- Initialize
  pi=get_variable('pi ')

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
!-----------------  end of tuneabt function --------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function tuneabt2(zn, ixy, initt, maxn, turns,dq)

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
  integer mft,nft,mf,ixy,initt,maxn,nftmax,npoint,nn,turns
  !     need dimensions?
  double precision zn(turns,6),dq(2*turns),ftmax,temp,step,cf1,cf2, &
       scra1,scra2,scra3,scra4,assk,co,si,p1,p2,pi,twopi,get_variable,   &
       zero,one,two,cf3
  parameter(zero=0d0,one=1d0,two=2d0)

  !---- Initialize
  pi=get_variable('pi ')
  twopi=get_variable('twopi ')

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
!-----------------  end of tuneabt2 function --------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function fitlyap(distvect, onelog, turnnumber,   &
     nturn)
  implicit none


  !----------------------------------------------------------------------*
  ! Purpose:
  !   Computes interpolated Lyapunov exponent.
  !   DISTVECT (normalized) distance between two companion particles
  !   NTURN is the number of turns.
  !----------------------------------------------------------------------*
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !     needs:
  !      slopexy
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer mf,n1,n2,n3,npoint,nturn
  double precision onelog(*),turnnumber(*),fitlyap1,fitlyap2,       &
       fitlyap3,deltalog1,deltalog2,deltalog3,distvect(*),slopexy,zero,  &
       one,dlmax
  parameter(zero=0d0,one=1d0,dlmax=1d-5)

  do mf = 1, nturn
     !        dq(ilogd+mf) = log(distvect(mf))
     !        dq(in+mf)    = mf
     !        dq(ilogn+mf) = log(dq(in+mf))
     if(distvect(mf).eq.zero) then
        onelog(mf) = zero
     else
        onelog(mf) = log(distvect(mf))
     endif
     turnnumber(mf)    = real(mf)
     !        turnlog(mf) = log(dble(turnnumber(mf)))
  enddo

  !---- Loglog fit over 3 subsequent periods of NPOINT = NTURN/4 TURNS
  !     starting at N1 = NPOINT + 1, i.e. at the second fourth
  npoint = int(nturn / 4)
  n1 = npoint + 1
  n2 = n1 + npoint
  n3 = n2 + npoint

  !---- DELTALOG = deviation from 1 of loglog slope.
  !
  deltalog1 = slopexy(turnnumber(n1), onelog(n1), npoint) - one
  deltalog2 = slopexy(turnnumber(n2), onelog(n2), npoint) - one
  deltalog3 = slopexy(turnnumber(n3), onelog(n3), npoint) - one

  if (deltalog1 .lt. dlmax  .and.  deltalog2 .lt. dlmax  .and.      &
       deltalog3 .lt. dlmax) then
     fitlyap = zero
  else
     fitlyap1 = slopexy(turnnumber(n1), onelog(n1), npoint)
     fitlyap2 = slopexy(turnnumber(n2), onelog(n2), npoint)
     fitlyap3 = slopexy(turnnumber(n3), onelog(n3), npoint)

     if (fitlyap1 .lt. fitlyap2) then
        fitlyap = fitlyap2
     else
        fitlyap = fitlyap1
     endif
     if (fitlyap .lt. fitlyap3) then
        fitlyap = fitlyap3
     endif
  endif

end function fitlyap
!-----------------  end of fitlyap function --------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine wmaxmin(track,eigen,znt)

  use wmaxmin0fi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes maximum and minimum betatron invariants during tracking.  *
  ! Input:                                                               *
  !   TRACK(6,*)(real)    Track coordinates: (X, PX, Y, PY, T, PT).      *
  !----------------------------------------------------------------------*
  integer kp,kq
  double precision track(6),wx,wxy,wy,znt(6),eigen(6,6),zero
  parameter(zero=0d0)

  !---- Copy track coordinates.
  !      do 10 i = 1, 6
  !        z(i) = track(i)
  !   10 continue

  !---- Convert to normalized values.
  do kq = 1, 5, 2
     kp = kq + 1
     znt(kq) = eigen(2,kp) * track(1) - eigen(1,kp) * track(2)       &
          + eigen(4,kp) * track(3) - eigen(3,kp) * track(4)                 &
          + eigen(6,kp) * track(5) - eigen(5,kp) * track(6)
     znt(kp) = eigen(1,kq) * track(2) - eigen(2,kq) * track(1)       &
          + eigen(3,kq) * track(4) - eigen(4,kq) * track(3)                 &
          + eigen(5,kq) * track(6) - eigen(6,kq) * track(5)
  enddo

  !---- Convert to amplitudes (and phases: not computed).
  wx = znt(1)**2 + znt(2)**2
  wy = znt(3)**2 + znt(4)**2
  wxy = wx + wy

  !---- Compare to and redefine WMIN and WMAX in TRDYNAP.
  if (wx.gt.wxmax) then
     wxmax = wx
  else if (wx.lt.wxmin.or.wxmin.eq.zero) then
     wxmin = wx
  endif
  if (wy.gt.wymax) then
     wymax = wy
  else if (wy.lt.wymin.or.wymin.eq.zero) then
     wymin = wy
  endif
  if (wxy.gt.wxymax) then
     wxymax = wxy
  else if (wxy.lt.wxymin.or.wxymin.eq.zero) then
     wxymin = wxy
  endif

end subroutine wmaxmin
!-----------------  end of wmaxmin subroutine --------------------------

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function slopexy(vectorx, vectory, nturn)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Computes slope a for linear fit of Y = a X + b.                    *
  ! VECTORX  is the array of abscissas X,                                *
  ! VECTORX  is the array of ordinates Y to interpolate,                 *
  ! NTURN    is their dimension.                                         *
  !----------------------------------------------------------------------*
  integer mf,nturn
  double precision vectorx(*),vectory(*),x2mean,xmean,xymean,ymean, &
       zero
  parameter(zero=0d0)

  xmean  = zero
  ymean  = zero
  x2mean = zero
  xymean = zero

  do mf = 1, nturn
     xmean  = xmean + vectorx(mf)
     ymean  = ymean + vectory(mf)
     x2mean = x2mean + vectorx(mf)**2
     xymean = xymean + vectorx(mf) * vectory(mf)
  enddo

  xmean  = xmean / nturn
  ymean  = ymean / nturn
  x2mean = x2mean / nturn
  xymean = xymean / nturn

  slopexy = (xymean - xmean * ymean) / (x2mean - xmean**2)

end function slopexy
!-----------------  end of slopexy function --------------------------
