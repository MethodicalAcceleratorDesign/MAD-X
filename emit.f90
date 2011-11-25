subroutine emit(deltap, tol, orbit0, disp0, rt, u0, emit_v,       &
     nemit_v, bmax, gmax, dismax, tunes, sig_v, pdamp)
  use bbfi
  use twiss0fi
  use emitfi
  implicit none


  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Compute emittances by A. Chao's method.                            *
  ! Input:                                                               *
  !   DELTAP     (real)   Average energy error desired.                  *
  !   tol        (real)   tolerance for distinction static/dynamic:
  !                       |e_val_5| < tol & |e-val_6| < tol: static,
  !                       default: 1.000001
  !   orbit0     (real)   closed orbit from TMCLOR
  !   disp0      (real)   dispersion from TMCLOR
  !   rt         (real)   one_turn matrix TMCLOR
  ! Output:
  !   u0         (real)   Radiation loss per turn in GeV
  !   emit_v     (real)   ex, ey, et (emittances)
  !   nemit_v    (real)   exn, eyn, etn (normalised emitt., MAD style)
  !   bmax       (real)   Maximum extents of modes.
  !   gmax       (real)   Maximum divergences of modes.
  !   dismax     (real)   Maximum dispersion.
  !   tunes      (real)   qx, qy, qs
  !   pdamp      (real)   damping partition numbers
  !   sig_v      (real)   sigx, sigy, sigt, sige
  !----------------------------------------------------------------------*
  !---- Communication area for radiation damping.
  double precision orbit0(6), orbit(6), orbit2(6), em(6,6), rd(6,6), reval(6)
  double precision aival(6), rt(6,6), orbit1(6), ek(6)
  double precision emit_v(3), nemit_v(3), tunes(3), sig_v(4)
  double precision u0, pdamp(3)
  double precision bmax(3,3), gmax(3,3), dismax(4)
  double precision em2(6,6), disp(6), disp0(6),al_errors(align_max)
  double precision tol, deltap, get_value, arad, suml, gammas, el
  double precision betas, bx, gx, re(6,6), te(6,6,6), tt(6,6,6)
  double precision get_variable, node_value
  integer i, j, j1, j2, k, k1, k2, eflag, restart_sequ, n_align
  integer node_al_errors, code, advance_node
  logical m66sta, fmap, stabx, staby, stabt, frad
  double precision zero, one, three, twopi
  parameter (zero = 0.0d0, one = 1.d0, three = 3.0d0)

  twopi = get_variable('twopi ')
  do i = 1, 6
     orbit(i) = orbit0(i)
     disp(i) = disp0(i)
  enddo
  call dzero(emit_v, 3)
  call dzero(nemit_v, 3)
  call dzero(bmax, 9)
  call dzero(gmax, 9)
  call dzero(dismax, 4)
  call dzero(tunes, 3)
  call dzero(sig_v, 4)
  call dzero(pdamp, 3)
  u0 = zero
  !---- Find eigenvectors at initial position.
  if (m66sta(rt)) then
     call laseig(rt, reval, aival, em)
     stabt = .false.
     !        print '('' Static map, eigenvalues:'',(/1X,2E15.8))',
     !     +    (reval(i), aival(i), i = 1, 4)
  else
     call ladeig(rt, reval, aival, em)
     stabt = reval(5)**2 + aival(5)**2 .le. tol  .and.               &
          reval(6)**2 + aival(6)**2 .le. tol
     !        print '('' Static map, eigenvalues:'',(/1X,2E15.8))',
     !     +    (reval(i), aival(i), i = 1, 6)
  endif
  stabx = reval(1)**2 + aival(1)**2 .le. tol  .and.                 &
       reval(2)**2 + aival(2)**2 .le. tol
  staby = reval(3)**2 + aival(3)**2 .le. tol  .and.                 &
       reval(4)**2 + aival(4)**2 .le. tol

  !---- Maximum extents.
  do j = 1, 3
     j1 = 2 * j -1
     j2 = 2 * j
     do k = 1, 3
        k1 = 2 * k - 1
        k2 = 2 * k
        bmax(j,k) = em(j1,k1) * em(j1,k1) + em(j1,k2) * em(j1,k2)
        gmax(j,k) = em(j2,k1) * em(j2,k1) + em(j2,k2) * em(j2,k2)
     enddo
  enddo
  arad = get_value('probe ','arad ')
  betas = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')
  cg = arad * gammas**3 / three
  frad = get_value('probe ','radiate ') .ne. zero
  !---- Initialize damping calculation.
  if (frad .and. stabt) then
     sum(1) = zero
     sum(2) = zero
     sum(3) = zero
     sumu0 = zero
     call m66one(rd)
  endif
  call dzero(tt,216)
  call m66one(rt)
  eflag = 0
  suml = zero
  bbd_cnt=0
  bbd_flag=1

  i = restart_sequ()
10 continue
  bbd_pos=i
  code = node_value('mad8_type ')
  if(code.eq.39) code=15
  if(code.eq.38) code=24
  el = node_value('l ')
  n_align = node_al_errors(al_errors)
  if (n_align.ne.0)  then
     call dcopy(orbit, orbit2, 6)
     call tmali1(orbit2,al_errors,betas,gammas,orbit,re)
     if (.not. stabt) call m66byv(re, disp, disp)
     call m66mpy(re, em, em)
     if (frad .and. stabt) call m66mpy(re, rd, rd)
  endif
  !*---- Keep orbit at entrance.
  call dcopy(orbit, orbit1, 6)
  !---- Element matrix and length.
  call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te)
  if (fmap) then
     !---- Advance dispersion.
     if (.not. stabt) then
        call m66byv(re, disp, disp)
        do j = 1, 4
           dismax(j) = max(abs(disp(j)),dismax(j))
        enddo
     endif
     !---- Radiation damping.
     call m66mpy(re, em, em2)
     if (frad .and. stabt) then
        call emdamp(code, deltap, em, em2, orbit1, orbit, re)
        call m66mpy(re, rd, rd)
     endif
     call dcopy(em2, em, 36)

     !---- Compute beta and gamma.
     do j = 1, 3
        j1 = 2 * j -1
        j2 = 2 * j
        do k = 1, 3
           k1 = 2 * k - 1
           k2 = 2 * k
           bx = em(j1,k1) * em(j1,k1) + em(j1,k2) * em(j1,k2)
           bmax(j,k) = max(bmax(j,k),bx)
           gx = em(j2,k1) * em(j2,k1) + em(j2,k2) * em(j2,k2)
           gmax(j,k) = max(gmax(j,k),gx)
        enddo
     enddo
     suml = suml + el
  endif
  if (n_align.ne.0)  then
     call dcopy(orbit, orbit2, 6)
     call tmali2(el,orbit2,al_errors,betas,gammas,orbit,re)
     if (.not. stabt) call m66byv(re, disp, disp)
     call m66mpy(re, em, em)
     if (frad .and. stabt) call m66mpy(re, rd, rd)
  endif
  if (advance_node().ne.0)  then
     i=i+1
     goto 10
  endif
  bbd_flag=0
  !---- Undamped tunes and beam extents.
  qx = atan2(aival(1), reval(1)) / twopi
  if (qx .lt. zero) qx = qx + one
  qy = atan2(aival(3), reval(3)) / twopi
  if (qy .lt. zero) qy = qy + one
  qs = atan2(aival(5), reval(5)) / twopi
  if (qs .lt. zero) qs = - qs
  !---- Summary output.
  call emsumm(rd,em,bmax,gmax,stabt,frad,u0,emit_v,nemit_v,tunes,   &
       sig_v,pdamp)
end subroutine emit
subroutine emsumm(rd,em,bmax,gmax,stabt,frad,u0,emit_v,nemit_v,   &
     tunes,sig_v,pdamp)
  use emitfi
  implicit none


  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Finish radiation damping calculation and print summary.            *
  ! Input:                                                               *
  !   RD(6,6)   (real)    Damped one-turn transfer matrix.               *
  !   EM(6,6)   (real)    Undamped eigenvectors.                         *
  !   BMAX(3,3) (real)    Maximum extents of modes.                      *
  !   GMAX(3,3) (real)    Maximum divergences of modes.                  *
  !   STABT     (logical) anti-STATIC flag
  !   FRAD      (logical) radiation flag
  ! Output:
  !   u0         (real)   Radiation loss per turn in GeV
  !   emit_v     (real)   ex, ey, et (emittances)
  !   nemit_v    (real)   exn, eyn, etn (normalised emitt., MAD style)
  !   tunes      (real)   qx, qy, qs
  !   pdamp      (real)   damping partition numbers
  !   sig_v      (real)   sigx, sigy, sigt, sige
  !----------------------------------------------------------------------*
  integer j,j1,j2,k,k1,k2,iqpr2
  double precision arad,gammas,clg,etpr,expr,eypr,f0,sal,en0
  double precision amass, clight, hbar, freq0, u0, betas
  double precision ten3p,tenp6,tenp9,three,twopi,one,two,zero,four
  double precision ex,ey,et,exn,eyn,sigx,sigy,sige,sigt
  double precision rd(6,6), em(6,6), bmax(3,3), gmax(3,3),pdamp(3)
  double precision get_variable,get_value
  double precision emit_v(3), nemit_v(3), tunes(3), sig_v(4)
  logical stabt, frad
  double precision reval(6), aival(6), alj(3), tau(3), tune(3)
  double precision sigma(6,6), bstar(3,3), gstar(3,3), dummy(6,6)
  parameter (zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.0d0)
  parameter (iqpr2 = 6, four = 4.d0)
  parameter (ten3p = 1.0d3, tenp6 = 1.0d6, tenp9 = 1.0d9)

  call dzero(sigma,36)
  ex=zero
  ey=zero
  et=zero
  twopi=get_variable('twopi ')
  clight = get_variable('clight ')
  hbar = get_variable('hbar ')
  arad = get_value('probe ','arad ')
  betas = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')
  amass = get_value('probe ','mass ')
  freq0 = get_value('probe ','freq0 ')
  !---- Synchrotron energy loss [GeV].
  if (stabt .and. frad) then
     u0 = sumu0

     !---- Tunes.
     call ladeig(rd, reval, aival, dummy)
     tune(1) = atan2(aival(1), reval(1)) / twopi
     if (tune(1) .lt. zero) tune(1) = tune(1) + one
     tune(2) = atan2(aival(3), reval(3)) / twopi
     if (tune(2) .lt. zero) tune(2) = tune(2) + one
     tune(3) = atan2(aival(5), reval(5)) / twopi
     if (tune(3) .lt. zero) tune(3) = - tune(3)

     !---- Damping constants per turn.
     alj(1) = - log(reval(1)**2 + aival(1)**2) / two
     alj(2) = - log(reval(3)**2 + aival(3)**2) / two
     alj(3) = - log(reval(5)**2 + aival(5)**2) / two

     !---- Damping partition numbers.
     en0 = get_value('probe ','energy ')
     sal = two * en0 / u0
     pdamp(1) = alj(1) * sal
     pdamp(2) = alj(2) * sal
     pdamp(3) = alj(3) * sal
     !---- Emittances.
     clg = ((55.d0 * hbar * clight) / (96.d0 * sqrt(three)))         &
          * ((arad * gammas**5) / amass)
     ex = clg * sum(1) / alj(1)
     ey = clg * sum(2) / alj(2)
     et = clg * sum(3) / alj(3)

     !---- Damping constants per second and damping times.
     f0 = freq0 * tenp6
     alj(1) = abs(alj(1) * f0)
     alj(2) = abs(alj(2) * f0)
     alj(3) = abs(alj(3) * f0)
     tau(1) = one / alj(1)
     tau(2) = one / alj(2)
     tau(3) = one / alj(3)
  endif

  !---- TRANSPORT sigma matrix.
  call emce2i(stabt, em, ex, ey, et, sigma)

  !---- Extents at interaction point.
  do j = 1, 3
     j1 = 2 * j -1
     j2 = 2 * j
     do k = 1, 3
        k1 = 2 * k - 1
        k2 = 2 * k
        bstar(j,k) = em(j1,k1) * em(j1,k1) + em(j1,k2) * em(j1,k2)
        gstar(j,k) = em(j2,k1) * em(j2,k1) + em(j2,k2) * em(j2,k2)
     enddo
  enddo

  exn = ex * four * betas * gammas
  eyn = ey * four * betas * gammas
  sigx = sqrt(abs(sigma(1,1)))
  sigy = sqrt(abs(sigma(3,3)))
  if (sigma(5,5) .gt. zero .or. sigma(6,6) .gt. zero)  then
     sigt = sqrt(abs(sigma(5,5)))
     sige = sqrt(abs(sigma(6,6)))
  else
     sigt = zero
     sige = zero
  endif
  tunes(1) = qx
  tunes(2) = qy
  tunes(3) = qs
  emit_v(1) = ex
  emit_v(2) = ey
  emit_v(3) = et
  nemit_v(1) = exn
  nemit_v(2) = eyn
  sig_v(1) = sigx
  sig_v(2) = sigy
  sig_v(3) = sigt
  sig_v(4) = sige
  !---- Summary output; header and global parameters.
  !---- Dynamic case.
  expr = ex * tenp6
  eypr = ey * tenp6
  etpr = et * tenp6
  if (stabt) then
     if (frad) write (iqpr2, 910) ten3p * u0
     write (iqpr2, 920) 1, 2, 3
     write (iqpr2, 930) qx, qy, qs
     if (frad) write (iqpr2, 940) tune
     write (iqpr2, 950) ((bstar(j,k), j = 1, 3), k = 1, 3),          &
          ((gstar(j,k), j = 1, 3), k = 1, 3),                               &
          ((bmax(j,k), j = 1, 3), k = 1, 3),                                &
          ((gmax(j,k), j = 1, 3), k = 1, 3)
     if (frad) then
        write (iqpr2, 960) pdamp, alj, (tau(j), j = 1, 3),            &
             expr, eypr, etpr
     endif
  else
     write (iqpr2, 920) 1, 2
     write (iqpr2, 930) qx, qy
     write (iqpr2, 970) ((bstar(j,k), j = 1, 2), k = 1, 2),          &
          ((gstar(j,k), j = 1, 2), k = 1, 2),                               &
          ((bmax(j,k), j = 1, 2), k = 1, 2),                                &
          ((gmax(j,k), j = 1, 2), k = 1, 2)
  endif

  !---- RF system.
  !      call enprrf

910 format(t6,'U0',t16,f14.6,' [MeV/turn]')
920 format(' '/' ',t42,3(9x,'M o d e',3x,i1:))
930 format(' Fractional tunes',t30,'undamped',t42,3f20.8)
940 format(' ',t30,'damped',t42,3f20.8)
950 format(' '/' beta* [m]',t30,'x',t42,3e20.8/t30,'y',t42,3e20.8/    &
       t30,'t',t42,3e20.8/                                               &
       ' '/' gamma* [1/m]',t30,'px',t42,3e20.8/t30,'py',t42,3e20.8/      &
       t30,'pt',t42,3e20.8/                                              &
       ' '/' beta(max) [m]',t30,'x',t42,3e20.8/t30,'y',t42,3e20.8/       &
       t30,'t',t42,3e20.8/                                               &
       ' '/' gamma(max) [1/m]',t30,'px',t42,3e20.8/t30,'py',t42,3e20.8/  &
       t30,'pt',t42,3e20.8)
960 format(' '/' Damping partition numbers',t42,3f20.8/               &
       ' Damping constants [1/s]',t46,3e20.8/                            &
       ' Damping times [s]',t46,3e20.8/                                  &
       ' Emittances [pi micro m]',t42,3e20.8)
970 format(' '/' beta* [m]',t30,'x',t42,2e20.8/t30,'y',t42,2e20.8/    &
       ' '/' gamma* [1/m]',t30,'px',t42,2e20.8/t30,'py',t42,2e20.8/      &
       ' '/' beta(max) [m]',t30,'x',t42,2e20.8/t30,'y',t42,2e20.8/       &
       ' '/' gamma(max) [1/m]',t30,'px',t42,2e20.8/t30,'py',t42,2e20.8)

end subroutine emsumm
subroutine emce2i(stabt, em, ex, ey, et, sigma)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Convert eigenvectors to internal sigma matrix form.                *
  ! Input:                                                               *
  !   EM(6,6)   (real)    Eigenvector matrix.                            *
  !   EX        (real)    Horizontal emittance.                          *
  !   EY        (real)    Vertical emittance.                            *
  !   ET        (real)    Longitudinal emittance.                        *
  ! Output:                                                              *
  !   SIGMA(6,6)(real)    Beam matrix in internal form.                  *
  !----------------------------------------------------------------------*
  integer j,k
  double precision et, ex, ey, em(6,6), sigma(6,6)
  logical stabt

  do j = 1, 6
     do k = 1, 6
        sigma(j,k) =                                                  &
             ex * (em(j,1) * em(k,1) + em(j,2) * em(k,2)) +                    &
             ey * (em(j,3) * em(k,3) + em(j,4) * em(k,4))
        if (stabt) then
           sigma(j,k) = sigma(j,k) +                                   &
                et * (em(j,5) * em(k,5) + em(j,6) * em(k,6))
        endif
     enddo
  enddo
end subroutine emce2i
subroutine emdamp(code, deltap, em1, em2, orb1, orb2, re)
  use twiss0fi
  use emitfi
  use twtrrfi
  implicit none


  !---------------------------------------------------------------------*
  ! Purpose:                                                            *
  !   Deal with radiation damping in an element.                        *
  ! Input:                                                              *
  !   code      (int)     MAD-8 element code                            *
  !   deltap    (real)    momentum error                                *
  !   EM1(6,6)  (real)    Matrix of eigenvectors at entrance.           *
  !   EM2(6,6)  (real)    Matrix of eigenvectors at exit.               *
  !   ORB1(6)   (real)    Orbit position at entrance.                   *
  !   ORB2(6)   (real)    Orbit position at exit.                       *
  ! Input/output:                                                       *
  !   RE(6,6)   (real)    Transfer matrix for the element; changed on   *
  !                       output to contain damping.                    *
  !---------------------------------------------------------------------*
  integer code
  double precision deltap
  double precision em1(6,6), em2(6,6), orb1(6), orb2(6), re(6,6)
  double precision ten3m, ten6p, zero, half, one, two, three, four
  double precision six, twelve
  parameter         (ten3m = 1.0d-3, ten6p = 1.0d+6)
  parameter         (zero  = 0.0d0,  half  = 0.5d0)
  parameter         (one   = 1.0d0,  two   = 2.0d0)
  parameter         (three = 3.0d0,  twelve = 12.d0)
  parameter         (four  = 4.0d0,  six   = 6.0d0)
  integer i, j, ir, ii, n, n_ferr, iord, nn, ns, nd, nord
  integer node_fd_errors
  double precision  rw(6,6), tw(6,6,6), ferror(2), rw0(6,6)
  double precision  normal(0:maxmul), skew(0:maxmul)
  double precision  vals(2,0:maxmul), field(2,0:maxmul)
  double precision  f_errors(0:50)
  double precision  o1(6), e1(6,6), o2(6), e2(6,6)
  double precision  x1, y1, t1, x2, y2, t2, px1, py1, pt1
  double precision  px2, py2, pt2
  equivalence       (x1, o1(1)), (px1, o1(2))
  equivalence       (y1, o1(3)), (py1, o1(4))
  equivalence       (t1, o1(3)), (pt1, o1(4))
  equivalence       (x2, o2(1)), (px2, o2(2))
  equivalence       (y2, o2(3)), (py2, o2(4))
  equivalence       (t2, o2(3)), (pt2, o2(4))
  double precision  el, node_value, tilt, bvk
  double precision  edg1, edg2, sk1, sk2, hgap, fint, sks, h, ct
  double precision  corr, hx, hy, hxx, hxy, hyy, h1, hcb1, hcbs1
  double precision  tedg1, fact1, fact1x, rfac1, rfac1x, rfac1y
  double precision  h2, hcb2, tedg2, fact2, fact2x, rfac2
  double precision  rfac2x, rfac2y, bi2gi2, betas, gammas
  double precision  get_value, e5sq1, e5sq2, e5sqs1, e5sqs2, x, y
  double precision  f1, f2, f1s, f2s, twon, str, st, pi, clight
  double precision  r1sq, r2sq, fh1, fh2, dr, di, drt, hcb
  double precision  rfv, rff, rfl, rfvlt, rffrq, rflag, time
  double precision  xkick, ykick, dpx, dpy, an, hyx, hcbs2,hbi
  double precision  sk3, rfac, rfacx, rfacy
  double precision  get_variable, fh
  !      double precision  sk0, sk0s
  if (code .eq. 8)  then
     !--- multipole
     el = node_value('lrad ')
  else
     el = node_value('l ')
  endif
  if (el .eq. zero.and.code .ne. 10) goto 500
  betas = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')
  !---- Prepare data.
  bvk = node_value('other_bv ')
  !---- Switch on element type.
  go to (500,  20,  30, 500,  50,  60,  70,  80, 500, 100,          &
       500, 500, 500, 140, 150, 160, 500, 500, 500, 500,                 &
       500, 500, 500, 500, 500, 500, 500, 500, 500, 500,                 &
       500, 500, 500, 500, 500, 500, 500, 500, 150, 500), code
  go to 500

  !---- Dipole.
20 continue
30 continue

  ! FRS 16.11.2004 This is still in the intermediate spirit of k0&k0s
  !      an = node_value('angle ')
  !      sk0 = an / el
  !      sk0s = node_value('k0s ')
  !      if (sk0s .eq. zero)  then
  !        tilt = zero
  !      else
  !        tilt = atan2(sk0s, sk0)
  !        sk0 = sqrt(sk0**2 + sk0s**2)
  !        an = sk0 * el
  !      endif
  !      an = bvk * an
  !      sk0 = bvk * sk0
  !      sk0s = bvk * sk0s
  ! FRS 16.11.2004 here the correction
  an = bvk * node_value('angle ') * el/node_value('l ')
  tilt = -node_value('tilt ')
  edg1 = bvk * node_value('e1 ')
  edg2 = bvk * node_value('e2 ')
  sk1 = bvk * node_value('k1 ')
  sk2 = bvk * node_value('k2 ')
  hgap = node_value('hgap ')
  fint = node_value('fint ')
  sks = zero
  h = an / el

  !---- Refer orbit and eigenvectors to magnet midplane.
  ct = cos(tilt)
  st = sin(tilt)
  o1(1) =   ct * orb1(1) + st * orb1(3)
  o1(2) =   ct * orb1(2) + st * orb1(4)
  o1(3) = - st * orb1(1) + ct * orb1(3)
  o1(4) = - st * orb1(2) + ct * orb1(4)
  o1(5) = orb1(5)
  o1(6) = orb1(6)
  o2(1) =   ct * orb2(1) + st * orb2(3)
  o2(2) =   ct * orb2(2) + st * orb2(4)
  o2(3) = - st * orb2(1) + ct * orb2(3)
  o2(4) = - st * orb2(2) + ct * orb2(4)
  o2(5) = orb2(5)
  o2(6) = orb2(6)
  do i = 1, 6
     e1(1,i) =   ct * em1(1,i) + st * em1(3,i)
     e1(2,i) =   ct * em1(2,i) + st * em1(4,i)
     e1(3,i) = - st * em1(1,i) + ct * em1(3,i)
     e1(4,i) = - st * em1(2,i) + ct * em1(4,i)
     e1(5,i) = em1(5,i)
     e1(6,i) = em1(6,i)
     e2(1,i) =   ct * em2(1,i) + st * em2(3,i)
     e2(2,i) =   ct * em2(2,i) + st * em2(4,i)
     e2(3,i) = - st * em2(1,i) + ct * em2(3,i)
     e2(4,i) = - st * em2(2,i) + ct * em2(4,i)
     e2(5,i) = em2(5,i)
     e2(6,i) = em2(6,i)
  enddo
  !---- Move through orbit through fringing field;
  !     Requested components of eigenvectors are not affected.
  corr = (h + h) * hgap * fint
  call m66one(rw)
  call dzero(tw,216)
  call tmfrng(.false.,h,sk1,edg1,zero,+one,corr,rw,tw)
  call m66byv(rw,o1,o1)
  call m66one(rw)
  call dzero(tw,216)
  call tmfrng(.false.,h,sk1,edg2,zero,-one,corr,rw,tw)
  call dcopy(rw,rw0,36)
  call m66inv(rw0,rw)
  call m66byv(rw,o2,o2)

  !---- Local curvature and its derivatives,
  !     Coefficients for damping matrix.
  hx = sk1*x1 + sks*y1 + h + half*sk2 * (x1**2 - y1**2)
  hy = sks*x1 - sk1*y1 - sk2*x1*y1
  hxx = sk1 + sk2*x1
  hxy = sks - sk2*y1
  hyx = hxy
  hyy = - hxx
  h1 = sqrt(hx**2 + hy**2)
  hcb1 = h1**3
  hcbs1 = three*h1 *                                                &
       (hx * (hxx*px1 + hxy*py1) + hy * (hxy*px1 + hyy*py1))

  tedg1  = tan(edg1)
  fact1  = (one + h*x1) * (one - tedg1*x1)
  fact1x = h - tedg1 - 2.0*h*tedg1*x1
  rfac1  = cg*el*h1**2*fact1
  rfac1x = cg*el * (two*(hx*hxx+hy*hyx)*fact1 + h1**2*fact1x)
  rfac1y = cg*el *  two*(hx*hxy+hy*hyy)*fact1

  hx = sk1*x2 + sks*y2 + h + half*sk2 * (x2**2 - y2**2)
  hy = sks*x2 - sk1*y2 - sk2*x2*y2
  hxx = sk1 + sk2*x2
  hxy = sks - sk2*y2
  hyx = hxy
  hyy = - hxx
  h2 = sqrt(hx**2 + hy**2)
  hcb2 = h2**3
  hcbs2 = three*h2 *                                                &
       (hx * (hxx*px2 + hxy*py2) + hy * (hxy*px2 + hyy*py2))

  tedg2  = tan(edg2)
  fact2  = (one + h*x2) * (one - tedg2*x2)
  fact2x = h - tedg2 - 2.0*h*tedg2*x2
  rfac2  = cg*el*h2**2*fact2
  rfac2x = cg*el * (two*(hx*hxx+hy*hyx)*fact2 + h2**2*fact2x)
  rfac2y = cg*el *  two*(hx*hxy+hy*hyy)*fact2

  !---- Cubic integration over h**3 * E(i,5) * conjg(E(i,5)).
  bi2gi2 = one / (betas * gammas)**2
  hbi = h / betas
  do i = 1, 3
     ir = 2 * i - 1
     ii = 2 * i

     !---- E(i,5) * conjg(E(i,5)) and its derivative w.r.t. S.
     e5sq1 = e1(5,ir)**2 + e1(5,ii)**2
     e5sq2 = e2(5,ir)**2 + e2(5,ii)**2
     e5sqs1 = two * (e1(5,ir) * (bi2gi2*e1(6,ir) - hbi*e1(1,ir))     &
          + e1(5,ii) * (bi2gi2*e1(6,ii) - hbi*e1(1,ii)))
     e5sqs2 = two * (e2(5,ir) * (bi2gi2*e2(6,ir) - hbi*e2(1,ir))     &
          + e2(5,ii) * (bi2gi2*e2(6,ii) - hbi*e2(1,ii)))

     !---- Integrand and its derivative w.r.t. S.
     f1 = hcb1 * e5sq1
     f2 = hcb2 * e5sq2
     f1s = hcbs1 * e5sq1 + hcb1 * e5sqs1
     f2s = hcbs2 * e5sq2 + hcb2 * e5sqs2

     !---- Actual integration.
     sum(i) = sum(i) + half * el * (f1 + f2) -                       &
          el**2 * (f2s - f1s) / twelve
  enddo
  go to 77

  !---- Quadrupole.
50 continue
  bvk = node_value('other_bv ')
  sk1 = bvk * node_value('k1 ')
  str  = sk1
  n    = 1
  twon = two
  go to 75

  !---- Sextupole.
60 continue
  bvk = node_value('other_bv ')
  sk2 = bvk * node_value('k2 ')
  str  = sk2 / two
  n    = 2
  twon = four
  go to 75

  !---- Octupole.
70 continue
  bvk = node_value('other_bv ')
  sk3 = bvk * node_value('k2 ')
  str  = sk3 / six
  n    = 3
  twon = six

  !---- Common to all pure multipoles.
75 continue
  call dcopy(orb1, o1, 6)
  call dcopy(orb2, o2, 6)
  call dcopy(em1, e1, 36)
  call dcopy(em2, e2, 36)

  !---- Local curvature.
  r1sq = orb1(1)**2 + orb1(3)**2
  r2sq = orb2(1)**2 + orb2(3)**2
  h1 = abs(str) * sqrt(r1sq)**n
  h2 = abs(str) * sqrt(r2sq)**n
  rfac = cg * str**2 * el
  rfac1 = rfac * r1sq**n
  rfac2 = rfac * r2sq**n
  rfac1x = twon * rfac * r1sq**(n-1) * x1
  rfac2x = twon * rfac * r1sq**(n-1) * x2
  rfac1y = twon * rfac * r1sq**(n-1) * y1
  rfac2y = twon * rfac * r1sq**(n-1) * y2

  !---- Trapezoidal integration over h**3 * E(k,5) * conjg(E(k,5)).
  fh1 = half * el * h1**3
  fh2 = half * el * h2**3
  sum(1) = sum(1) + fh1 * (e1(5,1)**2 + e1(5,2)**2)                 &
       + fh2 * (e2(5,1)**2 + e2(5,2)**2)
  sum(2) = sum(2) + fh1 * (e1(5,3)**2 + e1(5,4)**2)                 &
       + fh2 * (e2(5,3)**2 + e2(5,4)**2)
  sum(3) = sum(3) + fh1 * (e1(5,5)**2 + e1(5,6)**2)                 &
       + fh2 * (e2(5,5)**2 + e2(5,6)**2)

  !---- Damping matrices.
  !     Code common to bending magnet and pure multipoles.
77 call m66one(rw)
  rw(2,1) =     - rfac1x * (one + pt1) * px1
  rw(2,2) = one - rfac1  * (one + pt1)
  rw(2,3) =     - rfac1y * (one + pt1) * px1
  rw(2,6) =     - rfac1                * px1
  rw(4,1) =     - rfac1x * (one + pt1) * py1
  rw(4,3) =     - rfac1y * (one + pt1) * py1
  rw(4,4) = one - rfac1  * (one + pt1)
  rw(4,6) =     - rfac1                * py1
  rw(6,1) =     - rfac1x * (one + pt1)**2
  rw(6,3) =     - rfac1y * (one + pt1)**2
  rw(6,6) = one - two * rfac1 * (one + pt1)
  call m66mpy(re, rw, re)

  call m66one(rw)
  rw(2,1) =     - rfac2x * (one + pt2) * px2
  rw(2,2) = one - rfac2  * (one + pt2)
  rw(2,3) =     - rfac2y * (one + pt2) * px2
  rw(2,6) =     - rfac2                * px2
  rw(4,1) =     - rfac2x * (one + pt2) * py2
  rw(4,3) =     - rfac2y * (one + pt2) * py2
  rw(4,4) = one - rfac2  * (one + pt2)
  rw(4,6) =     - rfac2                * py2
  rw(6,1) =     - rfac2x * (one + pt2)**2
  rw(6,3) =     - rfac2y * (one + pt2)**2
  rw(6,6) = one - two * rfac2 * (one + pt2)
  call m66mpy(rw, re, re)
  go to 500

  !---- Thin multipoles, EL is the fictitious length for radiation.
80 continue
  !---- Multipole components.
  call dzero(f_errors,maxferr+1)
  n_ferr = node_fd_errors(f_errors)
  bvk = node_value('other_bv ')
  call dzero(normal,maxmul+1)
  call dzero(skew,maxmul+1)
  call get_node_vector('knl ',nn,normal)
  call get_node_vector('ksl ',ns,skew)
  call dzero(vals,2*(maxmul+1))
  do iord = 0, nn
     vals(1,iord) = normal(iord)
  enddo
  do iord = 0, ns
     vals(2,iord) = skew(iord)
  enddo

  !---- Field error vals.
  call dzero(field,2*(maxmul+1))
  if (n_ferr .gt. 0) then
     call dcopy(f_errors,field,n_ferr)
  endif
  nd = 2 * max(nn, ns, n_ferr/2-1)

  !---- Other components and errors.
  nord = 0
  do iord = 0, nd/2
     do j = 1, 2
        field(j,iord) = bvk * (vals(j,iord) + field(j,iord))          &
             / (one + deltap)
        if (field(j,iord) .ne. zero)  nord = iord
     enddo
  enddo

  !---- Track orbit.
  x = orb1(1)
  y = orb1(3)

  !---- Multipole kick.
  dr = zero
  di = zero
  do iord = nord, 0, -1
     drt = (dr * x - di * y) / float(iord+1) + field(1,iord)
     di  = (dr * y + di * x) / float(iord+1) + field(2,iord)
     dr  = drt
  enddo

  !---- H is local "curvature" due to multipole kick.
  h  = sqrt(dr**2 + di**2) / el
  hcb = half * el * h**3
  sum(1)  = sum(1) + hcb *                                          &
       (em1(5,1)**2 + em1(5,2)**2 + em2(5,1)**2 + em2(5,2)**2)
  sum(2)  = sum(2) + hcb *                                          &
       (em1(5,3)**2 + em1(5,4)**2 + em2(5,3)**2 + em2(5,4)**2)
  sum(3)  = sum(3) + hcb *                                          &
       (em1(5,5)**2 + em1(5,6)**2 + em2(5,5)**2 + em2(5,6)**2)

  !---- Damping matrix, is the same at both ends.
  rfac  = cg * (dr**2 + di**2) / el
  rfacx = cg * (- dr * re(2,1) + di * re(4,1)) / el
  rfacy = cg * (- dr * re(2,3) + di * re(4,3)) / el

  call m66one(rw)
  rw(2,1) = - rfacx * (one + orb1(6)) * orb1(2)
  rw(2,2) = one - rfac * (one + orb1(6))
  rw(2,3) = - rfacy * (one + orb1(6)) * orb1(2)
  rw(2,6) = - rfac * orb1(2)
  rw(4,1) = - rfacx * (one + orb1(6)) * orb1(4)
  rw(4,3) = - rfacy * (one + orb1(6)) * orb1(4)
  rw(4,4) = one - rfac * (one + orb1(6))
  rw(4,6) = - rfac * orb1(4)
  rw(6,1) = - rfacx * (one + orb1(6))
  rw(6,3) = - rfacy * (one + orb1(6))
  rw(6,6) = one - two * rfac * (one + orb1(6))
  call m66mpy(re, rw, re)
  call m66mpy(rw, re, re)
  go to 500

  !---- RF cavities.
100 continue
  rfv = node_value('volt ')
  rff = node_value('freq ')
  rfl = node_value('lag ')
  pi = get_variable('pi ')
  clight = get_variable('clight ')
  rfvlt = ten3m * rfv
  rffrq = rff * (ten6p * two * pi / clight)
  rflag = two * pi * rfl
  time = half * (orb1(5) + orb2(5))
  sumu0 = sumu0 + rfvlt * sin(rflag - rffrq * time)
  go to 500

  !---- Orbit correctors.
140 continue
150 continue
160 continue
  n_ferr = node_fd_errors(f_errors)
  do i = 1, 2
     ferror(i) = zero
  enddo
  if (n_ferr .gt. 0) call dcopy(f_errors, ferror, min(2, n_ferr))
  if(code.eq.14) then
     xkick=bvk*(node_value('kick ')+node_value('chkick ')+           &
          ferror(1))
     ykick=zero
  else if(code.eq.15.or.code.eq.39) then
     xkick=bvk*(node_value('hkick ')+node_value('chkick ')+          &
          ferror(1))
     ykick=bvk*(node_value('vkick ')+node_value('cvkick ')+          &
          ferror(2))
  else if(code.eq.16) then
     xkick=zero
     ykick=bvk*(node_value('kick ')+node_value('cvkick ')+           &
          ferror(2))
  else
     xkick=zero
     ykick=zero
  endif
  !---- Sum up total kicks.
  dpx = xkick / (one + deltap)
  dpy = ykick / (one + deltap)

  !---- Local curvature.
  hx = abs(dpx) / el
  hy = abs(dpy) / el
  rfac = cg * (hx**2 + hx**2) * el

  !---- Trapezoidal integration over h**3*E(k,5)*E*(k,5).
  fh = half * el * sqrt(hx**2 + hy**2)**3
  sum(1) = sum(1) + fh *                                            &
       (em1(5,1)**2 + em1(5,2)**2 + em2(5,1)**2 + em2(5,2)**2)
  sum(2) = sum(2) + fh *                                            &
       (em1(5,3)**2 + em1(5,4)**2 + em2(5,3)**2 + em2(5,4)**2)
  sum(3) = sum(3) + fh *                                            &
       (em1(5,5)**2 + em1(5,6)**2 + em2(5,5)**2 + em2(5,6)**2)

  !---- Damping matrices.
  call m66one(rw)
  rw(2,2) = one - rfac * (one + orb1(6))
  rw(2,6) = - rfac * orb1(2)
  rw(4,4) = one - rfac * (one + orb1(6))
  rw(4,6) = - rfac * orb1(4)
  rw(6,6) = one - two * rfac * (one + orb1(6))
  call m66mpy(re, rw, re)

  call m66one(rw)
  rw(2,2) = one - rfac * (one + orb2(6))
  rw(2,6) = - rfac * orb2(2)
  rw(4,4) = one - rfac * (one + orb2(6))
  rw(4,6) = - rfac * orb2(4)
  rw(6,6) = one - two * rfac * (one + orb2(6))
  call m66mpy(rw, re, re)
500 continue

end subroutine emdamp
