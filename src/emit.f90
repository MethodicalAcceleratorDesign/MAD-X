subroutine emit(deltap, tol, orbit0, disp0, rt, u0, emit_v, nemit_v, &
                bmax, gmax, dismax, tunes, sig_v, pdamp, updatebeam)
  use bbfi
  use twiss0fi
  use emitfi
  use matrices, only : EYE
  use math_constfi, only : zero, one, three, twopi
  use code_constfi
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
  !   nemit_v    (real)   exn, eyn, etn (normalised emittances)
  !   bmax       (real)   Maximum extents of modes.
  !   gmax       (real)   Maximum divergences of modes.
  !   dismax     (real)   Maximum dispersion.
  !   tunes      (real)   qx, qy, qs
  !   pdamp      (real)   damping partition numbers
  !   sig_v      (real)   sigx, sigy, sigt, sige
  !   updatebeam (logical) flag to trigger BEAM update upon return
  !----------------------------------------------------------------------*
  !---- Communication area for radiation damping.

  double precision, intent(IN) :: deltap, tol
  double precision, intent(IN) :: orbit0(6), disp0(6)
  double precision, intent(IN OUT) :: rt(6,6)
  double precision, intent(OUT) :: u0, emit_v(3), nemit_v(3)
  double precision, intent(OUT) :: bmax(3,3), gmax(3,3), dismax(4)
  double precision, intent(OUT) :: tunes(3), sig_v(4), pdamp(3)
  logical, intent(OUT) :: updatebeam

  double precision :: orbit(6), orbit1(6), orbit2(6)
  double precision :: em(6,6), em2(6,6), rd(6,6) ! eigenvalues and damping matrices
  double precision :: reval(6), aival(6), ek(6) ! re and im parts for damping and tunes
  double precision :: disp(6)
  double precision :: al_errors(align_max)
  double precision :: arad, suml, gammas, el
  double precision :: re(6,6), te(6,6,6), tt(6,6,6)
  double precision :: betas, bx, gx

  integer :: i, j, j1, j2, k, k1, k2, eflag, n_align, code
  logical :: fmap, stabx, staby, stabt, radiate

  integer, external :: restart_sequ, advance_node, node_al_errors
  logical, external :: m66sta
  double precision, external :: get_value, node_value

  ORBIT(:6) = ORBIT0(:6)
  DISP(:6) = DISP0(:6)

  EMIT_V(:3) = zero
  NEMIT_V(:3) = zero
  BMAX(:3,:3)  = zero
  GMAX(:3,:3)  = zero
  DISMAX(:4) = zero
  TUNES(:3) = zero
  SIG_V(:4) = zero
  PDAMP(:3) = zero
  u0 = zero

  !---- Find eigenvectors at initial position.
  if (m66sta(rt)) then
     call laseig(rt, reval, aival, em)
     stabt = .false.
  else
     call ladeig(rt, reval, aival, em)
     stabt = reval(5)**2 + aival(5)**2 .le. tol  .and.  &
             reval(6)**2 + aival(6)**2 .le. tol
  endif
  stabx = reval(1)**2 + aival(1)**2 .le. tol  .and.     &
          reval(2)**2 + aival(2)**2 .le. tol
  staby = reval(3)**2 + aival(3)**2 .le. tol  .and.     &
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
  radiate = get_value('probe ','radiate ') .ne. zero

  !---- Initialize damping calculation.
  if (radiate .and. stabt) then
     SUM(:3) = zero
     sumu0 = zero
     RD = EYE
  endif

  TT = zero
  RT = EYE  ! redefining RT

  eflag = 0
  suml = zero
  bbd_cnt=0
  bbd_flag=1

  i = restart_sequ()

10 continue
   bbd_pos=i
   el = node_value('l ')

   code = node_value('mad8_type ')
   if(code .eq. code_tkicker)     code = code_kicker
   if(code .eq. code_placeholder) code = code_instrument

   n_align = node_al_errors(al_errors)
   if (n_align .ne. 0)  then
      ORBIT2 = ORBIT
      call tmali1(orbit2,al_errors,betas,gammas,orbit,re)
      if (.not. stabt) DISP = matmul(RE,DISP)
      EM = matmul(RE,EM)
      if (radiate .and. stabt) RD = matmul(RE,RD)
   endif

  !---- Keep orbit at entrance.
  ORBIT1 = ORBIT

  !---- Element matrix and length.
  call tmmap(code,.true.,.true.,orbit,fmap,ek,re,te,.false.,el)

  if (fmap) then
     !---- Advance dispersion.
     if (.not. stabt) then
        DISP = matmul(RE,DISP)
        do j = 1, 4
           dismax(j) = max(abs(disp(j)),dismax(j))
        enddo
     endif

     !---- Radiation damping.
     EM2 = matmul(RE,EM)
     if (radiate .and. stabt) then
        call emdamp(code, deltap, em, em2, orbit1, orbit, re)
        RD = matmul(RE,RD)
     endif
     EM = EM2

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
     ORBIT2 = ORBIT
     call tmali2(el,orbit2,al_errors,betas,gammas,orbit,re)
     if (.not. stabt) DISP = matmul(RE,DISP)
     EM = matmul(RE,EM)
     if (radiate .and. stabt) RD = matmul(RE,RD)
  endif

  if (advance_node().ne.0)  then
     i=i+1
     goto 10
  endif

  bbd_flag=0

  !---- Undamped tunes and beam extents.
  qx = atan2(aival(1), reval(1)) / twopi ; if (qx .lt. zero) qx = qx + one
  qy = atan2(aival(3), reval(3)) / twopi ; if (qy .lt. zero) qy = qy + one
  qs = atan2(aival(5), reval(5)) / twopi ; if (qs .lt. zero) qs = - qs

  !---- Summary output.
  call emsumm(rd,em,bmax,gmax,stabt,radiate,u0,emit_v,nemit_v,tunes,sig_v,pdamp)

  updatebeam = radiate .and. stabt

end subroutine emit

subroutine emdamp(code, deltap, em1, em2, orb1, orb2, re)
  use twiss0fi
  use emitfi
  use twtrrfi
  use matrices
  use math_constfi, only : zero, one, two, three, four, six, twelve, pi, half, ten3m, ten6p
  use phys_constfi, only : clight
  use code_constfi
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
  integer :: code
  double precision :: deltap
  double precision :: em1(6,6), em2(6,6), orb1(6), orb2(6), re(6,6)

  integer :: i, j, ir, ii, n, n_ferr, iord, nn, ns, nd, nord
  double precision  :: rw(6,6), tw(6,6,6), rw0(6,6), ferror(2)

  double precision :: e1(6,6), e2(6,6)
  double precision :: normal(0:maxmul), skew(0:maxmul)
  double precision :: f_errors(0:maxferr)

  double precision :: o1(6), o2(6)
  double precision :: rot(6,6)
  double precision :: x1, y1, t1, px1, py1, pt1
  double precision :: x2, y2, t2, px2, py2, pt2

  double precision :: el, tilt, bvk
  double precision :: edg1, edg2, sk0, sk1, sk2, hgap, fint, sks, sksol, h, ct
  double precision :: corr, hx, hy, hxx, hxy, hyy, h1, hcb1, hcbs1
  double precision :: tedg1, fact1, dfact1_dx, rfac, rfac1, drfac1_dx, drfac1_dy
  double precision :: h2, hcb2, tedg2, fact2, dfact2_dx, rfac2
  double precision :: drfac2_dx, drfac2_dy, bi2gi2, betas, gammas
  double precision :: e5sq1, e5sq2, e5sqs1, e5sqs2, x, y
  double precision :: f1, f2, f1s, f2s, twon, str, st
  double precision :: r1sq, r2sq, fh1, fh2, dr, di, drt
  double precision :: rfv, rff, rfl, time
  double precision :: xkick, ykick, dpx, dpy, an, hyx, hcbs2,hbi
  double precision :: sk3, fh
  double precision :: drfac1_dpx, drfac1_dpy, drfac2_dpx, drfac2_dpy
  double precision :: denominator1, denominator2
  double precision :: bet1_sqr, bet2_sqr, dbet1_sqr_dpt, dbet2_sqr_dpt
  double precision :: p1, p2

  integer, external :: node_fd_errors
  double precision, external  :: node_value, get_value


  if (code .eq. code_multipole .or. code .eq. code_rfmultipole)  then
     !--- thin multipole and thin RF multipole
     el = node_value('lrad ')
  else
     el = node_value('l ')
  endif

  if (el.eq.zero .and. code.ne.code_rfcavity) return !- no damping
  ! RF cavities with zero length still accepted beyond this point

  betas  = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')

  !---- Prepare data.
  bvk = node_value('other_bv ')

  ! Switch based on element code for element specific damping
  select case (code)

     case (code_rbend, code_sbend) !---- DIPOLE
        an = bvk * node_value('angle ') * el/node_value('l ')
        tilt = -node_value('tilt ')
        edg1 = bvk * node_value('e1 ')
        edg2 = bvk * node_value('e2 ')
        sk0  = bvk * node_value('k0 ')
        sk1 = bvk * (node_value('k1 ') + node_value('k1tap ')) 
        sk2 = bvk * (node_value('k2 ') + node_value('k2tap '))
        hgap = node_value('hgap ')
        fint = node_value('fint ')
        sks = zero
        h = an / el

       if (sk0.ne.0)  h = sk0
       
        !---- Refer orbit and eigenvectors to magnet midplane.
        ct = cos(tilt)
        st = sin(tilt)

        ROT = reshape((/  ct, zero,   st, zero, zero, zero, &
                        zero,   ct, zero,   st, zero, zero, &
                         -st, zero,   ct, zero, zero, zero, &
                        zero,  -st, zero,   ct, zero, zero, &
                        zero, zero, zero, zero,  one, zero, &
                        zero, zero, zero, zero, zero,  one /), shape(ROT))

        O1 = matmul(ROT, ORB1)
        O2 = matmul(ROT, ORB2)

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

        !---- Move orbit through fringing field;
        !     Requested components of eigenvectors are not affected.
        corr = (h + h) * hgap * fint
        RW = EYE
        TW = zero
        call tmfrng(.false.,h,sk1,edg1,zero,+one,corr,rw,tw)
        O1 = matmul(RW, O1)

        RW = EYE
        TW = zero
        call tmfrng(.false.,h,sk1,edg2,zero,-one,corr,rw,tw)
        RW0 = RW
        RW = matmul(JMATT, matmul(transpose(RW0),JMAT)) !invert symplectic matrix
        O2 = matmul(RW, O2)

        !--- For better readibility of following equations
        x1 = o1(1); px1 = o1(2); y1 = o1(3); py1 = o1(4); t1 = o1(5); pt1 = o1(6)
        x2 = o2(1); px2 = o2(2); y2 = o2(3); py2 = o2(4); t2 = o2(5); pt2 = o2(6)

        !---- Local curvature and its derivatives,
        !     Coefficients for damping matrix.
        hx = sk1*x1 + sks*y1 + h + half*sk2 * (x1**2 - y1**2)
        hy = sks*x1 - sk1*y1 - sk2*x1*y1
        hxx = sk1 + sk2*x1
        hxy = sks - sk2*y1
        hyx = hxy
        hyy = -hxx
        h1 = sqrt(hx**2 + hy**2)
        hcb1 = h1**3
        hcbs1 = three*h1 * (hx * (hxx*px1 + hxy*py1) + hy * (hxy*px1 + hyy*py1))

        tedg1  = tan(edg1)
        fact1  = (one + h*x1) * (one - tedg1*x1)
        dfact1_dx = h - tedg1 - two*h*tedg1*x1
        rfac1  = cg*el*h1**2*fact1
        drfac1_dx = cg*el * (two*(hx*hxx+hy*hyx)*fact1 + h1**2*dfact1_dx)
        drfac1_dy = cg*el *  two*(hx*hxy+hy*hyy)*fact1

        hx = sk1*x2 + sks*y2 + h + half*sk2 * (x2**2 - y2**2)
        hy = sks*x2 - sk1*y2 - sk2*x2*y2
        hxx = sk1 + sk2*x2
        hxy = sks - sk2*y2
        hyx = hxy
        hyy = - hxx
        h2 = sqrt(hx**2 + hy**2)
        hcb2 = h2**3
        hcbs2 = three*h2 * (hx * (hxx*px2 + hxy*py2) + hy * (hxy*px2 + hyy*py2))

        tedg2  = tan(edg2)
        fact2  = (one + h*x2) * (one - tedg2*x2)
        dfact2_dx = h - tedg2 - two*h*tedg2*x2

        rfac2  = cg*el*h2**2*fact2
        drfac2_dx = cg*el * (two*(hx*hxx+hy*hyx)*fact2 + h2**2*dfact2_dx)
        drfac2_dy = cg*el *  two*(hx*hxy+hy*hyy)*fact2

        p1 = sqrt(pt1*pt1 + two*pt1/betas + one) ! delta + 1
        p2 = sqrt(pt2*pt2 + two*pt2/betas + one) ! delta + 1
        bet1_sqr = (pt1*pt1 + two*pt1/betas + one) / (one/betas + pt1)**2;
        bet2_sqr = (pt2*pt2 + two*pt2/betas + one) / (one/betas + pt2)**2;
        dbet1_sqr_dpt = (two/betas+two*pt1)/(pt1+one/betas)**2-(two*((pt1*two)/betas+pt1**2+one))/(pt1+one/betas)**3;
        dbet2_sqr_dpt = (two/betas+two*pt2)/(pt2+one/betas)**2-(two*((pt2*two)/betas+pt2**2+one))/(pt2+one/betas)**3;
        denominator1 = two*sqrt(((rfac1-two)*rfac1)/bet1_sqr+one);
        denominator2 = two*sqrt(((rfac2-two)*rfac2)/bet2_sqr+one);
        
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
           sum(i) = sum(i) + half*el * (f1 + f2) - el**2 * (f2s - f1s) / twelve
        enddo

        !---- Damping matrices.
        !     Code common to bending magnet and pure multipoles.
        RW = EYE
        rw(2,1) = (two*drfac1_dx*px1*rfac1-two*drfac1_dx*px1)/(bet1_sqr*denominator1);
        rw(2,2) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpx*px1*rfac1)-(two*drfac1_dpx*px1))/(bet1_sqr*denominator1);
        rw(2,3) = (two*drfac1_dy*px1*rfac1-two*drfac1_dy*px1)/(bet1_sqr*denominator1);
        rw(2,4) = (two*drfac1_dpy*px1*rfac1-two*drfac1_dpy*px1)/(bet1_sqr*denominator1);
        rw(2,6) = -(dbet1_sqr_dpt*px1*rfac1**2-two*dbet1_sqr_dpt*px1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(4,1) = (two*drfac1_dx*py1*rfac1-two*drfac1_dx*py1)/(bet1_sqr*denominator1);
        rw(4,2) = (two*drfac1_dpx*py1*rfac1-two*drfac1_dpx*py1)/(bet1_sqr*denominator1);
        rw(4,3) = (two*drfac1_dy*py1*rfac1-two*drfac1_dy*py1)/(bet1_sqr*denominator1);
        rw(4,4) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpy*py1*rfac1)-(two*drfac1_dpy*py1))/(bet1_sqr*denominator1);
        rw(4,6) = -(dbet1_sqr_dpt*py1*rfac1**2-two*dbet1_sqr_dpt*py1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(6,1) =     - drfac1_dx * p1**2
        rw(6,3) =     - drfac1_dy * p1**2
        rw(6,6) = one - two * rfac1 * (one + pt1)
        !! if the sixth variable was pt, the RW(6,*) equations would look like
        !rw(6,1) = -drfac1_dx*pt1-drfac1_dx/betas;
        !rw(6,2) = -drfac1_dpx*pt1-drfac1_dpx/betas;
        !rw(6,3) = -drfac1_dy*pt1-drfac1_dy/betas;
        !rw(6,4) = -drfac1_dpy*pt1-drfac1_dpy/betas;
        !rw(6,6) = one-rfac1;
        RE = matmul(RE,RW)

        RW = EYE
        rw(2,1) = (two*drfac2_dx*px2*rfac2-two*drfac2_dx*px2)/(bet2_sqr*denominator2);
        rw(2,2) = sqrt(rfac2**2/bet2_sqr-two*rfac2/bet2_sqr+one)+&
             ((two*drfac2_dpx*px2*rfac2)-(two*drfac2_dpx*px2))/(bet2_sqr*denominator2);
        rw(2,3) = (two*drfac2_dy*px2*rfac2-two*drfac2_dy*px2)/(bet2_sqr*denominator2);
        rw(2,4) = (two*drfac2_dpy*px2*rfac2-two*drfac2_dpy*px2)/(bet2_sqr*denominator2);
        rw(2,6) = -(dbet2_sqr_dpt*px2*rfac2**2-two*dbet2_sqr_dpt*px2*rfac2) / &
             (two*bet2_sqr**2*sqrt((rfac2**2-two*rfac2+bet2_sqr)/bet2_sqr));
        rw(4,1) = (two*drfac2_dx*py2*rfac2-two*drfac2_dx*py2)/(bet2_sqr*denominator2);
        rw(4,2) = (two*drfac2_dpx*py2*rfac2-two*drfac2_dpx*py2)/(bet2_sqr*denominator2);
        rw(4,3) = (two*drfac2_dy*py2*rfac2-two*drfac2_dy*py2)/(bet2_sqr*denominator2);
        rw(4,4) = sqrt(rfac2**2/bet2_sqr-two*rfac2/bet2_sqr+one)+&
             ((two*drfac2_dpy*py2*rfac2)-(two*drfac2_dpy*py2))/(bet2_sqr*denominator2);
        rw(4,6) = -(dbet2_sqr_dpt*py2*rfac2**2-two*dbet2_sqr_dpt*py2*rfac2) / &
             (two*bet2_sqr**2*sqrt((rfac2**2-two*rfac2+bet2_sqr)/bet2_sqr));
        rw(6,1) =     - drfac2_dx * p2**2
        rw(6,3) =     - drfac2_dy * p2**2
        rw(6,6) = one - two * rfac2 * (one + pt2)
        !! if the sixth variable was pt, the RW(6,:) equations would look like
        !rw(6,1) = -drfac2_dx*pt2-drfac2_dx/betas;
        !rw(6,2) = -drfac2_dpx*pt2-drfac2_dpx/betas;
        !rw(6,3) = -drfac2_dy*pt2-drfac2_dy/betas;
        !rw(6,4) = -drfac2_dpy*pt2-drfac2_dpy/betas;
        !rw(6,6) = one-rfac2;
        RE = matmul(RW,RE)
        
     case (code_quadrupole , code_sextupole, code_octupole, code_solenoid) !---- Common to all pure multipoles.
        sk1 = zero
        sk2 = zero
        sk3 = zero
        sksol = zero
        select case (code)
        case (code_quadrupole)  !---- Quadrupole
           sk1 = bvk * (node_value('k1 ') + node_value('k1tap ')) 
           str  = sk1
           n    = 1
           twon = two
        case (code_sextupole)   !---- Sextupole
           sk2 = bvk * (node_value('k2 ') + node_value('k2tap '))
           str  = sk2 / two
           n    = 2
           twon = four
        case (code_octupole)   !---- Octupole
           sk3 = bvk * node_value('k2 ')
           str  = sk3 / six
           n    = 3
           twon = six
        case (code_solenoid)  !---- Solenoid
           sksol = node_value('ks ');
           str   = zero
           n     = 0
           twon  = zero
        end select

        O1 = ORB1; O2 = ORB2
        E1 = EM1;  E2 = EM2

        !--- For better readibility of following equations
        x1 = o1(1); px1 = o1(2); y1 = o1(3); py1 = o1(4); t1 = o1(5); pt1 = o1(6)
        x2 = o2(1); px2 = o2(2); y2 = o2(3); py2 = o2(4); t2 = o2(5); pt2 = o2(6)

        !---- Local curvature.
        r1sq = x1**2 + y1**2
        r2sq = x2**2 + y2**2
        h1 = abs(str) * sqrt(r1sq)**n + sksol*(sksol*x1-py1) + sksol*(sksol*y1+px1)
        h2 = abs(str) * sqrt(r2sq)**n + sksol*(sksol*x2-py2) + sksol*(sksol*y2+px2)
        rfac = cg * str**2 * el
        rfac1 = rfac * r1sq**n + cg * sksol*(sksol*x1-py1) / el + cg * sksol*(sksol*y1+px1) / el;
        rfac2 = rfac * r2sq**n + cg * sksol*(sksol*x2-py2) / el + cg * sksol*(sksol*y2+px2) / el;
        drfac1_dx = twon * rfac * r1sq**(n-1) * x1 + (cg*sksol**2)/el
        drfac2_dx = twon * rfac * r1sq**(n-1) * x2 + (cg*sksol**2)/el
        drfac1_dy = twon * rfac * r1sq**(n-1) * y1 + (cg*sksol**2)/el
        drfac2_dy = twon * rfac * r1sq**(n-1) * y2 + (cg*sksol**2)/el
        drfac1_dpx = cg*sksol/el
        drfac2_dpx = cg*sksol/el
        drfac1_dpy = -cg*sksol/el
        drfac2_dpy = -cg*sksol/el

        !
        p1 = sqrt(pt1*pt1 + two*pt1/betas + one) ! delta + 1
        p2 = sqrt(pt2*pt2 + two*pt2/betas + one) ! delta + 1
        bet1_sqr = (pt1*pt1 + two*pt1/betas + one) / (one/betas + pt1)**2;
        bet2_sqr = (pt2*pt2 + two*pt2/betas + one) / (one/betas + pt2)**2;
        denominator1 = 2*sqrt(((rfac1-two)*rfac1)/bet1_sqr+one);
        denominator2 = 2*sqrt(((rfac2-two)*rfac2)/bet2_sqr+one);

        !---- Trapezoidal integration over h**3 * E(k,5) * conjg(E(k,5)).
        fh1 = half * el * h1**3
        fh2 = half * el * h2**3
        sum(1) = sum(1) + fh1 * (e1(5,1)**2 + e1(5,2)**2) + fh2 * (e2(5,1)**2 + e2(5,2)**2)
        sum(2) = sum(2) + fh1 * (e1(5,3)**2 + e1(5,4)**2) + fh2 * (e2(5,3)**2 + e2(5,4)**2)
        sum(3) = sum(3) + fh1 * (e1(5,5)**2 + e1(5,6)**2) + fh2 * (e2(5,5)**2 + e2(5,6)**2)

        !---- Damping matrices.
        !     Code common to bending magnet and pure multipoles.
        RW = EYE
        rw(2,1) = (two*drfac1_dx*px1*rfac1-two*drfac1_dx*px1)/(bet1_sqr*denominator1);
        rw(2,2) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpx*px1*rfac1)-(two*drfac1_dpx*px1))/(bet1_sqr*denominator1);
        rw(2,3) = (two*drfac1_dy*px1*rfac1-two*drfac1_dy*px1)/(bet1_sqr*denominator1);
        rw(2,4) = (two*drfac1_dpy*px1*rfac1-two*drfac1_dpy*px1)/(bet1_sqr*denominator1);
        rw(2,6) = -(dbet1_sqr_dpt*px1*rfac1**2-two*dbet1_sqr_dpt*px1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(4,1) = (two*drfac1_dx*py1*rfac1-two*drfac1_dx*py1)/(bet1_sqr*denominator1);
        rw(4,2) = (two*drfac1_dpx*py1*rfac1-two*drfac1_dpx*py1)/(bet1_sqr*denominator1);
        rw(4,3) = (two*drfac1_dy*py1*rfac1-two*drfac1_dy*py1)/(bet1_sqr*denominator1);
        rw(4,4) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpy*py1*rfac1)-(two*drfac1_dpy*py1))/(bet1_sqr*denominator1);
        rw(4,6) = -(dbet1_sqr_dpt*py1*rfac1**2-two*dbet1_sqr_dpt*py1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(6,1) =     - drfac1_dx * p1**2
        rw(6,3) =     - drfac1_dy * p1**2
        rw(6,6) = one - two * rfac1 * (one + pt1)
        !! if the sixth variable was pt, the RW(6,:) equations would look like
        !rw(6,1) = -drfac1_dx*pt1-drfac1_dx/betas;
        !rw(6,2) = -drfac1_dpx*pt1-drfac1_dpx/betas;
        !rw(6,3) = -drfac1_dy*pt1-drfac1_dy/betas;
        !rw(6,4) = -drfac1_dpy*pt1-drfac1_dpy/betas;
        !rw(6,6) = one-rfac1;
        RE = matmul(RE,RW)

        RW = EYE
        rw(2,1) = (two*drfac2_dx*px2*rfac2-two*drfac2_dx*px2)/(bet2_sqr*denominator2);
        rw(2,2) = sqrt(rfac2**2/bet2_sqr-two*rfac2/bet2_sqr+one)+&
             ((two*drfac2_dpx*px2*rfac2)-(two*drfac2_dpx*px2))/(bet2_sqr*denominator2);
        rw(2,3) = (two*drfac2_dy*px2*rfac2-two*drfac2_dy*px2)/(bet2_sqr*denominator2);
        rw(2,4) = (two*drfac2_dpy*px2*rfac2-two*drfac2_dpy*px2)/(bet2_sqr*denominator2);
        rw(2,6) = -(dbet2_sqr_dpt*px2*rfac2**2-two*dbet2_sqr_dpt*px2*rfac2) / &
             (two*bet2_sqr**2*sqrt((rfac2**2-two*rfac2+bet2_sqr)/bet2_sqr));
        rw(4,1) = (two*drfac2_dx*py2*rfac2-two*drfac2_dx*py2)/(bet2_sqr*denominator2);
        rw(4,2) = (two*drfac2_dpx*py2*rfac2-two*drfac2_dpx*py2)/(bet2_sqr*denominator2);
        rw(4,3) = (two*drfac2_dy*py2*rfac2-two*drfac2_dy*py2)/(bet2_sqr*denominator2);
        rw(4,4) = sqrt(rfac2**2/bet2_sqr-two*rfac2/bet2_sqr+one)+&
             ((two*drfac2_dpy*py2*rfac2)-(two*drfac2_dpy*py2))/(bet2_sqr*denominator2);
        rw(4,6) = -(dbet2_sqr_dpt*py2*rfac2**2-two*dbet2_sqr_dpt*py2*rfac2) / &
             (two*bet2_sqr**2*sqrt((rfac2**2-two*rfac2+bet2_sqr)/bet2_sqr));
        rw(6,1) =     - drfac2_dx * p2**2
        rw(6,3) =     - drfac2_dy * p2**2
        rw(6,6) = one - two * rfac2 * (one + pt2)
        !! if the sixth variable was pt, the RW(6,:) equations would look like
        !rw(6,1) = -drfac2_dx*pt2-drfac2_dx/betas;
        !rw(6,2) = -drfac2_dpx*pt2-drfac2_dpx/betas;
        !rw(6,3) = -drfac2_dy*pt2-drfac2_dy/betas;
        !rw(6,4) = -drfac2_dpy*pt2-drfac2_dpy/betas;
        !rw(6,6) = one-rfac2;
        RE = matmul(RW,RE)

     case (code_multipole) !---- Thin multipoles
        ! EL is ELRAD, the fictitious length for radiation.
        !---- Multipole components.
        F_ERRORS(0:maxferr) = zero ; n_ferr = node_fd_errors(F_ERRORS)

        NORMAL(0:maxmul) = zero ; call get_node_vector('knl ',nn,normal)
        SKEW(0:maxmul) = zero   ; call get_node_vector('ksl ',ns,skew)

        !---- Angle (bvk applied later)
        an = node_value('angle ')
        if (an .ne. 0) f_errors(0) = f_errors(0) + normal(0) - an

        !---- Other components and errors.
        nord = 0
        do i = 0, max(nn, ns, n_ferr/2-1)
           f_errors(2*i)   = bvk * (normal(i) + f_errors(2*i))   / (one + deltap)
           f_errors(2*i+1) = bvk * (skew(i)   + f_errors(2*i+1)) / (one + deltap)
           ! get the maximum effective order; loop runs over maximum of user given values
           if (f_errors(2*i) .ne. zero .or. f_errors(2*i+1) .ne. zero)  nord = i
        enddo

        !---- Multipole kick.
        dr = zero
        di = zero
        do i = nord, 0, -1
           drt = (dr * orb1(1) - di * orb1(3)) / float(i+1) + f_errors(2*i)
           di  = (dr * orb1(3) + di * orb1(1)) / float(i+1) + f_errors(2*i+1)
           dr  = drt
        enddo

        !---- H is local "curvature" due to multipole kick.
        h  = sqrt(dr**2 + di**2) / el
        sum(1) = sum(1) + half*el*h**3 * (em1(5,1)**2 + em1(5,2)**2 + em2(5,1)**2 + em2(5,2)**2)
        sum(2) = sum(2) + half*el*h**3 * (em1(5,3)**2 + em1(5,4)**2 + em2(5,3)**2 + em2(5,4)**2)
        sum(3) = sum(3) + half*el*h**3 * (em1(5,5)**2 + em1(5,6)**2 + em2(5,5)**2 + em2(5,6)**2)

        !---- Damping matrix, is the same at both ends.
        rfac1  = cg * (dr**2 + di**2) / el
        drfac1_dx = cg * (- dr * re(2,1) + di * re(4,1)) / el
        drfac1_dy = cg * (- dr * re(2,3) + di * re(4,3)) / el

        !---- Support variables
        x1 = orb1(1); px1 = orb1(2); y1 = orb1(3); py1 = orb1(4); t1 = orb1(5); pt1 = orb1(6)
        p1 = sqrt(pt1*pt1 + two*pt1/betas + one) ! delta + 1
        bet1_sqr = (pt1*pt1 + two*pt1/betas + one) / (one/betas + pt1)**2;
        dbet1_sqr_dpt = (two/betas+two*pt1)/(pt1+one/betas)**2-(two*((pt1*two)/betas+pt1**2+one))/(pt1+one/betas)**3;
        denominator1 = 2*sqrt(((rfac1-two)*rfac1)/bet1_sqr+one);
        
        RW = EYE
        rw(2,1) = (two*drfac1_dx*px1*rfac1-two*drfac1_dx*px1)/(bet1_sqr*denominator1);
        rw(2,2) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpx*px1*rfac1)-(two*drfac1_dpx*px1))/(bet1_sqr*denominator1);
        rw(2,3) = (two*drfac1_dy*px1*rfac1-two*drfac1_dy*px1)/(bet1_sqr*denominator1);
        rw(2,4) = (two*drfac1_dpy*px1*rfac1-two*drfac1_dpy*px1)/(bet1_sqr*denominator1);
        rw(2,6) = -(dbet1_sqr_dpt*px1*rfac1**2-two*dbet1_sqr_dpt*px1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(4,1) = (two*drfac1_dx*py1*rfac1-two*drfac1_dx*py1)/(bet1_sqr*denominator1);
        rw(4,2) = (two*drfac1_dpx*py1*rfac1-two*drfac1_dpx*py1)/(bet1_sqr*denominator1);
        rw(4,3) = (two*drfac1_dy*py1*rfac1-two*drfac1_dy*py1)/(bet1_sqr*denominator1);
        rw(4,4) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpy*py1*rfac1)-(two*drfac1_dpy*py1))/(bet1_sqr*denominator1);
        rw(4,6) = -(dbet1_sqr_dpt*py1*rfac1**2-two*dbet1_sqr_dpt*py1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(6,1) =     - drfac1_dx * p1**2
        rw(6,3) =     - drfac1_dy * p1**2
        rw(6,6) = one - two * rfac1 * (one + pt1)
        !! if the sixth variable was pt, the RW(6,:) equations would look like
        !rw(6,1) = -drfac1_dx*pt1-drfac1_dx/betas;
        !rw(6,2) = -drfac1_dpx*pt1-drfac1_dpx/betas;
        !rw(6,3) = -drfac1_dy*pt1-drfac1_dy/betas;
        !rw(6,4) = -drfac1_dpy*pt1-drfac1_dpy/betas;
        !rw(6,6) = one-rfac1;

        ! RE = RW * RE * RW
        RE = matmul(RW, matmul(RE,RW))


     case (code_rfcavity) !---- RF cavities.
        rfv = node_value('volt ') * ten3m ! MV but u0 is in GeV
        rff = node_value('freq ') * ten6p * two * pi / clight
        rfl = (node_value('lag ') + node_value('lagtap '))  * two * pi
        time = half * (orb1(5) + orb2(5))
        sumu0 = sumu0 + rfv * sin(rfl - rff * time)


     case (code_hkicker, code_kicker, code_vkicker, code_tkicker) !---- Orbit correctors.

        n_ferr = node_fd_errors(f_errors)

        FERROR(1:2) = zero

        if (n_ferr .gt. 0) FERROR(:2) = F_ERRORS(:min(2,n_ferr))

        select case (code)
        case (code_hkicker)
           xkick = bvk * (node_value('kick ') + node_value('chkick ') + ferror(1))
           ykick = zero
        case (code_kicker, code_tkicker)
           xkick = bvk * (node_value('hkick ') + node_value('chkick ') + ferror(1))
           ykick = bvk * (node_value('vkick ') + node_value('cvkick ') + ferror(2))
        case (code_vkicker)
           xkick = zero
           ykick = bvk * (node_value('kick ') + node_value('cvkick ') + ferror(2))
        case default
           xkick = zero
           ykick = zero
        end select

        !---- Sum up total kicks.
        dpx = xkick / (one + deltap)
        dpy = ykick / (one + deltap)

        !---- Local curvature.
        hx = abs(dpx) / el
        hy = abs(dpy) / el
        rfac = cg * (hx**2 + hx**2) * el
        drfac1_dx  = zero;
        drfac1_dy  = zero;
        drfac1_dpx = zero;
        drfac1_dpy = zero;
        drfac2_dx  = zero;
        drfac2_dy  = zero;
        drfac2_dpx = zero;
        drfac2_dpy = zero;

        !---- Trapezoidal integration over h**3*E(k,5)*E*(k,5).
        fh = half * el * sqrt(hx**2 + hy**2)**3
        sum(1) = sum(1) + fh * (em1(5,1)**2 + em1(5,2)**2 + em2(5,1)**2 + em2(5,2)**2)
        sum(2) = sum(2) + fh * (em1(5,3)**2 + em1(5,4)**2 + em2(5,3)**2 + em2(5,4)**2)
        sum(3) = sum(3) + fh * (em1(5,5)**2 + em1(5,6)**2 + em2(5,5)**2 + em2(5,6)**2)

        !---- Support variables
        x1 = orb1(1); px1 = orb1(2); y1 = orb1(3); py1 = orb1(4); t1 = orb1(5); pt1 = orb1(6)
        x2 = orb2(1); px2 = orb2(2); y2 = orb2(3); py2 = orb2(4); t2 = orb2(5); pt2 = orb2(6)
        
        p1 = sqrt(pt1*pt1 + two*pt1/betas + one)
        p2 = sqrt(pt2*pt2 + two*pt2/betas + one)
        bet1_sqr = (pt1*pt1 + two*pt1/betas + one) / (one/betas + pt1)**2;
        bet2_sqr = (pt2*pt2 + two*pt2/betas + one) / (one/betas + pt2)**2;
        dbet1_sqr_dpt = (two/betas+two*pt1)/(pt1+one/betas)**2-(two*((pt1*two)/betas+pt1**2+one))/(pt1+one/betas)**3;
        dbet2_sqr_dpt = (two/betas+two*pt2)/(pt2+one/betas)**2-(two*((pt2*two)/betas+pt2**2+one))/(pt2+one/betas)**3;
        denominator1 = 2*sqrt(((rfac1-two)*rfac1)/bet1_sqr+one);
        denominator2 = 2*sqrt(((rfac2-two)*rfac2)/bet2_sqr+one);
        
        !---- Damping matrices.
        RW = EYE
        rw(2,1) = (two*drfac1_dx*px1*rfac1-two*drfac1_dx*px1)/(bet1_sqr*denominator1);
        rw(2,2) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpx*px1*rfac1)-(two*drfac1_dpx*px1))/(bet1_sqr*denominator1);
        rw(2,3) = (two*drfac1_dy*px1*rfac1-two*drfac1_dy*px1)/(bet1_sqr*denominator1);
        rw(2,4) = (two*drfac1_dpy*px1*rfac1-two*drfac1_dpy*px1)/(bet1_sqr*denominator1);
        rw(2,6) = -(dbet1_sqr_dpt*px1*rfac1**2-two*dbet1_sqr_dpt*px1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(4,1) = (two*drfac1_dx*py1*rfac1-two*drfac1_dx*py1)/(bet1_sqr*denominator1);
        rw(4,2) = (two*drfac1_dpx*py1*rfac1-two*drfac1_dpx*py1)/(bet1_sqr*denominator1);
        rw(4,3) = (two*drfac1_dy*py1*rfac1-two*drfac1_dy*py1)/(bet1_sqr*denominator1);
        rw(4,4) = sqrt(rfac1**2/bet1_sqr-two*rfac1/bet1_sqr+one)+&
             ((two*drfac1_dpy*py1*rfac1)-(two*drfac1_dpy*py1))/(bet1_sqr*denominator1);
        rw(4,6) = -(dbet1_sqr_dpt*py1*rfac1**2-two*dbet1_sqr_dpt*py1*rfac1) / &
             (two*bet1_sqr**2*sqrt((rfac1**2-two*rfac1+bet1_sqr)/bet1_sqr));
        rw(6,1) =     - drfac1_dx * p1**2
        rw(6,3) =     - drfac1_dy * p1**2
        rw(6,6) = one - two * rfac1 * (one + pt1)
        !! if the sixth variable was pt, the RW(6,:) equations would look like
        !rw(6,1) = -drfac1_dx*pt1-drfac1_dx/betas;
        !rw(6,2) = -drfac1_dpx*pt1-drfac1_dpx/betas;
        !rw(6,3) = -drfac1_dy*pt1-drfac1_dy/betas;
        !rw(6,4) = -drfac1_dpy*pt1-drfac1_dpy/betas;
        !rw(6,6) = one-rfac1;
        RE = matmul(RE,RW)

        RW = EYE
        rw(2,1) = (two*drfac2_dx*px2*rfac2-two*drfac2_dx*px2)/(bet2_sqr*denominator2);
        rw(2,2) = sqrt(rfac2**2/bet2_sqr-two*rfac2/bet2_sqr+one)+&
             ((two*drfac2_dpx*px2*rfac2)-(two*drfac2_dpx*px2))/(bet2_sqr*denominator2);
        rw(2,3) = (two*drfac2_dy*px2*rfac2-two*drfac2_dy*px2)/(bet2_sqr*denominator2);
        rw(2,4) = (two*drfac2_dpy*px2*rfac2-two*drfac2_dpy*px2)/(bet2_sqr*denominator2);
        rw(2,6) = -(dbet2_sqr_dpt*px2*rfac2**2-two*dbet2_sqr_dpt*px2*rfac2) / &
             (two*bet2_sqr**2*sqrt((rfac2**2-two*rfac2+bet2_sqr)/bet2_sqr));
        rw(4,1) = (two*drfac2_dx*py2*rfac2-two*drfac2_dx*py2)/(bet2_sqr*denominator2);
        rw(4,2) = (two*drfac2_dpx*py2*rfac2-two*drfac2_dpx*py2)/(bet2_sqr*denominator2);
        rw(4,3) = (two*drfac2_dy*py2*rfac2-two*drfac2_dy*py2)/(bet2_sqr*denominator2);
        rw(4,4) = sqrt(rfac2**2/bet2_sqr-two*rfac2/bet2_sqr+one)+&
             ((two*drfac2_dpy*py2*rfac2)-(two*drfac2_dpy*py2))/(bet2_sqr*denominator2);
        rw(4,6) = -(dbet2_sqr_dpt*py2*rfac2**2-two*dbet2_sqr_dpt*py2*rfac2) / &
             (two*bet2_sqr**2*sqrt((rfac2**2-two*rfac2+bet2_sqr)/bet2_sqr));
        rw(6,1) =     - drfac2_dx * p2**2
        rw(6,3) =     - drfac2_dy * p2**2
        rw(6,6) = one - two * rfac2 * (one + pt2)
        !! if the sixth variable was pt, the RW(6,:) equations would look like
        !rw(6,1) = -drfac2_dx*pt2-drfac2_dx/betas;
        !rw(6,2) = -drfac2_dpx*pt2-drfac2_dpx/betas;
        !rw(6,3) = -drfac2_dy*pt2-drfac2_dy/betas;
        !rw(6,4) = -drfac2_dpy*pt2-drfac2_dpy/betas;
        !rw(6,6) = one-rfac2;
        RE = matmul(RW,RE)

     case (code_rfmultipole)  !---- thin RF multipole
        ! should provide a combination of damping as for thin multipole and
        ! energy gain as RF cavity

     case default
        ! nothing

     end select
     
end subroutine emdamp

subroutine emsumm(rd,em,bmax,gmax,stabt,radiate,u0,emit_v,nemit_v, &
                  tunes,sig_v,pdamp)
  use emitfi
  use math_constfi, only : zero, one, two, three, four, twopi
  use phys_constfi, only : clight, hbar
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
  !   nemit_v    (real)   exn, eyn, etn (normalised emittances)
  !   tunes      (real)   qx, qy, qs
  !   pdamp      (real)   damping partition numbers
  !   sig_v      (real)   sigx, sigy, sigt, sige
  !----------------------------------------------------------------------*
  double precision :: rd(6,6), em(6,6), bmax(3,3), gmax(3,3)
  logical :: stabt, radiate
  double precision :: u0
  double precision :: emit_v(3), nemit_v(3), tunes(3), sig_v(4), pdamp(3)

  integer :: j, j1, j2, k, k1, k2
  double precision :: arad, gammas, clg, en0
  double precision :: amass, freq0, betas
  double precision :: ex, ey, et, exn, eyn, sigx, sigy, sige, sigt
  double precision :: reval(6), aival(6), alj(3), tau(3), tune(3)
  double precision :: sigma(6,6), bstar(3,3), gstar(3,3), dummy(6,6)

  double precision, external :: get_value
  integer, parameter :: iqpr2 = 6
  double precision, parameter :: ten3p=1.0d3, tenp6=1.0d6, tenp9=1.0d9

  SIGMA(:6,:6) = zero
  ex=zero;  ey=zero;  et=zero

  arad   = get_value('probe ','arad ')
  betas  = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')
  amass  = get_value('probe ','mass ')
  freq0  = get_value('probe ','freq0 ')

  !---- Synchrotron energy loss [GeV].
  if (stabt .and. radiate) then
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
     PDAMP(:) = ALJ(:) * two * en0 / u0

     !---- Emittances.
     clg = ((55.d0*hbar*clight) / (96.d0*sqrt(three))) * ((arad * gammas**5) / amass)
     ex = clg * sum(1) / alj(1)
     ey = clg * sum(2) / alj(2)
     et = clg * sum(3) / alj(3)

     !---- Damping constants per second and damping times.
     ALJ(:) = abs(ALJ(:) * freq0 * tenp6)
     TAU(:) = one / ALJ(:)
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

  exn = ex * betas * gammas
  eyn = ey * betas * gammas

  sigx = sqrt(abs(sigma(1,1)))
  sigy = sqrt(abs(sigma(3,3)))
  if (sigma(5,5) .gt. zero .or. sigma(6,6) .gt. zero)  then
     sigt = sqrt(abs(sigma(5,5)));     sige = sqrt(abs(sigma(6,6)))
  else
     sigt = zero;     sige = zero
  endif
  tunes(1) = qx;    tunes(2) = qy;    tunes(3) = qs
  emit_v(1) = ex;   emit_v(2) = ey;   emit_v(3) = et
  nemit_v(1) = exn; nemit_v(2) = eyn
  sig_v(1) = sigx;  sig_v(2) = sigy;  sig_v(3) = sigt;  sig_v(4) = sige

  !---- Summary output; header and global parameters.

  if (stabt) then !---- Dynamic case.
     if (radiate) write (iqpr2, 910) ten3p * u0
     write (iqpr2, 920) 1, 2, 3
     write (iqpr2, 930) qx, qy, qs
     if (radiate) write (iqpr2, 940) tune
     write (iqpr2, 950) ((bstar(j,k), j = 1, 3), k = 1, 3), &
                        ((gstar(j,k), j = 1, 3), k = 1, 3), &
                        ((bmax(j,k), j = 1, 3), k = 1, 3),  &
                        ((gmax(j,k), j = 1, 3), k = 1, 3)
     if (radiate) then
        write (iqpr2, 960) pdamp, alj, (tau(j), j = 1, 3), &
                           ex*tenp6, ey*tenp6, et*tenp6
     endif
  else !---- Static case
     write (iqpr2, 920) 1, 2
     write (iqpr2, 930) qx, qy
     write (iqpr2, 970) ((bstar(j,k), j = 1, 2), k = 1, 2), &
                        ((gstar(j,k), j = 1, 2), k = 1, 2), &
                        ((bmax(j,k), j = 1, 2), k = 1, 2),  &
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
  logical :: stabt
  double precision :: ex, ey, et, em(6,6), sigma(6,6)

  integer :: j, k

  do j = 1, 6
     do k = 1, 6
        sigma(j,k) = ex * (em(j,1)*em(k,1) + em(j,2)*em(k,2)) + &
                     ey * (em(j,3)*em(k,3) + em(j,4)*em(k,4))
        if (stabt) &
             sigma(j,k) = sigma(j,k) + et * (em(j,5)*em(k,5) + em(j,6)*em(k,6))
     enddo
  enddo
end subroutine emce2i

subroutine getclor(orbit0, rt, tt, error)
  !----------------------------------------------------------------------*
  ! Purpose:
  !   Get periodic closed orbit (e.g. at start of Twiss),
  !   first + second order one-turn map
  ! Called from mad_emit.c
  ! Input:
  !   orbit0(6)   (real)  initial guess
  ! Output:
  !   rt(6,6)     (real)  one-turn matrix
  !   tt(6,6,6)   (real)  one-turn second-order map
  !   error       (int)   error flag (0: OK, else != 0)
  !----------------------------------------------------------------------*
  use twiss0fi
  use matrices, only : EYE
  use math_constfi, only : zero
  implicit none

  double precision :: orbit0(6), rt(6,6), tt(6,6,6)
  integer :: error

  double precision :: opt(fundim)

  RT  = EYE
  OPT = zero
  call tmclor(orbit0, .true., .true., opt, rt, tt, error)
end subroutine getclor
