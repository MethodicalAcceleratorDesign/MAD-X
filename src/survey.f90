!  Routines for the survey command in MADX / A. Verdier (started October 2001)
subroutine survey
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Execute SURVEY command.                                            *
  ! Attributes, must be given in this order in the dictionary:           *
  !   X0        (real)    Initial X position.                            *
  !   Y0        (real)    Initial Y position.                            *
  !   Z0        (real)    Initial Z position.                            *
  !   THETA0    (real)    Initial azimuthal angle.                       *
  !   PHI0      (real)    Initial elevation angle.                       *
  !   PSI0      (real)    Initial roll angle.                            *
  !----------------------------------------------------------------------*
  ! Modified: 01-APR-1999, M. Woodley (SLAC)                             *
  !   If we're doing tape file output and there are LCAVITY elements in  *
  !   the current beamline, initialize ENER1 (in COMMON /OPTIC1/) using  *
  !   ENERGY from BEAM common, and call TMLCAV for each one to update    *
  !   ENERGY                                                             *
  !----------------------------------------------------------------------*
  integer :: i, j, code, add_pass, passes, n_add_angle
  integer :: angle_count, node_count, node_ref(100)
  double precision :: dphi, dpsi, dtheta, phi, phi0, psi, psi0, theta, theta0
  double precision :: sums, el, suml, tilt, globaltilt
  double precision :: v(3), v0(3), ve(3), w(3,3), w0(3,3), we(3,3), tx(3)
  double precision :: add_angle(10), org_ang(100)

  integer, external :: restart_sequ, advance_node, set_cont_sequence
  double precision, external :: proxim, node_value, get_value

  !---- Retrieve command attributes.
  v0(1) =  get_value('survey ','x0 ')
  v0(2) =  get_value('survey ','y0 ')
  v0(3) =  get_value('survey ','z0 ')
  theta0 = get_value('survey ','theta0 ')
  phi0 =   get_value('survey ','phi0 ')
  psi0 =   get_value('survey ','psi0 ')

  !---- Initialise the angles
  theta = theta0
  phi = phi0
  psi = psi0

  !---- Set up initial V and W.
  call sumtrx(theta0, phi0, psi0, w0)
  V = V0
  W = W0

  suml = zero
  sums = zero

5 continue
  !---- loop over passes
  add_pass = get_value('sequence ','add_pass ')   ! multiple passes allowed
  do passes = 0, add_pass
     j = restart_sequ()
     angle_count = 0
     node_count = 0

10   continue
     !---- loop over elements
     node_count = node_count + 1
     if (passes .gt. 0)  then
        call get_node_vector('add_angle ',n_add_angle,add_angle)
        if (n_add_angle .gt. 0 .and. add_angle(passes) .ne. 0.) then
           if (passes .eq. 1) then
              angle_count = angle_count + 1
              node_ref(angle_count) = node_count
              org_ang(angle_count) = node_value('angle ')
           endif
           call store_node_value('angle ', add_angle(passes))
        endif
     endif

     code = node_value('mad8_type ')
     !if (code.eq.39) code=15 ! 2015-Aug-06  21:50:12  ghislain: not required here
     !if (code.eq.38) code=24
     !**** el is the arc length for all bends  ********
     ! LD: 2018.02.01, rbarc is computed by node_value (if needed)...
     el = node_value('l ')
     call suelem(el, ve, we, tilt)
     suml = suml + el
     !**  Compute the coordinates at each point
     !call sutrak(v, w, ve, we)
     V = V + matmul(W,VE)
     W = matmul(W,WE)
     !**  Compute globaltilt HERE : it's the value at the entrance
     globaltilt = psi + tilt
     !**  Compute the survey angles at each point
     call suangl(w, theta, phi, psi)
     !**  Fill the survey table
     call sufill(suml,v, theta, phi, psi,globaltilt)
     if (advance_node().ne.0)  goto 10
     !---- end of loop over elements  ***********************************
  enddo


  if (add_pass .gt. 0) then
     j = restart_sequ()
     angle_count = 1
     node_count = 0

20   continue
     !---- loop over elements to
     ! restore original angle to node if necessary
     node_count = node_count+1
     if (node_ref(angle_count) .eq. node_count)  then
        call store_node_value('angle ', org_ang(angle_count))
        angle_count = angle_count+1
     endif
     if (advance_node().ne.0)  goto 20

  endif
  if (set_cont_sequence() .ne. 0)  goto 5

  !---- Centre of machine. ! 2016.12.14 ldeniau: vars below are never used.
  TX = V - V0
  dtheta = theta - proxim(theta0, theta)
  dphi = phi - proxim(phi0, phi)
  dpsi = psi - proxim(psi0, psi)
end subroutine survey

subroutine suangl(w, theta, phi, psi)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Given a rotation matrix, compute the survey angles.                *
  ! Input:                                                               *
  !   W(3,3)    (real)    Rotation matrix.                               *
  ! Output:                                                              *
  !   THETA     (real)    Azimuthal angle.                               *
  !   PHI       (real)    Elevation angle.                               *
  !   PSI       (real)    Roll angle.                                    *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: w(3,3)
  double precision, intent(OUT) :: theta, phi, psi

  double precision :: arg
  double precision, external :: proxim

  arg = sqrt(w(2,1)**2 + w(2,2)**2)

  phi = atan2(w(2,3), arg)
  theta = proxim(atan2(w(1,3), w(3,3)), theta)
  psi = proxim(atan2(w(2,1), w(2,2)), psi)

end subroutine suangl

subroutine sumtrx(theta, phi, psi, w)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Given three survey angles, compute rotation matrix.                *
  ! Input:                                                               *
  !   THETA     (real)    Azimuthal angle.                               *
  !   PHI       (real)    Elevation angle.                               *
  !   PSI       (real)    Roll angle.                                    *
  ! Output:                                                              *
  !   W(3,3)    (real)    Rotation matrix.                               *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: theta, phi, psi
  double precision, intent(OUT) :: w(3,3)

  double precision :: cosphi, cospsi, costhe, sinphi, sinpsi, sinthe

  costhe = cos(theta);   cosphi = cos(phi);   cospsi = cos(psi)
  sinthe = sin(theta);   sinphi = sin(phi);   sinpsi = sin(psi)

  w(1,1) = + costhe * cospsi - sinthe * sinphi * sinpsi
  w(1,2) = - costhe * sinpsi - sinthe * sinphi * cospsi
  w(1,3) =                     sinthe * cosphi
  w(2,1) =                              cosphi * sinpsi
  w(2,2) =                              cosphi * cospsi
  w(2,3) =                              sinphi
  w(3,1) = - sinthe * cospsi - costhe * sinphi * sinpsi
  w(3,2) = + sinthe * sinpsi - costhe * sinphi * cospsi
  w(3,3) =                     costhe * cosphi

end subroutine sumtrx

subroutine suelem(el, ve, we, tilt)
  use twtrrfi
  use matrices, only : EYE
  use math_constfi, only : zero, one
  use code_constfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Compute Displacement and rotation for one element.                 *
  ! Output:                                                              *
  !   EL        (real)    Element length along design orbit.             *
  !   VE(3)     (real)    Displacement of exit w.r.t. entry.             *
  !   WE(3,3)   (real)    Rotation of exit w.r.t. entry.                 *
  ! Reference pointer used:                                              *
  !   LCELM     /REFER/   Current element bank.                          *
  ! Local links:                                                         *
  !   LSEQ                Beam lines sequence for a lump.                *
  !----------------------------------------------------------------------*
  ! Modified: 17-FEB-2005, F. Tecker                                     *
  !   Return tilt also for elements other than BENDs and MULTIPOLEs      *
  ! Modified: 28-DEC-1998, T. Raubenheimer (SLAC)                        *
  !   Added LCAVITY element at ISP 27                                    *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: el
  double precision, intent(OUT) :: ve(3), we(3,3), tilt

  integer :: code, nn, ns
  double precision :: angle, cospsi, costhe, sinpsi, sinthe, ds, dx, dy, bv, x_t, y_t, z_t
  double precision :: normal(0:maxmul), skew(0:maxmul)

  double precision, external :: node_value

  !---- Branch on subprocess code.
  angle = zero
  dx = zero
  ds = zero

  VE(1:2) = zero ; ve(3) = el
  WE = EYE(:3,:3)

  code = node_value('mad8_type ')
  bv   = node_value('other_bv ')
  tilt = node_value('tilt ') * bv

  select case (code)

     case (code_rbend, code_sbend) !---- RBEND, SBEND
        angle = node_value('angle ') * bv
        if (abs(angle) .ge. 1d-13) then
           dx = el * (cos(angle)-one)/angle
           ds = el * sin(angle)/angle
        else
           ds = el
        endif

        cospsi = cos(tilt);  sinpsi = sin(tilt)
        costhe = cos(angle); sinthe = sin(angle)

        ve(1) = dx * cospsi
        ve(2) = dx * sinpsi
        ve(3) = ds

        we(1,1) = costhe * cospsi*cospsi + sinpsi*sinpsi
        we(2,1) = (costhe - one) * cospsi * sinpsi
        we(3,1) = sinthe * cospsi
        we(1,2) = we(2,1)
        we(2,2) = costhe * sinpsi*sinpsi + cospsi*cospsi
        we(3,2) =  sinthe * sinpsi
        we(1,3) = - we(3,1)
        we(2,3) = - we(3,2)
        we(3,3) = costhe


     case (code_multipole) !---- MULTIPOLE (thin, no length)
        ! Must stay compatible with SBEND (makethin!), i.e. ignore ks0l
        ! LD 2017.11.20, attempt to add angle attribute precedence,
        angle = node_value('angle ')
        if (angle .eq. 0) then
          normal(0) = 0
          call get_node_vector('knl ', nn, normal)
          angle = normal(0)
        endif

        angle = angle * bv
        cospsi = cos(tilt);  sinpsi = sin(tilt)
        costhe = cos(angle); sinthe = sin(angle)

        ! VE is equal to default because length is zero

        we(1,1) = costhe * cospsi*cospsi + sinpsi*sinpsi
        we(2,1) = (costhe - one) * cospsi * sinpsi
        we(3,1) = sinthe * cospsi
        we(1,2) = we(2,1)
        we(2,2) = costhe * sinpsi*sinpsi + cospsi*cospsi
        we(3,2) = sinthe * sinpsi
        we(1,3) = - we(3,1)
        we(2,3) = - we(3,2)
        we(3,3) = costhe


     case (code_xrotation) !---- Rotation around X-axis.  QUESTIONABLE USEFULNESS  !!!!!!!!!!!!!
        dy = node_value('angle ') * bv
        we(2,2) =  cos(dy)
        we(3,2) =  sin(dy)
        we(2,3) = -sin(dy)
        we(3,3) =  cos(dy)


     case (code_yrotation) !---- Rotation around Y-axis.  QUESTIONABLE USEFULNESS  !!!!!!!!!!!!!
        dx = node_value('angle ') * bv
        we(1,1) =  cos(dx)
        we(3,1) = -sin(dx)
        we(1,3) =  sin(dx)
        we(3,3) =  cos(dx)

     case (code_srotation) !---- Rotation around S-axis. SPECIAL CASE
        tilt = node_value('angle ') * bv
        we(1,1) =  cos(tilt)
        we(2,1) =  sin(tilt) !should be - according to convention in MAD8 PhysG. or MADX manual?
        we(1,2) = -sin(tilt) !should be + according to convention in MAD8 PhysG. or MADX manual?
        we(2,2) =  cos(tilt)


     case(code_translation) !  Translation of the reference system.
        x_t = node_value('dx ')
        y_t = node_value('dy ')
        z_t = node_value('ds ')
        ve(1) =  x_t
        ve(2) =  y_t
        ve(3) =  z_t


     case default
       ! all straight elements and catch all; use default VE and WE

     end select

end subroutine suelem

subroutine sufill(suml, v, theta, phi, psi, globaltilt)
  use twtrrfi
  use math_constfi, only : zero
  use code_constfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   writes the survey output in the table "survey"                     *
  ! Output:                                                              *
  !   EL       (real)    Element length along design orbit.              *
  !   V(3)     (real)    Coordinate at the end of the element            *
  !   theta, phi, psi(real) : the survey angles                            *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: suml, v(3), theta, phi, psi, globaltilt

  integer :: code, nn, ns, i
  double precision :: ang, el, tmp, surv_vect(7)
  double precision :: normal(0:maxmul), skew(0:maxmul)

  double precision, external :: node_value

  el = node_value('l ')
  call string_to_table_curr('survey ', 'name ', 'name ')
  call string_to_table_curr('survey ', 'keyword ', 'base_name ')
  call string_to_table_curr('survey ', 'comments ', 'comments ')
  call double_to_table_curr('survey ', 's ',suml )
  call double_to_table_curr('survey ', 'l ',el )
  call double_to_table_curr('survey ', 'x ',v(1) )
  call double_to_table_curr('survey ', 'y ',v(2) )
  call double_to_table_curr('survey ', 'z ',v(3) )
  call double_to_table_curr('survey ', 'theta ',theta)
  call double_to_table_curr('survey ', 'phi ',phi)
  call double_to_table_curr('survey ', 'psi ',psi)
  call double_to_table_curr('survey ', 'globaltilt ',globaltilt)

  i = node_value('pass_flag ')
  if (i .eq. 0) then
     surv_vect(1) = v(1)
     surv_vect(2) = v(2)
     surv_vect(3) = v(3)
     surv_vect(4) = theta
     surv_vect(5) = phi
     surv_vect(6) = psi
     surv_vect(7) = suml
     tmp = 1
     call store_node_value('pass_flag ', tmp)
     i = 7
     call store_node_vector('surv_data ', i, surv_vect)
  endif

  code = node_value('mad8_type ')

  ang = zero
  if (code .eq. code_rbend .or. code .eq. code_sbend) then ! RBEND or SBEND
     ang = node_value('angle ') * node_value('other_bv ')
  else if (code .eq. code_multipole) then ! multipoles (LD 2014.10.15)
     normal(0) = zero; skew(0) = zero
     call get_node_vector('knl ',nn,normal)
     call get_node_vector('ksl ',ns,skew) ! process ks0l (LD 2014.10.15)
     if (nn.ne.0 .or. ns.ne.0) then
       ang = sqrt(normal(0)**2+skew(0)**2) * node_value('other_bv ')
     endif
  endif

  call double_to_table_curr('survey ', 'angle ',ang)

  tmp = node_value('slot_id ')
  call double_to_table_curr('survey ', 'slot_id ',tmp)

  tmp = node_value('assembly_id ')
  call double_to_table_curr('survey ', 'assembly_id ',tmp)

  tmp = node_value('mech_sep ')
  call double_to_table_curr('survey ', 'mech_sep ',tmp)

  tmp = node_value('v_pos ')
  call double_to_table_curr('survey ', 'v_pos ',tmp)

  tmp = node_value('tilt ')
  call double_to_table_curr('survey ', 'tilt ',tmp)


  call augment_count('survey ')
end subroutine sufill


subroutine survtest
  integer j, length, advance_node
  double precision vector(7)
  integer, external :: restart_sequ
  ! test routine for USE with option SURVEY
  call set_sequence('combseq ')
  print *, "# sequence combseq"
  j = restart_sequ()
10 continue
  call get_node_vector('surv_data ', length, vector)
  print *, vector
  if (advance_node() .ne. 0)  goto 10
  call set_sequence('new3 ')
  print *, "# sequence newcont"
  j = restart_sequ()
20 continue
  call get_node_vector('surv_data ', length, vector)
  print *, vector
  if (advance_node() .ne. 0)  goto 20
end subroutine survtest
!-----------------  end of surtest subroutine --------------------------
