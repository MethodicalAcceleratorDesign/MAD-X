!  Routines for the survey command in MADX / A. Verdier (started October 2001)
subroutine survey



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

  integer i,j,code,restart_sequ,advance_node,add_pass,passes,n_add_angle
  integer angle_count, node_count, node_ref(100),set_cont_sequence
  double precision dphi,dpsi,dtheta,phi,phi0,proxim,psi,psi0,sums,  &
       theta,theta0,v(3),v0(3),ve(3),w(3,3),w0(3,3),we(3,3),tx(3),       &
       node_value,el,suml,get_value,tilt,globaltilt,zero,add_angle(10), &
       org_ang(100)
  parameter(zero=0d0)

  !---- Retrieve command attributes.
  v0(1)=  get_value('survey ','x0 ')
  v0(2)=  get_value('survey ','y0 ')
  v0(3)=  get_value('survey ','z0 ')
  theta0 = get_value('survey ','theta0 ')
  phi0 =   get_value('survey ','phi0 ')
  psi0 =   get_value('survey ','psi0 ')

  !---- Initialise the angles
  theta = theta0
  phi = phi0
  psi = psi0

  !---- Set up initial V and W.
  suml = zero
  sums = zero
  call sumtrx(theta0, phi0, psi0, w0)
  !      theta = theta0
  !      phi =  phi0
  !      psi =  psi0


  !---- (replaces SUCOPY)
  do j = 1, 3
     v(j) = v0(j)
     do i = 1, 3
        w(i,j) = w0(i,j)
     enddo
  enddo
5 continue

  !---- loop over elements  NO SYMMETRIC SUPERPERIOD ANYMORE!   *******
  !      print *,"suml  length   theta(x)   phi(y)    psi(z)   coord."
  add_pass = get_value('sequence ','add_pass ')
  ! multiple passes allowed
  do passes = 0, add_pass
     j = restart_sequ()
     angle_count = 0
     node_count = 0
10   continue
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
     if(code.eq.39) code=15
     if(code.eq.38) code=24
     !      print *,"code   ", code
     !**** el is the arc length for all bends  ********
     el = node_value('l ')
     call suelem(el, ve, we,tilt)
     !      print *,"el, tilt", el, tilt
     suml = suml + el
     !**  Compute the coordinates at each point
     call sutrak(v, w, ve, we)
     !**  Compute globaltilt HERE : it's the value at the entrance
     globaltilt=psi+tilt
     !**  Compute the survey angles at each point
     call suangl(w, theta, phi, psi)
     !**  Fill the survey table
     call sufill(suml,v, theta, phi, psi,globaltilt)
     if (advance_node().ne.0)  goto 10
     !---- end of loop over elements  ***********************************
  enddo
  ! restore original angle to node if necessary
  if (add_pass .gt. 0) then
     j = restart_sequ()
     angle_count = 1
     node_count = 0
20   continue
     node_count = node_count+1
     if (node_ref(angle_count) .eq. node_count)  then
        call store_node_value('angle ', org_ang(angle_count))
        angle_count = angle_count+1
     endif
     if (advance_node().ne.0)  goto 20
  endif
  if (set_cont_sequence() .ne. 0)  goto 5

  !---- Centre of machine.
  do i = 1, 3
     tx(i) = v(i) - v0(i)
  enddo
  dtheta = theta - proxim(theta0, theta)
  dphi = phi - proxim(phi0, phi)
  dpsi = psi - proxim(psi0, psi)
  !      print *,v, theta, phi, psi
end subroutine survey
!-----------------  end of survey  subroutine -------------------------

!***********************************************************************
!  Subroutines necessary : suangl sutrak suelem
!**********************************************************************

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
  double precision arg,theta,phi,psi,w(3,3),proxim

  arg = sqrt(w(2,1)**2 + w(2,2)**2)
  phi = atan2(w(2,3), arg)
  !      print *,"SUANGL: phi =",phi," arg=",arg,"  w23 =",w(2,3),
  !      "  w22 =",w(2,2),"  w21 =",w(2,1)

  !*****  old procedure commented as incompatiblr with YROT
  !      if (arg .gt. 1.0e-20) then
  theta = proxim(atan2(w(1,3), w(3,3)), theta)
  psi = proxim(atan2(w(2,1), w(2,2)), psi)
  !      else
  !        theta = atan2(w(1,3), w(3,3))
  !        psi = proxim(atan2(-w(1,2), w(1,1))-theta, psi)
  !      endif
end subroutine suangl
!-----------------  end of suangl  subroutine -------------------------
!
!**********************************************************************


subroutine sutrak(v, w, ve, we)
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Update global position.                                            *
  ! Input:                                                               *
  !   V(3)      (real)    Global displacement before element.            *
  !   W(3,3)    (real)    Global rotation matrix before element.         *
  !   VE(3)     (real)    Displacement due to element.                   *
  !   WE(3,3)   (real)    Rotation due to element.                       *
  ! Output:                                                              *
  !   V(3)      (real)    Global displacement after element.             *
  !   W(3,3)    (real)    Global rotation matrix after element.          *
  !----------------------------------------------------------------------*
  integer i
  double precision v(3),ve(3),w(3,3),we(3,3),wt1,wt2,wt3

  do i = 1, 3
     v(i) = v(i) + w(i,1)*ve(1) + w(i,2)*ve(2) + w(i,3)*ve(3)
     wt1 = w(i,1)*we(1,1) + w(i,2)*we(2,1) + w(i,3)*we(3,1)
     wt2 = w(i,1)*we(1,2) + w(i,2)*we(2,2) + w(i,3)*we(3,2)
     wt3 = w(i,1)*we(1,3) + w(i,2)*we(2,3) + w(i,3)*we(3,3)
     w(i,1) = wt1
     w(i,2) = wt2
     w(i,3) = wt3
  enddo
end subroutine sutrak
!-----------------  end of sutrak subroutine --------------------------
!

!**********************************************************************
subroutine sufill(suml,v, theta, phi, psi,globaltilt)


  use twtrrfi
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   writes the survey output in the table "survey"                     *
  ! Output:                                                              *
  !   EL       (real)    Element length along design orbit.              *
  !   V(3)     (real)    Coordinate at the end of the element            *
  ! theta, phi, psi(real) : the survey angles                            *
  !----------------------------------------------------------------------*
  integer code,nn, i
  double precision ang,el,v(3),theta,phi,psi,node_value,suml,       &
       normal(0:maxmul),globaltilt,tmp, surv_vect(7)

  el = node_value('l ')
  call string_to_table('survey ', 'name ', 'name ')
  call double_to_table('survey ', 's ',suml )
  call double_to_table('survey ', 'l ',el )
  call double_to_table('survey ', 'x ',v(1) )
  call double_to_table('survey ', 'y ',v(2) )
  call double_to_table('survey ', 'z ',v(3) )
  call double_to_table('survey ', 'theta ',theta)
  call double_to_table('survey ', 'phi ',phi)
  call double_to_table('survey ', 'psi ',psi)
  call double_to_table('survey ', 'globaltilt ',globaltilt)
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
  if(code.eq.39) code=15
  if(code.eq.38) code=24
  if(code.eq.2.or.code.eq.3) then
     ang = node_value('angle ')*node_value('other_bv ')
  else if(code.eq.8) then
     call get_node_vector('knl ',nn,normal)
     if (nn .ne. 0) then
        ang = normal(0)
     else
        ang = 0d0
     endif
  else
     ang = 0d0
  endif
  call double_to_table('survey ', 'angle ',ang)

  !---------------- --------------------------
  !     copy over the attributes 'mech_sep' and 'assembly_id'
  !     FT 06.06.2008

  tmp = node_value('slot_id ')
  call double_to_table('survey ', 'slot_id ',tmp)

  tmp = node_value('assembly_id ')
  call double_to_table('survey ', 'assembly_id ',tmp)

  tmp = node_value('mech_sep ')
  call double_to_table('survey ', 'mech_sep ',tmp)

  !== jln dealt with the new property v_pos as for mech_sep
  tmp = node_value('v_pos ')
  call double_to_table('survey ', 'v_pos ',tmp)
  !==

  call augment_count('survey ')
end subroutine sufill
!-----------------  end of sufill subroutine --------------------------
subroutine survtest
  integer j, length, advance_node
  double precision vector(7)
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
