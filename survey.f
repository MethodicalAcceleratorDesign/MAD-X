!  Routines for the survey command in MADX / A. Verdier (October 2001)
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
      integer i,j,code,restart_sequ,advance_node
      double precision dphi,dpsi,dtheta,phi,phi0,proxim,psi,psi0,sums,  &
     &theta,theta0,v(3),v0(3),ve(3),w(3,3),w0(3,3),we(3,3),tx(3),       &
     &node_value,el,suml,get_value,costhe,sinthe,cosphi,sinphi,cospsi,  &
     &sinpsi,tilt

!---- Retrieve command attributes.
      v0(1)=  get_value('survey ','x0 ')
      v0(2)=  get_value('survey ','y0 ')
      v0(3)=  get_value('survey ','z0 ')
      theta0 = get_value('survey ','theta0 ')
      phi0 =   get_value('survey ','phi0 ')
      psi0 =   get_value('survey ','psi0 ')

!---- Initialise the angles
      theta = theta0
      phi =  phi0
      psi =  psi0

!---- Set up initial V and W.
      suml = 0
      sums = 0.
      costhe = cos(theta0)
      sinthe = sin(theta0)
      cosphi = cos(phi0)
      sinphi = sin(phi0)
      cospsi = cos(psi0)
      sinpsi = sin(psi0)
      w0(1,1) = + costhe * cospsi - sinthe * sinphi * sinpsi
      w0(1,2) = - costhe * sinpsi - sinthe * sinphi * cospsi
      w0(1,3) =                     sinthe * cosphi
      w0(2,1) =                              cosphi * sinpsi
      w0(2,2) =                              cosphi * cospsi
      w0(2,3) =                              sinphi
      w0(3,1) = - sinthe * cospsi - costhe * sinphi * sinpsi
      w0(3,2) = + sinthe * sinpsi - costhe * sinphi * cospsi
      w0(3,3) =                     costhe * cosphi
!---- (replaces SUCOPY)
      do j = 1, 3
        v(j) = v0(j)
        do i = 1, 3
          w(i,j) = w0(i,j)
        enddo
      enddo

!---- loop over elements  NO SYMMETRIC SUPERPERIOD ANYMORE!   *******
!      print *,"suml  length   theta(x)   phi(y)    psi(z)   coord."
      j = restart_sequ()
 10   continue
      code = node_value('mad8_type ')
!      print *,"code   ", code 
!**** el is the arc length for all bends  ********
      el = node_value('l ')
      call suelem(el, ve, we,tilt)
      suml = suml + el
!**  Compute the coordinates at each point
      call sutrak(v, w, ve, we)
!**  Compute the survey angles at each point
      call suangl(w, theta, phi, psi)
!**  Fill the survey table
      call sufill(suml,v, theta, phi, psi,tilt)
! Test :
!      print *,suml,"    ",el,"      ",theta,"      ", phi,"     ", psi,
!     +v(1),v(2),v(3)
      if (advance_node().ne.0)  goto 10
!---- end of loop over elements  ***********************************

!---- Centre of machine.
      do i = 1, 3
        tx(i) = v(i) - v0(i)
      enddo
      dtheta = theta - proxim(theta0, theta)
      dphi = phi - proxim(phi0, phi)
      dpsi = psi - proxim(psi0, psi)
!      print *,v, theta, phi, psi
      end
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
      double precision arg,theta,phi,psi,w(3,3),proxim,thetaint,psiint

      arg = sqrt(w(2,1)**2 + w(2,2)**2)
      phi = atan2(w(2,3), arg)
      psiint = proxim(atan2(-w(1,2), w(1,1))-theta, psi)
      if (arg .gt. 1.0e-20) then
        thetaint = proxim(atan2(w(1,3), w(3,3)), theta)
        psiint = proxim(atan2(w(2,1), w(2,2)), psi)
      else
        psiint = proxim(atan2(-w(1,2), w(1,1))-theta, psi)
        thetaint=theta
      endif
      theta=thetaint
      psi=psiint
      end
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
      end
!-----------------  end of sutrak subroutine --------------------------
!
!
!**********************************************************************
      subroutine suelem(el, ve, we,tilt)
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
! Modified: 28-DEC-1998, T. Raubenheimer (SLAC)                        *
!   Added LCAVITY element at ISP 27                                    *
!----------------------------------------------------------------------*
      integer code,nn,ns,maxmul
      parameter(maxmul=20)
      double precision angle,cospsi,costhe,ds,dx,sinpsi,sinthe,tilt,    &
     &ve(3),we(3,3),node_value,el,normal(0:maxmul),skew(0:maxmul),angv
     &,one,get_variable
      parameter(one=1d0)
!---- Branch on subprocess code.
      tilt = 0.0
      code = node_value('mad8_type ')
      go to ( 10,  20,  20,  40,  50,  60,  70,  80,  90, 100,          &
     &110, 120, 130, 140, 150, 160, 170, 180, 190, 200,                 &
     &210, 220, 230, 240, 250,  20, 270, 280, 290, 300,                 &
     &310, 310, 310, 310, 310, 310, 310, 310, 310, 310), code

!---- Drift space.
   10 continue
!---- Arbitrary matrix.
   40 continue
!---- Quadrupole.
   50 continue

!---- Sextupole.
   60 continue

!---- Octupole.
   70 continue

!---- Solenoid.
   90 continue

!---- RF cavity.
  100 continue

!---- Electrostatic separator.
  110 continue

!---- Kickers.
  140 continue
  150 continue
  160 continue

!---- Monitors.
  170 continue
  180 continue
  190 continue

!---- Apertures.
  200 continue
  210 continue

!---- Marker.
  250 continue

!---- Beam-beam.
  220 continue

!---- lcavity
  270 continue
!---- Reserved.
  280 continue
  290 continue
  300 continue

!---- Beam instrument.
  240 continue

!---- Lump.
  230 continue

!---- User-defined elements.
  310 continue

!****** end of straight elements ***************
      ve(1) = 0
      ve(2) = 0
      ve(3) = el
      we(1,1) = one
      we(2,1) = 0
      we(3,1) = 0
      we(1,2) = 0
      we(2,2) = one
      we(3,2) = 0
      we(1,3) = 0
      we(2,3) = 0
      we(3,3) = one
      go to 500

!---- multipoles , introduced  17.09.02 / AV
   80 continue
      call dzero(normal,maxmul+1)
      call dzero(skew,maxmul+1)
      call node_vector('knl ',nn,normal)
      call node_vector('ksl ',ns,skew)
!      print *,"mult ",code,"  angle",normal(0),"  skew ",ns
!     *,skew(0)
!-----  dipole_bv introduced to suppress SU in MADX input (AV  7.10.02)
      angle = normal(0)*node_value('dipole_bv ')
      angv = skew(0)
      if(angle.eq.0.0) then
      tilt = 0.0
        if(angv.ne.0.0) then
        tilt = get_variable('twopi ')*0.25
        endif
      else
      tilt = atan2(angv,angle)
      endif
! As el=0, there is no dx and no ds
        dx = 0.0
        ds = 0.0
      go to 490

!---- Any kind of  bend. 
   20 continue
!--------------  dipole_bv introduced to suppress SU (AV  7.10.02)
      angle = node_value('angle ')*node_value('dipole_bv ')
!      print *," BV = ",node_value('dipole_bv ')
      angv = node_value('k0s ')*el*node_value('dipole_bv ')
      if (angle .eq. 0.0) then
        dx = 0.0
        ds = el
        tilt = 0.0
           if(angv.ne.0.0) then
           tilt = get_variable('twopi ')*0.25
           dx = - el * (cos(angv)-one)/angv
           ds =  el * sin(angv)/angv
           angle = angv
           endif
      else
! el corrected 18.09.02 // identical to mad8(sector bend)
        tilt = atan2(angv,angle)
        dx = el * (cos(angle)-one)/angle
        ds = el * sin(angle)/angle
      endif
!      print *," *****  TILT = ",tilt,"   length= ",el
      go to 490

!---- Rotation around S-axis.
  120 continue
!        call ucopy(q(lcelm+meangr), angle, mwflt)
      we(1,1) = cos(angle)
      we(2,1) = sin(angle)
      we(1,2) = - we(2,1)
      we(2,2) = + we(1,1)
      go to 500

!---- Rotation around Y-axis.
  130 continue
      angle=0
      we(1,1) = cos(angle)
      we(3,1) = sin(angle)
      we(1,3) = - we(3,1)
      we(3,3) = + we(1,1)
      go to 500

!---- Common for bends and multipoles: Displacement and rotation matrix.
  490 continue
      cospsi = cos(tilt)
      sinpsi = sin(tilt)
      costhe = cos(angle)
      sinthe = sin(angle)
      ve(1) = dx * cospsi
      ve(2) = dx * sinpsi
      ve(3) = ds
      we(1,1) = costhe * cospsi*cospsi + sinpsi*sinpsi
      we(2,1) = (costhe - one) * cospsi * sinpsi
      we(3,1) = sinthe * cospsi
      we(1,2) = we(2,1)
      we(2,2) = costhe * sinpsi*sinpsi + cospsi*cospsi
      we(3,2) = - sinthe * sinpsi
      we(1,3) = - we(3,1)
      we(2,3) = - we(3,2)
      we(3,3) = costhe
  500 continue
      end
!-----------------  end of suelem subroutine --------------------------

!**********************************************************************
      subroutine sufill(suml,v, theta, phi, psi,tilt)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   writes the survey output in the table "survey"                     *
! Output:                                                              *
!   EL       (real)    Element length along design orbit.              *
!   V(3)     (real)    Coordinate at the end of the element            *
! theta, phi, psi(real) : the survey angles                            *
!----------------------------------------------------------------------*
      integer code,nn
      double precision ang,el,v(3),theta,phi,psi,node_value,suml,
     &normal(20),tilt,globaltilt

      globaltilt=psi+tilt
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

      code = node_value('mad8_type ')
      if(code.eq.2.or.code.eq.3) then
        ang = node_value('angle ')
      else if(code.eq.8) then
      call node_vector('knl ',nn,normal)
        ang = normal(1)
      else
        ang = 0
      endif
      call double_to_table('survey ', 'angle ',ang)

      call augment_count('survey ')
      end
!-----------------  end of sufill subroutine --------------------------

!**********************************************************************
      subroutine sumtrx(the, phi, psi, w)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Given three survey angles, compute rotation matrix.                *
! Input:                                                               *
!   THE       (real)    Azimuthal angle.                               *
!   PHI       (real)    Elevation angle.                               *
!   PSI       (real)    Roll angle.                                    *
! Output:                                                              *
!   W(3,3)    (real)    Rotation matrix.                               *
!----------------------------------------------------------------------*
      double precision cosphi,cospsi,costhe,phi,psi,sinphi,sinpsi,      &
     &sinthe,the,w(3,3)

      costhe = cos(the)
      sinthe = sin(the)
      cosphi = cos(phi)
      sinphi = sin(phi)
      cospsi = cos(psi)
      sinpsi = sin(psi)
      w(1,1) = + costhe * cospsi - sinthe * sinphi * sinpsi
      w(1,2) = - costhe * sinpsi - sinthe * sinphi * cospsi
      w(1,3) =                     sinthe * cosphi
      w(2,1) =                              cosphi * sinpsi
      w(2,2) =                              cosphi * cospsi
      w(2,3) =                              sinphi
      w(3,1) = - sinthe * cospsi - costhe * sinphi * sinpsi
      w(3,2) = + sinthe * sinpsi - costhe * sinphi * cospsi
      w(3,3) =                     costhe * cosphi

      end
!-----------------  end of sumtrx subroutine --------------------------

!**********************************************************************
      subroutine sutran(w, v, we)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Transform rotation W and displacement V from entrance to exit.     *
! Input:                                                               *
!   W(3,3)    (real)    Rotation matrix w.r.t. input system.           *
!   V(3)      (real)    Displacement w.r.t. input system.              *
!   WE(3,3)   (real)    Rotation matrix due to element.                *
! Output:                                                              *
!   W(3,3)    (real)    Rotation matrix w.r.t. output system.          *
!   V(3)      (real)    Displacement w.r.t. output system.             *
!----------------------------------------------------------------------*
      integer i,k
      double precision v(3),vt(3),w(3,3),we(3,3),wt(3,3)

!---- VT := transpose(WE) * V;
!     WT := transpose(WE) * W;
      do i = 1, 3
        vt(i) = we(1,i)*v(1) + we(2,i)*v(2) + we(3,i)*v(3)
        do k = 1, 3
          wt(i,k) = we(1,i)*w(1,k) + we(2,i)*w(2,k) + we(3,i)*w(3,k)
        enddo
      enddo

!---- V := VT       [= transpose(WE) * V];
!     W := WT * WE  [= transpose(WE) * W * WE];
      do i = 1, 3
        v(i) = vt(i)
        do k = 1, 3
          w(i,k) = wt(i,1)*we(1,k) + wt(i,2)*we(2,k) + wt(i,3)*we(3,k)
        enddo
      enddo
      end
!-----------------  end of sutran subroutine --------------------------
