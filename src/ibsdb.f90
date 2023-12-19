! **************************************************************************
! Note by F. Antoniou and F. Zimmermann March 2012
! Note that in this version the dispersion is corrected (multiplied by beta) 
! within the module as in the twiss table the disperion is given in the pt
! frame (dx/dpt) while for the ibs calculations the dx/dp is needed.
!
! 12/02/2019 A. Saa Hernandez. Formulas for coasting beam corrected
! 15/12/2023 F. Soubelet - Tried to document the code as much as possible
! ***************************************************************************

! --------------------------------------------------------------------------- !
! Here is a foreword about the logic of the IBS routine.                      !
!                                                                             !
! I will refer to the CERN-ATS-2012-066 note by F. Antoniou and F. Zimmermann !
! (found at https://cds.cern.ch/record/1445924/files/CERN-ATS-2012-066.pdf)   !
! in the comments in this file. When referring to an equation, a table, or to !
! "the note", it will be to this document.                                    !
!                                                                             !
! The IBS routine, essentially, loops through all elements in the lattice and !
! computes, at each element, the horizontal, vertical and longitudinal B&M    !
! integrals as they are defined in Eq (8) of the note. Beware: these          !
! intermediate values include the square bracket terms of Eq (8).             !
! Afterwards, these rates are averaged over the ring and multiplied with the  !
! common fraction of Eq (8), which results in the final values given back by  !
! the routine, and stored in the global variables 'ibs.tx', 'ibs.ty' and      !
! 'ibs.tl' later on.                                                          !
!                                                                             !
! The routine also computes "average" / "weighted" values, but these are only !
! displayed to the user via stdout.                                           !
! --------------------------------------------------------------------------- !

subroutine ibs
   use ibsdbfi
   use name_lenfi
   use math_constfi, only: zero, one, two, half
   use phys_constfi
   implicit none
   !----------------------------------------------------------------------*
   ! Purpose:                                                             *
   !   INTRABEAM SCATTERING, IBS Command                                  *
   !   These routines are a much reduced version of IBS as taken          *
   !   from the program ZAP, written by M. Zisman.                        *
   !   One should refer to the ZAP USERS MANUAL LBL-21270 UC-28.          *
   ! Attribute:                                                           *
   !   TABLE     (name)    Name of Twiss table.                           *
   !----------------------------------------------------------------------*

   ! Declare all variables to be used
   integer :: step, i, j, flag, testtype, range(2), n
   double precision :: alx, alxbar, alxwtd, aly, alybar, alywtd
   double precision :: betax, betay, beteff, beteffy, bxbar, bxinv, bybar, byinv, bywtd
   double precision :: ax1, ax2, ay1, ay2, bx1, bx2, by1, by2
   double precision :: dx, dx1, dx2, dxbar, dxwtd, dpx, dpx1, dpx2, dpxbr, dpxwtd
   double precision :: dy, dy1, dy2, dybar, dywtd, dpy, dpy1, dpy2, dpybr, dpywtd
   double precision :: taul, taux, tauy, tavl, tavlc, tavx, tavxc, tavy, tavyc
   double precision :: tlbar, tlidc, tlwtd, txbar, txidc, txwtd, tybar, tyidc, tywtd
   double precision :: salxb, salyb, sbxb, sbxinv, sbyb, sbyinv, sdpxb, sdxb, sdpyb, sdyb
   double precision :: hscrpt, hscwtd, hscrpty, hscwtdy
   double precision :: s1, s2, ss2, l1, l2, ll2, const, dels, wnorm, sdum, tol

   ! Functions that are provided by MAD-X to query / update the internal state 
   integer, external :: get_option, double_from_table_row, restart_sequ, advance_to_pos
   double precision, external :: get_value

   ! Query parameters from the common blocks
   charge = get_value('probe ', 'charge ')     ! particle charge, from the beam
   gammas = get_value('probe ', 'gamma ')      ! relativistic gamma, from the beam
   gamma  = get_value('probe ', 'gamma ')      ! relativistic gamma, from the beam
   en0    = get_value('probe ', 'energy ')     ! energy per particle, from the beam
   amass  = get_value('probe ', 'mass ')       ! particle rest mass, from the beam
   ex     = get_value('probe ', 'ex ')         ! geometric hor. emittance, from the beam
   ey     = get_value('probe ', 'ey ')         ! geometric ver. emittance, from the beam
   et     = get_value('probe ', 'et ')         ! geometric long. emittance, from the beam
   sigt   = get_value('probe ', 'sigt ')       ! bunch length, from the beam
   sige   = get_value('probe ', 'sige ')       ! relative energy spread, from the beam
   parnum = get_value('probe ', 'npart ')      ! number of particles per bunch, from the beam
   circ   = get_value('probe ', 'circ ')       ! machine circumference
   currnt = get_value('probe ', 'bcurrent ')   ! bunch current, from the beam
   betas  = get_value('probe ', 'beta ')       ! relativistic gamma, from the beam
   beta   = get_value('probe ', 'beta ')       ! relativistic gamma, from the beam
   arad   = get_value('probe ', 'arad ')       ! classical particle radium, from the beam
   alfa   = get_value('probe ', 'alfa ')       ! momentum compaction factor, from the beam
   freq0  = get_value('probe ', 'freq0 ')      ! revolution frequency, from the beam
   bunch  = get_value('probe ', 'kbunch ')     ! number of particles per bunch, from the beam
   step   = get_value('ibs ', 'steps ')        ! IBS step - not used
   tol    = get_value('ibs ', 'tolerance ')    ! IBS tolerance (to check for integral convergence)

   ! Computes dp/p from sigma_e (the dE/E), according to dp/p = (dE/E) / beta**2.
   ! For some reason the variable name is unchanged but it's good to remember it's dp/p now
   sige = sige/beta/beta
   print *, 'sige ', sige

   ! Initialize the following variables to accumulate weighted average lifetimes
   tavlc = zero; tavxc = zero; tavyc = zero                  ! these will hold / accumulate the intermediate B&M integrals
   sbxb = zero; salxb = zero; sdxb = zero; sdpxb = zero      ! will store the "accumulated" values of betx, alfx, dx and dpx
   sbyb = zero; salyb = zero; sdyb = zero; sdpyb = zero      ! will store the "accumulated" values of bety, alfy, dy and dpy
   sbxinv = zero; sbyinv = zero                              ! will store the "accumulated" values of 1/betx and 1/bety
   alxwtd = zero; dxwtd = zero; dpxwtd = zero                ! will store values to compute the "weighted" quantities (see later below)
   bywtd = zero; alywtd = zero; dywtd = zero; dpywtd = zero  ! will store values to compute the "weighted" quantities (see later below)
   hscwtd = zero; hscwtdy = zero                             ! will store values to compute the "weighted" quantities (see later below)
   wnorm = zero                                              ! will store the accumulated length of all elements with non-zero Dx

   ! Now we start a quick TWISS Table reading. The following block reads the value of various optics functions
   ! at the first element (Fortran is 1-indexed). The double_from_table_row returns a non-zero value in case of
   ! issue, which is caught by the check at end of line, and triggers the 'goto 102' statement which will stop
   ! the program and raise and error to the user. This is mostly to make sure TWISS was called before IBS.
   call table_range('twiss ', '#s/#e ', range)
   flag = double_from_table_row('twiss ', 's ', range(1), s1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'l ', range(1), l1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'betx ', range(1), bx1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'bety ', range(1), by1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'alfx ', range(1), ax1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'alfy ', range(1), ay1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'dx ', range(1), dx1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'dpx ', range(1), dpx1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'dy ', range(1), dy1); if (flag .ne. 0) goto 102
   flag = double_from_table_row('twiss ', 'dpy ', range(1), dpy1); if (flag .ne. 0) goto 102

   ! NOTE by F.A & F.Z
   ! ************************************************************************************
   ! Added 16.01.2012 to check if the twiss is taken at the center (testtype=2) or the
   ! exit (testtype=1) of the elements.
   ! If testtype=1 linear interpolation is used to calculate the twiss at the center of
   ! the elements.
   !*************************************************************************************

   ! We go back to the first element in the lattice, and loop until we get to an element of
   ! non-zero length, which we need just below to determine if TWISS was centered or not ****
   j = restart_sequ()
   do i = range(1) + 1, range(2)
      j = advance_to_pos('twiss ', i)
      flag = double_from_table_row('twiss ', 's ', i, ss2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'l ', i, ll2); if (flag .ne. 0) goto 102
      if (ll2 .gt. 0.0001) exit ! break loop
   end do

   ! The following checks if TWISS was computed at the center of elements. We are at the first element of non-zero
   ! length, and we compare s at next element minus s at first element of non-zero length (ss2 - s1) with the length
   ! of the second element. If it matches, then the delta_s corresponds to the length of the element and it means we
   ! are getting values at the exit of the elements: not centered. We set testtype to 1 (otherwise, set to 2).
   if ((ss2 - s1) .eq. ll2) then
      testtype = 1  ! this means twiss is not centered
      print *, 'Twiss was calculated at the exit of the elements.'
      print *, 'Twiss functions at the center of the elements are calculated through linear interpolation'
   else if ((ss2 - s1) .eq. (l1 + ll2)/2) then
      testtype = 2  ! this means twiss is centered
      print *, 'Twiss was calculated at the center of the elements. No interpolation is used'
   end if

   ! Check if "ibs_table" (I forced this to always be true if IBS command is called)
   n = get_option('ibs_table ')

   ! The following is a big loop over all elements in the TWISS table. At each step, we read optics functions
   ! for the given lattice element and use these to compute the growth rates at this element, as well as some
   ! other "average" and "weighted" values. Details are given below as we go.
   j = restart_sequ()  ! back to start of sequence.
   do i = range(1) + 1, range(2)
      ! Move to the next element and gets the optics functions.
      j = advance_to_pos('twiss ', i)
      ! The double_from_table_row returns a non-zero value in case of issue, which is caught by the check at end
      ! of line, and triggers the 'goto 102' statement which will stop the program and raise and error to the user.
      ! We get value in 'twiss' table, for given column, at row 'i', and store in the last provided variable (s2, l2 etc.).
      flag = double_from_table_row('twiss ', 's ', i, s2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'l ', i, l2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'betx ', i, bx2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'bety ', i, by2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'alfx ', i, ax2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'alfy ', i, ay2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'dx ', i, dx2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'dpx ', i, dpx2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'dy ', i, dy2); if (flag .ne. 0) goto 102
      flag = double_from_table_row('twiss ', 'dpy ', i, dpy2); if (flag .ne. 0) goto 102

      !  NOTE by F.A & F.Z
      ! ************************************************************************************
      ! Dispersion and Dispersion prime is multiplied by beta, in order to be in the deltap
      ! and not the pt frame. This correction is necessary for non-relativistic beams
      !*************************************************************************************

      ! For each element, we adjust the dx, dy, dpx and dpy values to be in deltap frame (multiplying by beta_rel).
      ! If TWISS was NOT at center of elements we interpolate linearly (see this block) to have values at center.
      if (testtype .eq. 1) then  ! testtype = 1 means twiss was not centered
         dels = s2 - s1                 ! interpolate delta_s
         sdum = half*(s2 + s1)          ! interpolate sum of s values
         betax = half*(bx2 + bx1)       ! interpolate betax
         betay = half*(by2 + by1)       ! interpolate betay
         alx = half*(ax2 + ax1)         ! interpolate alfx
         aly = half*(ay2 + ay1)         ! interpolate alfy
         dx = beta*half*(dx2 + dx1)     ! interpolate dx and do the adjustment to deltap frame
         dpx = beta*half*(dpx2 + dpx1)  ! interpolate dpx and do the adjustment to deltap frame
         dy = beta*half*(dy2 + dy1)     ! interpolate dy and do the adjustment to deltap frame
         dpy = beta*half*(dpy2 + dpy1)  ! interpolate dpy and do the adjustment to deltap frame

      ! Otherwise, if TWISS was at center of elements, we just query values from the table.
      ! The adjustment is still done to dx, dy, dpx and dpy in order to be in the deltap frame.
      else if (testtype .eq. 2) then
         dels = l2                      ! just take the length of the element
         sdum = s2                      ! just take the s value of the element
         betax = bx2                    ! just take the betax value of the element
         betay = by2                    ! just take the betay value of the element
         alx = ax2                      ! just take the alfx value of the element
         aly = ay2                      ! just take the alfy value of the element
         dx = beta*dx2                  ! do the adjustment to deltap frame
         dpx = beta*dpx2                ! do the adjustment to deltap frame
         dy = beta*dy2                  ! do the adjustment to deltap frame
         dpy = beta*dpy2                ! do the adjustment to deltap frame
      end if

      ! Add the queried values to "accumulated" variables that are used for average ring lifetimes calculations.
      sbxb = sbxb + betax*dels
      sbxinv = sbxinv + dels/betax
      sbyb = sbyb + betay*dels
      sbyinv = sbyinv + dels/betay
      salxb = salxb + alx*dels
      salyb = salyb + aly*dels
      sdxb = sdxb + dx*dels
      sdpxb = sdpxb + dpx*dels
      sdyb = sdyb + dy*dels
      sdpyb = sdpyb + dpy*dels

      ! Calculate weighted average values (?) in regions of non-zero Dx's.
      ! These values are used to calculate "average" and "weighted" ring lifetimes later on.
      if (dx .gt. zero) then
         wnorm = wnorm + dels                                                  ! This is the accumulated length of all elements with non-zero Dx
         dxwtd = dxwtd + dels*dx                                               ! This is the accumulated Dx*length of all elements with non-zero Dx
         dpxwtd = dpxwtd + dels*dpx                                            ! This is the accumulated Dpx*length of all elements with non-zero Dx
         dywtd = dywtd + dels*dy                                               ! This is the accumulated Dy*length of all elements with non-zero Dx
         dpywtd = dpywtd + dels*dpy                                            ! This is the accumulated Dpy*length of all elements with non-zero Dx
         bywtd = bywtd + dels/sqrt(betay)                                      ! This is the accumulated 1/sqrt(Bety)*length of all elements with non-zero Dx
         alxwtd = alxwtd + dels*alx                                            ! This is the accumulated alfx*length of all elements with non-zero Dx
         alywtd = alywtd + dels*aly                                            ! This is the accumulated alfy*length of all elements with non-zero Dx
         hscrpt = betax*dpx**2 + two*alx*dx*dpx + (one + alx**2)*dx**2/betax   ! This is the horizontal action variable
         hscrpty = betay*dpy**2 + two*aly*dy*dpy + (one + aly**2)*dy**2/betay  ! This is the vertical action variable
         hscwtd = hscwtd + dels*sqrt(hscrpt)                                   ! This is horizontal "emittance" that is incremented
         hscwtdy = hscwtdy + dels*sqrt(hscrpty)                                ! This is vertical "emittance" that is incremented
      end if

      ! We compute the B&M integrals at the given element. The computed results are stored into the
      ! txidc, tyidc and tlidc variables for the horizontal, vertical and longitudinal planes, respectively.
      call twsint(betax, betay, alx, aly, dx, dpx, dy, dpy, txidc, tyidc, tlidc)  ! for details, see the function

      ! We accumulate contributions. We will use these to do the averaging at the end of
      ! the loop, when we've gone and computed integrals at all elements in the lattice.
      tavlc = tavlc + tlidc*dels  ! accumulated longitudinal B&M integral
      tavxc = tavxc + txidc*dels  ! accumulated horizontal B&M integral
      tavyc = tavyc + tyidc*dels  ! accumulated vertical B&M integral

      ! Fill "ibs_table" (always): add the values computed above in the table (at current row).
      ! Empirically it seems it only gets created if asked to output to file (IBS, FILE=string;)
      if (n .ne. 0) then
         call string_to_table_curr('ibs ', 'name ', 'name ')
         call double_to_table_curr('ibs ', 's ', sdum)        ! Store s position
         call double_to_table_curr('ibs ', 'dels ', dels)     ! Store length difference between consecutive elements
         call double_to_table_curr('ibs ', 'tli ', tlidc)     ! Store longitudinal B&M integral at element
         call double_to_table_curr('ibs ', 'txi ', txidc)     ! Store horizontal B&M integral at element
         call double_to_table_curr('ibs ', 'tyi ', tyidc)     ! Store vertical B&M integral at element
         call double_to_table_curr('ibs ', 'betx ', betax)    ! Store horizontal beta function
         call double_to_table_curr('ibs ', 'alfx ', alx)      ! Store horizontal alpha function
         call double_to_table_curr('ibs ', 'dx ', dx)         ! Store horizontal dispersion function (adjusted for deltap frame)
         call double_to_table_curr('ibs ', 'dpx ', dpx)       ! Store horizontal dispersion prime function (adjusted for deltap frame)
         call double_to_table_curr('ibs ', 'bety ', betay)    ! Store vertical beta function
         call double_to_table_curr('ibs ', 'alfy ', aly)      ! Store vertical alpha function
         call double_to_table_curr('ibs ', 'dy ', dy)         ! Store vertical dispersion function (adjusted for deltap frame)
         call double_to_table_curr('ibs ', 'dpy ', dpy)       ! Store vertical dispersion prime function (adjusted for deltap frame)
         call augment_count('ibs ')                           ! Move to next row in the table
      endif

      ! *********** Make sure the following lines are not moved by the compiler ******
      ! Not sure what this is for
      s1 = s2
      bx1 = bx2
      by1 = by2
      ax1 = ax2
      ay1 = ay2
      dx1 = dx2
      dpx1 = dpx2
      dy1 = dy2
      dpy1 = dpy2

   end do
   ! Right now we are done looping: we've gone through the whole lattice.

   ! The block below computes the "average" value for each of these quantities across the ring.
   ! Considering we have iterated through the lattice with the s2 variable, its current value is the ring circumference.
   ! The variables divided contain "accumulated" values, so we divide the sum of betx (for instance) by the length of machine.
   bxbar = sbxb/s2
   bybar = sbyb/s2
   alxbar = salxb/s2
   alybar = salyb/s2
   dxbar = sdxb/s2
   dpxbr = sdpxb/s2
   dybar = sdyb/s2
   dpybr = sdpyb/s2
   bxinv = sbxinv/s2
   byinv = sbyinv/s2

   ! The same is done here for "effective" quantities. This time we divide "accumulated" variables by 'wnorm', which is
   ! the effective s (accumulated) for all places that have non-zero Dx. These values are like the "averages" just above,
   ! but only taking in consideration the parts of the lattice where we have horizontal dispersion.
   dxwtd = dxwtd/wnorm
   dpxwtd = dpxwtd/wnorm
   dywtd = dywtd/wnorm
   dpywtd = dpywtd/wnorm
   bywtd = bywtd/wnorm
   bywtd = one/bywtd**2
   alxwtd = alxwtd/wnorm
   alywtd = alywtd/wnorm
   hscwtd = (hscwtd/wnorm)**2
   beteff = dxwtd**2/hscwtd

   ! This either seems to be doing a correction to the effective beta_y based on the 'hscwtdy' variable.
   ! From above, it seems 'hscwtdy' is only non-zero (which forces the correction) if there is dispersion in x and y,
   ! because it is computed (and incremented) in the loop through the lattice from 'hscrpty' which depends on Dy, and
   ! the increment is only done if Dx is non-zero at the current element. Probably this is a correction because it is
   ! divided by the same in the calculation in this if statement block.
   if (hscwtdy .ne. 0.d0) then
      beteffy = dywtd**2/hscwtdy
   else
      beteffy = bywtd
   end if

   ! Compute beam sizes with average betas. This is straight forward: we get the average
   ! sigmas in horizontal and vertical from the average betas computed above.
   sigx = sqrt(ex*bxbar + (dx*sige)**2)
   sigy = sqrt(ey*bybar + (dy*sige)**2)

   ! These 3 calls print a whole bunch of info and values to stdout.
   call enprgl      ! print global information (flags, beam params etc.)
   call enprem      ! print emittances and sigmas
   call cavprt()    ! print some information for each of the RF elements in the sequence

   ! Here we compute the B&M integrals, but from the "average" quantities computed above. This is to later
   ! get the "average" growth rates across the ring. These values are only printed to stdout, not used.
   call twsint(bxbar, bybar, alxbar, alybar, dxbar, dpxbr, dybar, dpybr, txbar, tybar, tlbar)

   ! Here we compute the B&M integrals, but from the "effective" quantities computed above. This is to later
   ! get the "effective" growth rates across the ring. These values are only printed to stdout, not used.
   call twsint(beteff, beteffy, alxwtd, alywtd, dxwtd, dpxwtd, dywtd, dpywtd, txwtd, tywtd, tlwtd)

   ! We compute the constant term of all calculations (the first fraction in Eq (8)) to put into 'const'.
   ! The average betx, bety, dx, dy values are used (plus some global parameters accessed in the function).
   call twclog(bxbar, bybar, dxbar, dybar, const)

   ! The following block writes the (weighted) average values to stdout for the user.
   write (*, '(/a/)') " Ring average values (m) "
   write (*, '(5x,a,1pe13.5,4x,a,1pe13.5,4x,a,1pe12.5,4x,a,1pe12.5)') &
      "betx   = ", bxbar, "bety   = ", bybar, "Dx  = ", dxbar, "Dy  = ", dybar
   write (*, '(5x,a,1pe13.5,4x,a,1pe13.5,4x,a,1pe12.5,4x,a,1pe12.5)') &
      "alfx   = ", alxbar, "alfy   = ", alybar, "Dpx = ", dpxbr, "Dpy = ", dpybr
   write (*, '(5x,a,1pe13.5,4x,a,1pe13.5)') "1/betx = ", bxinv, "1/bety = ", byinv

   ! We now compute the final result: IBS growth rates! We multiply 'tavlc' / 'tavxc' / 'tavyc' by the
   ! constant term from just above, and we divide by s2 (the ring circumference) to get averages.
   tavl = tavlc*const/s2
   tavx = tavxc*const/s2
   tavy = tavyc*const/s2

   ! Compute the lifetimes too, it's simply the inverse of the growth rates.
   taul = one/tavl
   taux = one/tavx
   tauy = one/tavy

   ! And these are set as global variables in the MAD-X environment, to values guarded by the 'ibs' prefix.
   call set_variable('ibs.tx ', taux)
   call set_variable('ibs.ty ', tauy)
   call set_variable('ibs.tl ', taul)

   ! Some more info values are printed to stdout: the (weighted) average growth rates and lifetimes.
   write (*, '(/5x,a)') "(Weighted) average rates (1/sec):"
   write (*, '( 5x,a,1p,es15.6)') "Longitudinal= ", tavl
   write (*, '( 5x,a,   es15.6)') "Horizontal  = ", tavx
   write (*, '( 5x,a,   es15.6)') "Vertical    = ", tavy
   write (*, '(/5x,a)') "(Weighted) average lifetimes (sec):"
   write (*, '( 5x,a,1p,es15.6)') "Longitudinal= ", taul
   write (*, '( 5x,a,   es15.6)') "Horizontal  = ", taux
   write (*, '( 5x,a,   es15.6/)') "Vertical    = ", tauy

   return  ! this ends the 'ibs' routine call.

   ! This 102 is for the 'goto' statements when querying the TWISS table above in the routine.
   ! It is reached if a desired value is not in the table, and raises an error.
   102 call fort_fail('IBS: ', 'table value not found, rest skipped, program stops ')

end subroutine ibs

! This routine prints to stdout some global information (flags, beam params etc.)
subroutine enprgl
   use ibsdbfi
   use math_constfi, only: zero, one
   use code_constfi
   implicit none
   !----------------------------------------------------------------------*
   ! Purpose:                                                             *
   !   Print global data for machine.                                     *
   !----------------------------------------------------------------------*
   logical :: radiate
   double precision :: eta, gamtr, t0

   double precision, external :: get_value

   radiate = get_value('probe ', 'radiate ') .ne. 0

   !---- Global parameters.
   gamtr = zero
   if (alfa .ne. zero) gamtr = sign(one, alfa)*sqrt(one/abs(alfa))

   t0 = one/freq0
   eta = alfa - one/gamma**2

   write (*, '(/,a,/)') " Global parameters for the machine: "

   write (*, '(a,l1,a/)') "radiate = ", radiate, ":"

   write (*, '(t6,a,t16,f14.6,a,t46,a,t56,f14.6,a,t86,a,t96,f14.6,a)') &
      "C", circ, " m", "f0", freq0, " MHz", "T0", t0, " microseconds"
   write (*, '(t6,a,t16,e18.6,t46,a,t56,e18.6,t86,a,t96,f14.6)') &
      "alfa", alfa, "eta", eta, "gamma(tr)", gamtr
   write (*, '(t6,a,t16,f14.6,a,t46,a,t56,i6,t86,a,t96,e18.6,a)') &
      "Bcurrent", currnt, " A/bunch", "Kbunch", bunch, "Npart", parnum, " per bunch"
   write (*, '(t6,a,t16,f14.6,a,t46,a,t56,f14.6,t86,a,t96,f14.6)') &
      "E", en0, " GeV", "gamma", gamma, "beta", beta

end subroutine enprgl

! This routine prints to stdout the emittances and sigmas
subroutine enprem
   use ibsdbfi
   use math_constfi, only: ten6p, ten3p
   implicit none
   !----------------------------------------------------------------------*
   ! Purpose:                                                             *
   !   Print emittances and sigmas.                                       *
   !----------------------------------------------------------------------*

   write (*, '(/a/)') " Emittances:"
   write (*, '(t6,a,t16,e16.6,a,t48,a,t58,f14.6,a)') &
      "Ex", ten6p*ex, " pi*mm*mrad", "sigx", ten3p*sigx, " mm"
   write (*, '(t6,a,t16,e16.6,a,t48,a,t58,f14.6,a)') &
      "Ey", ten6p*ey, " pi*mm*mrad", "sigy", ten3p*sigy, " mm"
   write (*, '(t6,a,t16,e16.6,a,t48,a,t58,f14.6,a,t88,a,t96,f14.6,a/)') &
      "Et", ten6p*et, " pi*mm*mrad", "sigt", ten3p*sigt, " mm", &
      "sigE", ten3p*sige, " 1/1000"

end subroutine enprem

! This routine prints to stdout information for each of the RF elements in the sequence
subroutine cavprt()
   use name_lenfi
   use code_constfi
   implicit none

   integer :: i, lg
   double precision :: el, rfv, rff, rfl, deltap
   character(len=name_len) :: sequ_name, el_name

   integer, external :: get_string, restart_sequ, advance_node
   double precision, external ::  get_value, node_value

   lg = get_string('sequence ', 'name ', sequ_name)
   if (lg .gt. 0) write (*, '("sequence name: ",a/)') sequ_name(:lg)
   i = restart_sequ()

   do
      if (node_value('mad8_type ') .eq. code_rfcavity) then
         lg = get_string('element ', 'name ', el_name)
         el = node_value('l ')
         rfv = node_value('volt ')
         rff = node_value('freq ')
         rfl = node_value('lag ')
         deltap = get_value('probe ', 'deltap ')
         print '(a,5g14.6)', el_name(:lg), el, rfv, rff, rfl, deltap
      end if
      if (advance_node() .eq. 0) exit
   end do

end subroutine cavprt

! This routine computes the full common constant term (first fraction) of Eq (8), including the Coulomb logarithm.
subroutine twclog(bxbar, bybar, dxbar, dybar, const)
   use ibsdbfi
   use math_constfi, only: zero, two, four, eight, pi
   use phys_constfi, only: hbar, clight, qelect
   implicit none
   !----------------------------------------------------------------------*
   ! Purpose:                                                             *
   !   Calculation of Coulomb logarithm (and print)                       *
   !   based on the formulae in AIP physics vade mecum p.264 (1981)       *
   !   Afterwards, the complete constant term of eq. (8) is calculated.   *
   !   and stored into the provided 'const' variable.
   ! Input:                                                               *
   !   BXBAR     (real)    Average horizontal beta.                       *
   !   BYBAR     (real)    Average vertical beta.                         *
   ! Output:                                                              *
   !   CONST     (real)    Constant in eq. (IV.9.1), ZAP user's manual.   *
   !----------------------------------------------------------------------*
   double precision :: bxbar, bybar, dxbar, dybar, const

   logical :: fbch
   double precision :: bgam, cbunch, coulog
   double precision :: debyel, densty, etrans, pnbtot, qion, tempev, vol
   double precision :: rmax, rmin, rmincl, rminqm, sigtcm, sigxcm, sigycm

   double precision, external :: get_value

   double precision, parameter :: ot2 = 1d2, ft8 = 5d8, ot5 = 1d5, ttm3 = 2d-3
   double precision, parameter :: fac1 = 743.4d0, fac2 = 1.44d-7

   ! Determine if the beam is bunched or coasting. This leads to a correction
   ! factor in the coulomb logarithm calculation (see footnote 4, page 2 in the note).
   fbch = get_value('probe ', 'bunched ') .ne. 0

   ! Calculate transverse temperature as 2*P*X' (assume transverse energy is temperature/2)
   qion = abs(charge)
   etrans = ft8*(gammas*en0 - amass)*(ex/bxbar)
   tempev = two*etrans

   ! Calculate beam volume to get density (in cm**-3).
   sigxcm = ot2*sqrt(ex*bxbar + (dxbar*sige)**2)
   sigycm = ot2*sqrt(ey*bybar + (dybar*sige)**2)
   sigtcm = ot2*sigt
   if (fbch) then                        ! if bunched beam
      vol = eight*sqrt(pi**3)*sigxcm*sigycm*sigtcm
      densty = parnum/vol
   else                                  ! if coasting beam
      vol = four*pi*sigxcm*sigycm*ot2*circ
      densty = parnum/vol
   end if

   ! Calculate RMAX as smaller of SIGXCM and DEBYE length.
   debyel = fac1*sqrt(tempev/densty)/qion
   rmax = min(sigxcm, debyel)

   ! Calculate RMIN as larger of classical distance of closest approach
   ! or quantum mechanical diffraction limit from nuclear radius.
   rmincl = fac2*qion**2/tempev
   rminqm = hbar*clight*ot5/(two*sqrt(ttm3*etrans*amass))
   rmin = max(rmincl, rminqm)
   coulog = log(rmax/rmin)
   bgam = betas*gammas
   if (fbch) then                        ! if bunched beam
      const = parnum*coulog*arad**2*clight/ &
              (eight*pi*betas**3*gammas**4*ex*ey*sige*sigt)
      cbunch = qion*parnum*qelect*betas*clight/circ
   else                                  ! if coasting beam
      const = parnum*coulog*arad**2*clight/ &
              (four*sqrt(pi)*betas**3*gammas**4*ex*ey*sige*circ)
   end if

   ! Print Coulomb logarithm and some beam parameters
   write (*, '(/t6,a,1p,e14.6)') "CONST               = ", const
   write (*, '(/5x,a,f14.6,a)') "ENERGY              = ", en0, " GeV"
   write (*, '( 5x,a,f14.6)') "BETA                = ", betas
   write (*, '( 5x,a,f14.3)') "GAMMA               = ", gammas
   write (*, '( 5x,a,f14.3)') "COULOMB LOG         = ", coulog

   ! We set the coulomb log and the full constant as global variables in the MAD-X environment, guarded by the 'ibs' prefix.
   call set_variable('ibs.coulog ', coulog)
   call set_variable('ibs.const ', const)

   ! Print warning here if Coulomb logarithm gave bad results. Usually this
   ! error is due to a starting guess far from the equilibrium value.
   if (coulog .lt. zero) &
      call fort_warn('TWCLOG: ', 'Coulomb logarithm gives invalid result --- check input parameters.')

   write (*, '(/5x,a,1p,e14.6,a)') "X-emittance         = ", ex, " m*rad"
   write (*, '(5x,a,1p,e14.6,a/)') "Y-emittance         = ", ey, " m*rad"

   if (fbch) then                        ! if bunched beam
      write (*, '(5x,a,1p,e14.6)') "Momentum spread     = ", sige
      write (*, '(5x,a,0p,f14.6,a/)') "Bunch length        = ", sigt, " m"
      write (*, '(5x,a,1p,e14.6)') "Particles per bunch = ", parnum
      write (*, '(5x,a,1p,e14.6,a)') "Bunch current       = ", cbunch, " A"
   else                                  ! if coasting beam
      write (*, '(5x,a,1p,e14.6/)') "Momentum spread     = ", sige
      write (*, '(5x,a,0p,f14.6,a)') "Current             = ", currnt, " A"
   end if

end subroutine twclog


! This routine computes the B&M integral with the provided terms. It is called at each element by the main ibs routine.
! BEWARE: the computed values correspond to the integrals MULTIPLIED by the square brackets term in Eq (8) of the note.
! So, remember that the txi, tyi and tli values (later stored in the IBS table) are not the integrals: they are integrals * brackets.
subroutine twsint(betax, betay, alx, aly, dx, dpx, dy, dpy, txi, tyi, tli)
   use ibsdbfi
   use math_constfi, only: zero, one, two, three, four, six, ten
   implicit none
   !----------------------------------------------------------------------*
   ! Purpose:                                                             *
   !   Subroutine uses Simpson's rule integration                         *
   !   to calculate Bjorken/Mtingwa integrals (eqn. 3.4)                  *
   !   Particle Accelerators 13, 115 (1983)                               *
   !                                                                      *
   !   The expressions found in Conte/Martini                             *
   !   Particle Accelerators 17, 1 (1985) contain two false               *
   !   terms in the expression for tau_x which have been corrected        *
   !   in this version of MADX;                                           *
   !   contributions from vertical dispersion were also added;            *
   !   AB Note by Frank Zimmermann be published (2005)                    *
   !                                                                      *
   !   Integrals are broken into decades to optimize speed.               *
   !                                                                      *
   !   For the VAX, values may not exceed 10**33, therefore TSTLOG=33     *
   !   For the IBM, values may not exceed 10**74, therefore TSTLOG=74     *
   !   (PMG, March 1988)                                                  *
   !                                                                      *
   !   The integral is split into MAXDEC decades with NS steps /decade.   *
   !   TEST is used for testing convergence of the integral               *
   ! Input:                                                               *
   !   BETAX     (real)    Horizontal beta.                               *
   !   BETAY     (real)    Vertical beta.                                 *
   !   ALX       (real)    Horizontal alpha.                              *
   !   ALY       (real)    Vertical alpha.                                *
   !   DX        (real)    Horizontal dispersion.                         *
   !   DPX       (real)    Derivative of horizontal dispersion.           *
   !   DY        (real)    Vertical dispersion.                           *
   !   DPY       (real)    Derivative of vertical dispersion.             *
   ! Output:                                                              *
   !   TXI       (real)    Horizontal result (integral * const_x).        *
   !   TYI       (real)    Vertical result (integral * const_y).          *
   !   TLI       (real)    Longitudinal result (integral * const_l).      *
   !----------------------------------------------------------------------*
   double precision, intent(IN) :: betax, betay, alx, aly, dx, dpx, dy, dpy
   double precision, intent(OUT) :: txi, tyi, tli

   integer :: iiz, iloop, n
   integer, parameter :: maxdec = 30, ns = 50
   integer, external :: get_option
   
   double precision :: a, b, am, c1, c2, c3, c1y, c2y, chy, cx, cy, cl, r1
   double precision :: cscale, chklog, cprime, ccy
   double precision :: zintl, zintx, zinty
   double precision :: td1, td2, tl1, tl2, tx1, tx2, ty1, ty2
   double precision :: al(31), bl(30), h, aloop
   double precision :: term, func, polyl, polyx, polyy, suml, sumx, sumy
   double precision :: cof, f, alam, phi, phiy, tmpl, tmpx, tmpy

   double precision, parameter :: onetominus20 = 1d-20
   double precision, parameter :: tstlog = 74d0, power = -two/three, test = 1d-7

   double precision :: coeff(2)
   data coeff/2d0, 4d0/  ! coefficients for the polynomial, an array of [2.0, 4.0]

   ! We compute intermediate terms that are used later on. Most of these are specific recurrent
   ! parts of the terms defined in Table 1 of the note (a, b, c, ax, bx, al, bl, ay and by).
   ! Remember that 'gammas' is the relativistic gamma
   phi = dpx + (alx*dx/betax)                              ! This is Phi_x from Eq (6) of the note
   phiy = dpy + (aly*dy/betay)                             ! This is Phi_y from Eq (6) of the note
   am = one
   c1 = (gammas*dx)**2/(ex*betax)                          ! (gamma^2 Dx^2) / (eps_x betx) intermediate term
   c1y = (gammas*dy)**2/(ey*betay)                         ! (gamma^2 Dy^2) / (eps_y bety) intermediate term
   c3 = betax/ex                                           ! (betx) / (eps_x) intermediate term
   c2 = c3*(gammas*phi)**2                                 ! (betx gamma^2 phi_x^2) / (eps_x) intermediate term
   cx = c1 + c2                                            ! (gamma^2 Dx^2) / (eps_x betx) + (betx gamma^2 phi_x^2) / (eps_x) intermediate term
   cl = am*(gammas/sige)**2                                ! gamma^2 / sig_p^2 intermediate term
   cy = betay/ey                                           ! (bety) / (eps_y) intermediate term
   c2y = cy*(gammas*phiy)**2                               ! (bety gamma^2 phi_y^2) / (eps_y) intermediate term
   chy = c1y + c2y                                         ! (betx gamma^2 phi_x^2) / (eps_x) + (bety gamma^2 phi_y^2) / (eps_y) intermediate term
   r1 = three/cy                                           ! (3 bety) / (eps_y) intermediate term

   ! Now we compute actual terms from Table 1 of the note (a, b, c). We define CPRIME=C*CSCALE to try to keep the value
   ! small enough for the VAX in single precision or IBM in double precision. Test LOG(C) to see if it needs scaling.
   a = cx + cl + chy + c3 + cy                             ! The a term in Table 1 of the note
   b = (c3 + cy)*(c1 + cl + c1y) + cy*c2 + c3*c2y + c3*cy  ! The b term in Table 1 of the note
   cscale = one
   chklog = log10(c3) + log10(cy) + log10(c1 + cl)
   if (chklog .gt. tstlog) cscale = ten**(tstlog - chklog)
   cprime = c3*cy*cscale*(c1 + cl + c1y)                  ! The c term in Table 1 of the note

   ! The three variables below are used to store (accumulate) the intermediate results on integration intervals
   zintl = zero
   zintx = zero
   zinty = zero

   ! Constants for integration loop. To keep the numbers reasonable, the numerator is scaled by 1/CPRIME
   ! and the denominator by 1/CPRIME**2. The extra factor of CPRIME is accounted for after integrating.
   ccy = cprime**power
   td1 = (a - cy)*ccy
   td2 = one/(sqrt(ccy)*cscale*cy)
   tl1 = (two*a - three*cy - three*c3)/cprime                           ! Corresponds to al / c as they are defined in table 1 of the note
   tl2 = (b - three*c3*cy)/cprime                                       ! Corresponds to bl / c as they are defined in table 1 of the note
   ty1 = (-a + three*cy - chy - chy/cy*(c3 - two*gammas**2/sige**2) &   ! Corresponds to ay / c as they are defined in table 1 of the note
          + two*chy*(cx + chy)/cy + six*c2y)/cprime
   ty2 = (b - c1y*(c3 + cy) + chy*(cy + chy) &                          ! Corresponds to by / c as they are defined in table 1 of the note
          + chy*ey*(one/ey + betax/(betay*ex))*gammas**2/sige**2 &
          - chy*betax/ex*four + (one + (betax*ey)/(betay*ex))*cx*chy &
          + (chy**2)*(betax*ey)/(betay*ex) - chy*ey*c2*c3/betay &
          - c2y*(cy + c3 + chy) + three*c3*(two*c2y + c1y)) &
         /cprime - r1/cscale

   ! Here is a bit of trickery, done for the sake of the order of magnitude of the numbers involved. For the calculation
   ! of ax and bx below, the horizontal square bracket term of Eq (8) (gamma^2 * Hx / epsx) is included too.
   tx1 = (two*(a - c3 - cy)*(cx - c3) - cy*cx &                         ! Corresponds to  (gamma^2 * Hx / epsx) * ax / c as they are defined in table 1 of the note
          + c3*(c1y + six*c2 + c2y + two*c3 + cl - cy))/cprime
   tx2 = (c3 + cx)*((b - c1y*(c3 + cy))/cprime) &                       ! Corresponds to  (gamma^2 * Hx / epsx) * bx / c as they are defined in table 1 of the note
         - six/cscale + three*c3*cy*(cl/cprime) &
         + (six*c3*cy*c1y + (betay/ey + betax/ex)*chy*cx + &
            chy*(c3**2 - two*cy*c3) - c2y*cx*(cy + c3) + (two*cy*c3 - c3*c3)*c2y)/cprime

   al(1) = zero

   ! -------------------------------------------------------------------------- !
   ! The loop below corresponds to the calculation of the integrals themselves. !
   ! It is an iterative solving which follows the Simpson integration rule.     !
   ! The Simpson rule essentially goes as:                                      !
   ! 1. The integration integral is split into "decades", aka sub-intervals.    !
   !    Let's say we have n+1 equally spaced points, with n even, the points    !
   !    are: x0=a, x1, x2, ..., xn=b.                                           !
   ! 2. For each consecutive points pair (xi, xi+1), fit a quadratic polynomial !
   !    that passes through the three points (xi, f(xi)), (xi+1, f(xi+1)),      !
   !    and the midpoint (xi + xi+1)/2, f((xi + xi+1)/2). It is of the form     !
   !    P(x) = a(x−xi+1​/2)**2 + b(x−xi+1​/2) + c  (not the a, b, c of the note)  !
   ! 3. Integrate each quadratic polynomial over its interval [xi,xi+1] to get  !
   !    the approximate area under the curve for that interval.                 !
   ! 4. The formula for the approximation using Simpson's Rule is given by:     !
   !    int_a^b f(x) dx = h/3 * [f(x0​) + 4*f(x1​) + 2*f(x2​) + 4*f(x3​) + ...      !
   !    + 2*f(xn−2​) + 4*f(xn−1​) + f(xn​)]                                        !
   !    where h is the width of each decade: h = (b−a)/n​.                       !
   !                                                                            !
   ! The loop below essentially does the first and fourth steps above. Note:    !
   !   -> Powers of 10 are used for the "decades" limits                         !
   ! -------------------------------------------------------------------------- !
   do iloop = 1, maxdec
      bl(iloop) = ten**iloop             ! this is the end of this subinterval
      al(iloop + 1) = bl(iloop)          ! this is the start of this subinterval
      h = (bl(iloop) - al(iloop))/ns     ! the width of each decade in this subinterval (ns is the number of steps done)
      aloop = al(iloop)

      ! Evaluate Simpson's rule summation for one interval. The integrand is calculated in the loop itself.
      ! Here we evaluate what is essentially the first term in the expansion (step 4 in the description above).
      if (abs(cy + aloop) .gt. onetominus20) then
         term = sqrt((cy + aloop)*ccy)* &
                sqrt((aloop*ccy*aloop + td1*aloop + td2) + aloop*c2y*(c3 - cy)*ccy/(cy + aloop))
      else
         term = sqrt((cy + aloop)*ccy)*sqrt((aloop*ccy*aloop + td1*aloop + td2))
      end if
      func = sqrt(aloop)/term**3  ! we start accumulating contributions
      polyl = tl1*aloop + tl2     ! corresponds to ?
      polyx = tx1*aloop + tx2     ! corresponds to ?
      polyy = ty1*aloop + ty2     ! corresponds to ?
      suml = func*polyl           ! corresponds to ?
      sumx = func*polyx           ! corresponds to ?
      sumy = func*polyy           ! corresponds to ?

      ! We do the same for ns steps: these are the other components of step 4 in the description above.
      ! We have a check in the loop to see if the integral has converged (see lower down).
      do iiz = 1, ns
         alam = aloop + iiz*h
         cof = coeff(mod(iiz, 2) + 1)  ! this is the 4 or 2 in the Simpson's rule summation (starts with 4)
         if (abs(cy + alam) .gt. onetominus20) then
            term = sqrt((cy + alam)*ccy)* &
                   sqrt((alam*ccy*alam + td1*alam + td2) + alam*c2y*(c3 - cy)*ccy/(cy + alam))
         else
            term = sqrt((cy + alam)*ccy)*sqrt((alam*ccy*alam + td1*alam + td2))
         end if
         f = sqrt(alam)/term**3     ! same as above for the next term
         polyl = tl1*alam + tl2     ! same as above for the next term
         polyx = tx1*alam + tx2     ! same as above for the next term
         polyy = ty1*alam + ty2     ! same as above for the next term
         suml = suml + cof*f*polyl  ! same as above for the next term
         sumx = sumx + cof*f*polyx  ! same as above for the next term
         sumy = sumy + cof*f*polyy  ! same as above for the next term
      end do

      ! Accumulate contributions from this decade: compute the full summation (step 4 in the description above).
      suml = suml - f*polyl
      sumx = sumx - f*polyx
      sumy = sumy - f*polyy
      tmpl = (suml/three)*h  ! this is the h/3 * [f(x0​) + 4*f(x1​) + 2*f(x2​) + 4*f(x3​) + ... + 2*f(xn−2​) + 4*f(xn−1​) + f(xn​)] term
      tmpx = (sumx/three)*h  ! this is the h/3 * [f(x0​) + 4*f(x1​) + 2*f(x2​) + 4*f(x3​) + ... + 2*f(xn−2​) + 4*f(xn−1​) + f(xn​)] term
      tmpy = (sumy/three)*h  ! this is the h/3 * [f(x0​) + 4*f(x1​) + 2*f(x2​) + 4*f(x3​) + ... + 2*f(xn−2​) + 4*f(xn−1​) + f(xn​)] term
      zintl = zintl + tmpl   ! we increment the "accumulated" result with this decade's contribution
      zintx = zintx + tmpx   ! we increment the "accumulated" result with this decade's contribution
      zinty = zinty + tmpy   ! we increment the "accumulated" result with this decade's contribution

      ! We test to see if the integral has converged: if the contribution of this decade is lower than 1e-7 (for each
      ! plane) we goto 100 which essentially breaks out of the loop and goes directly to the final result calculation.
      if (abs(tmpl/zintl) .lt. test .and. &
          abs(tmpx/zintx) .lt. test .and. &
          abs(tmpy/zinty) .lt. test) goto 100
   end do

   ! We only get here if there was an issue in the integral calculation: namely if there was no convergence (we would have
   ! gone (goto) to the '100' statement below. So this block prints some values to the terminal and then raises an error.
   write (*, *) tmpl, zintl, tmpx, zintx, tmpy, zinty, test
   write (*, '(a,i3,a)') "Bjorken/Mtingwa integrals did not converge in ", maxdec, " decades."
   call fort_fail('TWSINT: ', 'Problem with TWSINT, program stopped ')

   100 continue

   ! We divide results by CPRIME to account for the above scaling or numerator and denominator. For the vertical and longitudinal
   ! terms, we multiply the integral result by their respective square bracket terms from Eq (8). Remember that for the horizontal
   ! term, tx1 and tx2 already include the square bracket bracket so no multiplication is done (or needed) here.
   txi = (zintx/cprime)     ! horizontal result: bracket_x * integral_x
   tli = cl*(zintl/cprime)  ! longitudinal result: bracket_l * integral_l
   tyi = cy*(zinty/cprime)  ! vertical result: bracket_y * integral_y
   ! The result is stored in the provided variables txi, tyi and tli which are provided by the caller.

end subroutine twsint
