      subroutine trrun(switch,turns,orbit0,rt,part_id,last_turn,        &
     &last_pos,z,dxt,dyt,last_orbit,eigen,coords,e_flag,code_buf,l_buf)

      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!          Interface RUN and DYNAP command to tracking routine         *
!                                                                      *
!-- Input:                                                             *
!   switch  (int)         1: RUN, 2: DYNAP fastune                     *
!   turns   (int)         number of turns to track                     *
!   orbit0  (dp. array)   start of closed orbit                        *
!   rt      (dp. matrix)  one-turn matrix                              *
!-- Output:                                                            *
!   eigen       dp(6,6)   eigenvector matrix                           *
!   coords      dp(6,0:turns,npart) (only switch > 1) particle coords. *
!   e_flag      (int)         (only switch > 1) 0: OK, else part. lost *
!-- Buffers:                                                           *
!   part_id     int(npart)                                             *
!   last_turn   int(npart)                                             *
!   last_pos    dp(npart)                                              *
!   z           dp(6,npart)                                            *
!   last_orbit  dp(6,npart)                                            *
!   code_buf    int(nelem)  local mad-8 code storage                   *
!   l_buf       dp(nelem)   local length storage                       *
!----------------------------------------------------------------------*
      include 'twiss0.fi'
      include 'name_len.fi'
      include 'track.fi'
      include 'bb.fi'
      logical onepass,onetable,last_out,info,aperflag
      integer j,code,restart_sequ,advance_node,node_al_errors,n_align,  &
     &nlm,jmax,j_tot,turn,turns,i,k,get_option,ffile,SWITCH,nint,ndble, &
     &nchar,part_id(*),last_turn(*),char_l,segment, e_flag, nobs,lobs,  &
     &int_arr(1),tot_segm,code_buf(*)
      double precision tmp_d,orbit0(6),orbit(6),el,re(6,6),rt(6,6),     &
     &al_errors(align_max),z(6,*),dxt(*),dyt(*),eigen(6,6),sum,         &
     &node_value,one,                                                   &
     &get_variable,last_pos(*),last_orbit(6,*),maxaper(6),get_value,    &
     &zero,obs_orb(6),coords(6,0:turns,*),l_buf(*)
      parameter(zero=0d0,one=1d0)
      character*12 tol_a, char_a
!hbu
      double precision spos
!hbu
      character*4 vec_names(7)
!hbu
      character*(name_len) el_name
      character*120 msg
      data tol_a,char_a / 'maxaper ', ' ' /
!hbu
      data vec_names / 'x', 'px', 'y', 'py', 't', 'pt','s' /


!---- AK 2006 04 23
!---- This version of trrun.F gets rid of all problems concerning delta_p
!---- by eliminating any delta_p dependence and using full 6D formulae only!!!
!---- Only the parts of the code that deal with radiation effects still use
!---- the quantity delta_p

!      print *,"madX::trrun.F"
!      print *," "
!      print *," AK special version 2006/04/23"
!      print *," ============================="
!      print *," Full 6D formulae internally only."      


      if(fsecarb) then
        write (msg,*) 'Second order terms of arbitrary Matrix not '//   &
     &'allowed for tracking.'
        call aawarn('trrun: ',msg)
        return
      endif
      aperflag = .false.
      e_flag = 0
! flag to avoid double entry of last line
      last_out = .false.
      onepass = get_option('onepass ') .ne. 0
      onetable = get_option('onetable ') .ne. 0
      info = get_option('info ') * get_option('warn ') .gt. 0
      if(onepass) call m66one(rt)
      call trbegn(rt,eigen)
      if (switch .eq. 1)  then
        ffile = get_value('run ', 'ffile ')
      else
        ffile = 1
      endif
! for one_table case
      segment = 0
      tot_segm = turns / ffile + 1
      if (mod(turns, ffile) .ne. 0) tot_segm = tot_segm + 1
      call trinicmd(switch,orbit0,eigen,jmax,z,turns,coords)
!--- jmax may be reduced by particle loss - keep number in j_tot
      j_tot = jmax
!--- get vector of six coordinate maxapers (both RUN and DYNAP)
      call comm_para(tol_a, nint, ndble, nchar, int_arr, maxaper,       &
     &char_a, char_l)
!--- set particle id
      do k=1,jmax
        part_id(k) = k
      enddo
!hbu--- init info for tables initial s position is 0
!hbu initial s position is 0
      spos=0
!hbu start of line, element 0
      nlm=0
!hbu
      el_name='start           '
!--- enter start coordinates in summary table
      do  i = 1,j_tot
        tmp_d = i
        call double_to_table('tracksumm ', 'number ', tmp_d)
        tmp_d = 0
        call double_to_table('tracksumm ', 'turn ', tmp_d)
        do j = 1, 6
          tmp_d = z(j,i) - orbit0(j)
          call double_to_table('tracksumm ', vec_names(j), tmp_d)
        enddo
!hbu add s
        call double_to_table('tracksumm ',vec_names(7),spos)
        call augment_count('tracksumm ')
      enddo
!--- enter first turn, and possibly eigen in tables
      if (switch .eq. 1)  then
        if (onetable)  then
          call track_pteigen(eigen)
!hbu add s, node id and name
          call tt_putone(jmax, 0, tot_segm, segment, part_id,           &
     &z, orbit0,spos,nlm,el_name)
        else
          do i = 1, jmax
!hbu
            call tt_puttab(part_id(i), 0, 1, z(1,i), orbit0,spos)
          enddo
        endif
      endif
!---- Initialize kinematics and orbit
      bet0  = get_value('beam ','beta ')
      betas = get_value('probe ','beta ')
      gammas= get_value('probe ','gamma ')
      bet0i = one / bet0
      beti   = one / betas
      dtbyds = get_value('probe ','dtbyds ')
      deltas = get_variable('track_deltap ')
      arad = get_value('probe ','arad ')
      dorad = get_value('probe ','radiate ') .ne. 0
      dodamp = get_option('damp ') .ne. 0
      dorand = get_option('quantum ') .ne. 0
      call dcopy(orbit0,orbit,6)
!--- loop over turns
      nobs = 0
      do turn = 1, turns
        j = restart_sequ()
        nlm = 0
        sum=zero
!--- loop over nodes
 10     continue
        bbd_pos=j
        if (turn .eq. 1)  then
          code = node_value('mad8_type ')
          el = node_value('l ')
          code_buf(nlm+1) = code
          l_buf(nlm+1) = el
!hbu get current node name
          call element_name(el_name,len(el_name))
        else
          el = l_buf(nlm+1)
          code = code_buf(nlm+1)
        endif
        if ((code .ne. 1 .and. code .ne. 4) .and. el .ne. zero) then
          print*," "
          print*,"code: ",code," el: ",el,"   THICK ELEMENT FOUND"
          sum = node_value('name ')
          print*," "
          print*,"Track dies nicely"
          print*,"Thick lenses will get nowhere"
          print*,"MAKETHIN will save you"
          print*," "
          print*," "
! Better to use stop. If use return, irrelevant tracking output is printed
          stop
!          call aafail('TRRUN',                                          &
!     &    '----element with length found : CONVERT STRUCTURE WITH '//   &
!     &    'MAKETHIN')
        endif
        if (switch .eq. 1) then
          nobs = node_value('obs_point ')
        endif
!--------  Misalignment at beginning of element (from twissfs.f)
        if (code .ne. 1)  then
          call dzero(al_errors, align_max)
          n_align = node_al_errors(al_errors)
          if (n_align .ne. 0)  then
            do i = 1, jmax
              call tmali1(z(1,i),al_errors, betas, gammas,z(1,i), re)
            enddo
          endif
        endif
!-------- Track through element  // suppress dxt 13.12.04
        call ttmap(code,el,z,jmax,dxt,dyt,sum,turn,part_id, last_turn,  &
     &last_pos, last_orbit,aperflag,maxaper,al_errors,onepass)
!--------  Misalignment at end of element (from twissfs.f)
        if (code .ne. 1)  then
          if (n_align .ne. 0)  then
            do i = 1, jmax
              call tmali2(el,z(1,i), al_errors, betas, gammas,          &
     &z(1,i), re)
            enddo
          endif
        endif
        nlm = nlm+1
        if (nobs .gt. 0)  then
          call dzero(obs_orb,6)
          call get_node_vector('obs_orbit ', lobs, obs_orb)
          if (lobs .lt. 6)                                              &
     &call aafail('TRACK', 'obs. point orbit not found')
          if (onetable)  then
!hbu
            spos=sum
!hbu get current node name
            call element_name(el_name,len(el_name))
!hbu
            call tt_putone(jmax, turn, tot_segm, segment, part_id,      &
     &z, obs_orb,spos,nlm,el_name)
          else
            do i = 1, jmax
!hbu add spos
              call tt_puttab(part_id(i), turn, nobs, z(1,i), obs_orb,   &
     &spos)
            enddo
          endif
        endif
        if (advance_node().ne.0)  then
          j=j+1
          go to 10
        endif
!--- end of loop over nodes
        if (switch .eq. 1)  then
          if (mod(turn, ffile) .eq. 0)  then
            if (turn .eq. turns)  last_out = .true.
            if (onetable)  then
!hbu
              spos=sum
!hbu spos added
              call tt_putone(jmax, turn, tot_segm, segment, part_id,    &
     &z, orbit0,spos,nlm,el_name)
            else
              do i = 1, jmax
!hbu
                call tt_puttab(part_id(i), turn, 1, z(1,i), orbit0,spos)
              enddo
            endif
          endif
        else
          do i = 1, jmax
            do j = 1, 6
              coords(j,turn,i) = z(j,i) - orbit0(j)
            enddo
          enddo
        endif
        if (jmax .eq. 0 .or. (switch .gt. 1 .and. jmax .lt. j_tot))     &
     &goto 20
        if (switch .eq. 2 .and. info) then
          if (mod(turn,100) .eq. 0) print *, 'turn :', turn
        endif
      enddo
!--- end of loop over turns
 20   continue
      if (switch .gt. 1 .and. jmax .lt. j_tot)  then
        e_flag = 1
        return
      endif
      do i = 1, jmax
        last_turn(part_id(i)) = min(turns, turn)
        last_pos(part_id(i)) = sum
        do j = 1, 6
          last_orbit(j,part_id(i)) = z(j,i)
        enddo
      enddo
      turn = min(turn, turns)
!--- enter last turn in tables if not done already
      if (.not. last_out)  then
        if (switch .eq. 1)  then
          if (onetable)  then
!hbu
            spos=sum
!hbu get current node name
            call element_name(el_name,len(el_name))
!hbu spos added
            call tt_putone(jmax, turn, tot_segm, segment, part_id,      &
     &z, orbit0,spos,nlm,el_name)
          else
            do i = 1, jmax
!hbu
              call tt_puttab(part_id(i), turn, 1, z(1,i), orbit0,spos)
            enddo
          endif
        endif
      endif
!--- enter last turn in summary table
      do  i = 1,j_tot
        tmp_d = i
        call double_to_table('tracksumm ', 'number ', tmp_d)
        tmp_d = last_turn(i)
        call double_to_table('tracksumm ', 'turn ', tmp_d)
        do j = 1, 6
          tmp_d = last_orbit(j,i) - orbit0(j)
          call double_to_table('tracksumm ', vec_names(j), tmp_d)
        enddo
!hbu
        spos=last_pos(i)
!hbu
        call double_to_table('tracksumm ',vec_names(7),spos)
        call augment_count('tracksumm ')
      enddo
 999  end
      subroutine ttmap(code,el,track,ktrack,dxt,dyt,sum,turn,part_id,   &
     &last_turn,last_pos, last_orbit,aperflag,maxaper,al_errors,        &
     &onepass)                                                          &
      implicit none

!----------------------------------------------------------------------*
! Purpose:                                                             *
!   -  Test apertures                                                  *
!   -  Track through thin lenses ONLY.                                 *
!                                                                      *
! Input/output:                                                        *
!   SUM       (double)    Accumulated length.                          *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   NUMBER(*) (integer) Number of current track.                       *
!   KTRACK    (integer) number of surviving tracks.                    *
!----------------------------------------------------------------------*
      include 'twiss0.fi'
      include 'name_len.fi'
      include 'twtrr.fi'
      logical aperflag,bend,fmap,onepass
      integer turn,code,ktrack,part_id(*),last_turn(*),nn,jtrk,         &
     &get_option
      double precision apx,apy,el,sum,node_value,track(6,*),            &
     &last_pos(*),last_orbit(6,*),parvec(26),get_value,ct,tmp,          &
     &aperture(maxnaper),one,maxaper(6), zero,al_errors(align_max),     &
     &st,theta,ek(6),re(6,6),te(6,6,6),craporb(6),dxt(*),dyt(*)
      character*(name_len) aptype
      parameter(zero = 0.d0, one=1d0)

      fmap=.false.
      call dzero(ek,6)
      call dzero(craporb,6)
      call dzero(re,36)
      call dzero(te,216)
!---- Drift space
      if(code.eq.1) then
        call ttdrf(el,track,ktrack)
        go to 502
      endif
!---- Translation of particles // AK 23.02.2006
!---- useful for combination of different transfer lines or rings
      if(code.eq.36) then
        if (onepass) then
          call tttrans(track,ktrack)
          go to 502
        endif
      endif
!---- Collimator with elliptic aperture.
      if(code.eq.20) then
        apx = node_value('xsize ')
        apy = node_value('ysize ')
        if(apx.eq.zero) then
          apx=maxaper(1)
        endif
        if(apy.eq.zero) then
          apy=maxaper(3)
        endif
        call trcoll(1, apx, apy, turn, sum, part_id, last_turn,         &
     &last_pos, last_orbit, track, ktrack,al_errors)
        go to 502
      endif
!---- Collimator with rectangular aperture.
      if(code.eq.21) then
        apx = node_value('xsize ')
        apy = node_value('ysize ')
        if(apx.eq.zero) then
          apx=maxaper(1)
        endif
        if(apy.eq.zero) then
          apy=maxaper(3)
        endif
        call trcoll(2, apx, apy, turn, sum, part_id, last_turn,         &
     &last_pos, last_orbit, track, ktrack,al_errors)
        go to 502
      endif
!---- Beam-beam,  standard 4D, no aperture
      if(code.eq.22) then
        parvec(5)=get_value('probe ', 'arad ')
        parvec(6)=node_value('charge ') * get_value('probe ', 'npart ')
        parvec(7)=get_value('probe ','gamma ')
        call ttbb(track, ktrack, parvec)
        go to 502
      endif

!____ Test aperture. ALL ELEMENTS BUT DRIFTS
      aperflag = get_option('aperture ') .ne. 0
      if(aperflag) then
        nn=name_len
        call node_string('apertype ',aptype,nn)
        call dzero(aperture,maxnaper)
        call get_node_vector('aperture ',nn,aperture)
!        print *, " TYPE ",aptype
!     &  "values  x y lhc",aperture(1),aperture(2),aperture(3)
!------------  ellipse case ----------------------------------
        if(aptype.eq.'ellipse') then
          apx = aperture(1)
          apy = aperture(2)
          call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!------------  circle case ----------------------------------
        else if(aptype.eq.'circle') then
          apx = aperture(1)
!        print *,"radius of circle in element",apx
          if(apx.eq.zero) then
            apx = maxaper(1)
!        print *,"radius of circle by default",apx
          endif
          apy = apx
!        print *,"circle, radius= ",apx
          call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!------------  rectangle case ----------------------------------
        else if(aptype.eq.'rectangle') then
          apx = aperture(1)
          apy = aperture(2)
          call trcoll(2, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!------------  LHC screen case ----------------------------------
        else if(aptype.eq.'lhcscreen') then
!        print *, "LHC screen start, Xrect= ",
!     &  aperture(1),"  Yrect= ",aperture(2),"  Rcirc= ",aperture(3)
          apx = aperture(3)
          apy = aperture(3)
!JMJ!     Making essential changes in AV's absence, 16/7/2003
!JMJ!     this tests whether the particle is outside the circumscribing
!JMJ!     circle.
          call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!JMJ!
!JMJ!     This tests whether particles are outside the space bounded by
!JMJ!     two or four lines intersecting the circle.
!JMJ!     The previous version checked that it was outside a rectangle,
!JMJ!     one of whose dimensions was zero (in the LHC) so particles
!JMJ!     were ALWAYS lost on the beam screen.
!JMJ!     The new version of the test works on the understanding that
!JMJ!     values of aperture(1) or aperture(2) greater than aperture(3)
!JMJ!     have no meaning and that specifying a zero value is equivalent.
!JMJ!     The most general aperture described by the present
!JMJ!     implementation of LHCSCREEN is a rectangle with rounded
!JMJ!     corners.
!JMJ!     N.B. apx and apy already have the value aperture(3); the "if"s
!JMJ!     ensure that they don't get set to zero.
          if(aperture(1).gt.0.) apx = aperture(1)
          if(aperture(2).gt.0.) apy = aperture(2)
          call trcoll(2, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!        print *, "LHC screen end"
!------------  marguerite case ----------------------------------
        else if(aptype.eq.'marguerite') then
          apx = aperture(1)
          apy = aperture(2)
          call trcoll(3, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!------------  rectellipse case ----------------------------------
        else if(aptype.eq.'rectellipse') then
!*****         test ellipse
          apx = aperture(3)
          apy = aperture(4)
          call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!*****         test rectangle
          apx = aperture(1)
          apy = aperture(2)
          call trcoll(2, apx, apy, turn, sum, part_id, last_turn,       &
     &last_pos, last_orbit, track, ktrack,al_errors)
!       print*, " test apertures"
!       print*, "      apx=",apx, " apy=",apy," apxell=",apxell,
!     &         " apyell=",apyell
!          call trcoll(3, apx, apy, turn, sum, part_id, last_turn,       &
!     &    last_pos, last_orbit, track, ktrack,al_errors)
        endif
      endif
!----  END OF IF(APERFLAG)  ----------------
      bend = code .eq. 2 .or. code .eq. 3
      if (.not. bend)  then
        theta = node_value('tilt ')
        if (theta .ne. zero)  then
!--- rotate trajectory before entry
          st = sin(theta)
          ct = cos(theta)
          do jtrk = 1,ktrack
            tmp = track(1,jtrk)
            track(1,jtrk) = ct * tmp + st * track(3,jtrk)
            track(3,jtrk) = ct * track(3,jtrk) - st * tmp
            tmp = track(2,jtrk)
            track(2,jtrk) = ct * tmp + st * track(4,jtrk)
            track(4,jtrk) = ct * track(4,jtrk) - st * tmp
          enddo
        endif
      endif
!-- switch on element type BUT DRIFT, COLLIMATORS, BEAM_BEAM / 13.03.03
!-- 500 has been specified at the relevant places in go to below
!-- code =1 for drift, treated above, go to 500 directly
!      print *,"   CODE    ",code
      go to ( 500,  20,  30,  40,  50,  60,  70,  80,  90, 100,         &
     &110, 120, 130, 140, 150, 160, 170, 180, 190, 500,                 &
     &500, 500, 230, 240, 250, 260, 270, 280, 290, 300,                 &
     &310, 320, 330, 500, 350, 360, 500, 500, 500, 500), code
!
!---- Make sure that nothing is execute if element is not known
      go to 500
!
!---- Bending magnet. OBSOLETE, to be kept for go to
   20 continue  ! RBEND
   30 continue  ! SBEND
      go to 500
!---- Arbitrary matrix. OBSOLETE, to be kept for go to
   40 continue
      call tmarb(.false.,.false.,craporb,fmap,el,ek,re,te)
      call tttrak(ek,re,te,track,ktrack)
      go to 500
!---- Quadrupole. OBSOLETE, to be kept for go to
   50 continue
      go to 500
!---- Sextupole. OBSOLETE, to be kept for go to
   60 continue
      go to 500
!---- Octupole. OBSOLETE, to be kept for go to
   70 continue
      go to 500
!---- Monitors, beam instrument., MONITOR, HMONITOR, VMONITOR
  170 continue
  180 continue
  190 continue
      go to 500
!---- Multipole
   80 continue
      call ttmult(track,ktrack,dxt,dyt)
      go to 500
!---- Solenoid.
   90 continue
      call trsol(track, ktrack)
      go to 500
!---- RF cavity.
  100 continue
      call ttrf(track,ktrack)
      go to 500
!---- Electrostatic separator.
  110 continue
!        call ttsep(el, track, ktrack)
      go to 500
!---- Rotation around s-axis.
  120 continue
      call ttsrot(track, ktrack)
      go to 500
!---- Rotation around y-axis.
  130 continue
!        call ttyrot(track, ktrack)
      go to 500
!---- Correctors.
  140 continue
  150 continue
  160 continue
      call ttcorr(el, track, ktrack)
      go to 500
!---- ECollimator, RCollimator, BeamBeam, Lump. ??
  200 continue
  210 continue
  220 continue
  230 continue
      go to 500
!---- Instrument
  240 continue
      go to 500
!---- Marker.
  250 continue
      go to 500
!---- General bend (dipole, quadrupole, and skew quadrupole).
  260 continue
      go to 500
!---- LCAV cavity.
  270 continue
!        call ttlcav(el, track, ktrack)
      go to 500
!---- Reserved. (PROFILE,WIRE,SLMONITOR,BLMONITOR,IMONITOR)
  280 continue
  290 continue
  300 continue
  310 continue
  320 continue
      go to 500
!---- Dipole edge
  330 continue
      call ttdpdg(track,ktrack)
      go to 500
!---- Changeref ???
  350 continue
      go to 500
!---- Translation
  360 continue
      go to 500


  500 continue
      if (.not. bend)  then
        if (theta .ne. zero)  then
!--- rotate trajectory at exit
          do jtrk = 1,ktrack
            tmp = track(1,jtrk)
            track(1,jtrk) = ct * tmp - st * track(3,jtrk)
            track(3,jtrk) = ct * track(3,jtrk) + st * tmp
            tmp = track(2,jtrk)
            track(2,jtrk) = ct * tmp - st * track(4,jtrk)
            track(4,jtrk) = ct * track(4,jtrk) + st * tmp
          enddo
        endif
      endif
!---- Accumulate length.
 502  sum = sum + el
      return
      end

      subroutine ttmult(track, ktrack,dxt,dyt)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!    Track particle through a general thin multipole.                  *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer) Number of surviving tracks.                    *
!   dxt       (double)  local buffer                                   *
!   dyt       (double)  local buffer                                   *
!----------------------------------------------------------------------*
      include 'twtrr.fi'
      include 'track.fi'
      logical first
      integer iord,jtrk,nd,nord,ktrack,j,n_ferr,nn,ns,node_fd_errors,   &
     &get_option
      double precision     const,curv,dbi,dbr,dipi,dipr,dx,dy,elrad,    &
     &pt,px,py,rfac,rpt1,rpt2,rpx1,rpx2,rpy1,rpy2,                      &
     &f_errors(0:maxferr),field(2,0:maxmul),vals(2,0:maxmul),           &
     &ordinv(maxmul),track(6,*),dxt(*),dyt(*),normal(0:maxmul),         &
     &skew(0:maxmul),bv0,bvk,node_value,zero,one,two,three,half,        &
     &ttt

      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,half=5d-1)
      save first,ordinv
      data first / .true. /

!---- Precompute reciprocals of orders.
      if (first) then
        do iord = 1, maxmul
          ordinv(iord) = one / float(iord)
        enddo
        first = .false.
      endif
      call dzero(f_errors, maxferr+1)
      n_ferr = node_fd_errors(f_errors)
      bv0 = node_value('dipole_bv ')
      bvk = node_value('other_bv ')
!---- Multipole length for radiation.
      elrad = node_value('lrad ')
!---- Multipole components.
      call dzero(normal,maxmul+1)
      call dzero(skew,maxmul+1)
      call get_node_vector('knl ',nn,normal)
      call get_node_vector('ksl ',ns,skew)
      nd = 2 * max(nn, ns)
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
!---- Dipole error.
!      dbr = bv0 * field(1,0) / (one + deltas)
!      dbi = bv0 * field(2,0) / (one + deltas)
      dbr = bv0 * field(1,0)
      dbi = bv0 * field(2,0)
!---- Nominal dipole strength.
!      dipr = bv0 * vals(1,0) / (one + deltas)
!      dipi = bv0 * vals(2,0) / (one + deltas)
      dipr = bv0 * vals(1,0)
      dipi = bv0 * vals(2,0)
!---- Other components and errors.
      nord = 0
      do iord = 1, nd/2
        do j = 1, 2
!          field(j,iord) = bvk * (vals(j,iord) + field(j,iord))          &
!     &/ (one + deltas)
          field(j,iord) = bvk * (vals(j,iord) + field(j,iord))          
          if (field(j,iord) .ne. zero)  nord = iord
        enddo
      enddo
!---- Pure dipole: only quadrupole kicks according to lrad.
      if (nord .eq. 0) then
        do jtrk = 1,ktrack
          dxt(jtrk) = zero
          dyt(jtrk) = zero
        enddo
!----------- introduction of dipole focusing
        if(elrad.gt.zero.and.get_option('thin_foc ').eq.1) then
          do jtrk = 1,ktrack
            dxt(jtrk) =  dipr*dipr*track(1,jtrk)/elrad
            dyt(jtrk) =  dipi*dipi*track(3,jtrk)/elrad
          enddo
        endif
!---- Accumulate multipole kick from highest multipole to quadrupole.
      else
        do jtrk = 1,ktrack
          dxt(jtrk) =                                                   &
     &field(1,nord)*track(1,jtrk) - field(2,nord)*track(3,jtrk)
          dyt(jtrk) =                                                   &
     &field(1,nord)*track(3,jtrk) + field(2,nord)*track(1,jtrk)
        enddo

        do iord = nord - 1, 1, -1
          do jtrk = 1,ktrack
            dx = dxt(jtrk)*ordinv(iord+1) + field(1,iord)
            dy = dyt(jtrk)*ordinv(iord+1) + field(2,iord)
            dxt(jtrk) = dx*track(1,jtrk) - dy*track(3,jtrk)
            dyt(jtrk) = dx*track(3,jtrk) + dy*track(1,jtrk)
          enddo
        enddo
!        do jtrk = 1,ktrack
!          dxt(jtrk) = dxt(jtrk) / (one + deltas)
!          dyt(jtrk) = dyt(jtrk) / (one + deltas)
!        enddo
        if(elrad.gt.zero.and.get_option('thin_foc ').eq.1) then
          do jtrk = 1,ktrack
            dxt(jtrk) = dxt(jtrk) + dipr*dipr*track(1,jtrk)/elrad
            dyt(jtrk) = dyt(jtrk) + dipi*dipi*track(3,jtrk)/elrad
          enddo
        endif
      endif

!---- Radiation loss at entrance.
      if (dorad .and. elrad .ne. 0) then
        const = arad * gammas**3 / three

!---- Full damping.
        if (dodamp) then
          do jtrk = 1,ktrack
            curv = sqrt((dipr + dxt(jtrk))**2 +                         &
     &(dipi + dyt(jtrk))**2) / elrad

            if (dorand) then
              call trphot(elrad,curv,rfac,deltas)
            else
              rfac = const * curv**2 * elrad
            endif

            px = track(2,jtrk)
            py = track(4,jtrk)
            pt = track(6,jtrk)
            track(2,jtrk) = px - rfac * (one + pt) * px
            track(4,jtrk) = py - rfac * (one + pt) * py
            track(6,jtrk) = pt - rfac * (one + pt) ** 2
          enddo

!---- Energy loss like for closed orbit.
        else

!---- Store energy loss on closed orbit.
          rfac = const * ((dipr + dxt(1))**2 + (dipi + dyt(1))**2)
          rpx1 = rfac * (one + track(6,1)) * track(2,1)
          rpy1 = rfac * (one + track(6,1)) * track(4,1)
          rpt1 = rfac * (one + track(6,1)) ** 2

          do jtrk = 1,ktrack
            track(2,jtrk) = track(2,jtrk) - rpx1
            track(4,jtrk) = track(4,jtrk) - rpy1
            track(6,jtrk) = track(6,jtrk) - rpt1
          enddo

        endif
      endif

!---- Apply multipole effect including dipole.
      do jtrk = 1,ktrack
!       Added for correct Ripken implementation of formulae      	
        ttt = sqrt(one+two*track(6,jtrk)*bet0i+track(6,jtrk)**2)
!        track(2,jtrk) = track(2,jtrk) -                                 &
!     &(dbr + dxt(jtrk) - dipr * (deltas + beti*track(6,jtrk)))
!        track(4,jtrk) = track(4,jtrk) +                                 &
!     &(dbi + dyt(jtrk) - dipi * (deltas + beti*track(6,jtrk)))
!        track(5,jtrk) = track(5,jtrk)                                   &
!     &- (dipr*track(1,jtrk) - dipi*track(3,jtrk)) * beti
        track(2,jtrk) = track(2,jtrk) -                                 &
     &   (dbr + dxt(jtrk) - dipr * (ttt - one))
        track(4,jtrk) = track(4,jtrk) +                                 &
     &   (dbi + dyt(jtrk) - dipi * (ttt - one))
        track(5,jtrk) = track(5,jtrk) -                                 &
     &   (dipr*track(1,jtrk) - dipi*track(3,jtrk)) *                    &
     &   ((one + bet0*track(6,jtrk))/ttt)*bet0i
      enddo

!---- Radiation loss at exit.
      if (dorad .and. elrad .ne. 0) then

!---- Full damping.
        if (dodamp) then
          do jtrk = 1,ktrack
            curv = sqrt((dipr + dxt(jtrk))**2 +                         &
     &(dipi + dyt(jtrk))**2) / elrad

            if (dorand) then
              call trphot(elrad,curv,rfac,deltas)
            else
              rfac = const * curv**2 * elrad
            endif

            px = track(2,jtrk)
            py = track(4,jtrk)
            pt = track(6,jtrk)
            track(2,jtrk) = px - rfac * (one + pt) * px
            track(4,jtrk) = py - rfac * (one + pt) * py
            track(6,jtrk) = pt - rfac * (one + pt) ** 2
          enddo

!---- Energy loss like for closed orbit.
        else

!---- Store energy loss on closed orbit.
          rfac = const * ((dipr + dxt(1))**2 + (dipi + dyt(1))**2)
          rpx2 = rfac * (one + track(6,1)) * track(2,1)
          rpy2 = rfac * (one + track(6,1)) * track(4,1)
          rpt2 = rfac * (one + track(6,1)) ** 2

          do jtrk = 1,ktrack
            track(2,jtrk) = track(2,jtrk) - rpx2
            track(4,jtrk) = track(4,jtrk) - rpy2
            track(6,jtrk) = track(6,jtrk) - rpt2
          enddo
        endif
      endif

      end
      subroutine trphot(el,curv,rfac,deltap)
      implicit none

!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Generate random energy loss for photons, using a look-up table to  *
!   invert the function Y.  Ultra-basic interpolation computed;        *
!   leads to an extrapolation outside the table using the two outmost  *
!   point on each side (low and high).                                 *
!   Assumes ultra-relativistic particles (beta = 1).                   *
! Author: Ghislain Roy                                                 *
! Input:                                                               *
!   EL     (double)       Element length.                              *
!   CURV   (double)       Local curvature of orbit.                    *
! Output:                                                              *
!   RFAC   (double)       Relative energy loss due to photon emissions.*
!----------------------------------------------------------------------*
!---- Generate pseudo-random integers in batches of NR.                *
!     The random integers are generated in the range [0, MAXRAN).      *
!---- Table definition: maxtab, taby(maxtab), tabxi(maxtab)            *
!----------------------------------------------------------------------*
      integer i,ierror,j,nphot,nr,maxran,maxtab
      parameter(nr=55,maxran=1000000000,maxtab=101)
      double precision amean,curv,dlogr,el,frndm,rfac,scalen,scaleu,    &
     &slope,ucrit,xi,deltap,hbar,clight,arad,pc,gamma,amass,real_am,    &
     &get_value,get_variable,tabxi(maxtab),taby(maxtab),zero,one,two,   &
     &three,five,twelve,fac1
      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,five=5d0,twelve=12d0,&
     &fac1=3.256223d0)
      character*20 text
      data (taby(i), i = 1, 52)                                         &
     &/ -1.14084005d0,  -0.903336763d0, -0.769135833d0, -0.601840854d0, &
     &-0.448812515d0, -0.345502228d0, -0.267485678d0, -0.204837948d0,   &
     &-0.107647471d0, -0.022640628d0,  0.044112321d0,  0.0842842236d0,  &
     &0.132941082d0,  0.169244036d0,  0.196492359d0,  0.230918407d0,    &
     &0.261785239d0,  0.289741248d0,  0.322174788d0,  0.351361096d0,    &
     &0.383441716d0,  0.412283719d0,  0.442963421d0,  0.472622454d0,    &
     &0.503019691d0,  0.53197819d0,   0.561058342d0,  0.588547111d0,    &
     &0.613393188d0,  0.636027336d0,  0.675921738d0,  0.710166812d0,    &
     &0.725589216d0,  0.753636241d0,  0.778558254d0,  0.811260045d0,    &
     &0.830520391d0,  0.856329501d0,  0.879087269d0,  0.905612588d0,    &
     &0.928626955d0,  0.948813677d0,  0.970829248d0,  0.989941061d0,    &
     &1.0097903d0,    1.02691281d0,   1.04411256d0,   1.06082714d0,     &
     &1.0750246d0,    1.08283985d0,   1.0899564d0,    1.09645379d0 /
      data (taby(i), i = 53, 101)                                       &
     &/  1.10352755d0,   1.11475027d0,   1.12564385d0,   1.1306442d0,   &
     &1.13513422d0,   1.13971806d0,   1.14379156d0,   1.14741969d0,     &
     &1.15103698d0,   1.15455759d0,   1.15733826d0,   1.16005647d0,     &
     &1.16287541d0,   1.16509759d0,   1.16718769d0,   1.16911888d0,     &
     &1.17075884d0,   1.17225218d0,   1.17350936d0,   1.17428589d0,     &
     &1.17558432d0,   1.17660713d0,   1.17741513d0,   1.17805469d0,     &
     &1.17856193d0,   1.17896497d0,   1.17928565d0,   1.17954147d0,     &
     &1.17983139d0,   1.1799767d0,    1.18014216d0,   1.18026078d0,     &
     &1.18034601d0,   1.1804074d0,    1.18045175d0,   1.1804837d0,      &
     &1.18051291d0,   1.18053186d0,   1.18054426d0,   1.18055236d0,     &
     &1.18055761d0,   1.18056166d0,   1.18056381d0,   1.1805656d0,      &
     &1.18056655d0,   1.18056703d0,   1.18056726d0,   1.1805675d0,      &
     &1.18056762d0 /
      data (tabxi(i), i = 1, 52)                                        &
     &/ -7.60090017d0,  -6.90775537d0,  -6.50229025d0,  -5.99146461d0,  &
     &-5.52146101d0,  -5.20300722d0,  -4.96184492d0,  -4.76768923d0,    &
     &-4.46540833d0,  -4.19970512d0,  -3.98998451d0,  -3.86323285d0,    &
     &-3.70908213d0,  -3.59356928d0,  -3.50655794d0,  -3.39620972d0,    &
     &-3.29683733d0,  -3.20645332d0,  -3.10109282d0,  -3.0057826d0,     &
     &-2.9004221d0,   -2.80511189d0,  -2.70306253d0,  -2.60369015d0,    &
     &-2.50103593d0,  -2.4024055d0 ,  -2.30258512d0,  -2.20727491d0,    &
     &-2.12026358d0,  -2.04022098d0,  -1.89712d0   ,  -1.7719568d0,     &
     &-1.71479833d0,  -1.60943794d0,  -1.51412773d0,  -1.38629436d0,    &
     &-1.30933332d0,  -1.20397282d0,  -1.10866261d0,  -0.99425226d0,    &
     &-0.89159810d0,  -0.79850775d0,  -0.69314718d0,  -0.59783697d0,    &
     &-0.49429631d0,  -0.40047753d0,  -0.30110508d0,  -0.19845095d0,    &
     &-0.10536054d0,  -0.05129330d0,   0.0d0,          0.048790119d0 /
      data (tabxi(i), i = 53, 101)                                      &
     &/  0.104360029d0,  0.198850885d0,  0.300104618d0,  0.350656837d0, &
     &0.398776114d0,  0.451075643d0,  0.500775278d0,  0.548121393d0,    &
     &0.598836541d0,  0.652325153d0,  0.69813472d0 ,  0.746687889d0,    &
     &0.802001595d0,  0.850150883d0,  0.900161386d0,  0.951657832d0,    &
     &1.00063193d0,   1.05082154d0,   1.09861231d0,   1.13140213d0,     &
     &1.1939224d0,    1.25276291d0,   1.3083328d0,    1.36097658d0,     &
     &1.4109869d0,    1.45861506d0,   1.50407743d0,   1.54756248d0,     &
     &1.60943794d0,   1.64865863d0,   1.70474803d0,   1.75785792d0,     &
     &1.80828881d0,   1.85629797d0,   1.90210748d0,   1.9459101d0,      &
     &2.0014801d0,    2.05412364d0,   2.10413408d0,   2.15176225d0,     &
     &2.19722462d0,   2.25129175d0,   2.29253483d0,   2.35137534d0,     &
     &2.40694523d0,   2.45100522d0,   2.501436d0,     2.60268974d0,     &
     &2.64617491d0 /

!Get constants
      clight = get_variable('clight ')
      hbar   = get_variable('hbar ')
      arad   = get_value('probe ','arad ')
      pc     = get_value('probe ','pc ')
      amass  = get_value('probe ','mass ')
      gamma  = get_value('probe ','gamma ')
      scalen = five / (twelve * hbar * clight)
      scaleu = hbar * three * clight / two

!---- AMEAN is the average number of photons emitted.,
!     NPHOT is the integer number generated from Poisson's law.
      amean = scalen * abs(arad*pc*(one+deltap)*el*curv) * sqrt(three)
      rfac = zero
      real_am = amean
      if (real_am .gt. zero) then
        call dpoissn(real_am, nphot, ierror)

        if (ierror .ne. 0) then
          write(text, '(1p,d20.12)') amean
          call aafail('TRPHOT: ','Fatal: Poisson input mean =' // text)
        endif

!---- For all photons, sum the radiated photon energy,
!     in units of UCRIT (relative to total energy).

        if (nphot .ne. 0) then
          ucrit = scaleu * gamma**2 * abs(curv) / amass
          xi = zero
          do i = 1, nphot

!---- Find a uniform random number in the range [ 0,3.256223 ].
!     Note that the upper limit is not exactly 15*sqrt(3)/8
!     because of imprecision in the integration of F.
            dlogr = log(fac1 * frndm())

!---- Now look for the energy of the photon in the table TABY/TABXI
            do j = 2, maxtab
              if (dlogr .le. taby(j) ) go to 20
            enddo

!---- Perform linear interpolation and sum up energy lost.
 20         slope = (dlogr - taby(j-1)) / (taby(j) - taby(j-1))
            xi = dexp(tabxi(j-1) + slope * (tabxi(j) - tabxi(j-1)))
            rfac = rfac + ucrit * xi
          enddo
        endif
      endif
      end

      subroutine ttsrot(track,ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Coordinate change due to rotation about s axis                     *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer) number of surviving tracks.        tmsrot      *
!----------------------------------------------------------------------*
      include 'track.fi'
      integer itrack,ktrack,j
      double precision track(6,*),psi,node_value,ct,st,trb(4)
      psi = node_value('angle ')
      ct = cos(psi)
      st = -sin(psi)
      do  itrack = 1, ktrack
        do  j = 1,4
          trb(j)=track(j,itrack)
        enddo
        track(1,itrack) = trb(1)*ct-trb(3)*st
        track(2,itrack) = trb(2)*ct-trb(4)*st
        track(3,itrack) = trb(1)*st+trb(3)*ct
        track(4,itrack) = trb(2)*st+trb(4)*ct
      enddo
      end


      subroutine ttyrot(track,ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Coordinate change due to rotation about s axis                     *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer) number of surviving tracks.        tmsrot      *
!----------------------------------------------------------------------*
      include 'track.fi'
      integer itrack,ktrack,j
      double precision track(6,*),theta,node_value,ct,st,trb(6)
      theta = node_value('angle ')
      ct = cos(theta)
      st = -sin(theta)
      do  itrack = 1, ktrack
        do  j = 1,6
          trb(j)=track(j,itrack)
        enddo
        track(1,itrack) = trb(1)*ct-trb(5)*st
        track(2,itrack) = trb(2)*ct-trb(6)*st
        track(5,itrack) = trb(1)*st+trb(5)*ct
        track(6,itrack) = trb(2)*st+trb(6)*ct
      enddo
      end

      subroutine ttdrf(el,track,ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Track a set of particle through a drift space.                     *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer) number of surviving tracks.                    *
! Output:                                                              *
!   EL        (double)    Length of quadrupole.                        *
!----------------------------------------------------------------------*
      include 'track.fi'
      integer itrack,ktrack
      double precision el,pt,px,py,track(6,*),ttt,one,two
      parameter(one=1d0,two=2d0)

! picked from trturn in madx.ss
      do  itrack = 1, ktrack
        px = track(2,itrack)
        py = track(4,itrack)
        pt = track(6,itrack)
        ttt = el/sqrt(one+two*pt*bet0i+pt**2 - px**2 - py**2)
        track(1,itrack) = track(1,itrack) + ttt*px
        track(3,itrack) = track(3,itrack) + ttt*py
!        track(5,itrack) = track(5,itrack)                               &
!     &+ el*(beti + pt * dtbyds) - (beti+pt)*ttt
!---- AK 20060413
!---- Ripken DESY-95-189 p.36
        track(5,itrack) = track(5,itrack)                               &
     &                  + bet0i*(el - (one + bet0*pt) * ttt)
      enddo
      end
      subroutine ttrf(track,ktrack)

      implicit none

!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Track a set of trajectories through a thin cavity (zero length).   *
!   The cavity is sandwiched between two drift spaces of half length.  *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer) number of surviving tracks.                    *
! Output:                                                              *
!   EL        (double)    Length of quadrupole.                        *
!----------------------------------------------------------------------*
! Modified: 06-JAN-1999, T. Raubenheimer (SLAC)                        *
!   Modified to allow wakefield tracking however no modification to    *
!   the logic and the nominal energy is not updated -- see routine     *
!   TTLCAV to change the nominal energy                                *
!----------------------------------------------------------------------*
      integer itrack,ktrack
      double precision bi2gi2,dl,el,omega,dtbyds,phirf,pt,rff,betas,    &
     &gammas,rfl,rfv,track(6,*),clight,twopi,vrf,deltap,deltas,pc,pc0,  &
     &get_variable,node_value,get_value,one,two,half,ten3m,ten6p
!      double precision px,py,ttt,beti,el1
      parameter(one=1d0,two=2d0,half=5d-1,ten3m=1d-3,ten6p=1d6)

!---- Initialize
      clight=get_variable('clight ')
      twopi=get_variable('twopi ')

!---- Fetch data.
!      el = node_value('l ')
!      el1 = node_value('l ')
      rfv = node_value('volt ')
      rff = node_value('freq ')
      rfl = node_value('lag ')
!      deltap = get_value('probe ','deltap ')
!      deltas = get_variable('track_deltap ')
!      pc = get_value('probe ','pc ')
      pc0 = get_value('beam ','pc ')
!      betas = get_value('probe ','beta ')
!      gammas= get_value('probe ','gamma ')
!      dtbyds = get_value('probe ','dtbyds ')

!      print *,"RF cav.  volt=",rfv, "  freq.",rff

!*---- Get the longitudinal wakefield filename (parameter #17).
!      if (iq(lcelm+melen+15*mcsiz-2) .eq. 61) then
!        lstr = iq(lcelm+melen+15*mcsiz)
!        call uhtoc(iq(lq(lcelm-17)+1), mcwrd, lfile, 80)
!      else
!        lfile = ' '
!      endif

!*---- Get the transverse wakefield filename (parameter #18).
!      if (iq(lcelm+melen+16*mcsiz-2) .eq. 61) then
!        lstr = iq(lcelm+melen+16*mcsiz)
!        call uhtoc(iq(lq(lcelm-18)+1), mcwrd, tfile, 80)
!      else
!        tfile = ' '
!      endif

!*---- If there are wakefields split the cavity.
!      if (lfile .ne. ' ' .or. tfile .ne. ' ') then
!        el1 = el / two
!        rfv = rfv / two
!        lwake = .true.
!      else
!        el1 = el
!        lwake = .false.
!      endif
!---- Set up.
      omega = rff * (ten6p * twopi / clight)
!      vrf   = rfv * ten3m / (pc * (one + deltas))
      vrf   = rfv * ten3m / pc0
      phirf = rfl * twopi
!      dl    = el * half
!      bi2gi2 = one / (betas * gammas) ** 2
!---- Loop for all particles.
!      do itrack = 1, ktrack
!        pt = track(6,itrack)
!        pt = pt + vrf * sin(phirf - omega * track(5,itrack))
!        track(6,itrack) = pt
!!      print *," pt ",pt, "   track(5,itrack)",track(5,itrack)
!      enddo
      do itrack = 1, ktrack
        pt = track(6,itrack)
        pt = pt + vrf * sin(phirf - omega * track(5,itrack))
        track(6,itrack) = pt
!      print *," pt ",pt, "   track(5,itrack)",track(5,itrack)
      enddo

!*---- If there were wakefields, track the wakes and then the 2nd half
!*     of the cavity.
!      if (lwake) then
!        call ttwake(two*el1, nbin, binmax, lfile, tfile, ener1, track,
!     +              ktrack)
!
!*---- Track 2nd half of cavity -- loop for all particles.
!      do 20 itrack = 1, ktrack
!
!*---- Drift to centre.
!         px = track(2,itrack)
!         py = track(4,itrack)
!         pt = track(6,itrack)
!         ttt = one/sqrt(one+two*pt*beti+pt**2 - px**2 - py**2)
!         track(1,itrack) = track(1,itrack) + dl*ttt*px
!         track(3,itrack) = track(3,itrack) + dl*ttt*py
!         track(5,itrack) = track(5,itrack)
!     +        + dl*(beti - (beti+pt)*ttt) + dl*pt*dtbyds
!
!*---- Acceleration.
!         pt = pt + vrf * sin(phirf - omega * track(5,itrack))
!         track(6,itrack) = pt
!
!*---- Drift to end.
!         ttt = one/sqrt(one+two*pt*beti+pt**2 - px**2 - py**2)
!         track(1,itrack) = track(1,itrack) + dl*ttt*px
!         track(3,itrack) = track(3,itrack) + dl*ttt*py
!         track(5,itrack) = track(5,itrack)
!     +        + dl*(beti - (beti+pt)*ttt) + dl*pt*dtbyds
! 20   continue
!      endif
      end
      subroutine ttcorr(el,track,ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Track particle through an orbit corrector.                         *
!   The corrector is sandwiched between two half-length drift spaces.  *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer) number of surviving tracks.                    *
! Output:                                                              *
!   EL        (double)    Length of quadrupole.                        *
!----------------------------------------------------------------------*
      include 'track.fi'
      include 'twtrr.fi'
!      logical dorad,dodamp,dorand
      integer itrack,ktrack,n_ferr,node_fd_errors,code,bv0,i,get_option
      double precision bi2gi2,bil2,curv,d,dpx,dpy,el,pt,px,py,rfac,rpt, &
     &rpx,rpy,track(6,*),xkick,ykick, !deltas,arad,betas,gammas,dtbyds, &
     &div,f_errors(0:maxferr),field(2),get_variable,get_value,          &
     &node_value,zero,one,two,three,half
      double precision dpxx,dpyy,ttt,ddd,xp,yp,pxx,pyy
      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,half=5d-1)

!---- Initialize.
      bv0 = node_value('dipole_bv ')
      deltas = get_variable('track_deltap ')
      arad = get_value('probe ','arad ')
      betas = get_value('probe ','beta ')
      gammas = get_value('probe ','gamma ')
      dtbyds = get_value('probe ','dtbyds ')
      code = node_value('mad8_type ')
      call dzero(f_errors, maxferr+1)
      n_ferr = node_fd_errors(f_errors)
      dorad = get_value('probe ','radiate ') .ne. 0
      dodamp = get_option('damp ') .ne. 0
      dorand = get_option('quantum ') .ne. 0
      if (el .eq. zero)  then
        div = one
      else
        div = el
      endif

      do i = 1, 2
        field(i) = zero
      enddo
      if (n_ferr .gt. 0) call dcopy(f_errors, field, min(2, n_ferr))
      if (code.eq.14) then
        xkick=bv0*(node_value('kick ')+node_value('chkick ')+           &
     &field(1)/div)
        ykick=zero
      else if (code.eq.15) then
        xkick=bv0*(node_value('hkick ')+node_value('chkick ')+          &
     &field(1)/div)
        ykick=bv0*(node_value('vkick ')+node_value('cvkick ')+          &
     &field(2)/div)
      else if(code.eq.16) then
        xkick=zero
        ykick=bv0*(node_value('kick ')+node_value('cvkick ')+           &
     &field(2)/div)
      else
        xkick=zero
        ykick=zero
      endif

!---- Sum up total kicks.
      dpx = xkick / (one + deltas)
      dpy = ykick / (one + deltas)
!---- Thin lens 6D version
      dpxx = xkick
      dpyy = ykick

      bil2 = el / (two * betas)
      bi2gi2 = one / (betas * gammas) ** 2

!---- Half radiation effects at entrance.
      if (dorad  .and.  el .ne. 0) then
        if (dodamp .and. dorand) then
          curv = sqrt(dpx**2 + dpy**2) / el
        else
          rfac = arad * gammas**3 * (dpx**2 + dpy**2) / (three * el)
        endif

!---- Full damping.
        if (dodamp) then
          do itrack = 1, ktrack
            if (dorand) call trphot(el,curv,rfac,deltas)
            px = track(2,itrack)
            py = track(4,itrack)
            pt = track(6,itrack)
            track(2,itrack) = px - rfac * (one + pt) * px
            track(4,itrack) = py - rfac * (one + pt) * py
            track(6,itrack) = pt - rfac * (one + pt) ** 2
          enddo

!---- Energy loss as for closed orbit.
        else
          rpx = rfac * (one + track(6,1)) * track(2,1)
          rpy = rfac * (one + track(6,1)) * track(4,1)
          rpt = rfac * (one + track(6,1)) ** 2

          do itrack = 1, ktrack
            track(2,itrack) = track(2,itrack) - rpx
            track(4,itrack) = track(4,itrack) - rpy
            track(6,itrack) = track(6,itrack) - rpt
          enddo
        endif
      endif

!---- Thick lens code taken out! AK 21.04.2006
!!---- Half kick at entrance.
!      do itrack = 1, ktrack
!        px = track(2,itrack) + half * dpx
!        py = track(4,itrack) + half * dpy
!        pt = track(6,itrack)
!
!!---- Drift through corrector.
!        d = (one - pt / betas) * el
!        track(1,itrack) = track(1,itrack) + px * d
!        track(3,itrack) = track(3,itrack) + py * d
!        track(5,itrack) = track(5,itrack) + d * bi2gi2 * pt -           &
!     &bil2 * (px**2 + py**2 + bi2gi2*pt**2) + el*pt*dtbyds
!
!!---- Half kick at exit.
!        track(2,itrack) = px + half * dpx
!        track(4,itrack) = py + half * dpy
!        track(6,itrack) = pt
!      enddo

!---- Kick at dipole corrector magnet
!     including PT-dependence 
      do itrack = 1, ktrack
        px = track(2,itrack)
        py = track(4,itrack)
        pt = track(6,itrack)
!        ttt = sqrt(one+two*pt*bet0i+pt**2 - px**2 - py**2)
!        ddd = sqrt(one+two*pt*bet0i+pt**2)
!        track(2,itrack) = px + dpxx*ttt/ddd
!        track(4,itrack) = py + dpyy*ttt/ddd
!        ttt = sqrt(one+two*pt*bet0i+pt**2 - dpxx**2 - dpyy**2)
!        ddd = sqrt(one+two*pt*bet0i+pt**2)
!        track(2,itrack) = px + dpxx*ttt/ddd
!        track(4,itrack) = py + dpyy*ttt/ddd

!        ttt = sqrt(one+two*pt*bet0i+pt**2 - px**2 - py**2)
!        ddd = sqrt(one+two*pt*bet0i+pt**2)
!
!        xp = px/ttt + dpxx/ddd  ! Apply kick in standard coordinates!
!        yp = py/ttt + dpyy/ddd
!
!        ttt = sqrt(one + xp**2 + yp**2)
!        
!        px = xp*ddd/ttt ! Transform back to canonical coordinates
!        py = yp*ddd/ttt
!        
!        track(2,itrack) = px
!        track(4,itrack) = py

        track(2,itrack) = px + dpxx
        track(4,itrack) = py + dpyy

!       Add time of flight effects (stolen from Ripken-Dipole)
!        track(5,itrack) = track(5,itrack) -                             &
!     &   (dpxx*track(1,itrack) - dpyy*track(3,itrack)) *                &
!     &   ((one + bet0*track(6,itrack))/ddd)*bet0i
        
      enddo

!---- Half radiation effects at exit.
!     If not random, use same RFAC as at entrance.
      if (dorad  .and.  el .ne. 0) then

!---- Full damping.
        if (dodamp) then
          do itrack = 1, ktrack
            if (dorand) call trphot(el,curv,rfac,deltas)
            px = track(2,itrack)
            py = track(4,itrack)
            pt = track(6,itrack)
            track(2,itrack) = px - rfac * (one + pt) * px
            track(4,itrack) = py - rfac * (one + pt) * py
            track(6,itrack) = pt - rfac * (one + pt) ** 2
          enddo

!---- Energy loss as for closed orbit.
        else
          rpx = rfac * (one + track(6,1)) * track(2,1)
          rpy = rfac * (one + track(6,1)) * track(4,1)
          rpt = rfac * (one + track(6,1)) ** 2

          do itrack = 1, ktrack
            track(2,itrack) = track(2,itrack) - rpx
            track(4,itrack) = track(4,itrack) - rpy
            track(6,itrack) = track(6,itrack) - rpt
          enddo
        endif
      endif

      end
      subroutine dpoissn (amu,n,ierror)
      implicit none
!----------------------------------------------------------------------*
!    POISSON GENERATOR                                                 *
!    CODED FROM LOS ALAMOS REPORT      LA-5061-MS                      *
!    PROB(N)=EXP(-AMU)*AMU**N/FACT(N)                                  *
!        WHERE FACT(N) STANDS FOR FACTORIAL OF N                       *
!    ON RETURN IERROR.EQ.0 NORMALLY                                    *
!              IERROR.EQ.1 IF AMU.LE.0.                                *
!----------------------------------------------------------------------*
      integer n, ierror
      double precision amu,ran,pir,frndm,grndm,expma,amax,zero,one,half
      parameter(zero=0d0,one=1d0,half=5d-1)
!    AMAX IS THE VALUE ABOVE WHICH THE NORMAL DISTRIBUTION MUST BE USED
      data amax/88d0/

      ierror= 0
      if(amu.le.zero) then
!    MEAN SHOULD BE POSITIVE
        ierror=1
        n = 0
        go to 999
      endif
      if(amu.gt.amax) then
!   NORMAL APPROXIMATION FOR AMU.GT.AMAX
        ran = grndm()
        n=ran*sqrt(amu)+amu+half
        goto 999
      endif
      expma=exp(-amu)
      pir=one
      n=-1
 10   n=n+1
      pir=pir*frndm()
      if(pir.gt.expma) go to 10
 999  end
      subroutine ttbb(track,ktrack,parvec)

      implicit none
!----------------------------------------------------------------------*
! purpose:                                                             *
!   track a set of particle through a beam-beam interaction region.    *
!   see mad physicist's manual for the formulas used.                  *
!input:                                                                *
! input/output:                                                        *
!   track(6,*)(double)  track coordinates: (x, px, y, py, t, pt).      *
!   ktrack    (integer) number of tracks.                              *
!----------------------------------------------------------------------*
      include 'bb.fi'
      logical bborbit
      integer ktrack,itrack,ipos,get_option
      double precision track(6,*),parvec(*),pi,sx,sy,xm,ym,sx2,sy2,xs,  &
     &ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,cbx,cby,      &
     &get_variable,node_value,zero,one,two,three,ten3m,explim
      parameter(zero=0d0,one=1d0,two=2d0,three=3d0,ten3m=1d-3,          &
     &explim=150d0)
!     if x > explim, exp(-x) is outside machine limits.

!---- initialize.
      bborbit = get_option('bborbit ') .ne. 0
      pi=get_variable('pi ')
      sx = node_value('sigx ')
      sy = node_value('sigy ')
      xm = node_value('xma ')
      ym = node_value('yma ')
      fk = two * parvec(5) * parvec(6) / parvec(7)
      if (fk .eq. zero)  return
      ipos = 0
      if (.not. bborbit)  then
!--- find position of closed orbit bb_kick
        do ipos = 1, bbd_cnt
          if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
        enddo
        ipos = 0
    1   continue
      endif
      sx2 = sx*sx
      sy2 = sy*sy
!---- limit formulae for sigma(x) = sigma(y).
      if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
        do itrack = 1, ktrack
          xs = track(1,itrack) - xm
          ys = track(3,itrack) - ym
          rho2 = xs * xs + ys * ys
          tk = rho2 / (two * sx2)
          if (tk .gt. explim) then
            phix = xs * fk / rho2
            phiy = ys * fk / rho2
          else if (rho2 .ne. zero) then
            phix = xs * fk / rho2 * (one - exp(-tk) )
            phiy = ys * fk / rho2 * (one - exp(-tk) )
          else
            phix = zero
            phiy = zero
          endif
          if (ipos .ne. 0)  then
!--- subtract closed orbit kick
            phix = phix - bb_kick(1,ipos)
            phiy = phiy - bb_kick(2,ipos)
          endif
          track(2,itrack) = track(2,itrack) + phix
          track(4,itrack) = track(4,itrack) + phiy
        enddo

!---- case sigma(x) > sigma(y).
      else if (sx2 .gt. sy2) then
        r2 = two * (sx2 - sy2)
        r  = sqrt(r2)
        rk = fk * sqrt(pi) / r
        do itrack = 1, ktrack
          xs = track(1,itrack) - xm
          ys = track(3,itrack) - ym
          xr = abs(xs) / r
          yr = abs(ys) / r
          call ccperrf(xr, yr, crx, cry)
          tk = (xs * xs / sx2 + ys * ys / sy2) / two
          if (tk .gt. explim) then
            phix = rk * cry
            phiy = rk * crx
          else
            xb = (sy / sx) * xr
            yb = (sx / sy) * yr
            call ccperrf(xb, yb, cbx, cby)
            phix = rk * (cry - exp(-tk) * cby)
            phiy = rk * (crx - exp(-tk) * cbx)
          endif
          track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
          track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
          if (ipos .ne. 0)  then
!--- subtract closed orbit kick
            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
          endif
        enddo

!---- case sigma(x) < sigma(y).
      else
        r2 = two * (sy2 - sx2)
        r  = sqrt(r2)
        rk = fk * sqrt(pi) / r
        do itrack = 1, ktrack
          xs = track(1,itrack) - xm
          ys = track(3,itrack) - ym
          xr = abs(xs) / r
          yr = abs(ys) / r
          call ccperrf(yr, xr, cry, crx)
          tk = (xs * xs / sx2 + ys * ys / sy2) / two
          if (tk .gt. explim) then
            phix = rk * cry
            phiy = rk * crx
          else
            xb  = (sy / sx) * xr
            yb  = (sx / sy) * yr
            call ccperrf(yb, xb, cby, cbx)
            phix = rk * (cry - exp(-tk) * cby)
            phiy = rk * (crx - exp(-tk) * cbx)
          endif
          track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
          track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
          if (ipos .ne. 0)  then
!--- subtract closed orbit kick
            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
          endif
        enddo
      endif
      end

      subroutine trkill(n, turn, sum, jmax, part_id,                    &
     &last_turn, last_pos, last_orbit, z,aptype)

!hbu--- kill particle:  print, modify part_id list
      implicit none
      include 'name_len.fi'
      integer i,j,n,turn,part_id(*),jmax,last_turn(*)
      double precision sum, z(6,*), last_pos(*), last_orbit(6,*)
      character*(name_len) aptype
!hbu
      character*(name_len) el_name

      last_turn(part_id(n)) = turn
      last_pos(part_id(n)) = sum
      do j = 1, 6
        last_orbit(j,part_id(n)) = z(j,n)
      enddo

!hbu
      call element_name(el_name,len(el_name))
!hbu
      write(6,'(''particle #'',i6,'' lost turn '',i6,''  at pos. s ='', &
     &f10.2,'' element='',a,'' aperture ='',a)')                        &
     &part_id(n),turn,sum,el_name,aptype
      print *,"   X=",z(1,n),"  Y=",z(3,n),"  T=",z(5,n)

      do i = n+1, jmax
        part_id(i-1) = part_id(i)
        do j = 1, 6
          z(j,i-1) = z(j,i)
        enddo
      enddo
      jmax = jmax - 1
      end

      subroutine tt_putone(npart,turn,tot_segm,segment,part_id,z,orbit0,&
     &spos,ielem,el_name)

!hbu added spos, ielem, el_name
      implicit none
!----------------------------------------------------------------------*
!--- purpose: enter all particle coordinates in one table              *
!    input:                                                            *
!    npart  (int)           number of particles                        *
!    turn   (int)           turn number                                *
!    tot_segm (int)         total (target) number of entries           *
!    segment(int)           current segment count                      *
!    part_id (int array)    particle identifiers                       *
!    z (double (6,*))       particle orbits                            *
!    orbit0 (double array)  reference orbit                            *
!----------------------------------------------------------------------*
      include 'name_len.fi'
      integer i,j,npart,turn,tot_segm,segment,part_id(*),length
      double precision z(6,*),orbit0(6),tmp,tt,ss
!hbu was *36 allow longer info
      character*80 table,comment
!hbu
      integer ielem
!hbu name of element
      character*(name_len) el_name
!hbu
      double precision spos
!hbu
      character*4 vec_names(7)
!hbu
      data vec_names / 'x', 'px', 'y', 'py', 't', 'pt','s' /
      data table / 'trackone' /

!hbu
      length = len(comment)
      segment = segment + 1
!hbu
      write(comment, '(''#segment'',4i8,1X,A)')                         &
     &segment,tot_segm,npart,ielem,el_name
      call comment_to_table(table, comment, length)
      tt = turn
      do i = 1, npart
        call double_to_table(table, 'turn ', tt)
        ss = part_id(i)
        call double_to_table(table, 'number ', ss)
        do j = 1, 6
          tmp = z(j,i) - orbit0(j)
          call double_to_table(table, vec_names(j), tmp)
        enddo
!hbu spos
        call double_to_table(table,vec_names(7),spos)
        call augment_count(table)
      enddo
      end
      subroutine tt_puttab(npart,turn,nobs,orbit,orbit0,spos)

!hbu added spos
      implicit none
!----------------------------------------------------------------------*
!--- purpose: enter particle coordinates in table                      *
!    input:                                                            *
!    npart  (int)           particle number                            *
!    turn   (int)           turn number                                *
!    nobs   (int)           observation point number                   *
!    orbit  (double array)  particle orbit                             *
!    orbit0 (double array)  reference orbit                            *
!----------------------------------------------------------------------*
      integer npart,turn,j,nobs
      double precision orbit(6),orbit0(6),tmp,tt,tn
      double precision energy, get_value
      character*36 table
!hbu
      double precision spos
!hbu
      character*4 vec_names(7)
!hbu
      data vec_names / 'x', 'px', 'y', 'py', 't', 'pt', 's' /
      data table / 'track.obs$$$$.p$$$$' /

      tt = turn 
      tn = npart
      
      write(table(10:13), '(i4.4)') nobs
      write(table(16:19), '(i4.4)') npart
      
      energy = get_value('probe ','energy ')
      
      call double_to_table(table, 'turn ', tt)   ! the number of the current turn
      call double_to_table(table, 'number ', tn) ! the number of the current particle
      call double_to_table(table, 'e ', energy)
      
      do j = 1, 6
        tmp = orbit(j) - orbit0(j)
        call double_to_table(table, vec_names(j), tmp)
      enddo
!hbu spos
      call double_to_table(table,vec_names(7),spos)
      call augment_count(table)
      end
      subroutine trcoll(flag, apx, apy, turn, sum, part_id, last_turn,  &
     &last_pos, last_orbit, z, ntrk,al_errors)

      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   test for collimator aperture limits.                               *
! input:                                                               *
!   flag      (integer)   aperture type flag:                          *
!                         1: elliptic, 2: rectangular                  *
!   apx       (double)    x aperture or half axis                      *
!   apy       (double)    y aperture or half axis                      *
!   turn      (integer)   current turn number.                         *
!   sum       (double)    accumulated length.                          *
! input/output:                                                        *
!   part_id   (int array) particle identification list                 *
!   last_turn (int array) storage for number of last turn              *
!   last_pos  (dp. array) storage for last position (= sum)            *
!   last_orbit(dp. array) storage for last orbit                       *
!   z(6,*)    (double)    track coordinates: (x, px, y, py, t, pt).    *
!   ntrk      (integer) number of surviving tracks.                    *
!----------------------------------------------------------------------*
      include 'twiss0.fi'
      include 'name_len.fi'
      integer flag,turn,part_id(*),last_turn(*),ntrk,i,n,nn
      double precision apx,apy,sum,last_pos(*),last_orbit(6,*),z(6,*),  &
     &one,al_errors(align_max)
      parameter(one=1d0)
      character*(name_len) aptype

      n = 1
 10   continue
      do i = n, ntrk
!---- Is particle outside aperture?
        if (flag .eq. 1.and.((z(1,i)-al_errors(11))/apx)**2             &
     &+((z(3,i)-al_errors(12)) / apy)**2 .gt. one) then
          go to 99
        else if(flag .eq. 2                                             &
     &.and. (abs(z(1,i)-al_errors(11)) .gt. apx                         &
     &.or. abs(z(3,i)-al_errors(12)) .gt. apy)) then
          go to 99
!***  Introduction of marguerite : two ellipses
        else if(flag .eq. 3.and. ((z(1,i)-al_errors(11)) / apx)**2      &
     &+ ((z(3,i)-al_errors(12)) / apy)**2 .gt. one .and.                &
     &((z(1,i)-al_errors(11)) / apy)**2 +                               &
     &((z(3,i)-al_errors(12)) / apx)**2 .gt. one) then
          go to 99
        endif
        go to 98
 99     n = i
        nn=name_len
        call node_string('apertype ',aptype,nn)
        call trkill(n, turn, sum, ntrk, part_id,                        &
     &last_turn, last_pos, last_orbit, z,aptype)
        goto 10
 98     continue
      enddo
      end

      subroutine trinicmd(switch,orbit0,eigen,jend,z,turns,coords)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Define initial conditions for all particles to be tracked          *
! input:                                                               *
!   switch (int)  1: run, 2: dynap fastune, 3: dynap aperture          *
!   orbit0(6) - closed orbit                                           *
!   x, px, y, py, t, deltap, fx, phix, fy, phiy, ft, phit              *
!             - raw coordinates from start list                        *
!   eigen     - Eigenvectors                                           *
! output:                                                              *
!   jend      - number of particles to track                           *
!   z(6,jend) - Transformed cartesian coordinates incl. c.o.           *
!   coords      dp(6,0:turns,npart) (only switch > 1) particle coords. *
!----------------------------------------------------------------------*
      include 'track.fi'
      include 'bb.fi'
      logical zgiv,zngiv
      integer j,jend,k,kp,kq,next_start,itype(23),switch,turns,jt
      double precision phi,track(12),zstart(12),twopi,z(6,1000),zn(6),  &
     &ex,ey,et,orbit0(6),eigen(6,6),x,px,y,py,t,deltae,fx,phix,fy,phiy, &
     &ft,phit,get_value,get_variable,zero,deltax,coords(6,0:turns,*)
      double precision deltap,one
      parameter(zero=0d0,one=1d0)
      character*120 msg(2)

      deltap = get_variable('track_deltap ')
!---- Initialise orbit, emittances and eigenvectors etc.
      j = 0
      twopi = get_variable('twopi ')
      ex = get_value('probe ','ex ')
      ey = get_value('probe ','ey ')
      et = get_value('probe ','et ')
!-----get x add-on for lyaponuv calculation from dynap table
        deltax = get_value('dynap ', 'lyapunov ')
 1    continue
      jt  =  next_start(x,px,y,py,t,deltae,fx,phix,fy,phiy,ft,phit)
      if (switch.gt.1) then
      j = 2*jt-1
      else
      j = jt
      endif
      if (j .ne. -1 .and. j.ne.0)  then
        if (switch .lt. 3 .or. j .eq. 1)  then
          jend = j
          if (switch.gt.1) jend=jend+1 
!---- Get start coordinates
          track(1) = x
          track(2) = px
          track(3) = y
          track(4) = py
          track(5) = t
!          track(6) = deltae
!---- Here we absorb deltap into the PT variable
          track(6) = deltae + sqrt((one+deltap)**2+bet0i**2-one)-bet0i
          track(7) = fx
          track(8) = phix
          track(9) = fy
          track(10) = phiy
          track(11) = ft
          track(12) = phit
          do k = 1,12
            if(abs(track(k)).ne.zero) then
              itype(k) = 1
            else
              itype(k) = 0
            endif
          enddo
!---- Normalized coordinates.
          do kq = 1, 5, 2
            kp = kq + 1
            phi = twopi * track(kq+7)
            zn(kq) =   track(kq+6) * cos(phi)
            zn(kp) = - track(kq+6) * sin(phi)
          enddo
!---- Transform to unnormalized coordinates and refer to closed orbit.
          zgiv = .false.
          zngiv = .false.
          do k = 1, 6
            if (itype(k) .ne. 0) zgiv = .true.
            if (itype(k+6) .ne. 0) zngiv = .true.
            zstart(k) = track(k)                                        &
     &+ sqrt(ex) * (eigen(k,1) * zn(1) + eigen(k,2) * zn(2))            &
     &+ sqrt(ey) * (eigen(k,3) * zn(3) + eigen(k,4) * zn(4))            &
     &+ sqrt(et) * (eigen(k,5) * zn(5) + eigen(k,6) * zn(6))
          enddo
          if (switch .gt. 1)  then
!--- keep initial coordinates for dynap
            do k = 1, 6
              coords(k,0,j) = zstart(k)
            enddo
          endif
!---- Warn user about possible data conflict.
          if (zgiv .and. zngiv) then
            msg(1) = 'Absolute and normalized coordinates given,'
            msg(2) = 'Superposition used.'
            call aawarn('START-1: ', msg(1))
            call aawarn('START-2: ', msg(2))
          endif
          do k = 1, 6
            z(k,j) = orbit0(k) + zstart(k)
          enddo
          if (switch.gt.1) then
          do k = 1, 6
            z(k,j+1) = z(k,j)
          enddo
            z(1,j+1) = z(1,j) + deltax
          endif
        endif
        goto 1
      endif
c      if (switch .eq. 3)  then
!--- create second particle with x add-on
c        deltax = get_value('dynap ', 'lyapunov ')
c        jend = 2
c        z(1,jend) = z(1,1) + deltax
c        do k = 2, 6
c          z(k,jend) = z(k,1)
c        enddo
c      endif
      end
      subroutine trbegn(rt,eigen)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Initialize tracking mode; TRACK command execution.                 *
!----------------------------------------------------------------------*
      logical m66sta,onepass
      integer get_option
      double precision rt(6,6),reval(6),aival(6),eigen(6,6),zero
      parameter(zero=0d0)

!---- Initialize
      onepass = get_option('onepass ') .ne. 0
!---- One-pass system: Do not normalize.
      if (onepass) then
        call m66one(eigen)
      else
!---- Find eigenvectors.
        if (m66sta(rt)) then
          call laseig(rt,reval,aival,eigen)
        else
          call ladeig(rt,reval,aival,eigen)
        endif
      endif
      end
      subroutine ttdpdg(track, ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose: computes the effect of the dipole edge                      *
! Input/Output:  ktrack , number of surviving trajectories             *
!                 track , trajectory coordinates
!----------------------------------------------------------------------*
      integer ktrack,itrack
      double precision fint,bv0,e1,h,hgap,corr,rw(6,6),tw(6,6,6),       &
     &track(6,*),node_value,zero,one
      parameter(zero=0d0,one=1d0)

      call m66one(rw)
!      call dzero(tw, 216)

      bv0 = node_value('dipole_bv ')
      e1 = bv0 * node_value('e1 ')
!--- bvk applied further down
      h = node_value('h ')
      hgap = node_value('hgap ')
      fint = node_value('fint ')
      corr = (h + h) * hgap * fint
!          print*,"------------------------------------------ "
!---- Fringe fields effects computed from the TWISS routine tmfrng
!     tmfrng returns the matrix elements rw(used) and tw(unused)
!     No radiation effects as it is a pure thin lens with no lrad
      call tmfrng(.false.,h,zero,e1,zero,zero,corr,rw,tw)
      do itrack = 1, ktrack
        track(2,itrack) = track(2,itrack) + rw(2,1)*track(1,itrack)
        track(4,itrack) = track(4,itrack) + rw(4,3)*track(3,itrack)
      enddo
      return
      end

      subroutine trsol(track,ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Track a set of particles through a thin solenoid.                  *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!   KTRACK    (integer)   number of surviving tracks.                  *
! Output:                                                              *
!   EL        (double)    Length of solenoid.                          *
!----------------------------------------------------------------------*
      integer itrack,ktrack
      double precision dtbyds,gamma,beta,deltap
      double precision get_value, node_value
      double precision el,track(6,*),one,two,elrad
      double precision sk,skl,sks,sksl,cosTh,sinTh,Q,R,Z
      double precision xf,yf,pxf,pyf,sigf,psigf
      double precision onedp,fpsig,fppsig
      parameter(one=1d0,two=2d0)
!
!---- Initialize.
!      dtbyds  = get_value('probe ','dtbyds ')
!      gamma   = get_value('probe ','gamma ')
      beta    = get_value('probe ','beta ')
!      deltap  = get_value('probe ','deltap ')
!
!---- Get solenoid parameters
!      elrad   = node_value('lrad ')
      sksl    = node_value('ksl ')
      sks     = node_value('ks ')
!
!---- Set up strengths
!      sk    = sks / two / (one + deltap)
      sk    = sks / two
      skl   = sksl / two
!
!---- Loop over particles
!
      do  itrack = 1, ktrack
!
!     Ripken formulae p.28 (3.35 and 3.36)
        xf   = track(1,itrack)
        yf   = track(3,itrack)
        psigf= track(6,itrack)/beta
!
!     We do not use a constant deltap!!!!! WE use full 6D formulae!
        onedp   = sqrt(one+2*psigf+(beta**2)*(psigf**2))
        fpsig   = onedp - one
        fppsig  = (one+(beta**2)*psigf)/onedp
!
!     Set up C,S, Q,R,Z
        cosTh = cos(skl/onedp)
        sinTh = sin(skl/onedp)
        Q = -skl*sk/onedp
        R = fppsig/(onedp**2)*skl*sk
        Z = fppsig/(onedp**2)*skl
!
!
        pxf  = track(2,itrack) + xf*Q
        pyf  = track(4,itrack) + yf*Q
        sigf = track(5,itrack)*beta - 0.5D0*(xf**2 + yf**2)*R

!       Ripken formulae p.29 (3.37)
!
        track(1,itrack) =  xf*cosTh  + yf*sinTh
        track(2,itrack) =  pxf*cosTh + pyf*sinTh
        track(3,itrack) = -xf*sinTh  + yf*cosTh
        track(4,itrack) = -pxf*sinTh + pyf*cosTh
        track(5,itrack) =  (sigf + (xf*pyf - yf*pxf)*Z)/beta
!        track(6,itrack) =  psigf*beta

      enddo
      end

      subroutine tttrans(track,ktrack)
      implicit none
!----------------------------------------------------------------------*
! Purpose:                                                             *
!   Translation of a set of particles.                                 *
! Input/output:                                                        *
!   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
!----------------------------------------------------------------------*
      integer itrack,ktrack
      double precision get_value, node_value
      double precision track(6,*)
      double precision t_x,t_px,t_y,t_py,t_sig,t_psig
!
!---- Initialize.
!
!---- Get translation parameters
      t_x    = node_value('x ')
      t_px   = node_value('px ')
      t_y    = node_value('y ')
      t_py   = node_value('py ')
      t_sig  = node_value('t ')
      t_psig = node_value('pt ')
!
!---- Loop over particles
!
      do  itrack = 1, ktrack
!
!     Add vector to particle coordinates
        track(1,itrack) = track(1,itrack) + t_x       
        track(2,itrack) = track(2,itrack) + t_px      
        track(3,itrack) = track(3,itrack) + t_y      
        track(4,itrack) = track(4,itrack) + t_py      
        track(5,itrack) = track(5,itrack) + t_sig      
        track(6,itrack) = track(6,itrack) + t_psig      
!
      enddo
      end
