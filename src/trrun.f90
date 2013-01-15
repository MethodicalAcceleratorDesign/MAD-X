LOGICAL FUNCTION  is_drift()
  double precision node_value, code
  code = node_value('mad8_type ')
  is_drift = code.eq.1;
END FUNCTION is_drift

LOGICAL FUNCTION  is_matrix()
  double precision node_value, code
  code = node_value('mad8_type ')
  is_matrix = code.eq.4;
END FUNCTION is_matrix

LOGICAL FUNCTION  is_quad()
  double precision node_value, code
  code = node_value('mad8_type ')
  is_quad = code.eq.5;
END FUNCTION is_quad

LOGICAL FUNCTION  is_thin()
  double precision node_value, el, zero
  parameter(zero=0d0)
  el = node_value('l ')
  is_thin = el.eq.zero;
END FUNCTION is_thin

subroutine trrun(switch,turns,orbit0,rt,part_id,last_turn,        &
     last_pos,z,dxt,dyt,last_orbit,eigen,coords,e_flag,code_buf,l_buf)

  use bbfi
  use twiss0fi
  use name_lenfi
  use trackfi
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
  logical onepass,onetable,last_out,info,aperflag,doupdate
  integer j,code,restart_sequ,advance_node,node_al_errors,n_align,       &
       nlm,jmax,j_tot,turn,turns,i,k,get_option,ffile,SWITCH,nint,ndble, &
       nchar,part_id(*),last_turn(*),char_l,segment, e_flag, nobs,lobs,  &
       int_arr(1),tot_segm,code_buf(*)
  double precision tmp_d,orbit0(6),orbit(6),el,re(6,6),rt(6,6),          &
       al_errors(align_max),z(6,*),zz(6),dxt(*),dyt(*),eigen(6,6),sum,   &
       node_value,one,                                                   &
       get_variable,last_pos(*),last_orbit(6,*),maxaper(6),get_value,    &
       zero,obs_orb(6),coords(6,0:turns,*),l_buf(*),deltap
  parameter(zero=0d0,one=1d0)
  character(12) tol_a, char_a
  !hbu
  double precision spos
  !hbu
  character(4) vec_names(7)
  !hbu
  character(name_len) el_name
  character(120) msg
  data tol_a,char_a / 'maxaper ', ' ' /
  !hbu
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt','s' /

  logical is_drift, is_thin, is_quad, is_matrix

  !-------added by Yipeng SUN 01-12-2008--------------
  deltap = get_value('probe ','deltap ')

  if(deltap.eq.0) then
     onepass = get_option('onepass ') .ne. 0
     if(onepass)   then
     else
        call trclor(orbit0)
     endif
  else
  endif
  !-------added by Yipeng SUN 01-12-2008--------------

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
          'allowed for tracking.'
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
       char_a, char_l)
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
     call double_to_table_curr('tracksumm ', 'number ', tmp_d)
     tmp_d = 0
     call double_to_table_curr('tracksumm ', 'turn ', tmp_d)
     do j = 1, 6
        tmp_d = z(j,i) - orbit0(j)
        call double_to_table_curr('tracksumm ', vec_names(j), tmp_d)
     enddo
     !hbu add s
     call double_to_table_curr('tracksumm ',vec_names(7),spos)
     call augment_count('tracksumm ')
  enddo
  !--- enter first turn, and possibly eigen in tables
  if (switch .eq. 1)  then
     if (onetable)  then
        call track_pteigen(eigen)
        !hbu add s, node id and name
        call tt_putone(jmax, 0, tot_segm, segment, part_id,           &
             z, orbit0,spos,nlm,el_name)
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

  doupdate = get_option('update ') .ne. 0

  !--- loop over turns
  nobs = 0
  do turn = 1, turns

     if (doupdate) call trupdate(turn)

     j = restart_sequ()
     nlm = 0
     sum=zero
     !--- loop over nodes
10   continue
     bbd_pos=j
     if (turn .eq. 1)  then
        code = node_value('mad8_type ')
        if(code.eq.39) code=15
        if(code.eq.38) code=24
        el = node_value('l ')
        code_buf(nlm+1) = code
        l_buf(nlm+1) = el
        !hbu get current node name
        call element_name(el_name,len(el_name))
        if (.not.(is_drift() .or. is_thin() .or. is_quad() .or. is_matrix())) then
           print*," "
           print*,"code: ",code," el: ",el,"   THICK ELEMENT FOUND"
           sum = node_value('name ')
           print*," "
           print*,"Track dies nicely"
           print*,"Thick lenses will get nowhere"
           print*,"MAKETHIN will save you"
           print*," "
           print*," "
           !           Better to use stop. If use return, irrelevant tracking output
           !           is printed
           stop
           !            call aafail('TRRUN',                                          &
           !           '----element with length found : CONVERT STRUCTURE WITH '//   &
           !           'MAKETHIN')
        endif
     else
        el = l_buf(nlm+1)
        code = code_buf(nlm+1)
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
              call dcopy(z(1,i),zz,6)
              call tmali1(zz,al_errors, betas, gammas,z(1,i), re)
           enddo
        endif
     endif
     !-------- Track through element  // suppress dxt 13.12.04
     call ttmap(code,el,z,jmax,dxt,dyt,sum,turn,part_id,             &
          last_turn,last_pos, last_orbit,aperflag,maxaper,al_errors,onepass)
     !--------  Misalignment at end of element (from twissfs.f)
     if (code .ne. 1)  then
        if (n_align .ne. 0)  then
           do i = 1, jmax
              call dcopy(z(1,i),zz,6)
              call tmali2(el,zz, al_errors, betas, gammas,z(1,i), re)
           enddo
        endif
     endif
     nlm = nlm+1
     if (nobs .gt. 0)  then
        call dzero(obs_orb,6)
        call get_node_vector('obs_orbit ', lobs, obs_orb)
        if (lobs .lt. 6)                                              &
             call aafail('TRACK', 'obs. point orbit not found')
        if (onetable)  then
           !hbu
           spos=sum
           !hbu get current node name
           call element_name(el_name,len(el_name))
           !hbu
           call tt_putone(jmax, turn, tot_segm, segment, part_id,      &
                z, obs_orb,spos,nlm,el_name)
        else
          if (mod(turn, ffile) .eq. 0)  then
            do i = 1, jmax
                !hbu add spos
                call tt_puttab(part_id(i), turn, nobs, z(1,i), obs_orb,   &
                    spos)
            enddo
           endif
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
                   z, orbit0,spos,nlm,el_name)
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
          goto 20
     if (switch .eq. 2 .and. info) then
        if (mod(turn,100) .eq. 0) print *, 'turn :', turn
     endif
  enddo
  !--- end of loop over turns
20 continue
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
                z, orbit0,spos,nlm,el_name)
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
     call double_to_table_curr('tracksumm ', 'number ', tmp_d)
     tmp_d = last_turn(i)
     call double_to_table_curr('tracksumm ', 'turn ', tmp_d)
     do j = 1, 6
        tmp_d = last_orbit(j,i) - orbit0(j)
        call double_to_table_curr('tracksumm ', vec_names(j), tmp_d)
     enddo
     !hbu
     spos=last_pos(i)
     !hbu
     call double_to_table_curr('tracksumm ',vec_names(7),spos)
     call augment_count('tracksumm ')
  enddo
end subroutine trrun


subroutine ttmap(code,el,track,ktrack,dxt,dyt,sum,turn,part_id,   &
     last_turn,last_pos, last_orbit,aperflag,maxaper,al_errors,        &
     onepass)

  use twtrrfi
  use twiss0fi
  use name_lenfi
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
  logical aperflag,fmap,onepass
  integer turn,code,ktrack,part_id(*),last_turn(*),nn,jtrk,         &
       get_option
  double precision apx,apy,apr,el,sum,node_value,track(6,*),        &
       last_pos(*),last_orbit(6,*),parvec(26),get_value,ct,tmp,          &
       aperture(maxnaper),one,maxaper(6), zero,al_errors(align_max),st,  &
       theta,ek(6),re(6,6),te(6,6,6),craporb(6),dxt(*),dyt(*),offset(2), &
       offx,offy
  character(name_len) aptype
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

  !---- Rotate trajectory before entry
  theta = node_value('tilt ')
  if (theta .ne. zero)  then
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

  !---- Translation of particles // AK 23.02.2006
  !---- useful for combination of different transfer lines or rings
  if(code.eq.36) then
     if (onepass) then
        call tttrans(track,ktrack)
        go to 500
     endif
  endif

  !---- Beam-beam,  standard 4D, no aperture
  if(code.eq.22) then
     parvec(5)=get_value('probe ', 'arad ')
     parvec(6)=node_value('charge ') * get_value('probe ', 'npart ')
     parvec(7)=get_value('probe ','gamma ')
     call ttbb(track, ktrack, parvec)
     go to 500
  endif

  !---- Special colllimator aperture check taken out AK 20071211
  !!---- Collimator with elliptic aperture.
  !      if(code.eq.20) then
  !        apx = node_value('xsize ')
  !        apy = node_value('ysize ')
  !        if(apx.eq.zero) then
  !          apx=maxaper(1)
  !        endif
  !        if(apy.eq.zero) then
  !          apy=maxaper(3)
  !        endif
  !        call trcoll(1, apx, apy, turn, sum, part_id, last_turn,         &
  !     last_pos, last_orbit, track, ktrack,al_errors)
  !        go to 500
  !      endif
  !!---- Collimator with rectangular aperture.
  !      if(code.eq.21) then
  !        apx = node_value('xsize ')
  !        apy = node_value('ysize ')
  !        if(apx.eq.zero) then
  !          apx=maxaper(1)
  !        endif
  !        if(apy.eq.zero) then
  !          apy=maxaper(3)
  !        endif
  !        call trcoll(2, apx, apy, turn, sum, part_id, last_turn,         &
  !     last_pos, last_orbit, track, ktrack,al_errors)
  !        go to 500
  !      endif


  !---- Test aperture. ALL ELEMENTS BUT DRIFTS
  aperflag = get_option('aperture ') .ne. 0
  if(aperflag) then
     nn=name_len
     call node_string('apertype ',aptype,nn)
     call dzero(aperture,maxnaper)
     call get_node_vector('aperture ',nn,aperture)
     call dzero(offset,2)
     call get_node_vector('aper_offset ',nn,offset)

     offx = offset(1)
     offy = offset(2)
     !        print *, " TYPE ",aptype &
     !       "values  x y lhc",aperture(1),aperture(2),aperture(3)
     !------------  ellipse case ----------------------------------
     if(aptype.eq.'ellipse') then
        apx = aperture(1)
        apy = aperture(2)
        call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
             last_pos, last_orbit, track, ktrack,al_errors,offx,offy)
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
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !------------  rectangle case ----------------------------------
     else if(aptype.eq.'rectangle') then
        apx = aperture(1)
        apy = aperture(2)
        call trcoll(2, apx, apy, turn, sum, part_id, last_turn,       &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !-------------Racetrack type , added by Yipeng SUN 21-10-2008---
     else if(aptype.eq.'racetrack') then
        apx = aperture(1)
        apy = aperture(2)
        apr = aperture(3)
        call trcoll1(4, apx, apy, turn, sum, part_id, last_turn,      &
             last_pos,last_orbit,track,ktrack,al_errors,apr,offx,offy)
        !------------  LHC screen case ----------------------------------
     else if(aptype.eq.'lhcscreen') then
        !        print *, "LHC screen start, Xrect= ",
        !       aperture(1),"  Yrect= ",aperture(2),"  Rcirc= ",aperture(3)
        apx = aperture(3)
        apy = aperture(3)
        !JMJ!     Making essential changes in AV's absence, 16/7/2003
        !JMJ!     this tests whether the particle is outside the circumscribing
        !JMJ!     circle.
        call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
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
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !        print *, "LHC screen end"
        !------------  marguerite case ----------------------------------
     else if(aptype.eq.'marguerite') then
        apx = aperture(1)
        apy = aperture(2)
        call trcoll(3, apx, apy, turn, sum, part_id, last_turn,       &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !------------  rectellipse case ----------------------------------
     else if(aptype.eq.'rectellipse') then
        !*****         test ellipse
        apx = aperture(3)
        apy = aperture(4)
        call trcoll(1, apx, apy, turn, sum, part_id, last_turn,       &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !*****         test rectangle
        apx = aperture(1)
        apy = aperture(2)
        call trcoll(2, apx, apy, turn, sum, part_id, last_turn,       &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !       print*, " test apertures"
        !       print*, "      apx=",apx, " apy=",apy," apxell=",apxell,
        !              " apyell=",apyell
        !          call trcoll(3, apx, apy, turn, sum, part_id, last_turn,       &
        !         last_pos, last_orbit, track,ktrack,al_errors)
     endif
     !      else
     !!---- Check on collimator xsize/ysize even if user did not use 'aperture'
     !!---- option in track command. This simulates roughly the old MAD-X
     !!---- behaviour.      	
     !!---- Collimator with elliptic aperture.
     !      if(code.eq.20) then
     !        apx = node_value('xsize ')
     !        apy = node_value('ysize ')
     !        if(apx.eq.zero) then
     !          apx=maxaper(1)
     !        endif
     !        if(apy.eq.zero) then
     !          apy=maxaper(3)
     !        endif
     !        call trcoll(1, apx, apy, turn, sum, part_id, last_turn,         &
     !     last_pos, last_orbit, track, ktrack,al_errors)
     !        go to 500
     !      endif
     !!---- Collimator with rectangular aperture.
     !      if(code.eq.21) then
     !        apx = node_value('xsize ')
     !        apy = node_value('ysize ')
     !        if(apx.eq.zero) then
     !          apx=maxaper(1)
     !        endif
     !        if(apy.eq.zero) then
     !          apy=maxaper(3)
     !        endif
     !        call trcoll(2, apx, apy, turn, sum, part_id, last_turn,         &
     !      last_pos, last_orbit, track, ktrack,al_errors)
     !        go to 500
     !      endif      	
  endif
  !----  END OF IF(APERFLAG)  ----------------


  !-- switch on element type BUT DRIFT, COLLIMATORS, BEAM_BEAM / 13.03.03
  !-- 500 has been specified at the relevant places in go to below
  !-- code =1 for drift, treated above, go to 500 directly
  !      print *,"   CODE    ",code
  go to ( 500,  20,  30,  40,  50,  60,  70,  80,  90, 100,         &
       110, 120, 130, 140, 150, 160, 170, 180, 190, 500,                 &
       500, 500, 230, 240, 250, 260, 270, 280, 290, 300,   310, 320,     &
  !     330, 500, 350, 360, 370,500,500,400,410,500, 500, 500, 500), code
  ! Use this line to enable non-linear thin lens
  !     330, 500, 350, 360, 370,500,500,400,410,420,500, 500, 500, 500), code
  ! Enable non-linear thin lens and RF-Multipole
       330, 500, 350, 360, 370,500,500,400,410,420,430, 500, 500, 500), code
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
  call tmarb(.false.,.false.,craporb,fmap,ek,re,te)
  call tttrak(ek,re,track,ktrack)
  go to 500
  !---- Quadrupole. OBSOLETE, to be kept for go to
  !---- AL: not so fast, if here it's a thick quadrupole
50 continue
  call tttquad(track,ktrack)
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
  call ttmult(track,ktrack,dxt,dyt,turn)
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
  call ttsep(el, track, ktrack)
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
  call ttcorr(el, track, ktrack, turn)
  go to 500
  !---- ECollimator, RCollimator, BeamBeam, Lump. ??
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
  !---- Crab cavity.
370 continue
  call ttcrabrf(track,ktrack,turn)
  go to 500
  !---- h ac dipole.
400 continue
  call tthacdip(track,ktrack,turn)
  go to 500
  !---- v ac dipole.
410 continue
  call ttvacdip(track,ktrack,turn)
  go to 500
  !---- nonlinear elliptical lens
420 continue
  call ttnllens(track,ktrack)
  go to 500
  !---- rf multipoles
430 continue
  call ttrfmult(track,ktrack,turn)
  go to 500

  !---- Solenoid.


500 continue

  !---- Rotate trajectory at exit
  if (theta .ne. zero)  then
     do jtrk = 1,ktrack
        tmp = track(1,jtrk)
        track(1,jtrk) = ct * tmp - st * track(3,jtrk)
        track(3,jtrk) = ct * track(3,jtrk) + st * tmp
        tmp = track(2,jtrk)
        track(2,jtrk) = ct * tmp - st * track(4,jtrk)
        track(4,jtrk) = ct * track(4,jtrk) + st * tmp
     enddo
  endif

  !---- Accumulate length.
502 sum = sum + el
  return
end subroutine ttmap

subroutine ttmult(track, ktrack,dxt,dyt,turn)

  use twtrrfi
  use name_lenfi
  use trackfi
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
  logical first
  integer iord,jtrk,nd,nord,ktrack,j,n_ferr,nn,ns,node_fd_errors,   &
       get_option,turn,noisemax,nn1,in
  double precision     const,curv,dbi,dbr,dipi,dipr,dx,dy,elrad,    &
       pt,px,py,rfac,rpt1,rpt2,rpx1,rpx2,rpy1,rpy2,                      &
       f_errors(0:maxferr),field(2,0:maxmul),vals(2,0:maxmul),           &
       ordinv(maxmul),track(6,*),dxt(*),dyt(*),normal(0:maxmul),         &
       skew(0:maxmul),bvk,node_value,zero,one,two,three,half,ttt,        &
       npeak(100), nlag(100), ntune(100), temp,noise

  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,half=5d-1)
  character(name_len) aptype

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
  bvk = node_value('other_bv ')
  !---- Multipole length for radiation.
  elrad = node_value('lrad ')
  noise = node_value('noise ')
  !---- Multipole components.
  call dzero(normal,maxmul+1)
  call dzero(skew,maxmul+1)
  call get_node_vector('knl ',nn,normal)
  call get_node_vector('ksl ',ns,skew)
  nd = 2 * max(nn, ns)
  call dzero(vals,2*(maxmul+1))

  if(noise .eq. 1)   then
     nn1=name_len
     noisemax = node_value('noisemax ')
     call dzero(npeak,noisemax)
     call dzero(ntune,noisemax)
     call dzero(nlag,noisemax)
     call get_node_vector('npeak ',nn1,npeak)
     call get_node_vector('ntune ',nn1,ntune)
     call get_node_vector('nlag ',nn1,nlag)

     temp = 0
     do in = 1, noisemax
        temp = temp + npeak(in) * sin(nlag(in) + ntune(in) * turn)
     enddo


     !   temp = npeak * sin(nlag + ntune * turn)
     do iord = 0, nn
        vals(1,iord) = normal(iord) * (1+temp)
     enddo
     do iord = 0, ns
        vals(2,iord) = skew(iord) * (1+temp)
     enddo
  else
     do iord = 0, nn
        vals(1,iord) = normal(iord)
     enddo
     do iord = 0, ns
        vals(2,iord) = skew(iord)
     enddo
  endif

  !  do iord = 0, nn
  !     vals(1,iord) = normal(iord)
  !  enddo
  !  do iord = 0, ns
  !     vals(2,iord) = skew(iord)
  !  enddo
  !---- Field error vals.
  call dzero(field,2*(maxmul+1))
  if (n_ferr .gt. 0) then
     call dcopy(f_errors,field,n_ferr)
  endif
  !-----added FrankS, 10-12-2008
  nd = 2 * max(nn, ns, n_ferr/2-1)
  !---- Dipole error.
  !      dbr = bvk * field(1,0) / (one + deltas)
  !      dbi = bvk * field(2,0) / (one + deltas)
  dbr = bvk * field(1,0)
  dbi = bvk * field(2,0)
  !---- Nominal dipole strength.
  !      dipr = bvk * vals(1,0) / (one + deltas)
  !      dipi = bvk * vals(2,0) / (one + deltas)
  dipr = bvk * vals(1,0)
  dipi = bvk * vals(2,0)
  !---- Other components and errors.
  nord = 0
  do iord = 1, nd/2
     do j = 1, 2
        !          field(j,iord) = bvk * (vals(j,iord) + field(j,iord))          &
        !     / (one + deltas)
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
             field(1,nord)*track(1,jtrk) - field(2,nord)*track(3,jtrk)
        dyt(jtrk) =                                                   &
             field(1,nord)*track(3,jtrk) + field(2,nord)*track(1,jtrk)
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
                (dipi + dyt(jtrk))**2) / elrad

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
     !     (dbr + dxt(jtrk) - dipr * (deltas + beti*track(6,jtrk)))
     !        track(4,jtrk) = track(4,jtrk) +                                 &
     !     (dbi + dyt(jtrk) - dipi * (deltas + beti*track(6,jtrk)))
     !        track(5,jtrk) = track(5,jtrk)                                   &
     !     - (dipr*track(1,jtrk) - dipi*track(3,jtrk)) * beti
     track(2,jtrk) = track(2,jtrk) -                                 &
          (dbr + dxt(jtrk) - dipr * (ttt - one))
     track(4,jtrk) = track(4,jtrk) +                                 &
          (dbi + dyt(jtrk) - dipi * (ttt - one))
     track(5,jtrk) = track(5,jtrk) -                                 &
          (dipr*track(1,jtrk) - dipi*track(3,jtrk)) *                       &
          ((one + bet0*track(6,jtrk))/ttt)*bet0i
  enddo

  !---- Radiation loss at exit.
  if (dorad .and. elrad .ne. 0) then

     !---- Full damping.
     if (dodamp) then
        do jtrk = 1,ktrack
           curv = sqrt((dipr + dxt(jtrk))**2 +                         &
                (dipi + dyt(jtrk))**2) / elrad

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

end subroutine ttmult

subroutine ttsrot(track,ktrack)

  use trackfi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Coordinate change due to rotation about s axis                     *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.        tmsrot      *
  !----------------------------------------------------------------------*
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
end subroutine ttsrot


subroutine ttyrot(track,ktrack)

  use trackfi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Coordinate change due to rotation about y axis                     *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.        tmsrot      *
  !----------------------------------------------------------------------*
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
end subroutine ttyrot

subroutine ttdrf(el,track,ktrack)

  use trackfi
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
     !     + el*(beti + pt * dtbyds) - (beti+pt)*ttt
     !---- AK 20060413
     !---- Ripken DESY-95-189 p.36
     track(5,itrack) = track(5,itrack)                               &
          + bet0i*(el - (one + bet0*pt) * ttt)
  enddo
end subroutine ttdrf
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
  double precision omega,phirf,pt,rff,rfl,rfv,track(6,*),clight,    &
       twopi,vrf,pc0,get_variable,node_value,get_value,one,two,half,bvk, &
       ten3m,ten6p
  !      double precision px,py,ttt,beti,el1
  parameter(one=1d0,two=2d0,half=5d-1,ten3m=1d-3,ten6p=1d6)

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- BV flag
  bvk = node_value('other_bv ')

  !---- Fetch data.
  !      el = node_value('l ')
  !      el1 = node_value('l ')
  rfv = bvk * node_value('volt ')
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
end subroutine ttrf
subroutine ttcrabrf(track,ktrack,turn)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of trajectories through a thin crab cavity (zero l.)   *
  !   The crab cavity is sandwiched between 2 drift spaces half length.  *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.                    *
  ! Output:                                                              *
  !   EL        (double)    Length of quadrupole.                        *
  !----------------------------------------------------------------------*
  ! Added: R. Calaga/F. Zimmermann (10/06, 11/07, 09/10)                 *
  !----------------------------------------------------------------------*
  integer itrack,ktrack,turn,t1,t2,t3,t4,p1,p2
  double precision omega,phirf,pt,rff,rfl,rfv,eph,track(6,*),clight,twopi,vrf,pc0,  &
       get_variable,node_value,get_value,one,two,half,ten3m,ten6p,px
  !      double precision px,py,ttt,beti,el1
  parameter(one=1d0,two=2d0,half=5d-1,ten3m=1d-3,ten6p=1d6)

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- Fetch data.
  rfv = node_value('volt ')
  rff = node_value('freq ')
  rfl = node_value('lag ')
  pc0 = get_value('beam ','pc ')

  !--- specify voltage ramp (relative)
  t1 = node_value('rv1 ')
  t2 = node_value('rv2 ') + t1
  t3 = node_value('rv3 ') + t2
  t4 = node_value('rv4 ') + t3

  !--- specify phase change in two steps
  p1 = node_value('rph1 ')
  p2 = node_value('rph2 ') + p1
  eph = node_value('lagf ')

  !---- Set up voltage ramp.
  omega = rff * (ten6p * twopi / clight)
  
  if (turn .lt. t1)  then
     vrf = 0.0
  else if (turn .ge. t1 .and. turn .lt. t2) then
     vrf = (turn-t1)*rfv*ten3m/pc0/(t2-t1)
  else if (turn .ge. t2 .and. turn .lt. t3) then
     vrf = rfv*ten3m/pc0
  else if (turn .ge. t3 .and. turn .lt. t4) then
     vrf = (t4-turn)*rfv*ten3m/pc0/(t4-t3)
  else
     vrf = 0.0
  endif

  !---- Set up phase ramp.
  if (turn .lt. p1)  then
     phirf = rfl*twopi
  else if (turn .ge. p1 .and. turn .lt. p2) then
     phirf = (turn-p1)*eph*twopi/(p2-p1)
  else
     phirf= eph*twopi
  endif

  !  print*," turn: ",turn, " phase: ", phirf*360/twopi


  do itrack = 1, ktrack
     px  = track(2,itrack)                                           &
          + vrf * sin(phirf - omega * track(5,itrack))
     pt = track(6,itrack)                                            &
          - omega* vrf * track(1,itrack) *                                  &
          cos(phirf - omega * track(5,itrack))

     !---- track(2,jtrk) = track(2,jtrk)
     !        pt = track(6,itrack)
     !        pt = pt + vrf * sin(phirf - omega * track(5,itrack))

     track(2,itrack) = px
     track(6,itrack) = pt
  enddo
end subroutine ttcrabrf
subroutine ttsep(el,track,ktrack)

  use twtrrfi
  use trackfi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  ! ttsep tracks particle through an electrostatic separator with an     *
  ! homogeneous, purely perpendicularly acting electric field. The       *
  ! element is given as a thick element and then sliced by makethin.     *
  ! Input:                                                               *
  !   EL         (double)     Length of the element                      *
  !   KTRACK:    (integer)    Number of surviving tracks                 *
  ! Input/Output:                                                        *
  !   TRACK(6,*) (double)     Track coordinantes: (X, PX, Y, PY, T, PT)  *
  !----------------------------------------------------------------------*
  integer ktrack,itrack
  double precision track(6,*),node_value,get_value
  double precision el,ex,ey,tilt,charge,mass,pc,beta,deltap,kick0
  double precision cos_tilt,sin_tilt,efieldx,efieldy,pt
  double precision one,ten3m
  parameter(one=1.0d0,ten3m=1.0d-3)
  !
  ex=node_value('ex_l ')
  ey=node_value('ey_l ')
  tilt=node_value('tilt ')
  cos_tilt=cos(tilt)
  sin_tilt=sin(tilt)

  charge=get_value('probe ','charge ')
  mass  =get_value('probe ','mass ')
  pc    =get_value('probe ','pc ')
  beta  =get_value('probe ','beta ')
  efieldx=ex*cos_tilt + ey*sin_tilt
  efieldy=-ex*sin_tilt + ey*cos_tilt
  do itrack = 1, ktrack
     pt=track(6,itrack)
     deltap=sqrt(one-one/beta/beta+(pt+one/beta)**2) - one
     kick0=charge*ten3m/pc/(one+deltap)/beta
     track(2,itrack) = track(2,itrack) + kick0*efieldx
     track(4,itrack) = track(4,itrack) + kick0*efieldy
  end do

end subroutine ttsep
subroutine ttcorr(el,track,ktrack,turn)

  use twtrrfi
  use trackfi
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
  !      logical dorad,dodamp,dorand
  integer itrack,ktrack,n_ferr,node_fd_errors,code,bvk,i,get_option
  integer sinkick,turn
  double precision bi2gi2,bil2,curv,dpx,dpy,el,pt,px,py,rfac,rpt,   &
       rpx,rpy,track(6,*),xkick,ykick,                                   &
       div,f_errors(0:maxferr),field(2),get_variable,get_value,          &
       node_value,zero,one,two,three,half,twopi
  double precision dpxx,dpyy
  double precision temp,sinpeak,sintune,sinphase
  parameter(zero=0d0,one=1d0,two=2d0,three=3d0,half=5d-1)

  !---- Initialize.
  twopi=get_variable('twopi ')
  bvk = node_value('other_bv ')
  deltas = get_variable('track_deltap ')
  arad = get_value('probe ','arad ')
  betas = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')
  dtbyds = get_value('probe ','dtbyds ')
  code = node_value('mad8_type ')
  if(code.eq.39) code=15
  if(code.eq.38) code=24
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
     xkick=bvk*(node_value('kick ')+node_value('chkick ')+           &
          field(1)/div)
     ykick=zero
  else if (code.eq.15) then
     xkick=bvk*(node_value('hkick ')+node_value('chkick ')+          &
          field(1)/div)
     ykick=bvk*(node_value('vkick ')+node_value('cvkick ')+          &
          field(2)/div)
  else if(code.eq.16) then
     xkick=zero
     ykick=bvk*(node_value('kick ')+node_value('cvkick ')+           &
          field(2)/div)
  else
     xkick=zero
     ykick=zero
  endif

  !---- Sinusoidal kick
  sinkick = node_value('sinkick ')
  if(sinkick.eq.1) then
     sinpeak = node_value('sinpeak ')
     sintune = node_value('sintune ')
     sinphase = node_value('sinphase ')
     temp = sinpeak * sin(sinphase + twopi * sintune * turn)
     if (code.eq.14) then
        xkick=xkick+temp
     else if (code.eq.15) then
        xkick=xkick+temp
        ykick=ykick+temp
     else if(code.eq.16) then
        ykick=ykick+temp
     endif
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
  !     bil2 * (px**2 + py**2 + bi2gi2*pt**2) + el*pt*dtbyds
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
     !        (dpxx*track(1,itrack) - dpyy*track(3,itrack)) *                &
     !        ((one + bet0*track(6,itrack))/ddd)*bet0i

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

end subroutine ttcorr
subroutine tt_putone(npart,turn,tot_segm,segment,part_id,z,orbit0,&
     spos,ielem,el_name)

  use name_lenfi
  implicit none

  !hbu added spos, ielem, el_name
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
  integer i,j,npart,turn,tot_segm,segment,part_id(*),length
  double precision z(6,*),orbit0(6),tmp,tt,ss
  !hbu was *36 allow longer info
  character(120) table,comment
  !hbu
  integer ielem
  !hbu name of element
  character(name_len) el_name
  !hbu
  double precision spos
  !hbu
  character(4) vec_names(7)
  !hbu
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt','s' /
  data table / 'trackone' /

  !hbu
  length = len(comment)
  segment = segment + 1
  !hbu
  write(comment, '(''#segment'',4i8,1X,A)')                         &
       segment,tot_segm,npart,ielem,el_name
  call comment_to_table_curr(table, comment, length)
  tt = turn
  do i = 1, npart
     call double_to_table_curr(table, 'turn ', tt)
     ss = part_id(i)
     call double_to_table_curr(table, 'number ', ss)
     do j = 1, 6
        tmp = z(j,i) - orbit0(j)
        call double_to_table_curr(table, vec_names(j), tmp)
     enddo
     !hbu spos
     call double_to_table_curr(table,vec_names(7),spos)
     call augment_count(table)
  enddo
end subroutine tt_putone



subroutine tt_puttab(npart,turn,nobs,orbit,orbit0,spos)

  implicit none

  !hbu added spos
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
  character(36) table
  !hbu
  double precision spos
  !hbu
  character(4) vec_names(7)
  !hbu
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt', 's' /
  data table / 'track.obs$$$$.p$$$$' /

  tt = turn
  tn = npart

  write(table(10:13), '(i4.4)') nobs
  write(table(16:19), '(i4.4)') npart

  energy = get_value('probe ','energy ')

  call double_to_table_curr(table, 'turn ', tt)   ! the number of the cur
  call double_to_table_curr(table, 'number ', tn) ! the number of the cur
  call double_to_table_curr(table, 'e ', energy)

  do j = 1, 6
     tmp = orbit(j) - orbit0(j)
     call double_to_table_curr(table, vec_names(j), tmp)
  enddo
  !hbu spos
  call double_to_table_curr(table,vec_names(7),spos)
  call augment_count(table)
end subroutine tt_puttab

subroutine trinicmd(switch,orbit0,eigen,jend,z,turns,coords)

  use bbfi
  use trackfi
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
  logical zgiv,zngiv
  integer j,jend,k,kp,kq,next_start,itype(23),switch,turns,jt
  double precision phi,track(12),zstart(12),twopi,z(6,*),zn(6),     &
       ex,ey,et,orbit0(6),eigen(6,6),x,px,y,py,t,deltae,fx,phix,fy,phiy, &
       ft,phit,get_value,get_variable,zero,deltax,coords(6,0:turns,*)
  double precision deltap,one
  parameter(zero=0d0,one=1d0)
  character(120) msg(2)

  deltap = get_variable('track_deltap ')
  !---- Initialise orbit, emittances and eigenvectors etc.
  j = 0
  twopi = get_variable('twopi ')
  ex = get_value('probe ','ex ')
  ey = get_value('probe ','ey ')
  et = get_value('probe ','et ')
  bet0  = get_value('beam ','beta ')
  bet0i = one / bet0
  !-----get x add-on for lyaponuv calculation from dynap table
  deltax = get_value('dynap ', 'lyapunov ')
1 continue
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
                + sqrt(ex) * (eigen(k,1) * zn(1) + eigen(k,2) * zn(2))            &
                + sqrt(ey) * (eigen(k,3) * zn(3) + eigen(k,4) * zn(4))            &
                + sqrt(et) * (eigen(k,5) * zn(5) + eigen(k,6) * zn(6))
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
              coords(k,0,j+1) = z(k,j+1)
           enddo
           z(1,j+1) = z(1,j) + deltax
           coords(1,0,j+1) = z(1,j+1)
        endif
     endif
     goto 1
  endif
  !      if (switch .eq. 3)  then
  !--- create second particle with x add-on
  !        deltax = get_value('dynap ', 'lyapunov ')
  !        jend = 2
  !        z(1,jend) = z(1,1) + deltax
  !        do k = 2, 6
  !          z(k,jend) = z(k,1)
  !        enddo
  !      endif
end subroutine trinicmd
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
end subroutine trbegn
subroutine ttdpdg(track, ktrack)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose: computes the effect of the dipole edge                      *
  ! Input/Output:  ktrack , number of surviving trajectories             *
  !                 track , trajectory coordinates                       *
  !----------------------------------------------------------------------*
  integer ktrack,itrack
  double precision fint,e1,h,hgap,corr,rw(6,6),tw(6,6,6),track(6,*),&
       node_value,zero,one
  parameter(zero=0d0,one=1d0)

  call m66one(rw)
  !      call dzero(tw, 216)

  e1 = node_value('e1 ')
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
end subroutine ttdpdg

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
  double precision beta
  double precision get_value, node_value
  double precision track(6,*),one,two
  double precision sk,skl,sks,sksl,cosTh,sinTh,Q,R,Z
  double precision xf,yf,pxf,pyf,sigf,psigf,bvk
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
  sksl    = node_value('ksi ')
  sks     = node_value('ks ')

  !---- BV flag
  bvk = node_value('other_bv ')
  sks = sks * bvk
  sksl = sksl * bvk

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
end subroutine trsol

subroutine tttrans(track,ktrack)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Translation of a set of particles.                                 *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !----------------------------------------------------------------------*
  integer itrack,ktrack
  double precision node_value
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
end subroutine tttrans

subroutine tttrak(ek,re,track,ktrack)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of particle with a given TRANSPORT map.                *
  ! Input:                                                               *
  !   EK(6)     (real)    Kick.                                          *
  !   RE(6,6)   (real)    First-order terms.                             *
  ! Input/output:                                                        *
  !   TRACK(6,*)(real)    Track coordinates: (X, PX, Y, PY, T, PT).      *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  integer i,itrack,ktrack
  double precision ek(6),re(36),temp(6),track(6,*)
  !      double precision se(36)

  do itrack = 1, ktrack
     !        do i = 1, 36
     !          se(i) = re(i)                                                 &
     !           + te(i,1)*track(1,itrack) + te(i,2) * track(2,itrack)       &
     !           + te(i,3)*track(3,itrack) + te(i,4) * track(4,itrack)       &
     !           + te(i,5)*track(5,itrack) + te(i,6) * track(6,itrack)
     !        enddo
     do i = 1, 6
        temp(i) = ek(i)                                               &
             + re(i)    * track(1,itrack) + re(i+ 6) * track(2,itrack)         &
             + re(i+12) * track(3,itrack) + re(i+18) * track(4,itrack)         &
             + re(i+24) * track(5,itrack) + re(i+30) * track(6,itrack)
     enddo
     do i = 1, 6
        track(i,itrack) = temp(i)
     enddo
  enddo
end subroutine tttrak

subroutine trupdate(turn)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Update global turnnumber-VARIABLE,                                 *
  !   which can be used to make other variables turn-dependent.          *
  ! Input/output:                                                        *
  !   turn     (integer)    Current turn number.                         *
  !----------------------------------------------------------------------*
  integer           turn
  character(50)      cmd

  !---- call pro_input('TR$TURN := turn;')
  write(cmd, '(''tr$turni := '',i8,'' ; '')') turn
  call pro_input(cmd)
  write(cmd, '(''exec, tr$macro($tr$turni) ; '')')
  call pro_input(cmd)
end subroutine trupdate

subroutine trclor(orbit0)

  use twiss0fi
  use name_lenfi
  use trackfi
  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   -  Find 6D closed orbit                                            *
  !   -  Use Newton approximation                                        *
  !                                                                      *
  ! Input/output:                                                        *
  !   ORBIT0(6) (double)  Closed Orbit Coordinates @ sequence start      *
  !                       (X, PX, Y, PY, T, PT)                          *
  !----------------------------------------------------------------------*
  double precision zero, one
  parameter(zero=0d0,one=1d0)
  double precision orbit0(6),z(6,7),zz(6),z0(6,7),z00(6,7),a(6,7),deltap,ddd(6)
  integer itra, itmax, j, bbd_pos, j_tot
  integer code
  double precision el,dxt(200),dyt(200)
  integer restart_sequ
  integer advance_node, get_option, node_al_errors
  double precision node_value, get_value, get_variable
  integer n_align
  double precision al_errors(align_max)
  double precision re(6,6)

  logical aperflag, onepass
  integer pmax, turn, turnmax, part_id(1),last_turn(1)
  double precision sum, orbit(6)
  double precision last_pos(6),last_orbit(6,1),maxaper(6)

  integer nint, ndble, nchar, int_arr(1), char_l
  character(12) char_a

  integer i,k, irank

  double precision cotol, err

  logical is_drift, is_thin, is_quad, is_matrix

  print *," "
  !      print *," AK special version 2007/12/13"
  !      print *," ============================="

  !      print *," Modified by Yipeng SUN 19/11/2008 "
  !      print *," ======================="
  print *," Full 6D closed orbit search."
  print*," Initial value of 6-D closed orbit from Twiss: "
  print*," orbit0 ",orbit0

  !---- Initialize needed variables
  turn    = 1
  turnmax = 1
  itmax   = 10
  pmax    = 7
  sum     = zero
  aperflag  = .false.
  onepass   = .true.

  do k=1,7
     call dcopy(orbit0,z(1,k),6)
  enddo

  !      ddd(1) = orbit0(1)/100000
  !      ddd(2) = orbit0(2)/100000
  !      ddd(3) = orbit0(3)/100000
  !      ddd(4) = orbit0(4)/100000
  !      ddd(5) = orbit0(5)/100000
  !      ddd(6) = orbit0(6)/100000

  ! for on momentum pt =0
  !      ddd(1) = 1e-15
  !      ddd(2) = 1e-15
  !      ddd(3) = 1e-15
  !      ddd(4) = 1e-15
  !      ddd(5) = 1e-15
  !      ddd(6) = 1e-15
  !--------------------------------------
  ddd(1) = 1e-15
  ddd(2) = 1e-15
  ddd(3) = 1e-15
  ddd(4) = 1e-15
  ddd(5) = 1e-15
  ddd(6) = 1e-15

  !      do k=1,6
  !        z(k,k+1) = z(k,k+1) + ddd(k)
  !      enddo

  do k=1,7
     call dcopy(z(1,k),z0(1,k),6)
     call dcopy(z(1,k),z00(1,k),6)
  enddo

  !--- jmax may be reduced by particle loss - keep number in j_tot
  j_tot = pmax
  !--- get vector of six coordinate maxapers (both RUN and DYNAP)
  call comm_para('maxaper ',nint,ndble,nchar,int_arr,maxaper,       &
       char_a, char_l)
  !--- set particle id




  cotol = get_variable('twiss_tol ')

  !---- Initialize kinematics and orbit
  bet0  = get_value('beam ','beta ')
  betas = get_value('probe ','beta ')
  gammas= get_value('probe ','gamma ')
  bet0i = one / bet0
  beti   = one / betas
  dtbyds = get_value('probe ','dtbyds ')
  deltas = get_variable('track_deltap ')
  deltap = get_value('probe ','deltap ')
  arad = get_value('probe ','arad ')
  dorad = get_value('probe ','radiate ') .ne. 0
  dodamp = get_option('damp ') .ne. 0
  dorand = get_option('quantum ') .ne. 0
  call dcopy(orbit0,orbit,6)


  !---- Iteration for closed orbit.
  do itra = 1, itmax
     !        print*," before tracking "
     !        print*," itra: ",itra, " z: ", z
     j = restart_sequ()

     !---- loop over nodes
10   continue


     bbd_pos = j
     code    = node_value('mad8_type ')
     if(code.eq.39) code=15
     if(code.eq.38) code=24
     el      = node_value('l ')
     if (itra .eq. 1)  then
        if (.not.(is_drift() .or. is_thin() .or. is_quad() .or. is_matrix())) then
           print*," "
           print*,"code: ",code," el: ",el,"   THICK ELEMENT FOUND"
           print*," "
           print*,"Track dies nicely"
           print*,"Thick lenses will get nowhere"
           print*,"MAKETHIN will save you"
           print*," "
           print*," "
           stop
        endif
     endif



     !--------  Misalignment at beginning of element (from twissfs.f)
     if (code .ne. 1)  then
        call dzero(al_errors, align_max)
        n_align = node_al_errors(al_errors)
        if (n_align .ne. 0)  then
           do i = 1, pmax
              call dcopy(z(1,i),zz,6)
              call tmali1(zz,al_errors, betas, gammas,z(1,i), re)
           enddo
        endif
     endif

     !-------- Track through element
     call ttmap(code,el,z,pmax,dxt,dyt,sum,turn,part_id,             &
          last_turn,last_pos, last_orbit,aperflag,maxaper,al_errors,onepass)

     !--------  Misalignment at end of element (from twissfs.f)
     if (code .ne. 1)  then
        if (n_align .ne. 0)  then
           do i = 1, pmax
              call dcopy(z(1,i),zz,6)
              call tmali2(el,zz, al_errors, betas, gammas,z(1,i), re)
           enddo
        endif
     endif


     if (advance_node().ne.0)  then
        j=j+1
        goto 10
     endif
     !---- end of loop over nodes


     !---- construct one-turn map
     do k=1,6
        do i=1,6
           a(i,k) = (z(i,k+1) - z(i,1))/ddd(i)
        enddo
     enddo
     !---- Solve for dynamic case.
     err = zero
     do i= 1,6
        a(i,i) = a(i,i) - one
        a(i,7) = z(i,1) - z0(i,1)
        err = max(abs(a(i,7)), err)
     enddo

     call solver(a,6,1,irank)
     if (irank.lt.6) go to 100
     do i = 1, 6
        z0(i,1) = z0(i,1) - a(i,7)
        z00(i,1) = z00(i,1) - a(i,7)
     enddo
     !---- Solve for dynamic case.
     !      do i = 1, 6	
     !      print*," a(i,7) ",a(i,7)
     !      enddo

     !      if (err.lt.cotol) then

     !      print *, ' '
     !      print '(''iteration: '',i3,'' error: '',1p,e14.6,'' deltap: ''    &
     !     ,                                                                 &
     !     1p,e14.6)',itra,err,deltap
     !          print '(''orbit: '', 1p,6e14.6)', z0

     !      goto 110

     !      endif



     do k=2,7
        call dcopy(z00(1,1),z0(1,k),6)
     enddo

     do k=1,6
        z0(k,k+1) = z0(k,k+1) + ddd(k)
     enddo

     do k=1,7
        call dcopy(z0(1,k),z(1,k),6)
     enddo


     !---- end of Iteration
  enddo

  call dcopy(z0(1,1),orbit0,6)
  !      print *,"madX::trrun.F::trclor()"
  !      print *," Iteration maximum reached. NO convergence."
  goto 110

100 continue
  print *," Singular matrix occurred during closed orbit search."

110 continue
  print *, ' '
  print *, '6D closed orbit found by subroutine trclor '
  print '(''iteration: '',i3,'' error: '',1p,e14.6,'' deltap: '',1p,e14.6)',itra,err,deltap
  print '(''orbit: '', 1p,6e14.6)', orbit0

end subroutine trclor

subroutine ttnllens(track,ktrack)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of trajectories through a thin nonlinear lens          *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  ! Added by Alexander Valishev 22 Oct 2010                              *
  !----------------------------------------------------------------------*
  integer itrack,ktrack
  double precision track(6,*),node_value,get_variable,    &
       pi,one,two,knll,cnll, &
       dd, u, v, dUu, dUv, dux, duy, dvx, dvy, x, y
  parameter(one=1d0,two=2d0)
  
  cnll=node_value('cnll ')
  knll=node_value('knll ')/cnll
  pi=get_variable('pi ')

  do itrack = 1, ktrack

    x = track(1,itrack)/cnll
    y = track(3,itrack)/cnll

    u=0.5*sqrt((x-1)**2+y**2)+0.5*sqrt((x+1)**2+y**2)
    v=0.5*sqrt((x+1)**2+y**2)-0.5*sqrt((x-1)**2+y**2)
    if (u.eq.1d0) then 
       dd=0 
    else 
       dd=u**2*log(u+sqrt(u*u-1))/sqrt(u**2-1)
    endif
    dUu=(u+log(u+sqrt(u*u-1))*sqrt(u**2-1)+dd)/(u**2-v**2) &
     -2*u*(u*log(u+sqrt(u*u-1))*sqrt(u**2-1) &
     +v*(acos(v)-0.5*pi)*sqrt(1-v**2)) /(u**2-v**2)**2
    dUv=2*v*(u*log(u+sqrt(u*u-1))*sqrt(u**2-1) &
         +v*(acos(v)-0.5*pi)*sqrt(1-v**2)) /(u**2-v**2)**2 &
         -(v-(acos(v)-0.5*pi)*sqrt(1-v**2)+v**2*(acos(v)-0.5*pi)/sqrt(1-v**2))&
         /(u**2-v**2)
    dux=0.5*(x-1)/sqrt((x-1)**2+y**2) +0.5*(x+1)/sqrt((x+1)**2+y**2)
    duy=0.5*y/sqrt((x-1)**2+y**2) +0.5*y/sqrt((x+1)**2+y**2)
    dvx=0.5*(x+1)/sqrt((x+1)**2+y**2) -0.5*(x-1)/sqrt((x-1)**2+y**2)
    dvy=0.5*y/sqrt((x+1)**2+y**2) -0.5*y/sqrt((x-1)**2+y**2)

    track(2,itrack)=track(2,itrack)+knll*(dUu*dux+dUv*dvx)
    track(4,itrack)=track(4,itrack)+knll*(dUu*duy+dUv*dvy)

  enddo
end subroutine ttnllens

subroutine ttrfmult(track, ktrack, turn)

  use twtrrfi
  use trackfi
  implicit none

  !--------------------*
  ! Andrea Latina 2012 *
  !--------------------*
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !    Track particle through a general thin rf-multipole.               *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer)   Number of surviving tracks.                  *
  !----------------------------------------------------------------------*

  double precision beta, angle, cangle, sangle, dtmp
  double precision node_value, get_value
  double precision bvk, deltap, elrad
  double precision f_errors(0:maxferr)
  double precision vals(2,0:maxmul)
  integer n_ferr, jtrk, iord, nord

  !--- AL: RF-multipole
  integer dummyi, nn, ns, node_fd_errors, turn, ktrack
  double precision normal(0:maxmul), skew(0:maxmul)
  double precision track(6,*), field(2,0:maxmul)
  double precision field_cos(2,0:maxmul)
  double precision field_sin(2,0:maxmul)
  double precision zero, one, two, three
  double precision pc, krf, rfac
  double precision pi, clight, ten3m
  double precision x, y, z, dpx, dpy, dpt
  double precision freq, volt, lag, harmon
  double precision pnl(0:maxmul), psl(0:maxmul)
  complex*16 ii, Cm2, Sm2, Cm1, Sm1, Cp0, Sp0, Cp1, Sp1
    
  parameter ( pi = 3.14159265358979d0 )
  parameter ( clight = 299792458d0 )
  parameter ( zero=0d0, one=1d0, two=2d0, three=3d0, ten3m=1d-3)
  parameter ( ii=(0d0,1d0) )
  
  !---- Zero the arrays
  call dzero(normal,maxmul+1)
  call dzero(skew,maxmul+1)
  call dzero(pnl,maxmul+1)
  call dzero(psl,maxmul+1)
  call dzero(f_errors,maxferr+1)
  call dzero(field,2*(maxmul+1))
  
  !---- Read-in the parameters
  freq = node_value('freq ');
  volt = node_value('volt ');
  lag = node_value('lag ');
  harmon = node_value('harmon ');
  bvk = node_value('other_bv ')
  elrad = node_value('lrad ')
  deltap = get_value('probe ', 'deltap ')
  dorad = get_value('probe ','radiate ') .ne. zero
  arad = get_value('probe ','arad ')
  gammas = get_value('probe ','gamma ')
  beta = get_value('probe ','beta ')
  pc = get_value('probe ','pc ')
  n_ferr = node_fd_errors(f_errors);
  call get_node_vector('knl ', nn, normal)
  call get_node_vector('ksl ', ns, skew)
  call get_node_vector('pnl ', dummyi, pnl);
  call get_node_vector('psl ', dummyi, psl);
 
  rfac = zero
  
  !---- Set-up some parameters
  krf = 2*pi*freq*1d6/clight;
  
  if (n_ferr.gt.0) then
     call dcopy(f_errors,field,n_ferr)
  endif
  nord = max(nn, ns, n_ferr/2-1);
      
  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1,ktrack
    x = track(1,jtrk);
    y = track(3,jtrk);
    z = track(5,jtrk);
    !---- Vector with strengths + field errors
    do iord = 0, nord;
      field_cos(1,iord) = bvk * (normal(iord) * cos(pnl(iord) * 2 * pi - krf * z) + field(1,iord));
      field_sin(1,iord) = bvk * (normal(iord) * sin(pnl(iord) * 2 * pi - krf * z));
      field_cos(2,iord) = bvk * (skew(iord)   * cos(psl(iord) * 2 * pi - krf * z) + field(2,iord));
      field_sin(2,iord) = bvk * (skew(iord)   * sin(psl(iord) * 2 * pi - krf * z));
    enddo
    Cm2 = 0d0;
    Sm2 = 0d0;
    Cm1 = 0d0;
    Sm1 = 0d0;
    Cp0 = 0d0;
    Sp0 = 0d0;
    Cp1 = 0d0;
    Sp1 = 0d0;
    do iord = nord, 0, -1
      if (iord.ge.2) then
        Cm2 = Cm2 * (x+ii*y) / (iord-1) + field_cos(1,iord)+ii*field_cos(2,iord);
        Sm2 = Sm2 * (x+ii*y) / (iord-1) + field_sin(1,iord)+ii*field_sin(2,iord);
      endif
      if (iord.ge.1) then
        Cm1 = Cm1 * (x+ii*y) / (iord)   + field_cos(1,iord)+ii*field_cos(2,iord);
        Sm1 = Sm1 * (x+ii*y) / (iord)   + field_sin(1,iord)+ii*field_sin(2,iord);
      endif
      Cp0 = Cp0 * (x+ii*y) / (iord+1)   + field_cos(1,iord)+ii*field_cos(2,iord);
      Sp0 = Sp0 * (x+ii*y) / (iord+1)   + field_sin(1,iord)+ii*field_sin(2,iord);
      Cp1 = Cp1 * (x+ii*y) / (iord+2)   + field_cos(1,iord)+ii*field_cos(2,iord);
      Sp1 = Sp1 * (x+ii*y) / (iord+2)   + field_sin(1,iord)+ii*field_sin(2,iord);
    enddo
    Sp1 = Sp1 * (x+ii*y);
    Cp1 = Cp1 * (x+ii*y);

    !---- The kick
    dpx = -REAL(Cp0) / (one + track(6,jtrk));
    dpy = AIMAG(Cp0) / (one + track(6,jtrk));
    dpt = (volt * ten3m * sin(lag * 2 * pi - krf * z) / pc - krf * REAL(Sp1)) / (one + track(6,jtrk));

    !---- Radiation effects at entrance.
    if (dorad  .and.  elrad .ne. zero) then
      rfac = arad * gammas**3 * (dpx**2+dpy**2) / (three*elrad)
      track(2,jtrk) = track(2,jtrk) - rfac * (one + track(6,jtrk)) * track(2,jtrk)
      track(4,jtrk) = track(4,jtrk) - rfac * (one + track(6,jtrk)) * track(4,jtrk)
      track(6,jtrk) = track(6,jtrk) - rfac * (one + track(6,jtrk)) ** 2
    endif

    !---- Apply the kick
    track(2,jtrk) = track(2,jtrk) + dpx
    track(4,jtrk) = track(4,jtrk) + dpy
    track(6,jtrk) = track(6,jtrk) + dpt

    !---- Radiation effects at exit.
    if (dorad  .and.  elrad .ne. zero) then
      track(2,jtrk) = track(2,jtrk) - rfac * (one + track(6,jtrk)) * track(2,jtrk)
      track(4,jtrk) = track(4,jtrk) - rfac * (one + track(6,jtrk)) * track(4,jtrk)
      track(6,jtrk) = track(6,jtrk) - rfac * (one + track(6,jtrk)) ** 2
    endif
  enddo   

end subroutine ttrfmult

subroutine tttquad(track, ktrack)
  
  use twtrrfi
  use trackfi
  implicit none
  
  !--------------------*
  ! Andrea Latina 2012 *
  !--------------------*
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !    Track particle through a general thick quadrupole.                *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer)   Number of surviving tracks.                  *
  !----------------------------------------------------------------------*
  
  double precision track(6,*)
  integer ktrack
  
  double precision node_value, get_value
  double precision k1, k1s, length, ksqrt, tmp
  double precision sx, cx, sy, cy, ct, st
  double precision x, px, y, py, z, pz
  logical skew, focusing
  integer jtrk
  
  double precision zero, one, two, sqrt2
  parameter ( zero=0d0, one=1d0, two=2d0 )
  parameter ( sqrt2=1.41421356237310d0 )
  
  !---- Read-in the parameters
  k1 = node_value('k1 ');
  k1s = node_value('k1s ');
  length = node_value('l ');
  
  if ((k1.ne.zero).and.(k1s.ne.zero)) then
     call aawarn('trrun: ',&
          'a quadrupole cannot have both K1 and K1S different than zero');
     return  
  endif

  if (k1s.ne.zero) then
     skew = .true.
     focusing = k1s.gt.zero
  else
     skew = .false.
     focusing = k1.gt.zero
  endif
  
  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1,ktrack
     !---- The particle position
     x  = track(1,jtrk);
     px = track(2,jtrk);
     y  = track(3,jtrk);
     py = track(4,jtrk);
     z  = track(5,jtrk);
     pz = track(6,jtrk);
  
!!$    !---- Radiation effects at entrance.
!!$    if (dorad  .and.  elrad .ne. zero) then
!!$      rfac = arad * gammas**3 * (dpx**2+dpy**2) / (three*elrad)
!!$      track(2,jtrk) = track(2,jtrk) - rfac * (one + track(6,jtrk)) * track(2,jtrk)
!!$      track(4,jtrk) = track(4,jtrk) - rfac * (one + track(6,jtrk)) * track(4,jtrk)
!!$      track(6,jtrk) = track(6,jtrk) - rfac * (one + track(6,jtrk)) ** 2
!!$    endif
     
     !---- If SKEW rotates by -45 degrees
     if (skew) then
        ct =  sqrt2 / two
        st = -sqrt2 / two
        tmp = x
        x = ct * tmp + st * y
        y = ct * y   - st * tmp
        tmp = px
        px = ct * tmp + st * py
        py = ct * py  - st * tmp
     endif
     
     !---- Computes the kick
     if (skew) then
        ksqrt = sqrt(abs(k1s / (one + track(6, jtrk))))
     else
        ksqrt = sqrt(abs(k1  / (one + track(6, jtrk))))
     endif
     if (focusing) then
        cx = cos(ksqrt*length)
        sx = sin(ksqrt*length)
        cy = cosh(ksqrt*length)
        sy = sinh(ksqrt*length)
        tmp = x
        x  = tmp * cx + px * sx / ksqrt
        px = -sx * ksqrt * tmp + px * cx
        tmp = y
        y  = tmp * cy + py * sy / ksqrt
        py =  sy * ksqrt * tmp + py * cy
     else
        cx = cosh(ksqrt*length)
        sx = sinh(ksqrt*length)
        cy = cos(ksqrt*length)
        sy = sin(ksqrt*length)
        tmp = x
        x  = tmp * cx + px * sx / ksqrt
        px =  sx * ksqrt * tmp + px * cx
        tmp = y
        y  = tmp * cy + py * sy / ksqrt
        py = -sy * ksqrt * tmp + py * cy
     endif
     
     !---- If SKEW rotates by +45 degrees
     if (skew) then
        ct = sqrt2 / two
        st = sqrt2 / two
        tmp = x
        x = ct * tmp + st * y
        y = ct * y   - st * tmp
        tmp = px
        px = ct * tmp + st * py
        py = ct * py  - st * tmp
     endif
     
     !---- Applies the kick
     track(1,jtrk) = x
     track(2,jtrk) = px
     track(3,jtrk) = y
     track(4,jtrk) = py
     
!!$    !---- Radiation effects at exit.
!!$    if (dorad  .and.  elrad .ne. zero) then
!!$      track(2,jtrk) = track(2,jtrk) - rfac * (one + track(6,jtrk)) * track(2,jtrk)
!!$      track(4,jtrk) = track(4,jtrk) - rfac * (one + track(6,jtrk)) * track(4,jtrk)
!!$      track(6,jtrk) = track(6,jtrk) - rfac * (one + track(6,jtrk)) ** 2
!!$    endif
     
  enddo
  
end subroutine tttquad
