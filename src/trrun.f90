LOGICAL FUNCTION  is_thin()
  double precision node_value, el
  el = node_value('l ')
  is_thin = el.eq.0d0;
END FUNCTION is_thin

LOGICAL FUNCTION  is_drift()
  double precision node_value, code
  code = node_value('mad8_type ')
  is_drift = code.eq.1;
END FUNCTION is_drift

LOGICAL FUNCTION  is_dipole()
  double precision node_value, code
  code = node_value('mad8_type ')
  is_dipole = code.eq.3;
END FUNCTION is_dipole

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

subroutine trrun(switch,turns,orbit0,rt,part_id,last_turn,        &
     last_pos,z,dxt,dyt,last_orbit,eigen,coords,e_flag,code_buf,l_buf)

  use twtrrfi
  use bbfi
  use time_varfi
  use spch_bbfi
  use twiss0fi
  use name_lenfi
  use trackfi
  use fasterror
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
  logical onepass,onetable,last_out,info,aperflag,doupdate,first,        &
       bb_sxy_update,virgin_state,emittance_update,checkpnt_restart,     &
       fast_error_func,exit_loss_turn
  integer j,code,restart_sequ,advance_node,                              &
       node_al_errors,n_align,nlm,jmax,j_tot,turn,turns,i,k,get_option,  &
       ffile,SWITCH,nint,ndble,nchar,part_id(*),last_turn(*),char_l,     &
       segment, e_flag, nobs,lobs,int_arr(1),tot_segm,code_buf(*),       &
       tot_turn,max_part
  parameter(max_part=20000)
  integer part_id_keep(max_part),last_turn_keep(max_part)
  double precision tmp_d,orbit0(6),orbit(6),el,rt(6,6),re(6,6),          &
       al_errors(align_max),z(6,*),zz(6),dxt(*),dyt(*),eigen(6,6),sum,   &
       node_value,get_variable,last_pos(*),last_orbit(6,*),              &
       maxaper(6),get_value,obs_orb(6),coords(6,0:turns,*),l_buf(*),     &
       betx_start,bety_start,alfx_start,alfy_start,gamx_start,gamy_start,&
       dx_start,dpx_start,dy_start,dpy_start,deltap,                     &
       N_ions_in_beam, Npart_gain, t_rms, pt_rms,                        &
       N_ions_ini, n_ions_macro, sigma_z_ini, z_factor,                  &
       N_ions_for_bb,z_keep(6,max_part),ex_rms0,ey_rms0,sigma_p0,sigma_z0
  data ex_rms0,ey_rms0,sigma_p0,sigma_z0 / 0d0, 0d0, 0d0, 0d0 /

  character(12) tol_a, char_a
  character(20) text

!VVK 20100321 -------------------------------------------------
      integer i_part                      ! local counter
      double precision   Summ_t_mean      ! local for mean value
      double precision   Summ_t_square    ! local for rms value
!-------------------------------------------------------------------

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
  save first
  data first / .true. /
  save tot_turn
  data tot_turn / 0 /
  save betx_start, bety_start,                             &
       alfx_start, alfy_start,                             &
       gamx_start, gamy_start,                             &
       dx_start,    dpx_start,                             &
       dy_start,    dpy_start,                             &
       N_ions_in_beam, Npart_gain, t_rms, pt_rms,          &
       N_ions_ini, n_ions_macro, sigma_z_ini, z_factor,    &
       N_ions_for_bb,z_keep,part_id_keep,                  &
       last_turn_keep,jmax,segment,                        &
       ex_rms0,ey_rms0,sigma_p0,sigma_z0

  data betx_start, bety_start,                             &
       alfx_start, alfy_start,                             &
       gamx_start, gamy_start,                             &
       dx_start,    dpx_start,                             &
       dy_start,    dpy_start / 1d0, 1d0, 0d0, 0d0,        &
       0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /

  logical is_drift, is_thin, is_quad, is_matrix, is_dipole

  !-------added by Yipeng SUN 01-12-2008--------------
  deltap = get_value('probe ','deltap ')

  if(deltap.eq.0d0) then
     onepass = get_option('onepass ') .ne. 0
     if(onepass)   then
     else
        call trclor(switch,orbit0)
     endif
  else
  endif

  !---- Set fast_error_func flag to use faster error function
  !---- including tables. Thanks to late G. Erskine
  fast_error_func = get_option('fast_error_func ') .ne. 0
  if(fast_error_func.and..not.fasterror_on) then
     call wzset
     fasterror_on = .true.
  endif

  if(switch.eq.1) then
     bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
!FRS open unit 90 filename=checkpoint_restart.dat to prolong run
     open(90,file='checkpoint_restart.dat',form='unformatted',status='unknown')
     checkpnt_restart = get_value('run ', 'checkpnt_restart ') .ne. 0d0
  else
     bb_sxy_update = .false.
     checkpnt_restart = .false.
  endif
  if(bb_sxy_update) then
     emittance_update = get_option('emittance_update ') .ne. 0
     virgin_state = get_value('run ', 'virgin_state ') .ne. 0d0
     if(virgin_state) first=.true.

     call table_input(                                 &
          betx_start, bety_start,                      &
          alfx_start, alfy_start,                      &
          gamx_start, gamy_start,                      &
          dx_start,    dpx_start,                      &
          dy_start,    dpy_start)

     if(first) call make_bb6d_ixy(turns)
     exit_loss_turn = get_option('exit_loss_turn ') .ne. 0
  endif

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
  if(first) then
     segment = 0
  endif
  tot_segm = turns / ffile + 1
  if (mod(turns, ffile) .ne. 0) tot_segm = tot_segm + 1

  if(first) then
     if(checkpnt_restart) then
        read(90,END=100) jmax
        read(90,END=100) Ex_rms0
        read(90,END=100) Ey_rms0
        do i = 1, jmax
           do j=1,6
              read(90,END=100) z(j,i)
           enddo
        enddo
     else
        call trinicmd(switch,orbit0,eigen,jmax,z,turns,coords)
     endif
     !--- set particle id
!$OMP PARALLEL PRIVATE(k)
!$OMP DO
     do k=1,jmax
        part_id(k) = k
     enddo
!$OMP END DO
!$OMP END PARALLEL
  else
     if(jmax.eq.0) then
        call aawarn('trrun: ',&
             'Particle Number equals zero: exit from trrun')
        return
     endif
     if(jmax.gt.max_part) then
        write(text, '(1p,d20.12)') max_part
        call aafail('TRRUN: ','Fatal: Maximum Particle Number exceeded =' // text)
     endif
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
     do i = 1, jmax
        last_turn(i)=last_turn_keep(i)
        part_id(i)=part_id_keep(i)
        do j=1,6
           z(j,i)=z_keep(j,i)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

  !--- jmax may be reduced by particle loss - keep number in j_tot
  j_tot = jmax
  !--- get vector of six coordinate maxapers (both RUN and DYNAP)
  call comm_para(tol_a, nint, ndble, nchar, int_arr, maxaper,       &
       char_a, char_l)
  !hbu--- init info for tables initial s position is 0
  !hbu initial s position is 0
  spos=0
  !hbu start of line, element 0
  nlm=0
  !hbu
  el_name='start           '

  !---- Initialize kinematics and orbit
  bet0  = get_value('beam ','beta ')
  betas = get_value('probe ','beta ')
  gammas= get_value('probe ','gamma ')
  bet0i = 1d0 / bet0
  beti   = 1d0 / betas
  dtbyds = get_value('probe ','dtbyds ')
  deltas = get_variable('track_deltap ')
  arad = get_value('probe ','arad ')
  dorad = get_value('probe ','radiate ') .ne. 0d0
  dodamp = get_option('damp ') .ne. 0
  dorand = get_option('quantum ') .ne. 0

  if(first) then
     !--- enter start coordinates in summary table
     if (switch .eq. 1)  then
        t_max=get_value('probe ','circ ')/get_value('run ', 'track_harmon ')/betas
        pt_max = get_value('run ', 'deltap_max ')
        pt_max=(sqrt((betas*(pt_max+1d0))**2+1d0/gammas**2)-1d0)*beti
     else
        t_max=1d20
        pt_max=1d20
     endif
     do  i = 1,j_tot
        if(abs(z(5,i)).gt.t_max) then
           write(text, '(1p,d13.5,a1,i6)') t_max,"p",i
           call aafail('TRACK_INITIAL: ','Fatal: T-Coordinate larger than' // text)
        endif
        if(abs(z(6,i)).gt.pt_max) then
           write(text, '(1p,d13.5,a1,i6)') pt_max,"p",i
           call aafail('TRACK_INITIAL: ','Fatal: PT-Coordinate larger than' // text)
        endif
        tmp_d = i
        call double_to_table_curr('tracksumm ', 'number ', tmp_d)
        tmp_d = tot_turn
        call double_to_table_curr('tracksumm ', 'turn ', tmp_d)
        do j = 1, 6
           tmp_d = z(j,i) - orbit0(j)
           call double_to_table_curr('tracksumm ', vec_names(j), tmp_d)
        enddo
        !hbu add s
        call double_to_table_curr('tracksumm ',vec_names(7),spos)
        call augment_count('tracksumm ')
     enddo
  endif
  !--- enter first turn, and possibly eigen in tables
  if (switch .eq. 1)  then
     if (onetable)  then
        if(first) then
           call track_pteigen(eigen)
           !hbu add s, node id and name
           call tt_putone(jmax, tot_turn, tot_segm, segment, part_id,           &
                z, orbit0,spos,nlm,el_name)
        endif
     else
        do i = 1, jmax
           !hbu
           call tt_puttab(part_id(i), 0, 1, z(1,i), orbit0,spos)
        enddo
     endif
  endif
  call dcopy(orbit0,orbit,6)

  doupdate = get_option('update ') .ne. 0

  !--- loop over turns
  nobs = 0

  if(bb_sxy_update) then
     trrun_nt=0

     if(first) then
        time_var_m_cnt=0
        time_var_p_cnt=0
        time_var_c_cnt=0
        !      <<N_macro_part_ini=N_macro_surv+N_macro_lost>>
        N_ions_in_beam=get_value('probe ', 'npart ') !BEAM->NPART
        if(N_ions_in_beam .LT. 0d0) call aafail('TRRUN: ','N_ions_in_beam .LE. 0.0')
        Npart_gain = get_value('run ', 'n_part_gain ')
        N_ions_ini=Npart_gain*N_ions_in_beam
        N_macro_surv=jmax    ! = number of START lines submitted
        n_ions_macro=N_ions_ini/N_macro_surv

        N_for_I=N_macro_surv ! at start (to be redefined in Ixy)
        if(N_macro_surv .GT. N_macro_max) call aafail('TRRUN: ',&
             'Number N_macro_surv exceeds N_macro_max (array size)')
        if(N_macro_surv .GT. N_macro_surv) call aafail('TRRUN: ',&
             'Number START-lines exceeds the initial number of macroparticles N_macro_surv')
        t_rms = get_value('run ', 'sigma_z ')*beti
        pt_rms = get_value('run ', 'deltap_rms ')
        pt_rms=(sqrt((betas*(pt_rms+1d0))**2+1d0/gammas**2)-1d0)*beti
        sigma_z_ini=t_rms !betas: BEAM->BETA
        sigma_z=sigma_z_ini !at start (to be redefined in Ixy)
        sigma_p=pt_rms       !default
        z_factor=1d0 !at start sigma_z_ini/sigma_z
        Ex_rms=get_value('probe ', 'ex ') !BEAM->Ex
        Ey_rms=get_value('probe ', 'ey ') !BEAM->Ey
        if(checkpnt_restart) then
           Ex_rms=Ex_rms0
           Ey_rms=Ey_rms0
        endif
!        write(8,'(4(g16.9,1x))') Ex_rms, Ey_rms,sigma_z,sigma_p
        first=.false.
     endif
  endif

  do turn = 1, turns

!--- Write checkpoint_restart data - disable for speed reasons
!     rewind 90
!     write(90) jmax
!     write(90) Ex_rms
!     write(90) Ey_rms
!     do i = 1, jmax
!        do j=1,6
!           write(90) z(j,i)
!        enddo
!     enddo

     if (doupdate) call trupdate(turn)

     j = restart_sequ()

     if(bb_sxy_update) then
        trrun_nt=turn
        time_var_m_lnt=0
        time_var_p_lnt=0
        time_var_c_lnt=0

        N_macro_surv=jmax
        i_spch = 0 !a special spch-update counter

        ex_rms0=ex_rms
        ey_rms0=ey_rms
        sigma_z0=sigma_z
!sigma_p0 = sigma_p !CM, 3/11/14
        !fill, table=Ixy_unsorted; column=i_macro_part, Ix, Iy, dpi, z_part;
!new on 3/31/14:
        if (emittance_update) then

        call ixy_calcs(betas, orbit0, z,                                &
             betx_start, bety_start,                      &
             alfx_start, alfy_start,                      &
             gamx_start, gamy_start,                      &
             dx_start,    dpx_start,                      &
             dy_start,    dpy_start)
        call ixy_fitting()

        call double_to_table_curr('bb6d_ixy ', 'turn ', dble(tot_turn+turn))
        call double_to_table_curr('bb6d_ixy ', 'n_macro_surv ', dble(n_macro_surv))
        call double_to_table_curr('bb6d_ixy ', 'n_for_i ', dble(n_for_i))
        call double_to_table_curr('bb6d_ixy ', 'ex_rms ', ex_rms)
        call double_to_table_curr('bb6d_ixy ', 'ey_rms ', ey_rms)
        call double_to_table_curr('bb6d_ixy ', 'sigma_p ', sigma_p)
        call double_to_table_curr('bb6d_ixy ', 'sigma_z ', sigma_z)
        call augment_count('bb6d_ixy ')

        if(sigma_p0.eq.0d0) sigma_p0=sigma_p
!new on 3/31/14:
        endif

        if(.not.emittance_update) then
           ex_rms=ex_rms0
           ey_rms=ey_rms0
           sigma_z=sigma_z0
           sigma_p=sigma_p0
        endif
!        write(8,'(4(g16.9,1x))') Ex_rms, Ey_rms,sigma_z,sigma_p
        if((sigma_z.GT.0d0).AND.(sigma_z_ini.GT.0d0)) then
           z_factor=sigma_z_ini/sigma_z
        else
           z_factor=1D0
        endif

        N_ions_for_bb=n_ions_macro*N_for_I*z_factor
        if(N_ions_in_beam.le.0d0) then
           rat_bb_n_ions=0d0
        else
           rat_bb_n_ions=N_ions_for_bb/N_ions_in_beam
        endif

        if(idnint(time_var_m_nt(time_var_m_cnt+1)).eq.tot_turn+turn) then
           time_var_m=.true.
        else
           time_var_m=.false.
        endif
        if(idnint(time_var_p_nt(time_var_p_cnt+1)).eq.tot_turn+turn) then
           time_var_p=.true.
        else
           time_var_p=.false.
        endif
        if(idnint(time_var_c_nt(time_var_c_cnt+1)).eq.tot_turn+turn) then
           time_var_c=.true.
        else
           time_var_c=.false.
        endif
     endif

     nlm = 0
     sum=0d0

     if((sigma_t.eq.0d0.or.emittance_update).and.bb_sxy_update) then
!VVK 20100321 -------- Find RMS-value of t ----------------------
! if we do 1-turn tracking, orbit0(5)=0 always
        Summ_t_mean=0d0
        do i_part=1, jmax
           if (abs(z(5,i_part)) .GE. 0d0) then
              Summ_t_mean=Summ_t_mean+z(5,i_part)
           else
              Print *, 'NaN z(5,i) ? :', i_part, z(5,i_part)
           endif
        enddo
        mean_t=Summ_t_mean/dble(jmax)

        Summ_t_square=0d0
        do i_part=1, jmax
           if (abs(z(5,i_part)) .GE. 0d0) &
                Summ_t_square=Summ_t_square+(z(5,i_part)-mean_t)**2
        enddo
        sigma_t=sqrt(Summ_t_square/dble(jmax))
        if(abs(sigma_t).eq.0d0) then
           sigma_t=t_max/2d0
           call aawarn('TTRUN Frozen SC: sigma_t = zero: ','sigma_t set to L/track_harmon/betas/2')
        endif
!-----------------------------------------------------------------
     endif

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
        if (.not.(is_drift() .or. is_thin() .or. is_quad() .or. is_dipole() .or. is_matrix())) then
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
           !            call aafail('TRRUN: Fatal ',                                 &
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
        !print *, "n_align = ", n_align
        !print *, "align_max = ", align_max
        if (n_align .ne. 0)  then
           do i = 1, jmax
              call dcopy(z(1,i),zz(1),6)
              call tmali1(zz(1),al_errors, betas, gammas,z(1,i), re(1,1))
           enddo
        endif
     endif
     !-------- Track through element  // suppress dxt 13.12.04
     call ttmap(switch,code,el,z,jmax,dxt,dyt,sum,tot_turn+turn,part_id,             &
          last_turn,last_pos, last_orbit,aperflag,maxaper,al_errors,onepass)
     !--------  Misalignment at end of element (from twissfs.f)
     if (bb_sxy_update .and. emittance_update .and. is_lost) then
        call ixy_calcs(betas, orbit0, z,                                &
             betx_start, bety_start,                      &
             alfx_start, alfy_start,                      &
             gamx_start, gamy_start,                      &
             dx_start,    dpx_start,                      &
             dy_start,    dpy_start)
        call ixy_fitting()
        is_lost = .false.
     endif
     if (code .ne. 1)  then
        if (n_align .ne. 0)  then
           do i = 1, jmax
              call dcopy(z(1,i),zz(1),6)
              call tmali2(el,zz(1),al_errors, betas, gammas,z(1,i),re(1,1))
           enddo
        endif
     endif
     nlm = nlm+1
     if (nobs .gt. 0)  then
        call dzero(obs_orb,6)
        call get_node_vector('obs_orbit ', lobs, obs_orb)
        if (lobs .lt. 6)                                              &
             call aafail('TRACK: Fatal ', 'obs. point orbit not found')
        if (onetable)  then
           !hbu
           spos=sum
           !hbu get current node name
           call element_name(el_name,len(el_name))
           !hbu
           call tt_putone(jmax, tot_turn+turn, tot_segm, segment, part_id,      &
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
              !2013-Nov-05  15:05:33  ghislain: get current node name
              call element_name(el_name,len(el_name))
              call tt_putone(jmax, tot_turn+turn, tot_segm, segment, part_id,    &
                   z, orbit0,spos,nlm,el_name)
           else
              do i = 1, jmax
                 !hbu
                 call tt_puttab(part_id(i), turn, 1, z(1,i), orbit0,spos)
              enddo
           endif
        endif
     else
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
        do i = 1, jmax
           do j = 1, 6
              coords(j,turn,i) = z(j,i) - orbit0(j)
           enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL
     endif
     if (jmax .eq. 0 .or. (switch .gt. 1 .and. jmax .lt. j_tot))     &
          goto 20
     if (switch .eq. 2 .and. info) then
        if (mod(turn,100) .eq. 0) print *, 'turn :', turn
     endif

     if(exit_loss_turn .and. lost_in_turn) then
        lost_in_turn = .false.
        write(text, '(i6)') turn
        call aawarn('TRRUN Frozen SC: ',&
             'Stop in Loss Turn: '//text)
        goto 20
     endif

  enddo
  !--- end of loop over turns
20 continue
  if (switch .gt. 1 .and. jmax .lt. j_tot)  then
     e_flag = 1
     return
  endif
  do i = 1, jmax
     last_turn(part_id(i)) = min(turns+tot_turn, turn+tot_turn)
     last_pos(part_id(i)) = sum
     do j = 1, 6
        last_orbit(j,part_id(i)) = z(j,i)
     enddo
  enddo
  turn = min(turn, turns)

  if(bb_sxy_update) then
     tot_turn=tot_turn+turn
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
     do i = 1, jmax
        part_id_keep(i)=part_id(i)
        last_turn_keep(i)=last_turn(i)
        do j=1,6
           z_keep(j,i)=z(j,i)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

  !--- enter last turn in tables if not done already
  if (.not. last_out)  then
     if (switch .eq. 1)  then
        if (onetable)  then
           !hbu
           spos=sum
           !hbu get current node name
           call element_name(el_name,len(el_name))
           !hbu spos added
           call tt_putone(jmax, tot_turn+turn, tot_segm, segment, part_id,      &
                z, orbit0,spos,nlm,el_name)
        else
           do i = 1, jmax
              !hbu
              call tt_puttab(part_id(i), turn, 1, z(1,i), orbit0,spos)
           enddo
        endif
     endif
  endif

!--- Write checkpoint_restart data
  rewind 90
  write(90) jmax
  write(90) Ex_rms
  write(90) Ey_rms
  do i = 1, jmax
     do j=1,6
        write(90) z(j,i)
     enddo
  enddo

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
  goto 101
100 call aafail('TRACK: Fatal ', 'checkpoint_restart file corrupted')
101 continue
end subroutine trrun


subroutine ttmap(switch,code,el,track,ktrack,dxt,dyt,sum,turn,part_id,   &
     last_turn,last_pos,last_orbit,aperflag,maxaper,al_errors, onepass)

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
  !   switch    (integer)   1: RUN, 2: DYNAP fastune                     *
  !   SUM       (double)    Accumulated length.                          *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   NUMBER(*) (integer) Number of current track.                       *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  logical aperflag,fmap,onepass
  integer turn,code,ktrack,part_id(*),last_turn(*),nn,jtrk,         &
       get_option, optiondebug, switch
  double precision apx,apy,apr,el,sum,node_value,track(6,*),        &
       last_pos(*),last_orbit(6,*),parvec(26),get_value,ct,tmp,          &
       aperture(maxnaper),one,maxaper(6),al_errors(align_max),st,  &
       theta,ek(6),re(6,6),te(6,6,6),craporb(6),dxt(*),dyt(*),offset(2), &
       offx,offy,min_double
  character(name_len) aptype
  parameter(one=1d0)
  parameter(min_double = 1.e-36)

  optiondebug = get_option('debug ')

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
  if (theta .ne. 0d0)  then
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

! !  ---- Special colllimator aperture check taken out AK 20071211
! !  ---- Collimator with elliptic aperture.
!        if(code.eq.20) then
!          apx = node_value('xsize ')
!          apy = node_value('ysize ')
!          if(apx.eq.0d0) then
!            apx=maxaper(1)
!          endif
!          if(apy.eq.0d0) then
!            apy=maxaper(3)
!          endif
!          call trcoll(1, apx, apy, apr, turn, sum, part_id, last_turn,    &
!               last_pos, last_orbit, track, ktrack,al_errors,offx,offy)
!          go to 500
!        endif
! !  ---- Collimator with rectangular aperture.
!        if(code.eq.21) then
!          apx = node_value('xsize ')
!          apy = node_value('ysize ')
!          if(apx.eq.0d0) then
!            apx=maxaper(1)
!          endif
!          if(apy.eq.0d0) then
!            apy=maxaper(3)
!          endif
!          call trcoll(2, apx, apy, apr, turn, sum, part_id, last_turn,    &
!               last_pos, last_orbit, track, ktrack,al_errors,offx,offy)
!          go to 500
!        endif


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

     if (optiondebug .ne. 0) then
        print *, " aperture type ",aptype
        print *, "          aperture ", aperture(1),aperture(2),aperture(3),aperture(4)
        print *, "          offsets  ", offx, offy
        print *, " "
     endif

     !------------  ellipse case ----------------------------------
     if(aptype.eq.'ellipse') then
        apx = aperture(1)
        apy = aperture(2)
        apr = 0.0d0
        if (optiondebug .ne. 0) print *, "ellipse: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(1, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track, ktrack,al_errors,offx,offy)

        !------------  circle case ----------------------------------
     else if(aptype.eq.'circle') then
        apx = aperture(1)
!        if(apx.eq.0d0) then     ! This FP comparison was triggering a bug in some conditions...
        if(apx.lt.min_double) then
           if (optiondebug .ne. 0) print *, " zero circle radius ", apx, " replaced by default ", maxaper(1)
           apx = maxaper(1)
        endif
        apy = apx
        apr = 0.0d0
        if (optiondebug .ne. 0) print *, "circle: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(1, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)

        !------------  rectangle case ----------------------------------
     else if(aptype.eq.'rectangle') then
        apx = aperture(1)
        apy = aperture(2)
        apr = 0.0d0
        if (optiondebug .ne. 0) print *, "rectangle: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(2, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)

        !-------------Racetrack type , added by Yipeng SUN 21-10-2008---
     else if(aptype.eq.'racetrack') then
        apx = aperture(1)
        apy = aperture(2)
        apr = aperture(3)
        if (optiondebug .ne. 0) print *, "racetrack: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(4, apx, apy, apr, turn, sum, part_id, last_turn, &
             last_pos,last_orbit,track,ktrack,al_errors,offx,offy)

        !------------  LHC screen case ----------------------------------
     else if(aptype.eq.'lhcscreen') then
        !        print *, "LHC screen start, Xrect= ",
        !       aperture(1),"  Yrect= ",aperture(2),"  Rcirc= ",aperture(3)
        apx = aperture(3)
        apy = aperture(3)
        apr = 0.0d0
        !JMJ!     Making essential changes in AV's absence, 16/7/2003
        !JMJ!     this tests whether the particle is outside the circumscribing
        !JMJ!     circle.
        if (optiondebug .ne. 0) print *, "lhcscreen 1/2: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(1, apx, apy, apr, turn, sum, part_id, last_turn,  &
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
        if(aperture(1).gt.0d0) apx = aperture(1)
        if(aperture(2).gt.0d0) apy = aperture(2)
        if (optiondebug .ne. 0) print *, "lhcscreen 2/2: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(2, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !        print *, "LHC screen end"

        !------------  marguerite case ----------------------------------
     else if(aptype.eq.'marguerite') then
        apx = aperture(1)
        apy = aperture(2)
        apr = 0.0d0
        if (optiondebug .ne. 0) print *, "marguerite: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(3, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)

        !------------  rectellipse case ----------------------------------
     else if(aptype.eq.'rectellipse') then
        !*****         test ellipse
        apx = aperture(3)
        apy = aperture(4)
        apr = 0.0d0
        if (optiondebug .ne. 0) print *, "rectellipse 1/2: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(1, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !*****         test rectangle
        apx = aperture(1)
        apy = aperture(2)
        apr = 0.0d0
        if (optiondebug .ne. 0) print *, "rectellipse 2/2: apx, apy, apr, offx, offy ", &
             apx, apy, apr, offx, offy
        call trcoll(2, apx, apy, apr, turn, sum, part_id, last_turn,  &
             last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
        !       print*, " test apertures"
        !       print*, "      apx=",apx, " apy=",apy," apxell=",apxell,
        !              " apyell=",apyell
        !          call trcoll(3, apx, apy, apr, turn, sum, part_id, last_turn,  &
        !         last_pos, last_orbit, track,ktrack,al_errors)

     endif

  !  else
     !  ---- 2013-May-22  11:35:54  ghislain: this section for backward compatibility
     !       with the MAD8 style definition of R/ECOLLIMATOR with xsize and ysize parameters,
     !       but ONLY IF apertype and aperture parameters were not specified.
     !
     !       The case where apertype and aperture were specified was
     !       already taken into account in preceding section under if(aperflag) where
     !       the eventual xsize and ysize parameters were simply ignored.

     !       Hence the new style definition of collimators
     !       COLL : RCOLLIMATOR, apertype=rectangle, aperture={1.e-3,1.e-3} ;
     !       has precedence over the apertype/aperture old style definition
     !       COLL : RCOLLIMATOR, xsize=1.e-3, ysize=1.e-3 ;

     !  NOTE: this calls for redefinition of the COLLIMATOR element...
     !        and eventually the demise of the old xsize/ysize syntax

     !!---- Collimator with elliptic aperture.
     ! if(code.eq.20) then
     !    apx = node_value('xsize ')
     !    apy = node_value('ysize ')
     !    if(apx.eq.0d0) then
     !       apx=maxaper(1)
     !    endif
     !    if(apy.eq.0d0) then
     !       apy=maxaper(3)
     !    endif
     !    call trcoll(1, apx, apy, apr, turn, sum, part_id, last_turn,    &
     !         last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
     ! endif

     !!---- Collimator with rectangular aperture.
     ! if(code.eq.21) then
     !    apx = node_value('xsize ')
     !    apy = node_value('ysize ')
     !    if(apx.eq.0d0) then
     !       apx=maxaper(1)
     !    endif
     !    if(apy.eq.0d0) then
     !       apy=maxaper(3)
     !    endif
     !    call trcoll(2, apx, apy, apr, turn, sum, part_id, last_turn,    &
     !         last_pos, last_orbit, track,ktrack,al_errors,offx,offy)
     ! endif

  endif
  !----  END OF IF(APERFLAG) ... ELSE ...  ----------------


  !-- switch on element type BUT DRIFT, COLLIMATORS, BEAM_BEAM / 13.03.03
  !-- 500 has been specified at the relevant places in go to below
  !-- code =1 for drift, treated above, go to 500 directly
  !      print *,"   CODE    ",code
  go to ( 500,  20,  30,  40,  50,  60,  70,  80,  90, 100,         &
       110, 120, 130, 140, 150, 160, 170, 180, 190, 500,                 &
       500, 500, 230, 240, 250, 260, 270, 280, 290, 300,   310, 320,     &
  !     330, 500, 350, 360, 370,500,500,400,410,      500, 500, 500, 500), code
  ! Use this line to enable non-linear thin lens
  !     330, 500, 350, 360, 370,500,500,400,410, 420, 500, 500, 500, 500), code
  ! Enable non-linear thin lens and RF-Multipole
        330, 500, 350, 360, 370,500,500,400,410, 420, 430, 500, 500, 500), code
  !
  !---- Make sure that nothing is executed if element is not known
  go to 500
  !
  !---- Bending magnet.
  !---- AL: use thick-dipole element
20 continue  ! RBEND
30 continue  ! SBEND
  call tttdipole(track,ktrack)
  go to 500
  !---- Arbitrary matrix. OBSOLETE, to be kept for go to
40 continue
  call tmarb(.false.,.false.,craporb,fmap,ek,re,te)
  call tttrak(ek,re,track,ktrack)
  go to 500
  !---- Quadrupole. OBSOLETE, to be kept for go to
  !---- AL: not so fast, here is the thick quadrupole
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
  if (switch .eq. 1) then
     call ttrfloss(turn, sum, part_id, last_turn, &
       last_pos,last_orbit,track,ktrack)
  endif
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
  if (theta .ne. 0d0)  then
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
  use time_varfi
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
  logical first,time_var
  integer iord,jtrk,nd,nord,ktrack,i,j,n_ferr,nn,ns,node_fd_errors,   &
       store_no_fd_err,get_option,turn,noisemax,nn1,in,mylen
  double precision     const,curv,dbi,dbr,dipi,dipr,dx,dy,elrad,    &
       pt,px,py,rfac,rpt1,rpt2,rpx1,rpx2,rpy1,rpy2,                      &
       f_errors(0:maxferr),field(2,0:maxmul),vals(2,0:maxmul),           &
       ordinv(maxmul),track(6,*),dxt(*),dyt(*),normal(0:maxmul),         &
       skew(0:maxmul),bvk,node_value,ttt,        &
       npeak(100), nlag(100), ntune(100), temp,noise

  character(name_len) name

  save first,ordinv
  data first / .true. /

  !---- Precompute reciprocals of orders.
  if (first) then
     do iord = 1, maxmul
        ordinv(iord) = 1d0 / dble(iord)
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
  nd = 2 * max(nn, ns, n_ferr/2-1)
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
  time_var = node_value('time_var ') .ne. 0d0
  if(time_var.and.time_var_m) then
     time_var_m_cnt=time_var_m_cnt+1
     time_var_m_lnt=time_var_m_lnt+1
     if(idnint(time_var_m_ind(time_var_m_cnt)).ne.time_var_m_lnt)    then
          call aafail('TTMULT: ', 'wrong index in Table: time_var_mul')
       endif
     call element_name(name,len(name))
     mylen=len_trim(name)
     if(time_var_m_ch(time_var_m_cnt)(:mylen).ne.name(:mylen))       &
          call aafail('TTMULT: ', 'wrong element name in Table: '//     &
          'time_var_mul')
     do i=maxmul,0,-1
        if(abs(myfield(time_var_m_cnt,1,i)).gt.0d0.or.               &
             abs(myfield(time_var_m_cnt,2,i)).gt.0d0) then
           n_ferr=i
           goto 101
        endif
     enddo
101  n_ferr=2*n_ferr+2
     do i=0,(n_ferr-2)/2
        f_errors(i)=myfield(time_var_m_cnt,1,i)
        f_errors(n_ferr/2+i)=myfield(time_var_m_cnt,2,i)
     enddo
     nd = 2 * max(nn, ns, n_ferr/2-1)
     call dcopy(f_errors,field,nd+2)
     n_ferr = store_no_fd_err(f_errors,n_ferr)
  else
     if (n_ferr .gt. 0) then
        call dcopy(f_errors,field,n_ferr)
     endif
  endif
  !-----added FrankS, 10-12-2008
  !nd = 2 * max(nn, ns, n_ferr/2-1)
  !---- Dipole error.
  !      dbr = bvk * field(1,0) / (1d0 + deltas)
  !      dbi = bvk * field(2,0) / (1d0 + deltas)
  dbr = bvk * field(1,0)
  dbi = bvk * field(2,0)
  !---- Nominal dipole strength.
  !      dipr = bvk * vals(1,0) / (1d0 + deltas)
  !      dipi = bvk * vals(2,0) / (1d0 + deltas)
  dipr = bvk * vals(1,0)
  dipi = bvk * vals(2,0)
  !---- Other components and errors.
  nord = 0
  do iord = 1, nd/2
     do j = 1, 2
        !          field(j,iord) = bvk * (vals(j,iord) + field(j,iord))          &
        !     / (1d0 + deltas)
        field(j,iord) = bvk * (vals(j,iord) + field(j,iord))
        if (field(j,iord) .ne. 0d0)  nord = iord
     enddo
  enddo
  !---- Pure dipole: only quadrupole kicks according to lrad.
  if (nord .eq. 0) then
     do jtrk = 1,ktrack
        dxt(jtrk) = 0d0
        dyt(jtrk) = 0d0
     enddo
     !----------- introduction of dipole focusing
     if(elrad.gt.0d0.and.get_option('thin_foc ').eq.1) then
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
     !          dxt(jtrk) = dxt(jtrk) / (1d0 + deltas)
     !          dyt(jtrk) = dyt(jtrk) / (1d0 + deltas)
     !        enddo
     if(elrad.gt.0d0.and.get_option('thin_foc ').eq.1) then
        do jtrk = 1,ktrack
           dxt(jtrk) = dxt(jtrk) + dipr*dipr*track(1,jtrk)/elrad
           dyt(jtrk) = dyt(jtrk) + dipi*dipi*track(3,jtrk)/elrad
        enddo
     endif
  endif

  !---- Radiation loss at entrance.
  if (dorad .and. elrad .ne. 0) then
     const = arad * gammas**3 / 3d0

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
           track(2,jtrk) = px - rfac * (1d0 + pt) * px
           track(4,jtrk) = py - rfac * (1d0 + pt) * py
           track(6,jtrk) = pt - rfac * (1d0 + pt) ** 2
        enddo

        !---- Energy loss like for closed orbit.
     else

        !---- Store energy loss on closed orbit.
        rfac = const * ((dipr + dxt(1))**2 + (dipi + dyt(1))**2)
        rpx1 = rfac * (1d0 + track(6,1)) * track(2,1)
        rpy1 = rfac * (1d0 + track(6,1)) * track(4,1)
        rpt1 = rfac * (1d0 + track(6,1)) ** 2

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
     ttt = sqrt(1d0+2d0*track(6,jtrk)*bet0i+track(6,jtrk)**2)
     !        track(2,jtrk) = track(2,jtrk) -                                 &
     !     (dbr + dxt(jtrk) - dipr * (deltas + beti*track(6,jtrk)))
     !        track(4,jtrk) = track(4,jtrk) +                                 &
     !     (dbi + dyt(jtrk) - dipi * (deltas + beti*track(6,jtrk)))
     !        track(5,jtrk) = track(5,jtrk)                                   &
     !     - (dipr*track(1,jtrk) - dipi*track(3,jtrk)) * beti
     track(2,jtrk) = track(2,jtrk) -                                 &
          (dbr + dxt(jtrk) - dipr * (ttt - 1d0))
     track(4,jtrk) = track(4,jtrk) +                                 &
          (dbi + dyt(jtrk) - dipi * (ttt - 1d0))
     track(5,jtrk) = track(5,jtrk) -                                 &
          (dipr*track(1,jtrk) - dipi*track(3,jtrk)) *                       &
          ((1d0 + bet0*track(6,jtrk))/ttt)*bet0i
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
           track(2,jtrk) = px - rfac * (1d0 + pt) * px
           track(4,jtrk) = py - rfac * (1d0 + pt) * py
           track(6,jtrk) = pt - rfac * (1d0 + pt) ** 2
        enddo

        !---- Energy loss like for closed orbit.
     else

        !---- Store energy loss on closed orbit.
        rfac = const * ((dipr + dxt(1))**2 + (dipi + dyt(1))**2)
        rpx2 = rfac * (1d0 + track(6,1)) * track(2,1)
        rpy2 = rfac * (1d0 + track(6,1)) * track(4,1)
        rpt2 = rfac * (1d0 + track(6,1)) ** 2

        do jtrk = 1,ktrack
           track(2,jtrk) = track(2,jtrk) - rpx2
           track(4,jtrk) = track(4,jtrk) - rpy2
           track(6,jtrk) = track(6,jtrk) - rpt2
        enddo
     endif
  endif

end subroutine ttmult
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
  double precision amean,curv,dlogr,el,frndm,rfac,scalen,scaleu,         &
       slope,ucrit,xi,deltap,hbar,clight,arad,pc,gamma,amass,real_am,    &
       get_value,get_variable,tabxi(maxtab),taby(maxtab),five,twelve,fac1
  parameter(five=5d0,twelve=12d0,fac1=3.256223d0)
  character(20) text
  data (taby(i), i = 1, 52)                                         &
       / -1.14084005d0,  -0.903336763d0, -0.769135833d0, -0.601840854d0, &
       -0.448812515d0, -0.345502228d0, -0.267485678d0, -0.204837948d0,   &
       -0.107647471d0, -0.022640628d0,  0.044112321d0,  0.0842842236d0,  &
       0.132941082d0,  0.169244036d0,  0.196492359d0,  0.230918407d0,    &
       0.261785239d0,  0.289741248d0,  0.322174788d0,  0.351361096d0,    &
       0.383441716d0,  0.412283719d0,  0.442963421d0,  0.472622454d0,    &
       0.503019691d0,  0.53197819d0,   0.561058342d0,  0.588547111d0,    &
       0.613393188d0,  0.636027336d0,  0.675921738d0,  0.710166812d0,    &
       0.725589216d0,  0.753636241d0,  0.778558254d0,  0.811260045d0,    &
       0.830520391d0,  0.856329501d0,  0.879087269d0,  0.905612588d0,    &
       0.928626955d0,  0.948813677d0,  0.970829248d0,  0.989941061d0,    &
       1.0097903d0,    1.02691281d0,   1.04411256d0,   1.06082714d0,     &
       1.0750246d0,    1.08283985d0,   1.0899564d0,    1.09645379d0 /
  data (taby(i), i = 53, 101)                                       &
       /  1.10352755d0,   1.11475027d0,   1.12564385d0,   1.1306442d0,   &
       1.13513422d0,   1.13971806d0,   1.14379156d0,   1.14741969d0,     &
       1.15103698d0,   1.15455759d0,   1.15733826d0,   1.16005647d0,     &
       1.16287541d0,   1.16509759d0,   1.16718769d0,   1.16911888d0,     &
       1.17075884d0,   1.17225218d0,   1.17350936d0,   1.17428589d0,     &
       1.17558432d0,   1.17660713d0,   1.17741513d0,   1.17805469d0,     &
       1.17856193d0,   1.17896497d0,   1.17928565d0,   1.17954147d0,     &
       1.17983139d0,   1.1799767d0,    1.18014216d0,   1.18026078d0,     &
       1.18034601d0,   1.1804074d0,    1.18045175d0,   1.1804837d0,      &
       1.18051291d0,   1.18053186d0,   1.18054426d0,   1.18055236d0,     &
       1.18055761d0,   1.18056166d0,   1.18056381d0,   1.1805656d0,      &
       1.18056655d0,   1.18056703d0,   1.18056726d0,   1.1805675d0,      &
       1.18056762d0 /
  data (tabxi(i), i = 1, 52)                                        &
       / -7.60090017d0,  -6.90775537d0,  -6.50229025d0,  -5.99146461d0,  &
       -5.52146101d0,  -5.20300722d0,  -4.96184492d0,  -4.76768923d0,    &
       -4.46540833d0,  -4.19970512d0,  -3.98998451d0,  -3.86323285d0,    &
       -3.70908213d0,  -3.59356928d0,  -3.50655794d0,  -3.39620972d0,    &
       -3.29683733d0,  -3.20645332d0,  -3.10109282d0,  -3.0057826d0,     &
       -2.9004221d0,   -2.80511189d0,  -2.70306253d0,  -2.60369015d0,    &
       -2.50103593d0,  -2.4024055d0 ,  -2.30258512d0,  -2.20727491d0,    &
       -2.12026358d0,  -2.04022098d0,  -1.89712d0   ,  -1.7719568d0,     &
       -1.71479833d0,  -1.60943794d0,  -1.51412773d0,  -1.38629436d0,    &
       -1.30933332d0,  -1.20397282d0,  -1.10866261d0,  -0.99425226d0,    &
       -0.89159810d0,  -0.79850775d0,  -0.69314718d0,  -0.59783697d0,    &
       -0.49429631d0,  -0.40047753d0,  -0.30110508d0,  -0.19845095d0,    &
       -0.10536054d0,  -0.05129330d0,   0.0d0,          0.048790119d0 /
  data (tabxi(i), i = 53, 101)                                      &
       /  0.104360029d0,  0.198850885d0,  0.300104618d0,  0.350656837d0, &
       0.398776114d0,  0.451075643d0,  0.500775278d0,  0.548121393d0,    &
       0.598836541d0,  0.652325153d0,  0.69813472d0 ,  0.746687889d0,    &
       0.802001595d0,  0.850150883d0,  0.900161386d0,  0.951657832d0,    &
       1.00063193d0,   1.05082154d0,   1.09861231d0,   1.13140213d0,     &
       1.1939224d0,    1.25276291d0,   1.3083328d0,    1.36097658d0,     &
       1.4109869d0,    1.45861506d0,   1.50407743d0,   1.54756248d0,     &
       1.60943794d0,   1.64865863d0,   1.70474803d0,   1.75785792d0,     &
       1.80828881d0,   1.85629797d0,   1.90210748d0,   1.9459101d0,      &
       2.0014801d0,    2.05412364d0,   2.10413408d0,   2.15176225d0,     &
       2.19722462d0,   2.25129175d0,   2.29253483d0,   2.35137534d0,     &
       2.40694523d0,   2.45100522d0,   2.501436d0,     2.60268974d0,     &
       2.64617491d0 /

  !Get constants
  clight = get_variable('clight ')
  hbar   = get_variable('hbar ')
  arad   = get_value('probe ','arad ')
  pc     = get_value('probe ','pc ')
  amass  = get_value('probe ','mass ')
  gamma  = get_value('probe ','gamma ')
  deltap = get_value('probe ','deltap ')
  scalen = five / (twelve * hbar * clight)
  scaleu = hbar * 3d0 * clight / 2d0

  !---- AMEAN is the average number of photons emitted.,
  !     NPHOT is the integer number generated from Poisson's law.
  amean = scalen * abs(arad*pc*(1d0+deltap)*el*curv) * sqrt(3d0)
  rfac = 0d0
  real_am = amean
  if (real_am .gt. 0d0) then
     call dpoissn(real_am, nphot, ierror)

     if (ierror .ne. 0) then
        write(text, '(1p,d20.12)') amean
        call aafail('TRPHOT: ','Fatal: Poisson input mean =' // text)
     endif

     !---- For all photons, sum the radiated photon energy,
     !     in units of UCRIT (relative to total energy).

     if (nphot .ne. 0) then
        ucrit = scaleu * gamma**2 * abs(curv) / amass
        xi = 0d0
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
end subroutine trphot

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
  !   EL        (double)    Length of drift.                             *
  !----------------------------------------------------------------------*
  integer itrack,ktrack
  double precision el,pt,px,py,track(6,*),ttt

  ! picked from trturn in madx.ss
!$OMP PARALLEL PRIVATE(itrack, px, py, pt, ttt)
!$OMP DO
  do  itrack = 1, ktrack
     px = track(2,itrack)
     py = track(4,itrack)
     pt = track(6,itrack)
     ttt = el/sqrt(1d0+2d0*pt*bet0i+pt**2 - px**2 - py**2)
     track(1,itrack) = track(1,itrack) + ttt*px
     track(3,itrack) = track(3,itrack) + ttt*py
     !        track(5,itrack) = track(5,itrack)                               &
     !     + el*(beti + pt * dtbyds) - (beti+pt)*ttt
     !---- AK 20060413
     !---- Ripken DESY-95-189 p.36
     track(5,itrack) = track(5,itrack)                               &
          + bet0i*(el - (1d0 + bet0*pt) * ttt)
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine ttdrf
subroutine ttrf(track,ktrack)

  use twtrrfi
  use name_lenfi
  use time_varfi
  use trackfi
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
  logical time_var
  integer itrack,ktrack,mylen
  double precision omega,phirf,pt,rff,rfl,rfv,track(6,*),clight,    &
       twopi,vrf,pc0,get_variable,node_value,get_value,bvk,    &
       ten3m,ten6p
  !      double precision px,py,ttt,beti,el1
  parameter(ten3m=1d-3,ten6p=1d6)
  character*(name_len) name

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- BV flag
  bvk = node_value('other_bv ')

  !---- Fetch data.
  !      el = node_value('l ')
  !      el1 = node_value('l ')
  rfv = bvk * node_value('volt ')
  time_var = node_value('time_var ') .ne. 0d0
  if(time_var.and.time_var_c) then
     time_var_c_cnt=time_var_c_cnt+1
     time_var_c_lnt=time_var_c_lnt+1
     if(idnint(time_var_c_ind(time_var_c_cnt)).ne.time_var_c_lnt)    &
          call aafail('TMARB: ', 'wrong index in Table: time_var_cav')
     call element_name(name,len(name))
     mylen=len_trim(name)
     if(time_var_c_ch(time_var_c_cnt)(:mylen).ne.name(:mylen))       &
          call aafail('TMARB: ', 'wrong element name in Table: '//   &
          'time_var_cav')
!---- Overwrite Volt
     rfv = cav_volt(time_var_c_cnt)
!---- Store Volt
     call store_node_value('volt ',rfv)
  endif

!      rfv = node_value('volt ')
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
  !        el1 = el / 2d0
  !        rfv = rfv / 2d0
  !        lwake = .true.
  !      else
  !        el1 = el
  !        lwake = .false.
  !      endif
  !---- Set up.
  omega = rff * (ten6p * twopi / clight)
  !      vrf   = rfv * ten3m / (pc * (1d0 + deltas))
  vrf   = rfv * ten3m / pc0
  phirf = rfl * twopi
  !      dl    = el * 5d-1
  !      bi2gi2 = 1d0 / (betas * gammas) ** 2
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
  !        call ttwake(2d0*el1, nbin, binmax, lfile, tfile, ener1, track,
  !     +              ktrack)
  !
  !*---- Track 2nd half of cavity -- loop for all particles.
  !      do 20 itrack = 1, ktrack
  !
  !*---- Drift to centre.
  !         px = track(2,itrack)
  !         py = track(4,itrack)
  !         pt = track(6,itrack)
  !         ttt = 1d0/sqrt(1d0+2d0*pt*beti+pt**2 - px**2 - py**2)
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
  !         ttt = 1d0/sqrt(1d0+2d0*pt*beti+pt**2 - px**2 - py**2)
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
       get_variable,node_value,get_value,ten3m,ten6p,px
  !      double precision px,py,ttt,beti,el1
  parameter(ten3m=1d-3,ten6p=1d6)

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
     vrf = 0d0
  else if (turn .ge. t1 .and. turn .lt. t2) then
     vrf = (turn-t1)*rfv*ten3m/pc0/(t2-t1)
  else if (turn .ge. t2 .and. turn .lt. t3) then
     vrf = rfv*ten3m/pc0
  else if (turn .ge. t3 .and. turn .lt. t4) then
     vrf = (t4-turn)*rfv*ten3m/pc0/(t4-t3)
  else
     vrf = 0d0
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
subroutine tthacdip(track,ktrack,turn)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of trajectories through a thin h ac dipole (zero l.)   *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  ! Added by Yipeng SUN on 11 Nov 2009                                   *
  !----------------------------------------------------------------------*
  integer itrack,ktrack,turn,turn1,turn2,turn3,turn4
  double precision omega,phirf,rff,rfl,rfv,track(6,*),clight,twopi,vrf, &
      pc0,get_variable,node_value,get_value,ten3m,ten6p,px
  !      double precision px,py,ttt,beti,el1
  parameter(ten3m=1d-3,ten6p=1d6)

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- Fetch data.
  rfv = node_value('volt ')
  rff = node_value('freq ')
  rfl = node_value('lag ')
  pc0 = get_value('beam ','pc ')

  turn1 = node_value('ramp1 ')
  turn2 = node_value('ramp2 ')
  turn3 = node_value('ramp3 ')
  turn4 = node_value('ramp4 ')
  !---- Set up.
  omega = rff * twopi
  vrf   = 300 * rfv * ten3m / pc0
  phirf = rfl * twopi

  if (turn .lt. turn1)  then
     vrf = 0
  else if (turn .ge. turn1 .and. turn .lt. turn2) then
     vrf = (turn-turn1) * vrf / (turn2-turn1)
  else if (turn .ge. turn2 .and. turn .lt. turn3) then
     vrf = vrf
  else if (turn .ge. turn3 .and. turn .lt. turn4) then
     vrf = (turn4-turn) * vrf / (turn4-turn3)
  else
     vrf = 0
  endif
  !      if (turn .le. 10 .or. turn .gt. 1990) then
  !       print*," turn: ",turn, " vrf: ", vrf
  !      endif
  do itrack = 1, ktrack
     px  = track(2,itrack)                                           &
          + vrf * sin(phirf + omega * turn)

     !---- track(2,jtrk) = track(2,jtrk)
     !        pt = track(6,itrack)
     !        pt = pt + vrf * sin(phirf - omega * track(5,itrack))

     track(2,itrack) = px
  enddo
end subroutine tthacdip
subroutine ttvacdip(track,ktrack,turn)

  implicit none

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Track a set of trajectories through a thin v ac dipole (zero l.)   *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  ! Added by Yipeng SUN on 11 Nov 2009                                   *
  !----------------------------------------------------------------------*
  integer itrack,ktrack,turn,turn1,turn2,turn3,turn4
  double precision omega,phirf,rff,rfl,rfv,track(6,*),clight,twopi,vrf, &
       pc0,get_variable,node_value,get_value,ten3m,ten6p,py
  !      double precision px,py,ttt,beti,el1
  parameter(ten3m=1d-3,ten6p=1d6)

  !---- Initialize
  clight=get_variable('clight ')
  twopi=get_variable('twopi ')

  !---- Fetch data.
  rfv = node_value('volt ')
  rff = node_value('freq ')
  rfl = node_value('lag ')
  pc0 = get_value('beam ','pc ')

  turn1 = node_value('ramp1 ')
  turn2 = node_value('ramp2 ')
  turn3 = node_value('ramp3 ')
  turn4 = node_value('ramp4 ')
  !---- Set up.
  omega = rff * twopi
  vrf   = 300 * rfv * ten3m / pc0
  phirf = rfl * twopi


  if (turn .lt. turn1)  then
     vrf = 0
  else if (turn .ge. turn1 .and. turn .lt. turn2) then
     vrf = (turn-turn1) * vrf / (turn2-turn1)
  else if (turn .ge. turn2 .and. turn .lt. turn3) then
     vrf = vrf
  else if (turn .ge. turn3 .and. turn .lt. turn4) then
     vrf = (turn4-turn) * vrf / (turn4-turn3)
  else
     vrf = 0
  endif
  !      if (turn .le. 10 .or. turn .gt. 1990) then
  !       print*," turn: ",turn, " vrf: ", vrf
  !      endif
  do itrack = 1, ktrack
     py  = track(4,itrack)                                           &
          + vrf * sin(phirf + omega * turn)

     !---- track(2,jtrk) = track(2,jtrk)
     !        pt = track(6,itrack)
     !        pt = pt + vrf * sin(phirf - omega * track(5,itrack))

     track(4,itrack) = py
  enddo
end subroutine ttvacdip
!FIXME Unused dummy argument 'el'
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
  parameter(one=1d0,ten3m=1d-3)
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
     deltap=sqrt(1d0-1d0/beta/beta+(pt+1d0/beta)**2) - 1d0
     kick0=charge*ten3m/pc/(1d0+deltap)/beta
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
       node_value,twopi
  double precision dpxx,dpyy
  double precision temp,sinpeak,sintune,sinphase

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
  dorad = get_value('probe ','radiate ') .ne. 0d0
  dodamp = get_option('damp ') .ne. 0
  dorand = get_option('quantum ') .ne. 0
  if (el .eq. 0d0)  then
     div = 1d0
  else
     div = el
  endif

  do i = 1, 2
     field(i) = 0d0
  enddo
  if (n_ferr .gt. 0) call dcopy(f_errors, field, min(2, n_ferr))
  if (code.eq.14) then
     xkick=bvk*(node_value('kick ')+node_value('chkick ')+           &
          field(1)/div)
     ykick=0d0
  else if (code.eq.15) then
     xkick=bvk*(node_value('hkick ')+node_value('chkick ')+          &
          field(1)/div)
     ykick=bvk*(node_value('vkick ')+node_value('cvkick ')+          &
          field(2)/div)
  else if(code.eq.16) then
     xkick=0d0
     ykick=bvk*(node_value('kick ')+node_value('cvkick ')+           &
          field(2)/div)
  else
     xkick=0d0
     ykick=0d0
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
  dpx = xkick / (1d0 + deltas)
  dpy = ykick / (1d0 + deltas)
  !---- Thin lens 6D version
  dpxx = xkick
  dpyy = ykick

  bil2 = el / (2d0 * betas)
  bi2gi2 = 1d0 / (betas * gammas) ** 2

  !---- Half radiation effects at entrance.
  if (dorad  .and.  el .ne. 0) then
     if (dodamp .and. dorand) then
        curv = sqrt(dpx**2 + dpy**2) / el
     else
        rfac = arad * gammas**3 * (dpx**2 + dpy**2) / (3d0 * el)
     endif

     !---- Full damping.
     if (dodamp) then
        do itrack = 1, ktrack
           if (dorand) call trphot(el,curv,rfac,deltas)
           px = track(2,itrack)
           py = track(4,itrack)
           pt = track(6,itrack)
           track(2,itrack) = px - rfac * (1d0 + pt) * px
           track(4,itrack) = py - rfac * (1d0 + pt) * py
           track(6,itrack) = pt - rfac * (1d0 + pt) ** 2
        enddo

        !---- Energy loss as for closed orbit.
     else
        rpx = rfac * (1d0 + track(6,1)) * track(2,1)
        rpy = rfac * (1d0 + track(6,1)) * track(4,1)
        rpt = rfac * (1d0 + track(6,1)) ** 2

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
  !        px = track(2,itrack) + 5d-1 * dpx
  !        py = track(4,itrack) + 5d-1 * dpy
  !        pt = track(6,itrack)
  !
  !!---- Drift through corrector.
  !        d = (1d0 - pt / betas) * el
  !        track(1,itrack) = track(1,itrack) + px * d
  !        track(3,itrack) = track(3,itrack) + py * d
  !        track(5,itrack) = track(5,itrack) + d * bi2gi2 * pt -           &
  !     bil2 * (px**2 + py**2 + bi2gi2*pt**2) + el*pt*dtbyds
  !
  !!---- Half kick at exit.
  !        track(2,itrack) = px + 5d-1 * dpx
  !        track(4,itrack) = py + 5d-1 * dpy
  !        track(6,itrack) = pt
  !      enddo

  !---- Kick at dipole corrector magnet
  !     including PT-dependence
  do itrack = 1, ktrack
     px = track(2,itrack)
     py = track(4,itrack)
     pt = track(6,itrack)
     !        ttt = sqrt(1d0+2d0*pt*bet0i+pt**2 - px**2 - py**2)
     !        ddd = sqrt(1d0+2d0*pt*bet0i+pt**2)
     !        track(2,itrack) = px + dpxx*ttt/ddd
     !        track(4,itrack) = py + dpyy*ttt/ddd
     !        ttt = sqrt(1d0+2d0*pt*bet0i+pt**2 - dpxx**2 - dpyy**2)
     !        ddd = sqrt(1d0+2d0*pt*bet0i+pt**2)
     !        track(2,itrack) = px + dpxx*ttt/ddd
     !        track(4,itrack) = py + dpyy*ttt/ddd

     !        ttt = sqrt(1d0+2d0*pt*bet0i+pt**2 - px**2 - py**2)
     !        ddd = sqrt(1d0+2d0*pt*bet0i+pt**2)
     !
     !        xp = px/ttt + dpxx/ddd  ! Apply kick in standard coordinates!
     !        yp = py/ttt + dpyy/ddd
     !
     !        ttt = sqrt(1d0 + xp**2 + yp**2)
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
     !        ((1d0 + bet0*track(6,itrack))/ddd)*bet0i

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
           track(2,itrack) = px - rfac * (1d0 + pt) * px
           track(4,itrack) = py - rfac * (1d0 + pt) * py
           track(6,itrack) = pt - rfac * (1d0 + pt) ** 2
        enddo

        !---- Energy loss as for closed orbit.
     else
        rpx = rfac * (1d0 + track(6,1)) * track(2,1)
        rpy = rfac * (1d0 + track(6,1)) * track(4,1)
        rpt = rfac * (1d0 + track(6,1)) ** 2

        do itrack = 1, ktrack
           track(2,itrack) = track(2,itrack) - rpx
           track(4,itrack) = track(4,itrack) - rpy
           track(6,itrack) = track(6,itrack) - rpt
        enddo
     endif
  endif

end subroutine ttcorr
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
  double precision amu,ran,pir,frndm,grndm,expma,amax
  !    AMAX IS THE VALUE ABOVE WHICH THE NORMAL DISTRIBUTION MUST BE USED
  data amax/88d0/

  ierror= 0
  if(amu.le.0d0) then
     !    MEAN SHOULD BE POSITIVE
     ierror=1
     n = 0
     go to 999
  endif
  if(amu.gt.amax) then
     !   NORMAL APPROXIMATION FOR AMU.GT.AMAX
     ran = grndm()
     n=ran*sqrt(amu)+amu+5d-1
     goto 999
  endif
  expma=exp(-amu)
  pir=1d0
  n=-1
10 n=n+1
  pir=pir*frndm()
  if(pir.gt.expma) go to 10
999 end subroutine dpoissn
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
  integer ktrack,beamshape,b_dir_int,get_option
  double precision track(6,*),parvec(*),fk,dp
  double precision gamma0,beta0,beta_dp,ptot,b_dir
  double precision q,q_prime,get_value,node_value,get_variable
  double precision cme32
  logical first,bb_ultra_relati
  save first
  data first / .true. /
  parameter(cme32=1d-32)
  !     if x > explim, exp(-x) is outside machine limits.

  !---- Calculate momentum deviation and according changes
  !     of the relativistic factor beta0
  dp  = get_variable('track_deltap ')
  q = get_value('probe ','charge ')
  q_prime = node_value('charge ')
  gamma0 = parvec(7)
  beta0 = sqrt(1d0-1d0/gamma0**2)
  ptot = beta0*gamma0*(1d0+dp)
  beta_dp = ptot / sqrt(1d0 + ptot**2)
  b_dir_int = node_value('bbdir ')
  b_dir=dble(b_dir_int)
  b_dir = b_dir/sqrt(b_dir*b_dir + cme32)
  !---- pre-factor, if zero, anything else does not need to be calculated
  bb_ultra_relati = get_option('bb_ultra_relati ') .ne. 0
  if(bb_ultra_relati) then
     fk = 2d0 * parvec(5) * parvec(6) / parvec(7)
  else
     fk = 2d0*parvec(5)*parvec(6)/parvec(7)/beta0/(1d0+dp)/q*          &
       (1d0-beta0*beta_dp*b_dir)/(beta_dp+5d-1*(b_dir-1d0)*b_dir*beta0)
  endif
  !
  if (fk .eq. 0d0)  return
  !---- choose beamshape: 1-Gaussian, 2-flattop=trapezoidal, 3-hollow-parabolic
  beamshape = node_value('bbshape ')
  if(beamshape.lt.1.or.beamshape.gt.3) then
     beamshape=1
     if(first) then
        first = .false.
        call aawarn('TTBB: ',                                         &
             'beamshape out of range, set to default=1')
     endif
  endif
  if(beamshape.eq.1) call ttbb_gauss(track,ktrack,fk)
  if(beamshape.eq.2) call ttbb_flattop(track,ktrack,fk)
  if(beamshape.eq.3) call ttbb_hollowparabolic(track,ktrack,fk)

end subroutine ttbb
subroutine ttbb_gauss(track,ktrack,fk)

  use name_lenfi
  use bbfi
  use spch_bbfi
  use fasterror
  implicit none

  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having a Gaussian perpendicular shape                       *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  logical bborbit,bb_sxy_update
  integer ktrack,itrack,ipos,get_option,mylen
  double precision track(6,*),pi,sx,sy,xm,ym,sx2,sy2,xs,                 &
       ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,cbx,cby,      &
       get_variable,node_value,ten3m,explim
  double precision xrv(ktrack), yrv(ktrack), crxv(ktrack), cryv(ktrack), &
       xbv(ktrack), ybv(ktrack), cbxv(ktrack), cbyv(ktrack),tkv(ktrack), &
       phixv(ktrack),phiyv(ktrack),xsv(ktrack),ysv(ktrack),rkv(ktrack)
  parameter(ten3m=1d-3,explim=150d0)
 !VVK 20100321 ------------------------------------------------------
      real  gauss_factor_t
      character*20 text
      character*(name_len) name
!-------------------------------------------------------------------

 !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  pi=get_variable('pi ')
  bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
  if(bb_sxy_update) then
     name=' '
     call element_name(name,len(name))
     mylen=len_trim(name)
     i_spch=i_spch+1
     if(i_spch.gt.N_spch) then
        write(text, '(1p,i8)') i_spch
        call aafail('TTBB: ', 'Table with too few BB elements: '//  &
             text)
     endif
     if(spch_bb_name(i_spch)(:mylen).ne.name(:mylen)) then
        call aafail('TTBB: ', 'wrong element name in Table: '//     &
             'spch_bb')
     endif

     sx=sqrt(betx_bb(i_spch)*Ex_rms+(dx_bb(i_spch)*sigma_p)**2)
     sy=sqrt(bety_bb(i_spch)*Ey_rms+(dy_bb(i_spch)*sigma_p)**2)
  else
     sx = node_value('sigx ')
     sy = node_value('sigy ')
  endif

  xm = node_value('xma ')
  ym = node_value('yma ')
  if(bb_sxy_update) fk = fk * rat_bb_n_ions !Ratio_for_bb_N_ions
  if (fk .eq. 0d0)  return
  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif
  sx2 = sx*sx
  sy2 = sy*sy
  !---- limit formulae for sigma(x) = sigma(y).
  if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
!$OMP PARALLEL PRIVATE(itrack, xs, ys, rho2, tk, gauss_factor_t, phix, phiy)
!$OMP DO
     do itrack = 1, ktrack
        xs = track(1,itrack) - xm
        ys = track(3,itrack) - ym
        rho2 = xs * xs + ys * ys
        tk = rho2 / (2d0 * sx2)

        if(bb_sxy_update) then
           gauss_factor_t=                                 &!VVK 20100321
                exp(-5d-1*(track(5,itrack)-mean_t)**2/sigma_t**2)!VVK 20100321

           if (tk .gt. explim) then
              phix = xs * fk / rho2 &
                   *gauss_factor_t !VVK 20100321
              phiy = ys * fk / rho2 &
                   *gauss_factor_t !VVK 20100321
           else if (rho2 .ne. 0d0) then
              phix = xs * fk / rho2 * (1d0 - exp(-tk) ) &
                   *gauss_factor_t !VVK 20100321
              phiy = ys * fk / rho2 * (1d0 - exp(-tk) ) &
                   *gauss_factor_t !VVK 20100321
           else
              phix = 0d0
              phiy = 0d0
           endif
        else
           if (tk .gt. explim) then
              phix = xs * fk / rho2
              phiy = ys * fk / rho2
           else if (rho2 .ne. 0d0) then
              phix = xs * fk / rho2 * (1d0 - exp(-tk) )
              phiy = ys * fk / rho2 * (1d0 - exp(-tk) )
           else
              phix = 0d0
              phiy = 0d0
           endif
        endif

        if (ipos .ne. 0)  then
           !--- subtract closed orbit kick
           phix = phix - bb_kick(1,ipos)
           phiy = phiy - bb_kick(2,ipos)
        endif
        track(2,itrack) = track(2,itrack) + phix
        track(4,itrack) = track(4,itrack) + phiy
     enddo
!$OMP END DO
!$OMP END PARALLEL

     !---- case sigma(x) > sigma(y).
  else if (sx2 .gt. sy2) then
     r2 = 2d0 * (sx2 - sy2)
     r  = sqrt(r2)
 !        rk = fk * sqrt(pi) / r                 !VVK 20100321
     rk = fk * sqrt(pi) / r
     if(fasterror_on) then
!$OMP PARALLEL PRIVATE(itrack, gauss_factor_t)
!$OMP DO
        do itrack = 1, ktrack

           if(bb_sxy_update) then
              gauss_factor_t= &                                !VVK 20100321
                   exp(-5d-1*(track(5,itrack)-mean_t)**2/sigma_t**2)!VVK 20100321
              rkv(itrack) = fk * sqrt(pi) / r*gauss_factor_t !VVK 20100321
           endif

           xsv(itrack) = track(1,itrack) - xm
           ysv(itrack) = track(3,itrack) - ym
           xrv(itrack) = abs(xsv(itrack)) / r
           yrv(itrack) = abs(ysv(itrack)) / r
           tkv(itrack) = (xsv(itrack) * xsv(itrack) / sx2 + ysv(itrack) * ysv(itrack) / sy2) / 2d0
           xbv(itrack) = (sy / sx) * xrv(itrack)
           ybv(itrack) = (sx / sy) * yrv(itrack)
        enddo
!$OMP END DO
!$OMP END PARALLEL
        call wzsubv(ktrack,xrv, yrv, crxv, cryv)
        call wzsubv(ktrack,xbv, ybv, cbxv, cbyv)
!        do itrack = 1, ktrack
!           if (tkv(itrack) .gt. explim) then
!              phixv(itrack) = rk * cryv(itrack)
!              phiyv(itrack) = rk * crxv(itrack)
!           endif
!        enddo
!$OMP PARALLEL PRIVATE(itrack)
!$OMP DO
        do itrack = 1, ktrack
           phixv(itrack) = rkv(itrack) * (cryv(itrack) - exp(-tkv(itrack)) * cbyv(itrack))
           phiyv(itrack) = rkv(itrack) * (crxv(itrack) - exp(-tkv(itrack)) * cbxv(itrack))
           track(2,itrack) = track(2,itrack) + phixv(itrack) * sign(1d0,xsv(itrack))
           track(4,itrack) = track(4,itrack) + phiyv(itrack) * sign(1d0,ysv(itrack))
           if (ipos .ne. 0)  then
              !--- subtract closed orbit kick
              track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
              track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
           endif
        enddo
!$OMP END DO
!$OMP END PARALLEL
     else
!$OMP PARALLEL PRIVATE(itrack, rk, xs, ys, xr, yr, tk, xb, yb, crx, cry, cbx, cby, gauss_factor_t, phix, phiy)
!$OMP DO
        do itrack = 1, ktrack

           if(bb_sxy_update) then
              gauss_factor_t= &                                !VVK 20100321
                   exp(-5d-1*(track(5,itrack)-mean_t)**2/sigma_t**2)!VVK 20100321
              rk = fk * sqrt(pi) / r*gauss_factor_t !VVK 20100321
           endif

           xs = track(1,itrack) - xm
           ys = track(3,itrack) - ym
           xr = abs(xs) / r
           yr = abs(ys) / r
           call ccperrf(xr, yr, crx, cry)
           tk = (xs * xs / sx2 + ys * ys / sy2) / 2d0
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
           track(2,itrack) = track(2,itrack) + phix * sign(1d0,xs)
           track(4,itrack) = track(4,itrack) + phiy * sign(1d0,ys)
           if (ipos .ne. 0)  then
              !--- subtract closed orbit kick
              track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
              track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
           endif
        enddo
!$OMP END DO
!$OMP END PARALLEL

     endif

     !---- case sigma(x) < sigma(y).
  else
     r2 = 2d0 * (sy2 - sx2)
     r  = sqrt(r2)
!        rk = fk * sqrt(pi) / r                 !VVK 20100321
     rk = fk * sqrt(pi) / r
!$OMP PARALLEL PRIVATE(itrack, rk, xs, ys, xr, yr, tk, xb, yb, crx, cry, cbx, cby, gauss_factor_t, phix, phiy)
!$OMP DO
     do itrack = 1, ktrack

        if(bb_sxy_update) then
           gauss_factor_t= &                                !VVK 20100321
                exp(-5d-1*(track(5,itrack)-mean_t)**2/sigma_t**2)!VVK 20100321
           rk = fk * sqrt(pi) / r*gauss_factor_t !VVK 20100321
        endif

        xs = track(1,itrack) - xm
        ys = track(3,itrack) - ym
        xr = abs(xs) / r
        yr = abs(ys) / r

        if(fasterror_on) then
           call wzsub(yr, xr, cry, crx)
        else
           call ccperrf(yr, xr, cry, crx)
        endif

        tk = (xs * xs / sx2 + ys * ys / sy2) / 2d0
        if (tk .gt. explim) then
           phix = rk * cry
           phiy = rk * crx
        else
           xb  = (sy / sx) * xr
           yb  = (sx / sy) * yr

           if(fasterror_on) then
              call wzsub(yb, xb, cby, cbx)
           else
              call ccperrf(yb, xb, cby, cbx)
           endif

           phix = rk * (cry - exp(-tk) * cby)
           phiy = rk * (crx - exp(-tk) * cbx)
        endif
        track(2,itrack) = track(2,itrack) + phix * sign(1d0,xs)
        track(4,itrack) = track(4,itrack) + phiy * sign(1d0,ys)
        if (ipos .ne. 0)  then
           !--- subtract closed orbit kick
           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
        endif
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif
end subroutine ttbb_gauss
subroutine ttbb_flattop(track,ktrack,fk)

  use bbfi
  implicit none

  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having an trapezoidal and, so flat top radial profile       *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  logical bborbit,first
  integer ktrack,itrack,ipos,get_option
  double precision track(6,*),pi,r0x,r0y,wi,wx,wy,xm,ym,                   &
       r0x2,r0y2,xs,ys,rho,rho2,fk,phir,phix,phiy,get_variable,            &
       node_value,four,six,eight,twelve,twentyfour,     &
       fortyeight,quarter,ten3m,explim,norm,r1,zz
  parameter(four=4d0,six=6d0,eight=8d0,                                    &
       twelve=12d0,twentyfour=24d0,fortyeight=48d0,quarter=25d-2,          &
       ten3m=1d-3,explim=150d0)
  save first
  data first / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  pi=get_variable('pi ')
  ! mean radii of the is given via variables sigx and sigy
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  wi = node_value('width ')
  xm = node_value('xma ')
  ym = node_value('yma ')
  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif
  r0x2 = r0x*r0x
  r0y2 = r0y*r0y
  wx = r0x*wi
  wy = r0y*wi
  !---- limit formulae for mean radius(x) = mean radius(y),
  !-----      preliminary the only case considered.
  !
  if (abs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
     zz = 5d-1*(r0x + r0y)
     r0x=zz
     r0y=zz
     r0x2=r0x*r0x
     r0y2=r0y*r0y
     if(first) then
        first=.false.
        call aawarn('TTBB_FLATTOP: ','beam is assumed to be circular')
     endif
  endif
  norm = (twelve*r0x**2 + wx**2)/twentyfour
  r1=r0x-wx/2d0
!$OMP PARALLEL PRIVATE(itrack, xs, ys, rho, rho2, phir, phix, phiy)
!$OMP DO
  do itrack = 1, ktrack
     xs = track(1,itrack) - xm
     ys = track(3,itrack) - ym
     rho2 = xs * xs + ys * ys
     rho  = sqrt(rho2)
     if(rho.le.r1) then
        phir = 5d-1/norm
        phix = phir*xs
        phiy = phir*ys
     else if(rho.gt.r1.and.rho.lt.r1+wx) then
        phir = ((r0x**2/four - r0x**3/six/wx - r0x*wx/eight +            &
             wx**2/fortyeight)/rho2 + quarter + 5d-1*r0x/wx -            &
             rho/3d0/wx)/norm
        phix = phir*xs
        phiy = phir*ys
     else if(rho.ge.r1+wx) then
        phir = 1d0/rho2
        phix = xs*phir
        phiy = ys*phir
     endif
     track(2,itrack) = track(2,itrack)+phix*fk
     track(4,itrack) = track(4,itrack)+phiy*fk
  end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine ttbb_flattop
subroutine ttbb_hollowparabolic(track,ktrack,fk)

  use bbfi
  implicit none

  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having a hollow-parabolic perpendicular shape               *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  logical bborbit,first
  integer ktrack,itrack,ipos,get_option
  double precision track(6,*),pi,r0x,r0y,wi,wx,wy,xm,ym,                   &
       r0x2,r0y2,xs,ys,rho,rho2,fk,phir,phix,phiy,get_variable,            &
       node_value,four,twelve,ten3m,  &
       explim,zz
  parameter(four=4d0,twelve=12d0,ten3m=1d-3,explim=150d0)
  save first
  data first / .true. /
  !     if x > explim, exp(-x) is outside machine limits.

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  pi=get_variable('pi ')
  ! mean radii of the is given via variables sigx and sigy
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  wi = node_value('width ')
  ! width is given as FWHM of parabolic density profile, but formulas were
  ! derived with half width at the bottom of the parabolic density profile
  wi = wi/sqrt(2d0)
  xm = node_value('xma ')
  ym = node_value('yma ')
  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif
  r0x2 = r0x*r0x
  r0y2 = r0y*r0y
  wx  = wi*r0x
  wy  = wi*r0y
  !---- limit formulae for mean radius(x) = mean radius(y),
  !-----      preliminary the only case considered.
  !
  if (abs(r0x2 - r0y2) .gt. ten3m * (r0x2 + r0y2)) then
     zz = 5d-1*(r0x + r0y)
     r0x=zz
     r0y=zz
     r0x2=r0x*r0x
     r0y2=r0y*r0y
     if(first) then
        first=.false.
        call aawarn('TTBB_HOLLOWPARABOLIC: ',                         &
             'beam is assumed to be circular')
     endif
  endif
!$OMP PARALLEL PRIVATE(itrack, xs, ys, rho, rho2, phir, phix, phiy)
!$OMP DO
  do itrack = 1, ktrack
     xs = track(1,itrack) - xm
     ys = track(3,itrack) - ym
     rho2 = xs * xs + ys * ys
     rho  = sqrt(rho2)
     if(rho.le.r0x-wx) then
        phix = 0d0
        phiy = 0d0
     else if(rho.gt.r0x-wx.and.rho.lt.r0x+wx) then
        phir=75d-2/wx/r0x/rho2*(r0x**4/twelve/wx**2 - r0x**2/2d0 +       &
             2d0*r0x*wx/3d0 - wx**2/four + rho2/2d0*(1d0 -                    &
             r0x**2/wx**2) + rho**3/3d0*2d0*r0x/wx**2 -                       &
             rho**4/four/wx**2)
        phix = phir*xs
        phiy = phir*ys
     else
        phir = 1d0/rho2
        phix = xs*phir
        phiy = ys*phir
     endif
     track(2,itrack) = track(2,itrack)+phix*fk
     track(4,itrack) = track(4,itrack)+phiy*fk
  end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine ttbb_hollowparabolic
!!$subroutine ttbb_old(track,ktrack,parvec)
!!$
!!$  use bbfi
!!$  implicit none
!!$
!!$  !----------------------------------------------------------------------*
!!$  ! purpose:                                                             *
!!$  !   track a set of particle through a beam-beam interaction region.    *
!!$  !   see mad physicist's manual for the formulas used.                  *
!!$  !input:                                                                *
!!$  ! input/output:                                                        *
!!$  !   track(6,*)(double)  track coordinates: (x, px, y, py, t, pt).      *
!!$  !   ktrack    (integer) number of tracks.                              *
!!$  !----------------------------------------------------------------------*
!!$  logical bborbit
!!$  integer ktrack,itrack,ipos,get_option
!!$  double precision track(6,*),parvec(*),pi,sx,sy,xm,ym,sx2,sy2,xs,  &
!!$       ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,cbx,cby,      &
!!$       get_variable,node_value,ten3m,explim
!!$  parameter(ten3m=1d-3,          &
!!$       explim=150d0)
!!$  !     if x > explim, exp(-x) is outside machine limits.
!!$
!!$  !---- initialize.
!!$  bborbit = get_option('bborbit ') .ne. 0
!!$  pi=get_variable('pi ')
!!$  sx = node_value('sigx ')
!!$  sy = node_value('sigy ')
!!$  xm = node_value('xma ')
!!$  ym = node_value('yma ')
!!$  fk = 2d0 * parvec(5) * parvec(6) / parvec(7)
!!$  if (fk .eq. 0d0)  return
!!$  ipos = 0
!!$  if (.not. bborbit)  then
!!$     !--- find position of closed orbit bb_kick
!!$     do ipos = 1, bbd_cnt
!!$        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
!!$     enddo
!!$     ipos = 0
!!$1    continue
!!$  endif
!!$  sx2 = sx*sx
!!$  sy2 = sy*sy
!!$  !---- limit formulae for sigma(x) = sigma(y).
!!$  if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
!!$     do itrack = 1, ktrack
!!$        xs = track(1,itrack) - xm
!!$        ys = track(3,itrack) - ym
!!$        rho2 = xs * xs + ys * ys
!!$        tk = rho2 / (2d0 * sx2)
!!$        if (tk .gt. explim) then
!!$           phix = xs * fk / rho2
!!$           phiy = ys * fk / rho2
!!$        else if (rho2 .ne. 0d0) then
!!$           phix = xs * fk / rho2 * (1d0 - exp(-tk) )
!!$           phiy = ys * fk / rho2 * (1d0 - exp(-tk) )
!!$        else
!!$           phix = 0d0
!!$           phiy = 0d0
!!$        endif
!!$        if (ipos .ne. 0)  then
!!$           !--- subtract closed orbit kick
!!$           phix = phix - bb_kick(1,ipos)
!!$           phiy = phiy - bb_kick(2,ipos)
!!$        endif
!!$        track(2,itrack) = track(2,itrack) + phix
!!$        track(4,itrack) = track(4,itrack) + phiy
!!$     enddo
!!$
!!$     !---- case sigma(x) > sigma(y).
!!$  else if (sx2 .gt. sy2) then
!!$     r2 = 2d0 * (sx2 - sy2)
!!$     r  = sqrt(r2)
!!$     rk = fk * sqrt(pi) / r
!!$     do itrack = 1, ktrack
!!$        xs = track(1,itrack) - xm
!!$        ys = track(3,itrack) - ym
!!$        xr = abs(xs) / r
!!$        yr = abs(ys) / r
!!$        call ccperrf(xr, yr, crx, cry)
!!$        tk = (xs * xs / sx2 + ys * ys / sy2) / 2d0
!!$        if (tk .gt. explim) then
!!$           phix = rk * cry
!!$           phiy = rk * crx
!!$        else
!!$           xb = (sy / sx) * xr
!!$           yb = (sx / sy) * yr
!!$           call ccperrf(xb, yb, cbx, cby)
!!$           phix = rk * (cry - exp(-tk) * cby)
!!$           phiy = rk * (crx - exp(-tk) * cbx)
!!$        endif
!!$        track(2,itrack) = track(2,itrack) + phix * sign(1d0,xs)
!!$        track(4,itrack) = track(4,itrack) + phiy * sign(1d0,ys)
!!$        if (ipos .ne. 0)  then
!!$           !--- subtract closed orbit kick
!!$           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
!!$           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
!!$        endif
!!$     enddo
!!$
!!$     !---- case sigma(x) < sigma(y).
!!$  else
!!$     r2 = 2d0 * (sy2 - sx2)
!!$     r  = sqrt(r2)
!!$     rk = fk * sqrt(pi) / r
!!$     do itrack = 1, ktrack
!!$        xs = track(1,itrack) - xm
!!$        ys = track(3,itrack) - ym
!!$        xr = abs(xs) / r
!!$        yr = abs(ys) / r
!!$        call ccperrf(yr, xr, cry, crx)
!!$        tk = (xs * xs / sx2 + ys * ys / sy2) / 2d0
!!$        if (tk .gt. explim) then
!!$           phix = rk * cry
!!$           phiy = rk * crx
!!$        else
!!$           xb  = (sy / sx) * xr
!!$           yb  = (sx / sy) * yr
!!$           call ccperrf(yb, xb, cby, cbx)
!!$           phix = rk * (cry - exp(-tk) * cby)
!!$           phiy = rk * (crx - exp(-tk) * cbx)
!!$        endif
!!$        track(2,itrack) = track(2,itrack) + phix * sign(1d0,xs)
!!$        track(4,itrack) = track(4,itrack) + phiy * sign(1d0,ys)
!!$        if (ipos .ne. 0)  then
!!$           !--- subtract closed orbit kick
!!$           track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
!!$           track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
!!$        endif
!!$     enddo
!!$  endif
!!$end subroutine ttbb_old
subroutine trkill(n, turn, sum, jmax, part_id,                    &
     last_turn, last_pos, last_orbit, z, aptype)

  use name_lenfi
  use spch_bbfi
  implicit none

  !hbu--- kill particle:  print, modify part_id list
  logical recloss, exit_loss_turn
  integer i,j,n,turn,part_id(*),jmax,last_turn(*),get_option
  double precision sum, z(6,*), last_pos(*), last_orbit(6,*),       &
       torb(6) !, theta, node_value, st, ct, tmp
  character(name_len) aptype
  !hbu
  character(name_len) el_name

  recloss = get_option('recloss ') .ne. 0
  exit_loss_turn = get_option('exit_loss_turn ') .ne. 0

  !!--- As elements might have a tilt we have to transform back
  !!--- into the original coordinate system!
  !      theta = node_value('tilt ')
  !      if (theta .ne. 0d0)  then
  !          st = sin(theta)
  !          ct = cos(theta)
  !!--- rotate trajectory (at exit)
  !            tmp = z(1,n)
  !            z(1,n) = ct * tmp - st * z(3,n)
  !            z(3,n) = ct * z(3,n) + st * tmp
  !            tmp = z(2,n)
  !            z(2,n) = ct * tmp - st * z(4,n)
  !            z(4,n) = ct * z(4,n) + st * tmp
  !      endif

  last_turn(part_id(n)) = turn
  last_pos(part_id(n)) = sum
  do j = 1, 6
     torb(j) = z(j,n)
     last_orbit(j,part_id(n)) = z(j,n)
  enddo

  !hbu
  call element_name(el_name,len(el_name))
  !hbu
  write(6,'(''particle #'',i6,'' lost turn '',i6,''  at pos. s ='', f10.2,'' element='',a,'' aperture ='',a)') &
       part_id(n),turn,sum,el_name,aptype
  print *,"   X=",z(1,n),"  Y=",z(3,n),"  T=",z(5,n)

  if(exit_loss_turn) then
     lost_in_turn = .true.
     is_lost = .true.
  endif
 
  if(recloss) then
     call tt_ploss(part_id(n),turn,sum,torb,el_name)
  endif

  do i = n+1, jmax
     part_id(i-1) = part_id(i)
     do j = 1, 6
        z(j,i-1) = z(j,i)
     enddo
  enddo
  jmax = jmax - 1

end subroutine trkill


subroutine tt_ploss(npart,turn,spos,orbit,el_name)

  use name_lenfi
  implicit none

  !hbu added spos
  !----------------------------------------------------------------------*
  !--- purpose: enter lost particle coordinates in table                 *
  !    input:                                                            *
  !    npart  (int)           particle number                            *
  !    turn   (int)           turn number                                *
  !    spos    (double)       s-coordinate when loss happens             *
  !    orbit  (double array)  particle orbit                             *
  !    orbit0 (double array)  reference orbit                            *
  !----------------------------------------------------------------------*
  integer npart,turn,j
  double precision orbit(6),tmp,tt,tn
  double precision energy,get_value
  character(120) table
  character(name_len) el_name
  !hbu
  double precision spos
  !hbu
  character(4) vec_names(7)
  !hbu
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt', 's' /
  data table / 'trackloss' /

  tn = npart
  tt = turn

  energy = get_value('probe ','energy ')


  ! the number of the current particle
  call double_to_table_curr(table, 'number ', tn)
  ! the number of the current turn
  call double_to_table_curr(table, 'turn ', tt)
  !hbu spos

  call double_to_table_curr(table,vec_names(7),spos)

  do j = 1, 6
     tmp = orbit(j)
     call double_to_table_curr(table, vec_names(j), tmp)
  enddo
  call double_to_table_curr(table, 'e ', energy)
  call string_to_table_curr(table, 'element ', el_name)

  call augment_count(table)
end subroutine tt_ploss


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
  logical first
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
  save first
  data first / .true. /

  !hbu
  length = len(comment)
  segment = segment + 1
  !hbu
  write(comment, '(''#segment'',4i8,1X,A)')                         &
       segment,tot_segm,npart,ielem,el_name
  if(first) call comment_to_table_curr(table, comment, length)
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

subroutine trcoll(flag, apx, apy, apr, turn, sum, part_id, last_turn,  &
     last_pos, last_orbit, z, ntrk,al_errors,offx,offy)

  use twiss0fi
  use name_lenfi
  use Inf_NaN_Detection
  implicit none

  ! 2013-May-22  11:46:15  ghislain: Merged the trcoll and trcoll1 subroutines,
  ! hence added the apr argument and the code for the type 4 case of
  ! racetrack aperture

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   test for collimator aperture limits.                               *
  ! input:                                                               *
  !   flag      (integer)   aperture type flag:                          *
  !                         1: elliptic, 2: rectangular, 3: marguerite   *
  !                         4: racetrack                                 *
  !   apx       (double)    x aperture or half axis                      *
  !   apy       (double)    y aperture or half axis                      *
  !   apr       (double)    radius of curvature for racetrack corners    *
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
  integer flag,turn,part_id(*),last_turn(*),ntrk,i,n,nn,  &
       get_option, optiondebug
  double precision apx,apy,apr,sum,last_pos(*),last_orbit(6,*),z(6,*), &
       al_errors(align_max),offx,offy
  character(name_len) aptype

  optiondebug = get_option('debug ')


  if (optiondebug .ne. 0) then
     print *, "trcoll called with parameters: flag, apx, apy, apr : ", flag, apx, apy, apr
  endif


  n = 1
10 continue
  do i = n, ntrk

     if(ISNAN(z(1,i)).or.ISNAN(z(3,i))) goto 99

     !*** case of ellipse
     if (flag .eq. 1) then
        if ( ((z(1,i)-al_errors(11)-offx)/apx)**2 + &
             ((z(3,i)-al_errors(12)-offy)/apy)**2 .gt. 1d0 ) then
           if (optiondebug .ne. 0) then
              print *, "trcoll ellipse: x, al_err_x, offx, apx : ", z(1,i), al_errors(11), offx, apx
              print *, "                y, al_err_y, offy, apy : ", z(3,i), al_errors(12), offy, apy
           endif
           go to 99
        endif
        go to 98

     !*** case of rectangle
     else if(flag .eq. 2) then
        if ( (abs(z(1,i)-al_errors(11)-offx) .gt. apx .or. &
             abs(z(3,i)-al_errors(12)-offy) .gt. apy) ) then
           if (optiondebug .ne. 0) then
              print *, "trcoll rectangle: x, al_err_x, offx, apx : ", z(1,i), al_errors(11), offx, apx
              print *, "                  y, al_err_y, offy, apy : ", z(3,i), al_errors(12), offy, apy
           endif
           go to 99
        endif
        go to 98

     !***  case of marguerite: two ellipses
     else if(flag .eq. 3) then
        if ( ((z(1,i)-al_errors(11)-offx)/apx)**2 + &
             ((z(3,i)-al_errors(12)-offy)/apy)**2 .gt. 1d0 .and. &
             ((z(1,i)-al_errors(11)-offx)/apy)**2 + &
             ((z(3,i)-al_errors(12)-offy)/apx)**2 .gt. 1d0 ) then
           if (optiondebug .ne. 0) then
              print *, "trcoll marguerite: x, al_err_x, offx, apx : ", z(1,i), al_errors(11), offx, apx
              print *, "                   y, al_err_y, offy, apy : ", z(3,i), al_errors(12), offy, apy
           endif
           go to 99
        endif
        go to 98

     !*** case of racetrack: test outer rectangle (apx+apr,apy+apr) then test radius for rounded corners.
     else if (flag .eq. 4) then
        if ( abs(z(1,i)-al_errors(11)-offx) .gt. (apr+apx) .or. &
             abs(z(3,i)-al_errors(12)-offy) .gt. (apy+apr) .or. &
             ( abs(z(1,i)-al_errors(11)-offx) .gt. apx .and. &
             abs(z(3,i)-al_errors(12)-offy) .gt. apy .and. &
             (abs(z(1,i)-al_errors(11)-offx)-apx)**2 + &
             (abs(z(3,i)-al_errors(12)-offy)-apy)**2 .gt. apr**2 ) ) then
           if (optiondebug .ne. 0) then
              print *, "trcoll racetrack: flag is", flag, "and should be 4"
              print *, "trcoll racetrack: x, al_err_x, offx, apx      : ", z(1,i), al_errors(11), offx, apx
              print *, "                  y, al_err_y, offy, apy, apr : ", z(3,i), al_errors(12), offy, apy, apr
           endif
           go to 99
        endif
        go to 98

     ! add error case when option is not in range 1..4 tested above; we should never get there.
     else
        call aawarn('trcoll:','called with unknown flag option. exit from trcoll')
        return
     endif

     ! lose particle if it is outside aperture
99   n = i
     nn=name_len
     call node_string('apertype ',aptype,nn)
     call trkill(n, turn, sum, ntrk, part_id,                        &
          last_turn, last_pos, last_orbit, z, aptype)
     if(ntrk.eq.0) then
        call aawarn('trcoll: ',&
             'Particle Number equals zero: exit from trcoll')
        return
     endif
     ! particle numbering has been reset by trkill, restart the loop with new paramaters.
     goto 10

98   continue
  enddo
end subroutine trcoll

subroutine ttrfloss(turn, sum, part_id, last_turn, last_pos, last_orbit, z, ntrk)

  use twiss0fi
  use name_lenfi
  use trackfi
  use Inf_NaN_Detection
  implicit none

  ! FRS: August 2013

  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Find particles outside the bucket.                                 *
  ! input:                                                               *
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
  integer turn,part_id(*),last_turn(*),ntrk,i,n,nn,  &
       get_option, optiondebug
  double precision sum,last_pos(*),last_orbit(6,*),z(6,*), &
       al_errors(align_max),offx,offy
  character(name_len) :: non_app="RF-Bucket"

  optiondebug = get_option('debug ')

  n = 1
10 continue
  do i = n, ntrk

     !---- Is particle outside the bucket?
     if(ISNAN(z(5,i)).or.ISNAN(z(6,i))) goto 99
     if(abs(z(5,i)).gt.t_max.or.abs(z(6,i)).gt.pt_max) goto 99
     ! break if particle is inside aperture and continue the loop
     go to 98

     ! lose particle if it is outside aperture
99   n = i
     nn=name_len
     call trkill(n, turn, sum, ntrk, part_id,                        &
          last_turn, last_pos, last_orbit, z, non_app)
     if(ntrk.eq.0) then
        call aawarn('trcoll: ',&
             'Particle Number equals zero: exit from trcoll')
        return
     endif
     goto 10

98   continue
  enddo
end subroutine ttrfloss


! subroutine trcoll1(flag, apx, apy, turn, sum, part_id, last_turn,  &
!      last_pos, last_orbit, z, ntrk,al_errors, apr,offx,offy)

! 2013-May-22  11:50:43  ghislain: commented out after merging the racetrack case inside the trcoll subroutine
! will be deleted after tests are done successfully

!   use twiss0fi
!   use name_lenfi
!   implicit none

!   !----------------------------------------------------------------------*
!   ! Similar with trcoll, for racetrack type aperture
!   !-------------Racetrack type , added by Yipeng SUN 21-10-2008---
!   !----------------------------------------------------------------------*
!   integer flag,turn,part_id(*),last_turn(*),ntrk,i,n,nn
!   double precision apx,apy,sum,last_pos(*),last_orbit(6,*),z(6,*),  &
!        one,al_errors(align_max),apr,offx,offy
!   parameter(one=1d0)
!   character(name_len) aptype
!
!   n = 1
! 10 continue
!   do i = n, ntrk
!
!      !---- Is particle outside aperture?
!      if (flag .eq. 4 .and. &
!          (abs(z(1,i)-al_errors(11)-offx) .gt. (apr+apx) .or. &
!           abs(z(3,i)-al_errors(12)-offy) .gt. (apy+apr) ) .or. &
!           (((abs(z(1,i)-al_errors(11)-offx)-apx)**2 + &
!             (abs(z(3,i)-al_errors(12)-offy)-apy)**2) .gt. apr**2) .and. &
!            abs(z(1,i)-al_errors(11)-offx) .gt. apx .and. &
!            abs(z(3,i)-al_errors(12)-offy) .gt. apy) then
!         go to 99
!      endif
!      go to 98
! 99   n = i
!      nn=name_len
!      call node_string('apertype ',aptype,nn)
!      call trkill(n, turn, sum, ntrk, part_id,                        &
!           last_turn, last_pos, last_orbit, z, aptype)
!      goto 10
! 98   continue
!   enddo
! end subroutine trcoll1

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
       ft,phit,get_value,get_variable,deltax,coords(6,0:turns,*)
  double precision deltap
  character(120) msg(2)

  deltap = get_variable('track_deltap ')
  !---- Initialise orbit, emittances and eigenvectors etc.
  j = 0
  twopi = get_variable('twopi ')
  ex = get_value('probe ','ex ')
  ey = get_value('probe ','ey ')
  et = get_value('probe ','et ')
  bet0  = get_value('beam ','beta ')
  bet0i = 1d0 / bet0
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
        track(6) = deltae + sqrt((1d0+deltap)**2+bet0i**2-1d0)-bet0i
        track(7) = fx
        track(8) = phix
        track(9) = fy
        track(10) = phiy
        track(11) = ft
        track(12) = phit
        do k = 1,12
           if(abs(track(k)).ne.0d0) then
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
  double precision rt(6,6),reval(6),aival(6),eigen(6,6)

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
  double precision fint,e1,h,hgap,corr,tilt,ek(6),rw(6,6),tw(6,6,6),track(6,*),&
       node_value

  call dzero(ek,6)
  call m66one(rw)
  call dzero(tw, 216)

  e1 = node_value('e1 ')
  h = node_value('h ')
  hgap = node_value('hgap ')
  fint = node_value('fint ')
  tilt = node_value('tilt ')
  corr = (h + h) * hgap * fint
  !          print*,"------------------------------------------ "
  !---- Fringe fields effects computed from the TWISS routine tmfrng
  !     tmfrng returns the matrix elements rw(used) and tw(unused)
  !     No radiation effects as it is a pure thin lens with no lrad
  call tmfrng(.false.,h,0d0,e1,0d0,0d0,corr,rw,tw)
  call tmtilt(.false.,tilt,ek,rw,tw)
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
  double precision track(6,*)
  double precision sk,skl,sks,sksl,cosTh,sinTh,Q,R,Z
  double precision xf,yf,pxf,pyf,sigf,psigf,bvk
  double precision onedp,fpsig,fppsig
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
  !      sk    = sks / 2d0 / (1d0 + deltap)
  sk    = sks / 2d0
  skl   = sksl / 2d0
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
     onedp   = sqrt(1d0+2*psigf+(beta**2)*(psigf**2))
     fpsig   = onedp - 1d0
     fppsig  = (1d0+(beta**2)*psigf)/onedp
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
     sigf = track(5,itrack)*beta - 5d-1*(xf**2 + yf**2)*R

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
!$OMP PARALLEL PRIVATE(itrack)
!$OMP DO
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
!$OMP END DO
!$OMP END PARALLEL
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

subroutine trclor(switch,orbit0)

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
  !   switch  (int)         1: RUN, 2: DYNAP fastune                     *
  !   ORBIT0(6) (double)  Closed Orbit Coordinates @ sequence start      *
  !                       (X, PX, Y, PY, T, PT)                          *
  !----------------------------------------------------------------------*
  double precision orbit0(6),z(6,7),zz(6),z0(6,7),z00(6,7),a(6,7),deltap,ddd(6)
  integer itra, itmax, j, bbd_pos, j_tot
  integer code,switch
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

  logical is_drift, is_thin, is_quad, is_matrix, is_dipole

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
  sum     = 0d0
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
  bet0i = 1d0 / bet0
  beti   = 1d0 / betas
  dtbyds = get_value('probe ','dtbyds ')
  deltas = get_variable('track_deltap ')
  deltap = get_value('probe ','deltap ')
  arad = get_value('probe ','arad ')
  dorad = get_value('probe ','radiate ') .ne. 0d0
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
        if (.not.(is_drift() .or. is_thin() .or. is_quad() .or. is_dipole() .or. is_matrix())) then
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
     call ttmap(switch,code,el,z,pmax,dxt,dyt,sum,turn,part_id,             &
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
     err = 0d0
     do i= 1,6
        a(i,i) = a(i,i) - 1d0
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

subroutine ixy_fitting()
  use name_lenfi
  use bbfi
  use spch_bbfi
  implicit none

  logical emittance_update
  integer i, iii, jjj,get_option
  ! n_lines_in_turn_table=1
  double precision Summ_dpi_square, Summ_z_part_square
  double precision Summ_x, Summ_y
  double precision Ix(N_macro_max), Iy(N_macro_max),                &
       dpi(N_macro_max), z_part(N_macro_max)
  double precision Ix_sorted(N_macro_max), Iy_sorted(N_macro_max)
  double precision Ix_min,Iy_min, Ix_min_last,Iy_min_last
  ! I_x/Ex_rms+I_y/Ey_rms <= I_div_E_sum_max
  ! limit for particles taken for Ix, Iy evaluations
  double precision    Ix_i, Iy_i, dpi_i, z_part_i, Ixy_rel_summ,    &
       N_for_I_dble
  double precision get_value
  double precision ce10
  double precision c_sumx, y_sumx, t_sumx, c_sumy, y_sumy, t_sumy, a_sum
  parameter(ce10=1d10)
  integer i_for_I

  !----------------------------------------------------------------------
  emittance_update = get_option('emittance_update ') .ne. 0

  i_for_I=0     ! I-evaluations
  I_div_E_sum_max = get_value('run ', 'i_div_e_sum_max ')
  DO i=1,N_macro_surv
     Ix_i=Ix_array(i)
     Iy_i=Iy_array(i)
     dpi_i=dpi_array(i)
     z_part_i=z_part_array(i)
! Lethal bug found by Valery!!!
!     Ixy_rel_summ=Ix_i/Ex_rms+Iy_i/Ex_rms
     Ixy_rel_summ=Ix_i/Ex_rms+Iy_i/Ey_rms
     if (Ixy_rel_summ .LE. I_div_E_sum_max) then
        i_for_I=i_for_I+1
        Ix(i_for_I)=Ix_i
        Iy(i_for_I)=Iy_i
        dpi(i_for_I)=dpi_i
        z_part(i_for_I)=z_part_i
     endif
  ENDDO
  if(i_for_I.eq.0) then
          call aawarn('trrun: ',&
          'the RMS emittances cannot be calculated: exit from IXY_FITTING');
     return
  endif
  N_for_I=i_for_I

  N_for_I_dble=dble(N_for_I)

  !     rms of the relative momentum
  Summ_dpi_square=0d0
  sigma_p_loop: DO iii=1, N_for_I
     Summ_dpi_square=Summ_dpi_square+dpi(iii)**2
  ENDDO sigma_p_loop
  Summ_dpi_square=Summ_dpi_square/N_for_I_dble
  if (Summ_dpi_square .GE. 0d0) then
     if(emittance_update) sigma_p=sqrt(Summ_dpi_square)
  else
     call aafail('IXY_FITTING: Fatal: ','Summ_dpi_square<0')
  endif

  !     rms of the bunch length
  Summ_z_part_square=0d0
  sigma_z_loop: DO iii=1, N_for_I
     Summ_z_part_square=Summ_z_part_square+z_part(iii)**2
  ENDDO sigma_z_loop
  Summ_z_part_square=Summ_z_part_square/N_for_I_dble
  if (Summ_z_part_square .GE. 0d0) then
     sigma_z=sqrt(Summ_z_part_square)
  else
     call aafail('IXY_FITTING: Fatal: ','Summ_z_part_square<0')
  endif



  !     SORTING Ix

  Ix_min_last=0d0
  iii_loop: DO iii=1, N_for_I

     Ix_min=ce10
!$OMP PARALLEL PRIVATE(jjj) SHARED(Ix_min)
!$OMP DO REDUCTION(MIN:Ix_min)
     jjj_loop: DO jjj=1,N_for_I
        if (Ix(jjj) < Ix_min .AND. Ix(jjj) > Ix_min_last) then
           Ix_min=Ix(jjj)
        endif
     ENDDO jjj_loop
!$OMP END DO
!$OMP END PARALLEL
     Ix_sorted(iii)=Ix_min
     Ix_min_last=Ix_min
  ENDDO iii_loop

  Iy_min_last=0d0
  iii_loop_y: DO iii=1, N_for_I

     Iy_min=ce10
!$OMP PARALLEL PRIVATE(jjj) SHARED(Iy_min)
!$OMP DO REDUCTION(MIN:Iy_min)
     jjj_loop_y: DO jjj=1,N_for_I
        if (Iy(jjj) < Iy_min .AND. Iy(jjj) > Iy_min_last) then
           Iy_min=Iy(jjj)
        endif
     ENDDO jjj_loop_y
!$OMP END DO
!$OMP END PARALLEL
     Iy_sorted(iii)=Iy_min
     Iy_min_last=Iy_min
  ENDDO iii_loop_y


  !     Summ of step-function for Ex/Ey evaluation
  !     Kahan summation algorithm
  Summ_x=0d0
  Summ_y=0d0
  c_sumx=0d0
  c_sumy=0d0
  alpha = get_value('run ', 'alpha ')
!!!!!$OMP PARALLEL PRIVATE(iii,a_sum,c_sumx,c_sumy,y_sumx,y_sumy,t_sumx,t_sumy)
!!!!!$OMP DO REDUCTION(+:Summ_x,Summ_y)
  Summ_loop: DO iii=1,N_for_I
     a_sum=Log(1d0-(alpha+dble(iii-1))/N_for_I_dble)
     y_sumx=a_sum/Ix_sorted(iii)-c_sumx
     y_sumy=a_sum/Iy_sorted(iii)-c_sumy
     t_sumx=Summ_x+y_sumx
     t_sumy=Summ_y+y_sumy
     c_sumx=(t_sumx-Summ_x)-y_sumx
     c_sumy=(t_sumy-Summ_y)-y_sumy
     Summ_x=t_sumx
     Summ_y=t_sumy
  ENDDO Summ_loop
!!!!!$OMP END DO
!!!!!$OMP END PARALLEL
  Ex_rms=-N_for_I_dble
  Ey_rms=Ex_rms/Summ_y
  Ex_rms=Ex_rms/Summ_x

  return

END subroutine ixy_fitting

subroutine ixy_calcs(betas, orbit, z,                             &
     betax_start, betay_start,                    &
     alfax_start, alfay_start,                    &
     gamax_start, gamay_start,                    &
     dx_start,    dpx_start,                      &
     dy_start,    dpy_start)
  use name_lenfi
  use bbfi
  use spch_bbfi
  implicit none

  integer i
  double precision betas,                                           &
       betax_start, betay_start, alfax_start, alfay_start,          &
       gamax_start, gamay_start,                                    &
       dx_start,    dpx_start,  dy_start,   dpy_start,              &
       orbit(6), z(6,N_macro_surv)
  double precision dpi,xi,pxi,yi,pyi
  !-----------------------------------------------------------------------

  !      DPI=pt_000KK/beta_part_ini; z_part=t_000KK*beta_part_ini;
  !      XI=x_000KK-dx_start*DPI; PXI=px_000KK-dpx_start*DPI;
  !      YI=y_000KK-dy_start*DPI; PYI=py_000KK-dpy_start*DPI;
  !      Ix=(gamax_start*XI*XI+2*alfax_start*XI*PXI+betax_start*PXI*PXI)/2;
  !      Iy=(gamay_start*YI*YI+2*alfay_start*YI*PYI+betay_start*PYI*PYI)/2;

  do i=1,N_macro_surv
     ! Exact formulation might be too computational time costly
     !        DPI=(sqrt((1d0+z(6,i)*betas)**2-gammas**(-2)))/betas-1d0
     DPI=(z(6,i)    - orbit(6)) / betas
     XI =z(1,i)     - orbit(1) - dx_start  * DPI
     PXI=z(2,i)     - orbit(2) - dpx_start * DPI
     YI =z(3,i)     - orbit(3) - dy_start  * DPI
     PYI=z(4,i)     - orbit(4) - dpy_start * DPI

     Ix_array(i)=(gamax_start*XI*XI+2d0*alfax_start*XI*PXI+         &
          betax_start*PXI*PXI)/2d0
     Iy_array(i)=(gamay_start*YI*YI+2d0*alfay_start*YI*PYI+         &
          betay_start*PYI*PYI)/2d0
     dpi_array(i)=DPI
     z_part_array(i)=(z(5,i)-orbit(5))*betas
  enddo

end subroutine ixy_calcs

subroutine table_input(                                           &
     betx_start, bety_start,                      &
     alfx_start, alfy_start,                      &
     gamx_start, gamy_start,                      &
     dx_start,    dpx_start,                      &
     dy_start,    dpy_start)

  use name_lenfi
  use bbfi
  use twtrrfi
  use time_varfi
  use spch_bbfi
  implicit none
  integer i,ii,j,flag,range(2),string_from_table_row,advance_to_pos,          &
       double_from_table_row
  double precision position,                                        &
       betx_start,bety_start,alfx_start,alfy_start,gamx_start,gamy_start,&
       dx_start,dpx_start,dy_start,dpy_start
  double precision cme10
  parameter(cme10=1d-10)
  character*(name_len) name
  character*20 text

  call dzero(myfield,n_time_var*2*(maxmul+1))
  call dzero(phase_tromb,n_time_var*36)
  call dzero(cav_volt,n_time_var)
  call dzero(time_var_m_ind,n_time_var)
  call dzero(time_var_p_ind,n_time_var)
  call dzero(time_var_c_ind,n_time_var)
  call dzero(time_var_m_nt,n_time_var)
  call dzero(time_var_p_nt,n_time_var)
  call dzero(time_var_c_nt,n_time_var)

  call table_range('spch_bb ', '#s/#e ', range)
  print *, 'Range for Table spch_bb : ', range(1), range(2)
  if(range(1).eq.0.and.range(2).eq.0) then
     print*," Info: Table spch_bb is empty "
  endif
  name=" "
  if(range(2).gt.bbd_max) then
     write(text, '(1p,i8)') bbd_max
     call aafail('TRRUN: Fatal: ',                                   &
          'overrun of the number of BB elements in table spch_bb =' //    &
          text)
  endif
  N_spch = range(2)-range(1)
  if(N_spch.lt.1) then
     call aafail('TRRUN: Fatal: ',                                  &
          'Table: spch_bb holds no BB elements')
  endif
  do i=range(1),range(2)
     j = advance_to_pos('spch_bb ', i)
     if(i.eq.range(1)) then
        flag = string_from_table_row('spch_bb ', 'name ', i, name)
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 's ',i, position)
        if(flag.ne.0) goto 98
        if(name(:10).ne."RING$START".and.position.ne.0d0) then
           write(text, '(1p,i8)') i
           call aafail('TRRUN: Fatal: ',                             &
                'Global TWISS not readable from table spch_bb'// text)
        endif
        flag = double_from_table_row('spch_bb ', 'betx ', i, betx_start)
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'bety ', i, bety_start)
        if(flag.ne.0) goto 98
        if(abs(betx_start).lt.cme10.or.abs(bety_start).lt.cme10) then
           write(text, '(1p,i8)') i
           call aafail('TRRUN: Fatal: ',                             &
                'start beta values from TWISS table smaller than '// &
                '1e-10, location: ' // text)
        endif
        flag = double_from_table_row('spch_bb ', 'alfx ', i, alfx_start)
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'alfy ', i, alfy_start)
        if(flag.ne.0) goto 98
        gamx_start = (1d0+alfx_start*alfx_start)/betx_start
        gamy_start = (1d0+alfy_start*alfy_start)/bety_start
        flag = double_from_table_row('spch_bb ', 'dx ', i,  dx_start)
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'dpx ', i, dpx_start)
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'dy ', i,  dy_start)
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'dpy ', i, dpy_start)
        if(flag.ne.0) goto 98
     else
        ii=i-1
        flag = string_from_table_row('spch_bb ', 'name ', i, spch_bb_name(ii))
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'betx ', i, betx_bb(ii))
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'bety ', i, bety_bb(ii))
        if(flag.ne.0) goto 98
        if(abs(betx_bb(ii)).lt.cme10.or.abs(bety_bb(ii)).lt.cme10) then
           write(text, '(1p,i8)') i
           call aafail('TRRUN: Fatal: ',                                &
                'BB beta values from TWISS table smaller '//            &
                'than 1e-10 at location:  '//text)
        endif
        flag = double_from_table_row('spch_bb ', 'alfx ', i, alfx_bb(ii))
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'alfy ', i, alfy_bb(ii))
        if(flag.ne.0) goto 98
        gamx_bb(ii) = (1d0+alfx_bb(ii)*alfx_bb(ii))/betx_bb(ii)
        gamy_bb(ii) = (1d0+alfy_bb(ii)*alfy_bb(ii))/bety_bb(ii)
        flag = double_from_table_row('spch_bb ', 'dx ', i,  dx_bb(ii))
        if(flag.ne.0) goto 98
        flag = double_from_table_row('spch_bb ', 'dy ', i,  dy_bb(ii))
        if(flag.ne.0) goto 98
     endif
  enddo
  goto 99
98 write(text, '(1p,i8)') i
  call aafail('TRRUN: Fatal: ',                                     &
       'Table: spch_bb corrupted at row =' // text)
99 continue

  call table_range('time_var_mul ', '#s/#e ', range)
  print *, 'Range for Table time_var_mul: ', range(1), range(2)
  if(range(1).eq.0.and.range(2).eq.0) then
     print*," Info: Table time_var_mul is empty "
     goto 102
  endif
  name=" "
  if(range(2).gt.n_time_var) then
     write(text, '(1p,i8)') n_time_var
     call aafail('TRRUN: Fatal: ',                                   &
          'overrun of time varying mult arrays =' // text)
  endif
  do i=range(1),range(2)
     j = advance_to_pos('time_var_mul ', i)
     flag = double_from_table_row('time_var_mul ', 'number ',i,          &
          time_var_m_ind(i))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'turn ',i,            &
          time_var_m_nt(i))
     if(flag.ne.0) goto 101
     flag = string_from_table_row('time_var_mul ', 'name ', i,              &
          time_var_m_ch(i))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn0 ',i,            &
          myfield(i,1,0))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks0 ',i,            &
          myfield(i,1,1))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn1 ',i,            &
          myfield(i,1,2))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks1 ',i,            &
          myfield(i,1,3))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn2 ',i,            &
          myfield(i,1,4))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks2 ',i,            &
          myfield(i,1,5))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn3 ',i,            &
          myfield(i,1,6))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks3 ',i,            &
          myfield(i,1,7))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn4 ',i,            &
          myfield(i,1,8))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks4 ',i,            &
          myfield(i,1,9))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn5 ',i,            &
          myfield(i,1,10))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks5 ',i,            &
          myfield(i,2,0))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn6 ',i,            &
          myfield(i,2,1))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks6 ',i,            &
          myfield(i,2,2))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn7 ',i,            &
          myfield(i,2,3))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks7 ',i,            &
          myfield(i,2,4))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn8 ',i,            &
          myfield(i,2,5))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks8 ',i,            &
          myfield(i,2,6))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn9 ',i,            &
          myfield(i,2,7))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks9 ',i,            &
          myfield(i,2,8))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dkn10 ',i,           &
          myfield(i,2,9))
     if(flag.ne.0) goto 101
     flag = double_from_table_row('time_var_mul ', 'dks10 ',i,           &
          myfield(i,2,10))
     if(flag.ne.0) goto 101
  enddo
  goto 102
101 write(text, '(1p,i8)') i
  call aafail('TRRUN: Fatal: ',                                     &
       'Table: time_var_mul corrupted at row =' // text)
102 continue

  call table_range('time_var_pha ', '#s/#e ', range)
  print *, 'Range for Table time_var_pha: ', range(1), range(2)
  if(range(1).eq.0.and.range(2).eq.0) then
     print*," Info: Table time_var_pha is empty "
     goto 104
  endif
  name=" "
  if(range(2).gt.n_time_var) then
     write(text, '(1p,i8)') n_time_var
     call aafail('TRRUN: Fatal: ',                                   &
          'overrun of time varying pha arrays =' // text)
  endif
  do i=range(1),range(2)
     j = advance_to_pos('time_var_pha ', i)
     flag = double_from_table_row('time_var_pha ', 'number ',i,          &
          time_var_p_ind(i))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'turn ',i,            &
          time_var_p_nt(i))
     if(flag.ne.0) goto 103
     flag = string_from_table_row('time_var_pha ', 'name ', i,              &
          time_var_p_ch(i))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm11 ',i,            &
          phase_tromb(i,1))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm12 ',i,            &
          phase_tromb(i,2))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm13 ',i,            &
          phase_tromb(i,3))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm14 ',i,            &
          phase_tromb(i,4))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm15 ',i,            &
          phase_tromb(i,5))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm16 ',i,            &
          phase_tromb(i,6))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm21 ',i,            &
          phase_tromb(i,7))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm22 ',i,            &
          phase_tromb(i,8))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm23 ',i,            &
          phase_tromb(i,9))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm24 ',i,            &
          phase_tromb(i,10))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm25 ',i,            &
          phase_tromb(i,11))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm26 ',i,            &
          phase_tromb(i,12))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm31 ',i,            &
          phase_tromb(i,13))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm32 ',i,            &
          phase_tromb(i,14))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm33 ',i,            &
          phase_tromb(i,15))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm34 ',i,            &
          phase_tromb(i,16))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm35 ',i,            &
          phase_tromb(i,17))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm36 ',i,            &
          phase_tromb(i,18))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm41 ',i,            &
          phase_tromb(i,19))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm42 ',i,            &
          phase_tromb(i,20))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm43 ',i,            &
          phase_tromb(i,21))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm44 ',i,            &
          phase_tromb(i,22))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm45 ',i,            &
          phase_tromb(i,23))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm46 ',i,            &
          phase_tromb(i,24))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm51 ',i,            &
          phase_tromb(i,25))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm52 ',i,            &
          phase_tromb(i,26))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm53 ',i,            &
          phase_tromb(i,27))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm54 ',i,            &
          phase_tromb(i,28))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm55 ',i,            &
          phase_tromb(i,29))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm56 ',i,            &
          phase_tromb(i,30))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm61 ',i,            &
          phase_tromb(i,31))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm62 ',i,            &
          phase_tromb(i,32))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm63 ',i,            &
          phase_tromb(i,33))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm64 ',i,            &
          phase_tromb(i,34))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm65 ',i,            &
          phase_tromb(i,35))
     if(flag.ne.0) goto 103
     flag = double_from_table_row('time_var_pha ', 'rm66 ',i,            &
          phase_tromb(i,36))
     if(flag.ne.0) goto 103
  enddo
  goto 104
103 write(text, '(1p,i8)') i
  call aafail('TRRUN: Fatal: ',                                     &
       'Table: time_var_pha corrupted at row =' // text)
104 continue

  call table_range('time_var_cav ', '#s/#e ', range)
  print *, 'Range for Table time_var_cav: ', range(1), range(2)
  if(range(1).eq.0.and.range(2).eq.0) then
     print*," Info: Table time_var_cav is empty "
     goto 106
  endif
  name=" "
  if(range(2).gt.n_time_var) then
     write(text, '(1p,i8)') n_time_var
     call aafail('TRRUN: Fatal: ',                                   &
          'overrun of time varying cav arrays =' // text)
  endif
  do i=range(1),range(2)
     j = advance_to_pos('time_var_cav ', i)
     flag = double_from_table_row('time_var_cav ', 'number ',i,          &
          time_var_c_ind(i))
     if(flag.ne.0) goto 105
     flag = double_from_table_row('time_var_cav ', 'turn ',i,            &
          time_var_c_nt(i))
     if(flag.ne.0) goto 105
     flag = string_from_table_row('time_var_cav ', 'name ', i,              &
          time_var_c_ch(i))
     if(flag.ne.0) goto 105
     flag = double_from_table_row('time_var_cav ', 'volt ',i,            &
          cav_volt(i))
     if(flag.ne.0) goto 105
  enddo
  goto 106
105 write(text, '(1p,i8)') i
  call aafail('TRRUN: Fatal: ',                                     &
       'Table: time_var_cav corrupted at row =' // text)
106 continue

end subroutine table_input

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
       pi,knll,cnll, &
       dd, u, v, dUu, dUv, dux, duy, dvx, dvy, x, y

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

!FIXME Unused dummy argument 'turn'
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

  double precision beta
  double precision node_value, get_value
  double precision bvk, deltap, elrad
  double precision f_errors(0:maxferr)
  integer n_ferr, jtrk, iord, nord

  !--- AL: RF-multipole
  integer dummyi, nn, ns, node_fd_errors, turn, ktrack
  double precision normal(0:maxmul), skew(0:maxmul)
  double precision track(6,*), field(2,0:maxmul)
  double precision field_cos(2,0:maxmul)
  double precision field_sin(2,0:maxmul)
  double precision pc, krf, rfac
  double precision twopi, clight, ten3m, get_variable
  double precision x, y, z, dpx, dpy, dpt
  double precision freq, volt, lag, harmon
  double precision pnl(0:maxmul), psl(0:maxmul)
  complex*16 ii, Cm2, Sm2, Cm1, Sm1, Cp0, Sp0, Cp1, Sp1

  parameter ( ten3m=1d-3)
  parameter ( ii=(0d0,1d0) )

  twopi=get_variable('twopi ')  
  clight=get_variable('clight ')

  !---- Zero the arrays
  call dzero(normal,maxmul+1)
  call dzero(skew,maxmul+1)
  call dzero(pnl,maxmul+1)
  call dzero(psl,maxmul+1)
  call dzero(f_errors,maxferr+1)
  call dzero(field,2*(maxmul+1))

  !---- Read-in the parameters
  freq = node_value('freq ');
  lag = node_value('lag ');
  harmon = node_value('harmon ');
  bvk = node_value('other_bv ')
  elrad = node_value('lrad ')
  deltap = get_value('probe ', 'deltap ')
  dorad = get_value('probe ','radiate ') .ne. 0d0
  arad = get_value('probe ','arad ')
  gammas = get_value('probe ','gamma ')
  beta = get_value('probe ','beta ')
  pc = get_value('probe ','pc ')
  n_ferr = node_fd_errors(f_errors);
  call get_node_vector('knl ', nn, normal)
  call get_node_vector('ksl ', ns, skew)
  call get_node_vector('pnl ', dummyi, pnl);
  call get_node_vector('psl ', dummyi, psl);

  rfac = 0d0

  !---- Set-up some parameters
  volt = bvk * node_value('volt ');
  krf = twopi*freq*1d6/clight;

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
      field_cos(1,iord) = bvk * (normal(iord) * cos(pnl(iord) * twopi - krf * z) + field(1,iord));
      field_sin(1,iord) = bvk * (normal(iord) * sin(pnl(iord) * twopi - krf * z));
      field_cos(2,iord) = bvk * (skew(iord)   * cos(psl(iord) * twopi - krf * z) + field(2,iord));
      field_sin(2,iord) = bvk * (skew(iord)   * sin(psl(iord) * twopi - krf * z));
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
    dpx = -REAL(Cp0);
    dpy = AIMAG(Cp0);
    dpt = (volt * ten3m * sin(lag * twopi - krf * z) / pc - krf * REAL(Sp1));

    !---- Radiation effects at entrance.
    if (dorad  .and.  elrad .ne. 0d0) then
      rfac = arad * gammas**3 * (dpx**2+dpy**2) / (3d0*elrad)
      track(2,jtrk) = track(2,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(2,jtrk)
      track(4,jtrk) = track(4,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(4,jtrk)
      track(6,jtrk) = track(6,jtrk) - rfac * (1d0 + track(6,jtrk)) ** 2
    endif

    !---- Apply the kick
    track(2,jtrk) = track(2,jtrk) + dpx
    track(4,jtrk) = track(4,jtrk) + dpy
    track(6,jtrk) = track(6,jtrk) + dpt

    !---- Radiation effects at exit.
    if (dorad  .and.  elrad .ne. 0d0) then
      track(2,jtrk) = track(2,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(2,jtrk)
      track(4,jtrk) = track(4,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(4,jtrk)
      track(6,jtrk) = track(6,jtrk) - rfac * (1d0 + track(6,jtrk)) ** 2
    endif
  enddo

end subroutine ttrfmult

subroutine tttquad(track, ktrack)

  use twtrrfi
  use trackfi
  implicit none

  !-------------------------*
  ! Andrea Latina 2012-2013 *
  !-------------------------*
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !    Track particle through a general thick quadrupole.                *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer)   Number of surviving tracks.                  *
  !----------------------------------------------------------------------*

  double precision track(6,*)
  integer ktrack

  double precision node_value
  double precision k1, k1s, length
  double precision kk0, kk, ksqrt ! kk0 for the design momentum, kk for this particle's momentum
  double precision x, px, y, py, z, pt
  double precision x_, px_, y_, py_, z_, pt_
  double precision C,  S,  ksqrt_S,  S_over_ksqrt
  double precision Ch, Sh, ksqrt_Sh, Sh_over_ksqrt
  double precision delta_p1
  double precision bet0sqr
  double precision ct, st
  double precision tmp
  logical skew
  integer jtrk

  double precision sqrt2
  parameter ( sqrt2=1.41421356237310d0 )

  !---- Read-in the parameters
  bet0sqr = bet0*bet0;
  k1 = node_value('k1 ');
  k1s = node_value('k1s ');
  length = node_value('l ');

  if ((k1.eq.0d0).and.(k1s.eq.0d0)) then
     call ttdrf(length,track,ktrack);
     return
  else if ((k1.ne.0d0).and.(k1s.ne.0d0)) then
     call aawarn('trrun: ',&
          'a quadrupole cannot have *both* K1 and K1S different than zero!');
     return
  endif

  if (k1s.ne.0d0) then
     kk0 = k1s;
     skew = .true.
  else
     kk0 = k1;
     skew = .false.
  endif

  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1,ktrack
     !---- The particle position
     x  = track(1,jtrk);
     px = track(2,jtrk);
     y  = track(3,jtrk);
     py = track(4,jtrk);
     z  = track(5,jtrk);
     pt = track(6,jtrk);

!!$    !---- Radiation effects at entrance.
!!$    if (dorad  .and.  elrad .ne. 0d0) then
!!$      rfac = arad * gammas**3 * (dpx**2+dpy**2) / (3d0*elrad)
!!$      track(2,jtrk) = track(2,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(2,jtrk)
!!$      track(4,jtrk) = track(4,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(4,jtrk)
!!$      track(6,jtrk) = track(6,jtrk) - rfac * (1d0 + track(6,jtrk)) ** 2
!!$    endif

     !---- If SKEW rotates by -45 degrees
     if (skew) then
        ct =  sqrt2 / 2d0
        st = -sqrt2 / 2d0
        tmp = x
        x = ct * tmp + st * y
        y = ct * y   - st * tmp
        tmp = px
        px = ct * tmp + st * py
        py = ct * py  - st * tmp
     endif

     !---- Computes 1+delta and kk
     delta_p1 = sqrt(pt*pt+2d0*pt/bet0+1d0);
     kk = kk0 / delta_p1;

     !---- Computes the kick
     if (kk.gt.0d0) then
        ksqrt = sqrt(kk);
        C = cos(ksqrt*length);
        S = sin(ksqrt*length);
        Ch = cosh(ksqrt*length);
        Sh = sinh(ksqrt*length);
        ksqrt_S  =  ksqrt*S;
        ksqrt_Sh = -ksqrt*Sh;
        S_over_ksqrt  =  S/ksqrt;
        Sh_over_ksqrt = Sh/ksqrt;
     else
        ksqrt = sqrt(-kk);
        C = cosh(ksqrt*length);
        S = sinh(ksqrt*length);
        Ch = cos(ksqrt*length);
        Sh = sin(ksqrt*length);
        ksqrt_S  = -ksqrt*S;
        ksqrt_Sh =  ksqrt*Sh;
        S_over_ksqrt  =  S/ksqrt;
        Sh_over_ksqrt = Sh/ksqrt;
     endif

     !---- Equations of motion
     !---- X
     x_  = C * x + S_over_ksqrt * px / delta_p1;
     px_ = -ksqrt_S * delta_p1 * x + C * px;
     !---- Y
     y_  = Ch * y + Sh_over_ksqrt * py / delta_p1;
     py_ = -ksqrt_Sh * delta_p1 * y + Ch * py;
     !---- Z
     z_ = z + pt*length*(1d0-bet0sqr)/bet0sqr - &
          (0.5) * (bet0*pt+1d0)/bet0/(delta_p1*delta_p1) * &
          (0.5 * kk0 * ((x*x)*(length-C*S_over_ksqrt) - (y*y)*(length-Ch*Sh_over_ksqrt)) + &
          0.5 * ((px*px)*(length+C*S_over_ksqrt) + (py*py)*(length+Ch*Sh_over_ksqrt)) / delta_p1 - &
          (x*px*(1d0-C*C)+y*py*(1d0-Ch*Ch)));
     !pt_ = pt; ! unchanged

     x = x_;
     y = y_;
     z = z_;
     px = px_;
     py = py_;
     !pt = pt_; ! unchanged

     !---- If SKEW rotates by +45 degrees
     if (skew) then
        ct = sqrt2 / 2d0
        st = sqrt2 / 2d0
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
     track(5,jtrk) = z
     !track(6,jtrk) = pt ! unchanged


!!$    !---- Radiation effects at exit.
!!$    if (dorad  .and.  elrad .ne. 0d0) then
!!$      track(2,jtrk) = track(2,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(2,jtrk)
!!$      track(4,jtrk) = track(4,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(4,jtrk)
!!$      track(6,jtrk) = track(6,jtrk) - rfac * (1d0 + track(6,jtrk)) ** 2
!!$    endif

  enddo

end subroutine tttquad

subroutine tttdipole(track, ktrack)
  
  use twtrrfi
  use trackfi
  implicit none
  
  !-------------------------*
  ! Andrea Latina 2013-2014 *
  !-------------------------*
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !    Track particle through a general thick dipole.                    *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer)   Number of surviving tracks.                  *
  !----------------------------------------------------------------------*
  
  double precision track(6,*)
  integer ktrack
  
  double precision node_value
  double precision L, angle, rho, h, k0
  double precision x, px, y, py, z, pt
  double precision x_, px_, y_, py_, z_, pt_
  double precision delta_plus_1, delta_plus_1_sqr, sqrt_delta_plus_1
  double precision sqrt_h_sqrt_k0, sqrt_h_div_sqrt_k0, sqrt_k0_div_sqrt_h
  double precision C, S, C_sqr
  double precision bet0sqr
  integer jtrk
  
  !---- Read-in the parameters
  bet0sqr = bet0*bet0;
  L = node_value('l ');
  angle = node_value('angle ');
  rho = abs(L/angle);
  h = angle/L;
  k0 = h;

  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1,ktrack
     !---- The particle position
     x  = track(1,jtrk);
     px = track(2,jtrk);
     y  = track(3,jtrk);
     py = track(4,jtrk);
     z  = track(5,jtrk);
     pt = track(6,jtrk);
  
!!$    !---- Radiation effects at entrance.
!!$    if (dorad  .and.  elrad .ne. 0d0) then
!!$      rfac = arad * gammas**3 * (dpx**2+dpy**2) / (3d0*elrad)
!!$      track(2,jtrk) = track(2,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(2,jtrk)
!!$      track(4,jtrk) = track(4,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(4,jtrk)
!!$      track(6,jtrk) = track(6,jtrk) - rfac * (1d0 + track(6,jtrk)) ** 2
!!$    endif

     delta_plus_1_sqr = pt*pt+2.0*pt/bet0+1;
     delta_plus_1 = sqrt(delta_plus_1_sqr);
     sqrt_delta_plus_1 = sqrt(delta_plus_1);
     sqrt_h_sqrt_k0 = sign(sqrt(h*k0),k0);
     sqrt_h_div_sqrt_k0 = sqrt(h/k0);
     sqrt_k0_div_sqrt_h = sqrt(k0/h);
     C=cos(sqrt_h_sqrt_k0*L/sqrt_delta_plus_1);
     S=sin(sqrt_h_sqrt_k0*L/sqrt_delta_plus_1);
     C_sqr = C*C;
     x_ = px*S/(sqrt_delta_plus_1*sqrt_h_sqrt_k0)+x*C-delta_plus_1*C/k0+C/h+delta_plus_1/k0-1.0/h;
     px_ = -sqrt_delta_plus_1*sqrt_h_sqrt_k0*x*S- &
          sqrt_delta_plus_1*sqrt_k0_div_sqrt_h*S+delta_plus_1*sqrt_delta_plus_1*sqrt_h_div_sqrt_k0*S+px*C;
     y_  = y + py * L / delta_plus_1; 
     py_ = py; 
     z_  = z + pt*L*(1.0-bet0sqr)/bet0sqr + &
          (-(0.5)*(bet0*pt+1.0)/bet0/(delta_plus_1**3) * &
          (x*x*delta_plus_1*(h*k0*L-sqrt_delta_plus_1*sqrt_h_sqrt_k0*C*S)*0.5 + &
          px*px*(sqrt_delta_plus_1*C*S/sqrt_h_sqrt_k0+L)*0.5 + &
          px*(-delta_plus_1_sqr*C_sqr/k0+delta_plus_1*C_sqr/h+delta_plus_1_sqr/k0-delta_plus_1/h) + &
          x*(-delta_plus_1**(3.0*0.5)*sqrt_k0_div_sqrt_h*C*S+ &
          delta_plus_1**(5.0*0.5)*sqrt_h_div_sqrt_k0*C*S+ &
          delta_plus_1*k0*L-delta_plus_1_sqr*h*L)+ &
          x*px*delta_plus_1*(C_sqr-1.0) + &
          py*py*L + &
          (-delta_plus_1**(3.0*0.5)*sqrt_k0_div_sqrt_h*C*S*0.5/h + &
          delta_plus_1**(5.0*0.5)*C*S/(sqrt_h_sqrt_k0)+ &
          (-delta_plus_1**(7.0*0.5)*sqrt_h_div_sqrt_k0*C*S*0.5/k0 + &
          delta_plus_1*L*(delta_plus_1*(delta_plus_1*h/k0*0.5-1.0)+k0/h*0.5)))));
     !pt_ = pt; ! unchanged
     
     x = x_;
     y = y_;
     z = z_;
     px = px_;
     py = py_;
     !pt = pt_; ! unchanged
     
     !---- Applies the kick
     track(1,jtrk) = x
     track(2,jtrk) = px
     track(3,jtrk) = y
     track(4,jtrk) = py
     track(5,jtrk) = z
     !track(6,jtrk) = pt ! unchanged
     
!!$    !---- Radiation effects at exit.
!!$    if (dorad  .and.  elrad .ne. 0d0) then
!!$      track(2,jtrk) = track(2,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(2,jtrk)
!!$      track(4,jtrk) = track(4,jtrk) - rfac * (1d0 + track(6,jtrk)) * track(4,jtrk)
!!$      track(6,jtrk) = track(6,jtrk) - rfac * (1d0 + track(6,jtrk)) ** 2
!!$    endif
     
  enddo
  
end subroutine tttdipole

subroutine wzset
  !  *********************************************************************
  !
  !  Adapted from SixTrack version 4.5.04 from 12.01.2014
  !  Authors: Eric McIntosh & Riccardo de Maria
  !
  !  This subroutine must be called before subroutine WZSUB can be used to
  !  compute values of the complex error function w(z).
  !
  !  Parameters xcut and ycut specify the opposite corners (xcut,0) and
  !  (0,ycut) of the rectangle inside which interpolation is to be used
  !  by subroutine WZSUB.
  !
  !  Parameter h is the side of the squares of the interpolation grid.
  !
  !  Parameters nx and ny must be set to the nearest integers to xcut/h
  !  and ycut/h respectively (or to larger values).
  !
  !  Calls MYWWERF new version of (CERN library) WWERF (C335)
  !
  !  (G.A.Erskine, 29.09.1995)
  !
  !  *********************************************************************
  use fasterror
  integer i,j,k
  double precision h,wi,wr,x,y
  parameter ( h = 1.d0/63.d0 )
  save
  !-----------------------------------------------------------------------
  hrecip = 1.d0/h
  kstep = nx+2
  k = 0
  do j=0,ny+1
     do i=0,nx+1
        k = k+1
        !hr05       x=i*h
        x=dble(i)*h                                                  !hr05
        !hr05       y=j*h
        y=dble(j)*h                                                  !hr05
        call mywwerf(x,y,wr,wi)
        wtreal(k)=wr
        wtimag(k)=wi
     end do
  end do
end subroutine wzset
subroutine mywwerf(x,y,wr,wi)
  !  *********************************************************************
  !
  !  Adapted from SixTrack version 4.5.04 from 12.01.2014
  !  Authors: Eric McIntosh & Riccardo de Maria
  !
  !  New version of (CERN library) WWERF (C335)
  !
  !  (G.A.Erskine, 29.09.1995)
  !
  !  *********************************************************************
  implicit none
  integer n
  double precision c,c1,c2,c3,c4,hf,p,rr,ri,sr0,sr,si,tr,ti,vi,vr,  &
       &wi,wr,x,xa,xl,y,ya,zhi,zhr,z1,z10
  parameter (z1=1,hf=z1/2d0,z10=10d0)
  parameter (c1=74d0/z10,c2=83d0/z10,c3=z10/32d0,c4=16d0/z10)
  !     parameter (c=1.12837916709551257d0,p=(2d0*c4)**33)
  parameter (c=1.12837916709551257d0,p=46768052394588893.3825d0)
  dimension rr(37),ri(37)
  save
  !-----------------------------------------------------------------------
  xa=abs(x)
  ya=abs(y)
  if(ya.lt.c1.and.xa.lt.c2) then
     !        zh=dcmplx(ya+c4,xa)
     zhr=ya+c4
     zhi=xa
     rr(37)=0d0
     ri(37)=0d0
     do n=36,1,-1
        !          t=zh+n*dconjg(r(n+1))
        !hr05     tr=zhr+n*rr(n+1)
        tr=zhr+dble(n)*rr(n+1)                                         !hr05
        !hr05     ti=zhi-n*ri(n+1)
        ti=zhi-dble(n)*ri(n+1)                                         !hr05
        !          r(n)=hf*t/(dreal(t)**2+dimag(t)**2)
        !hr05     rr(n)=hf*tr/(tr**2+ti**2)
        rr(n)=(hf*tr)/(tr**2+ti**2)                                    !hr05
        !hr05     ri(n)=hf*ti/(tr**2+ti**2)
        ri(n)=(hf*ti)/(tr**2+ti**2)                                    !hr05
     enddo
     xl=p
     sr=0d0
     si=0d0
     do n=33,1,-1
        xl=c3*xl
        !          s=r(n)*(s+xl)
        sr0=rr(n)*(sr+xl)-ri(n)*si
        si=rr(n)*si+ri(n)*(sr+xl)
        sr=sr0
     enddo
     !        v=c*s
     vr=c*sr
     vi=c*si
  else
     zhr=ya
     zhi=xa
     rr(1)=0d0
     ri(1)=0d0
     do n=9,1,-1
        !          t=zh+n*dconjg(r(1))
        !hr05     tr=zhr+n*rr(1)
        tr=zhr+dble(n)*rr(1)                                           !hr05
        !hr05     ti=zhi-n*ri(1)
        ti=zhi-dble(n)*ri(1)                                           !hr05
        !          r(1)=hf*t/(dreal(t)**2+dimag(t)**2)
        !hr05     rr(1)=hf*tr/(tr**2+ti**2)
        rr(1)=(hf*tr)/(tr**2+ti**2)                                    !hr05
        !hr05     ri(1)=hf*ti/(tr**2+ti**2)
        ri(1)=(hf*ti)/(tr**2+ti**2)                                    !hr05
     enddo
     !        v=c*r(1)
     vr=c*rr(1)
     vi=c*ri(1)
  endif
  !hr05 if(ya.eq.0) then
  if(ya.eq.0d0) then                                                 !hr05
     !        v=dcmplx(exp(-xa**2),dimag(v))
     !hr05   vr=exp_rn(-xa**2)
     vr=exp(-1d0*xa**2)                                              !hr05
  endif
  if(y.lt.0d0) then
     !        v=2*exp(-dcmplx(xa,ya)**2)-v
     !hr05   vr=2d0*exp_rn(ya**2-xa**2)*cos_rn(2d0*xa*ya)-vr
     vr=(2d0*exp(ya**2-xa**2))*cos((2d0*xa)*ya)-vr              !hr05
     vi=(-2d0*exp(ya**2-xa**2))*sin((2d0*xa)*ya)-vi             !hr05
     !hr05   if(x.gt.0) vi=-vi
     if(x.gt.0d0) vi=-1d0*vi                                          !hr05
  else
     !hr05   if(x.lt.0) vi=-vi
     if(x.lt.0d0) vi=-1d0*vi                                          !hr05
  endif
  wr=vr
  wi=vi
  return
end subroutine mywwerf
subroutine wzsubv(n,vx,vy,vu,vv)
  !  *********************************************************************
  !
  !  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
  !  where w(z) is the complex error function defined by formula 7.1.3 in
  !  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
  !  Washington, 1966].  The absolute error of the computed value is less
  !  than 1E-8.
  !
  !  *** Note.  Subroutine WZSET must have been called before this sub-
  !  routine can be used.
  !
  !  For (x,y) inside the rectangle with opposite corners (xcut,0) and
  !  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
  !  lation formula is used.  For (x,y) outside this rectangle, a two-
  !  term rational approximation is used.
  !
  !  (G.A.Erskine, 29.09.1997)
  !
  !  Vectorised for up to 64 argument values by E.McIntosh, 30.10.1997.
  !
  !
  !  Third-order divided-difference interpolation over the corners of a
  !  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
  !  (F.B.Hildebrand New York, 1957), but with complex nodes and
  !  function values].
  !
  !  In the interpolation formula the corners of the grid square contain-
  !  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
  !  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
  !  and 3 respectively, and a preceding 't' indicates twice the value.
  !
  !
  !  Two-term rational approximation to w(z) [Footnote to Table 7.9
  !  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
  !  I.A.Stegun, Washington, 1966), but with additional digits in
  !  the constants]:
  !              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
  !  Maximum absolute error:
  !        <1.E-6  for  x>=4.9  or  y>=4.4
  !        <1.E-7  for  x>=6.1  or  y>=5.7
  !        <1.E-8  for  x>=7.8  or  y>=7.5
  !
  !  *********************************************************************
  implicit none
  dimension vx(*),vy(*),vu(*),vv(*)
  integer i,j,k,n,vmu,vnu
  double precision a1,a2,b1,b2,vd12i,vd12r,vd23i,vd23r,             &
       &vd34i,vd34r,vp,vq,vqsq,vr,vsimag,vsreal,vt,vtdd13i,vtdd13r,       &
       &vtdd24i,vtdd24r,vtdddi,vtdddr,vti,vtr,vu,vusum,vusum3,vv,         &
       &vvsum,vvsum3,vw1i,vw1r,vw2i,vw2r,vw3i,vw3r,vw4i,vw4r,vx,          &
       &vxh,vxhrel,vy,vyh,vyhrel
  integer npart
  parameter(npart = 64)
  integer idim,kstep,nx,ny
  double precision h,half,hrecip,one,wtimag,wtreal,xcut,ycut
  parameter ( xcut = 7.77d0, ycut = 7.46d0 )
  parameter ( h = 1.d0/63.d0 )
  parameter ( nx = 490, ny = 470 )
  parameter ( idim = (nx+2)*(ny+2) )
  parameter ( half = 0.5d0, one = 1.d0 )
  common /wzcom1/ hrecip, kstep
  common /wzcom2/ wtreal(idim), wtimag(idim)
  parameter ( a1 = 0.5124242248d0, a2 = 0.0517653588d0 )
  parameter ( b1 = 0.2752551286d0, b2 = 2.7247448714d0 )
  double precision xm,xx,yy
  parameter (xm=1d120)
  !     temporary arrays to facilitate vectorisation
  integer in,out,ins,outs
  dimension ins(npart),outs(npart)
  !-----------------------------------------------------------------------
  integer ilo,ihi
  !$    integer nt,th,s
  !$    integer omp_get_num_threads,omp_get_thread_num
  !-----------------------------------------------------------------------
  ilo=1
  ihi=n
  !$OMP PARALLEL DEFAULT(PRIVATE), SHARED(n,vx,vy,vu,vv,hrecip,kstep,wtreal,wtimag)
  !$    nt=omp_get_num_threads()
  !$    th=omp_get_thread_num()
  !$    s=n/nt
  !$    ilo=1+s*th
  !$    ihi=ilo+s-1
  !$    if (th.eq.nt-1) ihi=n
  !-----------------------------------------------------------------------
  in=0
  out=0
  do i=ilo,ihi
     if (vx(i).ge.xcut.or.vy(i).ge.ycut) then
        out=out+1
        outs(out)=i
        if (out.eq.npart) then
           !     everything outside the rectangle so approximate
           !     write (*,*) 'ALL outside'
           !     write (*,*) 'i=',i
           do j=1,out
              xx=vx(outs(j))
              yy=vy(outs(j))
              if (xx.ge.xm) xx=xm
              if (yy.ge.xm) yy=xm
              vp=xx**2-yy**2
              vq=(2.d0*xx)*yy
              vqsq=vq**2
              !  First term.
              vt=vp-b1
              vr=a1/(vt**2+vqsq)
              vsreal=vr*vt
              vsimag=-vr*vq
              !  Second term
              vt=vp-b2
              vr=a2/(vt**2+vqsq)
              vsreal=vsreal+vr*vt
              vsimag=vsimag-vr*vq
              !  Multiply by i*z.
              vu(outs(j))=-(yy*vsreal+xx*vsimag)
              vv(outs(j))=xx*vsreal-yy*vsimag
           enddo
           out=0
        endif
     else
        in=in+1
        ins(in)=i
        if (in.eq.npart) then
           !     everything inside the square, so interpolate
           !     write (*,*) 'ALL inside'
           do j=1,in
              vxh = hrecip*vx(ins(j))
              vyh = hrecip*vy(ins(j))
              vmu = int(vxh)
              vnu = int(vyh)
              !  Compute divided differences.
              k = 2 + vmu + vnu*kstep
              vw4r = wtreal(k)
              vw4i = wtimag(k)
              k = k - 1
              vw3r = wtreal(k)
              vw3i = wtimag(k)
              vd34r = vw4r - vw3r
              vd34i = vw4i - vw3i
              k = k + kstep
              vw2r = wtreal(k)
              vw2i = wtimag(k)
              vd23r = vw2i - vw3i
              vd23i = vw3r - vw2r
              vtr = vd23r - vd34r
              vti = vd23i - vd34i
              vtdd24r = vti - vtr
              !hr05 vtdd24i(j) = - ( vtr(j) + vti(j) )
              vtdd24i = -1d0* ( vtr + vti )                             !hr05
              k = k + 1
              vw1r = wtreal(k)
              vw1i = wtimag(k)
              vd12r = vw1r - vw2r
              vd12i = vw1i - vw2i
              vtr = vd12r - vd23r
              vti = vd12i - vd23i
              vtdd13r = vtr + vti
              vtdd13i = vti - vtr
              vtdddr = vtdd13i - vtdd24i
              vtdddi = vtdd24r - vtdd13r
              !  Evaluate polynomial.
              vxhrel = vxh - dble(vmu)
              vyhrel = vyh - dble(vnu)
              vusum3=half*(vtdd13r+                                       &
                   &(vxhrel*vtdddr-vyhrel*vtdddi))
              vvsum3=half*(vtdd13i+                                       &
                   &(vxhrel*vtdddi+vyhrel*vtdddr))
              vyhrel = vyhrel - one
              vusum=vd12r+(vxhrel*vusum3-vyhrel*vvsum3)
              vvsum=vd12i+(vxhrel*vvsum3+vyhrel*vusum3)
              vxhrel = vxhrel - one
              vu(ins(j))=vw1r+(vxhrel*vusum-vyhrel*vvsum)
              vv(ins(j))=vw1i+(vxhrel*vvsum+vyhrel*vusum)
           enddo
           in=0
        endif
     endif
  enddo
  !     everything outside the rectangle so approximate
  !     write (*,*) 'ALL outside'
  !     write (*,*) 'i=',i
  do j=1,out
     xx=vx(outs(j))
     yy=vy(outs(j))
     if (xx.ge.xm) xx=xm
     if (yy.ge.xm) yy=xm
     vp=xx**2-yy**2
     vq=(2.d0*xx)*yy
     vqsq=vq**2
     !  First term.
     vt=vp-b1
     vr=a1/(vt**2+vqsq)
     vsreal=vr*vt
     vsimag=-vr*vq
     !  Second term
     vt=vp-b2
     vr=a2/(vt**2+vqsq)
     vsreal=vsreal+vr*vt
     vsimag=vsimag-vr*vq
     !  Multiply by i*z.
     vu(outs(j))=-(yy*vsreal+xx*vsimag)
     vv(outs(j))=xx*vsreal-yy*vsimag
  enddo
  !     everything inside the square, so interpolate
  !     write (*,*) 'ALL inside'
  do j=1,in
     vxh = hrecip*vx(ins(j))
     vyh = hrecip*vy(ins(j))
     vmu = int(vxh)
     vnu = int(vyh)
     !  Compute divided differences.
     k = 2 + vmu + vnu*kstep
     vw4r = wtreal(k)
     vw4i = wtimag(k)
     k = k - 1
     vw3r = wtreal(k)
     vw3i = wtimag(k)
     vd34r = vw4r - vw3r
     vd34i = vw4i - vw3i
     k = k + kstep
     vw2r = wtreal(k)
     vw2i = wtimag(k)
     vd23r = vw2i - vw3i
     vd23i = vw3r - vw2r
     vtr = vd23r - vd34r
     vti = vd23i - vd34i
     vtdd24r = vti - vtr
     !hr05 vtdd24i(j) = - ( vtr(j) + vti(j) )
     vtdd24i = -1d0* ( vtr + vti )                             !hr05
     k = k + 1
     vw1r = wtreal(k)
     vw1i = wtimag(k)
     vd12r = vw1r - vw2r
     vd12i = vw1i - vw2i
     vtr = vd12r - vd23r
     vti = vd12i - vd23i
     vtdd13r = vtr + vti
     vtdd13i = vti - vtr
     vtdddr = vtdd13i - vtdd24i
     vtdddi = vtdd24r - vtdd13r
     !  Evaluate polynomial.
     vxhrel = vxh - dble(vmu)
     vyhrel = vyh - dble(vnu)
     vusum3=half*(vtdd13r+                                       &
          &(vxhrel*vtdddr-vyhrel*vtdddi))
     vvsum3=half*(vtdd13i+                                       &
          &(vxhrel*vtdddi+vyhrel*vtdddr))
     vyhrel = vyhrel - one
     vusum=vd12r+(vxhrel*vusum3-vyhrel*vvsum3)
     vvsum=vd12i+(vxhrel*vvsum3+vyhrel*vusum3)
     vxhrel = vxhrel - one
     vu(ins(j))=vw1r+(vxhrel*vusum-vyhrel*vvsum)
     vv(ins(j))=vw1i+(vxhrel*vvsum+vyhrel*vusum)
  enddo
  !$OMP END PARALLEL
  return
end subroutine wzsubv
subroutine wzsub(x,y,u,v)
  !  *********************************************************************
  !
  !  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
  !  where w(z) is the complex error function defined by formula 7.1.3 in
  !  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
  !  Washington, 1966].  The absolute error of the computed value is less
  !  than 1E-8.
  !
  !  *** Note.  Subroutine WZSET must have been called before this sub-
  !  routine can be used.
  !
  !  For (x,y) inside the rectangle with opposite corners (xcut,0) and
  !  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
  !  lation formula is used.  For (x,y) outside this rectangle, a two-
  !  term rational approximation is used.
  !
  !  (G.A.Erskine, 29.09.1997)
  !
  !
  !  Third-order divided-difference interpolation over the corners of a
  !  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
  !  (F.B.Hildebrand New York, 1957), but with complex nodes and
  !  function values].
  !
  !  In the interpolation formula the corners of the grid square contain-
  !  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
  !  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
  !  and 3 respectively, and a preceding 't' indicates twice the value.
  !
  !  *********************************************************************
  use fasterror
  integer k,mu,nu
  double precision a1,a2,b1,b2,d12i,d12r,d23i,d23r,d34i,d34r,p,           &
       &q,qsq,r,simag,sreal,t,tdd13i,tdd13r,tdd24i,tdd24r,tdddi,tdddr,ti, &
       &tr,u,usum,usum3,v,vsum,vsum3,w1i,w1r,w2i,w2r,w3i,w3r,w4i,w4r,x,xh,&
       &xhrel,y,yh,yhrel
  double precision half,one,xcut,ycut
  parameter ( xcut = 7.77d0, ycut = 7.46d0 )
  parameter ( half = 0.5d0, one = 1.d0 )
  parameter ( a1 = 0.5124242248d0, a2 = 0.0517653588d0 )
  parameter ( b1 = 0.2752551286d0, b2 = 2.7247448714d0 )
  save
  !-----------------------------------------------------------------------
  if ( x.ge.xcut .or. y.ge.ycut ) goto 1000
  xh = hrecip*x
  yh = hrecip*y
  mu = int(xh)
  nu = int(yh)
  !  Compute divided differences.
  k = 2 + mu + nu*kstep
  w4r = wtreal(k)
  w4i = wtimag(k)
  k = k - 1
  w3r = wtreal(k)
  w3i = wtimag(k)
  d34r = w4r - w3r
  d34i = w4i - w3i
  k = k + kstep
  w2r = wtreal(k)
  w2i = wtimag(k)
  d23r = w2i - w3i
  d23i = w3r - w2r
  tr = d23r - d34r
  ti = d23i - d34i
  tdd24r = ti - tr
  !hr05 tdd24i = - ( tr + ti )
  tdd24i = -1d0* ( tr + ti )                                         !hr05
  k = k + 1
  w1r = wtreal(k)
  w1i = wtimag(k)
  d12r = w1r - w2r
  d12i = w1i - w2i
  tr = d12r - d23r
  ti = d12i - d23i
  tdd13r = tr + ti
  tdd13i = ti - tr
  tdddr = tdd13i - tdd24i
  tdddi = tdd24r - tdd13r
  !  Evaluate polynomial.
  xhrel = xh - dble(mu)
  yhrel = yh - dble(nu)
  usum3 = half*( tdd13r + ( xhrel*tdddr - yhrel*tdddi ) )
  vsum3 = half*( tdd13i + ( xhrel*tdddi + yhrel*tdddr ) )
  yhrel = yhrel - one
  usum = d12r + ( xhrel*usum3 - yhrel*vsum3 )
  vsum = d12i + ( xhrel*vsum3 + yhrel*usum3 )
  xhrel = xhrel - one
  u = w1r + ( xhrel*usum - yhrel*vsum )
  v = w1i + ( xhrel*vsum + yhrel*usum )
  return
  !
  !  Two-term rational approximation to w(z) [Footnote to Table 7.9
  !  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
  !  I.A.Stegun, Washington, 1966), but with additional digits in
  !  the constants]:
  !              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
  !  Maximum absolute error:
  !        <1.E-6  for  x>=4.9  or  y>=4.4
  !        <1.E-7  for  x>=6.1  or  y>=5.7
  !        <1.E-8  for  x>=7.8  or  y>=7.5
  !
1000 p=x**2-y**2
  !hr05 q=2.d0*x*y
  q=(2.d0*x)*y                                                       !hr05
  qsq=q**2
  !  First term.
  t=p-b1
  r=a1/(t**2+qsq)
  sreal=r*t
  !hr05 simag=-r*q
  simag=(-1d0*r)*q                                                   !hr05
  !  Second term
  t=p-b2
  r=a2/(t**2+qsq)
  sreal=sreal+r*t
  simag=simag-r*q
  !  Multiply by i*z.
  !hr05 u=-(y*sreal+x*simag)
  u=-1d0*(y*sreal+x*simag)                                           !hr05
  v=x*sreal-y*simag
  return
  !
end subroutine wzsub
