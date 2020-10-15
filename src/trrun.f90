subroutine trrun(switch, turns, orbit0, rt, part_id, last_turn, last_pos, &
                 z, dxt, dyt, last_orbit, eigen, coords, e_flag, code_buf, & 
                 l_buf)
  use twtrrfi
  use bbfi
  use time_varfi
  use spch_bbfi
  use twiss0fi
  use name_lenfi
  use trackfi
  use fasterror
  use matrices, only : EYE
  use math_constfi, only : zero, one, two
  use code_constfi
  use SpaceCharge
  use track_enums
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !          Interface RUN and DYNAP command to tracking routine         *
  !                                                                      *
  !-- Input:                                                             *
  !   switch  (int)         1: RUN, 2: DYNAP fastune            f         *
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
  integer, intent(IN)  :: switch, turns
  integer, intent(OUT) :: e_flag

  double precision, intent(IN) :: orbit0(6)
  double precision, intent(IN OUT) :: rt(6,6)
  double precision, intent(OUT) :: eigen(6,6), coords(6,0:turns,*)

  integer :: part_id(*), last_turn(*), code_buf(*)
  double precision :: last_pos(*), z(6,*), dxt(*), dyt(*)
  double precision :: last_orbit(6,*),  l_buf(*)
  double precision :: theta
  double precision, dimension (:), allocatable :: theta_buf   
  logical :: onepass, onetable, last_out, info, aperflag, doupdate, debug
  logical :: run=.false.,dynap=.false., thin_foc, onlyaver
  logical, save :: first=.true.
  logical :: fast_error_func
  integer :: i, j, k, code, ffile
  integer :: n_align, nlm, j_tot, turn, nobs, lobs
  integer :: nint, ndble, nchar, char_l, tot_segm, int_arr(1)
  
  double precision :: orbit(6), el, re(6,6), deltap, sum, spos
  double precision :: al_errors(align_max), zz(6), maxaper(6), obs_orb(6)
    
  character(len=12) :: tol_a='maxaper ', char_a=' '
  character(len=name_len) :: el_name
  character(len=4) :: vec_names(7)
  data vec_names /'x', 'px', 'y', 'py', 't', 'pt', 's'/
  character(len=20) :: text

  integer, external :: restart_sequ, advance_node, get_option, node_al_errors, get_nnodes
  double precision, external :: node_value, get_variable, get_value, node_obs_point

  external :: set_tt_attrib, alloc_tt_attrib, set_tt_multipoles, get_tt_multipoles

  ! 2015-Jul-08  19:16:53  ghislain: make code more readable
  run   = switch .eq. 1
  dynap = switch .eq. 2
  allocate ( theta_buf(get_nnodes()) )
  !--- Initialize
  deltap = get_value('probe ','deltap ')
  betas  = get_value('probe ','beta ')
  gammas = get_value('probe ','gamma ')
  dtbyds = get_value('probe ','dtbyds ')
  arad   = get_value('probe ','arad ')
  radiate  = get_value('probe ','radiate ') .ne. zero

  thin_cf =  get_option('thin_cf ').ne.zero
  bet0    =  get_value('beam ','beta ')

  bet0i   = one / bet0
  beti    = one / betas

  deltas  = get_variable('track_deltap ')

  damp = get_option('damp ') .ne. 0
  quantum = get_option('quantum ') .ne. 0

  debug = get_option('debug ') .ne. 0
  thin_foc = get_option('thin_foc ').eq.1

  onlyaver = get_option('only_average ') .ne. 0
  call init_elements()
  !-------added by Yipeng SUN 01-12-2008--------------
  if (deltap .eq. zero) then
     onepass = get_option('onepass ') .ne. 0
     if (.not.onepass) call trclor(switch, orbit0)
  endif

  !---- Set fast_error_func flag to use faster error function
  !---- including tables. Thanks to late G. Erskine
  fast_error_func = get_option('fast_error_func ') .ne. 0
  if (fast_error_func .and. .not.fasterror_on) then
     call wzset
     fasterror_on = .true.
  endif

  !---- get options for space charge variables
  call SC_Init(first, run, dynap, turns);
  
  if (fsecarb) then
     call fort_warn('TRRUN: ','Second order terms of arbitrary Matrix not allowed for tracking.')
     return
  endif

  aperflag = get_option('aperture ') .ne. 0
  e_flag = 0
  ! flag to avoid double entry of last line
  last_out = .false.
  onepass  = get_option('onepass ') .ne. 0
  onetable = get_option('onetable ') .ne. 0

  info = get_option('info ') * get_option('warn ') .gt. 0 ! why not just respect info and warn ?

  if (onepass) RT = EYE
  call trbegn(rt,eigen)
  if (run) then
     ffile = get_value('run ', 'ffile ')
  else if (dynap) then
     ffile = 1
  endif
  ! for one_table case
  if (first) segment = 0

  tot_segm = turns / ffile + 1
  if (mod(turns, ffile) .ne. 0) tot_segm = tot_segm + 1

  if (first) then
     if (checkpnt_restart) then
        read(90,END=100) jmax
        read(90,END=100) Ex_rms0
        read(90,END=100) Ey_rms0
        do i = 1, jmax
           do j=1,6
              read(90,END=100) z(j,i)
           enddo
        enddo
!frs on 07.06.2016 - fixing
!  longitudinal plane must be frozen too!
        read(90,END=100) sigma_t
        read(90,END=100) mean_t
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
     if (jmax .eq. 0) then
        call fort_warn('TRRUN: ','Particle Number equals zero: exit from trrun')
        return
     endif
     if (jmax .gt. max_part) then
        write(text,'(a,d20.12,a)') "Maximum number of particles (",max_part,") exceeded"
        call fort_fail('TRRUN: ',text)
     endif
!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
     do i = 1, jmax
        last_turn(i) = last_turn_keep(i)
        part_id(i) = part_id_keep(i)
        do j = 1,6
           z(j,i) = z_keep(j,i)
        enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif

  !--- jmax may be reduced by particle loss - keep number in j_tot
  j_tot = jmax
  !--- get vector of six coordinate maxapers (both RUN and DYNAP)
  call comm_para(tol_a, nint, ndble, nchar, int_arr, maxaper, char_a, char_l)
  !hbu--- init info for tables initial s position is 0
  spos=0 ; nlm=0 ; el_name='start           '

  if (first) then
     !--- enter start coordinates in summary table
     if (run) then
        t_max = get_value('probe ','circ ') / get_value('run ', 'track_harmon ') / betas
        pt_max = get_value('run ', 'deltap_max ')
        pt_max = (sqrt( (betas*(pt_max + one))**2 + one/gammas**2) - one) * beti
     else if (dynap) then
        t_max  = 1d20
        pt_max = 1d20
     endif

     do  i = 1, j_tot
        if (abs(z(5,i)) .gt. t_max) then
           write(text, '(1p,d13.5,a1,i6)') t_max,"p",i
           call fort_fail('TRACK_INITIAL: ','Fatal: T-Coordinate larger than' // text)
        endif
        if (abs(z(6,i)).gt.pt_max) then
           write(text, '(1p,d13.5,a1,i6)') pt_max,"p",i
           call fort_fail('TRACK_INITIAL: ','Fatal: PT-Coordinate larger than' // text)
        endif
        call double_to_table_curr('tracksumm ', 'number ', dble(i))
        call double_to_table_curr('tracksumm ', 'turn ', dble(tot_turn))
        do j = 1, 6
           call double_to_table_curr('tracksumm ', vec_names(j), z(j,i) - orbit0(j))
        enddo
        call double_to_table_curr('tracksumm ', vec_names(7), spos)
        call augment_count('tracksumm ')
     enddo
  endif

  !--- enter first turn, and possibly eigen in tables
  if (run) then
     if (onetable)  then
        if (first) then
           call track_pteigen(eigen)
           call tt_putone(jmax, tot_turn, tot_segm, segment, part_id, &
                z, orbit0, spos, nlm, el_name, onlyaver)
        endif
     else
        do i = 1, jmax
           call tt_puttab(part_id(i), 0, 1, z(1,i), orbit0, spos)
        enddo
     endif
  endif

  ORBIT = ORBIT0

  doupdate = get_option('update ') .ne. 0

  !--- loop over turns
  nobs = 0

  call BB_Init(first);
  
  turnloop: do turn = 1, turns

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

     call BB_Update(jmax, orbit0, z);

     nlm = 0
     sum = zero




     do !--- loop over nodes

        bbd_pos = j
        if (turn .eq. 1)  then

           code = node_value('mad8_type ')
!           if (code .eq. code_tkicker)     code = code_kicker
           if (code .eq. code_placeholder) code = code_instrument

           el = node_value('l ')
           theta = node_value('tilt ')
           theta_buf(nlm+1) = theta
           code_buf(nlm+1) = code
           l_buf(nlm+1) = el
           !param(nlm+1, enum_bvk) = 
           !param(nlm+1, enum_lrad) -= 
           !param(nlm+1, enum_bvk)
           !param(nlm+1, enum_bvk)

           if ((code.eq.code_sextupole .or. &
              code.eq.code_octupole .or. &
              code.eq.code_elseparator .or. &
              code.eq.code_crabcavity) .and. el.ne.zero) then
              !if (.not. (is_drift() .or. is_thin() .or. is_quad() .or. is_dipole() .or. is_matrix()) ) then
              call element_name(el_name,len(el_name))
              print *," "
              print *,el_name, "code: ",code," el: ",el,"   THICK ELEMENT FOUND"
              print *," "
              print *,"Track dies nicely"
              print *,"Thick lenses will get nowhere"
              print *,"MAKETHIN will save you"
              print *," "
              print *," "

              sum = node_value('name ')
              call fort_fail('TRRUN: Fatal ','----element with length found : CONVERT STRUCTURE WITH MAKETHIN')
           endif

        else
           el = l_buf(nlm+1)
           code = code_buf(nlm+1)
           theta = theta_buf(nlm+1)
        endif

        if (run) nobs = node_obs_point()

        !--------  Misalignment at beginning of element (from twissfs.f)
        if (code .ne. code_drift)  then
           AL_ERRORS = zero
           n_align = node_al_errors(al_errors)
           if (n_align .ne. 0)  then
              do i = 1, jmax
                 ZZ(:) = Z(:,i)
                 call tmali1(zz, al_errors, betas, gammas, z(1,i), re)
              enddo
           endif
        endif
        
        !-------- Track through element  // suppress dxt 13.12.04
        call ttmap(switch, code, el, z, jmax, dxt, dyt, sum, tot_turn+turn, part_id, &
             last_turn, last_pos, last_orbit, aperflag, maxaper, al_errors, onepass, debug, theta, thin_foc)

        !-------- Space Charge update
        call SC_Update(orbit0, z);
        
        !--------  Misalignment at end of element (from twissfs.f)
        if (code .ne. code_drift .and. n_align.ne.0)  then
           do i = 1, jmax
              ZZ(:) = Z(:,i)
              call tmali2(el, zz, al_errors, betas, gammas, z(1,i), re)
           enddo
        endif

        nlm = nlm + 1
        if (nobs .gt. 0)  then
           OBS_ORB = zero
           call get_node_vector('obs_orbit ', lobs, obs_orb)
           if (lobs .lt. 6) &
                call fort_fail('TRACK: Fatal ', 'obs. point orbit not found')
           if (onetable)  then
              spos = sum
              call element_name(el_name,len(el_name))
              call tt_putone(jmax, tot_turn+turn, tot_segm, segment, part_id, &
                   z, obs_orb,spos,nlm,el_name,onlyaver)
           else
              if (mod(turn, ffile) .eq. 0)  then
                 do i = 1, jmax
                    call tt_puttab(part_id(i), turn, nobs, z(1,i), obs_orb,spos)
                 enddo
              endif
           endif
        endif

        if (advance_node() .eq. 0)  exit

        j=j+1
     end do !--- end of loop over nodes

     if (run) then
        if (mod(turn, ffile) .eq. 0)  then
           if (turn .eq. turns)  last_out = .true.
           if (onetable)  then
              spos=sum
              call element_name(el_name,len(el_name))
              call tt_putone(jmax, tot_turn+turn, tot_segm, segment, part_id, &
                   z, orbit0,spos,nlm,el_name,onlyaver)
           else
              do i = 1, jmax
                 call tt_puttab(part_id(i), turn, 1, z(1,i), orbit0, spos)
              enddo
           endif
        endif
     else if (dynap) then
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

     if (jmax .eq. 0 .or. (switch .gt. 1 .and. jmax .lt. j_tot)) goto 20 ! break loop over turns
     if (dynap .and. info .and. mod(turn,100) .eq. 0) print *, 'turn :', turn

     if (exit_loss_turn .and. lost_in_turn) then
        lost_in_turn = .false.
        write(text, '(i6)') turn
        call fort_warn('TRRUN Frozen SC: ', 'Stop in Loss Turn: '//text)
        goto 20 ! break loop over turns
     endif

  end do turnloop !--- end of loop over turns

20 continue

  if (dynap .and. jmax.lt.j_tot)  then
     e_flag = 1
     return
  endif

  do i = 1, jmax
     last_turn(part_id(i)) = min(turns+tot_turn, turn+tot_turn)
     last_pos(part_id(i)) = sum
     last_orbit(1:6,part_id(i)) = z(1:6,i)
  enddo
  turn = min(turn, turns)

  call BB_Update2(jmax, orbit0, z, part_id, last_turn);
  
  !--- enter last turn in tables if not done already
  if (run .and. .not.last_out)  then
     if (onetable)  then
        spos = sum
        call element_name(el_name,len(el_name))
        call tt_putone(jmax, tot_turn+turn, tot_segm, segment, part_id, &
             z, orbit0,spos,nlm,el_name,onlyaver)
     else
        do i = 1, jmax
           call tt_puttab(part_id(i), turn, 1, z(1,i), orbit0,spos)
        enddo
     endif
  endif

  !--- Write checkpoint_restart data
  call BB_Write(jmax, orbit0, z);
  
  !--- enter last turn in summary table
  do  i = 1, j_tot
     call double_to_table_curr('tracksumm ', 'number ', dble(i))
     call double_to_table_curr('tracksumm ', 'turn ', dble(last_turn(i)))
     do j = 1, 6
        call double_to_table_curr('tracksumm ', vec_names(j), last_orbit(j,i) - orbit0(j))
     enddo
     spos = last_pos(i)
     call double_to_table_curr('tracksumm ',vec_names(7),spos)
     call augment_count('tracksumm ')
  enddo

  return

100 call fort_fail('TRACK: Fatal ', 'checkpoint_restart file corrupted')
end subroutine trrun

subroutine init_elements
  use track_enums
  use trackfi
  use twtrrfi
  use code_constfi
  implicit none
  logical:: aperflag
  integer:: j, code
  integer, external :: restart_sequ, advance_node, get_option
  double precision, external :: node_value
  external :: update_node_aperture
 
  aperflag = get_option('aperture ') .ne. 0
  
  j = restart_sequ()
  do !---- loop over nodes
    code    = node_value('mad8_type ')
  !        if (code .eq. code_tkicker)     code = code_kicker
    if(code .eq. code_multipole) then
     call alloc_tt_attrib(total_enums)
     call set_tt_attrib(enum_other_bv, node_value('other_bv '))
     call set_tt_attrib(enum_lrad, node_value('lrad '))
     call set_tt_attrib(enum_noise, node_value('noise '))
     call set_tt_attrib(enum_angle, node_value('angle '))
     call set_tt_attrib(enum_time_var, node_value('time_var '))
     call set_tt_multipoles(maxmul)
    endif

    if(code.eq.code_hkicker .or. code.eq.code_vkicker .or. &
      code.eq.code_kicker .or.  code.eq.code_tkicker) then
      call alloc_tt_attrib(total_enums)
      call set_tt_attrib(enum_other_bv, node_value('other_bv '))
      call set_tt_attrib(enum_sinkick, node_value('sinkick '))
      call set_tt_attrib(enum_kick, node_value('kick '))
      call set_tt_attrib(enum_chkick, node_value('chkick '))
      call set_tt_attrib(enum_cvkick, node_value('cvkick '))
      call set_tt_attrib(enum_hkick, node_value('hkick '))
      call set_tt_attrib(enum_vkick, node_value('vkick '))
    endif
    if(aperflag .and. code .ne. code_drift) then
       call update_node_aperture()
    endif
    
    if (advance_node() .eq. 0)  exit

  end do !--- end of loop over nodes to set upt things

end subroutine init_elements

subroutine ttmap(switch,code,el,track,ktrack,dxt,dyt,sum,turn,part_id, &
     last_turn,last_pos,last_orbit,aperflag,maxaper,al_errors,onepass, debug, theta, thin_foc)
  use twtrrfi
  use twiss0fi
  use name_lenfi
  use math_constfi, only : zero, one
  use code_constfi
  use aperture_enums
  use trackfi
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
  integer :: switch, code, ktrack
  integer :: turn, part_id(*), last_turn(*)
  double precision :: el, sum
  double precision :: track(6,*), dxt(*), dyt(*), last_pos(*)
  double precision :: last_orbit(6,*), maxaper(6), al_errors(align_max)
  logical :: aperflag, onepass, lost_global

  logical :: fmap, debug
  logical :: run=.false., dynap=.false., thin_foc
  integer :: i, nn, jtrk, apint
  double precision :: ct, tmp, st, theta
  double precision :: ap1, ap2, ap3, ap4, aperture(maxnaper)
  double precision :: offset(2), offx, offy
  double precision :: ek(6), re(6,6), te(6,6,6), craporb(6)
  character(len=name_len) :: aptype

  integer, external :: get_option, node_apertype
  double precision, external :: get_value, node_value

  double precision, parameter :: min_double = 1.e-36
  external :: node_aperture_vector, node_aperture_offset
  ! 2015-Jul-08  19:16:53  ghislain: make code more readable
  run   = switch .eq. 1
  dynap = switch .eq. 2

  fmap=.false.
 
  !---- Drift space; no rotation or aperture check, go straight to tracking and return
  if (code .eq. code_drift) then
     call ttdrf(el,track,ktrack)
     sum = sum + el
     return
  endif

  !---- Rotate trajectory before entry

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

  !---- Test aperture. ALL ELEMENTS BUT DRIFTS and BEAMBEAM
  !     print *, "apint", apint, "ap_notset", ap_notset
  if (aperflag .and. code.ne.code_beambeam) then
     nn=name_len
    
     apint=node_apertype()

     if(apint .eq. ap_notset) then
    ! make global check even if aperture is not defined
    lost_global =.false.
      do jtrk = 1,ktrack
        lost_global =  ISNAN(track(2,jtrk)) .or. ISNAN(track(4,jtrk)) .or. &
          ISNAN(track(5,jtrk)) .or. ISNAN(track(6,jtrk))                                .or. &
          abs(track(1, jtrk)) .gt. maxaper(1) .or.  abs(track(2, jtrk)) .gt. maxaper(2) .or. &
          abs(track(3, jtrk)) .gt. maxaper(3) .or.  abs(track(4, jtrk)) .gt. maxaper(4) .or. &
          abs(track(5, jtrk)) .gt. maxaper(5) .or.  abs(track(6, jtrk)) .gt. maxaper(6)
          if(lost_global) then
          APERTURE(:maxnaper) = zero
          call get_node_vector('aperture ',nn,aperture)

          OFFSET = zero
          call get_node_vector('aper_offset ',nn,offset)
          call trcoll(apint,  aperture, offset, al_errors,  maxaper, &
                turn, sum, part_id, last_turn, last_pos, last_orbit, track, ktrack, debug)
          EXIT ! They are anway checked against all the particles so no need to continue to loop
          endif 
        enddo 

     else
     !APERTURE(:maxnaper) = zero
     !call get_node_vector('aperture ',nn,aperture)
     call node_aperture_vector(aperture)
     call node_aperture_offset(offset)
     !OFFSET = zero
     !call get_node_vector('aper_offset ',nn,offset)

     if (debug) then
        print *, " aperture type       ", apint
        print *, "          parameters ", (aperture(i),i=1,4)
        print *, "          offsets    ", (offset(i),i=1,2)
        print *, " alignment errors    ", (al_errors(i),i=11,12)
        print *, " maximum aperture    ", (offset(i),i=1,2)
        print *, " "
     endif

     call trcoll(apint,  aperture, offset, al_errors,  maxaper, &
          turn, sum, part_id, last_turn, last_pos, last_orbit, track, ktrack, debug)
     endif
  endif ! Test aperture

  ! Switch based on element code for element specific tracking
  select case (code)

    case (code_rbend, code_sbend)
       call tttdipole(track,ktrack, code)

    case (code_matrix)
       EK(:6) = zero
       RE(:6,:6) = zero
       TE(:6,:6,:6) = zero
       CRAPORB(:6) = zero
       call tmarb(.false.,.false.,craporb,fmap,ek,re,te)
       call tttrak(ek,re,track,ktrack)

    case (code_quadrupole)
       call tttquad(track,ktrack)

    case (code_multipole)
      if(thin_cf .and. node_value('lrad ') .gt. zero ) then
        call ttmult_cf_mini(track,ktrack,dxt,dyt,turn,thin_foc)
      else
        call ttmult(track,ktrack,dxt,dyt,turn,thin_foc)
      endif

    case (code_solenoid)
       call trsol(track, ktrack,dxt,dyt)

    case (code_rfcavity)
       call ttrf(track,ktrack)
       if (run) &
            call ttrfloss(turn,sum,part_id,last_turn,last_pos,last_orbit,track,ktrack)

    case (code_elseparator)
       call ttsep(track, ktrack)

    case (code_srotation)
       call ttsrot(track, ktrack)

    case (code_yrotation)
       call ttyrot(track, ktrack)

    case (code_xrotation)
       call ttxrot(track, ktrack)

    case (code_hkicker, code_vkicker, code_kicker, code_tkicker)
       call ttcorr(el, track, ktrack, turn, code)

    !case (code_ecollimator)
    !   call fort_warn('TRRUN: ','found deprecated ECOLLIMATOR element; should be replaced by COLLIMATOR')

   ! case (code_rcollimator)
   !    call fort_warn('TRRUN: ','found deprecated RCOLLIMATOR element; should be replaced by COLLIMATOR')

    case (code_beambeam)
       call ttbb(track, ktrack)

    case (code_twcavity)
       ! call ttlcav(el, track, ktrack)

    case (code_dipedge)
       call ttdpdg(track,ktrack)

    case (code_translation)
       if (onepass) call tttrans(track,ktrack)

    case (code_crabcavity)
       call ttcrabrf(track,ktrack,turn)

    case (code_hacdipole)
       call tthacdip(track,ktrack,turn)

    case (code_vacdipole)
       call ttvacdip(track,ktrack,turn)

    case (code_nllens)
       call ttnllens(track,ktrack)

    case (code_rfmultipole)
       call ttrfmult(track,ktrack,turn)

    case(code_changerefp0)
      call ttchangep0(track,ktrack)

    case (code_hmonitor:code_rcollimator, code_instrument, &
        code_slmonitor:code_imonitor, code_placeholder, code_collimator)
        if(el .gt. 0) call ttdrf(el,track,ktrack)
    case default ! The rest: do nothing

  end select

  ! This is where we should Test Aperture at exit of elements
  ! by calling trcoll again

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
  sum = sum + el

  return
end subroutine ttmap

SUBROUTINE  ttmult_cf_mini(track,ktrack,dxt,dyt,turn, thin_foc)
  use twtrrfi
  use twissbeamfi, only : deltap, beta
  use math_constfi, only : zero, one, two, three
  use time_varfi
  use trackfi
  use time_varfi
  use track_enums

  implicit none 
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Computes thin-lens kick through combined-function magnet.        *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !     Detailed description:                                            *
  !     Phys. Rev. AccelBeams 19.054002 by M. Titze
  !----------------------------------------------------------------------*
  double precision :: track(6,*), dxt(*), dyt(*), ttt
  logical ::  time_var,thin_foc
  integer :: ktrack, turn
  logical :: fsec, ftrk, fmap
  integer :: nord, k, j, nn, ns, bvk, iord, n_ferr, jtrk, nd
  integer, external :: Factorial
  double precision :: dpx, dpy, tilt, kx, ky, elrad, bp1, h0
  double precision :: dipr, dipi, dbr, dbi, dtmp, an, angle, tilt2
  double precision :: gstr, sstr, x, px0, y, py0, orb50, orb60, deltapp
  double precision :: normal(0:maxmul), skew(0:maxmul), f_errors(0:maxferr)
  !double precision :: orbit(6),
  double complex :: kappa, barkappa, sum0, del_p_g, pkick, dxdpg, dydpg, &
                    dxx, dxy, dyy, rp, rm
  double complex :: lambda(0:maxmul)
  double complex :: g(0:maxmul, 0:maxmul)

  double precision, external :: node_value, get_tt_attrib, get_value
  integer, external :: node_fd_errors
  fmap = .true.
  ! Read magnetic field components & fill lambda's according to field
  ! components relative to given plane
  F_ERRORS(0:maxferr) = zero
  n_ferr = node_fd_errors(f_errors)

  bvk = get_tt_attrib(enum_other_bv)
    !---- Multipole length for radiation.
  elrad = get_tt_attrib(enum_lrad)
  an = get_tt_attrib(enum_angle)
  time_var = get_tt_attrib(enum_time_var) .ne. 0  
  
  !dbr = bvk * f_errors(0) !field(1,0)
  !dbi = bvk * f_errors(1) !field(2,0)

  !---- Nominal dipole strength.
  !      dipr = bvk * vals(1,0) / (one + deltas)
  !      dipi = bvk * vals(2,0) / (one + deltas)

  bet0   = get_value('probe ','beta ')
  !---- Multipole components.
  NORMAL(0:maxmul) = zero! ; call get_node_vector('knl ',nn,normal)
  SKEW(0:maxmul) = zero  ! ; call get_node_vector('ksl ',ns,skew)

  call get_tt_multipoles(nn,normal,ns,skew)

  dipr = bvk * normal(0) !vals(1,0)
  dipi = bvk * skew(0)  


     ! cf magnet with quadrupole & sextupole
     gstr = normal(1)/elrad
     sstr = normal(2)/elrad

    do jtrk = 1,ktrack
       x = track(1,jtrk)
       px0 = track(2,jtrk)
       y = track(3,jtrk)
       py0 = track(4,jtrk)
       orb50 = track(5,jtrk)
       orb60 = track(6,jtrk)

       ! get \Delta p/p + 1 (which we denote by the variable
       ! deltapp here) out of orbit(6), for the corresponding
       ! particle
       deltapp = bet0i*sqrt((1d0 + bet0*orb60)**2 - 1 + bet0**2)
       
       ! orbit transformation:
       ! attention: The following formulas constitute only the kick part
       ! of the CF in a drift-kick-drift decomposition.
       track(2, jtrk) = (elrad**2*(-gstr*x - 0.5*x**2*sstr +&
     &0.5*y**2*sstr) + elrad*px0 +&
     &elrad*(dipi*y**3*sstr/6.0d0 - dipr*gstr*x**2 +&
     &0.5*dipr*gstr*y**2 - 0.5*dipr*x**3*sstr + dipr*x*y**2*sstr +&
     &dipr*deltapp - dipr) -&
     &dipr*(-dipi*gstr*y**3/6.0d0 + dipi*y +&
     &dipr*x))/elrad

       track(4, jtrk) = (elrad**2*y*(gstr + x*sstr) + elrad*py0 +&
     &elrad*(0.5*dipi*gstr*y**2 + 0.5*dipi*x*y**2*sstr + dipi*deltapp&
     & - dipi + dipr*gstr*x*y + dipr*x**2*y*sstr -&
     &dipr*y**3*sstr/6.0d0) -&
     &dipi**2*gstr*y**3/6.0d0 - dipi**2*y +&
     &0.5*dipi*dipr*gstr*x*y**2 - dipi*dipr*x +&
     &dipr**2*gstr*y**3/6.0d0)/elrad

       track(5, jtrk) = (bet0*orb50*deltapp - (bet0*orb60 +&
     &1.0)*(dipi*y + dipr*x))/(bet0*deltapp)

    enddo

end SUBROUTINE ttmult_cf_mini

SUBROUTINE  ttmult_cf(track,ktrack,dxt,dyt,turn, thin_foc)
  use twtrrfi
  use twissbeamfi, only : beta
  use math_constfi, only : zero, one, two, three
  use time_varfi
  use trackfi
  use time_varfi
  use track_enums

  implicit none 
  !----------------------------------------------------------------------*
  !     Purpose:                                                         *
  !     Computes thin-lens kick through combined-function magnet.        *
  !     Input:                                                           *
  !     fsec      (logical) if true, return second order terms.          *
  !     ftrk      (logical) if true, track orbit.                        *
  !     Input/output:                                                    *
  !     orbit(6)  (double)  closed orbit.                                *
  !     Output:                                                          *
  !     fmap      (logical) if true, element has a map.                  *
  !     re(6,6)   (double)  transfer matrix.                             *
  !     te(6,6,6) (double)  second-order terms.                          *
  !     Detailed description:                                            *
  !     Phys. Rev. AccelBeams 19.054002 by M. Titze
  !----------------------------------------------------------------------*
  double precision :: track(6,*), dxt(*), dyt(*), ttt
  logical ::  time_var,thin_foc
  integer :: ktrack, turn
  logical :: fsec, ftrk, fmap
  integer :: nord, k, j, nn, ns, bvk, iord, n_ferr, jtrk, nd
  integer, external :: Factorial
  double precision :: dpx, dpy, tilt, kx, ky, elrad, bp1, h0
  double precision :: dtmp, an, angle, tilt2, etahat
  double precision :: normal(0:maxmul), skew(0:maxmul), f_errors(0:maxferr)
  double complex :: kappa, barkappa, sum0, del_p_g, pkick, dxdpg, dydpg, &
                    dxx, dxy, dyy, rp, rm
  double complex :: lambda(0:maxmul)
  double complex :: g(0:maxmul, 0:maxmul)

  double precision, external :: node_value, get_tt_attrib
  integer, external :: node_fd_errors
  fmap = .true.
  
  ! Read magnetic field components & fill lambda's according to field
  ! components relative to given plane
  F_ERRORS(0:maxferr) = zero
  n_ferr = node_fd_errors(f_errors)

  bvk = get_tt_attrib(enum_other_bv)
    !---- Multipole length for radiation.
  elrad = get_tt_attrib(enum_lrad)
  an = get_tt_attrib(enum_angle)
  time_var = get_tt_attrib(enum_time_var) .ne. 0  

  !---- Multipole components.
  NORMAL(0:maxmul) = zero! ; call get_node_vector('knl ',nn,normal)
  SKEW(0:maxmul) = zero  ! ; call get_node_vector('ksl ',ns,skew)
  tilt2 = 0
  call get_tt_multipoles(nn,normal,ns,skew)

  !---- Angle (no bvk in track)
  if (an .ne. 0) f_errors(0) = f_errors(0) + normal(0) - an

  !-----added FrankS, 10-12-2008
  nd = 2 * max(nn, ns, n_ferr/2-1)

  !Below here should not be commented output
  !---- Other components and errors.
  ! that loop should start at one since nominal dipole strength already taken into account above
  !needs to be here though
  nord = 0
  do iord = 1, nd/2
     f_errors(2*iord)   = bvk * (f_errors(2*iord) + normal(iord))
     f_errors(2*iord+1) = bvk * (f_errors(2*iord+1) + skew(iord))
     if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero) nord=iord
  enddo
  !Done with all the setting up... 

  lambda(0:maxmul) = 0

  if (elrad.gt.zero) then
    lambda(0) = (normal(0) + (0, 1)*skew(0))/elrad
     do k = 1, nord
        lambda(k) = (f_errors(2*k) + (0, 1)*f_errors(2*k+1))/elrad/Factorial(k)
     enddo
  else
     lambda = zero
  endif

  kx = real(lambda(0))    ! N.B. B_y |_{\varphi = tilt, r = 0} = kx
  ky = - aimag(lambda(0)) !      B_x |_{\varphi = tilt, r = 0} = -ky, see Eqs. (18) in 
                          ! Phys. Rev. AccelBeams 19.054002

  kappa = kx + (0, 1)*ky
  barkappa = conjg(kappa)

  ! Now fill up the g_{ij}'s for j = 0, ..., i and i = 0, ..., nord + 1.
  g(0, 0) = (0, 0)
  g(1, 0) = -lambda(0)
  g(1, 1) = conjg(g(1, 0))

  do k = 1, nord
     do j = 0, k - 1
        ! Eq. (6), in Ref. above
        g(k + 1, j + 1) = (barkappa*g(k, j + 1)*(j + one)*(j - k + three/two) +  &
             kappa*g(k, j)*(k - j)*(one/two - j))/(k - j)/(j + one)
     enddo
     ! Eq. (8) in Ref. above
     sum0 = 0
     do j = 1, k
       sum0 = sum0 - (k + 1 - j)*g(k + 1, j)*exp(-two*(0, 1)*j*tilt2)
     enddo
     g(k + 1, 0) = ( sum0 - two**k*exp(-(0, 1)*k*tilt2)*( lambda(k) &
                    + one/two*(barkappa*exp((0, 1)*tilt2) + kappa*exp(-(0, 1)*tilt2)) &
                    *lambda(k - 1) ) )/(k + one)
     g(k + 1, k + 1) = conjg(g(k + 1, 0))
  enddo

  do jtrk = 1, ktrack
     etahat = sqrt(two*track(6,jtrk)/beta + track(6,jtrk)**2 + one) - one ! etahat = deltap of individual particle
     h0 = sqrt((one + etahat)**2 - track(2,jtrk)**2 - track(4,jtrk)**2)

     rp = (track(1,jtrk) + (0, 1)*track(3,jtrk))/two
     rm = conjg(rp)

     ! Compute \partial_+ G using Eq. (7) in Ref. above     
     del_p_g = 0
     do k = 1, nord
        sum0 = 0
        do j = 0, k - 1
           sum0 = sum0 + (k - j)*g(k, j)*rp**(k - 1 - j)*rm**j
        enddo
        del_p_g = del_p_g + sum0
     enddo
     ! Now compute kick (Eqs. (38) in Ref. above)
     pkick = elrad*(barkappa*h0 + del_p_g)
     dpx = real(pkick)
     dpy = - aimag(pkick)
     track(1, jtrk) = track(1, jtrk) + elrad*(kx*track(1, jtrk) + ky*track(3, jtrk))*track(2, jtrk)/h0
     track(2, jtrk) = track(2, jtrk) + dpx
     track(3, jtrk) = track(3, jtrk) + elrad*(kx*track(1, jtrk) + ky*track(3, jtrk))*track(4, jtrk)/h0
     track(4, jtrk) = track(4, jtrk) + dpy
     ! N.B. orbit(5) = \sigma/beta and orbit(6) = beta*p_\sigma
     track(5, jtrk) = track(5, jtrk) - elrad*(kx*track(1, jtrk) + ky*track(3, jtrk)) &
                *(one + beta*track(6, jtrk))/beta/(one + etahat)
  enddo


end subroutine ttmult_cf

subroutine ttmult(track,ktrack,dxt,dyt,turn, thin_foc)
  use twtrrfi
  use name_lenfi
  use trackfi
  use time_varfi
  use math_constfi, only : zero, one, two, three
  use track_enums
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !    Track particle through a general thin multipole.                  *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) Number of surviving tracks.                    *
  !   dxt       (double)  local buffer                                   *
  !   dyt       (double)  local buffer                                   *
  !   el_num    (integer) elemenent number in sequence                   *
  !   para(*,20)(double)  matrix containing the values                   *
  !----------------------------------------------------------------------*
  double precision :: track(6,*), dxt(*), dyt(*)
  integer :: ktrack, turn

  logical, save :: first=.true.
  logical ::  time_var,thin_foc
  integer :: iord, jtrk, nd, nord, i, j, n_ferr, nn, ns, noisemax, nn1, in, mylen
  integer :: nnt, nst
  double precision :: curv, dbi, dbr, dipi, dipr, dx, dy, elrad
  double precision :: pt, px, py, rfac
  double precision :: f_errors(0:maxferr)
  double precision :: field(2,0:maxmul)
  !double precision :: vals(2,0:maxmul)
  double precision :: normal(0:maxmul), skew(0:maxmul),normalt(0:maxmul),skewt(0:maxmul),an
  double precision, save :: ordinv(maxmul), const
  double precision :: bvk, node_value, ttt
  double precision :: npeak(100), nlag(100), ntune(100), temp, noise
  character(len=name_len) name
  double precision :: beta_sqr, f_damp_t

  integer :: node_fd_errors, store_no_fd_err, get_option
  double precision , external:: get_tt_attrib  
  external:: get_tt_multipoles

  !---- Precompute reciprocals of orders and radiation constant
  if (first) then
     do iord = 1, maxmul
        ordinv(iord) = one / dble(iord)
     enddo
     const = arad * (betas * gammas)**3 / three
     first = .false.
  endif

  F_ERRORS(0:maxferr) = zero
  n_ferr = node_fd_errors(f_errors)

  bvk = get_tt_attrib(enum_other_bv)
    !---- Multipole length for radiation.
  elrad = get_tt_attrib(enum_lrad)
  noise = get_tt_attrib(enum_noise)
  an = get_tt_attrib(enum_angle)
  time_var = get_tt_attrib(enum_time_var) .ne. 0  

  
  !---- Multipole components.
  NORMAL(0:maxmul) = zero! ; call get_node_vector('knl ',nn,normal)
  SKEW(0:maxmul) = zero  ! ; call get_node_vector('ksl ',ns,skew)

  call get_tt_multipoles(nn,normal,ns,skew)


  nd = 2 * max(nn, ns, n_ferr/2-1)

  !---- Angle (no bvk in track)
  if (an .ne. 0) f_errors(0) = f_errors(0) + normal(0) - an

  !----
  if (noise .eq. 1)   then
     nn1 = name_len
     noisemax = node_value('noisemax ')
     ! 2015-Jun-11  12:37:29  ghislain: should insert a guard for noisemax < 100
     NPEAK(:noisemax) = zero ; call get_node_vector('npeak ',nn1,npeak)
     NTUNE(:noisemax) = zero ; call get_node_vector('ntune ',nn1,ntune)
     NLAG(:noisemax) = zero  ; call get_node_vector('nlag ',nn1,nlag)

     temp = 0
     do in = 1, noisemax
        temp = temp + npeak(in) * sin(nlag(in) + ntune(in) * turn)
     enddo
     NORMAL(:nn) = NORMAL(:nn) * (1+temp)
     SKEW(:nn)   = SKEW(:nn)   * (1+temp)
  endif

  !--- Time variation for fields in matrix, multipole or RF-cavity
  ! 2015-Jun-24  18:55:43  ghislain: DOC FIXME not documented!!!
 ! time_var = node_value('time_var ') .ne. zero

  if (time_var .and. time_var_m) then
     time_var_m_cnt = time_var_m_cnt + 1
     time_var_m_lnt = time_var_m_lnt + 1
     if (idnint(time_var_m_ind(time_var_m_cnt)) .ne. time_var_m_lnt)    &
          call fort_fail('TTMULT: ', 'wrong index in Table: time_var_mul')

     call element_name(name,len(name))
     mylen = len_trim(name)
     if (time_var_m_ch(time_var_m_cnt)(:mylen) .ne. name(:mylen)) &
          call fort_fail('TTMULT: ', 'wrong element name in Table: time_var_mul')

     !--- find maximum order for myfield
     do i = maxmul, 0, -1
        if (abs(myfield(time_var_m_cnt,1,i)) .gt. zero .or.   &
           abs(myfield(time_var_m_cnt,2,i)) .gt. zero) then
           n_ferr = i ! replacing the previous count of field-errors
           goto 101
        endif
     enddo

101  n_ferr = 2*n_ferr + 2
     do i=0,(n_ferr-2)/2
        f_errors(i)          = myfield(time_var_m_cnt,1,i)
        f_errors(n_ferr/2+i) = myfield(time_var_m_cnt,2,i)
     enddo
     nd = 2 * max(nn, ns, n_ferr/2-1)
     call dcopy(f_errors,field,nd+2)
     n_ferr = store_no_fd_err(f_errors,n_ferr)

  endif

  !-----added FrankS, 10-12-2008
  !nd = 2 * max(nn, ns, n_ferr/2-1)

  !---- Dipole error.
  !      dbr = bvk * field(1,0) / (one + deltas)
  !      dbi = bvk * field(2,0) / (one + deltas)
  dbr = bvk * f_errors(0) !field(1,0)
  dbi = bvk * f_errors(1) !field(2,0)

  !---- Nominal dipole strength.
  !      dipr = bvk * vals(1,0) / (one + deltas)
  !      dipi = bvk * vals(2,0) / (one + deltas)
  dipr = bvk * normal(0) !vals(1,0)
  dipi = bvk * skew(0)   !vals(2,0)

  !---- Other components and errors.
  nord = 0
  do iord = 1, nd/2
     f_errors(2*iord)   = bvk * (f_errors(2*iord) + normal(iord))
     f_errors(2*iord+1) = bvk * (f_errors(2*iord+1) + skew(iord))
     if (f_errors(2*iord).ne.zero .or. f_errors(2*iord+1).ne.zero) nord=iord
  enddo

  !---- Pure dipole: only quadrupole kicks according to lrad.
  if (nord .eq. 0) then
     dxt(:ktrack) = zero
     dyt(:ktrack) = zero
     !----------- introduction of dipole focusing
     if (elrad.gt.zero .and. thin_foc) then

        DXT(:ktrack) = dipr*dipr*TRACK(1,:ktrack)/elrad
        DYT(:ktrack) = dipi*dipi*TRACK(3,:ktrack)/elrad
     endif
  !---- Accumulate multipole kick from highest multipole to quadrupole.
  else
     DXT(:ktrack) = f_errors(2*nord)*TRACK(1,:ktrack) - f_errors(2*nord+1)*TRACK(3,:ktrack)
     DYT(:ktrack) = f_errors(2*nord)*TRACK(3,:ktrack) + f_errors(2*nord+1)*TRACK(1,:ktrack)

     do iord = nord - 1, 1, -1
        do jtrk = 1,ktrack
           dx = dxt(jtrk)*ordinv(iord+1) + f_errors(2*iord)
           dy = dyt(jtrk)*ordinv(iord+1) + f_errors(2*iord+1)
           dxt(jtrk) = dx*track(1,jtrk) - dy*track(3,jtrk)
           dyt(jtrk) = dx*track(3,jtrk) + dy*track(1,jtrk)
        enddo
     enddo
     if (elrad.gt.zero .and. thin_foc) then
        DXT(:ktrack) = DXT(:ktrack) + dipr*dipr*TRACK(1,:ktrack)/elrad
        DYT(:ktrack) = DYT(:ktrack) + dipi*dipi*TRACK(3,:ktrack)/elrad
     endif
  endif

  !---- Radiation loss at entrance.
  if (radiate .and. elrad .ne. 0) then
     !---- Full damping.
     if (damp) then
        do jtrk = 1,ktrack
           curv = sqrt((dipr + dxt(jtrk))**2 + (dipi + dyt(jtrk))**2) / elrad
           if (quantum) then
              call trphot(elrad,curv,rfac,pt)
           else
              rfac = const * curv**2 * elrad
           endif
           px = track(2,jtrk)
           py = track(4,jtrk)
           pt = track(6,jtrk)
           beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           track(2,jtrk) = px * f_damp_t;
           track(4,jtrk) = py * f_damp_t;
           track(6,jtrk) = pt * (one - rfac) - rfac / bet0;
        enddo
        !---- Energy loss like for closed orbit.
     else
        !---- Store energy loss on closed orbit.
        ! 2016-Mar-16  18:45:41  ghislain: track(i,1) is not the closed orbit but the first particle!!!
        rfac = const * ((dipr + dxt(1))**2 + (dipi + dyt(1))**2)
        pt = track(6,1)
        beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
        TRACK(2,:ktrack) = TRACK(2,:ktrack) * f_damp_t;
        TRACK(4,:ktrack) = TRACK(4,:ktrack) * f_damp_t;
        TRACK(6,:ktrack) = TRACK(6,:ktrack) * (one - rfac) - rfac / bet0;
     endif
  endif

  !---- Apply multipole effect including dipole.
  do jtrk = 1,ktrack
     !       Added for correct Ripken implementation of formulae
     ttt = sqrt( one + two*track(6,jtrk)*bet0i + track(6,jtrk)**2 )
     ! track(2,jtrk) = track(2,jtrk) - (dbr + dxt(jtrk) - dipr * (deltas + beti*track(6,jtrk)))
     ! track(4,jtrk) = track(4,jtrk) + (dbi + dyt(jtrk) - dipi * (deltas + beti*track(6,jtrk)))
     ! track(5,jtrk) = track(5,jtrk) - (dipr*track(1,jtrk) - dipi*track(3,jtrk)) * beti
     track(2,jtrk) = track(2,jtrk) - (dbr + dxt(jtrk) - dipr * (ttt - one))
     track(4,jtrk) = track(4,jtrk) + (dbi + dyt(jtrk) - dipi * (ttt - one))
     track(5,jtrk) = track(5,jtrk) - &
          (dipr*track(1,jtrk) - dipi*track(3,jtrk)) *   &
          ((one + bet0*track(6,jtrk))/ttt) * bet0i
  enddo

  !---- Radiation loss at exit.
  if (radiate .and. elrad .ne. 0) then
     !---- Full damping.
     if (damp) then
        do jtrk = 1,ktrack
           curv = sqrt((dipr + dxt(jtrk))**2 + (dipi + dyt(jtrk))**2) / elrad
           if (quantum) then
              call trphot(elrad,curv,rfac,pt)
           else
              rfac = const * curv**2 * elrad
           endif
           px = track(2,jtrk)
           py = track(4,jtrk)
           pt = track(6,jtrk)
           beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           track(2,jtrk) = px * f_damp_t;
           track(4,jtrk) = py * f_damp_t;
           track(6,jtrk) = pt * (one - rfac) - rfac / bet0;
        enddo

        !---- Energy loss like for closed orbit.
     else

        !---- Store energy loss on closed orbit.
        rfac = const * ((dipr + dxt(1))**2 + (dipi + dyt(1))**2)
        pt = track(6,1)
        beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
        TRACK(2,:ktrack) = TRACK(2,:ktrack) * f_damp_t;
        TRACK(4,:ktrack) = TRACK(4,:ktrack) * f_damp_t;
        TRACK(6,:ktrack) = TRACK(6,:ktrack) * (one - rfac) - rfac / bet0;

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
  double precision :: track(6,*)
  integer :: ktrack

  double precision :: psi, ct, st, trb(4)
  integer :: i

  double precision, external :: node_value

  psi = node_value('angle ')
  ct = cos(psi)
  st = -sin(psi)
  do  i = 1, ktrack
     TRB = TRACK(1:4,i)
     track(1,i) = trb(1)*ct - trb(3)*st
     track(2,i) = trb(2)*ct - trb(4)*st
     track(3,i) = trb(1)*st + trb(3)*ct
     track(4,i) = trb(2)*st + trb(4)*ct
  enddo
end subroutine ttsrot

subroutine ttxrot(track,ktrack)
  use trackfi
  use math_constfi, only : one, two
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Coordinate change due to rotation about y axis                     *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.        tmsrot      *
  !----------------------------------------------------------------------*
  double precision :: track(6,*)
  integer :: ktrack

  integer :: i
  double precision :: angle, ca, sa, ta
  double precision :: x, px, y, py, t, pt, pz, ptt
  double precision :: node_value

  angle = node_value('angle ')
  if (angle .eq. 0) return

  angle = angle * node_value('other_bv ')
  ca = cos(angle)
  sa = sin(angle)
  ta = tan(angle)
  do i = 1, ktrack
    x  = TRACK(1,i)
    px = TRACK(2,i)
    y  = TRACK(3,i)
    py = TRACK(4,i)
    t  = TRACK(5,i)
    pt = TRACK(6,i)
    pz = sqrt(one + two*pt/bet0i + pt**2 - px**2 - py**2)
    ptt = 1 - ta*py/pz

    track(3,i) = y/(ca*ptt)
    track(4,i) = ca*py + sa*pz
    track(1,i) = x + ta*y*px/(pz*ptt)
    track(5,i) = t - ta*y*(one/bet0i+pt)/(pz*ptt)
  enddo
end subroutine ttxrot

subroutine ttyrot(track,ktrack)
  use trackfi
  use math_constfi, only : one, two
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Coordinate change due to rotation about y axis                     *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !   KTRACK    (integer) number of surviving tracks.        tmsrot      *
  !----------------------------------------------------------------------*
  double precision :: track(6,*)
  integer :: ktrack

  integer :: i
  double precision :: angle, ca, sa, ta
  double precision :: x, px, y, py, t, pt, pz, ptt
  double precision :: node_value

  angle = node_value('angle ')
  if (angle .eq. 0) return

  angle = angle * node_value('other_bv ')
  ca = cos(angle)
  sa = sin(angle)
  ta = tan(angle)
  do i = 1, ktrack
    x  = TRACK(1,i)
    px = TRACK(2,i)
    y  = TRACK(3,i)
    py = TRACK(4,i)
    t  = TRACK(5,i)
    pt = TRACK(6,i)

    pz = sqrt(one + two*pt/bet0i + pt**2 - px**2 - py**2)
    ptt = 1 - ta*px/pz
    track(1,i) = x/(ca*ptt)
    track(2,i) = ca*px + sa*pz
    track(3,i) = y + ta*x*py/(pz*ptt)
    track(5,i) = t - ta*x*(one/bet0i+pt)/(pz*ptt)
  enddo
end subroutine ttyrot

subroutine ttdrf(el,track,ktrack)
  use trackfi
  use math_constfi, only : one, two
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
  double precision :: el, track(6,*)
  integer :: ktrack

  double precision :: pt, px, py, l_pz
  integer :: i

! picked from trturn in madx.ss
!$OMP PARALLEL PRIVATE(i, px, py, pt, l_pz)
!$OMP DO
  do  i = 1, ktrack
     px = track(2,i)
     py = track(4,i)
     pt = track(6,i)
     ! L/pz
     l_pz = el / sqrt( one + two*pt*bet0i + pt**2 - px**2 - py**2)
     track(1,i) = track(1,i) + l_pz*px
     track(3,i) = track(3,i) + l_pz*py
     ! track(5,i) = track(5,i) + el*(beti + pt*dtbyds) - (beti+pt)*l_pz
     !---- AK 20060413
     !---- Ripken DESY-95-189 p.36
     track(5,i) = track(5,i) + bet0i*(el - (one + bet0*pt) * l_pz)
  enddo
!$OMP END DO
!$OMP END PARALLEL
end subroutine ttdrf

subroutine ttrf(track,ktrack)
  use twtrrfi
  use name_lenfi
  use time_varfi
  use trackfi
  use math_constfi, only : zero, one, two, ten3m, ten6p, twopi, half,pi
  use phys_constfi, only : clight
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
  double precision :: track(6,*)
  integer :: ktrack, itrack, get_option
  logical :: time_var, Fringe
  integer :: i, mylen
  double precision :: omega, phirf, pt, rff, rfl, rfv, vrf, pc0, bvk
  double precision :: bucket_half_length, dummy, tcorr
  double precision :: el, V, c1(ktrack), s1(ktrack), jc
  !      double precision px,py,ttt,beti,el1
  character(len=name_len) :: name

  double precision :: get_variable, node_value, get_value

  !---- BV flag
  bvk = node_value('other_bv ')


  rfv = bvk * node_value('volt ')

  !---- Time Variation
  time_var = node_value('time_var ') .ne. zero
  if (time_var.and.time_var_c) then
     time_var_c_cnt = time_var_c_cnt + 1
     time_var_c_lnt = time_var_c_lnt + 1
     if (idnint(time_var_c_ind(time_var_c_cnt)) .ne. time_var_c_lnt)    &
          call fort_fail('TTRF: ', 'wrong index in Table: time_var_cav')
     call element_name(name,len(name))
     mylen = len_trim(name)
     if (time_var_c_ch(time_var_c_cnt)(:mylen) .ne. name(:mylen))       &
          call fort_fail('TTRF: ', 'wrong element name in Table: time_var_cav')
     !---- Overwrite Volt
     rfv = cav_volt(time_var_c_cnt)
     !---- Store Volt
     call store_node_value('volt ',rfv)
  endif

  rff = node_value('freq ')
  rfl = node_value('lag ')+node_value('lagtap ')
  pc0 = get_value('beam ','pc ')
  omega = rff * (ten6p * twopi / clight)

  ! vrf   = rfv * ten3m / (pc * (one + deltas))
  vrf   = rfv * ten3m
  phirf = rfl * twopi
  ! dl    = el / two
  ! bi2gi2 = one / (betas * gammas) ** 2
  el  = node_value('l ')
  fringe = node_value('fringe ') .gt. zero
  if(el .gt. zero) then
    if(fringe) then
      tcorr = el/(2*betas)
      jc = one
      V = jc*vrf/(pc0*el)
      s1=sin(phirf - omega*(TRACK(5,1:ktrack)+tcorr))
      c1=cos(phirf - omega*(TRACK(5,1:ktrack)+tcorr))

      TRACK(2,1:ktrack)=TRACK(2,1:ktrack)-V*S1*(TRACK(1,1:ktrack))*half
      TRACK(4,1:ktrack)=TRACK(4,1:ktrack)-V*S1*(TRACK(3,1:ktrack))*half
      TRACK(6,1:ktrack)=TRACK(6,1:ktrack)+0.25d0*(TRACK(1,1:ktrack)**2+TRACK(3,1:ktrack)**2)*V*c1*omega
    endif
    call ttdrf(el/2,track,ktrack)
  endif

  TRACK(6,1:ktrack) = TRACK(6,1:ktrack) +  vrf * sin(phirf - omega*TRACK(5,1:ktrack)) / pc0

  if(el .gt. zero) then
    call ttdrf(el/2,track,ktrack)
    if(fringe) then
      jc = -one
      V = jc*vrf/(pc0*el)
      s1=sin(phirf - (omega)*(TRACK(5,1:ktrack)+jc*tcorr))
      c1=cos(phirf - (omega)*(TRACK(5,1:ktrack)+jc*tcorr))

      TRACK(2,1:ktrack)=TRACK(2,1:ktrack)-V*S1*(TRACK(1,1:ktrack))*half
      TRACK(4,1:ktrack)=TRACK(4,1:ktrack)-V*S1*(TRACK(3,1:ktrack))*half
      TRACK(6,1:ktrack)=TRACK(6,1:ktrack)+0.25d0*(TRACK(1,1:ktrack)**2+TRACK(3,1:ktrack)**2)*V*c1*omega
    endif
  endif

   if(get_option('bucket_swap ').eq.1) then
    bucket_half_length = &
      get_value('probe ','circ ') / (two * get_value('probe ', 'beta ') &
                                     * node_value('harmon '));
    do itrack = 1, ktrack
      if (abs(track(5,itrack)) .gt. bucket_half_length) then
         dummy = mod((abs(track(5,itrack)) + bucket_half_length), bucket_half_length*2)
         track(5,itrack)=sign(1d0,track(5,itrack))*(dummy-bucket_half_length)
      endif
    enddo
  endif
  !! frs add-on end
 
end subroutine ttrf

subroutine ttchangep0(track,ktrack)
  use math_constfi, only : zero, two, one
  use phys_constfi, only : clight
  implicit none
  double precision :: track(6,*) 
  double precision :: get_value, bet0
  double precision :: pc0, px_, py_, pt_, onedp
  integer :: i, ktrack

  pc0 = get_value('beam ','pc ')
  bet0 = get_value('beam ','beta ')
  do i =1, ktrack
    px_ = track(1,i)
    py_ = track(3,i)
    pt_ = track(6,i) 
    
    onedp   = sqrt( one + two*pt_/bet0 + (pt_**2))

    TRACK(2,i) = TRACK(2,i)/onedp
    TRACK(4,i) = TRACK(4,i)/onedp
    TRACK(6,i) = zero
  end do

end subroutine ttchangep0



subroutine ttcrabrf(track,ktrack,turn)
  use math_constfi, only : zero, ten3m, ten6p, twopi
  use phys_constfi, only : clight
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
  double precision :: track(6,*)
  integer :: ktrack, turn

  integer :: i, t1, t2, t3, t4, p1, p2, bvk
  double precision :: omega, phirf, pt, rff, rfl, rfv, eph, vrf, pc0,  px
  !      double precision px,py,ttt,beti,el1

  double precision :: node_value, get_value

  !---- Initialize
  bvk = node_value('other_bv ')

  !---- Fetch data.
  rfv = bvk * node_value('volt ')
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
     vrf = zero
  else if (turn .lt. t2) then
     vrf = (turn-t1) * rfv * ten3m/pc0/(t2-t1)
  else if (turn .lt. t3) then
     vrf = rfv * ten3m/pc0
  else if (turn .lt. t4) then
     vrf = (t4-turn) * rfv * ten3m/pc0/(t4-t3)
  else
     vrf = zero
  endif

  !---- Set up phase ramp.
  if (turn .lt. p1)  then
     phirf = rfl * twopi
  else if (turn .lt. p2) then
     phirf = (turn-p1) * eph * twopi/(p2-p1)
  else
     phirf = eph * twopi
  endif

  TRACK(2,1:ktrack) = TRACK(2,1:ktrack) + vrf*sin(phirf - bvk*omega*TRACK(5,1:ktrack))

  TRACK(6,1:ktrack) = TRACK(6,1:ktrack) - &
       omega * vrf * TRACK(1,1:ktrack) * cos(phirf - bvk*omega*TRACK(5,1:ktrack))
  
 
end subroutine ttcrabrf

subroutine tthacdip(track,ktrack,turn)
  use math_constfi, only : zero, ten3m, ten6p, twopi
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
  double precision, intent(INOUT) :: track(6,*)
  integer, intent(IN) :: ktrack, turn

  integer :: turn1, turn2, turn3, turn4
  double precision :: omega, phirf, rff, rfl, rfv, vrf, pc0

  double precision :: node_value, get_value

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

  if (turn .lt. turn1)  then ! votage stable at zero
     vrf = zero
  else if (turn .lt. turn2) then ! ramping up the voltage
     vrf = (turn-turn1) * vrf / (turn2-turn1)
  else if (turn .lt. turn3) then ! voltage stable at maximum
     vrf = vrf
  else if (turn .lt. turn4) then ! ramping down the voltage
     vrf = (turn4-turn) * vrf / (turn4-turn3)
  else ! stable again at zero
     vrf = zero
  endif

  TRACK(2,:ktrack) = TRACK(2,:ktrack) + vrf * sin(phirf + omega * turn)

end subroutine tthacdip

subroutine ttvacdip(track,ktrack,turn)
  use math_constfi, only : zero, ten3m, ten6p, twopi
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
  double precision, intent(INOUT) :: track(6,*)
  integer, intent(IN) :: ktrack, turn

  integer :: turn1, turn2, turn3, turn4
  double precision :: omega, phirf, rff, rfl, rfv, vrf, pc0

  double precision :: node_value, get_value

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

  if (turn .lt. turn1)  then ! voltage stable at zero
     vrf = zero
  else if (turn .lt. turn2) then ! ramping up the voltage
     vrf = (turn-turn1) * vrf / (turn2-turn1)
  else if (turn .lt. turn3) then ! voltage stable at maximum
     vrf = vrf
  else if (turn .lt. turn4) then ! ramping down the voltage
     vrf = (turn4-turn) * vrf / (turn4-turn3)
  else ! stable again at zero
     vrf = zero
  endif

  TRACK(4,:ktrack) = TRACK(4,:ktrack) + vrf * sin(phirf + omega * turn)

end subroutine ttvacdip


subroutine ttsep(track,ktrack)
  use twtrrfi
  use trackfi
  use math_constfi, only : one, ten3m
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
  double precision :: track(6,*)
  integer :: ktrack

  integer :: i
  double precision :: el, ex, ey, tilt, charge, mass, pc, beta, deltap, kick0
  double precision :: pt

  double precision :: node_value, get_value

  ex = node_value('ex_l ')
  ey = node_value('ey_l ')

  charge = get_value('probe ','charge ')
  pc     = get_value('probe ','pc ')
  beta   = get_value('probe ','beta ')

  !Before there was rotation here. This made rotation not work because any rotation was later cancelled.

  do i = 1, ktrack
     deltap = sqrt(one - one/beta/beta + (track(6,i)+one/beta)**2) - one
     kick0 = charge * ten3m / pc / (one+deltap) / beta
     track(2,i) = track(2,i) + kick0*ex
     track(4,i) = track(4,i) + kick0*ey
  end do

end subroutine ttsep

subroutine ttcorr(el,track,ktrack,turn, code)
  use twtrrfi
  use trackfi
  use math_constfi, only : zero, one, two, three, twopi
  use code_constfi
  use track_enums
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
  double precision :: el, track(6,*)
  integer :: ktrack, turn

  integer :: i, n_ferr, code, bvk, sinkick
  double precision :: bi2gi2, bil2, curv, dpx, dpy, pt, px, py, rfac, rpt
  double precision :: rpx, rpy, xkick, ykick, div, temp
  double precision :: f_errors(0:maxferr), field(2)
  double precision :: dpxx, dpyy
  double precision :: sinpeak, sintune, sinphase
  double precision :: beta_sqr, f_damp_t

  integer :: node_fd_errors, get_option
  double precision :: get_variable, get_value, node_value
  double precision :: external, get_tt_attrib

  !---- Initialize.
        
      
  bvk = get_tt_attrib(enum_other_bv)
  sinkick = get_tt_attrib(enum_sinkick)
  

  
  !deltas = get_variable('track_deltap ')
  
  !arad = get_value('probe ','arad ')
  !betas = get_value('probe ','beta ')
  !gammas = get_value('probe ','gamma ')
  !dtbyds = get_value('probe ','dtbyds ')
  !radiate = get_value('probe ','radiate ') .ne. zero
  
  !damp = get_option('damp ') .ne. 0
  !quantum = get_option('quantum ') .ne. 0

  
!  if (code .eq. code_tkicker)      code = code_kicker
  !if (code .eq. code_placeholder) code = code_instrument

  F_ERRORS(0:maxferr) = zero ; n_ferr = node_fd_errors(f_errors)

  if (el .eq. zero)  then
     div = one
  else
     div = el
  endif

  FIELD(1:2) = zero

  select case (code)
    case (code_hkicker)
       xkick = bvk*(get_tt_attrib(enum_kick)+get_tt_attrib(enum_chkick)+field(1)/div)
       ykick = zero
    case (code_kicker, code_tkicker)
       xkick = bvk*(get_tt_attrib(enum_hkick)+get_tt_attrib(enum_chkick)+field(1)/div)
       ykick = bvk*(get_tt_attrib(enum_vkick)+get_tt_attrib(enum_cvkick)+field(2)/div)
    case (code_vkicker)
       xkick = zero
       ykick = bvk*(get_tt_attrib(enum_kick)+get_tt_attrib(enum_cvkick)+field(2)/div)
    case default
       xkick = zero
       ykick = zero
  end select

  !---- Sinusoidal kick (not supported by tkicker)
  if (sinkick .eq. 1) then
     sinpeak = node_value('sinpeak ')
     sintune = node_value('sintune ')
     sinphase = node_value('sinphase ')
     select case (code)
       case (code_hkicker)
          xkick = xkick + sinpeak * sin(sinphase + twopi * sintune * turn)
       case (code_kicker)
          xkick = xkick + sinpeak * sin(sinphase + twopi * sintune * turn)
          ykick = ykick + sinpeak * sin(sinphase + twopi * sintune * turn)
       case (code_vkicker)
          ykick = ykick + sinpeak * sin(sinphase + twopi * sintune * turn)
     end select
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
  if (radiate  .and.  el .ne. 0) then
     if (damp .and. quantum) then
        curv = sqrt(dpx**2 + dpy**2) / el
     else
        rfac = arad * (betas * gammas)**3 * (dpx**2 + dpy**2) / (three * el)
     endif

     if (damp) then !---- Full damping.
        do i = 1, ktrack
           px = track(2,i)
           py = track(4,i)
           pt = track(6,i)
           if (quantum) call trphot(el,curv,rfac,pt)
           beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           track(2,i) = px * f_damp_t;
           track(4,i) = py * f_damp_t;
           track(6,i) = pt * (one - rfac) - rfac / bet0;
        enddo
     else !---- Energy loss as for closed orbit.
        pt = track(6,1)
        beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
        track(2,:ktrack) = track(2,:ktrack) * f_damp_t;
        track(4,:ktrack) = track(4,:ktrack) * f_damp_t;
        track(6,:ktrack) = track(6,:ktrack) * (one - rfac) - rfac / bet0
     endif
  endif

  !---- Thick lens code taken out! AK 21.04.2006
  !!---- Half kick at entrance.
  !      do i = 1, ktrack
  !        px = track(2,i) + dpx/two
  !        py = track(4,i) + dpy/two
  !        pt = track(6,i)
  !
  !!---- Drift through corrector.
  !        d = (one - pt / betas) * el
  !        track(1,i) = track(1,i) + px * d
  !        track(3,i) = track(3,i) + py * d
  !        track(5,i) = track(5,i) + d * bi2gi2 * pt -           &
  !     bil2 * (px**2 + py**2 + bi2gi2*pt**2) + el*pt*dtbyds
  !
  !!---- Half kick at exit.
  !        track(2,i) = px + dpx/two
  !        track(4,i) = py + dpy/two
  !        track(6,i) = pt
  !      enddo

  !---- Kick at dipole corrector magnet
  !     including PT-dependence
  if(el.gt.zero) then
    call ttdrf(el/two,track,ktrack); !Tracks to the middle
  endif
  do i = 1, ktrack
     px = track(2,i)
     py = track(4,i)
     pt = track(6,i)
     !        ttt = sqrt(one+two*pt*bet0i+pt**2 - px**2 - py**2)
     !        ddd = sqrt(one+two*pt*bet0i+pt**2)
     !        track(2,i) = px + dpxx*ttt/ddd
     !        track(4,i) = py + dpyy*ttt/ddd
     !        ttt = sqrt(one+two*pt*bet0i+pt**2 - dpxx**2 - dpyy**2)
     !        ddd = sqrt(one+two*pt*bet0i+pt**2)
     !        track(2,i) = px + dpxx*ttt/ddd
     !        track(4,i) = py + dpyy*ttt/ddd

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
     !        track(2,i) = px
     !        track(4,i) = py

     track(2,i) = px + dpxx
     track(4,i) = py + dpyy

     !       Add time of flight effects (stolen from Ripken-Dipole)
     !        track(5,i) = track(5,i) -                             &
     !        (dpxx*track(1,i) - dpyy*track(3,i)) *                &
     !        ((one + bet0*track(6,i))/ddd)*bet0i

  enddo
  if(el .gt. zero) then
    call ttdrf(el/two,track,ktrack); !Tracks from the middle to the
  endif
  !---- Half radiation effects at exit.
  !     If not random, use same RFAC as at entrance.
  if (radiate  .and.  el .ne. 0) then
     if (damp) then !---- Full damping.
        do i = 1, ktrack
           pt = track(6,i)
           if (quantum) call trphot(el,curv,rfac,pt)
           beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           track(2,i) = track(2,i) * f_damp_t;
           track(4,i) = track(4,i) * f_damp_t;
           track(6,i) = track(6,i) * (one - rfac) - rfac / bet0;
        enddo
     else !---- Energy loss as for closed orbit.
        pt = track(6,1)
        beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
        f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
        track(2,:ktrack) = track(2,:ktrack) * f_damp_t;
        track(4,:ktrack) = track(4,:ktrack) * f_damp_t;
        track(6,:ktrack) = track(6,:ktrack) * (one - rfac) - rfac / bet0;
     endif
  endif

end subroutine ttcorr

subroutine ttbb(track,ktrack)
  use bbfi, only : explim
  use math_constfi, only : zero, one, two, half
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
  double precision :: track(6,*)
  integer :: ktrack

  integer :: beamshape, b_dir_int
  logical, save :: first=.true.
  logical :: bb_ultra_relati
  double precision :: parvec(26), fk, q, q_prime, dp
  double precision :: gamma0, beta0, beta_dp, ptot, b_dir

  integer :: get_option
  double precision :: get_value, node_value, get_variable
  double precision, parameter :: cme32=1d-32

  !---- Calculate momentum deviation and according changes
  !     of the relativistic factor beta0
  q = get_value('probe ','charge ')
  q_prime = node_value('charge ')
  parvec(5) = get_value('probe ', 'arad ')
  parvec(6) = node_value('charge ') * get_value('probe ', 'npart ')
  parvec(7) = get_value('probe ','gamma ')

  dp = get_variable('track_deltap ')
  gamma0 = parvec(7)
  beta0 = sqrt(one - one/gamma0**2)
  ptot = beta0*gamma0*(one+dp)
  beta_dp = ptot / sqrt(one + ptot**2)
  b_dir_int = node_value('bbdir ')
  b_dir = dble(b_dir_int)
  b_dir = b_dir/sqrt(b_dir*b_dir + cme32)

  !---- pre-factor, if zero, anything else does not need to be calculated
  bb_ultra_relati = get_option('bb_ultra_relati ') .ne. 0
  if (bb_ultra_relati) then
     fk = two * parvec(5) * parvec(6) / parvec(7)
  else
     fk = two*parvec(5)*parvec(6)/parvec(7)/beta0/(one+dp)/q*          &
       (one-beta0*beta_dp*b_dir)/(beta_dp+half*(b_dir-one)*b_dir*beta0)
  endif
  if (fk .eq. zero) return

  !---- choose beamshape: 1-Gaussian (default), 2-flattop=trapezoidal, 3-hollow-parabolic
  beamshape = node_value('bbshape ')
  select case (beamshape)
    case (1)
       call ttbb_gauss(track,ktrack,fk)
    case (2)
       call ttbb_flattop(track,ktrack,fk)
    case (3)
       call ttbb_hollowparabolic(track,ktrack,fk)
    case default
       if (first) then
         first = .false.
         call fort_warn('TTBB: ','beamshape out of range, set to default=1')
       endif
       call ttbb_gauss(track,ktrack,fk)
  end select

end subroutine ttbb

subroutine ttbb_gauss(track,ktrack,fk)
  use name_lenfi
  use bbfi
  use spch_bbfi
  use fasterror
  use math_constfi, only : zero, one, two, half, ten3m, pi
  implicit none
  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having a Gaussian perpendicular shape                       *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  double precision :: track(6,*), fk
  integer, intent(IN) :: ktrack

  integer :: i, ipos, mylen
  double precision :: sx, sy, xm, ym, sx2, sy2, xs, ys, rho2, tk
  double precision :: phix, phiy, rk, xb, yb, crx, cry, xr, yr, r, r2, cbx, cby
  double precision :: xrv(ktrack), yrv(ktrack), crxv(ktrack), cryv(ktrack)
  double precision :: xbv(ktrack), ybv(ktrack), cbxv(ktrack), cbyv(ktrack)
  double precision :: tkv(ktrack), phixv(ktrack), phiyv(ktrack)
  double precision :: xsv(ktrack), ysv(ktrack), rkv(ktrack)
  logical :: bborbit, bb_sxy_update, long_coup_on
  !VVK 20100321 ------------------------------------------------------
  double precision ::  gauss_factor_t ! promoted by GR from REAL
  character(len=20) :: text
  character(len=name_len) :: name
  !-------------------------------------------------------------------
  integer :: get_option
  double precision :: node_value

  !---- initialize.
  bborbit       = get_option('bborbit ') .ne. 0
  bb_sxy_update = get_option('bb_sxy_update ') .ne. 0
  long_coup_on  = get_option('long_coup_off ') .eq. 0
  if (bb_sxy_update) then
     name=' '
     call element_name(name,len(name))
     mylen = len_trim(name)
     i_spch = i_spch+1
     if (i_spch .gt. N_spch) then
        write(text, '(1p,i8)') i_spch
        call fort_fail('TTBB: ', 'Table with too few BB elements: '// text)
     endif
     if (spch_bb_name(i_spch)(:mylen) .ne. name(:mylen)) then
        call fort_fail('TTBB: ', 'wrong element name in Table: spch_bb')
     endif

     sx = sqrt(betx_bb(i_spch)*Ex_rms+(dx_bb(i_spch)*sigma_p)**2)
     sy = sqrt(bety_bb(i_spch)*Ey_rms+(dy_bb(i_spch)*sigma_p)**2)
  else
     sx = node_value('sigx ')
     sy = node_value('sigy ')
  endif

  xm = node_value('xma ')
  ym = node_value('yma ')
  if (bb_sxy_update) fk = fk * rat_bb_n_ions !Ratio_for_bb_N_ions

  if (fk .eq. zero)  return

  ipos = 0
  if (.not. bborbit)  then
     !--- find position of closed orbit bb_kick
     do ipos = 1, bbd_cnt
        if (bbd_loc(ipos) .eq. bbd_pos)  goto 1
     enddo
     ipos = 0
1    continue
  endif

!  print *, 'bborbit=', bborbit, &
!         ', ipos=', ipos, ', bbd_pos=', bbd_pos, ', bbd_cnt=', bbd_cnt, &
!         ', bb_kick_x=', bb_kick(1,ipos), ', bb_kick_y=', bb_kick(2,ipos)

  sx2 = sx*sx
  sy2 = sy*sy
  !---- limit formulae for sigma(x) = sigma(y).
  if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then
!$OMP PARALLEL PRIVATE(i, xs, ys, rho2, tk, gauss_factor_t, phix, phiy)
!$OMP DO
     do i = 1, ktrack
        xs = track(1,i) - xm
        ys = track(3,i) - ym
        rho2 = xs * xs + ys * ys
        tk = rho2 / (two * sx2)

        if (bb_sxy_update .and. long_coup_on) then
           gauss_factor_t = exp(-half*(track(5,i)-mean_t)**2/sigma_t**2)!VVK 20100321

           if (tk .gt. explim) then
              ! if x > explim, exp(-x) is outside machine limits.
              phix = xs * fk / rho2 * gauss_factor_t !VVK 20100321
              phiy = ys * fk / rho2 * gauss_factor_t !VVK 20100321
           else if (rho2 .ne. zero) then
              phix = xs * fk / rho2 * (one - exp(-tk) ) * gauss_factor_t !VVK 20100321
              phiy = ys * fk / rho2 * (one - exp(-tk) ) * gauss_factor_t !VVK 20100321
           else
              phix = zero
              phiy = zero
           endif
        else
           if (tk .gt. explim) then
              ! if x > explim, exp(-x) is outside machine limits.
              phix = xs * fk / rho2
              phiy = ys * fk / rho2
           else if (rho2 .ne. zero) then
              phix = xs * fk / rho2 * (one - exp(-tk) )
              phiy = ys * fk / rho2 * (one - exp(-tk) )
           else
              phix = zero
              phiy = zero
           endif
        endif

        if (ipos .ne. 0) then
           !--- subtract closed orbit kick
           phix = phix - bb_kick(1,ipos)
           phiy = phiy - bb_kick(2,ipos)
        endif
        track(2,i) = track(2,i) + phix
        track(4,i) = track(4,i) + phiy
     enddo
!$OMP END DO
!$OMP END PARALLEL

     !---- case sigma(x) > sigma(y).
  else if (sx2 .gt. sy2) then
     r2 = two * (sx2 - sy2)
     r  = sqrt(r2)
     ! rk = fk * sqrt(pi) / r                 !VVK 20100321
     rk = fk * sqrt(pi) / r
     if (fasterror_on) then
!$OMP PARALLEL PRIVATE(i, gauss_factor_t)
!$OMP DO
        do i = 1, ktrack

           if (bb_sxy_update) then
              gauss_factor_t = exp(-half*(track(5,i)-mean_t)**2/sigma_t**2)!VVK 20100321
              rkv(i) = fk * sqrt(pi) / r*gauss_factor_t !VVK 20100321
           endif

           xsv(i) = track(1,i) - xm
           ysv(i) = track(3,i) - ym
           xrv(i) = abs(xsv(i)) / r
           yrv(i) = abs(ysv(i)) / r
           tkv(i) = (xsv(i) * xsv(i) / sx2 + ysv(i) * ysv(i) / sy2) / two
           xbv(i) = (sy / sx) * xrv(i)
           ybv(i) = (sx / sy) * yrv(i)
        enddo
!$OMP END DO
!$OMP END PARALLEL
        call wzsubv(ktrack,xrv, yrv, crxv, cryv)
        call wzsubv(ktrack,xbv, ybv, cbxv, cbyv)
!        do i = 1, ktrack
!           if (tkv(i) .gt. explim) then
               ! if x > explim, exp(-x) is outside machine limits.
!              phixv(i) = rk * cryv(i)
!              phiyv(i) = rk * crxv(i)
!           endif
!        enddo
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
        do i = 1, ktrack
           phixv(i) = rkv(i) * (cryv(i) - exp(-tkv(i)) * cbyv(i))
           phiyv(i) = rkv(i) * (crxv(i) - exp(-tkv(i)) * cbxv(i))
           track(2,i) = track(2,i) + phixv(i) * sign(one,xsv(i))
           track(4,i) = track(4,i) + phiyv(i) * sign(one,ysv(i))
           if (ipos .ne. 0)  then
              !--- subtract closed orbit kick
              track(2,i) = track(2,i) - bb_kick(1,ipos)
              track(4,i) = track(4,i) - bb_kick(2,ipos)
           endif
        enddo
!$OMP END DO
!$OMP END PARALLEL
     else
!$OMP PARALLEL PRIVATE(i, rk, xs, ys, xr, yr, tk, xb, yb, crx, cry, cbx, cby, gauss_factor_t, phix, phiy)
!$OMP DO
        do i = 1, ktrack

           if (bb_sxy_update) then
              gauss_factor_t = exp(-half*(track(5,i)-mean_t)**2/sigma_t**2)!VVK 20100321
              rk = fk * sqrt(pi) / r*gauss_factor_t !VVK 20100321
           endif

           xs = track(1,i) - xm
           ys = track(3,i) - ym
           xr = abs(xs) / r
           yr = abs(ys) / r
           call ccperrf(xr, yr, crx, cry)
           tk = (xs * xs / sx2 + ys * ys / sy2) / two
           if (tk .gt. explim) then
              ! if x > explim, exp(-x) is outside machine limits.
              phix = rk * cry
              phiy = rk * crx
           else
              xb = (sy / sx) * xr
              yb = (sx / sy) * yr
              call ccperrf(xb, yb, cbx, cby)
              phix = rk * (cry - exp(-tk) * cby)
              phiy = rk * (crx - exp(-tk) * cbx)
           endif
           track(2,i) = track(2,i) + phix * sign(one,xs)
           track(4,i) = track(4,i) + phiy * sign(one,ys)
           if (ipos .ne. 0)  then
              !--- subtract closed orbit kick
              track(2,i) = track(2,i) - bb_kick(1,ipos)
              track(4,i) = track(4,i) - bb_kick(2,ipos)
           endif
        enddo
!$OMP END DO
!$OMP END PARALLEL

     endif

     !---- case sigma(x) < sigma(y).
  else
     r2 = two * (sy2 - sx2)
     r  = sqrt(r2)
     ! rk = fk * sqrt(pi) / r                 !VVK 20100321
     rk = fk * sqrt(pi) / r
!$OMP PARALLEL PRIVATE(i, rk, xs, ys, xr, yr, tk, xb, yb, crx, cry, cbx, cby, gauss_factor_t, phix, phiy)
!$OMP DO
     do i = 1, ktrack

        if (bb_sxy_update) then
           gauss_factor_t = exp(-half*(track(5,i)-mean_t)**2/sigma_t**2)!VVK 20100321
           rk = fk * sqrt(pi) / r*gauss_factor_t !VVK 20100321
        endif

        xs = track(1,i) - xm
        ys = track(3,i) - ym
        xr = abs(xs) / r
        yr = abs(ys) / r

        if (fasterror_on) then
           call wzsub(yr, xr, cry, crx)
        else
           call ccperrf(yr, xr, cry, crx)
        endif

        tk = (xs * xs / sx2 + ys * ys / sy2) / two
        if (tk .gt. explim) then
           ! if x > explim, exp(-x) is outside machine limits.
           phix = rk * cry
           phiy = rk * crx
        else
           xb  = (sy / sx) * xr
           yb  = (sx / sy) * yr

           if (fasterror_on) then
              call wzsub(yb, xb, cby, cbx)
           else
              call ccperrf(yb, xb, cby, cbx)
           endif

           phix = rk * (cry - exp(-tk) * cby)
           phiy = rk * (crx - exp(-tk) * cbx)
        endif
        track(2,i) = track(2,i) + phix * sign(one,xs)
        track(4,i) = track(4,i) + phiy * sign(one,ys)
        if (ipos .ne. 0)  then
           !--- subtract closed orbit kick
           track(2,i) = track(2,i) - bb_kick(1,ipos)
           track(4,i) = track(4,i) - bb_kick(2,ipos)
        endif
     enddo
!$OMP END DO
!$OMP END PARALLEL
  endif
end subroutine ttbb_gauss

subroutine ttbb_flattop(track,ktrack,fk)
  use bbfi
  use math_constfi, only : one, two, three, four, six, eight, twelve, half, ten3m, pi
  implicit none
  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having an trapezoidal and, so flat top radial profile       *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  double precision :: track(6,*), fk
  integer :: ktrack

  integer :: i, ipos
  double precision :: r0x, r0y, wi, wx, wy, xm, ym, r0x2, r0y2, xs, ys
  double precision :: rho, rho2, phir, phix, phiy, norm, r1, zz
  logical :: bborbit
  logical, save :: first= .true.

  integer :: get_option
  double precision :: node_value

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
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
  if (abs(r0x2 - r0y2) .gt. ten3m*(r0x2+r0y2)) then
     zz = half*(r0x + r0y)
     r0x = zz
     r0y = zz
     r0x2 = r0x*r0x
     r0y2 = r0y*r0y
     if (first) call fort_warn('TTBB_FLATTOP: ','beam is assumed to be circular')
     first=.false.
  endif

  norm = half*r0x**2 + wx**2/(six*four)
  r1 = r0x - wx/two
!$OMP PARALLEL PRIVATE(i, xs, ys, rho, rho2, phir, phix, phiy)
!$OMP DO
  do i = 1, ktrack
     xs = track(1,i) - xm
     ys = track(3,i) - ym
     rho2 = xs * xs + ys * ys
     rho  = sqrt(rho2)
     if (rho .le. r1) then
        phir = half/norm
        phix = phir*xs
        phiy = phir*ys
     else if (rho.gt.r1 .and. rho.lt.r1+wx) then
        phir = ( (r0x**2/four - r0x**3/six/wx - r0x*wx/eight + wx**2/(six*eight))/rho2 &
                + one/four + half*r0x/wx - rho/three/wx) / norm
        phix = phir*xs
        phiy = phir*ys
     else if (rho .ge. r1+wx) then
        phir = one/rho2
        phix = xs*phir
        phiy = ys*phir
     endif
     track(2,i) = track(2,i)+phix*fk
     track(4,i) = track(4,i)+phiy*fk
  end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine ttbb_flattop

subroutine ttbb_hollowparabolic(track,ktrack,fk)
  use bbfi
  use math_constfi, only : zero, one, two, three, four, twelve, half, ten3m, pi
  implicit none
  ! ---------------------------------------------------------------------*
  ! purpose: kicks the particles of the beam considered with a beam      *
  !          having a hollow-parabolic perpendicular shape               *
  ! input and output as those of subroutine ttbb                         *
  ! ---------------------------------------------------------------------*
  double precision :: track(6,*), fk
  integer :: ktrack

  integer :: i, ipos
  double precision :: r0x, r0y, wi, wx, wy, xm, ym, r0x2, r0y2, xs, ys
  double precision :: rho, rho2, phir, phix, phiy, norm, r1, zz
  logical :: bborbit
  logical, save :: first= .true.

  integer :: get_option
  double precision :: node_value

  !---- initialize.
  bborbit = get_option('bborbit ') .ne. 0
  ! mean radii of the is given via variables sigx and sigy
  r0x = node_value('sigx ')
  r0y = node_value('sigy ')
  wi = node_value('width ')
  ! width is given as FWHM of parabolic density profile, but formulas were
  ! derived with half width at the bottom of the parabolic density profile
  wi = wi/sqrt(two)
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
     zz = half*(r0x + r0y)
     r0x = zz
     r0y = zz
     r0x2 = r0x*r0x
     r0y2 = r0y*r0y
     if (first) call fort_warn('TTBB_HOLLOWPARABOLIC: ', 'beam is assumed to be circular')
     first = .false.
  endif

!$OMP PARALLEL PRIVATE(i, xs, ys, rho, rho2, phir, phix, phiy)
!$OMP DO
  do i = 1, ktrack
     xs = track(1,i) - xm
     ys = track(3,i) - ym
     rho2 = xs * xs + ys * ys
     rho  = sqrt(rho2)
     if (rho.le.r0x-wx) then
        phix = zero
        phiy = zero
     else if (rho.gt.r0x-wx .and. rho.lt.r0x+wx) then
        phir = three/four/wx/r0x/rho2 * (r0x**4/twelve/wx**2 - r0x**2/two + &
             two*r0x*wx/three - wx**2/four + rho2/two*(one - r0x**2/wx**2) + &
             rho**3/three*two*r0x/wx**2 - rho**4/four/wx**2)
        phix = phir*xs
        phiy = phir*ys
     else
        phir = one/rho2
        phix = xs*phir
        phiy = ys*phir
     endif
     track(2,i) = track(2,i) + phix*fk
     track(4,i) = track(4,i) + phiy*fk
  end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine ttbb_hollowparabolic

subroutine trkill(n, turn, sum, ntrk, part_id, &
     last_turn, last_pos, last_orbit, z, aptype)
  use name_lenfi
  use trackfi
  implicit none

  integer :: n, turn, ntrk, part_id(*), last_turn(*)
  double precision :: z(6,*), last_pos(*), last_orbit(6,*), sum
  character(len=name_len) :: aptype

  integer :: i, j
  double precision :: torb(6)
  logical :: recloss, exit_loss_turn
  character(len=name_len) :: el_name

  integer, external :: get_option

  recloss = get_option('recloss ') .ne. 0
  exit_loss_turn = get_option('exit_loss_turn ') .ne. 0

  !!--- As elements might have a tilt we have to transform back
  !!--- into the original coordinate system!
  !      theta = node_value('tilt ')
  !      if (theta .ne. zero)  then
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
  TORB = Z(:,n)
  LAST_ORBIT(:,part_id(n)) = Z(:,n)

  call element_name(el_name,len(el_name))

  write(6,'(''particle #'',i6,'' lost turn '',i6,''  at pos. s ='', f10.2,'' element='',a,'' aperture ='',a)') &
       part_id(n),turn,sum,el_name,aptype
  print *,"   X=",z(1,n),"  Y=",z(3,n),"  T=",z(5,n)

  if (exit_loss_turn) then
     lost_in_turn = .true.
     is_lost = .true.
  endif

  if (recloss) call tt_ploss(part_id(n),turn,sum,torb,el_name)

  do i = n+1, ntrk
     part_id(i-1) = part_id(i)
     Z(:,i-1) = Z(:,i)
  enddo
  ntrk = ntrk - 1

end subroutine trkill

subroutine tt_ploss(npart,turn,spos,orbit,el_name)
  use name_lenfi
  implicit none
  !----------------------------------------------------------------------*
  !--- purpose: enter lost particle coordinates in table                 *
  !    input:                                                            *
  !    npart  (int)           particle number                            *
  !    turn   (int)           turn number                                *
  !    spos    (double)       s-coordinate when loss happens             *
  !    orbit  (double array)  particle orbit                             *
  !    orbit0 (double array)  reference orbit                            *
  !----------------------------------------------------------------------*
  integer :: npart, turn
  double precision :: spos, orbit(6)
  character(len=name_len) :: el_name

  integer :: j
  double precision :: tmp, tt, tn, energy
  character(len=120) :: table='trackloss'
  character(len=4) :: vec_names(7)
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt', 's' /

  double precision, external :: get_value

  tn = npart
  tt = turn

  energy = get_value('probe ','energy ')

  ! the number of the current particle
  call double_to_table_curr(table, 'number ', tn)
  ! the number of the current turn
  call double_to_table_curr(table, 'turn ', tt)

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
                     spos,ielem,el_name, onlyaver)
  use name_lenfi
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
  integer :: npart, turn, tot_segm, segment, part_id(*), ielem
  double precision :: z(6,*), orbit0(6)
  character(len=name_len) :: el_name

  logical, save :: first = .true.
  logical :: onlyaver
  integer :: i, j, length
  double precision :: tmp, tt, ss, spos, tmp_v(6)
  character(len=120) :: table = 'trackone', comment
  character(len=4) :: vec_names(7)
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt','s' /

  length = len(comment)
  segment = segment + 1

  write(comment, '(''#segment'',4i8,1X,A)') segment, tot_segm, npart, ielem, el_name
  if (first) call comment_to_table_curr(table, comment, length)
  tt = turn
  if(onlyaver) then
    call double_to_table_curr(table, 'turn ', tt)
    ss = -1.0
    call double_to_table_curr(table, 'number ', ss)
    do j = 1, 6
    tmp = 0
      do i = 1, npart
        tmp = tmp + (z(j,i) - orbit0(j))
      enddo
      call double_to_table_curr(table, vec_names(j), tmp/npart)
    enddo
    
    call double_to_table_curr(table,vec_names(7),spos)
    call augment_count(table)
  else
    do i = 1, npart
       call double_to_table_curr(table, 'turn ', tt)
       ss = part_id(i)
       call double_to_table_curr(table, 'number ', ss)
       do j = 1, 6
          tmp = z(j,i) - orbit0(j)
          call double_to_table_curr(table, vec_names(j), tmp)
       enddo
       call double_to_table_curr(table,vec_names(7),spos)
       call augment_count(table)
    enddo
  endif

end subroutine tt_putone

subroutine tt_puttab(npart,turn,nobs,orbit,orbit0,spos)
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
  integer, intent(IN) :: npart, turn, nobs
  double precision, intent(IN) :: orbit(6), orbit0(6), spos

  integer :: i
  double precision :: energy, tmp, tt, tn
  character(len=36) :: table='track.obs$$$$.p$$$$'
  character(len=4) :: vec_names(7)
  data vec_names / 'x', 'px', 'y', 'py', 't', 'pt', 's' /

  double precision, external :: get_value

  tt = turn
  tn = npart

  write(table(10:13), '(i4.4)') nobs
  write(table(16:19), '(i4.4)') npart

  energy = get_value('probe ','energy ')

  call double_to_table_curr(table, 'turn ', tt)
  call double_to_table_curr(table, 'number ', tn)
  call double_to_table_curr(table, 'e ', energy)
  do i = 1, 6
     tmp = orbit(i) - orbit0(i)
     call double_to_table_curr(table, vec_names(i), tmp)
  enddo
  call double_to_table_curr(table,vec_names(7),spos)
  call augment_count(table)
end subroutine tt_puttab

subroutine trcoll(apint,  aperture, offset, al_errors, maxaper, &
                 turn, sum, part_id, last_turn, last_pos, last_orbit, z, ntrk, debug)
  use twiss0fi
  use name_lenfi
  use Inf_NaN_Detection
  use math_constfi, only : zero, one, pi
  use aperture_enums
  implicit none
  ! 2015-Feb-20  18:46:05  ghislain: rewrite of trcoll
  ! 2015-Mar-09  14:50:37  ghislain: adapted to new racetrack parameter definition
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   test for aperture limit.                                           *
  ! input:                                                               *
  !   apint  (integer)   aperture type among                           *
  !                       circle,rectangle,ellipse,rectcircle,lhcscreen, *
  !                       rectellipse,racetrack,octagon                  *
  !   aperture  (array double) aperture parameters,                      *
  !                         depending on aperture type                   *
  !   offset    (array double) aperture offsets                          *
  !   al_errors (array double) alignment errors                          *
  !   maxaper   (array double) maximum particle extensions               *
  !                                                                      *
  !   turn      (integer)   current turn number.                         *
  !   sum       (double)    accumulated length.                          *
  !                                                                      *
  ! input/output:                                                        *
  !   part_id   (int array) particle identification list                 *
  !   last_turn (int array) storage for number of last turn              *
  !   last_pos  (dp. array) storage for last position (= sum)            *
  !   last_orbit(dp. array) storage for last orbit                       *
  !   z(6,*)    (double)    track coordinates: (x, px, y, py, t, pt).    *
  !   ntrk      (integer) number of surviving tracks.                    *
  !----------------------------------------------------------------------*
  character(len=name_len) :: apertype
  double precision :: aperture(*), offset(2), al_errors(align_max), maxaper(6)
  double precision :: sum, last_pos(*), last_orbit(6,*), z(6,*)
  integer :: turn, part_id(*), last_turn(*), ntrk, apint

  integer :: i, n, nn, nna
  double precision :: ap1, ap2, ap3, ap4, x, y!, pi
  logical :: lost, debug, is_custom

  integer, external :: get_option, inside_userdefined_geometry, is_custom_set
  double precision, parameter :: min_double=1.d-36

 ! debug = get_option('debug ') .ne. 0

  !--- First check aperture parameters
  ap1 = zero ; ap2 = zero ; ap3 = zero ; ap4 = zero

  select case (apint)

     case (ap_circle)
        ap1 = aperture(1)
        if (ap1.lt.min_double) then
           if (debug) print *, " zero or negative circle radius ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif

     case (ap_ellipse)
        ap1 = aperture(1)
        ap2 = aperture(2)
        if (ap1.lt.min_double) then
           if (debug) print *, " zero half-axis ellipse  ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif
        if (ap2.lt.min_double) then
           if (debug) print *, " zero half-axis ellipse  ", ap2, " replaced by default ", maxaper(3)
           ap2 = maxaper(3)
        endif

     case (ap_rectangle)
        ap1 = aperture(1)
        ap2 = aperture(2)
        if (ap1.lt.min_double) then
           if (debug) print *, " zero or negative horizontal extent ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif
        if (ap2.lt.min_double) then
           if (debug) print *, " zero or negative vertical extent ", ap2, " replaced by default ", maxaper(3)
           ap2 = maxaper(3)
        endif

     case (ap_lhcscreen, ap_rectcircle)
        !*****         test circle
        ap3 = aperture(3)
        if (ap3.lt.min_double) then
           if (debug) print *, " zero or negative circle radius ", ap3, " replaced by default ", maxaper(1)
           ap3 = maxaper(1)
        endif
        !*****         test rectangle
        ap1 = aperture(1)
        ap2 = aperture(2)
        if (ap1.lt.min_double) then
           if (debug) print *, " zero or negative horizontal extent ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif
        if (ap2.lt.min_double) then
           if (debug) print *, " zero or negative vertical extent ", ap2, " replaced by default ", maxaper(3)
           ap2 = maxaper(3)
        endif

     case (ap_rectellipse)
        !*****         test ellipse
        ap3 = aperture(3)
        ap4 = aperture(4)
        if (ap3.lt.min_double) then
           if (debug) print *, " zero or negative ellipse half axis ", ap3, " replaced by default ", maxaper(1)
           ap3 = maxaper(1)
        endif
        if (ap4.lt.min_double) then
           if (debug) print *, " zero or negative ellipse half axis ", ap4, " replaced by default ", maxaper(3)
           ap4 = maxaper(3)
        endif
        !*****         test rectangle
        ap1 = aperture(1)
        ap2 = aperture(2)
       if (ap1.lt.min_double) then
           if (debug) print *, " zero or negative horizontal extent ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif
        if (ap2.lt.min_double) then
           if (debug) print *, " zero or negative vertical extent ", ap2, " replaced by default ", maxaper(3)
           ap2 = maxaper(3)
        endif

     case (ap_racetrack)
        ap1 = aperture(1)
        ap2 = aperture(2)
        ap3 = aperture(3)
        ap4 = aperture(4)
        if (ap1.lt.min_double) then
           if (debug) print *, " zero or negative horizontal extent ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif
        if (ap2.lt.min_double) then
           if (debug) print *, " zero or negative vertical extent ", ap2, " replaced by default ", maxaper(3)
           ap2 = maxaper(3)
        endif
        if (ap3.lt.zero) then ! zero extent of rounded corned is allowed
           if (debug) print *, " negative horizontal semi-axis ", ap3, &
                "; horizontal semi-axis reset to horizontal extent", ap1
           ap3 = ap1
        endif
        if (ap4.lt.zero) then ! zero extent of rounded corner is allowed
           if (debug) print *, " negative vertical semi-axis ", ap4, &
                "; horizontal semi-axis reset to horizontal extent", ap2
           ap4 = ap2
        endif
        if (ap3.gt.ap1) then
           if (debug) print *, " horizontal semi-axis ", ap3, " is larger than horizontal extent ", ap1, &
                "; horizontal semi-axis reset to horizontal extent"
           ap3 = ap1
        endif
        if (ap4.gt.ap2) then
           if (debug) print *, " vertical semi-axis", ap4, " is larger than vertical extent ", ap2, &
                "; vertical semi-axis reset to vertical extent"
           ap4 = ap2
        endif

     case (ap_octagon)
        ! 2015-Feb-20  18:42:26  ghislain: added octagon shape
        ap1 = aperture(1)
        ap2 = aperture(2)
        ap3 = aperture(3)
        ap4 = aperture(4)
        if (ap1.lt.min_double) then
           if (debug) print *, " zero or negative horizontal extent ", ap1, " replaced by default ", maxaper(1)
           ap1 = maxaper(1)
        endif
        if (ap2.lt.min_double) then
           if (debug) print *, " zero or negative vertical extent ", ap2, " replaced by default ", maxaper(3)
           ap2 = maxaper(3)
        endif
        if (ap3.lt.zero .or. ap3.gt.pi/2) then
           if (debug) print *, " first angle is not in first quadrant ", ap3, " replaced by default ", zero
           ap3 = zero
        endif
        if (ap4.lt.zero .or. ap4.gt.pi/2) then
           if (debug) print *, " second angle is not in first quadrant ", ap4, " replaced by default ", pi/2
           ap2 = pi/2
        endif
        if (ap3.gt.ap4) then
           call fort_warn('trcoll:','octagon aperture: first and second angles inverted. exit from trcoll')
           return
        endif
     case(ap_custom)
        ap1 = aperture(1)
        ap2 = aperture(2)
     case(ap_custom_inter)
     ! Intenitionaly left blank. 


     case default
        ! add error case for un-identified aperture type;
        ! this INCLUDES the case of aperture data given in file with input APERTYPE=filename!
        call fort_warn('trcoll:','called with unknown aperture type. Ignored')

  end select

  !--- Then track through

  n = 1
10 continue

  do i = n, ntrk

     if (ISNAN(z(1,i)) .or. ISNAN(z(3,i))) then
        lost = .true.
        goto 99 ! lost...
     endif



     lost = .false.

     x = abs(z(1,i) - al_errors(11) - offset(1))
     y = abs(z(3,i) - al_errors(12) - offset(2))

     select case (apint)
     case (ap_circle)
        lost = (x/ap1)**2 + (y/ap1)**2 .gt. one

     case (ap_ellipse)
        lost = (x/ap1)**2 + (y/ap2)**2 .gt. one

     case (ap_rectangle)
        lost =  x .gt. ap1 .or. y .gt. ap2

     case (ap_lhcscreen, ap_rectcircle)
        lost = x .gt. ap1 .or. y .gt. ap2 .or. (x/ap3)**2 + (y/ap3)**2 .gt. one

     case (ap_rectellipse)
       lost =  x .gt. ap1 .or. y .gt. ap2 .or. (x/ap3)**2 + (y/ap4)**2 .gt. one

     case (ap_racetrack)
        ! 2015-Mar-09  15:05:39  ghislain: adapted to new racetrack parameter definition
        !*** case of racetrack: test outer rectangle (ap1,ap2) first
        !    then test ellipse for corner part.
        lost =  x .gt. ap1 .or. y .gt. ap2 .or. &
             ( x .gt. ap1-ap3 .and. y .gt. ap2-ap4 .and. &
               ((x-(ap1-ap3)) / ap3)**2 + ((y-(ap2-ap4)) / ap4)**2 .gt. one )

     case (ap_octagon)
        ! 2015-Feb-20  18:42:26  ghislain: added octagon shape
        !*** case of octagon: test outer rectangle (ap1,ap2) then test cut corner.
        lost =  x .gt. ap1 .or. y .gt. ap2 .or. &
             (ap2*tan(pi/2 - ap4) - ap1)*(y - ap1*tan(ap3)) - (ap2 - ap1*tan(ap3))*(x - ap1) .lt. zero
     case (ap_custom)
        lost =  x .gt. ap1 .or. y .gt. ap2 ! First checks the user defined rectangle
        if(lost) then
          x = z(1,i) - al_errors(11) - offset(1)
          y = z(3,i) - al_errors(12) - offset(2)
          lost = inside_userdefined_geometry(x,y) .eq. 0
      endif
      case(ap_custom_inter)
        lost = .true.
     case default

     end select
     if(lost) then
       is_custom = is_custom_set() .eq. 1
       if(is_custom) then
          x = z(1,i) - al_errors(11) - offset(1)
          y = z(3,i) - al_errors(12) - offset(2)
          lost = inside_userdefined_geometry(x,y) .eq. 0  
       endif
     endif
    

     if (.not. lost) then
        lost =  ISNAN(z(2,i)) .or. ISNAN(z(4,i))                                .or. &
                ISNAN(z(5,i)) .or. ISNAN(z(6,i))                                .or. &
                abs(z(1, i)) .gt. maxaper(1) .or.  abs(z(2, i)) .gt. maxaper(2) .or. &
                abs(z(3, i)) .gt. maxaper(3) .or.  abs(z(4, i)) .gt. maxaper(4) .or. &
                abs(z(5, i)) .gt. maxaper(5) .or.  abs(z(6, i)) .gt. maxaper(6)
     endif

     ! lose particle if it is outside aperture
99   if (lost) then
        n = i
        nna=name_len
        call node_string('apertype ',apertype,nna)
        call trkill(n, turn, sum, ntrk, part_id, last_turn, last_pos, last_orbit, z, apertype)
        if (ntrk .eq. 0) then
           call fort_warn('trcoll: ','Particle Number equals zero: exit from trcoll')
           return
        endif
        ! particle numbering has been reset by trkill, restart the loop with new paramaters.
        goto 10
     endif

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
  integer :: turn, part_id(*), last_turn(*), ntrk
  double precision :: sum, last_pos(*), last_orbit(6,*), z(6,*)

  integer :: i, n, nn
  double precision :: al_errors(align_max), offx, offy
  character(len=name_len) :: non_app="RF-Bucket"

  n = 1
10 continue
  do i = n, ntrk

     !---- Is particle outside the bucket?
     if (ISNAN(z(5,i)) .or. ISNAN(z(6,i))) goto 99
     if (abs(z(5,i)).gt.t_max .or. abs(z(6,i)).gt.pt_max) goto 99
     ! if particle is inside aperture, stop treatment and continue the loop
     cycle

     ! lose particle if it is outside aperture
99   n = i
     nn = name_len
     call trkill(n, turn, sum, ntrk, part_id, last_turn, last_pos, last_orbit, z, non_app)
     if (ntrk .eq. 0) then
        call fort_warn('ttrfloss: ', 'Particle Number equals zero: exit from ttrfloss')
        return
     endif

     goto 10 ! restart loop starting where we had left off (n=i)

  enddo
end subroutine ttrfloss

subroutine trinicmd(switch,orbit0,eigen,jend,z,turns,coords)
  use bbfi
  use trackfi
  use math_constfi, only : zero, one, twopi
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
  integer :: switch, jend, turns
  double precision :: orbit0(6), eigen(6,6), z(6,*), coords(6,0:turns,*)

  logical :: run=.false., dynap=.false.
  logical :: zgiv, zngiv
  integer :: j, k, kp, kq,  itype(23), jt
  double precision :: track(12), zstart(12), zn(6)
  double precision :: ex, ey, et, x, px, y, py, t, deltae
  double precision :: fx, phix, fy, phiy, ft, phit
  double precision :: deltax, deltap, phi

  double precision, external :: get_value, get_variable
  integer, external :: next_start

  ! 2015-Jul-08  19:16:53  ghislain: make code more readable
  run   = switch .eq. 1
  dynap = switch .eq. 2

  deltap = get_variable('track_deltap ')
  !---- Initialise orbit, emittances and eigenvectors etc.
  ex = get_value('probe ','ex ')
  ey = get_value('probe ','ey ')
  et = get_value('probe ','et ')
  bet0  = get_value('beam ','beta ')
  bet0i = one / bet0
  !-----get x add-on for lyaponuv calculation from dynap table
  deltax = get_value('dynap ', 'lyapunov ')

  j = 0

  do ! loop over particles to be STARTed
     jt  =  next_start(x,px,y,py,t,deltae,fx,phix,fy,phiy,ft,phit)

     if (jt .eq. 0)  return

     j = jt   ; if (dynap) j = 2*jt - 1
     jend = j ; if (dynap) jend = jend + 1

     !---- Get start coordinates
     track(1) = x
     track(2) = px
     track(3) = y
     track(4) = py
     track(5) = t
     ! track(6) = deltae
     !---- Here we absorb deltap into the PT variable
     track(6) = deltae + sqrt((one+deltap)**2+bet0i**2-one)-bet0i
     track(7) = fx
     track(8) = phix
     track(9) = fy
     track(10) = phiy
     track(11) = ft
     track(12) = phit

     do k = 1, 12
        if (abs(track(k)) .ne. zero) then
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
     zgiv  = .false.
     zngiv = .false.
     do k = 1, 6
        if (itype(k) .ne. 0)   zgiv  = .true.
        if (itype(k+6) .ne. 0) zngiv = .true.
        zstart(k) = track(k) &
             + sqrt(ex) * (eigen(k,1) * zn(1) + eigen(k,2) * zn(2))  &
             + sqrt(ey) * (eigen(k,3) * zn(3) + eigen(k,4) * zn(4))  &
             + sqrt(et) * (eigen(k,5) * zn(5) + eigen(k,6) * zn(6))
     enddo

     !--- keep initial coordinates for dynap
     if (dynap) COORDS(1:6,0,j) = ZSTART(1:6)

     if (zgiv .and. zngiv) & !---- Warn user about possible data conflict.
          call fort_warn('START: ','Absolute and normalized coordinates given, superposition used.')

     Z(1:6,j) = ORBIT0(1:6) + ZSTART(1:6)

     if (dynap) then !--- build up paired partners with deltax difference in x
        z(1,j+1)   = z(1,j)  + deltax
        Z(2:6,j+1) = Z(2:6,j)
        COORDS(1:6,0,j+1) = Z(1:6,j+1)
     endif

  end do ! Loop broken with jt.eq.0

end subroutine trinicmd

subroutine trbegn(rt,eigen)
  use matrices, only : EYE
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Initialize tracking mode; TRACK command execution.                 *
  !----------------------------------------------------------------------*
  double precision :: rt(6,6),eigen(6,6)

  double precision :: reval(6), aival(6)
  logical :: onepass

  logical :: m66sta
  integer :: get_option

  !---- Initialize
  onepass = get_option('onepass ') .ne. 0
  !---- One-pass system: Do not normalize.
  if (onepass) then
     EIGEN = EYE
  else
     !---- Find eigenvectors.
     if (m66sta(rt)) then
        call laseig(rt,reval,aival,eigen)
     else
        call ladeig(rt,reval,aival,eigen)
     endif
  endif
end subroutine trbegn

subroutine ttdpdg_map(track, ktrack, e1, h, hgap, fint, tilt)
  use matrices, only : EYE
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose: computes the effect of the dipole edge                      *
  ! Input/Output:  ktrack , number of surviving trajectories             *
  !                 track , trajectory coordinates                       *
  !----------------------------------------------------------------------*
  integer :: ktrack
  double precision :: track(6,*), e1, h, hgap, fint, tilt

  double precision :: corr, ek(6), rw(6,6), tw(6,6,6)

  double precision, external :: node_value

  EK = zero
  RW = EYE
  TW = zero

  corr = (h + h) * hgap * fint

  !---- Fringe fields effects computed from the TWISS routine tmfrng
  !     tmfrng returns the matrix elements rw(used) and tw(unused)
  !     No radiation effects as it is a pure thin lens with no lrad
  call tmfrng(.false., h, zero, e1, zero, zero, corr, rw, tw)
!  call tmtilt(.false., tilt, ek, rw, tw) !! frs add-off

  TRACK(2,1:ktrack) = TRACK(2,1:ktrack) + RW(2,1) * TRACK(1,1:ktrack)
  TRACK(4,1:ktrack) = TRACK(4,1:ktrack) + RW(4,3) * TRACK(3,1:ktrack)

  return
end subroutine ttdpdg_map

subroutine ttdpdg(track, ktrack)
  use matrices, only : EYE
  use math_constfi, only : zero
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose: computes the effect of the dipole edge                      *
  ! Input/Output:  ktrack , number of surviving trajectories             *
  !                 track , trajectory coordinates                       *
  !----------------------------------------------------------------------*
  integer :: ktrack
  double precision :: track(6,*)

  double precision :: e1, h, hgap, fint, tilt, corr, bvk
  double precision :: ek(6), rw(6,6), tw(6,6,6)

  double precision :: node_value

  EK = zero
  RW = EYE
  TW = zero

  e1 = node_value('e1 ')
  h = node_value('h ')
  bvk = node_value('other_bv ') !! frs add-on
  hgap = node_value('hgap ')
  fint = node_value('fint ')
  tilt = node_value('tilt ')
  h = bvk * h                   !! frs add-on

  call ttdpdg_map(track, ktrack, e1, h, hgap, fint, tilt);

end subroutine ttdpdg

subroutine trsol(track,ktrack,dxt,dyt)
  use trackfi, only : radiate, arad, damp, quantum, gammas
  use math_constfi, only : zero, half, one, two, three, four
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
  double precision :: track(6,*)
  integer :: ktrack

  integer :: i, step
  double precision :: dxt(*), dyt(*)
  double precision :: bet0
  double precision :: sk, skl, cosTh, sinTh, Q, R, Z
  double precision :: xf, yf, pxf, pyf, sigf, psigf, bvk
  double precision :: onedp, fpsig, fppsig

  double precision :: get_value, node_value

  double precision :: omega, length
  double precision :: x_, y_, z_, px_, py_, pt_
  double precision :: pxf_, pyf_
  double precision :: bet, length_, elrad
  double precision :: curv, const, rfac
  double precision :: beta_sqr, f_damp_t

  !---- Initialize.
  bet0 = get_value('probe ','beta ')

  !---- Get solenoid parameters
  elrad   = node_value('lrad ')
  bvk = node_value('other_bv ')
  sk  = bvk * node_value('ks ') / two
  length = node_value('l ')

  if (length.eq.zero) then

     skl = bvk * node_value('ksi ') / two

     !---- Loop over particles
     do  i = 1, ktrack
        do step = 1, 3
           !     Ripken formulae p.28 (3.35 and 3.36)
           xf    = track(1,i)
           yf    = track(3,i)
           psigf = track(6,i) / bet0

           !     We do not use a constant deltap!!!!! WE use full 6D formulae!
           onedp   = sqrt( one + two*psigf + (bet0**2)*(psigf**2) )
           fpsig   = onedp - one
           fppsig  = ( one + (bet0**2)*psigf ) / onedp

           ! Set up C,S, Q,R,Z
           cosTh = cos(skl/onedp)
           sinTh = sin(skl/onedp)
           Q = -skl * sk / onedp
           R = fppsig / (onedp**2) * skl * sk
           Z = fppsig / (onedp**2) * skl

           pxf  = track(2,i) + xf*Q
           pyf  = track(4,i) + yf*Q
           sigf = track(5,i)*bet0 - half*(xf**2 + yf**2)*R

           ! For radiation calculations (initial angles)
           dxt(i) = track(2,i);
           dyt(i) = track(4,i);

           ! final angles after solenoid
           pxf_ =  pxf * cosTh  +  pyf * sinTh;
           pyf_ = -pxf * sinTh  +  pyf * cosTh;

           ! kick received by particle
           dxt(i) = pxf_ - track(2,i);
           dyt(i) = pyf_ - track(4,i);

           !---- Radiation loss at entrance (step.eq.1) and exit (step.eq.3)
           if ((step.eq.1).or.(step.eq.3)) then
              if (radiate .and. elrad .gt. zero) then
                 !---- Full damping.
                 if (damp) then
                    curv = sqrt(dxt(i)**2 + dyt(i)**2) / elrad;
                    if (quantum) then
                       call trphot(elrad,curv,rfac,track(6,i))
                    else
                       const = arad * (bet0 * gammas)**3 / three
                       rfac = const * curv**2 * elrad
                    endif
                    pt_ = track(6,i)
                    beta_sqr = (pt_*pt_ + two*pt_/bet0 + one) / (one/bet0 + pt_)**2;
                    f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
                    track(2,i) = track(2,i) * f_damp_t;
                    track(4,i) = track(4,i) * f_damp_t;
                    track(6,i) = track(6,i) * (one - rfac) - rfac / bet0;
                    !---- Energy loss like for closed orbit.
                 else
                    !---- Store energy loss on closed orbit.
                    rfac = const * (dxt(1)**2 + dyt(1)**2)
                    pt_ = track(6,1)
                    beta_sqr = (pt_*pt_ + two*pt_/bet0 + one) / (one/bet0 + pt_)**2;
                    f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
                    TRACK(2,i) = TRACK(2,i) * f_damp_t;
                    TRACK(4,i) = TRACK(4,i) * f_damp_t;
                    TRACK(6,i) = TRACK(6,i) * (one - rfac) - rfac / bet0;
                 endif
              endif
           else !   step.eq.2, body of the solenoid
              !       Ripken formulae p.29 (3.37)
              track(1,i) =  xf  * cosTh  +  yf  * sinTh
              track(2,i) =  pxf_
              track(3,i) = -xf  * sinTh  +  yf  * cosTh
              track(4,i) =  pyf_
              track(5,i) =  (sigf + (xf*pyf - yf*pxf)*Z) / bet0
              ! track(6,i) =  psigf*bet0
           endif
        enddo ! step
     enddo ! i

  else
     if (sk.ne.zero) then
        skl = sk*length

        !---- Loop over particles
        do  i = 1, ktrack
           do step = 1, 3
              ! initial phase space coordinates
              x_  = track(1,i)
              y_  = track(3,i)
              px_ = track(2,i)
              py_ = track(4,i)
              z_  = track(5,i)
              pt_ = track(6,i)

              ! set up constants
              onedp = sqrt(one + two*pt_/bet0 + pt_**2);

              ! set up constants
              cosTh = cos(two*skl/onedp)
              sinTh = sin(two*skl/onedp)
              omega = sk/onedp;

              ! Store the kick for radiation calculations
              pxf_ = (omega*((cosTh-one)*y_-sinTh*x_)+py_*sinTh+px_*(one+cosTh))/two;
              pyf_ = (omega*((one-cosTh)*x_-sinTh*y_)-px_*sinTh+py_*(one+cosTh))/two;
              dxt(i) = pxf_ - track(2,i);
              dyt(i) = pyf_ - track(4,i);

              if ((step.eq.1).or.(step.eq.3)) then
                 if (radiate .and. elrad .gt. zero) then
                    !---- Full damping.
                    if (damp) then
                       curv = sqrt(dxt(i)**2 + dyt(i)**2) / length;
                       if (quantum) then
                          call trphot(length,curv,rfac,track(6,i))
                       else
                          const = arad * (bet0 * gammas)**3 / three
                          rfac = const * curv**2 * length
                       endif
                       beta_sqr = (pt_*pt_ + two*pt_/bet0 + one) / (one/bet0 + pt_)**2;
                       f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
                       track(2,i) = track(2,i) * f_damp_t;
                       track(4,i) = track(4,i) * f_damp_t;
                       track(6,i) = track(6,i) * (one - rfac) - rfac / bet0;
                       !---- Energy loss like for closed orbit.
                    else
                       !---- Store energy loss on closed orbit.
                       rfac = const * (dxt(1)**2 + dyt(1)**2)
                       beta_sqr = (pt_*pt_ + two*pt_/bet0 + one) / (one/bet0 + pt_)**2;
                       f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
                       track(2,i) = track(2,i) * f_damp_t;
                       track(4,i) = track(4,i) * f_damp_t;
                       track(6,i) = track(6,i) * (one - rfac) - rfac / bet0;
                    endif
                 endif
              else
                 ! total path length traveled by the particle
                 bet = onedp / (one/bet0 + pt_);
                 length_ = length - half/(onedp**2)*(omega*(sinTh-two*length*omega)*(x_**2+y_**2)+&
                      two*(one-cosTh)*(px_*x_+py_*y_)-(sinTh/omega+two*length)*(px_**2+py_**2))/four;

                 ! Thick transport
                 track(1,i) = ((one+cosTh)*x_+sinTh*y_+(px_*sinTh-py_*(cosTh-one))/omega)/two;
                 track(3,i) = ((one+cosTh)*y_-sinTh*x_+(py_*sinTh+px_*(cosTh-one))/omega)/two;
                 track(2,i) = pxf_;
                 track(4,i) = pyf_;
                 track(5,i) = z_ + length/bet0 - length_/bet;
              endif
           enddo ! step

        enddo ! i
     else
        call ttdrf(length,track,ktrack);
     endif
  endif
end subroutine trsol

subroutine tttrans(track,ktrack)
  use trackfi, only : betas
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  !   Translation of a set of particles.                                 *
  ! Input/output:                                                        *
  !   TRACK(6,*)(double)    Track coordinates: (X, PX, Y, PY, T, PT).    *
  !----------------------------------------------------------------------*
  
  double precision :: track(6,*)
  integer :: ktrack

  integer :: i
  double precision :: t_x, t_y, t_z
  double precision :: node_value


  !---- Get translation parameters
  t_x    = node_value('dx ')
  t_y    = node_value('dy ')
  t_z    = node_value('ds ')

  !---- Loop over particles
  
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
  do  i = 1, ktrack
     ! Add vector to particle coordinates
     track(1,i) = track(1,i) - t_x
     track(3,i) = track(3,i) - t_y
     track(5,i) = track(5,i) - t_z/betas
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
  double precision :: ek(6), re(36), track(6,*)
  integer :: ktrack

  integer :: i, j
  double precision :: temp(6)

  do i = 1, ktrack
     do j = 1, 6
        temp(j) = ek(j)  &
             + re(j)    * track(1,i) + re(j+ 6) * track(2,i)  &
             + re(j+12) * track(3,i) + re(j+18) * track(4,i)  &
             + re(j+24) * track(5,i) + re(j+30) * track(6,i)
     enddo
     TRACK(1:6,i) = TEMP(1:6)
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
  integer       :: turn
  character(len=25) :: cmd1
  character(len=30) :: cmd2

  !---- call pro_input('TR$TURN := turn;')
  write(cmd1, '(''tr$turni := '',i8,'' ; '')') turn
  call pro_input(cmd1)
  write(cmd2, '(''exec, tr$macro($tr$turni) ; '')')
  call pro_input(cmd2)
  call init_elements() ! added since now temporary variables are used and need to update
end subroutine trupdate

subroutine trclor(switch,orbit0)
  use twiss0fi
  use name_lenfi
  use trackfi
  use matrices, only : EYE
  use math_constfi, only : zero, one
  use code_constfi
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
  integer :: switch
  double precision :: orbit0(6)

  logical :: aperflag, onepass, debug, thin_foc
  integer :: itra

  integer :: i, j, k, bbd_pos, j_tot, code, irank, n_align
  integer :: pmax, turn, turnmax, part_id(1), last_turn(1)
  integer :: nint, ndble, nchar, int_arr(1), char_l
  double precision :: re(6,6)
  double precision :: z(6,7), zz(6), z0(6,7), z00(6,7), a(6,7), ddd(6)
  double precision :: cotol, err, deltap, el, dxt(200), dyt(200)
  double precision :: al_errors(align_max)
  double precision :: sum, orbit(6), theta
  double precision :: last_pos(6), last_orbit(6,1), maxaper(6)
  character(len=12) :: char_a

  integer, parameter :: itmax=10

  integer, external :: restart_sequ, advance_node, get_option, node_al_errors
  double precision, external :: node_value, get_value, get_variable, get_length

  write (*,'(/a)')          'Full 6D closed orbit search.'
  write (*,'(a)')           'Initial value of 6-D closed orbit from Twiss: '
  write (*,'(a,1p,6e14.6)') 'orbit0 ', ORBIT0

  !---- Initialize variables
  turn    = 1
  turnmax = 1
  pmax    = 7
  sum     = zero
  aperflag  = .false.
  onepass   = .true.

  do k = 1, 7
     Z(:,k) = ORBIT0
  enddo

  DDD(1:6) = 1d-15

! How does it work without the code right after? i.e. A will always be singular!
!  do k = 1, 6
!     z(k,k+1) = z(k,k+1) + ddd(k)
!  enddo

  Z0  = Z
  Z00 = Z

  !--- jmax may be reduced by particle loss - keep number in j_tot
  j_tot = pmax

  !--- get vector of six coordinate maxapers (both RUN and DYNAP)
  call comm_para('maxaper ', nint, ndble, nchar, int_arr, maxaper, char_a, char_l)

  cotol = get_variable('twiss_tol ')

  !---- Initialize kinematics and orbit
  bet0    = get_value('beam ','beta ')
  betas   = get_value('probe ','beta ')
  gammas  = get_value('probe ','gamma ')
  bet0i   = one / bet0
  beti    = one / betas
  dtbyds  = get_value('probe ','dtbyds ')
  deltas  = get_variable('track_deltap ')
  deltap  = get_value('probe ','deltap ')
  arad    = get_value('probe ','arad ')
  radiate = get_value('probe ','radiate ') .ne. zero
  damp    = get_option('damp ') .ne. 0
  quantum = get_option('quantum ') .ne. 0
  debug   = get_option('debug ') .ne. 0 .and. get_option('trace ') .ne. 0

  ORBIT = ORBIT0

  !---- Iteration for closed orbit.
  debug = get_option('debug ') .ne. 0
  thin_foc = get_option('thin_foc ') .eq. 1
  do itra = 1, itmax
     j = restart_sequ()

     do !---- loop over nodes

        bbd_pos = j
        code    = node_value('mad8_type ')
!        if (code .eq. code_tkicker)     code = code_kicker
        if (code .eq. code_placeholder) code = code_instrument
        el      = node_value('l ')
        
        
      !  if (itra .eq. 1 .and. &
      !      code.ne.code_drift .and. &
      !      code.ne.code_quadrupole .and. &
      !      code.ne.code_rbend .and. &
      !      code.ne.code_sbend .and. &
      !      code.ne.code_matrix .and. el.ne.zero ) then
        if (itra .eq. 1 .and. &
            (code.eq.code_sextupole .or. &
            code.eq.code_octupole .or. &
            code.eq.code_elseparator .or. &
            code.eq.code_crabcavity) .and. el.ne.zero) then
           !  .not.(is_drift() .or. is_thin() .or. is_quad() .or. is_dipole() .or. is_matrix()) ) then
           print *,"\ncode: ",code," el: ",el,"   THICK ELEMENT FOUND\n"
           print *,"Track dies nicely"
           print *,"Thick lenses will get nowhere"
           print *,"MAKETHIN will save you\n\n"
           print *, code_sextupole, code_octupole
           call fort_fail('TRRUN: Fatal ','----element with length found : CONVERT STRUCTURE WITH MAKETHIN')
        endif

        !--------  Misalignment at beginning of element (from twissfs.f)
        if (code .ne. code_drift)  then
           AL_ERRORS(:align_max) = zero
           n_align = node_al_errors(al_errors)
           if (n_align .ne. 0)  then
              do i = 1, pmax
                 ZZ = Z(:,i)
                 call tmali1(zz, al_errors, betas, gammas, z(1,i), re)
              enddo
           endif
        endif
    
    theta = node_value('tilt ')
        !-------- Track through element
        call ttmap(switch,code,el,z,pmax,dxt,dyt,sum,turn,part_id, &
             last_turn,last_pos,last_orbit,aperflag,maxaper,al_errors,onepass, debug, theta, thin_foc)

        !--------  Misalignment at end of element (from twissfs.f)
        if (code .ne. code_drift .and. n_align .ne. 0)  then
           do i = 1, pmax
              ZZ = Z(:,i)
              call tmali2(el, zz, al_errors, betas, gammas, z(1,i), re)
           enddo
        endif

        if (advance_node() .eq. 0) exit

        j=j+1
     enddo !---- end of loop over nodes

     !---- construct one-turn map
     do k=1,6
        A(:,k) = ( Z(:,k+1) - Z(:,1) ) / DDD(:)
     enddo

     !---- Solve for dynamic case.
     A(:6,:6) = A(:6,:6) - EYE
     A(:6,7)  = Z(:6,1) - Z0(:6,1)
     err = maxval(abs(A(:,7)))

     call solver(a,6,1,irank)
     if (irank .lt. 6) goto 100

     Z0(:,1)  = Z0(:,1)  - A(:,7)
     Z00(:,1) = Z00(:,1) - A(:,7)

     do k = 2, 7
        Z0(:,k) = Z00(:,1)
     enddo

     do k = 1, 6
        z0(k,k+1) = z0(k,k+1) + ddd(k)
     enddo

    if (debug) then
      write (*,'(a,42e14.6)') 'Z= ', Z
      write (*,'(a,42e14.6)') 'A= ', A
      write (*,'(a,42e14.6)') 'Z0= ', Z0
    endif

     Z = Z0
     !---- end of Iteration
  enddo

  ORBIT0(1:6) = Z0(1:6,1)
  goto 110

100 continue
  write (*,'(a)') '  Singular matrix occurred during closed orbit search.'

110 continue
  write (*,'(/a)')                          '6D closed orbit found by subroutine trclor '
  write (*,'(a,i3,a,1p,e14.6,a,1p,e14.6)')  'iteration: ',itra,' error: ',err,' deltap: ',deltap
  write (*,'(a,1p,6e14.6)')                 'orbit: ', ORBIT0

end subroutine trclor

subroutine ttnllens(track,ktrack)
  use math_constfi, only : zero, one, two, half, pi
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
  double precision :: track(6,*)
  integer :: ktrack

  integer :: i
  double precision :: knll, cnll, dd, u, v, dUu, dUv, dux, duy, dvx, dvy, x, y
  double precision, external :: node_value

  cnll = node_value('cnll ')
  knll = node_value('knll ') / cnll

  do i = 1, ktrack

    x = track(1,i) / cnll
    y = track(3,i) / cnll

    u = half*sqrt((x-one)**2+y**2) + half*sqrt((x+one)**2+y**2)
    v = half*sqrt((x+one)**2+y**2) - half*sqrt((x-one)**2+y**2)
    if (u .eq. one) then
       dd = zero
    else
       dd = u**2 * log(u + sqrt(u*u-one))/sqrt(u**2-one)
    endif
    dUu = (u + log(u+sqrt(u*u-one))*sqrt(u**2-one) + dd)/(u**2-v**2) &
         - two*u*(u*log(u+sqrt(u*u-one))*sqrt(u**2-one) &
                + v*(acos(v)-pi/two)*sqrt(one-v**2)) /(u**2-v**2)**2
    dUv = two*v*(u*log(u+sqrt(u*u-one))*sqrt(u**2-one) &
               + v*(acos(v)-pi/two)*sqrt(one-v**2)) /(u**2-v**2)**2 &
         - (v - (acos(v)-pi/two)*sqrt(one-v**2) + v**2*(acos(v)-pi/two)/sqrt(one-v**2)) &
            / (u**2-v**2)
    dux = half*(x-one)/sqrt((x-one)**2+y**2) + half*(x+one)/sqrt((x+one)**2+y**2)
    duy = half*y/sqrt((x-one)**2+y**2) + half*y/sqrt((x+one)**2+y**2)
    dvx = half*(x+one)/sqrt((x+one)**2+y**2) - half*(x-one)/sqrt((x-one)**2+y**2)
    dvy = half*y/sqrt((x+one)**2+y**2) - half*y/sqrt((x-one)**2+y**2)

    track(2,i) = track(2,i) + knll*(dUu*dux+dUv*dvx)
    track(4,i) = track(4,i) + knll*(dUu*duy+dUv*dvy)

  enddo
end subroutine ttnllens

!FIXME Unused dummy argument 'turn'
subroutine ttrfmult(track, ktrack, turn)
  use twtrrfi
  use trackfi
  use math_constfi, only : zero, one, two, three, ten3m, ten6p, twopi
  use phys_constfi, only : clight
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
  ! BV-flag is applied as: P o inverse(M) o P

  double precision :: track(6,*)
  integer :: ktrack, turn

  integer :: n_ferr, jtrk, iord, nord
  !--- AL: RF-multipole
  integer :: dummyi, nn, ns
  double precision :: normal(0:maxmul), skew(0:maxmul)
  double precision :: f_errors(0:maxferr)
  double precision :: field(2,0:maxmul)
  double precision :: field_cos(2,0:maxmul)
  double precision :: field_sin(2,0:maxmul)
  double precision :: pnl(0:maxmul), psl(0:maxmul)
  double precision :: pc, krf, rfac, curv
  double precision :: x, y, z, px, py, pt, dpx, dpy, dpt
  double precision :: freq, volt, lag, harmon
  double precision :: beta, bvk, deltap, elrad
  double precision :: beta_sqr, f_damp_t

  ! LD: 2016.03.04 - use complex declaration compatible with Lahey
  integer, parameter:: dp=kind(0.d0)
  complex(kind=dp) :: Cm2, Sm2, Cm1, Sm1, Cp0, Sp0, Cp1, Sp1

  double precision, external :: node_value, get_value
  integer, external :: node_fd_errors
  complex(kind=dp), parameter :: ii=(0.0_dp,1.0_dp)

  !---- Zero the arrays
  NORMAL = zero
  SKEW = zero
  PNL = zero
  PSL = zero
  F_ERRORS = zero
  FIELD = zero

  !---- Read-in the parameters
  freq = node_value('freq ');
  lag = node_value('lag ');
  harmon = node_value('harmon ');
  bvk = node_value('other_bv ')
  elrad = node_value('lrad ')
  deltap = get_value('probe ', 'deltap ')
  radiate = get_value('probe ','radiate ') .ne. zero
  arad = get_value('probe ','arad ')
  gammas = get_value('probe ','gamma ')
  beta = get_value('probe ','beta ')
  pc = get_value('probe ','pc ')

  n_ferr = node_fd_errors(f_errors)
  call get_node_vector('knl ', nn, normal)
  call get_node_vector('ksl ', ns, skew)
  call get_node_vector('pnl ', dummyi, pnl)
  call get_node_vector('psl ', dummyi, psl)

  rfac = zero

  !---- Set-up some parameters
  volt = bvk * node_value('volt ')
  krf = twopi * freq * ten6p/clight

  if (n_ferr .gt. 0) call dcopy(f_errors,field,n_ferr)

  nord = max(nn, ns, n_ferr/2-1)

  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1,ktrack
    ! apply the transformation P: (-1, 1, 1, -1, -1, 1) * X
    x  = bvk * track(1,jtrk)
    px =       track(2,jtrk)
    y  =       track(3,jtrk)
    py = bvk * track(4,jtrk)
    z  = bvk * track(5,jtrk)
    pt =       track(6,jtrk)

    !---- Vector with strengths + field errors
    do iord = 0, nord
      field_cos(1,iord) = bvk * (normal(iord) * cos(pnl(iord) * twopi - krf * z) + field(1,iord))
      field_sin(1,iord) = bvk * (normal(iord) * sin(pnl(iord) * twopi - krf * z))
      field_cos(2,iord) = bvk * (skew(iord)   * cos(psl(iord) * twopi - krf * z) + field(2,iord))
      field_sin(2,iord) = bvk * (skew(iord)   * sin(psl(iord) * twopi - krf * z))
    enddo
    Cm2 = zero; Sm2 = zero; Cm1 = zero; Sm1 = zero;
    Cp0 = zero; Sp0 = zero; Cp1 = zero; Sp1 = zero;

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
    ! apply the transformation P: (-1, 1, 1, -1, -1, 1) * X
    dpx = -REAL(Cp0);
    dpy = AIMAG(Cp0);
    dpt = (volt * ten3m * sin(lag * twopi - krf * z) / pc - krf * REAL(Sp1));

    !---- Radiation effects at entrance.
    if (radiate  .and.  elrad .ne. zero) then
       if (quantum) then
          curv = sqrt(dpx**2+dpy**2) / elrad;
          call trphot(elrad,curv,rfac,pt)
       else
          rfac = arad * (betas * gammas)**3 * (dpx**2+dpy**2) / (three*elrad)
       endif
       if (damp) then
          beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
          f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
          px = px * f_damp_t;
          py = py * f_damp_t;
          pt = pt * (one - rfac) - rfac / bet0;
       else
          if (jtrk.eq.1) then
             beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
             f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
          endif
          px = px * f_damp_t;
          py = py * f_damp_t;
          pt = pt * (one - rfac) - rfac / bet0;
       endif
    endif

    !---- Apply the kick
    px = px + dpx
    py = py + dpy
    pt = pt + dpt

    !---- Radiation effects at exit.
    if (radiate  .and.  elrad .ne. zero) then
       if (quantum) then
          curv = sqrt(dpx**2+dpy**2) / elrad;
          call trphot(elrad,curv,rfac,pt)
       else
          rfac = arad * (betas * gammas)**3 * (dpx**2+dpy**2) / (three*elrad)
       endif
       if (damp) then
          beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
          f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
          px = px * f_damp_t;
          py = py * f_damp_t;
          pt = pt * (one - rfac) - rfac / bet0;
       else
          if (jtrk.eq.1) then
             beta_sqr = (pt*pt + two*pt/bet0 + one) / (one/bet0 + pt)**2;
             f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
          endif
          px = px * f_damp_t;
          py = py * f_damp_t;
          pt = pt * (one - rfac) - rfac / bet0;
       endif
    endif

    ! apply the transformation P: (-1, 1, 1, -1, -1, 1) * X
    track(1,jtrk) = bvk * x;
    track(2,jtrk) =       px;
    track(3,jtrk) =       y;
    track(4,jtrk) = bvk * py;
    track(5,jtrk) = bvk * z;
    track(6,jtrk) =       pt;

  enddo

end subroutine ttrfmult

subroutine ttcfd(x, px, y, py, z, pt, h, k0_, k1_, length)

  use trackfi
  use math_constfi, only : zero, half, one, two, three, six
  implicit none

  !-------------------------*
  ! Andrea Latina 2017      *
  !-------------------------*
  !----------------------------------------------------------------------------*
  ! Purpose:                                                                   *
  !    Track a particle through a combined function dipole.                    *
  ! Input:                                                                     *
  !   h (double)                     curvature of the reference orbit, 1/rho0  *
  !   k0 (double)                    bending strength (in a sbend k0 == h) 1/m *
  !   k1 (double)                    focusing strength 1/m^2                   *
  ! Input/output:                                                              *
  !   (X, PX, Y, PY, T, PT)(double)  Track coordinates.                        *
  !----------------------------------------------------------------------------*

  double precision :: x_, px_, y_, py_, z_, length_
  double precision :: x, px, y, py, z, pt, xp, yp
  double precision :: h, k0_, k1_, length, Kx, Ky
  double precision :: A, B, C, D, k0, k1
  double precision :: delta_plus_1, bet
  double precision :: Cx, Sx, Cy, Sy
  integer, parameter:: dp=kind(0.d0)
  complex(kind=dp) :: sqrt_Kx, sqrt_Ky

  delta_plus_1 = sqrt(pt*pt + two*pt/bet0 + one);
  bet = delta_plus_1/(one/bet0+pt);

  k0 = k0_ / delta_plus_1; ! 1/m
  k1 = k1_ / delta_plus_1; ! 1/m^2
  Kx = k0*h + k1; ! 1/m^2
  Ky =       -k1; ! 1/m^2
  if (Kx.ne.zero) then
     sqrt_Kx = cdsqrt(DCMPLX(Kx)); ! 1/m
     Sx = dreal(cdsin(sqrt_Kx*length) / sqrt_Kx); ! m
     Cx = dreal(cdcos(sqrt_Kx*length)); ! 1
  else
     sqrt_Kx = zero; ! 1/m
     Sx = length; ! m
     Cx = one; ! 1
  endif
  if (Ky.ne.zero) then ;
     sqrt_Ky = cdsqrt(DCMPLX(Ky)); ! 1/m
     Sy = dreal(cdsin(sqrt_Ky*length) / sqrt_Ky); ! m
     Cy = dreal(cdcos(sqrt_Ky*length)); ! 1
  else
     sqrt_Ky = zero; ! 1/m
     Sy = length; ! m
     Cy = one; ! 1
  endif

  xp = px/delta_plus_1;
  yp = py/delta_plus_1;

  ! useful constants
  A = -Kx*x-k0+h; ! 1/m
  B = xp;
  C = -Ky*y; ! 1/m
  D = yp;

  ! transverse map
  x_ = x*Cx + xp*Sx;
  y_ = y*Cy + yp*Sy;
  px_ = (A*Sx + B*Cx) * delta_plus_1;
  py_ = (C*Sy + D*Cy) * delta_plus_1;

  if (Kx.ne.zero) then
     x_ = x_ + (k0-h)*(Cx-one)/Kx;
  else
     x_ = x_ - (k0-h)*half*length**2;
  endif

  ! longitudinal map
  length_ = length; ! will be the total path length traveled by the particle
  if (Kx.ne.zero) then
     length_ = length_ - (h*((Cx-one)*xp+Sx*A+length*(k0-h)))/Kx;
     length_ = length_ + half*(-(A**2*Cx*Sx)/(two*Kx)+(B**2*Cx*Sx)/two+&
          (A**2*length)/(two*Kx)+(B**2*length)/two-(A*B*Cx**2)/Kx+(A*B)/Kx);
  else
     length_ = length_ + h*length*(three*length*xp+six*x-(k0-h)*length**2)/six;
     length_ = length_ + half*(B**2)*length;
  endif
  if (Ky.ne.zero) then
     length_ = length_ + half*(-(C**2*Cy*Sy)/(two*Ky)+(D**2*Cy*Sy)/two+&
          (C**2*length)/(two*Ky)+(D**2*length)/two-(C*D*Cy**2)/Ky+(C*D)/Ky);
  else
     length_ = length_ + half*(D**2)*length;
  endif
  z_ = z + length/bet0 - length_/bet;

  x  = x_;
  px = px_;
  y  = y_;
  py = py_;
  z  = z_;

end subroutine ttcfd

subroutine tttquad(track, ktrack)
  use twtrrfi
  use trackfi
  use twiss_elpfi
  use math_constfi, only : zero, one, two, three, half
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
  double precision :: track(6,*)
  integer :: ktrack

  double precision :: k1, k1s, length
  double precision :: x, px, y, py, z, pt
  double precision :: delta_plus_1
  double precision :: ct, st
  double precision :: tilt
  double precision :: tmp

  double precision :: hx, hy, rfac, gamma, beta, curv
  double precision :: beta_gamma, beta_sqr, f_damp_t
  integer :: jtrk, elpar_vl

  double precision, external :: node_value
  double precision, parameter ::  sqrt2=1.41421356237310d0

  double precision, external :: get_value

  double precision :: f_errors(0:maxferr)
  integer, external :: node_fd_errors, el_par_vector
  integer :: n_ferr
  
  !gamma = get_value('probe ','gamma ')
  !beta = get_value('probe ','beta ')

  !---- Read-in the parameters
  elpar_vl = el_par_vector(q_k1st, g_elpar)
  
  length = node_value('l ');
  tilt = g_elpar(q_tilt)

  f_errors = zero
  n_ferr = node_fd_errors(f_errors)
  k1  = g_elpar(q_k1)  + g_elpar(q_k1t)
  k1s = g_elpar(q_k1s) + g_elpar(q_k1st)
  
  !k1  = node_value('k1 ')
  !k1s = node_value('k1s ')


  if (length.ne.zero) then
     k1  = k1  + f_errors(2)/length
     k1s = k1s + f_errors(3)/length
  endif

  if (k1s.ne.zero) then
     tilt = -atan2(k1s, k1)/two ! + tilt
     k1 = sqrt(k1**2 + k1s**2)
  else
     tilt = zero
  endif

  if (k1.eq.zero) then
     call ttdrf(length,track,ktrack);
     return
  endif

  if (tilt.ne.zero)  then
     st = sin(tilt)
     ct = cos(tilt)
  endif

  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1, ktrack
     !---- The particle position
     x  = track(1,jtrk);
     px = track(2,jtrk);
     y  = track(3,jtrk);
     py = track(4,jtrk);
     z  = track(5,jtrk);
     pt = track(6,jtrk);

     !---  rotate orbit before entry
     if (tilt .ne. zero)  then
        tmp = x
        x = ct * tmp + st * y
        y = ct * y   - st * tmp
        tmp = px
        px = ct * tmp + st * py
        py = ct * py  - st * tmp
     endif

     !---- Computes 1+delta
     delta_plus_1 = sqrt(pt*pt + two*pt/beta + one);

     !---- Radiation effects at entrance
     if (radiate) then
        hx = (-k1*x) / delta_plus_1;
        hy = ( k1*y) / delta_plus_1;
        if (quantum) then
           curv = sqrt(hx**2+hy**2);
           call trphot(length,curv,rfac,pt)
        else
           beta_gamma = delta_plus_1 * gamma * beta;
           rfac = (arad * beta_gamma**3 * length / three) * (hx**2 + hy**2);
        endif
        if (damp) then
           beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        else
           if (jtrk.eq.1) then
              beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
              f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           endif
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        endif
     endif

     call ttcfd(x, px, y, py, z, pt, 0d0, 0d0, k1, length);

     !---- Radiation effects at exit
     if (radiate) then
        hx = (-k1*x) / delta_plus_1;
        hy = ( k1*y) / delta_plus_1;
        if (quantum) then
           curv = sqrt(hx**2+hy**2);
           call trphot(length,curv,rfac,pt)
        else
           beta_gamma = delta_plus_1 * gamma * beta;
           rfac = (arad * beta_gamma**3 * length / three) * (hx**2 + hy**2);
        endif
        if (damp) then
           beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        else
           if (jtrk.eq.1) then
              beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
              f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           endif
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        endif
     endif

     !---  rotate orbit at exit
     if (tilt .ne. zero)  then
        tmp = x
        x = ct * tmp - st * y
        y = ct * y   + st * tmp
        tmp = px
        px = ct * tmp - st * py
        py = ct * py  + st * tmp
     endif

     !---- Applies the kick
     track(1,jtrk) = x
     track(2,jtrk) = px
     track(3,jtrk) = y
     track(4,jtrk) = py
     track(5,jtrk) = z
     track(6,jtrk) = pt

  enddo

end subroutine tttquad

subroutine tttdipole(track, ktrack, code)
  use twtrrfi
  use trackfi
  use math_constfi, only : zero, one, two, three, half
  use code_constfi
  use twiss_elpfi
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

  double precision :: track(6,*)
  integer :: ktrack

  integer :: jtrk
  double precision :: length, angle, rho, h, k0, k1
  double precision :: x, px, y, py, z, pt, delta_plus_1
  double precision :: beta, gamma, hx, hy, rfac, curv
  double precision :: e1, e2, h1, h2, hgap, fint, fintx
  double precision :: beta_sqr, beta_gamma, f_damp_t

  double precision, external :: node_value, get_value
  double precision :: f_errors(0:maxferr)
  integer, external :: node_fd_errors, el_par_vector
  integer :: n_ferr, code, elpar_vl

  !code    = node_value('mad8_type ')
  !arad    = get_value('probe ','arad ')
  !beta    = get_value('probe ','beta ')
  !gamma   = get_value('probe ','gamma ')
  !radiate = get_value('probe ','radiate ') .ne. zero
  !All these were removed since they were global parameters. 
  
  elpar_vl = el_par_vector(b_k3s, g_elpar)
  !---- Read-in dipole edges angles
  !e1    = node_value('e1 ');
  !e2    = node_value('e2 ');
  e1 = g_elpar(b_e1)
  e2 = g_elpar(b_e2)
  !h1    = node_value('h1 ')
  !h2    = node_value('h2 ')
  h1 = g_elpar(b_h1)
  h2 = g_elpar(b_h2)
  !hgap  = node_value('hgap ')
  !fint  = node_value('fint ')
  !fintx = node_value('fintx ')
  hgap = g_elpar(b_hgap)
  fint = g_elpar(b_fint)
  fintx = g_elpar(b_fintx)
  
  length = node_value('l ')
  angle  = g_elpar(b_angle)

  rho = abs(length/angle)
  h = angle/length
  k0 = g_elpar(b_k0)
  k1 = g_elpar(b_k1)
  !k0 = node_value('k0 ') ! was h
  !k1 = node_value('k1 ')


  if (code .eq. code_rbend) then
     e1 = e1 + angle / two;
     e2 = e2 + angle / two;
  endif

  !---- Apply errors
  f_errors = zero
  n_ferr = node_fd_errors(f_errors)
  if (k0.ne.0) then
    f_errors(0) = f_errors(0) + k0*length - angle
  else
    k0 = h
  endif
  k0 = k0 + f_errors(0) / length ! dipole term
  k1 = k1 + f_errors(2) / length ! quad term

  if (k0.eq.zero .and. k1.eq.zero) then
     call ttdrf(length,track,ktrack);
     return
  endif
  !---- Apply entrance dipole edge effect
  if (node_value('kill_ent_fringe ') .eq. zero) &
       call ttdpdg_map(track, ktrack, e1, k0, hgap, fint, zero)

  !---- Prepare to calculate the kick and the matrix elements
  do jtrk = 1,ktrack
     !---- The particle position
     x  = track(1,jtrk);
     px = track(2,jtrk);
     y  = track(3,jtrk);
     py = track(4,jtrk);
     z  = track(5,jtrk);
     pt = track(6,jtrk);

     delta_plus_1 = sqrt(pt*pt + two*pt/beta + one);

     !---- Radiation effects at entrance.
     if (radiate) then
        !hx = (-k0 -k1*x + k1s*y  - k2*(x*x - y*y)/two + k2s*x*y) / delta_plus_1; ! if there were k1s k2 and k2s
        !hy = (     k1*y + k1s*x  + k2*x*y + k2s*(x*x - y*y)/two) / delta_plus_1;
        hx = (-k0 -k1*x) / delta_plus_1;
        hy = (     k1*y) / delta_plus_1;
        if (quantum) then
           curv = sqrt(hx**2+hy**2);
           call trphot(length * (one + h*x) - two * tan(e1)*x, curv, rfac, pt);
        else
           beta_gamma = delta_plus_1 * gammas * beta;
           rfac = (arad * beta_gamma**3 * two / three) * (hx**2 + hy**2) * (length / two * (one + h*x) - tan(e1)*x)
        endif
        if (damp) then
           beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        else
           if (jtrk.eq.1) then
              beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
              f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           endif
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        endif
     endif

     call ttcfd(x, px, y, py, z, pt, h, k0, k1, length);

     !---- Radiation effects at exit.
     if (radiate) then
        !hx = (-k0 -k1*x + k1s*y  - k2*(x*x - y*y)/two + k2s*x*y) / delta_plus_1; ! if there were k1s k2 and k2s
        !hy = (     k1*y + k1s*x  + k2*x*y + k2s*(x*x - y*y)/two) / delta_plus_1;
        hx = (-k0 -k1*x) / delta_plus_1;
        hy = (     k1*y) / delta_plus_1;
        if (quantum) then
           curv = sqrt(hx**2+hy**2);
           call trphot(length * (one + h*x) - two * tan(e2)*x, curv, rfac, pt);
        else
           beta_gamma = delta_plus_1 * gammas * beta;
           rfac = (arad * beta_gamma**3 * two / three) * (hx**2 + hy**2) * (length / two * (one + h*x) - tan(e2)*x)
        endif
        if (damp) then
           beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
           f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        else
           if (jtrk.eq.1) then
              beta_sqr = (pt*pt + two*pt/beta + one) / (one/beta + pt)**2;
              f_damp_t = sqrt(one + rfac*(rfac - two) / beta_sqr);
           endif
           px = px * f_damp_t;
           py = py * f_damp_t;
           pt = pt * (one - rfac) - rfac / beta;
        endif
     endif

     !---- Applies the kick
     track(1,jtrk) = x
     track(2,jtrk) = px
     track(3,jtrk) = y
     track(4,jtrk) = py
     track(5,jtrk) = z
     track(6,jtrk) = pt

  enddo

  !---- Apply exit dipole edge effect
  if (node_value('kill_exi_fringe ') .eq. zero) then
     if (fintx .lt. zero) fintx = fint
     call ttdpdg_map(track, ktrack, e2, k0, hgap, fintx, zero)
  endif

end subroutine tttdipole

subroutine trphot(el,curv,rfac,pt)
  use math_constfi, only : zero, one, two, three, five, twelve
  use phys_constfi, only : clight, hbar
  use trackfi, only : arad, gammas, betas
  implicit none
  !----------------------------------------------------------------------*
  ! Purpose:                                                             *
  ! option synrad=1                                                      *
  !   Generate random energy loss for photons, using a look-up table to  *
  !   invert the function Y.  Ultra-basic interpolation computed;        *
  !   leads to an extrapolation outside the table using the two outmost  *
  !   point on each side (low and high).                                 *
  !   Assumes ultra-relativistic particles (beta = 1).                   *
  !   Author: Ghislain Roy                                               *
  ! option synrad=2                                                      *
  !   Generate random energy loss for photons, using the                 *
  !   Synchrotron radiation spectrum generator                           *
  !   described in CERN-OPEN-2007-018 by Helmut Burkhardt                *                                                 *

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
  double precision :: el, curv, rfac, pt

  integer :: i, ierror, j, nphot
  double precision :: amean, dlogr, slope, ucrit, xi, sumxi, InvSynFracInt,rn_arg
  integer, parameter :: maxtab=101
  double precision :: tabxi(maxtab),taby(maxtab)
  double precision :: pc, gamma, amass
  character(len=20) text
  integer, external :: get_option
  integer :: synrad
  double precision :: get_value, frndm, delta_plus_1
  double precision, parameter :: fac1=3.256223d0
  ! integer, parameter :: maxran=1000000000

  data (taby(i), i = 1, 101)                                         &
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
       1.0750246d0,    1.08283985d0,   1.0899564d0,    1.09645379d0,     &
       1.10352755d0,   1.11475027d0,   1.12564385d0,   1.1306442d0,      &
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
  data (tabxi(i), i = 1, 101)                                        &
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
       -0.10536054d0,  -0.05129330d0,   0.0d0,          0.048790119d0,   &
       0.104360029d0,  0.198850885d0,  0.300104618d0,  0.350656837d0,    &
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
  pc     = get_value('probe ','pc ')
  amass  = get_value('probe ','mass ')

  gamma = (betas*pt + one)*gammas;
  delta_plus_1 = sqrt(pt*pt + two*pt/betas + one);

  !---- AMEAN is the average number of photons emitted.,
  !     NPHOT is the integer number generated from Poisson's law.
  !-AL- AMEAN implicitly takes el / 2 (half the element length)
  amean = five * sqrt(three) / (twelve * hbar * clight) * abs(arad * pc * delta_plus_1 * el * curv)
  ucrit = three/two * hbar * clight * gamma**3 * abs(curv)
  sumxi = zero
  if (amean.gt.3d-1) then
     print *,"More than 0.3 photons emitted in element."
     print *,"You might want to consider increasing the number of slices to reduce this number."
  endif

  if (amean .gt. zero) then
     call dpoissn(amean, nphot, ierror)

     if (ierror .ne. 0) then
        write(text, '(1p,d20.12)') amean
        call fort_fail('TRPHOT: ','Fatal: Poisson input mean =' // text)
     endif

     !---- For all photons, sum the radiated photon energy in units of UCRIT
     if (nphot .ne. 0) then
        synrad = get_option('synrad ')
        if (synrad .eq. 1) then
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
             sumxi = sumxi + xi
             ! write(60,*) xi ! 2016-Mar-16  18:22:30  ghislain: dump individual photons to file
          enddo
        else !-- synrad .eq. 2
          do i = 1, nphot
             rn_arg = frndm()
             xi = InvSynFracInt(rn_arg)
  !hbutest   print *,xi,rn_arg," InvSynFracIntdata"
             sumxi = sumxi + xi
          enddo
        endif
     endif
  endif

  ! normalize rfac to beam energy
  rfac = sumxi * ucrit / (gamma*amass)

end subroutine trphot

subroutine dpoissn (amu,n,ierror)
  use math_constfi, only : zero, one, half
  implicit none
  !----------------------------------------------------------------------*
  !    POISSON GENERATOR                                                 *
  !    CODED FROM LOS ALAMOS REPORT      LA-5061-MS                      *
  !    PROB(N)=EXP(-AMU)*AMU**N/FACT(N)                                  *
  !        WHERE FACT(N) STANDS FOR FACTORIAL OF N                       *
  !    ON RETURN IERROR.EQ.0 NORMALLY                                    *
  !              IERROR.EQ.1 IF AMU.LE.0.                                *
  !----------------------------------------------------------------------*
  double precision, intent(IN) :: amu
  integer, intent(OUT) :: n, ierror

  double precision :: ran, pir, expma
  double precision :: frndm, grndm
  ! AMAX IS THE VALUE ABOVE WHICH THE NORMAL DISTRIBUTION MUST BE USED
  double precision, parameter :: amax=88d0

  ierror= 0

  if (amu .le. zero) then
     ! MEAN SHOULD BE POSITIVE
     ierror = 1
     n = 0
     return
  endif

  if (amu .gt. amax) then
     ! NORMAL APPROXIMATION FOR AMU.GT.AMAX
     ran = grndm()
     n = ran * sqrt(amu) + amu + half
     return
  endif

   expma = exp(-amu)
   pir = one
   n = -1

10 n = n + 1
   pir = pir * frndm()
   if (pir .gt. expma) go to 10

end subroutine dpoissn

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
  integer :: i,j,k
  double precision :: wi, wr, x, y
  double precision, parameter :: h=1.d0/63.d0
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
  double precision :: x, y, wr, wi

  integer :: n
  double precision :: rr(37), ri(37), sr0, sr, si, tr, ti, vi, vr, xa, xl, ya, zhi, zhr

  double precision, parameter :: z1=1d0, hf=z1/2d0, z10=10d0
  double precision, parameter :: c1=74d0/z10, c2=83d0/z10, c3=z10/32d0, c4=16d0/z10
  !double precision, parameter :: c=1.12837916709551257d0, p=(2d0*c4)**33
  double precision, parameter :: c=1.12837916709551257d0
  double precision, parameter :: p=46768052394588893.3825d0

  save
  !-----------------------------------------------------------------------
  xa=abs(x)
  ya=abs(y)
  if (ya.lt.c1 .and. xa.lt.c2) then
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
  !hr05 if (ya .eq. 0) then
  if (ya .eq. 0d0) then                                                 !hr05
     !        v=dcmplx(exp(-xa**2),dimag(v))
     !hr05   vr=exp_rn(-xa**2)
     vr=exp(-1d0*xa**2)                                              !hr05
  endif
  if (y.lt.0d0) then
     !        v=2*exp(-dcmplx(xa,ya)**2)-v
     !hr05   vr=2d0*exp_rn(ya**2-xa**2)*cos_rn(2d0*xa*ya)-vr
     vr=(2d0*exp(ya**2-xa**2))*cos((2d0*xa)*ya)-vr              !hr05
     vi=(-2d0*exp(ya**2-xa**2))*sin((2d0*xa)*ya)-vi             !hr05
     !hr05   if (x.gt.0) vi=-vi
     if (x.gt.0d0) vi=-1d0*vi                                          !hr05
  else
     !hr05   if (x.lt.0) vi=-vi
     if (x.lt.0d0) vi=-1d0*vi                                          !hr05
  endif
  wr=vr
  wi=vi
  return
end subroutine mywwerf

subroutine wzsubv(n,vx,vy,vu,vv)
  use math_constfi, only : one, two, half
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
  integer :: n
  double precision :: vx(*), vy(*), vu(*), vv(*)

  integer, parameter :: npart=64, nx=490, ny=470, idim=(nx+2)*(ny+2)

  integer :: i, j, k, vmu, vnu, kstep
  double precision :: vd12i, vd12r, vd23i, vd23r, vd34i, vd34r
  double precision :: vp, vq, vqsq, vr, vsimag, vsreal
  double precision :: vt, vtdd13i, vtdd13r, vtdd24i, vtdd24r, vtdddi, vtdddr, vti, vtr
  double precision :: vusum, vusum3, vvsum, vvsum3
  double precision :: vw1i, vw1r, vw2i, vw2r, vw3i, vw3r, vw4i, vw4r
  double precision :: vxh, vxhrel, vyh, vyhrel
  double precision :: hrecip, wtimag(idim), wtreal(idim)
  double precision :: xx, yy

  double precision, parameter :: xcut=7.77d0, ycut=7.46d0
  double precision, parameter :: h=1.d0/63.d0
  double precision, parameter :: a1=0.5124242248d0, a2=0.0517653588d0
  double precision, parameter :: b1=0.2752551286d0, b2=2.7247448714d0
  double precision, parameter :: xm=1d120
  !     temporary arrays to facilitate vectorisation
  integer :: in, out, ins(npart), outs(npart)

  !common /wzcom1/ hrecip, kstep
  !common /wzcom2/ wtreal(idim), wtimag(idim)

  !-----------------------------------------------------------------------
  integer :: ilo, ihi
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
  !$    if (th .eq. nt-1) ihi=n
  !-----------------------------------------------------------------------
  in=0
  out=0
  do i=ilo,ihi
     if (vx(i).ge.xcut .or. vy(i).ge.ycut) then
        out=out+1
        outs(out)=i
        if (out .eq. npart) then
           !     everything outside the rectangle so approximate
           !     write (*,*) 'ALL outside'
           !     write (*,*) 'i=',i
           do j=1,out
              xx=vx(outs(j))
              yy=vy(outs(j))
              if (xx.ge.xm) xx=xm
              if (yy.ge.xm) yy=xm
              vp=xx**2-yy**2
              vq=(two*xx)*yy
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
        if (in .eq. npart) then
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
              vtdd24i = -one * ( vtr + vti )                             !hr05
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
              vusum3 = half * (vtdd13r + (vxhrel*vtdddr - vyhrel*vtdddi))
              vvsum3 = half * (vtdd13i + (vxhrel*vtdddi + vyhrel*vtdddr))
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
  do j= 1,out
     xx = vx(outs(j))
     yy = vy(outs(j))
     if (xx.ge.xm) xx=xm
     if (yy.ge.xm) yy=xm
     vp = xx**2 - yy**2
     vq = (two*xx)*yy
     vqsq = vq**2
     !  First term.
     vt = vp-b1
     vr = a1/(vt**2+vqsq)
     vsreal =  vr*vt
     vsimag = -vr*vq
     !  Second term
     vt = vp-b2
     vr = a2/(vt**2+vqsq)
     vsreal = vsreal +vr*vt
     vsimag = vsimag -vr*vq
     !  Multiply by i*z.
     vu(outs(j)) = -(yy*vsreal+xx*vsimag)
     vv(outs(j)) =   xx*vsreal-yy*vsimag
  enddo
  !     everything inside the square, so interpolate
  !     write (*,*) 'ALL inside'
  do j= 1,in
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
     vtdd24i = -one * ( vtr + vti )                             !hr05
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
     vusum3 = half * (vtdd13r + (vxhrel*vtdddr-vyhrel*vtdddi))
     vvsum3 = half * (vtdd13i + (vxhrel*vtdddi+vyhrel*vtdddr))
     vyhrel = vyhrel - one
     vusum = vd12r + (vxhrel*vusum3-vyhrel*vvsum3)
     vvsum = vd12i + (vxhrel*vvsum3+vyhrel*vusum3)
     vxhrel = vxhrel - one
     vu(ins(j)) = vw1r + (vxhrel*vusum-vyhrel*vvsum)
     vv(ins(j)) = vw1i + (vxhrel*vvsum+vyhrel*vusum)
  enddo
  !$OMP END PARALLEL
  return
end subroutine wzsubv

subroutine wzsub(x,y,u,v)
  use math_constfi, only : one, two, half
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
  double precision :: x, y, u, v

  integer :: k, mu, nu
  double precision :: d12i, d12r, d23i, d23r, d34i, d34r, p, q, qsq, r
  double precision :: simag, sreal, t, tdd13i, tdd13r, tdd24i, tdd24r, tdddi, tdddr, ti, tr
  double precision :: usum, usum3, vsum, vsum3, w1i, w1r, w2i, w2r, w3i, w3r, w4i, w4r
  double precision :: xh, xhrel, yh, yhrel

  double precision, parameter :: xcut=7.77d0, ycut=7.46d0
  double precision, parameter :: a1=0.5124242248d0, a2=0.0517653588d0
  double precision, parameter :: b1=0.2752551286d0, b2=2.7247448714d0

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
  q=(two*x)*y                                                       !hr05
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
