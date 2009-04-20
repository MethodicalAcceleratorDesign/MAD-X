module madx_ptc_trackline_module
  use madx_ptc_module
  use madx_ptc_intstate_module
  use madx_ptc_setcavs_module
  implicit none
  save
  public

  public                              :: ptc_trackline       ! subroutine inside the module
  public                              :: ptc_track_everystep

  ! flag for debugging ranges from 0 (no debug printout) to 10 (the most detailed)
  real(dp),allocatable :: Dismom(:,:)    ! <xnormal_(2*i-1)**(2j)>= dismon(i,j)*I_i**j

  !********************************************************************************************
  !********************************************************************************************
  !********************************************************************************************

contains

  subroutine ptc_track_everystep(nobs)
    ! subroutine that performs tracking with acceleration
    ! it is called as a result of ptc_trackline MAD-X command

    implicit none
    integer, intent (IN) :: nobs ! the maximum number of observation points >=1
    INTEGER, ALLOCATABLE :: observedelements(:)
    integer  :: charge    ! charge of an accelerated particle
    type(fibre), pointer :: p
    real (dp)            :: x(1:6)
    !    real (dp)            :: polarx(1:6)   ! track vector -
    real (dp)            :: xp, yp, pz, p0
    real (dp)            :: pathlegth = zero
    integer              :: npart = 1
    integer              :: n = 1
    integer              :: nturns = 1
    integer              :: t = 1
    logical(lp)          :: gcs
    logical(lp)          :: rplot
!    real (dp)            :: gposx, gposy, gposz
    integer              :: e, ni
    integer              :: obspointnumber ! observation point number in c-code
    integer              :: getnumberoftracks !function
    type(internal_state)  :: intstate
    real(kind(1d0))      :: get_value
    integer, external    :: restart_sequ, & !  restart beamline and return number of beamline node
         advance_node    !  advance to the next node in expanded sequence
    !                    !  =0 (end of range), =1 (else)
    REAL(KIND(1d0)), external :: node_value  !/*returns value for parameter par of current element */
    TYPE(BEAM) :: TheBEAM
    TYPE(INTEGRATION_NODE),POINTER :: CURR_SLICE,PREV_SLICE
    integer             :: mf

    !------------------------------------------------------
    !initialization
    npart = 1
    n = 1
    t = 1
    !------------------------------------------------------

    if(universe.le.0) then
       call fort_warn('return from ptc_trackline: ',' no universe created')
       print*,"Max number of nobs ", nobs
       return
    endif
    if(index_mad.le.0) then
       call fort_warn('return from ptc_trackline: ',' no layout created')
       return
    endif

    nturns = get_value('ptc_trackline ','turns ')
    if (getdebug() > 2) print *, 'ptc_trackline, nturns = ', nturns

    if ( (nturns > 1) .and. (my_ring%closed .eqv. .false.)) then
       call fort_warn('WARNING: You can not make more than one turn in a line!', &
            'Putting number of turns to 1!')
       nturns = 1
    endif


    gcs = get_value('ptc_trackline ','gcs ') .ne. 0

    rplot = get_value('ptc_trackline ','rootntuple ') .ne. 0

    intstate = getintstate()
    if (gcs .and.  intstate%TOTALPATH==1) then
       call fort_warn("ptc_trackline","Having global coordinates and totalpath for z is sensless")
       gcs = .false.
    endif


    allocate(observedelements(1:my_ring%n)); observedelements(:)=0 ! zero means that this element is not an obs. point

    c_%x_prime=.true.

    e=restart_sequ()
    p=>my_ring%start
    do e=1, my_ring%n

       obspointnumber=node_value('obs_point ')
       IF (e.eq.1) obspointnumber=1 ! node_value gives 0 for 1st (?)

       if (obspointnumber .gt. 0) then
          if (getdebug() > 0) print *,"Element ",e," is an observation point no. ",obspointnumber
          observedelements(e) = obspointnumber
       endif

       obspointnumber=advance_node() ! c-code go to the next node -> the passed value is never used, just to shut up a compiler
       p=>p%next
    enddo


    charge = get_value('beam ', "charge ");
    if (getdebug() > 3 ) print *, 'Read charge:', charge,' layout has charge ', my_ring%start%charge

    if (cavsareset .eqv. .false.) then
       call setcavities(my_ring,maxaccel)
    endif

    if (getdebug() > 0) print *, 'reading tracks starting posiotions from table ....'

    call gettrack(1,x(1),x(2),x(3),x(4),x(6),x(5))

    if (getdebug() > 0) print *, 'reading.... Done'

    if (getdebug() > 0) then
       print *, '###################################################'
       print *, '###################################################'
       print *, '######         TRACKING WITH PTC         ##########'
       print *, '###################################################'
       print *, '###################################################'
    endif

    if (rplot) then
       call newrplot()
    endif


    if(.not.associated(my_ring%t))  then
       CALL MAKE_node_LAYOUT(my_ring)
    endif

    n=1
    npart = getnumberoftracks()
    if (getdebug() > 0) print *, 'There is ', npart,' tracks'

    !     IF(.NOT.ASSOCIATED(TheBeam%N)) THEN
    CALL ALLOCATE_BEAM(TheBeam,npart)
    !     ELSEIF(TheBeam%N/=npart) THEN
    !        CALL KILL_BEAM(TheBeam)
    !        CALL ALLOCATE_BEAM(TheBeam,npart)
    !     ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!    READS DATA FROM MADX         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1, npart

       pathlegth = zero

       !if (getdebug() > 3 )
       print *, 'Getting track ',n

       call gettrack(n,TheBeam%X(n,1),TheBeam%X(n,2),TheBeam%X(n,3),TheBeam%X(n,4),TheBeam%X(n,6),TheBeam%X(n,5))

       !if (getdebug() > 0 )
       write(6,'(a10,1x,i8,1x,6(f9.6,1x))') 'Track ',n,TheBeam%X(n,1:6)

       TheBeam%X(n,7)=ZERO

       if( associated(TheBeam%POS(n)%NODE) ) then
          TheBeam%POS(n)%NODE=>my_ring%start%t1
       endif

    enddo !loop over tracks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!      TRACKING       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call kanalnummer(mf)

    open(unit=mf,file='thintracking_ptc.txt',POSITION='APPEND' , STATUS='UNKNOWN')

    do t=1, nturns
       print*, "TURN NUMBER ",t
       p=>my_ring%start

       PREV_SLICE => my_ring%start%T1
       CURR_SLICE => prev_slice%next
       e = 1

       !       print*,"Name of the first element ", my_ring%start%mag%name
       !       print*,"Position of the first element ", my_ring%start%T2%pos
       !       print*,"Name of the last element ", my_ring%end%mag%name
       !       print*,"Position of the last element ", my_ring%end%T2%pos

       do ni=1, my_ring%end%T2%pos

          if ( .not. associated(CURR_SLICE%PARENT_FIBRE, PREV_SLICE%PARENT_FIBRE) ) then
             e = e + 1
             p=>p%next
          endif


          !          call track_beam(my_ring,TheBeam,getintstate(), pos1=ni, pos2=ni+1)
          call track_beam(my_ring,TheBeam,getintstate(), node1=ni, node2=ni+1)
          pathlegth = curr_slice%s(3)

          if (getdebug() > 2 ) then
             write(6,*) e, 'l=',pathlegth
          endif

          do n=1, npart

             x = TheBeam%X(n,1:6)

             p0=(1+x(5))
             pz=sqrt(p0**2 - x(2)**2 - x(4)**2)
             p0 = p0*p%mag%p%p0c
             xp = x(2)/pz
             yp = x(4)/pz


             !             write(441,'(i8, 1x, a16, i4 ,1x, 2f8.4, 1x, 6f8.4)') ni, p%mag%name, e,&
             !                    pathlegth, TheBeam%X(n,7), &
             !     x(1), xp , x(3), yp , x(5), p0 , x(6)

             write(mf,'(i8,1x, a16, 1x, 3i4, 1x,2f8.4, 1x, 7f12.8)' ) ni, p%mag%name, e, n, t, &
                  pathlegth, TheBeam%X(n,7), &
                  x(1), xp , x(3), yp , x(5), p0 , x(6)


             if (rplot) then
                !For thin tracking I still do not know how to get global coordinates
                !   gcs = my_false !For thin tracking
                !   if (gcs) then
                !      !                write(6,'(a12,3f8.4)') "Magnet B ", p%mag%p%f%b(1), p%mag%p%f%b(2), p%mag%p%f%b(3)
                !      gposx = x(1)*p%chart%f%exi(1,1) + x(3)*p%chart%f%exi(1,2) + x(6)*p%chart%f%exi(1,3)
                !      gposy = x(1)*p%chart%f%exi(2,1) + x(3)*p%chart%f%exi(2,2) + x(6)*p%chart%f%exi(2,3)
                !      gposz = x(1)*p%chart%f%exi(3,1) + x(3)*p%chart%f%exi(3,2) + x(6)*p%chart%f%exi(3,3)
                !      !                write(6,'(a12,3f8.4)') " Rotated ", gposx,gposy,gposz
                !      gposx = gposx + p%chart%f%b(1)
                !      gposy = gposy + p%chart%f%b(2)
                !      gposz = gposz + p%chart%f%b(3)
                !
                !      write(6,'(a12, 2i6,3f8.4)') p%mag%name, n,e, gposx,gposy,gposz
                !      call plottrack(n, e, t, gposx, xp , gposy, yp , x(5), p0 , gposz)
                !   else
                call plottrack(n, e, t, x(1), xp , x(3), yp , x(5), p0 , x(6))
                !   endif
             endif

             if ( observedelements(e) .gt. 0 ) then
                if ( associated(CURR_SLICE, p%t2 ) ) then
                   print*, "Sending to table", n, e, pathlegth
                   call putintracktable(n,t,observedelements(e),x(1), xp , x(3), yp , x(6), x(5), pathlegth, p0)
                endif
             endif
             !fields in the table         "number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e"

          enddo

          if (associated(CURR_SLICE%next)) then
             PREV_SLICE => CURR_SLICE
             CURR_SLICE => CURR_SLICE%next
          else
             exit;
          endif

       enddo !over elements

    enddo !loop over turns

    close(mf)

    if (rplot) call rplotfinish()
    call deletetrackstrarpositions()

    c_%x_prime=.false.

    CALL KILL_BEAM(TheBeam)

    deallocate (observedelements)
    !==============================================================================
  end subroutine ptc_track_everystep
  !_________________________________________________________________________________


  subroutine putinstatustable (npart,turn,elno,elna,spos,stat,x,xini,e,mf)
    implicit none
    include 'name_len.fi'
 !   integer  :: npart,turn,nobs,stat,mf,elno
    integer  :: npart,turn,stat,mf,elno
!    real(kind(1d0)) :: tt
!    character*36 table_puttab
    character*36 table
    character(name_len) elna
    !hbu
    real (dp)            :: x(1:6)
    real (dp)            :: xini(1:6)
    real(dp) :: spos,e
!    integer :: get_option
    !hbu
    data table        / 'tracksumm        ' /
      
    write(mf,*) npart,spos,turn,elno,elna,xini, x, e,stat
    
    !"number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e",
    doublenum = npart
    call double_to_table(table, 'number ' , doublenum)
    doublenum = turn
    call double_to_table(table, 'turn ' , doublenum)
    doublenum = x(1)
    call double_to_table(table, 'x ' , doublenum)
    doublenum = x(2)
    call double_to_table(table, 'px ' , doublenum)
    doublenum = x(3)
    call double_to_table(table, 'y ' , doublenum)
    doublenum = x(4)
    call double_to_table(table, 'py ' , doublenum)
    doublenum = x(6)
    call double_to_table(table, 't ' , doublenum)
    doublenum = x(5)
    call double_to_table(table, 'pt ' , doublenum)
    doublenum = spos
    call double_to_table(table, 's ' , doublenum)
    doublenum = e
    call double_to_table(table, 'e ' , doublenum)

    call augment_count(table)
    
  end subroutine putinstatustable

  !_________________________________________________________________________________

  subroutine putintracktable (npart,turn,nobs,x,px,y,py,t,pt,spos,e)
    implicit none
    !--- purpose: enter particle coordinates in table                      *
    !    input:                                                            *
    !    npart  (int)           particle number                            *
    !    turn   (int)           turn number                                *
    !    nobs   (int)           observation point number                   *
    !----------------------------------------------------------------------*

    !vvk
    !      real(dp) :: tmp_coord_array(lnv), tmp_norm_array(lnv), tmp_norm
    integer  :: npart,turn,nobs
    real(kind(1d0)) :: tt
    character*36 table_puttab
    character*36 table
    !hbu
    real(dp) :: x,px,y,py,t,pt
    real(dp) :: spos,e
    integer :: get_option
    !hbu
    data table_puttab / 'track.obs$$$$.p$$$$' /
    data table        / 'trackone           ' /

    if ( get_option('onetable ') .ne. 0 ) then
       table='trackone'
       table(9:9)= achar(0)
    else
       table=table_puttab
       write(table(10:13), '(i4.4)') nobs
       write(table(16:19), '(i4.4)') npart
    endif


    tt = turn


    call double_to_table(table, 'turn ', tt)
    doublenum = x
    call double_to_table(table, 'x ' , doublenum)

    doublenum = px
    call double_to_table(table, 'px ', doublenum)

    doublenum = y
    call double_to_table(table, 'y ' , doublenum)

    doublenum = py
    call double_to_table(table, 'py ', doublenum)

    doublenum = t
    call double_to_table(table, 't ' , doublenum)

    doublenum = pt
    call double_to_table(table, 'pt ', doublenum)

    doublenum = spos
    call double_to_table(table, 's ' , doublenum)

    doublenum = e
    call double_to_table(table, 'e ' , doublenum)
    call augment_count(table)

  end subroutine putintracktable

  !_________________________________________________________________________________


  subroutine ptc_trackline(nobs)
    ! subroutine that performs tracking with acceleration
    ! it is called as a result of ptc_trackline MAD-X command

    implicit none
    integer, intent (IN) :: nobs ! the maximum number of observation points >=1
    INTEGER, ALLOCATABLE :: observedelements(:)
    integer  :: charge    ! charge of an accelerated particle
    type(fibre), pointer :: p
    real (dp)            :: x(1:6)
    real (dp)            :: xini(1:6)
    !    real (dp)            :: polarx(1:6)   ! track vector -
    real (dp)            :: xp, yp, pz, p0
    real (dp)            :: pathlegth = zero
    integer              :: npart = 1
    integer              :: n = 1
    integer              :: nturns = 1
    integer              :: t = 1
    logical(lp)          :: gcs
    logical(lp)          :: rplot
    real (dp)            :: gposx, gposy, gposz
    integer              :: e
    integer              :: apertflag
    character(200)       :: whymsg
    integer              :: why(9)
    !    integer              :: rplotno
    integer              :: obspointnumber ! observation point number in c-code
    integer              :: getnumberoftracks !function
    type(internal_state)  :: intstate
    real(kind(1d0))      :: get_value
    integer, external    :: restart_sequ, & !  restart beamline and return number of beamline node
         advance_node    !  advance to the next node in expanded sequence
    !                    !  =0 (end of range), =1 (else)
    REAL(KIND(1d0)), external :: node_value  !/*returns value for parameter par of current element */
    integer              :: mf
    type(work)           :: fen      ! Fibre ENergy
    !------------------------------------------------------
    !initialization
    npart = 1
    n = 1
    t = 1
    !------------------------------------------------------

    if(universe.le.0) then
       call fort_warn('return from ptc_trackline: ',' no universe created')
       print *, nobs
       return
    endif
    if(index_mad.le.0) then
       call fort_warn('return from ptc_trackline: ',' no layout created')
       return
    endif

    nturns = get_value('ptc_trackline ','turns ')
    if (getdebug() > 2) print *, 'ptc_trackline, nturns = ', nturns

    if ( (nturns > 1) .and. (my_ring%closed .eqv. .false.)) then
       call fort_warn('WARNING: You can not make more than one turn in a line!', &
            'Putting number of turns to 1!')
       nturns = 1
    endif


    gcs = get_value('ptc_trackline ','gcs ') .ne. 0

    rplot = get_value('ptc_trackline ','rootntuple ') .ne. 0

    intstate = getintstate()
    if (gcs .and.  intstate%TOTALPATH==1) then
       call fort_warn("ptc_trackline","Having global coordinates and totalpath for z is sensless")
       gcs = .false.
    endif


    allocate(observedelements(1:my_ring%n)); observedelements(:)=0 ! zero means that this element is not an obs. point

    c_%x_prime=.true.

    e=restart_sequ()
    p=>my_ring%start
    do e=1, my_ring%n

       obspointnumber=node_value('obs_point ')
       IF (e.eq.1) obspointnumber=1 ! node_value gives 0 for 1st (?)

       if (obspointnumber .gt. 0) then
          if (getdebug() > 0) print *,"Element ",e," is an observation point no. ",obspointnumber
          observedelements(e) = obspointnumber
       endif

       obspointnumber=advance_node() ! c-code go to the next node -> the passed value is never used, just to shut up a compiler
       p=>p%next
    enddo


    charge = get_value('beam ', "charge ");
    if (getdebug() > 3 ) print *, 'Read charge:', charge,' layout has charge ', my_ring%start%charge

    if (cavsareset .eqv. .false.) then
       call setcavities(my_ring,maxaccel)
    endif

    if (getdebug() > 0) print *, 'reading tracks starting posiotions from table ....'

    call gettrack(1,x(1),x(2),x(3),x(4),x(6),x(5))

    if (getdebug() > 0) print *, 'reading.... Done'

    if (getdebug() > 0) then
       print *, '###################################################'
       print *, '###################################################'
       print *, '######         TRACKING WITH PTC         ##########'
       print *, '###################################################'
       print *, '###################################################'
    endif

    if (rplot) then
       call newrplot()
    endif

    call kanalnummer(mf)
    open(unit=mf,file='ptctracklinestatus.txt',POSITION='APPEND' , STATUS='UNKNOWN')

    n=1
    npart = getnumberoftracks()
    if (getdebug() > 0) print *, 'There is ', npart,' tracks'
    do n=1, npart

       pathlegth = zero

       if (getdebug() > 3 ) print *, 'Getting track ',n

       call gettrack(n,x(1),x(2),x(3),x(4),x(6),x(5))

       if (getdebug() > 0 ) write(6,'(a10,1x,i8,1x,6(f9.6,1x))') 'Track ',n,x
       xini = x
       do t=1, nturns

          p=>my_ring%start

          do e=1, my_ring%n

             !print*, p%mag%name, p%mag%P%KILL_ENT_FRINGE, p%mag%P%KILL_EXI_FRINGE,  &
             !        p%mag%P%BEND_FRINGE,  p%mag%p%PERMFRINGE, p%mag%PERMFRINGE
             !write(6,'(a10,1x,i8,1x,6(f12.9,1x))') 'Track ',n,x

             call track(my_ring,x,e,e+1,getintstate())
             if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
                call fort_warn('ptc_trackline: ','DA got unstable')
                call seterrorflag(10,"ptc_trackline ","DA got unstable ");
                fen = p;
                
                call putinstatustable(n,t,e,p%MAG%name,pathlegth,3,x,xini,fen%energy,mf)
                goto 100 !for the time being lets try next particle, 
                         !but most probably we will need to stop tracking and reinit 
	     !goto 101
                
             endif
             
             pathlegth = pathlegth + p%mag%p%ld

             if (getdebug() > 2 ) then
                write(6,*) e, 'l=',pathlegth
                write(6,'(5f8.4, f16.8)') x(1),x(2),x(3),x(4),x(5),x(6)
             endif

             p0=(1+x(5))
             pz=sqrt(p0**2 - x(2)**2 - x(4)**2)
             p0 = p0*p%mag%p%p0c
             xp = x(2)/pz
             yp = x(4)/pz

             if (rplot) then
                if (gcs) then
                   !                write(6,'(a12,3f8.4)') "Magnet B ", p%mag%p%f%b(1), p%mag%p%f%b(2), p%mag%p%f%b(3)
                   gposx = x(1)*p%chart%f%exi(1,1) + x(3)*p%chart%f%exi(1,2) + x(6)*p%chart%f%exi(1,3)
                   gposy = x(1)*p%chart%f%exi(2,1) + x(3)*p%chart%f%exi(2,2) + x(6)*p%chart%f%exi(2,3)
                   gposz = x(1)*p%chart%f%exi(3,1) + x(3)*p%chart%f%exi(3,2) + x(6)*p%chart%f%exi(3,3)
                   !                write(6,'(a12,3f8.4)') " Rotated ", gposx,gposy,gposz
                   gposx = gposx + p%chart%f%b(1)
                   gposy = gposy + p%chart%f%b(2)
                   gposz = gposz + p%chart%f%b(3)

                   if (getdebug() > 3 ) write(6,'(a12, 2i6,3f8.4)') p%mag%name, n,e, gposx,gposy,gposz

                   call plottrack(n, e, t, gposx, xp , gposy, yp , x(5), p0 , gposz)
                else
                   call plottrack(n, e, t, x(1), xp , x(3), yp , x(5), p0 , x(6))
                endif
             endif

             if ( observedelements(e) .gt. 0) then
                call putintracktable(n,t,observedelements(e),x(1), xp , x(3), yp , x(6), x(5), pathlegth, p0)
             endif
             !fields in the table         "number", "turn", "x", "px", "y", "py", "t", "pt", "s", "e"

             call produce_aperture_flag(apertflag)
             if (apertflag/=0) then
                print *, 'Particle out of aperture!'

                call ANALYSE_APERTURE_FLAG(apertflag,why)
                Write(6,*) "ptc_trackline: APERTURE error for element: ",e," name: ",p%MAG%name
                Write(6,*) "Message: ",c_%message
                write(whymsg,*) 'APERTURE error: ',why
                call fort_warn('ptc_twiss: ',whymsg)
                call seterrorflag(10,"ptc_twiss: ",whymsg);
                
                fen = p;
                call putinstatustable(n,t,e,p%MAG%name,pathlegth,1,x,xini,fen%energy,mf)
                
                goto 100 !take next track
                
             endif
             p=>p%next
          enddo !over elements

          if (apertflag/=0) then
             exit; !goes to the next particle
          endif

       enddo !loop over turns
       !                     npart,turn,elno,elna,spos,stat,x,xini,e,mf
       fen = p%previous;
       t = t - 1
       print*, t 
       call putinstatustable(n,t,e,p%previous%MAG%name,pathlegth,0,x,xini,fen%energy,mf)
      
100    continue!take next track
       
    enddo !loop over tracks
    

101 continue!finish the program
    
    if (rplot) call rplotfinish()
    call deletetrackstrarpositions()

    close(mf)

    c_%x_prime=.false.

    deallocate (observedelements)
    !==============================================================================
  end subroutine ptc_trackline




  !_________________________________________________________________________________


end module madx_ptc_trackline_module


!              if (getdebug() > 3) then
!                 write(6,*) p%mag%name
!                 write(6,'(a12,3f8.4)') "Chart  B ", p%chart%f%b(1), p%chart%f%b(2), p%chart%f%b(3)
!                 write(6,'(a12,3f8.4)') "Magnet B ", p%mag%p%f%b(1), p%mag%p%f%b(2), p%mag%p%f%b(3)
!                 write(6,'(a12,3f8.4)') "Chart Exi1 ", p%chart%f%exi(1,1), p%chart%f%exi(1,2), p%chart%f%exi(1,3)
!                 write(6,'(a12,3f8.4)') "Chart Exi2 ", p%chart%f%exi(2,1), p%chart%f%exi(2,2), p%chart%f%exi(2,3)
!                 write(6,'(a12,3f8.4)') "Chart Exi2 ", p%chart%f%exi(3,1), p%chart%f%exi(3,2), p%chart%f%exi(3,3)
!                 write(6,'(a12,3f8.4)') "mag Exi1 ", p%mag%p%f%exi(1,1), p%mag%p%f%exi(1,2), p%mag%p%f%exi(1,3)
!                 write(6,'(a12,3f8.4)') "mag Exi2 ", p%mag%p%f%exi(2,1), p%mag%p%f%exi(2,2), p%mag%p%f%exi(2,3)
!                 write(6,'(a12,3f8.4)') "mag Exi2 ", p%mag%p%f%exi(3,1), p%mag%p%f%exi(3,2), p%mag%p%f%exi(3,3)
!              endif
!
