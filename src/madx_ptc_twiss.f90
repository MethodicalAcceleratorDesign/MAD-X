module madx_ptc_twiss_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! madx_ptc_distrib module
  ! Piotr K. Skowronski , Frank Schmidt (CERN)
  !
  ! This module contains service for twiss distributions with PTC
  use madx_ptc_module
  use madx_ptc_intstate_module, only : getdebug
  USE madx_ptc_setcavs_module
  USE madx_ptc_knobs_module
  USE madx_ptc_distrib_module

  implicit none

  save
  
  
    
  !============================================================================================
  !  PUBLIC INTERFACE
  public                         :: twiss,ptc_twiss



  !============================================================================================
  !  PRIVATE
  !    datta structures

  !PSk 2011.01.05 goes global to the modules so the slice tracking produces it for the summ table
  type(real_8)            :: theTransferMap(6)
  type(universal_taylor)  :: unimap(6)

  type twiss

     logical(lp) nf
     real(dp), dimension(3,3) ::  beta,alfa,gama
     real(dp), dimension(3,3) ::  beta_p,alfa_p,gama_p ! derivatives of the above w.r.t delta_p
     real(dp), dimension(3)   ::  mu
     real(dp), dimension(6)   ::  disp
     real(dp), dimension(6)   ::  disp_p ! derivatives of the dispersion w.r.t delta_p
     real(dp), dimension(6)   ::  disp_p2 ! second derivatives of the dispersion w.r.t delta_p
     real(dp), dimension(6)   ::  disp_p3 ! third order derivatives of dispersion w.r.t delta_p
     real(dp), dimension(3)   ::  tune
     real(dp), dimension(6,6) ::  eigen
  end type twiss

  interface assignment (=)
     module procedure equaltwiss
     module procedure zerotwiss
     module procedure normalform_normalform
  end interface assignment (=)

  interface alloc
     module procedure alloctwiss
  end interface alloc

  interface kill
     module procedure killtwiss
  end interface kill

  type(normalform)                      :: Normal
  real(dp), private, dimension(2,ndim2) :: angp
  real(dp), private, dimension(ndim2)   :: dicu
  integer,  private, parameter          :: ndd=ndim2
  integer,  private, dimension(4)       :: iia,icoast
  integer,  private                     :: np

  !new lattice function
  real(dp), private, dimension(3)       :: testold
  real(dp), private, dimension(3)       :: phase

  character(len=5), private, dimension(5), parameter :: str5 = (/'10000','01000','00100','00010','00001'/)

  integer, private, allocatable          :: J(:)
  integer, private, dimension(6)         :: j5 = (/0,0,0,0,1,0/)
  integer, private, dimension(6)         :: j6 = (/0,0,0,0,0,1/)
  integer, private, dimension(6,6)       :: fo = &
       reshape(          (/1,0,0,0,0,0,&
       0,1,0,0,0,0,&
       0,0,1,0,0,0,&
       0,0,0,1,0,0,&
       0,0,0,0,1,0,&
       0,0,0,0,0,1 /), &
       (/6,6/) )

  logical     :: slice_magnets, center_magnets, deltap_dependency, isRing

  logical     :: resetBetaExtrema;
  

  real(dp)    :: minBeta(3,3) ! jluc: to store extremas of Twiss functions (show-up in summary table
  real(dp)    :: maxBeta(3,3) ! jluc: to store extremas of Twiss functions (show-up in summary table)
  real(dp)    :: minBetX      ! Edwards-Teng betas
  real(dp)    :: maxBetX      ! Edwards-Teng betas
  real(dp)    :: minBetY      ! Edwards-Teng betas
  real(dp)    :: maxBetY      ! Edwards-Teng betas
  real(dp)    :: minDisp(4)  
  real(dp)    :: maxDisp(4)  

  logical     :: resetOrbitExtrema
  real(dp)    :: minOrbit(6)  
  real(dp)    :: maxOrbit(6)
  real(dp)    :: sum2Orbit(6) ! sum of squares for rms calculatios
  integer     :: nobsOrbit    ! counter of observation points for rms calculation
  
  ! slice tracking displays many points at the same position
  ! for rms only the last should be taken 
  real(dp)    :: prevOrbit(6)
  real(dp)    :: prevS(6)
  
  character(2000), private  :: whymsg
  
  
  !============================================================================================
  !  PRIVATE
  !    routines
  private zerotwiss,equaltwiss,killtwiss



contains
  !____________________________________________________________________________________________

  subroutine equaltwiss(s1,Y)
    implicit none
    type(twiss), intent(inout)::s1
    type(real_8), intent(in)::Y(6)
    integer jj,i,k, ndel, n
    real(dp) :: lat(0:6,6,3)
    real(dp) :: test, dph
    real(dp) :: epsil=1e-12  !
    integer  :: J(lnv)


    lat = zero

    n=3  ! 1 2 3 are tunes

    ndel=0
    if(c_%ndpt/=0 .and. isRing )  then
      ! in case of icase=56 the closed solution in long is not searched
       !print*, "We are at the mode 6D + nocav"
       ndel=1  !this is 6D without cavity (MADX icase=56) for a ring
    endif
    ! calculation of alpha, beta, gamma 
    J=0;
    
    !print*, "EqualTwiss nd2=",c_%nd2
    
    do i=1,c_%nd2-2*ndel
       do jj=i,c_%nd2-2*ndel
          do k=1,c_%nd-ndel
             n=n+1
             J(2*k-1)=1
             lat(i,jj,k)=              (Y(i)%t.sub.J)*(Y(jj)%t.sub.J)
             J(2*k-1)=0
             J(2*k)=1
             lat(i,jj,k)=lat(i,jj,k) + (Y(i)%t.sub.J)*(Y(jj)%t.sub.J)
             lat(jj,i,k)=lat(i,jj,k)
             J(2*k)=0
          enddo
       enddo
    enddo

    J=0
    !here ND2=4 and delta is present      nd2=6 and delta is a constant
    !      print*,"nv",c_%nv,"nd2",c_%nd2,"np",c_%np,"ndpt",c_%ndpt ,"=>",c_%nv-c_%nd2-c_%np
    if( (c_%npara==5)       .or.  (c_%ndpt/=0) ) then
       !when there is no cavity it gives us dispersions
       do i=1,4
          lat(0,i,1)=(Y(i)%t.sub.J5)
       enddo
    elseif (c_%nd2 == 6) then
       do i=1,4
          lat(0,i,1) =              (Y(i)%t.sub.J5)*(Y(6)%t.sub.J6)
          lat(0,i,1) = lat(0,i,1) + (Y(i)%t.sub.J6)*(Y(5)%t.sub.J5)
       enddo
    else
       do i=1,4
          lat(0,i,1)=zero
       enddo
    endif


    !!!!!!!!!!!!!!!!
    ! phase advance!
    !!!!!!!!!!!!!!!!

    k = 2
    if(c_%nd2==6.and.c_%ndpt==0) k = 3

    j=0
    do i=1, k
       jj=2*i -1
       TEST=ATAN2((Y(2*i -1).SUB.fo(2*i,:)),(Y(2*i-1).SUB.fo(2*i-1,:)))/TWOPI
       if (i == 3) then
          TEST = ATAN2((Y(6).SUB.fo(5,:)),(Y(6).SUB.fo(6,:)))/TWOPI
       endif


       IF(TEST<zero.AND.abs(TEST)>EPSIL)TEST=TEST+one
       DPH=TEST-TESTOLD(i)
       IF(DPH<zero.AND.abs(DPH)>EPSIL) DPH=DPH+one
       IF(DPH>half) DPH=DPH-one

       PHASE(i)=PHASE(i)+DPH
       TESTOLD(i)=TEST

    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,c_%nd
       do k=1,c_%nd
          s1%beta(k,i)= lat(2*k-1,2*k-1,i)
          s1%alfa(k,i)=-lat(2*k-1,2*k  ,i)
          s1%gama(k,i)= lat(2*k  ,2*k  ,i)
       enddo
    enddo

    ! --- derivatives of the Twiss parameters w.r.t delta_p
    if (deltap_dependency) then
       if( (c_%npara==5) .or. (c_%ndpt/=0) ) then ! condition to be checked
          call computeDeltapDependency(y,s1)
       endif
    endif
    ! ---


    !when there is no cavity it gives us dispersions
    do i=1,c_%nd2-2*ndel
       s1%disp(i)=lat(0,i,1)
    enddo

    if (c_%nd == 3) then
       do i=1,c_%nd
          test = s1%beta(3,i)
          s1%beta(3,i) = s1%gama(3,i)
          s1%gama(3,i) = test
       enddo
    endif


    s1%mu=phase

    do k=1,3
       do i=1,6
          s1%eigen(k*2-1,i) = Y(k*2-1).sub.fo(i,:)
          s1%eigen(k*2  ,i) = Y(k*2  ).sub.fo(i,:)
       enddo
    enddo

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
       write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost(:len_trim(messagelost))
       call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
       call seterrorflag(10,"ptc_twiss ",whymsg);
       return
    endif

  end subroutine equaltwiss



  subroutine  alloctwiss(s1)
    implicit none
    type (twiss),intent(inout)::s1

    s1=0
  end subroutine alloctwiss
  !_________________________________________________________________

  subroutine  killtwiss(s1)
    implicit none
    type (twiss),intent(inout)::s1

    s1%nf=.false.
    s1%beta(:,:)=zero
    s1%beta_p(:,:)=zero
    s1%alfa(:,:)=zero
    s1%alfa_p(:,:)=zero
    s1%gama(:,:)=zero
    s1%gama_p(:,:)=zero
    s1%mu(:)=zero
    s1%disp(:)=zero
    s1%disp_p(:)=zero
    s1%disp_p2(:)=zero
    s1%disp_p3(:)=zero
    s1%tune(:)=zero
    s1%eigen(:,:)=zero
  end subroutine killtwiss
  !_________________________________________________________________

  !________________________________________________________________________________

  subroutine zerotwiss(s1,i)
    implicit none
    type(twiss), intent(inout)::s1
    integer, intent(in)::i

    if(i==0) then

       call liepeek(iia,icoast)
       NP=iia(2)-c_%nd2

       s1%nf=.false.
       s1%beta(:,:)=zero
       s1%beta_p(:,:)=zero
       s1%alfa(:,:)=zero
       s1%alfa_p(:,:)=zero
       s1%gama(:,:)=zero
       s1%gama_p(:,:)=zero
       s1%mu(:)=zero
       s1%disp(:)=zero
       s1%disp_p(:)=zero
       s1%disp_p2(:)=zero
       s1%disp_p3(:)=zero
       s1%tune(:)=zero
       s1%eigen(:,:)=zero
       dicu(:)=zero
       angp(:,:)=zero
    endif

  end subroutine zerotwiss




  !_________________________________________________________________

  subroutine ptc_twiss(tab_name,summary_tab_name)
    use twissafi
    implicit none
    logical(lp)             :: closed_orbit,beta_flg, slice, goslice
    integer                 :: k,i,ii
    integer                 :: no,mynd2,npara,nda,icase,flag_index,why(9),my_nv,nv_min
    integer                 :: ioptfun,iii,restart_sequ,advance_node,mf1,mf2
    integer                 :: tab_name(*)
    integer                 :: summary_tab_name(*)
    real(dp)                :: x(6)
    real(dp)                :: deltap0,deltap ,d_val
    real(kind(1d0))         :: get_value,suml,s
    integer                 :: posstart, posnow
    integer                 :: geterrorflag !C function that returns errorflag value
    type(real_8)            :: y(6)
    type(twiss)             :: tw
    type(fibre), POINTER    :: current
    type(integration_node), pointer :: nodePtr, stopNode
    type(work)              :: startfen !Fibre energy at the start
    real(dp)                :: r,re(ndim2,ndim2),dt
    logical(lp)             :: initial_matrix_manual, initial_matrix_table, initial_map_manual
    logical(lp)             :: initial_distrib_manual, initial_ascript_manual, writetmap
    logical(lp)             :: maptable
    logical(lp)             :: ring_parameters  !! forces isRing variable to true, i.e. calclulation of closed solution
    integer                 :: rmatrix
    real(dp)                :: emi(3)
    logical(lp)             :: isputdata
    integer                 :: countSkipped
    character(48)           :: summary_table_name
    character(12)           :: tmfile='transfer.map'
    character(48)           :: charconv !routine
    
    if(universe.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_twiss: ',' no universe created')
       call seterrorflag(1,"ptc_twiss ","no universe created till now");
       return
    endif
    if(index_mad.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_twiss: ',' no layout created')
       call seterrorflag(2,"ptc_twiss ","no layout created till now");
       return
    endif

    call resetBetaExtremas()

    !skipnormalform = my_false
    countSkipped = 0
    
    !all zeroing
    testold = zero
    phase = zero
    do i=1,6
       unimap(i) = zero
    enddo

    if (getdebug() > 1) then
        print*,"ptc_twiss"
    endif
    call momfirstinit()

    !------------------------------------------------------------------------------
    table_name = charconv(tab_name)
    summary_table_name = charconv(summary_tab_name)

    if (getdebug() > 1) then
       print*,"ptc_twiss: Table name is ",table_name
       print*,"ptc_twiss: Summary table name is ", summary_table_name
    endif

    call cleartables() !defined in madx_ptc_knobs

    call kill_para(my_ring) !removes all the previous parameters

    nda = getnknobsall() !defined in madx_ptc_knobs
    suml=zero

    no = get_value('ptc_twiss ','no ')
    if ( no .lt. 1 ) then
       call fort_warn('madx_ptc_twiss.f90 <ptc_twiss>:','Order in twiss is smaller then 1')
       print*, "Order is ", no
       return
    endif

    icase = get_value('ptc_twiss ','icase ')

    deltap0 = get_value('ptc_twiss ','deltap ')

    rmatrix = get_value('ptc_twiss ','rmatrix ')

    deltap = zero

    call my_state(icase,deltap,deltap0)

    CALL UPDATE_STATES

    ! check that deltap-dependency selected if icase=5, which stands for
    ! 4 dimensions and deltap/p as an external parameter.
    ! indeed the simplified formulas we applied for derivation w.r.t deltap assume
    ! that deltap is an externally supplied constant parameter, which is no longer
    ! the case when icase = 6 as deltap then becomes a phase-space state-variable.
    ! and for icase=4, there is no dispersion.

    deltap_dependency = get_value('ptc_twiss ','deltap_dependency ') .ne. 0

    if (deltap_dependency) then
    
      if (.not.((icase.eq.5) .or. (icase.eq.56))) then
         call fort_warn('ptc_twiss: ','derivation w.r.t deltap assume deltap is fixed parameter')
         call fort_warn('ptc_twiss: ','derivation w.r.t deltap assume icase=5 (formula differentiation) or icase=56 (from map)')
      endif

     if (no < 2) then
       call fort_warn("ptc_twiss ","For dispersion's 1st order derivatives w.r.t delta-p there must be no>=2 ")
     endif

     if (no < 3) then
       call fort_warn("ptc_twiss ","For dispersion's 2nd order derivatives w.r.t delta-p there must be no>=3 ")
     endif

     if (no < 4) then
       call fort_warn("ptc_twiss ","For dispersion's 3rd order derivatives w.r.t delta-p there must be no>=4 ")
     endif
    
    endif
    
    !############################################################################
    !############################################################################
    !############################################################################

    slice_magnets = get_value('ptc_twiss ','slice_magnets ') .ne. 0
    center_magnets = get_value('ptc_twiss ','center_magnets ') .ne. 0

    slice = slice_magnets .or. center_magnets
    
    if ( slice) then
     call make_node_layout(my_ring) 
     call getBeamBeam()
    endif 

    !############################################################################
    !############################################################################
    !############################################################################

    
    x(:)=zero
    if(mytime) then
       call Convert_dp_to_dt (deltap, dt)
    else
       dt=deltap
    endif
    if(icase.eq.5 .or. icase.eq.56) x(5)=dt
    
    closed_orbit = get_value('ptc_twiss ','closed_orbit ') .ne. 0

    if( closed_orbit .and. (icav .gt. 0) .and. (my_ring%closed .eqv. .false.)) then
       call fort_warn('return from ptc_twiss: ',' Closed orbit requested on not closed layout.')
       call seterrorflag(3,"ptc_twiss ","Closed orbit requested on not closed layout.")
       return
    endif

    if(closed_orbit) then

       if ( .not. c_%stable_da) then
          call fort_warn('ptc_twiss: ','DA got unstable even before finding closed orbit')
          call seterrorflag(10,"ptc_twiss ","DA got unstable even before finding closed orbit")
          call aafail('ptc_twiss: ','DA got unstable even before finding closed orbit. program stops')
          !          return
       endif

       ! pass starting point for closed orbit search
       x(1)=get_value('ptc_twiss ','x ')
       x(2)=get_value('ptc_twiss ','px ')
       x(3)=get_value('ptc_twiss ','y ')
       x(4)=get_value('ptc_twiss ','py ')
       x(6)=get_value('ptc_twiss ','t ')
       x(5)=x(5)+get_value('ptc_twiss ','pt ')

       

       if (getdebug() > 2) then
         print*, "Looking for orbit"
         call print(default,6)
       endif
       
       
       if ( slice )  then
         call FIND_ORBIT_x(my_ring,x,default,c_1d_7)
       else
         call find_orbit(my_ring,x,1,default,c_1d_7)
       endif
       
       if ( .not. check_stable) then
          write(whymsg,*) 'DA got unstable during closed orbit search: PTC msg: ',messagelost(:len_trim(messagelost))
          call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
          call seterrorflag(10,"ptc_twiss ",whymsg);
          return
          !          return
       endif
       
      ! print*, "From closed orbit", w_p%nc
      ! if ( w_p%nc .gt. 0) then
      !   do i=1,w_p%nc
      !      call fort_warn('ptc_twiss: ',w_p%c(i))
      !      call seterrorflag(10,"ptc_twiss ",w_p%c(i));
      !   enddo
      ! endif    
       
      if (getdebug() > 1) then
         CALL write_closed_orbit(icase,x)
      endif

    elseif((my_ring%closed .eqv. .true.) .and. (getdebug() > 1)) then
       print*, "Closed orbit specified by the user!"
       !CALL write_closed_orbit(icase,x) at this position it isn't read
    endif

    mynd2 = 0
    npara = 0
    
    writetmap = get_value('ptc_twiss ','writetmap ') .ne. 0
    
    !this must be before initialization of the Bertz

    initial_distrib_manual = get_value('ptc_twiss ','initial_moments_manual ') .ne. 0
    if (initial_distrib_manual) then
       if (getdebug() > 1) then
          print*,"Initializing map with initial_moments_manual=true"
       endif
       call readinitialdistrib()
    endif


    call init(default,no,nda,BERZ,mynd2,npara)

    !This must be before init map
    call alloc(y)
    y=npara
    Y=X

    !    if (maxaccel .eqv. .false.) then
    !      cavsareset = .false.
    !    endif

    if ( (cavsareset .eqv. .false.) .and. (my_ring%closed .eqv. .false.) ) then

       call setcavities(my_ring,maxaccel)
       if (geterrorflag() /= 0) then
          return
       endif
    endif

    call setknobs(my_ring)

    
    
    call alloc(theTransferMap)
    theTransferMap = npara
    theTransferMap = X



    !############################################################################
    !############################################################################
    !############################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  INIT Y that is tracked          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call initmap(dt,slice)

    if (geterrorflag() /= 0) then
       !if arror occured then return
       return
    endif

    
    !############################################################################
    !############################################################################
    !############################################################################


    
    call alloc(tw)

    !Y

    !the initial twiss is needed to initialize propely calculation of some variables f.g. phase advance
    tw=y
    phase = zero !we have to do it after the very initial twiss params calculation above

    current=>MY_RING%start
    startfen = 0
    startfen = current!setting up start energy for record
    suml=zero
    iii=restart_sequ()
    print77=.false.
    read77=.false.

    if (getdebug() > 2) then
       call kanalnummer(mf1)
       open(unit=mf1,file='ptctwiss.txt')
       print *, "ptc_twiss: internal state is:"
       call print(default,6)
    endif

    call killsavedmaps() !delete all maps, if present
    mapsorder = 0 !it is set at the end, so we are sure the twiss was successful

    savemaps = get_value('ptc_twiss ','savemaps ') .ne. 0

    if (savemaps) then
       allocate(maps(MY_RING%n))
       do i=1,MY_RING%n
          do ii=1,6
             call alloc(maps(i)%unimap(ii),0,0)
             maps(i)%unimap(ii) = zero !this initializes and allocates the variables
          enddo
       enddo
    else
       nullify(maps) !assurance
    endif
    
    resetBetaExtrema = .true.;
    resetOrbitExtrema = .true.;
    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                              !
   !  T H E   M A I N   L O O P   !
   !                              !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,MY_RING%n

      if (getdebug() > 1) then
         write(6,*) ""
         write(6,*) "##########################################"
         write(6,'(i4, 1x,a, f10.6)') i,current%mag%name, suml
         write(6,'(a1,a,a1)') ">",current%mag%vorname,"<"
         write(6,'(a, f15.6, a)') "Ref Momentum ",current%mag%p%p0c," GeV/c"
         !          if (associated(current%mag%BN)) write(6,*) "k1=", current%mag%BN(2)
      endif
      
      
      ! Can not do this trick because beam beam can be defined within those elements
      ! so even stupid markers will occur at least twice in the twiss table
      ! skowron 2012.07.03
      !if (slice)  then
         goslice  = .true.
      !   if (current%mag%kind==kind0) then ! this is a MARKER
      !     goslice = .false.
      !   elseif(current%mag%kind==kind1) then ! this is a DRIFT, they go in one step anyway
      !     goslice = .false.
      !   elseif(current%mag%kind==kind11) then ! this is a MONITOR, they go in one step anyway
      !     goslice = .false.
      !   elseif(current%mag%kind==kind12) then ! this is a HMONITOR, they go in one step anyway
      !     goslice = .false.
      !   elseif(current%mag%kind==kind13) then ! this is a VMONITOR, they go in one step anyway
      !     goslice = .false.
      !   elseif(current%mag%kind==kind14) then ! this is a INSTRUMENT, they go in one step anyway
      !     goslice = .false.
      !   endif          
      !endif

      if (slice .and. goslice) then

        if (getdebug() > 1) then
           write(6,*) "##### SLICE MAGNETS"
        endif

        nodePtr => current%t1

        posstart = nodePtr%pos
        !posnow = posstart;

        if (associated(current%next)) then
          !write(6,*) "                       Next node exists"
          stopNode => current%next%t1
        else
          stopNode => current%t2
        endif  

        do while ( .not. (associated(nodePtr, stopNode)) ) 

          s = nodePtr%next%s(1) ! s(1) is the total arc-length, s(3) the total integration-distance

          !I do not know what JL meant here
          if ((s .eq. 0d0) .and. (nodePtr%pos .eq. (my_ring%t%n+posstart-1))) then
             s = nodePtr%s(1) + nodePtr%next%s(5) ! s of previous node + local offset
          endif

          if (getdebug() > 2) then
             write(6,*) "##### SLICE MAGNETS NODE ",&
                      & nodePtr%pos," => ",nodePtr%pos+1," s=",s
          endif

         if (nda > 0) then
            call track_probe_x(my_ring,y,+default, & ! +default in case of extra parameters !?
                 & node1=nodePtr%pos,node2=nodePtr%pos+1)
            call track_probe_x(my_ring,theTransferMap,+default, & ! +default in case of extra parameters !?
                 & node1=nodePtr%pos,node2=nodePtr%pos+1)
          else
            call track_probe_x(my_ring,y,default, &
                 & node1=nodePtr%pos,node2=nodePtr%pos+1)
            call track_probe_x(my_ring,theTransferMap,default, &
                 & node1=nodePtr%pos,node2=nodePtr%pos+1)
          endif

          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             
             write(whymsg,*) 'DA got unstable in tracking at s= ',s, &
                             ' magnet ',i,' ', current%mag%name,' ', current%mag%vorname, &
	         ' step ',nodePtr%pos,' PTC msg: ',messagelost(:len_trim(messagelost))
             call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
             call seterrorflag(10,"ptc_twiss ",whymsg);
             
             if (getdebug() > 2) close(mf1)
             return
          endif

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             write(whymsg,*) 'APERTURE error: ',why,' s=',s,'  element: ',i,' name: ',current%MAG%name
             call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
             call seterrorflag(10,"ptc_twiss: ",whymsg);
             !          Write(6,*) why ! See produce aperture flag routine in sd_frame
             goto 100
          endif

          isputdata = .false.

          !!!!!!!!!!!!!!!

          if (center_magnets ) then
            if ( associated(nodePtr,current%tm) ) then
              if (mod(current%mag%p%nst,2)/=0) then !checking if number od slices is even
                 !here it is odd (we should be in the middle)
                 if (current%mag%L==0) then
	! this is a zero-length marker or monitor or equivalent  keep it
	isputdata = .true.
                 elseif (current%mag%kind==31) then ! this is a DRIFT
	! write(0,*) 'SKIP a drift of length',fibrePtr%mag%L,'was cut into an odd number of slices: ',fibrePtr%mag%p%nst
	countSkipped = countSkipped + 1
                 else ! neither a zero-length element, nor a drift
	! write(0,*) 'SKIP element of name',fibrePtr%mag%name,'and kind',fibrePtr%mag%kind
	countSkipped = countSkipped + 1
                 endif
              else ! this element has an even number of slices
                 isputdata = .true.
              endif
            endif  
          else
            if (nodePtr%next%cas==case0) then
              !not center_magnets, take every reasonable node
              ! this an inner integration node i.e. neither an extremity nor a fringe node, both to be discarded
               isputdata = .true.
            else
              if (getdebug() > 2) then
                 write(6,*) "         Not Saving data CASE=",nodePtr%next%cas
              endif
            endif
          endif

          if (isputdata ) then

            if (getdebug() > 2) then
               write(6,*) "             Saving data CASE=",nodePtr%next%cas
            endif

            tw = y ! set the twiss parameters, with y being equal to the A_ phase advance
            suml = s; 

            call puttwisstable(theTransferMap)
            call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap, y)

          !else
          !  write(6,*) "                                                NOT Saving data"
          endif

          nodePtr => nodePtr%next
        enddo    

        if (isputdata .eqv. .false.) then ! always save the last point if it was not yet saved
          !write(6,*) "                                                  END OF ELEMENT Saving data"

          if (getdebug() > 2) then
             write(6,*) "               Saving anyway, it is the last node"
          endif

          tw = y ! set the twiss parameters, with y being equal to the A_ phase advance
          if (s > suml) then !work around against last element having s=0
            suml = s; 
          endif

          call puttwisstable(theTransferMap)
          call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap, y)

        endif

      else
        ! ELEMENT AT ONCE MODE
        if (nda > 0) then
           !         if (getnknobis() > 0) c_%knob = my_true
           !print*, "parametric",i,c_%knob
           call track(my_ring,y,i,i+1,+default)
           call track(my_ring,theTransferMap,i,i+1,+default)
        else
           call track(my_ring,y,i,i+1, default)
           call track(my_ring,theTransferMap,i,i+1,default)
        endif


        if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
           
           write(whymsg,*) 'DA got unstable in tracking at s= ',s, &
                           ' magnet ',i,' ', current%mag%name,' ', current%mag%vorname, &
	       ' PTC msg: ',messagelost(:len_trim(messagelost))
           call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
           call seterrorflag(10,"ptc_twiss ",whymsg);
           if (getdebug() > 2) close(mf1)
           return
        endif

        call PRODUCE_APERTURE_FLAG(flag_index)
        if(flag_index/=0) then
           call ANALYSE_APERTURE_FLAG(flag_index,why)
           write(whymsg,*) 'APERTURE error: ',why,' s=',s,'  element: ',i,' name: ',current%MAG%name
           call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
           call seterrorflag(10,"ptc_twiss: ",whymsg);
           goto 100
        endif

        if (getdebug() > 2) then
           write(mf1,*) "##########################################"
           write(mf1,'(i4, 1x,a, f10.6)') i,current%mag%name, suml
           call print(y,mf1)
        
           if (current%mag%kind==kind4) then
             print*,"CAVITY at s=",suml,"freq ", current%mag%freq,&
                    	"lag ", current%mag%lag, &
                    	"volt ", current%mag%volt
            endif
        endif
        
        suml=suml+current%MAG%P%ld

        if (savemaps) then
           do ii=1,6
              maps(i)%unimap(ii) = y(ii)
           enddo
           maps(i)%s = suml
           maps(i)%name = current%mag%name
        endif

        ! compute the Twiss parameters
        tw=y

        call puttwisstable(theTransferMap)
        
        call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap,y)


      endif

      iii=advance_node()
      current=>current%next
    enddo

100 continue


    ! relocated the following here to avoid side-effect
    print77=.false.
    read77=.false.

    
    if (writetmap) then
       !'===      TRANSFER MAP      ==='
       call kanalnummer(mf2)
       open(unit=mf2,file=tmfile)
       call print(theTransferMap,mf2)
       close(mf2)
    endif
    
    if (getdebug() > 2) then


       call kanalnummer(mf2)
       ! avoid file name conflict between slice.madx and center.madx
       ! which are located under same testing directory
       if (get_value('ptc_twiss ','center_magnets ').ne.0) then
          open(unit=mf2,file='end_center.map')
       else
          open(unit=mf2,file='end.map')
       endif

       call print(y,mf2)

       close(mf2)

    endif

   
    !must be after initmap that sets the isRing
    ring_parameters = get_value('ptc_twiss ','ring_parameters ') .ne. 0
    if (ring_parameters) then
      if (getdebug() > 1) then
        write(6,*) "User forces ring parameters calculation"
      endif
      isRing = .true.
    endif


    if(isRing .eqv. .true.) then
       call oneTurnSummary(theTransferMap ,y, x, suml)
    else
       print*, "Reduced SUMM Table (Inital parameters specified)"
       call onePassSummary(theTransferMap , x, suml)
    endif


    maptable = get_value('ptc_twiss ','maptable ') .ne. 0
    if(maptable) then
       call makemaptable(theTransferMap,no)
    endif


    call set_option('ptc_twiss_summary ',1)
    

    if (getdebug() > 1) then
       write(6,*) "##########################################"
       write(6,*) "##########################################"
       write(6,*) "###  END  OF  PTC_TWISS            #######"
       write(6,*) "##########################################"
       write(6,*) "##########################################"
       !          if (associated(current%mag%BN)) write(6,*) "k1=", current%mag%BN(2)
    endif

    c_%watch_user=.false.

    call kill(tw)
    CALL kill(y)

    CALL kill(theTransferMap)

    do i=1,6
       call kill(unimap(i))
    enddo

    call finishknobs()

    if (savemaps) then  !do it at the end, so we are sure the twiss was successful
       mapsorder = no
       mapsicase = icase
       if (getnmoments() > 0) call ptc_moments(no*2) !calcualate moments with the maximum available order
    endif


! f90flush is not portable, and useless...
!    call f90flush(20,my_false)

    if (getdebug() > 2) close(mf1)

    !****************************************************************************************
    !*********  E N D   O F   PTC_TWISS      ************************************************
    !****************************************************************************************
    !________________________________________________________________________________________

  contains  ! what follows are internal subroutines of ptc_twiss
    !____________________________________________________________________________________________

    subroutine initmap(dt,slice)
      implicit none
      integer     :: double_from_table_row
      integer     :: mman, mtab, mascr, mdistr !these variable allow to check if the user did not put too many options
      integer     :: mmap
      real(dp)    :: dt
      logical(lp) :: slice
      
      beta_flg = (get_value('ptc_twiss ','betx ').gt.0) .and. (get_value('ptc_twiss ','bety ').gt.0)

      mman  = get_value('ptc_twiss ','initial_matrix_manual ')
      mtab  = get_value('ptc_twiss ','initial_matrix_table ')
      mascr = get_value('ptc_twiss ','initial_ascript_manual ')
      mdistr = get_value('ptc_twiss ','initial_moments_manual ')
      mmap = get_value('ptc_twiss ','initial_map_manual ')


      isRing = .false. ! set to true in the case the map



      ! is calculated over a ring, about the closed orbit. Later-on in the same subroutine.


      initial_matrix_manual = mman .ne. 0
      initial_map_manual = mmap .ne. 0
      initial_matrix_table = mtab .ne. 0
      initial_ascript_manual = mascr .ne. 0


      if ( (mman + mtab + mascr + mdistr + mmap) > 1) then
         call seterrorflag(11,"ptc_twiss ","Ambigous command options");
         print*, "Only one of the following switches might be on:"
         print*, "initial_matrix_manual  = ", initial_matrix_manual
         print*, "initial_map_manual     = ", initial_map_manual
         print*, "initial_matrix_table   = ", initial_matrix_table
         print*, "initial_ascript_manual = ", initial_ascript_manual
         print*, "initial_moments_manual = ", mdistr
      endif

      !        print*, "initial_distrib_manual is ",initial_distrib_manual

      if(initial_matrix_table) then
         k = double_from_table_row("map_table ", "nv ", 1, doublenum)
         if(k.ne.-1) then
            call liepeek(iia,icoast)
            my_nv=int(doublenum)
            nv_min=min(c_%npara,my_nv)
            if (getdebug() > 2) then
              print*,"NV from table ", my_nv, " this command defined nv ",c_%npara," Using ", nv_min
            endif
         else
            if (getdebug() > 2) then
              print*,"Can not read NV from map_table. Exiting."
            endif
            
            call seterrorflag(10,"ptc_twiss initmap","Can not read NV from map_table. Exiting. ");
            
            return
         endif
      endif

      if(initial_matrix_table) then

         if (getdebug() > 1) then
            print*,"Initializing map with initial_matrix_table=true"
         endif
         call readmatrixfromtable()

         if (geterrorflag() /= 0) then
            return
         endif

      elseif(initial_ascript_manual) then

         if (getdebug() > 1) then
            print*,"Initializing map with initial_ascript_manual=true"
         endif
         call readinitialascript()
         if (geterrorflag() /= 0) then
            return
         endif

      elseif(initial_map_manual) then
         if (getdebug() > 1) then
            print*,"Initializing map with initial_map_manual=true"
         endif
         call readinitialmap()
         if (geterrorflag() /= 0) then
            return
         endif

      elseif(initial_matrix_manual) then

         if (getdebug() > 1) then
            print*,"Initializing map with initial_matrix_manual=true"
         endif
         call readinitialmatrix()

         if (geterrorflag() /= 0) then
            return
         endif
      elseif (initial_distrib_manual) then
         !matrix is already prepared beforehand
         if (getdebug() > 1) then
            print*,"Initializing map with initial_moments_manual=true"
            print*, "Initializing map from prepared UniTaylor"
         endif
         call readreforbit() !reads x

         do i=1, c_%nd2
            y(i) = unimap(i)
         enddo

         if (geterrorflag() /= 0) then
            return
         endif
      elseif(beta_flg) then

         if (getdebug() > 1) then
            print*,"Initializing map with initial twiss parameters"
         endif

         call readinitialtwiss(dt)

         if (geterrorflag() /= 0) then
            return
         endif
      else

         isRing = .true. ! compute momemtum compaction factor, tunes, chromaticies for ring

         if (getdebug() > 1) then
            print*,"Initializing map from one turn map: Start Map"
            call print(y,6)
         
            print*,"Tracking identity map to get closed solution. STATE:"
            call print(default,6)

         endif
         
         if (slice) then
           call track_probe_x(my_ring,y,default) !, MY_RING%start%t1,MY_RING%end%t2);
         else
           call track(my_ring,y,1,default)
         endif
                  
         if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
            write(whymsg,*) 'DA got unstable (one turn map production) at ', &
                             lost_fibre%mag%name, &
                            ' PTC msg: ',messagelost(:len_trim(messagelost))
            call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
            call seterrorflag(10,"ptc_twiss ",whymsg);
            return
         endif

         call PRODUCE_APERTURE_FLAG(flag_index)
         if(flag_index/=0) then
            call ANALYSE_APERTURE_FLAG(flag_index,why)

            write(whymsg,*) 'APERTURE unstable (one turn map production) - programs continues: ',why
            call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
            call seterrorflag(10,"ptc_twiss: ",whymsg);
            !          Write(6,*) "ptc_twiss unstable (map production)-programs continues "
            !          Write(6,*) why ! See produce aperture flag routine in sd_frame
            c_%watch_user=.false.
            CALL kill(y)
            return
         endif

         if (getdebug() > 1) then
            print*,"Initializing map from one turn map. One Turn Map"
            call print(y,6)
            
         endif

         call maptoascript()

         call reademittance()

      endif
      
    end subroutine initmap
    !____________________________________________________________________________________________

    function getdeltae()
      implicit none
      real(dp)   :: getdeltae
      type(work) :: cfen !current fibre energy

      cfen = 0 ! do not remove -> if it is removed energy is wrong because = adds energy to the previous value

      if ( (associated(current%next) .eqv. .false. ) .or. (associated( current%next, my_ring%start)) ) then
         ! if current is the last element in the sequence i.e.
         ! p%next == NULL (LINE) OR
         ! p%next points the first element (CIRCLE)
         cfen=current

         if (getdebug()>1) then
            !if it is the last element in the line
            print *, 'It is the last element  ', current%mag%name
            !(it is always marker, i.e element that does not change reference energy)
            print *, 'Its reference energy is ', cfen%p0c
         endif

         !take its reference energy
      else
         cfen=current%next      ! energy after passing this element
      endif

      getdeltae = cfen%energy/startfen%energy

      if (getdebug() > 2) then
         write(mf1,'(3(a, f12.6))') "Ref Momentum ",cfen%p0c," Energy ", cfen%energy," DeltaE ",getdeltae
      endif

    end function getdeltae
    !____________________________________________________________________________________________

    subroutine puttwisstable(transfermap)
      implicit none
      include "madx_ptc_knobs.inc"
      integer i1,i2,ii,i1a,i2a
      real(kind(1d0))   :: opt_fun(150),myx  ! opt_fun(72) -> opt_fun(81)
      ! increase to 150 to have extra space beyond what's needed to accomodate additional derivatives w.r.t. delta_p
      real(kind(1d0))   :: deltae
      type(real_8), target :: transfermap(6)
      ! added on 3 November 2010 to hold Edwards & Teng parametrization
      real(dp) :: betx,bety,alfx,alfy,R11,R12,R21,R22
      ! to convert between Ripken and Edwards-Teng parametrization
      real(dp) :: kappa,u,ax,ay,kx,ky,kxy2,usqrt,bx,by,cx,cy,cosvp,sinvp,cosvm,sinvm,cosv2,sinv2,cosv1,sinv1
      real(dp) :: deltaeValue

      if (getdebug() > 2) then
         write(mf1,*) "##########################################"
         write(mf1,*) ""
         write(mf1,'(i4, 1x,a, f10.6)') i,current%mag%name,suml
         write(mf1,*) ""
         call print(y,mf1)
      endif


      deltae = getdeltae()

      call double_to_table_curr(table_name, 's ', suml)

      doublenum = deltae * startfen%energy
      call double_to_table_curr(table_name, 'energy ', doublenum)


      opt_fun(:)=zero

      call liepeek(iia,icoast)
      allocate(j(c_%npara))
      j(:)=0
      do ii=1,c_%npara ! fish
         opt_fun(ii)=y(ii)%T.sub.j
      enddo

      call trackOrbitExtremaAndRms(opt_fun(1:6))

      ! swap 
      myx=opt_fun(6)
      opt_fun(6)=opt_fun(5)
      opt_fun(5)=myx
      deallocate(j)
      
      
      ioptfun=6
      call vector_to_table_curr(table_name, 'x ', opt_fun(1), ioptfun)


      opt_fun(1) = transfermap(1).sub.fo(1,:)
      opt_fun(2) = transfermap(1).sub.fo(2,:)
      opt_fun(3) = transfermap(1).sub.fo(3,:)
      opt_fun(4) = transfermap(1).sub.fo(4,:)
      opt_fun(5) = transfermap(1).sub.fo(6,:)
      opt_fun(6) = transfermap(1).sub.fo(5,:)


      opt_fun(7) = transfermap(2).sub.fo(1,:)
      opt_fun(8) = transfermap(2).sub.fo(2,:)
      opt_fun(9) = transfermap(2).sub.fo(3,:)
      opt_fun(10)= transfermap(2).sub.fo(4,:)
      opt_fun(11)= transfermap(2).sub.fo(6,:)
      opt_fun(12)= transfermap(2).sub.fo(5,:)

      opt_fun(13)= transfermap(3).sub.fo(1,:)
      opt_fun(14)= transfermap(3).sub.fo(2,:)
      opt_fun(15)= transfermap(3).sub.fo(3,:)
      opt_fun(16)= transfermap(3).sub.fo(4,:)
      opt_fun(17)= transfermap(3).sub.fo(6,:)
      opt_fun(18)= transfermap(3).sub.fo(5,:)

      opt_fun(19)= transfermap(4).sub.fo(1,:)
      opt_fun(20)= transfermap(4).sub.fo(2,:)
      opt_fun(21)= transfermap(4).sub.fo(3,:)
      opt_fun(22)= transfermap(4).sub.fo(4,:)
      opt_fun(23)= transfermap(4).sub.fo(6,:)
      opt_fun(24)= transfermap(4).sub.fo(5,:)


      opt_fun(25)= transfermap(6).sub.fo(1,:)
      opt_fun(26)= transfermap(6).sub.fo(2,:)
      opt_fun(27)= transfermap(6).sub.fo(3,:)
      opt_fun(28)= transfermap(6).sub.fo(4,:)
      opt_fun(29)= transfermap(6).sub.fo(6,:)
      opt_fun(30)= transfermap(6).sub.fo(5,:)


      opt_fun(31)= transfermap(5).sub.fo(1,:)
      opt_fun(32)= transfermap(5).sub.fo(2,:)
      opt_fun(33)= transfermap(5).sub.fo(3,:)
      opt_fun(34)= transfermap(5).sub.fo(4,:)
      opt_fun(35)= transfermap(5).sub.fo(6,:)
      opt_fun(36)= transfermap(5).sub.fo(5,:)

      ioptfun=36
      call vector_to_table_curr(table_name, 're11 ', opt_fun(1), ioptfun)
      
      deltae = deltae * (1.0_dp + y(5).sub.'0')

      opt_fun(beta11)= tw%beta(1,1) * deltae ! beta11=1
      opt_fun(beta12)= tw%beta(1,2) * deltae
      opt_fun(beta13)= tw%beta(1,3) * deltae
      opt_fun(beta21)= tw%beta(2,1) * deltae
      opt_fun(beta22)= tw%beta(2,2) * deltae
      opt_fun(beta23)= tw%beta(2,3) * deltae
      opt_fun(beta31)= tw%beta(3,1) * deltae
      opt_fun(beta32)= tw%beta(3,2) * deltae
      opt_fun(beta33)= tw%beta(3,3) * deltae

      opt_fun(alfa11)= tw%alfa(1,1) * deltae
      opt_fun(alfa12)= tw%alfa(1,2) * deltae
      opt_fun(alfa13)= tw%alfa(1,3) * deltae
      opt_fun(alfa21)= tw%alfa(2,1) * deltae
      opt_fun(alfa22)= tw%alfa(2,2) * deltae
      opt_fun(alfa23)= tw%alfa(2,3) * deltae
      opt_fun(alfa31)= tw%alfa(3,1) * deltae
      opt_fun(alfa32)= tw%alfa(3,2) * deltae
      opt_fun(alfa33)= tw%alfa(3,3) * deltae

      opt_fun(gama11)= tw%gama(1,1) * deltae
      opt_fun(gama12)= tw%gama(1,2) * deltae
      opt_fun(gama13)= tw%gama(1,3) * deltae
      opt_fun(gama21)= tw%gama(2,1) * deltae
      opt_fun(gama22)= tw%gama(2,2) * deltae
      opt_fun(gama23)= tw%gama(2,3) * deltae
      opt_fun(gama31)= tw%gama(3,1) * deltae
      opt_fun(gama32)= tw%gama(3,2) * deltae
      opt_fun(gama33)= tw%gama(3,3) * deltae


      ! --- derivatives of Twiss paramters w.r.t delta_p
      ! NOW why do we need to multiply by deltae, as for the other Twiss parameters?
      if (deltap_dependency) then
         opt_fun(beta11p)= tw%beta_p(1,1) * deltae
         opt_fun(beta12p)= tw%beta_p(1,2) * deltae
         opt_fun(beta13p)= tw%beta_p(1,3) * deltae
         opt_fun(beta22p)= tw%beta_p(2,1) * deltae
         opt_fun(beta22p)= tw%beta_p(2,2) * deltae
         opt_fun(beta23p)= tw%beta_p(2,3) * deltae
         opt_fun(beta32p)= tw%beta_p(3,2) * deltae
         opt_fun(beta33p)= tw%beta_p(3,3) * deltae

         opt_fun(alfa11p)= tw%alfa_p(1,1) * deltae
         opt_fun(alfa12p)= tw%alfa_p(1,2) * deltae
         opt_fun(alfa13p)= tw%alfa_p(1,3) * deltae
         opt_fun(alfa21p)= tw%alfa_p(2,1) * deltae
         opt_fun(alfa22p)= tw%alfa_p(2,2) * deltae
         opt_fun(alfa23p)= tw%alfa_p(2,3) * deltae
         opt_fun(alfa31p)= tw%alfa_p(3,1) * deltae
         opt_fun(alfa32p)= tw%alfa_p(3,2) * deltae
         opt_fun(alfa33p)= tw%alfa_p(3,3) * deltae

         opt_fun(gama11p)= tw%gama_p(1,1) * deltae
         opt_fun(gama12p)= tw%gama_p(1,2) * deltae
         opt_fun(gama13p)= tw%gama_p(1,3) * deltae
         opt_fun(gama21p)= tw%gama_p(2,1) * deltae
         opt_fun(gama22p)= tw%gama_p(2,2) * deltae
         opt_fun(gama23p)= tw%gama_p(2,3) * deltae
         opt_fun(gama31p)= tw%gama_p(3,1) * deltae
         opt_fun(gama32p)= tw%gama_p(3,2) * deltae
         opt_fun(gama33p)= tw%gama_p(3,3) * deltae
      endif
      ! --- end

      ! march 10th: do we need to multiply by deltae the following?
      opt_fun(mu1)=tw%mu(1) !* deltae
      opt_fun(mu2)=tw%mu(2) !* deltae
      opt_fun(mu3)=tw%mu(3) !* deltae

      ! write(0,*),"DEBUG = ", tw%mu(3) ! should give Qs

      opt_fun(disp1)=tw%disp(1) ! was 31 instead of 57
      opt_fun(disp2)=tw%disp(2) ! was 32 instead of 58
      opt_fun(disp3)=tw%disp(3) ! was 33 instead of 59
      opt_fun(disp4)=tw%disp(4) ! was 34 instead of 60

      ! 9 march 2009: add 4 items

      if (deltap_dependency) then
         opt_fun(disp1p) = tw%disp_p(1)
         opt_fun(disp2p) = tw%disp_p(2)
         opt_fun(disp3p) = tw%disp_p(3)
         opt_fun(disp4p) = tw%disp_p(4)
         ! 20 july 2009: add 8 items for second derivatives w.r.t deltap
         opt_fun(disp1p2) = tw%disp_p2(1)
         opt_fun(disp2p2) = tw%disp_p2(2)
         opt_fun(disp3p2) = tw%disp_p2(3)
         opt_fun(disp4p2) = tw%disp_p2(4)
         opt_fun(disp1p3) = tw%disp_p3(1)
         opt_fun(disp2p3) = tw%disp_p3(2)
         opt_fun(disp3p3) = tw%disp_p3(3)
         opt_fun(disp4p3) = tw%disp_p3(4)
      endif

      ! JLUC TODO
      ! opt_fun(61)=zero disp4 is now 61 in madx_ptc_knobs.inc
      ! jluc: left the following umodified, except 36->62
      opt_fun(62+4+8+6)=zero ! was 36 instead of 62 => on 9 march add 4 => on 3 July 2009 add 4 => on 3 November 2010 add 6
      do i1=1,6
         if(i1.le.4) then
            i1a=i1
         elseif(i1.eq.5) then
            i1a=6
         else
            i1a=5
         endif
         do i2=1,6
            if(i2.le.4) then
               i2a=i2
            elseif(i2.eq.5) then
               i2a=6
            else
               i2a=5
            endif
! idiotic counting reset (62+4+8+6) ==> 73
            ii=73+(i1a-1)*6+i2a ! was 36 instead of 62 => now (62+4) instead of 62 => 62+4+4
            opt_fun(ii)=tw%eigen(i1,i2) * deltae
            ! where do these eigen values go? no such dedicated column defined in madx_ptc_knobs.inc
            if(mytime.and.i2a.eq.6) opt_fun(ii)=-opt_fun(ii)
         enddo
      enddo

      if (getdebug() > 2)  then
         ! this part of the code is out of sync with the rest
         !j(:)=0
         write(6,'(a,1(f9.4,1x))') current%MAG%name,suml
         write(6,'(a,1(f11.7,1x))') "Delta E ", deltae
         write(6,'(a,3(i8.0,1x))')  "idxes   ", beta11,beta22,beta33
         write(6,'(a,3(f11.4,1x))')  "betas raw    ", tw%beta(1,1),tw%beta(2,2),tw%beta(3,3)
         write(6,'(a,3(f11.4,1x))')  "betas w/ener ", opt_fun(beta11),opt_fun(beta22),opt_fun(beta33)
         write(6,'(a,4(f11.4,1x))')  "dispersions  ", opt_fun(disp1),opt_fun(disp2),opt_fun(disp3),opt_fun(disp4)
         write(6,'(a,3(f11.4,1x))')  "phase adv.   ", tw%mu(1),tw%mu(2),tw%mu(3)
         write(6,'(a,4(f11.4,1x))')  "orbit transv.", y(1)%T.sub.'0',y(2)%T.sub.'0',y(3)%T.sub.'0',y(4)%T.sub.'0'
         write(6,'(a,2(f11.4,1x))')  "dp/p, T      ", y(5)%T.sub.'0',y(6)%T.sub.'0'
         
      endif

      ! the following works : we see the list of all elements in sequence - what about twiss_ptc_line & twiss_ptc_ring?
      ! jluc debug - begin
      !write(28,'(a,1(f8.4,1x))') current%MAG%name,suml
      ! jluc debug - end

      ! on July 3rd 2009, add another 8 for second/third derivatives of dispersions w.r.t. deltap
      ioptfun=81+4+8+6 !72->81 to accomodate additional derivatives w.r.t. delta_p => should one add 4 to this one, as above?
      ! actually 3*21+6 (???) elements from beta11 to include up to disp6p
      ! on november 3rd 2010, added 6

      ! overwrote the above for which I am not sure where the value comes from
      ioptfun = 79 + 36 ! 79 as for ntwisses in madx_ptc_knobs.inc + 36 eigenvalues

      ! LD: update columns access from twiss column description in mad_gcst.c
      ! WAS: fill contiguous data in one-go, from beta11 up to mu1, mu2, mu3
      ! NOW: fill contiguous data in one-go, from beta11 up to the end
      !      eigenvalue seems to be retrieved ~50 lines above...
      ioptfun = 18*3+16+3+36+1 ! = 110
      call vector_to_table_curr(table_name, 'beta11 ', opt_fun(1), ioptfun)

      ! convert between the Ripken and Edwards-Teng parametrization
      ! according to the formulas in "BETATRON MOTION WITH COUPLING OF HORIZONTAL AND VERTICAL DEGREES OF FREEDOM"
      ! from V. A. Lebedev    and  S. A. Bogacz

      deltaeValue = deltae ! equals 1.0 unless there is a cavity

      r11 = zero
      r12 = zero
      r21 = zero
      r22 = zero

      if (tw%beta(1,2)==zero .and. tw%beta(2,1)==zero) then

         ! in case there is absolutely no coupling kx and ky will be zero and u will be NaN
         ! and betx, bety, alfx, alfy will also evaluate as NaN if we apply the above formulae
         ! therefore we simply copy beta11 into betx and beta22 into bety in this case, so as
         ! to get the same values between twiss and ptc_twiss
         ! beta11, alfa11 etc... are multiplied by deltae before output
         ! hence we reflect this in the formula from Lebedev
         betx = tw%beta(1,1) * deltaeValue
         bety = tw%beta(2,2) * deltaeValue
         alfx = tw%alfa(1,1) * deltaeValue
         alfy = tw%alfa(2,2) * deltaeValue

      else

         kx=sqrt(tw%beta(1,2)/tw%beta(1,1)); ! multiplication by deltae in numerator and denominator
         ky=sqrt(tw%beta(2,1)/tw%beta(2,2));

         ! beta11, alfa11 etc... are multiplied by deltae before output
         ax=kx*tw%alfa(1,1) * deltaeValue -tw%alfa(1,2) * deltaeValue /kx;
         ! hence we reflect this in the formula from Lebedev
         ay=ky*tw%alfa(2,2) * deltaeValue -tw%alfa(2,1) * deltaeValue /ky;
         kxy2=kx*kx*ky*ky;
         
         
         if((abs(kx*kx-ky*ky).gt.TINY(ONE)).and.(abs(one-kxy2).gt.TINY(ONE))) then
            usqrt=kxy2*(one+(ax*ax-ay*ay)/(kx*kx-ky*ky)*(one-kxy2))
            if(usqrt.gt.TINY(ONE)) then
               usqrt=sqrt(usqrt)
!               u=(-kxy2+usqrt)/(one-kxy2)
               if(kxy2.le.usqrt) THEN
                  u=(-kxy2+usqrt)/(one-kxy2)
               else
                  u=(-kxy2-usqrt)/(one-kxy2)
               endif
            else
               u=-kxy2/(one-kxy2)
            endif

            ! betx, bety, alfx, alfy are the values computed by twiss with very good precision
            ! beta11, alfa11 etc... are multiplied by deltae before output
            ! hence we reflect this in the formula from Lebedev

            kappa=one-u

            betx = (tw%beta(1,1)/kappa) * deltaeValue
            bety = (tw%beta(2,2)/kappa) * deltaeValue
            alfx = (tw%alfa(1,1)/kappa) * deltaeValue
            alfy = (tw%alfa(2,2)/kappa) * deltaeValue

            bx = kx*kappa+u/kx
            by = ky*kappa-u/ky
            cx = kx*kappa-u/kx
            cy = ky*kappa+u/ky

            cosvp = (ax*ay-bx*cy)/(ay*ay+cy*cy)
            sinvp = (ax*cy+ay*bx)/(ay*ay+cy*cy)
            cosvm = (ax*ay+by*cx)/(ax*ax+cx*cx)
            sinvm = (ax*by-ay*cx)/(ax*ax+cx*cx)

            cosv2 =  sqrt((one+cosvp*cosvm-sinvp*sinvm)/two)
            sinv2 = -sqrt((one-cosvp*cosvm+sinvp*sinvm)/two)
            cosv1 = -sqrt((one+cosvp*cosvm+sinvp*sinvm)/two)
            sinv1 = -sqrt((one-cosvp*cosvm-sinvp*sinvm)/two)

            r11 = sqrt(tw%beta(2,2)/tw%beta(1,2))*(tw%alfa(1,2)*sinv2+u*cosv2)/kappa
            r12 = sqrt(tw%beta(1,1)*tw%beta(2,1))*sinv1/kappa
            r21 = (cosv2*(tw%alfa(1,2)*kappa-tw%alfa(2,2)*u)-sinv2*&
                 (u*kappa+tw%alfa(1,2)*tw%alfa(2,2)))/(kappa*sqrt(tw%beta(1,2)*tw%beta(2,2)))
            r22 = -sqrt(tw%beta(1,1)/tw%beta(2,1))*(u*cosv1+tw%alfa(2,1)*sinv1)/kappa
         else
            call fort_warn("ptc_twiss","Edwards-Teng beta function is set to regular beta because")
            call fort_warn("ptc_twiss","Argument of sqrt(kxy2*(1+(ax*ax-ay*ay)/(kx*kx-ky*ky)*(one-kxy2))) smaller than TINY")
            print*,"Argument of sqrt(kxy2*(1+(ax*ax-ay*ay)/(kx*kx-ky*ky)*(one-kxy2))) is: ",&
                 kxy2*(1+(ax*ax-ay*ay)/(kx*kx-ky*ky)*(one-kxy2))
                 
            print*, "kx=",kx
            print*, "ky=",ky
            print*, "kxy2=",kxy2
                 
            betx = tw%beta(1,1) * deltaeValue
            bety = tw%beta(2,2) * deltaeValue
            alfx = tw%alfa(1,1) * deltaeValue
            alfy = tw%alfa(2,2) * deltaeValue
         endif

      endif

      ! Edwards-Teng parameters go into betx, bety, alfx, alfy which are at the beginning of twiss_table_cols in madxl.h
      call double_to_table_curr(table_name, 'betx ', betx ) ! non contiguous with the above table entries
      call double_to_table_curr(table_name, 'bety ', bety ) ! hence we must store these values one by one
      call double_to_table_curr(table_name, 'alfx ', alfx )
      call double_to_table_curr(table_name, 'alfy ', alfy )
      call double_to_table_curr(table_name, 'r11 ', r11 )
      call double_to_table_curr(table_name, 'r12 ', r12 )
      call double_to_table_curr(table_name, 'r21 ', r21 )
      call double_to_table_curr(table_name, 'r22 ', r22 )

      call augment_count(table_name)

      !
      !--- track the Twiss functions' extremas
      call trackBetaExtrema(tw%beta,deltae,betx,bety,tw%disp)
      !---


    end subroutine puttwisstable
    !____________________________________________________________________________________________

    subroutine readrematrix
      !reads covariance matrix of the initial distribution
      implicit none

      re(1,1) = get_value('ptc_twiss ','re11 ')
      re(1,2) = get_value('ptc_twiss ','re12 ')
      re(1,3) = get_value('ptc_twiss ','re13 ')
      re(1,4) = get_value('ptc_twiss ','re14 ')
      re(1,5) = get_value('ptc_twiss ','re16 ')
      re(1,6) = get_value('ptc_twiss ','re15 ')
      re(2,1) = get_value('ptc_twiss ','re21 ')
      re(2,2) = get_value('ptc_twiss ','re22 ')
      re(2,3) = get_value('ptc_twiss ','re23 ')
      re(2,4) = get_value('ptc_twiss ','re24 ')
      re(2,5) = get_value('ptc_twiss ','re26 ')
      re(2,6) = get_value('ptc_twiss ','re25 ')
      re(3,1) = get_value('ptc_twiss ','re31 ')
      re(3,2) = get_value('ptc_twiss ','re32 ')
      re(3,3) = get_value('ptc_twiss ','re33 ')
      re(3,4) = get_value('ptc_twiss ','re34 ')
      re(3,5) = get_value('ptc_twiss ','re36 ')
      re(3,6) = get_value('ptc_twiss ','re35 ')
      re(4,1) = get_value('ptc_twiss ','re41 ')
      re(4,2) = get_value('ptc_twiss ','re42 ')
      re(4,3) = get_value('ptc_twiss ','re43 ')
      re(4,4) = get_value('ptc_twiss ','re44 ')
      re(4,5) = get_value('ptc_twiss ','re46 ')
      re(4,6) = get_value('ptc_twiss ','re45 ')
      re(5,1) = get_value('ptc_twiss ','re61 ')
      re(5,2) = get_value('ptc_twiss ','re62 ')
      re(5,3) = get_value('ptc_twiss ','re63 ')
      re(5,4) = get_value('ptc_twiss ','re64 ')
      re(5,5) = get_value('ptc_twiss ','re66 ')
      re(5,6) = get_value('ptc_twiss ','re65 ')
      re(6,1) = get_value('ptc_twiss ','re51 ')
      re(6,2) = get_value('ptc_twiss ','re52 ')
      re(6,3) = get_value('ptc_twiss ','re53 ')
      re(6,4) = get_value('ptc_twiss ','re54 ')
      re(6,5) = get_value('ptc_twiss ','re56 ')
      re(6,6) = get_value('ptc_twiss ','re55 ')

    end subroutine readrematrix
    !_________________________________________________________________

    subroutine readreforbit
      !reads covariance
      implicit none
      x(:)=zero
      x(1)=get_value('ptc_twiss ','x ')
      x(2)=get_value('ptc_twiss ','px ')
      x(3)=get_value('ptc_twiss ','y ')
      x(4)=get_value('ptc_twiss ','py ')
      x(5)=get_value('ptc_twiss ','pt ')
      x(6)=get_value('ptc_twiss ','t ')
    end subroutine readreforbit
    !_________________________________________________________________

    subroutine readinitialdistrib
      !reads covariance matrix of the initial distribution
      implicit none
      include 'madx_ptc_distrib.inc'
      type(taylor) ht
      type(pbfield) h
      type(damap) id
      type(normalform) norm
      real(dp) lam
      integer nd,nd_m
      logical fake_3
      integer jc(6)
      type(real_8) yy(6)
      integer :: dodo = 0
      real(dp) x(6)

      if(dodo==1) then
         x=zero
         call find_orbit(my_ring,x,1,default,c_1d_7)
         write(6,*) x
         call init(default,1,0,berz)
         call alloc(yy)
         call alloc(id)
         id=1
         yy=x+id
         call track(my_ring,yy,1,default)
         call print(yy,6)
         stop 999
      endif

      jc(1)=2
      jc(2)=1
      jc(3)=4
      jc(4)=3
      jc(5)=6
      jc(6)=5

      print*, "We are at initialization with moments of the distribution"

      call readrematrix()

      print*, re(1,:)
      print*, re(2,:)
      print*, re(3,:)
      print*, re(4,:)
      print*, re(5,:)
      print*, re(6,:)

      nd=2
      if(icase==6) nd=3
      nd_m=nd
      fake_3=.false.
      if (getdistrtype(3) /= distr_gauss.and.nd==3) then
         !here we have flat in delta
         nd_m=2
         fake_3=.true.
      endif

      call init(2,nd,0,0)

      call alloc(ht)
      call alloc(h)
      call alloc(id)
      call alloc(norm)

      do i = 1,nd_m*2
         do ii = 1,nd_m*2
            ht=ht+re(i,ii)*(-1)**(ii+i)*(1.0_dp.mono.jc(i))*(1.0_dp.mono.jc(ii))
         enddo
      enddo
      ht=-ht*pi

      lam=ten*full_abs(ht)
      ht=ht/lam
      if(fake_3) then !1959 is the yearh of birth of Etienne (just a number)
         ht=ht+(0.1959e0_dp.mono.'000020')+(0.1959e0_dp.mono.'000002')
      endif
      h=ht
      id=1
      id=texp(h,id)
      norm=id
      emi=zero
      emi(1:nd_m)=norm%tune(1:nd_m)*lam
      re=norm%a_t


      do i = 1,c_%nd2
         !call print(norm%a_t%v(i),6)
         unimap(i) = norm%a_t%v(i)
      enddo

      !      re=id
      call setemittances(emi(1),emi(2),emi(3))

      call kill(h)
      call kill(ht)
      call kill(id)
      call kill(norm)


    end subroutine readinitialdistrib
    !_________________________________________________________________

    subroutine readmatrixfromtable ! 26 april 2010: changed this routine
      ! because the format of the map_table has now changed to contain all terms,
      ! and not only the zeroth and first order ones...
      implicit none
      integer  :: double_from_table_row, table_length
      ! following added 26 april 2010
      integer :: order, row, nrows,i
      integer :: nx, nxp, ny, nyp, ndeltap, nt, index
      real(dp):: coeff
      !character(6) :: selector
      integer, dimension(6) :: jj ! 3 may 2010
      logical(lp) :: ignore_map_orbit
      real(dp),dimension(ndim2)::reval,aieval
      real(dp),dimension(ndim2,ndim2)::revec,aievec
      real(dp):: checkvalue
      
      
      order = get_value('ptc_twiss ','no ')

      row = 1 ! starts at one
      
      nrows = table_length("map_table ")

      do while(row .le. nrows) ! k=0 when read okay. k=-3 when the table has no row
         k = double_from_table_row("map_table ","coef ", row,doublenum)
         !write(0,*) 'k=',k
         !write(0,*) 'coef=',doublenum
         coeff=doublenum
         k = double_from_table_row("map_table ","n_vector ",row,doublenum)
         index = int(doublenum)
         k = double_from_table_row("map_table ","nx ",row,doublenum)
         nx = int(doublenum)
         !write(0,*) 'index=',index
         !write(0,*) 'nx=',nx
         k = double_from_table_row("map_table ","nxp ",row,doublenum)
         nxp = int(doublenum)
         !write(0,*) 'nxp=',nxp
         k = double_from_table_row("map_table ","ny ",row,doublenum)
         ny = int(doublenum)
         !write(0,*) 'ny=',ny
         k = double_from_table_row("map_table ","nyp ",row,doublenum)
         nyp = int(doublenum)
         !write(0,*) 'nyp=',nyp
         k = double_from_table_row("map_table ","ndeltap ",row,doublenum)
         ndeltap = int(doublenum)
         !write(0,*) 'ndeltap=',ndeltap
         k = double_from_table_row("map_table ","nt ",row,doublenum)
         nt = int(doublenum)
         !write(0,*) 'nt=',nt

         !      y(index)%T.sub.j=coeff ! are we able to do this? NO

         ! code proposed by piotr to achieve the equivalent of y(1).sub.'12345'=4.0
         !oldv = Y(1).sub.'12345'
         !newtoset = (4 - oldv).mono.'12345'
         !Y(1) = Y(1) + newtoset

         ! code proposed by etienne to achieve the same
         !call pok(y(1),j,4.d0)

         if (k.eq.0) then
            jj(1)=nx
            jj(2)=nxp
            jj(3)=ny
            jj(4)=nyp
            jj(5)=ndeltap
            jj(6)=nt
            call pok(y(index)%T,jj,coeff)
            ! the following gives the same result as the above
            !oldv = y(index).sub.jj
            !newtoset = (coeff - oldv).mono.jj ! mono for monomial
            !y(index)%t=y(index)%t+newtoset
            ! failed to compile the following work in one line
            !Y(index) = Y(index) + ((coeff-(Y(index).sub.jj)).mono.jj)
         endif

         row = row+1

      enddo

!      call daprint(y,28) ! to be compared with fort.18 created by ptc_normal
      
      ignore_map_orbit = get_value('ptc_twiss ','ignore_map_orbit ') .ne. 0
      
      if ( .not. ignore_map_orbit ) then
        do row=1,6
          x(row) = y(row).sub.'0'
        enddo  
      endif

      call maptoascript()

    end subroutine readmatrixfromtable

    !_________________________________________________________________
    subroutine readinitialmap ! from fort.18 file
      implicit none
      type(damap) :: map
      ! call readMapFromFort18(y)
      call alloc(map)
      call dainput(map,18)
      y = map
      call kill(map)
      call maptoascript()
      call reademittance()
    end subroutine readinitialmap
    !_________________________________________________________________

    !_________________________________________________________________

    subroutine readinitialmatrix
      !reads initial map elements from MAD-X ptc_twiss command parameters
      implicit none

      call readrematrix() !reads re
      call readreforbit() !reads x
      call initmapfrommatrix()
      call maptoascript()
      call reademittance()

    end subroutine readinitialmatrix
    !_________________________________________________________________

    subroutine readinitialascript
      !reads initial map elements from MAD-X ptc_twiss command parameters
      implicit none
      type(damap) :: map
      call alloc(map)
      call dainput(map,19)
      y = x+map
      call kill(map)
      call reademittance()

    end subroutine readinitialascript
    !_________________________________________________________________

    subroutine reademittance
      !initializes Y(6) from re(6,6)
      implicit none
      real(dp) :: emix,emiy,emiz

      emix = get_value('probe ','ex ')
      emiy = get_value('probe ','ey ')
      emiz = get_value('probe ','et ')

      call setemittances(emix,emiy,emiz)

    end subroutine reademittance
    !_________________________________________________________________

    subroutine initmapfrommatrix
      !initializes Y(6) from re(6,6)
      implicit none
      real(dp),dimension(ndim2)::reval,aieval
      real(dp),dimension(ndim2,ndim2)::revec,aievec
      real(dp):: checkvalue
      
      call liepeek(iia,icoast)
      if (getdebug() > 1) then
         write (6,'(8a8)')   "no","nv","nd","nd2","ndc","ndc2","ndt","ndpt"
         write (6,'(8i8)') iia(1),iia(2),iia(3),iia(4),icoast(1),icoast(2),icoast(3),icoast(4)
         print*, "c_%npara is ", c_%npara
      endif

      allocate(j(c_%nv))
      j(:)=0
      do i = 1,c_%npara
         do ii = 1,c_%npara
            j(ii)=1
            r=re(i,ii)-(y(i)%T.sub.j)
            y(i)%T=y(i)%T+(r.mono.j)
            j(ii)=0
         enddo
      enddo
      deallocate(j)

!      call daprint(y,29) ! to be compared with fort.18 created by ptc_normal and fort.28

      call eig6(re,reval,aieval,revec,aievec)
      do i=1,iia(4)-icoast(2)
         checkvalue = abs(reval(i)**2+aieval(i)**2 - one)
         if(checkvalue .gt.c_1d_10) then
            write(whymsg,*) "Provided matrix has eigenvalue more than 1e-10 off the unit circle ! plane = ",i, &
                            " r^2 = ", reval(i)**2+aieval(i)**2, " delta = ",checkvalue 
            call fort_warn("ptc_twiss",whymsg(:len_trim(whymsg)))
   
            if(checkvalue .gt.c_1d_8) then
            
              write(whymsg,*) "ERROR: Provided matrix has eigenvalue more than 1e-8 off the unit circle ! plane = ",i, &
                              " r^2 = ", reval(i)**2+aieval(i)**2, " delta = ",checkvalue 
              call seterrorflag(10,"ptc_twiss",whymsg);
              return;
            endif 
            
         endif
      enddo
    end subroutine initmapfrommatrix
    !_________________________________________________________________

    subroutine maptoascript
      !Performes normal form on a map, and plugs A_ in its place
      implicit none

      if (getdebug() > 2) then
         print*,"maptoascript: doing normal form"
      endif

      call alloc(normal)
      normal = y

      if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
         write(whymsg,*) 'DA got unstable during Normal Form: The closed solution does not exist. PTC msg: ', &
              messagelost(:len_trim(messagelost))
         call fort_warn('ptc_twiss::maptoascript: ',whymsg(:len_trim(whymsg)))
         if (icase == 6) then
           print*,""
           print*,"6D closed solution does not exist, you may try 4D or 5D (case = 4 or 5)"
           print*,"and if it works check setting of the cavities (LAG and VOLT)"
         endif

         call seterrorflag(10,"ptc_twiss::maptoascript ",whymsg);
         
         return
      endif

      if (getdebug() > 2) then
         print*,"maptoascript: normal form done"
      endif

      y = x + normal%a_t
      call kill(normal)

    end subroutine maptoascript
    !_________________________________________________________________

    subroutine readinitialtwiss(dt)
      !Reads initial twiss parameters from MAD-X command
      implicit none
      real(kind(1d0)) alpha(3),beta(3),disp(4),mu(3)
      type(real_8) al(3),be(3),di(4)
      type(pol_block_inicond) :: inicondknobs
      integer k_system
      real(dp)  sizept
      real(kind(1d0)) emiz
      real(dp) dt

      beta(1)  = get_value('ptc_twiss ','betx ')
      beta(2)  = get_value('ptc_twiss ','bety ')
      beta(3)  = get_value('ptc_twiss ','betz ')
      alpha(1) = get_value('ptc_twiss ','alfx ')
      alpha(2) = get_value('ptc_twiss ','alfy ')
      alpha(3) = get_value('ptc_twiss ','alfz ')
      disp(1)  = get_value('ptc_twiss ','dx ')
      disp(2)  = get_value('ptc_twiss ','dpx ')
      disp(3)  = get_value('ptc_twiss ','dy ')
      disp(4)  = get_value('ptc_twiss ','dpy ')
      mu(1)    = get_value('ptc_twiss ','mux ')
      mu(2)    = get_value('ptc_twiss ','muy ')
      mu(3)    = get_value('ptc_twiss ','muz ')

      if (getdebug() > 1) then
         print*,"TW: Ini betas ",beta
         print*,"TW: Ini alfas ",alpha
      endif

      if (c_%nd == 3) then
         if (beta(3) <= zero) then
            call fort_warn("ptc_twiss","Fatal Error: 6D requested and betz is smaller then or equal to 0!")
            call aafail("ptc_twiss","Fatal Error: 6D requested and betz is smaller then or equal to 0! program stops")
         endif

         beta(3) = (one+alpha(3)**2)/beta(3)
         alpha(3) =-alpha(3)

      endif


      x(:)=zero
      x(1)=get_value('ptc_twiss ','x ')
      x(2)=get_value('ptc_twiss ','px ')
      x(3)=get_value('ptc_twiss ','y ')
      x(4)=get_value('ptc_twiss ','py ')
      !      x(5)=get_value('ptc_twiss ','t ')
      !      x(6)=get_value('ptc_twiss ','pt ')
      x(6)=get_value('ptc_twiss ','t ')
      x(5)=get_value('ptc_twiss ','pt ')
      !frs plug deltap
      if(icase.eq.5 .or. icase.eq.56 ) x(5) = x(5) + dt


      if (getdebug() > 0 ) then
         CALL write_closed_orbit(icase,x)
      endif

      call reademittance()

      !Here we initialize Y(6)

      call alloc(be); call alloc(al); call alloc(di)

      !  code to power knows
      !


      do i=1,c_%nd
         !be(i)=beta(i)
         !al(i)=alpha(i)
         be(i)= beta(i)/(1_dp+x(5))
         al(i)=alpha(i)/(1_dp+x(5))
         !print*, 'be ', be(i), ' beta', beta(i)
      enddo

      do i=1,4
         di(i)=disp(i)
      enddo

      if (getnknobis() > 0) then
         !POWER KNOBS
         c_%knob = my_true

         k_system= (c_%npara - c_%nd2) + getnknobsm()
         inicondknobs = getknobinicond()

         do i=1,c_%nd

            if(inicondknobs%beta(i)/=0) then
               if (getdebug() > 1) then
                   print*,"Beta ",i," is knob no. ", inicondknobs%beta(i)
               endif
               call make_it_knob(be(i),k_system+inicondknobs%beta(i))
            endif

            if(inicondknobs%alfa(i)/=0) then
               if (getdebug() > 1) then
                   print*,"Alfa ",i," is knob no. ",  inicondknobs%alfa(i)
               endif
               call make_it_knob(al(i),k_system+inicondknobs%alfa(i))
            endif
         enddo

         do i=1,4
            if(inicondknobs%dispersion(i)/=0) then
               if (getdebug() > 1) then
                   print*,"Dispersion ",i," is knob no. ",  inicondknobs%dispersion(i)
               endif
               call make_it_knob(di(i),k_system+inicondknobs%dispersion(i))
            endif
         enddo

      endif


      y=x

      do i=1,c_%nd
                  !print*, " Beta(", i,")=", beta(i)
                  !call print(y(2*i-1),6)
                  !call print(y(2*i  ),6)

         y(2*i-1)= x(2*i-1) + sqrt(be(i)) * morph((one.mono.(2*i-1))    )
         y(2*i)= x(2*i) + one/sqrt(be(i)) * &
              (morph(  (one.mono.(2*i)) )-(al(i)) * morph((one.mono.(2*i-1))))

                  !call print(y(2*i-1),6)
                  !call print(y(2*i  ),6)
      enddo



      !--moments--!
      if( (c_%npara==5)   ) then
         
         !print*, "c_%npara ",c_%npara, " c_%ndpt ", c_%ndpt
         
         if ( (beta(3) .gt. zero) .and. (c_%ndpt/=0) ) then

            !Option one: sigma(5) is sqrt of emittance as in other two dimensions
            !          print*, "Init X5 with betaz ", beta(3)

            if (c_%nd < 3) then !otherwise it was already done, to be cleaned cause it is ugly and bug prone
               beta(3) = (one+alpha(3)**2)/beta(3)
               alpha(3) =-alpha(3)
            endif

            y(5) = x(5) +  sqrt( beta(3) )*morph((one.mono.5))
            y(6)= x(6) + one/sqrt(beta(3)) * (morph(  (one.mono.6) )-(alpha(3)) * morph(one.mono.5))

            emiz = get_value('probe ','et ')
            if ( emiz .le. 0  ) then
               sizept = get_value('probe ','sige ')
               emiz = sizept/sqrt(beta(3))
               !            print*, "Calculated Emittance ", emiz
            else
               emiz = sqrt(emiz)
               !            print*, "Read Emittance ", emiz
            endif

            call setsigma(5, emiz)
            call setsigma(6, emiz)

         else
            !by default we have no knowledge about longitudinal phase space, so init dp/p to ident
            !          print*, "Init X5 with ONE"
            !frs we need here the initial value of pt and t should not hurt
            y(5) = x(5) + morph((one.mono.5))
            y(6) = x(6) + morph((one.mono.5))
            call setsigma(5, get_value('probe ','sige '))
            call setsigma(6, get_value('probe ','sigt '))
         endif

      endif

      if ( (icase .gt. 5) .and. (get_value('probe ','et ') .le. 0) ) then !6 and 56
         !beta(3) is converted to gamma already (in 3rd coord the canonical planes are swapped)
         emiz = get_value('probe ','sige ')/sqrt(beta(3))
         print*, "icase=",icase," nd=",c_%nd, "sigma(5)=", emiz
         call setsigma(5, emiz)
         call setsigma(6, emiz)
      endif



      if(icase/=4) then
         do i=1,4
            y(i)= y(i) + di(i) * morph((one.mono.5))
         enddo
      endif


      if ( getdebug() > 2) then
         print*," Read the following BETA0 block in module ptc_twiss"
         print*," Twiss parameters:"
         write (6,'(6a8)')   "betx","alfx","bety","alfy","betz","alfz"
         write (6,'(6f8.4)')  beta(1),  alpha(1),  beta(2),  alpha(2),  beta(3),  alpha(3)
         write (6,'(4a8)')   "dx","dpx","dy","dpy"
         write (6,'(4f8.4)')  disp
         write (6,'(3a8)')   "mux","muy","muz"
         write (6,'(3f8.4)')  mu

         print*," Track:"
         write(6,'(6f8.4)') x
      endif


    end subroutine readinitialtwiss
    !____________________________________________________________________________________________

    subroutine onePassSummary(oneTurnMap,startorbit,suml)

      implicit none
      type(real_8),target :: oneTurnMap(6)
      real(dp),    target :: startorbit(6) ! six-dimensional phase-space state (usually referred-to as 'x')
      real(dp) :: suml ! cumulative length along the ring
      real(dp) :: rdp_mmilion ! float with zero (0)
      real(dp) :: deltap ! float with zero (0)
      
      rdp_mmilion= -1e6;
      
      call double_to_table_curr( summary_table_name, 'length ', suml ) ! total length of the machine

      call double_to_table_curr( summary_table_name, 'alpha_c ',    rdp_mmilion ) ! momemtum compaction factor
      call double_to_table_curr( summary_table_name, 'alpha_c_p ',  rdp_mmilion) ! derivative w.r.t delta-p/p
      call double_to_table_curr( summary_table_name, 'alpha_c_p2 ', rdp_mmilion) ! 2nd order derivative
      call double_to_table_curr( summary_table_name, 'alpha_c_p3 ', rdp_mmilion) ! 3rd order derivative
      call double_to_table_curr( summary_table_name, 'eta_c ',      rdp_mmilion) ! associated phase-slip factor
      call double_to_table_curr( summary_table_name, 'gamma_tr ',   rdp_mmilion) ! associated transition energy

      call double_to_table_curr( summary_table_name, 'q1 ', tw%mu(1))
      call double_to_table_curr( summary_table_name, 'q2 ', tw%mu(2))
      call double_to_table_curr( summary_table_name, 'dq1 ', rdp_mmilion)
      call double_to_table_curr( summary_table_name, 'dq2 ', rdp_mmilion)

      call double_to_table_curr( summary_table_name, 'qs ', rdp_mmilion)

      deltap = get_value('ptc_twiss ','deltap ')
      call double_to_table_curr( summary_table_name, 'deltap ', deltap)

      call putMinMaxRmses(summary_table_name,startorbit)

      call augment_count( summary_table_name ); ! only one row actually...


    end subroutine onePassSummary
    ! jluc
    ! compute momemtum-compaction factor in the same fashion it is carried-out in twiss.F

    subroutine oneTurnSummary(oneTurnMap,theAscript,startorbit,suml)

      implicit none
      type(real_8),target :: oneTurnMap(6),theAscript(6)
      real(dp),    target :: startorbit(6) 
      real(dp) :: suml ! cumulative length along the ring
      type(fibre), pointer :: fibrePtr
      real(dp) :: alpha_c, eta_c ! momentum-compaction factor & phase-slip factor
      real(dp) :: alpha_c_p ! first order derivative w.r.t delta-p/p
      real(dp) :: alpha_c_p2 ! second order derivative w.r.t delta-p/p
      real(dp) :: alpha_c_p3 ! third order derivative w.r.t delta-p/p
      real(dp) :: gamma_tr ! gamma_transition, or "transition energy" above which the particles' arrival time
      !  with respect to other particles is determined by its path length instead of by its velocity
      real(dp) :: deltap
      real(dp) :: betaRelativistic, gammaRelativistic
      real(dp) :: fractionalTunes(ndim)
      real(dp) :: chromaticities(2)
      integer :: i,j
      real(dp) :: sd ! as in twiss.F
      type(damap) :: yy ! added on November 6th to retreive momemtum compaction without differentiating the formula
      integer, dimension(6,6) :: coeffSelector = &
           reshape( (/1,0,0,0,0,0, &
           0,1,0,0,0,0, &
           0,0,1,0,0,0, &
           0,0,0,1,0,0, &
           0,0,0,0,1,0, &
           0,0,0,0,0,1/), (/6,6/))
      type(normalform) theNormalForm
      real(dp) :: dispersion(4)
      real(dp) :: dispersion_p(4) ! derivative of the dispersion w.r.t delta-p/p
      integer :: debugFiles
      integer :: icase
      integer :: order
      real(dp) :: rdp_mmilion 
      rdp_mmilion= -1e6;
    
      order = get_value('ptc_twiss ', 'no ')

      ! should end-up gracefully here in case the topology of the lattice is not those of a closed-ring

      debugFiles = 0 ! set it to one and fort.21, fort.22 and fort.23 are created


      ! 2. retreive the relativistic parameters beta and gamma
      ! (beta=v/c, gamma=E/mc^2 and gamma=1/sqrt(1-beta^2))
      betaRelativistic = get_value('probe ','beta ');
      gammaRelativistic = get_value('probe ','gamma ');
      icase = get_value('ptc_twiss ','icase ') ! mind the trailing space

      ! now retrieve the one-turn map's coefficients
      if (debugFiles .eq. 1) then
         do i=1,6
            do j=1,6
               write(21,*) "r(",i,j,")=",oneTurnMap(i).sub.coeffSelector(j,:)
            enddo
         enddo
      endif
     
      if (getdebug() > 1) then
        write(6,*) "Doing normal form ... "
      endif

      ! 4. retreive the dispersion coefficients
      ! (may be the coefficient of delta of the map?)
      ! decompose the map via a normal form to get the dispersion...
      call alloc(theNormalForm)
      theNormalForm = oneTurnMap

      if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
         write(whymsg,*) 'DA got unstable during Normal Form: The closed solution does not exist. PTC msg: ', &
              messagelost(:len_trim(messagelost))
         call fort_warn('ptc_twiss oneTurnSummary: ',whymsg(:len_trim(whymsg)))

         if (icase == 6) then
           print*,""
           print*,"6D closed solution does not exist, you may try 4D or 5D (case = 4 or 5)"
           print*,"and if it works check setting of the cavities (LAG and VOLT)"
         endif

         call seterrorflag(10,"ptc_twiss oneTurnSummary",whymsg)

         call kill(theNormalForm)
         return
      endif

      if (getdebug() > 1) then
        write(6,*) "Doing normal form ... Done"
      endif
      
      
      if (debugFiles .eq. 1) then
         call daprint(theNormalForm%A1,23) ! supposed to print dispersion's first and higher orders
         ! according to h_definition.f90: type normalform contains A1 as dispersion
         ! (would need to go through DHDJ to get the tune...)
      endif
      ! first order dispersions !?
      ! (at least checked that these values match those computed in twiss.F)
      dispersion(1) = theNormalForm%A1%v(1).sub.'000010'
      dispersion(2) = theNormalForm%A1%v(2).sub.'000010'
      dispersion(3) = theNormalForm%A1%v(3).sub.'000010'
      dispersion(4) = theNormalForm%A1%v(4).sub.'000010'

      if (debugFiles .eq. 1) then
         do i=1,4
            write(21,*) "dispersion(",i,")=", dispersion(i)
         enddo
      endif
       

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!                                       !!!!!!!!!!!!!!!!
      !!!!!!!!!!!! COMPACTION FACTOR, ITS HIGHER ORDERS  !!!!!!!!!!!!!!!!
      !!!!!!!!!!!! GAMMA TRANSITION, SLIP FACTOR         !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!                                       !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     

       
      ! in 4D we have no knowledge of longitudinal dynamics, it is not calculated
      ! the same 5D with time=false
      
      eta_c = rdp_mmilion
      alpha_c = rdp_mmilion
      gamma_tr = rdp_mmilion
      
      alpha_c_p = rdp_mmilion
      alpha_c_p2 = rdp_mmilion
      alpha_c_p3 = rdp_mmilion
      
      !call daprint(theNormalForm%dhdj,6) 

      
      if ( (icase.gt.4)  .and. (default%time) ) then
         
         ! Here R56 is dT/ddelta
         ! 5. apply formulas from twiss.F: 
         !sd = r(5,6)+r(5,1)*disp(1)+...+r(5,4)*disp(4)
         !print*,"ALPHA_C, GAMMA TR : TIME ON"

         sd = - oneTurnMap(6).sub.coeffSelector(5,:) ! 5/6 swap MADX/PTC
         !print*,'sd(0)',sd
         do i=1,4
            !print*, 'Disp',i,'=', dispersion(i)
            sd = sd - (oneTurnMap(6).sub.coeffSelector(i,:))*dispersion(i)
         enddo
         !print*,'sd(f)',sd

         eta_c = -sd * betaRelativistic**2 / suml 
         alpha_c = one / gammaRelativistic**2 + eta_c
         gamma_tr = one / sqrt(alpha_c) 

      elseif( (icase.eq.5)  .and. (default%time .eqv. .false.) ) then

         ! Here R56 is dL/ddelta 
         ! so we get alpha_c first from transfer matrix

         sd = +1.0*(oneTurnMap(6).sub.coeffSelector(5,:)) ! 5/6 swap MADX/PTC
         !print*,'sd(0)',sd
         do i=1,4
            !print*, 'Disp',i,'=', dispersion(i)
            sd = sd + (oneTurnMap(6).sub.coeffSelector(i,:))*dispersion(i)
         enddo
         !print*,'sd(f)',sd
         alpha_c = sd/suml
         eta_c = alpha_c - one / gammaRelativistic**2
         gamma_tr = one / sqrt(alpha_c)
         
      elseif( (icase.eq.56)  .and. (default%time .eqv. .false.) ) then

         !print*,"ALPHA_C, GAMMA TR : 56D TIME OFF"


         call alloc(yy)
         do i=1,c_%nd2 ! c_%nd2 is 6 when icase is 56 or 6 (but 4 when icase=5)
            yy%v(i) = oneTurnMap(i)%t
         enddo
         yy = theNormalForm%A1**(-1)*yy*theNormalForm%A1 ! takes away all dispersion dependency
         !write(0,*) 'for yy, c_%nd2 is ',c_%nd2 ! 0 is stderr
         alpha_c    = (yy%v(6).sub.'000010')/suml
         gamma_tr = one / sqrt(alpha_c)! overwrite the value obtained from the Twiss formula
         eta_c = alpha_c - one / gammaRelativistic**2

         alpha_c_p  = 2.0*(yy%v(6).sub.'000020')/suml

         if (order.ge.3) then
            alpha_c_p2 = 3.0*2.0*(yy%v(6).sub.'000030')/suml
         endif

         if (order.ge.4) then
            alpha_c_p3 = 4.0*3.0*2.0*(yy%v(6).sub.'000040')/suml
         endif


         call kill(yy)

      endif
      
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!! HIGHER ORDERS IN dP/P  !!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     if ( (icase.eq.5) .and. (default%time .eqv. .true.) ) then 

         !print*,"ALPHA_C dp/dp derivatives : 5D TIME ON"
         
         if (getdebug() > 2) then
           call kanalnummer(mf1)
           open(unit=mf1,file='oneTurnMap.5dt.txt')
           call daprint(oneTurnMap,mf1)
           close(mf1)
         endif
         
         ! compute delta-p/p dependency of alpha_c
         ! first order derivatives of the dispersions
         ! assuming icase=5
         !call daprint(oneTurnMap,88)

         ! always assume time=true, so that the fifth phase-space variable is deltap instead of pt!
         ! otherwise should issue a warning or an error!

         ! proceed slightly differently to retreive the dispersion's derivatives
         ! (so far only managed to get the dispersion on the normal form, but not its derivatives)
         
         ! note: should reuse information computed previously than redo this...
         ! if icase=5, Taylor series expansion disp = disp + delta_p
         do i=1,4
            ! apparently, we are up by a factor two
            dispersion_p(i) = 2.0*(theAscript(i)%t.sub.'000020') ! as usual in this file

         enddo

         ! compute derivative of the formula for eta_c

         sd = 0.0
         ! first dr65/ddeltap
         sd = sd - 1.0*(oneTurnMap(6).sub.'100010')*dispersion(1) &
              &  - 1.0*(oneTurnMap(6).sub.'010010')*dispersion(2) &
              &  - 1.0*(oneTurnMap(6).sub.'001010')*dispersion(3) &
              &  - 1.0*(oneTurnMap(6).sub.'000110')*dispersion(4) &
              &  - 2.0*(oneTurnMap(6).sub.'000020') &
	            ! new
              &  - 1.0*(oneTurnMap(6).sub.'000011')*(oneTurnMap(6).sub.'000001') &
	            ! then terms dr6i/ddeltap*dispersion(i)
              &  - 2.0*(oneTurnMap(6).sub.'200000')*dispersion(1)**2 &
              &  - (oneTurnMap(6).sub.'110000')*dispersion(1)*dispersion(2) &
              &  - (oneTurnMap(6).sub.'101000')*dispersion(1)*dispersion(3) &
              &  - (oneTurnMap(6).sub.'100100')*dispersion(1)*dispersion(4) &
	            ! new
              &  - (oneTurnMap(6).sub.'100010')*dispersion(1) &
	            ! dr62/ddeltap*disperion(2)
              &  - (oneTurnMap(6).sub.'110000')*dispersion(1)*dispersion(2) &
              &  - 2.0*(oneTurnMap(6).sub.'020000')*dispersion(2)**2 &
              &  - (oneTurnMap(6).sub.'011000')*dispersion(2)*dispersion(3) &
              &  - (oneTurnMap(6).sub.'010100')*dispersion(2)*dispersion(4) &
	            ! new
              &  - (oneTurnMap(6).sub.'010010')*dispersion(2) &
	            ! dr63/ddeltap*dispersion(3)
              &  - (oneTurnMap(6).sub.'101000')*dispersion(1)*dispersion(3) &
              &  - (oneTurnMap(6).sub.'011000')*dispersion(2)*dispersion(3) &
              &  - 2.0*(oneTurnMap(6).sub.'002000')*dispersion(3)**2 &
              &  - (oneTurnMap(6).sub.'001100')*dispersion(3)*dispersion(4) &
	            ! new
              &  - 1*(oneTurnMap(6).sub.'001010')*dispersion(3) &
	            ! dr64/ddeltap*dispersion(4)
              &  - (oneTurnMap(6).sub.'100100')*dispersion(1)*dispersion(4) &
              &  - (oneTurnMap(6).sub.'010100')*dispersion(2)*dispersion(4) &
              &  - (oneTurnMap(6).sub.'001100')*dispersion(3)*dispersion(4) &
              &  - 2.0*(oneTurnMap(6).sub.'000200')*dispersion(4)**2 &
	            ! new
              &  - 1*(oneTurnMap(6).sub.'000110')*dispersion(4)

         ! terms involving derivatives of the dispersions
         do i=1,4
!            print*, 'sd=',sd,' disp=',dispersion_p(i)
            sd = sd - (oneTurnMap(6).sub.coeffSelector(i,:))*dispersion_p(i)
         enddo

!         print*, 'sd=',sd
         alpha_c_p = sd * (betaRelativistic**2) / suml

         ! eventually, one could differentiate the above formula to obtain alpha_c_p2
         ! but for the time-being, expect icase to be 56 to compute alpha_c_p2 and alpha_c_p3.

         alpha_c_p2 = rdp_mmilion
         alpha_c_p3 = rdp_mmilion


             
     elseif ( (icase.eq.56) .and. (default%time .eqv. .true.) ) then ! here one may obtain the pathlength derivatives from the map

         !print*,"ALPHA_C dp/dp derivatives : 56D TIME ON"
         if (getdebug() > 2) then
           call kanalnummer(mf1)
           open(unit=mf1,file='oneTurnMap.56dt.txt')
           call daprint(oneTurnMap,mf1)
           close(mf1)
         endif
         
         call alloc(yy)
         do i=1,c_%nd2 ! c_%nd2 is 6 when icase is 56 or 6 (but 4 when icase=5)
            yy%v(i) = oneTurnMap(i)%t
         enddo
         yy = theNormalForm%A1**(-1)*yy*theNormalForm%A1 ! takes away all dispersion dependency

         alpha_c_p  = -2.0*(yy%v(6).sub.'000020')/suml
         
         if (order.ge.3) then
            alpha_c_p2 = -3.0*2.0*(yy%v(6).sub.'000030')/suml
         endif

         if (order.ge.4) then
            alpha_c_p3 = -4.0*3.0*2.0*(yy%v(6).sub.'000040')/suml
         endif

         call kill(yy)

      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!                                !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!    END OF COMPACTION FACTOR    !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!                                !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      ! also output the tune ...
      fractionalTunes = theNormalForm%tune
      ! the above is exactly equivalent to the following two lines (i.e. returns frac.tune)
      fractionalTunes(1) = theNormalForm%DHDJ%v(1).sub.'0000' ! as in So_fitting.f90
      fractionalTunes(2) = theNormalForm%DHDJ%v(2).sub.'0000' ! as in So_fitting.f90
      ! 26 november 2009
      if (icase.eq.6) then
         fractionalTunes(3) = - theNormalForm%DHDJ%v(3).sub.'0000' ! to enter here icase must be 6 but not only:
         ! there must be a cavity otherwise icase is set internally to 56 by my_state in ptc_module.f90
         !write(0,*) 'fractionTunes(3)=', fractionalTunes(3)
         ! in the above, inserted minus sign to match the 'phase' or 'tw%mu(3)'
         ! computed as atan2(ascript(6).sub.'000010',ascript(6).sub.'000001')/2*pi
      else
         fractionalTunes(3) = 0.0
         !write(0,*) 'nullify fractionalTunes(3), icase=',icase
      endif
      ! Q: is it possible to get the actual total tune, as returned by twiss.F?
      ! => no, not with a map...

      ! ... as well as the chromaticities
      if (icase.eq.5 .or. icase.eq.56) then
         chromaticities(1) = theNormalForm%DHDJ%v(1).sub.'00001' ! as in So_fitting.f90
         chromaticities(2) = theNormalForm%DHDJ%v(2).sub.'00001' ! as in So_fitting.f90
         ! to get chromaticities, went to higher order with above "call init_default(default,2,0)"
      else
         ! if icase = 6, delta_p is a phase-space variable and not an external parameter hence we can't compute chromaticies
         chromaticities(1) = 0.0
         chromaticities(2) = 0.0
      endif

      ! for debug: check the values by printing the map
      if (debugFiles .eq. 1) then
         call daprint(oneTurnMap,25) ! prints the one-turn map on file 25
         call daprint(theNormalForm%dhdj,26) ! print tunes, chromaticities and anharmonicities
         ! as done in madx_ptc_normal.f90
      endif

      call kill(theNormalForm)


      !write(6,*)
      !write(6,*) "Momentum compaction (momentum compaction factor & phase-slip factor):"
      !write(6,*) "alpha=", alpha_c, "eta=", eta_c
      !write(6,*) "fractional tunes=", fractionalTunes
      !write(6,*) "chromaticities=", chromaticities

      deltap = get_value('ptc_twiss ','deltap ')

      ! write the data into the ptc_twiss_summary table
      ! warning: must preserve the order EXACTLY
      call double_to_table_curr( summary_table_name, 'length ', suml ) ! total length of the machine
      call double_to_table_curr( summary_table_name, 'alpha_c ', alpha_c ) ! momemtum compaction factor
      call double_to_table_curr( summary_table_name, 'alpha_c_p ', alpha_c_p) ! derivative w.r.t delta-p/p
      call double_to_table_curr( summary_table_name, 'alpha_c_p2 ', alpha_c_p2) ! 2nd order derivative
      call double_to_table_curr( summary_table_name, 'alpha_c_p3 ', alpha_c_p3) ! 3rd order derivative
      call double_to_table_curr( summary_table_name, 'eta_c ', eta_c ) ! associated phase-slip factor
      call double_to_table_curr( summary_table_name, 'gamma_tr ', gamma_tr) ! associated transition energy
      call double_to_table_curr( summary_table_name, 'q1 ', fractionalTunes(1))
      call double_to_table_curr( summary_table_name, 'q2 ', fractionalTunes(2))
      call double_to_table_curr( summary_table_name, 'dq1 ', chromaticities(1))
      call double_to_table_curr( summary_table_name, 'dq2 ', chromaticities(2))
      ! 26 november 2009
      call double_to_table_curr( summary_table_name, 'qs ', fractionalTunes(3))
      ! write the extremas of the Twiss functions
      
      
      
      call double_to_table_curr( summary_table_name, 'deltap ', deltap)

      call putMinMaxRmses(summary_table_name, startorbit)


      call augment_count( summary_table_name ); ! only one row actually...

    end subroutine oneTurnSummary


  END subroutine ptc_twiss
  !____________________________________________________________________________________________
  
  
  subroutine putMinMaxRmses(summary_table_name,startorbit)
    implicit none
    character(48) :: summary_table_name
    real(dp) :: startorbit(6)
    real(dp) :: xrms(6)
    
      call double_to_table_curr( summary_table_name, 'beta_x_min ', minBetX)
      call double_to_table_curr( summary_table_name, 'beta_x_max ', maxBetX)
      call double_to_table_curr( summary_table_name, 'beta_y_min ', minBetY)
      call double_to_table_curr( summary_table_name, 'beta_y_max ', maxBetY)

      call double_to_table_curr( summary_table_name, 'beta_x_min ', minBetX)
      call double_to_table_curr( summary_table_name, 'beta_x_max ', maxBetX)
      call double_to_table_curr( summary_table_name, 'beta_y_min ', minBetY)
      call double_to_table_curr( summary_table_name, 'beta_y_max ', maxBetY)

      call double_to_table_curr( summary_table_name, 'beta11min ', minBeta(1,1))
      call double_to_table_curr( summary_table_name, 'beta12min ', minBeta(1,2))
      call double_to_table_curr( summary_table_name, 'beta13min ', minBeta(1,3))
      call double_to_table_curr( summary_table_name, 'beta21min ', minBeta(2,1))
      call double_to_table_curr( summary_table_name, 'beta22min ', minBeta(2,2))
      call double_to_table_curr( summary_table_name, 'beta23min ', minBeta(2,3))
      call double_to_table_curr( summary_table_name, 'beta31min ', minBeta(3,1))
      call double_to_table_curr( summary_table_name, 'beta32min ', minBeta(3,2))
      call double_to_table_curr( summary_table_name, 'beta33min ', minBeta(3,3))

      call double_to_table_curr( summary_table_name, 'beta11max ', maxBeta(1,1))
      call double_to_table_curr( summary_table_name, 'beta12max ', maxBeta(1,2))
      call double_to_table_curr( summary_table_name, 'beta13max ', maxBeta(1,3))
      call double_to_table_curr( summary_table_name, 'beta21max ', maxBeta(2,1))
      call double_to_table_curr( summary_table_name, 'beta22max ', maxBeta(2,2))
      call double_to_table_curr( summary_table_name, 'beta23max ', maxBeta(2,3))
      call double_to_table_curr( summary_table_name, 'beta31max ', maxBeta(3,1))
      call double_to_table_curr( summary_table_name, 'beta32max ', maxBeta(3,2))
      call double_to_table_curr( summary_table_name, 'beta33max ', maxBeta(3,3))


      call double_to_table_curr( summary_table_name, 'disp1min ', minDisp(1))
      call double_to_table_curr( summary_table_name, 'disp2min ', minDisp(2))
      call double_to_table_curr( summary_table_name, 'disp3min ', minDisp(3))
      call double_to_table_curr( summary_table_name, 'disp4min ', minDisp(4))

      call double_to_table_curr( summary_table_name, 'disp1max ', maxDisp(1))
      call double_to_table_curr( summary_table_name, 'disp2max ', maxDisp(2))
      call double_to_table_curr( summary_table_name, 'disp3max ', maxDisp(3))
      call double_to_table_curr( summary_table_name, 'disp4max ', maxDisp(4))

      !!Start orbit
      call double_to_table_curr( summary_table_name,'orbit_x ',  startorbit(1))
      call double_to_table_curr( summary_table_name,'orbit_px ', startorbit(2))
      call double_to_table_curr( summary_table_name,'orbit_y ',  startorbit(3))
      call double_to_table_curr( summary_table_name,'orbit_py ', startorbit(4))
      call double_to_table_curr( summary_table_name,'orbit_pt ', startorbit(5))
      call double_to_table_curr( summary_table_name,'orbit_-cT ',startorbit(6))

      xrms = sqrt(sum2Orbit / nobsOrbit)
      
       call double_to_table_curr(summary_table_name,'xcorms ', xrms(1))
       call double_to_table_curr(summary_table_name,'pxcorms ',xrms(2))
       call double_to_table_curr(summary_table_name,'ycorms ', xrms(3))
       call double_to_table_curr(summary_table_name,'pycorms ',xrms(4))
       call double_to_table_curr(summary_table_name,'ptcorms ',xrms(5))
       call double_to_table_curr(summary_table_name,'tcorms ', xrms(6))


       call double_to_table_curr(summary_table_name,'xcomin ' ,minOrbit(1))
       call double_to_table_curr(summary_table_name,'pxcomin ',minOrbit(2))
       call double_to_table_curr(summary_table_name,'ycomin ' ,minOrbit(3))
       call double_to_table_curr(summary_table_name,'pycomin ',minOrbit(4))
       call double_to_table_curr(summary_table_name,'ptcomin ',minOrbit(5))
       call double_to_table_curr(summary_table_name,'tcomin ' ,minOrbit(6))

       call double_to_table_curr(summary_table_name,'xcomax ' ,maxOrbit(1))
       call double_to_table_curr(summary_table_name,'pxcomax ',maxOrbit(2))
       call double_to_table_curr(summary_table_name,'ycomax ' ,maxOrbit(3))
       call double_to_table_curr(summary_table_name,'pycomax ',maxOrbit(4))
       call double_to_table_curr(summary_table_name,'ptcomax ',maxOrbit(5))
       call double_to_table_curr(summary_table_name,'tcomax ' ,maxOrbit(6))

     
      

  end subroutine putMinMaxRmses

  subroutine computeDeltapDependency(y,s1)
    implicit none
    type(real_8), intent(in)  :: y(6)
    type(twiss),  intent(inout)  :: s1
    integer :: k,i
    integer :: J(lnv) ! the map's coefficient selector, as usual
    integer :: Jderiv(lnv) ! to store the map's coefficient selector of the derivative w.r.t deltap
    real(kind(1d0)) :: get_value ! C-function
    integer :: no ! order must be at equal to 2 to be able to get terms of the form x*deltap
    ! required to evaluate the derivatives of Twiss parameters w.r.t deltap
    integer :: ndel ! as in subroutine 'equaltwiss'...
    real(dp) :: rdp_mmilion ! float with zero (0)

    rdp_mmilion= -1e6;

    ! in order to avoid this message, should prevent entering this subroutine
    ! in case none of the Twiss derivatives is selected...

    no = get_value('ptc_twiss ','no ')
    if ( no .lt. 2 ) then
       call fort_warn('madx_ptc_twiss.f90 <ptc_twiss>:','Order in computeDeltapDependency() is smaller then 2')
       print*, "Order is ", no
       return
    endif

    J=0
    do i=1,c_%nd ! c_%nd is global variable
       do k=1,c_%nd ! the same
          J(2*i-1)=1
          Jderiv = J ! vector copy
          Jderiv(5)=1 ! the delta_p coefficient
          s1%beta_p(k,i)= (Y(2*k-1).sub.Jderiv)*(Y(2*k-1).sub.J) + (Y(2*k-1).sub.J)*(Y(2*k-1).sub.Jderiv)
          s1%alfa_p(k,i)= -(Y(2*k-1).sub.Jderiv)*(Y(2*k).sub.J) - (Y(2*k-1).sub.J)*(Y(2*k).sub.Jderiv)
          s1%gama_p(k,i)= (Y(2*k).sub.Jderiv)*(Y(2*k).sub.J) + (Y(2*k).sub.J)*(Y(2*k).sub.Jderiv)
          J(2*i-1)=0
          J(2*i)=1
          Jderiv = J ! vector copy
          Jderiv(5)=1 ! the delta_p coefficient
          s1%beta_p(k,i) = s1%beta_p(k,i) &
               + (Y(2*k-1).sub.Jderiv)*(Y(2*k-1).sub.J) + (Y(2*k-1).sub.J)*(Y(2*k-1).sub.Jderiv)
          s1%alfa_p(k,i)= s1%alfa_p(k,i) &
               - (Y(2*k-1).sub.Jderiv)*(Y(2*k).sub.J) - (Y(2*k-1).sub.J)*(Y(2*k).sub.Jderiv)
          s1%gama_p(k,i)= s1%gama_p(k,i) &
               + (Y(2*k).sub.Jderiv)*(Y(2*k).sub.J)+(Y(2*k).sub.J)*(Y(2*k).sub.Jderiv)
          J(2*i)=0
       enddo
    enddo

    ! the computations above match the following formulas, obtained by derivation of the Twiss parameters using the chain-rule
    ! beta derivatives w.r.t delta_p
    !    beta11 = two * (y(1)%t.sub.'100000')*(y(1)%t.sub.'100010') + two * (y(1)%t.sub.'010000')*(y(1)%t.sub.'010010')
    !    beta12 = two * (y(1)%t.sub.'001000')*(y(1)%t.sub.'001010') + two * (y(1)%t.sub.'000100')*(y(1)%t.sub.'000110')
    !    beta21 = two * (y(3)%t.sub.'100000')*(y(3)%t.sub.'100010') + two * (y(3)%t.sub.'010000')*(y(3)%t.sub.'010010')
    !    beta22 = two * (y(3)%t.sub.'001000')*(y(3)%t.sub.'001010') + two * (y(3)%t.sub.'000100')*(y(3)%t.sub.'000110')
    ! alpha derivatives w.r.t delta_p
    !    alfa11 = -((y(1)%t.sub.'100010')*(y(2)%t.sub.'100000')+(y(1)%t.sub.'100000')*(y(2).sub.'100010')+&
    !         (y(1)%t.sub.'010010')*(y(2)%t.sub.'010000')+(y(1)%t.sub.'010000')*(y(2)%t.sub.'010010'))
    !    alfa12 = -((y(1)%t.sub.'001010')*(y(2)%t.sub.'001000')+(y(1)%t.sub.'001000')*(y(2)%t.sub.'001010')+&
    !         (y(1)%t.sub.'000110')*(y(2)%t.sub.'000100')+(y(1)%t.sub.'000100')*(y(2)%t.sub.'000110'))
    !    alfa21 = -((y(3)%t.sub.'100010')*(y(4)%t.sub.'100000')+(y(3)%t.sub.'100000')*(y(4)%t.sub.'100010')+&
    !         (y(3)%t.sub.'010010')*(y(4)%t.sub.'010000')+(y(3)%t.sub.'010000')*(y(4)%t.sub.'010010'))
    !    alfa22 = -((y(3)%t.sub.'001010')*(y(4)%t.sub.'001000')+(y(3)%t.sub.'001000')*(y(4)%t.sub.'001010')+&
    !         (y(3)%t.sub.'000110')*(y(4)%t.sub.'000100')+(y(3)%t.sub.'000100')*(y(4)%t.sub.'000110'))
    ! gamma derivatives w.r.t delta_p
    !    gama11 = (y(2)%t.sub.'100010')*(y(2)%t.sub.'100000')+(y(2)%t.sub.'100000')*(y(2)%t.sub.'100010')+&
    !         (y(2)%t.sub.'010010')*(y(2)%t.sub.'010000')+(y(2)%t.sub.'010000')*(y(2)%t.sub.'010010')
    !    gama12 = (y(2)%t.sub.'001010')*(y(2)%t.sub.'001000')+(y(2)%t.sub.'001000')*(y(2)%t.sub.'001010')+&
    !         (y(2)%t.sub.'000110')*(y(2)%t.sub.'000100')+(y(2)%t.sub.'000100')*(y(2)%t.sub.'000110')
    !    gama21 = (y(4)%t.sub.'100010')*(y(4)%t.sub.'100000')+(y(4)%t.sub.'100000')*(y(4)%t.sub.'100010')+&
    !         (y(4)%t.sub.'010010')*(y(4)%t.sub.'010000')+(y(4)%t.sub.'010000')*(y(4)%t.sub.'010010')
    !    gama22 = (y(4)%t.sub.'001010')*(y(4)%t.sub.'001000')+(y(4)%t.sub.'001000')*(y(4)%t.sub.'001010')+&
    !         (y(4)%t.sub.'000110')*(y(4)%t.sub.'000100')+(y(4)%t.sub.'000100')*(y(4)%t.sub.'000110')

    ! now compute deltap dependencies of the dispersion

    ! code to be differentiated, as in subroutine 'equaltwiss'
    !    J=0
    !    !here ND2=4 and delta is present      nd2=6 and delta is a constant
    !    !      print*,"nv",c_%nv,"nd2",c_%nd2,"np",c_%np,"ndpt",c_%ndpt ,"=>",c_%nv-c_%nd2-c_%np
    !    if( (c_%npara==5)       .or.  (c_%ndpt/=0) ) then
    !       !when there is no cavity it gives us dispersions
    !       do i=1,4
    !          lat(0,i,1)=(Y(i)%t.sub.J5)
    !       enddo
    !    elseif (c_%nd2 == 6) then
    !       do i=1,4
    !          lat(0,i,1) =              (Y(i)%t.sub.J5)*(Y(6)%t.sub.J6)
    !          lat(0,i,1) = lat(0,i,1) + (Y(i)%t.sub.J6)*(Y(5)%t.sub.J5)
    !       enddo
    !    else
    !       do i=1,4
    !          lat(0,i,1)=zero
    !       enddo
    !    endif
    !
    !    ...
    !
    !    !when there is no cavity it gives us dispersions
    !    do i=1,c_%nd2-2*ndel
    !       s1%disp(i)=lat(0,i,1)
    !    enddo



    ! differentiating the code that computes the dispersion...
    ndel=0
    ! as in subroutine 'equaltwiss'...
    if (c_%ndpt/=0) then
       ndel=1
    endif
    if ((c_%npara ==5) .or. (c_%ndpt /= 0)) then
       do i=1,c_%nd2-2*ndel ! should it be 1 to 4?
          ! in 4D+deltap case, coefficients are those of the Taylor series development
          ! of deltap (factorial)
          ! disp = disp0 + sigma (1/n!) disp_pn
          if (no.gt.2) then
             s1%disp_p(i) = 2.0*(y(i)%t.sub.'000020')
          else
             s1%disp_p(i) = rdp_mmilion
          endif
          if (no.gt.2) then
             s1%disp_p2(i) = 3.0*2.0*(y(i)%t.sub.'000030') ! assume at least order 3 for the map
          else
             s1%disp_p2(i) = rdp_mmilion
          endif
             !call fort_warn("ptc_twiss ","assume no>=3 for dispersion's 2nd order derivatives w.r.t delta-p")
          !endif
          if (no.gt.3) then
             s1%disp_p3(i) = 4.0*3.0*2.0*(y(i)%t.sub.'000040')
          else
          !   call fort_warn("ptc_twiss ","assume no>=4 for dispersion's 3rd order derivatives w.r.t delta-p")
             s1%disp_p3(i) = rdp_mmilion
          endif
       enddo
    elseif (c_%nd2 ==6) then
       do i=1,c_%nd2-2*ndel ! should it be 1 to 4?
          ! finally never enter here: these differentiation formulas turned-out not to apply to the case
          ! where deltap is a phase-space state-variable rather than an externally supplied
          ! parameter.
          ! u'v+uv'
          s1%disp_p(i) = &
               2.0*(y(i)%t.sub.'000020')*(y(6)%t.sub.'000001')+&
               (y(i)%t.sub.'000010')*(y(6)%t.sub.'000011')+&
               (y(i)%t.sub.'000011')*(y(5)%t.sub.'000010')+&
               (y(i)%t.sub.'000001')*2.0*(y(5)%t.sub.'000020')
       enddo
    else
       do i=1,c_%nd2-2*ndel ! should it be 1 to 4?
          s1%disp_p(i) = 0
       enddo
    endif
    ! write(6,*) "dispersion derivative (1) = ",s1%disp_p(1)
    ! checked the above yields non-zero value...

  end subroutine computeDeltapDependency
  !____________________________________________________________________________________________

  ! --- a set of routines to track extremas of Twiss functions
  subroutine resetBetaExtremas()
      resetBetaExtrema = .true.
  end subroutine resetBetaExtremas
  
  subroutine trackBetaExtrema(betas,den, betx, bety, disp)
    implicit none
    real(dp) :: betas(3,3)
    real(dp) :: den   ! delta energy from acceleration 
    real(dp) :: betx
    real(dp) :: bety
    real(dp) :: disp(4)
    integer  :: i,j ! iterators
    real(dp) :: value ! auxiliary variable
    
    if (resetBetaExtrema) then
    
       resetBetaExtrema = .false.
       do i=1,3
         do j=1,3
           minBeta(i,j) = betas(i,j)*den
           maxBeta(i,j) = betas(i,j)*den
         enddo
       enddo
       minBetX = betx
       maxBetX = betx
       minBetY = bety
       maxBetY = bety
       do i=1,4
         maxDisp(i) = disp(i)
         maxDisp(i) = disp(i)
       enddo
       return
    endif
    

    do i=1,3
      do j=1,3
        value = betas(i,j)*den
        if (minBeta(i,j) .gt. value) then
           minBeta(i,j) = value
        elseif (maxBeta(i,j) .lt. value) then
           maxBeta(i,j) = value
        endif
      enddo
    enddo
    
    if (minBetX .gt. betx) minBetX = betx
    if (maxBetX .lt. betx) maxBetX = betx
    if (minBetY .gt. bety) minBetY = bety
    if (maxBetY .lt. bety) maxBetY = bety
    
    do i=1,4
      if (minDisp(i) .gt. disp(i)) then
         minDisp(i) = disp(i)
      elseif (maxDisp(i) .lt. disp(i)) then
         maxDisp(i) = disp(i)
      endif
    enddo

  end subroutine trackBetaExtrema
  ! --- end of set of routines

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  
  subroutine trackOrbitExtremaAndRms(orbit)
    implicit none
    real(dp) :: orbit(6)
    integer  :: i
    
    
    if (resetOrbitExtrema) then
      
      resetOrbitExtrema = .false.
      nobsOrbit = 1

      do i=1,6
        minOrbit(i)  = orbit(i)
        maxOrbit(i)  = orbit(i)
        sum2Orbit(i) = orbit(i)*orbit(i)
      enddo
      
      return      
    endif

    nobsOrbit = nobsOrbit + 1

    do i=1,6
      if ( orbit(i) .lt. minOrbit(i) ) minOrbit(i)  = orbit(i)
      if ( orbit(i) .gt. maxOrbit(i) ) maxOrbit(i)  = orbit(i)
      sum2Orbit(i) = sum2Orbit(i) + orbit(i)*orbit(i)
    enddo

  end subroutine trackOrbitExtremaAndRms

  !____________________________________________________________________________________________


  subroutine getBeamBeam()
   implicit none
   integer                 :: i,e,elcode
   integer, external       :: restart_sequ, & !  restart beamline and return number of beamline node
                              advance_node    !  advance to the next node in expanded sequence
                                              !  =0 (end of range), =1 (else)
   REAL(KIND(1d0)), external :: node_value  !/*returns value for parameter par of current element */
   type(fibre), pointer    :: p
   real (dp)               :: fk
   
   TYPE(INTEGRATION_NODE),POINTER :: CURR_SLICE
   
   e=restart_sequ()
   p=>my_ring%start
   do e=1, my_ring%n !in slices e goes to nelem + 1 because the last slice is the fist one.

     elcode=node_value('mad8_type ')

     if (elcode .eq. 22) then

        if (getdebug() > 1 ) then
          write(6,*) " Beam-Beam position at element named >>",p%mag%name,"<<"
        endif

        CURR_SLICE => p%t1

        do while (.not. (CURR_SLICE%cas==case0.or.CURR_SLICE%cas==caset) )
          if (associated(CURR_SLICE,p%t2)) exit
          CURR_SLICE => CURR_SLICE%next
          !print*, CURR_SLICE%cas
        enddo

        !print *,  'BB Node Case NO: ',CURR_SLICE%cas

        if(((CURR_SLICE%cas==case0).or.(CURR_SLICE%cas==caset))) then !must be 0 or 3

          if(.not.associated(CURR_SLICE%BB)) call alloc(CURR_SLICE%BB)

          call getfk(fk)
          CURR_SLICE%bb%fk = fk
          CURR_SLICE%bb%sx = node_value('sigx ')
          CURR_SLICE%bb%sy = node_value('sigy ')
          CURR_SLICE%bb%xm = node_value('xma ')
          CURR_SLICE%bb%ym = node_value('yma ')
          CURR_SLICE%bb%PATCH=.true.
          if (getdebug() > 2 ) then
            print*, "BB fk=",CURR_SLICE%bb%fk
            print*, "BB sx=",CURR_SLICE%bb%sx
            print*, "BB sy=",CURR_SLICE%bb%sy
            print*, "BB xm=",CURR_SLICE%bb%xm
            print*, "BB ym=",CURR_SLICE%bb%ym
          endif

          do_beam_beam = .true.

        else
          call fort_warn('getBeamBeam: ','Bad node case for BeamBeam')
        endif 

      endif

      i=advance_node() ! c-code go to the next node -> the passed value is never used, just to shut up a compiler
      p=>p%next
    enddo

  end subroutine getBeamBeam


end module madx_ptc_twiss_module
