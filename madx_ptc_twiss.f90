module madx_ptc_twiss_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
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
  private

  !============================================================================================
  !  PUBLIC INTERFACE
  public                         :: twiss,ptc_twiss



  !============================================================================================
  !  PRIVATE
  !    data structures

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
  end interface

  interface alloc
     module procedure alloctwiss
  end interface

  interface kill
     module procedure killtwiss
  end interface

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

  logical :: slice_magnets, deltap_dependency, momentumCompactionToggle
 
  real(dp)                :: minBeta(3,3) ! jluc: to store extremas of Twiss functions (show-up in summary table
  real(dp)                :: maxBeta(3,3) ! jluc: to store extremas of Twiss functions (show-up in summary table) 

  !============================================================================================
  !  PRIVATE
  !    routines
  private zerotwiss,equaltwiss,killtwiss



contains
  !____________________________________________________________________________________________

  subroutine equaltwiss(s1,Y)
    implicit none
    type(twiss), intent(inout)::s1
    type(real_8), intent(in)::Y(ndd)
    integer jj,i,k, ndel, n
    real(dp) :: lat(0:6,6,3)
    real(dp) :: test, dph
    real(dp) :: epsil=1e-12  !
    integer  :: J(lnv)

    lat = zero

    n=3  ! 1 2 3 are tunes

    ndel=0
    if(c_%ndpt/=0)  then
       !       print*, "We are at the mode 6D + nocav"
       ndel=1  !this is 6D without cavity (MADX icase=56)
    endif

    J=0;
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


! in 4D we get following formulas for Twiss parameters. They can be derived to compute deltap dependency
! beta[ 1 1 ]=(Y(1)%t.sub.100000)*(Y(1)%t.sub.100000)+(Y(1)%t.sub.010000)*(Y(1)%t.sub.010000)
! alpha[ 1 1 ]=-(Y(1)%t.sub.100000)*(Y(2)%t.sub.100000)+(Y(1)%t.sub.010000)*(Y(2)%t.sub.010000)
! gamma[ 1 1 ]=(Y(2)%t.sub.100000)*(Y(2)%t.sub.100000)+(Y(2)%t.sub.010000)*(Y(2)%t.sub.010000)
! beta[ 2 1 ]=(Y(3)%t.sub.100000)*(Y(3)%t.sub.100000)+(Y(3)%t.sub.010000)*(Y(3)%t.sub.010000)
! alpha[ 2 1 ]=-(Y(3)%t.sub.100000)*(Y(4)%t.sub.100000)+(Y(3)%t.sub.010000)*(Y(4)%t.sub.010000)
! gamma[ 2 1 ]=(Y(4)%t.sub.100000)*(Y(4)%t.sub.100000)+(Y(4)%t.sub.010000)*(Y(4)%t.sub.010000)
! beta[ 1 2 ]=(Y(1)%t.sub.001000)*(Y(1)%t.sub.001000)+(Y(1)%t.sub.000100)*(Y(1)%t.sub.000100)
! alpha[ 1 2 ]=-(Y(1)%t.sub.001000)*(Y(2)%t.sub.001000)+(Y(1)%t.sub.000100)*(Y(2)%t.sub.000100)
! gamma[ 1 2 ]=(Y(2)%t.sub.001000)*(Y(2)%t.sub.001000)+(Y(2)%t.sub.000100)*(Y(2)%t.sub.000100)
! beta[ 2 2 ]=(Y(3)%t.sub.001000)*(Y(3)%t.sub.001000)+(Y(3)%t.sub.000100)*(Y(3)%t.sub.000100)
! alpha[ 2 2 ]=-(Y(3)%t.sub.001000)*(Y(4)%t.sub.001000)+(Y(3)%t.sub.000100)*(Y(4)%t.sub.000100)
! gamma[ 2 2 ]=(Y(4)%t.sub.001000)*(Y(4)%t.sub.001000)+(Y(4)%t.sub.000100)*(Y(4)%t.sub.000100)

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

    !--- track the Twiss functions' extremas
    call trackBetaExtrema(1,1,s1%beta(1,1))
    call trackBetaExtrema(2,2,s1%beta(2,2))
    !---
   
    s1%mu=phase

    do k=1,c_%nd
       do i=1,c_%nd2
          s1%eigen(k*2-1,i) = Y(k*2-1).sub.fo(i,:)
          s1%eigen(k*2  ,i) = Y(k*2  ).sub.fo(i,:)
       enddo
    enddo


    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
       call fort_warn('ptc_twiss: ','DA in twiss got unstable')
       call seterrorflag(10,"ptc_twiss ","DA got unstable in twiss ");
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
    implicit none
    include 'twissa.fi'
    logical(lp)             :: closed_orbit,beta_flg
    integer                 :: k,i,ii
    integer                 :: no,mynd2,npara,nda,icase,flag_index,why(9),my_nv,nv_min
    character(200)          :: whymsg
    integer                 :: ioptfun,iii,restart_sequ,advance_node,mf1,mf2
    integer                 :: tab_name(*)
    integer                 :: summary_tab_name(*)
    real(dp)                :: x(6)
    real(dp)                :: deltap0,deltap,d_val
    real(kind(1d0))         :: get_value,suml
    integer                 :: geterrorflag !C function that returns errorflag value
    type(real_8)            :: y(6), transfermap(6)
    type(twiss)             :: tw
    type(fibre), POINTER    :: current
    type(work)              :: startfen !Fibre energy at the start
    real(dp)                :: r,re(6,6),dt
    logical(lp)             :: initial_matrix_manual, initial_matrix_table, initial_map_manual
    logical(lp)             :: initial_distrib_manual, initial_ascript_manual
    integer                 :: row, rmatrix
    real(dp)                :: emi(3)
    logical(lp)             :: skipnormalform, tracktm
    character*48            :: summary_table_name
    
    skipnormalform = my_false

    !all zeroing
    testold = zero
    phase = zero
    do i=1,6
       unimap(i) = zero
    enddo

    if (getdebug() > 1) print*,"ptc_twiss"
    call momfirstinit()

    !------------------------------------------------------------------------------
    table_name = charconv(tab_name)
    summary_table_name = charconv(summary_tab_name)

    if (getdebug() > 1) then
       print*,"ptc_twiss: Table name is ",table_name
       print*,"ptc_twiss: Summary table name is", summary_table_name
    endif

    if(universe.le.0) then
       call fort_warn('return from ptc_twiss: ',' no universe created')
       call seterrorflag(1,"ptc_twiss ","no universe created till now");
       return
    endif
    if(index_mad.le.0) then
       call fort_warn('return from ptc_twiss: ',' no layout created')
       call seterrorflag(2,"ptc_twiss ","no layout created till now");
       return
    endif

    call cleartables() !defined in madx_ptc_knobs

    call kill_para(my_ring) !removes all the previous parameters

    nda = getnknobsall() !defined in madx_ptc_knobs
    suml=zero

    icase = get_value('ptc_twiss ','icase ')
    deltap0 = get_value('ptc_twiss ','deltap ')

    rmatrix = get_value('ptc_twiss ','rmatrix ')

    deltap = zero
    call my_state(icase,deltap,deltap0)
    CALL UPDATE_STATES

    ! check that deltap-dependency selected iff icase=5, which stands for
    ! 4 dimensions and deltap/p as an external parameter.
    ! indeed the simplified formulas we applied for derivation w.r.t deltap assume
    ! that deltap is an externally supplied constant parameter, which is no longer
    ! the case when icase = 6 as deltap then becomes a phase-space state-variable.
    ! and for icase=4, there is no dispersion.	
    
    deltap_dependency = get_value('ptc_twiss ','deltap_dependency ') .ne. 0
    
    if (deltap_dependency .and. .not.(icase.eq.5)) then
    	call fort_warn('ptc_twiss: ','derivation formulas w.r.t deltap assume deltap is fixed parameter')
	call fort_warn('ptc_twiss: ','derivation formulas w.r.t deltap therefore assume icase=5')
    endif
    
    x(:)=zero
    if(mytime) then
       call Convert_dp_to_dt (deltap, dt)
    else
       dt=deltap
    endif
    if(icase.eq.5) x(5)=dt

    closed_orbit = get_value('ptc_twiss ','closed_orbit ') .ne. 0

    if( closed_orbit .and. (icav .gt. 0) .and. (my_ring%closed .eqv. .false.)) then
       call fort_warn('return from ptc_twiss: ',' Closed orbit requested on not closed layout.')
       call seterrorflag(3,"ptc_twiss ","Closed orbit requested on not closed layout.");
       return
    endif

    if(closed_orbit) then

       if ( .not. c_%stable_da) then
          call fort_warn('ptc_twiss: ','DA got unstable even before finding orbit')
          call seterrorflag(10,"ptc_twiss ","DA got unstable even before finding orbit");
          stop
          return
       endif

       !print*, "Looking for orbit"
       call find_orbit(my_ring,x,1,default,c_1d_7)

       if ( .not. c_%stable_da) then
          call fort_warn('ptc_twiss: ','DA got unstable at attempt to find closed orbit')
          call seterrorflag(10,"ptc_twiss ","DA got unstable at attempt to find closed orbit");
          stop
          return
       endif
       CALL write_closed_orbit(icase,x)

    elseif(my_ring%closed .eqv. .true.) then
       print*, "Closed orbit specified by the user!"
       !CALL write_closed_orbit(icase,x) at this position it isn't read
    endif


    mynd2 = 0
    npara = 0
    no = get_value('ptc_twiss ','no ')
    if ( no .lt. 1 ) then
       call fort_warn('madx_ptc_twiss.f90 <ptc_twiss>:','Order in twiss is smaller then 1')
       print*, "Order is ", no
       return
    endif

    !this must be before initialization of the Bertz

    initial_distrib_manual = get_value('ptc_twiss ','initial_moments_manual ') .ne. 0
    if (initial_distrib_manual) then
       if (getdebug() > 1) then
          print*,"Initializing map with initial_moments_manual=true"
       endif
       call readinitialdistrib()
    endif


    call init(default,no,nda,BERZ,mynd2,npara)

    call alloc(y)
    y=npara
    Y=X

    call alloc(transfermap)
     
    if ( getnpushes() > 0 .or. rmatrix > 0) then
      tracktm = my_true
    else
      tracktm = my_false
    endif  
    
    
    if (tracktm) then
       transfermap = npara
       transfermap = X
    endif

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  INIT Y that is tracked          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call initmap(dt)

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

    slice_magnets = get_value('ptc_twiss ','slice_magnets ') .ne. 0
    

    if (.not. slice_magnets) then

       ! choice is to go one element after the other, instead of through the magnets' inner slices
       ! option is either absent or set to 'elements'

       do i=1,MY_RING%n

          if (getdebug() > 1) then
             write(6,*) "##########################################"
             write(6,'(i4, 1x,a, f10.6)') i,current%mag%name, suml
             write(6,'(a1,a,a1)') ">",current%mag%vorname,"<"
             write(6,'(a, f12.6, a)') "Ref Momentum ",current%mag%p%p0c," GeV/c"
             !          if (associated(current%mag%BN)) write(6,*) "k1=", current%mag%BN(2)
          endif

          if (nda > 0) then
             !         if (getnknobis() > 0) c_%knob = my_true
             !print*, "parametric",i,c_%knob
             call track(my_ring,y,i,i+1,+default)
             if (tracktm) then
                call track(my_ring,transfermap,i,i+1,+default)
             endif
          else
             call track(my_ring,y,i,i+1, default)
             if (tracktm) then
                call track(my_ring,transfermap,i,i+1,default)
             endif
          endif
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             call fort_warn('ptc_twiss: ','DA got unstable')
             call seterrorflag(10,"ptc_twiss ","DA got unstable ");
             if (getdebug() > 2) close(mf1)
             return
          endif

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) "ptc_twiss unstable (Twiss parameters) element: ",i," name: ",current%MAG%name,"-programs continues "
             write(whymsg,*) 'APERTURE error: ',why
             call fort_warn('ptc_twiss: ',whymsg)
             call seterrorflag(10,"ptc_twiss: ",whymsg);
             !          Write(6,*) why ! See produce aperture flag routine in sd_frame
             goto 100
          endif

          if (getdebug() > 2) then
             write(mf1,*) "##########################################"
             write(mf1,'(i4, 1x,a, f10.6)') i,current%mag%name, suml
             call print(y,mf1)
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
          
          call puttwisstable(transfermap) 
          
          if (tracktm) then
             call putusertable(i,current%mag%name,suml,getdeltae(),transfermap,y)
          endif



          iii=advance_node()
          current=>current%next
       enddo
100    continue

    else ! choice is to track along successive inner slices of the magnets
       call TrackAlongInnerSlices() ! for time-being, instead of the above
    endif


    ! relocated the following here to avoid side-effect
    print77=.false.
    read77=.false.

    call kanalnummer(mf2)
    open(unit=mf2,file='end.map')
    call print(y,mf2)
    close(mf2)
    ! relocated here to avoid side-effect


    if ( (momentumCompactionToggle .eqv. .true.)  .and. (getenforce6D() .eqv. .false.)) then
       ! only makes sense if the lattice is a ring (skipped for a line lattice)
       call oneTurnSummary()
       call set_option('ptc_twiss_summary ', 1)
    else
       call set_option('ptc_twiss_summary ',0) ! for time-being, do not support lines
    endif

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

  !  if (tracktm) then
       CALL kill(transfermap)
  !  endif
  
    do i=1,6
       call kill(unimap(i))
    enddo

    call finishknobs()

    if (savemaps) then  !do it at the end, so we are sure the twiss was successful
       mapsorder = no
       mapsicase = icase
       if (getnmoments() > 0) call ptc_moments(no*2) !calcualate moments with the maximum available order
    endif


    call f90flush(20,my_false)

    if (getdebug() > 2) close(mf1)

    !****************************************************************************************
    !*********  E N D   O F   PTC_TWISS      ************************************************
    !****************************************************************************************
    !________________________________________________________________________________________

  contains  ! what follows are internal subroutines of ptc_twiss
    !____________________________________________________________________________________________

    subroutine initmap(dt)
      implicit none
      integer  :: double_from_table
      integer  :: mman, mtab, mascr, mdistr !these variable allow to check if the user did not put too many options
      integer  :: mmap 
      real(dp) dt
      
      beta_flg = (get_value('ptc_twiss ','betx ').gt.0) .and. (get_value('ptc_twiss ','bety ').gt.0)

      mman  = get_value('ptc_twiss ','initial_matrix_manual ')
      mtab  = get_value('ptc_twiss ','initial_matrix_table ')
      mascr = get_value('ptc_twiss ','initial_ascript_manual ')
      mdistr = get_value('ptc_twiss ','initial_moments_manual ')
      mmap = get_value('ptc_twiss ','initial_map_manual ')

      
      momentumCompactionToggle = .false. ! set to true in the case the map
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
         k = double_from_table("map_table ", "nv ", 1, doublenum)
         if(k.ne.-1) then
            call liepeek(iia,icoast)
            my_nv=int(doublenum)
            nv_min=min(c_%npara,my_nv)
         else
            initial_matrix_table=.false.
         endif
      endif

      if(initial_matrix_table) then

         if (getdebug() > 1) then

            print*,"Initializing map with initial_matrix_table=true"
         endif
         call readmatrixfromtable()

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

         if (getdebug() > 1) then
            print*,"Initializing map from one turn map: Start Map"
            call print(y,6)
         endif

         momentumCompactionToggle = .true. ! compute momemtum compaction factor, tunes, chromaticies for ring

         call track(my_ring,y,1,default)
         if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
            call fort_warn('ptc_twiss: ','DA got unstable (one turn map production)')
            call seterrorflag(10,"ptc_twiss ","DA got unstable (one turn map production)");
            return
         endif

         call PRODUCE_APERTURE_FLAG(flag_index)
         if(flag_index/=0) then
            call ANALYSE_APERTURE_FLAG(flag_index,why)

            write(whymsg,*) 'APERTURE unstable (one turn map production) - programs continues: ',why
            call fort_warn('ptc_twiss: ',whymsg)
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
      type(work)              :: cfen !current fibre energy

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
         write(mf1,'(3(a, f10.6))') "Ref Momentum ",cfen%p0c," Energy ", cfen%energy," DeltaE ",getdeltae
      endif

    end function getdeltae
    !____________________________________________________________________________________________

    subroutine puttwisstable(transfermap)
      implicit none
      include "madx_ptc_knobs.inc"
      integer i1,i2,ii,i1a,i2a
      real(kind(1d0))   :: opt_fun(150),myx ! opt_fun(72) -> opt_fun(81) 
      ! increase to 150 to have extra space beyond what's needed to accomodate additional derivatives w.r.t. delta_p
      real(kind(1d0))   :: deltae
      type(real_8), target :: transfermap(6)

      if (getdebug() > 2) then
         write(mf1,*) "##########################################"
         write(mf1,*) ""
         write(mf1,'(i4, 1x,a, f10.6)') i,current%mag%name,suml
         write(mf1,*) ""
         call print(y,mf1)
      endif

      deltae = getdeltae()

      call double_to_table(table_name, 's ', suml)

      doublenum = deltae * startfen%energy
      call double_to_table(table_name, 'energy ', doublenum)


      opt_fun(:)=zero

      call liepeek(iia,icoast)
      allocate(j(c_%npara))
      j(:)=0
      do ii=1,c_%npara ! fish
         opt_fun(ii)=y(ii)%T.sub.j
      enddo

      myx=opt_fun(6)
      opt_fun(6)=opt_fun(5)
      opt_fun(5)=myx
      deallocate(j)

      ioptfun=6
      call vector_to_table(table_name, 'x ', ioptfun, opt_fun(1))

      
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
      call vector_to_table(table_name, 're11 ', ioptfun, opt_fun(1))



      opt_fun(beta11)= tw%beta(1,1) * deltae	! beta11=1
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
      opt_fun(62+4+8)=zero ! was 36 instead of 62 => on 9 march add 4 => on 3 July 2009 add 4
      do i1=1,c_%nd2
         if(i1.le.4) then
            i1a=i1
         elseif(i1.eq.5) then
            i1a=6
         else
            i1a=5
         endif
         do i2=1,c_%nd2
            if(i2.le.4) then
               i2a=i2
            elseif(i2.eq.5) then
               i2a=6
            else
               i2a=5
            endif
            ii=(62+4+8)+(i1a-1)*6+i2a ! was 36 instead of 62 => now (62+4) instead of 62 => 62+4+4
            opt_fun(ii)=tw%eigen(i1,i2) * deltae
	    ! where do these eigen values go? no such dedicated column defined in madx_ptc_knobs.inc
            if(mytime.and.i2a.eq.6) opt_fun(ii)=-opt_fun(ii)
         enddo
      enddo

      if (getdebug() > 2)  then
      	! this part of the code is out of sync with the rest
         write(6,'(a,1(f8.4,1x))') current%MAG%name,suml
         write(6,'(a,1(f10.8,1x))') "Delta E ", deltae
         write(6,'(a,3(i8.0,1x))')  "idxes ", beta11,beta22,beta33
         write(6,'(a,3(f8.4,1x))')  "betas raw   ", tw%beta(1,1),tw%beta(2,2),tw%beta(3,3)
         write(6,'(a,3(f8.4,1x))')  "betas w/ener", opt_fun(1),opt_fun(5),opt_fun(9)
         write(6,'(a,3(f8.4,1x))')  "disps       ", opt_fun(57),opt_fun(58),opt_fun(59),opt_fun(60)
         write(6,'(a,3(f8.4,1x))')  "tunes       ", tw%mu(1),tw%mu(2),tw%mu(3)
      endif

      ! the following works : we see the list of all elements in sequence - what about twiss_ptc_line & twiss_ptc_ring?
      ! jluc debug - begin
      !write(28,'(a,1(f8.4,1x))') current%MAG%name,suml
      ! jluc debug - end

      ! on July 3rd 2009, add another 8 for second/third derivatives of dispersions w.r.t. deltap
      ioptfun=81+4+8 !72->81 to accomodate additional derivatives w.r.t. delta_p => should one add 4 to this one, as above?
      ! actually 3*21+6 (???) elements from beta11 to include up to disp6p
      call vector_to_table(table_name, 'beta11 ', ioptfun, opt_fun(1))
      call augment_count(table_name)


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

    subroutine readmatrixfromtable
      implicit none
      integer  :: double_from_table

      x(:)=zero
      allocate(j(c_%npara))
      j(:)=0

      do i = 1,my_nv
         k   = double_from_table("map_table ", "coef ", i, doublenum)
         d_val=doublenum
         if(i.le.c_%npara) then
            x(i) = d_val-(y(i)%T.sub.j)
         endif
      enddo

      do i = 1,nv_min
         do ii = 1,nv_min
            j(ii)  = 1
            row    = i*my_nv+ii
            k   = double_from_table("map_table ", "coef ", row, doublenum)
            d_val=doublenum
            d_val  = d_val-(y(i)%T.sub.j)
            y(i)%T = y(i)%T + (d_val.mono.j)
            j(ii)=0
         enddo
      enddo

      call maptoascript()

      deallocate(j)
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
      call readrematrix() !reads re
      call readreforbit() !reads x
      call initmapfrommatrix()
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
      real(dp),dimension(6)::reval,aieval
      real(dp),dimension(6,6)::revec,aievec

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

      call eig6(re,reval,aieval,revec,aievec)
      do i=1,iia(4)-icoast(2)
         if(abs(reval(i)**2+aieval(i)**2 -one).gt.c_1d_10) then
            call fort_warn("ptc_twiss","Fatal Error: Eigenvalues off the unit circle!")
            stop
         endif
      enddo
    end subroutine initmapfrommatrix
    !_________________________________________________________________

    subroutine maptoascript
      !Performes normal form on a map, and plugs A_ in its place
      implicit none

      call alloc(normal)
      normal = y

      if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
         call fort_warn('ptc_twiss: ','Error: DA in twiss got unstable during Normal Form')
         call seterrorflag(10,"ptc_twiss ","DA in twiss got unstable during Normal Form");
         return
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
            stop
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
      if(icase.eq.5) x(5) = x(5) + dt


      if (getdebug() > 0 ) then
         CALL write_closed_orbit(icase,x)
      endif

      call reademittance()

      !Here we initialize Y(6)

      call alloc(be); call alloc(al); call alloc(di)

      !  code to power knows
      !


      do i=1,c_%nd
         be(i)=beta(i)
         al(i)=alpha(i)
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
               if (getdebug() > 1) print*,"Beta ",i," is knob no. ", inicondknobs%beta(i)
               call make_it_knob(be(i),k_system+inicondknobs%beta(i))
            endif

            if(inicondknobs%alfa(i)/=0) then
               if (getdebug() > 1) print*,"Alfa ",i," is knob no. ",  inicondknobs%alfa(i)
               call make_it_knob(al(i),k_system+inicondknobs%alfa(i))
            endif
         enddo

         do i=1,4
            if(inicondknobs%dispersion(i)/=0) then
               if (getdebug() > 1) print*,"Dispersion ",i," is knob no. ",  inicondknobs%dispersion(i)
               call make_it_knob(di(i),k_system+inicondknobs%dispersion(i))
            endif
         enddo

      endif


      y=x

      do i=1,c_%nd
         !         print*, " ", i, beta(i)
         !         call print(y(2*i-1),6)
         !         call print(y(2*i  ),6)

         y(2*i-1)= x(2*i-1) + sqrt(be(i)) * morph((one.mono.(2*i-1))    )
         y(2*i)= x(2*i) + one/sqrt(be(i)) * &
              (morph(  (one.mono.(2*i)) )-(al(i)) * morph((one.mono.(2*i-1))))

         !         call print(y(2*i-1),6)
         !         call print(y(2*i  ),6)
      enddo



      !--moments--!
      if( (c_%npara==5)       .or.  (c_%ndpt/=0) ) then

         if ( beta(3) .gt. zero ) then

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



    ! jluc
    ! compute momemtum-compaction factor in the same fashion it is carried-out in twiss.F

    subroutine oneTurnSummary()

      implicit none
      type(fibre), pointer :: fibrePtr
      real(dp) :: alpha_c, eta_c ! momentum-compaction factor & phase-slip factor
      real(dp) :: alpha_c_p ! first order derivative w.r.t delta-p/p
      real(dp) :: alpha_c_p2 ! second order derivative w.r.t delta-p/p
      real(dp) :: gamma_tr ! gamma_transition, or "transition energy" above which the particles' arrival time
      !  with respect to other particles is determined by its path length instead of by its velocity 
      real(dp) :: deltap
      real(dp) :: betaRelativistic, gammaRelativistic
      real(dp) :: suml ! cumulative length along the ring
      real(dp) :: fractionalTunes(3)
      real(dp) :: chromaticities(2)
      integer :: i,j
      real(dp) :: state(6) ! six-dimensional phase-space state (usually referred-to as 'x')
      real(dp) :: sd ! as in twiss.F
      type(real_8) :: oneTurnMap(6)
      type(real_8) :: theAscript(6) ! used here to compute dispersion's derivatives
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
      real(dp) :: dispersion_p2(4) ! second order of the dispersion w.r.t delta-p/p
      real(dp) :: dispersion_p3(4) ! third order of the dispersion w.r.t delta-p/p
      integer :: debugFiles
      integer :: icase
      integer :: order
      real(dp) :: dpathlength, dpathlength_p ! 6th coordinate and derivative w.r.t deltap
      
      order = get_value('ptc_twiss ', 'no ')
      
      ! should end-up gracefully here in case the topology of the lattice is not those of a closed-ring

      debugFiles = 0 ! set it to one and fort.21, fort.22 and fort.23 are created

      ! 1. track along the machine to get its cumulative length
      fibrePtr => my_ring%start
      suml = zero
      do i=1,my_ring%n
         suml = suml+fibrePtr%mag%p%ld
         fibrePtr => fibrePtr%next
         ! instead, could loop until we encounter the end-ground of the layout
      enddo

      ! 2. retreive the relativistic parameters beta and gamma
      ! (beta=v/c, gamma=E/mc^2 and gamma=1/sqrt(1-beta^2))
      betaRelativistic = get_value('probe ','beta ');
      gammaRelativistic = get_value('probe ','gamma ');

      ! 3. need the one-turn map's coefficients to reproduce twiss.F formulas for momentum compaction


      state(:)=zero

      ! added 9 march 2009 to take into dependency on deltap, as done in subroutine 'ptc_twiss'
      icase = get_value('ptc_twiss ','icase ') ! mind the trailing space
      deltap0 = get_value('ptc_twiss ','deltap ') ! mind the trailing space
      deltap = zero
      call my_state(icase,deltap,deltap0)
      if (mytime) then
         call Convert_dp_to_dt(deltap,dt)
      else
         dt = deltap
      endif
      if (icase.eq.5) then
         state(5) = dt
      endif
      ! end of added part

      call find_orbit(my_ring,state,1,default)

      ! write(6,*) 'NOW ORBIT IS ',state(:) => y, py, pt, -cT all ZERO !!! (idem in ptc_twiss)

      if (.not.c_%stable_da) then
         call fort_warn('ptc_twiss:','DA got unstable in momentum-compaction routine')
         call seterrorflag(10,"ptc_twiss","DA got unstable in momentum-compaction routine")
         stop
         return
      endif

      if (debugFiles .eq. 1) then
         write(21,*) "Closed orbit state=", state(1:6)
      endif


      ! need at least order 2 to get the chromaticities
      ! with an explicit call to init_default
      ! if we hadn't ned the chromaticities we could have simply discarded the following snippe
      call init_default(default,2,0)

      call alloc(oneTurnMap)
      oneTurnMap = npara ! ?
      oneTurnMap = state
      ! or should we create it as x+id?
      call track(my_ring,oneTurnMap,1,default)

      ! now retrieve the one-turn map's coefficients
      if (debugFiles .eq. 1) then
         do i=1,6
            do j=1,6
               write(21,*) "r(",i,j,")=",oneTurnMap(i).sub.coeffSelector(j,:)
            enddo
         enddo
      endif

      ! 4. retreive the dispersion coefficients
      ! (may be the coefficient of delta of the map?)
      ! decompose the map via a normal form to get the dispersion...   
      call alloc(theNormalForm)
      theNormalForm = oneTurnMap
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

      ! 5. apply formulas from twiss.F: sd = r(5,6)+r(5,1)*disp(1)+...+r(5,4)*disp(4)
      sd = - oneTurnMap(6).sub.coeffSelector(5,:) ! 5/6 swap MADX/PTC
      do i=1,4
         sd = sd - (oneTurnMap(6).sub.coeffSelector(i,:))*dispersion(i)
      enddo
      
      dpathlength = -sd ! added 18 august 2009 for latter reuse of the map's six coordinate
      ! dpathlength is the derivative of the pathlength w.r.t deltap
      
      eta_c = -sd * betaRelativistic**2 / suml
      alpha_c = one / gammaRelativistic**2 + eta_c
      gamma_tr = one / sqrt(alpha_c)

      ! compute delta-p/p dependency of alpha_c
      ! first order derivatives of the dispersions
      ! assuming icase=5
      !call daprint(oneTurnMap,88)
      if (icase.eq.5) then
      	
      	! always assume time=false, so that the fifth phase-space variable is deltap instead of pt!
	! otherwise should issue a warning or an error!
 	
      	call alloc(theAscript)
        ! proceed slightly differently to retreive the dispersion's derivatives
	! (so far only managed to get the dispersion on the normal form, but not its derivatives)
       	theAscript = state + theNormalForm%a_t
	
	! note: should reuse information computed previously than redo this...
	
      	! if icase=5, Taylor series expansion disp = disp + delta_p
	do i=1,4
		! apparently, we are up by a factor two
      		dispersion_p(i) = 2.0*(theAscript(i)%t.sub.'000020') ! as usual in this file
	enddo

	! compute derivative of the formula for eta_c

	sd = 0.0
	! first dr65/ddeltap
	sd = sd		- 1.0*(oneTurnMap(6).sub.'100010')*dispersion(1) &
	&		- 1.0*(oneTurnMap(6).sub.'010010')*dispersion(2) & 
	&		- 1.0*(oneTurnMap(6).sub.'001010')*dispersion(3) &
	&		- 1.0*(oneTurnMap(6).sub.'000110')*dispersion(4) &
	&		- 2.0*(oneTurnMap(6).sub.'000020') &
	! new
	&		- 1.0*(oneTurnMap(6).sub.'000011')*(oneTurnMap(6).sub.'000001') &
	! then terms dr6i/ddeltap*dispersion(i)
	&		- 2.0*(oneTurnMap(6).sub.'200000')*dispersion(1)**2 &
	&		- (oneTurnMap(6).sub.'110000')*dispersion(1)*dispersion(2) &
	&		- (oneTurnMap(6).sub.'101000')*dispersion(1)*dispersion(3) &
	&		- (oneTurnMap(6).sub.'100100')*dispersion(1)*dispersion(4) &
	! new
	&		- (oneTurnMap(6).sub.'100010')*dispersion(1) &
	! dr62/ddeltap*disperion(2)
	&		- (oneTurnMap(6).sub.'110000')*dispersion(1)*dispersion(2) &
	&		- 2.0*(oneTurnMap(6).sub.'020000')*dispersion(2)**2 &
	&		- (oneTurnMap(6).sub.'011000')*dispersion(2)*dispersion(3) &
	&		- (oneTurnMap(6).sub.'010100')*dispersion(2)*dispersion(4) &
	! new
	&		- (oneTurnMap(6).sub.'010010')*dispersion(2) &
	! dr63/ddeltap*dispersion(3)
	&		- (oneTurnMap(6).sub.'101000')*dispersion(1)*dispersion(3) &
	&		- (oneTurnMap(6).sub.'011000')*dispersion(2)*dispersion(3) &
	&		- 2.0*(oneTurnMap(6).sub.'002000')*dispersion(3)**2 &
	&		- (oneTurnMap(6).sub.'001100')*dispersion(3)*dispersion(4) &
	! new
	&		- 1*(oneTurnMap(6).sub.'001010')*dispersion(3) &
	! dr64/ddeltap*dispersion(4)
	&		- (oneTurnMap(6).sub.'100100')*dispersion(1)*dispersion(4) &
	&		- (oneTurnMap(6).sub.'010100')*dispersion(2)*dispersion(4) &
	&		- (oneTurnMap(6).sub.'001100')*dispersion(3)*dispersion(4) &
	&		- 2.0*(oneTurnMap(6).sub.'000200')*dispersion(4)**2 &
	! new
	&		- 1*(oneTurnMap(6).sub.'000110')*dispersion(4)

	! terms involving derivatives of the dispersions
	do i=1,4
		sd = sd - (oneTurnMap(6).sub.coeffSelector(i,:))*dispersion_p(i)
	enddo
	
	dpathlength_p = sd ! for later reuse of the sixth coordinate of the map
	dpathlength_p = 0.0
	
	alpha_c_p = sd * (betaRelativistic**2) / suml
	
	! now compute second order derivatives alpha_c_p2 (assumes map of order 3 at least)
	! to this end, we first need to compute second order derivatives of the dispersion
	! extrapolated from how we retreive the first order derivative...
	if (order.lt.3) then
		alpha_c_p2 = 0.0
         	call fort_warn('ptc_twiss:','order too low to compute alpha_c_p2 (should be >=3)')	
	else
		do i=1,4
      			dispersion_p2(i) = 3.0*2.0*(theAscript(i)%t.sub.'000030') ! to be checked
		enddo
		
		sd = 0.0
		! derivative of -2.0*oneTurnMap(6).sub.'100010'*dispersion(1)
		sd = 		- 1.0 * (&
		! negative factor out of the parenthesis
! DO NOT REMOVE TAG BELOW
!----- INSERT CODE --- START
! code start insert on Tue Sep  1 16:08:16 2009
! warning: did not check yet that the following formula is correct
! at the moment, there's a numerical discrepancy with the finite-difference approximation of this derivative
& +2.0*(oneTurnMap(6).sub.'200010')*dispersion(1)*dispersion(1) &
& +(oneTurnMap(6).sub.'110010')*dispersion(2)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'100020')*dispersion(1) &
& +(oneTurnMap(6).sub.'100010')*dispersion_p(1) &
& +(oneTurnMap(6).sub.'110010')*dispersion(1)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'020010')*dispersion(2)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'010020')*dispersion(2) &
& +(oneTurnMap(6).sub.'010010')*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'100020')*dispersion(1)*0.5*2.0 &
& +2.0*(oneTurnMap(6).sub.'010020')*dispersion(2)*0.5*2.0 &
& +6.0*(oneTurnMap(6).sub.'000030')*0.5*2.0 &
& +6.0*(oneTurnMap(6).sub.'300000')*dispersion(1)*0.5*2.0*dispersion(1)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'210000')*dispersion(2)*0.5*2.0*dispersion(1)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'200010')*0.5*2.0*dispersion(1)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'200000')*dispersion_p(1)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'210000')*dispersion(1)*dispersion(2)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'120000')*dispersion(2)*dispersion(2)*dispersion(1) &
& +(oneTurnMap(6).sub.'110010')*dispersion(2)*dispersion(1) &
& +(oneTurnMap(6).sub.'110000')*dispersion_p(2)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'200010')*dispersion(1)*dispersion(1) &
& +(oneTurnMap(6).sub.'110010')*dispersion(2)*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'100020')*dispersion(1) &
& +2.0*(oneTurnMap(6).sub.'200000')*dispersion(1)*dispersion_p(1) &
& +(oneTurnMap(6).sub.'110000')*dispersion(2)*dispersion_p(1) &
& +(oneTurnMap(6).sub.'100010')*dispersion_p(1) &
& +6.0*(oneTurnMap(6).sub.'300000')*dispersion(1)*0.5*dispersion(1)*2.0*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'210000')*dispersion(2)*0.5*dispersion(1)*2.0*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'200010')*0.5*dispersion(1)*2.0*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'200000')*dispersion_p(1)*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'210000')*dispersion(1)*dispersion(2)*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'120000')*dispersion(2)*dispersion(2)*dispersion_p(1) &
& +(oneTurnMap(6).sub.'110010')*dispersion(2)*dispersion_p(1) &
& +(oneTurnMap(6).sub.'110000')*dispersion_p(2)*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'200010')*dispersion(1)*dispersion_p(1) &
& +(oneTurnMap(6).sub.'110010')*dispersion(2)*dispersion_p(1) &
& +2.0*(oneTurnMap(6).sub.'100020')*dispersion_p(1) &
& +(oneTurnMap(6).sub.'100000')*dispersion_p2(1) &
& +2.0*(oneTurnMap(6).sub.'210000')*dispersion(1)*dispersion(1)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'120000')*dispersion(2)*dispersion(1)*dispersion(2) &
& +(oneTurnMap(6).sub.'110010')*dispersion(1)*dispersion(2) &
& +(oneTurnMap(6).sub.'110000')*dispersion_p(1)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'120000')*dispersion(1)*0.5*2.0*dispersion(2)*dispersion(2) &
& +6.0*(oneTurnMap(6).sub.'030000')*dispersion(2)*0.5*2.0*dispersion(2)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'020010')*0.5*2.0*dispersion(2)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'020000')*dispersion_p(2)*dispersion(2) &
& +(oneTurnMap(6).sub.'110010')*dispersion(1)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'020010')*dispersion(2)*dispersion(2) &
& +2.0*(oneTurnMap(6).sub.'010020')*dispersion(2) &
& +(oneTurnMap(6).sub.'110000')*dispersion(1)*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'020000')*dispersion(2)*dispersion_p(2) &
& +(oneTurnMap(6).sub.'010010')*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'210000')*dispersion(1)*dispersion(1)*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'120000')*dispersion(2)*dispersion(1)*dispersion_p(2) &
& +(oneTurnMap(6).sub.'110010')*dispersion(1)*dispersion_p(2) &
& +(oneTurnMap(6).sub.'110000')*dispersion_p(1)*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'120000')*dispersion(1)*0.5*dispersion(2)*2.0*dispersion_p(2) &
& +6.0*(oneTurnMap(6).sub.'030000')*dispersion(2)*0.5*dispersion(2)*2.0*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'020010')*0.5*dispersion(2)*2.0*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'020000')*dispersion_p(2)*dispersion_p(2) &
& +(oneTurnMap(6).sub.'110010')*dispersion(1)*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'020010')*dispersion(2)*dispersion_p(2) &
& +2.0*(oneTurnMap(6).sub.'010020')*dispersion_p(2) &
& +(oneTurnMap(6).sub.'010000')*dispersion_p2(2) &
! code end insert on Tue Sep  1 16:08:16 2009
!----- INSERT CODE --- END
! DO NOT REMOVE TAG ABOVE
		& 		)
		alpha_c_p2 = sd * betaRelativistic**2 / suml
	endif
			
	call kill(theAscript)
      else
      	alpha_c_p = 0.0 ! exactly zero means failure to compute the actual value
	alpha_c_p2 = 0.0
      endif ! if (icase.eq.5)
      
      ! also output the tune ...
      fractionalTunes = theNormalForm%tune
      ! the above is exactly equivalent to the following two lines (i.e. returns frac.tune)
      fractionalTunes(1) = theNormalForm%DHDJ%v(1).sub.'0000' ! as in So_fitting.f90
      fractionalTunes(2) = theNormalForm%DHDJ%v(2).sub.'0000' ! as in So_fitting.f90
      ! Q: is it possible to get the actual total tune, as returned by twiss.F?
      ! => no, not with a map...

      ! ... as well as the chromaticities
      chromaticities(1) = theNormalForm%DHDJ%v(1).sub.'00001' ! as in So_fitting.f90
      chromaticities(2) = theNormalForm%DHDJ%v(2).sub.'00001' ! as in So_fitting.f90
      ! to get chromaticities, went to higher order with above "call init_default(default,2,0)"

      ! for debug: check the values by printing the map
      if (debugFiles .eq. 1) then
         call daprint(oneTurnMap,25) ! prints the one-turn map on file 25
         call daprint(theNormalForm%dhdj,26) ! print tunes, chromaticities and anharmonicities
         ! as done in madx_ptc_normal.f90
      endif

      call kill(theNormalForm)

      call kill(oneTurnMap)


      !write(6,*)
      !write(6,*) "Momentum compaction (momentum compaction factor & phase-slip factor):"
      !write(6,*) "alpha=", alpha_c, "eta=", eta_c
      !write(6,*) "fractional tunes=", fractionalTunes
      !write(6,*) "chromaticities=", chromaticities

      deltap = get_value('ptc_twiss ','deltap ')

      ! write the data into the ptc_twiss_summary table
      ! warning: must preserve the order EXACTLY
      call double_to_table( summary_table_name, 'length ', suml ) ! total length of the machine
      call double_to_table( summary_table_name, 'alpha_c ', alpha_c ) ! momemtum compaction factor
      call double_to_table( summary_table_name, 'alpha_c_p ', alpha_c_p) ! derivative w.r.t delta-p/p
      call double_to_table( summary_table_name, 'alpha_c_p2 ', alpha_c_p2) ! 2nd order derivative
      call double_to_table( summary_table_name, 'eta_c ', eta_c ) ! associated phase-slip factor
      call double_to_table( summary_table_name, 'gamma_tr ', gamma_tr) ! associated transition energy
      call double_to_table( summary_table_name, 'q1 ', fractionalTunes(1))
      call double_to_table( summary_table_name, 'q2 ', fractionalTunes(2))
      call double_to_table( summary_table_name, 'dq1 ', chromaticities(1))
      call double_to_table( summary_table_name, 'dq2 ', chromaticities(2))
      ! write the extremas of the Twiss functions
      ! for the time-being, do not bother about the coupling terms
      call double_to_table( summary_table_name, 'beta_x_min ', minBeta(1,1))
      call double_to_table( summary_table_name, 'beta_x_max ', maxBeta(1,1))
      call double_to_table( summary_table_name, 'beta_y_min ', minBeta(2,2))
      call double_to_table( summary_table_name, 'beta_y_max ', maxBeta(2,2))
      call double_to_table( summary_table_name, 'deltap ', deltap)
      ! the 6-d closed orbit
      call double_to_table( summary_table_name,'orbit_x ',state(1))
      call double_to_table( summary_table_name,'orbit_px ', state(2))
      call double_to_table( summary_table_name,'orbit_y ', state(3))
      call double_to_table( summary_table_name,'orbit_py ', state(4))
      ! warning: if 'time=false', the last two phase-space state-variables
      ! should be deltap/p and path-length respectively
      call double_to_table( summary_table_name,'orbit_pt ', state(5))
      call double_to_table( summary_table_name,'orbit_-cT ', state(6))

      call augment_count( summary_table_name ); ! only one row actually...

    end subroutine oneTurnSummary


    subroutine TrackAlongInnerSlices()

      ! the ptc_twiss shall feature an additional flag to decide
      ! whether or not one evaluate the Twiss functions inside the elements
      ! or solely at the middle.

      implicit none
      integer :: initialThinLensPos, thinLensPos
      type(integration_node), pointer :: nodePtr
      type(fibre), pointer :: fibrePtr
      real(dp) :: s
      integer(dp) :: knobsNumber
      real(dp) :: state(6) ! six-dimensional phase-space state (usually referred-to as 'x')
      type(real_8) :: theTransferMap(6) ! transfer map between successive inner thin-lenses
      type(normalform) :: theNormalForm
      type(real_8) :: theAscript(6) ! the phase-advance
      integer :: returnedInteger

      knobsNumber = nda ! nda is a semi-global variable !!!

      call make_node_layout(my_ring) ! essential: the way to look inside the magnets

      state = 0.d0
      call find_orbit(my_ring,state,1,default,1.d-5) ! 1 for the first element

      call alloc(theTransferMap) ! transfer map between two successive inner slices
      ! is it the correct way to initialize, or should we rely on y=x+id??
      theTransferMap = npara ! to be checked later-on
      theTransferMap = state

      call track(my_ring, theTransferMap, 1, default)

      call alloc(theNormalForm)
      call alloc(theAscript)

      theNormalForm = theTransferMap ! decompose into a normal form
      theAscript = state + theNormalForm%A_t ! linear part of the normal form

      call kill(theNormalForm) ! theAscript only is of interest from now on...

      ! my_ring%nthin does return 0, instead of the actual value, which we can get from my_ring%t%n

      fibrePtr => my_ring%start
      initialThinLensPos = fibrePtr%t1%pos ! t1 points to the first integration node
      nodePtr => fibrePtr%t1 ! t1 points to the first integration node

      do thinLensPos = initialThinLensPos, my_ring%t%n+initialThinLensPos-1 
         ! "t" is the child thin-lens layout according to "type layout" definition

         s = nodePtr%next%s(3) ! s(1) is the total arc-length, s(3) the total integration-distance
         ! From Etienne:
         ! nodePtr%s(1) is the total arc length (LD)
         ! nodePtr%s(3) is the total integration distance.... ( LC)  LC=LD unless exact=true and true rbends are present.
         ! nodePtr%s(4) is the local LC. 
         ! nodePtr%s(5) is the local LD. 


         ! in  a closed-ring, s evaluates to zero when we reach my_ring%t%n+initialThinLensPos -1
         ! let's omit such final value in case of a closed-ring
         if ((s .eq. 0d0) .and. (thinLensPos .eq. (my_ring%t%n+initialThinLensPos-1))) then
            s = nodePtr%s(3) + nodePtr%next%s(5) ! s of previous node + local offset
         endif
       
         ! following part is cut/pasted from the equivalent reference
         ! code that successively evaluates Twiss at elements' ends
         if (knobsNumber > 0) then
            call track_probe_x(my_ring,theAscript,+default, & ! +default in case of extra parameters !?
                 & node1=thinLensPos,node2=thinLensPos+1)
            if (getnpushes() > 0) then ! try to understand this part later-on
               call track_probe_x(my_ring,theTransferMap,+default, & ! +default in case of extra parameters !?
                    & node1=thinLensPos,node2=thinLensPos+1)
            endif
         else
            ! do we really need track_probe_x, or simply track to go along inner slices?
            call track_probe_x(my_ring,theAscript,default, &
                 & node1=thinLensPos,node2=thinLensPos+1)
            !write(27,*) "beta new=", (theTransferMap(1).sub.'1')**2+(theTransferMap(1).sub.'01')**2
            if (getnpushes() > 0) then ! try to understand this part later-on
               call track_probe_x(my_ring,theTransferMap,default, &
                    & node1=thinLensPos,node2=thinLensPos+1)             
            endif
         endif

         if ((.not. check_stable) .or. (.not. c_%stable_da)) then
            call fort_warn('ptc_twiss:','DA got unstable during tracking through the thin-lens slices')
            call seterrorflag(10,"ptc_twiss","DA got unstable during tracking through the thin-lens slices")
            return ! try to continue to the next thin-lens position
         endif

         ! following part is cut/pasted from the equivalent reference
         ! code that successively evaluates Twiss at elements' ends
         call produce_aperture_flag(flag_index)
         if (flag_index/=0) then
            call analyse_aperture_flag(flag_index, why)
            write(6,*) "ptc_twiss unstable (Twiss inside magnets) element:",fibrePtr%mag%name,"-program continues"
            write(whymsg,*) 'APERTURE error:',why
            call fort_warn('ptc_twiss:', whymsg)
            call seterrorflag(10,"ptc_twiss:", whymsg)
            goto 200 ! skip loop altogether
         endif

         tw = theAscript ! set the twiss parameters, with y being equal to the A_ phase advance

         ! update some semi-global values before invoking puttwisstable!!!
         current => fibrePtr ! current defined outside this subroutine and required by puttwisstable!!!
         suml = s ! another global!!! to be updated later-on with the actual s within the magnet
         ! suml is used internally to puttwisstable to save the curvilign abciss...

         call puttwisstable(transfermap) ! writes the resulting tw above to an internal table

         if (getnpushes() > 0) then ! !writes user selected map coeffs to user created tables. See ptc_select 
                                      !also saves twiss parameters and used defined variables as Taylor series 
                                      !             in function of knobs. See ptc_knob
            ! presently, do not enter here anyway
            call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap, theAscript) 
         endif                                                                                 


         if (associated(nodePtr,fibrePtr%t1)) then
            write(24,*) s, thinLensPos, "located at the beginning of ", fibrePtr%mag%name
         endif
         if (associated(nodePtr,fibrePtr%tm)) then
            write(24,*) s, thinLensPos, "located at the middle of ", fibrePtr%mag%name
            ! if option is to evaluate Twiss parameters at the middle, then invoke computation here
         endif
         ! Note: sometime looks like beginning is immediately followed by the middle

         nodePtr => nodePtr%next

         if (associated(nodePtr,fibrePtr%t2)) then ! t2 is last integration node along the fibre
            write(24,*) s, thinLensPos, "located at the end of ", fibrePtr%mag%name
            fibrePtr => fibrePtr%next
            ! indicate the complete_twiss_table code in madxn.c that we moved to the next element
            ! so that the element name on the far left displays correctly
            returnedInteger = advance_node()
         endif

      enddo ! for the successive slices / thin-lenses in the magnets' sequence

200   continue

      call kill(theTransferMap)
      call kill(theAscript)

    end subroutine TrackAlongInnerSlices

    subroutine TrackAlongMagnets()
      ! should cut/paste the tracking part of the ptctwiss main subroutine
    end subroutine TrackAlongMagnets
  
  END subroutine ptc_twiss

  subroutine computeDeltapDependency(y,s1)
    implicit none
    type(real_8), intent(in)  :: y(ndd)
    type(twiss),  intent(inout)  :: s1
    integer :: k,i
    integer :: J(lnv) ! the map's coefficient selector, as usual
    integer :: Jderiv(lnv) ! to store the map's coefficient selector of the derivative w.r.t deltap
    real(kind(1d0)) :: get_value ! C-function
    integer :: no ! order must be at equal to 2 to be able to get terms of the form x*deltap
    ! required to evaluate the derivatives of Twiss parameters w.r.t deltap
    integer :: ndel ! as in subroutine 'equaltwiss'...

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
          s1%disp_p(i) = 2.0*(y(i)%t.sub.'000020')
	  if (no.gt.2) then
		  s1%disp_p2(i) = 3.0*2.0*(y(i)%t.sub.'000030') ! assume at least order 3 for the map
	  else
	  	call fort_warn("ptc_twiss ","assume no>=3 for dispersion's 2nd order derivatives w.r.t delta-p")
	  endif
	  if (no.gt.3) then
		  s1%disp_p3(i) = 4.0*3.0*2.0*(y(i)%t.sub.'000040') 
		  ! assume at least order 4 for the map
	  else
	  	call fort_warn("ptc_twiss ","assume no>=4 for dispersion's 3rd order derivatives w.r.t delta-p")
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
  
! --- a set of routines to track extremas of Twiss functions
subroutine trackBetaExtrema(i,j,value)
  implicit none
  logical, save :: firstTime(3,3) = .true.
  integer :: i,j
  real(dp) :: value
  if (firstTime(i,j)) then
     firstTime(i,j) = .false.
     minBeta(i,j) = value
     maxBeta(i,j) = value
!     write(80,*) 'first time in trackBetaExtrema for ',i,j
  else
     if (minBeta(i,j) .gt. value) then
        minBeta(i,j) = value
     elseif (maxBeta(i,j) .lt. value) then
        maxBeta(i,j) = value
     endif
  endif
end subroutine trackBetaExtrema
! --- end of set of routines


end module madx_ptc_twiss_module


