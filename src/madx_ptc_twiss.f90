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
  type(probe_8)            :: theTransferMap
  type(probe_8)            :: theRDTs
  type(universal_taylor)   :: unimap(6)

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
  integer,  private                     :: np, no

  !new lattice function
  real(dp), private, dimension(3)       :: testold
  real(dp), private, dimension(3)       :: phase
  real(dp), private, allocatable, dimension(:,:,:)        :: savedTM

  character(len=5), private, dimension(5), parameter :: str5 = (/'10000','01000','00100','00010','00001'/)
  integer, private, dimension(6,6,3 )    :: Iaa ! for i=1,3 Ia(2*i-1,2*i-1,i) =1.d0;  Ia(2*i,2*i,i)   = 1.d0;
  integer, private, allocatable          :: J(:)
  integer, private, dimension(6)         :: j5 = (/0,0,0,0,1,0/)
  integer, private, dimension(6)         :: j6 = (/0,0,0,0,0,1/)
  integer, private, dimension(6,6)       :: fo = &
	            reshape(   (/1,0,0,0,0,0,&
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

  character(48)           :: nl_table_name='nonlin'
  character(48)           :: rdt_table_name='twissrdt'

  !============================================================================================
  !  PRIVATE
  !    routines
  private zerotwiss,equaltwiss,killtwiss



contains
  !____________________________________________________________________________________________

  subroutine initIaaMatrix()
    implicit none
    integer i

    Iaa = 0;
    do i=1,3
      Iaa(2*i-1,2*i-1,i) =1;  Iaa(2*i,2*i,i)   = 1;
    enddo

  end subroutine initIaaMatrix

  subroutine dispersion6D(A_script,disp)
    implicit none
!    type(probe_8), intent(in)::A_script_probe
    type(real_8) ::A_script(6)
    real(dp), intent(out)  :: disp(4)
    real(dp)  :: amatrix(6,6)
    real(dp) ::amatrix_inv(6,6),Ha(6,6,3)
    integer i,j, inverr
    ! H based dispersion as in Chao-Sands paper

    do i=1,6
      do j=1,6
        amatrix(i,j) = A_script(i).sub.fo(j,:)
      enddo
    enddo

    inverr = 0
    call matinv(amatrix, amatrix_inv,6,6,inverr) ! amap**(-1) ! and invert it and copy to 6:6 matrix
    if (inverr .ne. 0) then
      call fort_warn('ptc_twiss::dispersion6D ',' Can not invert A_script linear')
      disp = zero;
    endif


    do i=1,c_%nd

      Ha(1:6,1:6,i)=matmul(matmul(amatrix,Iaa(1:6,1:6,i)),amatrix_inv)  ! (15b)

    enddo

   do i=1,4
     disp(i) = Ha(i,5,3)/Ha(5,5,3)
   enddo


  end subroutine dispersion6D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !version on polynomials to non-linear dispersion
  subroutine dispersion6Dp(A_script,disp)
    implicit none
    type(real_8) ::A_script(6)
    real(dp), intent(out)  :: disp(4)
    type(taylor)  dispT(4)
    type(damap)  :: amap, dispMap
    type(damap)  :: Iaa_map
    real(dp)  :: amatrix(6,6)
    real(dp) ::amatrix_inv(6,6),Ha(6,6,3), v
    integer i,j,k
    ! H based dispersion as in Chao-Sands paper


    call alloc(dispMap)
    call alloc(dispT)
    call alloc(amap)
    amap = A_script ! move the transformation to type damap

    print*,"++++++++++++++"
    print*, "amap"
    call print(amap,6)

    amatrix    =amap       ! now we can copy to simple 6:6 matrix
    amatrix_inv=amap**(-1) ! and invert it and copy to 6:6 matrixs

    call alloc(Iaa_map)

    do i=1,c_%nd

      Ha(1:6,1:6,i)=matmul(matmul(amatrix,Iaa(1:6,1:6,i)),amatrix_inv)  ! (15b)

    enddo

    do i=3,3 ! i=1,c_%nd
      Iaa_map = zero;
      do j=1,6
        print*, Iaa(j,:,i)
        do k=1,6
          v = Iaa(j,k,i)
          Iaa_map%v(j) = Iaa_map%v(j) + v*(1.0_dp.mono.fo(k,:))
        enddo
      enddo
      print*,"++++++++++++++"
      print*, "Iaa ", i
      call print(Iaa_map,6)

      dispMap = amap * Iaa_map * amap**(-1)


      !do ii=1,6
      !  print*,"        ",Ha(1:6,ii,i)
      !enddo

    enddo


   dispMap = dispMap * (1.0_dp /  Ha(5,5,3));

   ! print*,"++++++++++++++"
   ! print*, "dispMap"
   !call print(dispMap,6)


   do i=1,4
     disp(i) = Ha(i,5,3)/Ha(5,5,3)

     dispT(i) = dispMap%v(i).par.fo(5,:)

     if (getdebug() > 2) then
       print*, 'Disp ', i, ' R', disp(i) , ' Map ', dispMap%v(i).sub.fo(5,:)
       print*, "Taylor "
       call print(dispT(i),6)
     endif

   enddo

   call kill(amap)

  end subroutine dispersion6Dp

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine equaltwiss(s1,A_script)
    implicit none
    type(twiss), intent(inout)::s1
    type(real_8), intent(in) ::A_script(6)
    real(dp)  :: amatrix(6,6) ! first order A_script
    integer jj,i,k, ndel
    real(dp) :: lat(0:6,6,3)=0
    real(dp) :: test, dph
    real(dp) :: epsil=1e-12  !
    integer  :: J(lnv)
    real(dp)  :: beta(3)

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) ' check_stable ',check_stable,' c_%stable_da ',c_%stable_da,' PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('equaltwiss CHECK 0 : ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"equaltwiss CHECK 0 ",whymsg)
      return
    endif


    lat = zero

    ndel=0
    if(c_%ndpt/=0 .and. isRing )  then
      ! in case of icase=56 the closed solution in long is not searched
       !print*, "We are at the mode 6D + nocav"
       ndel=1  !this is 6D without cavity (MADX icase=56) for a ring
    endif
    ! calculation of alpha, beta, gamma
    J=0;

    !print*, "EqualTwiss nd=",c_%nd," nd2=",c_%nd2, " ndel=", ndel

    beta(1) = (A_script(1).sub.'100000')**2 + (A_script(1).sub.'010000')**2
    beta(2) = (A_script(3).sub.'001000')**2 + (A_script(3).sub.'000100')**2
    if ( c_%nd > 2 ) beta(3) = (A_script(6).sub.'000010')**2 + (A_script(6).sub.'000001')**2

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) ' check_stable ',check_stable,' c_%stable_da ',c_%stable_da,' PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('equaltwiss CHECK 1 : ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"equaltwiss CHECK 1 ",whymsg)
      return
    endif

    !print*, " Skowron 1 EqTwiss ", beta

    do i=1,c_%nd2-2*ndel
       do jj=i,c_%nd2-2*ndel
          do k=1,c_%nd-ndel
             J(2*k-1)=1
             lat(i,jj,k)=              (A_script(i).sub.J)*(A_script(jj).sub.J)
             J(2*k-1)=0
             J(2*k)=1
             lat(i,jj,k)=lat(i,jj,k) + (A_script(i).sub.J)*(A_script(jj).sub.J)
             lat(jj,i,k)=lat(i,jj,k)
             J(2*k)=0
          enddo
       enddo
    enddo

    !print*,"BETZ=",(A_script(6).sub.'000010')**2 + (A_script(6).sub.'000001')**2

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) ' check_stable ',check_stable,' c_%stable_da ',c_%stable_da,' PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('equaltwiss CHECK 2 : ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"equaltwiss CHECK 2 ",whymsg)
      return
    endif

    J=0
    !here ND2=4 and delta is present      nd2=6 and delta is a constant
    !      print*,"nv",c_%nv,"nd2",c_%nd2,"np",c_%np,"ndpt",c_%ndpt ,"=>",c_%nv-c_%nd2-c_%np
    if( (c_%npara==5)       .or.  (c_%ndpt/=0)  ) then
       !when there is no cavity it gives us dispersions
       do i=1,4
          lat(0,i,1)=(A_script(i).sub.J5)
       enddo
    elseif (c_%nd2 == 6) then

      call dispersion6D(A_script,lat(0,1:4,1))

    else
       do i=1,4
          lat(0,i,1)=zero
       enddo
    endif

    !when there is no cavity it gives us dispersions
    do i=1,c_%nd2-2*ndel
       s1%disp(i)=lat(0,i,1)
    enddo

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) ' check_stable ',check_stable,' c_%stable_da ',c_%stable_da,' PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('equaltwiss CHECK 3 : ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"equaltwiss CHECK 3 ",whymsg)
      return
    endif

    !!!!!!!!!!!!!!!!
    ! phase advance!
    !!!!!!!!!!!!!!!!

    k = 2
    if(c_%nd2==6.and.c_%ndpt==0) k = 3

    j=0
    do i=1, k
       jj=2*i -1
       if (i == 3) then
          TEST = ATAN2((A_script(6).SUB.fo(5,:)),(A_script(6).SUB.fo(6,:)))/TWOPI
       else
          TEST = ATAN2((A_script(2*i -1).SUB.fo(2*i,:)),(A_script(2*i-1).SUB.fo(2*i-1,:)))/TWOPI
       endif


       IF(TEST<zero.AND.abs(TEST)>EPSIL)TEST=TEST+one
       DPH=TEST-TESTOLD(i)
       IF(DPH<zero.AND.abs(DPH)>EPSIL) DPH=DPH+one
       IF(DPH>half) DPH=DPH-one

       PHASE(i)=PHASE(i)+DPH
       TESTOLD(i)=TEST

    enddo


    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) ' check_stable ',check_stable,' c_%stable_da ',c_%stable_da,' PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('equaltwiss CHECK 4 : ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"equaltwiss CHECK 4 ",whymsg)
      return
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*," "
    !print*, "beta 11 = ", lat(1,1,1)
    !print*, "beta 22 = ", lat(3,3,2)
    !print*, "beta 22 = ", lat(5,5,3)
    !print*,"  manual 22 ", (A_script(3).sub.'0010')**2 + (A_script(3).sub.'0001')**2

    do i=1,c_%nd
       do k=1,c_%nd
          s1%beta(k,i)= lat(2*k-1,2*k-1,i)
          s1%alfa(k,i)=-lat(2*k-1,2*k  ,i)
          s1%gama(k,i)= lat(2*k  ,2*k  ,i)
       enddo
    enddo

    !print*, " Skowron EqTwiss ", s1%beta
    ! --- derivatives of the Twiss parameters w.r.t delta_p
    if (deltap_dependency) then
       if( (c_%npara==5) .or. (c_%ndpt/=0) ) then ! condition to be checked
          call computeDeltapDependency(A_script,s1)
       endif
    endif
    ! ---

    !swap for longitudinal beta with gamma
    if (c_%nd == 3) then
       do i=1,c_%nd
          test = s1%beta(3,i)
          s1%beta(3,i) = s1%gama(3,i)
          s1%gama(3,i) = test
       enddo
    endif

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) ' check_stable ',check_stable,' c_%stable_da ',c_%stable_da,' PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('equaltwiss CHECK 5 : ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"equaltwiss CHECK 5 ",whymsg)
      return
    endif


    s1%mu=phase

    do k=1,3
       do i=1,6
          s1%eigen(k*2-1,i) = A_script(k*2-1).sub.fo(i,:)
          s1%eigen(k*2  ,i) = A_script(k*2  ).sub.fo(i,:)
       enddo
    enddo

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
       write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost(:len_trim(messagelost))
       call fort_warn('e qualtwiss: ',whymsg(:len_trim(whymsg)))
       call seterrorflag(10,"equaltwiss: ",whymsg);
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
    use s_extend_poly, only : MAPDUMP ! LD: 04.06.2019
    implicit none
    logical(lp)             :: closed_orbit,beta_flg, slice, goslice
    integer                 :: k,i,ii
    integer                 :: mynd2,npara,nda,icase,flag_index,why(9),my_nv,nv_min
    integer                 :: ioptfun,iii,restart_sequ,advance_node,mf1,mf2
    integer                 :: tab_name(*)
    integer                 :: summary_tab_name(*)
    real(dp)                :: deltap0,deltap ,d_val
    double precision         :: get_value,suml,s
    integer                 :: posstart, posnow
    integer                 :: geterrorflag !C function that returns errorflag value
    real(dp)                :: orbit(6)=0.d0
    type(probe)             :: orbit_probe
    type(probe_8)           :: A_script_probe ! A_script == A**(-1); oneTurnMap = A_script o R o A
    type(twiss)             :: tw
    type(fibre), POINTER    :: current
    type(integration_node), pointer :: nodePtr, stopNode
    type(work)              :: startfen !Fibre energy at the start
    real(dp)                :: relativisticBeta  !current beta0, set by getdeltae()
    real(dp)                :: r,re(ndim2,ndim2),dt
    logical(lp)             :: initial_matrix_manual, initial_matrix_table, initial_map_manual
    logical(lp)             :: initial_distrib_manual, initial_ascript_manual, writetmap
    logical(lp)             :: maptable
    logical(lp)             :: ring_parameters  !! forces isRing variable to true, i.e. calclulation of closed solution
    logical(lp)             :: doNormal         !! do normal form analysis
    logical(lp)             :: doRDTtracking    !!
    logical(lp)             :: isstochastic  !! tempurary veriable used in switching off stochastic in closed orbit search
    type(c_damap)           :: AscriptInPhasor, dummyMap  !! maps for RDTs calculations
    type(c_vector_field)    :: vectorField                !! defined here to avoid every step alloc and kill
    type(c_taylor)          :: theRDTs                    !!
    real(dp)                :: emi(3)
    logical(lp)             :: isputdata  ! in everystep mode (node by node) switch deciding if data are to be put in twiss table for a give node
    logical(lp)             :: rmatrix  ! flag to mark that transfer matrix should be saved (otherwise we might not track theTransferMap)
    logical(lp)             :: isTMsave ! flag that TM was recorded during pre-run for closed solution search
    logical(lp)             :: doTMtrack ! true if rmatrix==true .and. isRing==true . do not track theTransferMap and save time
                                           !      .or. already tracked form closed solution search
    logical(lp)             :: usertableActive = .false.  ! flag to mark that there was something requested with ptc_select
    integer                 :: countSkipped
    character(48)           :: summary_table_name
    character(12)           :: tmfile='transfer.map'
    character(48)           :: charconv !routine
    real(dp)                :: BETA0
    integer                 :: mapdumpbak ! LD: 04.06.2019

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

    if (getdebug() > 1) then
      if(.not.associated(MY_RING%t)) then
        print*,"ptc_twiss: NODE LAYOUT ALREADY CREATED"
      else
        print*,"ptc_twiss: NODE LAYOUT NOT YET CREATED"
      endif
    endif

    call resetBetaExtremas()
    call initIaaMatrix()

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

    rmatrix = get_value('ptc_twiss ','rmatrix ') .ne. 0
    if ( (getnknobsall() + getnpushes()) < 1) then
      usertableActive = .false.
    else
      usertableActive = .true.
      rmatrix = .true. ! for the time being force if ptc_select or ptc_knob was defined
    endif

    maptable = get_value('ptc_twiss ','maptable ') .ne. 0
    if (maptable) then
      rmatrix = .true. 
    endif
    
    isTMsave = .false.

    no = get_value('ptc_twiss ','no ')
    if ( no .lt. 1 ) then
       call fort_warn('madx_ptc_twiss.f90 <ptc_twiss>:','Order in twiss is smaller then 1')
       print*, "Order is ", no
       call tidy()
       return
    endif

    icase = get_value('ptc_twiss ','icase ')

    deltap0 = get_value('ptc_twiss ','deltap ')

    deltap = zero

    call my_state(icase,deltap,deltap0)
    if (getdebug() > 2) then
       print *, "ptc_twiss: internal state after my_state:"
       call print(default,6)
    endif

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

    call make_node_layout(my_ring)
    call getBeamBeam()


    !############################################################################
    !############################################################################
    !############################################################################


    orbit(:)=zero
    ! read the orbit
    ! if closed orbit is to be found pass it as the starting point for the searcher
    orbit(1)=get_value('ptc_twiss ','x ')
    orbit(2)=get_value('ptc_twiss ','px ')
    orbit(3)=get_value('ptc_twiss ','y ')
    orbit(4)=get_value('ptc_twiss ','py ')
    orbit(6)=-get_value('ptc_twiss ','t ') ! swap of t sign
    orbit(5)=get_value('ptc_twiss ','pt ')

    if(mytime) then
       call Convert_dp_to_dt (deltap, dt)
    else
       dt=deltap
    endif
    if(icase.eq.5 .or. icase.eq.56) orbit(5) = dt + orbit(5)

    closed_orbit = get_value('ptc_twiss ','closed_orbit ') .ne. 0

    if( closed_orbit .and. (icav .gt. 0) .and. (my_ring%closed .eqv. .false.)) then
       call fort_warn('return from ptc_twiss: ',' Closed orbit requested on not closed layout.')
       call seterrorflag(3,"ptc_twiss ","Closed orbit requested on not closed layout.")
       call tidy()
       return
    endif

    if(closed_orbit) then

       if ( .not. c_%stable_da) then
          call fort_warn('ptc_twiss: ','DA got unstable even before finding closed orbit')
          call seterrorflag(10,"ptc_twiss ","DA got unstable even before finding closed orbit")
          call aafail('ptc_twiss: ','DA got unstable even before finding closed orbit. program stops')
          !          return
       endif


       if (getdebug() > 2) then
         print*, "Looking for orbit"
         print*, "Init orbit ", orbit
         call print(default,6)
       endif

       ! disable stochastic for closed orbit seach
       isstochastic = default%stochastic
       default%stochastic = .false.

       current=>my_ring%start
       !global_verbose = .true.
       call FIND_ORBIT_x(orbit,default,c_1d_8,fibre1=current)
       !global_verbose = .false.

       default%stochastic = isstochastic

       if ( .not. check_stable) then
          write(whymsg,*) 'DA got unstable during closed orbit search: PTC msg: ',messagelost(:len_trim(messagelost))
          call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
          call seterrorflag(10,"ptc_twiss ",whymsg);
          call tidy()
          return
       endif

      if (getdebug() > 0) then
         CALL write_closed_orbit(icase,orbit)
      endif

    elseif((my_ring%closed .eqv. .true.) .and. (getdebug() > 1)) then
       print*, "Closed orbit specified by the user!"
       !CALL write_closed_orbit(icase,x) at this position it isn't read
    endif


    orbit_probe = orbit

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

    n_rf = nclocks

    ! old PTC
    !call init(default,no,nda,BERZ,mynd2,npara)

    !new complex PTC
    call init_all(default,no,nda,BERZ,mynd2,npara,nclocks) ! need to add number of clocks
    c_verbose=.false.

    i_piotr(:) = 0
    i_piotr(1)=1; i_piotr(2)=3; i_piotr(3)=5;

    c_normal_auto=.true.;


!    call init_all(default,no,nda)
    ! mynd2 and npara are outputs

    if (getdebug() > 2) then
       print *, "ptc_twiss: internal state after init:"
       call print(default,6)
       print*, "no,nda,BERZ,mynd2,npara"
       print*, no,nda,BERZ,mynd2,npara
       print*,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
       print*,""
    endif


    if ( (cavsareset .eqv. .false.) .and. (my_ring%closed .eqv. .false.) ) then

       call setcavities(my_ring,maxaccel)
       if (geterrorflag() /= 0) then
          call tidy()
          return
       endif
    endif

    call setknobs(my_ring)



    call alloc(A_script_probe)
    A_script_probe%u=my_false
    A_script_probe%x=npara
    A_script_probe%x=orbit


    !This must be before init map
    call alloc(theTransferMap)
    theTransferMap%u = .false.
    theTransferMap%x = npara
    theTransferMap%x = orbit

    !############################################################################
    !############################################################################
    !############################################################################
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  INIT A_script_probe that is tracked          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call initmap(dt,slice)

    if (geterrorflag() /= 0) then
       !if arror occured then return
       call tidy()
       return
    endif

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) 'DA got unstable during A_ initialization: The closed solution does not exist. PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"ptc_twiss INIT CHECK",whymsg)
      call tidy()
      return
    endif

    !############################################################################
    !############################################################################
    !############################################################################

    ! n_rf is a variable of PTC that must be set accordingly


    do i=1,nclocks
      A_script_probe%ac(i)%om = twopi*clocks(i)%tune * BETA0start / my_ring_length
      !omega of the the modulation
      !savedProbe%ac%om = twopi*clocks(1)%tune
      A_script_probe%ac(i)%x(1)  = one
      A_script_probe%ac(i)%x(2)  = zero  ! initial clock vector (sin like)
    enddo




       doRDTtracking = get_value('ptc_twiss ','trackrdts ') .ne. 0
       if (doRDTtracking) then
          call alloc(theRDTs)
          call alloc(vectorField)
          call alloc(AscriptInPhasor)
          call alloc(dummyMap)
       endif




    !############################################################################
    !############################################################################
    !############################################################################



    ! assume that we track the transfer map
    doTMtrack = .true.

    !and check now if we could skip it and save time
    if ((usertableActive .eqv. .false.) ) then
      !we do not need the full thing to fill usertable
      ! currently we save only the first order map (i.e. matrix)
      ! we can consider saving the map in universal taylor, if this still pays back

      if (getdebug() > 1) then
        print*, "There is no user tables (ptc_select or ptc_knobs), we can eventually skip the TM tracking"
        print*, "FLAGS: isTMsave  slice  rmatrix isRing"
        print*, isTMsave,  slice , rmatrix, isRing
        ! isTMsave - tells if Transfer Matrix was saved; true if A_ initized from closed solution
        ! slice -
        ! rmatrix - user requests rmatrix in the twiss table
        ! isRing - request closed solution to twiss table
      endif

      if ( isTMsave .and.  (slice .eqv. .false.) ) then
        !the matrix was tracked in the initmap, and slice is off (it is saved only for each element)
        ! in the future should do support for slicing
        doTMtrack = .false.
      endif


      !another independent condition
      if ( (isRing .eqv. .false.) .and. (rmatrix .eqv. .false.) ) then
        !we do not have to do the normal form at the end (isRing false)
        ! and
        !the user does not want transfer maps (rmatrix false)
        !
        doTMtrack = .false.
      endif


     endif

    if (doTMtrack) then
      !we will be tracking theTransferMap
      ! get clean initialization after initmap
      call kill(theTransferMap)
      call alloc(theTransferMap)
      theTransferMap%u = .false.
      theTransferMap%x = npara
      theTransferMap%x = orbit

      if (getdebug() > 1) then
        print*, "doTMtrack=true, theTransferMap is new"
      endif

    elseif (getdebug() > 1) then
        print*, "doTMtrack=false, theTransferMap stays as it was"

    endif


    !############################################################################
    !############################################################################
    !############################################################################



    call alloc(tw)

    !Y

    !the initial twiss is needed to initialize propely calculation of some variables f.g. phase advance
    tw = A_script_probe%x
    if (geterrorflag() /= 0) then
       call fort_warn('ptc_twiss: ','equaltwiss at the begining of the line returned with error')
       call tidy()
       return
    endif

    phase = zero !we have to do it after the very initial twiss params calculation above
    current=>MY_RING%start
    startfen = 0
!    print*,'Check 5'
!    call print(default,6)
!    print*,'my_state'
    startfen = current!setting up start energy for record
!    print*,'Check 6'
!    call print(default,6)
!    print*,'my_state'
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

    if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
      write(whymsg,*) 'DA got unstable during initialization: The closed solution does not exist. PTC msg: ', &
                       messagelost(:len_trim(messagelost))
      call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
      call seterrorflag(10,"ptc_twiss INIT CHECK",whymsg)
      call tidy()
      return
    endif

    allocate(j(c_%npara)) ! used in puttwisstable for defining monomials to peek


    if (geterrorflag() /= 0) then
       call fort_warn('ptc_twiss: ','equaltwiss at the beginning of the line returned with error')
       call tidy()
       return
    endif



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                              !
   !  T H E   M A I N   L O O P   !
   !                              !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,MY_RING%n

      if (getdebug() > 2) then
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
          !if ((s .eq. 0d0) .and. (nodePtr%pos .eq. (my_ring%t%n+posstart-1))) then
          !   s = nodePtr%s(1) + nodePtr%next%s(5) ! s of previous node + local offset
          !endif
          s = nodePtr%s(1)
          if (getdebug() > 2) then
             write(6,*) "##### SLICE MAGNETS NODE ",&
                      & nodePtr%pos," => ",nodePtr%pos+1," s=",s
          endif

          ! CALL PROPAGATE WITH PROPER OPTIONS
          call propagateswy()


          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then

             write(whymsg,*) 'DA got unstable in tracking at s= ',s, &
                             ' magnet ',i,' ', current%mag%name,' ', current%mag%vorname, &
	         ' step ',nodePtr%pos,' PTC msg: ',messagelost(:len_trim(messagelost))
             call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
             call seterrorflag(10,"ptc_twiss ",whymsg);

             if (getdebug() > 2) close(mf1)
             call tidy()
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

            tw = A_script_probe%x ! set the twiss parameters, with y being equal to the A_ phase advance
            suml = s;

            call puttwisstable(theTransferMap%x)
            if(doRDTtracking)   call putrdttable(current)
            if(usertableActive) call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap%x, A_script_probe%x)

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

          tw = A_script_probe%x ! set the twiss parameters
          if (s > suml) then !work around against last element having s=0
            suml = s;
          endif

          call puttwisstable(theTransferMap%x)
          if(doRDTtracking)   call putrdttable(current)
          if(usertableActive) call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap%x, A_script_probe%x)

        endif

      else
        ! ELEMENT AT ONCE MODE
        if (nda > 0) then
           if (mapdump .eq. 0 .or. mapdump .ge. 11) then ! mapdump = 0,11,12
             mapdumpbak = mapdump ; mapdump = modulo(mapdump, 10)
             call propagate(my_ring,A_script_probe,+default,fibre1=i,fibre2=i+1)
             mapdump = mapdumpbak
           endif
           if (doTMtrack .and. mapdump .ge. 0 .and. mapdump .le. 2) then ! mapdump = 0,1,2
             call propagate(my_ring,theTransferMap,+default,fibre1=i,fibre2=i+1)
           endif
        else
           if (mapdump .eq. 0 .or. mapdump .ge. 11) then ! mapdump = 0,11,12
             mapdumpbak = mapdump ; mapdump = modulo(mapdump, 10)
             call propagate(my_ring,A_script_probe,default, fibre1=i,fibre2=i+1)
             mapdump = mapdumpbak
           endif
           if (doTMtrack .and. mapdump .ge. 0 .and. mapdump .le. 2) then ! mapdump = 0,1,2
             call propagate(my_ring,theTransferMap,default,fibre1=i,fibre2=i+1)
           endif

        endif


        if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then

           write(whymsg,*) 'DA got unstable in tracking at s= ',s, &
                           ' magnet ',i,' ', current%mag%name,' ', current%mag%vorname, &
	       ' PTC msg: ',messagelost(:len_trim(messagelost))
           call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
           call seterrorflag(10,"ptc_twiss ",whymsg);
           if (getdebug() > 2) close(mf1)
           call tidy()
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
           call print(A_script_probe,mf1)

           if (current%mag%kind==kind4) then
             print*,"CAVITY at s=",suml,"freq ", current%mag%freq,&
                    	"lag ", current%mag%lag, &
                    	"volt ", current%mag%volt
            endif
        endif

        suml=suml+current%MAG%P%ld

        if (savemaps) then
           do ii=1,6
              maps(i)%unimap(ii) = A_script_probe%x(ii)
           enddo
           maps(i)%s = suml
           maps(i)%name = current%mag%name
        endif

        !print*,"Skowron 4 ", current%mag%name,  check_stable, c_%stable_da, A_script_probe%x(1).sub.'100000'

        ! compute the Twiss parameters
        tw = A_script_probe%x
        if (geterrorflag() /= 0) then
           call fort_warn('ptc_twiss: ','equaltwiss at ' // current%mag%name // ' returned with error')
           call tidy()
           return
        endif

        !print*,"Skowron 5 ", current%mag%name,  check_stable, c_%stable_da, A_script_probe%x(1).sub.'100000'

        if(isTMsave) then
          call puttwisstable(theTransferMap%x,transfermapSaved=savedTM(i,:,:))
        else
          call puttwisstable(theTransferMap%x)
        endif

        !print*,"Skowron 6 ", current%mag%name,  check_stable, c_%stable_da, A_script_probe%x(1).sub.'100000'

        if(doRDTtracking)   call putrdttable(current)
        if(usertableActive) call putusertable(i,current%mag%name,suml,getdeltae(),theTransferMap%x,A_script_probe%x)

        !print*,"Skowron 7 ", current%mag%name,  check_stable, c_%stable_da, A_script_probe%x(1).sub.'100000'

      endif

      if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
        write(whymsg,*) 'DA got unstable in ', current%mag%name, &
	    ' check_stable ',check_stable, ' c_%stable_da ',c_%stable_da, &
	    ' PTC msg: ',  messagelost(:len_trim(messagelost))
        call fort_warn('ptc_twiss: ',whymsg(:len_trim(whymsg)))
        call seterrorflag(10,"ptc_twiss end of loop check",whymsg)
        call tidy()
        return
      endif


      iii=advance_node()
      current=>current%next
    enddo

100 continue

   ! print*,"Skowron total tracking time", tsum

    ! relocated the following here to avoid side-effect
    print77=.false.
    read77=.false.


    if (writetmap) then
       !'===      TRANSFER MAP      ==='
       call kanalnummer(mf2)
       open(unit=mf2,file=tmfile)
       call print(theTransferMap%x,mf2)
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

       call print(A_script_probe,mf2)

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

    ! Normal
    doNormal = get_value('ptc_twiss ','normal ') .ne. 0

    if(isRing .eqv. .true.) then
       if (doNormal) call normalFormAnalysis(theTransferMap ,A_script_probe, orbit)
       call oneTurnSummary(theTransferMap ,A_script_probe%x, orbit, suml)
    else
      print*, "Reduced SUMM Table (Inital parameters specified)"
      call onePassSummary(theTransferMap%x , orbit, suml)
    endif


    if(maptable) then
       call makemaptable(theTransferMap%x,no)
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


    call finishknobs()

    if (savemaps) then  !do it at the end, so we are sure the twiss was successful
       mapsorder = no
       mapsicase = icase
       if (getnmoments() > 0) call ptc_moments(no*2) !calcualate moments with the maximum available order
    endif


    call tidy()

! f90flush is not portable, and useless...
!    call f90flush(20,my_false)

    if (getdebug() > 2) close(mf1)


    !****************************************************************************************
    !*********  E N D   O F   PTC_TWISS      ************************************************
    !****************************************************************************************
    !________________________________________________________________________________________

  contains  ! what follows are internal subroutines of ptc_twiss
    !____________________________________________________________________________________________

    subroutine propagateswy()
      implicit none

       if (nda > 0) then
          if (mapdump .eq. 0 .or. mapdump .ge. 11) then ! mapdump = 0,11,12
            mapdumpbak = mapdump ; mapdump = modulo(mapdump, 10)
            call propagate(my_ring,A_script_probe,+default,node1=nodePtr%pos,node2=nodePtr%pos+1)
            mapdump = mapdumpbak
          endif
          if (doTMtrack .and. mapdump .ge. 0 .and. mapdump .le. 2) then ! mapdump = 0,1,2
            call propagate(my_ring,theTransferMap,+default,node1=nodePtr%pos,node2=nodePtr%pos+1)
          endif
        else
          if (mapdump .eq. 0 .or. mapdump .ge. 11) then ! mapdump = 0,11,12
            mapdumpbak = mapdump ; mapdump = modulo(mapdump, 10)
            call propagate(my_ring,A_script_probe,default,node1=nodePtr%pos,node2=nodePtr%pos+1)
            mapdump = mapdumpbak
          endif
          if (doTMtrack .and. mapdump .ge. 0 .and. mapdump .le. 2) then ! mapdump = 0,1,2
            call propagate(my_ring,theTransferMap,default,node1=nodePtr%pos,node2=nodePtr%pos+1)
          endif
        endif

    end subroutine propagateswy

    !____________________________________________________________________________________________
    subroutine tidy()
      ! deallocates all the variables
      implicit none

      call kill(tw)
      CALL kill(A_script_probe)

      CALL kill(theTransferMap)

      do i=1,6
         call kill(unimap(i))
      enddo


      if (allocated(j)) deallocate(j)

      if (allocated(savedTM)) deallocate(savedTM)

    end subroutine tidy

    !____________________________________________________________________________________________
    subroutine initmap(dt,slice)
      implicit none
      integer     :: double_from_table_row
      integer     :: mman, mtab, mascr, mdistr !these variable allow to check if the user did not put too many options
      integer     :: mmap
      real(dp)    :: dt
      logical(lp) :: slice
      integer     :: mf
      integer     :: i,j,k

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
            A_script_probe%x(i) = unimap(i)
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
            call print(theTransferMap,6)

            print*,"Tracking identity map to get closed solution. STATE:"
            call print(default,6)

         endif

         if (getdebug() > 2) then
           print*, "printing the initial map"
           call print(theTransferMap,17)
         endif

         allocate(savedTM(my_ring%n,6,6))


         do i=1,MY_RING%n

           call propagate(my_ring,theTransferMap,default, fibre1=i,fibre2=i+1)

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
              return
           endif

           do j=1,6
             do k=1,6
               savedTM(i,j,k) = theTransferMap%x(j).sub.fo(k,:)
             enddo
           enddo

         enddo

         if (getdebug() > 2) then
            call print(theTransferMap,18)
            call kanalnummer(mf,file="theOneTurnMap.txt")
            call print(theTransferMap,mf)
            close(mf)

            print*,"Initializing map from One Turn Map."
         endif

        isTMsave = .true.

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

      relativisticBeta = cfen%beta0

      if (getdebug() > 2) then
         write(mf1,'(3(a, f12.6))') "Ref Momentum ",cfen%p0c," Energy ", cfen%energy," DeltaE ",getdeltae
      endif

    end function getdeltae
    !____________________________________________________________________________________________

    subroutine puttwisstable(transfermap,transfermapSaved)
      implicit none
      include "madx_ptc_knobs.inc"
      type(real_8), target :: transfermap(6)  !
      double precision , optional :: transfermapSaved(6,6)
      integer i1,i2,ii,i1a,i2a
      double precision  :: opt_fun(150), tmpa6(6),tmpa66(6,6) ,myx  ! opt_fun(72) -> opt_fun(81)
      ! increase to 150 to have extra space beyond what's needed to accomodate additional derivatives w.r.t. delta_p
      double precision    :: deltae ! for reference energy increase via acceleration
      double precision    :: deltap ! for deltap treatment
      ! added on 3 November 2010 to hold Edwards & Teng parametrization
      real(dp) :: betx,bety,alfx,alfy,R11,R12,R21,R22, pt_, onedp
      ! to convert between Ripken and Edwards-Teng parametrization
      real(dp) :: kappa,u,ax,ay,kx,ky,kxy2,usqrt,bx,by,cx,cy,cosvp,sinvp,cosvm,sinvm,cosv2,sinv2,cosv1,sinv1
      real(dp) :: deltaeValue

      if (getdebug() > 2) then
         write(mf1,*) "##########################################"
         write(mf1,*) ""
         write(mf1,'(i4, 1x,a, f10.6)') i,current%mag%name,suml
         write(mf1,*) ""
         call print(A_script_probe,mf1)
      endif


      deltae = getdeltae()

      call double_to_table_curr(table_name, 's ', suml)

      doublenum = deltae * startfen%energy
      call double_to_table_curr(table_name, 'energy ', doublenum)


      opt_fun(:)=zero

      call liepeek(iia,icoast)
      j(:)=0

      do ii=1,c_%npara ! fish
         opt_fun(ii)=A_script_probe%x(ii).sub.j
      enddo

      call trackOrbitExtremaAndRms(opt_fun(1:6))

      ! swap t and pt
      myx=opt_fun(6)
      opt_fun(6)=opt_fun(5)
      opt_fun(5)=-myx


      ioptfun=6
      call vector_to_table_curr(table_name, 'x ', opt_fun(1), ioptfun)
      opt_fun(:)=zero ! reset it

      if (rmatrix) then
        if (present(transfermapSaved)) then
          ! we have to swap 5 and 6, and to avoid confusion, I do it on a copy (original would be swapped, and if used again...)

          tmpa66 = transfermapSaved

          tmpa6 = tmpa66(:,6)
          tmpa66(:,6) = tmpa66(:,5)
          tmpa66(:,5) = tmpa6


          opt_fun( 1:6 ) = tmpa66(1,:)
          opt_fun( 7:12) = tmpa66(2,:)
          opt_fun(13:18) = tmpa66(3,:)
          opt_fun(19:24) = tmpa66(4,:)
          opt_fun(31:36) = tmpa66(5,:)
          opt_fun(25:30) = tmpa66(6,:)

        else

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

        endif

        ioptfun=36
        call vector_to_table_curr(table_name, 're11 ', opt_fun(1), ioptfun)
        opt_fun(:)=zero ! reset it

      endif

      !deltap = A_script_probe%x(5).sub.'0'
      !deltae = deltae * (1.0 + deltap)
      if(default%time) then
        pt_ = A_script_probe%x(5).sub.'0'
        onedp   = sqrt( one + two*pt_/relativisticBeta + (pt_**2))
      else
        onedp = one + A_script_probe%x(5).sub.'0'
      endif
      


      opt_fun(beta11)= tw%beta(1,1) * onedp ! beta11=1
      opt_fun(beta12)= tw%beta(1,2) * onedp
      opt_fun(beta13)= tw%beta(1,3) * onedp
      opt_fun(beta21)= tw%beta(2,1) * onedp
      opt_fun(beta22)= tw%beta(2,2) * onedp
      opt_fun(beta23)= tw%beta(2,3) * onedp
      opt_fun(beta31)= tw%beta(3,1) * onedp
      opt_fun(beta32)= tw%beta(3,2) * onedp
      opt_fun(beta33)= tw%beta(3,3) * onedp

      opt_fun(alfa11)= tw%alfa(1,1) 
      opt_fun(alfa12)= tw%alfa(1,2) 
      opt_fun(alfa13)= tw%alfa(1,3) 
      opt_fun(alfa21)= tw%alfa(2,1) 
      opt_fun(alfa22)= tw%alfa(2,2) 
      opt_fun(alfa23)= tw%alfa(2,3) 
      opt_fun(alfa31)= tw%alfa(3,1) 
      opt_fun(alfa32)= tw%alfa(3,2) 
      opt_fun(alfa33)= tw%alfa(3,3) 

      opt_fun(gama11)= tw%gama(1,1) / onedp
      opt_fun(gama12)= tw%gama(1,2) / onedp
      opt_fun(gama13)= tw%gama(1,3) / onedp
      opt_fun(gama21)= tw%gama(2,1) / onedp
      opt_fun(gama22)= tw%gama(2,2) / onedp
      opt_fun(gama23)= tw%gama(2,3) / onedp
      opt_fun(gama31)= tw%gama(3,1) / onedp
      opt_fun(gama32)= tw%gama(3,2) / onedp
      opt_fun(gama33)= tw%gama(3,3) / onedp


      ! --- derivatives of Twiss paramters w.r.t delta_p
      ! NOW why do we need to multiply by onedp, as for the other Twiss parameters?
      if (deltap_dependency) then
         opt_fun(beta11p)= tw%beta_p(1,1) * onedp
         opt_fun(beta12p)= tw%beta_p(1,2) * onedp
         opt_fun(beta13p)= tw%beta_p(1,3) * onedp
         opt_fun(beta21p)= tw%beta_p(2,1) * onedp
         opt_fun(beta22p)= tw%beta_p(2,2) * onedp
         opt_fun(beta23p)= tw%beta_p(2,3) * onedp
         opt_fun(beta32p)= tw%beta_p(3,2) * onedp
         opt_fun(beta33p)= tw%beta_p(3,3) * onedp

         opt_fun(alfa11p)= tw%alfa_p(1,1) 
         opt_fun(alfa12p)= tw%alfa_p(1,2) 
         opt_fun(alfa13p)= tw%alfa_p(1,3) 
         opt_fun(alfa21p)= tw%alfa_p(2,1) 
         opt_fun(alfa22p)= tw%alfa_p(2,2) 
         opt_fun(alfa23p)= tw%alfa_p(2,3) 
         opt_fun(alfa31p)= tw%alfa_p(3,1) 
         opt_fun(alfa32p)= tw%alfa_p(3,2) 
         opt_fun(alfa33p)= tw%alfa_p(3,3) 

         opt_fun(gama11p)= tw%gama_p(1,1) / onedp
         opt_fun(gama12p)= tw%gama_p(1,2) / onedp
         opt_fun(gama13p)= tw%gama_p(1,3) / onedp
         opt_fun(gama21p)= tw%gama_p(2,1) / onedp
         opt_fun(gama22p)= tw%gama_p(2,2) / onedp
         opt_fun(gama23p)= tw%gama_p(2,3) / onedp
         opt_fun(gama31p)= tw%gama_p(3,1) / onedp
         opt_fun(gama32p)= tw%gama_p(3,2) / onedp
         opt_fun(gama33p)= tw%gama_p(3,3) / onedp
      endif
      ! --- end

      ! march 10th: do we need to multiply by onedp the following?
      opt_fun(mu1)=tw%mu(1) !* deltae
      opt_fun(mu2)=tw%mu(2) !* deltae
      opt_fun(mu3)=tw%mu(3) !* deltae

      ! write(0,*),"DEBUG = ", tw%mu(3) ! should give Qs
      ! LD: 2019.12.18 removed after complains of FS
!      if ( default%time ) then
!        tw%disp(:) = tw%disp(:) * relativisticBeta
!      endif

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


      !skowron: 2017 Jan, twiss table has sigma matrix in this place
      !
      !     opt_fun(62+4+8+6)=zero ! was 36 instead of 62 => on 9 march add 4 => on 3 July 2009 add 4 => on 3 November 2010 add 6
      !     do i1=1,6
      !        if(i1.le.4) then
      !           i1a=i1
      !        elseif(i1.eq.5) then
      !           i1a=6
      !        else
      !           i1a=5
      !        endif
      !        do i2=1,6
      !          if(i2.le.4) then
      !	 i2a=i2
      !           elseif(i2.eq.5) then
      !	 i2a=6
      !           else
      !	 i2a=5
      !           endif
      !           ! idiotic counting reset (62+4+8+6) ==> 73
      !           ii=73+(i1a-1)*6+i2a ! was 36 instead of 62 => now (62+4) instead of 62 => 62+4+4
      !
      !           opt_fun(ii)=tw%eigen(i1,i2) * deltae
      !           ! where do these eigen values go? no such dedicated column defined in madx_ptc_knobs.inc
      !           if(mytime.and.i2a.eq.6) opt_fun(ii)=-opt_fun(ii)
      !
      !        enddo
      !     enddo


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
         write(6,'(a,4(f11.4,1x))')  "orbit transv.", A_script_probe%x(1).sub.'0',A_script_probe%x(2).sub.'0', &
                                                      A_script_probe%x(3).sub.'0',A_script_probe%x(4).sub.'0'
         write(6,'(a,2(f11.4,1x))')  "dp/p, T      ", A_script_probe%x(5).sub.'0',A_script_probe%x(6).sub.'0'

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
      opt_fun(:)=zero

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
         betx = tw%beta(1,1) * onedp
         bety = tw%beta(2,2) * onedp
         alfx = tw%alfa(1,1)
         alfy = tw%alfa(2,2)

      else

         kx=sqrt(tw%beta(1,2)/tw%beta(1,1)); ! multiplication by deltae in numerator and denominator
         ky=sqrt(tw%beta(2,1)/tw%beta(2,2));

         ! beta11, alfa11 etc... are multiplied by deltae before output
         ax=kx*tw%alfa(1,1) * onedp -tw%alfa(1,2) * onedp /kx;
         ! hence we reflect this in the formula from Lebedev
         ay=ky*tw%alfa(2,2) * onedp -tw%alfa(2,1) * onedp /ky;
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

            betx = (tw%beta(1,1)/kappa) * onedp
            bety = (tw%beta(2,2)/kappa) * onedp
            alfx = (tw%alfa(1,1)/kappa) 
            alfy = (tw%alfa(2,2)/kappa) 

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

            betx = tw%beta(1,1) * onedp
            bety = tw%beta(2,2) * onedp
            alfx = tw%alfa(1,1)
            alfy = tw%alfa(2,2)
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
      call trackBetaExtrema(tw%beta,onedp,betx,bety,tw%disp)
      !---

    end subroutine puttwisstable
    !____________________________________________________________________________________________

    subroutine putrdttable(fib)
      implicit none
      type(fibre), POINTER    :: fib
      complex(dp)   :: c_val
      real(dp)    :: im_val, re_val, d_val,  eps=1e-6
      integer     :: ind(10), i, mynres, order,rrr
      character(len=18):: nick

        dummyMap=A_script_probe%x
        call c_canonise(dummyMap,AscriptInPhasor)

        AscriptInPhasor=to_phasor() * AscriptInPhasor * from_phasor()
        call c_factor_map(AscriptInPhasor,dummyMap,vectorField,0)

        theRDTs = cgetpb(vectorField)

        call string_to_table_curr(rdt_table_name,"name ","name ")
        call double_to_table_curr(rdt_table_name, 's ', suml)


        call c_taylor_cycle(theRDTs,size=mynres)

        do rrr=1,mynres

            ind = 0
            call c_taylor_cycle(theRDTs,ii=rrr,value=c_val,j=ind(1:c_%nv))

            order = sum(ind(1:6))


            !print*,"GNFU ",ind(1:6)

            im_val = imag(c_val)
            re_val = real(c_val)
            d_val  = hypot(re_val, im_val)

            ! if amplitude is close to zero then it is not worth to output
            if (d_val .lt. eps) then
              if (getdebug()>2) print*,"putGnormaltable idx=",rrr," ",d_val," smaller then eps=",eps, " skipping "
              cycle
            endif

           write(nick,'(a4,6(a1,i1))') 'gnfa','_',ind(1),'_',ind(2),'_',ind(3), &
                    	      '_',ind(4),'_',ind(5),'_',ind(6)
           call double_to_table_curr2(rdt_table_name, nick, d_val )

           nick(4:4) = 'c'
           call double_to_table_curr2(rdt_table_name, nick, re_val )
           nick(4:4) = 's'
           call double_to_table_curr2(rdt_table_name, nick, im_val )


          ! write(*,*) nick, " = ", d_val

       enddo

       if (fib%mag%p%nmul > 0) then
       endif

       if (fib%mag%p%nmul > 1) then
         if (fib%mag%l > 0) then
           call double_to_table_curr2(rdt_table_name,'k1l ', fib%mag%bn(2)*fib%mag%l)
           call double_to_table_curr2(rdt_table_name,'k1sl ',fib%mag%an(2)*fib%mag%l)
         else
           call double_to_table_curr2(rdt_table_name,'k1l ', fib%mag%bn(2))
           call double_to_table_curr2(rdt_table_name,'k1sl ',fib%mag%an(2))
         endif
       endif

       if (fib%mag%p%nmul > 2) then
         if (fib%mag%l > 0) then
           call double_to_table_curr(rdt_table_name,'k2l ', fib%mag%bn(3)*fib%mag%l)
           call double_to_table_curr(rdt_table_name,'k2sl ',fib%mag%an(3)*fib%mag%l)
         else
           call double_to_table_curr(rdt_table_name,'k2l ', fib%mag%bn(3))
           call double_to_table_curr(rdt_table_name,'k2sl ',fib%mag%an(3))
         endif

       endif

       if (fib%mag%p%nmul > 3) then
         if (fib%mag%l > 0) then
           call double_to_table_curr(rdt_table_name,'k3l ', fib%mag%bn(4)*fib%mag%l)
           call double_to_table_curr(rdt_table_name,'k3sl ',fib%mag%an(4)*fib%mag%l)
         else
           call double_to_table_curr(rdt_table_name,'k3l ', fib%mag%bn(4))
           call double_to_table_curr(rdt_table_name,'k3sl ',fib%mag%an(4))
         endif
       endif

       call augment_count(rdt_table_name)
       ! write(*,*)

    end subroutine putrdttable

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
      orbit(:)=zero
      orbit(1)=get_value('ptc_twiss ','x ')
      orbit(2)=get_value('ptc_twiss ','px ')
      orbit(3)=get_value('ptc_twiss ','y ')
      orbit(4)=get_value('ptc_twiss ','py ')
      orbit(5)=get_value('ptc_twiss ','pt ')
      orbit(6)=-get_value('ptc_twiss ','t ')

      if(icase.eq.5 .or. icase.eq.56 ) orbit(5) = orbit(5) + dt

      orbit_probe = orbit
    end subroutine readreforbit
    !_________________________________________________________________
    logical(dp) function checksymmetric(r)
     implicit none
     real(dp)                :: r(6,6)
     integer i,j
     character(1024)  :: msg
     checksymmetric = .false.
      print*,"++++++++++++++"
      print*, r(1,:)
      print*, r(2,:)
      print*, r(3,:)
      print*, r(4,:)
      print*, r(5,:)
      print*, r(6,:)

     do i=1,6
       do j=i,6
         if ( abs( r(i,j) - r(j,i) ) > 1e-24 ) then
           print*,r(i,j),  r(j,i)
           write(msg,'(a,i1,a,i1,a,G10.2,a,i1,a,i1,a,G10.2)') &
              "checksymmetric re matrix re(",i,",",j,")=",r(i,j),&
              "   not equal to re("            ,j,",",i,")=",r(j,i)
           call fort_warn('ptc_twiss: ',msg(:len_trim(msg)))
           checksymmetric = .true.
           return
         endif
       enddo
     enddo


    end function checksymmetric

    subroutine readinitialdistrib
      !reads covariance matrix of the initial distribution
      implicit none
      include 'madx_ptc_distrib.inc'
      type(taylor) ht, ht0
      type(pbfield) h
      type(damap) id
      type(normalform) norm
      real(dp) lam
      integer nd,nd_m
      logical flat_longi
      integer jc(6)
      type(probe_8) yy
      real(dp) x(6)
!       integer :: dodo = 0
!       if(dodo==1) then
!          x=zero
!          call FIND_ORBIT_x(x,default,c_1d_7,fibre1=my_ring%start)
!          write(6,*) x
!          call init_all(default,1,0)
!          call alloc(yy)
!          call alloc(id)
!          id=1
!          yy%x=x+id
!          call propagate(my_ring,yy,default)
!          call print(yy,6)
!          stop 999
!       endif

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

      if ( checksymmetric(re(:6,:6)) ) then
        call seterrorflag(11,"ptc_twiss ","The covariance matrix is not symmetric, bailing out.");
        return;
      endif



      nd=2
      if(icase==6 .or. icase==56) nd=3
      nd_m=nd

      flat_longi=.false.
      if (getdistrtype(3) /= distr_gauss.and.nd==3) then
         !here we have flat in delta
         nd_m=2
         flat_longi=.true.
      endif

      call init_all(default,2,nd)

      call alloc(ht)
      call alloc(ht0)
      call alloc(h)
      call alloc(id)
      call alloc(norm)

      do i = 1,nd_m*2
         do ii = 1,nd_m*2

            if (abs(re(i,ii)) > 1e-24) then
               ht0 =  re(i,ii)* (-1)**(ii+i) *(1.0_dp.mono.jc(i))*(1.0_dp.mono.jc(ii))
               !print*,'(',i,ii,')'
               !call daprint(ht0,6)
               ht=ht + ht0
            endif
         enddo
      enddo
      ht=-ht*pi

      lam=ten*full_abs(ht)
      ht=ht/lam

      if(flat_longi) then !1959 is the yearh of birth of Etienne (just a number)
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
      print*,"Emittances: ", emi(1),emi(2),emi(3)

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


         if (k.eq.0) then
            jj(1)=nx
            jj(2)=nxp
            jj(3)=ny
            jj(4)=nyp
            jj(5)=ndeltap
            jj(6)=nt
            call pok(theTransferMap%x(index)%T,jj,coeff)
            ! the following gives the same result as the above
            !oldv = theTransferMap(index).sub.jj
            !newtoset = (coeff - oldv).mono.jj ! mono for monomial
            !theTransferMap(index)%t=theTransferMap(index)%t+newtoset

         endif

         row = row+1

      enddo

!      call daprint(A_script_probe,28) ! to be compared with fort.18 created by ptc_normal

      ignore_map_orbit = get_value('ptc_twiss ','ignore_map_orbit ') .ne. 0

      if ( .not. ignore_map_orbit ) then
        do row=1,6
          orbit(row) = theTransferMap%x(row).sub.'0'
        enddo

        orbit_probe = orbit
      endif

      call maptoascript()

    end subroutine readmatrixfromtable

    !_________________________________________________________________
    subroutine readinitialmap ! from fort.18 file
      implicit none
      type(damap) :: map
      integer :: row
      logical(lp) :: ignore_map_orbit
      ! call readMapFromFort18(y)
      call alloc(map)

      ! I don't know why, but between initialization, when these flags are set to false, and this point they are flipped to true
      read77  = .false.
      print77 = .false.
      !print*,"read77=", read77
      call dainput(map,18)
      close(18)

      theTransferMap%x = map
      call kill(map)


      ignore_map_orbit = get_value('ptc_twiss ','ignore_map_orbit ') .ne. 0

      if ( .not. ignore_map_orbit ) then
        do row=1,6
          orbit(row) = theTransferMap%x(row).sub.'0'
        enddo

        orbit_probe = orbit
      endif


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
      read77  = .false.
      print77 = .false.
      call dainput(map,19)
      close(19)
      A_script_probe%x = orbit + map
      call kill(map)
      call reademittance()

    end subroutine readinitialascript
    !_________________________________________________________________

    subroutine reademittance
      !initializes A_script_probe(6) from re(6,6)
      implicit none
      real(dp) :: emix,emiy,emiz

      emix = get_value('probe ','ex ')
      emiy = get_value('probe ','ey ')
      emiz = get_value('probe ','et ')

      call setemittances(emix,emiy,emiz)

    end subroutine reademittance
    !_________________________________________________________________

    subroutine initmapfrommatrix
      !initializes A_script_probe(6) from re(6,6)
      implicit none
      real(dp),dimension(ndim2)::reval,aieval
      real(dp),dimension(ndim2,ndim2)::revec,aievec
      real(dp):: checkvalue
      real(dp) :: orbit(6)
      type(taylor) :: t

      call alloc(t)

      allocate(j(6))
      j(:)=0

      do i=1,c_%npara
        orbit(i) = theTransferMap%x(i).sub.j
      enddo

  !    call kill(theTransferMap)
  !    call alloc(theTransferMap)
  !
     ! theTransferMap%u = .false.
     ! theTransferMap%x = c_%npara
     ! theTransferMap%x = orbit


      call liepeek(iia,icoast)
      if (getdebug() > 1) then
         write (6,'(8a8)')   "no","nv","nd","nd2","ndc","ndc2","ndt","ndpt"
         write (6,'(8i8)') iia(1),iia(2),iia(3),iia(4),icoast(1),icoast(2),icoast(3),icoast(4)
         print*, "c_%npara is ", c_%npara
      endif

      do i = 1,6

        t = orbit(i)

        do ii = 1,c_%npara
            j(ii)=1

            r=re(i,ii)
            t = t+(r.mono.j)

            if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
              write(*,*) 'DA got unstable during Setup of transfer map i=',i,' ii=',ii
            endif

            j(ii)=0
         enddo

         theTransferMap%x(i) = t
      enddo


      if (getdebug() > 2) then
        print*,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
        print*,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
        print*,"PRINT "
        call print(theTransferMap,6)
        print*,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
        print*,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
        print*,"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
      endif
      deallocate(j)

!      call daprint(A_script_probe,29) ! to be compared with fort.18 created by ptc_normal and fort.28


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
      type(c_normal_form) theNormalForm
      type(c_damap)  :: c_Map, a_cs
      integer :: mf

      if (getdebug() > 2) then
         print*,"maptoascript: doing normal form"
      endif

      call alloc(c_Map)
      c_Map = theTransferMap

      call alloc(theNormalForm)
      call  c_normal(c_Map,theNormalForm)       ! (4)

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
         call print(theNormalForm%a_t,19)

      endif

      !print*, "maptoascript: TUNES    NF ", theNormalForm%tune

      call kill(A_script_probe)
      call alloc(A_script_probe)

      A_script_probe%u=my_false

      !use Courant Snyder
      call alloc(a_cs)
      call c_full_canonise(theNormalForm%atot,a_cs)   ! (0)
      A_script_probe = orbit_probe +  a_cs

     ! A_script_probe =  orbit_probe + theNormalForm%a_t

      if (getdebug() > 2) then

        call kanalnummer(mf,file="NormalFormA_t.txt")
        call print(theNormalForm%a_t,mf)
        close(mf)

        call kanalnummer(mf,file="NormalFormA1.txt")
        call print(theNormalForm%a1,mf)
        close(mf)

        call kanalnummer(mf,file="Ascript_start.txt")
        call print(A_script_probe,mf)
        close(mf)

      endif

      call kill(theNormalForm)
      call kill(c_Map)
      call kill(a_cs)

    end subroutine maptoascript
    !_________________________________________________________________

    subroutine readinitialtwiss(dt)
      !Reads initial twiss parameters from MAD-X command
      implicit none
      double precision alpha(3),beta(3),disp(4),mu(3)
      type(real_8) al(3),be(3),di(4)
      type(pol_block_inicond) :: inicondknobs
      integer k_system
      real(dp)  sizept
      double precision emiz
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


      orbit(:)=zero

      call readreforbit()

      if (getdebug() > 0 ) then
         CALL write_closed_orbit(icase,orbit)
      endif

      call reademittance()

      !Here we initialize A_script_probe(6)

      call alloc(be); call alloc(al); call alloc(di)

      !  code to power knows
      !


      do i=1,c_%nd
         !be(i)=beta(i)
         !al(i)=alpha(i)
         be(i)= beta(i)/(1_dp+orbit(5))
         al(i)=alpha(i)/(1_dp+orbit(5))
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


      A_script_probe%x=orbit

      do i=1,c_%nd
                  !print*, " Beta(", i,")=", beta(i)
                  !call print(A_script_probe(2*i-1),6)
                  !call print(A_script_probe(2*i  ),6)

         A_script_probe%x(2*i-1)= orbit(2*i-1) + sqrt(be(i)) * morph((one.mono.(2*i-1))    )
         A_script_probe%x(2*i)  = orbit(2*i)   + (one/sqrt(be(i)) * &
              (morph(  (one.mono.(2*i)) )-(al(i)) * morph((one.mono.(2*i-1)))))

                  !call print(A_script_probe(2*i-1),6)
                  !call print(A_script_probe(2*i  ),6)
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

            A_script_probe%x(5) = orbit(5) +  sqrt( beta(3) )*morph((one.mono.5))
            A_script_probe%x(6)=  orbit(6) + one/sqrt(beta(3)) * (morph(  (one.mono.6) )-(alpha(3)) * morph(one.mono.5))

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
            !skowron: it is x
            A_script_probe%x(5) = orbit(5) + morph((one.mono.5))
            A_script_probe%x(6) = orbit(6) + morph((one.mono.5))
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
            A_script_probe%x(i)= A_script_probe%x(i) + di(i) * morph((one.mono.5))
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
         write(6,'(6f8.4)') orbit
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    subroutine oneTurnSummary(oneTurnMap,theAscript,startorbit,suml)

      implicit none
      type(real_8),target :: theAscript(6)
      type(probe_8),target :: oneTurnMap
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
      real(dp) :: b,t1,t2,t3,t4 ! stuff needed for calculating high order compaction factor out of time slip factor
      real(dp) :: dispersion(4)
      real(dp) :: dispersion_p(4) ! derivative of the dispersion w.r.t delta-p/p
      integer :: debugFiles, mf
      integer :: icase
      integer :: order
      real(dp) :: rdp_mmilion = -1e6
      type(c_damap) :: yy, theA_beta ! added on November 6th to retreive momemtum compaction without differentiating the formula
      type(c_vector_field) :: vf_kernel
      type(c_damap)        :: c_Map
      type(c_damap)        :: rotationMap, theA_CS, theA_Spin, theA0, theA1, theA2, theRot
      type(c_taylor)       :: thePhaseAndSlip(3), theNuSpin
      type(c_normal_form)  :: theNormalForm


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
               write(21,*) "r(",i,j,")=",oneTurnMap%x(i).sub.fo(j,:)
            enddo
         enddo
         flush(21)
      endif



      if (getdebug() > 1) then
        write(6,*) "Doing normal form ... "
      endif

      ! 4. retreive the dispersion coefficients
      ! (may be the coefficient of delta of the map?)
      ! decompose the map via a normal form to get the dispersion...

      !print*,"Skowron 3 ", check_stable, c_%stable_da

      call alloc(c_Map)

      !print*,"Skowron 4 ", check_stable, c_%stable_da

      c_Map = oneTurnMap

      !print*,"Skowron 5 ", check_stable, c_%stable_da

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

         call kill(c_Map)

         return
      endif

      call alloc(theNormalForm)
      call  c_normal(c_Map,theNormalForm)       ! (4)

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

         call kill(c_Map)
         call kill(theNormalForm)

         return
      endif

      if (getdebug() > 1) then
        write(6,*) "Doing normal form ... Done"
      endif


      if (debugFiles .eq. 1) then

        call kanalnummer(mf,file="NormalFormA_t.summ.txt")
        call print(theNormalForm%a_t,mf)
        close(mf)

        call kanalnummer(mf,file="NormalFormA1.summ.txt")
        call print(theNormalForm%a1,mf)
        close(mf)

      endif

      !********************************************!
      ! new all-in-ine algorithm to get normal form goodies
      ! unstable for 6D
!       call alloc(rotationMap)
!
!       rotationMap=theNormalForm%atot**(-1)*c_Map*theNormalForm%atot ! (7a)
!
!       call alloc(theA_CS)
!       call alloc(theA_Spin)
!       call alloc(theA0)
!       call alloc(theA1)
!       call alloc(theA2)
!       call alloc(theRot)
!
!       call alloc(thePhaseAndSlip)
!       call alloc(theNuSpin)
!
!       call c_full_canonise(rotationMap,theA_CS,theA_Spin,theA0,theA1,theA2,theRot,phase=thePhaseAndSlip,nu_spin=theNuSpin)
!
!       call print(thePhaseAndSlip(1),6)
!       call print(thePhaseAndSlip(2),6)
!       call print(thePhaseAndSlip(3),6)

      !********************************************!
      !skowron: old ptc
      !print*, "Cplx dispersion A1(1,5)", theNormalForm%A1%v(1).sub.'000010'


      if( (c_%npara==5)       .or.  (c_%ndpt/=0)  ) then
         !when there is no cavity it gives us dispersions
         dispersion(1) = real(theNormalForm%A_t%v(1).sub.'000010')
         dispersion(2) = real(theNormalForm%A_t%v(2).sub.'000010')
         dispersion(3) = real(theNormalForm%A_t%v(3).sub.'000010')
         dispersion(4) = real(theNormalForm%A_t%v(4).sub.'000010')
      elseif (c_%nd2 == 6) then
        !faster to use the A_script then normal form
        ! because the normal form A_t must be canonised to C-S and moved from cmplx to real_8
        call dispersion6D(theAscript,dispersion)

      else
         do i=1,4
            dispersion(i)=zero
         enddo
      endif


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

      gamma_tr = zero

      if ( icase.gt.4 ) then
        if ( default%time ) then
         ! ALL cases with time=true
         ! Here R56 is dT/ddelta
         ! 5. apply formulas from twiss.F:
         !sd = r(5,6)+r(5,1)*disp(1)+...+r(5,4)*disp(4)
         !print*,"ALPHA_C, GAMMA TR : TIME ON"

         sd = -1.0*(oneTurnMap%x(6).sub.fo(5,:)) ! 5/6 swap MADX/PTC
         !print*,'sd(0)',sd
         do i=1,4
            !print*, 'Disp',i,'=', dispersion(i), ' oneTurnMap%x(6).sub.fo(i,:) = ', oneTurnMap%x(6).sub.fo(i,:)
            sd = sd - (oneTurnMap%x(6).sub.fo(i,:))*dispersion(i)
            !print*,'sd(',i,')',sd
         enddo
         !print*,'sd(f)',sd

         eta_c = -sd * betaRelativistic**2 / suml

         alpha_c = one / gammaRelativistic**2 + eta_c

         if (icase.eq.56) then

           ! non lin way

           yy%n = 6
           call alloc(yy)

           do i=1,6
              yy%v(i) = oneTurnMap%x(i)%t
           enddo

           yy = theNormalForm%A_t**(-1) * yy * theNormalForm%A_t ! takes away all dispersion dependency


           t1 = yy%v(6).sub.'000010'


           b = betaRelativistic

          !   print*, 'First order alpha_c ',alpha_c

             alpha_c = (suml - b**2*suml + b**2*t1)/suml

           !  print*, 'New alpha_c ',alpha_c

           ! algorithm from F.Schmidt


           if (order.ge.2) then
             t2 = yy%v(6).sub.'000020'
             alpha_c_p = b**2*(3*(-1 + b**2)*suml - 3*(-1 + b**2)*t1 + 2*b*t2)
             alpha_c_p = alpha_c_p / suml
           endif

           if (order.ge.3) then
             t3 = yy%v(6).sub.'000030'
             alpha_c_p2 = 3*b**2*((-1 + 6*b**2 - 5*b**4)*suml + t1 - 6*b**2*t1 + 5*b**4*t1 + &
                          4*b*t2 - 4*b**3*t2 + 2*b**2*t3)
             alpha_c_p2 = alpha_c_p2 / suml
           endif

           if (order.ge.4) then
             t4 = yy%v(6).sub.'000040'
             alpha_c_p3 = 3*b**3*(35*b**5*(suml - t1) + 10*t2 + 30*b**4*t2 - &
                          10*b**3*(5*suml - 5*t1 + 2*t3) + 5*b*(3*suml - 3*t1 + 4*t3) + &
                          8*b**2*(-5*t2 + t4))
             alpha_c_p3 = alpha_c_p3 / suml
           endif



            ! l5= (-15*b**3*(63*b**7*(suml - t1) - 2*t2 + 56*b**6*t2 - &
            !     21*b**5*(5*suml - 5*t1 + 2*t3) - 3*b*(suml - t1 + 6*t3) + &
            !     12*b**2*(3*t2 - 2*t4) + b**4*(-90*t2 + 24*t4) + &
            !     b**3*(45*suml - 45*t1 + 60*t3 - 8*t5)))/120d0

           call kill(yy)

         endif

        else

          if (icase.eq.6 .or. icase.eq.5)  then
            sd = -1.0*(oneTurnMap%x(6).sub.fo(5,:)) ! 5/6 swap MADX/PTC
            do i=1,4
               sd = sd - (oneTurnMap%x(6).sub.fo(i,:))*dispersion(i)
            enddo

           alpha_c = -sd/suml
           eta_c = alpha_c - one / gammaRelativistic**2



          else !5D or 6D
           !!!!!!!!!!!!!!!!!!!!1
           !56D

             yy%n = 6
             call alloc(yy)

             do i=1,6
                yy%v(i) = oneTurnMap%x(i)%t
             enddo

             yy = theNormalForm%A_t**(-1) * yy * theNormalForm%A_t ! takes away all dispersion dependency

             alpha_c    = (yy%v(6).sub.'000010')/suml
             eta_c = alpha_c - one / gammaRelativistic**2


             if (order.ge.2) then
                alpha_c_p  = 2.0*(yy%v(6).sub.'000020')/suml
             endif

             if (order.ge.3) then
                alpha_c_p2 = 6.0*(yy%v(6).sub.'000030')/suml
             endif

             if (order.ge.4) then
                alpha_c_p3 = 24.0*(yy%v(6).sub.'000040')/suml
             endif

          !   print*,'New 5D algo a_c ', alpha_c,' alpha_c_p ', alpha_c_p ,' e_c ', eta_c

             call kill(yy)

           endif ! 5D or 56D
        endif !time=false
      endif !icase=4

      if (alpha_c .gt. zero) then
        gamma_tr = one / sqrt(alpha_c)
      endif

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!! HIGHER ORDERS IN dP/P in time true  !!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! this algorithm does not work, for low energy (beta0 << 1) it
      ! obtains 2nd order time slip factor not alpha_c and the conversion is not straight forward

 !    if ( (icase.eq.5) .and. (default%time .eqv. .true.) ) then
      !
      !    !print*,"ALPHA_C dp/dp derivatives : 5D TIME ON"
      !
      !    ! compute delta-p/p dependency of alpha_c
      !    ! first order derivatives of the dispersions
      !    ! assuming icase=5
      !    !call daprint(oneTurnMap%x,88)
      !
      !    ! always assume time=true, so that the fifth phase-space variable is deltap instead of pt!
      !    ! otherwise should issue a warning or an error!
      !
      !    ! proceed slightly differently to retreive the dispersion's derivatives
      !    ! (so far only managed to get the dispersion on the normal form, but not its derivatives)
      !
      !    ! note: should reuse information computed previously than redo this...
      !    ! if icase=5, Taylor series expansion disp = disp + delta_p
      !    do i=1,4
      !       ! apparently, we are up by a factor two
      !       dispersion_p(i) = 2.0*(theAscript(i).sub.'000020')  ! as usual in this file
      !
      !    enddo
      !
      !    ! compute derivative of the formula for eta_c
      !
      !    sd = 0.0
      !    ! first dr65/ddeltap
      !    sd = sd - 1.0*(oneTurnMap%x(6).sub.'100010')*dispersion(1) &
      !         &  - 1.0*(oneTurnMap%x(6).sub.'010010')*dispersion(2) &
      !         &  - 1.0*(oneTurnMap%x(6).sub.'001010')*dispersion(3) &
      !         &  - 1.0*(oneTurnMap%x(6).sub.'000110')*dispersion(4) &
      !         &  - 2.0*(oneTurnMap%x(6).sub.'000020') &
      !	            ! new
      !         &  - 1.0*(oneTurnMap%x(6).sub.'000011')*(oneTurnMap%x(6).sub.'000001') &
      !	            ! then terms dr6i/ddeltap*dispersion(i)
      !         &  - 2.0*(oneTurnMap%x(6).sub.'200000')*dispersion(1)**2 &
      !         &  - (oneTurnMap%x(6).sub.'110000')*dispersion(1)*dispersion(2) &
      !         &  - (oneTurnMap%x(6).sub.'101000')*dispersion(1)*dispersion(3) &
      !         &  - (oneTurnMap%x(6).sub.'100100')*dispersion(1)*dispersion(4) &
      !	            ! new
      !         &  - (oneTurnMap%x(6).sub.'100010')*dispersion(1) &
      !	            ! dr62/ddeltap*disperion(2)
      !         &  - (oneTurnMap%x(6).sub.'110000')*dispersion(1)*dispersion(2) &
      !         &  - 2.0*(oneTurnMap%x(6).sub.'020000')*dispersion(2)**2 &
      !         &  - (oneTurnMap%x(6).sub.'011000')*dispersion(2)*dispersion(3) &
      !         &  - (oneTurnMap%x(6).sub.'010100')*dispersion(2)*dispersion(4) &
      !	            ! new
      !         &  - (oneTurnMap%x(6).sub.'010010')*dispersion(2) &
      !	            ! dr63/ddeltap*dispersion(3)
      !         &  - (oneTurnMap%x(6).sub.'101000')*dispersion(1)*dispersion(3) &
      !         &  - (oneTurnMap%x(6).sub.'011000')*dispersion(2)*dispersion(3) &
      !         &  - 2.0*(oneTurnMap%x(6).sub.'002000')*dispersion(3)**2 &
      !         &  - (oneTurnMap%x(6).sub.'001100')*dispersion(3)*dispersion(4) &
      !	            ! new
      !         &  - 1*(oneTurnMap%x(6).sub.'001010')*dispersion(3) &
      !	            ! dr64/ddeltap*dispersion(4)
      !         &  - (oneTurnMap%x(6).sub.'100100')*dispersion(1)*dispersion(4) &
      !         &  - (oneTurnMap%x(6).sub.'010100')*dispersion(2)*dispersion(4) &
      !         &  - (oneTurnMap%x(6).sub.'001100')*dispersion(3)*dispersion(4) &
      !         &  - 2.0*(oneTurnMap%x(6).sub.'000200')*dispersion(4)**2 &
      !	            ! new
      !         &  - 1*(oneTurnMap%x(6).sub.'000110')*dispersion(4)
      !
      !    ! terms involving derivatives of the dispersions
      !    do i=1,4
      !        print*, 'sd=',sd,' disp=',dispersion_p(i)
      !       sd = sd - (oneTurnMap%x(6).sub.fo(i,:))*dispersion_p(i)
      !    enddo
      !
      !     print*, 'sd=',sd
      !    alpha_c_p = -sd * (betaRelativistic**2) / suml

      !    eventually, one could differentiate the above formula to obtain alpha_c_p2
      !     but for the time-being, expect icase to be 56 to compute alpha_c_p2 and alpha_c_p3.

  !   elseif ( (icase.eq.56) .and. (default%time .eqv. .true.) ) then ! here one may obtain the pathlength derivatives from the map
       ! here we get directly eta_c with all higher order terms
       ! converting them to alpha_c is not obvious. Run time false to get higher order alpha_c

  !   endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!                                !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!    END OF COMPACTION FACTOR    !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!                                !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call alloc(vf_kernel)
      vf_kernel=0
      call flatten_c_factored_lie(theNormalForm%ker,vf_kernel)

      ! also output the tunes ...
      fractionalTunes(1) = tuneFromComplexVF(vf_kernel%v(1).sub.'1000') ! as in So_fitting.f90
      fractionalTunes(2) = tuneFromComplexVF(vf_kernel%v(3).sub.'0010') ! as in So_fitting.f90


      ! 26 november 2009
      if (icase.eq.6) then
         fractionalTunes(3) = tuneFromComplexVF(vf_kernel%v(6).sub.'000001') ! to enter here icase must be 6 but not only:
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

      if (fractionalTunes(1) < 0) fractionalTunes(1) = fractionalTunes(1) + one
      if (fractionalTunes(2) < 0) fractionalTunes(2) = fractionalTunes(2) + one
      !if (fractionalTunes(3) < 0) fractionalTunes(3) = fractionalTunes(3) + one

      ! we can also take the tune directly from the normal form
      !print*, "TUNES    NF ", theNormalForm%tune
      !print*, "  vf_kernel ", fractionalTunes

      ! ... as well as the chromaticities
      if (icase.eq.5 .or. icase.eq.56) then
         chromaticities(1) = tuneFromComplexVF(vf_kernel%v(1).sub.'10001') ! as in So_fitting.f90
         chromaticities(2) = tuneFromComplexVF(vf_kernel%v(3).sub.'00101') ! as in So_fitting.f90

      ! LD: 2019.12.18 removed after complains of FS
!         if (default%time) then
!            chromaticities(:) = chromaticities(:) * betaRelativistic
!         endif
         ! to get chromaticities, went to higher order with above "call init_default(default,2,0)"
      else
         ! if icase = 6, delta_p is a phase-space variable and not an external parameter hence we can't compute chromaticies
         chromaticities(1) = zero
         chromaticities(2) = zero
      endif

      ! for debug: check the values by printing the map
      if (debugFiles .eq. 1) then

        call kanalnummer(mf,file="oneTurnMap.summ.txt")
        call print(oneTurnMap,mf)
        close(mf)

        call kanalnummer(mf,file="vf_kernel.summ.txt")
        call print(vf_kernel,mf)
        close(mf)

      endif

      call kill(vf_kernel)

      call kill(theNormalForm)
      call kill(c_Map)


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
      call double_to_table_curr( summary_table_name,'orbit_t ',-startorbit(6))

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

       call double_to_table_curr(summary_table_name,'xcomax ' ,maxOrbit(1))
       call double_to_table_curr(summary_table_name,'pxcomax ',maxOrbit(2))
       call double_to_table_curr(summary_table_name,'ycomax ' ,maxOrbit(3))
       call double_to_table_curr(summary_table_name,'pycomax ',maxOrbit(4))
       call double_to_table_curr(summary_table_name,'ptcomax ',maxOrbit(5))


       if( (default%totalpath .eq. 0) ) then
         call double_to_table_curr(summary_table_name,'tcomin ' ,-maxOrbit(6)) !change of sign to MADX
         call double_to_table_curr(summary_table_name,'tcomax ' ,-minOrbit(6)) !change of sign to MADX
       else
         !this breakes MADX notation that particles arriving earlier have negative T
         !but otherwise the total time/pathlengh gets megative
         ! to do it correctly we should use nominal particle time + diff
         call double_to_table_curr(summary_table_name,'tcomin ' ,minOrbit(6)) !change of sign to MADX
         call double_to_table_curr(summary_table_name,'tcomax ' ,maxOrbit(6)) !change of sign to MADX
       endif




  end subroutine putMinMaxRmses

  subroutine computeDeltapDependency(A_script,s1)
    implicit none
    type(real_8), intent(in)  :: A_script(6)
    type(twiss),  intent(inout)  :: s1
    integer :: k,i
    integer :: J(lnv) ! the map's coefficient selector, as usual
    integer :: Jderiv(lnv) ! to store the map's coefficient selector of the derivative w.r.t deltap
    double precision :: get_value ! C-function
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
          s1%beta_p(k,i)= (A_script(2*k-1).sub.Jderiv) * (A_script(2*k-1).sub.J) +  &
                          (A_script(2*k-1).sub.J)      * (A_script(2*k-1).sub.Jderiv)
          s1%alfa_p(k,i)= -(A_script(2*k-1).sub.Jderiv)* (A_script(2*k).sub.J) -    &
                           (A_script(2*k-1).sub.J)     * (A_script(2*k).sub.Jderiv)
          s1%gama_p(k,i)= (A_script(2*k).sub.Jderiv)*(A_script(2*k).sub.J) + &
                          (A_script(2*k).sub.J)     *(A_script(2*k).sub.Jderiv)
          J(2*i-1)=0
          J(2*i)=1
          Jderiv = J ! vector copy
          Jderiv(5)=1 ! the delta_p coefficient
          s1%beta_p(k,i) = s1%beta_p(k,i) &
               + (A_script(2*k-1).sub.Jderiv)*(A_script(2*k-1).sub.J) + (A_script(2*k-1).sub.J)*(A_script(2*k-1).sub.Jderiv)
          s1%alfa_p(k,i)= s1%alfa_p(k,i) &
               - (A_script(2*k-1).sub.Jderiv)*(A_script(2*k).sub.J) - (A_script(2*k-1).sub.J)*(A_script(2*k).sub.Jderiv)
          s1%gama_p(k,i)= s1%gama_p(k,i) &
               + (A_script(2*k).sub.Jderiv)*(A_script(2*k).sub.J)+(A_script(2*k).sub.J)*(A_script(2*k).sub.Jderiv)
          J(2*i)=0
       enddo
    enddo

    ! the computations above match the following formulas, obtained by derivation of the Twiss parameters using the chain-rule
    ! beta derivatives w.r.t delta_p
    !    beta11 = two * (A_script(1).sub.'100000')*(A_script(1).sub.'100010') + two * (A_script(1).sub.'010000')*(A_script(1).sub.'010010')
    !    beta12 = two * (A_script(1).sub.'001000')*(A_script(1).sub.'001010') + two * (A_script(1).sub.'000100')*(A_script(1).sub.'000110')
    !    beta21 = two * (A_script(3).sub.'100000')*(A_script(3).sub.'100010') + two * (A_script(3).sub.'010000')*(A_script(3).sub.'010010')
    !    beta22 = two * (A_script(3).sub.'001000')*(A_script(3).sub.'001010') + two * (A_script(3).sub.'000100')*(A_script(3).sub.'000110')
    ! alpha derivatives w.r.t delta_p
    !    alfa11 = -((A_script(1).sub.'100010')*(A_script(2).sub.'100000')+(A_script(1).sub.'100000')*(A_script(2).sub.'100010')+&
    !         (A_script(1).sub.'010010')*(A_script(2).sub.'010000')+(A_script(1).sub.'010000')*(A_script(2).sub.'010010'))
    !    alfa12 = -((A_script(1).sub.'001010')*(A_script(2).sub.'001000')+(A_script(1).sub.'001000')*(A_script(2).sub.'001010')+&
    !         (A_script(1).sub.'000110')*(A_script(2).sub.'000100')+(A_script(1).sub.'000100')*(A_script(2).sub.'000110'))
    !    alfa21 = -((A_script(3).sub.'100010')*(A_script(4).sub.'100000')+(A_script(3).sub.'100000')*(A_script(4).sub.'100010')+&
    !         (A_script(3).sub.'010010')*(A_script(4).sub.'010000')+(A_script(3).sub.'010000')*(A_script(4).sub.'010010'))
    !    alfa22 = -((A_script(3).sub.'001010')*(A_script(4).sub.'001000')+(A_script(3).sub.'001000')*(A_script(4).sub.'001010')+&
    !         (A_script(3).sub.'000110')*(A_script(4).sub.'000100')+(A_script(3).sub.'000100')*(A_script(4).sub.'000110'))
    ! gamma derivatives w.r.t delta_p
    !    gama11 = (A_script(2).sub.'100010')*(A_script(2).sub.'100000')+(A_script(2).sub.'100000')*(A_script(2).sub.'100010')+&
    !         (A_script(2).sub.'010010')*(A_script(2).sub.'010000')+(A_script(2).sub.'010000')*(A_script(2).sub.'010010')
    !    gama12 = (A_script(2).sub.'001010')*(A_script(2).sub.'001000')+(A_script(2).sub.'001000')*(A_script(2).sub.'001010')+&
    !         (A_script(2).sub.'000110')*(A_script(2).sub.'000100')+(A_script(2).sub.'000100')*(A_script(2).sub.'000110')
    !    gama21 = (A_script(4).sub.'100010')*(A_script(4).sub.'100000')+(A_script(4).sub.'100000')*(A_script(4).sub.'100010')+&
    !         (A_script(4).sub.'010010')*(A_script(4).sub.'010000')+(A_script(4).sub.'010000')*(A_script(4).sub.'010010')
    !    gama22 = (A_script(4).sub.'001010')*(A_script(4).sub.'001000')+(A_script(4).sub.'001000')*(A_script(4).sub.'001010')+&
    !         (A_script(4).sub.'000110')*(A_script(4).sub.'000100')+(A_script(4).sub.'000100')*(A_script(4).sub.'000110')

    ! now compute deltap dependencies of the dispersion

    ! code to be differentiated, as in subroutine 'equaltwiss'
    !    J=0
    !    !here ND2=4 and delta is present      nd2=6 and delta is a constant
    !    !      print*,"nv",c_%nv,"nd2",c_%nd2,"np",c_%np,"ndpt",c_%ndpt ,"=>",c_%nv-c_%nd2-c_%np
    !    if( (c_%npara==5)       .or.  (c_%ndpt/=0) ) then
    !       !when there is no cavity it gives us dispersions
    !       do i=1,4
    !          lat(0,i,1)=(A_script(i).sub.J5)
    !       enddo
    !    elseif (c_%nd2 == 6) then
    !       do i=1,4
    !          lat(0,i,1) =              (A_script(i).sub.J5)*(A_script(6).sub.J6)
    !          lat(0,i,1) = lat(0,i,1) + (A_script(i).sub.J6)*(A_script(5).sub.J5)
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
             s1%disp_p(i) = 2.0*(A_script(i).sub.'000020')
          else
             s1%disp_p(i) = rdp_mmilion
          endif
          if (no.gt.2) then
             s1%disp_p2(i) = 3.0*2.0*(A_script(i).sub.'000030') ! assume at least order 3 for the map
          else
             s1%disp_p2(i) = rdp_mmilion
          endif
             !call fort_warn("ptc_twiss ","assume no>=3 for dispersion's 2nd order derivatives w.r.t delta-p")
          !endif
          if (no.gt.3) then
             s1%disp_p3(i) = 4.0*3.0*2.0*(A_script(i).sub.'000040')
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
               2.0*(A_script(i).sub.'000020')*(A_script(6).sub.'000001')+&
               (A_script(i).sub.'000010')*(A_script(6).sub.'000011')+&
               (A_script(i).sub.'000011')*(A_script(5).sub.'000010')+&
               (A_script(i).sub.'000001')*2.0*(A_script(5).sub.'000020')
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


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  subroutine normalFormAnalysis(oneTurnMap,theAscript,startorbit)
    use resindexfi
    implicit none
    type(probe_8),target :: oneTurnMap,theAscript
    real(dp),    target :: startorbit(6)
    real(dp)  :: prec ! for printing in files
    real(dp) :: disp1stOrder(4) ! for 6D algo
    type(c_taylor) :: tempTaylor ! for 6D algo
    type(c_normal_form) theNormalForm
    integer     :: mf, filecode ! output file
    integer     :: i,o ! output file
    character(len=250)   :: fmt
    integer     	:: io,r, myn1,myn2,indexa(mnres,4),mynres,ind(10)
    type(c_damap)  :: c_Map, c_Map2, q_Map, a_cs, a_cs_1
    type(c_taylor)  :: nrmlzdPseudoHam, g_io
    type(c_vector_field) vf, vf_kernel

     !use_complex_in_ptc=my_true

     call alloc(c_Map)

     c_Map = oneTurnMap;

     call alloc(theNormalForm)

     call  c_normal(c_Map,theNormalForm)       ! (4)

     !theNormalForm=oneTurnMap !! HERE WE DO NORMAL FORM TYPE 1

     if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
        write(whymsg,*) 'DA got unstable in Normal Form: PTC msg: ',&
                         messagelost(:LEN_TRIM(messagelost))
        call fort_warn('ptc_normal: ',whymsg)
        call seterrorflag(10,"ptc_normal ",whymsg);
        return
     endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !!
    if (getdebug()>2) then
       call kanalnummer(mf)
       open(unit=mf,file='normal.tfs')

       write(mf,'(a2,a16,1x,a4,1x,a10,1x)') '@ ',ch16lft('NAME'),'%09s', 'PTC_NORMAL'
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NO'),  '%9d', c_%no
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NV'),  '%9d', c_%nv
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('ND'),  '%9d', c_%nd
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('ND2'), '%9d', c_%nd2
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NDPT'),'%9d', c_%ndpt
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NPARA'),    '%9d', c_%npara
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NPARA_FPP'),'%9d', c_%npara_fpp
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NP_POL'),   '%9d', c_%np_pol
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('NSPIN'),    '%9d', c_%nspin
       write(mf,'(a2,a16,1x,a4,1x,i10,1x)') '@ ',ch16lft('SPIN_POS'), '%9d', c_%SPIN_pos

       !! make space for knobs, to be completed with proper handling of T and PT if absent
       i=1+c_%nv
       if (i .lt. 7) i=7
       write (fmt,'(a,i1,a)')  '(a2,2(a16,1x),ES16.8,',i,'(1x,i16))'
       !write (fmt,'(a)')  '(a2,2(a16,1x),ES16.8,'

       write(mf,'(a2,a16,9(1x,a16))') '* ',ch16lft('NAME'),ch16lft('NICKNAME'), &
                                      'VALUE',  'ORDER', &
                                      'ORDER_X','ORDER_PX', &
                                      'ORDER_Y','ORDER_PY', &
                                      'ORDER_PT','ORDER_T'
       write(mf,'(a2,a16,9(1x,a16))') '$ ',ch16lft('%s'),ch16lft('%s'), &
                                      '%le','%le', &
                                      '%le','%le', &
                                      '%le','%le', &
                                      '%le','%le'
    endif


   !!!!!!!!!!!!!!!!!!!!!!!!!!
   !! TUNES

   if (no > 1) then
     ! this algorithm cuts one order of magnitude so for order 1 does not work at all

      call alloc(vf_kernel)
      vf_kernel=0
      call flatten_c_factored_lie(theNormalForm%ker,vf_kernel)

      call putQnormaltable(vf_kernel%v(1),1)
      call putQnormaltable(vf_kernel%v(3),2)
      if (c_%nd2 == 6) then
        call putQnormaltable(vf_kernel%v(5),3)
      endif

      call kill(vf_kernel)

    else
      ! This puts only the linear part
      call putTunesNormalTable(theNormalForm%tune)

    endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!
   !! DISPERSION

    if (c_%nd2 == 6) then
      ! to be implemented
      call dispersion6D(theAscript%x,disp1stOrder)
      call alloc(tempTaylor)
      do i=1,4
        tempTaylor = disp1stOrder(i)
        call putDnormaltable(tempTaylor,i)
      enddo

     ! call dispersion6Dp(theAscript%x,disp1stOrder)

      call kill(tempTaylor)

    else
      !this does not work with 6D dispersion
      do i=1,4
        call putDnormaltable(theNormalForm%A_t%V(i),i)
      enddo
    endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!
   !! EQUILIBRIUM EMITTANCE
   !! only if radiation present

   if ( default%radiation .and. default%envelope) then

     !Emittances
     call putEmiNormalTable(theNormalForm%s_ijr)
     !Equilibrium Beam Sizes (Sigmas)
     call putESigNormalTable(theNormalForm%s_ij0)
     !Damping decrements
     call putDampingNormalTable(theNormalForm%damping)

   endif



   !!!!!!!!!!!!!!!!!!!!!!!!!!
   !! EIGEN VALUES


    call alloc(a_CS)

    ! raw transformation are subject to random rotations due to numerical instabilities
    ! c_canonise fixes them straight to fit the Courant Snyder format
    call c_canonise(theNormalForm%atot,a_CS)

    do i=1,c_%nd2 !from damap type def: Ndim2=6 but allocated to nd2=2,4,6

      !call putEnormaltable(theNormalForm%A_t%V(i),i)
      call putEnormaltable(a_CS%V(i),i)

    enddo


    !!!!!!!!!!!!!!!!!!!!!!
    !Generating functions

    !n%g is the vecotor field for the transformation
    !from resonance basis (action angle coordinate system, (x+ipx),(x-ipx)) back to cartesion X,Y
    !the ndim polynomials need to be flattened to get RDT's

        ! from c_normal
        !    n%a1=a1
        !    n%a2=a2
        !    n%a_t=a1*a2*from_phasor()*texp(n%g)*from_phasor(-1)


    !!!!!!!!!!!!!!!!!!!
    !
    ! the OLD algoritm
    ! call alloc(vf);
    ! call alloc(g_io);
    ! call flatten_c_factored_lie(theNormalForm%G,vf)
    ! g_io =-cgetpb(vf)
    !
    !!!!!!!!!!!!!!!!!!!

    call alloc(vf)
    call alloc(g_io)
    call alloc(a_CS_1)



    a_CS=to_phasor()*a_CS*from_phasor()
    call c_factor_map(a_CS,a_CS_1,vf,0)

    g_io = cgetpb(vf)

    call putGnormaltable(g_io)


    !!!!!!!!!!!!!!!!!!!!!!
    !HAMILTONIAN
    !!!!!!!!!!!!!! Normalised Pseudo-Hamiltonian !!!!!!!!!!!!!!!

    call c_canonise(theNormalForm%a_t,a_CS)

    a_CS = a_CS.sub.1
    call c_canonise(a_CS,a_CS)

    call alloc(c_Map2)
    c_Map2 = a_CS**(-1)*c_Map*a_CS

    call alloc(q_Map)
    call c_factor_map(c_Map2,q_Map,vf,dir=1)


    vf=from_phasor()*vf

    call alloc(nrmlzdPseudoHam)
    nrmlzdPseudoHam = cgetpb(vf)
    call putHnormaltable(nrmlzdPseudoHam)

    !!!!!!!!!!!!!!!!!!!!!!
    !ONE TURN MAP

    call putMnormaltable(oneTurnMap%x)

    !!!!!!!!!!!!!!!!!!!!!!
    !CLEANING

    call kill(vf)
    call kill(g_io)

    call kill(c_Map)
    call kill(c_Map2)
    call kill(q_Map)
    call kill(a_CS)
    call kill(a_CS_1)

    call kill(nrmlzdPseudoHam)

    call kill(theNormalForm)
   !if (icase.eq.5 .or. icase.eq.56) then
   !
   ! endif


  contains

    !left just of max 16 char string (+1 for C null termination)
    function ch16lft(in)
     implicit none
     character(*) :: in
     character(len=17) :: ch16lft

     write(ch16lft,'(a16)') in
     ch16lft = adjustl(ch16lft)
    end function ch16lft

    !terminates string with null for C
    subroutine ch16cterm(in)
     implicit none
     character(len=17) :: in
     integer zeropos
      zeropos = len_trim(in)+1
      if (zeropos > 17) zeropos = 17
      in(zeropos:zeropos) = char(0)

    end subroutine ch16cterm

    subroutine puttonormaltable(name, nick, basevar ,d_val,order,ind)
      implicit none
      character(len=17):: name, nick, basevar
      real(dp)    :: d_val
      integer     :: ind(10), order

      !print*,"Nick in:  ", nick
      call ch16cterm(name)
      call ch16cterm(nick)
      call ch16cterm(basevar)
      !print*,"Nick out: ", nick
      !print*,"name nick basevar Nick out: ", name, nick, basevar

      call string_to_table_curr(nl_table_name, 'name ', name)
      call string_to_table_curr(nl_table_name, 'nickname ', nick)
      call string_to_table_curr(nl_table_name, 'basevariable ', basevar)

      call double_to_table_curr(nl_table_name, 'value ', d_val)
      d_val = order
      call double_to_table_curr(nl_table_name, 'order ', d_val)
      d_val = ind(1)
      call double_to_table_curr(nl_table_name, 'order_x ', d_val)
      d_val = ind(2)
      call double_to_table_curr(nl_table_name, 'order_px ', d_val)
      d_val = ind(3)
      call double_to_table_curr(nl_table_name, 'order_y ', d_val)
      d_val = ind(4)
      call double_to_table_curr(nl_table_name, 'order_py ', d_val)
      d_val = ind(5)
      call double_to_table_curr(nl_table_name, 'order_pt ', d_val)
      d_val = ind(6)
      call double_to_table_curr(nl_table_name, 'order_t ', d_val)
      call augment_count(nl_table_name)

    end  subroutine puttonormaltable

    subroutine putMnormaltable(vv)
      implicit none
      type(real_8) :: vv(6)
      integer planei !1...6, 1=dx, 2=dpx
      character(len=1) :: planec
      character(len=2) :: planel
      character(len=1) :: d='M'
      character(len=17):: parname
      character(len=17):: nn, nick,basevar
      integer     :: ind(10), cnv, illa, n,nw, order, i,j
      real(dp)    :: d_val


      !print*,"Putting one turn map inside "
      ind(:) = 0

      do j=1,c_%nv

        nw = 0
        i=1

        if (vv(j)%alloc) then


          call dacycle(vv(j)%t%i,i,d_val,n)

          !print*, "MAP j=",j," has ",n," coefs"

          do i=1,N
             call dacycle(vv(j)%t%I,i,d_val,illa,ind(1:c_%nv))

             order = sum(ind(1:c_%nv))

             !print*, 'M',j,'_',ind(1),'_',ind(2),'_',ind(3), &
             !              '_',ind(4),'_',ind(5),'_',ind(6), d_val

             write(basevar,'(a1,i1)') 'm',j
             write(nn,'(a1,i1,6(a1,i1))') 'm',j,'_',ind(1),'_',ind(2),'_',ind(3), &
                                                '_',ind(4),'_',ind(5),'_',ind(6)

             if (getdebug() > 2) then
               write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                                  d_val, order, ind(1:6)
             endif

             call puttonormaltable(nn,nn,basevar,d_val,order,ind)

          enddo
        else
          if (vv(j)%kind == 0) then
             write(nn,'(a4,i1,a12)') 'm',j,'_0_0_0_0_0_0'

          endif
        endif
      enddo

    end subroutine putMnormaltable

    ! when switched to real table from file mf, move it out of "contains"
    subroutine putDnormaltable(v,planei)
      implicit none
      type(c_taylor) :: v
      integer planei !1...6, 1=dx, 2=dpx
      character(len=1) :: planec
      character(len=17) :: planel
      character(len=4) :: d='disp'
      character(len=17):: parname
      character(len=17):: nn, ni
      integer     :: ind(10), cnv, illa, n,nw, order, i
      real(dp)    :: d_val
      complex(dp) :: c_val

     ! print*, 'putDnormaltable, plane=',planei

      cnv = c_%nv
      if (cnv .lt.6)  cnv=6
      if (cnv .gt.10) cnv=10
      ind(:) = 0

      write(planec,'(i1)') planei
      select case(planei)
       case(1)
         planec='1'
         planel='x'
       case(2)
         planec='2'
         planel='px'
       case(3)
         planec='3'
         planel='y'
       case(4)
         planec='4'
         planel='py'
      end select

      parname = ch16lft(d//planec)

      nw = 0
      i=1

      call c_taylor_cycle(v,size=N)

      !print*, "There is ", N ," coefficients of dispersion plane ", planei
      do i=1,N

         call c_taylor_cycle(v,ii=i,value=c_val,j=ind(1:c_%nv))
         d_val = real(c_val)

         if (abs(d_val) < 1.e-15_dp) cycle ! skip this component, it is basically zero

         order = sum(ind(1:cnv))

         !print*, 'putDnormaltable, plane=',planei,' i=',i,' order=',order,' ind(5)=',ind(5)

         if (order /= ind(5)) then
           !it is the first order term that stays in the calculation with 1.0 coeff
           cycle
         endif

         !print*, 'putDnormaltable, plane=',planei,' i=',i,' Value=',d_val,  ind(1:c_%nv)

         nn = parname
         ni = 'd'//trim(planel)

         if (order > 1) then
           nn = trim(nn) // '_p'
           ni = trim(ni) // '_p'
         endif
         if (order > 2) then
           write(nn,'(a,i1)') trim(nn),order-1
           write(ni,'(a,i1)') trim(ni),order-1
         endif

         !print*,nn,ni
         !write(mf,fmt) '  ',ch16lft(nn),  ch16lft(ni), &
         !              d_val, sum(ind(1:cnv)), ind(1:cnv)

         call puttonormaltable(nn,ni,planel,d_val,order,ind)

         nw = nw + 1
      enddo

      if (nw == 0) then
         !print*,"there was nothing written, it means there is no dispersion. Put zero to the table"
         ni = 'd'//trim(planel)
         ind(:)=0
         ind(5)=1

         if (getdebug() > 2) then
           write(mf,fmt) '  ',ch16lft(parname),  ch16lft(ni), &
                          0.0, 1, ind(1:cnv)
         endif

         d_val = 0.0
         call puttonormaltable(parname,ni,planel,d_val,1,ind)
      endif

      !print*, 'putDnormaltable DONE'

    end subroutine putDnormaltable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    subroutine putHnormaltable(nrmlzdPseudoHam)
      implicit none
      type(c_taylor) :: nrmlzdPseudoHam
      type(taylor) :: sinham,cosham
      integer      :: ind(10)
      integer      :: order
      character(len=17):: nn, nick, bv='HAMILTONIAN'
      character(len=17):: hamiltsin='hamilt_sin'
      character(len=17):: hamiltcos='hamilt_cos'
      character(len=17):: hamiltamp='hamilt_amp'
      integer     	:: i,r, myn1,myn2,indexa(mnres,4),mynres
      complex(dp)   :: c_val
      real(dp)    :: im_val, re_val, d_val, eps=1e-6
      integer     :: maxorder,o
      double precision :: get_value ! C-function

      maxorder = get_value('ptc_twiss ', 'no ')

      ind(:) = 0
      myn1 = 0
      myn2 = 0
      mynres = 0
      i=1
      call c_taylor_cycle(nrmlzdPseudoHam,size=mynres)

    do o=1,maxorder !print order by order, I don't know how to sort c_taylor (piotr)

      do r=1,mynres

        call c_taylor_cycle(nrmlzdPseudoHam,ii=r,value=c_val,j=ind(1:c_%nv))

        order = sum(ind(1:6))

        if ( order .ne. o) then
          cycle
        endif

        re_val = real(c_val)
        im_val = imag(c_val)
        d_val  = hypot(re_val, im_val)

        ! if amplitude is close to zero then it is not worth to output
        if (d_val .lt. eps) then

          if (getdebug()>2) then
            print*,"putHnormaltable idx=",r," ",d_val," smaller then eps=",eps, " skipping "
          endif

          cycle

        endif


        !print*,"HAML order ",order, ind(:)
        !print*,'       im=',im_val , ' re=',re_val  , ' amp=', d_val


        !!!!!!!!!!!!!!!!!!!!!!!!
        write(nn,'(a4,6(a1,i1))') 'hama','_',ind(1),'_',ind(2),'_',ind(3), &
                                       '_',ind(4),'_',ind(5),'_',ind(6)
        write(nick,'(a2,6(i1))') 'h_',ind(1),ind(2),ind(3), &
                                      ind(4),ind(5),ind(6)

        if (getdebug() > 2) then
         write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nick), &
	    d_val, order, ind(1:6)
        endif

        call puttonormaltable(nn,nick,hamiltamp,d_val,order,ind)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(nn,'(a4,6(a1,i1))') 'hams','_',ind(1),'_',ind(2),'_',ind(3), &
                                       '_',ind(4),'_',ind(5),'_',ind(6)
        write(nick,'(a4,3(a1,SP,i2))') 'hams','_',ind(1)-ind(2),'_',ind(3)-ind(4),'_',ind(5)-ind(6)

        write(nick,'(a2,6(i1),a3)') 'h_',ind(1),ind(2),ind(3),ind(4),ind(5),ind(6),'_IM'

        if (getdebug() > 2) then
         write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nick), &
	    re_val, order, ind(1:6)
        endif

        call puttonormaltable(nn,nick,hamiltsin,im_val,order,ind)


        write(nn,'(a4,6(a1,i1))') 'hamc','_',ind(1),'_',ind(2),'_',ind(3), &
                                       '_',ind(4),'_',ind(5),'_',ind(6)
        write(nick,'(a2,6(i1),a3)') 'h_',ind(1),ind(2),ind(3),ind(4),ind(5),ind(6),'_RE'

        if (getdebug() > 2) then
         write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nick), &
	    im_val, order, ind(1:6)
        endif

        call puttonormaltable(nn,nick,hamiltcos,re_val,order,ind)

      enddo
     enddo
    end subroutine putHnormaltable


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    subroutine putGnormaltable(gen)
    !gets generating function that are the linear part of A_t
      implicit none
      type(c_taylor) :: gen
      integer     :: order
      integer     :: ind(10), i
      character(len=17):: nn, nick
      character(len=17):: genfunsin='gen_fun_sin'
      character(len=17):: genfuncos='gen_fun_cos'
      character(len=17):: genfunamp='gen_fun_amp'
      logical skew
      integer     	:: r, myn1,myn2,indexa(mnres,4),mynres, illa
      complex(dp)   :: c_val
      real(dp)    :: im_val, re_val, d_val,  eps=1e-6
      integer     :: maxorder
      double precision :: get_value ! C-function

      maxorder = get_value('ptc_twiss ', 'no ')

      ind(:) = 0
      myn1 = 0
      myn2 = 0
      mynres = 0
      i=1
      call c_taylor_cycle(gen,size=mynres)


      do o=1,maxorder !print order by order, I don't know how to sort c_taylor (piotr)

        do r=1,mynres

          call c_taylor_cycle(gen,ii=r,value=c_val,j=ind(1:c_%nv))

          order = sum(ind(1:6))

          if ( order .ne. o) then
            cycle
          endif

          !print*,"GNFU ",ind(1:6)

          im_val = imag(c_val)
          re_val = real(c_val)
          d_val  = hypot(re_val, im_val)

          ! if amplitude is close to zero then it is not worth to output
          if (d_val .lt. eps) then
            if (getdebug()>2) print*,"putGnormaltable idx=",r," ",d_val," smaller then eps=",eps, " skipping "
            cycle
          endif


          write(nn,'(a4,6(a1,i1))') 'gnfa','_',ind(1),'_',ind(2),'_',ind(3), &
                                          '_',ind(4),'_',ind(5),'_',ind(6)

          write(nick,'(a2,6(i1))') 'f_',ind(1),ind(2),ind(3), &
                    	ind(4),ind(5),ind(6)
          !
          !write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
          !               d_val, order, ind(1:6)

          !write (fmt,'(a,i1,a)')  '(a2,2(a16,1x),ES16.8,',7,'(1x,i16))'
          !write(6,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
          !               d_val, order, ind(1:6)

          call puttonormaltable(nn,nick,genfunamp,d_val,order,ind)



          write(nn,'(a4,6(a1,i1))') 'gnfs','_',ind(1),'_',ind(2),'_',ind(3), &
                                          '_',ind(4),'_',ind(5),'_',ind(6)
          write(nick,'(a2,6(i1),a3)') 'f_',ind(1),ind(2),ind(3), &
                                           ind(4),ind(5),ind(6),'_im'
          !write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
          !               im_val, order, ind(1:6)
          call puttonormaltable(nn,nick,genfunsin,im_val,order,ind)

          write(nn,'(a4,6(a1,i1))') 'gnfc','_',ind(1),'_',ind(2),'_',ind(3), &
                                          '_',ind(4),'_',ind(5),'_',ind(6)
          write(nick,'(a2,6(i1),a3)') 'f_',ind(1),ind(2),ind(3), &
                                           ind(4),ind(5),ind(6),'_re'
          !write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
          !               re_val, order, ind(1:6)
          call puttonormaltable(nn,nick,genfuncos,re_val,order,ind)

        enddo
      enddo

      myn1 = 0
      myn2 = 0
      mynres = 0


    end subroutine putGnormaltable


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    subroutine putEmiNormalTable(va)
    !gets Emittances
      implicit none
      complex(dp) :: va(6,6), v
      integer     :: ind(10), i ! ind has 10 elements for extension to knobs and clocks
      real(dp)    :: d_val, emi(3)
      character(len=17):: nn
      character(len=17):: bv = 'emi'


      do i=1,3

        ind(:) = 0
        ind(2*i-1) = 1
        ind(2*i  ) = 1

        d_val = abs(real(va(2*i-1,2*i)))/2.0_dp

        write(nn,'(a4,i1)') 'emit',i

        if (getdebug() > 2) then
           write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                         d_val, 2, ind(1:6)
        endif


        call puttonormaltable(nn,nn,bv,d_val,2,ind)

        ind(:) = 0

      enddo

    end subroutine putEmiNormalTable

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    subroutine putESignormaltable(va)
    !gets Equilibrium Beam sizes
      implicit none
      complex(dp) :: va(6,6), v
      integer     :: ind(10), i,j ! ind has 10 elements for extension to knobs and clocks
      real(dp)    :: d_val
      character(len=17):: nn
      character(len=17):: bv = 'EquilSigma'


      do j=1,c_%nd2

      ind(:) = 0
      ind(j) = 1

        do i=1,c_%nv

          ind(i) = ind(i) + 1

          v= va(j,i)

          if (abs(v) < 1.e-15_dp) cycle ! skip this component, it is basically zero

          d_val = real(v)

          write(nn,'(a10,2i1)') 'EquilSigma',j,i

          if (getdebug() > 2) then
             write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                           d_val, 2, ind(1:6)
          endif


          call puttonormaltable(nn,nn,bv,d_val,2,ind)

          ind(:) = 0

        enddo
      enddo

    end subroutine putESignormaltable

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    subroutine putEnormaltable(v,planei)
    !gets eigenvectors that are the linear part of A_t
      implicit none
      type(c_taylor) :: v
      integer     :: planei
      integer     :: ind(10), i ! ind has 10 elements for extension to knobs and clocks
      real(dp)    :: d_val
      character(len=17):: nn
      character(len=17):: bv = 'a_t'

      ind(:) = 0
      do i=1,c_%nv
         ind(i)=1
         d_val = real(v.sub.ind(1:6))

        ! if (planei == 3 .and. i == 3) then
        !   call print(v,6)
        !   print*, "skowron 33 real 1 ", d_val
        !   print*, "skowron 33 complex 1 ", v.sub.'001000'
        ! endif

         write(nn,'(a4,2i1)') 'eign',planei,i

         if (getdebug() > 2) then
            write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                          d_val, 1, ind(1:6)
         endif


         call puttonormaltable(nn,nn,bv,d_val,1,ind)

         ind(i)=0

      enddo

    end subroutine putEnormaltable


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine putDampingNormalTable(va)
    !gets Emittances
      implicit none
      real(dp) :: va(3)
      integer     :: ind(10), i ! ind has 10 elements for extension to knobs and clocks
      real(dp)    :: d_val, emi(3)
      character(len=17):: nn
      character(len=17):: bv = 'damping'

      ind(:) = 0

      do i=1,c_%nd

        d_val = va(i)

        write(nn,'(a7,i1)') 'damping',i

        if (getdebug() > 2) then
           write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                         d_val, 0, ind(1:6)
        endif

        call puttonormaltable(nn,nn,bv,d_val,0,ind)

      enddo

    end subroutine putDampingNormalTable

    subroutine putTunesNormalTable(va)
    !gets Emittances
      implicit none
      real(dp) :: va(3)
      integer     :: ind(10), i ! ind has 10 elements for extension to knobs and clocks
      real(dp)    :: d_val, emi(3)
      character(len=17):: nn
      character(len=17):: bv = 'q'

      ind(:) = 0

      do i=1,c_%nd

        d_val = va(i)

        write(nn,'(a1,i1)') 'q',i

        if (getdebug() > 2) then
           write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nn), &
                         d_val, 0, ind(1:6)
        endif

        call puttonormaltable(nn,nn,bv,d_val,0,ind)

      enddo

    end subroutine putTunesNormalTable


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    ! when switched to real table from file mf, move it out of "contains"
    subroutine putQnormaltable(v,planei)
      implicit none
      type(c_taylor) :: v
      integer planei
      character(len=1) :: planec,planel
      character(len=1) :: q='q', u='_'
      character(len=17):: parname
      character(len=17):: nn, nick, bv
      integer     :: ind(10)
      real(dp)    :: d_val
      complex(dp) :: c_val
      integer     :: i, ii, ioa, n, order, orderX, orderY, orderPT, orderT
      integer     :: cnv

      cnv = c_%nv
      if (cnv .lt.6)  cnv=6
      if (cnv .gt.10) cnv=10
      ind(:) = 0

      write(planec,'(i1)') planei
      select case(planei)
       case(1)
         planec='1'
         planel='x'
       case(2)
         planec='2'
         planel='y'
       case(3)
         planec='3'
         planel='s'
      end select

      parname = ch16lft(q//planec)
      bv = ch16lft(q//planec)
     !!!!!!!!!!!!!!!!!!!!!

      N=-1

      i=1
      call c_taylor_cycle(v,size=N)

      do i=1,N
        call c_taylor_cycle(v,ii=i,value=c_val,j=ind(1:c_%nv))
         !print*, 'Value=',d_val,  ind(1:c_%nv)
         d_val = -aimag(c_val)/(2.*pi)
         ind(2*planei - 1) = ind(2*planei - 1) - 1 !kernel has extra exponent at the plane variable, q1=v(1).sub.'100000'

         if (mod(sum(ind(1:2)),2) /=  0 ) then
           !it should be always even
           call fort_warn('ptc_twiss: ',' strange dependence of tune on horizontal coordinates')
         endif

         if (mod(sum(ind(3:4)),2) /=  0 ) then
           !it should be always even
           call fort_warn('ptc_twiss: ',' strange dependence of tune on vertical coordinates')
         endif

         !PTC gives the same order in x and px to mark Jx and Jy, i.e.
         ! dQx/dJx has index 210000 ; after substracting the extra exponent in the plane of variable it is 110000
         ! dQy/dJy has index 002100
         ! we prefer to fix px and py to zero so order = sum(ind(:)) is consistent
         ind(2) = 0
         ind(4) = 0
         if (c_%nd2 == 6) then
          ! 6D case, we get Jt in longitudinal plane ind(5) == ind(6)
          if (mod(sum(ind(5:6)),2) /=  0 ) then
            !it should be always even
            call fort_warn('ptc_twiss: ',' strange dependence of tune on longitudinal coordinates')
          endif
          ind(5) = 0
         endif

         order = sum(ind(1:cnv))

         nn = parname

         if (order == 0 ) then
           ! O R D R E R   Z E R O

           nick = parname  ! tune  q1
           !tune sometimes comes negative, add one in this case
           if (d_val .lt. zero) d_val = d_val + one

         else


           orderX = ind(1)
           orderY = ind(3)
           orderPT= ind(5)
           orderT = ind(6)


           if (orderX > 0) nn = trim(nn) // '_jx'
           if (orderX > 1) write(nn,'(a,i1)') trim(nn),orderX
           if (orderY > 0) nn = trim(nn) // '_jy'
           if (orderY > 1) write(nn,'(a,i1)') trim(nn),orderY

           if (c_%nd2 == 6) then
            ! 6D case, we get Jt in longitudinal plane
            if (orderT > 0) nn = trim(nn) // '_jt'
            if (orderT > 1) write(nn,'(a,i1)') trim(nn),orderT


           else
             if (orderPT> 0) nn = trim(nn) // '_p'
             if (orderPT > 1) write(nn,'(a,i1)') trim(nn),orderPT
             if (orderT> 0) nn = trim(nn) // '_t'
             if (orderT > 1) write(nn,'(a,i1)') trim(nn),orderT
           endif

           nick = nn

           if (order == ind(5) ) then
             ! chroma or higher order in momentum
             if (order == 1) then
               nick = 'd'//trim(parname) ! chroma dqN
             else
               nick = nn
             endif

           endif

           if ( order == (sum( ind(2*planei-1:2*planei))) ) then
             nick = 'anh'//planel
             if (order > 1) then
                write(nick,'(a,a1,i1)') trim(nick),'_',order !qN_JxM
             endif
           else
             if ( (order==1) .and. (sum(ind(5:6)) == 0) ) then
                nick = 'anhc'
             endif
           endif

        endif !else order==0

         if (getdebug() > 2) then
            write(mf,fmt) '  ',ch16lft(nn),  ch16lft(nick), &
                          d_val, order, ind(1:cnv)
         endif

         call puttonormaltable(nn,nick,bv,d_val,order,ind)

      ENDDO

     !!!!!!!!!!!!!!!!!!!!!!



    end subroutine putQnormaltable



  end subroutine normalFormAnalysis


  real(dp) function tuneFromComplexVF(c_val)
   implicit none
   complex(dp) :: c_val

   tuneFromComplexVF = -aimag(c_val)/(2.*pi)

!   call print(c_val,6)

  end function tuneFromComplexVF

!function buildVariableName(base, ind )
!  implicit none
!  character(*) :: base
!
!end function buildVariableName



end module madx_ptc_twiss_module
