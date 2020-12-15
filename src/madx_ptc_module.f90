MODULE ptc_results
  USE madx_keywords
  implicit none
  public
  integer :: number_variables = 6
  integer :: order = 20
  character(len = 2), dimension(6) :: ptc_variables = (/'x ','xp','y ','yp','z ','dp'/)
  character(len = 2) :: ptc_var
  !type(normalform) n            ! with the new complex da it goes local to the modules
  !type (pbresonance) pbrg,pbrh  ! with the new complex da it goes local to the modules
END MODULE ptc_results

MODULE madx_ptc_module
  use S_fitting_new
  !  USE madx_keywords
  USE madx_ptc_setcavs_module
  USE madx_ptc_knobs_module
  use madx_ptc_intstate_module, only : getdebug

  implicit none
  public

  logical(lp) mytime
  integer icav

  integer :: universe=0,index_mad=0,EXCEPTION=0
  integer ipause
  integer,external :: mypause
  real(kind(1d0)) get_value,node_value
  type(layout),pointer :: my_ring=>null(), bmadl=>null()
  type(mad_universe),pointer ::  m_u=>null(),m_t=>null()
  integer, private, parameter :: mynreso=20
  integer, private, dimension(4) :: iia,icoast
  real(dp) :: mux_default=c_0_28, muy_default=c_0_31, muz_default=c_1d_3
  integer, private, allocatable :: J(:)
  logical(lp)             :: savemaps=.false.
  logical(lp) :: resplit,even
  real(dp) my_thin,my_xbend

  type mapbuffer
     type(universal_taylor) :: unimap(6)
     real(dp)               :: s
     character(nlp+1)       :: name
  end type mapbuffer

  type(mapbuffer), pointer  :: maps(:) !buffered maps from the last twiss
  integer                   :: mapsorder = 0  !order of the buffered maps, if 0 maps no maps buffered
  integer                   :: mapsicase = 0


  type fibreptr
    type(fibre), pointer    :: p => null()
  end type fibreptr


  integer, private, parameter:: maxelperclock = 10 ! maximum 10 ac dipols with given clock
  type clockdef
     real(dp)                :: tune = -1 ! negative means inactive, in fact it is tune, left like this for backward compatibility, can be changed during LS2
     real(dp)                :: lag = 0
     integer                 :: rampupstart = 0,  rampupstop = 0, rampdownstart = 0,  rampdownstop = 0
     integer                 :: nelements = 0
     type(fibreptr)          :: elements(maxelperclock)
  end type clockdef
  
  
  integer, private, parameter:: nmaxclocks = 3
  type(clockdef),  dimension(nmaxclocks) :: clocks ! 3 pointers
  integer                                :: nclocks = 0

  real(dp) :: beta0start
  real(dp) :: my_ring_length
  
  character(1000), private  :: whymsg
  external :: aafail, dcopy, get_node_vector, fort_warn
  external :: element_name, node_name, node_string
  external :: augment_count, string_to_table_curr, vector_to_table_curr
  external :: make_map_table, seterrorflag, comm_para, double_to_table_curr
CONTAINS

  subroutine resetclocks()
    implicit none
    integer i

    do i=1,nmaxclocks
      clocks(i)%tune = -1 ! negative means inactive
      clocks(i)%lag = 0
      clocks(i)%rampupstart = 0
      clocks(i)%rampupstop = 0
      clocks(i)%rampdownstart = 0
      clocks(i)%rampdownstop = 0
    enddo
    nclocks = 0
  end subroutine resetclocks

  subroutine ptc_create_universe()
    implicit none
    real(kind(1d0)) get_value
    integer maxnmul

    use_quaternion=.true.

    piotr_freq=.true. ! PTC flag in cavity tracking to have correct phasing with time=false

    check_longitudinal = .true. ! PTC flag to check stability of the closed orbit in longitudinal
                                ! to prevent finding unstable fixed point

    call set_aperture_all_case0(.true.)

    print77=.false.
    read77 =.false.
    lingyun_yang=get_value('ptc_create_universe ','ntpsa ').ne.0
    lielib_print(6)=get_value('ptc_create_universe ','symprint ')

    nullify(maps)

    if (getdebug()==0) global_verbose = .false.

    lielib_print =  (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
   !  lielib_print(1)=1   lieinit prints info
   !  lielib_print(2)=1   expflo warning if no convergence
   !  lielib_print(3)=1   Shows details in flofacg
   !  lielib_print(4)=1   tunes and damping
   !  lielib_print(5)=1  order in orbital normal form
   !  lielib_print(6)=1  symplectic condition
   !  lielib_print(7)=-1  go manual in normal form  (use auto command in fpp)
   !  lielib_print(8)=-1  To use nplane from FPP normalform%plane
   !  lielib_print(9)=1  print in checksymp(s1,norm) in j_tpsalie.f90
   !  lielib_print(10)=1  print lingyun's checks
   !  lielib_print(11)=1  print warning about Teng-Edwards
   !  lielib_print(12)=1  print info in make_node_layout
   !  lielib_print(13)=1  print info of normal form kernel into file kernel.txt and kernel_spin.txt
                     !  lielib_print(14)=1  print info about recutting
                     !  lielib_print(15)=1  print info during flat file reading and printing


    if (getdebug()>0) then
        lielib_print(9)=1 !prints symplecticity deviation
    endif

    if (getdebug()>1) then
        print*,"Now PTC"
        lielib_print =  (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
    endif

    sector_nmul_max = get_value('ptc_create_universe ','sector_nmul_max ')

    !    print*,">>ss1<< old sector_nmul",sector_nmul



    !    print*,">>ss1<< new sector_nmul",sector_nmul
    if (sector_nmul_max < 0) then
       maxnmul = getmaxnmul()

       if (maxnmul < 0) then
          call aafail('ptc_create_universe: ',&
                      'Seems that no sequence is in currently in use. Aborting.')
          return
       endif

       sector_nmul_max = maxnmul
       sector_nmul     = maxnmul
    else
       sector_nmul = get_value('ptc_create_universe ','sector_nmul ')
    endif


    if(sector_nmul_max.lt.sector_nmul) then
       call aafail('sector_nmul_max must be larger than sector_nmul: ',&
            'check your ptc_create_universe input')
    endif

    ! copy from Ss_fake_mad.f90:ptc_ini_no_append



    allocate(m_u)
    call set_up_universe(m_u)
    allocate(m_t)
    call set_up_universe(m_t)

    universe=universe+1

    allocate(bmadl)
    call set_up(bmadl)
    bmadl%NAME='BMAD REUSED FIBRE LAYOUT'
    call point_m_u(m_u,m_t)



  end subroutine ptc_create_universe
  !_________________________________________________________________

  subroutine ptc_create_layout()
    implicit none
    real(kind(1d0)) get_value



    if(universe.le.0 .or. EXCEPTION.ne.0) then
       call fort_warn('return from ptc_create_layout: ',' no universe created')
       return
    endif

    call append_empty_layout(m_u)
    call append_empty_layout(m_t)


    index_mad=index_mad+1
    my_ring=>m_u%end

    call ptc_input()

    if(EXCEPTION.eq.1) then
       call fort_warn('wrong magnet type KINDI which must be: ','1, 2, 3')
       return
    endif

    cavsareset = .false.
    mytime=get_value('ptc_create_layout ','time ').ne.0

    if(mytime) then
       default=getintstate()
       default=default+time
       call setintstate(default)
    endif


  end subroutine ptc_create_layout
  !_________________________________________________________________

  subroutine ptc_move_to_layout()
    implicit none
    real(kind(1d0)) get_value
    integer my_index

    if(universe.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_move_to_layout: ',' no universe created')
       return
    endif

    my_index = get_value('ptc_move_to_layout ','index ')

    if(my_index.gt.index_mad.or.my_index.le.0) then
       call fort_warn('return from ptc_move_to_layout: ',' layout outside allowed range')
       print*,"   Allowed range 0 < ",index_mad
       return
    endif

    call move_to_layout_i(m_u,my_ring,my_index)

  end subroutine ptc_move_to_layout
  !_________________________________________________________________

  subroutine ptc_input()
    use twtrrfi
    use twiss0fi
    use name_lenfi
    implicit none
    logical(lp) particle,doneit,isclosedlayout
    integer i,j,k,code,nt,icount,nn,ns,nd,mg,napoffset,get_string
    !    integer get_option
    integer double_from_table_row,table_cell_exists
    integer restart_sequ,advance_node,n_ferr,node_fd_errors
    integer, external :: get_userdefined_geometry, get_userdefined_geometry_len
    integer, parameter :: nt0=20000
    real(dp) l,l_machine,energy,kin,brho,beta0,p0c,pma,e0f,lrad,charge
    real(dp) f_errors(0:maxferr),aperture(maxnaper),normal(0:maxmul), apoffset(2)
    real(dp) patch_ang(3),patch_trans(3)
    real(dp) skew(0:maxmul),field(2,0:maxmul),fieldk(2),myfield(2*maxmul+2)
    real(dp) gamma,gamma2,gammatr2,freq,offset_deltap
    real(dp) modulationq
    real(dp) fint,fintx,div,muonfactor,edge,rhoi,hgap,corr,tanedg,secedg,psip
    real(dp) sk0,sk1,sk1s,sk2,sk2s,sk3,sk3s,tilt,dum1,dum2
    REAL(dp) ::  normal_0123(0:3), skew_0123(0:3) ! <= knl(1), ksl(1)
    real(dp) gammatr,ks,ksi,ex,ey ! LD
    real(kind(1d0)) get_value,node_value
    character(name_len) name
    character(name_len) aptype
    type(keywords) key
    character(20)       keymod0,keymod1
    character(1024)     msg ! warning message buffer
    character(name_len) magnet_name
    logical(lp)         exact0,no_cavity_totalpath
    integer             exact1
    integer             sector_nmul_max0,sector_nmul0
    integer             model
    integer             method,method0,method1
    integer             nst0,nst1,ord_max,kk
    REAL (dp) :: tempdp,bvk
    logical(lp):: ptcrbend,truerbend,errors_out
    !  Etienne helical
    character(nlp) heli(100)
    integer mheli,helit,ihelit
    type(fibre), pointer :: p => null()
    double precision, parameter :: zero=0.d0
    integer, parameter :: aplen=0
    REAL(DP), pointer, dimension (:) :: apx => null()
    REAL(DP), pointer, dimension (:) :: apy => null()
    


    !real :: tstart, tfinish, tsum
    !tsum = 0d0
    !call cpu_time(tstart)


    !---------------------------------------------------------------
    !---------------------------------------------------------------
    if (getdebug() > 1) then
       print *, '--------------------------------------------------------------'
       print *, '--------------------------------------------------------------'
       print *, '------    E X E C U T I N G     P T C     I N P U T   --------'
       print *, '--------------------------------------------------------------'
       print *, '--------------------------------------------------------------'
    endif

    energy=get_value('probe ','energy ')
    pma=get_value('probe ','mass ')
    charge=get_value('probe ','charge ')
    bvk=get_value('probe ','bv ')

    e0f=sqrt(ENERGY**2-pma**2)

    if (getdebug() > 0) then
       print *, 'MAD-X Beam Parameters'
       print '(a26, e13.6)', '      Energy :',energy
       print '(a26, e13.6)', '      Kinetic Energy :',energy-pma
       print '(a26, e13.6)', '      Particle Rest Mass :',pma
       print '(a26, e13.6)', '      Momentum :',e0f
    endif



    beta0=e0f/ENERGY


    if(abs(pma-pmae)/pmae<c_0_002) then
       if (getdebug() > 1) then
           print *,'Executing MAKE_STATES(TRUE), i.e. ELECTRON beam'
       endif
       particle=.true.
       CALL MAKE_STATES(PARTICLE)
    elseif(abs(pma-pmap)/pmap<c_0_002) then
       if (getdebug() > 1) then
           print *,'Executing MAKE_STATES(FALSE), i.e. PROTON beam'
       endif
       particle=.false.
       CALL MAKE_STATES(PARTICLE)
    else
       muonfactor=pma/pmae
       if (getdebug() > 1) then
           print '(a, f8.4, a)','Executing MAKE_STATES(',pma/pmae,'), i.e. PROTON beam'
       endif
       CALL MAKE_STATES(muonfactor)
    endif

    !the state is cleared at this stage
    call setintstate(default)
    !valid October 2002: oldscheme=.false.
    !!valid October 2002: oldscheme=.true.


    if (getdebug()==0) global_verbose = .false.

    !  with_external_frame=.false.
    !  with_internal_frame=.false.
    !  with_chart=.false.
    !  with_patch=.false.

    ! Global Keywords

    if (getdebug() > 1) then
       print *, '=============================================================='
       print *, 'INPUT PARAMETERS ARE:'
    endif

    sector_nmul_max0 = sector_nmul_max
    if (getdebug() > 1) then
        print*,'  Global max sector_nmul: ',sector_nmul_max0
    endif

    sector_nmul0 = sector_nmul
    if (getdebug() > 1) then
        print*,'  Global sector_nmul: ',sector_nmul0
    endif


    model = get_value('ptc_create_layout ','model ')
    if (getdebug() > 1) then
        print*,'  Global Model code is : ',model
    endif

    !*****************************
    !  MODEL Settings
    !*****************************
    select case(model)
    CASE(1)
       keymod0 = "DRIFT_KICK       "
    CASE(2)
       keymod0 = "MATRIX_KICK      "
    CASE(3)
       keymod0 = "DELTA_MATRIX_KICK"
    CASE DEFAULT
       PRINT *, 'EXCEPTION occured: Can not recognize model type ',model
       EXCEPTION=1
       index_mad=-1
       ipause=mypause(444)
       RETURN
    END SELECT

    if (getdebug() > 1) then
        print*,'  Global Model name (keymod0) is : ',keymod0
    endif

    method   = get_value('ptc_create_layout ','method ')
    if (getdebug() > 1) then
        print*,'  Global method is: ',method
    endif

    !*****************************
    !  METHOD Settings
    !*****************************
    select case(method)
    CASE(2)
       method0 = method
    CASE(4)
       method0 = method
    CASE(6)
       method0 = method
    CASE DEFAULT
       PRINT *, 'EXCEPTION occured: Can not recognize method order ',method
       EXCEPTION=1
       index_mad=-1
       ipause=mypause(444)
       RETURN
    END SELECT

    exact0    = get_value('ptc_create_layout ','exact ') .ne. 0
    if (getdebug() > 1) then
        print*,'  Global exact is: ',exact0
    endif

    nst0      = get_value('ptc_create_layout ','nst ')
    if (getdebug() > 1) then
        print*,'  Global Number of Integration Steps (nst) is: ',nst0
    endif

    ! MAD-X specials
    !    madlength = get_option('rbarc ') .eq. 0
    madlength = .false.
    if (getdebug() > 1) then
        print*,'  global rbend_length: ',madlength
    endif

    mad       = get_value('ptc_create_layout ','mad_mult ') .ne. 0
    if (getdebug() > 1) then
        print*,'  global mad_mult as in mad8: ',mad
    endif

    mad8      = get_value('ptc_create_layout ','mad8 ') .ne. 0
    if (getdebug() > 1) then
        print*,'  rbend as in mad8 (only global): ',mad8
    endif

    gamma     = get_value('probe ','gamma ')
    if (getdebug() > 1) then
        print*,'  gamma: ',gamma
    endif

    if (table_cell_exists('summ ', 'gammatr ', 1).ne.0) then
      k       = double_from_table_row('summ ','gammatr ',1,gammatr)
    else
      gammatr = 0
    endif

    if (getdebug() > 1) then
        print*,'  gammatr: ',gammatr
    endif

    gamma2    = gamma**2
    gammatr2  = gammatr**2

    if (getdebug() > 1) then
       print *, '=============================================================='
       print *, ''
    endif

    modulationtype = 1 ! simpler and faster modulation

    !  call Set_Up(MY_RING)

    if (getdebug() > 0) then
       print *, 'Setting MADx with '
       print *, '    energy        ',energy
       print *, '    method        ',method0
       print *, '    Num. of steps ',nst0
       print *, '    charge        ',charge
    endif
    ! etienne helical
    helit=0
    call kanalnummer(mheli)
    open(unit=mheli,file='helical.txt',status='OLD',err=1001)
    read(mheli,*) helit
    if(helit>100) then
       write(6,*) " too many helical dipole ",helit
       stop 99
    endif
    do ihelit=1,helit
       read(mheli,*) heli(ihelit)
       CALL CONTEXT(heli(ihelit))
    enddo
    close(mheli)
1001 continue
    helit=0
    call kanalnummer(mheli)
    open(unit=mheli,file='sixtrack_compatible.txt',status='OLD',err=1002)
    read(mheli,*) sixtrack_compatible
    close(mheli)
1002 continue
    ! end of etienne helical

    ! preliminary setting
    !    my_ring%charge=1
    initial_charge=1
    CALL SET_MADx(energy=energy,METHOD=method0,STEP=nst0)
    if (getdebug() > 1) then
        print *, 'MADx is set'
    endif

    call resetclocks() ! clocks for modulation

    icav=0
    j=restart_sequ()
    nt=0
    l_machine=zero

    errors_out = get_value('ptc_create_layout ','errors_out ').ne.0
    magnet_name=" "
    if(errors_out) mg = get_string('ptc_create_layout ','magnet_name ',magnet_name)

    ! it is safe for MADX because default for all magnets is 1m
    ! so if user defines it otherwise it means it knows what he is doing
    absolute_aperture = 1e6_dp
    t_aperture = 1e30;
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  ELEMENTS LOOP    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

10  continue
    nst1=node_value("nst ")
    if(nst1.gt.0) then
       nstd = nst1
    else
       nstd = nst0
    endif

    ord_max = -1

    call zero_key(key)

    !j=j+1
    nt=nt+1
    if(nt==nt0) then
       call fort_warn("Potential problem for very large structure: ","More than 20'000 elements found")
    endif
    icount=0
    l=zero
    l=node_value('l ')
    key%list%l=l
    l_machine=l_machine+l
    code=node_value('mad8_type ')
    if(code.eq.39) code=15
    if(code.eq.38) code=24
    call element_name(name,name_len)
    key%list%name=name

    call node_name(name,name_len)
    key%list%vorname=name

    !frs&piotr 18 Dec 2007: sector_nmul must stay global for the time being
    !local, if present, superseed global at current node


    !*****************************
    !  MODEL Settings
    !*****************************

    model = node_value('model ')
    keymod1 = " "
    select case(model)
    CASE(1)
       keymod1 = "DRIFT_KICK       "
    CASE(2)
       keymod1 = "MATRIX_KICK      "
    CASE(3)
       keymod1 = "DELTA_MATRIX_KICK"
    END SELECT


    if(keymod1.ne." ") then
       key%model=keymod1
    else
       key%model=keymod0
    endif
    method1=node_value("method ")
    if(method1.eq.2.or.method1.eq.4.or.method1.eq.6) then
       metd = method1
    else
       metd = method0
    endif

    exact1=node_value("exact ")

    if(exact1.eq.0.or.exact1.eq.1) then
       EXACT_MODEL = exact1 .ne. 0
    else
       EXACT_MODEL = exact0
    endif

    !special node keys
    key%list%permfringe=node_value("fringe ") ! transfer(node_value("fringe ") .ne. zero, key%list%permfringe)
    key%list%bend_fringe=node_value("bend_fringe ") .ne. zero
    key%list%kill_ent_fringe=node_value("kill_ent_fringe ") .ne. zero
    key%list%kill_exi_fringe=node_value("kill_exi_fringe ") .ne. zero

    nn=name_len
    call node_string('apertype ',aptype,nn)
    APERTURE = zero
    nn = 0
    call get_node_vector('aperture ',nn,aperture)
    
    apoffset = zero
    napoffset = 0
    call get_node_vector('aper_offset ',napoffset,apoffset)
    
    if (getdebug() > 2) then
       print*,' Aperture type: >',aptype,'< ',nn,' parameters:'
       do i=1,nn
         print*,'             ',i,' : ',aperture(i)
       enddo
       print*,'          offset: napoffset=', napoffset
       do i=1,napoffset
         print*,'             ',i,' : ',apoffset(i)
       enddo
    endif
    
    
    
    !print*, name,'madx_ptc_module: Got for aperture nn=',nn, aperture(1), aperture(2)

    if(.not.((aptype.eq."circle".and.aperture(1).eq.zero).or.aptype.eq." ")) then

      c_%APERTURE_FLAG=.true.
      select case(aptype)
        case("circle")
          key%list%aperture_on=.true.
          key%list%aperture_kind=1
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(1)
          if (getdebug() > 2) then
             print*,' Aperture: circle with ',aperture(1),' radius'
          endif
        case("ellipse")
          key%list%aperture_on=.true.
          key%list%aperture_kind=1
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(2)
          if (getdebug() > 2) then
             print*,' Aperture: ellipse with ',aperture(1),' ',aperture(2),' radii'
          endif
        case("rectangle")
          key%list%aperture_on=.true.
          key%list%aperture_kind=2
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
          if (getdebug() > 2) then
             print*,' Aperture: rectangle with ',aperture(1),' ',aperture(2),' XY'
          endif
        case("lhcscreen") ! 2015-Mar-10  14:28:41  ghislain: added
          key%list%aperture_on=.true.
          key%list%aperture_kind=3
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
          key%list%aperture_r(1)=aperture(3)
          key%list%aperture_r(2)=aperture(3)
        case("rectcircle") ! 2015-Mar-10  14:28:41  ghislain: added
          key%list%aperture_on=.true.
          key%list%aperture_kind=3
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
          key%list%aperture_r(1)=aperture(3)
          key%list%aperture_r(2)=aperture(3)
        case("rectellipse")
          key%list%aperture_on=.true.
          key%list%aperture_kind=3
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
          key%list%aperture_r(1)=aperture(3)
          key%list%aperture_r(2)=aperture(4)
        case("racetrack") ! 2015-Mar-10  14:25:24  ghislain: generalized racetrack
          key%list%aperture_on=.true.
          key%list%aperture_kind=5
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
          key%list%aperture_r(1)=aperture(3)
          key%list%aperture_r(2)=aperture(4)
        case("octagon") ! 2015-Mar-10  14:25:37  ghislain: added octagon
          key%list%aperture_on=.true.
          key%list%aperture_kind=7
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
          key%list%aperture_r(1)=aperture(3)
          key%list%aperture_r(2)=aperture(4)

          write(whymsg,*) 'Aperture: <<',aptype,'>> at magnet ',name(:len_trim(name)),' is not implemented by PTC'
          call fort_warn('ptc_createlayout: ',whymsg(:len_trim(whymsg)))
          call aafail('ptc_input:','Aperture type not implemented. Program stops')
          
        case("general") ! 2015-Mar-10  14:25:48  ghislain: kind was 6
          key%list%aperture_kind=8
          print*,"General aperture not implemented"
          call aafail('ptc_input:','General aperture not implemented. Program stops')
        
        case DEFAULT
          
          ! in case aperture is defined as file with arbitrary polygon points
          i = get_userdefined_geometry_len()
          if (i > 0) then
            if(nn > 1) then
              key%list%aperture_x=aperture(1)
              key%list%aperture_y=aperture(2)
            else
              key%list%aperture_x=0
              key%list%aperture_y=0
            endif
             
            allocate(apx(i))
            allocate(apy(i))
            
            i = get_userdefined_geometry(apx,apy,i)
            
            key%list%APERTURE_POLYGX => apx
            key%list%APERTURE_POLYGY => apy
            
            key%list%aperture_on=.true.
            key%list%aperture_kind=6
            
            if (getdebug()>1)  then
              print*, "Aperture defined as a polygon with ", i, " points "
            endif
          
          else
          
          
            write(whymsg,*) 'Aperture: <<',aptype,'>> at magnet ',name(:len_trim(name)),' is not recognized by PTC'
            call fort_warn('ptc_createlayout: ',whymsg(:len_trim(whymsg)))
            call aafail('ptc_input:','Aperture type not implemented. Program stops')
          endif
       end select
       
       
       key%list%aperture_dx=apoffset(1)
       key%list%aperture_dy=apoffset(2)
       
       
  !  else
  !   if( .not. ((code.eq.1) .or. (code.eq.4)) ) then
  !     write(*,'(a10,1x,a16,1x,a14,1x,6f10.6)') 'Aperture: ',aptype(1:16),'aperture pars:', aperture(1:6)
  !     write(whymsg,*) 'Aperture: ',aptype,' at magnet ',name(:len_trim(name)),' not supported by PTC'
  !     call fort_warn('ptc_createlayout: ',whymsg(:len_trim(whymsg)))
  !
  !   endif

    endif
    call append_empty(my_ring)

    if (getdebug() > 2) then
       print *,'Element ',key%list%name,' code is ',code
    endif


    select case(code)
    case(0,25)
       key%magnet="marker"
    case(4)
       call aafail('ptc_input:','PTC does not accept matrix elements. Program stops.')
    case(22)
       key%magnet="marker"
    !case(1,11)
    case(1,20,21,44)
       key%magnet="drift"
       CALL CONTEXT(key%list%name)

       do ihelit=1,helit
          IF(index(key%list%name,heli(ihelit)(1:len_trim(heli(ihelit))))/=0) then
             key%magnet="helicaldipole"
             write(6,*) " drift ",key%list%name, " became helical dipole in PTC "
          endif
       enddo
       ! end etienne Helical
    case(2) ! PTC accepts mults
       if(l.eq.zero) then
          key%magnet="marker"
          goto 100
       endif
       key%magnet="rbend"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       tempdp=sqrt(normal_0123(0)*normal_0123(0)+skew_0123(0)*skew_0123(0))
       key%list%b0=bvk*(node_value('angle ')+tempdp*l)

       !       print*, "RBEND: Angle: ", node_value('angle ')," tempdp ", tempdp, " l ", l
       !       print*, "RBEND: normal: ",normal_0123(0)," skew: ",skew_0123(0)

       key%list%k(2)=node_value('k1 ')+ key%list%k(2)
       key%list%k(3)=node_value('k2 ')+ key%list%k(3)
       key%list%k(4)=node_value('k3 ')+ key%list%k(4)

       key%list%ks(2)=node_value('k1s ')+ key%list%ks(2)
       key%list%ks(3)=node_value('k2s ')+ key%list%ks(3)
       key%list%ks(4)=node_value('k3s ')+ key%list%ks(4)

       if(EXACT_MODEL.and.(node_value('angle ').eq.zero)) then
          key%magnet="quadrupole"
          key%tiltd=node_value('tilt ')
       else

          ! Gymnastic needed since PTC expects MAD8 convention
          key%list%t1=node_value('e1 ')
          key%list%t2=node_value('e2 ')
          key%list%hgap=node_value('hgap ')
          !       key%list%fint=node_value('fint ')
          fint=node_value('fint ')
          fintx=node_value('fintx ')
          if((fintx.ne.fint).and.(fintx.gt.zero.and.fint.gt.zero)) then
             print*," The fint and fintx must be the same at each end or each might be zero"
             call aafail('ptc_input:','The fint and fintx must be the same at each end or each might be zero. Program stops.')
          endif
          if(fint.gt.zero) then
             key%list%fint=fint
             if(fintx.eq.zero) key%list%kill_exi_fringe=my_true
          else
             if(fintx.gt.zero) then
                key%list%fint=fintx
                key%list%kill_ent_fringe=my_true
             else
                key%list%fint=zero
             endif
          endif
          key%list%h1=node_value('h1 ')
          key%list%h2=node_value('h2 ')
          key%tiltd=node_value('tilt ')
          if(tempdp.gt.0) key%tiltd=key%tiltd + atan2(skew_0123(0),normal_0123(0))
          ptcrbend=node_value('ptcrbend ').ne.0
          if(ptcrbend) then
             call context(key%list%name)
             truerbend=node_value('truerbend ').ne.0
             if(truerbend) then
                key%magnet="TRUERBEND"
                if(key%list%t2/=zero) then
                   write(6,*) " The true parallel face bend "
                   write(6,*) " only accepts the total angle and e1 as an input "
                   write(6,*) " if e1=0, then the pipe angle to the entrance face is "
                   write(6,*) " angle/2. It is a normal rbend."
                   write(6,*) " If e1/=0, then the pipe angle to the entrance face is "
                   write(6,*) ' angle/2+e1 and the exit pipe makes an angle "angle/2-e1" '
                   write(6,*) " with the exit face."
                   write(6,*) " The offending non-zero t2 = (e2 - angle/2) is set to zero! "
                   write(6,*) " Make sure that this is what you want!!! "
                   !                write(6,*) " CHANGE YOUR LATTICE FILRE."
                   !                stop 666
                   key%list%t2=zero
                endif
             else
                key%magnet="WEDGRBEND"
             endif
          endif
       endif
       if(errors_out) then
          if(key%list%name(:len_trim(magnet_name)-1).eq. &
               magnet_name(:len_trim(magnet_name)-1)) then
             call string_to_table_curr('errors_dipole ', 'name ', key%list%name)
             call double_to_table_curr('errors_dipole ', 'k0l ', bvk*key%list%b0)
             call augment_count('errors_dipole ')
          endif
       endif
    case(3) ! PTC accepts mults watch out sector_nmul defaulted to 22

       if (getdebug()>2) print*,"Translating SBEND"
       if(l.eq.zero) then
          if (getdebug()>2) print*,"Length zero -> translating as MARKER"
          key%magnet="marker"
          goto 100
       endif
       key%magnet="sbend"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)


       if( (sector_nmul_max.lt.ord_max) .and. EXACT_MODEL) then

         sector_nmul_max =  ord_max;
         sector_nmul =  ord_max;

         call aafail('the order of multipoles in a sbend in exact mode cannot be ',&
            &'larger than sector_mul_max: check your ptc_create_universe input')

       endif

       tempdp=sqrt(normal_0123(0)*normal_0123(0)+skew_0123(0)*skew_0123(0))
       key%list%b0=bvk*(node_value('angle ')+ tempdp*l)

       key%list%k(2)=node_value('k1 ')+ key%list%k(2)
       key%list%k(3)=node_value('k2 ')+ key%list%k(3)
       key%list%k(4)=node_value('k3 ')+ key%list%k(4)

       key%list%ks(2)=node_value('k1s ')+ key%list%ks(2)
       key%list%ks(3)=node_value('k2s ')+ key%list%ks(3)
       key%list%ks(4)=node_value('k3s ')+ key%list%ks(4)

       key%list%t1=node_value('e1 ')
       key%list%t2=node_value('e2 ')
       key%list%hgap=node_value('hgap ')
       !       key%list%fint=node_value('fint ')
       fint=node_value('fint ')
       fintx=node_value('fintx ')
       if((fintx.ne.fint).and.(fintx.gt.zero.and.fint.gt.zero)) then
          print*," The fint and fintx must be the same at each end or each might be zero"
          call aafail('ptc_input:','The fint and fintx must be the same at each end or each might be zero. Program stops')
       endif
       if(fint.gt.zero) then
          key%list%fint=fint
          if(fintx.eq.zero) key%list%kill_exi_fringe=my_true
       else
          if(fintx.gt.zero) then
             key%list%fint=fintx
             key%list%kill_ent_fringe=my_true
          else
             key%list%fint=zero
          endif
       endif
       key%list%h1=node_value('h1 ')
       key%list%h2=node_value('h2 ')
       key%tiltd=node_value('tilt ')
       if(tempdp.gt.0) key%tiltd=key%tiltd + atan2(skew_0123(0),normal_0123(0))
       if(errors_out) then
          if(key%list%name(:len_trim(magnet_name)-1).eq. &
               magnet_name(:len_trim(magnet_name)-1)) then
             call string_to_table_curr('errors_dipole ', 'name ', key%list%name)
             call double_to_table_curr('errors_dipole ', 'k0l ', bvk*key%list%b0)
             call augment_count('errors_dipole ')
          endif
       endif

       if (getdebug()>2) then
         print*,"B0=", key%list%b0
         print*,"K=", key%list%k
         print*,"KS=", key%list%ks
         print*,"TILT=", key%tiltd
         print*,"T1=", key%list%t1 ! e1
         print*,"T2=", key%list%t2 ! e2
         print*,"H1=", key%list%h1
         print*,"H2=", key%list%h2
       endif

    case(5)
       key%magnet="quadrupole"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       ! Read data & fill %k(:), %ks(:) arrays which are
       ! summs of multipoles and errors

! LD: 19.06.2019
       sk0=node_value('k0 ')

       ! quadrupole components
       sk1= node_value('k1 ')+node_value('k1tap ')
       sk1s=node_value('k1s ')
       tilt=node_value('tilt ')
       dum1=key%list%k(2)-normal_0123(1)
       dum2=key%list%ks(2)-skew_0123(1)

       !print*,'normal_0123', normal_0123
       !print*,'skew_0123', skew_0123
       !print*,'sk1 sk1s dum1 dum2'
       !print*, sk1, sk1s, dum1, dum2

! LD: 19.06.2019
!       if(dum1.ne.zero.or.dum2.ne.zero) then                      !
          sk1= sk1 +dum1                                          !
          sk1s=sk1s+dum2                                          !
!       endif                                                      !

! LD: 19.06.2019
       if (sk1s .ne. zero) then
          if (ord_max .le. 2 .and. sk0 .eq. 0 .and. key%list%permfringe .eq. 0) then !
            tilt = -atan2(sk1s, sk1)/two + tilt
            sk1  = sqrt(sk1**2 + sk1s**2)
            sk1s = zero
          elseif (metd .lt. 4 .and. model .ne. 1) then
            call fort_warn('quadrupole with k1s and k0, k2, k2s or permfringe detected: ',&
                           'use method=4 or 6 for better results with model=2')
          endif
!          key%list%k(2) =sk1
!          key%list%ks(2)=zero
!          key%tiltd=tilt
       endif                                                      !
       key%list%k(2) =sk1                                         !
       key%list%ks(2) =sk1s                                        !
!       key%list%ks(2)=zero  ! added by VK                         !
       key%tiltd=tilt  !==========================================!

       !================================================================
       ! dipole component not active in MAD-X proper
! LD: 19.06.2019
       key%list%k(1)=key%list%k(1)+bvk*sk0

    case(6)
       key%magnet="sextupole"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       ! sextupole components
       sk2= node_value('k2 ') + node_value('k2tap ')
       sk2s=node_value('k2s ')
       tilt=node_value('tilt ')
       dum1=key%list%k(3)-normal_0123(2)
       dum2=key%list%ks(3)-skew_0123(2)

! LD: 19.06.2019
!       if(dum1.ne.zero.or.dum2.ne.zero) then                      !
          sk2= sk2 +dum1                                          !
          sk2s=sk2s+dum2                                          !
!       endif                                                      !
!       if (sk2s .ne. zero) then                                   !
!          tilt = -atan2(sk2s, sk2)/three + tilt                   !
!          sk2 =  sqrt(sk2**2 + sk2s**2)                           !
!       endif                                                      !
       key%list%k(3) =sk2                                         !
       key%list%ks(3) =sk2s                                         !
!       key%list%ks(3)=zero  ! added by VK                         !
       key%tiltd=tilt  !==========================================!

       !================================================================
       if(errors_out) then
          if(key%list%name(:len_trim(magnet_name)-1).eq. &
               magnet_name(:len_trim(magnet_name)-1)) then
             call string_to_table_curr('errors_total ', 'name ', key%list%name)
             myfield(:) = zero
             do kk=1,maxmul
                myfield(2*kk-1) = key%list%k(kk)
                myfield(2*kk)   = key%list%ks(kk)
             enddo
             call vector_to_table_curr('errors_total ', 'k0l ', myfield(1), i)
             call augment_count('errors_total ')
          endif
       endif

    case(7)
       key%magnet="octupole"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       ! octupole components
       sk3= node_value('k3 ')
       sk3s=node_value('k3s ')
       tilt=node_value('tilt ')
       dum1=key%list%k(4)-normal_0123(3)
       dum2=key%list%ks(4)-skew_0123(3)

! ! LD: 19.06.2019
!       if(dum1.ne.zero.or.dum2.ne.zero) then                      !
          sk3= sk3 +dum1                                          !
          sk3s=sk3s+dum2                                          !
!       endif                                                      !
!       if (sk3s .ne. zero) then                                   !
!          tilt = -atan2(sk3s, sk3)/four + tilt                    !
!          sk3 = sqrt(sk3**2 + sk3s**2)                            !
!       endif                                                      !
       key%list%k(4) =sk3                                         !
       key%list%ks(4)= sk3s
!       key%list%ks(4)=zero  ! added by VK                         !
       key%tiltd=tilt  !==========================================!

       !================================================================

    case(8)
       if (getdebug()>2) print*,"Translating MULTIPOLE"
       key%magnet="multipole"
       !---- Multipole components.
       F_ERRORS = zero
       n_ferr = node_fd_errors(f_errors)
       NORMAL = zero
       SKEW = zero
       call get_node_vector('knl ',nn,normal)
       call get_node_vector('ksl ',ns,skew)
       if(nn.ge.NMAX) nn=NMAX-1
       if(ns.ge.NMAX) ns=NMAX-1
       do i=1,NMAX
          key%list%k(i)=zero
          key%list%ks(i)=zero
       enddo
       skew(0)=-skew(0) ! frs error found 30.08.2008

       key%list%thin_h_angle=bvk*normal(0)
       key%list%thin_v_angle=bvk*skew(0)
       lrad=node_value('lrad ')
       if(lrad.gt.zero) then
          key%list%thin_h_foc=normal(0)*normal(0)/lrad
          key%list%thin_v_foc=skew(0)*skew(0)/lrad
       endif


       if(nn.gt.0) then
         ! should be like this but apparently for MADX compatibility it was left out
         ! key%list%k(1)=normal(0)

          do i=1,nn
             !print*, "multipole normal ", i, " = ", normal(i)
             key%list%k(i+1)=normal(i)
          enddo
       endif

       if(ns.gt.0) then
         ! should be like this but apparently for MADX compatibility it was left out
         ! key%list%ks(1)=skew(0)

          do i=1,ns
             !print*, "multipole skew ", i, " = ", skew(i)
             key%list%ks(i+1)=skew(i)
          enddo
       endif
       FIELD = zero
       if (n_ferr .gt. 0) then
          call dcopy(f_errors,field,n_ferr)
       endif
       nd = max(nn, ns, n_ferr/2)


       if(nd.ge.maxmul) nd=maxmul-1
       if(n_ferr.gt.0) then
          do i=0,nd
             key%list%k(i+1)=key%list%k(i+1)+field(1,i)
             key%list%ks(i+1)=key%list%ks(i+1)+field(2,i)
          enddo
       endif
       key%tiltd=node_value('tilt ')
       if(errors_out) then
          if(key%list%name(:len_trim(magnet_name)-1).eq. &
               magnet_name(:len_trim(magnet_name)-1)) then
             call string_to_table_curr('errors_field ', 'name ', key%list%name)
             call string_to_table_curr('errors_total ', 'name ', key%list%name)
             i=2*maxmul+2
             myfield(:) = zero
             do kk=1,nd+1
                myfield(2*kk-1) = field(1,kk-1)
                myfield(2*kk)   = field(2,kk-1)
             enddo
             call vector_to_table_curr('errors_field ', 'k0l ', myfield(1), i)
             myfield(:) = zero
             do kk=1,nd+1
                myfield(2*kk-1) = key%list%k(kk)
                myfield(2*kk)   = key%list%ks(kk)
             enddo
             call vector_to_table_curr('errors_total ', 'k0l ', myfield(1), i)
             call augment_count('errors_field ')
             call augment_count('errors_total ')
          endif
       endif

       if (getdebug()>2) then
         print*,"thin_h_angle=", key%list%thin_h_angle
         print*,"thin_v_angle=", key%list%thin_v_angle
         print*,"K=", key%list%k
         print*,"KS=", key%list%ks
         print*,"TILT=", key%tiltd
       endif


    case(9) ! PTC accepts mults
       key%magnet="solenoid"
       ks=node_value('ks ')
       if(l.ne.zero) then
          key%list%bsol=bvk*ks
       else
          ksi=node_value('ksi ')
          lrad=node_value('lrad ')
          if(lrad.eq.zero.and.ks.ne.zero) lrad=ksi/ks
          if(ksi.eq.zero.or.lrad.eq.zero) then
             key%magnet="marker"
             print*,"Thin solenoid: ",name," has no strength - set to marker"
          else
             key%list%bsol=bvk*ksi/lrad
             key%list%ls=lrad
          endif
       endif
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

    case(10)
       key%magnet="rfcavity"
       key%list%volt=bvk*node_value('volt ')
       freq=c_1d6*node_value('freq ')

       key%list%lag = -(node_value('lag ')+node_value('lagtap ') )*twopi


       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       if(offset_deltap.ne.zero) then

          default = getintstate()
          default=default+totalpath0
          call setintstate(default)
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one) !regular twiss values

       endif
       key%list%freq0=freq
       key%list%n_bessel=node_value('n_bessel ')
       key%list%harmon=one ! it is ignored by PTC because it does not know the circumference

       if(key%list%volt.ne.zero.and.key%list%freq0.ne.zero) then
         icav=1
         if (getdebug() > 2) then
            print*,"icav set to 1, RF Cavity detected with Volt ",key%list%volt, " and Freq ", key%list%freq0
         endif
       !else
       !  if (getdebug() > 2) then
       !     print*,"RF Cavity with zero voltage or frequency"
       !  endif
       endif

       !  case(11)
       !     key%magnet="elseparator"
       !     key%list%volt=node_value('ex ')
       !     key%list%lag=atan2(node_value('ey '),node_value('ex '))
       !     key%tiltd=node_value('tilt ')

       m_u%end%HARMONIC_NUMBER=node_value('harmon ')   ! etienne_harmon
       no_cavity_totalpath=node_value('no_cavity_totalpath ').ne.0
       if(no_cavity_totalpath) then
          key%list%cavity_totalpath=0
       else
          key%list%cavity_totalpath=1
          ! correction for time of flight through cavity
          ! we want particle with t=0 to be not accelerated
          key%list%lag = key%list%lag + twopi*freq*(l/2d0)/(clight*beta0)
       endif

       modulationq = node_value('modulationq ')
       if (abs(modulationq) .gt. 1e-12) then
         
         key%list%clockno_ac = getclockidx(modulationq)

         if (key%list%clockno_ac .lt. 0) then
           call aafail('ptc_input:', &
           'Too many AC Dipole clocks, PTC can accept max 3 clocks with given tune and ramp. Program stops.')
         endif

       
         key%list%n_ac = 1 
         
         key%list%d_volt = node_value('volterr ')
         key%list%d_phas = node_value('lagerr ')

         !print*,"RF Cavity modulation ON volt", key%list%d_volt, " lag ", key%list%d_phas
       endif
       !else
       !  print*,"RF Cavity modulation OFF"
       !endif

       ! LD: 09.04.2019
!       write (*,'(3(a,E25.16))') "@@ RF freq= ", freq," lag= ", key%list%lag, " lag= ", node_value('lag ')

!       print*,"madx_ptc_module::input volt: ", key%list%volt, &
!                                    " lag : ", key%list%lag, &
!                                    " harm: ", key%list%harmon, &
!                                    " freq: ", key%list%freq0

    case(11) ! LD: 04.07.2019
      key%magnet="elseparator"
      ex = node_value('ex ')
      ey = node_value('ey ')
      if (l .ne. 0) then
        ex = ex + node_value('ex_l ')/l
        ey = ey + node_value('ey_l ')/l
      endif
      key%list%volt=sqrt(ex**2 + ey**2)
      key%list%lag=atan2(ey,ex)

    case(12)
       ! actually our SROT element
       key%magnet="CHANGEREF"
       PATCH_ANG = zero
       PATCH_TRANS = zero
       patch_ang(3)=node_value('angle ')
       key%list%patchg=2
       do i=1,3
          key%list%ang(i)=patch_ang(i)
          key%list%t(i)=patch_trans(i)
       enddo
    case(13)
       ! actually our YROT element
       key%magnet="CHANGEREF"
       PATCH_ANG = zero
       PATCH_TRANS = zero
       patch_ang(2)=node_value('angle ')
       key%list%patchg=2
       do i=1,3
          key%list%ang(i)=patch_ang(i)
          key%list%t(i)=patch_trans(i)
       enddo
    case(14,15,16) ! PTC accepts mults
       ! kickers (corrector magnets)
       F_ERRORS = zero
       n_ferr = node_fd_errors(f_errors)
       do i=1,NMAX
          key%list%k(i)=zero
          key%list%ks(i)=zero
       enddo
       do i = 1, 2
          fieldk(i) = zero
       enddo
       if (n_ferr .gt. 0) call dcopy(f_errors, fieldk, min(2, n_ferr))
       if (l .eq. zero)  then
          div = one
       else
          div = l
       endif
       if(code.eq.14) then
          key%magnet="hkicker"
          key%list%k(1)=(node_value('kick ')+node_value('chkick ')+fieldk(1)/div)
      else if(code.eq.15) then
          key%magnet="kicker"
          key%list%k(1)=(node_value('hkick ')+node_value('chkick ')+fieldk(1)/div)
          key%list%ks(1)=(node_value('vkick ')+node_value('cvkick ')+fieldk(2)/div)
       else if(code.eq.16) then
          key%magnet="vkicker"
          key%list%ks(1)=(node_value('kick ')+node_value('cvkick ')+fieldk(2)/div)
       else
          key%magnet="marker"
       endif
       i=2*maxmul+2
       if(errors_out) then
          if(key%list%name(:len_trim(magnet_name)-1).eq. &
               magnet_name(:len_trim(magnet_name)-1)) then
             myfield(:) = zero
             myfield(1)=-key%list%k(1)
             myfield(2)=key%list%ks(1)
             call string_to_table_curr('errors_total ', 'name ', key%list%name)
             call vector_to_table_curr('errors_total ', 'k0l ', myfield(1), i)
             call augment_count('errors_total ')
          endif
       endif
       key%tiltd=node_value('tilt ')
    case(17)
       key%magnet="hmonitor"
    case(18)
       key%magnet="monitor"
    case(19)
       key%magnet="vmonitor"
       !  2015-Mar-05  13:13:35  ghislain: Warning !
       !                         ecollimator and rcollimator replaced by collimator in MAD-X, with code 44
       !  case(20)
       !     key%magnet="ecollimator"
       !     key%list%x_col=node_value('xsize ')
       !     key%list%y_col=node_value('ysize ')
       !     key%tiltd=node_value('tilt ')
       !  case(21)
       !     key%magnet="rcollimator"
       !     key%list%x_col=node_value('xsize ')
       !     key%list%y_col=node_value('ysize ')
       !     key%tiltd=node_value('tilt ')
    !case(20,21,44)
    !        key%magnet="madcollimator"
    !        key%list%x_col=1e3
    !        key%list%y_col=1e3
    !        key%tiltd=node_value('tilt ')
    case(33)
       !---- This is the dipedge element
       edge= node_value('e1 ')
       hgap= node_value('hgap ')
       rhoi= bvk * node_value('h ')
       fint= node_value('fint ')
       corr= 2 * rhoi * hgap * fint
       if(rhoi .ne. zero .and. ( edge .ne. zero .or. corr .ne. zero )) then
          key%magnet="multipole"
          tanedg = tan(edge)
          secedg = one / cos(edge)
          psip = edge - corr * secedg * (one + sin(edge)**2)
          key%list%hf= rhoi * tanedg
          key%list%vf= -rhoi * tan(psip)
          key%tiltd=node_value('tilt ') !! frs add-on
       else
          key%magnet="marker"
       endif
    case(24)
       key%magnet="instrument"
       key%tiltd=node_value('tilt ')
    case(27)
       key%magnet="twcavity"
       key%list%volt=bvk*node_value('volt ')
       freq=c_1d6*node_value('freq ')
       key%list%lag=-node_value('lag ')*twopi
       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       default=default+totalpath0 !fringe field calculation vitally relies on it!!!!
       if(offset_deltap.ne.zero) then
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one)
       endif
       key%list%freq0=freq
       key%list%dphas=node_value("delta_lag ")
       key%list%psi=node_value("psi ")
       key%list%harmon=one
       if(key%list%volt.ne.zero.and.key%list%freq0.ne.zero) icav=1
    case(34) ! XROTATION
       key%magnet="CHANGEREF"
       PATCH_ANG = zero
       PATCH_TRANS = zero
       patch_ang(1)=node_value('angle ')
       key%list%patchg=2
       do i=1,3
          key%list%ang(i)=patch_ang(i)
          key%list%t(i)=patch_trans(i)
       enddo
    case(35)
       key%magnet="CHANGEREF"
       PATCH_ANG = zero
       PATCH_TRANS = zero
       call get_node_vector('patch_ang ',3,patch_ang)
       call get_node_vector('patch_trans ',3,patch_trans)
       key%list%patchg=2
       do i=1,3
          key%list%ang(i)=patch_ang(i)
          key%list%t(i)=patch_trans(i)
       enddo
    case(36) ! TRANSLATION
       key%magnet="CHANGEREF"
       PATCH_ANG = zero
       patch_trans(1)=node_value('dx ')
       patch_trans(2)=node_value('dy ')
       patch_trans(3)=node_value('ds ')
       key%list%patchg=2
       do i=1,3
          key%list%ang(i)=patch_ang(i)
          key%list%t(i)=patch_trans(i)
       enddo
    case(37)!CRAB ??
       key%magnet="rfcavity"
       key%list%volt=zero
       do i=1,NMAX
          key%list%k(i)=zero
          key%list%ks(i)=zero
       enddo
       key%list%k(1)=node_value('volt ')*c_1d_3
       ! vertical crab
       ! maybe requires a flip of sign
       !       key%list%ks(1)= (+/-)  node_value('volt ')*c_1d_3
       !
       freq=c_1d6*node_value('freq ')
       key%list%lag=-node_value('lag ')*twopi+pih
       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       if(offset_deltap.ne.zero) then
          default = getintstate()
          default=default+totalpath0
          call setintstate(default)
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one)
       endif
       key%list%freq0=freq
       key%list%n_bessel=0
       key%list%harmon=one

       if(key%list%k(1).ne.zero.and.key%list%freq0.ne.zero) icav=1

    !RFMULTIPOLE, crab also falls here, but is made with special case where volt defines BN(1)

    case(40)

       key%magnet="hkicker"
       do i=1,NMAX
          key%list%k(i)=zero
          key%list%ks(i)=zero
       enddo

        key%list%n_ac = 1 ! only dipole
        ! need to convert voltage (E field) to corresponding B field
        if (L .gt. 0) then
          key%list%d_bn(1) =  0.3 * node_value('volt ')  / ( L * beta0 * get_value('beam ','pc '))
        else
          key%list%d_bn(1) =  0.3 * node_value('volt ')  / (beta0 * get_value('beam ','pc '))
        endif

        if (getdebug() > 1) then
          print*,"HACD bn(1)=", key%list%d_bn(1), "b0=",beta0, " pc=",get_value('beam ','pc '), " L=",L
        endif
        
        key%list%d_an(1) = zero

        key%list%D_ac = one ! extrac factor for amplitude; we use it for ramping

        ! parameters to modulate the nominal parameters. No modulation in MADX implemented.
        key%list%DC_ac = zero
        key%list%A_ac = zero
        key%list%theta_ac = -node_value('lag ') ! it is ignored with fast modulationtype = 1

        modulationq = node_value('freq ')
        key%list%clockno_ac = getclockidx(modulationq)

        if (key%list%clockno_ac .lt. 0) then
          call aafail('ptc_input:', &
          'Too many AC Dipole clocks, PTC can accept max 2 clocks with given tune and ramp. Program stops.')
        endif


    case(41)

       key%magnet="vkicker"
       do i=1,NMAX
          key%list%k(i)=zero
          key%list%ks(i)=zero
       enddo

        key%list%n_ac = 1 ! only dipole
        if (L .gt. 0) then
          key%list%d_an(1) =  0.3 * node_value('volt ') / ( L * beta0 * get_value('beam ','pc '))
        else
          key%list%d_an(1) =  0.3 * node_value('volt ')  / (beta0 * get_value('beam ','pc '))
        endif
        
        if (getdebug() > 1) then
          print*,"VACD bn(1)=", key%list%d_an(1), "b0=",beta0, " pc=",get_value('beam ','pc '), " L=",L
        endif
        
        key%list%d_bn(1) = zero

        key%list%D_ac = one ! extrac factor for amplitude; we use it for ramping

        ! parameters to modulate the nominal parameters. No modulation in MADX implemented.
        key%list%DC_ac = zero
        key%list%A_ac = zero
        key%list%theta_ac = -node_value('lag ')

        modulationq = node_value('freq ')
        key%list%clockno_ac = getclockidx(modulationq)

        if (key%list%clockno_ac .lt. 0) then
          call aafail('ptc_input:', &
          'Too many AC Dipole clocks, PTC can accept max 2 clocks with given tune and ramp. Program stops.')
        endif


    case(43)
       key%magnet="rfcavity" ! RFMULTIPOLE
       key%list%volt=bvk*node_value('volt ')
       freq=c_1d6*node_value('freq ')
       key%list%lag=-node_value('lag ')*twopi

       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       if(offset_deltap.ne.zero) then

          default = getintstate()
          default=default+totalpath0
          call setintstate(default)
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one) !regular twiss values

       endif
       key%list%freq0=freq
       key%list%n_bessel=node_value('n_bessel ')
       key%list%harmon=one ! it is ignored by PTC because it does not know the circumference

       key%list%k(:)=zero
       key%list%ks(:)=zero

       F_ERRORS=zero !call dzero(f_errors,maxferr+1)
       n_ferr = node_fd_errors(f_errors)
       NORMAL=zero !call dzero(normal,maxmul+1)
       SKEW=zero !call dzero(skew,maxmul+1)
       call get_node_vector('knl ',nn,normal)
       call get_node_vector('ksl ',ns,skew)
       if(nn.ge.NMAX) nn=NMAX-1
       if(ns.ge.NMAX) ns=NMAX-1

       skew(0)=-skew(0) ! frs error found 30.08.2008
       key%list%thin_h_angle=bvk*normal(0)
       key%list%thin_v_angle=bvk*skew(0)
       lrad=node_value('lrad ')
       if(lrad.gt.zero) then
          key%list%thin_h_foc=normal(0)*normal(0)/lrad
          key%list%thin_v_foc=skew(0)*skew(0)/lrad
       endif

       do i=0,nn
          key%list%k(i+1)=normal(i)
          if (normal(i) /= zero) icav=1
       enddo

       do i=0,ns
          key%list%ks(i+1)=skew(i)
          if (normal(i) /= zero) icav=1
       enddo

       FIELD=zero !call dzero(field,2*(maxmul+1))
       if (n_ferr .gt. 0) then
          call dcopy(f_errors,field,n_ferr)
       endif
       nd = max(nn, ns, n_ferr/2)


       if(nd.ge.maxmul) nd=maxmul-1
       if(n_ferr.gt.0) then
          do i=0,nd
             key%list%k(i+1)=key%list%k(i+1)+field(1,i)
             key%list%ks(i+1)=key%list%ks(i+1)+field(2,i)
          enddo
       endif

       key%tiltd=node_value('tilt ')
       if(errors_out) then
          if(key%list%name(:len_trim(magnet_name)-1).eq. &
               magnet_name(:len_trim(magnet_name)-1)) then
             call string_to_table_curr('errors_field ', 'name ', key%list%name)
             call string_to_table_curr('errors_total ', 'name ', key%list%name)
             i=2*maxmul+2
             myfield(:) = zero
             do kk=1,nd+1
                myfield(2*kk-1) = field(1,kk-1)
                myfield(2*kk)   = field(2,kk-1)
             enddo
             call vector_to_table_curr('errors_field ', 'k0l ', myfield(1), i)
             myfield(:) = zero
             do kk=1,nd+1
                myfield(2*kk-1) = key%list%k(kk)
                myfield(2*kk)   = key%list%ks(kk)
             enddo
             call vector_to_table_curr('errors_total ', 'k0l ', myfield(1), i)
             call augment_count('errors_field ')
             call augment_count('errors_total ')
          endif
       endif


       if(key%list%volt.ne.zero.and.key%list%freq0.ne.zero) icav=1

       m_u%end%HARMONIC_NUMBER=node_value('harmon ')   ! etienne_harmon
       no_cavity_totalpath=node_value('no_cavity_totalpath ').ne.0
       if(no_cavity_totalpath) then
          key%list%cavity_totalpath=0
       else
          key%list%cavity_totalpath=1
       endif

!       print*,"madx_ptc_module::input volt: ", key%list%volt, &
!                                    " lag : ", key%list%lag, &
!                                    " harm: ", key%list%harmon, &
!                                    " freq: ", key%list%freq0



    case default
       print*,"Element: ",name, "of type ",code," not implemented in PTC"
       call aafail('ptc_input:','Element not implemented in PTC. Program stops')
    end select
100 continue

    !apply BVK f;ag
    do i=1,NMAX
       key%list%k(i)=bvk*key%list%k(i)
       key%list%ks(i)=bvk*key%list%ks(i)
    enddo

   ! if (ord_max > 0) then
   !    key%list%nmul = ord_max
   ! endif

    call create_fibre(my_ring%end,key,EXCEPTION) !in ../libs/ptc/src/Sp_keywords.f90

    if(key%list%n_ac > 0 ) then
      !save pointer to the AC dipole element for ramping in tracking
      if (getdebug() > 1) then
         print*,"Adding Modulated Element: ",name, " of type ",code," to clock ",key%list%clockno_ac
      endif
      
      call addelementtoclock(my_ring%end,key%list%clockno_ac)
    endif


    if(advance_node().ne.0)  goto 10


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! END OF ELEMENTS LOOP    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (getdebug() > 0) then
       print*,' Length of machine: ',l_machine
    endif
    
    CALL GET_ENERGY(ENERGY,kin,BRHO,beta0,P0C)
    beta0start = beta0
    
    isclosedlayout=get_value('ptc_create_layout ','closed_layout ') .ne. 0

    if (getdebug() > 0) then
       if ( isclosedlayout .eqv. .true. ) then
          print *,'The machine is a RING'
       else
          print *,'The machine is a LINE'
       endif
    endif

    MY_RING%closed=isclosedlayout

    doneit=.true.
    call ring_l(my_ring,doneit)

    resplit=get_value('ptc_create_layout ','resplit ').ne.0
    if(resplit) then
       my_thin = get_value('ptc_create_layout ','thin ')
       my_xbend = get_value('ptc_create_layout ','xbend ')
       even = get_value('ptc_create_layout ','even ').ne.0
       resplit_cutting=2
       CALL THIN_LENS_resplit(my_ring,THIN=my_thin,even=even,xbend=my_xbend)
    endif

    if (getdebug() > 0) then
       write(6,*) "------------------------------------ PTC Survey ------------------------------------"
       write(6,*) "Before start: ",my_ring%start%chart%f%a
       write(6,*) "Before   end: ",my_ring%end%chart%f%b
    endif

    call survey(my_ring)

    if (getdebug() > 0) then
       write(6,*) "After  start: ",my_ring%start%chart%f%a
       write(6,*) "After    end: ",my_ring%end%chart%f%b
    endif

    call setintstate(default)
    
    call get_length(my_ring,l)
    my_ring_length = l
    if(my_ring%HARMONIC_NUMBER>0) then
       print*,"HARMONIC NUMBER defined in the ring: ", my_ring%HARMONIC_NUMBER
       

       j=restart_sequ()
       p=>my_ring%start

       do i=1,my_ring%n
          if(p%mag%kind==kind4) then
             if(p%mag%freq==zero) then

                tempdp = node_value('harmon ')

                p%mag%freq=clight*tempdp*BETA0/l
                p%magp%freq=p%mag%freq

                ! watch the msg buffer is 1024
                write(msg,*) " Cavity ",p%mag%name," defined with harmonic number ",tempdp,". Using SUM(LD) as ring length: ", l, &
                             " instead of real orbit length. Obtained freq. = ",p%mag%freq," Hz"

                call fort_warn("ptc_input",msg(:len_trim(msg)))
                if (p%mag%volt .ne. zero) icav=1
             endif
          endif
          p=>p%next
          j = advance_node()
       enddo
    endif

    if (getdebug() > 1) then
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^    F I N I S H E D      P T C     I N P U T    ^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
    endif

    return

  END subroutine ptc_input
  !_________________________________________________________________

  SUBROUTINE SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123, skew_0123,ord_max)
    use twtrrfi ! integer, maxmul,maxferr,maxnaper
    implicit none
    ! 1) read multipole coeff. and errors for a current thick element
    ! 2) fill the error and multiploes arrays of data-bases
    REAL(dp), INTENT(IN) :: l
    type(keywords), INTENT(INOUT) ::  key
    REAL(dp), INTENT(OUT) :: normal_0123(0:3), skew_0123(0:3) ! n/l;
    REAL(dp) :: normal(0:maxmul), skew  (0:maxmul), &
         f_errors(0:maxferr), field(2,0:maxmul), bk0
    INTEGER :: n_norm, n_skew, n_ferr ! number of terms in command line
    INTEGER :: n_max
    INTEGER :: node_fd_errors ! function
    integer :: i, i_count, n_dim_mult_err, ord_max

    double precision, parameter :: zero=0.d0
    real(kind(1d0))   ::  node_value

    !initialization
    normal_0123(:)=zero
    skew_0123(:)=zero
    do i=1,NMAX
       key%list%k(i)=zero
       key%list%ks(i)=zero
    enddo

    ! real(dp) f_errors(0:maxferr),normal(0:maxmul),skew(0:maxmul)
    ! Get multipole components on bench !-----------------------!
    NORMAL = zero ! make zero "normal"                          !
    SKEW = zero   ! make zero "skew"                            !
    !                                                           !
    ! madxdict.h: "knl = [r, {0}], "                            !
    !             "ksl = [r, {0}], "                            !
    ! Assign values from the command line                       !
    call get_node_vector('knl ',n_norm,normal)                  !
    call get_node_vector('ksl ',n_skew,skew)                    !
    skew(0)=-skew(0)                                            ! frs error found 30.08.2008
    if(n_norm.ge.maxmul) n_norm=maxmul-1                        !
    if(n_skew.ge.maxmul) n_skew=maxmul-1                        !
    ord_max=max(n_norm,n_skew)                                  !
    ! void get_node_vector(char*par,int*length,double* vector)  !
    ! /* returns vector for parameter par of current element */ !
    !                                                           !
    ! get errors                                                !
    F_ERRORS = zero                                             !
    n_ferr = node_fd_errors(f_errors) !                         !
    ! /* returns the field errors of a node */                  !
    FIELD = zero ! array to be zeroed.                          !
    if (n_ferr .gt. 0) then                                     !
       call dcopy(f_errors,field,n_ferr)                        !
    endif                                                       !
    !-----------------------------------------------------------!

    ! fill strength of ALL normal multipoles
    if(n_norm.gt.0) then  ! ===============================!
       do i_count=0,n_norm                                 !
          if(i_count.gt.0) then                            !
             if(l.ne.zero) then                            !
                key%list%k(i_count+1)=normal(i_count)/l    !
             else                                          !
                key%list%k(i_count+1)=normal(i_count)      !
             endif                                         !
          endif                                            !
          if (i_count.le.3) then                           !
             if(l.ne.zero) then                            !
                normal_0123(i_count)=normal(i_count)/l     !
             else                                          !
                normal_0123(i_count)=normal(i_count)       !
             endif                                         !
          endif                                            !
       enddo                                               !
    endif !================================================!

    ! fill strength of ALL skew multipoles
    if(n_skew.gt.0) then  ! ===============================!
       do i_count=0,n_skew                                 !
          if(i_count.gt.0) then                            !
             if(l.ne.zero) then                            !
                key%list%ks(i_count+1)=skew(i_count)/l     !
             else                                          !
                key%list%ks(i_count+1)=skew(i_count)       !
             endif                                         !
          endif                                            !
          if (i_count.le.3) then                           !
             if(l.ne.zero) then                            !
                skew_0123(i_count)=skew(i_count)/l         !
             else                                          !
                skew_0123(i_count)=skew(i_count)           !
             endif                                         !
          endif                                            !
       enddo                                               !
    endif !================================================!

    n_dim_mult_err = max(n_norm, n_skew, n_ferr/2) !===========!
    if(n_dim_mult_err.ge.maxmul) n_dim_mult_err=maxmul-1       !

    n_max = -1                                 !
    if(n_ferr.gt.0) then
       if (getdebug() > 2) then
         print*,"Reading errors ", n_ferr
       endif
       do i_count=n_dim_mult_err,0,-1                             !

          if (getdebug() > 2) then
            print*,"   >> error n ",i_count+1,"  kn = ", field(1,i_count), " ks = ", field(2,i_count)
          endif                                  !

          IF( (field(1,i_count)/=0.0_dp .or. field(2,i_count)/=0.0_dp) .and. (n_max < 0)) THEN
            n_max = i_count+1
          ENDIF

          if(l.ne.zero) then                                   !
             key%list%k(i_count+1)=key%list%k(i_count+1)+ &    !
                  field(1,i_count)/l                           !
             key%list%ks(i_count+1)=key%list%ks(i_count+1)+ &  !
                  field(2,i_count)/l                           !
          else                                                 !
             key%list%k(i_count+1)=key%list%k(i_count+1)+ &    !
                  field(1,i_count)                             !
             key%list%ks(i_count+1)=key%list%ks(i_count+1)+ &  !
                  field(2,i_count)                             !
          endif                                                !
       enddo                                                   !
    endif !====================================================!
    if (key%magnet == 'sbend' .or. key%magnet == 'rbend') then 
      bk0 = node_value('k0 ')
      if(bk0 .ne. 0) key%list%k(1) = key%list%k(1) + bk0 - node_value('angle ')/l
    endif

    ord_max=max(ord_max,n_max)

  END SUBROUTINE SUMM_MULTIPOLES_AND_ERRORS
  !----------------------------------------------------------------
 !_________________________________________________________________

  SUBROUTINE REFRESH_MULTIPOLES(l, normal_0123, skew_0123,ord_max,normal,skew)
    use twtrrfi ! integer, maxmul,maxferr,maxnaper
    implicit none
    ! 1) read multipole coeff. and errors for a current thick element
    ! 2) fill the error and multiploes arrays of data-bases
    REAL(dp), INTENT(IN) :: l
    REAL(dp), INTENT(OUT) :: normal_0123(0:3), skew_0123(0:3) ! n/l;
    REAL(dp) :: normal(0:maxmul), skew  (0:maxmul), &
         f_errors(0:maxferr), field(2,0:maxmul)
    INTEGER :: n_norm, n_skew, n_ferr ! number of terms in command line
    INTEGER :: node_fd_errors ! function
    integer :: i, i_count, n_dim_mult_err, ord_max

    double precision, parameter :: zero=0.d0

    !initialization
    normal_0123(:)=zero
    skew_0123(:)=zero

    ! real(dp) f_errors(0:maxferr),normal(0:maxmul),skew(0:maxmul)
    ! Get multipole components on bench !-----------------------!
    NORMAL = zero ! make zero "normal"                          !
    SKEW = zero   ! make zero "skew"                            !
    !                                                           !
    ! madxdict.h: "knl = [r, {0}], "                            !
    !             "ksl = [r, {0}], "                            !
    ! Assign values from the command line                       !
    call get_node_vector('knl ',n_norm,normal)                  !
    call get_node_vector('ksl ',n_skew,skew)                    !
    skew(0)=-skew(0)                                            ! frs error found 30.08.2008
    if(n_norm.ge.maxmul) n_norm=maxmul-1                        !
    if(n_skew.ge.maxmul) n_skew=maxmul-1                        !
    ord_max=max(n_norm,n_skew)                                  !
    if(l.ne.zero) then                            !
       NORMAL = NORMAL / l
       SKEW = SKEW / l
       !do i=0,maxmul
       !   normal(i)=normal(i)/l
       !   skew(i)=skew(i)/l
       !enddo
    endif
    NORMAL_0123(0:3) = NORMAL(0:3)
    SKEW_0123(0:3) = SKEW(0:3)
    !do i=0,3
    !   normal_0123(i)=normal(i)
    !   skew_0123(i)=skew(i)
    !enddo

    ! get errors                                                !
    F_ERRORS = zero                                             !
    n_ferr = node_fd_errors(f_errors) !                         !
    ! /* returns the field errors of a node */                  !
    FIELD = zero ! array to be zeroed.                          !
    if (n_ferr .gt. 0) then                                     !
       call dcopy(f_errors,field,n_ferr)                        !
    endif                                                       !
    !-----------------------------------------------------------!

    n_dim_mult_err = max(n_norm, n_skew, n_ferr/2) !===========!
    if(n_dim_mult_err.ge.maxmul) n_dim_mult_err=maxmul-1       !
    if(n_ferr.gt.0) then                                       !
       do i_count=0,n_dim_mult_err                             !
          if(l.ne.zero) then                                   !
             normal(i_count+1)=normal(i_count+1)+field(1,i_count)/l
             skew(i_count+1)=skew(i_count+1)+field(2,i_count)/l
          else                                                 !
             normal(i_count+1)=normal(i_count+1)+field(1,i_count)
             skew(i_count+1)=skew(i_count+1)+field(2,i_count)
          endif                                                !
       enddo                                                   !
    endif !====================================================!



  END SUBROUTINE REFRESH_MULTIPOLES
  !----------------------------------------------------------------

  subroutine ptc_getnfieldcomp(fibreidx, ncomp, nval)
    implicit none
    real(kind(1d0))      :: nval
    integer              :: fibreidx
    integer              :: ncomp
    type(fibre), pointer :: p
    integer              :: j

    p=>my_ring%start
    do j=1, fibreidx
       p=>p%next
    enddo

    ncomp = ncomp + 1
    nval = p%mag%BN(ncomp)

  end subroutine  ptc_getnfieldcomp
  !----------------------------------------------------------------

  subroutine ptc_getsfieldcomp(fibreidx, ncomp, nval)
    implicit none
    real(kind(1d0))      :: nval
    integer              :: fibreidx
    integer              :: ncomp
    type(fibre), pointer :: p
    integer              :: j

    p=>my_ring%start
    do j=1, fibreidx
       p=>p%next
    enddo

    ncomp = ncomp + 1

    nval = p%mag%AN(ncomp)
    print*, "Returning AN",nval," for ",p%mag%name


  end subroutine  ptc_getsfieldcomp
  !----------------------------------------------------------------

  subroutine ptc_setfieldcomp(fibreidx)
    implicit none
    integer              :: fibreidx
    type(fibre), pointer :: p
    integer              :: j, i
    integer              :: kn, ks
    real(dp)             :: v
    real(kind(1d0)) get_value

    if ( .not. associated(my_ring) ) then
       call fort_warn("ptc_setfieldcomp","No active PTC layout/period")
       return
    endif

    if (getdebug()>2) then
       print*, "I am in ptc_setfieldcomp: Element index is ", fibreidx
    endif

    if ( (fibreidx .lt. 1) .and. (fibreidx .gt. my_ring%n) ) then
       call fort_warn("ptc_setfieldcomp","element out of range of the current layout")
       return
    endif

    p=>my_ring%start
    do j=1, fibreidx
       p=>p%next
    enddo

    if (getdebug() > 1 ) then
       print*,"Found element no. ", fibreidx," named ", p%mag%name, &
            &" of kind ", p%mag%kind, mytype(p%mag%kind)
       print*,"Currently nmul is ", p%mag%p%nmul

       write(6,*) "BNs",p%mag%BN
       write(6,*) "ANs",p%mag%AN

       DO i=1,p%mag%p%nmul
          print*, "Polimorphic BN(",i,")"
          call print(p%mag%BN(i),6)
          print*, "Polimorphic AN(",i,")"
          call print(p%mag%AN(i),6)
       ENDDO

    endif

    kn = get_value('ptc_setfieldcomp ','kn ')
    v = get_value('ptc_setfieldcomp ','value ')

    if (kn >= 0) then
       kn = kn + 1

       if (getdebug() > 1) then
           print*,"Setting up KN ", kn, " from ", p%mag%BN(kn) ," to ", v
       endif

       call add(p%mag, kn,0,v)
       call add(p%magp,kn,0,v)


    else
       ks = get_value('ptc_setfieldcomp ','ks ')
       if (ks < 0) then
          call fort_warn("ptc_setfieldcomp","neither kn nor ks specified")
          return
       endif
       ks = ks + 1

       !      print*,"Setting up skew field component ", ks," to ", v

       if (getdebug() > 1) then
           print*,"Setting up KS ", ks, " from ", p%mag%AN(ks) ," to ", v
       endif
       call add(p%mag, -ks,0,v)
       call add(p%magp,-ks,0,v)

    endif

    if (getdebug() > 1 ) then
       write(6,*) "BNs",p%mag%BN
       write(6,*) "ANs",p%mag%AN
       write(6,*) ""
    endif
  end subroutine ptc_setfieldcomp
  !----------------------------------------------------------------

  subroutine ptc_align()
    use twiss0fi
    implicit none
    integer j,n_align,node_al_errors
    integer restart_sequ,advance_node
    real(dp) al_errors(align_max)
    type(fibre), pointer :: f
    !---------------------------------------------------------------


    j=restart_sequ()
    j=0
    f=>my_ring%start
10  continue

    j=j+1
    n_align = node_al_errors(al_errors)
    if (n_align.ne.0)  then
      if (getdebug() > 3) then
        write(6,*) " ----------------------------------------------- "
        write(6,*) f%mag%name," Translation Error "
        write(6,'(3f11.8)') al_errors(1:3)
        write(6,*) f%mag%name," Rotation Error "
        write(6,'(3f11.8)') al_errors(4:6)
        call print_elframes(f)
      endif

      ! this routine is buggy -> fixed by Etienne on 2018.01.29
      ! call mad_misalign_fibre(f,al_errors(1:6))

      ! this is PTC original, but it first rotates and than shifts (in frame after rotations)
      ! call misalign_fibre(f,al_errors(1:6))

      !our new routine
      call misalign_element(f,al_errors)

      !a workaround to handle misaligment of thin dipoles that is not handled by PTC
      call misalign_thindipole(f,al_errors)


      if (getdebug() > 3) then

        write(6,*) " vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv "
        call print_elframes(f)
      endif

    endif
    f=>f%next
    if(advance_node().ne.0)  goto 10


  END subroutine ptc_align
  !_________________________________________________________________

  subroutine print_elframes(f)
    implicit none
    TYPE(FIBRE),target,INTENT(INOUT):: f

        write(6,*) "Ac:", f%chart%f%a
        write(6,*) "Oc:", f%chart%f%o
        write(6,*) "Bc:", f%chart%f%b
        write(6,*)
        write(6,*) "Am:", f%mag%p%f%a
        write(6,*) "Om:", f%mag%p%f%o
        write(6,*) "Bm:", f%mag%p%f%b
        write(6,*)
        write(6,*) "entm(1,:) :", f%mag%p%f%ent(1,:)
        write(6,*) "entm(2,:) :", f%mag%p%f%ent(2,:)
        write(6,*) "entm(3,:) :", f%mag%p%f%ent(3,:)
        write(6,*)
        write(6,*) "midm(1,:) :", f%mag%p%f%mid(1,:)
        write(6,*) "midm(2,:) :", f%mag%p%f%mid(2,:)
        write(6,*) "midm(3,:) :", f%mag%p%f%mid(3,:)
        write(6,*)
        write(6,*) "exim(1,:) :", f%mag%p%f%exi(1,:)
        write(6,*) "exim(2,:) :", f%mag%p%f%exi(2,:)
        write(6,*) "exim(3,:) :", f%mag%p%f%exi(3,:)
        write(6,*)
        write(6,*) "ang_in    :", f%CHART%ang_in
        write(6,*) "ang_out   :", f%CHART%ang_out

  end subroutine print_elframes

  !_____________________________________________________
  ! Routine to handle thin dipole misalignments
  ! Etienne says that KIND3 is "highly illigal" and he does not support it
  ! therefore we need to tweak ut ourselves
  ! Patching from rotation errors is done using FIBRE%CHART%ANG_OUT and FIBRE%CHART%ANG_OUT
  ! Here we calculate these angles
  subroutine misalign_thindipole(f,al_errors)
    use twiss0fi, only: align_max
    TYPE(FIBRE),target,INTENT(INOUT):: f
    REAL(DP),INTENT(IN) :: al_errors(align_max)
    REAL(DP) :: Fent(3,3), F0ent(3,3), Fexi(3,3), F0exi(3,3), A
    REAL(DP) :: F1(3,3), F2(3,3), F3(3,3)

    if (f%mag%kind /= kind3 ) return

    !in MADX vertical bend is tilted horizontal bend, so thin_v_angle is always zero
    if ( abs(f%mag%K3%thin_h_angle) .lt. 1e-12 ) return;

    if (getdebug() > 2) then
      write(6,*) "misalign_thindipole angle = ",f%mag%K3%thin_h_angle
    endif


    A = f%mag%K3%thin_h_angle
    ! identity
    F0ent = zero
    F0ent(1,1) =  1
    F0ent(2,2) =  1
    F0ent(3,3) =  1

    !!!!!!!!!!!!!
    ! depends only on the bend angle
    F0exi = zero
    F0exi(1,1) =  cos(A)
    F0exi(1,3) =  sin(A)
    F0exi(2,2) =  1
    F0exi(3,1) = -sin(A)
    F0exi(3,3) =  cos(A)


    !!!!!!!!!!!!!
    ! Lumped error matrix
    f1 = zero
    f2 = zero
    f3 = zero

    f1(1,1) = 1
    f1(2,2) = cos(al_errors(4))
    f1(2,3) =-sin(al_errors(4))
    f1(3,2) = sin(al_errors(4))
    f1(3,3) = cos(al_errors(4))

    f2(1,1) = cos(al_errors(5))
    f2(1,3) =-sin(al_errors(5))
    f2(2,2) = 1
    f2(3,1) = sin(al_errors(5))
    f2(3,3) = cos(al_errors(5))

    ! change of sign
    f3(1,1) = cos(-al_errors(6))
    f3(1,2) =-sin(-al_errors(6))
    f3(2,1) = sin(-al_errors(6))
    f3(2,2) = cos(-al_errors(6))
    f3(3,3) = 1

    Fent = matmul(f2,  f1)
    Fent = matmul(f3,Fent)

    !____________________________________

    !lumped bend rot with alignment rotation
    Fexi = matmul(F0exi,Fent)

    ! Calculate athe angles needed for tracking (a PTC routine)
    CALL COMPUTE_ENTRANCE_ANGLE(F0ENT,FENT,f%CHART%ANG_IN)
    CALL COMPUTE_ENTRANCE_ANGLE(FEXI,F0EXI,f%CHART%ANG_OUT)

  end subroutine misalign_thindipole

  !_____________________________________________________
  subroutine misalign_element(f,al_errors)
    use twiss0fi, only: align_max
    TYPE(FIBRE),target,INTENT(INOUT):: f
    REAL(DP),INTENT(IN) :: al_errors(align_max)
    REAL(DP)             :: mis(align_max), omegat(3), basist(3,3)

      mis=0
      mis(4:6)=al_errors(4:6)
      mis(4:5)=-mis(4:5)

      OMEGAT=f%mag%p%f%a
      basist=f%mag%p%f%ent
      CALL MISALIGN_FIBRE(f,mis,OMEGAT,BASIST,ADD=.false.) ! false to remove previous alignment

      mis=0
      mis(1:3)=al_errors(1:3)
      CALL MISALIGN_FIBRE(f,mis,OMEGAT,BASIST,ADD=.true.)

  end subroutine misalign_element
  !_________________________________________________________________

  subroutine ptc_dumpmaps()
    !Dumps to file maps and/or matrixes (i.e. first order maps)
    implicit none
    type(fibre), pointer :: p
    type(damap)          :: id !identity map used for calculating maps for each element
    type(real_8)         :: y2(6)  !polimorphes array used for calculating maps for each element
    type(real_8)         :: yfull(6)  !polimorphes array used for calculating maps for each element
    real(dp)             :: xt(6)
    integer              :: i !iterators
    integer mf1,mf2
    character(200)       :: filename='ptcmaps.txt'
    character(200)       :: filenamefull='ptcmaps'
    integer              :: flag_index,why(9)
    real(kind(1d0))      :: suml=zero
    integer  geterrorflag !C function that returns errorflag value

    suml=zero

    if (cavsareset .eqv. .false.) then
       call setcavities(my_ring,maxaccel)
       if (geterrorflag() /= 0) then
          return
       endif
    endif

    if (getdebug() > 1) then
        print *, '<madx_ptc_module.f90 : ptc_dumpmaps> Maps are dumped to file ',filename
    endif
    call kanalnummer(mf1)
    open(unit=mf1,file=filename)

    !    write(filenamefull,*) filename,".",my_ring%start%mag%name,"-",my_ring%end%mag%name,".txt"
    filenamefull="ptcmaps.start-end.txt"
    print*, filenamefull
    call kanalnummer(mf2)
    open(unit=mf2,file=filenamefull)

    print*, "no=1"," mynd2=",c_%nd2," npara=",c_%npara
    call init(getintstate(),1,c_%np_pol,berz)

    call alloc(id);
    call alloc(y2);
    call alloc(yfull);

    xt(:) = zero
    id    = 1     ! making identity map

    yfull  = xt + id

    p=>my_ring%start
    do i=1,my_ring%n


       y2=xt+id ! we track identity map from the current position

       if( (p%mag%kind/=kind21) .and. (p%mag%kind/=kind4) ) then

          call track(my_ring,y2,i,i+1,getintstate())

          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps ",whymsg);
             close(mf1)
             close(mf2)
             return
          endif

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)

             Write(6,*) "ptc_dumpmaps: APERTURE error for element: ",i," name: ",p%MAG%name
             write(whymsg,*) 'APERTURE error: ',why
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps: ",whymsg);
             c_%watch_user=.false.
             close(mf1)
             close(mf2)
             return
          endif

          call track(my_ring,xt,i,i+1,getintstate())
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps ",whymsg);
             close(mf1)
             close(mf2)
             return
          endif

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) "ptc_dumpmaps: APERTURE error for element: ",i," name: ",p%MAG%name
             write(whymsg,*) 'APERTURE error: ',why
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps: ",whymsg);
             c_%watch_user=.false.
             close(mf1)
             close(mf2)
             return
          endif
       else
          if (getdebug() > 2) then
              print *, 'Track Cavity...'
          endif

          call track(my_ring,y2,i,i+2,getintstate())
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps ",whymsg);
             close(mf1)
             close(mf2)
             return
          endif

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) "ptc_dumpmaps: APERTURE error for element: ",i," name: ",p%MAG%name
             write(whymsg,*) 'APERTURE error: ',why
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps: ",whymsg);
             c_%watch_user=.false.
             close(mf1)
             close(mf2)
             return
          endif

          call track(my_ring,xt,i,i+2,getintstate())
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps ",whymsg);
             close(mf1)
             close(mf2)
             return
          endif

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) "ptc_dumpmaps: APERTURE error for element: ",i," name: ",p%MAG%name
             write(whymsg,*) 'APERTURE error: ',why
             call fort_warn('ptc_dumpmaps: ',whymsg)
             call seterrorflag(10,"ptc_dumpmaps: ",whymsg);
             c_%watch_user=.false.
             close(mf1)
             close(mf2)
             return
          endif
       endif


       call track(my_ring,yfull,i,i+1,getintstate())

       if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
          write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost
          call fort_warn('ptc_dumpmaps: ',whymsg)
          call seterrorflag(10,"ptc_dumpmaps ",whymsg);
          close(mf1)
          close(mf2)
          return
       endif

       call PRODUCE_APERTURE_FLAG(flag_index)
       if(flag_index/=0) then
          call ANALYSE_APERTURE_FLAG(flag_index,why)

          Write(6,*) "ptc_dumpmaps: APERTURE error for element: ",i," name: ",p%MAG%name
          write(whymsg,*) 'APERTURE error: ',why
          call fort_warn('ptc_dumpmaps: ',whymsg)
          call seterrorflag(10,"ptc_dumpmaps: ",whymsg);
          c_%watch_user=.false.
          close(mf1)
          close(mf2)
          return
       endif

       write(mf2,*) p%mag%name, suml,' m ==========================='
       call print(yfull,mf2)

       suml=suml+p%MAG%P%ld

       write(mf1,*) p%mag%name, suml,' m ==========================='
       if (c_%npara == 6) then
          call dump6dmap(y2, mf1)
       elseif (c_%npara == 5) then
          call dump5dmap(y2, mf1)
       elseif (c_%npara == 4) then
          call dump4dmap(y2, mf1)
       else
          call fort_warn("ptc_dumpmaps","c_%npara is neither 6,5 nor 4")
       endif
       p=>p%next
    enddo

    close(mf1)
    call kill(y2);
    call kill(id);

    !_________________________________________________________________
    !_________________________________________________________________
    !_________________________________________________________________

  contains
    !_________________________________________________________________
    subroutine dump4dmap(y2, fun)
      implicit none
      double precision a1000,a0100,a0010,a0001
      type(real_8) :: y2(6)  !polimorphes array used for calculating maps for each element
      integer      :: fun !file unit number
      integer      :: ii

      if (getdebug() > 1) then

      endif

      do ii=1,4
         a1000=y2(ii).sub.'1000'
         a0100=y2(ii).sub.'0100'
         a0010=y2(ii).sub.'0010'
         a0001=y2(ii).sub.'0001'
         write(fun,'(6f13.8)')  a1000, &
              &                 a0100, &
              &                 a0010, &
              &                 a0001
      enddo

    end subroutine dump4dmap
    !_________________________________________________________________

    subroutine dump5dmap(y2, fun)
      implicit none
      double precision a10000,a01000,a00100,a00010,a00001
      type(real_8) :: y2(6)  !polimorphes array used for calculating maps for each element
      integer      :: fun !file unit number
      integer      :: ii
      do ii=1,5
         a10000=y2(ii).sub.'10000'
         a01000=y2(ii).sub.'01000'
         a00100=y2(ii).sub.'00100'
         a00010=y2(ii).sub.'00010'
         a00001=y2(ii).sub.'00001'
         write(fun,'(6f13.8)')  a10000, &
              &                 a01000, &
              &                 a00100, &
              &                 a00010, &
              &                 a00001     !
      enddo

    end subroutine dump5dmap
    !_________________________________________________________________

    subroutine dump6dmap(y2, fun)
      implicit none
      double precision a100000,a010000,a001000,a000100,a000010,a000001
      type(real_8) :: y2(6)  !polimorphes array used for calculating maps for each element
      integer      :: fun !file unit number
      integer      :: ii

      do ii=1,4
         a100000=y2(ii).sub.'100000'
         a010000=y2(ii).sub.'010000'
         a001000=y2(ii).sub.'001000'
         a000100=y2(ii).sub.'000100'
         a000010=y2(ii).sub.'000010'
         a000001=y2(ii).sub.'000001'
         write(fun,'(6f13.8)')  a100000, &
              &                 a010000, &
              &                 a001000, &
              &                 a000100, &
              &                 a000001, & !madx format has dp/p at the last column
              &                 a000010    !
      enddo

      do ii=6,5,-1
         a100000=y2(ii).sub.'100000'
         a010000=y2(ii).sub.'010000'
         a001000=y2(ii).sub.'001000'
         a000100=y2(ii).sub.'000100'
         a000010=y2(ii).sub.'000010'
         a000001=y2(ii).sub.'000001'
         write(fun,'(6f13.8)')  a100000, &
              &                 a010000, &
              &                 a001000, &
              &                 a000100, &
              &                 a000001, & !madx format has dp/p at the last column
              &                 a000010    !
      enddo
    end subroutine dump6dmap


  end subroutine ptc_dumpmaps
  !________________________________________________________________________________________________

  RECURSIVE FUNCTION FACTORIAL (N) &
       RESULT (FACTORIAL_RESULT)
    INTEGER :: N, FACTORIAL_RESULT

    IF (N <= 0 ) THEN
       FACTORIAL_RESULT = 1
    ELSE
       FACTORIAL_RESULT = N * FACTORIAL (N-1)
    END IF
  END FUNCTION FACTORIAL
  !________________________________________________________________________________________________
  ! calculates product of factorials
  RECURSIVE FUNCTION FACTORIAL_PRODUCT (A,N) &
       RESULT (FACTORIAL_RESULT)
    INTEGER :: N, FACTORIAL_RESULT
    INTEGER :: A(6)
    integer :: i

    FACTORIAL_RESULT = 1

    do i=1,N
      FACTORIAL_RESULT = FACTORIAL_RESULT * FACTORIAL(A(i))
    enddo

  END FUNCTION FACTORIAL_PRODUCT
  !________________________________________________________________________________________________

  subroutine ptc_track()
    implicit none
    integer i,nint,ndble,nchar,int_arr(1),char_l,icase,turns,flag_index,why(9)
    integer j,next_start
    real(dp) x0(6),x(6),deltap0,deltap,dt
    real(dp)  xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx
    real(kind(1d0)) get_value
    logical(lp) closed_orbit
    character*12 char_a
    data char_a / ' ' /
    !------------------------------------------------------------------------------


    if(universe.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_track: ',' no universe created')
       return
    endif
    if(index_mad.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_track: ',' no layout created')
       return
    endif

    icase = get_value('ptc_track ','icase ')
    deltap0 = get_value('ptc_track ','deltap ')

    deltap = zero
    call my_state(icase,deltap,deltap0)

    if (getdebug() > 2) then
       print *, "ptc_track: internal state is:"
       call print(default,6)
    endif

    x0(:)=zero
    if(mytime) then
       call Convert_dp_to_dt (deltap, dt)
    else
       dt=deltap
    endif
    if(icase.eq.5) x0(5)=dt
    closed_orbit = get_value('ptc_track ','closed_orbit ') .ne. 0
    if(closed_orbit) then
       call find_orbit(my_ring,x0,1,default,c_1d_7)
       CALL write_closed_orbit(icase,x0)
    endif


    call comm_para('coord ',nint,ndble,nchar,int_arr,x,char_a,char_l)

    j  =  next_start(xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx)
    print*,"dat1",j,xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx
    j  =  next_start(xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx)
    print*,"dat2",j,xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx
    j  =  next_start(xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx)
    print*,"dat3",j,xx,pxx,yx,pyx,tx,deltaex,fxx,phixx,fyx,phiyx,ftx,phitx

    x(:)=x(:)+x0(:)
    print*,"  Initial Coordinates: ", x
    turns = get_value('ptc_track ','turns ')
    c_%watch_user=.true.
    do i=1,turns
       call track(my_ring,x,1,default)
       if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
          write(whymsg,*) 'DA got unstable: PTC msg: ',messagelost
          call fort_warn('ptc_track: ',whymsg)
          call seterrorflag(10,"ptc_track ",whymsg);
          return
       endif
       call PRODUCE_APERTURE_FLAG(flag_index)
       if(flag_index/=0) then
          call ANALYSE_APERTURE_FLAG(flag_index,why)
          Write(6,*) "ptc_track unstable (tracking)-programs continues "
          Write(6,*) why ! See produce aperture flag routine in sd_frame
          goto 100
       endif
    enddo
    c_%watch_user=.false.
    print*,"  End Coordinates: ",x
    return
100 continue
    c_%watch_user=.false.
    print*,"  Last Coordinates: ",x," after: ",i," turn(s)"

  END subroutine ptc_track
  !________________________________________________________________________________


  subroutine ptc_end()
    implicit none
    integer i

    if(universe.le.0.or.EXCEPTION.ne.0) then
       call fort_warn('return from ptc_end: ',' no universe can be killed')
       return
    endif

    call killsavedmaps() !module ptc_twiss -> kill buffered maps

    !    call killparresult()
    call resetknobs()  !remove the knobs

    call kill_map_cp()

    if ( associated(m_u%n) .eqv. .false. ) then
       print*, "We attempt to kill not initialized universe!"
    endif


    call kill_universe(m_u)
    nullify(my_ring)
    call kill_tpsa
!    do i=1,size(s_b)
!       call nul_coef(s_b(i))
!    enddo
!    deallocate(s_b)


    firsttime_coef=.true.

    universe=universe-1
  end subroutine ptc_end


  subroutine normalform_normalform(s1,s2)
    implicit none
    type (normalform),intent(inout)::s1
    type (normalform),intent(in)::s2
    integer i,j

    s1%a_t=s2%a_t
    s1%a1=s2%a1
    s1%a%constant(:)=s2%a%constant(:)
    s1%a%Linear=s2%a%Linear
    s1%a%nonlinear=s2%a%nonlinear
    s1%a%pb=s2%a%pb
    s1%normal%constant(:)=s2%normal%constant(:)
    s1%normal%Linear=s2%normal%Linear
    s1%normal%nonlinear=s2%normal%nonlinear
    s1%normal%pb=s2%normal%pb
    s1%DHDJ=s2%DHDJ
    do i=1,ndim
       s1%TUNE(i)=s2%TUNE(i)
       s1%damping(i)=s2%damping(i)
       s1%plane(i)=s2%plane(i)
       do j=1,mynreso
          s1%m(i,j)=s2%m(i,j)
       enddo
    enddo
    s1%nord=s2%nord
    s1%jtune=s2%jtune
    s1%nres=s2%nres
    s1%AUTO=s2%AUTO
  end subroutine normalform_normalform
  !_________________________________________________________________


  SUBROUTINE set_PARAMETERS(R,nt,iorder,IFAM,inda,scale)
    !Strength of Multipole of order iorder as parameter

    IMPLICIT NONE
    integer ipause, mypause
    logical(lp) ok
    INTEGER iorder,i,j,jj,k,lstr,IFAM,tot,nt,inda,min1
    INTEGER,parameter::ipara=100
    real(dp) scale(ipara),value
    character(20) str
    CHARACTER(3) STR1
    character(10),dimension(10)::multname
    type(layout) r
    type(fibre), POINTER :: current
    INTEGER,ALLOCATABLE,dimension(:)::DAFAM
    INTEGER,ALLOCATABLE,dimension(:,:)::FAM
    real(dp),ALLOCATABLE,dimension(:)::SFAM
    multname=(/"Dipole    ","Quadrupole","Sextupole ","Octupole  ","Decapole  ",&
         "Dodecapole","14-Pole   ","16-Pole   ","18-Pole   ","20-Pole   "/)

    ALLOCATE(FAM(IFAM,0:R%N),DAFAM(IFAM),SFAM(IFAM))

    min1=0
    if(iorder.lt.0) then
       min1=1
       iorder=-iorder
    endif
    DO I=1,IFAM
       OK=.TRUE.
       DO WHILE(OK)
          TOT=0
          if(min1.eq.0) WRITE(6,*) " Identify ",multname(iorder)
          if(min1.eq.1) WRITE(6,*) " Identify ","SKEW-"//multname(iorder)
          READ(5,*) STR
          STR=TRIM(ADJUSTL(STR))
          LSTR=LEN_TRIM (STR)
          current=>r%start
          DO J=1,R%N
             IF(current%MAG%NAME==STR.and.current%MAG%P%NMUL==iorder) THEN
                TOT=TOT+1
                FAM(I,TOT)=J
             ENDIF
             current=>current%next
          ENDDO
          WRITE(6,*) TOT," Is that OK? YES or NO?"
          READ (5,*) STR1
          STR1=TRIM(ADJUSTL(STR1))
          IF(STR1(1:1)=='Y'.OR.STR1(1:1)=='y') THEN
             OK=.FALSE.
             inda=inda+1
             if(inda.gt.100) then
                write(6,*) " Problem: Only ",ipara," Parameters allowed"
                ipause=mypause(2002)
             endif
             DAFAM(I)=inda
             WRITE(6,*) " Give Scaling Factor, '0' uses Default"
             read(5,*) value
             if(value==0) then
                WRITE(6,*) " Take Default Scaling Value : ",scale(inda)
                SFAM(I)=scale(inda)
             else
                SFAM(I)=value
             endif
          ENDIF
       ENDDO

       FAM(I,0)=TOT
       current=r%start
       DO JJ=1,FAM(I,0)
          J=FAM(I,JJ)
          ! ALLOCATION GYMNASTIC IF Multipole NOT YET ALLOCATED
          IF(current%MAGP%P%NMUL<iorder) THEN
             CALL KILL(current%MAGP%BN,current%MAGP%P%NMUL)
             CALL KILL(current%MAGP%AN,current%MAGP%P%NMUL)
             current%MAGP%P%NMUL=iorder
             DEALLOCATE(current%MAGP%BN)
             DEALLOCATE(current%MAGP%AN)
             CALL ALLOC(current%MAGP%BN,iorder)
             CALL ALLOC(current%MAGP%AN,iorder)
             ALLOCATE(current%MAGP%BN(iorder),current%MAGP%AN(iorder))
             DO K=1,current%MAG%P%NMUL
                current%MAGP%BN(K)=current%MAG%BN(K)
                current%MAGP%AN(K)=current%MAG%AN(K)
             ENDDO
             DEALLOCATE(current%MAG%BN)
             DEALLOCATE(current%MAG%AN)
             ALLOCATE(current%MAG%BN(iorder),current%MAG%AN(iorder))
             call equal(current%MAG,current%MAGP)
          ENDIF
          if(min1.eq.0) then
             current%MAGP%BN(iorder)%I=NT+I
             current%MAGP%BN(iorder)%KIND=3
          else
             current%MAGP%AN(iorder)%I=NT+I
             current%MAGP%AN(iorder)%KIND=3
          endif
          current=>current%next
       ENDDO
    ENDDO

    current=r%start
    DO I=1,IFAM
       DO JJ=1,1
          J=FAM(I,JJ)
          if(min1.eq.0) WRITE(6,*)  current%MAG%NAME,' ', current%MAG%BN(iorder)
          if(min1.eq.1) WRITE(6,*)  current%MAG%NAME,' ', current%MAG%AN(iorder)
          current=>current%next
       ENDDO
    ENDDO

    DEALLOCATE(FAM,STAT=I)
    !    WRITE(6,*) I
    DEALLOCATE(DAFAM,STAT=I)
    !    WRITE(6,*) I
    DEALLOCATE(SFAM,STAT=I)
    !    WRITE(6,*) I

  end subroutine set_PARAMETERS
  !______________________________________________________________________

  subroutine my_state(icase,deltap,deltap0, silent)
    implicit none
    integer icase,i
    real(dp) deltap0,deltap
    logical, optional :: silent
    logical :: verbose

    ! force no printout, ugly work around of ugly ptc_track
    if (present(silent)) then
      verbose = .not. silent
    else
      verbose = .true.
    endif


    default = getintstate()

    if (getdebug()>1 .and. verbose) then
       print*, "icase=",icase," deltap=",deltap," deltap0=",deltap0
    endif

    if (getdebug()>3 .and. verbose) then
       print*, "Input State"
       call print(default,6)
    endif

    deltap = zero
    select case(icase)
    CASE(4)
       if (getdebug()>1 .and. verbose) then
           print*, "my_state: Enforcing ONLY_4D+NOCAVITY and NO DELTA"
       endif
       default = default - delta0 + only_4d0 + NOCAVITY0
       i=4
    CASE(5)
       if (getdebug()>1 .and. verbose) then
           print*, "my_state: Enforcing DELTA"
       endif
       default = default + delta0
       deltap = deltap0
       i=5
    CASE(56)
       if (getdebug()>1 .and. verbose) then
           print*, "my_state: Enforcing coasting beam"
       endif
       default = default - delta0 - only_4d0 + NOCAVITY0
       deltap = deltap0
       i=56
    CASE(6)
       i=6
    CASE DEFAULT
       default = default + only_4d0 + NOCAVITY0
       i=4
    END SELECT

    if (i==6) then
       if (getdebug()>2 .and. verbose) then
         print*,"icav=",icav," my_ring%closed=",my_ring%closed," getenforce6D()=",getenforce6D()
       endif

       if ( (icav==0) .and. my_ring%closed .and. (getenforce6D() .eqv. .false.)) then
          default = default - delta0 - only_4d0 + NOCAVITY0
          call fort_warn('my_state: ',' no cavity - dimensionality reduced 6 -> 5 and 1/2')
          i=56
       else
          default = default - delta0 - only_4d0 - NOCAVITY0 !enforcing nocavity to false
       endif

    endif

    if (i==6 .and. (default%time .eqv. .false. )) then
      call fort_warn("my_state ",  &
            "TIME=false with RF cavities gives approximative results valid only for fully relativistic beams (beta~1)")
    endif


    call setintstate(default)
    CALL UPDATE_STATES

    if (getdebug()>0 .and. verbose) then
      !print*, "Resulting state"
      call print(default,6)
    endif

    !Update the global flag of ptc_module
    mytime = default%time

    icase = i

  end subroutine my_state
  !________________________________________________________________________________________________

  subroutine f90flush(i,option)
    implicit none
    integer i,ios
    logical(lp) ostat, fexist,option
    logical fexist1, ostat1
    character*20 faction,faccess,fform,fwstat,fposition
    character*255 fname
    inquire(err=1,iostat=ios,&
         unit=i,opened=ostat1,exist=fexist1,write=fwstat)
    fexist = fexist1
    ostat  = ostat1
    if (.not.ostat.or..not.fexist.or.fwstat.ne.'YES') return
    inquire(err=2,iostat=ios,&
         unit=i,action=faction,access=faccess,&
         form=fform,name=fname,position=fposition)
    close (unit=i,err=3)
    !     write (*,*) 'Re-opening ',i,' ',option,ios,faction,faccess,fform,fposition,fname
    if (option) then
       open(err=4,iostat=ios,&
            unit=i,action=faction,access=faccess,form=fform,&
            file=fname,status='old',position='append')
    else
       open(err=4,iostat=ios,&
            unit=i,action=faction,access=faccess,form=fform,&
            file=fname,status='replace',position='rewind')
    endif
    return
1   write (*,*)&
         ' F90FLUSH 1st INQUIRE FAILED with IOSTAT ',ios,' on UNIT ',i
    stop
2   write (*,*)&
         ' F90FLUSH 2nd INQUIRE FAILED with IOSTAT ', ios,' on UNIT ',i
    stop
3   write (*,*)&
         ' F90FLUSH CLOSE FAILED with IOSTAT ',ios,' on UNIT ',i
    stop
4   write (*,*)&
         ' F90FLUSH RE-OPEN FAILED with IOSTAT ',ios,' on UNIT ',i
    stop
  end subroutine f90flush
  !________________________________________________________________________________________________

  SUBROUTINE write_closed_orbit(icase,x)
    implicit none
    INTEGER,  INTENT(IN):: icase
    REAL (dp),INTENT(IN) :: x(6)
    if(icase.eq.4) then
       print*,"Closed orbit: ",x(1),x(2),x(3),x(4)
    elseif(icase.eq.5) then
       print*,"Closed orbit: ",x(1),x(2),x(3),x(4),x(5)
    elseif(icase.eq.6) then
       print*,"Closed orbit: ",x(1),x(2),x(3),x(4),-x(6),x(5)
    endif
  ENDSUBROUTINE write_closed_orbit
  !________________________________________________________________________________________________

  SUBROUTINE Convert_dp_to_dt(deltap, dt)
    implicit none
    ! convert deltap=(p-p0)/p0 to dt=deltaE/p0c
    REAL(dp), INTENT(IN)  :: deltap
    REAL(dp), INTENT(OUT) :: dt

    ! local
    real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet

    ! to get "energy" value
    Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)

    IF (beta0.gt.zero ) THEN
       dt=SQRT(deltap*(deltap+two)+one/beta0/beta0)-one/beta0
    ELSE  ! exculde devision by 0
       call aafail('SUBR. Convert_dp_to_dt: ',' CALL GET_ONE => beta0.LE.0')
    ENDIF

  END SUBROUTINE Convert_dp_to_dt
  !=============================================================================
  SUBROUTINE Convert_dt_to_dp(dt, deltap )
    implicit none
    ! convert deltap=(p-p0)/p0 to dt=deltaE/p0c
    REAL(dp), INTENT(IN)  :: dt
    REAL(dp), INTENT(OUT) :: deltap

    ! local
    real(dp) :: MASS_GeV, ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet

    ! to get "energy" value
    Call GET_ONE(MASS_GeV,ENERGY,KINETIC,BRHO,BETA0,P0C,gamma0I,gambet)

    IF (beta0.gt.zero ) THEN
       deltap=SQRT( one + two*dt/beta0 +dt*dt)-one
    ELSE  ! exculde devision by 0
       call aafail('SUBR. Convert_dp_to_dt: ',' CALL GET_ONE => beta0.LE.0')
    ENDIF

  END SUBROUTINE Convert_dt_to_dp
  !=============================================================================

  subroutine makemaptable(y,order)
    implicit none
    integer           :: order
    type(real_8):: y(6)
    integer,parameter :: i_map_coor=10
    integer           :: map_term, ja(6),i,ii,iii
    integer           :: i1,i2,i3,i4,i5,i6
    integer           :: no
    real(dp)          :: coef
    real(kind(1d0))   :: map_coor(i_map_coor)
    real(kind(1d0))   :: get_value
    character(len=10)  :: coeffCode="c0_000000"//CHAR(0) ! polynomial no + "_" + monomial code + space(end of line)


    !    type(universal_taylor) :: ut

    !    write(0,*) "MAP_TABLE"

    map_term=42
    call  make_map_table(map_term)


    call liepeek(iia,icoast)
    allocate(j(c_%npara))
    ja(:)    = 0
    j(:)     = 0

    ! note that the order in which the coefficients appear in the map_table slightly
    ! differ from the order in which they appear in fort.18
100 do i=1,c_%npara ! distribute exponents over 6 variables, knowing their sum
       do no=0,order
          if (c_%npara.eq.6) then
             do i1=no,0,-1
                do i2=no-i1,0,-1
                   do i3=no-i1-i2,0,-1
                      do i4=no-i1-i2-i3,0,-1
                         do i5=no-i1-i2-i3-i4,0,-1
                            do i6=no-i1-i2-i3-i4-i5,0,-1
                               if (i1+i2+i3+i4+i5+i6==no) then
                                  !write(0,'(6(i4))'), i1,i2,i3,i4,i5,i6
                                  j(1)=i1
                                  j(2)=i2
                                  j(3)=i3
                                  j(4)=i4
                                  j(5)=i5
                                  j(6)=i6
                                  coef = y(i)%T.sub.j
                                  if (coef.ne.zero) then
                                     map_coor(1)=coef
                                     map_coor(2)=i
                                     map_coor(3)=c_%npara
                                     map_coor(4)=no
                                     map_coor(5)=j(1)
                                     map_coor(6)=j(2)
                                     map_coor(7)=j(3)
                                     map_coor(8)=j(4)
                                     map_coor(9)=j(5)
                                     map_coor(10)=j(6)
                                     call vector_to_table_curr("map_table ", 'coef ', map_coor(1), i_map_coor)
                                     write(coeffCode,'(i1,a1,6i1)') i,'_',j(1:c_%npara)
                                     coeffCode = 'c'//coeffCode(1:8)//CHAR(0)
                                     call string_to_table_curr("map_table ", "name ", coeffCode);
                                     call augment_count("map_table ")
                                  endif
                                  !write(0,*) 'write coef', coef
                               endif
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif (c_%npara.eq.5) then ! distribute exponents over 5 variables, knowing their sum
             do i1=no,0,-1
                do i2=no-i1,0,-1
                   do i3=no-i1-i2,0,-1
                      do i4=no-i1-i2-i3,0,-1
                         do i5=no-i1-i2-i3-i4,0,-1
                            if (i1+i2+i3+i4+i5==no) then
                               j(1)=i1
                               j(2)=i2
                               j(3)=i3
                               j(4)=i4
                               j(5)=i5
                               coef = y(i)%T.sub.j
                               if (coef.ne.zero) then
                                 map_coor(1)=coef
                                 map_coor(2)=i
                                 map_coor(3)=c_%npara
                                 map_coor(4)=no
                                 map_coor(5)=j(1)
                                 map_coor(6)=j(2)
                                 map_coor(7)=j(3)
                                 map_coor(8)=j(4)
                                 map_coor(9)=j(5)
                                 map_coor(10) = 0
                                 call vector_to_table_curr("map_table ", 'coef ', map_coor(1), i_map_coor)
                                 write(coeffCode,'(i1,a1,6i1)') i,'_',j(1:5),0   !anyway fixed above
                                 coeffCode = 'c'//coeffCode(1:8)//CHAR(0)
                                 call string_to_table_curr("map_table ", "name ", coeffCode);
                                 call augment_count("map_table ")
                               endif
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif (c_%npara.eq.4) then ! distribute exponents over 4 variables, knowing their sum
             do i1=no,0,-1
                do i2=no-i1,0,-1
                   do i3=no-i1-i2,0,-1
                      do i4=no-i1-i2-i3,0,-1
                         if (i1+i2+i3+i4==no) then
                            j(1)=i1
                            j(2)=i2
                            j(3)=i3
                            j(4)=i4
                            coef = y(i)%T.sub.j
                            if (coef.ne.zero) then
                               map_coor(1)=coef
                               map_coor(2)=i
                               map_coor(3)=c_%npara
                               map_coor(4)=no
                               map_coor(5)=j(1)
                               map_coor(6)=j(2)
                               map_coor(7)=j(3)
                               map_coor(8)=j(4)
                               map_coor(9)=0
                               map_coor(10)=0
                               call vector_to_table_curr("map_table ", 'coef ', map_coor(1), i_map_coor)
                               write(coeffCode,'(i1,a1,6i1)') i,'_',j(1:4),0,0
                               coeffCode = 'c'//coeffCode(1:8)//CHAR(0)
                               call string_to_table_curr( "map_table ", "name ", coeffCode );
                               call augment_count("map_table ")
                            endif
                         endif
                      enddo
                   enddo
                enddo
             enddo
          else
             call fort_warn('makemaptable ','map output expects 4,5 or 6 variables')
          endif
       enddo
    enddo

    deallocate(j)


  end subroutine makemaptable
  !_________________________________________________________________

  subroutine killsavedmaps
    implicit none
    integer i,ii

    if(.not. savemaps) return

    if (.not. associated(maps)) then
       return
    endif

    do i=lbound(maps,1),ubound(maps,1)
       do ii=1,6
          call kill(maps(i)%unimap(ii))
       enddo
    enddo
    deallocate(maps)
    nullify(maps)

  end subroutine killsavedmaps
  !_________________________________________________________________


  SUBROUTINE ptc_read_errors()
    use twtrrfi
    use name_lenfi
    implicit none
    integer i,k,pos,nfac(maxmul),flag,string_from_table_row,double_from_table_row,l
    real(dp) d(2*maxmul),b(maxmul),a(maxmul),tilt,ab,bvk
    character(name_len) name,name2
    type(fibre),pointer :: p
    logical(lp) :: overwrite
    real(kind(1d0)) get_value
    character*4 :: mag_index1(10)=(/'k0l ','k1l ','k2l ','k3l ','k4l ','k5l ','k6l ','k7l ','k8l ','k9l '/)
    character*5 :: mag_index2(10)=(/'k0sl ','k1sl ','k2sl ','k3sl ','k4sl ','k5sl ','k6sl ','k7sl ','k8sl ','k9sl '/)
    character*5 :: mag_index3(11)=(/'k10l ','k11l ','k12l ','k13l ','k14l ','k15l ','k16l ','k17l ','k18l ','k19l ','k20l '/)
    character*6 :: mag_index4(11)=(/'k10sl ','k11sl ','k12sl ','k13sl ','k14sl ','k15sl ','k16sl ', &
         'k17sl ','k18sl ','k19sl ','k20sl '/)

    overwrite = get_value('ptc_read_errors ','overwrite ').ne.0
    bvk=get_value('probe ','bv ')

    nfac(1)=1
    do i=2,maxmul
       nfac(i)=nfac(i-1)*(i-1)
    enddo

    flag = string_from_table_row('errors_read ', 'name ',1,name)

    if(flag.ne.0) call aafail('fill_errors reports: ',' The >>> errors_read <<< table is empty ')
    i=0

    p=>my_ring%start
    do while(.true.)
       i=i+1
       b(:)=zero
       a(:)=zero
       d(:)=zero
       name2=" "
       flag = string_from_table_row('errors_read ', 'name ',i,name2)
       if(flag.ne.0) goto 100
       do k=1,maxmul
          if(k<=10) then
             flag = double_from_table_row('errors_read ',mag_index1(k),i,d(2*k-1))
             flag = double_from_table_row('errors_read ',mag_index2(k),i,d(2*k))
          else
             flag = double_from_table_row('errors_read ',mag_index3(k-10),i,d(2*k-1))
             flag = double_from_table_row('errors_read ',mag_index4(k-10),i,d(2*k))
          endif
       enddo
       if(flag.ne.0) goto 100
       do k=1,maxmul
          b(k)=d(2*k-1)/nfac(k)
          a(k)=d(2*k)/nfac(k)
       enddo
       name=" "
       name(:len_trim(name2)-1)=name2(:len_trim(name2)-1)
       call context(name)
       call move_to(my_ring,p,name,pos)

       !madxtilt =  get_orginal_madx_tilt(name)

       tilt=-(p%mag%p%tiltd)! - madxtilt)   ! here we should read tilt from MADX lattice and deduce back the automatic tilt from skew (+ normal)

       if(pos/=0) then
          if(p%mag%l/=zero) then
             do k=1,maxmul
                b(k)=b(k)/p%mag%l
                a(k)=a(k)/p%mag%l
             enddo
          endif
          do k=1,maxmul
             b(k)=bvk*b(k)
             a(k)=bvk*a(k)
          enddo
          if(tilt/=zero) then
             do k=1,maxmul
                ab=b(k)
                b(k)=b(k)*cos(tilt*k)+a(k)*sin(tilt*k)
                a(k)=-ab*sin(tilt*k)+a(k)*cos(tilt*k)
             enddo
          endif
          do k=NMAX,1,-1
             if(b(k)/=zero) then
                if(overwrite) then
                   call add(p,k,0,b(k))
                else
                   call add(p,k,1,b(k))
                endif
             endif
             if(a(k)/=zero) then
                if(overwrite) then
                   call add(p,-k,0,a(k))
                else
                   call add(p,-k,1,a(k))
                endif
             endif
          enddo
       else
          write(6,*) " name,pos, dir of dna ",name, p%mag%parent_fibre%dir
       endif
    enddo
100 continue
    return

  end SUBROUTINE ptc_read_errors
  !________________________________________________________________________________________________

  subroutine ptc_refresh_k()
    use twtrrfi
    use name_lenfi
    implicit none
    integer j,code,k,pos,nfac(maxmul)
    integer restart_sequ,advance_node
    type(fibre),pointer :: p
    real(dp) sk0,sk,sks,tilt,b(maxmul),a(maxmul),bvk
    real(kind(1d0))   :: get_value,node_value
    character(name_len) name
    logical(lp) :: overwrite
    !---------------------------------------------------------------

    overwrite = get_value('ptc_refresh_k ','overwrite ').ne.0
    bvk=get_value('probe ','bv ')

    nfac(1)=1
    do j=2,maxmul
       nfac(j)=nfac(j-1)*(j-1)
    enddo

    j=restart_sequ()
    j=0
    p=>my_ring%start
10  continue
    b(:)=zero
    a(:)=zero

    code=node_value('mad8_type ')
    if(code.ne.5.and.code.ne.6) goto 100
    if(code.eq.5) then
       ! quadrupole components code =  5
       k=2
       sk= node_value('k1 ')
       sks=node_value('k1s ')
       tilt=node_value('tilt ')
       b(k)=sk
! LD: 19.06.2019
       sk0= node_value('k0 ')
       if (sks .ne. zero .and. sk0 .eq. zero) then ! should also consider permfringe
          tilt = -atan2(sks, sk)/two + tilt
          b(k)=sqrt(sk**2+sks**2)/abs(sk)*sk
          ! bug: sks not updated
       endif
    elseif(code.eq.6) then
       ! sextupole components code = 6
       k=3
       sk= node_value('k2 ')+node_value('k2tap ')
       sks=node_value('k2s ')
       tilt=node_value('tilt ')
       b(k)=sk
! LD: 19.06.2019
!       if (sks .ne. zero) then
!          tilt = -atan2(sks, sk)/three + tilt
!          b(k)=sqrt(sk**2+sks**2)/abs(sk)*sk                           !
!       endif
    endif

    call element_name(name,name_len)
    call context(name)
    call move_to(my_ring,p,name,pos)
    if(pos/=0) then
       b(k)=b(k)/nfac(k)
       if(tilt/=zero) then
          a(k)=-b(k)*sin(tilt*k)
          b(k)=b(k)*cos(tilt*k)
       endif
       b(k)=bvk*b(k)
       a(k)=bvk*a(k)
       do j=1,maxmul
          if(overwrite) then
             call add(p,j,0,b(j))
             call add(p,-j,0,a(j))
          else
             call add(p,j,1,b(j))
             call add(p,-j,1,a(j))
          endif
       enddo
    else
       write(6,*) " name,pos, dir of dna ",name, p%mag%parent_fibre%dir
    endif

100 continue

    if(advance_node().ne.0)  goto 10

    return

  END subroutine ptc_refresh_k
  !________________________________________________________________________________________________

  subroutine getfk(fk)
  !returns FK factor for Beam-Beam effect
    implicit none
    real (dp) :: fk,dpp
    real (dp) :: gamma0,beta0,beta_dp,ptot,b_dir,arad,totch
    real (dp) :: q,q_prime
    integer   :: b_dir_int
    real(kind(1d0)) :: get_value
    real(kind(1d0)) :: get_variable
    integer         :: get_option
    REAL(KIND(1d0)) :: node_value  !/*returns value for parameter par of current element */

    !---- Calculate momentum deviation and according changes
    !     of the relativistic factor beta0

    gamma0 = get_value('probe ','gamma ')
    arad=get_value('probe ', 'arad ')
    totch=node_value('charge ') * get_value('probe ', 'npart ')

    if (getdebug()>1) then
      print*, 'getfk for beam-beam: charge npart ',node_value('charge '), get_value('probe ', 'npart ')
      print*, 'getfk for beam-beam: gamma0, arad, totch ',gamma0, arad, totch
    endif


    dpp  = get_variable('track_deltap ')
    q = get_value('probe ','charge ')
    q_prime = node_value('charge ')

    if (getdebug()>1) then
      print*, 'dpp q q_prime',dpp, q, q_prime
    endif

    beta0 = sqrt(one-one/gamma0**2)
    ptot = beta0*gamma0*(one+dpp)
    beta_dp = ptot / sqrt(one + ptot**2)
    b_dir_int = node_value('bbdir ')
    b_dir=dble(b_dir_int)
    b_dir = b_dir/sqrt(b_dir*b_dir + 1.0d-32)

    fk = two*arad*totch/gamma0 /beta0/(one+dpp)/q*          &
         (one-beta0*beta_dp*b_dir)/(beta_dp+0.5*(b_dir-one)*b_dir*beta0)

  end subroutine getfk

  !____________________________________________________________________________________________
  ! Configures beam-beam for every beambeam element defined in MADX lattice
  subroutine getBeamBeam()
   implicit none
   integer                 :: i,e,elcode
   integer, external       :: restart_sequ, & !  restart beamline and return number of beamline node
                              advance_node    !  advance to the next node in expanded sequence
                                              !  =0 (end of range), =1 (else)
   double precision, external :: node_value  !/*returns value for parameter par of current element */
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
          if (getdebug() > -2 ) then
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
  !________________________________________________________________________________________________
  ! this routine serves ptc_putbeambeam command
  ! that enables user to inject beambeam inside any element without splitting it
  subroutine putbeambeam()
    implicit none
    real (dp) :: fk, xma, yma, sigx, sigy,s
    integer :: elno
    logical(lp) :: found
    TYPE(INTEGRATION_NODE),POINTER :: node

    real(kind(1d0)), external :: get_value
    real(kind(1d0)), external :: get_variable
    integer, external         :: get_option
    REAL(KIND(1d0)), external :: node_value  !/*returns value for parameter par of current element */

    s    = get_value('ptc_putbeambeam ','global_s ')
    xma  = get_value('ptc_putbeambeam ','xma ')
    yma  = get_value('ptc_putbeambeam ','yma ')
    sigx = get_value('ptc_putbeambeam ','sigx ')
    sigy = get_value('ptc_putbeambeam ','sigy ')

    print*, 'Input   xma, yma, sigx, sigy, s'
    print*, 'Input', xma, yma, sigx, sigy, s

    if(.not.associated(my_ring%t))  then
       CALL MAKE_node_LAYOUT(my_ring)
    endif

    elno = 0
    call s_locate_beam_beam(my_ring,s,elno,node,found)

    if (.not. found) then
      print*,"could not find a node for beam-beam"
      return
    endif

    print*, 'Name of element in PTC: ', node%PARENT_FIBRE%mag%name

    write(6,*) 'node%a:',node%a            ! node entrance position
    write(6,*) 'node%ent(1,1:3):',node%ent(1,1:3)   ! node entrance e_1 vector
    write(6,*) 'node%ent(2,1:3):',node%ent(2,1:3)   ! node entrance e_2 vector
    write(6,*) 'node%ent(3,1:3):',node%ent(3,1:3)   ! node entrance e_3 vector
    write(6,*) " s variable of node and following node and "
    write(6,*) s, node%s(1),node%next%s(1), s - node%s(1), s - node%next%s(1)



    if(((node%cas==case0).or.(node%cas==caset))) then !must be 0 or 3

      if(.not.associated(node%BB)) call alloc(node%BB)

      call getfk(fk)
      node%bb%fk = fk
      node%bb%xm = xma
      node%bb%ym = yma
      node%bb%sx = sigx
      node%bb%sy = sigy
      node%bb%PATCH=.true.
      if (getdebug() > 2 ) then
        print*, "BB fk=",node%bb%fk
        print*, "BB sx=",node%bb%sx
        print*, "BB sy=",node%bb%sy
        print*, "BB xm=",node%bb%xm
        print*, "BB ym=",node%bb%ym
      endif

      do_beam_beam = .true.

    else
      call fort_warn('getBeamBeam: ','Bad node case for BeamBeam')
    endif

    !N.B. If nothing else is done, the beam-beam kick is placed at the entrance of the node.
    !The call FIND_PATCH(t%a,t%ent,o ,mid,D,ANG) needs to be invoked to place the beam-beam kick

  end subroutine putbeambeam
  !________________________________________________________________________________________________

  ! if clock with such freqency exists it returns its index
  ! if not, returns index of the next free slot
  ! if there is no free slots, returns -1
  integer function getclockidx(f)
    implicit none
    real(dp) f ! frequency to search
    integer i, r1, r2, r3, r4
    logical fits
    real(kind(1d0)) node_value
    getclockidx = -1

    ! frequency is in fact tune
    ! kept like this on Rogelio request not to break the codes before LS2
    ! afterwards "freq" should be changed to "tune" in definition of the AC_DIPOLE
    

    r1 = node_value('ramp1 ')
    r2 = node_value('ramp2 ')
    r3 = node_value('ramp3 ')
    r4 = node_value('ramp4 ')

    do i=1,nclocks

     if ( abs(clocks(i)%tune - f) .gt. c_1d_10 )   cycle
     if (r1 .ne.  clocks(i)%rampupstart )          cycle
     if (r2 .ne.  clocks(i)%rampupstop )           cycle
     if (r3 .ne.  clocks(i)%rampdownstart )        cycle
     if (r4 .ne.  clocks(i)%rampdownstop )         cycle

     ! this is the good clock
     getclockidx = i
     return

    enddo


    if (nclocks == nmaxclocks) then
      getclockidx = -1 ! repeated for code clarity
      return
    endif

    nclocks = nclocks + 1

    clocks(nclocks)%tune          = f
    clocks(nclocks)%rampupstart   = r1
    clocks(nclocks)%rampupstop    = r2
    clocks(nclocks)%rampdownstart = r3
    clocks(nclocks)%rampdownstop  = r4

    getclockidx = nclocks

    clocks(nclocks)%nelements = 0
    
    if (getdebug() > 1) then
      print*,"getclockidx: Created new clock. nclocks = ", nclocks
    endif


  end function getclockidx
  !________________________________________________________________________________________________

  subroutine addelementtoclock(p,c)
    implicit none
    type(fibre), pointer :: p
    integer c  ! clock index
    integer elidx

     if (clocks(c)%nelements .ge. maxelperclock) then
       call aafail('ptc_input:addelementtoclock:', &
        'Buffer for AC dipoles is too small. Contact MADX support to make it bigger.')
     endif

     clocks(c)%nelements = clocks(c)%nelements + 1
     elidx = clocks(c)%nelements

     clocks(c)%elements(elidx)%p=>p
     
     ! sets amplitude of modulation to maximum for ptc_twiss
     ! (in track this parameter is ramped up and down)
     p%magp%d_ac = 1

  end subroutine addelementtoclock

  !________________________________________________________________________________________________
  !
  subroutine acdipoleramping(t)
    implicit none
    !---------------------------------------    *
    !--- ramp up and down of the ac dipols      *
    !--- Adjust amplitudes in function of turns *
    integer  t
    integer  n,i
    real(dp) r
    type(fibre), pointer :: p
    
    !print*,"acdipoleramping t=",t
    
    do n=1,nclocks

      do i=1,clocks(n)%nelements

        p => clocks(n)%elements(i)%p

        !print*,"Setting ramp to clock ",n," element ", p%mag%name

        if (clocks(n)%rampupstop < 1) then
          ! no ramping, always full amplitude
          p%mag%d_ac = one
          cycle
        endif

        if (t < clocks(n)%rampupstart) then
          p%mag%d_ac = zero
          cycle
        endif

        if (t < clocks(n)%rampupstop) then
          r = (t - clocks(n)%rampupstart)
          p%mag%d_ac = r/(clocks(n)%rampupstop - clocks(n)%rampupstart)
          cycle
        endif

        if (t < clocks(n)%rampdownstart) then
          p%mag%d_ac = one
          cycle
        endif

        if (t < clocks(n)%rampdownstop) then
          r = (clocks(n)%rampdownstop - t)
          p%mag%d_ac = r/(clocks(n)%rampdownstop - clocks(n)%rampdownstart)
          cycle
        endif

        p%mag%d_ac = zero

      enddo
     

    enddo
    
    !print*,"acdipoleramping d_ac=",p%mag%d_ac
  end subroutine acdipoleramping

  !_________________________________________
  ! returns max multipole order of bends
  ! it is used to set
  function getmaxnmul()
    use twtrrfi, only: maxferr
    implicit none
    integer getmaxnmul
    integer i,j, maxnmul, maxk, code, n_ferr, max_n_ferr
    integer restart_sequ,advance_node,node_fd_errors
    integer n_norm, n_skew
    REAL(dp) :: tmp_0123(0:3), v
    REAL(dp) :: tmpmularr(0:maxferr), field(2,0:maxferr)
    real(kind(1d0)) node_value
    character (len = 3), dimension(3) :: kns = (/'k1 ','k2 ','k3 ' /)
    character (len = 4), dimension(3) :: kss = (/'k1s ','k2s ','k3s ' /)

    getmaxnmul = -1

    j=restart_sequ() ! returns -1 if error, +1 if OK


    do while (j .gt. 0)
      code=node_value('mad8_type ')
      !bend or rbend
      if (code/=2 .and. code/=3 ) then
        j = advance_node()  !returns 1 if OK, 0 otherhise
        cycle;
      endif

      maxk = 0
      do i=3,1,-1
        v = node_value(kns(i))
        if (v .ne. zero ) then
          maxk = i
          exit
        endif

        v = node_value(kss(i))
        if (v .ne. zero ) then
          maxk = i
          exit
        endif

      enddo


      call get_node_vector('knl ',n_norm,tmpmularr)
      call get_node_vector('ksl ',n_skew,tmpmularr)

      n_ferr = node_fd_errors(tmpmularr)
      field = zero ! array to be zeroed.                          !
      if (n_ferr .gt. 0) then                                     !
         call dcopy(tmpmularr,field,n_ferr)                        !
      endif                                                       !

      max_n_ferr = n_ferr/2
      n_ferr = -1
      do i=max_n_ferr,0,-1
      !  print*,"magnet =",j," err_n ", i, " ", field(1,i),  field(2,i)
        IF( (field(1,i)/=0.0_dp .or. field(2,i)/=0.0_dp) ) THEN
          n_ferr = i+1
          exit
        ENDIF
      enddo


      maxnmul = max(n_norm,n_skew)  ! max order between  knl and kns
      maxnmul = max(maxk,  maxnmul) ! max between the above and defined with k1,k2,k3, k1s,k2s,k3s
      maxnmul = max(n_ferr,maxnmul) ! max between the above and the errors

      if (maxnmul > getmaxnmul) getmaxnmul = maxnmul
      if (getdebug() > 2) then
         print*, "j   getmaxnmul   maxnmul  maxk  n_norm  n_skew  n_ferr"
         print*, j,getmaxnmul, maxnmul, maxk, n_norm,n_skew,n_ferr
      endif

      j = advance_node()

    enddo

    !
    getmaxnmul = getmaxnmul + 1

    if (getdebug() > 0) then
      write(6,'(a,i2)') "Determined SECTOR NMUL MAX : ", getmaxnmul
    endif


  end  function getmaxnmul


END MODULE madx_ptc_module