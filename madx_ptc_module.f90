MODULE ptc_results
  USE madx_keywords
  implicit none
  public
  integer :: number_variables = 6
  integer :: order = 20
  character(len = 2), dimension(6) :: ptc_variables = (/'x ','xp','y ','yp','z ','dp'/)
  character(len = 2) :: ptc_var
  type(real_8) y(6)
  type(normalform) n
  type (pbresonance) pbrg,pbrh
END MODULE ptc_results

MODULE madx_ptc_module
  USE madx_keywords
  USE madx_ptc_setcavs_module
  USE madx_ptc_tablepush_module
  use madx_ptc_intstate_module, only : getdebug

  implicit none
  public
  logical(lp) mytime
  integer icav
  integer :: universe=0,index=0,EXCEPTION=0
  integer ipause
  integer,external :: mypause
  real(kind(1d0)) get_value,node_value
  type(layout),pointer :: MY_RING
  type(mad_universe) m_u
  integer, private, parameter :: ndd=ndim2,mynreso=20
  integer, private, dimension(4) :: iia,icoast
  integer, private, allocatable :: J(:)
  integer, private :: NO,ND,ND2,NP,NDPT,NV,icount=0
  real(dp), private, dimension(ndim2,5) :: rdd
  real(dp), private, dimension(ndim2) :: dicu
  real(dp), private, dimension(2,ndim2) :: angp
  character(len=4), private, dimension(4), target :: str4 = (/'1000','0100','0010','0001'/)
  character(len=5), private, dimension(5), parameter :: str5 = (/'10000','01000','00100','00010','00001'/)
  character(len=6), private, dimension(6), target :: str6 = (/'100000','010000','001000','000100','000001','000010'/)
  private zerotwiss,equaltwiss,alloctwiss,killtwiss
  real(dp) :: mux_default=c_0_28, muy_default=c_0_31, muz_default=c_1d_3
  type twiss
     type(damap) a1
     type(damap) a_t
     type(damap) junk
     type(damap) junk1
     type(normalform) n
     logical(lp) nf
     real(dp), dimension(3,3) ::  beta,alfa,gama
     real(dp), dimension(3)   ::  mu
     real(dp), dimension(6)   ::  disp
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

CONTAINS

  subroutine ptc_create_universe()
    implicit none

    if (getdebug()>1) print*,"Now PTC"
    
    print77=.false.
    read77 =.false.

    if (getdebug()==0) global_verbose = .false.
    if (getdebug()>0) print*,"Now PTC"

    call set_up_universe(m_u)
    universe=universe+1
    
    
  end subroutine ptc_create_universe
  !_________________________________________________________________

  subroutine ptc_create_layout()
    implicit none
    real(kind(1d0)) get_value

    if(universe.le.0) then
       call fort_warn('return from ptc_create_layout: ',' no universe created')
       return
    endif

    call append_empty_layout(m_u)
    index=index+1
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

    if(universe.le.0) then
       call fort_warn('return from ptc_move_to_layout: ',' no universe created')
       return
    endif

    my_index = get_value('ptc_move_to_layout ','index ')

    if(my_index.gt.index.or.my_index.le.0) then
       call fort_warn('return from ptc_move_to_layout: ',' layout outside allowed range')
       print*,"   Allowed range 0 < ",index
       return
    endif

    call move_to_layout_i(m_u,my_ring,my_index)

  end subroutine ptc_move_to_layout
  !_________________________________________________________________

  subroutine ptc_input()
    implicit none
    include 'twtrr.fi'
    include 'name_len.fi'
    include 'twiss0.fi'
    logical(lp) particle,doneit,isclosedlayout
    integer i,j,k,code,nt,icount,nn,ns,nd
    integer get_option,double_from_table
    integer restart_sequ,advance_node,n_ferr,node_fd_errors
    integer, parameter :: nt0=20000,length=16
    real(dp) l,l_machine,energy,kin,brho,beta0,p0c,pma,e0f,lrad
    real(dp) f_errors(0:50),aperture(maxnaper),normal(0:maxmul)
    real(dp) patch_ang(3),patch_trans(3)
    real(dp) skew(0:maxmul),field(2,0:maxmul),fieldk(2)
    real(dp) gamma,gamma2,gammatr2,freq,offset_deltap
    real(dp) fint,fintx,div,muonfactor
    real(dp) sk1,sk1s,sk2,sk2s,sk3,sk3s,tilt
    REAL(dp) ::  normal_0123(0:3), skew_0123(0:3) ! <= knl(1), ksl(1)
    real(kind(1d0)) get_value,node_value,gammatr
    character(length) name
    character(name_len) aptype
    type(keywords) key
    character(20)       keymod0,keymod1
    logical(lp)         exact0
    integer             exact1
    integer             sector_nmul_max0,sector_nmul0,sector_nmul1
    integer             model
    integer             method0,method1
    integer             nst0,nst1
    REAL (dp) :: tempdp
    !---------------------------------------------------------------
    if (getdebug() > 1) then
       print *, '--------------------------------------------------------------'
       print *, '--------------------------------------------------------------'
       print *, '------    E X E C U T I N G     P T C     I N P U T   --------'
       print *, '--------------------------------------------------------------'
       print *, '--------------------------------------------------------------'
    endif

    energy=get_value('beam ','energy ')
    pma=get_value('beam ','mass ')
    e0f=sqrt(ENERGY**2-pma**2)

    if (getdebug() > 0) then
       print *, 'MAD-X Beam Parameters'
       print '(a20, f8.4)', '      Energy :',energy
       print '(a20, f8.4)', '      Kinetic Energy :',energy-pma
       print '(a20, f8.4)', '      Partice Rest Mass :',pma
       print '(a20, f8.4)', '      Momentum :',e0f
    endif



    beta0=e0f/ENERGY


    if(abs(pma-pmae)/pmae<c_0_002) then
       if (getdebug() > 1) print *,'Executing MAKE_STATES(TRUE), i.e. ELECTRON beam'
       particle=.true.
       CALL MAKE_STATES(PARTICLE)
    elseif(abs(pma-pmap)/pmap<c_0_002) then
       if (getdebug() > 1) print *,'Executing MAKE_STATES(FALSE), i.e. PROTON beam'
       particle=.false.
       CALL MAKE_STATES(PARTICLE)
    else
       if (getdebug() > 1) print '(a, f8.4, a)','Executing MAKE_STATES(',pma/pmae,'), i.e. PROTON beam'
       muonfactor=pma/pmae
       CALL MAKE_STATES(muonfactor)
    endif

    !valid October 2002: oldscheme=.false.
    !!valid October 2002: oldscheme=.true.

    if (getdebug() > 1) print '(a23, l7, a1)','Executing MAKE_STATES(',PARTICLE,')'

    CALL MAKE_STATES(PARTICLE)

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

    sector_nmul_max0 = get_value('ptc_create_layout ','sector_nmul_max ')
    if (getdebug() > 1) print*,'  Global max sector_nmul: ',sector_nmul_max0

    sector_nmul0 = get_value('ptc_create_layout ','sector_nmul ')
    if (getdebug() > 1) print*,'  Global sector_nmul: ',sector_nmul0


    model = get_value('ptc_create_layout ','model ')
    if (getdebug() > 1) print*,'  Global Model code is : ',model

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
       ipause=mypause(444)
       RETURN
    END SELECT



    if (getdebug() > 1) print*,'  Global Model name (keymod0) is : ',keymod0

    method0   = get_value('ptc_create_layout ','method ')
    if (getdebug() > 1) print*,'  Global method is: ',method0

    exact0    = get_value('ptc_create_layout ','exact ') .ne. 0
    if (getdebug() > 1) print*,'  Global exact is: ',exact0

    nst0      = get_value('ptc_create_layout ','nst ')
    if (getdebug() > 1) print*,'  Global Number of Integration Steps (nst) is: ',nst0

    ! MAD-X specials
    madlength = get_option('rbarc ') .eq. 0
    if (getdebug() > 1) print*,'  global rbend_length: ',madlength

    mad       = get_value('ptc_create_layout ','mad_mult ') .ne. 0
    if (getdebug() > 1) print*,'  global mad_mult as in mad8: ',mad

    mad8      = get_value('ptc_create_layout ','mad8 ') .ne. 0
    if (getdebug() > 1) print*,'  rbend as in mad8 (only global): ',mad8

    gamma     = get_value('beam ','gamma ')
    if (getdebug() > 1) print*,'  gamma: ',gamma

    k         = double_from_table('summ ','gammatr ',1,gammatr)
    if (getdebug() > 1) print*,'  gammatr: ',gammatr

    gamma2    = gamma**2
    gammatr2  = gammatr**2
    
    if (getdebug() > 1) then
      print *, '=============================================================='
      print *, ''
    endif

    !  call Set_Up(MY_RING)

    if (getdebug() > 0) then
       print *, 'Setting MADx with '
       print *, '    energy        ',energy
       print *, '    method        ',method0
       print *, '    Num. of steps ',nst0
    endif

    CALL SET_MADx(energy=energy,METHOD=method0,STEP=nst0)
    
    if (getdebug() > 1) print *, 'MADx is set'

    icav=0
    nt=0
    j=restart_sequ()
    j=0
    l_machine=zero
10  continue
    nst1=node_value("nst ")
    if(nst1.gt.0) then
       nstd = nst1
    else
       nstd = nst0
    endif

    call zero_key(key)
    j=j+1
    nt=nt+1
    if(nt==nt0) then
       print*,'More than the maximum number: ',nt0,' of elements in the structure==> Program stops'
       stop
    endif
    icount=0
    l=zero
    l=node_value('l ')
    key%list%l=l
    l_machine=l_machine+l
    code=node_value('mad8_type ')
    call element_name(name,length)
    key%list%name=name

    !local, if present, superseed global at current node
    sector_nmul1=node_value("sector_nmul ")
    if(sector_nmul1.gt.0) then
       sector_nmul = sector_nmul1
    else
       sector_nmul = sector_nmul0
    endif


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
    key%list%permfringe=node_value("permfringe ") .ne. zero
    key%list%kill_ent_fringe=node_value("kill_ent_fringe ") .ne. zero
    key%list%kill_exi_fringe=node_value("kill_exi_fringe ") .ne. zero
    key%list%bend_fringe=node_value("bend_fringe ") .ne. zero

    nn=name_len
    call node_string('apertype ',aptype,nn)
    call dzero(aperture,maxnaper)
    call get_node_vector('aperture ',nn,aperture)
    if(.not.((aptype.eq."circle".and.aperture(1).eq.zero).or.aptype.eq." ")) then
       c_%APERTURE_FLAG=.true.
       select case(aptype)
       case("circle")
          key%list%aperture_on=.true.
          key%list%aperture_kind=1
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(1)
       case("ellipse")
          key%list%aperture_on=.true.
          key%list%aperture_kind=1
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(2)
       case("rectangular")
          key%list%aperture_on=.true.
          key%list%aperture_kind=2
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
       case("lhcscreen")
          key%list%aperture_on=.true.
          key%list%aperture_kind=3
          key%list%aperture_x=aperture(3)
          key%list%aperture_y=aperture(3)
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(2)
       case("marguerite")
          key%list%aperture_on=.true.
          key%list%aperture_kind=4
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(2)
       case("general")
          key%list%aperture_kind=5
          print*,"General aperture not implemented"
          stop
       end select
    endif
    call append_empty(my_ring)


    !    print *,'____________________________________________________'
    !    print *,'Adding an element with code',code,' named ',name
    !    print *,'____________________________________________________'

    select case(code)
    case(0,4,25)
       key%magnet="marker"
    case(1,11,20,21)
       key%magnet="drift"
!       if (getdebug() > 2)  print *, 'This is a drift'
    case(2) ! PTC accepts mults
       if(l.eq.zero) then
          key%magnet="marker"
          goto 100
       endif
       key%magnet="rbend"
           !VK
           CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123) 

          tempdp=sqrt(normal_0123(0)*normal_0123(0)+skew_0123(0)*skew_0123(0))
          key%list%b0=node_value('angle ')+tempdp*l

          key%list%k(2)=node_value('k1 ')+ key%list%k(2) 
          key%list%k(3)=node_value('k2 ')+ key%list%k(3) 
          key%list%k(4)=node_value('k3 ')+ key%list%k(4)

          key%list%ks(2)=node_value('k1s ')+ key%list%ks(2) 
          key%list%ks(3)=node_value('k2s ')+ key%list%ks(3) 
          key%list%ks(4)=node_value('k3s ')+ key%list%ks(4)
 
       ! Gymnastic needed since PTC expects MAD8 convention
       key%list%t1=node_value('e1 ')-node_value('angle ')/two
       key%list%t2=node_value('e2 ')-node_value('angle ')/two
       key%list%hgap=node_value('hgap ')
       !       key%list%fint=node_value('fint ')
       fint=node_value('fint ')
       fintx=node_value('fintx ')
       if((fintx.ne.fint).and.(fintx.gt.zero.and.fint.gt.zero)) then
          print*," The fint and fintx must be the same at each end or each might be zero"
          stop
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
          if(tempdp.gt.0) key%tiltd=key%tiltd + asin(skew_0123(0)/tempdp)       
    case(3) ! PTC accepts mults watch out sector_nmul defaulted to 4
       if(l.eq.zero) then
          key%magnet="marker"
          goto 100
       endif
       key%magnet="sbend"
           !VK
           CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123) 

          tempdp=sqrt(normal_0123(0)*normal_0123(0)+skew_0123(0)*skew_0123(0))
          key%list%b0=node_value('angle ')+ tempdp*l

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
          stop
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
          if(tempdp.gt.0) key%tiltd=key%tiltd + asin(skew_0123(0)/tempdp)

   case(5) ! PTC accepts mults
       if (getdebug() > 9)  print *, 'This is a quadrupole'
       key%magnet="quadrupole"

       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123) 
       ! Read data & fill %k(:), %ks(:) arrays which are 
       ! summs of multipoles and errors

       ! quadrupole components 
       sk1=node_value('k1 ')
       sk1s=node_value('k1s ')

       ! A sum of quadrupole components from K1 & K1S and ========!
       ! from multipoles on the bench (without errors) defines    !
       ! a tilt angle of normal Q                                 !
       if(l.ne.0) then                                            !
         !sk1 = sk1 +  normal(1)/l                                !
          sk1 = sk1 +  normal_0123(1)                             !
         !sk1s = sk1s + skew(1)/l                                 !
          sk1s = sk1s + skew_0123(1)                              !
       endif                                                      !
       if (sk1s .eq. zero)  then                                  !
          tilt = zero                                             !
       else                                                       !
          tilt = asin(sk1s/sqrt(sk1**2 + sk1s**2)) / two          !
       endif                                                      ! 
                                                                  !
       if(l.ne.0) then                                            !
         !sk1  = sk1  + field(1,1)/l                              !
          sk1  = sk1  + (key%list%k(2)-normal_0123(1))            !
         !sk1s = sk1s + field(2,1)/l                              !
          sk1s = sk1s + (key%list%ks(2)-skew_0123(1))             !
       endif                                                      !
                                                                  !
       if (tilt .ne. zero) sk1 = sqrt(sk1**2 + sk1s**2)           !
       key%list%k(2)=sk1                                          !
       key%list%ks(2)=zero  ! added by VK                         !       
       key%tiltd=node_value('tilt ')+tilt  !======================!

   !================================================================

    case(6)
       key%magnet="sextupole"

       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123) !VK

       sk2=node_value('k2 ')
       sk2s=node_value('k2s ')

       ! A sum of sextupole components from K2 & K2S and  ========!
       ! from multipoles on the bench (without errors) defines    !
       ! a tilt angle of normal Sextupole                         !
       if(l.ne.0) then                                            !
         !sk2 = sk2 +  normal(2)/l                                !
          sk2 = sk2 +  normal_0123(2)                             !
         !sk2s = sk2s + skew(2)/l                                 !
          sk2s = sk2s + skew_0123(2)                              !
       endif                                                      !
       !                                                          !
       if (sk2s .eq. zero)  then                                  !
          tilt = zero                                             !
       else                                                       !
          tilt = asin(sk2s/sqrt(sk2**2 + sk2s**2)) / three        !
       endif                                                      !
       !                                                          !
       if(l.ne.0) then                                            !
         !sk2  = sk2 + field(1,2)/l                               !
          sk2  = sk2 + (key%list%k(3)-normal_0123(2))             !
         !sk2s = sk2s + field(2,2)/l                              !
          sk2s = sk2s + (key%list%ks(3)-skew_0123(2))             !
       endif                                                      !
       if (tilt .ne. zero) sk2 = sqrt(sk2**2 + sk2s**2)           ! 
       key%list%k(3)=sk2                                          !
       key%list%ks(3)=zero  ! added by VK                         ! 
       key%tiltd=node_value('tilt ')+tilt !-----------------------!

    case(7) ! PTC accepts mults
       key%magnet="octupole"

       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123) !VK

       sk3=node_value('k3 ')
       sk3s=node_value('k3s ')

       ! A sum of octupole components from K3 & K3S and ==========!
       ! from multipoles on the bench (without errors) defines    !
       ! a tilt angle of normal Octupole                          !
       if(l.ne.0) then                                            !
         !sk3 = sk3 +  normal(3)/l                                !
          sk3 = sk3 +  normal_0123(3)                             !
         !sk3s = sk3s + skew(3)/l                                 !
          sk3s = sk3s + skew_0123(3)                              !
       endif                                                      !
       !                                                          !
       if (sk3s .eq. zero)  then                                  !
          tilt = zero                                             !
       else                                                       !
          tilt = asin(sk3s/sqrt(sk3**2 + sk3s**2)) / four         !
       endif                                                      !
       !                                                          !
       if(l.ne.0) then                                            !
         !sk3 = sk3 + field(1,3)/l                                !
          sk3 = sk3 + (key%list%k(4)-normal_0123(3))              !
         !sk3s = sk3s + field(2,3)/l                              !
          sk3s = sk3s + (key%list%ks(3)-skew_0123(3))             !
       endif                                                      !
       if (tilt .ne. zero) sk3 = sqrt(sk3**2 + sk3s**2)           !
       key%list%k(4)=sk3                                          !
       key%list%ks(4)=zero  ! added by VK                         ! 
       key%tiltd=node_value('tilt ')+tilt !-----------------------!

    case(8)
       key%magnet="multipole"
       !---- Multipole components.
       call dzero(f_errors,maxferr+1)
       n_ferr = node_fd_errors(f_errors)
       call dzero(normal,maxmul+1)
       call dzero(skew,maxmul+1)
       call get_node_vector('knl ',nn,normal)
       call get_node_vector('ksl ',ns,skew)
       key%list%thin_h_angle=normal(0)
       key%list%thin_v_angle=skew(0)
       lrad=node_value('lrad ')
       if(lrad.gt.zero) then
          key%list%thin_h_foc=normal(0)*normal(0)/lrad
          key%list%thin_v_foc=skew(0)*skew(0)/lrad
       endif
       if(nn.gt.0) then
          do i=1,nn
             key%list%k(i+1)=normal(i)
          enddo
       endif
       if(ns.gt.0) then
          do i=1,ns
             key%list%ks(i+1)=skew(i)
          enddo
       endif
       call dzero(field,2*(maxmul+1))
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
    case(9) ! PTC accepts mults
       key%magnet="solenoid"
       key%list%bsol=node_value('ks ')
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123) !VK

    case(10)
       key%magnet="rfcavity"
       key%list%volt=node_value('volt ')
       freq=c_1d6*node_value('freq ')
       key%list%lag=node_value('lag ')*twopi
       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       if(offset_deltap.ne.zero) then
          default = getintstate()
          default=default+totalpath0
          call setintstate(default)
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one)
       endif
       key%list%freq0=freq
       key%list%n_bessel=node_value('n_bessel ')
       key%list%harmon=one
       if(key%list%volt.ne.zero.and.key%list%freq0.ne.zero) icav=1
       !  case(11)
       !     key%magnet="elseparator"
       !     key%list%volt=node_value('ex ')
       !     key%list%lag=atan2(node_value('ey '),node_value('ex '))
       !     key%tiltd=node_value('tilt ')
    case(14,15,16) ! PTC accepts mults
       call dzero(f_errors,maxferr+1)
       n_ferr = node_fd_errors(f_errors)
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
          key%list%k(1)=node_value('kick ')+node_value('chkick ')+fieldk(1)/div
       else if(code.eq.15) then
          key%magnet="kicker"
          key%list%k(1)=node_value('hkick ')+node_value('chkick ')+fieldk(1)/div
          key%list%ks(1)=node_value('vkick ')+node_value('cvkick ')+fieldk(2)/div
       else if(code.eq.16) then
          key%magnet="vkicker"
          key%list%ks(1)=node_value('kick ')+node_value('cvkick ')+fieldk(2)/div
       else
          key%magnet="marker"
       endif
       key%tiltd=node_value('tilt ')
    case(17)
       key%magnet="hmonitor"
    case(18)
       key%magnet="monitor"
    case(19)
       key%magnet="vmonitor"
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
    case(24)
       key%magnet="instrument"
       key%tiltd=node_value('tilt ')
    case(27)
!       if (getdebug() > 2)  print *, 'This is a twcavity'
       key%magnet="twcavity"
       key%list%volt=node_value('volt ')
       freq=c_1d6*node_value('freq ')
       key%list%lag=node_value('lag ')*twopi
       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       if(offset_deltap.ne.zero) then
          default = getintstate()
          default=default+totalpath0
          call setintstate(default)
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one)
       endif
       key%list%freq0=freq
       key%list%dphas=node_value("delta_lag ")
       key%list%psi=node_value("psi ")
       key%list%harmon=one
       if(key%list%volt.ne.zero.and.key%list%freq0.ne.zero) icav=1
    case(35)
       key%magnet="CHANGEREF"
       call dzero(patch_ang,3)
       call dzero(patch_trans,3)
       call get_node_vector('patch_ang ',3,patch_ang)
       call get_node_vector('patch_trans ',3,patch_trans)
       key%list%patchg=2
       do i=1,3
          key%list%ang(i)=patch_ang(i)
          key%list%t(i)=patch_trans(i)
       enddo
    case default
       print*,"Element: ",name," not implemented"
       stop
    end select
100 continue
    call create_fibre(my_ring%end,key,EXCEPTION)
    if(advance_node().ne.0)  goto 10

    if (getdebug() > 0) then
      print*,' Length of machine: ',l_machine
    endif
       
    CALL GET_ENERGY(ENERGY,kin,BRHO,beta0,P0C)

    isclosedlayout=get_value('ptc_create_layout ','closed_layout ') .ne. 0
    if (isclosedlayout .eqv. .true.) then
       print *,'The machine is a RING'
    else
       print *,'The machine is a LINE'
    endif

    MY_RING%closed=isclosedlayout

    doneit=.true.
    call ring_l(my_ring,doneit)

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
    
    if (getdebug() > 1) then 
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^    F I N I S H E D      P T C     I N P U T    ^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
    endif

  END subroutine ptc_input
  !_________________________________________________________________

  SUBROUTINE SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123, skew_0123)
  ! 1) read multipole coeff. and errors for a current thick element
  ! 2) fill the error and multiploes arrays of data-bases
    include 'twtrr.fi' ! integer, maxmul,maxferr,maxnaper
    REAL(dp), INTENT(IN) :: l
    type(keywords), INTENT(INOUT) ::  key    
    REAL(dp), INTENT(OUT) :: normal_0123(0:3), skew_0123(0:3) ! n/l;     
    REAL(dp) :: normal(0:maxmul), skew  (0:maxmul), &
                f_errors(0:50), field(2,0:maxmul)
    INTEGER :: n_norm, n_skew, n_ferr ! number of terms in command line
    INTEGER :: node_fd_errors ! function
    integer :: i_count, n_dim_mult_err

    !initialization
    normal_0123(:)=zero
    skew_0123(:)=zero

           ! real(dp) f_errors(0:50),normal(0:maxmul),skew(0:maxmul)
       ! Get multipole components on bench !-----------------------!
       call dzero(normal,maxmul+1) ! make zero "normal"            !                     
       call dzero(skew,maxmul+1)   ! make zero "skew"              !
                                                                   !
       ! madxdict.h: "knl = [r, {0}], "                            !
       !             "ksl = [r, {0}], "                            !
       ! Assign values from the command line                       !
       call get_node_vector('knl ',n_norm,normal)                  !    
       call get_node_vector('ksl ',n_skew,skew)                    !
       ! void get_node_vector(char*par,int*length,double* vector)  !
       ! /* returns vector for parameter par of current element */ !
                                                                   !
       ! get errors                                                !
       call dzero(f_errors,maxferr+1)                              !
       n_ferr = node_fd_errors(f_errors) !                         !
                ! /* returns the field errors of a node */         !
       call dzero(field,2*(maxmul+1)) ! array to be zeroed.        !
       if (n_ferr .gt. 0) then                                     !
          call dcopy(f_errors,field,n_ferr)                        !
            ! subroutine dcopy(in,out,n)                           !
            ! Purpose:   Copy arrays.                              !
       endif                                                       !
       !-----------------------------------------------------------!
 
       ! fill strength of ALL normal multipoles
       if(n_norm.gt.0) then  ! ========================!
          do i_count=0,n_norm                          !
             if(i_count.gt.0) &                        !
             key%list%k(i_count+1)=normal(i_count)/l   !
             if (i_count.le.3) &                       !
               normal_0123(i_count)=normal(i_count)/l  !
          enddo                                        !
       endif !=========================================!

       ! fill strength of ALL skew multipoles
       if(n_skew.gt.0) then !==========================! 
          do i_count=0,n_skew                          !
             if(i_count.gt.0) &                        !
             key%list%ks(i_count+1)=skew(i_count)/l    !
             if (i_count.le.3) &                       !
               skew_0123(i_count)=skew(i_count)/l      !
          enddo                                        !
       endif !=========================================!   
 
       n_dim_mult_err = max(n_norm, n_skew, n_ferr/2) !========!
       if(n_dim_mult_err.ge.maxmul) n_dim_mult_err=maxmul-1    !
       if(n_ferr.gt.0) then                                    !
          do i_count=0,n_dim_mult_err                          !
            key%list%k(i_count+1)=key%list%k(i_count+1)+ &     !
                                      field(1,i_count)/l       !
            key%list%ks(i_count+1)=key%list%ks(i_count+1)+ &   !
                                      field(2,i_count)/l       !
          enddo                                                !
       endif !=================================================!

  END SUBROUTINE SUMM_MULTIPOLES_AND_ERRORS
  !----------------------------------------------------------------


  subroutine ptc_align()
    implicit none
    include 'twiss0.fi'
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
       call mad_misalign_fibre(f,al_errors(1:6))
    endif
    f=>f%next
    if(advance_node().ne.0)  goto 10

  END subroutine ptc_align
  !_________________________________________________________________

  subroutine ptc_dumpmaps()
    !Dumps to file maps and/or matrixes (i.e. first order maps)
    implicit none
    type(fibre), pointer :: p
    type(damap)          :: id !identity map used for calculating maps for each element
    type(real_8)         :: y2(6)  !polimorphes array used for calculating maps for each element
    real(dp)             :: xt(6)
    integer              :: i, ii  !iterators
    character(200)       :: filename='ptcmaps.txt'
    integer              :: get_string
    real(kind(1d0))      :: get_value
    integer              :: flag_index,why(9)


    if (cavsareset .eqv. .false.) then
       call setcavities(my_ring,maxaccel)
    endif
    !      ii = get_string('ptc_dumpmaps ','file ',filename)
    !      if (ii < 1) then
    !        print *, '<madx_ptc_module.f90 : ptc_dumpmaps> Specified file name has 0 legth',filename
    !        return
    !      endif

    if (getdebug() > 1) print *, '<madx_ptc_module.f90 : ptc_dumpmaps> Maps are dumped to file ',filename
    open(unit=42,file=filename)

    call init(getintstate(),1,c_%np_pol,berz)

    call alloc(id);
    call alloc(y2);

    xt(:) = zero
    id    = 1     ! making identity map

    p=>my_ring%start
    do i=1,my_ring%n

       y2=xt+id ! we track identity map from the current position

       if(p%mag%kind/=kind21) then
          call track(my_ring,y2,i,i+1,getintstate())

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) " ptc_dumpmaps-1 "
             Write(6,*) why ! See produce aperture flag routine in sd_frame
             c_%watch_user=.false.
             return
          endif

          call track(my_ring,xt,i,i+1,getintstate())

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) " ptc_dumpmaps-2 "
             Write(6,*) why ! See produce aperture flag routine in sd_frame
             c_%watch_user=.false.
             return
          endif
       else
          if (getdebug() > 2) print *, 'Track Cavity...'

          call track(my_ring,y2,i,i+2,getintstate())

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) " ptc_dumpmaps-3 "
             Write(6,*) why ! See produce aperture flag routine in sd_frame
             c_%watch_user=.false.
             return
          endif

          call track(my_ring,xt,i,i+2,getintstate())

          call PRODUCE_APERTURE_FLAG(flag_index)
          if(flag_index/=0) then
             call ANALYSE_APERTURE_FLAG(flag_index,why)
             Write(6,*) " ptc_dumpmaps-4 "
             Write(6,*) why ! See produce aperture flag routine in sd_frame
             c_%watch_user=.false.
             return
          endif
       endif


       write(42,*) p%mag%name,' ==========================='
       do ii=1,4
          write(42,'(6f13.8)')  y2(ii).sub.'100000', &
               &                      y2(ii).sub.'010000', &
               &                      y2(ii).sub.'001000', &
               &                      y2(ii).sub.'000100', &
               &                      y2(ii).sub.'000001', & !madx format has dp/p at the last column
               &                      y2(ii).sub.'000010'    !
       enddo
       do ii=6,5,-1
          write(42,'(6f13.8)')  y2(ii).sub.'100000', &
               &                      y2(ii).sub.'010000', &
               &                      y2(ii).sub.'001000', &
               &                      y2(ii).sub.'000100', &
               &                      y2(ii).sub.'000001', & !madx format has dp/p at the last column
               &                      y2(ii).sub.'000010'    !
       enddo

       p=>p%next
    enddo

    close(42)
    call kill(y2);
    call kill(id);

  end subroutine ptc_dumpmaps

  !_________________________________________________________________

  subroutine ptc_twiss(tab_name)
    implicit none
    include 'twissa.fi'
    logical(lp) closed_orbit,beta_flg,betz_flg
    integer k,i,ii,no,mynd2,npara,nda,icase,flag_index,why(9),my_nv,nv_min
    integer inval,ioptfun,iii,restart_sequ,advance_node,get_option
    integer tab_name(*)
    real(dp) x(6),deltap0,deltap,betx,alfx,mux,bety,alfy,muy,betz,alfz,muz,dx,dpx,dy,dpy,d_val
    real(kind(1d0)) get_value,suml
    type(real_8) y(6)
    type(twiss) tw
    type(fibre), POINTER :: current
    type(work)   :: startfen !Fibre energy at the start
    real(dp) r,re(6,6),dt
    logical(lp) initial_matrix_manual, initial_matrix_table
    integer n_vector,order,nx,nxp,ny,nyp,nt,ndeltap
    integer row,double_from_table
    integer  :: charge    ! charge of an accelerated particle


    if (getdebug() > 1) print*,"ptc_twiss"
    !------------------------------------------------------------------------------
    table_name = charconv(tab_name)

    if (getdebug() > 1) print*,"ptc_twiss: Table name is ",table_name

    if(universe.le.0) then
       call fort_warn('return from ptc_twiss: ',' no universe created')
       return
    endif
    if(index.le.0) then
       call fort_warn('return from ptc_twiss: ',' no layout created')
       return
    endif

    call cleartables()

    nda=0
    suml=zero

    icase = get_value('ptc_twiss ','icase ')
    deltap0 = get_value('ptc_twiss ','deltap ')

    deltap = zero
    call my_state(icase,deltap,deltap0)
    CALL UPDATE_STATES
    
    
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
       return
    endif

    if(closed_orbit) then
       call find_orbit(my_ring,x,1,default,c_1d_7)
       CALL write_closed_orbit(icase,x)
    endif

    no = get_value('ptc_twiss ','no ')

    call init(default,no,nda,BERZ,mynd2,npara)
    call alloc(y)
    y=npara
    Y=X

    beta_flg = (get_value('ptc_twiss ','betx ').gt.0) .and. (get_value('ptc_twiss ','bety ').gt.0)
    betz_flg = (get_value('ptc_twiss ','betz ').gt.0) .and. (npara.eq.6)

    initial_matrix_manual = get_value('ptc_twiss ','initial_matrix_manual ') .ne. 0
    initial_matrix_table = get_value('ptc_twiss ','initial_matrix_table ') .ne. 0

    if(initial_matrix_table) then
       k = double_from_table("map_table ", "nv ", 1, doublenum)
       if(k.ne.-1) then
          call liepeek(iia,icoast)
          my_nv=int(doublenum)
          nv_min=min(iia(2),my_nv)
       else
          initial_matrix_table=.false.
       endif
    endif

    if(initial_matrix_table) then
       x(:)=zero
       allocate(j(iia(2)))
       j(:)=0
       do i = 1,my_nv
          k   = double_from_table("map_table ", "coef ", i, doublenum)
          d_val=doublenum
          if(i.le.iia(2)) then
             x(i)  = d_val-(y(i)%T.sub.j)
          endif
       enddo
       do i = 1,nv_min
          do ii = 1,nv_min
             j(ii)  = 1
             row    = i*my_nv+ii
             k   = double_from_table("map_table ", "coef ", row, doublenum)
             d_val=doublenum
             d_val  = d_val-(y(i)%T.sub.j)
             y(i)%T = y(i)%T+(d_val.mono.j)
             j(ii)=0
          enddo
       enddo
       deallocate(j)
    elseif(initial_matrix_manual) then
       call readinitialmatrix()
    elseif(beta_flg) then
       call readinitialtwiss()
    else
       c_%watch_user=.true.
       call track(my_ring,y,1,default)
       call PRODUCE_APERTURE_FLAG(flag_index)
       if(flag_index/=0) then
          call ANALYSE_APERTURE_FLAG(flag_index,why)
          Write(6,*) "ptc_twiss unstable (map production)-programs continues "
          Write(6,*) why ! See produce aperture flag routine in sd_frame
          c_%watch_user=.false.
          CALL kill(y)
          return
       endif
    endif


    if (cavsareset .eqv. .false.) then
       call setcavities(my_ring,maxaccel)
    endif


    !############################################################################
    !############################################################################
    !############################################################################

    call alloc(tw)
    tw=y
    y=npara
    Y=X
    current=>MY_RING%start
    startfen = 0
    startfen = current!setting up start energy for record
    suml=zero
    iii=restart_sequ()
    print77=.false.
    open(unit=21,file='ptctwiss.txt')

    if (getdebug() > 2) then
       print *, "ptc_twiss: internal state is:"
       call print(default,6)
    endif

    do i=1,MY_RING%n
       
       if (getdebug() > 1) then
          write(6,*) "##########################################"
          write(6,'(i4, 1x,a, f10.6)') i,current%mag%name, suml
          write(6,'(a, f9.6, a)') "Ref Momentum ",current%mag%p%p0c," GeV/c"
       endif

       call track(my_ring,y,i,i+1,default)

       call PRODUCE_APERTURE_FLAG(flag_index)
       if(flag_index/=0) then
          call ANALYSE_APERTURE_FLAG(flag_index,why)
          Write(6,*) "ptc_twiss unstable (Twiss parameters) element: ",i," name: ",current%MAG%name,"-programs continues "
          Write(6,*) why ! See produce aperture flag routine in sd_frame
          goto 100
       endif

       call putusertable(i,current%mag%name,y)

       suml=suml+current%MAG%P%ld
       tw=y
       call puttwisstable()

       iii=advance_node()
       current=>current%next
    enddo
100 continue
    c_%watch_user=.false.
    call kill(tw)
    CALL kill(y)
    call f90flush(20,.false.)

    close(21)

    !****************************************************************************************
    !*********  E N D   O F   PTC_TWISS      ************************************************
    !****************************************************************************************
    !________________________________________________________________________________________

  contains  ! what follows are internal subroutines of ptc_twiss
    !____________________________________________________________________________________________

    subroutine puttwisstable()
      implicit none
      include 'twissa.fi'
      integer i1,i2,ii,i1a,i2a
      real(kind(1d0))   :: opt_fun(72)
      real(dp)   :: deltae
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

      deltae = cfen%energy/startfen%energy

      write(21,*) "##########################################"
      write(21,*) ""
      write(21,'(i4, 1x,a, f10.6)') i,current%mag%name,suml
      write(21,'(3(a, f10.6))') "Ref Momentum ",cfen%p0c," Energy ", cfen%energy," DeltaE ",deltae
      write(21,*) ""
      call print(y,21)

      call double_to_table(table_name, 's ', suml)
      doublenum=current%mag%p%p0c
      call double_to_table(table_name, 'energy ', doublenum)

      opt_fun(:)=zero
      call liepeek(iia,icoast)
      allocate(j(iia(2)))
      j(:)=0
      do ii=1,iia(2) !iia(2)==nv
         opt_fun(ii)=y(ii)%T.sub.j
      enddo
      deallocate(j)

      ioptfun=6
      call vector_to_table(table_name, 'x ', ioptfun, opt_fun(1))

      opt_fun(1)=tw%beta(1,1) * deltae
      opt_fun(2)=tw%beta(1,2) * deltae
      opt_fun(3)=tw%beta(1,3) * deltae
      opt_fun(4)=tw%beta(2,1) * deltae
      opt_fun(5)=tw%beta(2,2) * deltae
      opt_fun(6)=tw%beta(2,3) * deltae
      opt_fun(7)=tw%beta(3,1) * deltae
      opt_fun(8)=tw%beta(3,2) * deltae
      opt_fun(9)=tw%beta(3,3) * deltae
      opt_fun(10)=tw%alfa(1,1) * deltae
      opt_fun(11)=tw%alfa(1,2) * deltae
      opt_fun(12)=tw%alfa(1,3) * deltae
      opt_fun(13)=tw%alfa(2,1) * deltae
      opt_fun(14)=tw%alfa(2,2) * deltae
      opt_fun(15)=tw%alfa(2,3) * deltae
      opt_fun(16)=tw%alfa(3,1) * deltae
      opt_fun(17)=tw%alfa(3,2) * deltae
      opt_fun(18)=tw%alfa(3,3) * deltae
      opt_fun(19)=tw%gama(1,1) * deltae
      opt_fun(20)=tw%gama(1,2) * deltae
      opt_fun(21)=tw%gama(1,3) * deltae
      opt_fun(22)=tw%gama(2,1) * deltae
      opt_fun(23)=tw%gama(2,2) * deltae
      opt_fun(24)=tw%gama(2,3) * deltae
      opt_fun(25)=tw%gama(3,1) * deltae
      opt_fun(26)=tw%gama(3,2) * deltae
      opt_fun(27)=tw%gama(3,3) * deltae
      opt_fun(28)=tw%mu(1) * deltae
      opt_fun(29)=tw%mu(2) * deltae
      opt_fun(30)=tw%mu(3) * deltae
      opt_fun(31)=tw%disp(1) * deltae
      opt_fun(32)=tw%disp(2) * deltae
      opt_fun(33)=tw%disp(3) * deltae
      opt_fun(34)=tw%disp(4) * deltae
      opt_fun(35)=tw%disp(5) * deltae
      opt_fun(36)=tw%disp(6) * deltae
      do i1=1,nd2
         if(i1.le.4) then
            i1a=i1
         elseif(i1.eq.5) then
            i1a=6
         else
            i1a=5
         endif
         do i2=1,nd2
            if(i2.le.4) then
               i2a=i2
            elseif(i2.eq.5) then
               i2a=6
            else
               i2a=5
            endif
            ii=37+(i1a-1)*6+i2a
            opt_fun(ii)=tw%eigen(i1,i2) * deltae
            if(mytime.and.i2a.eq.6) opt_fun(ii)=-opt_fun(ii)
         enddo
      enddo

      if (getdebug() > 2)  then
        write(6,'(a16,4f12.3)') 'b11,b12,b21,b22: ',&
             &opt_fun(1),opt_fun(2),opt_fun(4),opt_fun(5)
      endif

      ioptfun=72
      call vector_to_table(table_name, 'beta11 ', ioptfun, opt_fun(1))
      call augment_count(table_name)
      write(20,'(a,13(f9.6))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%mu(3),tw%beta(1,1),tw%beta(1,2),&
           tw%beta(2,1),tw%beta(2,2),tw%beta(3,1),tw%disp(1),tw%disp(3)
      !write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%mu(3),tw%beta(1,1),&
      !     tw%beta(2,1),tw%beta(2,2),&

    end subroutine puttwisstable
    !____________________________________________________________________________________________

    subroutine readinitialmatrix
      !reads initial map elements from MAD-X ptc_twiss command parameters
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
      x(:)=zero
      x(1)=get_value('ptc_twiss ','x ')
      x(2)=get_value('ptc_twiss ','px ')
      x(3)=get_value('ptc_twiss ','y ')
      x(4)=get_value('ptc_twiss ','py ')
      x(5)=get_value('ptc_twiss ','t ')
      x(6)=get_value('ptc_twiss ','pt ')

      call liepeek(iia,icoast)
      allocate(j(iia(2)))
      j(:)=0
      do i = 1,iia(2)
         do ii = 1,iia(2)
            j(ii)=1
            r=re(i,ii)-(y(i)%T.sub.j)
            y(i)%T=y(i)%T+(r.mono.j)
            j(ii)=0
         enddo
      enddo
      deallocate(j)

    end subroutine readinitialmatrix
    !_________________________________________________________________

    subroutine readinitialtwiss
      !Reads initial twiss parameters from MAD-X command
      implicit none
      integer get_option

      betx = get_value('ptc_twiss ','betx ')
      bety = get_value('ptc_twiss ','bety ')
      alfx = get_value('ptc_twiss ','alfx ')
      alfy = get_value('ptc_twiss ','alfy ')
      dx   = get_value('ptc_twiss ','dx ')
      dpx  = get_value('ptc_twiss ','dpx ')
      dy   = get_value('ptc_twiss ','dy ')
      dpy  = get_value('ptc_twiss ','dpy ')
      mux  = get_value('ptc_twiss ','mux ')
      muy  = get_value('ptc_twiss ','muy ')

      x(:)=zero
      x(1)=get_value('ptc_twiss ','x ')
      x(2)=get_value('ptc_twiss ','px ')
      x(3)=get_value('ptc_twiss ','y ')
      x(4)=get_value('ptc_twiss ','py ')
      x(5)=get_value('ptc_twiss ','t ')
      x(6)=get_value('ptc_twiss ','pt ')


      inval = get_option('twiss_inval ')
      if( (inval.eq.1) .and. (getdebug() > 1)) then
         print*," Read BETA0 block in module ptc_twiss"
      endif


      if(mux.eq.zero.or.fraction(mux).eq.zero) mux=mux_default
      if(muy.eq.zero.or.fraction(muy).eq.zero) muy=muy_default
      re(1,1) = cos(twopi*mux)+alfx*sin(twopi*mux)
      re(1,2) = betx*sin(twopi*mux)
      re(1,3) = zero
      re(1,4) = zero
      re(1,5) = (one-re(1,1))*dx-re(1,2)*dpx
      re(1,6) = zero
      re(2,1) = -(one+alfx**2)*sin(twopi*mux)/betx
      re(2,2) = cos(twopi*mux)-alfx*sin(twopi*mux)
      re(2,3) = zero
      re(2,4) = zero
      re(2,5) = -re(2,1)*dx+(1-re(2,2))*dpx
      re(2,6) = zero
      re(3,1) = zero
      re(3,2) = zero
      re(3,3) = cos(twopi*muy)+alfy*sin(twopi*muy)
      re(3,4) = bety*sin(twopi*muy)
      re(3,5) = (1-re(3,3))*dy-re(3,4)*dpy
      re(3,6) = zero
      re(4,1) = zero
      re(4,2) = zero
      re(4,3) = -(1+alfy**2)*sin(twopi*muy)/bety
      re(4,4) = cos(twopi*muy)-alfy*sin(twopi*muy)
      re(4,5) = -re(4,3)*dy+(1-re(4,4))*dpy
      re(4,6) = zero
      re(5,1) = zero
      re(5,2) = zero
      re(5,3) = zero
      re(5,4) = zero
      re(5,5) = one
      re(5,6) = zero
      re(6,1) = zero
      re(6,2) = zero
      re(6,3) = zero
      re(6,4) = zero
      re(6,5) = zero
      re(6,6) = one

      if(betz_flg) then
         betz = get_value('ptc_twiss ','betz ')
         alfz = get_value('ptc_twiss ','alfz ')
         muz  = get_value('ptc_twiss ','muz ')
         if(muz.eq.zero.or.fraction(muz).eq.zero) muz=muz_default
         re(1,5) = zero
         re(2,5) = zero
         re(3,5) = zero
         re(4,5) = zero
         re(6,6) = cos(twopi*muz)+alfz*sin(twopi*muz)
         re(6,5) = betz*sin(twopi*muz)
         re(5,6) = -(1+alfz**2)*sin(twopi*muz)/betz
         re(5,5) = cos(twopi*muz)-alfz*sin(twopi*muz)
      endif
      call liepeek(iia,icoast)
      allocate(j(iia(2)))
      j(:)=0
      do i = 1,iia(2)
         do ii = 1,iia(2)
            j(ii)=1
            r=re(i,ii)-(y(i)%T.sub.j)
            y(i)%T=y(i)%T+(r.mono.j)
            j(ii)=0
         enddo
      enddo
      deallocate(j)
      do i=1,iia(2)
         x(i) = 0
      enddo

    end subroutine readinitialtwiss

  END subroutine ptc_twiss
  !_________________________________________________________________


  FUNCTION double_from_ptc(var, column)
    USE ptc_results
    implicit none
    real(dp) double_from_ptc
    character (len = *) var
    integer :: column(*)
    integer j,k,n1,n2,nn,var_length,ptc_variable_length,ind(6)
    integer double_from_table, string_from_table
    logical(lp) var_check

    double_from_ptc = zero
    var_length = LEN(ptc_var)
    do j = 1 , number_variables
       ptc_variable_length = LEN(ptc_variables(j))
       if (var_length .NE. ptc_variable_length) THEN
          print *,"No such a variable"
          RETURN
       ENDIF
       var_check = .false.
       do k = 1, var_length
          if (ptc_variables(j)(k:k) .NE. ptc_var(k:k)) EXIT
          var_check = .true.
       enddo
       if (var_check) EXIT
    enddo
    if (.NOT. var_check) THEN
       print *,"No such a variable"
       RETURN
    ENDIF
    nn = 0
    do k = 1 , number_variables
       ind(k) = column(k)
       nn = nn + ind(k)
    enddo
    if (nn > order) THEN
       print *,"The order is larger than ",order
       RETURN
    ENDIF
    double_from_ptc = y(j)%t.sub.ind
  END FUNCTION double_from_ptc

  FUNCTION double_from_ptc_normal(name_var,row,icase)
    USE ptc_results
    implicit none
    logical(lp) name_l
    integer,intent(IN) ::  row,icase
    real(dp) double_from_ptc_normal, d_val, d_val1, d_val2
    integer idx,ii,i1,i2
    integer j,k,ind(6)
    integer double_from_table
    character(len = 4)  name_var
    character(len = 2)  name_var1
    character(len = 3)  name_var2

    name_l = .false.
    double_from_ptc_normal = zero

    name_var1 = name_var
    SELECT CASE (name_var1)
    CASE ("dx")
       ind(:)=0
       k = double_from_table("normal_results ", "order1 ", row, doublenum)
       ind(5) = int(doublenum)
       if (ind(5) == 0) ind(5) = 1
       ind(6) = 0
       d_val = n%A1%V(1).sub.ind
    CASE ('dy')
       ind(:)=0
       k = double_from_table("normal_results ", "order1 ", row, doublenum)
       ind(5) = int(doublenum)
       if (ind(5) == 0) ind(5) = 1
       ind(6) = 0
       d_val = n%A1%V(3).sub.ind
    CASE ('q1')
       ind(:)=0
       d_val = n%dhdj%V(3).sub.ind
    CASE ('q2')
       ind(:)=0
       d_val = n%dhdj%V(4).sub.ind
    CASE DEFAULT
       name_l = .true.
    END SELECT
    if (name_l) then
       name_l = .false.
       name_var2 = name_var
       SELECT CASE (name_var2)
       CASE ('dpx')
          ind(:)=0
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%A1%V(2).sub.ind
       CASE ('dpy')
          ind(:)=0
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%A1%V(4).sub.ind
       CASE ('dq1')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(:)=0
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%dhdj%V(3).sub.ind
       CASE ('dq2')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(:)=0
          ind(5) = int(doublenum)
          if (ind(5) == 0) ind(5) = 1
          ind(6) = 0
          d_val = n%dhdj%V(4).sub.ind
       CASE DEFAULT
          name_l = .true.
       END SELECT
    endif
    if (name_l) then
       SELECT CASE (name_var)
       CASE ('anhx')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          do j = 1,2
             ind(j) =  int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          do j = 3,4
             ind(j) = int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(5) = int(doublenum)
          ind(6) = 0
          d_val = n%dhdj%V(3).sub.ind
       CASE ('anhy')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          do j = 1,2
             ind(j) = int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          do j = 3,4
             ind(j) = int(doublenum)
          enddo
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(5) = int(doublenum)
          ind(6) = 0
          d_val = n%dhdj%V(4).sub.ind
       CASE ('hamc')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrh%cos%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('hams')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrh%sin%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('hama')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val1 = pbrh%cos%h.sub.ind
          d_val2 = pbrh%sin%h.sub.ind
          double_from_ptc_normal = SQRT(d_val1**2 + d_val2**2)
          RETURN
       CASE ('haml')
          double_from_ptc_normal = zero
          RETURN
       CASE ('gnfc')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrg%cos%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('gnfs')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val = pbrg%sin%h.sub.ind
          double_from_ptc_normal = d_val
          RETURN
       CASE ('gnfa')
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          ind(1) = int(doublenum)
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          ind(2) = int(doublenum)
          k = double_from_table("normal_results ", "order3 ", row, doublenum)
          ind(3) = int(doublenum)
          k = double_from_table("normal_results ", "order4 ", row, doublenum)
          ind(4) = int(doublenum)
          ind(5) = 0
          ind(6) = 0
          d_val1 = pbrg%cos%h.sub.ind
          d_val2 = pbrg%sin%h.sub.ind
          double_from_ptc_normal = SQRT(d_val1**2 + d_val2**2)
          RETURN
       CASE ('gnfu')
          double_from_ptc_normal = zero
          RETURN
       CASE ('eign')
          ii=(icase/2)*2
          k = double_from_table("normal_results ", "order1 ", row, doublenum)
          i1 = int(doublenum) 
          if(i1.gt.ii) call aafail('return from double_from_ptc_normal: ',' wrong # of eigenvectors')
          if(i1.eq.5) i1=6
          if(i1.eq.6) i1=5             
          k = double_from_table("normal_results ", "order2 ", row, doublenum)
          i2 = int(doublenum)
          if(i2.gt.ii) call aafail('return from double_from_ptc_normal: ',' eigenvectors too many components')
          if(i2.eq.5) i2=6
          if(i2.eq.6) i2=5             
          ind(:)=0
          ind(i2)=1
          double_from_ptc_normal = n%A_t%V(i1).sub.ind
          if(mytime.and.i2.eq.6) double_from_ptc_normal = -double_from_ptc_normal
          RETURN
       CASE DEFAULT
          print *,"--Error in the table normal_results-- Unknown input: ",name_var
       END SELECT
    endif
    double_from_ptc_normal = d_val*(factorial(ind(1))*factorial(ind(3))*factorial(ind(5)))
  END FUNCTION double_from_ptc_normal

  RECURSIVE FUNCTION FACTORIAL (N) &
       RESULT (FACTORIAL_RESULT)
    INTEGER :: N, FACTORIAL_RESULT

    IF (N <= 0 ) THEN
       FACTORIAL_RESULT = 1
    ELSE
       FACTORIAL_RESULT = N * FACTORIAL (N-1)
    END IF
  END FUNCTION FACTORIAL



  SUBROUTINE display_table_results()
    implicit none
    integer, external :: select_ptc_idx, string_from_table, double_from_table
    character(len = 4) name_var
    integer row,k
    integer :: ord(3)

    print *,"Variable name  Order 1  order 2  order 3        Value      "
    do row = 1 , select_ptc_idx()
       name_var=" "
       k = string_from_table("normal_results ", "name ", row, name_var)
       k = double_from_table("normal_results ", "order1 ", row, doublenum)
       ord(1) = int(doublenum)
       k = double_from_table("normal_results ", "order2 ", row, doublenum)
       ord(2) = int(doublenum)
       k = double_from_table("normal_results ", "order3 ", row, doublenum)
       ord(3) = int(doublenum)
       k = double_from_table("normal_results ", "value ", row, doublenum)
       WRITE(*,100) name_var,ord(1),ord(2),ord(3),doublenum
    enddo
100 FORMAT(3X,A4,14X,I1,8X,I1,8X,I1,5X,f25.18)
  END SUBROUTINE display_table_results
  !_________________________________________________________________

  SUBROUTINE ptc_normal()
    USE ptc_results
    USE madx_ptc_intstate_module
    implicit none
    logical(lp) closed_orbit,normal,maptable
    integer no,mynd2,npara,nda,icase,flag_index,why(9)
    integer i, ii, iii, j1, jj, ja(6), k, l, starti
    integer,parameter :: i_map_coor=10
    integer n_rows,row,n_haml,n_gnfu,nres,mynres,n1,n2,map_term
    integer,external :: select_ptc_idx, minimum_acceptable_order, &
         string_from_table, double_from_table, result_from_normal
    real(dp) x(6),deltap0,deltap,dt
    !type(real_8) y(6)
    integer :: column(6) = (/1,0,0,0,0,0/)
    integer :: ord(3), indexa(4)
    integer :: row_haml(101)
    integer :: index1(1000,2)
    real(dp) coef
    real(kind(1d0)) get_value,val_ptc,map_coor(i_map_coor)
    character(len = 4) name_var
    !------------------------------------------------------------------------------

    if(universe.le.0) then
       call fort_warn('return from ptc_normal: ',' no universe created')
       return
    endif
    if(index.le.0) then
       call fort_warn('return from ptc_normal: ',' no layout created')
       return
    endif

    nda=0

    icase = get_value('ptc_normal ','icase ')
    deltap0 = get_value('ptc_normal ','deltap ')
    
    deltap = zero
    call my_state(icase,deltap,deltap0)
    CALL UPDATE_STATES

    x(:)=zero
    if(mytime) then
       call Convert_dp_to_dt (deltap, dt)
    else
       dt=deltap
    endif
    if(icase.eq.5) x(5)=dt
    closed_orbit = get_value('ptc_normal ','closed_orbit ') .ne. 0
    if(closed_orbit) then
       call find_orbit(my_ring,x,1,default,c_1d_7)
       CALL write_closed_orbit(icase,x)
    endif

    no = get_value('ptc_normal ','no ')

    call init(default,no,nda,BERZ,mynd2,npara)
    call alloc(y)
    y=npara
    Y=X
    c_%watch_user=.true.
    call track(my_ring,y,1,default)
    call PRODUCE_APERTURE_FLAG(flag_index)
    if(flag_index/=0) then
       call ANALYSE_APERTURE_FLAG(flag_index,why)
       Write(6,*) "ptc_normal unstable (map production)-programs continues "
       Write(6,*) why ! See produce aperture flag routine in sd_frame
       CALL kill(y)
       c_%watch_user=.false.
       return
    endif
    c_%watch_user=.false.
    if (getdebug()>1) call daprint(y,18)

    maptable = get_value('ptc_normal ','maptable ') .ne. 0
    if(maptable) then
       map_term=42
       call  make_map_table(map_term)
       call liepeek(iia,icoast)
       allocate(j(iia(2)))
       ja(:)    = 0
       j(:)     = 0
       do iii=1,iia(2)
          coef = y(iii)%T.sub.j
          map_coor(1)=coef
          map_coor(2)=iii
          map_coor(3)=iia(2)
          map_coor(4)=0
          map_coor(5)=ja(1)
          map_coor(6)=ja(2)
          map_coor(7)=ja(3)
          map_coor(8)=ja(4)
          map_coor(9)=ja(5)
          map_coor(10)=ja(6)
          call vector_to_table("map_table ", 'coef ', i_map_coor, map_coor(1))
          call augment_count("map_table ")
       enddo
       do i = 1,iia(2)
          do ii = 1,iia(2)
             j(ii) = 1
             ja(ii) = j(ii)
             coef = y(i)%T.sub.j
             map_coor(1)=coef
             map_coor(2)=i
             map_coor(3)=iia(2)
             map_coor(4)=no
             map_coor(5)=ja(1)
             map_coor(6)=ja(2)
             map_coor(7)=ja(3)
             map_coor(8)=ja(4)
             map_coor(9)=ja(5)
             map_coor(10)=ja(6)
             call vector_to_table("map_table ", 'coef ', i_map_coor, map_coor(1))
             call augment_count("map_table ")
             j(:)  = 0
             ja(ii) = j(ii)
          enddo
       enddo
       deallocate(j)
    endif

    normal = get_value('ptc_normal ','normal ') .ne. 0
    if(normal) then
       call alloc(n)

       !------ Find the number of occurences of the attribute 'haml'

       n_rows = select_ptc_idx()
       n_haml = 0
       n_gnfu = 0
       if (n_rows > 0) then
          do row = 1,n_rows
             name_var=" "
             k = string_from_table("normal_results ", "name ", row, name_var)
             if (name_var .eq. 'gnfu') n_gnfu = n_gnfu + 1
             if (name_var .eq. 'haml') then
                n_haml = n_haml + 1
                row_haml(n_haml) = row
             endif
          enddo
          if (n_gnfu > 0) call alloc(pbrg)
          if (n_haml > 0) then
             call alloc(pbrh)
             do j1 =1,n_haml
                row = row_haml(j1)
                k = double_from_table("normal_results ", "value ", row, doublenum)
                mynres = int(doublenum)
                row = row_haml(j1) - 3*mynres + 2
                starti = 1
                if (j1 .eq. 1) then
                   k = double_from_table("normal_results ", "order1 ", row, doublenum)
                   indexa(1) = int(doublenum)
                   k = double_from_table("normal_results ", "order2 ", row, doublenum)
                   indexa(2) = int(doublenum)
                   k = double_from_table("normal_results ", "order3 ", row, doublenum)
                   indexa(3) = int(doublenum)
                   k = double_from_table("normal_results ", "order4 ", row, doublenum)
                   indexa(4) = int(doublenum)
                   index1(1,1) = indexa(1) - indexa(2)
                   index1(1,2) = indexa(3) - indexa(4)
                   n%m(1,1)= index1(1,1)
                   n%m(2,1)= index1(1,2)
                   if(ndim2.eq.6) n%m(3,1)= indexa(3)
                   nres = 1
                   starti = 2
                endif
                !============ nres is the number of resonances to be set

                if (mynres .ge. starti) then
                   do i = starti,mynres
                      ii = row + 3*(i-1)
                      k = double_from_table("normal_results ", "order1 ", ii, doublenum)
                      indexa(1) = int(doublenum)
                      k = double_from_table("normal_results ", "order2 ", ii, doublenum)
                      indexa(2) = int(doublenum)
                      k = double_from_table("normal_results ", "order3 ", ii, doublenum)
                      indexa(3) = int(doublenum)
                      k = double_from_table("normal_results ", "order4 ", ii, doublenum)
                      indexa(4) = int(doublenum)
                      n1 = indexa(1) - indexa(2)
                      n2 = indexa(3) - indexa(4)
                      do l = 1,nres
                         if (n1 .eq. index1(l,1) .and. n2 .eq. index1(l,2)) goto 100
                      enddo
                      nres = nres + 1
                      index1(nres,1) = n1
                      index1(nres,2) = n2
                      n%m(1,nres)= n1
                      n%m(2,nres)= n2
                      if(ndim2.eq.6) n%m(3,nres)= indexa(3)
100                   continue
                   enddo
                endif
             enddo
             n%nres = nres
          endif
       endif
       !------------------------------------------------------------------------

       n=y
       if (n_gnfu > 0) pbrg = n%a%pb
       if (n_haml > 0) pbrh = n%normal%pb
       write(19,'(/a/)') 'Dispersion, First and Higher Orders'
       if (getdebug()>1) call daprint(n%A1,19)

       !------ get values and store them in the table 'normal_results' ---------

       n_rows = select_ptc_idx()
       if (no .ne. minimum_acceptable_order()) then
          print *,"The minimum required order is ",minimum_acceptable_order(), &
               "while the 'no' attribute in the command ptc_normal is ",no
       endif
       if (no < minimum_acceptable_order()) then
          print *,"ptc_normal failed. MAD-X continues."
          stop
       endif
       
       if (n_rows > 0) then
          do row = 1,n_rows
             name_var=" "
             k = string_from_table("normal_results ", "name ", row, name_var)
             val_ptc = double_from_ptc_normal(name_var,row,icase)
             if (name_var .ne. 'haml'.and.name_var .ne. 'gnfu')    &
                  call double_to_table_row("normal_results ", "value ", row, val_ptc)
          enddo
       endif

       !------------------------------------------------------------------------

       write(19,'(/a/)') 'Tunes, Chromaticities and Anharmonicities'
       !  call daprint(n%A_t,19)
       !  call daprint(n%A,19)

       if (getdebug()>1) call daprint(n%dhdj,19)

       !       call daprint(pbrh,19)
       if (n_gnfu > 0) call kill(pbrg)
       if (n_haml > 0) call kill(pbrh)
       call kill(n)
    endif
    CALL kill(y)
    call f90flush(18,.false.)
    call f90flush(19,.false.)

  END subroutine ptc_normal
  !________________________________________________________________________________

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

    if(universe.le.0) then
       call fort_warn('return from ptc_track: ',' no universe created')
       return
    endif
    if(index.le.0) then
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

    if(universe.le.0) then
       call fort_warn('return from ptc_end: ',' no universe can be killed')
       return
    endif

    call kill_universe(m_u)
    call kill_tpsa
    do i=1,size(s_b)
       call nul_coef(s_b(i))
    enddo
    deallocate(s_b)
    firsttime_coef=.true.
  end subroutine ptc_end

  !________________________________________________________________________________

  subroutine zerotwiss(s1,i)
    implicit none
    type(twiss), intent(inout)::s1
    integer, intent(in)::i

    if(i==0) then

       call liepeek(iia,icoast)
       NO=iia(1)
       ND=iia(3)
       ND2=iia(3)*2
       NP=iia(2)-nd2
       NDPT=icoast(4)
       NV=iia(2)

       s1%a1=1
       s1%a_t=1
       s1%nf=.false.
       s1%beta(:,:)=zero
       s1%alfa(:,:)=zero
       s1%gama(:,:)=zero
       s1%mu(:)=zero
       s1%disp(:)=zero
       s1%tune(:)=zero
       s1%eigen(:,:)=zero
       dicu(:)=zero
       angp(:,:)=zero
    endif

  end subroutine zerotwiss

  !________________________________________________________________________________

  subroutine equaltwiss(s1,s2)
    implicit none
    type(twiss), intent(inout)::s1
    type(real_8), intent(in)::s2(ndd)
    integer i,j,j1,ii,ii1,ii2,iii,i2,i3,i1,ei2
    integer ind(6)
    real(dp) au(6,6),aui(2),sx,cx,dphi(3)
    character(len=nd2), dimension(:), pointer :: string
   
    icount=icount+1
    if(nd2.eq.4) string=>str4
    if(nd2.eq.6) string=>str6
    if(s2(1)%kind.eq.2) then

       if(.not.s1%nf) then
          s1%nf=.true.
          do i=1,nd2

             s1%junk%v(i)=s2(i)
          enddo
          s1%n=s1%junk
          s1%a1=s1%n%a1
          s1%a_t=s1%n%a_t
          s1%tune(:)=s1%n%tune(:)
          s1%junk=s1%n%a_t
          do i=1,nd2
             dicu(i)=s1%n%a1%v(i).sub.str5(5)
          enddo
       else
          do i=1,nd2
             s1%junk%v(i)=s2(i)
          enddo
          s1%junk=s1%junk*s1%a_t
       endif
    else
       print*,'A Taylor Map is needed as input for the Twiss calculation'
       return
    endif


    if(nd2.eq.4.and.np.ge.1) then
       do i=1,nd2
          do j=1,nd2
             rdd(i,j)=s2(i).sub.str5(j)
          enddo
          rdd(i,5)=s2(i).sub.str5(5)
       enddo
       do j1=1,2
          ii=2*j1

          iii=2*j1-1
          s1%disp(iii)=rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2)+rdd(ii-1,3)*dicu(3)+rdd(ii-1,4)*dicu(4)+rdd(ii-1,5)
          s1%disp(ii)=rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2)+rdd(ii,3)*dicu(3)+rdd(ii,4)*dicu(4)+rdd(ii,5)
       enddo
    endif

    if(nd.eq.3) then
       !       call daprint(s1%junk,6)
       s1%junk1=s1%junk**(-1)
    endif
    ind(:)=0
    do j=1,nd
       ii=2*j
       ii1=ii-1
       ii2=ii
       if(j.eq.1) then
          i2=4
          i3=6
       elseif(j.eq.2) then
          i2=2
          i3=6
       elseif(j.eq.3) then
          i2=2
          i3=4
          ii1=ii
          ii2=ii-1
       endif

       do ei2=1,nd2
          ind(ei2)=1
          s1%eigen(ii1,ei2)=s1%junk%V(ii1).sub.ind
          s1%eigen(ii,ei2)=s1%junk%V(ii).sub.ind
          ind(ei2)=0
       enddo
    
       angp(1,ii-1)=s1%junk%v(ii1).sub.string(ii-1)
       au(ii,ii-1)=s1%junk%v(ii2).sub.string(ii-1)
       angp(1,ii)=s1%junk%v(ii1).sub.string(ii)
       
       au(ii,ii)=s1%junk%v(ii2).sub.string(ii)
       au(i2-1,i2-1)=s1%junk%v(ii1).sub.string(i2-1)
       au(i2,i2-1)=s1%junk%v(ii2).sub.string(i2-1)
       au(i2-1,i2)=s1%junk%v(ii1).sub.string(i2)
       au(i2,i2)=s1%junk%v(ii2).sub.string(i2)

       s1%beta(1,j)=angp(1,ii-1)*angp(1,ii-1)+angp(1,ii)*angp(1,ii)
       s1%beta(2,j)=au(i2-1,i2-1)*au(i2-1,i2-1)+au(i2-1,i2)*au(i2-1,i2)
       s1%alfa(1,j)=-(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))
       s1%alfa(2,j)=-(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2))
       s1%gama(1,j)=au(ii,ii-1)*au(ii,ii-1)+au(ii,ii)*au(ii,ii)
       s1%gama(2,j)=au(i2,i2-1)*au(i2,i2-1)+au(i2,i2)*au(i2,i2)
       if(nd.eq.3) then
          au(i3-1,i3-1)=s1%junk%v(ii1).sub.string(i3-1)
          au(i3,i3-1)=s1%junk%v(ii2).sub.string(i3-1)
          au(i3-1,i3)=s1%junk%v(ii1).sub.string(i3)
          au(i3,i3)=s1%junk%v(ii2).sub.string(i3)
          s1%beta(3,j)=au(i3-1,i3-1)*au(i3-1,i3-1)+au(i3-1,i3)*au(i3-1,i3)
          s1%alfa(3,j)=-(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3))
          s1%gama(3,j)=au(i3,i3-1)*au(i3,i3-1)+au(i3,i3)*au(i3,i3)
          aui(1)=s1%junk1%v(6).sub.string(6)
          aui(2)=s1%junk1%v(5).sub.string(6)
          if(j.lt.3) then
             s1%disp(ii-1)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
             s1%disp(ii)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
          else
             s1%disp(5)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
             s1%disp(6)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
          endif
       endif

       sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
       cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
       if(abs(sx).gt.c_1d_15.or.abs(cx).gt.c_1d_15) then
          dphi(j)=atan2(sx,cx)*twopii
       else
          dphi(j)=zero
       endif
       s1%mu(j)=s1%mu(j)+dphi(j)
    enddo
    if(icount.gt.1) then
       do i=1,ndim2
          angp(2,i)=angp(1,i)
       enddo
    endif

  end subroutine equaltwiss
  !_________________________________________________________________

  subroutine  alloctwiss(s1)
    implicit none
    type (twiss),intent(inout)::s1

    call alloc(s1%a1)
    call alloc(s1%a_t)
    call alloc(s1%n)
    call alloc(s1%junk)
    call alloc(s1%junk1)

    s1=0
  end subroutine alloctwiss
  !_________________________________________________________________

  subroutine  killtwiss(s1)
    implicit none
    type (twiss),intent(inout)::s1

    call kill(s1%a1)
    call kill(s1%a_t)
    call kill(s1%n)
    call kill(s1%junk)
    call kill(s1%junk1)

    s1%nf=.false.
    s1%beta(:,:)=zero
    s1%alfa(:,:)=zero
    s1%gama(:,:)=zero
    s1%mu(:)=zero
    s1%disp(:)=zero
    s1%tune(:)=zero
    s1%eigen(:,:)=zero

  end subroutine killtwiss
  !_________________________________________________________________

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

  subroutine my_state(icase,deltap,deltap0)
    implicit none
    integer icase,i
    real(dp) deltap0,deltap

    default = getintstate()
    
    if (getdebug()>1) then
      print*, "icase=",icase," deltap=",deltap," deltap0=",deltap0
    endif
    
    deltap = zero
    select case(icase)
    CASE(4)
       if (getdebug()>1) print*, "my_state: Enforcing ONLY_4D+NOCAVITY"
       default=default+only_4d+NOCAVITY
       i=4
    CASE(5)
       if (getdebug()>1) print*, "my_state: Enforcing DELTA"
       default=default+delta
       deltap=deltap0
       i=5
    CASE(6)
       i=6
    CASE DEFAULT
       default=default+only_4d+NOCAVITY
       i=4
    END SELECT

    if (i==6) then
      if (icav==0) then
       default=default+only_4d+NOCAVITY
       call fort_warn('return mystate: ',' no cavity - dimensionality reduced 6 -> 4')
       i=4
      else
       default = default - NOCAVITY !enforcing nocavity to false
      endif 
    endif

    call setintstate(default)
    CALL UPDATE_STATES

    if (getdebug()>0) call print(default,6)

    icase = i

  end subroutine my_state

  !______________________________________________________________________

  subroutine f90flush(i,option)
    implicit none
    integer i,ios
    logical(lp) ostat, fexist,option
    character*20 faction,faccess,fform,fwstat
    character*255 fname
    inquire(err=1,iostat=ios,&
         unit=i,opened=ostat,exist=fexist,write=fwstat)
    if (.not.ostat.or..not.fexist.or.fwstat.ne.'YES') return
    inquire(err=2,iostat=ios,&
         unit=i,action=faction,access=faccess,&
         form=fform,name=fname)
    close (unit=i,err=3)
    !     write (*,*) 'Re-opening ',i,' ',faction,faccess,fform,fname
    if (option) then
       open(err=4,iostat=ios,&
            unit=i,action=faction,access=faccess,form=fform,&
            file=fname,status='old',position='append')
    else
       open(err=4,iostat=ios,&
            unit=i,action=faction,access=faccess,form=fform,&
            file=fname,status='old',position='rewind')
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

END MODULE madx_ptc_module
