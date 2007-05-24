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
  use ptc_spin
  !  USE madx_keywords
  USE madx_ptc_setcavs_module
  USE madx_ptc_knobs_module
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
  integer, private, parameter :: mynreso=20
  integer, private, dimension(4) :: iia,icoast
  integer, private :: NO,ND,ND2,NP,NDPT,NV
  real(dp) :: mux_default=c_0_28, muy_default=c_0_31, muz_default=c_1d_3
  integer, private, allocatable :: J(:)
  logical(lp)             :: savemaps=.false.

  type mapbuffer
     type(universal_taylor)  :: unimap(6)
     real(dp)                :: s
     character(nlp+1)        :: name
  end type mapbuffer

  type(mapbuffer), pointer  :: maps(:) !buffered maps from the last twiss
  integer                   :: mapsorder = 0  !order of the buffered maps, if 0 maps no maps buffered
  integer                   :: mapsicase = 0

CONTAINS

  subroutine ptc_create_universe()
    implicit none
    real(kind(1d0)) get_value

    print77=.false.
    read77 =.false.

    if (getdebug()==0) global_verbose = .false.
    if (getdebug()>0) print*,"Now PTC"
    sector_nmul_max = get_value('ptc_create_universe ','sector_nmul_max ')
    sector_nmul = get_value('ptc_create_universe ','sector_nmul ')
    if(sector_nmul_max.lt.sector_nmul) then
       call aafail('sector_nmul_max must be larger than sector_n: ',&
            'check your ptc_create_universe input')
    endif
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
    real(dp) l,l_machine,energy,kin,brho,beta0,p0c,pma,e0f,lrad,charge
    real(dp) f_errors(0:50),aperture(maxnaper),normal(0:maxmul)
    real(dp) patch_ang(3),patch_trans(3)
    real(dp) skew(0:maxmul),field(2,0:maxmul),fieldk(2)
    real(dp) gamma,gamma2,gammatr2,freq,offset_deltap
    real(dp) fint,fintx,div,muonfactor
    real(dp) sk1,sk1s,sk2,sk2s,sk3,sk3s,tilt,dum1,dum2
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
    integer             nst0,nst1,ord_max
    REAL (dp) :: tempdp
    logical(lp):: ptcrbend,truerbend
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
    charge=get_value('beam ','charge ')
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
    if (getdebug() > 1) print*,'  Global max sector_nmul: ',sector_nmul_max0

    sector_nmul0 = sector_nmul
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
    !    madlength = get_option('rbarc ') .eq. 0
    madlength = .false.
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
       print *, '    charge        ',charge
    endif

    my_ring%mass=pma
    ! preliminary setting
    my_ring%charge=1

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

    call node_name(name,length)
    key%list%vorname=name

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

    select case(code)
    case(0,4,25)
       key%magnet="marker"
    case(1,11,20,21)
       key%magnet="drift"
    case(2) ! PTC accepts mults
       if(l.eq.zero) then
          key%magnet="marker"
          goto 100
       endif
       key%magnet="rbend"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       tempdp=sqrt(normal_0123(0)*normal_0123(0)+skew_0123(0)*skew_0123(0))
       key%list%b0=node_value('angle ')+tempdp*l

       !       print*, "RBEND: Angle: ", node_value('angle ')," tempdp ", tempdp, " l ", l
       !       print*, "RBEND: normal: ",normal_0123(0)," skew: ",skew_0123(0)

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
                write(6,*) " CHANGE YOUR LATTICE FILE."
                stop 666
             endif
          else
             key%magnet="WEDGRBEND"
          endif
       endif
    case(3) ! PTC accepts mults watch out sector_nmul defaulted to 4
       if(l.eq.zero) then
          key%magnet="marker"
          goto 100
       endif
       key%magnet="sbend"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)
       if(sector_nmul_max.lt.ord_max.and.EXACT_MODEL) call aafail('the order of multipoles in a sbend in exact mode cannot be ',&
            &'larger than sector_mul_max: check your ptc_create_universe input')

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
       if(tempdp.gt.0) key%tiltd=key%tiltd + atan2(skew_0123(0),normal_0123(0))

    case(5)
       key%magnet="quadrupole"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       ! Read data & fill %k(:), %ks(:) arrays which are
       ! summs of multipoles and errors

       ! quadrupole components
       sk1= node_value('k1 ')
       sk1s=node_value('k1s ')
       tilt=node_value('tilt ')
       dum1=key%list%k(2)-normal_0123(1)
       dum2=key%list%ks(2)-skew_0123(1)

       if(dum1.ne.zero.or.dum2.ne.zero) then                      !
          sk1= sk1 +dum1                                          !
          sk1s=sk1s+dum2                                          !
       endif                                                      !
       if (sk1s .ne. zero) then                                   !
          tilt = -atan2(sk1s, sk1)/two + tilt                     !
          sk1 = sqrt(sk1**2 + sk1s**2)                            !
       endif                                                      !
       key%list%k(2) =sk1                                         !
       key%list%ks(2)=zero  ! added by VK                         !
       key%tiltd=tilt  !==========================================!

       !================================================================

    case(6)
       key%magnet="sextupole"
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

       ! sextupole components
       sk2= node_value('k2 ')
       sk2s=node_value('k2s ')
       tilt=node_value('tilt ')
       dum1=key%list%k(3)-normal_0123(2)
       dum2=key%list%ks(3)-skew_0123(2)

       if(dum1.ne.zero.or.dum2.ne.zero) then                      !
          sk2= sk2 +dum1                                          !
          sk2s=sk2s+dum2                                          !
       endif                                                      !
       if (sk2s .ne. zero) then                                   !
          tilt = -atan2(sk2s, sk2)/three + tilt                   !
          sk2 = sqrt(sk2**2 + sk2s**2)                            !
       endif                                                      !
       key%list%k(3) =sk2                                         !
       key%list%ks(3)=zero  ! added by VK                         !
       key%tiltd=tilt  !==========================================!

       !================================================================

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

       if(dum1.ne.zero.or.dum2.ne.zero) then                      !
          sk3= sk3 +dum1                                          !
          sk3s=sk3s+dum2                                          !
       endif                                                      !
       if (sk3s .ne. zero) then                                   !
          tilt = -atan2(sk3s, sk3)/four + tilt                    !
          sk3 = sqrt(sk3**2 + sk3s**2)                            !
       endif                                                      !
       key%list%k(4) =sk3                                         !
       key%list%ks(4)=zero  ! added by VK                         !

       key%tiltd=tilt  !==========================================!

       !================================================================

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
       if(l.ne.zero) then
          key%list%bsol=node_value('ks ')
       else
          print*,"Thin solenoid: ",name," not implemented in PTC"
          stop
       endif
       !VK
       CALL SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123,skew_0123,ord_max)

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
       key%magnet="twcavity"
       key%list%volt=node_value('volt ')
       freq=c_1d6*node_value('freq ')
       key%list%lag=node_value('lag ')*twopi
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
       print*,"Element: ",name," not implemented in PTC"
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

    if (getdebug() > 1) then
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^    F I N I S H E D      P T C     I N P U T    ^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
       print *, '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
    endif

  END subroutine ptc_input
  !_________________________________________________________________

  SUBROUTINE SUMM_MULTIPOLES_AND_ERRORS (l, key, normal_0123, skew_0123,ord_max)
    implicit none
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
    integer :: i_count, n_dim_mult_err, ord_max

    !initialization
    normal_0123(:)=zero
    skew_0123(:)=zero

    ! real(dp) f_errors(0:50),normal(0:maxmul),skew(0:maxmul)
    ! Get multipole components on bench !-----------------------!
    call dzero(normal,maxmul+1) ! make zero "normal"            !
    call dzero(skew,maxmul+1)   ! make zero "skew"              !
    !                                                           !
    ! madxdict.h: "knl = [r, {0}], "                            !
    !             "ksl = [r, {0}], "                            !
    ! Assign values from the command line                       !
    call get_node_vector('knl ',n_norm,normal)                  !
    call get_node_vector('ksl ',n_skew,skew)                    !
    if(n_norm.ge.maxmul) n_norm=maxmul-1       !
    if(n_skew.ge.maxmul) n_skew=maxmul-1       !
    ord_max=max(n_norm,n_skew)
    ! void get_node_vector(char*par,int*length,double* vector)  !
    ! /* returns vector for parameter par of current element */ !
    !                                                           !
    ! get errors                                                !
    call dzero(f_errors,maxferr+1)                              !
    n_ferr = node_fd_errors(f_errors) !                         !
    ! /* returns the field errors of a node */                  !
    call dzero(field,2*(maxmul+1)) ! array to be zeroed.        !
    if (n_ferr .gt. 0) then                                     !
       call dcopy(f_errors,field,n_ferr)                        !
       ! subroutine dcopy(in,out,n)                             !
       ! Purpose:   Copy arrays.                                !
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
    if(n_ferr.gt.0) then                                       !
       do i_count=0,n_dim_mult_err                             !
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



  END SUBROUTINE SUMM_MULTIPOLES_AND_ERRORS
  !----------------------------------------------------------------

  subroutine ptc_getnfieldcomp(fibreidx, ncomp, nval)
    implicit none
    include 'twissa.fi'
    real(kind(1d0))      :: nval
    integer              :: fibreidx
    integer              :: ncomp
    type(fibre), pointer :: p
    integer              :: j, i
    integer              :: kn, ks
    real(dp)             :: v
    real(kind(1d0))      :: tmpv
    real(kind(1d0)) get_value

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
    include 'twissa.fi'
    real(kind(1d0))      :: nval
    integer              :: fibreidx
    integer              :: ncomp
    type(fibre), pointer :: p
    integer              :: j, i
    integer              :: kn, ks
    real(dp)             :: v
    real(kind(1d0))      :: tmpv
    real(kind(1d0)) get_value

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
    include 'twissa.fi'
    integer              :: fibreidx
    type(fibre), pointer :: p
    integer              :: j, i
    integer              :: kn, ks
    real(dp)             :: v
    real(kind(1d0))      :: tmpv
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

       if (getdebug() > 1) print*,"Setting up KN ", kn, " from ", p%mag%BN(kn) ," to ", v

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

       if (getdebug() > 1) print*,"Setting up KS ", ks, " from ", p%mag%AN(ks) ," to ", v
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

  subroutine extendnmul(f,n)
    implicit none
    type(fibre),  pointer               :: f !fiber
    integer                             :: n !order
    real(dp)    , DIMENSION(:), POINTER :: ANR,BNR !real arrays for regular element
    type(real_8), DIMENSION(:), POINTER :: ANP,BNP !polimorphic arrays for polimorphic element
    integer                             :: i !iterator
    !P.Sk

    if (.not. associated(f)) then
       return
    endif

    ALLOCATE(ANR(n),BNR(n))
    ALLOCATE(ANP(n),BNP(n))
    CALL ALLOC(ANP,n)
    CALL ALLOC(BNP,n)

    DO I=1,f%mag%P%NMUL
       ANR(I)=f%mag%AN(I)
       BNR(I)=f%mag%BN(I)

       ANP(I)=f%magp%AN(I)
       ANP(I)=f%magp%BN(I)
    ENDDO

    DO I=f%mag%P%NMUL+1, n
       ANR(I)=zero
       BNR(I)=zero

       ANP(I)=zero
       ANP(I)=zero
    ENDDO

    call kill(f%magp%AN,f%magp%p%nmul)
    call kill(f%magp%BN,f%magp%p%nmul)
    deallocate(f%magp%AN,f%magp%BN)
    deallocate(f%mag%AN,f%mag%BN)

    f%mag%p%nmul  = n
    f%magp%p%nmul = n

    f%mag%AN=>ANR
    f%mag%BN=>BNR

    f%magp%AN=>ANP
    f%magp%BN=>BNP


  end subroutine extendnmul

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
       write(6,'(6f8.3)')  al_errors(1:6)
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
    type(real_8)         :: yfull(6)  !polimorphes array used for calculating maps for each element
    real(dp)             :: xt(6)
    integer              :: i, ii  !iterators
    character(200)       :: filename='ptcmaps.txt'
    character(200)       :: filenamefull='ptcmaps'
    integer              :: get_string
    real(kind(1d0))      :: get_value
    integer              :: flag_index,why(9)
    character(200)       :: whymsg
    real(kind(1d0))      :: suml=zero
    integer  geterrorflag !C function that returns errorflag value

    suml=zero

    if (cavsareset .eqv. .false.) then
       call setcavities(my_ring,maxaccel)
       if (geterrorflag() /= 0) then
          return
       endif
    endif

    if (getdebug() > 1) print *, '<madx_ptc_module.f90 : ptc_dumpmaps> Maps are dumped to file ',filename
    open(unit=42,file=filename)

    !    write(filenamefull,*) filename,".",my_ring%start%mag%name,"-",my_ring%end%mag%name,".txt"
    filenamefull="ptcmaps.start-end.txt"
    print*, filenamefull
    open(unit=43,file=filenamefull)

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

       if(p%mag%kind/=kind21) then

          call track(my_ring,y2,i,i+1,getintstate())

          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             call fort_warn('ptc_dumpmaps: ','DA got unstable')
             call seterrorflag(10,"ptc_dumpmaps ","DA got unstable ");
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
             return
          endif

          call track(my_ring,xt,i,i+1,getintstate())
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             call fort_warn('ptc_dumpmaps: ','DA got unstable')
             call seterrorflag(10,"ptc_dumpmaps ","DA got unstable ");
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
             return
          endif
       else
          if (getdebug() > 2) print *, 'Track Cavity...'

          call track(my_ring,y2,i,i+2,getintstate())
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             call fort_warn('ptc_dumpmaps: ','DA got unstable')
             call seterrorflag(10,"ptc_dumpmaps ","DA got unstable ");
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
             return
          endif

          call track(my_ring,xt,i,i+2,getintstate())
          if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
             call fort_warn('ptc_dumpmaps: ','DA got unstable')
             call seterrorflag(10,"ptc_dumpmaps ","DA got unstable ");
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
             return
          endif
       endif


       call track(my_ring,yfull,i,i+1,getintstate())

       if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
          call fort_warn('ptc_dumpmaps: ','DA got unstable')
          call seterrorflag(10,"ptc_dumpmaps ","DA got unstable ");
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
          return
       endif

       write(43,*) p%mag%name, suml,' m ==========================='
       call print(yfull,43)

       suml=suml+p%MAG%P%ld

       write(42,*) p%mag%name, suml,' m ==========================='
       if (c_%npara == 6) then
          call dump6dmap(y2, 42)
       elseif (c_%npara == 5) then
          call dump5dmap(y2, 42)
       elseif (c_%npara == 4) then
          call dump4dmap(y2, 42)
       else
          call fort_warn("ptc_dumpmaps","c_%npara is neither 6,5 nor 4")
       endif
       p=>p%next
    enddo

    close(42)
    call kill(y2);
    call kill(id);

    !_________________________________________________________________
    !_________________________________________________________________
    !_________________________________________________________________

  contains
    !_________________________________________________________________
    subroutine dump4dmap(y2, fun)
      implicit none
      type(real_8) :: y2(6)  !polimorphes array used for calculating maps for each element
      integer      :: fun !file unit number
      integer      :: ii

      if (getdebug() > 1) then

      endif

      do ii=1,4
         write(fun,'(4f13.8)')  y2(ii).sub.'1000', &
              &                      y2(ii).sub.'0100', &
              &                      y2(ii).sub.'0010', &
              &                      y2(ii).sub.'0001'
      enddo

    end subroutine dump4dmap
    !_________________________________________________________________

    subroutine dump5dmap(y2, fun)
      implicit none
      type(real_8) :: y2(6)  !polimorphes array used for calculating maps for each element
      integer      :: fun !file unit number
      integer      :: ii
      do ii=1,5
         write(fun,'(5f13.8)')  y2(ii).sub.'10000', &
              &                      y2(ii).sub.'01000', &
              &                      y2(ii).sub.'00100', &
              &                      y2(ii).sub.'00010', &
              &                      y2(ii).sub.'00001'    !
      enddo

    end subroutine dump5dmap
    !_________________________________________________________________

    subroutine dump6dmap(y2, fun)
      implicit none
      type(real_8) :: y2(6)  !polimorphes array used for calculating maps for each element
      integer      :: fun !file unit number
      integer      :: ii

      do ii=1,4
         write(fun,'(6f13.8)')  y2(ii).sub.'100000', &
              &                      y2(ii).sub.'010000', &
              &                      y2(ii).sub.'001000', &
              &                      y2(ii).sub.'000100', &
              &                      y2(ii).sub.'000001', & !madx format has dp/p at the last column
              &                      y2(ii).sub.'000010'    !
      enddo

      do ii=6,5,-1
         write(fun,'(6f13.8)')  y2(ii).sub.'100000', &
              &                      y2(ii).sub.'010000', &
              &                      y2(ii).sub.'001000', &
              &                      y2(ii).sub.'000100', &
              &                      y2(ii).sub.'000001', & !madx format has dp/p at the last column
              &                      y2(ii).sub.'000010'    !
      enddo
    end subroutine dump6dmap


  end subroutine ptc_dumpmaps


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


  RECURSIVE FUNCTION FACTORIAL (N) &
       RESULT (FACTORIAL_RESULT)
    INTEGER :: N, FACTORIAL_RESULT

    IF (N <= 0 ) THEN
       FACTORIAL_RESULT = 1
    ELSE
       FACTORIAL_RESULT = N * FACTORIAL (N-1)
    END IF
  END FUNCTION FACTORIAL

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
       if (( .not. check_stable ) .or. ( .not. c_%stable_da )) then
          call fort_warn('ptc_track: ','DA got unstable')
          call seterrorflag(10,"ptc_track ","DA got unstable ");
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

    if(universe.le.0) then
       call fort_warn('return from ptc_end: ',' no universe can be killed')
       return
    endif

    call killsavedmaps() !module ptc_twiss -> kill buffered maps

    !    call killparresult()
    call resetknobs()  !remove the knobs

    if ( associated(m_u%n) .eqv. .false. ) then
       print*, "We attempt to kill not initialized universe!"
    endif


    call kill_universe(m_u)
    nullify(my_ring)
    call kill_tpsa
    do i=1,size(s_b)
       call nul_coef(s_b(i))
    enddo
    deallocate(s_b)
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
       if (getdebug()>1) print*, "my_state: Enforcing ONLY_4D+NOCAVITY and NO DELTA"
       default=default-delta0
       default=default+only_4d0+NOCAVITY0
       i=4
    CASE(5)
       if (getdebug()>1) print*, "my_state: Enforcing DELTA"
       default=default+delta0
       deltap=deltap0
       i=5
    CASE(56)
       if (getdebug()>1) print*, "my_state: Enforcing coasting beam"
       default = default - delta0 - only_4d0
       default = default + NOCAVITY0
       deltap=deltap0
       i=56
    CASE(6)
       i=6
    CASE DEFAULT
       default=default+only_4d0+NOCAVITY0
       i=4
    END SELECT

    if (i==6) then
       if ( (icav==0) .and. my_ring%closed .and. (getenforce6D() .eqv. .false.)) then
          default = default - delta0 - only_4d0
          default=default +  NOCAVITY0
          call fort_warn('return mystate: ',' no cavity - dimensionality reduced 6 -> 5 and 1/2')
          i=56
       else
          default = default - delta0 - only_4d0 - NOCAVITY0 !enforcing nocavity to false
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
    logical fexist1, ostat1
    character*20 faction,faccess,fform,fwstat
    character*255 fname
    inquire(err=1,iostat=ios,&
         unit=i,opened=ostat1,exist=fexist1,write=fwstat)
    fexist = fexist1
    ostat  = ostat1
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

  subroutine makemaptable(y)
    implicit none
    type(real_8):: y(6)
    integer,parameter :: i_map_coor=10
    integer           :: map_term, ja(6),i,ii,iii
    real(dp)          :: coef
    real(kind(1d0))   :: map_coor(i_map_coor)
    !    type(universal_taylor) :: ut

    map_term=42
    call  make_map_table(map_term)
    call liepeek(iia,icoast)
    allocate(j(c_%npara))
    ja(:)    = 0
    j(:)     = 0
    do iii=1,c_%npara
       coef = y(iii)%T.sub.j
       map_coor(1)=coef
       map_coor(2)=iii
       map_coor(3)=c_%npara
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

    do i = 1,c_%npara

       do ii = 1,c_%npara
          j(ii) = 1
          ja(ii) = j(ii)
          coef = y(i)%T.sub.j
          map_coor(1)=coef
          map_coor(2)=i
          map_coor(3)=c_%npara! 29.06.2006 here was iia(2) - to be verified
          map_coor(4)=sum(ja(:))
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

       !           ut = y(i)
       !           do ii = 1,ut%n
       !              map_coor(1)=ut%c(ii) !coef
       !              map_coor(2)=i !index of taylor
       !              map_coor(3)=c_%npara
       !              map_coor(4)=sum(ut%j(i,:)) !order
       !              map_coor(5)=ut%j(ii,1)
       !              map_coor(6)=ut%j(ii,2)
       !              map_coor(7)=ut%j(ii,3)
       !              map_coor(8)=ut%j(ii,4)
       !              map_coor(9)=ut%j(ii,5)
       !              map_coor(10)=ut%j(ii,6)
       !           enddo


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

END MODULE madx_ptc_module
