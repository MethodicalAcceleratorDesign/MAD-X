MODULE madx_ptc_module
  USE madx_keywords
  implicit none
  integer icav
  integer :: universe=0,index=0,EXCEPTION=0
  integer ipause
  integer,external :: mypause
  real(kind(1d0)) get_value,node_value
  type(layout),pointer :: MY_RING
  type(mad_universe) m_u
  integer, private, parameter :: ndd=ndim2,mynreso=20
  integer, private, dimension(4) :: iia,icoast
  integer, private :: NO,ND,ND2,NP,NDPT,NV,icount=0
  real(dp), private, dimension(ndim2,5) :: rdd
  real(dp), private, dimension(ndim2) :: dicu
  real(dp), private, dimension(2,ndim2) :: angp
  character(len=4), private, dimension(4), target :: str4 = (/'1000','0100','0010','0001'/)
  character(len=5), private, dimension(5), parameter :: str5 = (/'10000','01000','00100','00010','00001'/)
  character(len=6), private, dimension(6), target :: str6 = (/'100000','010000','001000','000100','000001','000010'/)
  private zerotwiss,equaltwiss,alloctwiss,killtwiss

  type twiss
     type(damap) a1
     type(damap) a_t
     type(damap) junk
     type(damap) junk1
     type(normalform) n
     logical nf
     real(dp), dimension(3,3) ::  beta,alfa,gama
     real(dp), dimension(3)   ::  mu
     real(dp), dimension(6)   ::  disp
     real(dp), dimension(3)   ::  tune
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

    print77=.false.
    read77 =.false.
    print*,"Now PTC"

    call set_up_universe(m_u)
    universe=universe+1

  end subroutine ptc_create_universe

  subroutine ptc_create_layout()
    implicit none

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

  end subroutine ptc_create_layout

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

  subroutine ptc_input()
    implicit none
    include 'twtrr.fi'
    include 'name_len.fi'
    include 'twiss0.fi'
    logical(lp) particle,doneit
    integer i,j,k,code,nt,icount,nn,ns,nd
    integer get_option,double_from_table
    integer restart_sequ,advance_node,n_ferr,node_fd_errors
    integer, parameter :: imul=20,nt0=20000,length=16
    real(dp) l,l_machine,energy,kin,brho,beta0,p0c,pma,e0f,lrad
    real(dp) f_errors(0:50),aperture(100),normal(0:maxmul)
    real(dp) skew(0:maxmul),field(2,0:maxmul),fieldk(2)
    real(dp) gamma,gammatr,gamma2,gammatr2,freq,offset_deltap
    real(dp) fint,fintx,div
    real(kind(1d0)) get_value,node_value
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
    !---------------------------------------------------------------
    energy=get_value('beam ','energy ')
    pma=get_value('beam ','mass ')
    e0f=sqrt(ENERGY**2-pma**2)
    print*,'mad-X, energy, k. energy, pma, momentum: ',energy,energy-pma,pma,e0f
    beta0=e0f/ENERGY
    particle=.false.
    if(abs(pma-pmae)/pmae<02) particle=.true.
    !valid October 2002: oldscheme=.false.
    !!valid October 2002: oldscheme=.true.
    CALL MAKE_STATES(PARTICLE)

    !  with_external_frame=.false.
    !  with_internal_frame=.false.
    !  with_chart=.false.
    !  with_patch=.false.

    ! Global Keywords
    sector_nmul_max0 = get_value('ptc_create_layout ','sector_nmul_max ')
    sector_nmul0 = get_value('ptc_create_layout ','sector_nmul ')
    model = get_value('ptc_create_layout ','model ')
    select case(model)
    CASE(1)
       keymod0 = "DRIFT_KICK       "
    CASE(2)
       keymod0 = "MATRIX_KICK      "
    CASE(3)
       keymod0 = "DELTA_MATRIX_KICK"
    CASE DEFAULT
       EXCEPTION=1
       ipause=mypause(444)
       RETURN
    END SELECT
    method0   = get_value('ptc_create_layout ','method ')
    exact0    = get_value('ptc_create_layout ','exact ') .ne. 0
    nst0      = get_value('ptc_create_layout ','nst ')
    ! MAD-X specials
    madlength = get_option('rbarc ') .eq. 0
    mad       = get_value('ptc_create_layout ','mad_mult ') .ne. 0
    mad8      = get_value('ptc_create_layout ','mad8 ') .ne. 0
    gamma     = get_value('beam ','gamma ')
    k         = double_from_table('summ ','gammatr ',1,gammatr)
    gamma2    = gamma**2
    gammatr2  = gammatr**2

    print*,'global max sector_nmul: ',sector_nmul_max0
    print*,'global sector_nmul: ',sector_nmul0
    print*,'global model: ',keymod0
    print*,'global method: ',method0
    print*,'global exact: ',exact0
    print*,'global nst: ',nst0
    print*,'global rbend_length: ',madlength
    print*,'global mad_mult as in mad8: ',mad
    print*,'rbend as in mad8 (only global): ',mad8
    print*,'gamma: ',gamma
    print*,'gammatr: ',gammatr

    !  call Set_Up(MY_RING)

    CALL SET_MADx(energy=energy,METHOD=method0,STEP=nst0)

    icav=0
    nt=0

    j=restart_sequ()
    j=0
    l_machine=zero
10  continue
    call zero_key(key)
    j=j+1
    nt=nt+1
    if(nt==nt0) then
       print*,"More than the maximum number: ",nt0," of elements in the structure==> Program stops"
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
    nst1=node_value("nst ")
    if(nst1.gt.0) then
       nstd = nst1
    else
       nstd = nst0
    endif

    !special node keys
    key%list%permfringe=node_value("permfringe ") .ne. zero
    key%list%kill_ent_fringe=node_value("kill_ent_fringe ") .ne. zero
    key%list%kill_exi_fringe=node_value("kill_exi_fringe ") .ne. zero
    key%list%bend_fringe=node_value("bend_fringe ") .ne. zero

    nn=name_len
    call node_string('apertype ',aptype,nn)
    call get_node_vector('aperture ',nn,aperture)
    if(.not.(aptype.eq."circle".and.aperture(1).eq.zero)) then
       select case(aptype)
       case("circle")
          key%list%aperture_kind=1
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(1)
       case("ellipse") 
          key%list%aperture_kind=1
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(2)
       case("rectangular")
          key%list%aperture_kind=2
          key%list%aperture_x=aperture(1)
          key%list%aperture_y=aperture(2)
       case("lhcscreen")
          key%list%aperture_kind=3
          key%list%aperture_x=aperture(3)
          key%list%aperture_y=aperture(3)
          key%list%aperture_r(1)=aperture(1)
          key%list%aperture_r(2)=aperture(2)
       case("marguerite")
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
       key%list%b0=node_value('angle ')
       !     key%list%k(1)=node_value('k0 ') 
       key%list%k(2)=node_value('k1 ') 
       key%list%k(3)=node_value('k2 ')
       !     key%list%k(3)=node_value('k2 ') 
       !     key%list%ks(1)=node_value('k0s ') 
       !     key%list%ks(2)=node_value('k1s ') 
       !     key%list%ks(3)=node_value('k2s ') 
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
    case(3) ! PTC accepts mults watch out sector_nmul defaulted to 4
       if(l.eq.zero) then
          key%magnet="marker"
          goto 100
       endif
       key%magnet="sbend"
       key%list%b0=node_value('angle ')
       key%list%k(2)=node_value('k1 ')
       key%list%k(3)=node_value('k2 ')
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
    case(5) ! PTC accepts mults
       key%magnet="quadrupole"
       key%list%k(2)=node_value('k1 ')  ! other mul
       key%tiltd=node_value('tilt ')
    case(6)
       key%magnet="sextupole"
       key%list%k(3)=node_value('k2 ')
       key%tiltd=node_value('tilt ')
    case(7) ! PTC accepts mults
       key%magnet="octupole"
       key%list%k(4)=node_value('k3 ')
       key%tiltd=node_value('tilt ')
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
    case(10)
       key%magnet="rfcavity"
       key%list%volt=node_value('volt ')
       freq=c_1d6*node_value('freq ')
       key%list%lag=node_value('lag ')*twopi
       offset_deltap=get_value('ptc_create_layout ','offset_deltap ')
       if(offset_deltap.ne.zero) then 
          default=default+totalpath0
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
       if(offset_deltap.ne.zero) then 
          default=default+totalpath0
          freq=freq*((gammatr2-gamma2)*offset_deltap/gammatr2/gamma2+one)
       endif
       key%list%freq0=freq
       key%list%dphas=node_value("delta_lag ")
       key%list%psi=node_value("psi ")
       key%list%harmon=one
       if(key%list%volt.ne.zero.and.key%list%freq0.ne.zero) icav=1
    case default
       print*,"Element: ",name," not implemented"
       stop
    end select
100 continue
    call create_fibre(my_ring%end,key,EXCEPTION)
    if(advance_node().ne.0)  goto 10

    print*,' Length of machine: ',l_machine
    CALL GET_ENERGY(ENERGY,kin,BRHO,beta0,P0C)

    MY_RING%closed=.true.
    doneit=.true.
    call ring_l(my_ring,doneit)

    call survey(my_ring)

  END subroutine ptc_input

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

  subroutine ptc_twiss(tab_name)
    implicit none
    include 'twissa.fi'
    logical(lp) closed_orbit,initial_matrix
    integer i,ii,no,mynd2,npara,mynpa,nda,icase,flag_index,why(9)
    integer ioptfun,iii,restart_sequ,advance_node
    integer tab_name(*)
    real(dp) opt_fun(36)
    real(dp) x(6),suml,deltap0,deltap
    real(kind(1d0)) get_value
    type(real_8) y(6)
    type(twiss) tw
    type(fibre), POINTER :: current
    integer, allocatable :: J(:) 
    real(dp) r,re(6,6)
    !------------------------------------------------------------------------------
    table_name = charconv(tab_name)

    if(universe.le.0) then
       call fort_warn('return from ptc_twiss: ',' no universe created')
       return
    endif
    if(index.le.0) then
       call fort_warn('return from ptc_twiss: ',' no layout created')
       return
    endif

    nda=0
    suml=zero

    icase = get_value('ptc_twiss ','icase ')
    deltap0 = get_value('ptc_twiss ','deltap ')
    call my_state(icase,deltap,deltap0,mynpa)

    CALL UPDATE_STATES
!    call print(default,6)

    x(:)=zero
    if(icase.eq.5) x(5)=deltap
    closed_orbit = get_value('ptc_twiss ','closed_orbit ') .ne. 0
    if(closed_orbit) then
       call find_orbit(my_ring,x,1,default,1d-7)
       print*,"Closed orbit: ",x
    endif

    no = get_value('ptc_twiss ','no ')

    call init(default,no,nda,BERZ,mynd2,npara)
    call alloc(y)
    y=npara
    Y=X

    initial_matrix = get_value('ptc_twiss ','initial_matrix ') .ne. 0
    if(initial_matrix) then
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
    call alloc(tw)
    tw=y
    y=npara
    Y=X
    current=>MY_RING%start
    suml=zero
    iii=restart_sequ()
    do i=1,MY_RING%n
       call track(my_ring,y,i,i+1,default)
       call PRODUCE_APERTURE_FLAG(flag_index)
       if(flag_index/=0) then
          call ANALYSE_APERTURE_FLAG(flag_index,why)
          Write(6,*) "ptc_twiss unstable (Twiss parameters) element: ",i," name: ",current%MAG%name,"-programs continues "
          Write(6,*) why ! See produce aperture flag routine in sd_frame
          goto 100
       endif
       suml=suml+current%MAG%P%ld
       tw=y

       call double_to_table(table_name, 's ', suml)

       opt_fun(1)=tw%beta(1,1)
       opt_fun(2)=tw%beta(1,2)
       opt_fun(3)=tw%beta(1,3)
       opt_fun(4)=tw%beta(2,1)
       opt_fun(5)=tw%beta(2,2)
       opt_fun(6)=tw%beta(2,3)
       opt_fun(7)=tw%beta(3,1)
       opt_fun(8)=tw%beta(3,2)
       opt_fun(9)=tw%beta(3,3)
       opt_fun(10)=tw%alfa(1,1)
       opt_fun(11)=tw%alfa(1,2)
       opt_fun(12)=tw%alfa(1,3)
       opt_fun(13)=tw%alfa(2,1)
       opt_fun(14)=tw%alfa(2,2)
       opt_fun(15)=tw%alfa(2,3)
       opt_fun(16)=tw%alfa(3,1)
       opt_fun(17)=tw%alfa(3,2)
       opt_fun(18)=tw%alfa(3,3)
       opt_fun(19)=tw%gama(1,1)
       opt_fun(20)=tw%gama(1,2)
       opt_fun(21)=tw%gama(1,3)
       opt_fun(22)=tw%gama(2,1)
       opt_fun(23)=tw%gama(2,2)
       opt_fun(24)=tw%gama(2,3)
       opt_fun(25)=tw%gama(3,1)
       opt_fun(26)=tw%gama(3,2)
       opt_fun(27)=tw%gama(3,3)
       opt_fun(28)=tw%mu(1)
       opt_fun(29)=tw%mu(2)
       opt_fun(30)=tw%mu(3)
       opt_fun(31)=tw%disp(1)
       opt_fun(32)=tw%disp(2)
       opt_fun(33)=tw%disp(3)
       opt_fun(34)=tw%disp(4)
       opt_fun(35)=tw%disp(5)
       opt_fun(36)=tw%disp(6)

       ioptfun=36
       call vector_to_table(table_name, 'beta11 ', ioptfun, opt_fun(1))
       call augment_count(table_name)
       write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%mu(3),tw%beta(1,1),tw%beta(1,2),&
            tw%beta(2,1),tw%beta(2,2),tw%beta(3,1),tw%disp(1),tw%disp(3)
       !write(20,'(a,13(1x,1p,e21.14))') current%MAG%name,suml,tw%mu(1),tw%mu(2),tw%mu(3),tw%beta(1,1),&
       !     tw%beta(2,1),tw%beta(2,2),&
       iii=advance_node()
       current=>current%next
    enddo
100 continue
    c_%watch_user=.false.
    call kill(tw)
    CALL kill(y)
    call f90flush(20,.false.)

  END subroutine ptc_twiss

  subroutine ptc_normal()
    implicit none
    logical(lp) closed_orbit,normal
    integer no,mynd2,npara,mynpa,nda,icase,flag_index,why(9)
    real(dp) x(6),deltap0,deltap
    type(real_8) y(6)
    real(kind(1d0)) get_value
    type(normalform) n
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
    call my_state(icase,deltap,deltap0,mynpa)

    x(:)=zero
    if(icase.eq.5) x(5)=deltap
    closed_orbit = get_value('ptc_normal ','closed_orbit ') .ne. 0
    if(closed_orbit) then
       call find_orbit(my_ring,x,1,default,1d-7)
       print*,"Closed orbit: ",x
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
    call daprint(y,18)
    normal = get_value('ptc_normal ','normal ') .ne. 0
    if(normal) then
       call alloc(n)
       n=y
       write(19,'(/a/)') 'Disperion, First and Higher Orders'
       call daprint(n%A1,19)
       write(19,'(/a/)') 'Tunes, Chromaticities and Anharmonicities'
       !  call daprint(n%A_t,19)
       !  call daprint(n%A,19)
       call daprint(n%dhdj,19)
       call kill(n)
    endif
    CALL kill(y)
    call f90flush(18,.false.)
    call f90flush(19,.false.)

  END subroutine ptc_normal

  subroutine ptc_track()
    implicit none
    integer i,nint,ndble,nchar,int_arr(1),char_l,mynpa,icase,turns,flag_index,why(9)
    real(dp) x0(6),x(6),deltap0,deltap
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
    call my_state(icase,deltap,deltap0,mynpa)

    CALL UPDATE_STATES
!    call print(default,6)

    x0(:)=zero
    if(icase.eq.5) x0(5)=deltap
    closed_orbit = get_value('ptc_track ','closed_orbit ') .ne. 0
    if(closed_orbit) then
       call find_orbit(my_ring,x0,1,default,1d-7)
       print*,"Closed orbit: ",x0
    endif

    call comm_para('coord ',nint,ndble,nchar,int_arr,x,char_a,char_l)

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
       dicu(:)=zero
       angp(:,:)=zero
    endif

  end subroutine zerotwiss

  subroutine equaltwiss(s1,s2)
    implicit none
    type(twiss), intent(inout)::s1
    type(real_8), intent(in)::s2(ndd)
    integer i,j,j1,ii,ii1,ii2,iii,i2,i3
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
       s1%junk1=s1%junk**(-1)
    endif
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

  end subroutine killtwiss

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
                print*,current%MAG%P%NMUL,current%MAG%KIND
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

  subroutine my_state(icase,deltap,deltap0,mynpa)
    implicit none
    integer icase,mynpa
    real(dp) deltap0,deltap

    deltap = zero
    select case(icase)
    CASE(4)
       default=default+only_4d+NOCAVITY
       mynpa=4
    CASE(5)
       default=default+delta
       deltap=deltap0
       mynpa=5
    CASE(6)
       mynpa=6
    CASE DEFAULT
       default=default+only_4d+NOCAVITY
       mynpa=4
    END SELECT
    if(mynpa==6.and.icav==0) then
       default=default+delta
       mynpa=5
    endif

    default=default+time

    CALL UPDATE_STATES
    call print(default,6)

  end subroutine my_state

  subroutine f90flush(i,option)
    implicit none
    integer i,ios
    logical ostat, fexist,option
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

END MODULE madx_ptc_module
