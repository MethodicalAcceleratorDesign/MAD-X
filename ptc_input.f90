subroutine ptc_input(lhc,icav,EXCEPTION)
  use madx_keywords
  implicit none
  include 'twtrr.fi'
  logical(lp) particle,doneit,exact
  integer ipause, mypause
  integer i,j,code,nt,icount,icav,nn,ns,nd,kindi,EXCEPTION
  integer restart_sequ,advance_node,method,nst,n_ferr,node_fd_errors
  integer, parameter :: imul=20,nt0=20000,length=16
  real(dp) l,l_machine,energy,kin,brho,beta0,p0c,pma,e0f,lrad
  real(dp) f_errors(0:50),aperture(100),normal(0:maxmul),skew(0:maxmul),field(2,0:maxmul)
  real(kind(1d0)) get_value,node_value
  character(length) name
  character(24) aptype
  character(20) keymod
  type(layout) lhc
  type(keywords) key
  !---------------------------------------------------------------
  energy=get_value('beam ','energy ')
  pma=get_value('beam ','mass ')
  e0f=sqrt(ENERGY**2-pma**2)
  print*,'mad-X, energy, k. energy, pma, momentum: ',energy,energy-pma,pma,e0f
  beta0=e0f/ENERGY
  particle=.false.
  if(abs(pma-pmae)/pmae<02) particle=.true.
  oldscheme=.false.
  !oldscheme=.true.
  CALL MAKE_STATES(PARTICLE)

  with_external_frame=.false.
  with_internal_frame=.false.
  with_chart=.false.
  with_patch=.false.
  mad=.true. ! permanent
  madlength=.false.

  method = get_value('ptc ','method ')

  kindi = get_value('ptc ','kindi ')
  select case(kindi)
  CASE(1)
     keymod="DRIFT_KICK       "
  CASE(2)
     keymod="MATRIX_KICK      "
  CASE(3)
     keymod="DELTA_MATRIX_KICK"
     method=2
  CASE(4)
     keymod="BEND_KICK        "
  CASE DEFAULT
     EXCEPTION=1
     ipause=mypause(444)
     RETURN
  END SELECT

  exact = get_value('ptc ','exact ') .ne. 0
  if(exact.and.kindi.ne.3) EXACT_MODEL = .true.

  print*,'method: ',method
  nst = get_value('ptc ','nst ')
  print*,'nst: ',nst
  metd = method
  nstd = nst

  call Set_Up(LHC)
 
  CALL SET_MAD(energy=energy,METHOD=method,STEP=nst)

  icav=0
  nt=0

  j=restart_sequ()
  j=0
  l_machine=zero
10 continue
  call zero_key(key)
  key%mad8=.true.
  
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
  key%model=keymod

  nn=24
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
        print*,"Not implemented"
        stop
     end select
  endif
  call append_empty(lhc)
  select case(code)
  case(0,4,25)
     key%magnet="marker"
  case(1,11,20,21)
     key%magnet="drift"
  case(2) ! PTC accepts mults
     key%magnet="rbend"
     key%list%b0=node_value('angle ')
!     key%list%k(1)=node_value('k0 ') 
     key%list%k(2)=node_value('k1 ') 
!     key%list%k(3)=node_value('k2 ') 
!     key%list%ks(1)=node_value('k0s ') 
!     key%list%ks(2)=node_value('k1s ') 
!     key%list%ks(3)=node_value('k2s ') 
! Gymnastic needed since PTC expects MAD8 convention
     key%list%t1=node_value('e1 ')-node_value('angle ')/two
     key%list%t2=node_value('e2 ')-node_value('angle ')/two
     key%list%hgap=node_value('hgap ')
     key%list%fint=node_value('fint ')
     key%list%h1=node_value('h1 ')
     key%list%h2=node_value('h2 ')
     key%tiltd=node_value('tilt ')
  case(3) ! PTC accepts mults watch out sector_nmul defaulted to 4
     key%magnet="sbend"
     key%list%b0=node_value('angle ')
     key%list%k(2)=node_value('k1 ')
     key%list%t1=node_value('e1 ')
     key%list%t2=node_value('e2 ')
     key%list%hgap=node_value('hgap ')
     key%list%fint=node_value('fint ')
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
     if(abs(node_value('lag ')).eq.half) then
        key%list%freq0=c_1d6*node_value('freq ')
        key%list%lag=-(node_value('lag ')+half)
     elseif(abs(node_value('lag ')).eq.zero) then
        key%list%freq0=-c_1d6*node_value('freq ')
        key%list%lag=-(node_value('lag ')+half)
     endif
! Gymnastic needed since PTC expects MAD8 convention
     if(node_value('harmon ').ne.zero) key%list%freq0=key%list%freq0/node_value('harmon ')
!  case(11)
!     key%magnet="elseparator"
!     key%list%volt=node_value('ex ')
!     key%list%lag=atan2(node_value('ey '),node_value('ex '))
!     key%tiltd=node_value('tilt ')
  case(14,15,16) ! PTC accepts mults
     if(code.eq.14) then
        key%magnet="hkicker"
        key%list%k(1)=node_value('kick ')
     else if(code.eq.15) then
        key%magnet="kicker"
        key%list%k(1)=node_value('hkick ')
        key%list%ks(1)=node_value('vkick ')
     else if(code.eq.16) then
        key%magnet="vkicker"
        key%list%ks(1)=node_value('kick ')
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
  case default
     print*,"Not implemented"
     stop
  end select
  call create_fibre(lhc%end,key,EXCEPTION)
  if(advance_node().ne.0)  goto 10

  print*,' Length of machine: ',l_machine
  CALL GET_ENERGY(ENERGY,kin,BRHO,beta0,P0C)

  LHC%closed=.true.
  doneit=.true.
  call ring_l(lhc,doneit)

END subroutine ptc_input
