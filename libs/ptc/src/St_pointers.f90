module pointer_lattice
  use madx_ptc_module !etienne, my_state_mad=>my_state, my_ring_mad=>my_ring
use duan_zhe_map, probe_zhe=>probe,tree_element_zhe=>tree_element,dp_zhe=>dp, &
DEFAULT0_zhe=>DEFAULT0,TOTALPATH0_zhe=>TOTALPATH0,TIME0_zhe=>TIME0,ONLY_4d0_zhe=>ONLY_4d0,RADIATION0_zhe=>RADIATION0, &
NOCAVITY0_zhe=>NOCAVITY0,FRINGE0_zhe=>FRINGE0,STOCHASTIC0_zhe=>STOCHASTIC0,ENVELOPE0_zhe=>ENVELOPE0, &
DELTA0_zhe=>DELTA0,SPIN0_zhe=>SPIN0,MODULATION0_zhe=>MODULATION0,only_2d0_zhe=>only_2d0 , &
INTERNAL_STATE_zhe=>INTERNAL_STATE

  !  use madx_keywords
  USE gauss_dis
  implicit none
  public
  ! stuff for main program
  type(layout),pointer :: my_ering => null(), my_fring => null()
  type(internal_state),pointer :: my_estate => null()
!  type(internal_state),pointer :: my_old_state
  integer ,pointer :: my_start, MY_ORDER, MY_NP,MY_END,my_start_t
  real(dp), pointer :: my_fix(:),MY_DELTA
  real(dp), pointer ::my_target(:)
  type(pol_block), pointer :: my_pol_block(:)
  integer, pointer :: n_pol_block
  real(sp) my_scale_planar
  ! end of stuff for main program
  integer :: lat(0:6,0:6,3)=1,mres(3),mresp
  character(25) lat_name(81)
  real(dp) r_ap(2),x_ap,y_ap   ! default aperture
  integer :: kind_ap=2,file_ap=0                      ! Single aperture position and kind of aperture + file unit
  character(nlp) name_ap,namet(2)
  character(255) :: filename_ap = "Tracking.txt",datafile
  character(255), private :: file_zher,filezhe, name_zhe
  integer, private :: k_zhe, number_zhe_maps_local
  integer last_npara 
  integer :: i_layout=0,i_layout_t=1,pos_layout=1
  integer my_lost_position
  private thin
  real(dp) thin
  !  BEAM STUFF
  REAL(DP) SIG(6) 
  REAL(DP) ait(6,6)
 
  logical :: nophase=.false. 
  INTEGER :: N_BEAM=0,USE_BEAM=1
  logical, private :: m_u_t = .true.
  TYPE(REAL_8),private :: Y(6)
  TYPE(DAMAP),PRIVATE :: ID
  !  integer nd2,npara
  !  PRIVATE POWER_CAVITY,RADIA
  ! stuff from my fortran
  type(internal_state), target:: etat
  integer,private,target:: START ,FIN,ORDER,np,start_t
  real(dp),target:: xfix(6) ,DELT0
  integer :: logs_exp=30, num_iter = 20
  !logical :: absolute = .false.
 !private ind1,ind2,ipos,ireal,a_f,a_f0,yfit
  integer, allocatable:: ind1(:),ind2(:),ipos(:)
  logical, allocatable :: ireal(:)
  real(dp), allocatable :: a_f(:),a_f0(:),yfit(:),dyfit(:,:)
  integer sizeind1
  logical :: onefunc = .true.,skipzero=.false.,skipcomplex=.true.
 type(probe), pointer :: xs0g(:) => null()
 logical ::  use_hermite =.false.
 private get_polarisation
real(dp) n_ang(3),lm(2)
character(vp)snake
integer isnake,nlm,no_pol

  INTERFACE SCRIPT
     MODULE PROCEDURE read_ptc_command
  END INTERFACE


  type hermite
    integer n
    real(dp) :: h(2)
    real(dp) :: a(6,6)=0
    real(dp):: ai(6,6)=0
    real(dp):: m(6,6)=0
    real(dp):: mi(6,6)=0
    real(dp) ::f(6)=0
    real(dp) b(2,2)
    integer :: gen =0
    type(damap), pointer :: ms(:,:)
    type(probe_8), pointer :: p(:,:)
    real(dp), pointer ::  x0(:,:,:) => null()
    type(internal_state) state
    type(layout), pointer ::r
    integer pos
    integer noh
    integer no
    integer :: maxite =100, nint=10
    logical :: linear =.true.
    real(dp) :: eps=1.e-6_dp
    real(dp), pointer ::  he(:,:,:,:,:) => null()
  end type hermite





contains



  subroutine set_lattice_pointers
    implicit none

    print77=.true.
    read77 =.true.

    my_ering => m_u%start
m_u_t=.true.
if(m_t%n>0) then
    my_ering => m_t%end
    m_u_t=.false.
endif
    if(associated(my_estate)) then
    !  my_estate=>my_old_state
      etat=my_estate
      my_estate => etat
    else
     etat=DEFAULT !+nocavity0-time0
     my_estate => etat
    endif 
    my_start => START
    my_END => FIN
    my_fix=>xfix
    my_DELTA=> DELT0
    MY_ORDER=> ORDER
    MY_NP=> NP
    MY_ORDER=1
    my_start_t=>start_t
    start_t=1
    DELT0=0.D0
    START=1
    FIN=1
    ORDER=1
    NP=0
    xfix=0.d0
    my_scale_planar=100.d0
    my_fix(1:6)=.000d0

    write(6,*) " absolute_aperture  ", absolute_aperture
    write(6,*) " hyperbolic_aperture ", hyperbolic_aperture

    ! CALL create_p_ring

  end  subroutine set_lattice_pointers


  subroutine read_ptc_command(ptc_fichier)
    use madx_ptc_module !etienne , piotr_state=>my_state,piotr_my_ring=>my_ring
    implicit none
    CHARACTER*(120) com,COMT,filename,name_root,title,name_root_res,filetune,FILESMEAR
    character*(4) suffix,SUFFIX_res
    character(*) ptc_fichier
    integer i,ii,mf,i_layout_temp,IB,NO,i_layout_temp_t
    !  FITTING FAMILIES
    INTEGER NPOL,J,NMUL,K,ICN,N,np
    type(pol_block), ALLOCATABLE :: pol_(:)
    type(pol_block) :: pb
    CHARACTER*(NLP) NAME,VORNAME
    real(dp) targ_tune(2),targ_chrom(2),EPSF
    real(dp) targ_RES(4)
    !  END FITTING FAMILIES
    real(dp),pointer :: beta(:,:,:)
    REAL(DP) DBETA,tune(3),tunenew(2),CHROM(2),DEL,dtune(3)
    integer ntune(2)
    ! fitting and scanning tunes
    real(dp) tune_ini(2),tune_fin(2),dtu(2),fint,hgap
    integer nstep(2),i1,i2,I3,n_bessel,i11,i22,di12
    LOGICAL(LP) STRAIGHT,skip,fixp,skipcav,fact
    ! end
    ! TRACK 4D NORMALIZED
    INTEGER POS,NTURN,resmax,ngen
    real(dp) EMIT(6),APER(2),emit0(2),sca,nturns
    integer nscan,mfr,ITMAX,MRES(4)
    real(dp), allocatable :: resu(:,:)
    ! END
    ! RANDOM MULTIPOLE
    INTEGER addi
    REAL(DP) CUT,cn,cns
    LOGICAL(LP) integrated
    ! END
    ! LOCAL BEAM STUFF
    INTEGER NUMBER_OF_PARTICLE
    TYPE(DAMAP) MY_A
    INTEGER MY_A_NO
    ! APERTURE
    REAL(DP)  APER_R(2),APER_X,APER_Y
    INTEGER KINDAPER
    TYPE(integration_node), POINTER :: TL
    type(internal_state),target :: my_default
    TYPE(internal_state) tempstate
    ! DYN APERTURE
    REAL(DP) r_in,del_in,DLAM,ang_in,ang_out,dx,targ_tune_alex(2),sexr0
    INTEGER ITE,n_in,POSR
    logical(lp) found_it
    type(fibre),pointer ::p,f1,f2,ft
    ! TRACKING RAYS
    INTEGER IBN,N_name
    REAL(DP) X(6),DT(3),sc,NLAM,A1,B1,HPHA,B_TESLA,CUR1,CUR2
    REAL(DP)VOLT,PHASE,x_ref(6),x_ref0(6)
    type(internal_state) state0
    logical do_state0
    INTEGER HARMONIC_NUMBER
    ! changing magnet
    logical(lp) bend_like
    logical exists,noca
    ! remove_patches
    save my_default

    LOGICAL :: b_b,patchbb
    REAL(DP) xbend
    ! automatic track
!    type(internal_state),pointer :: my_old_state
    TYPE(WORK) W
    INTEGER   KINDA   ! 1,2,3,4
    REAL(DP) RA(2)
    REAL(DP) XA,YA,DXA,DYA, DC_ac,A_ac,theta_ac,D_ac
    real(dp), allocatable :: an(:),bn(:) !,n_co(:)
    real(dp) d_volt,d_phas
    integer icnmin,icnmax,n_ac,inode,icavv !,n_coeff
    logical :: longprintt,onemap
!    logical :: log_estate=.true. 
    integer :: mftune=6,nc
    real(dp), allocatable :: tc(:)
    type(integration_node), pointer  :: t

     longprintt=longprint 
     longprint=.true.
!    if(log_estate) then
!       nullify(my_estate)
!       log_estate=.false.
!    endif

    if(associated(my_estate)) then
!      my_old_state=>my_estate
      my_default=my_estate
    else
      my_default=default
    endif 
    my_estate=>my_default
    skip=.false.
    call kanalnummer(mf)
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(6,*) "$$$$$$$$$$         Write New read_ptc_command           $$$$$$$$$$"
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

    if(ptc_fichier/="screen") then
     open(unit=mf,file=ptc_fichier)
    else
     mf=5
    endif
    if(i_layout==0) then
             i_layout=1
             call move_to_layout_i(m_u,my_ering,i_layout)
             write(6,*) "Selected Layout in m_u",i_layout,"  called  ---> ",my_ering%name
             m_u_t=.true.
    endif
    do i=1,10000
       read(mf,'(a120)') comT
       


       COM=COMT
       call context(com)
       com=com(1:len_trim(com))
       if(index(com,'!')/=0) then
         if(index(com,'!')/=1) then
          comt=com(1:index(com,'!')-1)
          COM=COMT
          call context(com)
         endif
       endif


       if(com(1:1)==' ') THEN
          WRITE(6,*) ' '
          cycle
       ENDIF
  !     if(com(1:1)=='!'.and.com(2:2)/='!') THEN
       if(com(1:1)=='!') THEN
          cycle
       ENDIF
       if(com(1:5)=='PAUSE') THEN
          com=com(1:5)
       ENDIF
       if(.not.skip) then
          if(com(1:2)=="/*") then
             skip=.true.
             cycle
          endif
       endif
       if(skip) then !1
          if(com(1:2)=="*/") then
             skip=.false.
          endif
          cycle
       endif         ! 1
       write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       write(6,*) " "
       write(6,*) '            ',i,comT(1:LEN_TRIM(COMT))
       write(6,*) " "
       write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

       select case(com)
       case('SELECTLAYOUT','SELECTLATTICE')
          read(mf,*) i_layout_temp
          if(i_layout_temp>m_u%n) then
             write(6,*) " Universe Size ", m_u%n

             write(6,*) " Selected Layout does not exist "

          else

             i_layout=i_layout_temp
             call move_to_layout_i(m_u,my_ering,i_layout)
             write(6,*) "Selected Layout in m_u",i_layout,"  called  ---> ",my_ering%name
             m_u_t=.true.
          endif
       case('SELECTTRACKABLELAYOUT','SELECTTRACKABLELATTICE')
          read(mf,*) i_layout_temp
          if(i_layout_temp>m_t%n) then
             write(6,*) " Universe Size ", m_t%n

             write(6,*) " Selected Layout does not exist "

          else

             i_layout=i_layout_temp
             call move_to_layout_i(m_t,my_ering,i_layout)
             write(6,*) "Selected Layout in m_t",i_layout,"  called  ---> ",my_ering%name
             m_u_t=.false.
          endif
       case('SELECTLASTLAYOUT','SELECTLASTLATTICE')
         if(m_u_t) then
          i_layout_temp=m_u%n
          if(i_layout_temp>m_u%n) then
             write(6,*) " Universe Size ", m_u%n

             write(6,*) " Selected Layout does not exist "

          else

             i_layout=i_layout_temp
             call move_to_layout_i(m_u,my_ering,i_layout)
             write(6,*) "Selected Layout in m_u",i_layout,"  called  ---> ",my_ering%name
          endif
         else
          i_layout_temp=m_t%n
          if(i_layout_temp>m_t%n) then
             write(6,*) " Universe Size ", m_t%n

             write(6,*) " Selected Layout does not exist "

          else

             i_layout=i_layout_temp
             call move_to_layout_i(m_t,my_ering,i_layout)
             write(6,*) "Selected Layout in m_t",i_layout,"  called  ---> ",my_ering%name
          endif
         endif
       case('NAMELAYOUT')
          read(mf,*) NAME
       case('KILLLASTLAYOUT')
          read(mf,*) NAME
          call kill_last_layout(m_u)
          my_ering%NAME=NAME
!       case('RESTOREDEFAULT','RESTORE')
!          my_OLD_state=default
!          my_default=default
!          my_estate=>my_default
          ! Orbit stuff

       case('UPDATETWISS','UPDATETWISSFORORBIT')
           call update_twiss_for_orbit
       case('USEORBITUNITS')
          MY_ERING%t%ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=.true.
       case('DONOTUSEORBITUNITS')
          MY_ERING%t%ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=.false.
      case('FORCERESTORALLMAGNETS')
          CALL force_restore_ANBN(my_ering)
      case('RESTORALLMAGNETS')
          CALL restore_ANBN(my_ering)

       case('SETORBITRESTORE','RESTOREMAGNETFLAG')
          read(mf,*) restore_mag,restore_magp
!       case('SETORBITPHASOR')
!          read(mf,*) xsm%ac%x(1:2)
!          xsm0%ac%x(1:2)=xsm%ac%x(1:2)
!          ptc_node_old=-1
!          first_particle=my_false
       case('SETORBITOMEGA')   !!! TO CREATE A NODE
        read(mf,*) xbend
         if(associated(MY_ERING%t)) then
          if(associated(MY_ERING%t%ORBIT_LATTICE)) then
         write(6,*) " Old ORBIT_OMEGA = ",MY_ERING%t%ORBIT_LATTICE%ORBIT_OMEGA
        if(xbend<0) then 
          MY_ERING%t%ORBIT_LATTICE%ORBIT_OMEGA=abs(xbend)
        else
          MY_ERING%t%ORBIT_LATTICE%ORBIT_OMEGA=MY_ERING%t%ORBIT_LATTICE%ORBIT_OMEGA*xbend
        endif
         write(6,*) " New ORBIT_OMEGA = ",MY_ERING%t%ORBIT_LATTICE%ORBIT_OMEGA
        endif
        endif 
      case('SETORBITMARKER')
         READ(MF,*) i1
         allocate(orbitname(i1))
          do i1=1,size(orbitname)
           READ(MF,*) name
           call context(name)
           orbitname(i1)=name
          enddo
       case('USER1')
         call my_user_routine1
       case('SPINTWISSCAS')
           read(mf,*) no,noca   !  order and general canonise
         call SMALL_CODE_TWISS(my_ering,no,noca)
       case('SETORBITPHASORTIME','ORBITTIME')
          read(mf,*) xsmt
          xsm0t=xsmt
          x_orbit_sync=0.0_dp
          x_orbit_sync(6)=xsmt
          ptc_node_old=-1
          first_particle=my_false       
       case('MAKEORBITMARKER')
            extra_node=1
            if(associated(my_ering%t)) then
             if(associated(my_ering%t%ORBIT_LATTICE)) then
             write(6,*) "The orbit nodes have been already created "
             write(6,*) '"MAKE ORBIT MARKER" will not work!' 
             write(6,*) 'If running orbit put this command in the special ' 
             write(6,*) ' PTC script called pre_orbit_set.txt ' 
             write(6,*) ' Execution now interrupted ' 
               stop
             endif   
            endif    
       case('TIMEINUNITS','TIMEINSECONDS')
       
          read(mf,*) xsmt
!          write(6,*) " Using ",unit_time," seconds"
          xsmt=xsmt*clight  !*unit_time
          xsm0t=xsmt
          x_orbit_sync=0.0_dp
          x_orbit_sync(6)=xsmt
          ptc_node_old=-1
          first_particle=my_false
       case('INITIALTIMEINMYUNITS','TIMEINMYUNITS','READFROMFILEINITIALTIME')
         INQUIRE (FILE = INITIAL_setting, EXIST = exists)
          if(exists) then
             write(6,*) "file ", INITIAL_setting(1:len_trim(FINAL_setting)), &
                  " exists, interrupt execution if you do not want to overwrite!"
            call kanalnummer(i1,INITIAL_setting)
   !         read(i1,*) xsm%ac%t,unit_time,n_used_patch,include_patch
            read(i1,*) xsmt,n_used_patch,include_patch
            read(i1,*) nc
            if(nc/=0) then
                 allocate(tc(nc))
                  do i2=1,nc
                   read(i1,*) tc(i2) 
                  enddo
                 call set_all_tc_for_restarting(my_ering,tc,nc)
                 deallocate(tc)
            endif
            close(i1)
            if(include_patch) then
              call kanalnummer(i1,"time_patch.dat")
               read(i1,*) n_patch
               if(associated(my_ering%T%ORBIT_LATTICE%dt)) deallocate(my_ering%T%ORBIT_LATTICE%dt)
               allocate(my_ering%T%ORBIT_LATTICE%dt(n_patch))
                do i2=1,n_patch
                 read(i1,*) i3,my_ering%T%ORBIT_LATTICE%dt(i2)
                enddo
              close(i1)
!looking for element just before the cavity
   do i2=1,size(my_ering%T%ORBIT_LATTICE%ORBIT_NODES)
     T=>my_ering%T%ORBIT_LATTICE%ORBIT_NODES(i2)%NODE
  
    DO I1=1,my_ering%T%ORBIT_LATTICE%ORBIT_NODES(i2)%dpos
 

       if(t%parent_fibre%mag%kind==kind4) then
        my_ering%T%ORBIT_LATTICE%tp=>t%previous
        inode=i2
        goto 2222
       endif

       T=>T%NEXT
    ENDDO
   enddo
2222   write(6,*) "ptc mode # ",inode,"element ", my_ering%T%ORBIT_LATTICE%tp%parent_fibre%mag%name
  



            endif
           else
         read(mf,*) xsmt,n_used_patch,include_patch
!             read(mf,*) xsm%ac%t,unit_time,n_used_patch,include_patch
          endif
 !         write(6,*) " Using ",unit_time," seconds"
          xsmt=xsmt*clight 
          xsm0t=xsmt
          x_orbit_sync=0.0_dp
          x_orbit_sync(6)=xsmt
          ptc_node_old=-1
          first_particle=my_false
          write(6,*) "x_orbit_sync(6) = " , x_orbit_sync(6)
          
       case('FINALTIMEINMYUNITS')
      ! INQUIRE (FILE = FINAL_setting, EXIST = exists)
      !    if(exists) then
            call kanalnummer(i1,FINAL_setting)
             write(i1,*) x_orbit_sync(6)/clight,n_used_patch,include_patch
             write(6,*) "x_orbit_sync(6) = " , x_orbit_sync(6)
             write(6,*) "t_fin = " , x_orbit_sync(6)/clight
             
             call find_all_tc_for_restarting(my_ering,tc,nc)
             write(i1,*) nc
             do i2=1,nc
                write(i1,*) tc(i2)
             enddo
             deallocate(tc)
            close(i1)
      !    endif              

      
       case('SETORBITACCELERATION')
          accelerate=my_true          
       case('SETORBITNOACCELERATION')
          accelerate=my_false          
       case('SETORBITRAMPING','RAMPING')
          RAMP=my_true          
       case('SETORBITNORAMPING','NORAMPING')
          RAMP=my_false          
 !      case('SETORBITTIMEUNIT')
 !           read(mf,*) unit_time
       case('SETORBITSTATE')
          my_ORBIT_LATTICE%state=my_estate
       case('PUTORBITSTATE','USEORBITSTATE')
          my_estate=my_ORBIT_LATTICE%state
       case('NULLIFYACCELERATION')
        nullify(acc)
        nullify(ACCFIRST)
        nullify(paccfirst)
        nullify(paccthen)
       case('INITIALIZECAVITY','CAVITYTABLE')
        p=>my_ering%start
          READ(MF,*)n
          do i1=1,n
            read(mf,*) name,filename
            call move_to( my_ering,p,name,pos) 
            write(6,*) "Found cavity ",name," at position ",pos
            call lecture_fichier(p%mag,filename)        
          enddo
           paccthen%mag%c4%acc%next=>paccfirst
           paccfirst%mag%c4%acc%previous=>paccthen
       case('ENERGIZEORBITLATTICEATTIME')
        read(mf,*) xa
        call energize_ORBIT_lattice(xa)
      case('ENERGIZEORBITLATTICE')
        call energize_ORBIT_lattice
       case('SETALLRAMP')
        call set_all_ramp(my_ering)
        
       case('PRINTHERD','OPENHERD')
          IF(MF_HERD/=0) THEN
             WRITE(6,*) " CANNOT OPEN HERD FILE TWICE "
             STOP 55
          ENDIF
          call kanalnummer(MF_HERD)
          open(unit=MF_HERD,file=print_herd)
       case('CLOSEHERD')
          CLOSE(MF_HERD)
          MF_HERD=0
       case('PRINTPTCNODES','PRINTORBITNODES')
          call kanalnummer(ii)
          open(unit=ii,file=def_orbit_node)

          tL=>my_ORBIT_LATTICE%ORBIT_NODES(1)%NODE
          write(ii,*) " Number of PTC Nodes  from 0 to ",my_ORBIT_LATTICE%ORBIT_N_NODE-1

          do j=1,my_ORBIT_LATTICE%ORBIT_N_NODE
             write(ii,*) "****************************************************"
             write(ii,*) " Node Number ",j-1
             write(ii,*) " Number of PTC Integration Nodes",my_ORBIT_LATTICE%ORBIT_NODES(j)%dpos
             DO I1=1,my_ORBIT_LATTICE%ORBIT_NODES(j)%dpos
                write(ii,*) tL%S(1),tL%parent_fibre%mag%name,' ',tL%parent_fibre%pos,tL%cas
                TL=>TL%NEXT
             ENDDO
          enddo
          close(ii)

       case('SETORBITFORACCELERATION')
       !   CALL ptc_synchronous_set(-1)
       case('READORBITINTERMEDIATEDATA')
       !   CALL ptc_synchronous_set(-2)
       case('PRINTORBITINTERMEDIATEDATA')
       !   CALL ptc_synchronous_after(-2)
       case('FILEFORORBITINITIALSETTINGS')  ! ACCELERATION FILE
          READ(MF,*) initial_setting
          write(6,*) initial_setting

          !         INQUIRE (FILE = initial_setting, EXIST = exists)
          !        if(exists) then
          !         write(6,*) "file ", initial_setting(1:len_trim(def_orbit_node)), &
          !         " exists, interrupt execution if you do not want to overwrite!"
          !        endif
       case('FILEFORORBITFINALSETTINGS')  ! ACCELERATION FILE
          READ(MF,*) final_setting
          write(6,*) final_setting
          !          INQUIRE (FILE = final_setting, EXIST = exists)
          !         if(exists) then
          !          write(6,*) "file ", final_setting(1:len_trim(def_orbit_node)), &
          !          " exists, interrupt execution if you do not want to overwrite!"
          !         endif
       case('PTCNODEFILE','PTCNODE')  ! ACCELERATION FILE
          READ(MF,*) def_orbit_node
          write(6,*) def_orbit_node
          !        INQUIRE (FILE = def_orbit_node, EXIST = exists)
          !       if(exists) then
          !        write(6,*) "file ", def_orbit_node(1:len_trim(def_orbit_node)), &
          !        " exists, interrupt execution if you do not want to overwrite!"
          !    write(6,*) " you have 2 seconds to do so "
          !    CALL DATE_AND_TIME(values=temps)
          !    i1=temps(7)
          if(i1>=58) i1=i1-58
          !    do while(.true.)
          !     CALL DATE_AND_TIME(values=temps)
          !     if(temps(7)>i1+2) exit
          !    enddo
          !          endif
          ! end Orbit stuff
       case('RECORDLOSTPARTICLEINORBIT')  !
          wherelost=1
       case('RECORDLOSTPARTICLE')  ! 1= orbit , 2 thin lens PTC
          read(mf,*) wherelost
       case('NOLOSTPARTICLE')  !
          wherelost=0
       case('CLEANLOSTPARTICLE')  !
          if(associated(my_ering%T)) then
             TL=>my_ering%T%START
             DO j=1,my_ering%T%N
                tl%lost=0
                TL=>TL%NEXT
             ENDDO
          endif
       case('PRINTLOSTPARTICLE')  !
          read(mf,*) filename
          CALL  kanalnummer(mfr)
          OPEN(UNIT=MFR,FILE=FILENAME,status='REPLACE')
          write(mfr,'(a42)') "  s          , lost  , cas , pos ,  name "
          if(associated(my_ering%T)) then
             TL=>my_ering%T%START
             DO j=1,my_ering%T%N
                if(TL%lost>0) then
                   write(mfr,900) tl%s(1),TL%lost,tl%cas,TL%pos_in_fibre,tl%parent_fibre%mag%name
900                FORMAT(F13.6,1x,(i8,1x),2(i4,1x),a16)
                endif
                TL=>TL%NEXT
             ENDDO
          endif
          close(mfr)
       case('DEFAULT')
          my_estate=DEFAULT0
       case('RADIATION')
          my_estate=RADIATION          
       case('+NOCAVITY')
          my_estate=my_estate+NOCAVITY0
       case('-NOCAVITY')
          my_estate=my_estate-NOCAVITY0
       case('+ENVELOPE')
          my_estate=my_estate+ENVELOPE0
       case('-ENVELOPE')
          my_estate=my_estate-ENVELOPE0
       case('+CAVITY')
          my_estate=my_estate-NOCAVITY0
       case('-CAVITY')
          my_estate=my_estate+NOCAVITY0
       case('+FRINGE')
          my_estate=my_estate+FRINGE0
       case('-FRINGE')
          my_estate=my_estate-FRINGE0
       case('+TIME')
          my_estate=my_estate+TIME0
       case('-TIME')
          my_estate=my_estate-TIME0
       case('+TOTALPATH')
          my_estate=my_estate+TOTALPATH0
       case('-TOTALPATH')
          my_estate=my_estate-TOTALPATH0
       case('+ONLY_4D')
          my_estate=my_estate+ONLY_4D0
       case('-ONLY_4D')
          my_estate=my_estate-ONLY_4D0
       case('+DELTA')
          my_estate=my_estate+DELTA0
       case('-DELTA')
          my_estate=my_estate-DELTA0
       case('+RADIATION')
          my_estate=my_estate+RADIATION0
       case('-RADIATION')
          my_estate=my_estate-RADIATION0
       case('+MODULATION')
          my_estate=my_estate+MODULATION0
       case('-MODULATION',"+NOMODULATION")
          my_estate=my_estate-MODULATION0
       case('+SPIN')
          my_estate=my_estate+SPIN0

       case("PRINTSTATE")
         READ(MF,*) K
         CALL print(MY_ESTATE,K)
       case("PRINTSTATEONSCREEN")
         CALL print(MY_ESTATE,6)

       case("KILLPROBE","KILLPROBES")
        deallocate(xs0g)
        nullify(xs0g)
       case("READPROBE","READPROBES")
        read(mf,*) file_zher
        call  read_ptc_rays(file_zher)

       case('BERZ','GERMANIC','MARTIN')
          CALL change_package(2)
       case('YANG','CHINESE','LINGYUN')
          CALL change_package(1)
       case('ALLOCATEBETA')
          ALLOCATE(BETA(2,3,my_ering%N))

       case('DEALLOCATEBETA')
          DEALLOCATE(BETA)
       case('FILLBETA')
          READ(MF,*) IB,pos
          CALL FILL_BETA(my_ering,my_estate,pos,BETA,IB,DBETA,tune,tunenew)
       case('FITBENDS')
          CALL fit_all_bends(my_ering,my_estate)
       case('LIMITFORCUTTING')
          READ(MF,*) limit_int0
          WRITE(6,*) "limit_int0 =",limit_int0
       case('LMAX',"MAXIMUMDS")
          READ(MF,*) LMAX
          WRITE(6,*) "LMAX FOR SPACE CHARGE =",LMAX
       case('FUZZYLMAX')
          READ(MF,*) FUZZY_SPLIT
          WRITE(6,*) "FUZZY LMAX FOR SPACE CHARGE =",abs(LMAX),abs(LMAX*FUZZY_SPLIT)
       case('CUTTINGALGORITHM')
          READ(MF,*) resplit_cutting
          IF(resplit_cutting==0) WRITE(6,*) " RESPLITTING NORMALLY"
          IF(resplit_cutting==1) WRITE(6,*) " CUTTING DRIFT USING LMAX"
          IF(resplit_cutting==2) WRITE(6,*) " CUTTING EVERYTHING USING LMAX"
          IF(resplit_cutting==-2) WRITE(6,*) " CUTTING EVERYTHING USING LMAX EXCEPT DRIFTS"
          !       case('KIND7WITHMETHOD1')
          !          CALL PUT_method1_in_kind7(my_ering,1000)

       case('RADIATIONBENDSPLIT','RADIATIONBEND')
         read(mf,*) radiation_bend_split
       case('THINLENS=1')
          call THIN_LENS_restart(my_ering)
       case('MANUALTHINLENS')
          sagan_even=my_true
          THIN=-1
          CALL THIN_LENS_resplit(my_ering,THIN,lim=limit_int0)
       case('THINLENS')
         sagan_even=my_true
          READ(MF,*) THIN
          xbend=-1.0_dp
          sexr0=-1.d0
          if(thin<0) then
             READ(MF,*) sexr0,xbend
             if(xbend<0) then
                xbend=-xbend
                radiation_bend_split=MY_true
             endif
             thin=-thin
          endif
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(my_ering,THIN,lim=limit_int0,lmax0=lmax,sexr=sexr0,xbend=xbend)
          radiation_bend_split=MY_false
       case('EVENTHINLENS')
          READ(MF,*) THIN
          xbend=-1.0_dp
          sexr0=-1.d0
          if(thin<0) then
             READ(MF,*) sexr0,xbend
             if(xbend<0) then
                xbend=-xbend
                radiation_bend_split=MY_true
             endif
             thin=-thin
          endif
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(my_ering,THIN,EVEN=my_TRUE,lim=limit_int0,lmax0=lmax,sexr=sexr0,xbend=xbend)
       case('ODDTHINLENS')
          sagan_even=my_false
          READ(MF,*) THIN
          xbend=-1.0_dp
          sexr0=-1.d0
          if(thin<0) then
             READ(MF,*) sexr0,xbend
             if(xbend<0) then
                xbend=-xbend
                radiation_bend_split=MY_true
             endif
             thin=-thin
          endif
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(my_ering,THIN,EVEN=my_FALSE,lim=limit_int0,lmax0=lmax,sexr=sexr0,xbend=xbend)
          ! thin layout stuff

       case('RECUTKIND7NODRIFT','USEODDMETHODS')
          call RECUT_KIND7(my_ering,lmax,my_false)
       case('RECUTKIND7ANDDRIFT','USEODDMETHODSONDRIFT')
          call RECUT_KIND7(my_ering,lmax,my_true)
       case('MAKE_THIN_LAYOUT','MAKELAYOUT','MAKE_NODE_LAYOUT')

          if(.not.associated(my_ering%t)) CALL MAKE_node_LAYOUT(my_ering)
       case('SURVEY_THIN_LAYOUT','SURVEYLAYOUT','SURVEY_NODE_LAYOUT')

          IF(associated(my_ering%t)) THEN
             CALL fill_survey_data_in_NODE_LAYOUT(my_ering)
          ELSE
             WRITE(6,*) " NO NODE LAYOUT PRESENT "
          ENDIF


       case('CHECKKREIN')
          WRITE(6,*) "OLD CHECK_KREIN ",check_krein
          READ(MF,*) check_krein
          WRITE(6,*) "NEW CHECK_KREIN ",CHECK_KREIN

       case('KREINSIZE')
          WRITE(6,*) "OLD KREIN SIZE PARAMETER ",size_krein
          READ(MF,*) size_krein
          WRITE(6,*) "NEW KREIN SIZE PARAMETER",size_krein

       case('ABSOLUTEAPERTURE')
          WRITE(6,*) "OLD C_%ABSOLUTE_APERTURE ",C_%ABSOLUTE_APERTURE
          READ(MF,*) C_%ABSOLUTE_APERTURE
          WRITE(6,*) "NEW C_%ABSOLUTE_APERTURE ",C_%ABSOLUTE_APERTURE
          ! end of layout stuff
          ! random stuff
       case('GLOBALAPERTUREFLAG')
          WRITE(6,*) "OLD C_%APERTURE_FLAG ",C_%APERTURE_FLAG
          READ(MF,*) APERTURE_FLAG
          WRITE(6,*) "NEW C_%APERTURE_FLAG ",C_%APERTURE_FLAG
       case('TURNOFFONEAPERTURE','TOGGLEAPERTURE')
          READ(MF,*) POS
          CALL TURN_OFF_ONE_aperture(my_ering,pos)
       case('SETONEAPERTURE')
          read(MF,*) pos
          read(MF,*) kindaper, APER_R,APER_X,APER_Y,dxa,dya
          CALL assign_one_aperture(my_ering,pos,kindaper,APER_R,APER_X,APER_Y,dxa,dya)
          ! end of layout stuff
          ! random stuff

       case('GAUSSIANSEED')
          READ(MF,*) I1
          CALL gaussian_seed(i1)
 
       case('RANDOMMULTIPOLE')
          read(mf,*) i1, cut,STRAIGHT     !!!   seed, cut on sigma, print
          read(mf,*) name,fixp    !!!  name, full=> logical true full name compared
          read(mf,*) sc  !!!  percentage 

          if(sc<=0) then
            read(mf,*) addi,integrated
            read(mf,*) n,cns,cn      !!! n (negative means skew) bn=cns+cn*x
          endif
          call lattice_random_error_new(my_ering,name,fixp,i1,cut,n,addi,integrated,cn,cns,sc,STRAIGHT)
       case('MISALIGNEVERYTHING')
 
          read(mf,*) SIG(1:6),cut
          CALL MESS_UP_ALIGNMENT(my_ering,SIG,cut)
          ! end of random stuff
       case('MISALIGN','MISALIGNONE')
          read(mf,*) i1, cut,STRAIGHT
          read(mf,*) name,fixp
          read(mf,*) SIG(1:6) 
          call  MESS_UP_ALIGNMENT_name(my_ering,name,i1,fixp,sig,cut,STRAIGHT)
       case('ALWAYS_EXACTMIS',"ALWAYSEXACTMISALIGNMENTS")
          ! end of random stuff
          read(mf,*) ALWAYS_EXACTMIS
          if(ALWAYS_EXACTMIS) write(6,*) " ALWAYS EXACT MISALIGNMENTS "
          if(.NOT.ALWAYS_EXACTMIS) write(6,*) " EXACT MISALIGNMENTS SET USING STATE "
       case('KILLBEAMBEAM')
          if(associated(my_ering%T)) then
             TL=>my_ering%T%START
             DO j=1,my_ering%T%N
                if(associated(tl%BB)) then
                   write(6,*) tl%pos,tl%parent_fibre%mag%name,' killed'
                   call kill(tl%BB)
                endif
                TL=>TL%NEXT
             ENDDO
          endif
       case('E1CAS','E1WM')
       read(MF,*) sc
       if(sc/=e1_cas) then
         write(6,*) "E1 of the wonderful magnet company was ", e1_cas
         write(6,*) "                     |  "
         write(6,*) "                     |  "
         write(6,*) "                    \|/  "
           e1_cas=sc
         write(6,*) "E1 of the wonderful magnet company is now ", e1_cas     
        endif
       case('BEAMBEAM')
          READ(MF,*) SC,pos,patchbb
          read(mf,*) X_ref(1), X_ref(2), X_ref(3), X_ref(4)
          if(patchbb) then
           read(mf,*) x
          endif
          IF(.NOT.ASSOCIATED(my_ering%T)) THEN
             CALL MAKE_NODE_LAYOUT(my_ering)
          ENDIF
          ! s(1) total ld
          ! s(2) local integration distance
          !          SC=MOD(SC,MY_RING%T%END%S(1))
          b_b=.false.
          TL=>my_ering%T%START
          DO j=1,my_ering%T%N
             if(pos<1) then
                IF(TL%S(1)<=SC.AND.TL%NEXT%S(1)>SC) then
                   b_b=.true.
                   exit
                endif
             else
                if(j==pos) then
                   b_b=.true.
                   exit
                endif
             endif
             TL=>TL%NEXT
          ENDDO
          if(b_b.and.tl%cas/=case0) tl=>tl%next
          if(b_b.and.tl%cas==case0) then
             write(6,*) " Beam-Beam position at ",tl%parent_fibre%mag%name
             if(.not.associated(tl%BB)) call alloc(tl%BB)
             tl%bb%fk=X_ref(1)* X_ref(4)**2
             tl%bb%sx=X_ref(2)* X_ref(4)
             tl%bb%sy=X_ref(3)* X_ref(4)
             !           if(pos<1) tl%bb%ds=SC-TL%S(1)
             write(6,*) tl%pos,tl%parent_fibre%mag%name,' created'
             !              write(6,*) " ds = ",tl%bb%ds
             if(patchbb) then
              tl%bb%patch=patchbb
              tl%bb%a=x(1:3)
              tl%bb%d=x(4:6)              
             endif
          else
             write(6,*) " Beam-Beam position not found "
          endif

       case('TAPERING','TAPER')
           read(mf,*) i1    ! staircase for reachuing full radiation
           read(mf,*) filename   !!!  put "nofile" for no plot
           X_ref=0
           p=>my_ering%start
         call taper(p,X_ref,i1,my_estate,1.d-8,file=filename)  
       case('UNTAPERING','UNTAPER')
           p=>my_ering%start
           call untaper(p)  
          case('DOBEAMBEAM')
          do_beam_beam=my_true       
       case('NOBEAMBEAM')
          do_beam_beam=my_false
       case('SETFAMILIES','SETGFAMILIES')
          np=0
          READ(MF,*) NPOL
          ALLOCATE(pol_(NPOL))
          DO J=1,NPOL
             READ(MF,*) NMUL,NAME
             CALL CONTEXT(NAME)
             N_NAME=0
             IF(NAME(1:2)=='NO') THEN
                READ(MF,*) NAME
                call context(name)
                N_NAME=len_trim(name)
             ENDIF
             POL_(J)=0
             POL_(J)%NAME=NAME
             POL_(J)%N_NAME=N_NAME
             DO K=1,NMUL
                IF(COM=='SETGFAMILIES') THEN
                   READ(MF,*) N,ICN,POL_(J)%G,POL_(J)%NB,POL_(J)%NP  !     integer g,np,nb   !  group index  number of blocks
                ELSE
                   READ(MF,*) N,ICN
                ENDIF
                if(icn>np) np=icn
                IF(N>0) THEN
                   POL_(J)%IBN(N)=ICN
                ELSE
                   POL_(J)%IAN(-N)=ICN
                ENDIF
             ENDDO
          ENDDO
          Write(6,*) " Number of parameters = ",np
          Write(6,*) " Number of polymorphic blocks = ",NPOL

       case('SETGFAMILIESN')
          np=0
          READ(MF,*) NPOL
          ALLOCATE(pol_(NPOL))
          kindaper=0
          n_in=0
2010      continue
          READ(MF,*) kinda
          kindaper=kindaper+kinda
          n_in=n_in+1
          icnmin=100000
          icnmax=0
          DO J=kindaper-kinda+1,kindaper
             READ(MF,*) NMUL,NAME
             CALL CONTEXT(NAME)
             N_NAME=0
             IF(NAME(1:2)=='NO') THEN
                READ(MF,*) NAME
                call context(name)
                N_NAME=len_trim(name)
             ENDIF
             POL_(J)=0
             POL_(J)%NAME=NAME
             POL_(J)%N_NAME=N_NAME
             DO K=1,NMUL
                IF(COM=='SETGFAMILIES') THEN
                   READ(MF,*) N,ICN,POL_(J)%G,POL_(J)%NB,POL_(J)%NP  !     integer g,np,nb   !  group index  number of blocks
                ELSE
                   READ(MF,*) N,ICN
                   POL_(J)%NB=n_in
                ENDIF
                if(icn>np) np=icn
                if(icn<icnmin) icnmin=icn
                if(icn>icnmax) icnmax=icn
                IF(N>0) THEN
                   POL_(J)%IBN(N)=ICN
                ELSE
                   POL_(J)%IAN(-N)=ICN
                ENDIF
             ENDDO
          ENDDO
          DO J=kindaper-kinda+1,kindaper
             POL_(J)%g=icnmin
             POL_(J)%NP=icnmax-icnmin+1
          enddo
          if(kindaper<npol) goto 2010
          Write(6,*) " Number of parameters = ",np
          Write(6,*) " Number of polymorphic blocks = ",NPOL

       case('READFAMILIESANBN')
          read(mf,*) filename
          call kanalnummer(mfr,filename)
          READ(MFR,*) II
          DO J=1,II
             read(mfr,'(a120)') title
             READ(title,*) NAME,VORNAME
             call move_to( my_ering,p,name,VORNAME,pos)

             call context(title)
             do while(.true.)
                read(mfr,'(a120)') title
                name_root=title
                call context(name_root)
                if(name_root(1:3)=='END') exit
                read(title,*) n,r_in
                call add(p,n,0,r_in)
             enddo

          ENDDO

          close(mfr)

       case('FRANKCOMMAND','MADCOMMAND','COMMAND')
          read(mf,*) filename
          call read_mad_command77(filename)
          
       case('PRINTFAMILIESANBN')
          read(mf,*) filename
          call kanalnummer(mfr,filename)
          write(MFR,*) NPOL
          DO J=1,NPOL
             my_ering=pol_(J)
          ENDDO

          p=>my_ering%start

          do j=1,my_ering%n
             if(p%magp%knob) then
                write(mfr,*)  p%magp%name, p%magp%vorname
                do ii=1,p%magp%p%nmul
                   if(p%magp%bn(ii)%kind==3) then
                      write(mfr,*)ii, p%magp%bn(ii)%r
                   endif
                   if(p%magp%an(ii)%kind==3) then
                      write(mfr,*)-ii, p%magp%an(ii)%r
                   endif
                enddo
                write(mfr,*) "end of data"
             endif
             p=>p%next
          enddo
          call kill_para(my_ering)
          close(mfr)
       case('RAMP','RAMPMAGNET')
          READ(MF,*) NAME, filename, hgap

          CALL CONTEXT(NAME)
          N_NAME=0
          IF(NAME(1:2)=='NO') THEN
             READ(MF,*) NAME
             call context(name)
             N_NAME=len_trim(name)
          ENDIF
            DC_ac=1.0_dp
            A_ac=0.0_dp
            theta_ac=0.0_dp
            D_ac=0.0_dp
            n_ac=0
           ! n_coeff=0
          p=>my_ering%start
          do ii=1,my_ering%N
             found_it=MY_FALSE
             IF(N_NAME==0) THEN
                found_it=P%MAG%NAME==NAME
             ELSE
                found_it=P%MAG%NAME(1:N_NAME)==NAME(1:N_NAME)
             ENDIF

             IF(FOUND_IT) THEN
             
                write(6,*) " slow ramping magnet found ",P%MAG%name
                
                call reading_file(P%MAG,filename)
                p%mag%ramp%r=hgap
                p%magp%ramp%r=hgap
                
                i2=size(p%mag%ramp%table(0)%bn)
                
                if(.not.associated(P%MAG%DC_ac)) then
                   allocate(P%MAG%DC_ac)
                   allocate(P%MAG%A_ac)
                   allocate(P%MAG%theta_ac)
                   allocate(P%MAG%D_ac)

                   allocate(P%MAGP%DC_ac)
                   allocate(P%MAGP%A_ac)
                   allocate(P%MAGP%theta_ac)
                   CALL alloc(P%MAGP%DC_ac)
                   CALL alloc(P%MAGP%A_ac)
                   CALL alloc(P%MAGP%theta_ac)
                   allocate(P%MAGp%D_ac)
                   CALL alloc(P%MAGP%D_ac)



                   P%MAG%D_ac=D_ac
                   P%MAG%DC_ac=DC_ac
                   P%MAG%A_ac=A_ac
                   P%MAG%theta_ac=theta_ac*twopi
                   P%MAGP%D_ac=D_ac
                   P%MAGP%DC_ac=DC_ac
                   P%MAGP%A_ac=A_ac
                   P%MAGP%theta_ac=theta_ac*twopi
                   P%MAG%slow_ac=1
                   P%MAGP%slow_ac=1

                   if(i2>p%mag%p%nmul) then
                      CALL ADD(P,i2,0,0.0_dp)
                   endif
                   allocate(P%MAG%d_an(p%mag%p%nmul))
                   allocate(P%MAG%d_bn(p%mag%p%nmul))
                   allocate(P%MAGp%d_an(p%mag%p%nmul))
                   allocate(P%MAGp%d_bn(p%mag%p%nmul))
                   allocate(P%MAG%d0_an(p%mag%p%nmul))
                   allocate(P%MAG%d0_bn(p%mag%p%nmul))
                   allocate(P%MAGp%d0_an(p%mag%p%nmul))
                   allocate(P%MAGp%d0_bn(p%mag%p%nmul))


                   P%MAG%d_an=0.0_dp
                   P%MAG%d_bn=0.0_dp

                   call alloc(P%MAGp%d_an,p%mag%p%nmul)
                   call alloc(P%MAGp%d_bn,p%mag%p%nmul)
                   call alloc(P%MAGp%d0_an,p%mag%p%nmul)
                   call alloc(P%MAGp%d0_bn,p%mag%p%nmul)
                   do i1=1,p%mag%p%nmul
                      P%MAG%d0_bn(i1)=P%MAG%bn(i1)
                      P%MAG%d0_an(i1)=P%MAG%an(i1)
                      P%MAGp%d0_bn(i1)=P%MAG%bn(i1)
                      P%MAGp%d0_an(i1)=P%MAG%an(i1)
                   enddo


                else
                   write(6,*) " already associated --> change code "
                   stop 166
                ENDIF
             ENDIF


             p=>p%next
          enddo

       case('RFCLEANBMAD')
     !  read(mf,*)  IBN,HPHA,N_name ! bessel, phase, totalpath
           ibn=0
        hpha=-pi
        n_name=0
          p=>my_ering%start
          do ii=1,my_ering%N
           if(p%mag%kind==kind4) then
            p%mag%c4%N_BESSEL=ibn
            p%magp%c4%N_BESSEL=ibn
            p%mag%c4%phase0=HPHA
            p%magp%c4%phase0=HPHA
            p%mag%c4%CAVITY_TOTALPATH=n_name
            p%magp%c4%CAVITY_TOTALPATH=n_name
           endif
           P=>P%NEXT
          ENDDO
       case('MODULATE','ACMAGNET')
          READ(MF,*) NAME , posr

          CALL CONTEXT(NAME)
          N_NAME=0
          IF(NAME(1:2)=='NO') THEN
             READ(MF,*) NAME
             call context(name)
             N_NAME=len_trim(name)
          ENDIF
          READ(MF,*)  DC_ac,A_ac,theta_ac
          READ(MF,*)  D_ac,n_ac !,n_coeff
          i2=0
          if(n_ac>0) then        !n_ac>0
             allocate(an(n_ac))
             allocate(bn(n_ac))
             an=0.0_dp
             bn=0.0_dp
             do while(.true.)
                read(mf,*) i1,dtu(1:2)
                if(i1<=0) exit
                if(i1>i2) i2=i1
                an(i1)=dtu(2)
                bn(i1)=dtu(1)
             enddo
          endif                  !n_ac>0
      !    if(n_coeff>0) then        !n_ac>0
      !       allocate(n_co(n_coeff))
      !       n_co=zero
      !       read(mf,*) n_co
      !    endif                  !n_ac>0
          p=>my_ering%start
          do ii=1,my_ering%N
             found_it=MY_FALSE
             IF(N_NAME==0) THEN
                found_it=P%MAG%NAME==NAME
             ELSE
                found_it=P%MAG%NAME(1:N_NAME)==NAME(1:N_NAME)
             ENDIF

             IF(FOUND_IT) THEN
                write(6,*) " slow ac magnet found ",P%MAG%name
                if(.not.associated(P%MAG%DC_ac)) then
                   allocate(P%MAG%DC_ac)
                   allocate(P%MAG%A_ac)
                   allocate(P%MAG%theta_ac)
                   allocate(P%MAG%D_ac)

                   allocate(P%MAGP%DC_ac)
                   allocate(P%MAGP%A_ac)
                   allocate(P%MAGP%theta_ac)
                   CALL alloc(P%MAGP%DC_ac)
                   CALL alloc(P%MAGP%A_ac)
                   CALL alloc(P%MAGP%theta_ac)
                   allocate(P%MAGp%D_ac)
                   CALL alloc(P%MAGP%D_ac)



                   P%MAG%D_ac=D_ac
                   P%MAG%DC_ac=DC_ac
                   P%MAG%A_ac=A_ac
                   P%MAG%theta_ac=theta_ac*twopi
                   P%MAGP%D_ac=D_ac
                   P%MAGP%DC_ac=DC_ac
                   P%MAGP%A_ac=A_ac
                   P%MAGP%theta_ac=theta_ac*twopi
                 if(P%MAG%slow_ac/=0) then
                  write(6,*) P%MAG%name, " already modulated "
                  stop 180
                 endif
                   P%MAG%slow_ac=posr
                   P%MAGP%slow_ac=posr

                   if(i2>p%mag%p%nmul) then
                      CALL ADD(P,i2,0,0.0_dp)
                   endif
                   allocate(P%MAG%d_an(p%mag%p%nmul))
                   allocate(P%MAG%d_bn(p%mag%p%nmul))
                   allocate(P%MAGp%d_an(p%mag%p%nmul))
                   allocate(P%MAGp%d_bn(p%mag%p%nmul))
                   allocate(P%MAG%d0_an(p%mag%p%nmul))
                   allocate(P%MAG%d0_bn(p%mag%p%nmul))
                   allocate(P%MAGp%d0_an(p%mag%p%nmul))
                   allocate(P%MAGp%d0_bn(p%mag%p%nmul))


                   P%MAG%d_an=0.0_dp
                   P%MAG%d_bn=0.0_dp

                   call alloc(P%MAGp%d_an,p%mag%p%nmul)
                   call alloc(P%MAGp%d_bn,p%mag%p%nmul)
                   call alloc(P%MAGp%d0_an,p%mag%p%nmul)
                   call alloc(P%MAGp%d0_bn,p%mag%p%nmul)
                   do i1=1,p%mag%p%nmul
                      P%MAG%d0_bn(i1)=P%MAG%bn(i1)
                      P%MAG%d0_an(i1)=P%MAG%an(i1)
                      P%MAGp%d0_bn(i1)=P%MAG%bn(i1)
                      P%MAGp%d0_an(i1)=P%MAG%an(i1)
                   enddo

                   do i1=1,n_ac
                      P%MAG%d_an(i1) =an(i1)
                      P%MAGp%d_an(i1)=an(i1)
                      P%MAG%d_bn(i1) =bn(i1)
                      P%MAGp%d_bn(i1)=bn(i1)
                   enddo
                   !
                   !    IF(K%MODULATION) THEN
                   !       DV=(XS%AC%X(1)*COS(EL%theta_ac)-XS%AC%X(2)*SIN(EL%theta_ac))
                   !       V=EL%DC_ac+EL%A_ac*DV
                   !       DV=el%D_ac*DV
                   !    else
                   !       V=EL%DC_ac
                   !       DV=zero
                   !    endif
                   !                   P%MAG%DC_ac=DC_ac
                   !                   P%MAG%A_ac=A_ac
                   !                   P%MAG%theta_ac=theta_ac*twopi
                   !                   P%MAGP%DC_ac=DC_ac
                   !                   P%MAGP%A_ac=A_ac
                   !                   P%MAGP%theta_ac=theta_ac*twopi
                   !                   P%MAG%slow_ac=.true.
                   !                   P%MAGP%slow_ac=.true.
                else
                   write(6,*) " already associated --> change code "
                   stop 166
                ENDIF
             ENDIF


             p=>p%next
          enddo


          if(n_ac>0) then
             deallocate(an,bn)
          endif

      !    if(n_coeff>0) then
      !       deallocate(n_co)
      !    endif
                 case('MODULATERF','ACMAGNETRF')
          READ(MF,*) NAME , posr

          CALL CONTEXT(NAME)
          N_NAME=0
          IF(NAME(1:2)=='NO') THEN
             READ(MF,*) NAME
             call context(name)
             N_NAME=len_trim(name)
          ENDIF
          READ(MF,*)  DC_ac,A_ac,theta_ac
          READ(MF,*)  D_ac,n_ac !,n_coeff
          i2=0
          if(n_ac>0) then        !n_ac>0
             allocate(an(n_ac))
             allocate(bn(n_ac))
             an=0.0_dp
             bn=0.0_dp
             do while(.true.)
                read(mf,*) i1,dtu(1:2)
                if(i1<=0) exit
                if(i1>i2) i2=i1
                an(i1)=dtu(2)
                bn(i1)=dtu(1)
             enddo
          else
              n_ac=1
              allocate(an(n_ac))
             allocate(bn(n_ac))
             an=0.0_dp
             bn=0.0_dp
             i2=1             
          endif                  !n_ac>0

          p=>my_ering%start
          do ii=1,my_ering%N
             found_it=MY_FALSE
             IF(N_NAME==0) THEN
                found_it=P%MAG%NAME==NAME
             ELSE
                found_it=P%MAG%NAME(1:N_NAME)==NAME(1:N_NAME)
             ENDIF

             IF(FOUND_IT) THEN
                write(6,*) " slow ac magnet found ",P%MAG%name
                 if(associated(p%mag%volt)) then
                    write(6,*) " It is an electric element "
                else
                  write(6,*) " Not electric element "
                  stop
                endif
                READ(MF,*) d_volt , d_phas
                if(.not.associated(P%MAG%DC_ac)) then
                   allocate(P%MAG%DC_ac)
                   allocate(P%MAG%A_ac)
                   allocate(P%MAG%theta_ac)
                   allocate(P%MAG%D_ac)

                   allocate(P%MAGP%DC_ac)
                   allocate(P%MAGP%A_ac)
                   allocate(P%MAGP%theta_ac)
                   CALL alloc(P%MAGP%DC_ac)
                   CALL alloc(P%MAGP%A_ac)
                   CALL alloc(P%MAGP%theta_ac)
                   allocate(P%MAGp%D_ac)
                   CALL alloc(P%MAGP%D_ac)



                   P%MAG%D_ac=D_ac
                   P%MAG%DC_ac=DC_ac
                   P%MAG%A_ac=A_ac
                   P%MAG%theta_ac=theta_ac*twopi
                   P%MAGP%D_ac=D_ac
                   P%MAGP%DC_ac=DC_ac
                   P%MAGP%A_ac=A_ac
                   P%MAGP%theta_ac=theta_ac*twopi
                 if(P%MAG%slow_ac/=0) then
                  write(6,*) P%MAG%name, " already modulated "
                  stop 180
                 endif
                   P%MAG%slow_ac=posr
                   P%MAGP%slow_ac=posr

                   if(i2>p%mag%p%nmul) then
                      CALL ADD(P,i2,0,0.0_dp)
                   endif
                              
                   allocate(P%MAG%d_an(p%mag%p%nmul))
                   allocate(P%MAG%d_bn(p%mag%p%nmul))
                   allocate(P%MAGp%d_an(p%mag%p%nmul))
                   allocate(P%MAGp%d_bn(p%mag%p%nmul))
                   allocate(P%MAG%d0_an(p%mag%p%nmul))
                   allocate(P%MAG%d0_bn(p%mag%p%nmul))
                   allocate(P%MAGp%d0_an(p%mag%p%nmul))
                   allocate(P%MAGp%d0_bn(p%mag%p%nmul))


                   P%MAG%d_an=0.0_dp
                   P%MAG%d_bn=0.0_dp

                   call alloc(P%MAGp%d_an,p%mag%p%nmul)
                   call alloc(P%MAGp%d_bn,p%mag%p%nmul)
                   call alloc(P%MAGp%d0_an,p%mag%p%nmul)
                   call alloc(P%MAGp%d0_bn,p%mag%p%nmul)
                   do i1=1,p%mag%p%nmul
                      P%MAG%d0_bn(i1)=P%MAG%bn(i1)
                      P%MAG%d0_an(i1)=P%MAG%an(i1)
                      P%MAGp%d0_bn(i1)=P%MAG%bn(i1)
                      P%MAGp%d0_an(i1)=P%MAG%an(i1)
                   enddo

                   do i1=1,n_ac
                      P%MAG%d_an(i1) =an(i1)
                      P%MAGp%d_an(i1)=an(i1)
                      P%MAG%d_bn(i1) =bn(i1)
                      P%MAGp%d_bn(i1)=bn(i1)
                   enddo
                                
                            allocate(P%MAG%D0_Volt)         
                            allocate(P%MAGp%D0_Volt)         
                            allocate(P%MAG%D_Volt)         
                            allocate(P%MAGp%D_Volt) 
                            allocate(P%MAG%D0_phas)         
                            allocate(P%MAGp%D0_phas)         
                            allocate(P%MAG%D_phas)         
                            allocate(P%MAGp%D_phas) 
         
                       call alloc(P%MAGp%D0_Volt)             
                       call alloc(P%MAGp%D_Volt) 
                       call alloc(P%MAGp%D0_phas)             
                       call alloc(P%MAGp%D_phas) 
                       
                         P%MAG%D0_Volt= P%MAG%Volt
                         P%MAGp%D0_Volt= P%MAG%Volt
                         P%MAG%D_Volt= d_volt
                         P%MAGp%D_Volt= d_volt
                         P%MAG%D0_phas= P%MAG%phas
                         P%MAGp%D0_phas= P%MAG%phas
                         P%MAG%D_phas= d_phas
                         P%MAGp%D_phas= d_phas
             
                else
                   write(6,*) " already associated --> change code "
                   stop 166
                ENDIF
             ENDIF


             p=>p%next
          enddo


          if(n_ac>0) then
             deallocate(an,bn)
          endif

      !    if(n_coeff>0) then
      !       deallocate(n_co)
      !    endif

       case('MAPSFORZHE')
          READ(MF,*) i11,I22,number_zhe_maps_local ,hgap ! position  i1=i2 one turn map,  fact is precision of stochastic kick
          READ(MF,*) MY_A_NO,do_state0   ! ORDER OF THE MAP  
          READ(MF,*) filename
          state0=my_estate-radiation0-envelope0
          if(do_state0) then
           do_state0=my_estate%radiation.or.my_estate%envelope
          endif
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          if(i11==i22) i22=i11+my_ering%n
          di12=float(i22-i11)/number_zhe_maps_local
                        do k_zhe=1,number_zhe_maps_local
                             write(name_zhe,*) k_zhe
                           file_zher(1:len_trim(filename))=filename(1:len_trim(filename))
                           file_zher(1+len_trim(filename):len_trim(filename)+len_trim(name_zhe))=name_zhe(1:len_trim(name_zhe))
                           call context(file_zher)
                          write(6,*) "out file ",file_zher(1:len_trim(file_zher))
           
          i1=i11+(k_zhe-1)*di12
          i2=i1+di12
          if(i2>i22.or.k_zhe==number_zhe_maps_local) i2=i22
          write(6,*)" from to ", i1,i2,i22
!pause 873
           p=>my_ering%start
           f1=>p          
           do ii=2,i1
            p=>p%next
             f1=>p
           enddo
           if(I2==i1) then
            f2=>f1
           else
            p=>my_ering%start
            f2=>p 
            do ii=2,i2
              p=>p%next
              f2=>p
            enddo
            endif
             x_ref=0.0_dp


if(.not.my_estate%envelope) hgap=-1
            call FIND_ORBIT_x(x_ref,my_estate,1.d-7,fibre1=f1)
 if(do_state0)            call FIND_ORBIT_x(x_ref0,state0,1.d-7,fibre1=f1)

 if(do_state0)   then
              call fill_tree_element_line_zhe0(state0,my_estate,f1,f2,MY_A_NO,x_ref0,x_ref,file_zher,stochprec=hgap) 

 else
              call fill_tree_element_line_zhe(my_estate,f1,f2,MY_A_NO,x_ref,file_zher,stochprec=hgap) 
 endif

write(6,*) " State used "
          call print(my_estate,6)
write(6,*) " closed orbit at position ",i1
           write(6,*) x_ref(1:3)
           write(6,*) x_ref(4:6)
                  enddo

       case('MAPFORZHE')
          READ(MF,*) i1,I2,hgap  ! position  i1=i2 one turn map,  fact is precision of stochastic kick
          READ(MF,*) MY_A_NO   ! ORDER OF THE MAP  
          READ(MF,*) filename
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          
           p=>my_ering%start
           f1=>p          
           do ii=2,i1
            p=>p%next
             f1=>p
           enddo
           if(I2==i1) then
            f2=>f1
           else
            p=>my_ering%start
            f2=>p 
            do ii=2,i2
              p=>p%next
              f2=>p
            enddo
            endif
             x_ref=0.0_dp


if(.not.my_estate%envelope) hgap=-1
 if(my_a_no>0)            call FIND_ORBIT_x(x_ref,my_estate,1.d-7,fibre1=f1)

              call fill_tree_element_line_zhe(my_estate,f1,f2,iabs(MY_A_NO),x_ref,filename,stochprec=hgap) 


write(6,*) " State used "
          call print(my_estate,6)
write(6,*) " closed orbit at position ",i1
           write(6,*) x_ref(1:3)
           write(6,*) x_ref(4:6)
 

       case('MAKEMAPITOJ')
          READ(MF,*) i1,I2,i3  ! position
          READ(MF,*) MY_A_NO  ! ORDER OF THE MAP
          READ(MF,*) fixp,fact,noca  !  SYMPLECTIC , factored
          READ(MF,*) filename
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          
           
           p=>my_ering%start
           f1=>p          
           do ii=2,i1
            p=>p%next
             f1=>p
           enddo

           p=>my_ering%start
           f2=>p 
           do ii=2,i2
             p=>p%next
             f2=>p
           enddo
 
             x_ref=0.0_dp
             call MOVE_TO_LAYOUT_I(m_u,my_fring,i3)

             if(associated(my_ering,my_fring))then 
                ft=>f1
             else
                ft=>my_fring%start       
             endif

if(noca) then
             call FIND_ORBIT_x(x_ref,time0+nocavity0,1.d-7,fibre1=f1)
else
             call FIND_ORBIT_x(x_ref,time0,1.d-7,fibre1=f1)
endif  

             name_root=filename
             call context(name_root)
             if(name_root(1:2)=='NO') then
              call fill_tree_element_line(f1,f2,ft,i2,x_ref,fact,nocav=noca)
             else
              call fill_tree_element_line(f1,f2,ft,i2,x_ref,fact,nocav=noca,file=filename)
             endif
  


                    ft%mag%forward(3)%symptrack=FIXP
                    ft%magP%forward(3)%symptrack=FIXP
                    ft%mag%do1mapf=.true.
                    ft%magp%do1mapf=.true.
                    ft%mag%filef="one_turn_map.txt"
             if(associated(ft,f1))then
                p=>f1%next
               do i1=1,my_ering%n
                if(associated(f2,p)) exit
                 p%mag%skip_ptc_f=1
                 p%magp%skip_ptc_f=1
                p=>p%next 
               enddo
                
             endif 

       case('MAKEONETURNMAP','TRACKWITHONETURNMAP')
          READ(MF,*) i1,i3  ! position
          READ(MF,*) I2  ! ORDER OF THE MAP
          READ(MF,*) fixp,fact,noca  !  SYMPLECTIC , factored
          READ(MF,*) filename
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          
           p=>my_ering%start
           f1=>p          
           do ii=2,i1
            p=>p%next
             f1=>p
           enddo


             f2=>f1
             x_ref=0.0_dp
             call MOVE_TO_LAYOUT_I(m_u,my_fring,i3)

             if(associated(my_ering,my_fring))then 
                ft=>f1
             else
                ft=>my_fring%start       
             endif
if(noca) then
             call FIND_ORBIT_x(x_ref,time0+nocavity0,1.d-7,fibre1=f1)
else
             call FIND_ORBIT_x(x_ref,time0,1.d-7,fibre1=f1)
endif  
write(6,*) x_ref
       
             name_root=filename
             call context(name_root)
             if(name_root(1:2)=='NO') then
              call fill_tree_element_line(f1,f2,ft,i2,x_ref,fact,nocav=noca)
             else
              call fill_tree_element_line(f1,f2,ft,i2,x_ref,fact,nocav=noca,file=filename)
             endif
                    ft%mag%forward(3)%symptrack=FIXP
                    ft%magP%forward(3)%symptrack=FIXP
                    ft%mag%do1mapf=.true.
                    ft%magp%do1mapf=.true.
                    ft%mag%filef="one_turn_map.txt"
             if(associated(ft,f1))then
                f2=>f1%next
               do i1=1,my_ering%n
                if(associated(f2,f1)) exit
                 f2%mag%skip_ptc_f=1
                 f2%magp%skip_ptc_f=1
                f2=>f2%next 
               enddo
                
             endif 


       case('MAKEMAP','TRACKWITHMAP')
          READ(MF,*) NAME
          READ(MF,*) I1  ! ORDER OF THE MAP
          READ(MF,*)  onemap  ! use one map : no cutting
          READ(MF,*) fixp,fact  !  SYMPLECTIC , factored
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          x_ref=0.0_DP

          READ(MF,*) x_ref
          n_ac=0
          CALL CONTEXT(NAME)
          N_NAME=0
          IF(NAME(1:2)=='NO') THEN
             READ(MF,*) NAME
             call context(name)
             N_NAME=len_trim(name)
          ENDIF

          p=>my_ering%start
          do ii=1,my_ering%N
             found_it=MY_FALSE
             IF(N_NAME==0) THEN
                found_it=P%MAG%NAME==NAME
             ELSE
                found_it=P%MAG%NAME(1:N_NAME)==NAME(1:N_NAME)
             ENDIF

             IF(FOUND_IT) THEN
                write(6,*) "  magnet found FOR MAP REPLACEMENT ",P%MAG%name
                call fill_tree_element(p,I1,x_REF,onemap,fact)
                 n_ac=n_ac+1

                   IF(P%DIR==1) THEN
                    p%mag%forward(3)%symptrack=FIXP
                    p%magP%forward(3)%symptrack=FIXP
                do i2=1,3
                 write(p%mag%filef,*) "map",n_ac,".txt"
                 call context(p%mag%filef)
                enddo
                   ELSE
                do i2=1,3
                 write(p%mag%fileb,*) "map",n_ac,".txt"
                 call context(p%mag%fileb)
                enddo
                    p%mag%BACKward(3)%symptrack=FIXP
                    p%magP%BACKward(3)%symptrack=FIXP
                   ENDIF
             ENDIF

             p=>p%next
          enddo

       case('MAKEALLMAP','TRACKALLWITHMAP')
          READ(MF,*) I1  ! ORDER OF THE MAP
          READ(MF,*)  onemap  ! use one map : no cutting
          READ(MF,*) fixp,fact  !  SYMPLECTIC 
          READ(MF,*) skipcav  !  skip cavity 
          icnmin=0
          icnmax=0
          icavv=0
          x_ref=0.0_DP
           n_ac=0
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          p=>my_ering%start
          do ii=1,my_ering%N
   

             IF(p%mag%kind/=kind0) THEN
             icnmax=icnmax+1
              if(.not.skipcav.or.(p%mag%kind/=kind4.and.p%mag%kind/=kind21)) then
              !  write(6,*) "  magnet found FOR MAP REPLACEMENT ",P%MAG%name
                 n_ac=n_ac+1
                call fill_tree_element(p,I1,x_REF,onemap,fact)
                   IF(P%DIR==1) THEN
                    p%mag%forward(3)%symptrack=FIXP
                    p%magP%forward(3)%symptrack=FIXP
                do i2=1,3
                 write(p%mag%filef,*) "map",n_ac,".txt"
                 call context(p%mag%filef)
 
                enddo
                   ELSE
                do i2=1,3
                 write(p%mag%fileb,*) "map",n_ac,".txt"
                 call context(p%mag%fileb)

                enddo
                    p%mag%BACKward(3)%symptrack=FIXP
                    p%magP%BACKward(3)%symptrack=FIXP
                   ENDIF
              else
                icavv=icavv+1
              endif
             else
              icnmin=icnmin+1
             ENDIF


             p=>p%next
          enddo
         write(6,*) icnmax, " changed into Taylor maps "
         write(6,*) icnmin, " markers "
         write(6,*) icavv, " cavities left alone "
         write(6,*) my_ering%N, " total number of fibres "
       case('REMOVEALLMAPS')
       READ(MF,*) I1,I2  ! ORDER OF THE MAP

          p=>my_ering%start
          do ii=i1,i2-1
             p%MAG%skip_ptc_f=-p%MAG%skip_ptc_f
             p%MAGp%skip_ptc_f=-p%MAGp%skip_ptc_f
             p=>p%next
          enddo


          p=>my_ering%start
          do ii=1,my_ering%N
   
          if(associated(p%MAG%forward)) then
             p%MAG%usef=.false.
             p%MAGp%usef=.false.
          endif
           if(associated(p%MAG%backward)) then
             p%MAG%useb=.false.
             p%MAGp%useb=.false.
           endif
             p=>p%next
          enddo
    
       case('PUTBACKALLMAP')

       READ(MF,*) I1,I2  ! ORDER OF THE MAP

          p=>my_ering%start
          do ii=i1,i2-1
             p%MAG%skip_ptc_f=-p%MAG%skip_ptc_f
             p%MAGp%skip_ptc_f=-p%MAGp%skip_ptc_f
             p=>p%next
          enddo


          p=>my_ering%start
          do ii=1,my_ering%N
   
       if(associated(p%MAG%forward)) then
             p%MAG%usef=.true.
             p%MAGp%usef=.true.
       endif
       if(associated(p%MAG%backward)) then
             p%MAG%useb=.true.
             p%MAGp%useb=.true.
       endif
 
             p=>p%next
          enddo
    
       case('SETAPERTURE')
          READ(MF,*) KINDA,NAME

          CALL CONTEXT(NAME)
          N_NAME=0
          IF(NAME(1:2)=='NO') THEN
             READ(MF,*) NAME
             call context(name)
             N_NAME=len_trim(name)
          ENDIF
          READ(MF,*) XA,YA,DXA,DYA
          READ(MF,*) RA(1),RA(2)

          p=>my_ering%start
          do ii=1,my_ering%N
             found_it=MY_FALSE
             IF(N_NAME==0) THEN
                found_it=P%MAG%NAME==NAME
             ELSE
                found_it=P%MAG%NAME(1:N_NAME)==NAME(1:N_NAME)
             ENDIF

             IF(FOUND_IT) THEN
                IF(.NOT.ASSOCIATED(P%MAG%P%APERTURE)) THEN
                   CALL alloc(P%MAG%P%APERTURE)
                ENDIF
                IF(.NOT.ASSOCIATED(P%MAGP%P%APERTURE)) THEN
                   CALL alloc(P%MAGP%P%APERTURE)
                ENDIF

                P%MAG%P%APERTURE%KIND= KINDA  ! 1,2,3,4
                P%MAG%P%APERTURE%R= RA
                P%MAG%P%APERTURE%X= XA
                P%MAG%P%APERTURE%Y= YA
                P%MAG%P%APERTURE%DX= DXA
                P%MAG%P%APERTURE%DY= DYA
                P%MAGP%P%APERTURE%KIND= KINDA  ! 1,2,3,4
                P%MAGP%P%APERTURE%R= RA
                P%MAGP%P%APERTURE%X= XA
                P%MAGP%P%APERTURE%Y= YA
                P%MAGP%P%APERTURE%DX= DXA
                P%MAGP%P%APERTURE%DY= DYA

             ENDIF

             p=>p%next
          enddo




       case('PUTFAMILY','LAYOUT<=KNOB')
          do j=1,NPOL
             my_ering=pol_(j)
          enddo
       case('DEALLOCATEFAMILIES')
          call kill_para(my_ering)
          deallocate(POL_)
       case('PTCTWISSTENGEDWARDS','TWISSTENGEDWARDS')  !
          read(mf,*) filename, NAME, integrated
          read(mf,*) del

    !      call compute_twiss(my_ering,my_estate,filename,1,del,1,integrated,name,my_true,my_false)
       case('PTCTWISS','TWISS','PTCTWISSRIPKEN','TWISSRIPKEN')  !
          read(mf,*) filename, NAME, integrated
          read(mf,*) del

       !   call compute_twiss(my_ering,my_estate,filename,1,del,1,integrated,name,my_false,my_false)

       case('PTCTWISSSASHA','TWISSSASHA','PTCTWISSRIPKENSASHA','TWISSRIPKENSASHA')  !
          read(mf,*) filename, NAME, integrated
          read(mf,*) del

    !      call compute_twiss(my_ering,my_estate,filename,1,del,1,integrated,name,my_false,my_true)

       case('FITTUNESCAN','SEARCHAPERTUREX=Y')
          read(mf,*) epsf
          read(mf,*) ntune
          read(mf,*) tune(1:2)
          read(mf,*) dtune(1:2),targ_RES(3:4)
          if(com=='SEARCHAPERTUREX=Y') then
             READ(MF,*) r_in,del_in,dx,DLAM,fixp
             READ(MF,*) POS,NTURN,ITE,FILENAME,name
             call context(name)
             if(name(1:11)/='NONAMEGIVEN') then
                posr=pos
                call move_to( my_ering,p,name,posR,POS)
                if(pos==0) then
                   write(6,*) name, " not found "
                   stop
                endif
             endif
             call kanalnummer(mfr,filename)
          endif




          do i2=0,ntune(2)
             do i1=0,ntune(1)

                if(ntune(1)/=0) then
                   targ_tune(1)=tune(1)+((dtune(1)-tune(1))*i1)/(ntune(1))
                else
                   targ_tune(1)=tune(1)
                endif
                if(ntune(2)/=0) then
                   targ_tune(2)=tune(2)+((dtune(2)-tune(2))*i2)/(ntune(2))
                else
                   targ_tune(2)=tune(2)
                endif
                if(abs(targ_RES(3))>999) then
                   call lattice_fit_TUNE_gmap(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)
                else
                   targ_RES(1:2)=targ_tune(1:2)
                   call lattice_fit_tune_CHROM_gmap(my_ering,my_estate,EPSF,pol_,NPOL,targ_RES,NP)
                endif
                if(com=='SEARCHAPERTUREX=Y') then
                   ! write(mfr,*) targ_tune(1:2)
                   CALL dyn_aperalex(my_ering,r_in,del_in,dx,dlam,pos,nturn,ite,my_estate,MFR,targ_tune,fixp)
                endif
             enddo
          enddo
          if(com=='SEARCHAPERTUREX=Y') close(mfr)

       case('ALEXFITTUNESCAN','ALEXSEARCHAPERTUREX=Y')
          read(mf,*) epsf
          read(mf,*) ntune
          read(mf,*) tune(1:2)
          read(mf,*) dtune(1:2),targ_RES(3:4)
          if(com=='ALEXSEARCHAPERTUREX=Y') then
             READ(MF,*) r_in,del_in,dx,DLAM,fixp
             READ(MF,*) POS,NTURN,ITE,FILENAME,name
             call context(name)
             if(name(1:11)/='NONAMEGIVEN') then
                posr=pos
                call move_to( my_ering,p,name,posR,POS)
                if(pos==0) then
                   write(6,*) name, " not found "
                   stop
                endif
             endif
             call kanalnummer(mfr,filename)
          endif

          sc=1.d0


          do i2=0,ntune(2)
             do i1=0,ntune(1)

                if(ntune(1)/=0) then
                   targ_tune(1)=tune(1)+((dtune(1)-tune(1))*i1)/(ntune(1))
                else
                   targ_tune(1)=tune(1)
                endif
                if(ntune(2)/=0) then
                   targ_tune(2)=tune(2)+((dtune(2)-tune(2))*i2)/(ntune(2))
                else
                   targ_tune(2)=tune(2)
                endif
                targ_tune_alex(1)=22.0_dp+targ_tune(1)
                targ_tune_alex(2)=20.0_dp+targ_tune(2)

                if(abs(targ_RES(3))>999) then

                   CALL special_alex_main_ring_auto(my_ering,3,targ_tune_alex,sc,epsf)
                   call lattice_fit_TUNE_gmap(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)
                else
                   CALL special_alex_main_ring_auto(my_ering,3,targ_tune_alex,sc,epsf)
                   targ_RES(1:2)=targ_tune(1:2)
                   call lattice_fit_tune_CHROM_gmap(my_ering,my_estate,EPSF,pol_,NPOL,targ_RES,NP)
                endif
                if(com=='ALEXSEARCHAPERTUREX=Y') then
                   ! write(mfr,*) targ_tune(1:2)
                   CALL dyn_aperalex(my_ering,r_in,del_in,dx,dlam,pos,nturn,ite,my_estate,MFR,targ_tune,fixp)
                endif
             enddo
          enddo
          if(com=='ALEXSEARCHAPERTUREX=Y') close(mfr)

       case('FITTUNE')
          read(mf,*) epsf
          read(mf,*) targ_tune
          if(targ_tune(1)<=0.0_dp) targ_tune=tune(1:2)
          call lattice_fit_TUNE_gmap(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)

       case('FITTUNERAD')
          read(mf,*) epsf
          read(mf,*) targ_tune
          if(targ_tune(1)<=0.0_dp) targ_tune=tune(1:2)
          call lattice_fit_TUNE_gmap_rad(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)

       case('DELTAFITTUNE')
          read(mf,*) epsf
          read(mf,*) targ_tune
          tempstate=my_estate+nocavity0
          call lattice_GET_tune(my_ering,tempstate,mftune,tune)
          targ_tune(1:2)=targ_tune(1:2)+tune(1:2) 
          call lattice_fit_TUNE_gmap(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)
       case('FITTUNEAUTO')
          read(mf,*) epsf
          read(mf,*) targ_tune
          read(mf,*) namet(1), namet(2)
          if(targ_tune(1)<=0.0_dp) targ_tune=tune(1:2)
          call lattice_fit_TUNE_gmap_auto(my_ering,my_estate,EPSF,targ_tune,namet)
       case('SCANTUNE')
          STRAIGHT=.FALSE.
          read(mf,*) epsf
          read(mf,*) nstep
          read(mf,*) tune_ini,tune_fin
          read(mf,*) name_root,SUFFIX
          dtu=0.0_dp
          IF(NSTEP(2)/=0) THEN
             if(nstep(1)/=1) dtu(1)=(tune_fin(1)-tune_ini(1))/(nstep(1)-1)
             if(nstep(2)/=1) dtu(2)=(tune_fin(2)-tune_ini(2))/(nstep(2)-1)
          ELSE
             dtu(1)=(tune_fin(1)-tune_ini(1))/(nstep(1)-1)
             dtu(2)=(tune_fin(2)-tune_ini(2))/(nstep(1)-1)
             STRAIGHT=.TRUE.
             NSTEP(2)=1
          ENDIF

          I3=0
          do i1=0,nstep(1)-1
             do i2=0,nstep(2)-1
                IF(STRAIGHT) THEN
                   targ_tune(1)=tune_ini(1)+dtu(1)*i1
                   targ_tune(2)=tune_ini(2)+dtu(2)*I1
                ELSE
                   targ_tune(1)=tune_ini(1)+dtu(1)*i1
                   targ_tune(2)=tune_ini(2)+dtu(2)*i2
                ENDIF
                call lattice_fit_TUNE_gmap(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)
                write(title,*) 'tunes =',targ_tune
                I3=I3+1
                call create_name(name_root,i3,suffix,filename)
                write(6,*)" printing in "
                write(6,*)filename
                call print_bn_an(my_ering,2,title,filename)
                write(6,*)" PRINTED "
             enddo
          enddo
       case('FITSEX')
          read(mf,*) epsf
          read(mf,*) targ_chrom
          read(mf,*)i1
          call c_lattice_fit_CHROM_gmap1(my_ering,my_estate,EPSF,pol_,NPOL,targ_chrom,np,i1,mf)
       case('FITSEXLINEAR')
          read(mf,*) epsf
          read(mf,*) targ_chrom
          read(mf,*)i1
          call lattice_fit_CHROM_gmap2(my_ering,my_estate,EPSF,pol_,NPOL,targ_chrom,np,i1,mf)
       case('FITCHROMATICITY')
          read(mf,*) epsf
          read(mf,*) targ_chrom
          call lattice_fit_CHROM_gmap(my_ering,my_estate,EPSF,pol_,NPOL,targ_chrom,NP)
       case('FITTUNECHROMATICITY')
          read(mf,*) epsf
          read(mf,*) targ_RES
          call lattice_fit_tune_CHROM_gmap(my_ering,my_estate,EPSF,pol_,NPOL,targ_RES,NP)
       case('GETCHROMATICITY')
          call lattice_GET_CHROM(my_ering,my_estate,CHROM)
       case('OPENTUNEFILE')
          read(mf,*) filename
          call kanalnummer(mftune,filename)
!          write(mftune,*) " Time unit = ",unit_time ," seconds "
       case('CLOSETUNEFILE')
           close(mftune)
           mftune=6
       case('GETTUNE')
          call lattice_GET_tune(my_ering,my_estate,mftune,tune)
       case('STRENGTH','STRENGTHFILE')  !
          IF(mfpolbloc/=0) CLOSE(mfpolbloc)

          READ(MF,*) file_block_name

          if(file_block_name(1:7)/="noprint") then
             call kanalnummer(mfpolbloc)
             open(unit=mfpolbloc,file=file_block_name)
          endif

          WRITE(6,*) " KNOBS PRINTED IN ",file_block_name(1:LEN_TRIM(file_block_name))
       case('NOSTRENGTH','NOSTRENGTHFILE')  !
          file_block_name='noprint'
          IF(mfpolbloc/=0) CLOSE(mfpolbloc)
       case('FINALSETTING')  ! ACCELERATION FILE
          READ(MF,*) FINAL_setting
2000      INQUIRE (FILE = FINAL_setting, EXIST = exists)
          if(exists) then
             write(6,*) "file ", FINAL_setting(1:len_trim(FINAL_setting)), &
                  " exists, interrupt execution if you do not want to overwrite!"
             !     write(6,*) " you have 2 seconds to do so "
             !     CALL DATE_AND_TIME(values=temps)
             !     i1=temps(7)
             !      if(i1>=58) i1=i1-58
             !     do while(.true.)
             !      CALL DATE_AND_TIME(values=temps)
             !      if(temps(7)>i1+2) exit
             !     enddo
          endif
       case('INITIALSETTING') ! ACCELERATION FILE
          READ(MF,*) initial_setting
2001      INQUIRE (FILE = initial_setting, EXIST = exists)
          if(.not.exists) then
             write(6,*) "file ", initial_setting(1:len_trim(initial_setting)), &
                  " does not exist, please input now on screen "
             read(5,*) initial_setting
             GOTO 2001
          endif
       case('PAUSE')
          WRITE(6,*) " Type enter to continue execution "
          READ(5,*)
       case('POLARIZATION')
         read(mf,*) nturns, ngen,no_pol
         read(mf,*) n_ang,nlm,lm
         read(mf,*) snake,isnake
         read(mf,*) filename_ap,datafile
        call compute_polarisation(my_ering,no_pol,n_ang,nlm,lm,nturns,ngen,snake,isnake,filename_ap,datafile)
       case('PRINTONCE')
          print77=.true.
          read77=.true.
       CASE('4DMAP')
          READ(MF,*) NO
          READ(MF,*) POS, DEL
          READ(MF,*) filename
          CALL compute_map_4d(my_ering,my_estate,filename,pos,del,no)
       case('PRINTAINRESONANCE')
          READ(MF,*) MRES(1:4)
          READ(MF,*) NO,EMIT0
          READ(MF,*) FILENAME
          CALL lattice_PRINT_RES_FROM_A(my_ering,my_estate,NO,EMIT0,MRES,FILENAME)
       case('PRINTTWICE')
          print77=.false.
          read77=.false.
       case('PRINTBNAN','PRINTANBN','PRINTBN','PRINTAN')
          read(mf,*) title
          read(mf,*) nmul,filename
          call print_bn_an(my_ering,nmul,title,filename)
       case('READBNAN','READANBN','READBN','READAN')
          read(mf,*) filename
          call READ_bn_an(my_ering,filename)
       case('RETURN','EXIT','QUIT','END','STOP')
          goto 100
       case('PTCEND','ENDPTC','APOCALYSPE')
          CALL PTC_END()
          goto 100
       case('TRACK4DNORMALIZED')
          emit=0.0_dp
          read(mf,*) IB
          read(mf,*) POS,NTURN,ITMAX,resmax
          read(mf,*) EMIT(1:2),APER(1:2),emit(3),emit(6)
          read(mf,*) filename,filetune,FILESMEAR
          ! emit(3)=1.d38
          CALL track_aperture(my_ering,my_estate,beta,dbeta,tune,ib,ITMAX,emit,aper,pos,nturn,FILENAME,filetune,FILESMEAR,resmax)
       case('SCANTRACK4DNORMALIZED')
          emit=0.0_dp

          read(mf,*) POS,NTURN,ITMAX,resmax
          read(mf,*) EMIT0(1:2),APER,emit(3),emit(6)
          read(mf,*) nscan,name_root,SUFFIX
          read(mf,*) name_root_res,SUFFIX_res
          ib=1
          allocate(resu(nscan,4))
          do i1=1,nscan
             write(6,*) " CASE NUMBER ",i1
             call create_name(name_root,i1,suffix,filename)
             call READ_bn_an(my_ering,filename)
             call create_name(name_root_res,i1,SUFFIX_res,filename)
             COMT=name_root_res(1:len_trim(name_root_res))//"_resonance"
             call create_name(COMT,i1,SUFFIX_res,filetune)
             COMT=name_root_res(1:len_trim(name_root_res))//"_SMEAR"
             call create_name(COMT,i1,SUFFIX_res,FILESMEAR)
             if(i1/=1) ib=2
             emit(1:2)=EMIT0
             CALL track_aperture(my_ering,my_estate,beta,dbeta,tune,ib,ITMAX,emit,aper,pos,nturn,FILENAME,filetune,FILESMEAR,resmax)
             resu(i1,1:2)=emit(4:5)
             resu(i1,3:4)=emit(1:2)
          enddo
          filetune=name_root_res(1:len_trim(name_root_res))//'.'//SUFFIX_res
          call kanalnummer(mfr)

          OPEN(UNIT=mfr,FILE=filetune)
          write(mfr,*) " Precision = ",emit(3),emit(6)
          do i1=1,nscan
             write(mfr,205) i1,resu(i1,1:4)
          enddo
          deallocate(resu)
          close(mfr)
205       FORMAT(1x,i4,4(1X,D18.11))
       case('SEARCHAPERTURE')

          READ(MF,*) r_in,n_in,ang_in,ang_out,del_in,DLAM
          READ(MF,*) POS,NTURN,ITE,FILENAME,name
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ering,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif


          call kanalnummer(mfr)
          open(unit=mfr,file=filename)
          CALL dyn_aper(my_ering,r_in,n_in,ang_in,ang_out,del_in,dlam,pos,nturn,ite,my_estate,MFR)
          close(mfr)

       case('KNOB')

          READ(MF,*) POS, IBN

          if (pos > my_ering%n) stop

          p=>my_ering%start
          do ii=1,pos
             p=>p%next
          enddo


          write(6,*) "El name ", p%mag%name


          CALL INIT(default,3,1,BERZ)

          print*, "Npara is ", c_%NPARA

          pb = 0
          pb%name = p%mag%name
          write(6,*) "IBN ", IBN
          pb%ibn(ibn) = 1

          my_ering = pb

          CALL ALLOC(ID)
          CALL ALLOC(Y)
          x(:)=0
          ID=1
          Y=X+ID

          p=>my_ering%start
          do ii=1,my_ering%n
             write(6,*) "##########################################"
             write(6,'(i4, 1x,a, f10.6)') ii,p%mag%name
             write(6,'(a, f9.6, a)') "Ref Momentum ",p%mag%p%p0c," GeV/c"

             call track(my_ering,y,ii,ii+1,default)
             call daprint(y(1),6)
             p=>p%next

          enddo

       case('ALEXREMOVAL')
          call special_alex_main_ring_removal(my_ering)
       case('SASHASPECIALRCS')
          READ(MF,*) epsf,sca   ! aper scale >0 <=1
      ! call lattice_fit_bump_rcs(my_ering,epsf)
        call lattice_fit_bump_min_rcs(my_ering%next,my_ering,EPSF,pol_,NPOL,sca)
       case('PRINTFRAMES')

          READ(MF,*) FILENAME
          CALL print_frames(my_ering,filename)


       case('PRINTNEWFLATFILE')

          READ(MF,*) FILENAME

          call print_new_flat(my_ering,filename)

       case('READNEWFLATFILE')

          READ(MF,*) FILENAME
          call read_lattice_append(M_U,filename)
          WRITE(6,*) M_U%END%N, M_U%END%END%POS

       case('READFLATFILE')

          READ(MF,*) FILENAME
          CALL  READ_AND_APPEND_VIRGIN_general(M_U,filename)

          WRITE(6,*) M_U%END%N, M_U%END%END%POS

       case('PRINTFLATFILE')

          READ(MF,*) FILENAME
          CALL  print_COMPLEX_SINGLE_STRUCTURE(my_ering,filename,lmax0=lmax)

          WRITE(6,*) M_U%END%N, M_U%END%END%POS
       case('SKIPMARKER','SKIPMARKERS')
        print_marker=my_false
       case('INCLUDEMARKER','INCLUDEMARKERS')
        print_marker=my_true
        
       case('TOGGLEMARKER','TOGGLEMARKERS')
        print_marker=.not.print_marker
        if(print_marker) then
             Write(6,*)  'printing makers on flat file '
        else
            Write(6,*)  ' NOT printing makers on flat file '
        endif
        case('PERMFRINGEON')

          CALL  PUTFRINGE(my_ering,MY_TRUE)

       case('PERMFRINGEOFF')

          CALL  PUTFRINGE(my_ering,MY_FALSE)

       case('BENDFRINGEON')

          CALL  PUTbend_FRINGE(my_ering,MY_TRUE)

       case('BENDFRINGEOFF')

          CALL  PUTbend_FRINGE(my_ering,MY_FALSE)

       case('SOLENOID2','SOLENOIDFRINGE')

          write(6,*) " At the ends of all solenoids kick of the form "
          write(6,*) " X(2)=X(2)+X(1)*B/PZ      where pz=1+delta"
          write(6,*) " X(4)=X(4)+X(3)*B/PZ "
          write(6,*) " X(6)=X(6)-0.5_dp*B*(X(1)**2+X(3)**2)*TIME_FAC/PZ**2 "

          read(mf,*) fint,hgap

          p=>my_ering%start
          nscan=0

          do ii=1,my_ering%n
             if(associated(p%mag%s5)) then
                p%mag%s5%fint=fint
                p%mag%s5%hgap=hgap
                p%magp%s5%fint=fint
                p%magp%s5%hgap=hgap
                nscan=nscan+1
             endif
             p=>p%next
          enddo
          Write(6,*) " Found ",nscan," Solenoids"

       case('CAVITYBESSEL','CAVITYINNERFIELD')

          read(mf,*) n_bessel
          nscan=0
          p=>my_ering%start

          do ii=1,my_ering%n
             if(associated(p%mag%c4)) then
                p%mag%c4%n_bessel=n_bessel
                p%magp%c4%n_bessel=n_bessel
                nscan=nscan+1
             endif
             p=>p%next
          enddo

          Write(6,*) " Found ",nscan," Cavities"

       case('REVERSEBEAMLINE')

          CALL  REVERSE_BEAM_LINE(my_ering)

          !         WRITE(6,*) M_U%END%N, M_U%END%END%POS

       case('PSREXAMPLEOFPATCHING')

          call APPEND_EMPTY_LAYOUT(m_u)
          CALL remove_drifts(my_ering,m_u%END)
          m_u%end%name="psr_no_drift"
          call APPEND_EMPTY_LAYOUT(m_u)
          m_u%end%name="psr_quads_for_bends"
          CALL remove_drifts_bends(my_ering,m_u%END)

          WRITE(6,*) my_ering%N , m_u%END%N

       case('NORMALFORM')
          READ(MF,*)POS,name
          READ(MF,*) FILENAME

       case('TRANSLATELAYOUT')
          READ(MF,*)DT
          CALL TRANSLATE(my_ering,DT)
       case('TRANSLATEPARTOFLAYOUT')
          READ(MF,*)DT
          READ(MF,*) I1,I2
          CALL TRANSLATE(my_ering,DT,I1,I2)
          CALL MOVE_TO(my_ering,P,i1)
          CALL FIND_PATCH(P%PREVIOUS,P,NEXT=my_TRUE,ENERGY_PATCH=my_FALSE)
          CALL MOVE_TO(my_ering,P,i2)
          CALL FIND_PATCH(P,P%NEXT,NEXT=my_FALSE,ENERGY_PATCH=my_FALSE)
       case('ROTATEPARTOFLAYOUT')

          READ(MF,*)DT
          READ(MF,*) I1,I2
          CALL MOVE_TO(my_ering,P,i1)
          call ROTATE_LAYOUT(my_ering,P%mag%p%f%ent,DT,I1,I2)
          CALL MOVE_TO(my_ering,P,i1)
          CALL FIND_PATCH(P%PREVIOUS,P,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
          CALL MOVE_TO(my_ering,P,i2)
          CALL FIND_PATCH(P,P%NEXT,NEXT=MY_FALSE,ENERGY_PATCH=MY_FALSE)

       case('TRANSLATEFIBREANDPATCH')
          READ(MF,*)POS
          READ(MF,*)DT
          CALL MOVE_TO(my_ering,P,POS)
          CALL TRANSLATE_Fibre(P,DT,ORDER=1,BASIS=P%MAG%P%F%MID)
          CALL FIND_PATCH(P%PREVIOUS,P,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
          CALL FIND_PATCH(P,P%NEXT,NEXT=MY_FALSE,ENERGY_PATCH=MY_FALSE)
       case('DEBUG')
          read(mf,*) debug_flag,debug_acos
       case('VALISHEVON')
          valishev=my_true
          Write(6,*)"Valishev's multipoles are on"
       case('VALISHEVOFF')
          valishev=my_true
          Write(6,*)"Valishev's multipoles are off"
       case('POWERVALISHEV')
          valishev=my_true
          READ(MF,*)POS
          READ(MF,*)A1,B1
          CALL MOVE_TO(my_ering,P,POS)
          p%mag%VA=A1
          p%mag%VS=B1
          p%magp%VA=A1
          p%magp%VS=B1
       case('POWERVALISHEVQUADRUPOLE')
          valishev=my_true
          READ(MF,*)POS
          READ(MF,*)A1,B1
          CALL MOVE_TO(my_ering,P,POS)
          a1=p%mag%bn(2)
          cns=0.0_dp
          CALL ADD(P,2,0,cns)
          a1=0.5_dp*a1*b1**2
          p%mag%VA=A1
          p%mag%VS=B1
          p%magp%VA=A1
          p%magp%VS=B1
       case('VALISHEVALLQUADRUPOLE')
          valishev=my_true
          READ(MF,*)B1
          i2=0
          p=>my_ering%start
          do i1=1,my_ering%n
             if(.not.associated(p%mag%volt).and.associated(p%mag%bn)) then
                if(p%mag%p%nmul>=2) then
                   a1=p%mag%bn(2)
                   cns=0.0_dp
                   CALL ADD(P,2,0,cns)
                   a1=0.5_dp*a1*b1**2
                   p%mag%VA=A1
                   p%mag%VS=B1
                   p%magp%VA=A1
                   p%magp%VS=B1
                   i2=i2+1
                endif
             endif
             p=>p%next
          enddo
          write(6,*) i2, "Quadrupoles globally replaced by Valishev quadrupoles"
       case('UNDOVALISHEVQUADRUPOLE')
          valishev=my_false
          p=>my_ering%start
          i2=0
          do i1=1,my_ering%n
             if(.not.associated(p%mag%volt).and.associated(p%mag%bn)) then
                if(p%mag%p%nmul>=2) then
                   a1= p%magp%VA
                   B1=p%magp%VS
                   a1=a1*2.0_dp/b1**2
                   p%mag%VA=0.0_dp
                   p%mag%VS=0.0_dp
                   p%magp%VA=0.0_dp
                   p%magp%VS=0.0_dp
                   cns=a1
                   i2=i2+1
                   CALL ADD(P,2,1,cns)
                endif
             endif
             p=>p%next
          enddo
          write(6,*) i2,"Valishev quadrupoles globally replaced by normal quadrupoles"
       case('POWERVALISHEVNAME')
          valishev=my_true
          READ(MF,*)name   !,POS
          READ(MF,*)A1,B1
          !          call move_to(my_ering,p,name,POS)
          i2=0
          p=>my_ering%start
          do i1=1,my_ering%n
             if(.not.associated(p%mag%volt).and.associated(p%mag%bn).and.name==p%mag%name) then
                if(p%mag%p%nmul>=2) then
                   i2=i2+1
                   p%mag%VA=A1
                   p%mag%VS=B1
                   p%magp%VA=A1
                   p%magp%VS=B1
                endif
             endif
             p=>p%next
          enddo
          write(6,*) " found ",i2, " Valishev quadrupoles "
       case('POWERVALISHEVNAMEQUADRUPOLE')
          valishev=my_true
          READ(MF,*)name   !,POS
          READ(MF,*)A1,B1
          !          call move_to(my_ering,p,name,POS)
          p=>my_ering%start
          do i1=1,my_ering%n
             if(.not.associated(p%mag%volt).and.associated(p%mag%bn).and.name==p%mag%name) then
                if(p%mag%p%nmul>=2) then
                   i2=i2+1
                   a1=p%mag%bn(2)
                   cns=0.0_dp
                   CALL ADD(P,2,0,cns)
                   a1=0.5_dp*a1*b1**2
                   p%mag%VA=A1
                   p%mag%VS=B1
                   p%magp%VA=A1
                   p%magp%VS=B1
                endif
             endif
             p=>p%next
          enddo
          write(6,*) " found ",i2, " Valishev quadrupoles "
       case('POWERMULTIPOLE')
          READ(MF,*)POS
          READ(MF,*)n,cns, bend_like
          CALL MOVE_TO(my_ering,P,POS)
          CALL ADD(P,N,0,CNS)
          p%mag%p%bend_fringe=bend_like
          p%magp%p%bend_fringe=bend_like
       case('POWERMULTIPOLENAME')
          READ(MF,*)name   !,POS
          READ(MF,*)n,cns, bend_like
          call move_to(my_ering,p,name,POS)
          if(pos/=0) then
             CALL ADD(P,N,0,CNS)
             p%mag%p%bend_fringe=bend_like
             p%magp%p%bend_fringe=bend_like
          else
             write(6,*) name," Not found "
             stop 555
          endif
       case('SETALWAYSON')
         READ(MF,*)name,bend_like
          call move_to(my_ering,p,name,POS)
          if(pos/=0) then
            if(p%mag%kind==kind4) then
              p%mag%c4%always_on=bend_like
              p%magp%c4%always_on=bend_like
            else
             write(6,*) name," is not a cavity "
             stop 555
            endif
          else
             write(6,*) name," Not found "
             stop 555
          endif
         
           

       case('POWERHELICALDIPOLE','POWERHELICAL')

          !SRM DEBUG ... typo
          !B_TESLA --> B_TESLA

          READ(MF,*)name, VORNAME   !,POS
          READ(MF,*) A1,B1    !  A1,B1
          READ(MF,*) NLAM     ! NUMBER OF WAVES
          READ(MF,*) HPHA     ! PHASE
          READ(MF,*) CUR2,CUR1,B_TESLA   ! ACTUAL CURRENT, REF CURRENT, b FIELD
          READ(MF,*) I1,I2    ! NST,METHOD
          call move_to(my_ering,p,name,VORNAME,POS)
          if(pos/=0) then
             W=P
             B_TESLA=CUR2/CUR1*B_TESLA/w%brho

             p%mag%bn(1)=B_TESLA+B1
             p%magp%bn(1)=B_TESLA+B1
             p%mag%an(1)=-B_TESLA+A1
             p%magp%an(1)=-B_TESLA+A1

             !SRM DEBUG ... phase
             p%mag%phas=twopi*HPHA
             p%magp%phas=p%mag%phas

             !SRM DEBUG ... multiply NLAM
             !             p%mag%freq=twopi/p%mag%l/NLAM
             p%mag%freq=twopi/p%mag%l*NLAM
             p%magp%freq=p%mag%freq

             p%mag%p%nst=I1
             p%mag%p%method=I2
             p%magp%p%nst=I1
             p%magp%p%method=I2
          else
             write(6,*) name," Not found "
             stop 555
          endif
       case('COMPUTEMAP')
          READ(MF,*)POS,DEL,NO
          READ(MF,*) FILENAME

          CALL compute_map_general(my_ering,my_estate,filename,pos,del,no)

       case('ZEROSEXTUPOLES')
          call zero_sex(my_ering)
       case('POWERCAVITY')
          c_%CAVITY_TOTALPATH=0 ! fake pill box
          !          c_%phase0=0.0_dp ! because madx is crap anyway
          READ(MF,*)HARMONIC_NUMBER,VOLT,PHASE,c_%phase0,c_%CAVITY_TOTALPATH,epsf
          c_%phase0=c_%phase0*pi
          CALL power_cavity(my_ering,HARMONIC_NUMBER,VOLT,PHASE,epsf)


        case('CAVITYTOTALPATH')
        read(mf,*) pos

        if(pos/=0.and.pos/=1) then
            Write(6,*) "Cavity totalpath must 0 or 1"
            stop 665 
        endif

        call totalpath_cavity(my_ering,pos)

 
       case('NEWEQUILIBRIUMSIZES','EQUILIBRIUMSIZES')
          READ(MF,*) POS,FILENAME,NAME
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ering,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif
          call radia_new(my_ering,POS,my_estate,FILENAME)
       case('SPECIALALEX')
          READ(MF,*) I1
          READ(MF,*) targ_tune(1:2),sc

          CALL special_alex_main_ring(my_ering,i1,targ_tune,sc)

       case('SPECIALALEXAUTO')
          READ(MF,*) I1,epsf
          READ(MF,*) targ_tune(1:2),sc

          CALL special_alex_main_ring_auto(my_ering,i1,targ_tune,sc,epsf)


       case default

          write(6,*) " Command Ignored "

       end select


    enddo
100 continue
!    if(associated(my_old_state)) my_estate=>my_old_state

    write(6,*) " Exiting Command File ", ptc_fichier(1:len_trim(ptc_fichier))

     longprint=longprintt 

    if(mf/=5) close(mf)

  END subroutine read_ptc_command

!!!!!!!!!!!!!!!!   polarization scan  !!!!!!!!!!!!!!!
subroutine  compute_polarisation(r,no,ang,nlm,lm,nturns,ngen,snake,isnake,plotfile,datafile)
 implicit none
type(layout), pointer :: r
character(*) plotfile,datafile
CHARACTER(*) snake
integer isnake,nlm
!!!!!!  PTC stuff
real(dp) closed_orbit(6) , x(6),cut,energy,deltap,de,rotator,n_ave(3),xij,ang(3),dlm
real(dp)  n0i(3),spin_damp(6,6),e_ij(6,6) ,circum,lm(2) 

type(internal_state), target :: state,state_trackptc ,state0 
type(probe) ray   
type(probe_8) rayp
type(c_damap) M,ID,U_1,as,a0,a1,a2   
 integer i,k,no,mf,mfd,j,l,kprint,count,rplen,mfisf,kp,jj(6),i1,i2,i3,i4 !,n
type(fibre), pointer ::f  
type(integration_node), pointer :: it  
character(255) mapfile

type(tree_element_zhe)     t_olek_map(1:3) 
type(probe_zhe) xs0_zhe 
integer nbunch,je(6),ngen
type(bunch) bunch_zhe  
 
type(c_taylor) phase(4)
real(dp) :: nturns, xj,deb ,dpol 
 type(c_normal_form) nf
logical ip,doit,track_ptc,reload ,old_use_quaternion
real(dp) tune(1:3) , spin_tune 

TYPE(c_spinor) ISF  
     type(quaternion) q,q0
 type(q_linear)  q_c,q_ptc

    call get_length(r,circum)
    circum=circum/clight
 
 call kanalnummer(mfd,datafile) 
call kanalnummer(mf,plotfile)

! enforces quaternion rather than SO(3)
old_use_quaternion=use_quaternion
use_quaternion=.true.
 !!! removing some annoying prints 
c_verbose=.false.
lielib_print(4)=0
track_ptc=.false.

 reload=.true.
 
  if(track_ptc) then
kprint=nturns/1000
else
kprint=nturns/10
endif
nbunch=ngen**3

  call zhe_ini(use_quaternion)

write(mfd,*) nlm
dlm=(lm(2)-lm(1))/nlm

do i4=0,nlm-1
call alloc_bunch(bunch_zhe,nbunch)


 rotator=(lm(1)+i4*(dlm ))*twopi
 

 
!! initializes outside tracking

 

! PTC tracking state
state0=time0+spin0 +radiation0
 state=state0+envelope0    
state_trackptc=state0 +stochastic0


! the state used in my windows graphical interface. Not harmful.
 file_zhe="olek"

 
!!!  r is a lattice called layout in PTC. Linked list of so-called fibres: magnets for your purpose
r=>m_u%start
!! creates an integration node layout, necessary for spin and looking inside magnets
if(associated(r%t)) call make_node_layout(r)

 
!  calculation of beam size using envelope theory
call radia_new(r,1,state0,"radia.dat",e_ij=e_ij,spin_damp=spin_damp ,ngen=ngen,bunch_zhe=bunch_zhe  )


 
! "it" points to an integration node
it=>r%start%t1
 
cut=1.d0/sqrt(ang(1)**2+ang(2)**2+ang(3)**2)
k=0
f=>r%start
do i=1,r%n
if(f%mag%name(1:isnake)==snake(1:isnake)) then

f%patch%patch=5

f%patch%a_ang=ang*cut*rotator
write(6,'(a,a,a)') snake(1:isnake), " found at ",f%mag%name(1:isnake)
exit
endif
f=>f%next
enddo

number_zhe_maps = 1
!!! alloc an array of pointers to my fibres 


 it=>r%start%t1
 


!!!! finds the closed orbit at position s=0
closed_orbit=0
call FIND_ORBIT_x(closed_orbit,state0,1.d-8,node1=it)
 
    call GET_loss(r,energy,deltap)
 
    write(6,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap

IF(.NOT.check_stable) THEN
 WRITE(6,*) " UNSTABLE  IN ORBIT SEARCHER "
 stop
ENDIF
write(6,'(6(E20.13,1X))' ) closed_orbit
 
!!! First order maps for tests.
 
 ray=closed_orbit
call FIND_ORBIT_x(ray%x,state0,1.d-5,node1=it)


 
 
 

call alloc(M,ID,U_1,as,a0,a1,a2)
call alloc(rayp)
call alloc(nf)
call alloc(isf)
 
ID=1
rayp=ray + ID   ! closed orbit added to identity and thrown into a polymorphic ray. 
!See definition of C_damap and probe_8 in h_definition.f90
 
 
 call  propagate(rayp,state0,node1=it)
id=rayp
  
 
 !call kill(L_ns , N_pure_ns, L_s , N_s)

 call c_normal(id,nf,dospin=.true.)
 


tune=nf%tune(1:3)
spin_tune=nf%spin_tune

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
write(6,*) "tunes "
write(6,*) nf%tune(1:3)
write(6,*) "spin tune "
write(6,*) nf%spin_tune
 
call c_full_canonise(nf%atot,U_1,as,a0,a1,a2)

ISF=2   ! (Fa)
call makeso3(as)
ISF=as%s*ISF ! (Fb)
 

 
ray=closed_orbit
q0%x(0)=0
q0%x(1)=ISF%v(1)
q0%x(2)=ISF%v(2)
q0%x(3)=ISF%v(3)
n0i=q0%x(1:3)

cut=0
do i1=1,3
do i2=1,6
do i3=1,6
 
cut =  cut + (ISF%v(i1).d.i2)*e_ij(i2,i3)*(ISF%v(i1).d.i3)

enddo
enddo
enddo
 
 cut=-cut/2
 
write(6,*) " Polarization computed with dn/dz + final polarization after ",nturns," turns"
write(6,*)  cut,exp(nturns*cut)
deb=cut
dpol=exp(nturns*deb)

ray%q=q0

 

 call propagate(ray,state0,node1=it)

 

write(6,*) " Q0 and q(Q0)q^-1 "
 q=ray%q*q0*ray%q**(-1)
 
call kill(M,ID,U_1,as,a0,a1,a2)
call kill(rayp)
call kill(nf)
call kill(isf)
 
call init(state0,no,0)

call alloc(M,ID,U_1,as,a0,a1,a2)
call alloc(rayp)
call alloc(nf)
call alloc(isf)
call alloc(phase)
!!!!!!  generate distribution  !!!!!
 
  id=1
 
ray=closed_orbit
 
rayp=ray+id
   x=rayp%x
  call propagate(rayp,state,node1=it)
 
!!! create files "olek#" where #=1,2,3,....,n
!!! maps are  put in there
mapfile="olek"
 !!! probe_8 thrown into a C_dmap
  M=rayp
 !!! closed orbit saved
 ray=rayp
!!! initial map local closed orbit + identity  re-created for step i+1
 rayp=ray+ID
 !!! tracking object created for step i
!call print(m)
  
m%x0(1:6)=x
call print(m)

call fill_tree_element_line_zhe_outside_map(m ,mapfile,as_is=.false.,stochprec=1.d-8) 
 
 
call read_tree_zhe(t_olek_map(1:3),mapfile)


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
 
 !bunch_zhe=closed_orbit
 !!!!!!!!!!!!!!1
  
 do i=1,bunch_zhe%n
 
    bunch_zhe%xs(i)%q=1.d0
  
 enddo

xs0_zhe=closed_orbit  
write(6,*)  " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "
 xj=0
j=0

 
 
do while(.true.)
xj=xj+1
j=j+1
if(xj>nturns) exit

do k=1,nbunch
 
 if( bunch_zhe%stable(k)) then
  
!!! one-turn map
    call track_TREE_probe_complex_zhe(t_olek_map(1:3),bunch_zhe%xs(k),spin=.true.,rad=.true.,stoch=.true.)  !stoch=state%stochastic)
  if(bunch_zhe%xs(k)%u) exit

     if(bunch_zhe%xs(k)%u) then
     bunch_zhe%stable(k)=.false.
     bunch_zhe%turn(k)=xj
     bunch_zhe%r=bunch_zhe%r-1
     write(6,*) K, " lost at turn ",xj
     if(reload) then
       rplen=0
       do while(.not. (rplen>=1.and.rplen<=bunch_zhe%n.and.rplen/=k))  
         cut= RANF() *bunch_zhe%n
         rplen=nint(cut)
       enddo 
       bunch_zhe%xs(k)%u=.false.
       bunch_zhe%r=bunch_zhe%r+1
       bunch_zhe%turn(k)=0
        bunch_zhe%stable(k)=.true.
      bunch_zhe%xs(k)=bunch_zhe%xs(rplen)
      bunch_zhe%reloaded=bunch_zhe%reloaded+1
      write(6,*) K, " reloaded by ",rplen

    

    endif

 endif
endif



!write(mf,'(6(E20.13,1X))' ) bunch_zhe(k)%x-closed_orbit
enddo
 
 call get_polarisation(bunch_zhe,q0,n_ave,de)
 

if(mod(j,kprint)==0) then   !.or.xj>nturns-10) then
 write(mf,'(6(E20.13,1X))' ) xj,exp(xj*deb),de,n_ave
 ! write(6,*) xj
!  if(track_ptc)
 write(6,'(6(E20.13,1X))' ) xj,exp(xj*deb),de,n_ave

j=0

endif
enddo


 
write(mf,*)  " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "
 write( mf,*) "   survived        ", bunch_zhe%r," out of ",    bunch_zhe%n  
 write( mf,*) "   reloaded        ", bunch_zhe%reloaded
 write( mf,*) "   track_ptc        ", track_ptc     
    
 write( mf,*) " nturns             ", nturns          
 write( mf,*) " nbunch             ", nbunch          
 write( mf,*) " rotator/2pi  ", rotator/twopi       
 write( mf,*) " no                 ", no 
 write( mf,'((a20),(a))') " plot file          ", plotfile(1:len_trim(plotfile)) 
 write(mf,*) " q0 "
call print(q0,mf)
 write(mf,*)
 write(mf,*)
 
 write(mf,*) "computed depolarization = ",dpol
 write(mf,*) "tracked  depolarization = ",de
write(mf,*)  " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "

write(6,*)  " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "
write( 6,*) "   survived        ", bunch_zhe%r," out of ",    bunch_zhe%n  
write( 6,*) "   reloaded        ", bunch_zhe%reloaded
 write( 6,*) "   track_ptc        ", track_ptc     
 
 write( 6,*) " nturns             ", nturns          
 write( 6,*) " nbunch             ", nbunch          
 write( 6,*) " rotator/2pi  ", rotator/twopi       
 write( 6,*) " no                 ", no       
 write( 6,'((a20),(a))') " plot file          ", plotfile(1:len_trim(plotfile)) 
write(6,*) " q0 "
call print(q0)
 write(6,*)
 write(6,*)
 

 write(6,*) "computed depolarization = ",dpol
 write(6,*) "tracked  depolarization = ",de
write(6,*)  " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "

 write(mfd,'(11(E20.13,1X))' )tune(1:3),spin_tune,xj,circum*xj,exp(xj*deb),de,n_ave


call kill(M,ID,U_1,as,a0,a1,a2)
call kill(rayp)
call kill(nf)
call kill(isf)
call kill(phase)

call kill_bunch(bunch_zhe)

enddo
use_quaternion=old_use_quaternion

close(mf)

close(mfd)


end subroutine compute_polarisation

subroutine get_polarisation(b,q0,n_ave,de)
implicit none
 type(quaternion) q0,q
 real(dp) n_ave(3) ,de
 integer n,i,l
 type(bunch)  b
 
 
 n=b%n
  n_ave=0

 do i=1,n
  if(b%stable(i)) then
 do l=0,3
  q%x(l)= b%xs(i)%q%x(l)
 enddo
   q=q*q0*q**(-1)
  
  n_ave=n_ave + q%x(1:3)
 endif
 enddo

 n_ave=n_ave/b%r
de=sqrt(n_ave(1)**2+n_ave(2)**2+n_ave(3)**2)
 
end subroutine get_polarisation

   real(dp) function dis_gaussian(r)
    implicit none
    real(dp) r1,r2,x ,r


    R1 = -LOG(1.0_dp-r)
    R2 = 2.0_dp*PI*RANF()
    R1 = SQRT(2.0_dp*R1)
    X  = R1*COS(R2)
     dis_gaussian=x
 
    RETURN
  END function dis_gaussian

 SUBROUTINE radia_new(R,loc,estate,FILE1,fix,em,sij,sijr,tune,damping,e_ij,spin_damp,init_tpsa,ngen,bunch_zhe,file_bunch)
    implicit none
    TYPE(LAYOUT) R

    REAL(DP) X(6),m,energy,deltap
    CHARACTER(*), optional :: FILE1,file_bunch
    integer, optional :: ngen
    type(c_damap)  Id,a0 
    type(c_normal_form) normal
    integer  i,j ,i1,i2,i3,i4
    real(dp), optional :: fix(6), em(3),sij(6,6),tune(3),damping(3),e_ij(6,6),spin_damp(6,6)
    type(bunch), optional :: bunch_zhe
    logical, optional ::  init_tpsa
    complex(dp), optional :: sijr(6,6)   
    TYPE(INTERNAL_STATE) state
    TYPE(INTERNAL_STATE), target :: estate
    integer loc,mf1,mfg
    type(fibre), pointer :: p
    type(probe) xs0
    type(probe_8) xs
    character*48 fmd,fmd1
    real(dp) mat(6,6),matf(6,6),ki(6),ray(6),a1(3)
 

    if(present(FILE1)) then
    call kanalnummer(mf1)
    open(mf1,file=FILE1)
    endif
fmd= '(a12,1X,a3,I1,a3,i1,a4,D18.11,1x,D18.11)'
fmd1='(1X,a3,I1,a3,i1,a4,2(D18.11,1x),(f10.3,1x),a2)'


    if(present(init_tpsa)) then
     if(init_tpsa)     then 
        state=(estate-nocavity0)+radiation0
     else
     state =estate
    endif
    else 
        state=(estate-nocavity0)+radiation0
    endif
    x=0.d0

!do i=1,10
!radfac=float(i)/10.d0
 
    CALL FIND_ORBIT_x(R,X,STATE,1.0e-8_dp,fibre1=loc)
!write(6,format6) x
    if(.not.check_stable) then
      write(6,*) "Unstable in radia_new ",i
      stop
     endif
!enddo

    if(present(FILE1)) then
        WRITE(mf1,*) " CLOSED ORBIT AT LOCATION ",loc
        write(mf1,*) x
    endif
     if(present(fix)) fix=x


    call GET_loss(r,energy,deltap)
       if(present(FILE1)) then
    write(mf1,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap
    write(mf1,*) " stable closed orbit tracked "
    write(mf1,"(6(1X,D18.11))") x
    write(mf1,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap
    endif

   if(present(init_tpsa)) then
     if(init_tpsa)     then 
        CALL INIT(state,1,0)
    endif
    else 
        CALL INIT(state,1,0)
    endif
 
    CALL ALLOC(NORMAL)
    CALL ALLOC(ID,a0)
    call alloc(xs)

!    if(i1==1) normal%stochastic=my_true
    xs0=x
    ID=1
    xs=XS0+ID

    state=state+envelope0
 
    CALL TRACK_PROBE(r,xs,state, fibre1=loc)
     id=xs
 
    if(present(e_ij)) e_ij=xs%e_ij
   if(present(FILE1)) then
    write(mf1,*) " Full Map "    
     call print(id,mf1)
    write(mf1,*) " End of Full Map "
   endif
  call  c_normal(id,normal)             ! (8)
if(present(ngen)) then
matf=0
   if(present(bunch_zhe)) then
   if(ngen**3>bunch_zhe%n) then
    Write(6,*) " problems in radia "
    stop 388
   endif
   endif
   if(present(file_bunch)) call kanalnummer(mfg,file_bunch)
  a0=1
  a0%e_ij=normal%s_ij0 
  call c_stochastic_kick(a0,mat,ki,1.d-38)
  i4=0
  do i1=0,ngen-1
  do i2=0,ngen-1
  do i3=0,ngen-1
    i4=i4+1
    ray=0
     a1(1)=float(i1)/float(ngen)
     a1(2)=float(i2)/float(ngen)
     a1(3)=float(i3)/float(ngen)
    ray(1)= ki(1)*dis_gaussian(a1(1))
    ray(2)= ki(2)*dis_gaussian(a1(1))
    ray(3)= ki(3)*dis_gaussian(a1(2))
    ray(4)= ki(4)*dis_gaussian(a1(2))
    ray(5)= ki(5)*dis_gaussian(a1(3))
    ray(6)= ki(6)*dis_gaussian(a1(3))
    ray=matmul(mat,ray)
    if(present(file1)) then
    do i=1,6
    do j=i,6
    matf(i,j)=matf(i,j) + ray(i)*ray(j)
    enddo
    enddo
    endif
     ray=ray+x
     If(present(bunch_zhe)) bunch_zhe%xs(i4)%x=ray
     if(present(file_bunch))   write(mfg,format6) ray
  enddo
  enddo
  enddo
    if(present(file1)) then
     write(mf1,*) ' Generated  and calculated distributions using "ngen" '
    matf=matf/ngen**3
    do i=1,6
    do j=i,6
     write(mf1,*) i,j,matf(i,j),real(normal%s_ij0(i,j))
    enddo
    enddo
    endif
   if(present(file_bunch))   close(mfg)
endif
   if(present(FILE1)) then
    write(mf1,*)" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
    write(mf1,*)" Tunes "
    write(mf1,*) normal%tune
    write(mf1,*)" Damping "
    write(mf1,*) normal%damping
    write(mf1,*)" Emittances "
    write(mf1,*) normal%emittance
    write(mf1,*)" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
    write(mf1,*)"   "

    write(mf1,*)" Equilibrium Beam Sizes "
    do i=1,6
       do j=i,6
          write(mf1,*) i,j,normal%s_ij0(i,j)
       enddo
    enddo
write(mf1,*)
write(mf1,'(16X,a50)') "   Equilibrium moments in Phasors Basis           "
 do i=1,6
 do j=i,6 
  if(abs(normal%s_ijr(i,j))>1.d-20) then
   write(mf1,fmd) " Phasors -> ","<x_",i," x_",j,"> = ",  &   ! (14)
                    c_clean(normal%s_ijr(i,j),1.d-20)  
  endif
 enddo
 enddo 
    close(mf1)
endif
if(present(spin_damp)) then
mat=id
matf=mat

call furman_symp(matf)


id=0
id=matf
 
id=id**(-1)
 
matf=id 
 
matf=matmul(matf,mat)
do i=1,6
 matf(i,i)=matf(i,i)-1
enddo
 
 spin_damp=matmul(matmul(matf,normal%s_ij0),transpose(matf))
 
endif

if(present(em)) em=normal%emittance
if(present(tune)) tune=normal%tune(1:3)
if(present(damping)) damping=normal%damping(1:3)
if(present(sij)) then
do i=1,6
do j=1,6
 sij(i,j)=normal%s_ij0(i,j)
enddo
enddo
endif
if(present(sijr)) then
do i=1,6
do j=1,6
 sijr(i,j)=normal%s_ijr(i,j)
enddo
enddo
endif




    CALL KILL(NORMAL)
    CALL KILL(ID,a0)
    CALL KILL(xs)

  end subroutine radia_new

  subroutine totalpath_cavity(r,j)
    implicit none
    TYPE(LAYOUT), POINTER :: r
    type(fibre), pointer :: p
    integer i,j
 

    p=>r%start
    do i=1,r%n

       if(p%mag%kind==kind4) then
          write(6,*) " cavity found ",p%mag%name,p%mag%vorname
          p%mag%c4%CAVITY_TOTALPATH=mod(j,2)
          p%magp%c4%CAVITY_TOTALPATH=mod(j,2)
       endif

       p=>P%NEXT

    enddo
    c_%CAVITY_TOTALPATH=j
  end subroutine totalpath_cavity

  subroutine power_cavity(r,HARM,VOLT,PHAS,prec)
    implicit none
    TYPE(LAYOUT), POINTER :: r
    type(fibre), pointer :: p
    integer i,ip,HARM,j
    type(internal_state) state
    real(dp) closed(6),s,VOLT,PHAS,accuracy,circum,freq
    real(dp),optional :: prec

    call get_length(r,circum)

    write(6,*) " fiducial length ",circum

    state=my_estate+nocavity0-totalpath0

    accuracy=1.0e-10_dp
    if(present(prec)) accuracy=prec
    closed=0.0_dp
    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)

    closed(6)=0.d0
    call track(r,closed,1,state+totalpath0+time0)
    write(6,*) " CT = ",closed(6)

    freq=HARM*CLIGHT/closed(6)

    ip=0
    p=>r%start
    do i=1,r%n
       closed(6)=0.0_dp
       call track(r,closed,i,i+1,state)
       s=closed(6)
       circum=circum+s
       IF(ABS(s)>accuracy) THEN
          ip=ip+1
          p%patch%b_l=-s
          p%patch%b_t=-s/P%BETA0
          p%patch%time=2
       ELSE
          p%patch%time=0
       ENDIF
       p=>P%next
    enddo

    if(ip/=0) then
       WRITE(6,*)"?????????????????????????????????????????????????????????????"
       WRITE(6,*)"?????????????????????????????????????????????????????????????"
       WRITE(6,*)IP,"TIME PATCHES WERE NEEDED BECAUSE BENDS HAVE BEEN ADJUSTED "
       WRITE(6,*)"THEY WERE PUT AT THE END OF SEVERAL FIBRES "
       WRITE(6,*)"TIME PATCH MAYBE NEEDED IN IDEAL LATTICES IF BENDS ARE ADJUSTED "
       WRITE(6,*) " ACTUAL LENGTH ",CIRCUM
       WRITE(6,*)"?????????????????????????????????????????????????????????????"
       WRITE(6,*)"?????????????????????????????????????????????????????????????"
    endif


    p=>r%start
    do i=1,r%n

       if(p%mag%kind==kind4) then
          write(6,*) " Before "

          ! write(6,*) p%mag%name
          ! write(6,*) " volt    = ",p%mag%volt
          ! write(6,*) " freq    = ",p%mag%freq
          ! write(6,*) " phas    = ",p%mag%phas
          ! write(6,*) " ref p0c = ",p%mag%p%p0c
          ! write(6,*) "electron = ", c_%electron
          if(p%mag%l/=0.0_dp) then
             p%mag%volt=VOLT/p%mag%l
             p%magp%volt=p%mag%volt
          else
             p%mag%volt=VOLT   !/p%mag%l
             p%magp%volt=p%mag%volt
          endif
          p%mag%freq=freq      !CLIGHT*HARM/circum   ! harmonic=120.0_dp
          p%magp%freq=p%mag%freq
          p%mag%phas=PHAS
          p%magp%phas=p%mag%phas
          write(6,*) " After "

          write(6,*) p%mag%name
          write(6,*) " volt    = ",p%mag%volt
          write(6,*) " freq    = ",p%mag%freq
          write(6,*) " phas    = ",p%mag%phas
          write(6,*) " ref p0c = ",p%mag%p%p0c
          p%mag%c4%phase0=0.0_dp
          p%magp%c4%phase0=0.0_dp
          p%mag%c4%f(1)=1.0_dp
          p%magp%c4%f(1)=1.0_dp
          p%mag%c4%ph(1)=0.0_dp
          p%magp%c4%ph(1)=0.0_dp
          p%mag%c4%CAVITY_TOTALPATH=0
          p%magp%c4%CAVITY_TOTALPATH=0
          do j=2,p%mag%c4%nf
             p%mag%c4%f(j)=0.0_dp
             p%magp%c4%f(j)=0.0_dp
             p%mag%c4%ph(j)=0.0_dp
             p%magp%c4%ph(j)=0.0_dp
          enddo
          !          write(6,*) "electron = ", c_%electron
          !          write(6,*) " phase0  = ", c_%phase0
          if(c_%CAVITY_TOTALPATH==0) write(6,*) " fake cavity "
       endif
       p=>P%NEXT

    enddo

  end subroutine power_cavity

  subroutine zero_sex(r)
    implicit none
    TYPE(LAYOUT), POINTER :: r
    type(fibre), pointer :: p
    integer i

    p=>r%start
    do i=1,r%n
       if(associated(p%mag%bn)) then
          if(p%mag%p%nmul>=3) then
             call add(p,3,0,0.0_dp)
          endif
       endif
       p=>p%next
    enddo

  end subroutine zero_sex

  subroutine charge_dir(r)
    implicit none
    TYPE(LAYOUT), POINTER :: r
    type(fibre), pointer :: p
    integer i
    p=>r%start
    do i=1,r%n
       p%dir=-p%dir
       p=>p%next
    enddo

  end subroutine charge_dir


  SUBROUTINE remove_drifts(R,NR)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R,NR
    integer I
    type(fibre), pointer :: P
    logical(lp) doneit

    p=>r%start

    do i=1,r%n
       IF(P%MAG%KIND/=KIND0.AND.P%MAG%KIND/=KIND1) THEN

          !        call APPEND_EMPTY(NR)
          CALL APPEND( NR, P )

       ENDIF
       P=>P%NEXT
    ENDDO

    NR%closed=.true.
    doneit=.true.
    call ring_l(NR,doneit)


    write(6,*) " do you want patching ?"
    read(5,*) i
    if(i==0) return

    p=>nr%start

    do i=1,nr%n-1
       CALL FIND_PATCH(P,P%next,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)

       P=>P%NEXT
    ENDDO
    CALL FIND_PATCH(P,P%next,NEXT=my_false,ENERGY_PATCH=MY_FALSE)

    ! avoiding putting a patch on the very first fibre since survey does not allow it....


  end SUBROUTINE remove_drifts

  SUBROUTINE remove_drifts_bends(R,NR)  ! special example to be removed later
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R,NR
    integer I,IG
    type(fibre), pointer :: P,bend
    logical(lp) doneit
    real(dp) ent(3,3),a(3),ang(3),d(3)

    p=>r%start
    bend=>r%next%start%next   ! second layout in universe
    write(6,*) " using bend called ",bend%mag%name
    write(6,*) " 'USING SURVEY' TYPE 1 / 'USING GEOMETRY' TYPE 0 "
    READ(5,*) IG
    do i=1,r%n
       IF(P%MAG%KIND/=KIND0.AND.P%MAG%KIND/=KIND1.and.P%MAG%p%b0==0.0_dp) THEN

          CALL APPEND( NR, P )
       elseif(P%MAG%p%b0/=0.0_dp) then
          bend%mag%p%bend_fringe=.true.
          bend%magp%p%bend_fringe=.true.
          bend%mag%L=P%MAG%p%lc
          bend%magp%L=P%MAG%p%lc   ! give it correct arc length
          bend%mag%p%Lc=P%MAG%p%lc
          bend%magp%p%Lc=P%MAG%p%lc   ! give it correct arc length
          bend%mag%p%Ld=P%MAG%p%lc
          bend%magp%p%Ld=P%MAG%p%lc   ! give it correct arc length
          call add(bend,1,0,p%mag%bn(1))    ! Give a huge B field to quadrupole, i.e. looks like a kicker now
          CALL APPEND( NR, bend )
          ent=p%chart%f%mid     !  storing the bend location
          a=p%chart%f%a         !
          !     since we use a quadrupole, the entrance frame of this quad is the mid frame of the bend
          !     The  fibre bend must be rotated and translated into position
          ! easiest way is to survey it with initial condition correspounding to the actual position and orientation
          !
          IF(IG==1) THEN
             call SURVEY(nr%end,ENT,A)
          ELSE
             d=a-nr%end%chart%f%a
             CALL TRANSLATE_Fibre(nr%end,D)  ! translation in global frame
             CALL COMPUTE_ENTRANCE_ANGLE(nr%end%chart%f%ENT,ENT,ANG)
             CALL ROTATE_Fibre(nr%end,A,ang)  ! translation in global frame
          ENDIF


       ENDIF
       P=>P%NEXT
    ENDDO

    NR%closed=.true.
    doneit=.true.
    call ring_l(NR,doneit)



    p=>nr%start

    do i=1,nr%n-1
       CALL FIND_PATCH(P,P%next,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)

       P=>P%NEXT

    ENDDO
    CALL FIND_PATCH(P,P%next,NEXT=my_false,ENERGY_PATCH=MY_FALSE)

    ! avoiding putting a patch on the very first fibre since survey is not a self-check in that case


  end SUBROUTINE remove_drifts_bends




  SUBROUTINE print_frames(R,filename)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer I,mf
    type(fibre), pointer :: P
    character(*) filename
    call kanalnummer(mf)


    open(unit=mf,file=filename)
    write(mf,*) "Contains location of each fibre and the magnet within the fibre "
    write(mf,*) "N.B. Drifts and Markers are fibres in PTC "
    p=>r%start
    do i=1,r%n


       !   INTEGER(2), POINTER:: PATCH    ! IF TRUE, SPACIAL PATCHES NEEDED
       !   INTEGER, POINTER :: A_X1,A_X2   ! FOR ROTATION OF PI AT ENTRANCE = -1, DEFAULT = 1 ,
       !   INTEGER, POINTER :: B_X1,B_X2   ! FOR ROTATION OF PI AT EXIT = -1    , DEFAULT = 1
       !   REAL(DP),DIMENSION(:), POINTER:: A_D,B_D      !ENTRACE AND EXIT TRANSLATIONS  A_D(3)
       !   REAL(DP),DIMENSION(:), POINTER:: A_ANG,B_ANG   !ENTRACE AND EXIT ROTATIONS    A_ANG(3)
       !   INTEGER(2), POINTER:: ENERGY   ! IF TRUE, ENERGY PATCHES NEEDED
       !   INTEGER(2), POINTER:: TIME     ! IF TRUE, TIME PATCHES NEEDED
       !   REAL(DP), POINTER:: A_T,B_T     ! TIME SHIFT NEEDED SOMETIMES WHEN RELATIVE TIME IS USED
       write(mf,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       write(mf,*) "  "
       write(mf,*) "|| position = ", i,' || PTC kind = ', P%mag%kind," || name = ",P%mag%name, " ||"
       write(mf,*) "  "
       if(p%patch%patch==1.or.p%patch%patch==3) then
          write(mf,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
          write(mf,*) " Entrance geometrical Patch "
          write(mf,*) " Translations A_D(3) "
          write(mf,*) P%patch%a_d
          write(mf,*) " Rotations A_ANG(3) || PI rotations ->   ",p%patch%a_x1,p%patch%a_x2
          write(mf,*) P%patch%A_ANG
          write(mf,*) " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
          write(mf,*) "  "
       endif
       write(mf,*) " Fibre positioning or Ideal position in conventional parlance"
       write(mf,*) "  "
       write(mf,*) " Entrance origin A(3) "
       write(mf,*) P%chart%f%a
       write(mf,*) " Entrance frame (i,j,k) basis in the ent(3,3) array "
       write(mf,*) P%chart%f%ent(1,:)
       write(mf,*) P%chart%f%ent(2,:)
       write(mf,*) P%chart%f%ent(3,:)
       write(mf,*) " Middle origin O(3) "
       write(mf,*) P%chart%f%o
       write(mf,*) " Middle frame (i,j,k) basis in the ent(3,3) array "
       write(mf,*) P%chart%f%mid(1,:)
       write(mf,*) P%chart%f%mid(2,:)
       write(mf,*) P%chart%f%mid(3,:)
       write(mf,*) " Exit origin B(3) "
       write(mf,*) P%chart%f%B
       write(mf,*) " Exit frame (i,j,k) basis in the ent(3,3) array "
       write(mf,*) P%chart%f%exi(1,:)
       write(mf,*) P%chart%f%exi(2,:)
       write(mf,*) P%chart%f%exi(3,:)
       write(mf,*) "  "
       write(mf,*) " Actual magnet positioning  "
       write(mf,*) "  "
       write(mf,*) " Entrance origin A(3) "
       write(mf,*) P%mag%p%f%a
       write(mf,*) " Entrance frame (i,j,k) basis in the ent(3,3) array "
       write(mf,*) P%mag%p%f%ent(1,:)
       write(mf,*) P%mag%p%f%ent(2,:)
       write(mf,*) P%mag%p%f%ent(3,:)
       write(mf,*) " Middle origin O(3) "
       write(mf,*) P%mag%p%f%o
       write(mf,*) " Middle frame (i,j,k) basis in the ent(3,3) array "
       write(mf,*) P%mag%p%f%mid(1,:)
       write(mf,*) P%mag%p%f%mid(2,:)
       write(mf,*) P%mag%p%f%mid(3,:)
       write(mf,*) " Exit origin B(3) "
       write(mf,*) P%mag%p%f%B
       write(mf,*) " Exit frame (i,j,k) basis in the ent(3,3) array "
       write(mf,*) P%mag%p%f%exi(1,:)
       write(mf,*) P%mag%p%f%exi(2,:)
       write(mf,*) P%mag%p%f%exi(3,:)
       if(p%patch%patch==2.or.p%patch%patch==3) then
          write(mf,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
          write(mf,*) " Exit geometrical Patch "
          write(mf,*) " Translations B_D(3) "
          write(mf,*) P%patch%b_d
          write(mf,*) " Rotations B_ANG(3) || PI rotations ->   ",p%patch%b_x1,p%patch%b_x2
          write(mf,*) P%patch%B_ANG
          write(mf,*) " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
          write(mf,*) "  "
       endif
       write(mf,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

       P=>P%NEXT
    ENDDO

    close(mf)
  end SUBROUTINE print_frames

  subroutine printframes(filenameIA)
    use madx_ptc_module
    implicit none
    character*48 charconv
    !   include 'twissa.fi'
    integer   filenameIA(*)
    character(48) filename

    filename = charconv(filenameIA)
    call print_frames(my_ering,filename)

  end subroutine printframes




  subroutine Universe_max_n(n)
    !use build_lattice
    implicit none
    integer n,i
    type(layout), pointer :: L
    n=0

    l=>m_u%start
    do i=1,m_u%n
       n=n+l%n
       l=>l%next
    enddo

  end subroutine Universe_max_n

  subroutine Universe_max_node_n(n)
    !use build_lattice
    implicit none
    integer n,i
    type(layout), pointer :: L
    n=0

    l=>m_u%start
    do i=1,m_u%n
       if(associated(l%t) ) n=n+l%t%n
       l=>l%next
    enddo

  end subroutine Universe_max_node_n

  
 




subroutine read_ptc_rays(filename)
implicit none
character(*) filename
integer mf,nr,i
type(probe) r0
real(dp) X(6)

call kanalnummer(mf,filename)
r0=0
X=0
R0=X

read(mf,*) nr,r0%nac

allocate(xs0g(nr))

do i=1,r0%nac
 read(mf,*) r0%ac(i)%om
 read(mf,*) r0%ac(i)%x
enddo

do i=1,nr
 xs0g(i)=0
 xs0g(i)=r0
 read(mf,*) xs0g(i)%x
enddo
close(mf)

end subroutine read_ptc_rays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Oleksii's   Hermite   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mp=noh
! z(1:2)
! h(2)
! ord1  ord2  (0 no derivatives)
!  U is %he
! val is array of size (0:ord1,0:ord2)

subroutine eval_new(mp, z, h, ord1, ord2,  u, val)
    implicit none
    integer, intent(in) :: mp, ord1, ord2
    real(kind=dp), dimension(2), intent(in) :: z
    real(kind=dp), dimension(2), intent(in) :: h
    real(kind=dp), dimension(0:mp, 0:mp, 0:1, 0:1), intent(in) :: u
    real(kind=dp), intent(out) :: val(0:ord1, 0:ord2)
    integer :: i, j, k, idx, idy
    real(kind=dp), dimension(0:2*mp+1,0:2*mp+1) :: scofsmp, smcofsmp
    real(kind=dp), dimension(0:2*mp+1) :: mpdata
    real(kind=dp) :: zx, zy, xfact, yfact

    i = 0
    j = 0


    do idy = 0,mp
        mpdata(0:mp)        = u(:, idy, i  ,j)
        mpdata(mp+1:2*mp+1) = u(:, idy, i+1,j)
        call interpolate(mp, mpdata, smcofsmp(:,idy))

    end do
    do idy = 0,mp
        mpdata(0:mp)        = u(:,idy,i,   j+1)
        mpdata(mp+1:2*mp+1) = u(:,idy,i+1, j+1)
        call interpolate(mp, mpdata, smcofsmp(:,mp+1+idy))
    end do

    do idx = 0,2*mp+1
        call interpolate(mp, smcofsmp(idx,:), scofsmp(idx,:))
    end do



   do idy = 0, ord2
    do idx = 0, ord1
        val(idx, idy) = 0._dp
        zy = 1._dp
        do j = idy, 2*mp+1
        zx = 1._dp 
        do i = idx, 2*mp+1 ! - j   ! for total degree etienne
            xfact = 1._dp
            yfact = 1._dp
            do k = i-idx+1,i
               xfact = xfact * dble(k) / h(1)
            enddo
            do k = j-idy+1,j
               yfact = yfact * dble(k) / h(2)
            enddo
            val(idx, idy) = val(idx, idy) + xfact*yfact*scofsmp(i,j) * zx * zy
            zx = zx * z(1)
        enddo
        zy = zy * z(2)
        enddo
    enddo
    enddo
end subroutine eval_new






subroutine interpolate(m, u, cofs)
    integer, intent(in) :: m
    real(kind=dp), dimension(0:2*m+1), intent(in) :: u
    real(kind=dp), dimension(0:2*m+1), intent(out) :: cofs
    integer :: i,j
    real(kind=dp), dimension(0:2*m+1, 0:2*m+1) :: dd

    ! Set data
    do j = 0, m
        dd(0:m,j) = u(j)
        dd(m+1:2*m+1, j) = u(j+m+1)
    enddo
    do j = 1, m
        dd(m-j+1:m, j) = dd(m-j+2:m+1, j-1) - dd(m-j+1:m, j-1)
    enddo
    do j = m+1, 2*m+1
        dd(0:2*m+1-j, j) = dd(1:2*m+1-j+1, j-1) - dd(0:2*m+1-j, j-1)
    enddo

    cofs = dd(0,:)


    ! Dual Vander solve
    do i = 2*m,0,-1
        do j = i, 2*m
            if (i < m+1) then
                cofs(j) = cofs(j) + 0.5_dp * cofs(j+1)
            else
                cofs(j) = cofs(j) - 0.5_dp * cofs(j+1)
            endif
        enddo
    enddo
end subroutine interpolate

    subroutine interpolate_2D(mp, u, scofsmp)
       implicit none
       integer, intent(in) :: mp
       real(kind=dp), dimension(0:mp, 0:mp, 0:1, 0:1), intent(in) :: u
       integer :: i, j, k, idx, idy
       real(kind=dp), dimension(0:2*mp+1,0:2*mp+1) :: scofsmp, smcofsmp
       real(kind=dp), dimension(0:2*mp+1) :: mpdata

       i = 0; j = 0;

       do idy = 0,mp
       mpdata(0:mp)        = u(:, idy, i  ,j)
       mpdata(mp+1:2*mp+1) = u(:, idy, i+1,j)
       call interpolate(mp, mpdata, smcofsmp(:,idy))

       end do
       do idy = 0,mp
       mpdata(0:mp)        = u(:,idy,i,   j+1)
       mpdata(mp+1:2*mp+1) = u(:,idy,i+1, j+1)
       call interpolate(mp, mpdata, smcofsmp(:,mp+1+idy))
       end do

       do idx = 0,2*mp+1
       call interpolate(mp, smcofsmp(idx,:), scofsmp(idx,:))
       end do
    end subroutine

    subroutine eval_g(mp, z, h, ord1, ord2, u, val, d)
       implicit none
       integer, intent(in) :: mp, d, ord1, ord2
       real(kind=dp), dimension(2), intent(in) :: z
       real(kind=dp), dimension(2), intent(in) :: h
       real(kind=dp), dimension(0:mp, 0:mp, 0:1, 0:1, d), intent(in) :: u
       real(kind=dp), intent(out) :: val(0:ord1, 0:ord2)
       integer :: i, j, k, idx, idy
       real(kind=dp), dimension(0:2*mp+1,0:2*mp+1, d) :: scofsmp
       real(kind=dp), dimension(0:2*mp+2,0:2*mp+2) :: icofs
       real(kind=dp) :: zx, zy, xfact, yfact
       do i = 1,d
         call interpolate_2D(mp, u(:,:,:,:,i), scofsmp(:,:,i))
       end do

       icofs = 0.d0
       do idy = 0, 2*mp+1
         do idx = 0, 2*mp+1
           icofs(idx+1, idy) =  icofs(idx+1, idy) + scofsmp(idx, idy, 2) / (h(1)**idx * h(2)**idy) / (idx + idy + 1)
           icofs(idx, idy+1) =  icofs(idx, idy+1) + scofsmp(idx, idy, 1) / (h(1)**idx * h(2)**idy) / (idx + idy + 1)
         end do
       end do

       do idy = 0, ord2
       do idx = 0, ord1
       val(idx, idy) = 0._dp
       zy = 1._dp
       do j = idy, 2*mp+1
       zx = 1._dp
       do i = idx, 2*mp+1 ! - j   ! for total degree etienne
       xfact = 1._dp
       yfact = 1._dp
       do k = i-idx+1,i
       xfact = xfact * dble(k)
       enddo
       do k = j-idy+1,j
       yfact = yfact * dble(k)
       enddo
       val(idx, idy) = val(idx, idy) + xfact*yfact*icofs(i,j) * zx * zy
       zx = zx * z(1)
       enddo
       zy = zy * z(2)
       enddo
       enddo
       enddo
     end subroutine eval_g

subroutine track_hermite(mh,xs0)
implicit none
type(hermite), intent(in) :: mh
type(probe), intent(inout) :: xs0
real(dp) x(6), z0(6), xf(1,6) ,  val(0:0,0:0)
integer blk(6),nd2
 
!write(6,*) mh%gen
if(mh%gen==0) then
blk=0
nd2=2
 
  blk(1:nd2) = floor((xs0%x(1:nd2)-mh%x0(1:nd2,0,0)) * mh%n / (mh%b(:,2)))
 
  xs0%u= blk(1)>mh%n.or.blk(1)<-mh%n.or.blk(2)>mh%n.or.blk(2)<-mh%n
  if(xs0%u) then
    write(6,*) " grid 1"
    return
  endif
  z0(1:nd2) = (xs0%x(1:nd2) - mh%h*( 0.5d0)-mh%x0(1:nd2,blk(1),blk(2)) ) / mh%h
  
  !call  eval(mh%noh, 1, z0,mh%h, 0, 6, mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,1), xf(1,1))
  !call  eval(mh%noh, 1, z0,mh%h,0, 6, mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,2), xf(1,2))


   call  eval_new(mh%noh,  z0,mh%h, 0, 1, mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,1), val)
  xs0%x(1) = val(0,0)
   call  eval_new(mh%noh,  z0,mh%h, 0, 1, mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,2), val)
  xs0%x(2) = val(0,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Oleksii's   Hermite   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mp=noh
! z(1:2)
! h(2)
! ord1  ord2  (0 no derivatives)
!  U is %he
! val is array of size (0:ord1,0:ord2)

!subroutine eval_new(mp, z, h, ord1, ord2,  u, val)




  
 elseif(mh%gen==1) then
 stop 666
 elseif(mh%gen==2) then
   call track_hermite_linear(mh,xs0)! no Hermite
   call track_hermite_invert2(mh,xs0)   ! Hermite
 endif
end subroutine track_hermite



subroutine track_hermite_invert2(mh,xs0)
implicit none
type(hermite), intent(in) :: mh
type(probe), intent(inout) :: xs0
real(dp) z(1:15), z0(6), xf(1,6),dz(3),normb,norm ,pit,xi,val(0:1,0:1) 
integer blk(6),nd2
integer i
blk=0
xi=xs0%x(1)
pit=xs0%x(2)
!!! guess for z(1)
z(1)=xs0%x(1) 
z(2)=xs0%x(2) 
normb=1.d38
do i=1, mh%maxite
nd2=2
xf=0
  blk(1:nd2) = floor((z(1:nd2)-mh%x0(1:nd2,0,0)) * mh%n / (mh%b(:,2)))
 
   xs0%u= blk(1)>mh%n.or.blk(1)<-mh%n.or.blk(2)>mh%n.or.blk(2)<-mh%n
  if(xs0%u) then
    write(6,*) " grid 1",z(1:2)

    return
  endif
  z0(1:nd2) = (z(1:nd2) - mh%h*( 0.5d0)-mh%x0(1:nd2,blk(1),blk(2)) ) !/ mh%h
  
!!! value and derivative in the original not scaled 
! xf(1,1)  value
! xf(1,2)=0
! xf(1,3)=derivative
!subroutine eval_new(mp, z, h, ord1, ord2,  u, val)

 ! call  eval_new(mh%noh, z0, mh%h, 1,1,  mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,1), val)
   call  eval_g(mh%noh, z0, mh%h, 1, 1, mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,1:2), val, 2)
 dz=0
dz(1)=(xs0%x(1)-val(0,1))/val(1,1)
z(1)=z(1)+dz(1)
norm=abs(dz(1))  !+abs(dz(2))+abs(dz(3))
 
if(norm>mh%eps) then
  normb=norm
else
 if(norm>=normb) then
  xs0%x(1)=z(1)
z(2)= xs0%x(2)
 ! xf(1,2) in position 2
! xf(1,1) still preserved
  blk(1:nd2) = floor((z(1:nd2)-mh%x0(1:nd2,0,0)) * mh%n / (mh%b(:,2)))

   xs0%u= blk(1)>mh%n.or.blk(1)<-mh%n.or.blk(2)>mh%n.or.blk(2)<-mh%n
  if(xs0%u) then
    write(6,*) " grid 2"

    return
  endif
  z0(1:nd2) = (z(1:nd2) - mh%h*( 0.5d0)-mh%x0(1:nd2,blk(1),blk(2)) ) !/ mh%h

   call  eval_g(mh%noh, z0, mh%h, 1, 1, mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,1:2), val, 2)

  xs0%x(2)=val(1,0)
 !  z(1)=xs0%x(1)
 ! z(2)=pit
!  blk(1:nd2) = floor((z(1:nd2)-mh%x0(1:nd2,0,0)) * mh%n / (mh%b(:,2)))

 !  xs0%u= blk(1)>mh%n.or.blk(1)<-mh%n.or.blk(2)>mh%n.or.blk(2)<-mh%n
 ! z0(1,1:nd2) = (z(1:nd2) - mh%h*( 0.5d0)-mh%x0(1:nd2,blk(1),blk(2)) ) / mh%h

!  call  eval(mh%noh, 1, z0, mh%h,0, 6,  mh%he(:,:,blk(1):blk(1)+1,blk(2):blk(2)+1,1), xf(1,1))
 ! write(6,*) xf(1,1)
 ! write(6,*) xi
  return
 endif
  normb=norm
endif
 if(i>mh%maxite) then
   xs0%u=.true.
   write(6,*) " used too many iterations "
 endif
enddo

end subroutine track_hermite_invert2


  subroutine alloc_hermite(mh)
    implicit none
    type(hermite), intent(inout) :: mh
    integer i,j,k
    type(probe_8) xs
    type(probe)xs0
    type(c_damap) a,c_map,id
    type(c_normal_form) c_n
    mh%no=2*mh%noh
    allocate(mh%p(-mh%n:mh%n,-mh%n:mh%n))
    allocate(mh%ms(-mh%n:mh%n,-mh%n:mh%n))
    allocate(mh%x0(1:6,-mh%n:mh%n,-mh%n:mh%n))
    allocate(mh%he(0:mh%noh,0:mh%noh,-mh%n:mh%n,-mh%n:mh%n,2))
    mh%x0=0
    mh%he=0
    mh%f=0

    CALL FIND_ORBIT_x(mh%r,mh%f,mh%STATE,c_1d_5,fibre1=mh%pos)  ! (3)
    write(6,*) check_stable
    write(6,*) mh%f

    call init_all(mh%STATE,1,0)
    call alloc(a,c_map,id)
    call alloc(xs); call alloc(c_n);
    xs0=mh%f
    id=1
    xs=xs0+id
    CALL propagate(mh%r,XS,mh%STATE,FIBRE1=mh%pos)  ! (4)


    c_map=xs
    call c_normal(c_map,c_n)  ! (6)
    mh%a=c_n%a_t
    mh%ai=c_n%a_t**(-1)
    mh%m=c_map
    mh%mi=c_map**(-1)
    call kill(a,c_map,id)
    call kill(xs);
    call kill(c_n)

    do i=1,2
    mh%h(i)=(mh%b(i,2)-mh%b(i,1))/2/mh%n
    enddo

    do i=-mh%n,mh%n
    do  j = -mh%n,mh%n
    call alloc(mh%P(i,j))
 

    mh%x0(1,i,j)=mh%f(1)+i*mh%h(1)
    mh%x0(2,i,j)=mh%f(2)+j*mh%h(2)

    enddo
    enddo

  end subroutine alloc_hermite

  subroutine compute_hermite(mh)
    implicit none
    type(hermite), intent(inout) :: mh
    type(probe)xs0
    real(dp) x0(3),xf(6),cu
    type(c_damap) a,c_map,id
    type(damap) ms,idr
    integer i,j,k,js(6)
    type(c_normal_form) c_n

    call init_all(mh%STATE,mh%no,0)
    call alloc(a,c_map,id)
    call alloc(c_n);
    
    if(mh%gen==0) then
    id=1

    do i=1,c_%nd2
    id%v(i)=id%v(i)*mh%h(i)
    enddo

    do i=-mh%n,mh%n
    do  j = -mh%n,mh%n
    xs0=mh%x0(1:6,i,j)
    mh%p(i,j)=xs0+id

    CALL propagate(mh%r,mh%p(i,j),mh%STATE,FIBRE1=mh%pos)  ! (4)
 
    enddo
    enddo
     
    elseif(mh%gen==1) then
      call alloc(ms,idr)

    id=1
    idr=1
    do i=1,c_%nd2
      idr%v(i)=idr%v(i)*mh%h(i)
     enddo
 

    do i=-mh%n,mh%n
    do  j = -mh%n,mh%n
    CALL compute_partially_inverted_location(mh,i,j,x0)

    xs0=mh%x0(1:6,i,j)
    xs0%x(1)= x0(1)

    mh%p(i,j)=xs0+id
 
    call track_hermite_linear_inv_8(mh,mh%p(i,j))
    CALL propagate(mh%r,mh%p(i,j),mh%STATE,FIBRE1=mh%pos)  ! (4)
 !    mh%xf(1:6,i,j)=mh%p(i,j)%x
     
    ms=mh%p(i,j)%x    ! compute a damap
     cu=mh%p(i,j)%x(1)

 
    do k=1,c_%nd2
     ms%v(k)=ms%v(k)-(ms%v(k).sub.'0')
    enddo
 
      js=0
     js(1)=1;js(3)=1;js(5)=1;  
     
     ms=ms**(js(1:c_%nd2))  !partial inversion

    ms=ms*idr   ! scales

 

     do k=1,c_%nd2
       mh%p(i,j)%x(k)%t=ms%v(k) +mh%x0(k,i,j)
     enddo
 

    enddo
    enddo

      call kill(ms,idr)
    else
! gen = 2 

  mh%he=0
      call alloc(ms,idr)

    id=1
    idr=1
    do i=1,c_%nd2
      idr%v(i)=idr%v(i)*mh%h(i)
     enddo
 

    do i=-mh%n,mh%n
    do  j =-mh%n,mh%n
    call alloc(mh%ms(i,j))
    CALL compute_partially_inverted_location(mh,i,j,x0)

    xs0=mh%x0(1:6,i,j)
    xs0%x(1)= x0(1)
xf(1)=x0(1)

    mh%p(i,j)=xs0+id
 
    call track_hermite_linear_inv_8(mh,mh%p(i,j))
    CALL propagate(mh%r,mh%p(i,j),mh%STATE,FIBRE1=mh%pos)  ! (4)

xf(2)=mh%p(i,j)%x(2)
     
    ms=mh%p(i,j)%x    ! compute a damap

 
    do k=1,c_%nd2
     ms%v(k)=ms%v(k)-(ms%v(k).sub.'0')  
    enddo
 
      js=0
     js(1)=1;js(3)=1;js(5)=1;  
          
 
     ms=ms**(js(1:c_%nd2))  !partial inversion
 
      mh%ms(i,j)=ms


    do k=1,c_%nd2
     ms%v(k)=ms%v(k)+xf(k)
    enddo
     call fill_hermite_gen(mh,i,j,ms)

    call kill(mh%p(i,j))
 
      
 
 
 
 
     enddo
    enddo

      call kill(ms,idr)

      deallocate(mh%p)
    endif
    call kill(a,c_map,id)
    call kill(c_n)

  end subroutine compute_hermite

  subroutine compute_partially_inverted_location(mh,i,j,xf)
    implicit none
    type(hermite), intent(inout) :: mh
    integer, intent(in) :: i,j
    type(probe)xsf,xs1
    real(dp) mat(3,3) 
     real(dp) :: eps, xf(1:3),xfd(1:3), xfdm(1:3), x0(1:6),eps2,norm,normb,der(3),cu
    integer k
     
     eps=mh%h(1)/mh%nint/100
    eps2=eps
    x0=0
    x0(1:2)=mh%x0(1:2,i,j)
 
    normb=1.d8
     do k=1,mh%maxite   



    xs1=x0
    xs1%x(1)=xs1%x(1)-eps
 call track_hermite_linear_inv(mh,xs1)

    CALL propagate(mh%r,xs1,mh%STATE,FIBRE1=mh%pos)  ! (4)
    xfdm(1)=xs1%x(1)

    xs1=x0
    xs1%x(1)=xs1%x(1)+eps
 call track_hermite_linear_inv(mh,xs1)
    CALL propagate(mh%r,xs1,mh%STATE,FIBRE1=mh%pos)  ! (4)
     xfd(1)=xs1%x(1)

    xs1=x0
 call track_hermite_linear_inv(mh,xs1)
    CALL propagate(mh%r,xs1,mh%STATE,FIBRE1=mh%pos)  ! (4)
    xf(1)=xs1%x(1)


    mat(1,1)= (xfd(1)-xfdm(1))/eps/2
cu=xf(1)
  !!!!!!!!!!!!!!!!   invert mat

    mat(1,1)=1.d0/mat(1,1)
     der(1)=mat(1,1)*(mh%x0(1,i,j)-xf(1)) 

    x0(1)=x0(1)+der(1)
     norm=abs(der(1))
 
     if(norm> eps2) then
      normb=norm
     else
      if(normb<=norm) then
       xf(1) =x0(1)

        goto 111
       endif
         normb=norm
     endif


     
     enddo 
        write(6,'(4(1x,G21.14))') mh%x0(1,i,j),cu,der(1),eps2
       write(6,'(2(1x,G21.14))') norm,normb
     write(6,*) " failure ",i,j,k,x0(1)-mh%f(1)

stop


111 continue 
       write(6,'(4(1x,G21.14))') mh%x0(1,i,j),cu,der(1),eps2
       write(6,'(2(1x,G21.14))') norm,normb
     write(6,*) " success ",i,j,k,x0(1)-mh%f(1)
  
  end subroutine compute_partially_inverted_location

  subroutine track_hermite_linear_inv(mh,xs1)
  implicit none
type(probe) xs1
type(hermite) mh
     if(mh%linear) then
        xs1%x=xs1%x-mh%f
         xs1%x = matmul(mh%mi,xs1%x)
        xs1%x=  xs1%x+mh%f
       endif
end   subroutine track_hermite_linear_inv

  subroutine track_hermite_linear(mh,xs1)
  implicit none
type(probe) xs1
type(hermite) mh
     if(mh%linear) then
        xs1%x=xs1%x-mh%f
         xs1%x = matmul(mh%m,xs1%x)
        xs1%x=  xs1%x+mh%f
       endif
end   subroutine track_hermite_linear

  subroutine track_hermite_linear_inv_8(mh,xs1)
  implicit none
type(probe_8) xs1,xs0
type(hermite) mh
integer i,j,k

     if(mh%linear) then
call alloc(xs0)
       do i=1,c_%nd2
        xs1%x(i)=xs1%x(i)-mh%f(i)
        xs0%x(i)=0.0_dp
       enddo
       do i=1,c_%nd2
       do j=1,c_%nd2
           xs0%x(i)= mh%mi(i,j)*xs1%x(j)+xs0%x(i)
       enddo
       enddo
         
        xs1%x=xs0%x

       do i=1,c_%nd2
        xs1%x(i)=xs1%x(i)+mh%f(i)
       enddo
call kill(xs0)
       endif

end   subroutine track_hermite_linear_inv_8

  subroutine fill_hermite(mh)
    implicit none
    type(hermite) mh
    real(dp) x
    integer i,j,k,n,l
    integer, allocatable :: jc(:)
     if(mh%gen==2) return
    allocate(jc(c_%nd2))

    do i=-mh%n,mh%n
    do  j = -mh%n,mh%n
     do k=1,2
     call taylor_cycle(mh%P(i,j)%x(k)%t,size=n)
     do l=1,n
     call taylor_cycle(mh%P(i,j)%x(k)%t,ii=l,value=x,j=jc)
      if(jc(1)>mh%noh) cycle
      if(jc(2)>mh%noh) cycle
     mh%he(jc(1),jc(2),i,j,k)=x
    enddo
    enddo
      call kill(mh%P(i,j))
    enddo
    enddo
   
    deallocate(mh%p)
    deallocate(jc)

  end subroutine fill_hermite

  subroutine fill_hermite_gen(mh,i,j,ms)
    implicit none
    type(hermite) mh
    type(damap), intent(in):: ms 
    type(damap) c
    real(dp) x
    real(dp) xf(6)
    integer i,j,k,n,l
    integer, allocatable :: jc(:)


    allocate(jc(c_%nd2))



    do k=1,c_%nd2
     call taylor_cycle(ms%v(k),size=n)
    do l=1,n
     call taylor_cycle(ms%v(k),ii=l,value=x,j=jc)

     if(jc(1)<=mh%noh.and.jc(2)<=mh%noh) mh%he(jc(1),jc(2),i,j,k)=x * mh%h(1) ** jc(1) * mh%h(2) ** jc(2)

    enddo
    enddo
 

    deallocate(jc)

  end subroutine fill_hermite_gen

  subroutine kill_hermite(mh)
    implicit none
    type(hermite) mh
    integer i,j
    mh%n=0
    mh%h=0
    mh%a=0
    mh%f=0
    mh%b=0
    mh%pos=0
    mh%no=0
    if(associated(mh%p)) then
       do i=-mh%n,mh%n
       do  j = -mh%n,mh%n
          call kill(mh%P(i,j))
       enddo
       enddo
      deallocate(mh%p)
    endif
    if(associated(mh%ms)) then
       do i=-mh%n,mh%n
       do  j = -mh%n,mh%n
          call kill(mh%ms(i,j))
       enddo
       enddo
      deallocate(mh%ms)
    endif
    if(associated(mh%x0)) deallocate(mh%x0)
    if(associated(mh%he)) deallocate(mh%he)
    mh%r=> null()

 
  end subroutine kill_hermite

 



end module pointer_lattice




subroutine read_ptc_command77(ptc_fichier)
  use pointer_lattice
  implicit none
  character(*) ptc_fichier
  integer m

  if(ptc_fichier(1:len_trim(ptc_fichier))=='CPP') then
     call change_default_tpsa(1)
     return
  endif
  if(ptc_fichier(1:len_trim(ptc_fichier))=='FORTRAN') then
     call change_default_tpsa(2)
     return
  endif

  call kanalnummer(m)

  open(unit=m,file=ptc_fichier,status='OLD',err=2001)
  close(m)
  call read_ptc_command(ptc_fichier)
  return
2001 continue

  write(6,*) " Warning: command file does not exit "

end  subroutine read_ptc_command77



!Piotr says:
!I have implemented for you MAD-X command ptc_script.
!It has one parameter named "file". Hence, you call sth like
!ptc_script, file="directptctweaks.ptc"
!It executes subroutine execscript in file madx_ptc_script.f90.

subroutine gino_ptc_command77(gino_command)
  use pointer_lattice
  implicit none
  character(*) gino_command

  !  Etienne puts garbage here
  ! print*,"Fuck it works! ",gino_command

  call context(gino_command)
  call call_gino(gino_command)

end  subroutine gino_ptc_command77


subroutine read_mad_command77(ptc_fichier)
  use pointer_lattice
!  use pointer_lattice_frank
  implicit none
  character(*) ptc_fichier
  integer m

  call kanalnummer(m)

  open(unit=m,file=ptc_fichier,status='OLD',err=2001)
  close(m)
!  call read_mad_command(ptc_fichier)  ! inside pointer_lattice_frank
  return
2001 continue

  write(6,*) " Warning: mad command file does not exit "



end  subroutine read_mad_command77

subroutine my_user_routine1 
use madx_ptc_module
!use pointer_lattice
!use duan_zhe_map, probe_zhe=>probe,tree_element_zhe=>tree_element,dp_zhe=>dp, & 
!DEFAULT0_zhe=>DEFAULT0,TOTALPATH0_zhe=>TOTALPATH0,TIME0_zhe=>TIME0,ONLY_4d0_zhe=>ONLY_4d0,RADIATION0_zhe=>RADIATION0, &
!NOCAVITY0_zhe=>NOCAVITY0,FRINGE0_zhe=>FRINGE0,STOCHASTIC0_zhe=>STOCHASTIC0,ENVELOPE0_zhe=>ENVELOPE0, &
!DELTA0_zhe=>DELTA0,SPIN0_zhe=>SPIN0,MODULATION0_zhe=>MODULATION0,only_2d0_zhe=>only_2d0 , &
!INTERNAL_STATE_zhe=>INTERNAL_STATE
implicit none
 
type(layout),pointer :: fodo
type(c_normal_form) normal_form
integer pos,no,np,mf,i
real(dp) closed_orbit(6)
type(probe) xs0
type(probe_8) xs
type(c_damap) id,one_turn_map,a_cs,a0,a1,a2,a
type(c_taylor) x2div2,x2div2_f,phase(3),phase_one_turn_map(3)
type(internal_state) state
type(fibre), pointer :: f



 
!!!! reading the flat file produced by MAD-X
!call ptc_ini_no_append
call read_lattice_append(M_U,"../../files_for_cas/guido_fodo/flat.txt")
!call read_lattice_append(M_U,"C:\msys64\home\Etienne\MAD-X\files_for_cas\guido_fodo\flat.txt")
!call read_lattice_append(M_U,"guido_fodo_ubuntu.txt")
fodo=>m_u%end


call in_bmad_units   ! units similar to MAD-X
 

! finds the closed orbit at position 1 (should be (0,0,0,0,0,0))
pos=1
closed_orbit=0.d0;  
                                                
call find_orbit_x(fodo,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)  
 
state=nocavity   ! state that produces map with delta dependence
no=3
np=0
call init_all(state,no,np)
 

!   create these TPSA objects
call alloc(id,one_turn_map,a_cs,a0,a1,a2,a)
call alloc(xs)
call alloc(normal_form)
call alloc(x2div2,x2div2_f)
call alloc(phase)
call alloc(phase_one_turn_map)

xs0=closed_orbit   ! xs0 contains orbit and spin 
id=1   !    identity map
xs=id+xs0   !  xs is a probe_8 which can become a Taylor series
 
 
call propagate(fodo,xs,state,fibre1=pos) ! computes one turn map around closed orbit
 
one_turn_map=xs
 
 
call c_normal(one_turn_map,normal_form,phase=phase_one_turn_map)  ! one_turn_map= normal_form%atot o  rotation o  normal_form%atot^-1
 
call c_canonise(normal_form%atot,a_cs,a0,a1,a2) 

xs =  a_cs +xs0

call kanalnummer(mf,"twiss_from_guido.txt")

phase=0.0_dp
x2div2=2.0_dp*(1.0_dp.cmono.1)**2
call C_AVERAGE(x2div2,a1,x2div2_f)  

write(mf,*) " Phase advance in x "
call clean(phase,phase,prec=1.d-10)
call print(phase(1),mf)
write(mf,*) " Beginning of lattice "
call print(x2div2_f,mf)

f=>fodo%start
do i=1,fodo%n
 call propagate(fodo,xs,state,fibre1=i,fibre2=i+1)

write(mf,*) " end of Magnet ",f%mag%name

 a=xs  ! creates tracked canonical transformation
 xs0=xs
call c_canonise(a,a_cs,a0,a1,a2,phase) ;call clean(a1,a1,prec=1.d-10);
 
call C_AVERAGE(x2div2,a1,x2div2_f)  

write(mf,*) " Phase advance in x "
call clean(phase,phase,prec=1.d-10)
call print(phase(1),mf)
write(mf,*) " 2(x^2> ~ beta"
call print(x2div2_f,mf)

xs=xs0+a_cs

f=>f%next
enddo


write(mf,*) " Tune in x from one turn map"
call clean(phase_one_turn_map,phase_one_turn_map,prec=1.d-10)
call print(phase_one_turn_map(1),mf)

close(mf)

call kill(id,one_turn_map,a_cs,a0,a1,a2,a)
call kill(xs)
call kill(normal_form)
call kill(x2div2,x2div2_f)
call kill(phase)
call kill(phase_one_turn_map)


end  subroutine my_user_routine1

