module pointer_lattice
  use madx_ptc_module !etienne, my_state_mad=>my_state, my_ring_mad=>my_ring
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
  character(255) :: filename_ap = "Tracking.txt"
  integer last_npara 
  integer :: i_layout=0,i_layout_t=1
  integer my_lost_position
  private thin
  real(dp) thin
  !  BEAM STUFF
  REAL(DP) SIG(6) 
  REAL(DP) ait(6,6)
  private equal_spinor_fourier_spinor_fourier

  INTEGER :: N_BEAM=0,USE_BEAM=1
  logical, private :: m_u_t = .true.
  TYPE(REAL_8),private :: Y(6)
  TYPE(DAMAP),PRIVATE :: ID
  !  integer nd2,npara
  !  PRIVATE POWER_CAVITY,RADIA
  ! stuff from my fortran
  type(internal_state), target:: etat
  integer,target:: START ,FIN,ORDER,np,start_t
  real(dp),target:: xfix(6) ,DELT0
  integer :: logs_exp=30, num_iter = 20
  !logical :: absolute = .false.
 !private ind1,ind2,ipos,ireal,a_f,a_f0,yfit
  integer, allocatable:: ind1(:),ind2(:),ipos(:)
  logical, allocatable :: ireal(:)
  real(dp), allocatable :: a_f(:),a_f0(:),yfit(:),dyfit(:,:)
  integer sizeind1
  logical :: onefunc = .true.,skipzero=.false.,skipcomplex=.true.
 
   INTERFACE assignment (=)
        MODULE PROCEDURE     equal_spinor_fourier_spinor_fourier
        MODULE PROCEDURE     equal_matrix_fourier_matrix_fourier
   end INTERFACE 
   
  INTERFACE SCRIPT
     MODULE PROCEDURE read_ptc_command
  END INTERFACE

!!  new stuff on non-perturbative !

type  vector_field_fourier 
     real(dp) fix(6)    !   closed orbit of map
     integer :: ns(3)=0   ! integer for Fourier transform
     real(dp), pointer :: mr(:,:,:,:,:) =>null()   ! spin matrices produced by code (i,j,k,1:3,1:3)
     real(dp), pointer ::  x_i(:,:,:,:) =>null()   ! starting position in phase  x_i(i,j,k,1:6)=r%x(1:6)
     real(dp), pointer :: phis(:,:,:,:) =>null()   ! %phis(i,j,k,2)=j*dphi2
     type(spinor), pointer :: sp(:,:,:)            ! sp(:i,j,k)   spinor for all the matrices mr  
     integer n1,n2,n3,nd                           ! # fourier modes -n1:n1, etc... nd=degree of freedom
     real(dp)  mu(3),muf(3),em(3)                  ! tune and initial emitances
     complex(dp),  DIMENSION(:,:,:,:), POINTER :: f  !  vector field expansion f(1:3,-n1:n1,-n2:n2,-n3:n3)
end  type vector_field_fourier
   type(vector_field_fourier) af

type  spinor_fourier 
     integer n1,n2,n3
     complex(dp),  DIMENSION(:,:,:,:), POINTER :: s
end  type spinor_fourier 

type  matrix_fourier 
      real(dp) muf(3)
     type(spinor_fourier) v(3)
end  type matrix_fourier 

  type  explogs  
     integer n1,n2
     complex(dp),  DIMENSION(:,:,:), POINTER :: h
  end  type explogs


  type  logs 
     integer  m(3),ms  
     integer  ns,no
     type(explogs) h,a,n
     type(explogs), pointer ::  af(:)
     real(dp), pointer :: as(:,:,:)
     type(spinor), pointer :: sp(:,:)
     real(dp), pointer :: s(:,:,:,:)  !   spin matrices
     real(dp) em(2),mu(2),fix(6)
     real(dp), pointer :: x_i(:,:,:),phis(:,:,:)
  end  type logs


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
    integer nstep(2),i1,i2,I3,n_bessel
    LOGICAL(LP) STRAIGHT,skip,fixp,skipcav,fact
    ! end
    ! TRACK 4D NORMALIZED
    INTEGER POS,NTURN,resmax
    real(dp) EMIT(6),APER(2),emit0(2),sca
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
    REAL(DP) X(6),DT(3),x_ref(6),sc,NLAM,A1,B1,HPHA,B_TESLA,CUR1,CUR2
    REAL(DP)VOLT,PHASE
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
    integer icnmin,icnmax,n_ac,inode !,n_coeff
    logical :: log_estate=.true.,longprintt,onemap
    integer :: mftune=6,nc
    real(dp), allocatable :: tc(:)
    type(integration_node), pointer  :: t

     longprintt=longprint 
     longprint=.true.
    if(log_estate) then
       nullify(my_estate)
!       nullify(my_old_state)
       log_estate=.false.
    endif

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
!       case('SETORBITPHASORFREQUENCY')
!          read(mf,*) xsm%ac%om
!          xsm0%ac%om=xsm%ac%om
!          ptc_node_old=-1
!          first_particle=my_false
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
       case('+NOCAVITY')
          my_estate=my_estate+NOCAVITY0
       case('-NOCAVITY')
          my_estate=my_estate-NOCAVITY0
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
       case('-SPIN')
          my_estate=my_estate-SPIN0
          !  case('+EXACTMIS')
          !     my_estate=my_estate+EXACTMIS0
          !  case('-EXACTMIS')
          !     my_estate=my_estate-EXACTMIS0

          !       case('DEFAULTTPSA')
          !       read(mf,*) default_tpsa
          !       if(default_tpsa) then
          !        write(6,*) " Default TPSA is Chinese "
          !       else
          !        write(6,*) " Default TPSA is Germanic "
          !       endif

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
                   P%MAG%slow_ac=.true.
                   P%MAGP%slow_ac=.true.

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
          READ(MF,*) NAME

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
                   P%MAG%slow_ac=.true.
                   P%MAGP%slow_ac=.true.

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
       case('MAKEONETURNMAP','TRACKWITHONETURNMAP')
          READ(MF,*) i1,i3  ! position
          READ(MF,*) I2  ! ORDER OF THE MAP
          READ(MF,*) fixp,fact,noca  !  SYMPLECTIC , factored
    !      READ(MF,*) filename
          if(.not.associated(my_ering%t)) call make_node_layout(my_ering)
          
           p=>my_ering%start
           
           do ii=1,i1
             f1=>p
            p=>p%next
           enddo
    
             f2=>f1
             x_ref=0.0_dp
             call MOVE_TO_LAYOUT_I(m_u,my_fring,i3)

             if(associated(my_ering,my_fring))then 
                ft=>f1
             else
                ft=>my_fring%start       
             endif
             call FIND_ORBIT_x(x_ref,time0,1.d-7,fibre1=f1)
         !    name_root=filename
        !     call context(name_root)
       !      if(name_root(1:2)=='NO') then
              call fill_tree_element_line(f1,f2,ft,i2,x_ref,fact,nocav=noca)
       !      else
        !      call fill_tree_element_line(f1,f2,ft,i2,x_ref,fact,nocav=noca,file=filename)
      !       endif
                    ft%mag%forward(3)%symptrack=FIXP
                    ft%magP%forward(3)%symptrack=FIXP
                    ft%mag%do1mapf=.true.
                    ft%magp%do1mapf=.true.
                    ft%mag%filef="one_turn_map.txt"
             if(associated(ft,f1))then
                f2=>f1%next
               do i1=1,my_ering%n
                if(associated(f2,f1)) exit
                 f2%mag%skip_ptc_f=.true.
                 f2%magp%skip_ptc_f=.true.
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
  
              endif
             else
              icnmin=icnmin+1
             ENDIF


             p=>p%next
          enddo
         write(6,*) icnmax, " changed into Taylor maps "
         write(6,*) icnmin, " markers "
         write(6,*) my_ering%N, " total number of fibres "
       case('REMOVEALLMAP')

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
       case('FITTUNECOUPLING')
          read(mf,*) epsf
          read(mf,*) targ_tune
          call lattice_linear_res_gmap(my_ering,my_estate,epsf,pol_,NPOL,targ_tune,NP)
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
       case('PRINTSTATE')
          CALL PRINT(my_estate,6)
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

 SUBROUTINE radia_new(R,loc,estate,FILE1,fix,em,sij,sijr,tune,damping)
    implicit none
    TYPE(LAYOUT) R

    REAL(DP) X(6),m,energy,deltap
    CHARACTER(*), optional :: FILE1
    type(c_damap)  Id,a0,a_cs
    type(c_normal_form) normal
    integer  i,j 
    real(dp), optional :: fix(6), em(3),sij(6,6),tune(3),damping(3)
    complex(dp), optional :: sijr(6,6)   
    TYPE(INTERNAL_STATE) state
    TYPE(INTERNAL_STATE), target :: estate
    integer loc,mf1
    type(fibre), pointer :: p
    type(probe) xs0
    type(probe_8) xs
    character*48 fmd,fmd1

    if(present(FILE1)) then
    call kanalnummer(mf1)
    open(mf1,file=FILE1)
    endif
fmd= '(a12,1X,a3,I1,a3,i1,a4,D18.11,1x,D18.11)'
fmd1='(1X,a3,I1,a3,i1,a4,2(D18.11,1x),(f10.3,1x),a2)'



    state=(estate-nocavity0)+radiation0
    x=0.d0

    CALL FIND_ORBIT_x(R,X,STATE,1.0e-8_dp,fibre1=loc)
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
    CALL INIT(state,1,0)
    CALL ALLOC(NORMAL)
    CALL ALLOC(ID)
    call alloc(xs)

!    if(i1==1) normal%stochastic=my_true
    xs0=x
    ID=1
    xs=XS0+ID

    state=state+envelope0

    CALL TRACK_PROBE(r,xs,state, fibre1=loc)
     id=xs
   if(present(FILE1)) then
    write(mf1,*) " Full Map "    
     call print(id,mf1)
    write(mf1,*) " End of Full Map "
   endif
  call  c_normal(id,normal)             ! (8)

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
       do j=1,6
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

!    if(normal%stochastic) then
!       write(mf1,*) "Stochastic kicks"
!       write(mf1,*) normal%kick!

!       write(mf1,*)" Stochastic Transformation "
!       do i=1,6
!          do j=1,6
!             write(mf1,*) i,j,normal%STOCH(i,j)
!          enddo
!       enddo

!       write(mf1,*)" Inverse Stochastic Transformation "
!       do i=1,6
!          do j=1,6
!             write(mf1,*) i,j,normal%STOCH_inv(i,j)
!          enddo
!       enddo

!    endif



    CALL KILL(NORMAL)
    CALL KILL(ID)
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
          p%patch%b_t=-s
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

  
  !!!! new non perturbative  !!!
  
!!!!!!!!!!!!!!!!!!  new non-perturbative

subroutine make_matrix_fourier(af,a)
 implicit none
  type(vector_field_fourier)  af
 type(matrix_fourier) a
  type(spinor_fourier)  s
  integer i
  

  call alloc_spinor_fourier(s,h=af)
do i=1,3
 s%s(1:3,0,0,0)=0.d0
 s%s(i,0,0,0)=1.d0
 call exp_h_s(af,s,a%v(i))
enddo
  call kill_spinor_fourier(s)

end subroutine make_matrix_fourier

subroutine exp_h_m(af,a,b)
 implicit none
  type(vector_field_fourier)  af
 type(matrix_fourier) a,b
  type(spinor_fourier)  s
  integer i
  

  call alloc_spinor_fourier(s,h=af)

do i=1,3
 call exp_h_s(af,a%v(i),s)
 b%v(i)=s
enddo

  call kill_spinor_fourier(s)

end subroutine exp_h_m  

subroutine find_log_tpsa(a,af)
 implicit none
  type(vector_field_fourier)  af
 type(matrix_fourier) a,b
  type(spinor_fourier)  s
  type(vector_field_fourier)  bf
  real(dp)  norm,c
  integer i,mf
  
  call alloc_matrix_fourier(b,h=af)
  call alloc_vector_field_fourier(bf,h=af)

call norm_matrix_fourier(a,norm)

a%muf=0

b=a

b%v(1)%s(1,0,0,0)=b%v(1)%s(1,0,0,0)-1.0_dp
b%v(2)%s(2,0,0,0)=b%v(2)%s(2,0,0,0)-1.0_dp
b%v(3)%s(3,0,0,0)=b%v(3)%s(3,0,0,0)-1.0_dp

c=-1
do i=2,logs_exp*20
call  mul_matrix_fourier(a,b,b)
 a%v(1)%s=a%v(1)%s+(c/i)*b%v(1)%s
 a%v(2)%s=a%v(2)%s+(c/i)*b%v(2)%s
 a%v(3)%s=a%v(3)%s+(c/i)*b%v(3)%s
call norm_matrix_fourier(a,norm)
write(6,*) i,norm
c=-c
enddo

 call find_symetric_log(a,af)
call norm_vector_field_fourier(af,norm)
write(6,*) " symetric ",norm
 call find_asymetric_log(a,af)
write(6,*) " asymetric ",norm
call norm_vector_field_fourier(af,norm)

af%mu=0.0_dp
af%muf=0.0_dp


  call kill_matrix_fourier(b)
  call kill_vector_field_fourier(bf)

end subroutine find_log_tpsa

subroutine mul_matrix_fourier(a,b,c)
 implicit none
 type(matrix_fourier) a,b,c,t,bt

  integer i,j,k,n1,n2,n3
  integer i1,i2,i3,j1,j2,j3,k1,k2,k3

  call alloc_matrix_fourier(t,s=a%v(1))
  call alloc_matrix_fourier(bt,s=a%v(1))
 
n1=(size(a%v(1)%s,2)-1)/2
n2=(size(a%v(1)%s,3)-1)/2
n3=(size(a%v(1)%s,4)-1)/2

bt=a
call trans_mu_matrix(a%muf,bt)

 do i=1,3
 do j=1,3
 do k=1,3
 do j1=-n1,n1
 do k1=-n1,n1
     i1=j1+k1
    if(abs(i1)>n1) cycle
 do j2=-n2,n2
 do k2=-n2,n2
     i2=j2+k2
     if(abs(i2)>n2) cycle
 do j3=-n3,n3
 do k3=-n3,n3     
     i3=j3+k3
if(abs(i3)>n3) cycle


  t%v(k)%s(i,i1,i2,i3) = bt%v(j)%s(i,j1,j2,j3)*b%v(k)%s(j,k1,k2,k3) +  t%v(k)%s(i,i1,i2,i3)

 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo

c=t
  call kill_matrix_fourier(t)
  call kill_matrix_fourier(bt)

end subroutine mul_matrix_fourier


subroutine find_log(a,af)
 implicit none
  type(vector_field_fourier)  af
 type(matrix_fourier) a,b
  type(spinor_fourier)  s
  type(vector_field_fourier)  bf
  real(dp)  norm
  integer i,mf
  
call kanalnummer(mf,"junk.txt")
  call alloc_matrix_fourier(b,h=af)
  call alloc_vector_field_fourier(bf,h=af)

123 call norm_matrix_fourier(a,norm)
write(6,*) "1",norm
af%mu=-af%mu
af%f=-af%f
bf%f=0.d0
call exp_h_m(af,a,b)
af%mu=-af%mu
af%f=-af%f
call norm_matrix_fourier(b,norm)
write(6,*) "2",norm

call find_asymetric_log(b,bf)


!call find_log_tpsa(b,bf)

!stop

af%f=af%f+bf%f

!pause 

goto 123
  call kill_matrix_fourier(b)
  call kill_vector_field_fourier(bf)
close(mf)

end subroutine find_log  


subroutine norm_spinor_fourier(s,norm)
implicit none
type(spinor_fourier)  s
integer i,j,k,l
real(dp) norm

norm=0.0_dp

do i=1,3
do j=-s%n1,s%n1
do k=-s%n2,s%n2
do l=-s%n3,s%n3
 
norm=norm+abs(s%s(i,j,k,l))

enddo
enddo
enddo
enddo

end subroutine norm_spinor_fourier


subroutine norm_matrix_fourier(ma,norm)
implicit none
type(matrix_fourier)  ma
integer i,j,k,l
real(dp) norm,norm1

norm=0.0_dp

do i=1,3
 call norm_spinor_fourier(ma%v(i),norm1)
 norm=norm+norm1
enddo

norm=abs(norm-3)

end subroutine norm_matrix_fourier

subroutine norm_vector_field_fourier(af,norm)
implicit none
type(vector_field_fourier)  af
integer i,j,k,l
real(dp) norm,norm1


norm=0.0_dp

do i=1,3
do j=-af%n1,af%n1
do k=-af%n2,af%n2
do l=-af%n3,af%n3
 
norm=norm+abs(af%f(i,j,k,l))

enddo
enddo
enddo
enddo

norm=abs(norm-3)

end subroutine norm_vector_field_fourier


subroutine find_symetric_log(a,af)
implicit none
type(vector_field_fourier)  af
type(matrix_fourier)  a
integer n1,n2,n3

n1=af%n1
n2=af%n2
n3=af%n3
af%f=0.0_dp
 
af%f(3,-n1:n1,-n2:n2,-n3:n3) =(a%v(1)%s(2,-n1:n1,-n2:n2,-n3:n3)+a%v(2)%s(1,-n1:n1,-n2:n2,-n3:n3))/2.0_dp
af%f(1,-n1:n1,-n2:n2,-n3:n3) =(a%v(2)%s(3,-n1:n1,-n2:n2,-n3:n3)+a%v(3)%s(2,-n1:n1,-n2:n2,-n3:n3))/2.0_dp
af%f(2,-n1:n1,-n2:n2,-n3:n3) =(a%v(3)%s(1,-n1:n1,-n2:n2,-n3:n3)+a%v(1)%s(3,-n1:n1,-n2:n2,-n3:n3))/2.0_dp

end subroutine find_symetric_log

subroutine find_asymetric_log(a,af)
implicit none
type(vector_field_fourier)  af
type(matrix_fourier)  a
integer n1,n2,n3

n1=af%n1
n2=af%n2
n3=af%n3
af%f=0.0_dp
 
af%f(3,-n1:n1,-n2:n2,-n3:n3) =(a%v(1)%s(2,-n1:n1,-n2:n2,-n3:n3)-a%v(2)%s(1,-n1:n1,-n2:n2,-n3:n3))/2.0_dp
af%f(1,-n1:n1,-n2:n2,-n3:n3) =(a%v(2)%s(3,-n1:n1,-n2:n2,-n3:n3)-a%v(3)%s(2,-n1:n1,-n2:n2,-n3:n3))/2.0_dp
af%f(2,-n1:n1,-n2:n2,-n3:n3) =(a%v(3)%s(1,-n1:n1,-n2:n2,-n3:n3)-a%v(1)%s(3,-n1:n1,-n2:n2,-n3:n3))/2.0_dp

 
end subroutine find_asymetric_log


subroutine print_spinor_fourier(s,mf)
implicit none
type(spinor_fourier)  s
integer mf,i,j,k,l

do i=1,3
write(mf,*) " term ",i
do j=-s%n1,s%n1
do k=-s%n2,s%n2
do l=-s%n3,s%n3
write(mf,'(4(i4,1x),2(1x,E15.8))') i,j,k,l,s%s(i,j,k,l)

enddo
enddo
enddo
enddo


end subroutine print_spinor_fourier

subroutine print_vector_field(af,mf)
implicit none
type(vector_field_fourier)  af
integer mf,i,j,k,l

do i=1,3
write(mf,*) " term ",i

do j=-af%n1,af%n1
do k=-af%n2,af%n2
do l=-af%n3,af%n3

write(mf,'(4(i4,1x),2(1x,E15.8))') i,j,k,l,af%f(i,j,k,l)

enddo
enddo
enddo
enddo


end subroutine print_vector_field 

subroutine print_a_f0(mf)
implicit none
integer mf,i

do i=1,size(ind1)
 write(mf,'(L1,1x,3(i4,1x),1(1x,E15.8))') ireal(i),ipos(i),ind1(i),ind2(i),a_f0(i)
enddo

end subroutine print_a_f0
 
subroutine print_a_f(mf)
implicit none
integer mf,i

do i=1,size(ind1)
 write(mf,'(L1,1x,3(i4,1x),1(1x,E15.8))') ireal(i),ipos(i),ind1(i),ind2(i),a_f(i)
enddo

end subroutine print_a_f

subroutine alloc_matrix_fourier(a,n1,n2,n3,s,h)
 implicit none
 type(matrix_fourier) a
  type(spinor_fourier),optional :: s
  type(vector_field_fourier),optional :: h
 integer,optional :: n1,n2,n3
 integer i
a%muf=0.0_dp
do i=1,3
 call alloc_spinor_fourier(a%v(i),n1,n2,n3,s,h)
enddo

end subroutine alloc_matrix_fourier


subroutine kill_matrix_fourier(a)
 implicit none
 type(matrix_fourier) a
 integer i
a%muf=0.0_dp
do i=1,3
 call kill_spinor_fourier(a%v(i))
enddo

end subroutine kill_matrix_fourier


subroutine alloc_spinor_fourier(a,n1,n2,n3,s,h)
 implicit none
 type(spinor_fourier) a
  type(spinor_fourier),optional :: s
  type(vector_field_fourier),optional :: h
 integer,optional :: n1,n2,n3
 
if(present(n1)) then
a%n1=n1
a%n2=n2
a%n3=n3
elseif(present(h)) then
a%n1=h%n1
a%n2=h%n2
a%n3=h%n3    
elseif(present(s)) then   
a%n1=s%n1
a%n2=s%n2
a%n3=s%n3   
endif
    
allocate(a%s(3,-a%n1:a%n1,-a%n2:a%n2,-a%n3:a%n3) )
a%s=0.0_dp


end subroutine alloc_spinor_fourier


subroutine kill_spinor_fourier(a)
 implicit none
 type(spinor_fourier) a
a%n1=0
a%n2=0
a%n3=0
deallocate(a%s)
end subroutine kill_spinor_fourier

subroutine alloc_vector_field_fourier(a,n1,n2,n3,h,s)
 implicit none
 type(vector_field_fourier) a
type(spinor_fourier),optional :: s
type(vector_field_fourier),optional :: h
 integer,optional :: n1,n2,n3
integer i1,i2,k,i
 

if(present(n1)) then
a%n1=n1
a%n2=n2
a%n3=n3
elseif(present(h)) then
a%n1=h%n1
a%n2=h%n2
a%n3=h%n3    
elseif(present(s)) then   
a%n1=s%n1
a%n2=s%n2
a%n3=s%n3   
endif

allocate(a%f(3,-a%n1:a%n1,-a%n2:a%n2,-a%n3:a%n3))
a%f=0.0_dp
a%mu=0.0_dp  
a%muf=0.0_dp  
a%em=0.0_dp  


if(.not.allocated(ind1)) then
if(skipcomplex.and.skipzero) then
write(6,*) " skipcomplex and skipzero are true "
stop
endif

 k=1
 do  i1=-a%n1,a%n1
 do  i2=-a%n2,a%n2
 if(skipzero.and.i1==0.and.i2==0)goto 101
  k=k+1
 if(i1==0.and.i2==0)goto 101
 enddo

 enddo

101  continue
 k=k-1

if(skipcomplex) then
 allocate(ind1(6*k-3),ind2(6*k-3),a_f(6*k-3),a_f0(6*k-3),ipos(6*k-3),ireal(6*k-3))
else
 allocate(ind1(6*k),ind2(6*k),a_f(6*k),a_f0(6*k),ipos(6*k),ireal(6*k))
endif


a_f=0.0_dp
ireal=.true.
k=1
 do  i1=-a%n1,a%n1
 do  i2=-a%n2,a%n2
 if(skipzero.and.i1==0.and.i2==0)goto 102 
  ind1(k)=i1
  ind2(k)=i2
  k=k+1
 if(i1==0.and.i2==0)goto 102
 enddo

enddo

102  continue
k=k-1



if(skipcomplex) then
  ind1(k+1:2*k-1)=ind1(1:k-1)
  ind2(k+1:2*k-1)=ind2(1:k-1)
  ind1(2*k:3*k-1)=ind1(1:k)
  ind2(2*k:3*k-1)=ind2(1:k)
  ind1(3*k:4*k-2)=ind1(1:k-1)
  ind2(3*k:4*k-2)=ind2(1:k-1)
  ind1(4*k-1:5*k-2)=ind1(1:k)
  ind2(4*k-1:5*k-2)=ind2(1:k)
  ind1(5*k-1:6*k-3)=ind1(1:k-1)
  ind2(5*k-1:6*k-3)=ind2(1:k-1)
 ireal(k+1:2*k-1)=.false.
 ireal(3*k:4*k-2)=.false.
 ireal(5*k-1:6*k-3)=.false.
 ipos(1:2*k-1)=1
 ipos(2*k:4*k-2)=2
 ipos(4*k-1:6*k-3)=3
else
do i=1,k
  ind1(k+i)=ind1(i)
  ind2(k+i)=ind2(i)
  ind1(2*k+i)=ind1(i)
  ind2(2*k+i)=ind2(i)
  ind1(3*k+i)=ind1(i)
  ind2(3*k+i)=ind2(i)
  ind1(4*k+i)=ind1(i)
  ind2(4*k+i)=ind2(i)
  ind1(5*k+i)=ind1(i)
  ind2(5*k+i)=ind2(i)
enddo
 ireal(k+1:2*k)=.false.
 ireal(3*k+1:4*k)=.false.
 ireal(5*k+1:6*k)=.false.
 ipos(1:2*k)=1
 ipos(2*k+1:4*k)=2
 ipos(4*k+1:6*k)=3
endif




sizeind1=size(ind1)
a_f=0
a_f0=0
endif



end subroutine alloc_vector_field_fourier


subroutine set_a_f_a_f0
implicit none
integer i,k
if(skipcomplex) then

k=(size(ind1)+3)/6
do i=1,k
a_f(i)=real(af%f(1,ind1(i),ind2(i),0))
enddo
do i=k+1,2*k-1
a_f(i)=aimag(af%f(1,ind1(i),ind2(i),0))
enddo
do i=2*k,3*k-1
a_f(i)=real(af%f(2,ind1(i),ind2(i),0))
enddo
do i=3*k,4*k-2
a_f(i)=aimag(af%f(2,ind1(i),ind2(i),0))
enddo
do i=4*k-1,5*k-2
a_f(i)=real(af%f(3,ind1(i),ind2(i),0))
enddo
do i=5*k-1,6*k-3
a_f(i)=aimag(af%f(3,ind1(i),ind2(i),0))
enddo

else

k=size(ind1)/6
do i=1,k
a_f(i)=real(af%f(1,ind1(i),ind2(i),0))
enddo
do i=k+1,2*k
a_f(i)=aimag(af%f(1,ind1(i),ind2(i),0))
enddo
do i=2*k+1,3*k
a_f(i)=real(af%f(2,ind1(i),ind2(i),0))
enddo
do i=3*k+1,4*k
a_f(i)=aimag(af%f(2,ind1(i),ind2(i),0))
enddo
do i=4*k+1,5*k
a_f(i)=real(af%f(3,ind1(i),ind2(i),0))
enddo
do i=5*k+1,6*k
a_f(i)=aimag(af%f(3,ind1(i),ind2(i),0))
enddo
endif


a_f0=a_f

end subroutine set_a_f_a_f0


subroutine kill_vector_field_fourier(a)
 implicit none
 type(vector_field_fourier) a
 integer n1,n2,nd
 a%nd=0
 a%muf=0
 a%mu=0
a%n1=0
a%n2=0
a%n3=0
deallocate(a%f)
 
end subroutine kill_vector_field_fourier

subroutine mul_op_h1_h2(h1,h2,hf)
! h1.L h2  
implicit none
 type(vector_field_fourier) h1
 type(spinor_fourier) h2
  type(spinor_fourier) hf
complex(dp), allocatable :: ht(:,:,:,:)
integer k,n1,n2,n3
integer i1,i2,i3,j1,j2,j3,k1,k2,k3
integer i1k,i2k,i3k
complex(dp) c

k=size(h1%f,1)
n1=(size(h1%f,2)-1)/2
n2=(size(h1%f,3)-1)/2
n3=(size(h1%f,4)-1)/2

call alloc_h(ht,h1%f)

 do j1=-n1,n1
 do k1=-n1,n1
     i1=j1+k1
 !if(absolute) then
 !    i1k=abs(i1)+abs(k1)
 !else
     i1k=i1
 !endif
    if(abs(i1k)>n1) cycle
 do j2=-n2,n2
 do k2=-n2,n2
     i2=j2+k2
 !if(absolute) then
 !    i2k=abs(i2)+abs(k2)
 !else
     i2k=i2
 !endif
     if(abs(i2k)>n2) cycle     
 do j3=-n3,n3
 do k3=-n3,n3     
     i3=j3+k3
 !if(absolute) then
 !    i3k=abs(i3)+abs(k3)
 !else
     i3k=i3
 !endif
     if(abs(i3k)>n3) cycle         
ht(1,i1,i2,i3)=h1%f(2,j1,j2,j3)*h2%s(3,k1,k2,k3)-h1%f(3,j1,j2,j3)*h2%s(2,k1,k2,k3)+ht(1,i1,i2,i3)
ht(2,i1,i2,i3)=h1%f(3,j1,j2,j3)*h2%s(1,k1,k2,k3)-h1%f(1,j1,j2,j3)*h2%s(3,k1,k2,k3)+ht(2,i1,i2,i3)
ht(3,i1,i2,i3)=h1%f(1,j1,j2,j3)*h2%s(2,k1,k2,k3)-h1%f(2,j1,j2,j3)*h2%s(1,k1,k2,k3)+ht(3,i1,i2,i3)

 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 



 do k1=-n1,n1
 do k2=-n2,n2
 do k3=-n3,n3     

c=i_*(k1*h1%mu(1)+k2*h1%mu(2)+k3*h1%mu(3))
        
ht(1,k1,k2,k3)=c*h2%s(1,k1,k2,k3)+ht(1,k1,k2,k3)
ht(2,k1,k2,k3)=c*h2%s(2,k1,k2,k3)+ht(2,k1,k2,k3)
ht(3,k1,k2,k3)=c*h2%s(3,k1,k2,k3)+ht(3,k1,k2,k3)

 enddo
 enddo
 enddo

 hf%s=ht
 
deallocate(ht) 
end subroutine mul_op_h1_h2

subroutine trans_mu_s(mu,s)
! h1.L h2  
implicit none
  real(dp) mu(3)
  type(spinor_fourier) s
integer  k1,k2,k3,n1,n2,n3
complex(dp) c
 

 do k1=-n1,n1
 do k2=-n2,n2
 do k3=-n3,n3     

c=exp(i_*(k1*mu(1)+k2*mu(2)+k3*mu(3)))
        
s%s(1,k1,k2,k3)=c*s%s(1,k1,k2,k3) 
s%s(2,k1,k2,k3)=c*s%s(2,k1,k2,k3) 
s%s(3,k1,k2,k3)=c*s%s(3,k1,k2,k3) 

 enddo
 enddo
 enddo

end subroutine trans_mu_s


subroutine trans_mu_matrix(mu,a)
! h1.L h2  
implicit none
  real(dp) mu(3)
  type(matrix_fourier) a
  integer i

 do i=1,3
  call trans_mu_s(mu,a%v(i))
 enddo
 
end subroutine trans_mu_matrix



subroutine exp_h_s(h,s,sf)
implicit none
type(vector_field_fourier), intent(in):: h
type(spinor_fourier), intent(inout):: s,sf
type(spinor_fourier) st,stt
real(dp) fac
integer i

call alloc_spinor_fourier(st,s=s)
call alloc_spinor_fourier(stt,s=s)

st=s
stt=s
 
do i=1,logs_exp

!call mul_op_h1_h2(h%f,stt%s,stt%s)
call mul_op_h1_h2(h,stt,stt)
 
 
  stt%s=stt%s/i
  st%s=st%s + stt%s
 
enddo
 sf%s= st%s
 
call kill_spinor_fourier(st)
call kill_spinor_fourier(stt)

end subroutine exp_h_s

subroutine evaluate_matrix_fourier(phi,a,mc,mr)
implicit none
type(matrix_fourier), intent(in):: a
complex(dp), optional :: mc(3,3)
complex(dp) c,mt(3)
real(dp), optional ::  mr(3,3)
real(dp) phi(3)    
integer i

do i=1,3

call  evaluate_spinor_fourier(phi,a%v(i),mc(1:3,i),mr(1:3,i))

enddo
 
end subroutine evaluate_matrix_fourier

subroutine evaluate_spinor_fourier(phi,sf,mc,mr)
implicit none
type(spinor_fourier), intent(in):: sf
complex(dp), optional :: mc(3)
complex(dp) c,mt(3)
real(dp), optional ::  mr(3)
real(dp) phi(3)    
integer i,j,k,l



mt=0.0_dp
 


do j=-sf%n1,sf%n1
do k=-sf%n2,sf%n2
do l=-sf%n3,sf%n3
do i=1,3
 c=i_*(j*phi(1)+k*phi(2)+l*phi(3))
mt(i)=mt(i)+ sf%s(i,j,k,l)*exp(c) 
enddo
enddo
enddo
enddo
 
if(present(mc)) mc=mt
if(present(mr)) mr=mt
 
 
end subroutine evaluate_spinor_fourier


subroutine equal_spinor_fourier_spinor_fourier(s2,s1)
implicit none
type(spinor_fourier), intent(in):: s1
type(spinor_fourier), intent(inout):: s2

if(.not.associated(s2%s)) call alloc_spinor_fourier(s2,s=s1)
s2%n1=s1%n1
s2%n2=s1%n2
s2%n3=s1%n3
s2%s=s1%s

end subroutine equal_spinor_fourier_spinor_fourier

subroutine equal_matrix_fourier_matrix_fourier(a2,a1)
implicit none
type(matrix_fourier), intent(in):: a1
type(matrix_fourier), intent(inout):: a2
integer i
!if(.not.allocated(a2%v)) call alloc_matrix_fourier(a2,s=a1%v(1))

a2%muf=a1%muf
 
do i=1,3
 a2%v(i)=a1%v(i)
enddo

end subroutine equal_matrix_fourier_matrix_fourier

subroutine alloc_h(h,h1)
implicit none
complex(dp), intent (in) :: h1(:,:,:,:) 
complex(dp), allocatable :: h(:,:,:,:) 
integer k,n1,n2,n3

k=size(h1,1)
n1=(size(h1,2)-1)/2
n2=(size(h1,3)-1)/2
n3=(size(h1,4)-1)/2
allocate(h(k,-n1:n1,-n2:n2,-n3:n3))
h=0.0_dp


end subroutine alloc_h



 subroutine fourier_logs_new(spinmap,ns,U,fix,tune)  !,em,ns)
! part 1 just compute h before factorization
 implicit none
 type(vector_field_fourier) spinmap
 integer i,j,ns(3),k,i1,i2
 real(dp) dphi1,dphi2,x(6),mphi
type(c_damap), optional, intent(inout):: U
real(dp), optional, intent(inout):: fix(6),tune(:)

type(c_damap) c_map, c_id
type(c_normal_form) c_n
 type(probe) xs0
 type(probe_8) xs
type(internal_state) state
type(c_ray) r
type(c_spinor) n
type(spinor) nr


 allocate(spinmap%mr(0:ns(1)-1,0:ns(2)-1,0:ns(3)-1,3,3))
 allocate(spinmap%phis(0:ns(1)-1,0:ns(2)-1,0:ns(3)-1,3))
 allocate(spinmap%x_i(0:ns(1)-1,0:ns(2)-1,0:ns(3)-1,1:6))
 allocate(spinmap%sp(0:ns(1)-1,0:ns(2)-1,0:ns(3)-1))
if(onefunc) then
 if(.not.allocated(yfit)) allocate( yfit(1) )
else
 if(.not.allocated(yfit)) allocate( yfit(ns(1)*ns(2)*ns(3)) )

endif
 yfit=0
write(6,*) sizeind1
allocate(dyfit(size(yfit),sizeind1))
dyfit=0
!pause 76
spinmap%phis=0.0_dp
spinmap%mr=0.0_dp
spinmap%ns=ns
spinmap%f=0.0_dp

 
state=my_eSTATE-spin0
 
if(present(fix)) then
my_fix=fix
else
 my_fix=0.0_dp
 call find_orbit_x(my_ering,my_fix,my_eSTATE,1.e-8_dp,fibre1=start)  
endif

if(present(U))  then
spinmap%muf(1:spinmap%nd) = twopi*tune(1:spinmap%nd)
 call alloc(c_map,c_id); call alloc(c_n); call alloc(n)
c_id=U
 
else
 call init_all(my_estate,1,0)
 call alloc(c_map,c_id); call alloc(c_n); call alloc(n)

c_id=1
xs0=my_fix
xs=xs0+c_id
call propagate(my_ering,xs,my_estate,fibre1=start)  
c_id=xs
call c_normal(c_id,c_n)

spinmap%muf(1:2)=twopi*c_n%tune(1:2)

c_id=c_n%a_t
endif

dphi1=twopi/ns(1)
dphi2=twopi/ns(2)


r=0
do i=0,ns(1)-1
do j=0,ns(2)-1
r%x(1)= sqrt(spinmap%em(1))*cos(i*dphi1)
r%x(2)=-sqrt(spinmap%em(1))*sin(i*dphi1)
r%x(3)= sqrt(spinmap%em(2))*cos(j*dphi2)
r%x(4)=-sqrt(spinmap%em(2))*sin(j*dphi2)
spinmap%phis(i,j,0,1)=i*dphi1
spinmap%phis(i,j,0,2)=j*dphi2
 xs0=0
 r=c_id.o.r
 do k=1,6
   x(k)=r%x(k)+my_fix(k)
 enddo


 spinmap%x_i(i,j,0,1:6)=r%x(1:6)
 spinmap%fix=my_fix

 xs0=x
 call propagate(my_ering,xs0,my_estate,fibre1=start)  
 
 spinmap%mr(i,j,0,1:3,1:3)=xs0
 
call get_logs(xs0,spinmap%sp(i,j,0))
 
enddo
enddo

 
call fourier_trans_logs_new(spinmap) 
 
 call kill(c_map,c_id); call kill(c_n); call kill(n) 

 end  subroutine fourier_logs_new

!type  vector_field_fourier 
!     real(dp) fix(6)    !   closed orbit of map
!     integer :: ns(3)=0   ! integer for Fourier transform
!     real(dp), pointer :: mr(:,:,:,:,:) =>null()   ! spin matrices produced by code (i,j,k,1:3,1:3)
!     real(dp), pointer ::  x_i(:,:,:,:) =>null()   ! starting position in phase  x_i(i,j,k,1:6)=r%x(1:6)
!     real(dp), pointer :: phis(:,:,:,:) =>null()   !%phis(i,j,k,2)=j*dphi2
!     type(spinor), pointer :: sp(:,:,:)            !sp(:i,j,k)   spinor for all the matrices mr  
!     integer n1,n2,n3,nd                           ! # fourier modes -n1:n1, etc... nd=degree of freedom
!     real(dp)  mu(3),em(3)                         ! tune and initial emitances
!     complex(dp),  DIMENSION(:,:,:,:), POINTER :: f  !  vector field expansion f(1:3,-n1:n1,-n2:n2,-n3:n3)
!end  type vector_field_fourier

subroutine fourier_trans_logs_new(spinmap)
! fourier trnasforms all the matrices
implicit none
integer i1,i2,i3,i,j,k,l,ns(3)
real(dp) mphi,dphi1,dphi2,dphi3
type(vector_field_fourier) spinmap

do i=1,3
 ns(i)=size(spinmap%mr,i)
enddo
spinmap%f=0.0_dp

dphi1=1.0_dp/ns(1)
dphi2=1.0_dp/ns(2)
dphi3=1.0_dp/ns(3)

!allocate(a%f(3,-a%n1:a%n1,-a%n2:a%n2,-a%n3:a%n3))

do i1=-spinmap%n1,spinmap%n1
do i2=-spinmap%n2,spinmap%n2
do i3=-spinmap%n3,spinmap%n3

do i=0,ns(1)-1
do j=0,ns(2)-1
do k=0,ns(3)-1

mphi=(i1*i*dphi1+i2*j*dphi2+i3*k*dphi3)*twopi
do l=1,3
 spinmap%f(l,i1,i2,i3)=spinmap%f(l,i1,i2,i3) + exp(-i_*mphi)*spinmap%sp(i,j,k)%x(l)*dphi1*dphi2*dphi3
enddo

enddo
enddo
enddo
enddo
enddo
enddo
end subroutine fourier_trans_logs_new

subroutine get_logs(xs0,nr)
implicit none
type(probe) xs0
type(spinor) nr
type(c_spinor) n
type(c_spinmatrix) s
real(dp) x(3)

call alloc(s)
call alloc(n)


 s=xs0
n=log(s)
nr=n
x(2)=nr%x(2)
if(x(2)<0) then
x(1)=nr%x(1)
x(3)=nr%x(3)

x(1:3)=-x(1:3)*(twopi/sqrt(x(1)**2+x(2)**2+x(3)**2))

nr%x(1)=nr%x(1)+x(1)
nr%x(2)=nr%x(2)+x(2)
nr%x(3)=nr%x(3)+x(3)
 
endif 

call kill(s)
call kill(n)

end subroutine get_logs

subroutine change_vector_field_fourier(af,k,epsi,norm)
implicit none
type(vector_field_fourier) af
type(matrix_fourier) ma
integer i1,i2,i,k,j1,j2,pos,ifit,mf,mft
real(dp) epsi
complex(dp) mac(3,3),zer(3)
real(dp) phi(3),mar(3,3),norm,dmar(3,3),dnorm
 
if(abs(k)>size(ind1) ) then
 write(6,*) " k,size(ind1),size(ind2)",k,size(ind1),size(ind2)
 stop 100 
endif

zer= af%f(1:3,0,0,0) 
 

af%f=0
 
 
do i=1,sizeind1 
 
 if(ireal(i)) then
  af%f(ipos(i),ind1(i),ind2(i),0)  =  a_f(i)  +af%f(ipos(i),ind1(i),ind2(i),0)
  af%f(ipos(i),-ind1(i),-ind2(i),0)=  a_f(i)  +af%f(ipos(i),-ind1(i),-ind2(i),0)
 else
  af%f(ipos(i),ind1(i),ind2(i),0)  =  i_*a_f(i)  +af%f(ipos(i),ind1(i),ind2(i),0)
  af%f(ipos(i),-ind1(i),-ind2(i),0)=  -i_*a_f(i)  +af%f(ipos(i),-ind1(i),-ind2(i),0)
 endif
enddo

 
do i=1,3
 af%f(i,0,0,0)=af%f(i,0,0,0)/2
enddo

 

if(skipzero)  af%f(1:3,0,0,0)=zer

if(epsi>0) then
if(ireal(k)) then
j1=ind1(k);j2=ind2(k);pos=ipos(k);
 
 af%f(pos,j1,j2,0)=af%f(pos,j1,j2,0)+epsi
 af%f(pos,-j1,-j2,0)=af%f(pos,-j1,-j2,0)+epsi
else 
j1=ind1(k);j2=ind2(k);pos=ipos(k);
 
 af%f(pos,j1,j2,0)=af%f(pos,j1,j2,0)+i_*epsi
 af%f(pos,-j1,-j2,0)=af%f(pos,-j1,-j2,0)-i_*epsi
endif
endif

call alloc_matrix_fourier(ma,h=af)

call make_matrix_fourier(af,ma)

norm=0

!write(6,*) size(af%phis,1),size(af%phis,2)
yfit=0
ifit=0
do i1=0,size(af%phis,1)-1
do i2=0,size(af%phis,2)-1

phi=0
phi(1)=af%phis(i1,i2,0,1)
phi(2)=af%phis(i1,i2,0,2)








call evaluate_matrix_fourier(phi,ma,mc=mac,mr=mar)

!write(6,*) ' '
!do i=1,3
! WRITE(6,'(3(1x,g21.14))') mar(i,1:3)
!enddo

dmar=mar-af%mr(i1,i2,0,1:3,1:3)
call norm_spin_mat2(dmar,dnorm)
norm=dnorm+norm

if(.not.onefunc) then
ifit=ifit+1
yfit(ifit)=dnorm
endif

enddo
enddo

norm=norm/size(af%phis,1)/size(af%phis,2)
if(onefunc) yfit(1)=norm
!write(6,*) " norm ",norm

call kill_matrix_fourier(ma)
end subroutine change_vector_field_fourier

subroutine norm_spin_mat2(ma,norm)
implicit none
real(dp) ma(3,3),norm
integer i,j

norm=0.d0

do i=1,3
do j=1,3
norm=norm + ma(i,j)**2
enddo
enddo

norm=norm/9.0_dp

end subroutine norm_spin_mat2

subroutine merit_fourier(a1, yfit1, dyda1, status)
implicit none
  real(dp) eps,del1
  real(dp), intent(in) :: a1(:)
  real(dp), intent(out) :: yfit1(:)
  real(dp), intent(out) :: dyda1(:, :)

  integer status
  integer i,j,k
 
 del1=0.d0

a_f=a1
  eps=0.0_dp
 call change_vector_field_fourier(af,0,eps,del1) 
write(6,*) "norm  ",del1
  yfit1=yfit

  eps=1.d-8
  call set_af_eps(eps,yfit1)
  dyda1=dyfit

!write(6,*) size(yfit)
!write(6,*) size(yfit1)
!write(6,*) size(dyfit,1),size(dyfit,2)
!write(6,*) size(dyda1,1),size(dyda1,2)
!pause 6

!pause 7
!do i=1,size(dyfit,1)
!do j=1,size(dyfit,2)
!write(6,*) i,j,dyfit(i,j)
!enddo
!enddo
!pause 7

end  subroutine merit_fourier

subroutine merit_fourier1(a1, yfit1, dyda1, status)
implicit none
  real(dp) eps,del1
  real(dp), intent(in) :: a1(:)
  real(dp), intent(out) :: yfit1(:)
  real(dp), intent(out) :: dyda1(:, :)

  integer status
  integer i,j,k
 
 del1=0.d0


  yfit1(1)=sin(a1(1)+2*a1(2))-0.1d0

  eps=1.d-8

  dyda1(1,1)= (sin(a1(1)+eps+2*a1(2))-0.1d0-yfit1(1))/eps
  dyda1(1,2)= (sin(a1(1)+2*(a1(2)+eps))-0.1d0-yfit1(1))/eps


end  subroutine merit_fourier1

subroutine set_af_eps(eps,yfit1)
implicit none
integer k
real(dp) yfit1(:)
real(dp) del1,eps

do k=1,sizeind1
 call change_vector_field_fourier(af,k,eps,del1) 
dyfit(:,k)=(yfit-yfit1)/eps
enddo

end subroutine set_af_eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine alloc_logs(a,n1,n2,ns,em,no)
 implicit none
 type(logs) a
 integer i,n1,n2,no,ns
real(dp) em(2)

allocate(a%sp(0:ns-1,0:ns-1))
allocate(a%s(0:ns-1,0:ns-1,3,3))
allocate(a%as(num_iter,3,3))
allocate(a%af(num_iter))
allocate(a%x_i(0:ns-1,0:ns-1,1:6))
allocate(a%phis(0:ns-1,0:ns-1,1:2))

call alloc_explogs(a%h,n1,n2)
call alloc_explogs(a%a,n1,n2)
call alloc_explogs(a%n,n1,n2)

a%as=0.0_dp  
a%x_i=0.0_dp
a%phis=0.0_dp

do i =1,num_iter
 call alloc_explogs(a%af(i),n1,n2)
 a%as(i,1,1)=1.0_dp;a%as(i,2,2)=1.0_dp;a%as(i,3,3)=1.0_dp;
enddo

!(num_iter,3,3)

 a%m=0; a%ms=0;  
  a%ns=ns;a%no=no
 a%em=em
 a%mu=0

 a%s=0.0_dp

 end  subroutine alloc_logs


 subroutine alloc_explogs(a,n1,n2)
implicit none
type(explogs) a
integer n1,n2

   a%n1=n1 ; a%n2=n2;
   allocate(a%h(3,-n1:n1,-n2:n2))
   a%h = 0.0_dp
 end subroutine alloc_explogs

 subroutine kill_explogs(a)
implicit none
type(explogs) a

   deallocate(a%h)
 
 end subroutine kill_explogs

subroutine copy_logs(a,b)
 implicit none
 type(logs) a,b
 integer i,j

call  copy_explogs(a%h,b%h)
 b%ns=a%ns
 b%em=a%em

 do i=0,a%ns-1
 do j=0,a%ns-1
  b%sp(i,j)=a%sp(i,j)
 enddo
 enddo

end subroutine copy_logs


 subroutine copy_explogs(a,b)
 implicit none
 type(explogs) a,b

 b%h=a%h
 b%n1=a%n1
 b%n2=a%n2
 
 end subroutine copy_explogs

 subroutine fourier_logs(spinmap,U,fix,tune)  !,em,ns)
! part 1 just compute h before factorization
 implicit none
 type(logs) spinmap
 integer i,j,ns,k,i1,i2
 real(dp) dphi1,dphi2,x(6),em(2),mphi
type(c_damap), optional, intent(inout):: U
real(dp), optional, intent(inout):: fix(6),tune(:)

type(c_damap) c_map, c_id
type(c_normal_form) c_n
 type(probe) xs0
 type(probe_8) xs
type(internal_state) state
type(c_ray) r
type(c_spinor) n
type(spinor) nr
 
 
ns=spinmap%ns
em=spinmap%em
write(6,*) em
 
state=my_eSTATE-spin0
 
if(present(fix)) then
my_fix=fix
else
 my_fix=0.0_dp
 call find_orbit_x(my_ering,my_fix,my_eSTATE,1.e-8_dp,fibre1=start)  
 
endif

if(present(U))  then
spinmap%mu=twopi*tune(1:2)
 call alloc(c_map,c_id); call alloc(c_n); call alloc(n)
c_id=U
else
 call init_all(my_estate,1,0)
 call alloc(c_map,c_id); call alloc(c_n); call alloc(n)

c_id=1
xs0=my_fix
xs=xs0+c_id
call propagate(my_ering,xs,my_estate,fibre1=start)  
c_id=xs
call c_normal(c_id,c_n)

spinmap%mu=twopi*c_n%tune(1:2)

c_id=c_n%a_t
endif

dphi1=twopi/ns
dphi2=twopi/ns


r=0
do i=0,ns-1
do j=0,ns-1
r%x(1)= sqrt(em(1))*cos(i*dphi1)
r%x(2)=-sqrt(em(1))*sin(i*dphi1)
r%x(3)= sqrt(em(2))*cos(j*dphi2)
r%x(4)=-sqrt(em(2))*sin(j*dphi2)
spinmap%phis(i,j,1)=i*dphi1
spinmap%phis(i,j,2)=j*dphi2
 xs0=0
 r=c_id.o.r
 do k=1,6
   x(k)=r%x(k)+my_fix(k)
 enddo

 spinmap%x_i(i,j,1:6)=r%x(1:6)
 spinmap%fix=my_fix

 xs0=x
 call propagate(my_ering,xs0,my_estate,fibre1=start)  
 
 spinmap%s(i,j,1:3,1:3)=xs0
 
call get_logs(xs0,spinmap%sp(i,j))
 
!spinmap%sp(i,j)=nr
!write(6,*) i,j
!write(6,*) nr%x
enddo
enddo


call fourier_trans_logs(spinmap)

 call kill(c_map,c_id); call kill(c_n); call kill(n) 

 end  subroutine fourier_logs


subroutine fourier_trans_logs(spinmap)
! fourier trnasforms all the matrices
implicit none
integer i1,i2,i,j,k
real(dp) mphi,dphi1,dphi2
 type(logs) spinmap

spinmap%h%h=0.0_dp

dphi1=1.0_dp/spinmap%ns
dphi2=1.0_dp/spinmap%ns

do i1=-spinmap%h%n1,spinmap%h%n1
do i2=-spinmap%h%n2,spinmap%h%n2

do i=0,spinmap%ns-1
do j=0,spinmap%ns-1

mphi=(i1*i*dphi1+i2*j*dphi2)*twopi
do k=1,3
 spinmap%h%h(k,i1,i2)=spinmap%h%h(k,i1,i2) + exp(-i_*mphi)*spinmap%sp(i,j)%x(k)*dphi1*dphi2
enddo

enddo
enddo
enddo
enddo
end subroutine fourier_trans_logs

subroutine evaluate_logs(spinmap,i,j,s)
implicit none
integer i1,i2,i,j,k
real(dp) mphi,dphi1,dphi2
 type(logs) spinmap
type(spinor) s

dphi1=1.0_dp/spinmap%ns
dphi2=1.0_dp/spinmap%ns

s%x=0.0_dp

do i1=-spinmap%h%n1,spinmap%h%n1
do i2=-spinmap%h%n2,spinmap%h%n2



mphi=(i1*i*dphi1+i2*j*dphi2)*twopi
do k=1,3
 s%x(k)=exp(i_*mphi)*spinmap%h%h(k,i1,i2) + s%x(k)
enddo

 
enddo
enddo
end subroutine evaluate_logs


subroutine evaluate_all_logs(spinmap,i,j,spinmapn)
implicit none
integer i,j
 type(logs) spinmap,spinmapn
do i=0,spinmap%ns-1
do j=0,spinmap%ns-1
call evaluate_logs(spinmap,i,j,spinmapn%sp(i,j))
enddo
enddo

end subroutine evaluate_all_logs

subroutine multiply_logs_s(spinmap,s)
implicit none
integer i,j 
 type(logs) spinmap
 real(dp) s(3,3) 
 type(probe) xs0

do i=0,spinmap%ns-1
do j=0,spinmap%ns-1
 spinmap%s(i,j,1:3,1:3)=matmul(spinmap%s(i,j,1:3,1:3),s)

 xs0=spinmap%s(i,j,1:3,1:3)
 call get_logs(xs0,spinmap%sp(i,j))

enddo
enddo

 call fourier_trans_logs(spinmap)

end subroutine multiply_logs_s

subroutine multiply_s_logs(s,spinmap)
implicit none
integer i,j 
 type(logs) spinmap
 real(dp) s(3,3) 
 type(probe) xs0

do i=0,spinmap%ns-1
do j=0,spinmap%ns-1
 spinmap%s(i,j,1:3,1:3)=matmul(s,spinmap%s(i,j,1:3,1:3))

 xs0=spinmap%s(i,j,1:3,1:3)
 call get_logs(xs0,spinmap%sp(i,j))

enddo
enddo

 call fourier_trans_logs(spinmap)

end subroutine multiply_s_logs

!  type  logs    
!     integer n1,n2,ns,no
!     complex(dp),  DIMENSION(:,:,:), POINTER :: h
!     type(spinor), pointer :: sp(:,:)
!     real(dp), pointer :: s(:,:,:,:)  !   spin matrices
!     real(dp) em(2)
!  end  type logs



subroutine evaluate_constant_part(spinmap,s,si)
implicit none
integer i1,i2,i,j,k
 type(logs) spinmap
 real(dp) s(3,3),si(3,3)
 type(c_spinor) h
 type(probe) xs0
call alloc(h)

do k=1,3
 h%v(k)=spinmap%h%h(k,0,0)  
enddo

xs0=exp(h)
s=xs0

xs0=exp(-h)
si=xs0

call kill(h)

end subroutine evaluate_constant_part


subroutine normalise_logs_i(spinmap,m)
implicit none
type(logs), intent(inout) :: spinmap
type(c_spinor) spm
type(c_spinmatrix) sm,sm0,smi
real(dp) theta0
integer i,j,k,m
type(explogs)  a,b
type(probe) xs0
!  type  explogs  
!     integer n1,n2
!     complex(dp),  DIMENSION(:,:,:), POINTER :: h
!  end  type explogs


call alloc(spm)
call alloc(sm)
call alloc(sm0)
call alloc(smi)

call alloc_explogs(a,spinmap%h%n1,spinmap%h%n2)
call alloc_explogs(b,spinmap%h%n1,spinmap%h%n2)

!call exp_logs(a,spinmap%h,spinmap%h,spinmap%mu)

do k=1,3
 spm%v(k)=spinmap%h%h(k,0,0)
enddo

write(6,*)" before 0,0 "
write(6,*) spinmap%h%h(1,0,0) 
write(6,*) spinmap%h%h(2,0,0) 
write(6,*) spinmap%h%h(3,0,0) 

        call c_find_as(spm, sm0)
 !   MODULE PROCEDURE EQUAL_c_spinmatrix_3_by_3
 !   MODULE PROCEDURE EQUAL_3_by_3_probe
 !   MODULE PROCEDURE EQUAL_probe_3_by_3
 !   MODULE PROCEDURE EQUAL_probe_c_spinmatrix


sm0=sm0**(-1)

xs0=sm0
spinmap%as(m,:,:)=xs0

call simil_logs(sm0,spinmap%h)


theta0=spinmap%h%h(2,0,0)

write(6,*)" after 0,0 "
write(6,*) spinmap%h%h(1,0,0) 
write(6,*) spinmap%h%h(2,0,0) 
write(6,*) spinmap%h%h(3,0,0) 
write(6,*) " nu0 ",theta0/twopi
write(6,*) " mu(1) ",spinmap%mu(1)
write(6,*) " mu(2) ",spinmap%mu(2)

call n0_to_nr(spinmap%h,b)

do i=-b%n1,b%n1
do j=-b%n2,b%n2

 if(iabs(i)+iabs(j)/=0) then 
  b%h(2,i,j)= b%h(2,i,j)/i_/(i*spinmap%mu(1)+j*spinmap%mu(2))
  else
    b%h(2,i,j)=0
  endif
  b%h(1,i,j)= b%h(1,i,j)/(theta0 +(i*spinmap%mu(1)+j*spinmap%mu(2)))/i_
  b%h(3,i,j)= b%h(3,i,j)/(-theta0+(i*spinmap%mu(1)+j*spinmap%mu(2)))/i_

enddo
enddo

call nr_to_n0(b,a)

spinmap%af(m)%h=a%h

call exp_logs(a,spinmap%h,spinmap%h,spinmap%mu)

call kill_explogs(a)
call kill_explogs(b)

call kill(spm)
call kill(sm)
call kill(sm0)
call kill(smi)

end subroutine normalise_logs_i

subroutine evaluate_logs_a_all(spinmap,x,sm)
implicit none
type(logs) spinmap
integer i
real(dp) sm(3,3),sd(3,3),x(2)

sm=0.0_dp
do i=1,3
 sm(i,i)=1.0_dp
enddo

if(.true.) then
do i=1,num_iter

 sm=matmul(spinmap%as(i,1:3,1:3),sm)

 call evaluate_logs_a(spinmap,x,i,sd)

 sm=matmul(sd,sm)

enddo
else
do i=num_iter,1,-1

 call evaluate_logs_a(spinmap,x,i,sd)

 sm=matmul(sm,sd)

 sm=matmul(sm,spinmap%as(i,1:3,1:3))


enddo
endif
end subroutine evaluate_logs_a_all

subroutine evaluate_logs_a(spinmap,x,m,sm)
implicit none
integer i1,i2,i,j,k,m
real(dp) mphi,x(2),phi1,phi2
 type(logs) spinmap
type(spinor) s
type(c_spinor) sc
real(dp) sm(3,3)


phi1=x(1)
phi2=x(2)

s%x=0.0_dp

do i1=-spinmap%h%n1,spinmap%h%n1
do i2=-spinmap%h%n2,spinmap%h%n2

mphi=(i1*phi1+i2*phi2) 
do k=1,3
 s%x(k)=exp(i_*mphi)*spinmap%af(m)%h(k,i1,i2) + s%x(k)
enddo
 
enddo
enddo

call alloc(sc)

 sc=s

sm=exp(sc)

call kill(sc)

end subroutine evaluate_logs_a

  subroutine bracket_a_b_c(a,b,co)
    implicit none
    TYPE(explogs), INTENT(INout) :: a,b,co
    TYPE(explogs) c
    integer i1,i2,j1,j2
    call alloc_explogs(c,a%n1,a%n2)

do i1=-a%n1,a%n1
do j1=-a%n2,a%n2
do i2=-b%n1,b%n1
do j2=-b%n2,b%n2

  if(iabs(i1+i2)>a%n1) cycle
  if(iabs(j1+j2)>a%n2) cycle

  c%h(1,i1+i2,j1+j2)=c%h(1,i1+i2,j1+j2)+a%h(2,i1,j1)*b%h(3,i2,j2)-a%h(3,i1,j1)*b%h(2,i2,j2)
  c%h(2,i1+i2,j1+j2)=c%h(2,i1+i2,j1+j2)+a%h(3,i1,j1)*b%h(1,i2,j2)-a%h(1,i1,j1)*b%h(3,i2,j2)
  c%h(3,i1+i2,j1+j2)=c%h(3,i1+i2,j1+j2)+a%h(1,i1,j1)*b%h(2,i2,j2)-a%h(2,i1,j1)*b%h(1,i2,j2)

enddo
enddo
enddo
enddo

co%h=c%h

call kill_explogs(c)

end  subroutine bracket_a_b_c


subroutine exp_logs(a,b,c,mu)
implicit none
TYPE(explogs), INTENT(INout) :: a,b,c
TYPE(explogs) d,t
real(dp) fac,mu(:)
integer i,j

call alloc_explogs(d,a%n1,a%n2)
call alloc_explogs(t,a%n1,a%n2)


t%h=b%h
d%h=b%h
do i=1,logs_exp
 call bracket_a_b_c(a,d,d)
 fac=1.0_dp/i
 d%h=fac*d%h
 t%h=t%h+d%h
enddo

d%h=-a%h
     do i=-a%n1,a%n1
     do j=-a%n2,a%n2
         d%h(1:3,i,j)=i_*(i*mu(1)+j*mu(2))*d%h(1:3,i,j)
     enddo
     enddo 


t%h=t%h+d%h    


do i=2,logs_exp
call bracket_a_b_c(a,d,d)
fac=1.0_dp/i
d%h=fac*d%h
t%h=t%h+d%h
enddo

c%h=t%h

call kill_explogs(d)
call kill_explogs(t)

end subroutine exp_logs


subroutine simil_logs(s,h)
implicit none
type(explogs), intent(inout) :: h
type(c_spinor) spm
type(c_spinmatrix), intent(in) :: s 
integer i,j,k

call alloc(spm)
 


do i=-h%n1,h%n1
do j=-h%n2,h%n2
 
do k=1,3
 spm%v(k)=h%h(k,i,j)
enddo
 spm=s*spm
do k=1,3
 h%h(k,i,j)=spm%v(k)
enddo
 
enddo 
enddo

 

call kill(spm)
 
end subroutine simil_logs



  subroutine print_explogs(a,b,file)
    implicit none
    TYPE(explogs), INTENT(INout) :: a,b
    character(*) file
    integer mf,i,j,k
    real(dp) n1(3),n2(3)
    call kanalnummer(mf,file)

     do i=-a%n1,a%n1
     do j=-a%n2,a%n2

!     do i=-2,2
!     do j=-2,2
     if(iabs(i)+iabs(j)>2) cycle
!     write(mf,'(2(i4,1x),6(1x,E15.8))') i,j,a%h(1:3,i,j)
!     write(mf,'(2(5x),6(1x,E15.8))')    b%h(1:3,i,j)
   do k=1,3
    n1(k)=abs(a%h(k,i,j))
    n2(k)=abs(b%h(k,i,j))
   enddo
     write(mf,'(2(i4,1x),3(1x,E15.8))') i,j,n1(1:3)
     write(mf,'(2(5x),3(1x,E15.8))')    n2(1:3)

     enddo
     enddo
    close(mf)

end subroutine print_explogs

  subroutine n0_to_nr(n0,nr)
    implicit none
    TYPE(explogs), INTENT(INout) :: n0,nr
    complex(dp) nt(3)
    integer i,j
     do i=-n0%n1,n0%n1
     do j=-n0%n2,n0%n2
      nt(2)=n0%h(2,i,j)
      nt(1)=n0%h(1,i,j)-i_*n0%h(3,i,j) ! coefficient of  1/2(L_x + i L_z) 
      nt(3)=n0%h(1,i,j)+i_*n0%h(3,i,j) ! coefficient of  1/2(L_x - i L_z)
      nr%h(1:3,i,j)=nt
     enddo
     enddo  
    

  end subroutine n0_to_nr

  subroutine nr_to_n0(nr,n0) 
    implicit none
    TYPE(explogs), INTENT(INout) :: n0,nr
    complex(dp) nt(3)
    integer i,j

  
     do i=-nr%n1,nr%n1
     do j=-nr%n2,nr%n2
     nt(2)=nr%h(2,i,j)
     nt(1)=(nr%h(1,i,j)+nr%h(3,i,j))/2.0_dp    ! coefficient of L_x
     nt(3)=i_*(nr%h(1,i,j)-nr%h(3,i,j))/2.0_dp ! coefficient of L_z
      n0%h(1:3,i,j)=nt
      enddo
     enddo  

  end subroutine nr_to_n0

  subroutine norm_explogs(nr,norm,normphi) 
    implicit none
    TYPE(explogs), INTENT(INout) :: nr
    real(dp) norm,normphi
    integer i,j

     norm=0.0_dp
     normphi=0.0_dp  
     
     do i=-nr%n1,nr%n1
     do j=-nr%n2,nr%n2
     norm=abs(nr%h(1,i,j))+abs(nr%h(2,i,j))+abs(nr%h(3,i,j)) +norm
      if(iabs(i)+iabs(j)/=0) then
      normphi=abs(nr%h(1,i,j))+abs(nr%h(2,i,j))+abs(nr%h(3,i,j)) +normphi
      else
      normphi=abs(nr%h(3,i,j)) +normphi
      endif
      enddo
     enddo  

  end subroutine norm_explogs


!!!!!!!  stuff for demin  !!!!!!!

subroutine remove_energy_patches(r)
implicit none
type(layout), target :: r
type(fibre), pointer :: p
integer i
type(work) w

w=0

p=>r%start
w=p

do i=1,r%n
 p%patch%energy=0
 p=w
 p=>p%next
enddo
 
end subroutine remove_energy_patches

subroutine set_fcc(r,xr,norescale,state,del)
implicit none
type(layout), target :: r
type(fibre), pointer :: p
type(work) werk
real(dp) del0,deld,x(6),xb(6),del,xr(6)
integer nresc,i
logical norescale
type(internal_state) state

if(norescale) then
 call remove_energy_patches(r)
endif

x=xr
p=>r%start
werk=p
del0=werk%energy
deld=werk%p0c

 
 nresc=0
do i=1,r%n-1

xb=x
call track_probe_x(x,state,fibre1=p,fibre2=p%next)
!write(mf,*) p%t2%s(1),x(5)

if(norescale) then
 if(abs(x(5)-xb(5))>del) then
 nresc=nresc+1
  call no_rescale_fcc(r,p,x)  
 
  x=xb
  call track_probe_x(x,state,fibre1=p,fibre2=p%next)
 
  endif
else
 
  call rescale_fcc(r,p,x)  

endif
 

p=>p%next
enddo

p=>r%start
x=xr
 call track_probe_x(x,state,fibre1=p)
write(6,*) x

 werk=r%end
if(norescale) then
 write(6,*) " final energy = ",werk%energy
write(6,*) "delta/p0c = ", (werk%energy-del0)/deld
 x(5)=(werk%energy-del0)/deld
write(6,*) nresc, " # of energy patches "
else
 write(6,*)  " final energy = ",werk%energy+werk%p0c*x(5)
endif

xr=x
end subroutine set_fcc

subroutine no_rescale_fcc(r,f,x)
implicit none
type(layout), target ::r
type(fibre), target ::f
type(fibre), pointer ::p
real(dp) x5,x(6),e
integer i
type(work) w,we
type(internal_state) state

w=f
we=f
!write(6,*) w%energy,w%p0c,w%beta0

e=w%energy+w%p0c*x(5)
!write(6,*) e

call find_energy(w,energy=e)

!write(6,*) w%energy,w%p0c,w%beta0

p=>f%next
do i=1,r%n

p=w

if(associated(p,r%end)) exit
p=>p%next
enddo

!x(2)= we%p0c*x(2)/w%p0c
!x(4)= we%p0c*x(4)/w%p0c
!x(5)=0.d0

f%patch%energy=2

end subroutine no_rescale_fcc

subroutine rescale_fcc(r,f,x)
implicit none
type(layout), target ::r
type(fibre), target ::f
type(fibre), pointer ::p
real(dp) x5,x(6),e
integer i
type(work) w
type(internal_state) state

w=f
e=w%energy+w%p0c*x(5)
call find_energy(w,energy=e)
w%rescale=.true.
w%power=-1
f%next=w


end subroutine rescale_fcc

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

