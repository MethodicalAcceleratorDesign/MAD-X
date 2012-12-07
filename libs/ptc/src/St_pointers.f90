module pointer_lattice
  use madx_ptc_module !etienne, my_state_mad=>my_state, my_ring_mad=>my_ring
  !  use madx_keywords
  USE gauss_dis
  implicit none
  public
  ! stuff for main program
  type(layout),pointer :: my_ering
  type(internal_state),pointer :: my_estate
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
  integer :: i_layout=1, i_layout_t=1
  integer my_lost_position
  private thin
  real(dp) thin
  !  BEAM STUFF
  REAL(DP) SIG(6),ait(6,6)
  type(beam), allocatable:: my_beams(:)
  type(beam), pointer:: my_beam
  INTEGER :: N_BEAM=0,USE_BEAM=1

  TYPE(REAL_8),private :: Y(6)
  TYPE(DAMAP),PRIVATE :: ID
  !  integer nd2,npara
  !  PRIVATE POWER_CAVITY,RADIA
  ! stuff from my fortran
  type(internal_state), target:: etat
  integer,target:: START ,FIN,ORDER,np,start_t
  real(dp),target:: xfix(6) ,DELT0

  INTERFACE SCRIPT
     MODULE PROCEDURE read_ptc_command
  END INTERFACE

   

contains
  subroutine set_lattice_pointers
    implicit none

    print77=.true.
    read77 =.true.

    my_ering => m_u%start
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
    my_fix(1)=.000d0

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
    LOGICAL(LP) STRAIGHT,skip,fixp
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
    ! DYN APERTURE
    REAL(DP) r_in,del_in,DLAM,ang_in,ang_out,dx,targ_tune_alex(2),sexr0
    INTEGER ITE,n_in,POSR
    logical(lp) found_it
    type(fibre),pointer ::p
    ! TRACKING RAYS
    INTEGER IBN,N_name
    REAL(DP) X(6),DT(3),x_ref(6),sc,NLAM,A1,B1,HPHA,B_TESLA,CUR1,CUR2
    REAL(DP)VOLT,PHASE
    INTEGER HARMONIC_NUMBER
    ! changing magnet
    logical(lp) bend_like
    logical exists
    ! remove_patches
    save my_default
    integer :: limit_int(2) =(/4,18/)
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
    logical :: log_estate=.true.
    integer :: mftune=6,nc
    real(dp), allocatable :: tc(:)
    type(integration_node), pointer  :: t

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

    open(unit=mf,file=ptc_fichier)

    do i=1,10000
       read(mf,'(a120)') comT
       COM=COMT
       call context(com)



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
             write(6,*) "Selected Layout",i_layout,"  called  ---> ",my_ering%name
          endif
       case('SELECTTRACKABLELAYOUT','SELECTTRACKABLELATTICE')
          read(mf,*) i_layout_temp
          if(i_layout_temp>m_u%n) then
             write(6,*) " Universe Size ", m_u%n

             write(6,*) " Selected Layout does not exist "

          else

             i_layout=i_layout_temp
             call move_to_layout_i(m_u,my_ering,i_layout)
             write(6,*) "Selected Layout",i_layout,"  called  ---> ",my_ering%name
          endif
       case('SELECTLASTLAYOUT','SELECTLASTLATTICE')
          i_layout_temp=m_u%n
          if(i_layout_temp>m_u%n) then
             write(6,*) " Universe Size ", m_u%n

             write(6,*) " Selected Layout does not exist "

          else

             i_layout=i_layout_temp
             call move_to_layout_i(m_u,my_ering,i_layout)
             write(6,*) "Selected Layout",i_layout,"  called  ---> ",my_ering%name
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
       case('SETORBITRESTORE')
          read(mf,*) restore_mag,restore_magp
       case('SETORBITPHASOR')
          read(mf,*) xsm%ac%x(1:2)
          xsm0%ac%x(1:2)=xsm%ac%x(1:2)
          ptc_node_old=-1
          first_particle=my_false
       case('SETORBITMARKER')   !!! TO CREATE A NODE
         READ(MF,*) i1
         allocate(orbitname(i1))
          do i1=1,size(orbitname)
           READ(MF,*) name
           call context(name)
           orbitname(i1)=name
          enddo
       case('SETORBITPHASORFREQUENCY')
          read(mf,*) xsm%ac%om
          xsm0%ac%om=xsm%ac%om
          ptc_node_old=-1
          first_particle=my_false
       case('SETORBITPHASORTIME','ORBITTIME')
          read(mf,*) xsm%ac%t
          xsm0%ac%t=xsm%ac%t
          x_orbit_sync=0.0_dp
          x_orbit_sync(6)=xsm%ac%t
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
       
          read(mf,*) xsm%ac%t
!          write(6,*) " Using ",unit_time," seconds"
          xsm%ac%t=xsm%ac%t*clight  !*unit_time
          xsm0%ac%t=xsm%ac%t
          x_orbit_sync=0.0_dp
          x_orbit_sync(6)=xsm%ac%t
          ptc_node_old=-1
          first_particle=my_false
       case('INITIALTIMEINMYUNITS','TIMEINMYUNITS','READFROMFILEINITIALTIME')
         INQUIRE (FILE = INITIAL_setting, EXIST = exists)
          if(exists) then
             write(6,*) "file ", INITIAL_setting(1:len_trim(FINAL_setting)), &
                  " exists, interrupt execution if you do not want to overwrite!"
            call kanalnummer(i1,INITIAL_setting)
   !         read(i1,*) xsm%ac%t,unit_time,n_used_patch,include_patch
            read(i1,*) xsm%ac%t,n_used_patch,include_patch
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
         read(mf,*) xsm%ac%t,n_used_patch,include_patch
!             read(mf,*) xsm%ac%t,unit_time,n_used_patch,include_patch
          endif
 !         write(6,*) " Using ",unit_time," seconds"
          xsm%ac%t=xsm%ac%t*clight 
          xsm0%ac%t=xsm%ac%t
          x_orbit_sync=0.0_dp
          x_orbit_sync(6)=xsm%ac%t
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
       case('ALMOSTEXACT')
          almost_exact=.true.
          write(6,*) " Almost Exact = ",almost_exact
       case('TRUELYEXACT','EXACT')
          almost_exact=.false.
          write(6,*) " Almost Exact = ",almost_exact
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
          READ(MF,*) limit_int
          WRITE(6,*) "limit_int =",limit_int
       case('LMAX')
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
          !       case('KIND7WITHMETHOD1')
          !          CALL PUT_method1_in_kind7(my_ering,1000)
       case('THINLENS=1')
          call THIN_LENS_restart(my_ering)
       case('MANUALTHINLENS')
          THIN=-1
          CALL THIN_LENS_resplit(my_ering,THIN,lim=limit_int)
       case('THINLENS')
          READ(MF,*) THIN
          xbend=-1.0_dp
          if(thin<0) then
             READ(MF,*) sexr0,xbend
             if(xbend<0) then
                xbend=-xbend
                radiation_bend_split=MY_true
             endif
             thin=-thin
          endif
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(my_ering,THIN,lim=limit_int,lmax0=lmax,sexr=sexr0,xbend=xbend)
          radiation_bend_split=MY_false
       case('EVENTHINLENS')
          READ(MF,*) THIN
          xbend=-1.0_dp
          if(thin<0) then
             READ(MF,*) sexr0,xbend
             if(xbend<0) then
                xbend=-xbend
                radiation_bend_split=MY_true
             endif
             thin=-thin
          endif
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(my_ering,THIN,EVEN=my_TRUE,lim=limit_int,lmax0=lmax,sexr=sexr0,xbend=xbend)
       case('ODDTHINLENS')
          READ(MF,*) THIN
          xbend=-1.0_dp
          if(thin<0) then
             READ(MF,*) sexr0,xbend
             if(xbend<0) then
                xbend=-xbend
                radiation_bend_split=MY_true
             endif
             thin=-thin
          endif
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(my_ering,THIN,EVEN=my_FALSE,lim=limit_int,lmax0=lmax,sexr=sexr0,xbend=xbend)
          ! thin layout stuff

       case('ADDSURVEYINFO')

          call ADD_SURVEY_INFO(my_ering)

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


          ! BEAMS STUFF
       case('ALLOCATEBEAMS')
          READ(MF,*) N_BEAM
          ALLOCATE(MY_BEAMS(N_BEAM))
          CALL NULLIFY_BEAMS(MY_BEAMS)
       case('DEALLOCATEBEAMS')
          CALL KILL_BEAMS(MY_BEAMS)
          DEALLOCATE(MY_BEAMS)
       case('CREATEBEAM')
          READ(MF,*) USE_BEAM,NUMBER_OF_PARTICLE, FILENAME
          READ(MF,*) CUT,SIG(1:6)
          CALL CONTEXT(FILENAME)
          IF(FILENAME(1:2)=='NO'.OR.FILENAME(1:2)=='no') THEN
             CALL  create_PANCAKE(MY_BEAMS(USE_BEAM),NUMBER_OF_PARTICLE,CUT,SIG,my_ering%start%t1)
          ELSE
             CALL  kanalnummer(mfr)
             OPEN(UNIT=MFR,FILE=FILENAME)
             READ(MF,*) MY_A_NO
             CALL compute_A_4d(my_ering,my_estate,filename,pos,del,MY_A_NO,MY_A)
             CALL  create_PANCAKE(MY_BEAMS(USE_BEAM),NUMBER_OF_PARTICLE,CUT,SIG,my_ering%start%t1,MY_A)
             CALL KILL(MY_A)
             close(mfr)
          ENDIF
       case('COPYBEAM')
          READ(MF,*) I1,I2
          CALL COPY_BEAM(MY_BEAMS(I1),MY_BEAMS(I2))
       case('PRINTBEAM')
          READ(MF,*) i1,filename
          CALL  kanalnummer(mfr)
          OPEN(UNIT=MFR,FILE=FILENAME)
          call PRINT_beam_raw(MY_BEAMS(I1),MFr)
          CLOSE(MFR)
       case('BEAMSTATISTICS')
          READ(MF,*) USE_BEAM

          CALL Stat_beam_raw(MY_BEAMS(USE_BEAM),4,6)

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

          call compute_twiss(my_ering,my_estate,filename,1,del,1,integrated,name,my_true,my_false)
       case('PTCTWISS','TWISS','PTCTWISSRIPKEN','TWISSRIPKEN')  !
          read(mf,*) filename, NAME, integrated
          read(mf,*) del

          call compute_twiss(my_ering,my_estate,filename,1,del,1,integrated,name,my_false,my_false)

       case('PTCTWISSSASHA','TWISSSASHA','PTCTWISSRIPKENSASHA','TWISSRIPKENSASHA')  !
          read(mf,*) filename, NAME, integrated
          read(mf,*) del

          call compute_twiss(my_ering,my_estate,filename,1,del,1,integrated,name,my_false,my_true)

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
          call lattice_fit_CHROM_gmap1(my_ering,my_estate,EPSF,pol_,NPOL,targ_chrom,np,i1,mf)
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
          call lattice_GET_tune(my_ering,my_estate,mftune)
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

       case('EQUILIBRIUMSIZES')
          READ(MF,*) POS,FILENAME,fileTUNE, NAME
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ering,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif
          CALL radia(my_ering,POS,FILENAME,fileTUNE,my_estate)
       case('NEWEQUILIBRIUMSIZES')
          READ(MF,*) POS,FILENAME, i1,NAME
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ering,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif
          call radia_new(my_ering,POS,i1,FILENAME,my_estate)
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

    close(mf)

  END subroutine read_ptc_command

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


  SUBROUTINE radia(R,loc,FILE1,FILE2,estate,ast,asti,kick,mat,mat0,fixrad)
    implicit none
    TYPE(LAYOUT) R

    REAL(DP) X(6),m,as(6,6),energy,deltap
    TYPE(ENV_8) YS(6)   !
    type(beamenvelope) env
    CHARACTER(*) FILE1,FILE2
    type(normalform) normal
    integer nd2,npara,i,j,js(6),n1,n2
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    TYPE(INTERNAL_STATE) state
    TYPE(INTERNAL_STATE), target :: estate
    integer no,loc,mf1,mf2
    real(dp) av(6,6),e(3),mat0i(6,6)
    real(dp),optional :: ast(6,6),asti(6,6),kick(6),mat(6,6),mat0(6,6),fixrad(6)
    type(fibre), pointer :: p
    no=0
    call kanalnummer(mf1)
    open(mf1,file=FILE1)
    call kanalnummer(mf2)
    open(mf2,file=FILE2)


    if(present(mat0)) then
       state=(estate-nocavity0)+radiation0
       x=0.d0
       CALL FIND_ORBIT_x(R,X,STATE,1.0e-5_dp,fibre1=loc)
       CALL INIT(state,1,0)
       CALL ALLOC(Y)
       CALL ALLOC(NORMAL)
       call alloc(id)  ! ALLOCATE VARIABLES
       id=1
       y=id+x
       CALL TRACK_PROBE_x(r,y,state, fibre1=loc)
       id=y
       mat0=id
       normal=id
       mat0i=id**(-1)   ! mat0i has inverse
       write(6,*) "symplectic tunes "
       write(6,*) normal%tune

       CALL kill(Y)
       CALL kill(NORMAL)
       call kill(id)  ! ALLOCATE VARIABLES

    endif
    state=(estate-nocavity0)+radiation0
    x=0.d0

    CALL FIND_ORBIT_x(R,X,STATE,1.0e-5_dp,fibre1=loc)
    WRITE(6,*) " CLOSED ORBIT AT LOCATION ",loc
    write(6,*) x
    if(present(fixrad)) fixrad=x
    if(track_flag(r,x,loc,state)==0) then
       write(6,*) " stable closed orbit tracked "
    else
       write(6,*) " unstable closed orbit tracked "
       stop 333
    endif

    goto 100
    open(unit=30,file='junk.txt')
    call move_to(r,p,loc)
    do i=loc,loc+r%n
       call track(r,x,i,i+1,state)
       p=>p%next
       write(30,205) i,p%mag%name(1:8),x
    enddo
    close(30)
205 FORMAT(1x,i4,1x,a8,1x,6(1X,D18.11))
100 continue

    call GET_loss(r,energy,deltap)
    write(6,*) x
    write(6,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap

    write(mf1,*) " stable closed orbit tracked "
    write(mf1,"(6(1X,D18.11))") x
    write(mf1,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap

    CALL INIT(state,2,0,BERZ,ND2,NPARA)
    CALL ALLOC(Y)
    CALL ALLOC(NORMAL)
    call alloc(ys)
    call alloc(env)  ! ALLOCATE VARIABLES
    !Y=NPARA
    CALL ALLOC(ID)
    ID=1
    Y=X+ID
    ys=y

    CALL TRACK(R,YS,loc,state)
    call kanalnummer(npara,'junkys.txt')
    write(6,*) loc
    call print(ys,npara)
    close(npara)
    if(.not.check_stable) write(6,*) " unstable tracking envelope "
    env%stochastic=.true.
    env=ys
    if(.not.check_stable) write(6,*) " unstable in normalizing envelope "

    y=ys
    normal=y
    if(.not.check_stable) write(6,*) " unstable in normalizing map "
    as=normal%a_t
    goto 111
    id=y
    open(unit=66,file='crap.txt')
    call print(id,66)
    do i=1,6
       do j=1,6
          m=ys(i)%e(j)
          write(66,*) i,j,m
       enddo
    enddo
    close(66)
111 continue
    !  TYPE beamenvelope
    !     ! radiation normalization
    !     type (damap) transpose    ! Transpose of map which acts on polynomials
    !     type (taylor) bij         !  Represents the stochastic kick at the end of the turn  Env_f=M Env_f M^t + B
    !     TYPE (pbresonance) bijnr   !  Equilibrium beam sizes in resonance basis
    !     type (taylor) sij0  !  equilibrium beam sizes
    !     real(dp) emittance(3),tune(3),damping(3)
    !     logical AUTO,STOCHASTIC
    !     real(dp)  KICK(3)
    !     type (damap) STOCH
    !  END TYPE beamenvelope

    write(6,*) " Chao emittance "
    write(6,*) env%emittance
    write(6,*) " tunes "
    write(6,*) env%tune
    write(6,*) " damping decrements "
    write(6,*) env%damping
    write(mf1,*) " Chao emittance "
    write(mf1,*) env%emittance
    write(mf1,*) " tunes "
    write(mf1,*) env%tune
    write(mf1,*) " damping decrements "
    write(mf1,*) env%damping
    js=0
    do i=1,6
       do j=1,6
          js(j)=1
          m=env%STOCH%v(i).sub.js
          write(mf1,*) m,i,j
          js(j)=0
       enddo
    enddo
    js=0
    do i=1,6
       do j=1,6
          js(j)=1
          m=ys(i)%v.sub.js
          ! write(mf1,*) m,i,j
          js(j)=0
       enddo
    enddo
    if(present(kick)) then
       kick(1)=env%kick(1)
       kick(2)=env%kick(1)
       kick(3)=env%kick(2)
       kick(4)=env%kick(2)
       kick(5)=env%kick(3)
       kick(6)=env%kick(3)
    endif
    write(mf1,*) env%kick
    write(mf1,*) "B matrix"
    call print(env%STOCH,mf1)
    write(mf1,*) "m matrix"

    if(present(mat)) then
       id=ys%v
       mat=id
       mat=matmul(mat,mat0i)
    endif
    call print(ys%v,mf1)
    write(mf1,*) "equilibrium <X_i X_j>"
    call print(env%sij0,mf1)
    write(mf1,*) " Resonance Fluctuation "
    call print(env%bijnr,mf1)
    ! for etienne
    write(mf2,*) env%kick
    call print(env%STOCH,mf2)
    if(present(ast)) ast=env%STOCH
    env%STOCH=env%STOCH**(-1)
    if(present(asti))asti=env%STOCH
    call print(env%STOCH,mf2)
    call print(ys%v,mf2)
    write(mf2,*) " Damping  "
    write(mf2,*) env%damping
    write(mf2,*) " Stochastic Theta "
    call print(env%bij,mf2)

    js=0
    av=0.d0
    e=env%emittance
    write(mf1,*) " emmitances "
    write(mf1,*) e
    av(1,1)=(as(1,1)**2+as(1,2)**2)*e(1)+(as(1,3)**2+as(1,4)**2)*e(2)+ &
         (as(1,5)**2+as(1,6)**2)*e(3)
    av(3,3)=(as(3,1)**2+as(3,2)**2)*e(1)+(as(3,3)**2+as(3,4)**2)*e(2)+ &
         (as(3,5)**2+as(3,6)**2)*e(3)
    av(5,5)=(as(5,1)**2+as(5,2)**2)*e(1)+(as(5,3)**2+as(5,4)**2)*e(2)+ &
         (as(5,5)**2+as(5,6)**2)*e(3)
    av(3,5)=(as(3,1)*as(5,1)+as(3,2)*as(5,2))*e(1)+(as(3,3)*as(5,3)+as(3,4)*as(5,4))*e(2)+ &
         (as(3,5)*as(5,5)+as(3,6)*as(5,6))*e(3)
    n1=1
    n2=3
    av(n1,n2)=(as(n1,1)*as(n2,1)+as(n1,2)*as(n2,2))*e(1)+(as(n1,3)*as(n2,3)+as(n1,4)*as(n2,4))*e(2)+ &
         (as(n1,5)*as(n2,5)+as(n1,6)*as(n2,6))*e(3)
    n1=3
    n2=6
    av(n1,n2)=(as(n1,1)*as(n2,1)+as(n1,2)*as(n2,2))*e(1)+(as(n1,3)*as(n2,3)+as(n1,4)*as(n2,4))*e(2)+ &
         (as(n1,5)*as(n2,5)+as(n1,6)*as(n2,6))*e(3)

    m=env%sij0.sub.'2'

    write(mf1,*) " <X**2> exact and alex "
    write(mf1,*) m
    write(mf1,*) av(1,1)
    m=env%sij0.sub.'002'
    write(mf1,*) " <y**2> exact and alex "
    write(mf1,*) m
    write(mf1,*) av(3,3)
    m=env%sij0.sub.'00002'
    write(mf1,*) " <L**2> exact and alex "
    write(mf1,*) m
    write(mf1,*) av(5,5)
    m=env%sij0.sub.'001010'
    write(mf1,*) " <y delta> exact and alex "
    write(mf1,*) m/2.d0
    write(mf1,*) av(3,5)
    n1=1
    n2=3
    m=env%sij0.sub.'101000'
    write(mf1,*) " <x y> exact and alex "
    write(mf1,*) m/2.d0
    write(mf1,*) av(n1,n2)
    n1=3
    n2=6
    m=env%sij0.sub.'001001'
    write(mf1,*) " <y L> exact and alex "
    write(mf1,*) m/2.d0
    write(mf1,*) av(n1,n2)



    CALL KILL(Y)
    CALL KILL(Ys)
    CALL KILL(NORMAL)
    CALL KILL(env)

    ! compute map with radiation minus the cavity!
    ! cavity must be at the end and only one cavity

    if(no>0) then
       CALL INIT(STATE,no,0,BERZ,ND2,NPARA)
       CALL ALLOC(Y)  ! ALLOCATE VARIABLES
       write(17,*) x(1),x(2),x(3)
       write(17,*) x(4),x(5),x(6)
       Y=NPARA
       Y=X

       CALL TRACK(R,Y,1,r%n,STATE)

       call print(y,17)
       call kill(y)
    endif
    close(mf1)
    close(mf2)


    if(present(mat)) then
       call kanalnummer(mf1)
       open(mf1,file='barber_stochastic.txt')
       if(present(mat)) then
          do i=1,6
             do j=1,6
                write(mf1,*) i,j,mat(i,j)," mat"
             enddo
          enddo
       endif
       if(present(ast)) then
          do i=1,6
             do j=1,6
                write(mf1,*) i,j,ast(i,j)," ast"
             enddo
          enddo
       endif
       if(present(asti)) then
          do i=1,6
             do j=1,6
                write(mf1,*) i,j,asti(i,j)," asti"
             enddo
          enddo
       endif
       if(present(kick)) then
          do i=1,6
             write(mf1,*) i,kick(i)," kick"
          enddo
       endif

       close(mf1)
    endif

  end subroutine radia

  SUBROUTINE radia_new(R,loc,i1,FILE1,estate)
    implicit none
    TYPE(LAYOUT) R

    REAL(DP) X(6),m,energy,deltap
    CHARACTER(*) FILE1
    type(normal_spin) normal
    integer  i,j ,i1
    TYPE(damapspin) ID
    TYPE(INTERNAL_STATE) state
    TYPE(INTERNAL_STATE), target :: estate
    integer loc,mf1
    type(fibre), pointer :: p
    type(probe) xs0
    type(probe_8) xs

    call kanalnummer(mf1)
    open(mf1,file=FILE1)




    state=(estate-nocavity0)+radiation0
    x=0.d0

    CALL FIND_ORBIT_x(R,X,STATE,1.0e-5_dp,fibre1=loc)
    WRITE(6,*) " CLOSED ORBIT AT LOCATION ",loc
    write(6,*) x




    call GET_loss(r,energy,deltap)

    write(6,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap

    write(mf1,*) " stable closed orbit tracked "
    write(mf1,"(6(1X,D18.11))") x
    write(mf1,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap

    CALL INIT(state,1,0)
    CALL ALLOC(NORMAL)
    CALL ALLOC(ID)
    call alloc(xs)

    if(i1==1) normal%stochastic=my_true
    xs0=x
    ID=1
    xs=XS0+ID

    state=state+envelope0

    CALL TRACK_PROBE(r,xs,state, fibre1=loc)
    write(mf1,*) " Full Map "
     id=xs    
     call print(id,mf1)
    write(mf1,*) " End of Full Map "

    normal=id
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
    write(mf1,*)" Equilibrium Phasor Sizes "
    do i=1,6
       do j=1,6
          write(mf1,*) i,j,normal%s_ijr(i,j)
       enddo
    enddo

    if(normal%stochastic) then
       write(mf1,*) "Stochastic kicks"
       write(mf1,*) normal%kick

       write(mf1,*)" Stochastic Transformation "
       do i=1,6
          do j=1,6
             write(mf1,*) i,j,normal%STOCH(i,j)
          enddo
       enddo

       write(mf1,*)" Inverse Stochastic Transformation "
       do i=1,6
          do j=1,6
             write(mf1,*) i,j,normal%STOCH_inv(i,j)
          enddo
       enddo

    endif

    close(mf1)

    CALL KILL(NORMAL)
    CALL KILL(ID)
    CALL KILL(xs)

  end subroutine radia_new

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

