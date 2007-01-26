module pointer_lattice
  use madx_ptc_module, my_state_mad=>my_state, my_ring_mad=>my_ring
  !  use madx_keywords
  USE gauss_dis
  implicit none
  public
  ! stuff for main program
  type(layout),pointer :: my_ring
  type(internal_state),pointer :: my_state
  integer ,pointer :: my_start, MY_ORDER, MY_NP,MY_END,my_start_t
  real(dp), pointer :: my_fix(:),MY_DELTA
  real(dp), pointer ::my_target(:)
  type(pol_block), pointer :: my_pol_block(:)
  integer, pointer :: n_pol_block
  real(sp) my_scale_planar
  ! end of stuff for main program
  integer :: lat(0:6,0:6,3)=1,mres(3)
  character(25) lat_name(81)
  real(dp) r_ap(2),x_ap,y_ap   ! default aperture
  integer :: kind_ap=3,file_ap=0                      ! Single aperture position and kind of aperture + file unit
  character(nlp) name_ap
  character(255) :: filename_ap = "Tracking.txt"
  integer last_npara
  integer :: i_layout=1
  integer my_lost_position
  real(dp) thin
  !  BEAM STUFF
  REAL(DP) SIG(6)
  type(beam), allocatable:: my_beams(:)
  type(beam), pointer:: my_beam
  INTEGER :: N_BEAM=0,USE_BEAM=1

  TYPE(REAL_8),private :: Y(6)
  TYPE(DAMAP),PRIVATE :: ID
  integer nd2,npara
  !  PRIVATE POWER_CAVITY,RADIA
contains

  subroutine read_ptc_command(ptc_fichier)
    use madx_ptc_module , piotr_state=>my_state,piotr_my_ring=>my_ring
    implicit none
    CHARACTER*(120) com,COMT,filename,name_root,title,name_root_res,filetune,FILESMEAR
    character*(4) suffix,SUFFIX_res
    character(*) ptc_fichier
    integer i,ii,mf,i_layout_temp,LIM(2),IB,NO
    !  FITTING FAMILIES
    INTEGER NPOL,J,NMUL,K,ICN,N,np,MRESO(3)
    type(pol_block), ALLOCATABLE :: pol_(:)
    type(pol_block) :: pb
    CHARACTER*(NLP) NAME,flag
    real(dp) targ_tune(2),targ_chrom(2),EPSF
    real(dp) targ_RES(4)
    !  END FITTING FAMILIES
    real(dp),allocatable :: beta(:,:,:)
    REAL(DP) DBETA,tune(3),tunenew(2),CHROM(2),DEL
    ! fitting and scanning tunes
    real(dp) tune_ini(2),tune_fin(2),dtu(2)
    integer nstep(2),i1,i2,I3,neq,i4,i5,I6,it
    LOGICAL(LP) STRAIGHT,skip
    ! end
    ! TRACK 4D NORMALIZED
    INTEGER POS,NTURN,resmax
    real(dp) EMIT(6),APER(2),emit0(2)
    integer nscan,mfr,ITMAX,MRES(4)
    real(dp), allocatable :: resu(:,:)
    ! END
    ! RANDOM MULTIPOLE
    INTEGER iseed,nMULT,addi
    REAL(DP) CUT,cn,cns
    LOGICAL(LP) integrated
    ! END
    ! LOCAL BEAM STUFF
    INTEGER NUMBER_OF_PARTICLE,printmod
    TYPE(DAMAP) MY_A
    INTEGER MY_A_NO,MY_A_ND
    ! APERTURE
    REAL(DP)  APER_R,APER_X,APER_Y
    INTEGER KINDAPER
    TYPE(integration_node), POINTER :: TL
    type(internal_state),target :: my_default
    ! DYN APERTURE
    REAL(DP) r_in,del_in,DLAM,ang_in,ang_out
    INTEGER ITE,n_in,POSR
    logical(lp) found_it
    type(fibre),pointer ::p,p1
    ! TRACKING RAYS
    INTEGER NRAYS
    INTEGER IBN,N_name
    REAL(DP) X(6),DT(3),x_ref(6),sc
    REAL(DP)VOLT,PHASE
    INTEGER HARMONIC_NUMBER
    ! changing magnet
    logical(lp) bend_like
    integer              :: apertflag
    character(200)       :: whymsg
    integer              :: why(9)
    ! remove_patches
    logical(lp) do_not_remove,put_geo_patch
    real(dp) :: lmax=1.e38_dp
    save my_default
    integer :: limit_int(2) =(/4,18/)
    LOGICAL :: track_k,CLOSED_ORBIT=MY_FALSE
    REAL(DP) CLOSED(6),XP(6),te(6)
    ! automatic track
    integer lim_t(6,2),IEXT
    real(dp) dlim_t(6)
    type(internal_state),pointer :: my_old_state

    if(associated(my_state)) my_old_state=>my_state
    my_default=default
    my_state=>my_default
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
       if(com(1:1)=='!'.and.com(2:2)/='!') THEN
          cycle
       ENDIF
       if(com(1:5)=='PAUSE') THEN
          com=com(1:5)
       ENDIF
       if(.not.skip) then
          if(com(1:2)=="!!") then
             skip=.true.
             cycle
          endif
       endif
       if(skip) then !1
          if(com(1:2)=="!!") then
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
             call move_to_layout_i(m_u,my_ring,i_layout)
             write(6,*) "Selected Layout",i_layout,"  called  ---> ",my_ring%name
          endif
       case('NAMELAYOUT')
          read(mf,*) NAME
       case('KILLLASTLAYOUT')
          read(mf,*) NAME
          call kill_last_layout(m_u)
          MY_RING%NAME=NAME
       case('RESTOREMYSTATE')
          my_default=my_OLD_state
          my_state=>my_OLD_state
       case('DEFAULT')
          MY_STATE=DEFAULT0
       case('+NOCAVITY')
          MY_STATE=MY_STATE+NOCAVITY0
       case('-NOCAVITY')
          MY_STATE=MY_STATE-NOCAVITY0
       case('+FRINGE')
          MY_STATE=MY_STATE+FRINGE0
       case('-FRINGE')
          MY_STATE=MY_STATE-FRINGE0
       case('+TIME')
          MY_STATE=MY_STATE+TIME0
       case('-TIME')
          MY_STATE=MY_STATE-TIME0
       case('+TOTALPATH')
          MY_STATE=MY_STATE+TOTALPATH0
       case('-TOTALPATH')
          MY_STATE=MY_STATE-TOTALPATH0
       case('+ONLY_4D')
          MY_STATE=MY_STATE+ONLY_4D0
       case('-ONLY_4D')
          MY_STATE=MY_STATE-ONLY_4D0
       case('+DELTA')
          MY_STATE=MY_STATE+DELTA0
       case('-DELTA')
          MY_STATE=MY_STATE-DELTA0
       case('+RADIATION')
          MY_STATE=MY_STATE+RADIATION0
       case('-RADIATION')
          MY_STATE=MY_STATE-RADIATION0


       case('ALLOCATEBETA')
          ALLOCATE(BETA(2,2,MY_RING%N))

       case('DEALLOCATEBETA')
          DEALLOCATE(BETA)
       case('FILLBETA')
          READ(MF,*) IB,pos
          CALL FILL_BETA(MY_RING,MY_STATE,pos,BETA,IB,DBETA,tune,tunenew)
       case('FITBENDS')
          CALL fit_all_bends(MY_RING,MY_STATE)
       case('LIMITFORCUTTING')
          READ(MF,*) limit_int
          WRITE(6,*) "limit_int =",limit_int
       case('LMAX')
          READ(MF,*) LMAX
          WRITE(6,*) "LMAX FOR SPACE CHARGE =",LMAX
       case('FUZZYLMAX')
          READ(MF,*) FUZZY_SPLIT
          WRITE(6,*) "FUZZY LMAX FOR SPACE CHARGE =",LMAX,LMAX*FUZZY_SPLIT
       case('KIND7WITHMETHOD1')
          CALL PUT_method1_in_kind7(MY_ring,1000)
       case('THINLENS=1')
          call THIN_LENS_restart(MY_RING)
       case('THINLENS')
          READ(MF,*) THIN
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(MY_RING,THIN,lim=limit_int,lmax=lmax)
       case('EVENTHINLENS')
          READ(MF,*) THIN
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(MY_RING,THIN,EVEN=my_TRUE,lim=limit_int,lmax=lmax)
       case('ODDTHINLENS')
          READ(MF,*) THIN
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(MY_RING,THIN,EVEN=my_FALSE,lim=limit_int,lmax=lmax)
          ! thin layout stuff
       case('MAKE_THIN_LAYOUT','MAKELAYOUT','MAKE_NODE_LAYOUT')

          if(.not.associated(MY_RING%t)) CALL MAKE_node_LAYOUT(MY_RING)
       case('SURVEY_THIN_LAYOUT','SURVEYLAYOUT','SURVEY_NODE_LAYOUT')

          IF(associated(MY_RING%t)) THEN
             CALL fill_survey_data_in_NODE_LAYOUT(MY_RING)
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
             CALL  create_PANCAKE(MY_BEAMS(USE_BEAM),NUMBER_OF_PARTICLE,CUT,SIG,my_ring%start%t1)
          ELSE
             CALL  kanalnummer(mfr)
             OPEN(UNIT=MFR,FILE=FILENAME)
             READ(MF,*) MY_A_NO
             CALL compute_A_4d(my_ring,my_state,filename,pos,del,MY_A_NO,MY_A)
             CALL  create_PANCAKE(MY_BEAMS(USE_BEAM),NUMBER_OF_PARTICLE,CUT,SIG,my_ring%start%t1,MY_A)
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
       case('TRACKBEAM')
          READ(MF,*) USE_BEAM,POS, NTURN,printmod,filename
          IF(FILENAME(1:2)=='NO'.OR.FILENAME(1:2)=='no') THEN
             mfr=6
          else
             CALL  kanalnummer(mfr)
             OPEN(UNIT=MFR,FILE=FILENAME)
          endif
          DO I1=1,NTURN
             call track( my_ring,MY_BEAMS(USE_BEAM),MY_STATE,POS)
             if(mod(i1,printmod)==0) then
                write(6,*) " turn ",i1
                CALL Stat_beam_raw(MY_BEAMS(USE_BEAM),4,mfr)
             endif
          ENDDO
          if(mfr/=6) close(mfr)

          ! END BEAMS STUFF



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
       case('TURNOFFONEAPERTURE')
          READ(MF,*) POS
          CALL TURN_OFF_ONE_aperture(MY_RING,pos)
       case('SETONEAPERTURE')
          read(MF,*) pos
          read(MF,*) kindaper, APER_R,APER_X,APER_Y
          CALL assign_one_aperture(MY_RING,pos,kindaper,APER_R,APER_X,APER_Y)
          ! end of layout stuff
          ! random stuff
       case('TURNONBEAMBEAMKICK','TURNONBEAMBEAM')
          READ(MF,*) USE_BEAM
          MY_BEAMS(USE_BEAM)%BEAM_BEAM=my_true

       case('GAUSSIANSEED')
          READ(MF,*) I1
          CALL gaussian_seed(i1)
       case('RANDOMMULTIPOLE')
          read(mf,*) name
          read(mf,*) n,cns,cn
          read(mf,*) addi,integrated
          read(mf,*) cut
          write(6,*) " Distribution cut at ",cut," sigmas"
          call lattice_random_error(my_ring,name,cut,n,addi,integrated,cn,cns)
       case('MISALIGNEVERYTHING')
          read(mf,*) SIG(1:6),cut
          CALL MESS_UP_ALIGNMENT(my_ring,SIG,cut)
          ! end of random stuff
       case('SETFAMILIES')
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
                READ(MF,*) N,ICN
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
       case('DEALLOCATEFAMILIES')
          deallocate(POL_)

       case('FITTUNE')
          read(mf,*) epsf
          read(mf,*) targ_tune
          if(targ_tune(1)<=zero) targ_tune=tune(1:2)
          call lattice_fit_TUNE_gmap(my_ring,my_state,epsf,pol_,NPOL,targ_tune,NP)
       case('SCANTUNE')
          STRAIGHT=.FALSE.
          read(mf,*) epsf
          read(mf,*) nstep
          read(mf,*) tune_ini,tune_fin
          read(mf,*) name_root,SUFFIX
          dtu=zero
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
                call lattice_fit_TUNE_gmap(my_ring,my_state,epsf,pol_,NPOL,targ_tune,NP)
                write(title,*) 'tunes =',targ_tune
                I3=I3+1; call create_name(name_root,i3,suffix,filename);
                write(6,*)" printing in "
                write(6,*)filename
                call print_bn_an(my_ring,2,title,filename)
                write(6,*)" PRINTED "
             enddo
          enddo
       case('FITCHROMATICITY')
          read(mf,*) epsf
          read(mf,*) targ_chrom
          call lattice_fit_CHROM_gmap(my_ring,my_state,EPSF,pol_,NPOL,targ_chrom,NP)
       case('FITTUNECHROMATICITY')
          read(mf,*) epsf
          read(mf,*) targ_RES
          call lattice_fit_CHROM_gmap(my_ring,my_state,EPSF,pol_,NPOL,targ_RES,NP)
       case('GETCHROMATICITY')
          call lattice_GET_CHROM(my_ring,my_state,CHROM)
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
          CALL compute_map_4d(my_ring,my_state,filename,pos,del,no)
       case('PRINTAINRESONANCE')
          READ(MF,*) MRES(1:4)
          READ(MF,*) NO,EMIT0
          READ(MF,*) FILENAME
          CALL lattice_PRINT_RES_FROM_A(my_ring,my_state,NO,EMIT0,MRES,FILENAME)
       case('PRINTSTATE')
          CALL PRINT(MY_STATE,6)
       case('PRINTTWICE')
          print77=.false.
          read77=.false.
       case('PRINTBNAN','PRINTANBN','PRINTBN','PRINTAN')
          read(mf,*) title
          read(mf,*) nmul,filename
          call print_bn_an(my_ring,nmul,title,filename)
       case('READBNAN','READANBN','READBN','READAN')
          read(mf,*) filename
          call READ_bn_an(my_ring,filename)
       case('RETURN','EXIT','QUIT','END','STOP')
          goto 100
       case('PTCEND','ENDPTC','APOCALYSPE')
          CALL PTC_END()
          goto 100
       case('TRACK4DNORMALIZED')
          emit=zero
          read(mf,*) IB
          read(mf,*) POS,NTURN,ITMAX,resmax
          read(mf,*) EMIT(1:2),APER(1:2),emit(3),emit(6)
          read(mf,*) filename,filetune,FILESMEAR
          ! emit(3)=1.d38
          CALL track_aperture(my_ring,my_state,beta,dbeta,tune,ib,ITMAX,emit,aper,pos,nturn,FILENAME,filetune,FILESMEAR,resmax)
       case('SCANTRACK4DNORMALIZED')
          emit=zero

          read(mf,*) POS,NTURN,ITMAX,resmax
          read(mf,*) EMIT0(1:2),APER,emit(3),emit(6)
          read(mf,*) nscan,name_root,SUFFIX
          read(mf,*) name_root_res,SUFFIX_res
          ib=1
          allocate(resu(nscan,4))
          do i1=1,nscan
             write(6,*) " CASE NUMBER ",i1
             call create_name(name_root,i1,suffix,filename);
             call READ_bn_an(my_ring,filename)
             call create_name(name_root_res,i1,SUFFIX_res,filename);
             COMT=name_root_res(1:len_trim(name_root_res))//"_resonance"
             call create_name(COMT,i1,SUFFIX_res,filetune);
             COMT=name_root_res(1:len_trim(name_root_res))//"_SMEAR"
             call create_name(COMT,i1,SUFFIX_res,FILESMEAR);
             if(i1/=1) ib=2
             emit(1:2)=EMIT0
             CALL track_aperture(my_ring,my_state,beta,dbeta,tune,ib,ITMAX,emit,aper,pos,nturn,FILENAME,filetune,FILESMEAR,resmax)
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
       case('SEARCHAPERTUREALONGLINE')

          READ(MF,*) r_in,n_in,ang_in,ang_out,del_in,DLAM
          READ(MF,*) POS,NTURN,ITE,FILENAME,name
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ring,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif


          call kanalnummer(mfr)
          open(unit=mfr,file=filename)
          CALL dyn_aper(MY_RING,r_in,n_in,ang_in,ang_out,del_in,dlam,pos,nturn,ite,my_state,MFR)
          close(mfr)

       case('KNOB')

          READ(MF,*) POS, IBN

          if (pos > MY_RING%n) stop;

          p=>MY_RING%start
          do ii=1,pos
             p=>p%next
          enddo


          write(6,*) "El name ", p%mag%name


          CALL INIT(default,3,1,BERZ,ND2,NPARA)

          print*, "Npara is ", c_%NPARA

          pb = 0
          pb%name = p%mag%name
          write(6,*) "IBN ", IBN
          pb%ibn(ibn) = 1

          my_ring = pb

          CALL ALLOC(ID)
          CALL ALLOC(Y)
          x(:)=0
          ID=1
          Y=X+ID

          p=>MY_RING%start
          do ii=1,MY_RING%n
             write(6,*) "##########################################"
             write(6,'(i4, 1x,a, f10.6)') ii,p%mag%name
             write(6,'(a, f9.6, a)') "Ref Momentum ",p%mag%p%p0c," GeV/c"

             call track(my_ring,y,ii,ii+1,default)
             call daprint(y(1),6)
             p=>p%next

          enddo

       case('PREISNERAUTOTRACK')
          READ(MF,*) POS, NTURN, FILENAME,filesmear, name
          READ(MF,*) CLOSED_ORBIT
          read(mf,*) title
          do i1=1,6
             read(mf,*) LIM_T(i1,1:2),dlim_t(i1)
          enddo

          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ring,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          else
             call move_to(my_ring,p1,POS)
          endif

          call kanalnummer(mfr)
          if(title(1:1)=='A') then
             open(unit=mfr,file=filename,position='APPEND',status='OLD')
          ELSE
             open(unit=mfr,file=filename)
          endif

          call kanalnummer(ib)

          if(title(1:1)=='A') then
             open(unit=ib,file=filesmear,position='APPEND',status='OLD')
          ELSE
             open(unit=ib,file=filesmear)
          endif

          IEXT=0
          CLOSED=ZERO
          IF((.NOT.MY_STATE%NOCAVITY).AND.CLOSED_ORBIT) THEN
             IEXT=FIND_ORBIT_FLAG(MY_RING,CLOSED,POS,MY_STATE,c_1d_5)
          ENDIF

          IF(IEXT==0) THEN
             DO I5=LIM_T(5,1),LIM_T(5,2)

                IF(MY_STATE%NOCAVITY.AND.CLOSED_ORBIT) THEN
                   CLOSED(5) = I5*dlim_t(5)
                   IEXT=FIND_ORBIT_FLAG(MY_RING,CLOSED,POS,MY_STATE,c_1d_5)
                ENDIF
                IF(MY_STATE%NOCAVITY.AND.CLOSED_ORBIT) THEN
                   WRITE(6,*) " ENERGY X(5) = ",CLOSED(5)
                ELSE
                   WRITE(6,*) " ENERGY X(5) = ",I5*dlim_t(5)+CLOSED(5)
                ENDIF
                IF(IEXT/=0) THEN
                   WRITE(6,*) " 4D CLOSED ORBIT NOT FOUND "
                   CYCLE
                ENDIF
                DO I6=LIM_T(6,1),LIM_T(6,2)
                   DO I4=LIM_T(4,1),LIM_T(4,2)
                      DO I3=LIM_T(3,1),LIM_T(3,2)
                         DO I2=LIM_T(2,1),LIM_T(2,2)
                            DO I1=LIM_T(1,1),LIM_T(1,2)
                               XP(1)=I1*dlim_t(1)+CLOSED(1)
                               XP(2)=I2*dlim_t(2)+CLOSED(2)
                               XP(3)=I3*dlim_t(3)+CLOSED(3)
                               XP(4)=I4*dlim_t(4)+CLOSED(4)
                               XP(6)=I6*dlim_t(6)+CLOSED(6)
                               IF(MY_STATE%NOCAVITY.AND.CLOSED_ORBIT) THEN
                                  CLOSED(5) = I5*dlim_t(5)
                                  CALL FIND_ORBIT(MY_RING,CLOSED,POS,MY_STATE,c_1d_5)
                                  XP(5)=CLOSED(5)
                               ELSE
                                  XP(5)=I5*dlim_t(5)+CLOSED(5)
                               ENDIF
                               track_k=.true.
                               te(1)=I1*dlim_t(1)
                               te(2)=I2*dlim_t(2)
                               te(3)=I3*dlim_t(3)
                               te(4)=I4*dlim_t(4)
                               te(5)=I5*dlim_t(5)
                               te(6)=I6*dlim_t(6)
                               DO It=1,NTURN ! it
                                  p=>p1
                                  do ii=1,MY_RING%n  !ii

                                     call track(my_ring,XP,ii,ii+1,MY_STATE)


                                     call produce_aperture_flag(apertflag)
                                     if (apertflag/=0) then
                                        !     print *, 'Particle out of aperture!'

                                        call ANALYSE_APERTURE_FLAG(apertflag,why)
                                        !    Write(6,*) "PREISNER: paticle : ",i1, " lost"
                                        !    Write(6,301) te
                                        !     Write(6,*) "PREISNER: APERTURE error for element: ",ii," name: ",p%MAG%name
                                        !    Write(6,*) "Message: ",c_%message
                                        !                      write(whymsg,*) 'APERTURE error: ',why
                                        track_k=.false.
                                        goto 3001
                                        exit; !goes to the ne
                                     endif

                                     p=>p%next

                                  enddo ! ii


                               ENDDO  ! it
3001                           continue
                               if (track_k) then
                                  !                       Write(6,*) "PREISNER: particle  stable"
                                  !                       Write(6,301) te
                                  Write(mfr,302) te
                               else
                                  !                       Write(6,*) "PREISNER: particle  unstable"
                                  !                       Write(6,301) te
                                  Write(ib,302) te
                               endif


                            ENDDO  !1
                         ENDDO  !2
                      ENDDO  !3
                   ENDDO  !4
                ENDDO  !6
             ENDDO  !5
          ELSE
             WRITE(6,*) "NO TRACDKING DONE : 6D CLOSED ORBIT WITH CAVITY NOT FOUND "
          ENDIF ! IEXT
          close(ib)
          close(mfr)

       case('PREISNERTRACK')
          READ(MF,*) NRAYS,POS, NTURN, FILENAME,filesmear, name
          READ(MF,*) CLOSED_ORBIT
          read(mf,*) title
          print*,"PREISNER TRACK", NRAYS,POS, NTURN, FILENAME
          print*,"PREISNER TRACK name ",  name
          call context(name)
          call context(FILENAME)
          call context(filesmear)
          call context(title)
          CLOSED=0.D0
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ring,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          else
             call move_to(my_ring,p1,POS)
          endif
          call kanalnummer(mfr)
          if(title(1:1)=='A') then
             open(unit=mfr,file=filename,position='APPEND',status='OLD')
          ELSE
             open(unit=mfr,file=filename)
          endif

          call kanalnummer(ib)

          if(title(1:1)=='A') then
             open(unit=ib,file=filesmear,position='APPEND',status='OLD')
          ELSE
             open(unit=ib,file=filesmear)
          endif


          DO I1=1,NRAYS
             track_k=.true.
             READ(MF,*) XP
             te=xp
             WRITE(MFR,301) XP
             IF(CLOSED_ORBIT) THEN
                CLOSED(5)=XP(5)
                CALL FIND_ORBIT(MY_RING,CLOSED,POS,MY_STATE,c_1d_5)
                IF(MY_STATE%NOCAVITY) THEN
                   X(1:4)=XP(1:4)+CLOSED(1:4)
                   X(5)=CLOSED(5)
                   X(6)=0.D0
                ELSE
                   XP=X+CLOSED
                ENDIF
             ELSE
                X=XP
             ENDIF

             DO I2=1,NTURN
                p=>p1
                do ii=1,MY_RING%n

                   call track(my_ring,X,ii,ii+1,MY_STATE)


                   call produce_aperture_flag(apertflag)
                   if (apertflag/=0) then
                      print *, 'Particle out of aperture!'

                      call ANALYSE_APERTURE_FLAG(apertflag,why)
                      Write(6,*) "PREISNER: particle : ",i1, " lost"
                      Write(6,301) te
                      Write(6,*) "PREISNER: APERTURE error for element: ",ii," name: ",p%MAG%name
                      Write(6,*) "Message: ",c_%message
                      !                      write(whymsg,*) 'APERTURE error: ',why
                      track_k=.false.
                      goto 2001
                      exit; !goes to the ne
                   endif

                   p=>p%next

                enddo

                IF(CLOSED_ORBIT) THEN
                   IF(MY_STATE%NOCAVITY) THEN
                      XP(1:4)=X(1:4)-CLOSED(1:4)
                      XP(5:6)=X(5:6)
                   ELSE
                      XP=X-CLOSED
                   ENDIF
                ELSE
                   XP=X
                ENDIF
                WRITE(MFR,301) XP

             ENDDO
2001         continue
             if (track_k) then
                Write(6,*) "PREISNER: particle : ",i1, " stable"
                Write(6,301) te
                Write(ib,302) te,i1
             endif

             write(6,*) " particle",i1 ," done "
          ENDDO
          write(6,*) " tracking done "
          CLOSE(MFR)
          CLOSE(ib)
       case('TRACKRAYS')
          READ(MF,*) NRAYS,POS, NTURN, FILENAME, name
301       FORMAT(6(1X,E15.8))
302       FORMAT(6(1X,E15.8),1x,i4)
          print*,"TRACKRAYS", NRAYS,POS, NTURN, FILENAME
          print*,"TRACKRAYS name ",  name
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ring,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif

          call kanalnummer(mfr)
          open(unit=mfr,file=filename)

          MY_STATE = MY_STATE + EXACTMIS0 + NOCAVITY0 + TIME + FRINGE
          call print(MY_STATE,6)
          EXACT_MODEL = my_true


          DO I1=1,NRAYS
             READ(MF,*) X
             WRITE(MFR,301) X
             DO I2=1,NTURN
                p=>MY_RING%start
                do ii=1,MY_RING%n
                   write(6,*) "##########################################"
                   write(6,'(i4, 1x,a)') ii,p%mag%name

                   call track(my_ring,X,ii,ii+1,MY_STATE)

                   write(6,'(6E8.4)') x

                   write(MFR,'(i4, 1x,a)') ii,p%mag%name
                   WRITE(MFR,301) X


                   call produce_aperture_flag(apertflag)
                   write(6,*) "apertflag ", apertflag
                   if (apertflag/=0) then
                      print *, 'Particle out of aperture!'

                      call ANALYSE_APERTURE_FLAG(apertflag,why)
                      Write(6,*) "ptc_trackline: APERTURE error for element: ",ii," name: ",p%MAG%name
                      Write(6,*) "Message: ",c_%message
                      write(whymsg,*) 'APERTURE error: ',why

                      exit; !goes to the ne
                   endif

                   p=>p%next

                enddo

             ENDDO
          ENDDO

          CLOSE(MFR)

       case('PRINTFRAMES')

          READ(MF,*) FILENAME
          CALL print_frames(MY_RING,filename)

       case('READFLATFILE')

          READ(MF,*) FILENAME
          CALL  READ_AND_APPEND_VIRGIN_LAYOUT(M_U,filename)

          WRITE(6,*) M_U%END%N, M_U%END%END%POS
       case('PRINTFLATFILE')

          READ(MF,*) FILENAME
          CALL  print_LAYOUT(MY_RING,filename)

          WRITE(6,*) M_U%END%N, M_U%END%END%POS

       case('PSREXAMPLEOFPATCHING')

          call APPEND_EMPTY_LAYOUT(m_u)
          CALL remove_drifts(MY_RING,m_u%END)
          m_u%end%name="psr_no_drift"
          call APPEND_EMPTY_LAYOUT(m_u)
          m_u%end%name="psr_quads_for_bends"
          CALL remove_drifts_bends(MY_RING,m_u%END)

          WRITE(6,*) MY_RING%N , m_u%END%N

       case('REMOVEPATCHESFORPIOTR')
          read(mf,*)do_not_remove, put_geo_patch
          read(mf,*) x_ref
          read(mf,*) x
          read(mf,*) sc,FILENAME

          call remove_patches(my_ring,x_ref,x,sc,do_not_remove,put_geo_patch,FILENAME)

       case('NORMALFORM')
          READ(MF,*)POS,name
          READ(MF,*) FILENAME

       case('TRANSLATELAYOUT')
          READ(MF,*)DT
          CALL TRANSLATE(MY_RING,DT)
       case('TRANSLATEPARTOFLAYOUT')
          READ(MF,*)DT
          READ(MF,*) I1,I2
          CALL TRANSLATE(MY_RING,DT,I1,I2)
          CALL MOVE_TO(MY_RING,P,i1)
          CALL FIND_PATCH(P%PREVIOUS,P,NEXT=my_TRUE,ENERGY_PATCH=my_FALSE)
          CALL MOVE_TO(MY_RING,P,i2)
          CALL FIND_PATCH(P,P%NEXT,NEXT=my_FALSE,ENERGY_PATCH=my_FALSE)
       case('ROTATEPARTOFLAYOUT')

          READ(MF,*)DT
          READ(MF,*) I1,I2
          CALL MOVE_TO(MY_RING,P,i1)
          call ROTATE_LAYOUT(my_ring,P%mag%p%f%ent,DT,I1,I2)
          CALL MOVE_TO(MY_RING,P,i1)
          CALL FIND_PATCH(P%PREVIOUS,P,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
          CALL MOVE_TO(MY_RING,P,i2)
          CALL FIND_PATCH(P,P%NEXT,NEXT=MY_FALSE,ENERGY_PATCH=MY_FALSE)

       case('TRANSLATEFIBREANDPATCH')
          READ(MF,*)POS
          READ(MF,*)DT
          CALL MOVE_TO(MY_RING,P,POS)
          CALL TRANSLATE_Fibre(P,DT,ORDER=1,BASIS=P%MAG%P%F%MID)
          CALL FIND_PATCH(P%PREVIOUS,P,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
          CALL FIND_PATCH(P,P%NEXT,NEXT=MY_FALSE,ENERGY_PATCH=MY_FALSE)
       case('POWERMULTIPOLE')
          READ(MF,*)POS
          READ(MF,*)n,cns, bend_like
          CALL MOVE_TO(MY_RING,P,POS)
          CALL ADD(P,N,0,CNS)
          p%mag%p%bend_fringe=bend_like
          p%magp%p%bend_fringe=bend_like
       case('COMPUTEMAP')
          READ(MF,*)POS,DEL,NO
          READ(MF,*) FILENAME

          CALL compute_map_general(MY_RING,my_state,filename,pos,del,no)

       case('ZEROSEXTUPOLES')
          call zero_sex(my_ring)
       case('POWERCAVITY')
          c_%CAVITY_TOTALPATH=0 ! fake pill box
          c_%phase0=0.0_dp ! because madx is crap anyway
          READ(MF,*)HARMONIC_NUMBER,VOLT,PHASE,epsf

          CALL power_cavity(MY_RING,HARMONIC_NUMBER,VOLT,PHASE,epsf)
       case('EQUILIBRIUMSIZES')
          READ(MF,*) POS,FILENAME,fileTUNE, NAME
          call context(name)
          if(name(1:11)/='NONAMEGIVEN') then
             posr=pos
             call move_to( my_ring,p,name,posR,POS)
             if(pos==0) then
                write(6,*) name, " not found "
                stop
             endif
          endif
          CALL radia(MY_RING,POS,FILENAME,fileTUNE)

       case('THINEXAMPLE')
          READ(MF,*) I1,I2,FILENAME

          mfr=0
          IF(FILENAME(1:6)=='SCREEN'.OR.FILENAME(1:6)=='screen') mfr=6
          if(mfr==0) then
             call kanalnummer(mfr)
             open(unit=mfr,file=filename)
          endif
          CALL THIN_EXAMPLE(MY_RING,MY_BEAMS,I1,I2,MY_STATE,MFR)
          IF(MFR/=6) close(mfr)
       case default

          write(6,*) " Command Ignored "

       end select


    enddo
100 continue
    write(6,*) " Exiting Command File ", ptc_fichier(1:len_trim(ptc_fichier))

    close(mf)

  END subroutine read_ptc_command



  subroutine power_cavity(r,HARM,VOLT,PHAS,prec)
    implicit none
    TYPE(LAYOUT), POINTER :: r
    type(fibre), pointer :: p
    integer i,ip,HARM
    type(internal_state) state
    real(dp) closed(6),s,VOLT,PHAS,accuracy,circum
    real(dp),optional :: prec

    call get_length(r,circum)

    write(6,*) " fiducial length ",circum

    state=my_state+nocavity0-totalpath0

    accuracy=1.0e-10_dp
    if(present(prec)) accuracy=prec
    closed=0.0_dp
    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)




    ip=0
    p=>r%start
    do i=1,r%n
       closed(6)=zero
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
          write(6,*) p%mag%name
          write(6,*) p%mag%volt
          write(6,*) p%mag%freq
          write(6,*) p%mag%phas
          write(6,*) "ref energy =",p%mag%p%p0c
          write(6,*) "electron =", c_%electron
          p%mag%volt=VOLT   !/p%mag%l
          p%magp%volt=p%mag%volt
          p%mag%freq=CLIGHT*HARM/circum   ! harmonic=120.0_dp
          p%magp%freq=p%mag%freq
          p%mag%phas=PHAS
          p%magp%phas=p%mag%phas
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
    r%charge= r%charge
    p=>r%start
    do i=1,r%n
       p%dir=-p%dir
       p=>p%next
    enddo

  end subroutine charge_dir

  SUBROUTINE radia(R,loc,FILE1,FILE2)
    USE madx_keywords
    implicit none
    TYPE(LAYOUT) R

    REAL(DP) X(6),m,as(6,6),energy,deltap
    TYPE(ENV_8) YS(6)
    type(beamenvelope) env
    CHARACTER(*) FILE1,FILE2
    type(normalform) normal
    integer nd2,npara,i,j,js(6),n1,n2
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    TYPE(INTERNAL_STATE) state
    integer no,loc,mf1,mf2
    real(dp) av(6,6),e(3)
    type(fibre), pointer :: p
    no=0
    call kanalnummer(mf1)
    open(mf1,file=FILE1)
    call kanalnummer(mf2)
    open(mf2,file=FILE2)

    state=(my_state-nocavity0)+radiation0
    x=0.d0;
    CALL FIND_ORBIT(R,X,loc,STATE,1.0e-5_dp)
    WRITE(6,*) " CLOSED ORBIT AT LOCATION ",loc
    write(6,*) x
    if(track_flag(r,x,loc,state)==0) then
       write(6,*) " stable closed orbit tracked "
    else
       write(6,*) " unstable closed orbit tracked "
       stop 333
    endif

    open(unit=30,file='junk.txt')
    call move_to(r,p,loc)
    do i=loc,loc+r%n
       call track(r,x,i,i+1,state)
       p=>p%next
       write(30,205) i,p%mag%name(1:8),x
    enddo
    close(30)
205 FORMAT(1x,i4,1x,a8,1x,6(1X,D18.11))

    call GET_loss(r,energy,deltap)
    write(6,*) x
    write(6,*) "energy loss: GEV and DeltaP/p0c ",energy,deltap

    CALL INIT(state,2,0,BERZ,ND2,NPARA)
    CALL ALLOC(Y);CALL ALLOC(NORMAL);call alloc(ys);call alloc(env);  ! ALLOCATE VARIABLES
    !Y=NPARA
    CALL ALLOC(ID)
    ID=1
    Y=X+ID
    ys=y

    CALL TRACK(R,YS,loc,state)
    if(.not.check_stable) write(6,*) " unstable tracking envelope "
    env%stochastic=.true.
    env=ys
    if(.not.check_stable) write(6,*) " unstable in normalizing envelope "

    y=ys
    normal=y
    if(.not.check_stable) write(6,*) " unstable in normalizing map "
    as=normal%a_t

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
    write(mf1,*) env%kick
    write(mf1,*) "B matrix"
    call print(env%STOCH,mf1)
    write(mf1,*) "m matrix"
    call print(ys%v,mf1)
    write(mf1,*) "equilibrium <X_i X_j>"
    call print(env%sij0,mf1)
    write(mf1,*) " Resonance Fluctuation "
    call print(env%bijnr,mf1)
    ! for etienne
    write(mf2,*) env%kick
    call print(env%STOCH,mf2)
    env%STOCH=env%STOCH**(-1)
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
    n1=1;n2=3;
    av(n1,n2)=(as(n1,1)*as(n2,1)+as(n1,2)*as(n2,2))*e(1)+(as(n1,3)*as(n2,3)+as(n1,4)*as(n2,4))*e(2)+ &
         (as(n1,5)*as(n2,5)+as(n1,6)*as(n2,6))*e(3)
    n1=3;n2=6;
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
    n1=1;n2=3;
    m=env%sij0.sub.'101000'
    write(mf1,*) " <x y> exact and alex "
    write(mf1,*) m/2.d0
    write(mf1,*) av(n1,n2)
    n1=3;n2=6;
    m=env%sij0.sub.'001001'
    write(mf1,*) " <y L> exact and alex "
    write(mf1,*) m/2.d0
    write(mf1,*) av(n1,n2)



    CALL KILL(Y);CALL KILL(Ys);CALL KILL(NORMAL);CALL KILL(env);

    ! compute map with radiation minus the cavity!
    ! cavity must be at the end and only one cavity

    if(no>0) then
       CALL INIT(STATE,no,0,BERZ,ND2,NPARA)
       CALL ALLOC(Y);  ! ALLOCATE VARIABLES
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

  end subroutine radia

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
    logical(lp) doneit,first
    real(dp) ent(3,3),a(3),ang(3),d(3)

    p=>r%start
    bend=>r%next%start%next   ! second layout in universe
    write(6,*) " using bend called ",bend%mag%name
    write(6,*) " 'USING SURVEY' TYPE 1 / 'USING GEOMETRY' TYPE 0 "
    READ(5,*) IG
    do i=1,r%n
       IF(P%MAG%KIND/=KIND0.AND.P%MAG%KIND/=KIND1.and.P%MAG%p%b0==zero) THEN

          CALL APPEND( NR, P )
       elseif(P%MAG%p%b0/=zero) then
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

  SUBROUTINE remove_patches(R,x0,x,sc,do_not_remove,patch_everywhere,FILENAME)  ! special example to be removed later
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer I,mf
    real(dp) x(6),x0(6),sc
    type(fibre), pointer :: P
    type(INTEGRATION_NODE), pointer :: T
    logical(lp) do_not_remove,patch_everywhere
    type(three_d_info) v
    character(*) FILENAME

    if(do_not_remove) then
       p=>r%start
       do i=1,r%n
          p%patch=-1
          p%patch=0

          p=>p%next
       enddo

       if(patch_everywhere) then

          p=>r%start
          do i=1,r%n-1
             CALL FIND_PATCH(P,P%next,NEXT=MY_TRUE,ENERGY_PATCH=MY_FALSE)
             p=>p%next
          enddo

       endif
    endif
    CALL MAKE_NODE_LAYOUT(MY_RING)

    CALL fill_survey_data_in_NODE_LAYOUT(MY_RING)

    T=>MY_RING%T%start

    call alloc_three_d_info(v)

    V%reference_ray=x0
    V%X=x
    V%SCALE=sc
    call kanalnummer(mf)
    OPEN(UNIT=mf,FILE=FILENAME)
    t=>my_ring%t%start
    DO I=1,MY_RING%T%n
       CALL TRACK( my_ring,V,DEFAULT,POS1=I,POS2=I+1 )
       IF(V%U(1)) THEN
          WRITE(6,*) " UNSTABLE ",I
          goto 12
       ENDIF
       WRITE(mf,*) V%r(3),V%r(1),t%parent_fibre%mag%name,t%cas
       t=>t%next
    ENDDO
12  CLOSE(mf)


  end SUBROUTINE remove_patches




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
    call print_frames(MY_RING,filename)

  end subroutine printframes



end module pointer_lattice



subroutine read_ptc_command77(ptc_fichier)
  use pointer_lattice
  implicit none
  character(*) ptc_fichier
  integer m

  call kanalnummer(m)

  open(unit=m,file=ptc_fichier,status='old',err=2001)
  close(m)
  call read_ptc_command(ptc_fichier)
  return
2001 continue

  write(6,*) " Warning: command file does not exit "

end  subroutine read_ptc_command77



!Piotr says:
!I have implemented for you MAD-X command ptc_script.
!It has one parameter named "file". Hence, you call sth like
!ptc_script, file="directptctweaks.ptc";
!It executes subroutine execscript in file madx_ptc_script.f90.

subroutine gino_ptc_command77(gino_command)
  use pointer_lattice
  implicit none
  character(*) gino_command

  !  Etienne puts garbage here

  call context(gino_command)
  call call_gino(gino_command)

end  subroutine gino_ptc_command77
