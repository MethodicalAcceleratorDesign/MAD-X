module pointer_lattice
  use madx_keywords
  USE gauss_dis
  implicit none
  public
  ! stuff for main program
  type(layout),pointer :: my_ring
  type(internal_state),pointer :: my_state
  integer ,pointer :: my_start, MY_ORDER, MY_NP,MY_END
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
  INTEGER :: N_BEAM=0,USE_BEAM=1
contains

  subroutine read_ptc_command(ptc_fichier)
    use madx_ptc_module, piotr_state=>my_state
    implicit none
    CHARACTER*(120) com,COMT,filename,name_root,title,name_root_res,filetune,FILESMEAR
    character*(4) suffix,SUFFIX_res
    character(*) ptc_fichier
    integer i,mf,i_layout_temp,LIM(2),IB,NO
    !  FITTING FAMILIES
    INTEGER NPOL,J,NMUL,K,ICN,N,np,MRESO(3)
    type(pol_block), ALLOCATABLE :: pol_(:)
    CHARACTER*(NLP) NAME,flag
    real(dp) targ_tune(2),targ_chrom(2),epsf
    real(dp) targ_RES(4)
    !  END FITTING FAMILIES
    real(dp),allocatable :: beta(:,:,:)
    REAL(DP) DBETA,tune(3),tunenew(2),CHROM(2),DEL
    ! fitting and scanning tunes
    real(dp) tune_ini(2),tune_fin(2),dtu(2)
    integer nstep(2),i1,i2,I3,neq
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
    INTEGER NUMBER_OF_PARTICLE
    TYPE(DAMAP) MY_A
    INTEGER MY_A_NO,MY_A_ND

    type(internal_state),target :: my_default
    save my_default

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
             write(6,*) "Selected Layout",i_layout
             i_layout=i_layout_temp
             call move_to_layout_i(m_u,my_ring,i_layout)
          endif

       case('ALLOCATEBETA')
          ALLOCATE(BETA(2,2,MY_RING%N))

       case('DEALLOCATEBETA')
          DEALLOCATE(BETA)
       case('FILLBETA')
          READ(MF,*) IB,pos
          CALL FILL_BETA(MY_RING,MY_STATE,pos,BETA,IB,DBETA,tune,tunenew)
       case('FITBENDS')
          CALL fit_all_bends(MY_RING,MY_STATE)
       case('THINLENS')
          READ(MF,*) THIN
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(MY_RING,THIN)
       case('EVENTHINLENS')
          READ(MF,*) THIN
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(MY_RING,THIN,EVEN=.TRUE.)
       case('ODDTHINLENS')
          READ(MF,*) THIN
          WRITE(6,*) "THIN LENS FACTOR =",THIN
          CALL THIN_LENS_resplit(MY_RING,THIN,EVEN=.FALSE.)
          ! thin layout stuff
       case('MAKE_THIN_LAYOUT','MAKETHINLAYOUT')

          CALL MAKE_THIN_LAYOUT(MY_RING)


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
          IF(FILENAME(1:2)=='NO') THEN
             CALL  create_PANCAKE(MY_BEAMS(USE_BEAM),NUMBER_OF_PARTICLE,CUT,SIG,my_ring%start%t1)
          ELSE
             CALL  kanalnummer(mfr)
             OPEN(UNIT=MFR,FILE=FILENAME)
             READ(MFR,*) MY_A_NO,MY_A_ND
             CALL INIT(MY_A_NO,MY_A_ND,0,0)
             CALL ALLOC(MY_A)
             CALL READ(MY_A,MFR)
             CALL  create_PANCAKE(MY_BEAMS(USE_BEAM),N_BEAM,CUT,SIG,my_ring%start%t1,MY_A)
             CALL KILL(MY_A)
             close(mfr)
          ENDIF
       case('COPYBEAM')
          READ(MF,*) I1,I2
          CALL COPY_BEAM(MY_BEAMS(I1),MY_BEAMS(I2))
          ! END BEAMS STUFF

          ! end of layout stuff
          ! random stuff
       case('GAUSSIANSEED')
          READ(MF,*) I1
          CALL gaussian_seed(i1)
       case('RANDOMMULTIPOLE')
          read(mf,*) name
          read(mf,*) n,cns,cn
          read(mf,*) addi,integrated
          read(mf,*) cut
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
             POL_(J)=0
             POL_(J)%NAME=NAME
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
       case('FITSEXTUPOLERESONANCE')
          read(mf,*) epsf,neq
          call context(flag)
          READ(MF,*) MRESO
          read(mf,*) targ_RES
          call lattice_fit_SEXT_RES_from_a_gmap_vec(my_ring,my_state,EPSF,MRESO,pol_,NPOL,targ_RES,NP,neq)
       case('FITSEXTUPOLERESONANCE1')
          read(mf,*) epsf,neq
          call context(flag)
          READ(MF,*) MRESO
          read(mf,*) targ_RES
          call lattice_fit_SEXT_RES_from_a_gmap_vec1(my_ring,my_state,EPSF,MRESO,pol_,NPOL,targ_RES,NP,neq)
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

end module pointer_lattice

subroutine read_ptc_command77(ptc_fichier)
  use pointer_lattice
  implicit none
  character(*) ptc_fichier

  call read_ptc_command(ptc_fichier)


end  subroutine read_ptc_command77



!Piotr says:
!I have implemented for you MAD-X command ptc_script.
!It has one parameter named "file". Hence, you call sth like
!ptc_script, file="directptctweaks.ptc";
!It executes subroutine execscript in file madx_ptc_script.f90.
