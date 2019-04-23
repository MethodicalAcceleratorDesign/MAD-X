!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module madx_keywords
  use S_fitting
  implicit none
  public
  logical(lp)::mad8=my_false
  integer :: ifield_name=0
  logical(lp),private :: do_survey =my_false
  logical(lp) :: print_marker =my_true
  type(tree_element), private, allocatable :: t_e(:),t_ax(:),t_ay(:)
  real(dp), private :: a_(3),ent_(3,3), b_(3),exi_(3,3)
  logical :: old_name_vorname = .false.
  logical :: readingmaps = .true.
  type keywords
     character*20 magnet
     character*20 model
     logical(lp) FIBRE_flip
     INTEGER FIBRE_DIR
     integer method
     integer nstep
     logical(lp) exact
     logical(lp) madLENGTH
     logical(lp) mad8
     real(dp) tiltd
     type(el_list) LIST
  end type keywords

  type MADX_SURVEY
     REAL(DP) ALPHA,TILT,LD
     REAL(DP) PHI,THETA,PSI
     TYPE(CHART) CHART
  END type MADX_SURVEY

  include "a_namelists.inc"


 ! INTERFACE read_lattice_append
!     MODULE PROCEDURE read_universe_database 
 ! END  INTERFACE


contains


  subroutine create_fibre_append(append,mylat,key,EXCEPTION,magnet_only,br)  
    implicit none

!    type(mad_universe), target, intent(inout)  :: m_u
    type(layout), target, intent(inout)  :: mylat
    logical(lp), optional :: magnet_only
    type(keywords) key
    INTEGER EXCEPTION  !,NSTD0,METD0
    logical(lp) doneit,append
    type(fibre), pointer :: current
    type (taylor),optional, INTENT(INout):: br(:,:)

    if(append) then
     call append_empty(mylat)
    else
     if(associated(mylat%end)) then
      IF(ASSOCIATED(mylat%T)) THEN
         CALL kill_NODE_LAYOUT(mylat%T)  !  KILLING THIN LAYOUT
         nullify(mylat%T)
        if(lielib_print(12)==1) WRITE(6,*) " NODE LAYOUT HAS BEEN KILLED "
       ENDIF      
        mylat%end=-1
       else
        call append_empty(mylat)
     endif
    endif
     call  create_fibre(mylat%end,key,EXCEPTION,magnet_only,br)
     
    if(.not.append) then
     mylat%closed=my_true

     doneit=my_true
     call ring_l(mylat,doneit)

     call survey(mylat)
     call MAKE_NODE_LAYOUT( mylat)     
    endif
  end subroutine create_fibre_append


  subroutine create_fibre(el,key,EXCEPTION,magnet_only,br)
    implicit none
    integer ipause, mypause,i
    type(fibre), target, intent(inout)::el
    logical(lp), optional :: magnet_only
    type(keywords) key
    type(el_list) blank
    character*255 magnet
    character*17 MODEL
    INTEGER EXCEPTION  !,NSTD0,METD0
    LOGICAL(LP) EXACT0,magnet0
    logical(lp) FIBRE_flip0,MAD0
    logical(lp) :: t=my_true,f=my_false
    INTEGER FIBRE_DIR0,IL
    real(dp) e1_true,norm
    type (taylor),optional, INTENT(INout):: br(:,:)

    IL=15

    if(present(magnet_only)) then
       magnet0=magnet_only
    else
       magnet0=my_false
    endif

    blank=0
    magnet=key%magnet
    call context(magnet)
    model=key%model
    call context(model)

    CALL SET_MADX_(t,magnet0)


    select case(MODEL)
    CASE("DRIFT_KICK       ")
       MADTHICK=drift_kick_drift
    CASE("MATRIX_KICK      ")
       MADTHICK=matrix_kick_matrix
    CASE("DELTA_MATRIX_KICK")
       MADTHICK=kick_sixtrack_kick
    CASE DEFAULT
 
       EXCEPTION=1
       ipause=mypause(444)
       RETURN
    END SELECT
    !   MADTHICK=drift_kick_drift
    !    NSTD0=NSTD
    !    METD0=METD
    EXACT0=EXACT_MODEL
    FIBRE_FLIP0= FIBRE_FLIP
    FIBRE_DIR0=FIBRE_DIR
    MAD0=MAD

    KEY%LIST%nst=KEY%NSTEP
    KEY%LIST%method=KEY%METHOD
    EXACT_MODEL=KEY%EXACT
    FIBRE_FLIP = KEY%FIBRE_FLIP
    FIBRE_DIR  = KEY%FIBRE_DIR
    MADLENGTH=KEY%MADLENGTH

    !     real(dp) L,LD,LC,K(NMAX),KS(NMAX)
    !     real(dp) ang(3),t(3)
    !     real(dp) angi(3),ti(3)
    !     integer patchg
    !     real(dp) T1,T2,B0
    !     real(dp) volt,freq0,harmon,lag,DELTA_E,BSOL
    !     real(dp) tilt
    !     real(dp) FINT,hgap,h1,h2,X_COL,Y_COL
    !     real(dp) thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
    !     CHARACTER(120) file
    !     CHARACTER(120) file_rev
    !    CHARACTER(nlp) NAME
    !     CHARACTER(vp) VORNAME
    !     INTEGER KIND,nmul,nst,method
    !     LOGICAL(LP) APERTURE_ON
    !     INTEGER APERTURE_KIND
    !     REAL(DP) APERTURE_R(2),APERTURE_X,APERTURE_Y
    !     LOGICAL(LP) KILL_ENT_FRINGE,KILL_EXI_FRINGE,BEND_FRINGE,PERMFRINGE
    !     REAL(DP) DPHAS,PSI,dvds
    !     INTEGER N_BESSEL

    if(sixtrack_compatible) then
       EXACT_MODEL=my_false
       KEY%LIST%method=2
       MADTHICK=drift_kick_drift
    endif


    SELECT CASE(magnet(1:IL))
    CASE("DRIFT          ")
       BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("SUPERDRIFT     ")
       BLANK=SUPERDRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("SOLENOID       ")
       if(sixtrack_compatible) stop 1
       if(KEY%LIST%L/=0.0_dp) then
          BLANK=SOLENOID(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
          !  BLANK%bend_fringe=key%list%bend_fringe
       else
          write(6,*) "switch solenoid to dubious thin multipole "
          BLANK=MULTIPOLE_BLOCK(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       endif

    CASE("THICKMULTIPOLE ")
       BLANK=multipoleTILT(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe

    CASE("QUADRUPOLE     ")
       BLANK=QUADRUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("SEXTUPOLE     ")
       BLANK=SEXTUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("OCTUPOLE      ")
       BLANK=OCTUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("SBEND         ")
       BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("TRUERBEND     ")
       if(sixtrack_compatible) stop 2

       e1_true= KEY%LIST%b0/2.0_dp+ KEY%LIST%t1
       BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,e1=e1_true,list=KEY%LIST)

    CASE("WEDGRBEND     ")
       if(sixtrack_compatible) stop 3

       BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,e1=KEY%LIST%t1,e2=KEY%LIST%t2,list=KEY%LIST)

    CASE("RBEND         ")
       if(sixtrack_compatible) stop 4
       KEY%LIST%T1=KEY%LIST%T1+KEY%LIST%B0/2.0_dp
       KEY%LIST%T2=KEY%LIST%T2+KEY%LIST%B0/2.0_dp
       BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("KICKER         ","VKICKER        ","HKICKER        ")
       BLANK=KICKER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("MONITOR        ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       endif
    CASE("HMONITOR        ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND12;
       endif
    CASE("VMONITOR       ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND13;
       endif
    CASE("INSTRUMENT     ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND14;
       endif
    CASE("MARKER         ")
       BLANK=MARKER(KEY%LIST%NAME,list=KEY%LIST)
    CASE("CHANGEREF      ")
       if(sixtrack_compatible) stop 5
       BLANK=CHANGEREF(KEY%LIST%NAME,KEY%LIST%ANG,KEY%LIST%T,KEY%LIST%PATCHG)
    CASE("RFCAVITY       ")
       if(sixtrack_compatible) then
          If(KEY%LIST%L/=0.0_dp) stop 60
          If(KEY%LIST%N_BESSEL/=0.0_dp) stop 61
          norm=0.0_dp
          do i=1,nmax
             norm=norm+abs(KEY%LIST%k(i))+abs(KEY%LIST%ks(i))
          enddo
          norm=norm-abs(KEY%LIST%k(2))
          if(norm/=0.0_dp) then
             write(6,*) norm
             stop 62
          endif
       endif
       BLANK=RFCAVITY(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("TWCAVITY       ")
       if(sixtrack_compatible) stop 7
       BLANK=TWCAVITY(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("ELSEPARATOR    ")
       if(sixtrack_compatible) stop 8
       BLANK=ELSEPARATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("MULTIPOLE_BLOCK","MULTIPOLE      ")
       BLANK=MULTIPOLE_BLOCK(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("SMI            ","SINGLE_LENS    ")
       if(sixtrack_compatible) stop 9
       BLANK=SINGLE_LENS(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("RCOLLIMATOR    ")
       if(sixtrack_compatible) stop 10
       BLANK=RCOLLIMATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("ECOLLIMATOR    ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=ECOLLIMATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       endif
    CASE("WIGGLER        ")
       if(sixtrack_compatible) stop 12
       BLANK=WIGGLER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("HELICALDIPOLE  ")
       if(sixtrack_compatible) stop 13
       BLANK=HELICAL(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("PANCAKE        ")
       if(sixtrack_compatible) stop 13
       BLANK=pancake(KEY%LIST%NAME,KEY%LIST%file)
    CASE("ABELL_DRAGT    ")
       if(sixtrack_compatible) stop 13
       BLANK=abell_dragt(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("INTERNALPANCAKE")
       if(sixtrack_compatible) stop 13
       BLANK=pancake(KEY%LIST%NAME,br=br)

    CASE DEFAULT
       WRITE(6,*) " "
       WRITE(6,*) " THE MAGNET"
       WRITE(6,*) " "
       WRITE(6,*) "  --->   ",MAGNET(1:IL)
       WRITE(6,*) " "
       WRITE(6,*)  " IS NOT PERMITTED "
       STOP 666
    END SELECT

    BLANK%VORNAME = KEY%LIST%VORNAME
    CALL EL_Q_FOR_MADX(EL,BLANK)
    !  added 2007.07.09
    el%mag%parent_fibre =>el
    el%magp%parent_fibre=>el
    !  end of added 2007.07.09


    CALL SET_MADX_(f,f)

    !    NSTD=NSTD0
    !    METD=METD0
    EXACT_MODEL=EXACT0
    FIBRE_FLIP= FIBRE_FLIP0
    FIBRE_DIR=FIBRE_DIR0
    MAD=MAD0

    !    IF(ASSOCIATED(EL%PREVIOUS)) THEN
    !     if(.not.associated(EL%POS))allocate(EL%POS)
    !     EL%POS=EL%PREVIOUS%POS+1
    !    ELSE
    !     if(.not.associated(EL%POS))allocate(EL%POS)
    !     EL%POS=1
    !    ENDIF
    if(key%list%BEND_FRINGE) then
       el%mag%p%bend_fringe=my_true
       el%magp%p%bend_fringe=my_true
    endif

    if(el%mag%kind==kind4) then
       el%mag%c4%CAVITY_TOTALPATH=key%list%CAVITY_TOTALPATH
       el%magp%c4%CAVITY_TOTALPATH=key%list%CAVITY_TOTALPATH
    endif

  end subroutine create_fibre

  subroutine zero_key(key)
    implicit none

    type(keywords) , intent(out):: key
    key%magnet="CROTTE"
    select case(MADTHICK)
    CASE(drift_kick_drift)
       key%model="DRIFT_KICK       "
    CASE(matrix_kick_matrix)
       key%model="MATRIX_KICK      "
    CASE(kick_sixtrack_kick)
       key%model="DELTA_MATRIX_KICK"
    END SELECT


    key%FIBRE_flip=FIBRE_flip
    key%FIBRE_DIR=FIBRE_DIR
    key%method=METD
    key%nstep=NSTD
    key%exact=EXACT_MODEL
    key%madLENGTH=madLENGTH
    key%LIST%NMUL = 1
    key%mad8 = mad8
    key%tiltd=0.0_dp
    key%LIST=0

  end subroutine zero_key

  !  PRINTING FIBRES FOR FLAT FILES
  subroutine print_COMPLEX_SINGLE_STRUCTURE(L,FILENAME,LMAX0,NL)
    implicit none
    character(*) filename
    integer I,MF,N,n_l
    type(LAYOUT), TARGET :: L
    type(LAYOUT), pointer :: CL
    REAL(DP),OPTIONAL :: LMAX0
    integer,OPTIONAL :: NL

    n_l=0
    if(present(nl)) n_l=nl
    call kanalnummer(mf)
    open(unit=mf,file=filename)
    IF(ASSOCIATED(L%DNA)) THEN
       N=SIZE(L%DNA)
       Write(mf,*) N,N_L, " Number of pointers in the DNA array pointed at layouts"

       DO I=1,N
          L%DNA(I)%L%INDEX=I
          CALL print_LAYOUT(L%DNA(I)%L,FILENAME,LMAX0,MF)
       ENDDO
       !       ENDIF

       !       write(mf,*) " Beam Line DNA structure "
       !       do i=1,N
       !        ncon1=0
       !        ncon2=0
       !        if(associated(L%DNA(i)%L%con1)) then
       !         ncon1=size(L%DNA(i)%L%con1)
       !         ncon2=size(L%DNA(i)%L%con2)
       !          write(mf,*) ncon1,ncon2


       !          do j=1,max(ncon1,ncon2)
       !           if(j>ncon1) then
       !             write(mf,*) 0, L%DNA(i)%L%con2(j)%pos
       !            elseif(j>ncon2) then
       !             write(mf,*)  L%DNA(i)%L%con1(j)%pos,0
       !            else
       !             write(mf,*)  L%DNA(i)%L%con1(j)%pos, L%DNA(i)%L%con2(j)%pos
       !           endif
       !          enddo
       !        else
       !          write(mf,*) ncon1,ncon2
       !        endif
       !      enddo

       !       write(mf,*) " End of Beam Line DNA structure "

    ENDIF

    CL=>L


    if(n_l>0) then
       do i=1,n_l
          ! write(mf,*) " Beam Line DNA structure "
          ! write(mf,*) " End of Beam Line DNA structure "
          CALL print_LAYOUT(CL,FILENAME,LMAX0,MF)
          CL=>CL%NEXT
       enddo
    else
       CALL print_LAYOUT(L,FILENAME,LMAX0,MF)
    endif

    CLOSE(MF)
  END SUBROUTINE print_COMPLEX_SINGLE_STRUCTURE

  subroutine print_LAYOUT(L,FILENAME,LMAX0,MFF)
    implicit none
    character(*) filename
    integer I,MF,nmark
    INTEGER, OPTIONAL :: MFF
    type(LAYOUT), TARGET :: L
    type(FIBRE), pointer :: P
    REAL(DP),OPTIONAL :: LMAX0
    character*255 line
    logical(lp) print_temp
nmark=0
  if(.not.print_marker) then ! counting markers if not printing them
    P=>L%START
    DO I=1,L%N
       if(.not.(print_marker.or.p%mag%kind/=kind0.or.i==1)) then  
        nmark=nmark+1
       endif
       P=>P%NEXT
    ENDDO 
    if(L%START%mag%kind==kind0) then
        write(6,*) "Removing ",nmark, "markers (first 1.0_dp left in) "
    else
        write(6,*) "Removing ",nmark, "markers "
    endif
  endif  ! counting markers
  
    IF(PRESENT(MFF)) THEN
       MF=MFF
    ELSE
       call kanalnummer(mf)
       open(unit=mf,file=filename)
    ENDIF

    IF(PRESENT(LMAX0)) THEN
       WRITE(MF,*) L%N-nmark, LMAX0, " NUMBER OF FIBRES AND L_MAX  "
    ELSE
       WRITE(MF,*) L%N-nmark, 0, " NUMBER OF FIBRES AND L_MAX  "
    ENDIF
    if(l%name(1:1)/=' ') then
       write(MF,'(a17,a16)') " GLOBAL DATA FOR ",l%name
    else
       write(MF,*) " $$$$$$$$$ GLOBAL DATA  $$$$$$$$$"
    endif

    write(line,*) l%start%mass,L%START%mag%p%p0c,l%start%ag, " MASS, P0C, AG(spin)"
    write(MF,'(a255)') line
    write(line,*) phase0,compute_stoch_kick,l%start%charge, " PHASE0, compute_stoch_kick, CHARGE"
    write(MF,'(a255)') line
!     write(MF,*) phase0,compute_stoch_kick,l%start%charge, " PHASE0, compute_stoch_kick, CHARGE"
    write(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING, &
         "CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING"
    write(line,*) sector_nmul_max,SECTOR_NMUL,&
         OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE,&
         " SECTOR_NMUL_MAX,SECTOR_NMUL,OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE"
    write(mf,'(a255)')line
    write(line,*) wedge_coeff,valishev, " wedge_coeff", " Valishev Multipole "
    write(mf,'(a255)')line
    write(MF,*) MAD8_WEDGE, " MAD8_WEDGE"
    write(MF,*) " $$$$$$$ END GLOBAL DATA $$$$$$$$$$"


    P=>L%START
    DO I=1,L%N
       if(i==1) then
          print_temp=print_frame
          print_frame=my_true
       endif
       
       if(print_marker.or.p%mag%kind/=kind0.or.i==1) then  
        CALL print_FIBRE(P,mf)
       endif
        if(i==1) then
          print_frame=print_temp
       endif
       P=>P%NEXT
    ENDDO

    IF(.NOT.PRESENT(MFF)) CLOSE(MF)

  END subroutine print_LAYOUT


  subroutine READ_INTO_VIRGIN_LAYOUT(L,FILENAME,RING,LMAX0,mf1)
    implicit none
    character(*) filename
    integer mf,I,N,RES,se1,se2
    integer, optional :: mf1
    type(LAYOUT), TARGET :: L
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0
    LOGICAL(LP) RING_IT,doneit
    character*255 line
    character*255 lineg
    real(dp) p0c,MASSF,ag0
    type(internal_state) original
    logical ttt,uuu

    RING_IT=MY_TRUE

    IF(PRESENT(RING)) RING_IT=RING

    if(present(mf1)) then
       mf=mf1
    else
       call kanalnummer(mf)
       open(unit=mf,file=filename,status='OLD',err=2001)
    endif

    IF(PRESENT(LMAX0)) then
       READ(MF,*) N,LMAX0
    ELSE
       READ(MF,*) N
    ENDIF
    read(MF,'(a255)') line
    call context(line)

    if(index(line,"FOR")/=0) then
       l%name=line(index(line,"FOR")+3:index(line,"FOR")+2+nlp)
    endif
    read(MF,'(A255)') lineg
    res=INDEX (lineG, "AG(spin)")
    IF(RES==0) THEN
       read(lineg,*) MASSF,p0c
       IF(ABS(MASSF-pmap)/PMAP<0.01E0_DP) THEN
          A_PARTICLE=A_PROTON
       ELSEIF(ABS(MASSF-pmae)/pmae<0.01E0_DP) THEN
          A_PARTICLE=A_ELECTRON
       ELSEIF(ABS(MASSF-pmaMUON)/pmaMUON<0.01E0_DP) THEN
          A_PARTICLE=A_MUON
       ENDIF
    ELSE
       read(lineg,*) MASSF,p0c,A_PARTICLE
    ENDIF
    ag0=A_PARTICLE
    read(MF,*) phase0,compute_stoch_kick,initial_charge
    read(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING
    read(MF,*) se1,se2,OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE
    call input_sector(se2,se1)

    read(MF,'(A255)') lineg
    res=INDEX (lineG, "Valishev")
    IF(RES==0) THEN
       read(lineg,*) wedge_coeff
    ELSE
       read(lineg,*) wedge_coeff,valishev
    ENDIF

    read(MF,*) MAD8_WEDGE
    read(MF,'(a255)') line
    original=default
  !  if(allocated(s_b)) then
       firsttime_coef=my_true
  !     deallocate(s_b)
  !  endif
    !    L%MASS=MASSF
    MASSF=MASSF/pmae
    CALL MAKE_STATES(MASSF)
    A_PARTICLE=ag0
    default=original
    call Set_madx(p0c=p0c)
    DO I=1,N
       CALL APPEND_CLONE(L,muonfactor=massf,charge=initial_charge)
       CALL READ_FIBRE(L%END,mf)
       CALL COPY(L%END%MAG,L%END%MAGP)
    ENDDO

    if(.not.present(mf1)) CLOSE(MF)

    L%closed=RING_IT

    doneit=my_true
    call ring_l(L,doneit)
    ! if(do_survey) call survey(L)

    return

2001 continue

    Write(6,*) " File ",filename(1:len_trim(filename)) ," does not exist "

  END subroutine READ_INTO_VIRGIN_LAYOUT

  subroutine READ_AND_APPEND_VIRGIN_general(U,filename,RING,LMAX0)
    implicit none
    character(*) filename
    integer  mf
    type(MAD_UNIVERSE), TARGET :: U
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0
    character*120 line
    integer res

    res=0
    call kanalnummer(mf)
    open(unit=mf,file=filename,status='OLD',err=2001)
    read(mf,'(a120)') line
    res=INDEX (line, "DNA")
    if(res/=0) res=1
    close(mf)

    if(res==1) then
 
       call read_COMPLEX_SINGLE_STRUCTURE(U,filename,RING,LMAX0)
    else

       call APPEND_EMPTY_LAYOUT(U)
 
       CALL READ_INTO_VIRGIN_LAYOUT(U%END,FILENAME,RING,LMAX0)
       ! if(do_survey) call survey(u%end)
    endif
    do_survey=my_false
    return
2001 continue

    Write(6,*) " File ",filename(1:len_trim(filename)) ," does not exist "

  END subroutine READ_AND_APPEND_VIRGIN_general

  subroutine READ_AND_APPEND_VIRGIN_LAYOUT(U,filename,RING,LMAX0,mf)
    implicit none
    character(*) filename
    integer,optional:: mf
    type(MAD_UNIVERSE), TARGET :: U
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0


    call APPEND_EMPTY_LAYOUT(U)

    CALL READ_INTO_VIRGIN_LAYOUT(U%END,FILENAME,RING,LMAX0=LMAX0,mf1=MF)

  END subroutine READ_AND_APPEND_VIRGIN_LAYOUT

  subroutine print_FIBRE(m,mf)
    implicit none
    integer mf,siam_pos,siam_index,GIRD_POS,GIRD_index
    type(FIBRE), pointer :: m
    siam_pos=0
    siam_index=0
    GIRD_POS=0
    GIRD_index=0
    if(associated(m%mag%siamese)) then
       siam_index=m%mag%siamese%parent_fibre%parent_layout%index
       siam_pos=m%mag%siamese%parent_fibre%pos
    endif
    if(associated(m%mag%GIRDERS)) then
       GIRD_index=m%mag%GIRDERS%parent_fibre%parent_layout%index
       GIRD_POS=m%mag%GIRDERS%parent_fibre%pos
    endif
    WRITE(MF,*) " @@@@@@@@@@@@@@@@@@@@ FIBRE @@@@@@@@@@@@@@@@@@@@"
    if(siam_index==0.AND.GIRD_index==0) then
       WRITE(MF,'(A11,4(I4,1x))') " DIRECTION ", M%DIR, &
            m%mag%parent_fibre%parent_layout%index,m%mag%parent_fibre%pos, &
            m%mag%parent_fibre%parent_layout%n
    else

       WRITE(MF,'(A11,4(I4,1x),A16,4(I4,1x))') " DIRECTION ", M%DIR, &
            m%mag%parent_fibre%parent_layout%index,m%mag%parent_fibre%pos, &
            m%mag%parent_fibre%parent_layout%n," Siamese/Girder "         &
            ,siam_pos,siam_index,GIRD_POS,GIRD_index
    endif
    CALL print_chart(m%CHART,mf)
    CALL print_PATCH(m%PATCH,mf)
    CALL print_element(M,M%MAG,mf)
    WRITE(MF,*) " @@@@@@@@@@@@@@@@@@@@  END  @@@@@@@@@@@@@@@@@@@@"

  END subroutine print_FIBRE

  subroutine READ_FIBRE(m,mf)
    implicit none
    integer mf
    type(FIBRE), pointer :: m
    character*255 line
    READ(MF,*) LINE
    READ(MF,'(A11,I4)') LINE(1:11),M%DIR
    CALL READ_chart(m%CHART,mf)
    CALL READ_PATCH(m%PATCH,mf)
    CALL READ_element(m,M%MAG,mf)
    READ(MF,*) LINE

  END subroutine READ_FIBRE

  subroutine READ_FIBRE_2_lines(mf,DIR,index,pos,n,siam_index,siam_pos,gird_index,gird_pos)
    implicit none
    integer mf
    character*255 line

    integer DIR,index,pos,n,siam_index,siam_pos,gird_index,gird_pos
    READ(MF,*) LINE
    siam_index=0
    siam_pos=0
    gird_index=0
    gird_pos=0
    READ(MF,'(A11,4(I4,1x),A16,4(I4,1x))') LINE(1:11),DIR,index,pos,n, &
         LINE(12:27),siam_pos,siam_index,gird_index,gird_pos
    !    CALL READ_chart(m%CHART,mf)
    !    CALL READ_PATCH(m%PATCH,mf)
    !    CALL READ_element(M%MAG,mf)
    !    READ(MF,*) LINE

  END subroutine READ_FIBRE_2_lines

  subroutine print_PATCH(m,mf)
    implicit none
    integer mf,i1,i2,i3
    type(PATCH), pointer :: m
    character*255 line

    i1=M%PATCH
    i2=M%energy
    i3=M%time

    IF(IABS(i1)+iabs(i2)+iabs(i3)/=0) then
       WRITE(MF,*) " >>>>>>>>>>>>>>>>>> PATCH <<<<<<<<<<<<<<<<<<"
       WRITE(MF,*) M%PATCH,M%ENERGY,M%TIME," patch,energy,time"
       WRITE(MF,*) M%A_X1,M%A_X2,M%B_X1,M%B_X2," discrete 180 rotations"
       WRITE(LINE,*) M%A_D,M%A_ANG,"  a_d, a_ang "
       WRITE(MF,'(A255)') LINE
       WRITE(LINE,*) M%B_D,M%B_ANG,"  b_d, b_ang "
       WRITE(MF,'(A255)') LINE
       WRITE(MF,*) M%A_T,M%B_T,"  time patches a_t and b_t "
       WRITE(MF,*) M%A_L,M%B_L,"  patch length patches a_t and b_t "
       WRITE(MF,*) " >>>>>>>>>>>>>>>>>>  END  <<<<<<<<<<<<<<<<<<"
    else
       WRITE(MF,*) " NO PATCH "
    endif
  END subroutine print_PATCH

  subroutine READ_PATCH(m,mf)
    implicit none
    integer mf
    type(PATCH), pointer :: m
    character*255 line

    READ(MF,*)LINE
    if(index(line,"NO")==0) then
       READ(MF,*) M%PATCH,M%ENERGY,M%TIME
       READ(MF,*) M%A_X1,M%A_X2,M%B_X1,M%B_X2
       READ(MF,*) M%A_D,M%A_ANG
       READ(MF,*) M%B_D,M%B_ANG
       READ(MF,*) M%A_T,M%B_T !M%A_L,M%B_L
M%A_L=M%A_T
M%B_L=M%B_T
       READ(MF,*) LINE
    endif

  END subroutine READ_PATCH

  subroutine print_chart(m,mf)
    implicit none
    integer mf,I
    type(CHART), pointer :: m
    character*255 line
    real(dp) norm

    norm=0.0_dp
    do i=1,3
       norm=abs(M%D_IN(i))+norm
       norm=abs(M%ANG_IN(i))+norm
       norm=abs(M%ANG_OUT(i))+norm
       norm=abs(M%D_OUT(i))+norm
    enddo
    if(norm>0.0_dp.OR.print_frame) then
       write(mf,*) " THIS IS A CHART THIS IS A CHART THIS IS A CHART THIS IS A CHART "
       CALL print_magnet_frame(m%F,mf)
       WRITE(LINE,*) M%D_IN,M%ANG_IN
       WRITE(MF,'(A255)') LINE
       WRITE(LINE,*) M%D_OUT,M%ANG_OUT
       WRITE(MF,'(A255)') LINE
       write(mf,*) " END OF A CHART  END OF A CHART  END OF A CHART  END OF A CHART  "
    else
       write(mf,*) " NO CHART "
    endif
  end subroutine print_chart

  subroutine READ_chart(m,mf)
    implicit none
    integer mf
    type(CHART), pointer :: m
    character*60 line
    READ(mf,*) LINE
    if(index(line,"NO")==0) then
       CALL READ_magnet_frame(m%F,mf)
       READ(MF,*) M%D_IN,M%ANG_IN
       READ(MF,*) M%D_OUT,M%ANG_OUT
       READ(mf,*) LINE
    else
       do_survey=my_true
    endif
  end subroutine READ_chart


  subroutine READ_chart_fake(mf)
    implicit none
    integer mf
    character*60 line
    type(magnet_frame), pointer :: f
    real(dp) d1(3),d2(3)

    call alloc(f)

    READ(mf,*) LINE
    if(index(line,"NO")==0) then
       CALL READ_magnet_frame(F,mf)
       READ(MF,*) d1,d2
       READ(MF,*) d1,d2
       READ(mf,*) LINE
    else
       do_survey=my_true
    endif
    call kill(f)
  end subroutine READ_chart_fake


  subroutine print_element(P,m,mf)
    implicit none
    integer mf,I
    type(FIBRE), pointer :: P
    type(element), pointer :: m
    character*255 line
    integer f0
    logical perm_
    f0=1

    WRITE(MF,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ELEMENT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    if(m%vorname(1:1)==' ') then
       WRITE(MF,*) M%KIND,M%NAME, ' NOVORNAME'
    ELSE
       WRITE(MF,*) M%KIND,M%NAME,' ',M%VORNAME
    ENDIF
    perm_=.false.
    if(M%p%PERMFRINGE/=0) perm_=.true.
!    WRITE(MF,*) M%L,perm_,M%MIS , " L,PERMFRINGE,MIS "
    WRITE(MF,*) M%L,M%p%PERMFRINGE,M%MIS , " L,PERMFRINGE,MIS "
    WRITE(LINE,*) M%FINT,M%HGAP,M%H1,M%H2,m%va,m%vs, " FINT,HGAP,H1,H2,VA,VS "
    WRITE(MF,'(A255)') LINE
    WRITE(LINE,*) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, " no more mis"
    WRITE(MF,'(A255)') LINE
    IF(ASSOCIATED(M%DELTA_E).and.ASSOCIATED(M%FREQ)) THEN
       WRITE(MF,*) " CAVITY INFORMATION "
       WRITE(LINE,*) M%VOLT, M%FREQ,M%PHAS,M%DELTA_E,M%LAG,M%THIN, " VOLT,FREQ, PHAS, DELTA_E, LAG, THIN"
       WRITE(MF,'(A255)') LINE
    ELSEIF(.not.ASSOCIATED(M%DELTA_E).and.ASSOCIATED(M%FREQ)) THEN
       WRITE(MF,*) " HELICAL DIPOLE INFORMATION "
       WRITE(LINE,*) M%FREQ,M%PHAS, " K_Z, PHAS"
       WRITE(MF,'(A255)') LINE
    ELSEIF(ASSOCIATED(M%VOLT)) THEN
       WRITE(MF,*) " ELECTRIC SEPTUM INFORMATION "
       WRITE(MF,*) M%VOLT,M%PHAS, "VOLT, PHAS(rotation angle) "
    ELSE
       WRITE(MF,*) " NO ELECTRIC ELEMENT INFORMATION "
    ENDIF
    IF(ASSOCIATED(M%B_SOL)) THEN
       WRITE(MF,*)  " SOLENOID_PRESENT ",M%B_SOL, " B_SOL"," SADISTIC " 
    ELSE
       WRITE(MF,*) " NO_SOLENOID_PRESENT ",0.0_dp," SADISTIC" 
    ENDIF
    CALL print_magnet_chart(P,m%P,mf)
    if(p%MAG%KIND==KIND7) then
       f0=p%MAG%t7%f
    endif
    if(p%MAG%KIND==KIND2.and.p%MAG%p%method==2) then
       f0=p%MAG%k2%f
    endif
    if(associated(p%MAG%K16)) then
       if(p%MAG%K16%DRIFTKICK.and.p%MAG%p%method==2)  f0=p%MAG%K16%f
    endif
    if(associated(p%MAG%TP10)) then
       if(p%MAG%TP10%DRIFTKICK.and.p%MAG%p%method==2) f0=p%MAG%TP10%f
    endif
    !    if(f0>0) then
    !     Write(mf,*) f0," Internal Recutting "
    !    endif
    IF(ASSOCIATED(M%an)) THEN
       do i=1,m%p%NMUL
          write(line,*) m%bn(i),m%an(i),f0, "  BN AN %f ",I
          write(mf,'(a255)') line
       enddo
    endif
    call print_specific_element(m,mf)
    WRITE(MF,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   END   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
  end subroutine print_element

  subroutine print_pancake(el,mf)
    implicit none
    type(pancake), pointer :: el
    integer mf
    character*120 filename

    ifield_name=ifield_name+1
    filename(1:8)="fieldmap"
    write(filename(9:120),*) ifield_name
    call context(filename)
    filename=filename(1:len_trim(filename))//'.TXT'
    call context(filename)
    write(mf,*) filename
    call print_pancake_field(el,filename)
  end subroutine print_pancake


  subroutine print_pancake_field(el,filename)
    implicit none
    type(pancake), pointer :: el
    integer mf,nst,i,j
    character(*) filename
    real(dp) brho,cl
    type(real_8) b(3)


    call kanalnummer(mf)
    open(unit=mf,file=filename) !,recl=200)

    nst=2*el%p%nst+1

    cl=(clight/1e8_dp)
    BRHO=el%p%p0c*10.0_dp/cl

    call init(EL%B(1)%no,2)
    CALL ALLOC(B)

    write(mf,*) nst,el%p%ld,el%p%b0,EL%B(1)%no,my_false
    do i=1,nst
       B(1)=morph(1.0_dp.mono.1)
       B(2)=morph(1.0_dp.mono.2)
       B(3)=0.0_dp;
       CALL trackg(EL%B(i),B)
       do j=1,3
          b(j)=b(j)*brho
          call print(b(j),mf)
       enddo

    enddo
    CALL kill(B)

    close(mf)

  end subroutine print_pancake_field




  subroutine print_wig(el,mf)
    implicit none
    type(SAGAN), pointer :: el
    integer mf
    write(mf,*) el%internal
    call print_undu_R(el%w,mf)
  end subroutine print_wig

  subroutine read_wig(el,mf)
    implicit none
    type(SAGAN), pointer :: el
    integer mf
    if(.not.associated(el%internal)) allocate(el%internal(6))
    read(mf,*) el%internal
    call read_undu_R(el%w,mf)
  end subroutine read_wig

  subroutine print_undu_R(el,mf)
    implicit none
    type(undu_R), pointer :: el
    integer mf,n,i
    character*255 line

    write(mf,*) " Undulator internal type undu_R"
    n=size(EL%FORM)

    write(mf,*) n,EL%offset
    do i=1,n
       write(line,*) el%a(i),el%f(i),EL%FORM(i),EL%K(1:3,i)
       write(mf,'(a255)') line
    enddo

    write(mf,*) " End of Undulator internal type undu_R"

  end subroutine print_undu_R

  subroutine read_undu_R(el,mf)
    implicit none
    type(undu_R), pointer :: el
    integer mf,n,i,ne
    character*255 line
    real(dp) offset

    read(mf,'(a255)') line
    read(mf,*) n,offset
    ne=n
    call POINTERS_W(EL,N,ne)
    el%offset=offset
    do i=1,n
       read(mf,*) el%a(i),el%f(i),EL%FORM(i),EL%K(1:3,i)
    enddo

    read(mf,'(a255)') line

  end subroutine read_undu_R


  subroutine print_specific_element(el,mf)
    implicit none
    type(element), pointer :: el
    integer mf,i
    character*255 line

    select case(el%kind)
    CASE(KIND0,KIND1,kind2,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17,KIND22)
    case(kind3)
       WRITE(LINE,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle," patch_edge_ls ",&
            el%k3%patch,el%k3%hf,el%k3%vf,el%k3%ls 
       WRITE(MF,'(A255)') LINE
       WRITE(line,*) el%k3%dx,el%k3%dy,el%k3%pitch_x,el%k3%pitch_y, " SADISTIC" 
       WRITE(MF,'(A255)') LINE
    case(kind4)
       WRITE(line,*) el%c4%N_BESSEL, " HARMON ",el%c4%NF," constant&ripple ",el%c4%a,el%c4%r,el%c4%always_on
       WRITE(MF,'(A255)') LINE
       WRITE(MF,*) el%c4%t,el%c4%phase0,el%c4%CAVITY_TOTALPATH
       do i=1,el%c4%NF
          write(mf,*) el%c4%f(i),el%c4%ph(i)
       enddo
    case(kind5)
       WRITE(line,*) el%s5%dx,el%s5%dy,el%s5%pitch_x,el%s5%pitch_y," SADISTIC" 
       WRITE(MF,'(A255)') LINE
    case(kind10)
       WRITE(MF,*) el%tp10%DRIFTKICK,  " driftkick "
    case(kind16,kind20)
       WRITE(MF,*) el%k16%DRIFTKICK,el%k16%LIKEMAD, " driftkick,likemad"
    case(kind18)
!       WRITE(MF,*) " RCOLLIMATOR HAS AN INTRINSIC APERTURE "
!       CALL print_aperture(EL%RCOL18%A,mf)
    case(kind19)
 !      WRITE(MF,*) " ECOLLIMATOR HAS AN INTRINSIC APERTURE "
 !      CALL print_aperture(EL%ECOL19%A,mf)
    case(kind21)
       WRITE(MF,*) el%cav21%PSI,el%cav21%DPHAS,el%cav21%DVDS,el%cav21%always_on
    case(KINDWIGGLER)
       call print_wig(el%wi,mf)
    case(KINDpa)
       call print_pancake(el%pa,mf)
    case default
       write(6,*) " not supported in print_specific_element",el%kind
       stop 101
    end select

  end subroutine print_specific_element

  subroutine read_specific_element(el,mf,sad)
    implicit none
    type(element), pointer :: el
    integer mf,NB,NH,i,i1
    CHARACTER*6 HARMON
    CHARACTER*15 rip
    character*255 line
    real(dp) x1,x2,x3,x4
    logical(lp) always_on,sad
    select case(el%kind)
    CASE(KIND0,KIND1,kind2,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17,kind22)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
     case(kind5)
       CALL SETFAMILY(EL)
       if(sad) read(mf,*) el%s5%dx,el%s5%dy,el%s5%pitch_x,el%s5%pitch_y
    case(kind3)
       IF(.NOT.ASSOCIATED(el%B_SOL)) then
          ALLOCATE(el%B_SOL)
          el%B_SOL=0.0_dp
       endif
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       READ(MF,"(a255)") LINE
       if(index(line,"patch_edge_ls")/=0) then
          read(line,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle &
               ,rip,el%k3%patch,el%k3%hf,el%k3%vf,el%k3%ls
       elseif(index(line,"patch_edge")/=0) then
          read(line,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle &
               ,rip,el%k3%patch,el%k3%hf,el%k3%vf
       elseif(index(line,"patch")/=0) then
          read(line,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle &
               ,harmon,el%k3%patch
       else
          read(line,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle
       endif
       if(sad) read(mf,*) el%k3%dx,el%k3%dy,el%k3%pitch_x,el%k3%pitch_y
    case(kind4)
       NB=0
       NH=0
       x3=0.0_dp
       x4=1.0_dp
       READ(MF,'(a255)') LINE

       IF(INDEX(LINE, 'HARMON')/=0) THEN
          if(INDEX(LINE, 'ripple')/=0) then
             read(LINE,*) NB,HARMON,NH,rip,x3,x4,always_on
          else
             read(LINE,*) NB,HARMON,NH
          endif
          if(nh>N_CAV4_F) then
             N_CAV4_F=NH
          endif
          read(MF,*) x1,x2,i1
       ELSE
          read(LINE,*) NB
          x1=0.0_dp
          x2=c_%phase0
          i1=CAVITY_TOTALPATH
       ENDIF
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       el%c4%N_BESSEL=NB
       N_CAV4_F=1
       do i=1,nh
          read(mf,*) el%c4%f(i),el%c4%ph(i)
       enddo
       el%c4%t=x1
       el%c4%phase0=x2
       el%c4%a=x3
       el%c4%r=x4
       el%c4%CAVITY_TOTALPATH=i1
       el%c4%always_on=always_on
    case(kind10)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%tp10%DRIFTKICK
    case(kind16,kind20)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%k16%DRIFTKICK,el%k16%LIKEMAD
    case(kind18)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
  !     READ(MF,*) LINE
  !     CALL READ_aperture(EL%RCOL18%A,mf)
    case(kind19)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
    !   READ(MF,*) LINE
    !   CALL READ_aperture(EL%ECOL19%A,mf)
    case(kind21)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%cav21%PSI,el%cav21%DPHAS,el%cav21%DVDS,el%cav21%always_on
    case(KINDWIGGLER)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       call read_wig(el%wi,mf)
    case(KINDpa)
       call read_pancake(el,mf)  ! SET FAMILY DONE INSIDE
    case default
       write(6,*) " not supported in read_specific_element"
       stop 102
    end select

  end subroutine read_specific_element

  subroutine read_pancake(el,mf)
    implicit none
    type(ELEMENT), pointer :: el
    integer mf
    character*120 filename
    read(mf,*) filename
    call context(filename)

    call read_pancake_field(el,filename)
  end subroutine read_pancake


  subroutine read_pancake_field(el,filename)
    implicit none
    type(ELEMENT), pointer :: el
    integer mf,nst,ORDER,I
    real(dp)  L,hc,cl,BRHO
    logical(lp) REPEAT
    character(*) filename
    TYPE(TAYLOR) B(3)
    type(tree_element), allocatable :: t_e(:)

    cl=(clight/1e8_dp)
    BRHO=el%p%p0c*10.0_dp/cl


    call kanalnummer(mf,filename,old=.true.)
  !  open(unit=mf,file=filename)
    read(mf,*) nst,L,hc, ORDER,REPEAT
    CALL INIT(ORDER,2)
    CALL ALLOC(B)
    ALLOCATE(T_E(NST))
    DO I=1,NST
       CALL READ(B(1),mf);CALL READ(B(2),mf);CALL READ(B(3),mf);
       B(1)=B(1)/BRHO
       B(2)=B(2)/BRHO
       B(3)=B(3)/BRHO
       CALL SET_TREE_g(T_E(i),B)
    ENDDO
    close(mf)
    CALL KILL(B)
    CALL SETFAMILY(EL,t=T_E)  !,T_ax=T_ax,T_ay=T_ay)
    deallocate(T_E)

  end subroutine read_pancake_field

  subroutine READ_element(p,m,mf)
    implicit none
    integer mf,I,nmul
    type(fibre), pointer :: p
    type(element), pointer :: m
    character*255 line
    character*255 linet
    CHARACTER*21 SOL
    REAL(DP) B_SOL,r(3),d(3)
    integer f0
    logical perm_,ffl
    logical(lp) sad
    f0=0
    sad=.false.
 !   if(m%kind==kind5.or.m%kind==kind3) sad=.true.
    READ(MF,*) LINE
    READ(MF,*) M%KIND,M%NAME,M%VORNAME
    CALL CONTEXT(M%NAME);
    CALL CONTEXT(M%VORNAME);
    IF(M%VORNAME(1:9)=='NOVORNAME') M%VORNAME=' '

 !   READ(MF,*) M%L,perm_,M%MIS  !,M%EXACTMIS
    READ(MF,'(A120)') LINE
    linet=line
    CALL CONTEXT(LINEt)
ffl=(index(LINEt,"FFL")/=0).or.(index(LINEt,"TFL")/=0).or. &
     (index(LINEt,"FTL")/=0).or.(index(LINEt,"TTL")/=0)
     if(ffl) then
       READ(line,*) M%L,perm_,M%MIS
            if(perm_) then
             M%p%PERMFRINGE=1
            else
             M%p%PERMFRINGE=0     
            endif
       else
        READ(line,*) M%L,M%p%PERMFRINGE,M%MIS     
     endif

    READ(MF,'(a255)') LINE
!write(6,'(a255)') line 
 !   READ(MF,*) M%FINT,M%HGAP,M%H1,M%H2
    if(index(LINE,"VA")==0 ) then
    ! READ(line,*) M%FINT,M%HGAP,M%H1,M%H2
M%FINT=0.5d0
M%HGAP=0
M%H1=0
M%H2=0
     else
     READ(line,*) M%FINT,M%HGAP,M%H1,M%H2,m%VA,m%VS
    endif
    READ(MF,*) R,D
    READ(MF,'(a255)') LINE
    CALL CONTEXT(LINE)
    IF(LINE(1:1)=='C') THEN
       IF(.NOT.ASSOCIATED(M%VOLT)) ALLOCATE(M%VOLT)
       IF(.NOT.ASSOCIATED(M%FREQ)) ALLOCATE(M%FREQ)
       IF(.NOT.ASSOCIATED(M%PHAS)) ALLOCATE(M%PHAS)
       IF(.NOT.ASSOCIATED(M%DELTA_E))ALLOCATE(M%DELTA_E)
       IF(.NOT.ASSOCIATED(M%LAG))   ALLOCATE(M%LAG)
       IF(.NOT.ASSOCIATED(M%THIN))  ALLOCATE(M%THIN)
       READ(MF,*) M%VOLT, M%FREQ,M%PHAS,M%DELTA_E,M%LAG,M%THIN
    ELSEIF(LINE(1:1)=='H') THEN
       IF(.NOT.ASSOCIATED(M%FREQ)) ALLOCATE(M%FREQ)
       IF(.NOT.ASSOCIATED(M%PHAS)) ALLOCATE(M%PHAS)
       READ(MF,*) M%FREQ,M%PHAS
    ELSEIF(LINE(1:1)=='E') THEN
       IF(.NOT.ASSOCIATED(M%VOLT)) ALLOCATE(M%VOLT)
       IF(.NOT.ASSOCIATED(M%PHAS)) ALLOCATE(M%PHAS)
       READ(MF,*) M%VOLT, M%PHAS
    ENDIF
    READ(MF,'(a255)') LINE
     if(index(line,"SAD")/=0) sad=.true.
    READ(LINE,*) SOL,B_SOL
    CALL CONTEXT(SOL)
    IF(SOL(1:2)=='SO') THEN
       IF(.NOT.ASSOCIATED(M%B_SOL))ALLOCATE(M%B_SOL)
       M%B_SOL=B_SOL
    ENDIF
    CALL  READ_magnet_chart(p,m%P,mf)
    nmul = M%P%NMUL
    !     Write(mf,*) f0," Internal Recutting "
    if(m%kind==kind10) then
      M%P%NMUL=sector_nmul_max
    endif
    IF(M%P%NMUL/=0) THEN
       IF(.NOT.ASSOCIATED(M%AN)) THEN
          ALLOCATE(M%AN(M%P%NMUL))
          ALLOCATE(M%BN(M%P%NMUL))
       ELSE
          DEALLOCATE(M%AN)
          DEALLOCATE(M%BN)
          ALLOCATE(M%AN(M%P%NMUL))
          ALLOCATE(M%BN(M%P%NMUL))
       ENDIF
       !     write(6,*) M%KIND,M%NAME,M%VORNAME

       !     write(6,*) M%P%NMUL
       !          READ(MF,'(a120)') LINE
       !     write(6,'(a120)') line
       !     pause 1
       do i=1,nmul   !m%p%NMUL
          READ(MF,'(A255)') LINEt
          if(index(LINEt,"%f")==0 ) then
             READ(linet,*) m%bn(i),m%an(i)
          else
             READ(linet,*) m%bn(i),m%an(i),f0
          endif
          !          READ(mf,*) m%bn(i),m%an(i),f0
       enddo
    endif
    call read_specific_element(m,mf,sad)

    READ(MF,*) LINE

    if(f0>0) then
       if(p%MAG%KIND==KIND7) then
          p%MAG%t7%f=f0
       endif
       if(p%MAG%KIND==KIND2.and.p%MAG%p%method==2) then
          p%MAG%k2%f=f0
       endif
       if(associated(p%MAG%K16)) then
          if(p%MAG%K16%DRIFTKICK.and.p%MAG%p%method==2)  p%MAG%K16%f=f0
       endif
       if(associated(p%MAG%TP10)) then
          if(p%MAG%TP10%DRIFTKICK.and.p%MAG%p%method==2) p%MAG%TP10%f=f0
       endif
    endif

  end subroutine READ_element


  subroutine READ_fake_element(p,mf)
    implicit none
    integer mf,I
    type(fibre),pointer ::  p
    type(fibre),pointer ::  f
    type(element),pointer ::  m
    character*255 line
    CHARACTER*21 SOL
    REAL(DP) B_SOL,r(3),d(3)
    logical(lp) sad

    sad=.false.
 !   if(m%kind==kind5.or.m%kind==kind3) sad=.true.

    nullify(m)
    call alloc(f)
    m=>f%mag
    !    m=0
    READ(MF,'(a255)') LINE
    READ(MF,*) M%KIND,M%NAME,M%VORNAME
    CALL CONTEXT(M%NAME);
    CALL CONTEXT(M%VORNAME);
    IF(M%VORNAME(1:9)=='NOVORNAME') M%VORNAME=' '

    READ(MF,*) M%L,M%p%PERMFRINGE,M%MIS   !,M%EXACTMIS
 !   READ(MF,*) M%FINT,M%HGAP,M%H1,M%H2
    READ(MF,'(a255)') LINE
    if(index(LINE,"VA")==0 ) then
     READ(line,'(a255)') M%FINT,M%HGAP,M%H1,M%H2
     else
     READ(line,'(a255)') M%FINT,M%HGAP,M%H1,M%H2,m%VA,m%VS
    endif
    READ(MF,*) R,D
    READ(MF,'(a255)') LINE
    CALL CONTEXT(LINE)
    IF(LINE(1:1)=='C') THEN
       IF(.NOT.ASSOCIATED(M%VOLT)) ALLOCATE(M%VOLT)
       IF(.NOT.ASSOCIATED(M%FREQ)) ALLOCATE(M%FREQ)
       IF(.NOT.ASSOCIATED(M%PHAS)) ALLOCATE(M%PHAS)
       IF(.NOT.ASSOCIATED(M%DELTA_E))ALLOCATE(M%DELTA_E)
       IF(.NOT.ASSOCIATED(M%LAG))   ALLOCATE(M%LAG)
       IF(.NOT.ASSOCIATED(M%THIN))  ALLOCATE(M%THIN)
       READ(MF,*) M%VOLT, M%FREQ,M%PHAS,M%DELTA_E,M%LAG,M%THIN
    ELSEIF(LINE(1:1)=='E') THEN
       IF(.NOT.ASSOCIATED(M%VOLT)) ALLOCATE(M%VOLT)
       IF(.NOT.ASSOCIATED(M%PHAS)) ALLOCATE(M%PHAS)
       READ(MF,*) M%VOLT, M%PHAS
    ENDIF
    READ(mf,*) SOL,B_SOL
    CALL CONTEXT(SOL)
    IF(SOL(1:2)=='SO') THEN
       IF(.NOT.ASSOCIATED(M%B_SOL))ALLOCATE(M%B_SOL)
       M%B_SOL=B_SOL
    ENDIF
    CALL  READ_magnet_chart(p,m%P,mf)
    IF(M%P%NMUL/=0) THEN
       IF(.NOT.ASSOCIATED(M%AN)) THEN
          ALLOCATE(M%AN(M%P%NMUL))
          ALLOCATE(M%BN(M%P%NMUL))
       ELSE
          DEALLOCATE(M%AN)
          DEALLOCATE(M%BN)
          ALLOCATE(M%AN(M%P%NMUL))
          ALLOCATE(M%BN(M%P%NMUL))
       ENDIF
       !     write(6,*) M%KIND,M%NAME,M%VORNAME

       !     write(6,*) M%P%NMUL
       !          READ(MF,'(a120)') LINE
       !     write(6,'(a120)') line
       !     pause 1
       do i=1,m%p%NMUL
          READ(mf,*) m%bn(i),m%an(i)
       enddo
    endif
    call read_specific_element(m,mf,sad)

    READ(MF,*) LINE

    !          m=-1;
    call SUPER_zero_fibre(f,-1)
    !          deallocate(m);

  end subroutine READ_fake_element

  subroutine print_magnet_chart(P,m,mf)
    implicit none
    type(FIBRE), pointer :: P
    type(magnet_chart), pointer :: m
    integer mf
    character*200 line

    WRITE(MF,*) "MAGNET CHART MAGNET CHART MAGNET CHART MAGNET CHART MAGNET CHART MAGNET CHART "
    WRITE(MF,*) M%EXACT,M%METHOD,M%NST,M%NMUL, " EXACT METHOD NST NMUL"
    WRITE(line,*) M%LD, M%LC, M%B0,' TILT= ',M%TILTD, " LD LC B0 "
    WRITE(MF,'(A200)') LINE
    WRITE(LINE,*) P%BETA0,P%GAMMA0I, P%GAMBET, M%P0C, " BETA0 GAMMA0I GAMBET P0C"
    WRITE(MF,'(A200)') LINE
    WRITE(MF,*) M%EDGE, " EDGES"
    WRITE(MF,*) M%KILL_ENT_FRINGE,M%KILL_EXI_FRINGE,M%bend_fringe, " Kill_ent_fringe, kill_exi_fringe, bend_fringe "

    CALL print_magnet_frame(m%F,mf)
    CALL print_aperture(m%APERTURE,mf)
    write(mf,'(a68)') "END MAGNET CHART END MAGNET CHART END MAGNET CHART END MAGNET CHART "
  end subroutine print_magnet_chart

  subroutine READ_magnet_chart(p,m,mf)
    implicit none
    type(fibre), pointer :: p
    type(magnet_chart), pointer :: m
    integer mf
    character*200 line
    character*5 til
    real(dp) BETA0,GAMMA0I, GAMBET


    READ(MF,*) LINE
    READ(MF,*) M%EXACT,M%METHOD,M%NST,M%NMUL
    READ(MF,'(A200)') LINE

    if(index(line,"TILT=")/=0) then
       READ(LINE,*) M%LD, M%LC, M%B0,til,M%tiltd
    else
       READ(LINE,*) M%LD, M%LC, M%B0
       M%tiltd=0.0_dp
    endif


    READ(MF,*)BETA0,GAMMA0I, GAMBET, M%P0C
    READ(MF,*) M%EDGE
    READ(MF,*) M%KILL_ENT_FRINGE,M%KILL_EXI_FRINGE,M%bend_fringe

    CALL READ_magnet_frame(m%F,mf)
    CALL READ_aperture(m%APERTURE,mf)
    READ(MF,*) LINE
    p%BETA0=beta0
    p%GAMBET=GAMBET
    p%GAMMA0I=GAMMA0I
    !    p%P0C=M%P0C
    !    M%BETA0 =>p%BETA0
    !    M%GAMMA0I => p%GAMMA0I
    !    M%GAMBET => p%GAMBET

  end subroutine READ_magnet_chart

  subroutine print_magnet_frame(m,mf)
    implicit none
    type(magnet_frame), pointer :: m
    integer mf,i
    if(print_frame) then
       write(mf,'(a72)') "MAGNET FRAME MAGNET FRAME MAGNET FRAME MAGNET FRAME MAGNET FRAME MAGNET FRAME "
       WRITE(MF,*) m%a
       do i=1,3
          WRITE(MF,*) m%ent(i,1:3)
       enddo
       WRITE(MF,*) m%o
       do i=1,3
          WRITE(MF,*) m%mid(i,1:3)
       enddo
       WRITE(MF,*) m%b
       do i=1,3
          WRITE(MF,*) m%exi(i,1:3)
       enddo
       write(mf,'(a68)') "END MAGNET FRAME END MAGNET FRAME END MAGNET FRAME END MAGNET FRAME "
    else
       write(mf,'(a72)') " NO MAGNET FRAME NO MAGNET FRAME NO MAGNET FRAME NO MAGNET FRAME NO MAGNET    "
    endif
  end subroutine print_magnet_frame

  subroutine read_magnet_frame(m,mf)
    implicit none
    type(magnet_frame), pointer :: m
    integer mf,i
    character*120 line

    read(mf,'(a120)') line

    if(index(line,"NO")==0) then

       read(MF,*) m%a
       do i=1,3
          read(MF,*) m%ent(i,1:3)
       enddo
       read(MF,*) m%o
       do i=1,3
          read(MF,*) m%mid(i,1:3)
       enddo
       read(MF,*) m%b
       do i=1,3
          read(MF,*) m%exi(i,1:3)
       enddo
       read(mf,'(a120)') line
    else
       do_survey=my_true
    endif
  end subroutine read_magnet_frame

  subroutine print_aperture(m,mf)
    implicit none
    type(MADX_APERTURE), pointer :: m
    integer mf
    character*200 line

    IF(.NOT.ASSOCIATED(M)) THEN
       write(mf,'(a20)') " NO MAGNET APERTURE "
    ELSE
       write(mf,'(a21)') " HAS MAGNET APERTURE "
       WRITE(MF,*) m%KIND   ! 1,2,3,4
       WRITE(MF,*) m%R
       WRITE(line,*)  m%X,m%Y,' SHIFT ',m%dX,m%dY
       WRITE(MF,'(A200)') LINE
       write(mf,'(a23)')  " END OF MAGNET APERTURE"
    ENDIF

  end subroutine print_aperture


  subroutine READ_aperture(m,mf)
    implicit none
    type(MADX_APERTURE), pointer :: m
    integer mf
    character*200 line
    character*5 ch

    !    READ(mf,'(a120)') LINE(1:120)
    READ(mf,'(a120)') LINE

    CALL CONTEXT(LINE)

    IF(LINE(1:2)/='NO') THEN
       IF(.NOT.ASSOCIATED(M)) THEN
          CALL alloc(M)
       ENDIF

       READ(MF,*) m%KIND   ! 1,2,3,4
       READ(MF,*) m%R
       read(mf,'(a200)') line
       if(index(line,"SHIFT")==0) then
          READ(line,*) m%X,m%Y
       else
          READ(line,*) m%X,m%Y,ch,m%dX,m%dY
       endif
       READ(mf,'(a120)') LINE(1:120)
    ENDIF

  end subroutine READ_aperture

!!!!!

  SUBROUTINE change_fibre(p)
    IMPLICIT NONE
    INTEGER MF
    TYPE(FIBRE), POINTER :: P

    CALL KANALNUMMER(MF)

    OPEN(UNIT=MF,FILE='JUNK_CHANGE_FIBRE.TXT')

    CALL print_FIBRE(P,mf)
    P=-1
    REWIND MF
    CALL alloc_fibre( P )
    CALL READ_FIBRE(P,mf)


    CLOSE(MF)
  END SUBROUTINE change_fibre

!!!!!!!!!!  pointed at layout !!!!!!!!!!!!!!
  subroutine read_COMPLEX_SINGLE_STRUCTURE(U,filename,RING,LMAX0)
    implicit none
    character(*) filename
    integer mf,n,i,n_l,J
    type(MAD_UNIVERSE), TARGET :: U
    type(layout), pointer :: L
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0

    call kanalnummer(mf)
    open(unit=mf,file=filename,status='OLD',err=2001)


    read(mf,*) n,n_l
    write(6,*) n,n_l
    do i=1,n
       call READ_AND_APPEND_VIRGIN_LAYOUT(U,filename,RING,LMAX0,mf)
       if(i==1) L=>U%end
       write(6,*) " read layout ", i
       write(6,*) U%end%name
    enddo
    do i=1,n_l
       call APPEND_EMPTY_LAYOUT(U)
       allocate(U%END%DNA(N))
       U%END%DNA(1)%L=>L
       U%END%DNA(1)%counter=0
       DO J=2,N
          U%END%DNA(J)%L=>U%END%DNA(J-1)%L%NEXT
          U%END%DNA(j)%counter=0
       ENDDO

       !   read(mf,'(a120)') line
       WRITE(6,*) "LAYOUT ",I
       !   WRITE(6,*) LINE
       !       do k=1,N
       !        read(mf,*) ncon1,ncon2
       !        if(ncon1/=0) then
       !          allocate(U%END%DNA(k)%L%con1(ncon1))
       !          allocate(U%END%DNA(k)%L%con2(ncon2))
       !          U%END%DNA(k)%L%con1(1:ncon1)%POS=0
       !          U%END%DNA(k)%L%con2(1:ncon2)%POS=0
       !        else
       !          nullify(U%END%DNA(k)%L%con1)
       !          nullify(U%END%DNA(k)%L%con2)
       !       endif

       !          do j=1,max(ncon1,ncon2)
       !            READ(MF,*) POS1,POS2
       !           IF(POS1/=0) U%END%DNA(k)%L%con1(J)%POS=POS1
       !           IF(POS2/=0) U%END%DNA(k)%L%con2(J)%POS=POS2
       !         enddo

       !       enddo
       !   read(mf,'(a120)') line

       !   WRITE(6,*) I
       !   WRITE(6,*) LINE


       CALL READ_pointed_INTO_VIRGIN_LAYOUT(U%END,FILENAME,RING,LMAX0,mf1=MF)

       !     do k=1,N
       !      DO J=1,SIZE(U%END%DNA(k)%L%con1)
       !        CALL MOVE_TO(U%END,P,U%END%DNA(k)%L%CON1(J)%POS)
       !        U%END%DNA(k)%L%CON1(J)%P=>P
       !     ENDDO
       !     DO J=1,SIZE(U%END%DNA(k)%L%con2)
       !       CALL MOVE_TO(U%END,P,U%END%DNA(k)%L%CON2(J)%POS)
       !       U%END%DNA(k)%L%CON2(J)%P=>P
       !     ENDDO
       !    ENDDO

    enddo

    close(mf)



    return
2001 continue

    Write(6,*) " File ",filename(1:len_trim(filename)) ," does not exist "

  END subroutine read_COMPLEX_SINGLE_STRUCTURE

  ! MAKES NOT SENSE BECAUSE DNA ARRAY NOT PROVIDED!
  !  subroutine READ_pointed_AND_APPEND_VIRGIN_LAYOUT(U,filename,RING,mf)
  !    implicit none
  !    character(*) filename
  !    integer,optional :: mf
  !    type(MAD_UNIVERSE), TARGET :: U
  !    LOGICAL(LP), OPTIONAL :: RING


  !    call APPEND_EMPTY_LAYOUT(U)

  !    CALL READ_pointed_INTO_VIRGIN_LAYOUT(U%END,FILENAME,RING,mf)

  !  END subroutine READ_pointed_AND_APPEND_VIRGIN_LAYOUT

  subroutine READ_pointed_INTO_VIRGIN_LAYOUT(L,FILENAME,RING,LMAX0,mf1)
    implicit none
    character(*) filename
    integer I,mf,N,DIR,index_,pos,nt,siam_index,siam_pos,gird_index,gird_pos
    type(LAYOUT), TARGET :: L
    type(FIBRE), pointer :: P,current,siam,gird
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0
    LOGICAL(LP) RING_IT,doneit
    character*120 line
    character*255 lineg
    real(dp) p0c,MASSF,LMAX0t,ag0
    type(internal_state) original
    integer,optional :: mf1
    integer res
    integer se1,se2

    RING_IT=MY_TRUE

    IF(PRESENT(RING)) RING_IT=RING

    if(present(mf1) ) then
       mf=mf1
    else
       call kanalnummer(mf)
       open(unit=mf,file=filename,status='OLD',err=2001)
    endif

    READ(MF,*) N,LMAX0t
    write(6,*) N,LMAX0t
    IF(PRESENT(LMAX0)) then
       if(LMAX0t/=0.0_dp) LMAX0=LMAX0T
    ENDIF
    read(MF,'(a120)') line
    call context(line)
    write(6,*) line

    if(index(line,"FOR")/=0) then
       l%name=line(index(line,"FOR")+3:index(line,"FOR")+2+nlp)
    endif

    read(MF,'(A255)') lineg
    res=INDEX (lineG, "AG(spin)")
    IF(RES==0) THEN
       read(lineg,*) MASSF,p0c
       IF(ABS(MASSF-pmap)/PMAP<0.01E0_DP) THEN
          A_PARTICLE=A_PROTON
       ELSEIF(ABS(MASSF-pmae)/pmae<0.01E0_DP) THEN
          A_PARTICLE=A_ELECTRON
       ELSEIF(ABS(MASSF-pmaMUON)/pmaMUON<0.01E0_DP) THEN
          A_PARTICLE=A_MUON
       ENDIF
    ELSE
       read(lineg,*) MASSF,p0c,A_PARTICLE
    ENDIF
    ag0=A_PARTICLE

    !    read(MF,*) MASSF,p0c
    read(MF,*) phase0,compute_stoch_kick,initial_charge
    read(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING
    read(MF,*) se1,se2,OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE
    read(MF,'(A255)') lineg
    res=INDEX (lineG, "Valishev")
    IF(RES==0) THEN
       read(lineg,*) wedge_coeff
    ELSE
       read(lineg,*) wedge_coeff,valishev
    ENDIF
    read(MF,*) MAD8_WEDGE
    read(MF,'(a120)') line
    original=default
    call input_sector(se2,se1)
!    if(allocated(s_b)) then
       firsttime_coef=my_true
!       deallocate(s_b)
!    endif
    !    L%MASS=MASSF
    MASSF=MASSF/pmae
    CALL MAKE_STATES(MASSF)
    A_PARTICLE=ag0
    default=original
    call Set_madx(p0c=p0c)
    DO I=1,N
       call  READ_FIBRE_2_lines(mf,dir,index_,pos,nt,siam_index,siam_pos,gird_index,gird_pos)
       call  READ_chart_fake(mf)
       call move_to(L%DNA(index_)%L,p,pos)
       CALL APPEND_POINT(L,P)
       current=>l%end
       current%dir=dir

       if(siam_index/=0) then
          call move_to(L%DNA(siam_index)%L,siam,siam_pos)
          p%mag%siamese=>siam%mag
          write(6,*) p%mag%name,' is a siamese of ', siam%mag%name
       endif
       if(gird_index/=0) then
          call move_to(L%DNA(gird_index)%L,gird,gird_pos)
          p%mag%girderS=>gird%mag
          write(6,*) p%mag%name,' is on the girder of ', gird%mag%name
       endif
       !        if(pos==1) then
       !         if(associated(L%DNA(index)%l%con1)) then
       !           L%DNA(index)%counter=L%DNA(index)%counter+1
       !           kc=L%DNA(index)%counter
       !           L%DNA(index)%l%con1(kc)%p=>current
       !      write(6,*)  0,index,kc
       !          endif
       !        elseif(pos==nt) then
       !         if(associated(L%DNA(index)%l%con2)) then
       !           L%DNA(index)%l%con2(kc)%p=>current
       !     !       write(6,*)  1,index,kc
       !          endif
       !        endif
       CALL READ_PATCH(current%PATCH,mf)
       call READ_fake_element(current,mf)
       READ(MF,*) LINE
       !
       ! CALL READ_FIBRE(L%END,mf)
       !       CALL COPY(L%END%MAG,L%END%MAGP)
    ENDDO

    if(.not.present(mf1) ) CLOSE(MF)

    L%closed=RING_IT

    doneit=my_true
    call ring_l(L,doneit)
    if(do_survey) call survey(L)

    return

2001 continue

    Write(6,*) " File ",filename(1:len_trim(filename)) ," does not exist "

  END subroutine READ_pointed_INTO_VIRGIN_LAYOUT


!!!!!!!!!!!!  switching routines !!!!!!!!!!!!!
  SUBROUTINE switch_layout_to_cavity( L,name,sext,a,r,freq,t )  ! switch to kind7
    implicit none
    TYPE (layout), target :: L
    TYPE (fibre), pointer :: p
    character(*) name
    real(dp),OPTIONAL:: a,r,freq,t
    INTEGER I
    logical(lp) sext


    p=>L%start
    do i=1,L%n

       call  switch_to_cavity( p,name,sext,a,r,freq,t)

       p=>p%next
    enddo



  end SUBROUTINE switch_layout_to_cavity



  SUBROUTINE switch_to_cavity( el,name,sext,a,r,freq,t )  ! switch to kind7
    implicit none
    TYPE (fibre), target :: el
    character(*) name
    integer i,nm,EXCEPTION
    real(dp),OPTIONAL:: a,r,freq,t
    real(dp), allocatable :: an(:),bn(:)
    type(keywords) key
    logical(lp) sext
    ! This routines switches to cavity
    nm=len_trim(name)
    select case(el%mag%kind)
    case(kind10,kind16,kind2,kind7,kind6,KIND20)
       if(el%mag%name(1:nm)==name(1:nm)) then
          if(sext.and.el%mag%p%nmul>2)then
             write(6,*) el%mag%name
             call add(el,3,0,0.0_dp)
             call add(el,-3,0,0.0_dp)
          else
             write(6,*) el%mag%name, " not changed "
          endif
       endif
       if(el%mag%p%b0/=0.0_dp.or.el%mag%name(1:nm)==name(1:nm)) return
       write(6,*) el%mag%name
       if(el%mag%p%nmul/=size(el%mag%an)) then
          write(6,*) "error in switch_to_cavity "
          stop 666
       endif
       allocate(an(el%mag%p%nmul),bn(el%mag%p%nmul))
       an=el%mag%an*el%mag%p%p0c
       bn=el%mag%bn*el%mag%p%p0c
       call zero_key(key)
       key%magnet="rfcavity"
       key%list%volt=0.0_dp
       key%list%lag=0.0_dp
       key%list%freq0=0.0_dp
       IF(PRESENT(FREQ)) THEN
          key%list%freq0=FREQ
       ENDIF
       key%list%n_bessel=0
       key%list%harmon=1.0_dp
       key%list%l=el%mag%L
       key%list%name=el%mag%name
       key%list%vorname=el%mag%vorname
       EXACT_MODEL=el%mag%p%exact
       key%nstep=el%mag%p%nst
       key%method=el%mag%p%method


       el%mag=-1
       el%magp=-1
       el%mag=0
       el%magp=0
       call create_fibre(el,key,EXCEPTION,my_true)
       do i=size(an),1,-1
          call add(el,-i,0,an(i))
          call add(el,i,0,bn(i))
       enddo
       el%mag%c4%a=1.0_dp
       el%magp%c4%a=1.0_dp
       IF(PRESENT(a)) THEN
          el%mag%c4%a=a
          el%magp%c4%a=a
       ENDIF
       el%mag%c4%r=0.0_dp
       el%magp%c4%r=0.0_dp
       IF(PRESENT(r)) THEN
          el%mag%c4%r=r
          el%magp%c4%r=r
       ENDIF
       el%mag%c4%PHASE0=0.0_dp
       el%mag%c4%PHASE0=0.0_dp
       el%mag%c4%always_on=my_true
       el%magp%c4%always_on=my_true
       IF(PRESENT(FREQ)) THEN
          el%mag%FREQ=FREQ
          el%magP%FREQ=FREQ
       ENDIF
       IF(PRESENT(T).and.PRESENT(FREQ)) THEN
          el%mag%C4%T=T/(el%mag%C4%freq/CLIGHT)
          el%magP%C4%T=el%mag%C4%T
       ENDIF
       deallocate(an,bn)
    CASE(KIND4)
       IF(el%mag%c4%always_on) THEN
          IF(PRESENT(FREQ)) THEN
             el%mag%FREQ=FREQ
             el%magP%FREQ=FREQ
          ENDIF
          IF(PRESENT(T)) THEN
             el%mag%C4%T=T/(el%mag%C4%freq/CLIGHT)
             el%magP%C4%T=el%mag%C4%T
          ENDIF
          IF(PRESENT(r)) THEN
             el%mag%c4%r=r
             el%magp%c4%r=r
          ENDIF
          IF(PRESENT(a)) THEN
             el%mag%c4%a=a
             el%magp%c4%a=a
          ENDIF

       ENDIF

    case default
       return
    end select

  END SUBROUTINE switch_to_cavity

  SUBROUTINE switch_to_kind7( el )  ! switch to kind7
    implicit none
    TYPE (fibre), target :: el
    ! This routines switches to kind7 (not exact) from kind2,10,16
    select case(el%mag%kind)
    case(kind10,kind16,kind2,KIND20)
       el%magp=-1
       el%mag%L=el%mag%p%ld
       el%mag%p%lc=el%mag%p%ld
       el%mag%p%exact=my_false
       el%magp=0
    end select

    select case(el%mag%kind)
    case(kind10)
       EL%MAG%TP10=-1
       deallocate(EL%MAG%TP10)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    case(kind16,KIND20)
       EL%MAG%k16=-1
       deallocate(EL%MAG%k16)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    case(KIND2)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    end select

  END SUBROUTINE switch_to_kind7

!!!!!! New flat file for data base universe  M_T !!!!!!!!!!!!!!!!! 
  
subroutine  print_new_flat(ring,filename,last,com)

implicit none
type(layout), target :: ring
type(fibre), pointer :: f
logical(lp), optional ::last
character(6), optional ::com
character(*) filename
integer i,mf
character(120) line
logical(lp) fin
character (6) comt
comt='REWIND'
fin=my_true
if(present(last)) fin=last
!goto 1
if(present(com)) comt=com
call kanalnummer(mf)
open(unit=mf,file=filename,position=comt) !,recl=200) !comt could be append for a complex universe 

   write(mf,'(a120)') ring%name                        ! Sagan depedent line
   write(mf,*) highest_fringe  , " highest fringe "    !  Sagan depedent line DO NOT CHANGE
   write(mf,*) lmax, " Maximum Length for Orbit "
   write(MF,*) ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING  , "ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING "
   write(mf,*) SECTOR_NMUL,SECTOR_NMUL_MAX , " SECTOR_NMUL,SECTOR_NMUL_MAX "
    
 write(mf,*) " $$$$$$$$$$$$$$$$$ START OF LAYOUT $$$$$$$$$$$$$$$$$"


  
f=>ring%start

call Print_initial_chart(f,mf)

do i=1,ring%n
  call el_el0(f%mag,my_true,mf)
  call fib_fib0(f,my_true,mf)
  CALL MC_MC0(f%MAG%P,my_true,mf)
  CALL print_ElementLIST(f%mag,MY_TRUE,mf)
  if(f%patch%patch/=0.or.f%patch%time/=0.or.f%patch%energy/=0) call patch_patch0(f%patch,my_true,mf)
  if(f%mag%mis) call CHART_CHART0(f%chart,my_true,mf)
 write(mf,*) " $$$$$$$$$$$$$$$$$ END OF FIBRE $$$$$$$$$$$$$$$$$"
 f=>f%next    
enddo

 write(mf,*) "&ELENAME"
if(fin) then
 write(mf,*)  "ELE0%NAME_VORNAME='alldone','alldone',"
else
 write(mf,*)  "ELE0%NAME_VORNAME='endhere','endhere',"
endif 
write(mf,*) "/"


close(mf)


end subroutine print_new_flat

subroutine  print_universe(un,filename)
!call print_universe(m_u,'junk2.txt')
!call print_universe_pointed(m_u,m_t,'junk3.txt')
implicit none
type(mad_universe),target :: un
type(layout),pointer :: r
character(*) filename
integer i



r=>un%start


if(associated(r,un%end)) then
 call  print_new_flat(r,filename,last=my_true)
else
          do while(.not.associated(r,un%end))
            if(associated(r,un%start)) then
              call  print_new_flat(r,filename,last=my_false)  
            else
              call  print_new_flat(r,filename,last=my_false,com="APPEND")  
            endif
            r=>r%next
          enddo
              call  print_new_flat(r,filename,last=my_true,com="APPEND")  
 
endif
call print_universe_siamese(un,filename)
call print_universe_girders(un,filename)
end subroutine  print_universe

subroutine  print_universe_girders(un,filename)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,k,i1,i2,j1,j2,MF



call TIE_MAD_UNIVERSE(un)

p=>un%start%start
p0=>p
p=>p%next


k=0
do while(.not.associated(p0,p))


if(associated(p%mag%girders)) then
 
 if(associated(p%mag%GIRDER_FRAME)) then
 k=k+1

endif
endif
 
 p=>p%n
enddo

call kanalnummer(mf)
open(unit=mf,file=filename,position='APPEND')
write(MF,*) k, " girders in the universe "




p=>un%start%start
p0=>p
p=>p%next


i=0
do while(.not.associated(p0,p))
i=i+1

if(associated(p%mag%girders)) then
 
 
 if(associated(p%mag%GIRDER_FRAME)) then
  j1=0
  j2=0
  ps=>p
  call locate_in_universe(ps,i1,i2)
     write(MF,*) p%mag%girder_frame%a    
     write(MF,*) p%mag%girder_frame%ent(1,1:3)
     write(MF,*) p%mag%girder_frame%ent(2,1:3)
     write(MF,*) p%mag%girder_frame%ent(3,1:3)
     write(MF,*) p%mag%girder_frame%b    
     write(MF,*) p%mag%girder_frame%exi(1,1:3)
     write(MF,*) p%mag%girder_frame%exi(2,1:3)
     write(MF,*) p%mag%girder_frame%exi(3,1:3)



     write(MF,*) i1,i2,ps%loc,PS%MAG%NAME
  do while(j1/=i1.or.j2/=i2)
      ps=>ps%mag%girders%parent_fibre    
      call locate_in_universe(ps,j1,j2)
      if (j1/=i1.or.j2/=i2) then
        write(MF,*) j1,j2 ,ps%loc
      else
        write(MF,*) 0,0,0
      endif
  enddo
 endif

endif
 
 p=>p%n
enddo


CLOSE(MF)

end subroutine  print_universe_girders

subroutine  read_universe_girders(un,mf,ns)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
integer i,k,j,i1,i2,j1,MF,ns


call TIE_MAD_UNIVERSE(un)




do i=1,ns

call read_initial_chart(mf)

read(mf,*) i1,i2,j1

 p0=>un%start%start 
 do j=2,j1
  p0=>p0%n
 enddo
 ps=>p0

 do k=1,1000000 
  read(mf,*) i1,i2,j1
  if(i1==0) exit
  p=>un%start%start
  do j=2,j1
   p=>p%n
  enddo
  ps%mag%girders=>p%mag
  ps=>p
 enddo
  ps%mag%girders=>p0%mag

  allocate(p0%mag%girder_frame)
  call NULL_af(p0%mag%girder_frame)
  allocate(p0%mag%girder_frame%a(3))
  allocate(p0%mag%girder_frame%b(3))
  allocate(p0%mag%girder_frame%ent(3,3))
  allocate(p0%mag%girder_frame%exi(3,3))
  p0%mag%girder_frame%a=a_
  p0%mag%girder_frame%ent(1,1:3)=ent_(1,1:3)
  p0%mag%girder_frame%ent(2,1:3)=ent_(2,1:3)
  p0%mag%girder_frame%ent(3,1:3)=ent_(3,1:3)
  p0%mag%girder_frame%b=b_
  p0%mag%girder_frame%exi(1,1:3)=exi_(1,1:3)
  p0%mag%girder_frame%exi(2,1:3)=exi_(2,1:3)
  p0%mag%girder_frame%exi(3,1:3)=exi_(3,1:3)

enddo


end subroutine  read_universe_girders

subroutine  print_universe_siamese(un,filename)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,k,i1,i2,j1,j2,MF



call TIE_MAD_UNIVERSE(un)

p=>un%start%start
p0=>p
p=>p%next


k=0
do while(.not.associated(p0,p))


if(associated(p%mag%siamese)) then
 
 if(associated(p%mag%SIAMESE_FRAME)) then
 k=k+1

endif
endif
 
 p=>p%n
enddo

call kanalnummer(mf)
open(unit=mf,file=filename,position='APPEND') !,recl=200)
write(MF,*) k, " siamese in the universe "




p=>un%start%start
p0=>p
p=>p%next


i=0
do while(.not.associated(p0,p))
i=i+1

if(associated(p%mag%siamese)) then
 
 
 if(associated(p%mag%SIAMESE_FRAME)) then
  j1=0
  j2=0
  ps=>p
  call locate_in_universe(ps,i1,i2)
     write(MF,*) p%mag%SIAMESE_FRAME%ANGLE
     write(MF,*) p%mag%SIAMESE_FRAME%d
     write(MF,*) i1,i2,ps%loc,PS%MAG%NAME
  do while(j1/=i1.or.j2/=i2)
      ps=>ps%mag%siamese%parent_fibre    
      call locate_in_universe(ps,j1,j2)
      if (j1/=i1.or.j2/=i2) then
        write(MF,*) j1,j2 ,ps%loc
      else
        write(MF,*) 0,0,0
      endif
  enddo
 endif

endif
 
 p=>p%n
enddo


CLOSE(MF)

end subroutine  print_universe_siamese

subroutine  read_universe_siamese(un,mf,ns)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
integer i,k,j,i1,i2,j1,MF,ns
real(dp) a(3),d(3)


call TIE_MAD_UNIVERSE(un)




do i=1,ns

read(mf,*) a
read(mf,*) d

read(mf,*) i1,i2,j1

 p0=>un%start%start 
 do j=2,j1
  p0=>p0%n
 enddo
 ps=>p0

 do k=1,1000000 
  read(mf,*) i1,i2,j1
  if(i1==0) exit
  p=>un%start%start
  do j=2,j1
   p=>p%n
  enddo
  ps%mag%siamese=>p%mag
  ps=>p
 enddo
  ps%mag%siamese=>p0%mag

  allocate(p0%mag%siamese_frame)
  call NULL_af(p0%mag%siamese_frame)
  allocate(p0%mag%siamese_frame%angle(3))
  allocate(p0%mag%siamese_frame%d(3))
  p0%mag%siamese_frame%angle=a
  p0%mag%siamese_frame%d=d
enddo


end subroutine  read_universe_siamese

subroutine  read_universe_database(un,filename,arpent)
! the universes should be empty
!call read_universe_database(m_u,'junk2.txt',arpent=my_false)
!call read_universe_pointed(M_u,M_t,'junk3.txt')
!call create_dna(M_u,m_t)
!arpent= false => the database should not be surveyed.
! DNA is automatically created in create_dna
implicit none
type(mad_universe),target :: un
logical(lp), optional :: arpent
character(*) filename
integer mf,ns
ELE0%NAME_VORNAME(1)=' '
ELE0%NAME_VORNAME(2)=' '

 call kanalnummer(mf,filename(1:len_trim(filename)))
          do while(ELE0%NAME_VORNAME(1)/="alldone")
           call append_empty_layout(un)  
           call set_up(un%end)
           call read_lattice(un%end,filename,mf,arpent)

          enddo
read(mf,*) ns   ! number of siamese
 call read_universe_siamese(un,mf,ns)
read(mf,*) ns   ! number of girders
call read_universe_girders(un,mf,ns)
close(mf)
end subroutine  read_universe_database

subroutine  read_lattice_append(un,filename,arpent)
! the universes should be empty
!call read_universe_database(m_u,'junk2.txt',arpent=my_false)
!call read_universe_pointed(M_u,M_t,'junk3.txt')
!call create_dna(M_u,m_t)
!arpent= false => the databaseshould not be surveyed.
! DNA is automatically created in create_dna
implicit none
type(mad_universe),target :: un
logical(lp), optional :: arpent
character(*) filename
integer mf
ELE0%NAME_VORNAME(1)=' '
ELE0%NAME_VORNAME(2)=' '
 call kanalnummer(mf,filename(1:len_trim(filename)),old=.true.)

           call append_empty_layout(un)  
           call set_up(un%end)
           call read_lattice(un%end,filename,mf,arpent)
close(mf)
end subroutine  read_lattice_append



subroutine  read_lattice(r,filename,mfile,arpent)
implicit none
type(layout),target :: r
character(*) filename
logical(lp), optional :: arpent
integer , optional :: mfile
logical(lp) doneit,surv
character(120) line
type(fibre),pointer :: s22
type(element),pointer :: s2
type(elementp), pointer :: s2p
integer se2,se1

integer mf,n


if(present(mfile)) then
 mf=mfile
else
 call kanalnummer(mf,filename(1:len_trim(filename)),old=.true.)
endif
surv=my_true





!-----------------------------------


   read(mf,'(a120)') r%name
   read(mf,*) highest_fringe  
   read(mf,*) lmax  
   read(MF,*) ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING  
   read(mf,*) se2,se1
     call input_sector(se2,se1)
call make_states(my_false)
call set_mad(energy=2.0d0)


   
 read(mf,'(a120)') line


 
 read(mf,'(a120)') line
call read_initial_chart(mf)
 read(mf,'(a120)') line



n=0
do while(.true.) 
    call zero_ele0
   read(mf,NML=ELEname,end=999) !!  basic stuff : an,bn,L, b_sol, etc...
!write(6,*) ELE0%name_vorname


   if(ELE0%NAME_VORNAME(1)== "endhere".or.ELE0%NAME_VORNAME(1)=="alldone") goto 99
 !write(6,NML=ELEname)
   call zero_fib0
   read(mf,NML=FIBRENAME,end=999) !
!     real(dp) GAMMA0I_GAMBET_MASS_AG(4)  !GAMMA0I,GAMBET,MASS ,AG  BETA0 is computed
 !    real(dp)  CHARGE
 !    integer  DIR !DIR,CHARGE
 !    integer patch 
   call zero_MAGL0
   read(mf,NML=MAGLNAME,end=999) ! reads magnet frame: kill_fringe, nst, method, etc... 
 !write(6,NML=MAGLNAME)
 call read_ElementLIST(ELE0%kind,MF)  ! read individual elments and aperture at the end

                   
 if(fib0%patch/=0) then
  call zero_patch0
  read(mf,NML=patchname,end=999)    ! patch read if present
 endif

if(ele0%recut_even_electric_MIS(4)) then
 call zero_CHART0
 read(mf,NML=CHARTname)  ! reading misalignment
endif


 read(mf,'(a120)') line

 call append_empty(r)  ! appends a fibre
 s22=>r%end
  call  nullify_for_madx(s22)


  

    s2=>s22%mag;
    s2p=>s22%magp;



    call fib_fib0(s22,my_false)
IF(ELE0%kind==KIND10) THEN
 IF(MAGL0%METHOD_NST_NMUL_permfringe_highest(3)>SECTOR_NMUL_MAX) MAGL0%METHOD_NST_NMUL_permfringe_highest(3)=SECTOR_NMUL_MAX
ENDIF
     S2 = MAGL0%METHOD_NST_NMUL_permfringe_highest(3)    
     


!write(6,*) associated(s2%p)
 !pause 78
    call MC_MC0(s2%p,my_false)



 !pause 79
    call el_el0(s2,my_false)



    if(s2%kind/=kindpa) then
       CALL SETFAMILY(S2) 
    else
       CALL SETFAMILY(S2,t=T_E)  !,T_ax=T_ax,T_ay=T_ay)
       S2%P%METHOD=4
       deallocate(T_E,t_ax,t_ay)
    endif    

 

    call print_ElementLIST(s2,my_false)

    s2p=0   
 
 !pause 665

    call copy(s2,s2p)
    s22%mag%parent_fibre =>s22
    s22%magp%parent_fibre=>s22
       s22%CHART=0
       s22%PATCH=0
     if(fib0%patch/=0) then
       !s22%PATCH%patch=fib0%patch
      call patch_patch0(s22%patch,my_false)
    endif
   if(ele0%recut_even_electric_MIS(4)) call CHART_CHART0(s22%chart,my_false)





n=n+1
enddo


   100 CONTINUE


99 write(6,*) ELE0%NAME_VORNAME(1)
999 write(6,*) "Read ",n

 s22=>r%start
! if(s22%dir==1) then
  s22%chart%f%a=a_
  s22%chart%f%ent=ent_
! else
  s22%chart%f%b=b_
  s22%chart%f%exi=exi_
! endif
    r%closed=.true.

    doneit=.true.
    call ring_l(r,doneit)
 if(present(arpent)) surv=arpent

   if(surv) then


     call survey(r)
   endif

1000 continue

if(.not.present(mfile)) close(mf)

end subroutine read_lattice

  subroutine read_ElementLIST(kind,mf)
    implicit none

    integer mf,i,kind
    LOGICAL dir
    character*255 line


    select case(kind)
    CASE(KIND0,KIND1,kind2,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17)
  case(kind3)
     read(mf,NML=thin30name)
    case(kind4)
     read(mf,NML=CAVname)
    case(kind5)
     read(mf,NML=sol50name)
    case(kind10)
      read(mf,NML=tp100name)
    case(kindabell)
      read(mf,NML=ab0name)
    case(kind16,kind20)

     read(mf,NML=k160name)

    case(kind18)

    case(kind19)


    case(kindhel)
      read(mf,NML=helname)
    case(kind21)
     read(mf,NML=tCAVname)
    case(KINDWIGGLER)
      read(mf,NML=wigname)
    case(KINDpa)
 
   case default
       write(MF,*) " not supported in print_specific_element",kind
 !      stop 101
    end select
    
    CALL  r_ap_aplist(mf) 


  END SUBROUTINE read_ElementLIST

subroutine  fib_fib0(f,dir,mf)
implicit none
type(fibre), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
! fib0%t(1)=f%BETA0
 fib0%GAMMA0I_GAMBET_MASS_AG(1)=f%GAMMA0I
 fib0%GAMMA0I_GAMBET_MASS_AG(2)=f%GAMBET
 fib0%GAMMA0I_GAMBET_MASS_AG(3)=f%MASS
 fib0%GAMMA0I_GAMBET_MASS_AG(4)=f%AG
 fib0%DIR=f%DIR
 fib0%CHARGE=f%CHARGE
 fib0%patch=f%patch%patch+7*f%patch%energy+49*f%patch%time
 
    if(present(mf)) then
     write(mf,NML=fibrename)
    endif   
else
    if(present(mf)) then
     read(mf,NML=fibrename)
    endif   
    ! f%BETA0=fib0%t(1)
 f%GAMMA0I=fib0%GAMMA0I_GAMBET_MASS_AG(1)
 f%GAMBET=fib0%GAMMA0I_GAMBET_MASS_AG(2)
 f%MASS=fib0%GAMMA0I_GAMBET_MASS_AG(3)
 f%AG=fib0%GAMMA0I_GAMBET_MASS_AG(4)
 f%BETA0=sqrt(1.0_dp-f%GAMMA0I**2)   
 f%DIR=fib0%DIR 
 f%CHARGE=fib0%CHARGE
 !f%patch%patch=fib0%patch     ! f%patch%patch is not yet allocated
endif
endif

end subroutine fib_fib0

subroutine  patch_patch0(f,dir,mf)
implicit none
type(PATCH), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
! fib0%t(1)=f%BETA0
 
 patch0%A_X1=f%A_X1
 patch0%A_X2=f%A_X2
 patch0%B_X1=f%B_X1
 patch0%B_X2=f%B_X2
 patch0%A_D=f%A_D
 patch0%B_D=f%B_D
 patch0%A_ANG=f%A_ANG
 patch0%B_ANG=f%B_ANG
 patch0%A_T=f%A_T
 patch0%B_T=f%B_T
 patch0%A_L=f%A_L
 patch0%B_L=f%B_L
 patch0%ENERGY=f%ENERGY
 patch0%TIME=f%TIME
 patch0%geometry=f%patch

    if(present(mf)) then
     write(mf,NML=patchname)
    endif   
else
    if(present(mf)) then
     read(mf,NML=patchname)
    endif   

f%A_X1= patch0%A_X1
f%A_X2= patch0%A_X2
f%B_X1= patch0%B_X1
f%B_X2= patch0%B_X2

 f%A_D=patch0%A_D
 f%B_D=patch0%B_D
 f%A_ANG=patch0%A_ANG
 f%B_ANG=patch0%B_ANG
 f%A_T=patch0%A_T
 f%B_T=patch0%B_T
 f%A_L=patch0%A_L
 f%B_L=patch0%B_L
 f%ENERGY=patch0%ENERGY
 f%TIME=patch0%TIME
 f%patch=patch0%geometry

endif
endif
end subroutine patch_patch0


subroutine  CHART_CHART0(f,dir,mf)
implicit none
type(chart), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then    

 CHART0%D_IN=f%D_IN
 CHART0%D_OUT=f%D_OUT
 CHART0%ANG_IN=f%ANG_IN
 CHART0%ANG_OUT=f%ANG_OUT
 
    if(present(mf)) then
     write(mf,NML=CHARTname)
    endif   
else
    if(present(mf)) then
     read(mf,NML=CHARTname)
    endif   

 f%D_IN=CHART0%D_IN
 f%D_OUT=CHART0%D_OUT
 f%ANG_IN=CHART0%ANG_IN
 f%ANG_OUT=CHART0%ANG_OUT
 

endif
endif
end subroutine CHART_CHART0



subroutine  MC_MC0(f,dir,mf)
implicit none
type(MAGNET_CHART), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 MAGL0%LC_LD_B0_P0(1)=f%LC
 MAGL0%LC_LD_B0_P0(2)=f%LD
 MAGL0%LC_LD_B0_P0(3)=f%B0
 MAGL0%LC_LD_B0_P0(4)=f%P0C
 
 MAGL0%TILTD_EDGE(1)=f%TILTD
 MAGL0%TILTD_EDGE(2)=f%EDGE(1)
 MAGL0%TILTD_EDGE(3)=f%EDGE(2)

 MAGL0%KILL_SPIN(1)=f%KILL_ENT_SPIN
 MAGL0%KILL_SPIN(2)=f%KILL_EXI_SPIN

 MAGL0%KIN_KEX_BENDFRINGE_EXACT(1)=f%KILL_ENT_FRINGE
 MAGL0%KIN_KEX_BENDFRINGE_EXACT(2)=f%KILL_EXI_FRINGE
 MAGL0%KIN_KEX_BENDFRINGE_EXACT(3)=f%bend_fringe
 MAGL0%KIN_KEX_BENDFRINGE_EXACT(4)=f%EXACT

 MAGL0%METHOD_NST_NMUL_permfringe_highest(1)=f%METHOD
 MAGL0%METHOD_NST_NMUL_permfringe_highest(2)=f%NST
 MAGL0%METHOD_NST_NMUL_permfringe_highest(3)=f%NMUL
 MAGL0%METHOD_NST_NMUL_permfringe_highest(4)=f%permfringe
 MAGL0%METHOD_NST_NMUL_permfringe_highest(5)=f%highest_fringe

 if(present(mf)) then
     write(mf,NML=MAGLname)
    endif   
else
    if(present(mf)) then
     read(mf,NML=MAGLname)
    endif   
 f%LC=MAGL0%LC_LD_B0_P0(1)
 f%LD=MAGL0%LC_LD_B0_P0(2)
 f%B0=MAGL0%LC_LD_B0_P0(3)
 f%P0C=MAGL0%LC_LD_B0_P0(4)
 
 f%TILTD=MAGL0%TILTD_EDGE(1)
 f%EDGE(1)=MAGL0%TILTD_EDGE(2)
 f%EDGE(2)=MAGL0%TILTD_EDGE(3)

 f%KILL_ENT_SPIN=MAGL0%KILL_SPIN(1)
 f%KILL_EXI_SPIN=MAGL0%KILL_SPIN(2)

 f%KILL_ENT_FRINGE=MAGL0%KIN_KEX_BENDFRINGE_EXACT(1)
 f%KILL_EXI_FRINGE=MAGL0%KIN_KEX_BENDFRINGE_EXACT(2)
 f%bend_fringe=MAGL0%KIN_KEX_BENDFRINGE_EXACT(3)
 f%EXACT=MAGL0%KIN_KEX_BENDFRINGE_EXACT(4)

 f%METHOD=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(1)
 f%NST=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(2)
 f%NMUL=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(3)
 f%permFRINGE=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(4)
 f%highest_fringe=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(5)

endif
endif
end subroutine MC_MC0

subroutine  el_el0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
character(nlp+3) nname
integer n,np,no,inf,i

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 ELE0%KIND=F%KIND
! ELE0%name_vorname(1)="'"//f%name//"' "
! ELE0%name_vorname(2)="'"//f%vorname//"' "

 ELE0%name_vorname(1)= f%name 
 ELE0%name_vorname(2)= f%vorname 
  call context(ELE0%name_vorname(1))
  call context(ELE0%name_vorname(2))
 if(.not.old_name_vorname) then
  call context(ELE0%name_vorname(1),dollar=my_true)
  call context(ELE0%name_vorname(2),dollar=my_true)
 endif
 ele0%an=0.0_dp
 ele0%bn=0.0_dp
if(f%p%nmul>0) then
  ele0%an(1:f%p%nmul)=f%an(1:f%p%nmul)
  ele0%bn(1:f%p%nmul)=f%bn(1:f%p%nmul)
endif
 ele0%VOLT_FREQ_PHAS=0.0_dp
 ele0%B_SOL=0.0_dp
 
   ele0%fint_hgap_h1_h2_va_vs(1)=f%fint(1)
   ele0%fint_hgap_h1_h2_va_vs(2)=f%fint(2)
   ele0%fint_hgap_h1_h2_va_vs(3)=f%hgap(1)
   ele0%fint_hgap_h1_h2_va_vs(4)=f%hgap(2)
   ele0%fint_hgap_h1_h2_va_vs(5)=f%h1
   ele0%fint_hgap_h1_h2_va_vs(6)=f%h2
   ele0%fint_hgap_h1_h2_va_vs(7)=f%va
   ele0%fint_hgap_h1_h2_va_vs(8)=f%vs


   ele0%L=f%L
   IF(ASSOCIATED(f%B_SOL)) ele0%B_SOL=f%B_SOL
 
 if(associated(f%volt)) ele0%VOLT_FREQ_PHAS(1)=f%VOLT
 if(associated(f%FREQ)) ele0%VOLT_FREQ_PHAS(2)=f%FREQ
 if(associated(f%PHAS)) ele0%VOLT_FREQ_PHAS(3)=f%PHAS
 if(associated(f%THIN)) ele0%THIN=f%THIN


ele0%slow_ac= f%slow_ac
ele0%recut_even_electric_MIS(1) = f%recut
ele0%recut_even_electric_MIS(2) = f%even
ele0%recut_even_electric_MIS(3) = f%electric
ele0%recut_even_electric_MIS(4) = f%MIS
 ele0%usebf_do1bf(1)=f%useb
 ele0%usebf_do1bf(2)=f%usef 
 ele0%skipptcbf(1)=f%skip_ptc_b 
 ele0%skipptcbf(2)=f%skip_ptc_f 
 ele0%usebf_do1bf(3)=f%do1mapb 
 ele0%usebf_do1bf(4)=f%do1mapf
 ele0%filef=trim(f%filef)
 ele0%fileb=TRIM(f%fileb)
 



if(associated(f%forward)) then

  if(present(mf)) then
   call kanalnummer(inf,f%filef)
    call print_tree_elements(f%forward,inf)
   close(inf)
 endif
endif
if(associated(f%backward)) then
 if(present(mf)) then
   call kanalnummer(inf,f%fileb)
    call print_tree_elements(f%backward,inf)
   close(inf)
 endif
endif


    if(present(mf)) then
     write(mf,NML=ELEname)
    endif  


else



    if(present(mf)) then
     read(mf,NML=ELEname)
    endif   


 F%KIND=ELE0%KIND  
  call context(ELE0%name_vorname(1))
  call context(ELE0%name_vorname(2))
nname=ELE0%name_vorname(1)
 f%name=nname(1:nlp)
nname=ELE0%name_vorname(2)
 f%vorname=nname(1:nlp)
if(f%p%nmul>0) then
 f%an=0.0_dp
 f%bn=0.0_dp
 f%an(1:f%p%nmul)=ele0%an(1:f%p%nmul)
 f%bn(1:f%p%nmul)=ele0%bn(1:f%p%nmul)
endif

f%fint(1)= ele0%fint_hgap_h1_h2_va_vs(1)
f%fint(2)=ele0%fint_hgap_h1_h2_va_vs(2)
f%hgap(1)= ele0%fint_hgap_h1_h2_va_vs(3)
f%hgap(2)= ele0%fint_hgap_h1_h2_va_vs(4)
f%h1  = ele0%fint_hgap_h1_h2_va_vs(5)
f%h2  = ele0%fint_hgap_h1_h2_va_vs(6)
f%va  = ele0%fint_hgap_h1_h2_va_vs(7)
f%vs  = ele0%fint_hgap_h1_h2_va_vs(8)


if(f%kind==kind4.or.f%kind==kind21) then ! cavities
 if(.not.associated(f%volt)) then
  ALLOCATE(f%VOLT,f%FREQ,f%PHAS,f%DELTA_E,f%THIN,f%lag)
 endif
  f%VOLT=ele0%VOLT_FREQ_PHAS(1)
  f%FREQ=ele0%VOLT_FREQ_PHAS(2)
  f%PHAS=ele0%VOLT_FREQ_PHAS(3)
  f%THIN=ele0%THIN
  f%delta_e=0.0_dp
endif

if(f%kind==KIND15) then   ! electrip separetor
 if(.not.associated(f%volt)) then
  ALLOCATE(f%VOLT,f%PHAS)
 endif
  f%VOLT=ele0%VOLT_FREQ_PHAS(1)
  f%PHAS=ele0%VOLT_FREQ_PHAS(3)
endif

    if(f%kind==kind22) then  ! helical dipole
       if(.not.associated(f%freq)) ALLOCATE(f%FREQ,f%PHAS)
       f%FREQ=ele0%VOLT_FREQ_PHAS(2)
       f%PHAS=ele0%VOLT_FREQ_PHAS(3)
    endif

 f%slow_ac = ele0%slow_ac
 f%recut = ele0%recut_even_electric_MIS(1)
 f%even = ele0%recut_even_electric_MIS(2)
 f%electric = ele0%recut_even_electric_MIS(3)
 f%MIS = ele0%recut_even_electric_MIS(4)
 solve_electric=f%electric
   F%L=ele0%L



!     logical(lp) usebf_do1bf(4)!
!	 integer skipptcbf(2)

 f%useb=ele0%usebf_do1bf(1)
 f%usef=ele0%usebf_do1bf(2)
 f%skip_ptc_b=ele0%skipptcbf(1)
 f%skip_ptc_f=ele0%skipptcbf(2)
 f%do1mapb=ele0%usebf_do1bf(3) 
 f%do1mapf=ele0%usebf_do1bf(4)
 f%fileb= ele0%fileb
 f%filef= ele0%filef
if(ele0%filef/=' '.and.readingmaps) then
 if(.not.associated(f%forward)) then 
  allocate(f%forward(3))
 endif
 call kanalnummer(inf,ele0%filef)

  do i=1,3
    read(inf,*) n,np,no
    CALL ALLOC_TREE(f%forward(i),N,NP)
    f%forward(i)%N=n
    f%forward(i)%NP=np
    f%forward(i)%no=no
    call read_tree_element(f%forward(i),inf)
  enddo
close(inf)
endif

if(ele0%fileb/=' '.and.readingmaps) then
 if(.not.associated(f%backward)) then 
  allocate(f%backward(3))
 endif
 call kanalnummer(inf,ele0%fileb)

  do i=1,3
    read(inf,*) n,np,no
    CALL ALLOC_TREE(f%backward(i),N,NP)
    f%backward(i)%N=n
    f%backward(i)%NP=np
    f%backward(i)%no=no
    call read_tree_element(f%backward(i),inf)
  enddo
close(inf)
 endif
   
    if(f%kind==kind3.or.f%kind==kind5) then   
        IF(.not.ASSOCIATED(f%B_SOL)) ALLOCATE(f%B_SOL);
       F%B_SOL=ele0%B_SOL
    endif

   
endif
endif
end subroutine el_el0



  subroutine print_ElementLIST(el,dir,mf)
    implicit none
    type(element), pointer :: el
    integer i
     logical(lp),optional ::  dir
     integer,optional :: mf
    character*255 line


    select case(el%kind)
    CASE(KIND0,KIND1,kind2,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17)
  case(kind3)
     call thin3_thin30(el,dir,mf)
    case(kind4)
        call cav4_cav40(EL,dir,mf)
    case(kind5)
        call sol5_sol50(EL,dir,mf)
    case(kind10)
        call tp10_tp100(EL,dir,mf)
    case(kindabell)
        call ab_ab0(EL,dir,mf)

    case(kind16,kind20)
        call k16_k160(EL,dir,mf)

    case(kind18)
!       WRITE(MF,*) " RCOLLIMATOR HAS AN INTRINSIC APERTURE "
  !     CALL  ap_aplist(el,dir,mf) 
    case(kind19)
!       WRITE(MF,*) " ECOLLIMATOR HAS AN INTRINSIC APERTURE "
!       CALL print_aperture(EL%ECOL19%A,mf)
    case(kindhel)
        call hel_hel0(EL,dir,mf)
    case(kind21)
        call tcav4_tcav40(EL,dir,mf)
!       WRITE(MF,*) el%cav21%PSI,el%cav21%DPHAS,el%cav21%DVDS
    case(KINDWIGGLER)
        call wig_wig0(el,dir,mf)
 !      call print_wig(el%wi,mf)
    case(KINDpa)
 !      call print_pancake(el%pa,mf)
    case default
       write(MF,*) " not supported in print_specific_element",el%kind
 !      stop 101
    end select
 
    CALL  ap_aplist(el,dir,mf) 

 

  END SUBROUTINE print_ElementLIST

   
subroutine  cav4_cav40(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

 cav0%f=0.0_dp
 cav0%PH=0.0_dp   
 cav0%N_BESSEL=F%c4%N_BESSEL
 cav0%NF=F%c4%NF
 cav0%CAVITY_TOTALPATH=F%c4%CAVITY_TOTALPATH
 cav0%phase0=F%c4%phase0
 cav0%t=F%c4%t
 cav0%always_on=F%c4%always_on
 cav0%f(1:F%c4%NF)=F%c4%f
 cav0%PH(1:F%c4%NF)=F%c4%PH
 cav0%A=F%c4%A
 cav0%R=F%c4%R
    if(present(mf)) then
     write(mf,NML=CAVname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=CAVname)
    endif   
N_CAV4_F=cav0%NF
CALL SETFAMILY(f)
 F%c4%N_BESSEL=cav0%N_BESSEL
 F%c4%NF =cav0%NF


 F%c4%CAVITY_TOTALPATH=cav0%CAVITY_TOTALPATH
 F%c4%phase0=cav0%phase0
 F%c4%t=cav0%t
 F%c4%always_on=cav0%always_on
 F%c4%f=cav0%f(1:F%c4%NF)

 
 F%c4%PH=cav0%PH(1:F%c4%NF)
 F%c4%A=cav0%A
 F%c4%R=cav0%R  
endif
endif
end subroutine cav4_cav40

subroutine  hel_hel0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG


hel0%fake_shift=F%he22%fake_shift
hel0%N_BESSEL=F%he22%N_BESSEL

    if(present(mf)) then
     write(mf,NML=helname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=helname)
    endif   
 
CALL SETFAMILY(f)
 F%he22%N_BESSEL=hel0%N_BESSEL
 F%he22%fake_shift=hel0%fake_shift

endif
endif
end subroutine hel_hel0

subroutine  wig_wig0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
integer n,ne
n=0
ne=0
if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG


 wig0%internal=F%wi%internal
 wig0%offset=F%wi%w%offset
 wig0%ex=f%wi%w%ex
 wig0%ey=f%wi%w%ey
wig0%n=0
 if(associated(f%wi%w%a)) wig0%n=size(f%wi%w%a)
 n=wig0%n
wig0%a=0.0_dp
wig0%f=0.0_dp
wig0%form=0.0_dp
wig0%k=0.0_dp
if(n>0) then
 wig0%a(1:n)=f%wi%w%a(1:n)
 wig0%f(1:n)=f%wi%w%f(1:n)
 wig0%form(1:n)=f%wi%w%form(1:n)
 wig0%k(1:3,1:n)=f%wi%w%k(1:3,1:n)
else
 wig0%a(1:n)=0
 wig0%f(1:n)=0
 wig0%form(1:n)=0
 wig0%k(1:3,1:n)=0
endif

wig0%ne=0
 if(associated(f%wi%w%ae)) wig0%ne=size(f%wi%w%ae)
 ne=wig0%ne
if(ne>0) then
 wig0%ae(1:ne)=f%wi%w%ae(1:ne)
 wig0%fe(1:ne)=f%wi%w%fe(1:ne)
 wig0%forme(1:ne)=f%wi%w%forme(1:ne)
 wig0%ke(1:3,1:ne)=f%wi%w%ke(1:3,1:ne)
else
 wig0%ae(1:ne)=0
 wig0%fe(1:ne)=0
 wig0%forme(1:ne)=0
 wig0%ke(1:3,1:ne)=0
endif

    if(present(mf)) then
     write(mf,NML=wigname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=wigname)
    endif   
    if(.not.associated(f%wi%internal)) allocate(f%wi%internal(6))
  F%wi%internal=wig0%internal 
  N=wig0%N
  ne=wig0%Ne

    call pointers_w(f%wi%w,N,ne)

 F%wi%w%offset=wig0%offset
 F%wi%w%ex=wig0%ex
 F%wi%w%ey=wig0%ey

if(n>0) then
 F%wi%w%a(1:N)=wig0%a(1:N)
 F%wi%w%f(1:N)=wig0%f(1:N)
 F%wi%w%form(1:N)=wig0%form(1:N)
 F%wi%w%k(1:3,1:N)=wig0%k(1:3,1:N)
endif

if(ne>0) then
 F%wi%w%ae(1:Ne)=wig0%ae(1:Ne)
 F%wi%w%fe(1:Ne)=wig0%fe(1:Ne)
 F%wi%w%forme(1:Ne)=wig0%forme(1:Ne)
 F%wi%w%ke(1:3,1:Ne)=wig0%ke(1:3,1:Ne)
endif

endif
endif
end subroutine wig_wig0

subroutine  tcav4_tcav40(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG


 tcav0%PSI_DPHAS_DVDS(1)=F%cav21%psi
 tcav0%PSI_DPHAS_DVDS(2)=F%cav21%dphas
 tcav0%PSI_DPHAS_DVDS(3)=F%cav21%dvds
 tcav0%always_on=F%cav21%always_on
    if(present(mf)) then
     write(mf,NML=tCAVname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=tCAVname)
    endif   

 F%cav21%psi=tcav0%PSI_DPHAS_DVDS(1)
 F%cav21%dphas=tcav0%PSI_DPHAS_DVDS(2)
 F%cav21%dvds=tcav0%PSI_DPHAS_DVDS(3)
 F%cav21%always_on=tcav0%always_on
endif

endif
end subroutine tcav4_tcav40

subroutine  sol5_sol50(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

    


 sol50%dx_dy_pitchx_pitchy(1)=F%s5%dx
 sol50%dx_dy_pitchx_pitchy(2)=F%s5%dy
 sol50%dx_dy_pitchx_pitchy(3)=F%s5%pitch_x
 sol50%dx_dy_pitchx_pitchy(4)=F%s5%pitch_y

    if(present(mf)) then
     write(mf,NML=sol50name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=sol50name)
    endif   
       IF(.NOT.ASSOCIATED(f%B_SOL)) then
        ALLOCATE(f%B_SOL)
          f%B_SOL=0.0_dp
       endif
       CALL SETFAMILY(f) 

 F%s5%dx= sol50%dx_dy_pitchx_pitchy(1)
 F%s5%dy= sol50%dx_dy_pitchx_pitchy(2)
 F%s5%pitch_x= sol50%dx_dy_pitchx_pitchy(3)
 F%s5%pitch_y= sol50%dx_dy_pitchx_pitchy(4)
endif
endif
end subroutine sol5_sol50


subroutine  thin3_thin30(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

    

 thin30%thin_h_foc=F%k3%thin_h_foc
 thin30%thin_v_foc=F%k3%thin_v_foc
 thin30%thin_h_angle=F%k3%thin_h_angle
 thin30%thin_v_angle=F%k3%thin_v_angle
 thin30%hf=F%k3%hf
 thin30%vf=F%k3%vf
 thin30%patch=F%k3%patch
 thin30%ls=F%k3%ls 
  thin30%dx_dy_pitchx_pitchy(1)=F%k3%dx
  thin30%dx_dy_pitchx_pitchy(2)=F%k3%dy
  thin30%dx_dy_pitchx_pitchy(3)=F%k3%pitch_x
  thin30%dx_dy_pitchx_pitchy(4)=F%k3%pitch_y
    if(present(mf)) then
     write(mf,NML=thin30name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=thin30name)
    endif   
       IF(.NOT.ASSOCIATED(f%B_SOL)) then
        ALLOCATE(f%B_SOL)
          f%B_SOL=0.0_dp
       endif
       CALL SETFAMILY(f) 
 f%k3%thin_h_foc=thin30%thin_h_foc
 f%k3%thin_v_foc=thin30%thin_v_foc
 f%k3%thin_h_angle=thin30%thin_h_angle
 f%k3%thin_v_angle=thin30%thin_v_angle
 f%k3%hf=thin30%hf
 f%k3%vf=thin30%vf
 f%k3%patch=thin30%patch 
 f%k3%ls=thin30%ls
 F%k3%dx= thin30%dx_dy_pitchx_pitchy(1)
 F%k3%dy= thin30%dx_dy_pitchx_pitchy(2)
 F%k3%pitch_x= thin30%dx_dy_pitchx_pitchy(3)
 F%k3%pitch_y= thin30%dx_dy_pitchx_pitchy(4)
endif
endif
end subroutine thin3_thin30

subroutine  tp10_tp100(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
 

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

 tp100%ae=0.0_dp
 tp100%be=0.0_dp
 
 tp100%DRIFTKICK=F%tp10%DRIFTKICK
 if(f%electric) then
  tp100%ae(1:size(F%tp10%ae))=F%tp10%ae
  tp100%be(1:size(F%tp10%be))=F%tp10%be
 endif

     if(present(mf)) then
     write(mf,NML=tp100name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=tp100name)
    endif   
   CALL SETFAMILY(f)
 if(f%electric) then
  F%tp10%ae=tp100%ae(1:sector_nmul_max) 
  F%tp10%be=tp100%be(1:sector_nmul_max)
  call GETAEBE(f%TP10)
 endif

 F%tp10%DRIFTKICK=tp100%DRIFTKICK
endif
endif
end subroutine tp10_tp100


subroutine  ab_ab0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
 

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

 ab0%N_m=0 
 ab0%dz_t_te=0 
 ab0%b=0 
 ab0%e=0 
 ab0%scale_angc_xc_vc_dc_hc=0 
 
 ab0%n_m(1)= F%ab%n
 ab0%n_m(2)= F%ab%m
 ab0%b(1:ab0%n_m(2) ,1:ab0%n_m(2))= F%ab%b
 ab0%dz_t_te(2*ab0%n_m(2)+3:3*ab0%n_m(2)+3)=F%ab%te
ab0%dz_t_te(ab0%n_m(2)+2:2*ab0%n_m(2)+2)= F%ab%t
ab0%dz_t_te(1:ab0%n_m(2)+1)= F%ab%dz
ab0%scale_angc_xc_vc_dc_hc(1)=F%ab%SCALE 
 ab0%scale_angc_xc_vc_dc_hc(2)=F%ab%angc
 ab0%scale_angc_xc_vc_dc_hc(3)=F%ab%xc
 ab0%scale_angc_xc_vc_dc_hc(4)=F%ab%vc
 ab0%scale_angc_xc_vc_dc_hc(5)=F%ab%dc
 ab0%scale_angc_xc_vc_dc_hc(6)=F%ab%hc
 

! if(f%electric) then
!  tp100%ae(1:size(F%tp10%ae))=F%tp10%ae
!  tp100%be(1:size(F%tp10%be))=F%tp10%be
! endif

     if(present(mf)) then
     write(mf,NML=ab0name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=ab0name)
    endif  
   n_abell=ab0%n_m(1) 
   m_abell=ab0%n_m(2)
   CALL SETFAMILY(f)
 !if(f%electric) then
 ! F%tp10%ae=tp100%ae(1:sector_nmul_max) 
 ! F%tp10%be=tp100%be(1:sector_nmul_max)
 ! call GETAEBE(f%TP10)
 !endif
 F%ab%n=ab0%n_m(1) 
 F%ab%m=ab0%n_m(2) 
 F%ab%b=ab0%b(1:ab0%n_m(2) ,1:ab0%n_m(2))
 F%ab%te=ab0%dz_t_te(2*ab0%n_m(2)+3:3*ab0%n_m(2)+3)
 F%ab%t=ab0%dz_t_te(ab0%n_m(2)+2:2*ab0%n_m(2)+2)
 F%ab%dz=ab0%dz_t_te(1:ab0%n_m(2)+1)
F%ab%SCALE= ab0%scale_angc_xc_vc_dc_hc(1)
F%ab%angc= ab0%scale_angc_xc_vc_dc_hc(2)
F%ab%xc= ab0%scale_angc_xc_vc_dc_hc(3)
F%ab%vc= ab0%scale_angc_xc_vc_dc_hc(4)
F%ab%dc= ab0%scale_angc_xc_vc_dc_hc(5)
F%ab%hc= ab0%scale_angc_xc_vc_dc_hc(6)
 
endif
endif
end subroutine ab_ab0

subroutine  k16_k160(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 
 k160%DRIFTKICK=F%k16%DRIFTKICK
 k160%LIKEMAD=F%k16%LIKEMAD
     if(present(mf)) then
     write(mf,NML=k160name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=k160name)
    endif   
 F%k16%DRIFTKICK=k160%DRIFTKICK
 F%k16%LIKEMAD=k160%LIKEMAD
endif
endif
end subroutine k16_k160

subroutine  ap_aplist(f,dir,mf)
implicit none
type(element), target :: f
type(MADX_APERTURE), pointer :: a
logical(lp),optional ::  dir
integer,optional :: mf
logical(LP) :: here
CHARACTER(120) LINE

here=associated(f%p%APERTURE)

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 
if(here) then
 a=>f%p%aperture

    aplist%kind=a%kind
    aplist%r=a%r
    aplist%x=a%x
    aplist%y=a%y
    aplist%dx=a%dx
    aplist%dy=a%dy
     if(present(mf)) then
     Write(mf,*) " APERTURE "  
     write(mf,NML=aperturename)
    endif   
else
 Write(mf,*) " NO APERTURE "
endif
 
 else   ! dir=false
 !  if(present(mf)) then     
 !     READ(MF,'(a120)') LINE 
 !     CALL CONTEXT(LINE)
 !  ENDIF
    IF(aplist%on) THEN
       IF(.NOT.HERE) THEN
          CALL alloc(f%p%aperture)
           a=>f%p%aperture
       ENDIF
    
    
    if(present(mf)) then     
     read(mf,NML=aperturename)
    endif   
      a%kind=aplist%kind
      a%r=aplist%r
      a%x=aplist%x
      a%y=aplist%y
      a%dx=aplist%dx
      a%dy=aplist%dy      
    ENDIF
endif
endif
end subroutine ap_aplist


subroutine  r_ap_aplist(mf)
implicit none
integer  mf
CHARACTER(120) LINE
   
      READ(MF,'(a120)') LINE 
      CALL CONTEXT(LINE)
    IF(LINE(1:2)/='NO') THEN
        read(mf,NML=aperturename)
       aplist%on=.true.
   !    write(6,nml=aperturename)
    else
     aplist%on=.false.
    endif
    
end subroutine r_ap_aplist

!!!!!!!!!!!!!!!!!!!!!    pointed lattices !!!!!!!!!!!!!!!!!!!!!

subroutine  print_universe_pointed(ud,un,filename,com)
! call print_universe_pointed(m_u,m_t,'junk3.txt')
implicit none
type(mad_universe),target :: ud,un
type(layout),pointer :: r
character(6), optional ::com
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,j,i0,j0,i1,j1,jb,MF
character (6) comt
logical(lp) before,just
logical, allocatable :: dna(:)
character(nlp) name

allocate(dna(ud%n))

comt='REWIND'
if(present(com)) comt=com

call kanalnummer(mf)
open(unit=mf,file=filename,position=comt) !,recl=4000)

call TIE_MAD_UNIVERSE(ud)


r=>un%start

write(mf,*) un%n, ud%n," trackable and DNA Layouts"


do i=1,un%n
dna=.false.
write(mf,'(a120)') r%name
p0=>r%start
p=>p0
call locate_in_universe(p,i0,j0)
jb=j0
j1=j0
!   write(mf,*) i,r%n," New "
   write(mf,*)  r%n," Number elements in pointed layout "

dna(i0)=.true.

  call Print_initial_chart(p,mf)
 name=p%mag%name
 if(.not.old_name_vorname) then
  call context(name,dollar=my_true)
 endif
  write(mf,'(1x,i4,1x,i8,1x,i2,1x,a24,1x,i8)') i0,p%dir*j0,p%patch%patch,name, 1
  call fib_fib0(p,my_true,mf)
  before=my_false
  just=my_false
 if(p%patch%patch/=0) call patch_patch0(p%patch,my_true,mf)

 do j=2,r%n
  p=>p%next
    jb=j1
   call locate_in_universe(p,i1,j1)
name=p%mag%name
 if(.not.old_name_vorname) then
  call context(name,dollar=my_true)
 endif
 write(mf,'(1x,i4,1x,i8,1x,i2,1x,a24,1x,i8)') i1,p%dir*j1,p%patch%patch,name, j
    dna(i1)=my_true      
       just=my_false
       if(before) then 
        call fib_fib0(p,my_true,mf)
        just=my_true
        before=my_false
       endif
    if(p%patch%patch/=0) then
       if(.not.just) call fib_fib0(p,my_true,mf)  
       call patch_patch0(p%patch,my_true,mf)   
       before=my_true   
    endif
 enddo
 write(mf,*) " DNA SUMMARY "
do i0=1,ud%n
 write(mf,*) i0,dna(i0)
enddo
write(mf,*) " !!!!!!! End of Pointed Layout !!!!!!!"

 r=>r%next

enddo

close(mf)

deallocate(dna)

end subroutine  print_universe_pointed


subroutine  read_universe_pointed(ud,un,filename)
implicit none
type(mad_universe),target :: ud,un
type(layout),pointer :: r,rd
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,j,i0,MF,n,n_u,k(3),mypause,ipause,n_d,size_dna
integer pos
character(120) line,name1
logical(lp) doneit,first
character(nlp) name
logical, allocatable :: dna(:)
 
call kanalnummer(mf,filename,old=.true.)
!open(unit=mf,file=filename)


call TIE_MAD_UNIVERSE(ud)


read(mf,*)n_u,n_d

allocate(dna(n_d))

do i=1,n_u
   read(mf,'(a120)') name1 
!   read(mf,*) i0,n 
   read(mf,*)  n 

  read(mf,'(a120)') line
 call read_initial_chart(mf)
  read(mf,'(a120)') line

call append_empty_layout(un) 
!call set_up(un%end)  !
 
       R => un%end
r%name=name1

 do j=1,n 
  !     read(mf,'(1x,i4,1x,i8,1x,i2,1x,a24)' ) k    ,name
       read(mf,*) k    ,name
!        write(6,*)  k  
!        write(6,*)   name
!pause 234
if(j==1.or.k(3)>0) then
 read(mf,NML=FIBRENAME)
else
 if(r%end%patch%patch>0) read(mf,NML=FIBRENAME)    ! previous had a patch
endif

       call MOVE_TO_LAYOUT_I( ud,rd,k(1) )
       pos=iabs(k(2))   
       call move_to_p_safe( rd,p,POS)
       call append_point(r, p)  
       if( p%mag%name/=name) then
         write(6,*) j," serious error in read_universe_pointed "
         write(6,*)  k ,name
         write(6,*) i,p%mag%name ,pos
           read(mf,'(a120)') line
           write(6,'(a120)') line
           read(mf,'(a120)') line
           write(6,'(a120)') line
           read(mf,'(a120)') line
           write(6,'(a120)') line
          ipause=mypause(666)
         stop 666   
       endif

        if(k(3)/=0) then
        read(mf,NML=patchname)
        call patch_patch0(r%end%patch,my_false)
        endif
        r%end%patch%patch=k(3)

       call fib_fib0(r%end,my_false) 
       if(k(2)>0) then   ! because fib0 could have the wrong dir
        r%end%dir=1
       else
        r%end%dir=-1
       endif
 enddo

    write(6,*) r%index,ud%n,un%n
    r%closed=my_true
    doneit=my_true
    call ring_l(r,doneit)

  p=>r%start
  p%chart%f%a=a_
  p%chart%f%ent=ent_
! else
  p%chart%f%b=b_
  p%chart%f%exi=exi_
call survey(r)

read(mf,'(a120)') line

size_dna=0
do j=1,n_d
 read(mf,*) i0,dna(j)
 if(dna(j)) size_dna=size_dna+1
enddo

 read(mf,'(a10)') line(1:10)

! Setting up the dna
 allocate(r%dna(size_dna))

  i0=0
 do j=1,n_d
  if(dna(j))  then
      i0=i0+1
      call MOVE_TO_LAYOUT_I( ud,rd,j )    
     r%dna(i0)%L=>rd
  endif
 enddo 
 if(i0/=size_dna) then
  write(6,*) i0,size_dna, "error in read_universe_pointed " 
  stop 567
 endif

enddo



close(mf)

deallocate(dna)

end subroutine  read_universe_pointed


 subroutine Print_initial_chart(f,mf)
 implicit none
 type(fibre), target :: f
 integer mf

 write(mf,*) " $$$$$$$$$$$$$$$$$ INITIAL CHART $$$$$$$$$$$$$$$$$"
!  IF(F%DIR==1) THEN
   write(mf,*) f%chart%f%A
   write(mf,*) f%chart%f%ENT(1,1:3)
   write(mf,*) f%chart%f%ENT(2,1:3)
   write(mf,*) f%chart%f%ENT(3,1:3)
!  ELSE
   write(mf,*) f%chart%f%B
   write(mf,*) f%chart%f%EXI(1,1:3)
   write(mf,*) f%chart%f%EXI(2,1:3)
   write(mf,*) f%chart%f%EXI(3,1:3)
!  ENDIF
 write(mf,*) " $$$$$$$$$$$$$$$$$ END OF INITIAL CHART $$$$$$$$$$$$$$$$$"

 end subroutine Print_initial_chart

 subroutine read_initial_chart(mf)
 implicit none
 integer mf
 
   read(mf,*) A_
   read(mf,*) ENT_(1,1:3)
   read(mf,*) ENT_(2,1:3)
   read(mf,*) ENT_(3,1:3)

   read(mf,*) B_
   read(mf,*) EXI_(1,1:3)
   read(mf,*) EXI_(2,1:3)
   read(mf,*) EXI_(3,1:3)

 end subroutine read_initial_chart


subroutine create_dna(ud,ut)
implicit none
type(mad_universe), target :: ut,ud
type(layout), pointer :: r,rd
type(fibre), pointer :: ps
integer i,j,j1,j2,k,kn
logical(lp), allocatable :: here(:)



call TIE_MAD_UNIVERSE(ut)

allocate(here(ud%n))

r=>ut%start
do i=1,ut%n

here=my_false

ps => r%start

do j=1,r%n
 call locate_in_universe(ps,j1,j2)
 here(j1)=my_true
ps=>ps%next
enddo

kn=0
do k=1,size(here)
 if(here(k)) kn=kn+1
enddo

if(associated(r%DNA)) then 
deallocate(r%DNA)
write(6,*) "deallocated DNA"
endif
allocate(r%DNA(kn))

rd=>ud%start
kn=0
do k=1,size(here)
 if(here(k)) then
    kn=kn+1
    r%dna(kn)%L=>rd
    r%dna(kn)%counter=k
  endif
 rd=>rd%next
enddo


 r=>r%next

enddo

deallocate(here)

end subroutine create_dna

subroutine zero_ele0
implicit none

    ele0%name_vorname=' ' 
	ele0%L=0;ele0%B_SOL=0;
	ele0%an=0;ele0%bn=0;
    ele0%VOLT_FREQ_PHAS=0
    ele0%THIN=.false. 
    ele0%fint_hgap_h1_h2_va_vs=0
	ele0%recut_even_electric_MIS=.false.
    ele0%slow_ac=0
    ele0%usebf_do1bf=.false.
    ele0%skipptcbf=0
    ele0%filef=' '
    ele0%fileb=' '


!     logical(lp) usebf_do1bf(4)!
!	 integer skipptcbf(2)

end subroutine zero_ele0

subroutine zero_fib0
implicit none

    fib0%GAMMA0I_GAMBET_MASS_AG=0 
    fib0%CHARGE=0
    fib0%DIR=0
    fib0%patch=0

end subroutine zero_fib0

subroutine zero_CHART0
implicit none
    CHART0%D_IN=0; CHART0%D_OUT=0;CHART0%ANG_IN=0; CHART0%ANG_OUT=0;

end subroutine zero_CHART0

subroutine zero_MAGL0
implicit none
 
    MAGL0%LC_LD_B0_P0=0  ! LC LD B0 P0C
    MAGL0%TILTD_EDGE=0   ! TILTD EDGE
    MAGL0%KIN_KEX_BENDFRINGE_EXACT=.false. ! KILL_ENT_FRINGE, KILL_EXI_FRINGE, bend_fringe,EXACT
	MAGL0%METHOD_NST_NMUL_permfringe_highest=0  ! METHOD,NST,NMUL,permfringr, highest_fringe
    MAGL0%kill_spin=.false. 

end subroutine zero_MAGL0

subroutine zero_patch0
implicit none
 
 
     patch0%A_X1=0;patch0%A_X2=0;patch0%B_X1=0;patch0%B_X2=0
     patch0%A_D=0
     patch0%B_D=0
     patch0%A_ANG=0
     patch0%B_ANG=0 
     patch0%A_T=0
     patch0%B_T=0
     patch0%A_L=0
     patch0%B_L=0
     patch0%ENERGY=0
     patch0%TIME=0
     patch0%GEOMETRY=0

end subroutine zero_patch0

subroutine specify_element_type(f,mf)
implicit none
type(fibre), target ::f
integer mf

if(associated(f%mag%D0)) then

Write(mf,*) f%mag%kind, "TYPE(DRIFT1), POINTER :: D0"
return
endif
if(associated(f%mag%K2)) then

Write(mf,*) f%mag%kind,"TYPE(DKD2), POINTER :: K2 "
return
endif
if(associated(f%mag%K3)) then

Write(mf,*) f%mag%kind,"TYPE(KICKT3), POINTER :: K3"
return
endif

if(associated(f%mag%C4)) then

Write(mf,*) f%mag%kind,"TYPE(CAV4), POINTER :: C4 "
return
endif

if(associated(f%mag%S5)) then

Write(mf,*) f%mag%kind,"TYPE(SOL5), POINTER :: S5 "
return
endif

if(associated(f%mag%T6)) then

Write(mf,*) f%mag%kind,"TYPE(KTK), POINTER :: T6 "
return

endif
if(associated(f%mag%T6)) then

Write(mf,*) f%mag%kind,"TYPE(KTKP), POINTER :: T7 "
return
endif

if(associated(f%mag%T7)) then

Write(mf,*) f%mag%kind,"TYPE(TKTF), POINTER :: T7 "
return
endif

if(associated(f%mag%S8)) then

Write(mf,*) f%mag%kind,"TYPE(NSMI), POINTER :: S8"
return
endif

if(associated(f%mag%S9)) then

Write(mf,*) f%mag%kind,"TYPE(SSMI), POINTER :: S9"
return
endif

if(associated(f%mag%TP10)) then

Write(mf,*) f%mag%kind,"TYPE(TEAPOT), POINTER :: TP10 "
return
endif

if(associated(f%mag%MON14)) then

Write(mf,*) f%mag%kind,"TYPE(MON), POINTER :: MON14"
return
endif

if(associated(f%mag%SEP15)) then

Write(mf,*) f%mag%kind,"TYPE(ESEPTUM), POINTER :: SEP15"
return
endif

if(associated(f%mag%K16)) then

Write(mf,*) f%mag%kind,"TYPE(STREX), POINTER :: K16 "
return
endif

if(associated(f%mag%ENGE17)) then

Write(mf,*) f%mag%kind,"TYPE(ENGE), POINTER :: ENGE17"
return
endif

if(associated(f%mag%RCOL18)) then

Write(mf,*) f%mag%kind,"TYPE(RCOL), POINTER :: RCOL18"
return
endif

if(associated(f%mag%ECOL19)) then

Write(mf,*) f%mag%kind,"TYPE(ECOL), POINTER :: ECOL19"
return
endif

if(associated(f%mag%CAV21)) then

Write(mf,*) f%mag%kind,"TYPE(CAV_TRAV), POINTER :: CAV21"
return
endif

if(associated(f%mag%WI)) then

Write(mf,*) f%mag%kind,"TYPE(CAV_TRAV), POINTER :: WI"
return
endif
 
if(associated(f%mag%PA)) then

Write(mf,*) f%mag%kind,"TYPE(PANCAKE), POINTER :: PA"
return
endif

if(associated(f%mag%HE22)) then

Write(mf,*) f%mag%kind,"TYPE(HELICAL_DIPOLE), POINTER :: HE22"
return
endif

if(associated(f%mag%SDR)) then

Write(mf,*) f%mag%kind,"TYPE(HELICAL_DIPOLE), POINTER :: SDR"
return
endif
 
end subroutine specify_element_type 

end module madx_keywords

