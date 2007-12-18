!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module madx_keywords
  use S_fitting
  implicit none
  public
  logical(lp)::mad8=my_false
  integer :: ifield_name=0
  logical(lp),private :: do_survey =.false.

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




contains




  subroutine create_fibre(el,key,EXCEPTION,magnet_only)
    implicit none
    integer ipause, mypause
    type(fibre), target, intent(inout)::el
    logical(lp), optional :: magnet_only
    type(keywords) key
    type(el_list) blank
    character*255 magnet
    character*17 MODEL
    INTEGER EXCEPTION  !,NSTD0,METD0
    LOGICAL(LP) EXACT0,magnet0
    logical(lp) FIBRE_flip0,MAD0
    logical(lp) :: t=.true.,f=.false.
    INTEGER FIBRE_DIR0,IL,multipi
    real(dp) e1_true


    IL=15

    if(present(magnet_only)) then
       magnet0=magnet_only
    else
       magnet0=.false.
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

    SELECT CASE(magnet(1:IL))
    CASE("DRIFT          ")
       BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("SOLENOID       ")
       BLANK=SOLENOID(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("QUADRUPOLE     ")
       BLANK=QUADRUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("SEXTUPOLE     ")
       BLANK=SEXTUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("OCTUPOLE      ")
       BLANK=OCTUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("SBEND         ")
       BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("TRUERBEND     ")

       e1_true= KEY%LIST%b0/two+ KEY%LIST%t1
       BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,e1=e1_true,list=KEY%LIST)

    CASE("WEDGRBEND     ")

       BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,e1=KEY%LIST%t1,e2=KEY%LIST%t2,list=KEY%LIST)

    CASE("RBEND         ")
       KEY%LIST%T1=KEY%LIST%T1+KEY%LIST%B0/two
       KEY%LIST%T2=KEY%LIST%T2+KEY%LIST%B0/two
       BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("KICKER         ","VKICKER        ","HKICKER        ")
       BLANK=KICKER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("MONITOR        ")
       BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("HMONITOR        ")
       BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND12;
    CASE("VMONITOR       ")
       BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND13;
    CASE("INSTRUMENT     ")
       BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND14;
    CASE("MARKER         ")
       BLANK=MARKER(KEY%LIST%NAME)
    CASE("CHANGEREF      ")
       BLANK=CHANGEREF(KEY%LIST%NAME,KEY%LIST%ANG,KEY%LIST%T,KEY%LIST%PATCHG)
    CASE("RFCAVITY       ")
       BLANK=RFCAVITY(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("TWCAVITY       ")
       BLANK=TWCAVITY(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("ELSEPARATOR    ")
       BLANK=ELSEPARATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("MULTIPOLE_BLOCK","MULTIPOLE      ")
       BLANK=MULTIPOLE_BLOCK(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("SMI            ","SINGLE_LENS    ")
       BLANK=SINGLE_LENS(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("RCOLLIMATOR    ")
       BLANK=RCOLLIMATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("ECOLLIMATOR    ")
       BLANK=ECOLLIMATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("WIGGLER        ")
       BLANK=WIGGLER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       !    CASE("TAYLORMAP      ")
       !       IF(KEY%LIST%file/=' '.and.KEY%LIST%file_rev/=' ') THEN
       !          BLANK=TAYLOR_MAP(KEY%LIST%NAME,FILE=KEY%LIST%file,FILE_REV=KEY%LIST%file_REV,t=tilt.is.KEY%tiltd)
       !       ELSEIF(KEY%LIST%file/=' '.and.KEY%LIST%file_rev==' ') THEN
       !          BLANK=TAYLOR_MAP(KEY%LIST%NAME,FILE=KEY%LIST%file,t=tilt.is.KEY%tiltd)
       !       ELSEIF(KEY%LIST%file==' '.and.KEY%LIST%file_rev/=' ') THEN
       !          BLANK=TAYLOR_MAP(KEY%LIST%NAME,FILE_REV=KEY%LIST%file_REV,t=tilt.is.KEY%tiltd)
       !       ELSE
       !          BLANK=TAYLOR_MAP(KEY%LIST%NAME,t=tilt.is.KEY%tiltd)
       !       ENDIF
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
  end subroutine create_fibre

  subroutine zero_key(key)
    implicit none

    type(keywords) , intent(out):: key
    real h
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
    key%tiltd=ZERO
    key%LIST=0

  end subroutine zero_key

  !  PRINTING FIBRES FOR FLAT FILES
  subroutine print_COMPLEX_SINGLE_STRUCTURE(L,FILENAME,LMAX0,NL)
    implicit none
    character(*) filename
    integer I,MF,N,n_l,ncon1,ncon2,j
    type(LAYOUT), TARGET :: L
    type(LAYOUT), pointer :: CL
    type(FIBRE), pointer :: P
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
    integer I,MF
    INTEGER, OPTIONAL :: MFF
    type(LAYOUT), TARGET :: L
    type(FIBRE), pointer :: P
    REAL(DP),OPTIONAL :: LMAX0
    character*255 line
    logical(lp) print_temp

    IF(PRESENT(MFF)) THEN
       MF=MFF
    ELSE
       call kanalnummer(mf)
       open(unit=mf,file=filename)
    ENDIF

    IF(PRESENT(LMAX0)) THEN
       WRITE(MF,*) L%N, LMAX0, " NUMBER OF FIBRES AND L_MAX  "
    ELSE
       WRITE(MF,*) L%N, 0, " NUMBER OF FIBRES AND L_MAX  "
    ENDIF
    if(l%name(1:1)/=' ') then
       write(MF,'(a17,a16)') " GLOBAL DATA FOR ",l%name
    else
       write(MF,*) " $$$$$$$$$ GLOBAL DATA  $$$$$$$$$"
    endif

    write(MF,*) l%start%mass,L%START%mag%p%p0c, "MASS, P0C"
    write(MF,*) phase0,stoch_in_rec,l%start%charge, " PHASE0, STOCH_IN_REC, CHARGE"
    write(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING, &
         "CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING"
    write(line,*) SECTOR_NMUL_MAX,SECTOR_NMUL,&
         OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE,&
         " SECTOR_NMUL_MAX,SECTOR_NMUL,OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE"
    write(mf,'(a255)')line
    write(MF,*) wedge_coeff, " wedge_coeff"
    write(MF,*) MAD8_WEDGE, " MAD8_WEDGE"
    write(MF,*) " $$$$$$$ END GLOBAL DATA $$$$$$$$$$"


    P=>L%START
    DO I=1,L%N
       if(i==1) then
          print_temp=print_frame
          print_frame=my_true
       endif
       CALL print_FIBRE(P,mf)
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
    integer mf,I,N
    integer, optional :: mf1
    type(LAYOUT), TARGET :: L
    type(FIBRE), pointer :: P
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0
    LOGICAL(LP) RING_IT,doneit
    character*120 line
    real(dp) p0c,MASSF
    type(internal_state) original

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
    read(MF,'(a120)') line
    call context(line)

    if(index(line,"FOR")/=0) then
       l%name=line(index(line,"FOR")+3:index(line,"FOR")+2+nlp)
    endif
    read(MF,*) MASSF,p0c
    read(MF,*) phase0,stoch_in_rec,initial_charge
    read(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING
    read(MF,*) SECTOR_NMUL_MAX,SECTOR_NMUL,OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE
    read(MF,*) wedge_coeff
    read(MF,*) MAD8_WEDGE
    read(MF,'(a120)') line
    original=default
    if(allocated(s_b)) then
       firsttime_coef=.true.
       deallocate(s_b)
    endif
    !    L%MASS=MASSF
    MASSF=MASSF/pmae
    CALL MAKE_STATES(MASSF)
    default=original
    call Set_madx(p0c=p0c)
    DO I=1,N
       CALL APPEND_CLONE(L)
       CALL READ_FIBRE(L%END,mf)
       CALL COPY(L%END%MAG,L%END%MAGP)
    ENDDO

    if(.not.present(mf1)) CLOSE(MF)

    L%closed=RING_IT

    doneit=.true.
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
    do_survey=.false.
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
    integer mf,I,siam_pos,siam_index
    type(FIBRE), pointer :: m
    siam_pos=0
    siam_index=0
    if(associated(m%mag%siamese)) then
       siam_index=m%mag%siamese%parent_fibre%parent_layout%index
       siam_pos=m%mag%siamese%parent_fibre%pos
    endif
    WRITE(MF,*) " @@@@@@@@@@@@@@@@@@@@ FIBRE @@@@@@@@@@@@@@@@@@@@"
    if(siam_index==0) then
       WRITE(MF,'(A11,4(I4,1x))') " DIRECTION ", M%DIR, &
            m%mag%parent_fibre%parent_layout%index,m%mag%parent_fibre%pos, &
            m%mag%parent_fibre%parent_layout%n
    else
       WRITE(MF,'(A11,4(I4,1x),A9,2(I4,1x))') " DIRECTION ", M%DIR, &
            m%mag%parent_fibre%parent_layout%index,m%mag%parent_fibre%pos, &
            m%mag%parent_fibre%parent_layout%n," Siamese ",siam_pos,siam_index
    endif
    CALL print_chart(m%CHART,mf)
    CALL print_PATCH(m%PATCH,mf)
    CALL print_element(M,M%MAG,mf)
    WRITE(MF,*) " @@@@@@@@@@@@@@@@@@@@  END  @@@@@@@@@@@@@@@@@@@@"

  END subroutine print_FIBRE

  subroutine READ_FIBRE(m,mf)
    implicit none
    integer mf,I
    type(FIBRE), pointer :: m
    character*255 line
    READ(MF,*) LINE
    READ(MF,'(A11,I4)') LINE(1:11),M%DIR
    CALL READ_chart(m%CHART,mf)
    CALL READ_PATCH(m%PATCH,mf)
    CALL READ_element(m,M%MAG,mf)
    READ(MF,*) LINE

  END subroutine READ_FIBRE

  subroutine READ_FIBRE_2_lines(mf,DIR,index,pos,n,siam_index,siam_pos)
    implicit none
    integer mf,I
    character*255 line

    integer DIR,index,pos,n,siam_index,siam_pos
    READ(MF,*) LINE
    siam_index=0
    siam_pos=0

    READ(MF,'(A11,4(I4,1x),A9,2(I4,1x))') LINE(1:11),DIR,index,pos,n, &
         LINE(12:20),siam_pos,siam_index
    !    CALL READ_chart(m%CHART,mf)
    !    CALL READ_PATCH(m%PATCH,mf)
    !    CALL READ_element(M%MAG,mf)
    !    READ(MF,*) LINE

  END subroutine READ_FIBRE_2_lines

  subroutine print_PATCH(m,mf)
    implicit none
    integer mf,I
    type(PATCH), pointer :: m
    character*255 line

    IF(IABS(M%PATCH)+iabs(M%energy)+iabs(M%time)/=0) then
       WRITE(MF,*) " >>>>>>>>>>>>>>>>>> PATCH <<<<<<<<<<<<<<<<<<"
       WRITE(MF,*) M%PATCH,M%ENERGY,M%TIME," patch,energy,time"
       WRITE(MF,*) M%A_X1,M%A_X2,M%B_X1,M%B_X2," discrete 180 rotations"
       WRITE(LINE,*) M%A_D,M%A_ANG,"  a_d, a_ang "
       WRITE(MF,'(A255)') LINE
       WRITE(LINE,*) M%B_D,M%B_ANG,"  b_d, b_ang "
       WRITE(MF,'(A255)') LINE
       WRITE(MF,*) M%A_T,M%B_T,"  time patches a_t and b_t "
       WRITE(MF,*) " >>>>>>>>>>>>>>>>>>  END  <<<<<<<<<<<<<<<<<<"
    else
       WRITE(MF,*) " NO PATCH "
    endif
  END subroutine print_PATCH

  subroutine READ_PATCH(m,mf)
    implicit none
    integer mf,I
    type(PATCH), pointer :: m
    character*255 line

    READ(MF,*)LINE
    if(index(line,"NO")==0) then
       READ(MF,*) M%PATCH,M%ENERGY,M%TIME
       READ(MF,*) M%A_X1,M%A_X2,M%B_X1,M%B_X2
       READ(MF,*) M%A_D,M%A_ANG
       READ(MF,*) M%B_D,M%B_ANG
       READ(MF,*) M%A_T,M%B_T
       READ(MF,*) LINE
    endif

  END subroutine READ_PATCH

  subroutine print_chart(m,mf)
    implicit none
    integer mf,I
    type(CHART), pointer :: m
    character*255 line
    real(dp) norm

    norm=zero
    do i=1,3
       norm=abs(M%D_IN(i))+norm
       norm=abs(M%ANG_IN(i))+norm
       norm=abs(M%ANG_OUT(i))+norm
       norm=abs(M%D_OUT(i))+norm
    enddo
    if(norm>zero.OR.print_frame) then
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
    integer mf,I
    type(CHART), pointer :: m
    character*60 line
    READ(mf,*) LINE
    if(index(line,"NO")==0) then
       CALL READ_magnet_frame(m%F,mf)
       READ(MF,*) M%D_IN,M%ANG_IN
       READ(MF,*) M%D_OUT,M%ANG_OUT
       READ(mf,*) LINE
    else
       do_survey=.true.
    endif
  end subroutine READ_chart


  subroutine READ_chart_fake(mf)
    implicit none
    integer mf,I
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
       do_survey=.true.
    endif
    call kill(f)
  end subroutine READ_chart_fake


  subroutine print_element(P,m,mf)
    implicit none
    integer mf,I
    type(FIBRE), pointer :: P
    type(element), pointer :: m
    character*255 line

    WRITE(MF,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ELEMENT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    if(m%vorname(1:1)==' ') then
       WRITE(MF,*) M%KIND,M%NAME, ' NOVORNAME'
    ELSE
       WRITE(MF,*) M%KIND,M%NAME,' ',M%VORNAME
    ENDIF
    WRITE(MF,*) M%L,M%PERMFRINGE,M%MIS , " L,PERMFRINGE,MIS "
    WRITE(LINE,*) M%FINT,M%HGAP,M%H1,M%H2, " FINT,HGAP,H1,H2 "
    WRITE(MF,'(A255)') LINE
    WRITE(LINE,*) 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, " no more mis"
    WRITE(MF,'(A255)') LINE
    IF(ASSOCIATED(M%FREQ)) THEN
       WRITE(MF,*) " CAVITY INFORMATION "
       WRITE(LINE,*) M%VOLT, M%FREQ,M%PHAS,M%DELTA_E,M%LAG,M%THIN, " VOLT,FREQ, PHAS, DELTA_E, LAG, THIN"
       WRITE(MF,'(A255)') LINE
    ELSEIF(ASSOCIATED(M%VOLT)) THEN
       WRITE(MF,*) " ELECTRIC SEPTUM INFORMATION "
       WRITE(MF,*) M%VOLT,M%PHAS, "VOLT, PHAS(rotation angle) "
    ELSE
       WRITE(MF,*) " NO ELECTRIC ELEMENT INFORMATION "
    ENDIF
    IF(ASSOCIATED(M%B_SOL)) THEN
       WRITE(MF,*)  " SOLENOID_PRESENT ",M%B_SOL, " B_SOL"
    ELSE
       WRITE(MF,*) " NO_SOLENOID_PRESENT ",zero
    ENDIF
    CALL print_magnet_chart(P,m%P,mf)
    IF(ASSOCIATED(M%an)) THEN
       do i=1,m%p%NMUL
          write(mf,*) m%bn(i),m%an(i), " BN AN ",I
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
    open(unit=mf,file=filename)

    nst=2*el%p%nst+1

    cl=(clight/c_1d8)
    BRHO=el%p%p0c*ten/cl

    call init(EL%B(1)%no,2)
    CALL ALLOC(B)

    write(mf,*) nst,el%p%ld,el%p%b0,EL%B(1)%no,.false.
    do i=1,nst
       B(1)=morph(1.d0.mono.1)
       B(2)=morph(1.d0.mono.2)
       B(3)=ZERO;
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
    if(.not.associated(el%internal)) allocate(el%internal(3))
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
    integer mf,n,i
    character*255 line
    real(dp) offset

    read(mf,'(a255)') line
    read(mf,*) n,offset
    call INIT_SAGAN_POINTERS(EL,N)
    el%offset=offset
    do i=1,n
       read(mf,*) el%a(i),el%f(i),EL%FORM(i),EL%K(1:3,i)
    enddo

    read(mf,'(a255)') line

  end subroutine read_undu_R


  subroutine print_specific_element(el,mf)
    implicit none
    type(element), pointer :: el
    integer mf
    character*255 line

    select case(el%kind)
    CASE(KIND0,KIND1,kind2,kind5,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17)
    case(kind3)
       WRITE(LINE,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle
       WRITE(MF,'(A255)') LINE
    case(kind4)
       WRITE(MF,*) el%c4%N_BESSEL
    case(kind10)
       WRITE(MF,*) el%tp10%DRIFTKICK
    case(kind16,kind20)
       WRITE(MF,*) el%k16%DRIFTKICK,el%k16%LIKEMAD, " driftkick,likemad"
    case(kind18)
       WRITE(MF,*) " RCOLLIMATOR HAS AN INTRINSIC APERTURE "
       CALL print_aperture(EL%RCOL18%A,mf)
    case(kind19)
       WRITE(MF,*) " ECOLLIMATOR HAS AN INTRINSIC APERTURE "
       CALL print_aperture(EL%ECOL19%A,mf)
    case(kind21)
       WRITE(MF,*) el%cav21%PSI,el%cav21%DPHAS,el%cav21%DVDS
    case(KINDWIGGLER)
       call print_wig(el%wi,mf)
    case(KINDpa)
       call print_pancake(el%pa,mf)
    case default
       write(6,*) " not supported in print_specific_element",el%kind
       stop 101
    end select

  end subroutine print_specific_element

  subroutine read_specific_element(el,mf)
    implicit none
    type(element), pointer :: el
    integer mf
    character*255 line

    select case(el%kind)
    CASE(KIND0,KIND1,kind2,kind5,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
    case(kind3)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(mf,*) el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle
    case(kind4)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%c4%N_BESSEL
    case(kind10)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%tp10%DRIFTKICK
    case(kind16,kind20)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%k16%DRIFTKICK,el%k16%LIKEMAD
    case(kind18)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       READ(MF,*) LINE
       CALL READ_aperture(EL%RCOL18%A,mf)
    case(kind19)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       READ(MF,*) LINE
       CALL READ_aperture(EL%ECOL19%A,mf)
    case(kind21)
       CALL SETFAMILY(EL)   ! POINTERS MUST BE ESTABLISHED BETWEEN GENERIC ELEMENT M AND SPECIFIC ELEMENTS
       read(MF,*) el%cav21%PSI,el%cav21%DPHAS,el%cav21%DVDS
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

    cl=(clight/c_1d8)
    BRHO=el%p%p0c*ten/cl


    call kanalnummer(mf)
    open(unit=mf,file=filename)
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
    integer mf,I
    type(fibre), pointer :: p
    type(element), pointer :: m
    character*120 line
    CHARACTER*21 SOL
    REAL(DP) B_SOL,r(3),d(3)

    READ(MF,*) LINE
    READ(MF,*) M%KIND,M%NAME,M%VORNAME
    CALL CONTEXT(M%NAME);
    CALL CONTEXT(M%VORNAME);
    IF(M%VORNAME(1:9)=='NOVORNAME') M%VORNAME=' '

    READ(MF,*) M%L,M%PERMFRINGE,M%MIS  !,M%EXACTMIS
    READ(MF,*) M%FINT,M%HGAP,M%H1,M%H2
    READ(MF,*) R,D
    READ(MF,*) LINE
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
    call read_specific_element(m,mf)

    READ(MF,*) LINE


  end subroutine READ_element


  subroutine READ_fake_element(p,mf)
    implicit none
    integer mf,I
    type(fibre),pointer ::  p
    type(fibre),pointer ::  f
    type(element),pointer ::  m
    character*120 line
    CHARACTER*21 SOL
    REAL(DP) B_SOL,r(3),d(3)
    nullify(m)
    call alloc(f)
    m=>f%mag
    !    m=0
    READ(MF,*) LINE
    READ(MF,*) M%KIND,M%NAME,M%VORNAME
    CALL CONTEXT(M%NAME);
    CALL CONTEXT(M%VORNAME);
    IF(M%VORNAME(1:9)=='NOVORNAME') M%VORNAME=' '

    READ(MF,*) M%L,M%PERMFRINGE,M%MIS   !,M%EXACTMIS
    READ(MF,*) M%FINT,M%HGAP,M%H1,M%H2
    READ(MF,*) R,D
    READ(MF,*) LINE
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
    call read_specific_element(m,mf)

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
    real(dp) BETA0,GAMMA0I, GAMBET,P0C


    READ(MF,*) LINE
    READ(MF,*) M%EXACT,M%METHOD,M%NST,M%NMUL
    READ(MF,'(A200)') LINE

    if(index(line,"TILT=")/=0) then
       READ(LINE,*) M%LD, M%LC, M%B0,til,M%tiltd
    else
       READ(LINE,*) M%LD, M%LC, M%B0
       M%tiltd=zero
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
       do_survey=.true.
    endif
  end subroutine read_magnet_frame

  subroutine print_aperture(m,mf)
    implicit none
    type(MADX_APERTURE), pointer :: m
    integer mf
    IF(.NOT.ASSOCIATED(M)) THEN
       write(mf,'(a20)') " NO MAGNET APERTURE "
    ELSE
       write(mf,'(a20)') "    MAGNET APERTURE "
       WRITE(MF,*) m%KIND   ! 1,2,3,4
       WRITE(MF,*) m%R
       WRITE(MF,*) m%X,m%Y
       write(mf,'(a23)')  " END OF MAGNET APERTURE"
    ENDIF

  end subroutine print_aperture


  subroutine READ_aperture(m,mf)
    implicit none
    type(MADX_APERTURE), pointer :: m
    integer mf
    character*120 line

    READ(mf,'(a120)') LINE

    CALL CONTEXT(LINE)

    IF(LINE(1:2)/='NO') THEN
       IF(.NOT.ASSOCIATED(M)) THEN
          CALL alloc(M)
       ENDIF

       READ(MF,*) m%KIND   ! 1,2,3,4
       READ(MF,*) m%R
       READ(MF,*) m%X,m%Y
       READ(mf,'(a120)') LINE
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
    integer mf,n,i,n_l,J,ncon1,ncon2,k,m,POS1,POS2
    type(MAD_UNIVERSE), TARGET :: U
    type(layout), pointer :: L,DNA_END
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0
    character*120 line
    type(FIBRE), pointer :: P

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
    integer I,mf,N,DIR,index_,pos,nt,kc,siam_index,siam_pos
    type(LAYOUT), TARGET :: L
    type(FIBRE), pointer :: P,current,siam
    LOGICAL(LP), OPTIONAL :: RING
    REAL(DP), OPTIONAL :: LMAX0
    LOGICAL(LP) RING_IT,doneit
    character*120 line
    real(dp) p0c,MASSF,LMAX0t
    type(internal_state) original
    integer,optional :: mf1


    RING_IT=MY_TRUE

    IF(PRESENT(RING)) RING_IT=RING

    if(present(mf1) ) then
       mf=mf1
    else
       call kanalnummer(mf)
       open(unit=mf,file=filename,status='OLD',err=2001)
    endif

    READ(MF,*) N,LMAX0t
    IF(PRESENT(LMAX0)) then
       if(LMAX0t/=zero) LMAX0=LMAX0T
    ENDIF
    read(MF,'(a120)') line
    call context(line)

    if(index(line,"FOR")/=0) then
       l%name=line(index(line,"FOR")+3:index(line,"FOR")+2+nlp)
    endif
    read(MF,*) MASSF,p0c
    read(MF,*) phase0,stoch_in_rec,initial_charge
    read(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING
    read(MF,*) SECTOR_NMUL_MAX,SECTOR_NMUL,OLD_IMPLEMENTATION_OF_SIXTRACK,HIGHEST_FRINGE
    read(MF,*) wedge_coeff
    read(MF,*) MAD8_WEDGE
    read(MF,'(a120)') line
    original=default
    if(allocated(s_b)) then
       firsttime_coef=.true.
       deallocate(s_b)
    endif
    !    L%MASS=MASSF
    MASSF=MASSF/pmae
    CALL MAKE_STATES(MASSF)
    default=original
    call Set_madx(p0c=p0c)
    DO I=1,N
       call  READ_FIBRE_2_lines(mf,dir,index_,pos,nt,siam_index,siam_pos)
       call  READ_chart_fake(mf)
       call move_to(L%DNA(index_)%L,p,pos)
       CALL APPEND_POINT(L,P)
       current=>l%end
       current%dir=dir

       if(siam_index/=0) then
          call move_to(L%DNA(siam_index)%L,siam,siam_pos)
          p%mag%siamese=>siam%mag
          write(6,*) p%mag%name,' is connected to ', siam%mag%name
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

    doneit=.true.
    call ring_l(L,doneit)
    if(do_survey) call survey(L)

    return

2001 continue

    Write(6,*) " File ",filename(1:len_trim(filename)) ," does not exist "

  END subroutine READ_pointed_INTO_VIRGIN_LAYOUT

!!!!  SIXTRACK FLAT FILE !!!!!

  subroutine PRINT_LAYOUT_SIXTRACT(L,FILENAME)
    implicit none
    character(*) filename
    integer I,MF,DI,nt
    type(LAYOUT), TARGET :: L
    type(FIBRE), pointer :: P
    character*255 line
    REAL(DP) D(3),ANG(3),PREC,DS
    LOGICAL(LP) ENERGY_PATCH

    !    P=>L%START
    !    i=1
    !    nt=0
    !    do while(i<=l%n)

    !       call count_PATCH_sixtrack(di,p)
    !       i=i+1+DI
    !       nt=nt+1
    !    ENDDO

    nt=0

    PREC=1.D-10
    call kanalnummer(mf)
    open(unit=mf,file=filename)

    WRITE(MF,*) L%N, " NUMBER OF FIBRES   "
    if(l%name(1:1)/=' ') then
       write(MF,'(a17,a16)') " GLOBAL DATA FOR ",l%name
    else
       write(MF,*) " $$$$$$$$$ GLOBAL DATA $$$$$$$$$"
    endif

    write(MF,*) l%start%mass,L%START%mag%p%p0c, "MASS, P0C"
    write(MF,*) phase0,stoch_in_rec,l%start%charge, " PHASE0, STOCH_IN_REC, CHARGE"
    write(MF,*) CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING, &
         "CAVITY_TOTALPATH,ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING"
    write(MF,*) " $$$$$$$ END GLOBAL DATA $$$$$$$$$$"

    if(.not.associated(L%t)) call make_node_layout(L)
    CALL FILL_SURVEY_DATA_IN_NODE_LAYOUT(L)

    P=>L%START
    i=1
    do while(i<=l%n)

       CALL print_FIBRE_SIXTRACK(P,mf)
       p=>p%next
       i=i+1    !+DI
    ENDDO


  END subroutine PRINT_LAYOUT_SIXTRACT

  subroutine print_FIBRE_SIXTRACK(m,mf)
    implicit none
    integer mf,I
    type(FIBRE), pointer :: m

    WRITE(MF,*) " @@@@@@@@@@@@@@@@@@@@ FIBRE @@@@@@@@@@@@@@@@@@@@"
    WRITE(MF,*)  M%DIR, " DIRECTION "
    CALL print_PATCH_MIS_SIXTRACK(M,MF)
    CALL print_element_SIXTRACK(M,M%MAG,mf)
    WRITE(MF,*) " @@@@@@@@@@@@@@@@@@@@  END  @@@@@@@@@@@@@@@@@@@@"

  END subroutine print_FIBRE_SIXTRACK


  subroutine print_element_SIXTRACK(P,m,mf)
    implicit none
    integer mf,I
    type(FIBRE), pointer :: P
    type(element), pointer :: m
    character*255 line

    WRITE(MF,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ELEMENT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    if(m%vorname(1:1)==' ') then
       WRITE(MF,*) M%KIND,M%NAME, ' NOVORNAME'
    ELSE
       WRITE(MF,*) M%KIND,M%NAME,' ',M%VORNAME
    ENDIF
    WRITE(MF,*) M%L
    IF(ASSOCIATED(M%FREQ)) THEN
       WRITE(MF,*) " CAVITY INFORMATION "
       !HALF*EL%P%DIR*EL%P%CHARGE*EL%volt*c_1d_3*SIN(twopi*EL%freq*x(6)/CLIGHT+EL%PHAS+EL%phase0)/EL%P%P0C

       WRITE(LINE,*) c_1d_3*M%VOLT, twopi*M%FREQ/CLIGHT,M%PHAS+M%C4%phase0
       WRITE(MF,'(A255)') LINE
    ELSEIF(ASSOCIATED(M%VOLT)) THEN
       WRITE(MF,*) " ELECTRIC SEPTUM INFORMATION "
       WRITE(MF,*) M%VOLT,M%PHAS, "VOLT, PHAS(rotation angle) "
    ELSE
       WRITE(MF,*) " NO ELECTRIC ELEMENT INFORMATION "
    ENDIF
    IF(ASSOCIATED(M%B_SOL)) THEN
       WRITE(MF,*)  M%B_SOL, " B_SOL"
    ELSE
       WRITE(MF,*) zero
    ENDIF
    CALL print_magnet_chart_SIXTRACK(P,m%P,mf)
    IF(ASSOCIATED(M%an)) THEN
       do i=1,m%p%NMUL
          write(mf,*) m%bn(i),m%an(i), " BN AN ",I
       enddo
    endif
    !    call print_specific_element(m,mf)
    WRITE(MF,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   END   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
  end subroutine print_element_SIXTRACK

  subroutine print_magnet_chart_SIXTRACK(P,m,mf)
    implicit none
    type(FIBRE), pointer :: P
    type(magnet_chart), pointer :: m
    integer mf
    character*200 line

    WRITE(MF,*) "MAGNET CHART MAGNET CHART MAGNET CHART MAGNET CHART MAGNET CHART MAGNET CHART "
    WRITE(MF,*) M%METHOD,M%NST,M%NMUL, " METHOD NST NMUL"
    WRITE(line,*) M%B0, " B0"
    WRITE(MF,'(A200)') LINE
    WRITE(LINE,*) P%BETA0,P%GAMMA0I, P%GAMBET, M%P0C,m%tiltd, " BETA0 GAMMA0I GAMBET P0C TILT"
    WRITE(MF,'(A200)') LINE

    write(mf,'(a68)') "END MAGNET CHART END MAGNET CHART END MAGNET CHART END MAGNET CHART "
  end subroutine print_magnet_chart_SIXTRACK

  subroutine print_PATCH_MIS_SIXTRACK(m,mf)
    implicit none
    type(FIBRE), pointer :: m
    integer mf
    character*200 line

    write(mf,'(a10)')  " GEOMETRY "
    WRITE(MF,*) M%PATCH%A_X1,M%PATCH%A_X2,M%PATCH%B_X1,M%PATCH%B_X2
    WRITE(MF,*) M%PATCH%A_T, M%PATCH%B_T
    WRITE(LINE,*) M%PATCH%A_ANG,M%PATCH%A_D
    WRITE(MF,'(a200)') LINE
    WRITE(LINE,*) M%PATCH%B_ANG,M%PATCH%B_D
    WRITE(MF,'(a200)') LINE
    WRITE(LINE,*) M%CHART%ANG_IN,M%CHART%D_IN
    WRITE(MF,'(a200)') LINE
    WRITE(LINE,*) M%CHART%ANG_OUT,M%CHART%D_OUT
    WRITE(MF,'(a200)') LINE
    WRITE(LINE,*) M%MAG%P%TILTD

    write(mf,'(a14)') " END GEOMETRY "
  end subroutine print_PATCH_MIS_SIXTRACK


end module madx_keywords
