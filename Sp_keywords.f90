!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90

module madx_keywords
  use S_fitting
  implicit none
  logical(lp)::mad8=my_false

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
    type(fibre), intent(inout)::el
    logical(lp), optional :: magnet_only
    type(keywords) key
    type(el_list) blank
    character*255 magnet
    character*17 MODEL
    INTEGER EXCEPTION  !,NSTD0,METD0
    LOGICAL(LP) EXACT0,magnet0
    logical(lp) FIBRE_flip0,MAD0
    logical(lp) :: t=.true.,f=.false.
    INTEGER FIBRE_DIR0,IL
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
       BLANK=TRUERECTILT_MADX(KEY%LIST%NAME,tilt.is.KEY%tiltd,KEY%LIST)
    CASE("RBEND         ")
       IF(KEY%MAD8) THEN
          ! MADX does some massaging of e1 and e2 for compatibility with mad8
          ! PTC stays with MAD8
          IF(KEY%LIST%B0/=ZERO) THEN
             KEY%LIST%L=   KEY%LIST%B0*KEY%LIST%L /(two*SIN(KEY%LIST%B0/two))    !  ld is computed
          ENDIF
          KEY%LIST%T1=KEY%LIST%T1+KEY%LIST%B0/two
          KEY%LIST%T2=KEY%LIST%T2+KEY%LIST%B0/two
          BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       ELSE
          BLANK=RECTTILT_MADX(KEY%LIST%NAME,tilt.is.KEY%tiltd,KEY%LIST)
       ENDIF
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
    CASE("USER_1         ")
       BLANK=USER_1(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("USER_2         ")
       BLANK=USER_2(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("WIGGLER        ")
       BLANK=WIGGLER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("TAYLORMAP      ")
       IF(KEY%LIST%file/=' '.and.KEY%LIST%file_rev/=' ') THEN
          BLANK=TAYLOR_MAP(KEY%LIST%NAME,FILE=KEY%LIST%file,FILE_REV=KEY%LIST%file_REV,t=tilt.is.KEY%tiltd)
       ELSEIF(KEY%LIST%file/=' '.and.KEY%LIST%file_rev==' ') THEN
          BLANK=TAYLOR_MAP(KEY%LIST%NAME,FILE=KEY%LIST%file,t=tilt.is.KEY%tiltd)
       ELSEIF(KEY%LIST%file==' '.and.KEY%LIST%file_rev/=' ') THEN
          BLANK=TAYLOR_MAP(KEY%LIST%NAME,FILE_REV=KEY%LIST%file_REV,t=tilt.is.KEY%tiltd)
       ELSE
          BLANK=TAYLOR_MAP(KEY%LIST%NAME,t=tilt.is.KEY%tiltd)
       ENDIF
    CASE DEFAULT
       WRITE(6,*) " "
       WRITE(6,*) " THE MAGNET"
       WRITE(6,*) " "
       WRITE(6,*) "  --->   ",MAGNET(1:IL)
       WRITE(6,*) " "
       WRITE(6,*)  " IS NOT PERMITTED "
       STOP 666
    END SELECT

    CALL EL_Q_FOR_MADX(EL,BLANK)

    CALL SET_MADX_(f,f)

    !    NSTD=NSTD0
    !    METD=METD0
    EXACT_MODEL=EXACT0
    FIBRE_FLIP= FIBRE_FLIP0
    FIBRE_DIR=FIBRE_DIR0
    MAD=MAD0
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

  !  Printing universes

  subroutine print_universe(universe,style,time_patch_i,print_mis)
    implicit none
    type(mad_universe) universe
    INTEGER N,J,I,i_tot,i_L,POS, NU,nb,style,dirb,dira
    TYPE(LAYOUT), POINTER :: REC,rec1
    TYPE(FIBRE), POINTER :: P
    TYPE(MAGNET_CHART), POINTER :: MC
    REAL(DP) MASS,BRHO
    CHARACTER*11 FILER
    CHARACTER*100 NUT,nut0
    INTEGER MINIMAL,MINIMAL_RICH_AT_PATCHES,ROOT_ONLY,ROOT_PLUS_PATCHES
    LOGICAL(LP) FRAME,PATCH,TIME_PATCH,time_patch_i,check,printmis
    LOGICAL(LP), optional :: print_mis

    printmis=.false.

    if(present(print_mis)) printmis=print_mis

    TIME_PATCH=TIME_PATCH_i

    CALL GET_ONE(MASS=MASS)

    MINIMAL=0     !  PATCHES ONLY AT SKELETON MALFORMATION   ELSEWHERE      MAD8 SURVEY ASSUMED
    MINIMAL_RICH_AT_PATCHES=1  !  PATCHES AND NECESSARY FRAMES T SKELETON MALFORMATION  ELSEWHERE MAD8 IN BETWEEN
    ROOT_ONLY=2                !  NO PATCHES, ONLY FRAMES  EVERYWHERE
    ROOT_PLUS_PATCHES=3   !  EVERYTHING!

    FRAME=.FALSE.
    IF(STYLE>=1) FRAME=.TRUE.
    PATCH=.TRUE.
    IF(ROOT_ONLY/=2) PATCH=.FALSE.

    FILER='real_system'

    open(unit=40,file='real_magnet.txt')
    rec=>universe%start
    n=0

    write(40,195) " CHARGE = ",1," MASS = ",MASS

    do j=1,universe%SHARED

       p=>rec%start; MC=>P%MAG%P;
       do i=1,rec%n
          n=n+1
          write(40,210) 'Magnet # = ', n,' for etienne: layout # = ',j,' position in layout = ',i,' name = ',P%MAG%NAME
          write(40,211) "SKELETON ARC LENGTH = ",   MC%LD,' ANGLE = ',MC%LD*MC%B0,' ROLL = ',MC%TILTD,  ' P0C = ',MC%P0C
          if(p%mag%mis) then
             write(40,209) " MISALIGNEMENTS D AND R = ",p%mag%D,p%mag%R
          else
             write(40,*) " MISALIGNEMENTS NO PRESENT "
          endif
          write(40,*) " "
          p=>p%next;MC=>P%MAG%P;
       enddo
       rec=>rec%next

    enddo

    close(40)

209 FORMAT(A16,6(1x,E15.8))
210 FORMAT(a11,i4,A25,i4,a22,i4,a8,A16)
    !211 FORMAT(A40,E15.8 ,1x,E15.8 ,1x,E15.8 ,1x,E15.8 )
211 FORMAT(A22,E15.8 ,A9,E15.8 ,A8,E15.8, A7,E15.8 )

204 FORMAT(A15,i2)
202 FORMAT(A15,I2)
200 FORMAT(A16,3(E15.8,1X))
201 FORMAT(A50,2(I4,1X))
203 FORMAT(A23,E15.8)
199 FORMAT(A6,I4)
197 FORMAT(A16,I4,A24,A16)
198 FORMAT(A59,3(1X,I4))
196 FORMAT(A30,I4,A13,I4,A13,I4 )
195 FORMAT(A10,I4,A8,E15.8)
194 FORMAT(A36,3(1X,I2))


    WRITE(nut0,'(I4,a6,i4)') UNIVERSE%N-UNIVERSE%SHARED,'_style',style

    nb=0
    nut=' '
    do i=1,LEN_TRIM(NUT0)
       if(nut0(i:i)/=' ') then
          nb=nb+1
          nut(nb:nb)=nut0(i:i)
       endif
    enddo

    rec1=>rec
    check=.true.

    DO NU=1,UNIVERSE%N-UNIVERSE%SHARED
       p=>rec1%start;
       do i=1,rec1%n
          if(p%patch%patch+p%patch%energy/=0) then
             check=.false.
             exit
          endif
          p=>p%next
       enddo
       if(.not.check) exit
       rec1=>rec1%next
    ENDDO



    DO NU=1,UNIVERSE%N-UNIVERSE%SHARED

       open(unit=40,file=FILER//NUT(1:NB)//'.txt')
       WRITE(40,*) " STYLE OF FILE = ",STYLE
       p=>rec%start;
       IF(TIME_PATCH.or.check) THEN
          write(40,'(a72)') "THIS BEAM LINE    --> CAN <--         BE TRACKED WITH RELATIVE TIME OR PATH LENGTH"
          write(40,'(a72)') " BECAUSE THE TIME PATCHES ARE EITHER SET OR NOT NEEDED (STANDARD MAD8 SKELETON)   "
       ELSE
          write(40,'(a72)') "THIS BEAM LINE    --> CANNOT <--      BE TRACKED WITH RELATIVE TIME OR PATH LENGTH"
          write(40,'(a72)') " BECAUSE THE TIME PATCHES ARE NOT SET AND ARE NEEDED. (MAD8 SKELETON IS DEFORMED) "
       ENDIF
       write(40,*) "  "
       write(40,'(a40)') "POSITION OF THE ORIGIN OF THIS BEAM LINE"
       CALL PRINT_INITIAL_FRAME(P,40)

       write(40,'(a90)') "******************************************************************************************"
       write(40,'(a90)') "******************************************************************************************"
       write(40,'(a90)') "******************************************************************************************"


       do i=1,rec%n
          write(40,*) " "
          write(40,197) " FIBRE NUMBER = ", i," NAME OF MAGNET IN IT = ",p%mag%name
          call locate_in_universe(P,i_tot,i_L,POS)
          write(40,196)  " POSITION IN REAL MAGNET FILE ", i_tot, " BEAM LINE # ", I_L, " MAGNET # ",POS

          dirb=0;dira=0;
          if(associated(P%PREVIOUS)) dirb=P%PREVIOUS%DIR
          if(associated(P%NEXT))     dira=P%NEXT%DIR
          WRITE(40,198) " PROPAGATION DIRECTION OF PREVIOUS, PRESENT AND NEXT FIBRE ",dirb,P%DIR,dira

          write(40,194) " PATCHES FLAGS (ENERGY,TIME,FRAMES) ", p%patch%energy,p%patch%TIME,p%patch%patch


          if(p%patch%energy/=0) WRITE(40,204)   " ENERGY FLAG = ", P%PATCH%ENERGY
          ! TIME PATCHES
          select case(p%patch%TIME)
          case(1)
             IF(TIME_PATCH) THEN
                WRITE(40,204) "  TIME FLAG  = ",P%PATCH%TIME
                WRITE(40,203) " TIME PATCH ENTRANCE = ",P%PATCH%A_T
             ENDIF

          case(2)
             IF(TIME_PATCH) THEN
                WRITE(40,204) "  TIME FLAG  = ",P%PATCH%TIME
                WRITE(40,203) " TIME PATCH EXIT     = ",P%PATCH%B_T
             ENDIF
          case(3)
             IF(TIME_PATCH) THEN
                WRITE(40,204) "  TIME FLAG  = ",P%PATCH%TIME
                WRITE(40,203) " TIME PATCH ENTRANCE = ",P%PATCH%A_T
                WRITE(40,203) " TIME PATCH EXIT     = ",P%PATCH%B_T
             ENDIF
          END SELECT



          ! GEOMETRICAL PATCHES
          select case(p%patch%patch)
          CASE(0)

             IF(STYLE==ROOT_ONLY.OR.STYLE==ROOT_PLUS_PATCHES) THEN
                write(40,'(a90)')    "##########################################################################################"
                CALL PRINT_INITIAL_FRAME(P,40)
                write(40,'(a90)')    "##########################################################################################"
             ENDIF
          case(1)
             IF(FRAME) THEN
                write(40,*) "**************** EXIT FRAME OF THE PREVIOUS FIBRE ********************"
                CALL PRINT_INITIAL_FRAME(P%PREVIOUS,40)
                write(40,*) "**********************************************************************"
                write(40,*) "*********************** FRAME OF THE FIBRE ***************************"
                CALL PRINT_INITIAL_FRAME(P,40)
                write(40,*) "**********************************************************************"
             ENDIF
             IF(PATCH) THEN
                write(40,*)  "$$$$$$$$$$$$$$$$$$$ ALMOST MINIMALIST TYPE OF INFORMATION $$$$$$$$$$$$$$$$$$"
                WRITE(40,*)  " FRONTAL PATCH NEEDED"
                WRITE(40,200)   " TRANSLATION = ",P%PATCH%A_D
                WRITE(40,200)   "   ROTATION  = ",P%PATCH%A_ANG
                WRITE(40,201)   " DIRECTIONAL 180 DEGREES ROTATIONS X AND Y AXIS = ", P%PATCH%A_YZ,P%PATCH%A_XZ
                write(40,*)  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
             ENDIF
          case(2)
             IF(FRAME) THEN
                write(40,*) "*********************** FRAME OF THE FIBRE ***************************"
                CALL PRINT_INITIAL_FRAME(P,40)
                write(40,*) "**********************************************************************"
                write(40,*) "**************** ENTRANCE FRAME OF THE NEXT FIBRE ********************"
                CALL PRINT_INITIAL_FRAME(P%NEXT,40)
                write(40,*) "**********************************************************************"
             ENDIF
             IF(PATCH) THEN
                write(40,*)  "$$$$$$$$$$$$$$$$$$$ ALMOST MINIMALIST TYPE OF INFORMATION $$$$$$$$$$$$$$$$$$"
                WRITE(40,*)  " EXIT PATCH NEEDED"
                WRITE(40,200)   " TRANSLATION = ",P%PATCH%B_D
                WRITE(40,200)   "   ROTATION  = ",P%PATCH%B_ANG
                WRITE(40,201)   " DIRECTIONAL 180 DEGREES ROTATIONS X AND Y AXIS = ", P%PATCH%B_YZ,P%PATCH%B_XZ
                write(40,*)  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
             ENDIF
          case(3)
             IF(FRAME) THEN
                write(40,*) "**************** ENTRANCE FRAME OF THE NEXT FIBRE ********************"
                CALL PRINT_INITIAL_FRAME(P%PREVIOUS,40)
                write(40,*) "*********************** FRAME OF THE FIBRE ***************************"
                CALL PRINT_INITIAL_FRAME(P,40)
                write(40,*) "**************** ENTRANCE FRAME OF THE NEXT FIBRE ********************"
                CALL PRINT_INITIAL_FRAME(P%NEXT,40)
                write(40,*) "**********************************************************************"
             ENDIF
             IF(PATCH) THEN
                write(40,*)  "$$$$$$$$$$$$$$$$$$$ ALMOST MINIMALIST TYPE OF INFORMATION $$$$$$$$$$$$$$$$$$"
                WRITE(40,*)  " FRONTAL PATCH NEEDED"
                WRITE(40,200)   " TRANSLATION = ",P%PATCH%A_D
                WRITE(40,200)   "   ROTATION  = ",P%PATCH%A_ANG
                WRITE(40,201)   " DIRECTIONAL 180 DEGREES ROTATIONS X AND Y AXIS = ", P%PATCH%A_YZ,P%PATCH%A_XZ
                WRITE(40,*)  " EXIT PATCH NEEDED"
                WRITE(40,200)   " TRANSLATION = ",P%PATCH%B_D
                WRITE(40,200)   "   ROTATION  = ",P%PATCH%B_ANG
                WRITE(40,201)   " DIRECTIONAL 180 DEGREES ROTATIONS X AND Y AXIS = ", P%PATCH%B_YZ,P%PATCH%B_XZ
                write(40,*)  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
             ENDIF
          END SELECT
          if(printmis) then
             if(P%mag%mis) then
                CALL PRINT_INITIAL_FRAME_mag(p,40)
             endif
          endif




          p=>p%next
       enddo
       close(40)

       rec=>rec%next
    ENDDO

  END subroutine print_universe


  SUBROUTINE PRINT_INITIAL_FRAME(F,MF)
    implicit none
    INTEGER MF
    TYPE(FIBRE) F
    IF(F%DIR==1) THEN
       WRITE(MF,*) F%CHART%F%A
       WRITE(MF,*) F%CHART%F%ENT
    ELSE
       WRITE(MF,*) F%CHART%F%B
       WRITE(MF,*) F%CHART%F%EXI
    ENDIF

  END SUBROUTINE PRINT_INITIAL_FRAME

  SUBROUTINE PRINT_INITIAL_FRAME_mag(F,MF)
    implicit none
    INTEGER MF
    TYPE(FIBRE) F

209 FORMAT(A25,6(1x,E15.8))
    Write(mf,*) " &&&&&&&&&&&&&   Misalignment Information  &&&&&&&&&&&&&  "
    Write(mf,*) "  See Routine MIS_FIBR for actual ordering of Euclidean operators  "
    IF(F%DIR==1) THEN
       write(mf,209) " d(3),angle(3) to front = ",f%CHART%d_IN,f%CHART%ANG_IN
       Write(mf,*) " Front Frame "
       WRITE(MF,*) F%mag%p%F%A
       WRITE(MF,*) F%mag%p%F%ENT
       Write(mf,*) " Back Frame "
       WRITE(MF,*) F%mag%p%F%B
       WRITE(MF,*) F%mag%p%F%EXI
       write(mf,209) " d(3),angle(3) to front = ",f%CHART%d_out,f%CHART%ANG_out
    ELSE
       write(mf,209) " d(3),angle(3) to front = ",f%CHART%d_out,f%CHART%ANG_out
       Write(mf,*) " Back Frame (Reverse Propagator)"
       WRITE(MF,*) F%mag%p%F%B
       WRITE(MF,*) F%mag%p%F%EXI
       Write(mf,*) " Front Frame (Reverse Propagator) "
       WRITE(MF,*) F%mag%p%F%A
       WRITE(MF,*) F%mag%p%F%ENT
       write(mf,209) " d(3),angle(3) to front = ",f%CHART%d_IN,f%CHART%ANG_IN
    ENDIF
    Write(mf,*) " &&&&&&&&&&&&&   End Misalignment Information  &&&&&&&&&&&&&  "

  END SUBROUTINE PRINT_INITIAL_FRAME_mag

  SUBROUTINE SURVEY_nikolai(basis,OMEGA,b0,ld,roll,DIR,ent,a,mid,o,exi,b)
    ! SURVEYS A SINGLE ELEMENT FILLS IN CHART AND MAGNET_CHART; LOCATES ORIGIN AT THE ENTRANCE OR EXIT
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::DIR
    REAL(DP) ENT(3,3),EXI(3,3),HA,D(3),BASIS(3,3),OMEGA(3),Ang(3),N(3),out(3,3),a(3),b(3),mid(3,3),o(3)
    real(dp) angle,roll,lc,ld,b0

    INTEGER I,J


    angle=ld*b0
    HA=DIR*angle/TWO
    if(angle==zero) then
       lc=ld
    else
       lc=two*ld/angle*sin(angle/two)
    endif

    D=ZERO
    D(3)=DIR*LC/TWO
    IF(DIR==1) THEN
       Ang=ZERO;Ang(3)=roll  ;
       CALL GEO_ROT(basis,ENT    ,Ang  ,basis)
       Ang=ZERO;Ang(2)=HA ;
       CALL GEO_ROT(ENT     ,MID ,Ang     ,ENT)
       CALL GEO_ROT(MID,EXI     , Ang     ,MID)

       Ang=ZERO;Ang(3)=-roll ;
       CALL GEO_ROT(EXI     ,out ,Ang,EXI)

       O=OMEGA
       CALL GEO_TRA(O,MID,D,1)
       B=O
       CALL GEO_TRA(B,MID,D,1)

       a=OMEGA
       ent=basis
       exi=out
    ELSE
       Ang=ZERO;Ang(3)=roll  ;
       CALL GEO_ROT(basis,EXI      ,Ang  ,basis)
       Ang=ZERO;Ang(2)=HA ;
       CALL GEO_ROT(EXI     ,MID ,Ang     ,EXI)
       CALL GEO_ROT(MID,ENT     , Ang     ,MID)

       Ang=ZERO;A(3)=-roll ;
       CALL GEO_ROT(ENT     ,out ,Ang ,ENT)

       O=OMEGA
       CALL GEO_TRA(O,MID,D,1)
       A=O
       CALL GEO_TRA(A,MID,D,1)

       b=OMEGA
       exi=basis
       ent=out


    ENDIf

  END SUBROUTINE SURVEY_nikolai


end module madx_keywords
