module madx_keywords
  use S_fitting
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
    logical(lp) :: doneitt=.true.
    logical(lp) :: doneitf=.false.
    type(fibre), intent(inout)::el
    logical(lp), optional :: magnet_only
    type(keywords) key
    type(el_list) blank
    character*255 magnet
    character*17 MODEL
    INTEGER EXCEPTION,NSTD0,METD0
    LOGICAL(LP) EXACT0,magnet0
    logical(lp) FIBRE_flip0,MAD0
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

    CALL SET_MADX_(doneitt,magnet0)

    select case(MODEL)
    CASE("DRIFT_KICK       ")
       MADTHICK=drift_kick_drift
    CASE("MATRIX_KICK      ")
       MADTHICK=matrix_kick_matrix
    CASE("DELTA_MATRIX_KICK")
       MADTHICK=kick_sixtrack_kick
    CASE("BEND_KICK        ")
       MADTHICK=matrix_kick_matrix
    CASE DEFAULT
       EXCEPTION=1
       ipause=mypause(444)
       RETURN
    END SELECT


    NSTD0=NSTD
    METD0=METD
    EXACT0=EXACT_MODEL
    FIBRE_FLIP0= FIBRE_FLIP
    FIBRE_DIR0=FIBRE_DIR
    MAD0=MAD

    NSTD=KEY%NSTEP
    METD=KEY%METHOD
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

    CALL SET_MADX_(doneitf,doneitf)

    NSTD=NSTD0
    METD=METD0
    EXACT_MODEL=EXACT0
    FIBRE_FLIP= FIBRE_FLIP0
    FIBRE_DIR=FIBRE_DIR0
    MAD=MAD0
  end subroutine create_fibre

  subroutine zero_key(key)
    implicit none
    type(keywords) , intent(out):: key
    real h
    key%magnet="DRIFT"
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
    key%mad8 = .false.
    key%tiltd=ZERO
    key%LIST=0

  end subroutine zero_key

  ! MAD-X SURVEY ROUTINES AND INTERFACES TO PTC
  ! THESE ROUTINES DO NOT REQUIRE PTC LATTICES
  ! THEY REPLACED OLD MAD8 SHIT

  subroutine zero_MAD_SURVEY(FRAME)
    implicit none
    type(MADX_SURVEY) , intent(out):: FRAME

    FRAME%LD=ZERO
    FRAME%TILT=ZERO
    FRAME%ALPHA=ZERO
    FRAME%PSI=ZERO
    FRAME%PHI=ZERO
    FRAME%THETA=ZERO
    FRAME%CHART=1
  end subroutine zero_MAD_SURVEY


  subroutine drift_bending(el,F) ! this turns a drift into a bend for survey purposes
    implicit none
    type(fibre), intent(inout)::el
    TYPE(MADX_SURVEY), intent(in):: F

    ! THIS MUST BE USED ONLY IN MAD-X SURVEYS.
    ! TURN ON with_external_frame  TEMPORARILY IF OFF
    ! TURN ON with_internal_frame  TEMPORARILY IF OFF
    ! TURN ON with_chart           TEMPORARILY IF OFF
    ! TURN ON with_patch           TEMPORARILY IF OFF
    !
    EL%MAG%P%LD=F%LD
    EL%MAG%P%TILTD=F%TILT

    IF(F%ALPHA/=ZERO) THEN
       EL%MAG%P%LC=TWO*EL%MAG%P%LD*SIN(F%ALPHA/TWO)/F%ALPHA
       EL%MAG%L=EL%MAG%P%LC    ! NOT REALLY NECESSARY
       EL%MAG%P%B0=F%ALPHA/el%MAG%P%LD
    ELSE
       EL%MAG%L=EL%MAG%P%LD
       EL%MAG%P%LC=EL%MAG%P%LD
       EL%MAG%P%B0=ZERO
    ENDIF
    CALL SURVEY(EL)
  end subroutine drift_bending

  subroutine SURVEY_LEAP_FROG(D1,D2,F)
    IMPLICIT NONE
    type(fibre), intent(inout)::D1,D2
    TYPE(MADX_SURVEY), intent(inout)::F

    CALL drift_bending(D2,F)
    CALL ROTATE_FIBRE(D2,D1)
    CALL COPY(D2%MAG,D1%MAG)
    CALL COPY(D2%CHART,D1%CHART)
    CALL EXTRACT_MADX_FRAME(D1,F)


  END subroutine SURVEY_LEAP_FROG

  subroutine EXTRACT_MADX_FRAME(D1,F)
    IMPLICIT NONE
    type(fibre), intent(in)::D1
    TYPE(MADX_SURVEY), intent(inout)::F
    REAL(DP) AX,AY,AZ

    ! IMPORTANT LINE HERE, PTC ONLY CARES ABOUT FRAMES
    ! SINCE MEMORY IS CHEAP IN AD 2002
    F%CHART=D1%CHART

    ! MAD8 BULSHIT FOLLOWS I.E. (PHI,THETA,PSI)
    !COMPUTED BETWEEN -PI AND PI IN THE FORTRAN CONVENTION
    ! SAVING MEMORY I SUPPOSE IN MAD8

    IF(F%CHART%F%EXI(3,2)>ONE) F%CHART%F%EXI(3,2)=ONE
    IF(F%CHART%F%EXI(3,2)<-ONE) F%CHART%F%EXI(3,2)=-ONE
    F%PHI=ASIN(F%CHART%F%EXI(3,2))
    AX=F%CHART%F%EXI(3,1); AZ=F%CHART%F%EXI(3,3);
    F%THETA = ATAN2(AX,AZ)
    AX=F%CHART%F%EXI(1,2); AY=F%CHART%F%EXI(2,2);
    F%PSI = ATAN2(AX,AY)

  END subroutine EXTRACT_MADX_FRAME



end module madx_keywords
