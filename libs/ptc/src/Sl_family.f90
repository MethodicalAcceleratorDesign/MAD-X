!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

MODULE S_FAMILY
  USE S_FIBRE_BUNDLE
  IMPLICIT NONE
  public

  ! LINKED LIST
  PRIVATE SURVEY_EXIST_PLANAR_L_NEW ,SURVEY_EXIST_PLANAR_IJ,MISALIGN_FIBRE_EQUAL ,SURVEY_EXIST_PLANAR_I
  !,SURVEY_NO_PATCH
  PRIVATE COPY_LAYOUT,COPY_LAYOUT_I,KILL_PARA_L
  PRIVATE FIBRE_WORK,FIBRE_POL,FIBRE_BL,ADDP_ANBN,WORK_FIBRE,BL_FIBRE,layout_WORK
  PRIVATE COPY_LAYOUT_IJ,PUT_APERTURE_FIB,REMOVE_APERTURE_FIB,PUT_APERTURE_FIBt
  private copy_fibre
  ! old Sj_elements
  PRIVATE SURVEY_mag  !,TRANSLATE_magnet,rot_magnet
  !END old Sj_elements

  INTERFACE EL_TO_ELP
     !LINKED
     MODULE PROCEDURE EL_TO_ELP_L
  END INTERFACE

  INTERFACE ELP_TO_EL
     !LINKED
     MODULE PROCEDURE ELP_TO_EL_L
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_fibre
  END INTERFACE

  INTERFACE SURVEY
     ! LINK LIST
     !     MODULE PROCEDURE SURVEY_NO_PATCH            ! NO PATCH
     MODULE PROCEDURE SURVEY_EXIST_PLANAR_L_NEW  ! ORDINARY SURVEY STARTING AT POSITION 1
     MODULE PROCEDURE SURVEY_FIBRE
     MODULE PROCEDURE SURVEY_EXIST_PLANAR_IJ
     MODULE PROCEDURE SURVEY_EXIST_PLANAR_I
     ! old Sj_elements
     !     MODULE PROCEDURE SURVEY_mag  ! tracks a chart for survey
     !END old Sj_elements
  END INTERFACE

  !  INTERFACE TRACK
  ! LINK LIST
  !     MODULE PROCEDURE SURVEY_FIB
  !     MODULE PROCEDURE SURVEY_EXIST_PLANAR_IJ
  !     MODULE PROCEDURE SURVEY_EXIST_PLANAR_I
  ! old Sj_elements
  !     MODULE PROCEDURE SURVEY_mag  ! tracks a chart for survey
  !END old Sj_elements
  !  END INTERFACE

  INTERFACE KILL_PARA
     MODULE PROCEDURE KILL_PARA_L
  END INTERFACE

  INTERFACE ADD
     MODULE PROCEDURE ADDP_ANBN
  END INTERFACE

  INTERFACE PUT_APERTURE
     MODULE PROCEDURE PUT_APERTURE_FIB                               ! NEED UPGRADE
     MODULE PROCEDURE PUT_APERTURE_FIBt                               ! NEED UPGRADE
  END  INTERFACE

  INTERFACE REMOVE_APERTURE
     MODULE PROCEDURE REMOVE_APERTURE_FIB                               ! NEED UPGRADE
  END  INTERFACE



  INTERFACE ASSIGNMENT (=)
     ! LINKED
     MODULE PROCEDURE SCAN_FOR_POLYMORPHS
     MODULE PROCEDURE FIBRE_WORK
     MODULE PROCEDURE MISALIGN_FIBRE_EQUAL
     MODULE PROCEDURE FIBRE_POL
     MODULE PROCEDURE FIBRE_BL
     MODULE PROCEDURE BL_FIBRE
     MODULE PROCEDURE WORK_FIBRE
     MODULE PROCEDURE layout_WORK
  END  INTERFACE

  INTERFACE EQUAL
     ! LINKED
     MODULE PROCEDURE COPY_LAYOUT
  END  INTERFACE

  INTERFACE COPY
     ! LINKED
     MODULE PROCEDURE COPY_LAYOUT_I
     MODULE PROCEDURE COPY_LAYOUT_IJ
  END  INTERFACE

  INTERFACE TRANS
     ! LINKED
     MODULE PROCEDURE TRANSLATE_layout
     MODULE PROCEDURE TRANSLATE_fibre
     MODULE PROCEDURE TRANSLATE_frame
  END  INTERFACE

  INTERFACE TRANSLATE
     ! LINKED
     MODULE PROCEDURE TRANSLATE_layout
     MODULE PROCEDURE TRANSLATE_fibre
     MODULE PROCEDURE TRANSLATE_frame
     !     MODULE PROCEDURE TRANSLATE_magnet   !element input
  END  INTERFACE

  INTERFACE ROTATE
     ! LINKED
     MODULE PROCEDURE ROTATE_LAYOUT
     MODULE PROCEDURE ROTATE_FIBRE
     MODULE PROCEDURE ROTATE_FRAME
     MODULE PROCEDURE rotate_magnet
  END  INTERFACE

  INTERFACE ROTATION
     ! LINKED
     MODULE PROCEDURE ROTATE_LAYOUT
     MODULE PROCEDURE ROTATE_FIBRE
     MODULE PROCEDURE ROTATE_FRAME
  END  INTERFACE



CONTAINS

  ! old Sj_elements

  SUBROUTINE SURVEY_mag(C,el,dir,magnetframe,e_in) !  Tracks the chart through a magnet
    IMPLICIT NONE
    TYPE(CHART), TARGET ,optional, INTENT(INOUT):: C
    type(element),target,intent(inout) :: el
    !    TYPE(magnet_frame), OPTIONAL :: magnetframe
    TYPE(magnet_frame),target, OPTIONAL :: magnetframe
    INTEGER, intent(in):: dir
    TYPE(INNER_FRAME), OPTIONAL :: E_IN

    !  All PTC magnet have the same convention for the internal frame
    !  Show a user want to add a magnet with corckscrew survey
    !  his survey would have to be "caught" in this interface

    SELECT CASE(EL%KIND)
    case(kind0:kind22,KINDWIGGLER,kindpa,kindsuperdrift,kindabell)
       call SURVEY_chart(C,el%p,dir,magnetframe,E_IN)

       !    case(kind23)
       !       call SURVEY_CHART_layout(C,el,DIR,MAGNETFRAME,E_IN)
    case default
       write(6,*) el%kind," not supported SURVEY_mag in S_FAMILY"
    END SELECT

    ! RECURSIVE   SUBROUTINE SURVEY_EXIST_PLANAR_L_NEW(PLAN,ENT,A) ! CALLS ABOVE ROUTINE FROM FIBRE #1 TO #PLAN%N : STANDARD SURVEY

  end SUBROUTINE SURVEY_mag



  ! END old Sj_elements

  !NEW

  SUBROUTINE locate_mid_frame(R,mid,o,ld)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    TYPE(FIBRE), POINTER:: P
    real(dp), intent(out) :: mid(3,3),o(3),ld
    integer i


    if(mod(R%N,2)==0) then
       P=>R%START
       DO I=1,R%N/2
          P=>P%NEXT
       ENDDO
       mid=p%chart%f%ent
       o=p%chart%f%a
    else
       P=>R%START
       DO I=1,(R%N-1)/2
          P=>P%NEXT
       ENDDO
       mid=p%chart%f%mid
       o=p%chart%f%o
    endif

    call get_length(r,ld)





  END SUBROUTINE locate_mid_frame

  SUBROUTINE LOCATE_FIBRE(R,PIN,I)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    TYPE(FIBRE), POINTER:: P,PIN
    INTEGER, INTENT(INOUT) :: I
    P=>R%START
    DO I=1,R%N
       IF(ASSOCIATED(PIN,P) ) EXIT
       P=>P%NEXT
    ENDDO
  END SUBROUTINE LOCATE_FIBRE


  SUBROUTINE GET_FREQ(R,FREQ)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    REAL(DP), INTENT(OUT) :: FREQ
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    P=>R%START
    FREQ=0.0_dp
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             IF(FREQ==0.0_dp) THEN
                FREQ=P%MAG%FREQ
             ELSEIF(FREQ>P%MAG%FREQ) THEN
                FREQ=P%MAG%FREQ
             ENDIF
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE GET_FREQ

  SUBROUTINE GET_loss(R,energy,deltap)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    REAL(DP), INTENT(OUT) :: energy,deltap
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    P=>R%START
    energy=0.0_dp
    deltap=0.0_dp
    DO I=1,R%N
       IF(P%MAG%kind==kind4) THEN
          energy=energy+p%mag%delta_e
          if(lielib_print(10)==1) write(6,*) p%mag%name
       ENDIF
       P=>P%NEXT
    ENDDO
    P=>R%START
    deltap=energy/p%mag%p%p0c
  END SUBROUTINE GET_loss

  SUBROUTINE GET_ALL(R,FREQ,VOLT,PHAS)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    REAL(DP), INTENT(OUT) :: FREQ,VOLT,PHAS
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    P=>R%START
    FREQ=0.0_dp;VOLT=0.0_dp;PHAS=0.0_dp;
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             FREQ=TWOPI*P%MAG%FREQ/CLIGHT
             VOLT=-P%MAG%VOLT*volt_c/P%MAG%P%P0C
             PHAS=P%MAG%PHAS
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE GET_ALL

  SUBROUTINE GET_ALL_mad_like(R,FREQ,VOLT,PHAS)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    REAL(DP), INTENT(OUT) :: FREQ,VOLT,PHAS
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    P=>R%START
    FREQ=0.0_dp;VOLT=0.0_dp;PHAS=0.0_dp;
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             FREQ=P%MAG%FREQ
             VOLT=P%MAG%VOLT
             PHAS=-P%MAG%PHAS
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE GET_ALL_mad_like

  SUBROUTINE locate_next_cav(R,di,P)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT) :: R
    integer, INTENT(INout) :: di
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    di=1
    if(associated(p)) P=>P%NEXT
    DO I=1,R%N
       if(associated(p)) then
          IF(ASSOCIATED(P%MAG%FREQ)) THEN
             IF(P%MAG%FREQ/=0.0_dp) THEN
                exit
             ENDIF
          ENDIF
          di=di+1
          P=>P%NEXT
       endif
    ENDDO
  END SUBROUTINE locate_next_cav

  SUBROUTINE locate_all_cav(R,pos)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT) :: R
    TYPE(FIBRE), POINTER:: P
    INTEGER, pointer ::  pos(:)
    integer i,ic
    ic=0
    P=>r%start
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             ic=ic+1
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
    allocate(pos(ic))
    pos=0
    ic=0
    P=>r%start
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             ic=ic+1
             pos(ic)=i
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO

  END SUBROUTINE locate_all_cav


  SUBROUTINE SET_FREQ(R,FREQ)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT) :: R
    REAL(DP), INTENT(IN) :: FREQ
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    P=>R%START
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             P%MAG%FREQ=FREQ
             P%MAGP%FREQ=FREQ
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE SET_FREQ

  SUBROUTINE ADD_FREQ(R,FREQ)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT) :: R
    REAL(DP), INTENT(IN) :: FREQ
    TYPE(FIBRE), POINTER:: P
    INTEGER I
    P=>R%START
    DO I=1,R%N
       IF(ASSOCIATED(P%MAG%FREQ)) THEN
          IF(P%MAG%FREQ/=0.0_dp) THEN
             P%MAG%FREQ=P%MAG%FREQ+FREQ
             P%MAGP%FREQ=P%MAGP%FREQ+FREQ
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE ADD_FREQ

  !END NEW
  SUBROUTINE ADDP_ANBN(EL,NM,F,V,ELECTRIC) ! EXTENDS THE ADD ROUTINES FROM THE ELEMENT(P) TO THE FIBRE
    IMPLICIT NONE
    TYPE(FIBRE),target, INTENT(INOUT) ::EL
    REAL(DP), INTENT(IN) ::V
    INTEGER, INTENT(IN) ::NM,F
    LOGICAL(LP), OPTIONAL :: ELECTRIC

    CALL ADD(EL%MAG,NM,F,V,ELECTRIC)
    CALL ADD(EL%MAGP,NM,F,V,ELECTRIC)

  END SUBROUTINE ADDP_ANBN


  SUBROUTINE PUT_APERTURE_FIB(EL,KIND,R,X,Y)
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: R(2),X,Y
    INTEGER,INTENT(IN):: KIND
    TYPE(FIBRE),target,INTENT(INOUT):: EL

    CALL PUT_APERTURE(EL%MAG,KIND,R,X,Y,0.0_dp,0.0_dp)
    CALL PUT_APERTURE(EL%MAGP,KIND,R,X,Y,0.0_dp,0.0_dp)

  END  SUBROUTINE PUT_APERTURE_FIB

  SUBROUTINE PUT_APERTURE_FIBt(EL,KIND,R,X,Y,dx,dy)
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: R(2),X,Y,dx,dy
    INTEGER,INTENT(IN):: KIND
    TYPE(FIBRE),target,INTENT(INOUT):: EL

    CALL PUT_APERTURE(EL%MAG,KIND,R,X,Y,dx,dy)
    CALL PUT_APERTURE(EL%MAGP,KIND,R,X,Y,dx,dy)

  END  SUBROUTINE PUT_APERTURE_FIBt

  SUBROUTINE REMOVE_APERTURE_FIB(EL)
    IMPLICIT NONE
    TYPE(FIBRE),target,INTENT(INOUT):: EL

    CALL REMOVE_APERTURE_EL(EL%MAG)
    CALL REMOVE_APERTURE_ELP(EL%MAGP)

  END  SUBROUTINE REMOVE_APERTURE_FIB

  SUBROUTINE  layout_WORK(r,S1) ! CHANGES THE ENERGY OF THE FIBRE AND TURNS THE ENERGY PATCH ON
    IMPLICIT NONE
    TYPE (WORK),INTENT(IN):: S1
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(fibre),pointer :: p
    integer i

    p=> r%start

    do i=1,r%n
       p=s1

       p=>p%next
    enddo

  END SUBROUTINE layout_WORK


  SUBROUTINE  FIBRE_WORK(S2,S1) ! CHANGES THE ENERGY OF THE FIBRE AND TURNS THE ENERGY PATCH ON
    IMPLICIT NONE
    TYPE (WORK),INTENT(IN):: S1
    TYPE(FIBRE),target,INTENT(INOUT):: S2
    TYPE (WORK) w
    if(s2%mag%vorname=="RESCALE".and.force_rescale) then
     w=s1
     w%power=1
     w%rescale=.true.
     S2%MAG=w
     S2%MAGP=w
    else
     S2%MAG=S1
     S2%MAGP=S1
    endif
    if(S1%power/=-1) then       ! just rescaling  -1=ramping
       S2%mass=S1%mass
       S2%BETA0=S1%BETA0
       S2%GAMMA0I=S1%GAMMA0I
       S2%GAMBET=S1%GAMBET
    endif


  END SUBROUTINE FIBRE_WORK

  SUBROUTINE  WORK_FIBRE(S2,S1)  ! SUCKS THE ENERGY OUT OF A FIBRE BY LOKING AT ELEMENT
    IMPLICIT NONE
    TYPE (FIBRE),target,INTENT(IN):: S1
    TYPE(WORK),INTENT(INOUT):: S2
    electron=my_true
    muon=s1%mass/pmae
    s2%mass=s1%mass
    S2=S1%MAG
    IF(ABS(S1%MAG%P%P0C-S1%MAGP%P%P0C)>1e-10_dp) THEN
       !w_p=0
       !w_p%NC=3
       !w_p%FC='(2(1X,A72,/),(1X,A72))'
       write(6,*) " BEWARE : ELEMENT AND ELEMENTP SEEM TO HAVE "
       write(6,*) " DIFFERENT REFERENCE ENERGIES!"
       write(6,'(1X,G21.14,1X,g21.14)')  S1%MAG%P%P0C,S1%MAGP%P%P0C
       ! call !write_e(100)
    ENDIF

  END SUBROUTINE WORK_FIBRE
  !  S-aperture

  SUBROUTINE  alloc_s_aperture(S1,APERTURE)  ! copy full fibre
    IMPLICIT NONE
    TYPE(MADX_APERTURE), OPTIONAL :: APERTURE
    TYPE (FIBRE),target,INTENT(INout):: S1

    if(associated(S1%mag%p%a)) call kill(S1%mag%p%a)
    if(associated(S1%magp%p%a)) call kill(S1%magp%p%a)

    call alloc(S1%mag%p%a,S1%mag%p%nst+1,APERTURE)
    call alloc(S1%magp%p%a,S1%magp%p%nst+1,APERTURE)

  END SUBROUTINE alloc_s_aperture

  SUBROUTINE  kill_s_aperture(S1)  ! copy full fibre
    IMPLICIT NONE
    TYPE (FIBRE),INTENT(INout):: S1

    if(associated(S1%mag%p%a)) call kill(S1%mag%p%a)
    if(associated(S1%magp%p%a)) call kill(S1%magp%p%a)


  END SUBROUTINE kill_s_aperture


  SUBROUTINE  copy_fibre(S1,S2)  ! copy full fibre
    IMPLICIT NONE
    TYPE (FIBRE),target,INTENT(IN):: S1
    TYPE(FIBRE),target,INTENT(INOUT):: S2

    call copy(s1%mag,s2%mag)
    call copy(s1%mag,s2%magp)

  END SUBROUTINE copy_fibre

  SUBROUTINE  MISALIGN_FIBRE_EQUAL(S2,S1) ! MISALIGNS FULL FIBRE; FILLS IN CHART AND MAGNET_CHART
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    TYPE(FIBRE),target,INTENT(INOUT):: S2

    CALL MISALIGN_FIBRE(S2,S1)

  END SUBROUTINE  MISALIGN_FIBRE_EQUAL

  SUBROUTINE  FIND_AFFINE_SIAMESE(S2,CN,FOUND) !
    ! FIND THE ELEMENT (CN) ON WHICH THE AFFINE_FRAME IS
    IMPLICIT NONE
    TYPE(FIBRE),TARGET,INTENT(INOUT):: S2
    !    TYPE(AFFINE_FRAME),POINTER :: AF
    TYPE(ELEMENT), POINTER :: C,CN
    INTEGER K
    LOGICAL(LP),INTENT(INOUT)::FOUND

    !    NULLIFY(AF)
    !    NULLIFY(CN)


    FOUND=MY_FALSE
    K=0

    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN


       C=>S2%MAG
       CN=>S2%MAG%SIAMESE

       IF(ASSOCIATED(C%SIAMESE_FRAME)) THEN
          !          AF=>C%SIAMESE_FRAME
          CN=>C
          FOUND=MY_TRUE
          !          WRITE(6,*) " HERE 1 ", FOUND
          RETURN
       ENDIF

       DO WHILE(.NOT.ASSOCIATED(C,CN))
          IF(ASSOCIATED(CN%SIAMESE_FRAME)) THEN
             !             AF=>CN%SIAMESE_FRAME
             FOUND=MY_TRUE
             !             WRITE(6,*) " HERE 2 ", FOUND
             EXIT
          ENDIF
          CN=>CN%SIAMESE
          k=k+1
          IF(K>10000)THEN
             WRITE(6,*) " TOO MANY IN SIAMESE "
             STOP 666
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE  FIND_AFFINE_SIAMESE

  SUBROUTINE FIND_FRAME_SIAMESE(MAG,B,EXI,ADD) !
    ! ACTUALLY CONSTRUCTS THE AFFINE FRAME FOUND PREVIOUSLY ON MAG
    ! FROM  mag%SIAMESE_FRAME%D,mag%SIAMESE_FRAME%ANGLE
    ! IF ADD=FALSE, START FROM ORIGINAL FIBRE POSITION
    ! COMPUTE B,EXI
    IMPLICIT NONE
    TYPE(ELEMENT), POINTER :: MAG
    REAL(DP), INTENT(INOUT) :: B(3),EXI(3,3)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD
    LOGICAL(LP) ADDIN

    ADDIN=.FALSE.

    IF(PRESENT(ADD)) ADDIN=ADD

    IF(ADDIN) THEN
       CALL INVERSE_FIND_PATCH(mag%p%F%a,mag%p%F%ENT, &
            mag%SIAMESE_FRAME%D,mag%SIAMESE_FRAME%ANGLE,B,EXI)
    ELSE
       CALL INVERSE_FIND_PATCH(mag%PARENT_FIBRE%CHART%F%a,mag%PARENT_FIBRE%CHART%F%ENT, &
            mag%SIAMESE_FRAME%D,mag%SIAMESE_FRAME%ANGLE,B,EXI)
    ENDIF

  END SUBROUTINE FIND_FRAME_SIAMESE

  RECURSIVE SUBROUTINE  MISALIGN_SIAMESE(S2,S1,OMEGA,BASIS,ADD,preserve_girder)
    ! SAME AS MISALIGN_FIBRE: DEFAULT IS THE O,MID OF S2
    !   UNLESS IT FINDS  TYPE(AFFINE_FRAME), POINTER :: SIAMESE_FRAME ON S2%SIAMESE CHAIN
    ! ON ONE SIAMESE IN THE CHAIN

    ! THIS IS OVERWRITEN IF OMEGA AND BASIS ARE PRESENT


    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    TYPE(FIBRE),TARGET,INTENT(INOUT):: S2
    TYPE(ELEMENT), POINTER :: C,CN
    TYPE(fibre), POINTER :: P
    integer k
    REAL(DP) OMEGAT(3),BASIST(3,3),B(3),EXI(3,3)
    real(dp) a1(3),e1(3,3),a2(3),e2(3,3),dg1(3),ag1(3),mis(6)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD,preserve_girder
    LOGICAL(LP) ADDIN,pres
    !    TYPE(AFFINE_FRAME),POINTER :: AF
    LOGICAL(LP) FOUND
    LOGICAL(LP) FOUNDg
    type(element), pointer :: caf
    !    REAL(DP) D(3),ANG(3)

    FOUND=.FALSE.
    ADDIN=.FALSE.
    pres=.FALSE.
    if(present(preserve_girder)) pres=preserve_girder
    CALL FIND_AFFINE_SIAMESE(S2,CN,FOUND)  ! Looking for siamese WITH FRAME
    IF(FOUND) CALL FIND_FRAME_SIAMESE(CN,B,EXI,ADD) ! FIND ACTUAL FRAME

    FOUNDg=.false.
    if(pres.and.associated(s2%mag%girders).and.(.not.addin)) then
       CALL FIND_AFFINE_GIRDER(S2,CAF,FOUNDg)
       if(foundg) then
          a1=caf%girder_frame%a
          e1=caf%girder_frame%ent
          a2=caf%girder_frame%b
          e2=caf%girder_frame%exi
          !         call INVERSE_FIND_PATCH(a1,e1,dg1,ag1,a2,e2)
          call FIND_PATCH(a1,e1,a2,e2,dg1,ag1)
          mis=0.0_dp
          call MISALIGN_siamese(S2,MIS)
          MIS(1:3)=DG1
          MIS(4:6)=AG1

          call MISALIGN_SIAMESE(S2,MIS,OMEGA=a1,BASIS=e1)
          call MISALIGN_SIAMESE(S2,S1,OMEGA,BASIS,ADD=my_true,preserve_girder=my_false)
          return
       endif
    endif

    IF(PRESENT(ADD)) ADDIN=ADD

    IF(PRESENT(OMEGA)) THEN    ! Arbitrary Origin
       OMEGAT=OMEGA
    ELSE
       OMEGAT=S2%CHART%F%O   ! Centre of magnet otherwise
    ENDIF

    IF(PRESENT(BASIS)) THEN      ! Arbitrary Basis
       BASIST=BASIS
    ELSE
       BASIST=S2%CHART%F%MID  ! Centre of Magnet Otherwise
    ENDIF

    IF((.NOT.PRESENT(OMEGA)).AND.(.NOT.PRESENT(BASIS))) THEN
       IF(FOUND) THEN   ! If no special basis and no special origin
          OMEGAT=B         ! and siamese is found, then it uses the siamese basis
          BASIST=EXI        ! Notice that if ADD=true, the siamese frames move with the magnets
       ENDIF
    ENDIF

    CALL MISALIGN_FIBRE(S2,S1,OMEGAT,BASIST,ADD=ADDIN)
    k=1

    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN
       C=>S2%MAG
       CN=>S2%MAG%SIAMESE
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          P=>CN%PARENT_FIBRE
          CALL MISALIGN_FIBRE(P,S1,OMEGAT,BASIST,ADD=ADDIN)
          CN=>CN%SIAMESE
          k=k+1
       ENDDO
    ENDIF
    !    CALL MOVE_SIAMESE_FRAME(S2%MAG)
    if(global_verbose) write(6,*) k, " magnet misaligned "
  END SUBROUTINE  MISALIGN_SIAMESE

  SUBROUTINE  FIND_AFFINE_GIRDER(S2,CN,FOUND) !
    ! SAME AS FIND_AFFINE_SIAMESE
    ! LOCATES MAGNET CN WHERE C%GIRDER_FRAME IS
    IMPLICIT NONE
    TYPE(FIBRE),TARGET,INTENT(INOUT):: S2
    !    TYPE(AFFINE_FRAME),POINTER :: AF
    TYPE(ELEMENT), POINTER :: C,CN
    INTEGER K
    LOGICAL(LP),INTENT(INOUT)::FOUND

    !    NULLIFY(AF)
    !    NULLIFY(CN)


    FOUND=MY_FALSE
    K=0

    IF(ASSOCIATED(S2%MAG%GIRDERS)) THEN
       C=>S2%MAG
       CN=>S2%MAG%GIRDERS

       IF(ASSOCIATED(C%GIRDER_FRAME)) THEN
          !          AF=>C%GIRDER_FRAME
          CN=>C
          FOUND=MY_TRUE
          RETURN
       ENDIF

       DO WHILE(.NOT.ASSOCIATED(C,CN))
          IF(ASSOCIATED(CN%GIRDER_FRAME)) THEN
             !             AF=>CN%GIRDER_FRAME
             FOUND=MY_TRUE
             EXIT
          ENDIF
          CN=>CN%GIRDERS
          k=k+1
          IF(K>10000)THEN
             WRITE(6,*) " TOO MANY IN GIRDER "
             STOP 666
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE  FIND_AFFINE_GIRDER

  SUBROUTINE FIND_FRAME_GIRDER(MAG,B,EXI,ADD) !
    ! ACTUALLY LOCATES THE FRAME
    IMPLICIT NONE
    TYPE(ELEMENT), POINTER :: MAG
    REAL(DP), INTENT(INOUT) :: B(3),EXI(3,3)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD
    LOGICAL(LP) ADDIN

    ADDIN=.FALSE.

    IF(PRESENT(ADD)) ADDIN=ADD

    IF(.NOT.ADDIN) THEN
       mag%GIRDER_FRAME%EXI=mag%GIRDER_FRAME%ENT   ! ORIGINAL FRAME OF GIRDER
       mag%GIRDER_FRAME%B=mag%GIRDER_FRAME%A       ! ORIGINAL POSITION OF GIRDER
    ENDIF
    EXI=mag%GIRDER_FRAME%EXI
    B=mag%GIRDER_FRAME%B

  END SUBROUTINE FIND_FRAME_GIRDER

  SUBROUTINE  EXTRACT_GIRDER_FRAME(S2,A,ENT,FOUND) !
    ! USED IN PTC_GIRDER
    IMPLICIT NONE
    TYPE(element),TARGET,INTENT(INOUT):: S2
    TYPE(AFFINE_FRAME),POINTER :: AF
    TYPE(ELEMENT), POINTER :: C,CN
    INTEGER K
    LOGICAL(LP),INTENT(INOUT)::FOUND
    REAL(DP),INTENT(INOUT):: ENT(3,3),A(3)

    NULLIFY(AF)
    NULLIFY(CN)


    FOUND=MY_FALSE
    K=0

    IF(ASSOCIATED(S2%GIRDERS)) THEN
       C=>S2
       CN=>S2%GIRDERS

       IF(ASSOCIATED(C%GIRDER_FRAME)) THEN
          AF=>C%GIRDER_FRAME
          CN=>C
          FOUND=MY_TRUE
          RETURN
       ENDIF

       DO WHILE(.NOT.ASSOCIATED(C,CN))
          IF(ASSOCIATED(CN%GIRDER_FRAME)) THEN
             AF=>CN%GIRDER_FRAME
             FOUND=MY_TRUE
             EXIT
          ENDIF
          CN=>CN%GIRDERS
          k=k+1
          IF(K>10000)THEN
             WRITE(6,*) " TOO MANY IN GIRDER "
             STOP 666
          ENDIF
       ENDDO
    ENDIF

    IF(FOUND) THEN
       ENT=AF%ENT
       A=AF%A
    ENDIF
  END SUBROUTINE  EXTRACT_GIRDER_FRAME



  SUBROUTINE  MISALIGN_GIRDER(S2,S1,OMEGA,BASIS,ADD) !
    ! SIMILAR TO MISALIGN_SIAMESE
    ! COMMENT DIFFERENCES ONLY
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    TYPE(FIBRE),TARGET,INTENT(INOUT):: S2
    TYPE(ELEMENT), POINTER :: C,CN,CAF
    TYPE(fibre), POINTER :: P
    integer k
    REAL(DP) OMEGAT(3),BASIST(3,3),B(3),EXI(3,3),T_GLOBAL(3)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD
    LOGICAL(LP) ADDIN
    LOGICAL(LP) FOUND
    !    REAL(DP) D(3),ANG(3)
    TYPE(MAGNET_FRAME), POINTER :: F

    FOUND=.FALSE.
    ADDIN=.FALSE.
    CALL FIND_AFFINE_GIRDER(S2,CAF,FOUND)
    IF(FOUND) CALL FIND_FRAME_GIRDER(CAF,B,EXI,ADD)

    IF(PRESENT(ADD)) ADDIN=ADD

    IF(PRESENT(OMEGA)) THEN
       OMEGAT=OMEGA
    ELSE
       OMEGAT=S2%CHART%F%O
    ENDIF

    IF(PRESENT(BASIS)) THEN
       BASIST=BASIS
    ELSE
       BASIST=S2%CHART%F%MID
    ENDIF

    IF((.NOT.PRESENT(OMEGA)).AND.(.NOT.PRESENT(BASIS))) THEN
       IF(FOUND) THEN
          OMEGAT=B
          BASIST=EXI
       ENDIF
    ENDIF


    CALL MISALIGN_FIBRE(S2,S1,OMEGAT,BASIST,ADD=ADDIN)
    k=1

    IF(ASSOCIATED(S2%MAG%GIRDERS)) THEN
       C=>S2%MAG
       CN=>S2%MAG%GIRDERS
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          P=>CN%PARENT_FIBRE
          CALL MISALIGN_FIBRE(P,S1,OMEGAT,BASIST,ADD=ADDIN)
          CN=>CN%GIRDERS
          k=k+1
       ENDDO
    ENDIF



    IF(FOUND) THEN   !!! THE ORIGINAL GIRDER FRAME IS STILL GIRDER_FRAME%ENT AND GIRDER_FRAME%A
       !                    FINAL FRAME AFTER MISALIGNMENTS MUST BE COMPUTED
       call alloc(f)
       f%a=b
       f%ent=exi
       CALL ROTATE_FRAME(F,OMEGAT,S1(4:6),1,BASIS=BASIST)
       CALL   GEO_ROT(BASIST,S1(4:6),1)
       CALL CHANGE_BASIS(S1(1:3),BASIST,T_GLOBAL,GLOBAL_FRAME)
       F%A=F%A+T_GLOBAL
       CAF%GIRDER_FRAME%EXI=F%ent
       CAF%GIRDER_FRAME%B=F%A
       call kill(f)
    ENDIF

    if(global_verbose)     write(6,*) k, " magnet misaligned "
  END SUBROUTINE  MISALIGN_GIRDER



  recursive SUBROUTINE  MISALIGN_FIBRE(S2,S1,OMEGA,BASIS,ADD,preserve_girder)
    ! MISALIGNS FULL FIBRE; FILLS IN CHART AND MAGNET_CHART
    ! changed  add=true add extra misalignments TO EXISTING ONES
    ! O AND MID BY DEFAUTL, OTHERWISE OMEGA AND BASIS
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD,preserve_girder
    TYPE(FIBRE),target,INTENT(INOUT):: S2
    REAL(DP) ANGLE(3),T_GLOBAL(3),d(3),r(3)
    TYPE(MAGNET_FRAME), POINTER :: F,F0
    REAL(DP) D_IN(3),D_OUT(3),OMEGAT(3),BASIST(3,3)
    INTEGER I
    LOGICAL(LP) ADDIN,pres
    LOGICAL(LP) FOUNDg
    type(element), pointer :: caf
    real(dp) a1(3),e1(3,3),a2(3),e2(3,3),dg1(3),ag1(3),mis(6)



    ADDIN=.FALSE.
    pres=.FALSE.
    IF(PRESENT(ADD)) ADDIN=ADD
    if(present(preserve_girder)) pres=preserve_girder

    FOUNDg=.false.
    if(pres.and.associated(s2%mag%girderS).and.(.not.addin)) then
       CALL FIND_AFFINE_GIRDER(S2,CAF,FOUNDg)
       if(foundg) then
          a1=caf%girder_frame%a
          e1=caf%girder_frame%ent
          a2=caf%girder_frame%b
          e2=caf%girder_frame%exi
          !         call INVERSE_FIND_PATCH(a1,e1,dg1,ag1,a2,e2)
          call FIND_PATCH(a1,e1,a2,e2,dg1,ag1)
          mis=0.0_dp
          call MISALIGN_fibre(S2,MIS)
          MIS(1:3)=DG1
          MIS(4:6)=AG1

          call MISALIGN_fibre(S2,MIS,OMEGA=a1,BASIS=e1)
          call MISALIGN_fibre(S2,S1,OMEGA,BASIS,ADD=my_true,preserve_girder=my_false)
          return
       endif
    endif



    IF(ASSOCIATED(S2%CHART)) THEN
       !       IF(.NOT.ASSOCIATED(S2%MAG%D) ) ALLOCATE(S2%MAG%D(3))
       !       IF(.NOT.ASSOCIATED(S2%MAG%R) ) ALLOCATE(S2%MAG%R(3))
       !       IF(.NOT.ASSOCIATED(S2%MAGP%D)) ALLOCATE(S2%MAGP%D(3))
       !       IF(.NOT.ASSOCIATED(S2%MAGP%R)) ALLOCATE(S2%MAGP%R(3))
       !       DO I=1,3
       !          S2%MAG%D(I)=S1(I);   S2%MAGP%D(I)=S1(I);
       !          S2%MAG%R(I)=S1(3+I); S2%MAGP%R(I)=S1(3+I);
       !       ENDDO
       DO I=1,3
          D(I)=S1(I);   !D(I)=S1(I);
          R(I)=S1(3+I); !R(I)=S1(3+I);
       ENDDO
       S2%CHART%D_IN=0.0_dp;S2%CHART%D_OUT=0.0_dp;
       S2%CHART%ANG_IN=0.0_dp;S2%CHART%ANG_OUT=0.0_dp;
       S2%MAG%MIS=.TRUE.
       S2%MAGP%MIS=.TRUE.

       ! ADD CODE HERE
       CALL ALLOC(F)
       CALL ALLOC(F0)
       ! MOVE THE ORIGINAL INTERNAL CHART F
       IF(ADDIN) THEN
          F=S2%mag%p%f
          CALL SURVEY_NO_PATCH(S2,MAGNETFRAME=F0)
       ELSE
          CALL SURVEY_NO_PATCH(S2,MAGNETFRAME=F0)
          F=F0
       ENDIF

       ANGLE=r  !S2%MAG%R
       IF(PRESENT(BASIS)) THEN
          BASIST=BASIS
       ELSE
          BASIST=F%MID
       ENDIF
       IF(PRESENT(OMEGA)) THEN
          OMEGAT=OMEGA
       ELSE
          OMEGAT=F%O
       ENDIF

       CALL ROTATE_FRAME(F,OMEGAT,ANGLE,1,BASIS=BASIST)

       IF(PRESENT(BASIS)) THEN   ! MUST ROTATE THAT FRAME AS WELL FOR CONSISTENCY IN DEFINITION WHAT A MISALIGNMENT IS IN PTC
!          CALL   GEO_ROT(BASIST,ANGLE,1)
           e1=basist
          CALL   GEO_ROT(BASIST,e2,ANGLE,e1)   ! 2018 correction: agrees with successive rotations followed by rotations
           basist=e2
       ELSE
          BASIST=F%MID    ! ALREADY ROTATED
       ENDIF

       !       CALL CHANGE_BASIS(S2%MAG%D,BASIST,T_GLOBAL,GLOBAL_FRAME)
       CALL CHANGE_BASIS(D,BASIST,T_GLOBAL,GLOBAL_FRAME)


       F%A=F%A+T_GLOBAL
       F%O=F%O+T_GLOBAL
       F%B=F%B+T_GLOBAL



       CALL COMPUTE_ENTRANCE_ANGLE(F0%ENT,F%ENT,S2%CHART%ANG_IN)
       CALL COMPUTE_ENTRANCE_ANGLE(F%EXI,F0%EXI,S2%CHART%ANG_OUT)

       D_IN=F%A-F0%A
       D_OUT=F0%B-F%B

       !        WRITE(6,*) " IN GLOBAL BASIS D_IN AND D_OUT"

       !        WRITE(6,*) D_IN
       !        WRITE(6,*) D_OUT

       CALL CHANGE_BASIS(D_IN,GLOBAL_FRAME,S2%CHART%D_IN,F%ENT)
       CALL CHANGE_BASIS(D_OUT,GLOBAL_FRAME,S2%CHART%D_OUT,F0%EXI)

       !        WRITE(6,*) " IN LOCAL  BASIS D_IN AND D_OUT AS WELL AS ANGLES"
       !       WRITE(6,*) " ***************************************"
       !      WRITE(6,*) S2%CHART%ANG_IN
       !      WRITE(6,*) S2%CHART%ANG_OUT
       !      WRITE(6,*) S2%CHART%D_IN
       !      WRITE(6,*) S2%CHART%D_OUT

       !      WRITE(6,*) " ***************************************"

       CALL KILL(F)
       CALL KILL(F0)
       IF(ASSOCIATED(F)) deallocate(f)
       IF(ASSOCIATED(F0)) deallocate(f0)

       CALL SURVEY_NO_PATCH(S2)

       IF(ASSOCIATED(S2%T1)) THEN
          IF(ASSOCIATED(S2%T1%A)) THEN
             CALL fill_survey_ONE_FIBRE(S2)
          ENDIF
       ENDIF

    ELSE
       !w_p=0
       !w_p%NC=1
       !w_p%FC='((1X,A72))'
       WRITE(6,'(1X,A39,1X,A16)') " CANNOT MISALIGN THIS FIBRE: NO CHARTS ", S2%MAG%NAME
       ! call !write_e(100)
    ENDIF


  END SUBROUTINE MISALIGN_FIBRE


  SUBROUTINE  MAD_MISALIGN_FIBRE(S2,S1) ! MISALIGNS FULL FIBRE; FILLS IN CHART AND MAGNET_CHART
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    TYPE(FIBRE),target,INTENT(INOUT):: S2
    REAL(DP) ENT(3,3),ENT1(3,3),ENT2(3,3),T(3),MAD_ANGLE(3),T_GLOBAL(3),ANGLE(3),MIS(6)
    real(dp) r1(3,3),r2(3,3),r3(3,3)
    ent=S2%CHART%F%ent
    T(1)=S1(1);T(2)=S1(2);T(3)=S1(3);
    MAD_ANGLE(1)=-S1(4)
    MAD_ANGLE(2)=-S1(5)
    MAD_ANGLE(3)=S1(6)



    ANGLE=0.0_dp; ANGLE(1)=MAD_ANGLE(1)
    ent1=ent
    ent2=ent
    CALL GEO_ROT(ENT1,ENT,ANGLE,ENT2)
 
    ANGLE=0.0_dp; ANGLE(2)=MAD_ANGLE(2)
    ent1=ent
    ent2=ent
    CALL GEO_ROT(ENT1,ENT,ANGLE,ENT2)

    ANGLE=0.0_dp; ANGLE(3)=MAD_ANGLE(3)
    ent1=ent
    ent2=ent
    CALL GEO_ROT(ENT1,ENT,ANGLE,ENT2)
    




    CALL COMPUTE_ENTRANCE_ANGLE(S2%CHART%F%ent,ENT,ANGLE)

 
     ent1=S2%CHART%F%ent
     ent2=ent
     CALL CHANGE_BASIS(T,ENT1,T_GLOBAL,ent2)

 

    MIS(1:3)=T_GLOBAL
    MIS(4:6)=ANGLE

    ENT=S2%CHART%F%ent
    T=S2%CHART%F%A

        call MISALIGN_FIBRE(S2,MIS,T,ENT)

  END SUBROUTINE MAD_MISALIGN_FIBRE




  ! NEW ROUTINES TO CHANGE LAYOUT using only magnets!!!!

  SUBROUTINE TRANSLATE_girder(S2,D,ORDER,BASIS,PATCH,PREC) ! TRANSLATES A fibre
    IMPLICIT NONE
    TYPE (fibre),TARGET,INTENT(INOUT):: S2
    TYPE (element),pointer :: c,cn,caf
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    LOGICAL(LP),OPTIONAL :: PATCH
    REAL(DP),OPTIONAL :: PREC
    LOGICAL(LP) FOUND
    REAL(DP) exi(3,3), BASISt(3,3),b(3),t_global(3)
    integer k

    FOUND=.FALSE.
    CALL FIND_AFFINE_girder(S2,CAF,FOUND)
    IF(FOUND) CALL FIND_FRAME_girder(CAF,B,EXI,ADD=my_false)


    IF(PRESENT(BASIS)) THEN
       BASIST=BASIS
    ELSE
       if(found) then
          BASIST=exi
       else
          BASIST=global_frame
       endif
    ENDIF

    c=>s2%mag
    !    CALL TRANSLATE_magnet(c,D,ORDER,BASIST,PATCH=.false.)
    CALL TRANSLATE_magnet(c,D,ORDER,BASIST,PATCH,PREC)
    k=1

    IF(ASSOCIATED(S2%MAG%girderS)) THEN
       C=>S2%MAG
       CN=>S2%MAG%girderS
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          !       CALL TRANSLATE_magnet(cn,D,ORDER,BASIST,PATCH=.false.)
          CALL TRANSLATE_magnet(cn,D,ORDER,BASIST,PATCH,PREC)
          CN=>CN%girderS
          k=k+1
       ENDDO
    ENDIF

    if(global_verbose)     write(6,*) k, " magnets translated in girder "

    !    c=>s2%mag
    !         call patch_magnet(c,PATCH,PREC)
    !    k=1
    !
    !    IF(ASSOCIATED(S2%MAG%girder)) THEN
    !     C=>S2%MAG
    !     CN=>S2%MAG%girder
    !     DO WHILE(.NOT.ASSOCIATED(C,CN))
    !         call patch_magnet(cn,PATCH,PREC)
    !!      CN=>CN%girder
    !     k=k+1
    !    ENDDO
    !   ENDIF
    !if(global_verbose)     write(6,*) k, " magnets patched in translated in girder "


    IF(FOUND) THEN
       CALL CHANGE_BASIS(D,BASIST,T_GLOBAL,GLOBAL_FRAME)
       CAF%GIRDER_FRAME%a=CAF%GIRDER_FRAME%a+T_GLOBAL
       CAF%GIRDER_FRAME%B=CAF%GIRDER_FRAME%b+T_GLOBAL
    ENDIF
    !  ETIENNE DOES NOT UNDERSTAND THE CODE BELOW

    !     CALL FIND_AFFINE_girder(S2,CN,FOUND)

    !     IF(FOUND) CALL FIND_FRAME_girder(CN,B,EXI,ADD=my_false)

  END SUBROUTINE TRANSLATE_girder


  SUBROUTINE TRANSLATE_siamese(S2,D,ORDER,BASIS,PATCH,PREC) ! TRANSLATES A SIAMESE
    IMPLICIT NONE
    TYPE (fibre),TARGET,INTENT(INOUT):: S2
    TYPE (element),pointer :: c,cn
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    LOGICAL(LP),OPTIONAL :: PATCH
    REAL(DP),OPTIONAL :: PREC
    LOGICAL(LP) FOUND
    REAL(DP)  exi(3,3), BASISt(3,3),b(3)
    integer k

    FOUND=.FALSE.
    CALL FIND_AFFINE_SIAMESE(S2,CN,FOUND)
    IF(FOUND) CALL FIND_FRAME_SIAMESE(CN,B,EXI,ADD=my_false)


    IF(PRESENT(BASIS)) THEN
       BASIST=BASIS
    ELSE
       if(found) then
          BASIST=exi
       else
          BASIST=global_frame
       endif
    ENDIF

    c=>s2%mag
    CALL TRANSLATE_magnet(c,D,ORDER,BASIST,PATCH,PREC)
    !    CALL TRANSLATE_magnet(c,D,ORDER,BASIST,PATCH=.false.)
    k=1

    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN
       C=>S2%MAG
       CN=>S2%MAG%SIAMESE
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          CALL TRANSLATE_magnet(cn,D,ORDER,BASIST,PATCH,PREC)
          !       CALL TRANSLATE_magnet(cn,D,ORDER,BASIST,PATCH=.false.)
          CN=>CN%SIAMESE
          k=k+1
       ENDDO
    ENDIF

    !if(global_verbose)  write(6,*) k, " magnet translated "

    !    c=>s2%mag
    !         call patch_magnet(c,PATCH,PREC)
    !    k=1
    !
    !    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN
    !     C=>S2%MAG
    !     CN=>S2%MAG%SIAMESE
    !     DO WHILE(.NOT.ASSOCIATED(C,CN))
    !         call patch_magnet(cn,PATCH,PREC)
    !      CN=>CN%SIAMESE
    !      k=k+1
    !     ENDDO
    !    ENDIF
    !  if(global_verbose)   write(6,*) k, " magnets patched in translated in siamese "

  END SUBROUTINE TRANSLATE_siamese


  SUBROUTINE TRANSLATE_magnet(R,D,ORDER,BASIS,PATCH,PREC) ! TRANSLATES A fibre
    IMPLICIT NONE
    TYPE (element),TARGET,INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    TYPE(FIBRE), POINTER::P
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    LOGICAL(LP),OPTIONAL :: PATCH
    REAL(DP),OPTIONAL :: PREC
    LOGICAL(LP) PAT
    REAL(DP) PREC0
    type(fibre_appearance), pointer :: dk
    integer k

    PREC0=PUNY
    PAT=MY_FALSE

    IF(PRESENT(PATCH)) PAT=PATCH
    IF(PRESENT(PREC)) PREC0=PREC

    p=> r%parent_fibre

    call TRANSLATE(p,D,ORDER,BASIS)

    if(pat) then
       k=0
       if(associated(R%doko)) then  !!! PATCH TO DOKO'S  IF CREATED USING DNA I.E. APPEND_POINT
          dk=>r%doko
          do while(associated(dk))
             p=> dk%parent_fibre
             call FIND_PATCH(p,p%next,NEXT=my_false,ENERGY_PATCH=my_true,prec=PREC0)
             call FIND_PATCH(p%previous,p,NEXT=my_true,ENERGY_PATCH=my_true,prec=PREC0)
             k=k+1
             dk=>dk%next
          enddo
          if(global_verbose)    write(6,*) "in translate_magnet patched ",k,"times using doko"
       else    !!!   FOR COMPATIBILITY MODE   FIBRE=MAGNET
          call FIND_PATCH(p,p%next,NEXT=my_false,ENERGY_PATCH=my_true,prec=PREC0)
          call FIND_PATCH(p%previous,p,NEXT=my_true,ENERGY_PATCH=my_true,prec=PREC0)
       endif
    endif

  END SUBROUTINE TRANSLATE_magnet

  SUBROUTINE rotate_magnet(R,Ang,OMEGA,ORDER,BASIS,PATCH,PREC) ! TRANSLATES A fibre
    IMPLICIT NONE
    TYPE (element),TARGET,INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: ang(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    TYPE(FIBRE), POINTER::P
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    LOGICAL(LP),OPTIONAL :: PATCH
    REAL(DP),OPTIONAL :: PREC
    LOGICAL(LP) PAT
    REAL(DP) PREC0,omega(3)
    type(fibre_appearance), pointer :: dk
    integer k

    PREC0=PUNY
    PAT=MY_FALSE

    IF(PRESENT(PATCH)) PAT=PATCH
    IF(PRESENT(PREC)) PREC0=PREC

    p=> r%parent_fibre

    call rotate(p,OMEGA,Ang,ORDER,BASIS)

    if(pat) then
       k=0
       if(associated(R%doko)) then
          dk=>r% doko
          do while(associated(dk))  !!! PATCH TO DOKO'S  IF CREATED USING DNA I.E. APPEND_POINT
             p=> dk%parent_fibre
             call FIND_PATCH(p,p%next,NEXT=my_false,ENERGY_PATCH=my_true,prec=PREC0)
             call FIND_PATCH(p%previous,p,NEXT=my_true,ENERGY_PATCH=my_true,prec=PREC0)
             k=k+1
             dk=>dk%next
          enddo
          if(global_verbose)     write(6,*) "in rotate_magnet patched ",k,"times using doko"
       else     !!!   FOR COMPATIBILITY MODE   FIBRE=MAGNET
          call FIND_PATCH(p,p%next,NEXT=my_false,ENERGY_PATCH=my_true,prec=PREC0)
          call FIND_PATCH(p%previous,p,NEXT=my_true,ENERGY_PATCH=my_true,prec=PREC0)
       endif
    endif

  END SUBROUTINE rotate_magnet

  !  SUBROUTINE patch_magnet(R,PATCH,PREC) ! TRANSLATES A fibre
  !    IMPLICIT NONE
  !    TYPE (element),TARGET,INTENT(INOUT):: R
  !    TYPE(FIBRE), POINTER::P
  !    LOGICAL(LP),OPTIONAL :: PATCH
  !    REAL(DP),OPTIONAL :: PREC
  !    LOGICAL(LP) PAT
  !    REAL(DP) PREC0
  !    type(fibre_appearance), pointer :: dk
  !    integer k
  !
  !    PREC0=PUNY
  !    PAT=MY_FALSE
  !
  !    IF(PRESENT(PATCH)) PAT=PATCH
  !    IF(PRESENT(PREC)) PREC0=PREC
  !
  !    p=> r%parent_fibre
  !
  !
  !    if(pat) then
  !     k=0
  !     if(associated(R%doko)) then
  !      dk=>r%doko
  !       do while(associated(dk))
  !         p=> dk%parent_fibre
  !         call FIND_PATCH(p,p%next,NEXT=my_false,ENERGY_PATCH=my_true,prec=PREC0)
  !         call FIND_PATCH(p%previous,p,NEXT=my_true,ENERGY_PATCH=my_true,prec=PREC0)
  !         k=k+1
  !         dk=>dk%next
  !
  !       enddo
  ! if(global_verbose)      write(6,*) " patched ",k,"times using doko"
  !     else
  !      call FIND_PATCH(p,p%next,NEXT=my_false,ENERGY_PATCH=my_true,prec=PREC0)
  !      call FIND_PATCH(p%previous,p,NEXT=my_true,ENERGY_PATCH=my_true,prec=PREC0)
  !     endif
  !    endif
  !
  !  END SUBROUTINE patch_magnet

  SUBROUTINE rotate_siamese(S2,Ang,OMEGA,ORDER,BASIS,PATCH,PREC) !
    IMPLICIT NONE
    TYPE (fibre),TARGET,INTENT(INOUT):: s2
    REAL(DP),INTENT(IN):: ang(3)
    TYPE (element),pointer :: c,cn
    REAL(DP), OPTIONAL :: BASIS(3,3),omega(3)
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    LOGICAL(LP),OPTIONAL :: PATCH
    REAL(DP),OPTIONAL :: PREC
    LOGICAL(LP) FOUND
    REAL(DP) b(3),exi(3,3), BASISt(3,3),omegat(3)
    integer k

    FOUND=.FALSE.
    CALL FIND_AFFINE_SIAMESE(S2,CN,FOUND)
    IF(FOUND) CALL FIND_FRAME_SIAMESE(CN,B,EXI,ADD=my_false)

    !     write(6,*)found, b

    IF(PRESENT(BASIS)) THEN
       BASIST=BASIS
    ELSE
       if(found) then
          BASIST=exi
       else
          BASIST=global_frame
       endif
    ENDIF

    IF(PRESENT(OMEGA)) THEN
       OMEGAT=OMEGA
    ELSE
       if(found) then
          OMEGAT=b
       else
          OMEGAT=global_origin
       endif
    ENDIF



    c=>s2%mag
    CALL rotate_magnet(c,Ang,OMEGAt,ORDER,BASISt,PATCH,PREC)
    !    CALL rotate_magnet(c,Ang,OMEGAt,ORDER,BASISt,PATCH=.false.)
    k=1

    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN
       C=>S2%MAG
       CN=>S2%MAG%SIAMESE
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          !        CALL rotate_magnet(cn,Ang,OMEGAt,ORDER,BASISt,PATCH=.false.)
          CALL rotate_magnet(cn,Ang,OMEGAt,ORDER,BASISt,PATCH,PREC)
          CN=>CN%SIAMESE
          k=k+1
       ENDDO
    ENDIF

    !  if(global_verbose)   write(6,*) k, " magnets rotated in siamese"
    !    c=>s2%mag
    !         call patch_magnet(c,PATCH,PREC)
    !    k=1
    !
    !    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN
    !     C=>S2%MAG
    !     CN=>S2%MAG%SIAMESE
    !     DO WHILE(.NOT.ASSOCIATED(C,CN))
    !         call patch_magnet(cn,PATCH,PREC)
    !      CN=>CN%SIAMESE
    !      k=k+1
    !     ENDDO
    !    ENDIF
    !  if(global_verbose)   write(6,*) k, " magnets patched in rotated in siamese "

  END SUBROUTINE rotate_siamese

  SUBROUTINE rotate_girder(S2,Ang,OMEGA,ORDER,BASIS,PATCH,PREC) !
    IMPLICIT NONE
    TYPE (fibre),TARGET,INTENT(INOUT):: s2
    REAL(DP),INTENT(IN):: ang(3)
    TYPE (element),pointer :: c,cn,caf
    REAL(DP), OPTIONAL :: BASIS(3,3),omega(3)
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    LOGICAL(LP),OPTIONAL :: PATCH
    REAL(DP),OPTIONAL :: PREC
    LOGICAL(LP) FOUND
    REAL(DP) b(3),exi(3,3), BASISt(3,3),omegat(3)
    TYPE(MAGNET_FRAME), POINTER :: F
    integer k

    FOUND=.FALSE.
    CALL FIND_AFFINE_girder(S2,CAF,FOUND)
    IF(FOUND) CALL FIND_FRAME_girder(CAF,B,EXI,ADD=my_false)


    IF(PRESENT(BASIS)) THEN
       BASIST=BASIS
    ELSE
       if(found) then
          BASIST=exi
       else
          BASIST=global_frame
       endif
    ENDIF

    IF(PRESENT(OMEGA)) THEN
       OMEGAT=OMEGA
    ELSE
       if(found) then
          OMEGAT=b
       else
          OMEGAT=global_origin
       endif
    ENDIF



    c=>s2%mag
    !    CALL rotate_magnet(c,Ang,OMEGAt,ORDER,BASISt,PATCH=.false.)
    CALL rotate_magnet(c,Ang,OMEGAt,ORDER,BASISt,PATCH,PREC)
    k=1

    IF(ASSOCIATED(S2%MAG%girderS)) THEN
       C=>S2%MAG
       CN=>S2%MAG%girderS
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          !        CALL rotate_magnet(cn,Ang,OMEGAt,ORDER,BASISt,PATCH=.false.)
          CALL rotate_magnet(cn,Ang,OMEGAt,ORDER,BASISt,PATCH,PREC)
          CN=>CN%girderS
          k=k+1
       ENDDO
    ENDIF
    if(global_verbose)   write(6,*) k, " magnets rotated in girder "


    !    c=>s2%mag
    !         call patch_magnet(c,PATCH,PREC)
    !    k=1

    !    IF(ASSOCIATED(S2%MAG%girder)) THEN
    !     C=>S2%MAG
    !     CN=>S2%MAG%girder
    !     DO WHILE(.NOT.ASSOCIATED(C,CN))
    !         call patch_magnet(cn,PATCH,PREC)
    !      CN=>CN%girder
    !      k=k+1
    !     ENDDO
    !    ENDIF
    !if(global_verbose)     write(6,*) k, " magnets patched in rotated  girder "

!!! SAME GYMNASTICS AS IN TRANSLATE GIRDER
    IF(FOUND) THEN
       call alloc(f)
       f%a=CAF%GIRDER_FRAME%a
       f%ent=CAF%GIRDER_FRAME%ent
       f%b=CAF%GIRDER_FRAME%B
       f%exi=CAF%GIRDER_FRAME%exi
       CALL ROTATE_FRAME(F,OMEGAT,ang,ORDER,BASIS=BASIST)
       CAF%GIRDER_FRAME%ENT=F%ENT
       CAF%GIRDER_FRAME%A=F%A
       CAF%GIRDER_FRAME%EXI=F%EXI
       CAF%GIRDER_FRAME%B=F%B
       call kill(f)
    ENDIF




  END SUBROUTINE rotate_girder



  ! NEW ROUTINES TO CHANGE LAYOUT

  SUBROUTINE  TRANSLATE_layout(R,D,I1,I2,ORDER,BASIS) ! TRANSLATES A LAYOUT
    IMPLICIT NONE
    TYPE (LAYOUT),INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    TYPE(FIBRE), POINTER::P
    INTEGER I,I11,I22
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER,I1,I2
    ! THIS ROUTINE TRANSLATE THE ENTIRE LINE BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE D

    P=>R%START

    I11=1
    I22=R%N
    IF(PRESENT(I1))  I11=I1
    IF(PRESENT(I2))  I22=I2
    DO I=1,I11-1
       P=>P%NEXT
    ENDDO


    DO I=1,I22-I11+1
       CALL TRANSLATE_Fibre(P,D,ORDER,BASIS,dogirder=my_true)   ! DOGIRDER IS SET TO TRUE
       P=>P%NEXT
    ENDDO


  END SUBROUTINE TRANSLATE_layout

  SUBROUTINE TRANSLATE_Fibre(R,D,ORDER,BASIS,dogirder) ! TRANSLATES A fibre
    ! THIS ROUTINE TRANSLATE THE ENTIRE LINE BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE D
    IMPLICIT NONE
    TYPE (FIBRE),TARGET,INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    TYPE(FIBRE), POINTER::P,pp
    INTEGER IORDER
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    TYPE(INTEGRATION_NODE), POINTER :: T
    TYPE(element), POINTER::caf
    logical(lp), OPTIONAL, INTENT(IN) :: dogirder
    REAL(DP) DD(3)
    type(fibre_appearance), pointer :: dk

    logical(lp) dog
    dog=.false.
    if(present(dogirder)) dog=dogirder

    P=>R
    Pp=>R

    IORDER=1
    DD=D
    IF(PRESENT(ORDER)) IORDER=ORDER
    IF(PRESENT(BASIS)) THEN
       CALL CHANGE_BASIS(D,BASIS,DD,GLOBAL_FRAME)
    ENDIF



    !    IF(.NOT.ASSOCIATED(P%PARENT_CHART)) THEN  ! ONLY TRANSLATES ORIGINAL OTHERWISE
    ! THEY WILL TRANSLATE MORE THAN ONCE

    IF(ASSOCIATED(P%CHART)) THEN
       IF(ASSOCIATED(P%CHART%F)) THEN
          CALL TRANSLATE_FRAME(P%CHART%F,D,ORDER,BASIS)

          IF(ASSOCIATED(P%MAG%P%F)) THEN
             CALL TRANSLATE_FRAME(P%MAG%P%F,D,ORDER,BASIS)
             P%MAGP%P%F=P%MAG%P%F
          ENDIF
       ENDIF
    ENDIF

    !    ENDIF

    if(associated(R%mag%doko).and.associated(p,r%mag%parent_fibre)) then
       dk=>R%mag%doko
       do while(associated(dk))  !!! PATCH TO DOKO'S  IF CREATED USING DNA I.E. APPEND_POINT
          pP=> dk%parent_fibre
          IF(ASSOCIATED(PP%T1)) THEN
             IF(ASSOCIATED(PP%T1%A)) THEN
                T=>PP%T1
                DO WHILE(.NOT.ASSOCIATED(PP%T2,T))
                   CALL GEO_TRA(T%A,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
                   CALL GEO_TRA(T%B,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
                   T=>T%NEXT
                ENDDO
                CALL GEO_TRA(T%A,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
                CALL GEO_TRA(T%B,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
             ENDIF
          ENDIF

          dk=>dk%next
       enddo
    endif
    IF(ASSOCIATED(R%T1)) THEN
       IF(ASSOCIATED(R%T1%A)) THEN
          T=>P%T1
          DO WHILE(.NOT.ASSOCIATED(R%T2,T))
             CALL GEO_TRA(T%A,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
             CALL GEO_TRA(T%B,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
             T=>T%NEXT
          ENDDO
          CALL GEO_TRA(T%A,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
          CALL GEO_TRA(T%B,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
       ENDIF
    ENDIF

    if(dog) then    ! IF  DOGIRDER IS SET TO TRUE
       if(associated(r%mag%GIRDER_FRAME)) then
          caf=>r%mag
          CALL GEO_TRA(CAF%GIRDER_FRAME%a,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
          CALL GEO_TRA(CAF%GIRDER_FRAME%b,GLOBAL_FRAME,DD,IORDER)    ! A= A +I D*ENT
       endif
    endif
  END SUBROUTINE TRANSLATE_Fibre



  SUBROUTINE  ROTATE_LAYOUT(R,OMEGA,Ang,I1,I2,ORDER,BASIS) ! ROTATES A LAYOUT AROUND OMEGA BY A(3)  IN STANDARD PTC ORDER
    ! INVERSE => ORDER=-1   USING GLOBAL FRAME
    IMPLICIT NONE
    TYPE (LAYOUT),INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: OMEGA(3),Ang(3)
    TYPE(FIBRE), POINTER::P
    REAL(DP) OMEGAT(3)
    INTEGER I,IORDER,I11,I22
    INTEGER, OPTIONAL :: ORDER,I1,I2
    REAL(DP), OPTIONAL, INTENT(IN):: BASIS(3,3)
    real(dp) basist(3,3)
    ! THIS ROUTINE ROTATES THE ENTIRE LINE BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE THE ANGLES A(3) AND THE POINT OMEGA AROUND WHICH THE
    ! ROTATION HAPPENS
    ! OMEGA DEFINED IN THAT BASIS
    ! ANGLE AS WELL
    OMEGAT=OMEGA
    P=>R%START
    IORDER=1
    I11=1
    I22=R%N
    IF(PRESENT(ORDER)) IORDER=ORDER
    IF(PRESENT(I1))  I11=I1
    IF(PRESENT(I2))  I22=I2

    BASIST=GLOBAL_FRAME            ! NECESSARY SINCE BASIS CAN CHANGE DURING THE CALCULATION ASSUMING A POINTER IS PASSED
    IF(PRESENT(BASIS)) BASIST=BASIS


    DO I=1,I11-1
       P=>P%NEXT
    ENDDO


    DO I=1,I22-I11+1
       CALL ROTATE_FIBRE(P,OMEGA,Ang,ORDER,BASIST,dogirder=my_true) ! IF  DOGIRDER IS SET TO TRUE, DOES GIRDER
       P=>P%NEXT
    ENDDO

  END SUBROUTINE ROTATE_LAYOUT

  SUBROUTINE  ROTATE_FIBRE(R,OMEGA,Ang,ORDER,BASIS,dogirder) ! ROTATES A FIBRE AROUND OMEGA BY A(3)  IN STANDARD PTC ORDER
    ! INVERSE => ORDER=-1   USING GLOBAL FRAME IF BASIS NOT SPECIFIED
    IMPLICIT NONE
    TYPE (FIBRE),TARGET,INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: OMEGA(3),Ang(3)
    TYPE(FIBRE), POINTER::P,pp
    REAL(DP) OMEGAT(3)
    INTEGER IORDER
    INTEGER, OPTIONAL :: ORDER
    REAL(DP), OPTIONAL, INTENT(IN):: BASIS(3,3)
    real(dp) basist(3,3),D(3)
    TYPE(INTEGRATION_NODE), POINTER :: T
    TYPE(element), POINTER::caf
    logical(lp), OPTIONAL, INTENT(IN) :: dogirder
    logical(lp) dog
    type(fibre_appearance), pointer :: dk

    dog=.false.
    if(present(dogirder)) dog=dogirder


    OMEGAT=OMEGA
    P=>R
    Pp=>R
    IORDER=1
    IF(PRESENT(ORDER)) IORDER=ORDER
    BASIST=GLOBAL_FRAME            ! NECESSARY SINCE BASIS CAN CHANGE DURING THE CALCULATION ASSUMING A POINTER IS PASSED
    IF(PRESENT(BASIS)) BASIST=BASIS




    !    IF(.NOT.ASSOCIATED(P%PARENT_CHART)) THEN  ! ONLY ROTATES ORIGINAL OTHERWISE
    IF(ASSOCIATED(P%CHART)) THEN
       IF(ASSOCIATED(P%CHART%F)) THEN
          ! THEY WILL ROTATE MORE THAN ONCE
          CALL ROTATE_FRAME(P%CHART%F, OMEGAT,Ang,IORDER,BASIST)


          IF(ASSOCIATED(P%MAG%P%F)) THEN

             CALL ROTATE_FRAME(P%MAG%P%F, OMEGAT,Ang,IORDER,BASIST)
             P%MAGP%P%F=P%MAG%P%F
          ENDIF
       ENDIF
    ENDIF
    !    ENDIF

    if(associated(R%mag%doko).and.associated(p,r%mag%parent_fibre)) then
       dk=>R%mag%doko
       do while(associated(dk))  !!! PATCH TO DOKO'S  IF CREATED USING DNA I.E. APPEND_POINT
          pP=> dk%parent_fibre
          IF(ASSOCIATED(pp%T1)) THEN
             IF(ASSOCIATED(pp%T1%A)) THEN
                T=>pp%T1
                DO WHILE(.NOT.ASSOCIATED(pp%T2,T))
                   D=T%A-OMEGAT
                   CALL GEO_ROT(T%ENT,D,ANG,IORDER,BASIST)
                   T%A=OMEGAT+D
                   D=T%B-OMEGAT     ! ERROR BEFORE  2008.5.20
                   CALL GEO_ROT(T%EXI,D,ANG,IORDER,BASIST)
                   D=T%B-OMEGAT
                   T=>T%NEXT
                ENDDO
                D=T%A-OMEGAT
                CALL GEO_ROT(T%ENT,D,ANG,IORDER,BASIST)
                T%A=OMEGAT+D
                D=T%B-OMEGAT     ! ERROR BEFORE  2008.5.20
                CALL GEO_ROT(T%EXI,D,ANG,IORDER,BASIST)
                D=T%B-OMEGAT
             ENDIF
          ENDIF

          dk=>dk%next
       enddo
    endif
    IF(ASSOCIATED(R%T1)) THEN
       IF(ASSOCIATED(R%T1%A)) THEN
          T=>R%T1
          DO WHILE(.NOT.ASSOCIATED(R%T2,T))
             D=T%A-OMEGAT
             CALL GEO_ROT(T%ENT,D,ANG,IORDER,BASIST)
             T%A=OMEGAT+D
             D=T%B-OMEGAT     ! ERROR BEFORE  2008.5.20
             CALL GEO_ROT(T%EXI,D,ANG,IORDER,BASIST)
             D=T%B-OMEGAT
             T=>T%NEXT
          ENDDO
          D=T%A-OMEGAT
          CALL GEO_ROT(T%ENT,D,ANG,IORDER,BASIST)
          T%A=OMEGAT+D
          D=T%B-OMEGAT     ! ERROR BEFORE  2008.5.20
          CALL GEO_ROT(T%EXI,D,ANG,IORDER,BASIST)
          D=T%B-OMEGAT
       ENDIF
    ENDIF

    if(dog) then   ! IF  DOGIRDER IS SET TO TRUE, TRANSLATE GIRDER FRAME
       if(associated(r%mag%GIRDER_FRAME)) then
          caf=>r%mag
          D=CAF%GIRDER_FRAME%A-OMEGAT
          CALL GEO_ROT(CAF%GIRDER_FRAME%ENT,D,ANG,IORDER,BASIST)
          CAF%GIRDER_FRAME%A=OMEGAT+D
          D=CAF%GIRDER_FRAME%B-OMEGAT
          CALL GEO_ROT(CAF%GIRDER_FRAME%EXI,D,ANG,IORDER,BASIST)
          CAF%GIRDER_FRAME%B=OMEGAT+D
       endif
    endif


  END SUBROUTINE ROTATE_FIBRE






  SUBROUTINE  FIBRE_BL(S2,S1) ! PUTS A NEW MULTIPOLE BLOCK INTO FIBRE. EXTENDS ELEMENT(P) ROUTINES TO FIBRES
    IMPLICIT NONE
    TYPE (MUL_BLOCK),INTENT(IN):: S1
    TYPE(FIBRE),INTENT(INOUT):: S2

    S2%MAG=S1
    S2%MAGP=S1

  END   SUBROUTINE  FIBRE_BL

  SUBROUTINE  BL_FIBRE(S2,S1) ! SUCKS THE MULTIPOLE OUT LOOKING AT ELEMENT
    IMPLICIT NONE
    TYPE (FIBRE),INTENT(IN):: S1
    TYPE(MUL_BLOCK),INTENT(INOUT):: S2

    S2=S1%MAG


  END   SUBROUTINE  BL_FIBRE


  !RECURSIVE
  SUBROUTINE SURVEY_EXIST_PLANAR_IJ(PLAN,I1,I2,ENT,A) ! STANDARD SURVEY FROM FIBRE #I1 TO #I2
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: PLAN
    TYPE (FIBRE), POINTER :: C
    TYPE (PATCH), POINTER :: P
    REAL(DP),OPTIONAL, INTENT(INOUT) :: A(3),ENT(3,3)
    INTEGER , INTENT(IN)::I1,I2
    INTEGER I,i22
    REAL(DP) AT(3),ENTT(3,3),NORM
    LOGICAL(LP) SKIP

    SKIP=.FALSE.

    NULLIFY(C);

    if(i2>=i1) then
       i22=i2
    else
       i22=PLAN%n+i2
    endif


    CALL MOVE_TO(PLAN,C,MOD_N(I1,PLAN%N))


    IF((PRESENT(ENT).AND.(.NOT.PRESENT(A))).OR.(PRESENT(A).AND.(.NOT.PRESENT(ENT)))) THEN
       !w_p=0
       !w_p%NC=2
       !w_p%FC='(2(1X,A72,/),(1X,A72))'
       write(6,*) " BEWARE : ENT AND A  "
       write(6,*) " MUST BOTH BE PRESENT OR ABSENT"
       ! call !write_e(100)
    ELSEIF(PRESENT(ENT)) THEN
       ENTT=ENT
       AT=A
    ELSE
       IF(ASSOCIATED(C%CHART%F)) THEN
          IF(C%DIR==1) THEN
             ENTT=C%CHART%F%ENT
             AT=C%CHART%F%A
          ELSE
             ENTT=C%CHART%F%EXI
             AT=C%CHART%F%B
          ENDIF
       ELSE
          write(6,*) " No charts "
          STOP 888
       ENDIF
       IF(ASSOCIATED(C%PATCH)) THEN
          P=>C%PATCH
          IF(P%PATCH/=0) THEN
             NORM=0.0_dp
             DO I=1,3
                NORM=NORM+ABS(P%A_ANG(I))
                NORM=NORM+ABS(P%A_D(I))
             ENDDO
             NORM=NORM+ABS(P%A_X1-1)+ABS(P%A_X2-1)
             IF(NORM/=0.0_dp) THEN
                WRITE(6,*) " NORM IN SURVEY ", NORM
                WRITE(6,*) " THE SURVEY SKIPS THE FIRST PATCH"
                WRITE(6,*) " IT IS NOT A SELF-CHECK ANYMORE"
                SKIP=.TRUE.
             ENDIF
          ENDIF
       ENDIF

    ENDIF



    I=I1

    DO  WHILE(I<I22.AND.ASSOCIATED(C))

       CALL survey_FIBRE(C,ENTT,AT,SKIP)
       SKIP=.FALSE.
       C=>C%NEXT
       I=I+1
    ENDDO


    !    IF(PRESENT(ENT)) THEN
    !       ENT=ENTT
    !       A=AT
    !    ENDIF



  END SUBROUTINE SURVEY_EXIST_PLANAR_IJ


  !recursive
  SUBROUTINE SURVEY_FIBRE(C,ENT,A,nogirder,SKIP,E_IN)   !,MAGNETFRAME
    !changed
    ! SURVEYS A SINGLE ELEMENT FILLS IN CHART AND MAGNET_CHART; LOCATES ORIGIN AT THE ENTRANCE OR EXIT
    IMPLICIT NONE
    TYPE(FIBRE), TARGET , INTENT(INOUT):: C
    !    TYPE(MAGNET_FRAME), OPTIONAL :: MAGNETFRAME
    TYPE(INNER_FRAME), OPTIONAL :: E_IN
    logical(lp), OPTIONAL :: nogirder
    LOGICAL(LP), OPTIONAL :: SKIP
    REAL(DP), INTENT(INOUT)  :: ENT(3,3),A(3)
    REAL(DP) D(3),ANG(3),DT(3)
    LOGICAL(LP) SEL,SKIPT,dog
    TYPE (PATCH), POINTER :: P
    REAL(DP)dg1(3),ag1(3),ENT0(3,3)
    REAL(DP)dg2(3),ag2(3)
    dog=.true.
    if(present(nogirder)) dog=.not.nogirder     ! IF  DOGIRDER IS SET TO TRUE, MOVES GIRDER FRAME DURING SURVEY
!!!  RECORDS RELATIVE POSITION OF GIRDER
    IF(ASSOCIATED(C%MAG%GIRDER_FRAME).and.dog) THEN
       call FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            C%MAG%GIRDER_FRAME%a, C%MAG%GIRDER_FRAME%ent,dg1,ag1)
       call FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            C%MAG%GIRDER_FRAME%b, C%MAG%GIRDER_FRAME%exi,dg2,ag2)
    ENDIF

    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,-6)
    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,-5)
    SKIPT=.FALSE.

    IF(PRESENT(SKIP)) SKIPT=SKIP
    SEL=.FALSE.
    IF(ASSOCIATED(C%CHART)) THEN
       SEL=.FALSE.
       IF(ASSOCIATED(C%CHART%F)) SEL=.TRUE.
    ENDIF

    !        IF(.NOT.SEL) THEN !
    !        ENDIF
    IF(.NOT.SKIPT) THEN
       IF(ASSOCIATED(C%PATCH)) THEN
          P=>C%PATCH
          IF(P%PATCH==1.or.P%PATCH==2.or.P%PATCH==3) THEN
             ANG=0.0_dp
             ANG=P%A_ANG ;

             DT=P%A_D
             ! seems wrong !
             IF(P%A_X1*P%A_X2<0) ANG(1)=ANG(1)+PI
             IF(P%A_X2<0) then
                ANG(2)=-ANG(2)
                ANG(3)=-ANG(3)
                DT(2)=-DT(2)
                DT(3)=-DT(3)
             endif

             D=0.0_dp
             ent0=ent
             CALL GEO_ROT(ENT0,D,ANG,1,ENT)  ! for frank's flags
             ent=ent0

             CALL GEO_TRA(A,ENT,DT,1)

          ENDIF
       ENDIF

       IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,-4)
       IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,-3)
    ENDIF

    IF(C%DIR==1) THEN
       C%CHART%F%ENT=ENT
       C%CHART%F%A=A
    ELSE
       C%CHART%F%EXI=ENT
       C%CHART%F%B=A
    ENDIF


    CALL SURVEY_NO_PATCH(C,E_IN=E_IN)


    IF(C%DIR==1) THEN
       ENT=C%CHART%F%EXI
       A=C%CHART%F%B
    ELSE
       ENT=C%CHART%F%ENT
       A=C%CHART%F%A
    ENDIF

    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,E_IN%NST-3)
    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,E_IN%NST-2)


    IF(ASSOCIATED(C%PATCH)) THEN
       IF(P%PATCH==1.or.P%PATCH==2.or.P%PATCH==3) THEN
          ANG=0.0_dp
          ANG=P%B_ANG ;

          DT=P%B_D
          IF(P%B_X1*P%B_X2<0) ANG(1)=ANG(1)+PI
          IF(P%B_X2<0) then
             ANG(2)=-ANG(2)
             ANG(3)=-ANG(3)
             DT(2)=-DT(2)
             DT(3)=-DT(3)
          endif

          D=0.0_dp
          ent0=ent
          CALL GEO_ROT(ENT0,D,ANG,1,ENT)  ! for frank's flags
          ent=ent0

          CALL GEO_TRA(A,ENT,DT,1)

       ENDIF
    ENDIF
    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,E_IN%NST-1)
    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,E_IN%NST)
    IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT,A,-7)

!!!  PLACES INTO THE COMPUTED RELATIVE POSITION OF GIRDER
    IF(ASSOCIATED(C%MAG%GIRDER_FRAME).and.dog) THEN    ! IF  DOGIRDER IS SET TO TRUE, MOVES GIRDER FRAME DURING SURVEY
       call INVERSE_FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            dg1,ag1,C%MAG%GIRDER_FRAME%a, C%MAG%GIRDER_FRAME%ent)
       call INVERSE_FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            dg2,ag2,C%MAG%GIRDER_FRAME%b, C%MAG%GIRDER_FRAME%exi)
    ENDIF

  END SUBROUTINE SURVEY_FIBRE







  !RECURSIVE
  SUBROUTINE SURVEY_EXIST_PLANAR_L_NEW(PLAN,ENT,A) ! CALLS ABOVE ROUTINE FROM FIBRE #1 TO #PLAN%N : STANDARD SURVEY
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: PLAN
    REAL(DP),OPTIONAL, INTENT(INOUT) :: A(3),ENT(3,3)

    CALL survey(PLAN,1,ENT,A)

  END SUBROUTINE SURVEY_EXIST_PLANAR_L_NEW

  !RECURSIVE
  SUBROUTINE SURVEY_EXIST_PLANAR_I(PLAN,I1,ENT,A) ! STANDARD SURVEY FROM FIBRE #I1 TO #I2
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: PLAN
    REAL(DP),OPTIONAL, INTENT(INOUT) :: A(3),ENT(3,3)
    INTEGER , INTENT(IN)::I1
    INTEGER I2
    I2=PLAN%N+I1

    CALL survey(PLAN,I1,I2,ENT,A)

  END SUBROUTINE SURVEY_EXIST_PLANAR_I


  !RECURSIVE
  SUBROUTINE SURVEY_NO_PATCH(C,MAGNETFRAME,E_IN)
    !changed
    ! SURVEYS A SINGLE ELEMENT FILLS IN CHART AND MAGNET_CHART; LOCATES ORIGIN AT THE ENTRANCE OR EXIT
    IMPLICIT NONE
    TYPE(FIBRE), TARGET , INTENT(INOUT):: C
    TYPE(MAGNET_FRAME), TARGET, OPTIONAL :: MAGNETFRAME
    TYPE(INNER_FRAME), OPTIONAL :: E_IN
    LOGICAL(LP) SEL

    IF(.NOT.ASSOCIATED(C%CHART)) THEN
       RETURN
    ENDIF


    SEL=.FALSE.
    IF(ASSOCIATED(C%CHART)) THEN
       SEL=.FALSE.
       IF(ASSOCIATED(C%CHART%F)) SEL=.TRUE.
    ENDIF


    IF(SEL) THEN
       CALL SURVEY_mag(C%CHART,C%MAG,C%DIR,MAGNETFRAME,E_IN)
    ENDIF

    IF(ASSOCIATED(C%MAGP%P%F)) THEN
       C%MAGP%P%F=C%MAG%P%F
    ENDIF


    RETURN

  END SUBROUTINE SURVEY_NO_PATCH




  SUBROUTINE COPY_LAYOUT(R2,R1) ! COPY STANDARD LAYOUT ONLY
    IMPLICIT  NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: R1
    TYPE(LAYOUT),target, INTENT(INOUT):: R2
    TYPE (FIBRE), POINTER :: C
    LOGICAL(LP) :: DONEITT=.TRUE.
    INTEGER I   !  TGV
    NULLIFY(C)
    !    CALL LINE_L(R1,DONEIT)  !  TGV

    IF(ASSOCIATED(R2%N)) CALL KILL(R2)
    CALL SET_UP(R2)

    R2%CLOSED=.FALSE.
    R2%NTHIN=R1%NTHIN
    R2%THIN=R1%THIN
    R2%HARMONIC_NUMBER=R1%HARMONIC_NUMBER
    !    if(associated(r1%parent_universe)) R2%parent_universe=> r1%parent_universe
    C=> R1%START
    !    DO WHILE(ASSOCIATED(C))  !  TGV
    DO I=1,R1%N
       CALL APPEND(R2,C)
       C=>C%NEXT
    ENDDO
    R2%LASTPOS=R2%N
    R2%LAST=>R2%END

    R2%CLOSED=R1%CLOSED
    CALL RING_L(R2,DONEITT)
    !   CALL RING_L(R1,DONEIT)  !  TGV

  END SUBROUTINE COPY_LAYOUT

  SUBROUTINE COPY_LAYOUT_IJ(R1,I,J,R2)  ! COPY PIECES OF A STANDARD LAYOUT FROM FIBRE #I TO #J
    IMPLICIT  NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: R1
    TYPE(LAYOUT),target, INTENT(INOUT):: R2
    INTEGER, INTENT(IN):: I,J
    TYPE (FIBRE), POINTER :: C
    LOGICAL(LP) :: DONEITT=.TRUE.
    INTEGER K
    NULLIFY(C)

    !    CALL LINE_L(R1,DONEIT)   !TGV


    IF(ASSOCIATED(R2%N)) CALL KILL(R2)
    CALL SET_UP(R2)

    R2%CLOSED=.FALSE.
    R2%NTHIN=R1%NTHIN
    R2%THIN=R1%THIN
    !    if(associated(r1%parent_universe)) R2%parent_universe=> r1%parent_universe

    CALL MOVE_TO(R1,C,I)
    K=I
    !    DO WHILE(ASSOCIATED(C).AND.K<=J) !TGV
    DO K=I,J
       CALL APPEND(R2,C)
       !    CALL APPEND(R2,C%MAG)
       !    CALL EQUAL(R2%END%CHART,C%CHART)
       C=>C%NEXT
       !       K=K+1  !TGV
    ENDDO
    R2%LASTPOS=R2%N
    R2%LAST=>R2%END

    R2%CLOSED=R1%CLOSED
    CALL RING_L(R2,DONEITT)
    !   CALL RING_L(R1,DONEIT) !TGV

  END SUBROUTINE COPY_LAYOUT_IJ




  SUBROUTINE COPY_LAYOUT_I(R1,R2) ! COPIES IN THE COPY ORDER RATHER THAN THE LAYOUT ORDER
    IMPLICIT  NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: R1
    TYPE(LAYOUT),target, INTENT(INOUT):: R2

    CALL EQUAL(R2,R1)

  END SUBROUTINE COPY_LAYOUT_I

  SUBROUTINE KILL_PARA_L(R)  ! RESETS ALL THE PARAMETERS IN A LAYOUT : REMOVE POLYMORPHIC KNOBS
    IMPLICIT NONE
    TYPE(LAYOUT),target,INTENT(INOUT):: R
    TYPE (FIBRE), POINTER :: C
    INTEGER I
    c_%np_pol=0
    NULLIFY(C)

    !    CALL LINE_L(R,DONEIT)  ! TGV

    C=>R%START
    DO I=1,R%N
       !    DO WHILE(ASSOCIATED(C))  ! TGV
       if(mfpolbloc/=0)  call ELp_POL_print(C%MAGP)
       CALL RESET31(C%MAGP)
       C=>C%NEXT
    ENDDO
    !    CALL RING_L(R,DONEIT)  ! TGV
  END       SUBROUTINE KILL_PARA_L

  SUBROUTINE  FIBRE_POL(S2,S1)    !  SET POLYMORPH IN A FIBRE
    IMPLICIT NONE
    TYPE (POL_BLOCK),INTENT(IN):: S1
    TYPE(FIBRE),INTENT(INOUT):: S2
    S2%MAGP=S1
  END SUBROUTINE  FIBRE_POL

  SUBROUTINE  EL_POL_force(S2,S1)    !  SET POLYMORPH IN A FIBRE UNCONDITIONALLY
    IMPLICIT NONE
    TYPE (POL_BLOCK),INTENT(IN):: S1
    TYPE(FIBRE),INTENT(INOUT):: S2
    call ELp_POL_force(S2%MAGP,S1)
  END SUBROUTINE  EL_POL_force

  SUBROUTINE SCAN_FOR_POLYMORPHS(R,B)   !  SET POLYMORPH IN A FULL LAYOUT ONLY IF THE MAGNET IS A PRIMITIVE PARENT
    IMPLICIT  NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: R
    TYPE(POL_BLOCK), INTENT(IN):: B

    TYPE (FIBRE), POINTER :: C
    INTEGER I

    NULLIFY(C)
    !    CALL LINE_L(R,DONEIT)  ! TGV
    C=>R%START

    DO I=1,R%N
       !    DO WHILE(ASSOCIATED(C))  ! TGV
       !       IF(.NOT.ASSOCIATED(C%PARENT_MAG)) THEN
       C%MAGP=B
       !       ENDIF
       C=>C%NEXT
    ENDDO
    !    CALL RING_L(R,DONEIT)


  END SUBROUTINE SCAN_FOR_POLYMORPHS

  SUBROUTINE EL_TO_ELP_L(R)  ! COPY ALL PRIMITIVES ELEMENT INTO ELEMENTP
    IMPLICIT  NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: R
    TYPE (FIBRE), POINTER :: C
    INTEGER I

    NULLIFY(C)
    !    CALL LINE_L(R,DONEIT)  ! TGV

    C=>R%START
    DO I=1,R%N
       !    DO   WHILE(ASSOCIATED(C))  ! TGV
       !       IF(.NOT.ASSOCIATED(C%PARENT_MAG)) CALL COPY(C%MAG,C%MAGP)
       CALL COPY(C%MAG,C%MAGP)
       C=>C%NEXT
    ENDDO

    !    CALL RING_L(R,DONEIT)  ! TGV

  END SUBROUTINE EL_TO_ELP_L

  SUBROUTINE ELP_TO_EL_L(R) ! COPY ALL PRIMITIVES ELEMENTP INTO ELEMENT
    IMPLICIT  NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: R
    TYPE (FIBRE), POINTER :: C
    INTEGER I
    NULLIFY(C)

    !    CALL LINE_L(R,DONEIT) !  ! TGV
    C=>R%START
    !    DO   WHILE(ASSOCIATED(C))  ! TGV
    DO I=1,R%N
       !       IF(.NOT.ASSOCIATED(C%PARENT_MAG)) CALL COPY(C%MAGP,C%MAG)
       CALL COPY(C%MAGP,C%MAG)

       C=>C%NEXT
    ENDDO

    !    CALL RING_L(R,DONEIT)  ! TGV
  END SUBROUTINE ELP_TO_EL_L


  SUBROUTINE fill_survey_ONE_FIBRE(R)
    ! THIS SUBROUTINE ALLOCATES NODE FRAMES IF NEEDED
    ! IT SURVEYS THE NODES USING THE OLD REAL WORMS
    ! SHOULD BE CALLED AFTER MISALIGNMENTS OR MOVING PART OF LATTICE

    IMPLICIT NONE
    type(FIBRE),target:: r
    type(fibre), pointer ::c
    type(INTEGRATION_NODE), pointer ::t
    type(worm) vers
    integer ic,j
    real(dp) x(6),ent(3,3),a(3)
    !    INTEGER, TARGET :: CHARGE
    LOGICAL(LP) APER

    aper=APERTURE_FLAG
    APERTURE_FLAG=.FALSE.




    C=>R

    CALL ALLOC(vers,r)

    CALL XFRAME(vers%E,C%chart%f%ent,C%chart%f%A,-7)  ! initializes the survey part of worm
    vers%E%L(-1)=0.d0 !Starts beam line at z=0   fake distance along ld for cheap work

    !    do k=1,r%n
    x=0.0_dp
    CALL TRACK_FIBRE_RR(C,x,default,vers)

    t=>c%t1
    j=-6
    call gMID(vers,x,j)
    call G_FRAME(vers%e,ENT,A,j)
    t%ent=ent
    t%a=a

    t=>t%next
    if(t%cas/=case1) then
       write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%cas
       stop 665
    endif
    j=vers%POS(2)
    call gMID(vers,x,j)
    call G_FRAME(vers%e,ENT,A,j)
    t%ent=ent
    t%a=a
    t%previous%exi=ent
    t%previous%b=a
    t=>t%next
    ic=0
    DO J=vers%POS(2)+1,vers%POS(3)-1     ! pos(2)+1 to pos(3)-1 inside the magnet

       ic=ic+1

       call gMID(vers,x,j)
       call G_FRAME(vers%e,ENT,A,j)

       if(j/=vers%POS(2)+1) then
          t%previous%exi=ent
          t%previous%b=a
          if(t%previous%cas/=case0) then
             write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%previous%cas
             stop 666
          endif
       else
          t%previous%exi=ent
          t%previous%b=a
          if(t%previous%cas/=case1) then
             write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%previous%cas
             stop 664
          endif
       endif

       if(j/=vers%POS(3)-1) then
          t%ent=ent
          t%a=a
          if(t%cas/=case0) then
             write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%cas
             stop 666
          endif
       else
          t%ent=ent
          t%a=a
          if(t%cas/=case2) then
             write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%cas
             write(6,*)t%POS,T%PARENT_FIBRE%MAG%NAME
             write(6,*)T%PARENT_FIBRE%T1%POS,T%PARENT_FIBRE%T2%POS
             stop 668
          endif
       endif



       t=>t%next
    enddo
    j=vers%POS(3)
    call gMID(vers,x,j)
    call G_FRAME(vers%e,ENT,A,j)
    t%previous%exi=ent
    t%previous%b=a
    t%ent=ent
    t%a=a
    if(t%previous%cas/=case2) then
       write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%cas
       stop 669
    endif
    !      t=>t%next

    j=vers%nst
    call gMID(vers,x,j)
    call G_FRAME(vers%e,ENT,A,j)

    t%exi=ent
    t%b=a

    if(t%cas/=casep2) then
       write(6,*)" error in fill_survey_data_in_NODE_LAYOUT",j,t%cas
       stop 670
    endif


    if(ic/=c%mag%p%nst+1) then
       write(6,*)" error in fill_survey_data_in_NODE_LAYOUT"
       write(6,*) ic,c%mag%name,c%mag%p%nst
       stop 888
    endif
    c=>c%next
    !    enddo

    CALL kill(vers)

    APERTURE_FLAG=aper

  end  subroutine fill_survey_ONE_FIBRE

  SUBROUTINE TRACK_FIBRE_RR(C,X,K,X_IN)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: X(6)
    TYPE(WORM), OPTIONAL,INTENT(INOUT):: X_IN
    !    INTEGER,optional, target, INTENT(IN) :: CHARGE
    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER(2) PATCHT,PATCHG,PATCHE
    REAL(DP) ENT(3,3), A(3)
    !    integer,target :: charge1




    IF(PRESENT(X_IN)) then
       X_IN%F=>c ; X_IN%E%F=>C; X_IN%NST=>X_IN%E%NST;
    endif

    ! DIRECTIONAL VARIABLE
    C%MAG%P%DIR=>C%DIR
    !    if(present(charge)) then
    !       C%MAG%P%CHARGE=>CHARGE
    !    else
    !       charge1=1
    C%MAG%P%CHARGE=>C%charge
    !    endif
    !
    !    IF(.NOT.CHECK_STABLE) CHECK_STABLE=.TRUE.
    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-6)
       X_IN%POS(1)=X_IN%nst
    endif

    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-5)

    ! The chart frame of reference is located here implicitely
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-4)
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-3)

    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,-2)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,-1)
       X_IN%POS(2)=X_IN%nst
    endif

    CALL TRACK(C%MAG,X,K,X_IN)

    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(3)=X_IN%nst
    endif

    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)
    ! The magnet frame of reference is located here implicitely before misalignments
    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    IF(PRESENT(X_IN)) CALL XMID(X_IN,X,X_IN%nst+1)

    ! The CHART frame of reference is located here implicitely

    IF(PRESENT(X_IN)) then
       CALL XMID(X_IN,X,X_IN%nst+1)
       X_IN%POS(4)=X_IN%nst
    endif

    IF(PRESENT(X_IN))  THEN
       IF(X_IN%E%DO_SURVEY) THEN
          CALL G_FRAME(X_IN%E,ENT,A,-7)
          CALL  SURVEY(C,ENT,A,E_IN=X_IN%E)
       ELSE
          CALL SURVEY_INNER_MAG(X_IN%E)
       ENDIF
    ENDIF



    nullify(C%MAG%P%DIR)
    nullify(C%MAG%P%CHARGE)
  END SUBROUTINE TRACK_FIBRE_RR


END  MODULE        S_FAMILY
