!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
MODULE S_FAMILY
  use S_FIBRE_BUNDLE
  IMPLICIT NONE

  ! LINKED LIST
  PRIVATE SURVEY_EXIST_PLANAR_L_new ,SURVEY_EXIST_PLANAR_ij_new,survey_one
  PRIVATE COPY_LAYOUT,COPY_LAYOUT_I,kill_para_L
  private fibre_WORK,fibre_POL,fibre_BL,ADDP_ANBN,WORK_fibre,BL_fibre
  PRIVATE TRANS_D,rot_a,COPY_LAYOUT_ij


  INTERFACE EL_TO_ELP
     !LINKED
     MODULE PROCEDURE EL_TO_ELP_L
  END INTERFACE

  INTERFACE ELP_TO_EL
     !LINKED
     MODULE PROCEDURE ELP_TO_EL_L
  END INTERFACE

  INTERFACE SURVEY
     ! link list
     MODULE PROCEDURE SURVEY_one
     MODULE PROCEDURE SURVEY_EXIST_PLANAR_L_new
     MODULE PROCEDURE SURVEY_EXIST_PLANAR_IJ_new
  END INTERFACE


  INTERFACE kill_para
     MODULE PROCEDURE kill_para_L
  END INTERFACE

  INTERFACE ADD
     MODULE PROCEDURE ADDP_ANBN
  END INTERFACE

  INTERFACE assignment (=)
     ! linked
     MODULE PROCEDURE scan_for_polymorphs
     MODULE PROCEDURE fibre_WORK
     MODULE PROCEDURE misalign_fibre
     MODULE PROCEDURE fibre_POL
     MODULE PROCEDURE fibre_BL
     MODULE PROCEDURE BL_fibre
     MODULE PROCEDURE WORK_fibre
  end  INTERFACE

  INTERFACE EQUAL
     ! LINKED
     MODULE PROCEDURE COPY_LAYOUT
  end  INTERFACE

  INTERFACE copy
     ! LINKED
     MODULE PROCEDURE COPY_LAYOUT_I
     MODULE PROCEDURE COPY_LAYOUT_ij
  end  INTERFACE

  INTERFACE TRANS
     ! LINKED
     MODULE PROCEDURE TRANS_D
  end  INTERFACE

  INTERFACE TRANSLATE
     ! LINKED
     MODULE PROCEDURE TRANS_D
  end  INTERFACE

  INTERFACE rotate
     ! LINKED
     MODULE PROCEDURE rot_a
  end  INTERFACE

  INTERFACE rotation
     ! LINKED
     MODULE PROCEDURE rot_a
  end  INTERFACE



CONTAINS

  !NEW

  subroutine LOCATE_FIBRE(R,PIN,I)
    implicit none
    type(layout), intent(in) :: r
    type(fibre), pointer:: p,PIN
    integer, intent(inOUT) :: i
    p=>r%start
    do i=1,r%n
       IF(ASSOCIATED(PIN,P) ) EXIT
       p=>p%next
    enddo
  end subroutine LOCATE_FIBRE

  subroutine get_length(R,L)
    implicit none
    type(layout), intent(in) :: r
    real(dp), intent(out) :: L
    type(fibre), pointer:: p
    integer i
    p=>r%start
    L=0.d0
    do i=1,r%n
       L=L+p%mag%p%ld
       p=>p%next
    enddo
  end subroutine get_length

  subroutine get_freq(R,freq)
    implicit none
    type(layout), intent(in) :: r
    real(dp), intent(out) :: freq
    type(fibre), pointer:: p
    integer i
    p=>r%start
    freq=0.d0
    do i=1,r%n
       if(associated(p%mag%freq)) then
          if(p%mag%freq/=0.D0) THEN
             FREQ=p%mag%freq
          ENDIF
       endif
       p=>p%next
    enddo
  end subroutine get_freq

  subroutine get_ALL(R,freq,VOLT,PHAS)
    implicit none
    type(layout), intent(in) :: r
    real(dp), intent(out) :: freq,VOLT,PHAS
    type(fibre), pointer:: p
    integer i
    p=>r%start
    freq=0.d0;VOLT=0.D0;PHAS=0.D0;
    do i=1,r%n
       if(associated(p%mag%freq)) then
          if(p%mag%freq/=0.D0) THEN
             FREQ=twopi*p%mag%freq/CLIGHT
             VOLT=-P%MAG%volt*c_1d_3/P%MAG%P%P0C
             PHAS=P%MAG%PHAS
          ENDIF
       endif
       p=>p%next
    enddo
  end subroutine get_ALL

  subroutine set_freq(R,freq)
    implicit none
    type(layout), intent(inout) :: r
    real(dp), intent(in) :: freq
    type(fibre), pointer:: p
    integer i
    p=>r%start
    do i=1,r%n
       if(associated(p%mag%freq)) then
          if(p%mag%freq/=0.D0) THEN
             p%mag%freq=freq
             p%magp%freq=freq
          ENDIF
       endif
       p=>p%next
    enddo
  end subroutine set_freq

  subroutine add_freq(R,freq)
    implicit none
    type(layout), intent(inout) :: r
    real(dp), intent(in) :: freq
    type(fibre), pointer:: p
    integer i
    p=>r%start
    do i=1,r%n
       if(associated(p%mag%freq)) then
          if(p%mag%freq/=0.D0) THEN
             p%mag%freq=p%mag%freq+freq
             p%magp%freq=p%magp%freq+freq
          ENDIF
       endif
       p=>p%next
    enddo
  end subroutine add_freq

  !END NEW
  SUBROUTINE ADDP_ANBN(EL,NM,F,V) ! Extends the ADD routines from the Element(p) to the fibre
    IMPLICIT NONE
    TYPE(fibre), INTENT(INOUT) ::EL
    real(dp), INTENT(IN) ::V
    INTEGER, INTENT(IN) ::NM,F

    CALL ADD(EL%MAG,NM,F,V)
    CALL ADD(EL%MAGP,NM,F,V)

  END SUBROUTINE ADDP_ANBN

  SUBROUTINE  fibre_WORK(S2,S1) ! Changes the Energy of the fibre and turns the energy patch on
    implicit none
    type (WORK),INTENT(IN):: S1
    TYPE(fibre),INTENT(inOUT):: S2

    s2%mag=s1
    s2%magp=s1
    s2%PATCH%energy=.true.

  END SUBROUTINE fibre_WORK

  SUBROUTINE  WORK_fibre(S2,S1)  ! Sucks the energy out of a fibre by loking at ELEMENT
    implicit none
    type (fibre),INTENT(IN):: S1
    TYPE(work),INTENT(inOUT):: S2

    s2=s1%mag
    if(abs(s1%mag%p%p0c-s1%magp%p%p0c)>c_1d_10) then
       w_p=0
       w_p%nc=3
       w_p%fc='(2(1X,a72,/),(1X,a72))'
       w_p%c(1)=" Beware : Element and ElementP seem to have "
       w_p%c(2)=" different reference energies!"
       write(w_p%c(3),'(1x,g20.14,1x,g20.14)')  s1%mag%p%p0c,s1%magp%p%p0c
       call write_e(100)
    endif

  END SUBROUTINE WORK_fibre

  SUBROUTINE  misalign_fibre(S2,S1) ! Misaligns full fibre; fills in chart and Magnet_chart
    use rotation_mis
    implicit none
    real(dp),INTENT(IN):: S1(6)
    TYPE(fibre),INTENT(inOUT):: S2
    integer i

    if(ASSOCIATED(S2%CHART)) THEN
       IF(.NOT.ASSOCIATED(s2%mag%D) ) ALLOCATE(s2%mag%D(3))
       IF(.NOT.ASSOCIATED(s2%mag%R) ) ALLOCATE(s2%mag%R(3))
       IF(.NOT.ASSOCIATED(s2%magp%D)) ALLOCATE(s2%magp%D(3))
       IF(.NOT.ASSOCIATED(s2%magp%R)) ALLOCATE(s2%magp%R(3))
       DO I=1,3
          s2%mag%D(I)=S1(I);   s2%magp%D(I)=S1(I);
          s2%mag%R(I)=S1(3+I); s2%magp%R(I)=S1(3+I);
       ENDDO

       s2%mag%mis=.true.
       s2%magp%mis=.true.

       ! add code here
       CALL FACTORIZE_ROTATION(S1,S2%CHART%L,S2%CHART%ALPHA,S2%CHART%D_IN,   &
            & S2%CHART%ANG_IN,S2%CHART%D_OUT,S2%CHART%ANG_OUT)
       !
       CALL ADJUST_INTERNAL(S2)
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       write(w_p%c(1),'(1x,a39,1x,A16)') " CANNOT MISALIGN THIS FIBRE: NO CHARTS ", S2%MAG%NAME
       call write_e(100)
    ENDIF

  END SUBROUTINE misalign_fibre

  SUBROUTINE  ADJUST_INTERNAL(S2) ! Adjusts frames of magnet_chart on the basis of the misalignments
    IMPLICIT NONE
    TYPE(fibre),INTENT(inOUT):: S2
    if(associated(S2%mag%p%f)) then
       CALL COPY_INTERNAL(S2)
       CALL GEO_ROT(S2%mag%p%f%ent,S2%CHART%ANG_IN,1)
       CALL GEO_TRA(S2%mag%p%f%A,S2%CHART%f%ent,S2%CHART%D_IN,1)
       CALL GEO_TRA(S2%mag%p%f%B,S2%CHART%f%exi,S2%CHART%d_out,-1)
       CALL GEO_ROT(S2%mag%p%f%EXI,S2%CHART%ANG_out,-1)
       CALL GEO_AVE(S2%mag%p%f%ENT,S2%mag%p%f%A,S2%mag%p%f%EXI,   S2%mag%p%f%b,S2%mag%p%f%MID,S2%mag%p%f%o)

       CALL GEO_ROT(S2%magP%p%f%ent,S2%CHART%ANG_IN,1)
       CALL GEO_TRA(S2%magP%p%f%A,S2%CHART%f%ent,S2%CHART%D_IN,1)
       CALL GEO_TRA(S2%magP%p%f%B,S2%CHART%f%exi,S2%CHART%d_out,-1)
       CALL GEO_ROT(S2%magP%p%f%EXI,S2%CHART%ANG_out,-1)
       CALL GEO_AVE(S2%magP%p%f%ENT,S2%magp%p%f%A,S2%magP%p%f%EXI,S2%magp%p%f%b,S2%magP%p%f%MID,S2%magp%p%f%o)
    endif
  END SUBROUTINE  ADJUST_INTERNAL


  SUBROUTINE  COPY_INTERNAL(S2) ! Puts the frames of the chart into those of the magnet_chart
    IMPLICIT NONE
    TYPE(fibre),INTENT(inOUT):: S2

    if(associated(S2%mag%p%f)) then
       S2%mag%p%f%mid(:,:) = S2%CHART%f%mid(:,:)
       S2%mag%p%f%o(:)     = S2%CHART%f%o(:)
       S2%mag%p%f%ent(:,:) = S2%CHART%f%ent(:,:)
       S2%mag%p%f%a(:)     = S2%CHART%f%a(:)
       S2%mag%p%f%exi(:,:) = S2%CHART%f%exi(:,:)
       S2%mag%p%f%b(:)     = S2%CHART%f%b(:)
       S2%magP%p%f%mid(:,:) = S2%CHART%f%mid(:,:)
       S2%magP%p%f%o(:)     = S2%CHART%f%o(:)
       S2%magP%p%f%ent(:,:) = S2%CHART%f%ent(:,:)
       S2%magP%p%f%a(:)     = S2%CHART%f%a(:)
       S2%magP%p%f%exi(:,:) = S2%CHART%f%exi(:,:)
       S2%magP%p%f%b(:)     = S2%CHART%f%b(:)
    endif
  END SUBROUTINE  COPY_INTERNAL



  ! new routines to change layout

  SUBROUTINE  Trans_D(R,D) ! Translates a layout
    implicit none
    type (LAYOUT),INTENT(INOUT):: R
    real(dp),INTENT(IN):: D(3)
    TYPE(FIBRE), POINTER::P
    INTEGER I

    ! THIS ROUTINE TRANSLATE THE ENTIRE LINE BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE D

    P=>R%START

    DO I=1,R%N
       if(.not.associated(p%parent_chart)) then  ! Only translates original otherwise
          ! they will translate more than once

          IF(ASSOCIATED(P%CHART)) THEN
             if(associated(P%CHART%f)) then
                P%CHART%f%A=P%CHART%f%A+D
                P%CHART%f%B=P%CHART%f%B+D
                P%CHART%f%O=P%CHART%f%O+D

                if(associated(P%MAG%p%f)) then
                   P%MAG%p%f%A=P%MAG%p%f%A+D
                   P%MAG%p%f%B=P%MAG%p%f%B+D
                   P%MAG%p%f%O=P%MAG%p%f%O+D
                   P%MAGP%p%f%A=P%MAGP%p%f%A+D
                   P%MAGP%p%f%B=P%MAGP%p%f%B+D
                   P%MAGP%p%f%O=P%MAGP%p%f%O+D
                endif
             endif
          ENDIF

       endif
       P=>P%NEXT
    ENDDO


  END SUBROUTINE Trans_D

  SUBROUTINE  ROT_A(R,OMEGA,A) ! Rotates a layout around OMEGA by A(3)  in standard PTC order
    implicit none
    type (LAYOUT),INTENT(INOUT):: R
    real(dp),INTENT(IN):: OMEGA(3),A(3)
    TYPE(FIBRE), POINTER::P
    real(dp) D(3)
    INTEGER I
    ! THIS ROUTINE ROTATES THE ENTIRE LINE BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE THE ANGLES A(3) AND THE POINT OMEGA AROUND WHICH THE
    ! ROTATION HAPPENS

    P=>R%START

    DO I=1,R%N
       if(.not.associated(p%parent_chart)) then  ! Only rotates original otherwise
          IF(ASSOCIATED(P%CHART)) THEN
             if(associated(P%CHART%f)) then
                ! they will rotate more than once
                D=P%CHART%f%A-OMEGA
                CALL GEO_ROT(D,A)
                P%CHART%f%A=OMEGA+D

                D=P%CHART%f%O-OMEGA
                CALL GEO_ROT(D,A)
                P%CHART%f%O=OMEGA+D

                D=P%CHART%f%B-OMEGA
                CALL GEO_ROT(D,A)
                P%CHART%f%B=OMEGA+D


                CALL GEO_ROT(P%CHART%f%ENT,A,1)
                CALL GEO_ROT(P%CHART%f%MID,A,1)
                CALL GEO_ROT(P%CHART%f%EXI,A,1)


                if(associated(P%MAG%p%f)) then
                   D=P%MAG%p%f%A-OMEGA
                   CALL GEO_ROT(D,A)
                   P%MAG%p%f%A=OMEGA+D; P%MAGP%p%f%A=P%MAG%p%f%A
                   D=P%MAG%p%f%O-OMEGA
                   CALL GEO_ROT(D,A)
                   P%MAG%p%f%O=OMEGA+D; P%MAGP%p%f%O=P%MAG%p%f%O
                   D=P%MAG%p%f%B-OMEGA
                   CALL GEO_ROT(D,A)
                   P%MAG%p%f%B=OMEGA+D; P%MAGP%p%f%B=P%MAG%p%f%B

                   CALL GEO_ROT(P%MAG%p%f%ENT,A,1)
                   CALL GEO_ROT(P%MAG%p%f%MID,A,1)
                   CALL GEO_ROT(P%MAG%p%f%EXI,A,1)
                   CALL GEO_ROT(P%MAGP%p%f%ENT,A,1)
                   CALL GEO_ROT(P%MAGP%p%f%MID,A,1)
                   CALL GEO_ROT(P%MAGP%p%f%EXI,A,1)
                endif
             endif
          ENDIF
       endif

       P=>P%NEXT
    ENDDO

  END SUBROUTINE ROT_A

  SUBROUTINE  ROTATE_FIBRE(EL1,EL0)
    ! PUTS EL1 AT THE END OF EL0 FOR SURVEY PURPOSES. EL0 AND EL1 MUST HAVE SAME EL%DIR (Standard Survey)
    implicit none
    type (FIBRE),INTENT(INOUT):: EL1,EL0

    ! PUTS EL1 AT THE END OF EL0 FOR SURVEY PURPOSES. EL0 AND EL1 MUST HAVE SAME EL%DIR
    IF(ASSOCIATED(EL0%CHART).AND.ASSOCIATED(EL1%CHART)) THEN
       CALL ROTATE_C(EL1%CHART,EL0%CHART,EL0%DIR)
       CALL ADJUST_INTERNAL(EL1)
    ENDIF
  END SUBROUTINE  ROTATE_FIBRE

  !


  SUBROUTINE  fibre_BL(S2,S1) ! Puts a new multipole block into Fibre. Extends element(p) routines to fibres
    implicit none
    type (MUL_BLOCK),INTENT(IN):: S1
    TYPE(fibre),INTENT(inOUT):: S2

    s2%mag=s1
    s2%magp=s1

  END   SUBROUTINE  fibre_BL

  SUBROUTINE  BL_fibre(S2,S1) ! Sucks the multipole out looking at element
    implicit none
    type (fibre),INTENT(IN):: S1
    TYPE(MUL_BLOCK),INTENT(inOUT):: S2

    s2=s1%mag


  END   SUBROUTINE  BL_fibre


  SUBROUTINE SURVEY_EXIST_PLANAR_IJ_new(PLAN,I1,I2) ! standard survey from fibre #i1 to #i2
    IMPLICIT NONE
    TYPE(layout), INTENT(INOUT):: PLAN
    TYPE (fibre), POINTER :: C
    INTEGER , INTENT(IN)::I1,I2
    INTEGER I


    nullify(C);

    CALL move_to( plan,C,i1)

    !if(first) then
    ! c%CHART=1
    ! call survey(c)
    !endif


    c=>plan%START
    DO I=i1,i2-1
       CALL ROTATE_FIBRE(c%NEXT,c)
       c=>c%NEXT
    ENDDO



  END SUBROUTINE SURVEY_EXIST_PLANAR_IJ_new

  SUBROUTINE SURVEY_EXIST_PLANAR_L_new(PLAN) ! Calls above routine from fibre #1 to #plan%n : Standard Survey
    IMPLICIT NONE
    TYPE(layout), INTENT(INOUT):: PLAN

    call survey(plan,1,plan%n)

  END SUBROUTINE SURVEY_EXIST_PLANAR_L_new


  SUBROUTINE SURVEY_one(C) ! Surveys a single element fills in chart and magnet_chart; locates origin at the entrance or exit
    IMPLICIT NONE
    TYPE(FIBRE), TARGET , INTENT(INOUT):: C
    TYPE (CHART), POINTER :: BL,CL
    TYPE (CHART), TARGET :: CHART_i
    INTEGER J
    real(dp) NORM


    Nullify(BL);Nullify(CL);



    IF(.NOT.ASSOCIATED(C%CHART)) RETURN

    CHART_i=1

    CL=> C%CHART  ! CHART OF ELEMENT 1
    CL=1
    BL=> CHART_i


    if(c%dir==1) then
       CL%A_XY=C%MAG%P%tilTd  ! GO CHARTLY
       CL%L=C%MAG%P%LC
       CL%ALPHA=C%MAG%P%LD*C%MAG%P%B0
       if(associated(C%CHART%f)) then     !!!! doing survey
          !              CL%ALPHA=CL%A_XZ
          CALL ROTATE_V( BL,CL,1)

          CL%f%O(:)=half*(CL%f%A(:)+CL%f%B(:))

          CL%f%MID(:,:)=CL%f%ENT(:,:)
          DO J=1,3
             CL%f%MID(1,J)=(CL%f%ENT(1,J)+CL%f%EXI(1,J))
             CL%f%MID(3,J)=(CL%f%ENT(3,J)+CL%f%EXI(3,J))
          ENDDO
          NORM=zero
          DO J=1,3
             NORM=CL%f%MID(1,J)**2+NORM
          ENDDO
          NORM=SQRT(NORM)
          DO J=1,3
             CL%f%MID(1,J)=CL%f%MID(1,J)/NORM
          ENDDO
          NORM=zero
          DO J=1,3
             NORM=CL%f%MID(3,J)**2+NORM
          ENDDO
          NORM=SQRT(NORM)
          DO J=1,3
             CL%f%MID(3,J)=CL%f%MID(3,J)/NORM
          ENDDO

          CALL ADJUST_INTERNAL(C)
       endif  !!!! doing survey

    else
       CL%A_XY=C%MAG%P%tilTd  ! GO CHARTLY
       CL%L=C%MAG%P%LC
       CL%ALPHA=C%MAG%P%LD*C%MAG%P%B0
       !              CL%ALPHA=CL%A_XZ
       if(associated(C%CHART%f)) then     !!!! doing survey
          CALL ROTATE_V( BL,CL,-1)


          CL%f%O(:)=half*(CL%f%A(:)+CL%f%B(:))

          CL%f%MID(:,:)=CL%f%ENT(:,:)
          DO J=1,3
             CL%f%MID(1,J)=(CL%f%ENT(1,J)+CL%f%EXI(1,J))
             CL%f%MID(3,J)=(CL%f%ENT(3,J)+CL%f%EXI(3,J))
          ENDDO
          NORM=zero
          DO J=1,3
             NORM=CL%f%MID(1,J)**2+NORM
          ENDDO
          NORM=SQRT(NORM)
          DO J=1,3
             CL%f%MID(1,J)=CL%f%MID(1,J)/NORM
          ENDDO
          NORM=zero
          DO J=1,3
             NORM=CL%f%MID(3,J)**2+NORM
          ENDDO
          NORM=SQRT(NORM)
          DO J=1,3
             CL%f%MID(3,J)=CL%f%MID(3,J)/NORM
          ENDDO

          CALL ADJUST_INTERNAL(C)
       endif    !!!! doing survey

    endif




  END SUBROUTINE SURVEY_one


  SUBROUTINE COPY_LAYOUT(R2,R1) ! copy standard layout only
    IMPLICIT  NONE
    TYPE(layout), INTENT(INOUT):: R1
    TYPE(layout), INTENT(inOUT):: R2
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    logical(lp) :: doneitt=.true.
    nullify(C)

    CALL LINE_L(R1,doneit)


    IF(ASSOCIATED(R2%N)) CALL kill(R2)
    CALL SET_UP(R2)

    R2%CLOSED=.false.
    R2%NTHIN=R1%NTHIN
    R2%THIN=R1%THIN

    C=> R1%START
    DO WHILE(ASSOCIATED(C))
       CALL APPEND(R2,C)
       C=>C%NEXT
    ENDDO
    R2%lastpos=R2%N
    r2%LAST=>r2%end

    R2%CLOSED=R1%CLOSED
    CALL RING_L(R2,doneitt)
    CALL RING_L(R1,doneit)

  END SUBROUTINE COPY_LAYOUT

  SUBROUTINE COPY_LAYOUT_ij(R1,i,j,R2)  ! copy pieces of a standard layout from fibre #i to #j
    IMPLICIT  NONE
    TYPE(layout), INTENT(INOUT):: R1
    TYPE(layout), INTENT(inOUT):: R2
    integer, INTENT(in):: i,j
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    logical(lp) :: doneitt=.true.
    integer k
    nullify(C)

    CALL LINE_L(R1,doneit)


    IF(ASSOCIATED(R2%N)) CALL kill(R2)
    CALL SET_UP(R2)

    R2%CLOSED=.false.
    R2%NTHIN=R1%NTHIN
    R2%THIN=R1%THIN

    call move_to(r1,c,i)
    k=i
    DO WHILE(ASSOCIATED(C).and.k<=j)
       CALL APPEND(R2,C)
       !    CALL APPEND(R2,C%mag)
       !    CALL EQUAL(R2%end%CHART,C%CHART)
       C=>C%NEXT
       k=k+1
    ENDDO
    R2%lastpos=R2%N
    r2%LAST=>r2%end

    R2%CLOSED=R1%CLOSED
    CALL RING_L(R2,doneitt)
    CALL RING_L(R1,doneit)

  END SUBROUTINE COPY_LAYOUT_ij




  SUBROUTINE COPY_LAYOUT_I(R1,R2) ! Copies in the copy order rather than the Layout order
    IMPLICIT  NONE
    TYPE(layout), INTENT(INOUT):: R1
    TYPE(layout), INTENT(inOUT):: R2

    CALL EQUAL(R2,R1)

  END SUBROUTINE COPY_LAYOUT_I

  SUBROUTINE KILL_para_L(R)  ! Resets all the parameters in a layout : remove polymorphic knobs
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE (fibre), POINTER :: C
    logical(lp) doneit

    nullify(C)

    CALL LINE_L(R,doneit)

    C=>R%START

    DO WHILE(ASSOCIATED(C))
       CALL reset31(C%MAGP)
       C=>C%NEXT
    ENDDO
    CALL RING_L(R,doneit)
  END       SUBROUTINE KILL_para_L

  SUBROUTINE  fibre_POL(S2,S1)    !  Set polymorph in a fibre unconditionally
    implicit none
    type (POL_BLOCK),INTENT(IN):: S1
    TYPE(FIBRE),INTENT(inOUT):: S2
    S2%MAGP=S1
  END SUBROUTINE  fibre_POL

  SUBROUTINE scan_for_polymorphs(R,B)   !  Set polymorph in a full layout only if the magnet is a primitive parent
    IMPLICIT  NONE
    TYPE(layout), INTENT(inOUT):: R
    TYPE(POL_BLOCK), INTENT(IN):: B

    TYPE (fibre), POINTER :: C
    logical(lp) doneit

    nullify(C)
    CALL LINE_L(R,doneit)
    C=>R%START

    DO WHILE(ASSOCIATED(C))
       IF(.NOT.ASSOCIATED(C%PARENT_MAG)) THEN
          C%MAGP=b
       ENDIF
       C=>C%NEXT
    ENDDO
    CALL RING_L(R,doneit)

  END SUBROUTINE scan_for_polymorphs

  SUBROUTINE EL_TO_ELP_L(R)  ! Copy all primitives ELEMENT into ELEMENTP
    IMPLICIT  NONE
    TYPE(layout), INTENT(INOUT):: R
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    nullify(C)
    CALL LINE_L(R,doneit)

    C=>R%START
    do   WHILE(ASSOCIATED(C))
       IF(.NOT.ASSOCIATED(C%PARENT_MAG)) CALL COPY(C%MAG,C%MAGP)
       C=>C%NEXT
    ENDDO

    CALL RING_L(R,doneit)

  END SUBROUTINE EL_TO_ELP_L

  SUBROUTINE ELP_TO_EL_L(R) ! Copy all primitives ELEMENTP into ELEMENT
    IMPLICIT  NONE
    TYPE(layout), INTENT(INOUT):: R
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    nullify(C)

    CALL LINE_L(R,doneit)
    C=>R%START
    do   WHILE(ASSOCIATED(C))

       IF(.NOT.ASSOCIATED(C%PARENT_MAG)) CALL COPY(C%MAGP,C%MAG)

       C=>C%NEXT
    ENDDO

    CALL RING_L(R,doneit)
  END SUBROUTINE ELP_TO_EL_L


END  MODULE        S_FAMILY
