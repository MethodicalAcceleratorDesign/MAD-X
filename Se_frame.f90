!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
module S_FRAME
  use S_euclidean
  IMPLICIT NONE
  PRIVATE ZERO_CHART,COPY_CHART,COPY_CHART1,GEO_ROTA,GEO_ROTB,GEO_ROTV
  PRIVATE COPY_PATCH,COPY_PATCH1,ZERO_PATCH,FIND_PATCH_b
  Private ROTATE_E,alloc_f,dealloc_f,equal_f
  LOGICAL(lp),TARGET :: with_external_frame=.true.

  INTERFACE assignment (=)
     MODULE PROCEDURE ZERO_CHART
     MODULE PROCEDURE ZERO_PATCH
  end  INTERFACE

  INTERFACE EQUAL
     MODULE PROCEDURE COPY_CHART
     MODULE PROCEDURE COPY_PATCH
  end  INTERFACE
  INTERFACE copy
     MODULE PROCEDURE COPY_CHART1
     MODULE PROCEDURE COPY_PATCH1
  end  INTERFACE

  INTERFACE GEO_ROT
     MODULE PROCEDURE GEO_ROTA
     MODULE PROCEDURE GEO_ROTB
     MODULE PROCEDURE GEO_ROTV
  end  INTERFACE

  INTERFACE FIND_PATCH
     MODULE PROCEDURE FIND_PATCH_b
  end  INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE alloc_f
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE dealloc_f
  END INTERFACE

  INTERFACE assignment (=)
     MODULE PROCEDURE equal_f
  end  INTERFACE


  type magnet_frame
     real(dp), POINTER,dimension(:,:)::   ENT
     real(dp), POINTER,dimension(:)  ::   A
     real(dp), POINTER,dimension(:,:)::   EXI
     real(dp), POINTER,dimension(:)  ::   B
     real(dp), POINTER,dimension(:,:)::   MID
     real(dp), POINTER,dimension(:)  ::   O
  end type magnet_frame

  TYPE PATCH
     LOGICAL(lp), POINTER:: PATCH
     LOGICAL(lp), POINTER:: ENERGY
     LOGICAL(lp), POINTER:: TIME
     real(dp), POINTER:: A_T,B_T
     real(dp),dimension(:), POINTER:: A_D,B_D,A_ANG,B_ANG
  END TYPE PATCH

  TYPE CHART
     type(magnet_frame), pointer :: f
     real(dp), POINTER:: A_XY
     real(dp), POINTER:: L
     real(dp), POINTER:: ALPHA
     !  FIBRE MISALIGNMENTS   => Now always true!
     real(dp),dimension(:),  POINTER::   D_IN,ANG_IN
     real(dp),dimension(:),  POINTER::   D_OUT,ANG_OUT
  END TYPE CHART

CONTAINS

  SUBROUTINE  NULL_f(p)
    implicit none
    type (MAGNET_frame), pointer:: P

    nullify(P%O);nullify(P%MID);
    nullify(P%A);nullify(P%ENT);
    nullify(P%B);nullify(P%EXI);
  end subroutine NULL_f

  SUBROUTINE  alloc_f(p)
    implicit none
    type (MAGNET_frame), pointer:: P
    INTEGER I,J
    nullify(p)
    allocate(p)
    CALL NULL_f(p)
    ALLOCATE(P%O(3));ALLOCATE(P%MID(3,3));
    ALLOCATE(P%A(3));ALLOCATE(P%ENT(3,3));
    ALLOCATE(P%B(3));ALLOCATE(P%EXI(3,3));
    DO I=1,3
       P%O(i)=zero;P%A(i)=zero;P%B(i)=zero;
       DO J=1,3
          P%MID(i,j)=zero;P%ENT(i,j)=zero;P%EXI(i,j)=zero;
       ENDDO
       P%MID(i,i)=one;P%ENT(i,i)=one;P%EXI(i,i)=one;
    ENDDO
  end subroutine alloc_f

  SUBROUTINE  dealloc_f(p)
    implicit none
    type (MAGNET_frame), pointer:: P

    if(associated(p%a)) then
       DEALLOCATE(P%A);DEALLOCATE(P%ENT);
       DEALLOCATE(P%O);DEALLOCATE(P%MID);
       DEALLOCATE(P%B);DEALLOCATE(P%EXI);
    endif
  end SUBROUTINE  dealloc_f

  SUBROUTINE  equal_f(elp,el)
    implicit none
    type (MAGNET_frame),INTENT(inOUT)::elP
    type (MAGNET_frame),INTENT(IN)::el
    INTEGER I,J

    if(associated(ELP%O).and.associated(EL%O)) then
       DO I=1,3
          ELP%O(I)=EL%O(I)
          ELP%A(I)=EL%A(I)
          ELP%B(I)=EL%B(I)
          DO J=1,3
             ELP%MID(I,J)=EL%MID(I,J);ELP%ENT(I,J)=EL%ENT(I,J);ELP%EXI(I,J)=EL%EXI(I,J);
          ENDDO
       ENDDO
    endif
  end SUBROUTINE  equal_f

  SUBROUTINE COPY_PATCH(B,A)  ! B<-A
    IMPLICIT  NONE
    TYPE(PATCH), INTENT(inOUT):: B
    TYPE(PATCH), INTENT(IN):: A
    INTEGER I

    if(associated(b%PATCH).and.associated(a%PATCH)) then
       B%PATCH=A%PATCH
       B%ENERGY=A%ENERGY
       B%TIME=A%TIME

       B%A_T=A%A_T
       B%B_T=A%B_T

       DO I=1,3
          B%A_D(I)=A%A_D(I)
          B%B_D(I)=A%B_D(I)
          B%A_ANG(I)=A%A_ANG(I)
          B%B_ANG(I)=A%B_ANG(I)
       ENDDO
    endif

  END SUBROUTINE COPY_PATCH

  SUBROUTINE COPY_CHART(B,A)  ! B<-A
    IMPLICIT  NONE
    TYPE(CHART), INTENT(inOUT):: B
    TYPE(CHART), INTENT(IN):: A
    INTEGER I

    if(associated(b%A_XY).and.associated(a%A_XY)) then
       B%A_XY=A%A_XY
       B%L=A%L
       B%ALPHA=A%ALPHA

       DO I=1,3
          B%D_in(I)=A%D_in(I)
          B%ang_in(I)=A%ang_in(I)
          B%D_out(I)=A%D_out(I)
          B%ang_out(I)=A%ang_out(I)
       ENDDO
    endif
    if(associated(a%f)) then
       b%f=a%f
    endif


  END SUBROUTINE COPY_CHART

  SUBROUTINE COPY_CHART1(a,b)   ! a->b
    IMPLICIT  NONE
    TYPE(CHART), INTENT(inOUT):: B
    TYPE(CHART), INTENT(IN):: A

    call equal(b,a)

  END SUBROUTINE COPY_CHART1

  SUBROUTINE COPY_PATCH1(a,b)  ! a->b
    IMPLICIT  NONE
    TYPE(PATCH), INTENT(inOUT):: B
    TYPE(PATCH), INTENT(IN):: A

    call equal(b,a)

  END SUBROUTINE COPY_PATCH1




  SUBROUTINE ZERO_PATCH(F,R)   !R=0 nullifies and allocates ; R=-1 deallocates
    IMPLICIT  NONE
    TYPE(PATCH), INTENT(INOUT):: F
    INTEGER, INTENT(IN):: R
    INTEGER I

    IF(R==0.or.R==1) THEN
       NULLIFY(F%A_T,F%B_T,F%A_D,  F%B_D,   F%A_ANG,F%B_ANG)
       NULLIFY(F%TIME,F%ENERGY,F%PATCH)
       ALLOCATE(F%A_D(3),F%B_D(3),F%A_ANG(3),F%B_ANG(3))
       ALLOCATE(F%A_T,F%B_T)
       ALLOCATE(F%TIME,F%ENERGY,F%PATCH)
       F%A_T=zero
       F%B_T=zero
       DO I=1,3
          F%A_D(I)=zero
          F%B_D(I)=zero
          F%A_ANG(I)=zero
          F%B_ANG(I)=zero
       ENDDO
       f%patch=.false.
       f%ENERGY=.false.
       f%TIME=.false.
    ELSEIF(R==-1) THEN
       DEALLOCATE(F%A_D,F%B_D,F%A_ANG,F%B_ANG)
       DEALLOCATE(F%A_T,F%B_T)
       DEALLOCATE(F%TIME,F%ENERGY,F%PATCH)
       nullify(F%A_D,F%B_D,F%A_ANG,F%B_ANG)
       nullify(F%A_T,F%B_T)
       nullify(F%TIME,F%ENERGY,F%PATCH)
    ELSE
       w_p=1
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(a5,1x,i4,a30)') " R = ",R ," NOT DEFINED IN ZERO_CHART (1)"
       CALL WRITE_E(-1)
    ENDIF
  END SUBROUTINE ZERO_PATCH

  SUBROUTINE ZERO_CHART(F,R)   !R=0 nullifies and allocates ; R=-1 deallocates
    IMPLICIT  NONE
    TYPE(CHART), INTENT(INOUT):: F
    INTEGER, INTENT(IN):: R
    INTEGER I

    IF(R==0.or.R==1) THEN
       nullify(f%f)
       NULLIFY(    F%A_XY, F%L,    F%ALPHA)
       NULLIFY( F%d_in, F%ang_in,F%d_out,F%ang_out)
       if(with_external_frame) then
          call alloc(f%f)
       endif
       ALLOCATE(F%d_in(3),F%ang_in(3),F%d_out(3),F%ang_out(3))
       ALLOCATE(F%A_XY,F%L,F%ALPHA)
       IF(R/=1) then
          F%L=zero
          F%ALPHA=zero
       endif
       DO I=1,3
          F%Ang_in(I)=zero
          F%d_in(I)=zero
          F%Ang_out(I)=zero
          F%d_out(I)=zero
          IF(R==1.and.associated(f%f)) THEN
             F%f%ENT(I,I)=one
             F%f%EXI(I,I)=one
             F%f%MID(I,I)=one
          ENDIF
       ENDDO
       F%A_XY=zero
    ELSEIF(R==-1) THEN
       DEALLOCATE(F%d_in,F%ang_in,F%d_out,F%ang_out)
       DEALLOCATE(F%A_XY,F%L,F%ALPHA)
       if(associated(f%f)) then
          call kill(f%f)
          deallocate(f%f);
       endif
       NULLIFY(F%d_in,F%ang_in,F%d_out,F%ang_out)
       NULLIFY(F%A_XY,F%L,F%ALPHA,f%f)
    ELSE
       w_p=1
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       write(w_p%c(1),'(a5,1x,i4,a30)') " R = ",R ," NOT DEFINED IN ZERO_CHART (2)"
       CALL WRITE_E(-1)
    ENDIF

  END SUBROUTINE ZERO_CHART


  SUBROUTINE COPY_VECTOR(R,F)   ! R%ENT,R%MID,R%EXI -> F%ENT,F%MID,F%EXI
    IMPLICIT  NONE
    TYPE(CHART), INTENT(inOUT):: F
    TYPE(CHART), INTENT(IN):: R
    INTEGER I,J
    if(associated(r%f).and.associated(f%f)) then
       DO I=1,3
          DO J=1,3
             F%f%MID(I,J)=R%f%MID(I,J)
             F%f%EXI(I,J)=R%f%EXI(I,J)
             F%f%ENT(I,J)=R%f%ENT(I,J)
          ENDDO
       ENDDO
    endif
  END SUBROUTINE COPY_VECTOR




  SUBROUTINE GEO_TRA(A,ENT,D,I) ! Adds/subtracts D to A where D is expressed in the ENT frame
    implicit none
    real(dp), INTENT(INOUT):: A(3)
    real(dp), INTENT(IN):: ENT(3,3)
    real(dp), INTENT(IN):: D(3)
    INTEGER, INTENT(IN):: I
    real(dp) B(3)
    INTEGER J,K
    B=zero
    DO J=1,3
       DO K=1,3
          B(K)=D(J)*ENT(J,K)+B(K)
       ENDDO
    ENDDO
    A=A+I*B
  END SUBROUTINE GEO_TRA

  SUBROUTINE GEO_AVE(ENT,a,EXI,b,MID,o) ! Finds MID and O
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3),EXI(3,3),MID(3,3)
    real(dp), INTENT(INOUT):: a(3),b(3),o(3)
    INTEGER  I
    real(dp) N1,N3

    N1=zero;N3=zero;
    DO I=1,3
       MID(1,I)=ENT(1,I)+EXI(1,I)
       N1=N1+MID(1,I)**2
    ENDDO
    DO I=1,3
       MID(3,I)=ENT(3,I)+EXI(3,I)
       N3=N3+MID(3,I)**2
    ENDDO

    N1=SQRT(N1);N3=SQRT(N3);

    DO I=1,3
       MID(1,I)=MID(1,I)/N1 ; MID(3,I)=MID(3,I)/N3 ;
    ENDDO
    DO I=1,3
       o(I)=(a(I)+b(I))*half
    ENDDO

  END SUBROUTINE GEO_AVE


  SUBROUTINE GEO_ROTA(ENT,A,I) ! Rotates frame ENT by A(3) in the PTC or reverse PTC order
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3)
    real(dp), INTENT(IN):: A(3)
    INTEGER, INTENT(IN):: I

    IF(I==1) THEN
       CALL GEO_ROT(ENT,ENT,A(1),A(2),A(3))
    ELSE
       CALL GEO_ROT(ENT,ENT,zero,zero,-A(3))
       CALL GEO_ROT(ENT,ENT,zero,-A(2),zero)
       CALL GEO_ROT(ENT,ENT,-A(1),zero,zero)
    ENDIF

  END SUBROUTINE GEO_ROTA

  SUBROUTINE GEO_ROTB(ENT,EXI,A_XY,A_XZ,A_YZ) ! Rotates ENT into EXI
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3),EXI(3,3)
    real(dp), INTENT(IN):: A_XY,A_XZ,A_YZ
    real(dp) TEMP(3,3)
    INTEGER I
    ! Definition
    !  ENT(1,i) is the local x-vector
    !  ENT(2,i) is the local y-vector
    !  ENT(3,i) is the local z-vector
    ! Standardized to correspound to the standard PTC order
    ! Original order now changed was xy,xz,yz

    DO I=1,3
       TEMP(2,I)=  COS(A_YZ)*ENT(2,I)+SIN(A_YZ)*ENT(3,I)
       TEMP(3,I)= -SIN(A_YZ)*ENT(2,I)+COS(A_YZ)*ENT(3,I)
       TEMP(1,I)=  ENT(1,I)
    ENDDO

    DO I=1,3
       EXI(1,I)=  COS(A_XZ)*TEMP(1,I)+SIN(A_XZ)*TEMP(3,I)
       EXI(3,I)= -SIN(A_XZ)*TEMP(1,I)+COS(A_XZ)*TEMP(3,I)
       EXI(2,I)=  TEMP(2,I)
    ENDDO

    DO I=1,3
       TEMP(1,I)=  COS(A_XY)*EXI(1,I)+SIN(A_XY)*EXI(2,I)
       TEMP(2,I)= -SIN(A_XY)*EXI(1,I)+COS(A_XY)*EXI(2,I)
       TEMP(3,I)=  EXI(3,I)
    ENDDO


    DO I=1,3
       EXI(1,I)= TEMP(1,I)
       EXI(2,I)= TEMP(2,I)
       EXI(3,I)= TEMP(3,I)
    ENDDO


  END SUBROUTINE GEO_ROTB

  SUBROUTINE GEO_ROTV(V,A) ! Rotates V by angles A(3) in PTC order
    implicit none
    real(dp), INTENT(INOUT):: V(3)
    real(dp), INTENT(IN):: A(3)
    real(dp)  A_XY,A_XZ,A_YZ
    real(dp) W(3)
    INTEGER I
    ! Definition
    !  ENT(1,i) is the local x-vector
    !  ENT(2,i) is the local y-vector
    !  ENT(3,i) is the local z-vector
    ! Standardized to correspound to the standard PTC order
    ! Original order now changed was xy,xz,yz
    A_YZ=A(1);A_XZ=A(2);A_XY=A(3);

    DO I=1,3
       W(2)=  COS(A_YZ)*V(2)-SIN(A_YZ)*V(3)
       W(3)= +SIN(A_YZ)*V(2)+COS(A_YZ)*V(3)
       W(1)=  V(1)
    ENDDO

    DO I=1,3
       V(1)=  COS(A_XZ)*W(1)-SIN(A_XZ)*W(3)
       V(3)= +SIN(A_XZ)*W(1)+COS(A_XZ)*W(3)
       V(2)=  W(2)
    ENDDO

    DO I=1,3
       W(1)=  COS(A_XY)*V(1)-SIN(A_XY)*V(2)
       W(2)= +SIN(A_XY)*V(1)+COS(A_XY)*V(2)
       W(3)=  V(3)
    ENDDO


    V=W

  END SUBROUTINE GEO_ROTV



  SUBROUTINE ROTATE_V(E,F,dirver)  ! Rotates exit into entrance or vice versa (for survey)
    IMPLICIT NONE
    TYPE(CHART), INTENT(INOUT):: E
    TYPE(CHART), INTENT(INOUT):: F
    real(dp) V(3),A,ent(3,3)
    INTEGER I,dirver

    if(associated(E%f).and.associated(f%f)) then

       if(dirver==1) then
          !ROTATE EXIT OF LAST ONE IN X-Y
          !WRITE(16,*) "ENTERING ROTATE_V "
          !WRITE(16,*) " F%A_XY ",F%A_XY

          F%f%ENT=E%f%EXI
          CALL GEO_ROT(E%f%EXI,ENT,F%A_XY,zero,zero)



          !ROTATE ENTRANCE INTO EXIT

          CALL GEO_ROT(ENT,F%f%EXI,zero,F%Alpha,zero)


          !REMOVE 90% AT  EXIT

          CALL GEO_ROT(F%f%EXI,F%f%EXI,-F%A_XY,zero,zero)

          !WRITE(16,*) " F%B_XY ",F%B_XY


          !  SHOULD STILL BE OK!!!!
          A=F%Alpha/two
          DO I=1,3
             V(1)=   COS(A)*ENT(3,1)-SIN(A)*ENT(1,1)
             V(2)=   COS(A)*ENT(3,2)-SIN(A)*ENT(1,2)
             V(3)=   COS(A)*ENT(3,3)-SIN(A)*ENT(1,3)
          ENDDO

          DO I=1,3
             F%f%B(I)=F%L*V(I)+F%f%A(I)
          ENDDO
       else   !  reverse propagator
          !ROTATE ent OF LAST ONE IN X-Y
          F%f%exi=E%f%ent    ! ent is really exi here
          CALL GEO_ROT(E%f%ent,ent,F%A_XY,zero,zero)



          !ROTATE ENTRANCE INTO EXIT

          CALL GEO_ROT(ent,F%f%ent,zero,-F%Alpha,zero)


          !REMOVE 90% AT  EXIT

          CALL GEO_ROT(F%f%ent,F%f%ent,-F%A_XY,zero,zero)



          !  SHOULD STILL BE OK!!!!
          A=-F%Alpha/two
          !A=F%A_XZ/two
          DO I=1,3    ! ent is exi really!
             V(1)=   COS(A)*ent(3,1)-SIN(A)*ent(1,1)
             V(2)=   COS(A)*ent(3,2)-SIN(A)*ent(1,2)
             V(3)=   COS(A)*ent(3,3)-SIN(A)*ent(1,3)
          ENDDO

          DO I=1,3
             F%f%A(I)=-F%L*V(I)+F%f%B(I)
          ENDDO
       endif
    endif
  END SUBROUTINE ROTATE_V

  SUBROUTINE ROTATE_C(E,F,DIRVER) ! PUTS E AT THE END OF F FOR SURVEY PURPOSES. E AND F MUST HAVE SAME DIRVER
    IMPLICIT NONE
    TYPE(CHART), INTENT(INOUT):: E
    TYPE(CHART), INTENT(INOUT):: F
    real(dp) B(3),ROT(3,3)
    INTEGER,INTENT(IN):: dirver
    INTEGER I,J,K
    if(associated(E%f).and.associated(f%f)) then
       IF(DIRVER==1) THEN
          ROT=zero
          DO I=1,3 ;DO J=1,3 ;DO K=1,3 ;
             ROT(I,J)=ROT(I,J)+E%f%ENT(K,I)*F%f%EXI(K,J)
          ENDDO;ENDDO ;ENDDO;
          B=F%f%B
          CALL ROTATE_E(E%f%A,E%f%ENT,E%f%O,E%f%MID,E%f%B,E%f%EXI,B,ROT)
       ELSE
          ROT=zero
          DO I=1,3 ;DO J=1,3 ;DO K=1,3 ;
             ROT(I,J)=ROT(I,J)+E%f%EXI(K,I)*F%f%ENT(K,J)
          ENDDO;ENDDO ;ENDDO;
          B=F%f%A
          CALL ROTATE_E(E%f%B,E%f%EXI,E%f%O,E%f%MID,E%f%A,E%f%ENT,B,ROT)
       ENDIF
    endif
  END SUBROUTINE ROTATE_C


  SUBROUTINE ROTATE_E(A,ENT,O,MID,B,EXI,B_LAT,ROT) ! This routine is used in ROTATE_E

    IMPLICIT NONE
    real(dp), INTENT(INOUT):: ENT(3,3),MID(3,3),EXI(3,3),ROT(3,3)
    real(dp), INTENT(INOUT):: A(3),O(3),B(3),B_LAT(3)
    real(dp) TEMP(3,3),T(3)
    INTEGER I,J,K


    TEMP=zero
    DO I=1,3 ;DO J=1,3 ;DO K=1,3 ;
       TEMP(I,J)=TEMP(I,J)+ENT(I,K)*ROT(K,J)
    ENDDO;ENDDO ;ENDDO;
    ENT=TEMP
    TEMP=zero
    DO I=1,3 ;DO J=1,3 ;DO K=1,3 ;
       TEMP(I,J)=TEMP(I,J)+MID(I,K)*ROT(K,J)
    ENDDO;ENDDO ;ENDDO;
    MID=TEMP
    TEMP=zero
    DO I=1,3 ;DO J=1,3 ;DO K=1,3 ;
       TEMP(I,J)=TEMP(I,J)+EXI(I,K)*ROT(K,J)
    ENDDO;ENDDO ;ENDDO;
    EXI=TEMP

    O=O-A
    B=B-A
    A=B_LAT
    T=zero
    DO I=1,3 ; DO J=1,3;
       T(I)=ROT(J,I)*O(J)+T(I)
    ENDDO; ENDDO;
    O=T+A
    T=zero
    DO I=1,3 ; DO J=1,3;
       T(I)=ROT(J,I)*B(J)+T(I)
    ENDDO; ENDDO;
    B=T+A


  END SUBROUTINE ROTATE_E


  SUBROUTINE FIND_PATCH_b(A,ENT,B,EXI,D,ANG) ! Finds patch between ENT and EXI : Interfaced later for fibres
    USE rotation_mis
    IMPLICIT NONE
    real(dp), INTENT(INOUT):: ENT(3,3),EXI(3,3)
    real(dp), INTENT(INOUT):: A(3),B(3),D(3),ANG(3)
    !   THIS ROUTINE FINDS THE PATCH CONNECT A TO B IN THE FORWARD
    !   PROPAGATION  CASE
    TYPE(MATRIX_PTC) ROT
    INTEGER I,J,K


    ROT=0
    d=zero
    ANG=zero
    DO I=1,3 ;DO J=1,3 ;DO K=1,3 ;
       ROT%R(I,J)=ROT%R(I,J)+ENT(K,I)*EXI(K,J)   ! Vector rotation not component rotation
    ENDDO;ENDDO ;ENDDO;
    CALL factorize_m2(ROT,ANG)
    ANG(2)=-ANG(2)                             ! dynamical PTC convention
    rot%t=B-A
    ! express rot%t in the frame EXI : PTC convention
    do i=1,3
       do j=1,3
          d(i)=exi(i,j)*rot%t(j)+d(i)
       enddo
    enddo

  END SUBROUTINE FIND_PATCH_b


end module S_FRAME
