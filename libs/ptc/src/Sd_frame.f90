!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module S_FRAME
  use S_euclidean
  IMPLICIT NONE
  public
  PRIVATE ZERO_CHART,COPY_CHART,COPY_CHART1   !,GEO_ROTA,GEO_ROTB
  PRIVATE COPY_PATCH,COPY_PATCH1,ZERO_PATCH,FIND_PATCH_b,FIND_PATCH_bmad0
  Private alloc_f,dealloc_f,equal_f
  private make_rot_x,make_rot_y,make_rot_z   !,GEO_ROTAB_no_vec,GEO_ROTA_no_vec
  REAL(DP), public :: GLOBAL_FRAME(3,3)= RESHAPE((/1,0,0  ,0,1,0  ,0,0,1/),(/3,3/))
  REAL(DP), public :: GLOBAL_origin(3)= (/0,0,0/)
  integer :: ccc=0
  !  include "a_def_frame_patch_chart.inc"   ! sept 2007


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
     MODULE PROCEDURE GEO_ROTA      !(ENT,V,A,I)
     MODULE PROCEDURE GEO_ROTA_no_vec  !(ENT,A,I)
     MODULE PROCEDURE GEO_ROTAB_no_vec   ! (ENT,Exi,A_X2,A_X1,A_xy)
     MODULE PROCEDURE GEO_ROTB      ! (ENT,EXI,A,B,A_X2,A_X1,A_XY)
  end  INTERFACE



  INTERFACE FIND_PATCH
     MODULE PROCEDURE FIND_PATCH_b
  end  INTERFACE


  INTERFACE FIND_PATCH_bmad
     MODULE PROCEDURE FIND_PATCH_bmad0
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
    nullify(p)
    allocate(p)
    CALL NULL_f(p)
    ALLOCATE(P%O(3));ALLOCATE(P%MID(3,3));
    ALLOCATE(P%A(3));ALLOCATE(P%ENT(3,3));
    ALLOCATE(P%B(3));ALLOCATE(P%EXI(3,3));
    P%O=global_origin
    P%a=global_origin
    P%b=global_origin
    P%ENT=global_FRAME
    P%MID=global_FRAME
    P%EXI=global_FRAME

  end subroutine alloc_f

  SUBROUTINE  NULL_af(p)
    implicit none
    type (affine_frame), pointer:: P
    nullify(P%ANGLE);nullify(P%D);
    nullify(P%A);nullify(P%ENT);
    nullify(P%b);nullify(P%EXI);
  end subroutine NULL_af

  SUBROUTINE  NULL_gs(gs)
    implicit none
    integer i
    type (girder_siamese) :: gs(:)
    do i=1,size(gs)
       nullify(gs(i)%mag)
    enddo
  end subroutine NULL_gs

  SUBROUTINE  alloc_af(p,girder)
    implicit none
    type (affine_frame), pointer:: P
    logical(lp), optional :: girder
    logical(lp) gird
    gird=.false.
    if(present(girder)) gird=girder

    nullify(p)
    allocate(p)
    CALL NULL_af(p)
    if(gird) then
       ALLOCATE(P%A(3));ALLOCATE(P%ENT(3,3));
       ALLOCATE(P%B(3));ALLOCATE(P%EXI(3,3));
       P%a=global_origin
       P%ENT=global_FRAME
       P%B=global_origin
       P%EXI=global_FRAME
    else
       ALLOCATE(P%ANGLE(3));ALLOCATE(P%D(3));
       P%D=0.0_dp
       P%ANGLE=0.0_dp
    endif

  end subroutine alloc_af

  SUBROUTINE  kill_af(p)
    implicit none
    type (affine_frame), pointer:: P

    if(associated(p%D)) DEALLOCATE(P%D);
    if(associated(p%ANGLE)) DEALLOCATE(P%ANGLE);
    if(associated(p%a)) DEALLOCATE(P%A);
    if(associated(p%a)) DEALLOCATE(P%A);
    if(associated(p%ENT)) DEALLOCATE(P%ENT);
    if(associated(p%B)) DEALLOCATE(P%B);
    if(associated(p%EXI)) DEALLOCATE(P%EXI);

  end subroutine kill_af

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

    if(associated(ELP%O).and.associated(EL%O)) then
       ELP%O=EL%O
       ELP%A=EL%A
       ELP%B=EL%B
       ELP%MID=EL%MID;ELP%ENT=EL%ENT;ELP%EXI=EL%EXI;
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
       B%A_L=A%A_L
       B%B_L=A%B_L
       B%b0b=A%b0b
       B%p0b=A%p0b
       B%A_X1=A%A_X1
       B%B_X1=A%B_X1
       B%A_X2=A%A_X2
       B%B_X2=A%B_X2

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

    if(associated(b%f).and.associated(a%f)) then
       !   B%A_XY=A%A_XY
       !   B%L=A%L
       !   B%ALPHA=A%ALPHA

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
    !    if(r==1) then
    !    write(6,*) "i =0 "
    !    read(5,*) i
    !    i=1/i
    !   endif

    IF(R==0.or.R==1) THEN    !
       NULLIFY(F%A_T,F%B_T,F%A_L,F%B_L,F%A_D,  F%B_D,   F%A_ANG,F%B_ANG)
       NULLIFY(F%A_X1,F%A_X2,F%B_X1,F%B_X2,F%p0b,F%B0b)

       NULLIFY(F%TIME,F%ENERGY,F%PATCH)
       ALLOCATE(F%A_D(3),F%B_D(3),F%A_ANG(3),F%B_ANG(3))
       ALLOCATE(F%A_T,F%B_T,F%A_L,F%B_L)
       ALLOCATE(F%A_X1,F%A_X2,F%B_X1,F%B_X2)
       ALLOCATE(F%TIME,F%ENERGY,F%PATCH,F%p0b,F%B0b)
       F%A_T=0.0_dp
       F%B_T=0.0_dp
       F%A_L=0.0_dp
       F%B_L=0.0_dp
       F%A_X1=1; F%A_X2=1; F%B_X1=1; F%B_X2=1;
       F%A_D=0.0_dp
       F%B_D=0.0_dp
       F%A_ANG=0.0_dp
       F%B_ANG=0.0_dp
       F%p0b=0.0_dp
       F%B0b=0.0_dp
       f%patch=0
       f%ENERGY=0
       f%TIME=0
    ELSEIF(R==-1) THEN
       DEALLOCATE(F%A_D,F%B_D,F%A_ANG,F%B_ANG,F%p0b,F%B0b)
       DEALLOCATE(F%A_T,F%B_T,F%A_L,F%B_L)
       DEALLOCATE(F%A_X1,F%A_X2,F%B_X1,F%B_X2)
       DEALLOCATE(F%TIME,F%ENERGY,F%PATCH)
       nullify(F%A_D,F%B_D,F%A_ANG,F%B_ANG,F%p0b,F%B0b)
       nullify(F%A_T,F%B_T,F%A_L,F%B_L)
       nullify(F%A_X1,F%A_X2,F%B_X1,F%B_X2)
       nullify(F%TIME,F%ENERGY,F%PATCH)
    ELSE
 
       write(6,'(a5,1x,i4,a30)') " R = ",R ," NOT DEFINED IN ZERO_CHART (1)"
       ! call !write_e(-1)
    ENDIF
  END SUBROUTINE ZERO_PATCH



  SUBROUTINE ZERO_CHART(F,R)   !R=0 nullifies and allocates ; R=-1 deallocates
    IMPLICIT  NONE
    TYPE(CHART), INTENT(INOUT):: F
    INTEGER, INTENT(IN):: R
    !    if(r==1) then
    !    write(6,*) "i =0 "
    !    read(5,*) i
    !    i=1/i
    !   endif
    IF(R==0.or.R==1) THEN          !
       nullify(f%f)
       !       NULLIFY(    F%A_XY, F%L,    F%ALPHA)
       NULLIFY( F%d_in, F%ang_in,F%d_out,F%ang_out)
       call alloc(f%f)
       ALLOCATE(F%d_in(3),F%ang_in(3),F%d_out(3),F%ang_out(3))
       !      ALLOCATE(F%A_XY,F%L,F%ALPHA)
       !         F%L=zero
       !         F%ALPHA=zero
       !         F%A_XY=zero
       F%Ang_in=0.0_dp
       F%d_in=0.0_dp
       F%Ang_out=0.0_dp
       F%d_out=0.0_dp
       IF(associated(f%f)) THEN   ! R==1.and.
          F%f%ENT=global_frame
          F%f%EXI=global_frame
          F%f%MID=global_frame
          F%f%ENT=global_frame
          F%f%EXI=global_frame
          F%f%MID=global_frame
          F%f%A=GLOBAL_origin
          F%f%B=GLOBAL_origin
          F%f%O=GLOBAL_origin
       ENDIF
    ELSEIF(R==-1) THEN
       DEALLOCATE(F%d_in,F%ang_in,F%d_out,F%ang_out)
       !       DEALLOCATE(F%A_XY,F%L,F%ALPHA)
       if(associated(f%f)) then
          call kill(f%f)
          deallocate(f%f)
       endif
       NULLIFY(F%d_in,F%ang_in,F%d_out,F%ang_out,f%f)
       !       NULLIFY(F%A_XY,F%L,F%ALPHA)
    ELSEIF(R==2) THEN   ! SPECIAL TRYING TO FIX FRANK'S MEMORY IN EL_Q
       DEALLOCATE(F%d_in,F%ang_in,F%d_out,F%ang_out)
       !       DEALLOCATE(F%A_XY,F%L,F%ALPHA)
       if(associated(f%f)) then
          call kill(f%f)
          deallocate(f%f)
       endif
       NULLIFY(F%d_in,F%ang_in,F%d_out,F%ang_out,f%f)
       call alloc(f%f)
       ALLOCATE(F%d_in(3),F%ang_in(3),F%d_out(3),F%ang_out(3))
       !      ALLOCATE(F%A_XY,F%L,F%ALPHA)
       !         F%L=zero
       !         F%ALPHA=zero
       !         F%A_XY=zero
       F%Ang_in=0.0_dp
       F%d_in=0.0_dp
       F%Ang_out=0.0_dp
       F%d_out=0.0_dp
       IF(associated(f%f)) THEN   ! R==1.and.
          F%f%ENT=global_frame
          F%f%EXI=global_frame
          F%f%MID=global_frame
          F%f%ENT=global_frame
          F%f%EXI=global_frame
          F%f%MID=global_frame
          F%f%A=GLOBAL_origin
          F%f%B=GLOBAL_origin
          F%f%O=GLOBAL_origin
       ENDIF
    ELSE
 
       write(6,'(a5,1x,i4,a30)') " R = ",R ," NOT DEFINED IN ZERO_CHART (2)"
       ! call !write_e(-1)
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
  ! Geometry

  SUBROUTINE change_basis(A,ENT,B,EXI)
    implicit none
    real(dp), INTENT(IN):: A(3),ENT(3,3)
    real(dp), INTENT(INOUT):: EXI(3,3),B(3)
    real(dp) T(3)
    INTEGER L,K,N
    T=A
    B=0.0_dp
    DO L=1,3
       DO N=1,3
          DO K=1,3
             B(N)=T(L)*ENT(L,K)*EXI(N,K)+B(N)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE change_basis



  SUBROUTINE GEO_TRA(A,ENT,D,I)
    ! Adds/subtracts D to A where D is expressed in the ENT frame
    implicit none
    real(dp), INTENT(INOUT):: A(3)
    real(dp), INTENT(IN):: ENT(3,3)
    real(dp), INTENT(IN):: D(3)
    INTEGER, INTENT(IN):: I
    real(dp) B(3)
    INTEGER J,K
    B=0.0_dp
    DO J=1,3
       DO K=1,3
          B(K)=D(J)*ENT(J,K)+B(K)
       ENDDO
    ENDDO
    A=A+I*B
  END SUBROUTINE GEO_TRA


  SUBROUTINE GEO_ROTA_no_vec(ENT,ANG,I,basis)
    ! Rotates frame ENT by A(3) in the PTC or reverse
    !PTC order using global frame for angle definition
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3)
    real(dp), INTENT(IN):: ANG(3)
    INTEGER, INTENT(IN):: I
    real(dp) V(3)
    real(dp), optional, INTENT(IN):: basis(3,3)

    v=0.0_dp
    CALL GEO_ROT(ENT,V,ANG,i,basis)

  END SUBROUTINE GEO_ROTA_no_vec

  SUBROUTINE GEO_ROTA(ENT,A,ANG,I,basis)
    ! Rotates frame ENT by Ang(3) in the PTC or reverse PTC order
    !using global frame for angle definition
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3),A(3)
    real(dp), INTENT(IN):: ANG(3)
    INTEGER, INTENT(IN):: I
    real(dp), optional, INTENT(IN):: basis(3,3)
    REAL(DP) AT(3),AA(3),ent0(3,3)
    INTEGER J

    ! PTC order X first, than y and finally Z , here defined on global variables
    aa=a
    ent0=ent
    IF(I==1) THEN
       CALL GEO_ROT(ENT0,ENT,Aa,A,ANG,basis)
    ELSE
       DO J=1,3
          AT=0.0_dp;AT(J)=-ANG(J);
          CALL GEO_ROT(ENT0,ENT,Aa,A,AT,basis)
       ENDDO
    ENDIF

  END SUBROUTINE GEO_ROTA



  SUBROUTINE GEO_ROTAB_no_vec(ENT,Exi,Ang,basis) ! Used in survey stuff A_X2,A_X1,A_xy,
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3),exi(3,3)
    real(dp), INTENT(IN):: ANG(3)
    real(dp) v(3),vv(3)
    real(dp), optional, INTENT(IN):: basis(3,3)


    v=0.0_dp
    vv=0.0_dp

    CALL GEO_ROT(ENT,Exi,VV,V,ANG,basis)    !A_X2,A_X1,A_xy,

  END SUBROUTINE GEO_ROTAB_no_vec

  SUBROUTINE  ROTATE_FRAME(R,OMEGA,ANG,ORDER,BASIS)
    ! ROTATES A FRAME R AROUND OMEGA BY ANG(3)  IN STANDARD PTC ORDER
    ! INVERSE => ORDER=-1   USING GLOBAL FRAME
    IMPLICIT NONE
    TYPE (MAGNET_FRAME),TARGET,INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: OMEGA(3),ANG(3)
    TYPE(MAGNET_FRAME), POINTER::P
    REAL(DP) D(3),OMEGAT(3)
    INTEGER IORDER
    INTEGER, OPTIONAL :: ORDER
    REAL(DP), OPTIONAL, INTENT(IN):: BASIS(3,3)
    REAL(DP)  BASIST(3,3)
    OMEGAT=OMEGA
    P=>R
    IORDER=1
    IF(PRESENT(ORDER)) IORDER=ORDER

    BASIST=GLOBAL_FRAME            ! NECESSARY SINCE BASIS CAN CHANGE DURING THE CALCULATION ASSUMING A POINTER IS PASSED
    IF(PRESENT(BASIS)) BASIST=BASIS

    ! THEY WILL ROTATE MORE THAN ONCE
    D=P%A-OMEGAT
    CALL GEO_ROT(P%ENT,D,ANG,IORDER,BASIST)
    P%A=OMEGAT+D

    D=P%O-OMEGAT
    CALL GEO_ROT(P%MID,D,ANG,IORDER,BASIST)
    P%O=OMEGAT+D

    D=P%B-OMEGAT
    CALL GEO_ROT(P%EXI,D,ANG,IORDER,BASIST)
    P%B=OMEGAT+D


  END SUBROUTINE ROTATE_FRAME


  SUBROUTINE  TRANSLATE_FRAME(R,D,ORDER,BASIS) ! TRANSLATES A FRAME
    IMPLICIT NONE
    TYPE (MAGNET_FRAME),TARGET, INTENT(INOUT):: R
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    INTEGER IORDER
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    REAL(DP) DD(3)
    TYPE(MAGNET_FRAME), POINTER :: P
    ! THIS ROUTINE TRANSLATE THE ENTIRE LINE BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE D

    P=>R

    IORDER=1
    DD=D
    IF(PRESENT(ORDER)) IORDER=ORDER
    IF(PRESENT(BASIS)) THEN
       CALL CHANGE_BASIS(D,BASIS,DD,GLOBAL_FRAME)
    ENDIF



    P%A=P%A+IORDER*DD
    P%B=P%B+IORDER*DD
    P%O=P%O+IORDER*DD



  END SUBROUTINE TRANSLATE_FRAME

  SUBROUTINE  TRANSLATE_point(R,D,ORDER,BASIS) ! TRANSLATES A FRAME
    IMPLICIT NONE
    REAL(DP),TARGET, INTENT(INOUT):: R(3)
    REAL(DP),INTENT(IN):: D(3)
    REAL(DP), OPTIONAL :: BASIS(3,3)
    INTEGER IORDER
    INTEGER, OPTIONAL, INTENT(IN) :: ORDER
    REAL(DP) DD(3)

    ! THIS ROUTINE TRANSLATE r(3)BY A(3) IN STANDARD ORDER USING THE
    ! GLOBAL FRAME TO DEFINE D

    IORDER=1
    DD=D
    IF(PRESENT(ORDER)) IORDER=ORDER
    IF(PRESENT(BASIS)) THEN
       CALL CHANGE_BASIS(D,BASIS,DD,GLOBAL_FRAME)
    ENDIF

    r=r+IORDER*DD
 
  END SUBROUTINE TRANSLATE_point


  SUBROUTINE make_rot_x(r,a)
    implicit none
    real(dp), INTENT(INOUT):: r(3,3)
    real(dp), INTENT(IN):: a
    r=0.0_dp

    r(1,1)= 1.0_dp

    r(2,2)=cos(a)
    r(3,3)=cos(a)
    r(2,3)=sin(a)   ! defined acting on vectors, so y goes into z
    r(3,2)=-sin(a)

  end SUBROUTINE make_rot_x

  SUBROUTINE make_rot_y(r,a)
    implicit none
    real(dp), INTENT(INOUT):: r(3,3)
    real(dp), INTENT(IN):: a
    r=0.0_dp

    r(2,2)= 1.0_dp

    r(1,1)=cos(a)
    r(3,3)=cos(a)
    r(1,3)=sin(a)   ! defined acting on vectors, so x goes into z
    r(3,1)=-sin(a)

  end SUBROUTINE make_rot_y

  SUBROUTINE make_rot_z(r,a)
    implicit none
    real(dp), INTENT(INOUT):: r(3,3)
    real(dp), INTENT(IN):: a
    r=0.0_dp

    r(3,3)= 1.0_dp

    r(1,1)=cos(a)
    r(2,2)=cos(a)
    r(1,2)=sin(a)   ! defined acting on vectors, so x goes into y
    r(2,1)=-sin(a)

  end SUBROUTINE make_rot_z


  SUBROUTINE GEO_ROTB(ENT,EXI,A,B,ANG,BASIS)
    ! Rotates ENT into EXI using global frame otherwise uses basis
    implicit none
    real(dp), INTENT(INOUT):: ENT(3,3),EXI(3,3),A(3),B(3)
    real(dp), INTENT(IN):: ANG(3)
    real(dp) TEMP(3,3),exii(3,3),enti(3,3),T(3),R(3,3),basist(3,3)
    real(dp), optional, INTENT(IN):: basis(3,3)
    INTEGER I,j,k
    ! Definition
    !  ENT(1,i) is the local x-vector
    !  ENT(2,i) is the local y-vector
    !  ENT(3,i) is the local z-vector
    ! Standardized to correspound to the standard PTC order
    ! Original order now changed was xy,xz,yz
    ! TRANSPOSED OF GEO_ROTB

    exii(:,:)=0.0_dp
    enti(:,:)=0.0_dp
    exii(:,:)=0.0_dp
    temp(:,:)=0.0_dp


    call make_rot_z(enti,ANG(3))
    call make_rot_y(r,ANG(2))
    if(present(basis)) then
       basist=basis
    else
       basist=0.0_dp
       do i=1,3
          basist(i,i)=1.0_dp
       enddo
    endif

    DO I=1,3
       DO j=1,3
          DO k=1,3
             TEMP(i,k)=  basist(j,i)*enti(j,k)+TEMP(i,k)   ! basis^-1 * r(a(3))
          ENDDO
       ENDDO
    ENDDO
    enti=temp
    temp(:,:)=0.0_dp


    DO I=1,3
       DO j=1,3
          DO k=1,3
             TEMP(i,k)=  enti(i,j)*r(j,k)+TEMP(i,k)     ! basis^-1 * r(a(3)) *r(a(2))
          ENDDO
       ENDDO
    ENDDO

    call make_rot_x(r,ANG(1))

    DO I=1,3
       DO j=1,3
          DO k=1,3
             exii(i,k)=  TEMP(i,j)*r(j,k)+exii(i,k)   ! basis^-1 * r(a(3)) *r(a(2)) *r(a(1))
          ENDDO
       ENDDO
    ENDDO


    temp=0.0_dp
    DO I=1,3
       DO j=1,3
          DO k=1,3
             temp(i,k)=  exii(i,j)*basist(j,k)+temp(i,k)  ! basis^-1 * r(a(3)) *r(a(2)) *r(a(1)) *basis
          ENDDO
       ENDDO
    ENDDO
    exii=temp


    T(:)=0.0_dp

    do i=1,3; do j=1,3;
       T(i)=A(j)*exii(j,i)+T(i)
    enddo;enddo;

    B=T

    temp(:,:)=0.0_dp

    do i=1,3; do j=1,3;do k=1,3
       temp(i,k)=ent(i,j)*exii(j,k)+temp(i,k)
    enddo;enddo;enddo;

    exi=temp

    call check_frame(exi,t)

  END SUBROUTINE GEO_ROTB

  subroutine check_frame(ent,n)
    implicit none
    real(dp) ent(3,3),n(3),s,ss
    integer i,j

    n=0.0_dp
    do i=1,3
       do j=1,3
          n(i)=ent(i,j)**2+n(i)

       enddo

       n(i)=sqrt(n(i))
    enddo

    do i=1,3
       do j=1,3
          ent(i,j)=ent(i,j)/n(i)
       enddo
    enddo

    call COMPUTE_SCALAR(ent,1,ent,2,S)

    do j=1,3
       ent(2,j)=ent(2,j)-s*ent(1,j)
    enddo

    call COMPUTE_SCALAR(ent,2,ent,2,S)

    do j=1,3
       ent(2,j)=ent(2,j)/s
    enddo


    call COMPUTE_SCALAR(ent,1,ent,3,S)
    call COMPUTE_SCALAR(ent,2,ent,3,SS)

    do j=1,3
       ent(3,j)=ent(3,j)-s*ent(1,j)-SS*ENT(2,J)
    enddo

    call COMPUTE_SCALAR(ent,3,ent,3,S)

    do j=1,3
       ent(3,j)=ent(3,j)/s
    enddo



  end subroutine check_frame

  subroutine make_normal(n,s)
    implicit none
    real(dp) n(3),s
    integer i

    s=0.0_dp
    do i=1,3
       s=s+n(i)**2
    enddo

    if(s>EPS_FITTED) then
       s=sqrt(s)
       n=n/s
    else
       s=0.0_dp
    endif
  end subroutine make_normal

   SUBROUTINE COMPUTE_ENTRANCE_ANGLE(ENTL,ENTB,A)
    ! COMPUTES PTC'S ANGLES IN FRAME OF ENTL
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: ENTL(3,3), ENTB(3,3)
    REAL(DP), INTENT(OUT) :: A(3)
    REAL(DP) T1(3,3),T10(3,3),T2(3,3),S_IJ,S_JJ
    REAL(DP) AT(3)

    T1=ENTL
    T10=ENTL
    T2=ENTL


    CALL COMPUTE_SCALAR(T1,2,ENTB,3,S_IJ)
    CALL COMPUTE_SCALAR(T1,3,ENTB,3,S_JJ)

    if(S_IJ==0.0_dp.and.S_JJ==0.0_dp) then
       A(1)=0.0_dp
    else
       A(1)=ATAN2(-S_IJ,S_JJ)
    endif

    !   A(1)=ATAN2(-S_IJ,S_JJ)
    AT=0.0_dp;AT(1)=A(1);

    CALL GEO_ROT(T10,T1,AT,T2)
    T2=T1
    T10=T1

    CALL COMPUTE_SCALAR(T1,1,ENTB,3,S_IJ)
    CALL COMPUTE_SCALAR(T1,3,ENTB,3,S_JJ)

    if(S_IJ==0.0_dp.and.S_JJ==0.0_dp) then
       A(2)=0.0_dp
    else
       A(2)=ATAN2(-S_IJ,S_JJ)
    endif

    !    A(2)=ATAN2(-S_IJ,S_JJ)
    AT=0.0_dp;AT(2)=A(2);

    CALL GEO_ROT(T10,T1,AT,T2)
    T2=T1
    T10=T1

    CALL COMPUTE_SCALAR(T1,2,ENTB,1,S_IJ)
    CALL COMPUTE_SCALAR(T1,1,ENTB,1,S_JJ)

    if(S_IJ==0.0_dp.and.S_JJ==0.0_dp) then
       A(3)=0.0_dp
    else
       A(3)=ATAN2(S_IJ,S_JJ)
    endif
    !    A(3)=ATAN2(S_IJ,S_JJ)
   ! AT=0.0_dp;AT(3)=A(3);

   ! CALL GEO_ROT(T10,T1,AT,T2)
    !   T2=T1
    !   write(16,*) t2-entb


  END SUBROUTINE COMPUTE_ENTRANCE_ANGLE

  SUBROUTINE COMPUTE_ENTRANCE_ANGLE_bmad(ENTL,ENTB,A)
    ! COMPUTES BMAD'S ANGLES IN FRAME OF ENTL
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: ENTL(3,3), ENTB(3,3)
    REAL(DP), INTENT(OUT) :: A(3)
    REAL(DP) T1(3,3),T10(3,3),T2(3,3),S_IJ,S_JJ
    REAL(DP) AT(3)

    T1=ENTL
    T10=ENTL
    T2=ENTL


    CALL COMPUTE_SCALAR(T1,1,ENTB,3,S_IJ)
    CALL COMPUTE_SCALAR(T1,3,ENTB,3,S_JJ)

    if(S_IJ==0.0_dp.and.S_JJ==0.0_dp) then
       A(2)=0.0_dp
    else
       A(2)=ATAN2(-S_IJ,S_JJ)
    endif
    AT=0.0_dp;AT(2)=A(2);


    CALL GEO_ROT(T10,T1,AT,T2)
    T2=T1
    T10=T1


    CALL COMPUTE_SCALAR(T1,2,ENTB,3,S_IJ)
    CALL COMPUTE_SCALAR(T1,3,ENTB,3,S_JJ)

    if(S_IJ==0.0_dp.and.S_JJ==0.0_dp) then
       A(1)=0.0_dp
    else
       A(1)=ATAN2(-S_IJ,S_JJ)
    endif

 
    AT=0.0_dp;AT(1)=A(1);

    CALL GEO_ROT(T10,T1,AT,T2)
    T2=T1
    T10=T1


    CALL COMPUTE_SCALAR(T1,2,ENTB,1,S_IJ)
    CALL COMPUTE_SCALAR(T1,1,ENTB,1,S_JJ)

    if(S_IJ==0.0_dp.and.S_JJ==0.0_dp) then
       A(3)=0.0_dp
    else
       A(3)=ATAN2(S_IJ,S_JJ)
    endif

  !  AT=0.0_dp;AT(3)=A(3);

    a(1:2)=-a(1:2)




  END SUBROUTINE COMPUTE_ENTRANCE_ANGLE_bmad 

  SUBROUTINE  COMPUTE_SCALAR(ENTL,I,ENTB,J,S_IJ) ! Adjusts frames of magnet_chart on the basis of the misalignments
    IMPLICIT NONE
    real(dp),INTENT(in):: ENTL(3,3), ENTB(3,3)
    INTEGER,INTENT(in):: I,J
    real(dp),INTENT(inOUT) :: S_IJ
    INTEGER K

    S_IJ=0.0_dp
    DO K=1,3
       S_IJ=ENTL(I,K)*ENTB(J,K)+S_IJ
    ENDDO

  END SUBROUTINE  COMPUTE_SCALAR


  SUBROUTINE FIND_PATCH_BMAD0(A,ENT,B,EXI,D,ANG)
    ! FINDS PATCH BETWEEN ENT AND EXI : INTERFACED LATER FOR FIBRES
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT):: ENT(3,3),EXI(3,3)
    REAL(DP), INTENT(INOUT):: A(3),B(3),D(3),ANG(3)
    REAL(DP) dd(3)

    D=B-A
    dd=d
    CALL CHANGE_BASIS(Dd,GLOBAL_FRAME,D,ENT)

    CALL COMPUTE_ENTRANCE_ANGLE_bmad(ENT,EXI,ANG)



  END SUBROUTINE FIND_PATCH_BMAD0


  SUBROUTINE FIND_PATCH_B(A,ENT,B,EXI,D,ANG)
    ! FINDS PATCH BETWEEN ENT AND EXI : INTERFACED LATER FOR FIBRES
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT):: ENT(3,3),EXI(3,3)
    REAL(DP), INTENT(INOUT):: A(3),B(3),D(3),ANG(3)
    REAL(DP) dd(3)

    CALL COMPUTE_ENTRANCE_ANGLE(ENT,EXI,ANG)

    
    D=B-A
    dd=d
    CALL CHANGE_BASIS(Dd,GLOBAL_FRAME,D,EXI)


  END SUBROUTINE FIND_PATCH_B




  SUBROUTINE INVERSE_FIND_PATCH(A,ENT,D,ANG,B,EXI) !  used in misalignments of siamese
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT):: ENT(3,3),EXI(3,3)
    REAL(DP), INTENT(INOUT):: A(3),B(3),D(3),ANG(3)
    REAL(DP) DD(3)

    EXI=ENT
    CALL GEO_ROT(EXI,ANG,1,basis=ENT)
    !    CALL GEO_ROTA_no_vec(EXI,ANG,1,basis=ENT)
    CALL CHANGE_BASIS(D,EXI,DD,GLOBAL_FRAME)
    B=A+DD

  END SUBROUTINE INVERSE_FIND_PATCH


  SUBROUTINE INVERSE_FIND_PATCH_bmad(A,ENT,D,ANG,B,EXI) !  used in misalignments of siamese
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT):: ENT(3,3),EXI(3,3)
    REAL(DP), INTENT(INOUT):: A(3),B(3),D(3),ANG(3)
    REAL(DP) DD(3),an(3)


    CALL CHANGE_BASIS(D,ent,DD,GLOBAL_FRAME)
    B=A+DD
    EXI=ENT
    an=0
    an(3)=ang(3)
    CALL GEO_ROT(EXI,AN,1,basis=ENT)
    
    !ENT=EXI
    an=0
    an(1)=-ang(1)
    CALL GEO_ROT(EXI,AN,1,basis=ENT)

    !ENT=EXI
    an=0
    an(2)=-ang(2)
    CALL GEO_ROT(EXI,AN,1,basis=ENT)



  END SUBROUTINE INVERSE_FIND_PATCH_bmad

end module S_FRAME
