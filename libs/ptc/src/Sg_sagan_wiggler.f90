!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module sagan_WIGGLER
  use S_def_all_kinds
  implicit none
  public
  private INTR,INTP,ZEROr_SAGAN,ZEROp_SAGAN
  private ALLOC_SAGAN,KILL_SAGAN,POINTERS_SAGANR,POINTERS_SAGANP
  private copy_el_elp ,copy_elp_el ,copy_el_el
  private scale_SAGANR,scale_SAGANP,PRINT_,READ_
  PRIVATE DRIFTR,DRIFTP,DRIFT
  PRIVATE KICKPATHR,KICKPATHP,KICKPATH
  PRIVATE COMPX_R,COMPX_P,COMPY_R,COMPY_P,COMPZ_R,COMPZ_P,BF_R,BF_P
  PRIVATE COMPZ  !,B_FIELD
  PRIVATE KICKR,KICKP,KICK
  PRIVATE KILL_WIGGLER
  PRIVATE ALLOC_WIGGLER
  PRIVATE ZEROR_W,ZEROP_W,POINTERS_WP
  PRIVATE copy_W_WP,copy_WP_W,copy_W_W,INTR_SAGAN,INTp_SAGAN
  !  PRIVATE SET_R,SET_P,SET_W
  PRIVATE ADJUSTR_WI,ADJUSTP_WI,get_z_wiR,get_z_wiP
  integer :: wiggler_sagan=6
integer, parameter :: hyper_y_plane_y$ = 1, hyper_xy_plane_y$ = 2, hyper_x_plane_y$ = 3
integer, parameter :: hyper_y_plane_x$ = 4, hyper_xy_plane_x$ = 5, hyper_x_plane_x$ = 6
  integer :: limit_sag(2) =(/4,18/) 
 

  INTERFACE get_z_wi
     MODULE PROCEDURE get_z_wiR
     MODULE PROCEDURE get_z_wip
  END INTERFACE

  INTERFACE DRIFT
     MODULE PROCEDURE DRIFTR
     MODULE PROCEDURE DRIFTP
  END INTERFACE

  INTERFACE KICKPATH
     MODULE PROCEDURE KICKPATHR
     MODULE PROCEDURE KICKPATHP
  END INTERFACE

  INTERFACE KICK
     MODULE PROCEDURE KICKR
     MODULE PROCEDURE KICKP
  END INTERFACE

  !  INTERFACE SET_W
  !     MODULE PROCEDURE SET_R
  !     MODULE PROCEDURE SET_P
  !  END INTERFACE

  INTERFACE COMPX
     MODULE PROCEDURE COMPX_R
     MODULE PROCEDURE COMPX_P
  END INTERFACE

  INTERFACE COMPY
     MODULE PROCEDURE COMPY_R
     MODULE PROCEDURE COMPY_P
  END INTERFACE

  INTERFACE COMPZ
     MODULE PROCEDURE COMPZ_R
     MODULE PROCEDURE COMPZ_P
  END INTERFACE

  INTERFACE B_FIELD
     MODULE PROCEDURE BF_R
     MODULE PROCEDURE BF_P
  END INTERFACE

  INTERFACE TRACK
     MODULE PROCEDURE INTR
     MODULE PROCEDURE INTP
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_SAGAN
     MODULE PROCEDURE ALLOC_WIGGLER
  END INTERFACE

  INTERFACE POINTERS_SAGAN
     MODULE PROCEDURE POINTERS_SAGANR
     MODULE PROCEDURE POINTERS_SAGANP
  END INTERFACE

  INTERFACE POINTERS_W
     MODULE PROCEDURE INIT_SAGAN_POINTERS
     MODULE PROCEDURE POINTERS_WP
  END INTERFACE


  INTERFACE KILL
     MODULE PROCEDURE KILL_SAGAN
     MODULE PROCEDURE KILL_WIGGLER
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_el_elp
     MODULE PROCEDURE copy_elp_el
     MODULE PROCEDURE copy_el_el
     MODULE PROCEDURE copy_W_WP
     MODULE PROCEDURE copy_WP_W
     MODULE PROCEDURE copy_W_W
  END INTERFACE

  INTERFACE scale_SAGAN
     MODULE PROCEDURE scale_SAGANR
     MODULE PROCEDURE scale_SAGANP
  END INTERFACE

  INTERFACE PRINT_USER
     MODULE PROCEDURE PRINT_
  END INTERFACE

  INTERFACE READ_USER
     MODULE PROCEDURE READ_
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE ZEROr_SAGAN
     MODULE PROCEDURE ZEROp_SAGAN
     MODULE PROCEDURE ZEROr_W
     MODULE PROCEDURE ZEROp_W
  END INTERFACE

  INTERFACE TRACK_SLICE
     MODULE PROCEDURE INTR_SAGAN
     MODULE PROCEDURE INTP_SAGAN
  END INTERFACE

  INTERFACE ADJUST_WI
     MODULE PROCEDURE ADJUSTR_WI
     MODULE PROCEDURE ADJUSTP_WI
  END INTERFACE

contains

  SUBROUTINE ADJUSTR_WI(EL,X,k,J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(sagan),INTENT(INOUT):: EL
    TYPE(INTERNAL_STATE),OPTIONAL :: K    

    INTEGER, INTENT(IN) :: J


    IF(J==1.and.el%p%dir==-1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    elseIF(J==2.and.el%p%dir==1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    endif
  END SUBROUTINE ADJUSTR_WI

  SUBROUTINE ADJUSTP_WI(EL,X,k,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(saganP),INTENT(INOUT):: EL
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    INTEGER, INTENT(IN) :: J

    IF(J==1.and.el%p%dir==-1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    elseIF(J==2.and.el%p%dir==1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    endif


  END SUBROUTINE ADJUSTP_WI

  !  SUBROUTINE SET_R(EL)
  !    IMPLICIT NONE
  !    TYPE(undu_r),INTENT(INOUT):: EL
  !    INTEGER I,MYPAUSE
  !    IF(ASSOCIATED(EL%K)) THEN
  !       DO I=1,SIZE(EL%A)
  !          if(el%form(i)) then
  !             EL%K(2,I)=SQRT(EL%K(1,I)**2+EL%K(3,I)**2)
  !          else
  !             EL%K(2,I)=SQRT(-EL%K(1,I)**2+EL%K(3,I)**2)
  !          endif
  !       ENDDO
  !
  !
  !    ELSE
  !       WRITE(6,*) " WIGGLER FIELD NOT DEFINED "
  !       I=MYPAUSE(555)
  !    ENDIF
  !  END SUBROUTINE SET_R

  !  SUBROUTINE SET_P(EL)
  !    IMPLICIT NONE
  !    TYPE(undu_p),INTENT(INOUT):: EL
  !    INTEGER I,MYPAUSE
  !    IF(ASSOCIATED(EL%K)) THEN
  !       DO I=1,SIZE(EL%A)
  !          if(el%form(i)) then
  !             EL%K(2,I)=SQRT(EL%K(1,I)**2+EL%K(3,I)**2)
  !          else
  !             EL%K(2,I)=SQRT(-EL%K(1,I)**2+EL%K(3,I)**2)
  !          endif
  !       ENDDO
  !    ELSE
  !       WRITE(6,*) " WIGGLER FIELD NOT DEFINED "
  !       I=MYPAUSE(556)
  !    ENDIF
  !  END SUBROUTINE SET_P

  SUBROUTINE INTR(EL,X,k,mid)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(INOUT):: EL
    TYPE(WORM),OPTIONAL,INTENT(INOUT):: mid
    INTEGER I
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    !    CALL SET_W(EL%W)
    IF(PRESENT(MID)) CALL XMID(MID,X,0)

    DO I=1,EL%P%NST
       call track_slice(el,x,k,i)
       IF(PRESENT(MID)) CALL XMID(MID,X,i)
    ENDDO

    call ADJUST_WI(EL,X,k,2)

  END SUBROUTINE INTR

  SUBROUTINE get_z_wir(EL,i,z)
    IMPLICIT NONE
    TYPE(SAGAN),INTENT(INOUT):: EL
    integer i
    real(dp),INTENT(INOUT):: z
    real(dp) d

    D=EL%L/EL%P%NST
    IF(EL%P%DIR==1) THEN
       Z=(i-1)*d
    ELSE
       Z=EL%L-(i-1)*d
    ENDIF

  end SUBROUTINE get_z_wir

  SUBROUTINE get_z_wip(EL,i,z)
    IMPLICIT NONE
    TYPE(SAGANP),INTENT(INOUT):: EL
    integer i
    TYPE(REAL_8),INTENT(INOUT):: z
    TYPE(REAL_8) d

    CALL ALLOC(D)

    D=EL%L/EL%P%NST
    IF(EL%P%DIR==1) THEN
       Z=(i-1)*d
    ELSE
       Z=EL%L-(i-1)*d
    ENDIF

    CALL KILL(D)

  end SUBROUTINE get_z_wip

  SUBROUTINE INTR_SAGAN(EL,X,k,i)
    IMPLICIT NONE
    integer ipause, mypause
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(INOUT):: EL
    integer,INTENT(IN):: I
    real(dp) Z
    real(dp) D,DH
    real(dp) D1,D2,DK1,DK2
    real(dp) DF(4),DK(4)
    INTEGER J
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    !    CALL SET_W(EL%W)

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/2.0_dp/EL%P%NST
       D=EL%L/EL%P%NST
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       Z=Z+EL%P%DIR*DH
       CALL DRIFT(EL,DH,Z,1,X,k)
       CALL DRIFT(EL,DH,Z,2,X,k)
       CALL KICKPATH(EL,DH,X,k)
       CALL KICK(EL,D,Z,X,k)
       CALL KICKPATH(EL,DH,X,k)
       CALL DRIFT(EL,DH,Z,2,X,k)
       CALL DRIFT(EL,DH,Z,1,X,k)
       !       Z=Z+EL%P%DIR*DH

    CASE(4)
       D=EL%L/EL%P%NST

       DK1=D*FK1
       D1=DK1/2.0_dp
       DK2=D*FK2
       D2=DK2/2.0_dp
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       Z=Z+EL%P%DIR*D1
       CALL DRIFT(EL,D1,Z,1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL DRIFT(EL,D1,Z,1,X,k)
       Z=Z+EL%P%DIR*D1+D2
       CALL DRIFT(EL,D2,Z,1,X,k)
       CALL DRIFT(EL,D2,Z,2,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL KICK(EL,DK2,Z,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL DRIFT(EL,D2,Z,2,X,k)
       CALL DRIFT(EL,D2,Z,1,X,k)
       Z=Z+EL%P%DIR*(D1+D2)
       CALL DRIFT(EL,D1,Z,1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL DRIFT(EL,D1,Z,1,X,k)

    CASE(6)
       DO j =1,4
          DK(j)=EL%L*YOSK(J)/EL%P%NST
          DF(j)=DK(j)/2.0_dp
       ENDDO
       D=EL%L/EL%P%NST
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF
       DO J=4,1,-1
          Z=Z+EL%P%DIR*DF(J)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO
       DO J=2,4
          Z=Z+EL%P%DIR*DF(J)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO

    CASE DEFAULT
       WRITE(6,*) " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       ipause=mypause(357)
    END SELECT


  END SUBROUTINE INTR_SAGAN



  SUBROUTINE INTP(EL,X,k)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(INOUT):: EL
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    INTEGER I

    !    CALL SET_W(EL%W)


    DO I=1,EL%P%NST
       call track_slice(el,x,k,i)
    ENDDO

    call ADJUST_WI(EL,X,k,2)

  END SUBROUTINE INTP

  SUBROUTINE INTP_SAGAN(EL,X,k,i)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(INOUT):: EL
    integer,INTENT(IN):: I
    TYPE(REAL_8) Z
    TYPE(REAL_8) D,DH
    TYPE(REAL_8) D1,D2,DK1,DK2
    TYPE(REAL_8) DF(4),DK(4)
    INTEGER J
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    CALL ALLOC(Z,D,DH,D1,D2,DK1,DK2)
    CALL ALLOC(DF,4)
    CALL ALLOC(DK,4)

    !    CALL SET_W(EL%W)
    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/2.0_dp/EL%P%NST
       D=EL%L/EL%P%NST

       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       Z=Z+EL%P%DIR*DH
       CALL DRIFT(EL,DH,Z,1,X,k)
       CALL DRIFT(EL,DH,Z,2,X,k)
       CALL KICKPATH(EL,DH,X,k)
       CALL KICK(EL,D,Z,X,k)
       CALL KICKPATH(EL,DH,X,k)
       CALL DRIFT(EL,DH,Z,2,X,k)
       CALL DRIFT(EL,DH,Z,1,X,k)

    CASE(4)
       D=EL%L/EL%P%NST

       DK1=D*FK1
       D1=DK1/2.0_dp
       DK2=D*FK2
       D2=DK2/2.0_dp

       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       Z=Z+EL%P%DIR*D1
       CALL DRIFT(EL,D1,Z,1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL DRIFT(EL,D1,Z,1,X,k)
       Z=Z+EL%P%DIR*D1+D2
       CALL DRIFT(EL,D2,Z,1,X,k)
       CALL DRIFT(EL,D2,Z,2,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL KICK(EL,DK2,Z,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL DRIFT(EL,D2,Z,2,X,k)
       CALL DRIFT(EL,D2,Z,1,X,k)
       Z=Z+EL%P%DIR*(D1+D2)
       CALL DRIFT(EL,D1,Z,1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL DRIFT(EL,D1,Z,2,X,k)
       CALL DRIFT(EL,D1,Z,1,X,k)

    CASE(6)
       DO j =1,4
          DK(j)=EL%L*YOSK(J)/EL%P%NST
          DF(j)=DK(j)/2.0_dp
       ENDDO
       D=EL%L/EL%P%NST
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF
       DO J=4,1,-1
          Z=Z+EL%P%DIR*DF(J)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO
       DO J=2,4
          Z=Z+EL%P%DIR*DF(J)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL DRIFT(EL,DF(J),Z,2,X,k)
          CALL DRIFT(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO

    CASE DEFAULT
       WRITE(6,*) " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       ipause=mypause(357)
    END SELECT

    CALL KILL(Z,D,DH,D1,D2,DK1,DK2)
    CALL KILL(DF,4)
    CALL KILL(DK,4)



  END SUBROUTINE INTP_SAGAN

  SUBROUTINE ZEROR_SAGAN(EL,I)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%INTERNAL))  THEN
          DEALLOCATE(EL%INTERNAL)
          DEALLOCATE(EL%n_min)
          EL%W=-1
          DEALLOCATE(EL%W)
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%INTERNAL)
       NULLIFY(EL%n_min)
       NULLIFY(EL%W)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZEROR_SAGAN

  SUBROUTINE ZEROp_SAGAN(EL,I)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED KILL AND DEASSOCIATE
       IF(ASSOCIATED(EL%INTERNAL))  THEN
          CALL KILL(EL)              ! FPP DEALLOCATION FIRST OBVIOUSLY
          EL%W=-1
          DEALLOCATE(EL%INTERNAL)
          DEALLOCATE(EL%n_min)
          DEALLOCATE(EL%W)

       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%W)
       NULLIFY(EL%n_min)
       NULLIFY(EL%INTERNAL)  !!! was not there January 2014
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
       !
    endif

  END SUBROUTINE ZEROp_SAGAN

  SUBROUTINE ZEROR_W(EL,I)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%K))  THEN
          DEALLOCATE(EL%A)
          DEALLOCATE(EL%F,EL%x0,EL%y0)
          DEALLOCATE(EL%offset)
          DEALLOCATE(EL%FORM)
          DEALLOCATE(EL%K)
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%A)
       NULLIFY(EL%F,EL%x0,EL%y0)
       NULLIFY(EL%offset)
       NULLIFY(EL%FORM)
       NULLIFY(EL%K)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZEROR_W

  SUBROUTINE ZEROp_W(EL,I)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED KILL AND DEASSOCIATE
       IF(ASSOCIATED(EL%K))  THEN
          CALL KILL(EL)              ! FPP DEALLOCATION FIRST OBVIOUSLY
          DEALLOCATE(EL%A)
          DEALLOCATE(EL%F,EL%x0,EL%y0)
          DEALLOCATE(EL%offset)
          DEALLOCATE(EL%FORM)
          DEALLOCATE(EL%K)
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%A)
       NULLIFY(EL%F,EL%x0,EL%y0)
       NULLIFY(EL%offset)
       NULLIFY(EL%FORM)
       NULLIFY(EL%K)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
       !
    endif

  END SUBROUTINE ZEROp_W


  SUBROUTINE copy_el_elp(EL,ELP)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(in)::EL
    TYPE(SAGANP), INTENT(inout)::ELP

    ELP%INTERNAL(1)    =EL%INTERNAL(1)
    ELP%INTERNAL(2)    =EL%INTERNAL(2)
    ELP%INTERNAL(3)    =EL%INTERNAL(3)
    ELP%INTERNAL(4)    =EL%INTERNAL(4)
    ELP%INTERNAL(5)    =EL%INTERNAL(5)
    ELP%INTERNAL(6)    =EL%INTERNAL(6)
    ELP%n_min    =EL%n_min
    CALL COPY(EL%W,ELP%W)
    !  COPY CODING HERE NO ALLOCATION OF POINTERS OR POLYMORPH NEEDED
    !  IF DONE CORRECTLY

  END SUBROUTINE copy_el_elp

  SUBROUTINE copy_elp_el(EL,ELP)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(in)::EL
    TYPE(SAGAN), INTENT(inout)::ELP

    ELP%INTERNAL(1)    =EL%INTERNAL(1)
    ELP%INTERNAL(2)    =EL%INTERNAL(2)
    ELP%INTERNAL(3)    =EL%INTERNAL(3)
    ELP%INTERNAL(4)    =EL%INTERNAL(4)
    ELP%INTERNAL(5)    =EL%INTERNAL(5)
    ELP%INTERNAL(6)    =EL%INTERNAL(6)
    ELP%n_min    =EL%n_min
    CALL COPY(EL%W,ELP%W)

  END SUBROUTINE copy_elp_el

  SUBROUTINE copy_el_el(EL,ELP)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(in)::EL
    TYPE(SAGAN), INTENT(inout)::ELP

    ELP%INTERNAL(1)    =EL%INTERNAL(1)
    ELP%INTERNAL(2)    =EL%INTERNAL(2)
    ELP%INTERNAL(3)    =EL%INTERNAL(3)
    ELP%INTERNAL(4)    =EL%INTERNAL(4)
    ELP%INTERNAL(5)    =EL%INTERNAL(5)
    ELP%INTERNAL(6)    =EL%INTERNAL(6)
    ELP%n_min    =EL%n_min
    !  COPY CODING HERE NO ALLOCATION OF POINTERS
    CALL COPY(EL%W,ELP%W)

  END SUBROUTINE copy_el_el



  SUBROUTINE copy_W_W(EL,ELP)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(in)::EL
    TYPE(undu_r), INTENT(inout)::ELP

    INTEGER I,J

    if(associated(el%K)) then
       CALL POINTERS_W(ELP,SIZE(EL%A))
       DO I=1,3
          DO J=1,SIZE(EL%A)
             ELP%K(I,J)   =EL%K(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%A)
          ELP%A(I)    =EL%A(I)
          ELP%F(I)    =EL%F(I)
          ELP%x0(I)    =EL%x0(I)
          ELP%y0(I)    =EL%y0(I)
          ELP%FORM(I) =EL%FORM(I)
       ENDDO
       ELP%offset   =EL%offset
    endif

  END SUBROUTINE copy_W_W

  SUBROUTINE copy_W_WP(EL,ELP)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(in)::EL
    TYPE(undu_p), INTENT(inout)::ELP

    INTEGER I,J

    if(associated(el%K)) then
       CALL POINTERS_W(ELP,SIZE(EL%A))
       DO I=1,3
          DO J=1,SIZE(EL%A)
             ELP%K(I,J)   =EL%K(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%A)
          ELP%A(I)    =EL%A(I)
          ELP%F(I)    =EL%F(I)
          ELP%x0(I)    =EL%x0(I)
          ELP%y0(I)    =EL%y0(I)
          ELP%FORM(I) =EL%FORM(I)
       ENDDO
       ELP%offset   =EL%offset
    endif

  END SUBROUTINE copy_W_WP

  SUBROUTINE copy_WP_W(EL,ELP)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(in)::EL
    TYPE(undu_r), INTENT(inout)::ELP

    INTEGER I,J

    if(associated(el%K)) then
       CALL POINTERS_W(ELP,SIZE(EL%A))
       DO I=1,3
          DO J=1,SIZE(EL%A)
             ELP%K(I,J)   =EL%K(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%A)
          ELP%A(I)    =EL%A(I)
          ELP%x0(I)    =EL%x0(I)
          ELP%y0(I)    =EL%y0(I)
          ELP%F(I)    =EL%F(I)
          ELP%FORM(I) =EL%FORM(I)
       ENDDO
       ELP%offset   =EL%offset
    endif

  END SUBROUTINE copy_WP_W

  SUBROUTINE POINTERS_SAGANR(EL)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(INOUT)::EL

    ALLOCATE(EL%INTERNAL(6))
    EL%INTERNAL=0.0_dp
    allocate(EL%n_min)
     EL%n_min=wiggler_sagan
    ALLOCATE(EL%W)
    !CALL POINTERS_W(EL%W)
    el%w=0

    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_SAGANR

  SUBROUTINE POINTERS_SAGANP(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL

    ALLOCATE(EL%INTERNAL(6))
    allocate(EL%n_min)
     EL%n_min=wiggler_sagan
    ALLOCATE(EL%W)
    !CALL POINTERS_W(EL%W)
    el%w=0
    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_SAGANP

  SUBROUTINE INIT_SAGAN_POINTERS(EL,N)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(INOUT)::EL
    INTEGER, INTENT(IN)::N

    IF(ASSOCIATED(EL%A)) THEN
       EL=-1
    ENDIF
    EL=0
    ALLOCATE(EL%A(N))
    ALLOCATE(EL%F(N),EL%x0(N),EL%y0(N))
    ALLOCATE(EL%offset)
    ALLOCATE(EL%FORM(N))
    ALLOCATE(EL%K(3,N))
    EL%K=0.0_dp
    EL%x0=0.0_dp
    EL%y0=0.0_dp
    EL%A=0.0_dp
    EL%F=0.0_dp
    EL%offset=0.0_dp
    EL%FORM=0

  END SUBROUTINE INIT_SAGAN_POINTERS

  SUBROUTINE POINTERS_WP(EL,N)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER, INTENT(IN)::N

    IF(ASSOCIATED(EL%A)) THEN
       CALL KILL(EL)
       EL=-1
    ENDIF
    EL=0
    ALLOCATE(EL%A(N))
    ALLOCATE(EL%F(N),EL%x0(N),EL%y0(N))
    ALLOCATE(EL%offset)
    ALLOCATE(EL%K(3,N))
    ALLOCATE(EL%FORM(N))
    EL%FORM=0
    CALL ALLOC(EL)

    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_WP



  SUBROUTINE ALLOC_SAGAN(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL
    CALL ALLOC(EL%INTERNAL,6)
 
    ! CALL ALLOC(EL%W)
    ! ALLOC INTERNAL POLYMORPHS IF ANY
  END SUBROUTINE ALLOC_SAGAN

  SUBROUTINE ALLOC_WIGGLER(EL)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER I,J
    ! ALLOC INTERNAL POLYMORPHS IF ANY
    IF(ASSOCIATED(EL%K)) THEN  !DAVID
       DO I=1,3
          DO J=1,SIZE(EL%A)
             CALL ALLOC(EL%K(I,J));
          ENDDO
       ENDDO
       CALL ALLOC(EL%A,SIZE(EL%A));
       CALL ALLOC(EL%F,SIZE(EL%A));
       CALL ALLOC(EL%x0,SIZE(EL%A));
       CALL ALLOC(EL%y0,SIZE(EL%A));
       CALL ALLOC(EL%offset);
    ENDIF
  END SUBROUTINE ALLOC_WIGGLER


  SUBROUTINE KILL_SAGAN(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL

    CALL KILL(EL%INTERNAL,6)
    CALL KILL(EL%W)
    ! KILL INTERNAL POLYMORPHS IF ANY

  END SUBROUTINE KILL_SAGAN

  SUBROUTINE KILL_WIGGLER(EL)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER I,J
    ! KILL INTERNAL POLYMORPHS IF ANY
    IF(ASSOCIATED(EL%K)) then! DAVID
       DO I=1,3
          DO J=1,SIZE(EL%A)
             CALL KILL(EL%K(I,J));
          ENDDO
       ENDDO
       CALL KILL(EL%A,SIZE(EL%A));
       CALL KILL(EL%x0,SIZE(EL%A));
       CALL KILL(EL%y0,SIZE(EL%A));
       CALL KILL(EL%F,SIZE(EL%A));
       CALL KILL(EL%offset);
    ENDIF
  END SUBROUTINE KILL_WIGGLER

  SUBROUTINE reset_WI(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL

    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS

    CALL resetpoly_R31N(EL%INTERNAL,6)
    CALL reset_WIG(EL%W)

  END SUBROUTINE reset_WI

  SUBROUTINE reset_WIG(EL)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER I,J
    IF(ASSOCIATED(EL%K)) THEN ! DAVID
       DO I=1,3
          DO J=1,SIZE(EL%A)
             CALL resetpoly_R31(EL%K(I,J));
          ENDDO
       ENDDO

       CALL resetpoly_R31N(EL%A,SIZE(EL%A))
       CALL resetpoly_R31N(EL%F,SIZE(EL%A))
       CALL resetpoly_R31N(EL%x0,SIZE(EL%A))
       CALL resetpoly_R31N(EL%y0,SIZE(EL%A))
       CALL resetpoly_R31(EL%offset)
    ENDIF
  END SUBROUTINE reset_WIG


  SUBROUTINE  ELp_POL_SAGAN(S2,S2R,S1,DONEIT)
    implicit none
    integer ipause,mypause
    type (POL_BLOCK),INTENT(IN):: S1
    TYPE(SAGANp),INTENT(inOUT):: S2
    TYPE(SAGAN),INTENT(inOUT):: S2R
    LOGICAL(lp),INTENT(inOUT)::  DONEIT
    integer i

    ! ONE CAN LINK INTERNAL POLYMORPHS TO PART OF POL_BLOCK WHICH IS NOT USED
    ! HERE THE VARIABLE "INTERNAL" IS LINKED TO VOLT
    ! We also linked it to the pol_block2


    !   IF(S1%IVOLT>0) THEN
    !      s2%INTERNAL%I=S1%IVOLT+S1%NPARA
    !      s2%INTERNAL%S=S1%SVOLT
    !      s2%INTERNAL%KIND=3
    !!      DONEIT=.TRUE.
    !      IF(S1%SET_TPSAFIT) THEN
    !         s2%INTERNAL%R=s2%INTERNAL%R+s2%INTERNAL%S*s1%TPSAFIT(S1%IVOLT)
    !      ENDIF
    !    ENDIF
    !  or try
    DO I=1,6
       IF(S1%SAGAN%Iinternal(I)>0) THEN
          s2%INTERNAL(I)%I=S1%SAGAN%Iinternal(I)+S1%NPARA
          s2%INTERNAL(I)%S=S1%SAGAN%Sinternal(I)
          s2%INTERNAL(I)%KIND=3
          if(S1%SAGAN%Iinternal(I)>c_%np_pol) c_%np_pol=S1%SAGAN%Iinternal(I)
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN

             s2%INTERNAL(I)%R=s2%INTERNAL(I)%R+s2%INTERNAL(I)%S*s1%TPSAFIT(S1%SAGAN%Iinternal(I))
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2R%INTERNAL(I)=s2%INTERNAL(I)%R
          ENDIF
       ENDIF
    ENDDO
    if(size(s2%w%A)>size(S1%SAGAN%w%ia)) then
       write(6,*) " Pol_block for wiggler must be made bigger ",size(s2%w%A),size(S1%SAGAN%w%ia)
       ipause=mypause(121)
    endif

    do i=1,SIZE(s2%w%A)
       IF(S1%SAGAN%w%ia(i)>0) THEN
          s2%w%a(i)%I=S1%SAGAN%w%ia(i)+S1%NPARA
          s2%w%a(i)%S=S1%SAGAN%w%sa(i)
          s2%w%a(i)%KIND=3
          if(S1%SAGAN%w%ia(i)>c_%np_pol) c_%np_pol=S1%SAGAN%w%ia(i)
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN
             s2%w%a(i)%R=s2%w%a(i)%R+s2%w%a(i)%S*s1%TPSAFIT(S1%SAGAN%w%ia(i))
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2R%w%a(i)=s2%w%a(i)%R
          ENDIF
       ENDIF
    enddo

  end SUBROUTINE  ELp_POL_SAGAN

  SUBROUTINE  scale_SAGANR(S2,P0C_OLD,P0C_NEW,power)
    implicit none
    TYPE(SAGAN),INTENT(inOUT):: S2
    real(dp),INTENT(IN)::  P0C_OLD,P0C_NEW
    integer, INTENT(IN):: power
    INTEGER I
    ! EXAMPLE

    !    S2%INTERNAL= S2%INTERNAL*P0C_OLD/P0C_NEW
    DO I=1,SIZE(S2%W%A)
       S2%W%A(I)=S2%W%A(I)*(P0C_OLD/P0C_NEW)**power
    ENDDO

  end SUBROUTINE  scale_SAGANR

  SUBROUTINE  scale_SAGANP(S2,P0C_OLD,P0C_NEW,power)
    implicit none
    TYPE(SAGANp),INTENT(inOUT):: S2
    real(dp),INTENT(IN)::  P0C_OLD,P0C_NEW
    integer, INTENT(IN):: power
    INTEGER I
    ! EXAMPLE

    !    S2%INTERNAL= S2%INTERNAL*P0C_OLD/P0C_NEW
    DO I=1,SIZE(S2%W%A)
       S2%W%A(I)=S2%W%A(I)*(P0C_OLD/P0C_NEW)**power
    ENDDO

  end SUBROUTINE  scale_SAGANP

  ! split drifts
  SUBROUTINE DRIFTR(EL,L,Z,PLANE,X,k)
    IMPLICIT NONE
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(INOUT):: X(6)
    real(dp), INTENT(IN):: L,Z
    INTEGER, INTENT(IN)::PLANE
    real(dp) PZ,A,B,AP,BP
    TYPE(INTERNAL_STATE),OPTIONAL :: K


    IF(PLANE==1) THEN
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)-A
       X(4)=X(4)-AP
       if(k%TIME) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(1)=X(1)+L*X(2)/pz
          X(6)=X(6)+((X(2)*X(2))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(1)=X(1)+L*X(2)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(2)*X(2))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)+A
       X(4)=X(4)+AP
    ELSE
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)-BP
       X(4)=X(4)-B
       if(k%TIME) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(3)=X(3)+L*X(4)/pz
          X(6)=X(6)+((X(4)*X(4))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(3)=X(3)+L*X(4)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)+BP
       X(4)=X(4)+B
    ENDIF
  END SUBROUTINE DRIFTR

  SUBROUTINE DRIFTP(EL,L,Z,PLANE,X,k)
    IMPLICIT NONE
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8), INTENT(IN):: L,Z
    INTEGER, INTENT(IN)::PLANE
    TYPE(REAL_8) PZ,A,B,AP,BP
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    CALL ALLOC(PZ,A,B,AP,BP)
    IF(PLANE==1) THEN
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)-A
       X(4)=X(4)-AP
       if(k%TIME) then
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(1)=X(1)+L*X(2)/pz
          X(6)=X(6)+((X(2)*X(2))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(1)=X(1)+L*X(2)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(2)*X(2))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)+A
       X(4)=X(4)+AP
    ELSE
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)-BP
       X(4)=X(4)-B
       if(k%TIME) then
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(3)=X(3)+L*X(4)/pz
          X(6)=X(6)+((X(4)*X(4))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(3)=X(3)+L*X(4)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)+BP
       X(4)=X(4)+B
    ENDIF

    CALL KILL(PZ,A,B,AP,BP)

  END SUBROUTINE DRIFTP

  SUBROUTINE KICKPATHR(EL,L,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: L
    real(dp) PZ,PZ0,DPZ
    TYPE(INTERNAL_STATE),OPTIONAL :: K
    ! ETIENNE
    IF(EL%P%EXACT) THEN
       if(k%TIME) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2-X(2)**2-X(4)**2)
          PZ0=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp/EL%P%BETA0+X(5))/PZ +k%TOTALPATH*L/EL%P%BETA0
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)   ! WRONG NOT SYMPLECTIC 2015.8.11
          X(6)=X(6)-((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          PZ=ROOT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          PZ0=1.0_dp+X(5)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp+X(5))/PZ +k%TOTALPATH*L 
          PZ=ROOT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          X(6)=X(6)-(L/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
    ELSE
       if(k%TIME) then
          X(6)=X(6)+k%TOTALPATH*L/EL%P%BETA0
       else
          X(6)=X(6)+k%TOTALPATH*L
       endif
    ENDIF

  END SUBROUTINE KICKPATHR

  SUBROUTINE KICKPATHP(EL,L,X,k)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: L
    TYPE(REAL_8) PZ,PZ0,DPZ
    TYPE(INTERNAL_STATE),OPTIONAL :: K
    ! ETIENNE
    IF(EL%P%EXACT) THEN
       CALL ALLOC(PZ,PZ0,DPZ)
       if(k%TIME) then
          PZ=sqrt(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2-X(2)**2-X(4)**2)
          PZ0=sqrt(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp/EL%P%BETA0+X(5))/PZ +k%TOTALPATH*L/EL%P%BETA0
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)   ! WRONG NOT SYMPLECTIC 2015.8.11
          X(6)=X(6)-((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          PZ=sqrt((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          PZ0=1.0_dp+X(5)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp+X(5))/PZ +k%TOTALPATH*L 
          PZ=SQRT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          X(6)=X(6)-(L/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL KILL(PZ,PZ0,DPZ)
    ELSE
       if(k%TIME) then
          X(6)=X(6)+k%TOTALPATH*L/EL%P%BETA0
       else
          X(6)=X(6)+k%TOTALPATH*L
       endif
    ENDIF

  END SUBROUTINE KICKPATHP

  !   X_PLANE
  SUBROUTINE COMPX_R(EL,Z,X,A,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp


   DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =   EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))/2.0_dp)**2*EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2*0.5_dp &
                *COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
!          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_plane_x$) THEN
          A =  EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))*0.5_dp)**2 *EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2* 0.5_dp &
               *COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          A =  EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
!          B =  EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_x_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =  EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
!          A =  EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       endif
     enddo

    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 
    
  END SUBROUTINE COMPX_R

  SUBROUTINE COMPX_P(EL,Z,X,A,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp
! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X

   DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =   EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))/2.0_dp)**2*EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2*0.5_dp &
                *cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
!          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_plane_x$) THEN
          A =  EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))*0.5_dp)**2 *EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2* 0.5_dp &
               *cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          A =  EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
!          B =  EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_x_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =  EL%W%A(I)*EL%W%K(3,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
!          A =  EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       endif
     enddo

    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 
  END SUBROUTINE COMPX_P

  !   Y_PLANE
  SUBROUTINE COMPY_R(EL,Z,X,A,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp
 
    DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i))  &
               *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + B
!          A =  -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       elseif (EL%W%FORM(I) == hyper_xy_plane_y$) THEN
          A =  -EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2 *EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  -EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
!          B =  -EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_x_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) & 
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2*EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
!          B =   EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + B
       endif
ENDDO
 


    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 
  END SUBROUTINE COMPY_R

  SUBROUTINE COMPY_P(EL,Z,X,A,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp

    DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i))  &
               *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + B
!          A =  -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       elseif (EL%W%FORM(I) == hyper_xy_plane_y$) THEN
          A =  -EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2 *EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  -EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
!          B =  -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_x_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) & 
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2*EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
!          B =   EL%W%A(I)*EL%W%K(3,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + B
       endif
ENDDO


    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 

  END SUBROUTINE COMPY_P

  !   Z_PLANE

  SUBROUTINE COMPZ_R(EL,Z,X,A,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp

    DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_x$) THEN
          A =  -EL%W%A(I)*EL%W%K(1,i)**2*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =  EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))   &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  EL%W%A(I)*EL%W%K(1,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(3,i) + A
       elseif (EL%W%FORM(I) == hyper_x_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) * &
                sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
       ELSEif (EL%W%FORM(I) == hyper_y_plane_y$) THEN
          A =  -EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          B =  -EL%W%A(I)*EL%W%K(2,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)**2*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
                *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
 !         B =  -EL%W%A(I)*EL%W%K(2,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
 !               COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_x_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =   EL%W%A(I)*EL%W%K(2,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
       endif
    ENDDO


    a=a-EL%W%offset
    A=A*EL%P%CHARGE*EL%P%DIR 
    B=B*EL%P%CHARGE*EL%P%DIR 
  END SUBROUTINE COMPZ_R


  SUBROUTINE eval_thin_q(EL,Q,n)
    IMPLICIT NONE
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(OUT):: q
    INTEGER I,n

    q=0.0_dp
    n=el%n_min
    DO I=1,SIZE(EL%W%A)

        q = q + (EL%W%A(I)**2/EL%W%K(3,i)) * (EL%W%K(1,i)**2+EL%W%K(2,i)**2)/2.d0

    ENDDO

  END SUBROUTINE eval_thin_q

  SUBROUTINE COMPZ_P(EL,Z,X,A,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: A,B
    INTEGER I
 
    A=0.0_dp
    B=0.0_dp
    DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_x$) THEN
          A =  -EL%W%A(I)*EL%W%K(1,i)**2*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =  EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))   &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  EL%W%A(I)*EL%W%K(1,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(3,i) + A
       elseif (EL%W%FORM(I) == hyper_x_plane_x$) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) * &
                sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
       ELSEif (EL%W%FORM(I) == hyper_y_plane_y$) THEN
          A =  -EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          B =  -EL%W%A(I)*EL%W%K(2,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)**2*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
                *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
 !         B =  -EL%W%A(I)*EL%W%K(2,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
 !               COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_x_plane_y$) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =   EL%W%A(I)*EL%W%K(2,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
       endif
    ENDDO
    a=a-EL%W%offset
    A=A*EL%P%CHARGE*EL%P%DIR 
    B=B*EL%P%CHARGE*EL%P%DIR 
 
  END SUBROUTINE COMPZ_P

  SUBROUTINE BF_R(EL,Z,X,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: B(3)
    INTEGER I
    B=0.0_dp

! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X

 
    DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_x$) THEN


          B(1) = EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_plane_x$) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_plane_x$) THEN
          B(1) = EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = -EL%W%A(I)*EL%W%K(2,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
 
       ELSEif (EL%W%FORM(I) == hyper_y_plane_y$) THEN


          B(1) = -EL%W%A(I)*EL%W%K(1,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_plane_y$) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_plane_y$) THEN
          B(1) = EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
       else
          print *, 'ERROR IN BF_R: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
    ENDDO

       b(2)=b(2)+el%w%offset
 
  END SUBROUTINE BF_R

  SUBROUTINE BF_P(EL,Z,X,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: B(3)
    INTEGER I
    B(1)=0.0_dp;B(2)=0.0_dp;B(3)=0.0_dp;

! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X
 
 
    DO I=1,SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_plane_x$) THEN


          B(1) = EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_plane_x$) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_plane_x$) THEN
          B(1) = EL%W%A(I)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = -EL%W%A(I)*EL%W%K(2,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
 
 
       ELSEif (EL%W%FORM(I) == hyper_y_plane_y$) THEN


          B(1) = -EL%W%A(I)*EL%W%K(1,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_plane_y$) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_plane_y$) THEN
          B(1) = EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
       else
          print *, 'ERROR IN BF_R: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
     enddo

 
    b(2)=b(2)+el%w%offset
 

  END SUBROUTINE BF_P

  SUBROUTINE KICKR(EL,L,Z,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: L,Z
    real(dp) A,B
    TYPE(INTERNAL_STATE),OPTIONAL :: K







    CALL COMPZ(EL,Z,X,A,B)
    X(2)=X(2)+L*A
    X(4)=X(4)+L*B


  END SUBROUTINE KICKR

  SUBROUTINE KICKP(EL,L,Z,X,k)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: L,Z
    TYPE(REAL_8) A,B,AP,BP
    TYPE(INTERNAL_STATE),OPTIONAL :: K


    call alloc(A,B,AP,BP)


    CALL COMPZ(EL,Z,X,A,B)
    X(2)=X(2)+L*A
    X(4)=X(4)+L*B

    call KILL(A,B,AP,BP)

  END SUBROUTINE KICKP

  SUBROUTINE PRINT_(EL,MF)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(INOUT)::EL
    INTEGER MF,I

    !    WRITE(MF,*) EL%INTERNAL
    WRITE(MF,100) "NUMBER OF TERMS ",SIZE(EL%W%A)
    DO I=1,SIZE(EL%W%A)
       WRITE(MF,101) " A = ",EL%W%A(I), " K = ",EL%W%K(1,I),EL%W%K(2,I),EL%W%K(3,I), &
            & " PHASE = ",EL%W%F(I)," FORM = ",EL%W%FORM(I)
    ENDDO
    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS
100 FORMAT(A16,(1X,I4))
101 FORMAT(A5,(1X,g21.14),A5,3(1X,g21.14),A9,(1X,g21.14),A11,I3)

  END SUBROUTINE PRINT_

  SUBROUTINE READ_(EL,MF)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(INOUT)::EL
    INTEGER MF,I

    READ(MF,100) I

    CALL INIT_SAGAN_POINTERS(EL%W,I)
    DO I=1,SIZE(EL%W%A)
       READ(MF,101) EL%W%A(I), EL%W%K(1,I),EL%W%K(2,I),EL%W%K(3,I),EL%W%F(I),EL%W%FORM(I)
    ENDDO

    !   READ(MF,*) EL%INTERNAL
    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS
100 FORMAT(16X,(1X,I4))
101 FORMAT(5X,(1X,g21.14),5X,3(1X,g21.14),9X,(1X,g21.14),11X,I3)

  END SUBROUTINE READ_



end module sagan_WIGGLER
