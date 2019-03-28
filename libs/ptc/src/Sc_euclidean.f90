!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module S_euclidean
  use S_extend_poly
  use s_extend_poly, only : PRTP ! LD: 22.03.2019
  implicit none
  public

  PRIVATE       transr,transp !,TRANSP_P
  PRIVATE       ROT_XZR,ROT_XZP !,ROT_XZP_P
  PRIVATE       ROT_YZR,ROT_YZP !,ROT_YZP_P
  PRIVATE       ROT_XYR,ROT_XYP !,ROT_XYP_P
  private zero_r_z,zero_r_xy,zero_T_XYZ,zero_E_GENERAL,zero_E_GENERAL_s

  INTERFACE TRANS
     MODULE PROCEDURE transR
     MODULE PROCEDURE transp
     !     MODULE PROCEDURE TRANSP_P   ! New stuff
     !     MODULE PROCEDURE transS
     !     MODULE PROCEDURE TRANSS_S   ! New stuff
  END INTERFACE

  INTERFACE ROT_XZ
     MODULE PROCEDURE ROT_XZR
     MODULE PROCEDURE ROT_XZP
     !     MODULE PROCEDURE ROT_XZP_P
     !     MODULE PROCEDURE ROT_XZS
     !     MODULE PROCEDURE ROT_XZS_S
  END INTERFACE
  INTERFACE ROT_YZ
     MODULE PROCEDURE ROT_YZR
     MODULE PROCEDURE ROT_YZP
     !     MODULE PROCEDURE ROT_YZP_P
     !     MODULE PROCEDURE ROT_YZS
     !     MODULE PROCEDURE ROT_YZS_S
  END INTERFACE
  INTERFACE ROT_XY
     MODULE PROCEDURE ROT_XYR
     MODULE PROCEDURE ROT_XYP
     !     MODULE PROCEDURE ROT_XYP_P
     !     MODULE PROCEDURE ROT_XYS
     !     MODULE PROCEDURE ROT_XYS_S
  END INTERFACE

  TYPE R_XY
     REAL(DP) A(2)
  END TYPE R_XY

  TYPE R_Z
     REAL(DP) A
  END TYPE R_Z

  TYPE T_XYZ
     LOGICAL(LP) SIXTRACK
     REAL(DP) L_DESIGN,DL_SIXTRACK
     REAL(DP) D(3),DL
  END TYPE T_XYZ

  TYPE E_GENERAL
     INTEGER KIND
     TYPE(R_XY) T1
     TYPE(R_Z)  T2
     TYPE(T_XYZ) T3
  END TYPE E_GENERAL

  INTERFACE init
     MODULE PROCEDURE zero_r_xy
     MODULE PROCEDURE zero_r_z
     MODULE PROCEDURE zero_T_XYZ
     MODULE PROCEDURE zero_E_GENERAL
     MODULE PROCEDURE zero_E_GENERAL_s
  END INTERFACE


CONTAINS


  subroutine print_e_general(e,mf)
    IMPLICIT NONE
    type(e_general) e
    integer mf

    if(e%kind==1) then
      write(mf,*) " kind 1: x and y angle "
      write(mf,*) e%t1%a
    endif

    if(e%kind==2) then
      write(mf,*) " kind 2 : y angle "
      write(mf,*) e%t2%a
    endif

    if(e%kind==3) then
      write(mf,*) " kind 3 : dx,dy,dz  "
      write(mf,*) e%t3%d
      write(mf,*) " coeff of (1+delta)  "
      write(mf,*) e%t3%dl
      if(e%t3%SIXTRACK ) then
       write(mf,*) " L_DESIGN, DL_SIXTRACK  "
       write(mf,*) e%t3%L_DESIGN , e%t3%DL_SIXTRACK
      else
       write(mf,*) " L_DESIGN  "
       write(mf,*) e%t3%L_DESIGN
      endif
    endif

  end subroutine print_e_general


subroutine zero_r_xy(t)
 implicit none
 TYPE(R_XY) t
 t%a=0
end subroutine zero_r_xy

subroutine zero_r_z(t)
 implicit none
 TYPE(R_Z) t
 t%a=0
end subroutine zero_r_z

subroutine zero_T_XYZ(t)
 implicit none
 TYPE(T_XYZ) t
     t%SIXTRACK=.false.
      t%L_DESIGN=0
      t%DL_SIXTRACK=0
      t%DL=0
      t%D=0
end subroutine zero_T_XYZ


subroutine zero_E_GENERAL(t,i)
 implicit none
 TYPE(E_GENERAL) t
 integer i

  t%kind=i
   call init(t%t1)
   call init(t%t2)
   call init(t%t3)

end subroutine zero_E_GENERAL

subroutine zero_E_GENERAL_s(t)
 implicit none
 TYPE(E_GENERAL) t(:)
 integer i

do i=1,size(t)
   call init(t(i),0)
enddo

end subroutine zero_E_GENERAL_s


  SUBROUTINE TRANS_dl(A,dl,LD,X,b,ctime,DL_SIXTRACK,SIXTRACK)   ! for sixtrack recombination
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) PZ
    real(dp),INTENT(IN):: A(3),b,dl,LD,DL_SIXTRACK
    LOGICAL(lp),INTENT(IN):: ctime,SIXTRACK

    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    if(ctime) then     ! THIS IS SIXTRACK HERE
       PZ=sqrt(1.0_dp+2.0_dp*X(5)/b+x(5)**2)
       X(1)=X(1)+A(3)*X(2)/pz
       X(3)=X(3)+A(3)*X(4)/pz
       IF(SIXTRACK) THEN
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/b+x(5))*A(3)/pz-A(3)/B   ! SIXTRACK DRIFT  b=beta0
          X(6)=X(6)+ dl*(1.0_dp/b+x(5))/pz + DL_SIXTRACK/B     ! EXTRA
       ELSE
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/b+x(5))*A(3)/pz+ dl*(1.0_dp/b+x(5))/pz-LD/B
       ENDIF
    else
       X(1)=X(1)+A(3)*X(2)/(1.0_dp+X(5))
       X(3)=X(3)+A(3)*X(4)/(1.0_dp+X(5))
       IF(SIXTRACK) THEN
          X(6)=X(6)+(A(3)/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
          X(6)=X(6)+dl + DL_SIXTRACK
       ELSE
          X(6)=X(6)+(A(3)/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))+a(3)+ dl-LD
       ENDIF
    endif

  END SUBROUTINE TRANS_dl

  SUBROUTINE track_e_general_s(E,X,b,ctime)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: ctime
    TYPE(E_GENERAL),INTENT(IN):: E(:)
    integer i

    do i=1,size(e)
       call track_e_general(E(i),X,b,ctime)
    enddo

  END SUBROUTINE track_e_general_s


  SUBROUTINE track_e_general(E,X,b,ctime)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: ctime
    TYPE(E_GENERAL),INTENT(IN):: E

    IF(E%KIND==1) THEN
       CALL ROT_YZ(E%T1%A(2),X,b,my_FALSE,ctime)  ! inverted
       CALL ROT_XZ(E%T1%A(1),X,b,my_FALSE,ctime)
    ELSEIF(E%KIND==2) THEN
       CALL ROT_XY(E%T2%A,X)
    ELSE
       CALL TRANS_dl(E%T3%D,E%T3%dl,E%T3%L_DESIGN,X,b,ctime,E%T3%DL_SIXTRACK,E%T3%SIXTRACK)
    ENDIF

  END SUBROUTINE track_e_general

  SUBROUTINE recombine(S,E)
    IMPLICIT NONE
    TYPE(E_GENERAL),INTENT(INout):: S(:)
    TYPE(E_GENERAL),INTENT(INout)::E(3)
    integer i,j,K

    do j=1,size(s)
       do i=1,size(s)
          if(s(i)%kind==3) then
             do k=i,size(s)-1
                call commute_e(s(k),s(k+1))
             enddo
          endif
       enddo
    enddo

    do j=1,size(s)
       do i=1,size(s)
          if(s(i)%kind==2) then
             do k=i,size(s)-1
                if(s(k+1)%kind==3) exit
                call commute_e(s(k),s(k+1))
             enddo
          endif
       enddo
    enddo

    do i=1,size(s)
       write(6,*) s(i)%kind
    enddo

    E(1)%kind=1
    E(2)%kind=2
    E(3)%kind=3
    E(1)%T1%A=0.0_dp
    E(2)%T2%A=0.0_dp
    E(3)%T3%D=0.0_dp
    E(3)%T3%DL=0.0_dp
    E(3)%T3%DL_SIXTRACK=0.0_dp
    E(3)%T3%L_DESIGN=0.0_dp
    do i=1,size(s)
       if(S(I)%KIND==1) THEN
          E(1)%T1%A=E(1)%T1%A+S(I)%T1%A
       ENDIF
       if(S(I)%KIND==2) THEN
          E(2)%T2%A=E(2)%T2%A+S(I)%T2%A
       ENDIF
       if(S(I)%KIND==3) THEN
          E(3)%T3%D=E(3)%T3%D+S(I)%T3%D
          E(3)%T3%DL=E(3)%T3%DL+S(I)%T3%DL
          E(3)%T3%L_DESIGN=E(3)%T3%L_DESIGN+S(I)%T3%L_DESIGN
       ENDIF
    enddo
    E(3)%T3%DL_SIXTRACK=-E(3)%T3%L_DESIGN+E(3)%T3%D(3)



  END   SUBROUTINE recombine




  SUBROUTINE COMMUTE_E(E1,E2)
    IMPLICIT NONE
    TYPE(E_GENERAL),INTENT(INOUT):: E1,E2
    TYPE(E_GENERAL) ET
    real(dp) D(2),A

    IF(E1%KIND==1) THEN
       IF(E2%KIND==1) THEN
          ET=E1
          E1=E2
          E2=ET
       ELSEIF(E2%KIND==2) THEN
          ET=E1
          A=E2%T2%A
          D(1)=COS(A)*E1%T1%A(1)+SIN(A)*E1%T1%A(2)
          D(2)=COS(A)*E1%T1%A(2)-SIN(A)*E1%T1%A(1)
          E1=E2
          E2=ET
          E2%T1%A=D
       ELSEIF(E2%KIND==3) THEN
          ET=E1
          E2%T3%DL=E2%T3%DL+E2%T3%D(1)*E1%T1%A(1)+E2%T3%D(2)*E1%T1%A(2)
          E2%T3%DL=E2%T3%DL-E2%T3%D(3)*(E1%T1%A(1)**2+E1%T1%A(2)**2)/2.D0
          E2%T3%D(1:2)=E2%T3%D(1:2)-E2%T3%D(3)*E1%T1%A
          E1=E2
          E2=ET
       ENDIF
       return
    ENDIF

    IF(E1%KIND==2) THEN
       IF(E2%KIND==1) THEN
          ET=E2
          A=E1%T2%A
          D(1)=COS(A)*E2%T1%A(1)-SIN(A)*E2%T1%A(2)
          D(2)=COS(A)*E2%T1%A(2)+SIN(A)*E2%T1%A(1)
          E2=E1
          E1=ET
          E1%T1%A=D
       ELSEIF(E2%KIND==2) THEN
          ET=E1
          E1=E2
          E2=ET
       ELSEIF(E2%KIND==3) THEN
          ET=E2
          A=E1%T2%A
          D(1)=COS(A)*E2%T3%D(1)-SIN(A)*E2%T3%D(2)
          D(2)=COS(A)*E2%T3%D(2)+SIN(A)*E2%T3%D(1)
          ET%T3%D(1:2)=D
          E2=E1
          E1=ET
       ENDIF
       return
    ENDIF

    IF(E1%KIND==3) THEN
       IF(E2%KIND==1) THEN
          ET=E2
          E1%T3%DL=E1%T3%DL-E1%T3%D(1)*E2%T1%A(1)-E1%T3%D(2)*E2%T1%A(2)
          E1%T3%DL=E1%T3%DL-E1%T3%D(3)*(E2%T1%A(1)**2+E2%T1%A(2)**2)/2.D0
          E1%T3%D(1:2)=E1%T3%D(1:2)+E1%T3%D(3)*E2%T1%A
          E2=E1
          E1=ET
       ELSEIF(E2%KIND==2) THEN
          ET=E1
          A=E2%T2%A
          D(1)=COS(A)*E1%T3%D(1)+SIN(A)*E1%T3%D(2)
          D(2)=COS(A)*E1%T3%D(2)-SIN(A)*E1%T3%D(1)
          ET%T3%D(1:2)=D
          E1=E2
          E2=ET
       ELSEIF(E2%KIND==3) THEN
          ET=E1
          E1=E2
          E2=ET
       ENDIF
       return
    ENDIF
  END SUBROUTINE COMMUTE_E

  SUBROUTINE ROT_YZR(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) XN(6)
    real(dp),INTENT(IN):: A,b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    XN(1)=X(3)
    XN(2)=X(4)
    XN(3)=X(1)      ! Using a relabelling i.e. a symmetry with respect to plane z
    XN(4)=X(2)      ! Could have used an x-y rotation of angle 90 degrees (AIMIN)
    XN(5)=X(5)
    XN(6)=X(6)
    CALL ROT_XZ(A,XN,b,EXACT,ctime)
    X(1)=XN(3)
    X(2)=XN(4)
    X(3)=XN(1)
    X(4)=XN(2)
    X(5)=XN(5)
    X(6)=XN(6)
  END SUBROUTINE ROT_YZR

  SUBROUTINE ROT_YZP(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) XN(6)
    real(dp),INTENT(IN):: A,b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    call PRTP("ROT_YZ:0", X)

    CALL ALLOC(XN,6)
    XN(1)=X(3)
    XN(2)=X(4)
    XN(3)=X(1)
    XN(4)=X(2)
    XN(5)=X(5)
    XN(6)=X(6)
    CALL ROT_XZ(A,XN,b,EXACT,ctime)
    X(1)=XN(3)
    X(2)=XN(4)
    X(3)=XN(1)
    X(4)=XN(2)
    X(5)=XN(5)
    X(6)=XN(6)
    CALL KILL(XN,6)

    call PRTP("ROT_YZ:1", X)

  END SUBROUTINE ROT_YZP

  SUBROUTINE TRANSR(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) PZ
    real(dp),INTENT(IN):: A(3),b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    IF(EXACT) THEN
       if(ctime) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(1.0_dp/b+X(5))/PZ
       else
          PZ=ROOT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(1.0_dp+X(5))/PZ
       endif
    ELSE
       if(ctime) then     ! THIS IS SIXTRACK HERE
          PZ=sqrt(1.0_dp+2.0_dp*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(1.0_dp+X(5))
          X(3)=X(3)+A(3)*X(4)/(1.0_dp+X(5))
          X(6)=X(6)+(A(3)/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))+a(3)
       endif
    ENDIF


  END SUBROUTINE TRANSR

  SUBROUTINE TRANSP(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) PZ
    real(dp),INTENT(IN):: A(3),b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    call PRTP("TRANS:0", X)

    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    IF(EXACT) THEN
       CALL ALLOC(PZ)
       if(ctime) then
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(1.0_dp/b+X(5))/PZ
       else
          PZ=SQRT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(1.0_dp+X(5))/PZ
       endif
       CALL KILL(PZ)
    ELSE
       if(ctime) then
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(1.0_dp+X(5))
          X(3)=X(3)+A(3)*X(4)/(1.0_dp+X(5))
          X(6)=X(6)+(A(3)/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))+a(3)
       endif
    ENDIF

    call PRTP("TRANS:1", X)

  END SUBROUTINE TRANSP



  SUBROUTINE ROT_XYR(A,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) XN(4)
    real(dp),INTENT(IN):: A
    real(dp)           :: cosa, sina
    cosa = COS(A)
    sina = SIN(A)

    !    IF(EXACT) THEN
    XN(1)=COSA*X(1)+SINA*X(3)
    XN(3)=COSA*X(3)-SINA*X(1)
    XN(2)=COSA*X(2)+SINA*X(4)
    XN(4)=COSA*X(4)-SINA*X(2)
    X(1)=XN(1)
    X(2)=XN(2)
    X(3)=XN(3)
    X(4)=XN(4)
    !    ELSE
    !       X(1)=X(1)+A*X(3)
    !       X(4)=X(4)-A*X(2)
    !       X(2)=X(2)+A*X(4)
    !       X(3)=X(3)-A*X(1)
    !    ENDIF

  END SUBROUTINE ROT_XYR

  SUBROUTINE ROT_XYP(A,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) XN(4)
    real(dp),INTENT(IN):: A
    real(dp)           :: cosa, sina

    call PRTP("ROT_XY:0", X)

    cosa = COS(A)
    sina = SIN(A)

    !    IF(EXACT) THEN
    CALL ALLOC(XN,4)

    XN(1)=COSA*X(1)+SINA*X(3)
    XN(3)=COSA*X(3)-SINA*X(1)
    XN(2)=COSA*X(2)+SINA*X(4)
    XN(4)=COSA*X(4)-SINA*X(2)
    X(1)=XN(1)
    X(2)=XN(2)
    X(3)=XN(3)
    X(4)=XN(4)
    CALL KILL(XN,4)


    !    ELSE
    !       X(1)=X(1)+A*X(3)
    !       X(4)=X(4)-A*X(2)
    !       X(2)=X(2)+A*X(4)
    !       X(3)=X(3)-A*X(1)
    !    ENDIF

    call PRTP("ROT_XY:1", X)

  END SUBROUTINE ROT_XYP


  SUBROUTINE ROT_XZR(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) XN(6),PZ,PT
    real(dp),INTENT(IN):: A,b
    real(dp) sina, cosa, tana
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    IF(EXACT) THEN
       COSA = COS(A)
       SINA = SIN(A)
       TANA = TAN(A)
       if(ctime) then
          PZ=ROOT(1.0_dp+2.0_dp*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=1.0_dp-X(2)*TANA/PZ
          XN(1)=X(1)/COSA/PT
          XN(2)=X(2)*COSA+SINA*PZ
          XN(3)=X(3)+X(4)*X(1)*TANA/PZ/PT
          XN(6)=X(6)+X(1)*TANA/PZ/PT*(1.0_dp/b+x(5))
       else
          PZ=ROOT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          PT=1.0_dp-X(2)*TANA/PZ
          XN(1)=X(1)/COSA/PT
          XN(2)=X(2)*COSA+SINA*PZ
          XN(3)=X(3)+X(4)*X(1)*TANA/PZ/PT
          XN(6)=X(6)+(1.0_dp+X(5))*X(1)*TANA/PZ/PT
       endif
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    ELSE
       if(ctime) then   ! SIXTRACK
          PZ=sqrt(1.0_dp+2.0_dp*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(1.0_dp/b+x(5))/PZ
       else
          X(2)=X(2)+A*(1.0_dp+X(5))
          X(6)=X(6)+A*X(1)
       endif
    ENDIF
    !    CALL CHECK_STABILITY(X)
  END SUBROUTINE ROT_XZR

  SUBROUTINE ROT_XZP(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) XN(6),PZ,PT
    real(dp),INTENT(IN):: A,b
    real(dp) sina, cosa, tana
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    call PRTP("ROT_XZ:0", X)

    IF(EXACT) THEN
       COSA = COS(A)
       SINA = SIN(A)
       TANA = TAN(A)

       CALL ALLOC(XN,6)
       CALL ALLOC(PZ)
       CALL ALLOC(PT)
       if(ctime) then
          PZ=SQRT(1.0_dp+2.0_dp*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=1.0_dp-X(2)*TANA/PZ
          XN(1)=X(1)/COSA/PT
          XN(2)=X(2)*COSA+SINA*PZ
          XN(3)=X(3)+X(4)*X(1)*TANA/PZ/PT
          XN(6)=X(6)+X(1)*TANA/PZ/PT*(1.0_dp/b+x(5))
       else
          PZ=SQRT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          PT=1.0_dp-X(2)*TANA/PZ
          XN(1)=X(1)/COSA/PT
          XN(2)=X(2)*COSA+SINA*PZ
          XN(3)=X(3)+X(4)*X(1)*TANA/PZ/PT
          XN(6)=X(6)+(1.0_dp+X(5))*X(1)*TANA/PZ/PT
       endif

       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)

       CALL KILL(XN,6)
       CALL KILL(PZ)
       CALL KILL(PT)
    ELSE

       if(ctime) then
          CALL ALLOC(PZ)
          PZ=SQRT(1.0_dp+2.0_dp*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(1.0_dp/b+x(5))/PZ
          CALL KILL(PZ)
       else
          X(2)=X(2)+A*(1.0_dp+X(5))
          X(6)=X(6)+A*X(1)
       endif


    ENDIF

    call PRTP("ROT_XZ:1", X)

  END SUBROUTINE ROT_XZP


end module S_euclidean
