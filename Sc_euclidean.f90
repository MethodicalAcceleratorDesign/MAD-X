!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file A_SCRATCH_SIZE.F90
module S_euclidean   ! and module rotation_mis
  use S_extend_poly
  implicit none

  PRIVATE       transr,transp,TRANSP_P,transs,TRANSS_S
  PRIVATE       ROT_XZR,ROT_XZP,ROT_XZP_P,ROT_XZS,ROT_XZS_S
  PRIVATE       ROT_YZR,ROT_YZP,ROT_YZP_P,ROT_YZS,ROT_YZS_S
  PRIVATE       ROT_XYR,ROT_XYP,ROT_XYP_P,ROT_XYS,ROT_XYS_S
  INTERFACE TRANS
     MODULE PROCEDURE transR
     MODULE PROCEDURE transp
     MODULE PROCEDURE TRANSP_P   ! New stuff
     MODULE PROCEDURE transS
     MODULE PROCEDURE TRANSS_S   ! New stuff
  END INTERFACE

  INTERFACE ROT_XZ
     MODULE PROCEDURE ROT_XZR
     MODULE PROCEDURE ROT_XZP
     MODULE PROCEDURE ROT_XZP_P
     MODULE PROCEDURE ROT_XZS
     MODULE PROCEDURE ROT_XZS_S
  END INTERFACE
  INTERFACE ROT_YZ
     MODULE PROCEDURE ROT_YZR
     MODULE PROCEDURE ROT_YZP
     MODULE PROCEDURE ROT_YZP_P
     MODULE PROCEDURE ROT_YZS
     MODULE PROCEDURE ROT_YZS_S
  END INTERFACE
  INTERFACE ROT_XY
     MODULE PROCEDURE ROT_XYR
     MODULE PROCEDURE ROT_XYP
     MODULE PROCEDURE ROT_XYP_P
     MODULE PROCEDURE ROT_XYS
     MODULE PROCEDURE ROT_XYS_S
  END INTERFACE


CONTAINS

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
  END SUBROUTINE ROT_YZP

  SUBROUTINE ROT_YZP_P(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8),INTENT(IN)::  A
    TYPE(REAL_8) XN(6)
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime
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
  END SUBROUTINE ROT_YZP_P

  SUBROUTINE ROT_YZS(A,Y,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8) XN(6)
    TYPE(REAL_8)  X(6)
    real(dp),INTENT(IN):: A,b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime
    CALL ALLOC(XN,6)
    CALL ALLOC(X,6)
    X=Y
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
    Y=X
    CALL KILL(X,6)
    CALL KILL(XN,6)
  END SUBROUTINE ROT_YZS

  SUBROUTINE ROT_YZS_S(A,Y,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8),INTENT(INOUT):: A
    TYPE(REAL_8) XN(6)
    TYPE(REAL_8)  X(6)
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime
    CALL ALLOC(XN,6)
    CALL ALLOC(X,6)
    X=Y
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
    Y=X
    CALL KILL(X,6)
    CALL KILL(XN,6)
  END SUBROUTINE ROT_YZS_S

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
          PZ=ROOT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one/b+X(5))/PZ
       else
          PZ=ROOT((one+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one+X(5))/PZ
       endif
    ELSE
       if(ctime) then
          PZ=ROOT(one+two*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/two/pz**2+one)*(one/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(one+X(5))
          X(3)=X(3)+A(3)*X(4)/(one+X(5))
          X(6)=X(6)+(A(3)/(one+X(5)))*(X(2)*X(2)+X(4)*X(4))/two/(one+X(5))+a(3)
       endif
    ENDIF


  END SUBROUTINE TRANSR

  SUBROUTINE TRANSP(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) PZ
    real(dp),INTENT(IN):: A(3),b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    IF(EXACT) THEN
       CALL ALLOC(PZ)
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one/b+X(5))/PZ
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one+X(5))/PZ
       endif
       CALL KILL(PZ)
    ELSE
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/two/pz**2+one)*(one/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(one+X(5))
          X(3)=X(3)+A(3)*X(4)/(one+X(5))
          X(6)=X(6)+(A(3)/(one+X(5)))*(X(2)*X(2)+X(4)*X(4))/two/(one+X(5))+a(3)
       endif
    ENDIF

  END SUBROUTINE TRANSP

  SUBROUTINE TRANSP_P(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8),INTENT(IN):: A(3)
    TYPE(REAL_8) PZ
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    IF(EXACT) THEN
       CALL ALLOC(PZ)
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one/b+X(5))/PZ
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one+X(5))/PZ
       endif
       CALL KILL(PZ)
    ELSE
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/two/pz**2+one)*(one/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(one+X(5))
          X(3)=X(3)+A(3)*X(4)/(one+X(5))
          X(6)=X(6)+(A(3)/(one+X(5)))*(X(2)*X(2)+X(4)*X(4))/two/(one+X(5))+a(3)
       endif
    ENDIF

  END SUBROUTINE TRANSP_P

  SUBROUTINE TRANSS(A,Y,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8) X(6)
    TYPE(REAL_8) PZ
    real(dp),INTENT(IN):: A(3),b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime


    CALL ALLOC(X,6)
    X=Y
    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    IF(EXACT) THEN
       CALL ALLOC(PZ)
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one/b+X(5))/PZ
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one+X(5))/PZ
       endif
       CALL KILL(PZ)
    ELSE
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/two/pz**2+one)*(one/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(one+X(5))
          X(3)=X(3)+A(3)*X(4)/(one+X(5))
          X(6)=X(6)+(A(3)/(one+X(5)))*(X(2)*X(2)+X(4)*X(4))/two/(one+X(5))+a(3)
       endif
    ENDIF

    Y=X
    CALL KILL(X,6)

  END SUBROUTINE TRANSS

  SUBROUTINE TRANSS_S(A,Y,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8),INTENT(IN):: A(3)
    TYPE(REAL_8) X(6)
    TYPE(REAL_8) PZ
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime


    CALL ALLOC(X,6)
    X=Y
    X(1)=X(1)-A(1)
    X(3)=X(3)-A(2)
    IF(EXACT) THEN
       CALL ALLOC(PZ)
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one/b+X(5))/PZ
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          X(1)=X(1)+A(3)*X(2)/PZ
          X(3)=X(3)+A(3)*X(4)/PZ
          X(6)=X(6)+A(3)*(one+X(5))/PZ
       endif
       CALL KILL(PZ)
    ELSE
       if(ctime) then
          PZ=SQRT(one+two*X(5)/b+x(5)**2)
          X(1)=X(1)+A(3)*X(2)/pz
          X(3)=X(3)+A(3)*X(4)/pz
          X(6)=X(6)+((X(2)*X(2)+X(4)*X(4))/two/pz**2+one)*(one/b+x(5))*A(3)/pz
       else
          X(1)=X(1)+A(3)*X(2)/(one+X(5))
          X(3)=X(3)+A(3)*X(4)/(one+X(5))
          X(6)=X(6)+(A(3)/(one+X(5)))*(X(2)*X(2)+X(4)*X(4))/two/(one+X(5))+a(3)
       endif
    ENDIF

    Y=X
    CALL KILL(X,6)

  END SUBROUTINE TRANSS_S

  SUBROUTINE ROT_XYR(A,X,EXACT)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) XN(4)
    real(dp),INTENT(IN):: A
    LOGICAL(lp),INTENT(IN):: EXACT

    IF(EXACT) THEN
       XN(1)=COS(A)*X(1)+SIN(A)*X(3)
       XN(3)=COS(A)*X(3)-SIN(A)*X(1)
       XN(2)=COS(A)*X(2)+SIN(A)*X(4)
       XN(4)=COS(A)*X(4)-SIN(A)*X(2)
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(4)=XN(4)
    ELSE
       X(1)=X(1)+A*X(3)
       X(4)=X(4)-A*X(2)
       X(2)=X(2)+A*X(4)
       X(3)=X(3)-A*X(1)
    ENDIF

  END SUBROUTINE ROT_XYR

  SUBROUTINE ROT_XYP(A,X,EXACT)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) XN(4)
    real(dp),INTENT(IN):: A
    LOGICAL(lp),INTENT(IN):: EXACT

    IF(EXACT) THEN
       CALL ALLOC(XN,4)
       XN(1)=COS(A)*X(1)+SIN(A)*X(3)
       XN(3)=COS(A)*X(3)-SIN(A)*X(1)
       XN(2)=COS(A)*X(2)+SIN(A)*X(4)
       XN(4)=COS(A)*X(4)-SIN(A)*X(2)
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(4)=XN(4)
       CALL KILL(XN,4)
    ELSE
       X(1)=X(1)+A*X(3)
       X(4)=X(4)-A*X(2)
       X(2)=X(2)+A*X(4)
       X(3)=X(3)-A*X(1)
    ENDIF
  END SUBROUTINE ROT_XYP

  SUBROUTINE ROT_XYP_P(A,X,EXACT)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8),INTENT(IN):: A
    TYPE(REAL_8) XN(4)
    LOGICAL(lp),INTENT(IN):: EXACT

    IF(EXACT) THEN
       CALL ALLOC(XN,4)
       XN(1)=COS(A)*X(1)+SIN(A)*X(3)
       XN(3)=COS(A)*X(3)-SIN(A)*X(1)
       XN(2)=COS(A)*X(2)+SIN(A)*X(4)
       XN(4)=COS(A)*X(4)-SIN(A)*X(2)
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(4)=XN(4)
       CALL KILL(XN,4)
    ELSE
       X(1)=X(1)+A*X(3)
       X(4)=X(4)-A*X(2)
       X(2)=X(2)+A*X(4)
       X(3)=X(3)-A*X(1)
    ENDIF
  END SUBROUTINE ROT_XYP_P

  SUBROUTINE ROT_XYS(A,Y,EXACT)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8)  X(6)
    TYPE(REAL_8) XN(4)
    real(dp),INTENT(IN):: A
    LOGICAL(lp),INTENT(IN):: EXACT

    CALL ALLOC(X,6)
    X=Y
    IF(EXACT) THEN
       CALL ALLOC(XN,4)
       XN(1)=COS(A)*X(1)+SIN(A)*X(3)
       XN(3)=COS(A)*X(3)-SIN(A)*X(1)
       XN(2)=COS(A)*X(2)+SIN(A)*X(4)
       XN(4)=COS(A)*X(4)-SIN(A)*X(2)

       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(4)=XN(4)
       CALL KILL(XN,4)
    ELSE
       X(1)=X(1)+A*X(3)
       X(4)=X(4)-A*X(2)
       X(2)=X(2)+A*X(4)
       X(3)=X(3)-A*X(1)
    ENDIF
    Y=X
    CALL KILL(X,6)

  END SUBROUTINE ROT_XYS

  SUBROUTINE ROT_XYS_S(A,Y,EXACT)
    IMPLICIT NONE
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8),INTENT(IN):: A
    TYPE(REAL_8)  X(6)
    TYPE(REAL_8) XN(4)
    LOGICAL(lp),INTENT(IN):: EXACT

    CALL ALLOC(X,6)
    X=Y
    IF(EXACT) THEN
       CALL ALLOC(XN,4)
       XN(1)=COS(A)*X(1)+SIN(A)*X(3)
       XN(3)=COS(A)*X(3)-SIN(A)*X(1)
       XN(2)=COS(A)*X(2)+SIN(A)*X(4)
       XN(4)=COS(A)*X(4)-SIN(A)*X(2)

       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(4)=XN(4)
       CALL KILL(XN,4)
    ELSE
       X(1)=X(1)+A*X(3)
       X(4)=X(4)-A*X(2)
       X(2)=X(2)+A*X(4)
       X(3)=X(3)-A*X(1)
    ENDIF
    Y=X
    CALL KILL(X,6)

  END SUBROUTINE ROT_XYS_S

  SUBROUTINE ROT_XZR(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp) XN(6),PZ,PT
    real(dp),INTENT(IN):: A,b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    IF(EXACT) THEN
       if(ctime) then
          PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+X(1)*TAN(A)/PZ/PT*(one/b+x(5))
       else
          PZ=ROOT((one+X(5))**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+(one+X(5))*X(1)*TAN(A)/PZ/PT
       endif
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    ELSE
       if(ctime) then
          PZ=ROOT(one+two*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(one/b+x(5))/PZ
       else
          X(2)=X(2)+A*(one+X(5))
          X(6)=X(6)+A*X(1)
       endif
    ENDIF

  END SUBROUTINE ROT_XZR

  SUBROUTINE ROT_XZP(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8) XN(6),PZ,PT
    real(dp),INTENT(IN):: A,b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    IF(EXACT) THEN
       CALL ALLOC(XN,6)
       CALL ALLOC(PZ)
       CALL ALLOC(PT)
       if(ctime) then
          PZ=SQRT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+X(1)*TAN(A)/PZ/PT*(one/b+x(5))
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+(one+X(5))*X(1)*TAN(A)/PZ/PT
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
          PZ=SQRT(one+two*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(one/b+x(5))/PZ
          CALL KILL(PZ)
       else
          X(2)=X(2)+A*(one+X(5))
          X(6)=X(6)+A*X(1)
       endif
    ENDIF

  END SUBROUTINE ROT_XZP

  SUBROUTINE ROT_XZP_P(A,X,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8),INTENT(IN):: A
    TYPE(REAL_8) XN(6),PZ,PT
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    IF(EXACT) THEN
       CALL ALLOC(XN,6)
       CALL ALLOC(PZ)
       CALL ALLOC(PT)
       if(ctime) then
          PZ=SQRT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+X(1)*TAN(A)/PZ/PT*(one/b+x(5))
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+(one+X(5))*X(1)*TAN(A)/PZ/PT
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
          PZ=SQRT(one+two*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(one/b+x(5))/PZ
          CALL KILL(PZ)
       else
          X(2)=X(2)+A*(one+X(5))
          X(6)=X(6)+A*X(1)
       endif
    ENDIF

  END SUBROUTINE ROT_XZP_P

  SUBROUTINE ROT_XZS(A,Y,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8)  X(6)
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8) XN(6),PZ,PT
    real(dp),INTENT(IN):: A,b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime
    CALL ALLOC(X,6)
    X=Y

    IF(EXACT) THEN
       CALL ALLOC(XN,6)
       CALL ALLOC(PZ)
       CALL ALLOC(PT)
       if(ctime) then
          PZ=SQRT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+X(1)*TAN(A)/PZ/PT*(one/b+x(5))
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+(one+X(5))*X(1)*TAN(A)/PZ/PT
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
          PZ=SQRT(one+two*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(one/b+x(5))/PZ
          CALL KILL(PZ)
       else
          X(2)=X(2)+A*(one+X(5))
          X(6)=X(6)+A*X(1)
       endif
    ENDIF
    Y=X
    CALL KILL(X,6)

  END SUBROUTINE ROT_XZS

  SUBROUTINE ROT_XZS_S(A,Y,b,EXACT,ctime)
    IMPLICIT NONE
    TYPE(REAL_8)  X(6)
    TYPE(ENV_8),INTENT(INOUT):: Y(6)
    TYPE(REAL_8),INTENT(IN):: A
    TYPE(REAL_8) XN(6),PZ,PT
    real(dp),INTENT(IN):: b
    LOGICAL(lp),INTENT(IN):: EXACT,ctime

    CALL ALLOC(X,6)
    X=Y

    IF(EXACT) THEN
       CALL ALLOC(XN,6)
       CALL ALLOC(PZ)
       CALL ALLOC(PT)
       if(ctime) then
          PZ=SQRT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+X(1)*TAN(A)/PZ/PT*(one/b+x(5))
       else
          PZ=SQRT((one+X(5))**2-X(2)**2-X(4)**2)
          PT=one-X(2)*TAN(A)/PZ
          XN(1)=X(1)/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
          XN(3)=X(3)+X(4)*X(1)*TAN(A)/PZ/PT
          XN(6)=X(6)+(one+X(5))*X(1)*TAN(A)/PZ/PT
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
          PZ=SQRT(one+two*x(5)/b+X(5)**2)
          X(2)=X(2)+A*PZ
          X(6)=X(6)+A*X(1)*(one/b+x(5))/PZ
          CALL KILL(PZ)
       else
          X(2)=X(2)+A*(one+X(5))
          X(6)=X(6)+A*X(1)
       endif
    ENDIF
    Y=X
    CALL KILL(X,6)

  END SUBROUTINE ROT_XZS_S

end module S_euclidean
