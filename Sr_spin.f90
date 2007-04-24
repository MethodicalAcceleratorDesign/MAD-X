!   B~= (H B_x, H B_y,   B_z)
!   A~= (  A_x,   A_y, H A_z)
!   H = 1 + h x / RHO_0
!
module ptc_spin
  !use orbit_ptc
  !use beam_beam_ptc
  use accel_ptc
  implicit none
  PRIVATE get_fieldR !,get_field
  PRIVATE get_BfieldR,get_BfieldP,get_Bfield,get_fieldp
  private GETMULB_TEAPOTR,GETMULB_TEAPOTP,get_Bfieldt
  private B_PANCAkEr,B_PANCAkEp,B_PANCAkE
  PRIVATE DIRECTION_VR,DIRECTION_VP,DIRECTION_V
  PRIVATE  B_PARA_PERP_r,B_PARA_PERP_p,B_PARA_PERP
  PRIVATE get_omegaR,get_omegaP,get_omega
  PRIVATE PUSH_SPINR,PUSH_SPINP,PUSH_SPIN
  PRIVATE TRACK_FRINGE_spin_R,TRACK_FRINGE_spin_P,TRACK_FRINGE_spin
  PRIVATE TRACK_NODE_LAYOUT_FLAG_spin_R ,TRACK_NODE_LAYOUT_FLAG_spin_P
  PRIVATE GET_BE_CAVR,GET_BE_CAVP ,GET_BE_CAV
  private rot_spin_x,rot_spin_xr,rot_spin_xp,rot_spin_z,rot_spin_zr,rot_spin_zp
  private rot_spin_yr,rot_spin_yp,rot_spin_y
  private PATCH_SPINR,PATCH_SPINP,PATCH_SPIN
  private MIS_SPINR,MIS_SPINP,MIS_SPIN
  private DTILTR_SPIN,DTILTP_SPIN,DTILT_SPIN
  PRIVATE TRACK_SPIN_FRONTR,TRACK_SPIN_FRONTP,TRACK_SPIN_FRONT
  PRIVATE TRACK_SPIN_BACKR,TRACK_SPIN_BACKP,TRACK_SPIN_BACK
  private PUSH_SPIN_RAY8,TRACK_SPIN_BACK_RAY8,TRACK_SPIN_FRONT_ray8,TRACK_FRINGE_spin_ray8
  private radiate_2p,radiate_2r,radiate_2

  REAL(DP) :: AG=A_ELECTRON

  INTERFACE radiate_2
     MODULE PROCEDURE radiate_2r
     MODULE PROCEDURE radiate_2p
  END INTERFACE

  INTERFACE PATCH_SPIN
     MODULE PROCEDURE PATCH_SPINR
     MODULE PROCEDURE PATCH_SPINP
  END INTERFACE

  INTERFACE MIS_SPIN
     MODULE PROCEDURE MIS_SPINR
     MODULE PROCEDURE MIS_SPINP
  END INTERFACE

  INTERFACE DTILT_SPIN
     MODULE PROCEDURE DTILTR_SPIN
     MODULE PROCEDURE DTILTP_SPIN
  END INTERFACE

  INTERFACE TRACK_SPIN_FRONT
     MODULE PROCEDURE TRACK_SPIN_FRONTR
     MODULE PROCEDURE TRACK_SPIN_FRONTP
     MODULE PROCEDURE TRACK_SPIN_FRONT_ray8
  END INTERFACE

  INTERFACE TRACK_SPIN_BACK
     MODULE PROCEDURE TRACK_SPIN_BACKR
     MODULE PROCEDURE TRACK_SPIN_BACKP
     MODULE PROCEDURE TRACK_SPIN_BACK_RAY8
  END INTERFACE

  INTERFACE GET_BE_CAV
     MODULE PROCEDURE GET_BE_CAVR
     MODULE PROCEDURE GET_BE_CAVP
  END INTERFACE

  INTERFACE rot_spin_x
     MODULE PROCEDURE rot_spin_xr
     MODULE PROCEDURE rot_spin_xp
  END INTERFACE

  INTERFACE rot_spin_y
     MODULE PROCEDURE rot_spin_yr
     MODULE PROCEDURE rot_spin_yp
  END INTERFACE

  INTERFACE rot_spin_z
     MODULE PROCEDURE rot_spin_zr
     MODULE PROCEDURE rot_spin_zp
  END INTERFACE

  INTERFACE TRACK_FRINGE_spin
     MODULE PROCEDURE TRACK_FRINGE_spin_R
     MODULE PROCEDURE TRACK_FRINGE_spin_P
     MODULE PROCEDURE TRACK_FRINGE_spin_ray8
  END INTERFACE

  INTERFACE TRACK_NODE_LAYOUT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spin_R
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spin_P
  END INTERFACE


  INTERFACE PUSH_SPIN
     MODULE PROCEDURE PUSH_SPINR
     MODULE PROCEDURE PUSH_SPINP
     MODULE PROCEDURE PUSH_SPIN_RAY8
  END INTERFACE

  INTERFACE get_omega
     MODULE PROCEDURE get_omegaR
     MODULE PROCEDURE get_omegaP
  END INTERFACE

  INTERFACE B_PARA_PERP
     MODULE PROCEDURE B_PARA_PERP_r
     MODULE PROCEDURE B_PARA_PERP_p
  END INTERFACE

  INTERFACE DIRECTION_V
     MODULE PROCEDURE DIRECTION_Vr
     MODULE PROCEDURE DIRECTION_Vp
  END INTERFACE

  INTERFACE get_field
     MODULE PROCEDURE get_fieldr   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE get_fieldp
  END INTERFACE

  INTERFACE get_Bfield
     MODULE PROCEDURE get_BfieldR   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE get_BfieldP   ! MID DEFINED AS 1/2 L
  END INTERFACE


  INTERFACE get_Bfieldt
     MODULE PROCEDURE GETMULB_TEAPOTR   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE GETMULB_TEAPOTP
  END INTERFACE

  INTERFACE B_PANCAkE
     MODULE PROCEDURE B_PANCAkEr   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE B_PANCAkEp
  END INTERFACE

contains

  subroutine rot_spin_yr(s,ang)
    implicit none
    REAL(DP),INTENT(INOUT) ::  S(3)
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st

    CO =COS(ang)
    SI =sin(ang)

    ST=  CO *S(1)+SI *S(3)
    S(3)=CO *S(3)-SI *S(1)
    S(1)=ST

  END subroutine rot_spin_yr

  subroutine rot_spin_Xr(s,ang)
    implicit none
    REAL(DP),INTENT(INOUT) ::  S(3)
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st

    CO =COS(ang)
    SI =sin(ang)

    ST=  CO *S(2)+SI *S(3)
    S(3)=CO *S(3)-SI *S(2)
    S(2)=ST

  END subroutine rot_spin_Xr

  subroutine rot_spin_zr(s,ang)
    implicit none
    REAL(DP),INTENT(INOUT) ::  S(3)
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st

    CO =COS(ang)
    SI =sin(ang)

    ST=  CO *S(1)+SI *S(2)
    S(2)=CO *S(2)-SI *S(1)
    S(1)=ST

  END subroutine rot_spin_zr


  subroutine rot_spin_yp(s,ang)
    implicit none
    type(real_8),INTENT(INOUT) ::  S(3)
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    !type(real_8) co,si,st

    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    ST=  CO *S(1)+SI *S(3)
    S(3)=CO *S(3)-SI *S(1)
    S(1)=ST

    call kill(st)

  END subroutine rot_spin_yp

  subroutine rot_spin_xp(s,ang)
    implicit none
    type(real_8),INTENT(INOUT) ::  S(3)
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st

    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    ST=  CO *S(2)+SI *S(3)
    S(3)=CO *S(3)-SI *S(2)
    S(2)=ST

    call kill(st)

  END subroutine rot_spin_xp

  subroutine rot_spin_zp(s,ang)
    implicit none
    type(real_8),INTENT(INOUT) ::  S(3)
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st

    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    ST=  CO *S(1)+SI *S(2)
    S(2)=CO *S(2)-SI *S(1)
    S(1)=ST

    call kill(st)

  END subroutine rot_spin_zp

  subroutine PUSH_SPIN_RAY8(EL,DS,FAC,S,before,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    TYPE(REAL_8), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8) OM(3),CO(3),SI(3)
    TYPE(REAL_8) ST
    type(probe_8),INTENT(INOUT) ::S
    integer i,j
    LOGICAL(LP),intent(in) :: before
    if(EL%kind<=kind1) return



    call PUSH_SPIN(EL,DS,FAC,S%S%x,S%X,S%E_IJ,before,POS)

  end subroutine PUSH_SPIN_RAY8

  subroutine radiate_2r(EL,DS,FAC,X,b2,dlds,XP,before,POS)
    implicit none
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    real(dp),INTENT(INOUT) :: X(6),XP(2)
    real(dp), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    real(dp), intent(in):: B2,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    real(dp)  st


    if(EL%P%TIME) then
       ST=root(one+two*X(5)/EL%P%beta0+x(5)**2)-one
    else
       ST=X(5)
    endif

    ! X(5)=X(5)+B2*FAC*DS
    !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    X(5)=one/(one/(one+X(5))+CRADF(EL%P)*(one+X(5))**3*B2*DLDS*FAC*DS)-one

    if(el%kind/=kindpa) then
       if(EL%P%TIME) then
          X(2)=X(2)*root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
          X(4)=X(4)*root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
       else
          X(2)=X(2)*(one+X(5))/(one+ST)
          X(4)=X(4)*(one+X(5))/(one+ST)
       endif
    endif

  end subroutine radiate_2r

  subroutine PUSH_SPINR(EL,DS,FAC,S,X,before,POS)
    implicit none
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    REAL(DP),INTENT(INOUT) :: X(6),S(3)
    REAL(DP), INTENT(IN) :: DS,FAC
    REAL(DP) OM(3),CO(3),SI(3),B2,XP(2)
    REAL(DP) ST,dlds
    LOGICAL(LP),intent(in) :: BEFORE

    if(EL%kind<=kind1) return
    CALL get_omega(EL,OM,B2,dlds,XP,X,POS)

    if(el%p%radiation_new.AND.BEFORE) then
       call radiate_2(EL,DS,FAC,X,b2,dlds,XP,before,POS)
    endif


    if(EL%P%SPIN) then
       CO(1)=COS(FAC*DS*OM(1)/TWO)
       SI(1)=SIN(FAC*DS*OM(1)/TWO)
       CO(2)=COS(FAC*DS*OM(2)/TWO)
       SI(2)=SIN(FAC*DS*OM(2)/TWO)
       CO(3)=COS(FAC*DS*OM(3))
       SI(3)=SIN(FAC*DS*OM(3))

       ST=   CO(1)*S(2)-SI(1)*S(3)
       S(3)= CO(1)*S(3)+SI(1)*S(2)
       S(2)=ST
       ST=  CO(2)*S(1)+SI(2)*S(3)
       S(3)=CO(2)*S(3)-SI(2)*S(1)
       S(1)=ST
       ST=   CO(3)*S(1)-SI(3)*S(2)
       S(2)= CO(3)*S(2)+SI(3)*S(1)
       S(1)=ST
       ST=  CO(2)*S(1)+SI(2)*S(3)
       S(3)=CO(2)*S(3)-SI(2)*S(1)
       S(1)=ST
       ST=   CO(1)*S(2)-SI(1)*S(3)
       S(3)= CO(1)*S(3)+SI(1)*S(2)
       S(2)=ST
    endif

    if(el%p%radiation_new.AND.(.NOT.BEFORE)) then
       call radiate_2(EL,DS,FAC,X,b2,dlds,XP,before,POS)
    endif

  END subroutine PUSH_SPINR

  subroutine radiate_2p(EL,DS,FAC,X,E_IJ,b2,dlds,XP  ,before,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: X(6),XP(2)
    real(dp),INTENT(INOUT) :: E_IJ(6,6)
    TYPE(REAL_8), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8), intent(in):: B2,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    TYPE(REAL_8) st
    logical(lp) done_stoch
    real(dp) b30,x1,x3,denf,b3
    type(damap) xpmap
    integer i,j

    if(.not.before) then
       denf=(one+x(5))**5/SQRT((one+X(5))**2-Xp(1)**2-Xp(2)**2)
       b30=b2
       b30=b30**c_1_5
       b30=cflucf(el%p)*b30
       denf=denf*b30*FAC*DS

       call alloc(xpmap)
       xpmap%v(1)=x(1)
       xpmap%v(3)=x(3)
       xpmap%v(5)=x(5)
       xpmap%v(6)=x(6)
       xpmap%v(2)=xp(1)
       xpmap%v(4)=xp(2)
       xpmap=xpmap**(-1)
       do i=1,6
          do j=1,6
             X1=(xpmap%v(i)).sub.'000010'
             X3=(xpmap%v(j)).sub.'000010'
             E_IJ(i,j)=E_IJ(i,j)+denf*x1*x3
          enddo
       enddo
       call kill(xpmap)
    endif


    call alloc(st)

    if(EL%P%TIME) then
       ST=SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)-one
    else
       ST=X(5)
    endif

    ! X(5)=X(5)+B2*FAC*DS
    !   X(5)=one/(one/(one+X(5))-B2*FAC*DS)-one
    X(5)=one/(one/(one+X(5))+CRADF(EL%P)*(one+X(5))**3*B2*DLDS*FAC*DS)-one

    if(el%kind/=kindpa) then
       if(EL%P%TIME) then
          X(2)=X(2)*SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
          X(4)=X(4)*SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
       else
          X(2)=X(2)*(one+X(5))/(one+ST)
          X(4)=X(4)*(one+X(5))/(one+ST)
       endif
    endif

    call kill(st)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(before) then
       denf=(one+x(5))**5/SQRT((one+X(5))**2-Xp(1)**2-Xp(2)**2)
       b30=b2
       b30=b30**c_1_5
       b30=cflucf(el%p)*b30
       denf=denf*b30*FAC*DS

       call alloc(xpmap)
       xpmap%v(1)=x(1)
       xpmap%v(3)=x(3)
       xpmap%v(5)=x(5)
       xpmap%v(6)=x(6)
       xpmap%v(2)=xp(1)
       xpmap%v(4)=xp(2)
       xpmap=xpmap**(-1)
       do i=1,6
          do j=1,6
             X1=(xpmap%v(i)).sub.'000010'
             X3=(xpmap%v(j)).sub.'000010'
             E_IJ(i,j)=E_IJ(i,j)+denf*x1*x3
          enddo
       enddo
       call kill(xpmap)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine radiate_2p

  subroutine PUSH_SPINP(EL,DS,FAC,S,X,E_IJ,before,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: X(6),S(3)
    real(dp),INTENT(INOUT) :: E_IJ(6,6)
    TYPE(REAL_8), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8) OM(3),CO(3),SI(3),B2,XP(2)
    TYPE(REAL_8) ST,dlds
    LOGICAL(LP),intent(in) :: BEFORE

    if(EL%kind<=kind1) return

    CALL ALLOC(OM,3)
    CALL ALLOC(CO,3)
    CALL ALLOC(SI,3)
    CALL ALLOC(XP,2)
    CALL ALLOC(ST,B2,dlds)

    CALL get_omega(EL,OM,B2,dlds,XP,X,POS)
    if(el%p%radiation_new.AND.BEFORE) then
       call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,POS)
    endif

    if(EL%P%SPIN) then
       CO(1)=COS(FAC*DS*OM(1)/TWO)
       SI(1)=SIN(FAC*DS*OM(1)/TWO)
       CO(2)=COS(FAC*DS*OM(2)/TWO)
       SI(2)=SIN(FAC*DS*OM(2)/TWO)
       CO(3)=COS(FAC*DS*OM(3))
       SI(3)=SIN(FAC*DS*OM(3))

       ST=   CO(1)*S(2)-SI(1)*S(3)
       S(3)= CO(1)*S(3)+SI(1)*S(2)
       S(2)=ST
       ST=  CO(2)*S(1)+SI(2)*S(3)
       S(3)=CO(2)*S(3)-SI(2)*S(1)
       S(1)=ST
       ST=   CO(3)*S(1)-SI(3)*S(2)
       S(2)= CO(3)*S(2)+SI(3)*S(1)
       S(1)=ST
       ST=  CO(2)*S(1)+SI(2)*S(3)
       S(3)=CO(2)*S(3)-SI(2)*S(1)
       S(1)=ST
       ST=   CO(1)*S(2)-SI(1)*S(3)
       S(3)= CO(1)*S(3)+SI(1)*S(2)
       S(2)=ST
    endif
    if(el%p%radiation_new.AND.(.NOT.BEFORE)) then
       call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,POS)
    endif

    CALL KILL(OM,3)
    CALL KILL(CO,3)
    CALL KILL(SI,3)
    CALL KILL(XP,2)
    CALL KILL(ST,B2,dlds)

  END subroutine PUSH_SPINP

  subroutine get_omegar(EL,OM,B2,dlds,XP,X,POS)
    implicit none
    TYPE(ELEMENT), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    REAL(DP),INTENT(INOUT) :: X(6),OM(3),B2,XP(2),DLDS
    REAL(DP)  B(3),E(3),BPA(3),BPE(3),D1,D2,GAMMA,EB(3)
    REAL(DP) BETA0,GAMMA0I
    INTEGER I

    P=>EL%P
    ! DLDS IS  REALLY D(CT)/DS * (1/(ONE/BETA0+X(5)))
    OM(2)=ZERO
    EB=ZERO
    BPA=ZERO
    BPE=ZERO
    B=ZERO
    E=ZERO

    CALL get_field(EL,B,E,X,POS)

    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7) ! Straight for all practical purposes
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ENDIF

       OM(2)=P%b0
    case(KIND4) ! CAVITY
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP,E,EB)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ENDIF
    case(KIND16:kind17,KIND20)
       CALL B_PARA_PERP(EL,0,X,B,BPA,BPE,XP)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)
       ENDIF
    case(kind10)     ! TEAPOT real curvilinear
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)
       ENDIF
       OM(2)=P%b0
    case(KINDPA)     ! fitted field for real magnet
       STOP 123   !  PATCH PROBLEMS???? CONVERTING TO PX ????
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP)
       if(p%time) then
          beta0=p%beta0;GAMMA0I=p%GAMMA0I;
       else
          beta0=one;GAMMA0I=zero;
       endif
       d1=SQRT(x(2)**2+x(4)**2+(one+P%B0*x(1))**2)
       d2=one+two*x(5)/beta0+x(5)**2
       d2=gamma0I/beta0/d2
       DLDS=SQRT((ONE+d2**2))*d1/(ONE/BETA0+X(5))
       OM(2)=P%b0
    case(KINDWIGGLER)
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    CASE(KIND21)     ! travelling wave cavity
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    case default
       OM(1)=ZERO
       OM(2)=ZERO
       OM(3)=ZERO
    END SELECT

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( ONE/P%BETA0 + X(5) )

    OM(1)=-DLDS*( (ONE+AG*GAMMA)*BPE(1) + (ONE+AG)*BPA(1) )
    OM(2)=-DLDS*( (ONE+AG*GAMMA)*BPE(2) + (ONE+AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*( (ONE+AG*GAMMA)*BPE(3) + (ONE+AG)*BPA(3) )
    DO I=1,3
       OM(I)=OM(I)-DLDS*(AG*GAMMA+GAMMA/(ONE+GAMMA))*EB(I)
    ENDDO
    if(P%RADIATION_NEW) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF

  end subroutine get_omegar

  subroutine get_omegaP(EL,OM,B2,dlds,XP,X,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    TYPE(REAL_8), INTENT(INOUT) :: X(6),OM(3),B2,XP(2)
    TYPE(REAL_8)  B(3),E(3),BPA(3),BPE(3),DLDS,D1,D2,GAMMA,EB(3)
    REAL(DP) BETA0,GAMMA0I
    INTEGER I

    CALL ALLOC(B,3)
    CALL ALLOC(E,3)
    CALL ALLOC(EB,3)
    CALL ALLOC(BPA,3)
    CALL ALLOC(BPE,3)
    CALL ALLOC(D1,D2,GAMMA)

    P=>EL%P
    ! DLDS IS  REALLY D(CT)/DS * (1/(ONE/BETA0+X(5)))
    OM(2)=ZERO

    CALL get_field(EL,B,E,X,POS)
    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7) ! Straight for all practical purposes
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ENDIF
       OM(2)=P%b0
    case(KIND4) ! CAVITY
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP,E,EB)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ENDIF
    case(KIND16:kind17,KIND20)
       CALL B_PARA_PERP(el,0,X,B,BPA,BPE,XP)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)
       ENDIF
    case(kind10)     ! TEAPOT real curvilinear
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP)
       IF(P%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-X(2)**2-X(4)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)
       ENDIF
       OM(2)=P%b0
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PARA_PERP(EL,1,X,B,BPA,BPE,XP)
       if(p%time) then
          beta0=p%beta0;GAMMA0I=p%GAMMA0I;
       else
          beta0=one;GAMMA0I=zero;
       endif
       d1=SQRT(x(2)**2+x(4)**2+(one+P%B0*x(1))**2)
       d2=one+two*x(5)/beta0+x(5)**2
       d2=gamma0I/beta0/d2
       DLDS=SQRT((ONE+d2**2))*d1/(ONE/BETA0+X(5))
       OM(2)=P%b0
    case(KINDWIGGLER)
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "

    CASE(KIND21)     ! travelling wave cavity
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    case default
       OM(1)=ZERO
       OM(2)=ZERO
       OM(3)=ZERO
    END SELECT

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( ONE/P%BETA0 + X(5) )

    OM(1)=-DLDS*( (ONE+AG*GAMMA)*BPE(1) + (ONE+AG)*BPA(1) )
    OM(2)=-DLDS*( (ONE+AG*GAMMA)*BPE(2) + (ONE+AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*( (ONE+AG*GAMMA)*BPE(3) + (ONE+AG)*BPA(3) )
    DO I=1,3
       OM(I)=OM(I)-DLDS*(AG*GAMMA+GAMMA/(ONE+GAMMA))*EB(I)
    ENDDO
    if(P%RADIATION_NEW) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF
    CALL KILL(B,3)
    CALL KILL(E,3)
    CALL KILL(EB,3)
    CALL KILL(BPA,3)
    CALL KILL(BPE,3)
    CALL KILL(D1,D2,GAMMA)

  end subroutine get_omegaP

  subroutine get_fieldr(EL,B,E,X,POS)
    implicit none
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    REAL(DP),INTENT(INOUT) :: B(3),E(3)
    REAL(DP),INTENT(INOUT) :: X(6)
    INTEGER I
    B=ZERO
    E=ZERO

    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7,KIND16:kind17,KIND20) ! Straight for all practical purposes
       CALL get_Bfield(EL,B,X)
    case(kind10)     ! TEAPOT real curvilinear
       CALL get_Bfieldt(EL%TP10,B,X)
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PANCAkE(EL%PA,B,X,POS)
    case(KINDWIGGLER)
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    CASE(KIND4)      ! Pill box cavity
       CALL GET_BE_CAV(EL%C4,B,E,X)
    CASE(KIND21)     ! travelling wave cavity
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    case default
    END SELECT

    DO I=1,3
       B(I)=B(I)*EL%P%CHARGE
       E(I)=E(I)*EL%P%CHARGE
    ENDDO


  end subroutine get_fieldr

  subroutine get_fieldp(EL,B,E,X,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(INOUT) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: B(3),E(3)
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    INTEGER I

    DO I=1,3
       B(I)=ZERO
       E(I)=ZERO
    ENDDO

    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7,KIND16:kind17,KIND20) ! Straight for all practical purposes
       CALL get_Bfield(EL,B,X)
    case(kind10)     ! TEAPOT real curvilinear
       CALL get_Bfieldt(EL%TP10,B,X)
       !     WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PANCAkE(EL%PA,B,X,POS)
    case(KINDWIGGLER)
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    CASE(KIND4)      ! Pill box cavity
       CALL GET_BE_CAV(EL%C4,B,E,X)
    CASE(KIND21)     ! travelling wave cavity
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    case default
    END SELECT

    DO I=1,3
       B(I)=B(I)*EL%P%CHARGE
       E(I)=E(I)*EL%P%CHARGE
    ENDDO

  end subroutine get_fieldp

  SUBROUTINE get_BfieldR(EL,B,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENT),INTENT(IN):: EL
    real(dp) X1,X3,BBYTW,BBXTW,BBYTWT
    INTEGER J

    X1=X(1)
    X3=X(3)

    IF(EL%P%NMUL>=1) THEN
       BBYTW=EL%BN(EL%P%NMUL)
       BBXTW=EL%AN(EL%P%NMUL)


       DO  J=EL%P%NMUL-1,1,-1
          BBYTWT=X1*BBYTW-X3*BBXTW+EL%BN(J)
          BBXTW=X3*BBYTW+X1*BBXTW+EL%AN(J)
          BBYTW=BBYTWT
       ENDDO
    ELSE
       BBYTW=zero
       BBXTW=zero
    ENDIF
    B(1)=BBXTW;B(2)=BBYTW;
    IF(ASSOCIATED(EL%B_SOL)) THEN
       B(3)=EL%B_SOL;
    ELSE
       B(3)=ZERO
    ENDIF
  END SUBROUTINE get_BfieldR

  SUBROUTINE get_BfieldP(EL,B,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENTP),INTENT(IN):: EL
    TYPE(REAL_8)  X1,X3,BBYTW,BBXTW,BBYTWT
    INTEGER J

    CALL ALLOC(X1,X3,BBYTW,BBXTW,BBYTWT)
    X1=X(1)
    X3=X(3)

    IF(EL%P%NMUL>=1) THEN
       BBYTW=EL%BN(EL%P%NMUL)
       BBXTW=EL%AN(EL%P%NMUL)


       DO  J=EL%P%NMUL-1,1,-1
          BBYTWT=X1*BBYTW-X3*BBXTW+EL%BN(J)
          BBXTW=X3*BBYTW+X1*BBXTW+EL%AN(J)
          BBYTW=BBYTWT
       ENDDO
    ELSE
       BBYTW=zero
       BBXTW=zero
    ENDIF
    B(1)=BBXTW;B(2)=BBYTW;
    IF(ASSOCIATED(EL%B_SOL)) THEN
       B(3)=EL%B_SOL;
    ELSE
       B(3)=ZERO
    ENDIF
    CALL KILL(X1,X3,BBYTW,BBXTW,BBYTWT)
  END SUBROUTINE get_BfieldP

  SUBROUTINE GETMULB_TEAPOTR(EL,B,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3)
    TYPE(TEAPOT),INTENT(IN):: EL
    real(dp) X1,X3,BX,BY,BTX,BTY,BtYT
    INTEGER J,M,A,K

    X1=X(1)
    X3=X(3)

    BX=zero
    BY=zero

    k=0
    m=EL%P%nmul-1
    do a=m,1,-1
       BTX=zero
       BTY=zero
       do j=m-a,1,-1
          k=k+1
          !b%i(k)=a
          !b%j(k)=j
          BTX= (BTX+EL%BF_X(k))*X3  !x1
          BTY= (BTY+EL%BF_Y(k))*X3
       enddo

       k=k+1
       !  b%i(k)=a
       !  b%j(k)=0
       BTX= (BTX+EL%BF_X(k))
       BTY= (BTY+EL%BF_Y(k))
       BX= (BX+BTX)*X1
       BY= (BY+BTY)*X1
    enddo
    BTX=zero
    BTY=zero
    do j=m,1,-1
       k=k+1
       !  b%i(k)=0
       !  b%j(k)=j
       BTX= (BTX+EL%BF_X(k))*X3
       BTY= (BTY+EL%BF_Y(k))*X3
    enddo
    k=k+1
    !    b%i(k)=0
    !    b%j(k)=0
    BX= BX+BTX+EL%BF_X(k)  !+X3
    BY= BY+BTY+EL%BF_Y(k)  !+X3

    ! etienne
    IF(EL%P%NMUL>SECTOR_NMUL) THEN
       BtY=EL%BN(EL%P%NMUL)
       BtX=EL%AN(EL%P%NMUL)


       DO  J=EL%P%NMUL-1,SECTOR_NMUL+1,-1
          BtYT=X1*BtY-X3*BtX+EL%BN(J)
          BtX =X3*BtY+X1*BtX+EL%AN(J)
          BtY =BtYT
       ENDDO

       DO  J=SECTOR_NMUL, 1,-1
          BtYT=X1*BtY-X3*BtX
          BtX =X3*BtY+X1*BtX
          BtY =BtYT
       ENDDO

       BX= BX-BTX
       BY= BY+BTY

    ENDIF



    B(1)=BY/(one+EL%P%B0*X(1))
    B(2)=-BX/(one+EL%P%B0*X(1))
    B(3)=zero

  END SUBROUTINE GETMULB_TEAPOTR

  SUBROUTINE GETMULB_TEAPOTP(EL,B,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6),B(3)
    TYPE(TEAPOTP),INTENT(IN):: EL
    TYPE(REAL_8) X1,X3,BX,BY,BTX,BTY,BtYT
    INTEGER J,M,A,K

    CALL ALLOC(X1,X3,BX,BY,BTX,BTY,BtYT)

    X1=X(1)
    X3=X(3)

    BX=zero
    BY=zero

    k=0
    m=EL%P%nmul-1
    do a=m,1,-1
       BTX=zero
       BTY=zero
       do j=m-a,1,-1
          k=k+1
          !b%i(k)=a
          !b%j(k)=j
          BTX= (BTX+EL%BF_X(k))*X3  !x1
          BTY= (BTY+EL%BF_Y(k))*X3
       enddo

       k=k+1
       !  b%i(k)=a
       !  b%j(k)=0
       BTX= (BTX+EL%BF_X(k))
       BTY= (BTY+EL%BF_Y(k))
       BX= (BX+BTX)*X1
       BY= (BY+BTY)*X1
    enddo
    BTX=zero
    BTY=zero
    do j=m,1,-1
       k=k+1
       !  b%i(k)=0
       !  b%j(k)=j
       BTX= (BTX+EL%BF_X(k))*X3
       BTY= (BTY+EL%BF_Y(k))*X3
    enddo
    k=k+1
    !    b%i(k)=0
    !    b%j(k)=0
    BX= BX+BTX+EL%BF_X(k)  !+X3
    BY= BY+BTY+EL%BF_Y(k)  !+X3

    ! etienne
    IF(EL%P%NMUL>SECTOR_NMUL) THEN
       BtY=EL%BN(EL%P%NMUL)
       BtX=EL%AN(EL%P%NMUL)


       DO  J=EL%P%NMUL-1,SECTOR_NMUL+1,-1
          BtYT=X1*BtY-X3*BtX+EL%BN(J)
          BtX =X3*BtY+X1*BtX+EL%AN(J)
          BtY =BtYT
       ENDDO

       DO  J=SECTOR_NMUL, 1,-1
          BtYT=X1*BtY-X3*BtX
          BtX =X3*BtY+X1*BtX
          BtY =BtYT
       ENDDO

       BX= BX-BTX
       BY= BY+BTY

    ENDIF



    B(1)=BY/(one+EL%P%B0*X(1))
    B(2)=-BX/(one+EL%P%B0*X(1))
    B(3)=zero

    CALL KILL(X1,X3,BX,BY,BTX,BTY,BtYT)

  END SUBROUTINE GETMULB_TEAPOTP

  SUBROUTINE GET_BE_CAVR(EL,B,E,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3),E(3)
    TYPE(CAV4),INTENT(INOUT):: EL
    real(dp) DF,R2,F,DR2,O,VL
    INTEGER I

    E=ZERO
    B=ZERO
    IF(EL%P%NOCAVITY) RETURN

    O=twopi*EL%freq/CLIGHT
    VL=EL%volt*c_1d_3/EL%P%P0C

    DF=ZERO
    F=ONE
    R2=ONE

    DO I=1,EL%N_BESSEL
       R2=-R2*O**2/FOUR/(I+1)**2
       DR2=R2*I
       DF=DF+DR2*2
       R2=R2*(X(1)**2+X(3)**2)
       F=F+R2
    ENDDO

    !    EL%DELTA_E=x(5)

    IF(EL%N_BESSEL>0) THEN
       B(2)=-X(1)*DF*VL*COS(O*X(6)+EL%PHAS+phase0)/O
       B(1)= X(3)*DF*VL*COS(O*X(6)+EL%PHAS+phase0)/O
    ENDIF

    E(3)=-F*VL*SIN(O*x(6)+EL%PHAS+phase0)

  END SUBROUTINE GET_BE_CAVR

  SUBROUTINE GET_BE_CAVP(EL,B,E,X)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6),B(3),E(3)
    TYPE(CAV4p),INTENT(INOUT):: EL
    TYPE(REAL_8) DF,R2,F,DR2,O,VL
    INTEGER I

    DO I=1,3
       E(I)=ZERO
       B(I)=ZERO
    ENDDO

    IF(EL%P%NOCAVITY) RETURN
    CALL ALLOC(DF,R2,F,DR2,O,VL)

    O=twopi*EL%freq/CLIGHT
    VL=EL%volt*c_1d_3/EL%P%P0C

    DF=ZERO
    F=ONE
    R2=ONE

    DO I=1,EL%N_BESSEL
       R2=-R2*O**2/FOUR/(I+1)**2
       DR2=R2*I
       DF=DF+DR2*2
       R2=R2*(X(1)**2+X(3)**2)
       F=F+R2
    ENDDO

    !    EL%DELTA_E=x(5)

    IF(EL%N_BESSEL>0) THEN
       B(2)=-X(1)*DF*VL*COS(O*X(6)+EL%PHAS+phase0)/O
       B(1)= X(3)*DF*VL*COS(O*X(6)+EL%PHAS+phase0)/O
    ENDIF

    E(3)=-F*VL*SIN(O*x(6)+EL%PHAS+phase0)

    CALL KILL(DF,R2,F,DR2,O,VL)

  END SUBROUTINE GET_BE_CAVP

  subroutine B_PANCAkEr(EL,B,X,POS)
    IMPLICIT NONE
    real(dp), INTENT(INout) :: X(6)
    INTEGER, INTENT(INOUT) :: POS
    TYPE(PANCAKE),  INTENT(INOUT) :: EL
    real(dp) , INTENT(INOUT) ::B(3)

    B(1)=X(1);
    B(2)=X(3);
    B(3)=ZERO;

    CALL trackg(EL%B(POS),B)

    b(1)=EL%SCALE*b(1)
    b(2)=EL%SCALE*b(2)
    !    b(3)=EL%SCALE*el%p%charge*el%p%dir*b(3)
    b(3)=EL%SCALE*b(3)

  END subroutine B_PANCAkEr

  subroutine B_PANCAkEp(EL,B,X,POS)
    IMPLICIT NONE
    type(real_8), INTENT(INout) :: X(6)
    INTEGER, INTENT(INOUT) :: POS
    TYPE(PANCAKEP),  INTENT(INOUT) :: EL
    type(real_8), INTENT(INOUT) ::B(3)

    B(1)=X(1);
    B(2)=X(3);
    B(3)=ZERO;

    CALL trackg(EL%B(POS),B)

    b(1)=EL%SCALE*b(1)
    b(2)=EL%SCALE*b(2)
    !    b(3)=EL%SCALE*el%p%charge*el%p%dir*b(3)
    b(3)=EL%SCALE*b(3)

  END subroutine B_PANCAkEp

  subroutine B_PARA_PERP_R(EL,TEAPOT_LIKE,X,B,BPA,BPE,XP,EF,EFB)
    IMPLICIT NONE
    REAL(DP),  INTENT(INout) :: X(6)
    TYPE(ELEMENT),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    REAL(DP),  INTENT(INout) :: B(3),BPA(3),BPE(3),XP(2)
    REAL(DP),  OPTIONAL ::EF(3),EFB(3)
    INTEGER TEAPOT_LIKE,i
    REAL(DP) e(3),be

    P=>EL%P

    !  this routines gives us  B parallel and B perpendicular
    ! Also if EF is present, E perpendicular times beta is return

    call DIRECTION_V(EL,TEAPOT_LIKE,X,E,XP)

    be=b(1)*e(1)+b(2)*e(2)+b(3)*e(3)

    do i=1,3
       BPA(i)=be*e(i)
    enddo
    do i=1,3
       BPE(i)=B(i)-BPA(i)
    enddo

    IF(PRESENT(EF)) THEN
       BE=ROOT(ONE+TWO*X(5)+X(5)**2)/(ONE/P%BETA0+X(5)) !  BUGGGGGGG CHECK
       E=BE*E
       EFB(1)=EF(2)*E(3)-EF(3)*E(2)      ! BETA * E X e = Beta * E_perp
       EFB(2)=EF(3)*E(1)-EF(1)*E(3)
       EFB(3)=EF(1)*E(2)-EF(2)*E(1)
    ENDIF

  END subroutine B_PARA_PERP_R

  subroutine B_PARA_PERP_p(EL,TEAPOT_LIKE,X,B,BPA,BPE,XP,EF,EFB)
    IMPLICIT NONE
    type(real_8),  INTENT(INout) :: X(6)
    TYPE(ELEMENTP),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    type(real_8),  INTENT(INout) :: B(3),BPA(3),BPE(3),XP(2)
    type(real_8),  OPTIONAL ::EF(3),EFB(3)
    INTEGER TEAPOT_LIKE,i
    type(real_8) e(3),be

    P=>EL%P

    call alloc(e,3); call alloc(be);

    call DIRECTION_V(EL,TEAPOT_LIKE,X,E,XP)

    be=b(1)*e(1)+b(2)*e(2)+b(3)*e(3)

    do i=1,3
       BPA(i)=be*e(i)
    enddo
    do i=1,3
       BPE(i)=B(i)-BPA(i)
    enddo

    IF(PRESENT(EF)) THEN
       BE=SQRT(ONE+TWO*X(5)+X(5)**2)/(ONE/P%BETA0+X(5))
       E(1)=BE*E(1);E(2)=BE*E(2);E(3)=BE*E(3);
       EFB(1)=EF(2)*E(3)-EF(3)*E(2)
       EFB(2)=EF(3)*E(1)-EF(1)*E(3)
       EFB(3)=EF(1)*E(2)-EF(2)*E(1)
    ENDIF

    call kill(e,3); call kill(be);

  END subroutine B_PARA_PERP_p


  subroutine DIRECTION_VR(EL,TEAPOT_LIKE,X,E,XP)
    IMPLICIT NONE
    REAL(DP),  INTENT(INout) :: X(6),XP(2)
    TYPE(ELEMENT),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    REAL(DP),  INTENT(INOUT) ::E(3)
    REAL(DP) N,H,DP1
    INTEGER TEAPOT_LIKE

    P=>EL%P


    IF(EL%KIND/=KINDPA) THEN

       IF(ASSOCIATED(EL%B_SOL)) THEN  !SOLENOID
          IF(EL%P%TIME) THEN
             DP1=ROOT(ONE+TWO*X(5)/P%BETA0+X(5)**2)
          ELSE
             DP1=ONE+X(5)
          ENDIF

          N=ROOT(DP1**2-(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)**2-(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)**2)

          E(1)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)/DP1
          E(2)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)/DP1
          E(3)=N/DP1
       ELSE    ! NOT SOLENOID
          IF(EL%P%TIME) THEN
             DP1=ROOT(ONE+TWO*X(5)/P%BETA0+X(5)**2)
          ELSE
             DP1=ONE+X(5)
          ENDIF

          N=ROOT(DP1**2-X(2)**2-X(4)**2)

          E(1)=X(2)/DP1
          E(2)=X(4)/DP1
          E(3)=N/DP1
       ENDIF  ! END OF SOLENOID
       XP(1)=X(2)/N
       XP(2)=X(4)/N


    ELSE    ! NON CANONICAL VARIABLES
       H=ONE+P%B0*TEAPOT_LIKE*X(1)
       N=ROOT(H**2+X(2)**2+X(4)**2)
       E(1)=X(2)/N
       E(2)=X(4)/N
       E(3)=H/N
       XP(1)=X(2)
       XP(2)=X(4)
    ENDIF


  END subroutine DIRECTION_VR

  subroutine DIRECTION_VP(EL,TEAPOT_LIKE,X,E,XP)
    IMPLICIT NONE
    type(real_8), INTENT(INout) :: X(6),XP(2)
    TYPE(ELEMENTP),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    type(real_8), INTENT(INOUT) ::E(3)
    type(real_8) N,H,DP1
    INTEGER TEAPOT_LIKE

    P=>EL%P

    CALL ALLOC(N,H,DP1 )

    IF(EL%KIND/=KINDPA) THEN

       IF(ASSOCIATED(EL%B_SOL)) THEN  !SOLENOID
          IF(EL%P%TIME) THEN
             DP1=SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2)
          ELSE
             DP1=ONE+X(5)
          ENDIF

          N=SQRT(DP1**2-(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)**2-(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)**2)

          E(1)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)/DP1
          E(2)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)/DP1
          E(3)=N/DP1
       ELSE
          IF(EL%P%TIME) THEN
             DP1=SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2)
          ELSE
             DP1=ONE+X(5)
          ENDIF

          N=SQRT(DP1**2-X(2)**2-X(4)**2)

          E(1)=X(2)/DP1
          E(2)=X(4)/DP1
          E(3)=N/DP1
       ENDIF
       XP(1)=X(2)/N
       XP(2)=X(4)/N

    ELSE    ! NON CANONICAL VARIABLES
       H=ONE+P%B0*TEAPOT_LIKE*X(1)
       N=SQRT(H**2+X(2)**2+X(4)**2)
       E(1)=X(2)/N
       E(2)=X(4)/N
       E(3)=H/N
       XP(1)=X(2)
       XP(2)=X(4)

    ENDIF


    CALL KILL(N,H,DP1 )

  END subroutine DIRECTION_VP

!!!!!!!!!!!!   TRACKING ROUTINES    !!!!!!!!!!!!



  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spin_R(R,xs,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    type(probe), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C
    real(dp) fac,ds

    CALL RESET_APERTURE_FLAG

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       if(associated(c%bb)) call BBKICK(c%BB,XS%X)


       if(c%cas==0) then
          ds=c%parent_fibre%MAG%L/c%parent_fibre%MAG%p%nst
          fac=half
          if(k%radiation_new) then
          endif
          call PUSH_SPIN(c%parent_fibre%mag,ds,FAC,XS%S%X,XS%X,my_true)
          CALL TRACK_NODE_SINGLE(C,XS%X,K,R%CHARGE)
          call PUSH_SPIN(c%parent_fibre%mag,ds,FAC,XS%S%X,XS%X,my_false)
          if(k%radiation_new) then
          endif
       elseIF(c%cas==case1.or.c%cas==case2) then
          CALL TRACK_FRINGE_spin(C,XS%X,XS%S%X,K)
          CALL TRACK_NODE_SINGLE(C,XS%X,K,R%CHARGE)
       else
          IF(c%cas==caseP1) THEN
             CALL TRACK_SPIN_FRONT(C%PARENT_FIBRE,XS%S%X)
          ELSE
             CALL TRACK_SPIN_BACK(C%PARENT_FIBRE,XS%S%X)
          ENDIF
          CALL TRACK_NODE_SINGLE(C,XS%X,K,R%CHARGE)

       endif

       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spin_R


  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spin_P(R,XS,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    TYPE(probe_8), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C
    real(dp) fac
    TYPE(REAL_8)ds

    CALL ALLOC(DS)

    CALL RESET_APERTURE_FLAG

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       if(associated(c%bb)) call BBKICK(c%BB,XS%X)


       if(c%cas==0) then
          ds=c%parent_fibre%MAG%L/c%parent_fibre%MAG%p%nst
          fac=half

          call PUSH_SPIN(c%parent_fibre%magP,ds,FAC,XS,my_true)
          CALL TRACK_NODE_SINGLE(C,XS%X,K,R%CHARGE)
          call PUSH_SPIN(c%parent_fibre%magP,ds,FAC,XS,my_true)
       elseIF(c%cas==case1.or.c%cas==case2) then
          CALL TRACK_FRINGE_spin(C,XS,K)
          CALL TRACK_NODE_SINGLE(C,XS%X,K,R%CHARGE)
       else
          IF(c%cas==caseP1) THEN
             CALL TRACK_SPIN_FRONT(C%PARENT_FIBRE,XS)
          ELSE
             CALL TRACK_SPIN_BACK(C%PARENT_FIBRE,XS)
          ENDIF
          CALL TRACK_NODE_SINGLE(C,XS%X,K,R%CHARGE)
       endif




       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

    CALL KILL(DS)

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spin_P





  SUBROUTINE TRACK_FRINGE_spin_R(C,X,S,K)
    IMPLICIT NONE
    real(dp), INTENT(INOUT):: X(6),S(3)
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENT), POINTER :: EL

    el=>C%PARENT_FIBRE%MAG
    SELECT CASE(EL%KIND)
    CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
    case(KIND2)
    case(KIND4)
    case(KIND5)
    case(KIND6)
    case(KIND7)
    case(KIND10)
    case(KIND16)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%EDGE(1))
       ELSE
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%EDGE(2))
       ENDIF
    case(KIND20)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       ELSE
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       ENDIF

    case(KIND17)
    case(KIND21)
    case(KINDWIGGLER)
    case(KINDPA)
    CASE DEFAULT
       WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
       stop 666
    END SELECT


  END SUBROUTINE TRACK_FRINGE_spin_R

  SUBROUTINE TRACK_FRINGE_spin_P(C,X,S,K)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT):: X(6),S(3)
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENTP), POINTER :: EL

    el=>C%PARENT_FIBRE%MAGP
    SELECT CASE(EL%KIND)
    CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
    case(KIND2)
    case(KIND4)
    case(KIND5)
    case(KIND6)
    case(KIND7)
    case(KIND10)
    case(KIND16)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%EDGE(1))
       ELSE
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%EDGE(2))
       ENDIF
    case(KIND20)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       ELSE
          CALL rot_spin_y(s,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       ENDIF

    case(KIND17)
    case(KIND21)
    case(KINDWIGGLER)
    case(KINDPA)
    CASE DEFAULT
       WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
       stop 666
    END SELECT


  END SUBROUTINE TRACK_FRINGE_spin_P

  SUBROUTINE TRACK_FRINGE_spin_ray8(C,S,K)
    IMPLICIT NONE
    type(probe_8),INTENT(INOUT) ::S
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C

    integer i,j

    call TRACK_FRINGE_spin(C,S%X,S%S%X,K)



  END SUBROUTINE TRACK_FRINGE_spin_ray8


  SUBROUTINE TRACK_SPIN_FRONTR(C,S)
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    real(dp), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE

    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

    ! PUSHING BEAM
    !




    ! The chart frame of reference is located here implicitely
    IF((PATCHG==1).or.(PATCHG==3)) THEN
       CALL PATCH_SPIN(C,S,MY_TRUE)
    ENDIF


    CALL DTILT_SPIN(C%DIR,C%MAG%P%TILTD,1,S)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,S,DONEITT)
    ENDIF

  END SUBROUTINE TRACK_SPIN_FRONTR

  SUBROUTINE TRACK_SPIN_FRONTP(C,S)
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE

    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

    ! PUSHING BEAM
    !




    ! The chart frame of reference is located here implicitely
    IF((PATCHG==1).or.(PATCHG==3)) THEN
       CALL PATCH_SPIN(C,S,MY_TRUE)
    ENDIF


    CALL DTILT_SPIN(C%DIR,C%MAG%P%TILTD,1,S)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,S,DONEITT)
    ENDIF

  END SUBROUTINE TRACK_SPIN_FRONTP

  SUBROUTINE TRACK_SPIN_FRONT_RAY8(C,S)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    type(probe_8),INTENT(INOUT) ::S
    integer i,j
    type(real_8) sp(3)


    call TRACK_SPIN_front(C,S%S%X)


  end subroutine TRACK_SPIN_FRONT_RAY8


  ! back patch/misaglinments/tilt

  SUBROUTINE TRACK_SPIN_BACKR(C,S)
    implicit none
    logical(lp) :: doneitf=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    REAL(DP), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE



    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,S,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILT_SPIN(C%DIR,C%MAG%P%TILTD,2,S)


    IF((PATCHG==2).or.(PATCHG==3)) THEN
       CALL PATCH_SPIN(C,S,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely


  END SUBROUTINE TRACK_SPIN_BACKR

  SUBROUTINE TRACK_SPIN_BACK_RAY8(C,S)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    type(probe_8),INTENT(INOUT) ::S
    integer i,j


    call TRACK_SPIN_BACK(C,S%S%X)


  end subroutine TRACK_SPIN_BACK_RAY8



  SUBROUTINE TRACK_SPIN_BACKP(C,S)
    implicit none
    logical(lp) :: doneitf=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE



    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,S,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILT_SPIN(C%DIR,C%MAG%P%TILTD,2,S)


    IF((PATCHG==2).or.(PATCHG==3)) THEN
       CALL PATCH_SPIN(C,S,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely


  END SUBROUTINE TRACK_SPIN_BACKP




  SUBROUTINE PATCH_SPINR(C,s,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: s(3)
    logical(lp),INTENT(IN):: ENTERING
    real(dp) da
    IF(ENTERING) THEN
       da=C%PATCH%A_ANG(1)+((C%PATCH%A_X1-1)/2)*pi
       call rot_spin_x(s,da)
       call rot_spin_y(s,C%PATCH%A_ANG(2))
       da=C%PATCH%A_ANG(3)+((C%PATCH%A_X2-1)/2)*pi
       call rot_spin_z(s,da)
    ELSE
       da=C%PATCH%B_ANG(1)+((C%PATCH%B_X1-1)/2)*pi
       call rot_spin_x(s,da)
       call rot_spin_y(s,C%PATCH%A_ANG(2))
       da=C%PATCH%b_ANG(3)+((C%PATCH%B_X2-1)/2)*pi
       call rot_spin_z(s,da)
    ENDIF

  END SUBROUTINE PATCH_SPINR

  SUBROUTINE PATCH_SPINp(C,s,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    type(real_8), INTENT(INOUT):: s(3)
    logical(lp),INTENT(IN):: ENTERING
    real(dp) da

    IF(ENTERING) THEN
       da=C%PATCH%A_ANG(1)+((C%PATCH%A_X1-1)/2)*pi
       call rot_spin_x(s,da)
       call rot_spin_y(s,C%PATCH%A_ANG(2))
       da=C%PATCH%A_ANG(3)+((C%PATCH%A_X2-1)/2)*pi
       call rot_spin_z(s,da)
    ELSE
       da=C%PATCH%B_ANG(1)+((C%PATCH%B_X1-1)/2)*pi
       call rot_spin_x(s,da)
       call rot_spin_y(s,C%PATCH%A_ANG(2))
       da=C%PATCH%b_ANG(3)+((C%PATCH%B_X2-1)/2)*pi
       call rot_spin_z(s,da)
    ENDIF

  END SUBROUTINE PATCH_SPINp

  !   Misalignment routines
  SUBROUTINE MIS_SPINR(C,S,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    real(dp), INTENT(INOUT):: S(3)
    logical(lp),INTENT(IN):: ENTERING

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             call rot_spin_X(s,C%CHART%ANG_IN(1))
             call rot_spin_Y(s,C%CHART%ANG_IN(2))
             call rot_spin_Z(s,C%CHART%ANG_IN(3))
          ELSE
             call rot_spin_X(s,C%CHART%ANG_OUT(1))
             call rot_spin_Y(s,C%CHART%ANG_OUT(2))
             call rot_spin_Z(s,C%CHART%ANG_OUT(3))
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
             call rot_spin_Z(s,C%CHART%ANG_OUT(3))
             call rot_spin_Y(s,C%CHART%ANG_OUT(2))
             call rot_spin_X(s,C%CHART%ANG_OUT(1))
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
             call rot_spin_Z(s,C%CHART%ANG_IN(3))
             call rot_spin_Y(s,C%CHART%ANG_IN(2))
             call rot_spin_X(s,C%CHART%ANG_IN(1))
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_SPINR

  SUBROUTINE MIS_SPINP(C,S,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(REAL_8), INTENT(INOUT):: S(3)
    logical(lp),INTENT(IN):: ENTERING

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             call rot_spin_X(s,C%CHART%ANG_IN(1))
             call rot_spin_Y(s,C%CHART%ANG_IN(2))
             call rot_spin_Z(s,C%CHART%ANG_IN(3))
          ELSE
             call rot_spin_X(s,C%CHART%ANG_OUT(1))
             call rot_spin_Y(s,C%CHART%ANG_OUT(2))
             call rot_spin_Z(s,C%CHART%ANG_OUT(3))
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
             call rot_spin_Z(s,C%CHART%ANG_OUT(3))
             call rot_spin_Y(s,C%CHART%ANG_OUT(2))
             call rot_spin_X(s,C%CHART%ANG_OUT(1))
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
             call rot_spin_Z(s,C%CHART%ANG_IN(3))
             call rot_spin_Y(s,C%CHART%ANG_IN(2))
             call rot_spin_X(s,C%CHART%ANG_IN(1))
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_SPINP


  SUBROUTINE DTILTR_SPIN(DIR,TILTD,I,S)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: S(3)
    INTEGER,INTENT(IN):: I,DIR
    REAL(DP),INTENT(IN) :: TILTD
    real(dp) YS

    IF(TILTD==zero) RETURN
    IF(I==1) THEN
       YS=DIR*TILTD
       call rot_spin_Z(s,YS)
    ELSE
       YS=-DIR*TILTD
       call rot_spin_Z(s,YS)
    ENDIF

  END SUBROUTINE DTILTR_SPIN

  SUBROUTINE DTILTP_SPIN(DIR,TILTD,I,S)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: S(3)
    INTEGER,INTENT(IN):: I,DIR
    REAL(DP),INTENT(IN) :: TILTD
    real(dp) YS

    IF(TILTD==zero) RETURN
    IF(I==1) THEN
       YS=DIR*TILTD
       call rot_spin_Z(s,YS)
    ELSE
       YS=-DIR*TILTD
       call rot_spin_Z(s,YS)
    ENDIF

  END SUBROUTINE DTILTP_SPIN




end module ptc_spin
