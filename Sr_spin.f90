!   B~= (H B_x, H B_y,   B_z)
!   A~= (  A_x,   A_y, H A_z)
!   H = 1 + h x / RHO_0
!
module ptc_spin
  !use orbit_ptc
  !use beam_beam_ptc
  use accel_ptc
  implicit none
  public
  PRIVATE get_fieldR !,get_field
  PRIVATE get_BfieldR,get_BfieldP,get_Bfield,get_fieldp
  private GETMULB_TEAPOTR,GETMULB_TEAPOTP,get_Bfieldt
  private B_PANCAkEr,B_PANCAkEp,B_PANCAkE
  PRIVATE DIRECTION_VR,DIRECTION_VP,DIRECTION_V
  PRIVATE  B_PARA_PERP_r,B_PARA_PERP_p,B_PARA_PERP
  PRIVATE get_omega_spinR,get_omega_spinP,get_omega_spin
  PRIVATE PUSH_SPINR,PUSH_SPINP,PUSH_SPIN
  PRIVATE TRACK_FRINGE_spin_R,TRACK_FRINGE_spin_P,TRACK_FRINGE_spin
  PRIVATE TRACK_NODE_LAYOUT_FLAG_pr_s12_R,TRACK_NODE_LAYOUT_FLAG_pr_s12_P
  PRIVATE GET_BE_CAVR,GET_BE_CAVP ,GET_BE_CAV
  private rot_spin_x,rot_spin_xr,rot_spin_xp,rot_spin_z,rot_spin_zr,rot_spin_zp
  private rot_spin_yr,rot_spin_yp,rot_spin_y
  private PATCH_SPINR,PATCH_SPINP,PATCH_SPIN
  private MIS_SPINR,MIS_SPINP,MIS_SPIN
  private DTILTR_SPIN,DTILTP_SPIN,DTILT_SPIN
  PRIVATE TRACK_SPIN_FRONTR,TRACK_SPIN_FRONTP,TRACK_SPIN_FRONT
  PRIVATE TRACK_SPIN_BACKR,TRACK_SPIN_BACKP,TRACK_SPIN_BACK
  private PUSH_SPIN_RAY8,TRACK_SPIN_BACK_RAY8,TRACK_SPIN_FRONT_ray8
  private radiate_2p,radiate_2r,radiate_2
  private TRACK_NODE_FLAG_probe_R,TRACK_NODE_FLAG_probe_p,TRACK_NODE_LAYOUT_FLAG_spinr_x
  private FIND_ORBIT_LAYOUT_noda,TRACK_NODE_LAYOUT_FLAG_spinp_x
  PRIVATE get_Bfield_fringe_R,get_Bfield_fringe_p,get_Bfield_fringe
  private TRACK_LAYOUT_FLAG_spin12r_x,TRACK_LAYOUT_FLAG_spin12p_x
  PRIVATE TRACK_LAYOUT_FLAG_probe_spin12R,TRACK_LAYOUT_FLAG_probe_spin12P,TRACK_YS
  private PUSH_SPIN_fake_fringer,PUSH_SPIN_fake_fringep,PUSH_SPIN_fake_fringe
  PRIVATE TRACK_NODE_LAYOUT_FLAG_pr_t12_R,TRACK_NODE_LAYOUT_FLAG_pr_t12_P
  private TRACK_LAYOUT_FLAG_spint12r_x,TRACK_LAYOUT_FLAG_spint12p_x,alloc_temporal_beam
  private alloc_temporal_probe,FIND_ORBIT_LAYOUT_noda_spin
  !REAL(DP) :: AG=A_ELECTRON
  REAL(DP) :: bran_init=pi

  INTERFACE alloc
     MODULE PROCEDURE alloc_temporal_probe
     MODULE PROCEDURE alloc_temporal_beam
  END INTERFACE

  INTERFACE TRACK_PROBE2
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_s12_R  !#2  ! probe from node i1 to i2
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_s12_P
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_R  !#6 USING NODE1 TO NODE2 AS OBJECT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_P
  END INTERFACE

  INTERFACE TRACK_PROBE
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_PROBE_spin12R  !#3  ! probe from FIBRE
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_PROBE_spin12P    !
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_R  !#6 USING NODE1 TO NODE2 AS OBJECT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_P
  END INTERFACE

  INTERFACE TRACK_NODE_PROBE                  ! track probe in a node t
     MODULE PROCEDURE TRACK_NODE_FLAG_PROBE_R   ! #1
     MODULE PROCEDURE TRACK_NODE_FLAG_PROBE_p   ! #1p
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_R  !#6 USING NODE1 TO NODE2 AS OBJECT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_P
  END INTERFACE

  INTERFACE TRACK_node_x                !
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spinr_x  !#4 ! TRACK X THROUGH INTEGRATION_NODE T
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spinp_x
  END INTERFACE

  INTERFACE TRACK_node_v
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spin_v   !#5
  END INTERFACE

  INTERFACE TRACK
     MODULE PROCEDURE TRACK_YS
  END INTERFACE

  !     call TRACK_NODE_SINGLE(intnode,R%R,my_estate,my_ering%CHARGE)


  INTERFACE TRACK_probe_x
     !     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spinr_x  ! fibre i1 to i2
     !     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spinp_x
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spin12r_x    !#7
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spin12p_x
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spint12r_x
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spint12p_x
  END INTERFACE

  INTERFACE FIND_ORBIT_x
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
  END INTERFACE

  INTERFACE FIND_ORBIT_probe_x
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
  END INTERFACE


  INTERFACE FIND_ORBIT_probe_spin
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda_spin
  END INTERFACE


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
  END INTERFACE

  INTERFACE PUSH_SPIN
     MODULE PROCEDURE PUSH_SPINR
     MODULE PROCEDURE PUSH_SPINP
     MODULE PROCEDURE PUSH_SPIN_RAY8
  END INTERFACE

  INTERFACE get_omega_spin
     MODULE PROCEDURE get_omega_spinR
     MODULE PROCEDURE get_omega_spinP
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

  INTERFACE get_Bfield_fringe
     MODULE PROCEDURE get_Bfield_fringe_R   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE get_Bfield_fringe_p   ! MID DEFINED AS 1/2 L
  END INTERFACE

  INTERFACE PUSH_SPIN_fake_fringe
     MODULE PROCEDURE PUSH_SPIN_fake_fringer   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE PUSH_SPIN_fake_fringep   ! MID DEFINED AS 1/2 L
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

  subroutine PUSH_SPIN_RAY8(c,DS,FAC,S,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTP), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8) OM(3),CO(3),SI(3)
    TYPE(REAL_8) ST
    type(probe_8),INTENT(INOUT) ::S
    integer i,j
    LOGICAL(LP),intent(in) :: before
    type(internal_state) k
    el=>c%parent_fibre%magp
    if(EL%kind<=kind1) return
    call PUSH_SPIN(c,DS,FAC,S%S%x,S%X,S%E_IJ,before,k,POS)

  end subroutine PUSH_SPIN_RAY8

  subroutine radiate_2r(c,DS,FAC,X,b2,dlds,XP,before,k,POS)
    use gauss_dis
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    real(dp),INTENT(INOUT) :: X(6),XP(2)
    real(dp), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    real(dp), intent(in):: B2,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    real(dp)  st,z,av(3),t
    type(internal_state) k

    el=>c%parent_fibre%mag

    if(k%TIME) then
       ST=root(one+two*X(5)/EL%P%beta0+x(5)**2)-one
    else
       ST=X(5)
    endif

    ! X(5)=X(5)+B2*FAC*DS
    !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ! X(5)=one/(one/(one+X(5))+CRADF(EL%P)*(one+X(5))*B2*DLDS*FAC*DS)-one
    !        X(5)=X(5)-CRADF(EL%P)*(one+X(5))**3*B2*FAC*DS/SQRT((one+X(5))**2-X(2)**2-X(4)**2)
    X(5)=X(5)-CRADF(EL%P)*(one+X(5))**3*B2*FAC*DS*DLDS
    if(k%stochastic) then
       !         t=sqrt(12.e0_dp)*(bran(bran_init)-half)
       t=RANF()
       !         t=sqrt(12.d0)*(RANF()-half)
       if(t>half) then
          t=one
       else
          t=-one
       endif
       if(before) then
          x(5)=x(5)+t*c%delta_rad_in
       else
          x(5)=x(5)+t*c%delta_rad_out
       endif
    endif

    if(el%kind/=kindpa) then
       IF(ASSOCIATED(EL%B_SOL)) THEN
          if(k%TIME) then
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)*root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/two
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)*root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/two
          else
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)*(one+X(5))/(one+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/two
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)*(one+X(5))/(one+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/two
          endif
       ELSEif(el%kind==kind22) then

          IF(EL%HE22%P%DIR==1) THEN
             Z= pos*el%l/el%p%nst
          ELSE
             Z=EL%L-pos*el%l/el%p%nst
          ENDIF
          CALL compute_f4(EL%he22,X,Z,A=AV)
          if(k%TIME) then
             X(2)=(X(2)+EL%P%CHARGE*AV(1)) *root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2)) *root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          else
             X(2)=(X(2)+EL%P%CHARGE*AV(1))*(one+X(5))/(one+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2))*(one+X(5))/(one+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          endif


       ELSE
          if(k%TIME) then
             X(2)=X(2)*root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(4)=X(4)*root(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
          else
             X(2)=X(2)*(one+X(5))/(one+ST)
             X(4)=X(4)*(one+X(5))/(one+ST)
          endif
       ENDIF
    endif

    !       X(2)=X_MEC(2)*(one+X(5))/(one+X5)-EL%B_SOL*EL%P%CHARGE*X(3)/two
    !       X(4)=X_MEC(4)*(one+X(5))/(one+X5)+EL%B_SOL*EL%P%CHARGE*X(1)/two


  end subroutine radiate_2r

  subroutine PUSH_SPINR(c,DS,FAC,S,X,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    REAL(DP),INTENT(INOUT) :: X(6),S(3)
    REAL(DP), INTENT(IN) :: DS,FAC
    REAL(DP) OM(3),CO(3),SI(3),B2,XP(2)
    REAL(DP) ST,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k

    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    if(.not.(k%radiation.or.k%SPIN)) return
    el=>c%parent_fibre%mag
    if(EL%kind<=kind1) return

    CALL get_omega_spin(c,OM,B2,dlds,XP,X,POS,k)

    if(k%radiation.AND.BEFORE) then
       !if(el%p%radiation.AND.BEFORE) then
       call radiate_2(c,DS,FAC,X,b2,dlds,XP,before,k,POS)
    endif


    if(k%SPIN) then
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

    !if(el%p%radiation.AND.(.NOT.BEFORE)) then
    if(k%radiation.AND.(.NOT.BEFORE)) then
       call radiate_2(c,DS,FAC,X,b2,dlds,XP,before,k,POS)
    endif

  END subroutine PUSH_SPINR

  subroutine radiate_2p(c,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: X(6),XP(2)
    real(dp),INTENT(INOUT) :: E_IJ(6,6)
    TYPE(REAL_8), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8), intent(in):: B2,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    TYPE(REAL_8) st,av(3),z
    logical(lp) done_stoch
    real(dp) b30,x1,x3,denf,b3
    type(damap) xpmap
    integer i,j
    type(internal_state) k
    el=>c%parent_fibre%magp
    if(.not.before.and.stoch_in_rec) then

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
       c%delta_rad_out=root(denf)
       call kill(xpmap)
    endif


    call alloc(st)

    if(k%TIME) then
       ST=SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)-one
    else
       ST=X(5)
    endif

    ! X(5)=X(5)+B2*FAC*DS
    !   X(5)=one/(one/(one+X(5))-B2*FAC*DS)-one
    !   X(5)=one/(one/(one+X(5))+CRADF(EL%P)*(one+X(5))*B2*DLDS*FAC*DS)-one
    !          X(5)=X(5)-CRADF(EL%P)*(one+X(5))**3*B2*FAC*DS/SQRT((one+X(5))**2-X(2)**2-X(4)**2)
    X(5)=X(5)-CRADF(EL%P)*(one+X(5))**3*B2*FAC*DS*DLDS


    if(el%kind/=kindpa) then
       IF(ASSOCIATED(EL%B_SOL)) THEN
          if(k%TIME) then
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)*SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/two
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)*SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/two
          else
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)*(one+X(5))/(one+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/two
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)*(one+X(5))/(one+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/two
          endif
       ELSEif(el%kind==kind22) then
          call alloc(av,3)
          call alloc(z)

          IF(EL%HE22%P%DIR==1) THEN
             Z= pos*el%l/el%p%nst
          ELSE
             Z=EL%L-pos*el%l/el%p%nst
          ENDIF
          CALL compute_f4(EL%he22,X,Z,A=AV)
          if(k%TIME) then
             X(2)=(X(2)+EL%P%CHARGE*AV(1)) *sqrt(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2)) *sqrt(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          else
             X(2)=(X(2)+EL%P%CHARGE*AV(1))*(one+X(5))/(one+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2))*(one+X(5))/(one+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          endif
          call kill(av,3)
          call kill(z)
       ELSE
          if(k%TIME) then
             X(2)=X(2)*SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
             X(4)=X(4)*SQRT(one+two*X(5)/EL%P%beta0+x(5)**2)/(one+ST)
          else
             X(2)=X(2)*(one+X(5))/(one+ST)
             X(4)=X(4)*(one+X(5))/(one+ST)
          endif
       ENDIF
    endif

    call kill(st)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(before.and.stoch_in_rec) then
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
       c%delta_rad_in=root(denf)
       call kill(xpmap)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine radiate_2p

  subroutine PUSH_SPIN_fake_fringer(c,S,X,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    real(dp),INTENT(INOUT) :: X(6),S(3)

    real(dp) OM(3),CO(3),SI(3),B2,XP(2)
    real(dp) ST,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k

    if(.not.(k%radiation.or.k%SPIN)) return
    el=>c%parent_fibre%mag
    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    if(EL%kind<=kind1) return



    CALL get_omega_spin(c,OM,B2,dlds,XP,X,POS,k)
    !if(k%radiation.AND.BEFORE) then
    !if(el%p%radiation.AND.BEFORE) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif

    if(k%SPIN) then
       CO(1)=COS(OM(1)/TWO)
       SI(1)=SIN(OM(1)/TWO)
       CO(2)=COS(OM(2)/TWO)
       SI(2)=SIN(OM(2)/TWO)
       CO(3)=COS(OM(3))
       SI(3)=SIN(OM(3))

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
    !if(k%radiation.AND.(.NOT.BEFORE)) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif



  END subroutine PUSH_SPIN_fake_fringer

  subroutine PUSH_SPIN_fake_fringep(c,S,X,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTp), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: X(6),S(3)

    TYPE(REAL_8) OM(3),CO(3),SI(3),B2,XP(2)
    TYPE(REAL_8) ST,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k

    if(.not.(k%radiation.or.k%SPIN)) return
    el=>c%parent_fibre%magp
    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    if(EL%kind<=kind1) return

    CALL ALLOC(OM,3)
    CALL ALLOC(CO,3)
    CALL ALLOC(SI,3)
    CALL ALLOC(XP,2)
    CALL ALLOC(ST,B2,dlds)

    CALL get_omega_spin(c,OM,B2,dlds,XP,X,POS,k)
    !if(k%radiation.AND.BEFORE) then
    !if(el%p%radiation.AND.BEFORE) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif

    if(k%SPIN) then
       CO(1)=COS(OM(1)/TWO)
       SI(1)=SIN(OM(1)/TWO)
       CO(2)=COS(OM(2)/TWO)
       SI(2)=SIN(OM(2)/TWO)
       CO(3)=COS(OM(3))
       SI(3)=SIN(OM(3))

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
    !if(k%radiation.AND.(.NOT.BEFORE)) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif

    CALL KILL(OM,3)
    CALL KILL(CO,3)
    CALL KILL(SI,3)
    CALL KILL(XP,2)
    CALL KILL(ST,B2,dlds)

  END subroutine PUSH_SPIN_fake_fringep


  subroutine PUSH_SPINP(c,DS,FAC,S,X,E_IJ,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: X(6),S(3)
    real(dp),INTENT(INOUT) :: E_IJ(6,6)
    TYPE(REAL_8), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8) OM(3),CO(3),SI(3),B2,XP(2)
    TYPE(REAL_8) ST,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k
    if(.not.(k%radiation.or.k%SPIN)) return
    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    el=>c%parent_fibre%magp
    if(EL%kind<=kind1) return

    CALL ALLOC(OM,3)
    CALL ALLOC(CO,3)
    CALL ALLOC(SI,3)
    CALL ALLOC(XP,2)
    CALL ALLOC(ST,B2,dlds)

    CALL get_omega_spin(c,OM,B2,dlds,XP,X,POS,k)
    if(k%radiation.AND.BEFORE) then
       !if(el%p%radiation.AND.BEFORE) then
       call radiate_2(c,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)

    endif

    if(k%SPIN) then
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
    if(k%radiation.AND.(.NOT.BEFORE)) then
       !if(el%p%radiation.AND.(.NOT.BEFORE)) then
       call radiate_2(c,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    endif

    CALL KILL(OM,3)
    CALL KILL(CO,3)
    CALL KILL(SI,3)
    CALL KILL(XP,2)
    CALL KILL(ST,B2,dlds)

  END subroutine PUSH_SPINP

  subroutine get_omega_spinr(c,OM,B2,dlds,XP,X,POS,k)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    REAL(DP),INTENT(INOUT) :: X(6),OM(3),B2,XP(2),DLDS
    REAL(DP)  B(3),E(3),BPA(3),BPE(3),D1,D2,GAMMA,EB(3)
    REAL(DP) BETA0,GAMMA0I,XPA(2)
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    el=>c%parent_fibre%mag
    P=>EL%P
    P%DIR    => C%PARENT_FIBRE%DIR
    P%beta0  => C%PARENT_FIBRE%beta0
    P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    P%GAMBET => C%PARENT_FIBRE%GAMBET
    P%MASS => C%PARENT_FIBRE%MASS
    P%ag => C%PARENT_FIBRE%ag
    P%CHARGE=>C%PARENT_FIBRE%CHARGE
    ! DLDS IS  REALLY D(CT)/DS * (1/(ONE/BETA0+X(5)))
    OM(2)=ZERO
    EB=ZERO
    BPA=ZERO
    BPE=ZERO
    B=ZERO
    E=ZERO

    CALL get_field(EL,B,E,X,k,POS)

    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7) ! Straight for all practical purposes
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ENDIF

       OM(2)=P%b0
    case(KIND4) ! CAVITY
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,E,EB,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ENDIF
    case(KIND16:kind17,KIND20)
       CALL B_PARA_PERP(k,EL,0,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case(kind10)     ! TEAPOT real curvilinear
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ENDIF
       OM(2)=P%b0
    case(KINDPA)     ! fitted field for real magnet
       STOP 123   !  PATCH PROBLEMS???? CONVERTING TO PX ????
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,pos=POS)
       if(k%time) then
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
    case(KIND22)
       CALL B_PARA_PERP(k,EL,0,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case default
       OM(1)=ZERO
       OM(2)=ZERO
       OM(3)=ZERO
    END SELECT

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( ONE/P%BETA0 + X(5) )

    OM(1)=-DLDS*( (ONE+p%AG*GAMMA)*BPE(1) + (ONE+p%AG)*BPA(1) )
    OM(2)=-DLDS*( (ONE+p%AG*GAMMA)*BPE(2) + (ONE+p%AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*( (ONE+p%AG*GAMMA)*BPE(3) + (ONE+p%AG)*BPA(3) )
    DO I=1,3
       OM(I)=OM(I)-DLDS*(p%AG*GAMMA+GAMMA/(ONE+GAMMA))*EB(I)
    ENDDO
    if(k%RADIATION) then
       !      if(P%RADIATION) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF

  end subroutine get_omega_spinr

  subroutine get_omega_spinp(c,OM,B2,dlds,XP,X,POS,k)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTp), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8), INTENT(INOUT) :: X(6),OM(3),B2,XP(2)
    TYPE(REAL_8)  B(3),E(3),BPA(3),BPE(3),DLDS,D1,D2,GAMMA,EB(3),XPA(2)
    REAL(DP) BETA0,GAMMA0I
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    CALL ALLOC(B,3)
    CALL ALLOC(E,3)
    CALL ALLOC(EB,3)
    CALL ALLOC(BPA,3)
    CALL ALLOC(BPE,3)
    CALL ALLOC(XPA,2)
    CALL ALLOC(D1,D2,GAMMA)

    el=>c%parent_fibre%magp
    P=>EL%P
    P%DIR    => C%PARENT_FIBRE%DIR
    P%beta0  => C%PARENT_FIBRE%beta0
    P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    P%GAMBET => C%PARENT_FIBRE%GAMBET
    P%MASS => C%PARENT_FIBRE%MASS
    P%ag => C%PARENT_FIBRE%ag
    P%CHARGE=>C%PARENT_FIBRE%CHARGE
    ! DLDS IS  REALLY D(CT)/DS * (1/(ONE/BETA0+X(5)))
    OM(2)=ZERO

    CALL get_field(EL,B,E,X,k,POS)
    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7) ! Straight for all practical purposes
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ENDIF
       OM(2)=P%b0
    case(KIND4) ! CAVITY
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,E,EB,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-X(2)**2-X(4)**2)*(one+P%b0*X(1))
       ENDIF
    case(KIND16:kind17,KIND20)
       CALL B_PARA_PERP(k,el,0,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case(kind10)     ! TEAPOT real curvilinear
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)*(one+P%b0*X(1))
       ENDIF
       OM(2)=P%b0
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,pos=POS)
       if(k%time) then
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
    case(KIND22)
       CALL B_PARA_PERP(k,EL,0,X,B,BPA,BPE,XP,XPA,pos=POS)
       IF(k%TIME) THEN
          DLDS=ONE/SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=ONE/SQRT((ONE+X(5))**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case default
       OM(1)=ZERO
       OM(2)=ZERO
       OM(3)=ZERO
    END SELECT

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( ONE/P%BETA0 + X(5) )

    OM(1)=-DLDS*( (ONE+p%AG*GAMMA)*BPE(1) + (ONE+p%AG)*BPA(1) )
    OM(2)=-DLDS*( (ONE+p%AG*GAMMA)*BPE(2) + (ONE+p%AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*( (ONE+p%AG*GAMMA)*BPE(3) + (ONE+p%AG)*BPA(3) )
    !     DO I=1,3
    !       write(30,*) 'b, bpe and bpa',i,el%name,el%kind
    !       call print(b(i),30)
    !       call print(bpe(i),30)
    !       call print(bpa(i),30)
    !     enddo
    DO I=1,3
       OM(I)=OM(I)-DLDS*(p%AG*GAMMA+GAMMA/(ONE+GAMMA))*EB(I)
    ENDDO
    if(k%RADIATION) then
       !      if(P%RADIATION) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF
    CALL KILL(B,3)
    CALL KILL(E,3)
    CALL KILL(EB,3)
    CALL KILL(BPA,3)
    CALL KILL(BPE,3)
    CALL KILL(D1,D2,GAMMA)
    CALL KILL(XPA,2)

  end subroutine get_omega_spinp

  subroutine get_fieldr(EL,B,E,X,k,POS)
    implicit none
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    REAL(DP),INTENT(INOUT) :: B(3),E(3)
    REAL(DP),INTENT(INOUT) :: X(6)
    REAL(DP) Z
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    B=ZERO
    E=ZERO

    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7,KIND16:kind17,KIND20) ! Straight for all practical purposes
       if(present(pos)) then
          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,X,pos)   ! fringe effect
          ELSE
             call get_Bfield(EL,B,X)   ! fringe effect
          ENDIF
       else
          CALL get_Bfield(EL,B,X)
       endif
    case(kind10)     ! TEAPOT real curvilinear
       CALL get_Bfieldt(EL%TP10,B,X)
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PANCAkE(EL%PA,B,X,POS)
    case(KINDWIGGLER)
       CALL get_z_wi(EL%wi,POS,z)
       CALL B_FIELD(EL%wi,Z,X,B)

    CASE(KIND4)      ! Pill box cavity
       CALL GET_BE_CAV(EL%C4,B,E,X,k)
    CASE(KIND21)     ! travelling wave cavity
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    CASE(KIND22)     ! helical dipole
       IF(EL%HE22%P%DIR==1) THEN
          Z= pos*el%l/el%p%nst
       ELSE
          Z=EL%L-pos*el%l/el%p%nst
       ENDIF
       CALL compute_f4(EL%HE22,X,Z,B=B)    !    IF(EL%P%DIR==1) THEN
       !      write(6,*) z,z* EL%HE22%freq/twopi
       !      write(6,*) b
       !      pause  123
    case default


    END SELECT

    DO I=1,3
       B(I)=B(I)*EL%P%CHARGE
       E(I)=E(I)*EL%P%CHARGE
    ENDDO


  end subroutine get_fieldr

  subroutine get_fieldp(EL,B,E,X,k,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: B(3),E(3)
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    INTEGER I
    TYPE(REAL_8) z
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    DO I=1,3
       B(I)=ZERO
       E(I)=ZERO
    ENDDO

    SELECT CASE(EL%KIND)
    case(KIND2,kind5:kind7,KIND16:kind17,KIND20) ! Straight for all practical purposes
       if(present(pos)) then
          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,X,pos)   ! fringe effect
          ELSE
             call get_Bfield(EL,B,X)   ! fringe effect
          ENDIF
       else
          CALL get_Bfield(EL,B,X)
       endif
    case(kind10)     ! TEAPOT real curvilinear
       CALL get_Bfieldt(EL%TP10,B,X)
       !     WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PANCAkE(EL%PA,B,X,POS)
    case(KINDWIGGLER)
       call alloc(z)
       CALL get_z_wi(EL%wi,POS,z)
       CALL B_FIELD(EL%wi,Z,X,B)
       call kill(z)
    CASE(KIND4)      ! Pill box cavity
       CALL GET_BE_CAV(EL%C4,B,E,X,k)
    CASE(KIND21)     ! travelling wave cavity
       WRITE(6,*) EL%KIND,EL%NAME," NOT DONE "
    CASE(KIND22)     ! helical dipole
       IF(EL%HE22%P%DIR==1) THEN
          Z= pos*el%l/el%p%nst
       ELSE
          Z=EL%L-pos*el%l/el%p%nst
       ENDIF
       CALL compute_f4(EL%HE22,X,Z,B=B)    !    IF(EL%P%DIR==1) THEN
    case default
    END SELECT

    DO I=1,3
       B(I)=B(I)*EL%P%CHARGE
       E(I)=E(I)*EL%P%CHARGE
    ENDDO

  end subroutine get_fieldp

  SUBROUTINE get_Bfield_fringe_R(EL,B,X,pos)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENT),INTENT(IN):: EL
    INTEGER pos


    IF(ASSOCIATED(EL%B_SOL)) THEN
       B(1)=-(2*Pos+3)*EL%B_SOL*half*x(1);    ! POS =-2,-1  (ENT, EXIT)
       B(2)=-(2*Pos+3)*EL%B_SOL*half*x(3);
       B(3)=zero;
    else
       b(1)=zero
       b(2)=zero
       b(3)=zero
    ENDIF

  END SUBROUTINE get_Bfield_fringe_R


  SUBROUTINE get_Bfield_fringe_p(EL,B,X,pos)
    IMPLICIT NONE
    type(REAL_8),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENTP),INTENT(IN):: EL
    INTEGER pos


    IF(ASSOCIATED(EL%B_SOL)) THEN
       B(1)=-(2*Pos+3)*EL%B_SOL*half*x(1);    ! POS =-2,-1  (ENT, EXIT)
       B(2)=-(2*Pos+3)*EL%B_SOL*half*x(3);
       B(3)=zero;
    else
       b(1)=zero
       b(2)=zero
       b(3)=zero
    ENDIF

  END SUBROUTINE get_Bfield_fringe_p

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

  SUBROUTINE GET_BE_CAVR(EL,B,E,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3),E(3)
    TYPE(CAV4),INTENT(INOUT):: EL
    real(dp) DF,R2,F,DR2,O,VL
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    real(dp) BBYTWT,BBXTW,BBYTW,x1,x3
    integer J,dir,ko

    E=ZERO
    B=ZERO
    IF(k%NOCAVITY) RETURN

    O=twopi*EL%freq/CLIGHT
    VL=EL%volt*c_1d_3/EL%P%P0C
    do ko=1,el%nf

       DF=ZERO
       F=ONE
       R2=ONE

       DO I=1,EL%N_BESSEL
          R2=-R2*(ko*O)**2/FOUR/(I+1)**2
          DR2=R2*I
          DF=DF+DR2*2
          R2=R2*(X(1)**2+X(3)**2)
          F=F+R2
       ENDDO

       !    EL%DELTA_E=x(5)

       IF(EL%N_BESSEL>0) THEN
          B(2)=B(2)-EL%F(KO)*X(1)*DF*VL*COS(ko*O*X(6)+EL%PHAS+phase0)/(ko*O)
          B(1)=B(1)+EL%F(KO)*X(3)*DF*VL*COS(ko*O*X(6)+EL%PHAS+phase0)/(ko*O)
       ENDIF

       E(3)=E(3)-EL%F(KO)*F*VL*SIN(ko*O*x(6)+EL%PHAS+phase0)

       ! doing crabola

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

       ! multipole * cos(omega t+ phi)/p0c

       B(2)=B(2)+EL%F(KO)*BBYTW/EL%P%P0C*cos(ko*O*x(6)+EL%PHAS+EL%phase0)
       B(1)=B(1)+EL%F(KO)*BBXTW/EL%P%P0C*cos(ko*O*x(6)+EL%PHAS+EL%phase0)

       IF(EL%P%NMUL>=1) THEN
          BBYTW=-EL%BN(EL%P%NMUL)/EL%P%NMUL
          BBXTW=-EL%AN(EL%P%NMUL)/EL%P%NMUL


          DO  J=EL%P%NMUL,2,-1
             BBYTWT=X1*BBYTW-X3*BBXTW-EL%BN(J-1)/(J-1)
             BBXTW=X3*BBYTW+X1*BBXTW-EL%AN(J-1)/(J-1)
             BBYTW=BBYTWT
          ENDDO
          BBYTWT=X1*BBYTW-X3*BBXTW
          BBXTW=X3*BBYTW+X1*BBXTW
          BBYTW=BBYTWT
       ELSE
          BBYTW=zero
          BBXTW=zero
       ENDIF

       E(3)=E(3)+EL%F(KO)*ko*O*BBYTW/EL%P%P0C*sin(ko*O*x(6)+EL%PHAS+EL%phase0)
    enddo

  END SUBROUTINE GET_BE_CAVR

  SUBROUTINE GET_BE_CAVP(EL,B,E,X,k)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6),B(3),E(3)
    TYPE(CAV4p),INTENT(INOUT):: EL
    TYPE(REAL_8) DF,R2,F,DR2,O,VL
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    TYPE(REAL_8) BBYTWT,BBXTW,BBYTW,x1,x3
    integer J,dir,KO

    DO I=1,3
       E(I)=ZERO
       B(I)=ZERO
    ENDDO

    IF(k%NOCAVITY) RETURN
    CALL ALLOC(DF,R2,F,DR2,O,VL)
    CALL ALLOC(BBYTWT,BBXTW,BBYTW,x1,x3)

    O=twopi*EL%freq/CLIGHT
    VL=EL%volt*c_1d_3/EL%P%P0C
    do ko=1,el%nf

       DF=ZERO
       F=ONE
       R2=ONE

       DO I=1,EL%N_BESSEL
          R2=-R2*(ko*O)**2/FOUR/(I+1)**2
          DR2=R2*I
          DF=DF+DR2*2
          R2=R2*(X(1)**2+X(3)**2)
          F=F+R2
       ENDDO

       !    EL%DELTA_E=x(5)

       IF(EL%N_BESSEL>0) THEN
          B(2)=B(2)-EL%F(KO)*X(1)*DF*VL*COS(ko*O*X(6)+EL%PHAS+phase0)/(ko*O)
          B(1)=B(1)+EL%F(KO)*X(3)*DF*VL*COS(ko*O*X(6)+EL%PHAS+phase0)/(ko*O)
       ENDIF

       E(3)=E(3)-EL%F(KO)*F*VL*SIN(ko*O*x(6)+EL%PHAS+phase0)

       ! doing crabola

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

       ! multipole * cos(omega t+ phi)/p0c

       B(2)=B(2)+EL%F(KO)*BBYTW/EL%P%P0C*cos(ko*O*x(6)+EL%PHAS+EL%phase0)
       B(1)=B(1)+EL%F(KO)*BBXTW/EL%P%P0C*cos(ko*O*x(6)+EL%PHAS+EL%phase0)

       IF(EL%P%NMUL>=1) THEN
          BBYTW=-EL%BN(EL%P%NMUL)/EL%P%NMUL
          BBXTW=-EL%AN(EL%P%NMUL)/EL%P%NMUL


          DO  J=EL%P%NMUL,2,-1
             BBYTWT=X1*BBYTW-X3*BBXTW-EL%BN(J-1)/(J-1)
             BBXTW=X3*BBYTW+X1*BBXTW-EL%AN(J-1)/(J-1)
             BBYTW=BBYTWT
          ENDDO
          BBYTWT=X1*BBYTW-X3*BBXTW
          BBXTW=X3*BBYTW+X1*BBXTW
          BBYTW=BBYTWT
       ELSE
          BBYTW=zero
          BBXTW=zero
       ENDIF

       E(3)=E(3)+EL%F(KO)*ko*O*BBYTW/EL%P%P0C*sin(ko*O*x(6)+EL%PHAS+EL%phase0)
    enddo




    CALL KILL(BBYTWT,BBXTW,BBYTW,x1,x3)
    CALL KILL(DF,R2,F,DR2,O,VL)

  END SUBROUTINE GET_BE_CAVP

  subroutine B_PANCAkEr(EL,B,X,POS)
    IMPLICIT NONE
    real(dp), INTENT(INout) :: X(6)
    INTEGER, INTENT(IN) :: POS
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
    INTEGER, INTENT(IN) :: POS
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

  subroutine B_PARA_PERP_R(k,EL,TEAPOT_LIKE,X,B,BPA,BPE,XP,XPA,EF,EFB,POS)
    IMPLICIT NONE
    REAL(DP),  INTENT(INout) :: X(6)
    TYPE(ELEMENT),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    REAL(DP),  INTENT(INout) :: B(3),BPA(3),BPE(3),XP(2),XPA(2)
    REAL(DP),  OPTIONAL ::EF(3),EFB(3)
    integer, optional,intent(in) :: pos
    INTEGER TEAPOT_LIKE,i
    REAL(DP) e(3),be
    type(internal_state) k
    P=>EL%P

    !  this routines gives us  B parallel and B perpendicular
    ! Also if EF is present, E perpendicular times beta is return

    call DIRECTION_V(k,EL,TEAPOT_LIKE,X,E,XP,XPA,POS)

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

  subroutine B_PARA_PERP_p(k,EL,TEAPOT_LIKE,X,B,BPA,BPE,XP,XPA,EF,EFB,pos)
    IMPLICIT NONE
    type(real_8),  INTENT(INout) :: X(6)
    TYPE(ELEMENTP),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    type(real_8),  INTENT(INout) :: B(3),BPA(3),BPE(3),XP(2),XPA(2)
    type(real_8),  OPTIONAL ::EF(3),EFB(3)
    INTEGER TEAPOT_LIKE,i
    type(real_8) e(3),be
    type(internal_state) k
    integer, optional,intent(in) :: pos

    P=>EL%P

    call alloc(e,3); call alloc(be);

    call DIRECTION_V(k,EL,TEAPOT_LIKE,X,E,XP,XPA,POS)

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


  subroutine DIRECTION_VR(k,EL,TEAPOT_LIKE,X,E,XP,XPA,POS)
    IMPLICIT NONE
    REAL(DP),  INTENT(INout) :: X(6),XP(2),XPA(2)
    TYPE(ELEMENT),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    REAL(DP),  INTENT(INOUT) ::E(3)
    REAL(DP) N,H,DP1,A,AP,B,BP,z,AV(3)
    INTEGER TEAPOT_LIKE
    integer, optional,intent(in) :: pos
    type(internal_state) k

    P=>EL%P

    !    CALL COMPX(EL,Z,X,A,AP)
    !    X_MEC=zero
    !    X_MEC(2)=X(2)-A
    !    CALL COMPY(EL,Z,X,B,BP)
    !    X_MEC(4)=X(4)-B
    !    CALL B2PERP(EL%P,B_F,X_MEC,X5,B2)

    IF(k%TIME) THEN
       DP1=SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2)
    ELSE
       DP1=ONE+X(5)
    ENDIF

    IF(EL%KIND/=KINDPA) THEN

       IF(ASSOCIATED(EL%B_SOL)) THEN  !SOLENOID

          XPA(1)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)
          XPA(2)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N
       ELSEif(el%kind==kindwiggler) then
          CALL get_z_wi(EL%wi,POS,z)
          CALL COMPX(EL%wi,Z,X,A,AP)
          Xpa(1)=X(2)-A
          CALL COMPY(EL%wi,Z,X,B,BP)
          Xpa(2)=X(4)-B
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N

       ELSEif(el%kind==kind22) then

          IF(EL%HE22%P%DIR==1) THEN
             Z= pos*el%l/el%p%nst
          ELSE
             Z=EL%L-pos*el%l/el%p%nst
          ENDIF

          CALL compute_f4(EL%he22,X,Z,A=AV)
          Xpa(1)=X(2)-EL%P%CHARGE*AV(1)
          Xpa(2)=X(4)-EL%P%CHARGE*AV(2)
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N

       else


          N=SQRT(DP1**2-X(2)**2-X(4)**2)

          E(1)=X(2)/DP1
          E(2)=X(4)/DP1
          E(3)=N/DP1
          XPA(1)=X(2)
          XPA(2)=X(4)
          XP(1)=X(2)/N
          XP(2)=X(4)/N
       ENDIF

    ELSE    ! NON CANONICAL VARIABLES
       H=ONE+P%B0*TEAPOT_LIKE*X(1)
       N=SQRT(H**2+X(2)**2+X(4)**2)
       E(1)=X(2)/N
       E(2)=X(4)/N
       E(3)=H/N
       XPA(1)=X(2)
       XPA(2)=X(4)
       XP(1)=X(2)
       XP(2)=X(4)

    ENDIF

    E(1)=EL%P%dir*E(1)
    E(2)=EL%P%dir*E(2)
    E(3)=EL%P%dir*E(3)



  END subroutine DIRECTION_VR

  subroutine DIRECTION_VP(k,EL,TEAPOT_LIKE,X,E,XP,XPA,POS)
    IMPLICIT NONE
    type(real_8), INTENT(INout) :: X(6),XP(2),XPA(2)
    TYPE(ELEMENTP),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    type(real_8), INTENT(INOUT) ::E(3)
    type(real_8) N,H,DP1,A,AP,B,BP,z,AV(3)
    INTEGER TEAPOT_LIKE
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer, optional,intent(in) :: pos

    P=>EL%P

    CALL ALLOC(N,H,DP1 )

    IF(k%TIME) THEN
       DP1=SQRT(ONE+TWO*X(5)/P%BETA0+X(5)**2)
    ELSE
       DP1=ONE+X(5)
    ENDIF

    IF(EL%KIND/=KINDPA) THEN

       IF(ASSOCIATED(EL%B_SOL)) THEN  !SOLENOID

          XPA(1)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/two)
          XPA(2)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/two)
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N
       ELSEif(el%kind==kindwiggler) then
          call alloc(A,AP,B,BP,z)
          CALL get_z_wi(EL%wi,POS,z)
          CALL COMPX(EL%wi,Z,X,A,AP)
          Xpa(1)=X(2)-A
          CALL COMPY(EL%wi,Z,X,B,BP)
          Xpa(2)=X(4)-B
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N

          call kill(A,AP,B,BP,z)

       ELSEif(el%kind==kind22) then
          CALL ALLOC(AV,3)

          IF(EL%HE22%P%DIR==1) THEN
             Z= pos*el%l/el%p%nst
          ELSE
             Z=EL%L-pos*el%l/el%p%nst
          ENDIF
          CALL compute_f4(EL%he22,X,Z,A=AV)
          Xpa(1)=X(2)-EL%P%CHARGE*AV(1)
          Xpa(2)=X(4)-EL%P%CHARGE*AV(2)
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N
          CALL KILL(AV,3)

       else


          N=SQRT(DP1**2-X(2)**2-X(4)**2)

          E(1)=X(2)/DP1
          E(2)=X(4)/DP1
          E(3)=N/DP1
          XPA(1)=X(2)
          XPA(2)=X(4)
          XP(1)=X(2)/N
          XP(2)=X(4)/N
       ENDIF

    ELSE    ! NON CANONICAL VARIABLES
       H=ONE+P%B0*TEAPOT_LIKE*X(1)
       N=SQRT(H**2+X(2)**2+X(4)**2)
       E(1)=X(2)/N
       E(2)=X(4)/N
       E(3)=H/N
       XPA(1)=X(2)
       XPA(2)=X(4)
       XP(1)=X(2)
       XP(2)=X(4)

    ENDIF

    E(1)=EL%P%dir*E(1)
    E(2)=EL%P%dir*E(2)
    E(3)=EL%P%dir*E(3)


    CALL KILL(N,H,DP1 )

  END subroutine DIRECTION_VP

!!!!!!!!!!!!   GLOBAL TRACKING ROUTINES    !!!!!!!!!!!!

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_R(xs,k,fibre1,fibre2,node1,node2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    type(probe), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE),optional, POINTER :: node1,node2
    TYPE (fibre),optional, POINTER :: fibre1,fibre2
    TYPE (INTEGRATION_NODE), POINTER :: C,n1,n2
    !    INTEGER,TARGET :: CHARGE

    !    if(present(node1))CHARGE=NODE1%PARENT_FIBRE%CHARGE
    !    if(present(fibre1))CHARGE=fibre1%CHARGE

    CALL RESET_APERTURE_FLAG


    xs%u=my_false

    if(present(node1)) n1=>node1
    if(present(node2)) n2=>node2
    if(present(fibre1)) n1=>fibre1%t1
    if(present(fibre2)) n2=>fibre2%t1

    c=>n1



    DO  WHILE(.not.ASSOCIATED(C,n2))

       CALL TRACK_NODE_PROBE(C,XS,K)

       C=>C%NEXT
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_R

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_P(xs,k,fibre1,fibre2,node1,node2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    type(probe_8), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE),optional, POINTER :: node1,node2
    TYPE (fibre),optional, POINTER :: fibre1,fibre2
    TYPE (INTEGRATION_NODE), POINTER :: C,n1,n2
    !    INTEGER,TARGET :: CHARGE

    !    if(present(node1))CHARGE=NODE1%PARENT_FIBRE%CHARGE
    !    if(present(fibre1))CHARGE=fibre1%CHARGE

    CALL RESET_APERTURE_FLAG


    xs%u=my_false

    if(present(node1)) n1=>node1
    if(present(node2)) n2=>node2
    if(present(fibre1)) n1=>fibre1%t1
    if(present(fibre2)) n2=>fibre2%t1

    c=>n1



    DO  WHILE(.not.ASSOCIATED(C,n2))

       CALL TRACK_NODE_PROBE(C,XS,K)

       C=>C%NEXT
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_P



  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_R(R,xs,k,I1,I2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C

    CALL RESET_APERTURE_FLAG
    xs%u=my_false

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )


    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       CALL TRACK_NODE_PROBE(C,XS,K)

       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_R


  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_P(R,XS,k,I1,I2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(probe_8), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J   ,i22
    TYPE (INTEGRATION_NODE), POINTER :: C
    real(dp) fac


    CALL RESET_APERTURE_FLAG
    xs%u=my_false

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       CALL TRACK_NODE_PROBE(C,XS,K)  !,R%charge)

       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.


  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_P


  SUBROUTINE TRACK_LAYOUT_FLAG_probe_spin12r(r,xS,k,fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe),intent(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    integer i1,i2
    integer i11,i22
    type(fibre), pointer:: p

    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)
    i1=0
    i2=0
    i11=0
    i22=0
    if(present(node1)) i11=node1
    if(present(node2)) i22=node2
    if(present(fibre1)) then
       i1=fibre1
       CALL move_to( R,p,I1)
       i11=p%t1%pos
       if(fibre1>r%n) i11=i11+int(real(fibre1,kind=dp)/real(r%n,kind=dp))*r%t%n
    endif
    if(present(fibre2)) then
       i2=fibre2
       CALL move_to( R,p,I2)
       i22=p%t1%pos
       if(fibre2>r%n) i22=i22+int(real(fibre2,kind=dp)/real(r%n,kind=dp))*r%t%n
    endif

    IF(I22==0) then
       IF(R%CLOSED) THEN
          I22=I11+R%T%N
       ELSE
          I22=1+R%T%N
       ENDIF
    endif

    !     write(6,*) 'probe ',i11,i22
    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N
    !     write(6,*) 'probe ',i11,i22


    CALL TRACK_PROBE2(r,xs,K,i11,i22)


  END SUBROUTINE TRACK_LAYOUT_FLAG_probe_spin12r

  SUBROUTINE TRACK_LAYOUT_FLAG_probe_spin12P(r,xS,k,fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe_8),intent(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    integer i1,i2
    integer i11,i22
    type(fibre), pointer:: p

    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)
    i1=0
    i2=0
    i11=0
    i22=0
    if(present(node1)) i11=node1
    if(present(node2)) i22=node2
    if(present(fibre1)) then
       i1=fibre1
       CALL move_to( R,p,I1)
       i11=p%t1%pos
       if(fibre1>r%n) i11=i11+int(real(fibre1,kind=dp)/real(r%n,kind=dp))*r%t%n
    endif
    if(present(fibre2)) then
       i2=fibre2
       CALL move_to( R,p,I2)
       i22=p%t1%pos
       if(fibre2>r%n) i22=i22+int(real(fibre2,kind=dp)/real(r%n,kind=dp))*r%t%n
    endif

    IF(I22==0) then
       IF(R%CLOSED) THEN
          I22=I11+R%T%N
       ELSE
          I22=1+R%T%N
       ENDIF
    endif

    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE2(r,xs,K,i11,i22)


  END SUBROUTINE TRACK_LAYOUT_FLAG_probe_spin12P

  SUBROUTINE TRACK_LAYOUT_FLAG_spin12r_x(r,x,k,u,t, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe) xs
    real(dp),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    logical(lp), optional ::u
    type(integration_node),optional, pointer :: t

    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)

    xs%u=my_false
    XS=X
    if(present(t)) THEN
       !    ALLOCATE(xs%lost_node)
       !    t=>xs%lost_node
       nullify(t)
    ENDIF
    !    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE(r,xs,K, fibre1,fibre2,node1,node2)
    if(present(u)) u=xs%u
    X=XS%X
    if(present(t)) THEN
       t=>xs%lost_node
       ! deallocate(xs%lost_node)
       NULLIFY(xs%lost_node)
    ENDIF

  END SUBROUTINE TRACK_LAYOUT_FLAG_spin12r_x

  SUBROUTINE TRACK_LAYOUT_FLAG_spin12p_x(r,x,k,u,t, fibre1,fibre2,node1,node2)  ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe_8) xs
    type(real_8),target,intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    logical(lp), optional ::u
    type(integration_node),optional, pointer :: t

    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)
    call alloc(xs)
    xs%u=my_false
    XS%x=X
    if(present(t)) THEN
       !    ALLOCATE(xs%lost_node)
       !    t=>xs%lost_node
       nullify(t)
    ENDIF

    !    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE(r,xs,K, fibre1,fibre2,node1,node2)
    if(present(u)) u=xs%u
    if(present(t)) THEN
       t=>xs%lost_node
       !       deallocate(xs%lost_node)
       NULLIFY(xs%lost_node)
    ENDIF

    X=XS%X
    call kill(xs)
  END SUBROUTINE TRACK_LAYOUT_FLAG_spin12p_x

  SUBROUTINE TRACK_LAYOUT_FLAG_spint12r_x(x,k,u,t, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    type(probe) xs
    real(dp),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    type(fibre),optional,pointer:: fibre1,fibre2
    type(integration_node),optional,pointer:: node1,node2
    logical(lp), optional ::u
    type(integration_node),optional, pointer :: t


    xs%u=my_false
    XS=X
    if(present(t)) THEN
       !       ALLOCATE(xs%lost_node)
       !       t=>xs%lost_node
       nullify(t)
    ENDIF
    !    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE(xs,K, fibre1,fibre2,node1,node2)
    if(present(u)) u=xs%u
    X=XS%X
    if(present(t)) THEN
       t=>xs%lost_node
       !       deallocate(xs%lost_node)
       NULLIFY(xs%lost_node)
    ENDIF

  END SUBROUTINE TRACK_LAYOUT_FLAG_spint12r_x

  SUBROUTINE TRACK_LAYOUT_FLAG_spint12p_x(x,k,u,t, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    type(probe_8) xs
    type(real_8),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    type(fibre),optional,pointer:: fibre1,fibre2
    type(integration_node),optional,pointer:: node1,node2
    logical(lp), optional ::u
    type(integration_node),optional, pointer :: t


    call alloc(xs)
    xs%u=my_false
    XS%X=X
    if(present(t)) THEN
       !       ALLOCATE(xs%lost_node)
       !       t=>xs%lost_node
       nullify(t)
    ENDIF
    !    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE(xs,K, fibre1,fibre2,node1,node2)
    if(present(u)) u=xs%u
    X=XS%X
    if(present(t)) THEN
       t=>xs%lost_node
       !      deallocate(xs%lost_node)
       NULLIFY(xs%lost_node)
    ENDIF
    call kill(xs)

  END SUBROUTINE TRACK_LAYOUT_FLAG_spint12p_x

  SUBROUTINE TRACK_ys(r,ys,i1,k)  ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    integer,INTENT(IN):: i1
    type(probe_8) xs
    type(env_8),target,intent(INOUT) ::  ys(6)
    TYPE(INTERNAL_STATE) K
    integer i,j
    logical(lp) stoch
    stoch=stoch_in_rec
    stoch_in_rec=.true.
    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)
    call alloc(xs)
    xs%u=my_false

    do i=1,6
       XS%x(i)=ys(i)%v
       do j=1,6
          xs%e_ij(i,j)=ys(i)%e(j)
       enddo
    enddo

    CALL TRACK_PROBE(r,xs,K, fibre1=i1)

    do i=1,6
       ys(i)%v=XS%x(i)
       do j=1,6
          ys(i)%e(j)=xs%e_ij(i,j)
       enddo
    enddo
    stoch_in_rec=stoch

    call kill(xs)
  END SUBROUTINE TRACK_ys

  SUBROUTINE TRACK_fill_ref(r,fix,i1,k)  ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    integer,INTENT(IN):: i1
    real(dp), intent(INOUT) ::  fix(6)
    real(dp)   x(6)
    TYPE(INTERNAL_STATE) K
    integer i,ino1
    type(fibre), pointer :: p
    type(integration_node), pointer :: t


    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)

    x=fix

    call move_to(r,p,i1)
    ino1=p%t1%pos

    write(6,*) " Fibre ",i1, p%mag%name
    write(6,*) " Node ",ino1

    t=>p%t1
    do i=ino1,ino1+r%t%n
       t%ref(1)=x(1)
       t%ref(2)=x(3)
       CALL TRACK_PROBE_x(r,x,k, node1=i,node2=i+1)
       t%ref(3)=x(1)
       t%ref(4)=x(3)
       t=>t%next
    enddo

    write(6,*) " done "

  END SUBROUTINE TRACK_fill_ref




  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spin_v(T,v,k,ref) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE),pointer:: T
    type(probe) xs,XS_REF
    type(three_d_info),intent(INOUT) ::  v
    TYPE(INTERNAL_STATE) K
    REAL(DP) SC,x(6),reference_ray(6)
    TYPE(INTEGRATION_NODE),POINTER:: mag_in,mag_out
    logical(lp), optional :: ref
    logical(lp) ref0

    ref0=my_false
    if(present(ref)) ref0=ref

    IF(.NOT.ASSOCIATED(T%B)) THEN
       call FILL_SURVEY_DATA_IN_NODE_LAYOUT(t%parent_fibre%parent_LAYOUT)
       WRITE(6,*)  " SURVEY DONE FOR THIN LAYOUT IN TRACK_NODE_LAYOUT_FLAG_spin_v "
    ENDIF

    if(.not.check_stable) return
    xs%u=my_false
    xs_ref%u=my_false

    XS=V%X
    if(.not.ref0) then
       XS_REF=V%reference_ray
    endif

    X=V%X
    if(.not.ref0) then
       reference_ray=V%reference_ray
    else
       reference_ray=zero
       reference_ray(1)=t%ref(1)
       reference_ray(3)=t%ref(2)
    endif

    CALL TRACK_NODE_PROBE(T,xs,K)  !,t%parent_fibre%CHARGE)

    if(.not.ref0) CALL TRACK_NODE_PROBE(T,XS_REF,K)  !,t%parent_fibre%CHARGE)

    v%u(1)=XS%u

    if(.not.ref0) then
       v%u(2)=XS_REF%u
       v%reference_ray=XS_REF%x
    else
       v%u(2)=my_false
       v%reference_ray=zero
       v%reference_ray(1)=t%ref(3)
       v%reference_ray(3)=t%ref(4)
    endif

    v%x=XS%x

    IF(V%U(1).OR.V%U(2)) RETURN


    SC=ONE
    IF(v%SCALE/=zero) SC=v%SCALE
    !      t=>B%POS(1)%NODE%previous

    V%r0=t%A+(reference_ray(1)-SC*reference_ray(1))*t%ENT(1,1:3)+ SC*X(1)*t%ENT(1,1:3)
    V%r0=v%r0+(reference_ray(3)-SC*reference_ray(3))*t%ENT(2,1:3)+ SC*X(3)*t%ENT(2,1:3)

    V%r=t%B+(V%reference_ray(1)-SC*V%reference_ray(1))*t%EXI(1,1:3)+ SC*V%X(1)*t%EXI(1,1:3)
    V%r=v%r+(V%reference_ray(3)-SC*V%reference_ray(3))*t%EXI(2,1:3)+ SC*V%X(3)*t%EXI(2,1:3)
    mag_in=>t%parent_fibre%t1%next%next
    mag_out=>t%parent_fibre%t2%previous%previous
    v%a=mag_in%a
    v%ent=mag_in%ent
    v%b=mag_in%b
    v%exi=mag_in%exi
    v%o=t%B
    v%mid=t%exi


    IF(MAG_IN%PREVIOUS%CAS/=CASE1) STOP 201
    IF(MAG_OUT%NEXT%CAS/=CASE2) STOP 202

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spin_v

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spinr_x(T,x,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE),pointer:: T
    type(probe) xs
    real(dp),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K

    if(.not.check_stable) return
    xs%u=my_false
    XS=X

    CALL TRACK_NODE_PROBE(T,xs,K)  !t%parent_fibre%CHARGE)

    X=XS%X



  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spinr_x

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spinp_x(T,x,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE),pointer:: T
    type(probe_8) xs
    type(real_8),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K

    if(.not.check_stable) return
    xs%u=my_false
    XS%x=X

    CALL TRACK_NODE_PROBE(T,xs,K) !,t%parent_fibre%CHARGE)

    X=XS%X



  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_spinp_x




  !  FUNDAMENTAL TRACKING ROUTINES




  SUBROUTINE TRACK_NODE_FLAG_probe_R(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    REAL(DP) FAC,DS


    if(xs%u) return
    C%PARENT_FIBRE%MAG%P%DIR    => C%PARENT_FIBRE%DIR
    C%PARENT_FIBRE%MAG%P%beta0  => C%PARENT_FIBRE%beta0
    C%PARENT_FIBRE%MAG%P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    C%PARENT_FIBRE%MAG%P%GAMBET => C%PARENT_FIBRE%GAMBET
    C%PARENT_FIBRE%MAG%P%MASS => C%PARENT_FIBRE%MASS
    C%PARENT_FIBRE%MAG%P%ag => C%PARENT_FIBRE%ag
    C%PARENT_FIBRE%MAG%P%CHARGE=>C%PARENT_FIBRE%CHARGE
    !      ag=xs%s%g
    !      if(associated(c%bb)) call BBKICK(c%BB,XS%X)

    if(c%cas==0) then
       ds=c%parent_fibre%MAG%L/c%parent_fibre%MAG%p%nst
       fac=half
       !call PUSH_SPIN(c%parent_fibre%mag,ds,FAC,XS%S%X,XS%X,my_true,k,C%POS_IN_FIBRE-3)
       call PUSH_SPIN(c,ds,FAC,XS%S%X,XS%X,my_true,k,C%POS_IN_FIBRE-3)
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
       !call PUSH_SPIN(c%parent_fibre%mag,ds,FAC,XS%S%X,XS%X,my_false,k,C%POS_IN_FIBRE-2)
       call PUSH_SPIN(c,ds,FAC,XS%S%X,XS%X,my_false,k,C%POS_IN_FIBRE-2)
    elseIF(c%cas==case1.or.c%cas==case2) then
       CALL TRACK_FRINGE_spin(C,XS%X,XS%S%X,K)
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    else
       IF(c%cas==caseP1) THEN
          CALL TRACK_SPIN_FRONT(C%PARENT_FIBRE,XS%S%X)
       ELSE
          CALL TRACK_SPIN_BACK(C%PARENT_FIBRE,XS%S%X)
       ENDIF
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE

    endif
    xs%u=.not.check_stable
    if(xs%u) xs%lost_node=>c

  END SUBROUTINE TRACK_NODE_FLAG_probe_R

  SUBROUTINE TRACK_NODE_FLAG_probe_P(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe_8), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    REAL(DP) FAC
    type(real_8) ds

    if(xs%u) return
    C%PARENT_FIBRE%MAGp%P%DIR    => C%PARENT_FIBRE%DIR
    C%PARENT_FIBRE%MAGp%P%beta0  => C%PARENT_FIBRE%beta0
    C%PARENT_FIBRE%MAGp%P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    C%PARENT_FIBRE%MAGp%P%GAMBET => C%PARENT_FIBRE%GAMBET
    C%PARENT_FIBRE%MAGP%P%MASS => C%PARENT_FIBRE%MASS
    C%PARENT_FIBRE%MAGP%P%ag => C%PARENT_FIBRE%ag
    C%PARENT_FIBRE%MAGp%P%CHARGE=>C%PARENT_FIBRE%CHARGE
    !      ag=xs%s%g

    CALL ALLOC(DS)

    !      if(associated(c%bb)) call BBKICK(c%BB,XS%X)


    if(c%cas==0) then
       ds=c%parent_fibre%MAGp%L/c%parent_fibre%MAG%p%nst
       fac=half

       !call PUSH_SPIN(c%parent_fibre%magP,ds,FAC,XS,my_true,k,C%POS_IN_FIBRE-3)
       call PUSH_SPIN(c,ds,FAC,XS,my_true,k,C%POS_IN_FIBRE-3)
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
       call PUSH_SPIN(c,ds,FAC,XS,my_false,k,C%POS_IN_FIBRE-2)
    elseIF(c%cas==case1.or.c%cas==case2) then
       CALL TRACK_FRINGE_spin(C,XS%X,XS%S%X,K)
       !        CALL TRACK_FRINGE_spin(C,XS,K)
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    else
       IF(c%cas==caseP1) THEN
          CALL TRACK_SPIN_FRONT(C%PARENT_FIBRE,XS)
       ELSE
          CALL TRACK_SPIN_BACK(C%PARENT_FIBRE,XS)
       ENDIF
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    endif


    call kill(ds)
    xs%u=.not.check_stable
    if(xs%u) xs%lost_node=>c

  END SUBROUTINE TRACK_NODE_FLAG_probe_P



  SUBROUTINE TRACK_FRINGE_spin_R(C,X,S,K)
    IMPLICIT NONE
    real(dp), INTENT(INOUT):: X(6),S(3)
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENT), POINTER :: EL
    integer pos
    el=>C%PARENT_FIBRE%MAG
    SELECT CASE(EL%KIND)
    CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
    case(KIND2)
    case(KIND4)
    case(KIND5,KIND17)
       IF(C%CAS==CASE1) THEN
          pos=-2
          call PUSH_SPIN_fake_fringe(c,S,X,my_true,k,pos)
       elseif(C%CAS==CASE2) then
          pos=-1
          call PUSH_SPIN_fake_fringe(c,S,X,my_false,k,pos)
       endif
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

    case(KIND21,kind22)
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
    integer pos

    el=>C%PARENT_FIBRE%MAGP
    SELECT CASE(EL%KIND)
    CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
    case(KIND2)
    case(KIND4)
    case(KIND5,KIND17)
       IF(C%CAS==CASE1) THEN
          pos=-2
          call PUSH_SPIN_fake_fringe(c,S,X,my_true,k,pos)
       elseif(C%CAS==CASE2) then
          pos=-1
          call PUSH_SPIN_fake_fringe(c,S,X,my_false,k,pos)
       endif
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

    case(KIND21,kind22)
    case(KINDWIGGLER)
    case(KINDPA)
    CASE DEFAULT
       WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
       stop 666
    END SELECT


  END SUBROUTINE TRACK_FRINGE_spin_P

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
       call rot_spin_z(s,C%PATCH%A_ANG(3))
       da=((C%PATCH%A_X2-1)/2)*pi
       call rot_spin_x(s,da)
    ELSE
       da=C%PATCH%B_ANG(1)+((C%PATCH%B_X1-1)/2)*pi
       call rot_spin_x(s,da)
       call rot_spin_y(s,C%PATCH%A_ANG(2))
       call rot_spin_z(s,C%PATCH%b_ANG(3))
       da=((C%PATCH%B_X2-1)/2)*pi
       call rot_spin_x(s,da)
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
       call rot_spin_z(s,C%PATCH%A_ANG(3))
       da=((C%PATCH%A_X2-1)/2)*pi
       call rot_spin_x(s,da)
    ELSE
       da=C%PATCH%B_ANG(1)+((C%PATCH%B_X1-1)/2)*pi
       call rot_spin_x(s,da)
       call rot_spin_y(s,C%PATCH%A_ANG(2))
       call rot_spin_z(s,C%PATCH%b_ANG(3))
       da=((C%PATCH%B_X2-1)/2)*pi
       call rot_spin_x(s,da)
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

  SUBROUTINE FIND_ORBIT_LAYOUT_noda_spin(RING,FIX,STATE,eps,TURNS,fibre1,node1,theta0) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    type(probe) , intent(inOUT) :: FIX
    type(probe)  FIXT
    INTEGER , optional,intent(in) :: TURNS,node1,fibre1
    real(dp)  eps,x0(6),theta0r,s(3,3)
    real(dp), optional,intent(inout) :: theta0
    TYPE(INTERNAL_STATE), intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat
    type(probe_8) xs
    type(damapspin) ds
    integer i,turns0,J
    turns0=1
    theta0r=0.d0
    if(present(turns)) turns0=turns
    x0=FIX%x
    stat=state
    call find_orbit_x(RING,x0,STATE,eps,TURNS,fibre1,node1)

    stat=state+spin0+spin_only0

    FIX=x0

    call init(stat,1,0)

    call alloc(xs)
    call alloc(ds)

    !        ds=1
    !        xs= FIX + ds

    !      do i=1,TURNS0
    !       call track_probe(RING,xs,stat,fibre1,node1)
    !        enddo

    !        ds=xs

    DO J=1,3
       FIXT%X=x0
       FIXT%S%X=ZERO
       FIXT%S%X(J)=one

       do i=1,TURNS0
          call track_probe(RING,fixt,stat,fibre1,node1)
       enddo

       s(1:3,j)=FIXT%S%X
    ENDDO


    call get_spin_nx(S,theta0r,FIX%s%x)
    if(present(theta0)) theta0=theta0r

    call kill(xs)
    call kill(ds)



  END SUBROUTINE FIND_ORBIT_LAYOUT_noda_spin


  SUBROUTINE FIND_ORBIT_LAYOUT_noda(RING,FIX,STATE,eps,TURNS,fibre1,node1) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    real(dp) , intent(inOUT) :: FIX(6)
    INTEGER , optional,intent(in) :: TURNS,node1,fibre1
    real(dp)  eps
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat

    real(dp)  DIX(6),xdix,xdix0,tiny,freq
    real(dp) X(6),Y(6),MX(6,6),sxi(6,6),SX(6,6)
    integer NO1,ND2,I,IU,ITE,ier,j,ITEM,loct,LOC
    TYPE (fibre), POINTER :: C
    TYPE (integration_node), POINTER :: t
    logical(lp) APERTURE
    INTEGER TURNS0,trackflag
    !!    type(probe) xs
    if(.not.associated(RING%t)) call MAKE_NODE_LAYOUT(ring)
    !!    xs%x=zero
    !!    xs%s%x=zero

    TURNS0=1
    trackflag=0
    IF(PRESENT(TURNS)) TURNS0=TURNS
    freq=zero
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    !!    call move_to(ring,c,loc)
    !!    loct=c%t1%pos


    Nullify(C);

    if(.not.ring%closed) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=" This line is not ring : FIND_ORBIT_LAYOUT_noda "
       call write_e(100)
    endif
    dix(:)=zero
    tiny=c_1d_40
    xdix0=c_1d4*DEPS_tracking
    NO1=1
    if(.not.present(STATE)) then
       IF(default%NOCAVITY) THEN
          !    ND1=2
          stat=default+    only_4d
       ELSE
          !   ND1=3
          STAT=default
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72))'
          w_p%c(1)=" No Cavity in the Line "
          w_p%c(2)=" FIND_ORBIT_LAYOUT will crash "
          call write_e(101)
       ENDIF
    else
       IF(STATE%NOCAVITY) THEN
          ND2=4
          STAT=STATE+only_4d0
       ELSE
          ND2=6
          STAT=STATE
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72))'
          w_p%c(1)=" No Cavity in the Line "
          w_p%c(2)=" FIND_ORBIT_LAYOUT will crash "
          call write_e(112)
       ENDIF
    endif
101 continue


    if((stat%totalpath==1).and.(.not.stat%nocavity)) then
       C=>RING%START
       freq=zero
       i=1
       xdix=zero
       do while(i<=RING%n)
          if(associated(c%magp%freq)) then
             IF(FREQ==ZERO) THEN
                freq=c%magp%freq
             ELSEIF(c%magp%freq/=ZERO.AND.c%magp%freq<FREQ) THEN
                freq=c%magp%freq
             ENDIF
          endif
          IF(stat%TIME) THEN
             XDIX=XDIX+c%mag%P%LD/c%mag%P%BETA0
          ELSE
             XDIX=XDIX+c%mag%P%LD
          ENDIF


          c=>c%next
          i=i+1
       enddo
       if(freq==zero) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(113)
       endif
       IF(RING%HARMONIC_NUMBER>0) THEN
          FREQ=RING%HARMONIC_NUMBER*CLIGHT/FREQ
          STOP 476
       ELSE
          !   fancy way
          stat=stat+nocavity0
          x=fix
          if(present(fibre1)) then
             call FIND_ORBIT(RING,x,fibre1,stat,1.e-6_dp)
             x(6)=0.d0
             call track(RING,x,fibre1,stat)
             xdix=x(6)
          else
             CALL move_to_INTEGRATION_NODE( Ring%T,t,node1 )
             call FIND_ORBIT(RING,x,t%parent_fibre%pos,stat,1.e-6_dp)
             x(6)=0.d0
             call track(RING,x,t%parent_fibre%pos,stat)
             xdix=x(6)
          endif
          stat=stat-nocavity0
          !   end of fancy way
          XDIX=XDIX*FREQ/CLIGHT
          write(6,*) NINT(xdix)
          FREQ=NINT(XDIX)*CLIGHT/FREQ
          WRITE(6,*) " FREQ ",FREQ
       ENDIF
    endif




    ITEM=0
3   continue
    ITEM=ITEM+1
    X=FIX

    DO I=1,TURNS0
       !       CALL TRACK(RING,X,LOC,STAT)
       !       trackflag=TRACK_flag(RING,X,LOC,STAT)
       !!       xs%x=x
       call TRACK_probe_X(Ring,x,stat,fibre1=fibre1,node1=node1)
       !!       call TRACK_PROBE(Ring,xs,loct,loct+ring%t%n,stat)
       !!       x=xs%x
       !       if(trackflag/=0) then
       !         ITEM=MAX_FIND_ITER+100
       !       endif

    ENDDO
    !    write(6,*) x
    x(6)=x(6)-freq*turns0

    mx=zero
    DO J=1,ND2

       Y=FIX
       Y(J)=FIX(J)+EPS
       DO I=1,TURNS0
          !          CALL TRACK(RING,Y,LOC,STAT)
          !!       xs%x=y
          call TRACK_probe_X(Ring,Y,stat,fibre1=fibre1,node1=node1)
          !!       call TRACK_PROBE(Ring,xs,loct,loct+ring%t%n,stat)
          !!       y=xs%x
          !          if(.not.check_stable) then
          !             w_p=0
          !             w_p%nc=1
          !             w_p%fc='((1X,a72))'
          !             write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",3
          !             call write_i

          !             return
          !          endif
       ENDDO
       y(6)=y(6)-freq*turns0

       do i=1,ND2
          MX(I,J)=(Y(i)-X(i))/eps
       enddo

    ENDDO

    SX=MX;
    DO I=1,nd2   !  6 before
       SX(I,I)=MX(I,I)-one
    ENDDO

    DO I=1,ND2
       DIX(I)=FIX(I)-X(I)
    enddo

    CALL matinv(SX,SXI,ND2,6,ier)
    IF(IER==132)  then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=" Inversion failed in FIND_ORBIT_LAYOUT_noda"
       call write_e(333)
       return
    endif

    x=zero
    do i=1,nd2
       do j=1,nd2
          x(i)=sxi(i,j)*dix(j)+x(i)
       enddo
    enddo
    dix=x
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=zero
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo
    !    write(6,*) " Convergence Factor = ",nd2,xdix,deps_tracking
    !    pause 123321
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,a72))'
    write(w_p%c(1),'(a22,g20.14)') " Convergence Factor = ",xdix
    if(xdix.gt.deps_tracking) then
       ite=1
    else
       if(xdix.ge.xdix0.or.xdix<=tiny) then
          ite=0
       else
          ite=1
          xdix0=xdix
       endif
    endif

    if(iteM>=MAX_FIND_ITER)  then
       C_%stable_da=.FALSE.
       IF(iteM==MAX_FIND_ITER+100) THEN
          write(6,*) " Unstable in find_orbit without TPSA"
       ELSE
          write(6,*) " maximum number of iterations in find_orbit without TPSA"
       ENDIF
       ITE=0
    endif

    if(ite.eq.1)  then
       GOTO 3

    endif
    !    FIX(6)=FIX(6)+freq*turns0
    c_%APERTURE_FLAG=APERTURE

  END SUBROUTINE FIND_ORBIT_LAYOUT_noda

  !   backward compatible routine for radiation

  SUBROUTINE find_ENVELOPE(RING,YS,A1,FIX,LOC,STATE)
    ! Finds Envelope with TPSA in state supplied by user which may include parameters
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    TYPE(real_8),INTENT(INOUT)::A1(6)
    TYPE(ENV_8),INTENT(INOUT)::YS(6)
    real(dp),INTENT(INOUT)::FIX(6)
    INTEGER,INTENT(IN) :: LOC
    TYPE(INTERNAL_STATE),INTENT(IN) :: STATE
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    TYPE(NORMALFORM) NORMAL
    logical(lp) APERTURE
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.


    CALL ALLOC(Y)
    CALL ALLOC(ID)
    CALL ALLOC(NORMAL)

    y=6
    y=FIX

    CALL TRACK(RING,Y,LOC,STATE)
    if(.not.check_stable) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       write(w_p%c(1),'(a16)') " find_ENVELOPE 1"
       call write_i
       CALL KILL(NORMAL)
       CALL KILL(ID)
       CALL KILL(Y)

       return
    endif
    normal= y
    y=normal%a1+FIX

    a1=y
    ys=y
    CALL TRACK(RING,Ys,loc,STATE)
    if(.not.check_stable) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       write(w_p%c(1),'(a16)') " find_ENVELOPE 2"
       call write_i
       CALL KILL(NORMAL)
       CALL KILL(ID)
       CALL KILL(Y)

       return
    endif

    y=YS
    id=y
    id=(normal%a1**(-1))*id
    y=id+FIX
    ys=y


    CALL KILL(NORMAL)
    CALL KILL(ID)
    CALL KILL(Y)
    c_%APERTURE_FLAG=APERTURE


  END SUBROUTINE find_ENVELOPE

  ! time tracking
  SUBROUTINE TRACK_time(xT,DT,K)
    ! Tracks a single particle   of the beam for a time DT
    implicit none
    TYPE(INTEGRATION_NODE), POINTER:: T
    TYPE(temporal_probe),INTENT(INOUT):: xT
    REAL(DP), INTENT(IN) :: DT
    REAL(DP) XTT(6),DT0,YL,DT_BEFORE   !X(6),
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    type(element),pointer :: el
    LOGICAL(LP) END_OF_LINE
    INTEGER I
    END_OF_LINE=.FALSE.

    IF(K%TOTALPATH==0) then
       write(6,*) " Must used totalpath in tracking state "
       STOP 451
    endif
    IF(XT%xs%u) RETURN

    !    X=xs%x

    T=>XT%NODE
    !    T%PARENT_FIBRE%MAG=K
    !    T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR
    T%PARENT_FIBRE%MAG%P%CHARGE=>T%PARENT_FIBRE%CHARGE


    DT0=XT%xs%X(6)
    CALL DRIFT_BACK_TO_POSITION(T,XT%Ds,XT%xs%X,k)
    XT%Ds=ZERO
    DT0=XT%xs%X(6)-DT0

    YL=zero

    DO WHILE(DT0<=DT)
       XTT=XT%xs%X
       DT_BEFORE=DT0
       !         WRITE(6,*) " POS ",T%s(1),t%pos_in_fibre
       !         WRITE(6,*) " POS ",T%POS,T%CAS,T%PARENT_FIBRE%MAG%NAME
       ! putting spin and radiation
       CALL TRACK_NODE_PROBE(t,XT%XS,K) !,charge)
       ! no spin and no radiation
       !        CALL TRACK_NODE_SINGLE(t,XT%XS%X,K)  !,CHARGE

       !
       !
       DT0=DT0+(XT%xs%X(6)-XTT(6))
       T=>T%NEXT
       IF(.NOT.ASSOCIATED(T%NEXT)) THEN
          END_OF_LINE=.TRUE.
          EXIT
       ENDIF
    ENDDO

    IF(.NOT.END_OF_LINE) THEN
       IF(DT0/=DT) THEN
          XT%NODE=>T%PREVIOUS
          XT%xs%X=XTT
          DT0=DT-DT_BEFORE
          !           WRITE(6,*) " DT0 ", DT0
          CALL DRIFT_TO_TIME(T,YL,DT0,XT%xs%X,k)
       ELSE
          XT%NODE=>T
       ENDIF
    ELSE


       IF(DT0<DT) THEN
          XT%NODE=>T%PREVIOUS
          XT%xs%X=XTT
          DT0=DT-DT_BEFORE
          CALL DRIFT_TO_TIME(T,YL,DT0,XT%xs%X,k)
       ELSE
          XT%NODE=>T
       ENDIF

    ENDIF


    XT%Ds=yl
    !   xs%x=X
    !    XT%NODE=>T
    if(.not.CHECK_STABLE) then
       !       write(6,*) "unstable "
       CALL RESET_APERTURE_FLAG
       XT%xs%u=.true.
    endif

    call ptc_global_x_p(xt,k)

  END SUBROUTINE TRACK_time

  Subroutine ptc_global_x_p(xt,k)
    implicit none
    TYPE(temporal_probe),INTENT(INOUT):: xT
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer i
    real(dp) p(3),pz,betinv

    if(associated(XT%NODE%ENT)) then
       ! computing global position

       XT%pos=ZERO
       DO I=1,3
          XT%POS(i)=XT%POS(i) + XT%XS%X(1)*XT%NODE%ENT(1,I)     !
          XT%POS(i)=XT%POS(i) + XT%XS%X(3)*XT%NODE%ENT(2,I)     !
          XT%POS(i)=XT%POS(i) + XT%Ds*XT%NODE%PARENT_FIBRE%DIR*XT%NODE%ENT(3,I)     !
       ENDDO
       XT%pos(1:3) = XT%pos(1:3) + XT%NODE%A
       ! computing global momentum
       if(k%time) then
          betinv=one/XT%NODE%PARENT_FIBRE%beta0
       ELSE
          betinv=one
       ENDIF
       P(1)=XT%XS%X(2)*XT%NODE%PARENT_FIBRE%mag%p%p0c
       P(2)=XT%XS%X(4)*XT%NODE%PARENT_FIBRE%mag%p%p0c
       !       pz=XT%NODE%PARENT_FIBRE%mag%p%p0c* &
       !            (one+two*betinv*XT%XS%X(5)+XT%XS%X(5)**2)-xt%pos(4)**2-xt%pos(5)**2
       pz=(one+two*betinv*XT%XS%X(5)+XT%XS%X(5)**2)-XT%XS%X(2)**2-XT%XS%X(4)**2
       pz=XT%NODE%PARENT_FIBRE%mag%p%p0c*root(pz)*XT%NODE%PARENT_FIBRE%DIR
       P(3)=PZ
       DO I=1,3
          XT%POS(i+3)=XT%POS(i+3) + P(1)*XT%NODE%ENT(1,I)     !
          XT%POS(i+3)=XT%POS(i+3) + P(2)*XT%NODE%ENT(2,I)     !
          XT%POS(i+3)=XT%POS(i+3) + P(3)*XT%NODE%ENT(3,I)     !
       ENDDO
    else
       write(6,*) " FILL_SURVEY_DATA_IN_NODE_LAYOUT "
       write(6,*) " was not called, so no survey data on the nodes "
       write(6,*) " program will stop inside TRACK_time "
       stop

    endif

  end Subroutine ptc_global_x_p
  !  type probe
  !   real(dp) x(6)
  !   type(spinor) s
  !   logical u
  !    type(integration_node),pointer :: lost_node
  !  end type probe
  !type TEMPORAL_PROBE
  !    TYPE(probe)  XS
  !    TYPE(INTEGRATION_NODE), POINTER :: NODE
  !    real(DP)  DS,POS(6)
  !END type TEMPORAL_PROBE

  !type TEMPORAL_BEAM
  !    TYPE(TEMPORAL_PROBE), pointer :: TP(:)
  !    real(DP) a(3),ent(3,3),p0c,total_time
  !    integer n
  !    type(integration_node),pointer :: c   ! pointer close to a(3)
  !END type TEMPORAL_BEAM

  SUBROUTINE  alloc_temporal_beam(b,n,p0c) ! fibre i1 to i2
    IMPLICIT NONE
    type(temporal_beam), intent(inout):: b
    integer i,n
    real(dp) p0c


    allocate(b%tp(n))
    b%n=n
    b%a=GLOBAL_origin
    b%ent=global_frame
    b%total_time=zero
    b%p0c=p0c
    nullify(b%c)

    do i=1,n
       call alloc(b%tp(i))
    enddo

  END SUBROUTINE  alloc_temporal_beam

  SUBROUTINE  alloc_temporal_probe(p) ! fibre i1 to i2
    IMPLICIT NONE
    type(temporal_probe), intent(inout):: p
    integer i

    p%xs%u=my_false
    p%xs%x=zero
    p%xs%s%x=zero
    p%ds=zero
    p%pos=zero
    nullify(p%node)
    nullify(p%xs%lost_node)

  end SUBROUTINE  alloc_temporal_probe

  SUBROUTINE  position_temporal_beam(r,b,k) ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(temporal_beam),intent(INOUT) ::  b
    type(integration_node), pointer :: t,tw
    integer i
    real(dp) norm,d1,ds,p(3)
    type(internal_state) :: k

    if(k%totalpath/=1) then
       write(6,*) " One must use total path or time "
       write(6,*) " execution stopped in position_temporal_beam "
       stop 101
    endif
    b%state=k

    norm=mybig
    t=>r%t%start
    b%c=>t
    do i=1,r%t%n

       if(t%cas==case0) then ! 1
          d1=sqrt((t%a(1)-b%a(1))**2+(t%a(2)-b%a(2))**2+(t%a(3)-b%a(3))**2 )

          if(d1<norm) then
             norm=d1
             b%c=>t
          endif
       endif ! 1
       t=>t%next
    enddo

    do i=1,b%n
       b%tp(i)%pos(1:3)=zero
       !      GEO_TRA(A,ENT,D,I)    ! A= A +I D*ENT     I=1,-1
       ! puts the particle in the frame ent around the point a
       call GEO_TRA(b%tp(i)%pos(1:3),b%ent,b%tp(i)%xs%x(1:3),1)
       b%tp(i)%pos(1:3)=b%tp(i)%pos(1:3)+b%a
       !  In the array pos(1:3), the particle is in the global frame of PTC

       ! now we must find the local coordinates
       call locate_temporal_beam(r,b,tw,i)

    enddo


  end SUBROUTINE  position_temporal_beam

  SUBROUTINE  locate_temporal_beam(r,b,tw,j) !
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(temporal_beam),intent(INOUT) ::  b
    type(integration_node), pointer :: t,tw
    type(internal_state) :: k
    real(dp) a(3)
    integer i,j
    real(dp) norm,d1,ds,p(3),da(3),dal(3)


    ! locating the closest integration node
    ! pointing at it using tw
    norm=mybig
    t=> b%c
    a=b%tp(j)%pos(1:3)
    do i=1,r%t%n
       if(t%cas==case0) then ! 1
          d1=sqrt((t%a(1)-a(1))**2+(t%a(2)-a(2))**2+(t%a(3)-a(3))**2 )

          if(d1<norm) then
             norm=d1
             tw=>t
          endif
       endif ! 1
       t=>t%next
    enddo
    !
    !   Expressing the vector from tw to the particle
    ! in the entrance frame of tw
    da=a-tw%a
    call change_basis(DA,global_frame,dal,tw%ent)
    !
    !

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! some gymnastic if behind the integration node
    if(tw%parent_fibre%dir>0) then

       if(dal(3)<zero) then
          tw=>tw%previous
          do while(tw%cas/=case0)
             tw=>tw%previous
          enddo
          da=a-tw%a
          call change_basis(DA,global_frame,dal,tw%ent)
       endif

    else

       if(dal(3)>zero) then
          tw=>tw%previous
          do while(tw%cas/=case0)
             tw=>tw%previous
          enddo
          da=a-tw%a
          call change_basis(DA,b%ent,dal,tw%ent)
       endif
    endif
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    p=dal
    b%tp(j)%ds=p(3)
    b%tp(j)%node=>tw

    call original_p_to_ptc(b,j,p,tw)

  end SUBROUTINE  locate_temporal_beam

  SUBROUTINE original_p_to_ptc(b,j,p,tw) ! fibre i1 to i2
    IMPLICIT NONE
    type(temporal_beam),intent(INOUT) ::  b
    type(integration_node),pointer :: tw
    integer j
    real(dp) dal(3),d1,betinv,p(3)

    b%tp(j)%xs%x(4:6)=b%tp(j)%xs%x(4:6)*b%p0c/tw%parent_fibre%mag%p%p0c
    !    call change_basis(b%tp(j)%xs%x(4:6),global_frame,da,b%ent)
    ! expressing the momentum in the frame of the node tw
    call change_basis(b%tp(j)%xs%x(4:6),b%ent,dal,tw%ent)
!!!!!!!!!!!


    ! final setting of the PTC coordinates
    ! Does not take into acoount vector potential yet
    ! nor x', y' for pancake

    d1=dal(1)**2+dal(2)**2+dal(3)**2
    betinv=one/tw%parent_fibre%beta0

    if(b%state%time) then

       b%tp(j)%xs%x(5)=(d1-one)/ (sqrt(betinv**2-one +d1)+betinv)
       b%tp(j)%xs%x(2)=dal(1)
       b%tp(j)%xs%x(4)=dal(2)
       b%tp(j)%xs%x(1)=p(1)
       b%tp(j)%xs%x(3)=p(2)
       b%tp(j)%xs%x(6)=b%total_time

    else
       b%tp(j)%xs%x(5)=sqrt(d1)-one
       b%tp(j)%xs%x(2)=dal(1)
       b%tp(j)%xs%x(4)=dal(2)
       b%tp(j)%xs%x(1)=p(1)
       b%tp(j)%xs%x(3)=p(2)
       b%tp(j)%xs%x(6)=b%total_time
    endif

  end SUBROUTINE original_p_to_ptc

  SUBROUTINE TRACK_temporal_beam(b,dt,state) ! fibre i1 to i2
    IMPLICIT NONE
    type(temporal_beam),intent(INOUT) ::  b
    TYPE(INTERNAL_STATE), optional:: state
    real(dp) dt
    integer i
    TYPE(INTERNAL_STATE) K
    !BEAM_IN_X(B,I)
    !X_IN_BEAM(B,X,I,DL,T)
    k=b%state
    if(present(state)) k=state
    do i=1,b%n
       if(b%tp(i)%xs%u) cycle
       call TRACK_time(b%tp(i),DT,K)
    enddo
    !          call track_beam(my_ring,TheBeam,getintstate(), pos1=ni, pos2=ni+1)
  end SUBROUTINE TRACK_temporal_beam




  SUBROUTINE TRACK_beam(r,b,k,t, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    real(dp) x(6)
    type(beam),intent(INOUT) ::  b
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    integer i
    type(integration_node),optional, pointer :: t
    !BEAM_IN_X(B,I)
    !X_IN_BEAM(B,X,I,DL,T)

    do i=1,b%n
       if(b%u(i)) cycle
       x=BEAM_IN_X(B,I)
       call track_probe_x(r,x,k,b%u(i),t, fibre1,fibre2,node1,node2)
       call X_IN_BEAM(B,X,I)
    enddo
    !          call track_beam(my_ring,TheBeam,getintstate(), pos1=ni, pos2=ni+1)
  end SUBROUTINE TRACK_beam
end module ptc_spin
