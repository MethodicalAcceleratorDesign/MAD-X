!   B~= (H B_x, H B_y,   B_z)
!   A~= (  A_x,   A_y, H A_z)
!   H = 1 + h x / RHO_0
!
module ptc_spin
  !use orbit_ptc
  !use beam_beam_ptc
  use orbit_ptc
  implicit none
  public
  PRIVATE get_fieldR !,get_field
  PRIVATE get_BfieldR,get_BfieldP,get_Bfield,get_fieldp
 !,GETMULB_TEAPOT
  private B_PANCAkEr,B_PANCAkEp,B_PANCAkE
  PRIVATE DIRECTION_VR,DIRECTION_VP,DIRECTION_V
  PRIVATE  B_PARA_PERPr,B_PARA_PERPp,B_PARA_PERP
  PRIVATE get_omega_spinR,get_omega_spinP,get_omega_spin
  PRIVATE PUSH_SPINR,PUSH_SPINP !,PUSH_SPIN
  PRIVATE TRACK_FRINGE_spin_R,TRACK_FRINGE_spin_P,TRACK_FRINGE_spin
  PRIVATE TRACK_NODE_LAYOUT_FLAG_pr_s12_R,TRACK_NODE_LAYOUT_FLAG_pr_s12_P
  PRIVATE GET_BE_CAVR,GET_BE_CAVP ,GET_BE_CAV
  private rot_spin_xr,rot_spin_xp,rot_spin_zr,rot_spin_zp  !,rot_spin_z,rot_spin_x,
  private rot_spin_yr,rot_spin_yp   !,rot_spin_y 
  private PATCH_SPINR,PATCH_SPINP,PATCH_SPIN,superdrift_SPINR,superdrift_SPINp
  private MIS_SPINR,MIS_SPINP,MIS_SPIN,furman_step
  private DTILT_SPINR,DTILT_SPINP,DTILT_SPIN
  PRIVATE TRACK_SPIN_FRONTR,TRACK_SPIN_FRONTP,TRACK_SPIN_FRONT
  PRIVATE TRACK_SPIN_BACKR,TRACK_SPIN_BACKP,TRACK_SPIN_BACK
  !  private PUSH_SPIN_RAY8
  private radiate_2p,radiate_2r,radiate_2
  private TRACK_NODE_FLAG_probe_R,TRACK_NODE_FLAG_probe_p,TRACK_NODE_LAYOUT_FLAG_spinr_x

  PRIVATE get_Bfield_fringeR,get_Bfield_fringeP,get_Bfield_fringe,TRACK_NODE_LAYOUT_FLAG_spinp_x
  private TRACK_LAYOUT_FLAG_spin12r_x,TRACK_LAYOUT_FLAG_spin12p_x
  PRIVATE TRACK_LAYOUT_FLAG_probe_spin12R,TRACK_LAYOUT_FLAG_probe_spin12P
  private PUSH_SPIN_fake_fringer,PUSH_SPIN_fake_fringep,PUSH_SPIN_fake_fringe
  PRIVATE TRACK_NODE_LAYOUT_FLAG_pr_t12_R,TRACK_NODE_LAYOUT_FLAG_pr_t12_P
  private TRACK_LAYOUT_FLAG_spint12r_x,TRACK_LAYOUT_FLAG_spint12p_x,alloc_temporal_beam
  private alloc_temporal_probe,GET_BZ_fringe,GET_BZ_fringer,GET_BZ_fringep
  private TRACK_rotate_spin_r,TRACK_rotate_spin_p,TRACK_rotate_spin
  private TRACK_FRINGE_multipole_r,TRACK_FRINGE_multipole_p,TRACK_FRINGE_multipole
  private TRACK_wedge_spin_R,TRACK_wedge_spin_p,TRACK_wedge_spin, find_as,find_frac_r,find_n0
  !REAL(DP) :: AG=A_ELECTRON
  REAL(DP) :: bran_init=pi  
  logical :: locate_with_no_cavity = .false.,full_way=.true.
  integer  :: item_min=3,mfdebug
 
  !  INTEGER, PRIVATE :: ISPIN0P=0,ISPIN1P=3


  INTERFACE assignment (=)
     MODULE PROCEDURE equal_temporal
  end  INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE alloc_temporal_probe
     MODULE PROCEDURE alloc_temporal_beam
  END INTERFACE

  INTERFACE TRACK_PROBE2     ! semi private routine
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_s12_R  !#2  ! probe from node i1 to i2 (R,xs,k,I1,I2)
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_s12_P  ! Tracks probe from integer node i1 to i2 in state k
     !     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_R  !#6 USING NODE1 TO NODE2 AS OBJECT
     !     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_P  ! (xs,k,fibre1,fibre2,node1,node2)
  END INTERFACE

  INTERFACE TRACK_PROBE
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_PROBE_spin12R  !#3  ! probe from FIBRE
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_PROBE_spin12P  ! (r,xS,k,fibre1,fibre2,node1,node2) ! integer fibre i1 to i2
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_R  !#6 USING NODE1 TO NODE2 AS OBJECT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_P  ! (xs,k,fibre1,fibre2,node1,node2)
  END INTERFACE

  INTERFACE TRACK_NODE_PROBE                  ! (C,XS,K)  track probe in a node t
     MODULE PROCEDURE TRACK_NODE_FLAG_PROBE_R   ! #1
     MODULE PROCEDURE TRACK_NODE_FLAG_PROBE_p   ! #1p
  END INTERFACE

  INTERFACE TRACK_node_x                !
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spinr_x  !#4 ! TRACK X THROUGH INTEGRATION_NODE T
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spinp_x  !(T,x,k)
  END INTERFACE

  INTERFACE TRACK_node_v
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_spin_v   !#5
  END INTERFACE


  !     call TRACK_NODE_SINGLE(intnode,R%R,my_estate,my_ering%CHARGE)


  INTERFACE TRACK_probe_x
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spin12r_x    !#7
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spin12p_x  ! (r,x,k,u,t, fibre1,fibre2,node1,node2)  integer routine
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spint12r_x
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spint12p_x  !(x,k,u,t, fibre1,fibre2,node1,node2)  pointer routine
  END INTERFACE


  INTERFACE propagate
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spin12r_x    !#7
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spin12p_x  ! (r,x,k,u,t, fibre1,fibre2,node1,node2)  integer routine
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spint12r_x
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_spint12p_x  !(x,k,u,t, fibre1,fibre2,node1,node2)  pointer routine
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_PROBE_spin12R  !#3  ! probe from FIBRE
     MODULE PROCEDURE TRACK_LAYOUT_FLAG_PROBE_spin12P  ! (r,xS,k,fibre1,fibre2,node1,node2) ! integer fibre i1 to i2
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_R  !#6 USING NODE1 TO NODE2 AS OBJECT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_pr_t12_P  ! (xs,k,fibre1,fibre2,node1,node2)
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
     MODULE PROCEDURE DTILT_SPINR
     MODULE PROCEDURE DTILT_SPINP
  END INTERFACE

  INTERFACE TRACK_SPIN_FRONT
     MODULE PROCEDURE TRACK_SPIN_FRONTR
     MODULE PROCEDURE TRACK_SPIN_FRONTP
     !     MODULE PROCEDURE TRACK_SPIN_FRONT_ray8
  END INTERFACE

  INTERFACE TRACK_SPIN_BACK
     MODULE PROCEDURE TRACK_SPIN_BACKR
     MODULE PROCEDURE TRACK_SPIN_BACKP
     !     MODULE PROCEDURE TRACK_SPIN_BACK_RAY8
  END INTERFACE

  INTERFACE GET_BE_CAV
     MODULE PROCEDURE GET_BE_CAVR
     MODULE PROCEDURE GET_BE_CAVP
  END INTERFACE

  INTERFACE TRACK_rotate_spin
     MODULE PROCEDURE TRACK_rotate_spin_r
     MODULE PROCEDURE TRACK_rotate_spin_p
  END INTERFACE

  INTERFACE TRACK_wedge_spin
     MODULE PROCEDURE TRACK_wedge_spin_R
     MODULE PROCEDURE TRACK_wedge_spin_p
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



  INTERFACE TRACK_FRINGE_multipole
     MODULE PROCEDURE TRACK_FRINGE_multipole_r
     MODULE PROCEDURE TRACK_FRINGE_multipole_p
  END INTERFACE

  INTERFACE TRACK_FRINGE_spin
     MODULE PROCEDURE TRACK_FRINGE_spin_R
     MODULE PROCEDURE TRACK_FRINGE_spin_P
  END INTERFACE

  INTERFACE superdrift_SPIN
     MODULE PROCEDURE superdrift_SPINR
     MODULE PROCEDURE superdrift_SPINP
  END INTERFACE


  INTERFACE PUSH_SPIN
     MODULE PROCEDURE PUSH_SPINR
     MODULE PROCEDURE PUSH_SPINP
     !     MODULE PROCEDURE PUSH_SPIN_RAY8
  END INTERFACE

  INTERFACE get_omega_spin
     MODULE PROCEDURE get_omega_spinR
     MODULE PROCEDURE get_omega_spinP
  END INTERFACE

  INTERFACE B_PARA_PERP
     MODULE PROCEDURE B_PARA_PERPr
     MODULE PROCEDURE B_PARA_PERPp
  END INTERFACE

  INTERFACE DIRECTION_V
     MODULE PROCEDURE DIRECTION_Vr
     MODULE PROCEDURE DIRECTION_Vp
  END INTERFACE

  INTERFACE get_field
     MODULE PROCEDURE get_fieldr   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE get_fieldp
  END INTERFACE

  INTERFACE GET_BZ_fringe
     MODULE PROCEDURE GET_BZ_fringer   ! get fringe for multipoles
     MODULE PROCEDURE GET_BZ_fringep
  END INTERFACE

  INTERFACE get_Bfield
     MODULE PROCEDURE get_BfieldR   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE get_BfieldP   ! MID DEFINED AS 1/2 L
  END INTERFACE

  INTERFACE get_Bfield_fringe
     MODULE PROCEDURE get_Bfield_fringeR   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE get_Bfield_fringeP   ! MID DEFINED AS 1/2 L
  END INTERFACE

  INTERFACE PUSH_SPIN_fake_fringe
     MODULE PROCEDURE PUSH_SPIN_fake_fringer   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE PUSH_SPIN_fake_fringep   ! MID DEFINED AS 1/2 L
  END INTERFACE




  INTERFACE B_PANCAkE
     MODULE PROCEDURE B_PANCAkEr   ! MID DEFINED AS 1/2 L
     MODULE PROCEDURE B_PANCAkEp
  END INTERFACE


contains

  subroutine rot_spin_yr(P,ang)
    implicit none
    TYPE(PROBE),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st
    INTEGER I
    type(quaternion) dq
    if(p%use_q) then
     dq%x(0)=COS(ang/2)
     dq%x(2)=sin(ang/2)
     dq%x(1)=0
     dq%x(3)=0
     p%q=dq*p%q
    else
    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO
    endif
  END subroutine rot_spin_yr

  subroutine rot_spin_Xr(P,ang)
    implicit none
    TYPE(PROBE),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st
    INTEGER I
    type(quaternion) dq

    if(p%use_q) then

     dq%x(0)=COS(ang/2)
     dq%x(1)=-sin(ang/2)
     dq%x(2)=0
     dq%x(3)=0
     p%q=dq*p%q

    else
    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(2)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(2)
       P%S(I)%X(2)=ST
    ENDDO
   endif
  END subroutine rot_spin_Xr

  subroutine rot_spin_zr(P,ang)
    implicit none
    TYPE(PROBE),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st
    INTEGER I
    type(quaternion) dq

    if(p%use_q) then
  
     dq%x(0)=COS(ang/2)
     dq%x(3)=-sin(ang/2)
     dq%x(1)=0
     dq%x(2)=0
     p%q=dq*p%q
    else
    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(2)
       P%S(I)%X(2)=CO *P%S(I)%X(2)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO
    endif

  END subroutine rot_spin_zr


  subroutine rot_spin_yp(P,ang)
    implicit none
    type(PROBE_8),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    !type(real_8) co,si,st
    INTEGER I
    type(quaternion_8) dq
    if(p%use_q) then
     call alloc(dq)
     dq%x(0)=COS(ang/2)
     dq%x(2)=sin(ang/2)
     dq%x(1)=0.0_dp
     dq%x(3)=0.0_dp
     p%q=dq*p%q
     call kill(dq)
    else
    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0r,ISPIN1r
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO

    call kill(st)
    endif
  END subroutine rot_spin_yp

  subroutine rot_spin_xp(P,ang)
    implicit none
    type(PROBE_8),INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    INTEGER I
    type(quaternion_8) dq

    if(p%use_q) then
     call alloc(dq)
     dq%x(0)=COS(ang/2)
     dq%x(1)=-sin(ang/2)
     dq%x(2)=0.0_dp
     dq%x(3)=0.0_dp
     p%q=dq*p%q
     call kill(dq)
    else
    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(2)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(2)
       P%S(I)%X(2)=ST
    ENDDO

    call kill(st)
    endif
  END subroutine rot_spin_xp

  subroutine rot_spin_zp(P,ang)
    implicit none
    type(PROBE_8),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    INTEGER I

    type(quaternion_8) dq

    if(p%use_q) then
     call alloc(dq)
     dq%x(0)=COS(ang/2)
     dq%x(3)=-sin(ang/2)
     dq%x(1)=0.0_dp
     dq%x(2)=0.0_dp
     p%q=dq*p%q
     call kill(dq)
    else
    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(2)
       P%S(I)%X(2)=CO *P%S(I)%X(2)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO

    call kill(st)
    endif
  END subroutine rot_spin_zp

  subroutine radiate_2r(c,DS,FAC,X,b2,dlds,before,k,POS)
    use gauss_dis
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    real(dp),INTENT(INOUT) :: X(6)  !,XP(2)
    real(dp), INTENT(IN) :: DS
    REAL(DP), INTENT(IN) :: FAC
    real(dp), intent(in):: B2,dlds
    LOGICAL(LP),intent(in) :: BEFORE
    real(dp)  st,z,av(3),t
    type(internal_state) k

    IF(.NOT.CHECK_STABLE) return
    el=>c%parent_fibre%mag

    if(k%TIME) then
       ST=root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)-1.0_dp
    else
       ST=X(5)
    endif

    ! X(5)=X(5)+B2*FAC*DS
    !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ! X(5)=one/(one/(one+X(5))+CRADF(EL%P)*(one+X(5))*B2*DLDS*FAC*DS)-one
    !        X(5)=X(5)-CRADF(EL%P)*(one+X(5))**3*B2*FAC*DS/SQRT((one+X(5))**2-X(2)**2-X(4)**2)
    if(K%radiation) X(5)=X(5)-CRADF(EL%P)*(1.0_dp+X(5))**3*B2*FAC*DS*DLDS
    if(k%stochastic) then
       !         t=sqrt(12.e0_dp)*(bran(bran_init)-half)
       t=RANF()
       !         t=sqrt(12.d0)*(RANF()-half)
       if(t>0.5_dp) then
          t=1.0_dp
       else
          t=-1.0_dp
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
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp)*root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp)*root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp
          else
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp)*(1.0_dp+X(5))/(1.0_dp+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp)*(1.0_dp+X(5))/(1.0_dp+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp
          endif
       ELSEif(el%kind==kind22) then

          IF(EL%HE22%P%DIR==1) THEN
             Z= pos*el%l/el%p%nst
          ELSE
             Z=EL%L-pos*el%l/el%p%nst
          ENDIF
          CALL compute_f4(EL%he22,X,Z,A=AV)
          if(k%TIME) then
             X(2)=(X(2)+EL%P%CHARGE*AV(1)) *root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2)) *root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          else
             X(2)=(X(2)+EL%P%CHARGE*AV(1))*(1.0_dp+X(5))/(1.0_dp+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2))*(1.0_dp+X(5))/(1.0_dp+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          endif


       ELSE
          if(k%TIME) then
             X(2)=X(2)*root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(4)=X(4)*root(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
          else
             X(2)=X(2)*(1.0_dp+X(5))/(1.0_dp+ST)
             X(4)=X(4)*(1.0_dp+X(5))/(1.0_dp+ST)
          endif
       ENDIF
    endif

    !       X(2)=X_MEC(2)*(one+X(5))/(one+X5)-EL%B_SOL*EL%P%CHARGE*X(3)/two
    !       X(4)=X_MEC(4)*(one+X(5))/(one+X5)+EL%B_SOL*EL%P%CHARGE*X(1)/two


  end subroutine radiate_2r


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
    real(dp) b30,x1,x3,denf
    type(damap) xpmap
    integer i,j
    type(internal_state) k

    IF(.NOT.CHECK_STABLE) return
    el=>c%parent_fibre%magp
    if(.not.before.and.k%envelope) then

       denf=(1.0_dp+x(5))**5/SQRT((1.0_dp+X(5))**2-Xp(1)**2-Xp(2)**2)
       b30=b2
       b30=b30**1.5e0_dp
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
             X1=(xpmap%v(i)).sub.'000010'   ! Still works if BMAD units are used because xpmax**(-1) is needed!!!
             X3=(xpmap%v(j)).sub.'000010'
             E_IJ(i,j)=E_IJ(i,j)+denf*x1*x3 ! In a code internally using BMAD units '000001' is needed!!!
          enddo
       enddo
       if(compute_stoch_kick) c%delta_rad_out=root(denf)
       call kill(xpmap)
    endif


    call alloc(st)

    if(k%TIME) then
       ST=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)-1.0_dp
    else
       ST=X(5)
    endif

    ! X(5)=X(5)+B2*FAC*DS
    !   X(5)=one/(one/(one+X(5))-B2*FAC*DS)-one
    !   X(5)=one/(one/(one+X(5))+CRADF(EL%P)*(one+X(5))*B2*DLDS*FAC*DS)-one
    !          X(5)=X(5)-CRADF(EL%P)*(one+X(5))**3*B2*FAC*DS/SQRT((one+X(5))**2-X(2)**2-X(4)**2)
  if(K%radiation)  X(5)=X(5)-CRADF(EL%P)*(1.0_dp+X(5))**3*B2*FAC*DS*DLDS


    if(el%kind/=kindpa) then
       IF(ASSOCIATED(EL%B_SOL)) THEN
          if(k%TIME) then
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp)*SQRT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp)*SQRT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp
          else
             X(2)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp)*(1.0_dp+X(5))/(1.0_dp+ST)
             X(2)=X(2)-EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp
             X(4)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp)*(1.0_dp+X(5))/(1.0_dp+ST)
             X(4)=X(4)+EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp
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
             X(2)=(X(2)+EL%P%CHARGE*AV(1)) *sqrt(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2)) *sqrt(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          else
             X(2)=(X(2)+EL%P%CHARGE*AV(1))*(1.0_dp+X(5))/(1.0_dp+ST)
             X(2)=X(2)-EL%P%CHARGE*AV(1)
             X(4)=(X(4)-EL%P%CHARGE*AV(2))*(1.0_dp+X(5))/(1.0_dp+ST)
             X(4)=X(4)+EL%P%CHARGE*AV(2)
          endif
          call kill(av,3)
          call kill(z)
       ELSE
          if(k%TIME) then
             X(2)=X(2)*SQRT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
             X(4)=X(4)*SQRT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)/(1.0_dp+ST)
          else
             X(2)=X(2)*(1.0_dp+X(5))/(1.0_dp+ST)
             X(4)=X(4)*(1.0_dp+X(5))/(1.0_dp+ST)
          endif
       ENDIF
    endif

    call kill(st)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(before.and.k%envelope) then
       denf=(1.0_dp+x(5))**5/SQRT((1.0_dp+X(5))**2-Xp(1)**2-Xp(2)**2)

       b30=b2
       b30=b30**1.5e0_dp
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
       if(compute_stoch_kick) c%delta_rad_in=root(denf)
       call kill(xpmap)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine radiate_2p

  !  subroutine PUSH_SPIN_fake_fringer(c,p,before,k,POS)
  subroutine PUSH_SPIN_fake_fringer(c,p,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    !    real(dp),INTENT(INOUT) :: X(6),S(3)
    type(probe),INTENT(INOUT) :: p

    real(dp) OM(3),CO(3),SI(3),B2,XP(2)
    real(dp) ST,dlds,norm,stheta
    !    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k
    type(quaternion) dq
    INTEGER I

    IF(.NOT.CHECK_STABLE) return

    if(.not.((k%radiation.or.k%envelope).or.k%SPIN)) return
    el=>c%parent_fibre%mag
    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    if(EL%kind<=kind1) return



    CALL get_omega_spin(c,OM,B2,dlds,XP,p%X,POS,k)
    !if(k%radiation.AND.BEFORE) then
    !if(el%p%radiation.AND.BEFORE) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif

    if(k%SPIN) then
    if(p%use_q) then
        om=OM/2.0_dp
 
        norm=sqrt(om(1)**2+om(2)**2+om(3)**2)
        if(norm>0) then
        stheta=sin(norm)
        dq%x(0)=cos(norm)
        dq%x(1)=stheta*om(1)/norm
        dq%x(2)=stheta*om(2)/norm
        dq%x(3)=stheta*om(3)/norm
        p%q=dq*p%q
      endif
else
       CO(1)=COS(OM(1)/2.0_dp)
       SI(1)=SIN(OM(1)/2.0_dp)
       CO(2)=COS(OM(2)/2.0_dp)
       SI(2)=SIN(OM(2)/2.0_dp)
       CO(3)=COS(OM(3))
       SI(3)=SIN(OM(3))

       DO I=ISPIN0R,ISPIN1R
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%S(I)%X(2)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(3)*p%S(I)%X(1)-SI(3)*p%S(I)%X(2)
          p%S(I)%X(2)= CO(3)*p%S(I)%X(2)+SI(3)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%s(I)%X(2)=ST
       ENDDO
    endif
endif
    !if(k%radiation.AND.(.NOT.BEFORE)) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif



  END subroutine PUSH_SPIN_fake_fringer

  !  subroutine PUSH_SPIN_fake_fringep(c,p,before,k,POS)
  subroutine PUSH_SPIN_fake_fringep(c,p,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTp), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(probe_8),INTENT(INOUT) :: p
    !    TYPE(REAL_8),INTENT(INOUT) :: X(6),S(3)

    TYPE(REAL_8) OM(3),CO(3),SI(3),B2,XP(2)
    TYPE(REAL_8) ST,dlds,norm,stheta
    !    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k
    INTEGER I
    type(quaternion_8) dq

    IF(.NOT.CHECK_STABLE) return
    if(.not.((k%radiation.or.k%envelope).or.k%SPIN)) return
    el=>c%parent_fibre%magp
    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    if(EL%kind<=kind1) return

    IF(K%PARA_IN ) KNOB=.TRUE.


    CALL ALLOC(OM,3)
    CALL ALLOC(CO,3)
    CALL ALLOC(SI,3)
    CALL ALLOC(XP,2)
    CALL ALLOC(ST,B2,dlds)

    CALL get_omega_spin(c,OM,B2,dlds,XP,p%X,POS,k)
    !if(k%radiation.AND.BEFORE) then
    !if(el%p%radiation.AND.BEFORE) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif

    if(k%SPIN) then
    if(p%use_q) then
     call alloc(dq)
     call alloc(norm,stheta)
       do i=1,3
        om(i)=OM(i)/2.0_dp
       enddo
    
      norm=om(1)**2+om(2)**2+om(3)**2

        stheta=sin_quaternion(norm)
        dq%x(0)=cos_quaternion(norm)
        dq%x(1)=stheta*om(1)
        dq%x(2)=stheta*om(2)
        dq%x(3)=stheta*om(3)
        p%q=dq*p%q

       call kill(norm,stheta)
       call kill(dq)
   else
       CO(1)=COS(OM(1)/2.0_dp)
       SI(1)=SIN(OM(1)/2.0_dp)
       CO(2)=COS(OM(2)/2.0_dp)
       SI(2)=SIN(OM(2)/2.0_dp)
       CO(3)=COS(OM(3))
       SI(3)=SIN(OM(3))

       !       ST=   CO(1)*P%S%X(2)-SI(1)*P%S%X(3)
       !       P%S%X(3)= CO(1)*P%S%X(3)+SI(1)*P%S%X(2)
       !       P%S%X(2)=ST
       !       ST=  CO(2)*P%S%X(1)+SI(2)*P%S%X(3)
       !       P%S%X(3)=CO(2)*P%S%X(3)-SI(2)*P%S%X(1)
       !       P%S%X(1)=ST
       !       ST=   CO(3)*P%S%X(1)-SI(3)*P%S%X(2)
       !       P%S%X(2)= CO(3)*P%S%X(2)+SI(3)*P%S%X(1)
       !       P%S%X(1)=ST
       !       ST=  CO(2)*P%S%X(1)+SI(2)*P%S%X(3)
       !       P%S%X(3)=CO(2)*P%S%X(3)-SI(2)*P%S%X(1)
       !       P%S%X(1)=ST
       !       ST=   CO(1)*P%S%X(2)-SI(1)*P%S%X(3)
       !       P%S%X(3)= CO(1)*P%S%X(3)+SI(1)*P%S%X(2)
       !       P%S%X(2)=ST

       DO I=ISPIN0R,ISPIN1R
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%S(I)%X(2)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(3)*p%S(I)%X(1)-SI(3)*p%S(I)%X(2)
          p%S(I)%X(2)= CO(3)*p%S(I)%X(2)+SI(3)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%s(I)%X(2)=ST
       ENDDO

    endif
endif
    !if(k%radiation.AND.(.NOT.BEFORE)) then
    ! call radiate_2(EL,DS,FAC,X,E_IJ,b2,dlds,XP,before,k,POS)
    !endif

    CALL KILL(OM,3)
    CALL KILL(CO,3)
    CALL KILL(SI,3)
    CALL KILL(XP,2)
    CALL KILL(ST,B2,dlds)
    KNOB=.false.
  END subroutine PUSH_SPIN_fake_fringep

  !  subroutine PUSH_SPINR(c,DS,FAC,S,X,before,k,POS)
  subroutine PUSH_SPINR(c,DS,FAC,P,before,k,POS)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(PROBE), INTENT(INOUT) :: P
    !    REAL(DP),INTENT(INOUT) :: X(6),S(3)
    REAL(DP), INTENT(IN) :: DS,FAC
    REAL(DP) OM(3),CO(3),SI(3),B2,XP(2)
    REAL(DP) ST,dlds,norm,stheta
    type(quaternion) dq,mulq
    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k
    INTEGER I

    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    if(.not.(k%radiation.or.k%SPIN.or.k%envelope)) return
    IF(.NOT.CHECK_STABLE) return
    el=>c%parent_fibre%mag
    if(EL%kind<=kind1) return    ! should I prevent monitor here??? instead of xp=Px,y in get_omega_spin
    !    if(EL%kind>=kind11.and.EL%kind<=kind14) return    ! should I prevent monitor here??? instead of xp=Px,y in get_omega_spin
    !    if(EL%kind>=kind18.and.EL%kind<=kind19) return    ! should I prevent monitor here??? instead of xp=Px,y in get_omega_spin

    CALL get_omega_spin(c,OM,B2,dlds,XP,P%X,POS,k)
    if((k%radiation.or.k%envelope).AND.BEFORE) then
       !if(el%p%radiation.AND.BEFORE) then
       !       call radiate_2(c,DS,FAC,P%X,b2,dlds,XP,before,k,POS)
       call radiate_2(c,DS,FAC,P%X,b2,dlds,before,k,POS)
    endif
   if(k%spin) then
    if(p%use_q) then
      if(EL%kind/=kind3) then
        om=FAC*DS*OM/2.0_dp
      else
        om=FAC*OM/2.0_dp
      endif
        norm=sqrt(om(1)**2+om(2)**2+om(3)**2)
        if(norm>0) then
        stheta=sin(norm)
        dq%x(0)=cos(norm)
        dq%x(1)=stheta*om(1)/norm
        dq%x(2)=stheta*om(2)/norm
        dq%x(3)=stheta*om(3)/norm
        p%q=dq*p%q
 !         mulq%x=0.0_dp!

!          mulq%x(1)=dq%x(1)*p%q%x(1)-dq%x(2)*p%q%x(2)-dq%x(3)*p%q%x(3)-dq%x(4)*p%q%x(4)

 !        mulq%x(2)= dq%x(3)*p%q%x(4)-dq%x(4)*p%q%x(3)+ dq%x(1)*p%q%x(2)+ dq%x(2)*p%q%x(1)
 !        mulq%x(3)= dq%x(4)*p%q%x(2)-dq%x(2)*p%q%x(4)+ dq%x(1)*p%q%x(3)+ dq%x(3)*p%q%x(1)
!          p%q%x(4)= dq%x(2)*p%q%x(3)-dq%x(3)*p%q%x(2)+ dq%x(1)*p%q%x(4)+ dq%x(4)*p%q%x(1)
 !         p%q%x(1)=mulq%x(1)
  !        p%q%x(2)=mulq%x(2)
   !       p%q%x(3)=mulq%x(3)
      endif

    else
     if(EL%kind/=kind3) then
       CO(1)=COS(FAC*DS*OM(1)/2.0_dp)
       SI(1)=SIN(FAC*DS*OM(1)/2.0_dp)
       CO(2)=COS(FAC*DS*OM(2)/2.0_dp)
       SI(2)=SIN(FAC*DS*OM(2)/2.0_dp)
       CO(3)=COS(FAC*DS*OM(3))
       SI(3)=SIN(FAC*DS*OM(3))
    else
       CO(1)=COS(FAC*OM(1)/2.0_dp)
       SI(1)=SIN(FAC*OM(1)/2.0_dp)
       CO(2)=COS(FAC*OM(2)/2.0_dp)
       SI(2)=SIN(FAC*OM(2)/2.0_dp)
       CO(3)=COS(FAC*OM(3))
       SI(3)=SIN(FAC*OM(3))
    endif

       DO I=ISPIN0R,ISPIN1R
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%S(I)%X(2)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(3)*p%S(I)%X(1)-SI(3)*p%S(I)%X(2)
          p%S(I)%X(2)= CO(3)*p%S(I)%X(2)+SI(3)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%s(I)%X(2)=ST
       ENDDO
    endif
   endif
    !if(el%p%radiation.AND.(.NOT.BEFORE)) then
    if((k%radiation.or.k%envelope).AND.(.NOT.BEFORE)) then
       !       call radiate_2(c,DS,FAC,P%X,b2,dlds,XP,before,k,POS)
       call radiate_2(c,DS,FAC,P%X,b2,dlds,before,k,POS)
    endif

  END subroutine PUSH_SPINR

  subroutine PUSH_SPINP(c,DS,FAC,P,before,k,POS) !,E_IJ
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(PROBE_8),INTENT(INOUT) ::P
    !    TYPE(REAL_8),INTENT(INOUT) :: X(6),S(3)
    !    real(dp),INTENT(INOUT) :: E_IJ(6,6)
    TYPE(REAL_8), INTENT(INout) :: DS
    REAL(DP), INTENT(IN) :: FAC
    TYPE(REAL_8) OM(3),CO(3),SI(3),B2,XP(2)
    TYPE(REAL_8) ST,dlds,norm,stheta
    type(quaternion_8) dq
    LOGICAL(LP),intent(in) :: BEFORE
    type(internal_state) k
    INTEGER I
    if(.not.((k%radiation.or.k%envelope).or.k%SPIN)) return
    IF(.NOT.CHECK_STABLE) return
    !if(.not.(el%p%radiation.or.EL%P%SPIN)) return
    el=>c%parent_fibre%magp
    if(EL%kind<=kind1) return    ! should I prevent monitor here??? instead of xp=Px,y in get_omega_spin
    !    if(EL%kind>=kind11.and.EL%kind<=kind14) return    ! should I prevent monitor here??? instead of xp=Px,y in get_omega_spin
    !    if(EL%kind>=kind18.and.EL%kind<=kind19) return    ! should I prevent monitor here??? instead of xp=Px,y in get_omega_spin

    CALL ALLOC(OM,3)
    CALL ALLOC(CO,3)
    CALL ALLOC(SI,3)
    CALL ALLOC(XP,2)
    CALL ALLOC(ST,B2,dlds)


    IF(K%PARA_IN ) KNOB=.TRUE.
    CALL get_omega_spin(c,OM,B2,dlds,XP,P%X,POS,k)
    if((k%radiation.or.k%envelope).AND.BEFORE) then
       !if(el%p%radiation.AND.BEFORE) then
       call radiate_2(c,DS,FAC,P%X,P%E_IJ,b2,dlds,XP,before,k,POS)
       !       call radiate_2(c,DS,FAC,P%X,E_IJ,b2,dlds,XP,before,k,POS)

    endif

   if(k%spin) then
    if(p%use_q) then
     call alloc(dq)
     call alloc(norm,stheta)
      if(EL%kind/=kind3) then
       do i=1,3
        om(i)=FAC*DS*OM(i)/2.0_dp
       enddo
      else
       do i=1,3
        om(i)=FAC*OM(i)/2.0_dp
       enddo
      endif
      norm=om(1)**2+om(2)**2+om(3)**2

        stheta=sin_quaternion(norm)
        dq%x(0)=cos_quaternion(norm)
        dq%x(1)=stheta*om(1)
        dq%x(2)=stheta*om(2)
        dq%x(3)=stheta*om(3)
        p%q=dq*p%q

       call kill(norm,stheta)
       call kill(dq)

    else
     if(EL%kind/=kind3) then
       CO(1)=COS(FAC*DS*OM(1)/2.0_dp)
       SI(1)=SIN(FAC*DS*OM(1)/2.0_dp)
       CO(2)=COS(FAC*DS*OM(2)/2.0_dp)
       SI(2)=SIN(FAC*DS*OM(2)/2.0_dp)
       CO(3)=COS(FAC*DS*OM(3))
       SI(3)=SIN(FAC*DS*OM(3))
    else
       CO(1)=COS(FAC*OM(1)/2.0_dp)
       SI(1)=SIN(FAC*OM(1)/2.0_dp)
       CO(2)=COS(FAC*OM(2)/2.0_dp)
       SI(2)=SIN(FAC*OM(2)/2.0_dp)
       CO(3)=COS(FAC*OM(3))
       SI(3)=SIN(FAC*OM(3))
    endif

       DO I=ISPIN0R,ISPIN1R
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%S(I)%X(2)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(3)*p%S(I)%X(1)-SI(3)*p%S(I)%X(2)
          p%S(I)%X(2)= CO(3)*p%S(I)%X(2)+SI(3)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=  CO(2)*p%S(I)%X(1)+SI(2)*p%S(I)%X(3)
          p%S(I)%X(3)=CO(2)*p%S(I)%X(3)-SI(2)*p%S(I)%X(1)
          p%S(I)%X(1)=ST
          ST=   CO(1)*p%S(I)%X(2)-SI(1)*p%S(I)%X(3)
          p%S(I)%X(3)= CO(1)*p%S(I)%X(3)+SI(1)*p%S(I)%X(2)
          p%s(I)%X(2)=ST
       ENDDO
      endif
    endif
    if((k%radiation.or.k%envelope).AND.(.NOT.BEFORE)) then
       !if(el%p%radiation.AND.(.NOT.BEFORE)) then
       call radiate_2(c,DS,FAC,P%X,P%E_IJ,b2,dlds,XP,before,k,POS)
       !       call radiate_2(c,DS,FAC,P%X,E_IJ,b2,dlds,XP,before,k,POS)
    endif

    CALL KILL(OM,3)
    CALL KILL(CO,3)
    CALL KILL(SI,3)
    CALL KILL(XP,2)
    CALL KILL(ST,B2,dlds)
    knob=.false.
  END subroutine PUSH_SPINP

  subroutine get_omega_spinr(c,OM,B2,dlds,XP,X,POS,k)
    implicit none
    TYPE(integration_node), POINTER::c
    TYPE(ELEMENT), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    REAL(DP),INTENT(INOUT) :: X(6),OM(3),B2,XP(2),DLDS
    REAL(DP)  B(3),E(3),BPA(3),BPE(3),D1,D2,GAMMA,EB(3),EFD(3),beta,ed(3)
    REAL(DP) BETA0,GAMMA0I,XPA(2),phi,del,z
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(.NOT.CHECK_STABLE) return
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
    OM(2)=0.0_dp
    EB=0.0_dp
    BPA=0.0_dp
    BPE=0.0_dp
    B=0.0_dp
    E=0.0_dp
    EFD=0.0_dp
    phi=0.0_dp

    xp(1)=x(2)
    xp(2)=x(4)   !  to prevent a crash in monitors, etc... CERN june 2010
    dlds=0.0_dp
    del=x(5)
    CALL get_field(EL,B,E,phi,X,k,POS)

    SELECT CASE(EL%KIND) 
    case(KIND2,kind3,kind5:kind7,kindwiggler) ! Straight for all practical purposes
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
       if(pos>=0) OM(2)=p%dir*P%b0   ! not fake fringe
    case(KIND4) ! CAVITY
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-X(2)**2-X(4)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
    case(KIND16:kind17,KIND20)
       CALL B_PARA_PERP(k,el,0,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case(kind10)     ! TEAPOT real curvilinear
       x(5)=x(5)-phi*EL%P%CHARGE
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD,pos=POS)
       x(5)=x(5)+phi*EL%P%CHARGE

           DEL=x(5)-phi*EL%P%CHARGE
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
       if(pos>=0) OM(2)=p%dir*P%b0   ! not fake fringe
    case(kindabell)     ! TEAPOT real curvilinear

       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD,pos=POS)
       CALL get_z_ab(EL%ab,POS,z)
       call B_E_FIELD(EL%ab,x,Z,psie_in=phi)
           DEL=x(5)+phi*EL%P%CHARGE
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
       if(pos>=0) OM(2)=p%dir*P%b0   ! not fake fringe
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       if(k%time) then
          beta0=p%beta0;GAMMA0I=p%GAMMA0I;
       else
          beta0=1.0_dp;GAMMA0I=0.0_dp;
       endif
       d1=root(x(2)**2+x(4)**2+(1.0_dp+el%pa%hc*x(1))**2)
       d2=1.0_dp+2.0_dp*del/beta0+del**2
       d2=gamma0I/beta0/d2
       DLDS=root((1.0_dp+d2**2))*d1/(1.0_dp/BETA0+del)
       OM(2)=p%dir*el%pa%hc
    CASE(KIND21)     ! travelling wave cavity
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF

    case(KIND22)
       CALL B_PARA_PERP(k,EL,0,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case default
       OM(1)=0.0_dp
       OM(2)=0.0_dp
       OM(3)=0.0_dp
    END SELECT

    IF(.not.k%TIME) THEN
      del=(2*del+del**2)/(root(1.0_dp/p%beta0**2+2.0_dp*del+del**2)+1.0_dp/p%beta0)
    endif

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( 1.0_dp/P%BETA0 + del )

    OM(1)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(1) + (1.0_dp+p%AG)*BPA(1) )
    OM(2)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(2) + (1.0_dp+p%AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(3) + (1.0_dp+p%AG)*BPA(3) )


    beta=root(1.0_dp+2.0_dp*del/p%beta0+del**2)/(1.0_dp/P%BETA0 + del)  ! replaced


    DO I=1,3
       OM(I)=OM(I)+a_spin_scale*DLDS*beta*gamma*(p%AG+1.0_dp/(1.0_dp+GAMMA))*EB(I)
    ENDDO

   beta=root(1.0_dp+2.0_dp*x(5)/p%beta0+x(5)**2)*P%BETA0/P%GAMMA0I  ! replace  this

    om(1)=-DLDS*0.5_dp*e_muon*beta*(ed(2)*BPE(3)-ed(3)*BPE(2)) +  om(1)
    om(2)=-DLDS*0.5_dp*e_muon*beta*(ed(3)*BPE(1)-ed(1)*BPE(3)) +  om(2)
    om(3)=-DLDS*0.5_dp*e_muon*beta*(ed(1)*BPE(2)-ed(2)*BPE(1)) +  om(3)

    DO I=1,3
       OM(I)=OM(I)-DLDS*0.5_dp*e_muon*(GAMMA*E(I)+(1-GAMMA)*EFD(I))
    ENDDO

    !IF(.not.k%TIME) THEN
    !   del=(2.0_dp*del/p%beta0+del**2)/(sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)+1.0_dp)
    !endif

    if((k%radiation.or.k%envelope)) then
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
    TYPE(REAL_8)  B(3),E(3),BPA(3),BPE(3),DLDS,D1,D2,GAMMA,EB(3),efd(3),XPA(2),ed(3),beta,phi,del,z
    REAL(DP) BETA0,GAMMA0I
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    !  TESTBUG SATEESH
    !       TYPE(REAL_8) XS(6)
    !       TYPE(DAMAP) ID
    !       REAL(DP) CLO(6)

    !     CALL ALLOC(XS)
    !     CALL ALLOC(ID)
    !     CLO=X
    !     XS=X
    !     ID=1
    !     X=CLO+ID

    IF(.NOT.CHECK_STABLE) return

    CALL ALLOC(B,3)
    CALL ALLOC(E,3)
    CALL ALLOC(Ed,3)
    CALL ALLOC(efd,3)
    CALL ALLOC(beta,del,z)
    CALL ALLOC(EB,3)
    CALL ALLOC(BPA,3)
    CALL ALLOC(BPE,3)
    CALL ALLOC(XPA,2)
    CALL ALLOC(D1,D2,GAMMA,phi)

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
    OM(2)=0.0_dp
    xp(1)=x(2)
    xp(2)=x(4)   !  to prevent a crash in monitors, etc... CERN june 2010
    dlds=0.0_dp
    del=x(5)

    CALL get_field(EL,B,E,phi,X,k,POS)

    SELECT CASE(EL%KIND) 
    case(KIND2,kind3,kind5:kind7,kindwiggler) ! Straight for all practical purposes
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/SQRT(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/SQRT((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
       if(pos>=0) OM(2)=p%dir*P%b0   ! not fake fringe
    case(KIND4) ! CAVITY
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/SQRT(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/SQRT((1.0_dp+del)**2-X(2)**2-X(4)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
    case(KIND16:kind17,KIND20)
       CALL B_PARA_PERP(k,el,0,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/SQRT(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=1.0_dp/SQRT((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case(kind10)     ! TEAPOT real curvilinear
       x(5)=x(5)-phi*EL%P%CHARGE
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD,pos=POS)
       x(5)=x(5)+phi*EL%P%CHARGE

           DEL=x(5)-phi*EL%P%CHARGE
       IF(k%TIME) THEN
          DLDS=1.0_dp/SQRT(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/SQRT((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
       if(pos>=0) OM(2)=p%dir*P%b0   ! not fake fringe
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       if(k%time) then
          beta0=p%beta0;GAMMA0I=p%GAMMA0I;
       else
          beta0=1.0_dp;GAMMA0I=0.0_dp;
       endif
       d1=sqrt(x(2)**2+x(4)**2+(1.0_dp+el%pa%hc*x(1))**2)
       d2=1.0_dp+2.0_dp*del/beta0+del**2
       d2=gamma0I/beta0/d2
       DLDS=sqrt((1.0_dp+d2**2))*d1/(1.0_dp/BETA0+del)
       OM(2)=p%dir*el%pa%hc
    CASE(KIND21)     ! travelling wave cavity
       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/sqrt(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/sqrt((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
    case(KIND22)
       CALL B_PARA_PERP(k,EL,0,X,B,BPA,BPE,XP,XPA,ed,pos=POS)
       IF(k%TIME) THEN
          DLDS=1.0_dp/SQRT(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)
       ELSE
          DLDS=1.0_dp/SQRT((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)
       ENDIF
    case(kindabell)     ! TEAPOT real curvilinear

       CALL B_PARA_PERP(k,EL,1,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD,pos=POS)
       CALL get_z_ab(EL%ab,POS,z)
       call B_E_FIELD(EL%ab,x,Z,psie_in=phi)
           DEL=x(5)+phi*EL%P%CHARGE
       IF(k%TIME) THEN
          DLDS=1.0_dp/sqrt(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/sqrt((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF
       if(pos>=0) OM(2)=p%dir*P%b0   ! not fake fringe
    case default
       OM(1)=0.0_dp
       OM(2)=0.0_dp
       OM(3)=0.0_dp
    END SELECT

    IF(.not.k%TIME) THEN
      del=(2*del+del**2)/(sqrt(1.0_dp/p%beta0**2+2.0_dp*del+del**2)+1.0_dp/p%beta0)
    endif

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( 1.0_dp/P%BETA0 + del )

    OM(1)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(1) + (1.0_dp+p%AG)*BPA(1) )
    OM(2)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(2) + (1.0_dp+p%AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(3) + (1.0_dp+p%AG)*BPA(3) )


    beta=sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)/(1.0_dp/P%BETA0 + del)  ! replaced


    DO I=1,3
       OM(I)=OM(I)+a_spin_scale*DLDS*beta*gamma*(p%AG+1.0_dp/(1.0_dp+GAMMA))*EB(I)
    ENDDO

    e_muon_scale%r=e_muon
    beta=sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)*P%BETA0/P%GAMMA0I

    om(1)=-DLDS*0.5_dp*e_muon_scale*beta*(ed(2)*BPE(3)-ed(3)*BPE(2)) +  om(1)
    om(2)=-DLDS*0.5_dp*e_muon_scale*beta*(ed(3)*BPE(1)-ed(1)*BPE(3)) +  om(2)
    om(3)=-DLDS*0.5_dp*e_muon_scale*beta*(ed(1)*BPE(2)-ed(2)*BPE(1)) +  om(3)

    DO I=1,3
       OM(I)=OM(I)-DLDS*0.5_dp*e_muon_scale*(GAMMA*E(I)+(1-GAMMA)*EFD(I))
    ENDDO

    !IF(.not.k%TIME) THEN
    !   del=(2.0_dp*del/p%beta0+del**2)/(sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)+1.0_dp)
    !endif

    if((k%radiation.or.k%envelope)) then
       !      if(P%RADIATION) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF


    CALL KILL(B,3)
    CALL KILL(E,3)
    CALL KILL(EB,3)
    CALL KILL(BPA,3)
    CALL KILL(BPE,3)
    CALL KILL(D1,D2,GAMMA,phi,z)
    CALL KILL(XPA,2)
    CALL KILL(Ed,3)
    CALL KILL(efd,3)
    CALL KILL(beta,del)

    !  TESTBUG SATEESH
    !    CALL KILL(XS)
    !     CALL KILL(ID)

  end subroutine get_omega_spinp

  subroutine get_fieldr(EL,B,E,phi,X,k,POS)
    implicit none
    TYPE(ELEMENT), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    REAL(DP),INTENT(INOUT) :: B(3),E(3),phi
    REAL(DP),INTENT(INOUT) :: X(6)
    REAL(DP) Z,VM,a(3),ad(3)
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    B=0.0_dp
    E=0.0_dp
    a=0
    ad=0
    phi=0.0_dp
    SELECT CASE(EL%KIND)
    case(KIND2,kind3,kind5:kind7,KIND16:kind17,KIND20) ! Straight for all practical purposes

       if(present(pos)) then
          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,E,X,pos,k)   ! fringe effect
          ELSE
             call get_Bfield(EL,B,X)   ! fringe effect
          ENDIF
       else
          CALL get_Bfield(EL,B,X)
       endif
    case(kind10)     ! TEAPOT real curvilinear
       if(present(pos)) then
          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,E,X,pos,k)   ! fringe effect
          ELSE
!             if(EL%TP10%electric) then
              call GETELECTRIC(EL%TP10,E,phi,B,VM,X); E(3)=0.d0;
!             else
!              CALL GETMULB_TEAPOT(EL%TP10,B,VM,X)
!             endif
          ENDIF
       else
!             if(EL%TP10%electric) then
              call GETELECTRIC(EL%TP10,E,phi,B,VM,X); E(3)=0.d0;
!             else
!              CALL GETMULB_TEAPOT(EL%TP10,B,VM,X)
!             endif
       endif
    case(KINDPA)     ! fitted field for real magnet
       CALL B_PANCAkE(EL%PA,B,X,POS)
    case(KINDWIGGLER)
       CALL get_z_wi(EL%wi,POS,z)
       CALL B_FIELD(EL%wi,Z,X,B)

    CASE(KIND4)      ! Pill box cavity
      if(EL%C4%n_bessel/=-1) then
        CALL GET_BE_CAV(EL%C4,B,E,X,k)
      else
       IF(EL%c4%P%DIR==1) THEN
          Z= pos*el%l/el%p%nst
       ELSE
          Z=EL%L-pos*el%l/el%p%nst
       ENDIF
       call  Abmad_TRANS(EL%C4,Z,X,k,A,AD,B,E)
      endif

    CASE(KIND21)     ! travelling wave cavity
        call get_z_cav(EL%cav21,pos,z)

       call A_TRANS(EL%cav21,Z,X,k,A,AD,B,E)

    CASE(kindabell)     ! travelling wave cavity

       CALL get_z_ab(EL%ab,POS,z)


          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,E,X,pos,k)   ! fringe effect
          ELSE
              call B_E_FIELD(EL%ab,x,Z,E_in=E,B_in=b)
          ENDIF


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

  subroutine get_fieldp(EL,B,E,phi,X,k,POS)
    implicit none
    TYPE(ELEMENTP), POINTER::EL
    INTEGER,OPTIONAL,INTENT(IN) ::POS
    TYPE(REAL_8),INTENT(INOUT) :: B(3),E(3),phi
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    INTEGER I
    TYPE(REAL_8) z,VM,ad(3),a(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    CALL alloc(VM,Z)

    DO I=1,3
       B(I)=0.0_dp
       E(I)=0.0_dp
    ENDDO
    phi=0.0_dp
    SELECT CASE(EL%KIND)
    case(KIND2,kind3,kind5:kind7,KIND16:kind17,KIND20) ! Straight for all practical purposes
       if(present(pos)) then
          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,E,X,pos,k)   ! fringe effect
          ELSE
             call get_Bfield(EL,B,X)   ! fringe effect
          ENDIF
       else
          CALL get_Bfield(EL,B,X)
       endif
    case(kind10)     ! TEAPOT real curvilinear
       if(present(pos)) then
          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,E,X,pos,k)   ! fringe effect
          ELSE
!             if(EL%TP10%electric) then
              call GETELECTRIC(EL%TP10,E,phi,B,VM,X); E(3)=0.d0;
!             else
!              CALL GETMULB_TEAPOT(EL%TP10,B,VM,X)
!             endif
          ENDIF
       else
!             if(EL%TP10%electric) then
              call GETELECTRIC(EL%TP10,E,phi,B,VM,X); E(3)=0.d0;
!             else
!              CALL GETMULB_TEAPOT(EL%TP10,B,VM,X)
!             endif
       endif

    case(KINDPA)     ! fitted field for real magnet
       CALL B_PANCAkE(EL%PA,B,X,POS)
    case(KINDWIGGLER)
       CALL get_z_wi(EL%wi,POS,z)
       CALL B_FIELD(EL%wi,Z,X,B)
    CASE(KIND4)      ! Pill box cavity
      if(EL%C4%n_bessel/=-1) then
        CALL GET_BE_CAV(EL%C4,B,E,X,k)
      else
       call alloc(a,3)
       call alloc(ad,3)
       IF(EL%c4%P%DIR==1) THEN
          Z= pos*el%l/el%p%nst
       ELSE
          Z=EL%L-pos*el%l/el%p%nst
       ENDIF
       call  Abmad_TRANS(EL%C4,Z,X,k,A,AD,B,E)
       call kill(a,3)
       call kill(ad,3)
      endif
    CASE(KIND21)     ! travelling wave cavity
       call alloc(a,3)
       call alloc(ad,3)
        call get_z_cav(EL%cav21,pos,z)

       call A_TRANS(EL%cav21,Z,X,k,A,AD,B,E)
       call kill(a,3)
       call kill(ad,3)

    CASE(kindabell)     ! travelling wave cavity

          IF(POS<0) THEN
             call get_Bfield_fringe(EL,B,E,X,pos,k)   ! fringe effect
          ELSE
              call B_E_FIELD(EL%ab,x,Z,E_in=E,B_in=b)
          ENDIF

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

    CALL KILL(VM,Z)

  end subroutine get_fieldp

  SUBROUTINE get_Bfield_fringeR(EL,B,E,X,pos,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENT),INTENT(IN):: EL
    INTEGER pos
    real(dp) e(3),z,vm,phi
    TYPE(INTERNAL_STATE) K

e=0
    IF(ASSOCIATED(EL%B_SOL)) THEN
       B(1)=  EL%p%DIR*EL%p%CHARGE*(2*Pos+3)*EL%B_SOL*0.5_dp*x(1);    ! POS =-2,-1  (ENT, EXIT)
       B(2)=  EL%p%DIR*EL%p%CHARGE*(2*Pos+3)*EL%B_SOL*0.5_dp*x(3);
       B(3)=0.0_dp;
    else
       b(1)=0.0_dp
       b(2)=0.0_dp
       b(3)=0.0_dp
    ENDIF

    select case(el%kind)
    case(kind2,kind6,kind7,kind10) ! Not exact and Teapot (kind10)
       b(2)=-TAN(EL%p%EDGE(pos+3))*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(1)+b(2)
       if(.not.EL%p%exact) then
          b(1)=-TAN(EL%p%EDGE(pos+3))*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(3)+b(1)
       endif
    case(kind20)  !! likemad=true
       b(2)=-TAN(EL%p%EDGE(pos+3)-EL%p%b0*EL%p%LD/2.0_dp)*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(1)+b(2)
       b(1)=-TAN(EL%p%EDGE(pos+3)-EL%p%b0*EL%p%LD/2.0_dp)*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(3)+b(1)
    case(kind16)  ! likemad=false
    case(kindabell)
            if(pos==-2) then
             z=el%l*(1.0_dp-el%p%dir)/2.0_dp
              else
             z=el%l*(1.0_dp+el%p%dir)/2.0_dp
            endif
             call B_E_FIELD(EL%ab,X,Z,PSIE_in=phi,PSIM_in=vm)
            b(3)= (2*Pos+3)*vm  +b(3)     ! v here  e=-grad phi
            e(3)=-(2*Pos+3)*phi  +e(3)
    end select

    call GET_BZ_fringe(EL,X,B(3),e(3),pos,k)

  END SUBROUTINE get_Bfield_fringeR


  SUBROUTINE get_Bfield_fringeP(EL,B,e,X,pos,k)
    IMPLICIT NONE
    type(REAL_8),INTENT(INOUT):: X(6),B(3), e(3)
    TYPE(ELEMENTP),INTENT(IN):: EL
    INTEGER pos
    TYPE(INTERNAL_STATE) K
    type(REAL_8)vm,phi,z

call alloc(vm,phi,z)
    IF(ASSOCIATED(EL%B_SOL)) THEN
       B(1)= EL%p%DIR*EL%p%CHARGE*(2*Pos+3)*EL%B_SOL*0.5_dp*x(1);    ! POS =-2,-1  (ENT, EXIT)
       B(2)= EL%p%DIR*EL%p%CHARGE*(2*Pos+3)*EL%B_SOL*0.5_dp*x(3);
       B(3)=0.0_dp;
    else
       b(1)=0.0_dp
       b(2)=0.0_dp
       b(3)=0.0_dp
    ENDIF

    select case(el%kind)
    case(kind2,kind6,kind7,kind10) ! Not exact and Teapot (kind10)
       b(2)=-TAN(EL%p%EDGE(pos+3))*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(1)+b(2)
       if(.not.EL%p%exact) then
          b(1)=-TAN(EL%p%EDGE(pos+3))*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(3)+b(1)
       endif
    case(kind20)  !! likemad=true
       b(2)=-TAN(EL%p%EDGE(pos+3)-EL%p%b0*EL%p%LD/2.0_dp)*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(1)+b(2)
       b(1)=-TAN(EL%p%EDGE(pos+3)-EL%p%b0*EL%p%LD/2.0_dp)*EL%p%DIR*EL%p%CHARGE*el%BN(1)*X(3)+b(1)
    case(kind16)  ! likemad=false
    case(kindabell)
            if(pos==-2) then
             z=el%l*(1.0_dp-el%p%dir)/2.0_dp
              else
             z=el%l*(1.0_dp+el%p%dir)/2.0_dp
            endif
             call B_E_FIELD(EL%ab,X,Z,PSIE_in=phi,PSIM_in=vm)
            b(3)=(2*Pos+3)*vm  +b(3)     ! v here  e=-grad phi
            e(3)= -(2*Pos+3)*phi  +e(3)
    end select

    call GET_BZ_fringe(EL,X,B(3),e(3),pos,k)


call kill(vm,phi,z)

  END SUBROUTINE get_Bfield_fringeP

  SUBROUTINE GET_BZ_fringeR(EL,X,bz,ez,pos,k)
    IMPLICIT NONE
    TYPE(ELEMENT),INTENT(IN):: EL
    integer, intent(in) :: pos
    real(dp),INTENT(INOUT):: X(6),Bz,Ez
    real(dp) X1,X3,BBYTW,BBXTW,BBYTWT,E(3),phi,B(3),VM
    INTEGER J,jmax
    real(dp), allocatable :: an(:),bn(:)
    TYPE(INTERNAL_STATE) K

    if(el%electric.and.associated(el%tp10)) then
     call getelectric(EL%tp10,E,phi,B,VM,X)
            bz=(2*Pos+3)*vm       ! probably wrong : here  e=grad phi
            ez=(2*Pos+3)*phi
    else
    bz=0.0_dp
    IF(EL%P%BEND_FRINGE) then
       bz=-(2*Pos+3)*X(3)*EL%BN(1)
    endif

    IF(.not.(k%FRINGE.or.el%p%permfringe/=0)) return

    X1=X(1)
    X3=X(3)
    jmax=MIN(EL%p%NMUL,el%p%HIGHEST_FRINGE)+1

    allocate(an(jmax),bn(jmax))
    an(1)=0.0_dp
    bn(1)=0.0_dp
    do j=2,jmax
       IF(J==2.AND.EL%P%BEND_FRINGE) then
          an(j)=0.0_dp
          bn(j)= EL%AN(j-1)/(j-1)
       else
          an(j)=-EL%BN(j-1)/(j-1)
          bn(j)= EL%AN(j-1)/(j-1)
       endif
    enddo

    BBYTW=BN(jmax)
    BBXTW=AN(jmax)


    DO  J=jmax-1,1,-1
       BBYTWT=X1*BBYTW-X3*BBXTW+BN(J)
       BBXTW=X3*BBYTW+X1*BBXTW+AN(J)
       BBYTW=BBYTWT
    ENDDO

    BZ=-(2*Pos+3)*BBYTW+bz
    deallocate(an,bn)
    endif
  END SUBROUTINE GET_BZ_fringeR

  SUBROUTINE GET_BZ_fringep(EL,X,bz,ez,pos,k)
    IMPLICIT NONE
    TYPE(ELEMENTP),INTENT(IN):: EL
    integer, intent(in) :: pos
    type(real_8),INTENT(INOUT):: X(6),Bz,ez
    type(real_8) X1,X3,BBYTW,BBXTW,BBYTWT,E(3),phi,B(3),VM
    INTEGER J,jmax
    type(real_8), allocatable :: an(:),bn(:)
    TYPE(INTERNAL_STATE) K

    if(el%electric.and.associated(el%tp10)) then
     call alloc(phi,VM)
     call alloc(E);call alloc(b);
     call getelectric(EL%tp10,E,phi,B,VM,X)
            bz=(2*Pos+3)*vm       ! probably wrong : here  e=grad phi
            ez=(2*Pos+3)*phi
     call kill(phi,VM)
     call kill(E);call kill(b);
    else
    bz=0.0_dp
    IF(EL%P%BEND_FRINGE) then
       bz=-(2*Pos+3)*X(3)*EL%BN(1)
    endif

    IF(.not.(k%FRINGE.or.el%p%permfringe/=0)) return
    call alloc(X1,X3,BBYTW,BBXTW,BBYTWT)

    X1=X(1)
    X3=X(3)
    jmax=MIN(EL%p%NMUL,el%p%HIGHEST_FRINGE)+1

    allocate(an(jmax),bn(jmax))
    do j=1,jmax
       call alloc(an(j))
       call alloc(bn(j))
    enddo
    an(1)=0.0_dp
    bn(1)=0.0_dp
    do j=2,jmax
       IF(J==2.AND.EL%P%BEND_FRINGE) then
          an(j)=0.0_dp
          bn(j)= EL%AN(j-1)/(j-1)
       else
          an(j)=-EL%BN(j-1)/(j-1)
          bn(j)= EL%AN(j-1)/(j-1)
       endif
    enddo

    BBYTW=BN(jmax)
    BBXTW=AN(jmax)


    DO  J=jmax-1,1,-1
       BBYTWT=X1*BBYTW-X3*BBXTW+BN(J)
       BBXTW=X3*BBYTW+X1*BBXTW+AN(J)
       BBYTW=BBYTWT
    ENDDO

    BZ=-(2*Pos+3)*BBYTW+bz
    do j=1,jmax
       call kill(an(j))
       call kill(bn(j))
    enddo
    deallocate(an,bn)
    call kill(X1,X3,BBYTW,BBXTW,BBYTWT)
    endif
  END SUBROUTINE GET_BZ_fringep




  SUBROUTINE get_BfieldR(EL,B,X)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENT),INTENT(IN):: EL
    real(dp)  bsol
    INTEGER J



    IF(ASSOCIATED(EL%B_SOL)) THEN
       bsol=EL%B_SOL;
    ELSE
       bsol=0.0_dp
    ENDIF

   ! call GETNEWB(el%an,el%bn,el%b_sol,EL%P%NMUL,B,X)
    call GETNEWB(el%an,el%bn,bsol,EL%P%NMUL,B,X)

  END SUBROUTINE get_BfieldR

  SUBROUTINE get_BfieldP(EL,B,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6),B(3)
    TYPE(ELEMENTP),INTENT(IN):: EL
    TYPE(REAL_8)  bsol
    INTEGER J

    call alloc(bsol)

    IF(ASSOCIATED(EL%B_SOL)) THEN
       bsol=EL%B_SOL;
    ELSE
       bsol=0.0_dp
    ENDIF

!    call GETNEWB(el%an,el%bn,el%b_sol,EL%P%NMUL,B,X)
    call GETNEWB(el%an,el%bn,bsol,EL%P%NMUL,B,X)

       call kill(bsol)

  END SUBROUTINE get_BfieldP

  SUBROUTINE GET_BE_CAVR(EL,B,E,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6),B(3),E(3)
    TYPE(CAV4),INTENT(INOUT):: EL
    real(dp) DF,R2,F,DR2,O,VL
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    real(dp) BBYTWT,BBXTW,BBYTW,x1,x3
    integer J,ko

    E=0.0_dp
    B=0.0_dp
    IF(k%NOCAVITY) RETURN

       if(freq_redefine) then
        O=EL%freq
         else
        O=twopi*EL%freq/CLIGHT
       endif
    VL=EL%volt*volt_c/EL%P%P0C
    do ko=1,el%nf

       DF=0.0_dp
       F=1.0_dp
       R2=1.0_dp

       DO I=1,EL%N_BESSEL
          R2=-R2*(ko*O)**2/4.0_dp/(I+1)**2
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
          BBYTW=0.0_dp
          BBXTW=0.0_dp
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
          BBYTW=0.0_dp
          BBXTW=0.0_dp
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
    integer J,KO

    DO I=1,3
       E(I)=0.0_dp
       B(I)=0.0_dp
    ENDDO

    IF(k%NOCAVITY) RETURN
    CALL ALLOC(DF,R2,F,DR2,O,VL)
    CALL ALLOC(BBYTWT,BBXTW,BBYTW,x1,x3)

       if(freq_redefine) then
        O=EL%freq
         else
        O=twopi*EL%freq/CLIGHT
       endif
    VL=EL%volt*volt_c/EL%P%P0C
    do ko=1,el%nf

       DF=0.0_dp
       F=1.0_dp
       R2=1.0_dp

       DO I=1,EL%N_BESSEL
          R2=-R2*(ko*O)**2/4.0_dp/(I+1)**2
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
          BBYTW=0.0_dp
          BBXTW=0.0_dp
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
          BBYTW=0.0_dp
          BBXTW=0.0_dp
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
    real(dp) BE(nbe)
    Be(1)=X(1);
    Be(2)=X(3);
    Be(3)=0.0_dp;

    CALL trackg(EL%B(POS),BE)

    b(1)=EL%SCALE*be(1)
    b(2)=EL%SCALE*be(2)
    !    b(3)=EL%SCALE*el%p%charge*el%p%dir*b(3)
    b(3)=EL%SCALE*be(3)

  END subroutine B_PANCAkEr

  subroutine B_PANCAkEp(EL,B,X,POS)
    IMPLICIT NONE
    type(real_8), INTENT(INout) :: X(6)
    INTEGER, INTENT(IN) :: POS
    TYPE(PANCAKEP),  INTENT(INOUT) :: EL
    type(real_8), INTENT(INOUT) ::B(3)
    type(real_8) be(nbe)

    call alloc(be)
    Be(1)=X(1);
    Be(2)=X(3);
    Be(3)=0.0_dp;


    CALL trackg(EL%B(POS),Be)

    b(1)=EL%SCALE*be(1)
    b(2)=EL%SCALE*be(2)
    !    b(3)=EL%SCALE*el%p%charge*el%p%dir*b(3)
    b(3)=EL%SCALE*be(3)
    call kill(be)
  END subroutine B_PANCAkEp

  subroutine B_PARA_PERPr(k,EL,TEAPOT_LIKE,X,B,BPA,BPE,XP,XPA,e,EF,EFB,EFD,POS)
    IMPLICIT NONE
    REAL(DP),  INTENT(INout) :: X(6)
    TYPE(ELEMENT),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    REAL(DP),  INTENT(INout) :: B(3),BPA(3),BPE(3),XP(2),XPA(2),e(3)
    REAL(DP),  OPTIONAL ::EF(3),EFB(3),EFD(3)
    integer, optional,intent(in) :: pos
    INTEGER TEAPOT_LIKE,i
    REAL(DP) be
    type(internal_state) k
    P=>EL%P

    !  this routines gives us  B parallel and B perpendicular
    ! Also if EF is present, E perpendicular times beta is return

    call DIRECTION_V(k,EL,X,E,XP,XPA,POS)

    be=b(1)*e(1)+b(2)*e(2)+b(3)*e(3)

    do i=1,3
       BPA(i)=be*e(i)
    enddo
    do i=1,3
       BPE(i)=B(i)-BPA(i)
    enddo

    IF(PRESENT(EF)) THEN

       EFB(1)=-EF(2)*E(3)+EF(3)*E(2)      ! changed sign txE of Barber
       EFB(2)=-EF(3)*E(1)+EF(1)*E(3)
       EFB(3)=-EF(1)*E(2)+EF(2)*E(1)
       be=EF(1)*e(1)+EF(2)*e(2)+EF(3)*e(3)
       do i=1,3
         EFD(i)=be*e(i)
        enddo

    endif
  END subroutine B_PARA_PERPr

  subroutine B_PARA_PERPp(k,EL,TEAPOT_LIKE,X,B,BPA,BPE,XP,XPA,e,EF,EFB,EFD,pos)
    IMPLICIT NONE
    type(real_8),  INTENT(INout) :: X(6)
    TYPE(ELEMENTP),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    type(real_8),  INTENT(INout) :: B(3),BPA(3),BPE(3),XP(2),XPA(2),e(3)
    type(real_8),  OPTIONAL ::EF(3),EFB(3),EFD(3)
    INTEGER TEAPOT_LIKE,i
    type(real_8) be
    type(internal_state) k
    integer, optional,intent(in) :: pos

    P=>EL%P

     call alloc(be);

    call DIRECTION_V(k,EL,X,E,XP,XPA,POS)

    be=b(1)*e(1)+b(2)*e(2)+b(3)*e(3)

    do i=1,3
       BPA(i)=be*e(i)
    enddo
    do i=1,3
       BPE(i)=B(i)-BPA(i)
    enddo

    IF(PRESENT(EF)) THEN
       EFB(1)=-EF(2)*E(3)+EF(3)*E(2)      ! changed sign txE of Barber
       EFB(2)=-EF(3)*E(1)+EF(1)*E(3)
       EFB(3)=-EF(1)*E(2)+EF(2)*E(1)
       be=EF(1)*e(1)+EF(2)*e(2)+EF(3)*e(3)
       do i=1,3
         EFD(i)=be*e(i)
        enddo
    ENDIF

     call kill(be);

  END subroutine B_PARA_PERPp


  subroutine DIRECTION_VR(k,EL,X,E,XP,XPA,POS)
    IMPLICIT NONE
    REAL(DP),  INTENT(INout) :: X(6),XP(2),XPA(2)
    TYPE(ELEMENT),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    REAL(DP),  INTENT(INOUT) ::E(3)
    REAL(DP) N,H,DP1,A,AP,B,BP,z,ve,AV(3)
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
       DP1=root(1.0_dp+2.0_dp*X(5)/P%BETA0+X(5)**2)
    ELSE
       DP1=1.0_dp+X(5)
    ENDIF

    IF(EL%KIND/=KINDPA) THEN

       IF(ASSOCIATED(EL%B_SOL)) THEN  !SOLENOID

          XPA(1)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp)
          XPA(2)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp)
          N=root(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N
       ELSEif(el%kind==kindwiggler) then
          CALL get_z_wi(EL%wi,POS,z)

          if(el%wi%xprime) then
           Xpa(1)=X(2)
           Xpa(2)=X(4)
          else
          CALL COMPX(EL%wi,Z,X,A,AP)
          CALL COMPY(EL%wi,Z,X,B,BP)
           Xpa(1)=X(2)-A
           Xpa(2)=X(4)-B
          endif
          N=root(DP1**2-Xpa(1)**2-Xpa(2)**2)

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
          N=root(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N
       ELSEif(el%kind==kindabell) then



          if(el%ab%xprime) then
               H=1.0_dp+el%ab%hc*X(1)
               N=root(H**2+X(2)**2+X(4)**2)
               E(1)=X(2)/N
               E(2)=X(4)/N
               E(3)=H/N
               XPA(1)=X(2)
               XPA(2)=X(4)
               XP(1)=X(2)
               XP(2)=X(4)
          else
              CALL get_z_ab(EL%ab,POS,z)
              call B_E_FIELD(EL%ab,X,Z,psie_in=ve,A_in=av)

            IF(k%TIME) THEN
               DP1=root(1.0_dp+2.0_dp*(X(5)+EL%P%CHARGE*ve)/P%BETA0+(X(5)+ve)**2)
            ELSE
               DP1=1.0_dp+(X(5)+EL%P%CHARGE*ve)
            ENDIF

             Xpa(1)=X(2)-EL%P%CHARGE*AV(1)
             Xpa(2)=X(4)-EL%P%CHARGE*AV(2)
             N=root(DP1**2-Xpa(1)**2-Xpa(2)**2)

             E(1)=Xpa(1)/DP1
             E(2)=Xpa(2)/DP1
             E(3)=N/DP1
             XP(1)=XPA(1)/N
             XP(2)=XPA(2)/N

             endif
       else


          N=root(DP1**2-X(2)**2-X(4)**2)

          E(1)=X(2)/DP1
          E(2)=X(4)/DP1
          E(3)=N/DP1
          XPA(1)=X(2)
          XPA(2)=X(4)
          XP(1)=X(2)/N
          XP(2)=X(4)/N
       ENDIF

    ELSE    ! NON CANONICAL VARIABLES
       H=1.0_dp+el%pa%hc*X(1)
       N=root(H**2+X(2)**2+X(4)**2)
       E(1)=X(2)/N
       E(2)=X(4)/N
       E(3)=H/N
       XPA(1)=X(2)
       XPA(2)=X(4)
       XP(1)=X(2)
       XP(2)=X(4)

    ENDIF

!    E(1)=EL%P%dir*E(1)
!    E(2)=EL%P%dir*E(2)    etienne 2016_5_9
    E(3)=EL%P%dir*E(3)



  END subroutine DIRECTION_VR

  subroutine DIRECTION_VP(k,EL,X,E,XP,XPA,POS)
    IMPLICIT NONE
    type(real_8), INTENT(INout) :: X(6),XP(2),XPA(2)
    TYPE(ELEMENTP),  pointer :: EL
    TYPE(MAGNET_CHART),  pointer :: P
    type(real_8), INTENT(INOUT) ::E(3)
    type(real_8) N,H,DP1,A,AP,B,BP,z,AV(3),ve
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer, optional,intent(in) :: pos

    P=>EL%P

    CALL ALLOC(N,H,DP1,A,AP,B,BP,z,ve )
    CALL ALLOC(AV )

    IF(k%TIME) THEN
       DP1=SQRT(1.0_dp+2.0_dp*X(5)/P%BETA0+X(5)**2)
    ELSE
       DP1=1.0_dp+X(5)
    ENDIF

    IF(EL%KIND/=KINDPA) THEN

       IF(ASSOCIATED(EL%B_SOL)) THEN  !SOLENOID

          XPA(1)=(X(2)+EL%B_SOL*EL%P%CHARGE*X(3)/2.0_dp)
          XPA(2)=(X(4)-EL%B_SOL*EL%P%CHARGE*X(1)/2.0_dp)
          N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

          E(1)=Xpa(1)/DP1
          E(2)=Xpa(2)/DP1
          E(3)=N/DP1
          XP(1)=XPA(1)/N
          XP(2)=XPA(2)/N
       ELSEif(el%kind==kindwiggler) then
          CALL get_z_wi(EL%wi,POS,z)
          if(el%wi%xprime) then
           Xpa(1)=X(2)
           Xpa(2)=X(4)
          else
          CALL COMPX(EL%wi,Z,X,A,AP)
          CALL COMPY(EL%wi,Z,X,B,BP)
           Xpa(1)=X(2)-A
           Xpa(2)=X(4)-B
          endif
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
       ELSEif(el%kind==kindabell) then



          if(el%ab%xprime) then
               H=1.0_dp+el%ab%hc*X(1)
               N=SQRT(H**2+X(2)**2+X(4)**2)
               E(1)=X(2)/N
               E(2)=X(4)/N
               E(3)=H/N
               XPA(1)=X(2)
               XPA(2)=X(4)
               XP(1)=X(2)
               XP(2)=X(4)
          else
              CALL get_z_ab(EL%ab,POS,z)
              call B_E_FIELD(EL%ab,X,Z,psie_in=ve,A_in=av)

            IF(k%TIME) THEN
               DP1=SQRT(1.0_dp+2.0_dp*(X(5)+EL%P%CHARGE*ve)/P%BETA0+(X(5)+ve)**2)
            ELSE
               DP1=1.0_dp+(X(5)+EL%P%CHARGE*ve)
            ENDIF

             Xpa(1)=X(2)-EL%P%CHARGE*AV(1)
             Xpa(2)=X(4)-EL%P%CHARGE*AV(2)
             N=SQRT(DP1**2-Xpa(1)**2-Xpa(2)**2)

             E(1)=Xpa(1)/DP1
             E(2)=Xpa(2)/DP1
             E(3)=N/DP1
             XP(1)=XPA(1)/N
             XP(2)=XPA(2)/N

             endif
       else


          N=sqrt(DP1**2-X(2)**2-X(4)**2)

          E(1)=X(2)/DP1
          E(2)=X(4)/DP1
          E(3)=N/DP1
          XPA(1)=X(2)
          XPA(2)=X(4)
          XP(1)=X(2)/N
          XP(2)=X(4)/N
       ENDIF

    ELSE    ! NON CANONICAL VARIABLES
       H=1.0_dp+el%pa%hc*X(1)
       N=SQRT(H**2+X(2)**2+X(4)**2)
       E(1)=X(2)/N
       E(2)=X(4)/N
       E(3)=H/N
       XPA(1)=X(2)
       XPA(2)=X(4)
       XP(1)=X(2)
       XP(2)=X(4)

    ENDIF


!    E(1)=EL%P%dir*E(1)
!    E(2)=EL%P%dir*E(2)    etienne 2016_5_9
    E(3)=EL%P%dir*E(3)


     CALL kill(N,H,DP1,A,AP,B,BP,z,ve )
     CALL kill(AV )

  END subroutine DIRECTION_VP

!!!!!!!!!!!!   GLOBAL TRACKING ROUTINES    !!!!!!!!!!!!

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_R(xs,k,fibre1,fibre2,node1,node2) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    type(probe), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE),optional, POINTER :: node1,node2
    TYPE (fibre),optional, POINTER :: fibre1,fibre2
    TYPE (INTEGRATION_NODE), POINTER :: C,n1,n2,last
    logical donew
    real(dp) beta
    !    INTEGER,TARGET :: CHARGE

    !    if(present(node1))CHARGE=NODE1%PARENT_FIBRE%CHARGE
    !    if(present(fibre1))CHARGE=fibre1%CHARGE

    !    CALL RESET_APERTURE_FLAG
    nullify(n1)
    nullify(n2)


    xs%u=my_false

    if(present(node1)) n1=>node1
    if(present(node2)) n2=>node2
    if(present(fibre1)) n1=>fibre1%t1
    if(present(fibre2)) n2=>fibre2%t1

    c=>n1

    if(associated(n2)) then
       nullify(last)
    else
       if(n1%parent_fibre%parent_layout%closed) then
          last=>n1%previous
          n2  =>last
       else
          last=>n1%parent_fibre%parent_layout%t%end
          n2  =>n1%parent_fibre%parent_layout%t%end
       endif
    endif

    donew=(.not.(full_way.or.k%full_way)).and.(.not.present(node1)).and.(.not.present(node2))

    if(donew) then   ! actually calling old stuff pre-node
     call TRACK(xs%x,K,fibre1,fibre2=fibre2)
    else
     if(use_bmad_units.and.(.not.inside_bmad)) then 
       beta=C%PARENT_FIBRE%beta0
       if(C%PARENT_FIBRE%PATCH%ENERGY==4) beta=C%PARENT_FIBRE%PATCH%b0b
       call convert_bmad_to_ptc(xs,beta,k%time)
     endif
     DO  WHILE(.not.ASSOCIATED(C,n2))
        CALL TRACK_NODE_PROBE(C,XS,K)
        if(.not.check_stable) exit

        C=>C%NEXT
     ENDDO
     if(associated(last).and.check_stable) then
       CALL TRACK_NODE_PROBE(last,XS,K)
     endif
    if(use_bmad_units.and.(.not.inside_bmad)) then 
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==5) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_ptc_to_bmad(xs,beta,k%time)
    endif
    endif


    C_%STABLE_DA=.true.


    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_R

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_P(xs,k,fibre1,fibre2,node1,node2) ! Tracks double from i1 to i2 in state k
    use s_extend_poly, only : elem_name ! LD: 22.03.2019
    IMPLICIT NONE
    type(probe_8), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE),optional, POINTER :: node1,node2
    TYPE (fibre),optional, POINTER :: fibre1,fibre2
    TYPE (INTEGRATION_NODE), POINTER :: C,n1,n2,last
    logical donew
    real(dp) beta
    !    INTEGER,TARGET :: CHARGE

    !    if(present(node1))CHARGE=NODE1%PARENT_FIBRE%CHARGE
    !    if(present(fibre1))CHARGE=fibre1%CHARGE

    !    CALL RESET_APERTURE_FLAG
    nullify(n1)
    nullify(n2)


    xs%u=my_false

    if(present(node1)) n1=>node1
    if(present(node2)) n2=>node2
    if(present(fibre1)) n1=>fibre1%t1
    if(present(fibre2)) n2=>fibre2%t1

    c=>n1

    if(associated(n2)) then
       nullify(last)
    else
       if(n1%parent_fibre%parent_layout%closed) then
          last=>n1%previous
          n2  =>last
       else
          last=>n1%parent_fibre%parent_layout%t%end
          n2  =>n1%parent_fibre%parent_layout%t%end
       endif
    endif


 !   DO  WHILE(.not.ASSOCIATED(C,n2))

  !     CALL TRACK_NODE_PROBE(C,XS,K)
  !     if(.not.check_stable) exit!

!       C=>C%NEXT
!    ENDDO

!    if(associated(last).and.check_stable) then
!       CALL TRACK_NODE_PROBE(last,XS,K)
!    endif

    donew=(.not.(full_way.or.k%full_way)).and.(.not.present(node1)).and.(.not.present(node2))

    if(donew) then   ! actually calling old stuff pre-node
     call TRACK(xs%x,K,fibre1,fibre2=fibre2)
    else
     if(use_bmad_units.and.(.not.inside_bmad)) then 
       beta=C%PARENT_FIBRE%beta0
       if(C%PARENT_FIBRE%PATCH%ENERGY==4) beta=C%PARENT_FIBRE%PATCH%b0b
       call convert_bmad_to_ptc(xs,beta,k%time)
     endif
     DO  WHILE(.not.ASSOCIATED(C,n2))
        CALL TRACK_NODE_PROBE(C,XS,K)
        if(.not.check_stable) exit

        C=>C%NEXT
     ENDDO
     if(associated(last).and.check_stable) then
       elem_name = C%PARENT_FIBRE%MAGP%name  ! LD: 22.03.2019
       CALL TRACK_NODE_PROBE(last,XS,K)
     endif
    if(use_bmad_units.and.(.not.inside_bmad)) then 
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==5) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_ptc_to_bmad(xs,beta,k%time)
    endif
    endif


    C_%STABLE_DA=.true.

    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_t12_P



  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_R(R,xs,k,I1,I2) ! Tracks probe from integer node i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C
    real(dp) beta
    ! CALL RESET_APERTURE_FLAG
    xs%u=my_false

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )


    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    if(use_bmad_units.and.(.not.inside_bmad)) then 
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==4) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_bmad_to_ptc(xs,beta,k%time)
    endif

    DO  WHILE(J<I22.AND.ASSOCIATED(C))
       CALL TRACK_NODE_PROBE(C,XS,K)

       if(.not.check_stable) exit


       C=>C%NEXT
       J=J+1
    ENDDO

    if(use_bmad_units.and.(.not.inside_bmad)) then 
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==5) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_ptc_to_bmad(xs,beta,k%time)
    endif

    C_%STABLE_DA=.true.

    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_R


  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_P(R,XS,k,I1,I2) ! Tracks double from i1 to i2 in state k
    use s_extend_poly, only : elem_name ! LD: 22.03.2019
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(probe_8), INTENT(INOUT):: XS
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J   ,i22
    TYPE (INTEGRATION_NODE), POINTER :: C
    real(dp) beta

    !    CALL RESET_APERTURE_FLAG

    xs%u=my_false

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif


    J=I1

    if(use_bmad_units.and.(.not.inside_bmad)) then 
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==4) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_bmad_to_ptc(xs,beta,k%time)
    endif

    DO  WHILE(J<I22.AND.ASSOCIATED(C))
        elem_name = C%PARENT_FIBRE%MAGP%name  ! LD: 22.03.2019
        CALL TRACK_NODE_PROBE(C,XS,K)  !,R%charge)
        if(.not.check_stable) exit

       C=>C%NEXT
       J=J+1
    ENDDO

    if(use_bmad_units.and.(.not.inside_bmad)) then 
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==5) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_ptc_to_bmad(xs,beta,k%time)
    endif

    C_%STABLE_DA=.true.

    !    if(c_%watch_user) ALLOW_TRACKING=.FALSE.


  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_pr_s12_P


  SUBROUTINE TRACK_LAYOUT_FLAG_probe_spin12r(r,xS,k,fibre1,fibre2,node1,node2) ! integer fibre i1 to i2
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

  SUBROUTINE TRACK_LAYOUT_FLAG_spin12r_x(r,x,k, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe) xs
    real(dp),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    integer i1,i2
    logical donew
    !    logical(lp), optional ::u
    !    type(integration_node),optional, pointer :: t
      donew=(.not.(full_way.or.k%full_way)).and.(.not.present(node1)).and.(.not.present(node2))

    if(donew) then
      i1=fibre1
      if(present(fibre2) )THEN
         i2=FIBRE2
      else
         I2=r%n+i1
      endif
      if(i2<i1) then
       i2=r%n+i2
      endif
     CALL TRACK(r,x,I1,I2,K)
     else
       if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)

        xs%u=my_false
        XS=X
         CALL TRACK_PROBE(r,xs,K, fibre1,fibre2,node1,node2)
       X=XS%X
    endif
    !    if(present(u)) u=xs%u

    !    if(present(t)) THEN
    !       t=>xs%lost_node
    !       NULLIFY(xs%lost_node)
    !    ENDIF

  END SUBROUTINE TRACK_LAYOUT_FLAG_spin12r_x





  SUBROUTINE TRACK_LAYOUT_FLAG_spin12p_x(r,x,k, fibre1,fibre2,node1,node2)  ! fibre i1 to i2
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(probe_8) xs
    type(real_8),target,intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    integer,optional:: fibre1,fibre2,node1,node2
    integer i1,i2
    logical donew

      donew=(.not.(full_way.or.k%full_way)).and.(.not.present(node1)).and.(.not.present(node2))

    if(donew) then
      i1=fibre1
      if(present(fibre2) )THEN
         i2=FIBRE2
      else
         I2=r%n+i1
      endif
      if(i2<i1) then
       i2=r%n+i2
      endif
     CALL TRACK(r,x,I1,I2,K)
     else
       call alloc(xs)
       if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)

        xs%u=my_false
        XS%X=X
         CALL TRACK_PROBE(r,xs,K, fibre1,fibre2,node1,node2)
        X=XS%X
        call kill(xs)
    endif



  END SUBROUTINE TRACK_LAYOUT_FLAG_spin12p_x

  SUBROUTINE TRACK_LAYOUT_FLAG_spint12r_x(x,k, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    type(probe) xs
    real(dp),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    type(fibre),optional,pointer:: fibre1,fibre2
    type(integration_node),optional,pointer:: node1,node2
    !   logical(lp), optional ::u
    !   type(integration_node),optional, pointer :: t


    xs%u=my_false
    XS=X
    !   if(present(t)) THEN
    !       ALLOCATE(xs%lost_node)
    !       t=>xs%lost_node
    !     nullify(t)
    !  ENDIF
    !    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE(xs,K, fibre1,fibre2,node1,node2)
    !   if(present(u)) u=xs%u
    X=XS%X
    !   if(present(t)) THEN
    !      t=>xs%lost_node
    !       deallocate(xs%lost_node)
    !      NULLIFY(xs%lost_node)
    !   ENDIF

  END SUBROUTINE TRACK_LAYOUT_FLAG_spint12r_x

  SUBROUTINE TRACK_LAYOUT_FLAG_spint12p_x(x,k, fibre1,fibre2,node1,node2) ! fibre i1 to i2
    IMPLICIT NONE
    type(probe_8) xs
    type(real_8),intent(INOUT) ::  x(6)
    TYPE(INTERNAL_STATE) K
    type(fibre),optional,pointer:: fibre1,fibre2
    type(integration_node),optional,pointer:: node1,node2
    !    logical(lp), optional ::u
    !    type(integration_node),optional, pointer :: t


    call alloc(xs)
    xs%u=my_false
    XS%X=X
    !  if(present(t)) THEN
    !       ALLOCATE(xs%lost_node)
    !       t=>xs%lost_node
    !    nullify(t)
    ! ENDIF
    !    IF(I22==I11.AND.I2>I1) I22=I11+R%T%N

    CALL TRACK_PROBE(xs,K, fibre1,fibre2,node1,node2)
    !  if(present(u)) u=xs%u
    X=XS%X
    !  if(present(t)) THEN
    !     t=>xs%lost_node
    !      deallocate(xs%lost_node)
    !     NULLIFY(xs%lost_node)
    !  ENDIF
    call kill(xs)

  END SUBROUTINE TRACK_LAYOUT_FLAG_spint12p_x



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
    if(.not.check_stable) return

    ref0=my_false
    if(present(ref)) ref0=ref

    IF(.NOT.ASSOCIATED(T%B)) THEN
       call FILL_SURVEY_DATA_IN_NODE_LAYOUT(t%parent_fibre%parent_LAYOUT)
       WRITE(6,*)  " SURVEY DONE FOR THIN LAYOUT IN TRACK_NODE_LAYOUT_FLAG_spin_v "
    ENDIF

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
       reference_ray=0.0_dp
       reference_ray(1)=t%ref(1)
       reference_ray(3)=t%ref(2)
    endif

    CALL TRACK_NODE_PROBE(T,xs,K)  !,t%parent_fibre%CHARGE)

    if(.not.ref0.and.check_stable) CALL TRACK_NODE_PROBE(T,XS_REF,K)  !,t%parent_fibre%CHARGE)

    v%u(1)=XS%u

    if(.not.ref0) then
       v%u(2)=XS_REF%u
       v%reference_ray=XS_REF%x
    else
       v%u(2)=my_false
       v%reference_ray=0.0_dp
       v%reference_ray(1)=t%ref(3)
       v%reference_ray(3)=t%ref(4)
    endif

    v%x=XS%x

    IF(V%U(1).OR.V%U(2)) RETURN


    SC=1.0_dp
    IF(v%SCALE/=0.0_dp) SC=v%SCALE
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

    !    if(.not.check_stable) return
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

    !    if(.not.check_stable) return
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
    REAL(DP) FAC,DS,beta
    logical useptc,dofix0,dofix,doonemap
    type(tree_element), pointer :: arbre(:)
!    logical(lp) bmad
    IF(.NOT.CHECK_STABLE) then
       CALL RESET_APERTURE_FLAG
    endif

    !    if(xs%u) return
    C%PARENT_FIBRE%MAG%P%DIR    => C%PARENT_FIBRE%DIR
    C%PARENT_FIBRE%MAG%P%beta0  => C%PARENT_FIBRE%beta0
    C%PARENT_FIBRE%MAG%P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    C%PARENT_FIBRE%MAG%P%GAMBET => C%PARENT_FIBRE%GAMBET
    C%PARENT_FIBRE%MAG%P%MASS => C%PARENT_FIBRE%MASS
    C%PARENT_FIBRE%MAG%P%ag => C%PARENT_FIBRE%ag
    C%PARENT_FIBRE%MAG%P%CHARGE=>C%PARENT_FIBRE%CHARGE



    if(full_way.or.k%full_way) then
     useptc=.true.

     
  !  if(.not.(k%nocavity.and.(C%PARENT_FIBRE%MAG%kind==kind4.or.C%PARENT_FIBRE%MAG%kind==kind21))) then
     if(C%PARENT_FIBRE%dir==1) then
       if(C%PARENT_FIBRE%MAG%skip_ptc_f==1) return
       if(associated(C%PARENT_FIBRE%MAG%forward)) then
         if(C%PARENT_FIBRE%MAG%usef) useptc=.false.
          arbre=>C%PARENT_FIBRE%MAG%forward
          doonemap=C%PARENT_FIBRE%MAG%do1mapf
       endif
     else
       if(C%PARENT_FIBRE%MAG%skip_ptc_b==1) return
       if(associated(C%PARENT_FIBRE%MAG%backward)) then
         if(C%PARENT_FIBRE%MAG%useb) useptc=.false.
          arbre=>C%PARENT_FIBRE%MAG%backward
          doonemap=C%PARENT_FIBRE%MAG%do1mapb
       endif
     endif
 !   endif ! cavity
 

    if(use_bmad_units.and.inside_bmad) then
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==4) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_bmad_to_ptc(xs,beta,k%time)
    endif

    IF(K%MODULATION.and.xs%nac/=0) THEN !modulate
       if(c%parent_fibre%mag%slow_ac/=0) CALL MODULATE(C,XS,K) !modulate
       CALL TRACK_MODULATION(C,XS,K) !modulate
    ENDIF !modulate



    if(c%cas==0) then
       if(useptc) then
       ds=c%parent_fibre%MAG%L/c%parent_fibre%MAG%p%nst
       fac=0.5_dp
        call PUSH_SPIN(c,ds,FAC,XS,my_true,k,C%POS_IN_FIBRE-2)   ! -3 before....
        CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
        call PUSH_SPIN(c,ds,FAC,XS,my_false,k,C%POS_IN_FIBRE-2)
       elseif(doonemap) then
 
          if(C%POS_IN_FIBRE-2==1) then 
             dofix0=.true.;dofix=.true.
             call track_TREE_probe_complex(arbre,xs,dofix0,dofix,k)
          endif
       else
          dofix0=.false.;dofix=.false.
          if(C%POS_IN_FIBRE-2==1) dofix0=.true.
          if(C%POS_IN_FIBRE-C%PARENT_FIBRE%MAG%p%nst==2) dofix=.true.
        call track_TREE_probe_complex(arbre,xs,dofix0,dofix,k)
       endif
    elseIF(c%cas==case1.and.useptc) then
       CALL TRACK_FRINGE_spin(C,XS,K)
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    elseIF(c%cas==case2.and.useptc) then
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
       CALL TRACK_FRINGE_spin(C,XS,K)
    else
       IF(c%cas==caseP1) THEN

          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
          if(k%spin) then
            CALL TRACK_SPIN_FRONT(C%PARENT_FIBRE,XS)
            if(xs%use_q) xs%q%x=xs%q%x/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)

          endif
       
       ELSEif(c%cas==caseP2) THEN
          if(k%spin) then

                 CALL TRACK_SPIN_BACK(C%PARENT_FIBRE,XS)
             if(xs%use_q) xs%q%x=xs%q%x/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)

           endif
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
     ENDIF

    endif
    IF(K%MODULATION.and.xs%nac/=0.and.c%parent_fibre%mag%slow_ac/=0) then!modulate
       CALL restore_ANBN_SINGLE(C%PARENT_FIBRE%MAG,C%PARENT_FIBRE%MAGP)
    ENDIF  !modulate
  !  IF((K%MODULATION.or.ramp).and.c%parent_fibre%mag%slow_ac) THEN  !modulate
  !     CALL restore_ANBN_SINGLE(C%PARENT_FIBRE%MAG,C%PARENT_FIBRE%MAGP)
  !  ENDIF  !modulate
    if(use_bmad_units.and.inside_bmad) then
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==5) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_ptc_to_bmad(xs,beta,k%time)
    endif
 else ! full_way



    if(c%cas==0) then



        CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE


    elseIF(c%cas==case1) then
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    elseIF(c%cas==case2) then
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    else
       IF(c%cas==caseP1) THEN
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
       ELSEif(c%cas==caseP2) THEN
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE

     ENDIF

    endif

endif ! full_way
    xs%u=.not.check_stable
    if(xs%u) then
       lost_fibre=>c%parent_fibre
       lost_node=>c
       xlost=xs%x
    endif
       xs%last_node=>c
       xs%e=global_e
  END SUBROUTINE TRACK_NODE_FLAG_probe_R

  SUBROUTINE TRACK_NODE_FLAG_probe_P(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe_8), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    REAL(DP) FAC,beta
    type(real_8) ds
    logical(lp) CHECK_KNOB
    integer(2), pointer,dimension(:)::AN,BN
    integer ki
    logical useptc,dofix0,dofix,doonemap
    type(tree_element), pointer :: arbre(:)
!    logical(lp) bmad
    !   if(xs%u) return

    IF(.NOT.CHECK_STABLE) then
       CALL RESET_APERTURE_FLAG
    endif
    ki=c%parent_fibre%MAGp%kind
    C%PARENT_FIBRE%MAGp%P%DIR    => C%PARENT_FIBRE%DIR
    C%PARENT_FIBRE%MAGp%P%beta0  => C%PARENT_FIBRE%beta0
    C%PARENT_FIBRE%MAGp%P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    C%PARENT_FIBRE%MAGp%P%GAMBET => C%PARENT_FIBRE%GAMBET
    C%PARENT_FIBRE%MAGP%P%MASS => C%PARENT_FIBRE%MASS
    C%PARENT_FIBRE%MAGP%P%ag => C%PARENT_FIBRE%ag
    C%PARENT_FIBRE%MAGp%P%CHARGE=>C%PARENT_FIBRE%CHARGE


    if(full_way.or.k%full_way) then
    useptc=.true.
!    if(.not.(k%nocavity.and.(ki==kind4.or.ki==kind21))) then
     if(C%PARENT_FIBRE%dir==1) then
       if(C%PARENT_FIBRE%MAGp%skip_ptc_f==1) return
       if(associated(C%PARENT_FIBRE%MAGP%forward)) then
         if(C%PARENT_FIBRE%MAGP%usef) useptc=.false.
          arbre=>C%PARENT_FIBRE%MAGP%forward
          doonemap=C%PARENT_FIBRE%MAGp%do1mapf
       endif
     else
       if(C%PARENT_FIBRE%MAGp%skip_ptc_b==1) return
       if(associated(C%PARENT_FIBRE%MAGP%backward)) then
         if(C%PARENT_FIBRE%MAGP%useb) useptc=.false.
          arbre=>C%PARENT_FIBRE%MAGP%backward
          doonemap=C%PARENT_FIBRE%MAGp%do1mapb
       endif
     endif
!    endif
 
 

    if(use_bmad_units.and.inside_bmad) then
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==4) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_bmad_to_ptc(xs,beta,k%time)
    endif

    IF(K%MODULATION.and.xs%nac/=0) then
       if(c%parent_fibre%mag%slow_ac/=0)  CALL MODULATE(C,XS,K) !modulate
       CALL TRACK_MODULATION(C,XS,K) !modulate
    ENDIF !modulate




    CALL ALLOC(DS)

    !      if(associated(c%bb)) call BBKICK(c%BB,XS%X)

    if(c%cas==0) then
       if(useptc) then
       ds=c%parent_fibre%MAGp%L/c%parent_fibre%MAG%p%nst
       fac=0.5_dp
        if(ki==kind10)CALL MAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,K)
        call PUSH_SPIN(c,ds,FAC,XS,my_true,k,C%POS_IN_FIBRE-2)    ! -3 before
         if(ki==kind10)CALL UNMAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)
        CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
         if(ki==kind10)CALL MAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)
        call PUSH_SPIN(c,ds,FAC,XS,my_false,k,C%POS_IN_FIBRE-2)
         if(ki==kind10)CALL UNMAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)
       elseif(doonemap) then

          if(C%POS_IN_FIBRE-2==1) then
                     dofix0=.true.;dofix=.true.
 
           call track_TREE_probe_complex(arbre,xs,dofix0,dofix,k)  
 
          endif
       else
          dofix0=.false.;dofix=.false.
          if(C%POS_IN_FIBRE-2==1) dofix0=.true.
          if(C%POS_IN_FIBRE-C%PARENT_FIBRE%MAGp%p%nst==2) dofix=.true.
        call track_TREE_probe_complex(arbre,xs,dofix0,dofix,k)
       endif
    elseIF(c%cas==case1.and.useptc) then
if(ki==kind10)CALL MAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)
       CALL TRACK_FRINGE_spin(C,XS,K)
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
if(ki==kind10)CALL UNMAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)
    elseIF(c%cas==case2.and.useptc) then
if(ki==kind10)CALL MAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)
         CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
         CALL TRACK_FRINGE_spin(C,XS,K)
       !        CALL  (C,XS,K)
if(ki==kind10)CALL UNMAKEPOTKNOB(c%parent_fibre%MAGp%TP10,CHECK_KNOB,AN,BN,k)

    else
       IF(c%cas==caseP1) THEN
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
          if(k%spin) then

                 CALL TRACK_SPIN_FRONT(C%PARENT_FIBRE,XS)
       if(xs%use_q) then
           ds=1.0_dp/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)
           xs%q%x(0)=xs%q%x(0)*ds
           xs%q%x(1)=xs%q%x(1)*ds
           xs%q%x(2)=xs%q%x(2)*ds
           xs%q%x(3)=xs%q%x(3)*ds
        endif

          endif
       ELSEif(c%cas==caseP2) THEN
          if(k%spin) then

                 CALL TRACK_SPIN_BACK(C%PARENT_FIBRE,XS)
        if(xs%use_q) then
           ds=1.0_dp/sqrt(xs%q%x(1)**2+xs%q%x(2)**2+xs%q%x(3)**2+xs%q%x(0)**2)
           xs%q%x(0)=xs%q%x(0)*ds
           xs%q%x(1)=xs%q%x(1)*ds
           xs%q%x(2)=xs%q%x(2)*ds
           xs%q%x(3)=xs%q%x(3)*ds
        endif
           endif
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
     ENDIF


    endif


    IF(K%MODULATION.and.xs%nac/=0.and.c%parent_fibre%mag%slow_ac/=0) then
       CALL restore_ANBN_SINGLE(C%PARENT_FIBRE%MAG,C%PARENT_FIBRE%MAGP)
    ENDIF  !modulate


    call kill(ds)


    if(use_bmad_units.and.inside_bmad) then
      beta=C%PARENT_FIBRE%beta0
      if(C%PARENT_FIBRE%PATCH%ENERGY==5) beta=C%PARENT_FIBRE%PATCH%b0b
      call convert_ptc_to_bmad(xs,beta,k%time)
    endif
else


    if(c%cas==0) then
        CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    elseIF(c%cas==case1.or.c%cas==case2) then
       CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
    else
       IF(c%cas==caseP1) THEN
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
       ELSEif(c%cas==caseP2) THEN
          CALL TRACK_NODE_SINGLE(C,XS%X,K)  !,CHARGE
     ENDIF


    endif

endif
    xs%u=.not.check_stable
    if(xs%u) then
       lost_fibre=>c%parent_fibre
       lost_node=>c
       xlost=xs%x
    endif
       xs%last_node=>c
       xs%e=global_e
  END SUBROUTINE TRACK_NODE_FLAG_probe_P



  SUBROUTINE TRACK_FRINGE_spin_R(C,p,K)
    IMPLICIT NONE
    !    real(dp), INTENT(INOUT):: X(6),S(3)
    type(probe), INTENT(INOUT):: p
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    !    TYPE(ELEMENT), POINTER :: EL
    integer pos
    if(.not.(k%SPIN)) return

    !    el=>C%PARENT_FIBRE%MAG
    IF(.NOT.CHECK_STABLE) return
    IF(C%PARENT_FIBRE%MAG%KIND==KINDSUPERDRIFT) call superdrift_SPIN(c,P)
    if(C%PARENT_FIBRE%dir==1) then
       IF(C%CAS==CASE1) THEN
          call TRACK_rotate_spin(C,p,K)

          if(.not.C%parent_fibre%mag%p%kill_ent_spin) call TRACK_FRINGE_multipole(C,p,K)
          call TRACK_wedge_spin(C,p,K)
       else
          call TRACK_wedge_spin(C,p,K)

          if(.not.C%parent_fibre%mag%p%kill_exi_spin) call TRACK_FRINGE_multipole(C,p,K)

          call TRACK_rotate_spin(C,p,K)

       endif
    else
      ! write(6,*) " TRACK_FRINGE_spin_R "
       IF(C%CAS==CASE1) THEN
          call TRACK_rotate_spin(C,p,K)
          if(.not.C%parent_fibre%mag%p%kill_exi_spin) call TRACK_FRINGE_multipole(C,p,K)
          call TRACK_wedge_spin(C,p,K)
       else
          call TRACK_wedge_spin(C,p,K)
          if(.not.C%parent_fibre%mag%p%kill_ent_spin) call TRACK_FRINGE_multipole(C,p,K)
          call TRACK_rotate_spin(C,p,K)
       endif
     !  stop 888
    endif
  end SUBROUTINE TRACK_FRINGE_spin_R

  SUBROUTINE TRACK_FRINGE_spin_p(C,p,K)
    IMPLICIT NONE
    TYPE(probe_8), INTENT(INOUT)::p
    !    TYPE(REAL_8), INTENT(INOUT):: X(6),S(3)
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    !    TYPE(ELEMENTP), POINTER :: EL
    integer pos


    if(.not.(k%SPIN)) return
    !    el=>C%PARENT_FIBRE%MAGp
    IF(.NOT.CHECK_STABLE) return
    IF(C%PARENT_FIBRE%MAG%KIND==KINDSUPERDRIFT) call superdrift_SPIN(c,P)

    if(C%PARENT_FIBRE%dir==1) then
       IF(C%CAS==CASE1) THEN
          call TRACK_rotate_spin(C,p,K)
          if(.not.C%parent_fibre%magp%p%kill_ent_spin) call TRACK_FRINGE_multipole(C,p,K)
          call TRACK_wedge_spin(C,p,K)
       else
          call TRACK_wedge_spin(C,p,K)

          if(.not.C%parent_fibre%magp%p%kill_exi_spin) call TRACK_FRINGE_multipole(C,p,K)

          call TRACK_rotate_spin(C,p,K)

       endif
    else
       IF(C%CAS==CASE1) THEN
          call TRACK_rotate_spin(C,p,K)
          if(.not.C%parent_fibre%magp%p%kill_exi_spin) call TRACK_FRINGE_multipole(C,p,K)
          call TRACK_wedge_spin(C,p,K)
       else
          call TRACK_wedge_spin(C,p,K)
          if(.not.C%parent_fibre%magp%p%kill_ent_spin) call TRACK_FRINGE_multipole(C,p,K)
          call TRACK_rotate_spin(C,p,K)
       endif
      ! write(6,*) " TRACK_FRINGE_spin_p "
      ! stop 888
    endif

  end SUBROUTINE TRACK_FRINGE_spin_p

  SUBROUTINE TRACK_wedge_spin_R(C,p,K)
    IMPLICIT NONE
    !  this is a fake wedge.....
    type(probe), INTENT(INOUT):: p
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENT), POINTER :: EL
    integer pos
    el=>C%PARENT_FIBRE%MAG

    SELECT CASE(EL%KIND)
    case(KIND10)

       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(p,-C%PARENT_FIBRE%MAG%P%EDGE(1))
       ELSE
          CALL rot_spin_y(p,-C%PARENT_FIBRE%MAG%P%EDGE(2))
       ENDIF

       !       IF(C%CAS==CASE1) THEN
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(1))
       !       ELSE
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(2))
       !       ENDIF
       !       IF(C%CAS==CASE1) THEN
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(1))
       !       ELSE
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(2))
       !       ENDIF
       !    case(KIND20)
       !       IF(C%CAS==CASE1) THEN
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       !       ELSE
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       !       ENDIF

    END SELECT


  END SUBROUTINE TRACK_wedge_spin_R


  SUBROUTINE TRACK_wedge_spin_p(C,p,K)
    IMPLICIT NONE
    !  this is a fake wedge.....
    type(probe_8), INTENT(INOUT):: p
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENTP), POINTER :: EL
    integer pos
    el=>C%PARENT_FIBRE%MAGp

    SELECT CASE(EL%KIND)
    case(KIND10)

       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(p,-C%PARENT_FIBRE%MAGp%P%EDGE(1))
       ELSE
          CALL rot_spin_y(p,-C%PARENT_FIBRE%MAGp%P%EDGE(2))
       ENDIF

       !       IF(C%CAS==CASE1) THEN
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(1))
       !       ELSE
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(2))
       !       ENDIF
       !       IF(C%CAS==CASE1) THEN
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(1))
       !       ELSE
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%EDGE(2))
       !       ENDIF
       !    case(KIND20)
       !       IF(C%CAS==CASE1) THEN
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       !       ELSE
       !          CALL rot_spin_y(p,C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/TWO)
       !       ENDIF

    END SELECT


  END SUBROUTINE TRACK_wedge_spin_p


  SUBROUTINE TRACK_rotate_spin_r(C,p,K)
    IMPLICIT NONE
    !    real(dp), INTENT(INOUT):: X(6),S(3)
    type(probe), INTENT(INOUT):: p
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENT), POINTER :: EL
    integer pos
    el=>C%PARENT_FIBRE%MAG

    SELECT CASE(EL%KIND)
    case(KIND16,KIND10)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%EDGE(1))
       ELSE
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%EDGE(2))
       ENDIF
    case(KIND20)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/2.0_dp)
       ELSE
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/2.0_dp)
       ENDIF

    case(KINDPA)
       if(el%pa%hc==0.0_dp) then

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,el%p%dir*el%pa%angc)
            else
                CALL rot_spin_y(p,el%p%dir*el%pa%angc)
            endif
        else

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,-el%p%dir*el%pa%angc)
            else
                CALL rot_spin_y(p,-el%p%dir*el%pa%angc)
            endif

       endif
    case(KINDabell)
       if(el%ab%hc==0.0_dp) then

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,el%p%dir*el%ab%angc)
            else
                CALL rot_spin_y(p,el%p%dir*el%ab%angc)
            endif
        else

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,-el%p%dir*el%ab%angc)
            else
                CALL rot_spin_y(p,-el%p%dir*el%ab%angc)
            endif

       endif
    END SELECT


  END SUBROUTINE TRACK_rotate_spin_R

  SUBROUTINE TRACK_rotate_spin_p(C,p,K)
    IMPLICIT NONE
    TYPE(probe_8), INTENT(INOUT)::p
    !    TYPE(REAL_8), INTENT(INOUT):: X(6),S(3)
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENTP), POINTER :: EL
    integer pos
    el=>C%PARENT_FIBRE%MAGP

    SELECT CASE(EL%KIND)
    case(KIND16,KIND10)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%EDGE(1))
       ELSE
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%EDGE(2))
       ENDIF
    case(KIND20)
       IF(C%CAS==CASE1) THEN
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/2.0_dp)
       ELSE
          CALL rot_spin_y(p,el%p%dir*C%PARENT_FIBRE%MAG%P%B0*C%PARENT_FIBRE%MAG%P%LD/2.0_dp)
       ENDIF
    case(KINDPA)
       if(el%pa%hc==0.0_dp) then

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,el%p%dir*el%pa%angc)
            else
                CALL rot_spin_y(p,el%p%dir*el%pa%angc)
            endif
        else

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,-el%p%dir*el%pa%angc)
            else
                CALL rot_spin_y(p,-el%p%dir*el%pa%angc)
            endif

       endif
    case(KINDabell)
       if(el%ab%hc==0.0_dp) then

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,el%p%dir*el%ab%angc)
            else
                CALL rot_spin_y(p,el%p%dir*el%ab%angc)
            endif
        else

            IF(C%CAS==CASE1) THEN
                CALL rot_spin_y(p,-el%p%dir*el%ab%angc)
            else
                CALL rot_spin_y(p,-el%p%dir*el%ab%angc)
            endif

       endif

    END SELECT


  END SUBROUTINE TRACK_rotate_spin_p



  SUBROUTINE TRACK_FRINGE_multipole_R(C,p,K)
    IMPLICIT NONE
    !    real(dp), INTENT(INOUT):: X(6),S(3)
    type(probe), INTENT(INOUT):: p
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENT), POINTER :: EL
    integer pos
    el=>C%PARENT_FIBRE%MAG
    !    IF(.not.(k%FRINGE.or.el%p%permfringe)) return

    SELECT CASE(EL%KIND)
    CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       !    case(KIND2)
    case(KIND4)
    case(KIND2,kind5:kind7,KIND16:kind17,KIND20,KIND10,kindabell) ! Straight for all practical purposes
       IF(C%CAS==CASE1) THEN
          pos=-2
          !          call PUSH_SPIN_fake_fringe(c,p,my_true,k,pos)
          if(.not.el%P%KILL_ENT_spin) call PUSH_SPIN_fake_fringe(c,p,k,pos)
       elseif(C%CAS==CASE2) then
          pos=-1
          !          call PUSH_SPIN_fake_fringe(c,p,my_false,k,pos)

          if(.not.el%P%KILL_exi_spin) call PUSH_SPIN_fake_fringe(c,p,k,pos)

       endif
       !    case(KIND6)
       !    case(KIND7)
       !    case(KIND10)
       !    case(KIND16)
       !    case(KIND20)
    case(KIND21,kind22)
    case(KINDWIGGLER)
    case(KINDPA)
    case(kindsuperdrift)
    CASE DEFAULT
       WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
       stop 666
    END SELECT


  END SUBROUTINE TRACK_FRINGE_multipole_R

  SUBROUTINE TRACK_FRINGE_multipole_p(C,p,K)
    IMPLICIT NONE
    TYPE(probe_8), INTENT(INOUT)::p
    !    TYPE(REAL_8), INTENT(INOUT):: X(6),S(3)
    TYPE(INTERNAL_STATE) K
    TYPE (INTEGRATION_NODE), POINTER :: C
    TYPE(ELEMENTP), POINTER :: EL
    integer pos

    el=>C%PARENT_FIBRE%MAGP
    !    IF(.not.(k%FRINGE.or.el%p%permfringe)) return

    SELECT CASE(EL%KIND)
    CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       !    case(KIND2)
    case(KIND4)
    case(KIND2,kind5:kind7,KIND16:kind17,KIND20,KIND10,kindabell) ! Straight for all practical purposes
       IF(C%CAS==CASE1) THEN
          pos=-2
          !          call PUSH_SPIN_fake_fringe(c,p,my_true,k,pos)
          if(.not.el%P%KILL_ENT_spin) call PUSH_SPIN_fake_fringe(c,p,k,pos)
       elseif(C%CAS==CASE2) then
          pos=-1
          !          call PUSH_SPIN_fake_fringe(c,p,my_false,k,pos)

          if(.not.el%P%KILL_exi_spin) call PUSH_SPIN_fake_fringe(c,p,k,pos)

       endif
       !    case(KIND6)
       !    case(KIND7)
       !    case(KIND10)
       !    case(KIND16)
       !    case(KIND20)
    case(KIND21,kind22)
    case(KINDWIGGLER)
    case(KINDPA)
    case(kindsuperdrift)
    CASE DEFAULT
       WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
       stop 666
    END SELECT


  END SUBROUTINE TRACK_FRINGE_multipole_p

  SUBROUTINE TRACK_SPIN_FRONTR(C,P)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(PROBE), INTENT(INOUT) :: P
    !    real(dp), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE

    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

    ! PUSHING BEAM
    !




    ! The chart frame of reference is located here implicitely
    !    fake frontal spin snake PATCHG==5
    IF((PATCHG==1).or.(PATCHG==3).or.(PATCHG==5)) THEN
       CALL PATCH_SPIN(C,P,MY_TRUE)
    ENDIF


    CALL DTILT_SPIN(C%MAG%P%TILTD,1,P)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,P,my_true)
    ENDIF

  END SUBROUTINE TRACK_SPIN_FRONTR

  SUBROUTINE TRACK_SPIN_FRONTP(C,P)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(PROBE_8), INTENT(INOUT) :: P
    !    TYPE(REAL_8), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE

    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

    ! PUSHING BEAM
    !




    ! The chart frame of reference is located here implicitely
    !    fake frontal spin snake PATCHG==5
    IF((PATCHG==1).or.(PATCHG==3).or.(PATCHG==5)) THEN
       CALL PATCH_SPIN(C,P,MY_TRUE)
    ENDIF


    CALL DTILT_SPIN(C%MAG%P%TILTD,1,P)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,P,my_true)
    ENDIF

  END SUBROUTINE TRACK_SPIN_FRONTP

  !  SUBROUTINE TRACK_SPIN_FRONT_RAY8(C,S)
  !    implicit none
  !    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
  !    type(probe_8),INTENT(INOUT) ::S
  !    integer i,j
  !    type(real_8) sp(3)


  !    call TRACK_SPIN_front(C,S%S%X)


  !  end subroutine TRACK_SPIN_FRONT_RAY8


  ! back patch/misaglinments/tilt

  SUBROUTINE TRACK_SPIN_BACKR(C,P)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(PROBE), INTENT(INOUT) :: P
    !    REAL(DP), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE

    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,P,my_false)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILT_SPIN(C%MAG%P%TILTD,2,P)


    !    fake back spin snake PATCHG==6
    IF((PATCHG==2).or.(PATCHG==3).or.(PATCHG==6)) THEN
       CALL PATCH_SPIN(C,P,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely


  END SUBROUTINE TRACK_SPIN_BACKR

  !  SUBROUTINE TRACK_SPIN_BACK_RAY8(C,S)
  !    implicit none
  !    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
  !    type(probe_8),INTENT(INOUT) ::S
  !    integer i,j


  !    call TRACK_SPIN_BACK(C,S%S%X)
  !
  !
  !  end subroutine TRACK_SPIN_BACK_RAY8



  SUBROUTINE TRACK_SPIN_BACKP(C,P)
    implicit none
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    TYPE(PROBE_8), INTENT(INOUT) :: P
    !    TYPE(REAL_8), INTENT(INOUT) :: S(3)
    INTEGER(2) PATCHT,PATCHG,PATCHE



    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF


    IF(C%MAG%MIS) THEN
       CALL MIS_SPIN(C,P,my_false)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILT_SPIN(C%MAG%P%TILTD,2,P)


    !    fake back spin snake PATCHG==6
    IF((PATCHG==2).or.(PATCHG==3).or.(PATCHG==6)) THEN
       CALL PATCH_SPIN(C,P,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely


  END SUBROUTINE TRACK_SPIN_BACKP



  SUBROUTINE superdrift_SPINR(int,P)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(integration_node),target, INTENT(INOUT):: int
    TYPE(FIBRE),pointer:: c
    TYPE(PROBE), INTENT(INOUT):: P
    real(dp) da
    !    real(dp), INTENT(INOUT):: s(3)
      c=>int%parent_fibre
      if(c%dir==1.and.int%cas==case1) then
       da=C%mag%sdr%ANG(1)+((C%mag%sdr%A_X1-1)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_y(P,c%dir*C%mag%sdr%ANG(2)) ! 2016_5_9
       call rot_spin_z(P,C%mag%sdr%ANG(3))
       da=((C%mag%sdr%A_X2-1)/2)*pi
       call rot_spin_x(P,da)
      elseif(c%dir==-1.and.int%cas==case2) then
       da=((C%mag%sdr%A_X2-1)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_z(P,C%mag%sdr%ANG(3))
       call rot_spin_y(P,c%dir*C%mag%sdr%ANG(2))
       da=C%mag%sdr%ANG(1)+((C%mag%sdr%A_X1-1)/2)*pi
       call rot_spin_x(P,da)
      endif

  END SUBROUTINE superdrift_SPINR

  SUBROUTINE superdrift_SPINP(int,P)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(integration_node),target, INTENT(INOUT):: int
    TYPE(FIBRE),pointer:: c
    TYPE(PROBE_8), INTENT(INOUT):: P
    real(dp) da
       c=>int%parent_fibre
      if(c%dir==1.and.int%cas==case1) then
       da=C%mag%sdr%ANG(1)+((C%mag%sdr%A_X1-1)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_y(P,c%dir*C%mag%sdr%ANG(2)) ! 2016_5_9
       call rot_spin_z(P,C%mag%sdr%ANG(3))
       da=((C%mag%sdr%A_X2-1)/2)*pi
       call rot_spin_x(P,da)
      elseif(c%dir==-1.and.int%cas==case2) then
       da=((C%mag%sdr%A_X2-1)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_z(P,C%mag%sdr%ANG(3))
       call rot_spin_y(P,c%dir*C%mag%sdr%ANG(2))
       da=C%mag%sdr%ANG(1)+((C%mag%sdr%A_X1-1)/2)*pi
       call rot_spin_x(P,da)
      endif

  END SUBROUTINE superdrift_SPINP

  SUBROUTINE PATCH_SPINR(C,P,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(PROBE), INTENT(INOUT):: P
    !    real(dp), INTENT(INOUT):: s(3)
    logical(lp),INTENT(IN):: ENTERING
    real(dp) da
    IF(ENTERING) THEN
       da=C%PATCH%A_ANG(1)+((C%PATCH%A_X1-1)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_y(P,c%dir*C%PATCH%A_ANG(2)) ! 2016_5_9
       call rot_spin_z(P,C%PATCH%A_ANG(3))
       da=((C%PATCH%A_X2-1)/2)*pi
       call rot_spin_x(P,da)
    ELSE
       da=C%PATCH%B_ANG(1)+((C%PATCH%B_X1-1)/2)*pi
       call rot_spin_x(P,da)
       ! error etienne
       !      call rot_spin_y(P,C%PATCH%A_ANG(2))
       call rot_spin_y(P,c%dir*C%PATCH%b_ANG(2)) ! 2016_5_9
       call rot_spin_z(P,C%PATCH%b_ANG(3))
       da=((C%PATCH%B_X2-1)/2)*pi
       call rot_spin_x(P,da)
    ENDIF

  END SUBROUTINE PATCH_SPINR

  SUBROUTINE PATCH_SPINp(C,P,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    type(PROBE_8), INTENT(INOUT)::P
    !    type(real_8), INTENT(INOUT):: s(3)
    logical(lp),INTENT(IN):: ENTERING
    real(dp) da

    IF(ENTERING) THEN
       da=C%PATCH%A_ANG(1)+((C%PATCH%A_X1-1)/2)*pi

       call rot_spin_x(P,da)
       call rot_spin_y(P,c%dir*C%PATCH%A_ANG(2)) ! 2016_5_9
       call rot_spin_z(P,C%PATCH%A_ANG(3))


       da=((C%PATCH%A_X2-1)/2)*pi

       call rot_spin_x(P,da)
    ELSE
       da=C%PATCH%B_ANG(1)+((C%PATCH%B_X1-1)/2)*pi
       call rot_spin_x(P,da)
       ! error etienne
       !       call rot_spin_y(P,C%PATCH%A_ANG(2))
       call rot_spin_y(P,c%dir*C%PATCH%b_ANG(2)) ! 2016_5_9
       call rot_spin_z(P,C%PATCH%b_ANG(3))


       da=((C%PATCH%B_X2-1)/2)*pi

       call rot_spin_x(P,da)
    ENDIF


  END SUBROUTINE PATCH_SPINp

  !   Misalignment routines
  SUBROUTINE MIS_SPINR(C,P,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(PROBE), INTENT(INOUT):: P
    !    real(dp), INTENT(INOUT):: S(3)
    logical(lp),INTENT(IN):: ENTERING

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             call rot_spin_X(P,C%CHART%ANG_IN(1))
             call rot_spin_Y(P,C%CHART%ANG_IN(2))
             call rot_spin_Z(P,C%CHART%ANG_IN(3))
          ELSE
             call rot_spin_X(P,C%CHART%ANG_OUT(1))
             call rot_spin_Y(P,C%CHART%ANG_OUT(2))
             call rot_spin_Z(P,C%CHART%ANG_OUT(3))
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
             call rot_spin_Z(P,C%CHART%ANG_OUT(3))
             call rot_spin_Y(P,-C%CHART%ANG_OUT(2))   !2016_5_9
             call rot_spin_X(P,C%CHART%ANG_OUT(1))
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
             call rot_spin_Z(P,C%CHART%ANG_IN(3))
             call rot_spin_Y(P,-C%CHART%ANG_IN(2))   !2016_5_9
             call rot_spin_X(P,C%CHART%ANG_IN(1))
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_SPINR

  SUBROUTINE MIS_SPINP(C,P,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),INTENT(INOUT):: C
    TYPE(PROBE_8), INTENT(INOUT):: P
    logical(lp),INTENT(IN):: ENTERING

    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             call rot_spin_X(P,C%CHART%ANG_IN(1))
             call rot_spin_Y(P,C%CHART%ANG_IN(2))
             call rot_spin_Z(P,C%CHART%ANG_IN(3))
          ELSE
             call rot_spin_X(P,C%CHART%ANG_OUT(1))
             call rot_spin_Y(P,C%CHART%ANG_OUT(2))
             call rot_spin_Z(P,C%CHART%ANG_OUT(3))
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)


             call rot_spin_Z(P,C%CHART%ANG_OUT(3))
             call rot_spin_Y(P,-C%CHART%ANG_OUT(2))   !2016_5_9
             call rot_spin_X(P,C%CHART%ANG_OUT(1))

             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)

             call rot_spin_Z(P,C%CHART%ANG_IN(3))
             call rot_spin_Y(P,-C%CHART%ANG_IN(2))    !2016_5_9
             call rot_spin_X(P,C%CHART%ANG_IN(1))

             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_SPINP


  SUBROUTINE DTILT_SPINR(TILTD,I,P)
    IMPLICIT NONE
    TYPE(PROBE),INTENT(INOUT):: P
    !    real(dp),INTENT(INOUT):: S(3)
    INTEGER,INTENT(IN):: I
    REAL(DP),INTENT(IN) :: TILTD
    real(dp) YS

    IF(TILTD==0.0_dp) RETURN
    IF(I==1) THEN
       YS=TILTD
       call rot_spin_Z(P,YS)
    ELSE
       YS=-TILTD
       call rot_spin_Z(P,YS)
    ENDIF

  END SUBROUTINE DTILT_SPINR

  SUBROUTINE DTILT_SPINP(TILTD,I,P)
    IMPLICIT NONE
    TYPE(PROBE_8),INTENT(INOUT):: P
    !    TYPE(REAL_8),INTENT(INOUT):: S(3)
    INTEGER,INTENT(IN):: I
    REAL(DP),INTENT(IN) :: TILTD
    real(dp) YS

    IF(TILTD==0.0_dp) RETURN
    IF(I==1) THEN
       YS=TILTD
       call rot_spin_Z(P,YS)
    ELSE
       YS=-TILTD
       call rot_spin_Z(P,YS)
    ENDIF

  END SUBROUTINE DTILT_SPINP

 !   backward compatible routine for radiation

  SUBROUTINE stroboscopic_average(ring,xs0,xst,pos,mstate0,nturn,kp,n,mf)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    type(probe) , intent(inout) :: xs0,xst
    TYPE(INTERNAL_STATE) mstate0,mstate
    integer, intent(in) :: kp,pos,nturn
    integer, optional :: mf
    integer i,k,imax,nd2,mff
    type(spinor) n
    real(dp) norm,norm0,n0(3),theta0

    ! xs0 => ray being tracked
    ! xst = > average
    ! pos in layout
    ! nturn  = > number of turns
    ! kp => frequency of printing on screen while tracking
    !  spinor where info is stored for ISF

    mff=6
    if(present(mf)) mff=mf
    mstate=mstate0+spin0
    nd2=6
    if(mstate%nocavity) nd2=4

    write(mff,*); write(mff,*) " Results of Stroboscopic Averaging "
    write(mff,*) " every ",kp," turns "
    do k=1,nturn
       call track_probe(ring,xs0,mstate,node1=pos)  !,fibre2=3)
  if(use_quaternion) call probe_quaternion_to_matrix(xs0)
       do i=1,3
          xst%s(i)%x=xs0%s(i)%x+xst%s(i)%x  ! <---- Stroboscopic average
       enddo

       if(mod(k,kp)==0) then  ! kp
          write(mff,*) k,"#########################"
          do i=1,3
             norm=root(xst%s(1)%x(i)**2+xst%s(2)%x(i)**2+xst%s(3)%x(i)**2)
             if(norm>=0.0_dp) then
                norm=1.d0/norm
             endif
             write(mff,'(a14,4(1x,g20.13))') " Norm and ISF ", 1.d0/norm, xst%s(1)%x(i)*norm,xst%s(2)%x(i)*norm,xst%s(3)%x(i)*norm
          enddo
       endif ! kp
    enddo

    norm0=0.0_dp
    do i=1,3
       norm=root(xst%s(1)%x(i)**2+xst%s(2)%x(i)**2+xst%s(3)%x(i)**2)
       if(norm>=norm0) then
          imax=1
          n%x(1)=xst%s(1)%x(i)/norm
          n%x(2)=xst%s(2)%x(i)/norm
          if(abs(n%x(2))>=abs(n%x(1))) imax=2
          n%x(3)=xst%s(3)%x(i)/norm
          if(abs(n%x(3))>=abs(n%x(2))) imax=3
          if(n%x(imax)<0) n%x=-n%x
          norm0=norm
       endif
    enddo

  end SUBROUTINE stroboscopic_average


  ! time tracking


   SUBROUTINE TRACK_time(xT,DTt,K)
    ! Tracks a single particle   of the beam for a time DT
    implicit none
    TYPE(INTEGRATION_NODE), POINTER:: T
    TYPE(temporal_probe),INTENT(INOUT):: xT
    REAL(DP), INTENT(IN) :: DTt
    REAL(DP) DT0,DT_BEFORE,X(6),fac,dt
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    type(probe) b
    LOGICAL(LP) END_OF_LINE
    END_OF_LINE=.FALSE.

    CALL RESET_APERTURE_FLAG

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



    XT%r=0.0_dp
    DT0=0.d0
    dt=dtt+xt%dt0
!    XT%xs%X=XT%Xb

    FAC=0.0_dp

    DO WHILE(DT0<=DT)

       DT_BEFORE=DT0
       b=XT%xs
       !         WRITE(6,*) " POS ",T%s(1),t%pos_in_fibre
       !         WRITE(6,*) " POS ",T%POS,T%CAS,T%PARENT_FIBRE%MAG%NAME
       ! putting spin and radiation
       CALL TRACK_NODE_PROBE(t,XT%XS,K) !,charge)
       ! no spin and no radiation
       !        CALL TRACK_NODE_SINGLE(t,XT%XS%X,K)  !,CHARGE

       !
       !
       DT0=DT0+(XT%xs%X(6)-b%x(6))
!       write(6,*) dt0,dt,t%pos,t%pos_in_fibre
       T=>T%NEXT
       IF(.NOT.ASSOCIATED(T%NEXT)) THEN
          END_OF_LINE=.TRUE.
          EXIT
       ENDIF
       if(.not.check_stable.or.XT%xs%u) exit

    ENDDO
       XT%NODE=>T%previous

    IF(.NOT.END_OF_LINE.and.(.not.XT%xs%u)) THEN
       IF(DT0/=DT.and.check_stable) THEN
  !        XT%NODE=>T%previous
          XT%Dt0=dt-DT_BEFORE
          fac=(dt-DT_BEFORE)/(dt0-DT_BEFORE)
    !      X=fac*(XT%xs%X-XT%Xb)+XT%Xb
     !  yl=fac*XT%NODE%parent_fibre%mag%l/XT%NODE%parent_fibre%mag%p%nst
       ENDIF
    ELSE


       IF(DT0<DT.and.check_stable.and.(.not.XT%xs%u)) THEN
  !        XT%NODE=>T%previous
          XT%Dt0=dt-DT_BEFORE
          fac=(dt-DT_BEFORE)/(dt0-DT_BEFORE)
  !        X=fac*(XT%xs%X-XT%Xb)+XT%Xb
  !     yl=fac*XT%NODE%parent_fibre%mag%l/XT%NODE%parent_fibre%mag%p%nst
       ENDIF

    ENDIF


    XT%r=FAC

    if(.not.CHECK_STABLE.or.XT%xs%u) then
       lost_fibre=>t%parent_fibre
       lost_node=>t
       XT%xs%u=.true.
    endif

    call ptc_global_x_p(xt,b,k)

!call make_spinor_basis(xt%s(1),xt%s(2),xt%s(3))

!write(6,*) xt%s(1).dot.xt%s(1),xt%s(2).dot.xt%s(2),xt%s(3).dot.xt%s(3)
!write(6,*) xt%s(1).dot.xt%s(2),xt%s(2).dot.xt%s(3),xt%s(1).dot.xt%s(3)
    XT%t=b%x(6)+XT%dt0
    xt%xs=b
  END SUBROUTINE TRACK_time

  Subroutine ptc_global_x_p(xt,b,k)
    implicit none
    TYPE(temporal_probe),INTENT(INOUT):: xT
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer i,posi
    type(probe) b
    real(dp) p(3),pzb,pz,betinv,pos(3),dab(3),xp(2),xpa(2),e(3),eb(3),dal(3)

    if(associated(XT%NODE%ENT)) then
       ! computing global position
       call DIRECTION_V(k,xt%node%parent_fibre%mag,xt%xs%x,E,XP,XPA,xt%node%POS_IN_FIBRE-2)
       call DIRECTION_V(k,xt%node%parent_fibre%mag,b%x,EB,XP,XPA,xt%node%POS_IN_FIBRE-2)
       if(k%spin) call find_frac_r(xt,b)

       XT%pos=0.0_dp
       pos=0.0_dp
       DO I=1,3
          XT%POS(i)=XT%POS(i) + XT%XS%X(1)*XT%NODE%EXI(1,I)     !
          XT%POS(i)=XT%POS(i) + XT%XS%X(3)*XT%NODE%EXI(2,I)     !

          pos(i)=POS(i) + b%x(1)*XT%NODE%ent(1,I)
          pos(i)=POS(i) + b%x(3)*XT%NODE%ent(2,I)
       ENDDO
       dal(1:3)=pos(1:3)
       XT%pos(1:3) = XT%pos(1:3) + XT%NODE%B
       pos(1:3) = pos(1:3) + XT%NODE%A
       dab(1:3)=XT%pos(1:3)-pos(1:3)




       dab(1:3)=XT%r*dab(1:3)


       dal(1:3)=dal(1:3)+dab(1:3)
   ! pz=dal(1)*XT%NODE%ent(3,1)+dal(2)*XT%NODE%ent(3,2)+dal(3)*XT%NODE%ent(3,3)

       call change_basis(dal(1:3),global_frame,xt%ic,XT%NODE%ent)

!call ptc_print(XT%NODE%a,XT%NODE%ent,6)

       XT%pos(1:3)=pos(1:3)+dab(1:3)



       ! computing global momentum
       if(k%time) then
          betinv=1.0_dp/XT%NODE%PARENT_FIBRE%beta0
       ELSE
          betinv=1.0_dp
       ENDIF
       pz=sqrt(1.0_dp+2.0_dp*betinv*XT%XS%X(5)+XT%XS%X(5)**2)*XT%NODE%PARENT_FIBRE%mag%p%p0c
       pzb=sqrt(1.0_dp+2.0_dp*betinv*b%X(5)+b%X(5)**2)*XT%NODE%PARENT_FIBRE%mag%p%p0c


       pz=pzb+XT%r*(pz-pzb)
       e(1:3)=eb(1:3)+XT%r*(e(1:3)-eb(1:3))
       XT%POS(4:6)=e(1:3)*pz/sqrt(e(1)**2+e(2)**2+e(3)**2)

    else
       write(6,*) " FILL_SURVEY_DATA_IN_NODE_LAYOUT "
       write(6,*) " was not called, so no survey data on the nodes "
       write(6,*) " program will stop inside TRACK_time "
       stop

    endif

  end Subroutine ptc_global_x_p

  SUBROUTINE  locate_temporal_probe(r,tp,sc) !
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    type(temporal_probe),intent(INOUT) ::  tp
    real(dp), optional :: sc
    type(integration_node), pointer :: t,tw
    real(dp) a(3)
    integer i,j
    real(dp) norm,d1,p(3),da(3),dal(3)

    type(temporal_probe) tt

    ! locating the closest integration node
    ! pointing at it using tw
    norm=mybig
    t=> tp%node
    a=tp%pos(1:3)
    do i=1,r%t%n
       if(t%cas==case0) then ! 1
          d1=root((t%a(1)-a(1))**2+(t%a(2)-a(2))**2+(t%a(3)-a(3))**2 )

          if(d1<norm.and.t%parent_fibre%mag%l/=0.0_dp) then
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
    d1=da(1)*tw%ent(3,1)+da(2)*tw%ent(3,2)+da(3)*tw%ent(3,3)

    call change_basis(DA,global_frame,dal,tw%ent)
    !
    !

    !######################################
    ! some gymnastic if behind the integration node


   !    if(dal(3)<0.0_dp) then

       do while(dal(3)<0.0_dp)

          tw=>tw%previous
          do while(tw%cas/=case0)
             tw=>tw%previous
          enddo
          da=a-tw%a
          call change_basis(DA,global_frame,dal,tw%ent)
      enddo
 !   endif

   if(present(sc)) then
    sc=da(1)*tw%ent(3,1)+da(2)*tw%ent(3,2)+da(3)*tw%ent(3,3)
 !   write(6,*) " a.ent(3)",d1
    endif
    p=dal

 !   b%tp(j)%r=p(3)
    tp%node=>tw
!write(6,*) " frame in locate_temporal_probe "
!call ptc_print(tp%NODE%a,tp%NODE%ent,6)
  !  call original_p_to_ptc(b,j,p,tw)

  end SUBROUTINE  locate_temporal_probe

SUBROUTINE  fit_temporal_probe(tp,kk,done)
implicit none
type(temporal_probe) tp
type(internal_state) kk
logical done

if((tp%node%parent_fibre%mag%kind==kind4.and.(.not.kk%nocavity)).or.kk%radiation.and.(.not.locate_with_no_cavity)) then
 call fit_temporal_probe_cav(tp,kk,done)
else
 call fit_temporal_probe_nocav(tp,kk,done)
endif

end SUBROUTINE  fit_temporal_probe

SUBROUTINE  fit_temporal_probe_nocav(tp,kk,done)
implicit none
type(temporal_probe) tp,tt
type(internal_state) k,kk
real(dp) del,betinv,x(6),eps,d1,d2,r,ma(5,5),mai(5,5),r0,pos0(6),de(5),norm
real(dp) epsn,normb
integer i,j,ier,l,m,nit
type(work) w
logical done
nit=1000
done=.true.
eps=1.d-8
epsn=1.d-9
normb=1.d38
!write(6,*) " eps "
!read(5,*) eps
call alloc_temporal_probe(tt)
tt=tp

k=kk+nocavity0

del= (tp%pos(4)**2+tp%pos(5)**2+tp%pos(6)**2)/tp%node%parent_fibre%mag%p%p0c**2-1.0_dp
       if(k%time) then
          betinv=1.0_dp/tp%node%parent_fibre%beta0
       ELSE
          betinv=1.0_dp
       ENDIF

d1=0.0_dp;d2=0.0_dp;

do i=1,3
 d1=d1+(tp%node%a(i)-tp%pos(i))**2;  d2=d2+(tp%node%b(i)-tp%pos(i))**2;
enddo

!d1=sqrt(d1); d2=sqrt(d2);
!w=tp%node%parent_fibre
!r0=d1/(d1+d2)*tp%node%parent_fibre%mag%L/w%beta0/tp%node%parent_fibre%mag%p%nst
r0=0


!write(6,*) "dt initial ",r0
!write(6,'(3(1x,f20.15))') tp%node%a
!write(6,'(3(1x,f20.15))') tp%node%b
!pause 7
x=0.d0

del=del/( betinv+sqrt(betinv**2+del) )

x(5)=del

!write(6,*) del
do l=1,nit
!101 continue
!write(6,'(6(1x,f20.15))') tp%pos

!tt=tp
tt%xs%x=x
tt%dt0=0
r=r0
  call TRACK_time(tt,r,K)
 !call ptc_global_x_p(tt,b,k)
pos0=tt%pos



do i=1,4
 tt=tp
 tt%xs%x=x
 tt%xs%x(i)=x(i)+eps
 tt%dt0=0
 r=r0

! write(6,'(6(1x,f20.15))')tt%xs%x
 call TRACK_time(tt,r,K)

! write(6,'(6(1x,f20.15))')tt%xs%x
 ! call ptc_global_x_p(tt,b,k)

! write(6,'(6(1x,f20.15))') tt%pos
do j=1,5
 ma(j,i)=(tt%pos(j)-pos0(j))/eps
enddo

!pause
enddo

tt%xs%x=x
tt%dt0=0
r=r0+eps


!write(6,'(6(1x,f20.15))')tt%xs%x
  call TRACK_time(tt,r,K)


!write(6,'(6(1x,f20.15))')tt%xs%x
 !call ptc_global_x_p(tt,b,k)

!write(6,'(6(1x,f20.15))') tt%pos
do j=1,5
 ma(j,5)=(tt%pos(j)-pos0(j))/eps
enddo

!write(6,*) " matrix "
!do i=1,5
!write(6,'(5(1x,f20.15))')ma(i,1:5)
!enddo

call matinv(ma,mai,5,5,ier)
!write(6,*) " matrix inverse "
!do i=1,5
!write(6,'(5(1x,f20.15))')mai(i,1:5)
!enddo

!ma=matmul(ma,mai)
!write(6,*) " identity ? "
!do i=1,5
!write(6,'(5(1x,f20.15))')ma(i,1:5)
!enddo


de = matmul(mai,(tp%pos(1:5) - pos0(1:5) ))




!write(6,'(5(1x,f20.15))') de

x(1:4)=x(1:4)+de(1:4)
r0=r0+de(5)

norm=0
do m=1,6
 norm=abs(pos0(m)-tp%pos(m))+norm
enddo

if(norm>epsn.or.l<10) then
 normb=norm
else
 if(normb>=norm) exit
 normb=norm
endif

enddo

tt%xs%x=x
tt%dt0=0
r=r0
  call TRACK_time(tt,r,K)
 !call ptc_global_x_p(tt,b,k)
pos0=tt%pos

tt%xs%x(6)=tp%t-tt%dt0

!write(6,*) " norm ",norm,l
!write(6,*) " tp%dt0,tp%r ", tp%dt0,tp%r
!write(6,'(6(1x,f20.15))') pos0
!write(6,'(6(1x,f20.15))') tp%pos
tp=tt
if(l> nit-10) then
 write(6,*) " Could not converge in fit_temporal_probe_nocav. Norm = ",norm
 done=.false.
endif
   end SUBROUTINE  fit_temporal_probe_nocav

SUBROUTINE  fit_temporal_probe_cav(tp,k,done)
implicit none
type(temporal_probe) tp,tt
type(internal_state) k
real(dp) del,betinv,x(6),eps,d1,d2,r,ma(7,7),mai(7,7),r0,pos0(6),de(7),norm,t0
real(dp) epsn,normb
integer i,j,ier,l,m,nit
type(work) w
logical done
done=.true.
eps=1.d-8
epsn=1.d-9
normb=1.d38
 nit=1000

call alloc_temporal_probe(tt)
tt=tp



del= (tp%pos(4)**2+tp%pos(5)**2+tp%pos(6)**2)/tp%node%parent_fibre%mag%p%p0c**2-1.0_dp
       if(k%time) then
          betinv=1.0_dp/tp%node%parent_fibre%beta0
       ELSE
          betinv=1.0_dp
       ENDIF

d1=0.0_dp;d2=0.0_dp;

do i=1,3
 d1=d1+(tp%node%a(i)-tp%pos(i))**2;  d2=d2+(tp%node%b(i)-tp%pos(i))**2;
enddo

!d1=sqrt(d1); d2=sqrt(d2);
!w=tp%node%parent_fibre
!r0=d1/(d1+d2)*tp%node%parent_fibre%mag%L/w%beta0/tp%node%parent_fibre%mag%p%nst
r0=0

x=0.d0
del=del/( betinv+sqrt(betinv**2+del) )
x(6)=tp%t
x(5)=del


do l=1,nit



tt%xs%x=x
tt%dt0=0
r=r0
  call TRACK_time(tt,r,K)

pos0=tt%pos
t0=tt%t


do i=1,6
 tt=tp
 tt%xs%x=x
 tt%xs%x(i)=x(i)+eps
 tt%dt0=0
 r=r0


 call TRACK_time(tt,r,K)

do j=1,6

 ma(j,i)=(tt%pos(j)-pos0(j))/eps
enddo
ma(7,i)=(tt%t-t0)/eps

enddo

tt%xs%x=x
tt%dt0=0
r=r0+eps

  call TRACK_time(tt,r,K)

do j=1,6
 ma(j,7)=(tt%pos(j)-pos0(j))/eps
enddo
 ma(7,7)=(tt%t-t0)/eps

call matinv(ma,mai,7,7,ier)

de(1:6)=tp%pos(1:6) - pos0(1:6)
de(7)=tp%t-t0

de = matmul(mai,de)


x(1:6)=x(1:6)+de(1:6)
r0=r0+de(7)

norm=abs(t0-tp%t)
do m=1,6
 norm=abs(pos0(m)-tp%pos(m))+norm
enddo

if(norm>epsn.or.l<10) then
 normb=norm

else

 if(normb>=norm) exit
 normb=norm
endif

enddo

tt%xs%x=x
tt%dt0=0
r=r0
  call TRACK_time(tt,r,K)
 !call ptc_global_x_p(tt,b,k)
pos0=tt%pos


tp=tt
if(l> nit-10) then
 write(6,*) " Could not converge in fit_temporal_probe_cav. Norm = ",norm
 done=.true.
endif
   end SUBROUTINE  fit_temporal_probe_cav




Subroutine ptc_print(a,ent,mf)
implicit none
integer i,mf
real(dp) ent(3,3),a(3)

write(mf,*)"a"
 write(mf,'(3(1x,g20.13))') a
write(6,*)
write(mf,*)"ent"
do i=1,3
 write(mf,'(3(1x,g20.13))') ent(i,1:3)
enddo
write(mf,*)
end subroutine ptc_print


  Subroutine find_frac_r(xt,b)
    implicit none
    type(temporal_probe) xt
    type(probe) b
    type(spinor)n0
    real(dp) m1t(3,3),m2(3,3),r(3,3),a(3,3),ai(3,3),m1(3,3),ang,rat,m2out(3,3)
    integer i,j

   do i=1,3
   do j=1,3
     m1t(i,j)=b%s(i)%x(j)
     m1(i,j)=b%s(j)%x(i)
     m2(i,j)=xt%xs%s(j)%x(i)
   enddo
   enddo
     r=matmul(m1t,m2)


 if(abs(abs(r(1,1))+abs(r(2,2))+abs(r(3,3))-3.0_dp)>1.d-14 ) then


call  find_n0(r,n0)
call find_as(n0,a,ai)

    r=matmul(matmul(ai,r),a)


   ang=atan2(r(1,3),r(1,1))*xt%r
   r=0.0_dp

   r(1,1)=cos(ang);r(1,3)=sin(ang);
   r(3,3)=cos(ang);r(3,1)=-sin(ang);
   r(2,2)=1.0_dp


   r=matmul(matmul(matmul(a,r),ai),m1)
else
   r=m2
endif





   do i=1,3
   do j=1,3
     xt%s(j)%x(i)=r(i,j)
   enddo
   enddo


   end Subroutine find_frac_r


  subroutine find_as(n22,a,ai)
!#general : normal & manipulation
!# Find the c_spinmatrix "a" such that
!# e_y = (0,1,0)= a**(-1)*n0
!# because a*exp(theta n0.L)*a**(-1)= exp(theta (a*n0).L). See Sec.6.5.2.

    implicit none
    type(spinor), intent(inout) ::  n22
    real(dp)   a(3,3) ,  ai(3,3)
    type(spinor)  n1 ,n3,n2
     real(dp)    s,n
    real(dp) x
    integer i,is,j



    n2=n22

    x=n2.dot.n2

    x=sqrt(x)

    do i=1,3
     n2%x(i)=n2%x(i)/x
    enddo

    ! here we find smallest value of n2
    is=2
    if(abs(n2%x(1))< abs(n2%x(2))) is=1

    if(is==1) then
       if(abs(n2%x(3))<abs(n2%x(1))) is=3
    else
       if(abs(n2%x(3))<abs(n2%x(2))) is=3
    endif

    !  put n1 in along that value
    do i=1,3
       n1%x(i)=0.0_dp
    enddo
    n1%x(is)=1.0_dp

    s=n2%x(is)*n1%x(is)

    n=0.0_dp
    do i=1,3
       n1%x(i)=n1%x(i)-s*n2%x(i)
       n=n1%x(i)**2+n
    enddo
    do i=1,3
       n1%x(i)=n1%x(i)/sqrt(n)
    enddo

    n3%x(1)=n1%x(2)*n2%x(3)-n1%x(3)*n2%x(2)
    n3%x(2)=n1%x(3)*n2%x(1)-n1%x(1)*n2%x(3)
    n3%x(3)=n1%x(1)*n2%x(2)-n1%x(2)*n2%x(1)

    n=0.0_dp
    do i=1,3
       n=n3%x(i)**2+n
    enddo
    do i=1,3
       n3%x(i)=n3%x(i)/sqrt(n)
    enddo

 !   if(spin_normal_position==2) then
 if(abs(n1%x(1))>abs(n3%x(1))) then
       do i=1,3
          a(i,1)=n1%x(i)
          a(i,2)=n2%x(i)
          a(i,3)=n3%x(i)
       enddo
else
x=n3%x(1)    !.sub.'0'

if(x<0) then
       do i=1,3
          a(i,1)=-n3%x(i)
          a(i,2)=n2%x(i)
          a(i,3)=n1%x(i)
       enddo
else
       do i=1,3
          a(i,1)=n3%x(i)
          a(i,2)=n2%x(i)
          a(i,3)=-n1%x(i)
       enddo
endif

endif

   do i=1,3
   do j=1,3
    ai(i,j)=a(j,i)
   enddo
   enddo
  end subroutine find_as

 subroutine find_n0(s0,n0)

    implicit none
    real(dp) ,intent(inout) :: s0(3,3)
    type(spinor), intent(inout) :: n0
    real(dp)  norm0
    real(dp) det,detm
    real(dp) ss(3,3)

    integer i,is,j

   ss=0.0_dp

!    ss=s0


   do i=1,3
       do j=1,3

         ss(i,j)=s0(i,j)

       enddo
    enddo



    do i=1,3
       ss(i,i)=ss(i,i)-1.0_dp
    enddo

    det=(ss(2,2)*ss(3,3)-ss(2,3)*ss(3,2))

    is=1
    detm=(ss(1,1)*ss(3,3)-ss(1,3)*ss(3,1))

    if(abs(detm)>=abs(det)) then
       det=detm
       is=2
    endif

    detm=ss(1,1)*ss(2,2)-ss(1,2)*ss(2,1)
    if(abs(detm)>=abs(det)) then
       det=detm
       is=3
    endif


    n0%x(is)=1.0_dp
    if(is==1) then
       n0%x(2)=(-ss(3,3)*ss(2,1)+ss(2,3)*ss(3,1))/det
       n0%x(3)=(-ss(2,2)*ss(3,1)+ss(2,1)*ss(3,2))/det
    elseif(is==2) then
       n0%x(1)=(-ss(3,3)*ss(1,2)+ss(3,2)*ss(1,3))/det
       n0%x(3)=(-ss(1,1)*ss(3,2)+ss(1,2)*ss(3,1))/det
    else
       n0%x(1)=(-ss(2,2)*ss(1,3)+ss(2,3)*ss(1,2))/det
       n0%x(2)=(-ss(1,1)*ss(2,3)+ss(1,3)*ss(2,1))/det
    endif

     norm0=sqrt(n0%x(1)**2+n0%x(2)**2+n0%x(3)**2)



    do i=1,3
       n0%x(i)=n0%x(i)/norm0
    enddo



  end subroutine find_n0

subroutine equal_temporal(xtt,xt)
implicit none
type (temporal_probe),INTENT(IN)::xt
type (temporal_probe), INTENT(inOUT)::xtt
integer i
       xtt%r=xt%r
       xtt%xs=xt%xs
       xtt%dt0=xt%dt0
       xtt%pos=xt%pos
       xtt%ic=xt%ic
       xtt%t=xt%t
    do i=1,3
       xtt%s(i)%x=xt%s(i)%x
    enddo

xtt%node=>xt%node

end subroutine equal_temporal

  SUBROUTINE  alloc_temporal_beam(b,n,p0c) ! fibre i1 to i2
    IMPLICIT NONE
    type(temporal_beam), intent(inout):: b
    integer i,n
    real(dp) p0c


    allocate(b%tp(n))
    b%n=n
    b%a=GLOBAL_origin
    b%ent=global_frame
    b%total_time=0.0_dp
    b%p0c=p0c
    nullify(b%c)

    do i=1,n
       call alloc(b%tp(i))
    enddo

  END SUBROUTINE  alloc_temporal_beam

  SUBROUTINE  alloc_temporal_probe(p) ! fibre i1 to i2
    IMPLICIT NONE
    type(temporal_probe), intent(inout):: p

    p%xs%u=my_false
    p%xs%x=0.0_dp
    p%xs%s(1)=0
    p%xs%s(2)=0
    p%xs%s(3)=0
    p%r=0.0_dp
    p%pos=0.0_dp
    p%dt0=0.0_dp
    p%s(1)%x=0.0_dp
    p%s(2)%x=0.0_dp
    p%s(3)%x=0.0_dp
    p%IC=0.0_dp
    p%T=0.0_dp
    nullify(p%node)
    nullify(p%xs%last_node)

  end SUBROUTINE  alloc_temporal_probe


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


!!!! Routines to create maps for tracking out of PTC !!!!
!!! this fibre must be part of a layout and a thin layout as well



subroutine fill_tree_element(f,no,fix0,onemap,factor,file)   ! fix0 is the initial condition for the maps
implicit none
type(fibre), target :: f
type(layout), pointer :: r
TYPE(INTEGRATION_NODE),POINTER:: t1c,t2c
TYPE (NODE_LAYOUT), POINTER :: t
TYPE (tree_element), POINTER :: arbre(:)
type(internal_state) state
real(dp) fixr(6),fixs(6),fix(6),fix0(6),mat(6,6),e_ij(6,6),xn
type(probe) xs0
type(probe_8) xs
type(c_damap) m,mr
logical :: onemap,fact
logical,optional :: factor
integer no,i,mf
character(*), optional :: file 

fact=.false.

if(present(factor)) fact=factor


if(.not.associated(f%parent_layout)) then
 write(6,*) " parent layout not associated "
 stop
else
 r=>f%parent_layout
endif

if(.not.associated(f%parent_layout%t)) then
 write(6,*) " parent node layout not associated "
 stop
else
 t=>f%parent_layout%t
 t1c=>f%t1%next
 t2c=>f%t2   !%next
endif

! Classical radiation with stochastic envelope

state=radiation0+envelope0+time0

call init_all(state,1,0)

call alloc(xs);call alloc(m,mr)


! radiation

xs0=fix0
mr=1
xs=xs0+mr
call propagate(xs,state,node1=t1c,node2=t2c)


! For David
!!  mr: linear map with radiation would be read here instead of being computed,
!! and must be stored in fixr

fixr=xs%x    ! <---
mr=xs   ! <---

do i=1,6
 mr%v(i)=mr%v(i)-(mr%v(i).sub.0)
enddo

! For David
!!  The stochastic kicks are stored at e_ij

e_ij=xs%e_ij         ! <---
! no radiation


state=time0
xs0=fix0
m=1
xs=xs0+m
!write(6,*) t1c%parent_fibre%mag%name,t1c%parent_fibre%mag%p%nst
call propagate(xs,state,node1=t1c,node2=t2c)
fix=xs%x
! For David
!!  The same linear map is computed WITHOUT radiation : result put into m, the constant part is removed
!!
m=xs   ! <---
do i=1,6
 m%v(i)=m%v(i)-(m%v(i).sub.0)
enddo

m=m**(-1)*mr

if(.not.onemap) then
 call  nth_root(m,m,f%mag%p%nst)
endif

mat=m


call kill(xs);call kill(m);call kill(mr)

state=spin0+time0
call init_all(state,no,0)
call alloc(xs);call alloc(m)



xs0=fix0
 
m=1
xs=xs0+m
 

call propagate(xs,state,node1=t1c,node2=t2c)


! For David
!!  The full nonlinear map m is computed and the final orbit
!!
fix=xs%x  ! <---
m=xs  ! <---

m%e_ij=e_ij
do i=1,6
 m%v(i)=m%v(i)-(m%v(i).sub.0)
enddo


if(.not.onemap) then
 call  nth_root(m,m,f%mag%p%nst)
endif


if(f%dir==1) then
 if(.not.associated(f%mag%forward)) then
  allocate(f%mag%forward(3))
 ! allocate(f%mag%usef)
 else
  call KILL(f%mag%forward)
 endif
 
call SET_TREE_G_complex(f%mag%forward,m,fact)
 f%mag%do1mapf=onemap
 f%mag%usef=.true.
 arbre=>f%mag%forward
else
 if(.not.associated(f%mag%backward)) then
  allocate(f%mag%backward(3))
 ! allocate(f%mag%useb)
 else
  call KILL(f%mag%backward)
 endif
 call SET_TREE_G_complex(f%mag%backward,m,fact)
 f%mag%do1mapb=onemap
 f%mag%useb=.true.
 arbre=>f%mag%backward
endif

arbre(1)%rad=mat
arbre(1)%fix0(1:6)=fix0
arbre(1)%fixr(1:6)=fixr
arbre(1)%fix(1:6)=fix
if(onemap) then
 arbre(1)%ds=f%mag%p%ld
else
 arbre(1)%ds=f%mag%p%ld/f%mag%p%nst
endif
arbre(1)%beta0=f%beta0

if(f%dir==1) then
 if(.not.associated(f%magp%forward)) then
  allocate(f%magp%forward(3))
!  allocate(f%magp%usef)
 else
  call KILL(f%magp%forward)
 endif
 !call SET_TREE_G_complex(f%magp%forward,m)
do i=1,3
 call alloc_tree(f%magp%forward(i),f%mag%forward(i)%n,f%mag%forward(i)%np)
 call copy_tree(f%mag%forward(i),f%magp%forward(i))
enddo
 f%magp%do1mapf=onemap
 f%magp%usef=.true.
 arbre=>f%magp%forward
else

 if(.not.associated(f%magp%backward)) then
  allocate(f%magp%backward(3))
 ! allocate(f%magp%useb)
 else
  call KILL(f%magp%backward)
 endif
 !call SET_TREE_G_complex(f%magp%backward,m)
do i=1,3
 call alloc_tree(f%magp%backward(i),f%mag%backward(i)%n,f%mag%backward(i)%np)
 call copy_tree(f%mag%backward(i),f%magp%backward(i))
enddo
 f%magp%do1mapb=onemap
 f%magp%useb=.true.
 arbre=>f%magp%backward
endif

arbre(1)%rad=mat
arbre(1)%fix0(1:6)=fix0
arbre(1)%fixr(1:6)=fixr
arbre(1)%fix(1:6)=fix
if(onemap) then
 arbre(1)%ds=f%mag%p%ld
else
 arbre(1)%ds=f%mag%p%ld/f%mag%p%nst
endif
arbre(1)%beta0=f%beta0

call kill(xs);call kill(m)
 
 if(present(file)) then
  call kanalnummer(mf,file)
   call print_tree_elements(arbre,mf)
  close(mf)
 endif


end subroutine fill_tree_element

subroutine fill_tree_element_line(f1,f2,f,no,fix0,factor,nocav,file)   ! fix0 is the initial condition for the maps
implicit none
type(fibre), target :: f1,f2,f
type(layout), pointer :: r
TYPE(INTEGRATION_NODE),POINTER:: t1c,t2c
TYPE (NODE_LAYOUT), POINTER :: t
TYPE (tree_element), POINTER :: arbre(:)
type(internal_state) state
real(dp) fixr(6),fixs(6),fix(6),fix0(6),mat(6,6),e_ij(6,6),xn
type(probe) xs0
type(probe_8) xs
type(c_damap) m,mr
logical :: fact,noca
logical,optional :: factor,nocav
integer no,i,mf
type(fibre), pointer :: p
character(*), optional :: file 

fact=.false.
noca=.false.

if(present(factor)) fact=factor
if(present(factor)) noca=nocav

if(.not.associated(f1%parent_layout)) then
 write(6,*) " parent layout not associated "
 stop
else
 r=>f1%parent_layout
endif

if(.not.associated(f%parent_layout%t)) then
 write(6,*) " parent node layout not associated "
 stop
else
 t=>f1%parent_layout%t
 t1c=>f1%t1 !%next
 t2c=>f2%t1
endif

! Classical radiation with stochastic envelope

state=radiation0+envelope0+time0
state%NOCAVITY=noca

call init_all(state,1,0)

call alloc(xs);call alloc(m,mr)


! radiation

xs0=fix0
mr=1
xs=xs0+mr
if(associated(t1c,t2c)) then
 call propagate(xs,state,node1=t1c)
else
 call propagate(xs,state,node1=t1c,node2=t2c)
endif

! For David
!!  mr: linear map with radiation would be read here instead of being computed,
!! and must be stored in fixr

fixr=xs%x    ! <---
mr=xs   ! <---

do i=1,6
 mr%v(i)=mr%v(i)-(mr%v(i).sub.0)
enddo

! For David
!!  The stochastic kicks are stored at e_ij

e_ij=xs%e_ij         ! <---
! no radiation


state=time0
state%NOCAVITY=noca
xs0=fix0
m=1

if(associated(t1c,t2c)) then
 call propagate(xs,state,node1=t1c)
else
 call propagate(xs,state,node1=t1c,node2=t2c)
endif

fix=xs%x
! For David
!!  The same linear map is computed WITHOUT radiation : result put into m, the constant part is removed
!!
m=xs   ! <---
do i=1,6
 m%v(i)=m%v(i)-(m%v(i).sub.0)
enddo

m=m**(-1)*mr



mat=m


call kill(xs);call kill(m);call kill(mr)

state=spin0+time0
state%NOCAVITY=noca
call init_all(state,no,0)
call alloc(xs);call alloc(m)



xs0=fix0
m=1
xs=xs0+m
if(associated(t1c,t2c)) then
 call propagate(xs,state,node1=t1c)
else
 call propagate(xs,state,node1=t1c,node2=t2c)
endif


! For David
!!  The full nonlinear map m is computed and the final orbit
!!
fix=xs%x  ! <---
m=xs  ! <---

m%e_ij=e_ij
do i=1,6
 m%v(i)=m%v(i)-(m%v(i).sub.0)
enddo




if(f%dir==1) then
 if(.not.associated(f%mag%forward)) then
  allocate(f%mag%forward(3))
 ! allocate(f%mag%usef)
 else
  call KILL(f%mag%forward)
 endif

call SET_TREE_G_complex(f%mag%forward,m,fact)
 f%mag%do1mapf=.false.
 f%mag%usef=.true.
 arbre=>f%mag%forward

else
 if(.not.associated(f%mag%backward)) then
  allocate(f%mag%backward(3))
 ! allocate(f%mag%useb)
 else
  call KILL(f%mag%backward)
 endif
 call SET_TREE_G_complex(f%mag%backward,m,fact)
 f%mag%do1mapb=.false.
 f%mag%useb=.true.
 arbre=>f%mag%backward
endif

arbre(1)%rad=mat
arbre(1)%fix0(1:6)=fix0
arbre(1)%fixr(1:6)=fixr
arbre(1)%fix(1:6)=fix

 arbre(1)%ds=0.0_dp
 p=>f1
 do while(.not.associated(p,f2))
  arbre(1)%ds=p%mag%p%ld +arbre(1)%ds
  p=>p%next
 enddo

arbre(1)%beta0=f1%beta0




if(f%dir==1) then
 if(.not.associated(f%magp%forward)) then
  allocate(f%magp%forward(3))
!  allocate(f%magp%usef)
 else
  call KILL(f%magp%forward)
 endif
 !call SET_TREE_G_complex(f%magp%forward,m)
do i=1,3
 call alloc_tree(f%magp%forward(i),f%mag%forward(i)%n,f%mag%forward(i)%np)
 call copy_tree(f%mag%forward(i),f%magp%forward(i))
enddo
 f%magp%do1mapf=.false.
 f%magp%usef=.true.
 arbre=>f%magp%forward

else

 if(.not.associated(f%magp%backward)) then
  allocate(f%magp%backward(3))
 ! allocate(f%magp%useb)
 else
  call KILL(f%magp%backward)
 endif
 !call SET_TREE_G_complex(f%magp%backward,m)
do i=1,3
 call alloc_tree(f%magp%backward(i),f%mag%backward(i)%n,f%mag%backward(i)%np)
 call copy_tree(f%mag%backward(i),f%magp%backward(i))
enddo
 f%magp%do1mapb=.false.
 f%magp%useb=.true.
 arbre=>f%magp%backward
endif

arbre(1)%rad=mat
arbre(1)%fix0(1:6)=fix0
arbre(1)%fixr(1:6)=fixr
arbre(1)%fix(1:6)=fix
 arbre(1)%ds=0.0_dp
 p=>f1
 do while(.not.associated(p,f2))
  arbre(1)%ds=p%mag%p%ld +arbre(1)%ds
  p=>p%next
 enddo
arbre(1)%beta0=f1%beta0

call kill(xs);call kill(m)
 

 if(present(file)) then
  call kanalnummer(mf,file)
   call print_tree_elements(arbre,mf)
  close(mf)
 endif


end subroutine fill_tree_element_line

!!!!!!!!!!!!!!!!!!!!   stuff for Zhe  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_tree_element_line_zhe0(state_0,state,f1,f2,no,fix0_0,fix0,filef,stochprec,sagan_tree)   ! fix0 is the initial condition for the maps
implicit none
type(fibre), target :: f1,f2 
type(layout), pointer :: r
TYPE(INTEGRATION_NODE),POINTER:: t1c,t2c
TYPE (NODE_LAYOUT), POINTER :: t
type(internal_state), intent(in):: state,state_0
real(dp) fixr(6),fixs(6),fix(6),fix0(6),fix0_0(6),mat(6,6),xn,stoch,fix_0(6)
real(dp), optional :: stochprec
 
type(probe) xs0,xs0_0
type(probe_8) xs,xs_0
type(c_damap) m,mr,m_0
integer no,i,inf
type(fibre), pointer :: p
type(tree_element), pointer :: forward(:) =>null()
character(*),optional :: filef
type(tree_element),optional, target :: sagan_tree(3)

 

if(present(sagan_tree)) then
 forward=>sagan_tree
else
  allocate(forward(3))
endif
if(.not.associated(f1%parent_layout)) then
 write(6,*) " parent layout not associated "
 stop
else
 r=>f1%parent_layout
endif

 t=>f1%parent_layout%t
 t1c=>f1%t1 !%next
 t2c=>f2%t1


mat=0
do i=1,size(mat,1)
mat(i,i)=1
enddo

 
call init_all(state,no,0)
call alloc(xs);call alloc(m,m_0);call alloc(mr)
 call alloc(xs_0);


xs0=fix0
m=1
xs=xs0+m
if(associated(t1c,t2c)) then
 call propagate(xs,state,node1=t1c)
else
 call propagate(xs,state,node1=t1c,node2=t2c)
endif
 


xs0_0=fix0_0
m=1
xs_0=xs0_0+m
if(associated(t1c,t2c)) then
 call propagate(xs_0,state_0,node1=t1c)
else
 call propagate(xs_0,state_0,node1=t1c,node2=t2c)
endif
 
! For David
!!  The full nonlinear map m is computed and the final orbit
!!  
fix=xs%x  ! <---   
m=xs  ! <---   
 
do i=1,6
 m%v(i)=m%v(i)-(m%v(i).sub.0)
enddo 

fix_0=xs_0%x  ! <---   
m_0=xs_0 ! <---   
 
do i=1,6
 m_0%v(i)=m_0%v(i)-(m_0%v(i).sub.0)
enddo 

 

 
call SET_TREE_G_complex_zhe0(forward,m,m_0)

  stoch=-1.0_dp
if(present(stochprec)) stoch=stochprec

if(stoch>=0) then
  call c_stochastic_kick(m,forward(2)%rad,forward(2)%fix0,stoch)  
endif

 
forward(1)%rad=mat
forward(1)%fix0(1:6)=fix0
forward(1)%fixr(1:6)=fix
forward(1)%fix(1:6)=fix    ! always same fixed point
 
 
forward(3)%fix0(1:6)=fix0_0
forward(3)%fixr(1:6)=fix_0
forward(3)%fix(1:6)=fix_0    ! always same fixed point


 forward(1)%ds=0.0_dp
 p=>f1
 do while(.not.associated(p,f2))
  forward(1)%ds=p%mag%p%ld +forward(1)%ds
  p=>p%next
 enddo
forward(1)%beta0=f1%beta0

 if(present(filef)) then
  call kanalnummer(inf,filef)
    call print_tree_elements(forward,inf)
   close(inf)
  call KILL(forward)
  deallocate(forward)
endif

call kill(xs);call kill(m);call kill(mr)
 
end subroutine fill_tree_element_line_zhe0

subroutine fill_tree_element_line_zhe(state,f1,f2,no,fix0,filef,stochprec,sagan_tree)   ! fix0 is the initial condition for the maps
implicit none
type(fibre), target :: f1,f2
type(layout), pointer :: r
TYPE(INTEGRATION_NODE),POINTER:: t1c,t2c
TYPE (NODE_LAYOUT), POINTER :: t
type(internal_state), intent(in):: state
real(dp) fixr(6),fixs(6),fix(6),fix0(6),mat(6,6),xn,stoch
real(dp), optional :: stochprec
 
type(probe) xs0
type(probe_8) xs
type(c_damap) m,mr
integer no,i,inf
type(fibre), pointer :: p
type(tree_element), pointer :: forward(:) =>null()
character(*),optional :: filef
type(tree_element),optional, target :: sagan_tree(3)

 

if(present(sagan_tree)) then
 forward=>sagan_tree
else
  allocate(forward(3))
endif
if(.not.associated(f1%parent_layout)) then
 write(6,*) " parent layout not associated "
 stop
else
 r=>f1%parent_layout
endif

 t=>f1%parent_layout%t
 t1c=>f1%t1 !%next
 t2c=>f2%t1


mat=0
do i=1,size(mat,1)
mat(i,i)=1
enddo


call init_all(state,no,0)
call alloc(xs);call alloc(m);call alloc(mr)



xs0=fix0
m=1
xs=xs0+m
if(associated(t1c,t2c)) then
 call propagate(xs,state,node1=t1c)
else
 call propagate(xs,state,node1=t1c,node2=t2c)
endif


! For David
!!  The full nonlinear map m is computed and the final orbit
!!
fix=xs%x  ! <---
m=xs  ! <---

do i=1,6
 m%v(i)=m%v(i)-(m%v(i).sub.0)
enddo

 


call SET_TREE_G_complex_zhe(forward,m)

  stoch=-1.0_dp
if(present(stochprec)) stoch=stochprec

if(stoch>=0) then
  call c_stochastic_kick(m,forward(2)%rad,forward(2)%fix0,stoch)
endif


forward(1)%rad=mat
forward(1)%fix0(1:6)=fix0
forward(1)%fixr(1:6)=fix
forward(1)%fix(1:6)=fix    ! always same fixed point

 forward(1)%ds=0.0_dp
 p=>f1
 do while(.not.associated(p,f2))
  forward(1)%ds=p%mag%p%ld +forward(1)%ds
  p=>p%next
 enddo
forward(1)%beta0=f1%beta0

 if(present(filef)) then
  call kanalnummer(inf,filef)
    call print_tree_elements(forward,inf)
   close(inf)
  call KILL(forward)
  deallocate(forward)
endif

call kill(xs);call kill(m);call kill(mr)

end subroutine fill_tree_element_line_zhe

!!!!!!!!!!!!!!!!!!!!   tree tracking  for Zhe  : independent program

  SUBROUTINE SET_TREE_G_complex_zhe(T,Ma)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T(:)
    TYPE(c_damap), INTENT(INOUT) :: Ma
    INTEGER N,NP,i,k,j,kq
 
    real(dp) norm,mat(6,6)
    TYPE(taylor), ALLOCATABLE :: M(:), MG(:)
    TYPE(damap) ms
    integer js(6)
    type(c_damap) L_ns , N_pure_ns , N_s , L_s


 

    call alloc(L_ns , N_pure_ns , N_s , L_s)
    
    call symplectify_for_zhe(ma,L_ns , N_pure_ns, L_s , N_s )

!    np=ma%n+18
    if(ma%n/=6) then
     write(6,*) " you need a 6-d map in SET_TREE_G_complex for PTC "
     stop
    endif
    np=size_tree
! initialized in ptc ini
 !   ind_spin(1,1)=1+ma%n;ind_spin(1,2)=2+ma%n;ind_spin(1,3)=3+ma%n;
 !   ind_spin(2,1)=4+ma%n;ind_spin(2,2)=5+ma%n;ind_spin(2,3)=6+ma%n;
 !   ind_spin(3,1)=7+ma%n;ind_spin(3,2)=8+ma%n;ind_spin(3,3)=9+ma%n;
 !   k1_spin(1)=1;k2_spin(1)=1;
 !   k1_spin(2)=1;k2_spin(2)=2;
 !   k1_spin(3)=1;k2_spin(3)=3;
 !   k1_spin(4)=2;k2_spin(4)=1;
 !   k1_spin(5)=2;k2_spin(5)=2;
 !   k1_spin(6)=2;k2_spin(6)=3;
 !   k1_spin(7)=3;k2_spin(7)=1;
 !   k1_spin(8)=3;k2_spin(8)=2;
 !   k1_spin(9)=3;k2_spin(9)=3;


    ALLOCATE(M(NP))
    CALL ALLOC(M,NP)
    ALLOCATE(Mg(NP))
    CALL ALLOC(mg,NP)
    do i=1,np
     m(i)=0.e0_dp
     mg(i)=0.e0_dp
    enddo

      L_ns = L_ns*N_pure_ns

     do i=1,L_ns%n
      m(i)=L_ns%v(i)   ! orbital part
     enddo

    call c_full_norm_spin(Ma%s,k,norm)


if(use_quaternion) then
    call c_full_norm_quaternion(Ma%q,kq,norm)
    if(kq==-1) then
      do i=0,3
        m(ind_spin(1,1)+i)=ma%q%x(i)
      enddo
    elseif(kq/=-1) then
      m(ind_spin(1,1))=1.0_dp
      do i=ind_spin(1,1)+1,size_tree
        m(i)=0.0_dp
      enddo
    endif
else
    if(k==-1) then
      do i=1,3
      do j=1,3
        m(ind_spin(i,j))=ma%s%s(i,j)
      enddo
      enddo
    else
      do i=1,3
        m(ind_spin(i,i))=1.0e0_dp
      enddo
    endif
endif
      js=0
     js(1)=1;js(3)=1;js(5)=1; ! q_i(q_f,p_i) and p_f(q_f,p_i)
     call alloc(ms)


       ms=n_s



     ms=ms**js
!     do i=1,3
!      mg(i)=ms%v(2*i-1)   !  q_i(q_f,p_i)
!      mg(3+i)=ms%v(2*i)   !  p_f(q_f,p_i)
!     enddo
     do i=1,6
      mg(i)=ms%v(i)
     enddo
     do i=1,3
     do j=1,3
       mg(ind_spin(i,j))=ms%v(2*i-1).d.(2*j-1)  !   Jacobian for Newton search
     enddo
     enddo
     call kill(ms)  

     call SET_TREE_g(T(1),m(1:6))
 !    do i=1,ma%n
 !     m(i)=1.0_dp.cmono.i
 !    enddo 
 !    do i=ma%n+1,6
 !     m(i)=0.0_dp
 !    enddo
     call SET_TREE_g(T(2),m(7:15))
 
 !    call SET_TREE_g(T(2),m(1:size_tree))
     call SET_TREE_g(T(3),mg(1:size_tree))

!T(3)%ng=mul
!     write(6,*) " mul ",mul
      t(3)%rad=L_s


       mat=ma**(-1)
       t(1)%e_ij=ma%e_ij     !matmul(matmul(mat,ma%e_ij),transpose(mat))  not necessary I think

  

    call kill(m); call kill(mg);
    deallocate(M);    deallocate(Mg);
    call kill(L_ns , N_pure_ns , N_s , L_s)

  END SUBROUTINE SET_TREE_G_complex_zhe

  SUBROUTINE SET_TREE_G_complex_zhe0(T,Ma,ma_0)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T(:)
    TYPE(c_damap), INTENT(INOUT) :: Ma,ma_0
    INTEGER N,NP,i,k,j,kq
 
    real(dp) norm,mat(6,6)
    TYPE(taylor), ALLOCATABLE :: M(:), MG(:)
    TYPE(damap) ms
    integer js(6)
    type(c_damap) L_ns , N_pure_ns , N_s , L_s


 

    call alloc(L_ns , N_pure_ns , N_s , L_s)
    
    call symplectify_for_zhe0(ma,ma_0,L_ns , N_pure_ns, L_s , N_s )
    
!    np=ma%n+18
    if(ma%n/=6) then
     write(6,*) " you need a 6-d map in SET_TREE_G_complex for PTC "
     stop
    endif
    np=size_tree
! initialized in ptc ini
 !   ind_spin(1,1)=1+ma%n;ind_spin(1,2)=2+ma%n;ind_spin(1,3)=3+ma%n;
 !   ind_spin(2,1)=4+ma%n;ind_spin(2,2)=5+ma%n;ind_spin(2,3)=6+ma%n;
 !   ind_spin(3,1)=7+ma%n;ind_spin(3,2)=8+ma%n;ind_spin(3,3)=9+ma%n;    
 !   k1_spin(1)=1;k2_spin(1)=1;
 !   k1_spin(2)=1;k2_spin(2)=2;
 !   k1_spin(3)=1;k2_spin(3)=3;
 !   k1_spin(4)=2;k2_spin(4)=1;
 !   k1_spin(5)=2;k2_spin(5)=2;
 !   k1_spin(6)=2;k2_spin(6)=3;
 !   k1_spin(7)=3;k2_spin(7)=1;
 !   k1_spin(8)=3;k2_spin(8)=2;
 !   k1_spin(9)=3;k2_spin(9)=3;

   
    ALLOCATE(M(NP))
    CALL ALLOC(M,NP)
    ALLOCATE(Mg(NP))
    CALL ALLOC(mg,NP)
    do i=1,np
     m(i)=0.e0_dp
     mg(i)=0.e0_dp
    enddo
     
      L_ns = L_ns*N_pure_ns

     do i=1,L_ns%n
      m(i)=L_ns%v(i)   ! orbital part
     enddo

    call c_full_norm_spin(Ma%s,k,norm)


if(use_quaternion) then
    call c_full_norm_quaternion(Ma%q,kq,norm)
    if(kq==-1) then
      do i=0,3
        m(ind_spin(1,1)+i)=ma%q%x(i)
      enddo
    elseif(kq/=-1) then
      m(ind_spin(1,1))=1.0_dp
      do i=ind_spin(1,1)+1,size_tree
        m(i)=0.0_dp
      enddo
    endif
else
    if(k==-1) then
      do i=1,3
      do j=1,3
        m(ind_spin(i,j))=ma%s%s(i,j)
      enddo
      enddo
    else
      do i=1,3
        m(ind_spin(i,i))=1.0e0_dp
      enddo
    endif
endif
      js=0
     js(1)=1;js(3)=1;js(5)=1; ! q_i(q_f,p_i) and p_f(q_f,p_i)
     call alloc(ms)

 
       ms=n_s
 
 

     ms=ms**js
!     do i=1,3
!      mg(i)=ms%v(2*i-1)   !  q_i(q_f,p_i)
!      mg(3+i)=ms%v(2*i)   !  p_f(q_f,p_i)
!     enddo
     do i=1,6
      mg(i)=ms%v(i) 
     enddo
     do i=1,3
     do j=1,3
       mg(ind_spin(i,j))=ms%v(2*i-1).d.(2*j-1)  !   Jacobian for Newton search
     enddo
     enddo
          call kill(ms)  

     call SET_TREE_g(T(1),m(1:6))
 !    do i=1,ma%n
 !     m(i)=1.0_dp.cmono.i
 !    enddo 
 !    do i=ma%n+1,6
 !     m(i)=0.0_dp
 !    enddo
     call SET_TREE_g(T(2),m(7:15))
 
 !    call SET_TREE_g(T(2),m(1:size_tree))
     call SET_TREE_g(T(3),mg(1:size_tree))

!T(3)%ng=mul
!     write(6,*) " mul ",mul
      t(3)%rad=L_s
 

       mat=ma**(-1)
       t(1)%e_ij=ma%e_ij     !matmul(matmul(mat,ma%e_ij),transpose(mat))  not necessary I think

  

    call kill(m); call kill(mg);
    deallocate(M);    deallocate(Mg);
    call kill(L_ns , N_pure_ns , N_s , L_s)

  END SUBROUTINE SET_TREE_G_complex_zhe0



subroutine symplectify_for_zhe(m,L_ns , N_pure_ns , L_s, N_s )
implicit none
TYPE(c_damap),intent(inout):: m ,L_ns , N_pure_ns , N_s , L_s
type(c_vector_field) f,fs
complex(dp) v
type(c_taylor) t,dt
real(dp),allocatable::  mat(:,:)
integer i,j,k,n(11),nv,nd2,al,ii,a,mul
integer, allocatable :: je(:)
real(dp) dm,norm,normb,norma
TYPE(c_damap) mt
real(dp),allocatable::   S(:,:),id(:,:)

! m = L_ns o N_pure_ns o L_s o N_s
! d= = L_ns o N_pure_ns
! ms= L_s o N_s

allocate(S(m%n,m%n),id(m%n,m%n))

call c_get_indices(n,0)
nv=n(4)
nd2=n(3)

S=0
id=0
do i=1,nd2/2
 S(2*I-1,2*I)=1 ; S(2*I,2*I-1)=-1;
 Id(2*I-1,2*I-1)=1 ; id(2*I,2*I)=1;
enddo




 call alloc(f);call alloc(fs);
call alloc(t,dt);call alloc(mt);

allocate(mat(m%n,m%n))

mat=0



if(nv-nd2==0) then
mat=m.sub.1
else
write(6,*) " this map should not have parameters or modulated magnets "
stop 444
endif


! constructing Furman's contracting matrix from my review sec.3.8.2

call furman_symp(mat)
L_s=mat

mt=m*L_s**(-1)
L_ns=mt.sub.1

mt=L_ns**(-1)*mt

f=log(mt)

fs=0

! Integrating a symplectic operator using the hypercube's diagonal

allocate(je(nv))
je=0
do i=1,f%n

       j=1

        do while(.true.)

          call  c_cycle(f%v(i),j,v ,je); if(j==0) exit;
         dm=1
         do ii=1,nd2
          dm=dm+je(ii)
         enddo
        t=v.cmono.je
        do a=1,nd2
         dt=t.d.a
        do al=1,nd2
        do k=1,nd2
          fs%v(al)=fs%v(al)+s(a,al)*s(k,i)*(id(k,a)*t+(1.0_dp.cmono.k)*dt)/dm
        enddo ! k
        enddo ! al
        enddo ! a
        enddo

enddo
 


N_s=exp(fs)
N_pure_ns= mt*N_s**(-1)

N_s= L_s**(-1)*N_s*L_s 

!norma=1.d0/mul
!fs=norma*log(n_s)
!N_s=exp(fs)


deallocate(je);deallocate(s,id);
 call kill(f);call kill(fs);
call kill(t,dt);call kill(mt);
deallocate(mat)
end subroutine symplectify_for_zhe

subroutine symplectify_for_zhe0(m,m0,L_ns , N_pure_ns , L_s, N_s )
implicit none
TYPE(c_damap),intent(inout):: m ,m0,L_ns , N_pure_ns , N_s , L_s
type(c_vector_field) f,fs
complex(dp) v
type(c_taylor) t,dt
real(dp),allocatable::  mat(:,:)
integer i,j,k,n(11),nv,nd2,al,ii,a,mul
integer, allocatable :: je(:)
real(dp) dm,norm,normb,norma
TYPE(c_damap) mt
real(dp),allocatable::   S(:,:),id(:,:)

! m = L_ns o N_pure_ns o L_s o N_s
! d= = L_ns o N_pure_ns
! ms= L_s o N_s

allocate(S(m%n,m%n),id(m%n,m%n))

call c_get_indices(n,0)
nv=n(4)
nd2=n(3)

S=0
id=0
do i=1,nd2/2
 S(2*I-1,2*I)=1 ; S(2*I,2*I-1)=-1;
 Id(2*I-1,2*I-1)=1 ; id(2*I,2*I)=1;
enddo




 call alloc(f);call alloc(fs);
call alloc(t,dt);call alloc(mt);

allocate(mat(m%n,m%n))

mat=0



if(nv-nd2==0) then
mat=m0.sub.1
else
write(6,*) " this map should not have parameters or modulated magnets "
stop 444
endif


! constructing Furman's contracting matrix from my review sec.3.8.2

call furman_symp(mat)
L_s=mat

N_s=m0*L_s**(-1)

!goto 1111
f=log(N_s)
 
fs=0

! Integrating a symplectic operator using the hypercube's diagonal

allocate(je(nv))
je=0
do i=1,f%n

       j=1

        do while(.true.) 

          call  c_cycle(f%v(i),j,v ,je); if(j==0) exit;
         dm=1
         do ii=1,nd2
          dm=dm+je(ii)
         enddo
        t=v.cmono.je
        do a=1,nd2
         dt=t.d.a
        do al=1,nd2
        do k=1,nd2
          fs%v(al)=fs%v(al)+s(a,al)*s(k,i)*(id(k,a)*t+(1.0_dp.cmono.k)*dt)/dm
        enddo ! k
        enddo ! al
        enddo ! a
        enddo

enddo
 


N_s=exp(fs)


!1111 continue 


N_pure_ns= m*m0**(-1)

L_ns=N_pure_ns.sub.1

N_pure_ns=  L_ns**(-1)*N_pure_ns

N_s= L_s**(-1)*N_s*L_s 

 
 

deallocate(je);deallocate(s,id);
 call kill(f);call kill(fs);
call kill(t,dt);call kill(mt);
deallocate(mat)
end subroutine symplectify_for_zhe0



  SUBROUTINE furman_symp(r)
   implicit none
   real(dp)  r(:,:)
    real(dp), allocatable::rt(:,:)
    real(dp) eps,a,ab
   integer nmax,i,j,k,n
! Furmanizing the rotation
    n=size(r,1)
    allocate(rt(n,n))

    rt=0
    eps=1.d-8
    nmax=1000

    ab=1.d8
    do i=1,nmax
    ! rt=matmul(r,transpose(r))
    ! r= matmul((id-0.5e0_dp*rt),r)

      call furman_step(r,r,rt)

     a=-n
     do j=1,n
     do k=1,n
      a=a+abs(rt(j,k))
     enddo
     enddo
     a=abs(a)
     if(a<eps) then
      if(a>=ab) exit
      ab=a
     endif
    enddo
    if(i>nrmax-10) then
     write(6,*) i, a, "did not converge in orthonormalisep"
     ! stop
    endif
    deallocate(rt)
  end SUBROUTINE furman_symp

  SUBROUTINE furman_step(r,s,rt)
   implicit none
   real(dp)  r(:,:),s(:,:),rt(:,:)
   real(dp), allocatable :: id(:,:),ik(:,:), j(:,:),ji(:,:)
   integer i,n

    n=size(r,1)

allocate(id(n,n),ik(n,n), j(n,n),ji(n,n))
    id=0
    ik=0
     j=0
      ji=0
      do i=1,n/2
       j(2*i-1,2*i)=1
       j(2*i,2*i-1)=-1
      enddo
      ji=-j
      do i=1,n
       ik(i,i)=1.5e0_dp
      enddo


       id=matmul(ik-0.5_dp*  matmul( matmul(r,j),matmul(transpose(r),ji) )  , r)




        s=id

      rt=matmul(r,matmul(j,transpose(r)) )
deallocate(id,ik, j,ji)
end   SUBROUTINE furman_step

!!!!!!!!!!!!!!!!!!!!   stuff for Zhe  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! extract moments using initial moment and the probe_8 containing fluctuations and map
  subroutine extract_moments(p8,s_ij_in,s_ij_out)
    implicit none
    real(dp), intent(in) :: s_ij_in(6,6)
    real(dp), intent(out) :: s_ij_out(6,6)
    type(probe_8), intent(inout) ::p8
    type(damap) m
    real(dp) mat(6,6)

    call alloc(m)

    m=p8%x
    mat=m

    s_ij_out=0.0_dp
    s_ij_out=s_ij_in + p8%e_ij
    s_ij_out=matmul(mat,s_ij_out)
    mat=transpose(mat)
    s_ij_out=matmul(s_ij_out,mat)


    call kill(m)
  end subroutine extract_moments




end module ptc_spin

