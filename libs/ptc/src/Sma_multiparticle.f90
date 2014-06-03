
!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module ptc_multiparticle
  !  use S_TRACKING  !,FRINGE_=>FRINGE__MULTI !,FACE=>FACE_MULTI
  use beam_beam_ptc
  implicit none
  public
  CHARACTER*27 CASE_NAME(-2:3)
  PRIVATE fuzzy_eq,fuzzy_neq
  PRIVATE TRACK_FIBRE_FRONTR,TRACK_FIBRE_FRONTP
  PRIVATE TRACK_FIBRE_BACKR,TRACK_FIBRE_BACKP
  PRIVATE TRACKR_NODE_SINGLE,TRACKP_NODE_SINGLE,TRACKV_NODE_SINGLE
  private DRIFTr_BACK_TO_POSITION,DRIFTp_BACK_TO_POSITION  !,DRIFT_BACK_TO_POSITION
  private MAKE_NODE_LAYOUT_2 !,DRIFT_TO_TIME
  PRIVATE MODULATE_R,MODULATE_P
  PRIVATE TRACK_MODULATION_R,TRACK_MODULATION_P

  !  LOGICAL :: OLD_MOD=.TRUE.

  logical(lp),private, parameter :: dobb=.true.
  logical(lp),private, parameter :: aperture_all_case0=.false.
  type(probe) :: xsm,xsm0
  !real(dp) :: unit_time =1.0e-3_dp
  REAL(dp) :: x_orbit_sync(6)= 0.0_dp,dt_orbit_sync=0.0_dp
  
  INTERFACE TRACK_NODE_SINGLE
     MODULE PROCEDURE TRACKR_NODE_SINGLE     !@1  t,x,state,charge
     MODULE PROCEDURE TRACKP_NODE_SINGLE     !@1  t,y,state,charge
     MODULE PROCEDURE TRACKV_NODE_SINGLE     !@1  t,v,state,charge
  END INTERFACE




  INTERFACE DRIFT_BACK_TO_POSITION
     MODULE PROCEDURE DRIFTr_BACK_TO_POSITION
     MODULE PROCEDURE DRIFTp_BACK_TO_POSITION
  END INTERFACE

  INTERFACE TRACK_FIBRE_FRONT
     MODULE PROCEDURE TRACK_FIBRE_FRONTR
     MODULE PROCEDURE TRACK_FIBRE_FRONTP
  END INTERFACE

  INTERFACE TRACK_FIBRE_BACK
     MODULE PROCEDURE TRACK_FIBRE_BACKR
     MODULE PROCEDURE TRACK_FIBRE_BACKP
  END INTERFACE

  INTERFACE COPY
     MODULE PROCEDURE COPY_BEAM
  END INTERFACE

  INTERFACE OPERATOR (.feq.)
     MODULE PROCEDURE fuzzy_eq
  END INTERFACE
  INTERFACE OPERATOR (.fne.)
     MODULE PROCEDURE fuzzy_neq
  END INTERFACE

  INTERFACE MODULATE
     MODULE PROCEDURE MODULATE_R   ! PROPAGATE FAKE MODULATED (Q,P) ; 7TH AND 8TH VARIABLES
     MODULE PROCEDURE MODULATE_P   !
  END INTERFACE


  INTERFACE TRACK_MODULATION
     MODULE PROCEDURE TRACK_MODULATION_R   ! PROPAGATE FAKE MODULATED (Q,P) ; 7TH AND 8TH VARIABLES
     MODULE PROCEDURE TRACK_MODULATION_P   !
  END INTERFACE

  type three_d_info
     !   character(nlp),pointer :: name
     real(dp)  a(3),b(3)   ! Centre of entrance and exit faces
     real(dp)  ent(3,3),exi(3,3)  ! entrace and exit frames for drawing magnet faces
     real(dp)  wx,wy ! width of box for plotting purposes
     real(dp)  o(3),mid(3,3)   ! frames at the point of tracking
     real(dp)  reference_ray(6)  !
     real(dp) x(6)   ! ray tracked with reference_ray using a  type(beam)
     real(dp) r0(3),r(3)  ! ray position global returned
     real(dp) scale     !  magnification using reference_ray
     logical(lp) u(2)   ! unstable flag for both ray and reference_ray
  END type three_d_info

!  real :: ttime0,ttime1,dt1=0.0,dt2=0.0;

CONTAINS

  SUBROUTINE MODULATE_R(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    TYPE(ELEMENT),POINTER :: EL
    TYPE(ELEMENTP),POINTER :: ELp
    REAL(DP) v,dv



    EL=>C%PARENT_FIBRE%MAG
    ELP=>C%PARENT_FIBRE%MAGP

    IF(K%MODULATION) THEN

       DV=(XS%AC%X(1)*COS(EL%theta_ac)-XS%AC%X(2)*SIN(EL%theta_ac))
       V=EL%DC_ac+EL%A_ac*DV
       DV=el%D_ac*DV
     else
       V=0.0_dp
       DV=0.0_dp
     endif

    CALL transfer_ANBN(EL,ELP,VR=V,DVR=DV)


  END   SUBROUTINE MODULATE_R
  
   SUBROUTINE do_ramping_R(C,t,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    TYPE(INTERNAL_STATE) K
    TYPE(ELEMENT),POINTER :: EL
    TYPE(ELEMENTP),POINTER :: ELp
    REAL(DP) v,dv,t



    EL=>C%PARENT_FIBRE%MAG
    ELP=>C%PARENT_FIBRE%MAGP
    if(.not.associated(EL%ramp)) return

       V=EL%DC_ac
       DV=0.0_dp
       call set_ramp(C,t)

    CALL transfer_ANBN(EL,ELP,VR=V,DVR=DV)


  END   SUBROUTINE do_ramping_R
  
   SUBROUTINE DO_Ramping_p(C,t,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    TYPE(INTERNAL_STATE) K
    TYPE(ELEMENT),POINTER :: EL
    TYPE(ELEMENTP),POINTER :: ELP
    TYPE(REAL_8) V,DV
    real(dp) t

    EL=>C%PARENT_FIBRE%MAG
    ELP=>C%PARENT_FIBRE%MAGP
    
    if(.not.associated(EL%ramp)) return

    CALL ALLOC(V)
    CALL ALLOC(DV)


       V=elp%DC_ac
       DV=0.0_dp
         call set_ramp(C,t)

    CALL transfer_ANBN(EL,ELP,VP=V,DVP=DV)

    CALL KILL(V)
    CALL KILL(DV)


  END   SUBROUTINE DO_Ramping_p

   SUBROUTINE set_all_ramp(R)
    IMPLICIT NONE
    TYPE(layout), target  :: r
    TYPE(fibre), POINTER  :: p
    integer i
    REAL(DP) v,dv
    v=0.0_dp
    dv=0.0_dp     
     p=>r%start
     do i=1,r%n
     
      if(associated(p%mag%ramp)) then
       call set_ramp(p%t1,x_orbit_sync(6))
       CALL transfer_ANBN(p%mag,p%magp,VR=V,DVR=DV)

      endif 
     
      p=>p%next
     enddo      
     
    end SUBROUTINE set_all_ramp
     
  SUBROUTINE set_ramp(t,t0)
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE), POINTER  :: T
    integer i,it
    real(dp) r,ti,rat,dtot
    type(ramping), pointer :: a
    real(dp) an,bn,t0
    

 !   if(t%pos_in_fibre==1) return
    
    a=>t%parent_fibre%mag%ramp
    

    dtot=(a%table(a%n)%time-a%table(1)%time)  !/(a%n-1)
    
 !   ti=XSM0%ac%t/clight/a%unit_time    ! time in milliseconds
    ti=t0/clight    !/a%unit_time    ! time in milliseconds
    
!    if(ti>a%t_max.or.ti<a%table(1)%time) then
!    if(ti>a%table(a%n)%time.or.ti<a%table(1)%time) then
!     return
!    endif
     

     
    if(ti>=a%t_max.or.ti<a%table(1)%time) then
!    if(ti>a%table(a%n)%time.or.ti<a%table(1)%time) then
      if(ti>=a%t_max) then
           a%table(0)%bn=0.0_dp
           a%table(0)%an=0.0_dp
          do i=1,size(a%table(0)%bn)
           a%table(0)%bn(i)= a%table(a%n)%bn(i)*a%r
           a%table(0)%an(i)= a%table(a%n)%an(i)*a%r
          enddo 
            a%table(0)%b_t= a%table(a%n)%b_t
           a=>t%parent_fibre%magp%ramp
           a%table(0)%bn=0.0_dp
           a%table(0)%an=0.0_dp
          do i=1,size(a%table(0)%bn)
           a%table(0)%bn(i)= a%table(a%n)%bn(i)*a%r
           a%table(0)%an(i)= a%table(a%n)%an(i)*a%r
          enddo 
          a%table(0)%b_t= a%table(a%n)%b_t
      else
           a%table(0)%bn=0.0_dp
           a%table(0)%an=0.0_dp
          do i=1,size(a%table(0)%bn)
           a%table(0)%bn(i)= a%table(1)%bn(i)*a%r
           a%table(0)%an(i)= a%table(1)%an(i)*a%r
          enddo 
         a%table(0)%b_t= a%table(1)%b_t
          a=>t%parent_fibre%magp%ramp
           a%table(0)%bn=0.0_dp
           a%table(0)%an=0.0_dp
          do i=1,size(a%table(0)%bn)
           a%table(0)%bn(i)= a%table(1)%bn(i)*a%r
           a%table(0)%an(i)= a%table(1)%an(i)*a%r
          enddo 
          a%table(0)%b_t= a%table(1)%b_t
     
      endif

    else
    
         ti=ti-a%table(1)%time
         ti=mod(ti,dtot)+a%table(1)%time
          dtot=dtot/(a%n-1)
          ti=(ti-a%table(1)%time)/dtot+1
           
          it=int(ti)
!          it=idint(ti)
          
           rat=(ti-it)    
           

           a%table(0)%bn=0.0_dp
           a%table(0)%an=0.0_dp
          do i=1,size(a%table(0)%bn)
           a%table(0)%bn(i)=((a%table(it+1)%bn(i)-a%table(it)%bn(i))*rat + a%table(it)%bn(i))*a%r
           a%table(0)%an(i)= ((a%table(it+1)%an(i)-a%table(it)%an(i))*rat + a%table(it)%an(i))*a%r
          enddo 
           a%table(0)%b_t=((a%table(it+1)%b_t-a%table(it)%b_t)*rat + a%table(it)%b_t)

          a=>t%parent_fibre%magp%ramp
           a%table(0)%bn=0.0_dp
           a%table(0)%an=0.0_dp
          do i=1,size(a%table(0)%bn)
           a%table(0)%bn(i)=((a%table(it+1)%bn(i)-a%table(it)%bn(i))*rat + a%table(it)%bn(i))*a%r
           a%table(0)%an(i)= ((a%table(it+1)%an(i)-a%table(it)%an(i))*rat + a%table(it)%an(i))*a%r
          enddo 
           a%table(0)%b_t=((a%table(it+1)%b_t-a%table(it)%b_t)*rat + a%table(it)%b_t)
          
    endif
  end SUBROUTINE set_ramp
  
  SUBROUTINE MODULATE_P(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe_8), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    TYPE(ELEMENT),POINTER :: EL
    TYPE(ELEMENTP),POINTER :: ELP
    TYPE(REAL_8) V,DV

    EL=>C%PARENT_FIBRE%MAG
    ELP=>C%PARENT_FIBRE%MAGP

    CALL ALLOC(V)
    CALL ALLOC(DV)

    IF(K%MODULATION) THEN
       DV=(XS%AC%X(1)*COS(ELP%theta_ac)-XS%AC%X(2)*SIN(ELP%theta_ac))
       V=ELP%DC_ac+ELP%A_ac*DV
       DV=elp%D_ac*DV

       else  ! ramp
          V=0.0_dp
          DV=0.0_dp
       endif
 
    CALL transfer_ANBN(EL,ELP,VP=V,DVP=DV)

    CALL KILL(V)
    CALL KILL(DV)


  END   SUBROUTINE MODULATE_P

  
  SUBROUTINE TRACK_MODULATION_R(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    real(dp) xt
    real(dp),pointer :: beta0

    if(k%time) then
       beta0=>C%PARENT_FIBRE%beta0
       xs%ac%t=c%DS_AC/beta0+xs%ac%t
       xt = cos(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(1) + sin(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(2)
       XS%AC%X(2) = -sin(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(1) + cos(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(2)
       XS%AC%X(1) = xt
    else
       xt = cos(XS%AC%om * c%DS_AC) *XS%AC%X(1) + sin(XS%AC%om * c%DS_AC) *XS%AC%X(2)
       XS%AC%X(2) = -sin(XS%AC%om * c%DS_AC) *XS%AC%X(1) + cos(XS%AC%om * c%DS_AC) *XS%AC%X(2)
       XS%AC%X(1) = xt
       xs%ac%t=c%DS_AC+xs%ac%t
    endif

  END   SUBROUTINE TRACK_MODULATION_R

  SUBROUTINE TRACK_MODULATION_P(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe_8), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    TYPE(REAL_8) xt
    real(dp),pointer :: beta0

    CALL ALLOC(XT)

    if(k%time) then
       beta0=>C%PARENT_FIBRE%beta0
       xs%ac%t=c%DS_AC/beta0+xs%ac%t
       xt = cos(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(1) + sin(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(2)
       XS%AC%X(2) = -sin(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(1) + cos(XS%AC%om * c%DS_AC/beta0) *XS%AC%X(2)
       XS%AC%X(1) = xt
    else
       xt = cos(XS%AC%om * c%DS_AC) *XS%AC%X(1) + sin(XS%AC%om * c%DS_AC) *XS%AC%X(2)
       XS%AC%X(2) = -sin(XS%AC%om * c%DS_AC) *XS%AC%X(1) + cos(XS%AC%om * c%DS_AC) *XS%AC%X(2)
       XS%AC%X(1) = xt
       xs%ac%t=c%DS_AC+xs%ac%t
    endif

    CALL KILL(XT)

  END   SUBROUTINE TRACK_MODULATION_P



  FUNCTION fuzzy_eq( S1, S2 )
    implicit none
    logical(lp) fuzzy_eq
    real(dp), INTENT (IN) :: S1,S2
    fuzzy_eq=.false.

    if(abs(s1-s2)<=c_%eps_pos) fuzzy_eq=.true.

  end FUNCTION fuzzy_eq

  FUNCTION fuzzy_neq( S1, S2 )
    implicit none
    logical(lp) fuzzy_neq
    real(dp), INTENT (IN) :: S1,S2
    fuzzy_neq=.false.

    if(abs(s1-s2)>c_%eps_pos) fuzzy_neq=.true.

  end FUNCTION fuzzy_neq

  SUBROUTINE move_to_s( L,s,current,i,ds ) ! Moves position s
    implicit none
    TYPE (INTEGRATION_NODE), POINTER :: Current
    TYPE (NODE_LAYOUT) L
    real(dp) s,sp,ds
    integer i,k
    logical(lp) DOIT !,track_it

    !  track_it=.false.
    sp=mod(s,L%END%S(3))

    if(sp==0.0_dp.and.s/=0.0_dp) then
       current=>l%end
       i=l%n+1
       ds=0.0_dp
       !  track_it=.true.
       return
    endif

    if(sp==0.0_dp) then
       current=>l%start
       i=1
       ds=0.0_dp
       return
    endif


    nullify(current);
    Current => L%LAST

    k=L%LASTPOS
    I=K
    ds=0.0_dp

    IF(SP>CURRENT%S(3) ) then

       do i=k,l%n-1
          if(current%next%s(3)>=sp) exit
          current=>current%next
       enddo

       if(current%next%s(3)/=sp) ds=sp-current%s(3)
       if(current%next%s(3)==sp) THEN
          ds=0.0_dp
          CURRENT=>current%next
          I=I+1
       ENDIF

    elseif(SP<CURRENT%S(3)) then

       do i=k-1,1,-1
          current=>current%previous
          if(current%s(3)<=sp) exit
       enddo

       ds=sp-current%s(3)

    endif

    L%LASTPOS=I; L%LAST => Current;

    if(ds>0.0_dp) then
       if(CURRENT%S(4)-ds.feq.0.0_dp) then
          ds=0.0_dp
          current=>Current%next
          i=i+1
          L%LAST => Current;
       ELSEIF(ds.feq.0.0_dp) THEN
          DS=0.0_dp
       ENDIF
    endif


    DOIT=.TRUE.
    !DOIT=.FALSE.

    if(iabs(CURRENT%cas)==0.OR.iabs(CURRENT%cas)==1) then
       do while(DS==0.0_dp.AND.DOIT)    ! PUTS AT BEGINNING IF DS=ZERO
          CURRENT=>CURRENT%PREVIOUS
          IF(ASSOCIATED(CURRENT)) THEN

             IF((SP.FNE.CURRENT%S(3)).or.CURRENT%cas==-2) THEN
                CURRENT=>CURRENT%NEXT
                DOIT=.FALSE.
             ELSE
                I=I-1
             ENDIF

          ELSE
             CURRENT=>CURRENT%NEXT
             DOIT=.FALSE.
          ENDIF

       enddo
    elseif(iabs(CURRENT%cas)==2) then
       do while(DS==0.0_dp.AND.DOIT)    ! PUTS AT BEGINNING IF DS=ZERO
          CURRENT=>CURRENT%next
          IF(ASSOCIATED(CURRENT)) THEN

             IF((SP.FNE.CURRENT%S(3)).or.CURRENT%cas==1) THEN
                CURRENT=>CURRENT%previous
                DOIT=.FALSE.
             ELSE
                I=I+1
             ENDIF

          ELSE
             CURRENT=>CURRENT%previous
             DOIT=.FALSE.
          ENDIF

       enddo
    endif



    L%LASTPOS=I; L%LAST => Current;

    IF(I/=L%LAST%POS) THEN
       WRITE(6,*) " ERROR IN move_to_s ",I,L%LAST%POS
       STOP 999
    ENDIF

  END SUBROUTINE move_to_s

  ! tracking one steps in the body



  ! MULTIPARTICLE AT THE FIBRE LEVEL

  ! front patch/misaglinments/tilt
  SUBROUTINE TRACK_FIBRE_FRONTR(C,X,K)
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    !    TYPE(BEAM),TARGET,INTENT(INOUT):: B
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0

    !FRONTAL PATCH
    !    IF(ASSOCIATED(C%PATCH)) THEN
    PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    !    ELSE
    !       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    !    ENDIF

    ! PUSHING BEAM
    !



    IF(PATCHE/=0.AND.PATCHE/=2) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAG%P%P0C
             B0=>CN%MAG%P%BETA0

             X(2)=X(2)*P0/C%MAG%P%P0C
             X(4)=X(4)*P0/C%MAG%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=root(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/C%MAG%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAG%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAG%P%P0C-1.0_dp
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG==1.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%a_T/c%beta0
      else
       X(6)=X(6)-C%PATCH%a_T
      endif
    ENDIF

    CALL DTILTD(C%DIR,C%MAG%P%TILTD,1,X)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       ou = ALWAYS_EXACTMIS  !K%EXACTMIS.or.
       CALL MIS_FIB(C,X,k,OU,DONEITT)
    ENDIF

  END SUBROUTINE TRACK_FIBRE_FRONTR

  SUBROUTINE TRACK_FIBRE_FRONTP(C,X,K)
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    !    TYPE(BEAM),TARGET,INTENT(INOUT):: B
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0


    !FRONTAL PATCH
    !    IF(ASSOCIATED(C%PATCH)) THEN
    PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    !    ELSE
    !       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    !    ENDIF

    ! PUSHING BEAM
    !



    IF(PATCHE/=0.AND.PATCHE/=2) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN)) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAGP%P%P0C
             B0=>CN%MAGP%P%BETA0

             X(2)=X(2)*P0/C%MAGP%P%P0C
             X(4)=X(4)*P0/C%MAGP%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=SQRT(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAGP%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(SQRT(1.0_dp/C%MAGP%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAGP%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAGP%P%P0C-1.0_dp
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG==1.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%a_T/c%beta0
      else
       X(6)=X(6)-C%PATCH%a_T
      endif
    ENDIF

    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,1,X)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAGP%MIS) THEN
       ou = ALWAYS_EXACTMIS   !K%EXACTMIS.or.
       CALL MIS_FIB(C,X,k,OU,DONEITT)
    ENDIF


  END SUBROUTINE TRACK_FIBRE_FRONTP


  ! back patch/misaglinments/tilt

  SUBROUTINE TRACK_FIBRE_BACKR(C,X,K)
    implicit none
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    !    TYPE(BEAM),TARGET,INTENT(INOUT):: B
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0


    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAG%MIS) THEN
       ou = ALWAYS_EXACTMIS  !K%EXACTMIS.or.
       CALL MIS_FIB(C,X,k,OU,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%DIR,C%MAG%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%b_T/c%beta0
      else
       X(6)=X(6)-C%PATCH%b_T
      endif
    ENDIF

    IF(PATCHG==2.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely

    IF(PATCHE/=0.AND.PATCHE/=1) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       !       P0=>CN%MAG%P%P0C
       !       B0=>CN%MAG%P%BETA0
       P0=>CN%MAG%P%P0C
       B0=>CN%BETA0
       X(2)=X(2)*C%MAG%P%P0C/P0
       X(4)=X(4)*C%MAG%P%P0C/P0
       IF(k%TIME.or.recirculator_cheat)THEN
          X(5)=root(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW
          X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp
       ENDIF
    ENDIF




  END SUBROUTINE TRACK_FIBRE_BACKR

  SUBROUTINE TRACK_FIBRE_BACKP(C,X,K)
    implicit none
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    !    TYPE(BEAM),TARGET,INTENT(INOUT):: B
    type(real_8), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0


    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAGP%MIS) THEN
       ou = ALWAYS_EXACTMIS   !K%EXACTMIS.or.
       CALL MIS_FIB(C,X,k,OU,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%b_T/c%beta0
      else
       X(6)=X(6)-C%PATCH%b_T
      endif
    ENDIF

    IF(PATCHG==2.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely

    IF(PATCHE/=0.AND.PATCHE/=1) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       !       P0=>CN%MAGP%P%P0C
       !       B0=>CN%MAGP%P%BETA0
       P0=>CN%MAGP%P%P0C
       B0=>CN%BETA0
       X(2)=X(2)*C%MAGP%P%P0C/P0
       X(4)=X(4)*C%MAGP%P%P0C/P0
       IF(k%TIME.or.recirculator_cheat)THEN
          X(5)=sqrt(1.0_dp+2.0_dp*X(5)/C%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAGP%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW
          X(5)=(2.0_dp*X(5)+X(5)**2)/(sqrt(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAGP%P%P0C/P0-1.0_dp
       ENDIF
    ENDIF

  END SUBROUTINE TRACK_FIBRE_BACKP


  ! thin lens tracking

  SUBROUTINE TRACKV_NODE_SINGLE(T,V,K) !!
    implicit none
    TYPE(INTEGRATION_NODE),POINTER :: T
    TYPE(INTERNAL_STATE)  K
    REAL(DP) SC,reference_ray(6),x(6)
    type(three_d_info),intent(INOUT) ::  v
    TYPE(INTEGRATION_NODE),POINTER:: mag_in,mag_out

    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif

    if(abs(x(1))+abs(x(3))>absolute_aperture) then
       messageLOST="Sma_multiparticle.f90 TRACKV_NODE_SINGLE : exceed absolute_aperture in TRACKV_NODE_SINGLE"
       lost_node=>t
       lost_fibre=>t%parent_fibre
       xlost=x
       CHECK_STABLE=.false.
      
    endif

    x=V%X
    reference_ray=V%reference_ray

    CALL TRACK_NODE_SINGLE(T,V%X,K)

    IF(.NOT.CHECK_STABLE)      V%U(1)=.TRUE.

    CALL TRACK_NODE_SINGLE(T,V%reference_ray,K)

    IF(.NOT.CHECK_STABLE)   V%U(2)=.TRUE.
    IF(V%U(1).OR.V%U(2)) then
       write(6,*) " Unstable in TRACKV_NODE_SINGLE ", V%U
       RETURN
    endif
    IF(.NOT.ASSOCIATED(T%B)) THEN
       WRITE(6,*) " NO FRAMES IN INTEGRATION NODES "
       STOP 101
    ENDIF

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


  END SUBROUTINE TRACKV_NODE_SINGLE

  SUBROUTINE TRACKR_NODE_SINGLE(T,X,K) !!
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT):: T
    REAL(DP),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el

    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif

    if(abs(x(1))+abs(x(3))>absolute_aperture) then   !.or.(.not.CHECK_MADX_APERTURE)) then
       messageLOST="Sma_multiparticle.f90 TRACKR_NODE_SINGLE : exceed absolute_aperture in TRACKR_NODE_SINGLE"
       lost_node=>t
       lost_fibre=>t%parent_fibre
       xlost=x
       CHECK_STABLE=.false.
       
    endif


    !  T%PARENT_FIBRE%MAG=K
    T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR
    T%PARENT_FIBRE%MAG%P%beta0=>T%PARENT_FIBRE%beta0
    T%PARENT_FIBRE%MAG%P%GAMMA0I=>T%PARENT_FIBRE%GAMMA0I
    T%PARENT_FIBRE%MAG%P%GAMBET=>T%PARENT_FIBRE%GAMBET
    T%PARENT_FIBRE%MAG%P%MASS=>T%PARENT_FIBRE%MASS
    T%PARENT_FIBRE%MAG%P%CHARGE=>T%PARENT_FIBRE%CHARGE
       el=>T%PARENT_FIBRE%MAG

    !call cpu_time(ttime1)

    !dt1=ttime1-ttime0+dt1


    SELECT CASE(T%CAS)
    CASE(CASEP1)
       CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,X,K)
       if(associated(T%PARENT_FIBRE%MAG%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
    CASE(CASEP2)
       CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)

    CASE(CASE1,CASE2)
  !     el=>T%PARENT_FIBRE%MAG
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE.and.t%cas==case2) &
            call check_S_APERTURE_out(el%p,t%POS_IN_FIBRE-2,x)

       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19,kind22)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=X,k=k,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,X,k,1)
             CALL FRINGECAV(EL%C4,X,k=k,J=1)
          ELSE
             CALL FRINGECAV(EL%C4,X,k=k,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,X,k,2)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=X,k=k,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=X,k=k,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=X,k=k,J=T%CAS)
       case(KIND10)
          CALL FRINGE_teapot(EL%TP10,X,k=k,j=T%CAS)
       case(KIND16,KIND20)
          CALL fringe_STREX(EL%K16,X,k,T%CAS)
       case(KIND17)
          STOP 317
       case(KIND21)
          CALL FRINGE_CAV_TRAV(EL%CAV21,X=X,k=k,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,X,k,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)
          CALL ADJUST_WI(EL%WI,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,k,T%CAS)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT

    CASE(CASE0)
 !      el=>T%PARENT_FIBRE%MAG
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE)  &
            call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2,x)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif

       SELECT CASE(EL%KIND)
       CASE(KIND0)
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X,K)
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X,K,t%POS_IN_FIBRE-2)
       case(KIND3)
          CALL TRACK(EL%K3,X,K)
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X,K,t%POS_IN_FIBRE-2)
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X,K)
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X,K)
       case(KIND7)
          CALL TRACK_SLICE(EL%T7,X,K,t%POS_IN_FIBRE-2)
       case(KIND8)
          CALL TRACK(EL%S8,X,K)
       case(KIND9)
          CALL TRACK(EL%S9,X,K)
       case(KIND10)
          CALL TRACK_SLICE(EL%TP10,X,K,t%POS_IN_FIBRE-2)
       case(KIND11:KIND14)
          CALL MONTI(EL%MON14,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%MON14,X,K)
       case(KIND15)
          call SEPTTRACK(EL%SEP15,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%SEP15,X,K)
       case(KIND16,KIND20)
          CALL TRACK_SLICE(EL%K16,X,K,t%POS_IN_FIBRE-2)
       case(KIND17)
          STOP 317
       case(KIND18)
          call RCOLLIMATORI(EL%RCOL18,X,k,t%POS_IN_FIBRE-2)
       case(KIND19)
          CALL ECOLLIMATORI(EL%ECOL19,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%ECOL19,X,K)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,k,t%POS_IN_FIBRE-2)
       case(KIND22)
          CALL TRACK_SLICE(EL%he22,X,k,t%POS_IN_FIBRE-2)
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,k,t%POS_IN_FIBRE-2)
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,k,T%POS_IN_FIBRE-2)

       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT
       if(associated(T%PARENT_FIBRE%MAG%p%aperture).and.aperture_all_case0) &
            call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)


    case(CASET)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)
       endif
       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)
    case(CASETF1,CASETF2)

       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)


    END SELECT
    ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S


    !    T%PARENT_FIBRE%MAG=DEFAULT
    if(wherelost==2.and.(.not.check_stable)) then
       t%lost=t%lost+1
    endif
  END SUBROUTINE TRACKR_NODE_SINGLE


  SUBROUTINE TRACKP_NODE_SINGLE(T,X,K) !!
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT):: T
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(elementp),pointer :: el
    logical(lp) BN2,L
    logical(lp) CHECK_KNOB
    integer(2), pointer,dimension(:)::AN,BN

    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif

    if(abs(x(1))+abs(x(3))>absolute_aperture) then
       messageLOST="Sma_multiparticle.f90 TRACKP_NODE_SINGLE : exceed absolute_aperture in TRACKP_NODE_SINGLE"
       lost_node=>t
       lost_fibre=>t%parent_fibre
       xlost=x
       CHECK_STABLE=.false.
       
    endif

    !   T%PARENT_FIBRE%MAGP=K
    IF(K%PARA_IN ) KNOB=.TRUE.

    T%PARENT_FIBRE%MAGP%P%DIR=>T%PARENT_FIBRE%DIR
    T%PARENT_FIBRE%MAGP%P%beta0=>T%PARENT_FIBRE%beta0
    T%PARENT_FIBRE%MAGP%P%GAMMA0I=>T%PARENT_FIBRE%GAMMA0I
    T%PARENT_FIBRE%MAGP%P%GAMBET=>T%PARENT_FIBRE%GAMBET
    T%PARENT_FIBRE%MAGP%P%MASS=>T%PARENT_FIBRE%MASS
    T%PARENT_FIBRE%MAGP%P%CHARGE=>T%PARENT_FIBRE%CHARGE
       el=>T%PARENT_FIBRE%MAGP


    SELECT CASE(T%CAS)
    CASE(CASEP1)
       CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,X,K)
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,X)
    CASE(CASEP2)
       !    if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then ! new 2010
       !       CHECK_STABLE=.false.
       !    endif
       CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)

    CASE(CASE1,CASE2)
!       el=>T%PARENT_FIBRE%MAGP
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE.and.t%cas==case2) &
            call check_S_APERTURE_out(el%p,t%POS_IN_FIBRE-2,x)


       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19,kind22)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=X,k=k,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,X,k,1)
             CALL FRINGECAV(EL%C4,X,k=k,J=1)
          ELSE
             CALL FRINGECAV(EL%C4,X,k=k,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,X,k,2)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=X,k=k,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=X,k=k,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=X,k=k,J=T%CAS)
       case(KIND10)
          CALL FRINGE_teapot(EL%TP10,X,k,T%CAS)
       case(KIND16,KIND20)
          CALL fringe_STREX(EL%K16,X,k,T%CAS)
       case(KIND17)
          STOP 317
       case(KIND21)
          CALL FRINGE_CAV_TRAV(EL%CAV21,X=X,k=k,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,X,k,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)
          CALL ADJUST_WI(EL%WI,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,k,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT

    CASE(CASE0)
 !      el=>T%PARENT_FIBRE%MAGP
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE) &
            call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2,x)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
       SELECT CASE(EL%KIND)
       CASE(KIND0)
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X,K)
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X,K,t%POS_IN_FIBRE-2)
       case(KIND3)
          CALL TRACK(EL%K3,X,K)
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X,K,t%POS_IN_FIBRE-2)
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X,K)
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X,K)
       case(KIND7)
          IF((EL%T7%BN(2)%KIND==3.OR.EL%T7%L%KIND==3).AND.KNOB) THEN
             CALL GETMAT7(EL%T7)                                      ! RECOMPUTES ONLY IF KNOB (SPEED)
          ENDIF
          CALL TRACK_SLICE(EL%T7,X,K,t%POS_IN_FIBRE-2)
          IF(KNOB) THEN
             BN2=.FALSE.
             L=.FALSE.
             IF(EL%T7%BN(2)%KIND==3) THEN
                BN2=.TRUE.
             ENDIF
             IF(EL%T7%L%KIND==3) THEN
                L=.TRUE.
             ENDIF
             IF(BN2.OR.L) THEN
                EL%T7%BN(2)%KIND=1
                EL%T7%L%KIND=1
                CALL KILL(EL%T7)                               ! RECOMPUTES ONLY IF KNOB (SPEED)
                CALL ALLOC(EL%T7)                               ! KNOB IS REMOVED THE SLOW WAY(SPEED)
                CALL GETMAT7(EL%T7)
                IF(BN2) EL%T7%BN(2)%KIND=3
                IF(L)  EL%T7%L%KIND=3
             ENDIF
          ENDIF
       case(KIND8)
          CALL TRACK(EL%S8,X,K)
       case(KIND9)
          CALL TRACK(EL%S9,X,K)
       case(KIND10)
          CALL MAKEPOTKNOB(EL%TP10,CHECK_KNOB,AN,BN)
          CALL TRACK_SLICE(EL%TP10,X,K,t%POS_IN_FIBRE-2)
          CALL UNMAKEPOTKNOB(EL%TP10,CHECK_KNOB,AN,BN)
       case(KIND11:KIND14)
          CALL MONTI(EL%MON14,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%MON14,X,K)
       case(KIND15)
          call SEPTTRACK(EL%SEP15,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%SEP15,X,K)
       case(KIND16,KIND20)
          CALL TRACK_SLICE(EL%K16,X,K,t%POS_IN_FIBRE-2)
       case(KIND17)
          STOP 317
       case(KIND18)
          call RCOLLIMATORI(EL%RCOL18,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%RCOL18,X,K)
       case(KIND19)
          CALL ECOLLIMATORI(EL%ECOL19,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%ECOL19,X,K)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,k,t%POS_IN_FIBRE-2)
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,k,t%POS_IN_FIBRE-2)
       case(KIND22)
          CALL TRACK_SLICE(EL%he22,X,k,t%POS_IN_FIBRE-2)
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,k,T%POS_IN_FIBRE-2)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT
       if(associated(T%PARENT_FIBRE%MAG%p%aperture).and.aperture_all_case0) &
            call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)

    case(CASET)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)
       endif
       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)
    case(CASETF1,CASETF2)

       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)



    END SELECT
    ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S
    !    T%PARENT_FIBRE%MAGP=DEFAULT
    ! KNOB IS RETURNED TO THE PTC DEFAULT
    ! NEW STUFF WITH KIND=3
    KNOB=.FALSE.
    ! END NEW STUFF WITH KIND=3

  END SUBROUTINE TRACKP_NODE_SINGLE








  !  STUFF ABOUT LIST AND STRUCTURES

  SUBROUTINE MAKE_NODE_LAYOUT( R) !
    ! Creates the thin layout and puts in R%T by calling MAKE_NODE_LAYOUT_2
    ! At this point large patches would be problematic; i.e. |d(3)|>>>0 .
    implicit none
    TYPE (LAYOUT), TARGET :: R

    call MAKE_NODE_LAYOUT_2( R,R%T )

  end SUBROUTINE MAKE_NODE_LAYOUT

  SUBROUTINE sum_ds_ac( R,ds_ac_tot) !
    ! computes the total length for AC-modulation
    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE (NODE_LAYOUT), pointer :: L
    INTEGER I
    TYPE(INTEGRATION_NODE), POINTER :: T
    REAL(DP), intent(inout) :: ds_ac_tot

    t=>r%t%start
    ds_ac_tot=0.0_dp
    do i=1,r%t%n

       ds_ac_tot=ds_ac_tot+t%ds_ac

       t=>t%next
    enddo
  end SUBROUTINE sum_ds_ac

  SUBROUTINE MAKE_NODE_LAYOUT_2( R,L ) !
    ! Creates a thin layout.
    ! At this point large patches would be problematic; i.e. |d(3)|>>>0 .

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE (NODE_LAYOUT), pointer :: L
    TYPE(FIBRE), POINTER :: P
    INTEGER I,J,k,TEAPOT_LIKE
    REAL(DP) S,DLD,DL,LI,SL
    LOGICAL(LP) CIRCULAR
    TYPE(INTEGRATION_NODE), POINTER :: T1,T2,TM

    CASE_NAME(CASEP1)="THE ENTRANCE PATCH"
    CASE_NAME(CASEP2)="THE EXIT PATCH  "
    CASE_NAME(CASE1)= "THE ENTRANCE FRINGE"
    CASE_NAME(CASE2)="THE EXIT FRINGE"
    CASE_NAME(CASE0)="A STEP OF INTEGRATION BODY"
    !    CASE_NAME(CASE3)="BEAM BEAM PANCAKE"

    if(associated(L)) then
       CALL kill_NODE_LAYOUT(L)
       NULLIFY(L);
    endif

    allocate(L)
    CALL Set_Up_NODE_LAYOUT( L )
    S=0.0_dp
    SL=0.0_dp  !  INTEGRATION LENGTH
    P=>R%START
    k=1
    DO I=1,R%N

       TEAPOT_LIKE=0
       IF(P%MAG%P%B0/=0.0_dp) TEAPOT_LIKE=1
       IF(P%MAG%KIND==KIND16.OR.P%MAG%KIND==KIND16)TEAPOT_LIKE=0
       IF(P%MAG%KIND==KIND0.AND.P%MAG%P%NST/=1) THEN
          WRITE(6,*) "MARKER SHOULD HAVE NST=1 OTHERWISE PROBLEMS "
          WRITE(6,*) "WILL OCCUR WITH THE WORM AND THE NODE_LAYOUT SURVEY "
          STOP 500
       ENDIF
       IF(P%DIR==1) THEN
          LI=0.0_dp;
       ELSE
          LI=P%MAG%L;
       ENDIF
       DL=P%DIR*P%MAG%L/P%MAG%P%NST
       DLD=P%MAG%P%LD/P%MAG%P%NST
       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=0.0_dp;    ! s(1) total ld
       L%END%S(5)=0.0_dp;L%END%DS_AC=0.0_dp;
       T1=>L%END                            ! s(2) local integration distance
       ! s(3) total integration distance
       L%END%CAS=CASEP1                       ! s(4) end of step =  DL
       L%END%pos_in_fibre=1
       L%END%pos=k;k=k+1;
       L%END%PARENT_NODE_LAYOUT=>L
       L%END%PARENT_FIBRE=>P

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=0.0_dp;L%END%S(5)=0.0_dp;L%END%DS_AC=0.0_dp;
       L%END%CAS=CASE1
       L%END%pos_in_fibre=2
       L%END%pos=k;k=k+1;
       L%END%PARENT_NODE_LAYOUT=>L
       L%END%PARENT_FIBRE=>P

       DO J=1,P%MAG%P%NST
          CALL APPEND_EMPTY_THIN( L )
          L%END%TEAPOT_LIKE=TEAPOT_LIKE
          L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=DL;L%END%S(5)=DLD;L%END%ds_ac=DLD;
          L%END%CAS=CASE0
          L%END%pos_in_fibre=J+2
          L%END%pos=k;k=k+1;
          L%END%PARENT_NODE_LAYOUT=>L
          L%END%PARENT_FIBRE=>P
          S=S+DLD
          LI=LI+DL
          SL=SL+P%DIR*DL
          IF(MOD(P%MAG%P%NST,2)==0) THEN
             IF(J==P%MAG%P%NST/2) TM=>L%END    !+1
          ELSE
             IF(J==1) TM=>T1
          ENDIF
       ENDDO

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=0.0_dp;L%END%S(5)=0.0_dp;L%END%DS_AC=0.0_dp;
       L%END%CAS=CASE2
       L%END%pos_in_fibre=P%MAG%P%NST+3
       L%END%pos=k;k=k+1;
       L%END%PARENT_FIBRE=>P

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=0.0_dp;L%END%S(5)=0.0_dp;L%END%DS_AC=0.0_dp;
       L%END%CAS=CASEP2
       L%END%pos_in_fibre=P%MAG%P%NST+4
       L%END%pos=k;k=k+1;
       L%END%PARENT_NODE_LAYOUT=>L
       L%END%PARENT_FIBRE=>P
       T2=>L%END

       P%T1=>T1
       P%T2=>T2
       P%TM=>TM

       P=>P%NEXT
    ENDDO
    L%N=k-1

    L%PARENT_LAYOUT=>R

    IF(R%CLOSED) THEN
       l%closed=.true.
       CIRCULAR=.TRUE.
       CALL RING_L_THIN(L,CIRCULAR)
    ENDIF

   if(lielib_print(12)==1)  call stat_NODE_LAYOUT(l)

  END SUBROUTINE MAKE_NODE_LAYOUT_2


  SUBROUTINE stat_NODE_LAYOUT( L )
    implicit none
    TYPE (NODE_LAYOUT), pointer :: L

    WRITE(6,*)  " PARENT LAYOUT NAME :", L%PARENT_LAYOUT%NAME(1:len_trim(L%PARENT_LAYOUT%NAME))
    WRITE(6,*) " NUMBER OF ORIGINAL LAYOUT ELEMENTS :", L%PARENT_LAYOUT%N
    WRITE(6,*) " NUMBER OF THIN OBJECTS :", L%N
    WRITE(6,*) " TOTAL IDEAL LENGTH OF STRUCTURE :", L%END%S(1)
    WRITE(6,*) " TOTAL INTEGRATION LENGTH OF STRUCTURE (mad8 style survey) :", L%END%S(3)

  end SUBROUTINE stat_NODE_LAYOUT


  SUBROUTINE DRIFT_TO_TIME(T,YL,DT,X,k)
    ! Drifts to a given time using either a regular or TEAPOT drift. (Cylindrical coordinate drift)
    ! The amount of the drift YL is computed to achieve a time DT.

    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(INOUT):: YL
    real(dp),INTENT(IN):: DT
    TYPE(INTEGRATION_NODE), pointer :: T
    TYPE(magnet_chart), pointer :: p
    real(dp) PZ
    real(dp)  b
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    P=>T%PARENT_FIBRE%MAG%P

    IF(k%TIME) then
       B=P%BETA0
    ELSE
       B=1.0_dp
    ENDIF


    !   if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then
    !

    !       PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
    !       R=one/P%B0
    !       A=(X(1)+R)*(one/b+x(5))/PZ
    !       A=ATAN(DT/(A+DT*X(2)/PZ) )
    !       YL=A/P%B0
    !       PT=one-X(2)*TAN(A)/PZ
    !       XN(1)=(X(1)+R)/COS(A)/PT-R
    !       XN(2)=X(2)*COS(A)+SIN(A)*PZ
    !       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
    !       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
    !       X(1)=XN(1)
    !       X(2)=XN(2)
    !       X(3)=XN(3)
    !       X(6)=XN(6)
    !    else


    PZ=ROOT(1.0_dp+2.0_dp*X(5)/b+x(5)**2-X(2)**2-X(4)**2)

    YL=DT/(1.0_dp/b+X(5))*PZ

    X(1)=X(1)+YL*X(2)/PZ
    X(3)=X(3)+YL*X(4)/PZ
    X(6)=X(6)+YL*(1.0_dp/b+X(5))/PZ



    !    endif



  END SUBROUTINE DRIFT_TO_TIME

  SUBROUTINE DRIFTr_BACK_TO_POSITION(T,YL,X,k)
    ! This is a regular drift
    ! It is used in time tracking to project back to the beginning of the thin lens
    ! and it is used in S tracking to drift in the middle of a step.

    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: YL
    TYPE(INTEGRATION_NODE), pointer :: T
    !    TYPE(magnet_chart), pointer :: p
    !    TYPE(magnet_chart), pointer :: p
    real(dp) PZ
    real(dp)  b
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    !    P=>T%PARENT_FIBRE%MAG%P

    IF(k%TIME) then
       B=T%PARENT_FIBRE%BETA0
    ELSE
       B=1.0_dp
    ENDIF


    !    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then



    !      PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
    !      R=one/P%B0
    !      A=-YL*P%B0

    !       PT=one-X(2)*TAN(A)/PZ

    !       XN(1)=(X(1)+R*(two*sin(a/two)**2+X(2)*sin(A)/PZ))/COS(A)/PT
    !       XN(2)=X(2)*COS(A)+SIN(A)*PZ
    !       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
    !       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
    !          WRITE(6,*) "XN(6)-X(6)-DT , DT",DT,XN(6)-X(6)-DT
    !      X(1)=XN(1)
    !      X(2)=XN(2)
    !      X(3)=XN(3)
    !      X(6)=XN(6)
    !   else
    !       CALL DRIFT(YL,DL,P%beta0,1,P%EXACT,P%TIME,X)
    PZ=ROOT(1.0_dp+2.0_dp*X(5)/b+x(5)**2-X(2)**2-X(4)**2)


    X(1)=X(1)-YL*X(2)/PZ
    X(3)=X(3)-YL*X(4)/PZ
    X(6)=X(6)-YL*(1.0_dp/b+X(5))/PZ

    !  endif



  END SUBROUTINE DRIFTr_BACK_TO_POSITION


  SUBROUTINE DRIFTp_BACK_TO_POSITION(T,YL,X,k)
    ! This is a regular drift
    ! It is used in time tracking to project back to the beginning of the thin lens
    ! and it is used in S tracking to drift in the middle of an step.

    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: YL
    TYPE(INTEGRATION_NODE), pointer :: T
    TYPE(magnet_chart), pointer :: p
    type(real_8) XN(6),PZ,PT
    real(dp)  b
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    P=>T%PARENT_FIBRE%MAG%P

    IF(k%TIME) then
       B=P%BETA0
    ELSE
       B=1.0_dp
    ENDIF

    call alloc(xn,6); call alloc(pz,pt);

    !    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then



    !       PZ=sqrt(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
    !       R=one/P%B0
    !       A=-YL*P%B0

    !      PT=one-X(2)*TAN(A)/PZ

    !      XN(1)=(X(1)+R*(two*sin(a/two)**2+X(2)*sin(A)/PZ))/COS(A)/PT
    !      XN(2)=X(2)*COS(A)+SIN(A)*PZ
    !      XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
    !      XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
    !          WRITE(6,*) "XN(6)-X(6)-DT , DT",DT,XN(6)-X(6)-DT
    !      X(1)=XN(1)
    !      X(2)=XN(2)
    !      X(3)=XN(3)
    !      X(6)=XN(6)
    !   else
    !       CALL DRIFT(YL,DL,P%beta0,1,P%EXACT,P%TIME,X)
    PZ=sqrt(1.0_dp+2.0_dp*X(5)/b+x(5)**2-X(2)**2-X(4)**2)


    X(1)=X(1)-YL*X(2)/PZ
    X(3)=X(3)-YL*X(4)/PZ
    X(6)=X(6)-YL*(1.0_dp/b+X(5))/PZ

    !   endif

    call kill(xn,6); call kill(pz,pt);



  END SUBROUTINE DRIFTp_BACK_TO_POSITION


  !  Survey still worm like

  SUBROUTINE FILL_SURVEY_DATA_IN_NODE_LAYOUT(R)
    ! THIS SUBROUTINE ALLOCATES NODE FRAMES IF NEEDED
    ! IT SURVEYS THE NODES USING THE OLD REAL WORMS
    ! SHOULD BE CALLED AFTER MISALIGNMENTS OR MOVING PART OF LATTICE

    IMPLICIT NONE
    type(layout),target:: r
    type(fibre), pointer ::c
    type(INTEGRATION_NODE), pointer ::t
    type(worm) vers
    integer k,ic,j
    real(dp) x(6),ent(3,3),a(3)
    LOGICAL(LP) APER
    aper=APERTURE_FLAG
    APERTURE_FLAG=.FALSE.

    if(.not.associated(r%t)) call MAKE_NODE_LAYOUT(r)

    CALL  allocate_node_frame( R)

    call survey(r)

    CALL ALLOC(vers,r)
    C=>r%START
    CALL XFRAME(vers%E,C%chart%f%ent,C%chart%f%A,-7)  ! initializes the survey part of worm
    vers%E%L(-1)=0.d0 !Starts beam line at z=0   fake distance along ld for cheap work

    do k=1,r%n
       x=0.0_dp
       CALL TRACK(r,x,k,k+1,default,vers)
       if(.not.check_stable) then
          Write(6,*) " fake instability at ",c%mag%name, " ",k
          check_stable=.true.
          CALL RESET_APERTURE_FLAG
       endif

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


          !  omega(1)= a(1)+scale*(xr(1)*ent(1,1)+xr(3)*ent(2,1))
          !  omega(2)= a(2)+scale*(xr(1)*ent(1,2)+xr(3)*ent(2,2))
          !  omega(3)= a(3)+scale*(xr(1)*ent(1,3)+xr(3)*ent(2,3))
          !  r1(1)=omega0(3)
          !  r1(2)=omega0(1)
          !  r2(1)=omega(3)
          !  r2(2)=omega(1)
          !  if(abs(r1(1))>1.d6) r1=r2
          !  call gMoveTo2D(r1(1),r1(2))
          !  call  gDrawLineTo2D(r2(1),r2(2))
          !  omega0=omega
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
          write(6,*)k, ic,c%mag%name,c%mag%p%nst
          stop 888
       endif
       c=>c%next
    enddo

    CALL kill(vers)

    APERTURE_FLAG=aper

  end  subroutine fill_survey_data_in_NODE_LAYOUT


  ! BEAM STUFF

  subroutine create_beam(B,N,CUT,SIG,A,C,T)
    USE gauss_dis
    implicit none
    INTEGER N,I,J,K
    REAL(DP) CUT,SIG(6),X
    TYPE(BEAM) B
    REAL(DP), OPTIONAL :: A(6,6),C(6)
    REAL(DP)  AS(6,6),CL(6),XT(6)
    TYPE (INTEGRATION_NODE),optional,target::  T

    IF(.NOT.ASSOCIATED(B%N)) THEN
       CALL ALLOCATE_BEAM(B,N)
    ELSEIF(B%N/=N) THEN
       CALL KILL_BEAM(B)
       CALL ALLOCATE_BEAM(B,N)
    ENDIF

    CL=0.0_dp; AS=0.0_dp;
    DO I=1,6
       AS(I,I)=1.0_dp
    ENDDO

    IF(PRESENT(A)) AS=A
    IF(PRESENT(C)) CL=C

    DO I=1,N
       DO J=1,6
          CALL GRNF(X,cut)
          XT(J)=X*SIG(J)
       ENDDO
       B%X(I,1:6)=CL(:)
       DO J=1,6
          DO K=1,6
             B%X(I,J)=AS(J,K)*XT(K)+B%X(I,J)
          ENDDO

       ENDDO
    ENDDO


    if(present(t)) then
       DO I=1,N
          ! if(associated(B%POS(I)%NODE))then
          B%POS(I)%NODE=>T
          ! endif
       ENDDO

    endif

  end    subroutine create_beam

  subroutine create_PANCAKE(B,N,CUT,SIG,T,A)
    USE gauss_dis
    implicit none
    INTEGER N,I,J
    REAL(DP) CUT,SIG(6),X,Y(LNV),beta(2)
    TYPE(BEAM) B
    TYPE (INTEGRATION_NODE),optional,target::  T
    TYPE (DAMAP),OPTIONAL :: A
    TYPE (tree) monkey

    IF(.NOT.ASSOCIATED(B%N)) THEN
       CALL ALLOCATE_BEAM(B,N)
    ELSEIF(B%N/=N) THEN
       CALL KILL_BEAM(B)
       CALL ALLOCATE_BEAM(B,N)
    ENDIF
    write(6,*) n," particles created"
    Y=0.0_dp
    IF(.not.PRESENT(A)) THEN

       DO I=1,N
          DO J=1,6
             CALL GRNF(X,cut)
             B%X(I,J)=X*SIG(J)
          ENDDO
          B%X(I,7)=0.0_dp
       enddo
    ELSE
       call alloc(monkey)
       beta(1)=(a%v(1).sub.'1')**2+(a%v(1).sub.'01')**2
       beta(2)=(a%v(3).sub.'001')**2+(a%v(3).sub.'0001')**2
       write(6,*) " Betas in create_PANCAKE ",beta
       monkey=A
       DO I=1,N
          DO J=1,C_%ND
             CALL GRNF(X,cut)
             Y(2*j-1)=X*sqrt(SIG(j)/2.0_dp)
             CALL GRNF(X,cut)
             Y(2*j)=X*sqrt(SIG(j)/2.0_dp)
          ENDDO
          y=monkey*Y
          B%X(I,1:C_%ND2)=y(1:c_%nd2)

          DO J=C_%ND2+1,6
             CALL GRNF(X,cut)
             B%X(I,J)=X*SIG(J)
          ENDDO
          B%X(I,7)=0.0_dp
       enddo
       CALL KILL(MONKEY)
    ENDIF




    if(present(t)) then
       DO I=1,N
          if(associated(B%POS(I)%NODE))then
             B%POS(I)%NODE=>T
          endif
       ENDDO

    endif
  end    subroutine create_PANCAKE

  subroutine copy_beam(B1,B2)
    implicit none
    INTEGER I
    TYPE(BEAM), INTENT(INOUT) :: B1,B2

    IF(.NOT.ASSOCIATED(B2%N)) THEN
       CALL ALLOCATE_BEAM(B2,B1%N)
    ELSEIF(B1%N/=B2%N) THEN
       CALL KILL_BEAM(B2)
       CALL ALLOCATE_BEAM(B2,B1%N)
    ENDIF

    B2%X=B1%X
    B2%U=B1%U
    B2%N=B1%N
    !    B2%CHARGE=B1%CHARGE
    B2%LOST=B1%LOST
    DO I=0,B1%N
       if(associated(B1%POS(I)%NODE))then
          B2%POS(I)%NODE=>B1%POS(I)%NODE
       endif
    ENDDO

  END subroutine copy_beam

  subroutine READ_beam_raw(B,MF)
    implicit none
    INTEGER k,mf
    TYPE(BEAM), INTENT(IN):: B

    DO K=1,b%n
       IF(.not.B%U(K)) THEN
          if(associated(b%pos(k)%NODE)) then
             WRITE(MF,100) B%X(K,1:6),b%pos(k)%NODE%s(3)+B%X(K,7)
          else
             WRITE(MF,100) B%X(K,1:6),B%X(K,7)
          endif
       ENDIF
    ENDDO
100 FORMAT(7(1x,e13.6))
  END subroutine READ_beam_raw

  subroutine PRINT_beam_raw(B,MF)
    implicit none
    INTEGER k,mf
    TYPE(BEAM), INTENT(IN):: B

    DO K=1,b%n
       IF(.not.B%U(K)) THEN
          if(associated(b%pos(k)%NODE)) then
             WRITE(MF,100) B%X(K,1:6),b%pos(k)%NODE%s(3)+B%X(K,7)
          else
             WRITE(MF,100) B%X(K,1:6),B%X(K,7)
          endif
       ENDIF
    ENDDO
100 FORMAT(7(1x,e13.6))
  END subroutine PRINT_beam_raw

  subroutine stat_beam_raw(B,n,MF,xm)
    implicit none
    INTEGER i,j,k,mf,NOTlost,N
    TYPE(BEAM), INTENT(IN):: B
    real(dp), optional :: xm(6)
    real(dp), allocatable :: av(:,:)
    real(dp) em(2),beta(2),xma(6)
    allocate(av(n,n))
    av=0.0_dp
    notlost=0
    xma(:)=-1.0_dp
    DO K=1,b%n
       IF(.not.B%U(K)) THEN
          do i=1,6
             if(abs(b%x(k,i))>xma(i)) xma(i)=abs(b%x(k,i))
          enddo

          do i=1,n
             do j=i,n
                av(i,j)= b%x(k,i)*b%x(k,j)+av(i,j)
             enddo
          enddo
          notlost=notlost+1
       ENDIF
    ENDDO
    IF(NOTLOST==0) THEN
       if(mf/=6) then
          WRITE(mf,*) " ALL PARTICLES ARE LOST "
          WRITE(mf,*) " NO STATISTICS "
       else
          WRITE(6,*) " ALL PARTICLES ARE LOST "
          WRITE(6,*) " NO STATISTICS "
       endif
       deallocate(av)
       RETURN
    ENDIF
    if(notlost/=b%n-b%lost) then
       Write(6,*) " Error keeping track of lost particles "
       stop 999
    endif

    WRITE(MF,*) " NUMBER LEFT ",B%N-B%LOST
    if(mf/=6)WRITE(6,*) " NUMBER LEFT ",B%N-B%LOST
    WRITE(MF,*) " LOST ",B%LOST
    if(mf/=6)WRITE(6,*) " LOST ",B%LOST
    av=av/notlost
    em(1)=2.0_dp*sqrt(av(1,1)*av(2,2)-av(1,2)**2)
    em(2)=2.0_dp*sqrt(av(3,3)*av(4,4)-av(3,4)**2)
    beta(1)=2.0_dp*av(1,1)/em(1)
    beta(2)=2.0_dp*av(3,3)/em(2)

    write(mf,*) " average arrays "
    write(mf,*) "betas ",beta
    write(mf,*) "emittances ",em
    if(mf/=6) then
       write(6,*) " average arrays "
       write(6,*) "betas ",beta
       write(6,*) "emittances ",em
    endif
    write(6,*) " limits "
    write(6,*) xma(1:2)
    write(6,*) xma(3:4)
    write(6,*) xma(5:6)
    if(present(xm)) xm=xma

100 FORMAT(7(1x,e13.6))

    deallocate(av)
  END subroutine stat_beam_raw


  subroutine PRINT_beam(B,MF,I)
    implicit none
    INTEGER K,MF,I1,I2
    INTEGER,OPTIONAL:: I
    TYPE(BEAM), INTENT(IN):: B
    TYPE(INTEGRATION_NODE),POINTER::T
    TYPE(FIBRE),POINTER::F

    I1=1
    I2=B%N

    IF(PRESENT(I)) THEN
       I1=I
       I2=I
    ENDIF
    !    IF(B%TIME_INSTEAD_OF_S) THEN
    !       WRITE(MF,*) "____________________________ TIME TRACKED BEAM __________________________________"
    !    ELSE
    WRITE(MF,*) "_________________ POSITION TRACKED BEAM (AS IN PTC PROPER)_______________________"
    !    ENDIF

    DO K=I1,I2
       IF(B%U(K)) THEN
          WRITE(MF,*) " PARTICLE # ",K, " IS LOST "
       ELSE
          T=>B%POS(K)%NODE
          F=>T%PARENT_FIBRE
          WRITE(MF,*) "_________________________________________________________________________"
          WRITE(MF,*) " PARTICLE # ",K, " IS LOCATED AT SLICE # ",T%POS," IN FIBRE  ",F%MAG%NAME
          WRITE(MF,*) " IN THE FIBRE POSITION  ",T%pos_in_fibre
          WRITE(MF,*) " IN ",CASE_NAME(T%CAS)
          IF(T%CAS==CASE0)WRITE(MF,*) " AT THE STEP NUMBER ",T%pos_in_fibre-2

          WRITE(MF,*) "........................................................................."
          !          IF(B%TIME_INSTEAD_OF_S) THEN
          !             WRITE(MF,*) " TIME AND POSITION AFTER THIN SLICE = ",B%X(K,6:7)
          !          ELSE
          WRITE(MF,*) " TIME AND POSITION  = ",B%X(K,6:7)
          !          ENDIF
          WRITE(MF,*) " X,Y = ",B%X(K,1),B%X(K,3)
          WRITE(MF,*) " PX,PY = ",B%X(K,2),B%X(K,4)
          WRITE(MF,*) " ENERGY VARIABLE = ",B%X(K,5)
       ENDIF
       WRITE(MF,*) "_________________________________________________________________________"
    ENDDO

  END subroutine PRINT_beam



  SUBROUTINE NULLIFY_BEAM(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    NULLIFY(B%N,B%LOST)
    !    NULLIFY(B%Y)
    NULLIFY(B%X)
    NULLIFY(B%U)
    NULLIFY(B%POS)
    !    NULLIFY(B%CHARGE)
    !    NULLIFY(B%TIME_INSTEAD_OF_S)
    !    NULLIFY(B%SIGMA)
    !    NULLIFY(B%DX,B%ORBIT)
    !    NULLIFY(B%BBPAR,B%BEAM_BEAM,B%BBORBIT)
  END SUBROUTINE NULLIFY_BEAM

  SUBROUTINE NULLIFY_BEAMS(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B(:)
    INTEGER I
    DO I=1,SIZE(B)
       CALL NULLIFY_BEAM(B(i))
    ENDDO

  END SUBROUTINE NULLIFY_BEAMS

  subroutine alloc_three_d_info(v)
    IMPLICIT NONE
    TYPE(three_d_info) , INTENT (INOUT) :: V
    v%a=0.0_dp
    v%b=0.0_dp
    v%o=0.0_dp
    v%ent=global_frame
    v%exi=global_frame
    v%mid=global_frame
    v%reference_ray=0.0_dp
    v%r0=0.0_dp
    v%r=0.0_dp
    v%x=0.0_dp
    v%scale=1.0_dp
    v%u=my_false
    v%wx=0.1_dp
    v%wy=0.1_dp

  end subroutine alloc_three_d_info

  SUBROUTINE ALLOCATE_BEAM(B,N)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    INTEGER , INTENT (IN) :: N
    INTEGER I

    ALLOCATE(B%N,B%LOST)

    B%N=N
    B%LOST=0
    !    NULLIFY(B%Y)
    !    IF(PRESENT(POLYMORPH)) THEN
    !       IF(POLYMORPH) then
    !          ALLOCATE(B%Y(6))
    !          CALL ALLOC(B%Y)
    !       endif
    !    ENDIF
    ALLOCATE(B%X(N,7))
    ALLOCATE(B%U(0:N))
    ALLOCATE(B%POS(0:N))
    !    ALLOCATE(B%SIGMA(6))
    !    ALLOCATE(B%DX(3))
    !    ALLOCATE(B%ORBIT(6))
    !    ALLOCATE(B%BBPAR,B%BEAM_BEAM,B%BBORBIT)
    DO I=0,N
       NULLIFY(B%POS(i)%NODE)
    ENDDO
    !   ALLOCATE(B%CHARGE)
    !   ALLOCATE(B%TIME_INSTEAD_OF_S)

    B%X  = 0.0_dp
    B%U  = .FALSE.
    !    B%CHARGE=1
    !    B%TIME_INSTEAD_OF_S=.FALSE.

    !    B%SIGMA=ZERO
    !    B%DX=ZERO
    !    B%BBPAR=ZERO
    !    B%ORBIT=ZERO
    !    B%BEAM_BEAM=MY_FALSE
    !    B%BBORBIT=MY_FALSE
  END SUBROUTINE ALLOCATE_BEAM

  SUBROUTINE KILL_BEAM(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    !    IF(ASSOCIATED(B%Y)) THEN
    !       CALL KILL(B%Y)
    !       DEALLOCATE(B%Y)
    !    ENDIF
    IF(ASSOCIATED(B%N)) THEN
       DEALLOCATE(B%N,B%LOST,B%X,B%U,B%POS)
       !       DEALLOCATE(B%N,B%LOST,B%X,B%U,B%POS,B%CHARGE,B%TIME_INSTEAD_OF_S)
       !      DEALLOCATE(B%SIGMA,B%DX,B%BBPAR,B%ORBIT,B%BEAM_BEAM,B%BBORBIT)
    ENDIF
  END SUBROUTINE KILL_BEAM

  SUBROUTINE KILL_BEAMS(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B(:)
    INTEGER I
    DO I=1,SIZE(B)
       CALL KILL_BEAM(B(i))
    ENDDO
  END SUBROUTINE KILL_BEAMS


  FUNCTION BEAM_IN_X(B,I)
    IMPLICIT NONE
    REAL(DP) BEAM_IN_X(6)
    TYPE(BEAM), INTENT(INOUT) ::B
    INTEGER, INTENT(IN) :: I

    BEAM_IN_X=B%X(I,1:6)

  END  FUNCTION BEAM_IN_X

  SUBROUTINE X_IN_BEAM(B,X,I,DL,T)
    IMPLICIT NONE
    REAL(DP),OPTIONAL:: X(6)
    REAL(DP),OPTIONAL:: DL
    TYPE(BEAM), INTENT(INOUT) ::B
    TYPE(INTEGRATION_NODE),OPTIONAL,POINTER :: T
    INTEGER, INTENT(IN) :: I

    if(PRESENT(X)) B%X(I,1:6)=X(1:6)
    IF(PRESENT(DL)) B%X(I,7)=DL
    IF(PRESENT(T)) B%POS(I)%NODE=>T
    if(.not.CHECK_STABLE) then
       !       write(6,*) "unstable "
       CALL RESET_APERTURE_FLAG
       b%u(I)=.true.
       B%LOST=B%LOST+1
    endif

  END  SUBROUTINE X_IN_BEAM

  !  Beam Beam stuff

  subroutine s_locate_beam_beam(my_ering,sc,pos,tl,b_b)
    implicit none
    type(layout), target :: my_ering
    type(integration_node), pointer :: tl
    real(dp) sc
    integer pos,j
    logical(lp) b_b


    IF(.NOT.ASSOCIATED(my_ering%T)) THEN
       CALL MAKE_NODE_LAYOUT(my_ering)
    ENDIF
    ! s(1) total ld
    ! s(2) local integration distance
    !          SC=MOD(SC,MY_RING%T%END%S(1))
    b_b=.false.
    TL=>my_ering%T%START
    DO j=1,my_ering%T%N
       if(pos<1) then
          IF(TL%S(1)<=SC.AND.TL%NEXT%S(1)>SC) then
             b_b=.true.
             exit
          endif
       else
          if(j==pos) then
             b_b=.true.
             exit
          endif
       endif
       TL=>TL%NEXT
    ENDDO
    if(b_b.and.(tl%cas==case0.or.tl%cas==caset)) then
       write(6,*) " Beam-Beam position at ",tl%parent_fibre%mag%name
       if(.not.associated(tl%BB)) call alloc(tl%BB)
       write(6,*) tl%pos,tl%parent_fibre%mag%name,' created'
       b_b=my_true
    else
       b_b=my_false
       write(6,*) " Beam-Beam position not found "
    endif
  end  subroutine s_locate_beam_beam

  subroutine locate_beam_beam(my_ering,a,ent,tl,b_b)
    implicit none
    type(layout), target :: my_ering
    type(integration_node), pointer :: tl,tmin,tp
    real(dp) dmin,a(3),ent(3,3),d,cos
    integer pos,j
    logical(lp) b_b

    dmin=mybig

    IF(.NOT.ASSOCIATED(my_ering%T)) THEN
       CALL MAKE_NODE_LAYOUT(my_ering)
       call FILL_SURVEY_DATA_IN_NODE_LAYOUT(my_ering)
    ENDIF
    IF(.NOT.ASSOCIATED(my_ering%T%start%B)) THEN
       call FILL_SURVEY_DATA_IN_NODE_LAYOUT(my_ering)
    ENDIF

    b_b=.false.
    TL=>my_ering%T%START
    tmin=>tl
    tp=>tl
    DO j=1,my_ering%T%N
       if(tl%cas==case0.or.tl%cas==caset) then
          d=sqrt((a(1)-tl%a(1))**2+(a(2)-tl%a(2))**2+(a(3)-tl%a(3))**2)
          if(d<=dmin) then
             dmin=d
             tmin=>tl
             tp=>tl%parent_fibre%previous%t2%previous%previous
          endif
       endif
       TL=>TL%NEXT
    ENDDO
    if((tmin%cas==case0.or.tmin%cas==caset)) then

       write(6,*) " Tentative Beam-Beam position at ",tl%parent_fibre%mag%name
       write(6,*) tmin%pos,tmin%parent_fibre%mag%name,' created'
       b_b=my_true

       cos=0.0_dp
       do j=1,3
          cos=cos+ (a(j)-tmin%a(j))*tmin%ent(1,j)
       enddo
       tl=>tp
       if(cos<0.0_dp) then
          write(6,*) " Beam-Beam position replaced at ",tl%parent_fibre%mag%name,tl%cas
          write(6,*) tl%pos,tl%parent_fibre%mag%name,' created'
       endif
    else
       b_b=my_false
       write(6,*) " Beam-Beam position not found "
    endif
    if(b_b.and.(.not.associated(tl%BB))) call alloc(tl%BB)

  end  subroutine locate_beam_beam

end module ptc_multiparticle
