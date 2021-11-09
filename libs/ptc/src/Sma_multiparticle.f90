
!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module ptc_multiparticle
  !  use S_TRACKING  !,FRINGE_=>FRINGE__MULTI !,FACE=>FACE_MULTI
  USE S_FIBRE_BUNDLE

  use beam_beam_ptc
  implicit none
  public
  CHARACTER*27 CASE_NAME(-2:3)
  PRIVATE fuzzy_eq,fuzzy_neq
  PRIVATE TRACK_FIBRE_FRONTR,TRACK_FIBRE_FRONTP
  PRIVATE TRACK_FIBRE_BACKR,TRACK_FIBRE_BACKP
  PRIVATE TRACK_NODE_SINGLER,TRACK_NODE_SINGLEP,TRACK_NODE_SINGLEV
  private DRIFTr_BACK_TO_POSITION,DRIFTp_BACK_TO_POSITION  !,DRIFT_BACK_TO_POSITION
  private MAKE_NODE_LAYOUT_2 !,DRIFT_TO_TIME
  PRIVATE MODULATE_R,MODULATE_P
  PRIVATE SURVEY_EXIST_PLANAR_L_NEW,SURVEY_FIBRE_new,SURVEY_EXIST_PLANAR_IJ_new
  private survey_integration_layout
  PRIVATE TRACK_MODULATION_R,TRACK_MODULATION_P,FIND_PATCH_0_survey
  logical :: old_survey=.true.
  private TRACK_NODE_SINGLE_quar,TRACK_NODE_SINGLE_quaP
!!!!!!!  Old Survey   !!!!!!!!
private MISALIGN_FIBRE_EQUAL


 LOGICAL :: no_mis=.TRUE. 
  !  LOGICAL :: OLD_MOD=.TRUE.

  logical(lp),private, parameter :: dobb=.true.
  logical(lp),private            :: aperture_all_case0=.false.
 ! type(probe) :: xsm,xsm0
  real(dp) :: xsm0t=0.0_dp,xsmt=0.0_dp
  !real(dp) :: unit_time =1.0e-3_dp
  REAL(dp) :: x_orbit_sync(6)= 0.0_dp,dt_orbit_sync=0.0_dp



  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE MISALIGN_FIBRE_EQUAL
  END  INTERFACE

 ! INTERFACE survey_new

 ! END INTERFACE

  INTERFACE TRACK_NODE_SINGLE
     MODULE PROCEDURE TRACK_NODE_SINGLE_quar    !@1  t,x,state,charge
     MODULE PROCEDURE TRACK_NODE_SINGLE_quaP    !@1  t,x,state,charge
  END INTERFACE TRACK_NODE_SINGLE



  INTERFACE TRACK_NODE_SINGLE

     MODULE PROCEDURE TRACK_NODE_SINGLER     !@1  t,x,state,charge
     MODULE PROCEDURE TRACK_NODE_SINGLEP     !@1  t,y,state,charge
     MODULE PROCEDURE TRACK_NODE_SINGLEV     !@1  t,v,state,charge
  END INTERFACE

  INTERFACE convert_bmad_to_ptc
     MODULE PROCEDURE convert_bmad_to_ptcr
     MODULE PROCEDURE convert_bmad_to_ptcp
     MODULE PROCEDURE convert_bmad_to_ptcar
     MODULE PROCEDURE convert_bmad_to_ptcap
  END INTERFACE

  INTERFACE convert_ptc_to_bmad
     MODULE PROCEDURE convert_ptc_to_bmadr
     MODULE PROCEDURE convert_ptc_to_bmadp
     MODULE PROCEDURE convert_ptc_to_bmadar
     MODULE PROCEDURE convert_ptc_to_bmadap
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

  INTERFACE survey 
    MODULE PROCEDURE SURVEY_EXIST_PLANAR_L_NEW
!    MODULE PROCEDURE SURVEY_FIBRE_new
    MODULE PROCEDURE SURVEY_EXIST_PLANAR_IJ_new

   MODULE PROCEDURE survey_integration_layout
  end INTERFACE survey 

  INTERFACE FIND_PATCH_with_survey 
    MODULE PROCEDURE FIND_PATCH_0_survey

  end INTERFACE FIND_PATCH_with_survey 


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
    integer(2) n

    EL=>C%PARENT_FIBRE%MAG
    ELP=>C%PARENT_FIBRE%MAGP

    IF(K%MODULATION) THEN
      n=el%slow_ac
      
      if (modulationtype == 1) then
         V=zero
         DV=el%D_ac*XS%AC(n)%X(2)
      else 

         DV=(XS%AC(n)%X(1)*COS(EL%theta_ac)-XS%AC(n)%X(2)*SIN(EL%theta_ac))
         V=EL%DC_ac+EL%A_ac*DV
         DV=el%D_ac*DV
      endif   
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
    CALL transfer_ANBN(EL,ELP,VR=V,DVR=DV,t=T)


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
           a%table(0)%bn(i)= ((a%table(it+1)%bn(i)-a%table(it)%bn(i))*rat + a%table(it)%bn(i))*a%r
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
    integer(2) n

    EL=>C%PARENT_FIBRE%MAG
    ELP=>C%PARENT_FIBRE%MAGP

    CALL ALLOC(V)
    CALL ALLOC(DV)

    IF(K%MODULATION) THEN
      n=el%slow_ac
      if (modulationtype == 1) then
        V=zero
        DV=el%D_ac*XS%AC(n)%X(2)
      else
        DV=(XS%AC(n)%X(1)*COS(ELP%theta_ac)-XS%AC(n)%X(2)*SIN(ELP%theta_ac))
        V=ELP%DC_ac+ELP%A_ac*DV
        DV=elp%D_ac*DV
      endif
      
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
    integer(2) n
    if(xs%nac==0) return
    do n=1,xs%nac
      if(k%time) then
         beta0=>C%PARENT_FIBRE%beta0
         xs%ac%t=c%DS_AC/beta0+xs%ac(n)%t
         xt = cos(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(1) + sin(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(2)
         XS%AC(n)%X(2) = -sin(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(1) + cos(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(2)
         XS%AC(n)%X(1) = xt
      else
         xt = cos(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(1) + sin(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(2)
         XS%AC(n)%X(2) = -sin(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(1) + cos(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(2)
         XS%AC(n)%X(1) = xt
         xs%ac(n)%t=c%DS_AC+xs%ac(n)%t
      endif
    enddo
  END   SUBROUTINE TRACK_MODULATION_R

  SUBROUTINE TRACK_MODULATION_P(C,XS,K)
    IMPLICIT NONE
    type(INTEGRATION_NODE), pointer :: C
    type(probe_8), INTENT(INOUT) :: xs
    TYPE(INTERNAL_STATE) K
    TYPE(REAL_8) xt
    real(dp),pointer :: beta0
    integer(2) n
    if(xs%nac==0) return
    CALL ALLOC(XT)
    do n=1,xs%nac
        if(k%time) then
           beta0=>C%PARENT_FIBRE%beta0
           xs%ac(n)%t=c%DS_AC/beta0+xs%ac(n)%t
           xt = cos(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(1) + sin(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(2)
           XS%AC(n)%X(2) = -sin(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(1) + cos(XS%AC(n)%om * c%DS_AC/beta0) *XS%AC(n)%X(2)
           XS%AC(n)%X(1) = xt
        else
           xt = cos(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(1) + sin(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(2)
           XS%AC(n)%X(2) = -sin(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(1) + cos(XS%AC(n)%om * c%DS_AC) *XS%AC(n)%X(2)
           XS%AC(n)%X(1) = xt
           xs%ac(n)%t=c%DS_AC+xs%ac(n)%t
        endif
    enddo
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

    real(dp), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    real(dp) b1

    !FRONTAL PATCH
    !    IF(ASSOCIATED(C%PATCH)) THEN
    PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    !    ELSE
    !       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    !    ENDIF

    ! PUSHING BEAM
    !

    b1=C%BETA0

    IF(PATCHE/=0.AND.PATCHE/=2.AND.PATCHE/=5) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN).and.PATCHE/=4) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0.or.CN%PATCH%ENERGY==1.or.CN%PATCH%ENERGY==4) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAG%P%P0C
             B0=>CN%beta0  !    CN%MAG%P%BETA0   date 2021.7.1
 
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
       else  ! associated   
             P0=>C%PATCH%P0b
             B0=>C%PATCH%B0b
 
             X(2)=X(2)*P0/C%MAG%P%P0C
             X(4)=X(4)*P0/C%MAG%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=root(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/C%MAG%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAG%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAG%P%P0C-1.0_dp
             ENDIF           
       ENDIF ! ASSOCIATED
 
    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG==1.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%a_T  !/c%beta0
      else
       X(6)=X(6)-C%PATCH%a_L 
      endif
    ENDIF

    CALL DTILTD(C%MAG%P%TILTD,1,X)
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
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    real(dp) b1

    !FRONTAL PATCH
    !    IF(ASSOCIATED(C%PATCH)) THEN
    PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    !    ELSE
    !       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    !    ENDIF

    ! PUSHING BEAM
    !
    b1=C%BETA0


    IF(PATCHE/=0.AND.PATCHE/=2.AND.PATCHE/=5) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%PREVIOUS
       IF(ASSOCIATED(CN).and.PATCHE/=4) THEN ! ASSOCIATED
          !          IF(.NOT.CN%PATCH%ENERGY) THEN     ! No need to patch IF PATCHED BEFORE
          IF(CN%PATCH%ENERGY==0.or.CN%PATCH%ENERGY==1.or.CN%PATCH%ENERGY==4) THEN     ! No need to patch IF PATCHED BEFORE
             P0=>CN%MAGP%P%P0C
             B0=>CN%beta0  !    CN%MAG%P%BETA0   date 2021.7.1

           !  B0=>CN%MAGP%P%BETA0
 
 
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
       else  ! associated 
             P0=>C%PATCH%P0b
             B0=>C%PATCH%B0b
 
             X(2)=X(2)*P0/C%MAGP%P%P0C
             X(4)=X(4)*P0/C%MAGP%P%P0C
             IF(k%TIME.or.recirculator_cheat)THEN
                X(5)=SQRT(1.0_dp+2.0_dp*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAGP%P%P0C-1.0_dp !X(5) = DP/P0C_NEW
                X(5)=(2.0_dp*X(5)+X(5)**2)/(SQRT(1.0_dp/C%MAGP%P%BETA0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/C%MAGP%P%BETA0)
             ELSE
                X(5)=(1.0_dp+X(5))*P0/C%MAGP%P%P0C-1.0_dp
             ENDIF           
       ENDIF ! ASSOCIATED
 
    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG==1.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%a_T    !/c%beta0
      else
       X(6)=X(6)-C%PATCH%a_L
      endif
    ENDIF

    CALL DTILTD(C%MAGP%P%TILTD,1,X)
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
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    real(dp) b1

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
    CALL DTILTD(C%MAG%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%b_T   !/c%beta0
      else
       X(6)=X(6)-C%PATCH%b_L
      endif
    ENDIF

    IF(PATCHG==2.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely
    b1=C%BETA0

    IF(PATCHE/=0.AND.PATCHE/=1.AND.PATCHE/=4) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT 
!       IF(.NOT.ASSOCIATED(CN)) CN=>C
       IF(ASSOCIATED(CN).AND.PATCHE/=5) then
       !       P0=>CN%MAG%P%P0C
       !       B0=>CN%MAG%P%BETA0
       P0=>CN%MAG%P%P0C
       B0=>CN%BETA0
       b1=b0
       X(2)=X(2)*C%MAG%P%P0C/P0
       X(4)=X(4)*C%MAG%P%P0C/P0
       IF(k%TIME.or.recirculator_cheat)THEN
          X(5)=root(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW
          X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp
       ENDIF
    
    else
             P0=>C%PATCH%P0b
             B0=>C%PATCH%B0b
             b1=b0
             X(2)=X(2)*C%MAG%P%P0C/P0   ! 8/31/2016
             X(4)=X(4)*C%MAG%P%P0C/P0   ! 8/31/2016
             IF(k%TIME.or.recirculator_cheat)THEN   ! 8/31/2016
              X(5)=root(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD   ! 8/31/2016
              X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW   ! 8/31/2016
              X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)   ! 8/31/2016
             ELSE
               X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp   ! 8/31/2016
             ENDIF      
    endif
ENDIF

 
  
  END SUBROUTINE TRACK_FIBRE_BACKR

  SUBROUTINE TRACK_FIBRE_BACKP(C,X,K)
    implicit none
    logical(lp) :: doneitf=.false.
    TYPE(FIBRE),TARGET,INTENT(INOUT):: C
    type(real_8), INTENT(INOUT) :: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    logical(lp) ou,patch
    INTEGER(2) PATCHT,PATCHG,PATCHE
    TYPE (fibre), POINTER :: CN
    real(dp), POINTER :: P0,B0
    real(dp) b1

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
    CALL DTILTD(C%MAGP%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(K%TOTALPATH==0)) THEN
      if(K%time) then
       X(6)=X(6)-C%PATCH%b_T   !/c%beta0
      else
       X(6)=X(6)-C%PATCH%b_L
      endif
    ENDIF

    IF(PATCHG==2.or.PATCHG==3) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,k,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely
    b1=C%BETA0
    IF(PATCHE/=0.AND.PATCHE/=1.AND.PATCHE/=4) THEN

       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
!       IF(.NOT.ASSOCIATED(CN)) CN=>C
       IF(ASSOCIATED(CN).and.PATCHE/=5) then
       !       P0=>CN%MAGP%P%P0C
       !       B0=>CN%MAGP%P%BETA0
       P0=>CN%MAGP%P%P0C
       B0=>CN%BETA0
       b1=b0
       X(2)=X(2)*C%MAGP%P%P0C/P0
       X(4)=X(4)*C%MAGP%P%P0C/P0
       IF(k%TIME.or.recirculator_cheat)THEN
          X(5)=sqrt(1.0_dp+2.0_dp*X(5)/C%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAGP%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW
          X(5)=(2.0_dp*X(5)+X(5)**2)/(sqrt(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)
       ELSE
          X(5)=(1.0_dp+X(5))*C%MAGP%P%P0C/P0-1.0_dp
       ENDIF
    
    else
             P0=>C%PATCH%P0b
             B0=>C%PATCH%B0b
             b1=b0
             X(2)=X(2)*C%MAG%P%P0C/P0   ! 8/31/2016
             X(4)=X(4)*C%MAG%P%P0C/P0   ! 8/31/2016
             IF(k%TIME.or.recirculator_cheat)THEN   ! 8/31/2016
              X(5)=sqrt(1.0_dp+2.0_dp*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD   ! 8/31/2016
              X(5)=X(5)*C%MAG%P%P0C/P0-1.0_dp !X(5) = DP/P0C_NEW   ! 8/31/2016
              X(5)=(2.0_dp*X(5)+X(5)**2)/(sqrt(1.0_dp/B0**2+2.0_dp*X(5)+X(5)**2)+1.0_dp/B0)   ! 8/31/2016
             ELSE
               X(5)=(1.0_dp+X(5))*C%MAG%P%P0C/P0-1.0_dp   ! 8/31/2016
             ENDIF               
    ENDIF
    
endif
 

  END SUBROUTINE TRACK_FIBRE_BACKP


  ! thin lens tracking

  SUBROUTINE TRACK_NODE_SINGLEV(T,V,K) !!
    implicit none
    TYPE(INTEGRATION_NODE),POINTER :: T
    TYPE(INTERNAL_STATE)  K
    REAL(DP) SC,reference_ray(6),x(6)
    type(three_d_info),intent(INOUT) ::  v
    TYPE(INTEGRATION_NODE),POINTER:: mag_in,mag_out

    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif

    x=V%X
    if(abs(x(1))+abs(x(3))>absolute_aperture.or.abs(x(6))>t_aperture) then
       messageLOST="exceed absolute_aperture in TRACKV_NODE_SINGLE"
       lost_node=>t
       lost_fibre=>t%parent_fibre
       xlost=x
       CHECK_STABLE=.false.
    endif


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


  END SUBROUTINE TRACK_NODE_SINGLEV

  SUBROUTINE TRACK_NODE_SINGLE_quar(T,P,K) !!
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    type(probe), intent(inout) :: p
    TYPE(INTEGRATION_NODE), pointer, INTENT(INOUT):: T
     TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    LOGICAL TA 
    type(work) w,we
    IF(.NOT.CHECK_STABLE) return
 
    !       CALL RESET_APERTURE_FLAG
    !    endif

    if(abs(p%x(1))+abs(p%x(3))>absolute_aperture.or.abs(p%x(6))>t_aperture) then   !.or.(.not.CHECK_MADX_APERTURE)) then
       messageLOST="exceed absolute_aperture in TRACKR_NODE_SINGLE"
       lost_node=>t
       lost_fibre=>t%parent_fibre
       xlost=p%x
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

    if(k%stochastic) call kick_stochastic_before(t,p)

    SELECT CASE(T%CAS)
    CASE(CASEP1)

      CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,p%x,K)
      if(associated(T%PARENT_FIBRE%MAG%p%aperture)) then
        TA=T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos==-1 .OR. &
           T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos== 0
        if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,p%X)
      endif
      
      global_e= p%x(5)*el%p%p0c
      
    CASE(CASEP2)
      if(associated(T%PARENT_FIBRE%MAG%p%aperture)) then
        TA=T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos==1 .OR. &
           T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos==0
        if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,p%X)
      endif

      CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,p%x,K)
      global_e= p%x(5)*el%p%p0c
    
    CASE(CASE1,CASE2)
  !     el=>T%PARENT_FIBRE%MAG
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE.and.t%cas==case2) &
            call check_S_APERTURE_out(el%p,t%POS_IN_FIBRE-2,p%x)
 
       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=p%X,k=k,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,p%x,k,1)
             CALL FRINGECAV(EL%C4,p%x,k=k,J=1)
          ELSE
             CALL FRINGECAV(EL%C4,p%X,k=k,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,p%x,k,2)
          ENDIF
       case(KINDhel)
          IF(T%CAS==CASE2) THEN
            call fringe_hel(el%he22,p%x,2)
            call fake_shift(el%he22,p%x)
           else
            call fringe_hel(el%he22,p%x,1)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=p%X,k=k,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=p%X,k=k,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=p%X,k=k,J=T%CAS)
       case(KIND10)
          CALL FRINGE_teapot(EL%TP10,p%X,k=k,j=T%CAS)
       case(KIND16,KIND20)
          CALL fringe_STREX(EL%K16,p%x,k,T%CAS)
       case(KIND17)
          STOP 317
       case(KIND21)
          CALL FRINGE_CAV_TRAV(EL%CAV21,X=p%X,k=k,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,p%x,k,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)

          IF(T%CAS==CASE1) THEN
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,1)
          else
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,2)
          endif
          ELSE
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,2)
          else
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,1)
          endif
          CALL ADJUST_WI(EL%WI,p%x,k,T%CAS) 
          ENDIF

  ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,p%x,k,T%CAS)
       case(KINDabell)
          CALL ADJUST_abell(EL%ab,p%x,k,T%CAS)
       case(kindsuperdrift)
        if(el%p%dir==1.and.t%cas==case1) call  PATCH_drift(el%sdr,p%x,k,el%p%exact,1)
        if(el%p%dir==-1.and.t%cas==case2) call  PATCH_drift(el%sdr,p%x,k,el%p%exact,-1)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT
        global_e= p%x(5)*el%p%p0c
    CASE(CASE0)
 !      el=>T%PARENT_FIBRE%MAG
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE)  &
            call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2,p%x)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
 
          call BBKICK(t%bb,p%x,EL%p%BETA0,EL%P%EXACT,k%time)
 

          if(t%bb%patch)call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
 
       SELECT CASE(EL%KIND)
       CASE(KIND0)
         global_e= p%x(5)*el%p%p0c
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND2)
!          CALL TRACK_SLICE(EL%K2,p%x,K,t%POS_IN_FIBRE-2)
            CALL TRACK_SLICE_dkd2(p,k,T)
 
         global_e= p%x(5)*el%p%p0c
       case(KIND3)
          CALL TRACK(EL%K3,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND4)
          CALL TRACK_SLICE_CAV4(p,K,t)
 
    !      CALL TRACK_SLICE(EL%C4,p%x,K,t%POS_IN_FIBRE-2)
          global_e= p%x(5)*el%p%p0c
       case(KIND5)
         ! CALL TRACK_SLICE(EL%S5,p%x,K)
            CALL TRACK_SLICE_sol5(p,k,T)
           global_e= p%x(5)*el%p%p0c
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,p%x,K)
          global_e= p%x(5)*el%p%p0c
       case(KIND7)
      !    CALL TRACK_SLICE(EL%T7,p%x,K,t%POS_IN_FIBRE-2)
          CALL TRACK_SLICE_TKTF(p,K,t,t%POS_IN_FIBRE-2)
           global_e= p%x(5)*el%p%p0c
       case(KIND8)
          CALL TRACK(EL%S8,p%x,K)
          global_e= p%x(5)*el%p%p0c
       case(KIND9)
          CALL TRACK(EL%S9,p%x,K)
          global_e= p%x(5)*el%p%p0c
       case(KIND10)
 !         CALL TRACK_SLICE_TEAPOT_OLD(EL%TP10,p%x,K,t%POS_IN_FIBRE-2)
          CALL TRACK_SLICE_TEAPOT(p,K,t)  !,t%POS_IN_FIBRE-2)
           if(.not.el%electric)  global_e= p%x(5)*el%p%p0c
       case(KIND11:KIND14)
          CALL MONTI(EL%MON14,p%x,k,t%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
          !          CALL TRACK_SLICE(EL%MON14,p%x,K)
       case(KIND15)
          call SEPTTRACK(EL%SEP15,p%x,k,t%POS_IN_FIBRE-2)
         !   global_e= p%x(5)*el%p%p0c done inside
          !          CALL TRACK_SLICE(EL%SEP15,p%x,K)
       case(KIND16,KIND20)
            CALL TRACK_SLICE_strex(p,k,T,t%POS_IN_FIBRE-2)
      !    CALL TRACK_SLICE(EL%K16,p%x,K,t%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
       case(KIND17)
          STOP 317
       case(KIND18)
          call RCOLLIMATORI(EL%RCOL18,p%x,k,t%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
       case(KIND19)
          CALL ECOLLIMATORI(EL%ECOL19,p%x,k,t%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
          !          CALL TRACK_SLICE(EL%ECOL19,p%x,K)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,p%x,k,t%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
       case(KIND22)
          CALL TRACK_SLICE(EL%he22,p%x,k,t%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,p%x,k,t%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,p%x,k,T%POS_IN_FIBRE-2)
       global_e= p%x(5)*el%p%p0c
       case(KINDabell)
          CALL TRACK_SLICE(EL%ab,p%x,k,T%POS_IN_FIBRE-2)
 !      global_e= p%x(5)*el%p%p0c treat like electric
       case(kindsuperdrift)
          call track_slice(EL%sdr,p%x,k)
       global_e= p%x(5)*el%p%p0c
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT
       if(associated(T%PARENT_FIBRE%MAG%p%aperture).and.aperture_all_case0) &
            call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,p%X)


    case(CASET)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,p%x,EL%p%BETA0,EL%P%EXACT,k%time)
          if(t%bb%patch)call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)
       endif
!       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)
    case(CASETF1,CASETF2)

!       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)


    END SELECT
    ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S

     if(k%stochastic) call kick_stochastic_after(t,p)

    !    T%PARENT_FIBRE%MAG=DEFAULT
    if(wherelost==2.and.(.not.check_stable)) then
       t%lost=t%lost+1
    endif
  END SUBROUTINE TRACK_NODE_SINGLE_quaR


!kick_stochastic_after(c,p)

  SUBROUTINE TRACK_NODE_SINGLE_quaP(T,P,K) !!
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), pointer, INTENT(INOUT):: T
    TYPE(PROBE_8),INTENT(INOUT):: P
!    TYPE(REAL_8)  X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(elementp),pointer :: el
    logical(lp) BN2,L
    logical(lp) CHECK_KNOB
    integer(2), pointer,dimension(:)::AN,BN
     logical TA
    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif
!    CALL alloc(X)
!    X=P%X
 
    if(abs(p%x(1))+abs(p%x(3))>absolute_aperture.or.abs(p%x(6))>t_aperture) then
       messageLOST="exceed absolute_aperture in TRACKP_NODE_SINGLE"
       lost_node=>t
       lost_fibre=>t%parent_fibre
       xlost=p%x
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
       CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,p%x,K)
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture)) then
          TA=T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==-1 .OR.  &
             T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==0
          if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,p%x)
       endif
          global_e= p%x(5)*el%p%p0c
    CASE(CASEP2)
    
  
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture)) then
          TA=T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==1 .OR.  &
                 T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==0
          if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,p%X)
       endif

       CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,p%x,K)
       global_e= p%x(5)*el%p%p0c

  
    CASE(CASE1,CASE2)
!       el=>T%PARENT_FIBRE%MAGP
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE.and.t%cas==case2) &
            call check_S_APERTURE_out(el%p,t%POS_IN_FIBRE-2,p%x)


       SELECT CASE(EL%KIND)
 
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=p%X,k=k,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,p%x,k,1)
             CALL FRINGECAV(EL%C4,p%x,k=k,J=1)
          ELSE
             CALL FRINGECAV(EL%C4,p%x,k=k,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,p%x,k,2)
          ENDIF
       case(KINDhel)
          IF(T%CAS==CASE2) THEN
            call fake_shift(el%he22,p%x)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=p%X,k=k,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=p%X,k=k,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=p%X,k=k,J=T%CAS)
       case(KIND10)
          CALL FRINGE_teapot(EL%TP10,p%x,k,T%CAS)
       case(KIND16,KIND20)
          CALL fringe_STREx(EL%K16,p%x,k,T%CAS)
       case(KIND17)
          STOP 317
       case(KIND21)
          CALL FRINGE_CAV_TRAV(EL%CAV21,X=p%X,k=k,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,p%x,k,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)

          IF(T%CAS==CASE1) THEN
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,1)
          else
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,2)
          endif
          ELSE
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,2)
          else
            call ADJUST_LIKE_ABELL(EL%wi,p%x,k,1)
          endif
          CALL ADJUST_WI(EL%WI,p%x,k,T%CAS) 
          ENDIF

       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,p%x,k,T%CAS)    
       case(KINDabell)
          CALL ADJUST_ABELL(EL%AB,p%x,k,T%CAS)
 !      global_e= p%x(5)*el%p%p0c treat like electric
       case(kindsuperdrift)
        if(el%p%dir==1.and.t%cas==case1) call  PATCH_drift(el%sdr,p%x,k,el%p%exact,1)
        if(el%p%dir==-1.and.t%cas==case2) call  PATCH_drift(el%sdr,p%x,k,el%p%exact,-1)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT
        global_e= p%x(5)*el%p%p0c
    CASE(CASE0)

 !      el=>T%PARENT_FIBRE%MAGP
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE) &
            call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2,p%x)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,p%x,EL%p%BETA0,EL%P%EXACT,k%time)
          if(t%bb%patch)call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
       SELECT CASE(EL%KIND)
       CASE(KIND0)
         global_e= p%x(5)*el%p%p0c
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND2)
!          CALL TRACK_SLICE(EL%K2,p%x,K,t%POS_IN_FIBRE-2)
            CALL TRACK_SLICE_dkd2(p,k,T)
          global_e= p%x(5)*el%p%p0c
       case(KIND3)
          CALL TRACK(EL%K3,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND4)
          CALL TRACK_SLICE_CAV4(p,K,t)
     !      CALL TRACK_SLICE(EL%C4,p%x,K,t%POS_IN_FIBRE-2)
          global_e= p%x(5)*el%p%p0c
       case(KIND5)
 
            CALL TRACK_SLICE_sol5(p,k,T)
 
         global_e= p%x(5)*el%p%p0c
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND7)
          IF((EL%T7%BN(2)%KIND==3.OR.EL%T7%L%KIND==3).AND.KNOB) THEN
             CALL GETMAT7(EL%T7)                                      ! RECOMPUTES ONLY IF KNOB (SPEED)
          ENDIF
!          CALL TRACK_SLICE(EL%T7,p%x,K,t%POS_IN_FIBRE-2)
          CALL TRACK_SLICE_TKTF(p,K,t,t%POS_IN_FIBRE-2)
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
         global_e= p%x(5)*el%p%p0c
       case(KIND8)
          CALL TRACK(EL%S8,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND9)
          CALL TRACK(EL%S9,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND10)
          CALL MAKEPOTKNOB(EL%TP10,CHECK_KNOB,AN,BN)
!          CALL TRACK_SLICE_TEAPOT_OLD(EL%TP10,p%x,K,t%POS_IN_FIBRE-2)
          CALL TRACK_SLICE_TEAPOT(p,K,t)  !,t%POS_IN_FIBRE-2)
           CALL UNMAKEPOTKNOB(EL%TP10,CHECK_KNOB,AN,BN)
          if(.not.el%electric)  global_e= p%x(5)*el%p%p0c
       case(KIND11:KIND14)
          CALL MONTI(EL%MON14,p%x,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%MON14,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND15)
          call SEPTTRACK(EL%SEP15,p%x,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%SEP15,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND16,KIND20)
            CALL TRACK_SLICE_strex(p,k,T,t%POS_IN_FIBRE-2)
        !   CALL TRACK_SLICE(EL%K16,p%x,K,t%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
       case(KIND17)
          STOP 317
       case(KIND18)
          call RCOLLIMATORI(EL%RCOL18,p%x,k,t%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
          !          CALL TRACK_SLICE(EL%RCOL18,p%x,K)
       case(KIND19)
          CALL ECOLLIMATORI(EL%ECOL19,p%x,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%ECOL19,p%x,K)
         global_e= p%x(5)*el%p%p0c
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,p%x,k,t%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,p%x,k,t%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
       case(KIND22)
          CALL TRACK_SLICE(EL%he22,p%x,k,t%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,p%x,k,T%POS_IN_FIBRE-2)
         global_e= p%x(5)*el%p%p0c
       case(KINDabell)
          CALL TRACK_SLICE(EL%ab,p%x,k,T%POS_IN_FIBRE-2)
       case(kindsuperdrift)
          call track_slice(EL%sdr,p%x,k)
         global_e= p%x(5)*el%p%p0c
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture).and.aperture_all_case0) &
            call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,p%X)

    case(CASET)

       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,p%x,EL%p%BETA0,EL%P%EXACT,k%time)
          if(t%bb%patch)call PATCH_BB(t%bb,p%x,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)
       endif
 !      IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)
    case(CASETF1,CASETF2)

 !      IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)



    END SELECT
    ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S
    !    T%PARENT_FIBRE%MAGP=DEFAULT
    ! KNOB IS RETURNED TO THE PTC DEFAULT
    ! NEW STUFF WITH KIND=3
    KNOB=.FALSE.
    ! END NEW STUFF WITH KIND=3
 
 
  END SUBROUTINE TRACK_NODE_SINGLE_quaP


  SUBROUTINE TRACK_NODE_SINGLER(T,X,K) !!
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT):: T
    REAL(DP),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    LOGICAL TA
    type(work) w,we
    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif

    if(abs(x(1))+abs(x(3))>absolute_aperture.or.abs(x(6))>t_aperture) then   !.or.(.not.CHECK_MADX_APERTURE)) then
       messageLOST="exceed absolute_aperture in TRACKR_NODE_SINGLE"
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
      if(associated(T%PARENT_FIBRE%MAG%p%aperture)) then
        TA=T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos==-1 .OR. &
           T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos== 0
        if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
      endif
      
      global_e= x(5)*el%p%p0c
      
    CASE(CASEP2)
      if(associated(T%PARENT_FIBRE%MAG%p%aperture)) then
        TA=T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos==1 .OR. &
           T%PARENT_FIBRE%MAG%p%dir*T%PARENT_FIBRE%MAG%p%aperture%pos==0
        if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
      endif

      CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)
      global_e= x(5)*el%p%p0c
    
    CASE(CASE1,CASE2)
  !     el=>T%PARENT_FIBRE%MAG
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE.and.t%cas==case2) &
            call check_S_APERTURE_out(el%p,t%POS_IN_FIBRE-2,x)
 
       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
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
       case(KINDhel)
          IF(T%CAS==CASE2) THEN
            call fringe_hel(el%he22,x,2)
            call fake_shift(el%he22,x)
           else
            call fringe_hel(el%he22,x,1)
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

          IF(T%CAS==CASE1) THEN
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,X,k,1)
          else
            call ADJUST_LIKE_ABELL(EL%wi,X,k,2)
          endif
          ELSE
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,X,k,2)
          else
            call ADJUST_LIKE_ABELL(EL%wi,X,k,1)
          endif
          CALL ADJUST_WI(EL%WI,X,k,T%CAS) 
          ENDIF

  ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,k,T%CAS)
       case(KINDabell)
          CALL ADJUST_abell(EL%ab,X,k,T%CAS)
       case(kindsuperdrift)
        if(el%p%dir==1.and.t%cas==case1) call  PATCH_drift(el%sdr,X,k,el%p%exact,1)
        if(el%p%dir==-1.and.t%cas==case2) call  PATCH_drift(el%sdr,X,k,el%p%exact,-1)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT
        global_e= x(5)*el%p%p0c
    CASE(CASE0)
 !      el=>T%PARENT_FIBRE%MAG
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE)  &
            call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2,x)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
 
          call BBKICK(t%bb,X,EL%p%BETA0,EL%P%EXACT,k%time)
 

          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
 
       SELECT CASE(EL%KIND)
       CASE(KIND0)
         global_e= x(5)*el%p%p0c
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND2)
          CALL TRACK_SLICE_dkd2_OLD(EL%K2,X,K,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND3)
          CALL TRACK(EL%K3,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND4)
          CALL TRACK_SLICE_CAV4_OLD(EL%C4,X,K,t%POS_IN_FIBRE-2)
          global_e= x(5)*el%p%p0c
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X,K)
          global_e= x(5)*el%p%p0c
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X,K)
          global_e= x(5)*el%p%p0c
       case(KIND7)
          CALL TRACK_SLICE(EL%T7,X,K,t%POS_IN_FIBRE-2)
          global_e= x(5)*el%p%p0c
       case(KIND8)
          CALL TRACK(EL%S8,X,K)
          global_e= x(5)*el%p%p0c
       case(KIND9)
          CALL TRACK(EL%S9,X,K)
          global_e= x(5)*el%p%p0c
       case(KIND10)
          CALL TRACK_SLICE_TEAPOT_OLD(EL%TP10,X,K,t%POS_IN_FIBRE-2)
          if(.not.el%electric)  global_e= x(5)*el%p%p0c
       case(KIND11:KIND14)
          CALL MONTI(EL%MON14,X,k,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
          !          CALL TRACK_SLICE(EL%MON14,X,K)
       case(KIND15)
          call SEPTTRACK(EL%SEP15,X,k,t%POS_IN_FIBRE-2)
         !   global_e= x(5)*el%p%p0c done inside
          !          CALL TRACK_SLICE(EL%SEP15,X,K)
       case(KIND16,KIND20)
          CALL TRACK_SLICE_STREX_OLD(EL%K16,X,K,t%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
       case(KIND17)
          STOP 317
       case(KIND18)
          call RCOLLIMATORI(EL%RCOL18,X,k,t%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
       case(KIND19)
          CALL ECOLLIMATORI(EL%ECOL19,X,k,t%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
          !          CALL TRACK_SLICE(EL%ECOL19,X,K)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,k,t%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
       case(KIND22)
          CALL TRACK_SLICE(EL%he22,X,k,t%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,k,t%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,k,T%POS_IN_FIBRE-2)
       global_e= x(5)*el%p%p0c
       case(KINDabell)
          CALL TRACK_SLICE(EL%ab,X,k,T%POS_IN_FIBRE-2)
 !      global_e= x(5)*el%p%p0c treat like electric
       case(kindsuperdrift)
          call track_slice(EL%sdr,X,k)
       global_e= x(5)*el%p%p0c
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT
       if(associated(T%PARENT_FIBRE%MAG%p%aperture).and.aperture_all_case0) &
            call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)


    case(CASET)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X,EL%p%BETA0,EL%P%EXACT,k%time)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)
       endif
!       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)
    case(CASETF1,CASETF2)

!       IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)


    END SELECT
    ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S


    !    T%PARENT_FIBRE%MAG=DEFAULT
    if(wherelost==2.and.(.not.check_stable)) then
       t%lost=t%lost+1
    endif
  END SUBROUTINE TRACK_NODE_SINGLER


  SUBROUTINE TRACK_NODE_SINGLEP(T,X,K) !!
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
     logical TA
    IF(.NOT.CHECK_STABLE) return
    !       CALL RESET_APERTURE_FLAG
    !    endif
    
    if(abs(x(1))+abs(x(3))>absolute_aperture.or.abs(x(6))>t_aperture) then
       messageLOST="exceed absolute_aperture in TRACKP_NODE_SINGLE"
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
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture)) then
          TA=T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==-1 .OR.  &
             T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==0
          if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,X)
       endif
          global_e= x(5)*el%p%p0c
    CASE(CASEP2)
    
  
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture)) then
          TA=T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==1 .OR.  &
                 T%PARENT_FIBRE%MAGP%p%dir*T%PARENT_FIBRE%MAGP%p%aperture%pos==0
          if(TA) call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,X)
       endif

       CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)
       global_e= x(5)*el%p%p0c

  
    CASE(CASE1,CASE2)
!       el=>T%PARENT_FIBRE%MAGP
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE.and.t%cas==case2) &
            call check_S_APERTURE_out(el%p,t%POS_IN_FIBRE-2,x)


       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
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
       case(KINDhel)
          IF(T%CAS==CASE2) THEN
            call fake_shift(el%he22,x)
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

          IF(T%CAS==CASE1) THEN
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,X,k,1)
          else
            call ADJUST_LIKE_ABELL(EL%wi,X,k,2)
          endif
          ELSE
          if(el%p%dir==1) then
            call ADJUST_LIKE_ABELL(EL%wi,X,k,2)
          else
            call ADJUST_LIKE_ABELL(EL%wi,X,k,1)
          endif
          CALL ADJUST_WI(EL%WI,X,k,T%CAS) 
          ENDIF

       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,k,T%CAS)    
       case(KINDabell)
          CALL ADJUST_ABELL(EL%AB,X,k,T%CAS)
 !      global_e= x(5)*el%p%p0c treat like electric
       case(kindsuperdrift)
        if(el%p%dir==1.and.t%cas==case1) call  PATCH_drift(el%sdr,X,k,el%p%exact,1)
        if(el%p%dir==-1.and.t%cas==case2) call  PATCH_drift(el%sdr,X,k,el%p%exact,-1)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT
        global_e= x(5)*el%p%p0c
    CASE(CASE0)

 !      el=>T%PARENT_FIBRE%MAGP
       if(s_aperture_CHECK.and.associated(el%p%A).AND.CHECK_MADX_APERTURE) &
            call check_S_APERTURE(el%p,t%POS_IN_FIBRE-2,x)
       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X,EL%p%BETA0,EL%P%EXACT,k%time)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
       SELECT CASE(EL%KIND)
       CASE(KIND0)
         global_e= x(5)*el%p%p0c
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND2)
          CALL TRACK_SLICE_dkd2_OLD(EL%K2,X,K,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND3)
          CALL TRACK(EL%K3,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND4)
          CALL TRACK_SLICE_CAV4_OLD(EL%C4,X,K,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X,K)
         global_e= x(5)*el%p%p0c
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
         global_e= x(5)*el%p%p0c
       case(KIND8)
          CALL TRACK(EL%S8,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND9)
          CALL TRACK(EL%S9,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND10)
          CALL MAKEPOTKNOB(EL%TP10,CHECK_KNOB,AN,BN)
          CALL TRACK_SLICE_TEAPOT_OLD(EL%TP10,X,K,t%POS_IN_FIBRE-2)
          CALL UNMAKEPOTKNOB(EL%TP10,CHECK_KNOB,AN,BN)
          if(.not.el%electric)  global_e= x(5)*el%p%p0c
       case(KIND11:KIND14)
          CALL MONTI(EL%MON14,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%MON14,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND15)
          call SEPTTRACK(EL%SEP15,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%SEP15,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND16,KIND20)
          CALL TRACK_SLICE_STREX_OLD(EL%K16,X,K,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND17)
          STOP 317
       case(KIND18)
          call RCOLLIMATORI(EL%RCOL18,X,k,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
          !          CALL TRACK_SLICE(EL%RCOL18,X,K)
       case(KIND19)
          CALL ECOLLIMATORI(EL%ECOL19,X,k,t%POS_IN_FIBRE-2)
          !          CALL TRACK_SLICE(EL%ECOL19,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,k,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,k,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND22)
          CALL TRACK_SLICE(EL%he22,X,k,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,k,T%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KINDabell)
          CALL TRACK_SLICE(EL%ab,X,k,T%POS_IN_FIBRE-2)
       case(kindsuperdrift)
          call track_slice(EL%sdr,X,k)
         global_e= x(5)*el%p%p0c
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture).and.aperture_all_case0) &
            call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,X)

    case(CASET)

       if(associated(t%bb).and.dobb.and.do_beam_beam) then

          if(t%bb%patch) call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_true)
          call BBKICK(t%bb,X,EL%p%BETA0,EL%P%EXACT,k%time)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)
       endif
 !      IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)
    case(CASETF1,CASETF2)

 !      IF(ASSOCIATED(T%T)) CALL TRACK(T%T,X)



    END SELECT
    ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S
    !    T%PARENT_FIBRE%MAGP=DEFAULT
    ! KNOB IS RETURNED TO THE PTC DEFAULT
    ! NEW STUFF WITH KIND=3
    KNOB=.FALSE.
    ! END NEW STUFF WITH KIND=3

  END SUBROUTINE TRACK_NODE_SINGLEP








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
    INTEGER I,J,k,TEAPOT_LIKE,j0,dir
    REAL(DP) S,DLD,DL,LI,SL,u0(2),ub(2),uj(2),n0,nj,r0
    LOGICAL(LP) CIRCULAR,doit
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

doit=p%mag%kind==kind16.and.p%mag%p%b0/=0.0_dp
   if(doit) then
    j0=0
    dir=1
    if(p%dir==-1) dir=2
     r0=1.0_dp/p%mag%p%b0
    u0(1)=-sin(p%mag%p%edge(dir))*r0
    u0(2)=cos(p%mag%p%edge(dir))*r0
    ub=u0
   endif

       DO J=1,P%MAG%P%NST
  if(doit) then
        j0=j0+1
        LI=(j0*p%mag%L)/(p%mag%p%nst)
        uj=[li,0.0_dp]
        uj=uj+u0
         
        uj(2)=uj(2)+sqrt(r0**2-uj(1)**2)-u0(2)
        n0=sqrt(ub(1)**2+ub(2)**2)
        nj=sqrt(uj(1)**2+uj(2)**2)
        DLD=acos( ((ub(1)*uj(1)+ub(2)*uj(2)) /n0/nj)   )*r0
         ub=uj
 endif
!  2017.4.27
          S=S+DLD
          LI=LI+DL
          SL=SL+P%DIR*DL
          CALL APPEND_EMPTY_THIN( L )
          L%END%TEAPOT_LIKE=TEAPOT_LIKE
          L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=DL;L%END%S(5)=DLD;L%END%ds_ac=DLD;
          L%END%CAS=CASE0
          L%END%pos_in_fibre=J+2
          L%END%pos=k;k=k+1;
          L%END%PARENT_NODE_LAYOUT=>L
          L%END%PARENT_FIBRE=>P


!!!!!!!!!
!          S=S+DLD
!          LI=LI+DL
!          SL=SL+P%DIR*DL
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
   call allocate_node_frame( R )
 
 

  END SUBROUTINE MAKE_NODE_LAYOUT_2


  SUBROUTINE stat_NODE_LAYOUT( L )
    implicit none
    TYPE (NODE_LAYOUT), pointer :: L

  if(lielib_print(4)==1) then
    WRITE(6,*)  " PARENT LAYOUT NAME :", L%PARENT_LAYOUT%NAME(1:len_trim(L%PARENT_LAYOUT%NAME))
    WRITE(6,*) " NUMBER OF ORIGINAL LAYOUT ELEMENTS :", L%PARENT_LAYOUT%N
    WRITE(6,*) " NUMBER OF THIN OBJECTS :", L%N
    WRITE(6,*) " TOTAL IDEAL LENGTH OF STRUCTURE :", L%END%S(1)
    WRITE(6,*) " TOTAL INTEGRATION LENGTH OF STRUCTURE (mad8 style survey) :", L%END%S(3)
  endif
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
        call SURVEY(my_ering)
    ENDIF
    IF(.NOT.ASSOCIATED(my_ering%T%start%B)) THEN
        call SURVEY(my_ering)
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

    subroutine convert_bmad_to_ptcar(z,b1,time)
    IMPLICIT NONE
    real(dp),target,intent(INOUT) ::  z(6)
    real(dp) b0,t,b1
     logical(lp)  time

    if(time) then
     b0=b1
    else
     b0=1
    endif 

     if(tangent) then
      t=sqrt(1.0_dp+2.0_dp*z(5)/b0+z(5)**2)/sqrt(1.0_dp+z(2)**2+z(4)**2) 
      z(2)=z(2)*t
      z(4)=z(4)*t
     else
     t=z(6)
     z(6)=-z(5)*sqrt(1.0_dp/b0**2+2*t+t**2)/(1.0_dp+t)
     z(5)=sqrt(1.0_dp/b0**2+2*t+t**2)-1.0_dp/b0
     endif
     end subroutine convert_bmad_to_ptcar

    subroutine convert_bmad_to_ptcap(z,b1,time)
    IMPLICIT NONE
    type(real_8),target,intent(INOUT) ::  z(6)
    type(real_8) t
    real(dp) b0,b1
     logical(lp)  time

    if(time) then
     b0=b1
    else
     b0=1
    endif 

     call alloc(t)

     if(tangent) then
      t=sqrt(1.0_dp+2.0_dp*z(5)/b0+z(5)**2)/sqrt(1.0_dp+z(2)**2+z(4)**2) 
      z(2)=z(2)*t
      z(4)=z(4)*t
     else
     t=z(6)

     z(6)=-z(5)*sqrt(1.0_dp/b0**2+2*t+t**2)/(1.0_dp+t)
     z(5)=sqrt(1.0_dp/b0**2+2*t+t**2)-1.0_dp/b0
     endif
     call kill(t)

     end subroutine convert_bmad_to_ptcap 

    subroutine convert_ptc_to_bmadar(z,b1,time,LD)
    IMPLICIT NONE
    real(dp),target,intent(INOUT) :: z(6)
    real(dp), optional :: LD
    real(dp) b0,t,b1,l
     logical(lp)  time
    l=0
    if(present(ld)) l=ld
    if(time) then
     b0=b1 
     l=l/b1
    else
     b0=1
    endif 

     if(tangent) then
      t=1.0_dp/sqrt(1.0_dp+2.0_dp*z(5)/b0+z(5)**2-z(2)**2-z(4)**2) 
      z(2)=z(2)*t
      z(4)=z(4)*t
     else
      t=z(5)
      z(5)=-(z(6)-l)*sqrt(1.d0 +2*t/b0+t**2)/(1.d0/b0+t)
      z(6)=sqrt(1.0_dp +2*t/b0+t**2)-1.d0 
     endif
     end subroutine convert_ptc_to_bmadar   


    subroutine convert_ptc_to_bmadap(z,b1,time,ld)
    IMPLICIT NONE
    type(real_8),target,intent(INOUT) ::  z(6)
    type(real_8) t 
     real(dp), optional :: LD
     real(dp) b0,b1,l
     logical(lp)  time

     l=0
    if(present(ld)) l=ld
    if(time) then
     b0=b1
     l=l/b1
    else
     b0=1
    endif 

    
     call alloc(t)
     if(tangent) then
      t=1.0_dp/sqrt(1.0_dp+2.0_dp*z(5)/b0+z(5)**2-z(2)**2-z(4)**2) 
      z(2)=z(2)*t
      z(4)=z(4)*t
     else
      t=z(5)
      z(5)=-(z(6)-l)*sqrt(1.0_dp+2*t/b0+t**2)/(1.0_dp/b0+t)
      z(6)=sqrt(1.0_dp +2*t/b0+t**2)-1.0_dp 
     endif
     call kill(t)

     end subroutine convert_ptc_to_bmadap  


    subroutine convert_bmad_to_ptcr(z,b1,time)
    IMPLICIT NONE
    type(probe),target,intent(INOUT) ::  z
    real(dp) b0,t,b1
     logical(lp) time

    if(time) then
     b0=b1
    else
     b0=1
    endif 
     if(tangent) then
      t=sqrt(1.0_dp+2.0_dp*z%x(5)/b0+z%x(5)**2)/sqrt(1.0_dp+z%x(2)**2+z%x(4)**2) 
      z%x(2)=z%x(2)*t
      z%x(4)=z%x(4)*t
     else
      t=z%x(6)
      z%x(6)=-z%x(5)*sqrt(1.0_dp/b0**2+2*t+t**2)/(1.0_dp+t)
      z%x(5)=sqrt(1.0_dp/b0**2+2*t+t**2)-1.0_dp/b0
     endif
     end subroutine convert_bmad_to_ptcr   


    subroutine convert_bmad_to_ptcp(z,b1,time)
    IMPLICIT NONE
    type(probe_8),target,intent(INOUT) ::  z
    type(real_8) t
    real(dp) b0,b1
    logical(lp)  time
    if(time) then
     b0=b1
    else
     b0=1
    endif 
     call alloc(t)

     if(tangent) then
      t=sqrt(1.0_dp+2.0_dp*z%x(5)/b0+z%x(5)**2)/sqrt(1.0_dp+z%x(2)**2+z%x(4)**2) 
      z%x(2)=z%x(2)*t
      z%x(4)=z%x(4)*t
     else
     t=z%x(6)

     z%x(6)=-z%x(5)*sqrt(1.d0/b0**2+2*t+t**2)/(1.0_dp+t)
     z%x(5)=sqrt(1.0_dp/b0**2+2*t+t**2)-1.0_dp/b0
    endif 

     call kill(t)
     call convert_bmad_to_ptc(z%x0,b1,time)
     end subroutine convert_bmad_to_ptcp   

    subroutine convert_ptc_to_bmadr(z,b1,time,LD)
    IMPLICIT NONE
    type(probe),target,intent(INOUT) :: z
    real(dp) b0,t,b1,l
    logical(lp)  time
     real(dp), optional :: LD

     l=0
    if(present(ld)) l=ld
    if(time) then
     b0=b1
     l=l/b1
    else
     b0=1
    endif 
     if(tangent) then
      t=1.0_dp/sqrt(1.0_dp+2.0_dp*z%x(5)/b0+z%x(5)**2-z%x(2)**2-z%x(4)**2) 
      z%x(2)=z%x(2)*t
      z%x(4)=z%x(4)*t
     else
     t=z%x(5)
      z%x(5)=-(z%x(6)-l)*sqrt(1.0_dp +2*t/b0+t**2)/(1.0_dp/b0+t)
      z%x(6)=sqrt(1.0_dp+2*t/b0+t**2)-1.0_dp 
     endif
     end subroutine convert_ptc_to_bmadr   


   subroutine convert_ptc_to_bmadp(z,b1,time,LD)
    IMPLICIT NONE
    type(probe_8),target,intent(INOUT) ::  z
    type(real_8) t
    real(dp) b0,b1,l
    logical(lp)  time
     real(dp), optional :: LD

     l=0
    if(present(ld)) l=ld
    if(time) then
     b0=b1
     l=l/b1
    else
     b0=1
    endif 
     call alloc(t)

     if(tangent) then
      t=1.0_dp/sqrt(1.0_dp+2.0_dp*z%x(5)/b0+z%x(5)**2-z%x(2)**2-z%x(4)**2) 
      z%x(2)=z%x(2)*t
      z%x(4)=z%x(4)*t
     else
      t=z%x(5)
      z%x(5)=-(z%x(6)-l)*sqrt(1.0_dp +2*t/b0+t**2)/(1.0_dp/b0+t)
      z%x(6)=sqrt(1.0_dp +2*t/b0+t**2)-1.e0_dp 
     endif
     call kill(t)
 
           call convert_ptc_to_bmad(z%x0,b1,time,LD)

     end subroutine convert_ptc_to_bmadp 



     subroutine in_noncanonical_units
     implicit none  
      use_bmad_units=.true.
      tangent=.true.
      ndpt_bmad=0
     end subroutine in_noncanonical_units

     subroutine in_canonical_units
     implicit none  
      use_bmad_units=.false.
      tangent=.false.
      ndpt_bmad=0
     end subroutine in_canonical_units



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   New Survey Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
subroutine SURVEY_EXIST_PLANAR_IJ_new(L,pi,fi,ent,a)
implicit none
type(layout) L
integer :: pi,i
integer ,target, optional:: fi
type(fibre), pointer :: p
type(fibre),pointer :: f
real(dp), optional, intent(in):: a(3),ent(3,3)
 
if(.not.associated(l%t)) then
 call MAKE_NODE_LAYOUT(l)
endif


p=>L%start
do i=1,pi-1
p=>p%next
enddo

if(present(fi)) then
 f=>p
 do i=pi,fi-1
  f=>f%next
 enddo
 call  survey_integration_layout(p,f,ent=ent,a=a)
else
 call  survey_integration_layout(p,ent=ent,a=a)
endif

end subroutine SURVEY_EXIST_PLANAR_IJ_new


subroutine survey_integration_layout(p,f,ent,a)
implicit none
type(fibre), target :: p
type(fibre),target, optional:: f
type(fibre),pointer :: p1,p2
real(dp), optional, intent(in):: a(3),ent(3,3)
real(dp)  a0(3),ent0(3,3)
 
if(present(a)) then
 a0=a
else
 !a0=global_origin
!if(p%dir==1) then
   a0=p%t1%a
 ! a0=p%chart%f%a          !global_origin
!else
!  a0=p%chart%f%b          !global_origin
!   a0=p%t2%b
!endif
endif
if(present(ent)) then
 ent0=ent
else
 !ent0=global_FRAME
!if(p%dir==1) then
 !ent0=p%chart%f%ent !global_FRAME
   ent0=p%t1%ent

!else
 !ent0=p%chart%f%exi !global_FRAME
!   ent0=p%t2%ent

!endif
endif
 
 
call survey_integration_fibre(p,ent0,a0)

 
 
p1=>p%next
if(present(f)) then
 p2=>f
else
 p2=>p
endif
 
do while(.not.associated(p2,p1).and.associated(p1))
!if(p1%previous%dir==1) then 
 call survey_integration_fibre(p1,p1%previous%t2%exi,p1%previous%t2%b)
!else
! call survey_integration_fibre(p1,p1%previous%t2%ent,p1%previous%t2%a)
!endif 
!call survey_integration_fibre(p1,p1%previous%chart%f%b,p1%previous%chart%f%exi)
 
p1=>p1%next
enddo
 
end subroutine survey_integration_layout

  SUBROUTINE SURVEY_FIBRE_new(C,ENT,A,nogirder)   !,MAGNETFRAME
    !changed
    ! SURVEYS A SINGLE ELEMENT FILLS IN CHART AND MAGNET_CHART; LOCATES ORIGIN AT THE ENTRANCE OR EXIT
    IMPLICIT NONE
    TYPE(FIBRE), TARGET , INTENT(INOUT):: C
    !    TYPE(MAGNET_FRAME), OPTIONAL :: MAGNETFRAME
    logical(lp), OPTIONAL :: nogirder
    REAL(DP), INTENT(INOUT)  :: ENT(3,3),A(3)
    REAL(DP) D(3),ANG(3),DT(3)
    LOGICAL(LP) dog

    REAL(DP)dg1(3),ag1(3),ENT0(3,3)
    REAL(DP)dg2(3),ag2(3)
    dog=.true.
    if(present(nogirder)) dog=.not.nogirder     ! IF  DOGIRDER IS SET TO TRUE, MOVES GIRDER FRAME DURING SURVEY
!!!  RECORDS RELATIVE POSITION OF GIRDER
    IF(ASSOCIATED(C%MAG%GIRDER_FRAME).and.dog) THEN
       call FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            C%MAG%GIRDER_FRAME%a, C%MAG%GIRDER_FRAME%ent,dg1,ag1)
       call FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            C%MAG%GIRDER_FRAME%b, C%MAG%GIRDER_FRAME%exi,dg2,ag2)
    ENDIF


     call survey_integration_fibre(C,ent,a)


!!!  PLACES INTO THE COMPUTED RELATIVE POSITION OF GIRDER
    IF(ASSOCIATED(C%MAG%GIRDER_FRAME).and.dog) THEN    ! IF  DOGIRDER IS SET TO TRUE, MOVES GIRDER FRAME DURING SURVEY
       call INVERSE_FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            dg1,ag1,C%MAG%GIRDER_FRAME%a, C%MAG%GIRDER_FRAME%ent)
       call INVERSE_FIND_PATCH(C%chart%F%a,C%chart%F%ENT, &
            dg2,ag2,C%MAG%GIRDER_FRAME%b, C%MAG%GIRDER_FRAME%exi)
    ENDIF

  END SUBROUTINE SURVEY_FIBRE_new



  SUBROUTINE SURVEY_EXIST_PLANAR_L_NEW(PLAN,ENT,A) ! CALLS ABOVE ROUTINE FROM FIBRE #1 TO #PLAN%N : STANDARD SURVEY
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(INOUT):: PLAN
    REAL(DP),OPTIONAL, INTENT(INOUT) :: A(3),ENT(3,3)
  
    CALL survey(PLAN,pi=1,ent=ENT,a=A)
 
  END SUBROUTINE SURVEY_EXIST_PLANAR_L_NEW


 

subroutine survey_integration_fibre(p,exi,b,previous)
implicit none
type(fibre), target :: p
type(integration_node), pointer :: t
integer i,j
type(layout), pointer  :: r
real(dp),optional, intent(in):: b(3),exi(3,3)
!real(dp)  b0(3),exi0(3,3)
real(dp) a0(3),ent0(3,3),ang(3) 
 logical, optional :: previous
logical prev

prev=.true.

if(present(previous)) then
 prev=previous
endif

!a0=b0
!ent0=exi0
r=>p%parent_layout
if(.not.associated(r%t)) then
 call make_node_layout(r)
 call survey(r)
    CALL  allocate_node_frame( R)
! call FILL_SURVEY_DATA_IN_NODE_LAYOUT(r)
endif
if(.not.associated(p%t1%a))     CALL  allocate_node_frame( R)   !call FILL_SURVEY_DATA_IN_NODE_LAYOUT(r)
 
 IF(present(b).and.present(exi)) then
   a0=b
   ent0=exi
ELSE
if(associated(p%previous)) then
 if(prev) then
   ! if(p%previous%dir==1) then 
      ent0=p%previous%t2%exi
      a0=p%previous%t2%b
   ! else
    !  ent0=p%previous%t2%ent
    !  a0=p%previous%t2%a 
   ! endif
 else
  !  if(p%previous%dir==1) then 
      ent0=p%t1%ent
      a0=p%t1%a
 !   else
  !    ent0=p%t1%exi
  !    a0=p%t1%b
  !  endif
 endif
else
 !if(p%dir==1) then 
  ent0=p%t1%ent
  a0=p%t1%a
 !else
 !   ent0=p%t1%exi
 !   a0=p%t1%b
 !endif 
endif
ENDIF
 

t=>p%t1
 j=1
if(p%dir==-1) j=2
! write(6,*) p%mag%name
! write(6,*) a0
call survey_integration_node_p(t,ent0,a0,1)

 
t=>t%next
call survey_integration_fringe(t,ent0,a0,j)
 
do i=1,p%mag%p%nst
 t=>t%next

 call survey_integration_node_case0(t,ent0,a0)
 
enddo
t=>t%next
 j=2
if(p%dir==-1) j=1
call survey_integration_fringe(t,ent0,a0,j)
t=>t%next
call survey_integration_node_p(t,ent0,a0,2)
 
 
  
 

 
!!! entrance chart
 CALL COMPUTE_ENTRANCE_ANGLE(p%chart%f%ent,p%chart%f%exi,ANG)
p%chart%f%mid=p%chart%f%ent
 
ang=ang/2
p%chart%f%o=0.5_dp*(p%chart%f%a+p%chart%f%b)
CALL GEO_ROT(p%chart%f%mid,ANG,1,basis=p%chart%f%ent)
 
CALL COMPUTE_ENTRANCE_ANGLE(p%mag%p%f%ent,p%mag%p%f%exi,ANG)
p%mag%p%f%mid=p%mag%p%f%ent
ang=ang/2
p%mag%p%f%o=0.5_dp*(p%mag%p%f%a+p%mag%p%f%b)
CALL GEO_ROT(p%mag%p%f%mid,ANG,1,basis=p%mag%p%f%ent)

if(.not.p%mag%mis) then
p%magp%p%f%mid=p%mag%p%f%mid
p%magp%p%f%o=p%mag%p%f%o
endif
 


end subroutine survey_integration_fibre

subroutine survey_integration_node_p(t,ent0,a0,k)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
real(dp) pix1(3),pix2(3) ,a0(3),exi0(3,3),ent0(3,3)
logical(lp) :: ENTERING=my_true
 integer k
f=>t%parent_fibre
 

if(k==1) then
f=>t%parent_fibre
pix1=0.0_dp;pix2=0.0_dp;

!if(f%dir==1) then
 t%a=a0
 t%ent=ent0
 exi0=t%ent
!else
! t%b=a0
! t%exi=ent0
! exi0=t%exi
!endif

if(f%patch%A_X1==-1) pix1(1)=pi
if(f%patch%A_X2==-1) pix2(1)=pi

!
call GEO_ROT(exi0,pix1,1, ent0)
call GEO_ROT(exi0,f%patch%a_ang,1, exi0)
call TRANSLATE_point(a0,f%patch%A_D,1,exi0)  
call GEO_ROT(exi0,pix2,1, exi0)

!!!  These are chart frames so depends on dir=1 or -1

!!! entrance chart
if(f%dir==1) then
f%chart%f%ent=exi0
f%chart%f%a=a0
else
f%chart%f%exi=exi0
f%chart%f%b=a0
endif

pix1=0.0_dp
!if(f%dir==1) then
 pix1(3)=f%MAG%P%TILTD
!else
!  pix1(3)=f%MAG%P%TILTD
!endif
 call GEO_ROT(exi0,pix1,1, exi0)


! previous prot location




    IF(f%MAG%MIS) THEN
 
        call MIS_survey(a0,exi0,f,a0,exi0,ENTERING)
 
    ENDIF
 
!!!  These are magnet frames so depends on dir=1 or -1
 if(f%dir==1) then
  f%mag%p%f%ent=exi0
  f%mag%p%f%a=a0
  if(.not.f%mag%mis) then
  f%magp%p%f%ent=exi0
  f%magp%p%f%a=a0
  endif
else
  f%mag%p%f%exi=exi0
  f%mag%p%f%b=a0
  if(.not.f%mag%mis) then
  f%magp%p%f%exi=exi0
  f%magp%p%f%b=a0
 endif
endif

 ! new prot location
 




!if(f%dir==1) then
 t%b=a0
 t%exi=exi0
 ent0=exi0
!else
! t%a=a0
! t%ent=exi0
! ent0=exi0
!endif
 
 
else ! k


f=>t%parent_fibre

 
!if(f%dir==1) then
 t%a=a0
 t%ent=ent0
 exi0=t%ent
!else
! t%b=a0
! t%exi=ent0
! exi0=t%exi
!endif


!!!  These are magnet frames so depends on dir=1 or -1

 if(f%dir==-1) then
  f%mag%p%f%ent=ent0
  f%mag%p%f%a=a0
  if(.not. f%mag%mis) then
  f%magp%p%f%ent=ent0
  f%magp%p%f%a=a0
  endif
else
  f%mag%p%f%exi=ent0
  f%mag%p%f%b=a0
  if(.not. f%mag%mis) then
  f%magp%p%f%exi=ent0
  f%magp%p%f%b=a0
 endif
endif

 
        call MIS_survey(a0,exi0,f,a0,exi0,.not.ENTERING)
 


pix1=0.0_dp
!if(f%dir==1) then
 pix1(3)=-f%MAG%P%TILTD
!else  
!  pix1(3)=-f%MAG%P%TILTD    ! strange   
!endif

 call GEO_ROT(exi0,pix1,1, exi0)

!!!  These are chart frames so depends on dir=1 or -1

if(f%dir==1) then
 f%chart%f%exi=exi0
 f%chart%f%b=a0
else
 f%chart%f%ent=exi0
 f%chart%f%a=a0
endif

 

pix1=0.0_dp;pix2=0.0_dp;
if(f%patch%B_X1==-1) pix1(1)=pi
if(f%patch%B_X2==-1) pix2(1)=pi

!
call GEO_ROT(exi0,pix1,1, ent0)
call GEO_ROT(exi0,f%patch%B_ang,1, exi0)
call TRANSLATE_point(a0,f%patch%b_D,1,exi0)  
call GEO_ROT(exi0,pix2,1, exi0)


!!!!  all missing
t%b=a0
t%exi=exi0
ent0=exi0
 

!if(f%dir==1) then
 t%b=a0
 t%exi=exi0
 ent0=exi0
!else
! t%a=a0
! t%ent=exi0
! ent0=exi0
!endif


!       X(3)=C%PATCH%B_X1*X(3);X(4)=C%PATCH%B_X1*X(4);
!       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAG%P%BETA0,PATCH,k%TIME)
!       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAG%P%BETA0,PATCH,k%TIME)
!       CALL ROT_XY(C%PATCH%B_ANG(3),X)  !,PATCH)
!       CALL TRANS(C%PATCH%B_D,X,C%MAG%P%BETA0,PATCH,k%TIME)
!       X(3)=C%PATCH%B_X2*X(3);X(4)=C%PATCH%B_X2*X(4);


endif

 


end subroutine survey_integration_node_p

subroutine survey_integration_fringe(t,ent0,a0,j)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
real(dp) a0(3),ent0(3,3),pix1(3)
integer j

f=>t%parent_fibre

 !if(f%dir==1) then
 t%a=a0
 t%ent=ent0
!else
! t%b=a0
! t%exi=ent0
!endif




if(associated(t%parent_fibre%mag%sdr)) then
 call survey_integration_special_superdrift(t,ent0,a0)
else


 pix1=0.0_dp
if(f%mag%p%exact.and.(f%mag%kind/=kind10)) then

if(j==1) then
if(f%dir==1) then
 pix1(2)=f%MAG%P%edge(1)
else
 pix1(2)=-f%MAG%P%edge(2)
endif
else
if(f%dir==1) then
 pix1(2)=f%MAG%P%edge(2)
else
 pix1(2)=-f%MAG%P%edge(1)   
endif

endif

 call GEO_ROT(ent0,pix1,1, ent0)
endif
 

 
 



if(f%mag%kind==kindpa) then

call ADJUST_PANCAKE_frame(f%mag%pa,ent0,a0,j)
!
!write(6,*) " I am here in  1"
endif 

if(f%mag%kind==kind4) then


call ADJUST_cav_frame(f%mag%c4,ent0,a0,j,f%dir)
!
!write(6,*) " I am here in   2"
endif

if(f%mag%kind==kindabell) then

call ADJUST_abell_frame(f%mag%ab,ent0,a0,j)
!
!write(6,*) " I am here in   3"
endif  


 
!if(f%dir==1) then
  t%b=a0
  t%exi=ent0
  a0=t%b
  ent0=t%exi
! else
!  t%a=a0
!  t%ent=ent0
!  a0=t%a
!  ent0=t%ent
! endif
 
endif

end subroutine survey_integration_fringe




subroutine survey_integration_node_case0(t,ent0,b0)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
type(element), pointer :: m
type(magnet_chart), pointer :: p
real(dp) h,d(3),ang(3),b0(3),exi0(3,3),ent0(3,3)
integer dir
f=>t%parent_fibre
m=>f%mag
p=>m%p



 dir=f%dir

! if(f%dir==1) then
t%a=b0
t%ent=ent0
exi0=t%ent
!else
! t%b=b0
!  t%exi=ent0
! exi0=t%exi
!endif



select case(m%kind) 

CASE(KIND0,KIND1,KIND3:KIND5,KIND8:KIND9,KIND11:KIND15,KIND17:KIND22,kindwiggler,kindsuperdrift)
   h=dir*p%lc/p%nst
   d=(/0.0_dp,0.0_dp,h/)

   call geo_tra(b0,exi0,d,1)
CASE(KIND2,KIND6:KIND7,KIND10)
   h=dir*p%ld/p%nst
  if(p%b0==0.0_dp) then
   d=(/0.0_dp,0.0_dp,h/)

   call geo_tra(b0,exi0,d,1)
  else
  ang=0.0_dp
  ang(2)=h*p%b0/2
  h=2*sin(ang(2))/p%b0
  d=(/0.0_dp,0.0_dp,h/)
  call geo_rot(exi0,exi0,ang,exi0)
  call geo_tra(b0,exi0,d,1)
  call geo_rot(exi0,exi0,ang,exi0)

  endif

CASE(KINDPA)
   h=dir*m%l/p%nst
  if(m%pa%hc==0.0_dp) then
   d=(/0.0_dp,0.0_dp,h/)

   call geo_tra(b0,exi0,d,1)
  else
  ang=0.0_dp
  ang(2)=h*m%pa%hc/2
  h=2*sin(ang(2))/m%pa%hc
  d=(/0.0_dp,0.0_dp,h/)
  call geo_rot(exi0,exi0,ang,exi0)
  call geo_tra(b0,exi0,d,1)
  call geo_rot(exi0,exi0,ang,exi0)
  endif

CASE(KINDabell)
   h=dir*m%l/p%nst
  if(m%ab%hc==0.0_dp) then
   d=(/0.0_dp,0.0_dp,h/)

   call geo_tra(b0,exi0,d,1)
  else
  ang=0.0_dp
  ang(2)=h*m%ab%hc/2



  h=2*sin(ang(2))/m%ab%hc
  d=(/0.0_dp,0.0_dp,h/)
  call geo_rot(exi0,exi0,ang,exi0)
  call geo_tra(b0,exi0,d,1)
  call geo_rot(exi0,exi0,ang,exi0)
  endif

CASE(KIND16)
   h=dir*m%l/p%nst
   d=(/0.0_dp,0.0_dp,h/)
   call geo_tra(b0,exi0,d,1)

CASE default
 write(6,*) " not supported in survey_integration_node_case0 "
 stop

end select


! if(f%dir==1) then
t%b=b0
t%exi=exi0
!else
!t%a=b0
!t%ent=exi0
!endif
ent0=exi0

end subroutine survey_integration_node_case0





subroutine survey_integration_special_superdrift(t,ent0,a0)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
type(superdrift),pointer :: el
real(dp) pix1(3) ,pix2(3) ,a0(3),exi0(3,3),ent0(3,3)
logical(lp) :: ENTERING=my_false

 

f=>t%parent_fibre
el=>f%mag%sdr

pix1=0.0_dp; pix2=0;

t%a=a0
t%ent=ent0
exi0=t%ent


if(t%cas==case1.and.f%dir==1) then
!if(entering) then
! f%chart%f%ent=ent0
! f%chart%f%a=a0
!endif
!  f%mag%p%f%ent=ent0
!  f%mag%p%f%a=a0
!  f%magp%p%f%ent=ent0
!  f%magp%p%f%a=a0


if(el%A_X1==-1) pix1(1)=pi
if(el%A_X2==-1) pix2(1)=pi
call GEO_ROT(exi0,pix1,1, ent0) ! new
 pix1=0
 pix1(1)=el%ang(1)
! call GEO_ROT(exi0,pix1,1, ent0)
 call GEO_ROT(exi0,pix1,1, exi0)  !new
 pix1=0
 pix1(2)=el%ang(2)
 call GEO_ROT(exi0,pix1,1, exi0)
 pix1=0
 pix1(3)=el%ang(3)
 call GEO_ROT(exi0,pix1,1, exi0)
call TRANSLATE_point(a0,el%D,1,exi0)  
call GEO_ROT(exi0,pix2,1, exi0)  ! new
!if(.not.entering) then
! f%chart%f%ent=exi0
! f%chart%f%a=a0
!endif

!elseif(t%cas==case1.and.f%dir==-1) then
!  f%chart%f%exi=ent0
!  f%chart%f%b=a0
!  f%mag%p%f%exi=ent0
!  f%mag%p%f%b=a0
!  f%magp%p%f%exi=ent0
!  f%magp%p%f%b=a0
endif

if(t%cas==case2.and.f%dir==-1) then
!if(entering) then
! f%chart%f%ent=ent0
! f%chart%f%a=a0
!endif
!  f%mag%p%f%ent=ent0
!  f%mag%p%f%a=a0
!  f%magp%p%f%ent=ent0
!  f%magp%p%f%a=a0

if(el%A_X2==-1) pix2(1)=pi

call GEO_ROT(exi0,pix2,1, ent0)  ! new

 el%D(1)=-el%D(1)
 el%D(2)=-el%D(2)
 !  call TRANSLATE_point(a0,el%D,1,ent0)  
  call TRANSLATE_point(a0,el%D,1,exi0)  
 el%D(1)=-el%D(1)
 el%D(2)=-el%D(2)
 pix1=0
 pix1(3)=-el%ang(3)
 call GEO_ROT(exi0,pix1,1, exi0)

!call GEO_ROT(ent0,pix1,1, exi0)
 pix1=0
 pix1(2)=el%ang(2)
 call GEO_ROT(exi0,pix1,1, exi0)
 pix1=0
 pix1(1)=el%ang(1)
 call GEO_ROT(exi0,pix1,1, exi0)
pix1=0
if(el%A_X1==-1) pix1(1)=pi
 call GEO_ROT(exi0,pix1,1, exi0)

!if(.not.entering) then
!  f%chart%f%ent=exi0
!  f%chart%f%a=a0
!endif

!elseif(t%cas==case2.and.f%dir==1) then
!  f%chart%f%exi=ent0
!  f%chart%f%b=a0
!  f%mag%p%f%exi=ent0
!  f%mag%p%f%b=a0
!  f%magp%p%f%exi=ent0
 ! f%magp%p%f%b=a0
endif

t%b=a0
!t%ent=ent0   ! mistake????
t%exi=exi0
ent0=exi0
!t%next%a=t%b
!t%next%ent=t%exi

end subroutine survey_integration_special_superdrift


 SUBROUTINE ADJUST_cav_frame(EL,exi0,a0,J,dir)
    IMPLICIT NONE
    real(dp), target :: a0(3),exi0(3,3)
    TYPE(CAV4),INTENT(INOUT):: EL
    INTEGER, INTENT(IN) :: J,dir
    real(dp) d(3),ang(3)
    d=0
    ang=0

    IF(J==1) then
     d(3)= dir*el%H1
        call GEO_ROT(exi0,ang,1, exi0)
        call TRANSLATE_point(a0,D,1,exi0)  
    else
    d(3)= dir*el%H2
 
        call TRANSLATE_point(a0,D,1,exi0)  
        call GEO_ROT(exi0,ang,1, exi0)
    endif
 
  END SUBROUTINE ADJUST_cav_frame

 SUBROUTINE ADJUST_PANCAKE_frame(EL,exi0,a0,J)
    IMPLICIT NONE
    real(dp), target :: a0(3),exi0(3,3)
    TYPE(PANCAKE),INTENT(INOUT):: EL
    INTEGER, INTENT(IN) :: J
    real(dp) d(3),ang(3)
    d=0
    ang=0
    if(el%hc==0.0_dp) then  !<------ Rectangular geometry

    IF(J==1) then
    d(1)=el%xc; d(3)=el%dc; d(2)=el%vc; 
        ang(2)=el%angc
        call GEO_ROT(exi0,ang,1, exi0)
        call TRANSLATE_point(a0,D,1,exi0)  
    else
    d(1)=-el%xc ;d(3)=el%dc;d(2)=-el%vc;
        ang(2)=el%angc
        call TRANSLATE_point(a0,D,1,exi0)  
        call GEO_ROT(exi0,ang,1, exi0)
    endif
    else  !<------ Sector geometry
    IF(J==1) then
    d(1)=el%xc; d(3)=el%dc;d(2)=el%vc;
        ang(2)=el%angc
        call TRANSLATE_point(a0,D,1,exi0)  
        call GEO_ROT(exi0,ang,1, exi0)
    else
    d(1)=-el%xc; d(3)=el%dc;d(2)=-el%vc;
         ang(2)=el%angc
        call GEO_ROT(exi0,ang,1, exi0)
        call TRANSLATE_point(a0,D,1,exi0)  
    endif
    endif
  END SUBROUTINE ADJUST_PANCAKE_frame

 SUBROUTINE ADJUST_abell_frame(EL,exi0,a0,J)
    IMPLICIT NONE
    real(dp), target :: a0(3),exi0(3,3)
    TYPE(abell),INTENT(INOUT):: EL
    INTEGER, INTENT(IN) :: J
    real(dp) d(3),ang(3)
    d=0
    ang=0
    if(el%hc==0.0_dp) then  !<------ Rectangular geometry

    IF(J==1) then
    d(1)=el%xc; d(3)=el%dc; d(2)=el%vc; 
        ang(2)=el%angc
        call GEO_ROT(exi0,ang,1, exi0)
        call TRANSLATE_point(a0,D,1,exi0)  
    else
    d(1)=-el%xc ;d(3)=el%dc;d(2)=-el%vc;
        ang(2)=el%angc
        call TRANSLATE_point(a0,D,1,exi0)  
        call GEO_ROT(exi0,ang,1, exi0)
    endif
    else  !<------ Sector geometry
    IF(J==1) then
    d(1)=el%xc; d(3)=el%dc;d(2)=el%vc;
        ang(2)=el%angc
        call TRANSLATE_point(a0,D,1,exi0)  
        call GEO_ROT(exi0,ang,1, exi0)
    else
    d(1)=-el%xc; d(3)=el%dc;d(2)=-el%vc;
         ang(2)=el%angc
        call GEO_ROT(exi0,ang,1, exi0)
        call TRANSLATE_point(a0,D,1,exi0)  
    endif
    endif
  END SUBROUTINE ADJUST_abell_frame



 SUBROUTINE MIS_survey(a0,ent0,C,b0,exi0,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(FIBRE),target,INTENT(INOUT):: C
    real(dp), target :: a0(3),ent0(3,3),b0(3),exi0(3,3)
    real(dp) ang(3),d(3)
    logical(lp),INTENT(IN)::  ENTERING
    exi0=ent0
    b0=a0
    IF(ASSOCIATED(C%CHART)) THEN
       IF(C%DIR==1) THEN   ! FORWARD PROPAGATION
          IF(ENTERING) THEN
             call GEO_ROT(exi0,c%chart%ANG_IN,1, exi0)
  !           CALL ROT_YZ(C%CHART%ANG_IN(1),X,C%MAG%P%BETA0,OU,k%TIME)   ! ROTATIONS
  !           CALL ROT_XZ(C%CHART%ANG_IN(2),X,C%MAG%P%BETA0,OU,k%TIME)
  !           CALL ROT_XY(C%CHART%ANG_IN(3),X)  !,OU)
            call TRANSLATE_point(b0,c%chart%D_IN,1,exi0)  
  !           CALL TRANS(C%CHART%D_IN,X,C%MAG%P%BETA0,OU,k%TIME)         ! TRANSLATION
          ELSE
             call GEO_ROT(exi0,c%chart%ANG_OUT,1, exi0)
  !           CALL ROT_YZ(C%CHART%ANG_OUT(1),X,C%MAG%P%BETA0,OU,k%TIME)  ! ROTATIONS
  !           CALL ROT_XZ(C%CHART%ANG_OUT(2),X,C%MAG%P%BETA0,OU,k%TIME)
  !           CALL ROT_XY(C%CHART%ANG_OUT(3),X)  !,OU)
  !           CALL TRANS(C%CHART%D_OUT,X,C%MAG%P%BETA0,OU,k%TIME)        ! TRANSLATION
            call TRANSLATE_point(b0,c%chart%D_OUT,1,exi0) 
          ENDIF
       ELSE
          IF(ENTERING) THEN  ! BACKWARD PROPAGATION
              d=c%chart%D_OUT
              ang=C%CHART%ANG_OUT
    
            !   d(1:2)=-d(1:2)   
             call TRANSLATE_point(b0,d,-1,exi0)  
               d=ang
               ang=0
               ang(3)=-d(3)  
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(2)=-d(2)
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(1)=-d(1)
             call GEO_ROT(exi0,ang,1, exi0)
          ELSE
              d=C%CHART%D_IN
              ang=C%CHART%ANG_IN
                 
          !     d(1:2)=-d(1:2)   
             call TRANSLATE_point(b0,d,-1,exi0)  
               d=ang
               ang=0
               ang(3)=-d(3)    !!! strange
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(2)=-d(2)      !!! strange
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(1)=-d(1)
             call GEO_ROT(exi0,ang,1, exi0)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_survey

  subroutine set_aperture_all_case0(flag)
    implicit none
    logical flag 
    
    aperture_all_case0 = flag
    
  end subroutine set_aperture_all_case0


!!!!!!!!!!!!!!!!!!   Old Survey    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE  MISALIGN_FIBRE_EQUAL(S2,S1) ! MISALIGNS FULL FIBRE; FILLS IN CHART AND MAGNET_CHART
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    TYPE(FIBRE),target,INTENT(INOUT):: S2

    CALL MISALIGN_FIBRE(S2,S1)

  END SUBROUTINE  MISALIGN_FIBRE_EQUAL



SUBROUTINE MOVE_FRAMES(S2,s1,OMEGA,BASIS)
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    TYPE(FIBRE),target,INTENT(INOUT):: S2
    real(dp) d(3),ANGLE(3),OMEGAT(3),BASIST(3,3),e1(3,3),e2(3,3),T_GLOBAL(3),d_in(3),d_out(3)
    type(integration_node), pointer :: t,t2
    integer i  !,k
    TYPE (MAGNET_FRAME), pointer :: F0,F,fmag,fbasis

       CALL ALLOC(F)
       CALL ALLOC(F0)
       CALL ALLOC(fmag)
       CALL ALLOC(fbasis)
 
       DO I=1,3
          D(I)=S1(I);   !D(I)=S1(I);
          ANGLE(I)=S1(3+I); !R(I)=S1(3+I);
       ENDDO
      
  
       IF(PRESENT(BASIS)) THEN
          BASIST=BASIS
       ELSE
          BASIST=S2%mag%p%f%mid
       ENDIF
       IF(PRESENT(OMEGA)) THEN
          OMEGAT=OMEGA
       ELSE
          OMEGAT=S2%mag%p%f%o
       ENDIF
       f%mid=global_frame
       f%o=0
       t=>s2%t1%next
!k=0
       do while (.not.(associated(t,s2%t2)) )
        F%ent=t%ent
        F%exi=t%exi
        F%a=t%a
        F%b=t%b
        fmag=S2%mag%p%f
!      write(6,*) ANGLE
       CALL ROTATE_FRAME(F,OMEGAT,ANGLE,1,BASIS=BASIST)
       CALL ROTATE_FRAME(Fmag,OMEGAT,ANGLE,1,BASIS=BASIST)

       IF(PRESENT(BASIS)) THEN   ! MUST ROTATE THAT FRAME AS WELL FOR CONSISTENCY IN DEFINITION WHAT A MISALIGNMENT IS IN PTC
              fbasis%mid=basist
              fbasis%o=omegat
              CALL ROTATE_FRAME(fbasis,OMEGAT,ANGLE,1,BASIS=BASIST)
              e1=fbasis%mid
       !    e1=basist
        !  CALL   GEO_ROT(BASIST,e2,ANGLE,e1)   ! 2018 correction: agrees with successive rotations followed by rotations
        !   basist=e2
       ELSE
          e1=Fmag%MID    ! ALREADY ROTATED
       ENDIF

       CALL CHANGE_BASIS(D,e1,T_GLOBAL,GLOBAL_FRAME)


        F%A=F%A+T_GLOBAL
        F%O=F%O+T_GLOBAL
        F%B=F%B+T_GLOBAL
 
  !      k=k+1
 
        t%ent=f%ent
        t%exi=f%exi
        t%a=f%a
        t%b=f%b

        t=>t%next
       enddo 
!!!!  old fibre magnet attached frames  !!!!
 

        Fmag%A=Fmag%A+T_GLOBAL
        Fmag%O=Fmag%O+T_GLOBAL
        Fmag%B=Fmag%B+T_GLOBAL
        f0=S2%magp%p%f
        S2%mag%p%f=Fmag
        f=fmag


       CALL COMPUTE_ENTRANCE_ANGLE(F0%ENT,F%ENT,S2%CHART%ANG_IN)
       CALL COMPUTE_ENTRANCE_ANGLE(F%EXI,F0%EXI,S2%CHART%ANG_OUT)
 

       D_IN=F%A-F0%A
       D_OUT=F0%B-F%B

 

       CALL CHANGE_BASIS(D_IN,GLOBAL_FRAME,S2%CHART%D_IN,F%ENT)
       CALL CHANGE_BASIS(D_OUT,GLOBAL_FRAME,S2%CHART%D_OUT,F0%EXI)
 
       s2%mag%mis=.true.
       s2%magp%mis=.true.

       CALL kill(F)
       CALL kill(F0)
       CALL kill(fmag)
       CALL kill(fbasis)


    END SUBROUTINE MOVE_FRAMES






  recursive SUBROUTINE  MISALIGN_FIBRE(S2,S1,OMEGA,BASIS,ADD,preserve_girder)
    ! MISALIGNS FULL FIBRE; FILLS IN CHART AND MAGNET_CHART
    ! changed  add=true add extra misalignments TO EXISTING ONES
    ! O AND MID BY DEFAUTL, OTHERWISE OMEGA AND BASIS
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD,preserve_girder
    TYPE(FIBRE),target,INTENT(INOUT):: S2
    REAL(DP) ANGLE(3)   !,T_GLOBAL(3)  ,d(3),r(3)
 
    REAL(DP) D_IN(3),D_OUT(3),OMEGAT(3),BASIST(3,3)
    INTEGER I
    LOGICAL(LP) ADDIN,pres
    LOGICAL(LP) FOUNDg
    type(element), pointer :: caf
    real(dp) a1(3),e1(3,3),a2(3),e2(3,3),dg1(3),ag1(3),mis(6)



    ADDIN=.FALSE.
    pres=.FALSE.
    IF(PRESENT(ADD)) ADDIN=ADD
    if(present(preserve_girder)) pres=preserve_girder

    FOUNDg=.false.
    if(pres.and.associated(s2%mag%girderS).and.(.not.addin)) then
       CALL FIND_AFFINE_GIRDER(S2,CAF,FOUNDg)
       if(foundg) then
          a1=caf%girder_frame%a
          e1=caf%girder_frame%ent
          a2=caf%girder_frame%b
          e2=caf%girder_frame%exi
          !         call INVERSE_FIND_PATCH(a1,e1,dg1,ag1,a2,e2)
          call FIND_PATCH(a1,e1,a2,e2,dg1,ag1)
          mis=0.0_dp
          call MISALIGN_fibre(S2,MIS)
          MIS(1:3)=DG1
          MIS(4:6)=AG1

          call MISALIGN_fibre(S2,MIS,OMEGA=a1,BASIS=e1)
          call MISALIGN_fibre(S2,S1,OMEGA,BASIS,ADD=my_true,preserve_girder=my_false)
          return
       endif
    endif



    IF(ASSOCIATED(S2%CHART)) THEN
       !       IF(.NOT.ASSOCIATED(S2%MAG%D) ) ALLOCATE(S2%MAG%D(3))
       !       IF(.NOT.ASSOCIATED(S2%MAG%R) ) ALLOCATE(S2%MAG%R(3))
       !       IF(.NOT.ASSOCIATED(S2%MAGP%D)) ALLOCATE(S2%MAGP%D(3))
       !       IF(.NOT.ASSOCIATED(S2%MAGP%R)) ALLOCATE(S2%MAGP%R(3))
       !       DO I=1,3
       !          S2%MAG%D(I)=S1(I);   S2%MAGP%D(I)=S1(I);
       !          S2%MAG%R(I)=S1(3+I); S2%MAGP%R(I)=S1(3+I);
       !       ENDDO
     !  DO I=1,3
     !     D(I)=S1(I);   !D(I)=S1(I);
     !     R(I)=S1(3+I); !R(I)=S1(3+I);
     !  ENDDO
       S2%CHART%D_IN=0.0_dp;S2%CHART%D_OUT=0.0_dp;
       S2%CHART%ANG_IN=0.0_dp;S2%CHART%ANG_OUT=0.0_dp;
       S2%MAG%MIS=.TRUE.
       S2%MAGP%MIS=.TRUE.

       ! ADD CODE HERE
 
       ! MOVE THE ORIGINAL INTERNAL CHART F
       IF(ADDIN) THEN

                  call MOVE_FRAMES(S2,s1,OMEGA,BASIS)
            !      call survey_integration_fibre(s2)
                  call survey_integration_fibre(s2,s2%t1%ent,s2%t1%a) 

       ELSE
          !call survey_integration_fibre(s2)  !,previous=.false.)  
          call survey_integration_fibre(s2,s2%t1%ent,s2%t1%a) 
 
          call MOVE_FRAMES(S2,s1,OMEGA,BASIS)
     !     call survey_integration_fibre(s2)
          call survey_integration_fibre(s2,s2%t1%ent,s2%t1%a) 

       ENDIF
 


    ELSE
       !w_p=0
       !w_p%NC=1
       !w_p%FC='((1X,A72))'
       WRITE(6,'(1X,A39,1X,A16)') " CANNOT MISALIGN THIS FIBRE: NO CHARTS ", S2%MAG%NAME
       ! call !write_e(100)
    ENDIF


  END SUBROUTINE MISALIGN_FIBRE


  subroutine print_magnet_framet(m,mf)
    implicit none
    type(magnet_frame), target :: m
    integer mf,i
 
       write(mf,'(a72)') "MAGNET FRAME MAGNET FRAME MAGNET FRAME MAGNET FRAME MAGNET FRAME MAGNET FRAME "
       WRITE(MF,*) m%a
       WRITE(MF,*) m%o
       WRITE(MF,*) m%b
       write(mf,'(a)') " ENTRANCE "

       do i=1,3
          WRITE(MF,*) m%ent(i,1:3)
       enddo
       write(mf,'(a)') " MIDDLE "

       do i=1,3
          WRITE(MF,*) m%mid(i,1:3)
       enddo
       write(mf,'(a)') " EXIT "

       do i=1,3
          WRITE(MF,*) m%exi(i,1:3)
       enddo
       write(mf,'(a68)') "END MAGNET FRAME END MAGNET FRAME END MAGNET FRAME END MAGNET FRAME "
 
  end subroutine print_magnet_framet

  SUBROUTINE  MAD_MISALIGN_FIBRE(S2,S1) ! MISALIGNS FULL FIBRE; FILLS IN CHART AND MAGNET_CHART
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    TYPE(FIBRE),target,INTENT(INOUT):: S2
    REAL(DP) ENT(3,3),ENT1(3,3),ENT2(3,3),T(3),MAD_ANGLE(3),T_GLOBAL(3),ANGLE(3),MIS(6)
    real(dp) r1(3,3),r2(3,3),r3(3,3)
    ent=S2%CHART%F%ent
    T(1)=S1(1);T(2)=S1(2);T(3)=S1(3);
    MAD_ANGLE(1)=-S1(4)
    MAD_ANGLE(2)=-S1(5)
    MAD_ANGLE(3)=S1(6)



    ANGLE=0.0_dp; ANGLE(1)=MAD_ANGLE(1)
    ent1=ent
    ent2=ent
    CALL GEO_ROT(ENT1,ENT,ANGLE,ENT2)
 
    ANGLE=0.0_dp; ANGLE(2)=MAD_ANGLE(2)
    ent1=ent
    ent2=ent
    CALL GEO_ROT(ENT1,ENT,ANGLE,ENT2)

    ANGLE=0.0_dp; ANGLE(3)=MAD_ANGLE(3)
    ent1=ent
    ent2=ent
    CALL GEO_ROT(ENT1,ENT,ANGLE,ENT2)
    




    CALL COMPUTE_ENTRANCE_ANGLE(S2%CHART%F%ent,ENT,ANGLE)

 
     ent1=S2%CHART%F%ent
     ent2=ent
     CALL CHANGE_BASIS(T,ENT1,T_GLOBAL,ent2)

 

    MIS(1:3)=T_GLOBAL
    MIS(4:6)=ANGLE

    ENT=S2%CHART%F%ent
    T=S2%CHART%F%A

        call MISALIGN_FIBRE(S2,MIS,T,ENT)

  END SUBROUTINE MAD_MISALIGN_FIBRE

  SUBROUTINE  MISALIGN_GIRDER(S2,S1,OMEGA,BASIS,ADD) !
    ! SIMILAR TO MISALIGN_SIAMESE
    ! COMMENT DIFFERENCES ONLY
    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    TYPE(FIBRE),TARGET,INTENT(INOUT):: S2
    TYPE(ELEMENT), POINTER :: C,CN,CAF
    TYPE(fibre), POINTER :: P
    integer k
    REAL(DP) OMEGAT(3),BASIST(3,3),B(3),EXI(3,3),T_GLOBAL(3)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD
    LOGICAL(LP) ADDIN
    LOGICAL(LP) FOUND
    !    REAL(DP) D(3),ANG(3)
    TYPE(MAGNET_FRAME), POINTER :: F

    FOUND=.FALSE.
    ADDIN=.FALSE.
    CALL FIND_AFFINE_GIRDER(S2,CAF,FOUND)
    IF(FOUND) CALL FIND_FRAME_GIRDER(CAF,B,EXI,ADD)

    IF(PRESENT(ADD)) ADDIN=ADD

    IF(PRESENT(OMEGA)) THEN
       OMEGAT=OMEGA
    ELSE
       OMEGAT=S2%CHART%F%O
    ENDIF

    IF(PRESENT(BASIS)) THEN
       BASIST=BASIS
    ELSE
       BASIST=S2%CHART%F%MID
    ENDIF

    IF((.NOT.PRESENT(OMEGA)).AND.(.NOT.PRESENT(BASIS))) THEN
       IF(FOUND) THEN
          OMEGAT=B
          BASIST=EXI
       ENDIF
    ENDIF


    CALL MISALIGN_FIBRE(S2,S1,OMEGAT,BASIST,ADD=ADDIN)
    k=1

    IF(ASSOCIATED(S2%MAG%GIRDERS)) THEN
       C=>S2%MAG
       CN=>S2%MAG%GIRDERS
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          P=>CN%PARENT_FIBRE
          CALL MISALIGN_FIBRE(P,S1,OMEGAT,BASIST,ADD=ADDIN)
          CN=>CN%GIRDERS
          k=k+1
       ENDDO
    ENDIF



    IF(FOUND) THEN   !!! THE ORIGINAL GIRDER FRAME IS STILL GIRDER_FRAME%ENT AND GIRDER_FRAME%A
       !                    FINAL FRAME AFTER MISALIGNMENTS MUST BE COMPUTED
       call alloc(f)
       f%a=b
       f%ent=exi
       CALL ROTATE_FRAME(F,OMEGAT,S1(4:6),1,BASIS=BASIST)
       CALL   GEO_ROT(BASIST,S1(4:6),1)
       CALL CHANGE_BASIS(S1(1:3),BASIST,T_GLOBAL,GLOBAL_FRAME)
       F%A=F%A+T_GLOBAL
       CAF%GIRDER_FRAME%EXI=F%ent
       CAF%GIRDER_FRAME%B=F%A
       call kill(f)
    ENDIF

    if(global_verbose)     write(6,*) k, " magnet misaligned "
  END SUBROUTINE  MISALIGN_GIRDER


  RECURSIVE SUBROUTINE  MISALIGN_SIAMESE(S2,S1,OMEGA,BASIS,ADD,preserve_girder)
    ! SAME AS MISALIGN_FIBRE: DEFAULT IS THE O,MID OF S2
    !   UNLESS IT FINDS  TYPE(AFFINE_FRAME), POINTER :: SIAMESE_FRAME ON S2%SIAMESE CHAIN
    ! ON ONE SIAMESE IN THE CHAIN

    ! THIS IS OVERWRITEN IF OMEGA AND BASIS ARE PRESENT


    IMPLICIT NONE
    REAL(DP),INTENT(IN):: S1(6)
    REAL(DP), OPTIONAL, INTENT(IN) :: OMEGA(3),BASIS(3,3)
    TYPE(FIBRE),TARGET,INTENT(INOUT):: S2
    TYPE(ELEMENT), POINTER :: C,CN
    TYPE(fibre), POINTER :: P
    integer k
    REAL(DP) OMEGAT(3),BASIST(3,3),B(3),EXI(3,3)
    real(dp) a1(3),e1(3,3),a2(3),e2(3,3),dg1(3),ag1(3),mis(6)
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: ADD,preserve_girder
    LOGICAL(LP) ADDIN,pres
    !    TYPE(AFFINE_FRAME),POINTER :: AF
    LOGICAL(LP) FOUND
    LOGICAL(LP) FOUNDg
    type(element), pointer :: caf
    !    REAL(DP) D(3),ANG(3)

    FOUND=.FALSE.
    ADDIN=.FALSE.
    pres=.FALSE.
    if(present(preserve_girder)) pres=preserve_girder
    CALL FIND_AFFINE_SIAMESE(S2,CN,FOUND)  ! Looking for siamese WITH FRAME
    IF(FOUND) CALL FIND_FRAME_SIAMESE(CN,B,EXI,ADD) ! FIND ACTUAL FRAME

    FOUNDg=.false.
    if(pres.and.associated(s2%mag%girders).and.(.not.addin)) then
       CALL FIND_AFFINE_GIRDER(S2,CAF,FOUNDg)
       if(foundg) then
          a1=caf%girder_frame%a
          e1=caf%girder_frame%ent
          a2=caf%girder_frame%b
          e2=caf%girder_frame%exi
          !         call INVERSE_FIND_PATCH(a1,e1,dg1,ag1,a2,e2)
          call FIND_PATCH(a1,e1,a2,e2,dg1,ag1)
          mis=0.0_dp
          call MISALIGN_siamese(S2,MIS)
          MIS(1:3)=DG1
          MIS(4:6)=AG1

          call MISALIGN_SIAMESE(S2,MIS,OMEGA=a1,BASIS=e1)
          call MISALIGN_SIAMESE(S2,S1,OMEGA,BASIS,ADD=my_true,preserve_girder=my_false)
          return
       endif
    endif

    IF(PRESENT(ADD)) ADDIN=ADD

    IF(PRESENT(OMEGA)) THEN    ! Arbitrary Origin
       OMEGAT=OMEGA
    ELSE
       OMEGAT=S2%CHART%F%O   ! Centre of magnet otherwise
    ENDIF

    IF(PRESENT(BASIS)) THEN      ! Arbitrary Basis
       BASIST=BASIS
    ELSE
       BASIST=S2%CHART%F%MID  ! Centre of Magnet Otherwise
    ENDIF

    IF((.NOT.PRESENT(OMEGA)).AND.(.NOT.PRESENT(BASIS))) THEN
       IF(FOUND) THEN   ! If no special basis and no special origin
          OMEGAT=B         ! and siamese is found, then it uses the siamese basis
          BASIST=EXI        ! Notice that if ADD=true, the siamese frames move with the magnets
       ENDIF
    ENDIF

    CALL MISALIGN_FIBRE(S2,S1,OMEGAT,BASIST,ADD=ADDIN)
    k=1

    IF(ASSOCIATED(S2%MAG%SIAMESE)) THEN
       C=>S2%MAG
       CN=>S2%MAG%SIAMESE
       DO WHILE(.NOT.ASSOCIATED(C,CN))
          P=>CN%PARENT_FIBRE
          CALL MISALIGN_FIBRE(P,S1,OMEGAT,BASIST,ADD=ADDIN)
          CN=>CN%SIAMESE
          k=k+1
       ENDDO
    ENDIF
    !    CALL MOVE_SIAMESE_FRAME(S2%MAG)
    if(global_verbose) write(6,*) k, " magnet misaligned "
  END SUBROUTINE  MISALIGN_SIAMESE

  SUBROUTINE FIND_PATCH_0_survey(EL1,EL2_NEXT,NEXT,ENERGY_PATCH,PREC,patching) ! COMPUTES PATCHES
    IMPLICIT NONE
    TYPE (FIBRE),target :: EL1
    TYPE (FIBRE),TARGET,OPTIONAL, INTENT(INOUT) :: EL2_NEXT
    TYPE (FIBRE),POINTER :: EL2
    REAL(DP), OPTIONAL :: PREC
    LOGICAL(LP), OPTIONAL, INTENT(IN) ::  NEXT,ENERGY_PATCH
    LOGICAL(LP), OPTIONAL, INTENT(out) ::  patching
    LOGICAL(LP) patch
     call find_patch(EL1,EL2_NEXT,NEXT,ENERGY_PATCH,PREC,patching=patch)

     if(patch) then
       call survey_integration_fibre(el1)
       if (.not.present(el2_next)) call survey_integration_fibre(el1%next)
     endif
     if(present(patching) ) patching=patch

    end SUBROUTINE FIND_PATCH_0_survey

subroutine convert_mis_to_patch(ptc_fibre,use_ptc)
implicit none
type(fibre), pointer :: ptc_fibre
real(dp) ent(3,3),aent(3),d(3),ang(3)
real(dp) ent0(3,3),aent0(3),exi0f(3,3),bexi0f(3)
real(dp) exi0(3,3),pix(3),EXI(3,3),PREC
logical use_ptc
PREC=1.d-38

    ptc_fibre%mag%p%tiltd=0
    ptc_fibre%magp%p%tiltd=0
    ptc_fibre%chart%d_in=0
    ptc_fibre%chart%d_out=0
    ptc_fibre%chart%ang_in=0
    ptc_fibre%chart%ang_out=0
    ptc_fibre%mag%mis=.false.
    ptc_fibre%magp%mis=.false.
    ptc_fibre%PATCH%A_X1=1
    ptc_fibre%PATCH%B_X1=1
    ptc_fibre%PATCH%A_X2=1
    ptc_fibre%PATCH%B_X2=1
    ptc_fibre%PATCH%PATCH=3



ent=ptc_fibre%t1%next%ent
aent=ptc_fibre%t1%next%a


ent0=ptc_fibre%t1%ent
aent0=ptc_fibre%t1%a
exi0f=ptc_fibre%t2%exi
bexi0f =ptc_fibre%t2%b
 
call survey_integration_fibre(ptc_fibre,ptc_fibre%t1%ent,ptc_fibre%t1%a) 
 
D=aent-ptc_fibre%t1%next%a
 
CALL TRANSLATE_Fibre(ptc_fibre,D)
CALL COMPUTE_ENTRANCE_ANGLE(ptc_fibre%t1%next%ent,ent,ANG)
 

CALL ROTATE_FIBRE(ptc_fibre,aent,Ang,BASIS=ptc_fibre%t1%next%ent)  


if(associated(ptc_fibre%previous).and.use_ptc) then
   call find_patch(ptc_fibre%previous,ptc_fibre,NEXT=.true.,ENERGY_PATCH=.true.,PREC=prec)
else
    pix=0.0_dp
    if(ptc_fibre%dir==-1) then
     ptc_fibre%PATCH%A_X1=-1
     ptc_fibre%PATCH%B_X1=-1
     pix(1)=pi
   endif 
 
             EXI=ent0
             exi0=exi
             call geo_rot(exi,pix,1,basis=exi0)
             ENT=ptc_fibre%t1%next%ent
             exi0=ENT
             call geo_rot(ent,pix,1,basis=exi0)

 call find_patch(aent0,exi,ptc_fibre%t1%next%a,ent,ptc_fibre%patch%a_d,ptc_fibre%patch%a_ang)
 
endif
 

if(associated(ptc_fibre%next).and.use_ptc) then
   call find_patch(ptc_fibre,ptc_fibre%next,NEXT=.false.,ENERGY_PATCH=.true.,PREC=prec)
else

    pix=0.0_dp
    if(ptc_fibre%dir==-1) then
    ptc_fibre%PATCH%A_X2=-1
    ptc_fibre%PATCH%B_X2=-1
     pix(1)=pi
   endif 
             EXI=exi0f
             exi0=exi
             call geo_rot(exi,pix,1,basis=exi0)
             ENT=ptc_fibre%t2%previous%exi
             exi0=ENT
             call geo_rot(ent,pix,1,basis=exi0)

 
 call find_patch(ptc_fibre%t2%previous%b,ent,bexi0f,exi,ptc_fibre%patch%b_d,ptc_fibre%patch%b_ang)
endif
 
 call survey_integration_fibre(ptc_fibre,ent0,aent0) 


end subroutine convert_mis_to_patch
end module ptc_multiparticle
