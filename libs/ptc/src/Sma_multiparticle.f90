
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
  PRIVATE TRACK_NODE_SINGLER,TRACK_NODE_SINGLEP,TRACK_NODE_SINGLEV
  private DRIFTr_BACK_TO_POSITION,DRIFTp_BACK_TO_POSITION  !,DRIFT_BACK_TO_POSITION
  private MAKE_NODE_LAYOUT_2 !,DRIFT_TO_TIME
  PRIVATE MODULATE_R,MODULATE_P
  PRIVATE TRACK_MODULATION_R,TRACK_MODULATION_P
 LOGICAL :: no_mis=.TRUE. 
  !  LOGICAL :: OLD_MOD=.TRUE.

  logical(lp),private, parameter :: dobb=.true.
  logical(lp),private            :: aperture_all_case0=.false.
 ! type(probe) :: xsm,xsm0
  real(dp) :: xsm0t=0.0_dp,xsmt=0.0_dp
  !real(dp) :: unit_time =1.0e-3_dp
  REAL(dp) :: x_orbit_sync(6)= 0.0_dp,dt_orbit_sync=0.0_dp
    logical(lp) :: use_bmad_units=.false.,inside_bmad=.false.

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
          call BBKICK(t%bb,X)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
 
       SELECT CASE(EL%KIND)
       CASE(KIND0)
         global_e= x(5)*el%p%p0c
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X,K,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND3)
          CALL TRACK(EL%K3,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X,K,t%POS_IN_FIBRE-2)
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
          CALL TRACK_SLICE(EL%TP10,X,K,t%POS_IN_FIBRE-2)
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
          CALL TRACK_SLICE(EL%K16,X,K,t%POS_IN_FIBRE-2)
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
          call BBKICK(t%bb,X)
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
          call BBKICK(t%bb,X)
          if(t%bb%patch)call PATCH_BB(t%bb,X,k,EL%p%BETA0,ALWAYS_EXACT_PATCHING.or.EL%P%EXACT,my_false)

       endif
       SELECT CASE(EL%KIND)
       CASE(KIND0)
         global_e= x(5)*el%p%p0c
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X,K,t%POS_IN_FIBRE-2)
         global_e= x(5)*el%p%p0c
       case(KIND3)
          CALL TRACK(EL%K3,X,K)
         global_e= x(5)*el%p%p0c
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X,K,t%POS_IN_FIBRE-2)
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
          CALL TRACK_SLICE(EL%TP10,X,K,t%POS_IN_FIBRE-2)
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
          CALL TRACK_SLICE(EL%K16,X,K,t%POS_IN_FIBRE-2)
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
          call BBKICK(t%bb,X)
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
      z%x(6)=sqrt(1.0_dp +2*t/b0+t**2)-1.d0 
     endif
     call kill(t)


     end subroutine convert_ptc_to_bmadp 

     subroutine in_bmad_units
     implicit none  
      use_bmad_units=.true.
      ndpt_bmad=1
     end subroutine in_bmad_units

     subroutine in_ptc_units
     implicit none  
      use_bmad_units=.false.
      ndpt_bmad=0
     end subroutine in_ptc_units

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


subroutine survey_integration_layout(p,f,a,ent)
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
 a0=p%chart%f%a          !global_origin
endif
if(present(ent)) then
 ent0=ent
else
 !ent0=global_FRAME
 ent0=p%chart%f%ent !global_FRAME
endif

 
call survey_integration_fibre(p,a0,ent0)
 
 
p1=>p%next
if(present(f)) then
 p2=>f
else
 p2=>p
endif
 
do while(.not.associated(p2,p1))
 
 
call survey_integration_fibre(p1,p1%previous%t2%b,p1%previous%t2%exi)
 
!call survey_integration_fibre(p1,p1%previous%chart%f%b,p1%previous%chart%f%exi)
 
p1=>p1%next
enddo
 
end subroutine survey_integration_layout


subroutine survey_integration_fibre(p,b0,exi0)
implicit none
type(fibre), target :: p
type(integration_node), pointer :: t
integer i
type(layout), pointer  :: r
real(dp),intent(in):: b0(3),exi0(3,3)
real(dp) a0(3),ent0(3,3),ang(3) 
a0=b0
ent0=exi0
r=>p%parent_layout
if(.not.associated(r%t)) then
 call make_node_layout(r)
 call survey(r)
    CALL  allocate_node_frame( R)
! call FILL_SURVEY_DATA_IN_NODE_LAYOUT(r)
endif
if(.not.associated(p%t1%a))     CALL  allocate_node_frame( R)   !call FILL_SURVEY_DATA_IN_NODE_LAYOUT(r)
 
 
t=>p%t1
 
 
call survey_integration_node_p1(t,a0,ent0)


t=>t%next
call survey_integration_fringe(t,a0,ent0)

do i=1,p%mag%p%nst
 t=>t%next
 call survey_integration_node_case0(t,a0,ent0)
enddo
t=>t%next
call survey_integration_fringe(t,a0,ent0)
t=>t%next
call survey_integration_node_p2(t,a0,ent0)

 

 
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
p%magp%p%f%mid=p%mag%p%f%mid
p%magp%p%f%o=p%mag%p%f%o

 


end subroutine survey_integration_fibre



subroutine survey_integration_fringe(t,a0,ent0)
implicit none
type(integration_node), target :: t
real(dp) a0(3),ent0(3,3)

if(associated(t%parent_fibre%mag%sdr)) then
 call survey_integration_special_superdrift(t,a0,ent0)
else
 t%a=a0
 t%ent=ent0
 !t%parent_fibre%chart%f%a=a0
 !t%parent_fibre%chart%f%ent=ent0
 !t%parent_fibre%chart%f%b=a0
 !t%parent_fibre%chart%f%exi=ent0
 t%b=t%a
 t%exi=t%ent
 a0=t%b
 ent0=t%exi
! t%next%b=t%a
! t%next%ent=t%exi
endif
end subroutine survey_integration_fringe




subroutine survey_integration_node_case0(t,b0,ent0)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
type(element), pointer :: m
type(magnet_chart), pointer :: p
real(dp) h,d(3),ang(3),b0(3),exi0(3,3),ent0(3,3)

f=>t%parent_fibre
m=>f%mag
p=>m%p

t%a=b0
t%ent=ent0
exi0=t%ent
 
select case(m%kind) 

CASE(KIND0,KIND1,KIND3:KIND5,KIND8:KIND9,KIND11:KIND15,KIND17:KIND22,kindwiggler,kindsuperdrift)
   h=p%lc/p%nst
   d=(/0.0_dp,0.0_dp,h/)

   call geo_tra(b0,exi0,d,1)
CASE(KIND2,KIND6:KIND7,KIND10)
   h=p%ld/p%nst
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
   h=m%l/p%nst
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
   h=m%l/p%nst
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
   h=m%l/p%nst
   d=(/0.0_dp,0.0_dp,h/)
   call geo_tra(b0,exi0,d,1)

CASE default
 write(6,*) " not supported in survey_integration_node_case0 "
 stop

end select
t%b=b0
t%exi=exi0

ent0=exi0
!t%next%a=t%b
!t%next%ent=t%exi

end subroutine survey_integration_node_case0

subroutine survey_integration_node_p1(t,a0,ent0)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
real(dp) pix1(3),pix2(3) ,a0(3),exi0(3,3),ent0(3,3)
logical(lp) :: ENTERING=my_true

f=>t%parent_fibre
pix1=0.0_dp;pix2=0.0_dp;

t%a=a0
t%ent=ent0
exi0=t%ent

if(f%patch%A_X1==-1) pix1(1)=pi
if(f%patch%A_X2==-1) pix2(1)=pi

!
call GEO_ROT(exi0,pix1,1, ent0)
call GEO_ROT(exi0,f%patch%a_ang,1, exi0)
call TRANSLATE_point(a0,f%patch%A_D,1,exi0)  
call GEO_ROT(exi0,pix2,1, exi0)

!!! entrance chart
f%chart%f%ent=exi0
f%chart%f%a=a0


pix1=0.0_dp
pix1(3)=f%MAG%P%TILTD
 call GEO_ROT(exi0,pix1,1, exi0)




pix1=0.0_dp
if(f%mag%p%exact.and.(f%mag%kind/=kind10)) then
if(f%dir==1) then
 pix1(2)=f%MAG%P%edge(1)
else
 pix1(2)=f%MAG%P%edge(2)
endif

 call GEO_ROT(exi0,pix1,1, exi0)
endif


    IF(f%MAG%MIS) THEN
      call MIS_survey(a0,exi0,f,a0,exi0,ENTERING)
    ENDIF
 
if(f%mag%kind==kindpa) then

call ADJUST_PANCAKE_frame(f%mag%pa,a0,exi0,1)
!
write(6,*) " I am here in survey_integration_node_p1 "
endif
if(f%mag%kind==kindabell) then

call ADJUST_abell_frame(f%mag%ab,a0,exi0,1)
!
write(6,*) " I am here in survey_integration_node_p1 "
endif  
t%b=a0
!t%ent=ent0   ! mistake????
t%exi=exi0


ent0=exi0

 if(f%dir==1) then
  f%mag%p%f%ent=ent0
  f%mag%p%f%a=a0
  f%magp%p%f%ent=ent0
  f%magp%p%f%a=a0
else
  f%mag%p%f%exi=ent0
  f%mag%p%f%b=a0
  f%magp%p%f%exi=ent0
  f%magp%p%f%b=a0
endif


!t%next%a=t%b
!t%next%ent=t%exi
 

end subroutine survey_integration_node_p1


subroutine survey_integration_special_superdrift(t,a0,ent0)
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



 SUBROUTINE ADJUST_PANCAKE_frame(EL,a0,exi0,J)
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

 SUBROUTINE ADJUST_abell_frame(EL,a0,exi0,J)
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

subroutine survey_integration_node_p2(t,a0,ent0)
implicit none
type(integration_node), target :: t
type(fibre), pointer :: f
real(dp) pix1(3),pix2(3),ent0(3,3),a0(3),exi0(3,3)
logical(lp) :: ENTERING=my_FALSE

f=>t%parent_fibre

 
t%a=a0
t%ent=ent0
exi0=t%ent

 if(f%dir==-1) then
  f%mag%p%f%ent=ent0
  f%mag%p%f%a=a0
  f%magp%p%f%ent=ent0
  f%magp%p%f%a=a0
else
  f%mag%p%f%exi=ent0
  f%mag%p%f%b=a0
  f%magp%p%f%exi=ent0
  f%magp%p%f%b=a0
endif


pix1=0.0_dp
if(f%mag%p%exact.and.(f%mag%kind/=kind10)) then
if(f%dir==1) then
 pix1(2)=f%MAG%P%edge(2)
else
 pix1(2)=f%MAG%P%edge(1)
endif
 call GEO_ROT(exi0,pix1,1, exi0)
endif

if(f%mag%kind==kindpa) then

call ADJUST_PANCAKE_frame(f%mag%pa,a0,exi0,2)
!
write(6,*) " I am here in survey_integration_node_p1 "
endif 

if(f%mag%kind==kindabell) then

call ADJUST_abell_frame(f%mag%ab,a0,exi0,2)
!
write(6,*) " I am here in survey_integration_node_p1 "
endif  

    IF(f%MAG%MIS) THEN
      call MIS_survey(a0,exi0,f,a0,exi0,ENTERING)
    ENDIF


pix1=0.0_dp
pix1(3)=-f%MAG%P%TILTD
 call GEO_ROT(exi0,pix1,1, exi0)

f%chart%f%exi=exi0
f%chart%f%b=a0



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
!t%ent=ent0   ! mistake????
t%exi=exi0

ent0=exi0
!t%next%a=t%b
!t%next%ent=t%exi



!       X(3)=C%PATCH%B_X1*X(3);X(4)=C%PATCH%B_X1*X(4);
!       CALL ROT_YZ(C%PATCH%B_ANG(1),X,C%MAG%P%BETA0,PATCH,k%TIME)
!       CALL ROT_XZ(C%PATCH%B_ANG(2),X,C%MAG%P%BETA0,PATCH,k%TIME)
!       CALL ROT_XY(C%PATCH%B_ANG(3),X)  !,PATCH)
!       CALL TRANS(C%PATCH%B_D,X,C%MAG%P%BETA0,PATCH,k%TIME)
!       X(3)=C%PATCH%B_X2*X(3);X(4)=C%PATCH%B_X2*X(4);

end subroutine survey_integration_node_p2

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
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)

             call TRANSLATE_point(b0,d,-1,exi0)  
!             CALL TRANS(C%CHART%D_OUT,X,C%MAG%P%BETA0,OU,k%TIME)        ! TRANSLATION
!             CALL ROT_XY(C%CHART%ANG_OUT(3),X)  !,OU)
!             CALL ROT_XZ(C%CHART%ANG_OUT(2),X,C%MAG%P%BETA0,OU,k%TIME)
!             CALL ROT_YZ(C%CHART%ANG_OUT(1),X,C%MAG%P%BETA0,OU,k%TIME)  ! ROTATIONS
             d=ang
             ang=0.d0
             ang(3)=-d(3)
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(2)=-d(2)
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(1)=-d(1)
             call GEO_ROT(exi0,ang,1, exi0)
             C%CHART%D_OUT(1)=-C%CHART%D_OUT(1)
             C%CHART%D_OUT(2)=-C%CHART%D_OUT(2)
             C%CHART%ANG_OUT(3)=-C%CHART%ANG_OUT(3)
          ELSE
              d=C%CHART%D_IN
              ang=C%CHART%ANG_IN
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
             call TRANSLATE_point(b0,d,-1,exi0)
!             CALL TRANS(C%CHART%D_IN,X,C%MAG%P%BETA0,OU,k%TIME)         ! TRANSLATION
!             CALL ROT_XY(C%CHART%ANG_IN(3),X)  !,OU)
!             CALL ROT_XZ(C%CHART%ANG_IN(2),X,C%MAG%P%BETA0,OU,k%TIME)
!             CALL ROT_YZ(C%CHART%ANG_IN(1),X,C%MAG%P%BETA0,OU,k%TIME)   ! ROTATIONS
             d=ang
             ang=0.d0
             ang(3)=-d(3)
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(2)=-d(2)
             call GEO_ROT(exi0,ang,1, exi0)  
             ang=0.d0
             ang(1)=-d(1)
             call GEO_ROT(exi0,ang,1, exi0)
             C%CHART%D_IN(1)=-C%CHART%D_IN(1)
             C%CHART%D_IN(2)=-C%CHART%D_IN(2)
             C%CHART%ANG_IN(3)=-C%CHART%ANG_IN(3)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE MIS_survey

  subroutine set_aperture_all_case0(flag)
    implicit none
    logical flag 
    
    aperture_all_case0 = flag
    
  end subroutine set_aperture_all_case0

!!!!!!!!!!!!!!!!!!!


subroutine survey_integration_fibre_no_patch(p,b0,exi0)
implicit none
type(fibre), target :: p
type(integration_node), pointer :: t
integer i
type(layout), pointer  :: r
real(dp),intent(in):: b0(3),exi0(3,3) 
real(dp) a0(3),ent0(3,3),ang(3)
a0=b0
ent0=exi0
r=>p%parent_layout
if(.not.associated(r%t)) then
 call make_node_layout(r)
 call survey(r)
    CALL  allocate_node_frame( R)
! call FILL_SURVEY_DATA_IN_NODE_LAYOUT(r)
endif
if(.not.associated(p%t1%a))     CALL  allocate_node_frame( R)   !call FILL_SURVEY_DATA_IN_NODE_LAYOUT(r)
 


t=>p%t1
!call survey_integration_node_p1(t,a0,ent0)
t=>t%next
call survey_integration_fringe(t,a0,ent0)

do i=1,p%mag%p%nst
 t=>t%next
 call survey_integration_node_case0(t,a0,ent0)
enddo
t=>t%next
call survey_integration_fringe(t,a0,ent0)
!t=>t%next
!call survey_integration_node_p2(t,a0,ent0)

 
!!! entrance chart
 CALL COMPUTE_ENTRANCE_ANGLE(p%chart%f%ent,p%chart%f%exi,ANG)
p%chart%f%mid=p%chart%f%ent
!write(6,*) p%mag%name
!write(6,*) ang
ang=ang/2
p%chart%f%o=0.5_dp*(p%chart%f%a+p%chart%f%b)
CALL GEO_ROT(p%chart%f%mid,ANG,1,basis=p%chart%f%ent)
 
CALL COMPUTE_ENTRANCE_ANGLE(p%mag%p%f%ent,p%mag%p%f%exi,ANG)
p%mag%p%f%mid=p%mag%p%f%ent
ang=ang/2
p%mag%p%f%o=0.5_dp*(p%mag%p%f%a+p%mag%p%f%b)
CALL GEO_ROT(p%mag%p%f%mid,ANG,1,basis=p%mag%p%f%ent)
p%magp%p%f%mid=p%mag%p%f%mid
p%magp%p%f%o=p%mag%p%f%o

 
! I am here 
!etienne

end subroutine survey_integration_fibre_no_patch


end module ptc_multiparticle
