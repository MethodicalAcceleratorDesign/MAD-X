!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_def_all_kinds
  use S_status
  implicit none
  public
  private XMIDR,GMIDR,ALLOC_FIBRE
  include "a_def_worm.inc"
  !  include "a_def_all_kind.inc"
  !  include "a_def_sagan.inc"
  !  include "a_def_user1.inc"
  !!  include "a_def_arbitrary.inc"
  !  include "a_def_user2.inc"
  !  include "a_def_element_fibre_layout.inc"
  private ALLOC_midr,KILL_midr


  INTERFACE XMID
     MODULE PROCEDURE XMIDR
  END  INTERFACE

  INTERFACE GMID
     MODULE PROCEDURE GMIDR
  END  INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_midr
     MODULE PROCEDURE ALLOC_FIBRE
  END  INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILL_midr
  END  INTERFACE


contains

  !  RECURSIVE
  SUBROUTINE GET_LENGTH(R,L)
    IMPLICIT NONE
    TYPE(LAYOUT),target, INTENT(IN) :: R
    REAL(DP), INTENT(OUT) :: L
    TYPE(FIBRE), POINTER:: P

    INTEGER I
    P=>R%START
    L=0.0_dp
    DO I=1,R%N
       IF(P%MAG%KIND/=KIND23) THEN
          L=L+P%MAG%P%LD
          !       ELSE
          !          CALL GET_LENGTH(P%MAG%G23,LG)
          !          L=L+LG
       ENDIF
       P=>P%NEXT
    ENDDO
  END SUBROUTINE GET_LENGTH

  SUBROUTINE XFRAME(E_IN,ENT,A,I)
    IMPLICIT NONE
    TYPE(INNER_FRAME), INTENT(INOUT):: E_IN
    REAL(DP),optional, INTENT(IN):: ENT(3,3),A(3)
    INTEGER, INTENT(IN):: i
    INTEGER J,K

    IF(I<=SIZE(E_IN%ORIGIN,2)) THEN
       DO J=1,3
          ! write(6,*) i,SIZE(E_IN%ORIGIN,1),SIZE(E_IN%ORIGIN,2)
          ! write(6,*) LBOUND(E_IN%ORIGIN,DIM=1),uBOUND(E_IN%ORIGIN,DIM=1),LBOUND(E_IN%ORIGIN,DIM=2),uBOUND(E_IN%ORIGIN,DIM=2)
          if(PRESENT(A)) E_IN%ORIGIN(J,I)=A(J)
          DO K=1,3
             if(PRESENT(ENT)) E_IN%FRAME(J,K,I)=ENT(J,K)
          ENDDO
       ENDDO
    ELSE
       WRITE(6,*) I
    ENDIF
  END SUBROUTINE XFRAME


  SUBROUTINE G_FRAME(E_IN,ENT,A,I)
    IMPLICIT NONE
    TYPE(INNER_FRAME), INTENT(IN):: E_IN
    REAL(DP), INTENT(INOUT):: ENT(3,3),A(3)
    INTEGER, INTENT(IN):: i
    INTEGER J,K

    IF(I<=SIZE(E_IN%ORIGIN,2)) THEN
       DO J=1,3
          A(J)=E_IN%ORIGIN(J,I)
          DO K=1,3
             ENT(J,K)=E_IN%FRAME(J,K,I)
          ENDDO
       ENDDO
    ELSE
       WRITE(6,*) "ERROR IN GFRAME "
       WRITE(6,*) I,SIZE(E_IN%ORIGIN,2)
       STOP 345
    ENDIF
  END SUBROUTINE G_FRAME

  SUBROUTINE XMIDR(X_IN,X,I)
    IMPLICIT NONE
    TYPE(worm), INTENT(INOUT):: X_IN
    REAL(DP), INTENT(IN) :: X(6)
    INTEGER, INTENT(IN):: i
    INTEGER J

    X_IN%nst=i
    IF(I<=SIZE(X_IN%RAY,2)) THEN
       DO J=1,6
          X_IN%RAY(J,I)=X(J)
       ENDDO
    ELSE
       WRITE(6,*) I
       STOP 8
    ENDIF
  END SUBROUTINE XMIDR



  SUBROUTINE gMIDR(X_IN,X,I)
    IMPLICIT NONE
    TYPE(worm), INTENT(IN):: X_IN
    REAL(DP), INTENT(INOUT) :: X(6)
    INTEGER, INTENT(IN):: i
    INTEGER J

    IF(I<=SIZE(X_IN%RAY,2)) THEN
       DO J=1,6
          X(J)=  X_IN%RAY(J,I)
       ENDDO
    ELSE
       WRITE(6,*) I
       STOP 10
    ENDIF
  END SUBROUTINE gMIDR



  SUBROUTINE ALLOC_midr(X_IN,R)
    IMPLICIT NONE
    TYPE(worm), INTENT(INOUT):: X_IN
    TYPE(LAYOUT), INTENT(IN):: R
    INTEGER I
    TYPE(FIBRE), POINTER:: P

    P=>R%START
    allocate(x_in%nst)
    X_IN%NST=3
    DO I=1,R%N
       X_IN%nst= MAX(P%MAG%p%NST,X_IN%NST)
       P=>P%NEXT
    ENDDO

    allocate(x_in%RAY(6,-6:X_IN%nst+6))
    allocate(x_in%E)
    allocate(x_in%E%L(-1:X_IN%nst))
    allocate(x_in%POS(4))

    allocate(x_in%E%FRAME(3,3,-7:X_IN%nst+6))
    allocate(x_in%E%ORIGIN(3,-7:X_IN%nst+6))
    ALLOCATE(x_in%E%DO_SURVEY)
    x_in%E%DO_SURVEY=.TRUE.
    x_in%nst=0
    x_in%POS=0
    x_in%RAY=0.0_dp
    x_in%E%L=0.0_dp
    x_in%E%FRAME=0.0_dp
    x_in%E%ORIGIN=0.0_dp
    x_in%e%nst=>x_in%nst

  END SUBROUTINE ALLOC_midr

  SUBROUTINE ALLOC_FIBRE(X_IN,P)
    IMPLICIT NONE
    TYPE(worm), INTENT(INOUT):: X_IN
    TYPE(FIBRE),TARGET, INTENT(INOUT):: P

    allocate(x_in%nst)
    X_IN%NST=3
    X_IN%nst= MAX(P%MAG%p%NST,X_IN%NST)


    allocate(x_in%RAY(6,-6:X_IN%nst+6))
    allocate(x_in%E)
    allocate(x_in%E%L(-1:X_IN%nst))
    allocate(x_in%POS(4))

    allocate(x_in%E%FRAME(3,3,-7:X_IN%nst+6))
    allocate(x_in%E%ORIGIN(3,-7:X_IN%nst+6))
    ALLOCATE(x_in%E%DO_SURVEY)
    x_in%E%DO_SURVEY=.TRUE.
    x_in%nst=0
    x_in%POS=0
    x_in%RAY=0.0_dp
    x_in%E%L=0.0_dp
    x_in%E%FRAME=0.0_dp
    x_in%E%ORIGIN=0.0_dp
    x_in%e%nst=>x_in%nst

  END SUBROUTINE ALLOC_FIBRE

  SUBROUTINE KILL_midr(X_IN)
    IMPLICIT NONE
    TYPE(worm), INTENT(INOUT):: X_IN

    DEallocate(x_in%nst)

    DEallocate(x_in%RAY)
    DEallocate(x_in%POS)

    DEallocate(x_in%E%FRAME)
    DEallocate(x_in%E%ORIGIN)
    DEallocate(x_in%E%L)
    DEallocate(x_in%E%DO_SURVEY)
    DEallocate(x_in%E)

  END SUBROUTINE KILL_midr


  SUBROUTINE SURVEY_CHART(C,P,DIR,MAGNETFRAME,E_IN)
    !changed
    ! SURVEYS A SINGLE ELEMENT FILLS IN CHART AND MAGNET_CHART; LOCATES ORIGIN AT THE ENTRANCE OR EXIT
    IMPLICIT NONE
    TYPE(CHART), TARGET ,OPTIONAL, INTENT(INOUT):: C
    TYPE(MAGNET_CHART), TARGET,INTENT(INOUT) :: P
    TYPE (CHART), POINTER :: CL
    TYPE(MAGNET_FRAME), TARGET, OPTIONAL :: MAGNETFRAME
    TYPE(INNER_FRAME), OPTIONAL :: E_IN
    INTEGER, INTENT(IN) ::DIR
    TYPE(MAGNET_FRAME), POINTER :: F
    REAL(DP) ENT(3,3),EXI(3,3),HA,D(3),BASIS(3,3),OMEGA(3),A(3),N(3)
    INTEGER I,J
    CALL ALLOC(F)


    CL=> C  ! CHART OF ELEMENT 1

    HA=DIR*P%LD*P%B0/2.0_dp
    D=0.0_dp
    D(3)=DIR*P%LC/2.0_dp
    IF(ASSOCIATED(CL%F)) THEN     !!!! DOING SURVEY
       IF(DIR==1) THEN
          A=0.0_dp;A(3)=P%TILTD  ;
          CALL GEO_ROT(CL%F%ENT,ENT      ,A  ,CL%F%ENT)
          IF(PRESENT(E_IN) ) THEN

             CALL XFRAME(E_IN,ENT,CL%F%A,-2)
             !          WRITE(6,*) "E_IN%NST ",E_IN%NST
             !          WRITE(6,*) ENT
             !          WRITE(6,*) E_IN%FRAME(:,:,-2)

             !         PAUSE 123
          ENDIF
          A=0.0_dp;A(2)=HA ;
          CALL GEO_ROT(ENT     ,CL%F%MID ,A     ,ENT)
          CALL GEO_ROT(CL%F%MID,EXI     , A     ,CL%F%MID)

          IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT=EXI,I=E_IN%nst-4)
          A=0.0_dp;A(3)=-P%TILTD  ;
          CALL GEO_ROT(EXI     ,CL%F%EXI ,A,EXI)

          CL%F%O=CL%F%A
          CALL GEO_TRA(CL%F%O,CL%F%MID,D,1)
          CL%F%B=CL%F%O
          CALL GEO_TRA(CL%F%B,CL%F%MID,D,1)

          IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,A=CL%F%B,I=E_IN%nst-4)

       ELSE
          A=0.0_dp;A(3)=P%TILTD  ;
          CALL GEO_ROT(CL%F%EXI,EXI      ,A  ,CL%F%EXI)
          IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,EXI,CL%F%B,-2)
          A=0.0_dp;A(2)=HA ;
          CALL GEO_ROT(EXI     ,CL%F%MID ,A     ,EXI)
          CALL GEO_ROT(CL%F%MID,ENT     , A     ,CL%F%MID)

          IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,ENT=ENT,I=E_IN%nst-4)
          A=0.0_dp;A(3)=-P%TILTD  ;
          CALL GEO_ROT(ENT     ,CL%F%ENT ,A ,ENT)

          CL%F%O=CL%F%B
          CALL GEO_TRA(CL%F%O,CL%F%MID,D,1)
          CL%F%A=CL%F%O
          CALL GEO_TRA(CL%F%A,CL%F%MID,D,1)

          IF(PRESENT(E_IN) ) CALL XFRAME(E_IN,A=CL%F%A,I=E_IN%nst-4)



       ENDIF

       !  CAN BE THE SAME FOR DIR =1 OR DIR= -1
       F=CL%F

       F%ENT=ENT
       F%EXI=EXI
       BASIS=F%ENT
       OMEGA=F%A
       IF(PRESENT(MAGNETFRAME) )THEN
          MAGNETFRAME=F
       ENDIF

       !          DO I=1,3
       !          A=ZERO
       !          A(I)=C%ANG_IN(I)

       A=C%ANG_IN
       call ROTATE_FRAME(f,OMEGA,a,1,BASIS)

       !       D=F%A-OMEGA
       !       CALL GEO_ROT(F%ENT,D,A,1,BASIS)
       !       F%A=OMEGA+D

       !       D=F%O-OMEGA
       !       CALL GEO_ROT(F%MID,D,A,1,BASIS)
       !       F%O=OMEGA+D

       !       D=F%B-OMEGA
       !       CALL GEO_ROT(F%EXI,D,A,1,BASIS)
       !       F%B=OMEGA+D

       BASIS=F%ENT

       !          ENDDO

       CALL GEO_TRA(F%A,BASIS,C%D_IN,1)
       CALL GEO_TRA(F%B,BASIS,C%D_IN,1)
       CALL GEO_TRA(F%O,BASIS,C%D_IN,1)

       IF(PRESENT(E_IN) ) THEN
          IF(DIR==1) THEN
             CALL XFRAME(E_IN,F%ENT,F%A,-1)
             CALL XFRAME(E_IN,F%EXI,F%B,E_IN%nst-5)
          ELSE
             CALL XFRAME(E_IN,F%EXI,F%B,-1)
             CALL XFRAME(E_IN,F%ENT,F%A,E_IN%nst-5)
          ENDIF

          call SURVEY_inner_mag(E_IN)
       ENDIF


       ! CHECKING HERE THE CONSISTANCY
       ENT=F%EXI
       BASIS=F%EXI
       A=C%ANG_OUT
       OMEGA=F%B
       D=F%B-OMEGA   ! D=0 OF COURSE
       CALL GEO_ROT(ENT,D,A,1,BASIS)
       OMEGA=OMEGA+D

       CALL GEO_TRA(OMEGA,ENT,C%D_OUT,1)
       ENT=ENT-EXI
       OMEGA=OMEGA-CL%F%B
       N=0.0_dp
       DO I=1,3
          N(2)=ABS(OMEGA(I))+N(2)
          DO J=1,3
             N(1)=ABS(ENT(I,J))+N(1)
          ENDDO
       ENDDO

       !       IF(N(2)<EPS_FITTED) N(2)=N(2)/( ABS(CL%F%B(1))+ABS(CL%F%B(2))+ABS(CL%F%B(3)) )

       IF(N(1)>EPS_FITTED.OR.N(2)>EPS_FITTED) THEN
          WRITE(6,*) "INCONSISTANCY IN SURVEY_CHART "
          WRITE(6,*) N(1),N(2)
       ENDIF
       !


       IF(ASSOCIATED(P%F)) THEN
          P%F=F
       ENDIF



    ENDIF  !!!! DOING SURVEY


    CALL KILL(F)
    IF(ASSOCIATED(F)) deallocate(f)

  END SUBROUTINE SURVEY_CHART

  SUBROUTINE SURVEY_INNER_MAG(e_in) !  Tracks the chart through a magnet
    IMPLICIT NONE
    TYPE(INNER_FRAME), INTENT(INOUT):: e_in
    REAL(DP) ENT(3,3),A(3),MID(3,3),O(3),D(3)
    LOGICAL(LP) DONE
    INTEGER NST,I,start
    REAL(DP) LH,HA,ANG(3),ANGH,RHO
    TYPE(MAGNET_CHART), POINTER :: P

    NST=E_IN%NST-6
    DONE=.FALSE.
    start=nst*(1-(1+E_IN%F%dir)/2)


    P=>E_IN%F%MAG%P
    !E_IN%L

    IF(E_IN%DO_SURVEY) THEN   !  DOING THE SURVEY

       !     CALL GFRAME(E_IN,ENT,A,-1)
       IF(ASSOCIATED(E_IN%F%CHART)) THEN
          IF(ASSOCIATED(E_IN%F%CHART%F)) THEN
             MID=E_IN%F%CHART%F%MID
             O=E_IN%F%CHART%F%O
             DONE=.TRUE.
          ENDIF
       ENDIF

       IF(ASSOCIATED(E_IN%F%MAG%P%F)) THEN
          MID= P%F%MID
          O  = P%F%O
          DONE=.TRUE.
       ENDIF

       IF(.NOT.DONE) THEN
          WRITE(6,*) "ERROR IN SURVEY_INNER_MAG, NO FRAME WHATSOEVER "
          STOP 330
       ENDIF

       SELECT CASE(E_IN%F%MAG%KIND)

          !       CASE(KIND0)
          !          CALL XFRAME(E_IN,MID,O,0)
          !          E_IN%L(0)=zero +E_IN%L(-1)
          !          IF(NST/=0) THEN
          !             WRITE(6,*) "ERROR IN SURVEY_INNER_MAG at kind23"
          !             STOP 330
          !          ENDIF

          !       CASE(kind23)                 ! kind 23 layout
          !          call  GET_LENGTH(E_IN%F%mag%g23,Lh)
          !
          !          E_IN%L(start)=start*lh/nst  +E_IN%L(-1)
          !
          !          if(E_IN%F%dir==1) then
          !             CALL XFRAME(E_IN,P%F%ent,P%F%a,start)
          !          else
          !             CALL XFRAME(E_IN,P%F%exi,P%F%b,start)
          !          endif
          !
          !          start=start+E_IN%F%dir
          !
          !          E_IN%L(start)=start*lh/nst  +E_IN%L(-1)
          !
          !          if(E_IN%F%dir==1) then
          !             CALL XFRAME(E_IN,P%F%exi,P%F%b,start)
          !          else
          !             CALL XFRAME(E_IN,P%F%ent,P%F%a,start)
          !          endif

          !          IF(NST/=1) THEN
          !             WRITE(6,*) "ERROR IN SURVEY_INNER_MAG "
          !             STOP 331
          !          ENDIF
       CASE(KIND0,KIND1,KIND3:KIND5,KIND8:KIND9,KIND11:KIND15,KIND17:KIND22,kindwiggler)
          LH=P%LC/2.0_dp
          A=O
          D=0.0_dp;D(3)=-LH
          CALL GEO_TRA(A,MID,D,1)
          CALL XFRAME(E_IN,MID,A,start)

          HA=P%LC/NST
          E_IN%L(start)=start*P%LD/nst  +E_IN%L(-1)
          D=0.0_dp;D(3)=HA
          DO I=1,NST
             start=start+E_IN%F%dir
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
             CALL GEO_TRA(A,MID,D,1)
             CALL XFRAME(E_IN,MID,A,start)
          ENDDO
       CASE(KIND2,KIND6:KIND7,KIND10,KINDPA)
          IF(P%B0==0.0_dp) THEN
             LH=P%LC/2.0_dp
             A=O
             D=0.0_dp;D(3)=-LH
             CALL GEO_TRA(A,MID,D,1)
             CALL XFRAME(E_IN,MID,A,start)
             HA=P%LC/NST
             D=0.0_dp;D(3)=HA
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
             DO I=1,NST
                start=start+E_IN%F%dir
                E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
                CALL GEO_TRA(A,MID,D,1)
                CALL XFRAME(E_IN,MID,A,start)
             ENDDO
          ELSE
             RHO=1.0_dp/P%B0
             ANG=0.0_dp; D=0.0_dp;
             LH=P%LC/2.0_dp
             A=O
             D(3)=-LH
             ANGH=P%LD*P%B0/2.0_dp
             ANG(2)=-ANGH
             CALL GEO_TRA(A,MID,D,1)
             O=A
             CALL GEO_ROT(MID,ENT      ,ANG  ,MID)
             CALL XFRAME(E_IN,ENT,A,start)
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)

             ANG(2)=2.0_dp*ANGH/NST
             DO I=1,NST
                start=start+E_IN%F%dir
                E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
                HA=ANGH-I*ANG(2)
                CALL GEO_ROT(ENT,ENT      ,ANG  ,MID)
                D=0.0_dp
                D(1)=RHO*(COS(ha)-COS(ANGH))
                D(3)=P%LC/2.0_dp-sin(ha)*rho
                A=O
                CALL GEO_TRA(A,MID,D,1)
                CALL XFRAME(E_IN,ENT,A,start)
             ENDDO

          ENDIF

       CASE(KIND16)
          ANGH=P%LD*P%B0/2.0_dp
          LH=P%LC/2.0_dp
          A=O
          D=0.0_dp;D(3)=-LH
          CALL GEO_TRA(A,MID,D,1)
          ANG=0.0_dp;ANG(2)=-(ANGH-P%EDGE(1))
          CALL GEO_ROT(MID,MID      ,ANG  ,MID)

          CALL XFRAME(E_IN,MID,A,start)
          HA=E_IN%F%MAG%L/NST
          E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
          D=0.0_dp;D(3)=HA
          DO I=1,NST
             start=start+E_IN%F%dir
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
             CALL GEO_TRA(A,MID,D,1)
             CALL XFRAME(E_IN,MID,A,start)
          ENDDO


       CASE DEFAULT



          Write(6,*)"KIND = ",  E_IN%F%MAG%KIND," NOT SUPPORTED IN SURVEY_INNER_MAG "
          STOP 778

       END  SELECT



    ELSE    ! NOT  DOING THE SURVEY




       SELECT CASE(E_IN%F%MAG%KIND)

          !       CASE(KIND0)
          !          E_IN%L(0)=zero +E_IN%L(-1)
          !          IF(NST/=0) THEN
          !             WRITE(6,*) "ERROR IN SURVEY_INNER_MAG "
          !             STOP 330
          !          ENDIF
       CASE(KIND0,KIND1,KIND3:KIND5,KIND8:KIND9,KIND11:KIND15,KIND17:KIND22,kindwiggler)
          E_IN%L(start)=start*P%LD/nst  +E_IN%L(-1)
          DO I=1,NST
             start=start+E_IN%F%dir
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
          ENDDO
       CASE(KIND2,KIND6:KIND7,KIND10,KINDPA)
          IF(P%B0==0.0_dp) THEN
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
             DO I=1,NST
                start=start+E_IN%F%dir
                E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
             ENDDO
          ELSE
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)

             DO I=1,NST
                start=start+E_IN%F%dir
                E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
             ENDDO

          ENDIF

       CASE(KIND16)
          E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
          DO I=1,NST
             start=start+E_IN%F%dir
             E_IN%L(start)=start*P%LD/nst    +E_IN%L(-1)
          ENDDO

          !       CASE(kind23)                 ! kind 23 layout
          !          call  GET_LENGTH(E_IN%F%mag%g23,Lh)
          !
          !          E_IN%L(start)=start*lh/nst  +E_IN%L(-1)
          !
          !
          !          start=start+E_IN%F%dir
          !
          !          E_IN%L(start)=start*lh/nst  +E_IN%L(-1)

          !          IF(NST/=1) THEN
          !             WRITE(6,*) "ERROR IN SURVEY_INNER_MAG "
          !             STOP 331
          !          ENDIF
          !
       CASE DEFAULT



          Write(6,*)"KIND = ",  E_IN%F%MAG%KIND," NOT SUPPORTED IN SURVEY_INNER_MAG "
          STOP 778

       END  SELECT


    ENDIF



    E_IN%L(-1)=E_IN%L(-1)+P%Ld




  end SUBROUTINE SURVEY_INNER_MAG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   New Survey Routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine survey_integration_node_case0(a0,ent0,t,b0,exi0)
implicit none
type(integration_node), target :: t
real(dp), target :: a0(3),ent0(3,3),b0(3),exi0(3,3),ang(3)
type(fibre), pointer :: f
type(element), pointer :: m
type(magnet_chart), pointer :: p
real(dp) h,d(3)

f=>t%parent_fibre
m=>f%mag
p=>m%p
b0=a0
exi0=ent0
write(6,*) " kind ",kind2,m%kind
select case(m%kind) 

CASE(KIND0,KIND1,KIND3:KIND5,KIND8:KIND9,KIND11:KIND15,KIND17:KIND22,kindwiggler,kindsuper1)
   h=p%lc/p%nst
   d=(/0.0_dp,0.0_dp,h/)

   call geo_tra(b0,exi0,d,1)
CASE(KIND2,KIND6:KIND7,KIND10,KINDPA)
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

CASE(KIND16)
   h=m%l/p%nst
   d=(/0.0_dp,0.0_dp,h/)
   call geo_tra(b0,exi0,d,1)

CASE default
 write(6,*) " not supported in survey_integration_node_case0 "
 stop

end select

end subroutine survey_integration_node_case0

subroutine survey_integration_node_p1(a0,ent0,t,b0,exi0)
implicit none
type(integration_node), target :: t
real(dp), target :: a0(3),ent0(3,3),b0(3),exi0(3,3)
type(fibre), pointer :: f
real(dp) pix1(3),pix2(3) 
logical(lp) :: ENTERING=my_true

f=>t%parent_fibre
pix1=0.0_dp;pix2=0.0_dp;
b0=a0
exi0=ent0



if(f%patch%A_X1==-1) pix1(1)=pi
if(f%patch%A_X2==-1) pix2(1)=pi

!
call GEO_ROT(exi0,pix1,1, ent0)
call GEO_ROT(exi0,f%patch%a_ang,1, exi0)
call TRANSLATE_point(b0,f%patch%A_D,1,exi0)  
call GEO_ROT(exi0,pix2,1, exi0)

pix1=0.0_dp
pix1(3)=f%MAG%P%TILTD
 call GEO_ROT(exi0,pix1,1, exi0)

    IF(f%MAG%MIS) THEN
      call MIS_survey(b0,exi0,f,b0,exi0,ENTERING)
    ENDIF

!    IF(PATCHG==1.or.PATCHG==3) THEN
!       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
!       CALL PATCH_FIB(C,X,k,PATCH,MY_TRUE)
!    ENDIF
!       X(3)=C%PATCH%A_X1*X(3);X(4)=C%PATCH%A_X1*X(4);
!       CALL ROT_YZ(C%PATCH%A_ANG(1),X,C%MAG%P%BETA0,PATCH,k%TIME)
!       CALL ROT_XZ(C%PATCH%A_ANG(2),X,C%MAG%P%BETA0,PATCH,k%TIME)
!       CALL ROT_XY(C%PATCH%A_ANG(3),X)  !,PATCH)
!       CALL TRANS(C%PATCH%A_D,X,C%MAG%P%BETA0,PATCH,k%TIME)
!       X(3)=C%PATCH%A_X2*X(3);X(4)=C%PATCH%A_X2*X(4);

end subroutine survey_integration_node_p1


subroutine survey_integration_node_p2(a0,ent0,t,b0,exi0)
implicit none
type(integration_node), target :: t
real(dp), target :: a0(3),ent0(3,3),b0(3),exi0(3,3)
type(fibre), pointer :: f
real(dp) pix1(3),pix2(3)
logical(lp) :: ENTERING=my_FALSE

f=>t%parent_fibre

b0=a0
exi0=ent0




    IF(f%MAG%MIS) THEN
      call MIS_survey(b0,exi0,f,b0,exi0,ENTERING)
    ENDIF


pix1=0.0_dp
pix1(3)=-f%MAG%P%TILTD
 call GEO_ROT(exi0,pix1,1, exi0)

pix1=0.0_dp;pix2=0.0_dp;
if(f%patch%B_X1==-1) pix1(1)=pi
if(f%patch%B_X2==-1) pix2(1)=pi

!
call GEO_ROT(exi0,pix1,1, ent0)
call GEO_ROT(exi0,f%patch%B_ang,1, exi0)
call TRANSLATE_point(b0,f%patch%b_D,1,exi0)  
call GEO_ROT(exi0,pix2,1, exi0)


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


end module S_def_all_kinds
