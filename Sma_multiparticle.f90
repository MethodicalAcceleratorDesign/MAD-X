module ptc_multiparticle
  use S_TRACKING  !,FRINGE_=>FRINGE__MULTI !,FACE=>FACE_MULTI

  implicit none
  public
  private index
  PRIVATE INTER_CAV_TRAV,ADJUST_WI
  PRIVATE INTER_DRIFT1,INTER_dkd2,INTER_KICKT3,INTER_CAV4,INTER_SOL5,INTER_KTK,INTER_TEAPOT &
       ,INTER_STREX,INTER_SOLT,INTER_NSMI,INTER_SSMI,INTER_MON_11_14,INTER_ESEPTUM,INTER_RCOL18 &
       ,INTER_ECOL19,INTER_WI,INTER_TKTF
  PRIVATE INTEP_DRIFT1,INTEP_dkd2,INTEP_KICKT3,INTEP_CAV4,INTEP_SOL5,INTEP_KTK,INTEP_TEAPOT &
       ,INTEP_STREX,INTEP_SOLT,INTEP_NSMI,INTEP_SSMI,INTEP_MON_11_14,INTEP_ESEPTUM,INTEP_RCOL18 &
       ,INTEP_ECOL19,INTEP_WI,INTEp_TKTF
  PRIVATE TRACK_NODE_LAYOUT_S1,TRACK_NODE_LAYOUT_S12,TRACK_LAYOUT_ONE_12,TRACK_NODE_LAYOUT_S
  INTEGER, PARAMETER :: CASE1=1,CASE2=2, CASE0=0, CASEP1=-1,CASEP2=-2  !,CASE3=3
  INTEGER, PRIVATE :: TOTALPATH_FLAG
  integer :: index =0
  CHARACTER*27 CASE_NAME(-2:3)
  PRIVATE fuzzy_eq,fuzzy_neq
  PRIVATE ADJUSTR_WI,ADJUSTP_WI,TRACK_FIBRE_FRONTR,TRACK_FIBRE_FRONTP
  PRIVATE TRACK_FIBRE_BACKR,TRACK_FIBRE_BACKP
  PRIVATE FRINGER_CAV4,FRINGEP_CAV4,fringeR_STRAIGHT,fringeP_STRAIGHT
  PRIVATE fringer_STREX,fringep_STREX,fringer_TEAPOT,fringeP_TEAPOT
  PRIVATE FRINGER_CAV_TRAV,FRINGEP_CAV_TRAV
  PRIVATE TRACKR_NODE_SINGLE,TRACKP_NODE_SINGLE,TRACKV_NODE_SINGLE
  PRIVATE ADJUST_TIME_CAV4,ADJUSTR_TIME_CAV4,ADJUSTP_TIME_CAV4
  private DRIFTr_BACK_TO_POSITION,DRIFTp_BACK_TO_POSITION,DRIFT_BACK_TO_POSITION,DRIFT_BEAM_BACK_TO_POSITION
  PRIVATE TRACK_LAYOUT_XR_12,TRACK_LAYOUT_XP_12
  PRIVATE ADJUSTR_PANCAKE,ADJUSTP_PANCAKE,ADJUST_PANCAKE,INTER_PANCAKE,INTEP_PANCAKE
  private MAKE_NODE_LAYOUT_2,DRIFT_TO_TIME
  PRIVATE TRACK_NODE_LAYOUT_FLAG_R

  INTERFACE ADJUST_WI
     MODULE PROCEDURE ADJUSTR_WI
     MODULE PROCEDURE ADJUSTP_WI
  END INTERFACE

  INTERFACE DRIFT_BACK_TO_POSITION
     MODULE PROCEDURE DRIFTr_BACK_TO_POSITION
     MODULE PROCEDURE DRIFTp_BACK_TO_POSITION
  END INTERFACE

  INTERFACE ADJUST_TIME_CAV4
     MODULE PROCEDURE ADJUSTR_TIME_CAV4
     MODULE PROCEDURE ADJUSTP_TIME_CAV4
  END INTERFACE

  INTERFACE ADJUST_TIME_CAV_TRAV_OUT
     MODULE PROCEDURE ADJUSTR_TIME_CAV_TRAV_OUT
     MODULE PROCEDURE ADJUSTP_TIME_CAV_TRAV_OUT
  END INTERFACE

  INTERFACE ADJUST_PANCAKE
     MODULE PROCEDURE ADJUSTR_PANCAKE
     MODULE PROCEDURE ADJUSTP_PANCAKE
  END INTERFACE



  INTERFACE TRACK_FIBRE_FRONT
     MODULE PROCEDURE TRACK_FIBRE_FRONTR
     MODULE PROCEDURE TRACK_FIBRE_FRONTP
  END INTERFACE

  INTERFACE TRACK_FIBRE_BACK
     MODULE PROCEDURE TRACK_FIBRE_BACKR
     MODULE PROCEDURE TRACK_FIBRE_BACKP
  END INTERFACE


  INTERFACE TRACK_SLICE
     MODULE PROCEDURE INTER_DRIFT1
     MODULE PROCEDURE INTER_dkd2
     MODULE PROCEDURE INTER_KICKT3
     MODULE PROCEDURE INTER_CAV4
     MODULE PROCEDURE INTER_SOL5
     MODULE PROCEDURE INTER_KTK
     MODULE PROCEDURE INTER_TKTF
     MODULE PROCEDURE INTER_NSMI
     MODULE PROCEDURE INTER_SSMI
     MODULE PROCEDURE INTER_TEAPOT
     MODULE PROCEDURE INTER_STREX
     MODULE PROCEDURE INTER_SOLT
     MODULE PROCEDURE INTER_MON_11_14
     MODULE PROCEDURE INTER_ESEPTUM
     MODULE PROCEDURE INTER_CAV_TRAV
     MODULE PROCEDURE INTER_RCOL18
     MODULE PROCEDURE INTER_ECOL19
     MODULE PROCEDURE INTER_WI
     MODULE PROCEDURE INTER_PANCAKE

     MODULE PROCEDURE INTEP_DRIFT1
     MODULE PROCEDURE INTEP_dkd2
     MODULE PROCEDURE INTEP_KICKT3
     MODULE PROCEDURE INTEP_CAV4
     MODULE PROCEDURE INTEP_SOL5
     MODULE PROCEDURE INTEP_KTK
     MODULE PROCEDURE INTEP_TKTF
     MODULE PROCEDURE INTEP_NSMI
     MODULE PROCEDURE INTEP_SSMI
     MODULE PROCEDURE INTEP_TEAPOT
     MODULE PROCEDURE INTEP_STREX
     MODULE PROCEDURE INTEP_SOLT
     MODULE PROCEDURE INTEP_MON_11_14
     MODULE PROCEDURE INTEP_ESEPTUM
     MODULE PROCEDURE INTEP_CAV_TRAV
     MODULE PROCEDURE INTEP_RCOL18
     MODULE PROCEDURE INTEP_ECOL19
     MODULE PROCEDURE INTEP_WI
     MODULE PROCEDURE INTEP_PANCAKE
  END INTERFACE

  INTERFACE TRACK_FRINGE
     MODULE PROCEDURE fringeR_STRAIGHT   ! GOOD FOR DKD2,KTK, AND TKTF
     MODULE PROCEDURE fringeP_STRAIGHT   ! GOOD FOR DKD2,KTK, AND TKTF
     MODULE PROCEDURE FRINGER_CAV4
     MODULE PROCEDURE FRINGEP_CAV4
     MODULE PROCEDURE fringer_TEAPOT
     MODULE PROCEDURE fringeP_TEAPOT
     MODULE PROCEDURE fringer_strex
     MODULE PROCEDURE fringeP_STREX
     MODULE PROCEDURE FRINGER_CAV_TRAV
     MODULE PROCEDURE FRINGEP_CAV_TRAV
  END INTERFACE

  INTERFACE TRACK_NODE_SINGLE
     MODULE PROCEDURE TRACKR_NODE_SINGLE
     MODULE PROCEDURE TRACKP_NODE_SINGLE
     MODULE PROCEDURE TRACKV_NODE_SINGLE
  END INTERFACE

  !  INTERFACE TRACK
  !     MODULE PROCEDURE TRACKR_NODE_SINGLE
  !     MODULE PROCEDURE TRACKP_NODE_SINGLE
  !  END INTERFACE

  INTERFACE TRACK_NODE_LAYOUT_S
     MODULE PROCEDURE TRACK_NODE_LAYOUT_S12
     MODULE PROCEDURE TRACK_NODE_LAYOUT_S1
  END INTERFACE

  INTERFACE TRACK_NODE_LAYOUT
     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_R
     !     MODULE PROCEDURE TRACK_NODE_LAYOUT_FLAG_P
  END INTERFACE


  INTERFACE TRACK_LAYOUT_USING_NODE_S
     MODULE PROCEDURE TRACK_LAYOUT_XV_12
     MODULE PROCEDURE TRACK_LAYOUT_XR_12
     MODULE PROCEDURE TRACK_LAYOUT_XP_12
     MODULE PROCEDURE TRACK_LAYOUT_ONE_12
     !     MODULE PROCEDURE TRACK_LAYOUT_S12
  END INTERFACE

  INTERFACE TRACK
     MODULE PROCEDURE TRACK_LAYOUT_XV_12
     MODULE PROCEDURE TRACK_LAYOUT_XR_12
     MODULE PROCEDURE TRACK_LAYOUT_XP_12
     MODULE PROCEDURE TRACK_LAYOUT_ONE_12
     !     MODULE PROCEDURE TRACK_LAYOUT_S12
     MODULE PROCEDURE TRACK_NODE_T
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

  real :: ttime0,ttime1,dt1=0.0,dt2=0.0;

CONTAINS

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

    if(sp==zero.and.s/=zero) then
       current=>l%end
       i=l%n+1
       ds=zero
       !  track_it=.true.
       return
    endif

    if(sp==zero) then
       current=>l%start
       i=1
       ds=zero
       return
    endif


    nullify(current);
    Current => L%LAST

    k=L%LASTPOS
    I=K
    ds=zero

    IF(SP>CURRENT%S(3) ) then

       do i=k,l%n-1
          if(current%next%s(3)>=sp) exit
          current=>current%next
       enddo

       if(current%next%s(3)/=sp) ds=sp-current%s(3)
       if(current%next%s(3)==sp) THEN
          ds=ZERO
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

    if(ds>zero) then
       if(CURRENT%S(4)-ds.feq.zero) then
          ds=zero
          current=>Current%next
          i=i+1
          L%LAST => Current;
       ELSEIF(ds.feq.zero) THEN
          DS=ZERO
       ENDIF
    endif


    DOIT=.TRUE.
    !DOIT=.FALSE.

    if(iabs(CURRENT%cas)==0.OR.iabs(CURRENT%cas)==1) then
       do while(DS==ZERO.AND.DOIT)    ! PUTS AT BEGINNING IF DS=ZERO
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
       do while(DS==ZERO.AND.DOIT)    ! PUTS AT BEGINNING IF DS=ZERO
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



  ! tracking fringe areas

  SUBROUTINE fringer_STRAIGHT(EL,EL5,EL6,EL7,EL17,X,J)
    IMPLICIT NONE
    TYPE(DKD2),OPTIONAL,INTENT(IN):: EL
    TYPE(SOL5),OPTIONAL,INTENT(INOUT):: EL5
    TYPE(KTK),OPTIONAL,INTENT(INOUT):: EL6
    TYPE(TKTF),OPTIONAL,INTENT(INOUT):: EL7
    TYPE(SOLT),OPTIONAL,INTENT(INOUT):: EL17
    !      TYPE(BEAM), INTENT(INOUT) ::B
    integer,INTENT(IN):: J
    real(dp), INTENT(INOUT) :: X(6)
    integer  i
    ! J=1 front
    IF(PRESENT(EL)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL%P%DIR==1) THEN
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ELSE
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF
       else
          if(EL%P%DIR==1) THEN
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
          ENDIF
       ENDIF
       !         CALL X_IN_BEAM(X,B,I)

    ELSEIF(PRESENT(EL6)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL6%P%DIR==1) THEN
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,1,X)
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,1,X)
          ELSE
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,2,X)
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,2,X)
          ENDIF
       else
          if(EL6%P%DIR==1) THEN
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,2,X)
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,2,X)
          ELSE
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,1,X)
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,1,X)
          ENDIF
       ENDIF
       !         CALL X_IN_BEAM(X,B,I)

    ELSEIF(PRESENT(EL5)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL5%P%DIR==1) THEN
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,1,X)
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,1,X)
          ELSE
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,2,X)
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,2,X)
          ENDIF
       else
          if(EL5%P%DIR==1) THEN
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,2,X)
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,2,X)
          ELSE
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,1,X)
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,1,X)
          ENDIF
       ENDIF

    ELSEIF(PRESENT(EL7)) THEN

       if(J==1) then
          if(EL7%P%DIR==1) THEN
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,1,X)
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,1,X)
          ELSE
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,2,X)
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,2,X)
          ENDIF
       else
          if(EL7%P%DIR==1) THEN
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,2,X)
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,2,X)
          ELSE
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,1,X)
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,1,X)
          ENDIF
       ENDIF

    ELSEIF(PRESENT(EL17)) THEN

       if(J==1) then
          if(EL17%P%DIR==1) THEN
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,1,X)
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,1,X)
          ELSE
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,2,X)
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,2,X)
          ENDIF
       else
          if(EL17%P%DIR==1) THEN
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,2,X)
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,2,X)
          ELSE
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,1,X)
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,1,X)
          ENDIF
       ENDIF

    ENDIF
  END SUBROUTINE fringer_STRAIGHT


  SUBROUTINE fringeP_STRAIGHT(EL,EL5,EL6,EL7,EL17,X,J)
    IMPLICIT NONE
    TYPE(DKD2P),OPTIONAL,INTENT(IN):: EL
    TYPE(SOL5P),OPTIONAL,INTENT(INOUT):: EL5
    TYPE(KTKP),OPTIONAL,INTENT(INOUT):: EL6
    TYPE(TKTFP),OPTIONAL,INTENT(INOUT):: EL7
    TYPE(SOLTP),OPTIONAL,INTENT(INOUT):: EL17
    !      TYPE(BEAM), INTENT(INOUT) ::B
    integer,INTENT(IN):: J
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    integer  i
    ! J=1 front
    IF(PRESENT(EL)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL%P%DIR==1) THEN
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ELSE
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF
       else
          if(EL%P%DIR==1) THEN
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL EDGE(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
          ENDIF
       ENDIF
       !         CALL X_IN_BEAM(X,B,I)

    ELSEIF(PRESENT(EL6)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL6%P%DIR==1) THEN
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,1,X)
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,1,X)
          ELSE
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,2,X)
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,2,X)
          ENDIF
       else
          if(EL6%P%DIR==1) THEN
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,2,X)
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,2,X)
          ELSE
             IF(EL6%P%FRINGE) CALL MULTIPOLE_FRINGE(EL6%P,EL6%AN,EL6%BN,1,X)
             CALL EDGE(EL6%P,EL6%BN,EL6%H1,EL6%H2,EL6%FINT,EL6%HGAP,1,X)
          ENDIF
       ENDIF
       !         CALL X_IN_BEAM(X,B,I)

    ELSEIF(PRESENT(EL5)) THEN
       !       DO I=1,B%N
       !        IF(B%U(i)) CYCLE
       !         X=BEAM_IN_X(B,I)

       if(J==1) then
          if(EL5%P%DIR==1) THEN
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,1,X)
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,1,X)
          ELSE
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,2,X)
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,2,X)
          ENDIF
       else
          if(EL5%P%DIR==1) THEN
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,2,X)
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,2,X)
          ELSE
             IF(EL5%P%FRINGE) CALL MULTIPOLE_FRINGE(EL5%P,EL5%AN,EL5%BN,1,X)
             CALL EDGE(EL5%P,EL5%BN,EL5%H1,EL5%H2,EL5%FINT,EL5%HGAP,1,X)
          ENDIF
       ENDIF

    ELSEIF(PRESENT(EL7)) THEN

       if(J==1) then
          if(EL7%P%DIR==1) THEN
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,1,X)
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,1,X)
          ELSE
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,2,X)
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,2,X)
          ENDIF
       else
          if(EL7%P%DIR==1) THEN
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,2,X)
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,2,X)
          ELSE
             IF(EL7%P%FRINGE) CALL MULTIPOLE_FRINGE(EL7%P,EL7%AN,EL7%BN,1,X)
             CALL EDGE(EL7%P,EL7%BN,EL7%H1,EL7%H2,EL7%FINT,EL7%HGAP,1,X)
          ENDIF
       ENDIF

    ELSEIF(PRESENT(EL17)) THEN

       if(J==1) then
          if(EL17%P%DIR==1) THEN
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,1,X)
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,1,X)
          ELSE
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,2,X)
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,2,X)
          ENDIF
       else
          if(EL17%P%DIR==1) THEN
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,2,X)
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,2,X)
          ELSE
             IF(EL17%P%FRINGE) CALL MULTIPOLE_FRINGE(EL17%P,EL17%AN,EL17%BN,1,X)
             CALL EDGE(EL17%P,EL17%BN,EL17%H1,EL17%H2,EL17%FINT,EL17%HGAP,1,X)
          ENDIF
       ENDIF

    ENDIF
  END SUBROUTINE fringeP_STRAIGHT

  SUBROUTINE FRINGER_CAV4(EL,X,J)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) ::  X(6)
    TYPE(CAV4),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i,JC
    REAL(DP) C1,S1,V,O

    JC=-2*J+3

    IF(EL%P%NOCAVITY) RETURN
    IF(.NOT.EL%P%FRINGE) RETURN
    IF(EL%THIN) RETURN
    IF(jC==1.AND.EL%P%KILL_ENT_FRINGE) RETURN
    IF(jC==-1.AND.EL%P%KILL_EXI_FRINGE) RETURN


    O=EL%freq*twopi/CLIGHT
    V=jC*EL%P%CHARGE*EL%volt*c_1d_3/EL%P%P0C



    C1=COS(O*(x(6))+EL%PHAS+phase0)
    S1=SIN(O*(x(6))+EL%PHAS+phase0)



    X(2)=X(2)+V*S1*X(1)
    X(4)=X(4)+V*S1*X(3)
    x(5)=x(5)-HALF*(X(1)**2+X(3)**2)*V*C1*O




  END SUBROUTINE FRINGER_CAV4

  SUBROUTINE FRINGEP_CAV4(EL,X,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) ::  X(6)
    TYPE(CAV4P),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i,JC
    TYPE(REAL_8) C1,S1,V,O

    CALL ALLOC(C1,S1,V,O)

    JC=-2*J+3

    IF(EL%P%NOCAVITY) RETURN
    IF(.NOT.EL%P%FRINGE) RETURN
    IF(EL%THIN) RETURN
    IF(jC==1.AND.EL%P%KILL_ENT_FRINGE) RETURN
    IF(jC==-1.AND.EL%P%KILL_EXI_FRINGE) RETURN

    O=EL%freq*twopi/CLIGHT
    V=jC*EL%P%CHARGE*EL%volt*c_1d_3/EL%P%P0C

    C1=COS(O*(x(6))+EL%PHAS+phase0)
    S1=SIN(O*(x(6))+EL%PHAS+phase0)

    X(2)=X(2)+V*S1*X(1)
    X(4)=X(4)+V*S1*X(3)
    x(5)=x(5)-HALF*(X(1)**2+X(3)**2)*V*C1*O

    CALL KILL(C1,S1,V,O)

  END SUBROUTINE FRINGEP_CAV4

  SUBROUTINE ADJUSTR_TIME_CAV4(EL,X,J)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT)::  X(6)
    TYPE(CAV4),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i

    IF(J==1) THEN
       EL%DELTA_E=X(5)
       IF(EL%P%NOCAVITY) RETURN

       IF(EL%THIN) THEN
          CALL CAVITY(EL,X)
          EL%DELTA_E=(X(5)-EL%DELTA_E)*EL%P%P0C
          RETURN
       ENDIF

    ELSE
       IF(EL%THIN) RETURN

       if(EL%P%TIME) then
          X(6)=X(6)-(CAVITY_TOTALPATH-EL%P%TOTALPATH)*EL%P%LD/EL%P%BETA0
       else
          X(6)=X(6)-(CAVITY_TOTALPATH-EL%P%TOTALPATH)*EL%P%LD
       endif

       EL%DELTA_E=(X(5)-EL%DELTA_E)*EL%P%P0C
    ENDIF

  END SUBROUTINE ADJUSTR_TIME_CAV4

  SUBROUTINE ADJUSTP_TIME_CAV4(EL,X,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT)::  X(6)
    TYPE(CAV4P),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i

    IF(J==1) THEN
       EL%DELTA_E=X(5)
       IF(EL%P%NOCAVITY) RETURN

       IF(EL%THIN) THEN
          CALL CAVITY(EL,X)
          EL%DELTA_E=(X(5)-EL%DELTA_E)*EL%P%P0C
          RETURN
       ENDIF

    ELSE
       IF(EL%THIN) RETURN

       if(EL%P%TIME) then
          X(6)=X(6)-(CAVITY_TOTALPATH-EL%P%TOTALPATH)*EL%P%LD/EL%P%BETA0
       else
          X(6)=X(6)-(CAVITY_TOTALPATH-EL%P%TOTALPATH)*EL%P%LD
       endif

       EL%DELTA_E=(X(5)-EL%DELTA_E)*EL%P%P0C
    ENDIF

  END SUBROUTINE ADJUSTP_TIME_CAV4

  SUBROUTINE fringer_TEAPOT(EL,X,J)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(TEAPOT),INTENT(IN):: EL
    integer,INTENT(IN):: J
    INTEGER I

    IF(EL%P%DIR==1) THEN

       IF(J==1) THEN

          IF(EL%P%EDGE(1)/=zero) THEN
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) then
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
                x(2)=x(2)+EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)-EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)+EL%P%EDGE(1)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)-EL%P%EDGE(1)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL WEDGE(-EL%P%EDGE(1),X,EL2=EL)
          ELSE
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ENDIF

       ELSE ! J=2

          IF(EL%P%EDGE(2)/=zero) THEN
             CALL WEDGE(-EL%P%EDGE(2),X,EL2=EL)
             IF(EL%P%FRINGE) then
                x(2)=x(2)+EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)-EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)+EL%P%EDGE(2)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)-EL%P%EDGE(2)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
          ENDIF
       ENDIF ! J=2
    ELSE


       IF(J==1) THEN

          IF(EL%P%EDGE(2)/=zero) THEN
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) then
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
                x(2)=x(2)-EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)+EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)-EL%P%EDGE(2)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)+EL%P%EDGE(2)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL WEDGE(-EL%P%EDGE(2),X,EL2=EL)
          ELSE
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF

       ELSE ! J=2

          IF(EL%P%EDGE(1)/=zero) THEN
             CALL WEDGE(-EL%P%EDGE(1),X,EL2=EL)
             IF(EL%P%FRINGE) then
                x(2)=x(2)-EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)+EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)-EL%P%EDGE(1)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)+EL%P%EDGE(1)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
          ENDIF

       ENDIF

    ENDIF

  END SUBROUTINE fringer_TEAPOT

  SUBROUTINE fringeP_TEAPOT(EL,X,J)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(TEAPOTP),INTENT(IN):: EL
    integer,INTENT(IN):: J
    INTEGER I

    IF(EL%P%DIR==1) THEN

       IF(J==1) THEN

          IF(EL%P%EDGE(1)/=zero) THEN
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) then
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
                x(2)=x(2)+EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)-EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)+EL%P%EDGE(1)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)-EL%P%EDGE(1)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL WEDGE(-EL%P%EDGE(1),X,EL2=EL)
          ELSE
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ENDIF

       ELSE ! J=2

          IF(EL%P%EDGE(2)/=zero) THEN
             CALL WEDGE(-EL%P%EDGE(2),X,EL2=EL)
             IF(EL%P%FRINGE) then
                x(2)=x(2)+EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)-EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)+EL%P%EDGE(2)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)-EL%P%EDGE(2)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
          ENDIF
       ENDIF ! J=2
    ELSE


       IF(J==1) THEN

          IF(EL%P%EDGE(2)/=zero) THEN
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) then
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
                x(2)=x(2)-EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)+EL%P%EDGE(2)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)-EL%P%EDGE(2)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)+EL%P%EDGE(2)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL WEDGE(-EL%P%EDGE(2),X,EL2=EL)
          ELSE
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF

       ELSE ! J=2

          IF(EL%P%EDGE(1)/=zero) THEN
             CALL WEDGE(-EL%P%EDGE(1),X,EL2=EL)
             IF(EL%P%FRINGE) then
                x(2)=x(2)-EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(1)*x(1)**2-wedge_coeff(2)*x(3)**2*half)
                x(4)=x(4)+EL%P%EDGE(1)*el%bn(2)*(wedge_coeff(2)*x(1)*x(3))
                CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             ELSEIF(MAD8_WEDGE) THEN
                x(2)=x(2)-EL%P%EDGE(1)*el%bn(2)*(x(1)**2-x(3)**2)
                x(4)=x(4)+EL%P%EDGE(1)*el%bn(2)*(TWO*x(1)*x(3))
             endif
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
          ENDIF

       ENDIF

    ENDIF

  END SUBROUTINE fringeP_TEAPOT

  SUBROUTINE  fringer_STREX(EL,X,J)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    integer,INTENT(IN):: J
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(STREX),INTENT(IN):: EL
    real(dp) ANGH
    INTEGER I
    ! J=1 front

    IF(EL%P%DIR==1) THEN

       IF(J==1) THEN

          IF(EL%LIKEMAD) THEN

             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(1)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL  WEDGE(ANGH,X,EL1=EL)

          ELSE
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ENDIF

       ELSE  ! J==2


          IF(EL%LIKEMAD) THEN
             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(2)
             CALL  WEDGE(ANGH,X,EL1=EL)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
          ENDIF

       ENDIF ! J

    ELSE

       IF(J==1) THEN

          IF(EL%LIKEMAD) THEN

             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(2)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL  WEDGE(ANGH,X,EL1=EL)
          ELSE
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF


       ELSE ! J==2

          IF(EL%LIKEMAD) THEN
             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(1)
             CALL  WEDGE(ANGH,X,EL1=EL)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
          ENDIF

       ENDIF ! J

    ENDIF

  END SUBROUTINE fringer_STREX

  SUBROUTINE  fringeP_STREX(EL,X,J)
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    integer,INTENT(IN):: J
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(STREXP),INTENT(IN):: EL
    real(dp) ANGH
    INTEGER I
    ! J=1 front

    IF(EL%P%DIR==1) THEN

       IF(J==1) THEN

          IF(EL%LIKEMAD) THEN

             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(1)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL  WEDGE(ANGH,X,EL1=EL)

          ELSE
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
          ENDIF

       ELSE  ! J==2


          IF(EL%LIKEMAD) THEN
             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(2)
             CALL  WEDGE(ANGH,X,EL1=EL)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
          ENDIF

       ENDIF ! J

    ELSE

       IF(J==1) THEN

          IF(EL%LIKEMAD) THEN

             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(2)
             CALL ROT_XZ(EL%P%EDGE(2),X,EL%P%BETA0,DONEITT,EL%P%TIME)
             CALL FACE(EL%P,EL%BN,EL%H2,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
             CALL  WEDGE(ANGH,X,EL1=EL)
          ELSE
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,2,X)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,2,X)
          ENDIF


       ELSE ! J==2

          IF(EL%LIKEMAD) THEN
             ANGH=EL%P%B0*EL%P%LD*half-EL%P%EDGE(1)
             CALL  WEDGE(ANGH,X,EL1=EL)
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL FRINGE_(EL%P,EL%BN,EL%FINT,EL%HGAP,1,X)
             CALL FACE(EL%P,EL%BN,EL%H1,X)
             CALL ROT_XZ(EL%P%EDGE(1),X,EL%P%BETA0,DONEITT,EL%P%TIME)
          ELSE
             IF(EL%P%FRINGE) CALL MULTIPOLE_FRINGE(EL%P,EL%AN,EL%BN,1,X)
             CALL EDGE_TRUE_PARALLEL(EL%P,EL%BN,EL%H1,EL%H2,EL%FINT,EL%HGAP,1,X)
          ENDIF

       ENDIF ! J

    ENDIF

  END SUBROUTINE fringeP_STREX

  SUBROUTINE  FRINGER_CAV_TRAV(EL,X,J)
    IMPLICIT NONE
    !      TYPE(BEAM), INTENT(INOUT) ::B
    integer,INTENT(IN):: J
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAV),INTENT(INOUT):: EL
    INTEGER I
    ! J=1 front
    IF(J==1) THEN

       CALL FRINGECAV_TRAV(EL,EL%P%DIR,X)

    ELSE

       CALL FRINGECAV_TRAV(EL,-EL%P%DIR,X)

    ENDIF

  END SUBROUTINE FRINGER_CAV_TRAV

  SUBROUTINE  FRINGEP_CAV_TRAV(EL,X,J)
    IMPLICIT NONE
    !      TYPE(BEAM), INTENT(INOUT) ::B
    integer,INTENT(IN):: J
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAVP),INTENT(INOUT):: EL
    INTEGER I
    ! J=1 front
    IF(J==1) THEN

       CALL FRINGECAV_TRAV(EL,EL%P%DIR,X)

    ELSE

       CALL FRINGECAV_TRAV(EL,-EL%P%DIR,X)

    ENDIF

  END SUBROUTINE FRINGEP_CAV_TRAV

  SUBROUTINE ADJUSTR_TIME_CAV_TRAV_OUT(EL,X,J)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAV),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i

    IF(J==1) RETURN

    if(EL%P%TIME) then
       X(6)=X(6)-(1-EL%P%TOTALPATH)*EL%P%LD/EL%P%BETA0
    else
       X(6)=X(6)-(1-EL%P%TOTALPATH)*EL%P%LD
    endif

  END SUBROUTINE ADJUSTR_TIME_CAV_TRAV_OUT

  SUBROUTINE ADJUSTP_TIME_CAV_TRAV_OUT(EL,X,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAVP),INTENT(INOUT):: EL
    integer,INTENT(IN):: J
    integer i

    IF(J==1) RETURN

    if(EL%P%TIME) then
       X(6)=X(6)-(1-EL%P%TOTALPATH)*EL%P%LD/EL%P%BETA0
    else
       X(6)=X(6)-(1-EL%P%TOTALPATH)*EL%P%LD
    endif

  END SUBROUTINE ADJUSTP_TIME_CAV_TRAV_OUT

  ! tracking one steps in the body


  SUBROUTINE INTER_DRIFT1(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(DRIFT1),INTENT(IN):: EL

    real(dp) DH,DD

    SELECT CASE(EL%P%METHOD)
    CASE(2,4,6)
       DH=EL%L/EL%P%NST
       DD=EL%P%LD/EL%P%NST




       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_DRIFT1

  SUBROUTINE INTEP_DRIFT1(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(DRIFT1P),INTENT(IN):: EL
    TYPE(REAL_8) DH
    real(dp) DD

    CALL ALLOC(DH)
    SELECT CASE(EL%P%METHOD)
    CASE(2,4,6)
       DH=EL%L/EL%P%NST
       DD=EL%P%LD/EL%P%NST




       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    CALL KILL(DH)

  END SUBROUTINE INTEP_DRIFT1

  SUBROUTINE INTER_dkd2 (EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(DKD2),INTENT(IN):: EL

    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST




       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,D,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(4)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST




       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK2,X)

       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(6)
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)

       DO J=4,2,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICK (EL,DK(J),X)
       ENDDO
       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK(1),X)

       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       DO J=2,4
          CALL KICK (EL,DK(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO




    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_dkd2

  SUBROUTINE INTEP_dkd2 (EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(DKD2P),INTENT(IN):: EL

    real(dp) DD
    real(dp) DD1,DD2
    real(dp) DDF(4)
    TYPE(REAL_8) DH,D,D1,D2,DK1,DK2,DF(4),DK(4)

    INTEGER I,J

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       CALL ALLOC(DH,D)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST


       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,D,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CALL KILL(DH,D)

    CASE(4)
       CALL ALLOC(D1,D2,DK1,DK2)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST



       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK2,X)

       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CALL KILL(D1,D2,DK1,DK2)


    CASE(6)
       CALL ALLOC(DF,4);CALL ALLOC(DK,4);
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)

       DO J=4,2,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICK (EL,DK(J),X)
       ENDDO
       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICK (EL,DK(1),X)

       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       DO J=2,4
          CALL KICK (EL,DK(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO


       CALL KILL(DF,4);CALL KILL(DK,4);


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTEP_dkd2

  SUBROUTINE INTER_KICKT3(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(KICKT3),INTENT(IN):: EL

    CALL TRACK(EL,X)

  END SUBROUTINE INTER_KICKT3

  SUBROUTINE INTEP_KICKT3(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(KICKT3P),INTENT(IN):: EL

    CALL TRACK(EL,X)

  END SUBROUTINE INTEP_KICKT3

  SUBROUTINE INTER_NSMI(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(NSMI),INTENT(IN):: EL

    CALL TRACK(EL,X)

  END SUBROUTINE INTER_NSMI

  SUBROUTINE INTEP_NSMI(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(NSMIP),INTENT(IN):: EL

    CALL TRACK(EL,X)

  END SUBROUTINE INTEP_NSMI

  SUBROUTINE INTER_SSMI(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(SSMI),INTENT(IN):: EL
    INTEGER I

    CALL TRACK(EL,X)

  END SUBROUTINE INTER_SSMI

  SUBROUTINE INTEP_SSMI(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(SSMIP),INTENT(IN):: EL
    INTEGER I

    CALL TRACK(EL,X)

  END SUBROUTINE INTEP_SSMI

  SUBROUTINE INTER_CAV4(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  X(6)

    TYPE(CAV4),INTENT(INOUT):: EL
    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J


    TOTALPATH_FLAG=EL%P%TOTALPATH
    EL%P%TOTALPATH=CAVITY_TOTALPATH



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,D,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)


    CASE(4)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK2,X)

       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(6)
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       DO J=4,2,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKCAV (EL,DK(J),X)
       ENDDO
       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK(1),X)

       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       DO J=2,4
          CALL KICKCAV(EL,DK(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO




    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    EL%P%TOTALPATH=TOTALPATH_FLAG


  END SUBROUTINE INTER_CAV4

  SUBROUTINE INTEP_CAV4(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) ::  X(6)

    TYPE(CAV4P),INTENT(INOUT):: EL
    real(dp) DD
    real(dp) DD1,DD2
    real(dp) DDF(4)
    TYPE(REAL_8) DH,D,D1,D2,DK1,DK2,DF(4),DK(4)
    INTEGER I,J


    TOTALPATH_FLAG=EL%P%TOTALPATH
    EL%P%TOTALPATH=CAVITY_TOTALPATH



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       CALL ALLOC(DH,D)

       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,D,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CALL KILL(DH,D)

    CASE(4)
       CALL ALLOC(D1,D2,DK1,DK2)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK2,X)

       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CALL KILL(D1,D2,DK1,DK2)


    CASE(6)
       CALL ALLOC(DF,4);CALL ALLOC(DK,4);
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       DO J=4,2,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKCAV (EL,DK(J),X)
       ENDDO
       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL KICKCAV (EL,DK(1),X)

       CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       DO J=2,4
          CALL KICKCAV(EL,DK(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO


       CALL KILL(DF,4);CALL KILL(DK,4);


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    EL%P%TOTALPATH=TOTALPATH_FLAG


  END SUBROUTINE INTEP_CAV4

  SUBROUTINE INTER_SOL5(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::  X(6)
    TYPE(SOL5),INTENT(IN):: EL

    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2,D2H
    real(dp) dd1,dd2,DK(4),DF(4),DDF(4)
    INTEGER I,J




    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=(EL%P%LD)/two/EL%P%NST


       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,DH,X)
       CALL KICK_SOL(EL,DH,X)
       CALL KICKMUL(EL,D,X)

       CALL KICK_SOL(EL,DH,X)
       CALL SOL_ROT (EL,DH,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(4)

       D=EL%L/EL%P%NST
       D1=D*FD1
       D2=D*FD2
       DK1=D*FK1
       DK2=D*FK2
       D2H=DK2/two
       DD1=(EL%P%LD)/EL%P%NST*FD1
       DD2=(EL%P%LD)/EL%P%NST*FD2

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL KICKMUL(EL,DK1,X)

       CALL KICK_SOL(EL,D1,X)
       CALL SOL_ROT (EL,D1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D2H,X)
       CALL KICK_SOL(EL,D2H,X)
       CALL KICKMUL(EL,DK2,X)
       CALL KICK_SOL(EL,D2H,X)
       CALL SOL_ROT (EL,D2H,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL KICKMUL(EL,DK1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL SOL_ROT (EL,D1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)



    CASE(6)
       DO I =1,4
          DK(I)=EL%L*YOSK(I)/EL%P%NST
          DF(I)=DK(I)/two
          DDF(I)=EL%P%LD*YOSK(I)/two/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       DO J=4,1,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL KICKMUL(EL,DK(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO

       DO J=2,4
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL KICKMUL(EL,DK(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO


    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT



  END SUBROUTINE INTER_SOL5

  SUBROUTINE INTEP_SOL5(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) ::  X(6)
    TYPE(SOL5P),INTENT(IN):: EL

    real(dp) DD
    real(dp) DD1,DD2
    real(dp) DDF(4)
    TYPE(REAL_8) DH,D,D1,D2,DK1,DK2,DF(4),DK(4),D2H
    INTEGER I,J




    SELECT CASE(EL%P%METHOD)
    CASE(2)
       CALL ALLOC(DH,D)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=(EL%P%LD)/two/EL%P%NST


       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,DH,X)
       CALL KICK_SOL(EL,DH,X)
       CALL KICKMUL(EL,D,X)

       CALL KICK_SOL(EL,DH,X)
       CALL SOL_ROT (EL,DH,X)
       CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CALL KILL(DH,D)


    CASE(4)
       CALL ALLOC(D1,D2,DK1,DK2,D2H)

       D=EL%L/EL%P%NST
       D1=D*FD1
       D2=D*FD2
       DK1=D*FK1
       DK2=D*FK2
       D2H=DK2/two
       DD1=(EL%P%LD)/EL%P%NST*FD1
       DD2=(EL%P%LD)/EL%P%NST*FD2

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL KICKMUL(EL,DK1,X)

       CALL KICK_SOL(EL,D1,X)
       CALL SOL_ROT (EL,D1,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D2H,X)
       CALL KICK_SOL(EL,D2H,X)
       CALL KICKMUL(EL,DK2,X)
       CALL KICK_SOL(EL,D2H,X)
       CALL SOL_ROT (EL,D2H,X)
       CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       CALL SOL_ROT (EL,D1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL KICKMUL(EL,DK1,X)
       CALL KICK_SOL(EL,D1,X)
       CALL SOL_ROT (EL,D1,X)
       CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CALL KILL(D1,D2,DK1,DK2,D2H)


    CASE(6)
       CALL ALLOC(DF,4);CALL ALLOC(DK,4);
       DO I =1,4
          DK(I)=EL%L*YOSK(I)/EL%P%NST
          DF(I)=DK(I)/two
          DDF(I)=EL%P%LD*YOSK(I)/two/EL%P%NST
       ENDDO

       !       DO I=1,B%N

       !        X=BEAM_IN_X(B,I)
       DO J=4,1,-1
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL KICKMUL(EL,DK(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO

       DO J=2,4
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL KICKMUL(EL,DK(J),X)
          CALL KICK_SOL(EL,DF(J),X)
          CALL SOL_ROT (EL,DF(J),X)
          CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
       ENDDO

       CALL KILL(DF,4);CALL KILL(DK,4);

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTEP_SOL5


  SUBROUTINE INTER_KTK(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(KTK),INTENT(INOUT):: EL

    INTEGER I
    real(dp) DK,DK2,DK6,DK4,DK5



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DK2=EL%L/EL%P%NST
       DK=DK2/two

       IF(OLD_IMPLEMENTATION_OF_SIXTRACK) THEN

          CALL GETMAT(EL,X)
          CALL KICKKTK(EL,DK,X)  ! NEW
          CALL KICKPATH(EL,DK,X)
          CALL PUSHKTK(EL,X)

          CALL KICKPATH(EL,DK,X)
          if(EL%P%RADIATION) CALL GETMAT(EL,X)
          CALL KICKKTK(EL,DK,X)  ! NEW
       ELSE  ! OLD_IMPLEMENTATION_OF_SIXTRACK
          CALL GETMAT(EL,X)
          CALL PUSHKTK(EL,X)
          CALL KICKPATH(EL,DK,X)
          CALL KICKKTK(EL,DK2,X)  ! NEW
          CALL KICKPATH(EL,DK,X)
          CALL PUSHKTK(EL,X)
          if(EL%P%RADIATION) CALL GETMAT(EL,X)

       ENDIF ! OLD_IMPLEMENTATION_OF_SIXTRACK

    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two
       CALL GETMAT(EL,X)
       CALL KICKKTK(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKKTK(EL,DK6,X)
       CALL KICKPATH(EL,DK2,X)
       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKKTK(EL,DK,X) ! NEW
       if(EL%P%RADIATION) CALL GETMAT(EL,X)
    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two


       CALL GETMAT(EL,X)
       CALL KICKKTK(EL,DK,X)   ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHKTK(EL,X)

       CALL KICKKTK(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHKTK(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKKTK(EL,DK6,X)   ! SYMMETRY POINT
       CALL KICKPATH(EL,DK5,X)

       CALL PUSHKTK(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKKTK(EL,DK4,X)

       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKKTK(EL,DK,X)   ! NEW
       if(EL%P%RADIATION) CALL GETMAT(EL,X)

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_KTK

  SUBROUTINE INTEP_KTK(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(KTKP),INTENT(INOUT):: EL

    INTEGER I
    TYPE(REAL_8) DK,DK2,DK6,DK4,DK5



    SELECT CASE(EL%P%METHOD)
    CASE(2)
       CALL ALLOC(DK2,DK)
       DK2=EL%L/EL%P%NST
       DK=DK2/two

       IF(OLD_IMPLEMENTATION_OF_SIXTRACK) THEN



          CALL GETMAT(EL,X)
          CALL KICKKTK(EL,DK,X)  ! NEW
          CALL KICKPATH(EL,DK,X)
          CALL PUSHKTK(EL,X)

          CALL KICKPATH(EL,DK,X)
          if(EL%P%RADIATION) CALL GETMAT(EL,X)
          CALL KICKKTK(EL,DK,X)  ! NEW
       ELSE  ! OLD_IMPLEMENTATION_OF_SIXTRACK
          CALL GETMAT(EL,X)
          CALL PUSHKTK(EL,X)
          CALL KICKPATH(EL,DK,X)
          CALL KICKKTK(EL,DK2,X)  ! NEW
          CALL KICKPATH(EL,DK,X)
          CALL PUSHKTK(EL,X)
          if(EL%P%RADIATION) CALL GETMAT(EL,X)

       ENDIF ! OLD_IMPLEMENTATION_OF_SIXTRACK
       CALL KILL(DK2,DK)

    CASE(4)
       CALL ALLOC(DK,DK2,DK6)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two
       CALL GETMAT(EL,X)
       CALL KICKKTK(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKKTK(EL,DK6,X)
       CALL KICKPATH(EL,DK2,X)
       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKKTK(EL,DK,X) ! NEW
       if(EL%P%RADIATION) CALL GETMAT(EL,X)
       CALL KILL(DK,DK2,DK6)
    CASE(6)
       CALL ALLOC(DK,DK2,DK4,DK5,DK6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two


       CALL GETMAT(EL,X)
       CALL KICKKTK(EL,DK,X)   ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHKTK(EL,X)

       CALL KICKKTK(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHKTK(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKKTK(EL,DK6,X)   ! SYMMETRY POINT
       CALL KICKPATH(EL,DK5,X)

       CALL PUSHKTK(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKKTK(EL,DK4,X)

       CALL PUSHKTK(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKKTK(EL,DK,X)   ! NEW
       if(EL%P%RADIATION) CALL GETMAT(EL,X)
       CALL KILL(DK,DK2,DK4,DK5,DK6)

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTEP_KTK

  SUBROUTINE INTER_TKTF(EL,X,pos_in_fibre)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(TKTF),INTENT(INOUT):: EL
    INTEGER I,pos_in_fibre
    real(dp) DK,DK2,DK6,DK4,DK5



    SELECT CASE(EL%P%METHOD)
    CASE(1)
       DK=EL%L/EL%P%NST
       DK2=two*dk


       if(mod(pos_in_fibre,2)==1) then
          CALL PUSHTKT7(EL,X)
          CALL KICKPATH(EL,DK,X)
          CALL KICKTKT7(EL,DK2,X)
          CALL KICKPATH(EL,DK,X)
       else
          CALL PUSHTKT7(EL,X)
       endif

    CASE(2)
       DK2=EL%L/EL%P%NST
       DK=DK2/two

       !       CALL KICKTKT7(EL,DK2,X)
       !       CALL KICKPATH(EL,DK2,X)
       !       CALL PUSHTKT7(EL,X)
       !       CALL KICKPATH(EL,DK2,X)
       !       CALL KICKTKT7(EL,DK2,X)

       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK2,X)
       CALL KICKPATH(EL,DK,X)
       CALL PUSHTKT7(EL,X)

    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two

       CALL KICKTKT7(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKTKT7(EL,DK6,X)

       CALL KICKPATH(EL,DK2,X)
       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK,X) ! NEW

    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two

       CALL KICKTKT7(EL,DK,X)  ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKTKT7(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKTKT7(EL,DK6,X)   ! SYMMETRY POINT

       CALL KICKPATH(EL,DK5,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKTKT7(EL,DK4,X)

       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK,X)  ! NEW

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_TKTF

  SUBROUTINE INTEP_TKTF(EL,X,pos_in_fibre)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(TKTFP),INTENT(INOUT):: EL
    INTEGER I,pos_in_fibre
    TYPE(REAL_8) DK,DK2,DK6,DK4,DK5



    SELECT CASE(EL%P%METHOD)
    CASE(1)
       CALL ALLOC(DK,DK2)
       DK=EL%L/EL%P%NST
       DK2=two*dk


       if(mod(pos_in_fibre,2)==1) then
          CALL PUSHTKT7(EL,X)
          CALL KICKPATH(EL,DK,X)
          CALL KICKTKT7(EL,DK2,X)
          CALL KICKPATH(EL,DK,X)
       else
          CALL PUSHTKT7(EL,X)
       endif

       CALL KILL(DK,DK2)
    CASE(2)
       CALL ALLOC(DK,DK2)

       DK2=EL%L/EL%P%NST
       DK=DK2/two

       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK2,X)
       CALL KICKPATH(EL,DK,X)
       CALL PUSHTKT7(EL,X)


       CALL KILL(DK,DK2)

    CASE(4)
       CALL ALLOC(DK,DK2,DK6)

       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two

       CALL KICKTKT7(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKTKT7(EL,DK6,X)

       CALL KICKPATH(EL,DK2,X)
       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK,X) ! NEW

       CALL KILL(DK,DK2,DK6)

    CASE(6)
       CALL ALLOC(DK,DK2,DK6,DK4,DK5)

       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two

       CALL KICKTKT7(EL,DK,X)  ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKTKT7(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKTKT7(EL,DK6,X)   ! SYMMETRY POINT

       CALL KICKPATH(EL,DK5,X)

       CALL PUSHTKT7(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKTKT7(EL,DK4,X)

       CALL PUSHTKT7(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKTKT7(EL,DK,X)  ! NEW

       CALL KILL(DK,DK2,DK6,DK4,DK5)

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTEP_TKTF

  SUBROUTINE INTER_TEAPOT(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(TEAPOT),INTENT(IN):: EL
    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J


    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST

       CALL SSECH1(EL,DH,DD,X)
       CALL SKICK(EL,D,X)
       CALL SSECH1(EL,DH,DD,X)
    CASE(4)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST

       CALL SSECH1(EL,D1,DD1,X)
       CALL SKICK (EL,DK1,X)
       CALL SSECH1(EL,D2,DD2,X)
       CALL SKICK (EL,DK2,X)

       CALL SSECH1(EL,D2,DD2,X)
       CALL SKICK (EL,DK1,X)
       CALL SSECH1(EL,D1,DD1,X)



    CASE(6)
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       DO J=4,2,-1
          CALL SSECH1(EL,DF(J),DDF(J),X)
          CALL SKICK (EL,DK(J),X)
       ENDDO
       CALL SSECH1(EL,DF(1),DDF(1),X)
       CALL SKICK (EL,DK(1),X)

       CALL SSECH1(EL,DF(1),DDF(1),X)
       DO J=2,4
          CALL SKICK (EL,DK(J),X)
          CALL SSECH1(EL,DF(J),DDF(J),X)
       ENDDO



    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_TEAPOT

  SUBROUTINE INTEP_TEAPOT(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(TEAPOTP),INTENT(IN):: EL
    real(dp) DD
    real(dp) DD1,DD2
    real(dp) DDF(4)
    TYPE(REAL_8) DH,D,D1,D2,DK1,DK2,DF(4),DK(4)
    INTEGER I,J


    SELECT CASE(EL%P%METHOD)
    CASE(2)

       CALL ALLOC(DH,D)

       DH=EL%L/two/EL%P%NST
       D=EL%L/EL%P%NST
       DD=EL%P%LD/two/EL%P%NST

       CALL SSECH1(EL,DH,DD,X)
       CALL SKICK(EL,D,X)
       CALL SSECH1(EL,DH,DD,X)

       CALL KILL(DH,D)

    CASE(4)
       CALL ALLOC(D1,D2,DK1,DK2)
       D1=EL%L*FD1/EL%P%NST
       D2=EL%L*FD2/EL%P%NST
       DD1=EL%P%LD*FD1/EL%P%NST
       DD2=EL%P%LD*FD2/EL%P%NST
       DK1=EL%L*FK1/EL%P%NST
       DK2=EL%L*FK2/EL%P%NST

       CALL SSECH1(EL,D1,DD1,X)
       CALL SKICK (EL,DK1,X)
       CALL SSECH1(EL,D2,DD2,X)
       CALL SKICK (EL,DK2,X)

       CALL SSECH1(EL,D2,DD2,X)
       CALL SKICK (EL,DK1,X)
       CALL SSECH1(EL,D1,DD1,X)

       CALL KILL(D1,D2,DK1,DK2)

    CASE(6)
       CALL ALLOC(DF,4);CALL ALLOC(DK,4);
       DO I =1,4
          DF(I)=EL%L*YOSD(I)/EL%P%NST
          DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
          DK(I)=EL%L*YOSK(I)/EL%P%NST
       ENDDO

       DO J=4,2,-1
          CALL SSECH1(EL,DF(J),DDF(J),X)
          CALL SKICK (EL,DK(J),X)
       ENDDO
       CALL SSECH1(EL,DF(1),DDF(1),X)
       CALL SKICK (EL,DK(1),X)

       CALL SSECH1(EL,DF(1),DDF(1),X)
       DO J=2,4
          CALL SKICK (EL,DK(J),X)
          CALL SSECH1(EL,DF(J),DDF(J),X)
       ENDDO

       CALL KILL(DF,4);CALL KILL(DK,4);

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTEP_TEAPOT

  SUBROUTINE INTER_STREX(EL,X)
    IMPLICIT NONE
    TYPE(STREX),INTENT(IN):: EL
    real(dp), INTENT(INOUT) ::  X(6)
    real(dp) D,DH,DD
    real(dp) D1,D2,DK1,DK2
    real(dp) DD1,DD2
    real(dp) DF(4),DK(4),DDF(4)
    INTEGER I,J

    IF(EL%DRIFTKICK) THEN

       SELECT CASE(EL%P%METHOD)
       CASE(2)
          DH=EL%L/two/EL%P%NST
          D=EL%L/EL%P%NST
          DD=EL%P%LD/two/EL%P%NST

          CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,D,X)
          CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CASE(4)
          D1=EL%L*FD1/EL%P%NST
          D2=EL%L*FD2/EL%P%NST
          DD1=EL%P%LD*FD1/EL%P%NST
          DD2=EL%P%LD*FD2/EL%P%NST
          DK1=EL%L*FK1/EL%P%NST
          DK2=EL%L*FK2/EL%P%NST

          CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK1,X)
          CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK2,X)
          CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK1,X)
          CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)

       CASE(6)
          DO I =1,4
             DF(I)=EL%L*YOSD(I)/EL%P%NST
             DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
             DK(I)=EL%L*YOSK(I)/EL%P%NST
          ENDDO

          DO J=4,2,-1
             CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
             CALL KICKEX (EL,DK(J),X)
          ENDDO
          CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK(1),X)
          CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          DO J=2,4
             CALL KICKEX (EL,DK(J),X)
             CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          ENDDO

       CASE DEFAULT
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A72))'
          WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
          call write_e(357)
       END SELECT
    ELSE
       SELECT CASE(EL%P%METHOD)
       CASE(2)
          DH=EL%L/two/EL%P%NST
          D=EL%L/EL%P%NST
          DD=EL%P%LD/two/EL%P%NST

          CALL SPAR(EL,DH,DD,X)
          CALL KICKEX (EL,D,X)
          CALL SPAR(EL,DH,DD,X)

       CASE(4)
          D1=EL%L*FD1/EL%P%NST
          D2=EL%L*FD2/EL%P%NST
          DD1=EL%P%LD*FD1/EL%P%NST
          DD2=EL%P%LD*FD2/EL%P%NST
          DK1=EL%L*FK1/EL%P%NST
          DK2=EL%L*FK2/EL%P%NST

          CALL SPAR(EL,D1,DD1,X)
          CALL KICKEX (EL,DK1,X)
          CALL SPAR(EL,D2,DD2,X)
          CALL KICKEX (EL,DK2,X)
          CALL SPAR(EL,D2,DD2,X)
          CALL KICKEX (EL,DK1,X)
          CALL SPAR(EL,D1,DD1,X)

       CASE(6)
          DO I =1,4
             DF(I)=EL%L*YOSD(I)/EL%P%NST
             DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
             DK(I)=EL%L*YOSK(I)/EL%P%NST
          ENDDO

          DO J=4,2,-1
             CALL SPAR(EL,DF(J),DDF(J),X)
             CALL KICKEX (EL,DK(J),X)
          ENDDO
          CALL SPAR(EL,DF(1),DDF(1),X)
          CALL KICKEX (EL,DK(1),X)
          CALL SPAR(EL,DF(1),DDF(1),X)
          DO J=2,4
             CALL KICKEX (EL,DK(J),X)
             CALL SPAR(EL,DF(J),DDF(J),X)
          ENDDO

       CASE DEFAULT
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A72))'
          WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
          call write_e(357)
       END SELECT

    ENDIF


  END SUBROUTINE INTER_STREX

  SUBROUTINE INTEP_STREX(EL,X)
    IMPLICIT NONE
    TYPE(STREXP),INTENT(IN):: EL
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    real(dp) DD
    real(dp) DD1,DD2
    real(dp) DDF(4)
    TYPE(REAL_8) DH,D,D1,D2,DK1,DK2,DF(4),DK(4)
    INTEGER I,J

    IF(EL%DRIFTKICK) THEN

       SELECT CASE(EL%P%METHOD)
       CASE(2)
          CALL ALLOC(DH,D)
          DH=EL%L/two/EL%P%NST
          D=EL%L/EL%P%NST
          DD=EL%P%LD/two/EL%P%NST

          CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,D,X)
          CALL DRIFT(DH,DD,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KILL(DH,D)

       CASE(4)
          CALL ALLOC(D1,D2,DK1,DK2)
          D1=EL%L*FD1/EL%P%NST
          D2=EL%L*FD2/EL%P%NST
          DD1=EL%P%LD*FD1/EL%P%NST
          DD2=EL%P%LD*FD2/EL%P%NST
          DK1=EL%L*FK1/EL%P%NST
          DK2=EL%L*FK2/EL%P%NST

          CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK1,X)
          CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK2,X)
          CALL DRIFT(D2,DD2,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK1,X)
          CALL DRIFT(D1,DD1,EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KILL(D1,D2,DK1,DK2)

       CASE(6)
          CALL ALLOC(DF,4);CALL ALLOC(DK,4);
          DO I =1,4
             DF(I)=EL%L*YOSD(I)/EL%P%NST
             DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
             DK(I)=EL%L*YOSK(I)/EL%P%NST
          ENDDO

          DO J=4,2,-1
             CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
             CALL KICKEX (EL,DK(J),X)
          ENDDO
          CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          CALL KICKEX (EL,DK(1),X)
          CALL DRIFT(DF(1),DDF(1),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          DO J=2,4
             CALL KICKEX (EL,DK(J),X)
             CALL DRIFT(DF(J),DDF(J),EL%P%beta0,EL%P%TOTALPATH,EL%P%EXACT,EL%P%TIME,X)
          ENDDO
          CALL KILL(DF,4);CALL KILL(DK,4);

       CASE DEFAULT
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A72))'
          WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
          call write_e(357)
       END SELECT
    ELSE
       SELECT CASE(EL%P%METHOD)
       CASE(2)
          CALL ALLOC(DH,D)
          DH=EL%L/two/EL%P%NST
          D=EL%L/EL%P%NST
          DD=EL%P%LD/two/EL%P%NST

          CALL SPAR(EL,DH,DD,X)
          CALL KICKEX (EL,D,X)
          CALL SPAR(EL,DH,DD,X)
          CALL KILL(DH,D)

       CASE(4)
          CALL ALLOC(D1,D2,DK1,DK2)
          D1=EL%L*FD1/EL%P%NST
          D2=EL%L*FD2/EL%P%NST
          DD1=EL%P%LD*FD1/EL%P%NST
          DD2=EL%P%LD*FD2/EL%P%NST
          DK1=EL%L*FK1/EL%P%NST
          DK2=EL%L*FK2/EL%P%NST

          CALL SPAR(EL,D1,DD1,X)
          CALL KICKEX (EL,DK1,X)
          CALL SPAR(EL,D2,DD2,X)
          CALL KICKEX (EL,DK2,X)
          CALL SPAR(EL,D2,DD2,X)
          CALL KICKEX (EL,DK1,X)
          CALL SPAR(EL,D1,DD1,X)
          CALL KILL(D1,D2,DK1,DK2)

       CASE(6)
          CALL ALLOC(DF,4);CALL ALLOC(DK,4);
          DO I =1,4
             DF(I)=EL%L*YOSD(I)/EL%P%NST
             DDF(I)=EL%P%LD*YOSD(I)/EL%P%NST
             DK(I)=EL%L*YOSK(I)/EL%P%NST
          ENDDO

          DO J=4,2,-1
             CALL SPAR(EL,DF(J),DDF(J),X)
             CALL KICKEX (EL,DK(J),X)
          ENDDO
          CALL SPAR(EL,DF(1),DDF(1),X)
          CALL KICKEX (EL,DK(1),X)
          CALL SPAR(EL,DF(1),DDF(1),X)
          DO J=2,4
             CALL KICKEX (EL,DK(J),X)
             CALL SPAR(EL,DF(J),DDF(J),X)
          ENDDO
          CALL KILL(DF,4);CALL KILL(DK,4);

       CASE DEFAULT
          w_p=0
          w_p%nc=1
          w_p%fc='(1(1X,A72))'
          WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
          call write_e(357)
       END SELECT

    ENDIF


  END SUBROUTINE INTEP_STREX

  SUBROUTINE INTER_SOLT(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(SOLT),INTENT(INOUT):: EL

    INTEGER I
    real(dp) DK,DK2,DK4,DK5,DK6,DKT



    DKT=EL%L/EL%P%NST

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DK2=EL%L/EL%P%NST
       DK=DK2/two

       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)  ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHSOL(EL,X)

       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DK2,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X)  ! NEW

    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two

       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKMUL(EL,DK6,X)
       CALL KICKPATH(EL,DK2,X)
       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X) ! NEW
       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DKT,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two

       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)   ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHSOL(EL,X)

       CALL KICKMUL(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHSOL(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKMUL(EL,DK6,X)   ! SYMMETRY POINT
       CALL KICKPATH(EL,DK5,X)

       CALL PUSHSOL(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKMUL(EL,DK4,X)

       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X)   ! NEW
       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DKT,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

  END SUBROUTINE INTER_SOLT

  SUBROUTINE INTEP_SOLT(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(SOLTP),INTENT(INOUT):: EL

    INTEGER I
    TYPE(REAL_8) DK,DK2,DK4,DK5,DK6,DKT

    CALL ALLOC(DK,DK2,DK4,DK5,DK6,DKT)

    DKT=EL%L/EL%P%NST

    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DK2=EL%L/EL%P%NST
       DK=DK2/two

       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)  ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHSOL(EL,X)

       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DK2,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X)  ! NEW

    CASE(4)
       DK2=EL%L/EL%P%NST/three
       DK6=two*DK2
       DK=DK2/two

       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)    ! NEW
       CALL KICKPATH(EL,DK,X)
       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK2,X)
       CALL KICKMUL(EL,DK6,X)
       CALL KICKPATH(EL,DK2,X)
       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X) ! NEW
       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DKT,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

    CASE(6)
       DK2=c_14*EL%L/EL%P%NST/c_90
       DK4=c_32*EL%L/EL%P%NST/c_90
       DK6=twelve*EL%L/EL%P%NST/c_90
       DK5=DK6/two
       DK=DK2/two

       CALL GETMATSOL(EL,X)
       CALL KICKMUL(EL,DK,X)   ! NEW
       CALL KICKPATH(EL,DK,X)

       CALL PUSHSOL(EL,X)

       CALL KICKMUL(EL,DK4,X)
       CALL KICKPATH(EL,DK4,X)

       CALL PUSHSOL(EL,X)

       CALL KICKPATH(EL,DK5,X)
       CALL KICKMUL(EL,DK6,X)   ! SYMMETRY POINT
       CALL KICKPATH(EL,DK5,X)

       CALL PUSHSOL(EL,X)

       CALL KICKPATH(EL,DK4,X)
       CALL KICKMUL(EL,DK4,X)

       CALL PUSHSOL(EL,X)
       CALL KICKPATH(EL,DK,X)
       CALL KICKMUL(EL,DK,X)   ! NEW
       if(EL%P%RADIATION) THEN
          CALL KICK_SOL(EL,DKT,X)  ! RADIATION
          CALL GETMATSOL(EL,X)
       ENDIF

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    CALL KILL(DK,DK2,DK4,DK5,DK6,DKT)

  END SUBROUTINE INTEP_SOLT

  SUBROUTINE INTER_MON_11_14(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(mon),INTENT(INOUT):: EL
    INTEGER I

    CALL MONTI(EL,X,I)

  END SUBROUTINE INTER_MON_11_14

  SUBROUTINE INTEP_MON_11_14(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(monP),INTENT(INOUT):: EL
    INTEGER I

    CALL MONTI(EL,X,I)

  END SUBROUTINE INTEP_MON_11_14

  SUBROUTINE INTER_ESEPTUM(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(ESEPTUM),INTENT(INOUT):: EL

    INTEGER I

    CALL SEPTTRACK(EL,X,0)

  END SUBROUTINE INTER_ESEPTUM

  SUBROUTINE INTEP_ESEPTUM(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(ESEPTUMP),INTENT(INOUT):: EL

    INTEGER I

    CALL SEPTTRACK(EL,X,0)

  END SUBROUTINE INTEP_ESEPTUM

  SUBROUTINE INTER_PANCAKE(EL,X,POS)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(PANCAKE),INTENT(INOUT):: EL
    INTEGER I,IS,POS
    real(dp) ti,h

    H=el%L/el%p%NST

    SELECT CASE(EL%P%METHOD)
    CASE(4)
       IF(EL%P%DIR==1) THEN
          IS=-5+2*POS    ! POS=3 BEGINNING
          call rk4_m(IS,h,el,X)
       else
          IS=el%p%NST-2*POS+6
          call rk4_m(IS,h,el,X)
       ENDIF

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT


  END SUBROUTINE INTER_PANCAKE

  SUBROUTINE INTEP_PANCAKE(EL,X,POS)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(PANCAKEP),INTENT(INOUT):: EL
    INTEGER I,IS,POS
    TYPE(REAL_8) ti,h

    CALL ALLOC(TI,H)

    H=el%L/el%p%NST

    SELECT CASE(EL%P%METHOD)
    CASE(4)
       IF(EL%P%DIR==1) THEN
          IS=-5+2*POS    ! POS=3 BEGINNING
          call rk4_m(IS,h,el,X)
       else
          IS=el%p%NST-2*POS+6
          call rk4_m(IS,h,el,X)
       ENDIF

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    CALL KILL(TI,H)

  END SUBROUTINE INTEP_PANCAKE

  SUBROUTINE INTER_CAV_TRAV(EL,X,Z)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAV),INTENT(INOUT):: EL

    real(dp) , INTENT(IN) :: Z
    real(dp) D1
    INTEGER I
    REAL(DP) Z0
    INTEGER TOTALPATH

    Z0=Z

    TOTALPATH=EL%P%TOTALPATH
    EL%P%TOTALPATH=1
    SELECT CASE(EL%P%METHOD)
    CASE(2)
       D1=EL%L/EL%P%NST

       call rk2_cav(z0,d1,el,x)

    CASE(4)

       D1=EL%L/EL%P%NST

       call rk4_cav(z0,d1,el,x)

    CASE(6)

       D1=EL%L/EL%P%NST

       call rk6_cav(z0,d1,el,x)

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    !    IF(EL%P%FRINGE)

    EL%P%TOTALPATH=TOTALPATH

  END SUBROUTINE INTER_CAV_TRAV

  SUBROUTINE INTEP_CAV_TRAV(EL,X,Z)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(CAV_TRAVP),INTENT(INOUT):: EL

    !    TYPE(REAL_8), INTENT(IN) :: Z
    REAL(DP), INTENT(IN) :: Z
    INTEGER I
    TYPE(REAL_8) Z0,D1
    INTEGER TOTALPATH

    CALL ALLOC(Z0,D1)

    Z0=Z

    TOTALPATH=EL%P%TOTALPATH
    EL%P%TOTALPATH=1
    SELECT CASE(EL%P%METHOD)
    CASE(2)
       D1=EL%L/EL%P%NST

       call rk2_cav(z0,d1,el,x)

    CASE(4)

       D1=EL%L/EL%P%NST

       call rk4_cav(z0,d1,el,x)

    CASE(6)

       D1=EL%L/EL%P%NST

       call rk6_cav(z0,d1,el,x)

    CASE DEFAULT
       w_p=0
       w_p%nc=1
       w_p%fc='(1(1X,A72))'
       WRITE(w_p%c(1),'(a12,1x,i4,1x,a17)') " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       call write_e(357)
    END SELECT

    !    IF(EL%P%FRINGE)

    EL%P%TOTALPATH=TOTALPATH

    CALL KILL(Z0,D1)


  END SUBROUTINE INTEP_CAV_TRAV

  SUBROUTINE INTER_RCOL18(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(RCOL),INTENT(INOUT):: EL
    INTEGER I

    CALL RCOLLIMATORI(EL,X,I)
  END SUBROUTINE INTER_RCOL18

  SUBROUTINE INTEP_RCOL18(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(RCOLP),INTENT(INOUT):: EL
    INTEGER I

    CALL RCOLLIMATORI(EL,X,I)
  END SUBROUTINE INTEP_RCOL18


  SUBROUTINE INTER_ECOL19(EL,X)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(ECOL),INTENT(INOUT):: EL

    INTEGER I

    CALL ECOLLIMATORI(EL,X,I)

  END SUBROUTINE INTER_ECOL19

  SUBROUTINE INTEP_ECOL19(EL,X)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT) :: X(6)
    TYPE(ECOLP),INTENT(INOUT):: EL

    INTEGER I

    CALL ECOLLIMATORI(EL,X,I)

  END SUBROUTINE INTEP_ECOL19

  SUBROUTINE INTER_WI(EL,X,Z)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) ::X(6)
    TYPE(sagan),INTENT(INOUT):: EL

    real(dp) , INTENT(IN) :: Z
    INTEGER I

    call INTR_SAGAN(EL,X,Z)

  END SUBROUTINE INTER_WI

  SUBROUTINE INTEP_WI(EL,X,Z)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) ::X(6)
    TYPE(saganP),INTENT(INOUT):: EL

    real(dp) , INTENT(IN) :: Z
    INTEGER I

    call INTP_SAGAN(EL,X,Z)

  END SUBROUTINE INTEP_WI

  SUBROUTINE ADJUSTR_PANCAKE(EL,X,J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(PANCAKE),INTENT(INOUT):: EL

    INTEGER, INTENT(IN) :: J
    INTEGER I

    IF(J==1) then
       call conv_to_xp(el,x)
    else
       call conv_to_px(el,x)
    endif

  END SUBROUTINE ADJUSTR_PANCAKE

  SUBROUTINE ADJUSTP_PANCAKE(EL,X,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(PANCAKEP),INTENT(INOUT):: EL

    INTEGER, INTENT(IN) :: J
    INTEGER I

    IF(J==1) then
       call conv_to_xp(el,x)
    else
       call conv_to_px(el,x)
    endif

  END SUBROUTINE ADJUSTP_PANCAKE


  SUBROUTINE ADJUSTR_WI(EL,X,J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(sagan),INTENT(INOUT):: EL

    INTEGER, INTENT(IN) :: J
    INTEGER I

    IF(J==1) RETURN

    X(1)=X(1)-EL%INTERNAL(1)
    X(2)=X(2)-EL%INTERNAL(2)

  END SUBROUTINE ADJUSTR_WI

  SUBROUTINE ADJUSTP_WI(EL,X,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(saganP),INTENT(INOUT):: EL

    INTEGER, INTENT(IN) :: J
    INTEGER I

    IF(J==1) RETURN

    X(1)=X(1)-EL%INTERNAL(1)
    X(2)=X(2)-EL%INTERNAL(2)

  END SUBROUTINE ADJUSTP_WI





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

    ! DIRECTIONAL VARIABLE
    !    C%MAG%P%DIR=>C%DIR
    !    if(present(charge)) then
    !       C%MAG%P%CHARGE=>CHARGE
    !    else
    !       charge1=1
    !       C%MAG%P%CHARGE=>CHARGE1
    !    endif
    !

    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

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
             IF(C%MAG%P%TIME)THEN
                X(5)=root(one+two*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAG%P%P0C-one !X(5) = DP/P0C_NEW
                X(5)=(two*X(5)+X(5)**2)/(root(one/C%MAG%P%BETA0**2+two*X(5)+X(5)**2)+one/C%MAG%P%BETA0)
             ELSE
                X(5)=(one+X(5))*P0/C%MAG%P%P0C-one
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG/=0.AND.PATCHG/=2) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%a_T
    ENDIF

    CALL DTILTD(C%DIR,C%MAG%P%TILTD,1,X)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAG%MIS) THEN
       ou = K%EXACTMIS.or.C%MAG%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
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

    ! DIRECTIONAL VARIABLE
    !    C%MAG%P%DIR=>C%DIR
    !    if(present(charge)) then
    !       C%MAG%P%CHARGE=>CHARGE
    !    else
    !       charge1=1
    !       C%MAG%P%CHARGE=>CHARGE1
    !    endif
    !


    !FRONTAL PATCH
    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF

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
             IF(C%MAGP%P%TIME)THEN
                X(5)=SQRT(one+two*X(5)/B0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
                X(5)=X(5)*P0/C%MAGP%P%P0C-one !X(5) = DP/P0C_NEW
                X(5)=(two*X(5)+X(5)**2)/(SQRT(one/C%MAGP%P%BETA0**2+two*X(5)+X(5)**2)+one/C%MAGP%P%BETA0)
             ELSE
                X(5)=(one+X(5))*P0/C%MAGP%P%P0C-one
             ENDIF
          ENDIF ! No need to patch
       ENDIF ! ASSOCIATED

    ENDIF

    ! The chart frame of reference is located here implicitely
    IF(PATCHG/=0.AND.PATCHG/=2) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_TRUE)
    ENDIF

    IF(PATCHT/=0.AND.PATCHT/=2.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%a_T
    ENDIF

    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,1,X)
    ! The magnet frame of reference is located here implicitely before misalignments

    !      CALL TRACK(C,X,EXACTMIS=K%EXACTMIS)
    IF(C%MAGP%MIS) THEN
       ou = K%EXACTMIS.or.C%MAGP%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITT)
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
    INTEGER I


    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAG%MIS) THEN
       ou = K%EXACTMIS.or.C%MAG%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%DIR,C%MAG%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%b_T
    ENDIF

    IF(PATCHG/=0.AND.PATCHG/=1) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAG%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely

    IF(PATCHE/=0.AND.PATCHE/=1) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAG%P%P0C
       B0=>CN%MAG%P%BETA0
       X(2)=X(2)*C%MAG%P%P0C/P0
       X(4)=X(4)*C%MAG%P%P0C/P0
       IF(C%MAG%P%TIME)THEN
          X(5)=root(one+two*X(5)/C%MAG%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAG%P%P0C/P0-one !X(5) = DP/P0C_NEW
          X(5)=(two*X(5)+X(5)**2)/(root(one/B0**2+two*X(5)+X(5)**2)+one/B0)
       ELSE
          X(5)=(one+X(5))*C%MAG%P%P0C/P0-one
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
    INTEGER I


    IF(ASSOCIATED(C%PATCH)) THEN
       PATCHT=C%PATCH%TIME ;PATCHE=C%PATCH%ENERGY ;PATCHG=C%PATCH%PATCH;
    ELSE
       PATCHT=0 ; PATCHE=0 ;PATCHG=0;
    ENDIF



    IF(C%MAGP%MIS) THEN
       ou = K%EXACTMIS.or.C%MAGP%EXACTMIS
       CALL MIS_FIB(C,X,OU,DONEITF)
    ENDIF
    ! The magnet frame of reference is located here implicitely before misalignments
    CALL DTILTD(C%DIR,C%MAGP%P%TILTD,2,X)

    IF(PATCHT/=0.AND.PATCHT/=1.AND.(.NOT.K%TOTALPATH)) THEN
       X(6)=X(6)+C%PATCH%b_T
    ENDIF

    IF(PATCHG/=0.AND.PATCHG/=1) THEN
       patch=ALWAYS_EXACT_PATCHING.or.C%MAGP%P%EXACT
       CALL PATCH_FIB(C,X,PATCH,MY_FALSE)
    ENDIF

    ! The CHART frame of reference is located here implicitely

    IF(PATCHE/=0.AND.PATCHE/=1) THEN
       NULLIFY(P0);NULLIFY(B0);
       CN=>C%NEXT
       IF(.NOT.ASSOCIATED(CN)) CN=>C
       P0=>CN%MAGP%P%P0C
       B0=>CN%MAGP%P%BETA0
       X(2)=X(2)*C%MAGP%P%P0C/P0
       X(4)=X(4)*C%MAGP%P%P0C/P0
       IF(C%MAGP%P%TIME)THEN
          X(5)=sqrt(one+two*X(5)/C%MAGP%P%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
          X(5)=X(5)*C%MAGP%P%P0C/P0-one !X(5) = DP/P0C_NEW
          X(5)=(two*X(5)+X(5)**2)/(sqrt(one/B0**2+two*X(5)+X(5)**2)+one/B0)
       ELSE
          X(5)=(one+X(5))*C%MAGP%P%P0C/P0-one
       ENDIF
    ENDIF

  END SUBROUTINE TRACK_FIBRE_BACKP


  ! thin lens tracking

  SUBROUTINE TRACKV_NODE_SINGLE(T,V,K,CHARGE)
    implicit none
    TYPE(INTEGRATION_NODE),POINTER :: T
    TYPE(INTERNAL_STATE)  K
    INTEGER, optional, TARGET :: CHARGE
    REAL(DP) SC,reference_ray(6),x(6)
    type(three_d_info)  v
    TYPE(INTEGRATION_NODE),POINTER:: mag_in,mag_out

    x=V%X
    reference_ray=V%reference_ray

    CALL TRACK_NODE_SINGLE(T,V%X,K,CHARGE)

    IF(.NOT.CHECK_STABLE)      V%U(1)=.TRUE.

    CALL TRACK_NODE_SINGLE(T,V%reference_ray,K,CHARGE)

    IF(.NOT.CHECK_STABLE)   V%U(2)=.TRUE.
    IF(V%U(1).OR.V%U(2)) RETURN

    IF(.NOT.ASSOCIATED(T%B)) THEN
       WRITE(6,*) " NO FRAMES IN INTEGRATION NODES "
       STOP 101
    ENDIF

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


  END SUBROUTINE TRACKV_NODE_SINGLE

  SUBROUTINE TRACKR_NODE_SINGLE(T,X,K,CHARGE)
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT):: T
    REAL(DP),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    INTEGER, optional, TARGET :: CHARGE
    integer, target :: CHARGE1
    INTEGER I

    ! call cpu_time(ttime0)
    if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
       CHECK_STABLE=.false.
    endif


    T%PARENT_FIBRE%MAG=K
    T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR
    if(present(charge))  then
       T%PARENT_FIBRE%MAG%P%CHARGE=>CHARGE
    else
       charge1=1
       T%PARENT_FIBRE%MAG%P%CHARGE=>CHARGE1
    endif

    !call cpu_time(ttime1)

    !dt1=ttime1-ttime0+dt1


    SELECT CASE(T%CAS)
    CASE(CASEP1)
       CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,X,K)
    CASE(CASEP2)
       CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)

    CASE(CASE1,CASE2)
       el=>T%PARENT_FIBRE%MAG
       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=X,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,X,1)
             CALL TRACK_FRINGE(EL%C4,X,J=1)
          ELSE
             CALL TRACK_FRINGE(EL%C4,X,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,X,2)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=X,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=X,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=X,J=T%CAS)
       case(KIND10)
          CALL TRACK_FRINGE(EL%TP10,X,T%CAS)
       case(KIND16,KIND20)
          CALL TRACK_FRINGE(EL%K16,X,T%CAS)
       case(KIND17)
          CALL TRACK_FRINGE(EL17=EL%S17,X=X,J=T%CAS)
       case(KIND21)
          CALL TRACK_FRINGE(EL%CAV21,X=X,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)
          CALL ADJUST_WI(EL%WI,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT

    CASE(CASE0)
       if(associated(T%PARENT_FIBRE%MAG%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
       el=>T%PARENT_FIBRE%MAG
       SELECT CASE(EL%KIND)
       CASE(KIND0)
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X)
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X)
       case(KIND3)
          CALL TRACK_SLICE(EL%K3,X)
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X)
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X)
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X)
       case(KIND7)
          CALL TRACK_SLICE(EL%T7,X,t%pos_in_fibre)
       case(KIND8)
          CALL TRACK_SLICE(EL%S8,X)
       case(KIND9)
          CALL TRACK_SLICE(EL%S9,X)
       case(KIND10)
          CALL TRACK_SLICE(EL%TP10,X)
       case(KIND11:KIND14)
          CALL TRACK_SLICE(EL%MON14,X)
       case(KIND15)
          CALL TRACK_SLICE(EL%SEP15,X)
       case(KIND16,KIND20)
          CALL TRACK_SLICE(EL%K16,X)
       case(KIND17)
          CALL TRACK_SLICE(EL%S17,X)
       case(KIND18)
          CALL TRACK_SLICE(EL%RCOL18,X)
       case(KIND19)
          CALL TRACK_SLICE(EL%ECOL19,X)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,T%S(2))
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,T%S(2))
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,T%POS_IN_FIBRE)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT

       ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S


    END SELECT
    T%PARENT_FIBRE%MAG=DEFAULT

  END SUBROUTINE TRACKR_NODE_SINGLE

  SUBROUTINE TRACKR_NODE_SINGLE_orbit(T,X,K,CHARGE)
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT):: T
    REAL(DP),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    INTEGER, TARGET :: CHARGE
    integer, target :: CHARGE1
    INTEGER I

    ! call cpu_time(ttime0)


    !    T%PARENT_FIBRE%MAG=K
    !!   T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR    only standard lattice
    !    if(present(charge))  then
    !!    T%PARENT_FIBRE%MAG%P%CHARGE=>CHARGE            only standard lattice

    !    else
    !       charge1=1
    !      T%PARENT_FIBRE%MAG%P%CHARGE=>CHARGE1
    !    endif

    !call cpu_time(ttime1)

    !dt1=ttime1-ttime0+dt1


    SELECT CASE(T%CAS)
    CASE(CASEP1)
       if(abs(x(1))+abs(x(3))>absolute_aperture.or.(.not.CHECK_MADX_APERTURE)) then
          CHECK_STABLE=.false.
       endif
       if((T%PARENT_FIBRE%mag%mis)) CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,X,K)
    CASE(CASEP2)
       if((T%PARENT_FIBRE%mag%mis)) CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)

    CASE(CASE1,CASE2)
       el=>T%PARENT_FIBRE%MAG
       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=X,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,X,1)
             CALL TRACK_FRINGE(EL%C4,X,J=1)
          ELSE
             CALL TRACK_FRINGE(EL%C4,X,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,X,2)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=X,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=X,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=X,J=T%CAS)
       case(KIND10)
          CALL TRACK_FRINGE(EL%TP10,X,T%CAS)
       case(KIND16,KIND20)
          CALL TRACK_FRINGE(EL%K16,X,T%CAS)
       case(KIND17)
          CALL TRACK_FRINGE(EL17=EL%S17,X=X,J=T%CAS)
       case(KIND21)
          CALL TRACK_FRINGE(EL%CAV21,X=X,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)
          CALL ADJUST_WI(EL%WI,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT

    CASE(CASE0)
       if(associated(T%PARENT_FIBRE%MAG%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAG%p%aperture,X)
       el=>T%PARENT_FIBRE%MAG
       SELECT CASE(EL%KIND)
       CASE(KIND0)
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X)
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X)
       case(KIND3)
          CALL TRACK_SLICE(EL%K3,X)
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X)
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X)
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X)
       case(KIND7)
          CALL TRACK_SLICE(EL%T7,X,t%pos_in_fibre)
       case(KIND8)
          CALL TRACK_SLICE(EL%S8,X)
       case(KIND9)
          CALL TRACK_SLICE(EL%S9,X)
       case(KIND10)
          CALL TRACK_SLICE(EL%TP10,X)
       case(KIND11:KIND14)
          CALL TRACK_SLICE(EL%MON14,X)
       case(KIND15)
          CALL TRACK_SLICE(EL%SEP15,X)
       case(KIND16,KIND20)
          CALL TRACK_SLICE(EL%K16,X)
       case(KIND17)
          CALL TRACK_SLICE(EL%S17,X)
       case(KIND18)
          CALL TRACK_SLICE(EL%RCOL18,X)
       case(KIND19)
          CALL TRACK_SLICE(EL%ECOL19,X)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,T%S(2))
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,T%S(2))
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,T%POS_IN_FIBRE)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT

       ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S


    END SELECT
    !    T%PARENT_FIBRE%MAG=DEFAULT

  END SUBROUTINE TRACKR_NODE_SINGLE_orbit

  SUBROUTINE TRACKP_NODE_SINGLE(T,X,K,CHARGE)
    ! This routines tracks a single thin lens
    ! it is supposed to reproduce plain PTC
    implicit none
    TYPE(INTEGRATION_NODE), TARGET, INTENT(INOUT):: T
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(elementp),pointer :: el
    INTEGER, optional, TARGET :: CHARGE
    integer, target :: CHARGE1
    INTEGER I
    logical(lp) BN2,L
    logical(lp) CHECK_KNOB
    logical(lp), ALLOCATABLE,dimension(:)::AN,BN


    T%PARENT_FIBRE%MAGP=K
    IF(K%PARA_IN ) KNOB=.TRUE.

    T%PARENT_FIBRE%MAGP%P%DIR=>T%PARENT_FIBRE%DIR

    if(present(charge))  then
       T%PARENT_FIBRE%MAGP%P%CHARGE=>CHARGE
    else
       charge1=1
       T%PARENT_FIBRE%MAGP%P%CHARGE=>CHARGE1
    endif
    !    T%PARENT_FIBRE%MAGP%P%CHARGE=>CHARGE


    SELECT CASE(T%CAS)
    CASE(CASEP1)
       CALL TRACK_FIBRE_FRONT(T%PARENT_FIBRE,X,K)
    CASE(CASEP2)
       CALL TRACK_FIBRE_BACK(T%PARENT_FIBRE,X,K)

    CASE(CASE1,CASE2)
       el=>T%PARENT_FIBRE%MAGP
       SELECT CASE(EL%KIND)
       CASE(KIND0:KIND1,KIND3,KIND8:KIND9,KIND11:KIND15,KIND18:KIND19)
       case(KIND2)
          CALL TRACK_FRINGE(EL=EL%K2,X=X,J=T%CAS)
       case(KIND4)
          IF(T%CAS==CASE1) THEN
             CALL ADJUST_TIME_CAV4(EL%C4,X,1)
             CALL TRACK_FRINGE(EL%C4,X,J=1)
          ELSE
             CALL TRACK_FRINGE(EL%C4,X,J=2)
             CALL ADJUST_TIME_CAV4(EL%C4,X,2)
          ENDIF
       case(KIND5)
          CALL TRACK_FRINGE(EL5=EL%S5,X=X,J=T%CAS)
       case(KIND6)
          CALL TRACK_FRINGE(EL6=EL%T6,X=X,J=T%CAS)
       case(KIND7)
          CALL TRACK_FRINGE(EL7=EL%T7,X=X,J=T%CAS)
       case(KIND10)
          CALL TRACK_FRINGE(EL%TP10,X,T%CAS)
       case(KIND16,KIND20)
          CALL TRACK_FRINGE(EL%K16,X,T%CAS)
       case(KIND17)
          CALL TRACK_FRINGE(EL17=EL%S17,X=X,J=T%CAS)
       case(KIND21)
          CALL TRACK_FRINGE(EL%CAV21,X=X,J=T%CAS)
          CALL ADJUST_TIME_CAV_TRAV_OUT(EL%CAV21,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDWIGGLER)
          CALL ADJUST_WI(EL%WI,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       case(KINDPA)
          CALL ADJUST_PANCAKE(EL%PA,X,T%CAS)   ! ONLY DOES SOMETHING IF J==2
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 666
       END SELECT

    CASE(CASE0)
       if(associated(T%PARENT_FIBRE%MAGP%p%aperture)) call CHECK_APERTURE(T%PARENT_FIBRE%MAGP%p%aperture,X)
       el=>T%PARENT_FIBRE%MAGP
       SELECT CASE(EL%KIND)
       CASE(KIND0)
       case(KIND1)
          CALL TRACK_SLICE(EL%D0,X)
       case(KIND2)
          CALL TRACK_SLICE(EL%K2,X)
       case(KIND3)
          CALL TRACK_SLICE(EL%K3,X)
       case(KIND4)
          CALL TRACK_SLICE(EL%C4,X)
       case(KIND5)
          CALL TRACK_SLICE(EL%S5,X)
       case(KIND6)
          CALL TRACK_SLICE(EL%T6,X)
       case(KIND7)
          IF((EL%T7%BN(2)%KIND==3.OR.EL%T7%L%KIND==3).AND.KNOB) THEN
             CALL GETMAT7(EL%T7)                                      ! RECOMPUTES ONLY IF KNOB (SPEED)
          ENDIF
          CALL TRACK_SLICE(EL%T7,X,t%pos_in_fibre)
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
          CALL TRACK_SLICE(EL%S8,X)
       case(KIND9)
          CALL TRACK_SLICE(EL%S9,X)
       case(KIND10)
          IF(KNOB) THEN
             CALL CHECKPOTKNOB(EL%TP10,CHECK_KNOB) ! RECOMPUTES ONLY IF KNOB (SPEED)
             IF(CHECK_KNOB) THEN
                ALLOCATE(AN(EL%TP10%P%NMUL),BN(EL%TP10%P%NMUL))
                DO I=1,EL%TP10%P%NMUL
                   BN(I)=.FALSE.
                   AN(I)=.FALSE.
                   IF(EL%TP10%BN(I)%KIND==3) BN(I)=.TRUE.
                   IF(EL%TP10%AN(I)%KIND==3) AN(I)=.TRUE.
                ENDDO
                call GETANBN(EL%TP10)
             ENDIF
          ENDIF
          CALL TRACK_SLICE(EL%TP10,X)
       case(KIND11:KIND14)
          CALL TRACK_SLICE(EL%MON14,X)
       case(KIND15)
          CALL TRACK_SLICE(EL%SEP15,X)
       case(KIND16,KIND20)
          CALL TRACK_SLICE(EL%K16,X)
       case(KIND17)
          CALL TRACK_SLICE(EL%S17,X)
       case(KIND18)
          CALL TRACK_SLICE(EL%RCOL18,X)
       case(KIND19)
          CALL TRACK_SLICE(EL%ECOL19,X)
       case(KIND21)
          CALL TRACK_SLICE(EL%CAV21,X,T%S(2))
       case(KINDWIGGLER)
          CALL TRACK_SLICE(EL%WI,X,T%S(2))
       case(KINDPA)
          CALL TRACK_SLICE(EL%PA,X,T%POS_IN_FIBRE)
       CASE DEFAULT
          WRITE(6,*) "NOT IMPLEMENTED ",EL%KIND
          stop 999
       END SELECT

       ! CASE(CASE100)  ! FAKE BEAM BEAM CAKE AT SOME S


    END SELECT

    T%PARENT_FIBRE%MAGP=DEFAULT
    ! KNOB IS RETURNED TO THE PTC DEFAULT
    ! NEW STUFF WITH KIND=3
    KNOB=.FALSE.
    ! END NEW STUFF WITH KIND=3

  END SUBROUTINE TRACKP_NODE_SINGLE


  SUBROUTINE TRACK_NODE_SINGLE_FOR_time(B,I,DT,K)
    ! Tracks a single particle "I" of the beam for a time DT
    ! The particle is a location defined by the thin lens B%POS(I)%NODE
    ! and located B%X(I,7) metres in from of that thin lens
    implicit none
    TYPE(INTEGRATION_NODE), POINTER:: T
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP), INTENT(IN) :: DT
    REAL(DP) X(6),XT(6),DT0,YL,DT_BEFORE
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    LOGICAL(LP) END_OF_LINE
    INTEGER I
    END_OF_LINE=.FALSE.

    IF(.NOT.K%TOTALPATH) STOP 451
    IF(B%U(i)) RETURN

    X=BEAM_IN_X(B,I)
    T=>B%POS(I)%NODE
    T%PARENT_FIBRE%MAG=K
    T%PARENT_FIBRE%MAG%P%DIR=>T%PARENT_FIBRE%DIR
    T%PARENT_FIBRE%MAG%P%CHARGE=>B%CHARGE


    DT0=X(6)
    CALL DRIFT_BACK_TO_POSITION(T,B%X(I,7),X)
    B%X(I,7)=ZERO
    DT0=X(6)-DT0

    YL=zero

    DO WHILE(DT0<=DT)
       XT=X
       DT_BEFORE=DT0
       !         WRITE(6,*) " POS ",T%s(1),t%pos_in_fibre
       !         WRITE(6,*) " POS ",T%POS,T%CAS,T%PARENT_FIBRE%MAG%NAME
       CALL TRACK_NODE_SINGLE(T,X,K,B%CHARGE)
       DT0=DT0+(X(6)-XT(6))
       T=>T%NEXT
       IF(.NOT.ASSOCIATED(T%NEXT)) THEN
          END_OF_LINE=.TRUE.
          EXIT
       ENDIF
    ENDDO

    IF(.NOT.END_OF_LINE) THEN
       IF(DT0/=DT) THEN
          B%POS(I)%NODE=>T%PREVIOUS
          X=XT
          DT0=DT-DT_BEFORE
          !           WRITE(6,*) " DT0 ", DT0
          CALL DRIFT_TO_TIME(T,YL,DT0,X)
       ELSE
          B%POS(I)%NODE=>T
       ENDIF
    ELSE


       IF(DT0<DT) THEN
          B%POS(I)%NODE=>T%PREVIOUS
          X=XT
          DT0=DT-DT_BEFORE
          CALL DRIFT_TO_TIME(T,YL,DT0,X)
       ELSE
          B%POS(I)%NODE=>T
       ENDIF

    ENDIF



    CALL X_IN_BEAM(B,X,I,DL=YL)



    B%TIME_INSTEAD_OF_S=.TRUE.


  END SUBROUTINE TRACK_NODE_SINGLE_FOR_time

  SUBROUTINE TRACK_NODE_SINGLE_FOR_MAP(T,B,K)
    ! Tracks a full beam across a map representing several thin lenses

    implicit none
    TYPE(INTEGRATION_NODE),POINTER :: T
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP) X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    INTEGER I

    !    IF(ASSOCIATED(T%BT)) THEN
    !       IF(B%BEAM_BEAM) THEN
    !          CALL BBKICK(b,t)
    !       ENDIF
    !    ENDIF

    DO I=1,B%N
       IF(B%U(i)) CYCLE
       X=BEAM_IN_X(B,I)
       X=X-T%ORBIT
       CALL TRACK(T%TPSA_MAP,X)

       CALL X_IN_BEAM(B,X,I,DL=ZERO,T=T%INTEGRATION_NODE_AFTER_MAP)
       !       B%POS(I)%NODE=>T%INTEGRATION_NODE_AFTER_MAP

    ENDDO

    IF(ASSOCIATED(B%Y)) THEN
       DO I=1,6
          B%Y(I)=B%Y(I)-T%ORBIT(I)
       ENDDO
       CALL TRACK(T%TPSA_MAP,B%Y)
       CALL X_IN_BEAM(B,I=0,T=T%INTEGRATION_NODE_AFTER_MAP)
       !       B%POS(0)%NODE=>T%INTEGRATION_NODE_AFTER_MAP
    ENDIF

    B%TIME_INSTEAD_OF_S=.FALSE.

  END SUBROUTINE TRACK_NODE_SINGLE_FOR_MAP



  SUBROUTINE TRACK_NODE_SINGLE_FOR_S(T,B,K)
    ! Tracks a full beam across a single thin lens
    ! Position of macroparticle put at next thin lens

    implicit none
    TYPE(INTEGRATION_NODE),POINTER :: T
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP) X(6)
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    type(element),pointer :: el
    INTEGER I

    !    IF(ASSOCIATED(T%BT)) THEN
    !       IF(B%BEAM_BEAM) THEN
    !          CALL BBKICK(b,t)
    !       ENDIF
    !    ENDIF

    DO I=1,B%N
       IF(B%U(i)) CYCLE
       X=BEAM_IN_X(B,I)

       CALL TRACK_NODE_SINGLE(T,X,K,B%CHARGE)
       CALL X_IN_BEAM(B,X,I,DL=ZERO,T=T%NEXT)
       !       B%POS(I)%NODE=>T%NEXT

    ENDDO

    IF(ASSOCIATED(B%Y)) THEN
       CALL TRACK_NODE_SINGLE(T,B%Y,K,B%CHARGE)
       CALL X_IN_BEAM(B,I=0,T=T%NEXT)
       !       B%POS(0)%NODE=>T%NEXT
    ENDIF

    B%TIME_INSTEAD_OF_S=.FALSE.

  END SUBROUTINE TRACK_NODE_SINGLE_FOR_S

  SUBROUTINE TRACK_NODE_T(B,DT,K)
    ! Tracks to full beam for a time DT
    ! All the particles are at different locations
    ! Notice that the layout is hidden: this is consistant with time tracking
    ! Magnets are not ontological objects
    ! calls TRACK_NODE_SINGLE_FOR_time

    implicit none
    TYPE(BEAM),INTENT(INOUT):: B
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    real(dp),INTENT(IN):: DT

    INTEGER I
    DO I=1,B%N
       IF(B%U(i)) CYCLE
       call TRACK_NODE_SINGLE_FOR_time(B,I,DT,K)
    ENDDO
    B%TIME_INSTEAD_OF_S=.TRUE.

  END SUBROUTINE TRACK_NODE_T

  SUBROUTINE TRACK_LAYOUT_XV_12( R,V,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE )
    ! Tracks through the thin lens structure R%T of the layout R if it exists.
    ! Several posibilities:
    !1) Pos1 and Pos2 are given :  tracks from thin lens position pos1 to thin lens position pos2
    !2) Thin lens points T1 and T2 are given: Same as above pos1 and pos2 are derived from T1 and T2
    !3) P1 and P2 are fibres: Tracks from the first thin lens of P1 to the first of P2. Results should agree
    !   with plain PTC.
    !4) P1, IN_P1 and P2 , IN_P2 are give: moves to IN_P1 thin lens in P1 and tracks to IN_P2 position in p2
    !  For example, if the input is (P1,1) and (P2,1) the results is same as item #3 (same as plain PTC).
    !4) POS1_FIBRE and POS2_FIBRE are given: same as standard PTC

    ! In all the above cases one can elect to give only the first input. Then it tracks one turn around as in
    ! plain PTC.
    ! interfaced as TRACK_LAYOUT_USING_NODE_S

    ! V is of type three_d_info: position in 3d space is provided

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE(BEAM) B
    INTEGER, OPTIONAL,INTENT(IN):: POS1,POS2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE
    TYPE(FIBRE), OPTIONAL,POINTER :: P1,P2
    TYPE(INTEGRATION_NODE),POINTER, OPTIONAL :: T1,T2
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    REAL(DP) SC,reference_ray(6),x(6)
    type(three_d_info)  v
    TYPE(INTEGRATION_NODE),POINTER:: t,mag_in,mag_out

    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF

    CALL ALLOCATE_BEAM(B,2,MY_FALSE)
    B%X=ZERO
    B%X(1,1:6)=V%X
    B%X(2,1:6)=V%reference_ray
    x=V%X
    reference_ray=V%reference_ray
    CALL TRACK_LAYOUT_ONE_12( R,B,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE )

    V%X=B%X(1,1:6)
    V%reference_ray=B%X(2,1:6)

    IF(.NOT.ASSOCIATED(B%POS(1)%NODE%B)) THEN
       WRITE(6,*) " NO FRAMES IN INTEGRATION NODES "
       STOP 101
    ENDIF
    SC=ONE
    IF(v%SCALE/=zero) SC=v%SCALE
    t=>B%POS(1)%NODE%previous

    V%r0=t%A+(reference_ray(1)-SC*reference_ray(1))*t%ENT(1,1:3)+ SC*X(1)*t%ENT(1,1:3)
    V%r0=v%r0+(reference_ray(3)-SC*reference_ray(3))*t%ENT(2,1:3)+ SC*X(3)*t%ENT(2,1:3)

    V%r=t%B+(V%reference_ray(1)-SC*V%reference_ray(1))*t%EXI(1,1:3)+ SC*V%X(1)*t%EXI(1,1:3)
    V%r=v%r+(V%reference_ray(3)-SC*V%reference_ray(3))*t%EXI(2,1:3)+ SC*V%X(3)*t%EXI(2,1:3)
    mag_in=>t%previous%parent_fibre%t1%next%next
    mag_out=>t%previous%parent_fibre%t2%previous%previous
    v%a=mag_in%a
    v%ent=mag_in%ent
    v%b=mag_in%b
    v%exi=mag_in%exi
    v%o=t%B
    v%mid=t%exi
    v%U=B%U(1:2)



    IF(MAG_IN%PREVIOUS%CAS/=CASE1) STOP 201
    IF(MAG_OUT%NEXT%CAS/=CASE2) STOP 202

    CALL KILL_BEAM(B)

  END SUBROUTINE TRACK_LAYOUT_XV_12

  SUBROUTINE TRACK_LAYOUT_XR_12( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE,S1,S2 )
    ! Tracks through the thin lens structure R%T of the layout R if it exists.
    ! Several posibilities:
    !1) Pos1 and Pos2 are given :  tracks from thin lens position pos1 to thin lens position pos2
    !2) Thin lens points T1 and T2 are given: Same as above pos1 and pos2 are derived from T1 and T2
    !3) P1 and P2 are fibres: Tracks from the first thin lens of P1 to the first of P2. Results should agree
    !   with plain PTC.
    !4) P1, IN_P1 and P2 , IN_P2 are give: moves to IN_P1 thin lens in P1 and tracks to IN_P2 position in p2
    !  For example, if the input is (P1,1) and (P2,1) the results is same as item #3 (same as plain PTC).
    !4) POS1_FIBRE and POS2_FIBRE are given: same as standard PTC

    ! In all the above cases one can elect to give only the first input. Then it tracks one turn around as in
    ! plain PTC.
    ! interfaced as TRACK_LAYOUT_USING_NODE_S

    ! Here a usual trajectory is computed     REAL(DP),INTENT(INOUT):: X(6)

    implicit none
    TYPE (LAYOUT), TARGET :: R
    REAL(DP),INTENT(INOUT):: X(6)
    TYPE(BEAM) B
    INTEGER, OPTIONAL,INTENT(IN):: POS1,POS2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE
    TYPE(FIBRE), OPTIONAL,POINTER :: P1,P2
    TYPE(INTEGRATION_NODE),POINTER, OPTIONAL :: T1,T2
    REAL(DP), OPTIONAL :: S1,S2
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    LOGICAL(LP) U

    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF

    CALL ALLOCATE_BEAM(B,1,MY_FALSE)
    B%X=ZERO
    B%X(1,1:6)=X

    CALL TRACK_LAYOUT_ONE_12( R,B,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE,S1,S2)

    X=B%X(1,1:6)
    U=B%U(1)


    CALL KILL_BEAM(B)

  END SUBROUTINE TRACK_LAYOUT_XR_12

  SUBROUTINE TRACK_LAYOUT_XP_12( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    ! Tracks through the thin lens structure R%T of the layout R if it exists.
    ! Several posibilities:
    !1) Pos1 and Pos2 are given :  tracks from thin lens position pos1 to thin lens position pos2
    !2) Thin lens points T1 and T2 are given: Same as above pos1 and pos2 are derived from T1 and T2
    !3) P1 and P2 are fibres: Tracks from the first thin lens of P1 to the first of P2. Results should agree
    !   with plain PTC.
    !4) P1, IN_P1 and P2 , IN_P2 are give: moves to IN_P1 thin lens in P1 and tracks to IN_P2 position in p2
    !  For example, if the input is (P1,1) and (P2,1) the results is same as item #3 (same as plain PTC).
    !4) POS1_FIBRE and POS2_FIBRE are given: same as standard PTC

    ! In all the above cases one can elect to give only the first input. Then it tracks one turn around as in
    ! plain PTC.
    ! interfaced as TRACK_LAYOUT_USING_NODE_S

    ! Here a usual polymorphic trajectory is computed  TYPE(REAL_8),TARGET,INTENT(INOUT):: X(6)

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE(REAL_8),TARGET,INTENT(INOUT):: X(6)
    TYPE(BEAM) B
    INTEGER, OPTIONAL,INTENT(IN):: POS1,POS2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE
    TYPE(FIBRE), OPTIONAL,POINTER :: P1,P2
    TYPE(INTEGRATION_NODE),POINTER, OPTIONAL :: T1,T2
    TYPE(INTERNAL_STATE) K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    LOGICAL(LP) U

    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF
    CALL ALLOCATE_BEAM(B,1,MY_FALSE)
    B%X=ZERO
    B%X(1,1:6)=X
    B%Y=>X

    CALL TRACK_LAYOUT_ONE_12( R,B,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE )

    U=B%U(0)

    NULLIFY(B%Y)

    CALL KILL_BEAM(B)

  END SUBROUTINE TRACK_LAYOUT_XP_12


  SUBROUTINE TRACK_LAYOUT_ONE_12( R,B,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE,S1,S2 )
    ! Tracks through the thin lens structure R%T of the layout R if it exists.
    ! Several posibilities:
    !1) Pos1 and Pos2 are given :  tracks from thin lens position pos1 to thin lens position pos2
    !2) Thin lens points T1 and T2 are given: Same as above pos1 and pos2 are derived from T1 and T2
    !3) P1 and P2 are fibres: Tracks from the first thin lens of P1 to the first of P2. Results should agree
    !   with plain PTC.
    !4) P1, IN_P1 and P2 , IN_P2 are give: moves to IN_P1 thin lens in P1 and tracks to IN_P2 position in p2
    !  For example, if the input is (P1,1) and (P2,1) the results is same as item #3 (same as plain PTC).
    !4) POS1_FIBRE and POS2_FIBRE are given: same as standard PTC

    ! In all the above cases one can elect to give only the first input. Then it tracks one turn around as in
    ! plain PTC.
    ! interfaced as TRACK_LAYOUT_USING_NODE_S

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE(BEAM),INTENT(INOUT):: B
    INTEGER, OPTIONAL,INTENT(IN):: POS1,POS2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE
    TYPE(FIBRE), OPTIONAL,POINTER :: P1,P2
    TYPE(INTEGRATION_NODE),POINTER, OPTIONAL :: T1,T2
    REAL(DP), OPTIONAL :: S1,S2
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER I1,I2
    TYPE(FIBRE), POINTER :: C

    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF
    I1=-2*R%T%N;I2=-2*R%T%N;
    IF(PRESENT(POS1)) I1=POS1
    IF(PRESENT(POS2)) I2=POS2
    IF(PRESENT(T1))   I1=T1%POS
    IF(PRESENT(T2))   I2=T2%POS
    IF(PRESENT(POS1_FIBRE)) THEN
       call move_to(r,c,POS1_FIBRE)
       I1=C%T1%POS
    ENDIF
    IF(PRESENT(POS2_FIBRE)) THEN
       call move_to(r,c,POS2_FIBRE)
       I2=C%T1%POS
    ENDIF
    IF(PRESENT(P1))   I1=P1%T1%POS
    IF(PRESENT(P2))   I2=P2%T1%POS
    IF(PRESENT(IN_P1))I1=I1-1+IN_P1
    IF(PRESENT(IN_P2))I2=I2-1+IN_P2
    IF(I1<=0) THEN
       WRITE(6,*) " I1 AND I2 =   ",I1,I2
       WRITE(6,*) " CHECK INPUT"
       RETURN
    ENDIF

    IF(I2>0) THEN
       IF(I2<I1) I2=I2+R%T%N
       CALL TRACK_NODE_LAYOUT_S( R%T,B,I1,I2,K )
    ELSEIF(I1>0) THEN
       CALL TRACK_NODE_LAYOUT_S( R%T,B,I1,K )
    ELSEIF(PRESENT(S1).AND.PRESENT(S2)) THEN
       CALL TRACK_LAYOUT_S12( R,B,K,S1,S2 )
    ELSE
       CALL TRACK_LAYOUT_S12( R,B,K,S1 )
    ENDIF
    B%TIME_INSTEAD_OF_S=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_ONE_12

  SUBROUTINE TRACK_LAYOUT_S12( R,B,K,S1,S2 )
    ! Tracks through the thin lens structure from position S1 to position S2 (defined as the S(3) variables
    ! of the thin lens.
    ! The final position is stored as in time tracking but is obviously the same for all the particles.
    ! interfaced as TRACK_LAYOUT_USING_NODE_S

    implicit none
    TYPE (LAYOUT), TARGET :: R
    TYPE(BEAM),INTENT(INOUT):: B
    REAL(DP),INTENT(IN) :: S1
    REAL(DP),OPTIONAL :: S2
    REAL(DP) DS1,DS2,SFINAL
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER I1,I2,NTURN,I
    TYPE(INTEGRATION_NODE), POINTER :: LAST
    IF(.NOT.ASSOCIATED(R%T)) THEN
       WRITE(6,*) " NO THIN LAYOUT: TRACKING IMPOSSIBLE "
       RETURN
    ENDIF
    DS2=ZERO
    LAST=>R%T%LAST
    CALL move_to_s( R%T,s1,LAST,i1,ds1 )
    NTURN=0

    I2=-2*R%T%N;
    IF(I1<=0) THEN
       WRITE(6,*) " I1  =   ",I1
       WRITE(6,*) " CHECK INPUT FOR S1"
       RETURN
    ENDIF

    IF(DS1/=ZERO) THEN
       CALL DRIFT_BEAM_BACK_TO_POSITION(LAST,DS1,B)
    ENDIF

    IF(PRESENT(S2)) THEN
       CALL move_to_s( R%T,s2,LAST,i2,ds2 )
       IF(R%CLOSED) THEN
          NTURN=INT((S2-DS2-S1+DS1)/R%T%END%S(3))
       ENDIF
       SFINAL=S2
    ELSE
       SFINAL=S1
    ENDIF

    DO I=1,NTURN
       CALL TRACK_NODE_LAYOUT_S( R%T,B,I1,K )
    ENDDO

    IF(I2>0) THEN
       if(i2<i1) then
          i2=r%t%n+i2
       endif
       CALL TRACK_NODE_LAYOUT_S( R%T,B,I1,I2,K )
    ELSE
       CALL TRACK_NODE_LAYOUT_S( R%T,B,I1,K )
    ENDIF

    IF(DS2/=ZERO) THEN
       DS2=-DS2
       CALL DRIFT_BEAM_BACK_TO_POSITION(LAST,DS2,B)
    ENDIF
    B%X(:,7:7)=SFINAL
    B%TIME_INSTEAD_OF_S=.FALSE.
  END SUBROUTINE TRACK_LAYOUT_S12



  SUBROUTINE TRACK_NODE_LAYOUT_S12( R,B,I1,I2,K) ! T
    implicit none
    TYPE (NODE_LAYOUT), TARGET :: R
    TYPE (INTEGRATION_NODE), POINTER ::  C
    TYPE(BEAM),INTENT(INOUT):: B
    INTEGER, INTENT(IN):: I1,I2
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER I,J,J2


    call move_to_INTEGRATION_NODE(r,c,I1)


    if(i2>i1) then
       J=I1

       DO  WHILE(J<I2.AND.ASSOCIATED(C))

          IF(C%USE_TPSA_MAP) THEN !1

             CALL TRACK_NODE_SINGLE_FOR_MAP(C,B,K)
             J2=C%INTEGRATION_NODE_AFTER_MAP%POS
             IF(J2<C%POS) THEN
                J2=J2+R%N
             ENDIF
             IF(J2>I2) STOP 481
             C=>C%INTEGRATION_NODE_AFTER_MAP
             J=J2

          ELSE  !1
             CALL TRACK_NODE_SINGLE_FOR_S(C,B,K)
             C=>C%NEXT
             J=J+1

          ENDIF  !1
       ENDDO
    ELSEIF(I1>I2) THEN
       WRITE(6,*) " BACKWARDS FORBIDDEN IN TRACK_NODE_LAYOUT"
       STOP 666
    ENDIF
    B%TIME_INSTEAD_OF_S=.FALSE.
  end SUBROUTINE TRACK_NODE_LAYOUT_S12

  SUBROUTINE TRACK_NODE_LAYOUT_S1( R,B,I1,K ) ! TRACKS LAYOUT WITHOUT ANY COLLECTIVE FRIVOLITES
    implicit none
    TYPE (NODE_LAYOUT), TARGET :: R
    TYPE(BEAM),INTENT(INOUT):: B
    INTEGER, INTENT(IN):: I1
    TYPE(INTERNAL_STATE)  K
    !    TYPE(INTERNAL_STATE), INTENT(IN) :: K
    INTEGER II1,II2

    II1=I1
    IF(R%CLOSED) THEN
       II2=II1+R%N
    ELSE
       II2=R%N+1
    ENDIF

    CALL TRACK_NODE_LAYOUT_s(R,B,II1,II2,k)

    B%TIME_INSTEAD_OF_S=.FALSE.

  end SUBROUTINE TRACK_NODE_LAYOUT_S1

  !  STUFF ABOUT LIST AND STRUCTURES

  SUBROUTINE MAKE_NODE_LAYOUT( R) !
    ! Creates the thin layout and puts in R%T by calling MAKE_NODE_LAYOUT_2
    ! At this point large patches would be problematic; i.e. |d(3)|>>>0 .
    implicit none
    TYPE (LAYOUT), TARGET :: R

    call MAKE_NODE_LAYOUT_2( R,R%T )

  end SUBROUTINE MAKE_NODE_LAYOUT



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
    TYPE(INTEGRATION_NODE), POINTER :: T1,T2

    CASE_NAME(CASEP1)="THE ENTRANCE PATCH"
    CASE_NAME(CASEP2)="THE EXIT PATCH  "
    CASE_NAME(CASE1)= "THE ENTRANCE FRINGE"
    CASE_NAME(CASE2)="THE EXIT FRINGE"
    CASE_NAME(CASE0)="A STEP OF INTEGRATION BODY"
    !    CASE_NAME(CASE3)="BEAM BEAM PANCAKE"

    if(associated(L)) then
       CALL kill_NODE_LAYOUT(L)
       DEALLOCATE(L);
       NULLIFY(L);
    endif

    allocate(L)
    CALL Set_Up_NODE_LAYOUT( L )
    S=zero
    SL=ZERO  !  INTEGRATION LENGTH
    P=>R%START
    k=1
    DO I=1,R%N

       TEAPOT_LIKE=0
       IF(P%MAG%P%B0/=ZERO) TEAPOT_LIKE=1
       IF(P%MAG%KIND==KIND16.OR.P%MAG%KIND==KIND16)TEAPOT_LIKE=0
       IF(P%MAG%KIND==KIND0.AND.P%MAG%P%NST/=1) THEN
          WRITE(6,*) "MARKER SHOULD HAVE NST=1 OTHERWISE PROBLEMS "
          WRITE(6,*) "WILL OCCUR WITH THE WORM AND THE NODE_LAYOUT SURVEY "
          STOP 500
       ENDIF
       IF(P%DIR==1) THEN
          LI=zero;
       ELSE
          LI=P%MAG%L;
       ENDIF
       DL=P%DIR*P%MAG%L/P%MAG%P%NST
       DLD=P%MAG%P%LD/P%MAG%P%NST
       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;    ! s(1) total ld
       T1=>L%END                                                    ! s(2) local integration distance
       ! s(3) total integration distance
       L%END%CAS=CASEP1                                             ! s(4) end of step =  DL
       L%END%pos_in_fibre=1
       L%END%pos=k;k=k+1;
       L%END%PARENT_NODE_LAYOUT=>L
       L%END%PARENT_FIBRE=>P

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       L%END%CAS=CASE1
       L%END%pos_in_fibre=2
       L%END%pos=k;k=k+1;
       L%END%PARENT_NODE_LAYOUT=>L
       L%END%PARENT_FIBRE=>P

       DO J=1,P%MAG%P%NST
          CALL APPEND_EMPTY_THIN( L )
          L%END%TEAPOT_LIKE=TEAPOT_LIKE
          L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=DL;
          L%END%CAS=CASE0
          L%END%pos_in_fibre=J+2
          L%END%pos=k;k=k+1;
          L%END%PARENT_NODE_LAYOUT=>L
          L%END%PARENT_FIBRE=>P
          S=S+DLD
          LI=LI+DL
          SL=SL+P%DIR*DL
       ENDDO

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       L%END%CAS=CASE2
       L%END%pos_in_fibre=P%MAG%P%NST+3
       L%END%pos=k;k=k+1;
       L%END%PARENT_FIBRE=>P

       CALL APPEND_EMPTY_THIN( L )
       L%END%TEAPOT_LIKE=TEAPOT_LIKE
       L%END%S(1)=S;L%END%S(2)=LI;L%END%S(3)=SL;L%END%S(4)=zero;
       L%END%CAS=CASEP2
       L%END%pos_in_fibre=P%MAG%P%NST+4
       L%END%pos=k;k=k+1;
       L%END%PARENT_NODE_LAYOUT=>L
       L%END%PARENT_FIBRE=>P
       T2=>L%END

       P%T1=>T1
       P%T2=>T2

       P=>P%NEXT
    ENDDO
    L%N=k-1

    L%PARENT_LAYOUT=>R

    IF(R%CLOSED) THEN
       l%closed=.true.
       CIRCULAR=.TRUE.
       CALL RING_L_THIN(L,CIRCULAR)
    ENDIF

    call stat_NODE_LAYOUT(l)

  END SUBROUTINE MAKE_NODE_LAYOUT_2


  SUBROUTINE stat_NODE_LAYOUT( L )
    implicit none
    TYPE (NODE_LAYOUT), pointer :: L

    WRITE(6,*)  " PARENT LAYOUT NAME :", L%PARENT_LAYOUT%NAME(1:len_trim(L%PARENT_LAYOUT%NAME))
    WRITE(6,*) " NUMBER OF ORIGINAL LAYOUT ELEMENTS :", L%PARENT_LAYOUT%N
    WRITE(6,*) " NUMBER OF THIN OBJECTS :", L%N
    WRITE(6,*) " TOTAL IDEAL LENGTH OF STRUCTURE :", L%END%S(1)
    WRITE(6,*) " TOTAL INTEGERATION LENGTH OF STRUCTURE (mad8 style survey) :", L%END%S(3)

  end SUBROUTINE stat_NODE_LAYOUT


  SUBROUTINE DRIFT_TO_TIME(T,YL,DT,X)
    ! Drifts to a given time using either a regular or TEAPOT drift. (Cylindrical coordinate drift)
    ! The amount of the drift YL is computed to achieve a time DT.

    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(INOUT):: YL
    real(dp),INTENT(IN):: DT
    TYPE(INTEGRATION_NODE), pointer :: T
    TYPE(magnet_chart), pointer :: p
    real(dp) XN(6),PZ,PT
    real(dp)  A,b,R

    P=>T%PARENT_FIBRE%MAG%P

    IF(P%TIME) then
       B=P%BETA0
    ELSE
       B=one
    ENDIF


    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then


       PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
       R=one/P%B0
       A=(X(1)+R)*(one/b+x(5))/PZ
       A=ATAN(DT/(A+DT*X(2)/PZ) )
       YL=A/P%B0
       PT=one-X(2)*TAN(A)/PZ
       XN(1)=(X(1)+R)/COS(A)/PT-R
       XN(2)=X(2)*COS(A)+SIN(A)*PZ
       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    else


       PZ=ROOT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)

       YL=DT/(one/b+X(5))*PZ

       X(1)=X(1)+YL*X(2)/PZ
       X(3)=X(3)+YL*X(4)/PZ
       X(6)=X(6)+YL*(one/b+X(5))/PZ



    endif



  END SUBROUTINE DRIFT_TO_TIME

  SUBROUTINE DRIFTr_BACK_TO_POSITION(T,YL,X)
    ! This is a regular drift
    ! It is used in time tracking to project back to the beginning of the thin lens
    ! and it is used in S tracking to drift in the middle of a step.

    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: YL
    TYPE(INTEGRATION_NODE), pointer :: T
    TYPE(magnet_chart), pointer :: p
    real(dp) XN(6),PZ,PT
    real(dp)  A,b,R

    P=>T%PARENT_FIBRE%MAG%P

    IF(P%TIME) then
       B=P%BETA0
    ELSE
       B=one
    ENDIF


    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then



       PZ=ROOT(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
       R=one/P%B0
       A=-YL*P%B0

       PT=one-X(2)*TAN(A)/PZ
       !       XN(1)=(X(1)+R)/COS(A)/PT-R
       XN(1)=(X(1)+R*(two*sin(a/two)**2+X(2)*sin(A)/PZ))/COS(A)/PT
       XN(2)=X(2)*COS(A)+SIN(A)*PZ
       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
       !          WRITE(6,*) "XN(6)-X(6)-DT , DT",DT,XN(6)-X(6)-DT
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    else
       !       CALL DRIFT(YL,DL,P%beta0,1,P%EXACT,P%TIME,X)
       PZ=ROOT(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)


       X(1)=X(1)-YL*X(2)/PZ
       X(3)=X(3)-YL*X(4)/PZ
       X(6)=X(6)-YL*(one/b+X(5))/PZ

    endif



  END SUBROUTINE DRIFTr_BACK_TO_POSITION


  SUBROUTINE DRIFTp_BACK_TO_POSITION(T,YL,X)
    ! This is a regular drift
    ! It is used in time tracking to project back to the beginning of the thin lens
    ! and it is used in S tracking to drift in the middle of an step.

    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    real(dp),INTENT(IN):: YL
    TYPE(INTEGRATION_NODE), pointer :: T
    TYPE(magnet_chart), pointer :: p
    type(real_8) XN(6),PZ,PT
    real(dp)  A,b,R

    P=>T%PARENT_FIBRE%MAG%P

    IF(P%TIME) then
       B=P%BETA0
    ELSE
       B=one
    ENDIF

    call alloc(xn,6); call alloc(pz,pt);

    if(P%B0/=zero.AND.T%TEAPOT_LIKE==1) then



       PZ=sqrt(one+two*x(5)/b+X(5)**2-X(2)**2-X(4)**2)
       R=one/P%B0
       A=-YL*P%B0

       PT=one-X(2)*TAN(A)/PZ
       !       XN(1)=(X(1)+R)/COS(A)/PT-R
       XN(1)=(X(1)+R*(two*sin(a/two)**2+X(2)*sin(A)/PZ))/COS(A)/PT
       XN(2)=X(2)*COS(A)+SIN(A)*PZ
       XN(3)=X(3)+X(4)*(X(1)+R)*TAN(A)/PZ/PT
       XN(6)=X(6)+(X(1)+R)*TAN(A)/PZ/PT*(one/b+x(5))
       !          WRITE(6,*) "XN(6)-X(6)-DT , DT",DT,XN(6)-X(6)-DT
       X(1)=XN(1)
       X(2)=XN(2)
       X(3)=XN(3)
       X(6)=XN(6)
    else
       !       CALL DRIFT(YL,DL,P%beta0,1,P%EXACT,P%TIME,X)
       PZ=sqrt(one+two*X(5)/b+x(5)**2-X(2)**2-X(4)**2)


       X(1)=X(1)-YL*X(2)/PZ
       X(3)=X(3)-YL*X(4)/PZ
       X(6)=X(6)-YL*(one/b+X(5))/PZ

    endif

    call kill(xn,6); call kill(pz,pt);



  END SUBROUTINE DRIFTp_BACK_TO_POSITION

  SUBROUTINE DRIFT_BEAM_BACK_TO_POSITION(T,YL,B)
    ! This is a regular drift
    ! It is used in time tracking to project back to the beginning of the thin lens
    ! and it is used in S tracking to drift in the middle of an step.
    IMPLICIT NONE
    real(dp)  X(6)
    TYPE(BEAM), INTENT(INOUT) :: B
    real(dp),INTENT(IN):: YL
    TYPE(INTEGRATION_NODE), pointer :: T
    INTEGER I

    DO I=1,B%N
       IF(B%U(i)) CYCLE
       X=BEAM_IN_X(B,I)

       CALL DRIFT_BACK_TO_POSITION(T,YL,X)

       !       B%POS(I)%NODE=>T

       CALL X_IN_BEAM(B,X,I,DL=ZERO,T=T)


    ENDDO

    IF(ASSOCIATED(B%Y)) THEN
       CALL DRIFT_BACK_TO_POSITION(T,YL,B%Y)
       CALL X_IN_BEAM(B,I=0,T=T)
       !       B%POS(0)%NODE=>T%NEXT    LOOKS LIKE ERROR
    ENDIF



  END SUBROUTINE DRIFT_BEAM_BACK_TO_POSITION

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
    integer k,my_start,ic,j
    real(dp) x(6),ent(3,3),a(3)

    CALL  allocate_node_frame( R)

    call survey(r)

    CALL ALLOC(vers,r)
    C=>r%START
    CALL XFRAME(vers%E,C%chart%f%ent,C%chart%f%A,-7)  ! initializes the survey part of worm
    vers%E%L(-1)=0.d0 !Starts beam line at z=0   fake distance along ld for cheap work

    do k=1,r%n
       x=zero
       CALL TRACK(r,x,k,k+1,default,vers)


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


  end  subroutine fill_survey_data_in_NODE_LAYOUT

  ! BEAM STUFF

  subroutine create_beam(B,N,CUT,SIG,T)
    USE gauss_dis
    implicit none
    INTEGER N,I,J
    REAL(DP) CUT,SIG,X
    TYPE(BEAM) B
    TYPE (INTEGRATION_NODE),optional,target::  T

    IF(.NOT.ASSOCIATED(B%N)) THEN
       CALL ALLOCATE_BEAM(B,N)
    ELSEIF(B%N/=N) THEN
       CALL KILL_BEAM(B)
       CALL ALLOCATE_BEAM(B,N)
    ENDIF

    DO I=1,N
       DO J=1,6
          CALL GRNF(X,cut)
          B%X(I,J)=X*SIG
       ENDDO
    ENDDO


    if(present(t)) then
       DO I=1,N
          if(associated(B%POS(I)%NODE))then
             B%POS(I)%NODE=>T
          endif
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
    Y=ZERO
    IF(.not.PRESENT(A)) THEN

       DO I=1,N
          DO J=1,6
             CALL GRNF(X,cut)
             B%X(I,J)=X*SIG(J)
          ENDDO
          B%X(I,7)=ZERO
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
             Y(2*j-1)=X*sqrt(SIG(j)/two)
             CALL GRNF(X,cut)
             Y(2*j)=X*sqrt(SIG(j)/two)
          ENDDO
          y=monkey*Y
          B%X(I,1:C_%ND2)=y(1:c_%nd2)

          DO J=C_%ND2+1,6
             CALL GRNF(X,cut)
             B%X(I,J)=X*SIG(J)
          ENDDO
          B%X(I,7)=ZERO
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
    B2%CHARGE=B1%CHARGE
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
    TYPE(INTEGRATION_NODE),POINTER::T
    TYPE(FIBRE),POINTER::F

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
    TYPE(INTEGRATION_NODE),POINTER::T
    TYPE(FIBRE),POINTER::F

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

  subroutine stat_beam_raw(B,n,MF)
    implicit none
    INTEGER i,j,k,mf,NOTlost,N
    TYPE(BEAM), INTENT(IN):: B
    TYPE(INTEGRATION_NODE),POINTER::T
    TYPE(FIBRE),POINTER::F
    real(dp), allocatable :: av(:,:)
    real(dp) em(2),beta(2)
    allocate(av(n,n))
    av=zero
    notlost=0
    DO K=1,b%n
       IF(.not.B%U(K)) THEN
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
    em(1)=two*sqrt(av(1,1)*av(2,2)-av(1,2)**2)
    em(2)=two*sqrt(av(3,3)*av(4,4)-av(3,4)**2)
    beta(1)=two*av(1,1)/em(1)
    beta(2)=two*av(3,3)/em(2)

    write(mf,*) " average arrays "
    write(mf,*) "betas ",beta
    write(mf,*) "emittances ",em
    if(mf/=6) then
       write(6,*) " average arrays "
       write(6,*) "betas ",beta
       write(6,*) "emittances ",em
    endif
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
    IF(B%TIME_INSTEAD_OF_S) THEN
       WRITE(MF,*) "____________________________ TIME TRACKED BEAM __________________________________"
    ELSE
       WRITE(MF,*) "_________________ POSITION TRACKED BEAM (AS IN PTC PROPER)_______________________"
    ENDIF

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
          IF(B%TIME_INSTEAD_OF_S) THEN
             WRITE(MF,*) " TIME AND POSITION AFTER THIN SLICE = ",B%X(K,6:7)
          ELSE
             WRITE(MF,*) " TIME AND POSITION  = ",B%X(K,6:7)
          ENDIF
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
    NULLIFY(B%Y)
    NULLIFY(B%X)
    NULLIFY(B%U)
    NULLIFY(B%POS)
    NULLIFY(B%CHARGE)
    NULLIFY(B%TIME_INSTEAD_OF_S)
    NULLIFY(B%SIGMA)
    NULLIFY(B%DX,B%ORBIT)
    NULLIFY(B%BBPAR,B%BEAM_BEAM,B%BBORBIT)
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
    v%a=zero
    v%b=zero
    v%o=zero
    v%ent=global_frame
    v%exi=global_frame
    v%mid=global_frame
    v%reference_ray=zero
    v%r0=zero
    v%r=zero
    v%x=zero
    v%scale=one
    v%u=my_false
    v%wx=0.1_dp
    v%wy=0.1_dp

  end subroutine alloc_three_d_info

  SUBROUTINE ALLOCATE_BEAM(B,N,POLYMORPH)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    LOGICAL(LP), OPTIONAL, INTENT(IN) :: POLYMORPH
    INTEGER , INTENT (IN) :: N
    INTEGER I

    ALLOCATE(B%N,B%LOST)

    B%N=N
    B%LOST=0
    NULLIFY(B%Y)
    IF(PRESENT(POLYMORPH)) THEN
       IF(POLYMORPH) then
          ALLOCATE(B%Y(6))
          CALL ALLOC(B%Y)
       endif
    ENDIF
    ALLOCATE(B%X(N,7))
    ALLOCATE(B%U(0:N))
    ALLOCATE(B%POS(0:N))
    ALLOCATE(B%SIGMA(6))
    ALLOCATE(B%DX(3))
    ALLOCATE(B%ORBIT(6))
    ALLOCATE(B%BBPAR,B%BEAM_BEAM,B%BBORBIT)
    DO I=0,N
       NULLIFY(B%POS(i)%NODE)
    ENDDO
    ALLOCATE(B%CHARGE)
    ALLOCATE(B%TIME_INSTEAD_OF_S)

    B%X  = ZERO
    B%U  = .FALSE.
    B%CHARGE=1
    B%TIME_INSTEAD_OF_S=.FALSE.

    B%SIGMA=ZERO
    B%DX=ZERO
    B%BBPAR=ZERO
    B%ORBIT=ZERO
    B%BEAM_BEAM=MY_FALSE
    B%BBORBIT=MY_FALSE
  END SUBROUTINE ALLOCATE_BEAM

  SUBROUTINE KILL_BEAM(B)
    IMPLICIT NONE
    TYPE(BEAM) , INTENT (INOUT) :: B
    IF(ASSOCIATED(B%Y)) THEN
       CALL KILL(B%Y)
       DEALLOCATE(B%Y)
    ENDIF
    IF(ASSOCIATED(B%N)) THEN
       DEALLOCATE(B%N,B%LOST,B%X,B%U,B%POS,B%CHARGE,B%TIME_INSTEAD_OF_S)
       DEALLOCATE(B%SIGMA,B%DX,B%BBPAR,B%ORBIT,B%BEAM_BEAM,B%BBORBIT)
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

  SUBROUTINE TRACK_NODE_LAYOUT_FLAG_R(R,X,I1,I2,k) ! Tracks double from i1 to i2 in state k
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: R
    real(dp), INTENT(INOUT):: X(6)
    TYPE(INTERNAL_STATE) K
    INTEGER, INTENT(IN):: I1,I2
    INTEGER J,i22
    TYPE (INTEGRATION_NODE), POINTER :: C


    CALL RESET_APERTURE_FLAG

    CALL move_to_INTEGRATION_NODE( R%T,C,I1 )



    if(i2>=i1) then
       i22=i2
    else
       i22=r%T%n+i2
    endif

    J=I1

    DO  WHILE(J<I22.AND.ASSOCIATED(C))

       CALL TRACK_NODE_SINGLE(C,X,K,R%CHARGE)

       C=>C%NEXT
       J=J+1
    ENDDO

    if(c_%watch_user) ALLOW_TRACKING=.FALSE.

  END SUBROUTINE TRACK_NODE_LAYOUT_FLAG_R


end module ptc_multiparticle
