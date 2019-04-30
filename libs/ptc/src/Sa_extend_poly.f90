!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_extend_poly
  USE tree_element_MODULE
  IMPLICIT NONE
  public
  logical(lp), target :: ALWAYS_knobs=.false.
  type(real_8) e_muon_scale,a_spin_scale

  ! LD: 22.03.2019 (see Sc_euclidean.f90, Sh_def_kinf.f90 and Sr_spin.f90)
  character(len=150) :: ELEM_NAME = "UNKNOWN"
  integer            :: MAPDUMP = 0 ! 0 or 1 enable/disable PRTP

CONTAINS

  ! LD: 03.04.2019
  SUBROUTINE PRTP1(S, X)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN):: S
    TYPE(REAL_8), INTENT(IN):: X

    ! cancel all PRTP
    if (MAPDUMP .eq. 0) return

    ! @@ + elem + func + 7 columns
    WRITE(*, '(a,a15,a,a15,7E25.16)') '@@ ', ELEM_NAME, ' ', S, X.sub.'000000'&
                              , X.sub.'100000', X.sub.'010000', X.sub.'001000'&
                              , X.sub.'000100',-X.sub.'000001', X.sub.'000010'
  END SUBROUTINE PRTP1

  ! LD: 22.03.2019
  SUBROUTINE PRTP(S, X)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN):: S
    TYPE(REAL_8), OPTIONAL, INTENT(IN):: X(6)

    ! cancel all PRTP
    if (MAPDUMP .eq. 0) return

    ! special case: display only string without X
    if (.not. PRESENT(X)) then
      WRITE(*, '(a,a)') '@@ ', S
      return
    endif

    ! @@ + elem + func + 6 columns
    if (MAPDUMP .eq. 1) then
      WRITE(*, '(a,a15,a,a15,6E25.16)') '@@ ', ELEM_NAME, ' ', S &
        , X(1).sub.'000000', X(2).sub.'000000', X(3).sub.'000000', X(4).sub.'000000',-X(6).sub.'000000', X(5).sub.'000000'
      return
    endif

    ! @@ + elem + func + 42 columns
    WRITE(*, '(a,a15,a,a15,42E25.16)') '@@ ', ELEM_NAME, ' ', S &
      , X(1).sub.'000000', X(2).sub.'000000', X(3).sub.'000000', X(4).sub.'000000',-X(6).sub.'000000', X(5).sub.'000000'&
      , X(1).sub.'100000', X(1).sub.'010000', X(1).sub.'001000', X(1).sub.'000100',-X(1).sub.'000001', X(1).sub.'000010'&
      , X(2).sub.'100000', X(2).sub.'010000', X(2).sub.'001000', X(2).sub.'000100',-X(2).sub.'000001', X(2).sub.'000010'&
      , X(3).sub.'100000', X(3).sub.'010000', X(3).sub.'001000', X(3).sub.'000100',-X(3).sub.'000001', X(3).sub.'000010'&
      , X(4).sub.'100000', X(4).sub.'010000', X(4).sub.'001000', X(4).sub.'000100',-X(4).sub.'000001', X(4).sub.'000010'&
      ,-X(6).sub.'100000',-X(6).sub.'010000',-X(6).sub.'001000',-X(6).sub.'000100', X(6).sub.'000001',-X(6).sub.'000010'&
      , X(5).sub.'100000', X(5).sub.'010000', X(5).sub.'001000', X(5).sub.'000100',-X(5).sub.'000001', X(5).sub.'000010'
  END SUBROUTINE PRTP

  SUBROUTINE ANALYSE_APERTURE_FLAG(I,R)
    IMPLICIT NONE
    INTEGER I,B,K
    INTEGER :: R(:)

    K=I
    B=1
    r=-1
    DO WHILE (K>0.AND.B<=SIZE(R))
       R(B)=MOD(K,2)
       IF(MOD(K,2)==1) THEN
          K=(K-1)/2
       ELSE
          K=K/2
       ENDIF
       B=B+1
    ENDDO

  END   SUBROUTINE ANALYSE_APERTURE_FLAG




  REAL(DP) FUNCTION  SINEHX_X(X) ! REPLACES SINH(X)/X
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       sinehx_x=1.0_dp
       return
    endif

    IF((ABS(X)>hyperbolic_aperture).AND.ROOT_CHECK) THEN
       SINEHX_X=0.0_dp
       CHECK_STABLE=.FALSE.
       messagelost="Sa_extend_poly.f90 SINEXHX_X : argument out of range" !CERN
    ELSEIF(ABS(X)<=hyperbolic_aperture) THEN
       sinehx_x = sinhx_x(x)
    ELSE      !  IF X IS NOT A NUMBER
       sinehx_x=1.0_dp
       CHECK_STABLE=.FALSE.
       messagelost="Sa_extend_poly.f90 SINEXHX_X : should never happen" !CERN
    ENDIF

  END FUNCTION SINEHX_X


  ! Some polymorphism


  ! End of Some polymorphism

end module S_extend_poly



