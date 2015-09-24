!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_extend_poly
  USE tree_element_MODULE
  IMPLICIT NONE
  public
  logical(lp), target :: ALWAYS_knobs=.false.
  type(real_8) e_muon_scale,a_spin_scale



CONTAINS




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
    ELSEIF(ABS(X)<=hyperbolic_aperture) THEN
       sinehx_x = sinhx_x(x)
    ELSE      !  IF X IS NOT A NUMBER
       sinehx_x=1.0_dp
       CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION SINEHX_X


  ! Some polymorphism


  ! End of Some polymorphism

end module S_extend_poly



