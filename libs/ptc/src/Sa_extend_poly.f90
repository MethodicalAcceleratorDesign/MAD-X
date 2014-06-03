!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_extend_poly
  USE tree_element_MODULE
  IMPLICIT NONE
  public
  integer,private,parameter::ndd=6
  private PRINTenv,env_8map,env_8benv
  private ALLOCenv6,killenv6,REAL6env_8,env_8t,tenv_8
  logical(lp), target :: ALWAYS_knobs=.false.
  type(real_8) e_muon_scale
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE ENV_8MAP
     MODULE PROCEDURE REAL6env_8
     MODULE PROCEDURE env_8t
     MODULE PROCEDURE tenv_8
     MODULE PROCEDURE env_8benv
  END  INTERFACE


  INTERFACE PRINT
     MODULE PROCEDURE PRINTenv
  END  INTERFACE

  INTERFACE DAPRINT
     MODULE PROCEDURE PRINTenv
  END  INTERFACE



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


  SUBROUTINE  printenv(S1,mf)
    implicit none
    type (env_8),INTENT(INout)::S1(ndd)
    integer        mf,i,j
    do i=1,ndd
       do j=1,ndd
          write(mf,*) "Sigma0 ",i,j
          call print(s1(i)%sigma0(j),mf)
       enddo
    enddo

    do i=1,ndd
       do j=1,ndd
          write(mf,*) "Sigmaf ",i,j
          call print(s1(i)%sigmaf(j),mf)
       enddo
    enddo

    write(mf,*) "Map "
    do i=1,ndd
       call print(s1(i)%v,mf)
    enddo

    do i=1,ndd
       do j=1,ndd
          write(mf,*) "dB ",i,j
          call print(s1(i)%e(j),mf)
       enddo
    enddo


  END SUBROUTINE printenv

  SUBROUTINE  allocenv6(S1)
    implicit none
    type (env_8),INTENT(INout)::S1(ndd)

    call alloc(s1,ndd)

  END SUBROUTINE allocenv6

  SUBROUTINE  killenv6(S1)

    implicit none
    type (env_8),INTENT(INout)::S1(ndd)

    call kill(s1,ndd)

  END SUBROUTINE killenv6






  SUBROUTINE  REAL6env_8(S2,S1)
    implicit none
    real(dp),INTENT(inout)::S2(ndd)
    type (env_8),INTENT(in)::S1(ndd)
    integer i


    do i=1,ndd
       s2(i)=s1(i)%v
    enddo
  END SUBROUTINE REAL6env_8

  SUBROUTINE  env_8map(S1,S2)
    implicit none
    type (damap),INTENT(in)::S2
    type (env_8),INTENT(inOUT)::S1(ndd)
    integer i


    do i=1,ndd
       s1(i)%v=s2%v(i)
    enddo
  END SUBROUTINE env_8map

  SUBROUTINE  env_8benv(S1,S2)
    implicit none
    type (beamenvelope),INTENT(in)::S2
    type (env_8),INTENT(inOUT)::S1(ndd)

    s1=s2%sij0

  END SUBROUTINE env_8benv

  SUBROUTINE  env_8t(S1,S2)
    implicit none
    type (taylor),INTENT(in)::S2
    type (env_8),INTENT(inOUT)::S1(ndd)
    integer i,j
    character(6) ind,ind0
    real(dp) fac

    ind0='000000'
    do i=1,ndd
       do j=1,ndd
          ind  = ind0
          if(i==j) then
             ind(i:i)='2'
             fac=1.0_dp
          else
             fac=0.5_dp
             ind(i:i)='1'
             ind(j:j)='1'
          endif
          s1(i)%sigma0(j)=(s2.par.ind)*fac
          !        s1(i)%e(j)=(s2.par.ind)*fac
       enddo
    enddo
  END SUBROUTINE env_8t

  SUBROUTINE  tenv_8(S2,S1)
    implicit none
    type (taylor),INTENT(inout)::S2
    type (env_8),INTENT(in)::S1(ndd)
    type(damap) id
    type(real_8) x(6),s
    integer i,j

    call alloc(id)
    call alloc(x,6)
    call alloc(s)
    id=1
    x=id
    do i=1,ndd
       do j=1,ndd
          s=s1(i)%sigmaf(j)*x(i)*x(j)+s
       enddo
    enddo
    s2=s
    call kill(id)
    call kill(s)
    call kill(x,6)
  END SUBROUTINE tenv_8

  ! End of Some polymorphism

end module S_extend_poly



