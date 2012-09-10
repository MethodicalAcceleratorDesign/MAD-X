!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tree_element_MODULE
  USE polymorphic_complextaylor, mon1=>mon
  IMPLICIT NONE
  public
  integer,private,parameter::ndd=6

  PRIVATE track_TREE,track_TREEP,KILL_TREE,KILL_TREE_N,SET_TREE
  PRIVATE track_TREE_G,track_TREEP_g
  PRIVATE ALLOC_SPINOR_8,ALLOC_probe_8
  PRIVATE KILL_SPINOR_8,KILL_probe_8,KILL_DASPIN
  PRIVATE EQUAL_SPINOR8_SPINOR8,EQUAL_IDENTITY_SPINOR_8 !,EQUAL_SPINOR8_RAY8,EQUAL_RAY8_SPINOR8,
  PRIVATE  EQUAL_IDENTITY_probe_8,ALLOC_daspin,EQUAL_DAMAPSPIN_RAY8 ,EQUAL_RAY8_DAMAPSPIN
  private alloc_normal_spin,kill_normal_spin,EQUAL_IDENTITY_probe
  private READ_DASPIN,PRINT_DASPIN,EQUAL_IDENTITY_SPINOR,EQUAL_PROBE_REAL6
  PRIVATE EQUAL_SPINOR8_SPINOR,EQUAL_PROBE8_PROBE,EQUAL_PROBE8_REAL6
  private EQUAL_IDENTITY_SPINOR_8_r3 ,EQUAL_SPINOR_SPINOR8,exp_spinor_8
  private INTO_RES_SPIN8_eq,INTO_SPIN8_from_RES_eq,ALLOC_rf_phasor_8,KILL_rf_phasor_8

  private find_ar
  private find_ap,PRINT_spinor_8,dot_spinor_8,dot_spinor,print_res_spinor_8
  private  find_n_thetar,find_n_thetap,read_spinor_8
  !  private smatp,smatmulp
  private inv_asr,inv_asp !,inv_as
  PRIVATE EQUAL_DAmapSPIN_int,daddsc,scdadd,EQUAL_PROBE8_PROBE8,PRINT_probe8
  PRIVATE concat,assmap,EQUAL_damapspin,CUTORDER,assprobe_8,POWMAP,cmul,addm,mmul,spin8_mul_map
  private read_probe8,ALLOC_33t,ALLOC_33p,KILL_33t,KILL_33p,zero_33t
  private purge_transverse
  private get_spin_nx_r,get_spin_nx_t,get_spin_nx_rd,get_spin_nx_probe,get_spin_nx_spinor_8
  private scdaddo,daddsco,damapspin_spinor8_mul,damapspin_spinor_mul,eval_spinor_8
  private real_8REAL6,REAL6real_8,real_8REAL_8,PRINT6
  private check_fix,test_jc,A_OPT_damap,K_OPT_damap,factor_am,factor_as,concatxp
  private find_axisp,spin8_scal8_map,add_spin8_spin8,sub_spin8_spin8,mul_spin8_spin8
  private find_perp_basisp,find_exponentp
  private alloc_res_SPINOR_8,KILL_res_SPINOR_8,full_abst,EQUAL_RF8_RF8,extract_envelope_probe8
  PRIVATE EQUAL_RF8_RF,EQUAL_RF_RF8,print_rf_phasor_8,concat_envelope,extract_envelope_damap
  private EQUAL_DAMAP_RAY8,spimat_spinmat,EQUAL_damapspin_smat,EQUAL_smat_damapspin
  private flip,dlie,dphase,phase_shift  ! flip in lielib
  integer, target :: spin_extra_tpsa = 0 ,n0_normal= 2
  integer :: clockwise=1
  logical(lp) :: force_positive=.false.
  logical(lp) :: use_ptc_ac_position=.false.
  integer, private :: nd_used,i_phase,i_plane, decal
  logical :: onelie = my_false
  logical :: firstfac=.true.
  integer, private, parameter :: nfac=20
  real(dp), private :: fac(0:nfac)

  INTERFACE assignment (=)
     !
     MODULE PROCEDURE REAL_8REAL6
     MODULE PROCEDURE REAL6REAL_8
     MODULE PROCEDURE real_8REAL_8
     module PROCEDURE  spimat_spinmat
     !
     MODULE PROCEDURE EQUAL_IDENTITY_probe   ! probe_8=0 p%s=Id  p%x=x_i (tpsa) i=1,npara_fpp
     MODULE PROCEDURE EQUAL_IDENTITY_probe_8 ! probe=0 p%s=Id  p%x=0.d0
     MODULE PROCEDURE EQUAL_IDENTITY_SPINOR_8 ! SPINOR_8%x(r)=1 if SPINOR_8=r otherwise zero
     MODULE PROCEDURE EQUAL_IDENTITY_SPINOR  ! spinor%x(r)=1 if spinor=r otherwise zero
     MODULE PROCEDURE EQUAL_PROBE_REAL6
     MODULE PROCEDURE EQUAL_PROBE8_REAL6
     MODULE PROCEDURE EQUAL_PROBE8_PROBE

     MODULE PROCEDURE EQUAL_SPINOR8_SPINOR8
     MODULE PROCEDURE EQUAL_SPINOR8_SPINOR
     MODULE PROCEDURE EQUAL_DAMAPSPIN_RAY8
     MODULE PROCEDURE EQUAL_RAY8_DAMAPSPIN
     MODULE PROCEDURE EQUAL_DAMAP_RAY8

     MODULE PROCEDURE EQUAL_DAmapSPIN_int  ! damapspin = Identity (or zero)
     MODULE PROCEDURE EQUAL_PROBE8_PROBE8
     MODULE PROCEDURE      EQUAL_damapspin ! damapspin <= real spin matrix
     MODULE PROCEDURE EQUAL_damapspin_smat ! damapspin => real spin matrix
     MODULE PROCEDURE EQUAL_smat_damapspin
     MODULE PROCEDURE     EQUAL_SPINOR_SPINOR8
     MODULE PROCEDURE     EQUAL_PROBE_PROBE8
     MODULE PROCEDURE     normalise_spin
     MODULE PROCEDURE     INTO_RES_SPIN8_eq !
     MODULE PROCEDURE     INTO_SPIN8_from_RES_eq
     MODULE PROCEDURE     EQUAL_RF8_RF8
     MODULE PROCEDURE     EQUAL_RF8_RF
     MODULE PROCEDURE     EQUAL_RF_RF8
     !    MODULE PROCEDURE     EQUAL_nn_damap
  end  INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE concat  ! damapspin= damapspin o damapspin
     MODULE PROCEDURE cmul    ! damapspin%s= real(dp) * damapspin%s     orbital unchanged
     MODULE PROCEDURE mmul  ! damapspin%s=damapspin%s o damap
     MODULE PROCEDURE damapspin_spinor_mul  ! spinor_8%s(i)= damapspin%s(i,j) * spinor%s(j)
     MODULE PROCEDURE damapspin_spinor8_mul ! spinor_8%s(i)= damapspin%s(i,j) * spinor_8%s(j)
     MODULE PROCEDURE spin8_mul_map   ! spinor_8=spinor_8 o damap
     MODULE PROCEDURE eval_spinor_8   ! spinor=spinor_8 * xp(lnv)
     MODULE PROCEDURE concatxp   !  damapspin=damapspin*xp(lnv)
     MODULE PROCEDURE spin8_scal8_map  ! real_8*spinor_8
     MODULE PROCEDURE mul_spin8_spin8  !  mul_spin8_spin8= spin8 X spin8  cross product
  END  INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POWMAP
  END  INTERFACE

  INTERFACE OPERATOR (.cut.)
     MODULE PROCEDURE CUTORDER
  END  INTERFACE

  INTERFACE OPERATOR (.dot.)
     MODULE PROCEDURE dot_spinor
     MODULE PROCEDURE dot_spinor_8
  END  INTERFACE

  INTERFACE operator (+)
     MODULE PROCEDURE scdadd
     MODULE PROCEDURE daddsc
     MODULE PROCEDURE scdaddo
     MODULE PROCEDURE daddsco
     MODULE PROCEDURE addm    !    (damapspin1%m, damapspin1%s  + damapspin2%s)
     MODULE PROCEDURE add_spin8_spin8   ! spinor_8 = spinor_8 + spinor_8
  END  INTERFACE

  INTERFACE operator (-)
     MODULE PROCEDURE sub_spin8_spin8
  END  INTERFACE

  INTERFACE extract_beam_sizes   !
     MODULE PROCEDURE extract_envelope_damap  !
     MODULE PROCEDURE extract_envelope_probe8  !
  END INTERFACE

  INTERFACE exp   !  damapspin = exp(spinor_8)
     MODULE PROCEDURE exp_spinor_8  !
  END INTERFACE

  INTERFACE texp   !  damapspin = exp(spinor_8)
     MODULE PROCEDURE exp_spinor_8  !
  END INTERFACE

  INTERFACE find_n0   ! (s0(3,3),n0(3))
     MODULE PROCEDURE find_n_thetar  ! finds n0 the naive way
     MODULE PROCEDURE find_n_thetap !(s0(3,3),n0(3), spinor_8 optional)
  END INTERFACE

  INTERFACE find_axis   ! (ds,spinor_8)
     MODULE PROCEDURE find_axisp    ! s= exp(spinor_8.dot.L)
  END INTERFACE

  INTERFACE find_perp_basis   !
     MODULE PROCEDURE find_perp_basisp    !
  END INTERFACE

  INTERFACE find_exponent   !
     MODULE PROCEDURE find_exponentp    !
  END INTERFACE

  INTERFACE find_a    ! (n(3),a(3,3))  if you have n you get a(3,3)
     MODULE PROCEDURE find_ar  ! such that ai.s.a = exp(theta L_y)
     MODULE PROCEDURE find_ap
  END INTERFACE

  INTERFACE inv_as    ! gets inverse of spin matrix by transposing
     MODULE PROCEDURE inv_asr
     MODULE PROCEDURE inv_asp
  END INTERFACE

  INTERFACE get_spin_n0  ! (S,theta0,n0)
     MODULE PROCEDURE get_spin_nx_r  ! real s(3,3) = exp(theta0 n0.L)
     MODULE PROCEDURE get_spin_nx_rd  ! takes DS (damapspin but returns real)
     MODULE PROCEDURE get_spin_nx_t   ! takes DS (damapspin but returns real_8)
     MODULE PROCEDURE get_spin_nx_probe ! takes probe (damapspin but returns real)
     MODULE PROCEDURE get_spin_nx_spinor_8
  END INTERFACE

  INTERFACE find_exponent_jet  !
     MODULE PROCEDURE find_exponent_jet_p
  end interface

  INTERFACE full_abs  !
     MODULE PROCEDURE full_abst
  END INTERFACE   !



  INTERFACE PRINT
     MODULE PROCEDURE PRINT6
!!!!
     MODULE PROCEDURE PRINT_DASPIN
     MODULE PROCEDURE PRINT_probe8
     MODULE PROCEDURE PRINT_spinor_8
     MODULE PROCEDURE print_res_spinor_8
     MODULE PROCEDURE print_rf_phasor_8
  END INTERFACE

  INTERFACE daPRINT
     MODULE PROCEDURE PRINT6
!!!!
  END INTERFACE

  INTERFACE READ
     MODULE PROCEDURE READ_DASPIN
     MODULE PROCEDURE read_probe8   ! a bit illegal : reading polymorphs as taylor...
     MODULE PROCEDURE read_spinor_8
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_SPINOR_8
     MODULE PROCEDURE ALLOC_probe_8
     MODULE PROCEDURE ALLOC_daspin
     MODULE PROCEDURE A_OPT_damap
     MODULE PROCEDURE alloc_normal_spin
     MODULE PROCEDURE alloc_res_SPINOR_8
     MODULE PROCEDURE ALLOC_rf_phasor_8
     MODULE PROCEDURE SET_TREE
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILL_SPINOR_8
     MODULE PROCEDURE KILL_probe_8
     MODULE PROCEDURE KILL_DASPIN
     MODULE PROCEDURE K_OPT_damap
     MODULE PROCEDURE KILL_normal_spin
     MODULE PROCEDURE KILL_res_SPINOR_8
     MODULE PROCEDURE KILL_rf_phasor_8
  END INTERFACE

  INTERFACE ALLOC_33
     MODULE PROCEDURE ALLOC_33t
     MODULE PROCEDURE ALLOC_33p
  END INTERFACE

  INTERFACE ALLOC_nn
     MODULE PROCEDURE ALLOC_33t
     MODULE PROCEDURE ALLOC_33p
  END INTERFACE

  INTERFACE matmul_nn
     MODULE PROCEDURE matmul_33
  END INTERFACE

  INTERFACE zero_33
     MODULE PROCEDURE zero_33t
     MODULE PROCEDURE zero_33p
  END INTERFACE

  INTERFACE zero_nn
     MODULE PROCEDURE zero_33t
     MODULE PROCEDURE zero_33p
  END INTERFACE

  INTERFACE KILL_nn
     MODULE PROCEDURE KILL_33t
     MODULE PROCEDURE KILL_33p
  END INTERFACE

  INTERFACE KILL_33
     MODULE PROCEDURE KILL_33t
     MODULE PROCEDURE KILL_33p
  END INTERFACE

  INTERFACE track
     MODULE PROCEDURE track_TREE
     MODULE PROCEDURE track_TREEP
  END INTERFACE


  INTERFACE trackg
     MODULE PROCEDURE track_TREE_G
     MODULE PROCEDURE track_TREEP_g
  END INTERFACE


  INTERFACE KILL
     MODULE PROCEDURE KILL_TREE
     MODULE PROCEDURE KILL_TREE_N
     MODULE PROCEDURE KILL_33t
     MODULE PROCEDURE KILL_33p
  END INTERFACE

  INTERFACE factor
     MODULE PROCEDURE factor_am
     MODULE PROCEDURE factor_as
  END INTERFACE

  INTERFACE ass
     MODULE PROCEDURE assmap
     MODULE PROCEDURE assprobe_8
  END INTERFACE


CONTAINS
  SUBROUTINE  A_OPT_damap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (damapspin),INTENT(INout)::S1,S2
    type (damapspin),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call alloc(s1)
    call alloc(s2)
    if(present(s3)) call alloc(s3)
    if(present(s4)) call alloc(s4)
    if(present(s5)) call alloc(s5)
    if(present(s6)) call alloc(s6)
    if(present(s7)) call alloc(s7)
    if(present(s8)) call alloc(s8)
    if(present(s9)) call alloc(s9)
    if(present(s10))call alloc(s10)
  END SUBROUTINE A_opt_damap

  SUBROUTINE  K_OPT_damap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (damapspin),INTENT(INout)::S1,S2
    type (damapspin),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call KILL(s1)
    call KILL(s2)
    if(present(s3)) call KILL(s3)
    if(present(s4)) call KILL(s4)
    if(present(s5)) call KILL(s5)
    if(present(s6)) call KILL(s6)
    if(present(s7)) call KILL(s7)
    if(present(s8)) call KILL(s8)
    if(present(s9)) call KILL(s9)
    if(present(s10))call KILL(s10)
  END SUBROUTINE K_OPT_damap

!!! use to be in extend_poly

  FUNCTION daddsco( S1, S2 )
    implicit none
    TYPE (real_8) daddsco(ndd)
    TYPE (damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    integer localmaster,iia(4),ico(4),nd2,i
    call liepeek(iia,ico)
    nd2=iia(4)


    do i=1,nd2
       localmaster=master
       call ass(daddsco(i))
       daddsco(i)=s1%v(i)+s2(i)-(s1%v(i).sub.'0')
       master=localmaster
    enddo
    do i=nd2+1,ndd
       localmaster=master
       call ass(daddsco(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5) then
          daddsco(i)=s2(i)+(1.0_dp.mono.'00001')
       else
          daddsco(i)=s2(i)
       endif
       master=localmaster
    enddo

  END FUNCTION daddsco

  FUNCTION scdaddo( S2,S1  )
    implicit none
    TYPE (real_8) scdaddo(ndd)
    TYPE (damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    integer localmaster,iia(4),ico(4),nd2,i
    call liepeek(iia,ico)
    nd2=iia(4)

    do i=1,nd2
       localmaster=master
       call ass(scdaddo(i))
       scdaddo(i)=s1%v(i)+s2(i)-(s1%v(i).sub.'0')
       master=localmaster
    enddo
    do i=nd2+1,ndd
       localmaster=master
       call ass(scdaddo(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5) then
          scdaddo(i)=s2(i)+(1.0_dp.mono.'00001')
       else
          scdaddo(i)=s2(i)
       endif
       master=localmaster
    enddo


  END FUNCTION scdaddo

  SUBROUTINE  REAL6real_8(S2,S1)
    implicit none
    real(dp),INTENT(inOUT)::S2(ndd)
    type (real_8),INTENT(IN)::S1(ndd)
    integer i


    do i=1,ndd
       s2(i)=s1(i)          !%t
    enddo
  END SUBROUTINE REAL6real_8

  SUBROUTINE  real_8REAL_8(S1,S2)
    implicit none
    type (real_8),INTENT(in)::S2(ndd)
    type (real_8),INTENT(inOUT)::S1(ndd)
    integer i


    do i=1,ndd
       s1(i)=s2(i)
    enddo
  END SUBROUTINE real_8REAL_8

  SUBROUTINE  real_8REAL6(S1,S2)
    implicit none
    real(dp),INTENT(in)::S2(ndd)
    type (real_8),INTENT(inOUT)::S1(ndd)
    integer i


    do i=1,ndd
       s1(i)=s2(i)
    enddo
  END SUBROUTINE real_8REAL6



  SUBROUTINE  spimat_spinmat(S1,S2)
    implicit none
    type (spinmatrix),INTENT(in)::S2
    type (spinmatrix),INTENT(inOUT)::S1
    integer i,j

    do i=1,3
       do j=1,3
          s1%s(i,j)=s2%s(i,j)
       enddo
    enddo
  END SUBROUTINE spimat_spinmat

  SUBROUTINE  print6(S1,mf)
    implicit none
    type (real_8),INTENT(INout)::S1(ndd)
    integer        mf,i

    do i=1,ndd
       call print(s1(i),mf)
    enddo

  END SUBROUTINE print6

!!! end of "use to be in extend_poly"

  SUBROUTINE COPY_TREE(T,U)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: U

    IF(.NOT.ASSOCIATED(T%CC).OR..NOT.ASSOCIATED(U%CC) ) RETURN
    IF(SIZE(T%CC)/=SIZE(U%CC) ) STOP 888

    U%CC=T%CC
    U%JL=T%JL
    U%JV=T%JV
    U%N=T%N
    U%ND2=T%ND2
    U%no=T%no

  END SUBROUTINE COPY_TREE

  SUBROUTINE COPY_TREE_N(T,U)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T(:)
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: U(:)
    INTEGER I

    DO I=1,SIZE(T)
       CALL COPY_TREE(T(I),U(I))
    ENDDO

  END SUBROUTINE COPY_TREE_N



  SUBROUTINE NULL_TREE(T)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T

    NULLIFY(T%CC,T%JL,T%JV,T%N,T%ND2,T%no)

  END SUBROUTINE NULL_TREE


  SUBROUTINE ALLOC_TREE(T,N,ND2)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T
    INTEGER , INTENT(IN) :: N,ND2

    !IF(N==0) RETURN

    ALLOCATE(T%CC(N),T%fix(nd2),T%JL(N),T%JV(N),T%N,T%ND2,T%no)
    T%N=N
    T%ND2=ND2
    T%no=0
    T%fix=0.0_dp
  END SUBROUTINE ALLOC_TREE

  SUBROUTINE SET_TREE(T,MA)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T
    TYPE(DAMAP), INTENT(INOUT) :: MA
    INTEGER N
    TYPE(DAMAP) M

    CALL ALLOC(M)

    CALL  mtree(MA%v%I,C_%ND2,M%v%I,C_%ND2)
    CALL ppushGETN(M%v%I,C_%ND2,N)


    CALL ALLOC_TREE(T,N,C_%ND2)
    T%no=c_%no

    CALL ppushstore(M%v%I,C_%ND2,T%CC,T%JL,T%JV)

    CALL KILL(M)

  END SUBROUTINE SET_TREE

  ! FOR FAST B FIELD IN PACKAGE OF PTC
  SUBROUTINE SET_TREE_G(T,MA)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T
    TYPE(taylor), INTENT(INOUT) :: MA(:)
    INTEGER N,NP
    TYPE(taylor), ALLOCATABLE :: M(:)

    NP=SIZE(MA)

    ALLOCATE(M(NP))
    CALL ALLOC(M,NP)

    CALL  mtree(MA%I,NP,M%I,NP)
    CALL ppushGETN(M%I,NP,N)


    CALL ALLOC_TREE(T,N,NP)
    T%no=c_%no

    CALL ppushstore(M%I,NP,T%CC,T%JL,T%JV)

    CALL KILL(M,NP)
    deallocate(M)
  END SUBROUTINE SET_TREE_G

  SUBROUTINE track_TREE_G(T,XI)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    REAL(DP), INTENT(INOUT) :: XI(:)
    REAL(DP) XT(lno),XF(lnv),XM(lno+1),XX
    INTEGER JC,I,IV

    XT=0.0_dp
    XF=0.0_dp
    XM=0.0_dp

    do i=1,T%ND2
       xt(i)=xi(i)
    enddo
    do i=1,T%ND2
       xf(i) = T%cc(i)
    enddo

    XM(1) = 1.0_dp
    JC=T%ND2
    do i=1,(T%N-T%ND2)/T%ND2
       !
       xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%ND2
          jc=jc+1
          xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo
    xi=xf


  END SUBROUTINE track_TREE_G



  SUBROUTINE track_TREEP_g(T,XI)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    TYPE(REAL_8), INTENT(INOUT) :: XI(:)
    TYPE(REAL_8) XT(lno),XF(lnv),XM(lno+1),XX
    INTEGER JC,I,IV

    CALL ALLOC(XT,lno)
    CALL ALLOC(XF,lnv)
    CALL ALLOC(XM,lno+1)
    CALL ALLOC(XX)




    do i=1,T%ND2
       xt(i)=xi(i)
    enddo
    do i=1,T%ND2
       xf(i) = T%cc(i)
    enddo

    XM(1) = 1.0_dp
    JC=T%ND2

    do i=1,(T%N-T%ND2)/T%ND2
       !
       xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%ND2
          jc=jc+1
          xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo

    do i=1,T%ND2
       xI(i)=xF(i)
    enddo

    CALL KILL(XT,lno)
    CALL KILL(XF,lnv)
    CALL KILL(XM,lno+1)
    CALL KILL(XX)

  END SUBROUTINE track_TREEP_g





  ! END OF FAST B FIELD IN PACKAGE OF PTC

  !  type  tree_element
  !     real(dp) ,  DIMENSION(:), POINTER :: CC
  !     integer,  DIMENSION(:), POINTER :: JL,JV
  !     INTEGER,POINTER :: N,ND2
  !  end  type tree_element


  SUBROUTINE KILL_TREE(T)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T


    IF(ASSOCIATED(T%CC))   DEALLOCATE(T%CC,T%fix,T%JL,T%JV,T%N,T%ND2,T%No)


  END SUBROUTINE KILL_TREE

  SUBROUTINE KILL_TREE_N(T)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T(:)
    INTEGER I

    DO I=1,SIZE(T)
       CALL KILL(T(I))
    ENDDO

  END SUBROUTINE KILL_TREE_N




  SUBROUTINE track_TREE(T,XI,n)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    REAL(DP), INTENT(INOUT) :: XI(6)
    integer, optional, INTENT(IN) :: n
    integer n1,k
    REAL(DP) XT(lno),XF(6),XM(lno+1),XX
    INTEGER JC,I,IV
    xi(1:t%nd2)=xi(1:t%nd2)-t%fix
    n1=1
    if(present(n)) n1=n
    do k=1,n1
       if(.not.c_%CHECK_STABLE) return
       XT=0.0_dp
       XF=0.0_dp
       XM=0.0_dp

       do i=1,T%ND2
          xt(i)=xi(i)
       enddo
       do i=1,T%ND2
          xf(i) = T%cc(i)
       enddo

       XM(1) = 1.0_dp
       JC=T%ND2
       do i=1,(T%N-T%ND2)/T%ND2
          !
          xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
          xm(T%jl(JC+1)+1) = xx
          !
          do iv=1,T%ND2
             jc=jc+1
             xf(iv) = xf(iv) + t%cc(jc) * xx
          enddo
       enddo
       xi=xf

       if(abs(xi(1))>c_%absolute_aperture.or.abs(xi(3))>c_%absolute_aperture) then
          c_%CHECK_STABLE=.FALSE.
          xlost=xi
       endif
    enddo

  END SUBROUTINE track_TREE

  SUBROUTINE track_TREEP(T,XI,n)
    use da_arrays
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    TYPE(REAL_8), INTENT(INOUT) :: XI(6)
    integer, optional, INTENT(IN) :: n
    integer n1,k
    TYPE(REAL_8) XT(lno),XF(6),XM(lno+1),XX
    INTEGER JC,I,IV

    n1=1
    if(present(n)) n1=n
    do k=1,n1

       CALL ALLOC(XT,lno)
       CALL ALLOC(XF,6)
       CALL ALLOC(XM,lno+1)
       CALL ALLOC(XX)




       do i=1,T%ND2
          xt(i)=xi(i)
       enddo
       do i=1,T%ND2
          xf(i) = T%cc(i)
       enddo

       XM(1) = 1.0_dp
       JC=T%ND2

       do i=1,(T%N-T%ND2)/T%ND2
          !
          xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
          xm(T%jl(JC+1)+1) = xx
          !
          do iv=1,T%ND2
             jc=jc+1
             xf(iv) = xf(iv) + t%cc(jc) * xx
          enddo
       enddo

       do i=1,T%ND2
          xI(i)=xF(i)
       enddo

       CALL KILL(XT,lno)
       CALL KILL(XF,6)
       CALL KILL(XM,lno+1)
       CALL KILL(XX)

    enddo
  END SUBROUTINE track_TREEP

  !  READING COSY MAPS AND UNIVERSAL TAYLORS
  SUBROUTINE  dainput_SPECIAL6(S1,MFILE,COSY)
    implicit none
    INTEGER,INTENT(in)::MFILE,COSY
    type (damap),INTENT(INOUT)::S1
    INTEGER I,js(6),k1,k2,l,ncoef
    real(dp) x
    character*200 line
    S1=0
    SELECT CASE(COSY)
    CASE(0:1)
       call dainput(s1,mfile)
    CASE(2)  ! COSY INFINITY
       do i=1,c_%nd2
          read(mfile,'(a200)') line
          read(mfile,'(a200)') line
          do while(line(14:16)/='---')
             read(line,*)  k1,x,k2,js
             s1%v(i)=s1%v(i)+(x.mono.js)
             read(mfile,'(a200)') line
          enddo
       enddo

    CASE(-1)  ! sagan
       do i=1,c_%nd2
          read(mfile,*) ncoef
          do l=1,ncoef
             read(mfile,*)  x,js
             s1%v(i)=s1%v(i)+(x.mono.js)
          enddo
       enddo

    CASE DEFAULT

       WRITE(6,*) " NOT SUPPORTED IN DAREADMAP_SPECIAL"
       STOP 111

    END SELECT


  END SUBROUTINE dainput_SPECIAL6

  ! SYMPLECTIFY A MAP NEAR THE IDENTITY
  SUBROUTINE symplectic(m,eps,nst)
    IMPLICIT NONE
    type(damap), INTENT(INOUT) :: m
    type(onelieexponent) uno
    type(damap) id
    real(dp), optional :: eps
    integer nst

    call alloc(uno); call alloc(id);

    if(present(eps))then
       if(eps>0.0_dp) uno%eps=eps
    endif
    uno=m
    id=1
    uno%pb%h=uno%pb%h/nst
    m=texp(uno%pb,id)

    call kill(uno); call kill(id);


  end SUBROUTINE symplectic

  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  SPIN STUFF IS HERE




  subroutine assmap(s1)
    implicit none
    TYPE (damapspin) s1
    integer i,j

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)
       w_p=0
       w_p%nc=1
       w_p=(/" cannot indent anymore "/)
       w_p%fc='(1((1X,A72),/))'
       ! call !write_e(100)
    end select

    do i=1,3
       do j=1,3
          call assp_no_master(s1%s%s(i,j))
          !          call ass0(s1%s(i,j)%t)   ! taylor
          !          s1%s(i,j)%alloc=my_true
          !          s1%s(i,j)%kind=1
          !          s1%s(i,j)%i=0
       enddo
    enddo
    do i=1,c_%nd2
       call ass0(s1%m%v(i))
    enddo

  end subroutine assmap

  subroutine assprobe_8(s1)
    implicit none
    TYPE (probe_8) s1
    integer i,j

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)
       w_p=0
       w_p%nc=1
       w_p=(/" cannot indent anymore "/)
       w_p%fc='(1((1X,A72),/))'
       ! call !write_e(100)
    end select

    !    if(C_%SPIN_POS/=0) then
    !       do i=1,3
    !          call ass0(s1%s%x(i)%t)
    !          s1%s%x(i)%alloc=my_true
    !          s1%s%x(i)%kind=2
    !          s1%s%x(i)%i=0
    !       enddo
    do j=1,3
       do i=1,3
          call assp_no_master(s1%s(j)%x(i))
          !          call ass0(s1%s(j)%x(i)%t)
          !          s1%s(j)%x(i)%alloc=my_true
          !          s1%s(j)%x(i)%kind=1
          !          s1%s(j)%x(i)%i=0
       enddo
    enddo

    !    endif
    do i=1,6  !
       call assp_no_master(s1%x(i))
       !       call ass0(s1%x(i)%t)
       !       s1%x(i)%alloc=my_true
       !       s1%x(i)%kind=1
       !       s1%x(i)%i=0
    enddo

    do i=1,2  !
       call assp_no_master(s1%AC%x(i))
    enddo


  end subroutine assprobe_8




  FUNCTION FULL_ABST( S1 )
    implicit none
    REAL(DP) FULL_ABST
    TYPE (damapspin), INTENT (IN) :: S1
    integer i,j

    FULL_ABST=0.0_dp

    do i=1,3
       do j=1,3
          FULL_ABST=FULL_ABS(S1%s%s(i,j))+FULL_ABST
       enddo
    enddo


  END FUNCTION FULL_ABST

  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (damapspin) CUTORDER
    TYPE (damapspin), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer localmaster
    integer i,j

    localmaster=master
    call ass(CUTORDER)

    CUTORDER=0

    do i=1,3
       do j=1,3
          CUTORDER%s%s(i,j)=s1%s%s(i,j).cut.(s2-1)
       enddo
    enddo

    CUTORDER%m=s1%m.cut.s2

    master=localmaster

  END FUNCTION CUTORDER

  subroutine inv_damapspin(r,s)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: S
    TYPE(damapspin), INTENT(IN) :: r
    TYPE(damapspin) t
    INTEGER I,j

    call alloc(t)

    t=r
    s%m=t%m**(-1)

    do i=1,3
       do j=1,3
          s%s%s(i,j)=t%s%s(j,i)
          !          if(s%s(i,j)%kind==2) then   ! taylor
          s%s%s(i,j)=s%s%s(i,j)*s%m
          !           s%s(i,j)=s%s(i,j)%t*s%m ! taylor
          !          endif   ! taylor
       enddo
    enddo

    call kill(t)

  END    subroutine inv_damapspin


  !  type damapspin
  !   REAL(DP) X(6)
  !   type(damap) M
  !   type(real_8) s(3,3)
  !  end type damapspin




  subroutine EQUAL_damapspin(S,R)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: S
    TYPE(damapspin), INTENT(IN) :: r
    INTEGER I,j

    do i=1,3
       do j=1,3
          s%s%s(i,j)=r%s%s(i,j)
       enddo
    enddo

    s%m=r%m
    s%e_ij=r%e_ij
    !    s%s0=r%s0

    !     do i=1,6
    !     s%x(i)=r%x(i)
    !     enddo

  END    subroutine EQUAL_damapspin


  subroutine EQUAL_damapspin_smat(S,R)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: S
    real(dp), INTENT(IN) :: r(3,3)
    INTEGER I,j

    do i=1,3
       do j=1,3
          s%s%s(i,j)=r(i,j)
       enddo
    enddo


  END    subroutine EQUAL_damapspin_smat

  subroutine EQUAL_smat_damapspin(R,S)
    implicit none
    real(dp), INTENT(INout) :: r(3,3)
    TYPE(damapspin), INTENT(IN) :: S
    INTEGER I,j

    do i=1,3
       do j=1,3
          r(i,j)=s%s%s(i,j)
       enddo
    enddo


  END    subroutine EQUAL_smat_damapspin


  subroutine EQUAL_IDENTITY_SPINOR_8(S,R)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    INTEGER, INTENT(IN) :: R
    INTEGER I

    !     S%G=A_PARTICLE
    IF(R==1.or.r==2.or.r==3) THEN
       DO I=1,3    !C_%NSPIN
          S%X(I)=0.0_dp
       enddo
       S%X(r)=1.0_dp
       !        write(6,*) " in EQUAL_IDENTITY_SPINOR_8"
       !        stop 888
       !       DO I=1,C_%NSPIN
       !          S%X(I)=ONE.MONO.(c_%npara_fpp-C_%NSPIN+I)
       !       ENDDO
    ELSEIF(R==0) THEN
       DO I=1,3    !C_%NSPIN
          S%X(I)=0.0_dp
       enddo
    ELSE
       write(6,*) " stopped in EQUAL_IDENTITY_SPINOR_8"
       STOP 100
    ENDIF

  END    subroutine EQUAL_IDENTITY_SPINOR_8


  subroutine EQUAL_IDENTITY_SPINOR_8_r3(S,R)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    real(dp), INTENT(IN) :: R(3)
    INTEGER I

    !     S%G=A_PARTICLE

    DO I=1,3    !C_%NSPIN
       S%X(I)=r(i)   !+(ONE.MONO.(c_%npara_fpp-C_%NSPIN+I))
    ENDDO

  END    subroutine EQUAL_IDENTITY_SPINOR_8_r3

  subroutine EQUAL_IDENTITY_SPINOR (S,R)
    implicit none
    TYPE(SPINOR), INTENT(INOUT) :: S
    INTEGER, INTENT(IN) :: R
    INTEGER I

    !     S%G=A_PARTICLE
    IF(R==1.or.r==2.or.r==3) THEN

       DO I=1,3    !C_%NSPIN
          S%X(I)=0.0_dp
       enddo
       S%X(r)=1.0_dp
    ELSEIF(R==0) THEN
       DO I=1,3   !C_%NSPIN
          S%X(I)=0.0_dp
       enddo
    ELSE
       STOP 100
    ENDIF

  END    subroutine EQUAL_IDENTITY_SPINOR

  subroutine EQUAL_PROBE_REAL6 (P,X)
    implicit none
    TYPE(PROBE), INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: X(6)
    INTEGER I
    P%u=my_false
    !       P%s(0)%x=zero
    !       P%s(0)%x(n0_normal)=one

    DO I=1,ISPIN1R
       P%s(i)%x=0.0_dp
       P%s(i)%x(i)=1.0_dp
    enddo
    P%X=X
    P%ac%t=0.0_dp

  END    subroutine EQUAL_PROBE_REAL6

  subroutine EQUAL_PROBE8_REAL6 (P,X)
    implicit none
    TYPE(PROBE_8), INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: X(6)
    INTEGER I

    P%u=my_false
    !    P%S=0

    !       P%s(0)=0
    !       P%s(0)%x(n0_normal)=one

    DO I=1,3
       P%s(i)=0
       P%s(i)%x(i)=1.0_dp
    enddo

    DO I=1,6
       P%X(i)=X(i)
    enddo

    P%AC%X(1)=0.0_dp
    P%AC%X(2)=0.0_dp
    P%AC%t=0.0_dp
    p%e_ij=0.0_dp
  END    subroutine EQUAL_PROBE8_REAL6

  subroutine EQUAL_PROBE8_PROBE8(P8,P)
    implicit none
    TYPE(PROBE_8), INTENT(IN) :: P
    TYPE(PROBE_8), INTENT(INOUT) :: P8
    INTEGER I,J

    DO I=1,6
       P8%X(I)=P%X(I)
    ENDDO

    DO I=1,6
       DO j=1,6
          P8%E_ij(I,j)=P%E_ij(I,j)
       ENDDO
    ENDDO

    !    P8%S=P%S
    DO I=1,3
       P8%S(i)=P%S(i)
    enddo
    P8%AC=P%AC

    P8%u=P%u

  END subroutine EQUAL_PROBE8_PROBE8



  subroutine EQUAL_RF8_RF8(P8,P)
    implicit none
    TYPE(rf_phasor_8), INTENT(IN) :: P
    TYPE(rf_phasor_8), INTENT(INOUT) :: P8
    INTEGER I

    DO I=1,2
       P8%X(I)=P%X(I)
    ENDDO
    P8%om=P%om
    P8%t=P%t

  END subroutine EQUAL_RF8_RF8

  subroutine EQUAL_RF8_RF(P8,P)
    implicit none
    TYPE(rf_phasor_8), INTENT(INOUT) :: P8
    TYPE(rf_phasor), INTENT(IN) :: P
    INTEGER I

    DO I=1,2
       P8%X(I)=P%X(I)
    ENDDO
    P8%om=P%om
    P8%t=P%t

  END subroutine EQUAL_RF8_RF

  subroutine EQUAL_RF_RF8(P,P8)
    implicit none
    TYPE(rf_phasor), INTENT(INOUT) :: P
    TYPE(rf_phasor_8), INTENT(IN) :: P8
    INTEGER I

    DO I=1,2
       P%X(I)=P8%X(I)
    ENDDO
    P%om=P8%om
    P%t=P8%t

  END subroutine EQUAL_RF_RF8



  subroutine EQUAL_PROBE8_PROBE (P8,P)
    implicit none
    TYPE(PROBE), INTENT(IN) :: P
    TYPE(PROBE_8), INTENT(INOUT) :: P8
    INTEGER I
    DO I=1,6
       P8%X(I)=P%X(I)
    ENDDO
    !    P8%S=P%S
    !    do i=1,3
    !     P8%S(i)=0
    !     P8%S(i)%x(i)=one
    !    enddo
    !     P8%S(0)=P%s(0)
    do I=1,ISPIN1R
       P8%S(I)=P%s(I)
    enddo


    P8%u=P%u
    P8%e_ij=0.0_dp
  END subroutine EQUAL_PROBE8_PROBE

  subroutine EQUAL_PROBE_PROBE8 (P,P8)
    implicit none
    TYPE(PROBE), INTENT(INOUT) :: P
    TYPE(PROBE_8), INTENT(IN) :: P8
    INTEGER I
    DO I=1,6
       P%X(I)=P8%X(I)
    ENDDO
    !     P%S(0)=P8%S(0)
    DO I=1,ISPIN1R
       P%S(I)=P8%S(I)
    ENDDO
    P%u=P8%u

  END subroutine EQUAL_PROBE_PROBE8

  subroutine EQUAL_IDENTITY_probe(R,S)
    implicit none
    TYPE(probe), INTENT(INOUT) :: R
    INTEGER, INTENT(IN) :: S
    INTEGER I

    R%S(1)=0
    R%S(2)=0
    R%S(3)=0
    DO I=1,6
       R%X(I)=0.0_dp
    ENDDO
    IF(S==1) THEN
       !      R%S(0)%x(n0_normal)=one
       R%S(1)=1
       R%S(2)=2
       R%S(3)=3
    ELSEif(s==0) then
    else
       STOP 100
    ENDIF

  END    subroutine EQUAL_IDENTITY_probe

  subroutine EQUAL_IDENTITY_probe_8(R,S)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    INTEGER, INTENT(IN) :: S
    INTEGER I

    !    R%S=S
    !      R%S(0)=0
    R%S(1)=0
    R%S(2)=0
    R%S(3)=0
    DO I=1,6  !-C_%NSPIN
       R%X(I)=0.0_dp
    enddo
    IF(S==1) THEN
       DO I=1,c_%npara_fpp  !-C_%NSPIN
          R%X(I)=1.0_dp.MONO.I
       ENDDO
       !      R%S(0)%x(n0_normal)=one
       R%S(1)=1
       R%S(2)=2
       R%S(3)=3
    ELSEIF(S==0) THEN

       !       DO I=1,6  !-C_%NSPIN
       !        R%X(I)=ZERO!
       !       enddo
    ELSE
       STOP 100
    ENDIF
    R%e_ij=0.0_dp
  END    subroutine EQUAL_IDENTITY_probe_8


  subroutine EQUAL_SPINOR8_SPINOR8(R,S)
    implicit none
    TYPE(SPINOR_8), INTENT(IN) :: S
    TYPE(SPINOR_8), INTENT(INOUT) :: R
    INTEGER I

    DO I=1,3
       R%X(I)=S%X(I)
    ENDDO
    !    R%G=S%G
  END    subroutine EQUAL_SPINOR8_SPINOR8

  subroutine EQUAL_SPINOR_SPINOR8(R,S)
    implicit none
    TYPE(SPINOR_8), INTENT(IN) :: S
    TYPE(SPINOR), INTENT(INOUT) :: R
    INTEGER I

    DO I=1,3
       R%X(I)=S%X(I)
    ENDDO
    !    R%G=S%G
  END    subroutine EQUAL_SPINOR_SPINOR8

  subroutine EQUAL_SPINOR8_SPINOR(R,S)
    implicit none
    TYPE(SPINOR), INTENT(IN) :: S
    TYPE(SPINOR_8), INTENT(INOUT) :: R
    INTEGER I

    DO I=1,3
       R%X(I)=S%X(I)
    ENDDO
    !    R%G=S%G
  END    subroutine EQUAL_SPINOR8_SPINOR

  subroutine EQUAL_DAMAP_RAY8(DS,R)
    implicit none
    TYPE(probe_8), INTENT(IN) :: R
    TYPE(damap), INTENT(INOUT) :: DS
    INTEGER I,nd2t

    nd2t=C_%ND2
    if(doing_ac_modulation_in_ptc) then
       nd2t=C_%ND2-2
    endif



    DO I=1,nd2t
       DS%V(I)=R%X(I)
    ENDDO
    DO I=nd2t+1,C_%ND2
       DS%V(I)=R%ac%x(i-nd2t)
    ENDDO

  END subroutine EQUAL_DAMAP_RAY8


  subroutine EQUAL_DAMAPSPIN_RAY8(DS,R)
    implicit none
    TYPE(probe_8), INTENT(IN) :: R
    TYPE(damapspin), INTENT(INOUT) :: DS
    real(dp) s(3,3)
    INTEGER I,J,nd2t

    nd2t=C_%ND2
    if(doing_ac_modulation_in_ptc) then
       nd2t=C_%ND2-2
    endif



    DO I=1,nd2t
       DS%M%V(I)=R%X(I)
    ENDDO
    DO I=nd2t+1,C_%ND2
       DS%M%V(I)=R%ac%x(i-nd2t)
    ENDDO

    DO I=1,3
       DO J=1,3
          DS%S%s(I,J)=R%S(J)%X(I)
       ENDDO
    ENDDO
    !          DS%S(I,1)=R%SX%X(I)
    !          DS%S(I,2)=R%SY%X(I)
    !          DS%S(I,3)=R%SZ%X(I)
    ds%e_ij=r%e_ij

    !  DO I=1,3
    !    !     n(i)=R%S(0)%X(I)
    !     DO J=1,3
    !        S(I,J)=DS%S(I,J)
    !     ENDDO
    !  ENDDO

    !    call inv_as(s,s)

    !    n= matmul(s,n)

    !    ds%s0=n



    !     DS%G=R%S%G
  END subroutine EQUAL_DAMAPSPIN_RAY8

  subroutine EQUAL_RAY8_DAMAPSPIN(R,DS)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    TYPE(damapspin), INTENT(IN) :: DS
    INTEGER I,J

    !  if(C_%ND2/=6) then
    !       write(6,*) " not implemented in EQUAL_DASPIN_RAY8"
    !       stop
    !      endif




    DO I=1,C_%ND2
       R%X(I)= MORPH(DS%M%V(I))   !+R%X(I)
    ENDDO

    DO J=1,3
       DO I=1,3
          R%S(J)%X(I)=DS%S%s(I,J)                  ! taylor MORPH(DS%S(I,J))
       ENDDO
    ENDDO

    r%e_ij=ds%e_ij



    !    R%S%G=DS%G

  END subroutine EQUAL_RAY8_DAMAPSPIN

  subroutine EQUAL_DAmapSPIN_int(DS,I1)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: DS
    INTEGER, INTENT(IN) :: I1
    INTEGER I

    !       type daspin
    !   REAL(DP) X(6)
    !   type(damap) M
    !   type(real_8) s(3,3)
    !  end type daspin

    !    DS%X=ZERO
    DS%M=I1
    CALL zero_33(DS%S%s)

    DO I=1,3
       if(i1==1) then
          DS%S%s(I,I)=1.0_dp
       endif
    ENDDO
    DS%e_ij=0.0_dp

  END SUBROUTINE EQUAL_DAmapSPIN_int

  FUNCTION scdadd( S2,S1  )
    implicit none
    TYPE (probe_8) scdadd
    TYPE (damapspin), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j,nd2t,dc
    type(taylor) d


    scdadd%u=my_false
    scdadd%E_ij=0.0_dp

    if(doing_ac_modulation_in_ptc) then
       dc=2
    else
       dc=0
    endif
    if(c_%ndpt/=0) then                       ! delta-l without cavity
       nd2t=6
    elseif(c_%ndpt==0.and.c_%nd2==4+dc) then     !  no delta-l maybe delta as a parameter
       nd2t=4
    elseif(c_%ndpt==0.and.c_%nd2==6+dc) then      !   delta-l with cavity
       nd2t=6
    else
       write(6,*) " internal PTC error in scdadd o_tree_element.f90 "
       stop 666
    endif

    call alloc(d)
    do i=1,nd2t                 !   from 1-4 or 1-6 (if ndpt=0)
       localmaster=master
       call ass(scdadd%x(i))
       !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
       scdadd%x(i)=s1%M%V(i)+s2%x(i)-(s1%M%V(i).sub.'0')
       master=localmaster
    enddo

    do i=nd2t+1,6
       localmaster=master
       call ass(scdadd%x(i))
       if((c_%npara==5+dc).AND.I==5) then   ! npr
          scdadd%x(i)=s2%x(i)+(1.0_dp.mono.c_%npara)
       else
          scdadd%x(i)=s2%x(i)
       endif
       master=localmaster
    enddo

    !    if(doing_ac_modulation_in_ptc) then   ! useless now
    localmaster=master
    call ass(scdadd%AC%x(1))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    scdadd%ac%x(1)=s1%M%V(C_%ND2-1) +s2%AC%x(1)
    master=localmaster
    localmaster=master
    call ass(scdadd%AC%x(2))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    scdadd%ac%x(2)=s1%M%V(C_%ND2) +s2%AC%x(2)
    master=localmaster
    localmaster=master
    call ass(scdadd%AC%om)
!    call ass(scdadd%AC%t)
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    scdadd%AC%om=s2%AC%om
    scdadd%AC%t=s2%AC%t
    master=localmaster
    !    endif

    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(scdadd%s(J)%x(i))
          scdadd%s(J)%x(i)=S1%S%s(I,J)
          master=localmaster
       ENDDO

    ENDDO

    scdadd%e_ij=s1%e_ij
  END FUNCTION scdadd


  FUNCTION daddsc( S1,S2  )
    implicit none
    TYPE (probe_8) daddsc
    TYPE (damapspin), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j,nd2t,dc
    type(taylor) d




    !   call ass(daddsc)
    daddsc%u=my_false
    daddsc%E_ij=0.0_dp

    if(doing_ac_modulation_in_ptc) then
       dc=2
    else
       dc=0
    endif
    if(c_%ndpt/=0) then                       ! delta-l without cavity
       nd2t=6
    elseif(c_%ndpt==0.and.c_%nd2==4+dc) then     !  no delta-l maybe delta as a parameter
       nd2t=4
    elseif(c_%ndpt==0.and.c_%nd2==6+dc) then      !   delta-l with cavity
       nd2t=6
    else
       write(6,*) " internal PTC error in daddsc o_tree_element.f90 "
       stop 666
    endif

    call alloc(d)
    do i=1,nd2t                 !   from 1-4 or 1-6 (if ndpt=0)
       localmaster=master
       call ass(daddsc%x(i))
       !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
       daddsc%x(i)=s1%M%V(i)+s2%x(i)-(s1%M%V(i).sub.'0')
       master=localmaster
    enddo

    do i=nd2t+1,6
       localmaster=master
       call ass(daddsc%x(i))
       if((c_%npara==5+dc).AND.I==5) then   ! npr
          daddsc%x(i)=s2%x(i)+(1.0_dp.mono.c_%npara)
       else
          daddsc%x(i)=s2%x(i)
       endif
       master=localmaster
    enddo

    !    if(doing_ac_modulation_in_ptc) then   ! useless now
    localmaster=master
    call ass(daddsc%AC%x(1))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    daddsc%ac%x(1)=s1%M%V(C_%ND2-1) +s2%AC%x(1)
    master=localmaster
    localmaster=master
    call ass(daddsc%AC%x(2))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    daddsc%ac%x(2)=s1%M%V(C_%ND2) +s2%AC%x(2)
    master=localmaster
    localmaster=master
    call ass(daddsc%AC%om)
!    call ass(daddsc%AC%t)
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    daddsc%AC%om=s2%AC%om
    daddsc%AC%t=s2%AC%t
    master=localmaster
    !    endif

    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(daddsc%s(J)%x(i))
          daddsc%s(J)%x(i)=S1%S%s(I,J)
          master=localmaster
       ENDDO

    ENDDO

    daddsc%e_ij=s1%e_ij


  END FUNCTION daddsc

  subroutine print_DASPIN(DS,MF,prec)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: DS
    real(dp), optional :: prec
    INTEGER MF,I,J
    logical(lp) spin_in,rad_in

    !    WRITE(MF,*) " ORBIT "
    !     WRITE(MF,*) DS%X(1:3)
    !     WRITE(MF,*) DS%X(4:6)
    WRITE(MF,*) c_%nd2 ," DIMENSIONAL ORBITAL MAP "
    CALL PRINT(DS%M,MF,prec)
    call check_spin(DS,spin_in)
    if(spin_in) then
       WRITE(MF,*) " SPIN MATRIX "
       DO I=1,3
          DO J=1,3
             WRITE(MF,*) " SPIN MATRIX COMPONENT ",I,J
             CALL PRINT(DS%S%s(I,J),MF,prec)
          ENDDO
       ENDDO
    else

       WRITE(MF,*) "NO SPIN MATRIX OR IDENTITY  "

    endif
    call check_rad(DS%e_ij,rad_in)
    if(rad_in) then
       WRITE(MF,*) " STOCHASTIC KICK  "
       DO I=1,6
          DO J=1,6
             WRITE(MF,*) "  STOCHASTIC KICK COMPONENT ",I,J
             write(mf,*) ds%e_ij(i,j)
          ENDDO
       ENDDO
    else
       WRITE(MF,*) "NO STOCHASTIC KICK  "
    endif

  END subroutine print_DASPIN

  subroutine print_probe8(DS,MF)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: DS
    INTEGER MF,I,j
    logical(lp) rad_in

    WRITE(MF,*) " ORBIT "
    do i=1,6
       write(mf,*) ' Variable ',i
       call print(ds%x(i),mf)
    enddo
    !    WRITE(MF,*) " SPIN 0 "
    !    do i=1,3
    !       write(mf,*) ' Spin Variable ',i
    !       call print(ds%s(0)%x(i),mf)
    !    enddo
    WRITE(MF,*) " SPIN X "
    call print(ds%s(1),mf)
    !    do i=1,3
    !       write(mf,*) ' Spin Variable ',i
    !       call print(ds%s(1)%x(i),mf)
    !    enddo
    WRITE(MF,*) " SPIN Y "
    call print(ds%s(2),mf)
    !    do i=1,3
    !       write(mf,*) ' Spin Variable ',i
    !       call print(ds%s(2)%x(i),mf)
    !    enddo
    WRITE(MF,*) " SPIN Z "
    call print(ds%s(3),mf)
    !    do i=1,3
    !       write(mf,*) ' Spin Variable ',i
    !       call print(ds%s(3)%x(i),mf)
    !    enddo


    call check_rad(DS%e_ij,rad_in)



    if(rad_in) then

       WRITE(MF,*) " STOCHASTIC KICK  "

       do i=1,6
          do j=1,6
             write(6,*) i,j,ds%e_ij(i,j)
          enddo
       enddo
    else

       WRITE(MF,*) "NO STOCHASTIC KICK  "

    endif
    if(doing_ac_modulation_in_ptc) then
       call print(ds%ac,mf)
    else
       WRITE(MF,*) "NO MODULATION  "

    endif

  END subroutine print_probe8

  subroutine print_rf_phasor_8(S,MF)
    implicit none
    TYPE(rf_phasor_8), INTENT(INOUT) :: s
    INTEGER MF,I

    write(mf,*) ' AC INFORMATION '
    call print(s%om,mf)
    call print(s%t,mf)
    do i=1,2
       call print(s%x(i),mf)
    enddo

  END subroutine print_rf_phasor_8


  subroutine print_spinor_8(S,MF)
    implicit none
    TYPE(spinor_8), INTENT(INOUT) :: s
    INTEGER MF,I

    do i=1,3
       write(mf,*) ' Spin Variable ',i
       call print(s%x(i),mf)
    enddo

  END subroutine print_spinor_8

  subroutine print_res_spinor_8(S,MF)
    implicit none
    TYPE(res_spinor_8), INTENT(INOUT) :: s
    INTEGER MF,I

    do i=1,3
       write(mf,*) ' Spin Variable ',i
       call print(s%x(i),mf)
    enddo

  END subroutine print_res_spinor_8

  subroutine read_spinor_8(S,MF)
    implicit none
    TYPE(spinor_8), INTENT(INOUT) :: s
    INTEGER MF,I
    character*255 line
    type(taylor) t
    call alloc(t)

    do i=1,3
       read(mf,'(a255)') line
       call read(t,mf)
       s%x(i)=morph(t)
    enddo

    call kill(t)

  END subroutine read_spinor_8

  subroutine READ_DASPIN(DS,MF,file)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: DS
    INTEGER MF1,I,J,nd2
    CHARACTER*20 LINE
    CHARACTER(*), optional :: file
    INTEGER, optional :: MF
    TYPE(TAYLOR) T

    if(present(mf)) mf1=mf
    if(present(file)) then
       call kanalnummer(mf1)
       open(unit=mf1,file=file)
    endif
    CALL ALLOC(T)

    !    READ(MF,*) LINE
    !     READ(MF,*) DS%X(1:3)
    !     READ(MF,*) DS%X(4:6)
    READ(MF1,*) nd2, LINE(1:7)
    if(nd2>=c_%nd2) then
       CALL READ(DS%M,MF1)
       do i=c_%nd2+1,nd2
          call rea(t,mf1)
       enddo
    else
       do i=1,nd2
          call read(ds%m%v(i),mf1)
       enddo
       do i=nd2+1,c_%nd2
          ds%m%v(i)=1.0_dp.mono.i
       enddo
    endif
    READ(MF1,*) LINE
    DO I=1,3
       DO J=1,3
          READ(MF1,*) LINE
          CALL READ(T,MF1)
          DS%S%s(I,J)=T     !MORPH(T)
       ENDDO
    ENDDO
    READ(MF1,*) LINE
    DO I=1,3
       DO J=1,3
          READ(MF1,*) LINE
          read(mf,*) ds%e_ij(i,j)
       ENDDO
    ENDDO
    if(present(file)) close(mf1)

    CALL KILL(T)
  END subroutine READ_DASPIN

  subroutine read_probe8(DS,MF)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: DS
    INTEGER MF,I
    character*120 line
    type(taylor) t
    call alloc(t)

    read(mf,*) line
    do i=1,6
       read(mf,*) line
       call read(t,mf)
       ds%x(i)=morph(t)
    enddo
    do i=ISPIN0R,ISPIN1R
       read(mf,'(a120)') line
       call read(ds%s(i),mf)
    enddo

    call kill(t)

  END subroutine read_probe8



!!! NOT USED NOW
  subroutine CHECK_RES_ORBIT(J,NRES,M,SKIP)
    implicit none
    INTEGER M(:,:),NRES
    LOGICAL(lp) SKIP,SKIP1,SKIP2
    INTEGER I,K,J(:)

    SKIP=.FALSE.
    IF(NRES==0) RETURN

    DO I=1,NRES
       SKIP1=.TRUE.
       SKIP2=.TRUE.
       DO K=1,C_%ND
          SKIP1=((J(2*K)-J(2*K-1))==M(k,i)).AND.SKIP1
          SKIP2=((J(2*K)-J(2*K-1))==-M(k,i)).AND.SKIP2
       ENDDO
       IF(SKIP1.OR.SKIP2) THEN
          SKIP=.TRUE.
          RETURN
       ENDIF
    ENDDO

  END subroutine CHECK_RES_ORBIT



  subroutine find_n_thetar(s0,n0)
    implicit none
    real(dp),intent(in) :: s0(3,3)
    real(dp)  theta0,n0(3)
    real(dp)  det,ss(3,3),detm
    integer i,is,j





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
    detm=(ss(1,1)*ss(2,2)-ss(1,2)*ss(2,1))
    if(abs(detm)>=abs(det)) then
       det=detm
       is=3
    endif

    n0(is)=1.0_dp
    if(is==1) then
       n0(2)=(-ss(3,3)*ss(2,1)+ss(2,3)*ss(3,1))/det
       n0(3)=(-ss(2,2)*ss(3,1)+ss(2,1)*ss(3,2))/det
    elseif(is==2) then
       n0(1)=(-ss(3,3)*ss(1,2)+ss(3,2)*ss(1,3))/det
       n0(3)=(-ss(1,1)*ss(3,2)+ss(1,2)*ss(3,1))/det
    else
       n0(1)=(-ss(2,2)*ss(1,3)+ss(2,3)*ss(1,2))/det
       n0(2)=(-ss(1,1)*ss(2,3)+ss(1,3)*ss(2,1))/det
    endif


    theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)

    do i=1,3
       n0(i)=n0(i)/theta0
    enddo




  end subroutine find_n_thetar

  subroutine find_n_thetap(s0,n0)
    implicit none
    type(real_8),intent(in) :: s0(3,3)

    type(real_8)  theta0,n0(3)
    type(real_8)  det,ss(3,3),detm
    integer i,is,j


    call alloc(det)
    call alloc(detm)
    call alloc(theta0)
    call alloc_33(ss)


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



    n0(is)=1.0_dp
    if(is==1) then
       n0(2)=(-ss(3,3)*ss(2,1)+ss(2,3)*ss(3,1))/det
       n0(3)=(-ss(2,2)*ss(3,1)+ss(2,1)*ss(3,2))/det
    elseif(is==2) then
       n0(1)=(-ss(3,3)*ss(1,2)+ss(3,2)*ss(1,3))/det
       n0(3)=(-ss(1,1)*ss(3,2)+ss(1,2)*ss(3,1))/det
    else
       n0(1)=(-ss(2,2)*ss(1,3)+ss(2,3)*ss(1,2))/det
       n0(2)=(-ss(1,1)*ss(2,3)+ss(1,3)*ss(2,1))/det
    endif


    theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)


    do i=1,3
       n0(i)=n0(i)/theta0
    enddo


    call kill(det)
    call kill(detm)
    call kill(theta0)
    call kill_33(ss)

  end subroutine find_n_thetap

  !  find exponent of rotation routines
  subroutine find_axisp(ds,H_axis)  !
    implicit none
    type(damapspin),intent(in) :: ds
    type(spinor_8),intent(inout) :: H_axis

    type(real_8)  theta0,n0(3)
    type(real_8)  det,ss(3,3),detm
    integer i,is,j


    call alloc(det)
    call alloc(detm)
    call alloc(theta0)
    call alloc_33(ss)
    call alloc(n0)


    do i=1,3
       do j=1,3
          ss(i,j)=ds%s%s(i,j)
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



    n0(is)=1.0_dp
    if(is==1) then
       n0(2)=(-ss(3,3)*ss(2,1)+ss(2,3)*ss(3,1))/det
       n0(3)=(-ss(2,2)*ss(3,1)+ss(2,1)*ss(3,2))/det
    elseif(is==2) then
       n0(1)=(-ss(3,3)*ss(1,2)+ss(3,2)*ss(1,3))/det
       n0(3)=(-ss(1,1)*ss(3,2)+ss(1,2)*ss(3,1))/det
    else
       n0(1)=(-ss(2,2)*ss(1,3)+ss(2,3)*ss(1,2))/det
       n0(2)=(-ss(1,1)*ss(2,3)+ss(1,3)*ss(2,1))/det
    endif

    theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)

    do i=1,3
       n0(i)=n0(i)/theta0
    enddo

    h_axis%x(1)=n0(1);  h_axis%x(2)=n0(2) ; h_axis%x(3)=n0(3);





    call kill(det)
    call kill(detm)
    call kill(theta0)
    call kill_33(ss)
    call kill(n0)

  end subroutine find_axisp

  subroutine find_perp_basisp(y_axis,x_axis,z_axis)  !
    implicit none
    type(spinor_8),intent(inout) :: y_axis,x_axis,z_axis
    integer i,is
    type(real_8) norm

    call alloc(norm)

    is=1
    if(abs(y_axis%x(is))>abs(y_axis%x(2))) then
       is=2
    endif
    if(abs(y_axis%x(is))>abs(y_axis%x(3))) then
       is=3
    endif

    !  now is = smallest

    x_axis=is

    x_axis=x_axis-(x_axis.dot.y_axis)*y_axis

    norm=sqrt(x_axis.dot.x_axis)
    x_axis=(1.d0/norm)*x_axis



    z_axis=x_axis*y_axis
    norm=sqrt(z_axis.dot.z_axis)
    z_axis=(1.d0/norm)*z_axis



    call kill(norm)

  end subroutine find_perp_basisp

  subroutine find_exponent_jet_p(ds,h_axis)  !
    implicit none
    type(damapspin),intent(inout) :: ds
    type(spinor_8),intent(inout) :: h_axis
    type(damapspin) s,sf,st,SA
    integer i,nmax
    real(dp) eps1,EPS0
    LOGICAL DONOT
    nmax=1000
    call alloc(s,sf,st,SA)

    DONOT=.FALSE.
    EPS0=MYBIG

    st=1
    s=ds
    s%m=1
    sf=0
    do i=1,3
       s%s%s(i,i)=s%s%s(i,i)-1.0_dp
    enddo

    st=s
    do i=1,nmax
       SA=SF
       sf=sf + (1.0_dp/i)*st
       st= s*st
       st=(-1.0_dp)*st
       SA=SA+((-1.0_dp)*SF)
       eps1=FULL_ABS(SA)
       IF(EPS1>deps_tracking) THEN
          EPS0=EPS1
       ELSE
          IF(EPS1>=EPS0.or.EPS1<=PUNY) THEN
             DONOT=.TRUE.
          ELSE
             EPS0=EPS1
          ENDIF
       ENDIF

       IF(DONOT) EXIT

    enddo

    if(i>nmax-2) then
       write(6,*) "Did not converge in find_exponent_jet_p"
       stop 666
    endif

    h_axis%x(1)=sf%s%s(3,2)
    h_axis%x(2)=sf%s%s(1,3)
    h_axis%x(3)=sf%s%s(2,1)

    !etienne


    call kill(s,sf,st,SA)
  end subroutine find_exponent_jet_p

  subroutine find_exponentp(ds,y_axis,x_axis,z_axis,h_axis,spin_tune)  !
    implicit none
    type(damapspin),intent(inout) :: ds
    type(spinor_8),intent(inout) :: y_axis,x_axis,z_axis,h_axis
    type(spinor_8) s
    type(real_8) cos,sin,theta
    type(real_8) , optional :: spin_tune

    call alloc(s)
    call alloc(cos,sin,theta)

    call find_axis(ds,y_axis)
    !  etienne

    call find_perp_basis(y_axis,x_axis,z_axis)



    s=ds*x_axis


    cos=s.dot.x_axis
    sin=-(s.dot.z_axis)

    theta=clockwise*atan2(sin,cos)
    !if(force_positive.and.theta<zero) theta = theta + twopi    !!!! allow negative theta

    h_axis=(clockwise*theta)*y_axis

    if(present(spin_tune)) then
       spin_tune=(1.0_dp/twopi)*theta
    endif

    call kill(cos,sin,theta)
    call kill(s)
  end subroutine find_exponentp

  subroutine find_exponent_only(ds,h_axis,spin_tune)  !
    implicit none
    type(damapspin),intent(inout) :: ds
    type(spinor_8),intent(inout) :: h_axis
    type(spinor_8)  y_axis,x_axis,z_axis
    type(spinor_8) s
    type(real_8) cos,sin,theta
    type(real_8) , optional :: spin_tune

    call alloc(s)
    call alloc(cos,sin,theta)
    call alloc(y_axis)
    call alloc(x_axis)
    call alloc(z_axis)
    call find_axis(ds,y_axis)
    !  etienne

    call find_perp_basis(y_axis,x_axis,z_axis)



    s=ds*x_axis


    cos=s.dot.x_axis
    sin=-(s.dot.z_axis)

    theta=clockwise*atan2(sin,cos)
    !if(force_positive.and.theta<zero) theta = theta + twopi    !!!! allow negative theta

    h_axis=(clockwise*theta)*y_axis

    if(present(spin_tune)) then
       spin_tune=(1.0_dp/twopi)*theta
    endif

    call kill(y_axis)
    call kill(x_axis)
    call kill(z_axis)
    call kill(cos,sin,theta)
    call kill(s)
  end subroutine find_exponent_only

  subroutine remove_y_rot(as_xyz,r_y)
    implicit none
    type(damapspin), intent(inout) :: as_xyz
    type(damapspin), optional , intent(inout) ::r_y
    type(damapspin) temp,as_y,as_nl,rot_y
    type(spinor_8) n_expo,n_tune
    type(taylorresonance) tr
    type(taylor) t
    integer i,j
    integer  nmax
    real(dp) c,eps,norm1,norm2
    logical check
!!!  original as_xyz = as_xyz*r_y = a_y*a_nl*r_y  on exit
    check=.true.
    eps=1.d-1
    nmax=1000

    call alloc(n_expo)
    call alloc(n_tune)
    call alloc(temp,as_y,as_nl,rot_y)
    call alloc(tr)
    call alloc(t)

    as_y=as_xyz

    norm1=mybig
    rot_y=1

    check=.true.
    norm1=mybig
    do i=1,nmax
       call find_exponent_only(as_y,n_expo)
       !     call dalog_spinor_8(as_y,n_expo)
       norm2=0.0_dp
       n_tune%x(1)=0.0_dp
       n_tune%x(3)=0.0_dp
       ! do j=1,3
       j=2     ! do loop was a theory error
       t=n_expo%x(j)
       tr=t
       call cfu(tr%cos,phase_shift,tr%cos)
       tr%sin=0.0_dp
       t=tr
       n_tune%x(j)=t
       !  enddo
       temp=exp(n_tune)

       rot_y=temp*rot_y
       call inv_as(temp%s%s,temp%s%s)
       as_y=as_y*temp
       if(check) then
          if(norm2<eps) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2

    enddo
    if(i>nmax-10) then
       write(6,*) "no convergence in remove_y_rot "
       stop 1067
    endif
    as_xyz=as_y
    if(present(r_y)) r_y=rot_y

    call kill(n_expo)
    call kill(temp,as_y,as_nl,rot_y)
    call kill(tr)
    call kill(n_tune)
    call kill(t)

  end subroutine remove_y_rot

  subroutine remove_y_rot0(as_xyz,a_y,a_nl,r_y)
    implicit none
    type(damapspin), intent(inout) :: as_xyz
    type(damapspin), optional , intent(inout) :: a_y,a_nl,r_y
    type(damapspin) temp,as_y,as_nl,rot_y
    type(spinor_8) n_expo,n_tune
    type(taylorresonance) tr
    type(taylor) t
    integer i,j
    integer  nmax
    real(dp) c,eps,norm1,norm2
    logical check
!!!  original as_xyz = as_xyz*r_y = a_y*a_nl*r_y  on exit
    check=.true.
    eps=1.d-5
    nmax=1000

    call alloc(n_expo)
    call alloc(n_tune)
    call alloc(temp,as_y,as_nl,rot_y)
    call alloc(tr)
    call alloc(t)

    call factor_parameter_dependent_s0(as_xyz,as_y,as_nl,n_expo,1)

    norm1=mybig
    rot_y=1
    do i=1,nmax
       call find_exponent_only(as_y,n_expo)
       n_expo%x(1)=0.0_dp
       n_expo%x(3)=0.0_dp
       temp=exp(n_expo)
       rot_y=temp*rot_y
       call inv_as(temp%s%s,temp%s%s)
       as_y=as_y*temp
       call norm_spinor_8(n_expo,norm2)
       if(check) then
          if(norm2<eps) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2
    enddo
    if(i>nmax-10) then
       write(6,*) "no convergence in daexplogp "
       stop 1066
    endif
    as_nl=rot_y*as_nl*rot_y**(-1)

    check=.true.
    norm1=mybig
    do i=1,nmax
       call dalog_spinor_8(as_nl,n_expo)
       norm2=0.0_dp
       n_tune%x(1)=0.0_dp
       n_tune%x(3)=0.0_dp
       ! do j=1,3
       j=2     ! do loop was a theory error
       t=n_expo%x(j)
       tr=t
       call cfu(tr%cos,phase_shift,tr%cos)
       tr%sin=0.0_dp
       t=tr
       n_tune%x(j)=t
       norm2=norm2+full_abs(t)
       !  enddo
       temp=exp(n_tune)

       rot_y=temp*rot_y
       call inv_as(temp%s%s,temp%s%s)
       as_nl=as_nl*temp
       if(check) then
          if(norm2<eps) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2

    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in daexplogp "
       stop 1067
    endif
    as_xyz=as_y*as_nl
    if(present(r_y)) r_y=rot_y
    if(present(a_y)) a_y=as_y
    if(present(a_nl)) a_nl=as_nl

    call kill(n_expo)
    call kill(temp,as_y,as_nl,rot_y)
    call kill(tr)
    call kill(n_tune)
    call kill(t)

  end subroutine remove_y_rot0

  subroutine daexplogp(h_axis,ds)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(spinor_8), INTENT(IN) :: h_axis
    integer  nmax
    integer i
    TYPE(damapspin) dh,dhn,dr
    real(dp) c,eps,norm1,norm2
    logical check

    check=.true.
    eps=1.d-5
    nmax=1000

    call alloc(dh)
    call alloc(dhn)
    call alloc(dr)

    !  this only works with a da-map
    ds=1
    dh%m=1
    dh%s%s(2,1)=h_axis%x(3)
    dh%s%s(1,3)=h_axis%x(2)
    dh%s%s(3,2)=h_axis%x(1)
    dh%s%s(1,2)=-h_axis%x(3)
    dh%s%s(3,1)=-h_axis%x(2)
    dh%s%s(2,3)=-h_axis%x(1)
    dhn=1
    c=1.0_dp
    norm1=mybig
    do i=1,nmax
       dhn=dhn*dh
       c=c/i
       dr=ds
       ds=ds+c*dhn
       dr=ds+(-1.0_dp)*dr
       call norm_damapspin(dr,norm2)
       if(check) then
          if(norm2<eps) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2
    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in daexplogp "
       stop 1066
    endif

    call kill(dh)
    call kill(dhn)
    call kill(dr)

  end subroutine daexplogp

  subroutine norm_damapspin(ds,norm)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    real(dp) norm
    integer i,j

    norm=0.0_dp

    do i=1,3
       do j=1,3
          norm=norm+full_abs( ds%s%s(i,j) )
       enddo
    enddo

  end subroutine norm_damapspin

  subroutine norm_spinor_8(s,norm)
    implicit none
    TYPE(spinor_8), INTENT(INout) :: s
    real(dp) norm
    integer i

    norm=0.0_dp

    do i=1,3
       norm=norm+full_abs( s%x(i) )
    enddo

  end subroutine norm_spinor_8

  ! end of  find exponent of rotation routines


  subroutine find_ar(n2,a)
    implicit none
    real(dp)  n2(3),n1(3),n3(3)
    real(dp) a(3,3),s,n
    integer i,is

    ! here we find smallest value of n2
    is=2
    if(abs(n2(1))< abs(n2(2))) is=1

    if(is==1) then
       if(abs(n2(3))<abs(n2(1))) is=3
    else
       if(abs(n2(3))<abs(n2(2))) is=3
    endif

    !  put n1 in along that value
    n1=0.0_dp
    n1(is)=1.0_dp

    s=n2(is)*n1(is)

    n=0.0_dp
    do i=1,3
       n1(i)=n1(i)-s*n2(i)
       n=n1(i)**2+n
    enddo
    n1=n1/sqrt(n)

    n3(1)=n1(2)*n2(3)-n1(3)*n2(2)
    n3(2)=n1(3)*n2(1)-n1(1)*n2(3)
    n3(3)=n1(1)*n2(2)-n1(2)*n2(1)

    n=0.0_dp
    do i=1,3
       n=n3(i)**2+n
    enddo
    n3=n3/sqrt(n)
    ! spin_normal_position

    !    a=zero
    if(spin_normal_position==2) then
       a(:,1)=n1
       a(:,2)=n2
       a(:,3)=n3
    elseif(spin_normal_position==3) then
       a(:,2)=n1
       a(:,3)=n2
       a(:,1)=n3
    else
       a(:,3)=n1
       a(:,1)=n2
       a(:,2)=n3
    endif


  end subroutine find_ar

  subroutine find_ap(n2,a)
    implicit none
    type(real_8)  n2(3),n1(3),n3(3)
    type(real_8)  a(3,3),s,n
    real(dp) aa(3,3)
    integer j
    integer i,is
    call alloc(n1,3)
    call alloc(n3,3)
    call alloc(s,n)
    ! here we find smallest value of n2
    is=2
    if(abs(n2(1))< abs(n2(2))) is=1

    if(is==1) then
       if(abs(n2(3))<abs(n2(1))) is=3
    else
       if(abs(n2(3))<abs(n2(2))) is=3
    endif

    !  put n1 in along that value
    do i=1,3
       n1(i)=0.0_dp
    enddo
    n1(is)=1.0_dp

    s=n2(is)*n1(is)

    n=0.0_dp
    do i=1,3
       n1(i)=n1(i)-s*n2(i)
       n=n1(i)**2+n
    enddo
    do i=1,3
       n1(i)=n1(i)/sqrt(n)
    enddo

    n3(1)=n1(2)*n2(3)-n1(3)*n2(2)
    n3(2)=n1(3)*n2(1)-n1(1)*n2(3)
    n3(3)=n1(1)*n2(2)-n1(2)*n2(1)

    n=0.0_dp
    do i=1,3
       n=n3(i)**2+n
    enddo
    do i=1,3
       n3(i)=n3(i)/sqrt(n)
    enddo

    if(spin_normal_position==2) then
       do i=1,3
          a(i,1)=n1(i)
          a(i,2)=n2(i)
          a(i,3)=n3(i)
       enddo
    elseif(spin_normal_position==3) then
       do i=1,3
          a(i,2)=n1(i)
          a(i,3)=n2(i)
          a(i,1)=n3(i)
       enddo
    else
       do i=1,3
          a(i,3)=n1(i)
          a(i,1)=n2(i)
          a(i,2)=n3(i)
       enddo
    endif


    call kill(n1,3)
    call kill(n3,3)
    call kill(s,n)
  end subroutine find_ap



  subroutine copy_damap_matrix(mi,a)
    implicit none
    type(taylor), intent(inout) :: a(:,:)
    type(damap), intent(in) :: mi
    type(damap) m

    integer i,j,nt
    logical doflip
    integer, allocatable :: jl(:)

    call alloc(m)

    m=mi

    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       perform_flip=.false.
       call flip_damap(m,m)
       doflip=.true.
    else
       doflip=.false.
    endif
    nt=c_%nd2
    if(c_%ndpt/=0)nt=nt-2

    allocate(jl(nt))
    jl=0
    do i=1,min(nt,size(a,dim=1))
       do j=1,min(nt,size(a,dim=2))
          jl(j)=1
          a(i,j)=m%v(i).par.jl
          jl(j)=0
       enddo
    enddo

    if(doflip) then
       call flip_damap(m,m)
       perform_flip=.true.
    endif

    call kill(m)
    deallocate(jl)

  end subroutine copy_damap_matrix

  subroutine copy_matrix_matrix(ma,a)
    implicit none
    type(taylor), intent(inout) :: a(:,:)
    type(taylor), intent(in) :: ma(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          a(i,j)=ma(i,j)
       enddo
    enddo

  end subroutine copy_matrix_matrix



  subroutine alloc_33t(a)
    implicit none
    type(taylor) a(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          call alloc(a(i,j))
       enddo
    enddo

  end subroutine alloc_33t

  subroutine kill_33t(a)
    implicit none
    type(taylor) a(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          call kill(a(i,j))
       enddo
    enddo

  end subroutine kill_33t

  subroutine zero_33t(a,r)
    implicit none
    type(taylor) a(:,:) ! taylor
    !    type(taylor) a(3,3)
    integer i,j
    real(dp), optional :: r

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          a(i,j)=0.0_dp
       enddo
       if(present(r)) a(i,i)=1.0_dp
    enddo
  end subroutine zero_33t

  subroutine zero_33p(a,r)
    implicit none
    type(real_8) a(:,:) ! taylor
    !    type(taylor) a(3,3)
    integer i,j
    real(dp), optional :: r

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          a(i,j)=0.0_dp
       enddo
       if(present(r)) a(i,i)=1.0_dp
    enddo

  end subroutine zero_33p

  subroutine alloc_33p(a)
    implicit none
    type(real_8) a(:,:) ! taylor
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          call alloc(a(i,j))
       enddo
    enddo

  end subroutine alloc_33p

  subroutine kill_33p(a)
    implicit none
    type(real_8) a(:,:)
    integer i,j

    do i=1,size(a,dim=1)
       do j=1,size(a,dim=2)
          call kill(a(i,j))
       enddo
    enddo

  end subroutine kill_33p

  subroutine matmul_33(m,n,mo,sc)
    implicit none
    type(taylor) m(:,:),n(:,:),mo(:,:)
    type(taylor), allocatable :: a(:,:)
    real(dp), optional :: sc
    real(dp) sc0
    integer i,j,k
    sc0=1.0_dp
    allocate(a(size(m,dim=1),size(n,dim=2)))

    call alloc_33(a)

    do i=1,size(m,dim=1)
       do j=1,size(m,dim=2)
          do k=1,size(n,dim=2)
             a(i,k)=m(i,j)*n(j,k)+a(i,k)
          enddo
       enddo
    enddo

    if(present(sc)) sc0=sc
    do i=1,size(mo,dim=1)
       do j=1,size(mo,dim=2)
          mo(i,j)=sc0*a(i,j)
       enddo
    enddo
    call kill_33(a)
    deallocate(a)
  end subroutine matmul_33

  subroutine matmul_3344(m,n,mo)
    implicit none
    type(taylor) m(:,:),n(:,:),mo(:,:)
    type(taylor), allocatable :: a(:,:)
    integer i,j,k

    allocate(a(size(m,dim=1),size(n,dim=2)))

    call alloc_33(a)

    do i=1,size(m,dim=1)
       do j=1,size(m,dim=2)
          do k=1,size(n,dim=2)
             a(i,k)=m(i,j)*n(j,k)+a(i,k)
          enddo
       enddo
    enddo


    do i=1,size(mo,dim=1)
       do j=1,size(mo,dim=2)
          mo(i,j)=a(i,j)
       enddo
    enddo
    call kill_33(a)
    deallocate(a)
  end subroutine matmul_3344

  subroutine matmulp(m,n,mo)
    implicit none
    type(real_8) m(3,3),n(3,3),mo(3,3),a(3,3)
    integer i,j,k


    call alloc_33(a)

    do i=1,3
       do j=1,3
          do k=1,3
             a(i,k)=m(i,j)*n(j,k)+a(i,k)
          enddo
       enddo
    enddo


    do i=1,3
       do j=1,3
          mo(i,j)=a(i,j)
       enddo
    enddo
    call kill_33(a)

  end subroutine matmulp

  subroutine smatmulp(sc,m,n,mo)
    implicit none
    type(real_8) m(3,3),n(3,3),mo(3,3),a(3,3)
    integer i,j
    real(dp) sc

    call alloc_33(a)

    do i=1,3
       do j=1,3
          a(i,j)=m(i,j)*sc+n(i,j)
       enddo
    enddo


    do i=1,3
       do j=1,3
          mo(i,j)=a(i,j)
       enddo
    enddo
    call kill_33(a)

  end subroutine smatmulp

  subroutine smatp(sc,m,mo)
    implicit none
    type(real_8) m(3,3),mo(3,3),a(3,3)
    integer i,j
    real(dp) sc

    call alloc_33(a)

    do i=1,3
       do j=1,3
          a(i,j)=m(i,j)*sc
       enddo
    enddo


    do i=1,3
       do j=1,3
          mo(i,j)=a(i,j)
       enddo
    enddo
    call kill_33(a)

  end subroutine smatp


  subroutine make_unitary_p(m)
    !  inverse 3x3 rrotation
    implicit none
    integer i, nmax
    type(real_8) m(3,3)
    type(real_8) mt(3,3),id(3,3)
    real(dp) r,eps,rb
    logical doit

    doit=.true.
    nmax=100
    eps=1.e-4_dp
    call alloc_33(mt)
    call alloc_33(id)

    do i=1,3
       id(i,i)=1.5e0_dp
    enddo

    rb=mybig

    do i=1,nmax
       call transpose_p(m,mt)

       call matmulp(m,mt,mt)

       call smatp(-0.5_dp,mt,mt)

       call smatmulp(1.0_dp,mt,id,mt)

       call matmulp(mt,m,m)

       call transpose_p(m,mt)
       call matmulp(mt,m,mt)
       call smatmulp(-1.0_dp/1.5e0_dp,id,mt,mt)
       call absolute_p(mt,r)
       if(r<eps.and.doit) then
          doit=.false.
          rb=r
       else
          if(rb<=r) exit
          rb=r
       endif
    enddo

    if(i>=nmax) then
       write(6,*) " Too many iterations "
       write(6,*) " Perhaps did not converged in make_unitary_p "

    endif

    call kill_33(id)
    call kill_33(mt)

  end subroutine make_unitary_p

  subroutine check_unitary_p(m,r,rt,cc)
    !  inverse 3x3 rrotation
    implicit none
    integer i, j,c
    type(real_8) m(3,3)
    type(real_8) mt(3,3)
    real(dp) r,rt
    integer, optional :: cc

    c=c_%no+1
    if(present(cc)) c=cc

    call alloc_33(mt)

    r=0.0_dp
    rt=0.0_dp
    call transpose_p(m,mt)

    call matmulp(m,mt,mt)

    do i=1,3
       do j=1,3

          r=r+abs(mt(i,j))
          rt=rt+full_abs(mt(i,j).cut.c)
       enddo
    enddo



    call kill_33(mt)

  end subroutine check_unitary_p


  subroutine absolute_p(m,r)
    implicit none
    integer i,j
    type(real_8) m(3,3)
    real(dp) r

    r=0.0_dp

    do i=1,3
       do j=1,3
          r=r+full_abs(m(i,j))
       enddo
    enddo

  end subroutine absolute_p

  subroutine transpose_p(m,mi)
    !  inverse 3x3 rrotation
    implicit none
    integer i,j
    type(real_8), INTENT(INOUT) :: m(3,3),mi(3,3)
    type(real_8) mt(3,3)

    call alloc_33(mt)

    call smatp(1.0_dp,m,mt)

    do i=1,3
       do j=1,3
          mi(j,i)=mt(i,j)
       enddo
    enddo

    call kill_33(mt)

  end subroutine transpose_p

  subroutine inv_asr(m,mi)
    !  inverse 3x3 rrotation
    implicit none
    integer i,j
    real(dp) m(3,3),mi(3,3)
    real(dp)  n(3,3)


    do i=1,3
       do j=1,3
          n(j,i)=m(i,j)
       enddo
    enddo

    do i=1,3
       do j=1,3
          mi(i,j)=n(i,j)
       enddo
    enddo

  end subroutine inv_asr


  subroutine inv_asp(m,mi)
    !  inverse 3x3 rrotation
    implicit none
    integer i,j
    type(real_8) m(3,3),mi(3,3)
    type(real_8)  n(3,3)

    call alloc_33(n)

    do i=1,3
       do j=1,3
          n(j,i)=m(i,j)
       enddo
    enddo

    do i=1,3
       do j=1,3
          mi(i,j)=n(i,j)
       enddo
    enddo

    call kill_33(n)
  end subroutine inv_asp

  subroutine  trans_mat(m,ri,mi)
    implicit none
    integer i,j
    type(real_8) m(3,3),mi(3,3)
    type(real_8)  n(3,3)
    type(damap)  ri
    type(real_8)  r

    call alloc_33(n)
    call alloc(r)

    do i=1,3
       do j=1,3
          if(c_%nd2/=0) then
             r=(m(i,j)*ri)
          else
             r=m(i,j)
          endif
          n(i,j)=r    !morph(r)
          ! advances by r**(-1)
       enddo
    enddo


    do i=1,3
       do j=1,3
          mi(i,j)=n(i,j)
       enddo
    enddo

    call kill(r)
    call kill_33(n)
  end subroutine trans_mat






  subroutine ALLOC_DASPIN(D)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: D
    INTEGER I,J
    !    D%s0=ZERO
    CALL ALLOC(D%M)
    DO I=1,3
       DO J=1,3
          CALL ALLOC(D%S%s(I,J))
          d%e_ij(i,j)=0.0_dp
       ENDDO
    ENDDO

  END    subroutine ALLOC_DASPIN

  subroutine KILL_DASPIN(D)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: D
    INTEGER I,J

    !    D%s0=ZERO
    CALL KILL(D%M)
    DO I=1,3
       DO J=1,3
          CALL KILL(D%S%s(I,J))
          d%e_ij(i,j)=0.0_dp
       ENDDO
    ENDDO

  END    subroutine KILL_DASPIN




  subroutine ALLOC_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S

    CALL ALLOC(S%X,3)
    !     S%G=A_PARTICLE
  END    subroutine ALLOC_SPINOR_8

  subroutine ALLOC_res_SPINOR_8(S)
    implicit none
    TYPE(RES_SPINOR_8), INTENT(INOUT) :: S

    CALL ALLOC(S%X,3)
    !     S%G=A_PARTICLE
  END    subroutine ALLOC_res_SPINOR_8

  subroutine ALLOC_probe_8(R)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    INTEGER I !,J

    !    CALL ALLOC(R%S)
    DO I=1,3
       CALL ALLOC(R%S(I))
    ENDDO
    CALL ALLOC(R%X,6)
    !      R%S(0)%X(N0_NORMAL)=ONE
    DO I=1,3
       R%S(I)=0
    ENDDO
    CALL ALLOC(R%ac)
    r%e_ij=0.0_dp

  END    subroutine ALLOC_probe_8

  subroutine ALLOC_rf_phasor_8(R)
    implicit none
    TYPE(rf_phasor_8), INTENT(INOUT) :: R
    INTEGER I  !,J

    DO I=1,2
       CALL alloc(R%X(I))
    ENDDO
    CALL alloc(R%om)
!    CALL alloc(R%t)

  END    subroutine ALLOC_rf_phasor_8

  subroutine KILL_res_SPINOR_8(S)
    implicit none
    TYPE(RES_SPINOR_8), INTENT(INOUT) :: S

    CALL KILL(S%X,3)
    !     S%G=A_PARTICLE
  END    subroutine KILL_res_SPINOR_8

  subroutine KILL_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S

    CALL KILL(S%X,3)


  END    subroutine KILL_SPINOR_8

  subroutine KILL_probe_8(R)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    INTEGER I  !,J

    !    CALL KILL(R%S)
    DO I=1,3
       CALL KILL(R%S(I))
    ENDDO
    CALL KILL(R%X,6)

    CALL KILL(R%ac)
    r%e_ij=0.0_dp
  END    subroutine KILL_probe_8

  subroutine kill_rf_phasor_8(R)
    implicit none
    TYPE(rf_phasor_8), INTENT(INOUT) :: R
    INTEGER I  !,J

    DO I=1,2
       CALL KILL(R%X(I))
    ENDDO
    CALL KILL(R%om)
!    CALL KILL(R%t)

  END    subroutine kill_rf_phasor_8




!!!!!!!!!!!!!!!!   new stuff
  subroutine get_spin_nx_spinor_8(DS,theta0,n0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    type(real_8), intent(inout) :: theta0
    type(spinor_8), intent(inout) :: n0

    call get_spin_n0(DS,theta0,n0%x)  !get_spin_nx_t

  end subroutine get_spin_nx_spinor_8


  subroutine get_spin_nx_t(DS,theta0,n0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    type(damapspin) a1i,ds0
    type(real_8), intent(inout) :: theta0,n0(3)
    type(real_8) s(3,3),a(3,3),ai(3,3)
    type(real_8) a11,a13

    call alloc(a1i)
    call alloc(ds0)
    call alloc_33(s)
    call alloc_33(a)
    call alloc_33(ai)
    call alloc(a11,a13)

    call find_n0(ds%s%s,n0)

    !     call print(theta0,6)


    call find_a(n0,a)

    call inv_as(a,ai)



    call matmulp(ai,ds%s%s,ai)    !
    call matmulp(ai,a,s)    !

    a11=(s(1,1))
    a13=(s(1,3))

    theta0=clockwise*atan2(a13,a11)
    if(force_positive.and.theta0<0.0_dp) theta0 = theta0 + twopi    !!!! allow negative theta0

    ! at this stage the spin map is:

    !  exp( theta0 n0 L )   where n0 is not the real n but it depends on (x,p)


    !     call print(theta0,6)

    !    theta0=theta0*nsc

    !     call print(theta0,6)

    !     write(6,*) mod((theta0.sub.'0'),twopi)+twopi

    call kill(a1i)
    call kill(ds0)
    call kill_33(s)
    call kill_33(a)
    call kill_33(ai)

    call kill(a11,a13)

  end subroutine get_spin_nx_t

  subroutine get_spin_nx_r(S,theta0,n0)
    implicit none
    real(dp), intent(inout) :: theta0,n0(3)
    real(dp) a(3,3),ai(3,3)
    real(dp),  INTENT(INout) :: S(3,3)

    call find_n0(s,n0)





    call find_a(n0,a)
    call inv_as(a,ai)
    ai=matmul(ai,s)    !
    s=matmul(ai,a)    !

    theta0=clockwise*atan2(s(1,3),s(1,1))
    if(force_positive.and.theta0<0.0_dp)  theta0 = theta0 + twopi   !!!! allow negative theta0


  end subroutine get_spin_nx_r

  subroutine get_spin_nx_probe(xs0,theta0,n0)
    implicit none
    type(probe), intent(in) :: xs0
    real(dp), intent(inout) :: theta0,n0(3)
    integer i,j
    real(dp)  S(3,3)

    do i=1,3
       do j=1,3
          s(i,j)=xs0%s(j)%x(i)
       enddo
    enddo

    call  get_spin_n0(S,theta0,n0)


  end subroutine get_spin_nx_probe


  subroutine get_spin_nx_rd(dS,theta0,n0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    real(dp), intent(inout) :: theta0,n0(3)
    real(dp) a(3,3),ai(3,3),S(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          s(i,j)=ds%s%s(i,j)
       enddo
    enddo

    call find_n0(s,n0)

    call find_a(n0,a)
    call inv_as(a,ai)
    ai=matmul(ai,s)    !
    s=matmul(ai,a)    !

    theta0=clockwise*atan2(s(1,3),s(1,1))
    if(force_positive.and.theta0<0.0_dp)  theta0 = theta0 + twopi   !!!! allow negative theta0



  end subroutine get_spin_nx_rd



!!!!!!!!!!!!!!
  FUNCTION POWMAP( S1, R2 )
    implicit none
    TYPE (damapspin) POWMAP
    TYPE (damapspin), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (damapspin) S11
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    !    call checkdamap(s1)

    call ass(POWMAP)

    call alloc(s11)

    s11=1


    R22=IABS(R2)
    DO I=1,R22
       s11=s1*s11
    ENDDO

    IF(R2.LT.0) THEN
       ! if(old) then

       s11%m=s11%m**(-1)
       call trans_mat(s11%s%s,s11%m,s11%s%s)
       call inv_as(s11%s%s,s11%s%s)

       !      CALL etinv(S11%v%i,S11%v%i)
       !       else
       !          CALL newetinv(S11%v%j,S11%v%j)
       !       endif
    ENDIF

    powmap=s11


    ! powmap=junk
    call kill(s11)

    master=localmaster

  END FUNCTION POWMAP

  FUNCTION dot_spinor( S1, S2 )
    implicit none
    real(dp) dot_spinor
    TYPE (SPINOR), INTENT (IN) :: S1,S2

    INTEGER I

    dot_spinor=0.0_dp

    DO I=1,3
       dot_spinor=dot_spinor+s1%x(i)*s2%x(i)
    ENDDO



  END FUNCTION dot_spinor

  FUNCTION dot_spinor_8( S1, S2 )
    implicit none
    TYPE (real_8) dot_spinor_8
    TYPE (SPINOR_8), INTENT (IN) :: S1,S2

    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    !    call checkdamap(s1)

    call ass(dot_spinor_8)

    dot_spinor_8=0.0_dp

    DO I=1,3
       dot_spinor_8=dot_spinor_8+s1%x(i)*s2%x(i)
    ENDDO


    master=localmaster

  END FUNCTION dot_spinor_8

  FUNCTION concat(S2,S1)
    implicit none
    TYPE (damapspin) concat,t2
    TYPE (damapspin), INTENT (IN) :: S1, S2
    integer i,j,k
    integer localmaster
    type(real_8) s(3,3)


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call ass(concat)
    call alloc(t2);
    call alloc_33(s);

    t2=s2
    concat=s2


    do i=1,3
       do j=1,3
          !     if(s2%s(i,j)%kind==2) then
          !          if(c_%nd2/=0) then
          t2%s%s(i,j)=s2%s%s(i,j)*s1%m
          !          else
          !           t2%s(i,j)=s2%s(i,j)
          !          endif
          !     else
          !      t2%s(i,j)=s2%s(i,j)
          !     endif
       enddo
    enddo

    if(c_%nd2/=0) concat%m=s2%m*s1%m

    !    do i=1,6
    !     concat%x(i)=s2%x(i)
    !    enddo
    do i=1,3
       do j=1,3
          do k=1,3
             s(i,j)= t2%s%s(i,k)*s1%s%s(k,j)+ s(i,j)
          enddo
       enddo
    enddo
    call smatp(1.0_dp,s,concat%s%s)
    !    concat%s0=s1%s0

    call concat_envelope(S2,S1,concat)


    call kill_33(s);
    call kill(t2);
    master=localmaster

  END FUNCTION concat

  subroutine concat_envelope(S2,S1,S3)
    implicit none
    TYPE (damapspin), INTENT (IN) :: S1
    TYPE (damapspin), INTENT (IN) :: S2
    TYPE (damapspin), INTENT (INout) :: s3
    real(dp) s1mi(ndim2,ndim2)
    real(dp) e(ndim2,ndim2)

    s1mi=0.0_dp
    e=0.0_dp

    s1mi=(s1%m.sub.1)**(-1)

    e(1:c_%nd2,1:c_%nd2)=s2%e_ij
    e=matmul(matmul(s1mi,e),transpose(s1mi))

    s3%e_ij=e(1:c_%nd2,1:c_%nd2)+s1%e_ij

  end subroutine concat_envelope

  subroutine extract_envelope_damap(S1,E0_ij,E_ij)
    implicit none
    TYPE (damapspin), INTENT (IN) :: S1
    real(dp), INTENT (in) :: E0_ij(6,6)
    real(dp), INTENT (out) :: E_ij(6,6)
    real(dp) s1m(ndim2,ndim2)
    real(dp) e(ndim2,ndim2)

    s1m=0.0_dp
    e=0.0_dp


    s1m=s1%m.sub.1


    e(1:c_%nd2,1:c_%nd2)=s1%e_ij+E0_ij
    e=matmul(matmul(s1m,e),transpose(s1m))

    E_ij=e(1:c_%nd2,1:c_%nd2)

  end subroutine extract_envelope_damap

  subroutine extract_envelope_probe8(p,E0_ij,E_ij)
    implicit none
    TYPE (probe_8), INTENT (IN) :: p
    TYPE (damapspin)  S1
    real(dp), INTENT (in) :: E0_ij(6,6)
    real(dp), INTENT (out) :: E_ij(6,6)
    real(dp) s1m(ndim2,ndim2)
    real(dp) e(ndim2,ndim2)

    call alloc(s1)
    s1=p

    call extract_envelope_damap(s1,E0_ij,E_ij)

    call kill(s1)
  end subroutine extract_envelope_probe8


  FUNCTION cmul(S2,S1)   ! multiply spin part with real(dp) s1
    implicit none
    TYPE (damapspin) cmul
    TYPE (damapspin), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2
    integer i,j
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call ass(cmul)

    if(c_%nd2/=0) cmul%m=s1%m



    do i=1,3
       do j=1,3
          cmul%s%s(i,j)=s2*s1%s%s(i,j)
       enddo
    enddo
    master=localmaster

  END FUNCTION cmul

  FUNCTION mmul(S1,S2)  ! transform spin part with damap s2
    implicit none
    TYPE (damapspin) mmul
    TYPE (damapspin), INTENT (IN) :: S1
    type(damap), INTENT (IN) :: s2
    !    integer i,j,k
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call ass(mmul)

    mmul=s1
    call trans_mat( mmul%s%s,s2, mmul%s%s)

    !    do i=1,3
    !       do j=1,3
    !     if(s2%s(i,j)%kind==2) then
    !           mmul%s(i,j)=s1%s(i,j)*s2
    !       enddo
    !    enddo

    master=localmaster

  END FUNCTION mmul

  FUNCTION damapspin_spinor8_mul(S1,S2) ! transform spin part with spinor_8 s2 and returns spinor_8
    implicit none
    type(spinor_8) damapspin_spinor8_mul
    TYPE (damapspin), INTENT (IN) :: S1
    type(spinor_8), INTENT (IN) :: S2
    integer i,j
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(damapspin_spinor8_mul%x(1))
    call assp_no_master(damapspin_spinor8_mul%x(2))
    call assp_no_master(damapspin_spinor8_mul%x(3))

    damapspin_spinor8_mul%x(1)=0.0_dp
    damapspin_spinor8_mul%x(2)=0.0_dp
    damapspin_spinor8_mul%x(3)=0.0_dp

    do i=1,3
       do j=1,3
          damapspin_spinor8_mul%x(i)=(s1%s%s(i,j))*S2%x(j)+damapspin_spinor8_mul%x(i)
       enddo
    enddo

    master=localmaster

  END FUNCTION damapspin_spinor8_mul


  FUNCTION damapspin_spinor_mul(S1,S2) ! transform spin part with spinor_8 s2 and returns spinor_8
    implicit none
    type(spinor_8) damapspin_spinor_mul
    TYPE (damapspin), INTENT (IN) :: S1
    type(spinor), INTENT (IN) :: S2
    integer i,j
    integer localmaster



    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(damapspin_spinor_mul%x(1))
    call assp_no_master(damapspin_spinor_mul%x(2))
    call assp_no_master(damapspin_spinor_mul%x(3))

    damapspin_spinor_mul%x(1)=0.0_dp
    damapspin_spinor_mul%x(2)=0.0_dp
    damapspin_spinor_mul%x(3)=0.0_dp

    do i=1,3
       do j=1,3
          damapspin_spinor_mul%x(i)=(s1%s%s(i,j))*S2%x(j)+damapspin_spinor_mul%x(i)
       enddo
    enddo

    master=localmaster


  END FUNCTION damapspin_spinor_mul

  FUNCTION exp_spinor_8(S1)  ! transform spin part with damap s2
    implicit none
    TYPE (damapspin) exp_spinor_8
    TYPE (spinor_8), INTENT (IN) :: S1
    !    integer i,j,k
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call ass(exp_spinor_8)

    call daexplogp(S1,exp_spinor_8)

    master=localmaster

  END FUNCTION exp_spinor_8


  FUNCTION spin8_mul_map(S2,S1)  ! transforms spinor_8 s2 with damap s1
    implicit none
    type(spinor_8) spin8_mul_map
    TYPE (damap), INTENT (IN) :: S1
    type(spinor_8), INTENT (IN) :: S2
    integer i
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(spin8_mul_map%x(1))
    call assp_no_master(spin8_mul_map%x(2))
    call assp_no_master(spin8_mul_map%x(3))

    spin8_mul_map%x(1)=0.0_dp
    spin8_mul_map%x(2)=0.0_dp
    spin8_mul_map%x(3)=0.0_dp

    do i=1,3
       if(S2%x(i)%kind==2) then
          spin8_mul_map%x(i)=S2%x(i)%t*s1
       else
          spin8_mul_map%x(i)=S2%x(i)
       endif
    enddo

    master=localmaster

  END FUNCTION spin8_mul_map

  FUNCTION spin8_scal8_map(S1,S2)  ! transforms spinor_8 s2 with damap s1
    implicit none
    type(spinor_8) spin8_scal8_map
    TYPE (real_8), INTENT (IN) :: S1
    type(spinor_8), INTENT (IN) :: S2
    integer i
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(spin8_scal8_map%x(1))
    call assp_no_master(spin8_scal8_map%x(2))
    call assp_no_master(spin8_scal8_map%x(3))

    spin8_scal8_map%x(1)=0.0_dp
    spin8_scal8_map%x(2)=0.0_dp
    spin8_scal8_map%x(3)=0.0_dp

    do i=1,3
       spin8_scal8_map%x(i)=s1*S2%x(i)
    enddo

    master=localmaster

  END FUNCTION spin8_scal8_map


  FUNCTION add_spin8_spin8(S1,S2)  ! transforms spinor_8 s2 with damap s1
    implicit none
    type(spinor_8) add_spin8_spin8
    type(spinor_8), INTENT (IN) :: S1,S2
    integer i
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(add_spin8_spin8%x(1))
    call assp_no_master(add_spin8_spin8%x(2))
    call assp_no_master(add_spin8_spin8%x(3))

    add_spin8_spin8%x(1)=0.0_dp
    add_spin8_spin8%x(2)=0.0_dp
    add_spin8_spin8%x(3)=0.0_dp

    do i=1,3
       add_spin8_spin8%x(i)=s1%x(i)+S2%x(i)
    enddo

    master=localmaster

  END FUNCTION add_spin8_spin8

  FUNCTION mul_spin8_spin8(S1,S2)  ! transforms spinor_8 s2 with damap s1
    implicit none
    type(spinor_8) mul_spin8_spin8
    type(spinor_8), INTENT (IN) :: S1,S2
    integer i
    integer localmaster

    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(mul_spin8_spin8%x(1))
    call assp_no_master(mul_spin8_spin8%x(2))
    call assp_no_master(mul_spin8_spin8%x(3))

    mul_spin8_spin8%x(1)=0.0_dp
    mul_spin8_spin8%x(2)=0.0_dp
    mul_spin8_spin8%x(3)=0.0_dp

    mul_spin8_spin8%x(1)=s1%x(2)*S2%x(3)-s1%x(3)*S2%x(2)
    mul_spin8_spin8%x(2)=s1%x(3)*S2%x(1)-s1%x(1)*S2%x(3)
    mul_spin8_spin8%x(3)=s1%x(1)*S2%x(2)-s1%x(2)*S2%x(1)

    master=localmaster

  END FUNCTION mul_spin8_spin8

  FUNCTION sub_spin8_spin8(S1,S2)  ! transforms spinor_8 s2 with damap s1
    implicit none
    type(spinor_8) sub_spin8_spin8
    type(spinor_8), INTENT (IN) :: S1,S2
    integer i
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call assp_master(sub_spin8_spin8%x(1))
    call assp_no_master(sub_spin8_spin8%x(2))
    call assp_no_master(sub_spin8_spin8%x(3))

    sub_spin8_spin8%x(1)=0.0_dp
    sub_spin8_spin8%x(2)=0.0_dp
    sub_spin8_spin8%x(3)=0.0_dp

    do i=1,3
       sub_spin8_spin8%x(i)=s1%x(i)-S2%x(i)
    enddo

    master=localmaster

  END FUNCTION sub_spin8_spin8

  function eval_spinor_8(s,x)
    implicit none
    TYPE(spinor) eval_spinor_8
    TYPE(spinor_8), INTENT(IN) :: s
    real(dp), INTENT(IN) :: x(lnv)
    integer i


    do i=1,3
       if(s%x(i)%kind==2) then
          eval_spinor_8%x(i)=s%x(i)%t*x
       endif
    enddo

  end function eval_spinor_8

  FUNCTION concatxp(S2,x)
    implicit none
    TYPE (damapspin) concatxp
    TYPE (damapspin), INTENT (IN) :: S2
    real(dp), INTENT (IN) :: x(lnv)
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call ass(concatxp)

    call eval_spin_matrix(S2,x,concatxp)

    master=localmaster

  END FUNCTION concatxp


  FUNCTION addm(S2,S1)  ! adds spin part of s1 and s2
    implicit none
    TYPE (damapspin) addm
    TYPE (damapspin), INTENT (IN) :: S1, S2
    integer i,j
    integer localmaster


    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call ass(addm)

    addm=s2

    do i=1,3
       do j=1,3
          addm%s%s(i,j)=s2%s%s(i,j)+s1%s%s(i,j)
       enddo
    enddo


    master=localmaster

  END FUNCTION addm


  real(dp) function purge_transverse(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    if(.not.c_%stable_da) return

    purge_transverse=1.0_dp


    do i=1,c_%nd2

       if(i/=c_%ndpt) then
          if(j(i)/=0) then
             purge_transverse=0.0_dp
             exit
          endif
       endif
    enddo

    return
  end function purge_transverse

  subroutine clean_orbital_33(s,sf)
    implicit none
    TYPE (real_8), intent(inout) :: s(3,3),sf(3,3)
    integer i,j
    type(taylor) t
    logical doflip
    call alloc(t)

    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       perform_flip=.false.
       do i=1,3
          do j=1,3
             if(s(i,j)%kind==2) call fliptaylor(s(i,j)%t,s(i,j)%t,1)
          enddo
       enddo
       doflip=.true.
    else
       doflip=.false.
    endif

    do i=1,3
       do j=1,3
          if(s(i,j)%kind==2) then
             call cfu(s(i,j)%t,purge_transverse,t)
             sf(i,j)=t
          else
             sf(i,j)=s(i,j)
          endif
       enddo
    enddo
    call kill(t)

    if(doflip) then
       do i=1,3
          do j=1,3
             if(s(i,j)%kind==2) call fliptaylor(s(i,j)%t,s(i,j)%t,-1)
             if(s(i,j)%kind==2.and.s(i,j)%t%i/=sf(i,j)%t%i) then
                call fliptaylor(sf(i,j)%t,sf(i,j)%t,-1)
             endif
          enddo
       enddo
       perform_flip=.true.
    endif


  end  subroutine clean_orbital_33

!!!!!!!!!!  Normal form !!!!!!!!!!!!!!!

  subroutine alloc_normal_spin(D)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: D

    CALL alloc(D%n0,3)
    CALL alloc(D%theta0)
    CALL alloc(D%n)
    CALL alloc(D%a_t)
    CALL alloc(D%as)
    CALL alloc(D%a1)
    CALL alloc(D%ar)
    D%NRES=0
    D%M=0
    D%Ms=0
    D%s_ij0=0.0_dp
    D%s_ijr=0.0_dp
    D%emittance=0.0_dp
    D%tune=0.0_dp
    D%damping=0.0_dp
    D%AUTO=my_true
    D%STOCHASTIC=my_false
    D%STOCH=0.0_dp
    D%STOCH_inv=0.0_dp
    D%nu=0.0_dp

  END    subroutine alloc_normal_spin

  subroutine KILL_normal_spin(D)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: D

    CALL KILL(D%n0,3)
    CALL KILL(D%theta0)
    CALL KILL(D%n)
    CALL KILL(D%a_t)
    CALL KILL(D%as)
    CALL KILL(D%a1)
    CALL KILL(D%ar)

  END    subroutine KILL_normal_spin

!!!!!!!!!!!!!!!!   new stuff
  subroutine Go_to_closed(ns,DS,spin_in)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: ns
    TYPE(damapspin), INTENT(INout) :: DS
    logical(lp), INTENT(IN) :: spin_in
    logical(lp) rad_in
    type(damapspin) a1i,ds0
    type(damapspin)s ,a ,ai
    type(taylor) nn
    type(real_8) a11,a13
    integer i,j,jj(lnv)
    real(dp) ss(3,3)
!!!!!!!!!!   wrapping for radiation  !!!!!!!!!!!!!!!!!!!
    type(radtaylor) ys(ndim2)
    type(beamenvelope) env


    call alloc(a1i)
    call alloc(ds0)
    call alloc(s)
    call alloc(a)
    call alloc(ai)
    call alloc(nn)
    call alloc(a11,a13)

    ns%as=1   ! damapspin to fix point

    ns%n=ds%m   ! ORBITAL NORMALIZATION

    ns%as%m=ns%n%a1    ! fix point map a1 of orbital normal form

    a1i=ns%as**(-1) ! fix point map a1 of orbital normal form


    DS0=a1i*ds*ns%as
    !    DS0=a1i*ds*ns%as  ! at this stage ds0 is the map around the fixed point
    ! but spin matrix unchanged except for its dependence on transverse

!!!!!!!!!!   wrapping for radiation  !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!   later to be change changed  !!!!!!!!!!!!!!!
    call check_rad(ds%e_ij,rad_in)
    if(rad_in) then
       if(c_%no==1) then
          call normalize_envelope(ns,DS)
          Write(6,*) "New envelope calculation attempted with NO=1 "
       else
          call alloc(ys,6)
          call alloc(env)
          jj=0
          do i=1,6
             ys(i)%v=ds%m%v(i)
             do j=1,6
                ys(i)%e(j)=ds%e_ij(i,j)
             enddo
          enddo
          env%stochastic=ns%stochastic
          env%auto=ns%auto
          env=ys

          ns%s_ij0      =  env%s_ij0
          ns%emittance  =  env%emittance
          ns%KICK       =  env%KICK
          do i=1,6
             do j=1,6
                jj(j)=1
                ns%STOCH(i,j)  =  env%STOCH%v(i).sub.jj
                jj(j)=0
             enddo
          enddo

          call kill(env)
          call kill(ys,6)
       endif
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    s=1
    a=1
    ai=1

    call clean_orbital_33(ds0%s%s,s%s%s)
    ! at this stage s is the spin map  without dependence on orbital

    if(spin_in) then           !!!!  spin in
       call find_n0(s%s%s,ns%n0)

       call find_a(ns%n0,a%s%s)

       ai=a**(-1)

       s=(ai*s)*a

       !write(6,*) " diagonal "
       !    do i=1,3
       !    do j=1,3
       !    ss(i,j)=s%s(i,j)
       !    enddo
       !    write(6,*) ss(i,1:3)
       !    enddo

       !pause 777

       a11=s%s%s(1,1)
       a13=s%s%s(1,3)

       ns%theta0=clockwise*atan2(a13,a11)
       ns%nu=ns%theta0/twopi
       if(force_positive.and.ns%theta0<0.0_dp)  ns%theta0 = ns%theta0 + twopi  !!!! allow negative theta0

       ns%as=ns%as*a

       ns%ar=1
       ns%ar%m=ns%as%m**(-1)*ns%n%a_t
    else
       ns%ar=1
       ns%ar%m=ns%as%m**(-1)*ns%n%a_t

       ns%theta0=0.0_dp  ! arbitrary junk
       ns%n0(1)=0.0_dp   ! arbitrary junk
       ns%n0(2)=1.0_dp    ! arbitrary junk
       ns%n0(3)=0.0_dp   ! arbitrary junk
    endif

    call kill(a1i)
    call kill(ds0)
    call kill(s)
    call kill(a)
    call kill(ai)

    call kill(a11,a13)
    call kill(nn)

  end subroutine Go_to_closed

  subroutine normalize_envelope(norm_spin,m_spin)
    implicit none
    type(normal_spin) norm_spin
    type(damapspin) m_spin
    type(normalform) norm
    type(damap) m
    integer i,j,i1,i2
    real(dp) a(6,6),ai(6,6),ait(6,6),mat(6,6), sigma_inf(6,6),at(6,6),br(6,6)
    complex(dp) c(6,6),ci(6,6),cit(6,6), b(6,6),ct(6,6)
    complex(dp) coef,ba(6,6),b_phasor(6,6)
    complex(dp) R(6,6),r_phasor(6,6), sigma_inf_phasor(6,6)
    real(dp) xj(6,6),mj(6,6),xn,jb(6,6),bs(6,6)

    b=m_spin%e_ij
    bs=m_spin%e_ij

    c=0.0_dp
    ci=0.0_dp
    do i=1,3
       do j=1,3
          xj(2*i,2*i-1)=-1.0_dp
          xj(2*i-1,2*i)=1.0_dp
          c(2*i-1,2*i-1)=0.5_dp
          c(2*i-1,2*i)=0.5_dp
          c(2*i,2*i-1)=0.5_dp/i_
          c(2*i,2*i)=-0.5_dp/i_
          ci(2*i-1,2*i-1)=1.0_dp
          ci(2*i-1,2*i)=i_
          ci(2*i,2*i-1)=1.0_dp
          ci(2*i,2*i)=-i_
       enddo
    enddo

    mat=m_spin%m
    a=norm_spin%n%a_t
    ai=norm_spin%n%a_t**(-1)
    ait=transpose(ai)
    at=transpose(a)
    cit=transpose(ci)
    ct=transpose(c)

    R=matmul(matmul(ai,b),ait)

    ba=matmul(matmul(ai,b),ait)

    b_phasor=matmul(matmul(ci,ba),cit)

    r=matmul(matmul(ai,mat),a)
    r_phasor=matmul(matmul(ci,r),c)

    do i=1,6
       do j=1,6
          sigma_inf_phasor(i,j)= r_phasor(i,i)*r_phasor(j,j)/(1.0_dp-r_phasor(i,i)*r_phasor(j,j))*b_phasor(i,j)
       enddo
    enddo
    do i=1,3
       norm_spin%emittance(i)=sigma_inf_phasor(2*i-1,2*i)/2.0_dp
    enddo

    sigma_inf=matmul(matmul(c,sigma_inf_phasor),ct)
    sigma_inf=matmul(matmul(a,sigma_inf),at)

    norm_spin%s_ij0=sigma_inf
    norm_spin%s_ijr=sigma_inf_phasor

    norm_spin%tune=norm_spin%n%tune(1:3)
    norm_spin%damping=norm_spin%n%damping(1:3)



    if(norm_spin%STOCHASTIC) then
       call diagonalise_envelope_a(bs,br,a,ai,norm_spin%kick)   !diagonalise_envelope_a(b,br,a,ai,kick)
       norm_spin%STOCH=a
       norm_spin%STOCH_inv=ai
    endif

  end subroutine normalize_envelope

  subroutine fetch_s0(DS,s)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(damapspin), INTENT(INout) :: s
    TYPE(damapspin) s0
    integer i,j

    call alloc(s0)

    s0=1

    do i=1,3
       do j=1,3
          s0%s%s(i,j)=ds%s%s(i,j).sub.'0'
       enddo
    enddo
    s=s0
    call kill(s0)

  end subroutine fetch_s0


  subroutine dalog_spinor_8(DS,n)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(spinor_8), INTENT(INout) :: n
    TYPE(taylor) om(3)
    integer i
    call alloc(om)

    call dalog(DS,om)

    do i=1,3
       n%x(i)=morph(om(i))
    enddo

    call kill(om)
  end subroutine dalog_spinor_8

  subroutine factor_parameter_dependent_s0(DS,s0,NS,N_AXIS,DIR)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS,s0,NS
    TYPE(damapspin)A_F,A_S
    TYPE(spinor_8), INTENT(INout) :: N_AXIS
    INTEGER DIR

    CALL ALLOC(A_F,A_S)

    A_f=DS
    a_f%m=1
    CALL clean_orbital_33(A_f%S%s,A_f%S%s)
    A_S=DS
    a_s%m=1
    IF(DIR==1) THEN
       a_S=a_f**(-1)*a_s   !   =s0*ns
    ELSE
       a_S=a_s*a_f**(-1)   !   =ns*s0    of course s0 is different
    ENDIF
    CALL dalog_spinor_8(a_S,N_AXIS)   !  S=exp(N_AXIS)
    S0=A_F
    NS=A_S


    CALL KILL(A_F,A_S)

  END subroutine factor_parameter_dependent_s0


  subroutine factor_s0(DS,s0,NS,N_AXIS,DIR)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS,s0,NS
    TYPE(damapspin)A_F,A_S
    TYPE(spinor_8), INTENT(INout) :: N_AXIS
    INTEGER DIR

    CALL ALLOC(A_F,A_S)

    A_f=DS
    a_f%m=1
    CALL fetch_s0(A_f,A_f)
    A_S=DS
    a_s%m=1
    IF(DIR==1) THEN
       a_S=a_f*a_s**(-1)   ! angle has no constant part (DA) BUT PARAMETERS
    ELSE
       a_S=a_f**(-1)*a_s   ! angle has no constant part (DA)
    ENDIF
    CALL dalog_spinor_8(a_S,N_AXIS)   !  S=exp(N_AXIS)
    S0=A_F
    NS=A_S


    CALL KILL(A_F,A_S)

  END subroutine factor_s0


  subroutine dalog(DS,om)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(taylor), INTENT(INout) :: om(3)
    TYPE(damapspin) h
    integer i
    TYPE(damapspin) dh,dhn
    real(dp) c,cl
    call alloc(h)
    call alloc(dh)
    call alloc(dhn)
    !  this only works with a da-map

    dh=ds
    dh%m=1
    h=0
    dhn=1
    dh%s%s(1,1)=dh%s%s(1,1)-1.0_dp
    dh%s%s(2,2)=dh%s%s(2,2)-1.0_dp
    dh%s%s(3,3)=dh%s%s(3,3)-1.0_dp

    c=1.0_dp
    do i=1,c_%no
       dhn=dhn*dh
       cl=c/i
       h=h+cl*dhn
       c=-c
    enddo

    om(1)=h%s%s(3,2)
    om(2)=h%s%s(1,3)
    om(3)=h%s%s(2,1)

    call kill(h)
    call kill(dh)
    call kill(dhn)

  end subroutine dalog

  subroutine daexplog(om,ds)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(taylor), INTENT(INout) :: om(3)
    integer i
    TYPE(damapspin) dh,dhn
    real(dp) c
    call alloc(dh)
    call alloc(dhn)
    !  this only works with a da-map
    ds=1
    dh%m=1
    dh%s%s(2,1)=om(3)
    dh%s%s(1,3)=om(2)
    dh%s%s(3,2)=om(1)
    dh%s%s(1,2)=-om(3)
    dh%s%s(3,1)=-om(2)
    dh%s%s(2,3)=-om(1)
    dhn=1
    c=1.0_dp
    do i=1,c_%no
       dhn=dhn*dh
       c=c/i
       ds=ds+c*dhn
    enddo


    call kill(dh)
    call kill(dhn)

  end subroutine daexplog

  subroutine get_kernel(ns,om,oma)
    implicit none
    TYPE(normal_spin), INTENT(INout) :: ns
    TYPE(taylor), INTENT(INout) :: om(3),oma(3)
    type(taylorresonance) t
    TYPE(complextaylor) omr(3),omc(3)
    integer N,i,nd,j
    integer, allocatable :: jc(:)
    real(dp) value,ang,tune(4)
    complex(dp) denom
    logical doit,doflip

    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       perform_flip=.false.
       call fliptaylor(om(1),om(1),1)
       call fliptaylor(om(2),om(2),1)
       call fliptaylor(om(3),om(3),1)
       call flip_real_array(ns%n%tune,ns%n%tune,1)
       if(use_ptc_ac_position) then
          call flip_resonance(ns%m,ns%m,1)
       endif
       doflip=.true.
    else
       doflip=.false.
    endif


    nd=c_%nd2/2
    if(c_%ndpt/=0) nd=nd-1

    call alloc(omc,3)
    call alloc(omr,3)
    call alloc(t)
    !    write(6,*) c_%nd2,c_%ndpt,nd
    !    write(6,*) ns%n%tune


    omc(1)= om(1)-i_*om(3)  ! coeff of lambdda
    omc(3)= om(1)+i_*om(3)  ! coeff of lambdda*
    omc(2)= om(2)

    allocate(jc(c_%nv))


    !   om(2)    is always real!!!
    t=omc(2)%r
    !   cos part of omc(2) !!!

    call taylor_cycle(t%cos,size=n)

    do i=1,n
       call taylor_cycle(t%cos,ii=i,value=value,j=jc)

       call test_jc(ns,jc,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-1.0_dp)
          omr(2)=omr(2)+ (denom.mono.jc)
       endif
    enddo

    !  sin  part of om(2)
    call taylor_cycle(t%sin,size=N)

    do i=1,n
       call taylor_cycle(t%sin,ii=i,value=value,j=jc)

       call test_jc(ns,jc,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-1.0_dp)
          omr(2)=omr(2)+ i_*(denom.mono.jc)
       endif
    enddo
    doit=.true.
    !  real part of om(1)
    t=omc(1)%r
    call taylor_cycle(t%cos,N)

    do i=1,n
       call taylor_cycle(t%cos,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(clockwise*i_*ns%theta0))
          omr(1)=omr(1)+ (denom.mono.jc)
       endif
    enddo

    call taylor_cycle(t%SIN,N)

    do i=1,n
       call taylor_cycle(t%SIN,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(clockwise*i_*ns%theta0))
          omr(1)=omr(1)+ I_*(denom.mono.jc)
       endif
    enddo

    !  imaginary part of om(1)
    t=omc(1)%i
    call taylor_cycle(t%cos,N)

    do i=1,n
       call taylor_cycle(t%cos,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(clockwise*i_*ns%theta0))
          omr(1)=omr(1)+ i_*(denom.mono.jc)
       endif
    enddo

    call taylor_cycle(t%SIN,N)

    do i=1,n
       call taylor_cycle(t%SIN,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(clockwise*i_*ns%theta0))
          omr(1)=omr(1)-(denom.mono.jc)
       endif
    enddo

    !  real part of om(3)
    t=omc(3)%r
    call taylor_cycle(t%cos,N)

    do i=1,n
       call taylor_cycle(t%cos,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,-clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(-clockwise*i_*ns%theta0))
          omr(3)=omr(3)+ (denom.mono.jc)
       endif
    enddo

    call taylor_cycle(t%SIN,N)

    do i=1,n
       call taylor_cycle(t%SIN,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,-clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(-clockwise*i_*ns%theta0))
          omr(3)=omr(3)+ I_*(denom.mono.jc)
       endif
    enddo

    !  imaginary part of om(3)
    t=omc(3)%i
    call taylor_cycle(t%cos,N)

    do i=1,n
       call taylor_cycle(t%cos,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,-clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(-clockwise*i_*ns%theta0))
          omr(3)=omr(3)+ i_*(denom.mono.jc)
       endif
    enddo

    call taylor_cycle(t%SIN,N)

    do i=1,n
       call taylor_cycle(t%SIN,ii=i,value=value,j=jc)

       call test_jc_spin(ns,jc,-clockwise,nd,doit)

       if(doit) then
          ang=0.0_dp
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-exp(-clockwise*i_*ns%theta0))
          omr(3)=omr(3)-(denom.mono.jc)
       endif
    enddo
1111 continue
    call kill(omc,3)
    call alloc(omc,3)

    omc(1)= (omr(3)+omr(1))/2.0_dp
    omc(3)= (omr(3)-omr(1))/2.0_dp/i_
    omc(2)=omr(2)

    do i=1,3
       t%cos=omc(i)%r
       t%sin=omc(i)%i
       oma(i)=t
    enddo

    if(doflip) then
       call flip_real_array(ns%n%tune,ns%n%tune,-1)
       call fliptaylor(om(1),om(1),-1)
       call fliptaylor(om(2),om(2),-1)
       call fliptaylor(om(3),om(3),-1)
       call fliptaylor(oma(1),oma(1),-1)
       call fliptaylor(oma(2),oma(2),-1)
       call fliptaylor(oma(3),oma(3),-1)
       if(use_ptc_ac_position) then
          call flip_resonance(ns%m,ns%m,-1)
       endif
       perform_flip=.true.
    endif

    deallocate(jc)
    call kill(t)
    call kill(omr,3)
    call kill(omc,3)

  end subroutine get_kernel

  subroutine test_jc(ns,jc,nd,doit)
    implicit none
    logical doit
    integer i,nd,k,l,j
    integer  jc(:)
    TYPE(normal_spin), INTENT(INout) :: ns

    doit=.true.

    k=0
    do i=1,nd
       k=iabs(jc(i*2-1)-jc(i*2))+k
    enddo
    if(k==0) doit=.false.

    if(.not.doit) return

    do j=1,ns%n%nres
       k=0
       l=0
       do i=1,nd
          k=iabs(jc(i*2-1)-jc(i*2)-ns%n%m(i,j))+k
          l=iabs(jc(i*2-1)-jc(i*2)+ns%n%m(i,j))+l
       enddo
       if(k==0.or.l==0) then
          doit=.false.
          exit
       endif
    enddo


  end subroutine test_jc

  subroutine test_jc_spin(ns,jc,is,nd,doit)
    implicit none
    logical doit
    integer i,nd,k,l,j,is
    integer  jc(:)
    TYPE(normal_spin), INTENT(INout) :: ns

    doit=.true.

    !    write(6,*) jc(1:4),is
    !    write(6,*) ns%m(1,1),ns%m(2,1),ns%ms(1)

    do j=1,ns%nres
       k=0
       l=0
       do i=1,nd
          k=iabs(jc(i*2-1)-jc(i*2)-ns%m(i,j))+k
          l=iabs(jc(i*2-1)-jc(i*2)+ns%m(i,j))+l
       enddo
       k=k+iabs(is-ns%ms(j))
       l=l+iabs(is+ns%ms(j))
       if(k==0.or.l==0) then
          doit=.false.
          exit
       endif
    enddo


  end subroutine test_jc_spin

  subroutine check_spin(DS,spin_in)
    implicit none
    logical(lp), INTENT(INout) :: spin_in
    TYPE(damapspin), INTENT(IN) :: DS
    integer i,j
    real(dp) norm

    spin_in=.true.
    norm=0.0_dp
    do i=1,3
       do j=1,3
          norm=norm+full_abs(DS%s%s(i,j))
       enddo
    enddo
    norm=abs(norm-3.0_dp)

    if(norm<=eps_tpsalie) then
       write(6,*) " Spin Map is identity : not normalized "
       spin_in=.false.
    endif

  end  subroutine check_spin

  subroutine check_rad(e_ij,rad_in)
    implicit none
    logical(lp), INTENT(INout) :: rad_in
    real(dp) e_ij(6,6)
    integer i,j
    real(dp) norm

    rad_in=.true.
    norm=0.0_dp
    do i=1,6
       do j=1,6
          norm=norm+abs(e_ij(i,j))
       enddo
    enddo

    if(norm==0.0_dp) then
       if(global_verbose) write(6,*) " Radiation Envelope  is 0.0_dp : not printed "
       rad_in=.false.
    endif

  end  subroutine check_rad

  subroutine normalise_spin(ns,DS_in)
    implicit none
    TYPE(normal_spin), INTENT(INout) :: ns
    TYPE(damapspin), INTENT(IN) :: DS_in
    TYPE(damapspin)  ds,ds0,dst,dsi
    TYPE(damap)  r0
    type(taylor) om(3),oma(3)
    integer i
    logical(lp) spin_in

    call alloc(ds)
    call alloc(dst)
    call alloc(dsi)
    call alloc(ds0)
    call alloc(r0)
    call alloc(om)
    call alloc(oma)
    spin_in=.false.

    call check_spin(DS_in,spin_in)
    ds=ds_in
    ! step 1
    !  Normalize to the parameter dependent n0 of the theory
    call Go_to_closed(ns,DS,spin_in)

    ! step 2
    ! Apply the transformation found in step 1 to ds
    ds=ns%as**(-1)*ds*ns%as
    ! Apply the transformation found in step 1 to ds
    ds=ns%ar**(-1)*ds*ns%ar   ! around fixed point and transversely normalised

    call fetch_s0(DS,ds0)   ! ds0=(I,exp(theta0 L_y) )


    ns%a_t=1
    r0=ns%n%normal
    if(spin_in) then

       do i=1,c_%no+1

          dst=ds*ds0**(-1)    ! ds=exp(small)*ds0
          ! dst=exp(small)

          call dalog(DSt,om)  ! small=log(dst)  finite # of steps

          call get_kernel(ns,om,oma)

          call daexplog(oma,dsi)     ! dsi=exp(oma(i)*L_i)
          ns%a_t=ns%a_t*dsi          ! updated a_t= a_t*dsi

          ds=ds*dsi                  ! ds=dsi**(-1) ds * dsi
          dsi=dsi**(-1)              ! eventually ds is normalized
          ds=dsi*ds
       enddo                       ! keep going

    endif
    ns%a_t=ns%as*ns%ar*ns%a_t


    dsi=ns%a_t

    ns%a1=1
    ns%a1%m=ns%n%a1

    dsi=ns%a1**(-1)*ns%a_t
    ns%as=dsi
    ns%as%m=1
    ns%as=ns%as*dsi%m**(-1)
    ns%ar=1
    ns%ar%m=dsi%m       ! a_t= a1 o A_s o A_r
    call kill(ds)
    call kill(dst)
    call kill(dsi)
    call kill(ds0)
    call kill(r0)
    call kill(om)
    call kill(oma)
  end subroutine normalise_spin

  subroutine eval_spin_matrix(a_in,x,a_out)
    implicit none
    TYPE(damapspin), INTENT(IN) :: a_in
    TYPE(damapspin), INTENT(out) :: a_out
    TYPE(damapspin)  a_temp
    real(dp) x(lnv)
    integer i,j

    call alloc(a_temp)

    do i=c_%nv+1,lnv
       x(i)=0.0_dp
    enddo

    a_temp=a_in

    do i=1,3
       do j=1,3
          if(a_in%s%s(i,j)%kind==2) then
             a_temp%s%s(i,j)=a_in%s%s(i,j)%t*x
          endif
       enddo
    enddo
    a_out=a_temp

    call kill(a_temp)

  end subroutine eval_spin_matrix


  subroutine factor_as(a_t,a_f,a_s,a_l,a_nl,DR,R_TE,CS_TE,COSLIKE,s0,s_nl)
    implicit none
    TYPE(damapspin), INTENT(INout) :: a_t,a_f,a_s,a_l,a_nl
    logical(lp) lagrange0,factor_spin0
    TYPE(damapspin), optional, intent(inout) ::R_TE,CS_TE,DR,s0,s_nl
    logical(lp) , optional, intent(inout)  :: COSLIKE
    TYPE(damapspin) rot_y,temp,tempi

    factor_spin0=my_false
    if(present(s0).and.present(s_nl)) then
       factor_spin0=my_true
    elseif(.not.(present(s0)).and.(.not.present(s_nl))) then
       ! do nothing
    else
       write(6,*) " error in factor_as "
       stop 333
    endif
    a_f=1
    a_l=1
    a_nl=1
    lagrange0=my_false
    if(present(dr))  then
       lagrange0=my_true
       dr=1
    endif
    if(present(R_TE).and.present(CS_TE)) then
       R_TE=1
       CS_TE=1
    elseif((.not.present(R_TE)).and.(.not.present(CS_TE))) then
       ! nothing to be done
    else
       write(6,*) " factor_as has an error: R_TE and CS_TE must be both present or absent "
       stop 1068

    endif

    call factor(a_t%m,a_f%m,a_l%m,a_nl%m,DR%m,R_TE%m,CS_TE%m,COSLIKE)

    if(lagrange0) then
       a_s=1
       a_s%s=a_t%s
       a_s=a_s*dr**(-1)

       call alloc(rot_y,temp,tempi)
       if(factor_spin0) then
          call remove_y_rot0(a_s,s0,s_nl,r_y=rot_y)
       else
          call remove_y_rot(a_s,r_y=rot_y)
       endif
       if(present(dr)) dr=dr*rot_y
       call kill(rot_y,temp,tempi)

       a_t%s=a_s%s
    endif



    a_s=a_f**(-1)*a_t
    a_s%m=1
    a_s=a_s*a_t%m**(-1)




!!! this creates

    !! (a_t%m,a_t%s) = (a_f%m, I ) o (I ,a_s%s) o (a_l%m,I) o (a_nl%m,I)

    !
    !if(present(a_s0)) then
    !a_s0=1
    !call clean_orbital_33(a_s%s,a_s0%s)

    !a_s=a_s0**(-1)*a_s
    !endif



  end subroutine factor_as

  subroutine CANONIZE( a_t,A_cs,PHASE_ADVANCE,R_TE,CS_TE,COSLIKE )
    implicit none
    TYPE(damap), INTENT(INout) :: a_t,A_cs
    type(taylor),optional, INTENT(INout) ::  PHASE_ADVANCE(:)
    TYPE(damap), optional, intent(inout) ::R_TE,CS_TE 
    logical(lp) , optional, intent(inout)  :: COSLIKE
    TYPE(damap) a_f,a_l,a_nl,dr1,a_tt
    type(onelieexponent) uno
    logical(lp) doflip
    integer i


    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       !    write(6,*) " flipping ",c_%ndpt,c_%nd2-1
       !    pause 7123

       perform_flip=.false.
       call flip_damap(a_t,a_t)
       doflip=.true.
    else
       doflip=.false.
    endif

    call alloc(a_f,a_l,a_nl,dr1,a_tt)
    call alloc(uno)
    dr1=1
    a_tt=a_t
    call factor(a_tt,a_f,a_l,a_nl,DR1,R_TE,CS_TE,COSLIKE)

    A_cs=a_f*a_l*a_nl
    if(present(PHASE_ADVANCE)) then
    uno=dr1

     if(c_%ndpt==0) then
      do i=1,c_%nd
        PHASE_ADVANCE(i)=PHASE_ADVANCE(i)+((uno%VECTOR%v(2*i-1)).k.(2*i))/twopi
      enddo
     else
       if(c_%ndpt>c_%nd2-2) then
        do i=1,c_%nd-1
         PHASE_ADVANCE(i)=PHASE_ADVANCE(i)+((uno%VECTOR%v(2*i-1)).k.(2*i))/twopi
        enddo
       else
        do i=1,c_%nd-2
         PHASE_ADVANCE(i)=PHASE_ADVANCE(i)+((uno%VECTOR%v(2*i-1)).k.(2*i))/twopi
        enddo
         i=c_%nd-1
         PHASE_ADVANCE(i+1)=PHASE_ADVANCE(i)+((uno%VECTOR%v(2*i-1)).k.(2*i))/twopi
       endif
     endif
      
    endif
    

    if(doflip) then
       call flip_damap(a_t,a_t)
       call flip_damap(a_cs,a_cs)
       call flip_damap(dr1,dr1)
     if(present(PHASE_ADVANCE)) then
       do i=1,c_%nd
          call flip_taylor(PHASE_ADVANCE(i),PHASE_ADVANCE(i),-1)
       enddo
     endif  
       perform_flip=.true.
    endif
    call kill(uno)
    call kill(a_f,a_l,a_nl,dr1,a_tt)

  end subroutine CANONIZE

  subroutine factor_am(a_t,a_f,a_l,a_nl,DR,R_TE,CS_TE,COSLIKE)
    implicit none
    TYPE(damap), INTENT(INout) :: a_t,a_nl,a_l,a_f
    integer i,n,k,nt,j
    integer, allocatable :: jc(:)
    real(dp) value,alpha0,sip
    logical doit
    TYPE(damap) atemp,s1,s1i,m1
    TYPE(taylor)a12,a11,p(ndim)
    logical(lp) doflip,dote,lagrange0
    TYPE(damap), optional, intent(inout) ::R_TE,CS_TE,DR
    logical(lp) , optional, intent(inout)  :: COSLIKE
    type(taylor) m(ndim2,ndim2)
    type(taylor) at(2,2),bt(2,2),ct(2,2),dt(2,2),ati(2,2),bti(2,2),alpha,det
    type(gmap) g
    type(onelieexponent) un
    type(reversedragtfinn) rdf
    type(vecresonance) vr
    logical(lp) t_e
    
    t_e=my_true
    lagrange0=my_false
    if(present(dr)) lagrange0=my_true

    call alloc(atemp,s1,s1i,m1)
    call alloc(a12,a11)
    call alloc(p,ndim)
    call alloc_nn(m)
    call alloc_nn(at)
    call alloc_nn(bt)
    call alloc_nn(ct)
    call alloc_nn(dt)
    call alloc_nn(ati)
    call alloc_nn(bti)
    call alloc(alpha,det)
    allocate(jc(c_%nv))



    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       !    write(6,*) " flipping ",c_%ndpt,c_%nd2-1
       !    pause 7123

       perform_flip=.false.
       call flip_damap(a_t,a_t)
       doflip=.true.
    else
       doflip=.false.
    endif

    atemp=1
    do k=1,c_%nd2
       call taylor_cycle(a_t%v(k),N)

       do i=1,n
          call taylor_cycle(a_t%v(k),ii=i,value=value,j=jc)
          call check_fix(jc,0,doit)
          if(doit) atemp%v(k)= atemp%v(k) + (value.mono.jc)
       enddo
    enddo

    !   in case we have ndpt/=0

    if(c_%ndpt==c_%nd2-1) then   ! ptc convention
       atemp%v(c_%nd2-1)=atemp%v(c_%nd2-1)-(1.0_dp.mono.(c_%nd2-1))
       atemp%v(c_%nd2)=atemp%v(c_%nd2)-(1.0_dp.mono.(c_%nd2))
       do i=1,c_%nd-1
          atemp%v(c_%ndpt+1)=atemp%v(c_%ndpt+1) &
               +(atemp%v(2*i).d.c_%ndpt)*(1.0_dp.mono.2*i-1)-(atemp%v(2*i-1).d.c_%ndpt)*(1.0_dp.mono.2*i)
       enddo
    elseif(c_%ndpt==c_%nd2) then      ! Marylie convention
       atemp%v(c_%nd2-1)=atemp%v(c_%nd2-1)-(1.0_dp.mono.(c_%nd2-1))
       atemp%v(c_%nd2)=atemp%v(c_%nd2)-(1.0_dp.mono.(c_%nd2))
       do i=1,c_%nd-1
          atemp%v(c_%ndpt-1)=atemp%v(c_%ndpt-1) &
               -(atemp%v(2*i).d.c_%ndpt)*(1.0_dp.mono.2*i-1)+(atemp%v(2*i-1).d.c_%ndpt)*(1.0_dp.mono.2*i)
       enddo
    endif


    a_f=atemp   !!!  not exact ! temporal part could be wrong!!!

    a_l=a_f**(-1)*a_t

    atemp=0
    do k=1,c_%nd2
       call taylor_cycle(a_l%v(k),N)

       do i=1,n
          call taylor_cycle(a_l%v(k),ii=i,value=value,j=jc)
          call check_fix(jc,1,doit)
          if(doit) atemp%v(k)= atemp%v(k) + (value.mono.jc)
       enddo
    enddo


    !   in case we have ndpt/=0

    if(c_%ndpt/=0) then   ! ptc convention
       atemp%v(c_%nd2-1)=(1.0_dp.mono.(c_%nd2-1))
       atemp%v(c_%nd2)= (1.0_dp.mono.(c_%nd2))

       if(c_%ndpt==c_%nd2) then
          k=c_%nd2-1
       else
          k=c_%nd2
       endif

       call taylor_cycle(a_l%v(k),N)

       do i=1,n
          call taylor_cycle(a_l%v(k),ii=i,value=value,j=jc)
          call check_fix(jc,2,doit)
          if(doit) atemp%v(k)= atemp%v(k) + (value.mono.jc)
       enddo

    endif

    a_nl=atemp**(-1)*a_l
    a_l=atemp

!!! a_t=a_f*a_l*a_nl at this stage

    !   if(present(lagrange)) then
    if(lagrange0) then
       !   enforce Teng-Edwards to all orders in  parameters i.e. A_12=0 and A_3_4=0 etc....
       nt=c_%nd
       if(C_%ndpt/=0) nt=nt-1

       s1=1
       s1i=1

       if(C_%ndpt/=0) then
          s1%v(C_%ndpt)=1.0_dp.mono.C_%ndpt
          s1i%v(C_%ndpt)=1.0_dp.mono.C_%ndpt
       endif

       do i=1,nt
          a11=a_l%v(2*i-1).d.(2*i-1)
          a12=a_l%v(2*i-1).d.(2*i)
          p(i)=-a12/a11
          p(i)=atan(p(i))
          sip=a11*cos(p(i))-a12*sin(p(i))
          if(sip<0) p(i)=p(i)+pi
          if(C_%ndpt/=0) then
             if(mod(C_%ndpt,2)==1) then
                s1%v(C_%ndpt+1)=s1%v(C_%ndpt+1)-(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
                s1i%v(C_%ndpt+1)=s1i%v(C_%ndpt+1)+(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
             else
                s1%v(C_%ndpt-1)=s1%v(C_%ndpt-1)+(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
                s1i%v(C_%ndpt-1)=s1i%v(C_%ndpt-1)-(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
             endif
          endif
          s1%v(2*i-1)=COS(p(i))*(1.0_dp.mono.(2*i-1))+SIN(p(i))*(1.0_dp.mono.(2*i))
          s1%v(2*i)  =COS(p(i))*(1.0_dp.mono.(2*i))-SIN(p(i))*(1.0_dp.mono.(2*i-1))
          s1i%v(2*i-1)=COS(p(i))*(1.0_dp.mono.(2*i-1))-SIN(p(i))*(1.0_dp.mono.(2*i))
          s1i%v(2*i)  =COS(p(i))*(1.0_dp.mono.(2*i))+SIN(p(i))*(1.0_dp.mono.(2*i-1))
       enddo
       if(.not.courant_snyder) then
        s1i=1
        s1=1
       endif
       a_nl=s1i*a_nl*s1
       a_l=a_l*s1
       dr=s1i
       !       a_t_cs= a_f*a_l_cs* a_nl     ! a_nl is not yet " standard "
!!!!!   nonlinear part if s1i !!!!!!
       if(c_%no>1.and.courant_snyder) then
          if(onelie) then
             nd_used=c_%nd2/2
             if(c_%ndpt/=0) nd_used=nd_used-1
             call alloc(vr)
             call alloc(un)

             s1i=a_nl   ! preserve it
             s1=1
             do i=1,c_%no   ! i=1,c_%no
                un%eps=-c_%no
                un=s1i
                vr=un%vector
                decal=1   ! not coasting plane
                do k=1,nd_used*2
                   i_phase=k
                   if(mod(k,2)==0) then
                      i_plane=k/2
                   else
                      i_plane=(k+1)/2
                   endif
                   call cfu(vr%cos%v(k),dphase,vr%cos%v(k))
                   call cfu(vr%sin%v(k),dphase,vr%sin%v(k))
                enddo

                do k=nd_used*2+1,c_%nd2   !  fix coasting later
                   vr%cos%v(k)=0.0_dp
                   vr%sin%v(k)=0.0_dp
                enddo

                un%vector=vr
                s1=texp(un%vector,s1)
                s1i=texp(un%vector,s1i)
             enddo    ! i=1,c_%no

!!!! now fix coasting !!!!!
             if( nd_used == c_%nd) then
                a_nl=s1i
             else  !  coasting
                un%eps=-c_%no
                un=s1
                alpha0=1.0_dp
                if(mod(c_%ndpt,2)==0) alpha0=-1.0_dp
                do i=1,c_%nd2-2
                   un%vector%v(i)=alpha0*(un%vector%v(i).d.c_%ndpt)
                enddo
                call int_partial(un%vector,un%pb,2)

                if(mod(c_%ndpt,2)==0) then
                   s1%v(c_%ndpt-1)=s1%v(c_%ndpt-1)+un%pb%h
                else
                   s1%v(c_%ndpt+1)=s1%v(c_%ndpt+1)+un%pb%h
                endif

                a_nl=a_nl*s1
             endif   ! end of nd_used == c_%nd

!!!! end of now fix coasting !!!!!

             dr=dr*s1**(-1)
             call kill(vr)
             call kill(un)
          else  !  reverse Dragt-Finn
             nd_used=c_%nd2/2
             if(c_%ndpt/=0) nd_used=nd_used-1
             call alloc(vr)
             call alloc(rdf)

             s1i=a_nl   ! preserve it
             s1=1
             do i=1,c_%no   ! i=1,c_%no
                rdf=s1i
                vr=rdf%nonlinear
                decal=1   ! not coasting plane
                do k=1,nd_used*2
                   i_phase=k
                   if(mod(k,2)==0) then
                      i_plane=k/2
                   else
                      i_plane=(k+1)/2
                   endif
                   call cfu(vr%cos%v(k),dphase,vr%cos%v(k))
                   call cfu(vr%sin%v(k),dphase,vr%sin%v(k))
                enddo

                do k=nd_used*2+1,c_%nd2   !  fix coasting later
                   vr%cos%v(k)=0.0_dp
                   vr%sin%v(k)=0.0_dp
                enddo

                rdf%nonlinear=vr
                s1=texp(rdf%nonlinear,s1)
                s1i=texp(rdf%nonlinear,s1i)
             enddo    ! i=1,c_%no

!!!! now fix coasting !!!!!
             if( nd_used == c_%nd) then
                a_nl=s1i
             else  !  coasting
                rdf=s1
                alpha0=1.0_dp
                if(mod(c_%ndpt,2)==0) alpha0=-1.0_dp
                do i=1,c_%nd2-2
                   rdf%nonlinear%v(i)=alpha0*(rdf%nonlinear%v(i).d.c_%ndpt)
                enddo
                call int_partial(rdf%nonlinear,rdf%pb,2)

                if(mod(c_%ndpt,2)==0) then
                   s1%v(c_%ndpt-1)=s1%v(c_%ndpt-1)+rdf%pb%h
                else
                   s1%v(c_%ndpt+1)=s1%v(c_%ndpt+1)+rdf%pb%h
                endif

                a_nl=a_nl*s1
             endif   ! end of nd_used == c_%nd

!!!! end of now fix coasting !!!!!

             dr=dr*s1**(-1)
             call kill(vr)
             call kill(rdf)
          endif
       endif   ! no>1


!!!!
    endif   ! it (te)
    !       endif   ! present(te)


    dote=present(R_TE).and.present(CS_TE)
    if(doing_ac_modulation_in_ptc) then
       dote=dote.and.(c_%nd2==6.or.c_%ndpt>=5)
    else
       dote=dote.and.(c_%nd2==4.or.c_%ndpt>=5)
    endif
    if(dote) then

       call  copy_damap_matrix(a_l,m)
       call  copy_matrix_matrix(m(1:2,1:2),at)
       call  copy_matrix_matrix(m(1:2,3:4),ct)
       call  copy_matrix_matrix(m(3:4,1:2),dt)
       call  copy_matrix_matrix(m(3:4,3:4),bt)
        
       call invert_22(at,ati)
       call invert_22(bt,bti)
       if(.not.c_%STABLE_DA) then
        t_e=my_false
       endif 


       call matmul_nn(dt,ati,ati,sc=-1.0_dp)
       call matmul_nn(ati,ct,ct)
       call matmul_nn(ct,bti,ct)
       if(.not.c_%STABLE_DA) then
        t_e=my_false
        goto 888
       endif 


       alpha=ct(1,1)
       alpha0=alpha

       if(alpha0<=-1.0_dp) then
        t_e=my_false
        goto 888
       endif
        
       det=sqrt(1.0_dp/(1.0_dp+alpha))




       if(alpha0>=0.0_dp) then
          COSLIKE=my_true
       else
          ! det=sqrt(one/(one-alpha))
          COSLIKE=my_false
       endif
 
        CS_TE=0
       do i=1,2
          do j=1,2
             CS_TE%v(i)=at(i,j)*(1.0_dp.mono.j)/det+CS_TE%v(i)
             CS_TE%v(i+2)=bt(i,j)*(1.0_dp.mono.(j+2))/det+CS_TE%v(i+2)
          enddo
       enddo
 

        
       !  The rotation matrix is created but it may not have the correct path length
       !dependence
       if(c_%ndpt/=0.and.t_e) then
          call alloc(g,c_%nv)
          call alloc(un)
          ! write(6,*) " epseone "
          ! read(5,*) un%eps
          un%eps=-c_%no
          do i=1,c_%nv
             g%v(i)=1.0_dp.mono.i
          enddo
          g%v(c_%ndpt)=0.0_dp

          do i=1,c_%nd2
             s1%v(i) = CS_TE%v(i)*g
          enddo
!          s1%v(c_%ndpt)=one.mono.c_%ndpt
!!!!!new
           s1%v(c_%nd2)=1.0_dp.mono.c_%nd2
           s1%v(c_%nd2-1)=1.0_dp.mono.c_%nd2-1


           CS_TE%v(c_%nd2)=1.0_dp.mono.c_%nd2
           CS_TE%v(c_%nd2-1)=1.0_dp.mono.c_%nd2-1
!!!!!
          s1i=s1**(-1)
          s1i=s1i*CS_TE    ! s1i is completely nonlinear.
          un=s1i

          alpha0=1.0_dp
          if(mod(c_%ndpt,2)==0) alpha0=-1.0_dp
          do i=1,c_%nd2-2
             un%vector%v(i)=alpha0*(un%vector%v(i).d.c_%ndpt)
          enddo
          call int_partial(un%vector,un%pb,2)
          !  un%pb=un%vector  ! this is the longitudinal part
          if(mod(c_%ndpt,2)==0) then
             s1i%v(c_%ndpt-1)=s1i%v(c_%ndpt-1)+un%pb%h
          else
             s1i%v(c_%ndpt+1)=s1i%v(c_%ndpt+1)+un%pb%h
          endif

   
          CS_TE=s1*s1i



          call kill(un)
          call kill(g)
          !endif !eps_te
       endif
       
         888 continue
       if(.not.t_e) then       
        c_%STABLE_DA=my_true
        cs_te=0
        R_TE=0
        write(6,*) " Teng-Edwards is crap !"
       else       
        R_TE=a_l*cs_TE**(-1)
       endif  

    endif ! end of T-E  done

    if(lagrange0) then
       a_t=a_f*a_l*a_nl
    endif

    if(doing_ac_modulation_in_ptc.and.present(CS_TE)) then   ! removing useless tiny numbers
       CS_TE%v(c_%nd2-1)=1.0_dp.mono.(c_%nd2-1)
       CS_TE%v(c_%nd2)=1.0_dp.mono.(c_%nd2)
    endif

    if(doflip) then
       call flip_damap(a_t,a_t)
       call flip_damap(a_nl,a_nl)
       call flip_damap(a_l,a_l)
       call flip_damap(a_f,a_f)
       if(present(dr))  call flip_damap(dr,dr)
       if(present(CS_TE))  call flip_damap(CS_TE,CS_TE)
       if(present(R_TE))   call flip_damap(R_TE,R_TE)
       perform_flip=.true.
    endif

    deallocate(jc)
    call kill(atemp,s1,s1i,m1)
    call kill(a12,a11)
    call kill(p,ndim)
    call kill_nn(m)
    call kill_nn(at)
    call kill_nn(bt)
    call kill_nn(ct)
    call kill_nn(dt)
    call kill_nn(ati)
    call kill_nn(bti)
    call kill(alpha,det)


  end subroutine factor_am

  subroutine factor_am_special(a_t,a_f,a_l,a_nl,DR)
    implicit none
    TYPE(damap), INTENT(INout) :: a_t,a_nl,a_l,a_f
    integer i,n,k,nt,j
    integer, allocatable :: jc(:)
    real(dp) value,alpha0
    logical doit
    TYPE(damap) atemp,s1,s1i,m1
    TYPE(taylor)a12,a11,p(ndim)
    logical(lp) doflip,dote,lagrange0
    TYPE(damap), optional, intent(inout) ::DR
 !   type(taylor) m(ndim2,ndim2)
 !   type(taylor) at(2,2),bt(2,2),ct(2,2),dt(2,2),ati(2,2),bti(2,2),alpha,det
 !   type(gmap) g
 !   type(onelieexponent) un
 !   type(reversedragtfinn) rdf
 !   type(vecresonance) vr


    lagrange0=my_false
    if(present(dr)) lagrange0=my_true

    call alloc(atemp,s1,s1i,m1)
    call alloc(a12,a11)
    call alloc(p,ndim)
!    call alloc_nn(m)
!    call alloc_nn(at)
!    call alloc_nn(bt)
!    call alloc_nn(ct)
!    call alloc_nn(dt)
!    call alloc_nn(ati)
!    call alloc_nn(bti)
!    call alloc(alpha,det)
    allocate(jc(c_%nv))



    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       !    write(6,*) " flipping ",c_%ndpt,c_%nd2-1
       !    pause 7123

       perform_flip=.false.
       call flip_damap(a_t,a_t)
       doflip=.true.
    else
       doflip=.false.
    endif

    atemp=1
    do k=1,c_%nd2
       call taylor_cycle(a_t%v(k),N)

       do i=1,n
          call taylor_cycle(a_t%v(k),ii=i,value=value,j=jc)
          call check_fix(jc,0,doit)
          if(doit) atemp%v(k)= atemp%v(k) + (value.mono.jc)
       enddo
    enddo

    !   in case we have ndpt/=0

    if(c_%ndpt==c_%nd2-1) then   ! ptc convention
       atemp%v(c_%nd2-1)=atemp%v(c_%nd2-1)-(1.0_dp.mono.(c_%nd2-1))
       atemp%v(c_%nd2)=atemp%v(c_%nd2)-(1.0_dp.mono.(c_%nd2))
       do i=1,c_%nd-1
          atemp%v(c_%ndpt+1)=atemp%v(c_%ndpt+1) &
               +(atemp%v(2*i).d.c_%ndpt)*(1.0_dp.mono.2*i-1)-(atemp%v(2*i-1).d.c_%ndpt)*(1.0_dp.mono.2*i)
       enddo
    elseif(c_%ndpt==c_%nd2) then      ! Marylie convention
       atemp%v(c_%nd2-1)=atemp%v(c_%nd2-1)-(1.0_dp.mono.(c_%nd2-1))
       atemp%v(c_%nd2)=atemp%v(c_%nd2)-(1.0_dp.mono.(c_%nd2))
       do i=1,c_%nd-1
          atemp%v(c_%ndpt-1)=atemp%v(c_%ndpt-1) &
               -(atemp%v(2*i).d.c_%ndpt)*(1.0_dp.mono.2*i-1)+(atemp%v(2*i-1).d.c_%ndpt)*(1.0_dp.mono.2*i)
       enddo
    endif


    a_f=atemp   !!!  not exact ! temporal part could be wrong!!!

    a_l=a_f**(-1)*a_t

    atemp=0
    do k=1,c_%nd2
       call taylor_cycle(a_l%v(k),N)

       do i=1,n
          call taylor_cycle(a_l%v(k),ii=i,value=value,j=jc)
          call check_fix(jc,1,doit)
          if(doit) atemp%v(k)= atemp%v(k) + (value.mono.jc)
       enddo
    enddo


    !   in case we have ndpt/=0

    if(c_%ndpt/=0) then   ! ptc convention
       atemp%v(c_%nd2-1)=(1.0_dp.mono.(c_%nd2-1))
       atemp%v(c_%nd2)= (1.0_dp.mono.(c_%nd2))

       if(c_%ndpt==c_%nd2) then
          k=c_%nd2-1
       else
          k=c_%nd2
       endif

       call taylor_cycle(a_l%v(k),N)

       do i=1,n
          call taylor_cycle(a_l%v(k),ii=i,value=value,j=jc)
          call check_fix(jc,2,doit)
          if(doit) atemp%v(k)= atemp%v(k) + (value.mono.jc)
       enddo

    endif

    a_nl=atemp**(-1)*a_l
    a_l=atemp

!!! a_t=a_f*a_l*a_nl at this stage

    !   if(present(lagrange)) then
    if(lagrange0) then
       !   enforce Teng-Edwards to all orders in  parameters i.e. A_12=0 and A_3_4=0 etc....
       nt=c_%nd
       if(C_%ndpt/=0) nt=nt-1

       s1=1
       s1i=1

       if(C_%ndpt/=0) then
          s1%v(C_%ndpt)=1.0_dp.mono.C_%ndpt
          s1i%v(C_%ndpt)=1.0_dp.mono.C_%ndpt
       endif

       do i=1,nt
          a11=a_l%v(2*i-1).d.(2*i-1)
          a12=a_l%v(2*i-1).d.(2*i)
          p(i)=-a12/a11
          p(i)=atan(p(i))
          if(C_%ndpt/=0) then
             if(mod(C_%ndpt,2)==1) then
                s1%v(C_%ndpt+1)=s1%v(C_%ndpt+1)-(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
                s1i%v(C_%ndpt+1)=s1i%v(C_%ndpt+1)+(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
             else
                s1%v(C_%ndpt-1)=s1%v(C_%ndpt-1)+(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
                s1i%v(C_%ndpt-1)=s1i%v(C_%ndpt-1)-(p(i).d.C_%ndpt)*((1.0_dp.mono.(2*i))**2+(1.0_dp.mono.(2*i-1))**2)*0.5_dp
             endif
          endif
          s1%v(2*i-1)=COS(p(i))*(1.0_dp.mono.(2*i-1))+SIN(p(i))*(1.0_dp.mono.(2*i))
          s1%v(2*i)  =COS(p(i))*(1.0_dp.mono.(2*i))-SIN(p(i))*(1.0_dp.mono.(2*i-1))
          s1i%v(2*i-1)=COS(p(i))*(1.0_dp.mono.(2*i-1))-SIN(p(i))*(1.0_dp.mono.(2*i))
          s1i%v(2*i)  =COS(p(i))*(1.0_dp.mono.(2*i))+SIN(p(i))*(1.0_dp.mono.(2*i-1))
       enddo
       if(.not.courant_snyder) then
        s1i=1
        s1=1
       endif
       a_nl=s1i*a_nl*s1
       a_l=a_l*s1
       dr=s1i
       !       a_t_cs= a_f*a_l_cs* a_nl     ! a_nl is not yet " standard "
!!!!!   nonlinear part if s1i !!!!!!



!!!!
    endif   ! it (te)
 



    if(lagrange0) then
       a_t=a_f*a_l*a_nl
    endif


    if(doflip) then
       call flip_damap(a_t,a_t)
       call flip_damap(a_nl,a_nl)
       call flip_damap(a_l,a_l)
       call flip_damap(a_f,a_f)
       if(present(dr))  call flip_damap(dr,dr)

       perform_flip=.true.
    endif

    deallocate(jc)
    call kill(atemp,s1,s1i,m1)
    call kill(a12,a11)
    call kill(p,ndim)
!    call kill_nn(m)
!    call kill_nn(at)
!    call kill_nn(bt)
!    call kill_nn(ct)
!    call kill_nn(dt)
!    call kill_nn(ati)
!    call kill_nn(bti)
!    call kill(alpha,det)


  end subroutine factor_am_special

  subroutine int_partial(v,h,nd0)
    implicit none
    ! IF SCA=-one
    !     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:
    !
    ! IF SCA=one
    !     \VEC{V}.GRAD  = GRAD H . GRAD
    integer i,nd0
    type(vecfield) v
    type(pbfield) h
    type(taylor) b4,b3,b2,b1
    type(damap) x
    logical doflip
    if(.not.c_%stable_da) return

    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       !    write(6,*) " flipping ",c_%ndpt,c_%nd2-1
       !    pause 7123

       perform_flip=.false.
       call flip_vecfield(v,v,1)
       call flip_taylor(h%h,h%h,1)
       doflip=.true.
    else
       doflip=.false.
    endif

    nd_used=nd0
    call alloc(x)
    call alloc(b4,b3,b2,b1)

    x=1

    do i=1,nd_used
       call cfu(v%v(2*i-1),dlie,b3)
       call cfu(v%v(2*i),dlie,b1)
       b2=b1*x%v(2*i-1)
       b1=b3*x%v(2*i)
       b3=b2-b1
       b2=b3+b4
       b4=b2
    enddo
    h%h=b4


    call kill(b4,b3,b2,b1)
    call kill(x)

    if(doflip) then
       call flip_vecfield(v,v,-1)
       call flip_taylor(h%h,h%h,-1)
       perform_flip=.true.
    endif
  end subroutine int_partial

  real(dp) function dlie(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    if(.not.c_%stable_da) return

    dlie=0.0_dp
    do i=1,nd_used
       dlie=REAL(j(2*i-1)+j(2*i),kind=DP)+dlie
    enddo
    dlie=dlie+1.0_dp
    dlie=1.0_dp/dlie
    return
  end function dlie


  real(dp) function dphase(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    integer t,tu
    if(.not.c_%stable_da) return

    t=-decal
    tu=0
    dphase=0.0_dp
    do i=1,nd_used
       t=abs(j(2*i-1)-j(2*i))+t
    enddo
    if(t==0.and.decal/=0) then
       tu=(j(2*i_plane-1)-j(2*i_plane))
       if(mod(i_phase,2)==1) then
          tu=tu-1
       else
          tu=tu+1
       endif
       if(tu==0) dphase=-1.0_dp
    else
       if(t==0) dphase=-1.0_dp
    endif

    return
  end function dphase

  real(dp) function phase_shift(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    integer nd,t
    if(.not.c_%stable_da) return
    nd=c_%nd2/2
    if(c_%ndpt/=0) nd=nd-1


    phase_shift=0.0_dp
    t=0
    do i=1,nd
       t=abs(j(2*i-1)-j(2*i))+t
    enddo
    if(t==0) phase_shift=1.0_dp

    return
  end function phase_shift


  subroutine invert_22(a,ai)
    implicit none
    type(taylor) a(2,2),ai(2,2),t(2,2)
    type(taylor) det
    call alloc_nn(t)
    call alloc(det)

    det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    t(1,1)=a(2,2)/det
    t(2,2)=a(1,1)/det
    t(1,2)=-a(1,2)/det
    t(2,1)=-a(2,1)/det

    ai(1,1)=t(1,1)
    ai(1,1)=t(1,1)
    ai(1,1)=t(1,1)
    ai(1,1)=t(1,1)

    call copy_matrix_matrix(t,ai)

    call kill_nn(t)
    call kill(det)

  end subroutine invert_22


  subroutine check_fix(jc,ord,doit)
    implicit none
    logical doit
    integer i,nd2,l,ord
    integer  jc(:)

    nd2=c_%nd2
    if(c_%ndpt/=0) nd2=c_%nd2-2

    doit=.false.


    l=0
    do i=1,nd2
       l=jc(i)+l
    enddo
    if(l==ord) doit=.true.

  end subroutine check_fix


  subroutine INTO_RES_SPIN8_eq(H_RES,H)
    implicit none
    TYPE(spinor_8), INTENT(in) :: h
    TYPE(res_spinor_8), INTENT(inout) :: H_RES
    call   INTO_RES_SPIN8(H,H_RES)
  end subroutine INTO_RES_SPIN8_eq

  subroutine INTO_RES_SPIN8(H,H_RES)
    implicit none
    TYPE(spinor_8), INTENT(in) :: h
    TYPE(res_spinor_8), INTENT(inout) :: H_RES
    type(taylorresonance) tr,ti

    call alloc( tr)
    call alloc( ti)


    H_RES%X(1)= H%X(1)-i_*H%X(3)
    H_RES%X(2)= H%X(1)+i_*H%X(3)
    H_RES%X(3)= H%X(2)

    tr=H_RES%X(3)%t%r

    H_RES%X(3)= tr%cos + i_*tr%sin

    tr=H_RES%X(1)%t%r
    ti=H_RES%X(1)%t%i

    H_RES%X(1)=tr%cos + i_*tr%sin + i_* (ti%cos + i_*ti%sin)

    tr=H_RES%X(2)%t%r
    ti=H_RES%X(2)%t%i

    H_RES%X(2)=tr%cos + i_*tr%sin + i_* (ti%cos + i_*ti%sin)


    call kill( tr)
    call kill( ti)


  END SUBROUTINE INTO_RES_SPIN8

  subroutine INTO_SPIN8_from_RES_eq(H,H_RES)
    implicit none
    TYPE(spinor_8), INTENT(inout) :: h
    TYPE(res_spinor_8), INTENT(in) :: H_RES
    type(taylorresonance) tr
    type(taylor) t
    type(complextaylor) c

    call alloc( tr)
    call alloc( t)
    call alloc( c)

    tr%cos=H_RES%X(3)%t%r
    tr%sin=H_RES%X(3)%t%i
    t=tr
    H%X(2)=morph(t)
    c=(H_RES%X(1)+H_RES%X(2))/2.0_dp

    tr%cos=c%r
    tr%sin=c%i
    t=tr
    H%X(1)=morph(t)

    c=i_*(H_RES%X(1)-H_RES%X(2))/2.0_dp

    tr%cos=c%r
    tr%sin=c%i
    t=tr
    H%X(3)=morph(t)

    call kill( c)
    call kill( t)
    call kill( tr)


  END SUBROUTINE INTO_SPIN8_from_RES_eq

!!!!!!!!!!    normal form on theta into theta(H) !!!!!!!!!!!!!!!!

  subroutine normal_thetaH(ds,ns)
    implicit none
    TYPE(damapspin), INTENT(INout) :: ds
    TYPE(normal_spin), INTENT(INout) :: ns
    TYPE(onelieexponent) uno
    TYPE(pbfield) h
    TYPE(damapspin) nc
    TYPE(normalform) nf
    TYPE(taylorresonance) tr,theta0r
    TYPE(taylor) theta0
    TYPE(complextaylor) gam
    real(dp) muc,epsy,ath,vvb
    complex(dp) dely,bth,beta,alpha


    muc=0.25e0_dp*twopi
    call alloc(h)
    call alloc(uno)
    call alloc(nc)
    call alloc(tr)
    call alloc(theta0r)
    call alloc(theta0)
    call alloc(gam)
    call alloc(nf)

    nc=ns%a_t**(-1)*ds*ns%a_t

    h%h=muc*((1.0_dp.mono.'002')+(1.0_dp.mono.'0002'))/2.0_dp
    h%h=h%h+muc*((1.0_dp.mono.'2')+(1.0_dp.mono.'02'))/2.0_dp

    nc%m=texp(h,nc%m)


    !     call print(nc%m,6)
    !     nf=nc%m

    uno=nc%m
    uno%pb%h=uno%pb%h-muc*((1.0_dp.mono.'2')+(1.0_dp.mono.'02'))/2.0_dp
    tr=uno%pb%h

    call print(ns%theta0,6)

    !    a11=(s(1,1))
    !    a13=(s(1,3))

    theta0=atan2(nc%s%s(1,3),nc%s%s(1,1))
    if((theta0.sub.'0')<0.0_dp) theta0 = theta0 + twopi

    theta0r=theta0

!!    call taylor_clean(theta0r%cos,1.d-1)
!    call taylor_clean(theta0r%sin,1.d-1)
!    call taylor_clean(tr%cos,1.d-5)
!    call taylor_clean(tr%sin,1.d-5)

    call print(theta0r%cos,6)
    call print(tr%cos,6)

    call print(theta0r%sin,6)
    call print(tr%sin,6)

    epsy=tr%cos.sub.'0011'
    dely=(tr%cos.sub.'0040')+i_*(tr%sin.sub.'0040')
    ath=theta0r%cos.sub.'0011'
    bth=(theta0r%cos.sub.'0040')+i_*(theta0r%sin.sub.'0040')
    write(6,*) epsy
    write(6,*) dely
    write(6,*) ath
    write(6,*) bth
    vvb=(0.2d-3)**2

    alpha=(ath*dely/epsy-bth)*i_/8.d0/dely
    write(6,*) alpha
    alpha=abs(alpha)
    write(6,*) alpha
    alpha=0
    beta=((ath*dely/epsy-bth)*i_/8.d0-dely*alpha)/epsy
    beta=abs(beta)
    write(6,*) beta
    alpha=0.1d0
    beta=((ath*dely/epsy-bth)*i_/8.d0-dely*alpha)/epsy
    beta=abs(beta)
    write(6,*) beta
    ! alpha=beta*(dely/epsy)*vvb


    call kill(uno)
    call kill(h)
    call kill(nc)
    call kill(tr)
    call kill(theta0r)
    call kill(theta0)
    call kill(nf)
    call kill(gam)

  end subroutine normal_thetaH


!!! Some useful routines

  subroutine AVERAGE(F,A,F_floquet,F_xp,use_J)
    implicit none
    type(damap) A
    TYPE(TAYLOR), intent(inout):: F
    TYPE(TAYLOR), intent(inout):: F_floquet
    TYPE(TAYLOR), optional, intent(inout):: F_xp
    type(taylor) tt,tt_xp
    TYPE(taylorresonance) fq
    real(dp) value,valuexp
    integer, allocatable :: jc(:)
    logical, optional :: use_J
    logical usej
    integer i,n,j,it,nd,iu
    logical doflip,uj

    if(perform_flip.and.new_ndpt.and.c_%ndpt/=0) then
       perform_flip=.false.
       call fliptaylor(F,F,1)
       call flip_damap(A,A)
       doflip=.true.
    else
       doflip=.false.
    endif

    uj=.false.
    if(present(use_j)) uj=use_j

    nd=c_%nd2/2
    if(c_%ndpt/=0) nd=nd-1

    allocate(jc(c_%nv))
    call alloc(tt,tt_xp)
    call alloc(fq)
    !! USE_J is false by default


    fq=F*A
    fq%sin=0.0_dp

    call taylor_cycle(fq%cos,size=n)

    do i=1,n
       call taylor_cycle(fq%cos,ii=i,value=value,j=jc)

       it=0
       iu=0
       do j=1,nd
          it=iabs(jc(j*2-1)-jc(j*2))+it
          iu=iabs(jc(j*2-1)+jc(j*2))+iu
       enddo
       if(it==0) then
          iu=iu/2
          valuexp=value
          if(uj) then
             value=valuexp*2.0_dp**iu
          endif
          tt=(value.mono.jc)+tt
          tt_xp=(valuexp.mono.jc)+tt_xp
       endif

    enddo

    fq%cos=tt
    F_floquet=tt

    if(present(F_xp)) then
       fq%cos=tt_xp
       F_xp=fq
       F_xp=F_xp*A**(-1)
    endif

    deallocate(jc)
    call kill(tt,tt_xp)
    call kill(fq)

    if(doflip) then
       call fliptaylor(F,F,-1)
       call flip_damap(A,A)
       call fliptaylor(F_floquet,F_floquet,-1)
       if(present(F_xp)) call fliptaylor(F_xp,F_xp,-1)
       perform_flip=.true.
    endif

  end subroutine AVERAGE

  ! remove small numbers

  SUBROUTINE  clean_damapspin(S1,S2,prec)
    implicit none
    type (damapspin),INTENT(INOUT)::S2
    type (damapspin), intent(INOUT):: s1
    real(dp) prec
    integer i,j

    call clean_damap(s1%m,s2%m,prec)
    do i=1,3
       do j=1,3
          call clean_real_8(s1%s%s(i,j),s2%s%s(i,j),prec)
       enddo
    enddo


  END SUBROUTINE clean_damapspin

  SUBROUTINE  clean_spinor_8(S1,S2,prec)
    implicit none
    type (spinor_8),INTENT(INOUT)::S2
    type (spinor_8), intent(INOUT):: s1
    real(dp) prec
    integer i

    do i=1,3
       call clean_real_8(s1%x(i),s2%x(i),prec)
    enddo

  END SUBROUTINE clean_spinor_8

  SUBROUTINE  clean_res_spinor_8(S1,S2,prec)
    implicit none
    type (res_spinor_8),INTENT(INOUT)::S2
    type (res_spinor_8), intent(INOUT):: s1
    real(dp) prec
    integer i

    do i=1,3
       call clean_double_complex(s1%x(i),s2%x(i),prec)
    enddo

  END SUBROUTINE clean_res_spinor_8


  ! TABLE STUFF FOR FUTURE DA

  integer function number_mon(n,m)
    implicit none
    integer i,n,m

    number_mon=1


    do i=n+m,max(n,m)+1,-1

       number_mon=number_mon*i
    enddo

    do i=2,min(n,m)
       number_mon=number_mon/i
    enddo

  end  function number_mon

 integer function mul_fac(ju)
    implicit none
    integer ju(:),nv
    integer i,k,no
    
    mul_fac=1.0_dp
    if(firstfac) then
     call make_fac
     firstfac=.false.
    endif
    nv=size(ju)
    
    k=0
    do i=1,nv
     k=k+ju(i)
    enddo


    do i=1,nv
     if(ju(i)==0) cycle
     
     mul_fac=(fac(k)/fac(k-ju(i))/fac(ju(i)))*mul_fac
     k=k-ju(i)
     
    enddo
    

end function mul_fac

subroutine make_fac()
    implicit none
    integer i

    fac(0)=1.0_dp
    do i=1,nfac
    fac(i)=i*fac(i-1)
    enddo

end subroutine make_fac


  integer function pos_mon(ju,nomax,nv)
    implicit none
    integer ju(:),no,nv,nomax
    integer i,k,nk

    pos_mon=0
    no=0
    do i=1,nv
       no=no+ju(i)
    enddo

    nk=no

    if(nk>nomax) then
       pos_mon=0
       return
    endif

    do k=1,nv-1
       if(ju(k)/=0) then
          pos_mon=pos_mon+number_mon(nk,nv-k)-number_mon(nk-ju(k),nv-k)
          nk=nk-ju(k)
       endif
    enddo
    pos_mon=pos_mon+1

    if(no>0) pos_mon=pos_mon+number_mon(no-1,nv)

  end  function pos_mon

  integer function pos_no(no,nomax,nv)
    implicit none
    integer ju(lnv),no,nv,nomax
 
    if(no==0) then
     pos_no=0
    endif
   if(no>nomax) then
     pos_no=-1
    endif
    ju=0
    ju(nv)=no
    pos_no=pos_mon(ju,nomax,nv)
  end  function pos_no

  subroutine find_exp(p,ju,no,nv)
    implicit none
    integer ju(:),no,nv
    integer i,nk,nvk,p,p0,p1,pg

    ju=0

    if(p==1) then
       return
    endif

    if(p<=nv+1) then
       ju(nv-p+2)=1
       return
    endif


    do i=1,no
       p1=number_mon(i,nv)
       if(p1>=p) then
          nk=i
          exit
       endif
       p0=p1
    enddo



    nvk=nv-1
    pg=p0
    do while(nvk>0)

       p1=pg
       do i=0,nk
          p1=number_mon(nk-i,nvk-1)+p1
          if(p1>=p) then
             nk=nk-i
             ju(nv-nvk)=i
             nvk=nvk-1
             exit
          endif
          pg=p1
       enddo
    enddo
    ju(nv)=nk
  end subroutine find_exp



end module tree_element_MODULE
