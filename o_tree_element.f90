!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tree_element_MODULE
  USE polymorphic_complextaylor, mon1=>mon
  IMPLICIT NONE
  public
  integer,private,parameter::ndd=6

  PRIVATE track_TREE,track_TREEP,KILL_TREE,KILL_TREE_N
  PRIVATE track_TREE_G,track_TREEP_g
  PRIVATE ALLOC_SPINOR_8,ALLOC_probe_8
  PRIVATE KILL_SPINOR_8,KILL_probe_8,KILL_DASPIN
  PRIVATE EQUAL_SPINOR8_SPINOR8,EQUAL_IDENTITY_SPINOR_8 !,EQUAL_SPINOR8_RAY8,EQUAL_RAY8_SPINOR8,
  PRIVATE  EQUAL_IDENTITY_probe_8,ALLOC_daspin,EQUAL_DASPIN_RAY8 ,EQUAL_RAY8_DASPIN
  private alloc_normal_spin,kill_normal_spin,EQUAL_IDENTITY_probe
  private READ_DASPIN,PRINT_DASPIN,EQUAL_IDENTITY_SPINOR,EQUAL_PROBE_REAL6
  PRIVATE EQUAL_SPINOR8_SPINOR,EQUAL_PROBE8_PROBE,EQUAL_PROBE8_REAL6
  private EQUAL_IDENTITY_SPINOR_8_r3 ,EQUAL_SPINOR_SPINOR8

  private find_ar
  private find_ap,PRINT_spinor_8,dot_spinor_8,dot_spinor
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

  integer, target :: spin_extra_tpsa = 0 ,n0_normal= 2
  integer :: clockwise=1
  logical(lp) :: force_positive=.false.

  INTERFACE assignment (=)
     !
     MODULE PROCEDURE REAL_8REAL6
     MODULE PROCEDURE REAL6REAL_8
     MODULE PROCEDURE real_8REAL_8
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
     MODULE PROCEDURE EQUAL_DASPIN_RAY8
     MODULE PROCEDURE EQUAL_RAY8_DASPIN

     MODULE PROCEDURE EQUAL_DAmapSPIN_int  ! damapspin = Identity (or zero)
     MODULE PROCEDURE EQUAL_PROBE8_PROBE8
     MODULE PROCEDURE      EQUAL_damapspin
     MODULE PROCEDURE     EQUAL_SPINOR_SPINOR8
     MODULE PROCEDURE     EQUAL_PROBE_PROBE8
     MODULE PROCEDURE     normalise_spin

  end  INTERFACE
  INTERFACE OPERATOR (*)
     MODULE PROCEDURE concat  !! damapspin= damapspin o damapspin
     MODULE PROCEDURE cmul
     MODULE PROCEDURE mmul  ! damapspin%s=damapspin%s o damap
     MODULE PROCEDURE damapspin_spinor_mul  ! spinor_8%s(i)= damapspin%s(i,j) * spinor%s(j)
     MODULE PROCEDURE damapspin_spinor8_mul ! spinor_8%s(i)= damapspin%s(i,j) * spinor_8%s(j)
     MODULE PROCEDURE spin8_mul_map   ! spinor_8=spinor_8 o damap
     MODULE PROCEDURE eval_spinor_8   ! spinor=spinor_8 * xp(lnv)
     MODULE PROCEDURE concatxp   !  damapspin=damapspin*xp(lnv)
     MODULE PROCEDURE spin8_scal8_map  ! real_8*spinor_8
     MODULE PROCEDURE mul_spin8_spin8
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
     MODULE PROCEDURE addm
     MODULE PROCEDURE add_spin8_spin8
  END  INTERFACE

  INTERFACE operator (-)
     MODULE PROCEDURE sub_spin8_spin8
  END  INTERFACE

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


  INTERFACE PRINT
     MODULE PROCEDURE PRINT6
!!!!
     MODULE PROCEDURE PRINT_DASPIN
     MODULE PROCEDURE PRINT_probe8
     MODULE PROCEDURE PRINT_spinor_8
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
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILL_SPINOR_8
     MODULE PROCEDURE KILL_probe_8
     MODULE PROCEDURE KILL_DASPIN
     MODULE PROCEDURE K_OPT_damap
     MODULE PROCEDURE KILL_normal_spin
  END INTERFACE

  INTERFACE ALLOC_33
     MODULE PROCEDURE ALLOC_33t
     MODULE PROCEDURE ALLOC_33p
  END INTERFACE

  INTERFACE zero_33
     MODULE PROCEDURE zero_33t
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
          daddsco(i)=s2(i)+(one.mono.'00001')
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
          scdaddo(i)=s2(i)+(one.mono.'00001')
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

    ALLOCATE(T%CC(N),T%JL(N),T%JV(N),T%N,T%ND2,T%no)
    T%N=N
    T%ND2=ND2
    T%no=0

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

    XT=zero
    XF=zero
    XM=zero

    do i=1,T%ND2
       xt(i)=xi(i)
    enddo
    do i=1,T%ND2
       xf(i) = T%cc(i)
    enddo

    XM(1) = one
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

    XM(1) = one
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


    IF(ASSOCIATED(T%CC))   DEALLOCATE(T%CC,T%JL,T%JV,T%N,T%ND2,T%No)


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

    n1=1
    if(present(n)) n1=n
    do k=1,n1
       if(.not.c_%CHECK_STABLE) return
       XT=zero
       XF=zero
       XM=zero

       do i=1,T%ND2
          xt(i)=xi(i)
       enddo
       do i=1,T%ND2
          xf(i) = T%cc(i)
       enddo

       XM(1) = one
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

       XM(1) = one
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
       if(eps>zero) uno%eps=eps
    endif
    uno=m
    id=1
    uno%pb%h=uno%pb%h/nst
    m=texp(uno%pb,id)

    call kill(uno); call kill(id);


  end SUBROUTINE symplectic

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
       CALL WRITE_E(100)
    end select

    do i=1,3
       do j=1,3
          call assp_no_master(s1%s(i,j))
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
       CALL WRITE_E(100)
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


  end subroutine assprobe_8




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
          CUTORDER%s(i,j)=s1%s(i,j).cut.s2
       enddo
    enddo

    CUTORDER%m=s1%m.cut.(s2+1)

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
          s%s(i,j)=t%s(j,i)
          !          if(s%s(i,j)%kind==2) then   ! taylor
          s%s(i,j)=s%s(i,j)*s%m
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
          s%s(i,j)=r%s(i,j)
       enddo
    enddo

    s%m=r%m
    !    s%s0=r%s0

    !     do i=1,6
    !     s%x(i)=r%x(i)
    !     enddo

  END    subroutine EQUAL_damapspin




  subroutine EQUAL_IDENTITY_SPINOR_8(S,R)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    INTEGER, INTENT(IN) :: R
    INTEGER I

    !     S%G=A_PARTICLE
    IF(R==1.or.r==2.or.r==3) THEN
       DO I=1,3    !C_%NSPIN
          S%X(I)=ZERO
       enddo
       S%X(r)=one
       !        write(6,*) " in EQUAL_IDENTITY_SPINOR_8"
       !        stop 888
       !       DO I=1,C_%NSPIN
       !          S%X(I)=ONE.MONO.(c_%npara_fpp-C_%NSPIN+I)
       !       ENDDO
    ELSEIF(R==0) THEN
       DO I=1,3    !C_%NSPIN
          S%X(I)=ZERO
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
          S%X(I)=ZERO
       enddo
       S%X(r)=one
    ELSEIF(R==0) THEN
       DO I=1,3   !C_%NSPIN
          S%X(I)=ZERO
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
       P%s(i)%x=zero
       P%s(i)%x(i)=one
    enddo
    P%X=X

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
       P%s(i)%x(i)=one
    enddo

    DO I=1,6
       P%X(i)=X(i)
    enddo
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

    P8%u=P%u

  END subroutine EQUAL_PROBE8_PROBE8

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
    P8%e_ij=zero
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
       R%X(I)=zero
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
       R%X(I)=ZERO
    enddo
    IF(S==1) THEN
       DO I=1,c_%npara_fpp  !-C_%NSPIN
          R%X(I)=ONE.MONO.I
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

  subroutine EQUAL_DASPIN_RAY8(DS,R)
    implicit none
    TYPE(probe_8), INTENT(IN) :: R
    TYPE(damapspin), INTENT(INOUT) :: DS
    real(dp) s(3,3)
    INTEGER I,J
    ! if(C_%ND2/=6) then
    !  write(6,*) " not implemented in EQUAL_DASPIN_RAY8"
    !  stop
    ! endif

    DO I=1,C_%ND2
       DS%M%V(I)=R%X(I)
    ENDDO

    DO I=1,3
       DO J=1,3
          DS%S(I,J)=R%S(J)%X(I)
       ENDDO
    ENDDO
    !          DS%S(I,1)=R%SX%X(I)
    !          DS%S(I,2)=R%SY%X(I)
    !          DS%S(I,3)=R%SZ%X(I)


    DO I=1,3
       !     n(i)=R%S(0)%X(I)
       DO J=1,3
          S(I,J)=DS%S(I,J)
       ENDDO
    ENDDO

    !    call inv_as(s,s)

    !    n= matmul(s,n)

    !    ds%s0=n



    !     DS%G=R%S%G
  END subroutine EQUAL_DASPIN_RAY8

  subroutine EQUAL_RAY8_DASPIN(R,DS)
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
          R%S(J)%X(I)=DS%S(I,J)                  ! taylor MORPH(DS%S(I,J))
       ENDDO
    ENDDO





    !    R%S%G=DS%G

  END subroutine EQUAL_RAY8_DASPIN

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
    CALL zero_33(DS%S)

    DO I=1,3
       if(i1==1) then
          DS%S(I,I)=ONE
       endif
    ENDDO


  END SUBROUTINE EQUAL_DAmapSPIN_int

  FUNCTION scdadd( S2,S1  )
    implicit none
    TYPE (probe_8) scdadd
    TYPE (damapspin), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j
    type(taylor) d




    !   call ass(scdadd)
    scdadd%u=my_false
    scdadd%E_ij=zero

    call alloc(d)
    do i=1,c_%nd2
       localmaster=master
       call ass(scdadd%x(i))
       !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
       scdadd%x(i)=s1%M%V(i)+s2%x(i)-(s1%M%V(i).sub.'0')
       master=localmaster
    enddo

    do i=c_%nd2+1,6
       localmaster=master
       call ass(scdadd%x(i))
       if(c_%nd2==4.and.(c_%npara==5.or.c_%npara==8).AND.I==5) then   ! npr
          scdadd%x(i)=s2%x(i)+(one.mono.'00001')
       else
          scdadd%x(i)=s2%x(i)
       endif
       master=localmaster
    enddo


    !    if(C_%SPIN_POS/=0) then
    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(scdadd%s(J)%x(i))
          scdadd%s(J)%x(i)=S1%S(I,J)
          master=localmaster
       ENDDO
       !          localmaster=master
       !          call ass(scdadd%s(0)%x(i))
       !          scdadd%s(0)%x(i)=S2%S(0)%X(I)
       !          master=localmaster

    ENDDO
    !    endif



  END FUNCTION scdadd


  FUNCTION daddsc( S1,S2  )
    implicit none
    TYPE (probe_8) daddsc
    TYPE (damapspin), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j
    type(taylor) d




    !   call ass(daddsc)
    daddsc%u=my_false
    daddsc%E_ij=zero

    call alloc(d)
    do i=1,c_%nd2
       localmaster=master
       call ass(daddsc%x(i))
       !       daddsc%x(i)=s1%m%v(i)+s2%x(i)
       daddsc%x(i)=s1%M%V(i)+s2%x(i)-(s1%M%V(i).sub.'0')
       master=localmaster
    enddo

    do i=c_%nd2+1,6
       localmaster=master
       call ass(daddsc%x(i))
       if(c_%nd2==4.and.(c_%npara==5.or.c_%npara==8).AND.I==5) then   ! npr
          daddsc%x(i)=s2%x(i)+(one.mono.'00001')
       else
          daddsc%x(i)=s2%x(i)
       endif
       master=localmaster
    enddo


    !    if(C_%SPIN_POS/=0) then
    DO I=1,3
       !          call ass(daddsc%s%x(i))

       DO J=1,3
          localmaster=master
          call ass(daddsc%s(J)%x(i))
          daddsc%s(J)%x(i)=S1%S(I,J)
          master=localmaster
       ENDDO
       !          localmaster=master
       !          call ass(daddsc%s(0)%x(i))
       !          daddsc%s(0)%x(i)=S2%S(0)%X(I)
       !          master=localmaster


    ENDDO
    !    endif



  END FUNCTION daddsc

  subroutine print_DASPIN(DS,MF,prec)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: DS
    real(dp), optional :: prec
    INTEGER MF,I,J

    !    WRITE(MF,*) " ORBIT "
    !     WRITE(MF,*) DS%X(1:3)
    !     WRITE(MF,*) DS%X(4:6)
    WRITE(MF,*) c_%nd2 ," DIMENSIONAL ORBITAL MAP "
    CALL PRINT(DS%M,MF,prec)
    WRITE(MF,*) " SPIN MATRIX "
    DO I=1,3
       DO J=1,3
          WRITE(MF,*) " SPIN MATRIX COMPONENT ",I,J
          CALL PRINT(DS%S(I,J),MF,prec)
       ENDDO
    ENDDO

  END subroutine print_DASPIN

  subroutine print_probe8(DS,MF)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: DS
    INTEGER MF,I

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


  END subroutine print_probe8

  subroutine print_spinor_8(S,MF)
    implicit none
    TYPE(spinor_8), INTENT(INOUT) :: s
    INTEGER MF,I

    do i=1,3
       write(mf,*) ' Spin Variable ',i
       call print(s%x(i),mf)
    enddo

  END subroutine print_spinor_8

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
          ds%m%v(i)=one.mono.i
       enddo
    endif
    READ(MF1,*) LINE
    DO I=1,3
       DO J=1,3
          READ(MF1,*) LINE
          CALL READ(T,MF1)
          DS%S(I,J)=T     !MORPH(T)
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
       ss(i,i)=ss(i,i)-one
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

    n0(is)=one
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
       ss(i,i)=ss(i,i)-one
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



    n0(is)=one
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
          ss(i,j)=ds%s(i,j)
       enddo
    enddo



    do i=1,3
       ss(i,i)=ss(i,i)-one
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



    n0(is)=one
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

    norm=sqrt(z_axis.dot.z_axis)
    z_axis=x_axis*y_axis
    z_axis=(1.d0/norm)*z_axis



    call kill(norm)

  end subroutine find_perp_basisp

  subroutine find_exponentp(ds,y_axis,x_axis,z_axis,h_axis)  !
    implicit none
    type(damapspin),intent(inout) :: ds
    type(spinor_8),intent(inout) :: y_axis,x_axis,z_axis,h_axis
    type(spinor_8) s
    type(real_8) cos,sin,theta

    call alloc(s)
    call alloc(cos,sin,theta)

    call find_axis(ds,y_axis)
    call find_perp_basis(y_axis,x_axis,z_axis)

    s=ds*x_axis

    cos=s.dot.x_axis
    sin=-(s.dot.z_axis)

    theta=clockwise*atan2(sin,cos)
    !if(force_positive.and.theta<zero) theta = theta + twopi    !!!! allow negative theta

    h_axis=(clockwise*theta)*y_axis

    call kill(cos,sin,theta)
    call kill(s)
  end subroutine find_exponentp

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
    n1=zero
    n1(is)=one

    s=n2(is)*n1(is)

    n=zero
    do i=1,3
       n1(i)=n1(i)-s*n2(i)
       n=n1(i)**2+n
    enddo
    n1=n1/sqrt(n)

    n3(1)=n1(2)*n2(3)-n1(3)*n2(2)
    n3(2)=n1(3)*n2(1)-n1(1)*n2(3)
    n3(3)=n1(1)*n2(2)-n1(2)*n2(1)

    n=zero
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
       n1(i)=zero
    enddo
    n1(is)=one

    s=n2(is)*n1(is)

    n=zero
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

    n=zero
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



  subroutine alloc_33t(a)
    implicit none
    type(taylor) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call alloc(a(i,j))
       enddo
    enddo

  end subroutine alloc_33t

  subroutine kill_33t(a)
    implicit none
    type(taylor) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call kill(a(i,j))
       enddo
    enddo

  end subroutine kill_33t

  subroutine zero_33t(a)
    implicit none
    type(real_8) a(3,3) ! taylor
    !    type(taylor) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          a(i,j)=zero
       enddo
    enddo

  end subroutine zero_33t

  subroutine alloc_33p(a)
    implicit none
    type(real_8) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call alloc(a(i,j))
       enddo
    enddo

  end subroutine alloc_33p

  subroutine kill_33p(a)
    implicit none
    type(real_8) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call kill(a(i,j))
       enddo
    enddo

  end subroutine kill_33p

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

       call smatp(-half,mt,mt)

       call smatmulp(one,mt,id,mt)

       call matmulp(mt,m,m)

       call transpose_p(m,mt)
       call matmulp(mt,m,mt)
       call smatmulp(-one/1.5e0_dp,id,mt,mt)
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

    r=zero
    rt=zero
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

    r=zero

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

    call smatp(one,m,mt)

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
          CALL ALLOC(D%S(I,J))
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
          CALL KILL(D%S(I,J))
       ENDDO
    ENDDO

  END    subroutine KILL_DASPIN




  subroutine ALLOC_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S

    CALL ALLOC(S%X,3)
    !     S%G=A_PARTICLE
  END    subroutine ALLOC_SPINOR_8

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


  END    subroutine ALLOC_probe_8

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


    !  DO I=1,6
    !  DO J=1,6
    !   CALL KILL(R%E_IJ(I,J))
    !  ENDDO
    !  ENDDO

  END    subroutine KILL_probe_8




!!!!!!!!!!!!!!!!   new stuff
  subroutine get_spin_nx_spinor_8(DS,theta0,n0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    type(real_8), intent(inout) :: theta0
    type(spinor_8), intent(inout) :: n0

    call get_spin_n0(DS,theta0,n0%x)

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

    call find_n0(ds%s,n0)

    !     call print(theta0,6)


    call find_a(n0,a)

    call inv_as(a,ai)



    call matmulp(ai,ds%s,ai)    !
    call matmulp(ai,a,s)    !

    a11=(s(1,1))
    a13=(s(1,3))

    theta0=clockwise*atan2(a13,a11)
    if(force_positive.and.theta0<zero) theta0 = theta0 + twopi    !!!! allow negative theta0

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
    if(force_positive.and.theta0<zero)  theta0 = theta0 + twopi   !!!! allow negative theta0


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
          s(i,j)=ds%s(i,j)
       enddo
    enddo

    call find_n0(s,n0)

    call find_a(n0,a)
    call inv_as(a,ai)
    ai=matmul(ai,s)    !
    s=matmul(ai,a)    !

    theta0=clockwise*atan2(s(1,3),s(1,1))
    if(force_positive.and.theta0<zero)  theta0 = theta0 + twopi   !!!! allow negative theta0



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
       call trans_mat(s11%s,s11%m,s11%s)
       call inv_as(s11%s,s11%s)

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

    dot_spinor=ZERO

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

    dot_spinor_8=ZERO

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
          t2%s(i,j)=s2%s(i,j)*s1%m
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
             s(i,j)= t2%s(i,k)*s1%s(k,j)+ s(i,j)
          enddo
       enddo
    enddo
    call smatp(one,s,concat%s)
    !    concat%s0=s1%s0

    call kill_33(s);
    call kill(t2);
    master=localmaster

  END FUNCTION concat

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



    do i=1,3
       do j=1,3
          cmul%s(i,j)=s2*s1%s(i,j)
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
    call trans_mat( mmul%s,s2, mmul%s)

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

    damapspin_spinor8_mul%x(1)=zero
    damapspin_spinor8_mul%x(2)=zero
    damapspin_spinor8_mul%x(3)=zero

    do i=1,3
       do j=1,3
          damapspin_spinor8_mul%x(i)=(s1%s(i,j))*S2%x(j)+damapspin_spinor8_mul%x(i)
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

    damapspin_spinor_mul%x(1)=zero
    damapspin_spinor_mul%x(2)=zero
    damapspin_spinor_mul%x(3)=zero

    do i=1,3
       do j=1,3
          damapspin_spinor_mul%x(i)=(s1%s(i,j))*S2%x(j)+damapspin_spinor_mul%x(i)
       enddo
    enddo

    master=localmaster


  END FUNCTION damapspin_spinor_mul

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

    spin8_mul_map%x(1)=zero
    spin8_mul_map%x(2)=zero
    spin8_mul_map%x(3)=zero

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

    spin8_scal8_map%x(1)=zero
    spin8_scal8_map%x(2)=zero
    spin8_scal8_map%x(3)=zero

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

    add_spin8_spin8%x(1)=zero
    add_spin8_spin8%x(2)=zero
    add_spin8_spin8%x(3)=zero

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

    mul_spin8_spin8%x(1)=zero
    mul_spin8_spin8%x(2)=zero
    mul_spin8_spin8%x(3)=zero

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

    sub_spin8_spin8%x(1)=zero
    sub_spin8_spin8%x(2)=zero
    sub_spin8_spin8%x(3)=zero

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
          addm%s(i,j)=s2%s(i,j)+s1%s(i,j)
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

    purge_transverse=one


    do i=1,c_%nd2

       if(i/=c_%ndpt) then
          if(j(i)/=0) then
             purge_transverse=zero
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
    call alloc(t)
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

  end  subroutine clean_orbital_33

  subroutine clean_orbital_3(s,sf)
    implicit none
    TYPE (taylor), intent(inout) :: s(3),sf(3)
    integer i

    do i=1,3
       call cfu(s(i),purge_transverse,sf(i))
    enddo

  end  subroutine clean_orbital_3

  subroutine clean_orbital(s,sf)
    implicit none
    TYPE (taylor), intent(inout) :: s,sf

    call cfu(s,purge_transverse,sf)

  end  subroutine clean_orbital

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
  subroutine Go_to_closed(ns,DS)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: ns
    TYPE(damapspin), INTENT(INout) :: DS
    type(damapspin) a1i,ds0
    type(damapspin)s ,a ,ai
    type(taylor) nn
    !    type(taylor) s(3,3),a(3,3),ai(3,3)
    type(real_8) a11,a13
    integer i,j
    real(dp) ss(3,3)

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


    s=1
    a=1
    ai=1

    call clean_orbital_33(ds0%s,s%s)
    ! at this stage s is the spin map  without dependence on orbital

    call find_n0(s%s,ns%n0)

    call find_a(ns%n0,a%s)

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

    a11=s%s(1,1)
    a13=s%s(1,3)

    ns%theta0=clockwise*atan2(a13,a11)
    if(force_positive.and.ns%theta0<zero)  ns%theta0 = ns%theta0 + twopi  !!!! allow negative theta0

    ns%as=ns%as*a

    ns%ar=1
    ns%ar%m=ns%as%m**(-1)*ns%n%a_t

!!!!   as=(A1,A_spin)
!!!!   ar=(Ar,I)

    call kill(a1i)
    call kill(ds0)
    call kill(s)
    call kill(a)
    call kill(ai)

    call kill(a11,a13)
    call kill(nn)

  end subroutine Go_to_closed

  subroutine fetch_s0(DS,s0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(damapspin), INTENT(INout) :: s0
    integer i,j

    s0=1

    do i=1,3
       do j=1,3
          s0%s(i,j)=ds%s(i,j).sub.'0'
       enddo
    enddo



  end subroutine fetch_s0

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
    dh%s(1,1)=dh%s(1,1)-one
    dh%s(2,2)=dh%s(2,2)-one
    dh%s(3,3)=dh%s(3,3)-one

    c=one
    do i=1,c_%no
       dhn=dhn*dh
       cl=c/i
       h=h+cl*dhn
       c=-c
    enddo

    om(1)=h%s(3,2)
    om(2)=h%s(1,3)
    om(3)=h%s(2,1)

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
    dh%s(2,1)=om(3)
    dh%s(1,3)=om(2)
    dh%s(3,2)=om(1)
    dh%s(1,2)=-om(3)
    dh%s(3,1)=-om(2)
    dh%s(2,3)=-om(1)
    dhn=1
    c=one
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
    real(dp) value,ang
    complex(dp) denom
    logical doit

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
          ang=zero
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-one)
          omr(2)=omr(2)+ (denom.mono.jc)
       endif
    enddo

    !  sin  part of om(2)
    call taylor_cycle(t%sin,size=N)

    do i=1,n
       call taylor_cycle(t%sin,ii=i,value=value,j=jc)

       call test_jc(ns,jc,nd,doit)

       if(doit) then
          ang=zero
          do j=1,nd
             ang=(jc(j*2-1)-jc(j*2))*twopi*ns%n%tune(j)+ang
          enddo
          denom=value/(exp(-i_*ang)-one)
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
          ang=zero
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
          ang=zero
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
          ang=zero
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
          ang=zero
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
          ang=zero
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
          ang=zero
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
          ang=zero
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
          ang=zero
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

    omc(1)= (omr(3)+omr(1))/two
    omc(3)= (omr(3)-omr(1))/two/i_
    omc(2)=omr(2)

    do i=1,3
       t%cos=omc(i)%r
       t%sin=omc(i)%i
       oma(i)=t
    enddo

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
    if(.not.doit) then
       write(6,*) jc(1:4),is
       write(6,*) ns%m(1,1),ns%m(2,1),ns%ms(1)
    endif

    !    write(6,*) doit

    !    pause 313

  end subroutine test_jc_spin

  subroutine normalise_spin(ns,DS_in)
    implicit none
    TYPE(normal_spin), INTENT(INout) :: ns
    TYPE(damapspin), INTENT(IN) :: DS_in
    TYPE(damapspin)  ds,ds0,dst,dsi
    TYPE(damap)  r0
    type(taylor) om(3),oma(3)
    integer i

    call alloc(ds)
    call alloc(dst)
    call alloc(dsi)
    call alloc(ds0)
    call alloc(r0)
    call alloc(om)
    call alloc(oma)

    ds=ds_in
    ! step 1
    !  Normalize to the parameter dependent n0 of the theory
    call Go_to_closed(ns,DS)

    ! step 2
    ! Apply the transformation found in step 1 to ds
    ds=ns%as**(-1)*ds*ns%as
    ! Apply the transformation found in step 1 to ds
    ds=ns%ar**(-1)*ds*ns%ar

    call fetch_s0(DS,ds0)


    ns%a_t=1
    r0=ns%n%normal

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
       x(i)=zero
    enddo

    a_temp=a_in

    do i=1,3
       do j=1,3
          if(a_in%s(i,j)%kind==2) then
             a_temp%s(i,j)=a_in%s(i,j)%t*x
          endif
       enddo
    enddo
    a_out=a_temp

    call kill(a_temp)

  end subroutine eval_spin_matrix


  subroutine factor_as(a_t,a_f,a_s,a_l,a_nl)
    implicit none
    TYPE(damapspin), INTENT(INout) :: a_t,a_f,a_s,a_l,a_nl
    a_f=1
    a_l=1
    a_nl=1

    call factor(a_t%m,a_f%m,a_l%m,a_nl%m)

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




  subroutine factor_am(a_t,a_f,a_l,a_nl)
    implicit none
    TYPE(damap), INTENT(INout) :: a_t,a_nl,a_l,a_f
    integer i,n,k
    integer, allocatable :: jc(:)
    real(dp) value
    logical doit
    TYPE(damap) atemp

    call alloc(atemp)
    allocate(jc(c_%nv))

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
       atemp%v(c_%nd2-1)=atemp%v(c_%nd2-1)-(one.mono.(c_%nd2-1))
       atemp%v(c_%nd2)=atemp%v(c_%nd2)-(one.mono.(c_%nd2))
       do i=1,c_%nd-1
          atemp%v(c_%ndpt+1)=atemp%v(c_%ndpt+1) &
               +(atemp%v(2*i).d.c_%ndpt)*(one.mono.2*i-1)-(atemp%v(2*i-1).d.c_%ndpt)*(one.mono.2*i)
       enddo
    elseif(c_%ndpt==c_%nd2) then      ! Marylie convention
       atemp%v(c_%nd2-1)=atemp%v(c_%nd2-1)-(one.mono.(c_%nd2-1))
       atemp%v(c_%nd2)=atemp%v(c_%nd2)-(one.mono.(c_%nd2))
       do i=1,c_%nd-1
          atemp%v(c_%ndpt-1)=atemp%v(c_%ndpt-1) &
               -(atemp%v(2*i).d.c_%ndpt)*(one.mono.2*i-1)+(atemp%v(2*i-1).d.c_%ndpt)*(one.mono.2*i)
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
       atemp%v(c_%nd2-1)=(one.mono.(c_%nd2-1))
       atemp%v(c_%nd2)= (one.mono.(c_%nd2))

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


    deallocate(jc)
    call kill(atemp)


  end subroutine factor_am

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

    h%h=muc*((1.d0.mono.'002')+(1.d0.mono.'0002'))/two
    h%h=h%h+muc*((1.d0.mono.'2')+(1.d0.mono.'02'))/two

    nc%m=texp(h,nc%m)


    !     call print(nc%m,6)
    !     nf=nc%m

    uno=nc%m
    uno%pb%h=uno%pb%h-muc*((1.d0.mono.'2')+(1.d0.mono.'02'))/two
    tr=uno%pb%h

    call print(ns%theta0,6)

    !    a11=(s(1,1))
    !    a13=(s(1,3))

    theta0=atan2(nc%s(1,3),nc%s(1,1))
    if((theta0.sub.'0')<zero) theta0 = theta0 + twopi

    theta0r=theta0

    call taylor_clean(theta0r%cos,1.d-1)
    call taylor_clean(theta0r%sin,1.d-1)
    call taylor_clean(tr%cos,1.d-5)
    call taylor_clean(tr%sin,1.d-5)

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
          call clean_real_8(s1%s(i,j),s2%s(i,j),prec)
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



end module tree_element_MODULE
