!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tree_element_MODULE
  USE polymorphic_complextaylor, mon1=>mon
  IMPLICIT NONE
  public
  integer,private,parameter::ndd=6

  PRIVATE track_TREE,track_TREEP,KILL_TREE,KILL_TREE_N   !,SET_TREE
  PRIVATE track_TREE_G,track_TREEP_g
  PRIVATE ALLOC_SPINOR_8,ALLOC_probe_8,r_AVERAGE
  PRIVATE KILL_SPINOR_8,KILL_probe_8
  PRIVATE EQUAL_SPINOR8_SPINOR8,EQUAL_IDENTITY_SPINOR_8 !,EQUAL_SPINOR8_RAY8,EQUAL_RAY8_SPINOR8,
  PRIVATE  EQUAL_IDENTITY_probe_8
  private EQUAL_IDENTITY_probe
  private EQUAL_IDENTITY_SPINOR,EQUAL_PROBE_REAL6
  PRIVATE EQUAL_SPINOR8_SPINOR,EQUAL_PROBE8_PROBE,EQUAL_PROBE8_REAL6
  private EQUAL_IDENTITY_SPINOR_8_r3 ,EQUAL_SPINOR_SPINOR8
  private ALLOC_rf_phasor_8,KILL_rf_phasor_8,realdp_spinor,cross_real
  private sub_spinor
  private EQUAL_PROBE_PROBE
  private dot_spinor_8,dot_spinor,dot_real
  private  read_spinor_8
  !  private smatp,smatmulp

  PRIVATE EQUAL_PROBE8_PROBE8,PRINT_probe8,PRINT_probe

  private read_probe8


  private scdaddo,daddsco
  private real_8REAL6,REAL6real_8,real_8REAL_8
  private probe_quaternion_to_matrixr,probe_quaternion_to_matrixp


  private EQUAL_RF8_RF8 !,extract_envelope_probe8
  PRIVATE EQUAL_RF8_RF,EQUAL_RF_RF8,print_rf_phasor_8 !,extract_envelope_damap
  private EQUAL_DAMAP_RAY8,cross_spinor,cross_spinor8
  private flip  ! flip in lielib
  integer, target :: spin_extra_tpsa = 0 ,n0_normal= 2
  logical(lp) :: force_positive=.false.
  logical(lp) :: use_ptc_ac_position=.false.
  integer, private :: nd_used,i_phase,i_plane, decal
  logical :: onelie = my_false, phase_nonlinear=my_true
  logical :: firstfac=.true.
  integer, private, parameter :: nfac=20
  real(dp), private :: fac(0:nfac)
  integer :: nbe=8
  integer :: n_rf=0  !number of modulation clocks in the simulation
  integer :: modulationtype=0 ! 0 is the full blown and internal anf externa field, 1 is simple one on external field only without cos(theta)
  
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
     MODULE PROCEDURE EQUAL_DAMAP_RAY8
     MODULE PROCEDURE EQUAL_PROBE8_PROBE8
     MODULE PROCEDURE     EQUAL_SPINOR_SPINOR8
     MODULE PROCEDURE     EQUAL_PROBE_PROBE8
     MODULE PROCEDURE     EQUAL_PROBE_PROBE


     MODULE PROCEDURE     EQUAL_RF8_RF8
     MODULE PROCEDURE     EQUAL_RF8_RF
     MODULE PROCEDURE     EQUAL_RF_RF8
     !    MODULE PROCEDURE     EQUAL_nn_damap
  end  INTERFACE

  INTERFACE probe_quaternion_to_matrix
     MODULE PROCEDURE probe_quaternion_to_matrixr
     MODULE PROCEDURE probe_quaternion_to_matrixp
  END INTERFACE
 
  INTERFACE OPERATOR (.dot.)
     MODULE PROCEDURE dot_real
     MODULE PROCEDURE dot_spinor
     MODULE PROCEDURE dot_spinor_8
  END  INTERFACE

  INTERFACE OPERATOR (.cross.)
     MODULE PROCEDURE cross_real
  END  INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE cross_spinor
     MODULE PROCEDURE realdp_spinor
     MODULE PROCEDURE cross_spinor8
  END  INTERFACE

  INTERFACE operator (+)
     MODULE PROCEDURE scdaddo
     MODULE PROCEDURE daddsco
  END  INTERFACE

  INTERFACE operator (-)
     MODULE PROCEDURE sub_spinor
  END  INTERFACE


  INTERFACE PRINT

!!!
     MODULE PROCEDURE PRINT_probe
     MODULE PROCEDURE PRINT_probe8
     MODULE PROCEDURE PRINT_spinor_8
     MODULE PROCEDURE print_rf_phasor_8
  END INTERFACE



  INTERFACE READ
     MODULE PROCEDURE read_probe8   ! a bit illegal : reading polymorphs as taylor...
     MODULE PROCEDURE read_spinor_8
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_SPINOR_8
     MODULE PROCEDURE ALLOC_probe_8
     MODULE PROCEDURE ALLOC_rf_phasor_8
     MODULE PROCEDURE SET_TREE
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILL_SPINOR_8
     MODULE PROCEDURE KILL_probe_8
     MODULE PROCEDURE KILL_rf_phasor_8
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

  INTERFACE AVERAGE
     MODULE PROCEDURE r_AVERAGE
  END INTERFACE



CONTAINS
 

 
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
       if(nd2<=4.and.(c_%npara==3.or.c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
          If(ndpt_bmad==0) then
           if(nd2==4) daddsco(i)=s2(i)+(1.0_dp.mono.'00001')
           if(nd2==2) daddsco(i)=s2(i)+(1.0_dp.mono.'001')
          endif
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
       if(nd2<=4.and.(c_%npara==3.or.c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
         if(nd2==4)  scdaddo(i)=s2(i)+(1.0_dp.mono.'00001')
         if(nd2==2)  scdaddo(i)=s2(i)+(1.0_dp.mono.'001')
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
    U%NP=T%NP
    U%no=T%no
    U%FIXr=T%FIXr
    U%ds=T%ds
    U%beta0=T%beta0
    U%FIX=T%FIX
    U%FIX0=T%FIX0
    U%e_ij=T%e_ij
    U%rad=T%rad
!    U%file=T%file
    U%eps=T%eps
    U%symptrack=T%symptrack
    U%usenonsymp=T%usenonsymp
    U%factored=T%factored
   ! U%ng=T%ng
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

    NULLIFY(T%CC,T%JL,T%JV,T%N,T%NP,T%no,t%fixr,t%fix,t%fix0,t%beta0,t%e_ij,t%rad,t%ds)  !,t%file)

  END SUBROUTINE NULL_TREE


  SUBROUTINE ALLOC_TREE(T,N,np)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T
    INTEGER , INTENT(IN) :: N,np
    integer i
    !IF(N==0) RETURN


    ALLOCATE(T%CC(N),T%fix0(6),T%fix(6),T%fixr(6),T%JL(N),T%JV(N),T%N,T%ds,T%beta0,T%np,T%no, & 
  !  t%e_ij(c_%nd2,c_%nd2),T%rad(c_%nd2,c_%nd2),t%usenonsymp, t%symptrack, t%eps)  !,t%file)
     t%e_ij(6,6),T%rad(6,6),t%usenonsymp, t%symptrack, t%eps,t%factored) !,t%ng)  !,t%file)
    t%cc=0
    t%jl=0
    t%jv=0
    T%N=N
    T%np=np
    T%no=0
    T%fix=0.0_dp
    T%fix0=0.0_dp
    T%fixr=0.0_dp
    T%e_ij=0.0_dp
    T%ds=0.0_dp
    T%beta0=0.0_dp
    T%rad=0.0_dp
!    T%file=' '
    do i=1,6
     T%rad(i,i)=1.0_dp
    enddo
    t%eps=1.d-7
    t%symptrack=.false.
    t%usenonsymp=.false.
    t%factored=.false.
  !  t%ng=1
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

    do i=1,T%np
       xt(i)=xi(i)
    enddo
    do i=1,T%np
       xf(i) = T%cc(i)
    enddo

    XM(1) = 1.0_dp
    JC=T%np
    do i=1,(T%N-T%np)/T%np
       !
       xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%np
          jc=jc+1
          xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo
    do i=1,size(xi)
       xI(i)=xF(i)
    enddo

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




    do i=1,T%np
       xt(i)=xi(i)
    enddo
    do i=1,T%np
       xf(i) = T%cc(i)
    enddo

    XM(1) = 1.0_dp
    JC=T%np

    do i=1,(T%N-T%np)/T%np
       !
       xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%np
          jc=jc+1
          xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo

    do i=1,size(xi)
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


     IF(ASSOCIATED(T%CC))DEALLOCATE(T%CC,T%fix0,T%fix,T%fixr,t%ds,t%beta0,T%JL,T%JV,T%N,T%NP, &
    T%No,t%e_ij,t%rad,t%eps,t%symptrack,t%usenonsymp,t%factored)  !,t%file)


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
  !  xi(1:t%np)=xi(1:t%np)-t%fix
    n1=1
    if(present(n)) n1=n
    do k=1,n1
       if(.not.c_%CHECK_STABLE) return
       XT=0.0_dp
       XF=0.0_dp
       XM=0.0_dp

       do i=1,T%np
          xt(i)=xi(i)
       enddo
       do i=1,T%np
          xf(i) = T%cc(i)
       enddo

       XM(1) = 1.0_dp
       JC=T%np
       do i=1,(T%N-T%np)/T%np
          !
          xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
          xm(T%jl(JC+1)+1) = xx
          !
          do iv=1,T%np
             jc=jc+1
             xf(iv) = xf(iv) + t%cc(jc) * xx
          enddo
       enddo
       xi=xf

       if(abs(xi(1))>c_%absolute_aperture.or.abs(xi(3))>c_%absolute_aperture) then
          c_%CHECK_STABLE=.FALSE.
          xlost=xi
          messagelost="o_tree_element.f90 track_tree : aperture exeeded"
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




       do i=1,T%np
          xt(i)=xi(i)
       enddo
       do i=1,T%np
          xf(i) = T%cc(i)
       enddo

       XM(1) = 1.0_dp
       JC=T%np

       do i=1,(T%N-T%np)/T%np
          !
          xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
          xm(T%jl(JC+1)+1) = xx
          !
          do iv=1,T%np
             jc=jc+1
             xf(iv) = xf(iv) + t%cc(jc) * xx
          enddo
       enddo

       do i=1,T%np
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
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER I
    P%u=my_false
    !       P%s(0)%x=zero
    !       P%s(0)%x(n0_normal)=one

    DO I=1,ISPIN1R
       P%s(i)%x=0.0_dp
       P%s(i)%x(i)=1.0_dp
    enddo
! quaternion
    p%q=1.0_dp
    p%use_q=use_quaternion
    DO I=1,size(x,1)
       P%X(i)=X(i)
    enddo
    P%ac%t=0.0_dp
    p%nac=n_rf
  END    subroutine EQUAL_PROBE_REAL6

  subroutine EQUAL_PROBE8_REAL6 (P,X)
    implicit none
    TYPE(PROBE_8), INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: X(:)
    INTEGER I

    P%u=my_false
    !    P%S=0

    !       P%s(0)=0
    !       P%s(0)%x(n0_normal)=one

    DO I=1,3
       P%s(i)=0
       P%s(i)%x(i)=1.0_dp
    enddo
! quaternion
    p%q=1.0_dp
    DO I=1,size(x,1)
       P%X(i)=X(i)
    enddo
    p%use_q=use_quaternion
    p%nac=n_rf
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
! quaternion
    p8%q=p%q
    !    P8%S=P%S
    DO I=1,3
       P8%S(i)=P%S(i)
    enddo
     p8%q=p%q
     P8%nac=P%nac
    do i=1,P%nac
    P8%AC(i)=P%AC(i)
    enddo
    P8%u=P%u
    p8%use_q=P%use_q
    P8%e=P%e
    P8%x0=P%x0


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
! quaternion
    P8%q=P%q
    P8%nac=P%nac
!!!new 2018.1.4
    do I=1,P%nac
       P8%ac(i)=P%ac(i)
    enddo

!!!!
    P8%e=P%e
    P8%u=P%u
    P8%e_ij=0.0_dp
    p8%use_q=P%use_q
    p8%x0=P%x
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
    P%e=P8%e
! quaternion
    P%q=P8%q
    p%use_q=P8%use_q

!!!new 2018.1.4
    p%nac=P8%nac
    do I=1,P8%nac
       P%ac(i)=P8%ac(i)
    enddo
  END subroutine EQUAL_PROBE_PROBE8

  subroutine EQUAL_PROBE_PROBE (P,P8)
    implicit none
    TYPE(PROBE), INTENT(INOUT) :: P
    TYPE(PROBE), INTENT(IN) :: P8
    INTEGER I
    DO I=1,6
       P%X(I)=P8%X(I)
    ENDDO
    !     P%S(0)=P8%S(0)
    DO I=1,ISPIN1R
       P%S(I)=P8%S(I)
    ENDDO
    P%u=P8%u
    P%e=P8%e
! quaternion
    P%q=P8%q
    p%use_q=P8%use_q
!!!new 2018.1.4
    p%nac=P8%nac
    do I=1,P8%nac
       P%ac(i)=P8%ac(i)
    enddo

  END subroutine EQUAL_PROBE_PROBE

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
       R%S(1)=0
       R%S(2)=0
       R%S(3)=0
    else
       STOP 100
    ENDIF
! quaternion
    r%q=1.0_dp
    r%u=.false.
    r%use_q=use_quaternion
    r%e=0
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
       r%q=1.0_dp
    ELSEIF(S==0) THEN

       !       DO I=1,6  !-C_%NSPIN
       !        R%X(I)=ZERO!
       !       enddo
       r%q=0.0_dp
    ELSE
       STOP 100
    ENDIF
    R%e_ij=0.0_dp
    r%u=.false.
    r%use_q=use_quaternion
    r%e=0
    r%x0=0
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
    INTEGER I,nd2t,j

    nd2t=C_%ND2
    if(doing_ac_modulation_in_ptc) then
       nd2t=C_%ND2-2
    endif



    DO I=1,nd2t
       DS%V(I)=R%X(I)
    ENDDO
    j=1
    DO I=nd2t+1,C_%ND2,2
       DS%V(I)=R%ac(j)%x(1)
       DS%V(I+1)=R%ac(j)%x(2)
       j=j+1
    ENDDO

  END subroutine EQUAL_DAMAP_RAY8
 

   subroutine  probe_quaternion_to_matrixr(p)
    implicit none
    TYPE(probe), INTENT(INOUT) :: p
    type(quaternion) s,sf
    integer i
    do i=1,3
     s=0.0_dp
     s%x(i)=1.0_dp
     sf=p%q*s*p%q**(-1)
     p%s(i)%x=sf%x(1:3)
    enddo

    end subroutine  probe_quaternion_to_matrixr 

   subroutine  probe_quaternion_to_matrixp(p)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: p
    type(quaternion_8) s,sf
    integer i,j
     call ALLOC(s)
     call ALLOC(sf)
    do i=1,3
     s=0.0_dp
     s%x(i)=1.0_dp
     sf=p%q*s*p%q**(-1)
     do j=1,3
       p%s(i)%x(j)=sf%x(j+1)
     enddo
    enddo
     call kill(s)
     call kill(sf)
    end subroutine  probe_quaternion_to_matrixp 

   subroutine print_probe(DS,MF)
    implicit none
    TYPE(probe), INTENT(INOUT) :: DS
    INTEGER I,mfi
     integer,optional :: mf

     mfi=6
     if(present(mf)) mfi=mf

    WRITE(MFi,*) " ORBIT "
    do i=1,6
       write(mfi,*) ' Variable ',i
       write(mfi,'(6(1X,G20.13))') ds%x(i) 
    enddo
   if(ds%use_q) then
    WRITE(MFi,*) " quaternion "
     call print(ds%q,mfi)
    else
    WRITE(MFi,*) " SPIN X "
       write(mfi,'(3(1X,G20.13))') ds%s(1)%x 
 
    WRITE(MFi,*) " SPIN Y "
       write(mfi,'(3(1X,G20.13))') ds%s(2)%x 
 
    WRITE(MFi,*) " SPIN Z "
       write(mfi,'(3(1X,G20.13))') ds%s(3)%x 
   endif
 

  END subroutine print_probe

 
  subroutine print_probe8(DS,MF)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: DS
    INTEGER MFi,I,j
    logical(lp) rad_in
    integer,optional :: mf
     mfi=6
     if(present(mf)) mfi=mf

    WRITE(MFi,*) " ORBIT "
    do i=1,6
       write(mfi,*) ' Variable ',i
       call print(ds%x(i),mfi)
    enddo

   if(ds%use_q) then
    WRITE(MFi,*) " quaternion "
     call print(ds%q,mfi)
    else

    WRITE(MFi,*) " SPIN X "
    call print(ds%s(1),mfi)
 
    WRITE(MFi,*) " SPIN Y "
    call print(ds%s(2),mfi)
 
    WRITE(MFi,*) " SPIN Z "
    call print(ds%s(3),mfi)
 
    endif

    call check_rad(DS%e_ij,rad_in)



    if(rad_in) then

       WRITE(MFi,*) " STOCHASTIC KICK  "

       do i=1,6
          do j=1,6
             write(mfi,*) i,j,ds%e_ij(i,j)
          enddo
       enddo
    else

       WRITE(MFi,*) "NO STOCHASTIC KICK  "

    endif
    if(doing_ac_modulation_in_ptc) then
      write(mfi,*) ds%nac, " clocks "
       do i=1,ds%nac
        call print(ds%ac(i),mfi)
       enddo
    else
       WRITE(MFi,*) "NO MODULATION  "

    endif

  END subroutine print_probe8

  subroutine print_rf_phasor_8(S,MF)
    implicit none
    TYPE(rf_phasor_8), INTENT(INOUT) :: s
    INTEGER MFi,I
    integer,optional :: mf
     mfi=6
     if(present(mf)) mfi=mf

    write(mfi,*) ' AC INFORMATION : omega, pseudo-time, hands of the clock'
    call print(s%om,mfi)
    call print(s%t,mfi)
    do i=1,2
       call print(s%x(i),mfi)
    enddo

  END subroutine print_rf_phasor_8


  subroutine print_spinor_8(S,MF)
    implicit none
    TYPE(spinor_8), INTENT(INOUT) :: s
    INTEGER MFi,I
    integer,optional :: mf
     mfi=6
     if(present(mf)) mfi=mf

    do i=1,3
       write(mfi,*) ' Spin Variable ',i
       call print(s%x(i),mfi)
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
       read(mf,'(a255)') line
       call read(t,mf)
       s%x(i)=morph(t)
    enddo

    call kill(t)

  END subroutine read_spinor_8

  

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



!



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


 



  subroutine ALLOC_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S

    CALL ALLOC(S%X,3)
    !     S%G=A_PARTICLE
  END    subroutine ALLOC_SPINOR_8


  subroutine ALLOC_probe_8(R,n)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    INTEGER I,nac
    integer, optional :: n
     nac=1
     if(present(n)) nac=n
    !    CALL ALLOC(R%S)
    DO I=1,3
       CALL ALLOC(R%S(I))
    ENDDO
    CALL ALLOC(R%X,6)
    CALL ALLOC(R%q)
    !      R%S(0)%X(N0_NORMAL)=ONE
    DO I=1,3
       R%S(I)=0
    ENDDO
    r%nac=nac
    do i=1,nac
     CALL ALLOC(R%ac(i))
    enddo
    r%e_ij=0.0_dp
    r%u=.false.
    r%use_q=use_quaternion
    r%e=0
    r%x0=0
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
    CALL KILL(R%q)
    do i=1,r%nac
     CALL KILL(R%ac(i))
    enddo
    r%nac=0
    r%e_ij=0.0_dp
    r%u=.false.
    r%e=0
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
  

  FUNCTION dot_real( S1, S2 )
    implicit none
    real(dp) dot_real
    real(dp), INTENT (IN) :: S1(:),S2(:)

    INTEGER I

    dot_real=0.0_dp
   
    DO I=1,min(size(s1),size(s2))
       dot_real=dot_real+s1(i)*s2(i)
    ENDDO

  END FUNCTION dot_real



  subroutine make_spinor_basis(s1,s2,s3)
  implicit none
    TYPE (SPINOR) s1,s2,s3
    
     s1=(1.0_dp/sqrt(s1.dot.s1))*s1

     s2=s2 - (s2.dot.s1)*s1
     
     s2=(1.0_dp/sqrt(s2.dot.s2))*s2

     s3=s1*s2

  end subroutine make_spinor_basis


  FUNCTION realdp_spinor( S1, S2 )
    implicit none
    TYPE (SPINOR) realdp_spinor
    real(dp), intent(in) :: s1
    TYPE (SPINOR), INTENT (IN) :: S2

 
 

       realdp_spinor%x(1)= s1*s2%x(1)
       realdp_spinor%x(2)= s1*s2%x(2)
       realdp_spinor%x(3)= s1*s2%x(3)
    

  END FUNCTION realdp_spinor


  FUNCTION sub_spinor( S1, S2 )
    implicit none
    TYPE (SPINOR) sub_spinor
    TYPE (SPINOR), INTENT (IN) :: S1,S2

 
 

       sub_spinor%x(1)= s1%x(1)-s2%x(1)
       sub_spinor%x(2)= s1%x(2)-s2%x(2)
       sub_spinor%x(3)= s1%x(3)-s2%x(3)
    

  END FUNCTION sub_spinor

  FUNCTION cross_spinor( S1, S2 )
    implicit none
    TYPE (SPINOR) cross_spinor
    TYPE (SPINOR), INTENT (IN) :: S1,S2

       cross_spinor%x(1)= s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
       cross_spinor%x(2)= -s1%x(1)*s2%x(3)+s1%x(3)*s2%x(1)
       cross_spinor%x(3)= s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)
    
  END FUNCTION cross_spinor

  FUNCTION cross_real( S1, S2 )
    implicit none
    real(dp) cross_real(3)
    real(dp), INTENT (IN) :: S1(3),S2(3)

       cross_real(1)= s1(2)*s2(3)-s1(3)*s2(2)
       cross_real(2)=-s1(1)*s2(3)+s1(3)*s2(1)
       cross_real(3)= s1(1)*s2(2)-s1(2)*s2(1)
    
  END FUNCTION cross_real

  FUNCTION cross_spinor8( S1, S2 )
    implicit none
    TYPE (SPINOR_8) cross_spinor8
    TYPE (SPINOR_8), INTENT (IN) :: S1,S2
    integer localmaster

    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call ass(cross_spinor8%x(1))
    call ass(cross_spinor8%x(2))
    call ass(cross_spinor8%x(3))
 

       cross_spinor8%x(1)= s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
       cross_spinor8%x(2)= -s1%x(1)*s2%x(3)+s1%x(3)*s2%x(1)
       cross_spinor8%x(3)= s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)
    

    master=localmaster

  END FUNCTION cross_spinor8

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

 
!!! Some useful routines

  subroutine r_AVERAGE(F,A,F_floquet,F_xp,use_J)
    implicit none
    type(damap) A
    TYPE(TAYLOR), intent(inout):: F
    TYPE(TAYLOR), intent(inout):: F_floquet
    TYPE(TAYLOR), optional, intent(inout):: F_xp
    type(taylor) tt,tt_xp
    TYPE(taylorresonance) fq
    real(dp) value,valuexp
    integer, allocatable :: jc(:)
    logical(lp), optional :: use_J
    integer i,n,j,it,nd,iu
    logical uj



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
          tt_xp=(valuexp.mono.jc)+tt_xp
          if(uj) then
             value=valuexp*2.0_dp**iu
               do j=1,nd
                  jc(j*2)=0
               enddo           
          endif
          tt=(value.mono.jc)+tt

       endif

    enddo

  !  fq%cos=tt
    F_floquet=tt

    if(present(F_xp)) then
       fq%cos=tt_xp
       F_xp=fq
       F_xp=F_xp*A**(-1)
    endif

    deallocate(jc)
    call kill(tt,tt_xp)
    call kill(fq)


  end subroutine r_AVERAGE

  ! remove small numbers

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
    integer i,k 
    
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
