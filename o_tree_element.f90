!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tree_element_MODULE
  USE polymorphic_complextaylor
  IMPLICIT NONE
  public

  PRIVATE track_TREE,track_TREEP,KILL_TREE,KILL_TREE_N
  PRIVATE track_TREE_G,track_TREEP_g
  PRIVATE ALLOC_SPINOR_8,ALLOC_RAY_8
  PRIVATE KILL_SPINOR_8,KILL_RAY_8,KILL_DASPIN
  PRIVATE EQUAL_RAY8_SPINOR8,EQUAL_SPINOR8_SPINOR8,EQUAL_IDENTITY_SPINOR_8,EQUAL_SPINOR8_RAY8
  PRIVATE  EQUAL_IDENTITY_RAY_8,ALLOC_daspin,EQUAL_DASPIN_RAY8,EQUAL_RAY8_DASPIN
  private alloc_normal_spin,kill_normal_spin,EQUAL_NORMAL_DASPIN
  private READ_DASPIN,PRINT_DASPIN

  private power_rotr,find_ar,norm_matr,anti_matr
  private power_rotp,find_ap,norm_matp,anti_matp
  private  find_n_thetar,find_n_thetap
  !  private smatp,smatmulp
  private exp_n_thetar,exp_n_thetap,inv_asr,inv_asp !,inv_as


  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL_IDENTITY_RAY_8
     MODULE PROCEDURE EQUAL_IDENTITY_SPINOR_8
     MODULE PROCEDURE EQUAL_RAY8_SPINOR8
     MODULE PROCEDURE EQUAL_SPINOR8_RAY8
     MODULE PROCEDURE EQUAL_SPINOR8_SPINOR8
     MODULE PROCEDURE EQUAL_DASPIN_RAY8
     MODULE PROCEDURE EQUAL_RAY8_DASPIN
     MODULE PROCEDURE EQUAL_NORMAL_DASPIN
  end  INTERFACE

  INTERFACE exp_n_theta
     MODULE PROCEDURE exp_n_thetar
     MODULE PROCEDURE exp_n_thetap
  END INTERFACE

  INTERFACE find_n_theta
     MODULE PROCEDURE find_n_thetar
     MODULE PROCEDURE find_n_thetap
  END INTERFACE

  INTERFACE power_rot
     MODULE PROCEDURE power_rotr
     MODULE PROCEDURE power_rotp
  END INTERFACE

  INTERFACE find_a
     MODULE PROCEDURE find_ar
     MODULE PROCEDURE find_ap
  END INTERFACE

  INTERFACE norm_mat
     MODULE PROCEDURE norm_matr
     MODULE PROCEDURE norm_matp
  END INTERFACE

  INTERFACE anti_mat
     MODULE PROCEDURE anti_matr
     MODULE PROCEDURE anti_matp
  END INTERFACE

  INTERFACE inv_as
     MODULE PROCEDURE inv_asr
     MODULE PROCEDURE inv_asp
  END INTERFACE

  INTERFACE PRINT
     MODULE PROCEDURE PRINT_DASPIN
  END INTERFACE

  INTERFACE READ
     MODULE PROCEDURE READ_DASPIN
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_SPINOR_8
     MODULE PROCEDURE ALLOC_RAY_8
     MODULE PROCEDURE ALLOC_daspin
     MODULE PROCEDURE alloc_normal_spin
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILL_SPINOR_8
     MODULE PROCEDURE KILL_RAY_8
     MODULE PROCEDURE KILL_DASPIN
     MODULE PROCEDURE KILL_normal_spin
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




CONTAINS

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
    INTEGER JC,I,IV,ng

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
    integer i,k,nk,nvk,p,p0,p1,pg

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


  subroutine EQUAL_IDENTITY_SPINOR_8(S,R)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    INTEGER, INTENT(IN) :: R
    INTEGER I,j

    DO I=1,3
       DO J=1,3
          S%M(I,J)=ZERO
       ENDDO
    ENDDO
    IF(R==1) THEN

       DO I=1,C_%NSPIN
          S%M(I,I)=ONE
          S%X(I)=ONE.MONO.(c_%npara_fpp-C_%NSPIN+I)
       ENDDO
    ELSEIF(R==0) THEN
       DO I=1,C_%NSPIN
          S%X(I)=ZERO
       enddo
    ELSE
       STOP 100
    ENDIF

  END    subroutine EQUAL_IDENTITY_SPINOR_8

  subroutine EQUAL_IDENTITY_RAY_8(R,S)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    INTEGER, INTENT(IN) :: S
    INTEGER I

    R%S=S
    IF(S==1) THEN
       DO I=1,c_%npara_fpp-C_%NSPIN
          R%X(I)=ONE.MONO.I
       ENDDO
    ELSEIF(S==0) THEN
       R%X(I)=ZERO
    ELSE
       STOP 100
    ENDIF

  END    subroutine EQUAL_IDENTITY_RAY_8

  subroutine EQUAL_SPINOR8_RAY8(S,R)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    TYPE(probe_8), INTENT(IN) :: R

    S=R%S

  END    subroutine EQUAL_SPINOR8_RAY8

  subroutine EQUAL_RAY8_SPINOR8(R,S)
    implicit none
    TYPE(SPINOR_8), INTENT(IN) :: S
    TYPE(probe_8), INTENT(INOUT) :: R
    R%S=S

  END    subroutine EQUAL_RAY8_SPINOR8

  subroutine EQUAL_SPINOR8_SPINOR8(R,S)
    implicit none
    TYPE(SPINOR_8), INTENT(IN) :: S
    TYPE(SPINOR_8), INTENT(INOUT) :: R
    INTEGER I,J

    DO I=1,3
       R%X(I)=S%X(I)
       DO J=1,3
          R%M(I,J)=S%M(I,J)
       ENDDO
    ENDDO

  END    subroutine EQUAL_SPINOR8_SPINOR8

  subroutine EQUAL_DASPIN_RAY8(DS,R)
    implicit none
    TYPE(probe_8), INTENT(IN) :: R
    TYPE(DASPIN), INTENT(INOUT) :: DS
    INTEGER I,J

    DO I=1,6
       DS%X(I)=R%X(I)
    ENDDO

    DO I=1,C_%ND2
       DS%M%V(I)=R%X(I)
    ENDDO

    IF(C_%track_spint_mat) THEN
       DO I=1,3
          DO J=1,3
             DS%S(I,J)=R%S%M(I,J)
          ENDDO
       ENDDO
    ELSE
       DO I=1,3
          DO J=1,3
             DS%S(I,J)=(R%S%X(I)%T).D.(C_%SPIN_POS+J-1)
          ENDDO
       ENDDO
    ENDIF
  END subroutine EQUAL_DASPIN_RAY8

  subroutine EQUAL_RAY8_DASPIN(R,DS)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    TYPE(DASPIN), INTENT(IN) :: DS
    INTEGER I,J

    DO I=1,6
       R%X(I)=DS%X(I)
    ENDDO

    DO I=1,C_%ND2
       R%X(I)= MORPH(DS%M%V(I))+R%X(I)
    ENDDO

    IF(C_%track_spint_mat) THEN
       DO I=1,3
          DO J=1,3
             R%S%M(I,J)=DS%S(I,J)
          ENDDO
       ENDDO
    ELSE
       DO I=1,3
          DO J=1,3
             R%S%X(I)=DS%S(I,J)*MORPH(ONE.MONO.(C_%SPIN_POS+J-1))
          ENDDO
       ENDDO
    ENDIF


  END subroutine EQUAL_RAY8_DASPIN

  subroutine print_DASPIN(DS,MF)
    implicit none
    TYPE(DASPIN), INTENT(INOUT) :: DS
    INTEGER MF,I,J

    WRITE(MF,*) " ORBIT "
    WRITE(MF,*) DS%X(1:3)
    WRITE(MF,*) DS%X(4:6)
    WRITE(MF,*) " ORBITAL MAP "
    CALL PRINT(DS%M,MF)
    WRITE(MF,*) " SPIN MATRIX "
    DO I=1,3
       DO J=1,3
          WRITE(MF,*) " SPIN MATRIX COMPONENT ",I,J
          CALL PRINT(DS%S(I,J),MF)
       ENDDO
    ENDDO

  END subroutine print_DASPIN

  subroutine READ_DASPIN(DS,MF)
    implicit none
    TYPE(DASPIN), INTENT(INOUT) :: DS
    INTEGER MF,I,J
    CHARACTER*20 LINE
    TYPE(TAYLOR) T

    CALL ALLOC(T)

    READ(MF,*) LINE
    READ(MF,*) DS%X(1:3)
    READ(MF,*) DS%X(4:6)
    READ(MF,*) LINE
    CALL READ(DS%M,MF)
    READ(MF,*) LINE
    DO I=1,3
       DO J=1,3
          READ(MF,*) LINE
          CALL READ(T,MF)
          DS%S(I,J)=MORPH(T)
       ENDDO
    ENDDO

    CALL KILL(T)
  END subroutine READ_DASPIN



  subroutine EQUAL_NORMAL_DASPIN(R,DS)
    implicit none
    integer ipause, mypause
    TYPE(normal_spin), INTENT(INOUT) :: R
    TYPE(DASPIN), INTENT(IN) :: DS
    real(dp) s00(3,3),tunes(4)
    integer i,j
    type(real_8) n0(3),theta0,a(3,3),ai(3,3),s0(3,3),s0i(3,3),b(3)
    type(real_8) s1(3,3),s1i(3,3)
    type(damap) ri
    type(taylor) t
    call alloc_33(s1)
    call alloc_33(s1i)
    call alloc_33(s0)
    call alloc_33(s0i)
    call alloc_33(a)
    call alloc_33(ai)
    call alloc(n0,3)
    call alloc(b,3)
    call alloc(theta0)
    call alloc(ri)
    call alloc(t)

    !   normalize orbital map
    r%n=ds%m
    ri=r%n%normal
    ri=ri**(-1)
    tunes(1:3)=twopi*r%n%tune
    !
    !   do i=1,3
    !   do j=1,3
    !    r%ns%s(i,j)=morph(ds%s(i,j)%t*r%n%a_t)  ! go into normalized variables
    !    r%ns%s(i,j)=morph(ds%s(i,j)%t*ri)  ! advances by r**(-1)
    !   enddo
    !   enddo
    call trans_mat(ds%s,r%n%a_t,r%ns) ! go into normalized variables
    call trans_mat(r%ns,ri,r%ns)      ! advances by r**(-1)

    do i=1,3
       do j=1,3
          s00(i,j)=r%ns(i,j)
          s0(i,j)=s00(i,j)
       enddo
    enddo

    call find_n_theta(s0,theta0,n0)  ! s0=exp(theta0*n0.L)
    tunes(4)=theta0
    !    finds a
    call find_a(n0,a)
    call inv_as(a,ai)  ! a0**(-1) * s0 * a0 = exp(theta0*Ly)

    call matmulp(ai,r%ns,r%ns) ! a0**(-1) * r%ns%s * a0
    call matmulp(r%ns,a,r%ns) ! put result in r%ns%s

    do i=1,3
       do j=1,3
          s00(i,j)=r%ns(i,j)
          s0(i,j)=s00(i,j)
          ! call print(s0(i,j),6)
          !pause
       enddo
    enddo

    call smatp(one,a,r%as)

    call inv_as(s0,s0i)
    call matmulp(r%ns,s0i,r%ns)    !

    do i=1,c_%no+2

       call find_n_theta(r%ns,theta0,n0,my_false)  ! s0=exp(theta0*n0.L)
       !    call print(theta0,6)

       do j=1,3
          write(6,*) i,j
          t=(n0(j)%t).cut.c_%no
          call print(t,6)
          ipause=mypause(0)
       enddo

       call res_bas_spin(R,tunes,n0,b)
       theta0=one
       call exp_n_theta(theta0,b,s1)
       call inv_as(s1,s1i)
       call matmulp(r%as,s1,r%as)     ! building A
       call trans_mat(s1,ri,s1)

       call matmulp(s1i,r%ns,r%ns)
       call matmulp(r%ns,s0,r%ns)
       call matmulp(r%ns,s1,r%ns)
       call matmulp(r%ns,s0i,r%ns)

    enddo

    ri=(r%n%a_t)**(-1)
    call trans_mat(r%as,ri,r%as)

    ! checking

    call matmulp(ds%s,r%as,s0)
    call trans_mat(r%as,ds%m,s0i)
    call inv_as(s0i,s1)
    call matmulp(s1,s0,s1i)

    do i=1,3
       do j=1,3

          write(6,*) "final ",i,j
          t=s1i(i,j)%t
          !     t=s0(i,j)%t-s0i(i,j)%t
          t=t.cut.c_%no
          call print(t,6)
          ipause=mypause(0)

       enddo
    enddo

    call trans_mat(s1i,r%n%a_t,s1i)
    do i=1,3
       do j=1,3

          write(6,*) "final normal",i,j
          t=s1i(i,j)%t
          !     t=s0(i,j)%t-s0i(i,j)%t
          t=t.cut.c_%no
          call print(t,6)
          ipause=mypause(0)

       enddo
    enddo
    call kill_33(s1)
    call kill_33(s1i)
    call kill_33(s0)
    call kill_33(s0i)
    call kill_33(a)
    call kill_33(ai)
    call kill(n0,3)
    call kill(b,3)
    call kill(theta0)
    call kill(ri)
    call kill(t)


  END subroutine EQUAL_NORMAL_DASPIN

  subroutine CHECK_RES_ORBIT(J,NRES,M,SKIP)
    implicit none
    INTEGER M(:,:),NRES
    LOGICAL SKIP,SKIP1,SKIP2
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

  subroutine res_bas_spin(R,tunes,n0,b)
    implicit none
    type(real_8) n0(3),b(3)
    real(dp) tunes(4)
    type(NORMAL_SPIN) R
    type(taylorresonance) tr
    type(complextaylor) nc(2),ni(2),nr(2),ax(2),az(2)
    integer n,i,k,ndc,kk
    integer, allocatable :: j(:)
    real(dp) val
    complex(dp) r1,r2,ex
    logical(lp) reso,SKIP
    type(taylor) t

    ndc=0
    if(c_%ndpt/=0) ndc=1
    !    is=2
    allocate(j(c_%nv))
    call alloc(t)
    call alloc(tr)
    call alloc(nc,2)
    call alloc(ni,2)
    call alloc(nr,2)
    call alloc(ax,2)
    call alloc(az,2)

    tr=n0(1)%t

    nc(1)=zero
    nc(2)=zero
    call taylor_cycle(tr%cos,N)
    kk=1
    do i=1,n
       call taylor_cycle(tr%cos,kk,val,J)
       ex=one
       do k=1,c_%nd-ndc
          ex=ex*exp(i_*(j(2*k-1)-j(2*k))*tunes(k)  )
       enddo

       r1=one/(one- ex*exp(i_*tunes(4)) )
       r2=one/(one- ex*exp(-i_*tunes(4)) )
       nc(1)=nc(1)+((r1*val).mono.j)
       nc(2)=nc(2)+((r2*val).mono.j)
    enddo

    ni(1)=zero
    ni(2)=zero
    call taylor_cycle(tr%sin,N)
    kk=1
    do i=1,n
       call taylor_cycle(tr%sin,kk,val,J)

       ex=one
       do k=1,c_%nd-ndc
          ex=ex*exp(i_*(j(2*k-1)-j(2*k))*tunes(k)  )
       enddo
       r1=one/(one- ex*exp(i_*tunes(4)) )
       r2=one/(one- ex*exp(-i_*tunes(4)) )
       ni(1)=ni(1)+((r1*val).mono.j)
       ni(2)=ni(2)+((r2*val).mono.j)
    enddo

    ax(1)=nc(1)+i_*ni(1)
    ax(2)=nc(2)+i_*ni(2)

    tr=n0(3)%t

    nc(1)=zero
    nc(2)=zero
    call taylor_cycle(tr%cos,N)
    kk=1
    do i=1,n
       call taylor_cycle(tr%cos,kk,val,J)
       ex=one
       do k=1,c_%nd-ndc
          ex=ex*exp(i_*(j(2*k-1)-j(2*k))*tunes(k)  )
       enddo
       r1=one/(one- ex*exp(i_*tunes(4)) )
       r2=one/(one- ex*exp(-i_*tunes(4)) )
       nc(1)=nc(1)+((r1*val).mono.j)
       nc(2)=nc(2)+((r2*val).mono.j)
    enddo

    ni(1)=zero
    ni(2)=zero
    call taylor_cycle(tr%sin,N)
    kk=1
    do i=1,n
       call taylor_cycle(tr%sin,kk,val,J)
       ex=one
       do k=1,c_%nd-ndc
          ex=ex*exp(i_*(j(2*k-1)-j(2*k))*tunes(k)  )
       enddo
       r1=one/(one- ex*exp(i_*tunes(4)) )
       r2=one/(one- ex*exp(-i_*tunes(4)) )
       ni(1)=ni(1)+((r1*val).mono.j)
       ni(2)=ni(2)+((r2*val).mono.j)
    enddo
    az(1)=nc(1)+i_*ni(1)
    az(2)=nc(2)+i_*ni(2)


    nr(1)=zero
    nr(1)=(ax(1)+ax(2))/two +i_*(az(2)-az(1))/two
    nr(2)=zero
    nr(2)=i_*(ax(1)-ax(2))/two +(az(1)+az(2))/two

    tr%cos=nr(1)%r
    tr%sin=nr(1)%i
    t=tr
    b(1)=morph(t)
    tr%cos=nr(2)%r
    tr%sin=nr(2)%i
    t=tr
    b(3)=morph(t)

    tr=n0(2)%t

    nc(1)=zero
    call taylor_cycle(tr%cos,N)
    kk=1
    do i=1,n
       call taylor_cycle(tr%cos,kk,val,J)
       CALL CHECK_RES_ORBIT(J,R%N%NRES,R%N%M,SKIP)
       IF(.NOT.SKIP) THEN ! SKIP
          ex=one
          reso=.true.
          do k=1,c_%nd-ndc
             ex=ex*exp(i_*(j(2*k-1)-j(2*k))*tunes(k)  )
             reso=reso.and.(j(2*k-1)==j(2*k))
          enddo
          r1=(one- ex )
          if(.not.reso) then
             nc(1)=nc(1)+((val/r1).mono.j)
          endif
       ENDIF ! SKIP
    enddo

    ni(1)=zero
    call taylor_cycle(tr%sin,N)
    kk=1
    do i=1,n
       call taylor_cycle(tr%sin,kk,val,J)
       CALL  CHECK_RES_ORBIT(J,R%N%NRES,R%N%M,SKIP)
       IF(.NOT.SKIP) THEN ! SKIP
          ex=one
          reso=.true.
          do k=1,c_%nd-ndc
             ex=ex*exp(i_*(j(2*k-1)-j(2*k))*tunes(k)  )
             reso=reso.and.(j(2*k-1)==j(2*k))
          enddo
          r1=(one- ex )
          if(.not.reso) then
             ni(1)=ni(1)+((val/r1).mono.j)
          endif

       ENDIF ! SKIP
    enddo

    az(1)=nc(1)+i_*ni(1)

    tr%cos=az(1)%r
    tr%sin=az(1)%i


    t=tr
    b(2)=morph(t)


    deallocate(j)
    call kill(t)
    call kill(tr)
    call kill(nc,2)
    call kill(ni,2)
    call kill(nr,2)
    call kill(ax,2)
    call kill(az,2)


  end subroutine res_bas_spin

  subroutine exp_n_thetar(theta0,n0,s0)
    implicit none
    real(dp),intent(inout):: s0(3,3),theta0,n0(3)
    real(dp) om(3,3),xx(3,3)
    real(dp) deps,nb,n,dn,dnb
    integer i,j

    deps=1.d-7

    om=zero

    om(1,2)=-theta0*n0(3)
    om(2,1)=-om(1,2)
    om(1,3)=theta0*n0(2)
    om(3,1)=-om(1,3)
    om(2,3)=-theta0*n0(1)
    om(3,2)=-om(2,3)


    s0=zero
    xx=zero

    do i=1,3
       s0(i,i)=one
       xx(i,i)=one
    enddo

    dnb=1.e38_dp
    nb=1.e38_dp

    do i=1,1000
       xx= matmul(xx,om)/i
       s0=s0+xx
       call norm_mat(s0,n)
       !write(6,*) s0
       !pause 451
       dn=abs(n)
       if(i>10.and.(abs(dn-dnb)<deps)) then
          if(dn>=dnb) exit
       endif
       dnb=dn
    enddo

    if(i>950) then
       write(6,*) " Did not converge in exp_n_theta"
    endif


  end subroutine exp_n_thetar

  subroutine exp_n_thetap(theta0,n0,s0)
    implicit none
    type(real_8),intent(inout):: s0(3,3),theta0,n0(3)
    type(real_8) om(3,3),xx(3,3)
    real(dp) deps,nb,n,dn,dnb,sc
    integer i,j

    deps=1.d-7

    call alloc_33(om)
    call alloc_33(xx)

    om(1,2)=-theta0*n0(3)
    om(2,1)=-om(1,2)
    om(1,3)=theta0*n0(2)
    om(3,1)=-om(1,3)
    om(2,3)=-theta0*n0(1)
    om(3,2)=-om(2,3)


    !    s0=zero
    !    xx=zero
    do i=1,3
       do j=1,3
          s0(i,j)=zero
       enddo
    enddo

    do i=1,3
       s0(i,i)=one
       xx(i,i)=one
    enddo

    dnb=1.e38_dp
    nb=1.e38_dp

    do i=1,1000
       sc=one/i
       call matmulp(xx,om,xx)
       call smatp(sc,xx,xx)
       call smatmulp(one,xx,s0,s0)
       !      xx= matmul(xx,om)/i
       !      s0=s0+xx
       call norm_mat(s0,n)
       !write(6,*) s0
       !pause 451
       dn=abs(n)
       if(i>10.and.(abs(dn-dnb)<deps)) then
          if(dn>=dnb) exit
       endif
       dnb=dn
    enddo

    if(i>950) then
       write(6,*) " Did not converge in exp_n_theta"
    endif

    call kill_33(om)
    call kill_33(xx)


  end subroutine exp_n_thetap


  subroutine find_n_thetar(s0,theta0,n0)
    implicit none
    real(dp) s0(3,3),theta0,n0(3)
    real(dp) ss(3,3),sc,ln(3,3),xx(3,3),an,anb
    integer i,j
    sc=0.1_dp

    ss=s0
    call power_rot(ss,sc)

    do i=1,3
       ss(i,i)=ss(i,i)-one
    enddo

    anb=1.e38_dp
    ln=zero
    xx=ss
    do i=1,1000
       ln= xx/i+ln
       xx=-matmul(ss,xx)
       if(i>10) then
          call anti_mat(ln,an)
          if(an>=anb) exit
          anb=an
       endif
    enddo
    if(i>950) then
       write(6,*) " Did not converge in find_n_theta"
       write(6,*) i,an
    endif
    !    write(6,*) ln(1,:)
    !    write(6,*) ln(2,:)
    !    write(6,*) ln(3,:)
    n0(1)=ln(3,2)
    n0(2)=ln(1,3)
    n0(3)=ln(2,1)


    theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)
    n0=n0/theta0

    if(n0(spin_normal_position)<zero) then
       n0=-n0
       theta0=-theta0
    endif

    write(6,*) " theta0 = ",theta0
    write(6,*) " n0 = ",n0

  end subroutine find_n_thetar

  subroutine find_n_thetap(s0,theta0,n0,powerflag)
    implicit none
    type(real_8),intent(in) :: s0(3,3)
    type(real_8)  theta0,n0(3)
    type(real_8)  ss(3,3),ln(3,3),xx(3,3)
    real(dp) sc,an,anb
    integer i,j
    logical(lp), optional :: powerflag
    logical(lp):: powerf

    powerf=my_true
    call alloc_33(ss)
    call alloc_33(ln)
    call alloc_33(xx)

    sc=0.1_dp

    do i=1,3
       do j=1,3
          ln(i,j)=zero
          ss(i,j)=s0(i,j)
       enddo
    enddo

    if(present(powerflag)) powerf=powerflag
    if(powerf) call power_rot(ss,sc)

    do i=1,3
       ss(i,i)=ss(i,i)-one
    enddo

    anb=1.e38_dp
    do i=1,3
       do j=1,3
          ln(i,j)=zero
          xx(i,j)=ss(i,j)
       enddo
    enddo

    do i=1,1000

       sc=one/i
       call  smatmulp(sc,xx,ln,ln)
       !     ln= xx/i+ln
       !      xx=-matmul(ss,xx)
       call matmulp(ss,xx,xx)
       call smatp(-one,xx,xx)
       if(i>10) then
          call anti_mat(ln,an)
          if(an>=anb) exit
          anb=an
       endif
    enddo
    if(i>950) then
       write(6,*) " Did not converge in find_n_theta"
       write(6,*) i,an
    endif
    !    write(6,*) ln(1,:)
    !    write(6,*) ln(2,:)
    !    write(6,*) ln(3,:)
    n0(1)=ln(3,2)
    n0(2)=ln(1,3)
    n0(3)=ln(2,1)
    theta0=one
    if(powerf) then

       theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)
       do i=1,3
          n0(i)=n0(i)/theta0
       enddo


       if(n0(spin_normal_position)<zero) then
          do i=1,3
             n0(i)=-n0(i)
          enddo
          theta0=-theta0
       endif

    endif

    call kill_33(ss)
    call kill_33(ln)
    call kill_33(xx)
  end subroutine find_n_thetap

  subroutine find_ar(n2,a)
    implicit none
    real(dp)  n2(3),n1(3),n3(3)
    real(dp) a(3,3),s,n
    integer i,j,is

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

    write(6,*) n1
    write(6,*) n2
    write(6,*) n3
    ! spin_normal_position

    a=zero
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
    integer i,j,is
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

    !    write(6,*) n1
    !    write(6,*) n2
    !    write(6,*) n3
    ! spin_normal_position

    !a=zero
    if(spin_normal_position==2) then
       do i=1,3
          a(i,1)=n1(i)
          a(i,2)=n2(i)
          a(i,3)=n3(i)
       enddo
    elseif(spin_normal_position==3) then
       a(i,2)=n1(i)
       a(i,3)=n2(i)
       a(i,1)=n3(i)
    else
       a(i,3)=n1(i)
       a(i,1)=n2(i)
       a(i,2)=n3(i)
    endif

    call kill(n1,3)
    call kill(n3,3)
    call kill(s,n)
  end subroutine find_ap


  subroutine power_rotr(m,deps)
    implicit none
    real(dp) m(3,3),a(3,3),deps,n
    integer i
    a=m

    do i=1,1000
       a=matmul(a,m)
       call norm_mat(a,n)
       if(abs(n-3.0_dp)<deps) exit
    enddo

    !     write(6,*) i,abs(n-3.0_dp)
    if(i>999) then
       write(6,*) i,abs(n-3.0_dp)
       stop 999
    endif

  end subroutine power_rotr

  subroutine alloc_33(a)
    implicit none
    type(real_8) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call alloc(a(i,j))
       enddo
    enddo

  end subroutine alloc_33

  subroutine kill_33(a)
    implicit none
    type(real_8) a(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call kill(a(i,j))
       enddo
    enddo

  end subroutine kill_33

  subroutine power_rotp(m,deps)
    implicit none
    type(real_8) m(3,3),a(3,3)
    real(dp) deps,n
    integer i,j

    call alloc_33(a)

    do i=1,3
       do j=1,3
          a(i,j)=m(i,j)
       enddo
    enddo

    do i=1,1000
       call matmulp(a,m,a)
       call norm_mat(a,n)
       if(abs(n-3.0_dp)<deps) exit
    enddo

    !     write(6,*) i,abs(n-3.0_dp)
    if(i>999) then
       write(6,*) i,abs(n-3.0_dp)
       stop 1999
    endif


    call kill_33(a)

  end subroutine power_rotp

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

  subroutine norm_matr(m,n)
    implicit none
    integer i,j
    real(dp) m(3,3),n

    n=zero
    do i=1,3
       do j=1,3
          n=n+abs(m(i,j))
       enddo
    enddo

  end subroutine norm_matr

  subroutine norm_matp(m,n1)
    implicit none
    integer i,j
    type(real_8) m(3,3),n
    real(dp) n1
    call alloc(n)

    n=zero
    do i=1,3
       do j=1,3
          n=n+abs(m(i,j))
       enddo
    enddo

    if(n%kind==1) then
       n1=abs(n)
    else
       n1=full_abs(n%t)
    endif

    call kill(n)

  end subroutine norm_matp

  subroutine anti_matr(m,n)
    implicit none
    integer i,j
    real(dp) m(3,3),n

    n=zero
    do i=1,3
       do j=1,3
          n=n+(m(i,j))
       enddo
    enddo

  end subroutine anti_matr


  subroutine anti_matp(m,n1)
    implicit none
    integer i,j
    type(real_8) m(3,3)
    real(dp) n1
    type(real_8)  n

    call alloc(n)

    n=zero
    do i=1,3
       do j=1,3
          n=n+(m(i,j))
       enddo
    enddo

    !    if(n%kind==1) then
    n1=abs(n)
    !    else
    !     n1=full_abs(n%t)
    !    endif

    call kill(n)
  end subroutine anti_matp

  subroutine inv_asr(m,mi)
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

    call alloc_33(n)

    do i=1,3
       do j=1,3
          n(i,j)=morph(m(i,j)%t*ri)  ! advances by r**(-1)
       enddo
    enddo


    do i=1,3
       do j=1,3
          mi(i,j)=n(i,j)
       enddo
    enddo

    call kill_33(n)
  end subroutine trans_mat



  !  type daspin
  !   REAL(DP) X(6)
  !   type(damap) M
  !   type(real_8) s(3,3)
  !  end type daspin

  !  type normal_spin
  !   type(normalform) N
  !  end type normal_spin



  subroutine ALLOC_DASPIN(D)
    implicit none
    TYPE(DASPIN), INTENT(INOUT) :: D
    INTEGER I,J

    D%X=ZERO
    CALL ALLOC(D%M)
    DO I=1,3
       DO J=1,3
          CALL ALLOC(D%S(I,J))
       ENDDO
    ENDDO

  END    subroutine ALLOC_DASPIN

  subroutine KILL_DASPIN(D)
    implicit none
    TYPE(DASPIN), INTENT(INOUT) :: D
    INTEGER I,J

    D%X=ZERO
    CALL KILL(D%M)
    DO I=1,3
       DO J=1,3
          CALL KILL(D%S(I,J))
       ENDDO
    ENDDO

  END    subroutine KILL_DASPIN

  subroutine alloc_normal_spin(D)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: D

    CALL alloc(D%n)
    CALL alloc_33(D%ns)
    CALL alloc_33(D%as)
    D%NRES=0
    D%M=0

  END    subroutine alloc_normal_spin

  subroutine KILL_normal_spin(D)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: D

    CALL KILL(D%n)
    CALL KILL_33(D%ns)
    CALL KILL_33(D%as)

  END    subroutine KILL_normal_spin


  subroutine ALLOC_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    INTEGER I,J
    CALL ALLOC(S%X,3)
    DO I=1,3
       DO J=1,3
          CALL ALLOC(S%M(I,J))
       ENDDO
    ENDDO

  END    subroutine ALLOC_SPINOR_8

  subroutine ALLOC_RAY_8(R)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R

    CALL ALLOC(R%S)
    CALL ALLOC(R%X,6)

  END    subroutine ALLOC_RAY_8

  subroutine KILL_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    INTEGER I,J

    CALL KILL(S%X,3)

    DO I=1,3
       DO J=1,3
          CALL KILL(S%M(I,J))
       ENDDO
    ENDDO

  END    subroutine KILL_SPINOR_8

  subroutine KILL_RAY_8(R)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R

    CALL KILL(R%S)
    CALL KILL(R%X,6)

  END    subroutine KILL_RAY_8



  subroutine init_SPIN(NO1,ND1,NP1,NSPIN1,NDPT1,PACKAGE)
    implicit none
    integer NO1,ND1,NP1,NSPIN1
    integer,optional ::  NDPT1
    logical(lp),optional :: PACKAGE
    !  write(6,*) NO1,ND1,NP1,NSPIN1,NDPT1
    NSPIN=0
    IF(NSPIN1>1) NSPIN=NSPIN1
    IF(ND1>0) THEN
       call init(NO1,ND1,NP1+NSPIN,ndpt1,PACKAGE)
    ELSE
       call init(NO1,NP1+NSPIN,PACKAGE)
       c_%NDPT=0
       c_%ND=0
       c_%ND2=0
    ENDIF
    c_%NSPIN=NSPIN1
    C_%SPIN_POS=C_%ND2+1   ! CAN BE OVERWRITTEN BY PTC
  end subroutine  init_SPIN



end module tree_element_MODULE
