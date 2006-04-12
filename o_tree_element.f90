module tree_element_MODULE
  USE polymorphic_complextaylor
  IMPLICIT NONE
  public

  PRIVATE track_TREE,track_TREEP,KILL_TREE,KILL_TREE_N
  PRIVATE track_TREE_G,track_TREEP_g



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

    NULLIFY(T%CC,T%JL,T%JV,T%N,T%ND2)

  END SUBROUTINE NULL_TREE


  SUBROUTINE ALLOC_TREE(T,N,ND2)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T
    INTEGER , INTENT(IN) :: N,ND2

    !IF(N==0) RETURN

    ALLOCATE(T%CC(N),T%JL(N),T%JV(N),T%N,T%ND2)
    T%N=N
    T%ND2=ND2

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


    IF(ASSOCIATED(T%CC))   DEALLOCATE(T%CC,T%JL,T%JV,T%N,T%ND2)


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




end module tree_element_MODULE
