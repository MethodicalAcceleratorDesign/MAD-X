module tree_element_MODULE
  USE polymorphic_complextaylor
  IMPLICIT NONE

  PRIVATE track_TREE,track_TREEP,KILL_TREE

  type  tree_element
     real(dp) ,  DIMENSION(:), POINTER :: CC
     integer,  DIMENSION(:), POINTER :: JL,JV
     INTEGER,POINTER :: N,ND2
  end  type tree_element


  INTERFACE track
     MODULE PROCEDURE track_TREE
     MODULE PROCEDURE track_TREEP
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILL_TREE
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
    CALL ppushGETN(M%v%I,C_%ND,N)


    CALL ALLOC_TREE(T,N,C_%ND2)

    CALL ppushstore(M%v%I,C_%ND,T%CC,T%JL,T%JV)

    CALL KILL(M)

  END SUBROUTINE SET_TREE


  SUBROUTINE KILL_TREE(T)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(INOUT) :: T


    IF(ASSOCIATED(T%CC))   DEALLOCATE(T%CC,T%JL,T%JV,T%N,T%ND2)


  END SUBROUTINE KILL_TREE


  SUBROUTINE track_TREE(T,XI)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    REAL(DP), INTENT(INOUT) :: XI(6)
    REAL(DP) XT(6),XF(6),XM(7),XX
    INTEGER JC,I,IV

    XT=0.0_DP
    XF=0.0_DP
    XM=0.0_DP

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
  END SUBROUTINE track_TREE

  SUBROUTINE track_TREEP(T,XI)
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    TYPE(REAL_8), INTENT(INOUT) :: XI(6)
    TYPE(REAL_8) XT(6),XF(6),XM(7),XX
    INTEGER JC,I,IV

    CALL ALLOC(XT,6)
    CALL ALLOC(XF,6)
    CALL ALLOC(XM,7)
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

    CALL KILL(XT,6)
    CALL KILL(XF,6)
    CALL KILL(XM,7)
    CALL KILL(XX)

  END SUBROUTINE track_TREEP


  SUBROUTINE  dainput_SPECIAL6(S1,MFILE,COSY)
    implicit none
    INTEGER,INTENT(in)::MFILE,COSY
    type (damap),INTENT(INOUT)::S1
    INTEGER I,js(6),k1,k2
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
          do while(line(4:6)/='---')
             read(line,*)  k1,x,k2,js
             s1%v(i)=s1%v(i)+(x.mono.js)
             read(mfile,'(a200)') line
          enddo
       enddo

    CASE DEFAULT

       WRITE(6,*) " NOT SUPPORTED IN DAREADMAP_SPECIAL"
       STOP 111

    END SELECT


  END SUBROUTINE dainput_SPECIAL6


end module tree_element_MODULE








