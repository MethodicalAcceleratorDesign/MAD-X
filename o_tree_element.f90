!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tree_element_MODULE
  USE polymorphic_complextaylor, mon1=>mon
  IMPLICIT NONE
  public
  integer,private,parameter::ndd=6

  PRIVATE track_TREE,track_TREEP,KILL_TREE,KILL_TREE_N
  PRIVATE track_TREE_G,track_TREEP_g
  PRIVATE ALLOC_SPINOR_8,ALLOC_RAY_8
  PRIVATE KILL_SPINOR_8,KILL_RAY_8,KILL_DASPIN
  PRIVATE EQUAL_RAY8_SPINOR8,EQUAL_SPINOR8_SPINOR8,EQUAL_IDENTITY_SPINOR_8,EQUAL_SPINOR8_RAY8
  PRIVATE  EQUAL_IDENTITY_RAY_8,ALLOC_daspin,EQUAL_DASPIN_RAY8,EQUAL_RAY8_DASPIN
  private alloc_normal_spin,kill_normal_spin !,EQUAL_NORMAL_DASPIN
  private READ_DASPIN,PRINT_DASPIN,EQUAL_IDENTITY_SPINOR,EQUAL_PROBE_REAL6
  PRIVATE EQUAL_SPINOR8_SPINOR,EQUAL_PROBE8_PROBE,EQUAL_PROBE8_REAL6
  private EQUAL_IDENTITY_SPINOR_8_r3 ,EQUAL_SPINOR_SPINOR8

  private power_rotr,find_ar,norm_matr,anti_matr
  private power_rotp,find_ap,norm_matp,anti_matp
  private  find_n_thetar,find_n_thetap
  !  private smatp,smatmulp
  private exp_n_thetar,exp_n_thetap,inv_asr,inv_asp !,inv_as
  PRIVATE EQUAL_DAmapSPIN_int,daddsc,scdadd,EQUAL_PROBE8_PROBE8,PRINT_probe8
  PRIVATE concat,assmap,EQUAL_damapspin,CUTORDER,assprobe_8,POWMAP
  private read_probe8,ALLOC_33t,ALLOC_33p,KILL_33t,KILL_33p
  private purge_transverse,norm_matsr,norm_matsp
  private get_spin_nx_r,get_spin_nx_t,get_spin_nx_rd
  private scdaddo,daddsco
  private real_8REAL6,REAL6real_8,real_8REAL_8,PRINT6
  integer, target :: spin_extra_tpsa = 0

  INTERFACE assignment (=)
     !
     MODULE PROCEDURE REAL_8REAL6
     MODULE PROCEDURE REAL6REAL_8
     MODULE PROCEDURE real_8REAL_8
     !
     MODULE PROCEDURE EQUAL_IDENTITY_RAY_8
     MODULE PROCEDURE EQUAL_IDENTITY_SPINOR_8
     MODULE PROCEDURE EQUAL_IDENTITY_SPINOR
     MODULE PROCEDURE EQUAL_PROBE_REAL6
     MODULE PROCEDURE EQUAL_PROBE8_REAL6
     MODULE PROCEDURE EQUAL_PROBE8_PROBE
     MODULE PROCEDURE EQUAL_RAY8_SPINOR8
     MODULE PROCEDURE EQUAL_SPINOR8_RAY8
     MODULE PROCEDURE EQUAL_SPINOR8_SPINOR8
     MODULE PROCEDURE EQUAL_SPINOR8_SPINOR
     MODULE PROCEDURE EQUAL_DASPIN_RAY8
     MODULE PROCEDURE EQUAL_RAY8_DASPIN
     MODULE PROCEDURE EQUAL_IDENTITY_SPINOR_8_r3
     MODULE PROCEDURE EQUAL_DAmapSPIN_int
     MODULE PROCEDURE EQUAL_PROBE8_PROBE8
     MODULE PROCEDURE      EQUAL_damapspin
     MODULE PROCEDURE     EQUAL_SPINOR_SPINOR8
     MODULE PROCEDURE     EQUAL_PROBE_PROBE8
     !     MODULE PROCEDURE EQUAL_NORMAL_DASPIN
  end  INTERFACE
  INTERFACE OPERATOR (*)
     MODULE PROCEDURE concat
  END  INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POWMAP
  END  INTERFACE

  INTERFACE OPERATOR (.cut.)
     MODULE PROCEDURE CUTORDER
  END  INTERFACE

  INTERFACE operator (+)
     MODULE PROCEDURE scdadd
     MODULE PROCEDURE daddsc
     MODULE PROCEDURE scdaddo
     MODULE PROCEDURE daddsco
  END  INTERFACE


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

  INTERFACE norm_mats
     MODULE PROCEDURE norm_matsr
     MODULE PROCEDURE norm_matsp
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

  INTERFACE get_spin_nx
     MODULE PROCEDURE get_spin_nx_r
     MODULE PROCEDURE get_spin_nx_rd  ! takes DS (damapsin but return real)
     MODULE PROCEDURE get_spin_nx_t
  END INTERFACE


  INTERFACE PRINT
     MODULE PROCEDURE PRINT6
!!!!
     MODULE PROCEDURE PRINT_DASPIN
     MODULE PROCEDURE PRINT_probe8
  END INTERFACE

  INTERFACE daPRINT
     MODULE PROCEDURE PRINT6
!!!!
  END INTERFACE

  INTERFACE READ
     MODULE PROCEDURE READ_DASPIN
     MODULE PROCEDURE read_probe8   ! a bit illegal : reading polymorphs as taylor...
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

  INTERFACE ALLOC_33
     MODULE PROCEDURE ALLOC_33t
     MODULE PROCEDURE ALLOC_33p
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

  INTERFACE ass
     MODULE PROCEDURE assmap
     MODULE PROCEDURE assprobe_8
  END INTERFACE


  type normal_spin
     type(normalform) N
     type(damapspin) a1   !!! full analysis of n0
     type(damapspin) ar   !!! full analysis of n0
     type(damapspin) a_t
     type(taylor) N0(3)
     type(taylor) theta0


     real(dp) tune
     integer NRES,M(NDIM)

  end type normal_spin



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
          call ass0(s1%s(i,j))
          !       s1%s(i,j)%alloc=my_true
          !       s1%s(i,j)%kind=2
          !       s1%s(i,j)%i=0
       enddo
    enddo
    do i=1,c_%nd2
       call ass0(s1%m%v(i))
    enddo

  end subroutine assmap

  subroutine assprobe_8(s1)
    implicit none
    TYPE (probe_8) s1
    integer i

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

    if(C_%SPIN_POS/=0) then
       do i=1,3
          call ass0(s1%s%x(i)%t)
          s1%s%x(i)%alloc=my_true
          s1%s%x(i)%kind=2
          s1%s%x(i)%i=0
       enddo
    endif
    do i=1,6  !
       call ass0(s1%x(i)%t)
       s1%x(i)%alloc=my_true
       s1%x(i)%kind=1
       s1%x(i)%i=0
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
          !          if(s%s(i,j)%kind==2) then
          s%s(i,j)=s%s(i,j)*s%m
          !          endif
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
    s%s0=r%s0

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
    IF(R==1) THEN

       DO I=1,C_%NSPIN
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


  subroutine EQUAL_IDENTITY_SPINOR_8_r3(S,R)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S
    real(dp), INTENT(IN) :: R(3)
    INTEGER I

    !     S%G=A_PARTICLE

    DO I=1,C_%NSPIN
       S%X(I)=r(i)+(ONE.MONO.(c_%npara_fpp-C_%NSPIN+I))
    ENDDO

  END    subroutine EQUAL_IDENTITY_SPINOR_8_r3

  subroutine EQUAL_IDENTITY_SPINOR (S,R)
    implicit none
    TYPE(SPINOR), INTENT(INOUT) :: S
    INTEGER, INTENT(IN) :: R
    INTEGER I

    !     S%G=A_PARTICLE
    IF(R==1) THEN

       DO I=1,C_%NSPIN
          S%X(I)=ZERO
       ENDDO
    ELSEIF(R==0) THEN
       DO I=1,C_%NSPIN
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

    P%u=my_false
    P%S=0
    P%X=X
    P%S%x(2)=one

  END    subroutine EQUAL_PROBE_REAL6

  subroutine EQUAL_PROBE8_REAL6 (P,X)
    implicit none
    TYPE(PROBE_8), INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: X(6)
    INTEGER I

    P%u=my_false
    P%S=0
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
    DO I=1,3
       P8%S%X(I)=P%S%X(I)
    ENDDO

    DO I=1,6
       DO j=1,6
          P8%E_ij(I,j)=P%E_ij(I,j)
       ENDDO
    ENDDO

    P8%S=P%S
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
    P8%S=P%S
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
    P%S=P8%S
    P%u=P8%u

  END subroutine EQUAL_PROBE_PROBE8

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
    !     S%G=R%S%G

  END    subroutine EQUAL_SPINOR8_RAY8

  subroutine EQUAL_RAY8_SPINOR8(R,S)
    implicit none
    TYPE(SPINOR_8), INTENT(IN) :: S
    TYPE(probe_8), INTENT(INOUT) :: R
    R%S=S
    !     R%S%G=S%G

  END    subroutine EQUAL_RAY8_SPINOR8

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
    type(gmap) g
    real(dp) s(3,3),n(3)
    INTEGER I,J
    ! if(C_%ND2/=6) then
    !  write(6,*) " not implemented in EQUAL_DASPIN_RAY8"
    !  stop
    ! endif
    call alloc(g)
    do i=1,C_%ND2
       g%v(i)=one.mono.i
    enddo

    if(C_%ND2==4.and.(c_%npara==5.or.c_%npara==8)) then
       g%v(5)=one.mono.5
    endif

    do i=C_%SPIN_POS+3,C_%nv
       g%v(i)=one.mono.(i-3)
    enddo
    do i=C_%nv-2,c_%NV
       g%v(C_%SPIN_POS+I-C_%nv+2)=one.mono.i
    enddo

    DO I=1,C_%ND2
       DS%M%V(I)=R%X(I)
       DS%M%V(I)=DS%M%V(I)*G
    ENDDO

    DO I=1,3
       DO J=1,3
          DS%S(I,J)=((R%S%X(I)%T).D.(C_%SPIN_POS+J-1))*G
       ENDDO
    ENDDO


    DO I=1,3
       n(i)=R%S%X(I)
       DO J=1,3
          S(I,J)=DS%S(I,J)
       ENDDO
    ENDDO

    call inv_as(s,s)

    n= matmul(s,n)

    ds%s0=n


    call kill(g)

    !     DS%G=R%S%G
  END subroutine EQUAL_DASPIN_RAY8

  subroutine EQUAL_RAY8_DASPIN(R,DS)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    TYPE(damapspin), INTENT(IN) :: DS
    INTEGER I,J
    type(taylor) d
    type(gmap) g

    !  if(C_%ND2/=6) then
    !       write(6,*) " not implemented in EQUAL_DASPIN_RAY8"
    !       stop
    !      endif

    call alloc(g)
    call alloc(d)


    call alloc(g)
    do i=1,C_%ND2
       g%v(i)=one.mono.i
    enddo
    if(C_%ND2==4.and.(c_%npara==5.or.c_%npara==8)) then
       g%v(5)=one.mono.5
    endif


    do i=C_%SPIN_POS+3,C_%nv
       g%v(i)=one.mono.(i-3)
    enddo
    do i=C_%nv-2,c_%NV
       g%v(C_%SPIN_POS+I-C_%nv+2)=one.mono.i
    enddo
    g=g**(-1)

    DO I=1,C_%ND2
       R%X(I)= MORPH(DS%M%V(I))   !+R%X(I)
    ENDDO


    DO I=1,3
       R%S%X(I)=zero
       DO J=1,3
          d=DS%S(I,J)*g
          !       call print(DS%S(I,J),6)
          !       call print(d,6)
          R%S%X(I)=morph( d* ((ONE.MONO.(C_%SPIN_POS+J-1))+spin_extra_tpsa*ds%s0(j))  )+R%S%X(I)

       ENDDO
    ENDDO

    DO I=1,C_%ND2
       R%X(I)%t=R%X(I)%t*g
    ENDDO



    call kill(d)
    call kill(g)

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
    CALL ALLOC_33(DS%S)

    DO I=1,3
       DS%S(I,I)=ONE
    ENDDO


  END SUBROUTINE EQUAL_DAmapSPIN_int

  FUNCTION scdadd( S2,S1  )
    implicit none
    TYPE (probe_8) scdadd
    TYPE (damapspin), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,iia(4),ico(4),nd2,i,j
    type(taylor) d
    TYPE (probe_8) ds1
    integer      spin_extra_tpsat

    call alloc(ds1)

    spin_extra_tpsat=spin_extra_tpsa
    spin_extra_tpsa=0
    ds1=s1
    spin_extra_tpsa=spin_extra_tpsat

    call liepeek(iia,ico)
    nd2=iia(4)

    !   call ass(scdadd)
    scdadd%u=my_false
    scdadd%E_ij=zero

    call alloc(d)
    do i=1,nd2
       localmaster=master
       call ass(scdadd%x(i))
       !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
       scdadd%x(i)=ds1%x(i)+s2%x(i)-(ds1%x(i).sub.'0')
       master=localmaster
    enddo

    do i=nd2+1,6
       localmaster=master
       call ass(scdadd%x(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).AND.I==5) then   ! npr
          scdadd%x(i)=s2%x(i)+(one.mono.'00001')
       else
          scdadd%x(i)=s2%x(i)
       endif
       master=localmaster
    enddo


    if(C_%SPIN_POS/=0) then
       DO I=1,3
          localmaster=master
          call ass(scdadd%s%x(i))
          d=S2%s%x(i)
          DO J=1,3
             d=S1%S(I,J)*(ONE.MONO.(C_%SPIN_POS+J-1))+d
          ENDDO
          scdadd%s%x(i)=d
          !         scdadd%s%x(i)=ds1%s%x(i)
          master=localmaster
       ENDDO
    endif


    call kill(ds1)
    call kill(d)
  END FUNCTION scdadd

  FUNCTION daddsc( S1,S2  )
    implicit none
    TYPE (probe_8) daddsc
    TYPE (damapspin), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,iia(4),ico(4),nd2,i,j
    type(taylor) d
    TYPE (probe_8) ds1
    integer      spin_extra_tpsat

    call alloc(ds1)

    spin_extra_tpsat=spin_extra_tpsa
    spin_extra_tpsa=0
    ds1=s1
    spin_extra_tpsa=spin_extra_tpsat


    call liepeek(iia,ico)
    nd2=iia(4)
    !   call ass(daddsc)
    daddsc%u=my_false
    daddsc%E_ij=zero

    call alloc(d)
    do i=1,nd2
       localmaster=master
       call ass(daddsc%x(i))
       daddsc%x(i)=ds1%x(i)+s2%x(i)-(ds1%x(i).sub.'0')
       !       daddsc%x(i)=s1%m%v(i)+s2%x(i)
       master=localmaster
    enddo
    do i=nd2+1,6
       localmaster=master
       call ass(daddsc%x(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).AND.I==5) then   ! npr
          daddsc%x(i)=s2%x(i)+(one.mono.'00001')
       else
          daddsc%x(i)=s2%x(i)
       endif
       master=localmaster
    enddo


    if(C_%SPIN_POS/=0) then
       DO I=1,3
          localmaster=master
          call ass(daddsc%s%x(i))
          d=S2%s%x(i)
          DO J=1,3
             d=S1%S(I,J)*(ONE.MONO.(C_%SPIN_POS+J-1))+d
          ENDDO
          daddsc%s%x(i)=d
          !         daddsc%s%x(i)=ds1%s%x(i)

          master=localmaster
       ENDDO
    endif

    !      call alloc(d)
    !       d=sqrt(daddsc%s%x(1)**2+daddsc%s%x(2)**2+daddsc%s%x(3)**2)

    !       do i=1,3
    !     daddsc%s%x(i)=daddsc%s%x(i)/d
    !       enddo

    call kill(ds1)
    call kill(d)
  END FUNCTION daddsc


  subroutine print_DASPIN(DS,MF)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: DS
    INTEGER MF,I,J

    !    WRITE(MF,*) " ORBIT "
    !     WRITE(MF,*) DS%X(1:3)
    !     WRITE(MF,*) DS%X(4:6)
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

  subroutine print_probe8(DS,MF)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: DS
    INTEGER MF,I

    WRITE(MF,*) " ORBIT "
    do i=1,6
       write(mf,*) ' Variable ',i
       call print(ds%x(i),mf)
    enddo
    WRITE(MF,*) " SPIN "
    do i=1,3
       write(mf,*) ' Spin Variable ',i
       call print(ds%s%x(i),mf)
    enddo


  END subroutine print_probe8

  subroutine READ_DASPIN(DS,MF,file)
    implicit none
    TYPE(damapspin), INTENT(INOUT) :: DS
    INTEGER MF1,I,J
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
    READ(MF1,*) LINE
    CALL READ(DS%M,MF1)
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
    read(mf,*) line
    do i=1,3
       read(mf,*) line
       call read(t,mf)
       ds%s%x(i)=morph(t)
    enddo

    call kill(t)

  END subroutine read_probe8




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

  subroutine res_bas_spin(R,tunes,n0,b)
    implicit none
    type(taylor) n0(3),b(3)
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

    tr=n0(1)

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

    tr=n0(3)

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

    tr=n0(2)

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
    integer i

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
    type(taylor),intent(inout):: s0(3,3),theta0,n0(3)
    type(taylor) om(3,3),xx(3,3)
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


  subroutine find_n_thetar(s0,n0)
    implicit none
    real(dp),intent(in) :: s0(3,3)
    real(dp)  theta0,n0(3)
    real(dp)  det,ss(3,3)
    real(dp) sc,scmax
    integer i,j,is





    do i=1,3
       do j=1,3
          ss(i,j)=s0(i,j)
       enddo
    enddo



    scmax=zero
    do i=1,3
       sc=ss(i,i)
       sc=abs(sc)
       if(sc>scmax) then
          is=i
          scmax=sc
       endif
       ss(i,i)=ss(i,i)-one
    enddo

    if(is==1) then
       n0(is)=one
       det=ss(2,2)*ss(3,3)-ss(2,3)*ss(3,2)
       n0(2)=(-ss(3,3)*ss(2,1)+ss(2,3)*ss(3,1))/det
       n0(3)=(-ss(1,1)*ss(3,1)+ss(3,2)*ss(2,1))/det
    elseif(is==2) then
       n0(is)=one
       det=ss(1,1)*ss(3,3)-ss(1,3)*ss(3,1)
       n0(1)=(-ss(3,3)*ss(1,2)+ss(3,2)*ss(1,3))/det
       n0(3)=(-ss(1,1)*ss(3,2)+ss(1,2)*ss(3,1))/det
    else
       n0(is)=one
       det=ss(1,1)*ss(2,2)-ss(1,2)*ss(2,1)
       n0(1)=(-ss(2,2)*ss(1,3)+ss(3,1)*ss(1,2))/det
       n0(2)=(-ss(1,1)*ss(2,3)+ss(1,3)*ss(2,1))/det
    endif

    theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)

    do i=1,3
       n0(i)=n0(i)/theta0
    enddo



    if(n0(spin_normal_position)<zero) then
       do i=1,3
          n0(i)=-n0(i)
       enddo
    endif

  end subroutine find_n_thetar

  subroutine find_n_thetap(s0,n0)
    implicit none
    type(taylor),intent(in) :: s0(3,3)
    type(taylor)  theta0,n0(3)
    type(taylor)  det,ss(3,3)
    real(dp) sc,scmax
    integer i,j,is



    call alloc(det)
    call alloc(theta0)
    call alloc_33(ss)


    do i=1,3
       do j=1,3
          ss(i,j)=s0(i,j)
       enddo
    enddo



    scmax=zero
    do i=1,3
       sc=ss(i,i)
       sc=abs(sc)
       if(sc>scmax) then
          is=i
          scmax=sc
       endif
       ss(i,i)=ss(i,i)-one
    enddo

    if(is==1) then
       n0(is)=one
       det=ss(2,2)*ss(3,3)-ss(2,3)*ss(3,2)
       n0(2)=(-ss(3,3)*ss(2,1)+ss(2,3)*ss(3,1))/det
       n0(3)=(-ss(1,1)*ss(3,1)+ss(3,2)*ss(2,1))/det
    elseif(is==2) then
       n0(is)=one
       det=ss(1,1)*ss(3,3)-ss(1,3)*ss(3,1)
       n0(1)=(-ss(3,3)*ss(1,2)+ss(3,2)*ss(1,3))/det
       n0(3)=(-ss(1,1)*ss(3,2)+ss(1,2)*ss(3,1))/det
    else
       n0(is)=one
       det=ss(1,1)*ss(2,2)-ss(1,2)*ss(2,1)
       n0(1)=(-ss(2,2)*ss(1,3)+ss(3,1)*ss(1,2))/det
       n0(2)=(-ss(1,1)*ss(2,3)+ss(1,3)*ss(2,1))/det
    endif

    theta0=sqrt(n0(1)**2+n0(2)**2+n0(3)**2)

    do i=1,3
       n0(i)=n0(i)/theta0
    enddo



    if((n0(spin_normal_position).sub.'0')<zero) then
       do i=1,3
          n0(i)=-n0(i)
       enddo
    endif

    call kill(det)
    call kill(theta0)
    call kill_33(ss)

  end subroutine find_n_thetap

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

    !     write(6,*) n1
    !     write(6,*) n2
    !    write(6,*) n3
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
    type(taylor)  n2(3),n1(3),n3(3)
    type(taylor)  a(3,3),s,n
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


  subroutine power_rotr(m,deps,k)
    implicit none
    real(dp) m(3,3),a(3,3),deps,n
    integer i,k
    a=m

    k=1
    do i=1,100000
       a=matmul(a,m)
       k=k+1
       call norm_mats(a,n)
       if(abs(n)<deps) exit
    enddo

    !     write(6,*) i,abs(n-3.0_dp)
    if(i>100000-1) then
       write(6,*) i,abs(n-3.0_dp)
       stop 999
    endif
    m=a

  end subroutine power_rotr

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

  subroutine power_rotp(m,deps,k)
    implicit none
    type(taylor) m(3,3),a(3,3)
    real(dp) deps,n
    integer i,j,k

    call alloc_33(a)

    do i=1,3
       do j=1,3
          a(i,j)=m(i,j)
       enddo
    enddo

    k=1
    do i=1,100000
       call matmulp(a,m,a)
       k=k+1
       call norm_mats(a,n)
       if(abs(n)<deps) exit
    enddo

    !     write(6,*) i,abs(n-3.0_dp)
    if(i>100000-1) then
       write(6,*) i,abs(n-3.0_dp)
       stop 1999
    endif

    do i=1,3
       do j=1,3
          m(i,j)=a(i,j)
       enddo
    enddo

    call kill_33(a)

  end subroutine power_rotp

  subroutine matmulvp(m,n,mo)
    implicit none
    type(taylor) m(3,3),n(3),mo(3),a(3)
    integer i,j


    call alloc(a,3)

    do i=1,3
       do j=1,3
          a(i)=m(i,j)*n(j)+a(i)
       enddo
    enddo


    do i=1,3
       mo(i)=a(i)
    enddo
    call kill(a,3)

  end subroutine matmulvp

  subroutine matmulp(m,n,mo)
    implicit none
    type(taylor) m(3,3),n(3,3),mo(3,3),a(3,3)
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
    type(taylor) m(3,3),n(3,3),mo(3,3),a(3,3)
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
    type(taylor) m(3,3),mo(3,3),a(3,3)
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

  subroutine norm_matsr(m,n)
    implicit none
    integer i,j
    real(dp) m(3,3),n

    n=zero
    do i=1,3
       do j=1,3
          if(i==j) then
             n=n+abs(m(i,j)-one)
          else
             n=n+abs(m(i,j))
          endif
       enddo
    enddo

  end subroutine norm_matsr

  subroutine norm_matsp(m,n)
    implicit none
    integer i,j
    type(taylor) m(3,3)
    real(dp) n

    n=zero
    do i=1,3
       do j=1,3
          if(i==j) then
             n=n+abs(m(i,j)-one)
          else
             n=n+abs(m(i,j))
          endif
       enddo
    enddo

  end subroutine norm_matsp

  subroutine norm_matp(m,n1)
    implicit none
    integer i,j
    type(taylor) m(3,3),n
    real(dp) n1
    call alloc(n)

    n=zero
    do i=1,3
       do j=1,3
          n=n+abs(m(i,j))
       enddo
    enddo

    !    if(n%kind==1) then
    !     n1=abs(n)
    !    else
    n1=full_abs(n)
    !     n1=abs(n)
    !    endif

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
    type(taylor) m(3,3)
    real(dp) n1
    type(taylor)  n

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

  subroutine make_unitary_p(m,mi)
    !  inverse 3x3 rrotation
    implicit none
    integer i, nmax
    type(taylor) m(3,3),mi(3,3)
    type(taylor) mt(3,3),id(3,3)
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
    type(taylor) m(3,3)
    type(taylor) mt(3,3)
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
    type(taylor) m(3,3)
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
    type(taylor), INTENT(INOUT) :: m(3,3),mi(3,3)
    type(taylor) mt(3,3)

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
    type(taylor) m(3,3),mi(3,3)
    type(taylor)  n(3,3)

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
    type(taylor) m(3,3),mi(3,3)
    type(taylor)  n(3,3)
    type(damap)  ri
    type(taylor)  r

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
    TYPE(damapspin), INTENT(INOUT) :: D
    INTEGER I,J
    D%s0=ZERO
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

    D%s0=ZERO
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

  subroutine ALLOC_RAY_8(R)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    ! INTEGER I,J

    CALL ALLOC(R%S)
    CALL ALLOC(R%X,6)

    !  DO I=1,6
    !  DO J=1,6
    !   CALL ALLOC(R%E_IJ(I,J))
    !  ENDDO
    !  ENDDO

  END    subroutine ALLOC_RAY_8

  subroutine KILL_SPINOR_8(S)
    implicit none
    TYPE(SPINOR_8), INTENT(INOUT) :: S

    CALL KILL(S%X,3)


  END    subroutine KILL_SPINOR_8

  subroutine KILL_RAY_8(R)
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    ! INTEGER I,J

    CALL KILL(R%S)
    CALL KILL(R%X,6)

    !  DO I=1,6
    !  DO J=1,6
    !   CALL KILL(R%E_IJ(I,J))
    !  ENDDO
    !  ENDDO

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


!!!!!!!!!!!!!!!!   new stuff
  subroutine get_spin_nx_t(DS,theta0,n0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    type(damapspin) a1i,ds0
    type(taylor), intent(inout) :: theta0,n0(3)
    type(taylor) s(3,3),a(3,3),ai(3,3)
    type(real_8) a11,a13
    integer i,j,nsc

    call alloc(a1i)
    call alloc(ds0)
    call alloc_33(s)
    call alloc_33(a)
    call alloc_33(ai)
    call alloc(a11,a13)

    call find_n_theta(ds%s,n0)

    !     call print(theta0,6)


    call find_a(n0,a)

    call inv_as(a,ai)



    call matmulp(ai,ds%s,ai)    !
    call matmulp(ai,a,s)    !

    a11=morph(s(1,1))
    a13=morph(s(1,3))

    theta0=atan2(a13,a11)


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
    !    TYPE(damapspin), INTENT(INout) :: DS
    real(dp), intent(inout) :: theta0,n0(3)
    real(dp) a(3,3),ai(3,3)
    integer i,j,nsc
    real(dp),  INTENT(INout) :: S(3,3)

    !     do i=1,3
    !     do j=1,3
    !      s(i,j)=ds%s(i,j)
    !      enddo
    !      enddo

    call find_n_theta(s,n0)

    call find_a(n0,a)
    call inv_as(a,ai)
    ai=matmul(ai,s)    !
    s=matmul(ai,a)    !

    theta0=atan2(s(1,3),s(1,1))



  end subroutine get_spin_nx_r

  subroutine get_spin_nx_rd(dS,theta0,n0)
    implicit none
    TYPE(damapspin), INTENT(INout) :: DS
    real(dp), intent(inout) :: theta0,n0(3)
    real(dp) a(3,3),ai(3,3),S(3,3)
    integer i,j,nsc

    do i=1,3
       do j=1,3
          s(i,j)=ds%s(i,j)
       enddo
    enddo

    call find_n_theta(s,n0)

    call find_a(n0,a)
    call inv_as(a,ai)
    ai=matmul(ai,s)    !
    s=matmul(ai,a)    !

    theta0=atan2(s(1,3),s(1,1))



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

  FUNCTION concat(S2,S1)
    implicit none
    TYPE (damapspin) concat,t2
    TYPE (damapspin), INTENT (IN) :: S1, S2
    integer i,j,k
    integer localmaster
    type(taylor) s(3,3)


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
          if(c_%nd2/=0) then
             t2%s(i,j)=s2%s(i,j)*s1%m
          else
             t2%s(i,j)=s2%s(i,j)
          endif
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
    concat%s0=s1%s0

    call kill_33(s);
    call kill(t2);
    master=localmaster

  END FUNCTION concat

  real(dp) function purge_transverse(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    if(.not.c_%stable_da) return

    purge_transverse=one

    do i=1,c_%nd2
       if(j(i)/=0) then
          purge_transverse=zero
          exit
       endif
    enddo

    return
  end function purge_transverse

  subroutine clean_orbital_33(s,sf)
    implicit none
    TYPE (taylor), intent(inout) :: s(3,3),sf(3,3)
    integer i,j

    do i=1,3
       do j=1,3
          call cfu(s(i,j),purge_transverse,sf(i,j))
       enddo
    enddo

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
  subroutine NORMAL_DASPIN_new(R,DS)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: R
    TYPE(damapspin), INTENT(INout) :: DS
    !    real(dp) s00(3,3),tunes(4)
    !    integer i,j
    !   type(taylor) n0(3),a(3,3),ai(3,3),s0(3,3),s0i(3,3),b(3),theta0
    !   type(taylor) s1(3,3),s1i(3,3)
    !  type(damap) ri
    !   type(taylor) t
    !   type(taylorresonance) tr


    !   normalize orbital map
    global_verbose=.true.
    r%n%auto=.false.

    r%n=ds%m
    !    ri=r%n%normal
    !    ri=ri**(-1)
    !    tunes(1:3)=twopi*r%n%tune


  END subroutine NORMAL_DASPIN_new

  subroutine alloc_normal_spin(D)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: D

    CALL alloc(D%n0,3)
    CALL alloc(D%theta0)
    CALL alloc(D%n)
    CALL alloc(D%a_t)
    CALL alloc(D%a1)
    CALL alloc(D%ar)
    D%NRES=0
    D%M=0


  END    subroutine alloc_normal_spin

  subroutine KILL_normal_spin(D)
    implicit none
    TYPE(normal_spin), INTENT(INOUT) :: D

    CALL KILL(D%n0,3)
    CALL KILL(D%theta0)
    CALL KILL(D%n)
    CALL KILL(D%a_t)
    CALL KILL(D%a1)
    CALL KILL(D%ar)

  END    subroutine KILL_normal_spin

!!!!!!!!!!!!!!!!   new stuff
  subroutine Go_to_closed(ns,DS,a1)
    implicit none
    integer ipause, mypause
    TYPE(normal_spin), INTENT(INOUT) :: ns
    TYPE(damapspin), INTENT(INout) :: DS
    TYPE(damapspin), INTENT(INout) :: a1
    type(damapspin) a1i,ds0
    type(damapspin)s ,a ,ai
    type(taylor) nn
    !    type(taylor) s(3,3),a(3,3),ai(3,3)
    real(dp) theta0r,n0r(3),r,rt !
    type(real_8) a11,a13
    integer i,j

    call alloc(a1i)
    call alloc(ds0)
    call alloc(s)
    call alloc(a)
    call alloc(ai)
    call alloc(nn)
    !    call alloc_33(s)
    !    call alloc_33(a)
    !    call alloc_33(ai)
    call alloc(a11,a13)

    A1=1   ! damapspin to fix point


    a1%m=ns%n%a1    ! fix point map a1 of orbital normal form

    a1i=a1**(-1) ! fix point map a1 of orbital normal form

    DS0=a1i*ds*a1  ! at this stage ds0 is the map around the fixed point
    ! but spin matrix unchanged except for its dependence on transverse


    s=1
    a=1
    ai=1

    call clean_orbital_33(ds0%s,s%s)
    ! at this stage s is the spin map  without dependence on orbital

    call find_n_theta(s%s,ns%n0)

    call find_a(ns%n0,a%s)

    ai=a**(-1)
    !    call inv_as(a%s,ai%s)

    s=(ai*s)*a
    do i=1,3
       call print(ns%n0(i),6)
    enddo
    !    call find_n_theta(s%s,ns%n0)
    !do i=1,3
    ! call print(ns%n0(i),6)
    !enddo
    ipause=mypause(200)
    !    call matmulp(ai,s,ai)    !
    !    call matmulp(ai,a,s)    !

    a11=morph(s%s(1,1))
    a13=morph(s%s(1,3))

    ns%theta0=atan2(a13,a11)

    a1=a1*a

    call kill(a1i)
    call kill(ds0)
    !    call kill_33(s)
    !    call kill_33(a)
    !    call kill_33(ai)
    call kill(s)
    call kill(a)
    call kill(ai)

    call kill(a11,a13)
    call kill(nn)

  end subroutine Go_to_closed

end module tree_element_MODULE
