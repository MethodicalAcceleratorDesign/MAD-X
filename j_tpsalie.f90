!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size
module tpsalie
  use tpsa
  implicit none
  private  ASSVEC,ASSMAP,ASSPB,asstaylor,explieflo,expliepb,expflot,exppb
  private  checkmap,checkpb,checkvec,checktaylor,MAPmatrixr,matrixMAPr,TREEMAP
  private  TAYLORSMAP,DPEKMAP,DPOKMAP,zeroEQUALMAP,IdentityEQUALMAP
  private DABSMAP,EQUALMAP,EQUALVEC    !,EQUALMAPVEC,EQUALVECMAP
  private EQUALvecpb,EQUALpbpb,EQUALpbvec,EQUALpbda,EQUALdapb,CUTORDER,CUTORDERPB,CUTORDERVEC
  private GETORDERVEC,GETORDERMAP,GETORDERPB,concator,pushtree,pushmap,concat,pushmatrixr,push1polslow
  private trxflow,trxpb,trxtaylor !,DMULMAPsc,MULMAPsc,IMULMAPsc,DMULVECsc,MULVECsc,IMULVECsc
  !  private DMULpbsc,MULpbsc,IMULpbsc
  private scDMULMAP,scMULMAP,scIMULMAP !,scDMULVEC,scMULVEC,scIMULVEC,scDMULpb,scMULpb,scIMULpb
  private ADDMAP   !,VECMAP,MAPVEC,VECPB,PBVEC
  private SUBMAP,POWMAP,DAREADTAYLORS,DAREADMAP,DAREADVEC,DAREApb,DAREADTAYLOR
  PRIVATE DAPRINTTAYLORS,DAPRINTMAP
  private DAPRINTVEC,DAPRINTTAYLOR,DAPRINTpb,allocmap,allocvec,allocpb,alloctree,allocrad,allocrads
  private KILLmap,KILLvec,KILLpb,KILLtree,killrad,killrads
  private A_OPT_damap,k_OPT_damap,A_OPT_vecfield,K_OPT_vecfield,A_OPT_pbfield,K_OPT_pbfield
  private A_OPT_tree,K_OPT_tree
  !private push1pol,allocTAYLORS,KILLTAYLORS
  integer,private::NO,ND,ND2,NP,NDPT,NV
  logical(lp),private::old

  !frs real(dp),private::eps=c_1e_9
  integer::nrmax=400
  private mul_PBf_t,mul_VECf_t,mul_VECf_MAP,mul_PBf_MAP

  TYPE DAMAP
     TYPE (TAYLOR) V(NDIM2)    ! Ndim2=6 but allocated to nd2=2,4,6
  END TYPE DAMAP

  TYPE vecfield
     type (taylor) v(ndim2)
     integer ifac
  END TYPE vecfield

  TYPE pbfield
     type (taylor) h
     integer ifac
  END TYPE pbfield



  TYPE tree
     type (taylor) branch(ndim2)
  END TYPE tree

  !Radiation
  TYPE radtaylor
     type (taylor) v
     type (taylor) e(ndim2)
  END TYPE radtaylor


  INTERFACE ASSDAMAP
     MODULE PROCEDURE ASSVEC
     MODULE PROCEDURE ASSMAP
     MODULE PROCEDURE ASSPB
     MODULE PROCEDURE asstaylor
  END INTERFACE

  INTERFACE texp
     MODULE PROCEDURE explieflo ! flow on maps
     MODULE PROCEDURE expliepb  ! pb on maps
     MODULE PROCEDURE expflot    ! flow on taylor
     MODULE PROCEDURE exppb     ! pb on taylor
  END INTERFACE

  INTERFACE exp
     MODULE PROCEDURE explieflo ! flow on maps
     MODULE PROCEDURE expliepb  ! pb on maps
     MODULE PROCEDURE expflot    ! flow on taylor
     MODULE PROCEDURE exppb     ! pb on taylor
  END INTERFACE

  INTERFACE checkdamap
     MODULE PROCEDURE checkmap
     MODULE PROCEDURE checkpb
     MODULE PROCEDURE checkvec
     MODULE PROCEDURE checktaylor
  END INTERFACE


  INTERFACE assignment (=)
     MODULE PROCEDURE EQUALMAP
     MODULE PROCEDURE MAPTAYLORS
     MODULE PROCEDURE TAYLORSMAP
     MODULE PROCEDURE IdentityEQUALMAP
     MODULE PROCEDURE zeroEQUALMAP
     MODULE PROCEDURE MAPmatrixr
     MODULE PROCEDURE matrixMAPr
     !     MODULE PROCEDURE DABSMAP
     !     MODULE PROCEDURE ABSMAP
     MODULE PROCEDURE DPEKMAP
     MODULE PROCEDURE DPOKMAP
     MODULE PROCEDURE EQUALVEC
     MODULE PROCEDURE EQUALpbpb
     MODULE PROCEDURE EQUALpbda
     MODULE PROCEDURE EQUALdapb
     MODULE PROCEDURE EQUALvecpb
     MODULE PROCEDURE EQUALpbvec
     MODULE PROCEDURE TREEMAP
     !radiation
     MODULE PROCEDURE radEQUAL
     MODULE PROCEDURE EQUALrad
  end  INTERFACE

  INTERFACE full_abs
     MODULE PROCEDURE DABSMAP
  END INTERFACE

  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDERVEC
     MODULE PROCEDURE GETORDERMAP
     MODULE PROCEDURE GETORDERPB
  end  INTERFACE

  INTERFACE OPERATOR (.CUT.)
     MODULE PROCEDURE CUTORDER
     MODULE PROCEDURE CUTORDERPB
     MODULE PROCEDURE CUTORDERVEC
  END INTERFACE

  INTERFACE OPERATOR (.o.)
     MODULE PROCEDURE concator
     MODULE PROCEDURE trxtaylorc
     MODULE PROCEDURE trxpbc
     MODULE PROCEDURE trxflowc
  end  INTERFACE

  INTERFACE OPERATOR (*)
     !  Move just below  please : here it works
     MODULE PROCEDURE mul_PBf_MAP  !   Lines to be moved
     MODULE PROCEDURE mul_PBf_t    !   Lines to be moved
     MODULE PROCEDURE mul_VECf_t   !   Lines to be moved
     MODULE PROCEDURE mul_VECf_MAP !   Lines to be moved
     MODULE PROCEDURE pushmap
     MODULE PROCEDURE pushmatrixr
     MODULE PROCEDURE push1polslow
     MODULE PROCEDURE pushtree

     ! DA concatenation
     MODULE PROCEDURE concat
     MODULE PROCEDURE trxflow
     MODULE PROCEDURE trxpb
     MODULE PROCEDURE trxtaylor


     MODULE PROCEDURE DMULMAPsc
     MODULE PROCEDURE MULMAPsc
     MODULE PROCEDURE IMULMAPsc
     MODULE PROCEDURE scDMULMAP
     MODULE PROCEDURE scMULMAP
     MODULE PROCEDURE scIMULMAP
     !  Move just below  please
  END INTERFACE




  INTERFACE OPERATOR (+)
     MODULE PROCEDURE ADDMAP
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE SUBMAP
  END INTERFACE


  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POWMAP
  END INTERFACE


  INTERFACE DAinput
     MODULE PROCEDURE DAREADTAYLORS
     MODULE PROCEDURE DAREADMAP
     MODULE PROCEDURE DAREADVEC
     MODULE PROCEDURE DAREApb
     MODULE PROCEDURE DAREADTAYLOR
  END INTERFACE

  INTERFACE read
     MODULE PROCEDURE DAREADTAYLORS
     MODULE PROCEDURE DAREADMAP
     MODULE PROCEDURE DAREADVEC
     MODULE PROCEDURE DAREApb
     MODULE PROCEDURE DAREADTAYLOR
  END INTERFACE



  INTERFACE DAprint
     MODULE PROCEDURE DAPRINTTAYLORS
     MODULE PROCEDURE DAPRINTMAP
     MODULE PROCEDURE DAPRINTVEC
     MODULE PROCEDURE DAPRINTTAYLOR
     MODULE PROCEDURE DAPRINTpb
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE DAPRINTTAYLORS
     MODULE PROCEDURE DAPRINTMAP
     MODULE PROCEDURE DAPRINTVEC
     MODULE PROCEDURE DAPRINTTAYLOR
     MODULE PROCEDURE DAPRINTpb
  END INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE A_OPT_damap
     MODULE PROCEDURE A_OPT_vecfield
     MODULE PROCEDURE A_OPT_pbfield
     MODULE PROCEDURE A_OPT_tree
     MODULE PROCEDURE allocmap
     MODULE PROCEDURE allocvec
     MODULE PROCEDURE allocpb
     MODULE PROCEDURE alloctree
     !radiation
     MODULE PROCEDURE allocrad
     MODULE PROCEDURE allocrads
  END INTERFACE

  INTERFACE allocTPSA
     MODULE PROCEDURE allocmap
     MODULE PROCEDURE allocvec
     MODULE PROCEDURE allocpb
     MODULE PROCEDURE alloctree
     !radiation
     MODULE PROCEDURE allocrad
     MODULE PROCEDURE allocrads
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE k_OPT_damap
     MODULE PROCEDURE k_OPT_vecfield
     MODULE PROCEDURE k_OPT_pbfield
     MODULE PROCEDURE k_OPT_tree
     MODULE PROCEDURE KILLmap
     MODULE PROCEDURE KILLvec
     MODULE PROCEDURE KILLpb
     MODULE PROCEDURE KILLtree
     !radiation
     MODULE PROCEDURE KILLrad
     MODULE PROCEDURE KILLrads
  END INTERFACE
  INTERFACE KILLTPSA
     MODULE PROCEDURE KILLmap
     MODULE PROCEDURE KILLvec
     MODULE PROCEDURE KILLpb
     MODULE PROCEDURE KILLtree
     !radiation
     MODULE PROCEDURE KILLrad
     MODULE PROCEDURE KILLrads
  END INTERFACE



contains
  ! new Poisson stuff
  FUNCTION mul_PBf_t( S1, S2 )   ! Computes  s1 s2
    implicit none
    TYPE (taylor) mul_PBf_t
    TYPE (PBfield), INTENT (IN) :: S1
    TYPE (taylor)  , INTENT (IN) :: S2
    integer localmaster,i,J
    localmaster=master
    call ass(mul_PBf_t)

    mul_PBf_t=zero

    mul_PBf_t=s1%h.pb.s2

    master=localmaster
  END FUNCTION mul_PBf_t

  FUNCTION mul_VECf_t( S1, S2 )   ! Computes  s1 s2
    implicit none
    TYPE (taylor) mul_VECf_t
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (taylor)  , INTENT (IN) :: S2
    integer localmaster,i
    localmaster=master
    call ass(mul_VECf_t)

    mul_VECf_t=zero

    do i=1,c_%nd2
       mul_VECf_t=mul_VECf_t+s1%v(i)*(s2.d.i)
    enddo
    master=localmaster
  END FUNCTION mul_VECf_t


  FUNCTION mul_VECf_MAP( S1, S2 )   ! Computes  s1 s2
    implicit none
    TYPE (DAMAP) mul_VECf_MAP
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (DAMAP)  , INTENT (IN) :: S2
    integer localmaster,i
    localmaster=master

    call assDAMAP(mul_VECf_MAP)
    mul_VECf_MAP=0

    DO i=1,c_%nd2
       mul_VECf_MAP%V(I)=S1*s2%V(I)
    ENDDO


    master=localmaster
  END FUNCTION mul_VECf_MAP

  FUNCTION mul_PBf_MAP( S1, S2 )   ! Computes  s1 s2
    implicit none
    TYPE (DAMAP) mul_PBf_MAP
    TYPE (PBfield), INTENT (IN) :: S1
    TYPE (DAMAP)  , INTENT (IN) :: S2
    integer localmaster,i
    localmaster=master

    call assDAMAP(mul_PBf_MAP)
    mul_PBf_MAP=0

    DO I=1,C_%ND2
       mul_PBf_MAP%V(I)=S1*S2%V(I)
    ENDDO
    master=localmaster
  END FUNCTION mul_PBf_MAP

  ! end new Poisson stuff


  subroutine set_in_tpsalie( NO1,ND1,ND21,NP1,NDPT1,NV1,log)
    implicit none
    integer NO1,ND1,ND21,NP1,NDPT1,NV1
    logical(lp) log
    old=log
    NO=NO1
    ND=ND1
    ND2=ND21
    NP=NP1
    NDPT=NDPT1
    NV=NV1
  end  subroutine set_in_tpsalie



  SUBROUTINE  alloctree(S1)
    implicit none
    type (tree),INTENT(INOUT)::S1
    if(old) then
       call etall(s1%branch%i,nd2)
    else
       call NEWetall(s1%branch%j,nd2)
    endif
  END SUBROUTINE alloctree


  SUBROUTINE  allocmap(S1)
    implicit none
    type (damap),INTENT(INOUT)::S1

    if(old) then
       call etall(s1%v%i,nd2)
    else
       call NEWetall(s1%v%j,nd2)
    endif

  END SUBROUTINE allocmap


  SUBROUTINE  allocvec(S1)
    implicit none
    type (vecfield),INTENT(INOUT)::S1

    s1%ifac=0
    if(old) then
       call etall(s1%v%i,nd2)
    else
       call NEWetall(s1%v%j,nd2)
    endif
  END SUBROUTINE allocvec

  SUBROUTINE  allocpb(S1)
    implicit none
    type (pbfield),INTENT(INOUT)::S1
    !     call alloc(s1%h)
    s1%ifac=0
    if(old) then
       call etall1(s1%h%i)
    else
       call NEWetall(s1%h%j,1)
    endif
  END SUBROUTINE allocpb


  SUBROUTINE  KILLtree(S1)
    implicit none
    type (tree),INTENT(INOUT)::S1
    if(old) then
       call DADAL(s1%branch%i,nd2)
    else
       call newDADAL(s1%branch%j,nd2)
    endif
  END SUBROUTINE KILLtree


  SUBROUTINE  KILLmap(S1)
    implicit none
    type (damap),INTENT(INOUT)::S1
    if(old) then
       call DADAL(s1%v%i,nd2)
    else
       call newDADAL(s1%v%j,nd2)
    endif
  END SUBROUTINE KILLmap

  SUBROUTINE  KILLvec(S1)
    implicit none
    type (vecfield),INTENT(INOUT)::S1
    s1%ifac=0
    if(old) then
       call DADAL(s1%v%i,nd2)
    else
       call newDADAL(s1%v%j,nd2)
    endif
  END SUBROUTINE KILLvec

  SUBROUTINE  KILLpb(S1)
    implicit none
    type (pbfield),INTENT(INOUT)::S1
    s1%ifac=0
    if(old) then
       call DADAL1(s1%h%i)
    else
       call newDADAL(s1%h%j,1)
    endif
  END SUBROUTINE KILLpb

  SUBROUTINE  A_OPT_damap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (damap),INTENT(INout)::S1,S2
    type (damap),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
    type (damap),INTENT(INout)::S1,S2
    type (damap),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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

  SUBROUTINE  A_OPT_vecfield(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (vecfield),INTENT(INout)::S1,S2
    type (vecfield),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE A_OPT_vecfield

  SUBROUTINE  K_OPT_vecfield(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (vecfield),INTENT(INout)::S1,S2
    type (vecfield),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE K_OPT_vecfield

  SUBROUTINE  A_OPT_pbfield(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (pbfield),INTENT(INout)::S1,S2
    type (pbfield),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE A_OPT_pbfield

  SUBROUTINE  K_OPT_pbfield(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (pbfield),INTENT(INout)::S1,S2
    type (pbfield),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE K_OPT_pbfield

  SUBROUTINE  A_OPT_tree(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (tree),INTENT(INout)::S1,S2
    type (tree),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE A_OPT_tree

  SUBROUTINE  K_OPT_tree(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (tree),INTENT(INout)::S1,S2
    type (tree),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE K_OPT_tree



  SUBROUTINE  DAREADTAYLORS(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (TAYLOR),INTENT(IN)::S1(NDIM2)
    INTEGER I

    DO I=1,ND2
       CALL REA(s1(I),MFILE)
    ENDDO

  END SUBROUTINE DAREADTAYLORS

  SUBROUTINE  DAREADTAYLOR(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (TAYLOR),INTENT(IN)::S1

    CALL REA(s1,MFILE)

  END SUBROUTINE DAREADTAYLOR

  SUBROUTINE  DAREADMAP(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (damap),INTENT(in)::S1
    INTEGER I

    DO I=1,ND2
       CALL REA(s1%V(I),MFILE)
    ENDDO

  END SUBROUTINE DAREADMAP

  SUBROUTINE  DAREADVEC(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (VECFIELD),INTENT(inout)::S1
    INTEGER I

    read(mfile,*) i
    s1%ifac=i
    DO I=1,ND2
       CALL REA(s1%V(I),MFILE)
    ENDDO

  END SUBROUTINE DAREADVEC

  SUBROUTINE  DAREAPB(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (PBFIELD),INTENT(inout)::S1
    integer i
    read(mfile,*) i
    s1%ifac = i
    CALL REA(s1%H,MFILE)

  END SUBROUTINE DAREAPB

  SUBROUTINE  DAPRINTMAP(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (damap),INTENT(IN)::S1
    INTEGER I

    DO I=1,ND2
       CALL PRI(s1%V(I),MFILE)
    ENDDO
  END SUBROUTINE DAPRINTMAP

  SUBROUTINE  DAPRINTTAYLORS(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (TAYLOR),INTENT(IN)::S1(NDIM2)
    INTEGER I

    DO I=1,ND2
       CALL PRI(s1(i),MFILE)
    ENDDO
  END SUBROUTINE DAPRINTTAYLORS

  SUBROUTINE  DAPRINTVEC(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (VECFIELD),INTENT(IN)::S1
    INTEGER I

    write(mfile,*) s1%ifac,' Type of factorization represented'
    DO I=1,ND2
       CALL PRI(s1%V(I),MFILE)
    ENDDO

  END SUBROUTINE DAPRINTVEC


  SUBROUTINE  DAPRINTPB(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (PBFIELD),INTENT(IN)::S1

    write(mfile,*) s1%ifac,' Type of factorization represented'
    CALL PRI(s1%H,MFILE)

  END SUBROUTINE DAPRINTPB

  SUBROUTINE  DAPRINTTAYLOR(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (TAYLOR),INTENT(IN)::S1

    CALL PRI(s1,MFILE)

  END SUBROUTINE DAPRINTTAYLOR

  function  DABSMAP(S1)
    implicit none
    real(dp) DABSMAP
    type (damap),INTENT(IN)::S1
    real(dp) R1
    INTEGER I

    DABSMAP=zero
    R1=zero
    if(old) then
       if(s1%V(1)%i==0) call crap1("DABSMAP 1")  !call etall(s1%V%i,ND2)

       DO i=1,ND2
          !          R1=s1%V(I) !2002.10.17
          R1=full_abs(s1%V(I)) !2002.10.17
          DABSMAP=DABSMAP+R1
       ENDDO
    else
       if(.NOT. ASSOCIATED(s1%V(1)%j%r))call crap1("DABSMAP 2")  ! call newetall(s1%V%j,ND2)

       DO i=1,ND2
          !          R1=s1%V(I)
          R1=full_abs(s1%V(I))
          DABSMAP=DABSMAP+R1
       ENDDO
    endif

  END function DABSMAP


  SUBROUTINE  DPEKMAP(S2,S1)
    implicit none
    real(dp),INTENT(inOUT),dimension(:)::S2
    type (damap),INTENT(IN)::S1
    call check_snake
    if(old) then
       CALL DAPEK0(S1%V%I,S2,nd2)
    else
       CALL newDAPEK0(S1%V%J,S2,nd2)
    endif
  END SUBROUTINE DPEKMAP

  SUBROUTINE  DPOKMAP(S1,S2)
    implicit none
    real(dp),INTENT(IN),dimension(:)::S2
    type (damap),INTENT(inOUT)::S1
    if(old) then
       if(s1%V(1)%i==0) call crap1("DPOKMAP 1") !call allocw_old(s1%V(1))   !call etall(s1%V%i,ND2)
       CALL DAPOK0(S1%V%I,S2,nd2)
    else
       if(.NOT. ASSOCIATED(s1%V(1)%j%r))  call crap1("DPOKMAP 2") !  !call allocw_old(s1%V(1))  !call newetall(s1%V%j,ND2)
       CALL NEWDAPOK0(S1%V%J,S2,nd2)
    endif
  END SUBROUTINE DPOKMAP


  SUBROUTINE  TREEMAP(S1,S2)
    implicit none
    type (damap),INTENT(IN)::S2
    type (TREE),INTENT(inOUT)::S1
    integer i
    call check_snake
    do i=1,nd2
       call maketree(S2%v(i),S1%branch(i))
    enddo

  END SUBROUTINE TREEMAP



  SUBROUTINE  matrixMAPr(S2,S1)
    implicit none
    real(dp),INTENT(inOUT)::S2(:,:)            !(ndim2,ndim2)
    type (damap),INTENT(IN)::S1
    integer i,j,JL(lnv)
    real(dp) m(ndim2,ndim2)
    call check_snake

    do i=1,lnv
       JL(i)=0
    enddo

    if(old) then
       do i=1,nd2
          do j=1,nd2
             JL(j)=1
             call dapek(S1%v(i)%i,JL,m(i,j))
             JL(j)=0
          enddo
       enddo
       do i=1,nd2
          do j=1,nd2
             s2(i,j)=m(i,j)
          enddo
       enddo
    else
       do i=1,nd2
          do j=1,nd2
             jL(j)=1
             call NEWdapek(S1%v(i)%J,jL,m(i,j))
             jL(j)=0
          enddo
       enddo
       do i=1,nd2
          do j=1,nd2
             s2(i,j)=m(i,j)
          enddo
       enddo
    endif
  END SUBROUTINE matrixMAPr

  SUBROUTINE  MAPmatrixr(S1,S2)
    implicit none
    real(dp),INTENT(in)::S2(:,:)       !    (ndim2,ndim2)
    type (damap),INTENT(inout)::S1
    integer i,j,JL(lnv)

    do i=1,lnv
       JL(i)=0
    enddo

    s1=0
    if(old) then
       do i=1,nd2
          do j=1,nd2
             JL(j)=1
             call dapok(S1%v(i)%i,JL,s2(i,j))
             JL(j)=0
          enddo
       enddo
    else
       do i=1,nd2
          do j=1,nd2
             jL(j)=1
             call NEWdapok(S1%v(i)%J,jL,s2(i,j))
             jL(j)=0
          enddo
       enddo
    endif
  END SUBROUTINE MAPmatrixr


  SUBROUTINE  MAPTAYLORS(S2,S1)
    implicit none
    integer i
    type (damap),INTENT(inOUT)::S2
    type (TAYLOR),INTENT(IN),dimension(:)::S1
    if(old) then
       if(s2%V(1)%i==0) call crap1("MAPTAYLORS 1")  !call etall(s2%V%i,ND2)
       ! CALL DACOPD(S1%I,S2%v%I)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%j%r)) call crap1("MAPTAYLORS 2")  !call newetall(s2%V%j,ND2)
       ! CALL NEWDACOPD(S1%J,S2%v%J)
    endif

    do i=1,nd2
       s2%v(i)=s1(i)
    enddo
  END SUBROUTINE MAPTAYLORS


  SUBROUTINE  TAYLORSMAP(S1,S2)
    implicit none
    integer i
    type (damap),INTENT(IN)::S2
    type (TAYLOR),INTENT(inOUT),dimension(:)::S1
    call check_snake


    if(old) then
       if(s1(1)%i==0) call crap1("TAYLORSMAP 1")  !call allocw_old(s1(1))
       ! CALL DACOPD(S2%v%I,S1%I)
    else
       if(.NOT. ASSOCIATED(s1(1)%j%r)) call crap1("TAYLORSMAP 2")  !call allocw_old(s1(1))
       ! CALL NEWDACOPD(S2%v%J,S1%J)
    endif

    do i=1,nd2
       s1(i)=s2%v(i)
    enddo


  END SUBROUTINE TAYLORSMAP

  SUBROUTINE  EQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    integer i
    call check_snake

    if(old) then
       if(s2%V(1)%i==0) call crap1("EQUALMAP 1")  !call etall(s2%V%i,ND2)
       ! CALL DACOPD(S1%V%I,S2%v%I)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("EQUALMAP 2")  !call newetall(s2%V%j,ND2)
       ! CALL NEWDACOPD(S1%V%J,S2%v%J)
    endif
    do i=1,nd2
       s2%v(i)=s1%v(i)
    enddo
  END SUBROUTINE EQUALMAP



  SUBROUTINE  IdentityEQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    if(old) then
       if(s2%V(1)%i==0) call crap1("IdentityEQUALMAP 1") !call etall(s2%V%i,ND2)
       IF(S1.EQ.1) CALL etini(S2%V%I)
       IF(S1.EQ.0) CALL DACLRD(S2%V%I)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%J%r))call crap1("IdentityEQUALMAP 2") ! call newetall(s2%V%j,ND2)

       IF(S1.EQ.1) CALL NEWetini(S2%V%J)
       IF(S1.EQ.0) CALL NEWDACLRD(S2%V%J)
    endif
  END SUBROUTINE IdentityEQUALMAP


  SUBROUTINE  zeroEQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    real(dp) zero_(ndim2)
    integer i
    do i=1,ndim2
       zero_(i)=zero
    enddo
    if(old) then
       if(s2%V(1)%i==0) call crap1("zeroEQUALMAP 1") !call etall(s2%V%i,ND2)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("zeroEQUALMAP 2") !call newetall(s2%V%j,ND2)
    endif
    IF(S1.EQ.zero) s2=zero_
  END SUBROUTINE zeroEQUALMAP


  SUBROUTINE  EQUALVEC(S2,S1)
    implicit none
    type (vecfield),INTENT(inOUT)::S2
    type (vecfield),INTENT(IN)::S1
    call check_snake
    if(old) then
       if(s2%V(1)%i==0) call crap1("EQUALVEC 1") !call etall(s2%V%i,ND2)
       CALL DACOPD(S1%V%I,S2%v%I)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("EQUALVEC 2") !call newetall(s2%V%j,ND2)
       CALL NEWDACOPD(S1%V%J,S2%v%J)
    endif
    s2%ifac=s1%ifac

  END SUBROUTINE EQUALVEC

  SUBROUTINE  EQUALvecpb(S2,S1)
    implicit none
    type (vecfield),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1
    call check_snake

    if(old) then
       if(s2%V(1)%i==0) call crap1("EQUALvecpb 1")  !call etall(s2%V%i,ND2)
       CALL DIFD(S1%h%i,S2%v%i,-one)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("EQUALvecpb 2")  !call newetall(s2%V%j,ND2)
       CALL NEWDIFD(S1%h%J,S2%v%J,-one)
    endif
    s2%ifac=s1%ifac

  END SUBROUTINE EQUALvecpb

  SUBROUTINE  EQUALpbvec(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (vecfield),INTENT(IN)::S1

    call check_snake

    if(old) then
       if(s2%h%i==0) call crap1("EQUALpbvec 1")  !call etall1(s2%h%i)
       CALL intd(S1%v%i,s2%h%i,-one)
    else
       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALpbvec 2")  !call newetall(s2%h%J,1)
       CALL NEWintd(S1%v%J,s2%h%J,-one)
    endif
    s2%ifac=s1%ifac
  END SUBROUTINE EQUALpbvec


  SUBROUTINE  EQUALpbpb(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1
    call check_snake
    if(old) then
       if(s2%h%i==0) call crap1("EQUALpbpb 1")  ! call etall1(s2%h%i)
       CALL dacop(S1%h%i,S2%h%i)
    else
       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALpbpb 2")  !call newetall(s2%h%J,1)
       CALL NEWdacop(S1%h%J,S2%h%J)
    endif
    s2%ifac=s1%ifac
  END SUBROUTINE EQUALpbpb


  SUBROUTINE  EQUALpbda(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (taylor),INTENT(IN)::S1
    call check_snake

    if(old) then
       if(s2%h%i==0)  call crap1("EQUALpbda 1")  !call allocw_old(s2%h )
    else
       if(.NOT. ASSOCIATED(s2%h%J%r))  call crap1("EQUALpbda 2")  !call allocw_old(s2%h )
    endif
    s2%h=s1
  END SUBROUTINE EQUALpbda

  SUBROUTINE  EQUALdapb(S2,S1)
    implicit none
    type (taylor),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1

    call check_snake

    if(old) then
       if(s2%i==0) call crap1("EQUALdapb 1")  ! call allocw_old(s2)
       !      CALL dacop(S1%h%i,S2%i)
    else
       if(.NOT. ASSOCIATED(s2%J%r)) call crap1("EQUALdapb 2")  ! call allocw_old(s2)
       !      CALL NEWdacop(S1%h%J,S2%J)
    endif
    s2=s1%h

  END SUBROUTINE EQUALdapb


  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (damap) CUTORDER
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(CUTORDER)

    DO I=1,ND2
       CUTORDER%V(I)=(S1%V(I)).cut.S2
    ENDDO

    master=localmaster

  END FUNCTION CUTORDER

  FUNCTION CUTORDERPB( S1, S2 )
    implicit none
    TYPE (PBFIELD) CUTORDERPB
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(CUTORDERPB)

    CUTORDERPB%H=(S1%H).cut.S2
    master=localmaster

  END FUNCTION CUTORDERPB

  FUNCTION CUTORDERVEC( S1, S2 )
    implicit none
    TYPE (VECFIELD) CUTORDERVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(CUTORDERVEC)

    DO I=1,ND2
       CUTORDERVEC%V(I)=(S1%V(I)).cut.S2
    ENDDO
    master=localmaster

  END FUNCTION CUTORDERVEC

  FUNCTION GETORDERMAP( S1, S2 )
    implicit none
    TYPE (damap) GETORDERMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(GETORDERMAP)

    DO I=1,ND2
       GETORDERMAP%V(I)=(S1%V(I)).SUB.S2
    ENDDO

    master=localmaster

  END FUNCTION GETORDERMAP

  FUNCTION GETORDERVEC( S1, S2 )
    implicit none
    TYPE (VECFIELD) GETORDERVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(GETORDERVEC)

    DO I=1,ND2
       GETORDERVEC%V(I)=(S1%V(I)).SUB.S2
    ENDDO

    master=localmaster

  END FUNCTION GETORDERVEC

  FUNCTION GETORDERPB( S1, S2 )
    implicit none
    TYPE (PBFIELD) GETORDERPB
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    integer localmaster
    localmaster=master


    call checkdamap(s1)
    call assdamap(GETORDERPB)


    GETORDERPB%H=(S1%H).SUB.S2

    master=localmaster

  END FUNCTION GETORDERPB


  FUNCTION pushtree( S1, S2 )
    implicit none
    TYPE (tree), INTENT (IN) :: S1
    real(dp), intent(in),dimension(:)::s2
    real(dp) pushtree(lnv) ,junk(lnv)
    integer i

    do i=1,nd2
       junk(i)=push1pol( S1%branch(i), S2 )
    enddo


    do i=1,nd2
       pushtree(I)=junk(i)
    enddo
    do i=nd2+1,size(s2)
       pushtree(I)=S2(i)
    enddo

  END FUNCTION pushtree

  FUNCTION pushmatrixr( S1, S2 )
    implicit none
    real(dp), INTENT (IN) :: S1(ndim2,ndim2)
    real(dp), intent(in)  :: s2(ndim2)
    real(dp) pushmatrixr(ndim2) ,junk(ndim2)
    integer i,j

    do i=1,nd2
       junk(i)=zero
    enddo
    do i=1,nd2
       do j=1,nd2
          junk(i)=s1(i,j)*s2(j)+junk(i)
       enddo
    enddo


    do i=1,nd2
       pushmatrixr(I)=junk(i)
    enddo

  END FUNCTION pushmatrixr

  FUNCTION pushmap( S1, S2 )
    implicit none
    TYPE (damap), INTENT (IN) :: S1
    TYPE (tree)  arbre
    real(dp), intent(in),dimension(:)::s2
    real(dp) pushmap(lnv),junk(lnv)
    integer i

    call alloc(arbre)
    arbre=s1
    do i=1,nv
       junk(i)=s2(i)
    enddo
    pushmap=arbre*junk
    call kill(arbre)
  END FUNCTION pushmap

  FUNCTION push1pol( S1, S2 )
    implicit none
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), intent(in),dimension(:)::s2
    real(dp) push1pol

    if(old) then
       call ppush1(S1%i,s2,push1pol)
    else
       call newppush1(S1%j,s2,push1pol)
    endif


  END FUNCTION push1pol

  FUNCTION push1polslow( S1, S2 )
    implicit none
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), intent(in),dimension(:)::s2
    real(dp) push1polslow
    TYPE (taylor) t

    call alloc(t)

    call maketree(s1,t)
    if(old) then
       call ppush1(t%i,s2,push1polslow)
    else
       call newppush1(t%j,s2,push1polslow)
    endif

    call kill(t)

  END FUNCTION push1polslow






  FUNCTION concat(S1,S2)
    implicit none
    TYPE (damap) concat,t1,t2,tempnew
    TYPE (damap), INTENT (IN) :: S1, S2
    real(dp) v1(ndim2),zero_(ndim2)
    integer i
    integer localmaster
    localmaster=master
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(concat)
    call alloc(t1);call alloc(t2);call alloc(tempnew);


    do i=1,ndim2
       v1(i)=zero
       zero_(i)=zero
    enddo
    t1=s1;t2=s2;
    v1=s2
    t1=zero_;t2=zero_;
    if(old) then
       call etcct(t1%v%i,t2%v%i,tempnew%v%i)
       call dacopd(tempnew%v%i,concat%v%i)
    else
       call NEWetcct(t1%v%J,t2%v%J,tempnew%v%j)
       call NEWdacopd(tempnew%v%j,concat%v%J)
    endif
    concat=v1

    call kill(t1);call kill(t2);call kill(tempnew);
    master=localmaster

  END FUNCTION concat

  FUNCTION concator( S1, S2 )
    implicit none
    TYPE (damap) concator
    TYPE (damap), INTENT (IN) :: S1, S2
    TYPE (damap) tempnew
    integer localmaster
    localmaster=master
    call alloc(tempnew)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(concator)

    if(old) then
       call etcct(s1%v%i,s2%v%i,tempnew%v%i)
       call dacopd(tempnew%v%i,concator%v%i)
    else
       call NEWetcct(s1%v%J,s2%v%J,tempnew%v%j)
       call NEWdacopd(tempnew%v%j,concator%v%J)
    endif

    master=localmaster
    call kill(tempnew)

  END FUNCTION concator


  FUNCTION trxflow(S2,S1)
    implicit none
    TYPE (vecfield) trxflow
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    type(damap) s22,tempnew
    real(dp) zero_(ndim2)
    integer i
    integer localmaster
    localmaster=master

    call alloc(s22,tempnew)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxflow)

    do i=1,nd2
       zero_(i)=zero
    enddo
    s22=s2
    s22=zero_
    if(old) then
       call trxflo(s1%v%i,tempnew%v%i,s22%v%i)
       call dacopd(tempnew%v%i,trxflow%v%i)
    else
       call NEWtrxflo(s1%v%J,tempnew%v%J,s22%v%J)
       call NEWdacopd(tempnew%v%J,trxflow%v%J)
    endif

    call kill(s22,tempnew)
    master=localmaster
  END FUNCTION trxflow

  FUNCTION trxflowc(S1, S2 )
    implicit none
    TYPE (vecfield) trxflowc
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap) tempnew
    integer localmaster
    localmaster=master

    call alloc(tempnew)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxflowc)

    if(old) then
       call trxflo(s1%v%i,tempnew%v%i,s2%v%i)
       call dacopd(tempnew%v%i,trxflowc%v%i)
    else
       call NEWtrxflo(s1%v%J,tempnew%v%J,s2%v%J)
       call NEWdacopd(tempnew%v%J,trxflowc%v%J)
    endif

    master=localmaster
    call kill(tempnew)

  END FUNCTION trxflowc


  FUNCTION trxpb( S2, S1 )
    implicit none
    TYPE (pbfield) trxpb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap) S22
    real(dp) zero_(ndim2)
    integer i
    integer localmaster
    localmaster=master

    call alloc(s22)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxpb)

    do i=1,nd2
       zero_(i)=zero
    enddo
    s22=s2
    s22=zero_

    if(old) then
       call trx(s1%h%i,temp,s22%v%i)
       call dacop(temp,trxpb%h%i)
    else
       call NEWtrx(s1%h%J,tempL,s22%v%J)
       call NEWdacop(tempL,trxpb%h%J)
    endif

    call kill(s22)
    master=localmaster

  END FUNCTION trxpb

  FUNCTION trxpbc( S2, S1 )
    implicit none
    TYPE (pbfield) trxpbc
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxpbc)

    if(old) then
       call trx(s1%h%i,temp,s2%v%i)
       call dacop(temp,trxpbc%h%i)
    else
       call NEWtrx(s1%h%J,tempL,s2%v%J)
       call NEWdacop(tempL,trxpbc%h%J)
    endif

    master=localmaster

  END FUNCTION trxpbc

  FUNCTION trxtaylor( S1, S2 )
    implicit none
    TYPE (taylor) trxtaylor
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap)  S22
    real(dp) zero_(ndim2)
    integer i
    integer localmaster
    localmaster=master

    do i=1,nd2
       zero_(i)=zero
    enddo


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxtaylor)

    call alloc(s22)

    s22=s2
    s22=zero_
    if(old) then
       call trx(s1%i,temp,s22%v%i)
       call dacop(temp,trxtaylor%i)
    else
       call NEWtrx(s1%J,tempL,s22%v%J)
       call NEWdacop(tempL,trxtaylor%J)
    endif

    call kill(s22)
    master=localmaster

  END FUNCTION trxtaylor

  FUNCTION trxtaylorc( S1, S2 )
    implicit none
    TYPE (taylor) trxtaylorc
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxtaylorc)

    if(old) then
       call trx(s1%i,temp,s2%v%i)
       call dacop(temp,trxtaylorc%i)
    else
       call NEWtrx(s1%J,tempL,s2%v%J)
       call NEWdacop(tempL,trxtaylorc%J)
    endif

    master=localmaster

  END FUNCTION trxtaylorc



  FUNCTION texpdf( S1, S2,NRMIN,NRMAX,SCA,IFAC )
    implicit none
    TYPE (damap) texpdf
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    integer NRMIN,NRMAX,IFAC
    TYPE (damap) tempnew
    real(dp) sca
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(texpdf)
    call alloc(tempnew)

    if(old) then
       !write(6,*) NRMIN,NRMAX,SCA,IFAC
       call  FACFLOD(s1%v%i,s2%v%i,tempnew%v%i,NRMIN,NRMAX,SCA,IFAC )
       call dacopd(tempnew%v%i,texpdf%v%i)
    else
       call newFACFLOD(s1%v%j,s2%v%j,tempnew%v%j,NRMIN,NRMAX,SCA,IFAC )
       call NEWdacopd(tempnew%v%j,texpdf%v%J)
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION texpdf

  FUNCTION texpdft( S1, S2,NRMIN,NRMAX,SCA,IFAC )
    implicit none
    TYPE (taylor) texpdft
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (taylor), INTENT (IN) ::  S2
    integer NRMIN,NRMAX,IFAC
    real(dp) sca
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(texpdft)

    if(old) then
       !write(6,*) NRMIN,NRMAX,SCA,IFAC
       call  FACFLO(s1%v%i,s2%i,temp,NRMIN,NRMAX,SCA,IFAC )
       call dacop(temp,texpdft%i)
    else
       call newFACFLO(s1%v%j,s2%j,tempL,NRMIN,NRMAX,SCA,IFAC )
       call NEWdacop(tempL,texpdft%J)
    endif
    master=localmaster

  END FUNCTION texpdft


  FUNCTION explieflo( S1, S2 )
    implicit none
    TYPE (damap) explieflo
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap) tempnew

    integer no1
    integer localmaster
    localmaster=master
    no1=no

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(explieflo)
    call alloc(tempnew)
    if(s1%ifac/=0) then
       explieflo=texpdf( S1, S2,2,no1,one,s1%ifac )
       !         write(6,*) "Improper usage: map is factorized "
       !     write(6,*) "Sorry, we will implement in later version "
       !     write(6,*) "Use a Dragt-Finn or Reverse-Dragt-Finn scratch variable "
    else
       if(old) then
          call EXPFLOD(s1%v%i,s2%v%i,tempnew%v%i,eps_tpsalie,nrmax)
          call dacopd(tempnew%v%i,explieflo%v%i)
       else
          call NEWEXPFLOD(s1%v%J,s2%v%J,tempnew%v%J,eps_tpsalie,nrmax)
          call NEWdacopd(tempnew%v%J,explieflo%v%J)
       endif
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION explieflo

  FUNCTION expflot( S1, S2 )
    implicit none
    TYPE (taylor) expflot
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (taylor), INTENT (IN) ::  S2
    integer no1
    integer localmaster
    localmaster=master

    no1=no
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(expflot)

    if(s1%ifac/=0) then
       expflot=texpdft( S1, S2,2,no1,one,s1%ifac )
       !         write(6,*) "Improper usage: map is factorized "
       !    write(6,*) "Sorry, we will implement in later version "
       !     write(6,*) "Use a Dragt-Finn or Reverse-Dragt-Finn scratch variable "
    else
       if(old) then
          call EXPFLO(s1%v%i,s2%i,temp,eps_tpsalie,nrmax)
          call dacop(temp,expflot%i)
       else
          call NEWEXPFLO(s1%v%J,s2%J,tempL,eps_tpsalie,nrmax)
          call NEWdacop(tempL,expflot%J)
       endif
    endif
    master=localmaster

  END FUNCTION expflot

  FUNCTION expliepb( S1, S2 )
    implicit none
    TYPE (damap) expliepb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    type (vecfield) T
    integer localmaster
    localmaster=master

    call alloc(T)

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(expliepb)
    T=S1
    expliepb=explieflo(T,S2)

    !   if(s1%ifac/=0) then
    !    write(6,*) "Improper usage: map is factorized "
    !    write(6,*) "Sorry, we will implement in later version "
    !    write(6,*) "Use a Dragt-Finn or Reverse-Dragt-Finn scratch variable "
    !    endif
    !    if(old) then
    !     call EXPnd2(s1%h%i,s2%v%i,temp,eps_tpsalie,nrmax)
    !     call dacopd(temp,expliepb%v%i)
    !    else
    !     call NEWEXPnd2(s1%h%J,s2%v%J,tempLD,eps_tpsalie,nrmax)
    !     call NEWdacopd(tempLD,expliepb%v%J)
    !    endif
    CALL KILL(T)
    master=localmaster

  END FUNCTION expliepb

  FUNCTION exppb( S1, S2 )
    implicit none
    TYPE (taylor) exppb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (taylor), INTENT (IN) ::  S2
    type (vecfield) T
    integer localmaster
    localmaster=master

    call alloc(T)

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(exppb)
    T=S1
    exppb=expfloT(T,S2)

    !  if(s1%ifac/=0) then
    !   write(6,*) "Improper usage: map is factorized "
    !   write(6,*) "Sorry, we will implement in later version "
    !   write(6,*) "Use a Dragt-Finn or Reverse-Dragt-Finn scratch variable "
    !     endif
    !    if(old) then
    !      call EXP1D(s1%h%i,s2%i,temp,eps_tpsalie,nrmax)
    !      call dacop(temp,exppb%i)
    !     else
    !      call NEWEXP1D(s1%h%J,s2%J,tempL,eps_tpsalie,nrmax)
    !      call NEWdacop(tempL,exppb%J)
    !     endif
    CALL KILL(T)
    master=localmaster

  END FUNCTION exppb


  FUNCTION ADDMAP( S1, S2 )
    implicit none
    TYPE (damap) ADDMAP
    TYPE (damap), INTENT (IN) :: S1, S2
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(ADDMAP)

    if(old) then
       call DALIND(s1%v%i,one,s2%v%i,one,ADDMAP%v%i)
    else
       call newDALIND(s1%v%j,one,s2%v%j,one,ADDMAP%v%j)
    endif

    master=localmaster

  END FUNCTION ADDMAP

  FUNCTION SUBMAP( S1, S2 )
    implicit none
    TYPE (damap) SUBMAP
    TYPE (damap), INTENT (IN) :: S1, S2
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(SUBMAP)

    if(old) then
       call DALIND(s1%v%i,one,s2%v%i,-one,SUBMAP%v%i)
    else
       call newDALIND(s1%v%j,one,s2%v%j,-one,SUBMAP%v%j)
    endif

    master=localmaster

  END FUNCTION SUBMAP

  FUNCTION DMULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  DMULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I
    integer localmaster
    TYPE (damap)  tempnew
    localmaster=master

    call checkdamap(s1)
    call assdamap(DMULMAPsc)
    call alloc(tempnew)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,SC,tempnew%v(I)%i)
       ENDDO
       call dacopd(tempnew%v%i,DMULMAPsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,SC,tempnew%v(I)%J)
       ENDDO
       call NEWdacopd(tempnew%v%J,DMULMAPsc%v%J)
    endif

    master=localmaster
    call kill(tempnew)

  END FUNCTION DMULMAPsc

  FUNCTION MULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  MULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I
    TYPE (damap)  tempnew
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(MULMAPsc)
    call alloc(tempnew)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew%v(I)%i)
       ENDDO
       call dacopd(tempnew%v%i,MULMAPsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),tempnew%v(I)%j)
       ENDDO
       call NEWdacopd(tempnew%v%j,MULMAPsc%v%J)
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION MULMAPsc

  FUNCTION IMULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  IMULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I
    TYPE (damap)  tempnew
    integer localmaster
    localmaster=master

    call checkdamap(s1)
    call assdamap(IMULMAPsc)
    call alloc(tempnew)
    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew%v(I)%i)
       ENDDO
       call dacopd(tempnew%v%i,IMULMAPsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),tempnew%v(I)%j)
       ENDDO
       call NEWdacopd(tempnew%v%j,IMULMAPsc%v%J)
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION IMULMAPsc

  FUNCTION scDMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scDMULMAP
    TYPE (damap), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I
    TYPE (damap)  tempnew
    integer localmaster
    localmaster=master
    call checkdamap(s1)
    call assdamap(scDMULMAP)

    call alloc(tempnew)
    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,SC,tempnew%v(I)%i)
       ENDDO
       call dacopd(tempnew%v%i,scDMULMAP%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,SC,tempnew%v(i)%j)
       ENDDO
       call NEWdacopd(tempnew%v%j,scDMULMAP%v%J)
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION scDMULMAP

  FUNCTION scMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scMULMAP
    TYPE (damap), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I
    integer localmaster
    TYPE (damap)  tempnew
    localmaster=master

    call checkdamap(s1)
    call assdamap(scMULMAP)
    call alloc(tempnew)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew%v(I)%i)
       ENDDO
       call dacopd(tempnew%v%i,scMULMAP%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),tempnew%v(I)%j)
       ENDDO
       call NEWdacopd(tempnew%v%j,scMULMAP%v%J)
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION scMULMAP

  FUNCTION scIMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scIMULMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I
    integer localmaster
    TYPE (damap)  tempnew
    localmaster=master

    call checkdamap(s1)
    call assdamap(scIMULMAP)
    call alloc(tempnew)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew%v(I)%i)
       ENDDO
       call dacopd(tempnew%v%i,scIMULMAP%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),tempnew%v(I)%j)
       ENDDO
       call NEWdacopd(tempnew%v%j,scIMULMAP%v%J)
    endif
    master=localmaster
    call kill(tempnew)

  END FUNCTION scIMULMAP



  FUNCTION POWMAP( S1, R2 )
    implicit none
    TYPE (damap) POWMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (damap) S11
    INTEGER I,R22
    integer localmaster
    localmaster=master

    call checkdamap(s1)

    call assdamap(POWMAP)

    call alloc(s11)

    s11=1


    R22=IABS(R2)
    DO I=1,R22
       s11=s11*s1
    ENDDO

    IF(R2.LT.0) THEN
       if(old) then
          CALL etinv(S11%v%i,S11%v%i)
       else
          CALL newetinv(S11%v%j,S11%v%j)
       endif
    ENDIF

    powmap=s11


    ! powmap=junk
    call kill(s11)

    master=localmaster

  END FUNCTION POWMAP


  subroutine checkmap(s1)
    implicit none
    TYPE (damap) s1
    integer i
    if(old) then
       do i=1,nd2
          if(s1%v(i)%i==0) then
             w_p=0
             w_p%nc=1
             w_p=(/"Should not be here: Assign variables"/)
             w_p%fc='(1((1X,A72),/))'
             CALL WRITE_E(200)
             !             s1%v(i)%i=dummymap(i)
          endif
       enddo
    else
       do i=1,nd2
          if(.NOT. ASSOCIATED(s1%v(i)%J%r)) then
             !             s1%v(i)%j=dummymapl(i)
             w_p=0
             w_p%nc=1
             w_p=(/"Should not be here: Assign variables"/)
             w_p%fc='(1((1X,A72),/))'
             CALL WRITE_E(201)
          endif
       enddo
    endif

  end subroutine checkmap

  subroutine checkvec(s1)
    implicit none
    TYPE (vecfield) s1
    integer i

    if(old) then
       do i=1,nd2
          if(s1%v(i)%i==0) then
             !             s1%v(i)%i=dummymap(i)
             w_p=0
             w_p%nc=1
             w_p=(/"Should not be here: Assign variables"/)
             w_p%fc='(1((1X,A72),/))'
             CALL WRITE_E(202)
          endif
       enddo
    else
       do i=1,nd2
          if(.NOT. ASSOCIATED(s1%v(i)%J%r)) then
             w_p=0
             w_p%nc=1
             w_p=(/"Should not be here: Assign variables"/)
             w_p%fc='(1((1X,A72),/))'
             CALL WRITE_E(203)
             !             s1%v(i)%j=dummymapl(i)
          endif
       enddo
    endif

  end subroutine checkvec

  subroutine checkpb(s1)
    implicit none
    TYPE (pbfield) s1

    if(old) then
       if(s1%h%i==0) then
          w_p=0
          w_p%nc=1
          w_p=(/"Should not be here: Assign variables"/)
          w_p%fc='(1((1X,A72),/))'
          CALL WRITE_E(204)
          !          s1%h%i=dummy
       endif
    else
       if(.NOT. ASSOCIATED(s1%h%J%r)) then
          !          s1%h%j=dummyl
          w_p=0
          w_p%nc=1
          w_p=(/"Should not be here: Assign variables"/)
          w_p%fc='(1((1X,A72),/))'
          CALL WRITE_E(205)
       endif
    endif
  end subroutine checkpb

  subroutine checktaylor(s1)
    implicit none
    TYPE (taylor) s1
    if(old) then
       if(s1%i==0) then
          w_p=0
          w_p%nc=1
          w_p=(/"Should not be here: Assign variables"/)
          w_p%fc='(1((1X,A72),/))'
          CALL WRITE_E(206)
          !          s1%i=dummy
       endif
    else
       if(.NOT. ASSOCIATED(s1%J%r)) then
          w_p=0
          w_p%nc=1
          w_p=(/"Should not be here: Assign variables"/)
          w_p%fc='(1((1X,A72),/))'
          CALL WRITE_E(207)
          !          s1%j=dummyl
       endif
    endif
  end subroutine checktaylor



  subroutine assPB(s1)
    implicit none
    TYPE (PBFIELD) s1

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

    call ass0(s1%h)

  end subroutine assPB


  subroutine asstaylor(s1)
    implicit none
    TYPE (taylor) s1

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

    call ass0(s1)

  end subroutine asstaylor

  subroutine assvec(s1)
    implicit none
    TYPE (vecfield) s1
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

    do i=1,nd2
       call ass0(s1%v(i))
    enddo

  end subroutine assvec

  subroutine assmap(s1)
    implicit none
    TYPE (damap) s1
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

    do i=1,nd2
       call ass0(s1%v(i))
    enddo

  end subroutine assmap

  !Radiation
  SUBROUTINE  allocrad(S1)
    implicit none
    type (radtaylor),INTENT(INOUT)::S1
    if(old) then
       call etall1(s1%v%i)
       call etall(s1%e%i,nd2)
    else
       call newetall(s1%v%j,1)
       call newetall(s1%e%j,nd2)
    endif
  END SUBROUTINE allocrad

  SUBROUTINE  KILLrad(S1)
    implicit none
    type (radtaylor),INTENT(INOUT)::S1
    if(old) then
       call dadal1(s1%v%i)
       call dadal(s1%e%i,nd2)
    else
       call newdadal(s1%v%j,1)
       call newdadal(s1%e%j,nd2)
    endif

  END SUBROUTINE KILLrad

  SUBROUTINE  allocrads(S1,n)
    implicit none
    type (radtaylor),INTENT(INOUT),dimension(:)::S1
    integer i ,n
    do i=1,n
       call allocrad(s1(i))
    enddo

  END SUBROUTINE allocrads

  SUBROUTINE  KILLrads(S1,n)
    implicit none
    type (radtaylor),INTENT(INOUT),dimension(:)::s1
    integer i,n
    do i=1,n
       call killrad(s1(i))
    enddo


  END SUBROUTINE KILLrads



  SUBROUTINE  EQUALrad(S2,S1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    type (radTAYLOR),INTENT(IN)::S1
    !     if(iass0>ndum) then
    !     call ndum_warning
    !     endif
    !      iass0=0
    if(old) then
       if(s2%i==0) call crap1("EQUALrad 1")  ! call allocw_old(s2)
       CALL DACOP(S1%v%I,S2%I)
    else
       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUALrad 2")  !call allocw_old(s2)
       call newdacop(S1%v%j,S2%j)
    endif
  END SUBROUTINE EQUALrad


  SUBROUTINE  radEQUAL(S2,S1)
    implicit none
    type (radTAYLOR),INTENT(inOUT)::S2
    type (TAYLOR),INTENT(IN)::S1
    integer localmaster
    call check_snake

    !     if(iass0>ndum) then
    !     call ndum_warning
    !     endif
    !      iass0=0
    if(old) then
       if(s2%v%i==0) call crap1("radEQUAL 1")  ! call allocw_old(s2%v)
       CALL DACOP(S1%I,S2%v%I)
    else
       IF (.NOT. ASSOCIATED(s2%v%j%r)) call crap1("radEQUAL 2")  !call allocw_old(s2%v)
       call newdacop(S1%j,S2%v%j)
    endif
  END SUBROUTINE radEQUAL

end module tpsalie
