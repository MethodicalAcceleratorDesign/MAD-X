!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tpsalie
  use tpsa,  difd_lielib => difd, intd_lielib=> intd !
  implicit none
  public
  private  ASSVEC,ASSMAP,ASSPB,asstaylor,explieflo,expliepb,expflot,exppb
  private  checkmap,checkpb,checkvec,checktaylor,MAPmatrixr,matrixMAPr,TREEMAP
  private  TAYLORSMAP,DPEKMAP,DPOKMAP,zeroEQUALMAP,IdentityEQUALMAP
  private DABSMAP,EQUALMAP,EQUALVEC    !,EQUALMAPVEC,EQUALVECMAP
  private EQUALvecpb,EQUALpbpb,EQUALpbvec,EQUALpbda,EQUALdapb,CUTORDER,CUTORDERPB,CUTORDERVEC
  private GETORDERVEC,GETORDERMAP,GETORDERPB,concator,pushtree,concat,pushmatrixr,push1polslow
  private pushmap
  private trxflow,trxpb,trxtaylor !,DMULMAPsc,MULMAPsc,IMULMAPsc,DMULVECsc,MULVECsc,IMULVECsc
  !  private DMULpbsc,MULpbsc,IMULpbsc
  private scDMULMAP,scMULMAP,scIMULMAP !,scDMULVEC,scMULVEC,scIMULVEC,scDMULpb,scMULpb,scIMULpb
  private ADDMAP   !,VECMAP,MAPVEC,VECPB,PBVEC
  private SUBMAP,POWMAP,POWMAP_INV,DAREADTAYLORS,DAREADMAP,DAREADVEC,DAREApb,DAREADTAYLOR
  PRIVATE DAPRINTTAYLORS,DAPRINTMAP
  private DAPRINTVEC,DAPRINTTAYLOR,DAPRINTpb,allocmap,allocvec,allocpb,alloctree !,allocrad,allocrads
  private KILLmap,KILLvec,KILLpb,KILLtree  !,killrad,killrads
  private A_OPT_damap,k_OPT_damap,A_OPT_vecfield,K_OPT_vecfield,A_OPT_pbfield,K_OPT_pbfield
  private A_OPT_tree,K_OPT_tree
  !private push1pol,allocTAYLORS,KILLTAYLORS
  integer,private::NO,ND,ND2,NP,NDPT,NV
  logical(lp),private::old
  !frs real(dp),private::eps=c_1e_9
  integer::nrmax=400
  private mul_PBf_t,mul_VECf_t,mul_VECf_MAP,mul_PBf_MAP


  private A_OPT_gmap,k_OPT_gmap,allocgmap,KILLgmap,EQUALgMAP,IdentityEQUALgMAP,DAPRINTgMAP,concatorg
  private assgmap,concatg,DPEKgMAP,DPOKgMAP,gPOWMAP,trxgtaylorc,trxgtaylor,gPOWMAPtpsa,GETORDERgMAP,CUTORDERg
  private matrixtMAPr,EQUALgMAPdamap,PRINT_for_bmad_parsem,EQUALdamapgmap

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUALMAP
     MODULE PROCEDURE EQUALgMAP
     MODULE PROCEDURE EQUALgMAPdamap
     MODULE PROCEDURE EQUALdamapgmap
     MODULE PROCEDURE MAPTAYLORS
     MODULE PROCEDURE TAYLORSMAP
     MODULE PROCEDURE IdentityEQUALMAP
     MODULE PROCEDURE IdentityEQUALgMAP
     MODULE PROCEDURE zeroEQUALMAP
     MODULE PROCEDURE MAPmatrixr
     MODULE PROCEDURE matrixMAPr
     MODULE PROCEDURE matrixtMAPr   ! Taylor matrix =  damap
     !     MODULE PROCEDURE DABSMAP
     !     MODULE PROCEDURE ABSMAP
     MODULE PROCEDURE DPEKMAP
     MODULE PROCEDURE DPOKMAP
     MODULE PROCEDURE DPEKgMAP
     MODULE PROCEDURE DPOKgMAP
     MODULE PROCEDURE EQUALVEC
     MODULE PROCEDURE EQUALpbpb
     MODULE PROCEDURE EQUALpbda
     MODULE PROCEDURE EQUALdapb
     MODULE PROCEDURE EQUALvecpb
     MODULE PROCEDURE EQUALpbvec
     MODULE PROCEDURE TREEMAP
     !radiation
!     MODULE PROCEDURE radEQUAL
!     MODULE PROCEDURE EQUALrad
  end  INTERFACE

  INTERFACE OPERATOR (*)
     !  Move just below  please : here it works
     MODULE PROCEDURE mul_PBf_MAP  !   Lines to be moved
     MODULE PROCEDURE mul_PBf_t    !   Lines to be moved
     MODULE PROCEDURE mul_VECf_t   !   Lines to be moved
     MODULE PROCEDURE mul_VECf_MAP !   Lines to be moved
     MODULE PROCEDURE pushmap   ! slow lnv
     MODULE PROCEDURE pushmatrixr
     MODULE PROCEDURE push1polslow
     MODULE PROCEDURE pushtree

     ! DA concatenation
     MODULE PROCEDURE concat
     MODULE PROCEDURE concatg
     MODULE PROCEDURE trxflow
     MODULE PROCEDURE trxpb
     MODULE PROCEDURE trxtaylor
     MODULE PROCEDURE trxgtaylor


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
     MODULE PROCEDURE gPOWMAP
     MODULE PROCEDURE POWMAP_INV
  END INTERFACE

  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDERVEC
     MODULE PROCEDURE GETORDERMAP
     MODULE PROCEDURE GETORDERgMAP
     MODULE PROCEDURE GETORDERPB
  end  INTERFACE

  INTERFACE OPERATOR (.CUT.)
     MODULE PROCEDURE CUTORDER
     MODULE PROCEDURE CUTORDERg
     MODULE PROCEDURE CUTORDERPB
     MODULE PROCEDURE CUTORDERVEC
  END INTERFACE

  INTERFACE OPERATOR (.o.)
     MODULE PROCEDURE concator
     MODULE PROCEDURE concatorg
     MODULE PROCEDURE trxtaylorc
     MODULE PROCEDURE trxgtaylorc
     MODULE PROCEDURE push1polslow

     !     MODULE PROCEDURE trxpbc
     !     MODULE PROCEDURE trxflowc
  end  INTERFACE

  INTERFACE OPERATOR (.oo.)
     MODULE PROCEDURE gPOWMAPtpsa
  end  INTERFACE

  ! i/o

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
     MODULE PROCEDURE DAPRINTgMAP
     MODULE PROCEDURE DAPRINTVEC
     MODULE PROCEDURE DAPRINTTAYLOR
     MODULE PROCEDURE DAPRINTpb
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE DAPRINTTAYLORS
     MODULE PROCEDURE DAPRINTMAP
     MODULE PROCEDURE DAPRINTgMAP
     MODULE PROCEDURE DAPRINTVEC
     MODULE PROCEDURE DAPRINTTAYLOR
     MODULE PROCEDURE DAPRINTpb
  END INTERFACE


  INTERFACE print_for_bmad
     MODULE PROCEDURE print_for_bmad_parsem
  END INTERFACE

  ! Exponential of Lie Operators

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

  INTERFACE full_abs
     MODULE PROCEDURE DABSMAP
  END INTERFACE


  ! Constructors and Destructors

  INTERFACE alloc
     MODULE PROCEDURE A_OPT_damap
     MODULE PROCEDURE A_OPT_gmap
     MODULE PROCEDURE A_OPT_vecfield
     MODULE PROCEDURE A_OPT_pbfield
     MODULE PROCEDURE A_OPT_tree
     MODULE PROCEDURE allocmap
     MODULE PROCEDURE allocgmap
     MODULE PROCEDURE allocvec
     MODULE PROCEDURE allocpb
     MODULE PROCEDURE alloctree
     !radiation
!     MODULE PROCEDURE allocrad
!     MODULE PROCEDURE allocrads
  END INTERFACE

  INTERFACE allocTPSA
     MODULE PROCEDURE allocmap
     MODULE PROCEDURE allocvec
     MODULE PROCEDURE allocpb
     MODULE PROCEDURE alloctree
     !radiation
 !    MODULE PROCEDURE allocrad
 !    MODULE PROCEDURE allocrads
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE k_OPT_damap
     MODULE PROCEDURE k_OPT_gmap
     MODULE PROCEDURE k_OPT_vecfield
     MODULE PROCEDURE k_OPT_pbfield
     MODULE PROCEDURE k_OPT_tree
     MODULE PROCEDURE KILLmap
     MODULE PROCEDURE KILLgmap
     MODULE PROCEDURE KILLvec
     MODULE PROCEDURE KILLpb
     MODULE PROCEDURE KILLtree
     !radiation
   !  MODULE PROCEDURE KILLrad
   !  MODULE PROCEDURE KILLrads
  END INTERFACE

  INTERFACE KILLTPSA
     MODULE PROCEDURE KILLmap
     MODULE PROCEDURE KILLvec
     MODULE PROCEDURE KILLpb
     MODULE PROCEDURE KILLtree
     !radiation
   !  MODULE PROCEDURE KILLrad
   !  MODULE PROCEDURE KILLrads
  END INTERFACE

  ! Management routines

  INTERFACE ASSDAMAP
     MODULE PROCEDURE ASSVEC
     MODULE PROCEDURE ASSMAP
     MODULE PROCEDURE ASSgMAP
     MODULE PROCEDURE ASSPB
     MODULE PROCEDURE asstaylor
  END INTERFACE

  ! Checking routines

  INTERFACE checkdamap
     MODULE PROCEDURE checkmap
     MODULE PROCEDURE checkpb
     MODULE PROCEDURE checkvec
     MODULE PROCEDURE checktaylor
  END INTERFACE


contains
  ! new Poisson stuff
  FUNCTION mul_PBf_t( S1, S2 )   ! Computes  s1 s2
    implicit none
    TYPE (taylor) mul_PBf_t
    TYPE (PBfield), INTENT (IN) :: S1
    TYPE (taylor)  , INTENT (IN) :: S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       mul_PBf_t%i=0
      RETURN
    endif
    localmaster=master
    call ass(mul_PBf_t)

    mul_PBf_t=0.0_dp

    mul_PBf_t=s1%h.pb.s2

    master=localmaster
  END FUNCTION mul_PBf_t

  FUNCTION mul_VECf_t( S1, S2 )   ! Computes  s1 s2
    implicit none
    TYPE (taylor) mul_VECf_t
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (taylor)  , INTENT (IN) :: S2
    integer localmaster,i
    IF(.NOT.C_%STABLE_DA) then
       mul_VECf_t%i=0
      RETURN
    endif
    localmaster=master

    call ass(mul_VECf_t)

    mul_VECf_t=0.0_dp

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
    IF(.NOT.C_%STABLE_DA) then
       mul_VECf_MAP%v%i=0
      RETURN
    endif
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
    IF(.NOT.C_%STABLE_DA) then
       mul_PBf_MAP%v%i=0
      RETURN
    endif
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
     integer i
   do i=1,nd2
    call alloc(s1%branch(i))
   enddo

  END SUBROUTINE alloctree

 ! SUBROUTINE  damap_clean(S1,value)
 !   implicit none
 !   type (damap),INTENT(INOUT)::S1
 !   real(dp),INTENT(INOUT)::value
 !   INTEGER I

!    DO I=1,ND2
!       CALL taylor_clean(S1%V(I),value)
!    ENDDO

!  END SUBROUTINE damap_clean

  SUBROUTINE  allocmap(S1)
    implicit none
    type (damap),INTENT(INOUT)::S1
    INTEGER I

    DO I=1,ND2
       CALL ALLOC(S1%V(I))
    ENDDO
    ! if(old) then
    ! call etall(s1%v%i,nd2)
    !    else
    !       call NEWetall(s1%v%j,nd2)
    !    endif

  END SUBROUTINE allocmap

  SUBROUTINE  allocgmap(S1,n)
    implicit none
    type (gmap),INTENT(INOUT)::S1
    integer, optional :: n
    integer m,i

    m=nv
    if(present(n)) m=n
    s1%n=m

    ! if(old) then
    do i=1,s1%n
     call alloc(s1%v(i))
    enddo
    !    else
    !       call NEWetall(s1%v%j,s1%n)
    !    endif

  END SUBROUTINE allocgmap


  SUBROUTINE  allocvec(S1)
    implicit none
    type (vecfield),INTENT(INOUT)::S1
    integer i
    s1%ifac=0
    ! if(old) then
    do i=1,nd2
     call alloc(s1%v(i))
    enddo
    !    else
    !       call NEWetall(s1%v%j,nd2)
    !    endif
  END SUBROUTINE allocvec

  SUBROUTINE  allocpb(S1)
    implicit none
    type (pbfield),INTENT(INOUT)::S1
    !     call alloc(s1%h)
    s1%ifac=0
    ! if(old) then
    call alloc(s1%h)
    !    else
    !       call NEWetall(s1%h%j,1)
    !    endif
  END SUBROUTINE allocpb


  SUBROUTINE  KILLtree(S1)
    implicit none
    type (tree),INTENT(INOUT)::S1
    integer i
    ! if(old) then
     do i=1,nd2
    call kill(s1%branch(i))
    enddo
    !    else
    !       call newDADAL(s1%branch%j,nd2)
    !    endif
  END SUBROUTINE KILLtree


  SUBROUTINE  KILLmap(S1)
    implicit none
    type (damap),INTENT(INOUT)::S1
    INTEGER I
    ! if(old) then
    DO I=1,ND2
       CALL KILL(s1%v(I))
    ENDDO
    !    call DADAL(s1%v%i,nd2)
    !    else
    !       call newDADAL(s1%v%j,nd2)
    !    endif
  END SUBROUTINE KILLmap

  SUBROUTINE  KILLgmap(S1)
    implicit none
    type (gmap),INTENT(INOUT)::S1
    integer i
    ! if(old) then
    do i=1,s1%n
     call kill(s1%v(i))
    enddo
    !    else
    !       call newDADAL(s1%v%j,s1%n)
    !    endif
  END SUBROUTINE KILLgmap

  SUBROUTINE  KILLvec(S1)
    implicit none
    type (vecfield),INTENT(INOUT)::S1
     integer i
    s1%ifac=0
    ! if(old) then
    do i=1,nd2
     call kill(s1%v(i))
    enddo

    !    else
    !       call newDADAL(s1%v%j,nd2)
    !    endif
  END SUBROUTINE KILLvec

  SUBROUTINE  KILLpb(S1)
    implicit none
    type (pbfield),INTENT(INOUT)::S1
    s1%ifac=0
    ! if(old) then
    call kill(s1%h)
    !    else
    !       call newDADAL(s1%h%j,1)
    !    endif
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

  SUBROUTINE  A_OPT_gmap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10,n)
    implicit none
    type (gmap),INTENT(INout)::S1,S2
    type (gmap),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    integer, optional :: n
    integer m

    m=nv
    if(present(n)) m=n

    call alloc(s1,n=m)
    call alloc(s2,n=m)
    if(present(s3)) call alloc(s3,n=m)
    if(present(s4)) call alloc(s4,n=m)
    if(present(s5)) call alloc(s5,n=m)
    if(present(s6)) call alloc(s6,n=m)
    if(present(s7)) call alloc(s7,n=m)
    if(present(s8)) call alloc(s8,n=m)
    if(present(s9)) call alloc(s9,n=m)
    if(present(s10))call alloc(s10,n=m)
  END SUBROUTINE A_OPT_gmap

  SUBROUTINE  K_OPT_gmap(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (gmap),INTENT(INout)::S1,S2
    type (gmap),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE K_OPT_gmap

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
    type (TAYLOR),INTENT(INOUT)::S1(NDIM2)
    INTEGER I

    DO I=1,ND2
       CALL REA(s1(I),MFILE)
    ENDDO

  END SUBROUTINE DAREADTAYLORS

  SUBROUTINE  DAREADTAYLOR(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (TAYLOR),INTENT(INOUT)::S1

    CALL REA(s1,MFILE)

  END SUBROUTINE DAREADTAYLOR

  SUBROUTINE  DAREADMAP(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (damap),INTENT(INOUT)::S1
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

  SUBROUTINE  DAPRINTMAP(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (damap),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    INTEGER I

    DO I=1,ND2
       CALL PRI(s1%V(I),MFILE,PREC,ind=i)
    ENDDO
  END SUBROUTINE DAPRINTMAP

  SUBROUTINE  PRINT_for_bmad_parsem(S1,MFILE,ref0,ref1,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (damap),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    REAL(DP),OPTIONAL,INTENT(IN)::ref0(6),ref1(6)
    INTEGER I,mfi
    type(damap) t,idf
    character(255) line
     mfi=6
     if(present(mfile)) mfi=mfile

    call alloc(t,idf)
     t=s1
     if(present(ref0)) then
 !  adding orbit for bmad
       do i=1,nd2
        idf%v(i)=(1.d0.mono.i)-ref0(i)
       enddo
       do i=1,nd2
        t%v(i)=t%v(i)-(t%v(i).sub.'0')
       enddo
       t=t.o.idf
       do i=1,nd2
        t%v(i)=t%v(i)+ref1(i)
       enddo
       write(line,*) "ref_orbit=(",ref0(1),",",ref0(2),",",ref0(3),",",ref0(4),",",ref0(5),",",ref0(6),"),"
     endif
       call context(line)
       line=adjustl(line)
       write(mfi,'(a255)') line
    DO I=1,ND2
       CALL PRINT_for_bmad_parser(t%V(I),MFILE,PREC,ind=i)
    ENDDO
    call KILL(t,idf)
  END SUBROUTINE PRINT_for_bmad_parsem

  SUBROUTINE  DAPRINTgMAP(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (gmap),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    INTEGER I

    DO I=1,s1%n
       CALL PRI(s1%V(I),MFILE,PREC)
    ENDDO
  END SUBROUTINE DAPRINTgMAP

  SUBROUTINE  DAPRINTTAYLORS(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (TAYLOR),INTENT(IN)::S1(:)
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    INTEGER I,mfi
     mfi=6
     if(present(mfile)) mfi=mfile

    DO I=1,size(S1)
       if(s1(i)%i>0) then
          if(size(S1)>1) write(mfi,*) "Taylor #",i
          CALL PRI(s1(i),MFILE,PREC)
       endif
    ENDDO
  END SUBROUTINE DAPRINTTAYLORS

  SUBROUTINE  DAPRINTVEC(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (VECFIELD),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    INTEGER I,mfi
     mfi=6
     if(present(mfile)) mfi=mfile

    write(mfi,*) s1%ifac,' Factorization represented'
    DO I=1,ND2
       CALL PRI(s1%V(I),MFILE,PREC)
    ENDDO

  END SUBROUTINE DAPRINTVEC


  SUBROUTINE  DAPRINTPB(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (PBFIELD),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    integer mfi
     mfi=6
     if(present(mfile)) mfi=mfile

    write(mfi,*) s1%ifac,' Factorization represented'
    CALL PRI(s1%H,MFILE,PREC)

  END SUBROUTINE DAPRINTPB

  SUBROUTINE  DAPRINTTAYLOR(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (TAYLOR),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC

    CALL PRI(s1,MFILE,PREC)

  END SUBROUTINE DAPRINTTAYLOR

  function  DABSMAP(S1)
    implicit none
    real(dp) DABSMAP
    type (damap),INTENT(IN)::S1
    real(dp) R1
    INTEGER I
    IF(.NOT.C_%STABLE_DA) then
       R1=0
      RETURN
    endif

    DABSMAP=0.0_dp
    R1=0.0_dp
    ! if(old) then
  !  if(s1%V(1)%i==0) call crap1("DABSMAP 1")  !call etall(s1%V%i,ND2)

    DO i=1,ND2
       !          R1=s1%V(I) !2002.10.17
       R1=full_abs(s1%V(I)) !2002.10.17
       DABSMAP=DABSMAP+R1
    ENDDO
    !    else
    !       if(.NOT. ASSOCIATED(s1%V(1)%j%r))call crap1("DABSMAP 2")  ! call newetall(s1%V%j,ND2)
    !
    !       DO i=1,ND2
    !          !          R1=s1%V(I)
    !          R1=full_abs(s1%V(I))
    !          DABSMAP=DABSMAP+R1
    !       ENDDO
    !    endif

  END function DABSMAP


  SUBROUTINE  DPEKMAP(S2,S1)
    implicit none
    real(dp),INTENT(inOUT),dimension(:)::S2
    type (damap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then


     call PEK0(S1%V,S2,nd2)

    !    else
    !       CALL newDAPEK0(S1%V%J,S2,nd2)
    !    endif
  END SUBROUTINE DPEKMAP

  SUBROUTINE  DPEKgMAP(S2,S1)
    implicit none
    real(dp),INTENT(inOUT),dimension(:)::S2
    type (gmap),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
     call PEK0(S1%V,S2,s1%n)

     !    else
    !       CALL newDAPEK0(S1%V%J,S2,s1%n)
    !    endif
  END SUBROUTINE DPEKgMAP

  SUBROUTINE  DPOKMAP(S1,S2)
    implicit none
    real(dp),INTENT(IN),dimension(:)::S2
    type (damap),INTENT(inOUT)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
 !   if(s1%V(1)%i==0) call crap1("DPOKMAP 1") !call allocw_old(s1%V(1))   !call etall(s1%V%i,ND2)
!
    call pok0(s1%v,S2,nd2)
  !  CALL DAPOK0(S1%V%I,S2,nd2)
    !    else
    !       if(.NOT. ASSOCIATED(s1%V(1)%j%r))  call crap1("DPOKMAP 2") !  !call allocw_old(s1%V(1))  !call newetall(s1%V%j,ND2)
    !       CALL NEWDAPOK0(S1%V%J,S2,nd2)
    !    endif
  END SUBROUTINE DPOKMAP

  SUBROUTINE  DPOKgMAP(S1,S2)
    implicit none
    real(dp),INTENT(IN),dimension(:)::S2
    type (gmap),INTENT(inOUT)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
 !   if(s1%V(1)%i==0) call crap1("DPOKMAP 1") !call allocw_old(s1%V(1))   !call etall(s1%V%i,ND2)
    CALL POK0(S1%V,S2,s1%n)
    !    else
    !       if(.NOT. ASSOCIATED(s1%V(1)%j%r))  call crap1("DPOKMAP 2") !  !call allocw_old(s1%V(1))  !call newetall(s1%V%j,ND2)
    !       CALL NEWDAPOK0(S1%V%J,S2,s1%n)
    !    endif
  END SUBROUTINE DPOKgMAP

  SUBROUTINE  TREEMAP(S1,S2)
    implicit none
    type (damap),INTENT(IN)::S2
    type (TREE),INTENT(inOUT)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    do i=1,nd2
       call maketree(S2%v(i),S1%branch(i))
    enddo

  END SUBROUTINE TREEMAP

  SUBROUTINE  matrixtMAPr(S2,S1)
    implicit none
    type(taylor),INTENT(inOUT)::S2(:,:)            !(ndim2,ndim2)
    type (damap),INTENT(IN)::S1
    integer i,j
    type(taylor) m(ndim2,ndim2)
    integer, allocatable :: jl(:)

    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake


    allocate(JL(nd2))
    do i=1,nd2
       do j=1,nd2
          call alloc(m(i,j))
       enddo
    enddo
    ! if(old) then
    do i=1,nd2
       do j=1,nd2
          JL(j)=1
          m(i,j)=S1%v(i).par.jl
          JL(j)=0
       enddo
    enddo
    do i=1,nd2
       do j=1,nd2
          s2(i,j)=m(i,j)
       enddo
    enddo

    do i=1,nd2
       do j=1,nd2
          call kill(m(i,j))
       enddo
    enddo

    deallocate(jl)

  END SUBROUTINE matrixtMAPr



  SUBROUTINE  matrixMAPr(S2,S1)
    implicit none
    real(dp),INTENT(inOUT)::S2(:,:)            !(ndim2,ndim2)
    type (damap),INTENT(IN)::S1
    integer i,j,JL(lnv)
    real(dp) m(ndim2,ndim2)
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    do i=1,lnv
       JL(i)=0
    enddo

    ! if(old) then
    do i=1,nd2
       do j=1,nd2
          JL(j)=1
          call pek(S1%v(i),JL,m(i,j))
          JL(j)=0
       enddo
    enddo
    do i=1,nd2
       do j=1,nd2
          s2(i,j)=m(i,j)
       enddo
    enddo

  END SUBROUTINE matrixMAPr

  SUBROUTINE  MAPmatrixr(S1,S2)
    implicit none
    real(dp),INTENT(in)::S2(:,:)       !    (ndim2,ndim2)
    type (damap),INTENT(inout)::S1
    integer i,j,JL(lnv)
    IF(.NOT.C_%STABLE_DA) RETURN

    do i=1,lnv
       JL(i)=0
    enddo

    s1=0
    ! if(old) then
    do i=1,nd2
       do j=1,nd2
          JL(j)=1
          call pok(S1%v(i),JL,s2(i,j))
          JL(j)=0
       enddo
    enddo


  END SUBROUTINE MAPmatrixr


  SUBROUTINE  MAPTAYLORS(S2,S1)
    implicit none
    integer i
    type (damap),INTENT(inOUT)::S2
    type (TAYLOR),INTENT(IN),dimension(:)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
   ! if(s2%V(1)%i==0) call crap1("MAPTAYLORS 1")  !call etall(s2%V%i,ND2)
    ! CALL DACOPD(S1%I,S2%v%I)
    !    else
    !       if(.NOT. ASSOCIATED(s2%V(1)%j%r)) call crap1("MAPTAYLORS 2")  !call newetall(s2%V%j,ND2)
    !    endif

    do i=1,nd2
       s2%v(i)=s1(i)
    enddo
  END SUBROUTINE MAPTAYLORS


  SUBROUTINE  TAYLORSMAP(S1,S2)
    implicit none
    integer i
    type (damap),INTENT(IN)::S2
    type (TAYLOR),INTENT(inOUT),dimension(:)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake


    ! if(old) then
   ! if(s1(1)%i==0) call crap1("TAYLORSMAP 1")  !call allocw_old(s1(1))
    ! CALL DACOPD(S2%v%I,S1%I)
    !    else
    !       if(.NOT. ASSOCIATED(s1(1)%j%r)) call crap1("TAYLORSMAP 2")  !call allocw_old(s1(1))
    !    endif

    do i=1,nd2
       s1(i)=s2%v(i)
    enddo


  END SUBROUTINE TAYLORSMAP

  SUBROUTINE  EQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    ! if(old) then
   ! if(s2%V(1)%i==0) call crap1("EQUALMAP 1")  !call etall(s2%V%i,ND2)
    ! CALL DACOPD(S1%V%I,S2%v%I)
    !    else
    !       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("EQUALMAP 2")  !call newetall(s2%V%j,ND2)
    !    endif
    do i=1,nd2
       s2%v(i)=s1%v(i)
    enddo
  END SUBROUTINE EQUALMAP

  SUBROUTINE  EQUALgMAP(S2,S1)
    implicit none
    type (gmap),INTENT(inOUT)::S2
    type (gmap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    do i=1,s1%n
       s2%v(i)=s1%v(i)
    enddo
  END SUBROUTINE EQUALgMAP

  SUBROUTINE  EQUALgMAPdamap(S2,S1)
    implicit none
    type (gmap),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    do i=1,c_%nd2
       s2%v(i)=s1%v(i)
    enddo
  END SUBROUTINE EQUALgMAPdamap

  SUBROUTINE  EQUALdamapgmap(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    type (gmap),INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    do i=1,c_%nd2
       s2%v(i)=s1%v(i)
    enddo
  END SUBROUTINE EQUALdamapgmap



  SUBROUTINE  IdentityEQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
     integer i
    IF(.NOT.C_%STABLE_DA) RETURN

    ! if(old) then
   ! if(s2%V(1)%i==0) call crap1("IdentityEQUALMAP 1") !call etall(s2%V%i,ND2)
    IF(S1.EQ.1) then
    ! CALL etini(S2%V%I)
     do i=1,nd2
     s2%v(i)=1.0_dp.mono.i
     enddo
    endif
    IF(S1.EQ.0) then
!     CALL DACLRD(S2%V%I)
     do i=1,nd2
     s2%v(i)=0.0_dp
     enddo
endif
    !    else
    !       if(.NOT. ASSOCIATED(s2%V(1)%J%r))call crap1("IdentityEQUALMAP 2") ! call newetall(s2%V%j,ND2)
    !
    !       IF(S1.EQ.1) CALL NEWetini(S2%V%J)
    !       IF(S1.EQ.0) CALL NEWDACLRD(S2%V%J)
    !    endif
  END SUBROUTINE IdentityEQUALMAP

  SUBROUTINE  IdentityEQUALgMAP(S2,S1)
    implicit none
    type (gmap),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN

    do i=1,s2%n
       if(s1==1) then
          s2%v(i)=1.0_dp.mono.i
       else
          s2%v(i)=0.0_dp

       endif
    enddo

  END SUBROUTINE IdentityEQUALgMAP


  SUBROUTINE  zeroEQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    real(dp) zero_(ndim2)
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    do i=1,ndim2
       zero_(i)=0.0_dp
    enddo
    ! if(old) then
   ! if(s2%V(1)%i==0) call crap1("zeroEQUALMAP 1") !call etall(s2%V%i,ND2)
    !    else
    !       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("zeroEQUALMAP 2") !call newetall(s2%V%j,ND2)
    !    endif
    IF(S1.EQ.0.0_dp) s2=zero_
  END SUBROUTINE zeroEQUALMAP


  SUBROUTINE  EQUALVEC(S2,S1)
    implicit none
    type (vecfield),INTENT(inOUT)::S2
    type (vecfield),INTENT(IN)::S1
    integer i
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
  !  if(s2%V(1)%i==0) call crap1("EQUALVEC 1") !call etall(s2%V%i,ND2)
    do i=1,nd2
!     CALL DACOPD(S1%V%I,S2%v%I)
      s2%v(i)=s1%v(i)
    enddo
    !    else
    !       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("EQUALVEC 2") !call newetall(s2%V%j,ND2)
    !       CALL NEWDACOPD(S1%V%J,S2%v%J)
    !    endif
    s2%ifac=s1%ifac

  END SUBROUTINE EQUALVEC

  SUBROUTINE  EQUALvecpb(S2,S1)
    implicit none
    type (vecfield),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    ! if(old) then
!    if(s2%V(1)%i==0) call crap1("EQUALvecpb 1")  !call etall(s2%V%i,ND2)
    CALL DIFD(S1%h,S2%v,-1.0_dp)
    !    else
    !       if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call crap1("EQUALvecpb 2")  !call newetall(s2%V%j,ND2)
    !       CALL NEWDIFD(S1%h%J,S2%v%J,-one)
    !    endif
    s2%ifac=s1%ifac

  END SUBROUTINE EQUALvecpb

  subroutine difd(h1,v,sca)
    implicit none
    ! INVERSE OF INTD ROUTINE
    type(taylor) , intent(in):: h1
    type(taylor) h,b1
    integer i
    type(taylor) , intent(inout):: v(:)

    real(dp),intent (in) :: sca
    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(h)
    h=h1

    do i=1,nd
        v(2*i)=h.d.(2*i-1)
        b1=h.d.(2*i)
        v(2*i-1)=sca*b1
!      call dader(2*i,h,b1)
!      call   dacmu(b1,sca,v(2*i-1))
    enddo

    call kill(h)
    call kill(b1)
    return
  end subroutine difd

  SUBROUTINE  EQUALpbvec(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (vecfield),INTENT(IN)::S1

    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    ! if(old) then
  !  if(s2%h%i==0) call crap1("EQUALpbvec 1")  !call etall1(s2%h%i)
    CALL intd(S1%v,s2%h,-1.0_dp)
    !    else
    !       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALpbvec 2")  !call newetall(s2%h%J,1)
    !       CALL NEWintd(S1%v%J,s2%h%J,-one)
    !    endif
    s2%ifac=s1%ifac
  END SUBROUTINE EQUALpbvec

  real(dp) function dliet(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    if(.not.c_%stable_da) return

    dliet=0.0_dp
    do i=1,nd
       dliet=REAL(j(2*i-1)+j(2*i),kind=DP)+dliet
    enddo
    dliet=dliet+1.0_dp
    dliet=1.0_dp/dliet
    return
  end function dliet

  subroutine intd(v,h,sca)
    implicit none
    ! IF SCA=-one
    !     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:
    !
    ! IF SCA=one
    !     \VEC{V}.GRAD  = GRAD H . GRAD
    integer i
    type(taylor), intent(inout) :: h
    type(taylor),intent(in) ::v(ndim2)
    type(taylor) x(ndim2)

    type(taylor) b1,b2,b3,b4
    real(dp) sca


    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(b2)
    call alloc(b3)
    call alloc(b4)
    call alloc(x,nd2)

    do i=1,nd2
      x(i)=1.0_dp.mono.i
    enddo

    do i=1,nd
       call cfu(v(2*i-1),dliet,b3)
       call cfu(v(2*i),dliet,b1)
       b2=b1*x(2*i-1)
       b1=b3*x(2*i)
       b3=b2+sca*b1
       b4=b3+b4
    enddo
    h=b4

    call kill(x,nd2)
    call kill(b4)
    call kill(b3)
    call kill(b2)
    call kill(b1)

    return
  end subroutine intd


  SUBROUTINE  EQUALpbpb(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
 !   if(s2%h%i==0) call crap1("EQUALpbpb 1")  ! call etall1(s2%h%i)
!    CALL dacop(S1%h%i,S2%h%i)
      s2%h=s1%h
    !    else
    !       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALpbpb 2")  !call newetall(s2%h%J,1)
    !       CALL NEWdacop(S1%h%J,S2%h%J)
    !    endif
    s2%ifac=s1%ifac
  END SUBROUTINE EQUALpbpb


  SUBROUTINE  EQUALpbda(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (taylor),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    ! if(old) then
!    if(s2%h%i==0)  call crap1("EQUALpbda 1")  !call allocw_old(s2%h )
    !    else
    !       if(.NOT. ASSOCIATED(s2%h%J%r))  call crap1("EQUALpbda 2")  !call allocw_old(s2%h )
    !    endif
    s2%h=s1
  END SUBROUTINE EQUALpbda

  SUBROUTINE  EQUALdapb(S2,S1)
    implicit none
    type (taylor),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1

    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    ! if(old) then
!    if(s2%i==0) call crap1("EQUALdapb 1")  ! call allocw_old(s2)
    !      CALL dacop(S1%h%i,S2%i)
    !    else
    !       if(.NOT. ASSOCIATED(s2%J%r)) call crap1("EQUALdapb 2")  ! call allocw_old(s2)
    !    endif
    s2=s1%h

  END SUBROUTINE EQUALdapb


  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (damap) CUTORDER
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(CUTORDER)

    DO I=1,ND2
       CUTORDER%V(I)=(S1%V(I)).cut.S2
    ENDDO

    master=localmaster

  END FUNCTION CUTORDER

  FUNCTION CUTORDERg( S1, S2 )
    implicit none
    TYPE (gmap) CUTORDERg
    TYPE (gmap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call assdamap(CUTORDERg)

    DO I=1,ND2
       CUTORDERg%V(I)=(S1%V(I)).cut.S2
    ENDDO

    master=localmaster

  END FUNCTION CUTORDERg

  FUNCTION CUTORDERPB( S1, S2 )
    implicit none
    TYPE (PBFIELD) CUTORDERPB
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(CUTORDERPB)

    CUTORDERPB%H=(S1%H).cut.S2
    master=localmaster
    CUTORDERPB%ifac=S1%ifac

  END FUNCTION CUTORDERPB

  FUNCTION CUTORDERVEC( S1, S2 )
    implicit none
    TYPE (VECFIELD) CUTORDERVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(CUTORDERVEC)

    DO I=1,ND2
       CUTORDERVEC%V(I)=(S1%V(I)).cut.S2
    ENDDO
    CUTORDERVEC%ifac=S1%ifac
    master=localmaster

  END FUNCTION CUTORDERVEC

  FUNCTION GETORDERMAP( S1, S2 )
    implicit none
    TYPE (damap) GETORDERMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(GETORDERMAP)

    DO I=1,ND2
       GETORDERMAP%V(I)=(S1%V(I)).SUB.S2
    ENDDO

    master=localmaster

  END FUNCTION GETORDERMAP

  FUNCTION GETORDERgMAP( S1, S2 )
    implicit none
    TYPE (gmap) GETORDERgMAP
    TYPE (gmap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call assdamap(GETORDERgMAP)

    DO I=1,ND2
       GETORDERgMAP%V(I)=(S1%V(I)).SUB.S2
    ENDDO

    master=localmaster

  END FUNCTION GETORDERgMAP

  FUNCTION GETORDERVEC( S1, S2 )
    implicit none
    TYPE (VECFIELD) GETORDERVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
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
    IF(.NOT.C_%STABLE_DA) RETURN
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

    pushmatrixr=0.0_dp

    do i=1,nd2
       junk(i)=0.0_dp
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
    real(dp) push1pol,push1polt(1)
    real(dp) s2t(lnv)
    integer si
    type(taylor) s1t(1)
     if(old_package) then
        call ppush1(S1%i,s2,push1pol)
    else
        call alloc(s1t(1))
        s1t(1)=s1
        s2t=0
        s2t(1:size(s2))=s2
! GTPSA REMOVED !           call mad_tpsa_eval(1,S1t%j,nv,s2t,push1polt)
        call kill(s1t(1))
        push1pol=push1polt(1)
    endif


  END FUNCTION push1pol

  FUNCTION push1polslow( S1, S2 )
    implicit none
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), intent(in),dimension(:)::s2
    real(dp) push1polslow,push1polt(1)
    TYPE (taylor) t
    real(dp) s2t(lnv)
    integer si
    type(taylor) s1t(1)
     if(old_package) then

    call alloc(t)

    call maketree(s1,t)

    call ppush1(t%i,s2,push1polslow)
    else
        call alloc(s1t(1))
        s1t(1)=s1
        s2t=0
        s2t(1:size(s2))=s2
! GTPSA REMOVED !           call mad_tpsa_eval(1,S1t%j,nv,s2t,push1polt)
        call kill(s1t(1))
        push1polslow=push1polt(1)
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
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(concat)
    call alloc(t1);call alloc(t2);call alloc(tempnew);


    do i=1,ndim2
       v1(i)=0.0_dp
       zero_(i)=0.0_dp
    enddo
    t1=s1;t2=s2;
    v1=s1     ! change oct 2004.10
    t1=zero_;t2=zero_;
     if(old_package) then
    call etcct(t1%v%i,t2%v%i,tempnew%v%i)
    call dacopd(tempnew%v%i,concat%v%i)
        else
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(nd2,t1%v,nd2,t2%v,t2%v)
concat=t2
    !       call NEWetcct(t1%v%J,t2%v%J,tempnew%v%j)
    !       call NEWdacopd(tempnew%v%j,concat%v%J)
        endif
    concat=v1

    call kill(t1);call kill(t2);call kill(tempnew);
    master=localmaster

  END FUNCTION concat

  FUNCTION concatg(S1,S2)
    implicit none
    TYPE (gmap) concatg,t1,t2,tempnew
    TYPE (gmap), INTENT (IN) :: S1, S2
    real(dp) v1(lnv),zero_(lnv)
    integer i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    concatg%n=s1%n
    call assdamap(concatg)
    call alloc(t1,S1%n);call alloc(t2,S1%n);call alloc(tempnew,S1%n);


    v1=0.0_dp
    zero_=0.0_dp

    t1=s1;t2=s2;
    v1=s2
    t1=zero_;t2=zero_;
     if(old_package) then
    call getcct(t1%v%i,t2%v%i,tempnew%v%i,s1%n)
    do i=1,s1%n
       call dacop(tempnew%v(i)%i,concatg%v(i)%i)
    enddo
        else
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(s1%n,t1%v,s2%n,t2%v,concatg%v)
        endif
    concatg=v1

    call kill(t1);call kill(t2);call kill(tempnew);
    master=localmaster

  END FUNCTION concatg

    subroutine fpp_mad_tpsa_compose(n1,t1,n2,t2,t3)
    implicit none
    type(taylor), intent(in) :: t1(:),t2(:)
    type(taylor), intent(inout) :: t3(:)
    type(taylor) t2t(lnv)
    integer i,n1,n2

    call alloc(t2t(1:nv),nv)

     do i=1,min(size(t2),nv)
       t2t(i)=t2(i)
     enddo

     do i=min(size(t2),nv)+1,nv
       t2t(i)=1.0_dp.mono.i
     enddo
! GTPSA REMOVED !     call   mad_tpsa_compose(n1,t1(1:n1)%j,nv,t2t(1:nv)%j,t3(1:n1)%j)

    call KILL(t2t(1:nv),nv)

    end subroutine fpp_mad_tpsa_compose

  FUNCTION concator( S1, S2 )
    implicit none
    TYPE (damap) concator
    TYPE (damap), INTENT (IN) :: S1, S2
    TYPE (damap) tempnew
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call alloc(tempnew)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(concator)

     if(old_package) then
    call etcct(s1%v%i,s2%v%i,tempnew%v%i)
    call dacopd(tempnew%v%i,concator%v%i)
        else
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(nd2,s1%v,nd2,s2%v,concator%v)
        endif

    master=localmaster
    call kill(tempnew)

  END FUNCTION concator

  FUNCTION concatorg( S1, S2 )
    implicit none
    TYPE (gmap) concatorg
    TYPE (gmap), INTENT (IN) :: S1, S2
    TYPE (gmap) tempnew
    integer localmaster,i
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call alloc(tempnew)
    concatorg%n=s1%n
    call assdamap(concatorg)

   if(old_package) then
    call getcct(s1%v%i,s2%v%i,tempnew%v%i,s1%n)
    do i=1,s1%n
       call dacop(tempnew%v(i)%i,concatorg%v(i)%i)
    enddo
        else
    !       call NEWgetcct(s1%v%J,s2%v%J,tempnew%v%j,s1%n)
! GTPSA REMOVED ! !       call mad_tpsa_compose(s1%n,s1%v%j,s1%n,s2%v%j,concatorg%v%j)
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(s1%n,s1%v,s2%n,s2%v,concatorg%v)
    !       do i=1,s1%n
    !          call newdacop(tempnew%v(i)%j,concatorg%v(i)%j)
    !       enddo
        endif

    master=localmaster
    call kill(tempnew)

  END FUNCTION concatorg


  FUNCTION trxflow(S2,S1)
    implicit none
    TYPE (vecfield) trxflow
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    type(damap) s22,tempnew
    real(dp) zero_(ndim2)
    integer i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call alloc(s22,tempnew)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxflow)

    do i=1,nd2
       zero_(i)=0.0_dp
    enddo
    s22=s2
    s22=zero_
     if(old_package) then
    call trxflo(s1%v%i,tempnew%v%i,s22%v%i)
    call dacopd(tempnew%v%i,trxflow%v%i)
        else
    call trxflo_g(s1%v,tempnew%v,s22%v)
    call dacopd_g(tempnew%v,trxflow%v)
        endif

    call kill(s22,tempnew)
    master=localmaster
    trxflow%IFAC=S1%IFAC
  END FUNCTION trxflow

  subroutine trxflo_g(h,rh,y)
    implicit none
    !  *RH* = Y *H* Y^-1  CHANGE OF A VECTOR FLOW OPERATOR
    integer j,k
    type(taylor) b1,b2
    type(taylor),dimension(:)::h,rh,y
    type(taylor),dimension(ndim2)::yi,ht
    type(taylor) s22g(2)
    if(.not.c_%stable_da) return

    call alloc(yi,nd2)  !  ,'YI        ')
    call alloc(ht,nd2)  !  ,'HT        ')
    call alloc(b1 )
    call alloc(b2 )
    call alloc(s22g,2 )

    !call etinv(y,yi)
! GTPSA REMOVED !       call fpp_mad_tpsa_minv(nd2,y,yi)
    !----- HT= H o Y
 !   call etcct(h,y,ht)
! GTPSA REMOVED !    call fpp_mad_tpsa_compose(nd2,h,nd2,y,ht)
    !----
    do j=1,nd2
    rh(j)=0.0_dp
    enddo
 !   call daclrd(rh)
    do j=1,nd2
       do k=1,nd2
        b1=yi(j).d.k
       !   call dader(k,yi(j),b1)
       s22g(1:1)=b1
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(1,s22g(1:1),nd2,y,s22g(2:2))
         ! call trx(b1,b2,y)
            b1=b2*ht(k)
       !   call damul(b2,ht(k),b1)
            b2=b1+rh(j)
      !    call daadd(b1,rh(j),b2)
            rh(j)=b2
      !    call dacop(b2,rh(j))
       enddo
    enddo

    call kill(s22g,2 )
    call kill(b2)
    call kill(b1)
    call kill(ht,nd2)
    call kill(yi,nd2)
    return
  end subroutine trxflo_g

  FUNCTION trxpb( S2, S1 )
    implicit none
    TYPE (pbfield) trxpb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    !    TYPE (damap) S22
    !    real(dp) zero_(ndim2)
    !    integer i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    !    call alloc(s22)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxpb)


    trxpb%h=s1%h*s2

    master=localmaster
    trxpb%ifac=S1%ifac

  END FUNCTION trxpb


  FUNCTION trxtaylor( S1, S2 )
    implicit none
    TYPE (taylor) trxtaylor
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap)  S22,s22g
    real(dp) zero_(ndim2)
    integer i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    do i=1,nd2
       zero_(i)=0.0_dp
    enddo


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxtaylor)

    call alloc(s22,s22g)

    s22=s2
    s22=zero_
     if(old) then
        call trx(s1%i,temp%i,s22%v%i)
        call dacop(temp%i,trxtaylor%i)
    else
        s22g%v(1)=s1
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(1,s22g%v(1:1),nd2,s22%v,s22g%v(2:2))
        trxtaylor=s22g%v(2)
    endif

    call kill(s22,s22g)
    master=localmaster

  END FUNCTION trxtaylor

  FUNCTION trxgtaylor( S1, S2 )
    implicit none
    TYPE (taylor) trxgtaylor
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (gmap), INTENT (IN) ::  S2
    TYPE (gmap)  S22,s22g
    real(dp),allocatable:: zero_(:)
    integer i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master



    call assdamap(trxgtaylor)

    call alloc(s22,s2%n)

    allocate(zero_(S22%n))

    do i=1,S22%n
       zero_(i)=0.0_dp
    enddo


    s22=s2
    s22=zero_
    if(old_package) then
    call gtrx(s1%i,temp%i,s22%v%i,s2%n)
    call dacop(temp%i,trxgtaylor%i)
     else
        call alloc(s22g,s2%n)
        s22g%v(1)=s1
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(1,s22g%v(1:1),s2%n,s22%v,s22g%v(2:2))
        trxgtaylor=s22g%v(2)
        call kill(s22g)
    endif
    call kill(s22)
    deallocate(zero_)
    master=localmaster

  END FUNCTION trxgtaylor

  FUNCTION trxtaylorc( S1, S2 )
    implicit none
    TYPE (taylor) trxtaylorc
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap) s22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxtaylorc)

    if(old_package) then
    call trx(s1%i,temp%i,s2%v%i)
    call dacop(temp%i,trxtaylorc%i)
        else
        call alloc(s22)
        s22%v(1)=s1
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(1,s22%v(1:1),nd2,s2%v,s22%v(2:2))
        trxtaylorc=s22%v(2)
        call kill(s22)
       endif

    master=localmaster

  END FUNCTION trxtaylorc

  FUNCTION trxgtaylorc( S1, S2 )
    implicit none
    TYPE (taylor) trxgtaylorc
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (gmap), INTENT (IN) ::  S2
    TYPE (damap) s22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call assdamap(trxgtaylorc)

    if(old_package) then
     call gtrx(s1%i,temp%i,s2%v%i,s2%n)
     call dacop(temp%i,trxgtaylorc%i)
    else
        call alloc(s22)
        s22%v(1)=s1
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(1,s22%v(1:1),s2%n,s2%v,s22%v(2:2))
        trxgtaylorc=s22%v(2)
        call kill(s22)
    endif

    master=localmaster

  END FUNCTION trxgtaylorc



  FUNCTION texpdf( S1, S2,NRMIN,NRMAX,SCA,IFAC )
    implicit none
    TYPE (damap) texpdf
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    integer NRMIN,NRMAX,IFAC
    TYPE (damap) tempnew
    real(dp) sca
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(texpdf)
    call alloc(tempnew)

     if(old_package) then
    !write(6,*) NRMIN,NRMAX,SCA,IFAC
    call  FACFLOD(s1%v%i,s2%v%i,tempnew%v%i,NRMIN,NRMAX,SCA,IFAC )
    call dacopd(tempnew%v%i,texpdf%v%i)
        else
         call  FACFLOD_g(s1%v,s2%v,texpdf%v,NRMIN,NRMAX,SCA,IFAC )
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
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(texpdft)

     if(old_package) then
    !write(6,*) NRMIN,NRMAX,SCA,IFAC
    call  FACFLO(s1%v%i,s2%i,temp%i,NRMIN,NRMAX,SCA,IFAC )
    call dacop(temp%i,texpdft%i)
         else
     call  FACFLO_g(s1%v ,s2 ,texpdft ,NRMIN,NRMAX,SCA,IFAC )
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
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    no1=no

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(explieflo)
    call alloc(tempnew)
    if(old_package) then
      if(s1%ifac/=0) then
         explieflo=texpdf( S1, S2,2,no1,1.0_dp,s1%ifac )
      else
         call EXPFLOD(s1%v%i,s2%v%i,tempnew%v%i,eps_tpsalie,nrmax)
         call dacopd(tempnew%v%i,explieflo%v%i)

      endif
    else
      if(s1%ifac/=0) then
         explieflo=texpdf( S1, S2,2,no1,1.0_dp,s1%ifac )
      else

         call EXPFLOD_g(s1%v,s2%v,tempnew%v,eps_tpsalie,nrmax)
         call dacopd_g(tempnew%v,explieflo%v)

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
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    no1=no
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(expflot)

    if(old_package) then
      if(s1%ifac/=0) then
         expflot=texpdft( S1, S2,2,no1,1.0_dp,s1%ifac )
      else
         ! if(old) then
         call EXPFLO(s1%v%i,s2%i,temp%i,eps_tpsalie,nrmax)
         call dacop(temp%i,expflot%i)
      endif
    else
      if(s1%ifac/=0) then
         expflot=texpdft( S1, S2,2,no1,1.0_dp,s1%ifac )
      else
         ! if(old) then
         call EXPFLO_g(s1%v,s2,temp,eps_tpsalie,nrmax)
         expflot=temp
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
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call alloc(T)

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(expliepb)
    T=S1
    expliepb=explieflo(T,S2)

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
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call alloc(T)

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(exppb)
    T=S1
    exppb=expfloT(T,S2)

    CALL KILL(T)
    master=localmaster

  END FUNCTION exppb


  FUNCTION ADDMAP( S1, S2 )
    implicit none
    TYPE (damap) ADDMAP
    TYPE (damap), INTENT (IN) :: S1, S2
    integer localmaster,i
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(ADDMAP)


       do i=1,nd2
         ADDMAP%v(i)=s1%v(i)+s2%v(i)
       enddo



    master=localmaster

  END FUNCTION ADDMAP

  FUNCTION SUBMAP( S1, S2 )
    implicit none
    TYPE (damap) SUBMAP
    TYPE (damap), INTENT (IN) :: S1, S2
    integer localmaster,i
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(SUBMAP)



       do i=1,nd2
         SUBMAP%v(i)=s1%v(i)-s2%v(i)
       enddo



    master=localmaster

  END FUNCTION SUBMAP

  FUNCTION DMULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  DMULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I
    integer localmaster

    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(DMULMAPsc)


       do i=1,nd2
         DMULMAPsc%v(i)=sc*s1%v(i)
       enddo


    master=localmaster


  END FUNCTION DMULMAPsc

  FUNCTION MULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  MULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I

    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(MULMAPsc)



       do i=1,nd2
         MULMAPsc%v(i)=sc*s1%v(i)
       enddo

    master=localmaster


  END FUNCTION MULMAPsc

  FUNCTION IMULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  IMULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(IMULMAPsc)



       do i=1,nd2
         IMULMAPsc%v(i)=sc*s1%v(i)
       enddo

    master=localmaster


  END FUNCTION IMULMAPsc

  FUNCTION scDMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scDMULMAP
    TYPE (damap), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I

    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master

    call checkdamap(s1)
    call assdamap(scDMULMAP)




       do i=1,nd2
         scDMULMAP%v(i)=sc*s1%v(i)
       enddo

    master=localmaster


  END FUNCTION scDMULMAP

  FUNCTION scMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scMULMAP
    TYPE (damap), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I
    integer localmaster

    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(scMULMAP)



        do i=1,nd2
         scMULMAP%v(i)=sc*s1%v(i)
       enddo

    master=localmaster


  END FUNCTION scMULMAP

  FUNCTION scIMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scIMULMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I
    integer localmaster

    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)
    call assdamap(scIMULMAP)


        do i=1,nd2
         scIMULMAP%v(i)=sc*s1%v(i)
       enddo

    master=localmaster


  END FUNCTION scIMULMAP

    subroutine fpp_mad_tpsa_minv(ns,tpsa_a,tpsa_r)
    implicit none
    type(taylor) :: tpsa_a(:),tpsa_r(:)
    type(taylor)  a(lnv),b(lnv)
    integer i,n,ns

    call alloc(a,nv)
    call alloc(b,nv)

    do i=1,ns
     a(i)=tpsa_a(i)
    enddo
    do i=ns+1,nv
        a(i)=1.0_dp.mono.i
    enddo

! GTPSA REMOVED !    call mad_tpsa_minv(nv,a(1:nv)%j,b(1:nv)%j)
    do i=1,ns
     tpsa_r(i)=b(i)
    enddo
    call kill(a,nv)
    call kill(b,nv)

    end subroutine fpp_mad_tpsa_minv

    subroutine fpp_mad_tpsa_pminv(ns,tpsa_a,tpsa_r,select)
    implicit none
    type(taylor) :: tpsa_a(:),tpsa_r(:)
    type(taylor)  a(lnv),b(lnv)
    integer :: select(:)
    integer i,n,ns,j(lnv)

    call alloc(a,nv)
    call alloc(b,nv)

    do i=1,ns
     a(i)=tpsa_a(i)
    enddo
    do i=ns+1,nv
        a(i)=1.0_dp.mono.i
    enddo
    j(1:size(select))=select

! GTPSA REMOVED !    call mad_tpsa_pminv(nv,a(1:nv)%j,b(1:nv)%j,j(1:nv))
    do i=1,ns
     tpsa_r(i)=b(i)
    enddo
    call kill(a,nv)
    call kill(b,nv)

    end subroutine fpp_mad_tpsa_pminv

  FUNCTION POWMAP( S1, R2 )
    implicit none
    TYPE (damap) POWMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (damap) S11
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    call checkdamap(s1)

    call assdamap(POWMAP)

    call alloc(s11)

    s11=1


    R22=IABS(R2)
    DO I=1,R22
       s11=s1*s11
    ENDDO

    IF(R2.LT.0) THEN
        if(old_package) then
       CALL etinv(S11%v%i,S11%v%i)

               else
! GTPSA REMOVED !       call fpp_mad_tpsa_minv(nd2,S11%v,S11%v)
               endif
    ENDIF

    powmap=s11


    ! powmap=junk
    call kill(s11)

    master=localmaster

  END FUNCTION POWMAP

  FUNCTION gPOWMAP( S1, R2 )
    implicit none
    TYPE (gmap) gPOWMAP
    TYPE (gmap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (gmap) S11
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    gPOWMAP%n=s1%n
    call assdamap(gPOWMAP)

    call alloc(s11,s1%n)

    s11=1


    R22=IABS(R2)
    DO I=1,R22
       s11=s11*s1
    ENDDO

    IF(R2.LT.0) THEN
        if(old_package) then
       CALL getinv(S11%v%i,S11%v%i,s11%n)
              else
! GTPSA REMOVED !              call fpp_mad_tpsa_minv(s11%n,S11%v,S11%v)
              endif
    ENDIF

    gpowmap=s11


    ! powmap=junk
    call kill(s11)

    master=localmaster

  END FUNCTION gPOWMAP

  FUNCTION gPOWMAPtpsa( S1, R2 )
    implicit none
    TYPE (gmap) gPOWMAPtpsa
    TYPE (gmap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (gmap) S11,s0
    INTEGER I,R22
    integer localmaster
    real(dp), allocatable :: v(:)

    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    gPOWMAPtpsa%n=s1%n
    call assdamap(gPOWMAPtpsa)

    call alloc(s11,s1%n)
    call alloc(s0,s1%n)

    s11=1


    R22=IABS(R2)
    DO I=1,R22
       s11=s11.o.s1
    ENDDO
    allocate(v(s1%n))

    DO I=1,s1%n
       v(i)=s11%v(i).sub.'0'
    ENDDO


    IF(R2.LT.0) THEN
        if(old_package) then
       CALL getinv(S11%v%i,S11%v%i,s11%n)
              else
! GTPSA REMOVED !              call fpp_mad_tpsa_minv(s11%n,S11%v,S11%v)
              endif
    ENDIF

    do i=1,s1%n
       s0%v(i)=(1.0_dp.mono.i)-v(i)
       !s11%v(i)=s11%v(i)-(s11%v(i).sub.'0')
    enddo

    !     s0=s11.o.s0
    s11=s11.o.s0

    do i=1,s1%n
       gPOWMAPtpsa%v(i)=s11%v(i) !+ (s0%v(i).sub.'0')
    enddo






    ! powmap=junk
    call kill(s11,s0)
    deallocate(v)

    master=localmaster

  END FUNCTION gPOWMAPtpsa

  FUNCTION POWMAP_INV( S1, R2 )
    implicit none
    TYPE (damap) POWMAP_INV
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2(:)
    TYPE (damap) S11
    INTEGER I,jn(lnv)
    integer localmaster
    IF(.NOT.C_%STABLE_DA) RETURN
    localmaster=master


    do i=1,lnv
       jn(i)=0
    enddo
    do i=1,nd2
       jn(i)=R2(I)
    enddo
    call checkdamap(s1)

    call assdamap(POWMAP_INV)

    call alloc(s11)

     if(old_package) then
    if(s1%v(1)%i==0) call crap1("POWMAP_INV 2")  !call etall(s2%m%v%i,nd2)
    call etpin(S1%V%i,S11%v%i,jn)
        else
! GTPSA REMOVED !       call fpp_mad_tpsa_pminv(nd2,S1%V,S11%V,jn)
       endif

    POWMAP_INV=s11


    ! powmap=junk
    call kill(s11)

    master=localmaster

  END FUNCTION POWMAP_INV

  subroutine checksymp(s1,norm,orthogonal)
    implicit none
    TYPE (damap) s1
    real(dp)  norm1,mat(8,8),xj(8,8)
    real(dp), optional :: norm
    integer i,j
    logical(lp), optional :: orthogonal
    logical(lp) nn
    ! checks symplectic conditions on linear map
    nn=.not.present(norm)
    mat=0.d0
    mat=s1
    xj=0.d0
    if(present(orthogonal)) then
       if(orthogonal) then
          do i=1,nd
             xj(2*i-1,2*i-1)=1.0_dp
             xj(2*i,2*i)=1.0_dp
          enddo
       else
          do i=1,nd
             xj(2*i-1,2*i)=1.0_dp
             xj(2*i,2*i-1)=-1.0_dp
          enddo
       endif
    else
       do i=1,nd
          xj(2*i-1,2*i)=1.0_dp
          xj(2*i,2*i-1)=-1.0_dp
       enddo
    endif
    xj= MATMUL( transpose(mat),MATMUL(xj,mat))

    norm1=0.d0
    do i=1,nd2
       if(lielib_print(9)==1.or.nn)     write(6,'(6(1x,E15.8))') xj(i,1:nd2)
       do j=1,nd2

          norm1=norm1+abs(xj(i,j))
       enddo
    enddo
    norm1=abs(norm1-nd2)
    if(lielib_print(9)==1.or.nn) write(6,'(a29,(1x,E15.8))')"deviation from symplecticity ", norm1
    if(present(norm)) norm=abs(norm1)

  end subroutine checksymp


  subroutine checkmap(s1)
    implicit none
    TYPE (damap) s1
    integer i
     if(old_package) then
    do i=1,nd2
       if(s1%v(i)%i==0) then

            write(6,*) "Should not be here: checkmap"

          ! call !write_e(200)
          !             s1%v(i)%i=dummymap(i)
       endif
    enddo
         else

       endif

  end subroutine checkmap

  subroutine checkvec(s1)
    implicit none
    TYPE (vecfield) s1
    integer i

     if(old_package) then
    do i=1,nd2
       if(s1%v(i)%i==0) then

            write(6,*)  "Should not be here: checkvec"

       endif
    enddo
        else
    !       do i=1,nd2
    !          if(.NOT. ASSOCIATED(s1%v(i)%J%r)) then
    !             w_p=0
    !             w_p%nc=1
    !             w_p=(/"Should not be here: Assign variables"/)
    !             w_p%fc='(1((1X,A72),/))'
    !             ! call !write_e(203)
    !             !             s1%v(i)%j=dummymapl(i)
    !          endif
    !       enddo
        endif

  end subroutine checkvec

  subroutine checkpb(s1)
    implicit none
    TYPE (pbfield) s1

     if(old_package) then
    if(s1%h%i==0) then

         write(6,*) "Should not be here: Assign variables checkbp"

    endif
       else
    !       if(.NOT. ASSOCIATED(s1%h%J%r)) then
    !          !          s1%h%j=dummyl
    !          w_p=0
    !          w_p%nc=1
    !          w_p=(/"Should not be here: Assign variables"/)
    !          w_p%fc='(1((1X,A72),/))'
    !          ! call !write_e(205)
    !       endif
       endif
  end subroutine checkpb

  subroutine checktaylor(s1)
    implicit none
    TYPE (taylor) s1
    ! if(old) then
    if(old_package) then
         if(s1%i==0) then

            write(6,*) "Should not be here: Assign variables checktaylor "

         endif
        else

      endif
  end subroutine checktaylor



  subroutine assPB(s1)
    implicit none
    TYPE (PBFIELD) s1

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)

         write(6,*) " cannot indent anymore "

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

         write(6,*) "cannot indent anymore : asstaylor"

       ! call !write_e(100)
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

         write(6,*) " cannot indent anymore assvec"

       ! call !write_e(100)
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

         write(6,*) " cannot indent anymore assmap "

       ! call !write_e(100)
    end select

    do i=1,nd2
       call ass0(s1%v(i))
    enddo

  end subroutine assmap

  subroutine assgmap(s1)
    implicit none
    TYPE (gmap) s1
    integer i

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)

         write(6,*) " cannot indent anymore assgmap"

       ! call !write_e(100)
    end select

    do i=1,s1%n
       call ass0(s1%v(i))
    enddo

  end subroutine assgmap

  !Radiation
!  SUBROUTINE  allocrad(S1)
!    implicit none
!    type (radtaylor),INTENT(INOUT)::S1
!    ! if(old) then
!    call etall1(s1%v%i)
!    call etall(s1%e%i,nd2)
!    !    else
!    !       call newetall(s1%v%j,1)
!    !       call newetall(s1%e%j,nd2)
!    !   endif
!  END SUBROUTINE allocrad
!
!  SUBROUTINE  KILLrad(S1)
!    implicit none
!    type (radtaylor),INTENT(INOUT)::S1
!    ! if(old) then
!    call dadal1(s1%v%i)
!    call dadal(s1%e%i,nd2)
!    !    else
!    !       call newdadal(s1%v%j,1)
!    !       call newdadal(s1%e%j,nd2)
!    !   endif
!
!  END SUBROUTINE KILLrad

!  SUBROUTINE  allocrads(S1,n)
!    implicit none
!    type (radtaylor),INTENT(INOUT),dimension(:)::S1
!    integer i ,n
!    do i=1,n
!       call allocrad(s1(i))
!    enddo
!
!  END SUBROUTINE allocrads
!
!  SUBROUTINE  KILLrads(S1,n)
!    implicit none
!    type (radtaylor),INTENT(INOUT),dimension(:)::s1
!    integer i,n
!    do i=1,n
!       call killrad(s1(i))
!    enddo
!
!
!  END SUBROUTINE KILLrads



!  SUBROUTINE  EQUALrad(S2,S1)
!    implicit none
!    type (TAYLOR),INTENT(inOUT)::S2
!    type (radTAYLOR),INTENT(IN)::S1
!    !     if(iass0>ndum) then
!    !     call ndum_warning
!    !     endif
!    !      iass0=0
!    ! if(old) then
!    if(s2%i==0) call crap1("EQUALrad 1")  ! call allocw_old(s2)
!    CALL DACOP(S1%v%I,S2%I)
!    !    else
!    !       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUALrad 2")  !call allocw_old(s2)
!    !       call newdacop(S1%v%j,S2%j)
!    !   endif
!  END SUBROUTINE EQUALrad


!  SUBROUTINE  radEQUAL(S2,S1)
!    implicit none
!    type (radTAYLOR),INTENT(inOUT)::S2
!    type (TAYLOR),INTENT(IN)::S1
!    call check_snake
!
!    !     if(iass0>ndum) then
!    !     call ndum_warning
!    !     endif
!    !      iass0=0
!    ! if(old) then
!    if(s2%v%i==0) call crap1("radEQUAL 1")  ! call allocw_old(s2%v)
!    CALL DACOP(S1%I,S2%v%I)
!    !    else
!    !       IF (.NOT. ASSOCIATED(s2%v%j%r)) call crap1("radEQUAL 2")  !call allocw_old(s2%v)
!    !       call newdacop(S1%j,S2%v%j)
!    !   endif
!  END SUBROUTINE radEQUAL
!
!  ! remove small numbers


!!!!! Laurent GTPSA various 'FLO' routines

  subroutine daflo_g(h,x,y)
    implicit none
    ! LIE EXPONENT ROUTINES WITH FLOW OPERATORS
    !
    !     \VEC{H}.GRAD X =Y
    integer i
    type(taylor) x,y,b1,b2,b3
    type(taylor),dimension(:)::h
    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(b2)
    call alloc(b3)


    do i=1,nd2
        b2=x.d.i
        b3=b2*h(i)
        b2=b3+b1
        b1=b2
    enddo

    y=b1
    call kill(b3)
    call kill(b2)
    call kill(b1)
    return
  end subroutine daflo_g

  subroutine daflod_g(h,x,y)
    implicit none
    integer i
    type(taylor),dimension(:)::h,x,y
    type(taylor),dimension(ndim2)::b1,b2
    if(.not.c_%stable_da) return

    call alloc(b1,nd2)
    call alloc(b2,nd2)
    b1=h
    b2=x


    do i=1,nd2
       call daflo_g(b1,b2(i),y(i))
    enddo

    call kill(b1,nd2)
    call kill(b2,nd2)
    return
  end subroutine daflod_g

  subroutine expflo_g(h,x,y,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y
    logical(lp) more
    integer i,nrmax
     type(taylor) x,y,b1,b2,b3,b4
    type(taylor) ,dimension(:)::h
    real(dp) coe,eps,r,rbefore
    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(b2)
    call alloc(b3)
    call alloc(b4)
    b4=x
    b1=x

    more=.true.
    rbefore=1e30_dp
    do i=1,nrmax
       coe=1.0_dp/REAL(i,kind=DP)
        b2=b1*coe
       call daflo_g(h,b2,b1)
        b3=b1+b4
        r=full_abs(b1)

       if(more) then
          if(r.gt.eps) then
             rbefore=r
             goto 100
          else
             rbefore=r
             more=.false.
          endif
       else
          if(r.ge.rbefore) then
            y=b3

             call kill(b4)
             call kill(b3)
             call kill(b2)
             call kill(b1)
             return
          endif
          rbefore=r
       endif
100    continue
        b4=b3

    enddo
    if(lielib_print(2)==1) then
       write(6,'(a6,1x,G21.14,1x,a25)') ' NORM ',eps,' NEVER REACHED IN EXPFLO '
    endif
     y=b3
    call kill(b4)
    call kill(b3)
    call kill(b2)
    call kill(b1)
    return
  end subroutine expflo_g

  subroutine facflod_g(h,x,w,nrmin,nrmax,sca,ifac)
    implicit none
    ! IFAC=1
    ! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX )  \VEC{X}= \VEC{Y}
    ! IFAC=-1
    ! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) \VEC{X}= \VEC{Y}
    integer i,ifac,nrmax,nrmin
    type(taylor),dimension(:)::x,w,h
    real(dp) sca
    if(.not.c_%stable_da) return

    do i=1,nd2
       call facflo_g(h,x(i),w(i),nrmin,nrmax,sca,ifac)
    enddo

    return
  end subroutine facflod_g

  subroutine facflo_g(h,x,w,nrmin,nrmax,sca,ifac)
    implicit none
    ! IFAC=1
    ! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX ) X= Y
    ! IFAC=-1
    ! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) X= Y
    integer i,ifac,nmax,nrmax,nrmin
    type(taylor) x,w,v
    type(taylor) ,dimension(:)::h
    type(taylor) ,dimension(ndim2)::bm,b0
    real(dp) eps,sca
    if(.not.c_%stable_da) return

    call alloc(bm,nd2) !  ,'BM        ')
    call alloc(b0,nd2) !  ,'B0        ')
    call alloc(v)
    v=x
   ! call dacop(x,v)

    !    eps=-one
    !    call daeps(eps)
    eps=epsflo
    nmax=100
    !
    ! IFAC =1 ---> V = EXP(:SCA*H(NRMAX):)...EXP(:SCA*H(NRMIN):)X
    if(ifac.eq.1) then
       do i=nrmax,nrmin,-1
          call taked_g(h,i,b0)

          call dacmud_g(b0,sca,bm)

          call expflo_g(bm,v,b0(1),eps,nmax)
        v=b0(1)
!          call dacop(b0(1),v)
       enddo
    else
       ! IFAC =-1 ---> V = EXP(:SCA*H(NRMIN):)...EXP(:SCA*H(NRMAX):)X
       do i=nrmin,nrmax
          call taked_g(h,i,b0)
          call dacmud_g(b0,sca,bm)

          call expflo_g(bm,v,b0(1),eps,nmax)
            v=b0(1)
!          call dacop(b0(1),v)
       enddo
    endif
!    call dacop(v,w)
        w=v
    call kill(v)
    call kill(b0,nd2)
    call kill(bm,nd2)
    return
  end subroutine facflo_g

  subroutine taked_g(h,m,ht)
    implicit none
    !  \VEC{HT}= \VEC{H_M}  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
    integer i,m
    type(taylor) b1,b2
    type(taylor),dimension(:)::h,ht
    type(taylor),dimension(ndim2)::x
    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(b2)
    call alloc(x,nd2) !  ,'X         ')



    do   i=1,nd2
       ht(i)=h(i).sub.m
 !      call take(h(i),m,ht(i))
    enddo
    call kill(x,nd2)
    call kill(b2)
    call kill(b1)
    return
  end subroutine taked_g

  subroutine dacmud_g(h,sca,ht)
    implicit none
    integer i
    type(taylor),dimension(:)::h,ht
    real(dp) sca
    if(.not.c_%stable_da) return

    do i=1,nd2
       ht(i)=sca*h(i)  ! call dacmu(h(i),sca,ht(i))
    enddo
    return
  end subroutine dacmud_g

  subroutine expflod_g(h,x,w,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) \VEC{X} = \VEC{Y}
    integer j,nrmax
    type(taylor) b0
    type(taylor),dimension(:)::x,w,h
    type(taylor),dimension(ndim2)::v
    real(dp) eps
    if(.not.c_%stable_da) return

    call alloc(b0 )
    call alloc(v,nd2) !  ,'V         ')
    v=x
  !  call dacopd(x,v)
    do j=1,nd2
       call expflo_g(h,v(j),b0,eps,nrmax)
        v(j)=b0
 !      call dacop(b0,v(j))
    enddo
    call dacopd_g(v,w)
    call kill(v,nd2)
    call kill(b0)
    return
  end subroutine expflod_g


  subroutine dacopd_g(h,ht)
    implicit none
    !    H goes into HT  (nd2 array)
    integer i
    type(taylor),dimension(:)::h,ht
    if(.not.c_%stable_da) return

    do i=1,nd2
       ht(i)=h(i)
    enddo
    return
  end subroutine dacopd_g
end module tpsalie
