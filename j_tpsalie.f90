!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size
module tpsalie
  use tpsa
  implicit none
  private  ASSVEC,ASSMAP,ASSPB,asstaylor,explieflo,expliepb,expflot,exppb
  private  checkmap,checkpb,checkvec,checktaylor,MAPmatrixr,matrixMAPr,TREEMAP
  private  TAYLORSMAP,DPEKMAP,DPOKMAP,zeroEQUALMAP,IdentyEQUALMAP
  private DABSMAP,ABSMAP,EQUALMAP,EQUALVEC    !,EQUALMAPVEC,EQUALVECMAP
  !private EQUALvectaylor,EQUALtaylorvec
  !private DPEKFIELD,DPOKFIELD
  !private MAPmatrix,matrixMAP
  !PRIVATE DADDVECFIELD,DVECFIELDADD,DADDMAP,DMAPADD
  private EQUALvecpb,EQUALpbpb,EQUALpbvec,EQUALpbda,EQUALdapb
  private GETORDERVEC,GETORDERMAP,GETORDERPB,concator,pushtree,pushmap,concat,pushmatrixr,push1polslow
  private trxflow,trxpb,trxtaylor,DMULMAPsc,MULMAPsc,IMULMAPsc,DMULVECsc,MULVECsc,IMULVECsc
  private DMULpbsc,MULpbsc,IMULpbsc
  private scDMULMAP,scMULMAP,scIMULMAP,scDMULVEC,scMULVEC,scIMULVEC,scDMULpb,scMULpb,scIMULpb
  private ADDMAP,VECMAP,MAPVEC,VECPB,PBVEC
  private SUBMAP,POWMAP,DAREADTAYLORS,DAREADMAP,DAREADVEC,DAREApb,DAREADTAYLOR
  PRIVATE DAPRINTTAYLORS,DAPRINTMAP
  private DAPRINTVEC,DAPRINTTAYLOR,DAPRINTpb,allocTAYLORS,allocmap,allocvec,allocpb,alloctree,allocrad,allocrads
  private KILLTAYLORS,KILLmap,KILLvec,KILLpb,KILLtree,killrad,killrads
  !private push1pol
  integer,private::NO,ND,ND2,NP,NDPT,NV
  logical(lp),private::old

  !frs real(dp),private::eps=c_1e_9
  integer::nrmax=400
  integer,dimension(ndim2)::dummymap,tempnew
  TYPE (TAYLORLOW),dimension(ndim2)::tempLD,dummymapl


  integer iassmapdol

  type(dummapping) dummap(ndum)

  TYPE damap
     type (taylor) v(ndim2)
  END TYPE damap

  TYPE vecfield
     type (taylor) v(ndim2)
     integer ifac
  END TYPE vecfield

  TYPE pbfield
     type (taylor) h
     integer ifac
  END TYPE pbfield


  !TYPE matrix
  !real(dp) M(ndim2,ndim2)     ! no really used
  !END TYPE matrix

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
     MODULE PROCEDURE IdentyEQUALMAP
     MODULE PROCEDURE zeroEQUALMAP
     MODULE PROCEDURE DABSMAP
     MODULE PROCEDURE ABSMAP
     MODULE PROCEDURE EQUALVEC
     MODULE PROCEDURE EQUALpbpb
     !MODULE PROCEDURE EQUALMAPVEC   ! obsolete
     !MODULE PROCEDURE EQUALVECMAP ! obsolete
     MODULE PROCEDURE EQUALpbda
     MODULE PROCEDURE EQUALdapb
     !MODULE PROCEDURE EQUALvectaylor
     !MODULE PROCEDURE EQUALtaylorvec
     MODULE PROCEDURE DPEKMAP
     MODULE PROCEDURE DPOKMAP
     !MODULE PROCEDURE DPEKFIELD
     !MODULE PROCEDURE DPOKFIELD
     MODULE PROCEDURE EQUALvecpb
     MODULE PROCEDURE EQUALpbvec
     !MODULE PROCEDURE MAPmatrix
     !MODULE PROCEDURE matrixMAP
     MODULE PROCEDURE MAPmatrixr
     MODULE PROCEDURE matrixMAPr
     MODULE PROCEDURE TREEMAP
     !radiation
     MODULE PROCEDURE radEQUAL
     MODULE PROCEDURE EQUALrad
  end  INTERFACE

  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDERVEC
     MODULE PROCEDURE GETORDERMAP
     MODULE PROCEDURE GETORDERPB
  end  INTERFACE

  !INTERFACE OPERATOR (.i.)
  !MODULE PROCEDURE GETindexMAP    ! Useless out perhaps permanently
  !MODULE PROCEDURE GETindexvec
  !end  INTERFACE


  INTERFACE OPERATOR (.o.)
     MODULE PROCEDURE concator
     MODULE PROCEDURE trxtaylorc
     MODULE PROCEDURE trxpbc
     MODULE PROCEDURE trxflowc
  end  INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE pushmap
     MODULE PROCEDURE pushmatrixr
     !MODULE PROCEDURE push1pol
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

     MODULE PROCEDURE DMULVECsc
     MODULE PROCEDURE MULVECsc
     MODULE PROCEDURE IMULVECsc
     MODULE PROCEDURE scDMULVEC
     MODULE PROCEDURE scMULVEC
     MODULE PROCEDURE scIMULVEC

     MODULE PROCEDURE DMULpbsc
     MODULE PROCEDURE MULpbsc
     MODULE PROCEDURE IMULpbsc
     MODULE PROCEDURE scDMULpb
     MODULE PROCEDURE scMULpb
     MODULE PROCEDURE scIMULpb
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE ADDMAP
     MODULE PROCEDURE VECMAP
     MODULE PROCEDURE MAPVEC
     MODULE PROCEDURE VECPB
     MODULE PROCEDURE PBVEC
     !MODULE PROCEDURE ADDID
     ! MODULE PROCEDURE IDADD      ! Adding identity Useless
     ! MODULE PROCEDURE DADDMAP
     ! MODULE PROCEDURE DMAPADD
     ! MODULE PROCEDURE DADDVECFIELD
     ! MODULE PROCEDURE DVECFIELDADD
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
     MODULE PROCEDURE allocTAYLORS
     MODULE PROCEDURE allocmap
     MODULE PROCEDURE allocvec
     MODULE PROCEDURE allocpb
     MODULE PROCEDURE alloctree
     !radiation
     MODULE PROCEDURE allocrad
     MODULE PROCEDURE allocrads
  END INTERFACE

  INTERFACE allocTPSA
     MODULE PROCEDURE allocTAYLORS
     MODULE PROCEDURE allocmap
     MODULE PROCEDURE allocvec
     MODULE PROCEDURE allocpb
     !  MODULE PROCEDURE allocgen
     MODULE PROCEDURE alloctree
     !radiation
     MODULE PROCEDURE allocrad
     MODULE PROCEDURE allocrads
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILLTAYLORS
     MODULE PROCEDURE KILLmap
     MODULE PROCEDURE KILLvec
     MODULE PROCEDURE KILLpb
     MODULE PROCEDURE KILLtree
     !radiation
     MODULE PROCEDURE KILLrad
     MODULE PROCEDURE KILLrads
  END INTERFACE
  INTERFACE KILLTPSA
     MODULE PROCEDURE KILLTAYLORS
     MODULE PROCEDURE KILLmap
     MODULE PROCEDURE KILLvec
     MODULE PROCEDURE KILLpb
     !  MODULE PROCEDURE KILLgen
     MODULE PROCEDURE KILLtree
     !radiation
     MODULE PROCEDURE KILLrad
     MODULE PROCEDURE KILLrads
  END INTERFACE



contains

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



  SUBROUTINE  allocTAYLORS(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1(ndim2)
    if(old) then
       call etall(s1%i,nd2)
    else
       call NEWetall(s1%j,nd2)
    endif
  END SUBROUTINE allocTAYLORS

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



  SUBROUTINE  KILLTAYLORS(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1(NDIM2)
    if(old) then
       call DADAL(s1%i,nd2)
    else
       call newDADAL(s1%j,nd2)
    endif
  END SUBROUTINE KILLTAYLORS

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

  SUBROUTINE  DABSMAP(R2,S1)
    implicit none
    real(dp),INTENT(inOUT)::R2
    type (damap),INTENT(IN)::S1
    real(dp) R1
    INTEGER I
    R2=zero
    R1=zero
    if(old) then
       if(s1%V(1)%i==0) call crap1("DABSMAP 1")  !call etall(s1%V%i,ND2)

       DO i=1,ND2
          !          R1=s1%V(I) !2002.10.17
          R1=abs(s1%V(I)) !2002.10.17
          R2=R2+R1
       ENDDO
    else
       if(.NOT. ASSOCIATED(s1%V(1)%j%r))call crap1("DABSMAP 2")  ! call newetall(s1%V%j,ND2)

       DO i=1,ND2
          !          R1=s1%V(I)
          R1=abs(s1%V(I))
          R2=R2+R1
       ENDDO
    endif

  END SUBROUTINE DABSMAP

  SUBROUTINE  ABSMAP(R2,S1)
    implicit none
    REAL(SP),INTENT(inOUT)::R2
    type (damap),INTENT(IN)::S1
    REAL(sp) R1
    INTEGER I
    R2=zero
    R1=zero
    if(old) then
       if(s1%V(1)%i==0) call crap1("ABSMAP 1")  !call etall(s1%V%i,ND2)

       DO i=1,ND2
          !          R1=s1%V(I) !2002.10.17
          R1=abs(s1%V(I)) !2002.10.17
          R2=R2+R1
       ENDDO
    else
       if(.NOT. ASSOCIATED(s1%V(1)%j%r)) call crap1("ABSMAP 2")  !call newetall(s1%V%j,ND2)

       DO i=1,ND2
          !          R1=s1%V(I)
          R1=abs(s1%V(I))
          R2=R2+R1
       ENDDO
    endif

  END SUBROUTINE ABSMAP

  SUBROUTINE  DPEKMAP(S2,S1)
    implicit none
    real(dp),INTENT(inOUT),dimension(:)::S2
    type (damap),INTENT(IN)::S1
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

  !SUBROUTINE  DPEKFIELD(S2,S1)
  !    real(dp),INTENT(OUT)::S2(NDIM2)
  !    type (VECFIELD),INTENT(IN)::S1
  !     if(old) then
  !      CALL DAPEK0(S1%V%I,S2,nd2)
  !     else
  !      CALL NEWDAPEK0(S1%V%J,S2,nd2)
  !     endif
  !END SUBROUTINE DPEKFIELD

  !SUBROUTINE  DPOKFIELD(S1,S2)
  !    real(dp),INTENT(IN)::S2(NDIM2)
  !    type (VECFIELD),INTENT(OUT)::S1
  !     if(old) then
  !      if(s1%V(1)%i==0) call etall(s1%V%i,ND2)
  !      CALL DAPOK0(S1%V%I,S2,nd2)
  !     else
  !      if(.NOT. ASSOCIATED(s1%V(1)%j%r)) call newetall(s1%V%j,ND2)
  !      CALL NEWDAPOK0(S1%V%J,S2,nd2)
  !     endif
  !END SUBROUTINE DPOKFIELD

  SUBROUTINE  TREEMAP(S1,S2)
    implicit none
    type (damap),INTENT(IN)::S2
    type (TREE),INTENT(inOUT)::S1
    integer i

    do i=1,nd2
       call maketree(S2%v(i),S1%branch(i))
    enddo

  END SUBROUTINE TREEMAP

  !SUBROUTINE  matrixMAP(S2,S1)
  !    type (matrix),INTENT(OUT)::S2
  !    type (damap),INTENT(IN)::S1
  !     integer jj(40),i,j,JL(100)
  !     real(dp) m(ndim2,ndim2)
  !     if(old) then
  !      do i=1,nd2
  !      do j=1,nd2
  !      jj(j)=1
  !      call dapek(S1%v(i)%i,jj,m(i,j))
  !      jj(j)=0
  !      enddo
  !      enddo
  !      do i=1,nd2
  !      do j=1,nd2
  !      s2%m(i,j)=m(i,j)
  !      enddo
  !      enddo
  !     else
  !      do i=1,nd2
  !      do j=1,nd2
  !      jL(j)=1
  !      call NEWdapek(S1%v(i)%J,jL,m(i,j))
  !      jL(j)=0
  !      enddo
  !      enddo
  !      do i=1,nd2
  !      do j=1,nd2
  !      s2%m(i,j)=m(i,j)
  !      enddo
  !      enddo
  !     endif
  !END SUBROUTINE matrixMAP

  !SUBROUTINE  MAPmatrix(S1,S2)
  !    type (matrix),INTENT(in)::S2
  !    type (damap),INTENT(out)::S1
  !     integer jj(40),i,j,JL(100)
  !     s1=0
  !     if(old) then
  !      do i=1,nd2
  !      do j=1,nd2
  !      jj(j)=1
  !      call dapok(S1%v(i)%i,jj,s2%m(i,j))
  !      jj(j)=0
  !      enddo
  !      enddo
  !     else
  !      do i=1,nd2
  !      do j=1,nd2
  !      jL(j)=1
  !      call NEWdapok(S1%v(i)%J,jL,s2%m(i,j))
  !      jL(j)=0
  !      enddo
  !      enddo
  !     endif
  !END SUBROUTINE MAPmatrix


  SUBROUTINE  matrixMAPr(S2,S1)
    implicit none
    real(dp),INTENT(inOUT)::S2(ndim2,ndim2)
    type (damap),INTENT(IN)::S1
    integer i,j,JL(lnv)
    real(dp) m(ndim2,ndim2)

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
    real(dp),INTENT(in)::S2(ndim2,ndim2)
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


  SUBROUTINE  IdentyEQUALMAP(S2,S1)
    implicit none
    type (damap),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    if(old) then
       if(s2%V(1)%i==0) call crap1("IdentyEQUALMAP 1") !call etall(s2%V%i,ND2)
       IF(S1.EQ.1) CALL etini(S2%V%I)
       IF(S1.EQ.0) CALL DACLRD(S2%V%I)
    else
       if(.NOT. ASSOCIATED(s2%V(1)%J%r))call crap1("IdentyEQUALMAP 2") ! call newetall(s2%V%j,ND2)

       IF(S1.EQ.1) CALL NEWetini(S2%V%J)
       IF(S1.EQ.0) CALL NEWDACLRD(S2%V%J)
    endif
  END SUBROUTINE IdentyEQUALMAP


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


  !SUBROUTINE  EQUALVECMAP(S2,S1)
  !implicit none
  !    type (vecfield),INTENT(OUT)::S2
  !    type (damap),INTENT(IN)::S1
  !     if(old) then
  !      if(s2%V(1)%i==0) call etall(s2%V%i,ND2)
  !      CALL DACOPD(S1%V%I,S2%v%I)
  !     else
  !      if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call newetall(s2%V%j,ND2)
  !      CALL NEWDACOPD(S1%V%J,S2%v%J)
  !     endif
  !END SUBROUTINE EQUALVECMAP

  !SUBROUTINE  EQUALMAPVEC(S2,S1)
  !implicit none
  !    type (damap),INTENT(OUT)::S2
  !    type (vecfield),INTENT(IN)::S1
  !     if(old) then
  !      if(s2%V(1)%i==0) call etall(s2%V%i,ND2)
  !      CALL DACOPD(S1%V%I,S2%v%I)
  !     else
  !      if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call newetall(s2%V%j,ND2)
  !      CALL NEWDACOPD(S1%V%J,S2%v%J)
  !     endif
  !END SUBROUTINE EQUALMAPVEC

  SUBROUTINE  EQUALVEC(S2,S1)
    implicit none
    type (vecfield),INTENT(inOUT)::S2
    type (vecfield),INTENT(IN)::S1
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
    if(old) then
       if(s2%h%i==0) call crap1("EQUALpbvec 1")  !call etall1(s2%h%i)
       CALL intd(S1%v%i,s2%h%i,-one)
    else
       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALpbvec 2")  !call newetall(s2%h%J,1)
       CALL NEWintd(S1%v%J,s2%h%J,-one)
    endif
    s2%ifac=s1%ifac
  END SUBROUTINE EQUALpbvec

  !SUBROUTINE  EQUALvectaylor(S2,S1)
  !    type (vecfield),INTENT(OUT)::S2
  !    type (taylor),INTENT(IN)::S1
  !     if(old) then
  !      if(s2%V(1)%i==0) call etall(s2%V%i,ND2)
  !      CALL DIFD(S1%i,S2%v%i,-one)
  !     else
  !      if(.NOT. ASSOCIATED(s2%V(1)%J%r)) call newetall(s2%V%j,ND2)
  !      CALL NEWDIFD(S1%J,S2%v%J,-one)
  !     endif
  !END SUBROUTINE EQUALvectaylor

  !SUBROUTINE  EQUALtaylorvec(S2,S1)
  !    type (taylor),INTENT(OUT)::S2
  !    type (vecfield),INTENT(IN)::S1
  !     if(old) then
  !      if(s2%i==0) call etall1(s2%i)
  !      CALL intd(S1%v%i,s2%i,-one)
  !     else
  !      if(.NOT. ASSOCIATED(s2%J%r)) call newetall(s2%J,1)
  !      CALL NEWintd(S1%v%J,s2%J,-one)
  !     endif
  !END SUBROUTINE EQUALtaylorvec

  SUBROUTINE  EQUALpbpb(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (pbfield),INTENT(IN)::S1
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
    !     if(old) then
    !      if(s2%h%i==0) call etall1(s2%h%i)
    !      CALL dacop(S1%i,S2%h%i)
    !     else
    !      if(.NOT. ASSOCIATED(s2%h%J%r)) call newetall(s2%h%J,1)
    !      CALL NEWdacop(S1%J,S2%h%J)
    !     endif

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
    if(old) then
       if(s2%i==0) call crap1("EQUALdapb 1")  ! call allocw_old(s2)
       !      CALL dacop(S1%h%i,S2%i)
    else
       if(.NOT. ASSOCIATED(s2%J%r)) call crap1("EQUALdapb 2")  ! call allocw_old(s2)
       !      CALL NEWdacop(S1%h%J,S2%J)
    endif
    s2=s1%h
  END SUBROUTINE EQUALdapb


  !FUNCTION GETindexMAP( S1, S2 )
  !implicit none
  !  TYPE (taylor) GETindexMAP
  !  TYPE (damap), INTENT (IN) :: S1
  !  INTEGER, INTENT (IN) :: S2
  !
  !  call checkdamap(s1)
  !  call ass(GETindexMAP)
  !
  !      GETindexMAP=S1%V(S2)
  !
  !END FUNCTION GETindexMAP

  !FUNCTION GETindexvec( S1, S2 )
  !implicit none
  !  TYPE (taylor) GETindexvec
  !  TYPE (VECFIELD), INTENT (IN) :: S1
  !  INTEGER, INTENT (IN) :: S2
  !
  !  call checkdamap(s1)
  !  call ass(GETindexvec)
  !
  !      GETindexvec=S1%V(S2)
  !
  !END FUNCTION GETindexvec


  FUNCTION GETORDERMAP( S1, S2 )
    implicit none
    TYPE (damap) GETORDERMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I

    call checkdamap(s1)
    call assdamap(GETORDERMAP)

    DO I=1,ND2
       GETORDERMAP%V(I)=(S1%V(I)).SUB.S2
    ENDDO


  END FUNCTION GETORDERMAP

  FUNCTION GETORDERVEC( S1, S2 )
    implicit none
    TYPE (VECFIELD) GETORDERVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I

    call checkdamap(s1)
    call assdamap(GETORDERVEC)

    DO I=1,ND2
       GETORDERVEC%V(I)=(S1%V(I)).SUB.S2
    ENDDO


  END FUNCTION GETORDERVEC

  FUNCTION GETORDERPB( S1, S2 )
    implicit none
    TYPE (PBFIELD) GETORDERPB
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2


    call checkdamap(s1)
    call assdamap(GETORDERPB)


    GETORDERPB%H=(S1%H).SUB.S2


  END FUNCTION GETORDERPB

  FUNCTION GETdiffMAP( S1, S2 )
    implicit none
    TYPE (damap) GETdiffMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I

    call checkdamap(s1)
    call assdamap(GETdiffMAP)

    DO I=1,ND2
       GETdiffMAP%V(I)=S1%V(I).d.S2
    ENDDO


  END FUNCTION GETdiffMAP

  FUNCTION GETdiffVEC( S1, S2 )
    implicit none
    TYPE (VECFIELD) GETdiffVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2
    INTEGER I

    call checkdamap(s1)
    call assdamap(GETdiffVEC)

    DO I=1,ND2
       GETdiffVEC%V(I)=S1%V(I).d.S2
    ENDDO


  END FUNCTION GETdiffVEC

  FUNCTION GETdiffPB( S1, S2 )
    implicit none
    TYPE (PBFIELD) GETdiffPB
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: S2


    call checkdamap(s1)
    call assdamap(GETdiffPB)


    GETdiffPB%H=S1%H.d.S2


  END FUNCTION GETdiffPB


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
    TYPE (damap) concat,t1,t2
    TYPE (damap), INTENT (IN) :: S1, S2
    real(dp) v1(ndim2),zero_(ndim2)
    integer i
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(concat)
    call alloc(t1);call alloc(t2);


    do i=1,ndim2
       v1(i)=zero
       zero_(i)=zero
    enddo
    t1=s1;t2=s2;
    v1=s2
    t1=zero_;t2=zero_;
    if(old) then
       call etcct(t1%v%i,t2%v%i,tempnew)
       call dacopd(tempnew,concat%v%i)
    else
       call NEWetcct(t1%v%J,t2%v%J,tempLd)
       call NEWdacopd(tempLd,concat%v%J)
    endif
    concat=v1

    call kill(t1);call kill(t2);

  END FUNCTION concat

  FUNCTION concator( S1, S2 )
    implicit none
    TYPE (damap) concator
    TYPE (damap), INTENT (IN) :: S1, S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(concator)

    if(old) then
       call etcct(s1%v%i,s2%v%i,tempnew)
       call dacopd(tempnew,concator%v%i)
    else
       call NEWetcct(s1%v%J,s2%v%J,tempLd)
       call NEWdacopd(tempLd,concator%v%J)
    endif


  END FUNCTION concator


  FUNCTION trxflow(S1,S2)
    implicit none
    TYPE (vecfield) trxflow
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    type(damap) s22
    real(dp) zero_(ndim2)
    integer i

    call alloc(s22)
    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxflow)
    do i=1,nd2
       zero_(i)=zero
    enddo
    s22=s2
    s22=zero_
    if(old) then
       call trxflo(s1%v%i,tempnew,s22%v%i)
       call dacopd(tempnew,trxflow%v%i)
    else
       call NEWtrxflo(s1%v%J,tempLd,s22%v%J)
       call NEWdacopd(tempLd,trxflow%v%J)
    endif

    call kill(s22)
  END FUNCTION trxflow

  FUNCTION trxflowc(S1, S2 )
    implicit none
    TYPE (vecfield) trxflowc
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(trxflowc)

    if(old) then
       call trxflo(s1%v%i,tempnew,s2%v%i)
       call dacopd(tempnew,trxflowc%v%i)
    else
       call NEWtrxflo(s1%v%J,tempLd,s2%v%J)
       call NEWdacopd(tempLd,trxflowc%v%J)
    endif


  END FUNCTION trxflowc


  FUNCTION trxpb( S1, S2 )
    implicit none
    TYPE (pbfield) trxpb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap) S22
    real(dp) zero_(ndim2)
    integer i

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
  END FUNCTION trxpb

  FUNCTION trxpbc( S1, S2 )
    implicit none
    TYPE (pbfield) trxpbc
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2

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


  END FUNCTION trxpbc

  FUNCTION trxtaylor( S1, S2 )
    implicit none
    TYPE (taylor) trxtaylor
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    TYPE (damap)  S22
    real(dp) zero_(ndim2)
    integer i

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

  END FUNCTION trxtaylor

  FUNCTION trxtaylorc( S1, S2 )
    implicit none
    TYPE (taylor) trxtaylorc
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2

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


  END FUNCTION trxtaylorc



  FUNCTION texpdf( S1, S2,NRMIN,NRMAX,SCA,IFAC )
    implicit none
    TYPE (damap) texpdf
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    integer NRMIN,NRMAX,IFAC
    real(dp) sca

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(texpdf)

    if(old) then
       !write(6,*) NRMIN,NRMAX,SCA,IFAC
       call  FACFLOD(s1%v%i,s2%v%i,tempnew,NRMIN,NRMAX,SCA,IFAC )
       call dacopd(tempnew,texpdf%v%i)
    else
       call newFACFLOD(s1%v%j,s2%v%j,tempLD,NRMIN,NRMAX,SCA,IFAC )
       call NEWdacopd(tempLD,texpdf%v%J)
    endif

  END FUNCTION texpdf

  FUNCTION texpdft( S1, S2,NRMIN,NRMAX,SCA,IFAC )
    implicit none
    TYPE (taylor) texpdft
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (taylor), INTENT (IN) ::  S2
    integer NRMIN,NRMAX,IFAC
    real(dp) sca

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

  END FUNCTION texpdft


  FUNCTION explieflo( S1, S2 )
    implicit none
    TYPE (damap) explieflo
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    integer no1
    no1=no

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(explieflo)

    if(s1%ifac/=0) then
       explieflo=texpdf( S1, S2,2,no1,one,s1%ifac )
       !         write(6,*) "Improper usage: map is factorized "
       !     write(6,*) "Sorry, we will implement in later version "
       !     write(6,*) "Use a Dragt-Finn or Reverse-Dragt-Finn scratch variable "
    else
       if(old) then
          call EXPFLOD(s1%v%i,s2%v%i,tempnew,eps_tpsalie,nrmax)
          call dacopd(tempnew,explieflo%v%i)
       else
          call NEWEXPFLOD(s1%v%J,s2%v%J,tempLD,eps_tpsalie,nrmax)
          call NEWdacopd(tempLD,explieflo%v%J)
       endif
    endif

  END FUNCTION explieflo

  FUNCTION expflot( S1, S2 )
    implicit none
    TYPE (taylor) expflot
    TYPE (vecfield), INTENT (IN) :: S1
    TYPE (taylor), INTENT (IN) ::  S2
    integer no1
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

  END FUNCTION expflot

  FUNCTION expliepb( S1, S2 )
    implicit none
    TYPE (damap) expliepb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) ::  S2
    type (vecfield) T

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

  END FUNCTION expliepb

  FUNCTION exppb( S1, S2 )
    implicit none
    TYPE (taylor) exppb
    TYPE (pbfield), INTENT (IN) :: S1
    TYPE (taylor), INTENT (IN) ::  S2
    type (vecfield) T

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

  END FUNCTION exppb


  FUNCTION ADDMAP( S1, S2 )
    implicit none
    TYPE (damap) ADDMAP
    TYPE (damap), INTENT (IN) :: S1, S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(ADDMAP)

    if(old) then
       call DALIND(s1%v%i,one,s2%v%i,one,ADDMAP%v%i)
    else
       call newDALIND(s1%v%j,one,s2%v%j,one,ADDMAP%v%j)
    endif


  END FUNCTION ADDMAP

  FUNCTION SUBMAP( S1, S2 )
    implicit none
    TYPE (damap) SUBMAP
    TYPE (damap), INTENT (IN) :: S1, S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(SUBMAP)

    if(old) then
       call DALIND(s1%v%i,one,s2%v%i,-one,SUBMAP%v%i)
    else
       call newDALIND(s1%v%j,one,s2%v%j,-one,SUBMAP%v%j)
    endif


  END FUNCTION SUBMAP

  !FUNCTION IDADD( S1, S2 )
  !implicit none
  !  TYPE (damap) IDADD
  !  TYPE (damap), INTENT (IN) :: S2
  !  INTEGER, INTENT (IN) :: S1
  !  integer junk
  !
  !  call checkdamap(s2)
  !  call assdamap(IDADD)
  !  ! To remove warning
  !junk=s1
  !
  !
  !     if(old) then
  !      CALL ETINI(tempnew)
  !      call DALIND(s2%v%i,one,tempnew,one,IDADD%v%i)
  !     else
  !      CALL newETINI(TEMPLD)
  !      call newDALIND(s2%v%J,one,TEMPLD,one,IDADD%v%J)
  !     endif
  !
  !END FUNCTION IDADD

  !FUNCTION ADDID( S2, S1 )
  !implicit none
  !  TYPE (damap) ADDID
  !  TYPE (damap), INTENT (IN) :: S2
  !  INTEGER, INTENT (IN) :: S1
  !  integer junk
  !
  !  call checkdamap(s2)
  !  call assdamap(ADDID)
  !   ! To remove warning
  !junk=s1
  !   !
  !     if(old) then
  !      CALL ETINI(tempnew)
  !      call DALIND(s2%v%i,one,tempnew,one,ADDID%v%i)
  !     else
  !      CALL NEWETINI(TEMPLD)
  !      call NEWDALIND(s2%v%J,one,TEMPLD,one,ADDID%v%J)
  !     endif
  !
  !END FUNCTION ADDID



  FUNCTION VECMAP( S1, S2 )
    implicit none
    TYPE (VECFIELD) VECMAP
    TYPE (VECFIELD), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) :: S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(VECMAP)

    if(old) then
       call DALIND(s1%v%i,one,s2%v%i,one,tempnew)
       call dacopd(tempnew,VECMAP%v%i)
    else
       call NEWDALIND(s1%v%J,one,s2%v%J,one,tempLD)
       call NEWdacopd(tempLD,VECMAP%v%J)
    endif


  END FUNCTION VECMAP

  FUNCTION MAPVEC( S2, S1 )
    implicit none
    TYPE (VECFIELD) MAPVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    TYPE (damap), INTENT (IN) :: S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(MAPVEC)

    if(old) then
       call DALIND(s1%v%i,one,s2%v%i,one,tempnew)
       call dacopd(tempnew,MAPVEC%v%i)
    else
       call NEWDALIND(s1%v%J,one,s2%v%J,one,tempLD)
       call NEWdacopd(tempLD,MAPVEC%v%J)
    endif


  END FUNCTION MAPVEC

  FUNCTION VECPB( S1, S2 )
    implicit none
    TYPE (VECFIELD) VECPB
    TYPE (VECFIELD), INTENT (IN) :: S1
    TYPE (PBFIELD), INTENT (IN) :: S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(VECPB)

    if(old) then
       CALL DIFD(s2%h%i,tempnew,-one)
       call DALIND(s1%v%i,one,tempnew,one,VECPB%v%i)
    else
       CALL NEWDIFD(s2%h%J,TEMPLD,-one)
       call NEWDALIND(s1%v%J,one,TEMPLD,one,VECPB%v%J)
    endif
  END FUNCTION VECPB

  FUNCTION PBVEC( S2, S1 )
    implicit none
    TYPE (VECFIELD) PBVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    TYPE (PBFIELD), INTENT (IN) :: S2

    call checkdamap(s1)
    call checkdamap(s2)
    call assdamap(PBVEC)

    if(old) then
       CALL DIFD(s2%h%i,tempnew,-one)
       call DALIND(s1%v%i,one,tempnew,one,PBVEC%v%i)
    else
       CALL NEWDIFD(s2%h%J,TEMPLD,-one)
       call NEWDALIND(s1%v%J,one,TEMPLD,one,PBVEC%v%J)
    endif


  END FUNCTION PBVEC


  !FUNCTION DADDVECFIELD( S1, S2 )
  !  TYPE (VECFIELD) DADDVECFIELD
  !  TYPE (VECFIELD), INTENT (IN) :: S2
  !  real(dp), INTENT (IN) :: S1(NDIM2)
  !  INTEGER I
  !
  !  call checkdamap(s2)
  !  call assdamap(DADDVECFIELD)
  !
  !      DO I=1,ND2
  !      DADDVECFIELD%v(I)=s2%v(I)+S1(I)
  !      ENDDO
  !
  !END FUNCTION DADDVECFIELD

  !FUNCTION DVECFIELDADD(  S2,S1 )
  !  TYPE (VECFIELD) DVECFIELDADD
  !  TYPE (VECFIELD), INTENT (IN) :: S2
  !  real(dp), INTENT (IN) :: S1(NDIM2)
  !  INTEGER I
  !
  !  call checkdamap(s2)
  !  call assdamap(DVECFIELDADD)
  !
  !   DO I=1,ND2
  !   DVECFIELDADD%v(I)=s2%v(I)+S1(I)
  !   ENDDO
  !
  !END FUNCTION DVECFIELDADD

  ! FUNCTION DADDMAP( S1, S2 )
  !   implicit none
  !   TYPE (DAMAP) DADDMAP
  !   TYPE (DAMAP), INTENT (IN) :: S2
  !   real(dp), INTENT (IN),dimension(:)::S1
  !   INTEGER I
  !
  !   call checkdamap(s2)
  !   call assdamap(DADDMAP)
  !
  !   DO I=1,ND2
  !      DADDMAP%v(I)=s2%v(I)+S1(I)
  !   ENDDO
  !
  !END FUNCTION DADDMAP

  ! FUNCTION DMAPADD(  S2,S1 )
  !   implicit none
  !   TYPE (DAMAP) DMAPADD
  !   TYPE (DAMAP), INTENT (IN) :: S2
  !   real(dp), INTENT (IN),dimension(:)::S1
  !   INTEGER I
  !
  !   call checkdamap(s2)
  !   call assdamap(DMAPADD)
  !
  !   DO I=1,ND2
  !      DMAPADD%v(I)=s2%v(I)+S1(I)
  !   ENDDO
  !
  ! END FUNCTION DMAPADD

  FUNCTION DMULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  DMULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(DMULMAPsc)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,SC,tempnew(I))
       ENDDO
       call dacopd(tempnew,DMULMAPsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,SC,TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,DMULMAPsc%v%J)
    endif


  END FUNCTION DMULMAPsc

  FUNCTION MULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  MULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(MULMAPsc)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,MULMAPsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,MULMAPsc%v%J)
    endif

  END FUNCTION MULMAPsc

  FUNCTION IMULMAPsc( S1, Sc )
    implicit none
    TYPE (damap)  IMULMAPsc
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(IMULMAPsc)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,IMULMAPsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,IMULMAPsc%v%J)
    endif

  END FUNCTION IMULMAPsc

  FUNCTION scDMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scDMULMAP
    TYPE (damap), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(scDMULMAP)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,SC,tempnew(I))
       ENDDO
       call dacopd(tempnew,scDMULMAP%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,SC,TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,scDMULMAP%v%J)
    endif

  END FUNCTION scDMULMAP

  FUNCTION scMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scMULMAP
    TYPE (damap), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(scMULMAP)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,scMULMAP%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,scMULMAP%v%J)
    endif

  END FUNCTION scMULMAP

  FUNCTION scIMULMAP( Sc , S1 )
    implicit none
    TYPE (damap)  scIMULMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(scIMULMAP)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,scIMULMAP%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,scIMULMAP%v%J)
    endif

  END FUNCTION scIMULMAP


  FUNCTION DMULVECsc( S1, Sc )
    implicit none
    TYPE (VECFIELD)  DMULVECsc
    TYPE (VECFIELD), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(DMULVECsc)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,SC,tempnew(I))
       ENDDO
       call dacopd(tempnew,DMULVECsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,SC,TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,DMULVECsc%v%J)
    endif

  END FUNCTION DMULVECsc

  FUNCTION MULVECsc( S1, Sc )
    implicit none
    TYPE (VECFIELD)  MULVECsc
    TYPE (VECFIELD), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(MULVECsc)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,MULVECsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,MULVECsc%v%J)
    endif


  END FUNCTION MULVECsc

  FUNCTION IMULVECsc( S1, Sc )
    implicit none
    TYPE (VECFIELD)  IMULVECsc
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(IMULVECsc)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,IMULVECsc%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,IMULVECsc%v%J)
    endif

  END FUNCTION IMULVECsc

  FUNCTION scDMULVEC( Sc, S1 )
    implicit none
    TYPE (VECFIELD)  scDMULVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(scDMULVEC)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,SC,tempnew(I))
       ENDDO
       call dacopd(tempnew,scDMULVEC%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,SC,TEMPLD(I))
       ENDDO
       call NEWdacopd(tempLD,scDMULVEC%v%J)
    endif

  END FUNCTION scDMULVEC

  FUNCTION scMULVEC( Sc, S1 )
    implicit none
    TYPE (VECFIELD)  scMULVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(scMULVEC)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,scMULVEC%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(TEMPLD,scMULVEC%v%J)
    endif

  END FUNCTION scMULVEC

  FUNCTION scIMULVEC( Sc, S1 )
    implicit none
    TYPE (VECFIELD)  scIMULVEC
    TYPE (VECFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC
    INTEGER I

    call checkdamap(s1)
    call assdamap(scIMULVEC)

    if(old) then
       DO I=1,ND2
          CALL DACMU(s1%v(I)%i,REAL(SC,kind=DP),tempnew(I))
       ENDDO
       call dacopd(tempnew,scIMULVEC%v%i)
    else
       DO I=1,ND2
          CALL NEWDACMU(s1%v(I)%J,REAL(SC,kind=DP),TEMPLD(I))
       ENDDO
       call NEWdacopd(TEMPLD,scIMULVEC%v%J)
    endif

  END FUNCTION scIMULVEC


  FUNCTION DMULpbsc( S1, Sc )
    implicit none
    TYPE (PBFIELD)  DMULpbsc
    TYPE (PBFIELD), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC

    call checkdamap(s1)
    call assdamap(DMULpbsc)

    if(old) then
       CALL DACMU(s1%H%i,SC,temp)
       call dacop(temp,DMULpbsc%H%i)
    else
       CALL NEWDACMU(s1%H%J,SC,tempL)
       call NEWdacop(tempL,DMULpbsc%H%J)
    endif

  END FUNCTION DMULpbsc

  FUNCTION MULpbsc( S1, Sc )
    implicit none
    TYPE (PBFIELD)  MULpbsc
    TYPE (PBFIELD), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC

    call checkdamap(s1)
    call assdamap(MULpbsc)

    if(old) then
       CALL DACMU(s1%H%i,REAL(SC,kind=DP),temp)
       call dacop(temp,MULpbsc%H%i)
    else
       CALL NEWDACMU(s1%H%J,REAL(SC,kind=DP),tempL)
       call NEWdacop(tempL,MULpbsc%H%J)
    endif

  END FUNCTION MULpbsc

  FUNCTION IMULpbsc( S1, Sc )
    implicit none
    TYPE (PBFIELD)  IMULpbsc
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC

    call checkdamap(s1)
    call assdamap(IMULpbsc)

    if(old) then
       CALL DACMU(s1%H%i,REAL(SC,kind=DP),temp)
       call dacop(temp,IMULpbsc%H%i)
    else
       CALL NEWDACMU(s1%H%J,REAL(SC,kind=DP),tempL)
       call NEWdacop(tempL,IMULpbsc%H%J)
    endif

  END FUNCTION IMULpbsc


  FUNCTION scDMULpb( SC, S1 )
    implicit none
    TYPE (PBFIELD)  scDMULpb
    TYPE (PBFIELD), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: SC

    call checkdamap(s1)
    call assdamap(scDMULpb)

    if(old) then
       CALL DACMU(s1%H%i,SC,temp)
       call dacop(temp,scDMULpb%H%i)
    else
       CALL NEWDACMU(s1%H%J,SC,tempL)
       call NEWdacop(tempL,scDMULpb%H%J)
    endif


  END FUNCTION scDMULpb

  FUNCTION scMULpb( SC, S1 )
    implicit none
    TYPE (PBFIELD)  scMULpb
    TYPE (PBFIELD), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: SC

    call checkdamap(s1)
    call assdamap(scMULpb)

    if(old) then
       CALL DACMU(s1%H%i,REAL(SC,kind=DP),temp)
       call dacop(temp,scMULpb%H%i)
    else
       CALL NEWDACMU(s1%H%J,REAL(SC,kind=DP),tempL)
       call NEWdacop(tempL,scMULpb%H%J)
    endif

  END FUNCTION scMULpb

  FUNCTION scIMULpb( SC, S1 )
    implicit none
    TYPE (PBFIELD)  scIMULpb
    TYPE (PBFIELD), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: SC

    call checkdamap(s1)
    call assdamap(scIMULpb)

    if(old) then
       CALL DACMU(s1%H%i,REAL(SC,kind=DP),temp)
       call dacop(temp,scIMULpb%H%i)
    else
       CALL NEWDACMU(s1%H%J,REAL(SC,kind=DP),tempL)
       call NEWdacop(tempL,scIMULpb%H%J)
    endif

  END FUNCTION scIMULpb


  FUNCTION POWMAP( S1, R2 )
    implicit none
    TYPE (damap) POWMAP
    TYPE (damap), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    TYPE (damap) S11
    INTEGER I,R22

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
  END FUNCTION POWMAP


  subroutine checkmap(s1)
    implicit none
    TYPE (damap) s1
    integer i
    if(old) then
       do i=1,nd2
          if(s1%v(i)%i==0) then
             s1%v(i)%i=dummymap(i)
          endif
       enddo
    else
       do i=1,nd2
          if(.NOT. ASSOCIATED(s1%v(i)%J%r)) then
             s1%v(i)%j=dummymapl(i)
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
             s1%v(i)%i=dummymap(i)
          endif
       enddo
    else
       do i=1,nd2
          if(.NOT. ASSOCIATED(s1%v(i)%J%r)) then
             s1%v(i)%j=dummymapl(i)
          endif
       enddo
    endif

  end subroutine checkvec

  subroutine checkpb(s1)
    implicit none
    TYPE (pbfield) s1

    if(old) then
       if(s1%h%i==0) then
          s1%h%i=dummy
       endif
    else
       if(.NOT. ASSOCIATED(s1%h%J%r)) then
          s1%h%j=dummyl
       endif
    endif
  end subroutine checkpb

  subroutine checktaylor(s1)
    implicit none
    TYPE (taylor) s1
    if(old) then
       if(s1%i==0) then
          s1%i=dummy
       endif
    else
       if(.NOT. ASSOCIATED(s1%J%r)) then
          s1%j=dummyl
       endif
    endif
  end subroutine checktaylor

  subroutine assmap(s1)
    implicit none
    TYPE (damap) s1
    integer i

    iassmapdol=iassmapdol+1
    if(iassmapdol.eq.ndum+1) iassmapdol=1

    if(old) then
       do i=1,nd2
          s1%v(i)%i=dummap(iassmapdol)%d(i)
       enddo
    else
       do i=1,nd2
          s1%v(i)%j=dummap(iassmapdol)%e(i)
       enddo
    endif

  end subroutine assmap

  subroutine assVEC(s1)
    implicit none
    TYPE (VECFIELD) s1
    integer i

    iassmapdol=iassmapdol+1
    if(iassmapdol.eq.ndum+1) iassmapdol=1

    if(old) then
       do i=1,nd2
          s1%v(i)%i=dummap(iassmapdol)%d(i)
       enddo
    else
       do i=1,nd2
          s1%v(i)%j=dummap(iassmapdol)%e(i)
       enddo
    endif

  end subroutine assVEC


  subroutine assPB(s1)
    implicit none
    TYPE (PBFIELD) s1

    call ass0(s1%h)
    !iassdol=iassdol+1
    !if(old) then
    !  s1%H%i=dum(iassdol)
    !else
    !  s1%H%j=duml(iassdol)
    !endif

  end subroutine assPB


  subroutine asstaylor(s1)
    implicit none
    TYPE (taylor) s1

    call ass0(s1)
    !iassdol=iassdol+1
    !if(iassdol.eq.ndum+1) iassdol=1
    !if(old) then
    !  s1%i=dum(iassdol)
    !else
    !  s1%j=duml(iassdol)
    !endif

  end subroutine asstaylor

  subroutine ASSIGNMAP
    implicit none
    integer i
    iassmapdol=0
    if(old) then
       CALL ETALL(DUMMYmap,nd2)
       CALL ETALL(tempnew,nd2)
       do i=1,ndum
          CALL ETALL(DUMmap(i)%d,nd2)
       enddo
    else
       CALL newETALL(DUMMYmapl,nd2)
       CALL newETALL(templd,nd2)
       do i=1,ndum
          CALL newETALL(DUMmap(i)%e,nd2)
       enddo
    endif
    !WRITE(6,*) "ASSIGNING DUMMYMAP"
  end subroutine ASSIGNMAP

  subroutine DEASSIGNMAP
    implicit none
    integer i
    iassmapdol=0
    if(old) then
       if(DUMMYmap(1)/=0) then
          CALL DADAL(DUMMYmap,nd2)
          CALL DADAL(tempnew,nd2)
          do i=1,ndum
             CALL DADAL(DUMmap(i)%d,nd2)
          enddo
       endif

    else
       CALL newDADAL(DUMMYmapl,nd2)
       CALL newDADAL(templd,nd2)
       do i=1,ndum
          CALL newDADAL(DUMmap(i)%e,nd2)
       enddo
    endif
    !WRITE(6,*) "DEASSIGNING DUMMYMAP"
  end subroutine DEASSIGNMAP

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
