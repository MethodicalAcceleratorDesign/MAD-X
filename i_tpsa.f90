!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size

MODULE TPSA
  !use newda
  use definition
  use file_handler
  IMPLICIT NONE
  integer,private::ndel ,nd2par
  integer,private,dimension(lnv)::jfil

  private equal,DAABSEQUAL,Dequaldacon ,equaldacon ,Iequaldacon  !,AABSEQUAL 2002.10.17
  private pow,powr,powr8,dlogt, GETORDER,CUTORDER,getchar,GETint
  private getdiff,getdATRA  ,mul,dmulsc,dscmul
  private mulsc,scmul,imulsc,iscmul
  private div,ddivsc,dscdiv,divsc,scdiv,idivsc,iscdiv
  private unaryADD,add,daddsc,dscadd,addsc,scadd,iaddsc,iscadd
  private  unarySUB,subs,dsubsc,dscsub,subsc,scsub,isubsc,iscsub
  private allocda,KILLda,dexpt,dcost,dsint,dsqrtt,datant,dtant
  PRIVATE GETCHARnd2,GETintnd2,dputchar,dputint, filter,check_snake,dsinHt,dCOSHt
  PRIVATE DEQUAL,REQUAL
  !  PUBLIC VAR,ASS
  private pbbra,full_absT
  PRIVATE null_0,ALLOC_U,FILL_N,REFILL_N
  private FILL_R ! new sagan

  private NO,ND,ND2,NP,NDPT,NV
  integer NP,NO,ND,ND2,NDPT,NV
  private old
  logical(lp) old
  logical(lp) :: FIRSTSCHEME=.true.
  logical(lp),target  :: real_warning =.true.

  PRIVATE null_it,Set_Up,de_Set_Up,LINE_L,RING_L,kill_DALEVEL,dealloc_DASCRATCH,set_up_level
  private insert_da,append_da


  !    scratch levels of DA using linked list

  type dascratch
     type(taylor), pointer :: t
     TYPE (dascratch),POINTER :: PREVIOUS
     TYPE (dascratch),POINTER :: NEXT
  end type dascratch

  TYPE dalevel
     INTEGER,  POINTER :: N     ! TOTAL ELEMENT IN THE CHAIN
     !
     logical(lp),POINTER ::CLOSED
     TYPE (dascratch), POINTER :: PRESENT
     TYPE (dascratch), POINTER :: END
     TYPE (dascratch), POINTER :: START
     TYPE (dascratch), POINTER :: START_GROUND ! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
     TYPE (dascratch), POINTER :: END_GROUND ! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
  END TYPE dalevel

  type(dalevel) scratchda(ndumt)

  !   end of  scratch levels of DA using linked list


  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL
     !     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
     !     MODULE PROCEDURE AABSEQUAL   ! remove 2002.10.17
     MODULE PROCEDURE DEQUAL  ! added 2002.10.17    ! check2002.10.17
     MODULE PROCEDURE REQUAL   ! added 2002.10.17   ! check2002.10.17
     MODULE PROCEDURE Dequaldacon
     MODULE PROCEDURE equaldacon
     MODULE PROCEDURE Iequaldacon
     ! UNIVERSAL_TAYLOR
     MODULE PROCEDURE null_0
     MODULE PROCEDURE FILL_N
     MODULE PROCEDURE FILL_R  ! new sagan
     MODULE PROCEDURE REFILL_N
  end  INTERFACE

  INTERFACE abs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
  END INTERFACE

  INTERFACE dabs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
  END INTERFACE


  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW
     MODULE PROCEDURE POWR
     MODULE PROCEDURE POWR8
  END INTERFACE

  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDER
     MODULE PROCEDURE getchar
     MODULE PROCEDURE GETint
  END INTERFACE


  INTERFACE OPERATOR (.CUT.)
     MODULE PROCEDURE CUTORDER
  END INTERFACE

  INTERFACE OPERATOR (.mono.)
     MODULE PROCEDURE dputchar
     MODULE PROCEDURE dputint
  END INTERFACE

  INTERFACE OPERATOR (.d.)
     MODULE PROCEDURE getdiff
  END INTERFACE

  INTERFACE OPERATOR (.K.)
     MODULE PROCEDURE getdATRA    ! Used internally primarily
  END INTERFACE

  INTERFACE OPERATOR (.pb.)
     MODULE PROCEDURE pbbra
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul
     MODULE PROCEDURE dmulsc
     MODULE PROCEDURE dscmul
     MODULE PROCEDURE mulsc
     MODULE PROCEDURE scmul
     MODULE PROCEDURE imulsc
     MODULE PROCEDURE iscmul
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div
     MODULE PROCEDURE ddivsc
     MODULE PROCEDURE dscdiv
     MODULE PROCEDURE divsc
     MODULE PROCEDURE scdiv
     MODULE PROCEDURE idivsc
     MODULE PROCEDURE iscdiv
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryADD
     MODULE PROCEDURE add
     MODULE PROCEDURE daddsc
     MODULE PROCEDURE dscadd
     MODULE PROCEDURE addsc
     MODULE PROCEDURE scadd
     MODULE PROCEDURE iaddsc
     MODULE PROCEDURE iscadd
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE unarySUB
     MODULE PROCEDURE subs
     MODULE PROCEDURE dsubsc
     MODULE PROCEDURE dscsub
     MODULE PROCEDURE subsc
     MODULE PROCEDURE scsub
     MODULE PROCEDURE isubsc
     MODULE PROCEDURE iscsub
  END INTERFACE


  INTERFACE pek
     MODULE PROCEDURE pek000 ! not private
  END INTERFACE

  INTERFACE pok
     MODULE PROCEDURE pok000  ! not private
  END INTERFACE

  INTERFACE shiftda
     MODULE PROCEDURE shift000  ! not private
  END INTERFACE

  INTERFACE var
     MODULE PROCEDURE var000  ! not private
     MODULE PROCEDURE var001  ! not private
  END INTERFACE

  INTERFACE cfu
     MODULE PROCEDURE cfu000  ! not private
  END INTERFACE


  INTERFACE ass
     MODULE PROCEDURE ass0
  END INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE allocda
     MODULE PROCEDURE allocdas
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILLda
     MODULE PROCEDURE KILLdas
  END INTERFACE

  INTERFACE alloctpsa
     MODULE PROCEDURE allocda
  END INTERFACE

  INTERFACE KILLtpsa
     MODULE PROCEDURE KILLda
  END INTERFACE



  INTERFACE full_abs
     MODULE PROCEDURE full_absT
  END INTERFACE

  INTERFACE dexp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE exp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE cexp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE cdexp
     MODULE PROCEDURE dexpt
  END INTERFACE

  INTERFACE cdcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE dcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE cos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE ccos
     MODULE PROCEDURE dcost
  END INTERFACE

  INTERFACE dcosH
     MODULE PROCEDURE dcosHt
  END INTERFACE
  INTERFACE cosH
     MODULE PROCEDURE dcosHt
  END INTERFACE

  INTERFACE cdsin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE ccsin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE dsin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE sin
     MODULE PROCEDURE dsint
  END INTERFACE

  INTERFACE dsinH
     MODULE PROCEDURE dsinHt
  END INTERFACE
  INTERFACE sinH
     MODULE PROCEDURE dsinHt
  END INTERFACE

  INTERFACE dlog
     MODULE PROCEDURE dlogt
  END INTERFACE
  INTERFACE log
     MODULE PROCEDURE dlogt
  END INTERFACE
  INTERFACE cdlog
     MODULE PROCEDURE dlogt
  END INTERFACE
  INTERFACE clog
     MODULE PROCEDURE dlogt
  END INTERFACE

  INTERFACE dsqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE
  INTERFACE sqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE

  ! NEW
  INTERFACE OPERATOR (.PAR.)
     MODULE PROCEDURE getcharnd2
     MODULE PROCEDURE GETintnd2
  END INTERFACE


  INTERFACE datan
     MODULE PROCEDURE datant
  END INTERFACE
  INTERFACE atan
     MODULE PROCEDURE datant
  END INTERFACE

  INTERFACE dtan
     MODULE PROCEDURE dtant
  END INTERFACE
  INTERFACE tan
     MODULE PROCEDURE dtant
  END INTERFACE


  INTERFACE SET_UP0
     MODULE PROCEDURE SET_UP
  END INTERFACE




CONTAINS



  subroutine set_in_tpsa(NO1,ND1,ND21,NP1,NDPT1,NV1,log)
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
  end  subroutine set_in_tpsa

  FUNCTION unaryADD( S1 )
    implicit none
    TYPE (TAYLOR) unaryADD
    TYPE (TAYLOR), INTENT (IN) :: S1

    call check(s1)
    call ass(unaryADD)

    unaryADD=s1


  END FUNCTION unaryADD

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (TAYLOR) unarySUB
    TYPE (TAYLOR), INTENT (IN) :: S1

    call check(s1)
    call ass(unarySUB)

    ! unarySUB=(-one)*s1
    if(old) then
       call dacmu(s1%i,-one,temp)
       call dacop(temp,unarySUB%i)
    else
       call newdacmu(s1%j,-one,unarySUB%j)
       !  call newdacmu(s1%j,-one,templ)
       !  call newdacop(templ,unarySUB%j)
    endif



  END FUNCTION unarySUB

  SUBROUTINE  maketree(S1,s2)
    implicit none
    type (TAYLOR),INTENT(IN)::S1
    type (TAYLOR),INTENT(inOUT):: s2
    if(old) then
       call mtree((/s1%i/),1,(/s2%i/),1)
    else
       call newdacop(s1%j,s2%j)
    endif
  END SUBROUTINE maketree

  SUBROUTINE  allocda(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1
    IF(first_time) THEN
       w_p=0
       w_p%nc=1
       w_p=(/" No TPSA package ever initialized "/)
       w_p%fc='(1((1X,A72),/))'
       CALL WRITE_E(111)
    ENDIF
    if(old) then
       s1%i=0
       call etall1(s1%i)
    else
       call allocnewda(s1%j)
    endif
  END SUBROUTINE allocda


  SUBROUTINE  allocdaS(S1,I)
    implicit none
    type (TAYLOR),INTENT(INOUT),dimension(:)::S1
    INTEGER,INTENT(IN)::I
    INTEGER J
    DO   J=1,I
       CALL ALLOCDA(S1(j))
    ENDDO
  END SUBROUTINE allocdaS

  SUBROUTINE  KILLda(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1
    if(old) then
       call DADAL1(s1%i)
    else
       call KILLNEWDAs(s1%j)
    endif

  END SUBROUTINE KILLda

  SUBROUTINE  KILLDAS(S1,I)
    implicit none
    type (TAYLOR),INTENT(INOUT),dimension(:)::S1
    INTEGER,INTENT(IN)::I
    INTEGER J
    DO   J=1,I
       CALL KILLDA(S1(j))
    ENDDO

  END SUBROUTINE KILLDAS


  SUBROUTINE  EQUAL(S2,S1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    type (TAYLOR),INTENT(IN)::S1

    call check_snake
    if(old) then
       if(s2%i==0) then
          call crap1("EQUAL 1 in tpsa") !call allocw(s2)
       endif
       if(s1%i==0) call crap1("EQUAL 2") ! call allocw(s1)
       CALL DACOP(S1%I,S2%I)
    else
       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUAL 3") !call allocw(s2)
       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("EQUAL 4") !call allocw(s1)
       call newdacop(S1%j,S2%j)
    endif
  END SUBROUTINE EQUAL



  !  SUBROUTINE  DAABSEQUAL(R1,S2)
  !    implicit none
  !    type (TAYLOR),INTENT(IN)::S2
  !    real(dp), INTENT(inOUT)::R1
  !    call check_snake
  !
  !    if(old) then
  !       CALL DAABS(S2%I,R1)
  !    else
  !       CALL newDAABS(S2%j,R1)
  !    endif
  !  END SUBROUTINE DAABSEQUAL
  !
  !  SUBROUTINE  AABSEQUAL(R1,S2)
  !    implicit none
  !    type (TAYLOR),INTENT(IN)::S2
  !    REAL(SP), INTENT(inOUT)::R1
  !    real(dp) R2
  !    if(real_warning) call real_stop
  !    call check_snake
  !
  !    if(real_warning) call real_stop
  !    r2=s2
  !    R1=R2  ! 2002.10.17
  !
  !  END SUBROUTINE AABSEQUAL

  SUBROUTINE  DEQUAL(R1,S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    real(dp), INTENT(inOUT)::R1
    call check_snake

    R1=S2.SUB.'0'
  END SUBROUTINE DEQUAL

  SUBROUTINE  REQUAL(R1,S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    REAL(SP), INTENT(inOUT)::R1

    if(real_warning) call real_stop
    call check_snake

    R1=S2.SUB.'0'

  END SUBROUTINE REQUAL

  function  DAABSEQUAL(S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    real(dp) DAABSEQUAL

    if(old) then
       CALL DAABS(S2%I,DAABSEQUAL)
    else
       CALL newDAABS(S2%j,DAABSEQUAL)
    endif
  END function DAABSEQUAL


  SUBROUTINE  DEQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    real(dp), INTENT(IN)::R1

    call check_snake

    if(old) then
       if(s2%i==0)  call crap1("DEQUALDACON 1") !call allocw(s2)
       CALL DACON(S2%I,R1)
    else
       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("DEQUALDACON 2") !call allocw(s2)
       CALL newDACON(S2%j,R1)
    endif
  END SUBROUTINE DEQUALDACON

  SUBROUTINE  EQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    REAL(SP), INTENT(IN)::R1
    real(dp) R2
    if(real_warning) call real_stop
    call check_snake

    if(real_warning) call real_stop
    if(old) then
       if(s2%i==0) call crap1("EQUALDACON 1") !call allocw(s2)
    else
       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUALDACON 2") !call allocw(s2)
    endif
    r2=REAL(r1,kind=DP)
    s2=r2
  END SUBROUTINE EQUALDACON

  SUBROUTINE  IEQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    INTEGER, INTENT(IN)::R1
    real(dp) r2
    call check_snake


    if(old) then
       if(s2%i==0) call crap1("IEQUALDACON 1") !call allocw(s2)
    else
       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("IEQUALDACON 2") !call allocw(s2)
    endif
    r2=REAL(r1,kind=DP)
    s2=r2
  END SUBROUTINE IEQUALDACON

  FUNCTION dexpt( S1 )
    implicit none
    TYPE (taylor) dexpt
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dexpt)

    if(old) then
       call dafun('EXP ',s1%i,temp)
       call dacop(temp,dexpt%i)
    else
       call newdafun('EXP ',s1%j,dexpt%j)
       !call newdafun('EXP ',s1%j,templ)
       !  call newdacop(templ,dexpt%j)
    endif


  END FUNCTION dexpt

  FUNCTION FULL_ABST( S1 )
    implicit none
    real(dp) FULL_ABST
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)

    if(old) then
       CALL DAABS(S1%I,FULL_ABST)
    else
       CALL newDAABS(S1%j,FULL_ABST)
    endif

  END FUNCTION FULL_ABST




  FUNCTION dtant( S1 )
    implicit none
    TYPE (taylor) dtant
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dtant)

    if(old) then
       call dafun('SIN ',s1%i,temp)
       call dacop(temp,dtant%i)
       call dafun('COS ',s1%i,temp)
       call dadiv(dtant%i,temp,dtant%i)
    else
       call newdafun('SIN ',s1%j,templ)
       call newdacop(templ,dtant%j)
       call newdafun('COS ',s1%j,templ)
       call newdadiv(dtant%j,templ,dtant%j)
    endif


  END FUNCTION dtant

  FUNCTION dcost( S1 )
    implicit none
    TYPE (taylor) dcost
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dcost)

    if(old) then
       call dafun('COS ',s1%i,temp)
       call dacop(temp,dcost%i)
    else
       call newdafun('COS ',s1%j,dcost%j)
       !call newdafun('COS ',s1%j,templ)
       !  call newdacop(templ,dcost%j)
    endif


  END FUNCTION dcost

  FUNCTION dsint( S1 )
    implicit none
    TYPE (taylor) dsint
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dsint)
    if(old) then
       call dafun('SIN ',s1%i,temp)
       call dacop(temp,dsint%i)
    else
       call newdafun('SIN ',s1%j,dsint%j)
       !call newdafun('SIN ',s1%j,templ)
       !  call newdacop(templ,dsint%j)
    endif


  END FUNCTION dsint

  FUNCTION dsinHt( S1 )
    implicit none
    TYPE (taylor) dsinHt
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dsinHt)
    if(old) then
       call dafun('SINH',s1%i,temp)
       call dacop(temp,dsinHt%i)
    else
       call newdafun('SINH',s1%j,dsinHt%j)
    endif

  END FUNCTION dsinHt

  FUNCTION DCOSHT( S1 )
    implicit none
    TYPE (taylor) DCOSHT
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(DCOSHT)
    if(old) then
       call dafun('COSH',s1%i,temp)
       call dacop(temp,DCOSHT%i)
    else
       call newdafun('COSH',s1%j,DCOSHT%j)
    endif


  END FUNCTION DCOSHT


  FUNCTION dlogt( S1 )
    implicit none
    TYPE (taylor) dlogt
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dlogt)
    if(old) then
       call dafun('LOG ',s1%i,temp)
       call dacop(temp,dlogt%i)
    else
       call newdafun('LOG ',s1%j,dlogt%j)
       !call newdafun('LOG ',s1%j,templ)
       !  call newdacop(templ,dlogt%j)
    endif


  END FUNCTION dlogt

  FUNCTION dsqrtt( S1 )
    implicit none
    TYPE (taylor) dsqrtt
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dsqrtt)

    if(old) then
       call dafun('SQRT',s1%i,temp)
       call dacop(temp,dsqrtt%i)
    else
       call newdafun('SQRT',s1%j,dsqrtt%j)
       !call newdafun('SQRT',s1%j,templ)
       !call newdacop(templ,dsqrtt%j)
    endif


  END FUNCTION dsqrtt

  FUNCTION dATANT( S1 )
    implicit none
    TYPE (taylor) dATANT
    TYPE (taylor), INTENT (IN) :: S1

    call check(s1)
    call ass(dATANT)
    if(old) then
       call dafun('ATAN ',s1%i,temp)
       call dacop(temp,dATANT%i)
    else
       call newdafun('ATAN ',s1%j,dATANT%j)
       !call newdafun('ATAN ',s1%j,templ)
       !call newdacop(templ,dATANT%j)
    endif


  END FUNCTION dATANT


  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (taylor) mul
    TYPE (taylor), INTENT (IN) :: S1, S2

    call check(s1)
    call check(s2)
    call ass(mul)

    if(old) then
       call damul(s1%i,s2%i,temp)
       call dacop(temp,mul%i)
    else
       call newdamul(s1%j,s2%j,mul%j)
       ! call newdamul(s1%j,s2%j,templ)
       !  call newdacop(templ,mul%j)
    endif


  END FUNCTION mul

  FUNCTION pbbra( S1, S2 )
    implicit none
    TYPE (taylor) pbbra
    TYPE (taylor), INTENT (IN) :: S1, S2

    call check(s1)
    call check(s2)
    call ass(pbbra)

    if(old) then
       call DAPOI(s1%i,s2%i,temp,nd)
       call dacop(temp,pbbra%i)
    else
       call newDAPOI(s1%j,s2%j,templ,nd)
       call newdacop(templ,pbbra%j)
    endif


  END FUNCTION pbbra

  FUNCTION GETORDER( S1, S2 )
    implicit none
    TYPE (taylor) GETORDER
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    call check(s1)
    call ass(GETORDER)

    if(old) then
       CALL TAKE(S1%I,S2,TEMP)
       call dacop(temp,GETORDER%i)
    else
       CALL NEWTAKE(S1%J,S2,TEMPL)
       call NEWdacop(tempL,GETORDER%J)
    endif

  END FUNCTION GETORDER




  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (taylor) CUTORDER
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    INTEGER I
    call check(s1)
    call ass(CUTORDER)

    if(old) then
       call dacop(S1%I,CUTORDER%i)
       DO I=S2,NO
          CALL TAKE(CUTORDER%I,I,TEMP)
          CALL DASUB(CUTORDER%I,TEMP,CUTORDER%I)
       ENDDO
    else
       call NEWdacop(S1%J,CUTORDER%J)
       DO I=S2,NO
          CALL NEWTAKE(CUTORDER%J,I,TEMPL)
          CALL NEWDASUB(CUTORDER%J,TEMPL,CUTORDER%J)
       ENDDO
    endif

  END FUNCTION CUTORDER



  FUNCTION dputchar( S1, S2 )
    implicit none
    TYPE (taylor) dputchar
    real(dp), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i

    call ass(dputchar)


    resul = trim(ADJUSTL (s2))

    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= len(trim(ADJUSTL (s2)))
    !frs    do i=1,len(trim(ADJUSTL (s2)))
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),J(I))
       if(i>nv) then
          if(j(i)>0) then
             call var(dputchar,zero,0)
             return
          endif
       endif
    enddo



    call var(dputchar,zero,0)
    CALL pok(dputchar,j,s1)

  END FUNCTION dputchar

  FUNCTION dputint( S1, S2 )
    implicit none
    TYPE (taylor) dputint
    real(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    integer j(lnv),i

    call ass(dputint)



    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= size(s2)
    !frs    do i=1,len(trim(ADJUSTL (s2)))
    do i=1,nd2par
       j(i)=s2(i)
    enddo
    do i=1,nd2par
       if(i>nv) then
          if(j(i)>0) then
             call var(dputint,zero,0)
             return
          endif
       endif
    enddo



    call var(dputint,zero,0)
    CALL pok(dputint,j,s1)

  END FUNCTION dputint

  FUNCTION GETchar( S1, S2 )
    implicit none
    real(dp) GETchar,r1
    TYPE (taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i

    call check(s1)


    resul = trim(ADJUSTL (s2))

    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= len(trim(ADJUSTL (s2)))
    !frs    do i=1,len(trim(ADJUSTl (s2)))
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),J(I))
    enddo


    if(old) then
       CALL dapek(S1%I,j,r1)
    else
       CALL newdapek(S1%j,j,r1)
    endif
    GETchar=r1

  END FUNCTION GETchar

  FUNCTION GETint( S1, S2 )
    implicit none
    real(dp) GETint,r1
    TYPE (taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2(:)
    integer j(lnv),i

    call check(s1)



    do i=1,lnv
       j(i)=0
    enddo

    !frs get around compiler problem
    nd2par= size(s2)
    !frs    do i=1,len(trim(ADJUSTl (s2)))
    do i=1,nd2par
       J(I)=s2(i)
    enddo


    if(old) then
       CALL dapek(S1%I,j,r1)
    else
       CALL newdapek(S1%j,j,r1)
    endif
    GETint=r1

  END FUNCTION GETint


  SUBROUTINE CHARINT(A,I)
    IMPLICIT NONE
    INTEGER I
    CHARACTER(1) A

    i=-1
    IF(A=='1') I=1
    IF(A=='2') I=2
    IF(A=='3') I=3
    IF(A=='4') I=4
    IF(A=='5') I=5
    IF(A=='6') I=6
    IF(A=='7') I=7
    IF(A=='8') I=8
    IF(A=='9') I=9
    IF(A=='0') I=0
    if(i==-1) ndel=1
    IF(A=='a') I=1
    IF(A=='b') I=2
    IF(A=='c') I=3
    IF(A=='d') I=4
    IF(A=='e') I=5
    IF(A=='f') I=6
    IF(A=='g') I=7
    IF(A=='h') I=8
    IF(A=='i') I=9
    IF(A==' ') I=0
    IF(A=='o') I=0
    IF(A=='A') I=1
    IF(A=='B') I=2
    IF(A=='C') I=3
    IF(A=='D') I=4
    IF(A=='E') I=5
    IF(A=='F') I=6
    IF(A=='G') I=7
    IF(A=='H') I=8
    IF(A=='I') I=9
    IF(A=='O') I=0



  END SUBROUTINE CHARINT



  FUNCTION GETdiff( S1, S2 )
    implicit none
    TYPE (taylor) GETdiff
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    call check(s1)
    call ass(GETdiff)

    if(old) then
       CALL dader(S2,S1%I,TEMP)
       call dacop(temp,GETdiff%i)
    else
       CALL NEWdader(S2,S1%J,TEMPL)
       call NEWdacop(tempL,GETdiff%J)
    endif

  END FUNCTION GETdiff

  FUNCTION GETdatra( S1, S2 )
    implicit none
    TYPE (taylor) GETdatra
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    call check(s1)
    call ass(GETdatra)

    if(old) then
       CALL datra(S2,S1%I,TEMP)
       call dacop(temp,GETdatra%i)
    else
       CALL NEWdatra(S2,S1%J,TEMPL)
       call NEWdacop(tempL,GETdatra%J)
    endif

  END FUNCTION GETdatra

  FUNCTION POW( S1, R2 )
    implicit none
    TYPE (taylor) POW
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22

    call check(s1)
    call ass(POW)

    if(old) then
       CALL DACON(TEMP,one)

       R22=IABS(R2)
       DO I=1,R22
          CALL DAMUL(TEMP,S1%I,TEMP)
       ENDDO
       IF(R2.LT.0) THEN
          CALL DADIC(TEMP,one,TEMP)
       ENDIF
       call dacop(temp,POW%i)
    ELSE

       CALL newDACON(TEMPl,one)

       R22=IABS(R2)
       DO I=1,R22
          CALL newDAMUL(TEMPl,S1%j,TEMPl)
       ENDDO
       IF(R2.LT.0) THEN
          CALL newDADIC(TEMPl,one,TEMPl)
       ENDIF
       call newdacop(templ,POW%j)
    endif
  END FUNCTION POW

  FUNCTION POWR8( S1, R2 )
    implicit none
    TYPE (taylor) POWR8
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: R2

    call check(s1)
    call ass(POWR8)

    if(old) then
       CALL DAFUN('LOG ',S1%I,TEMP)
       CALL DACMU(TEMP,R2,TEMP)
       CALL DAFUN('EXP ',TEMP,TEMP)
       call dacop(temp,POWR8%i)
    ELSE
       CALL NEWDAFUN('LOG ',S1%J,TEMPL)
       CALL NEWDACMU(TEMPL,R2,TEMPL)
       CALL NEWDAFUN('EXP ',TEMPL,POWR8%J)
       !  CALL NEWDAFUN('EXP ',TEMPL,TEMPL)
       !  call newdacop(TEMPL,POWR8%J)
    endif
  END FUNCTION POWR8

  FUNCTION POWR( S1, R2 )
    implicit none
    TYPE (taylor) POWR
    TYPE (taylor), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: R2

    if(real_warning) call real_stop
    call check(s1)
    call ass(POWR)

    if(old) then
       CALL DAFUN('LOG ',S1%I,TEMP)
       CALL DACMU(TEMP,REAL(R2,kind=DP),TEMP)
       CALL DAFUN('EXP ',TEMP,TEMP)
       call dacop(temp,POWR%i)
    ELSE
       CALL NEWDAFUN('LOG ',S1%J,TEMPL)
       CALL NEWDACMU(TEMPL,REAL(R2,kind=DP),TEMPL)
       CALL NEWDAFUN('EXP ',TEMPL,POWR%J)
       !  CALL NEWDAFUN('EXP ',TEMPL,TEMPL)
       !  call newdacop(TEMPL,POWR%J)
    endif
  END FUNCTION POWR



  FUNCTION dmulsc( S1, sc )
    implicit none
    TYPE (taylor) dmulsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(dmulsc)

    if(old) then
       call dacmu(s1%i,sc,temp)
       call dacop(temp,dmulsc%i)
    else
       call newdacmu(s1%j,sc,dmulsc%j)
       !  call newdacmu(s1%j,sc,templ)
       ! call newdacop(templ,dmulsc%j)
    endif


  END FUNCTION dmulsc

  FUNCTION mulsc( S1, sc )
    implicit none
    TYPE (taylor) mulsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(mulsc)

    if(old) then
       call dacmu(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,mulsc%i)
    else
       call newdacmu(s1%j,REAL(sc,kind=DP),mulsc%j)
       ! call newdacmu(s1%j,REAL(sc,kind=DP),templ)
       ! call newdacop(templ,mulsc%j)
    endif



  END FUNCTION mulsc

  FUNCTION imulsc( S1, sc )
    implicit none
    TYPE (taylor) imulsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(imulsc)


    if(old) then
       call dacmu(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,imulsc%i)
    else
       call newdacmu(s1%j,REAL(sc,kind=DP),imulsc%j)
       !  call newdacmu(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,imulsc%j)
    endif


  END FUNCTION imulsc

  FUNCTION dscmul( sc,S1 )
    implicit none
    TYPE (taylor) dscmul
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(dscmul)

    if(old) then
       call dacmu(s1%i,sc,temp)
       call dacop(temp,dscmul%i)
    else
       call newdacmu(s1%j,sc,dscmul%j)
       !  call newdacmu(s1%j,sc,templ)
       !  call newdacop(templ,dscmul%j)
    endif


  END FUNCTION dscmul

  FUNCTION scmul( sc,S1 )
    implicit none
    TYPE (taylor) scmul
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(scmul)


    if(old) then
       call dacmu(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,scmul%i)
    else
       call newdacmu(s1%j,REAL(sc,kind=DP),scmul%j)
       !  call newdacmu(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,scmul%j)
    endif


  END FUNCTION scmul

  FUNCTION iscmul( sc,S1 )
    implicit none
    TYPE (taylor) iscmul
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(iscmul)

    if(old) then
       call dacmu(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,iscmul%i)
    else
       call newdacmu(s1%j,REAL(sc,kind=DP),iscmul%j)
       !  call newdacmu(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,iscmul%j)
    endif


  END FUNCTION iscmul

  FUNCTION div( S1, S2 )
    implicit none
    TYPE (taylor) div
    TYPE (taylor), INTENT (IN) :: S1, S2

    call check(s1)
    call check(s2)
    call ass(div)

    if(old) then
       call dadiv(s1%i,s2%i,temp)
       call dacop(temp,div%i)
    else
       ! call newdadiv(s1%j,s2%j,div%j)
       call newdadiv(s1%j,s2%j,templ)
       call newdacop(templ,div%j)
    endif



  END FUNCTION div

  FUNCTION dscdiv( sc,S1 )
    implicit none
    TYPE (taylor) dscdiv
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(dscdiv)

    if(old) then
       call dadic(s1%i,sc,temp)
       call dacop(temp,dscdiv%i)
    else
       call newdadic(s1%j,sc,dscdiv%j)
       !  call newdadic(s1%j,sc,templ)
       !  call newdacop(templ,dscdiv%j)
    endif


  END FUNCTION dscdiv

  FUNCTION scdiv( sc,S1 )
    implicit none
    TYPE (taylor) scdiv
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(scdiv)


    if(old) then
       call dadic(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,scdiv%i)
    else
       call newdadic(s1%j,REAL(sc,kind=DP),scdiv%j)
       ! call newdadic(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,scdiv%j)
    endif


  END FUNCTION scdiv

  FUNCTION iscdiv( sc,S1 )
    implicit none
    TYPE (taylor) iscdiv
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(iscdiv)

    if(old) then
       call dadic(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,iscdiv%i)
    else
       call newdadic(s1%j,REAL(sc,kind=DP),iscdiv%j)
       !  call newdadic(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,iscdiv%j)
    endif



  END FUNCTION iscdiv

  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (taylor) ddivsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(ddivsc)


    if(old) then
       call dacdi(s1%i,sc,temp)
       call dacop(temp,ddivsc%i)
    else
       call newdacdi(s1%j,sc,ddivsc%j)
       !  call newdacdi(s1%j,sc,templ)
       !  call newdacop(templ,ddivsc%j)
    endif

  END FUNCTION ddivsc

  FUNCTION divsc( S1, sc )
    implicit none
    TYPE (taylor) divsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(divsc)

    if(old) then
       call dacdi(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,divsc%i)
    else
       call newdacdi(s1%j,REAL(sc,kind=DP),divsc%j)
       !  call newdacdi(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,divsc%j)
    endif


  END FUNCTION divsc


  FUNCTION idivsc( S1, sc )
    implicit none
    TYPE (taylor) idivsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(idivsc)


    if(old) then
       call dacdi(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,idivsc%i)
    else
       call newdacdi(s1%j,REAL(sc,kind=DP),idivsc%j)
       !  call newdacdi(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,idivsc%j)
    endif


  END FUNCTION idivsc


  FUNCTION add( S1, S2 )
    implicit none
    TYPE (taylor) add
    TYPE (taylor), INTENT (IN) :: S1, S2

    call check(s1)
    call check(s2)
    call ass(add)


    if(old) then
       call daadd(s1%i,s2%i,add%i)
       !  call dacop(temp,add%i)
       !  call daadd(s1%i,s2%i,temp)
       !  call dacop(temp,add%i)
    else
       call newdaadd(s1%j,s2%j,add%j)
       !  call newdaadd(s1%j,s2%j,templ)
       !  call newdacop(templ,add%j)
    endif



  END FUNCTION add

  FUNCTION daddsc( S1, sc )
    implicit none
    TYPE (taylor) daddsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(daddsc)

    if(old) then
       call dacad(s1%i,sc,temp)
       call dacop(temp,daddsc%i)
    else
       call newdacad(s1%j,sc,daddsc%j)
       !  call newdacad(s1%j,sc,templ)
       !  call newdacop(templ,daddsc%j)
    endif


  END FUNCTION daddsc

  FUNCTION addsc( S1, sc )
    implicit none
    TYPE (taylor) addsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(addsc)


    if(old) then
       call dacad(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,addsc%i)
    else
       call newdacad(s1%j,REAL(sc,kind=DP),addsc%j)
       !  call newdacad(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,addsc%j)
    endif


  END FUNCTION addsc

  FUNCTION iaddsc( S1, sc )
    implicit none
    TYPE (taylor) iaddsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(iaddsc)

    if(old) then
       call dacad(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,iaddsc%i)
    else
       call newdacad(s1%j,REAL(sc,kind=DP),iaddsc%j)
       !  call newdacad(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,iaddsc%j)
    endif


  END FUNCTION iaddsc

  FUNCTION dscadd( sc,S1)
    implicit none
    TYPE (taylor) dscadd
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(dscadd)

    if(old) then
       call dacad(s1%i,sc,temp)
       call dacop(temp,dscadd%i)
    else
       call newdacad(s1%j,sc,dscadd%j)
       !  call newdacad(s1%j,sc,templ)
       !  call newdacop(templ,dscadd%j)
    endif



  END FUNCTION dscadd

  FUNCTION scadd( sc,S1)
    implicit none
    TYPE (taylor) scadd
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(scadd)

    if(old) then
       call dacad(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,scadd%i)
    else
       call newdacad(s1%j,REAL(sc,kind=DP),scadd%j)
       ! call newdacad(s1%j,REAL(sc,kind=DP),templ)
       ! call newdacop(templ,scadd%j)
    endif



  END FUNCTION scadd

  FUNCTION iscadd( sc,S1)
    implicit none
    TYPE (taylor) iscadd
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(iscadd)


    if(old) then
       call dacad(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,iscadd%i)
    else
       call newdacad(s1%j,REAL(sc,kind=DP),iscadd%j)
       ! call newdacad(s1%j,REAL(sc,kind=DP),templ)
       ! call newdacop(templ,iscadd%j)
    endif


  END FUNCTION iscadd

  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (taylor) subs
    TYPE (taylor), INTENT (IN) :: S1, S2

    call check(s1)
    call check(s2)
    call ass(subs)



    if(old) then
       call dasub(s1%i,s2%i,temp)
       call dacop(temp,subs%i)
    else
       call newdasub(s1%j,s2%j,subs%j)
       !call newdasub(s1%j,s2%j,templ)
       !call newdacop(templ,subs%j)
    endif


  END FUNCTION subs

  FUNCTION dsubsc( S1, sc )
    implicit none
    TYPE (taylor) dsubsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(dsubsc)

    if(old) then
       call dacsu(s1%i,sc,temp)
       call dacop(temp,dsubsc%i)
    else
       call newdacsu(s1%j,sc,dsubsc%j)
       !  call newdacsu(s1%j,sc,templ)
       !  call newdacop(templ,dsubsc%j)
    endif



  END FUNCTION dsubsc

  FUNCTION subsc( S1, sc )
    implicit none
    TYPE (taylor) subsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(subsc)

    if(old) then
       call dacsu(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,subsc%i)
    else
       call newdacsu(s1%j,REAL(sc,kind=DP),subsc%j)
       !  call newdacsu(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,subsc%j)
    endif


  END FUNCTION subsc

  FUNCTION isubsc( S1, sc )
    implicit none
    TYPE (taylor) isubsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(isubsc)

    if(old) then
       call dacsu(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,isubsc%i)
    else
       call newdacsu(s1%j,REAL(sc,kind=DP),isubsc%j)
       !  call newdacsu(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,isubsc%j)
    endif

  END FUNCTION isubsc

  FUNCTION dscsub( sc,S1)
    implicit none
    TYPE (taylor) dscsub
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc

    call check(s1)
    call ass(dscsub)

    if(old) then
       call dasuc(s1%i,sc,temp)
       call dacop(temp,dscsub%i)
    else
       call newdasuc(s1%j,sc,dscsub%j)
       !  call newdasuc(s1%j,sc,templ)
       !  call newdacop(templ,dscsub%j)
    endif


  END FUNCTION dscsub

  FUNCTION scsub( sc,S1)
    implicit none
    TYPE (taylor) scsub
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc

    if(real_warning) call real_stop
    call check(s1)
    call ass(scsub)

    if(old) then
       call dasuc(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,scsub%i)
    else
       call newdasuc(s1%j,REAL(sc,kind=DP),scsub%j)
       ! call newdasuc(s1%j,REAL(sc,kind=DP),templ)
       ! call newdacop(templ,scsub%j)
    endif


  END FUNCTION scsub

  FUNCTION iscsub( sc,S1)
    implicit none
    TYPE (taylor) iscsub
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc

    call check(s1)
    call ass(iscsub)

    if(old) then
       call dasuc(s1%i,REAL(sc,kind=DP),temp)
       call dacop(temp,iscsub%i)
    else
       call newdasuc(s1%j,REAL(sc,kind=DP),iscsub%j)
       !  call newdasuc(s1%j,REAL(sc,kind=DP),templ)
       !  call newdacop(templ,iscsub%j)
    endif



  END FUNCTION iscsub

  subroutine check(s1)
    implicit none
    TYPE (taylor) s1
    if(.not.checkass) return
    if(old) then
       if(s1%i==0) then
          w_p=0
          w_p%nc=1
          w_p=(/"Should not be here: Assign variables"/)
          w_p%fc='(1((1X,A72),/))'
          CALL WRITE_E(100)
       endif
    else
       IF (.NOT. ASSOCIATED(s1%j%r)) then
          w_p=0
          w_p%nc=1
          w_p=(/"Should not be here: Assign variables"/)
          w_p%fc='(1((1X,A72),/))'
          CALL WRITE_E(101)
       endif
    endif
  end subroutine check

  subroutine ASSIGN
    implicit none
    integer i,j
    do i=1,ndumt
       iassdoluser(i)=0
       iass0user(i)=0
    enddo
    if(old) then
       CALL ETALL1(DUMMY)
       call etall1(temp)
    else
       CALL allocnewda(DUMMYl)
       call allocnewda(templ)
    endif
    IF(OLDSCHEME) THEN
       if(old) then
          do i=1,ndumt
             if(ndumuser(i)>ndummax) then
                w_p=0
                w_p%nc=1
                write(w_p%c(1),'(a10,1x,i4,1x,a10)') "ndumuser(", int(i),  ")>ndummax "
                w_p%fc='(1((1X,A72),/))'
                CALL WRITE_E(111)
             endif
             do j=1,ndumuser(i)
                CALL ETALL1(DUMuser(i,j))
             enddo
          enddo
       else
          do i=1,ndumt
             do j=1,ndumuser(i)
                CALL allocnewda(DUMluser(i,j))
             enddo
          enddo
       endif
    ELSE
       CALL set_up_level
    ENDIF
  end subroutine ASSIGN

  subroutine DEASSIGN
    implicit none
    integer i,j
    do i=1,ndumt
       iassdoluser(i)=0
       iass0user(i)=0
    enddo
    if(old) then
       CALL DADAL1(DUMMY)
       call DADAL1(temp)
    else
       CALL KILLnewdaS(DUMMYl)
       call KILLnewdaS(templ)
    endif
    IF(FIRSTSCHEME) THEN
       if(old) then
          if(dummy/=0) then
             do i=1,ndumt
                do j=1,ndumuser(i)
                   CALL dadal1(DUMuser(i,j))
                enddo
             enddo
          endif
       else
          do i=1,ndumt
             do j=1,ndumuser(i)
                CALL KILLnewdaS(DUMluser(i,j))
             enddo
          enddo
       endif
    ELSE  ! NEW SCHEME
       do i=1,ndumt
          CALL kill_DALEVEL(scratchda(I))
       ENDDO
    ENDIF
  end subroutine DEASSIGN


  subroutine ass0(s1)
    implicit none
    TYPE (taylor) s1
    if(.not.no_ndum_check) iass0user(master)=iass0user(master)+1

    if(oldscheme) then
       iassdoluser(master)=iassdoluser(master)+1
       if(iassdoluser(master).eq.ndumuser(master)+1) iassdoluser(master)=1
       if(old) then
          s1%i=dumuser(master,iassdoluser(master))
       else
          s1%j=dumluser(master,iassdoluser(master))
       endif
    else
       if(iass0user(master)>scratchda(master)%n) then
          call INSERT_DA( scratchda(master) )
       ELSE
          scratchda(master)%PRESENT=>scratchda(master)%PRESENT%NEXT
       ENDIF
       if(old) then
          s1%i=scratchda(master)%PRESENT%T%i
       else
          s1%j=scratchda(master)%PRESENT%T%j
       endif

    endif

  end subroutine ASS0

  SUBROUTINE  ndum_warning_user
    implicit none
    integer ipause,II(0:1)

    if(oldscheme) then
       w_p=0
       w_p%nc=13
       w_p%fc='(13((1X,A72),/))'

       write(w_p%c(1),'(a18,1x,i4)') "Master variable = ", master
       write(w_p%c(2),'(a20,1x,i4)') "Number of Scratch = ", NDUMUSER(master)

       w_p%c(3)= " *****************************************************************"
       w_p%c(4)= " *  In User's new overloaded code                                *"
       w_p%c(5)= " *  line maybe be too big in overloaded code                     *"
       w_p%c(6)= " *  you can disconnect this check if you think things are fine   *"
       w_p%c(7)= " *  Set global variable no_ndum_check=.true.                     *"
       w_p%c(8)= " *****************************************************************"
       w_p%c(9)= "  "
       w_p%c(10)= " You may also increase the number "
       write(w_p%c(11),'(a31,1x,i4)') " of scratch variables at level ",master
       write(w_p%c(12),'(a4,1x,i4)') " to ",iass0user(master)
       w_p%c(13)= "    "
       CALL WRITE_E(0)
    else
       w_p=0
       w_p%nc=3
       w_p%fc='(3((1X,A72),/))'
       w_p%c(1)=  " *****************************************************************"
       w_p%c(2)=  " *  Should never be here in New Linked List Scheme               *"
       w_p%c(3)=  " *****************************************************************"
    endif
    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A72),/))'
    w_p%c(1)= " do you want a crash? "
    call write_e
    call read(ipause)
    ii(2000*ipause)=0

  end SUBROUTINE  ndum_warning_user


  !  These are new general TPSA-Routines

  SUBROUTINE  VAR000(S1,R1,I1)
    implicit none
    INTEGER,INTENT(IN)::I1
    real(dp),INTENT(IN)::R1
    type (taylor),INTENT(INOUT)::S1

    if(old) then
       if(s1%i==0) call crap1("VAR000  1" )  !call  etall1(s1%i)
       if(i1.ne.0) then
          CALL DAVAR(s1%i,R1,I1)
       else
          CALL DACON(s1%i,R1)
       endif
    else
       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("VAR000  2" ) !call newetall(s1%j,1)

       if(i1.ne.0) then
          CALL NEWDAVAR(s1%J,R1,I1)
       else
          CALL NEWDACON(s1%j,R1)
       endif
    endif

  END SUBROUTINE VAR000

  SUBROUTINE  VAR001(S1,R1,R2,I1)
    implicit none
    INTEGER,INTENT(IN)::I1
    real(dp),INTENT(IN)::R1,R2
    type (taylor),INTENT(INOUT)::S1

    if(old) then
       if(s1%i==0)  call crap1("VAR001  1" ) !etall1(s1%i)
       if(i1.ne.0) then
          CALL DAVAR(s1%i,zero,I1)
          CALL DACMU(s1%i,R2,TEMP)
          CALL DACAD(TEMP,R1,s1%i)
       else
          CALL DACON(s1%i,R1)
       endif
    else
       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("VAR001  2" )  !call newetall(s1%j,1)

       if(i1.ne.0) then
          CALL NEWDAVAR(s1%J,zero,I1)
          CALL NEWDACMU(s1%J,R2,TEMPL)
          CALL NEWDACAD(TEMPL,R1,s1%J)
       else
          CALL NEWDACON(s1%j,R1)
       endif
    endif

  END SUBROUTINE VAR001

  SUBROUTINE  shift000(S1,S2,s)
    implicit none
    INTEGER,INTENT(IN)::s
    type (taylor),INTENT(IN)::S1
    type (taylor),INTENT(inout)::S2

    if(old) then
       if(s2%i==0) call crap1("shift000  1" )  !call etall1(s2%i)
       CALL DAshift(s1%i,s2%i,s)
    else
       if(.NOT. ASSOCIATED(s2%j%r))call crap1("shift000  2" )   ! call newetall(s2%j,1)

       CALL NEWDAshift(s1%j,s2%j,s)
    endif

  END SUBROUTINE shift000


  SUBROUTINE  pek000(S1,J,R1)
    implicit none
    INTEGER,INTENT(IN),dimension(:)::j
    real(dp),INTENT(inOUT)::R1
    type (taylor),INTENT(IN)::S1
    if(old) then
       if(s1%i==0) call crap1("pek000  1" )  !call etall1(s1%i)
       CALL DApek(s1%i,j,r1)
    else
       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("pek000  2" ) ! newetall(s1%j,1)

       CALL newDApek(s1%j,j,r1)
    endif

  END SUBROUTINE pek000

  SUBROUTINE  pok000(S1,J,R1)
    implicit none
    INTEGER,INTENT(in),dimension(:)::j
    real(dp),INTENT(in)::R1
    type (taylor),INTENT(inout)::S1
    if(old) then
       if(s1%i==0) call crap1("pok000 1" )  ! call etall1(s1%i)
       CALL DApok(s1%i,j,r1)
    else
       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("pok000  2" )  ! call newetall(s1%j,1)

       CALL newDApok(s1%j,j,r1)
    endif

  END SUBROUTINE pok000


  SUBROUTINE  tran(S1,r1,R2)
    implicit none
    real(dp),INTENT(in)::R1
    real(dp),INTENT(inout)::R2
    type (taylor),INTENT(inout)::S1
    if(old) then
       if(s1%i==0) call crap1("tran  1" )  ! call etall1(s1%i)
       call daran(s1%i,r1,R2)
    else
       if(.NOT. ASSOCIATED(s1%j%r))call crap1("tran  2" )  !  call newetall(s1%j,1)

       call newdaran(s1%j,r1,R2)
    endif

  END SUBROUTINE tran



  SUBROUTINE  CFU000(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    real(dp) FUN
    EXTERNAL FUN
    if(old) then
       if(s1%i==0) call crap1("CFU000  1" )  !  call etall1(s1%i)
       CALL DACFU(s2%i,FUN,s1%i)
    else
       if(.NOT. ASSOCIATED(s1%j%r))call crap1("CFU000  2" )  !  call newetall(s1%j,1)
       CALL NEWDACFU(s2%J,FUN,s1%J)
    endif

  END SUBROUTINE CFU000

  SUBROUTINE  CFUR(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    complex(dp) FUN
    EXTERNAL FUN
    if(old) then
       if(s1%i==0) call crap1("CFUR  1" )  ! call etall1(s1%i)
       CALL DACFUR(s2%i,FUN,s1%i)
    else
       if(.NOT. ASSOCIATED(s1%j%r))call crap1("CFUR  2" )  !  call newetall(s1%j,1)
       CALL NEWDACFUR(s2%J,FUN,s1%J)
    endif

  END SUBROUTINE CFUR

  SUBROUTINE  CFUI(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    complex(dp) FUN
    EXTERNAL FUN
    if(old) then
       if(s1%i==0)call crap1("CFUI  1" ) ! call etall1(s1%i)
       CALL DACFUI(s2%i,FUN,s1%i)
    else
       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("CFUI  2" ) !call newetall(s1%j,1)
       CALL NEWDACFUI(s2%J,FUN,s1%J)
    endif

  END SUBROUTINE CFUI

  SUBROUTINE  tpsaeps(r1)
    implicit none
    real(dp),INTENT(INOUT)::r1
    if(old) then
       CALL DAeps(r1)
    else
       CALL newDAeps(r1)
    endif

  END SUBROUTINE tpsaeps


  SUBROUTINE  pri(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (TAYLOR),INTENT(IN)::S1

    if(old) then
       if(print77) then
          CALL DAPRI77(s1%i,MFILE)
       else
          CALL DAPRI(s1%i,MFILE)
       endif
    else
       if(newprint) then
          CALL newDAPRI(s1%j,MFILE)
       else
          if(print77) then
             CALL oldDAPRI77(s1%j,MFILE)
          else
             CALL oldDAPRI(s1%j,MFILE)
          endif
       endif
    endif

  END SUBROUTINE pri

  SUBROUTINE  REA(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (TAYLOR),INTENT(IN)::S1

    if(old) then
       if(s1%i==0)call crap1("REA  1" ) !  call etall1(s1%i)

       if(read77) then
          CALL DAREA77(s1%i,MFILE)
       else
          CALL DAREA(s1%i,MFILE)
       endif
    else
       if(.NOT. ASSOCIATED(s1%j%r))call crap1("REA  2" ) ! call newetall(s1%j,1)
       if(newread) then
          CALL newDAREA(s1%j,MFILE)
       else
          if(read77) then
             CALL oldDAREA77(s1%j,MFILE)
          else
             CALL oldDAREA(s1%j,MFILE)
          endif
       endif
    endif

  END SUBROUTINE REA


  function filter(j)
    implicit none
    real(dp) filter
    integer i
    integer,dimension(:)::j

    filter=one
    !do i=1,nd2+ndel
    do i=1,nd2par
       if(jfil(i)/=j(i)) filter=zero
    enddo

  end  function filter

  FUNCTION GETCHARnd2( S1, S2 )
    implicit none
    TYPE (taylor) GETCHARnd2,junk
    TYPE (taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer i,k
    ndel=0
    call check(s1)
    call ass(GETCHARnd2)

    call alloc(junk)
    resul = trim(ADJUSTR (s2))

    do i=1,lnv
       jfil(i)=0
    enddo

    nd2par= len(trim(ADJUSTR (s2)))

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),Jfil(I))
       if(i>nv) then
          if(Jfil(i)>0) then
             GETCHARnd2=zero
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2par+1,nv
       if(jfil(i)/=0) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72),/))'
          w_p%c(1)=" error in getchar for .para. "
          call write_e(0)
          stop
       endif
    enddo

    call cfu(s1,filter,junk)

    !DO I=1,ND2+ndel
    DO I=1,ND2par
       DO K=1,JFIL(I)
          JUNK=JUNK.K.I
       ENDDO
    ENDDO

    GETCHARnd2=junk

    call kill(junk)

  END FUNCTION GETCHARnd2

  FUNCTION GETintnd2( S1, S2 )
    implicit none
    TYPE (taylor) GETintnd2,junk
    TYPE (taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2(:)
    integer i,k,nt
    call check(s1)
    call ass(GETintnd2)

    call alloc(junk)

    do i=1,lnv
       jfil(i)=0
    enddo
    nd2par=size(s2)
    ndel=0

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2par
       Jfil(I)=s2(i)
       if(i>nv) then
          if(Jfil(i)>0) then
             GETintnd2=zero
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2par+1,nv
       if(jfil(i)/=0) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72),/))'
          w_p%c(1)=" error in GETintnd2 for .para. "
          call write_e(0)
          stop
       endif
    enddo

    call cfu(s1,filter,junk)

    !DO I=1,ND2+ndel
    DO I=1,ND2par
       DO K=1,JFIL(I)
          JUNK=JUNK.K.I
       ENDDO
    ENDDO

    GETintnd2=junk

    call kill(junk)

  END FUNCTION GETintnd2



  SUBROUTINE  null_0(S2,S1)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    integer, intent(in):: s1
    IF(S1==0) THEN
       NULLIFY(S2%N,S2%NV,S2%C,S2%J)
    ELSEIF(S1==-1) THEN
       DEALLOCATE(S2%N,S2%NV,S2%C,S2%J)
    ENDIF
  END SUBROUTINE null_0

  SUBROUTINE  ALLOC_U(S2,N,NV)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    integer, intent(in):: N,NV
    ALLOCATE(S2%N,S2%NV)
    if(N==0) then
       allocate(S2%C(1),S2%J(1,NV));S2%C(1)=zero;S2%J(:,:)=0;
    else
       allocate(S2%C(N),S2%J(N,NV))
    endif
    S2%N=N
    S2%NV=NV
  END SUBROUTINE ALLOC_U

  SUBROUTINE  FILL_R(S2,S1)  !new sagan
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    real (dp), intent(in):: s1
    INTEGER inoc,invc,ipoc,k,n,I,J(LNV)
    logical(lp) DOIT


    IF(ASSOCIATED(S2%N)) S2=-1
    S2=0
    CALL ALLOC_U(S2,1,nv)
    J=0
    DO N=1,S2%NV
       S2%J(1,N)=J(N)
    ENDDO
    S2%C(1)=S1

  END SUBROUTINE FILL_R

  SUBROUTINE  FILL_N(S2,S1)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    type (TAYLOR), intent(in):: s1
    INTEGER inoc,invc,ipoc,k,n,I,J(LNV)
    logical(lp) DOIT

    if(old) then
       if(s1%i==0)  call crap1("FILL_N 1")
    else
       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("FILL_N 2")
    endif


    IF(ASSOCIATED(S2%N)) S2=-1
    S2=0
    IF(OLD) THEN
       CALL dainf(S1%I,inoc,invc,ipoc,k,N)
       CALL ALLOC_U(S2,N,invc)
       do i=1,N
          k=ipoC+i-1
          CALL GET_C_J(k,S2%C(I),J)
          if(no==1) then
             j=0
             k=k-ipoc
             if(k/=0) j(k)=1
          endif
          DO N=1,S2%NV
             S2%J(i,N)=J(N)
          ENDDO
       ENDDO
    else
       N=0
       DO i=1,SIZE(S1%J%R)
          IF(PACKING) THEN
             IF(S1%J%YES(I)) N=N+1
          ELSE
             IF(S1%J%R(I)/=zero)  N=N+1
          ENDIF
       ENDDO
       CALL ALLOC_U(S2,N,nv)
       N=0
       DO i=1,SIZE(S1%J%R)
          DOIT=.FALSE.
          IF(PACKING) THEN
             IF(S1%J%YES(I)) DOIT=.TRUE.
          ELSE
             IF(S1%J%R(I)/=zero)  DOIT=.TRUE.
          ENDIF

          IF(DOIT) THEN

             N=N+1
             CALL newdancd(I,j)
             DO K=1,S2%NV
                S2%J(N,K)=J(K)
             ENDDO
             S2%C(N)=S1%J%R(I)
          ENDIF
       ENDDO

    ENDIF
  END SUBROUTINE FILL_N

  SUBROUTINE  REFILL_N(S1,S2)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(IN)::S2
    type (TAYLOR), intent(inOUT):: s1
    INTEGER IPAUSE,MYPAUSE,I,K,J(LNV)
    logical(lp) DOIT

    if(old) then
       if(s1%i==0)  call crap1("REFILL_N 1")
    else
       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("REFILL_N 2")
    endif


    S1=zero

    IF(.not.ASSOCIATED(S2%N)) THEN
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72),/))'
       w_p%c(1)=" ERROR IN REFILL_N: UNIVERSAL_TAYLOR DOES NOT EXIST"
       call write_e(123)
    ENDIF
    J=0
    DO I=1,S2%N
       DOIT=.TRUE.
       IF(S2%NV>NV) THEN
          K=NV
          DO WHILE(DOIT.AND.K<=S2%NV)
             IF(S2%J(I,K)/=0) DOIT=.FALSE.
             K=K+1
          ENDDO
       ENDIF

       IF(DOIT) THEN
          DO K=1,NV
             J(K)=S2%J(I,K)
          ENDDO
          CALL POK(S1,J,S2%C(I))
       ENDIF
    ENDDO

  END SUBROUTINE REFILL_N

  subroutine check_snake    !Checks equal signs only in defined in TPSA.f90
    implicit none
    select case (master)
    case(1)
       if(oldscheme) then
          if(iass0user(master)>ndumuser(master)) then
             call ndum_warning_user
          endif
          iass0user(master)=0
       else
          if(iass0user(master)>scratchda(master)%n.or.scratchda(master)%n>newscheme_max) then
             w_p=0
             w_p%nc=1
             w_p%fc='(1((1X,A72),/))'
             w_p%fi='(3((1X,i4)))'
             w_p%c(1)= "iass0user(master),scratchda(master)%n,newscheme_max"
             w_p=(/iass0user(master),scratchda(master)%n,newscheme_max/)
             call write_e
             call ndum_warning_user
          endif
          iass0user(master)=0
       endif

    end select
  end  subroutine check_snake

  subroutine crap1(STRING)
    implicit none
    CHARACTER(*) STRING

    w_p=0
    w_p%nc=2
    w_p%fc='((1X,A72,/),(1X,A72))'
    w_p%c(1)= "ERROR IN :"
    w_p%c(2)= STRING
    call write_e(3478)

  end subroutine crap1


  ! linked list

  SUBROUTINE Set_Up( L ) ! Sets up a layout: gives a unique negative index
    implicit none
    TYPE (dalevel) L
    call null_it(L)
    ALLOCATE(L%n);
    ALLOCATE(L%CLOSED);
    L%closed=.FALSE.
    L%N=0
  END SUBROUTINE Set_Up

  SUBROUTINE de_Set_Up( L ) ! deallocates layout content
    implicit none
    TYPE (dalevel) L
    deallocate(L%closed);
    deallocate(L%n);
  END SUBROUTINE de_Set_Up



  SUBROUTINE null_it( L ) ! Nullifies layout content
    implicit none
    TYPE (dalevel), intent(inout) :: L
    nullify(L%N )
    nullify(L%CLOSED )
    nullify(L%PRESENT )
    !
    nullify(L%END )
    nullify(L%START )
    nullify(L%START_GROUND )! STORE THE GROUNDED VALUE OF START DURING CIRCULAR SCANNING
    nullify(L%END_GROUND )! STORE THE GROUNDED VALUE OF END DURING CIRCULAR SCANNING
  END SUBROUTINE null_it

  SUBROUTINE LINE_L(L,doneit) ! makes into line temporarily
    implicit none
    TYPE (DALEVEL) L
    logical(lp) doneit
    doneit=.false.
    if(L%closed)  then
       if(associated(L%end%next)) then
          L%end%next=>L%start_ground
          doneit=.true.
       endif
       if(associated(L%start%previous)) then
          L%start%previous=>L%end_ground
       endif
    endif
  END SUBROUTINE LINE_L

  SUBROUTINE RING_L(L,doit) ! Brings back to ring if needed
    implicit none
    TYPE (DALEVEL) L
    logical(lp) doit
    if(L%closed.and.doit)  then
       if(.NOT.(associated(L%end%next))) then
          L%start_ground=>L%end%next      ! saving grounded pointer
          L%end%next=>L%start
       endif
       if(.NOT.(associated(L%start%previous))) then
          L%end_ground=>L%start%previous  ! saving grounded pointer
          L%start%previous=>L%end
       endif
    endif
  END SUBROUTINE RING_L

  SUBROUTINE APPEND_DA( L ) ! Standard append that clones everything
    implicit none
    TYPE (dascratch), POINTER :: Current
    TYPE (DALEVEL), TARGET,intent(inout):: L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    L%N=L%N+1
    nullify(current);ALLOCATE(current);

    call alloc_DA(current)

    if(L%N==1) current%next=> L%start
    Current % previous => L % end  ! point it to next fibre
    if(L%N>1)  THEN
       L % end % next => current      !
    ENDIF

    L % end => Current
    if(L%N==1) L%start=> Current
    L%PRESENT=>CURRENT    ! ALWAYS IF APPENDING
    CALL RING_L(L,doneit)
  END SUBROUTINE APPEND_DA

  SUBROUTINE INSERT_DA( L ) ! Standard append that clones everything
    implicit none
    logical(lp) :: doneitt=.true.
    TYPE (dascratch), POINTER :: Current
    TYPE (DALEVEL), TARGET,intent(inout):: L
    IF(L%N>1.AND.(.NOT.ASSOCIATED(L%PRESENT,L%END))) THEN

       L%N=L%N+1
       nullify(current);ALLOCATE(current);

       call alloc_DA(current)

       Current % previous => L % PRESENT   ! 2P -> 2
       Current % NEXT => L % PRESENT%NEXT  ! 2P -> 3
       L%PRESENT%NEXT=> CURRENT            ! 2  -> 2P
       Current % NEXT%PREVIOUS => CURRENT  ! 3  -> 2P
       L%PRESENT=>CURRENT                  ! 2P BECOMES 3

    ELSE

       CALL APPEND_DA( L )
       if(L%N==1) THEN
          L%CLOSED=.TRUE.
          CALL RING_L(L,doneitt)
       ENDIF

    ENDIF
  END SUBROUTINE INSERT_DA

  SUBROUTINE alloc_DA( c ) ! Does the full allocation of fibre and initialization of internal variables
    implicit none
    type(dascratch),pointer:: c
    ALLOCATE(C%T)
    CALL ALLOC(C%T)
    NULLIFY(C%NEXT)
    NULLIFY(C%PREVIOUS)

  end SUBROUTINE alloc_DA

  SUBROUTINE kill_DALEVEL( L )  ! Destroys a layout
    implicit none
    TYPE (DASCRATCH), POINTER :: Current
    TYPE (DALEVEL) L
    logical(lp) doneit
    CALL LINE_L(L,doneit)
    nullify(current)
    Current => L % end      ! end at the end
    DO WHILE (ASSOCIATED(L % end))
       L % end => Current % previous  ! update the end before disposing
       call dealloc_DASCRATCH(Current)
       Current => L % end     ! alias of last fibre again
       L%N=L%N-1
    END DO
    call de_set_up(L)
  END SUBROUTINE kill_DALEVEL

  SUBROUTINE dealloc_DASCRATCH( c ) ! destroys internal data  if it is not pointing (i.e. not a parent)
    implicit none
    type(DASCRATCH),pointer :: c
    IF(ASSOCIATED(C)) THEN
       CALL KILL(C%T)
       IF(ASSOCIATED(C%T)) DEALLOCATE(C%T)
       !       IF(ASSOCIATED(C%NEXT)) DEALLOCATE(C%NEXT)
       !       IF(ASSOCIATED(C%PREVIOUS)) DEALLOCATE(C%PREVIOUS)
       deallocate(c);
    ENDIF
  end SUBROUTINE dealloc_DASCRATCH

  SUBROUTINE set_up_level
    implicit none
    integer i
    do i=1,ndumt
       call set_up(scratchda(i))
       !    do j=1,n
       !      call INSERT_da(scratchda(i))
       !    enddo
       !    scratchda(i)%CLOSED=.TRUE.
       !    CALL RING_L(scratchda(i),.TRUE.)
    enddo

  end   SUBROUTINE set_up_level

  SUBROUTINE report_level
    implicit none
    integer i
    if((.not.oldscheme).and.associated(scratchda(1)%n)) then
       do i=1,ndumt
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72)))'
          write(w_p%c(1),'(a6,1x,i4,a5,1x,i4,1x,a7)') "Level ",i, " has ",scratchda(i)%n, "Taylors"
          call write_e
       enddo
    else
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72)))'
       w_p%c(1)="There is nothing to report on the levels"
       call write_e
    endif
  END   SUBROUTINE report_level

  SUBROUTINE real_stop
    implicit none
    integer i(1),j

    w_p=0
    w_p%nc=3
    w_p%c(1)=" You are using a kind(one) "
    w_p%c(2)=" set real_warning to false to permit this "
    w_p%c(3)=" write 1 to continue or -1 for a crash "
    w_p%fc='(3((1X,A60),/))'
    CALL WRITE_E
    call read(j)
    i(j)=0
    real_warning=.false.

  END   SUBROUTINE real_stop

END MODULE  tpsa
