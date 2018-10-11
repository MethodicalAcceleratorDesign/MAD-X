!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module tpsalie_analysis
  use tpsalie
  implicit none
  public
  private allocvecres,allocpbres,allocONELIE,allocdf,allocfd,allocnormal
  private normalMAP,MAPrevdf,revdfMAP,MAPdf,dfMAP,respb,pbres,resovec,vecreso
  private ONEEXPMAP,killnormal,KILLdf,KILLREVdf,killpbres,KILLvecres,killONELIE
  private DAPRINTvecres,DAPRINTpbres,DAPRINTdf,DAPRINTrevdf
  private df_map,fd_map,map_df,map_fd,one_map,map_one,MAPnormal
  private  resta,tares,alloctares,killtares
  private allocgen,KILLgen   !,ASSgen   ,ADDID,IDADD
  private EQUALgenMAP,EQUALMAPgen,pushgen
  !private rotf,rotsymp,equalgengen
  integer,private::NO,ND,ND2,NP,NDPT,NV
  logical(lp),private::old
  logical(lp) imaxflag
  !
  !
  !
  private NRESO

  INTERFACE init_tpsalie
     MODULE PROCEDURE init_map
     MODULE PROCEDURE init_tpsa
  END INTERFACE


  INTERFACE assignment (=)
     MODULE PROCEDURE normalMAP
     MODULE PROCEDURE MAPnormal
     MODULE PROCEDURE MAPrevdf
     MODULE PROCEDURE revdfMAP
     MODULE PROCEDURE MAPdf
     MODULE PROCEDURE dfMAP
     MODULE PROCEDURE respb
     MODULE PROCEDURE pbres
     MODULE PROCEDURE resta  !polynomial in resonance basis
     MODULE PROCEDURE tares
     MODULE PROCEDURE resovec
     MODULE PROCEDURE vecreso
     MODULE PROCEDURE MAPONEEXP
     MODULE PROCEDURE ONEEXPMAP
     !     MODULE PROCEDURE EQUALgengen !  not ready
     MODULE PROCEDURE EQUALMAPgen !  not ready
     MODULE PROCEDURE EQUALgenMAP  !  not ready
     !radiation
!     MODULE PROCEDURE beamrad
  end  INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE df_map
     MODULE PROCEDURE fd_map
     MODULE PROCEDURE map_df
     MODULE PROCEDURE map_fd
     MODULE PROCEDURE one_map
     MODULE PROCEDURE map_one
     MODULE PROCEDURE pushgen
  end  INTERFACE

  INTERFACE alloc
     MODULE PROCEDURE allocvecres
     MODULE PROCEDURE allocpbres
     MODULE PROCEDURE allocONELIE
     MODULE PROCEDURE allocdf
     MODULE PROCEDURE allocfd
     MODULE PROCEDURE allocnormal
     MODULE PROCEDURE allocgen
     MODULE PROCEDURE alloctares
     !radiation stochastic
!    MODULE PROCEDURE allocbeamenvelope
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE killONELIE
     MODULE PROCEDURE killnormal
     MODULE PROCEDURE KILLdf
     MODULE PROCEDURE KILLREVdf
     MODULE PROCEDURE killpbres
     MODULE PROCEDURE KILLvecres
     MODULE PROCEDURE KILLgen
     MODULE PROCEDURE killtares
     !radiation stochastic
!     MODULE PROCEDURE killbeamenvelope
  END INTERFACE



  INTERFACE DAprint
     MODULE PROCEDURE DAPRINTvecres
     MODULE PROCEDURE DAPRINTpbres
     MODULE PROCEDURE DAPRINTdf
     MODULE PROCEDURE DAPRINTrevdf
     MODULE PROCEDURE DAPRINTonelie
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE DAPRINTvecres
     MODULE PROCEDURE DAPRINTpbres
     MODULE PROCEDURE DAPRINTdf
     MODULE PROCEDURE DAPRINTrevdf
     MODULE PROCEDURE DAPRINTonelie
  END INTERFACE

  INTERFACE DAinput
     MODULE PROCEDURE DAREADvecres
     MODULE PROCEDURE DAREADpbres
     MODULE PROCEDURE DAREADdf
     MODULE PROCEDURE DAREADrevdf
     MODULE PROCEDURE DAREADonelie
  END INTERFACE

  INTERFACE read
     MODULE PROCEDURE DAREADvecres
     MODULE PROCEDURE DAREADpbres
     MODULE PROCEDURE DAREADdf
     MODULE PROCEDURE DAREADrevdf
     MODULE PROCEDURE DAREADonelie
  END INTERFACE


  !INTERFACE rotation
  !MODULE PROCEDURE rotf
  !MODULE PROCEDURE rotsymp
  !END INTERFACE

contains

  FUNCTION one_map(S1,S2)
    implicit none
    TYPE (damap) one_map ,junk
    TYPE (damap), INTENT (IN) :: S2
    TYPE (onelieexponent), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       one_map%v%i=0
      RETURN
    endif
    localmaster=master

    call checkdamap(s2)
    call assdamap(one_map)

    call alloc(junk)

    junk=1

    junk=s1

    one_map=junk*s2

    call kill(junk)
    master=localmaster

  END FUNCTION one_map

  FUNCTION map_one(S2,S1)
    implicit none
    TYPE (damap) map_one,junk
    TYPE (damap), INTENT (IN) :: S2
    TYPE (onelieexponent), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       map_one%v%i=0
      RETURN
    endif
    localmaster=master

    call checkdamap(s2)
    call assdamap(map_one)

    call alloc(junk)

    junk=1

    junk=s1

    map_one=s2*junk   ! 2002.12.6

    call kill(junk)
    master=localmaster

  END FUNCTION map_one


  FUNCTION df_map(S1,S2)
    implicit none
    TYPE (damap) df_map ,junk
    TYPE (damap), INTENT (IN) :: S2
    TYPE (dragtfinn), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       df_map%v%i=0
      RETURN
    endif
    localmaster=master

    call checkdamap(s2)
    call assdamap(df_map)

    call alloc(junk)

    junk=1

    junk=s1

    df_map=junk*s2

    call kill(junk)
    master=localmaster

  END FUNCTION df_map

  FUNCTION map_df(S2,S1)
    implicit none
    TYPE (damap) map_df,junk
    TYPE (damap), INTENT (IN) :: S2
    TYPE (dragtfinn), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       map_df%v%i=0
      RETURN
    endif
    localmaster=master

    call checkdamap(s2)
    call assdamap(map_df)

    call alloc(junk)

    junk=1

    junk=s1

    map_df=s2*junk    ! 2002.12.6

    call kill(junk)
    master=localmaster

  END FUNCTION map_df

  FUNCTION fd_map(S1,S2)
    implicit none
    TYPE (damap) fd_map ,junk
    TYPE (damap), INTENT (IN) :: S2
    TYPE (reversedragtfinn), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       fd_map%v%i=0
      RETURN
    endif
    localmaster=master

    call checkdamap(s2)
    call assdamap(fd_map)

    call alloc(junk)

    junk=1

    junk=s1

    fd_map=junk*s2

    call kill(junk)
    master=localmaster

  END FUNCTION fd_map

  FUNCTION map_fd(S2,S1)
    implicit none
    TYPE (damap) map_fd ,junk
    TYPE (damap), INTENT (IN) :: S2
    TYPE (reversedragtfinn), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       map_fd%v%i=0
      RETURN
    endif
    localmaster=master

    call checkdamap(s2)
    call assdamap(map_fd)

    call alloc(junk)

    junk=1

    junk=s1

    map_fd=s2*junk   ! 2002.12.6

    call kill(junk)
    master=localmaster

  END FUNCTION map_fd

  SUBROUTINE  go_to_fix_point(m,a1,a1i,nord)
    implicit none
    type (damap),INTENT(inOUT)::m,a1,a1i
    integer nord

    IF(.NOT.C_%STABLE_DA) RETURN

    call gofix(m%v%i,a1%v%i,a1i%v%i,nord)


  END SUBROUTINE go_to_fix_point


  SUBROUTINE  normalMAP(S2,S1)
    implicit none
    type (normalform),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    type (damap) JUNK
    real(dp) ZERO_(NDIM2)
    logical(lp) glo

    IF(.NOT.C_%STABLE_DA) RETURN
    glo=global_verbose
    call check_snake
    call alloc(junk)
    zero_(:)=0.0_dp
    !    if(old) then
    if(s2%normal%linear%V(1)%i==0)  call crap1("normalMAP 1") !call allocw(s2%normal%linear%V(1))  ! changed
    lielib_print(7)=0
    if(.not.S2%auto) then
       lielib_print(7)=-1
    endif
    ! global_verbose=.true.
    call setidpr(S2%plane)
    CALL INPUTRES(S2%M,S2%NRES)
    JUNK=S1
    S2%a%CONSTANT=JUNK
    JUNK=ZERO_
    if(s2%nord.le.0.or.s2%nord.gt.no) s2%nord=no
    if(s2%jtune.lt.0) s2%jtune=0
    CALL MAPNORMF(junk%V%i,s2%a%nonlinear%v%i,s2%a%linear%v%i,&
         & s2%A1%v%i,s2%normal%linear%v%i,s2%normal%nonlinear%v%i,s2%NORD,s2%jtune)
    !S2%a1=S2%a%CONSTANT
    !     S2%a%linear=S2%a%CONSTANT

    s2%normal%linear=(s2%normal%linear).sub.1
    S2%normal%pb=S2%normal%nonlinear
    S2%a%pb=S2%a%nonlinear

    CALL DHDJFLO(S2%NORMAL%NONLINEAR%V%i,S2%DHDJ%V%i)
    call GETTURA(s2%TUNE,S2%DAMPING)
    !    else
    !       if(.NOT.ASSOCIATED(s2%normal%linear%V(1)%j%r))  call crap1("normalMAP 2") ! call alloc(s2%normal%linear%V(1))
    !       if(S2%auto) then
    !          !     call setidpr(-100,S2%plane)
    !          call idprset(-101)
    !       else
    !          call setidpr(0,S2%plane)
    !       endif
    !       CALL INPUTRES(S2%M,S2%NRES)
    !       JUNK=S1
    !       !     S2%a%CONSTANT=JUNK
    !       JUNK=ZERO_
    !       if(s2%nord.le.0.or.s2%nord.gt.no) s2%nord=no
    !       if(s2%jtune.lt.0) s2%jtune=0
    !       CALL newMAPNORMF(junk%V%j,s2%a%nonlinear%v%j,s2%a%linear%v%j,&
    !            & s2%A1%v%j,s2%normal%linear%v%j,s2%normal%nonlinear%v%j,s2%NORD,s2%jtune)
    !       !     S2%a1=S2%a%CONSTANT
    !       !     S2%a%linear=S2%a%CONSTANT
    !       s2%normal%linear=(s2%normal%linear).sub.1
    !       S2%normal%pb=S2%normal%nonlinear
    !       S2%a%pb=S2%a%nonlinear
    !       CALL NEWDHDJFLO(S2%NORMAL%NONLINEAR%V%J,S2%DHDJ%V%J)
    !       call GETTURA(s2%TUNE,S2%DAMPING)
    !    endif

    s2%a_t=S2%a
    s2%a_t=s2%a1*s2%a_t
    S2%A_T=ZERO_
    CALL KILL(JUNK)
    global_verbose=glo


  END SUBROUTINE normalMAP

  SUBROUTINE  MAPnormal(S1,S2)
    implicit none
    type (normalform),INTENT(in)::S2
    type (damap),INTENT(INOUT)::S1
    type (damap) JUNK,id
    IF(.NOT.C_%STABLE_DA) RETURN

    call check_snake

    call alloc(JUNK,id)

    id=1
    id=texp(S2%normal%nonlinear,id)
    junk=S2%a1*S2%a

    s1=junk**(-1)*id*junk

    call kill(JUNK,id)

  END SUBROUTINE MAPnormal


  SUBROUTINE  ONEEXPMAP(S2,S1)
    implicit none
    type (ONELIEEXPONENT),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    type (damap) JUNK
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    call alloc(junk)
    !    if(old) then
    if(s2%VECTOR%V(1)%i==0) call crap1("ONEEXPMAP 1") ! call allocw(s2%VECTOR%V(1))
    JUNK=S1
    !     S2%CONSTANT=JUNK
    !     DO I=1,ND2
    !     ZERO_(I)=zero
    !     ENDDO
    !     JUNK=ZERO_
    IF(s2%eps==0.0_dp) THEN
       S2%EPS=1e3_dp*FULL_ABS(S1)
       call FLOFACG(JUNK%V%i,s2%VECTOR%V%i,s2%eps)
       S2%EPS=0.0_dp
    ELSE
       call FLOFACG(JUNK%V%i,s2%VECTOR%V%i,s2%eps)
    ENDIF
    S2%pb=S2%VECTOR
    CALL KILL(JUNK)
    !    else
    !       if(.NOT.ASSOCIATED(s2%VECTOR%V(1)%j%r)) THEN
    !          call crap1("ONEEXPMAP 2")  !call allocw(s2%VECTOR%V(1))
    !       ENDIF
    !       JUNK=S1
    !       !     S2%CONSTANT=JUNK
    !       !     DO I=1,ND2
    !       !     ZERO_(I)=zero
    !       !     ENDDO
    !       !     JUNK=ZERO_
    !       IF(s2%eps==ZERO) THEN
    !          S2%EPS=c_1d3*FULL_ABS(S1)
    !          call newFLOFACG(JUNK%V%j,s2%VECTOR%V%j,s2%eps)
    !          S2%EPS=ZERO
    !       ELSE
    !          call newFLOFACG(JUNK%V%j,s2%VECTOR%V%j,s2%eps)
    !       ENDIF
    !       S2%pb=S2%VECTOR
    !       CALL KILL(JUNK)
    !    endif
    !
  END SUBROUTINE ONEEXPMAP


  SUBROUTINE  MAPONEEXP(S1,S2)
    implicit none
    type (ONELIEEXPONENT),INTENT(IN)::S2
    type (damap),INTENT(inout)::S1
    type (damap) JUNK
    IF(.NOT.C_%STABLE_DA) RETURN
    call alloc(junk)
    call check_snake
    junk=1
    s1=texp(s2%vector,junk)
    !     s1=s2%constant
    call kill(junk)
  END SUBROUTINE MAPONEEXP




  SUBROUTINE  revdfMAP(S2,S1)
    implicit none
    type (reversedragtfinn),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    type (damap) JUNK
    real(dp) ZERO_(NDIM2)
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    call alloc(junk)
    ZERO_(:)=0.0_dp
    !    if(old) then
    if(s2%linear%V(1)%i==0)  call crap1("revdfMAP 1")  !call allocw(s2%linear%V(1))
    JUNK=S1
    S2%CONSTANT=JUNK
    JUNK=ZERO_
    junk=junk**(-1)
    call FLOFAC(JUNK%V%i,s2%linear%V%i,s2%nonlinear%V%i)
    s2%linear=s2%linear**(-1)
    call dacmud(s2%nonlinear%V%i,-1.0_dp,s2%nonlinear%V%i)
    !    S2%nonlinear=S2%nonlinear*S2%linear
    S2%pb=S2%nonlinear
    !    else
    !       if(.NOT.ASSOCIATED(s2%linear%V(1)%j%r))  call crap1("revdfMAP 2")  !call allocw(s2%linear%V(1))
    !       JUNK=S1
    !       S2%CONSTANT=JUNK
    !       JUNK=ZERO_
    !       junk=junk**(-1)
    !       call newFLOFAC(JUNK%V%j,s2%linear%V%j,s2%nonlinear%V%j)
    !       s2%linear=s2%linear**(-1)
    !       call newdacmud(s2%nonlinear%V%j,-one,s2%nonlinear%V%j)
    !       !    S2%nonlinear=S2%nonlinear*S2%linear
    !       S2%pb=S2%nonlinear
    !    endif
    CALL KILL(JUNK)

  END SUBROUTINE revdfMAP


  SUBROUTINE  dfMAP(S2,S1)
    implicit none
    type (dragtfinn),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    type (damap) JUNK
    real(dp) ZERO_(NDIM2)
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    call alloc(junk)
    ZERO_(:)=0.0_dp
    !    if(old) then
    if(s2%linear%V(1)%i==0)  call crap1("dfMAP 1")  !call allocw(s2%linear%V(1))
    JUNK=S1
    S2%CONSTANT=JUNK
    JUNK=ZERO_
    call FLOFAC(JUNK%V%i,s2%linear%V%i,s2%nonlinear%V%i)
    S2%pb=S2%nonlinear
    !    else
    !       if(.NOT.ASSOCIATED(s2%linear%V(1)%j%r))  call crap1("dfMAP 2")  !call allocw(s2%linear%V(1))
    !       JUNK=S1
    !       S2%CONSTANT=JUNK
    !       JUNK=ZERO_
    !       call newFLOFAC(JUNK%V%j,s2%linear%V%j,s2%nonlinear%V%j)
    !       S2%pb=S2%nonlinear
    !    endif
    CALL KILL(JUNK)
  END SUBROUTINE dfMAP

  SUBROUTINE  resovec(S1,S2)
    implicit none
    type (vecfield),INTENT(IN)::S2
    type (vecresonance),INTENT(inOUT)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    !    if(old) then
    if(S1%cos%v(1)%i==0) call crap1("resovec 1")  !call allocw(s1%cos%v(1))
    call ctorflo(s2%v%i,S1%cos%v%i,S1%sin%v%i)
    !    else
    !       if(.NOT.ASSOCIATED(S1%cos%v(1)%j%r)) call crap1("resovec 2")  !call allocw(s1%cos%v(1))
    !       call newctorflo(s2%v%j,S1%cos%v%j,S1%sin%v%j)
    !    endif
    s1%ifac=s2%ifac
    s1%cos%ifac=s2%ifac
    s1%sin%ifac=s2%ifac

  END SUBROUTINE resovec

  SUBROUTINE  vecreso(S2,S1)
    implicit none
    type (vecfield),INTENT(inout)::S2
    type (vecresonance),INTENT(in)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
    if(S2%v(1)%i==0) call crap1("vecreso 1")  !call allocw(s2%v(1))
    call RTOCflo(S1%cos%v%i,S1%sin%v%i,s2%v%i)
    !    else
    !       if(.not.associated(S2%v(1)%j%r)) call crap1("vecreso 2")  !call allocw(s2%v(1))
    !       call newRTOCflo(S1%cos%v%j,S1%sin%v%j,s2%v%j)
    !    endif
    s2%ifac=s1%ifac
  END SUBROUTINE vecreso

  SUBROUTINE  respb(S1,S2)
    implicit none
    type (pbfield),INTENT(IN)::S2
    type (pbresonance),INTENT(inOUT)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
    if(S1%cos%h%i==0) call crap1("respb 1")  !call allocw(s1%cos%h)
    call ctor(s2%h%i,S1%cos%h%i,S1%sin%h%i)
    !    else
    !       if(.not.associated(S1%cos%h%j%r)) call crap1("respb 1")
    !       call newctor(s2%h%j,S1%cos%h%j,S1%sin%h%j)
    !    endif
    s1%ifac=s2%ifac
    s1%cos%ifac=s2%ifac
    s1%sin%ifac=s2%ifac


  END SUBROUTINE respb

  SUBROUTINE  pbres(S2,S1)
    implicit none
    type (pbfield),INTENT(inOUT)::S2
    type (pbresonance),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
    if(s2%h%i==0) call crap1("respb 1")
    call RTOC(S1%cos%h%i,S1%sin%h%i,s2%h%i)
    !    else
    !       if(.not.associated(s2%h%j%r)) call crap1("respb 2")  !call allocw(s2%h)
    !       call newRTOC(S1%cos%h%j,S1%sin%h%j,s2%h%j)
    !    endif
    s2%ifac=s1%ifac
  END SUBROUTINE pbres

  SUBROUTINE  resta(S1,S2)
    implicit none
    type (taylor),INTENT(IN)::S2
    type (taylorresonance),INTENT(inOUT)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
    if(S1%cos%i==0)  call crap1("resta 1")  !call allocw(s1%cos)
    call ctor(s2%i,S1%cos%i,S1%sin%i)
    !    else
    !       if(.not.associated(S1%cos%j%r))  call crap1("resta 2")  !call allocw(s1%cos)
    !       call newctor(s2%j,S1%cos%j,S1%sin%j)
    !    endif
  END SUBROUTINE resta

  SUBROUTINE  tares(S2,S1)
    implicit none
    type (taylor),INTENT(inOUT)::S2
    type (taylorresonance),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    ! if(old) then
    if(s2%i==0)  call crap1("tares 1")  !call allocw(s2)
    call RTOC(S1%cos%i,S1%sin%i,s2%i)
    !   else
    !      if(.not.associated(s2%j%r))  call crap1("tares 2")  !call allocw(s2)
    !      call newRTOC(S1%cos%j,S1%sin%j,s2%j)
    !   endif
  END SUBROUTINE tares

  SUBROUTINE  MAPdf(S1,S2)
    implicit none
    type (dragtfinn),INTENT(IN)::S2
    type (damap),INTENT(inOUT)::S1
    TYPE (DAMAP) ID
    integer no1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    no1=no
    call alloc(ID)
    ! if(old) then
    if(s1%V(1)%i==0) call crap1("MAPdf 1")  !call allocw(s1%V(1))
    !    else
    !       if(.not.associated(s1%V(1)%j%r)) call crap1("MAPdf 2")  !call allocw(s1%V(1))
    !    endif
    ID=1
    ID=texpdf( S2%NONLINEAR, ID,2,NO1,1.0_dp,1 )
    S1=ID*S2%LINEAR
    S1=S2%CONSTANT
    CALL KILL(ID)
  END SUBROUTINE MAPdf

  SUBROUTINE  MAPrevdf(S1,S2)
    implicit none
    type (reversedragtfinn),INTENT(IN)::S2
    type (damap),INTENT(inOUT)::S1
    TYPE (DAMAP) ID
    integer no1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    no1=no
    call alloc(ID)
    ! if(old) then
    if(s1%V(1)%i==0) call crap1("MAPrevdf 1")  !call allocw(s1%V(1))
    !    else
    !       if(.NOT.ASSOCIATED(s1%V(1)%j%r)) call crap1("MAPrevdf 2")  !call allocw(s1%V(1))
    !    endif
    ID=1
    ID=texpdf( S2%NONLINEAR, ID,2,NO1,1.0_dp,-1 )
    S1=S2%LINEAR*ID
    S1=S2%CONSTANT
    CALL KILL(ID)
  END SUBROUTINE MAPrevdf

  SUBROUTINE  allocgen(S1)
    implicit none
    type (genfield),INTENT(INOUT)::S1
    integer i,j
    call alloc(s1%h)

    s1%eps=1e-6_dp
    s1%ifac=1
    s1%linear_in=.false.
    s1%imax=1000
    s1%constant(:)=0.0_dp
    imaxflag=.false.
    s1%no_cut=no+1
    do i=1,nd
       do j=1,nd
          call alloc(s1%d(i,j))
       enddo
    enddo
    do i=1,nd2
       call alloc(s1%linear%v(i))
       call alloc(s1%lineart%v(i))
       call alloc(s1%m%v(i))
       call alloc(s1%mt%v(i))
    enddo



  END SUBROUTINE allocgen

  SUBROUTINE  KILLgen(S1)
    implicit none
    type (genfield),INTENT(INOUT)::S1
    integer i,j
    s1%ifac=0
    s1%imax=0
    s1%no_cut=0
    do i=1,nd
       do j=1,nd
          call kill(s1%d(i,j))
       enddo
    enddo


    do i=1,nd2
       call kill(s1%lineart%v(i))
       call kill(s1%linear%v(i))
       call kill(s1%mt%v(i))
       call kill(s1%m%v(i))
    enddo
    ! if(old) then
    call DADAL1(s1%h%i)
    !    else
    !       call newDADAL(s1%h%j,1)
    !    endif
  END SUBROUTINE KILLgen

  SUBROUTINE  EQUALgenMAP(S2,S1)
    implicit none
    type (genfield),INTENT(inOUT)::S2
    type (damap),INTENT(IN)::S1
    type (damap)  w
    type(taylor) t
    type(onelieexponent) one_
    real(dp) zero_(ndim2)
    integer i,j,jn(lnv),k

    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    call alloc(w)
    call alloc(t)
    call alloc(one_)

    s2%constant=s1
    do i=1,ndim2
       zero_(i)=0.0_dp
    enddo
    do i=1,lnv
       jn(i)=0
    enddo
    do i=1,nd
       jn(2*i)=1
    enddo
    w=s1
    w=zero_
    s2%linear=w.sub.1
    !
    if(s2%linear_in) then
       s2%linear=1
    endif
    !
    w=w*s2%linear**(-1)
    do i=1,nd2
       call maketree(s2%linear%v(i),s2%lineart%v(i))
    enddo

    if(s2%ifac>1) then
       one_=w
       do k=1,nd2
          one_%vector%v(k)=one_%vector%v(k)/s2%ifac
       enddo
       w=one_
    endif


    ! if(old) then
    if(s2%h%i==0) call crap1("EQUALgenMAP 1")  ! call etall1(s2%h%i)
    if(s2%m%v(1)%i==0) call crap1("EQUALgenMAP 2")  !call etall(s2%m%v%i,nd2)
    call etpin(w%v%i,s2%m%v%i,jn)
    call intd(s2%m%v%i,s2%h%i,1.0_dp)
    !   else
    !       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALgenMAP 3")  !call newetall(s2%h%j,1)
    !       if(.NOT. ASSOCIATED(s2%m%v(1)%J%r)) call crap1("EQUALgenMAP 4")  !call newetall(s2%m%v%J,nd2)
    !       call newetpin(w%v%j,s2%m%v%j,jn)
    !       call newintd(s2%m%v%j,s2%h%j,one)
    !    endif
    ! Truncating to order no-1: not necessary completly but for self-conssitancy
    s2%m=s2%m.cut.s2%no_cut
    do i=1,nd
       do j=1,nd
          t=  (s2%m%v(2*i)).d.(2*j)
          call maketree(t,s2%d(i,j))
       enddo
    enddo

    do i=1,nd2
       call maketree(s2%m%v(i),s2%mt%v(i))
    enddo


    call kill(one_)
    call kill(t)
    call kill(w)

  END SUBROUTINE EQUALgenMAP

  SUBROUTINE  EQUALMAPgen(S1,S2)
    implicit none
    type (genfield),INTENT(in)::S2
    type (damap),INTENT(inout)::S1
    type (damap)  w
    integer i,jn(lnv)
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    call alloc(w)

    do i=1,lnv
       jn(i)=0
    enddo
    do i=1,nd
       jn(2*i)=1
    enddo


    ! if(old) then
    if(s1%v(1)%i==0)call crap1("EQUALMAPgen 1")  ! call etall(s1%v%i,nd2)
    call difd(s2%h%i,w%v%i,1.0_dp)
    call etpin(w%v%i,s1%v%i,jn)
    !   else
    !      if(.NOT. ASSOCIATED(s1%v(1)%J%r))call crap1("EQUALMAPgen 2")  ! call newetall(s1%v%J,nd2)
    !      call newdifd(s2%h%j,w%v%j,one)
    !      call newetpin(w%v%j,s1%v%j,jn)
    !   endif

    s1=s1*s2%linear

    s1=s2%constant
    call kill(w)
  END SUBROUTINE EQUALMAPgen

  !  SUBROUTINE  EQUALgengen(S2,S1)
  !    implicit none
  !    type (genfield),INTENT(inOUT)::S2
  !    type (genfield),INTENT(IN)::S1
  !    call check_snake
  !    if(old) then
  !       if(s2%h%i==0) call crap1("EQUALgengen 1")  ! call etall1(s2%h%i)
  !       CALL dacop(S1%h%i,S2%h%i)
  !    else
  !       if(.NOT. ASSOCIATED(s2%h%J%r)) call crap1("EQUALgengen 2")  ! call newetall(s2%h%J,1)
  !       CALL NEWdacop(S1%h%J,S2%h%J)
  !    endif
  !    s2%ifac=s1%ifac
  !  END SUBROUTINE EQUALgengen


  FUNCTION pushgen( S1, S2 )
    implicit none
    ! integer ipause, mypause
    TYPE (genfield), INTENT (IN) :: S1
    real(dp), intent(in),dimension(:)::s2
    real(dp) pushgen(lnv),s2t(nv)
    real(dp) junk(lnv),junk2(lnv),et(ndim),de(ndim),mt(ndim,ndim),MI(ndim,ndim)
    real(dp) eb,e
    logical(lp) more
    integer i,j,k,imax,ier,ifac

    IF(.NOT.C_%STABLE_DA) then
       pushgen=0
      RETURN
    endif

    if(.not.imaxflag) then
       junk(:)=0
       do i=1,nd2
          junk(i)= push1pol(s1%lineart%v(i),s2)
       enddo
       do i=nd2+1,nv   ! added
          junk(i)=s2(i)
       enddo

       imax=s1%imax

       do ifac=1,s1%ifac   ! ifac

          more=.true.

          eb=1.1111e5_dp

          do i=1,nv
             S2t(i)= junk(i)
          enddo


          do i=1,imax

             de(:)=0.0_dp



             do k=1,nd
                de(k)=s2t(2*k)-push1pol(s1%mt%v(2*k),junk)
             enddo

             do j=1,nd
                do k=1,nd
                   mt(j,k)=push1pol(s1%d(j,k),junk)
                enddo
             enddo
             call matinv(mt,mI,nd,ndim,ier)
             et(:)=0.0_dp
             do j=1,nd
                do k=1,nd
                   et(j)=mI(j,k)*de(k)+et(j)
                enddo
             enddo

             e=0.0_dp
             do k=1,nd
                junk(2*k)=junk(2*k)+et(k)
                e=abs(et(k))+e
             enddo

             if(more) then
                if(e>S1%eps) then
                   eb=e
                else
                   eb=e
                   more=.false.
                endif
             else
                if(e>=eb) then
                   goto 100
                else
                   eb=e
                endif
             endif

          enddo

          !   write(6,*) " Imax reached in generating function "
          imaxflag=.true.
          !   ipause=mypause(472)
          goto 1000
100       continue
          do k=1,nd
             junk2(2*k-1)=push1pol(s1%mt%v(2*k-1),junk)
          enddo

          do k=1,nd
             junk(2*k-1)=junk2(2*k-1)
          enddo
       enddo ! ifac

       do k=1,nV
          pushgen(K)=junk(K)
       enddo


    endif

1000 continue

  END FUNCTION pushgen

  SUBROUTINE  allocpbres(S2)
    implicit none
    type (pbresonance),INTENT(INOUT)::S2
    call alloc(s2%cos)
    call alloc(s2%sin)
    s2%ifac=0
  END SUBROUTINE allocpbres

  SUBROUTINE  alloctares(S2)
    implicit none
    type (taylorresonance),INTENT(INOUT)::S2
    call alloc(s2%cos)
    call alloc(s2%sin)
  END SUBROUTINE alloctares



  SUBROUTINE  allocvecres(S2)
    implicit none
    type (vecresonance),INTENT(INOUT)::S2
    !     call alloc(s2%cos)
    !     call alloc(s2%sin)
    call allocTPSA(s2%cos)
    call allocTPSA(s2%sin)
    s2%ifac=0
  END SUBROUTINE allocvecres

  SUBROUTINE  allocdf(S2)
    implicit none
    type (dragtfinn),INTENT(INOUT)::S2
    !     call alloc(s2%linear)
    !     call alloc(s2%nonlinear)
    !     call alloc(s2%pb)
    call allocTPSA(s2%linear)
    call allocTPSA(s2%nonlinear)
    call allocTPSA(s2%pb)
    s2%nonlinear%ifac=1
    s2%pb%ifac=1
    s2%constant(:)=0.0_dp
  END SUBROUTINE allocdf

  SUBROUTINE  allocONELIE(S2)
    implicit none
    type (ONELIEEXPONENT)  S2
    call allocTPSA(s2%VECTOR)
    call allocTPSA(s2%pb)
    S2%EPS=0.0_dp
  END SUBROUTINE allocONELIE


  SUBROUTINE  allocnormal(S2)
    implicit none
    type (normalform),INTENT(INOUT)::S2
    integer i,j
    call allocTPSA(s2%a_t)
    call allocTPSA(s2%a%Linear)
    call allocTPSA(s2%a%nonlinear)
    call allocTPSA(s2%a%pb)
    call allocTPSA(s2%normal%Linear)
    call allocTPSA(s2%normal%nonlinear)
    call allocTPSA(s2%normal%pb)
    s2%normal%constant(:)=0.0_dp
    s2%a%constant(:)=0.0_dp
    s2%NORMAL%nonlinear%ifac=1
    s2%NORMAL%PB%ifac=1
    s2%A%nonlinear%ifac=-1
    s2%A%PB%ifac=-1
    s2%nord=no
    s2%auto=.true.
    s2%nres=0
    s2%jtune=0
    do i=1,ndim
       s2%tune(i)=0.0_dp
       s2%DAMPING(i)=0.0_dp
       s2%PLANE(i)=2*i-1
       do j=1,nreso
          s2%M(i,j)=0
       enddo
    enddo
    call allocTPSA(s2%a1)
    call allocTPSA(s2%DHDJ)
  END SUBROUTINE allocnormal

  SUBROUTINE  allocfd(S2)
    implicit none
    type (reversedragtfinn),INTENT(INOUT)::S2
    call allocTPSA(s2%linear)
    call allocTPSA(s2%nonlinear)
    call allocTPSA(s2%pb)
    s2%nonlinear%ifac=-1
    s2%pb%ifac=-1
    s2%constant(:)=0.0_dp
  END SUBROUTINE allocfd

  SUBROUTINE  killnormal(S2)
    implicit none
    type (normalform),INTENT(INOUT)::S2
    call killTPSA(s2%a_t)
    call killTPSA(s2%a%Linear)
    call killTPSA(s2%a%nonlinear)
    call killTPSA(s2%a%pb)
    call killTPSA(s2%normal%Linear)
    call killTPSA(s2%normal%nonlinear)
    call killTPSA(s2%normal%pb)
    call killTPSA(s2%a1)
    call killTPSA(s2%DHDJ)
  END SUBROUTINE killnormal

  SUBROUTINE  killONELIE(S2)
    implicit none
    type (ONELIEEXPONENT)  S2
    call killTPSA(s2%VECTOR)
    call killTPSA(s2%pb)
  END SUBROUTINE killONELIE

  SUBROUTINE  killpbres(S2)
    implicit none
    type (pbresonance),INTENT(INOUT)::S2
    !     call kill(s2%cos)
    !     call kill(s2%sin)
    call killTPSA(s2%cos)
    call killTPSA(s2%sin)
  END SUBROUTINE killpbres

  SUBROUTINE  killtares(S2)
    implicit none
    type (taylorresonance),INTENT(INOUT)::S2
    !     call kill(s2%cos)
    !     call kill(s2%sin)
    call killTPSA(s2%cos)
    call killTPSA(s2%sin)
  END SUBROUTINE killtares

  SUBROUTINE  killvecres(S2)
    implicit none
    type (vecresonance),INTENT(INOUT)::S2
    !     call kill(s2%cos)
    !     call kill(s2%sin)
    call killTPSA(s2%cos)
    call killTPSA(s2%sin)
  END SUBROUTINE killvecres


  SUBROUTINE  killdf(S2)
    implicit none
    type (dragtfinn),INTENT(INOUT)::S2
    !     call kill(s2%linear)
    !     call kill(s2%nonlinear)
    !     call kill(s2%pb)
    call killTPSA(s2%linear)
    call killTPSA(s2%nonlinear)
    call killTPSA(s2%pb)
    s2%nonlinear%ifac=0
    s2%pb%ifac=0
  END SUBROUTINE killdf

  SUBROUTINE  killREVdf(S2)
    implicit none
    type (reversedragtfinn),INTENT(INOUT)::S2
    !     call kill(s2%linear)
    !     call kill(s2%nonlinear)
    !     call kill(s2%pb)
    call killTPSA(s2%linear)
    call killTPSA(s2%nonlinear)
    call killTPSA(s2%pb)
    s2%nonlinear%ifac=0
    s2%pb%ifac=0
  END SUBROUTINE killREVdf



  SUBROUTINE  DAPRINTonelie(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (ONELIEEXPONENT),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC
    integer mfi
     mfi=6
     if(present(mfile)) mfi=mfile

    write(mfi,*) s1%eps, " Convergence Test Number"
    !do i=1,nd2
    !write(mfile,*) s1%constant(i)
    !enddo

    CALL DAPRINT(s1%VECTOR,MFILE,PREC)
    CALL DAPRINT(s1%pb,MFILE,PREC)

  END SUBROUTINE DAPRINTonelie

  SUBROUTINE  DAreadonelie(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (ONELIEEXPONENT),INTENT(inout)::S1

    read(mfile,*) s1%eps
    !do i=1,nd2
    !     read(mfile,*) s1%constant(i)
    !     enddo

    CALL DAinput(s1%VECTOR,MFILE)
    CALL DAinput(s1%pb,MFILE)

  END SUBROUTINE DAreadonelie


  SUBROUTINE  DAPRINTvecres(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (vecresonance),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC


    CALL DAPRINT(s1%cos,MFILE,PREC)
    CALL DAPRINT(s1%sin,MFILE,PREC)

  END SUBROUTINE DAPRINTvecres

  SUBROUTINE  DAreadvecres(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (vecresonance),INTENT(inout)::S1
    CALL dainput(s1%cos,MFILE)
    CALL dainput(s1%sin,MFILE)

  END SUBROUTINE DAreadvecres

  SUBROUTINE  DAPRINTdf(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (dragtfinn),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC
    integer i,mfi   
     mfi=6
     if(present(mfile)) mfi=mfile

    do i=1,nd2
       write(mfi,*) s1%constant(i)
    enddo
    CALL DAPRINT(s1%linear,MFILE,PREC)
    CALL DAPRINT(s1%nonlinear,MFILE,PREC)

  END SUBROUTINE DAPRINTdf

  SUBROUTINE  dareaddf(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (dragtfinn),INTENT(inout)::S1
    integer i
    do i=1,nd2
       read(mfile,*) s1%constant(i)
    enddo
    call dainput(s1%linear,MFILE)
    CALL dainput(s1%nonlinear,MFILE)

  END SUBROUTINE dareaddf

  SUBROUTINE  DAPRINTrevdf(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (reversedragtfinn),INTENT(inout)::S1
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC
    integer i,mfi
     mfi=6
     if(present(mfile)) mfi=mfile
    do i=1,nd2
       write(mfi,*) s1%constant(i)
    enddo
    CALL DAPRINT(s1%linear,MFILE,PREC)
    CALL DAPRINT(s1%nonlinear,MFILE,PREC)

  END SUBROUTINE DAPRINTrevdf

  SUBROUTINE  DAreadrevdf(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (reversedragtfinn),INTENT(inout)::S1
    integer i
    do i=1,nd2
       read(mfile,*) s1%constant(i)
    enddo
    CALL dainput(s1%linear,MFILE)
    CALL dainput(s1%nonlinear,MFILE)

  END SUBROUTINE DAreadrevdf

  SUBROUTINE  DAPRINTpbres(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (pbresonance),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC
    CALL DAPRINT(s1%cos,MFILE,PREC)
    CALL DAPRINT(s1%sin,MFILE,PREC)

  END SUBROUTINE DAPRINTpbres

  SUBROUTINE  DAreadpbres(S1,MFILE)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    type (pbresonance),INTENT(inout)::S1
    CALL dainput(s1%cos,MFILE)
    CALL dainput(s1%sin,MFILE)

  END SUBROUTINE DAreadpbres

  subroutine init_map(NO1,ND1,NP1,NDPT1,log1)
    implicit none
    integer NO1,ND1,NP1,NDPT1,i
    LOGICAL(lp) log1,present_tpsa

    nb_=0
    !    if(first_time) then
    !       first_time=.false.
    !       w_p=0
    !       w_p%nc=1
    !       w_p=(/" Welcome to TPSA Overloaded"/)
    !       w_p%fc='(1((1X,A72),/))'
    !       ! call ! WRITE_I
    !    else
    !       call DATERMINATE
    !       call kill(varf1)
    !       call kill(varf2)
    !    endif

    !    if(.not.first_time) then
    present_tpsa=lingyun_yang
    if(last_tpsa==1) then
       lingyun_yang=.true.
       call DATERMINATE
       call kill(varf1)
       call kill(varf2)
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       call DATERMINATE
       call kill(varf1)
       call kill(varf2)
    endif
    lingyun_yang=present_tpsa
    !    endif


    master=0  !  master=1   2002.12.25

    call RESET_APERTURE_FLAG
    old=log1
    NO=NO1
    ND=ND1
    ND2=ND1*2
    NP=NP1
    NDPT=NDPT1
    NV=ND2+NP
    newprint=.false.
    ! if(old) then
    call LIEINIT(NO1,NV,ND1,NDPT1)   !,0
 !   w_p=0
 !   w_p%nc=1
 !   w_p=(/" Berz's Package  "/)
 !   w_p%fc='(1((1X,A72),/))'
    !       ! call ! WRITE_I
    !    else
    !       if(no1>3) then
    !          w_p=0
    !          w_p%nc=1
    !          w_p=(/" No1 is too big: run old=.true."/)
    !          w_p%fc='(1((1X,A72)))'
    !          ! call !write_e(-1)
    !       endif
    !       w_p=0
    !       w_p%nc=1
    !       w_p=(/" Etienne's Experimental Package  "/)
    !       w_p%fc='(1((1X,A72)))'
    !       ! call ! WRITE_I
    !       call newLIEINIT(NO1,NV,ND1,NDPT1,0)
    !    endif
    !
    call set_in_tpsa( NO,ND,ND2,NP,NDPT,NV,old)
    call set_in_tpsalie( NO,ND,ND2,NP,NDPT,NV,old)



!    w_p=0
!    w_p%nc=1
!    w_p=(/"          NO          ND         ND2          NP        NDPT          NV"/)
!    w_p%fc='(1((1X,A72)))'
!    w_p%fi='(1x,6(6x,i6))'
!    w_p=(/NO,ND,ND2,NP,NDPT,NV/)
    ! call ! WRITE_I
    CALL ASSIGN
    !    CALL ASSIGNMAP
    call alloc(varf1)
    call alloc(varf2)
    npara_fpp=nd2
    NPARA_original=npara_fpp



  end subroutine init_map



  subroutine KILL_fpp()
    implicit none
    logical present_tpsa
    present_tpsa=lingyun_yang
    if(last_tpsa==1) then
       lingyun_yang=.true.
       call kill(varf1)
       call kill(varf2)
       call DATERMINATE  ! IN THIS MODULE
       CALL dealloc_all    ! IN DABNEW
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       call kill(varf1)
       call kill(varf2)
       call DATERMINATE  ! IN THIS MODULE
       CALL dealloc_all    ! IN DABNEW
    endif
    lingyun_yang=present_tpsa
  END  subroutine KILL_fpp


  subroutine init_tpsa(NO1,NP1,log1)
    implicit none
    integer NO1,ND1,NP1,NDPT1,i
    LOGICAL(lp) log1,present_tpsa

    !    if(first_time) then
    !       first_time=.false.
    !       w_p=0
    !       w_p%nc=1
    !       w_p=(/" Welcome to TPSA Overloaded"/)
    !       w_p%fc='(1((1X,A72),/))'
    !       ! call ! WRITE_I
    !    else
    !       call DATERMINATE
    !       call kill(varf1)
    !       call kill(varf2)
    !    endif

    !    if(.not.first_time) then
    nb_=0
    present_tpsa=lingyun_yang
    if(last_tpsa==1) then
       lingyun_yang=.true.
       call DATERMINATE
       call kill(varf1)
       call kill(varf2)
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       call DATERMINATE
       call kill(varf1)
       call kill(varf2)
    endif
    lingyun_yang=present_tpsa
    !    endif

    master=0  !  master=1   2002.12.25

    old=log1
    nd1=0
    ndpt1=0
    NO=NO1
    ND=ND1
    ND2=ND1*2
    NP=NP1
    NDPT=NDPT1
    NV=ND2+NP
    newprint=.false.

    ! if(old) then
    call RESET_APERTURE_FLAG
    call LIEINIT(NO1,NV,ND1,NDPT1)   !,0
   
    call set_in_tpsa( NO,ND,ND2,NP,NDPT,NV,old)
    call set_in_tpsalie( NO,ND,ND2,NP,NDPT,NV,old)


    CALL ASSIGN
    !    CALL ASSIGNMAP
    call alloc(varf1)
    call alloc(varf2)
    npara_fpp=0
    NPARA_original=0
 

  end subroutine init_tpsa




  subroutine DATERMINATE()
    implicit none
    logical present_tpsa
    present_tpsa=lingyun_yang
    if(last_tpsa==1) then
       lingyun_yang=.true.
       CALL DEASSIGN
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       CALL DEASSIGN
    endif
    lingyun_yang=present_tpsa

    !    CALL DEASSIGNMAP
    !   IF(.NOT.OLD) then
    !      CALL DE_initialize_da
    !   endif

  end subroutine DATERMINATE


end module tpsalie_analysis






