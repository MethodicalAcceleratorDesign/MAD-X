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
  private df_map,fd_map,map_df,map_fd,MAPnormal  !,one_map,map_one
  private  resta,tares,alloctares,killtares
  private allocgen,KILLgen   !,ASSgen   ,ADDID,IDADD
  private EQUALgenMAP,EQUALMAPgen,pushgen
  !private rotf,rotsymp,equalgengen
  integer,private::NO,ND,ND2,NP,NDPT,NV
  logical(lp),private::old
  logical(lp) imaxflag

  integer, private :: ndc,ndc2,ndt,iref,itu,iflow,jtune,nres !,idpr
  integer, private,dimension(ndim)::nplane,idsta,ista
  real(dp), private,dimension(ndim)::dsta,sta,angle,rad,ps,rads
  real(dp), private,dimension(ndim,nreso)::mx
  real(dp), private :: stmem(ndim)
  private REXT
  private XGAM,XGBM,COMCFU
  private DFILT,FILT,reelflo_g
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
   !  MODULE PROCEDURE one_map
!     MODULE PROCEDURE map_one
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

!  FUNCTION one_map(S1,S2)
!    implicit none
!    TYPE (damap) one_map ,junk
!    TYPE (damap), INTENT (IN) :: S2
!    TYPE (onelieexponent), INTENT (IN) :: S1
!    integer localmaster
!    IF(.NOT.C_%STABLE_DA) then
!       one_map%v%i=0
!      RETURN
!    endif
!    localmaster=master
!
!    call checkdamap(s2)
!    call assdamap(one_map)
!
!    call alloc(junk)
!
!    junk=1
!
!    junk=s1
!    one_map=junk*s2
!
!    call kill(junk)
!    master=localmaster
!
!  END FUNCTION one_map

!  FUNCTION map_one(S2,S1)
!    implicit none
!    TYPE (damap) map_one,junk
!    TYPE (damap), INTENT (IN) :: S2
!    TYPE (onelieexponent), INTENT (IN) :: S1
!    integer localmaster
!    IF(.NOT.C_%STABLE_DA) then
!       map_one%v%i=0
!      RETURN
!    endif
!    localmaster=master
!
!    call checkdamap(s2)
!    call assdamap(map_one)
!
!    call alloc(junk)
!
!    junk=1
!
!    junk=s1
!
!    map_one=s2*junk   ! 2002.12.6
!
!    call kill(junk)
!    master=localmaster
!
!  END FUNCTION map_one


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

    if(old_package) then
     call gofix(m%v%i,a1%v%i,a1i%v%i,nord)
    else
     call gofix_g(m%v,a1%v,a1i%v,nord)
    endif

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
  !  if(s2%normal%linear%V(1)%i==0)  call crap1("normalMAP 1") !call allocw(s2%normal%linear%V(1))  ! changed
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
        if(old_package) then

    CALL MAPNORMF(junk%V%i,s2%a%nonlinear%v%i,s2%a%linear%v%i,&
         & s2%A1%v%i,s2%normal%linear%v%i,s2%normal%nonlinear%v%i,s2%NORD,s2%jtune)

    s2%normal%linear=(s2%normal%linear).sub.1
    S2%normal%pb=S2%normal%nonlinear
    S2%a%pb=S2%a%nonlinear

    CALL DHDJFLO(S2%NORMAL%NONLINEAR%V%i,S2%DHDJ%V%i)
        else
     CALL MAPNORMF_g(junk%V ,s2%a%nonlinear%v ,s2%a%linear%v ,&
         & s2%A1%v ,s2%normal%linear%v ,s2%normal%nonlinear%v ,s2%NORD,s2%jtune)

    s2%normal%linear=(s2%normal%linear).sub.1
    S2%normal%pb=S2%normal%nonlinear
    S2%a%pb=S2%a%nonlinear

    CALL DHDJFLO_g(S2%NORMAL%NONLINEAR%V,S2%DHDJ%V)

        endif

    call GETTURA(s2%TUNE,S2%DAMPING)

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
    call intd(s2%m%v,s2%h,1.0_dp)
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
    call difd(s2%h,w%v,1.0_dp)
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
    integer(1) ko_
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

!     old=old_package
    old=log1
    old_package=old
    NO=NO1
    ND=ND1
    ND2=ND1*2
    NP=NP1
    NDPT=NDPT1
    NV=ND2+NP
    newprint=.false.
     if(old) then
    call LIEINIT(NO1,NV,ND1,NDPT1)   !,0
    else
    call LIEINIT(NO1,NV,ND1,NDPT1,da_init=.false.)
     if(allocated(vo_berz)) deallocate(vo_berz)

     if(allocated(mo_gtpsa)) deallocate(mo_gtpsa)
     allocate(vo_berz(1:nv))
     allocate(mo_gtpsa(1:nv))
     mo_gtpsa=0
     do i=1,nv
      vo_berz(i)=no
     enddo
      ko_=0
! GTPSA REMOVED !      d_berz=mad_desc_newv(nv,  vo_berz, 0,ko_)
    endif
    ! if(old) then
 !   call LIEINIT(NO1,NV,ND1,NDPT1)   !,0
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
    integer(1) ko_

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
!     old=old_package
    old=log1
    old_package=old
    nd1=0
    ndpt1=0
    NO=NO1
    ND=ND1
    ND2=ND1*2
    NP=NP1
    NDPT=NDPT1
    NV=ND2+NP
    newprint=.false.

    call RESET_APERTURE_FLAG
     if(old) then
    call LIEINIT(NO1,NV,ND1,NDPT1)   !,0
    else
    call LIEINIT(NO1,NV,ND1,NDPT1,da_init=.false.)
     if(allocated(vo_berz)) deallocate(vo_berz)
     if(allocated(mo_gtpsa)) deallocate(mo_gtpsa)
     allocate(vo_berz(1:nv))
     allocate(mo_gtpsa(1:nv))
     mo_gtpsa=0
     do i=1,nv
      vo_berz(i)=no
     enddo
      ko_=0
! GTPSA REMOVED !      d_berz=mad_desc_newv(nv, vo_berz, 0,ko_)
    endif
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


!!!!!!!!!!!!!!!!!!! lielib for gtpsa

  subroutine mapnormf_g(x,ft,a2,a1,xy,h,nord,isi)
    implicit none
    integer ij,isi,nord
    type(taylor),dimension(ndim2)::a1i,a2i
    type(taylor),dimension(:)::x,a1,a2,ft,xy,h
    real(dp),dimension(ndim)::angle,rad,st,p
    if(.not.c_%stable_da) return

    call alloc(a1i,nd2) !  ,'A1I       ')
    call alloc(a2i,nd2) !  ,'A2I       ')
    !     frank/etienne
    do itu=1,ndim
       angle(itu)=0.0_dp
       p(itu)=0.0_dp
       st(itu)=0.0_dp
       rad(itu)=0.0_dp
       ps(itu)=0.0_dp
       rads(itu)=0.0_dp
    enddo
    jtune=isi
    call dacopd_g(x,xy)
    ! go to fix point in the parameters + pt to order nord>=1
    if(nv>nd2.or.ndc==1) then
       call gofix_g(xy,a1,a1i,nord)
       call simil_g(a1i,xy,a1,xy)
    else    !  this "if" was added to remove crashes when y-=plane is nearly identity
       call etini_g(a1)   ! in stochastic kick calculations
       call etini_g(a1i)  ! 2002.10.20
    endif
    ! linear part
    call midbflo_g(xy,a2,a2i,angle,rad,st)
    do ij=1,nd-ndc
       p(ij)=angle(ij)*(st(ij)*(twopii-1.0_dp)+1.0_dp)
    enddo
    stmem=st
    if(ndc.eq.1) p(nd)=angle(nd)


    do ij=1,nd       !  -ndc    Frank
       ps(ij)=p(ij)
       rads(ij)=rad(ij)
    enddo
    call initpert(st,angle,rad)
    call simil_g(a2i,xy,a2,xy)
    call dacopd_g(xy,a2i)
    call orderflo_g(h,ft,xy,angle,rad)
    do ij=1,nd-ndc
       p(ij)=angle(ij)
       !       if(angle(ij).gt.pi.and.st(ij).gt.zero)  then  !.and.itu.eq.1)then
       !          p(ij)=angle(ij)-twopi
       !          w_p=0
       !          w_p%nc=1
       !          w_p%fc='((1X,A72))'
       !          write(w_p%c(1),'(i4,a27,g21.14)') ij,' TH TUNE MODIFIED IN H2 TO ',p(ij)*twopii
       !          CALL WRITE_a
       !       endif
    enddo
    call h2pluflo_g(h,p,rad)
    !      CALL TAKED(A2I,1,XY)
    call taked_g(a2i,1,a1i)
! GTPSA REMOVED !call fpp_mad_tpsa_compose(nd2,xy,nd2,a1i,xy)
  !  call etcct(xy,a1i,xy)


    call kill(a2i,nd2)
    call kill(a1i,nd2)
    return
  end subroutine mapnormf_g

  subroutine orderflo_g(h,ft,x,ang,ra)
    implicit none
    !-   NONLINEAR NORMALIZATION PIECE OF MAPNORMF
    integer k
    type(taylor),dimension(ndim2)::w,v,rel,roi,b1,b5,b6,b9
    type(taylor),dimension(:)::x,h,ft
    real(dp),dimension(ndim)::ang,ra
    if(.not.c_%stable_da) return

    call alloc(w,nd2) !  ,'W         ')
    call alloc(v,nd2) !  ,'V         ')
    call alloc(rel,nd2) !  ,'REL       ')
    call alloc(roi,nd2) !  ,'ROI       ')
    call alloc(b1,nd2) !  ,'B1        ')
    call alloc(b5,nd2) !  ,'B5        ')
    call alloc(b6,nd2) !  ,'B6        ')
    call alloc(b9,nd2) !  ,'B9        ')

    call rotiflo_g(roi,ang,ra)
     call etini_g(rel)
    call kill(h)
    call alloc(ft,nd2)
! GTPSA REMOVED !call fpp_mad_tpsa_compose(nd2,x,nd2,roi,x)
!    call etcct(x,roi,x)
    do k=2,no
       ! IF K>2 V = H(K)^-1 X(K)
       call facflod_g(h,x,v,2,k-1,-1.0_dp,-1)
       ! EXTRACTING K TH DEGREE OF V ----> W
       call taked_g(v,k,w)
       !  write(16,*) "$$$$$$$$  K  $$$$$$$$$$", k
       ! W = EXP(B5) + ...
       call dacopd_g(w,b5)
       !      CALL INTD(W,B5,-one)
       ! B5 ON EXIT IS THE NEW CONTRIBUTION TO H
       ! B6 IS THE NEW CONTRIBUTION TO FT
       call nuanaflo_g(b5,b6)
       call dalind_g(b5,1.0_dp,h,1.0_dp,b1)
       call dacopd_g(b1,h)
       ! EXP(B9) = EXP( : ROTI B6 :)
       !call trxflo(b6,b9,roi)
      call trxflo_g(b6,b9,roi)

       ! V = EXP(-B6) REL
       call facflod_g(b6,rel,v,k,k,-1.0_dp,1)
       ! W = V o X
  !     call etcct(v,x,w)
! GTPSA REMOVED !call fpp_mad_tpsa_compose(nd2,v,nd2,x,w)
       if(lielib_print(5)==1) then

          write(6,'(a13,i4)') ' ORDERFLO K= ', k

       endif
       ! X = EXP(B9) W
       call facflod_g(b9,w,x,k,k,1.0_dp,1)
       ! B6 IS THE NEW CONTRIBUTION TO FT
       call dalind_g(b6,1.0_dp,ft,1.0_dp,b1)
       call dacopd_g(b1,ft)
    enddo
    call kill(b9,nd2)
    call kill(b6,nd2)
    call kill(b5,nd2)
    call kill(b1,nd2)
    call kill(roi,nd2)
    call kill(rel,nd2)
    call kill(v,nd2)
    call kill(w,nd2)
    return
  end subroutine orderflo_g

  subroutine nuanaflo_g(h,ft)
    implicit none
    ! RESONANCE DENOMINATOR OPERATOR (1-R^-1)^-1
    integer i
    type(taylor),dimension(:)::h,ft
    type(taylor),dimension(ndim2)::br,bi,c,ci
    type(taylor),dimension(2)::t1,t2
    if(.not.c_%stable_da) return

    call alloc(br,nd2)
    call alloc(bi,nd2)
    call alloc(c,nd2)
    call alloc(ci,nd2)
    call alloc(t1,2)
    call alloc(t2,2)

    call ctorflo_g(h,br,bi)

    ! FILTERING RESONANCES AND TUNE SHIFTS
    ! ASSUMING REALITY I.E. B(2*I-1)=CMPCJG(B(2*I))

    do i=1,nd2
       iflow=i
       call cfu(br(i),filt,c(i))
       call cfu(bi(i),filt,ci(i))
    enddo
    call rtocflo_g(c,ci,h)

    do i=1,nd2

       iflow=i
       call cfu(br(i),dfilt,br(i))
       call cfu(bi(i),dfilt,bi(i))
    enddo
    !  NOW WE MUST REORDER C AND CI TO SEPARATE THE REAL AND IMAGINARY PART
    ! THIS IS NOT NECESSARY WITH :H: OPERATORS

    do i=1,nd2
       t1(1)=br(i)
       t1(2)=bi(i)
       t2(1)=c(i)
       t2(2)=ci(i)
       iflow=i
       call comcfu(t1,xgam,xgbm,t2)
       c(i)=t2(1)
       ci(i)=t2(2)
    enddo

    call rtocflo_g(c,ci,ft)

    call kill(br,nd2)
    call kill(bi,nd2)
    call kill(c,nd2)
    call kill(ci,nd2)
    call kill(t1,2)
    call kill(t2,2)
    return
  end subroutine nuanaflo_g

  subroutine comcfu(b,f1,f2,c)
    implicit none
    ! Complex dacfu
    type(taylor),dimension(2)::b,c
    type(taylor),dimension(4)::t
    real(dp),external::f1,f2

    if(.not.c_%stable_da) return



    call alloc(t,4)

    call cfu(b(1),f1,t(1))
    call cfu(b(1),f2,t(2))
    call cfu(b(2),f1,t(3))
    call cfu(b(2),f2,t(4))
    c(1)=t(1)-t(4)
   ! call dasub(t(1),t(4),c(1))
    c(2)=t(2)-t(4)
   ! call daadd(t(2),t(3),c(2))

    call kill(t,4)

    return
  end subroutine comcfu

  real(dp) function xgam(j)
    implicit none
    ! XGAM AND XGBM ARE THE EIGENVALUES OF THE OPERATOR NEWANAFLO
    integer i,ic,ij,ik
    !      INTEGER J(NTT),JJ(NDIM),JP(NDIM)
    integer,dimension(:)::j
    integer,dimension(ndim)::jj,jp
    real(dp) ad,ans,as,ex,exh
    if(.not.c_%stable_da) return

    xgam=0.0_dp
    ad=0.0_dp
    as=0.0_dp
    ic=0
    do i=1,nd-ndc
       ik=2*i-1
       ij=2*i
       jp(i)=j(ik)+j(ij)
       jj(i)=j(ik)-j(ij)
       if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
          jp(i)=jp(i)-1
       endif
       ic=ic+iabs(jj(i))
    enddo

    do i=1,nd-ndc
       ad=dsta(i)*REAL(jj(i),kind=DP)*angle(i)-REAL(jp(i),kind=DP)*rad(i)+ad
       as=sta(i)*REAL(jj(i),kind=DP)*angle(i)+as
    enddo

    exh=EXP(ad/2.0_dp)
    ex=exh**2
    ans=4.0_dp*ex*(SINH(ad/2.0_dp)**2+SIN(as/2.0_dp)**2)
    if(ans.eq.0.0_dp) then
       print*,"NormalForm makes no sense!"
       print*,"no,nv,nd,nd2",no,nv,nd,nd2
       print*,"ndc,ndc2,ndt,ndpt",ndc,ndc2,ndt,ndpt
       stop
    endif
    xgam=2.0_dp*(-exh*SINH(ad/2.0_dp)+ex*SIN(as/2.0_dp)**2)/ans

    return
  end function xgam

  real(dp) function xgbm(j)
    implicit none
    integer i,ic,ij,ik
    real(dp) ad,ans,as,ex,exh
    !      INTEGER J(NTT),JJ(NDIM),JP(NDIM)
    integer,dimension(:)::j
    integer,dimension(ndim)::jj,jp
    if(.not.c_%stable_da) return

    xgbm=0.0_dp
    ad=0.0_dp
    as=0.0_dp
    ic=0
    do i=1,nd-ndc
       ik=2*i-1
       ij=2*i
       jp(i)=j(ik)+j(ij)
       jj(i)=j(ik)-j(ij)
       if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
          jp(i)=jp(i)-1
       endif
       ic=ic+iabs(jj(i))
    enddo

    do i=1,nd-ndc
       ad=dsta(i)*REAL(jj(i),kind=DP)*angle(i)-REAL(jp(i),kind=DP)*rad(i)+ad
       as=sta(i)*REAL(jj(i),kind=DP)*angle(i)+as
    enddo

    exh=EXP(ad/2.0_dp)
    ex=exh**2
    ans=4.0_dp*ex*(SINH(ad/2.0_dp)**2+SIN(as/2.0_dp)**2)
    if(ans.eq.0.0_dp) then
       print*,"NormalForm makes no sense!"
       print*,"no,nv,nd,nd2",no,nv,nd,nd2
       print*,"ndc,ndc2,ndt,ndpt",ndc,ndc2,ndt,ndpt
       stop
    endif
    xgbm=SIN(as)*ex/ans

    return
  end function xgbm

  real(dp) function filt(j)
    implicit none
    !  PROJECTION FUNCTIONS ON THE KERNEL ANMD RANGE OF (1-R^-1)
    !-  THE KERNEL OF (1-R^-1)
    integer i,ic,ic1,ic2,ij,ik,ji
    !      INTEGER J(NTT),JJ(NDIM)
    integer,dimension(:)::j
    integer,dimension(ndim)::jj
    if(.not.c_%stable_da) return

    filt=1.0_dp

    ic=0
    do i=1,nd-ndc
       ik=2*i-1
       ij=2*i
       jj(i)=j(ik)-j(ij)
       if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
       endif
       ic=ic+iabs(jj(i))
    enddo

    if(ic.eq.0.and.jtune.eq.0) return

    do i=1,nres
       ic1=1
       ic2=1
       do ji=1,nd-ndc
          if(mx(ji,i).ne.jj(ji)) ic1=0
          if(mx(ji,i).ne.-jj(ji)) ic2=0
          if(ic1.eq.0.and.ic2.eq.0) goto 3
       enddo
       return
3      continue
    enddo

    filt=0.0_dp
    return
  end function filt

  real(dp) function dfilt(j)
    implicit none
    !-  THE RANGE OF (1-R^-1)^1
    !- CALLS FILT AND EXCHANGES 1 INTO 0 AND 0 INTO 1.
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    real(dp) fil
    if(.not.c_%stable_da) return

    fil=filt(j)
    if(fil.gt.0.5_dp) then
       dfilt=0.0_dp
    else
       dfilt=1.0_dp
    endif
    return
  end function dfilt

  subroutine reelflo_g(c,ci,f,fi)
    implicit none
    ! DOES THE SPINOR PART IN RTOCFLO
    integer i
!    integer,dimension(:)::c,ci,f,fi
!    integer,dimension(ndim2)::e,ei
    type(taylor),dimension(:)::c,ci,f,fi
    type(taylor),dimension(ndim2)::e,ei

    if(.not.c_%stable_da) return


    call alloc(e,nd2)
    call alloc(ei,nd2)

    do i=1,nd-ndc
     e(2*i-1)=c(2*i-1)*0.5_dp+c(2*i)*0.5_dp
    !   call dalin(c(2*i-1),0.5_dp,c(2*i),0.5_dp,e(2*i-1))
     ei(2*i-1)=ci(2*i-1)*0.5_dp+ci(2*i)*0.5_dp
    !   call dalin(ci(2*i-1),0.5_dp,ci(2*i),0.5_dp,ei(2*i-1))
       if(ista(i).eq.1) then
        e(2*i)=ci(2*i-1)*0.5_dp-ci(2*i)*0.5_dp
    !      call dalin(ci(2*i-1),0.5_dp,ci(2*i),-0.5_dp,e(2*i))
          ei(2*i)=-c(2*i-1)*0.5_dp+c(2*i)*0.5_dp
      !    call dalin(c(2*i-1),-0.5_dp,c(2*i),0.5_dp,ei(2*i))
       else
            ei(2*i)=ci(2*i-1)*0.5_dp-ci(2*i)*0.5_dp
       !   call dalin(ci(2*i-1),0.5_dp,ci(2*i),-0.5_dp,ei(2*i))
            e(2*i)=c(2*i-1)*0.5_dp-c(2*i)*0.5_dp
       !   call dalin(c(2*i-1),0.5_dp,c(2*i),-0.5_dp,e(2*i))
       endif
    enddo

    do i=nd2-ndc2+1,nd2
    e(i)=c(i)
     !  call dacop(c(i),e(i))
    ei(i)=ci(i)
     !  call dacop(ci(i),ei(i))
    enddo

    call dacopd_g(e,f)
    call dacopd_g(ei,fi)

    call kill(e,nd2)
    call kill(ei,nd2)

    return
  end subroutine reelflo_g

  subroutine ctord_g(c,cr,ci)
    implicit none
    ! ROUTINES USED IN THE INTERMEDIATE STEPS OF CTORFLO AND RTOCFLO
    ! SAME AS CTOR  OVER ARRAYS CONTAINING ND2 COMPONENTS
    ! ROUTINE USEFUL IN INTERMEDIATE FLOW CHANGE OF BASIS
    integer i
    type(taylor),dimension(:)::c,ci,cr
    if(.not.c_%stable_da) return

    do i=1,nd2
       call ctor_g(c(i),cr(i),ci(i))
    enddo
    return
  end subroutine ctord_g

  subroutine ctor_g(c1,r2,i2)
    implicit none
    ! CHANGES OF BASIS
    !   C1------> R2+I R1
    type(taylor) c1,r2,i2,b1,b2
    type(taylor),dimension(ndim2)::x

    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(b2)
    call alloc(x,nd2) !  ,'X         ')

    call ctoi_g(c1,b1)
    call etcjg_g(x)
    call trx_g(b1,b2,x)
    r2=b1*0.5_dp+b2*0.5_dp
!    call dalin(b1,0.5_dp,b2,0.5_dp,r2)
     i2=b1*0.5_dp-b2*0.5_dp
!    call dalin(b1,0.5_dp,b2,-0.5_dp,i2)



    call alloc(x,nd2)
    call alloc(b2)
    call alloc(b1)
    return
  end subroutine ctor_g


  subroutine ctoi_g(f1,f2) ! no flip needed
    implicit none
    type(taylor) f1,f2,b1
    type(taylor),dimension(ndim2)::x
    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(x,nd2) !  ,'X         ')

    call cpart_g(f1,b1)
    call etctr_g(x)
    call trx_g(b1,f2,x)
    call kill(x,nd2)
    call kill(b1)
    return
  end subroutine ctoi_g


  subroutine etctr_g(x) ! no flip needed
    implicit none
    integer i
    type(taylor),dimension(:)::x
    type(taylor),dimension(ndim2)::rel
    if(.not.c_%stable_da) return

    call alloc(rel,nd2) !  ,'REL       ')

    call etini_g(rel)
    call etini_g(x)
    do i=1,nd-ndc
       call dalin_g(rel(2*i-1),0.5_dp,rel(2*i),0.5_dp,x(2*i-1))
       call dalin_g(rel(2*i-1),0.5_dp,rel(2*i),-0.5_dp,x(2*i))
    enddo
    call kill(rel,nd2)
    return
  end subroutine etctr_g

  subroutine dalin_g(x1,r1,x2,r2,x3) ! no flip needed
    implicit none
    type(taylor) x1,x2,x3
    real(dp) r1,r2
    if(.not.c_%stable_da) return

     x3=x1*r1+x2*r2

  end subroutine dalin_g


  subroutine etcjg_g(x) ! no flip needed
    implicit none
    integer i
    type(taylor),dimension(:)::x
    type(taylor),dimension(ndim2)::rel
    if(.not.c_%stable_da) return

    call alloc(rel,nd2) !  ,'REL       ')

    call etini_g(rel)
    call etini_g(x)
    do i=1,nd-ndc
       if(ista(i).eq.1) then
          call dacop_g(rel(2*i-1),x(2*i))
          call dacop_g(rel(2*i),x(2*i-1))
       else
          call dacop_g(rel(2*i-1),x(2*i-1))
          call dacop_g(rel(2*i),x(2*i))
       endif
    enddo
    call kill(rel,nd2)
    return
  end subroutine etcjg_g

  subroutine dacop_g(a1,a2)
    implicit none
    type(taylor) a1,a2
   a2=a1
  end   subroutine dacop_g

  subroutine trx_g(a1,a2,v)
    implicit none
     type(taylor) v(:)
    type(taylor) a1,a2
    type(taylor) t(2)

     call alloc(t,2)

       t(1)=a1
       t(2)=a2
! GTPSA REMOVED !        call fpp_mad_tpsa_compose(1,t(1:1),nd2,v,t(2:2))
       a2=t(2)

     call kill(t,2)

 end   subroutine trx_g

  subroutine rtoc_g(r1,i1,c2)
    implicit none
    !  INVERSE OF CTOR
    type(taylor) c2,r1,i1,b1

    if(.not.c_%stable_da) return

    call alloc(b1)
    b1=r1+i1
!    call daadd(r1,i1,b1)
    call itoc_g(b1,c2)

    call kill(b1)

    return
  end subroutine rtoc_g

  subroutine itoc_g(f1,f2) ! no flip needed
    implicit none
    type(taylor) f1,f2,b1
    type(taylor),dimension(ndim2)::x
    if(.not.c_%stable_da) return

    call alloc(b1)
    call alloc(x,nd2)

    call etrtc_g(x)
    call trx_g(f1,b1,x)
    call cpart_g(b1,f2)

    call kill(x,nd2)
    call kill(b1)
    return
  end subroutine itoc_g

  real(dp) function rext(j) ! no flip needed
    implicit none
    integer i,lie,mo
    integer,dimension(:)::j
    if(.not.c_%stable_da) return

    lie=0
    do i=1,nd-ndc
       lie=ista(i)*j(2*i)+lie
    enddo
    mo=mod(lie,4)+1
    !frs
    select case(mo)
    case(1,4)
       rext = 1.0_dp
    case(2,3)
       rext = -1.0_dp
    end select
    return

  end function rext

  subroutine cpart_g(h,ch)  ! no flip needed
    implicit none
    type(taylor) h,ch
    if(.not.c_%stable_da) return

    call cfu(h,rext,ch)
    return
  end subroutine cpart_g

  subroutine etrtc_g(x) ! no flip needed
    implicit none
    integer i
    type(taylor),dimension(:)::x
    type(taylor),dimension(ndim2)::rel
    if(.not.c_%stable_da) return

    call alloc(rel,nd2) !  ,'REL       ')

    call etini_g(rel)
    call etini_g(x)

    do i=1,nd-ndc
    x(2*i-1)=rel(2*i-1)+rel(2*i)
   !    call daadd(rel(2*i-1),rel(2*i),x(2*i-1))
    x(2*i)=rel(2*i-1)-rel(2*i)
   !    call dasub(rel(2*i-1),rel(2*i),x(2*i))
    enddo

    call kill(rel,nd2)
    return
  end subroutine etrtc_g

  subroutine dalind_g(h,rh,ht,rt,hr)
    implicit none
    integer i
    type(taylor),dimension(:)::h,ht,hr
    type(taylor),dimension(ndim2)::b
    real(dp) rh,rt
    if(.not.c_%stable_da) return

    call alloc(b,nd2) !  ,'B         ')

    do i=1,nd2
     b(i)=h(i)*rh+ht(i)*rt
  !     call dalin(h(i),rh,ht(i),rt,b(i))
    enddo
    call dacopd_g(b,hr)
    call kill(b,nd2)
    return
  end subroutine dalind_g

  subroutine etini_g(x)
    implicit none
    !  X=IDENTITY
    integer i
    type(taylor),dimension(:)::x
    if(.not.c_%stable_da) return

    do i=1,nd2
      x(i)=1.0_dp.mono.i
     !  call davar(x(i),0.0_dp,i)
    enddo
    return
  end subroutine etini_g


  subroutine ctorflo_g(c,dr,di)
    implicit none
    ! FLOW CTOR
    type(taylor),dimension(:)::dr,di,c
    if(.not.c_%stable_da) return

    call ctord_g(c,dr,di)
    call resvec_g(dr,di)

    return
  end subroutine ctorflo_g

  subroutine rtocflo_g(dr,di,c)
    implicit none
    ! FLOW RTOC
    type(taylor),dimension(:)::dr,di,c
    type(taylor),dimension(ndim2)::er,ei
    if(.not.c_%stable_da) return

    call alloc(er,nd2)
    call alloc(ei,nd2)

    call reelflo_g(dr,di,er,ei)
    call rtocd_g(er,ei,c)

    call kill(er,nd2)
    call kill(ei,nd2)

    return
  end subroutine rtocflo_g

  subroutine daclrd_g(h)
    implicit none
    ! clear a map : a vector of nd2 polynomials
    integer i
    type(taylor),dimension(:)::h
    if(.not.c_%stable_da) return

    do i=1,nd2
        h(i)=0.0_dp
    enddo
    return
  end subroutine daclrd_g


  subroutine simil_g(a,x,ai,y)
    implicit none
    !  Y= AoXoAI
    type(taylor),dimension(:)::x,y,a,ai
    type(taylor),dimension(ndim2)::w,v
    if(.not.c_%stable_da) return

    call alloc(w,nd2) ! ,'W         ')
    call alloc(v,nd2) ! ,'V         ')

!    call etcct(a,x,w)
!    call etcct(w,ai,v)
! GTPSA REMOVED !call fpp_mad_tpsa_compose(nd2,a,nd2,x,w)
! GTPSA REMOVED !call fpp_mad_tpsa_compose(nd2,w,nd2,ai,v)

    call dacopd_g(v,y)

    call kill(v,nd2)
    call kill(w,nd2)
    return
  end subroutine simil_g

  subroutine gofix_g(xy,a1,a1i,nord)
    implicit none
    ! GETTING TO THE FIXED POINT AND CHANGING TIME APPROPRIATELY IN THE
    ! COASTING BEAM CASE
    !****************************************************************
    ! X = A1 XY A1I WHERE X IS TO THE FIXED POINT TO ORDER NORD
    ! for ndpt not zero, works in all cases. (coasting beam: eigenvalue
    !1 in Jordan form)
    !****************************************************************
    integer i,nord
    type(taylor),dimension(:)::xy,a1,a1i
    type(taylor),dimension(ndim2)::x,w,v,rel
    real(dp) xic
    if(.not.c_%stable_da) return



    call alloc(x,nd2) !  ,  'X         ')
    call alloc(w,nd2) !  ,  'W         ')
    call alloc(v,nd2) !  ,  'V         ')
    call alloc(rel,nd2) !  ,'REL       ')


    ! COMPUTATION OF A1 AND A1I USING DAINV



    call etini_g(rel)

    !    call danot(nord)

    call etini_g(v)

    do i=1,nd2-ndc2
       !       call dacop(xy(i),x(i))
        x(i)=x(i).cut.(nord+1)
     !  call datrunc(xy(i),nord+1,x(i))
       v(i)=x(i)-rel(i)
     !  call dalin(x(i),1.0_dp,rel(i),-1.0_dp,v(i))
    enddo
! GTPSA REMOVED !call fpp_mad_tpsa_minv(nd2,v,w)
!    call etinv(v,w)
    call datruncd_g(w,nord+1,w)
    call daclrd_g(x)
    if(ndc.eq.1) then
        x(ndpt)=1.0_dp.mono.ndpt
!       call davar(x(ndpt),0.0_dp,ndpt)
    endif
!    call etcct(w,x,v)
! GTPSA REMOVED !call fpp_mad_tpsa_compose(nd2,w,nd2,x,v)
    if(ndc.eq.1) then
        v(nd2)=0.0_dp
        v(nd2-ndc)=0.0_dp
 !      call daclr(v(nd2))
 !      call daclr(v(nd2-ndc))
    endif
    call dalind_g(rel,1.0_dp,v,1.0_dp,a1)
    call dalind_g(rel,1.0_dp,v,-1.0_dp,a1i)

    if(ndpt.ne.0) then

       !  CORRECTIONS
       call daclrd_g(w)
       call daclrd_g(v)
       call daclrd_g(x)

       do i=1,nd2-ndc2
            w(i)=a1(i)-rel(i)
!          call dalin(a1(i),1.0_dp,rel(i),-1.0_dp,w(i))
       enddo

       !      COMPUTE Deta/Ddelta
       call dacopd_g(w,a1)

       do i=1,nd2-ndc2
            w(i)=w(i).d.ndpt
      !    call dader(ndpt,w(i),w(i))
       enddo
       !      COMPUTE J*Deta/dDELTA

       do i=1,nd-ndc
        v(2*i-1) = w(2*i)
        v(2*i) = -w(2*i-1)
!          call dacmu(w(2*i),1.0_dp,v(2*i-1) )
!          call dacmu(w(2*i-1),-1.0_dp,v(2*i) )
       enddo

       xic=(-1)**(ndt)

       do i=1,nd2-ndc2
        x(1)=v(i)*rel(i)
        w(ndt)=x(1)+w(ndt)
        w(i)=a1(i)
     !     call damul(v(i),rel(i),x(1))
     !     call daadd(x(1),w(ndt),w(ndt))
     !     call dacop(a1(i),w(i))
       enddo

        w(ndt)=xic*w(ndt)
!       call dacmu(w(ndt),xic,w(ndt))

       call expflod_g(w,rel,a1,1e-7_dp,10000)
       ! END OF  CORRECTIONS

       call datruncd_g(a1,nord+1,a1)
   !    call etinv_g(a1,a1i)
! GTPSA REMOVED !      call fpp_mad_tpsa_minv(nd2,a1,a1i)
       call datruncd_g(a1i,nord+1,a1i)
    endif


    !    call danot(no)

    call kill(rel,nd2)
    call kill(v,nd2)
    call kill(w,nd2)
    call kill(x,nd2)
    return
  end subroutine gofix_g

  subroutine datruncd_g(h,io,ht)
    implicit none
    !    H goes into HT  (nd2 array)
    integer i
    type(taylor),dimension(:)::h,ht
    integer io
    if(.not.c_%stable_da) return

    do i=1,nd2
        ht(i)=h(i).cut.io
   !    call datrunc(h(i),io,ht(i))
    enddo
    return
  end subroutine datruncd_g

  subroutine rotiflo_g(roi,ang,ra)
    implicit none
    ! CREATES  R^-1
    integer i
    real(dp) ch,sh,sim,simv,xx
    integer,dimension(lnv)::j
    type(taylor),dimension(:)::roi
    real(dp),dimension(ndim)::co,si,ang,ra
    if(.not.c_%stable_da) return

    !    do i=1,10
    j=0
    !    enddo

    call daclrd_g(roi)
    do i=1,nd-ndc
       xx=EXP(-ra(i))
       if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
          si(i)=-sh*xx
       else
          co(i)=COS(ang(i))*xx
          si(i)=SIN(ang(i))*xx
       endif
    enddo
    do i=1,nd-ndc
       if(ista(i).eq.0)then
          sim=si(i)
       else
          sim=-si(i)
       endif
       j(2*i-1)=1
       call pok(roi(2*i-1),j,co(i))
       simv=-sim
       call pok(roi(2*i),j,simv)
       j(2*i-1)=0
       j(2*i)=1
       simv=-si(i)
       call pok(roi(2*i),j,co(i))
       call pok(roi(2*i-1),j,simv)
       j(2*i)=0
    enddo

    if(ndc.eq.1) then
       j(ndt)=1
       call pok(roi(ndt),j,1.0_dp)
       call pok(roi(ndpt),j,0.0_dp)
       j(ndt)=0
       j(ndpt)=1
       call pok(roi(ndt),j,-ang(nd))
       call pok(roi(ndpt),j,1.0_dp)
       j(ndpt)=0
    endif

    return
  end subroutine rotiflo_g

  subroutine rtocd_g(cr,ci,c)
    implicit none
    !  INVERSE OF CTORD
    integer i
    type(taylor),dimension(:)::c,ci,cr
    if(.not.c_%stable_da) return

    do i=1,nd2
       call rtoc_g(cr(i),ci(i),c(i))
    enddo
    return
  end subroutine rtocd_g

  subroutine h2pluflo_g(h,ang,ra)
    implicit none
    ! POKES IN \VEC{H}  ANGLES AND DAMPING COEFFFICIENTS
    !
    integer i
    integer,dimension(lnv)::j
    type(taylor),dimension(:)::h
    real(dp) r1,r2
    real(dp),dimension(ndim)::ang,ra,st
    if(.not.c_%stable_da) return

    do i=1,nd
       st(i)=2.0_dp*sta(i)-1.0_dp
    enddo

    do i=1,lnv
       j(i)=0
    enddo

    do i=1,nd-ndc
       j(2*i-1)=1
       r1=-ang(i)
       !-----
       call pok(h(2*i),j,r1)

       r2=ra(i)
       call pok(h(2*i-1),j,r2)
       j(2*i-1)=0

       j(2*i)=1
       r1=ang(i)*st(i)
       call pok(h(2*i-1),j,r1)
       call pok(h(2*i),j,r2)
       j(2*i)=0

    enddo

    if(ndpt.eq.nd2-1) then
       j(ndpt)=1
       call pok(h(ndt),j,ang(nd))
    elseif(ndpt.eq.nd2) then
       j(ndpt)=1
       call pok(h(ndt),j,ang(nd))     !!!! correct 2014.2.12 with David Sagan and Chris Mayes
    endif
    return
  end subroutine h2pluflo_g

  subroutine resvec_g(cr,ci)
    implicit none
    ! DOES THE SPINOR PART IN CTORFLO
    integer i
    type(taylor),dimension(:)::ci,cr
    type(taylor),dimension(2)::tr,ti

    if(.not.c_%stable_da) return


    call alloc(tr,2)
    call alloc(ti,2)

    do i=1,nd-ndc
       if(ista(i).eq.1) then

       tr(1)=cr(2*i-1)-ci(2*i)
       ti(1)=ci(2*i-1)+cr(2*i)
       tr(2)=cr(2*i-1)+ci(2*i)
       ti(2)=ci(2*i-1)-cr(2*i)

        cr(2*i-1)=tr(1)
        cr(2*i)  =tr(2)
        ci(2*i-1)=ti(1)
        ci(2*i)  =ti(2)

      !    call dasub(cr(2*i-1),ci(2*i),tr(1))
      !    call daadd(ci(2*i-1),cr(2*i),ti(1))
      !   call daadd(cr(2*i-1),ci(2*i),tr(2))
      !    call dasub(ci(2*i-1),cr(2*i),ti(2))
       !   call dacop(tr(1),cr(2*i-1))
       !   call dacop(tr(2),cr(2*i))
       !   call dacop(ti(1),ci(2*i-1))
      !    call dacop(ti(2),ci(2*i))
       else
        tr(1)=cr(2*i-1)+cr(2*i)
        ti(1)=ci(2*i-1)+ci(2*i)
        tr(2)=cr(2*i-1)-cr(2*i)
        ti(2)=ci(2*i-1)-ci(2*i)
        cr(2*i-1) =tr(1)
        cr(2*i)   =tr(2)
        ci(2*i-1) =ti(1)
        ci(2*i)   =ti(2)

     !     call daadd(cr(2*i-1),cr(2*i),tr(1))
     !     call daadd(ci(2*i-1),ci(2*i),ti(1))
     !     call dasub(cr(2*i-1),cr(2*i),tr(2))
     !     call dasub(ci(2*i-1),ci(2*i),ti(2))
     !     call dacop(tr(1),cr(2*i-1))
     !     call dacop(tr(2),cr(2*i))
     !     call dacop(ti(1),ci(2*i-1))
     !     call dacop(ti(2),ci(2*i))
       endif
    enddo

    !    do i=nd2-ndc2+1,nd2
    !       call dacop(cr(i),dr(i))
    !       call dacop(ci(i),di(i))
    !    enddo


    call kill(tr,2)
    call kill(ti,2)
    return
  end subroutine resvec_g

  subroutine midbflo_g(c,a2,a2i,q,a,st)
    implicit none
    ! LINEAR EXACT NORMALIZATION USING EIGENVALUE PACKAGE OF NERI
    !
    integer i,j
    integer,dimension(lnv)::jx
    type(taylor),dimension(:)::c,a2,a2i
    real(dp) ch,r,shm
    real(dp),dimension(ndim2,ndim2)::cr,sa,sai,cm
    real(dp),dimension(ndim)::st,q,a
    if(.not.c_%stable_da) return

    do i=1,lnv
       jx(i)=0
    enddo

    !     frank/etienne
    do i=1,ndim
       st(i)=0.0_dp
       q(i)=0.0_dp
       a(i)=0.0_dp
    enddo
    !     frank/etienne
    do i=1,ndim2
       !     frank/etienne
       do j=1,ndim2
          sai(i,j)=0.0_dp
          sa(i,j)=0.0_dp
          cm(i,j)=0.0_dp
          cr(i,j)=0.0_dp
       enddo
    enddo

    do i=1,nd2
       do j=1,nd2
          jx(j)=1
          call  pek(c(i),jx,r)
          jx(j)=0
          cm(i,j)=r
       enddo
    enddo

    call mapflol(sa,sai,cr,cm,st)
    do i=1,nd-ndc
       if(st(i)+1e-3_dp.gt.1.0_dp) then
          a(i)=sqrt(cr(2*i-1,2*i-1)**2+cr(2*i-1,2*i)**2)
          q(i)=ARCCOS_lielib(cr(2*i-1,2*i-1)/a(i))
          a(i)=LOGE_lielib(a(i))
          if(cr(2*i-1,2*i).lt.0.0_dp) q(i)=twopi-q(i)
       else
          a(i)=sqrt(cr(2*i-1,2*i-1)**2-cr(2*i-1,2*i)**2)
          ch=cr(2*i-1,2*i-1)/a(i)
          shm=cr(2*i-1,2*i)/a(i)
          !       CH=CH+SQRT(CH**2-one)
          !       q(i)=LOG(CH)
          q(i)=-LOGE_lielib(ch+shm)   ! half integer ???? blows up
          !       IF(cr(2*i-1,2*i).gt.zero) Q(I)=-Q(I)
          a(i)=LOGE_lielib(a(i))
       endif
    enddo
    if(ndc.eq.0) then
       if(time_plane>0) then
          if(new_ndpt) then
             !       do i=1,nd
             !         write(6,*) i,q(i)/twopi
             !        if(st(i)+c_1d_3.gt.one.and.q(i).gt.pi) q(i)=q(i)-twopi
             !          write(6,*) i,q(i)/twopi
             !          pause 77
             !      enddo
             if(st(time_plane)+1e-3_dp.gt.1.0_dp.and.nd.ge.3.and.q(time_plane).gt.pi) q(time_plane)=q(time_plane)-twopi
          else
             if(st(time_plane)+1e-3_dp.gt.1.0_dp.and.nd.ge.3.and.q(time_plane).gt.pi) q(time_plane)=q(time_plane)-twopi
          endif
       endif
    else
       if(new_ndpt) then
          !       do i=1,nd-1
          !        if(st(i)+c_1d_3.gt.one.and.q(i).gt.pi) q(i)=q(i)-twopi
          !       enddo
       endif
       q(nd)=cr(ndt,ndpt)
    endif

    call daclrd_g(a2)
    call daclrd_g(a2i)

    do i=1,nd2
       do j=1,nd2
          jx(j)=1
          r=sa(i,j)
          if(r.ne.0.0_dp)call  pok(a2(i),jx,r)
          jx(j)=1
          r=sai(i,j)
          if(r.ne.0.0_dp)call  pok(a2i(i),jx,r)
          jx(j)=0
       enddo
    enddo

    return
  end subroutine midbflo_g

  subroutine dhdjflo_g(h,t)
    implicit none
    ! CONVENIENT TUNE SHIFT FINDED FOR SYMPLECTIC CASE (NU,DL)(H)=T
    integer i,j1,j2
    type(taylor),dimension(:)::h,t
    type(taylor),dimension(ndim2)::b1,b2
    type(taylor) bb1,bb2

    if(.not.c_%stable_da) return

    call alloc(b1,nd2)
    call alloc(b2,nd2)
    call alloc(bb1)
    call alloc(bb2)

    call ctorflo_g(h,b1,b2)

    do i=1,nd-ndc
        bb1=b2(2*i).k.(2*i)
!       call datra(2*i,b2(2*i),bb1)
        t(i+nd)=bb1*twopii
!       call dacmu(bb1,twopii,t(i+nd))
        bb1=t(i+nd)
!       call dacop(t(i+nd),bb1)
        bb1=0.0_dp
!       call daclr(bb2)
       call rtoc_g(bb1,bb2,bb1)
        t(i)=bb1
 !      call dacop(bb1,t(i))
    enddo

    if(ndpt.ne.0) then
        t(nd)=h(ndt)
  !     call dacop(h(ndt),t(nd))
        t(nd2)=b1(ndt)
 !      call dacop(b1(ndt),t(nd2))
    endif

    call kill(bb2)
    call kill(bb1)
    call kill(b2,nd2)
    call kill(b1,nd2)
    return
  end subroutine dhdjflo_g

end module tpsalie_analysis






