!The Full Polymorphic Package
!Copyright (C) Etienne Forest


MODULE TPSA
  !use newda
  use definition
  use file_handler
    use precision_constants
  IMPLICIT NONE
  public


  private factorial, poly_eval
  integer,private::ndel ,nd2par,nd2part,nd2partt,nd2parfilter
  integer,private,dimension(lnv)::jfil,jfilt

  private equal,DAABSEQUAL,Dequaldacon ,equaldacon ,Iequaldacon  !,AABSEQUAL 2002.10.17
  private pow,powr,powr8,dlogt, GETORDER,CUTORDER,getchar,GETint
  private getdiff,getdATRA  ,mul,dmulsc,dscmul
  private mulsc,scmul,imulsc,iscmul
  private div,ddivsc,dscdiv,divsc,scdiv,idivsc,iscdiv
  private unaryADD,add,daddsc,dscadd,addsc,scadd,iaddsc,iscadd
  private unarySUB,subs,dsubsc,dscsub,subsc,scsub,isubsc,iscsub
  private allocda,KILLda,A_OPT,K_opt
  private dexpt,dcost,dsint,dsqrtt,dtant,datanht,dtanht,c_exp_quaternion
  PRIVATE GETCHARnd2,GETintnd2,dputchar,dputint, filter,check_j,dsinHt,dCOSHt
  private GETintnd2t,print_for_bmad_parse,cdivq,cmulq,csubq,caddq
  PRIVATE DEQUAL,REQUAL,varf,varf001,dputint0  !,CHARINT
  !  PUBLIC VAR,ASS
  private pbbra,full_absT,asstaylor,getcharnd2s,GETintnd2s,GETintk
  private shiftda,shift000,cunaryADDq,cunarySUBq,cinvq,cabsq,cabsq2
  !PRIVATE null_0,ALLOC_U,FILL_N,REFILL_N
  !  public, alloc_uni, null_uni, fill_uni, refill_uni
  private rcmulq,cmulqr,ccmulq,cmulqc
  private fill_uni_r,taylor_cycle_dabnew,taylor_cycle_gtpsa ! new sagan
! quaternion
   private subq,unarysubq,addq,unaryADDq,absq,absq2,mulq,divq,EQUALq,EQUALqr,EQUALqi,powq,printq ,invq
   private EQUALcq,cEQUALqr,cEQUALqi,cPOWq,EQUALq_cq,EQUALcq_q
  private NO,ND,ND2,NP,NDPT,NV
  integer NP,NO,ND,ND2,NDPT,NV
  integer, TARGET :: NSPIN=0
  integer, TARGET :: SPIN_pos=0
  private old
  logical(lp) old
  logical(lp),target  :: real_warning =.true.

  PRIVATE null_it,Set_Up,de_Set_Up,LINE_L,RING_L,kill_DALEVEL,dealloc_DASCRATCH,set_up_level
  private insert_da,append_da,GETINTegrate
INTEGER, private, PARAMETER :: I4B = SELECTED_INT_KIND(9)
!INTEGER, private, PARAMETER :: DP = KIND(1.0D0)

logical:: switch_bessel=.true.
private norm_bessel_Ir,nbit,nbittr,nbitrt,etienne_bessel_Ir,etienne_bessel_It,etienne_bessel_Itr,etienne_bessel_Irt
private nbitreal,nbittaylor,nbittaylorrt,nbittaylortr,CFU000_new
  type(dalevel) scratchda(ndumt)   !scratch levels of DA using linked list
  real(dp), pointer :: tn0(:)=>null()







 INTERFACE nbi
     MODULE PROCEDURE nbitreal
     MODULE PROCEDURE nbittaylor
     MODULE PROCEDURE nbittaylorrt
     MODULE PROCEDURE nbittaylortr
  END INTERFACE

 INTERFACE nbi_etienne
     MODULE PROCEDURE etienne_bessel_Itr
     MODULE PROCEDURE etienne_bessel_Irt
     MODULE PROCEDURE etienne_bessel_Ir
     MODULE PROCEDURE etienne_bessel_It
  END INTERFACE



 INTERFACE nbi_david
     MODULE PROCEDURE norm_bessel_Ir
     MODULE PROCEDURE nbit
     MODULE PROCEDURE nbittr
     MODULE PROCEDURE nbitrt
  END INTERFACE

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL
     MODULE PROCEDURE EQUALq
     MODULE PROCEDURE EQUALcq
     MODULE PROCEDURE EQUALq_cq
     MODULE PROCEDURE EQUALcq_q
     MODULE PROCEDURE EQUALqi
     MODULE PROCEDURE EQUALqr
     MODULE PROCEDURE cEQUALqr
     MODULE PROCEDURE cEQUALqi
     !     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
     !     MODULE PROCEDURE AABSEQUAL   ! remove 2002.10.17
     MODULE PROCEDURE DEQUAL  ! added 2002.10.17    ! check2002.10.17
     MODULE PROCEDURE REQUAL   ! added 2002.10.17   ! check2002.10.17
     MODULE PROCEDURE Dequaldacon
     MODULE PROCEDURE equaldacon
     MODULE PROCEDURE Iequaldacon
     ! UNIVERSAL_TAYLOR

     MODULE PROCEDURE fill_uni_r
     MODULE PROCEDURE null_uni
     MODULE PROCEDURE fill_uni  ! new sagan
     MODULE PROCEDURE refill_uni
  end  INTERFACE

  INTERFACE clean
     MODULE PROCEDURE clean_taylor
     MODULE PROCEDURE clean_pbfield
     MODULE PROCEDURE clean_pbresonance
     MODULE PROCEDURE clean_damap
     MODULE PROCEDURE clean_vecfield
     MODULE PROCEDURE clean_vecresonance
     MODULE PROCEDURE clean_onelieexponent
     MODULE PROCEDURE clean_complextaylor
     MODULE PROCEDURE clean_gmap
  END INTERFACE

  INTERFACE print_for_bmad_parser
     MODULE PROCEDURE print_for_bmad_parse
  END INTERFACE


  INTERFACE print
     MODULE PROCEDURE printunitaylor
     MODULE PROCEDURE printq
     MODULE PROCEDURE cprintq
  END INTERFACE


  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryADD  !@2 This is a unary operation
     MODULE PROCEDURE add
     MODULE PROCEDURE unaryADDq  !@2 This is a unary operation
     MODULE PROCEDURE cunaryADDq
     MODULE PROCEDURE addq
     MODULE PROCEDURE caddq
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
     MODULE PROCEDURE unarySUBq
     MODULE PROCEDURE cunarySUBq
     MODULE PROCEDURE subq
     MODULE PROCEDURE csubq
     MODULE PROCEDURE dsubsc
     MODULE PROCEDURE dscsub
     MODULE PROCEDURE subsc
     MODULE PROCEDURE scsub
     MODULE PROCEDURE isubsc
     MODULE PROCEDURE iscsub
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul
     MODULE PROCEDURE mulq
     MODULE PROCEDURE cmulq
     MODULE PROCEDURE dmulsc
     MODULE PROCEDURE dscmul
     MODULE PROCEDURE mulsc
     MODULE PROCEDURE scmul
     MODULE PROCEDURE imulsc
     MODULE PROCEDURE iscmul
     MODULE PROCEDURE rcmulq
     MODULE PROCEDURE cmulqr
     MODULE PROCEDURE ccmulq
     MODULE PROCEDURE cmulqc
  END INTERFACE



  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div
     MODULE PROCEDURE divq
     MODULE PROCEDURE cdivq
     MODULE PROCEDURE ddivsc
     MODULE PROCEDURE dscdiv
     MODULE PROCEDURE divsc
     MODULE PROCEDURE scdiv
     MODULE PROCEDURE idivsc
     MODULE PROCEDURE iscdiv
  END INTERFACE


  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW
     MODULE PROCEDURE POWq
     MODULE PROCEDURE cPOWq
     MODULE PROCEDURE POWR
     MODULE PROCEDURE POWR8
  END INTERFACE


  ! New Operators

  INTERFACE OPERATOR (.mono.)
     MODULE PROCEDURE dputint0   !@1 &nbsp; single integer
     MODULE PROCEDURE dputint    !@1 &nbsp; Accepts J(nv)
     MODULE PROCEDURE dputchar   !@1 &nbsp; Accepts String such as '12'
  END INTERFACE

  INTERFACE OPERATOR (.var.)
     MODULE PROCEDURE varf        !@1 &nbsp; replaces var (overloads DAVAR)
     MODULE PROCEDURE varf001       !@1 replaces var001
  END INTERFACE

  INTERFACE OPERATOR (.d.)
     MODULE PROCEDURE getdiff    !@1 takes derivatives
  END INTERFACE

  INTERFACE OPERATOR (.i.)
     MODULE PROCEDURE GETINTegrate    !@1 takes anti-derivatives
  END INTERFACE


  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDER
     MODULE PROCEDURE getchar
     MODULE PROCEDURE GETint
  END INTERFACE

  INTERFACE OPERATOR (.PAR.)
     MODULE PROCEDURE getcharnd2
     MODULE PROCEDURE GETintnd2
  END INTERFACE

  INTERFACE OPERATOR (.part.)
     MODULE PROCEDURE GETintnd2t
  END INTERFACE


  INTERFACE OPERATOR (<=)
     MODULE PROCEDURE getcharnd2s
     MODULE PROCEDURE GETintnd2s
     MODULE PROCEDURE GETintk
  END INTERFACE

  INTERFACE OPERATOR (.CUT.)
     MODULE PROCEDURE CUTORDER
  END INTERFACE

  INTERFACE OPERATOR (.K.)
     MODULE PROCEDURE getdATRA    ! Used internally primarily
  END INTERFACE

  INTERFACE OPERATOR (.pb.)
     MODULE PROCEDURE pbbra
  END INTERFACE

  ! intrisic functions overloaded

  INTERFACE abs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
     MODULE PROCEDURE absq
     MODULE PROCEDURE cabsq
  END INTERFACE


  INTERFACE abs_square
     MODULE PROCEDURE absq2
     MODULE PROCEDURE cabsq2
  END INTERFACE

  INTERFACE dabs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
  END INTERFACE

  INTERFACE exp
     MODULE PROCEDURE dexpt
     MODULE PROCEDURE c_exp_quaternion
  END INTERFACE

  INTERFACE dexp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE cexp
     MODULE PROCEDURE dexpt
  END INTERFACE
  INTERFACE cdexp
     MODULE PROCEDURE dexpt
  END INTERFACE

  INTERFACE cos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE cdcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE dcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE ccos
     MODULE PROCEDURE dcost
  END INTERFACE

  INTERFACE cosH
     MODULE PROCEDURE dcosHt
  END INTERFACE
  INTERFACE dcosH
     MODULE PROCEDURE dcosHt
  END INTERFACE

  INTERFACE sin
     MODULE PROCEDURE dsint
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

  INTERFACE sinH
     MODULE PROCEDURE dsinHt
  END INTERFACE
  INTERFACE dsinH
     MODULE PROCEDURE dsinHt
  END INTERFACE

  INTERFACE log
     MODULE PROCEDURE dlogt
  END INTERFACE
  INTERFACE dlog
     MODULE PROCEDURE dlogt
  END INTERFACE
  INTERFACE cdlog
     MODULE PROCEDURE dlogt
  END INTERFACE
  INTERFACE clog
     MODULE PROCEDURE dlogt
  END INTERFACE



  INTERFACE sqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE
  INTERFACE dsqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE

  INTERFACE atanh
     MODULE PROCEDURE datanht
  END INTERFACE



  INTERFACE tan
     MODULE PROCEDURE dtant
  END INTERFACE

  INTERFACE tanh
     MODULE PROCEDURE dtanht
  END INTERFACE


  INTERFACE dtan
     MODULE PROCEDURE dtant
  END INTERFACE

  ! Non-intrisic Functions

  INTERFACE pek
     MODULE PROCEDURE pek000 ! not private
  END INTERFACE

  INTERFACE pok
     MODULE PROCEDURE pok000  ! not private
  END INTERFACE

  INTERFACE shiftda
     MODULE PROCEDURE shift000  ! not private
  END INTERFACE

  !  INTERFACE var
  !     MODULE PROCEDURE var000  ! not private
  !     MODULE PROCEDURE var001  ! not private
  !  END INTERFACE

  INTERFACE cfu
     MODULE PROCEDURE CFU000_new  ! not private
  END INTERFACE

  INTERFACE full_abs
     MODULE PROCEDURE full_absT
  END INTERFACE

  !  INTERFACE daread
  !     MODULE PROCEDURE rea
  !  END INTERFACE

  !  INTERFACE read
  !     MODULE PROCEDURE rea
  !  END INTERFACE

  !  INTERFACE daprint
  !     MODULE PROCEDURE pri
  !  END INTERFACE



  !  INTERFACE print
  !     MODULE PROCEDURE pri
  !  END INTERFACE


  ! Constructors and Destructors

  INTERFACE alloc
     MODULE PROCEDURE allocda
     MODULE PROCEDURE A_OPT
     MODULE PROCEDURE allocdas
     MODULE PROCEDURE alloc_u
  END INTERFACE

  INTERFACE KILL
     MODULE PROCEDURE KILLda
     MODULE PROCEDURE KILLdas
     MODULE PROCEDURE K_opt
     MODULE PROCEDURE kill_uni
  END INTERFACE

  INTERFACE alloctpsa
     MODULE PROCEDURE allocda
  END INTERFACE

  INTERFACE KILLtpsa
     MODULE PROCEDURE KILLda
  END INTERFACE


  ! management routines

  INTERFACE ass
     MODULE PROCEDURE asstaylor   !2000.12.25
  END INTERFACE



CONTAINS




  SUBROUTINE  change_default_tpsa(i)
    implicit none
    INTEGER, intent(in) :: I
    if(last_tpsa==0) then
       if(i==1) then
          default_tpsa=.true.
          if(i==1.and.lingyun_yang )write(6,*) " Default TPSA is CPP package of Yang"
          call change_package(i)
       else
          default_tpsa=.false.
          call change_package(i)
          if(i==2.and.(.not.lingyun_yang) )write(6,*) " Default TPSA is FORTRAN package of Berz (LBNL)"
       endif
    else
       write(6,*) " You could not change default TPSA here "
       write(6,*) " Only prior to any call to TPSA or PTC or after a PTC_END "
       stop 666
    endif
  end   SUBROUTINE  change_default_tpsa


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

  subroutine count_taylor(n,ns,ne)
    implicit none
    integer n,ns,ne,i

    call count_da(n)
    ns=0
    do i=1,ndumt
       ns=scratchda(i)%n+ns
    enddo
    ne=n-ns
  end subroutine count_taylor




  FUNCTION unaryADDq( S1 )
    implicit none
    TYPE (quaternion) unaryADDq
    TYPE (quaternion), INTENT (IN) :: S1





    unaryADDq=s1

  END FUNCTION unaryADDq

  FUNCTION unarySUBq( S1 )
    implicit none
    TYPE (quaternion) unarySUBq
    TYPE (quaternion), INTENT (IN) :: S1


         unarySUBq%x= -s1%x

  END FUNCTION unarySUBq


  FUNCTION invq( S1 )
    implicit none
    TYPE (quaternion) invq
    TYPE (quaternion), INTENT (IN) :: S1
    real(dp) norm
     integer i

              invq=s1
              do i=1,3
                invq%x(i)=-invq%x(i)
              enddo
                norm=abs_square(invq)
              do i=0,3
                invq%x(i)=invq%x(i)/norm
              enddo

  END FUNCTION invq


  FUNCTION absq( S1 )
    implicit none
    real(dp) absq
    TYPE (quaternion), INTENT (IN) :: S1
    integer i



     absq=sqrt(abs_square(s1))
  END FUNCTION absq


  FUNCTION absq2( S1 )
    implicit none
    real(dp) absq2
    TYPE (quaternion), INTENT (IN) :: S1
    integer i


           absq2=0
       do i=0,3
         absq2 = s1%x(i)**2+absq2
       enddo
  END FUNCTION absq2


! complex quaternion

  FUNCTION cunaryADDq( S1 )
    implicit none
    TYPE (complex_quaternion) cunaryADDq
    TYPE (complex_quaternion), INTENT (IN) :: S1




    cunaryADDq=s1

  END FUNCTION cunaryADDq

  FUNCTION cunarySUBq( S1 )
    implicit none
    TYPE (complex_quaternion) cunarySUBq
    TYPE (complex_quaternion), INTENT (IN) :: S1


         cunarySUBq%x= -s1%x

  END FUNCTION cunarySUBq


  FUNCTION cinvq( S1 )
    implicit none
    TYPE (complex_quaternion) cinvq
    TYPE (complex_quaternion), INTENT (IN) :: S1
    real(dp) norm
     integer i

              cinvq=s1
              do i=1,3
                cinvq%x(i)=-cinvq%x(i)
              enddo
                norm=abs_square(cinvq)
              do i=0,3
                cinvq%x(i)=cinvq%x(i)/norm
              enddo

  END FUNCTION cinvq


  FUNCTION cabsq( S1 )
    implicit none
    real(dp) cabsq
    TYPE (complex_quaternion), INTENT (IN) :: S1
    integer i



     cabsq=sqrt(abs_square(s1))
  END FUNCTION cabsq


  FUNCTION cabsq2( S1 )
    implicit none
    complex(dp) cabsq2
    TYPE (complex_quaternion), INTENT (IN) :: S1
    integer i


           cabsq2=0
       do i=0,3
         cabsq2 = abs(s1%x(i))**2+cabsq2
       enddo
  END FUNCTION cabsq2


 subroutine log_complex_quaternion(qn0,logqn0,epso)
!#restricted: normal
!# This routine takes the log assuming no orbital map
!#
  implicit none

  type(complex_quaternion), intent (inout) :: qn0,logqn0
  type(complex_quaternion)  q1,qs,dr
  complex(dp) s
  real(dp) sn,normb,norma,eps
  real(dp), optional :: epso
  complex(dp) cn
  logical diff
  integer i,k

  eps=1.d-3
  if(present(epso)) eps=epso
  sn=sqrt(qn0%x(1)**2+qn0%x(2)**2+qn0%x(3)**2)

if(sn>eps) then
  q1=qn0
  q1%x(0)=0.0_dp
  qs=0.0_dp
  s=sqrt(q1%x(1)**2+q1%x(2)**2+q1%x(3)**2)
  qs%x(0)=1.0_dp/s
  logqn0=q1*qs   ! q1=n
  s= sqrt(1.0_dp-s**2) + i_* s
  s=-i_*log(s)

  logqn0%x(0)=0.0_dp
  do i=1,3
   logqn0%x(i)=s*logqn0%x(i)
  enddo
goto 1
else
  diff=.false.
!write(6,*) "small quaternion "

 q1=qn0
 q1%x(0)=q1%x(0)-1.0_dp
 logqn0=0.0_dp
 qs=q1
normb=1.d38
 do i=1,100
  cn=-(-1.0_dp)**i/i
dr=logqn0
  logqn0=logqn0+cn*qs
dr=dr-logqn0
  qs=qs*q1

! call c_full_norm_quaternion(dr,k,norma)
norma=abs(dr)
    if(diff) then
      if(normb>=norma) goto 1
      normb=norma
    else
     if(norma<1.d-10) then
      diff=.true.
      normb=norma
     endif
    endif
 enddo
endif

write(6,*) " no convergence  in  log_quaternion"

1 continue


 end  subroutine log_complex_quaternion

  function c_exp_quaternion(h_axis,ds) ! spin routine
    implicit none
    TYPE(complex_quaternion) c_exp_quaternion
    TYPE(complex_quaternion),optional, INTENT(INout) :: DS
    TYPE(complex_quaternion), INTENT(IN) :: h_axis
    integer  nmax
    integer i,localmaster,k
    TYPE(complex_quaternion) dh,dhn,dr,dst
    real(dp) eps,norm1,norm2
    complex(dp) c
    logical check




    check=.true.
    eps=1.d-5
    nmax=1000

     c_exp_quaternion=1.0_dp

    dh=h_axis


    dhn=1.0_dp
    c=1.0_dp
    norm1=mybig
    do i=1,nmax
       dhn=dhn*dh
       c=1.0_dp/i
       dhn=c*dhn

       dr=c_exp_quaternion

       c_exp_quaternion=c_exp_quaternion+dhn

       dr=c_exp_quaternion+(-1.0_dp,0.0_dp)*dr

       norm2=abs(dr)


       if(check) then
          if(norm2<eps.and.i>10) then
             check=.false.
          endif
       else
          if(norm2>=norm1) exit
       endif
       norm1=norm2
    enddo

    if(i>nmax-10) then
       write(6,*) "no convergence in c_exp_quaternion, enter 0 to stop "
       read(5,*) norm1
       if(norm1==0)  stop 1066
    endif
    if(present(ds)) c_exp_quaternion=c_exp_quaternion*ds


  end   function c_exp_quaternion



!!!!!!!!!!!!!!!!!!!!!

  FUNCTION unaryADD( S1 )
    implicit none
    TYPE (TAYLOR) unaryADD
    TYPE (TAYLOR), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_%STABLE_DA) then
     unaryADD%i=0
     RETURN
    endif

    localmaster=master

    !    call check(s1)
    call ass(unaryADD)

    unaryADD=s1

    master=localmaster

  END FUNCTION unaryADD

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (TAYLOR) unarySUB
    TYPE (TAYLOR), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_%STABLE_DA) then
     unarysub%i=0
     RETURN
    endif
    localmaster=master

    call ass(unarySUB)
    if(old_package) then
     call dacmu(s1%i,-1.0_dp,temp%i)
     call dacop(temp%i,unarySUB%i)
     else
! GTPSA REMOVED !      call mad_tpsa_scl(s1%j,-1.0_dp,unarySUB%j)
     endif
    master=localmaster

  END FUNCTION unarySUB


  SUBROUTINE  maketree(S1,s2)
    implicit none
    type (TAYLOR),INTENT(IN)::S1
    type (TAYLOR),INTENT(inOUT):: s2
    IF(.NOT.C_%STABLE_DA) RETURN

        if(old_package) then
    call mtree((/s1%i/),1,(/s2%i/),1)
       else
          s2=s1
       endif
  END SUBROUTINE maketree

  SUBROUTINE  allocda(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1

    !    IF(first_time) THEN
    IF(last_tpsa==0.and.old_package) THEN

       write(6,*) " No TPSA package ever initialized "


    ENDIF
        if(old_package) then
    s1%i=0
    call etall1(s1%i)
        else
! GTPSA REMOVED !  s1%j=mad_tpsa_newd(d_berz , mad_tpsa_default)
        !   call nullnewda(s1%j)
        !   call allocnewda(s1%j)
        endif
  END SUBROUTINE allocda

  SUBROUTINE  A_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (taylor),INTENT(INout)::S1,S2
    type (taylor),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call allocda(s1)
    call allocda(s2)
    if(present(s3)) call allocda(s3)
    if(present(s4)) call allocda(s4)
    if(present(s5)) call allocda(s5)
    if(present(s6)) call allocda(s6)
    if(present(s7)) call allocda(s7)
    if(present(s8)) call allocda(s8)
    if(present(s9)) call allocda(s9)
    if(present(s10))call allocda(s10)
  END SUBROUTINE A_opt

  SUBROUTINE  K_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (taylor),INTENT(INout)::S1,S2
    type (taylor),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
    call KILLDA(s1)
    call KILLDA(s2)
    if(present(s3)) call KILLDA(s3)
    if(present(s4)) call KILLDA(s4)
    if(present(s5)) call KILLDA(s5)
    if(present(s6)) call KILLDA(s6)
    if(present(s7)) call KILLDA(s7)
    if(present(s8)) call KILLDA(s8)
    if(present(s9)) call KILLDA(s9)
    if(present(s10))call KILLDA(s10)
  END SUBROUTINE K_opt


  SUBROUTINE  ALLOCDAS(S1,k)
    implicit none
    type (TAYLOR),INTENT(INOUT),dimension(:)::S1
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S1,DIM=1)
       N=LBOUND(S1,DIM=1)+K-1
    else
       I=LBOUND(S1,DIM=1)
       N=UBOUND(S1,DIM=1)
    endif

    DO   J=I,N
       CALL allocDA(S1(j))
    ENDDO

  END SUBROUTINE ALLOCDAS

  SUBROUTINE  KILLda(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1
        if(old_package) then
    call DADAL1(s1%i)
        else
! GTPSA REMOVED !             call mad_tpsa_del  (S1%j);
! GTPSA REMOVED !              S1%j=c_null
        endif

  END SUBROUTINE KILLda

  SUBROUTINE  KILLDAS(S1,k)
    implicit none
    type (TAYLOR),INTENT(INOUT),dimension(:)::S1
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S1,DIM=1)
       N=LBOUND(S1,DIM=1)+K-1
    else
       I=LBOUND(S1,DIM=1)
       N=UBOUND(S1,DIM=1)
    endif

    DO   J=I,N
       CALL KILLDA(S1(j))
    ENDDO

  END SUBROUTINE KILLDAS


  SUBROUTINE  EQUAL(S2,S1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    type (TAYLOR),INTENT(IN)::S1
    IF(.NOT.C_%STABLE_DA) RETURN

    call check_snake
    !    if(old) then
    if(s2%i==0.and.old_package) then
       call crap1("EQUAL 1 in tpsa") !call allocw(s2)
! GTPSA REMOVED !    elseif(.not.c_associated(s2%j).and.(.not.(old_package))) then
       call crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
     if(old_package) then
      CALL DACOP(S1%I,S2%I)
     else
! GTPSA REMOVED !      call mad_tpsa_copy(S1%j,S2%j)
     endif
    !    else
    !      IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUAL 3") !call allocw(s2)
    !      IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("EQUAL 4") !call allocw(s1)
    !      call newdacop(S1%j,S2%j)
    !   endif
  END SUBROUTINE EQUAL

  SUBROUTINE  EQUALq(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    type (quaternion),INTENT(IN)::S1
    integer i


    do i=0,3
    s2%x(i)=s1%x(i)
    enddo

  end SUBROUTINE  EQUALq


  SUBROUTINE  EQUALcq(S2,S1)
    implicit none
    type (complex_quaternion),INTENT(inOUT)::S2
    type (complex_quaternion),INTENT(IN)::S1
    integer i


    do i=0,3
    s2%x(i)=s1%x(i)
    enddo

  end SUBROUTINE  EQUALcq


  SUBROUTINE  EQUALcq_q(S2,S1)
    implicit none
    type (complex_quaternion),INTENT(inOUT)::S2
    type (quaternion),INTENT(IN)::S1
    integer i

    do i=0,3
    s2%x(i)=s1%x(i)
    enddo

  end SUBROUTINE  EQUALcq_q

  SUBROUTINE  EQUALq_cq(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    type (complex_quaternion),INTENT(IN)::S1
    integer i


    do i=0,3
    s2%x(i)=s1%x(i)
    enddo

  end SUBROUTINE  EQUALq_cq

  SUBROUTINE  EQUALqr(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    integer i


    do i=0,3
    s2%x(i)=0
    enddo
    s2%x(0)=s1
  end SUBROUTINE  EQUALqr


  SUBROUTINE  cEQUALqr(S2,S1)
    implicit none
    type (complex_quaternion),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    integer i


    do i=0,3
    s2%x(i)=0
    enddo
    s2%x(0)=s1
  end SUBROUTINE  cEQUALqr

  SUBROUTINE  EQUALqi(S2,S1)
    implicit none
    type (quaternion),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i


    do i=0,3
    s2%x(i)=0
    enddo
    s2%x(s1)=1
  end SUBROUTINE  EQUALqi


  SUBROUTINE  cEQUALqi(S2,S1)
    implicit none
    type (complex_quaternion),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i

    do i=0,3
    s2%x(i)=0
    enddo
    s2%x(s1)=1
  end SUBROUTINE  cEQUALqi

  SUBROUTINE  DEQUAL(R1,S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    real(dp), INTENT(inOUT)::R1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake
    if(old_package ) then
    R1=S2.SUB.'0'
    else
! GTPSA REMOVED !     r1=mad_tpsa_get0(s2%j)
    endif
  END SUBROUTINE DEQUAL

  SUBROUTINE  REQUAL(R1,S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    REAL(SP), INTENT(inOUT)::R1
    IF(.NOT.C_%STABLE_DA) RETURN

    if(real_warning) call real_stop
    call check_snake

    if(old_package ) then
    R1=S2.SUB.'0'
    else
! GTPSA REMOVED !     r1=mad_tpsa_get0(s2%j)
    endif

  END SUBROUTINE REQUAL

  function  DAABSEQUAL(S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    real(dp) DAABSEQUAL
    DAABSEQUAL=0
    IF(.NOT.C_%STABLE_DA) RETURN

    if(old_package ) then
    DAABSEQUAL=abs(S2.sub.'0')
    else
! GTPSA REMOVED !     DAABSEQUAL=abs(mad_tpsa_get0(s2%j))
    endif


  END function DAABSEQUAL


  SUBROUTINE  DEQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    real(dp), INTENT(IN)::R1
    IF(.NOT.C_%STABLE_DA) RETURN

    !    if(old) then
    if(s2%i==0.and.old_package) then
       call crap1("EQUAL 1 in DEQUALDACON") !call allocw(s2)
! GTPSA REMOVED !    elseif(.not.c_associated(s2%j).and.(.not.(old_package))) then
       call crap1("EQUAL 1 in DEQUALDACON") !call allocw(s2)
    endif
    if(old_package ) then
     CALL DACON(S2%I,R1)
    else
! GTPSA REMOVED !     call  mad_tpsa_setvar(S2%j,r1,0,0.0_dp)
    endif

    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("DEQUALDACON 2") !call allocw(s2)
    !       CALL newDACON(S2%j,R1)
    !    endif
  END SUBROUTINE DEQUALDACON

  SUBROUTINE  EQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    REAL(SP), INTENT(IN)::R1
    real(dp) R2
    IF(.NOT.C_%STABLE_DA) RETURN
    if(real_warning) call real_stop
    call check_snake

    if(real_warning) call real_stop
    !    if(old) then
  !  if(s2%i==0) call crap1("EQUALDACON 1") !call allocw(s2)
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUALDACON 2") !call allocw(s2)
    !    endif

     r2=REAL(r1,kind=DP)
     s2=r2

  END SUBROUTINE EQUALDACON

  SUBROUTINE  IEQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    INTEGER, INTENT(IN)::R1
    real(dp) r2
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake


    ! if(old) then
  !  if(s2%i==0) call crap1("IEQUALDACON 1") !call allocw(s2)
    !    else
    !       IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("IEQUALDACON 2") !call allocw(s2)
    !    endif
    r2=REAL(r1,kind=DP)
    s2=r2


  END SUBROUTINE IEQUALDACON

  FUNCTION dexpt( S1 )
    implicit none
    TYPE (taylor) dexpt
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
     dexpt%i=0
     RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(dexpt)

     if(old_package) then
    call dafun('EXP ',s1%i,temp%i)
    call dacop(temp%i,dexpt%i)
        else
! GTPSA REMOVED !          call  mad_tpsa_exp(s1%j,dexpt%j)
        endif

    master=localmaster

  END FUNCTION dexpt

  FUNCTION FULL_ABST( S1 )
    implicit none
    real(dp) FULL_ABST
    TYPE (taylor), INTENT (IN) :: S1
    FULL_ABST=0
    IF(.NOT.C_%STABLE_DA) RETURN
    !    call check(s1)

     if(old_package ) then
      CALL DAABS(S1%I,FULL_ABST)
        else
! GTPSA REMOVED !           FULL_ABST= mad_tpsa_nrm1(S1%j,c_null)
        endif

  END FUNCTION FULL_ABST




  FUNCTION dtant( S1 )
    implicit none
    TYPE (taylor) dtant
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster

    IF(.NOT.C_%STABLE_DA) then
      dtant%i=0
     RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(dtant)

     if(old_package) then
    call dafun('SIN ',s1%i,temp%i)
    call dacop(temp%i,dtant%i)
    call dafun('COS ',s1%i,temp%i)
    call dadiv(dtant%i,temp%i,dtant%i)
        else
! GTPSA REMOVED !    call mad_tpsa_tan(s1%j,dtant%j)
        endif

    master=localmaster

  END FUNCTION dtant

  FUNCTION datanht( S1 )
    implicit none
    TYPE (taylor) datanht
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      datanht%i=0
     RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(datanht)

    if(old_package ) then
     datanht=log((1+s1)/sqrt(1-s1))/2.0_dp
    else
! GTPSA REMOVED !     call mad_tpsa_atanh(S1%j,datanht%j)
    endif
    master=localmaster

  END FUNCTION datanht

  FUNCTION dcost( S1 )
    implicit none
    TYPE (taylor) dcost
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      dcost%i=0
     RETURN
    endif
    localmaster=master



    !    call check(s1)
    call ass(dcost)

     if(old_package) then
    call dafun('COS ',s1%i,temp%i)
    call dacop(temp%i,dcost%i)
        else
! GTPSA REMOVED !    call mad_tpsa_cos(s1%j,dcost%j)
        endif

    master=localmaster

  END FUNCTION dcost

  FUNCTION dsint( S1 )
    implicit none
    TYPE (taylor) dsint
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      dsint%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dsint)
     if(old_package) then
    call dafun('SIN ',s1%i,temp%i)
    call dacop(temp%i,dsint%i)
        else
! GTPSA REMOVED !         call  mad_tpsa_sin(s1%j,dsint%j)
        endif

    master=localmaster

  END FUNCTION dsint

  FUNCTION dsinHt( S1 )
    implicit none
    TYPE (taylor) dsinHt
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      dsinHt%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dsinHt)
     if(old_package) then
    call dafun('SINH',s1%i,temp%i)
    call dacop(temp%i,dsinHt%i)
        else
! GTPSA REMOVED !           call mad_tpsa_sinh(s1%j,dsinHt%j)
        endif
    master=localmaster

  END FUNCTION dsinHt

  FUNCTION DCOSHT( S1 )
    implicit none
    TYPE (taylor) DCOSHT
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      DCOSHT%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(DCOSHT)
     if(old_package) then
    call dafun('COSH',s1%i,temp%i)
    call dacop(temp%i,DCOSHT%i)
        else
! GTPSA REMOVED !           call mad_tpsa_cosh(s1%j,DCOSHT%j)
        endif

    master=localmaster

  END FUNCTION DCOSHT

  FUNCTION dtanht( S1 )
    implicit none
    TYPE (taylor) dtanht
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      dtanht%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dtanht)
     if(old_package) then

    dtanht=sinh(s1)/cosh(s1)
        else
! GTPSA REMOVED !           call mad_tpsa_tanh(s1%j,dtanht%j)
        endif

    master=localmaster

  END FUNCTION dtanht

  FUNCTION dlogt( S1 )
    implicit none
    TYPE (taylor) dlogt
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      dlogt%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dlogt)
     if(old_package) then
    call dafun('LOG ',s1%i,temp%i)
    call dacop(temp%i,dlogt%i)
        else
! GTPSA REMOVED !           call mad_tpsa_log(s1%j,dlogt%j)
        endif

    master=localmaster

  END FUNCTION dlogt

  FUNCTION dsqrtt( S1 )
    implicit none
    TYPE (taylor) dsqrtt
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      dsqrtt%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dsqrtt)

     if(old_package) then
    call dafun('SQRT',s1%i,temp%i)
    call dacop(temp%i,dsqrtt%i)
        else
! GTPSA REMOVED !           call mad_tpsa_sqrt(s1%j,dsqrtt%j)
        endif

    master=localmaster

  END FUNCTION dsqrtt

  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (taylor) mul
    TYPE (taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      mul%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    !    call check(s2)
    call ass(mul)

     if(old_package) then
    call damul(s1%i,s2%i,temp%i)
    call dacop(temp%i,mul%i)
        else
! GTPSA REMOVED !           call mad_tpsa_mul(s1%j,s2%j,mul%j)
        endif

    master=localmaster

  END FUNCTION mul

  FUNCTION pbbra( S1, S2 )
    implicit none
    TYPE (taylor) pbbra
    TYPE (taylor), INTENT (IN) :: S1, S2
    integer localmaster
    integer i
    IF(.NOT.C_%STABLE_DA) then
      pbbra%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    !    call check(s2)
    call ass(pbbra)

     if(old_package ) then
    pbbra=0.0_dp
    do i=1,nd
       pbbra=(s1.d.(2*i-1))*(s2.d.(2*i))-(s2.d.(2*i-1))*(s1.d.(2*i))+pbbra
    enddo

        else
! GTPSA REMOVED !          call mad_tpsa_poisson(s1%j,s2%j,pbbra%j,nd2)
        endif

    master=localmaster

  END FUNCTION pbbra

  FUNCTION GETORDER( S1, S2 )
    implicit none
    TYPE (taylor) GETORDER
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      GETORDER%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(GETORDER)

     if(old_package ) then
    CALL TAKE(S1%I,S2,GETORDER%i)
 !   call dacop(temp%i,GETORDER%i)
        else
! GTPSA REMOVED !       call  mad_tpsa_cutord(S1%j,temp%j,s2+1)
! GTPSA REMOVED !       call  mad_tpsa_cutord(S1%j,GETORDER%j,s2)
! GTPSA REMOVED !       call mad_tpsa_sub(temp%j,GETORDER%j,GETORDER%j)
        endif
    master=localmaster

  END FUNCTION GETORDER




  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (taylor) CUTORDER
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,s
    IF(.NOT.C_%STABLE_DA) then
       CUTORDER%i=0
      RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(CUTORDER)

     if(old_package ) then
       if(s2>0) then
      call datrunc(S1%I,s2,CUTORDER%i)
      else
       s=-s2+1
       call datrunc(S1%I,s,CUTORDER%i)
        call dasub(S1%I,CUTORDER%i,CUTORDER%i)
      endif
    else
! GTPSA REMOVED !       call  mad_tpsa_cutord(S1%j,CUTORDER%j,s2)
   endif

    master=localmaster

  END FUNCTION CUTORDER

  FUNCTION dputchar( S1, S2 )
    implicit none
    TYPE (taylor) dputchar
    real(dp), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dputchar%i=0
      RETURN
    endif

    localmaster=master


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
   !    if(i>nv) then
   !       if(j(i)>0) then
   !          dputchar=0.0_dp
   !          !             call var(dputchar,zero,0)
   ! master=localmaster
   !          return
   !       endif
!       endif
    enddo



    dputchar=0.0_dp
    !    call var(dputchar,zero,0)
    CALL pok(dputchar,j,s1)

    master=localmaster

  END FUNCTION dputchar

  FUNCTION dputint( S1, S2 )
    implicit none
    TYPE (taylor) dputint
    real(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    integer j(lnv),i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dputint%i=0
      RETURN
    endif


    localmaster=master


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
             !             call var(dputint,zero,0)
             dputint=0.0_dp
    master=localmaster
             return
          endif
       endif
    enddo


    dputint=0.0_dp
    !    call var(dputint,zero,0)
    CALL pok(dputint,j,s1)

    master=localmaster

  END FUNCTION dputint

  FUNCTION dputint0( S1, S2 )
    implicit none
    TYPE (taylor) dputint0
    real(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer j(lnv)
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dputint0%i=0
      RETURN
    endif
if(check_gtpsa(i=s2)) then
 else
    localmaster=master
    call ass(dputint0)
   dputint0=S1
    master=localmaster
       return
!  write(6,*) " check_gtpsa failed dputint0 "
!  stop 999
 endif

    localmaster=master


    call ass(dputint0)

    if(s2>nv) then
       dputint0=S1
    master=localmaster
       return
    endif

if(old_package ) then
    j=0



    dputint0=0.0_dp
    !    call var(dputint0,zero,s2)

    j(s2)=1
    CALL pok(dputint0,j,s1)
else

! GTPSA REMOVED ! call mad_tpsa_clear(dputint0%j)
!  bug laurent
! GTPSA REMOVED !  call  mad_tpsa_setsm(dputint0%j,2,[s2,1],0.0_dp,s1)

endif
    master=localmaster

  END FUNCTION dputint0


  FUNCTION GETCHARnd2s( S1, S2 )
    implicit none
    TYPE (taylor) GETCHARnd2s
    TYPE (taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2

    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETCHARnd2s%i=0
      RETURN
    endif

    localmaster=master


    call ass(GETCHARnd2s)


    GETCHARnd2s=s1.par.s2
    call  shiftda(GETCHARnd2s,GETCHARnd2s, len(trim(ADJUSTR (s2) )))

    master=localmaster


  END FUNCTION GETCHARnd2s

  FUNCTION GETintnd2s( S1, S2 )
    implicit none
    TYPE (taylor) GETintnd2s
    TYPE (taylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)

    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETintnd2s%i=0
      RETURN
    endif
    localmaster=master


    call ass(GETintnd2s)

    GETintnd2s=s1.par.s2

    call  shiftda(GETintnd2s,GETintnd2s, size(s2) )

    master=localmaster


  END FUNCTION GETintnd2s

  FUNCTION GETintk( S1, S2 )
    implicit none
    TYPE (taylor) GETintk
    TYPE (taylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2

    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETintk%i=0
      RETURN
    endif
    localmaster=master


    call ass(GETintk)


    call  shiftda(s1,GETintk, s2 )

    master=localmaster


  END FUNCTION GETintk



  FUNCTION GETchar( S1, S2 )

    implicit none
    real(dp) GETchar,r1
    TYPE (taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer j(lnv),i,c

    getchar=0
    IF(.NOT.C_%STABLE_DA) RETURN


    resul = s2
    call context(resul)

    do i=1,lnv
       j(i)=0
    enddo


    nd2par= len_trim(resul)





    do i=1,nd2par
       CALL  CHARINT(RESUL(I:I),J(I))
    enddo

c=0
    do i=c_%nv+1,lnv
       c=j(i)+c
    enddo





if(c>0) then
r1=0.0_dp
else
    CALL pek(S1,j,r1)
endif

     GETchar=r1

  END FUNCTION GETchar




  FUNCTION GETint( S1, S2 )
    implicit none
    real(dp) GETint,r1
    TYPE (taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2(:)
    integer j(lnv),i,c

     getint=0
    IF(.NOT.C_%STABLE_DA) RETURN



    do i=1,lnv
       j(i)=0
    enddo


    nd2par= size(s2)

    do i=1,nd2par
       J(I)=s2(i)
    enddo

    c=0
    do i=c_%nv+1,lnv
       c=j(i)+c
    enddo

if(c>0) then
r1=0.0_dp
else
    CALL pek(S1,j,r1)
endif

    GETint=r1

  END FUNCTION GETint




  FUNCTION GETdiff( S1, S2 )
    implicit none
    TYPE (taylor) GETdiff
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster

    IF(.NOT.C_%STABLE_DA) then
       GETdiff%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(GETdiff)

     if(old_package ) then
    CALL dader(S2,S1%I,temp%i)
    call dacop(temp%i,GETdiff%i)
         else
! GTPSA REMOVED !      call mad_tpsa_deriv(S1%j,GETdiff%j,s2)
        endif
    master=localmaster

  END FUNCTION GETdiff

  FUNCTION GETINTegrate( S1, S2 )
    implicit none
    TYPE (taylor) GETINTegrate
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,n,i
    type(taylor) t,x
    real(dp) value
    integer, allocatable :: jc(:)
    IF(.NOT.C_%STABLE_DA) then
       GETINTegrate%i=0
      RETURN
    endif
    localmaster=master

!if(old_package ) then
    allocate(jc(c_%nv))
    jc=0
    !    call check(s1)
    call ass(GETINTegrate)
    call alloc(t,x)
    t=s1
    x=0
    call taylor_cycle(t,size=n)

    do i=1,n
       call taylor_cycle(t,ii=i,value=value,j=jc)
         x=((value/(jc(s2)+1)).mono.jc)*(1.0_dp.mono.s2)+x

    enddo

    GETINTegrate=x

    call kill(t,x)
    deallocate(jc)
!else
!endif
    master=localmaster

  END FUNCTION GETINTegrate

 function check_gtpsa(i,J)
 implicit none
 logical check_gtpsa
 integer,optional :: i
 integer,optional :: J(:)
 integer ii,k

  check_gtpsa=.true.

 if(present(i)) then
   if(.not.(i>=1.and.i<=c_%nv)) then
    check_gtpsa=.false.
    write(6,*)"i,c_%nv", i,c_%nv
  endif
 else
  k=0
 do ii=1,min(c_%nv,size(j))
  k=k+j(ii)
 enddo
 if(k>c_%no) then
   check_gtpsa=.false.
 !   write(6,*)"k>c_%no", k,c_%no
  endif
 k=0
  do ii=c_%nv+1,size(j)
  k=k+j(ii)
  enddo
  if(k>0) then
   check_gtpsa=.false.
    write(6,*)"=c_%nv+1,size(j) ", j(c_%nv+1:size(j))
  endif
 endif

 end function check_gtpsa

  FUNCTION GETdatra( S1, S2 )
    implicit none
    TYPE (taylor) GETdatra
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster,j(lnv)
     type(taylor) junk
      integer i,mf
    real(dp) value
    IF(.NOT.C_%STABLE_DA) then
       GETdatra%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(GETdatra)

     if(old_package ) then
    CALL datra(S2,S1%I,temp%i)
    call dacop(temp%i,GETdatra%i)
        else
    call alloc(junk)
   junk=0.0_dp
   call taylor_cycle(S1,size=mf)

   do i=1,mf
   call taylor_cycle(s1,ii=i,value=VALUE,j=J)
!   write(6,*) i,value

   if(value/=0.d0) then
    if(j(S2)/=0) then
     j(S2)=j(S2)-1
     junk=junk +  (VALUE.mono.J)
    endif
   endif
  enddo
    GETdatra=junk
    call kill(junk)
        endif
    master=localmaster

  END FUNCTION GETdatra


  FUNCTION POWq( S1, R2 )
    implicit none
    TYPE (quaternion) POWq,temp
    TYPE (quaternion), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster

    temp=1.0_dp

    R22=IABS(R2)
    DO I=1,R22
       temp=temp*s1
    ENDDO
    IF(R2.LT.0) THEN
       temp=invq(temp)
    ENDIF
     powq=temp

  END FUNCTION POWq


  FUNCTION cPOWq( S1, R2 )
    implicit none
    TYPE (complex_quaternion) cPOWq,temp
    TYPE (complex_quaternion), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster

    temp=1.0_dp

    R22=IABS(R2)
    DO I=1,R22
       temp=temp*s1
    ENDDO
    IF(R2.LT.0) THEN
       temp=cinvq(temp)
    ENDIF
     cpowq=temp

  END FUNCTION cPOWq

  FUNCTION POW( S1, R2 )
    implicit none
    TYPE (taylor) POW
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       POW%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(POW)

     if(old_package) then
    CALL DACON(temp%i,1.0_dp)

    R22=IABS(R2)
    DO I=1,R22
       CALL DAMUL(temp%i,S1%I,temp%i)
    ENDDO
    IF(R2.LT.0) THEN
       CALL DADIC(temp%i,1.0_dp,temp%i)
    ENDIF
    call dacop(temp%i,POW%i)
    else
    temp=1.0_dp
 !   CALL DACON(temp%i,1.0_dp)

    R22=IABS(R2)
    DO I=1,R22
! GTPSA REMOVED !      call mad_tpsa_mul(temp%j,s1%j,temp%j)
!       CALL DAMUL(temp%i,S1%I,temp%i)
    ENDDO
    IF(R2.LT.0) THEN
       CALL DADIC(temp%i,1.0_dp,temp%i)
    ENDIF
      pow=temp
 !   call dacop(temp%i,POW%i)
    endif
    master=localmaster
  END FUNCTION POW

  FUNCTION POWR8( S1, R2 )
    implicit none
    TYPE (taylor) POWR8
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: R2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       POWR8%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(POWR8)

    ! if(old) then
    CALL DAFUN('LOG ',S1%I,temp%i)
    CALL DACMU(temp%i,R2,temp%i)
    CALL DAFUN('EXP ',temp%i,temp%i)
    call dacop(temp%i,POWR8%i)
    !    ELSE
    !       CALL NEWDAFUN('LOG ',S1%J,TEMPL)
    !       CALL NEWDACMU(TEMPL,R2,TEMPL)
    !       CALL NEWDAFUN('EXP ',TEMPL,POWR8%J)
    !    endif
    master=localmaster
  END FUNCTION POWR8

  FUNCTION POWR( S1, R2 )
    implicit none
    TYPE (taylor) POWR
    TYPE (taylor), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: R2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       POWR%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(POWR)

    ! if(old) then
    CALL DAFUN('LOG ',S1%I,temp%i)
    CALL DACMU(temp%i,REAL(R2,kind=DP),temp%i)
    CALL DAFUN('EXP ',temp%i,temp%i)
    call dacop(temp%i,POWR%i)
    !    ELSE
    !       CALL NEWDAFUN('LOG ',S1%J,TEMPL)
    !       CALL NEWDACMU(TEMPL,REAL(R2,kind=DP),TEMPL)
    !       CALL NEWDAFUN('EXP ',TEMPL,POWR%J)
    !    endif
    master=localmaster
  END FUNCTION POWR



  FUNCTION dmulsc( S1, sc )
    implicit none
    TYPE (taylor) dmulsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dmulsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dmulsc)

     if(old_package) then
    call dacmu(s1%i,sc,temp%i)
    call dacop(temp%i,dmulsc%i)
        else
! GTPSA REMOVED !          call mad_tpsa_scl(s1%j,sc,dmulsc%j)
         endif

    master=localmaster
  END FUNCTION dmulsc

  FUNCTION mulsc( S1, sc )
    implicit none
    TYPE (taylor) mulsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       mulsc%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(mulsc)

     if(old_package) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,mulsc%i)
        else
! GTPSA REMOVED !              call mad_tpsa_scl(s1%j,REAL(sc,kind=DP),mulsc%j)
        endif
    master=localmaster
  END FUNCTION mulsc

  FUNCTION imulsc( S1, sc )
    implicit none
    TYPE (taylor) imulsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       imulsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(imulsc)


     if(old_package) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,imulsc%i)
        else
! GTPSA REMOVED !    call mad_tpsa_scl(s1%j,REAL(sc,kind=DP),imulsc%j)
        endif

    master=localmaster
  END FUNCTION imulsc

  FUNCTION dscmul( sc,S1 )
    implicit none
    TYPE (taylor) dscmul
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dscmul%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dscmul)

     if(old) then
    call dacmu(s1%i,sc,temp%i)
    call dacop(temp%i,dscmul%i)
        else
! GTPSA REMOVED !           call mad_tpsa_scl(s1%j,REAL(sc,kind=DP),dscmul%j)
        endif

    master=localmaster

  END FUNCTION dscmul

  FUNCTION scmul( sc,S1 )
    implicit none
    TYPE (taylor) scmul
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       scmul%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(scmul)


     if(old_package) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scmul%i)
        else
! GTPSA REMOVED !           call mad_tpsa_scl(s1%j,REAL(sc,kind=DP),scmul%j)
        endif

    master=localmaster

  END FUNCTION scmul

  FUNCTION iscmul( sc,S1 )
    implicit none
    TYPE (taylor) iscmul
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       iscmul%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(iscmul)

     if(old_package) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscmul%i)
        else
! GTPSA REMOVED !           call mad_tpsa_scl(s1%j,REAL(sc,kind=DP),iscmul%j)
        endif

    master=localmaster

  END FUNCTION iscmul

  FUNCTION div( S1, S2 )
    implicit none
    TYPE (taylor) div
    TYPE (taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       div%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    !    call check(s2)
    call ass(div)

     if(old_package) then
    call dadiv(s1%i,s2%i,temp%i)
    call dacop(temp%i,div%i)
        else
! GTPSA REMOVED !           call mad_tpsa_div(s1%j,s2%j,div%j)
        endif

    master=localmaster
  END FUNCTION div

  FUNCTION dscdiv( sc,S1 )
    implicit none
    TYPE (taylor) dscdiv
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dscdiv%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dscdiv)

     if(old_package) then
    call dadic(s1%i,sc,temp%i)
    call dacop(temp%i,dscdiv%i)
        else
        temp=sc
! GTPSA REMOVED !         call mad_tpsa_div(temp%j,s1%j,dscdiv%j)
        endif

    master=localmaster

  END FUNCTION dscdiv

  FUNCTION scdiv( sc,S1 )
    implicit none
    TYPE (taylor) scdiv
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       scdiv%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(scdiv)


     if(old_package) then
    call dadic(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scdiv%i)
        else
        temp=REAL(sc,kind=DP)
! GTPSA REMOVED !         call mad_tpsa_div(temp%j,s1%j,scdiv%j)    !    endif
    endif
    master=localmaster
  END FUNCTION scdiv

  FUNCTION iscdiv( sc,S1 )
    implicit none
    TYPE (taylor) iscdiv
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       iscdiv%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(iscdiv)

     if(old_package) then
    call dadic(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscdiv%i)
        else
        temp=REAL(sc,kind=DP)
! GTPSA REMOVED !         call mad_tpsa_div(temp%j,s1%j,iscdiv%j)    !    endif
    endif

    master=localmaster
  END FUNCTION iscdiv

  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (taylor) ddivsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    real(dp) v
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       ddivsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(ddivsc)


     if(old_package) then
    call dacdi(s1%i,sc,temp%i)
    call dacop(temp%i,ddivsc%i)
        else
        v=1.0_dp/sc
! GTPSA REMOVED !       !  call mad_tpsa_div(s1%j,temp%j,ddivsc%j)    !    endif
! GTPSA REMOVED !          call mad_tpsa_scl(s1%j,v,ddivsc%j)
     endif

    master=localmaster

  END FUNCTION ddivsc

  FUNCTION divsc( S1, sc )
    implicit none
    TYPE (taylor) divsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    real(dp) v
    IF(.NOT.C_%STABLE_DA) then
       divsc%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(divsc)

     if(old_package) then
    call dacdi(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,divsc%i)
        else
        v=1.0_dp/REAL(sc,kind=DP)
! GTPSA REMOVED !       !  call mad_tpsa_div(s1%j,temp%j,ddivsc%j)    !    endif
! GTPSA REMOVED !          call mad_tpsa_scl(s1%j,v,divsc%j)
     endif

    master=localmaster

  END FUNCTION divsc


  FUNCTION idivsc( S1, sc )
    implicit none
    TYPE (taylor) idivsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    real(dp) v
    IF(.NOT.C_%STABLE_DA) then
       idivsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(idivsc)


     if(old_package) then
    call dacdi(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,idivsc%i)
        else
        v=1.0_dp/REAL(sc,kind=DP)
! GTPSA REMOVED !       !  call mad_tpsa_div(s1%j,temp%j,ddivsc%j)    !    endif
! GTPSA REMOVED !          call mad_tpsa_scl(s1%j,v,idivsc%j)
     endif

    master=localmaster

  END FUNCTION idivsc


  FUNCTION add( S1, S2 )
    implicit none
    TYPE (taylor) add
    TYPE (taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       add%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    !    call check(s2)
    call ass(add)


     if(old) then
    call daadd(s1%i,s2%i,add%i)
       else
! GTPSA REMOVED !           call mad_tpsa_add(s1%j,s2%j,add%j)
        endif

    master=localmaster

  END FUNCTION add


  FUNCTION addq( S1, S2 )
    implicit none
    TYPE (quaternion) addq
    TYPE (quaternion), INTENT (IN) :: S1, S2


       addq%x=s1%x+s2%x
  END FUNCTION addq


  FUNCTION subq( S1, S2 )
    implicit none
    TYPE (quaternion) subq
    TYPE (quaternion), INTENT (IN) :: S1, S2

          subq%x=s1%x-s2%x

  END FUNCTION subq

  FUNCTION mulq( S1, S2 )
    implicit none
    TYPE (quaternion) mulq
    TYPE (quaternion), INTENT (IN) :: S1, S2
    integer i


          mulq=0.0_dp

          mulq%x(0)=s1%x(0)*s2%x(0)-s1%x(1)*s2%x(1)-s1%x(2)*s2%x(2)-s1%x(3)*s2%x(3)

         mulq%x(1)=  s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
         mulq%x(2)=  s1%x(3)*s2%x(1)-s1%x(1)*s2%x(3)
         mulq%x(3)=  s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)

        do i=1,3
         mulq%x(i)= mulq%x(i) + s1%x(0)*s2%x(i)+ s1%x(i)*s2%x(0)
        enddo

  END FUNCTION mulq

  FUNCTION divq( S1, S2 )
    implicit none
    TYPE (quaternion) divq
    TYPE (quaternion), INTENT (IN) :: S1, S2



       divq=s1*invq(s2)

  END FUNCTION divq

!!!! complex_quaternion


  FUNCTION caddq( S1, S2 )
    implicit none
    TYPE (complex_quaternion) caddq
    TYPE (complex_quaternion), INTENT (IN) :: S1, S2


       caddq%x=s1%x+s2%x
  END FUNCTION caddq


  FUNCTION csubq( S1, S2 )
    implicit none
    TYPE (complex_quaternion) csubq
    TYPE (complex_quaternion), INTENT (IN) :: S1, S2


          csubq%x=s1%x-s2%x

  END FUNCTION csubq

  FUNCTION cmulq( S1, S2 )
    implicit none
    TYPE (complex_quaternion) cmulq
    TYPE (complex_quaternion), INTENT (IN) :: S1, S2
    integer i


          cmulq=0.0_dp

          cmulq%x(0)=s1%x(0)*s2%x(0)-s1%x(1)*s2%x(1)-s1%x(2)*s2%x(2)-s1%x(3)*s2%x(3)

         cmulq%x(1)=  s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
         cmulq%x(2)=  s1%x(3)*s2%x(1)-s1%x(1)*s2%x(3)
         cmulq%x(3)=  s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)

        do i=1,3
         cmulq%x(i)= cmulq%x(i) + s1%x(0)*s2%x(i)+ s1%x(i)*s2%x(0)
        enddo

  END FUNCTION cmulq

  FUNCTION cmulqc( S1, c )
    implicit none
    TYPE (complex_quaternion) cmulqc
    TYPE (complex_quaternion), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: c
    integer i


          cmulqc%x= s1%x*c

  END FUNCTION cmulqc

  FUNCTION ccmulq( c,s1 )
    implicit none
    TYPE (complex_quaternion) ccmulq
    TYPE (complex_quaternion), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: c
    integer i

          ccmulq%x= s1%x*c

  END FUNCTION ccmulq

  FUNCTION cmulqr( S1, c )
    implicit none
    TYPE (complex_quaternion) cmulqr
    TYPE (complex_quaternion), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: c
    integer i


          cmulqr%x= s1%x*c

  END FUNCTION cmulqr

  FUNCTION rcmulq( c,s1 )
    implicit none
    TYPE (complex_quaternion) rcmulq
    TYPE (complex_quaternion), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: c
    integer i



          rcmulq%x= s1%x*c

  END FUNCTION rcmulq


  FUNCTION cdivq( S1, S2 )
    implicit none
    TYPE (complex_quaternion) cdivq
    TYPE (complex_quaternion), INTENT (IN) :: S1, S2



       cdivq=s1*cinvq(s2)

  END FUNCTION cdivq

!!!!!

  FUNCTION daddsc( S1, sc )
    implicit none
    TYPE (taylor) daddsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       daddsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(daddsc)

     if(old_package) then
    call dacad(s1%i,sc,temp%i)
    call dacop(temp%i,daddsc%i)
        else
         daddsc=s1
! GTPSA REMOVED !           call mad_tpsa_set0(daddsc%j,1.0_dp,sc) ! t[0] = a*t[0]+b
        endif
    master=localmaster

  END FUNCTION daddsc

  FUNCTION addsc( S1, sc )
    implicit none
    TYPE (taylor) addsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       addsc%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(addsc)


     if(old_package) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,addsc%i)
        else
         addsc=s1
! GTPSA REMOVED !           call mad_tpsa_set0(addsc%j,1.0_dp,REAL(sc,kind=DP)) ! t[0] = a*t[0]+b
        endif
    master=localmaster

  END FUNCTION addsc

  FUNCTION iaddsc( S1, sc )
    implicit none
    TYPE (taylor) iaddsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       iaddsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(iaddsc)

     if(old_package) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iaddsc%i)
        else
         iaddsc=s1
! GTPSA REMOVED !           call mad_tpsa_set0(iaddsc%j,1.0_dp,REAL(sc,kind=DP)) ! t[0] = a*t[0]+b
        endif
    master=localmaster

  END FUNCTION iaddsc

  FUNCTION dscadd( sc,S1)
    implicit none
    TYPE (taylor) dscadd
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dscadd%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dscadd)

     if(old_package) then
    call dacad(s1%i,sc,temp%i)
    call dacop(temp%i,dscadd%i)
        else
         dscadd=s1
! GTPSA REMOVED !           call mad_tpsa_set0(dscadd%j,1.0_dp,sc) ! t[0] = a*t[0]+b
        endif
    master=localmaster

  END FUNCTION dscadd

  FUNCTION scadd( sc,S1)
    implicit none
    TYPE (taylor) scadd
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       scadd%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(scadd)

     if(old_package) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scadd%i)
        else
         scadd=s1
! GTPSA REMOVED !           call mad_tpsa_set0(scadd%j,1.0_dp,REAL(sc,kind=DP)) ! t[0] = a*t[0]+b
        endif

    master=localmaster

  END FUNCTION scadd

  FUNCTION iscadd( sc,S1)
    implicit none
    TYPE (taylor) iscadd
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       iscadd%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(iscadd)


     if(old_package) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscadd%i)
        else
         iscadd=s1
! GTPSA REMOVED !           call mad_tpsa_set0(iscadd%j,1.0_dp,REAL(sc,kind=DP)) ! t[0] = a*t[0]+b
        endif
    master=localmaster

  END FUNCTION iscadd

  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (taylor) subs
    TYPE (taylor), INTENT (IN) :: S1, S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       subs%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    !    call check(s2)
    call ass(subs)



     if(old_package) then
    call dasub(s1%i,s2%i,temp%i)
    call dacop(temp%i,subs%i)
        else
! GTPSA REMOVED !       call mad_tpsa_sub(s1%j,s2%j,subs%j)
        endif
    master=localmaster

  END FUNCTION subs

  FUNCTION dsubsc( S1, sc )
    implicit none
    TYPE (taylor) dsubsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dsubsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dsubsc)

      if(old_package) then
    call dacsu(s1%i,sc,temp%i)
    call dacop(temp%i,dsubsc%i)
    else
       temp=sc
! GTPSA REMOVED !       call mad_tpsa_sub(s1%j,temp%j,dsubsc%j)
    endif

    master=localmaster


  END FUNCTION dsubsc

  FUNCTION subsc( S1, sc )
    implicit none
    TYPE (taylor) subsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       subsc%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(subsc)

     if(old_package) then
    call dacsu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,subsc%i)
   else
       temp=REAL(sc,kind=DP)
! GTPSA REMOVED !       call mad_tpsa_sub(s1%j,temp%j,subsc%j)
    endif
    master=localmaster

  END FUNCTION subsc

  FUNCTION isubsc( S1, sc )
    implicit none
    TYPE (taylor) isubsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       isubsc%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(isubsc)

     if(old_package) then
    call dacsu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,isubsc%i)
    else
       temp=REAL(sc,kind=DP)
! GTPSA REMOVED !       call mad_tpsa_sub(s1%j,temp%j,isubsc%j)
    endif
    master=localmaster

  END FUNCTION isubsc

  FUNCTION dscsub( sc,S1)
    implicit none
    TYPE (taylor) dscsub
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       dscsub%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(dscsub)

     if(old_package) then
    call dasuc(s1%i,sc,temp%i)
    call dacop(temp%i,dscsub%i)
    else
       temp=REAL(sc,kind=DP)
! GTPSA REMOVED !       call mad_tpsa_sub(temp%j,s1%j,dscsub%j)
    endif
    master=localmaster

  END FUNCTION dscsub

  FUNCTION scsub( sc,S1)
    implicit none
    TYPE (taylor) scsub
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       scsub%i=0
      RETURN
    endif
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(scsub)

     if(old_package) then
    call dasuc(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scsub%i)
     else
       temp=REAL(sc,kind=DP)
! GTPSA REMOVED !       call mad_tpsa_sub(temp%j,s1%j,scsub%j)
    endif
    master=localmaster

  END FUNCTION scsub

  FUNCTION iscsub( sc,S1)
    implicit none
    TYPE (taylor) iscsub
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       iscsub%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(iscsub)

     if(old_package) then
    call dasuc(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscsub%i)
     else
       temp=REAL(sc,kind=DP)
! GTPSA REMOVED !       call mad_tpsa_sub(temp%j,s1%j,iscsub%j)
    endif
    master=localmaster

  END FUNCTION iscsub


  !  These are new general TPSA-Routines



  FUNCTION varf( S1, S2 )
    implicit none
    TYPE (taylor) varf
    real(dp), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       varf%i=0
      RETURN
    endif
    localmaster=master


    call ass(varf)

    varf=S1 + (1.0_dp.mono.S2)

    master=localmaster

  END FUNCTION varf

  FUNCTION varf001( S1, S2 )
    implicit none
    TYPE (taylor) varf001
    real(dp), INTENT (IN) :: S1(2)
    integer  , INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       varf001%i=0
      RETURN
    endif
    localmaster=master


    call ass(varf001)

    varf001=S1(1) + (s1(2).mono.S2)

    master=localmaster

  END FUNCTION varf001




  SUBROUTINE  shift000(S1,S2,s)
    implicit none
    INTEGER,INTENT(IN)::s
    type (taylor),INTENT(IN)::S1
    type (taylor),INTENT(inout)::S2
    type(taylor) junk
    integer j(lnv),jj(lnv),i,mf
    real(dp) value
    IF(.NOT.C_%STABLE_DA) RETURN


!    if(s2%i==0) call crap1("shift000  1" )  !call etall1(s2%i)
     if(old_package) then
        CALL DAshift(s1%i,s2%i,s)
       else
    j=0
    call alloc(junk)
   junk=0.0_dp
   call taylor_cycle(S1,size=mf)

   do i=1,mf
   call taylor_cycle(s1,ii=i,value=VALUE,j=JJ)

   if(value/=0.d0) then
    j=0
    j(1:c_%nv-s)  =   jj(s+1:c_%nv)
    junk=junk + (VALUE.mono.J)
   endif
  enddo
    s2=junk
    call kill(junk)

   endif

  END SUBROUTINE shift000


  SUBROUTINE  pek0(S1,R1,j)
    implicit none
    INTEGER,INTENT(IN)::j
    real(dp),INTENT(inOUT)::R1(:)
    type (taylor),INTENT(IN)::S1(:)
    integer i
 !   integer k
    IF(.NOT.C_%STABLE_DA) RETURN
      do i=1,j
       r1(i)=s1(i)
      enddo
    !
  END SUBROUTINE pek0

  SUBROUTINE  pok0(S1,R1,j)
    implicit none
    INTEGER,INTENT(IN)::j
    real(dp),INTENT(in)::R1(:)
    type (taylor),INTENT(INout)::S1(:)
    integer i,jj(lnv)
 !   integer k
      jj=0
    IF(.NOT.C_%STABLE_DA) RETURN
      do i=1,j
       call pok(s1(i),jj,r1(i))
      enddo
    !
  END SUBROUTINE pok0


  SUBROUTINE  pek000(S1,J,R1)
    implicit none
    INTEGER,INTENT(IN),dimension(:)::j
    real(dp),INTENT(inOUT)::R1
    type (taylor),INTENT(IN)::S1
 !   integer k
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
   ! if(s1%i==0) call crap1("pek000  1" )  !call etall1(s1%i)
 !   k=s1%i
!    write(6,*) r1,k
     if(.not.check_gtpsa(j=j)) then  
         r1=0
     !   write(6,*) "error in pek000"
     !    write(6,*) j
     !    write(6,*)"c_%nv,c_%no",c_%nv,c_%no
       return
     endif
    if(old_package) then
       CALL DApek(s1%i,j,r1)
        else
         mo_gtpsa=0
         mo_gtpsa=j(1:size(j))
! GTPSA REMOVED !        r1=mad_tpsa_getm(s1%j,c_%nv,mo_gtpsa)
    endif
    !
  END SUBROUTINE pek000

  SUBROUTINE  pok000(S1,J,R1)
    implicit none
    INTEGER,INTENT(in),dimension(:)::j
    real(dp),INTENT(in)::R1
    type (taylor),INTENT(inout)::S1
    integer i,siz
    integer(1) jj(lnv)
    IF(.NOT.C_%STABLE_DA) RETURN

if(check_gtpsa(j=j)) then
 else
 ! write(6,*) " check_gtpsa failed pok000 "
  return
!  stop 999
 endif
!    if(check_j(j)/=0) return
    ! if(old) then
    !if(s1%i==0) call crap1("pok000 1" )  ! call etall1(s1%i)
    if(old_package) then
      CALL DApok(s1%i,j,r1)
    else

     jj=0
     jj(1:size(j))=j
     siz=min(size(j),nv)

! GTPSA REMOVED !      call  mad_tpsa_setm(s1%j,siz,jj ,0.0_dp,r1)
    endif

  END SUBROUTINE pok000

  SUBROUTINE  TAYLOR_ran(S1,r1,R2)
    implicit none
    real(dp),INTENT(in)::R1
    real(dp),INTENT(inout)::R2
    integer j(lnv),mf,i
    type (taylor),INTENT(inout)::S1
    type(taylor) junk
    real(dp) value
    IF(.NOT.C_%STABLE_DA) RETURN

    !
    !     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
    !     FOR R1 > 0, THE VECTOR IS FILLED WITH REALS,
    !     FOR R1 < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
    !     ABS(R1) IS THE FILLING FACTOR
    if(old) then
    if(s1%i==0) call crap1("tAYLOR_ran  1" )  ! call etall1(s1%i)
    call daran(s1%i,r1,R2)
        else

    j=0
    call alloc(junk)
   junk=0.0_dp
   call taylor_cycle(S1,size=mf)

   do i=1,mf
   call taylor_cycle(s1,ii=i,value=VALUE,j=J)
         if(r1.gt.0.0_dp) then
             VALUE = bran(R2)
             if(VALUE.gt.r1) VALUE = 0.0_dp
          elseif(r1.lt.0.0_dp) then
             VALUE = int(1+10*bran(R2))
             if(VALUE.gt.-ten*r1) VALUE = 0.0_dp
          endif

     junk=junk + (VALUE.mono.J)

  enddo
    s1=junk
    call kill(junk)

       endif
    !
  END SUBROUTINE TAYLOR_ran

! obsolete
!  SUBROUTINE  intd_taylor(S1,S2,factor)
!    implicit none
!    type (taylor),INTENT(inOUT)::S2
!   type (taylor),INTENT(IN)::S1(:)
!    real(dp),INTENT(IN):: factor
!    IF(.NOT.C_%STABLE_DA) RETURN

!    CALL intd(S1%i,s2%i,factor)

!  END SUBROUTINE intd_taylor

!  SUBROUTINE  DIFd_taylor(S2,S1,factor)
!    implicit none
!    type (taylor),INTENT(in)::S2
!    type (taylor),INTENT(INOUT)::S1(:)
!    real(dp),INTENT(IN):: factor
!    IF(.NOT.C_%STABLE_DA) RETURN
!
!    CALL DIFD(S2%i,s1%i,factor)
!
!  END SUBROUTINE DIFd_taylor


  SUBROUTINE  CFU000(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    real(dp) FUN
    EXTERNAL FUN
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
    if(s1%i==0) call crap1("CFU000  1" )  !  call etall1(s1%i)
    CALL DACFU(s2%i,FUN,s1%i)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r))call crap1("CFU000  2" )  !  call newetall(s1%j,1)
    !       CALL NEWDACFU(s2%J,FUN,s1%J)
    !    endif

  END SUBROUTINE CFU000

  SUBROUTINE  DACFU_GTPSA(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    real(dp) VALUE
    INTEGER MF,I,J(LNV)
    type (taylor)   junk
    interface
       function fun(abc)
         use precision_constants
         implicit none
         real(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface

    IF(.NOT.C_%STABLE_DA) RETURN
    j=0
    call alloc(junk)
   junk=0.0_dp
   call taylor_cycle(S2,size=mf)

   do i=1,mf
   call taylor_cycle(s2,ii=i,value=VALUE,j=J)
!   write(6,*) i,value

   if(value/=0.d0) then
    junk=junk + fun(j)*(VALUE.mono.J)
   endif
  enddo
    s1=junk
    call kill(junk)
  END SUBROUTINE DACFU_GTPSA

  SUBROUTINE  CFU000_new (S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    real(dp) FUN
    EXTERNAL FUN

    IF(.NOT.C_%STABLE_DA) RETURN
      if(old_package) then
          call CFU000 (S2,FUN,S1)
         else
            call DACFU_GTPSA(S2,FUN,S1)
         endif

  END SUBROUTINE CFU000_new

! obsolete
!  SUBROUTINE  CFUR(S2,FUN,S1)
!    implicit none
!    type (taylor),INTENT(INOUT)::S1
!    type (taylor),INTENT(IN)::S2
!    complex(dp) FUN
!    EXTERNAL FUN
!    IF(.NOT.C_%STABLE_DA) RETURN
!
!    if(s1%i==0) call crap1("CFUR  1" )  ! call etall1(s1%i)
!    CALL DACFUR(s2%i,FUN,s1%i)
!
!  END SUBROUTINE CFUR

! obsolete
!  SUBROUTINE  CFUI(S2,FUN,S1)
!    implicit none
!    type (taylor),INTENT(INOUT)::S1
!    type (taylor),INTENT(IN)::S2
!    complex(dp) FUN
!    EXTERNAL FUN
!    IF(.NOT.C_%STABLE_DA) RETURN
!    ! if(old) then
!    if(s1%i==0)call crap1("CFUI  1" ) ! call etall1(s1%i)
!    CALL DACFUI(s2%i,FUN,s1%i)
!
!  END SUBROUTINE CFUI

  SUBROUTINE  taylor_eps(r1)
    implicit none
    real(dp),INTENT(INOUT)::r1
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
    CALL DAeps(r1)
    !   else
    !      CALL newDAeps(r1)
    !   endif

  END SUBROUTINE taylor_eps



  FUNCTION GETCHARnd2( S1, S2 )
    implicit none
    TYPE (taylor) GETCHARnd2,junk
    TYPE (taylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    CHARACTER (LEN = LNV)  resul
    integer i,k
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETCHARnd2%i=0
      RETURN
    endif
    localmaster=master

    ndel=0
    !    call check(s1)
    call ass(GETCHARnd2)


    resul = trim(ADJUSTR (s2))

    do i=1,lnv
       jfil(i)=0
    enddo

    nd2parfilter= len(trim(ADJUSTR (s2)))

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2parfilter
       CALL  CHARINT(RESUL(I:I),Jfil(I))
       if(i>nv) then
          if(Jfil(i)>0) then
             GETCHARnd2=0.0_dp
          write(6,*) " error jfil(i)/=0 for .para. "
               stop  101

          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2parfilter+1,nv
       if(jfil(i)/=0) then

         write(6,*) " error in getchar for .para. "

          stop
       endif
    enddo
    call alloc(junk)
    call cfu(s1,filter,junk)


    DO I=1,nd2parfilter
       DO K=1,JFIL(I)
          JUNK=JUNK.K.I
       ENDDO
    ENDDO

    GETCHARnd2=junk

    call kill(junk)
    master=localmaster

  END FUNCTION GETCHARnd2

  FUNCTION GETintnd2( S1, S2 )
    implicit none
    TYPE (taylor) GETintnd2,junk
    TYPE (taylor), INTENT (IN) :: S1
    integer , INTENT (IN) ::  S2(:)
    integer i,k
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETintnd2%i=0
      RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(GETintnd2)

    call alloc(junk)

    do i=1,lnv
       jfil(i)=0
    enddo
    nd2parfilter=size(s2)
    ndel=0

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2parfilter
       Jfil(I)=s2(i)
       if(i>nv) then
          if(Jfil(i)>0) then
             GETintnd2=0.0_dp
          write(6,*) " error in i>nv for .para. "
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2parfilter+1,nv
       if(jfil(i)/=0) then
          write(6,*) " error in jfil(i)/=0 for .para. "
          stop 102
       endif
    enddo

    call cfu(s1,filter,junk)

    !DO I=1,ND2+ndel
    DO I=1,nd2parfilter
       DO K=1,JFIL(I)
          JUNK=JUNK.K.I
       ENDDO
    ENDDO

    GETintnd2=junk

    call kill(junk)
    master=localmaster

  END FUNCTION GETintnd2

  FUNCTION GETintnd2t( S1, S22 )
    implicit none
    TYPE (taylor) GETintnd2t,junk
    TYPE (taylor), INTENT (IN) :: S1
    type(sub_taylor), INTENT (IN) :: S22
    integer s2(lnv)
    integer i
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETintnd2t%i=0
      RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(GETintnd2t)

    call alloc(junk)

    do i=1,lnv
       jfilt(i)=0
    enddo
    s2=s22%j
    nd2part=s22%min
    nd2partt=s22%max
    ndel=0
    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=nd2part,nd2partt
       jfilt(I)=s2(i)
       if(i>nv) then
          if(jfilt(i)>0) then
             GETintnd2t=0.0_dp
             return
          endif
       endif

    enddo

    !do i=nd2+ndel+1,nv
    do i=nd2partt+1,nv
       if(jfilt(i)/=0) then

            write(6,*) " error in GETintnd2t for .part_taylor. "
          ! call !write_e(0)
          stop
       endif
    enddo

    call cfu(s1,filter_part,junk)

    !DO I=1,ND2+ndel
    !    DO I=1,ND2par
    !       DO K=1,jfilt(I)
    !          JUNK=JUNK.K.I
    !       ENDDO
    !    ENDDO

    GETintnd2t=junk

    call kill(junk)
    master=localmaster

  END FUNCTION GETintnd2t

  SUBROUTINE  taylor_cycle(S1,size,ii,VALUE,J)
    implicit none
    type (taylor),INTENT(IN)::S1
    integer,optional, intent(inout):: size
    integer,optional, intent(in):: ii
    integer,optional, intent(inout)::J(:)
    real(dp), OPTIONAL, intent(inout):: value
     if(old_package) then
      call taylor_cycle_dabnew(S1,size,ii,VALUE,J)
     else
      call taylor_cycle_gtpsa(S1,size,ii,VALUE,J)
     endif
  end SUBROUTINE  taylor_cycle


  SUBROUTINE  taylor_cycle_dabnew(S1,size,ii,VALUE,J)
    implicit none
    type (taylor),INTENT(IN)::S1
    integer,optional, intent(inout):: size
    integer,optional, intent(in):: ii
    integer,optional, intent(inout)::J(:)
    real(dp), OPTIONAL, intent(inout):: value
    INTEGER ipresent,ILLA
    real(dp) VALUE0
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) THEN
    IF(PRESENT(J).AND.PRESENT(VALUE).and.present(ii)) THEN
       call dacycle(S1%i,ii,value,illa,J)
    ELSEif(present(size)) then
       call dacycle(S1%i,ipresent,value0,size)
    else
       write(6,*) "error in taylor_cycle"
       stop 888
    ENDIF

  END SUBROUTINE taylor_cycle_dabnew

  SUBROUTINE  taylor_cycle_gtpsa(S1,size,ii,VALUE,J)
    implicit none
    type (taylor),INTENT(IN)::S1
    integer,optional, intent(inout):: size
    integer,optional, intent(in):: ii
    integer,optional, intent(inout)::J(:)
    real(dp), OPTIONAL, intent(inout):: value
    INTEGER err,ii_gtpsa
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) THEN
    IF(PRESENT(J).AND.PRESENT(VALUE).and.present(ii)) THEN
   !    call dacycle(S1%i,ii,value,illa,J)
! GTPSA REMOVED !        ii_gtpsa=ii-1
! GTPSA REMOVED !   VALUE=mad_tpsa_geti(S1%j,ii_gtpsa)
! GTPSA REMOVED !  err=mad_tpsa_mono(S1%j,c_%nv,mo_gtpsa,ii_gtpsa)
! GTPSA REMOVED !    j(1:c_%nv)=mo_gtpsa
! GTPSA REMOVED !  if(err<0) then
! GTPSA REMOVED !     write(6,*) " error in taylor_cycle_gtpsa "
! GTPSA REMOVED !   stop
! GTPSA REMOVED !  endif

    ELSEif(present(size)) then
! GTPSA REMOVED !       size=mad_desc_maxlen(d_berz)
    else
       write(6,*) "error in taylor_cycle"
       stop 888
    ENDIF

  END SUBROUTINE taylor_cycle_gtpsa

 ! SUBROUTINE  taylor_clean(S1,VALUE)
 !   implicit none
 !   type (taylor),INTENT(INout)::S1
 !   real(dp) value
 !   call daclean(S1%i,value)
 ! END SUBROUTINE taylor_clean

  subroutine check_snake()
    implicit none
    master=master+1
    select case (master)
    case(1:ndumt)
       if(iass0user(master)>scratchda(master)%n.or.scratchda(master)%n>newscheme_max) then
          call ndum_warning_user
       endif
       iass0user(master)=0
    case(ndumt+1:)

         write(6,*) "Should not be here in check_snake"

       ! call !write_e(101)
    end select
    master=master-1
  end subroutine check_snake

  ! functions used inside other routines

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


  function check_j(j)
    implicit none
    integer check_j
    INTEGER,INTENT(in),dimension(:)::j
    integer i,no

    IF(.NOT.C_%STABLE_DA) then
      check_j=0
     RETURN
    endif

    check_j=0

    no=0
    do i=1,size(j)
       no=j(i)+no
    enddo

    if(no>c_%no) then
       check_j=no
       return
    endif

    do i=c_%nv+1,size(j)
       if(j(i)/=0) then
          check_j=-i
       endif
    enddo
  end function check_j



  function filter(j)
    implicit none
    real(dp) filter
    integer i
    integer,dimension(:)::j

    filter=1.0_dp
    !do i=1,nd2+ndel
    do i=1,nd2parfilter
       if(jfil(i)/=j(i)) filter=0.0_dp
    enddo

  end  function filter

  function filter_part(j)
    implicit none
    real(dp) filter_part
    integer i
    integer,dimension(:)::j
    !    WRITE(6,*) jfilt(1:4)
    !    WRITE(6,*)nd2part,nd2partt
    filter_part=1.0_dp
    !do i=1,nd2+ndel
    do i=nd2part,nd2partt
       if(jfilt(i)/=j(i)) filter_part=0.0_dp
    enddo

  end  function filter_part

  !  i/o routines

  SUBROUTINE  printq(S1,MFILE,PREC)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (quaternion),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    INTEGER I,mfi
     mfi=6
     if(present(mfile)) mfi=mfile
      write(mfi,*) " real quaternion "
    DO I=0,3
      write(mfi,*) s1%x(i)
    ENDDO
  END SUBROUTINE printq


  SUBROUTINE  cprintq(S1,MFILE,PREC,pr)
    implicit none
    INTEGER,OPTIONAL,INTENT(IN)::MFILE
    type (complex_quaternion),INTENT(IN)::S1
    REAL(DP),OPTIONAL,INTENT(IN)::PREC
    logical,OPTIONAL,INTENT(INout)::pr
    real(dp) norm
    INTEGER I,mfi
     mfi=6
     if(present(mfile)) mfi=mfile
     if(present(prec) ) then
     norm=0
      do i=0,3
       norm=norm+abs(s1%x(i))
      enddo
     if(norm>prec) then
     if(present(pr))pr=.true.
     if(mfi/=0) then
     write(mfi,*) " complex_quaternion "
       DO I=0,3
         write(mfi,*) s1%x(i)
       ENDDO
      endif
      else
            if(present(pr))pr=.false.
      endif
     else

      write(mfi,*) " complex_quaternion "
    DO I=0,3
      write(mfi,*) s1%x(i)
    ENDDO
   endif
  END SUBROUTINE cprintq

  SUBROUTINE  print_for_bmad_parse(S1,MFILE,prec,ind)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(IN)::prec
    integer ,OPTIONAL,INTENT(IN)::ind
    type (TAYLOR),INTENT(IN)::S1

     bmadparser=1
    call pri(S1,MFILE,prec,ind)
     bmadparser=0

  end SUBROUTINE  print_for_bmad_parse

  SUBROUTINE  pri(S1,MFILE,prec,ind)
    implicit none
        INTEGER,OPTIONAL,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(IN)::prec
    integer ,OPTIONAL,INTENT(IN)::ind
    type (TAYLOR),INTENT(IN)::S1
    REAL(DP) PREC1,depst,value
    integer i,j,it,n,indo,k,kt,kl,mfi
    integer, allocatable :: jc(:)
    character(255) line,line0
    mfi=6
    if(present(mfile)) mfi=mfile
    IF(PRESENT(prec)) THEN
       PREC1=-1.0_dp
       depst=prec
       CALL taylor_eps(PREC1)
       CALL taylor_eps(depst)
    ENDIF

 if(bmadparser>0) then
    IF(PRESENT(ind)) THEN
     indo=ind
    else
     indo=0
    endif

    kt=0
    allocate(jc(c_%nv))
     call taylor_cycle(s1,size=n)
    do i=1,n
       call taylor_cycle(s1,ii=i,value=value,j=jc)

       it=0
       do j=c_%nd2+1,c_%nv
          it=jc(i)+it
       enddo
       if(it==0.and.abs(value)>depst) then
        kt=kt+1
       endif

    enddo

    deallocate(jc)

    allocate(jc(c_%nv))
    kl=0
     call taylor_cycle(s1,size=n)
    do i=1,n
       call taylor_cycle(s1,ii=i,value=value,j=jc)

       it=0
       do j=c_%nd2+1,c_%nv
          it=jc(i)+it
       enddo
       if(it==0.and.abs(value)>depst) then
        kl=kl+1
        write(line,*) "{",indo,":",value,","
        call context(line)
        do j=1,c_%nd2
         write(line(len_trim(line)+1:255),*)jc(j),"&"
         call context(line)
        enddo
        if(kl==kt.and.indo==c_%nd2) then
          write(line(len_trim(line)+1:255),*)"}"
        else
          write(line(len_trim(line)+1:255),*)"},"
        endif
         call context(line)
         k=0
         line0=' '
         do j=1,len_trim(line)
          if(line(j:j)/=' ') then
           !line(j:j)=' '
           !else
          k=k+1
           line0(k:k)=line(j:j)
          endif
         enddo
         do j=1,len_trim(line0)
          if(line0(j:j)=='&') line0(j:j)=' '
         enddo
        write(mfi,*) line0(1:len_trim(line0))
       endif

    enddo

    deallocate(jc)


   else




     if(old_package) then
    if(print77) then
       CALL DAPRI77(s1%i,mfi)
    else
       CALL DAPRI(s1%i,mfi)
    endif
    else
      PREC1=0.0_dp
      if(present(prec)) prec1=prec
     if(mfi==6) then
! GTPSA REMOVED !       call mad_tpsa_print(s1%j, "no name"//c_eos, prec1, 0, c_null);  ! 1 enleve le header nv, etc..
       else
     endif
    endif
    !
endif
    IF(PRESENT(prec))  CALL taylor_eps(PREC1)

  END SUBROUTINE pri

  SUBROUTINE  REA(S1,MFILE)
    implicit none
    INTEGER,INTENT(in)::MFILE
    type (TAYLOR),INTENT(IN)::S1

    ! if(old) then
    if(s1%i==0)call crap1("REA  1" ) !  call etall1(s1%i)

    if(read77) then
       CALL DAREA77(s1%i,MFILE)
    else
       CALL DAREA(s1%i,MFILE)
    endif
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r))call crap1("REA  2" ) ! call newetall(s1%j,1)
    !       if(newread) then
    !          CALL newDAREA(s1%j,MFILE)
    !       else
    !          if(read77) then
    !             CALL oldDAREA77(s1%j,MFILE)
    !          else
    !             CALL oldDAREA(s1%j,MFILE)
    !          endif
    !       endif
    !    endif

  END SUBROUTINE REA


  ! Universal Taylor Routines   (Sagan's Stuff)

  SUBROUTINE  kill_uni(S2)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2

    DEALLOCATE(S2%N,S2%NV,S2%C,S2%J)
    NULLIFY(S2%N,S2%NV,S2%C,S2%J)

  END SUBROUTINE  kill_uni

  SUBROUTINE  null_uni(S2,S1)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    integer, intent(in):: s1
    IF(S1==0) THEN
       NULLIFY(S2%N,S2%NV,S2%C,S2%J)
    ELSEIF(S1==-1) THEN
       DEALLOCATE(S2%N,S2%NV,S2%C,S2%J)
       NULLIFY(S2%N,S2%NV,S2%C,S2%J)
    ENDIF
  END SUBROUTINE null_uni


  SUBROUTINE  ALLOC_U(S2,N,NV)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    integer, intent(in):: N,NV
    ALLOCATE(S2%N,S2%NV)
    if(N==0) then
       allocate(S2%C(1),S2%J(1,NV));S2%C(1)=0.0_dp;S2%J(:,:)=0;
    else
       allocate(S2%C(N),S2%J(N,NV))
    endif
    S2%N=N
    S2%NV=NV
  END SUBROUTINE ALLOC_U

  SUBROUTINE  fill_uni_r(S2,S1)  !new sagan
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    real (dp), intent(in):: s1
    INTEGER n,J(LNV)


    IF(ASSOCIATED(S2%N)) S2=-1
    S2=0
    CALL ALLOC_U(S2,1,nv)
    J=0
    DO N=1,S2%NV
       S2%J(1,N)=J(N)
    ENDDO
    S2%C(1)=S1

  END SUBROUTINE fill_uni_r

  SUBROUTINE  FILL_UNI(S2,S1)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(INOUT)::S2
    type (TAYLOR), intent(in):: s1
    INTEGER ipresent,k,n,I
    real(dp) value
    INTEGER, allocatable :: j(:)
    call check_snake

    ! if(old) then
   ! if(s1%i==0)  call crap1("FILL_N 1")
    !    else
    !       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("FILL_N 2")
    !    endif


    IF(ASSOCIATED(S2%N)) S2=-1
    S2=0
    ipresent=1
    call taylor_cycle(S1,n)
    CALL ALLOC_U(S2,N,c_%nv)
    allocate(j(c_%nv))

    do i=1,N
       call taylor_cycle(S1,ii=i,value=value,j=j)
       S2%C(I)=value
       DO k=1,S2%NV
          S2%J(i,k)=J(k)
       ENDDO
    ENDDO

    deallocate(j)

  END SUBROUTINE FILL_UNI




  SUBROUTINE  REFILL_UNI(S1,S2)
    implicit none
    type (UNIVERSAL_TAYLOR),INTENT(IN)::S2
    type (TAYLOR), intent(inOUT):: s1
    INTEGER I,K,J(LNV)
    logical(lp) DOIT

    ! if(old) then
    !if(s1%i==0)  call crap1("REFILL_N 1")
    !    else
    !       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("REFILL_N 2")
    !    endif


    S1=0.0_dp

    IF(.not.ASSOCIATED(S2%N)) THEN

         write(6,*) " ERROR IN REFILL_N: UNIVERSAL_TAYLOR DOES NOT EXIST"
       ! call !write_e(123)
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

  END SUBROUTINE REFILL_UNI


  !_________________________________________________________________________________


  subroutine printunitaylor(ut,iunit)
    implicit none
    type(universal_taylor) :: ut
    integer, optional               :: iunit
    integer                :: i,ii
    integer inuit0
    inuit0=6
    if(present(iunit)) inuit0=iunit
    if (.not. associated(ut%n)) then
       write(iunit,'(A)') '    UNIVERSAL_TAYLOR IS EMPTY (NOT ASSOCIATED)'
       write(6,'(A)') '    UNIVERSAL_TAYLOR IS EMPTY (NOT ASSOCIATED)'
       return
    endif

    write(iunit,'(/1X,A,I5,A,I5,A/1X,A/)') 'UNIV_TAYLOR   NO =',ut%n,', NV =',ut%nv,', INA = unita',&
         '*********************************************'
    if(ut%n /= 0) then
       write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
    else
       write(iunit,'(A)') '   ALL COMPONENTS 0.0_dp '
    endif

    do i = 1,ut%n
       write(iunit,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') i,ut%c(i),sum(ut%j(i,:)),(ut%j(i,ii),ii=1,ut%nv)
       if( .not. print77) then
          write(iunit,*)  ut%c(i)
       endif
    enddo

    write(iunit,'(A)') '                                      '

  end subroutine printunitaylor


  ! End of Universal Taylor Routines




  ! Warning Routines

  subroutine crap1(STRING)
    implicit none
    CHARACTER(*) STRING


      write(6,*) "ERROR IN :"
      write(6,*) STRING


  end subroutine crap1

  SUBROUTINE real_stop()
    implicit none
    integer i(1),j


    write(6,*) " You are using a kind(1.0_dp) "
    write(6,*)" set real_warning to false to permit this "
    write(6,*)" write 1 to continue or -1 for a crash "
    call read(j)
    i(j)=0
    real_warning=.false.

  END   SUBROUTINE real_stop


  SUBROUTINE  ndum_warning_user()
    implicit none
    integer ipause,II(0:1)




      write(6,*)  " *  Should never be here in New Linked List Scheme               *"

    call read(ipause)
    ii(2000*ipause)=0

  end SUBROUTINE  ndum_warning_user

  ! End of  Warning Routines

  ! linked list of da for scratch levels

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

  SUBROUTINE set_up_level()
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

  SUBROUTINE report_level()
    implicit none
    integer i
    if(associated(scratchda(1)%n)) then
       do i=1,ndumt

          write(6,'(a6,1x,i4,a5,1x,i4,1x,a7)') "Level ",i, " has ",scratchda(i)%n, "Taylors"
          !          write(w_p%c(1),'(a6,1x,i4,a5,1x,i4,1x,a7)') "Level ",i, " has ",scratchda(i)%n, "Taylors"
          !          ! call !write_e
       enddo
    endif
  END   SUBROUTINE report_level

  ! end linked list of da for scratch levels

  ! Assignments Routines

  subroutine ASSIGN()
    implicit none
    integer i
    do i=1,ndumt
       iassdoluser(i)=0
       iass0user(i)=0
    enddo
    ! if(old) then
!    CALL ETALL1(DUMMY)
 !   call etall1(temp)
    call alloc(temp)
    !    else
    !       CALL allocnewda(DUMMYl)
    !       call allocnewda(templ)
    !    endif
    CALL set_up_level
  end subroutine ASSIGN

  subroutine DEASSIGN()
    implicit none
    integer i
    do i=1,ndumt
       iassdoluser(i)=0
       iass0user(i)=0
    enddo
    ! if(old) then
!    CALL DADAL1(DUMMY)
    call kill(temp)
!    call DADAL1(temp)
    !    else
    !       CALL KILLnewdaS(DUMMYl)
    !       call KILLnewdaS(templ)
    !    endif
    do i=1,ndumt
       CALL kill_DALEVEL(scratchda(I))
    ENDDO
  end subroutine DEASSIGN

  subroutine ASStaylor(s1)
    implicit none
    TYPE (taylor) s1
    !  lastmaster=master  ! 2002.12.13

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)
       write(6,*) " cannot indent anymore ",ndumt

       master=sqrt(-dble(master))
    end select
    !    write(26,*) "   taylor ",master
    call ass0(s1)

! GTPSA REMOVED !    if(.not.old_package) call mad_tpsa_clear(s1%j)

  end subroutine ASStaylor

  subroutine ass0(s1)
    implicit none
    integer ipause, mypause
    TYPE (taylor) s1

    IF(MASTER>NDUMT.or.master==0) THEN
       WRITE(6,*) "more scratch level needed ",master,NDUMT
       ipause=mypause(123)
       write(6,*) 1/sqrt(-dble(1000+master))
       stop 123
    ENDIF

    if(.not.no_ndum_check) iass0user(master)=iass0user(master)+1
    if(iass0user(master)>scratchda(master)%n) then
       call INSERT_DA( scratchda(master) )
    ELSE
       scratchda(master)%PRESENT=>scratchda(master)%PRESENT%NEXT
    ENDIF
     if(old_package) then
    s1%i=scratchda(master)%PRESENT%T%i
        else
! GTPSA REMOVED !           s1%j=scratchda(master)%PRESENT%T%j
        endif


  end subroutine ASS0


  ! remove small numbers

  SUBROUTINE  clean_taylor(S1,S2,prec)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S2
    type (TAYLOR), intent(INOUT):: s1
    real(dp) prec
    INTEGER ipresent,n,I,illa
    real(dp) value
    INTEGER, allocatable :: j(:)
    type (TAYLOR) t

    call alloc(t)
    t=0.0_dp

     call taylor_cycle(s1,size=n)
    allocate(j(c_%nv))

    do i=1,N
       call taylor_cycle(S1,ii=i,value=value,j=j)
       if(abs(value)>prec) then
          t=t+(value.mono.j)
       endif
    ENDDO
    s2=t
    deallocate(j)
    call kill(t)

  END SUBROUTINE clean_taylor


  SUBROUTINE  clean_pbfield(S1,S2,prec)
    implicit none
    type (pbfield),INTENT(INOUT)::S2
    type (pbfield), intent(INOUT):: s1
    real(dp) prec

    call clean_taylor(s1%h,s2%h,prec)

  END SUBROUTINE clean_pbfield

  SUBROUTINE  clean_pbresonance (S1,S2,prec)
    implicit none
    type (pbresonance),INTENT(INOUT)::S2
    type (pbresonance), intent(INOUT):: s1
    real(dp) prec

    call clean_pbfield(s1%cos,s2%cos,prec)
   call clean_pbfield(s1%sin,s2%sin,prec)

  END SUBROUTINE clean_pbresonance

  SUBROUTINE  clean_damap(S1,S2,prec)
    implicit none
    type (damap),INTENT(INOUT)::S2
    type (damap), intent(INOUT):: s1
    real(dp) prec
    integer i

    do i=1,c_%nd2
       call clean_taylor(s1%v(i),s2%v(i),prec)
    enddo


  END SUBROUTINE clean_damap

  SUBROUTINE  clean_vecfield(S1,S2,prec)
    implicit none
    type (vecfield),INTENT(INOUT)::S2
    type (vecfield), intent(INOUT):: s1
    real(dp) prec
    integer i

    do i=1,c_%nd2
       call clean_taylor(s1%v(i),s2%v(i),prec)
    enddo


  END SUBROUTINE clean_vecfield

  SUBROUTINE  clean_vecresonance(S1,S2,prec)
    implicit none
    type (vecresonance),INTENT(INOUT)::S2
    type (vecresonance), intent(INOUT):: s1
    real(dp) prec



       call clean_vecfield(s1%cos,s2%cos,prec)
       call clean_vecfield(s1%sin,s2%sin,prec)



  END SUBROUTINE clean_vecresonance



  SUBROUTINE  clean_onelieexponent(S1,S2,prec)
    implicit none
    type (onelieexponent),INTENT(INOUT)::S2
    type (onelieexponent), intent(INOUT):: s1
    real(dp) prec



       call clean_vecfield(s1%vector,s2%vector,prec)
       call clean_pbfield(s1%pb,s2%pb,prec)



  END SUBROUTINE clean_onelieexponent


  ! remove small numbers

  SUBROUTINE  clean_complextaylor(S1,S2,prec)
    implicit none
    type (complextaylor),INTENT(INOUT)::S2
    type (complextaylor), intent(INOUT):: s1
    real(dp) prec

    call clean_taylor(S1%r,S2%r,prec)
    call clean_taylor(S1%i,S2%i,prec)


  END SUBROUTINE clean_complextaylor

  SUBROUTINE  clean_gmap(S1,s2,prec)
    implicit none
    type (gmap),INTENT(INOUT)::S1
    type (gmap),INTENT(INOUT)::S2
    real(dp) prec
    INTEGER I

    DO I=1,s1%n
       CALL clean_taylor(S1%V(I),S2%V(I),prec)
    ENDDO

  END SUBROUTINE clean_gmap

function etienne_bessel_Ir(n, x, y,km) result (value)

implicit none

real(dp) x, y, value,dvalo
real(dp) r2, rr2, denom, r2k, r, scale, dval,eps
integer, optional :: km
integer n, k,km0
integer, parameter :: nk_max = 1000
logical done
!
eps=1.e-8_dp
done=.false.
r2 = x**2 + y**2
km0=15
if(present(km)) km0=km
scale = 1.0_dp/2**n !/ (2**n * factorial(n))


do k=1,n
scale=scale/k
enddo



! Close to origin case.
! Crossover point is hurestically derived for n <= 30

!if (r2 < 2.28 * (n+7)) then
  value = 1

  rr2 = r2 / 4
  denom = 1.0_dp
  r2k = 1
  dvalo=1.d38
  do k = 1, nk_max
    r2k = r2k * rr2
    denom = denom * k * (n + k)
    dval = r2k / denom
    value = value + dval
    dvalo=dval
  if(done) then
     if(dvalo>=dval) exit
    else
    if (k>km0.and.dval < eps * value) then
      done=.true.
    endif
  endif
    if (k == nk_max) then
      print *, 'Internal error in norm_bessel_I: No convergence!'
      stop
    endif
  enddo
if(present(km)) write(6,*) k
  value = scale * value


end function etienne_bessel_Ir


function etienne_bessel_It(n, x, y,km) result (value)

implicit none

type(taylor) x, y
type(taylor)   r2, rr2, denom, r2k, r, scale, dval, value,dvalo
real(dp) eps
integer, optional :: km
integer n, k,km0
integer, parameter :: nk_max = 1000
logical done
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      value%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(value)
call alloc(r2, rr2, denom, r2k, r, scale, dval,dvalo)
eps=1.e-8_dp
done=.false.
r2 = x**2 + y**2
km0=15
if(present(km)) km0=km
scale = 1.0_dp/2**n !/ (2**n * factorial(n))


do k=1,n
scale=scale/k
enddo



! Close to origin case.
! Crossover point is hurestically derived for n <= 30

!if (r2 < 2.28 * (n+7)) then
  value = 1

  rr2 = r2 / 4
  denom = 1.0_dp
  r2k = 1
  dvalo=1.d38
  do k = 1, nk_max
    r2k = r2k * rr2
    denom = denom * k * (n + k)
    dval = r2k / denom
    value = value + dval
    dvalo=dval
  if(done) then
     if(full_abs(dvalo)>=full_abs(dval)) exit
    else
    if (k>km0.and.full_abs(dval) < eps * full_abs(value)) then
      done=.true.
    endif
  endif
    if (k == nk_max) then
      print *, 'Internal error in norm_bessel_I: No convergence!'
      stop
    endif
  enddo
if(present(km)) write(6,*) k

  value = scale * value
call kill(r2, rr2, denom, r2k, r, scale, dval,dvalo)

    master=localmaster
end function etienne_bessel_It


function etienne_bessel_Itr(n, x, y,km) result (value)

implicit none

type(taylor) x
type(taylor)   r2, rr2, denom, r2k, r, scale, dval, value,dvalo
real(dp) eps,y
integer, optional :: km
integer n, k,km0
integer, parameter :: nk_max = 1000
logical done
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      value%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(value)
call alloc(r2, rr2, denom, r2k, r, scale, dval,dvalo)
eps=1.e-8_dp
done=.false.
r2 = x**2 + y**2
km0=15
if(present(km)) km0=km
scale = 1.0_dp/2**n !/ (2**n * factorial(n))


do k=1,n
scale=scale/k
enddo



! Close to origin case.
! Crossover point is hurestically derived for n <= 30

!if (r2 < 2.28 * (n+7)) then
  value = 1

  rr2 = r2 / 4
  denom = 1.0_dp
  r2k = 1
  dvalo=1.d38
  do k = 1, nk_max
    r2k = r2k * rr2
    denom = denom * k * (n + k)
    dval = r2k / denom
    value = value + dval
    dvalo=dval
  if(done) then
     if(full_abs(dvalo)>=full_abs(dval)) exit
    else
    if (k>km0.and.full_abs(dval) < eps * full_abs(value)) then
      done=.true.
    endif
  endif
    if (k == nk_max) then
      print *, 'Internal error in norm_bessel_I: No convergence!'
      stop
    endif
  enddo
if(present(km)) write(6,*) k

  value = scale * value
call kill(r2, rr2, denom, r2k, r, scale, dval,dvalo)

    master=localmaster
end function etienne_bessel_Itr


function etienne_bessel_Irt(n, x, y,km) result (value)

implicit none

type(taylor)  y
type(taylor)   r2, rr2, denom, r2k, r, scale, dval, value,dvalo
real(dp) x,eps
integer, optional :: km
integer n, k,km0
integer, parameter :: nk_max = 1000
logical done
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      value%i=0
     RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(value)
call alloc(r2, rr2, denom, r2k, r, scale, dval,dvalo)
eps=1.e-8_dp
done=.false.
r2 = x**2 + y**2
km0=15
if(present(km)) km0=km
scale = 1.0_dp/2**n !/ (2**n * factorial(n))


do k=1,n
scale=scale/k
enddo



! Close to origin case.
! Crossover point is hurestically derived for n <= 30

!if (r2 < 2.28 * (n+7)) then
  value = 1

  rr2 = r2 / 4
  denom = 1.0_dp
  r2k = 1
  dvalo=1.d38
  do k = 1, nk_max
    r2k = r2k * rr2
    denom = denom * k * (n + k)
    dval = r2k / denom
    value = value + dval
    dvalo=dval
  if(done) then
     if(full_abs(dvalo)>=full_abs(dval)) exit
    else
    if (k>km0.and.full_abs(dval) < eps * full_abs(value)) then
      done=.true.
    endif
  endif
    if (k == nk_max) then
      print *, 'Internal error in norm_bessel_I: No convergence!'
      stop
    endif
  enddo
if(present(km)) write(6,*) k

  value = scale * value
call kill(r2, rr2, denom, r2k, r, scale, dval,dvalo)

    master=localmaster
end function etienne_bessel_Irt

  !!! bessel   !!!!!!!!!!  2017

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!+
! Function norm_bessel_I(n, x, y) result (value)
!
! Routine to return the normalized bessel function defined to be
!   I_n(r) / r^n
! where I_n is the standard modified Bessel function of the first kind of order n and
!   r^2 = x^2 + y^2
!
! Input:
!   n     -- integer: Bessel order.
!   x, y  -- real(dp): Values to evaluate at.
!
! Output:
!   value -- real(rp): Normalized Bessel value.
!-

function norm_bessel_Ir(n, x, y) result (value)

implicit none

real(dp) x, y, value
real(dp) r2, rr2, denom, r2k, r, scale, dval

integer n, k
integer, parameter :: nk_max = 100

!

r2 = x**2 + y**2

scale = 1.0_dp / (2**n * factorial(n))

if (r2 == 0) then
  value = scale
  return
endif

if (n > 30) then
  print *, 'Error in norm_bessel_I: Bessel order greater than 30: ', n
  stop
endif

! Close to origin case.
! Crossover point is hurestically derived for n <= 30

if (r2 < 2.28 * (n+7)) then
  value = 1

  rr2 = r2 / 4
  denom = 1.0_dp
  r2k = 1

  do k = 1, nk_max
    r2k = r2k * rr2
    denom = denom * k * (n + k)
    dval = r2k / denom
    value = value + dval
    if (dval < 1d-16 * value) exit
    if (k == nk_max) then
      print *, 'Internal error in norm_bessel_I: No convergence!'
      stop
    endif
  enddo

  value = scale * value
endif

! Far from origin case

r = sqrt(r2)

select case (n)
case (0)
  value = bessel_i0(r)

case (1)
  value = bessel_i1(r) / r

case default
  value = bessel_I(n,r) / r**n
end select

end function norm_bessel_Ir

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------


function factorial(n) result (fact)

implicit none

real(8), parameter :: twopi=6.283185307179586476925286766559005768394_8
real(8) fact, lnn
real(8), parameter :: f(0:30) = [ &
      1.000000000000D+00, 1.000000000000D+00, 2.000000000000D+00, 6.000000000000D+00, &
      2.400000000000D+01, 1.200000000000D+02, 7.200000000000D+02, 5.040000000000D+03, &
      4.032000000000D+04, 3.628800000000D+05, 3.628800000000D+06, 3.991680000000D+07, &
      4.790016000000D+08, 6.227020800000D+09, 8.717829120000D+10, 1.307674368000D+12, &
      2.092278988800D+13, 3.556874280960D+14, 6.402373705728D+15, 1.216451004088D+17, &
      2.432902008177D+18, 5.109094217171D+19, 1.124000727778D+21, 2.585201673888D+22, &
      6.204484017332D+23, 1.551121004333D+25, 4.032914611266D+26, 1.088886945042D+28, &
      3.048883446117D+29, 8.841761993740D+30, 2.652528598122D+32]

integer n, i

! Use Sterling's formula if n is very large

if (n > ubound(f, 1)) then
  if (n > 170) stop
  fact = exp(n * log(real(n, 8)) - n + log(twopi * n) / 2)

else
  fact = f(n)
endif

end function

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

	FUNCTION bessel_I0(x)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessel_I0
	REAL(dp) :: ax
	REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
		3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
		0.45813e-2_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
		0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
		-0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
		0.392377e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessel_I0=poly_eval(real((x/3.75_dp)**2,dp),p)
	else
		bessel_I0=(exp(ax)/sqrt(ax))*poly_eval(real(3.75_dp/ax,dp),q)
	end if
	END FUNCTION bessel_I0

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

	FUNCTION bessel_I1(x)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessel_I1
	REAL(dp) :: ax
	REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
		0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
		0.301532e-2_dp,0.32411e-3_dp/)
	REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
		-0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
		0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
		-0.420059e-2_dp/)
	ax=abs(x)
	if (ax < 3.75) then
		bessel_I1=ax*poly_eval(real((x/3.75_dp)**2,dp),p)
	else
		bessel_I1=(exp(ax)/sqrt(ax))*poly_eval(real(3.75_dp/ax,dp),q)
	end if
	if (x < 0.0) bessel_I1=-bessel_I1
	END FUNCTION bessel_I1

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

	FUNCTION bessel_I(n,x)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: n
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessel_I
	INTEGER, PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER :: j,m
	REAL(dp) :: bi,bim,bip,tox
	if (n < 2) stop
	bessel_I=0.0
	if (x*x <= 8.0_dp*tiny(x)) RETURN
	tox=2.0_dp/abs(x)
	bip=0.0
	bi=1.0
	m=2*((n+int(sqrt(real(IACC*n,dp)))))
	do j=m,1,-1
		bim=bip+j*tox*bi
		bip=bi
		bi=bim
		if (exponent(bi) > IEXP) then
			bessel_I=scale(bessel_I,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end if
		if (j == n) bessel_I=bip
	end do
	bessel_I=bessel_I*bessel_I0(x)/bi
	if (x < 0.0 .and. mod(n,2) == 1) bessel_I=-bessel_I
	END FUNCTION bessel_I

!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------

	FUNCTION poly_eval(x,coeffs)
	REAL(dp), INTENT(IN) :: x
	REAL(dp), DIMENSION(:), INTENT(IN) :: coeffs
	REAL(dp) :: poly_eval
	REAL(dp) :: pow
	REAL(dp), DIMENSION(:), ALLOCATABLE :: vec
	INTEGER :: i,n,nn
	INTEGER, PARAMETER :: NPAR_POLY=8
	n=size(coeffs)
	if (n <= 0) then
		poly_eval=0.0_dp
	else if (n < NPAR_POLY) then
		poly_eval=coeffs(n)
		do i=n-1,1,-1
			poly_eval=x*poly_eval+coeffs(i)
		end do
	else
		allocate(vec(n+1))
		pow=x
		vec(1:n)=coeffs
		do
			vec(n+1)=0.0_dp
			nn=ishft(n+1,-1)
			vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
			if (nn == 1) exit
			pow=pow*pow
			n=nn
		end do
		poly_eval=vec(1)
		deallocate(vec)
	end if
	END FUNCTION poly_eval

  FUNCTION nbittaylortr( n,x,y )
    implicit none
    TYPE (taylor) nbittaylortr
    TYPE (taylor), INTENT (IN) :: x
    real(dp) y
    integer, INTENT (IN) :: n
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      nbittaylortr%i=0
     RETURN
    endif
    localmaster=master
    call ass(nbittaylortr)
    if(switch_bessel) then

     nbittaylortr=nbi_etienne(n,x,y)
    else

     nbittaylortr=nbi_david(n,x,y)

    endif
    master=localmaster
    END FUNCTION nbittaylortr

  FUNCTION nbittaylorrt( n,x,y )
    implicit none
    TYPE (taylor) nbittaylorrt
    TYPE (taylor), INTENT (IN) :: y
    real(dp) x
    integer, INTENT (IN) :: n
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      nbittaylorrt%i=0
     RETURN
    endif
    localmaster=master
    call ass(nbittaylorrt)
    if(switch_bessel) then

     nbittaylorrt=nbi_etienne(n,x,y)
    else

     nbittaylorrt=nbi_david(n,x,y)

    endif
    master=localmaster

    END FUNCTION nbittaylorrt

  FUNCTION nbittaylor( n,x,y )
    implicit none
    TYPE (taylor) nbittaylor
    TYPE (taylor), INTENT (IN) :: x,y
    integer, INTENT (IN) :: n
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
      nbittaylor%i=0
     RETURN
    endif
    localmaster=master
    call ass(nbittaylor)
    if(switch_bessel) then

     nbittaylor=nbi_etienne(n,x,y)
    else

     nbittaylor=nbi_david(n,x,y)

    endif
    master=localmaster
    END FUNCTION nbittaylor

  FUNCTION nbitreal( n,x,y )
    implicit none
    real(dp) nbitreal
    real(dp),  INTENT (IN) :: x,y
    integer, INTENT (IN) :: n

    if(switch_bessel) then

     nbitreal=nbi_etienne(n,x,y)
    else

     nbitreal=nbi_david(n,x,y)

    endif

    END FUNCTION nbitreal

  FUNCTION nbit( n,x,y )
    implicit none
    TYPE (taylor) nbit
    TYPE (taylor), INTENT (IN) :: x,y
    integer, INTENT (IN) :: n
    integer localmaster,i
    TYPE (taylor) t,dx,tn
    real(dp) x0,y0,fac,div
    IF(.NOT.C_%STABLE_DA) then
      nbit%i=0
     RETURN
    endif
    localmaster=master

    call ass(nbit)
    call alloc(dx,tn)

     x0=x
     y0=y
     dx=x**2+y**2-(x0**2+y0**2)
     tn=1.0_dp
     nbit=nbi_david(n,x0,y0)
     fac=1.0_dp
     do i=1,no
      tn=tn*dx
      fac=fac/i/2.0_dp
      div=nbi_david(n+i,x0,y0)*fac
      nbit=nbit+tn*div
     enddo

    call kill(dx,tn)
    master=localmaster

  END FUNCTION nbit

FUNCTION nbitrt( n,xx,y )
    implicit none
    TYPE (taylor) nbitrt
    TYPE (taylor), INTENT (IN) :: y
    real(dp), INTENT (IN) :: xx
    integer, INTENT (IN) :: n
    integer localmaster,i
    TYPE (taylor) x

    IF(.NOT.C_%STABLE_DA) then
      nbitrt%i=0
     RETURN
    endif
    localmaster=master

    call ass(nbitrt)
    call alloc(x)
      x=xx
      nbitrt=nbit( n,x,y )

    call kill(x)
    master=localmaster

  END FUNCTION nbitrt

FUNCTION nbittr( n,x,yy )
    implicit none
    TYPE (taylor) nbittr
    TYPE (taylor), INTENT (IN) :: x
    real(dp), INTENT (IN) :: yy
    integer, INTENT (IN) :: n
    integer localmaster,i
    TYPE (taylor) y

    IF(.NOT.C_%STABLE_DA) then
      nbittr%i=0
     RETURN
    endif
    localmaster=master

    call ass(nbittr)
    call alloc(y)
      y=yy
      nbittr=nbit( n,x,y )

    call kill(y)
    master=localmaster

  END FUNCTION nbittr

END MODULE  tpsa
