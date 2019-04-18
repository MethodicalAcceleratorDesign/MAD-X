!The Full Polymorphic Package
!Copyright (C) Etienne Forest
! Taylor polymorphism at execution is based on an idea
! and C++ prototype developed  by J. Bengtsson circa 1990

module polymorphic_taylor
  use complex_taylor
  implicit none
  public
  logical(lp),private,parameter::t=.true.,f=.false.
  integer,private::NO,ND,ND2,NP,NDPT,NV          !,lastmaster   ! 2002.12.13
  INTEGER,PRIVATE::NMAX_pol=100
  !  real(dp),PRIVATE::EPS=c_1d_6
  ! integer ent,exi
  integer,private,parameter::m1=mmmmmm1,m2=mmmmmm2,m3=mmmmmm3,ms=mmmmmm4
  integer,private,parameter:: m11=m1+ms*m1,m12=m1+ms*m2,m13=m1+ms*m3,  &
       m21=m2+ms*m1,m22=m2+ms*m2,m23=m2+ms*m3,                         &
       m31=m3+ms*m1,m32=m3+ms*m2,m33=m3+ms*m3
  logical(lp),private::old
  private resetpoly,resetpolyn ,resetpoly0,resetpolyn0 ,allocpoly,allocpolyn,resetpoly_R,resetpoly_Rn
  PRIVATE K_OPT,A_OPT
  private equal,Dequaldacon,equaldacon,iequaldacon ,realEQUAL,singleequal
  private taylorEQUAL,EQUALtaylor,complexreal_8
  private real_8univ,univreal_8
  PRIVATE real_8MAP,MAPREAL_8
  private unaryADD,add,daddsc,dscadd,addsc,scadd,iaddsc,iscadd
  private unarySUB,subs,dscsub,dsubsc,scsub,subsc,iscsub,isubsc
  private mul,dmulsc,dscmul,mulsc,scmul,imulsc,iscmul
  private div,ddivsc,dscdiv,divsc,scdiv,idivsc,iscdiv
  private POW,POWR,POWR8,POWq
  private dexpt,dcost,dcosdt,dsindt,dsint,dlogt,dsqrtt
  private greaterthan,greatereq
  private lessthan,dlessthansc,dsclessthan,lessthansc,sclessthan,ilessthansc,isclessthan
  private igreatersc,iscgreater,dgreatersc,dscgreater,greatersc,scgreater
  private dgreatereqsc,dscgreatereq,greatereqsc,scgreatereq,igreatereqsc,iscgreatereq
  private abst,full_abst,print6
  private lesseq,dlesseqsc,dsclesseq,lesseqsc,sclesseq,ilesseqsc,isclesseq
  private eq,deqsc,dsceq,eqsc,sceq,ieqsc,isceq
  private neq,dneqsc,dscneq,neqsc,scneq,ineqsc,iscneq
  PRIVATE GETCHARnd2,GETintnd2,getchar,GETint,GETORDER,CUTORDER
  private  normal_p !,beamENV_8
  private absoftdatan2dr,absoftdcosdr,absoftdsindr,absoftdtandr,absoftdatandr
  !complex stuff
  private datant,datanDt,datan2t,dasint,dacost,dtant,dtandt,datanht,datan2tt
  private dcosht,dsinht,dtanht,SINX_XT,SINHX_XT,polymorpht
  ! PRIVATE EQUAL1D,EQUAL2D
  ! end complex stuff
  private printpoly,printdouble,printsingle,dmulmapconcat,nbiP
  private line,Mequaldacon,cos_quaternionr,cos_quaternionp,sin_quaternionr,sin_quaternionp
  private make_it_knobr,kill_knobr
  character(120) line
  !integer npol
  !parameter (npol=20)
  !INTEGER iasspol
  !type(real_8)  DUMMYpol,DUMpol(ndum)
  !integer iasspol0
  !private assp
  PRIVATE SINH_HR
  PRIVATE SIN_HR
  ! PROBE_8 STUFF
  real(dp) :: sinhx_x_min=1e-4_dp
  real(dp) :: sinhx_x_minp=1.0_dp  !  1.e-9  !c_0_0001
 
INTEGER, private, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  private allocquaternion,allocquaternionn,A_OPT_quaternion
  private k_opt_quaternion,killquaternionn,killquaternion
  private equalq,mulq,subq,addq,EQUALq_8_r,EQUALq_r_8,printpolyq,absq,absq2,invq,EQUALq_r

 

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL   ! 2002.10.9
     MODULE PROCEDURE EQUALq 
     MODULE PROCEDURE EQUALq_8_r   ! 2002.10.9
     MODULE PROCEDURE EQUALq_r_8
     MODULE PROCEDURE EQUALq_r
     !     MODULE PROCEDURE EQUAL1D  ! 2004.7.10
     !     MODULE PROCEDURE EQUAL2D  ! 2004.7.10
     MODULE PROCEDURE complexreal_8
     !zgoubi
     MODULE PROCEDURE realEQUAL   !
     MODULE PROCEDURE singleequal
     MODULE PROCEDURE taylorEQUAL    !taylor=real_8
     MODULE PROCEDURE EQUALtaylor    !real_8= taylor    ! here 2002.10.9
     MODULE PROCEDURE Mequaldacon   ! dequaldacon on array
     MODULE PROCEDURE Dequaldacon  !
     MODULE PROCEDURE equaldacon   !
     ! For resetting
     MODULE PROCEDURE Iequaldacon
     MODULE PROCEDURE Iequaldaconn
     ! For universal_taylor
     MODULE PROCEDURE real_8univ
     MODULE PROCEDURE univreal_8   ! new sagan
     ! for damap
     MODULE PROCEDURE MAPreal_8
     MODULE PROCEDURE real_8MAP


     !  allowing normalization of polymorphic arrays directly
!     MODULE PROCEDURE beamENV_8
     MODULE PROCEDURE normal_p

  end  INTERFACE


  !@   <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber2" width="400" height="135">
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">+</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">
  !@         Real(dp)</font></span></td>
  !@         <td width="78" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@         <td width="56" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase; font-weight:700">
  !@         <font face="Times New Roman" size="1">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#ADD">add</a></font></span></td>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase; font-weight:700">
  !@         <font face="Times New Roman" size="1">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#DADDSC">daddsc</a></font></span></td>
  !@         <td width="78" height="20" align="center"><b>
  !@         <font size="1" face="Times New Roman">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#ADDSC">ADDSC</a></font></b></td>
  !@         <td width="56" height="20" align="center"><b>
  !@         <font size="1" face="Times New Roman">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#IADDSC">IADDSC</a></font></b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">
  !@         Real(dp)</font></span></td>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase; font-weight:700">
  !@         <font face="Times New Roman" size="1">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#DSCADD">dscadd</a></font></span></td>
  !@         <td width="77" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@         <td width="78" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@         <td width="56" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@         <td width="77" height="20" align="center"><b>
  !@         <font size="1" face="Times New Roman">
  !@         <a style="text-decoration: none" href="m_real_polymorph#SCADD">SCADD</a></font></b></td>
  !@         <td width="77" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@         <td width="78" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@         <td width="56" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">Integer</font></span></td>
  !@         <td width="77" height="20" align="center"><b>
  !@         <font size="1" face="Times New Roman">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#ISCADD">ISCADD</a></font></b></td>
  !@         <td width="77" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@         <td width="78" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@         <td width="56" height="20" align="center">
  !@         <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@       </tr>
  !@     </table>

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryADD  !
     MODULE PROCEDURE add   !
     MODULE PROCEDURE addq   
     MODULE PROCEDURE daddsc  !
     MODULE PROCEDURE dscadd  !
     MODULE PROCEDURE addsc  !
     MODULE PROCEDURE scadd  !
     MODULE PROCEDURE iaddsc  !
     MODULE PROCEDURE iscadd  !
  END INTERFACE

  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="400" height="135">
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">-</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="78" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@        <td width="56" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Integer</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#SUBS">SUBS</a></font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#DSUBSC">dSUBsc</a></font></span></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#SUBSC">SUBSC</a></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#ISUBSC">ISUBSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#DSCSUB">dscSUB</a></font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#SCSUB">SCSUB</a></font></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Integer</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#ISCSUB">ISCSUB</a></font></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@    </table>

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE unarySUB   !
     MODULE PROCEDURE subs    !
     MODULE PROCEDURE subq   !
     MODULE PROCEDURE dscsub   !
     MODULE PROCEDURE dsubsc   !
     MODULE PROCEDURE subsc   !
     MODULE PROCEDURE scsub   !
     MODULE PROCEDURE isubsc   !
     MODULE PROCEDURE iscsub   !
  END INTERFACE

  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="400" height="134">
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">*</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="78" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@        <td width="56" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Integer</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#MUL">MUL</a></font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#DMULSC">dMULsc</a></font></span></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#MULSC">MULSC</a></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#IMULSC">IMULSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="19" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="77" height="19" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#DSCMUL">dscMUL</a></font></span></td>
  !@        <td width="77" height="19" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="19" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="19" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#SCMUL">SCMUL</a></font></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Integer</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a style="text-decoration: none; font-weight:700" href="m_real_polymorph.htm#ISCMUL">ISCMUL</a></font></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@    </table>

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul    !
     MODULE PROCEDURE mulq   !
     MODULE PROCEDURE dmulsc   !
     MODULE PROCEDURE dscmul   !
     MODULE PROCEDURE mulsc   !
     MODULE PROCEDURE scmul   !
     MODULE PROCEDURE imulsc   !
     MODULE PROCEDURE iscmul   !
     MODULE PROCEDURE      dmulmapconcat  ! concatenation real_8*damap
  END INTERFACE

  !@       <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="400" height="135">
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">/</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Integer</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">REAL_8</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#DIV">div</a></font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#DDIVSC">dDIVsc</a></font></span></td>
  !@        <td width="0" height="20" align="center"><font size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#DIVSC">DIVSC</a></font></td>
  !@        <td width="0" height="20" align="center"><font size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#IDIVSC">IDIVSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#DSCDIV">dscDIV</a></font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Real(sp)</font></span></td>
  !@        <td width="0" height="20" align="center"><font size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#SCDIV">SCDIV</a></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Integer</font></span></td>
  !@        <td width="0" height="20" align="center"><font size="1">
  !@        <a style="text-decoration: none; font-weight: 700" href="m_real_polymorph.htm#ISCDIV">ISCDIV</a></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@    </table>



  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div    !
     MODULE PROCEDURE ddivsc   !
     MODULE PROCEDURE dscdiv   !
     MODULE PROCEDURE divsc   !
     MODULE PROCEDURE scdiv   !
     MODULE PROCEDURE idivsc   !
     MODULE PROCEDURE iscdiv   !
  END INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW     !
     MODULE PROCEDURE POWq     !
     MODULE PROCEDURE POWR     !
     MODULE PROCEDURE POWR8    !
  END INTERFACE

  !  Logical Operators

  INTERFACE OPERATOR (<)
     MODULE PROCEDURE lessthan
     MODULE PROCEDURE dlessthansc
     MODULE PROCEDURE dsclessthan
     MODULE PROCEDURE lessthansc
     MODULE PROCEDURE sclessthan
     MODULE PROCEDURE ilessthansc
     MODULE PROCEDURE isclessthan
  END INTERFACE

  INTERFACE OPERATOR (>)
     MODULE PROCEDURE greaterthan
     MODULE PROCEDURE dgreatersc
     MODULE PROCEDURE dscgreater
     MODULE PROCEDURE greatersc
     MODULE PROCEDURE scgreater
     MODULE PROCEDURE igreatersc
     MODULE PROCEDURE iscgreater
  END INTERFACE

  INTERFACE OPERATOR (>=)
     MODULE PROCEDURE greatereq
     MODULE PROCEDURE dgreatereqsc
     MODULE PROCEDURE dscgreatereq
     MODULE PROCEDURE greatereqsc
     MODULE PROCEDURE scgreatereq
     MODULE PROCEDURE igreatereqsc
     MODULE PROCEDURE iscgreatereq
  END INTERFACE

  INTERFACE OPERATOR (<=)
     MODULE PROCEDURE lesseq
     MODULE PROCEDURE dlesseqsc
     MODULE PROCEDURE dsclesseq
     MODULE PROCEDURE lesseqsc
     MODULE PROCEDURE sclesseq
     MODULE PROCEDURE ilesseqsc
     MODULE PROCEDURE isclesseq
  END INTERFACE

  INTERFACE OPERATOR (==)
     MODULE PROCEDURE eq
     MODULE PROCEDURE deqsc
     MODULE PROCEDURE dsceq
     MODULE PROCEDURE eqsc
     MODULE PROCEDURE sceq
     MODULE PROCEDURE ieqsc
     MODULE PROCEDURE isceq
  END INTERFACE

  INTERFACE OPERATOR (/=)
     MODULE PROCEDURE neq
     MODULE PROCEDURE dneqsc
     MODULE PROCEDURE dscneq
     MODULE PROCEDURE neqsc
     MODULE PROCEDURE scneq
     MODULE PROCEDURE ineqsc
     MODULE PROCEDURE iscneq
  END INTERFACE



  ! New Operators

  INTERFACE OPERATOR (.SUB.)
     MODULE PROCEDURE GETORDER
     MODULE PROCEDURE getchar
     MODULE PROCEDURE GETint
  END INTERFACE

  INTERFACE OPERATOR (.CUT.)
     MODULE PROCEDURE CUTORDER
  END INTERFACE

  INTERFACE OPERATOR (.PAR.)
     MODULE PROCEDURE getcharnd2
     MODULE PROCEDURE GETintnd2
  END INTERFACE

  ! Intrinsic function  (some are Fortran extension not always supported)

  INTERFACE exp
     MODULE PROCEDURE dexpt
  END INTERFACE

  INTERFACE invqq
     MODULE PROCEDURE invq
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
  INTERFACE ccos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE cdcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE dcos
     MODULE PROCEDURE dcost
  END INTERFACE

  INTERFACE atanh
     MODULE PROCEDURE datanht
  END INTERFACE

  INTERFACE tan
     MODULE PROCEDURE dtant
  END INTERFACE


  INTERFACE dtan
     MODULE PROCEDURE dtant
  END INTERFACE

  INTERFACE tand
     MODULE PROCEDURE dtandt
     !     MODULE PROCEDURE absoftdtandr  ! for non PC
  END INTERFACE
  INTERFACE dtand
     MODULE PROCEDURE dtandt
     !    MODULE PROCEDURE absoftdtandr   ! for non PC absoft
  END INTERFACE


  INTERFACE cosd
     MODULE PROCEDURE dcosdt
     !     MODULE PROCEDURE absoftdcosdr  ! for non PC
  END INTERFACE
  INTERFACE dcosd
     MODULE PROCEDURE dcosdt
     !    MODULE PROCEDURE absoftdcosdr  ! for non PC absoft
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

  INTERFACE sind
     MODULE PROCEDURE dsindt
     !     MODULE PROCEDURE absoftdsindr  ! for non PC
  END INTERFACE
  INTERFACE dsind
     MODULE PROCEDURE dsindt
     !    MODULE PROCEDURE absoftdsindr  ! for non PC absoft
  END INTERFACE

 INTERFACE norm_bessel_I
     MODULE PROCEDURE nbiP
  END INTERFACE

 INTERFACE nbi
     MODULE PROCEDURE nbiP
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

  INTERFACE abs
     MODULE PROCEDURE abst
     MODULE PROCEDURE absq
  END INTERFACE


  INTERFACE abs_square
     MODULE PROCEDURE absq2
  END INTERFACE

  INTERFACE dabs
     MODULE PROCEDURE abst
  END INTERFACE

  ! complex stuff
  INTERFACE atan
     MODULE PROCEDURE datant
  END INTERFACE
  INTERFACE datan
     MODULE PROCEDURE datant
  END INTERFACE

  INTERFACE atan2
     MODULE PROCEDURE datan2t
     MODULE PROCEDURE datan2tt
  END INTERFACE
  INTERFACE datan2
     MODULE PROCEDURE datan2t
     MODULE PROCEDURE datan2tt
  END INTERFACE

  INTERFACE atan2d
     MODULE PROCEDURE datan2dt
     !   MODULE PROCEDURE absoftdatan2dr      ! for non PC
  END INTERFACE
  INTERFACE datan2d
     MODULE PROCEDURE datan2dt
     !    MODULE PROCEDURE absoftdatan2dr      ! for non PC absoft
  END INTERFACE

  INTERFACE atand
     MODULE PROCEDURE datandt
     !     MODULE PROCEDURE absoftdatandr   ! for non PC
  END INTERFACE
  INTERFACE datand
     MODULE PROCEDURE datandt
     !    MODULE PROCEDURE absoftdatandr    ! for non PC  absoft
  END INTERFACE

  INTERFACE asin
     MODULE PROCEDURE dasint
  END INTERFACE
  INTERFACE dasin
     MODULE PROCEDURE dasint
  END INTERFACE

  INTERFACE acos
     MODULE PROCEDURE dacost
  END INTERFACE
  INTERFACE dacos
     MODULE PROCEDURE dacost
  END INTERFACE

  INTERFACE cosh
     MODULE PROCEDURE dcosht
  END INTERFACE
  INTERFACE dcosh
     MODULE PROCEDURE dcosht
  END INTERFACE

  INTERFACE sinh
     MODULE PROCEDURE dsinht
  END INTERFACE
  INTERFACE dsinh
     MODULE PROCEDURE dsinht
  END INTERFACE

  INTERFACE tanh
     MODULE PROCEDURE dtanht
  END INTERFACE
  INTERFACE dtanh
     MODULE PROCEDURE dtanht
  END INTERFACE

  ! end Intrinsic function

  ! Non Intrinsic function

  INTERFACE full_abs
     MODULE PROCEDURE full_abst
  END INTERFACE

  INTERFACE morph
     MODULE PROCEDURE polymorpht
  END INTERFACE

  INTERFACE SINHX_X
     MODULE PROCEDURE SINH_HR
     MODULE PROCEDURE SINHX_XT
  END INTERFACE

  INTERFACE SINX_X
     MODULE PROCEDURE SIN_HR
     MODULE PROCEDURE SINX_XT
  END INTERFACE

  INTERFACE cos_quaternion
     MODULE PROCEDURE cos_quaternionr
     MODULE PROCEDURE cos_quaternionp
   end INTERFACE 

  INTERFACE sin_quaternion
     MODULE PROCEDURE sin_quaternionr
     MODULE PROCEDURE sin_quaternionp
   end INTERFACE 
  ! i/o

  INTERFACE daprint
     MODULE PROCEDURE printpoly !
     MODULE PROCEDURE printdouble
     MODULE PROCEDURE printsingle
     MODULE PROCEDURE printpolyq
     MODULE PROCEDURE print6
  END INTERFACE
  INTERFACE print
     MODULE PROCEDURE printpoly   !
     MODULE PROCEDURE printdouble
     MODULE PROCEDURE printsingle
     MODULE PROCEDURE printpolyq
     MODULE PROCEDURE print6
  END INTERFACE



  ! constructors and destructors

  INTERFACE alloc
     MODULE PROCEDURE A_OPT    !allocpoly   !
     MODULE PROCEDURE allocpolyn  !
     MODULE PROCEDURE  allocquaternionn
     MODULE PROCEDURE  A_OPT_quaternion

  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE k_opt
     MODULE PROCEDURE killquaternionn
     MODULE PROCEDURE k_opt_quaternion
     MODULE PROCEDURE resetpolyn0  !
     MODULE PROCEDURE resetpoly_R
     MODULE PROCEDURE resetpoly_RN
  END INTERFACE

  INTERFACE reset
     MODULE PROCEDURE resetpoly   !
     MODULE PROCEDURE resetpolyn  !
  END INTERFACE

interface make_it_knob
 module procedure make_it_knobr
end INTERFACE

interface kill_knob
module procedure kill_knobr
end INTERFACE


  ! Managing

  INTERFACE ass
     MODULE PROCEDURE assp
  END INTERFACE

contains

  subroutine make_it_knobr(k,i,s)
    implicit none
    TYPE (real_8), intent(inout) :: k
    real(dp), optional :: s
    integer, intent(in) :: i
    if(i==0) return
    k%s=1.0_dp
    if(present(s)) k%s=s
    k%i=i
    k%kind=3
  end subroutine make_it_knobr

  subroutine kill_knobr(k)
    implicit none
    TYPE (real_8), intent(inout) :: k
    k%s=1.0_dp
    k%i=0
    k%kind=1
  end subroutine kill_knobr


  FUNCTION polymorpht( S1 )
    implicit none
    TYPE (real_8) polymorpht
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster

    localmaster=master
    call ass(polymorpht)
    polymorpht%t= s1
    master=localmaster

  END FUNCTION polymorpht


  FUNCTION GETCHARnd2( S1, S2 )
    implicit none
    TYPE (real_8) GETCHARnd2
    TYPE (real_8), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    integer localmaster
    type(taylor) t

    localmaster=master
    call ass(GETCHARnd2)
    call alloc(t)
    t=s1

    t=t.par.s2
    GETCHARnd2%t=t



    !       GETCHARnd2%kind=m2
    !    if(s1%kind==m2) then
    !       GETCHARnd2%t=s1%t.par.s2
    !    else
    !       GETCHARnd2%t=zero
    !    endif

    call kill(t)
    master=localmaster

  END FUNCTION GETCHARnd2

  FUNCTION GETintnd2( S1, S2 )
    implicit none
    TYPE (real_8) GETintnd2
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) ::  S2(:)
    integer localmaster
    type(taylor) t

    localmaster=master
    call ass(GETintnd2)
    call alloc(t)
    t=s1
    !  if(s1%kind==m2) then
    !      GETintnd2%t=s1%t.par.s2
    !   else
    !       GETintnd2%kind=m2
    !       GETintnd2%t=zero
    !
    !  endif
    t=t.par.s2
    GETintnd2%t=t


    call kill(t)

    master=localmaster

  END FUNCTION GETintnd2

  FUNCTION GETchar( S1, S2 )
    implicit none
    real(dp) GETchar
    TYPE (real_8), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    !  integer localmaster

    GETchar=0.0_dp
    if(s1%kind==m2) then
       ! GETchar%t=s1%t.sub.s2   !  OLD
       GETchar=s1%t.sub.s2   !  CHANGE

    endif

  END FUNCTION GETchar

  FUNCTION GETint( S1, S2 )
    implicit none
    real(dp) GETint
    TYPE (real_8), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    integer i

    GETint=0.0_dp
    if(s1%kind==m2) then
       GETint=s1%t.sub.s2   !  CHANGE
    elseif(s1%kind==m1) then
       !zgoubi
       GETint=s1
       do i=1,size(s2)
          if(S2(i)/=0) then
             GETint=0.0_dp
             exit
          endif
       enddo

    endif

  END FUNCTION GETint

  FUNCTION GETORDER( S1, S2 )
    implicit none
    TYPE (real_8) GETORDER
    TYPE (real_8), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer localmaster
    localmaster=master
    call ass(GETORDER)

    if(s1%kind==m2) then
       GETORDER%t=s1%t.sub.s2
    else
       GETORDER%kind=m1
       GETORDER%r=0.0_dp
       if(s2==0) GETORDER%r=s1%r
    endif

    master=localmaster

  END FUNCTION GETORDER




  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (real_8) CUTORDER
    TYPE (real_8), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2
    integer localmaster

    localmaster=master
    call ass(CUTORDER)

    CUTORDER=0.0_dp
    if(s1%kind==m2) then
       cutorder%kind=m2
       CUTORDER%t=s1%t.CUT.s2
    elseif(s1%kind==m1) then
       if(s2>=1) CUTORDER%r=s1%r
       cutorder%kind=m1
    endif

    master=localmaster

  END FUNCTION CUTORDER




  FUNCTION greatereq( S1, S2 )
    implicit none
    logical(lp) greatereq
    TYPE (real_8), INTENT (IN) :: S1,S2

    greatereq=.FALSE.
    select case(s1%kind+ms*s2%kind)
    case(m11)
       IF(S1%r>=S2%r) greatereq=.TRUE.
    case(m12,m21,m22)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          IF((S1%t.SUB.'0')>=S2%r) greatereq=.TRUE.
       case(m12)
          IF(S1%r>=(S2%t.SUB.'0')) greatereq=.TRUE.
       case(m22)
          IF((S1%t.SUB.'0')>=(S2%t.SUB.'0')) greatereq=.TRUE.
       end select
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31,m13,m33)
          IF(S1%r>=S2%r) greatereq=.TRUE.
       case(m32)
          IF(S1%r>=(S2%t.SUB.'0')) greatereq=.TRUE.
       case(m23)
          IF((S1%t.SUB.'0')>=S2%r) greatereq=.TRUE.
       end select
    case default
  
         write(6,*) " trouble in greatereq "
         write(6,*) "s1%kind ,s2%kind ",s1%kind ,s2%kind 
 
       ! call !write_e(0)
    end select

  END FUNCTION greatereq

  FUNCTION dgreatereqsc( S1, S2 )
    implicit none
    logical(lp) dgreatereqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dgreatereqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r>=S2) dgreatereqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')>=S2) dgreatereqsc=.TRUE.
    case default
 
        write(6,*) " trouble in dgreatereqsc "
        write(6,*) "s1%kind ",s1%kind
 
       ! call !write_e(0)
    end select

  END FUNCTION dgreatereqsc

  FUNCTION dscgreatereq( S2,S1  )
    implicit none
    logical(lp) dscgreatereq
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dscgreatereq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2>=S1%r) dscgreatereq=.TRUE.
    case(m2)
       IF(S2>=(S1%t.SUB.'0')) dscgreatereq=.TRUE.
    case default

 
      write(6,*)" trouble in dscgreatereq "
      write(6,*)"s1%kind ",s1%kind 
 
    end select

  END FUNCTION dscgreatereq

  FUNCTION greatereqsc( S1, S2 )
    implicit none
    logical(lp) greatereqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    greatereqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r>=S2) greatereqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')>=S2) greatereqsc=.TRUE.
    case default
 
       write(6,*) " trouble in greatereqsc "
       write(6,*)"s1%kind ",s1%kind
 
    end select

  END FUNCTION greatereqsc

  FUNCTION scgreatereq( S2,S1  )
    implicit none
    logical(lp) scgreatereq
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    scgreatereq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2>=S1%r) scgreatereq=.TRUE.
    case(m2)
       IF(S2>=(S1%t.SUB.'0')) scgreatereq=.TRUE.
    case default
 
       write(6,*)" trouble in scgreatereq "
       write(6,*)"s1%kind ",s1%kind
   
    end select

  END FUNCTION scgreatereq

  FUNCTION igreatereqsc( S1, S2 )
    implicit none
    logical(lp) igreatereqsc
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    igreatereqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r>=S2) igreatereqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')>=S2) igreatereqsc=.TRUE.
    case default
 
       write(6,*)" trouble in igreatereqsc "
       write(6,*)"s1%kind ",s1%kind
 
    end select

  END FUNCTION igreatereqsc

  FUNCTION iscgreatereq( S2,S1  )
    implicit none
    logical(lp) iscgreatereq
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    iscgreatereq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2>=S1%r) iscgreatereq=.TRUE.
    case(m2)
       IF(S2>=(S1%t.SUB.'0')) iscgreatereq=.TRUE.
    case default

        write(6,*)" trouble in iscgreatereq "
        write(6,*)"s1%kind ",s1%kind
 
       ! call !write_e(0)
    end select

  END FUNCTION iscgreatereq

  FUNCTION greaterthan( S1, S2 )
    implicit none
    logical(lp) greaterthan
    TYPE (real_8), INTENT (IN) :: S1,S2

    greaterthan=.FALSE.
    select case(s1%kind+ms*s2%kind)
    case(m11)
       IF(S1%r>S2%r) greaterthan=.TRUE.
    case(m12,m21,m22)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          IF((S1%t.SUB.'0')>S2%r) greaterthan=.TRUE.
       case(m12)
          IF(S1%r>(S2%t.SUB.'0')) greaterthan=.TRUE.
       case(m22)
          IF((S1%t.SUB.'0')>(S2%t.SUB.'0')) greaterthan=.TRUE.
       end select
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31,m13,m33)
          IF(S1%r>S2%r) greaterthan=.TRUE.
       case(m32)
          IF(S1%r>(S2%t.SUB.'0')) greaterthan=.TRUE.
       case(m23)
          IF((S1%t.SUB.'0')>S2%r) greaterthan=.TRUE.
       end select
    case default
 
 
       write(6,*) "s1%kind ,s2%kind ",s1%kind ,s2%kind 
 
       ! call !write_e(0)
    end select

  END FUNCTION greaterthan

  FUNCTION igreatersc( S1, S2 )
    implicit none
    logical(lp) igreatersc
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    igreatersc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r>S2) igreatersc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')>S2) igreatersc=.TRUE.
    case default
 
       write(6,*) " trouble in igreatersc "
       write(6,*)  "s1%kind ",s1%kind
 
 
    end select

  END FUNCTION igreatersc

  FUNCTION iscgreater( S2,S1  )
    implicit none
    logical(lp) iscgreater
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    iscgreater=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2>S1%r) iscgreater=.TRUE.
    case(m2)
       IF(S2>(S1%t.SUB.'0')) iscgreater=.TRUE.
    case default
 
       write(6,*)" trouble in iscgreater "
        write(6,*)"s1%kind ",s1%kind

 
    end select

  END FUNCTION iscgreater

  FUNCTION dgreatersc( S1, S2 )
    implicit none
    logical(lp) dgreatersc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dgreatersc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r>S2) dgreatersc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')>S2) dgreatersc=.TRUE.
    case default
 
      write(6,*) " trouble in dgreatersc "
       write(6,*)"s1%kind ",s1%kind
 
    end select

  END FUNCTION dgreatersc

  FUNCTION dscgreater( S2,S1  )
    implicit none
    integer ipause, mypause
    logical(lp) dscgreater
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dscgreater=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2>S1%r) dscgreater=.TRUE.
    case(m2)
       IF(S2>(S1%t.SUB.'0')) dscgreater=.TRUE.
    case default
 
           write(6,*)" trouble in dscgreater "
          write(6,*) "s1%kind ",s1%kind
 
    end select

  END FUNCTION dscgreater

  FUNCTION greatersc( S1, S2 )
    implicit none
    logical(lp) greatersc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    greatersc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r>S2) greatersc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')>S2) greatersc=.TRUE.
    case default
 
        write(6,*)" trouble in greatersc "
        write(6,*)"s1%kind   ",s1%kind
 
    end select

  END FUNCTION greatersc

  FUNCTION scgreater( S2,S1  )
    implicit none
    logical(lp) scgreater
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    scgreater=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2>S1%r) scgreater=.TRUE.
    case(m2)
       IF(S2>(S1%t.SUB.'0')) scgreater=.TRUE.
    case default
 
       write(6,*)" trouble in scgreater "
       write(6,*)"s1%kind   ",s1%kind 
 
       ! call !write_e(0)
    end select

  END FUNCTION scgreater

  FUNCTION lessthan( S1, S2 )
    implicit none
    logical(lp) lessthan
    TYPE (real_8), INTENT (IN) :: S1,S2

    lessthan=.FALSE.
    select case(s1%kind+ms*s2%kind)
    case(m11)
       IF(S1%r<S2%r) lessthan=.TRUE.
    case(m12,m21,m22)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          IF((S1%t.SUB.'0')<S2%r) lessthan=.TRUE.
       case(m12)
          IF(S1%r<(S2%t.SUB.'0')) lessthan=.TRUE.
       case(m22)
          IF((S1%t.SUB.'0')<(S2%t.SUB.'0')) lessthan=.TRUE.
       end select
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31,m13,m33)
          IF(S1%r<S2%r) lessthan=.TRUE.
       case(m32)
          IF(S1%r<(S2%t.SUB.'0')) lessthan=.TRUE.
       case(m23)
          IF((S1%t.SUB.'0')<S2%r) lessthan=.TRUE.
       end select
    case default
 
         write(6,*)" trouble in lessthan "
         write(6,*) "s1%kind ,s2%kind ",s1%kind ,s2%kind
 
       ! call !write_e(0)
    end select

  END FUNCTION lessthan

  FUNCTION dlessthansc( S1, S2 )
    implicit none
    logical(lp) dlessthansc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dlessthansc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r<S2) dlessthansc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')<S2) dlessthansc=.TRUE.
    case default

        write(6,*) " trouble in dlessthansc "
        write(6,*) "s1%kind   ",s1%kind
 
    end select

  END FUNCTION dlessthansc


  FUNCTION dsclessthan(  S2,S1 )
    implicit none
    logical(lp) dsclessthan
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dsclessthan=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2<S1%r) dsclessthan=.TRUE.
    case(m2)
       IF(S2<(S1%t.SUB.'0')) dsclessthan=.TRUE.
    case default

       write(6,*)" trouble in dsclessthan "
      write(6,*)"s1%kind   ",s1%kind 
 
    end select

  END FUNCTION dsclessthan

  FUNCTION lessthansc( S1, S2 )
    implicit none
    logical(lp) lessthansc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    lessthansc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r<S2) lessthansc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')<S2) lessthansc=.TRUE.
    case default
       !w_p=0

       write(6,*)" trouble in lessthansc "
       write(6,*)"s1%kind   ",s1%kind
 
    end select

  END FUNCTION lessthansc


  FUNCTION sclessthan(  S2,S1 )
    implicit none
    logical(lp) sclessthan
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    sclessthan=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2<S1%r) sclessthan=.TRUE.
    case(m2)
       IF(S2<(S1%t.SUB.'0')) sclessthan=.TRUE.
    case default
 
        write(6,*) " trouble in sclessthan "
        write(6,*)"s1%kind   ",s1%kind
  
    end select

  END FUNCTION sclessthan

  FUNCTION ilessthansc( S1, S2 )
    implicit none
    logical(lp) ilessthansc
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    ilessthansc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r<S2) ilessthansc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')<S2) ilessthansc=.TRUE.
    case default

       write(6,*) " trouble in ilessthansc "
       write(6,*) "s1%kind ,s2  ",s1%kind ,s2 
 
 
    end select

  END FUNCTION ilessthansc


  FUNCTION isclessthan(  S2,S1 )
    implicit none
    logical(lp) isclessthan
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    isclessthan=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2<S1%r) isclessthan=.TRUE.
    case(m2)
       IF(S2<(S1%t.SUB.'0')) isclessthan=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in isclessthan "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION isclessthan

  FUNCTION lesseq( S1, S2 )
    implicit none
    logical(lp) lesseq
    TYPE (real_8), INTENT (IN) :: S1,S2

    lesseq=.FALSE.
    select case(s1%kind+ms*s2%kind)
    case(m11)
       IF(S1%r<=S2%r) lesseq=.TRUE.
    case(m12,m21,m22)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          IF((S1%t.SUB.'0')<=S2%r) lesseq=.TRUE.
       case(m12)
          IF(S1%r<=(S2%t.SUB.'0')) lesseq=.TRUE.
       case(m22)
          IF((S1%t.SUB.'0')<=(S2%t.SUB.'0')) lesseq=.TRUE.
       end select
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31,m13,m33)
          IF(S1%r<=S2%r) lesseq=.TRUE.
       case(m32)
          IF(S1%r<=(S2%t.SUB.'0')) lesseq=.TRUE.
       case(m23)
          IF((S1%t.SUB.'0')<=S2%r) lesseq=.TRUE.
       end select
    case default

           write(6,*) " trouble in lesseq "
           write(6,*)  "s1%kind ,s2%kind ",s1%kind ,s2%kind
 
    end select


  END FUNCTION lesseq

  FUNCTION dlesseqsc( S1, S2 )
    implicit none
    logical(lp) dlesseqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dlesseqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r<=S2) dlesseqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')<=S2) dlesseqsc=.TRUE.
    case default
 
        write(6,*) " trouble in dlesseqsc "
        write(6,*) "s1%kind   ",s1%kind   
 
    end select

  END FUNCTION dlesseqsc

  FUNCTION dsclesseq( S2,S1  )
    implicit none
    logical(lp) dsclesseq
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dsclesseq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2<=S1%r) dsclesseq=.TRUE.
    case(m2)
       IF(S2<=(S1%t.SUB.'0')) dsclesseq=.TRUE.
    case default
 
         write(6,*) " trouble in dsclesseq "
         write(6,*) "s1%kind   ",s1%kind
  
    end select

  END FUNCTION dsclesseq

  FUNCTION lesseqsc( S1, S2 )
    implicit none
    logical(lp) lesseqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    lesseqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r<=S2) lesseqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')<=S2) lesseqsc=.TRUE.
    case default
  
      write(6,*)" trouble in lesseqsc "
      write(6,*) "s1%kind   ",s1%kind 
 
    end select

  END FUNCTION lesseqsc

  FUNCTION sclesseq( S2,S1  )
    implicit none
    logical(lp) sclesseq
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    sclesseq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2<=S1%r) sclesseq=.TRUE.
    case(m2)
       IF(S2<=(S1%t.SUB.'0')) sclesseq=.TRUE.
    case default
 
       write(6,*) " trouble in sclesseq "
       write(6,*) "s1%kind   ",s1%kind
 
    end select

  END FUNCTION sclesseq

  FUNCTION ilesseqsc( S1, S2 )
    implicit none
    logical(lp) ilesseqsc
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    ilesseqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r<=S2) ilesseqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')<=S2) ilesseqsc=.TRUE.
    case default

      write(6,*) " trouble in ilesseqsc "
      write(6,*) "s1%kind   ",s1%kind
 
       ! call !write_e(0)
    end select

  END FUNCTION ilesseqsc

  FUNCTION isclesseq( S2,S1  )
    implicit none
    logical(lp) isclesseq
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    isclesseq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2<=S1%r) isclesseq=.TRUE.
    case(m2)
       IF(S2<=(S1%t.SUB.'0')) isclesseq=.TRUE.
    case default
 
       write(6,*)" trouble in isclesseq "
       write(6,*)"s1%kind   ",s1%kind 
 
    end select

  END FUNCTION isclesseq

  FUNCTION eq( S1, S2 )
    implicit none
    logical(lp) eq
    TYPE (real_8), INTENT (IN) :: S1,S2

    eq=.FALSE.
    select case(s1%kind+ms*s2%kind)
    case(m11)
       IF(S1%r==S2%r) eq=.TRUE.
    case(m12,m21,m22)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          IF((S1%t.SUB.'0')==S2%r) eq=.TRUE.
       case(m12)
          IF(S1%r==(S2%t.SUB.'0')) eq=.TRUE.
       case(m22)
          IF((S1%t.SUB.'0')==(S2%t.SUB.'0')) eq=.TRUE.
       end select
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31,m13,m33)
          IF(S1%r==S2%r) eq=.TRUE.
       case(m32)
          IF(S1%r==(S2%t.SUB.'0')) eq=.TRUE.
       case(m23)
          IF((S1%t.SUB.'0')==S2%r) eq=.TRUE.
       end select
    case default

       write(6,*)" trouble in eq "
       write(6,*) "s1%kind ,s2%kind ",s1%kind ,s2%kind 
 
    end select


  END FUNCTION eq

  FUNCTION deqsc( S1, S2 )
    implicit none
    logical(lp) deqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    deqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r==S2) deqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')==S2) deqsc=.TRUE.
    case default
 
       write(6,*)  " trouble in deqsc "
       write(6,*) "s1%kind   ",s1%kind
 
    end select

  END FUNCTION deqsc

  FUNCTION dsceq( S2,S1  )
    implicit none
    logical(lp) dsceq
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dsceq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2==S1%r) dsceq=.TRUE.
    case(m2)
       IF(S2==(S1%t.SUB.'0')) dsceq=.TRUE.
    case default
 
         write(6,*) " trouble in dsceq "
         write(6,*) "s1%kind   "
 
    end select

  END FUNCTION dsceq

  FUNCTION eqsc( S1, S2 )
    implicit none
    logical(lp) eqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    eqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r==S2) eqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')==S2) eqsc=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in eqsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION eqsc

  FUNCTION sceq( S2,S1  )
    implicit none
    logical(lp) sceq
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    sceq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2==S1%r) sceq=.TRUE.
    case(m2)
       IF(S2==(S1%t.SUB.'0')) sceq=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in sceq "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION sceq


  FUNCTION ieqsc( S1, S2 )
    implicit none
    logical(lp) ieqsc
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    ieqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r==S2) ieqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')==S2) ieqsc=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in ieqsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION ieqsc

  FUNCTION isceq( S2,S1  )
    implicit none
    logical(lp) isceq
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    isceq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2==S1%r) isceq=.TRUE.
    case(m2)
       IF(S2==(S1%t.SUB.'0')) isceq=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in isceq "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION isceq

  FUNCTION neq( S1, S2 )
    implicit none
    logical(lp) neq
    TYPE (real_8), INTENT (IN) :: S1,S2

    neq=.FALSE.
    select case(s1%kind+ms*s2%kind)
    case(m11)
       IF(S1%r/=S2%r) neq=.TRUE.
    case(m12,m21,m22)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          IF((S1%t.SUB.'0')/=S2%r) neq=.TRUE.
       case(m12)
          IF(S1%r/=(S2%t.SUB.'0')) neq=.TRUE.
       case(m22)
          IF((S1%t.SUB.'0')/=(S2%t.SUB.'0')) neq=.TRUE.
       end select
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31,m13,m33)
          IF(S1%r/=S2%r) neq=.TRUE.
       case(m32)
          IF(S1%r/=(S2%t.SUB.'0')) neq=.TRUE.
       case(m23)
          IF((S1%t.SUB.'0')/=S2%r) neq=.TRUE.
       end select
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in neq "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select


  END FUNCTION neq

  FUNCTION dneqsc( S1, S2 )
    implicit none
    logical(lp) dneqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dneqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r/=S2) dneqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')/=S2) dneqsc=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dneqsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION dneqsc

  FUNCTION dscneq( S2,S1  )
    implicit none
    logical(lp) dscneq
    TYPE (real_8), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: S2

    dscneq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2/=S1%r) dscneq=.TRUE.
    case(m2)
       IF(S2/=(S1%t.SUB.'0')) dscneq=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dscneq "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION dscneq

  FUNCTION neqsc( S1, S2 )
    implicit none
    logical(lp) neqsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    neqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r/=S2) neqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')/=S2) neqsc=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in neqsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION neqsc

  FUNCTION scneq( S2,S1  )
    implicit none
    logical(lp) scneq
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2

    if(real_warning) call real_stop
    scneq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2/=S1%r) scneq=.TRUE.
    case(m2)
       IF(S2/=(S1%t.SUB.'0')) scneq=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in scneq "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION scneq

  FUNCTION ineqsc( S1, S2 )
    implicit none
    logical(lp) ineqsc
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    ineqsc=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S1%r/=S2) ineqsc=.TRUE.
    case(m2)
       IF((S1%t.SUB.'0')/=S2) ineqsc=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in ineqsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION ineqsc

  FUNCTION iscneq( S2,S1  )
    implicit none
    logical(lp) iscneq
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2

    iscneq=.FALSE.
    select case(s1%kind)
    case(m1,m3)
       IF(S2/=S1%r) iscneq=.TRUE.
    case(m2)
       IF(S2/=(S1%t.SUB.'0')) iscneq=.TRUE.
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in iscneq "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION iscneq

  FUNCTION dexpt( S1 )
    implicit none
    TYPE (real_8) dexpt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dexpt%r=exp(s1%r)
       dexpt%kind=1
    case(m2)
       localmaster=master
       call ass(dexpt)
       dexpt%t= exp(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dexpt)

          call varfk1(S1)
          dexpt%t= exp(varf1)
          master=localmaster
       else
          dexpt%r= exp(S1%r)
          dexpt%kind=1
       endif

    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dexpt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dexpt


  FUNCTION abst( S1 )
    implicit none
    real(dp) abst
    TYPE (real_8), INTENT (IN) :: S1


    select case(s1%kind)
    case(m1,m3)
       abst=abs(s1%r)
    case(m2)
       abst=abs(s1%t.sub.'0')
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in abst "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION abst

  FUNCTION Pabs( S1 )
    implicit none
    TYPE (real_8) Pabs
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       Pabs%R=abs(s1%r)
       Pabs%kind=1
    case(m2)
       localmaster=master
       call ass(Pabs)
       IF((s1%t.sub.'0')<0) THEN
          Pabs%T=-s1%t
       ELSE
          Pabs%T=s1%t
       ENDIF
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(Pabs)
          call varfk1(S1)
          IF((s1%t.sub.'0')<0) THEN
             Pabs%T=-varf1
          ELSE
             Pabs%T=varf1
          ENDIF
          master=localmaster
       else
          PABS%r= ABS(S1%r)
          PABS%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in Pabs "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION Pabs



  FUNCTION full_abst( S1 )
    implicit none
    real(dp) full_abst
    TYPE (real_8), INTENT (IN) :: S1


    select case(s1%kind)
    case(m1,m3)
       full_abst=abs(s1%r)
    case(m2)
       full_abst=full_abs(s1%t)    ! 2002.10.17
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in full_abst "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION full_abst

  FUNCTION dtant( S1 )
    implicit none
    TYPE (real_8) dtant
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dtant%r=tan(s1%r)
       dtant%kind=1
    case(m2)
       localmaster=master
       call ass(dtant)
       dtant%t= tan(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dtant)
          call varfk1(S1)
          dtant%t= tan(varf1)
          master=localmaster
       else
          dtant%r= tan(S1%r)
          dtant%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dtant "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dtant

  FUNCTION dtandt( S1 )
    implicit none
    TYPE (real_8) dtandt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dtandt%r=tan(s1%r*DEG_TO_RAD_)
       dtandt%kind=1
    case(m2)
       localmaster=master
       call ass(dtandt)
       dtandt%t= tan(s1%t*DEG_TO_RAD_)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dtandt)
          call varfk1(S1)
          dtandt%t= tan(varf1*DEG_TO_RAD_)
          master=localmaster
       else
          dtandt%r= tan(S1%r*DEG_TO_RAD_)
          dtandt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dtandt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dtandt

  FUNCTION absoftdtandr( S1 )
    implicit none
    real(dp) absoftdtandr
    real(dp), INTENT (IN) :: S1
    absoftdtandr=tan(s1*DEG_TO_RAD_)
  END FUNCTION absoftdtandr

  FUNCTION dcost( S1 )
    implicit none
    TYPE (real_8) dcost
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dcost%r=cos(s1%r)
       dcost%kind=1
    case(m2)
       localmaster=master
       call ass(dcost)
       dcost%t= cos(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dcost)
          call varfk1(S1)
          dcost%t= cos(varf1)
          master=localmaster
       else
          dcost%r= cos(S1%r)
          dcost%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dcost "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dcost

  FUNCTION dcosdt( S1 )
    implicit none
    TYPE (real_8) dcosdt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dcosdt%r=cos(s1%r*DEG_TO_RAD_)
       dcosdt%kind=1
    case(m2)
       localmaster=master
       call ass(dcosdt)
       dcosdt%t= cos(s1%t*DEG_TO_RAD_)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dcosdt)
          call varfk1(S1)
          dcosdt%t= cos(varf1*DEG_TO_RAD_)
          master=localmaster
       else
          dcosdt%r= cos(S1%r)
          dcosdt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dcosdt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dcosdt

  FUNCTION absoftdcosdr( S1 )
    implicit none
    real(dp) absoftdcosdr
    real(dp) , INTENT (IN) :: S1
    absoftdcosdr=cos(s1*DEG_TO_RAD_ )
  END FUNCTION absoftdcosdr


  FUNCTION dsint( S1 )
    implicit none
    TYPE (real_8) dsint
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dsint%r=sin(s1%r)
       dsint%kind=1
    case(m2)
       localmaster=master
       call ass(dsint)
       dsint%t= sin(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dsint)
          call varfk1(S1)
          dsint%t= sin(varf1)
          master=localmaster
       else
          dsint%r= sin(S1%r)
          dsint%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dsint "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dsint

  FUNCTION dsindt( S1 )
    implicit none
    TYPE (real_8) dsindt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dsindt%r=sin(s1%r*DEG_TO_RAD_ )
       dsindt%kind=1
    case(m2)
       localmaster=master
       call ass(dsindt)
       dsindt%t= sin(s1%t*DEG_TO_RAD_)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dsindt)
          call varfk1(S1)
          dsindt%t= sin(varf1*DEG_TO_RAD_)
          master=localmaster
       else
          dsindt%r= sin(S1%r*DEG_TO_RAD_)
          dsindt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dsindt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dsindt

  FUNCTION absoftdsindr( S1 )
    implicit none
    real(dp) absoftdsindr
    real(dp) , INTENT (IN) :: S1
    absoftdsindr=sin(s1*DEG_TO_RAD_ )
  END FUNCTION absoftdsindr


  FUNCTION dlogt( S1 )
    implicit none
    TYPE (real_8) dlogt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dlogt%r=loge(s1%r)
       dlogt%kind=1
    case(m2)
       localmaster=master
       call ass(dlogt)
       dlogt%t= log(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dlogt)
          call varfk1(S1)
          dlogt%t= log(varf1)
          master=localmaster
       else
          dlogt%r= log(S1%r)
          dlogt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dlogt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dlogt

  FUNCTION dsqrtt( S1 )
    implicit none
    TYPE (real_8) dsqrtt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       !       dsqrtt%r=SQRT(s1%r)
       dsqrtt%r=root(s1%r)
       dsqrtt%kind=1
    case(m2)
       localmaster=master
       call ass(dsqrtt)
       dsqrtt%t= SQRT(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dsqrtt)
          call varfk1(S1)
          dsqrtt%t= SQRT(varf1)
          master=localmaster
       else
          !          dsqrtt%r= SQRT(S1%r)
          dsqrtt%r= root(S1%r)
          dsqrtt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dsqrtt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dsqrtt



  FUNCTION POW( S1, S2 )
    implicit none
    TYPE (real_8) POW
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       POW%r=s1%r**s2
       POW%kind=1
    case(m2)
       localmaster=master
       call ass(POW)
       POW%t= s1%t**s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(POW)
          call varfk1(S1)
          POW%t= varf1**s2
          master=localmaster
       else
          POW%r= S1%r**s2
          POW%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in POW "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION POW

  FUNCTION POWR( S1, S2 )
    implicit none
    TYPE (real_8) POWR
    TYPE (real_8), INTENT (IN) :: S1
    REAL(SP) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       POWR%r=s1%r**REAL(s2,kind=DP)
       POWR%kind=1
    case(m2)
       localmaster=master
       call ass(POWR)
       POWR%t= s1%t**REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(POWR)
          call varfk1(S1)
          POWR%t= varf1**REAL(s2,kind=DP)
          master=localmaster
       else
          POWR%r= S1%r**REAL(s2,kind=DP)
          POWR%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in POWR "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION POWR

  FUNCTION POWR8( S1, S2 )
    implicit none
    TYPE (real_8) POWR8
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       POWR8%r=s1%r**s2
       POWR8%kind=1
    case(m2)
       localmaster=master
       call ass(POWR8)
       POWR8%t= s1%t**s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(POWR8)
          call varfk1(S1)
          POWR8%t= varf1**s2
          master=localmaster
       else
          POWR8%r= S1%r**s2
          POWR8%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in POWR8 "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION POWR8


  FUNCTION unaryADD( S1 )
    implicit none
    TYPE (real_8) unaryADD
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       unaryADD%r=s1%r
       unaryADD%kind=1
    case(m2)
       localmaster=master
       call ass(unaryADD)
       unaryADD%t= s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(unaryADD)
          call varfk1(S1)
          unaryADD%t= varf1
          master=localmaster
       else
          unaryADD%r= S1%r
          unaryADD%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in POWR8 "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION unaryADD

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (real_8) unarySUB
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       unarySUB%r=-s1%r
       unarySUB%kind=1
    case(m2)
       localmaster=master
       call ass(unarySUB)
       unarySUB%t= -s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(unarySUB)
          call varfk1(S1)
          unarySUB%t= -varf1
          master=localmaster
       else
          unarySUB%r= -S1%r
          unarySUB%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in unarySUB "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION unarySUB

  FUNCTION nbip( n,S1, S2 )
    implicit none
    TYPE (real_8) nbip
    TYPE (real_8), INTENT (IN) :: S1, S2
    integer, intent (in) :: n
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       nbip%r=nbi(n,s1%r,s2%r)
       nbip%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(nbip)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          nbip%t= nbi(n,s1%t,s2%r) 
       case(m12)
          nbip%t=  nbi(n,s1%r,s2%t) 
       case(m22)
          nbip%t= nbi(n,s1%t,s2%t) 
       end select
 
       master=localmaster

    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(nbip)
             call varfk1(S1)
             nbip%t= nbi(n,varf1,s2%r)
             master=localmaster
          else
             nbip%r= nbi(n,s1%r,s2%r)
             nbip%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(nbip)
             call varfk1(S2)
             nbip%t= nbi(n,s1%r,varf1)  
             master=localmaster
          else
             nbip%r= nbi(n,s1%r,s2%r)
             nbip%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(nbip)
          if(knob) then
             call varfk1(S1)
             nbip%t=  nbi(n,varf1,s2%t)  
          else
             nbip%t= nbi(n,s1%r,s2%t)    
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(nbip)
          if(knob) then
             call varfk1(S2)
             nbip%t=  nbi(n,s1%t,varf1)  
          else
             nbip%t=  nbi(n,s1%t,s2%r)    
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(nbip)
             call varfk1(S1)
             call varfk2(S2)
             nbip%t=  nbi(n,varf1,varf2)     
             master=localmaster
          else
             nbip%r=nbi(n,s1%r,s2%r)   
             nbip%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in nbip "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION nbip

  FUNCTION add( S1, S2 )
    implicit none
    TYPE (real_8) add
    TYPE (real_8), INTENT (IN) :: S1, S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       add%r=s1%r+s2%r
       add%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(add)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          add%t= s1%t+s2%r
       case(m12)
          add%t= s1%r+s2%t
       case(m22)
          add%t= s1%t+s2%t
       end select
       master=localmaster

    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(add)
             call varfk1(S1)
             add%t= varf1+s2%r
             master=localmaster
          else
             add%r= s1%r+s2%r
             add%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(add)
             call varfk1(S2)
             add%t= s1%r+varf1
             master=localmaster
          else
             add%r= s1%r+s2%r
             add%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(add)
          if(knob) then
             call varfk1(S1)
             add%t= varf1+s2%t
          else
             add%t= s1%r+s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(add)
          if(knob) then
             call varfk1(S2)
             add%t= s1%t+varf1
          else
             add%t= s1%t+s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(add)
             call varfk1(S1)
             call varfk2(S2)
             add%t= varf1+varf2
             master=localmaster
          else
             add%r= s1%r+s2%r
             add%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in add "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION add

  FUNCTION addq( S1, S2 )
    implicit none
    TYPE (quaternion_8) addq
    TYPE (quaternion_8), INTENT (IN) :: S1, S2
    type(real_8) temp
    integer i,localmaster


       call alloc(temp)
       do i=0,3
        temp=s1%x(i)+s2%x(i)
       if(temp%kind==2) then
             localmaster=master
        call ass(addq%x(i))
         addq%x(i)=temp
          master=localmaster
        else
         addq%x(i)%r=temp%r
         addq%x(i)%kind=1
        endif
       enddo
       call kill(temp)

  END FUNCTION addq

  FUNCTION absq2( S1 )
    implicit none
    real(dp) absq2
    TYPE (quaternion_8), INTENT (IN) :: S1
    integer i

    IF(.NOT.C_%STABLE_DA) then
     absq2=0
     RETURN
    endif
           absq2=0
       do i=0,3
         absq2 = s1%x(i)**2+absq2
       enddo
  END FUNCTION absq2

  FUNCTION absq( S1 )
    implicit none
    real(dp) absq
    TYPE (quaternion_8), INTENT (IN) :: S1
    integer i

    IF(.NOT.C_%STABLE_DA) then
     absq=0
     RETURN
    endif
           absq=sqrt(abs_square(s1))
            
  END FUNCTION absq

  FUNCTION subq( S1, S2 )
    implicit none
    TYPE (quaternion_8) subq
    TYPE (quaternion_8), INTENT (IN) :: S1, S2
    type(real_8) temp
    integer i,localmaster
 
 
       call alloc(temp)
       do i=0,3
        temp=s1%x(i)-s2%x(i)
       if(temp%kind==2) then
             localmaster=master
        call ass(subq%x(i))
         subq%x(i)%t=temp%t
          master=localmaster
        else

         subq%x(i)%r=temp%r
         subq%x(i)%kind=1
        endif
       enddo
       call kill(temp)

  END FUNCTION subq

  FUNCTION invq( S1 )
    implicit none
    TYPE (quaternion_8) invq
    TYPE (quaternion_8), INTENT (IN) :: S1
    type(real_8) norm
    type(real_8) temp(0:3)
    integer i,localmaster

    IF(.NOT.C_%STABLE_DA) then
     invq%x(1)=0
     RETURN
    endif
 
       call alloc(temp)
      call alloc(norm)

              norm=s1%x(0)**2+s1%x(1)**2+s1%x(2)**2+s1%x(3)**2

            temp(1)=s1%x(1)
              do i=1,3
                temp(i)=-s1%x(i)
              enddo
              do i=0,3
                temp(i)=temp(i)/norm
              enddo

    do i=0,3
       if(temp(i)%kind==2) then
             localmaster=master
        call ass(invq%x(i))
         invq%x(i)=temp(i)
          master=localmaster
        else
         invq%x(i)%r=temp(i)%r
         invq%x(i)%kind=1
        endif
    enddo

          call kill(norm)
          call kill(temp)

 
  END FUNCTION invq

  FUNCTION POWq( S1, R2 )
    implicit none
    TYPE (quaternion_8) POWq,qtemp
    TYPE (quaternion_8), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       POWq=0.0_dp
      RETURN
    endif

 
     call alloc(qtemp)    
     qtemp=1.0_dp

    R22=IABS(R2)
    DO I=1,R22
       qtemp=qtemp*s1
    ENDDO
    IF(R2.LT.0) THEN
       qtemp=invq(qtemp)
    ENDIF

     do i=0,3
       if(qtemp%x(i)%kind==2) then
             localmaster=master
        call ass(powq%x(i))
         powq%x(i)=qtemp%x(i)
        master=localmaster
        else
         powq%x(i)%r=qtemp%x(i)%r
         powq%x(i)%kind=1
        endif

     enddo
      call kill(qtemp)
        master=localmaster
  END FUNCTION POWq

  FUNCTION mulq( S1, S2 )
    implicit none
    TYPE (quaternion_8) mulq
    TYPE (quaternion_8), INTENT (IN) :: S1, S2
    type(real_8) temp(0:3)
    integer i,localmaster

 !       call ass_quaternion_8(mulq)
        call alloc(temp)
 
         temp(0)=s1%x(0)*s2%x(0)-s1%x(1)*s2%x(1)-s1%x(2)*s2%x(2)-s1%x(3)*s2%x(3)

         temp(1)=s1%x(2)*s2%x(3)-s1%x(3)*s2%x(2)
         temp(2)=s1%x(3)*s2%x(1)-s1%x(1)*s2%x(3)
         temp(3)=s1%x(1)*s2%x(2)-s1%x(2)*s2%x(1)

        do i=1,3
         temp(i)= temp(i) + s1%x(0)*s2%x(i)+ s1%x(i)*s2%x(0)
        enddo

        do i=0,3
            if(temp(i)%kind==2) then
             localmaster=master
               call ass(mulq%x(i))

                  mulq%x(i)=temp(i)
                 master=localmaster
             else
              mulq%x(i)%r=temp(i)%r
              mulq%x(i)%kind=1
             endif      
       enddo
       call kill(temp)


  END FUNCTION mulq



  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (real_8) subs
    TYPE (real_8), INTENT (IN) :: S1, S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       subs%r=s1%r-s2%r
       subs%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(subs)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          subs%t= s1%t-s2%r
       case(m12)
          subs%t= s1%r-s2%t
       case(m22)
          subs%t= s1%t-s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)

          if(knob) then
             localmaster=master
             call ass(subs)
             call varfk1(S1)
             subs%t= varf1-s2%r
             master=localmaster
          else
             subs%r= s1%r-s2%r
             subs%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(subs)
             call varfk1(S2)
             subs%t= s1%r-varf1
             master=localmaster
          else
             subs%r= s1%r-s2%r
             subs%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(subs)
          if(knob) then
             call varfk1(S1)
             subs%t= varf1-s2%t
          else
             subs%t= s1%r-s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(subs)
          if(knob) then
             call varfk1(S2)
             subs%t= s1%t-varf1
          else
             subs%t= s1%t-s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(subs)
             call varfk1(S1)
             call varfk2(S2)
             subs%t= varf1-varf2
             master=localmaster
          else
             subs%r= s1%r-s2%r
             subs%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in subs "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION subs

  FUNCTION daddsc( S1, S2 )
    implicit none
    TYPE (real_8) daddsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       daddsc%r=s1%r+s2
       daddsc%kind=1
    case(m2)
       localmaster=master
       call ass(daddsc)
       daddsc%t= s1%t+s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(daddsc)
          call varfk1(S1)
          daddsc%t= varf1+s2
          master=localmaster
       else
          daddsc%r= s1%r+s2
          daddsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in daddsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION daddsc

  FUNCTION dscadd( S2, S1 )
    implicit none
    TYPE (real_8) dscadd
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dscadd%r=s1%r+s2
       dscadd%kind=1
    case(m2)
       localmaster=master
       call ass(dscadd)
       dscadd%t= s1%t+s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dscadd)
          call varfk1(S1)
          dscadd%t= varf1+s2
          master=localmaster
       else
          dscadd%r= S1%r+s2
          dscadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dscadd "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dscadd

  FUNCTION dsubsc( S1, S2 )
    implicit none
    TYPE (real_8) dsubsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dsubsc%r=s1%r-s2
       dsubsc%kind=1
    case(m2)
       localmaster=master
       call ass(dsubsc)
       dsubsc%t= s1%t-s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dsubsc)
          call varfk1(S1)
          dsubsc%t= varf1-s2
          master=localmaster
       else
          dsubsc%r=s1%r-s2
          dsubsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dsubsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dsubsc

  FUNCTION dscsub( S2, S1 )
    implicit none
    TYPE (real_8) dscsub
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dscsub%r=s2-s1%r
       dscsub%kind=1
    case(m2)
       localmaster=master
       call ass(dscsub)
       dscsub%t=s2-s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dscsub)
          call varfk1(S1)
          dscsub%t=s2-varf1
          master=localmaster
       else
          dscsub%r=s2-s1%r
          dscsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dscsub "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dscsub

  FUNCTION addsc( S1, S2 )
    implicit none
    TYPE (real_8) addsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       addsc%r=s1%r+REAL(s2,kind=DP)
       addsc%kind=1
    case(m2)
       localmaster=master
       call ass(addsc)
       addsc%t= s1%t+REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(addsc)
          call varfk1(S1)
          addsc%t= varf1+REAL(s2,kind=DP)
          master=localmaster
       else
          addsc%r=s1%r+REAL(s2,kind=DP)
          addsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in addsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION addsc

  FUNCTION scadd( S2, S1 )
    implicit none
    TYPE (real_8) scadd
    TYPE (real_8), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       scadd%r=s1%r+REAL(s2,kind=DP)
       scadd%kind=1
    case(m2)
       localmaster=master
       call ass(scadd)
       scadd%t= s1%t+REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(scadd)
          call varfk1(S1)
          scadd%t= varf1+REAL(s2,kind=DP)
          master=localmaster
       else
          scadd%r=s1%r+REAL(s2,kind=DP)
          scadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in addsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION scadd

  FUNCTION subsc( S1, S2 )
    implicit none
    TYPE (real_8) subsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       subsc%r=s1%r-REAL(s2,kind=DP)
       subsc%kind=1
    case(m2)
       localmaster=master
       call ass(subsc)
       subsc%t= s1%t-REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(subsc)
          call varfk1(S1)
          subsc%t= varf1-REAL(s2,kind=DP)
          master=localmaster
       else
          subsc%r=s1%r-REAL(s2,kind=DP)
          subsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in subsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION subsc

  FUNCTION scsub( S2, S1 )
    implicit none
    TYPE (real_8) scsub
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       scsub%r=REAL(s2,kind=DP)-s1%r
       scsub%kind=1
    case(m2)
       localmaster=master
       call ass(scsub)
       scsub%t=REAL(s2,kind=DP)-s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(scsub)
          call varfk1(S1)
          scsub%t=REAL(s2,kind=DP)-varf1
          master=localmaster
       else
          scsub%r=REAL(s2,kind=DP)-s1%r
          scsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in scsub "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION scsub

  FUNCTION iaddsc( S1, S2 )
    implicit none
    TYPE (real_8) iaddsc
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       iaddsc%r=s1%r+REAL(s2,kind=DP)
       iaddsc%kind=1
    case(m2)
       localmaster=master
       call ass(iaddsc)
       iaddsc%t= s1%t+REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(iaddsc)
          call varfk1(S1)
          iaddsc%t= varf1+REAL(s2,kind=DP)
          master=localmaster
       else
          iaddsc%r=s1%r+REAL(s2,kind=DP)
          iaddsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in iaddsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION iaddsc

  FUNCTION iscadd( S2, S1 )
    implicit none
    TYPE (real_8) iscadd
    TYPE (real_8), INTENT (IN) :: S1
    integer, INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       iscadd%r=s1%r+REAL(s2,kind=DP)
       iscadd%kind=1
    case(m2)
       localmaster=master
       call ass(iscadd)
       iscadd%t= s1%t+REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(iscadd)
          call varfk1(S1)
          iscadd%t= varf1+REAL(s2,kind=DP)
          master=localmaster
       else
          iscadd%r=s1%r+REAL(s2,kind=DP)
          iscadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in iscadd "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION iscadd

  FUNCTION isubsc( S1, S2 )
    implicit none
    TYPE (real_8) isubsc
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       isubsc%r=s1%r-REAL(s2,kind=DP)
       isubsc%kind=1
    case(m2)
       localmaster=master
       call ass(isubsc)
       isubsc%t= s1%t-REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(isubsc)
          call varfk1(S1)
          isubsc%t= varf1-REAL(s2,kind=DP)
          master=localmaster
       else
          isubsc%r=s1%r-REAL(s2,kind=DP)
          isubsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in isubsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION isubsc

  FUNCTION iscsub( S2, S1 )
    implicit none
    TYPE (real_8) iscsub
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       iscsub%r=REAL(s2,kind=DP)-s1%r
       iscsub%kind=1
    case(m2)
       localmaster=master
       call ass(iscsub)
       iscsub%t=REAL(s2,kind=DP)-s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(iscsub)
          call varfk1(S1)
          iscsub%t=REAL(s2,kind=DP)-varf1
          master=localmaster
       else
          iscsub%r=REAL(s2,kind=DP)-s1%r
          iscsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in iscsub "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION iscsub

  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (real_8) mul
    TYPE (real_8), INTENT (IN) :: S1, S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       mul%r=s1%r*s2%r
       mul%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(mul)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          mul%t= s1%t*s2%r
       case(m12)
          mul%t= s1%r*s2%t
       case(m22)
          mul%t= s1%t*s2%t
       end select
       master=localmaster

    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)

          if(knob) then
             localmaster=master
             call ass(mul)
             call varfk1(S1)
             mul%t= varf1*s2%r
             master=localmaster
          else
             mul%r= s1%r*s2%r
             mul%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(mul)
             call varfk1(S2)
             mul%t= s1%r*varf1
             master=localmaster
          else
             mul%r= s1%r*s2%r
             mul%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(mul)
          if(knob) then
             call varfk1(S1)
             mul%t= varf1*s2%t
          else
             mul%t= s1%r*s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(mul)
          if(knob) then
             call varfk1(S2)
             mul%t= s1%t*varf1
          else
             mul%t= s1%t*s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(mul)
             call varfk1(S1)
             call varfk2(S2)
             mul%t= varf1*varf2
             master=localmaster
          else
             mul%r= s1%r*s2%r
             mul%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in mul "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION mul

  FUNCTION div( S1, S2 )
    implicit none
    TYPE (real_8) div
    TYPE (real_8), INTENT (IN) :: S1, S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       div%r=s1%r/s2%r
       div%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(div)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          div%t= s1%t/s2%r
       case(m12)
          div%t= s1%r/s2%t
       case(m22)
          div%t= s1%t/s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)

          if(knob) then
             localmaster=master
             call ass(div)
             call varfk1(S1)
             div%t= varf1/s2%r
             master=localmaster
          else
             div%r= s1%r/s2%r
             div%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(div)
             call varfk1(S2)
             div%t= s1%r/varf1
             master=localmaster
          else
             div%r= s1%r/s2%r
             div%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(div)
          if(knob) then
             call varfk1(S1)
             div%t= varf1/s2%t
          else
             div%t= s1%r/s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(div)
          if(knob) then
             call varfk1(S2)
             div%t= s1%t/varf1
          else
             div%t= s1%t/s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(div)
             call varfk1(S1)
             call varfk2(S2)
             div%t= varf1/varf2
             master=localmaster
          else
             div%r= s1%r/s2%r
             div%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in div "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION div

  FUNCTION dmulmapconcat( S1, S2 )
    implicit none
    TYPE (real_8) dmulmapconcat
    TYPE (real_8), INTENT (IN) :: S1
    type(damap) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dmulmapconcat%r=s1%r
       dmulmapconcat%kind=1
    case(m2)
       localmaster=master
       call ass(dmulmapconcat)
       dmulmapconcat%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dmulmapconcat)
          call varfk1(S1)
          dmulmapconcat%t= varf1*s2
          master=localmaster
       else
          dmulmapconcat%r=s1%r
          dmulmapconcat%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dmulmapconcat "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dmulmapconcat


  FUNCTION dmulsc( S1, S2 )
    implicit none
    TYPE (real_8) dmulsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dmulsc%r=s1%r*s2
       dmulsc%kind=1
    case(m2)
       localmaster=master
       call ass(dmulsc)
       dmulsc%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dmulsc)
          call varfk1(S1)
          dmulsc%t= varf1*s2
          master=localmaster
       else
          dmulsc%r=s1%r*s2
          dmulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dmulsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dmulsc

  FUNCTION dscmul( S2, S1 )
    implicit none
    TYPE (real_8) dscmul
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dscmul%r=s1%r*s2
       dscmul%kind=1
    case(m2)
       localmaster=master
       call ass(dscmul)
       dscmul%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dscmul)
          call varfk1(S1)
          dscmul%t= varf1*s2
          master=localmaster
       else
          dscmul%r=s1%r*s2
          dscmul%kind=1
       endif
    end select
  END FUNCTION dscmul

  FUNCTION ddivsc( S1, S2 )
    implicit none
    TYPE (real_8) ddivsc
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       ddivsc%r=s1%r/s2
       ddivsc%kind=1
    case(m2)
       localmaster=master
       call ass(ddivsc)
       ddivsc%t= s1%t/s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(ddivsc)
          call varfk1(S1)
          ddivsc%t= varf1/s2
          master=localmaster
       else
          ddivsc%r=s1%r/s2
          ddivsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in ddivsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION ddivsc

  FUNCTION dscdiv( S2, S1 )
    implicit none
    TYPE (real_8) dscdiv
    TYPE (real_8), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       dscdiv%r=s2/s1%r
       dscdiv%kind=1
    case(m2)
       localmaster=master
       call ass(dscdiv)
       dscdiv%t= s2/s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dscdiv)
          call varfk1(S1)
          dscdiv%t= s2/varf1
          master=localmaster
       else
          dscdiv%r=s2/s1%r
          dscdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in ddivsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dscdiv

  FUNCTION mulsc( S1, S2 )
    implicit none
    TYPE (real_8) mulsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       mulsc%r=s1%r*REAL(s2,kind=DP)
       mulsc%kind=1
    case(m2)
       localmaster=master
       call ass(mulsc)
       mulsc%t= s1%t*REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(mulsc)
          call varfk1(S1)
          mulsc%t= varf1*REAL(s2,kind=DP)
          master=localmaster
       else
          mulsc%r=s1%r*REAL(s2,kind=DP)
          mulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in mulsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION mulsc

  FUNCTION scmul( S2, S1 )
    implicit none
    TYPE (real_8) scmul
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       scmul%r=s1%r*REAL(s2,kind=DP)
       scmul%kind=1
    case(m2)
       localmaster=master
       call ass(scmul)
       scmul%t= s1%t*REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(scmul)
          call varfk1(S1)
          scmul%t= varf1*REAL(s2,kind=DP)
          master=localmaster
       else
          scmul%r=s1%r*REAL(s2,kind=DP)
          scmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in scmul "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION scmul

  FUNCTION divsc( S1, S2 )
    implicit none
    TYPE (real_8) divsc
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       divsc%r=s1%r/REAL(s2,kind=DP)
       divsc%kind=1
    case(m2)
       localmaster=master
       call ass(divsc)
       divsc%t= s1%t/REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(divsc)
          call varfk1(S1)
          divsc%t= varf1/REAL(s2,kind=DP)
          master=localmaster
       else
          divsc%r=s1%r/REAL(s2,kind=DP)
          divsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in divsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION divsc

  FUNCTION scdiv( S2, S1 )
    implicit none
    TYPE (real_8) scdiv
    TYPE (real_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
    integer localmaster

    if(real_warning) call real_stop
    select case(s1%kind)
    case(m1)
       scdiv%r=REAL(s2,kind=DP)/s1%r
       scdiv%kind=1
    case(m2)
       localmaster=master
       call ass(scdiv)
       scdiv%t= REAL(s2,kind=DP)/s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(scdiv)
          call varfk1(S1)
          scdiv%t= REAL(s2,kind=DP)/varf1
          master=localmaster
       else
          scdiv%r=REAL(s2,kind=DP)/s1%r
          scdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in scdiv "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION scdiv

  FUNCTION imulsc( S1, S2 )
    implicit none
    TYPE (real_8) imulsc
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       imulsc%r=s1%r*REAL(s2,kind=DP)
       imulsc%kind=1
    case(m2)
       if(s2/=0) then
        localmaster=master
        call ass(imulsc)
        imulsc%t= s1%t*REAL(s2,kind=DP)
        master=localmaster
       else
        imulsc%r=0.0_dp
        imulsc%kind=1
       endif
    case(m3)
       if(knob) then
          localmaster=master
          call ass(imulsc)
          call varfk1(S1)
          imulsc%t= varf1*REAL(s2,kind=DP)
          master=localmaster
       else
          imulsc%r=s1%r*REAL(s2,kind=DP)
          imulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in imulsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION imulsc

  FUNCTION iscmul( S2, S1 )
    implicit none
    TYPE (real_8) iscmul
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       iscmul%r=s1%r*REAL(s2,kind=DP)
       iscmul%kind=1
    case(m2)
       if(s2/=0) then
        localmaster=master
        call ass(iscmul)
        iscmul%t= s1%t*REAL(s2,kind=DP)
        master=localmaster
       else
        iscmul%r=0.0_dp
        iscmul%kind=1
       endif
    case(m3)
       if(knob) then
          localmaster=master
          call ass(iscmul)
          call varfk1(S1)
          iscmul%t= varf1*REAL(s2,kind=DP)
          master=localmaster
       else
          iscmul%r=s1%r*REAL(s2,kind=DP)
          iscmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in iscmul "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION iscmul

  FUNCTION idivsc( S1, S2 )
    implicit none
    TYPE (real_8) idivsc
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       idivsc%r=s1%r/REAL(s2,kind=DP)
       idivsc%kind=1
    case(m2)
       localmaster=master
       call ass(idivsc)
       idivsc%t= s1%t/REAL(s2,kind=DP)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(idivsc)
          call varfk1(S1)
          idivsc%t= varf1/REAL(s2,kind=DP)
          master=localmaster
       else
          idivsc%r=s1%r/REAL(s2,kind=DP)
          idivsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in idivsc "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION idivsc

  FUNCTION iscdiv( S2, S1 )
    implicit none
    TYPE (real_8) iscdiv
    TYPE (real_8), INTENT (IN) :: S1
    integer , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       iscdiv%r=REAL(s2,kind=DP)/s1%r
       iscdiv%kind=1
    case(m2)
       localmaster=master
       call ass(iscdiv)
       iscdiv%t= REAL(s2,kind=DP)/s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(iscdiv)
          call varfk1(S1)
          iscdiv%t= REAL(s2,kind=DP)/varf1
          master=localmaster
       else
          iscdiv%r=REAL(s2,kind=DP)/s1%r
          iscdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in iscdiv "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION iscdiv

  SUBROUTINE  printpolyq(S2,mf,prec)
    implicit none
    integer ipause, mypauses,i,k
    type (quaternion_8),INTENT(INOUT)::S2
    integer,optional :: mf
    real(dp), optional :: prec
    i=6
    if(present(mf)) i=mf
        write(i,*) " quaternion_8 " 
    do k=0,3
      call print(s2%x(k),i,prec)
    enddo

  END SUBROUTINE printpolyq

  SUBROUTINE  printpoly(S2,mf,prec)
    implicit none
    integer ipause, mypauses,i
    type (real_8),INTENT(INOUT)::S2
    integer,optional :: mf
    real(dp), optional :: prec
    i=6
    if(present(mf)) i=mf
    if(s2%kind/=0) then

       select  case (s2%kind)
       case(m1)
          write(i,*)  s2%r
       case(m2)
          call pri(S2%t,i,prec)
       case(m3)

          if(s2%i>0) then
             write(line,*) s2%r,"  +",s2%s,"  (x_",s2%i,")"

         else
          write(line,*) s2%r
          endif
             call context(line,maj=.false.)
             write(i,'(a)') adjustr(line(1:len_trim(line)))
          if(s2%alloc) then
             write(line,'(a41)')  " weird Taylor part should be deallocated "
             ipause=mypauses(0,line)
          endif
       end   select
    else

       line=  "Warning not defined in Printpoly (real polymorph)"
       ipause=mypauses(1,line)

    endif

  END SUBROUTINE printpoly

  SUBROUTINE  print6(S1,mf)
    implicit none
    type (real_8),INTENT(INout)::S1(:)
    integer,optional :: mf
    integer        i
    
 !   if(size(s1)==6) then
 !    do i=1,ndd
 !       call print(s1(i),mf)
 !    enddo
 !   else
     do i=lbound(s1,1),ubound(s1,1)
        call print(s1(i),mf)
     enddo
 !   endif
  END SUBROUTINE print6


  SUBROUTINE  printdouble(S2,mf)
    implicit none
    real(dp),INTENT(INOUT)::S2
    integer,optional :: mf
    integer i
 
    i=6
    if(present(mf)) i=mf
    write(mf,*)  s2

  END SUBROUTINE printdouble

  SUBROUTINE  printsingle(S2,mf)
    implicit none
    real(sp),INTENT(INOUT)::S2
    integer,optional :: mf
    integer i
 
    i=6
    if(present(mf)) i=mf
    write(mf,*)  s2

  END SUBROUTINE printsingle

  SUBROUTINE  resetpoly(S2)
    implicit none
    type (real_8),INTENT(INOUT)::S2

    if(s2%alloc) call killtpsa(s2%t)
    s2%alloc=f
    s2%kind=0
    s2%r=0.0_dp
    !s2%s=one
    !s2%i=0

  END SUBROUTINE resetpoly


  SUBROUTINE  resetpolyn(S2,K)
    implicit none
    type (real_8),INTENT(INOUT),dimension(:)::S2
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S2,DIM=1)
       N=LBOUND(S2,DIM=1)+K-1
    else
       I=LBOUND(S2,DIM=1)
       N=UBOUND(S2,DIM=1)
    endif

    DO   J=I,N
       call resetpoly(s2(j))
    enddo


  END SUBROUTINE resetpolyn


  SUBROUTINE  resetpoly_R(S2,FL)  !   STAYS REAL
    implicit none
    type (real_8),INTENT(INOUT)::S2
    logical(lp),INTENT(IN)::FL

    if(s2%alloc) call killtpsa(s2%t)
    s2%alloc=f
    s2%kind=1
    s2%r=0.0_dp

    IF(.NOT.FL) THEN
       s2%i=0
       s2%s=1.0_dp
    ENDIF

  END SUBROUTINE resetpoly_R

  SUBROUTINE  resetpoly_RN(S2,FL,K)
    implicit none
    type (real_8),INTENT(INOUT),dimension(:)::S2
    logical(lp),INTENT(IN)::FL

    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S2,DIM=1)
       N=LBOUND(S2,DIM=1)+K-1
    else
       I=LBOUND(S2,DIM=1)
       N=UBOUND(S2,DIM=1)
    endif

    DO   J=I,N
       call resetpoly_R(s2(j),FL)
    enddo


  END SUBROUTINE resetpoly_RN

  SUBROUTINE  resetpoly_R31(S2)  !   STAYS REAL FOR PTC
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(INOUT)::S2

    if(s2%kind==3) then
       if(s2%alloc) then
          line=  "Allocated in resetpoly_R31"
          ipause=mypauses(2,line)
       endif
       s2%kind=1
       s2%i=0
       s2%s=1.0_dp
    ENDIF
  END SUBROUTINE resetpoly_R31

  SUBROUTINE  resetpoly_R31N(S2,K)
    implicit none
    type (real_8),INTENT(INOUT),dimension(:)::S2

    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S2,DIM=1)
       N=LBOUND(S2,DIM=1)+K-1
    else
       I=LBOUND(S2,DIM=1)
       N=UBOUND(S2,DIM=1)
    endif

    DO   J=I,N
       call resetpoly_R31(s2(j))
    enddo


  END SUBROUTINE resetpoly_R31N


  SUBROUTINE  resetpoly0(S2)
    implicit none
    type (real_8),INTENT(INOUT)::S2

    if(s2%alloc) call killtpsa(s2%t)
    s2%alloc=f
    s2%kind=0
    s2%r=0.0_dp

    s2%i=0
    s2%s=1.0_dp

  END SUBROUTINE resetpoly0

  SUBROUTINE  K_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (real_8),INTENT(INout)::S1
    type (real_8),optional, INTENT(INout):: S2,s3,s4,s5,s6,s7,s8,s9,s10
    call resetpoly0(s1)
    if(present(s2)) call resetpoly0(s2)
    if(present(s3)) call resetpoly0(s3)
    if(present(s4)) call resetpoly0(s4)
    if(present(s5)) call resetpoly0(s5)
    if(present(s6)) call resetpoly0(s6)
    if(present(s7)) call resetpoly0(s7)
    if(present(s8)) call resetpoly0(s8)
    if(present(s9)) call resetpoly0(s9)
    if(present(s10))call resetpoly0(s10)
  END SUBROUTINE K_OPT








  !  SUBROUTINE  kill_c(k)
  !    implicit none
  !    integer ipause, mypauses
  !    integer k

  !    call dallsta(exi)   ! Etienne checking routine
  !    if(exi/=ent) then
  !       write(line,'(3(1x,i8))')  ent,exi,k
  !       ipause=mypauses(1999,line)
  !       k=k*10000
  !    endif

  !  END SUBROUTINE kill_c

  SUBROUTINE  resetpolyn0(S2,K)
    implicit none
    type (real_8),INTENT(INOUT),dimension(:)::S2
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S2,DIM=1)
       N=LBOUND(S2,DIM=1)+K-1
    else
       I=LBOUND(S2,DIM=1)
       N=UBOUND(S2,DIM=1)
    endif

    DO   J=I,N
       call resetpoly0(s2(j))
    enddo


  END SUBROUTINE resetpolyn0


  SUBROUTINE  allocpoly(S2)
    implicit none
    type (real_8),INTENT(INOUT)::S2

    !if(s2%alloc) call killtpsa(s2%t)     ADDED ETIENNE
    s2%alloc=f
    s2%kind=1
    s2%r=0.0_dp
    s2%i=0
!    s2%g=0
!    s2%nb=0
    s2%s=1.0_dp

  END SUBROUTINE allocpoly

  SUBROUTINE  allocquaternion(S2)
    implicit none
    type (quaternion_8),INTENT(INOUT)::S2
    integer i
     do i=0,3
      call allocpoly(s2%x(i))
    enddo

  END SUBROUTINE allocquaternion

  SUBROUTINE  killquaternion(S2)
    implicit none
    type (quaternion_8),INTENT(INOUT)::S2
    integer i
     do i=0,3
      call resetpoly0(s2%x(i))
    enddo

  END SUBROUTINE killquaternion

  SUBROUTINE  allocquaternionn(S2)
    implicit none
    type (quaternion_8),INTENT(INOUT)::S2(:)
    integer i
     do i=1,size(s2)
      call alloc(s2(i))
    enddo

  END SUBROUTINE allocquaternionn
 
  SUBROUTINE  killquaternionn(S2)
    implicit none
    type (quaternion_8),INTENT(INOUT)::S2(:)
    integer i
     do i=1,size(s2)
      call kill(s2(i))
    enddo

  END SUBROUTINE killquaternionn

  SUBROUTINE  A_OPT_quaternion(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (quaternion_8),INTENT(INout)::S1
    type (quaternion_8),optional, INTENT(INout):: S2,s3,s4,s5,s6,s7,s8,s9,s10
    call allocquaternion(s1)
    if(present(s2)) call allocquaternion(s2)
    if(present(s3)) call allocquaternion(s3)
    if(present(s4)) call allocquaternion(s4)
    if(present(s5)) call allocquaternion(s5)
    if(present(s6)) call allocquaternion(s6)
    if(present(s7)) call allocquaternion(s7)
    if(present(s8)) call allocquaternion(s8)
    if(present(s9)) call allocquaternion(s9)
    if(present(s10))call allocquaternion(s10)

  END SUBROUTINE A_OPT_quaternion

  SUBROUTINE  K_OPT_quaternion(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (quaternion_8),INTENT(INout)::S1
    type (quaternion_8),optional, INTENT(INout):: S2,s3,s4,s5,s6,s7,s8,s9,s10
    call killquaternion(s1)
    if(present(s2)) call killquaternion(s2)
    if(present(s3)) call killquaternion(s3)
    if(present(s4)) call killquaternion(s4)
    if(present(s5)) call killquaternion(s5)
    if(present(s6)) call killquaternion(s6)
    if(present(s7)) call killquaternion(s7)
    if(present(s8)) call killquaternion(s8)
    if(present(s9)) call killquaternion(s9)
    if(present(s10))call killquaternion(s10)

  END SUBROUTINE K_OPT_quaternion

  SUBROUTINE  A_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (real_8),INTENT(INout)::S1
    type (real_8),optional, INTENT(INout):: S2,s3,s4,s5,s6,s7,s8,s9,s10
    call allocpoly(s1)
    if(present(s2)) call allocpoly(s2)
    if(present(s3)) call allocpoly(s3)
    if(present(s4)) call allocpoly(s4)
    if(present(s5)) call allocpoly(s5)
    if(present(s6)) call allocpoly(s6)
    if(present(s7)) call allocpoly(s7)
    if(present(s8)) call allocpoly(s8)
    if(present(s9)) call allocpoly(s9)
    if(present(s10))call allocpoly(s10)
  END SUBROUTINE A_opt

  ! SUBROUTINE  alloc_c
  !   implicit none
  ! Etienne checking routine
  !   call dallsta(ent)


  !  END SUBROUTINE alloc_c

  SUBROUTINE  allocpolyn(S2,K)
    implicit none
    type (real_8),INTENT(INOUT),dimension(:)::S2
    INTEGER,optional,INTENT(IN)::k
    INTEGER J,i,N

    if(present(k)) then
       I=LBOUND(S2,DIM=1)
       N=LBOUND(S2,DIM=1)+K-1
    else
       I=LBOUND(S2,DIM=1)
       N=UBOUND(S2,DIM=1)
    endif

    DO   J=I,N
       call allocpoly(s2(j))
    enddo


  END SUBROUTINE allocpolyn


  subroutine init_map_p(NO1,ND1,NP1,NDPT1,log)
    implicit none
    integer NO1,ND1,NP1,NDPT1
    logical(lp) log
    call init_map_c(NO1,ND1,NP1,NDPT1,log)
    call set_in_polyp(log)
  end subroutine  init_map_p

  subroutine init_tpsa_p(NO1,NP1,log)
    implicit none
    integer NO1,NP1
    logical(lp) log
    call init_tpsa_c(NO1,NP1,log)
    call set_in_polyp(log)
  end subroutine  init_tpsa_p


  subroutine set_in_polyp(log)
    implicit none
    logical(lp) log
    integer iia(4),icoast(4)
    call liepeek(iia,icoast)
    old=log
    NO=iia(1)
    ND=iia(3)
    ND2=iia(3)*2
    NP=iia(2)-nd2
    NDPT=icoast(4)
    NV=iia(2)
    !    i_ =cmplx(zero,one,kind=dp)
  end  subroutine set_in_polyp

  SUBROUTINE  EQUAL2D(S2,S1)
    implicit none
    type (real_8),INTENT(inOUT)::S2(:,:)
    type (real_8),INTENT(IN)::S1(:,:)
    integer i,J,I1(2),I2(2)

    I1(1)=LBOUND(S1,DIM=1)
    I1(2)=LBOUND(S1,DIM=2)
    I2(1)=LBOUND(S2,DIM=1)
    I2(2)=LBOUND(S2,DIM=2)
    do i=I1(1),UBOUND(S1,DIM=1)
       do J=I1(2),UBOUND(S1,DIM=2)
          S2(I-I1(1)+I2(1),J-I1(2)+I2(2))=S1(I,J)
       ENDDO
    ENDDO

  end SUBROUTINE  EQUAL2D

  SUBROUTINE  EQUAL1D(S2,S1)
    implicit none
    type (real_8),INTENT(inOUT)::S2(:)
    type (real_8),INTENT(IN)::S1(:)
    integer i,I1,I2

    I1=LBOUND(S1,DIM=1)
    I2=LBOUND(S2,DIM=1)
    do i=I1,UBOUND(S1,DIM=1)
       S2(I-I1+I2)=S1(I)
    ENDDO

  end SUBROUTINE  EQUAL1D

  SUBROUTINE  EQUALq(S2,S1)
    implicit none
    integer ipause, mypauses
    type (quaternion_8),INTENT(inOUT)::S2
    type (quaternion_8),INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=s1%x(i)
      enddo

 end   SUBROUTINE  EQUALq

  SUBROUTINE  EQUALq_r_8(S2,S1)
    implicit none
    integer ipause, mypauses
    type (quaternion),INTENT(inOUT)::S2
    type (quaternion_8),INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=s1%x(i)
      enddo

 end   SUBROUTINE  EQUALq_r_8

  SUBROUTINE  EQUALq_8_r(S2,S1)
    implicit none
    integer ipause, mypauses
    type (quaternion_8),INTENT(inOUT)::S2
    type (quaternion),INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=s1%x(i)
      enddo

 end   SUBROUTINE  EQUALq_8_r

  SUBROUTINE  EQUALq_r(S2,S1)
    implicit none
    integer ipause, mypauses
    type (quaternion_8),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1
    integer i

      do i=0,3
        s2%x(i)=0.0_dp
      enddo
        s2%x(1)=s1

 end   SUBROUTINE  EQUALq_r

  SUBROUTINE  EQUAL(S2,S1)
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(inOUT)::S2
    type (real_8),INTENT(IN)::S1
    !    integer localmaster


    if(s1%kind==0) then
       line=  " You are putting kind=0  into something"
       ipause=mypauses(10,line)
    endif
    if(s2%kind==3.and.(.not.setknob)) then
       line=  " You are putting something  into a knob kind=3"
       ipause=mypauses(11,line)
    endif

    if (s2%kind>0) then       !   S2 exist
       if(S2%kind==S1%kind) then
          select case(S1%kind)
          case(m1)
             !        localmaster=master  ! added because of knob problems 2001.9.19
             !        call check_snake    ! added because of knob problems 2001.9.19
             !        master=0
             S2%R=S1%R
             !        master=localmaster  ! added because of knob problems 2001.9.19
          case(m2)
             !             localmaster=master
             call check_snake
             !  2002.12.26 master=0
             S2%t=S1%t
             !             master=localmaster
          case(m3)
             s2%r=S1%r  ! Knob stays a knob and real stays real 2002.10.9
             ! line=  " You are putting kind=3 (TPSA) into another kind=3"
             ! ipause=mypauses(13,line)
          end select
       elseif(S2%kind>S1%kind ) then
          if(S1%kind/=2) then
             s2%r=S1%r
             if(s2%kind/=3) s2%kind=1 !Knob stays a knob and real stays real 2002.10.9
          else
             s2%r=S1%t.sub.'0'  ! setting a kind=3
          endif
       else
          select case(s1%kind)
          case(m2)
             if(.not.s2%alloc) then
                call alloc(s2%t)
                s2%alloc=t
             endif
             s2%kind=2
             !             localmaster=master
             call check_snake
             !  2002.12.26 master=0

             S2%t=S1%t

             !             master=localmaster
          case(m3)
             if(.not.s2%alloc) then
                call alloc(s2%t)
                s2%alloc=t
             endif
             s2%kind=2
             !             localmaster=master
             call check_snake
             !  2002.12.26 master=0
             if(knob) then
                call varfk1(S1)
                S2%t=varf1
             else
                S2%R=S1%R
                s2%kind=1
             endif
             !             master=localmaster
          end select
       endif

    else        !   S2 does not exist



       if(S1%kind==1) then  ! what is s1
          if(s2%i==0) then
             S2%R=S1%R
             s2%kind=1
          elseif(s2%i>0.and.s2%i<=nv)  then
             call alloc(s2%t)
             s2%t=(/S1%R,S2%S/).var.s2%i
             !             call var(s2%t,S1%R,S2%S,s2%i)
             !      s2%i=0
             s2%kind=2
             s2%alloc=t
          else
             write(6,*) "EQUAL IN m_POLYMORPH ",s2%i,S2%KIND,S2%ALLOC
             WRITE(6,*) " I "
             READ(5,*) s2%i
             WRITE(6,*) SQRT(FLOAT(s2%i))
             stop 777
          endif

       else
          if(insane_PTC) then
             if(.not.s2%alloc) then
                call alloc(s2%t)
             endif
             !             localmaster=master
             call check_snake
             !  2002.12.26 master=0
             S2%t=S1%t
             !             master=localmaster
             s2%kind=2
             s2%alloc=t
          else
             !w_p=0
             !w_p%nc=5
             !w_p%fc='(4(1X,A72,/),(1X,A72))'
             write(6,'(A23,I4,A19)') " You are putting kind= ", s1%kind," (TPSA) in a kind=0"
               write(6,*)  " We do not allow that anymore for safety reasons"
             !w_p%c(3)=  " If you insist on it, modify real_polymorph and complex_polymorph"
             !w_p%c(4)= " at your own insane risk "
             !w_p%c(5)= " Etienne Forest/Frank Schmidt"
             ! call !write_e(778)
             !w_p%nc=sqrt(-float(!w_p%nc))
          endif
       endif    ! end of what is s1

    endif   ! S2  does not exist

  END SUBROUTINE EQUAL

  SUBROUTINE  complexreal_8(S2,S1)
    implicit none
    complex(dp) ,INTENT(inOUT)::S2
    type (real_8),INTENT(IN)::S1
    !    integer localmaster


    select case(S1%kind)
    case(m1)
       S2=S1%R
    case(m2)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%t.sub.'0'
       !       master=localmaster
    case(m3)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%r
       !       master=localmaster
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in complexreal_8 "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END SUBROUTINE complexreal_8



  SUBROUTINE  realEQUAL(S2,S1)
    implicit none
    real(dp),INTENT(inOUT)::S2
    type (real_8),INTENT(IN)::S1
    !    integer localmaster


    select case(S1%kind)
    case(m1)
       S2=S1%R
    case(m2)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%t.sub.'0'
       !       master=localmaster
    case(m3)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%r
       !       master=localmaster
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in realEQUAL "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END SUBROUTINE realEQUAL

  SUBROUTINE  singleEQUAL(S2,S1)
    implicit none
    real(sp),INTENT(inOUT)::S2
    type (real_8),INTENT(IN)::S1
    !    integer localmaster


    select case(S1%kind)
    case(m1)
       S2=S1%R
    case(m2)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%t.sub.'0'
       !       master=localmaster
    case(m3)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%r
       !       master=localmaster
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in realEQUAL "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END SUBROUTINE singleEQUAL



  SUBROUTINE  taylorEQUAL(S2,S1)
    implicit none
    type (taylor),INTENT(inOUT)::S2
    type (real_8),INTENT(IN)::S1
    !    integer localmaster


    select case(S1%kind)
    case(m1)
       S2=S1%R
    case(m2)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%t
       !       master=localmaster
    case(m3)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       if(knob) then
          call varfk1(S1)
          S2=varf1
       else
          S2=S1%R
       endif
       !       master=localmaster
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in taylorEQUAL "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END SUBROUTINE taylorEQUAL

  SUBROUTINE  EQUALtaylor(S2,S1)
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(inOUT)::S2
    type (taylor),INTENT(IN)::S1
    !    integer localmaster

    IF(S2%KIND==3.and.(.not.setknob)) THEN  ! 2002.10.9
       line=  "Forbidden in EQUALtaylor: s2 is a knob "
       ipause=mypauses(0,line)
    ENDIF

    !    localmaster=master
    call check_snake
    !  2002.12.26 master=0
    if(s2%kind/=3) then
       if(.not.s2%alloc) then
          call alloc(s2%t)
          s2%alloc=t
       endif
       S2%t=S1
       s2%KIND=2
    else
       s2%r=S1.sub.'0' ! 2002.10.9
    endif
    !    master=localmaster

  END SUBROUTINE EQUALtaylor

  SUBROUTINE  real_8univ(S2,S1)
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(inOUT)::S2
    type (universal_taylor),INTENT(IN)::S1
    !    integer localmaster

    IF(S2%KIND==M3) THEN
       line=  "Forbidden in real_8univ: s2 is a knob "
       ipause=mypauses(0,line)
    ENDIF

    !    localmaster=master
    !    call check_snake
    !  2002.12.26 master=0
    if(.not.s2%alloc) then
       call alloc(s2%t)
       s2%alloc=t
    endif
    S2%t=S1
    s2%KIND=2

    !    master=localmaster
  END SUBROUTINE real_8univ

  SUBROUTINE  univreal_8(S2,S1)   ! new sagan
    implicit none
    type (universal_taylor),INTENT(inOUT)::S2
    type (real_8),INTENT(IN)::S1
    !    integer localmaster

    select case(S1%kind)
    case(m1)
       S2=S1%R
    case(m2)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       S2=S1%t
       !       master=localmaster
    case(m3)
       !       localmaster=master
       call check_snake
       !  2002.12.26 master=0
       if(knob) then
          call varfk1(S1)
          S2=varf1
       else
          S2=S1%R
       endif
       !       master=localmaster
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in univreal_8 "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END SUBROUTINE univreal_8

  SUBROUTINE  Mequaldacon(S2,R1)
    implicit none
    type (real_8),INTENT(inOUT)::S2(:)
    real(dp),INTENT(IN)::R1
    integer i
    do i=1,size(s2)
      s2(i)=0.0_dp
    enddo
    end SUBROUTINE  Mequaldacon
  
  SUBROUTINE  Dequaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::R1

    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=r1
          return
       else
          line=  "Forbidden in Dequaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF


    !localmaster=master
    ! call check_snake
    !  master=0
    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1.or.s2%kind==3) then
          S2%R=R1
       else
          !     S2%t=R1
          s2%r=r1
          s2%kind=1
          !       call kill(s2%t)
          !       s2%alloc=f
       endif

    else        !   S2 does not exist



       if(s2%i==0) then
          S2%R=R1
          s2%kind=1
       elseif(s2%i>0.and.s2%i<=nv)  then
          call alloc(s2%t)
          s2%t=(/R1,s2%S/).var.s2%i
          !          call var(s2%t,R1,s2%S,s2%i)
          !      s2%i=0
          s2%kind=2
          s2%alloc=t
       else
          line=  "trouble in Dequaldacon in Real_polymorph"
          ipause=mypauses(-788,line)
       endif


    endif   ! S2 not allocated
    !   master=localmaster
  END SUBROUTINE Dequaldacon

  SUBROUTINE  iequaldaconn(S2,R1)
    implicit none
    type (real_8),INTENT(inOUT),dimension(:)::S2
    integer,INTENT(IN)::R1
    integer i
    !    integer localmaster

    !    localmaster=master
    !    call check_snake
    !  2002.12.26 master=0
    do i=1,r1
       call resetpoly(s2(i))
       s2(i)%i=i
    enddo
    !    master=localmaster
  END SUBROUTINE iequaldaconn


  SUBROUTINE  equaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(inOUT)::S2
    real(sp),INTENT(IN)::R1

    if(real_warning) call real_stop
    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=REAL(R1,kind=DP)
          return
       else
          line=  "Forbidden in equaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF
    !localmaster=master
    ! call check_snake
    !  master=0
    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1.or.s2%kind==3) then
          S2%R=REAL(R1,kind=DP)
       else
          !     S2%t=REAL(R1,kind=DP)
          !     S2%t=R1
          s2%r=REAL(R1,kind=DP)
          s2%kind=1
          !       call kill(s2%t)
          !       s2%alloc=f
       endif

    else        !   S2 does not exist



       if(s2%i==0) then
          S2%R=REAL(R1,kind=DP)
          s2%kind=1
       elseif(s2%i>0.and.s2%i<=nv)  then
          call alloc(s2%t)
          s2%t=(/REAL(R1,kind=DP),S2%S/).var.s2%i
          !          call var(s2%t,REAL(R1,kind=DP),S2%S,s2%i)
          !      s2%i=0
          s2%kind=2
          s2%alloc=t
       else
          stop 779
       endif


    endif   ! S2 not allocated
    !   master=localmaster
  END SUBROUTINE equaldacon

  SUBROUTINE  iequaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (real_8),INTENT(inOUT)::S2
    integer,INTENT(IN)::R1

    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=r1
          return
       else
          line=  "Forbidden in iequaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF
    !localmaster=master
    ! call check_snake
    !  master=0
    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1.or.s2%kind==3) then
          S2%R=REAL(R1,kind=DP)
       else
          !   S2%t=REAL(R1,kind=DP)
          s2%r=REAL(R1,kind=DP)
          s2%kind=1
          !       call kill(s2%t)
          !       s2%alloc=f
       endif

    else        !   S2 does not exist



       if(s2%i==0) then
          S2%R=REAL(R1,kind=DP)
          s2%kind=1
       elseif(s2%i>0.and.s2%i<=nv)  then
          call alloc(s2%t)
          s2%t=(/REAL(R1,kind=DP),S2%S/).var.s2%i
          !          call var(s2%t,REAL(R1,kind=DP),S2%S,s2%i)
          !     s2%i=0
          s2%kind=2
          s2%alloc=t
       else
          stop 776
       endif


    endif   ! S2 not allocated
    !   master=localmaster
  END SUBROUTINE iequaldacon

  subroutine assp(s1)
    implicit none
    TYPE (real_8) s1
    integer ipause,mypauses
    ! lastmaster=master  ! 2002.12.13

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)
       line=  " cannot indent anymore "
       ipause=mypauses(0,line)
    end select
    !    write(26,*) " real polymorph  ",master

    call ass0(s1%t)
    s1%alloc=t
    s1%kind=2
    s1%i=0

  end subroutine ASSp


  subroutine assp_no_master(s1)
    implicit none
    TYPE (real_8) s1


    call ass0(s1%t)
    s1%alloc=t
    s1%kind=1
    s1%i=0

  end subroutine assp_no_master

 

  FUNCTION datant( S1 )
    implicit none
    TYPE (real_8) datant
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complextaylor) w
    integer localmaster

    select case(s1%kind)
    case(m1)
       datant%r=atan(s1%r)
       datant%kind=1
    case(m2)
       localmaster=master
       call ass(datant)
       call alloc(w)
       w%r=s1%t
       w= atan(w)
       datant%t=w%r
       call kill(w)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(datant)
          call alloc(w)
          call varfk1(S1)
          w%r=varf1
          w= atan(w)
          datant%t=w%r
          call kill(w)
          master=localmaster
       else
          datant%r=atan(s1%r)
          datant%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in datant "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION datant

  FUNCTION datanht( S1 )
    implicit none
    TYPE (real_8) datanht
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       datanht%r=log((1+s1%r)/(1-s1%r))/2.0_dp
       datanht%kind=1
    case(m2)
       localmaster=master
       call ass(datanht)
       datanht%t=atanh(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(datanht)
          call varfk1(S1)
          datanht%t=atanh(varf1)
          master=localmaster
       else
          datanht%r=log((1+s1%r)/sqrt(1-s1%r))/2.0_dp
          datanht%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in datant "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION datanht





  FUNCTION datanDt( S1 )
    implicit none
    TYPE (real_8) datanDt
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complextaylor) w
    integer localmaster

    select case(s1%kind)
    case(m1)
       datanDt%r=atan(s1%r)*RAD_TO_DEG_
       datanDt%kind=1
    case(m2)
       localmaster=master
       call ass(datanDt)
       call alloc(w)
       w%r=s1%t
       w= atan(w)
       datanDt%t=w%r*RAD_TO_DEG_
       call kill(w)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(datanDt)
          call alloc(w)
          call varfk1(S1)
          w%r=varf1
          w= atan(w)
          datanDt%t=w%r*RAD_TO_DEG_
          call kill(w)
          master=localmaster
       else
          datanDt%r=atan(s1%r)*RAD_TO_DEG_
          datanDt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in datanDt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION datanDt

  FUNCTION absoftdatanDr( S1 )
    implicit none
    real(dp) absoftdatanDr
    real(dp), INTENT (IN) :: S1
    absoftdatanDr=atan(s1)*RAD_TO_DEG_
  END FUNCTION absoftdatanDr

  FUNCTION datan2t( S2,S1 )
    implicit none
    integer ipause, mypauses
    TYPE (real_8) datan2t
    TYPE (real_8), INTENT (IN) :: S2,S1
    TYPE (complextaylor) w
    real(dp) ANG
    integer localmaster ,si
    si=1.0_dp
    if(s2%kind==0.or.S1%kind==0) then
       line=  " Problems in datan2t "
       ipause=mypauses(0,line)
    endif
    if(S2%kind==S1%kind) then  ! 1
       select case(S1%kind)
       case(m1)
          DATAN2T%R=ATAN2(S2%R,S1%R)
          datan2t%kind=1
       case(m2)
          localmaster=master
          call ass(datan2t)
          call alloc(w)
          ANG=ATAN2(S2%T.SUB.'0',S1%T.SUB.'0')
          if(pil>abs(ANG).or.abs(ANG)>pim) then
             w%r=s2%t/s1%t
             w=atan(w)
             datan2t%t=w%r
             datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG
          else
             if((s2%t.sub.'0')<0.0_dp) si=-1.0_dp
             w%r=s1%t/SQRT(s2%t**2+s1%t**2)
             w=acos(w)
             datan2t%t=w%r
             datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
          endif
          call kill(w)
          master=localmaster
       case(m3)
          localmaster=master
          if(knob) then        !knob 1
             call ass(datan2t)
             call alloc(w)
             call varfk1(S1)
             call varfk2(S2)
             ANG=ATAN2(varf2.SUB.'0',varf1.SUB.'0')
             if(pil>abs(ANG).or.abs(ANG)>pim) then
                w%r=varf2/varf1
                w=atan(w)
                datan2t%t=w%r
                datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG
             else
                if((varf2.sub.'0')<0.0_dp) si=-1.0_dp
                w%r=varf1/SQRT(varf2**2+varf1**2)
                w=acos(w)
                datan2t%t=w%r
                datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
             endif
             call kill(w)
          else             !knob 1
             DATAN2T%R=ATAN2(S2%R,S1%R)
             datan2t%kind=1
          endif              !knob 1
          master=localmaster
       end select
    elseif(S2%kind>S1%kind) then  !1
       localmaster=master
       !       call ass(datan2t)
       call alloc(w)
       if(s2%kind==2) then !
          call ass(datan2t)
          ANG=ATAN2(S2%T.SUB.'0',S1%R)
          if(pil>abs(ANG).or.abs(ANG)>pim) then
             w%r=s2%t/s1%r
             w=atan(w)
             datan2t%t=w%r
             datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG

          else
             if((s2%t.sub.'0')<0.0_dp) si=-1.0_dp
             w%r=s1%t/SQRT(s2%t**2+s1%r**2)
             w=acos(w)
             datan2t%t=w%r
             datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
          endif
       else   !
          if(knob) then ! knob
             call ass(datan2t)
             call varfk2(S2)
             if(s1%kind==1) then !!
                ANG=ATAN2(varf2.SUB.'0',S1%R)
                if(pil>abs(ANG).or.abs(ANG)>pim) then  !!!
                   w%r=varf2/s1%r
                   w=atan(w)
                   datan2t%t=w%r
                   datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG

                else    !!!
                   if((varf2.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=s1%R/SQRT(varf2**2+s1%R**2)   !   etienne error here????
                   w=acos(w)
                   datan2t%t=w%r
                   datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
                endif !!!

             else  !!
                ANG=ATAN2(varf2.SUB.'0',S1%t.sub.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then  !!!
                   w%r=varf2/s1%t
                   w=atan(w)
                   datan2t%t=w%r
                   datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG

                else    !!!
                   if((varf2.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=s1%t/SQRT(varf2**2+s1%t**2)
                   w=acos(w)
                   datan2t%t=w%r
                   datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
                endif !!!
             endif !!
          else ! knob

             if(s1%kind==1) then !!
                DATAN2T%R=ATAN2(S2%R,S1%R)
                datan2t%kind=1

             else  !!
                call ass(datan2t)
                ANG=ATAN2(S2%R,S1%t.sub.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then  !!!
                   w%r=S2%R/s1%t
                   w=atan(w)
                   datan2t%t=w%r
                   datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG

                else    !!!
                   if((S2%R)<0.0_dp) si=-1.0_dp
                   w%r=s1%t/SQRT(S2%R**2+s1%t**2)
                   w=acos(w)
                   datan2t%t=w%r
                   datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
                endif !!!
             endif !!


          endif ! knob
       endif  !

       call kill(w)
       master=localmaster
    ELSE   !  1
       localmaster=master
       !       call ass(datan2t)
       call alloc(w)
       if(s1%kind==2) then !
          call ass(datan2t)
          ANG=ATAN2(S2%R,S1%T.SUB.'0')
          if(pil>abs(ANG).or.abs(ANG)>pim) then
             w%r=s2%r/s1%t
             w=atan(w)
             datan2t%t=w%r
             datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG
          else
             if(s2%r<0.0_dp) si=-1.0_dp
             w%r=s1%t/SQRT(s2%r**2+s1%t**2)
             w=acos(w)
             datan2t%t=w%r
             datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
          endif
       else !
          IF(KNOB) THEN ! KNOB 2
             call ass(datan2t)
             if(s2%kind==1) then  !!
                call varfk1(S1)
                ANG=ATAN2(S2%R,varf1.SUB.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then !!!
                   w%r=s2%r/varf1
                   w=atan(w)
                   datan2t%t=w%r
                   datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG
                else  !!!
                   if(s2%r<0.0_dp) si=-1.0_dp
                   w%r=varf1/SQRT(s2%r**2+varf1**2)
                   w=acos(w)
                   datan2t%t=w%r
                   datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
                endif  !!!
             else    !! s2%kind ==2
                call varfk1(S1)
                ANG=ATAN2(S2%t.sub.'0',varf1.SUB.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then !!!
                   w%r=s2%t/varf1
                   w=atan(w)
                   datan2t%t=w%r
                   datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG
                else  !!!
                   if((S2%t.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=varf1/SQRT(s2%t**2+varf1**2)
                   w=acos(w)
                   datan2t%t=w%r
                   datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
                endif  !!!
             endif   !!

          ELSE   !KNOB 2
             if(s2%kind==1) then  !!
                DATAN2T%R=ATAN2(S2%R,S1%R)
                datan2t%kind=1
             else    !! s2%kind ==2
                call ass(datan2t)
                ANG=ATAN2(S2%t.sub.'0',S1%r)
                if(pil>abs(ANG).or.abs(ANG)>pim) then !!!
                   w%r=s2%t/S1%r
                   w=atan(w)
                   datan2t%t=w%r
                   datan2t%t=datan2t%t-(datan2t%t.SUB.'0')+ ANG
                else  !!!
                   if((S2%t.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=S1%r/SQRT(s2%t**2+S1%r**2)
                   w=acos(w)
                   datan2t%t=w%r
                   datan2t%t=si*(datan2t%t-(datan2t%t.SUB.'0'))+ ANG
                endif  !!!
             endif   !!




          ENDIF  ! KNOB 2


       endif !

       call kill(w)
       master=localmaster
    endif


  END FUNCTION datan2t

  FUNCTION datan2tt( S2,S1 )
    implicit none
 
    TYPE (taylor) datan2tt
    TYPE (taylor), INTENT (IN) :: S2,S1
    TYPE (real_8) temp,s1_8,s2_8
    integer localmaster 
    

    localmaster=master
    call ass(datan2tt)
     call alloc(temp,s1_8,s2_8)   

     s1_8=morph(S1)
     s2_8=morph(S2)

     temp=atan2(s2_8,s1_8)
    
     datan2tt=temp%t
     
    master=localmaster
    call kill(temp,s1_8,s2_8)

  END FUNCTION datan2tt

  FUNCTION datan2dt( S2,S1 )
    implicit none
    integer ipause, mypauses
    TYPE (real_8) datan2dt
    TYPE (real_8), INTENT (IN) :: S2,S1
    TYPE (complextaylor) w
    real(dp) ANG
    integer localmaster ,si
    si=1.0_dp
    if(s2%kind==0.or.S1%kind==0) then
       line=  " Problems in datand2t "
       ipause=mypauses(0,line)
    endif
    if(S2%kind==S1%kind) then
       select case(S1%kind)
       case(m1)
          datan2dt%R=ATAN2(S2%R,S1%R)*RAD_TO_DEG_
          datan2dt%kind=1
       case(m2)
          localmaster=master
          call ass(datan2dt)
          call alloc(w)
          ANG=ATAN2(S2%T.SUB.'0',S1%T.SUB.'0')
          if(pil>abs(ANG).or.abs(ANG)>pim) then
             w%r=s2%t/s1%t
             w=atan(w)
             datan2dt%t=w%r
             datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
             datan2dt%t=datan2dt%t*RAD_TO_DEG_
          else
             if((s2%t.sub.'0')<0.0_dp) si=-1.0_dp
             w%r=s1%t/SQRT(s2%t**2+s1%t**2)
             w=acos(w)
             datan2dt%t=w%r
             datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
             datan2dt%t=datan2dt%t*RAD_TO_DEG_
          endif
          call kill(w)
          master=localmaster
       case(m3)
          localmaster=master
          if(knob) then   ! knob 1
             call ass(datan2dt)
             call alloc(w)
             call varfk1(S1)
             call varfk2(S2)
             ANG=ATAN2(varf2.SUB.'0',varf1.SUB.'0')  !etienne
             if(pil>abs(ANG).or.abs(ANG)>pim) then
                w%r=varf2/varf1
                w=atan(w)
                datan2dt%t=w%r
                datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                datan2dt%t=datan2dt%t*RAD_TO_DEG_
             else
                if((varf2.sub.'0')<0.0_dp) si=-1.0_dp
                w%r=varf1/SQRT(varf2**2+varf1**2)
                w=acos(w)
                datan2dt%t=w%r
                datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                datan2dt%t=datan2dt%t*RAD_TO_DEG_
             endif
             call kill(w)

          else  ! knob 1

             datan2dt%R=ATAN2(S2%R,S1%R)*RAD_TO_DEG_
             datan2dt%kind=1

          endif  ! knob 1

          master=localmaster
       end select
    elseif(S2%kind>S1%kind) then
       localmaster=master
       !       call ass(datan2dt)
       call alloc(w)
       if(s2%kind==2) then !
          call ass(datan2dt)
          ANG=ATAN2(S2%T.SUB.'0',S1%R)
          if(pil>abs(ANG).or.abs(ANG)>pim) then
             w%r=s2%t/s1%r
             w=atan(w)
             datan2dt%t=w%r
             datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
             datan2dt%t=datan2dt%t*RAD_TO_DEG_
          else
             if((s2%t.sub.'0')<0.0_dp) si=-1.0_dp
             w%r=s1%t/SQRT(s2%t**2+s1%r**2)
             w=acos(w)
             datan2dt%t=w%r
             datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
             datan2dt%t=datan2dt%t*RAD_TO_DEG_
          endif
       else !
          if(knob) then      !knob 2
             call ass(datan2dt)
             call varfk2(S2)
             if(s1%kind==1) then !!
                ANG=ATAN2(varf2.SUB.'0',S1%R)
                if(pil>abs(ANG).or.abs(ANG)>pim) then   !!!
                   w%r=varf2/s1%r
                   w=atan(w)
                   datan2dt%t=w%r
                   datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                else !!!
                   if((varf2.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=s1%t/SQRT(varf2**2+s1%t**2)
                   w=acos(w)
                   datan2dt%t=w%r
                   datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                endif   !!!
             else !!
                ANG=ATAN2(varf2.SUB.'0',S1%t.sub.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then   !!!
                   w%r=varf2/s1%t
                   w=atan(w)
                   datan2dt%t=w%r
                   datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                else !!!
                   if((varf2.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=s1%t/SQRT(varf2**2+s1%t**2)
                   w=acos(w)
                   datan2dt%t=w%r
                   datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                endif   !!!
             endif   !!

          else     !knob 2

             if(s1%kind==1) then !!
                datan2dt%R=ATAN2(S2%R,S1%R)*RAD_TO_DEG_
                datan2dt%kind=1
             else !!
                call ass(datan2dt)
                ANG=ATAN2(S2%R,S1%t.sub.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then   !!!
                   w%r=S2%R/s1%t
                   w=atan(w)
                   datan2dt%t=w%r
                   datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                else !!!
                   if((S2%R)<0.0_dp) si=-1.0_dp
                   w%r=s1%t/SQRT(S2%R**2+s1%t**2)
                   w=acos(w)
                   datan2dt%t=w%r
                   datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                endif   !!!
             endif   !!

          endif     !knob 2


       endif !
       call kill(w)
       master=localmaster
    ELSE  !1
       localmaster=master
       !       call ass(datan2dt)
       call alloc(w)
       if(s1%kind==2) then !
          call ass(datan2dt)
          ANG=ATAN2(S2%R,S1%T.SUB.'0')
          if(pil>abs(ANG).or.abs(ANG)>pim) then
             w%r=s2%r/s1%t
             w=atan(w)
             datan2dt%t=w%r
             datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
             datan2dt%t=datan2dt%t*RAD_TO_DEG_
          else
             if(s2%r<0.0_dp) si=-1.0_dp
             w%r=s1%t/SQRT(s2%r**2+s1%t**2)
             w=acos(w)
             datan2dt%t=w%r
             datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
             datan2dt%t=datan2dt%t*RAD_TO_DEG_
          endif
       else !
          if(knob) then  ! knob 3
             call ass(datan2dt)
             if(s2%kind==1) then  !!
                call varfk1(S1)
                ANG=ATAN2(S2%R,varf1.SUB.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then
                   w%r=s2%r/varf1
                   w=atan(w)
                   datan2dt%t=w%r
                   datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                else
                   if(s2%r<0.0_dp) si=-1.0_dp
                   w%r=varf1/SQRT(s2%r**2+varf1**2)
                   w=acos(w)
                   datan2dt%t=w%r
                   datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                endif
             else    !! s2%kind ==2
                call varfk1(S1)
                ANG=ATAN2(S2%t.sub.'0',varf1.SUB.'0')
                if(pil>abs(ANG).or.abs(ANG)>pim) then
                   w%r=s2%T/varf1
                   w=atan(w)
                   datan2dt%t=w%r
                   datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                else
                   if((S2%t.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=varf1/SQRT(s2%T**2+varf1**2)
                   w=acos(w)
                   datan2dt%t=w%r
                   datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                endif
             endif   !!
          else  ! knob 3
             if(s2%kind==1) then  !!
                datan2dt%R=ATAN2(S2%R,S1%R)*RAD_TO_DEG_
                datan2dt%kind=1
             else    !! s2%kind ==2
                call ass(datan2dt)
                ANG=ATAN2(S2%t.sub.'0',S1%r)
                if(pil>abs(ANG).or.abs(ANG)>pim) then
                   w%r=s2%T/S1%r
                   w=atan(w)
                   datan2dt%t=w%r
                   datan2dt%t=datan2dt%t-(datan2dt%t.SUB.'0')+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                else
                   if((S2%t.sub.'0')<0.0_dp) si=-1.0_dp
                   w%r=S1%r/SQRT(s2%T**2+S1%r**2)
                   w=acos(w)
                   datan2dt%t=w%r
                   datan2dt%t=si*(datan2dt%t-(datan2dt%t.SUB.'0'))+ ANG
                   datan2dt%t=datan2dt%t*RAD_TO_DEG_
                endif
             endif   !!
          endif ! knob 3
       endif !

       call kill(w)
       master=localmaster
    endif

  END FUNCTION datan2dt

  FUNCTION absoftdatan2dr( S2,S1 )
    implicit none
    real(dp) absoftdatan2dr
    real(dp), INTENT (IN) :: S2,S1
    absoftdatan2dr=ATAN2(S2,S1)*RAD_TO_DEG_
  END FUNCTION absoftdatan2dr

  FUNCTION dasint( S1 )
    implicit none
    TYPE (real_8) dasint
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (taylor) w
    integer localmaster

    select case(s1%kind)
    case(m1)
       dasint%r=arcsin(s1%r)
       !       dasint%r=asin(s1%r)
       dasint%kind=1
    case(m2)
       localmaster=master
       call ass(dasint)
       call alloc(w)
       w=s1%t
       w= asin(w)
       dasint%t=w
       call kill(w)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dasint)
          call alloc(w)
          call varfk1(S1)
          w=varf1
          w= asin(w)
          dasint%t=w
          call kill(w)
          master=localmaster
       else
          dasint%r=asin(s1%r)
          dasint%kind=1
       endif
       !       call unvarkind3(S1)
    end select
  END FUNCTION dasint

  FUNCTION dacost( S1 )
    implicit none
    TYPE (real_8) dacost
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (taylor) w
    integer localmaster

    select case(s1%kind)
    case(m1)
       dacost%r=arccos(s1%r)
       dacost%kind=1
    case(m2)
       localmaster=master
       call ass(dacost)
       call alloc(w)
       w=s1%t
       w= acos(w)
       dacost%t=w
       call kill(w)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dacost)
          call alloc(w)
          call varfk1(S1)
          w=varf1
          w= acos(w)
          dacost%t=w
          call kill(w)
          master=localmaster
       else
          dacost%r=acos(s1%r)
          dacost%kind=1
       endif

    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dacost "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dacost


  FUNCTION dcosht( S1 )
    implicit none
    TYPE (real_8) dcosht
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complextaylor) w
    integer localmaster

    select case(s1%kind)
    case(m1)
       dcosht%r=cosh(s1%r)
       dcosht%kind=1
    case(m2)
       localmaster=master
       call ass(dcosht)
       call alloc(w)
       w%r=s1%t
       w= cosh(w)
       dcosht%t=w%r
       call kill(w)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dcosht)
          call alloc(w)
          call varfk1(S1)
          w%r=varf1
          w= cosh(w)
          dcosht%t=w%r
          call kill(w)
          master=localmaster
       else
          dcosht%r=cosh(s1%r)
          dcosht%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dcosht "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dcosht

  FUNCTION dsinht( S1 )
    implicit none
    TYPE (real_8) dsinht
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complextaylor) w
    integer localmaster

    select case(s1%kind)
    case(m1)
       dsinht%r=sinh(s1%r)
       dsinht%kind=1
    case(m2)
       localmaster=master
       call ass(dsinht)
       call alloc(w)
       w%r=s1%t
       w= sinh(w)
       dsinht%t=w%r
       call kill(w)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dsinht)
          call alloc(w)
          call varfk1(S1)
          w%r=varf1
          w= sinh(w)
          dsinht%t=w%r
          call kill(w)
          master=localmaster
       else
          dsinht%r=sinh(s1%r)
          dsinht%kind=1
       endif

    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dsinht "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dsinht

  FUNCTION dtanht( S1 )
    implicit none
    TYPE (real_8) dtanht
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dtanht%r=tanh(s1%r)
       dtanht%kind=1
    case(m2)
       localmaster=master
       call ass(dtanht)
       dtanht%t=tanh(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(dtanht)
          call varfk1(S1)
          dtanht%t=tanh(varf1)
          master=localmaster
       else
          dtanht%r=tanh(s1%r)
          dtanht%kind=1
       endif

    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dtanht "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dtanht

  SUBROUTINE  MAPreal_8(S2,S1)
    implicit none
    type (damap),INTENT(inout)::S2
    type (real_8),INTENT(IN),dimension(:)::S1
    integer i

    call check_snake

    do i=1,nd2
       s2%v(i)=s1(i)          !%t
    enddo
  END SUBROUTINE MAPreal_8

  SUBROUTINE  normal_p(S2,S1)
    implicit none
    type (normalform),INTENT(inOUT)::S2
    type (real_8),INTENT(IN),dimension(:)::S1
    type (damap)  t
    integer i
    call alloc(t)

    call check_snake

    do i=1,nd2
       t%v(i)=s1(i)          !%t
    enddo
    s2=t

    call kill(t)
  END SUBROUTINE normal_p


  SUBROUTINE  real_8MAP(S1,S2)
    implicit none
    type (damap),INTENT(in)::S2
    type (real_8),INTENT(inout),dimension(:)::S1
    integer i

    call check_snake

    do i=1,nd2
       s1(i)=s2%v(i)
    enddo

  END SUBROUTINE real_8MAP




  SUBROUTINE  varfk1(S2)
    implicit none
    !    type (real_8),INTENT(INOUT)::S2
    type (real_8),INTENT(IN)::  S2

    if(knob) then
    !   if(nb_==0) then
          varf1=(/S2%R,S2%S/).var.(s2%i+npara_fpp)
    !   elseif(s2%nb==nb_) then
    !      varf1=(/S2%R,S2%S/).var.(s2%i+npara_fpp-s2%g+1)
    !   else
    !      varf1=S2%R
    !   endif
    else ! Not a knob
       stop 333
       varf1=(/S2%R,0.0_dp/).var.0  ! this is a buggy line never used
    endif


  end SUBROUTINE  varfk1

  SUBROUTINE  varfk2(S2)
    implicit none
    !    type (real_8),INTENT(INOUT)::S2
    type (real_8),INTENT(IN)::  S2

    if(knob) then
  !     if(nb_==0) then
          varf2=(/S2%R,S2%S/).var.(s2%i+npara_fpp)
  !     elseif(s2%nb==nb_) then
  !        varf2=(/S2%R,S2%S/).var.(s2%i+npara_fpp-s2%g+1)
  !     else
  !        varf2=S2%R
  !     endif
    else ! Not a knob
       stop 334
       varf2=(/S2%R,0.0_dp/).var.0   ! this is a buggy line never used
    endif


  end SUBROUTINE  varfk2


  ! SINH(X)/X DONE partly COSY-INFINITY WISE FOR TPSA
  function SINH_HR(X)
    implicit none
    real(dp), INTENT (IN) :: X
    real(dp) SINH_HR,Y,NORM0,NORM,SINH_HR0
    logical(lp) NOTDONE,CHECK
    INTEGER I,ipause,mypauses

    if(abs(x)<sinhx_x_min) then
       NOTDONE=.TRUE.
       CHECK=.TRUE.

       SINH_HR=1.0_dp
       Y=1.0_dp
       I=1
       NORM0=1e5_dp
       DO WHILE(I<NMAX_pol.AND.NOTDONE)
          Y=Y*X**2/REAL(I+1,kind=DP)/REAL(I+2,kind=DP)
          SINH_HR0=SINH_HR+Y
          NORM=ABS(SINH_HR-SINH_HR0)
          IF(NORM<=EPS_real_poly.AND.CHECK) THEN
             NORM0=NORM
             CHECK=.FALSE.
          ELSE
             IF(NORM>=NORM0) THEN
                NOTDONE=.FALSE.
             ELSE
                NORM0=NORM
             ENDIF
          ENDIF
          SINH_HR=SINH_HR0
          I=I+2
       ENDDO
       IF(I==NMAX_pol) THEN
          line="NO CONVERGENCE IN SINH_HR"
          ipause=mypauses(NMAX_pol,line)
       ENDIF
    else
       SINH_HR=SINH(X)/X
    endif
    return
  end function SINH_HR


  ! SIN(X)/X DONE partly COSY-INFINITY WISE FOR TPSA
  function SIN_HR(X)
    implicit none
    real(dp), INTENT (IN) :: X
    real(dp) SIN_HR,Y,NORM0,NORM,SINH_HR0
    logical(lp) NOTDONE,CHECK
    INTEGER I,ipause,mypauses

    if(abs(x)<sinhx_x_min) then
       NOTDONE=.TRUE.
       CHECK=.TRUE.

       SIN_HR=1.0_dp
       Y=1.0_dp
       I=1
       NORM0=1e5_dp
       DO WHILE(I<NMAX_pol.AND.NOTDONE)
          Y=-Y*X**2/REAL(I+1,kind=DP)/REAL(I+2,kind=DP)
          SINH_HR0=SIN_HR+Y
          NORM=ABS(SIN_HR-SINH_HR0)
          IF(NORM<=EPS_real_poly.AND.CHECK) THEN
             NORM0=NORM
             CHECK=.FALSE.
          ELSE
             IF(NORM>=NORM0) THEN
                NOTDONE=.FALSE.
             ELSE
                NORM0=NORM
             ENDIF
          ENDIF
          SIN_HR=SINH_HR0
          I=I+2
       ENDDO
       IF(I==NMAX_pol) THEN
          line="NO CONVERGENCE IN SINH_HR"
          ipause=mypauses(NMAX_pol,line)
       ENDIF
    else
       SIN_HR=SIN(X)/X
    endif
    return
  end function SIN_HR

  FUNCTION sinX_Xt( S1 )
    implicit none
    integer ipause, mypauses
    TYPE (real_8) sinX_Xt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster
    TYPE(TAYLOR)SIN_HP
    TYPE(TAYLOR)Y,SINH_HP0
    real(dp) NORM0,NORM
    logical(lp) NOTDONE,CHECK
    INTEGER I

    select case(s1%kind)
    case(m1)
       sinX_Xt%r=SIN_HR(s1%r)
       sinX_Xt%kind=1
    case(m2)
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
       localmaster=master
       call ass(sinX_Xt)

       if(ABS(S1%T.SUB.'0')<sinhx_x_minp) then
          NOTDONE=.TRUE.
          CHECK=.TRUE.

          SIN_HP=1.0_dp
          Y=1.0_dp
          I=1
          NORM0=1e5_dp
          DO WHILE(I<NMAX_pol.AND.NOTDONE)
             Y=-Y*S1%T**2/REAL(I+1,kind=DP)/REAL(I+2,kind=DP)
             SINH_HP0=SIN_HP+Y
             NORM=full_abs(SIN_HP-SINH_HP0)
             IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                NORM0=NORM
                CHECK=.FALSE.
             ELSE
                IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                   NOTDONE=.FALSE.
                ELSE
                   NORM0=NORM
                ENDIF
             ENDIF
             SIN_HP=SINH_HP0
             I=I+2
          ENDDO
          IF(I==NMAX_pol) THEN
             line="NO CONVERGENCE IN SIN_HP"
             ipause=mypauses(NMAX_pol,line)
          ENDIF
       else
          SIN_HP=SIN(S1%T)/S1%T
       endif
       sinX_Xt%t= SIN_HP

       master=localmaster
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);
    case(m3)
       if(knob) then
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
          localmaster=master
          call ass(sinX_Xt)
          call varfk1(S1)
          if(ABS(varf1.SUB.'0')<sinhx_x_minp) then
             NOTDONE=.TRUE.
             CHECK=.TRUE.

             SIN_HP=1.0_dp
             Y=1.0_dp
             I=1
             NORM0=1e5_dp
             DO WHILE(I<NMAX_pol.AND.NOTDONE)
                Y=-Y*varf1**2/REAL(I+1,kind=DP)/REAL(I+2,kind=DP)
                SINH_HP0=SIN_HP+Y
                NORM=full_abs(SIN_HP-SINH_HP0)
                IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                   NORM0=NORM
                   CHECK=.FALSE.
                ELSE
                   IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                      NOTDONE=.FALSE.
                   ELSE
                      NORM0=NORM
                   ENDIF
                ENDIF
                SIN_HP=SINH_HP0
                I=I+2
             ENDDO
             IF(I==NMAX_pol) THEN
                line="NO CONVERGENCE IN SINH_HP"
                ipause=mypauses(NMAX_pol,line)
             ENDIF
       CALL kill(SIN_HP); CALL kill(Y); CALL kill(SINH_HP0);

          else
             SIN_HP=SIN(varf1)/varf1
          endif
          sinX_Xt%t= SIN_HP
          master=localmaster
       else
          sinX_Xt%r= SIN_HR(S1%r)
          sinX_Xt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in sinX_XT "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION sinX_XT

  ! cosine for quaternion DONE partly COSY-INFINITY WISE FOR TPSA
  function sin_quaternionr(X)
    implicit none
    real(dp), INTENT (IN) :: X
    real(dp) sin_quaternionr,Y,NORM0,NORM,SINH_HR0
    logical(lp) NOTDONE,CHECK
    INTEGER I,ipause,mypauses

    if(abs(x)<.1d0) then
       NOTDONE=.TRUE.
       CHECK=.TRUE.
       sin_quaternionr=1.0_dp
       Y=1.0_dp
       I=2
       NORM0=1e5_dp
       DO WHILE(I<NMAX_pol.AND.NOTDONE)
          Y=-Y*X/REAL(I,kind=DP)/REAL(I+1,kind=DP)
          SINH_HR0=sin_quaternionr+Y
          NORM=ABS(sin_quaternionr-SINH_HR0)
          IF(NORM<=EPS_real_poly.AND.CHECK) THEN
             NORM0=NORM
             CHECK=.FALSE.
          ELSE
             IF(NORM>=NORM0) THEN
                NOTDONE=.FALSE.
             ELSE
                NORM0=NORM
             ENDIF
          ENDIF
          sin_quaternionr=SINH_HR0
          I=I+2
       ENDDO
       IF(I==NMAX_pol) THEN
          line="NO CONVERGENCE IN SINH_HR"
          ipause=mypauses(NMAX_pol,line)
       ENDIF
    else
       sin_quaternionr=sin(sqrt(x))/sqrt(x)
    endif
    return
  end function sin_quaternionr

FUNCTION sin_quaternionp( S1 )
    implicit none
    integer ipause, mypauses
    TYPE (real_8) sin_quaternionp
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster
    TYPE(TAYLOR)SIN_HP
    TYPE(TAYLOR)Y,SINH_HP0
    real(dp) NORM0,NORM
    logical(lp) NOTDONE,CHECK
    INTEGER I

    select case(s1%kind)
    case(m1)
       sin_quaternionp%r=sin_quaternionr(s1%r)
       sin_quaternionp%kind=1
    case(m2)
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
       localmaster=master
       call ass(sin_quaternionp)

       if(ABS(S1%T.SUB.'0')<.1d0) then
          NOTDONE=.TRUE.
          CHECK=.TRUE.

          SIN_HP=1.0_dp
          Y=1.0_dp
          I=2
          NORM0=1e5_dp
          DO WHILE(I<NMAX_pol.AND.NOTDONE)
             Y=-Y*S1%T/REAL(I,kind=DP)/REAL(I+1,kind=DP)
             SINH_HP0=SIN_HP+Y
             NORM=full_abs(SIN_HP-SINH_HP0)
             IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                NORM0=NORM
                CHECK=.FALSE.
             ELSE
                IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                   NOTDONE=.FALSE.
                ELSE
                   NORM0=NORM
                ENDIF
             ENDIF
             SIN_HP=SINH_HP0
             I=I+2
          ENDDO
          IF(I==NMAX_pol) THEN
             line="NO CONVERGENCE IN SIN_HP"
             ipause=mypauses(NMAX_pol,line)
          ENDIF
       else
          SIN_HP=sin(sqrt(S1%T)) / sqrt(S1%T)
       endif
       sin_quaternionp%t= SIN_HP

       master=localmaster
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);
    case(m3)
       if(knob) then
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
          localmaster=master
          call ass(sin_quaternionp)
          call varfk1(S1)
          if(ABS(varf1.SUB.'0')<.1d0) then
             NOTDONE=.TRUE.
             CHECK=.TRUE.

             SIN_HP=1.0_dp
             Y=1.0_dp
             I=2
             NORM0=1e5_dp
             DO WHILE(I<NMAX_pol.AND.NOTDONE)
                Y=-Y*varf1/REAL(I,kind=DP)/REAL(I+1,kind=DP)
                SINH_HP0=SIN_HP+Y
                NORM=full_abs(SIN_HP-SINH_HP0)
                IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                   NORM0=NORM
                   CHECK=.FALSE.
                ELSE
                   IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                      NOTDONE=.FALSE.
                   ELSE
                      NORM0=NORM
                   ENDIF
                ENDIF
                SIN_HP=SINH_HP0
                I=I+2
             ENDDO
             IF(I==NMAX_pol) THEN
                line="NO CONVERGENCE IN SINH_HP"
                ipause=mypauses(NMAX_pol,line)
             ENDIF
          else
             SIN_HP=sin(sqrt(varf1))/sqrt(varf1)
          endif
          sin_quaternionp%t= SIN_HP
          master=localmaster
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);
       else
          sin_quaternionp%r= sin_quaternionr(S1%r)
          sin_quaternionp%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in sinX_XT "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION sin_quaternionp

  ! cosine for quaternion DONE partly COSY-INFINITY WISE FOR TPSA
  function cos_quaternionr(X)
    implicit none
    real(dp), INTENT (IN) :: X
    real(dp) cos_quaternionr,Y,NORM0,NORM,SINH_HR0
    logical(lp) NOTDONE,CHECK
    INTEGER I,ipause,mypauses

    if(abs(x)<.1d0) then
       NOTDONE=.TRUE.
       CHECK=.TRUE.
       cos_quaternionr=1.0_dp
       Y=1.0_dp
       I=1
       NORM0=1e5_dp
       DO WHILE(I<NMAX_pol.AND.NOTDONE)
          Y=-Y*X/REAL(I,kind=DP)/REAL(I+1,kind=DP)
          SINH_HR0=cos_quaternionr+Y
          NORM=ABS(cos_quaternionr-SINH_HR0)
          IF(NORM<=EPS_real_poly.AND.CHECK) THEN
             NORM0=NORM
             CHECK=.FALSE.
          ELSE
             IF(NORM>=NORM0) THEN
                NOTDONE=.FALSE.
             ELSE
                NORM0=NORM
             ENDIF
          ENDIF
          cos_quaternionr=SINH_HR0
          I=I+2
       ENDDO
       IF(I==NMAX_pol) THEN
          line="NO CONVERGENCE IN SINH_HR"
          ipause=mypauses(NMAX_pol,line)
       ENDIF
    else
       cos_quaternionr=cos(sqrt(x))
    endif
    return
  end function cos_quaternionr

  FUNCTION cos_quaternionp( S1 )
    implicit none
    integer ipause, mypauses
    TYPE (real_8) cos_quaternionp
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster
    TYPE(TAYLOR)SIN_HP
    TYPE(TAYLOR)Y,SINH_HP0
    real(dp) NORM0,NORM
    logical(lp) NOTDONE,CHECK
    INTEGER I

    select case(s1%kind)
    case(m1)
       cos_quaternionp%r=cos_quaternionr(s1%r)
       cos_quaternionp%kind=1
    case(m2)
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
       localmaster=master
       call ass(cos_quaternionp)

       if(ABS(S1%T.SUB.'0')<.1d0) then
          NOTDONE=.TRUE.
          CHECK=.TRUE.

          SIN_HP=1.0_dp
          Y=1.0_dp
          I=1
          NORM0=1e5_dp
          DO WHILE(I<NMAX_pol.AND.NOTDONE)
             Y=-Y*S1%T/REAL(I,kind=DP)/REAL(I+1,kind=DP)
             SINH_HP0=SIN_HP+Y
             NORM=full_abs(SIN_HP-SINH_HP0)
             IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                NORM0=NORM
                CHECK=.FALSE.
             ELSE
                IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                   NOTDONE=.FALSE.
                ELSE
                   NORM0=NORM
                ENDIF
             ENDIF
             SIN_HP=SINH_HP0
             I=I+2
          ENDDO
          IF(I==NMAX_pol) THEN
             line="NO CONVERGENCE IN SIN_HP"
             ipause=mypauses(NMAX_pol,line)
          ENDIF
       else
          SIN_HP=cos(sqrt(S1%T)) 
       endif
       cos_quaternionp%t= SIN_HP

       master=localmaster
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);
    case(m3)
       if(knob) then
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
          localmaster=master
          call ass(cos_quaternionp)
          call varfk1(S1)
          if(ABS(varf1.SUB.'0')<.1d0) then
             NOTDONE=.TRUE.
             CHECK=.TRUE.

             SIN_HP=1.0_dp
             Y=1.0_dp
             I=1
             NORM0=1e5_dp
             DO WHILE(I<NMAX_pol.AND.NOTDONE)
                Y=-Y*varf1/REAL(I,kind=DP)/REAL(I+1,kind=DP)
                SINH_HP0=SIN_HP+Y
                NORM=full_abs(SIN_HP-SINH_HP0)
                IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                   NORM0=NORM
                   CHECK=.FALSE.
                ELSE
                   IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                      NOTDONE=.FALSE.
                   ELSE
                      NORM0=NORM
                   ENDIF
                ENDIF
                SIN_HP=SINH_HP0
                I=I+2
             ENDDO
             IF(I==NMAX_pol) THEN
                line="NO CONVERGENCE IN SINH_HP"
                ipause=mypauses(NMAX_pol,line)
             ENDIF
          else
             SIN_HP=cos(sqrt(varf1))
          endif
          cos_quaternionp%t= SIN_HP
          master=localmaster
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);
       else
          cos_quaternionp%r= cos_quaternionr(S1%r)
          cos_quaternionp%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in sinX_XT "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION cos_quaternionp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION sinHX_Xt( S1 )
    implicit none
    integer ipause, mypauses
    TYPE (real_8) sinHX_Xt
    TYPE (real_8), INTENT (IN) :: S1
    integer localmaster
    TYPE(TAYLOR)SIN_HP
    TYPE(TAYLOR)Y,SINH_HP0
    real(dp) NORM0,NORM
    logical(lp) NOTDONE,CHECK
    INTEGER I

    select case(s1%kind)
    case(m1)
       sinHX_Xt%r=SINH_HR(s1%r)
       sinHX_Xt%kind=1
    case(m2)
       CALL ALLOC(SIN_HP); CALL ALLOC(Y); CALL ALLOC(SINH_HP0);
       localmaster=master
       call ass(sinHX_Xt)

       if(ABS(S1%T.SUB.'0')<sinhx_x_minp) then
          NOTDONE=.TRUE.
          CHECK=.TRUE.

          SIN_HP=1.0_dp
          Y=1.0_dp
          I=1
          NORM0=1e5_dp
          DO WHILE(I<NMAX_pol.AND.NOTDONE)
             Y=Y*S1%T**2/REAL(I+1,kind=DP)/REAL(I+2,kind=DP)
             SINH_HP0=SIN_HP+Y
             NORM=full_abs(SIN_HP-SINH_HP0)
             IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                NORM0=NORM
                CHECK=.FALSE.
             ELSE
                IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                   NOTDONE=.FALSE.
                ELSE
                   NORM0=NORM
                ENDIF
             ENDIF
             SIN_HP=SINH_HP0
             I=I+2
          ENDDO
          IF(I==NMAX_pol) THEN
             line="NO CONVERGENCE IN SINH_HP"
             ipause=mypauses(NMAX_pol,line)
          ENDIF
       else

          SIN_HP=SINH(S1%T)/S1%T
       endif
       sinHX_Xt%t= SIN_HP

       master=localmaster
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);
    case(m3)
       if(knob) then
       CALL alloc(SIN_HP); CALL alloc(Y); CALL alloc(SINH_HP0);

          localmaster=master
          call ass(sinHX_Xt)
          call varfk1(S1)
          if(ABS(varf1.SUB.'0')<sinhx_x_minp) then
             NOTDONE=.TRUE.
             CHECK=.TRUE.

             SIN_HP=1.0_dp
             Y=1.0_dp
             I=1
             NORM0=1e5_dp
             DO WHILE(I<NMAX_pol.AND.NOTDONE)
                Y=Y*varf1**2/REAL(I+1,kind=DP)/REAL(I+2,kind=DP)
                SINH_HP0=SIN_HP+Y
                NORM=full_abs(SIN_HP-SINH_HP0)
                IF(NORM<=EPS_real_poly.AND.CHECK) THEN
                   NORM0=NORM
                   CHECK=.FALSE.
                ELSE
                   IF(NORM>=NORM0.and.(.not.check)) THEN  ! sagan
                      NOTDONE=.FALSE.
                   ELSE
                      NORM0=NORM
                   ENDIF
                ENDIF
                SIN_HP=SINH_HP0
                I=I+2
             ENDDO
             IF(I==NMAX_pol) THEN
                line="NO CONVERGENCE IN SINH_HP"
                ipause=mypauses(NMAX_pol,line)
             ENDIF
          else
             SIN_HP=SINH(varf1)/varf1
          endif
          sinHX_Xt%t= SIN_HP
          master=localmaster
       else
          sinHX_Xt%r= SINH_HR(S1%r)
          sinHX_Xt%kind=1
       endif
       CALL KILL(SIN_HP); CALL KILL(Y); CALL KILL(SINH_HP0);

    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in sinHX_Xt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION sinHX_Xt

  ! remove small numbers

  SUBROUTINE  clean_real_8(S1,S2,prec)
    implicit none
    type (real_8),INTENT(INOUT)::S2
    type (real_8), intent(INOUT):: s1
    real(dp) prec
    type(real_8) t

    call alloc(t)
    t=s1

    select case(s1%kind)
    case(m1)
       if(abs(t%r)<prec) t%r=0.0_dp
    case(m2)
       call clean_taylor(t%t,t%t,prec)
    case(m3)
       Write(6,*) " cannot clean a knob "
       stop 601
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in clean_real_8 "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
    s2=t
    call kill(t)


  END SUBROUTINE clean_real_8

  SUBROUTINE  flip_real_8(S1,S2,i)
    implicit none
    type (real_8),INTENT(INOUT)::S2
    type (real_8), intent(INOUT):: s1
    integer i
    if(s1%kind==2) then
       s2=s1  ! to make s2 taylor!
       call flip_taylor(S1%t,S2%t,i)
    endif
  end SUBROUTINE  flip_real_8


 

end module  polymorphic_taylor
