!The Full Polymorphic Package
!Copyright (C) Etienne Forest
! Taylor polymorphism at execution is based on an idea
! and C++ prototype developed  by J. Bengtsson circa 1990

module polymorphic_complextaylor
  use polymorphic_taylor
  implicit none
  public
  integer,private,parameter::m1=1,m2=2,m3=3,ms=4
  integer,private,parameter:: m11=m1+ms*m1,m12=m1+ms*m2,m13=m1+ms*m3,  &
       m21=m2+ms*m1,m22=m2+ms*m2,m23=m2+ms*m3,                         &
       m31=m3+ms*m1,m32=m3+ms*m2,m33=m3+ms*m3
  logical(lp),private,parameter::t=.true.,f=.false.
  integer,target,private::NO,ND,ND2,NP,NDPT,NV
  logical(lp),private,TARGET::old
  private set_in_poly,init_map_cp,init_tpsa_cp
  private equal ,Dequaldacon,cequaldacon,equaldacon,iequaldacon
  private EQUALRP,complexEQUAL,EQUALcomplext,complextEQUAL
  private resetpoly,resetpolyn,resetpoly0,resetpolyn0,k_opt
  private printpoly,allocpoly,allocpolyn,resetpoly_R,resetpoly_Rn,A_OPT
  private add,unaryADD,daddsc,caddsc,dscadd,cscadd,addsc,scadd,iaddsc,iscadd
  private subs,unarySUB,cscsub,csubsc,dscsub,dsubsc,subsc,scsub,isubsc,iscsub
  private mul,dmulsc,dscmul,cmulsc,cscmul,mulsc,scmul,imulsc,iscmul
  private pmul,mulp,padd,addp,psub,subp,pdiv,divp
  private cpmulsc,cpscmul,cpaddsc,cpscadd,cpsubsc,cpscsub ,cpdivsc,cpscdiv
  private div,cdivsc,cscdiv,ddivsc,dscdiv,divsc,scdiv,idivsc,iscdiv
  private POW,POWR,POWR8
  private dexpt,dcost,dsint,dlogt,dsqrtt,abst,CONJGT,print6
  PRIVATE dimagt,drealt,dcmplxt,GETint,GETORDER,CUTORDER,getchar,GETCHARnd2,GETintnd2
  !
  private asscp,make_it_knobc,kill_knobc
  private line
  character(120) line

  ! Starts the Taylor series Package
  INTERFACE init
     MODULE PROCEDURE init_tpsa_cp    !@1 &nbsp; Initializes a pure TPSA package
     MODULE PROCEDURE init_map_cp     !@1 &nbsp; Initializes a Lie/Differential Algebraic Package
  END INTERFACE


  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL
     MODULE PROCEDURE complexEQUAL
     MODULE PROCEDURE EQUALcomplext
     MODULE PROCEDURE complextEQUAL
     MODULE PROCEDURE EQUALRP  !2002.10.9
     MODULE PROCEDURE RPEQUAL !2002.10.9
     MODULE PROCEDURE Dequaldacon
     MODULE PROCEDURE equaldacon
     MODULE PROCEDURE iequaldacon
     MODULE PROCEDURE cequaldacon
  end  INTERFACE

  ! Operators
  !@     <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="684" height="445">
  !@       <tr>
  !@         <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; +</b>&nbsp;&nbsp;&nbsp; </font></span></td>
  !@         <td width="176" height="34" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">&nbsp;COMPLEX</font></span></td>
  !@         <td width="306" height="42" align="center" colspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="94" height="84" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="84" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@         <td width="92" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="105" height="42" align="center">
  !@         <font size="2">REAL_8</font></td>
  !@         <td width="102" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@         <td width="99" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="108" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX</font></span></td>
  !@         <td width="39" height="53" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@         <td width="84" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#ADD">ADD</a></font></td>
  !@         <td width="92" height="53" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CADDSC">CADDSC</a></span></font></td>
  !@         <td width="105" height="53" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#ADDP">ADDP</a></span></font></td>
  !@         <td width="102" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DADDSC">DADDSC</a></font></td>
  !@         <td width="99" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#ADDSC">ADDSC</a></font></td>
  !@         <td width="94" height="53" align="center"><font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#IADDSC">IADDSC</a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="50" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="84" height="50" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CSCADD">CSCADD</a></span></font></td>
  !@         <td width="92" height="50" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="50" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CPSCADD">CPSCADD</a></span></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@            <td width="78" height="20" align="center">
  !@            <b>F90</b></td>
  !@            <td width="56" height="20" align="center">
  !@            <b>&nbsp;&nbsp;&nbsp; F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="167" align="center" rowspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL_8</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#PADD">PADD</a></span></font></td>
  !@         <td width="92" height="55" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CPADDSC">CPADDSC</a></span></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <font color="#FF0000" size="2"><i>
  !@         <a href="m_real_polymorph.htm#ADD" style="text-decoration: none">
  !@         <font color="#FF0000">ADD</font></a></i></font></td>
  !@         <td width="102" height="55" align="center">
  !@         <i><font size="2" color="#FF0000">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#DADDSC">
  !@         <font color="#FF0000">DADDSC</font></a></font></i></td>
  !@         <td width="99" height="55" align="center">
  !@         <i><font size="2" color="#FF0000">
  !@         <a href="m_real_polymorph.htm#ADDSC" style="text-decoration: none">
  !@         <font color="#FF0000">ADDSC</font></a></font></i></td>
  !@         <td width="94" height="55" align="center"><i>
  !@         <font size="2" color="#FF0000">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#IADDSC">
  !@         <font color="#FF0000">IADDSC</font></a></font></i></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DSCADD">DSCADD</a></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="55" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a href="m_real_polymorph.htm#DSCADD" style="text-decoration: none">
  !@         <font color="#FF0000">DSCADD</font></a></i></font></td>
  !@         <td width="102" height="55" align="center"><b>F90</b></td>
  !@         <td width="99" height="55" align="center"><b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="52" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@         <td width="84" height="52" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#SCADD">SCADD</a></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="52" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a href="m_real_polymorph.htm#SCADD" style="text-decoration: none">
  !@         <font color="#FF0000">SCADD</font></a></i></font></td>
  !@         <td width="102" height="52" align="center"><b>F90</b></td>
  !@         <td width="99" height="52" align="center"><b>F90</b></td>
  !@         <td width="94" height="52" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="78" height="61" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@         <td width="84" height="61" align="center"><font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#ISCADD">ISCADD</a></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="61" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a href="m_real_polymorph.htm#ISCADD" style="text-decoration: none">
  !@         <font color="#FF0000">ISCADD</font></a></i></font></td>
  !@         <td width="102" height="61" align="center"><b>F90</b></td>
  !@         <td width="99" height="61" align="center"><b>F90</b></td>
  !@         <td width="94" height="61" align="center"><b>F90</b></td>
  !@       </tr>
  !@     </table>



  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add
     MODULE PROCEDURE padd
     MODULE PROCEDURE addp
     MODULE PROCEDURE cpscadd
     MODULE PROCEDURE cpaddsc
     MODULE PROCEDURE unaryADD
     MODULE PROCEDURE daddsc
     MODULE PROCEDURE caddsc
     MODULE PROCEDURE dscadd
     MODULE PROCEDURE cscadd
     MODULE PROCEDURE addsc
     MODULE PROCEDURE scadd
     MODULE PROCEDURE iaddsc
     MODULE PROCEDURE iscadd
  end  INTERFACE

  !@     <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber4" width="684" height="445">
  !@       <tr>
  !@         <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </b>
  !@         -&nbsp;&nbsp; </font></span></td>
  !@         <td width="176" height="34" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">&nbsp;COMPLEX</font></span></td>
  !@         <td width="306" height="42" align="center" colspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="94" height="84" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="84" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@         <td width="92" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="105" height="42" align="center">
  !@         <font size="2">REAL_8</font></td>
  !@         <td width="102" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@         <td width="99" height="45" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="108" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX</font></span></td>
  !@         <td width="39" height="53" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@         <td width="84" height="53" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#SUBS">SUSB</a></span></font></td>
  !@         <td width="92" height="53" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CSUBSC">CSUBSC</a></span></font></td>
  !@         <td width="105" height="53" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#SUBP">SUBP</a></span></font></td>
  !@         <td width="102" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DSUBSC">DSUBSC</a></font></td>
  !@         <td width="99" height="53" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#SUBSC">
  !@         SUBSC</a></span></font></td>
  !@         <td width="94" height="53" align="center"><font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#ISUBSC">ISUBSC</a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="50" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="84" height="50" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CSCSUB">
  !@         CSCSUB</a></span></font></td>
  !@         <td width="92" height="50" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="50" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CPSCSUB">
  !@         CPSCSUB</a></span></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@            <td width="78" height="20" align="center">
  !@            <b>F90</b></td>
  !@            <td width="56" height="20" align="center">
  !@            <b>&nbsp;&nbsp;&nbsp; F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="167" align="center" rowspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL_8</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#PSUB">PSUB</a></span></font></td>
  !@         <td width="92" height="55" align="center">
  !@         <font size="2"><span style="font-weight: 700">
  !@         <a style="text-decoration: none" href="n_complex_polymorph.htm#CPSUBSC">CPSUBSC</a></span></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#SUBS">
  !@         <font color="#FF0000">SUBS</font></a></i></font></td>
  !@         <td width="102" height="55" align="center">
  !@         <i><font size="2" color="#FF0000">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#DSUBSC">
  !@         <font color="#FF0000">DSUBSC</font></a></font></i></td>
  !@         <td width="99" height="55" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#SUBSC">
  !@         <font color="#FF0000">SUBSC</font></a></i></font></td>
  !@         <td width="94" height="55" align="center"><i>
  !@         <font size="2" color="#FF0000">
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#ISUBSC">
  !@         <font color="#FF0000">ISUBSC</font></a></font></i></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DSCSUB">
  !@         DSCSUB</a></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="55" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#DSCSUB">
  !@         <font color="#FF0000">DSCSUB</font></a></i></font></td>
  !@         <td width="102" height="55" align="center"><b>F90</b></td>
  !@         <td width="99" height="55" align="center"><b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="52" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@         <td width="84" height="52" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#SCSUB">
  !@         SCSUB</a></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="52" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#SCSUB">
  !@         <font color="#FF0000">SCSUB</font></a></i></font></td>
  !@         <td width="102" height="52" align="center"><b>F90</b></td>
  !@         <td width="99" height="52" align="center"><b>F90</b></td>
  !@         <td width="94" height="52" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="78" height="61" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@         <td width="84" height="61" align="center"><font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#ISCSUB">
  !@         ISCSUB</a></font></td>
  !@            <td width="77" height="20" align="center">
  !@            <b>F90</b></td>
  !@         <td width="105" height="61" align="center">
  !@         <font size="2" color="#FF0000"><i>
  !@         <a style="text-decoration: none" href="m_real_polymorph.htm#ISCSUB">
  !@         <font color="#FF0000">ISCSUB</font></a></i></font></td>
  !@         <td width="102" height="61" align="center"><b>F90</b></td>
  !@         <td width="99" height="61" align="center"><b>F90</b></td>
  !@         <td width="94" height="61" align="center"><b>F90</b></td>
  !@       </tr>
  !@     </table>

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE unarySUB
     MODULE PROCEDURE subs
     MODULE PROCEDURE cpscsub
     MODULE PROCEDURE cpsubsc
     MODULE PROCEDURE psub
     MODULE PROCEDURE subp
     MODULE PROCEDURE cscsub
     MODULE PROCEDURE csubsc
     MODULE PROCEDURE dscsub
     MODULE PROCEDURE dsubsc
     MODULE PROCEDURE subsc
     MODULE PROCEDURE scsub
     MODULE PROCEDURE isubsc
     MODULE PROCEDURE iscsub
  END INTERFACE

  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber5" width="684" height="445">
  !@      <tr>
  !@        <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </b>
  !@        *&nbsp;&nbsp;&nbsp; </font></span></td>
  !@        <td width="176" height="34" align="center" colspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">&nbsp;COMPLEX</font></span></td>
  !@        <td width="306" height="42" align="center" colspan="3">
  !@        <font size="2">REAL</font></td>
  !@        <td width="94" height="84" align="center" rowspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">Integer</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="84" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@        <td width="92" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@        <td width="105" height="42" align="center">
  !@        <font size="2">REAL_8</font></td>
  !@        <td width="102" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@        <td width="99" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="108" align="center" rowspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">COMPLEX</font></span></td>
  !@        <td width="39" height="53" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@        <td width="84" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#MUL">MUL</a></span></font></td>
  !@        <td width="92" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CMULSC">CMULSC</a></span></font></td>
  !@        <td width="105" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#MULP">MULP</a></span></font></td>
  !@        <td width="102" height="53" align="center">
  !@        <font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DMULSC">DMULSC</a></font></td>
  !@        <td width="99" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#MULSC">
  !@        MULSC</a></span></font></td>
  !@        <td width="94" height="53" align="center"><font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#IMULSC">IMULSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="50" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@        <td width="84" height="50" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CSCMUL">
  !@        CSCMUL</a></span></font></td>
  !@        <td width="92" height="50" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="50" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CPSCMUL">
  !@        CPSCMUL</a></span></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@           <td width="78" height="20" align="center">
  !@           <b>F90</b></td>
  !@           <td width="56" height="20" align="center">
  !@           <b>&nbsp;&nbsp;&nbsp; F90</b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="167" align="center" rowspan="3">
  !@        <font size="2">REAL</font></td>
  !@        <td width="39" height="55" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL_8</font></span></td>
  !@        <td width="84" height="55" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#PMUL">PMUL</a></span></font></td>
  !@        <td width="92" height="55" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CPMULSC">CPMULSC</a></span></font></td>
  !@        <td width="105" height="55" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#MUL">
  !@        <font color="#FF0000">MUL</font></a></i></font></td>
  !@        <td width="102" height="55" align="center">
  !@        <i><font size="2" color="#FF0000">
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#DMULSC">
  !@        <font color="#FF0000">DMULSC</font></a></font></i></td>
  !@        <td width="99" height="55" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#MULSC">
  !@        <font color="#FF0000">MULSC</font></a></i></font></td>
  !@        <td width="94" height="55" align="center"><i>
  !@        <font size="2" color="#FF0000">
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#IMULSC">
  !@        <font color="#FF0000">IMULSC</font></a></font></i></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="55" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@        <td width="84" height="55" align="center">
  !@        <font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DSCMUL">
  !@        DSCMUL</a></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="55" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#DSCMUL">
  !@        <font color="#FF0000">DSCMUL</font></a></i></font></td>
  !@        <td width="102" height="55" align="center"><b>F90</b></td>
  !@        <td width="99" height="55" align="center"><b>F90</b></td>
  !@        <td width="94" height="55" align="center"><b>F90</b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="52" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@        <td width="84" height="52" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#SCMUL">
  !@        SCMUL</a></span></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="52" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#SCMUL">
  !@        <font color="#FF0000">SCMUL</font></a></i></font></td>
  !@        <td width="102" height="52" align="center"><b>F90</b></td>
  !@        <td width="99" height="52" align="center"><b>F90</b></td>
  !@        <td width="94" height="52" align="center"><b>F90</b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="78" height="61" align="center" colspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">Integer</font></span></td>
  !@        <td width="84" height="61" align="center"><font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#ISCMUL">
  !@        ISCMUL</a></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="61" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#ISCMUL">
  !@        <font color="#FF0000">ISCMUL</font></a></i></font></td>
  !@        <td width="102" height="61" align="center"><b>F90</b></td>
  !@        <td width="99" height="61" align="center"><b>F90</b></td>
  !@        <td width="94" height="61" align="center"><b>F90</b></td>
  !@      </tr>
  !@    </table>


  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul
     MODULE PROCEDURE pmul
     MODULE PROCEDURE mulp
     MODULE PROCEDURE cmulsc
     MODULE PROCEDURE cscmul
     MODULE PROCEDURE cpmulsc
     MODULE PROCEDURE cpscmul
     MODULE PROCEDURE dmulsc
     MODULE PROCEDURE dscmul
     MODULE PROCEDURE mulsc
     MODULE PROCEDURE scmul
     MODULE PROCEDURE imulsc
     MODULE PROCEDURE iscmul
  END INTERFACE

  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber6" width="684" height="445">
  !@      <tr>
  !@        <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  !@        /</b> &nbsp;&nbsp;&nbsp; </font></span></td>
  !@        <td width="176" height="34" align="center" colspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">&nbsp;COMPLEX</font></span></td>
  !@        <td width="306" height="42" align="center" colspan="3">
  !@        <font size="2">REAL</font></td>
  !@        <td width="94" height="84" align="center" rowspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">Integer</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="84" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@        <td width="92" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@        <td width="105" height="42" align="center">
  !@        <font size="2">REAL_8</font></td>
  !@        <td width="102" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@        <td width="99" height="45" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="108" align="center" rowspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">COMPLEX</font></span></td>
  !@        <td width="39" height="53" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">DOUBLE COMPLEX</font></span></td>
  !@        <td width="84" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#DIV">DIV</a></span></font></td>
  !@        <td width="92" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CDIVSC">CDIVSC</a></span></font></td>
  !@        <td width="105" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#DIVP">DIVP</a></span></font></td>
  !@        <td width="102" height="53" align="center">
  !@        <font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DDIVSC">DDIVSC</a></font></td>
  !@        <td width="99" height="53" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#DIVSC">
  !@        DIVSC</a></span></font></td>
  !@        <td width="94" height="53" align="center"><font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#IDIVSC">IDIVSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="50" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@        <td width="84" height="50" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CSCDIV">
  !@        CSCDIV</a></span></font></td>
  !@        <td width="92" height="50" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="50" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CPSCDIV">
  !@        CPSCDIV</a></span></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@           <td width="78" height="20" align="center">
  !@           <b>F90</b></td>
  !@           <td width="56" height="20" align="center">
  !@           <b>&nbsp;&nbsp;&nbsp; F90</b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="167" align="center" rowspan="3">
  !@        <font size="2">REAL</font></td>
  !@        <td width="39" height="55" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL_8</font></span></td>
  !@        <td width="84" height="55" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#PDIV">PDIV</a></span></font></td>
  !@        <td width="92" height="55" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#CPDIVSC">CPDIVSC</a></span></font></td>
  !@        <td width="105" height="55" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#DIV">
  !@        <font color="#FF0000">DIV</font></a></i></font></td>
  !@        <td width="102" height="55" align="center">
  !@        <i><font size="2" color="#FF0000">
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#DDIVSC">
  !@        <font color="#FF0000">DDIVSC</font></a></font></i></td>
  !@        <td width="99" height="55" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#DIVSC">
  !@        <font color="#FF0000">DIVSC</font></a></i></font></td>
  !@        <td width="94" height="55" align="center"><i>
  !@        <font size="2" color="#FF0000">
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#IDIVSC">
  !@        <font color="#FF0000">IDIVSC</font></a></font></i></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="55" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@        <td width="84" height="55" align="center">
  !@        <font size="2">
  !@        <a style="text-decoration: none; font-weight: 700" href="n_complex_polymorph.htm#DSCDIV">
  !@        DSCDIV</a></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="55" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#DSCDIV">
  !@        <font color="#FF0000">DSCDIV</font></a></i></font></td>
  !@        <td width="102" height="55" align="center"><b>F90</b></td>
  !@        <td width="99" height="55" align="center"><b>F90</b></td>
  !@        <td width="94" height="55" align="center"><b>F90</b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="39" height="52" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@        <td width="84" height="52" align="center">
  !@        <font size="2"><span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#SCDIV">
  !@        SCDIV</a></span></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="52" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#SCDIV">
  !@        <font color="#FF0000">SCDIV</font></a></i></font></td>
  !@        <td width="102" height="52" align="center"><b>F90</b></td>
  !@        <td width="99" height="52" align="center"><b>F90</b></td>
  !@        <td width="94" height="52" align="center"><b>F90</b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="78" height="61" align="center" colspan="2">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="2">Integer</font></span></td>
  !@        <td width="84" height="61" align="center"><font size="2">
  !@        <span style="font-weight: 700">
  !@        <a style="text-decoration: none" href="n_complex_polymorph.htm#ISCDIV">
  !@        ISCDIV</a></span></font></td>
  !@           <td width="77" height="20" align="center">
  !@           <b>F90</b></td>
  !@        <td width="105" height="61" align="center">
  !@        <font size="2" color="#FF0000"><i>
  !@        <a style="text-decoration: none" href="m_real_polymorph.htm#ISCDIV">
  !@        <font color="#FF0000">ISCDIV</font></a></i></font></td>
  !@        <td width="102" height="61" align="center"><b>F90</b></td>
  !@        <td width="99" height="61" align="center"><b>F90</b></td>
  !@        <td width="94" height="61" align="center"><b>F90</b></td>
  !@      </tr>
  !@    </table>

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div
     MODULE PROCEDURE cpdivsc
     MODULE PROCEDURE cpscdiv
     MODULE PROCEDURE pdiv
     MODULE PROCEDURE divp
     MODULE PROCEDURE cdivsc
     MODULE PROCEDURE cscdiv
     MODULE PROCEDURE ddivsc
     MODULE PROCEDURE dscdiv
     MODULE PROCEDURE divsc
     MODULE PROCEDURE scdiv
     MODULE PROCEDURE idivsc
     MODULE PROCEDURE iscdiv
  END INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW
     MODULE PROCEDURE POWR
     MODULE PROCEDURE POWR8
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




  ! Intrinsic Functions

  INTERFACE aimag
     MODULE PROCEDURE dimagt
  END INTERFACE
  INTERFACE dimag
     MODULE PROCEDURE dimagt
  END INTERFACE


  INTERFACE real
     MODULE PROCEDURE drealt
  END INTERFACE
  INTERFACE dreal
     MODULE PROCEDURE drealt
  END INTERFACE
  INTERFACE dble
     MODULE PROCEDURE drealt
  END INTERFACE


  INTERFACE cmplx
     MODULE PROCEDURE dcmplxt
  END INTERFACE
  INTERFACE dcmplx
     MODULE PROCEDURE dcmplxt
  END INTERFACE


  INTERFACE exp
     MODULE PROCEDURE dexpt
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
  INTERFACE cdsqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE

  INTERFACE abs
     MODULE PROCEDURE abst
  END INTERFACE
  INTERFACE dabs
     MODULE PROCEDURE abst
  END INTERFACE

  ! Non Intrinsic Functions
  INTERFACE morph
     MODULE PROCEDURE polymorpht
  END INTERFACE

  INTERFACE CONJG
     MODULE PROCEDURE CONJGT
  END INTERFACE

  ! i/o

  INTERFACE daprint
     MODULE PROCEDURE printpoly
     MODULE PROCEDURE print6
  END INTERFACE
  INTERFACE print
     MODULE PROCEDURE print6
     MODULE PROCEDURE printpoly
  END INTERFACE


  ! Constructors and Destructors
  INTERFACE alloc
     MODULE PROCEDURE A_OPT !allocpoly
     MODULE PROCEDURE allocpolyn
  END INTERFACE

  INTERFACE kill
     MODULE PROCEDURE K_OPT  !resetpoly0
     MODULE PROCEDURE resetpolyn0
     MODULE PROCEDURE resetpoly_R  !
     MODULE PROCEDURE resetpoly_Rn
  END INTERFACE

  INTERFACE reset
     MODULE PROCEDURE resetpoly
     MODULE PROCEDURE resetpolyn
  END INTERFACE


interface make_it_knob
module procedure make_it_knobc
end INTERFACE

interface kill_knob
module procedure kill_knobc
end INTERFACE


  ! end Constructors and Destructors


  ! managing

  INTERFACE ass
     MODULE PROCEDURE asscp
  END INTERFACE


contains


  FUNCTION polymorpht( S1 )
    implicit none
    TYPE (complex_8) polymorpht
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster

    localmaster=master
    call ass(polymorpht)
    polymorpht%t= s1
    master=localmaster
  END FUNCTION polymorpht


  subroutine set_da_pointers()
    use da_arrays
    implicit none
    c_%total_da_size => total_da_size
    c_%lda_used => lda_used
    c_%real_warning => real_warning
    c_%check_da => check_da
    c_%stable_da => stable_da
    ! c_%stable_da => CHECK_STABLE
    c_%no => no
    c_%nv => nv
    c_%nd => nd
    c_%nd2 => nd2
    c_%np => np
    c_%NPARA => NPARA_original
    c_%NSPIN => NSPIN
    c_%SPIN_POS => SPIN_POS
    c_%ndpt => ndpt
    c_%npara_fpp => npara_fpp
    c_%knob => knob
    c_%setknob => setknob
    c_%OLD => OLD
    c_%global_verbose => global_verbose
    c_%da_absolute_aperture=>da_absolute_aperture

    c_%ROOT_CHECK => ROOT_CHECK
    c_%CHECK_STABLE => CHECK_STABLE
    c_%CHECK_MADX_APERTURE => CHECK_MADX_APERTURE
    c_%ROOT_CHECK => ROOT_CHECK
    c_%APERTURE_FLAG => APERTURE_FLAG
    c_%s_aperture_CHECK => s_aperture_CHECK
    c_%absolute_aperture => absolute_aperture

    c_%hyperbolic_aperture => hyperbolic_aperture
    c_%WATCH_USER => WATCH_USER
    c_%no_hyperbolic_in_normal_form => no_hyperbolic_in_normal_form

  end subroutine set_da_pointers

  subroutine init_map_cp(NO1,ND1,NP1,NDPT1,PACKAGE)
    implicit none
    integer NO1,ND1,NP1
    integer,optional ::  NDPT1
    logical(lp),optional :: PACKAGE
    logical(lp) PACKAGE1
    integer ndptt,i
 
     if(associated(dz_8)) then
      call kill(dz_8)
      deallocate(dz_8)
      nullify(dz_8)
     endif
     if(associated(dz_t)) then
      call kill(dz_t)
      deallocate(dz_t)
      nullify(dz_t)
     endif

    package1=.true.
    ndptt=0
    if(present(PACKAGE)) PACKAGE1=PACKAGE
    if(present(ndpt1)) ndptt=ndpt1
   ! !w_p=>W_I                   ! default output, comment out if necessary
    call set_da_pointers
    call init_map_p(NO1,ND1,NP1,ndptt,PACKAGE1)
    call set_in_poly(PACKAGE1)
    call set_in_polyp(PACKAGE1)

    allocate(dz_8(nv))
    call alloc(dz_8)
    allocate(dz_t(nv))
    call alloc(dz_t)

    do i=1,nv
     dz_8(i)=morph(1.0_dp.mono.i)   
    enddo
    do i=1,nv
     dz_t(i)=1.0_dp.mono.i   
    enddo

  end subroutine  init_map_cp

  subroutine init_tpsa_cp(NO1,NV1,PACKAGE)
    implicit none
    integer NO1,NV1,i
    logical(lp),optional :: PACKAGE
    logical(lp) PACKAGE1
     if(associated(dz_8)) then
      call kill(dz_8)
      deallocate(dz_8)
      nullify(dz_8)
     endif
     if(associated(dz_t)) then
      call kill(dz_t)
      deallocate(dz_t)
      nullify(dz_t)
     endif
    package1=.true.
    if(present(PACKAGE)) PACKAGE1=PACKAGE
    !w_p=>W_I                  ! default output, comment out if necessary
    call set_da_pointers
    call init_tpsa_p(NO1,NV1,PACKAGE1)
    call set_in_poly(PACKAGE1)
    call set_in_polyp(PACKAGE1)

    allocate(dz_8(nv))
    call alloc(dz_8)
    allocate(dz_t(nv))
    call alloc(dz_t)

    do i=1,nv
     dz_8(i)=morph(1.0_dp.mono.i)   
    enddo
    do i=1,nv
     dz_t(i)=1.0_dp.mono.i   
    enddo
  end subroutine  init_tpsa_cp


  subroutine set_in_poly(log)
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
    !i_ =cmplx(zero,one,kind=dp)
  end  subroutine set_in_poly

  SUBROUTINE  resetpoly(S2)
    implicit none
    type (complex_8),INTENT(INOUT)::S2

    if(s2%alloc) call killcomplex(s2%t)
    s2%alloc=f
    s2%kind=0
    s2%r=0.0_dp
    !s2%s=one
    !s2%i=0
    !s2%j=0

  END SUBROUTINE resetpoly

  SUBROUTINE  resetpolyn(S2,K)
    implicit none
    type (complex_8),INTENT(INOUT),dimension(:)::S2
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

  SUBROUTINE  resetpoly0(S2)
    implicit none
    type (complex_8),INTENT(INOUT)::S2

    if(s2%alloc) call killcomplex(s2%t)
    s2%alloc=f
    s2%kind=0
    s2%r=0.0_dp
    s2%i=0
    s2%j=0
    s2%s=1.0_dp

  END SUBROUTINE resetpoly0

  !  FUNCTION GETchar( S1, S2 )
  !    implicit none
  !    complex(dp) GETchar
  !    TYPE (complex_8), INTENT (IN) :: S1
  !    CHARACTER(*)  , INTENT (IN) ::  S2
  !
  !    if(s1%kind==m2) then
  !       GETchar=s1%t.sub.s2   !  CHANGE
  !    else
  !       GETchar=s1
  !    endif
  !
  !  END FUNCTION GETchar
  !
  !
  !  FUNCTION GETint( S1, S2 )
  !    implicit none
  !    complex(dp) GETint
  !    TYPE (complex_8), INTENT (IN) :: S1
  !    integer  , INTENT (IN) ::  S2(:)
  !    !  integer localmaster
  !
  !    if(s1%kind==m2) then
  !       ! GETchar%t=s1%t.sub.s2   !  OLD
  !       GETint=s1%t.sub.s2   !  CHANGE
  !    else
  !       GETint=s1
  !    endif
  !
  !  END FUNCTION GETint
  !
  !  FUNCTION GETORDER( S1, S2 )
  !    implicit none
  !    TYPE (complex_8) GETORDER
  !    TYPE (complex_8), INTENT (IN) :: S1
  !    integer  , INTENT (IN) ::  S2
  !
  !    if(s1%kind==m2) then
  !       localmaster=master
  !       call ass(GETORDER)
  !       GETORDER%t=s1%t.sub.s2
  !       master=localmaster
  !    else
  !       GETORDER=s1
  !    endif
  !
  !  END FUNCTION GETORDER
  !
  !
  !  FUNCTION CUTORDER( S1, S2 )
  !    implicit none
  !    TYPE (complex_8) CUTORDER
  !    TYPE (complex_8), INTENT (IN) :: S1
  !    integer  , INTENT (IN) ::  S2
  !
  !    if(s1%kind==m2) then
  !       localmaster=master
  !       call ass(CUTORDER)
  !       CUTORDER%t=s1%t.CUT.s2
  !       master=localmaster
  !    else
  !       CUTORDER=s1
  !    endif
  !
  !  END FUNCTION CUTORDER

  FUNCTION GETCHARnd2( S1, S2 )
    implicit none
    TYPE (complex_8) GETCHARnd2
    TYPE (complex_8), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    integer localmaster
    type(complextaylor) t

    localmaster=master
    call ass(GETCHARnd2)
    call alloc(t)
    t=s1

    t=t.par.s2
    GETCHARnd2%t=t



    call kill(t)
    master=localmaster

    !    if(s1%kind==m2) then
    !       localmaster=master
    !       call ass(GETCHARnd2)
    !       GETCHARnd2%t=s1%t.par.s2
    !       master=localmaster
    !    endif

  END FUNCTION GETCHARnd2

  FUNCTION GETintnd2( S1, S2 )
    implicit none
    TYPE (complex_8) GETintnd2
    TYPE (complex_8), INTENT (IN) :: S1
    integer, INTENT (IN) ::  S2(:)
    integer localmaster
    type(complextaylor) t

    !    if(s1%kind==m2) then
    !       localmaster=master
    !       call ass(GETintnd2)
    !       GETintnd2%t=s1%t.par.s2
    !       master=localmaster
    !    endif


    localmaster=master
    call ass(GETintnd2)
    call alloc(t)
    t=s1
    t=t.par.s2
    GETintnd2%t=t


    call kill(t)

    master=localmaster

  END FUNCTION GETintnd2


  FUNCTION GETchar( S1, S2 )
    implicit none
    complex(dp) GETchar
    TYPE (complex_8), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    integer i,j
    !  integer localmaster

    GETchar=0.0_dp
    if(s1%kind==m2) then
       ! GETchar%t=s1%t.sub.s2   !  OLD
       GETchar=s1%t.sub.s2   !  CHANGE
    elseif(s1%kind==m1) then
       GETchar=s1
       do i=1,len_trim(s2)
          CALL  CHARINT(s2(i:i),j)
          if(j/=0) then
             GETchar=0.0_dp
             exit
          endif
       enddo

    endif

  END FUNCTION GETchar

  FUNCTION GETint( S1, S2 )
    implicit none
    complex(dp) GETint
    TYPE (complex_8), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    integer i

    GETint=0.0_dp
    if(s1%kind==m2) then
       GETint=s1%t.sub.s2   !  CHANGE
    elseif(s1%kind==m1) then
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
    TYPE (complex_8) GETORDER
    TYPE (complex_8), INTENT (IN) :: S1
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
    TYPE (complex_8) CUTORDER
    TYPE (complex_8), INTENT (IN) :: S1
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




  SUBROUTINE  A_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (complex_8),INTENT(INout)::S1
    type (complex_8),optional, INTENT(INout):: S2,s3,s4,s5,s6,s7,s8,s9,s10
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

  SUBROUTINE  K_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (complex_8),INTENT(INout)::S1
    type (complex_8),optional, INTENT(INout):: S2,s3,s4,s5,s6,s7,s8,s9,s10
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

  SUBROUTINE  resetpolyn0(S2,K)
    implicit none
    type (complex_8),INTENT(INOUT),dimension(:)::S2
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

  SUBROUTINE  resetpoly_R(S2,FL)  !   STAYS REAL
    implicit none
    type (complex_8),INTENT(INOUT)::S2
    logical(lp),INTENT(IN)::FL

    if(s2%alloc) call killcomplex(s2%t)
    s2%alloc=f
    s2%kind=1
    s2%r=0.0_dp
    IF(.NOT.FL) THEN
       s2%i=0
       s2%j=0
       s2%s=1.0_dp
    ENDIF

  END SUBROUTINE resetpoly_R

  SUBROUTINE  resetpoly_RN(S2,FL,K)
    implicit none
    type (complex_8),INTENT(INOUT),dimension(:)::S2
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



  SUBROUTINE  allocpoly(S2)
    implicit none
    type (complex_8),INTENT(INOUT)::S2

    !  if(s2%alloc) call killcomplex(s2%t)
    s2%alloc=f
    s2%kind=1
    s2%r=0.0_dp
    s2%i=0
    s2%j=0
!    s2%g=0
    s2%s=1.0_dp

  END SUBROUTINE allocpoly

  SUBROUTINE  allocpolyn(S2,K)
    implicit none
    type (complex_8),INTENT(INOUT),dimension(:)::S2
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

  SUBROUTINE  printpoly(S2,mf)
    implicit none
    type (complex_8),INTENT(INOUT)::S2
    integer ipause,mypauses
    integer,optional :: mf
    integer i
    character(255) line
    i=6
    if(present(mf)) i=mf
          write(i,*) " printing a complex polymorph (complex_8)"
    if(s2%kind/=0) then

       select  case (s2%kind)
       case(m1)
          write(i,*)  s2%r
       case(m2)
          call printcomplex(S2%t,i)
       case(m3)

          if(s2%i>0.and.s2%j>0) then
             write(line,*) s2%r,"  +",s2%s,"  (x_",s2%i,"+ i","*x_",s2%j,")"
          elseif(s2%i>0)then
             write(line,*) s2%r,"  +",s2%s,"  (x_",s2%i,")"

           elseif(s2%j>0)then
             write(line,*) s2%r,"  +",s2%s,"  ( ","i*x_",s2%j,")"

           else
            write(line,*) s2%r
          endif
             call context(line,maj=.false.)
             write(i,'(a)') adjustr(line(1:len_trim(line)))
        
       end   select
    else

       line=" Warning not defined "
       ipause=mypauses(0,line)

    endif

  END SUBROUTINE printpoly

  SUBROUTINE  print6(S1,mf)
    implicit none
    type (complex_8),INTENT(INout)::S1(:)
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

  SUBROUTINE EQUAL(S2,S1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    type (complex_8),INTENT(IN)::S1
    !    integer localmaster

    if(s1%kind==0) then
       line=" You are putting kind=0  into something"
       ipause=mypauses(0,line)
    endif
    if(s2%kind==3.and.(.not.setknob)) then
       line=" You are putting something  into a knob kind=3"
       ipause=mypauses(0,line)
    endif

    if (s2%kind>0) then       !   S2 exist
       if(S2%kind==S1%kind) then
          select case(S1%kind)
          case(m1)
             S2%R=S1%R
          case(m2)
             !             localmaster=master
             call check_snake
             !2002.12.26 master=0
             S2%t=S1%t
             !             master=localmaster
          case(m3)
             s2%r=S1%r ! Knob stays a knob and real stays real 2002.10.9
             !             !w_p=0
             !             !w_p%nc=2
             !             !w_p%fc='((1X,A72,/,1x,a72))'
             !               write(6,*) " You are putting kind=3 (TPSA) into another kind=3"
             !               write(6,*) " The left handside is not a knob!"
             !             ! call !write_e(0)
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
             !2002.12.26 master=0
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
             !2002.12.26 master=0
             if(knob) then
                call varck1(s1)
                S2%t=varc1
             else
                s2%r=S1%r
                s2%kind=1
             endif
             !         call varkind3(S1)
             !         S2%t=S1%t
             !         call unvarkind3(S1)
             !             master=localmaster
          end select
       endif

    else        !   S2 does not exist



       if(S1%kind==1) then  ! what is s1
          if(s2%i==0) then
             S2%R=S1%R
             s2%kind=1
          elseif((s2%i>0.and.s2%i<=nv).and.(s2%j>0.and.s2%j<=nv))  then
             call alloc(s2%t)
             s2%t=(/S1%R,S2%S/).var.(/s2%i,s2%j/)
             !             call var(s2%t,S1%R,S2%S,s2%i,s2%j)
             !      s2%i=0
             !      s2%j=0
             s2%kind=2
             s2%alloc=t
          else
             stop 777
          endif

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
       endif    ! end of what is s1

    endif   ! S2  does not exist


  END SUBROUTINE EQUAL

  SUBROUTINE  EQUALRP(S2,S1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    type (REAL_8),INTENT(IN)::S1
    !    integer localmaster

    if(s1%kind==0) then
       line=" You are putting kind=0  into something"
       ipause=mypauses(0,line)
    endif
    if(s2%kind==3.and.(.not.setknob)) then
       line=" You are putting something  into a knob kind=3"
       ipause=mypauses(0,line)
    endif

    if (s2%kind>0) then       !   S2 exist
       if(S2%kind==S1%kind) then
          select case(S1%kind)
          case(m1)
             S2%R=S1%R
          case(m2)
             !             localmaster=master
             call check_snake
             !2002.12.26 master=0
             S2%t=S1%t
             !             master=localmaster
          case(m3)
             s2%r=S1%r ! Knob stays a knob and real stays real 2002.10.9
             !             !w_p=0
             !             !w_p%nc=2
             !             !w_p%fc='((1X,A72,/,1x,a72))'
             !               write(6,*) " You are putting kind=3 (TPSA) into another kind=3"
             !               write(6,*) " The left handside is not a knob!"
             !             ! call !write_e(0)
          end select
       elseif(S2%kind>S1%kind ) then
          if(S1%kind/=2) then
             s2%r=S1%r
             if(s2%kind/=3) s2%kind=1 !Knob stays a knob and real stays real 2002.10.9
          else
             s2%r=S1%t.sub.'0'  ! setting a kind=3
          endif
          !          !     S2%t=S1%r    all removed 2002.10.9
          !          s2%r=S1%r
          !          s2%kind=1
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
             !2002.12.26 master=0
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
             !2002.12.26 master=0
             if(knob) then
                call varfk1(s1)
                S2%t=varf1
             else
                s2%r=S1%r
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
          elseif((s2%i>0.and.s2%i<=nv).and.(s2%j>0.and.s2%j<=nv))  then
             call alloc(s2%t)
             s2%t=(/cmplx(S1%R,0.0_dp,kind=dp),S2%S/).var.(/s2%i,s2%j/)
             !             call var(s2%t,cmplx(S1%R,zero,kind=dp),S2%S,s2%i,s2%j)
             !      s2%i=0
             !      s2%j=0
             s2%kind=2
             s2%alloc=t
          else
             stop 777
          endif

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
       endif    ! end of what is s1

    endif   ! S2  does not exist


  END SUBROUTINE EQUALRP

  SUBROUTINE  RPEQUAL(S2,S1)
    implicit none
    integer ipause, mypauses
    type (REAL_8),INTENT(inOUT)::S2
    type (complex_8),INTENT(IN)::S1
    !    integer localmaster

    if(s1%kind==0) then
       line=" You are putting kind=0  into something"
       ipause=mypauses(0,line)
    endif
    if(s2%kind==3.and.(.not.setknob)) then
       line= " You are putting something  into a knob kind=3"
       ipause=mypauses(0,line)
    endif

    if (s2%kind>0) then       !   S2 exist
       if(S2%kind==S1%kind) then
          select case(S1%kind)
          case(m1)
             S2%R=S1%R
          case(m2)
             !             localmaster=master
             call check_snake
             !2002.12.26 master=0
             S2%t=S1%t
             !             master=localmaster
          case(m3)
             s2%r=S1%r ! Knob stays a knob and real stays real 2002.10.9
             !             !w_p=0
             !             !w_p%nc=2
             !             !w_p%fc='((1X,A72,/,1x,a72))'
             !               write(6,*) " You are putting kind=3 (TPSA) into another kind=3"
             !               write(6,*) " The left handside is not a knob!"
             !             ! call !write_e(0)
          end select
       elseif(S2%kind>S1%kind ) then
          if(S1%kind/=2) then
             s2%r=S1%r
             if(s2%kind/=3) s2%kind=1 !Knob stays a knob and real stays real 2002.10.9
          else
             s2%r=S1%t.sub.'0'  ! setting a kind=3
          endif
          !          !     S2%t=S1%r
          !          s2%r=S1%r
          !          s2%kind=1
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
             !2002.12.26 master=0
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
             !2002.12.26 master=0
             if(knob) then
                call varck1(s1)
                S2%t=varc1
             else
                S2%r=s1%r
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
          elseif((s2%i>0.and.s2%i<=nv))  then
             call alloc(s2%t)
             s2%t=(/REAL(S1%R,kind=DP),S2%S/).var.s2%i
             !             call var(s2%t,REAL(S1%R,kind=DP),S2%S,s2%i)
             !      s2%i=0
             s2%kind=2
             s2%alloc=t
          else
             stop 777
          endif

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
       endif    ! end of what is s1

    endif   ! S2  does not exist


  END SUBROUTINE RPEQUAL

  FUNCTION drealt( S1 )
    implicit none
    TYPE (real_8) drealt
    TYPE (complex_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       drealt%r=REAL(s1%r,kind=DP)
       drealt%kind=1
    case(m2)
       localmaster=master
       call assp(drealt)
       drealt%t=s1%t%r
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call assp(drealt)
          call varck1(s1)
          drealt%t=varc1%r
          master=localmaster
       else
          drealt%r=REAL(s1%r,kind=DP)
          drealt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in drealt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END FUNCTION drealt

  FUNCTION dimagt( S1 )
    implicit none
    TYPE (real_8) dimagt
    TYPE (complex_8), INTENT (IN) :: S1
    integer localmaster
    select case(s1%kind)
    case(m1)
       dimagt%r=aimag(s1%r)
       dimagt%kind=1
    case(m2)
       localmaster=master
       call assp(dimagt)
       dimagt%t=s1%t%i
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call assp(dimagt)
          call varck1(s1)
          dimagt%t=varc1%i
          master=localmaster
       else
          dimagt%r=aimag(s1%r)
          dimagt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dimagt "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
  END FUNCTION dimagt

  FUNCTION dcmplxt( S1,s2 )
    implicit none
    TYPE (complex_8) dcmplxt
    TYPE (real_8), INTENT (IN) :: S1,s2
    integer localmaster
    select case(s1%kind+ms*s2%kind)
    case(m11)
       dcmplxt%r=cmplx(s1%r,s2%r,kind=dp)
       dcmplxt%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(dcmplxt)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          dcmplxt%t%r=s1%t
          dcmplxt%t%i=s2%r
       case(m12)
          dcmplxt%t%r=s1%r
          dcmplxt%t%i=s2%t
       case(m22)
          dcmplxt%t%r=s1%t
          dcmplxt%t%i=s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(dcmplxt)
             call varfk1(s1)
             dcmplxt%t%r=varf1
             dcmplxt%t%i=s2%r
             master=localmaster
          else
             dcmplxt%r=cmplx(s1%r,s2%r,kind=dp)
             dcmplxt%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(dcmplxt)
             call varfk2(s2)
             dcmplxt%t%i=varf2
             dcmplxt%t%r=s1%r
             master=localmaster
          else
             dcmplxt%r=cmplx(s1%r,s2%r,kind=dp)
             dcmplxt%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(dcmplxt)
          if(knob) then
             call varfk1(s1)
             dcmplxt%t%r=varf1
             dcmplxt%t%i=s2%t
          else
             dcmplxt%t%r=s1%r
             dcmplxt%t%i=s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(dcmplxt)
          if(knob) then
             call varfk2(s2)
             dcmplxt%t%i=varf2
             dcmplxt%t%r=s1%t
          else
             dcmplxt%t%i=s2%r
             dcmplxt%t%r=s1%t
          endif
          master=localmaster

       case(m33)
          if(knob) then
             localmaster=master
             call ass(dcmplxt)
             call varFk1(s1)
             call varFk2(s2)
             dcmplxt%t%r=varf1
             dcmplxt%t%i=varf2
             master=localmaster
          else
             dcmplxt%r=cmplx(s1%r,s2%r,kind=dp)
             dcmplxt%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in dcmplxt "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dcmplxt


  SUBROUTINE  complexEQUAL(S2,S1)
    implicit none
    complex(dp) ,INTENT(inout)::S2
    type (complex_8),INTENT(IN)::S1
    !    integer localmaster


    select case(S1%kind)
    case(m1)
       S2=S1%R
    case(m2)
       !       localmaster=master
       call check_snake
       !2002.12.26 master=0
       S2=S1%t.sub.'0'
       !       master=localmaster
    case(m3)
       S2=S1%r
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in complexEQUAL "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select

  END SUBROUTINE complexEQUAL

  SUBROUTINE  EQUALcomplext(S2,S1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    type (complextaylor),INTENT(IN)::S1
    !    integer localmaster

    if(s2%kind==3.and.(.not.setknob)) then
       line= " You are putting something  into a knob kind=3"
       ipause=mypauses(0,line)
    endif
    !    localmaster=master
    call check_snake
    !2002.12.26 master=0
    if(s2%kind/=3) then
       if(.not.s2%alloc) then
          call alloc(s2%t)
          s2%alloc=t
       endif
       s2%kind=2
       s2%t=s1
    else
       s2%r=S1.sub.'0' ! 2002.10.9
    endif
    !    master=localmaster


  END SUBROUTINE EQUALcomplext

  SUBROUTINE  complextEQUAL(S1,S2)
    implicit none
    type (complex_8),INTENT(in)::S2
    type (complextaylor),INTENT(inout)::S1
    !    integer localmaster



    select case(S2%kind)
    case(m1)
       S1=S2%R
    case(m2)
       !       localmaster=master
       call check_snake
       !2002.12.26 master=0
       S1=S2%t
       !       master=localmaster
    case(m3)
       !       localmaster=master
       call check_snake
       !2002.12.26 master=0
       if(knob) then
          call varck2(s2)
          S1=varc2
       else
          S1=S2%R
       endif
       !       master=localmaster
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in complexEQUAL "
         write(6,*) "s2%kind   "
       !w_p=(/s2%kind  /)
       ! call !write_e(0)
    end select

  END SUBROUTINE complextEQUAL

  SUBROUTINE  Dequaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::R1

    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=r1
          return
       else
          line="Forbidden in Dequaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF

    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1) then
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
       elseif((s2%i>0.and.s2%i<=nv).and.(s2%j>0.and.s2%j<=nv))  then
          call alloc(s2%t)
          s2%t=(/cmplx(R1,0.0_dp,kind=dp),s2%s/).var.(/s2%i,s2%j/)
          !          call var(s2%t,cmplx(R1,zero,kind=dp),s2%s,s2%i,s2%j)
          !      s2%i=0
          !      s2%j=0
          s2%kind=2
          s2%alloc=t
       else
          stop 777
       endif


    endif   ! S2 not allocated
  END SUBROUTINE Dequaldacon

  SUBROUTINE  equaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    real(sp),INTENT(IN)::R1

    if(real_warning) call real_stop
    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=r1
          return
       else
          line="Forbidden in equaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF

    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1) then
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
       elseif((s2%i>0.and.s2%i<=nv).and.(s2%j>0.and.s2%j<=nv))  then
          call alloc(s2%t)
          s2%t=(/cmplx(REAL(R1,kind=DP),0.0_dp,kind=dp),S2%S/).var.(/s2%i,s2%j/)
          !          call var(s2%t,cmplx(REAL(R1,kind=DP),zero,kind=dp),S2%S,s2%i,s2%j)
          !      s2%i=0
          !      s2%j=0
          s2%kind=2
          s2%alloc=t
       else
          stop 777
       endif


    endif   ! S2 not allocated
  END SUBROUTINE equaldacon

  SUBROUTINE  iequaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    integer,INTENT(IN)::R1

    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=r1
          return
       else
          line= "Forbidden in iequaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF

    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1) then
          S2%R=REAL(R1,kind=DP)
       else
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
       elseif((s2%i>0.and.s2%i<=nv).and.(s2%j>0.and.s2%j<=nv))  then
          call alloc(s2%t)
          s2%t=(/cmplx(REAL(R1,kind=DP),0.0_dp,kind=dp),s2%s/).var.(/s2%i,s2%j/)
          !          call var(s2%t,cmplx(REAL(R1,kind=DP),zero,kind=dp),s2%s,s2%i,s2%j)
          !      s2%i=0
          !      s2%j=0
          s2%kind=2
          s2%alloc=t
       else
          stop 777
       endif


    endif   ! S2 not allocated
  END SUBROUTINE iequaldacon


  SUBROUTINE  cequaldacon(S2,R1)
    implicit none
    integer ipause, mypauses
    type (complex_8),INTENT(inOUT)::S2
    complex(dp),INTENT(IN)::R1

    IF(S2%KIND==M3) THEN
       if(setknob) then
          s2%r=r1
          return
       else
          line="Forbidden in cequaldacon: s2 is a knob "
          ipause=mypauses(0,line)
       endif
    ENDIF

    if (s2%kind/=0) then       !   S2 exist
       if(S2%kind==1) then
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
       elseif((s2%i>0.and.s2%i<=nv).and.(s2%j>0.and.s2%j<=nv))  then
          call alloc(s2%t)
          s2%t=(/R1,s2%s/).var.(/s2%i,s2%j/)
          !          call var(s2%t,R1,s2%s,s2%i,s2%j)
          !      s2%i=0
          !      s2%j=0
          s2%kind=2
          s2%alloc=t
       else
          stop 777
       endif


    endif   ! S2 not allocated
  END SUBROUTINE cequaldacon

  FUNCTION add( S1, S2 )
    implicit none
    TYPE (complex_8) add
    TYPE (complex_8), INTENT (IN) :: S1, S2
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
             call varck1(s1)
             add%t= varc1+s2%r
             master=localmaster
          else
             add%r=s1%r+s2%r
             add%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(add)
             call varck2(s2)
             add%t= s1%r+varc2
             master=localmaster
          else
             add%r=s1%r+s2%r
             add%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(add)
          if(knob) then
             call varck1(s1)
             add%t= varc1+s2%t
          else
             add%t=s1%r+s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(add)
          if(knob) then
             call varck2(s2)
             add%t= s1%t+varc2
          else
             add%t= s1%t+s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(add)
             call varck1(s1)
             call varck2(s2)
             add%t= varc1+varc2
             master=localmaster
          else
             add%r=s1%r+s2%r
             add%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in add "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION add

  FUNCTION unaryADD( S1 )
    implicit none
    TYPE (complex_8) unaryADD
    TYPE (complex_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       unaryADD%r=s1%r
       unaryADD%kind=1
    case(m2)
       localmaster=master
       call ass(unaryADD)
       unaryADD%t=s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(unaryADD)
          call varck1(s1)
          unaryADD%t=varc1
          master=localmaster
       else
          unaryADD%r=s1%r
          unaryADD%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in unaryADD "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION unaryADD

  FUNCTION daddsc( S1, S2 )
    implicit none
    TYPE (complex_8) daddsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          daddsc%t= varc1+s2
          master=localmaster
       else
          daddsc%r=s1%r+s2
          daddsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in daddsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION daddsc

  FUNCTION caddsc( S1, S2 )
    implicit none
    TYPE (complex_8) caddsc
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       caddsc%r=s1%r+s2
       caddsc%kind=1
    case(m2)
       localmaster=master
       call ass(caddsc)
       caddsc%t= s1%t+s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(caddsc)
          call varck1(s1)
          caddsc%t= varc1+s2
          master=localmaster
       else
          caddsc%r=s1%r+s2
          caddsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in caddsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION caddsc

  FUNCTION dscadd( S2, S1 )
    implicit none
    TYPE (complex_8) dscadd
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dscadd%t= varc1+s2
          master=localmaster
       else
          dscadd%r=s1%r+s2
          dscadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dscadd "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dscadd

  FUNCTION cscadd( S2, S1 )
    implicit none
    TYPE (complex_8) cscadd
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cscadd%r=s1%r+s2
       cscadd%kind=1
    case(m2)
       localmaster=master
       call ass(cscadd)
       cscadd%t= s1%t+s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cscadd)
          call varck1(s1)
          cscadd%t= varc1+s2
          master=localmaster
       else
          cscadd%r=s1%r+s2
          cscadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cscadd "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cscadd

  FUNCTION addsc( S1, S2 )
    implicit none
    TYPE (complex_8) addsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          addsc%t= varc1+REAL(s2,kind=DP)
          master=localmaster
       else
          addsc%r=s1%r+REAL(s2,kind=DP)
          addsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in addsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION addsc

  FUNCTION scadd( S2, S1  )
    implicit none
    TYPE (complex_8) scadd
    TYPE (complex_8), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: S2
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
          call varck1(s1)
          scadd%t= varc1+REAL(s2,kind=DP)
          master=localmaster
       else
          scadd%r=s1%r+REAL(s2,kind=DP)
          scadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in scadd "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION scadd

  FUNCTION iaddsc( S1, S2 )
    implicit none
    TYPE (complex_8) iaddsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          iaddsc%t= varc1+REAL(s2,kind=DP)
          master=localmaster
       else
          iaddsc%r=s1%r+REAL(s2,kind=DP)
          iaddsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in iaddsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION iaddsc

  FUNCTION iscadd( S2, S1 )
    implicit none
    TYPE (complex_8) iscadd
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          iscadd%t= varc1+REAL(s2,kind=DP)
          master=localmaster
       else
          iscadd%r=s1%r+REAL(s2,kind=DP)
          iscadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in iscadd "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION iscadd

  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (complex_8) subs
    TYPE (complex_8), INTENT (IN) :: S1, S2
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
             call varck1(s1)
             subs%t= varc1-s2%r
             master=localmaster
          else
             subs%r=s1%r-s2%r
             subs%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(subs)
             call varck2(s2)
             subs%t= s1%r-varc2
             master=localmaster
          else
             subs%r=s1%r-s2%r
             subs%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(subs)
          if(knob) then
             call varck1(s1)
             subs%t= varc1-s2%t
          else
             subs%t=s1%r-s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(subs)
          if(knob) then
             call varck2(s2)
             subs%t= s1%t-varc2
          else
             subs%t= s1%t-s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(subs)
             call varck1(s1)
             call varck2(s2)
             subs%t= varc1-varc2
             master=localmaster
          else
             subs%r=s1%r-s2%r
             subs%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in subs "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION subs

  FUNCTION unarySUB( S1 )
    implicit none
    TYPE (complex_8) unarySUB
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          unarySUB%t= -varc1
          master=localmaster
       else
          unarySUB%r=-s1%r
          unarySUB%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in unarySUB "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION unarySUB

  FUNCTION dsubsc( S1, S2 )
    implicit none
    TYPE (complex_8) dsubsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dsubsc%t= varc1-s2
          master=localmaster
       else
          dsubsc%r=s1%r-s2
          dsubsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dsubsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dsubsc

  FUNCTION dscsub( S2, S1 )
    implicit none
    TYPE (complex_8) dscsub
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dscsub%t=s2-varc1
          master=localmaster
       else
          dscsub%r=s2-s1%r
          dscsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dscsub "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dscsub

  FUNCTION csubsc( S1, S2 )
    implicit none
    TYPE (complex_8) csubsc
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       csubsc%r=s1%r-s2
       csubsc%kind=1
    case(m2)
       localmaster=master
       call ass(csubsc)
       csubsc%t= s1%t-s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(csubsc)
          call varck1(s1)
          csubsc%t= varc1-s2
          master=localmaster
       else
          csubsc%r=s1%r-s2
          csubsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in csubsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION csubsc

  FUNCTION cscsub( S2, S1 )
    implicit none
    TYPE (complex_8) cscsub
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cscsub%r=s2-s1%r
       cscsub%kind=1
    case(m2)
       localmaster=master
       call ass(cscsub)
       cscsub%t=s2-s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cscsub)
          call varck1(s1)
          cscsub%t=s2-varc1
          master=localmaster
       else
          cscsub%r=s2-s1%r
          cscsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cscsub "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cscsub

  FUNCTION subsc( S1, S2 )
    implicit none
    TYPE (complex_8) subsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          subsc%t= varc1-REAL(s2,kind=DP)
          master=localmaster
       else
          subsc%r=s1%r-REAL(s2,kind=DP)
          subsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in subsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION subsc

  FUNCTION scsub( S2, S1 )
    implicit none
    TYPE (complex_8) scsub
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          scsub%t=REAL(s2,kind=DP)-varc1
          master=localmaster
       else
          scsub%r=REAL(s2,kind=DP)-s1%r
          scsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in scsub "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION scsub

  FUNCTION isubsc( S1, S2 )
    implicit none
    TYPE (complex_8) isubsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          isubsc%t= varc1-REAL(s2,kind=DP)
          master=localmaster
       else
          isubsc%r=s1%r-REAL(s2,kind=DP)
          isubsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in isubsc "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION isubsc

  FUNCTION iscsub( S2, S1 )
    implicit none
    TYPE (complex_8) iscsub
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          iscsub%t=REAL(s2,kind=DP)-varc1
          master=localmaster
       else
          iscsub%r=REAL(s2,kind=DP)-s1%r
          iscsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in iscsub "
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION iscsub

  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (complex_8) mul
    TYPE (complex_8), INTENT (IN) :: S1, S2
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
             call varck1(s1)
             mul%t= varc1*s2%r
             master=localmaster
          else
             mul%r=s1%r*s2%r
             mul%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(mul)
             call varck2(s2)
             mul%t= s1%r*varc2
             master=localmaster
          else
             mul%r=s1%r*s2%r
             mul%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(mul)
          if(knob) then
             call varck1(s1)
             mul%t= varc1*s2%t
          else
             mul%t=s1%r*s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(mul)
          if(knob) then
             call varck2(s2)
             mul%t= s1%t*varc2
          else
             mul%t= s1%t*s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(mul)
             call varck1(s1)
             call varck2(s2)
             mul%t= varc1*varc2
             master=localmaster
          else
             mul%r=s1%r*s2%r
             mul%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in mul "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION mul

  FUNCTION pmul( S1, S2 )
    implicit none
    TYPE (complex_8) pmul
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       pmul%r=s1%r*s2%r
       pmul%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(pmul)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          pmul%t= s1%t*s2%r
       case(m12)
          pmul%t= s1%r*s2%t
       case(m22)
          pmul%t= s1%t*s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(pmul)
             call varfk1(s1)
             pmul%t= varf1*s2%r
             master=localmaster
          else
             pmul%r=s1%r*s2%r
             pmul%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(pmul)
             call varck2(s2)
             pmul%t= s1%r*varc2
             master=localmaster
          else
             pmul%r=s1%r*s2%r
             pmul%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(pmul)
          if(knob) then
             call varfk1(s1)
             pmul%t= varf1*s2%t
          else
             pmul%t=s1%r*s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(pmul)
          if(knob) then
             call varck2(s2)
             pmul%t= s1%t*varc2
          else
             pmul%t=s1%t*s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(pmul)
             call varfk1(s1)
             call varck2(s2)
             pmul%t= varf1*varc2
             master=localmaster
          else
             pmul%r=s1%r*s2%r
             pmul%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in pmul "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION pmul

  FUNCTION mulp( S1, S2 )
    implicit none
    TYPE (complex_8) mulp
    TYPE (real_8), INTENT (IN) :: S2
    TYPE (complex_8), INTENT (IN) ::  S1
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       mulp%r=s1%r*s2%r
       mulp%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(mulp)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          mulp%t= s1%t*s2%r
       case(m12)
          mulp%t= s1%r*s2%t
       case(m22)
          mulp%t= s1%t*s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(mulp)
             call varck1(s1)
             mulp%t= varc1*s2%r
             master=localmaster
          else
             mulp%r=s1%r*s2%r
             mulp%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(mulp)
             call varfk2(s2)
             mulp%t= s1%r*varf2
             master=localmaster
          else
             mulp%r=s1%r*s2%r
             mulp%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(mulp)
          if(knob) then
             call varck1(s1)
             mulp%t= varc1*s2%t
          else
             mulp%t= s1%r*s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(mulp)
          if(knob) then
             call varfk2(s2)
             mulp%t= s1%t*varf2
          else
             mulp%t= s1%t*s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(mulp)
             call varck1(s1)
             call varfk2(s2)
             mulp%t= varc1*varf2
             master=localmaster
          else
             mulp%r=s1%r*s2%r
             mulp%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in mulp "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION mulp

  FUNCTION padd( S1, S2 )
    implicit none
    TYPE (complex_8) padd
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       padd%r=s1%r+s2%r
       padd%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(padd)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          padd%t= s1%t+s2%r
       case(m12)
          padd%t= s1%r+s2%t
       case(m22)
          padd%t= s1%t+s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(padd)
             call varfk1(s1)
             padd%t= varf1+s2%r
             master=localmaster
          else
             padd%r=s1%r+s2%r
             padd%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(padd)
             call varck2(s2)
             padd%t= s1%r+varc2
             master=localmaster
          else
             padd%r=s1%r+s2%r
             padd%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(padd)
          if(knob) then
             call varfk1(s1)
             padd%t= varf1+s2%t
          else
             padd%t=s1%r+s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(padd)
          if(knob) then
             call varck2(s2)
             padd%t= s1%t+varc2
          else
             padd%t=s1%t+s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(padd)
             call varfk1(s1)
             call varck2(s2)
             padd%t= varf1+varc2
             master=localmaster
          else
             padd%r=s1%r+s2%r
             padd%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in padd "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION padd

  FUNCTION addp( S2, S1 )
    implicit none
    TYPE (complex_8) addp
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       addp%r=s1%r+s2%r
       addp%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(addp)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          addp%t= s1%t+s2%r
       case(m12)
          addp%t= s1%r+s2%t
       case(m22)
          addp%t= s1%t+s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(addp)
             call varfk1(s1)
             addp%t= varf1+s2%r
             master=localmaster
          else
             addp%r=s1%r+s2%r
             addp%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(addp)
             call varck2(s2)
             addp%t= s1%r+varc2
             master=localmaster
          else
             addp%r=s1%r+s2%r
             addp%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(addp)
          if(knob) then
             call varfk1(s1)
             addp%t= varf1+s2%t
          else
             addp%t=s1%r+s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(addp)
          if(knob) then
             call varck2(s2)
             addp%t= s1%t+varc2
          else
             addp%t=s1%t+s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(addp)
             call varfk1(s1)
             call varck2(s2)
             addp%t= varf1+varc2
             master=localmaster
          else
             addp%r=s1%r+s2%r
             addp%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in addp "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION addp

  FUNCTION psub( S1, S2 )
    implicit none
    TYPE (complex_8) psub
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       psub%r=s1%r-s2%r
       psub%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(psub)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          psub%t= s1%t-s2%r
       case(m12)
          psub%t= s1%r-s2%t
       case(m22)
          psub%t= s1%t-s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(psub)
             call varfk1(s1)
             psub%t= varf1-s2%r
             master=localmaster
          else
             psub%r=s1%r-s2%r
             psub%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(psub)
             call varck2(s2)
             psub%t= s1%r-varc2
             master=localmaster
          else
             psub%r=s1%r-s2%r
             psub%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(psub)
          if(knob) then
             call varfk1(s1)
             psub%t= varf1-s2%t
          else
             psub%t=s1%r-s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(psub)
          if(knob) then
             call varck2(s2)
             psub%t= s1%t-varc2
          else
             psub%t=s1%t-s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(psub)
             call varfk1(s1)
             call varck2(s2)
             psub%t= varf1-varc2
             master=localmaster
          else
             psub%r=s1%r-s2%r
             psub%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in psub "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION psub

  FUNCTION subp(S2 , S1 )
    implicit none
    TYPE (complex_8) subp
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       subp%r=s2%r-s1%r
       subp%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(subp)
       select case(s1%kind+ms*s2%kind)
       case(m12)
          subp%t= s2%t-s1%r
       case(m21)
          subp%t= s2%r-s1%t
       case(m22)
          subp%t= s2%t-s1%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(subp)
             call varfk1(s1)
             subp%t= s2%r-varf1
             master=localmaster
          else
             subp%r=s2%r-s1%r
             subp%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(subp)
             call varck2(s2)
             subp%t= s1%r-varc2
             master=localmaster
          else
             subp%r=s2%r-s1%r
             subp%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(subp)
          if(knob) then
             call varfk1(s1)
             subp%t= s2%t-varf1
          else
             subp%t=s2%r-s1%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(subp)
          if(knob) then
             call varck2(s2)
             subp%t= varc2-s1%t
          else
             subp%t=s2%t-s1%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(subp)
             call varfk1(s1)
             call varck2(s2)
             subp%t= varc2-varf1
             master=localmaster
          else
             subp%r=s2%r-s1%r
             subp%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in subp "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION subp

  FUNCTION pdiv( S1, S2 )
    implicit none
    TYPE (complex_8) pdiv
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       pdiv%r=s1%r/s2%r
       pdiv%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(pdiv)
       select case(s1%kind+ms*s2%kind)
       case(m21)
          pdiv%t= s1%t/s2%r
       case(m12)
          pdiv%t= s1%r/s2%t
       case(m22)
          pdiv%t= s1%t/s2%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(pdiv)
             call varfk1(s1)
             pdiv%t= varf1/s2%r
             master=localmaster
          else
             pdiv%r=s1%r/s2%r
             pdiv%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(pdiv)
             call varck2(s2)
             pdiv%t= s1%r/varc2
             master=localmaster
          else
             pdiv%r=s1%r/s2%r
             pdiv%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(pdiv)
          if(knob) then
             call varfk1(s1)
             pdiv%t= varf1/s2%t
          else
             pdiv%t=s1%r/s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(pdiv)
          if(knob) then
             call varck2(s2)
             pdiv%t= s1%t/varc2
          else
             pdiv%t=s1%t/s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(pdiv)
             call varfk1(s1)
             call varck2(s2)
             pdiv%t= varf1/varc2
             master=localmaster
          else
             pdiv%r=s1%r/s2%r
             pdiv%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in pdiv "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION pdiv

  FUNCTION divp(S2 , S1 )
    implicit none
    TYPE (complex_8) divp
    TYPE (real_8), INTENT (IN) :: S1
    TYPE (complex_8), INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind+ms*s2%kind)
    case(m11)
       divp%r=s2%r/s1%r
       divp%kind=1
    case(m12,m21,m22)
       localmaster=master
       call ass(divp)
       select case(s1%kind+ms*s2%kind)
       case(m12)
          divp%t= s2%t/s1%r
       case(m21)
          divp%t= s2%r/s1%t
       case(m22)
          divp%t= s2%t/s1%t
       end select
       master=localmaster
    case(m13,m31,m32,m23,m33)
       select case(s1%kind+ms*s2%kind)
       case(m31)
          if(knob) then
             localmaster=master
             call ass(divp)
             call varfk1(s1)
             divp%t= s2%r/varf1
             master=localmaster
          else
             divp%r=s2%r/s1%r
             divp%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(divp)
             call varck2(s2)
             divp%t= varc2/s1%r
             master=localmaster
          else
             divp%r=s2%r/s1%r
             divp%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(divp)
          if(knob) then
             call varfk1(s1)
             divp%t= s2%t/varf1
          else
             divp%t=s2%t/s1%r
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(divp)
          if(knob) then
             call varck2(s2)
             divp%t= varc2/s1%t
          else
             divp%t=s2%r/s1%t
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(divp)
             call varfk1(s1)
             call varck2(s2)
             divp%t= varc2/varf1
             master=localmaster
          else
             divp%r=s2%r/s1%r
             divp%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in divp "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION divp


  FUNCTION cmulsc( S1, S2 )
    implicit none
    TYPE (complex_8) cmulsc
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cmulsc%r=s1%r*s2
       cmulsc%kind=1
    case(m2)
       localmaster=master
       call ass(cmulsc)
       cmulsc%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cmulsc)
          call varck1(s1)
          cmulsc%t= varc1*s2
          master=localmaster
       else
          cmulsc%r=s1%r*s2
          cmulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cmulsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cmulsc

  FUNCTION cscmul( S2, S1 )
    implicit none
    TYPE (complex_8) cscmul
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cscmul%r=s1%r*s2
       cscmul%kind=1
    case(m2)
       localmaster=master
       call ass(cscmul)
       cscmul%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cscmul)
          call varck1(s1)
          cscmul%t= varc1*s2
          master=localmaster
       else
          cscmul%r=s1%r*s2
          cscmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cscmul"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cscmul

  FUNCTION cpmulsc( S1, S2 )
    implicit none
    TYPE (complex_8) cpmulsc
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpmulsc%r=s1%r*s2
       cpmulsc%kind=1
    case(m2)
       localmaster=master
       call ass(cpmulsc)
       cpmulsc%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpmulsc)
          call varfk1(s1)
          cpmulsc%t= varf1*s2
          master=localmaster
       else
          cpmulsc%r=s1%r*s2
          cpmulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpmulsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpmulsc

  FUNCTION cpscmul(S2 , S1 )
    implicit none
    TYPE (complex_8) cpscmul
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpscmul%r=s1%r*s2
       cpscmul%kind=1
    case(m2)
       localmaster=master
       call ass(cpscmul)
       cpscmul%t= s1%t*s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpscmul)
          call varfk1(s1)
          cpscmul%t= varf1*s2
          master=localmaster
       else
          cpscmul%r=s1%r*s2
          cpscmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpscmul"
         write(6,*) "s1%kind ",s1%kind
         read(5,*) localmaster
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpscmul

  FUNCTION cpaddsc( S1, S2 )
    implicit none
    TYPE (complex_8) cpaddsc
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpaddsc%r=s1%r+s2
       cpaddsc%kind=1
    case(m2)
       localmaster=master
       call ass(cpaddsc)
       cpaddsc%t= s1%t+s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpaddsc)
          call varfk1(s1)
          cpaddsc%t= varf1+s2
          master=localmaster
       else
          cpaddsc%r=s1%r+s2
          cpaddsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpaddsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpaddsc

  FUNCTION cpscadd( S2  ,S1)
    implicit none
    TYPE (complex_8) cpscadd
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpscadd%r=s1%r+s2
       cpscadd%kind=1
    case(m2)
       localmaster=master
       call ass(cpscadd)
       cpscadd%t= s1%t+s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpscadd)
          call varfk1(s1)
          cpscadd%t= varf1+s2
          master=localmaster
       else
          cpscadd%r=s1%r+s2
          cpscadd%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpscadd"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpscadd

  FUNCTION cpsubsc( S1, S2 )
    implicit none
    TYPE (complex_8) cpsubsc
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpsubsc%r=s1%r-s2
       cpsubsc%kind=1
    case(m2)
       localmaster=master
       call ass(cpsubsc)
       cpsubsc%t= s1%t-s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpsubsc)
          call varfk1(s1)
          cpsubsc%t= varf1-s2
          master=localmaster
       else
          cpsubsc%r=s1%r-s2
          cpsubsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpsubsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpsubsc

  FUNCTION cpscsub( S2  ,S1)
    implicit none
    TYPE (complex_8) cpscsub
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpscsub%r=s2-s1%r
       cpscsub%kind=1
    case(m2)
       localmaster=master
       call ass(cpscsub)
       cpscsub%t= s2-s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpscsub)
          call varfk1(s1)
          cpscsub%t= s2-varf1
          master=localmaster
       else
          cpscsub%r=s2-s1%r
          cpscsub%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpscsub"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpscsub

  FUNCTION cpdivsc( S1, S2 )
    implicit none
    TYPE (complex_8) cpdivsc
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpdivsc%r=s1%r/s2
       cpdivsc%kind=1
    case(m2)
       localmaster=master
       call ass(cpdivsc)
       cpdivsc%t= s1%t/s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpdivsc)
          call varfk1(s1)
          cpdivsc%t= varf1/s2
          master=localmaster
       else
          cpdivsc%r=s1%r/s2
          cpdivsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpdivsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpdivsc

  FUNCTION cpscdiv( S2  ,S1)
    implicit none
    TYPE (complex_8) cpscdiv
    TYPE (real_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cpscdiv%r=s2/s1%r
       cpscdiv%kind=1
    case(m2)
       localmaster=master
       call ass(cpscdiv)
       cpscdiv%t= s2/s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cpscdiv)
          call varfk1(s1)
          cpscdiv%t= s2/varf1
          master=localmaster
       else
          cpscdiv%r=s2/s1%r
          cpscdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cpscdiv"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cpscdiv

  FUNCTION dmulsc( S1, S2 )
    implicit none
    TYPE (complex_8) dmulsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dmulsc%t= varc1*s2
          master=localmaster
       else
          dmulsc%r=s1%r*s2
          dmulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dmulsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dmulsc

  FUNCTION dscmul( S2, S1 )
    implicit none
    TYPE (complex_8) dscmul
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dscmul%t= varc1*s2
          master=localmaster
       else
          dscmul%r=s1%r*s2
          dscmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dscmul"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dscmul

  FUNCTION mulsc( S1, S2 )
    implicit none
    TYPE (complex_8) mulsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          mulsc%t= varc1*REAL(s2,kind=DP)
          master=localmaster
       else
          mulsc%r=s1%r*REAL(s2,kind=DP)
          mulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in mulsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION mulsc

  FUNCTION scmul( S2, S1 )
    implicit none
    TYPE (complex_8) scmul
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          scmul%t= varc1*REAL(s2,kind=DP)
          master=localmaster
       else
          scmul%r=s1%r*REAL(s2,kind=DP)
          scmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in scmul"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION scmul

  FUNCTION imulsc( S1, S2 )
    implicit none
    TYPE (complex_8) imulsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          imulsc%t= varc1*REAL(s2,kind=DP)
          master=localmaster
       else
          imulsc%r=s1%r*REAL(s2,kind=DP)
          imulsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in imulsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION imulsc

  FUNCTION iscmul( S2, S1 )
    implicit none
    TYPE (complex_8) iscmul
    TYPE (complex_8), INTENT (IN) :: S1
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
        iscmul%     kind=1
       endif
    case(m3)
       if(knob) then
          localmaster=master
          call ass(iscmul)
          call varck1(s1)
          iscmul%t= varc1*REAL(s2,kind=DP)
          master=localmaster
       else
          iscmul%r=s1%r*REAL(s2,kind=DP)
          iscmul%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in iscmul"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION iscmul

  FUNCTION div( S1, S2 )
    implicit none
    TYPE (complex_8) div
    TYPE (complex_8), INTENT (IN) :: S1, S2
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
             call varck1(s1)
             div%t= varc1/s2%r
             master=localmaster
          else
             div%r=s1%r/s2%r
             div%kind=1
          endif
       case(m13)
          if(knob) then
             localmaster=master
             call ass(div)
             call varck2(s2)
             div%t= s1%r/varc2
             master=localmaster
          else
             div%r=s1%r/s2%r
             div%kind=1
          endif
       case(m32)
          localmaster=master
          call ass(div)
          if(knob) then
             call varck1(s1)
             div%t= varc1/s2%t
          else
             div%t=s1%r/s2%t
          endif
          master=localmaster
       case(m23)
          localmaster=master
          call ass(div)
          if(knob) then
             call varck2(s2)
             div%t= s1%t/varc2
          else
             div%t= s1%t/s2%r
          endif
          master=localmaster
       case(m33)
          if(knob) then
             localmaster=master
             call ass(div)
             call varck1(s1)
             call varck2(s2)
             div%t= varc1/varc2
             master=localmaster
          else
             div%r=s1%r/s2%r
             div%kind=1
          endif
       end select
    case default
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(2((1X,A72)))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in div "
         write(6,*) "s1%kind ,s2%kind "
       !w_p=(/s1%kind ,s2%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION div

  FUNCTION ddivsc( S1, S2 )
    implicit none
    TYPE (complex_8) ddivsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          ddivsc%t= varc1/s2
          master=localmaster
       else
          ddivsc%r=s1%r/s2
          ddivsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in ddivsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION ddivsc

  FUNCTION dscdiv( S2, S1 )
    implicit none
    TYPE (complex_8) dscdiv
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dscdiv%t= s2/varc1
          master=localmaster
       else
          dscdiv%r=s2/s1%r
          dscdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dscdiv"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dscdiv

  FUNCTION cdivsc( S1, S2 )
    implicit none
    TYPE (complex_8) cdivsc
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cdivsc%r=s1%r/s2
       cdivsc%kind=1
    case(m2)
       localmaster=master
       call ass(cdivsc)
       cdivsc%t= s1%t/s2
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cdivsc)
          call varck1(s1)
          cdivsc%t= varc1/s2
          master=localmaster
       else
          cdivsc%r=s1%r/s2
          cdivsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in cdivsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cdivsc

  FUNCTION cscdiv( S2, S1 )
    implicit none
    TYPE (complex_8) cscdiv
    TYPE (complex_8), INTENT (IN) :: S1
    complex(dp) , INTENT (IN) :: S2
    integer localmaster

    select case(s1%kind)
    case(m1)
       cscdiv%r=s2/s1%r
       cscdiv%kind=1
    case(m2)
       localmaster=master
       call ass(cscdiv)
       cscdiv%t= s2/s1%t
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(cscdiv)
          call varck1(s1)
          cscdiv%t= s2/varc1
          master=localmaster
       else
          cscdiv%r=s2/s1%r
          cscdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in "
         write(6,*) "s1%kind cscdiv"
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION cscdiv

  FUNCTION divsc( S1, S2 )
    implicit none
    TYPE (complex_8) divsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          divsc%t= varc1/REAL(s2,kind=DP)
          master=localmaster
       else
          divsc%r=s1%r/REAL(s2,kind=DP)
          divsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in "
         write(6,*) "s1%kind divsc"
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION divsc

  FUNCTION scdiv( S2, S1 )
    implicit none
    TYPE (complex_8) scdiv
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          scdiv%t= REAL(s2,kind=DP)/varc1
          master=localmaster
       else
          scdiv%r=REAL(s2,kind=DP)/s1%r
          scdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in "
         write(6,*) "s1%kind scdiv"
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION scdiv

  FUNCTION idivsc( S1, S2 )
    implicit none
    TYPE (complex_8) idivsc
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          idivsc%t= varc1/REAL(s2,kind=DP)
          master=localmaster
       else
          idivsc%r=s1%r/REAL(s2,kind=DP)
          idivsc%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in idivsc"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION idivsc

  FUNCTION iscdiv( S2, S1 )
    implicit none
    TYPE (complex_8) iscdiv
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          iscdiv%t= REAL(s2,kind=DP)/varc1
          master=localmaster
       else
          iscdiv%r=REAL(s2,kind=DP)/s1%r
          iscdiv%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in iscdiv"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION iscdiv

  FUNCTION POW( S1, S2 )
    implicit none
    TYPE (complex_8) POW
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          POW%t= varc1**s2
          master=localmaster
       else
          POW%r=s1%r**s2
          POW%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in POW"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION POW

  FUNCTION POWR( S1, S2 )
    implicit none
    TYPE (complex_8) POWR
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          POWR%t= varc1**REAL(s2,kind=DP)
          master=localmaster
       else
          POWR%r=s1%r**REAL(s2,kind=DP)
          POWR%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in POWR"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION POWR

  FUNCTION POWR8( S1, S2 )
    implicit none
    TYPE (complex_8) POWR8
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          POWR8%t= varc1**s2
          master=localmaster
       else
          POWR8%r=s1%r**s2
          POWR8%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in POWR8"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION POWR8

  FUNCTION dexpt( S1 )
    implicit none
    TYPE (complex_8) dexpt
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dexpt%t= exp(varc1)
          master=localmaster
       else
          dexpt%r=exp(s1%r)
          dexpt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dexpt"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dexpt

  FUNCTION CONJGT( S1 )
    implicit none
    TYPE (complex_8) CONJGT
    TYPE (complex_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       CONJGT%r=CONJG(s1%r)
       CONJGT%kind=1
    case(m2)
       localmaster=master
       call ass(CONJGT)
       CONJGT%t= CONJG(s1%t)
       master=localmaster
    case(m3)
       if(knob) then
          localmaster=master
          call ass(CONJGT)
          call varck1(s1)
          CONJGT%t= CONJG(varc1)
          master=localmaster
       else
          CONJGT%r=CONJG(s1%r)
          CONJGT%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in CONJGT"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION CONJGT

  FUNCTION dcost( S1 )
    implicit none
    TYPE (complex_8) dcost
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dcost%t= cos(varc1)
          master=localmaster
       else
          dcost%r=cos(s1%r)
          dcost%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dcost"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dcost

  FUNCTION dsint( S1 )
    implicit none
    TYPE (complex_8) dsint
    TYPE (complex_8), INTENT (IN) :: S1
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
          call varck1(s1)
          dsint%t= sin(varc1)
          master=localmaster
       else
          dsint%r=sin(s1%r)
          dsint%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in "
         write(6,*) "s1%kind dsint"
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dsint


  FUNCTION dlogt( S1 )
    implicit none
    TYPE (complex_8) dlogt
    TYPE (complex_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dlogt%r=log(s1%r)
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
          call varck1(s1)
          dlogt%t= log(varc1)
          master=localmaster
       else
          dlogt%r=log(s1%r)
          dlogt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dlogt"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dlogt

  FUNCTION dsqrtt( S1 )
    implicit none
    TYPE (complex_8) dsqrtt
    TYPE (complex_8), INTENT (IN) :: S1
    integer localmaster

    select case(s1%kind)
    case(m1)
       dsqrtt%r=sqrt(s1%r)
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
          call varck1(s1)
          dsqrtt%t= SQRT(varc1)
          master=localmaster
       else
          dsqrtt%r=sqrt(s1%r)
          dsqrtt%kind=1
       endif
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in dsqrtt"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION dsqrtt

  FUNCTION abst( S1 )
    implicit none
    real(dp) abst
    TYPE (complex_8), INTENT (IN) :: S1



    select case(s1%kind)
    case(m1,m3)
       abst=SQRT(REAL(s1%r,kind=DP)**2+aimag(s1%r)**2)
    case(m2)
       abst=SQRT(REAL(s1%t.sub.'0',kind=DP)**2+aimag(s1%t.sub.'0')**2)
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(1((1X,i4)))'
         write(6,*) " trouble in abst"
         write(6,*) "s1%kind "
       !w_p=(/s1%kind/)
       ! call !write_e(0)
    end select
  END FUNCTION abst


  subroutine asscp(s1)
    implicit none
    TYPE (complex_8) s1
    integer ipause,mypauses

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt+1)
       line=" cannot indent anymore "
       ipause=mypauses(0,line)
    end select

    !    write(26,*) "  complex polymorph ",master
    call ass0(s1%t%r)
    call ass0(s1%t%i)
    !call assc(s1%t)
    !call ass1(s1%t%i)
    !call ass0(s1%t%i)
    s1%alloc=t
    s1%kind=2
    s1%i=0
    s1%j=0


  end subroutine asscp

  !  subroutine asscp0
  !    implicit none
  !    integer ipause,mypauses!

  !    select case(master)
  !    case(0:ndumt-1)
  !       master=master+1
  !    case(ndumt+1)
  !       line=" cannot indent anymore "
  !       ipause=mypauses(0,line)
  !    end select

  !  end subroutine asscp0

  subroutine make_it_knobc(k,i,j,s)
    implicit none
    TYPE (complex_8), intent(inout) :: k
    complex(dp), optional :: s
    integer, intent(in) :: i,j
    if(i==0.and.j==0) return
    k%s=1.0_dp
    if(present(s)) k%s=s
    k%i=i
    k%j=j
    k%kind=3
  end subroutine make_it_knobc

  subroutine kill_knobc(k)
    implicit none
    TYPE (complex_8), intent(inout) :: k
    k%s=1.0_dp
    k%i=0
    k%j=0
    k%kind=1
  end subroutine kill_knobc

  SUBROUTINE  varck1(S2)
    implicit none
    type (complex_8)  S2

    if(knob) then
    !   if(nb_==0) then
          varc1=(/S2%R,S2%S/).var.(/s2%i+npara_fpp,s2%j+npara_fpp/)
     !  elseif(s2%nb==nb_) then
     !     varc1=(/S2%R,S2%S/).var.(/s2%i+npara_fpp-s2%g+1,s2%j+npara_fpp-s2%g+1/)
     !  else
      !    varc1=S2%R
     !  endif
    else ! Not a knob
       stop 3330   ! buggy never used
    !   varc1=(/S2%R,S2%S/).var.(/0,0/)
    endif

  end SUBROUTINE  varck1

  SUBROUTINE  varck2(S2)
    implicit none
    type (complex_8)  S2


    if(knob) then
   !    if(nb_==0) then
          varc2=(/S2%R,S2%S/).var.(/s2%i+npara_fpp,s2%j+npara_fpp/)
    !   elseif(s2%nb==nb_) then
    !      varc2=(/S2%R,S2%S/).var.(/s2%i+npara_fpp-s2%g+1,s2%j+npara_fpp-s2%g+1/)
    !   else
    !      varc2=S2%R
    !   endif
    else ! Not a knob
       stop 3331   ! buggy never used
    !   varc2=(/S2%R,S2%S/).var.(/0,0/)
    endif

  end SUBROUTINE  varck2

  ! remove small numbers

  SUBROUTINE  clean_complex_8(S1,S2,prec)
    implicit none
    type (complex_8),INTENT(INOUT)::S2
    type (complex_8), intent(INOUT):: s1
    real(dp) prec
 
    type(complex_8) t

    call alloc(t)
    t=s1

    select case(s1%kind)
    case(m1)
       if(abs(t%r)<prec) t%r=0.0_dp
    case(m2)
       call clean_complextaylor(t%t,t%t,prec)
    case(m3)
       Write(6,*) " cannot clean a knob "
       stop 601
    case default
       !w_p=0
       !w_p%nc=2
       !w_p%fc='((1X,A72,/,1x,a72))'
       !w_p%fi='(2((1X,i4)))'
         write(6,*) " trouble in clean_complex_8 "
         write(6,*) "s1%kind   "
       !w_p=(/s1%kind  /)
       ! call !write_e(0)
    end select
    s2=t
    call kill(t)


  END SUBROUTINE clean_complex_8



end module polymorphic_complextaylor
