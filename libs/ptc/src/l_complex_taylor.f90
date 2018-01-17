!The Full Polymorphic Package
!Copyright (C) Etienne Forest

module complex_taylor
  use tpsalie_analysis
  implicit none
  public
  private mul,cscmul,cmulsc,dscmul,dmulsc,mulsc,scmul,imulsc,iscmul
  private ctmul, cmult,ctadd,caddt,ctsub,csubt,ctdiv,cdivt
  private add,cscadd,dscadd,caddsc,daddsc,unaryADD,addsc,scadd,iaddsc,iscadd
  private tadd,addt,tmul,mult,tsub,subt,tdiv,divt
  private inv,div,dscdiv,cscdiv,cdivsc,ddivsc,divsc,scdiv,idivsc,iscdiv
  private subs,cscsub,dscsub,csubsc,dsubsc,iscsub,isubsc,scsub,subsc,unarySUB
  private EQUAL,cequaldacon,Dequaldacon,equaldacon,Iequaldacon,ctEQUAL,tcEQUAL
  private pow,powr,POWR8  !,DAABSEQUAL,AABSEQUAL 2002.10.17
  private alloccomplex,A_OPT,K_OPT   !,printcomplex ,killcomplex
  private logtpsat,exptpsat,abstpsat,dcost,dsint,datant,tant,dasint,dacost
  private dcosht,dsinht,dtanht,dsqrtt
  private getdiff,getdATRA,GETORDER,CUTORDER,getchar ,dputchar,dputint
  private set_in_complex   !, assc  !check,
  private dimagt,drealt,dcmplxt,CEQUAL,DEQUAL,REQUAL,CONJGT
  private GETCHARnd2,GETintnd2,GETint,getcharnd2s,GETintnd2s,GETintk
  private CFUC,CFURES,varco,varco1
  !  completing tpsa.f90
  private datantt,dasintt,dacostt,full_abstpsat
  integer,private::NO,ND,ND2,NP,NDPT,NV           !,lastmaster 2002.12.13
  logical(lp),private::old
  logical :: debug_flag =.false.
  logical :: debug_acos=.false.

  INTERFACE assignment (=)
     MODULE PROCEDURE EQUAL
     MODULE PROCEDURE ctEQUAL
     MODULE PROCEDURE tcEQUAL
     MODULE PROCEDURE CEQUAL
     MODULE PROCEDURE DEQUAL
     MODULE PROCEDURE REQUAL
     !     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
     !     MODULE PROCEDURE AABSEQUAL    ! remove 2002.10.17
     MODULE PROCEDURE cequaldacon
     MODULE PROCEDURE Dequaldacon
     MODULE PROCEDURE equaldacon
     MODULE PROCEDURE Iequaldacon
  end  INTERFACE

  !@     <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="684" height="445">
  !@       <tr>
  !@         <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; +</b>&nbsp;&nbsp;&nbsp; </font></span></td>
  !@         <td width="193" height="40" align="center" colspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="288" height="40" align="center" colspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="95" height="82" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="88" height="39" align="center">
  !@         <p>
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></p>
  !@         </td>
  !@         <td width="105" height="42" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(dp)</font></span></td>
  !@         <td width="88" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="101" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@         <td width="100" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="109" align="center" rowspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></td>
  !@         <td width="84" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#ADD" style="text-decoration: none; font-weight: 700">ADD</a></font></td>
  !@         <td width="105" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CADDSC" style="text-decoration: none; font-weight: 700">CADDSC</a></font></td>
  !@         <td width="92" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#ADDT" style="text-decoration: none; font-weight: 700">
  !@         ADDT</a></font></td>
  !@         <td width="102" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DADDSC" style="text-decoration: none; font-weight: 700">DADDSC</a></font></td>
  !@         <td width="99" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#ADDSC" style="text-decoration: none; font-weight: 700">ADDSC</a></font></td>
  !@         <td width="94" height="53" align="center"><font size="2">
  !@         <a href="l_complex_taylor.htm#IADDSC" style="text-decoration: none; font-weight: 700">
  !@         IADDSC</a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CSCADD" style="text-decoration: none; font-weight: 700">CSCADD</a></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="92" height="55" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CADDT" style="text-decoration: none; font-weight: 700">
  !@         CADDT</a></font></td>
  !@         <td width="102" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="99" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="171" align="center" rowspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="84" height="50" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#TADD" style="text-decoration: none; font-weight: 700">TADD</a></font></td>
  !@         <td width="105" height="50" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CTADD" style="text-decoration: none; font-weight: 700">CTADD</a></font></td>
  !@         <td width="92" height="50" align="center">
  !@            <span style="text-transform: uppercase; ">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#ADD" style="text-decoration: none; font-style:italic">
  !@         <font color="#FF0000">add</font></a></font></span></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase; ">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#DADDSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">daddsc</font></a></font></span></td>
  !@            <td width="78" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#ADDSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">ADDSC</font></a></font></td>
  !@            <td width="56" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <font color="#FF0000">&nbsp;&nbsp; </font>
  !@            <a href="i_tpsa.htm#IADDSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">IADDSC</font></a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DSCADD" style="text-decoration: none; font-weight: 700">DSCADD</a></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase; ">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#DSCADD" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">dscadd</font></a></font></span></td>
  !@         <td width="102" height="55" align="center"><b>F90</b></td>
  !@         <td width="99" height="55" align="center"><b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@         <td width="84" height="52" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#SCADD" style="text-decoration: none; font-weight: 700">
  !@         SCADD</a></font></td>
  !@         <td width="105" height="52" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#SCADD" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">SCADD</font></a></font></td>
  !@         <td width="102" height="52" align="center"><b>F90</b></td>
  !@         <td width="99" height="52" align="center"><b>F90</b></td>
  !@         <td width="94" height="52" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="78" height="56" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@         <td width="84" height="61" align="center"><font size="2">
  !@         <a href="l_complex_taylor.htm#ISCADD" style="text-decoration: none; font-weight: 700">ISCADD</a></font></td>
  !@         <td width="105" height="61" align="center"><b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#ISCADD" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">ISCADD</font></a></font></td>
  !@         <td width="102" height="61" align="center"><b>F90</b></td>
  !@         <td width="99" height="61" align="center"><b>F90</b></td>
  !@         <td width="94" height="61" align="center"><b>F90</b></td>
  !@       </tr>
  !@     </table>
  !@

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add
     MODULE PROCEDURE tadd
     MODULE PROCEDURE addt
     MODULE PROCEDURE cscadd
     MODULE PROCEDURE dscadd
     MODULE PROCEDURE ctadd
     MODULE PROCEDURE caddt
     MODULE PROCEDURE caddsc
     MODULE PROCEDURE daddsc
     MODULE PROCEDURE unaryADD
     MODULE PROCEDURE addsc
     MODULE PROCEDURE scadd
     MODULE PROCEDURE iaddsc
     MODULE PROCEDURE iscadd
  END INTERFACE

  !@     <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber2" width="684" height="445">
  !@       <tr>
  !@         <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  !@         -</b>&nbsp;&nbsp; </font></span></td>
  !@         <td width="193" height="40" align="center" colspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="288" height="40" align="center" colspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="95" height="82" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="88" height="39" align="center">
  !@         <p>
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></p>
  !@         </td>
  !@         <td width="105" height="42" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(dp)</font></span></td>
  !@         <td width="88" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="101" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@         <td width="100" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="109" align="center" rowspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></td>
  !@         <td width="84" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="l_complex_taylor.htm#SUBS">SUBS</a></font></td>
  !@         <td width="105" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="l_complex_taylor.htm#CSUBSC">CSUBSC</a></font></td>
  !@         <td width="92" height="53" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="l_complex_taylor.htm#SUBT">SUBT</a></font></td>
  !@         <td width="102" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DSUBSC" style="text-decoration: none; font-weight:700">DSUBSC</a></font></td>
  !@         <td width="99" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#SUBSC" style="text-decoration: none; font-weight:700">SUBSC</a></font></td>
  !@         <td width="94" height="53" align="center"><font size="2">
  !@         <a href="l_complex_taylor.htm#ISUBSC" style="text-decoration: none; font-weight:700">ISUBSC</a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@          <ahref="l_complex_taylor.htm#CSCSUB" style="text-decoration: none; font-weight:700">
  !@          <font size="2">
  !@         <a href="l_complex_taylor.htm#CSCSUB" style="text-decoration: none; font-weight: 700">CSCSUB</a></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="92" height="55" align="center">
  !@         <ahref="l_complex_taylor.htm#CSUBT" style="text-decoration: none">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CSUBT" style="text-decoration: none"><b>CSUBT</b></a></a></font></td>
  !@         <td width="102" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="99" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="171" align="center" rowspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="84" height="50" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#TSUB" style="text-decoration: none; font-weight:700">TSUB</a></font></td>
  !@         <td width="105" height="50" align="center">
  !@         <font size="2">
  !@         <a style="text-decoration: none; font-weight: 700" href="l_complex_taylor.htm#CTSUB">
  !@         CTSUB</a></font></td>
  !@         <td width="92" height="50" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#SUBS">
  !@         <font color="#FF0000">SUBS</font></a></font></span></td>
  !@            <td width="77" height="20" align="center">
  !@             <span style="text-transform: uppercase">
  !@             <font face="Times New Roman" size="2">
  !@            <a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#DSUBSC">
  !@            <font color="#FF0000">dSUBsc</font></a></font></span></td>
  !@            <td width="78" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#SUBSC">
  !@             <font color="#FF0000">SUBSC</font></a></font></td>
  !@            <td width="56" height="20" align="center">
  !@             <font size="2" face="Times New Roman">
  !@            <font color="#FF0000">&nbsp;&nbsp;
  !@            <a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#ISUBSC">&nbsp;</a></font><a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#ISUBSC"><font color="
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <ahref="l_complex_taylor.htm#DSCSUB" style="text-decoration: none; font-weight:700">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DSCSUB" style="text-decoration: none; font-weight: 700">DSCSUB</a></a></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#DSCSUB">
  !@            <font color="#FF0000">dscSUB</font></a></font></span></td>
  !@         <td width="102" height="55" align="center"><b>F90</b></td>
  !@         <td width="99" height="55" align="center"><b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@         <td width="84" height="52" align="center">
  !@         <ahref="l_complex_taylor.htm#SCSUB" style="text-decoration: none; font-weight:700">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#SCSUB" style="text-decoration: none; font-weight: 700">SCSUB</a></a></font></td>
  !@         <td width="105" height="52" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a style="text-decoration: none; font-style:italic" href="i_tpsa.htm#SCSUB">
  !@            <font color="#FF0000">SCSUB</font></a></font></td>
  !@         <td width="102" height="52" align="center"><b>F90</b></td>
  !@         <td width="99" height="52" align="center"><b>F90</b></td>
  !@         <td width="94" height="52" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="78" height="56" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@         <td width="84" height="61" align="center">
  !@         <ahref="l_complex_taylor.htm#ISCSUB" style="text-decoration: none; font-weight:700"><font size="2">
  !@          <a href="l_complex_taylor.htm#ISCSUB" style="text-decoration: none; font-weight: 700">ISCSUB</a></font></td>
  !@         <td width="105" height="61" align="center"><b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#ISCSUB" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">ISCSUB</font></a></font></td>
  !@         <td width="102" height="61" align="center"><b>F90</b></td>
  !@         <td width="99" height="61" align="center"><b>F90</b></td>
  !@         <td width="94" height="61" align="center"><b>F90</b></td>
  !@       </tr>
  !@     </table>


  INTERFACE OPERATOR (-)
     MODULE PROCEDURE unarySUB
     MODULE PROCEDURE subs
     MODULE PROCEDURE ctsub
     MODULE PROCEDURE csubt
     MODULE PROCEDURE tsub
     MODULE PROCEDURE subt
     MODULE PROCEDURE cscsub
     MODULE PROCEDURE dscsub
     MODULE PROCEDURE csubsc
     MODULE PROCEDURE dsubsc
     MODULE PROCEDURE subsc
     MODULE PROCEDURE scsub
     MODULE PROCEDURE isubsc
     MODULE PROCEDURE iscsub
  END INTERFACE

  !@     <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber3" width="684" height="445">
  !@       <tr>
  !@         <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </b>
  !@         *&nbsp;&nbsp; </font></span></td>
  !@         <td width="193" height="40" align="center" colspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="288" height="40" align="center" colspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="95" height="82" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="88" height="39" align="center">
  !@         <p>
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></p>
  !@         </td>
  !@         <td width="105" height="42" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(dp)</font></span></td>
  !@         <td width="88" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="101" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@         <td width="100" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="109" align="center" rowspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></td>
  !@         <td width="84" height="53" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#MUL" style="text-decoration: none">MUL</a></font></b></td>
  !@         <td width="105" height="53" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#CMULSC" style="text-decoration: none">CMULSC</a></font></b></td>
  !@         <td width="92" height="53" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#MULT" style="text-decoration: none">MULT</a></font></b></td>
  !@         <td width="102" height="53" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#DMULSC" style="text-decoration: none">DMULSC</a></font></b></td>
  !@         <td width="99" height="53" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#MULSC" style="text-decoration: none">MULSC</a></font></b></td>
  !@         <td width="94" height="53" align="center"><b><font size="2">
  !@         <a href="l_complex_taylor.htm#IMULSC" style="text-decoration: none">IMULSC</a></font></b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#CSCMUL" style="text-decoration: none">CSCMUL</a></font></b></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="92" height="55" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#CMULT" style="text-decoration: none">CMULT</a></font></b></td>
  !@         <td width="102" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="99" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="171" align="center" rowspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="84" height="50" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#TMUL" style="text-decoration: none">TMUL</a></font></b></td>
  !@         <td width="105" height="50" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#CTMUL" style="text-decoration: none">CTMUL</a></font></b></td>
  !@         <td width="92" height="50" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#MUL" style="text-decoration: none; font-style:italic">
  !@         <font color="#FF0000">MUL</font></a></font></span></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#DMULSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">dMULsc</font></a></font></span></td>
  !@            <td width="78" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#MULSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">MULSC</font></a></font></td>
  !@            <td width="56" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <font color="#FF0000">&nbsp;&nbsp; </font>
  !@            <a href="i_tpsa.htm#IMULSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">IMULSC</font></a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#DSCMUL" style="text-decoration: none">DSCMUL</a></font></b></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#DSCMUL" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">dscMUL</font></a></font></span></td>
  !@         <td width="102" height="55" align="center"><b>F90</b></td>
  !@         <td width="99" height="55" align="center"><b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@         <td width="84" height="52" align="center">
  !@         <b>
  !@         <font size="2"><a href="l_complex_taylor.htm#SCMUL" style="text-decoration: none">SCMUL</a></font></b></td>
  !@         <td width="105" height="52" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#SCMUL" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">SCMUL</font></a></font></td>
  !@         <td width="102" height="52" align="center"><b>F90</b></td>
  !@         <td width="99" height="52" align="center"><b>F90</b></td>
  !@         <td width="94" height="52" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="78" height="56" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@         <td width="84" height="61" align="center"><b><font size="2">
  !@         <a href="l_complex_taylor.htm#ISCMUL" style="text-decoration: none">ISCMUL</a></font></b></td>
  !@         <td width="105" height="61" align="center"><b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2" face="Times New Roman">
  !@            <a href="i_tpsa.htm#ISCMUL" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">ISCMUL</font></a></font></td>
  !@         <td width="102" height="61" align="center"><b>F90</b></td>
  !@         <td width="99" height="61" align="center"><b>F90</b></td>
  !@         <td width="94" height="61" align="center"><b>F90</b></td>
  !@      </tr>
  !@    </table>


  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul
     MODULE PROCEDURE tmul
     MODULE PROCEDURE mult
     MODULE PROCEDURE  cscmul
     MODULE PROCEDURE  ctmul
     MODULE PROCEDURE  dscmul
     MODULE PROCEDURE  cmulsc
     MODULE PROCEDURE  cmult
     MODULE PROCEDURE  dmulsc
     MODULE PROCEDURE mulsc
     MODULE PROCEDURE scmul
     MODULE PROCEDURE imulsc
     MODULE PROCEDURE iscmul
  END INTERFACE

  !@     <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber4" width="684" height="445">
  !@       <tr>
  !@         <td width="78" height="84" align="center" rowspan="2" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="4"><b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; </b>
  !@         /&nbsp; </font></span></td>
  !@         <td width="193" height="40" align="center" colspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="288" height="40" align="center" colspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="95" height="82" align="center" rowspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="88" height="39" align="center">
  !@         <p>
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></p>
  !@         </td>
  !@         <td width="105" height="42" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(dp)</font></span></td>
  !@         <td width="88" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="101" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">rEAL(DP)</font></span></td>
  !@         <td width="100" height="39" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="109" align="center" rowspan="2">
  !@         <font size="2">COMPLEX</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX TAYLOR</font></span></td>
  !@         <td width="84" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DIV" style="text-decoration: none; font-weight:700">DIV</a></font></td>
  !@         <td width="105" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CDIVSC" style="text-decoration: none; font-weight:700">CDIVSC</a></font></td>
  !@         <td width="92" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DIVT" style="text-decoration: none; font-weight:700">DIVT</a></font></td>
  !@         <td width="102" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DDIVSC" style="text-decoration: none; font-weight:700">DDIVSC</a></font></td>
  !@         <td width="99" height="53" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DIVSC" style="text-decoration: none; font-weight:700">DIVSC</a></font></td>
  !@         <td width="94" height="53" align="center"><font size="2">
  !@         <a href="l_complex_taylor.htm#IDIVSC" style="text-decoration: none; font-weight:700">IDIVSC</a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="55" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">COMPLEX(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@          <font size="2">
  !@          <a href="l_complex_taylor.htm#CSCDIV" style="text-decoration: none; font-weight:700">CSCDIV</a></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="92" height="55" align="center">
  !@          <font size="2">
  !@          <a href="l_complex_taylor.htm#CDIVT" style="text-decoration: none; font-weight:700">CDIVT</a></font></td>
  !@         <td width="102" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="99" height="55" align="center">
  !@         <b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="171" align="center" rowspan="3">
  !@         <font size="2">REAL</font></td>
  !@         <td width="39" height="54" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">TAYLOR</font></span></td>
  !@         <td width="84" height="50" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#TDIV" style="text-decoration: none; font-weight:700">TDIV</a></font></td>
  !@         <td width="105" height="50" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#CTDIV" style="text-decoration: none; font-weight:700">CTDIV</a></font></td>
  !@         <td width="92" height="50" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#DIV" style="text-decoration: none; font-style:italic">
  !@         <font color="#FF0000">div</font></a></font></span></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase">
  !@             <font face="Times New Roman" size="2">
  !@             <a href="i_tpsa.htm#DDIVSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">dDIVsc</font></a></font></span></td>
  !@            <td width="78" height="20" align="center">
  !@            <font size="2">
  !@            <a href="i_tpsa.htm#DIVSC" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">DIVSC</font></a></font></td>
  !@            <td width="56" height="20" align="center">
  !@            <font size="2">
  !@             <font color="#FF0000">&nbsp;&nbsp;&nbsp;&nbsp; </font>
  !@             <a href="i_tpsa.htm#IDIVSC" style="text-decoration: none; font-style:italic">
  !@             <font color="#FF0000">IDIVSC</font></a></font></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(DP)</font></span></td>
  !@         <td width="84" height="55" align="center">
  !@         <font size="2">
  !@         <a href="l_complex_taylor.htm#DSCDIV" style="text-decoration: none; font-weight:700">DSCDIV</a></font></td>
  !@         <td width="105" height="55" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <span style="text-transform: uppercase">
  !@            <font face="Times New Roman" size="2">
  !@            <a href="i_tpsa.htm#DSCDIV" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">dscDIV</font></a></font></span></td>
  !@         <td width="102" height="55" align="center"><b>F90</b></td>
  !@         <td width="99" height="55" align="center"><b>F90</b></td>
  !@         <td width="94" height="55" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="39" height="56" align="center">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">REAL(SP)</font></span></td>
  !@         <td width="84" height="52" align="center">
  !@          <font size="2">
  !@          <a href="l_complex_taylor.htm#SCDIV" style="text-decoration: none; font-weight:700">SCDIV</a></font></td>
  !@         <td width="105" height="52" align="center">
  !@         <b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2">
  !@             <a href="i_tpsa.htm#SCDIV" style="text-decoration: none; font-style:italic">
  !@            <font color="#FF0000">SCDIV</font></a></font></td>
  !@         <td width="102" height="52" align="center"><b>F90</b></td>
  !@         <td width="99" height="52" align="center"><b>F90</b></td>
  !@         <td width="94" height="52" align="center"><b>F90</b></td>
  !@       </tr>
  !@       <tr>
  !@         <td width="78" height="56" align="center" colspan="2">
  !@         <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="2">Integer</font></span></td>
  !@         <td width="84" height="61" align="center"><font size="2">
  !@         <a href="l_complex_taylor.htm#ISCDIV" style="text-decoration: none; font-weight:700">
  !@         ISCDIV</a></font></td>
  !@         <td width="105" height="61" align="center"><b>F90</b></td>
  !@            <td width="77" height="20" align="center">
  !@            <font size="2">
  !@             <a href="i_tpsa.htm#ISCDIV" style="text-decoration: none; font-style:italic">
  !@             <font color="#FF0000">ISCDIV</font></a></font></td>
  !@         <td width="102" height="61" align="center"><b>F90</b></td>
  !@         <td width="99" height="61" align="center"><b>F90</b></td>
  !@         <td width="94" height="61" align="center"><b>F90</b></td>
  !@       </tr>
  !@     </table>


  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div
     MODULE PROCEDURE divt
     MODULE PROCEDURE tdiv
     MODULE PROCEDURE ctdiv
     MODULE PROCEDURE cdivt
     MODULE PROCEDURE ddivsc
     MODULE PROCEDURE cdivsc
     MODULE PROCEDURE dscdiv
     MODULE PROCEDURE cscdiv
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

  INTERFACE OPERATOR (.var.)
     MODULE PROCEDURE varco
     MODULE PROCEDURE varco1
  END INTERFACE

  INTERFACE OPERATOR (.mono.)
     MODULE PROCEDURE dputint     !@1 Accepts J(nv) </br>
     MODULE PROCEDURE dputchar  !@1 Accepts String such as '12 </br>
  END INTERFACE

  INTERFACE OPERATOR (.d.)
     MODULE PROCEDURE getdiff
  END INTERFACE

  INTERFACE OPERATOR (.K.)
     MODULE PROCEDURE getdATRA
  END INTERFACE

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

  INTERFACE OPERATOR (<=)
     MODULE PROCEDURE getcharnd2s
     MODULE PROCEDURE GETintnd2s
     MODULE PROCEDURE GETintk
  END INTERFACE

  ! Intrinsic Routines

  INTERFACE aimag
     MODULE PROCEDURE dimagt
  END INTERFACE
  INTERFACE dimag
     MODULE PROCEDURE dimagt
  END INTERFACE

  INTERFACE dble
     MODULE PROCEDURE drealt
  END INTERFACE
  INTERFACE dreal
     MODULE PROCEDURE drealt
  END INTERFACE

  INTERFACE cmplx
     MODULE PROCEDURE dcmplxt
  END INTERFACE
  INTERFACE dcmplx
     MODULE PROCEDURE dcmplxt
  END INTERFACE

  INTERFACE CONJG
     MODULE PROCEDURE CONJGT
  END INTERFACE

  INTERFACE abs
     MODULE PROCEDURE abstpsat
  END INTERFACE
  INTERFACE dabs
     MODULE PROCEDURE abstpsat
  END INTERFACE

  INTERFACE log
     MODULE PROCEDURE logtpsat
  END INTERFACE
  INTERFACE dlog
     MODULE PROCEDURE logtpsat
  END INTERFACE
  INTERFACE clog
     MODULE PROCEDURE logtpsat
  END INTERFACE
  INTERFACE cdlog
     MODULE PROCEDURE logtpsat
  END INTERFACE


  INTERFACE atan
     MODULE PROCEDURE datant
     MODULE PROCEDURE datantt
  END INTERFACE
  INTERFACE datan
     MODULE PROCEDURE datant
     MODULE PROCEDURE datantt
  END INTERFACE

  INTERFACE asin
     MODULE PROCEDURE dasint
     MODULE PROCEDURE dasintt
  END INTERFACE
  INTERFACE dasin
     MODULE PROCEDURE dasint
     MODULE PROCEDURE dasintt
  END INTERFACE

  INTERFACE acos
     MODULE PROCEDURE dacost
     MODULE PROCEDURE dacostt
  END INTERFACE
  INTERFACE dacos
     MODULE PROCEDURE dacost
     MODULE PROCEDURE dacostt
  END INTERFACE

  INTERFACE tan
     MODULE PROCEDURE tant
  END INTERFACE
  INTERFACE dtan
     MODULE PROCEDURE tant
  END INTERFACE



  INTERFACE cos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE cdcos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE ccos
     MODULE PROCEDURE dcost
  END INTERFACE
  INTERFACE dcos
     MODULE PROCEDURE dcost
  END INTERFACE


  INTERFACE sin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE cdsin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE csin
     MODULE PROCEDURE dsint
  END INTERFACE
  INTERFACE dsin
     MODULE PROCEDURE dsint
  END INTERFACE


  INTERFACE exp
     MODULE PROCEDURE exptpsat
  END INTERFACE
  INTERFACE dexp
     MODULE PROCEDURE exptpsat
  END INTERFACE
  INTERFACE cexp
     MODULE PROCEDURE exptpsat
  END INTERFACE
  INTERFACE cdexp
     MODULE PROCEDURE exptpsat
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

  INTERFACE sqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE
  INTERFACE dsqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE
  INTERFACE cdsqrt
     MODULE PROCEDURE dsqrtt
  END INTERFACE

  ! End Intrinsic Routines

  ! Non-intrisic Functions

  !  INTERFACE var
  !     MODULE PROCEDURE varc
  !     MODULE PROCEDURE varcC
  !  END INTERFACE
  !
  !  INTERFACE shiftda
  !     MODULE PROCEDURE shiftc
  !  END INTERFACE

  INTERFACE pok
     MODULE PROCEDURE pokc
  END INTERFACE

  INTERFACE pek
     MODULE PROCEDURE pekc
  END INTERFACE

  INTERFACE CFU
     MODULE PROCEDURE CFUC
     MODULE PROCEDURE CFURES
  END INTERFACE

  INTERFACE full_abs
     MODULE PROCEDURE full_abstpsat
  END INTERFACE


  ! i/o

  INTERFACE daprint
     MODULE PROCEDURE printcomplex
  END INTERFACE

  INTERFACE print
     MODULE PROCEDURE printcomplex
  END INTERFACE

  INTERFACE read
     MODULE PROCEDURE inputcomplex
  END INTERFACE

  INTERFACE dainput
     MODULE PROCEDURE inputcomplex
  END INTERFACE

  ! end of /o

  ! Constructors  and Destructors

  INTERFACE alloc
     MODULE PROCEDURE alloccomplex
     MODULE PROCEDURE a_opt
     MODULE PROCEDURE alloccomplexn
  END INTERFACE


  INTERFACE kill
     MODULE PROCEDURE killcomplex
     MODULE PROCEDURE k_opt
     MODULE PROCEDURE killcomplexn
  END INTERFACE

  ! end Constructors  and Destructors

  ! managing

  INTERFACE ass
     MODULE PROCEDURE assc
  END INTERFACE

  ! end managing


contains

  FUNCTION dimagt( S1 )
    implicit none
    TYPE (taylor) dimagt
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dimagt)    !2002.12.25

    dimagt=s1%i

    master=localmaster
  END FUNCTION dimagt

  FUNCTION drealt( S1 )
    implicit none
    TYPE (taylor) drealt
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(drealt)    !2002.12.25

    drealt=s1%r

    master=localmaster
  END FUNCTION drealt

  FUNCTION GETCHARnd2( S1, S2 )
    implicit none
    TYPE (complextaylor) GETCHARnd2
    TYPE (complextaylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2

    integer localmaster
    localmaster=master

    call ass(GETCHARnd2)


    GETCHARnd2%r=s1%r.par.s2
    GETCHARnd2%i=s1%i.par.s2

    master=localmaster


  END FUNCTION GETCHARnd2

  FUNCTION GETintnd2( S1, S2 )
    implicit none
    TYPE (complextaylor) GETintnd2
    TYPE (complextaylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)

    integer localmaster
    localmaster=master

    call ass(GETintnd2)


    GETintnd2%r=s1%r.par.s2
    GETintnd2%i=s1%i.par.s2

    master=localmaster


  END FUNCTION GETintnd2


  FUNCTION GETCHARnd2s( S1, S2 )
    implicit none
    TYPE (complextaylor) GETCHARnd2s
    TYPE (complextaylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2

    integer localmaster
    localmaster=master

    call ass(GETCHARnd2s)


    GETCHARnd2s%r=s1%r<=s2
    GETCHARnd2s%i=s1%i<=s2

    master=localmaster


  END FUNCTION GETCHARnd2s

  FUNCTION GETintnd2s( S1, S2 )
    implicit none
    TYPE (complextaylor) GETintnd2s
    TYPE (complextaylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)

    integer localmaster
    localmaster=master

    call ass(GETintnd2s)


    GETintnd2s%r= s1%r<=s2
    GETintnd2s%i= s1%i<=s2

    master=localmaster


  END FUNCTION GETintnd2s

  FUNCTION dputchar( S1, S2 )
    implicit none
    TYPE (complextaylor) dputchar
    complex(dp) , INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2
    !    CHARACTER (LEN = LNV)  resul
    !    integer j(lnv),i,nd2par

    integer localmaster
    localmaster=master

    call ass(dputchar)

    !    resul = trim(ADJUSTL (s2))
    !
    !    do i=1,lnv
    !       j(i)=0
    !    enddo
    !
    !    nd2par= len(trim(ADJUSTL (s2)))
    !    !frs    do i=1,len(trim(ADJUSTL (s2)))
    !    do i=1,nd2par
    !       CALL  CHARINT(RESUL(I:I),J(I))
    !       if(i>nv) then
    !          if(j(i)>0) then
    !             call var(dputchar,cmplx(zero,zero,kind=dp),0,0)
    !             return
    !          endif
    !       endif
    !   enddo
    !
    !
    !    call var(dputchar,cmplx(zero,zero,kind=dp),0,0)
    !    call pok(dputchar,j,s1)
    !
    dputchar%r= real(S1,kind=dp).mono.S2
    dputchar%i= aimag(S1).mono.S2


    master=localmaster

  END FUNCTION dputchar

  FUNCTION dputint( S1, S2 )
    implicit none
    TYPE (complextaylor) dputint
    complex(dp) , INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    !   integer j(lnv),i,nd2par

    integer localmaster
    localmaster=master

    call ass(dputint)

    !    do i=1,lnv
    !       j(i)=0
    !    enddo!
    !
    !    nd2par=size(s2)
    !    do i=1,nd2par
    !       J(I)=s2(i)
    !    enddo

    !frs    do i=1,len(trim(ADJUSTL (s2)))
    !    do i=1,nd2par
    !       if(i>nv) then
    !          if(j(i)>0) then
    !             call var(dputint,cmplx(zero,zero,kind=dp),0,0)
    !             return
    !          endif
    !!       endif
    !    enddo
    !


    !   call var(dputint,cmplx(zero,zero,kind=dp),0,0)
    !    call pok(dputint,j,s1)


    dputint%r= real(S1,kind=dp).mono.S2
    dputint%i= aimag(S1).mono.S2


    master=localmaster

  END FUNCTION dputint



  FUNCTION varco(s1,s2)
    implicit none
    TYPE (complextaylor) varco
    complex(dp) , INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(2)

    integer localmaster
    localmaster=master

    call ass(varco)

    varco%r=REAL(s1,kind=DP) + (1.0_dp.mono.s2(1))
    varco%i=aimag(s1) + (1.0_dp.mono.s2(2))


    !varco%r=REAL(s1,kind=DP).var.s2(1)
    !varco%i=aimag(s1).var.s2(2)

    master=localmaster

  END FUNCTION varco


  FUNCTION varco1(s1,s2)
    implicit none
    TYPE (complextaylor) varco1
    complex(dp) , INTENT (IN) :: S1(2)
    integer  , INTENT (IN) ::  S2(2)

    integer localmaster
    localmaster=master

    call ass(varco1)



    varco1=s1(1)+s1(2)*((1.0_dp.mono.s2(1))+i_*(1.0_dp.mono.s2(2)))

    master=localmaster

  END FUNCTION varco1




  FUNCTION GETORDER( S1, S2 )
    implicit none
    TYPE (complextaylor) GETORDER
    TYPE (complextaylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    integer localmaster
    localmaster=master
    call ass(GETORDER)

    GETORDER%r=S1%r.sub.s2
    GETORDER%i=S1%i.sub.s2


    master=localmaster

  END FUNCTION GETORDER

  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (complextaylor) CUTORDER
    TYPE (complextaylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    integer localmaster
    localmaster=master
    call ass(CUTORDER)

    CUTORDER%r=S1%r.CUT.s2
    CUTORDER%i=S1%i.CUT.s2


    master=localmaster

  END FUNCTION CUTORDER

  FUNCTION GETchar( S1, S2 )
    implicit none
    complex(dp) GETchar
    real(dp) r1,r2
    TYPE (complextaylor), INTENT (IN) :: S1
    CHARACTER(*)  , INTENT (IN) ::  S2

    integer localmaster
    localmaster=master

    r1=S1%r.sub.s2
    r2=S1%i.sub.s2

    GETchar=cmplx(r1,r2,kind=dp)

    master=localmaster
  END FUNCTION GETchar

  FUNCTION GETint( S1, S2 )   ! 2002.12.20
    implicit none
    complex(dp) GETint
    TYPE (complextaylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2(:)
    real(dp) r1,r2


    r1=S1%r.sub.s2
    r2=S1%i.sub.s2

    GETint=cmplx(r1,r2,kind=dp)


  END FUNCTION GETint


  FUNCTION POW( S1, R2 )
    implicit none
    TYPE (complextaylor) POW,temp
    TYPE (complextaylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) :: R2
    INTEGER I,R22

    integer localmaster
    localmaster=master
    call ass(pow)

    call alloc(temp)

    TEMP=1.0_dp


    R22=IABS(R2)
    DO I=1,R22
       temp=temp*s1
    ENDDO
    IF(R2.LT.0) THEN
       temp=1.0_dp/temp
    ENDIF

    POW=temp

    call kill(temp)
    master=localmaster


  END FUNCTION POW

  FUNCTION POWR( S1, R2 )
    implicit none
    TYPE (complextaylor) POWR,temp
    TYPE (complextaylor), INTENT (IN) :: S1
    REAL(SP), INTENT (IN) :: R2
    integer localmaster

    if(real_warning) call real_stop
    localmaster=master
    call ass(POWR)
    call alloc(temp)

    temp=log(s1)
    temp=temp*(r2)
    temp=exp(temp)
    POWR=temp

    call kill(temp)
    master=localmaster

  END FUNCTION POWR

  FUNCTION POWR8( S1, R2 )
    implicit none
    TYPE (complextaylor) POWR8,temp
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: R2

    integer localmaster
    localmaster=master
    call ass(powr8)
    call alloc(temp)

    temp=log(s1)
    temp=temp*r2
    temp=exp(temp)
    POWR8=temp

    call kill(temp)
    master=localmaster

  END FUNCTION POWR8

  FUNCTION GETdiff( S1, S2 )
    implicit none
    TYPE (complextaylor) GETdiff
    TYPE (complextaylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    integer localmaster
    localmaster=master
    call ass(GETdiff)

    getdiff%r=s1%r.d.s2
    getdiff%i=s1%i.d.s2

    master=localmaster
  END FUNCTION GETdiff

  FUNCTION GETdatra( S1, S2 )
    implicit none
    TYPE (complextaylor) GETdatra
    TYPE (complextaylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2

    integer localmaster
    localmaster=master
    call ass(GETdatra)

    GETdatra%r=s1%r.k.s2
    GETdatra%i=s1%i.k.s2

    master=localmaster
  END FUNCTION GETdatra

  SUBROUTINE  alloccomplex(S2)
    implicit none
    type (complextaylor),INTENT(INOUT)::S2
    call alloctpsa(s2%r)
    call alloctpsa(s2%i)
  END SUBROUTINE alloccomplex

  SUBROUTINE  alloccomplexn(S2,K)
    implicit none
    type (complextaylor),INTENT(INOUT),dimension(:)::S2
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
       call alloctpsa(s2(j)%r)
       call alloctpsa(s2(j)%i)
    enddo

  END SUBROUTINE alloccomplexn

  SUBROUTINE  A_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (complextaylor),INTENT(INout)::S1,S2
    type (complextaylor),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE A_opt

  SUBROUTINE  K_OPT(S1,S2,s3,s4,s5,s6,s7,s8,s9,s10)
    implicit none
    type (complextaylor),INTENT(INout)::S1,S2
    type (complextaylor),optional, INTENT(INout):: s3,s4,s5,s6,s7,s8,s9,s10
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
  END SUBROUTINE K_opt

  SUBROUTINE  printcomplex(S2,i,PREC)
    implicit none
    type (complextaylor),INTENT(INOUT)::S2
    integer,optional :: i
    REAL(DP),OPTIONAL,INTENT(INOUT)::PREC

    call daprint(s2%r,i,PREC)
    call daprint(s2%i,i,PREC)
  END SUBROUTINE printcomplex

  SUBROUTINE  inputcomplex(S2,i)
    implicit none
    type (complextaylor),INTENT(INOUT)::S2
    integer i
    call dainput(s2%r,i)
    call dainput(s2%i,i)
  END SUBROUTINE inputcomplex


  SUBROUTINE  killcomplex(S2)
    implicit none
    type (complextaylor),INTENT(INOUT)::S2
    call killTPSA(s2%r)
    call killTPSA(s2%i)
  END SUBROUTINE killcomplex

  SUBROUTINE  killcomplexn(S2,K)
    implicit none
    type (complextaylor),INTENT(INOUT),dimension(:)::S2
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
       call killtpsa(s2(j)%r)
       call killtpsa(s2(j)%i)
    enddo

  END SUBROUTINE killcomplexn

  FUNCTION mul( S1, S2 )
    implicit none
    TYPE (complextaylor) mul
    TYPE (complextaylor), INTENT (IN) :: S1, S2
    integer localmaster
    localmaster=master
    call ass(mul)
    mul%r=s1%r*s2%r-s1%i*s2%i
    mul%i=s1%r*s2%i+s1%i*s2%r
    master=localmaster
  END FUNCTION mul

  FUNCTION div( S1, S2 )
    implicit none
    TYPE (complextaylor) div ,t
    TYPE (complextaylor), INTENT (IN) :: S1, S2
    integer localmaster
    localmaster=master

    call ass(div)

    call alloc(t)
    call inv(s2,t)
    div%r=s1%r*t%r-s1%i*t%i
    div%i=s1%r*t%i+s1%i*t%r
    call kill(t)
    master=localmaster
  END FUNCTION div

  FUNCTION cscdiv( S1, S2 )
    implicit none
    TYPE (complextaylor) cscdiv ,t
    TYPE (complextaylor), INTENT (IN) ::  S2
    complex(dp), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(cscdiv)
    call alloc(t)
    call inv(s2,t)
    cscdiv%r=REAL(s1,kind=DP)*t%r-aimag(s1)*t%i
    cscdiv%i=REAL(s1,kind=DP)*t%i+aimag(s1)*t%r
    call kill(t)
    master=localmaster
  END FUNCTION cscdiv

  FUNCTION dscdiv( S1, S2 )
    implicit none
    TYPE (complextaylor) dscdiv ,t
    TYPE (complextaylor), INTENT (IN) ::  S2
    real(dp), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(dscdiv)
    call alloc(t)
    call inv(s2,t)
    dscdiv%r=s1*t%r
    dscdiv%i=s1*t%i
    call kill(t)
    master=localmaster
  END FUNCTION dscdiv

  FUNCTION scdiv( S1, S2 )
    implicit none
    TYPE (complextaylor) scdiv ,t
    TYPE (complextaylor), INTENT (IN) ::  S2
    real(sp), INTENT (IN) :: S1
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(scdiv)
    call alloc(t)
    call inv(s2,t)
    scdiv%r=s1*t%r
    scdiv%i=s1*t%i
    call kill(t)
    master=localmaster
  END FUNCTION scdiv

  FUNCTION iscdiv( S1, S2 )
    implicit none
    TYPE (complextaylor) iscdiv ,t
    TYPE (complextaylor), INTENT (IN) ::  S2
    integer, INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(iscdiv)
    call alloc(t)
    call inv(s2,t)
    iscdiv%r=s1*t%r
    iscdiv%i=s1*t%i
    call kill(t)
    master=localmaster
  END FUNCTION iscdiv

  FUNCTION idivsc( S2,S1  )
    implicit none
    TYPE (complextaylor) idivsc
    TYPE (complextaylor), INTENT (IN) ::  S2
    integer, INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(idivsc)
    idivsc%r=(1.0_dp/s1)*s2%r
    idivsc%i=(1.0_dp/s1)*s2%i
    master=localmaster
  END FUNCTION idivsc


  FUNCTION divsc( S2,S1  )
    implicit none
    TYPE (complextaylor) divsc
    TYPE (complextaylor), INTENT (IN) ::  S2
    real(sp), INTENT (IN) :: S1
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(divsc)
    divsc%r=(1.0_dp/s1)*s2%r
    divsc%i=(1.0_dp/s1)*s2%i
    master=localmaster
  END FUNCTION divsc


  FUNCTION ddivsc( S2,S1  )
    implicit none
    TYPE (complextaylor) ddivsc
    TYPE (complextaylor), INTENT (IN) ::  S2
    real(dp), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(ddivsc)
    ddivsc%r=(1.0_dp/s1)*s2%r
    ddivsc%i=(1.0_dp/s1)*s2%i
    master=localmaster
  END FUNCTION ddivsc

  FUNCTION cdivsc( S2,S1  )
    implicit none
    TYPE (complextaylor) cdivsc
    TYPE (complextaylor), INTENT (IN) ::  S2
    complex(dp), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(cdivsc)
    cdivsc%r=REAL((1.0_dp/s1),kind=DP)*s2%r-aimag((1.0_dp/s1))*s2%i
    cdivsc%i=REAL((1.0_dp/s1),kind=DP)*s2%i+aimag((1.0_dp/s1))*s2%r
    master=localmaster
  END FUNCTION cdivsc

  FUNCTION cscmul( sc,S1 )
    implicit none
    TYPE (complextaylor) cscmul
    TYPE (complextaylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(cscmul)

    cscmul%r=REAL(sc,kind=DP)*s1%r-aimag(sc)*s1%i
    cscmul%i=REAL(sc,kind=DP)*s1%i+aimag(sc)*s1%r
    master=localmaster
  END FUNCTION cscmul

  FUNCTION ctmul( S1,sc )
    implicit none
    TYPE (complextaylor) ctmul
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(ctmul)

    ctmul%r=REAL(sc,kind=DP)*s1
    ctmul%i=aimag(sc)*s1
    master=localmaster
  END FUNCTION ctmul

  FUNCTION cmult( sc,S1 )
    implicit none
    TYPE (complextaylor) cmult
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(cmult)

    cmult%r=REAL(sc,kind=DP)*s1
    cmult%i=aimag(sc)*s1
    master=localmaster
  END FUNCTION cmult

  FUNCTION caddt( sc,S1 )
    implicit none
    TYPE (complextaylor) caddt
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(caddt)

    caddt%r=REAL(sc,kind=DP)+s1
    caddt%i=aimag(sc)
    master=localmaster
  END FUNCTION caddt

  FUNCTION ctadd(S1, sc )
    implicit none
    TYPE (complextaylor) ctadd
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(ctadd)

    ctadd%r=REAL(sc,kind=DP)+s1
    ctadd%i=aimag(sc)
    master=localmaster
  END FUNCTION ctadd

  FUNCTION csubt( sc,S1 )
    implicit none
    TYPE (complextaylor) csubt
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(csubt)

    csubt%r=REAL(sc,kind=DP)-s1
    csubt%i=aimag(sc)
    master=localmaster
  END FUNCTION csubt

  FUNCTION ctsub( S1,sc )
    implicit none
    TYPE (complextaylor) ctsub
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(ctsub)

    ctsub%r=s1-REAL(sc,kind=DP)
    ctsub%i=-aimag(sc)
    master=localmaster
  END FUNCTION ctsub

  FUNCTION cdivt( sc,S1 )
    implicit none
    TYPE (complextaylor) cdivt
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(cdivt)

    cdivt%r=REAL(sc,kind=DP)/s1
    cdivt%i=aimag(sc)/s1
    master=localmaster
  END FUNCTION cdivt

  FUNCTION ctdiv( S1,sc )
    implicit none
    TYPE (complextaylor) ctdiv
    TYPE (taylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    complex(dp) w
    integer localmaster
    localmaster=master
    call ass(ctdiv)
    w=1.0_dp/sc
    ctdiv%r=s1*REAL(w,kind=DP)
    ctdiv%i=s1*aimag(w)
    master=localmaster
  END FUNCTION ctdiv

  FUNCTION dscmul( sc,S1 )
    implicit none
    TYPE (complextaylor) dscmul
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(dscmul)

    dscmul%r=sc*s1%r
    dscmul%i=sc*s1%i
    master=localmaster

  END FUNCTION dscmul

  FUNCTION scmul( sc,S1 )
    implicit none
    TYPE (complextaylor) scmul
    TYPE (complextaylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(scmul)

    scmul%r=sc*s1%r
    scmul%i=sc*s1%i
    master=localmaster

  END FUNCTION scmul

  FUNCTION iscmul( sc,S1 )
    implicit none
    TYPE (complextaylor) iscmul
    TYPE (complextaylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(iscmul)

    iscmul%r=sc*s1%r
    iscmul%i=sc*s1%i
    master=localmaster

  END FUNCTION iscmul


  FUNCTION cmulsc( S1,sc )
    implicit none
    TYPE (complextaylor) cmulsc
    TYPE (complextaylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(cmulsc)

    cmulsc%r=REAL(sc,kind=DP)*s1%r-aimag(sc)*s1%i
    cmulsc%i=REAL(sc,kind=DP)*s1%i+aimag(sc)*s1%r
    master=localmaster
  END FUNCTION cmulsc



  FUNCTION dmulsc( S1, sc)
    implicit none
    TYPE (complextaylor) dmulsc
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(dmulsc)

    dmulsc%r=sc*s1%r
    dmulsc%i=sc*s1%i
    master=localmaster

  END FUNCTION dmulsc

  FUNCTION mulsc( S1, sc)
    implicit none
    TYPE (complextaylor) mulsc
    TYPE (complextaylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(mulsc)

    mulsc%r=sc*s1%r
    mulsc%i=sc*s1%i
    master=localmaster

  END FUNCTION mulsc

  FUNCTION imulsc( S1, sc)
    implicit none
    TYPE (complextaylor) imulsc
    TYPE (complextaylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(imulsc)

    imulsc%r=sc*s1%r
    imulsc%i=sc*s1%i
    master=localmaster

  END FUNCTION imulsc

  SUBROUTINE  EQUAL(S2,S1)
    implicit none
    type (complextaylor),INTENT(inOUT)::S2
    type (complextaylor),INTENT(IN)::S1
    call check_snake
    !    master=0
    S2%R=S1%R
    S2%I=S1%I

  END SUBROUTINE EQUAL

  SUBROUTINE  ctEQUAL(S2,S1)
    implicit none
    type (complextaylor),INTENT(inOUT)::S2
    type (taylor),INTENT(IN)::S1
    call check_snake

    S2%R=S1
    S2%I=0.0_dp

  END SUBROUTINE ctEQUAL

  SUBROUTINE  tcEQUAL(S1,S2)
    implicit none
    type (complextaylor),INTENT(in)::S2
    type (taylor),INTENT(inout)::S1
    call check_snake
    !    master=0
    S1=S2%R

  END SUBROUTINE tcEQUAL


  SUBROUTINE  CEQUAL(R1,S2)          ! 2002.12.22
    implicit none
    type (complextaylor),INTENT(IN)::S2
    COMPLEX(dp), INTENT(inOUT)::R1
    call check_snake

    R1=S2.SUB.'0'
  END SUBROUTINE CEQUAL

  SUBROUTINE  DEQUAL(R1,S2)     ! 2002.12.22
    implicit none
    type (complextaylor),INTENT(IN)::S2
    real(dp), INTENT(inOUT)::R1
    call check_snake

    R1=S2.SUB.'0'
  END SUBROUTINE DEQUAL

  SUBROUTINE  REQUAL(R1,S2)   ! 2002.12.22
    implicit none
    type (complextaylor),INTENT(IN)::S2
    REAL(SP), INTENT(inOUT)::R1

    if(real_warning) call real_stop
    call check_snake

    R1=S2.SUB.'0'

  END SUBROUTINE REQUAL


  SUBROUTINE  CEQUALDACON(S2,R1)
    implicit none
    type (complextaylor),INTENT(inout)::S2
    complex(dp), INTENT(IN)::R1
    call check_snake

    !    master=0

    S2%R=REAL(R1,kind=DP)
    S2%I=aimag(R1)

  END SUBROUTINE CEQUALDACON

  SUBROUTINE  dEQUALDACON(S2,R1)
    implicit none
    type (complextaylor),INTENT(inout)::S2
    real(dp) , INTENT(IN)::R1
    call check_snake


    S2%R=R1
    S2%I=0.0_dp

  END SUBROUTINE dEQUALDACON

  SUBROUTINE  EQUALDACON(S2,R1)
    implicit none
    type (complextaylor),INTENT(inout)::S2
    real(sp) , INTENT(IN)::R1
    if(real_warning) call real_stop

    call check_snake

    S2%R=REAL(R1,kind=DP)
    S2%I=0.0_dp

  END SUBROUTINE EQUALDACON

  SUBROUTINE  iEQUALDACON(S2,R1)
    implicit none
    type (complextaylor),INTENT(inout)::S2
    integer , INTENT(IN)::R1
    call check_snake


    S2%R=REAL(R1,kind=DP)
    S2%I=0.0_dp

  END SUBROUTINE iEQUALDACON

  FUNCTION add( S1, S2 )
    implicit none
    TYPE (complextaylor) add
    TYPE (complextaylor), INTENT (IN) :: S1, S2
    integer localmaster
    localmaster=master
    call ass(add)
    add%r=s1%r+s2%r
    add%i=s1%i+s2%i
    master=localmaster
  END FUNCTION add

  FUNCTION tadd( S1, S2 )
    implicit none
    TYPE (complextaylor) tadd
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(tadd)
    tadd%r=s1+s2%r
    tadd%i=s2%i
    master=localmaster
  END FUNCTION tadd

  FUNCTION addt(  S2,S1 )
    implicit none
    TYPE (complextaylor) addt
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(addt)
    addt%r=s1+s2%r
    addt%i=s2%i
    master=localmaster
  END FUNCTION addt

  FUNCTION tsub( S1, S2 )
    implicit none
    TYPE (complextaylor) tsub
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(tsub)
    tsub%r=s1-s2%r
    tsub%i=-s2%i
    master=localmaster
  END FUNCTION tsub

  FUNCTION subt(  S2,S1 )
    implicit none
    TYPE (complextaylor) subt
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(subt)
    subt%r=s2%r-s1
    subt%i=s2%i
    master=localmaster
  END FUNCTION subt

  FUNCTION tmul( S1, S2 )
    implicit none
    TYPE (complextaylor) tmul
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(tmul)
    tmul%r=s1*s2%r
    tmul%i=s1*s2%i
    master=localmaster
  END FUNCTION tmul

  FUNCTION mult(  S2,S1 )
    implicit none
    TYPE (complextaylor) mult
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(mult)
    mult%r=s1*s2%r
    mult%i=s1*s2%i
    master=localmaster
  END FUNCTION mult

  FUNCTION tdiv( S1, S2 )
    implicit none
    TYPE (complextaylor) tdiv,temp
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(tdiv)
    call alloc(temp)
    temp=1.0_dp/s2
    tdiv%r=s1*temp%r
    tdiv%i=s1*temp%i
    master=localmaster
    call kill(temp)
  END FUNCTION tdiv

  FUNCTION divt(S2 , S1 )
    implicit none
    TYPE (complextaylor) divt
    type (taylor) temp
    TYPE (complextaylor), INTENT (IN) ::  S2
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(divt)
    call alloc(temp)
    temp=1.0_dp/s1
    divt%r=temp*s2%r
    divt%i=temp*s2%i
    master=localmaster
    call kill(temp)
  END FUNCTION divt

  FUNCTION csubSC( S1,sc )
    implicit none
    TYPE (complextaylor) csubSC
    TYPE (complextaylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master

    call ass(csubSC)

    csubSC%r=-REAL(sc,kind=DP)+s1%r
    csubSC%i=-aimag(sc)+s1%i
    master=localmaster
  END FUNCTION csubSC

  FUNCTION DsubSC(S1,sc)
    implicit none
    TYPE (complextaylor) DsubSC
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(DsubSC)
    DsubSC%r=s1%r-sc
    DsubSC%i=s1%i
    master=localmaster
  END FUNCTION DsubSC

  FUNCTION subSC(S1,sc)
    implicit none
    TYPE (complextaylor) subSC
    TYPE (complextaylor), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: sc
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(subSC)
    subSC%r=s1%r-sc
    subSC%i=s1%i
    master=localmaster
  END FUNCTION subSC

  FUNCTION isubSC(S1,sc)
    implicit none
    TYPE (complextaylor) isubSC
    TYPE (complextaylor), INTENT (IN) :: S1
    integer , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(isubSC)
    isubSC%r=s1%r-sc
    isubSC%i=s1%i
    master=localmaster
  END FUNCTION isubSC

  FUNCTION cSCsub( sc,S1  )
    implicit none
    TYPE (complextaylor) cSCsub
    TYPE (complextaylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(cSCsub)

    cSCsub%r=REAL(sc,kind=DP)-s1%r
    cSCsub%i=aimag(sc)-s1%i
    master=localmaster
  END FUNCTION cSCsub


  FUNCTION DSCsub( sc,S1 )
    implicit none
    TYPE (complextaylor) DSCsub
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(DSCsub)
    DSCsub%r=sc-s1%r
    DSCsub%i=-s1%i
    master=localmaster
  END FUNCTION DSCsub

  FUNCTION SCsub( sc,S1 )
    implicit none
    TYPE (complextaylor) SCsub
    TYPE (complextaylor), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: sc
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(SCsub)
    SCsub%r=sc-s1%r
    SCsub%i=-s1%i
    master=localmaster
  END FUNCTION SCsub

  FUNCTION iSCsub( sc,S1 )
    implicit none
    TYPE (complextaylor) iSCsub
    TYPE (complextaylor), INTENT (IN) :: S1
    integer , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(iSCsub)
    iSCsub%r=sc-s1%r
    iSCsub%i=-s1%i
    master=localmaster
  END FUNCTION iSCsub

  FUNCTION unarysub( S1 )
    implicit none
    TYPE (complextaylor) unarysub
    TYPE (complextaylor), INTENT (IN) :: S1

    integer localmaster
    localmaster=master
    call ass(unarysub)
    unarysub%r=-s1%r
    unarysub%i=-s1%i
    master=localmaster
  END FUNCTION unarysub

  FUNCTION subs( S1, S2 )
    implicit none
    TYPE (complextaylor) subs
    TYPE (complextaylor), INTENT (IN) :: S1, S2
    integer localmaster
    localmaster=master
    call ass(subs)
    subs%r=s1%r-s2%r
    subs%i=s1%i-s2%i
    master=localmaster
  END FUNCTION subs

  FUNCTION cSCADD( sc,S1  )
    implicit none
    TYPE (complextaylor) cSCADD
    TYPE (complextaylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(cSCADD)

    cSCADD%r=REAL(sc,kind=DP)+s1%r
    cSCADD%i=aimag(sc)+s1%i
    master=localmaster
  END FUNCTION cSCADD

  FUNCTION cADDSC(  S1,sc )
    implicit none
    TYPE (complextaylor) cADDSC
    TYPE (complextaylor), INTENT (IN) :: S1
    complex(dp), INTENT (IN) :: sc
    integer localmaster
    localmaster=master

    call ass(cADDSC)

    cADDSC%r=REAL(sc,kind=DP)+s1%r
    cADDSC%i=aimag(sc)+s1%i
    master=localmaster
  END FUNCTION cADDSC

  FUNCTION DADDSC( S1, sc )
    implicit none
    TYPE (complextaylor) DADDSC
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(DADDSC)
    DADDSC%r=sc+s1%r
    DADDSC%i=s1%i
    master=localmaster
  END FUNCTION DADDSC

  FUNCTION ADDSC( S1, sc )
    implicit none
    TYPE (complextaylor) ADDSC
    TYPE (complextaylor), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: sc
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(ADDSC)
    ADDSC%r=sc+s1%r
    ADDSC%i=s1%i
    master=localmaster
  END FUNCTION ADDSC

  FUNCTION iADDSC( S1, sc )
    implicit none
    TYPE (complextaylor) iADDSC
    TYPE (complextaylor), INTENT (IN) :: S1
    integer , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(iADDSC)
    iADDSC%r=sc+s1%r
    iADDSC%i=s1%i
    master=localmaster
  END FUNCTION iADDSC

  FUNCTION DSCADD( sc,S1 )
    implicit none
    TYPE (complextaylor) DSCADD
    TYPE (complextaylor), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(DSCADD)
    DSCADD%r=sc+s1%r
    DSCADD%i=s1%i
    master=localmaster
  END FUNCTION DSCADD

  FUNCTION SCADD( sc,S1 )
    implicit none
    TYPE (complextaylor) SCADD
    TYPE (complextaylor), INTENT (IN) :: S1
    real(sp) , INTENT (IN) :: sc
    integer localmaster
    if(real_warning) call real_stop
    localmaster=master
    call ass(SCADD)
    SCADD%r=sc+s1%r
    SCADD%i=s1%i
    master=localmaster
  END FUNCTION SCADD

  FUNCTION iSCADD( sc,S1 )
    implicit none
    TYPE (complextaylor) iSCADD
    TYPE (complextaylor), INTENT (IN) :: S1
    integer , INTENT (IN) :: sc
    integer localmaster
    localmaster=master
    call ass(iSCADD)
    iSCADD%r=sc+s1%r
    iSCADD%i=s1%i
    master=localmaster
  END FUNCTION iSCADD

  FUNCTION unaryADD( S1 )
    implicit none
    TYPE (complextaylor) unaryADD
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master
    call ass(unaryADD)
    unaryADD%r=s1%r
    unaryADD%i=s1%i
    master=localmaster
  END FUNCTION unaryADD

  subroutine inv( S1, s2 )
    implicit none
    TYPE (complextaylor), INTENT (IN) :: S1
    TYPE (complextaylor), INTENT (inout) :: S2
    TYPE (complextaylor)  s,ss
    complex(dp) d1
    real(dp) r1 ,i1
    integer i

    call alloc(s)
    call alloc(ss)


    r1=s1%r.sub.'0'
    i1=s1%i.sub.'0'
    s=s1
    d1=cmplx(r1,i1,kind=dp)
    d1=1.0_dp/d1

    s=d1*s1

    s=s-1.0_dp
    s=(-1.0_dp)*s

    ss=cmplx(1.0_dp,0.0_dp,kind=dp)
    s2=cmplx(1.0_dp,0.0_dp,kind=dp)

    do i=1,no
       ss=ss*s
       s2=s2+ss
    enddo

    s2=d1*s2

    call kill(s)
    call kill(ss)
  END subroutine inv

  FUNCTION logtpsat( S1 )
    implicit none
    TYPE (complextaylor) logtpsat
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master

    call ass(logtpsat)

    call logtpsa(s1,logtpsat)
    master=localmaster
  END FUNCTION logtpsat

  subroutine logtpsa( S1, s2 )
    implicit none
    TYPE (complextaylor) S1
    TYPE (complextaylor)  S2
    TYPE (complextaylor)  s,ss
    complex(dp) d1
    real(dp) r1 ,i1
    integer i


    call alloc(s)
    call alloc(ss)

    r1=s1%r.sub.'0'
    i1=s1%i.sub.'0'

    d1=cmplx(r1,i1,kind=dp)
    s=(1.0_dp/d1)*s1-1.0_dp

    s2=s


    ss=s

    do i=2,no
       ss=cmplx(-1.0_dp,0.0_dp,kind=dp)*ss*s
       s2=s2+ss/REAL(i,kind=DP)
    enddo

    s2=log(d1)+s2

    call kill(s)
    call kill(ss)
  END subroutine logtpsa

  FUNCTION full_abstpsat( S1 )
    implicit none
    real(dp) full_abstpsat ,r1,r2
    TYPE (complextaylor), INTENT (IN) :: S1


    r1=full_abs(s1%r) ! 2002.10.17
    r2=full_abs(s1%i)
    full_abstpsat=r1+r2
    !abstpsat=SQRT((s1%r.sub.'0')**2+(s1%i.sub.'0')**2)

  END FUNCTION full_abstpsat

  FUNCTION abstpsat( S1 )
    implicit none
    real(dp) abstpsat ,r1,r2
    TYPE (complextaylor), INTENT (IN) :: S1

    r1=abs(s1%r) ! 2002.10.17    etienne crap
    r2=abs(s1%i)
    abstpsat=SQRT(r1**2+r2**2)

  END FUNCTION abstpsat

  FUNCTION dcmplxt( S1,s2 )
    implicit none
    TYPE (complextaylor) dcmplxt
    TYPE (taylor), INTENT (IN) :: S1,s2
    integer localmaster
    localmaster=master

    call ass(dcmplxt)

    dcmplxt%r=s1
    dcmplxt%i=s2

    master=localmaster
  END FUNCTION dcmplxt

  FUNCTION CONJGT( S1 )
    implicit none
    TYPE (complextaylor) CONJGT
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master

    call ass(CONJGT)

    CONJGT%r=s1%R
    CONJGT%i=-S1%I

    master=localmaster
  END FUNCTION CONJGT




  FUNCTION datant( S1 )
    implicit none
    TYPE (complextaylor) datant,temp
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master

    call ass(datant)
    call alloc(temp)

    temp=(1.0_dp+s1*i_)

    temp=temp/(1.0_dp-s1*i_)

    temp=log(temp)
    datant=temp/2.0_dp/i_

    call kill(temp)

    master=localmaster
  END FUNCTION datant

  FUNCTION datantt( S1 )
    implicit none
    TYPE (taylor) datantt
    TYPE (complextaylor) temp
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master

    call ass(datantt)
    call alloc(temp)

    temp%r=s1
    temp=atan(temp)
    datantt=temp%r
    call kill(temp)

    master=localmaster
  END FUNCTION datantt

  FUNCTION dasintt( S1 )
    implicit none
    TYPE (taylor) dasintt
    TYPE (complextaylor) temp
    TYPE (taylor), INTENT (IN) :: S1
    integer localmaster
    real(dp) a0
    localmaster=master
    call ass(dasintt)
    call alloc(temp)

    temp%r=s1
    a0=abs(temp%r)
    if(a0>1.0_dp) then
       dasintt%i=0
       check_stable=.false.
       stable_da=.false.
       messagelost= "l_complex_taylor.f90 dasintt : abs(x)>1 Not defined in dasintt of complex_taylor "
    endif

    temp=asin(temp)
    dasintt=temp%r
    call kill(temp)

    master=localmaster
  END FUNCTION dasintt

  FUNCTION dasint( S1 )
    implicit none
    TYPE (complextaylor) dasint,temp
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master

    call ass(dasint)
    call alloc(temp)

    temp=(1.0_dp-s1**2)
    temp=temp**(0.5_dp)

    temp=i_*s1+ temp
    dasint=-i_*log(temp)

    call kill(temp)

    master=localmaster
  END FUNCTION dasint




  FUNCTION dacostt( S1 )
    implicit none
    TYPE (taylor) dacostt
    TYPE (complextaylor) temp
    TYPE (taylor), INTENT (IN) :: S1
    TYPE (taylor) t
    integer localmaster
    real(dp) a0
    localmaster=master

    call ass(dacostt)
    call alloc(temp)
    call alloc(t)

    a0=abs(s1)
    if(a0>1.0_dp) then
       dacostt%i=0
       check_stable=.false.
       stable_da=.false.
       messagelost= "l_complex_taylor.f90 dacostt : abs(x)>1 Not defined in dacostt of complex_taylor "
    endif

    !   if(debug_flag) then
    !    if(debug_acos) then
    temp%r=s1
    temp=acos(temp)
    !    else
    !     temp=-one+s1**2
    !     temp=temp**(half)

    !      temp=(s1+ temp)
    !      temp=-i_*log(temp)
    !     endif
    !    else
    !     t=sqrt(one-s1**2)

    !     temp=(s1+ i_*t)
    !     temp=-i_*log(temp)
    !    endif
    !   call print(temp%r,10)
    dacostt=temp%r

    call kill(t)
    call kill(temp)

    master=localmaster
  END FUNCTION dacostt



  FUNCTION dacost( S1 )
    implicit none
    TYPE (complextaylor) dacost,temp
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master

    call ass(dacost)
    call alloc(temp)

    temp=1.0_dp-s1**2
    temp=temp**(0.5_dp)

    temp=(s1+ i_*temp)
    dacost=-i_*log(temp)

    call kill(temp)

    master=localmaster
  END FUNCTION dacost




  FUNCTION tant( S1 )
    implicit none
    TYPE (complextaylor) tant, temp
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(tant)

    call alloc(temp)

    temp=exp(i_*s1)
    temp=temp-exp(-i_*s1)
    tant=exp(i_*s1)
    tant=tant+exp(-i_*s1)
    tant=temp/tant/i_

    call kill(temp)

    master=localmaster
  END FUNCTION tant

  FUNCTION dtanht( S1 )
    implicit none
    TYPE (complextaylor) dtanht, temp
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dtanht)

    call alloc(temp)

    temp=exp(s1)
    temp=temp-exp(-s1)
    dtanht=exp(s1)
    dtanht=dtanht+exp(-s1)
    dtanht=temp/dtanht

    call kill(temp)

    master=localmaster
  END FUNCTION dtanht

  FUNCTION dcost( S1 )
    implicit none
    TYPE (complextaylor) dcost
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dcost)

    dcost=exp(i_*s1)
    dcost=dcost+exp(-i_*s1)
    dcost=dcost/2.0_dp

    master=localmaster
  END FUNCTION dcost

  FUNCTION dcosht( S1 )
    implicit none
    TYPE (complextaylor) dcosht
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dcosht)

    dcosht=exp(s1)
    dcosht=dcosht+exp(-s1)
    dcosht=dcosht/2.0_dp

    master=localmaster
  END FUNCTION dcosht


  FUNCTION dsint( S1 )
    implicit none
    TYPE (complextaylor) dsint
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dsint)

    dsint=exp(i_*s1)
    dsint=dsint-exp(-i_*s1)
    dsint=dsint/2.0_dp/i_

    master=localmaster
  END FUNCTION dsint

  FUNCTION dsinht( S1 )
    implicit none
    TYPE (complextaylor) dsinht
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dsinht)

    dsinht=exp(s1)
    dsinht=dsinht-exp(-s1)
    dsinht=dsinht/2.0_dp

    master=localmaster
  END FUNCTION dsinht


  FUNCTION exptpsat( S1 )
    implicit none
    TYPE (complextaylor) exptpsat
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(exptpsat)
    call exptpsa(s1,exptpsat)
    master=localmaster
  END FUNCTION exptpsat

  FUNCTION dsqrtt( S1 )
    implicit none
    TYPE (complextaylor) dsqrtt
    TYPE (complextaylor), INTENT (IN) :: S1
    integer localmaster
    localmaster=master


    call ass(dsqrtt)
    dsqrtt= S1**0.5_dp

    master=localmaster
  END FUNCTION dsqrtt

  subroutine exptpsa( S1, s2 )
    implicit none
    TYPE (complextaylor), INTENT (IN) :: S1
    TYPE (complextaylor), INTENT (inout) :: S2
    TYPE (complextaylor)  s,ss
    complex(dp) d1
    real(dp) r1 ,i1
    integer i
    call alloc(s)
    call alloc(ss)


    r1=s1%r.sub.'0'
    i1=s1%i.sub.'0'

    d1=cmplx(r1,i1,kind=dp)
    s=s1

    s=s1-d1


    ss=cmplx(1.0_dp,0.0_dp,kind=dp)
    s2=cmplx(1.0_dp,0.0_dp,kind=dp)

    do i=1,no
       ss=ss*s
       ss=ss/REAL(i,kind=DP)
       s2=s2+ss
    enddo

    s2=exp(d1)*s2

    call kill(s)
    call kill(ss)

  END subroutine exptpsa

  subroutine assc(s1)
    implicit none
    TYPE (complextaylor) s1
    !  lastmaster=master  ! 2002.12.13

    select case(master)
    case(0:ndumt-1)
       master=master+1
    case(ndumt)
 
         write(6,*) " cannot indent anymore assc" 
 
    end select
    !    write(26,*) " complex  taylor ",master

    call ass0(s1%r)
    call ass0(s1%i)


  end subroutine ASSc



  subroutine KILL_TPSA()
    IMPLICIT NONE
    logical present_tpsa

    present_tpsa=lingyun_yang
    !    if(.not.first_time) then
    if(last_tpsa==1) then
       lingyun_yang=.true.
       call KILL(varc1)
       call KILL(varc2)
       CALL KILL_fpp   ! IN TPSALIE_ANALISYS
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       call KILL(varc1)
       call KILL(varc2)
       CALL KILL_fpp   ! IN TPSALIE_ANALISYS
    endif
    lingyun_yang=default_tpsa
    last_tpsa=0
    !    endif
    !    first_time=.true.

  END subroutine KILL_TPSA

  subroutine init_map_c(NO1,ND1,NP1,NDPT1,log)
    implicit none
    integer NO1,ND1,NP1,NDPT1
    LOGICAL(lp) log,present_tpsa
    present_tpsa=lingyun_yang
    !    if(.not.first_time) then
    if(last_tpsa==1) then
       lingyun_yang=.true.
       call kill(varc1)
       call kill(varc2)
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       call kill(varc1)
       call kill(varc2)
    endif
    lingyun_yang=present_tpsa
    !    endif

    call init_map(NO1,ND1,NP1,NDPT1,log)
    call set_in_complex(log)
    call alloc(varc1)
    call alloc(varc2)

  end subroutine  init_map_c

  subroutine init_tpsa_c(NO1,NP1,log)
    implicit none
    integer NO1,NP1
    LOGICAL(lp) log,present_tpsa
    present_tpsa=lingyun_yang
    !    if(.not.first_time) then
    if(last_tpsa==1) then
       lingyun_yang=.true.
       call kill(varc1)
       call kill(varc2)
    elseif(last_tpsa==2) then
       lingyun_yang=.false.
       call kill(varc1)
       call kill(varc2)
    endif
    lingyun_yang=present_tpsa
    !    endif
    call init_tpsa(NO1,NP1,log)
    call set_in_complex(log)
    call alloc(varc1)
    call alloc(varc2)

  end subroutine  init_tpsa_c


  subroutine set_in_complex(log)
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
  end  subroutine set_in_complex

  !  SUBROUTINE  VARcC(S1,R1,R2,I1,I2)
  !    implicit none
  !    INTEGER,INTENT(IN)::I1,I2
  !    complex(dp),INTENT(IN)::R1
  !    complex(dp),INTENT(IN)::R2                ! big change
  !    type (complextaylor),INTENT(INOUT)::S1
  !    integer localmaster
  !    localmaster=master
  !
  !      s1=r1+r2*((one.mono.i1) + i_*    (one.mono.i2) )
  !!     s1%r=(/REAL(R1,kind=DP),R2/).var.i1
  !!     s1%i=(/aimag(R1),R2/).var.i2
  !!    call var001(s1%r,REAL(R1,kind=DP),R2,i1)
  !!    call var001(s1%i,aimag(R1),R2,i2)
  !    master=localmaster
  !
  !
  !
  !  END SUBROUTINE VARcC
  !
  !  SUBROUTINE  VARc(S1,R1,I1,I2)
  !    implicit none
  !    INTEGER,INTENT(IN)::I1,I2
  !    complex(dp),INTENT(IN)::R1
  !    type (complextaylor),INTENT(INOUT)::S1
  !
  !    integer localmaster
  !    localmaster=master
  !
  !     s1%r=REAL(R1,kind=DP).var.i1
  !     s1%i=aimag(R1).var.i2
  !!    call var000(s1%r,REAL(R1,kind=DP),i1)
  !!    call var000(s1%i,aimag(R1),i2)
  !    master=localmaster
  !
  !
  !  END SUBROUTINE VARc


  !  SUBROUTINE  shiftc(S1,S2,s)
  !    implicit none
  !    INTEGER,INTENT(IN)::s
  !    type (complextaylor),INTENT(IN)::S1
  !    type (complextaylor),INTENT(inout)::S2
  !
  !    call shift000(S1%r,S2%r,s)
  !    call shift000(S1%i,S2%i,s)

  !  END SUBROUTINE shiftc

  FUNCTION GETintk( S1, S2 )
    implicit none
    TYPE (complextaylor) GETintk
    TYPE (complextaylor), INTENT (IN) :: S1
    integer  , INTENT (IN) ::  S2

    integer localmaster
    localmaster=master

    call ass(GETintk)

    GETintk%r=S1%r<=s2
    GETintk%i=S1%i<=s2

    !    call  shiftda(GETintk,GETintk, s2 )

    master=localmaster


  END FUNCTION GETintk



  SUBROUTINE  pekc(S1,J,R1)
    implicit none
    INTEGER,INTENT(IN),dimension(:)::j
    complex(dp),INTENT(inout)::R1
    type (complextaylor),INTENT(IN)::S1
    real(dp) xr,xi

    call pek000(s1%r,j,xr)
    call pek000(s1%i,j,xi)

    r1=cmplx(xr,xi,kind=dp)

  END SUBROUTINE pekc


  SUBROUTINE  pokc(S1,J,R1)
    implicit none
    INTEGER,INTENT(in),dimension(:)::j
    complex(dp),INTENT(in)::R1
    type (complextaylor),INTENT(inout)::S1

    call pok000(s1%r,J,REAL(r1,kind=DP))
    call pok000(s1%i,J,aimag(r1))
  END SUBROUTINE pokc

  SUBROUTINE  CFUC(S2,FUN,S1)!
    implicit none
    type (complextaylor),INTENT(INOUT)::S1
    type (complextaylor),INTENT(IN)::S2
    type (complextaylor) T
    type (taylor) W
    complex(dp) FUN
    EXTERNAL FUN
    CALL ALLOC(T)
    CALL ALLOC(W)

    CALL CFUR(S2%R,FUN,W)
    T%R=W
    CALL CFUI(S2%I,FUN,W)
    T%R=T%R-W
    CALL CFUR(S2%I,FUN,W)
    T%I=W
    CALL CFUI(S2%R,FUN,W)
    T%I=T%I+W
    S1=T
    CALL KILL(T)
    CALL KILL(W)

  END SUBROUTINE CFUC

  SUBROUTINE  CFURES(S2,FUN,S1)!
    implicit none
    type (pbresonance),INTENT(INOUT)::S1
    type (pbresonance),INTENT(IN)::S2
    type (complextaylor) T
    type (taylor) W
    complex(dp) FUN
    EXTERNAL FUN
    CALL ALLOC(T)
    CALL ALLOC(W)

    CALL CFUR(S2%COS%H,FUN,W)
    T%R=W
    CALL CFUI(S2%SIN%H,FUN,W)
    T%R=T%R-W
    CALL CFUR(S2%SIN%H,FUN,W)
    T%I=W
    CALL CFUI(S2%COS%H,FUN,W)
    T%I=T%I+W
    S1%COS%H=T%R
    S1%SIN%H=T%I
    CALL KILL(T)
    CALL KILL(W)

  END SUBROUTINE CFURES


end module  complex_taylor
