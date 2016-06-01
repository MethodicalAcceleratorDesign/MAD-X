!The Full Polymorphic Package
!Copyright (C) Etienne Forest


MODULE TPSA
  !use newda
  use definition
  use file_handler
  IMPLICIT NONE
  public
  integer,private::ndel ,nd2par,nd2part,nd2partt
  integer,private,dimension(lnv)::jfil,jfilt

  private equal,DAABSEQUAL,Dequaldacon ,equaldacon ,Iequaldacon  !,AABSEQUAL 2002.10.17
  private pow,powr,powr8,dlogt, GETORDER,CUTORDER,getchar,GETint
  private getdiff,getdATRA  ,mul,dmulsc,dscmul
  private mulsc,scmul,imulsc,iscmul
  private div,ddivsc,dscdiv,divsc,scdiv,idivsc,iscdiv
  private unaryADD,add,daddsc,dscadd,addsc,scadd,iaddsc,iscadd 
  private unarySUB,subs,dsubsc,dscsub,subsc,scsub,isubsc,iscsub
  private allocda,KILLda,A_OPT,K_opt
  private dexpt,dcost,dsint,dsqrtt,dtant,datanht,dtanht
  PRIVATE GETCHARnd2,GETintnd2,dputchar,dputint, filter,check_j,dsinHt,dCOSHt
  private GETintnd2t,print_for_bmad_parse
  PRIVATE DEQUAL,REQUAL,varf,varf001  !,CHARINT
  !  PUBLIC VAR,ASS
  private pbbra,full_absT,asstaylor,getcharnd2s,GETintnd2s,GETintk
  private shiftda,shift000
  !PRIVATE null_0,ALLOC_U,FILL_N,REFILL_N
  !  public, alloc_uni, null_uni, fill_uni, refill_uni

  private fill_uni_r ! new sagan

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
private bessi_se,bessi0_se,poly_e,bessi1_se,I_nt,I_nr,In_nt,In_enz
real(dp) :: switch_bessel=0.001d0

  type(dalevel) scratchda(ndumt)   !scratch levels of DA using linked list

  INTERFACE I_ns
     MODULE PROCEDURE I_nr
     MODULE PROCEDURE I_nt
  END INTERFACE

 INTERFACE I_n
     MODULE PROCEDURE In_enz
     MODULE PROCEDURE In_nt
  END INTERFACE

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

     MODULE PROCEDURE fill_uni_r
     MODULE PROCEDURE null_uni
     MODULE PROCEDURE fill_uni  ! new sagan
     MODULE PROCEDURE refill_uni
  end  INTERFACE



  INTERFACE print_for_bmad_parser
     MODULE PROCEDURE print_for_bmad_parse
  END INTERFACE


  INTERFACE print
     MODULE PROCEDURE printunitaylor
  END INTERFACE




  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber2" width="400" height="135">
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">+</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@         <font face="Times New Roman" size="1">Taylor</font></span></td>
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
  !@         <font face="Times New Roman" size="1">Taylor</font></span></td>
  !@         <td width="77" height="20" align="center">
  !@         <span style="text-transform: uppercase; font-weight:700">
  !@         <font face="Times New Roman" size="1">
  !@         <a href="i_tpsa.htm#ADD" style="text-decoration: none">add</a></font></span></td>
  !@         <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase; font-weight:700">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DADDSC" style="text-decoration: none">daddsc</a></font></span></td>
  !@        <td width="78" height="20" align="center"><b>
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#ADDSC" style="text-decoration: none">ADDSC</a></font></b></td>
  !@        <td width="56" height="20" align="center"><b>
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#IADDSC" style="text-decoration: none">
  !@        IADDSC</a></font></b></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase; font-weight:700">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DSCADD" style="text-decoration: none">dscadd</a></font></span></td>
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
  !@        <td width="77" height="20" align="center"><b>
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#SCADD" style="text-decoration: none">SCADD</a></font></b></td>
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
  !@        <td width="77" height="20" align="center"><b>
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#ISCADD" style="text-decoration: none">
  !@        ISCADD</a></font></b></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@     </table>
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unaryADD  !@2 This is a unary operation
     MODULE PROCEDURE add
     MODULE PROCEDURE daddsc
     MODULE PROCEDURE dscadd
     MODULE PROCEDURE addsc
     MODULE PROCEDURE scadd
     MODULE PROCEDURE iaddsc
     MODULE PROCEDURE iscadd
  END INTERFACE




  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="400" height="135">
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">-</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Taylor</font></span></td>
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
  !@        <font face="Times New Roman" size="1">Taylor</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#SUBS" style="text-decoration: none; font-weight: 700">SUBS</a></font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DSUBSC" style="text-decoration: none; font-weight: 700">dSUBsc</a></font></span></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#SUBSC" style="text-decoration: none; font-weight: 700">SUBSC</a></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#ISUBSC" style="text-decoration: none; font-weight: 700">ISUBSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DSCSUB" style="text-decoration: none; font-weight: 700">dscSUB</a></font></span></td>
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
  !@        <a href="i_tpsa.htm#SCSUB" style="text-decoration: none; font-weight: 700">
  !@        SCSUB</a></font></td>
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
  !@        <a href="i_tpsa.htm#ISCSUB" style="text-decoration: none; font-weight: 700">
  !@        ISCSUB</a></font></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@     </table>
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



  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="400" height="134">
  !@      <tr>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">*</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Taylor</font></span></td>
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
  !@        <font face="Times New Roman" size="1">Taylor</font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#MUL" style="text-decoration: none; font-weight:700">MUL</a></font></span></td>
  !@        <td width="77" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DMULSC" style="text-decoration: none; font-weight:700">dMULsc</a></font></span></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#MULSC" style="text-decoration: none; font-weight:700">MULSC</a></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="1" face="Times New Roman">
  !@        <a href="i_tpsa.htm#IMULSC" style="text-decoration: none; font-weight:700">IMULSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="77" height="19" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="77" height="19" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DSCMUL" style="text-decoration: none; font-weight:700">dscMUL</a></font></span></td>
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
  !@        <a href="i_tpsa.htm#SCMUL" style="text-decoration: none; font-weight:700">
  !@        SCMUL</a></font></td>
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
  !@        <a href="i_tpsa.htm#ISCMUL" style="text-decoration: none; font-weight:700">
  !@        ISCMUL</a></font></td>
  !@        <td width="77" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="78" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="56" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@     </table>

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul
     MODULE PROCEDURE dmulsc
     MODULE PROCEDURE dscmul
     MODULE PROCEDURE mulsc
     MODULE PROCEDURE scmul
     MODULE PROCEDURE imulsc
     MODULE PROCEDURE iscmul
  END INTERFACE

  !@    <table border="4" cellspacing="1" bordercolor="#000000" id="AutoNumber1" width="400" height="135">
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">/</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">Taylor</font></span></td>
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
  !@        <font face="Times New Roman" size="1">Taylor</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DIV" style="text-decoration: none; font-weight: 700">div</a></font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DDIVSC" style="text-decoration: none; font-weight: 700">dDIVsc</a></font></span></td>
  !@        <td width="0" height="20" align="center"><font size="1">
  !@        <a href="i_tpsa.htm#DIVSC" style="text-decoration: none; font-weight: 700">DIVSC</a></font></td>
  !@        <td width="0" height="20" align="center"><font size="1">
  !@        <a href="i_tpsa.htm#IDIVSC" style="text-decoration: none; font-weight: 700">
  !@        IDIVSC</a></font></td>
  !@      </tr>
  !@      <tr>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        Real(dp)</font></span></td>
  !@        <td width="0" height="20" align="center">
  !@        <span style="text-transform: uppercase">
  !@        <font face="Times New Roman" size="1">
  !@        <a href="i_tpsa.htm#DSCDIV" style="text-decoration: none; font-weight: 700">dscDIV</a></font></span></td>
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
  !@        <a href="i_tpsa.htm#SCDIV" style="text-decoration: none; font-weight: 700">
  !@        SCDIV</a></font></td>
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
  !@        <a href="i_tpsa.htm#ISCDIV" style="text-decoration: none; font-weight: 700">
  !@        ISCDIV</a></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@        <td width="0" height="20" align="center">
  !@        <font size="2" face="Times New Roman"><b>F90</b></font></td>
  !@      </tr>
  !@    </table>

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div
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
  END INTERFACE
  INTERFACE dabs
     MODULE PROCEDURE DAABSEQUAL  ! remove 2002.10.17
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
     MODULE PROCEDURE cfu000  ! not private
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

  subroutine fliptaylor(xy,xyf,i)
    implicit none
    type(taylor), intent(inout) ::  xy,xyf
    integer i
    call flip_i(xy%i,xyf%i,i)
  end subroutine fliptaylor


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

    call dacmu(s1%i,-1.0_dp,temp%i)
    call dacop(temp%i,unarySUB%i)

    master=localmaster

  END FUNCTION unarySUB

  SUBROUTINE  maketree(S1,s2)
    implicit none
    type (TAYLOR),INTENT(IN)::S1
    type (TAYLOR),INTENT(inOUT):: s2
    IF(.NOT.C_%STABLE_DA) RETURN

    !    if(old) then
    call mtree((/s1%i/),1,(/s2%i/),1)
    !   else
    !      call newdacop(s1%j,s2%j)
    !   endif
  END SUBROUTINE maketree

  SUBROUTINE  allocda(S1)
    implicit none
    type (TAYLOR),INTENT(INOUT)::S1

    !    IF(first_time) THEN
    IF(last_tpsa==0) THEN
       w_p=0
       w_p%nc=1
       w_p=(/" No TPSA package ever initialized "/)
       w_p%fc='(1((1X,A72),/))'
       ! call !write_e(111)
    ENDIF
    !    if(old) then
    s1%i=0
    call etall1(s1%i)
    !    else
    !       call nullnewda(s1%j)
    !       call allocnewda(s1%j)
    !    endif
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
    !    if(old) then
    call DADAL1(s1%i)
    !    else
    !       call KILLNEWDAs(s1%j)
    !    endif

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
    if(s2%i==0) then
       call crap1("EQUAL 1 in tpsa") !call allocw(s2)
    endif
    if(s1%i==0) call crap1("EQUAL 2") ! call allocw(s1)
    CALL DACOP(S1%I,S2%I)
    !    else
    !      IF (.NOT. ASSOCIATED(s2%j%r)) call crap1("EQUAL 3") !call allocw(s2)
    !      IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("EQUAL 4") !call allocw(s1)
    !      call newdacop(S1%j,S2%j)
    !   endif
  END SUBROUTINE EQUAL

  SUBROUTINE  DEQUAL(R1,S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    real(dp), INTENT(inOUT)::R1
    IF(.NOT.C_%STABLE_DA) RETURN
    call check_snake

    R1=S2.SUB.'0'
  END SUBROUTINE DEQUAL

  SUBROUTINE  REQUAL(R1,S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    REAL(SP), INTENT(inOUT)::R1
    IF(.NOT.C_%STABLE_DA) RETURN

    if(real_warning) call real_stop
    call check_snake

    R1=S2.SUB.'0'

  END SUBROUTINE REQUAL

  function  DAABSEQUAL(S2)
    implicit none
    type (TAYLOR),INTENT(IN)::S2
    real(dp) DAABSEQUAL
    DAABSEQUAL=0
    IF(.NOT.C_%STABLE_DA) RETURN


    DAABSEQUAL=abs(S2.sub.'0')

  END function DAABSEQUAL


  SUBROUTINE  DEQUALDACON(S2,R1)
    implicit none
    type (TAYLOR),INTENT(inOUT)::S2
    real(dp), INTENT(IN)::R1
    IF(.NOT.C_%STABLE_DA) RETURN

    !    if(old) then
    if(s2%i==0)  call crap1("DEQUALDACON 1") !call allocw(s2)
    CALL DACON(S2%I,R1)
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
    if(s2%i==0) call crap1("EQUALDACON 1") !call allocw(s2)
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
    if(s2%i==0) call crap1("IEQUALDACON 1") !call allocw(s2)
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

    ! if(old) then
    call dafun('EXP ',s1%i,temp%i)
    call dacop(temp%i,dexpt%i)
    !    else
    !       call newdafun('EXP ',s1%j,dexpt%j)
    !    endif

    master=localmaster

  END FUNCTION dexpt

  FUNCTION FULL_ABST( S1 )
    implicit none
    real(dp) FULL_ABST
    TYPE (taylor), INTENT (IN) :: S1
    FULL_ABST=0
    IF(.NOT.C_%STABLE_DA) RETURN
    !    call check(s1)

    ! if(old) then
    CALL DAABS(S1%I,FULL_ABST)
    !    else
    !       CALL newDAABS(S1%j,FULL_ABST)
    !    endif

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

    ! if(old) then
    call dafun('SIN ',s1%i,temp%i)
    call dacop(temp%i,dtant%i)
    call dafun('COS ',s1%i,temp%i)
    call dadiv(dtant%i,temp%i,dtant%i)
    !    else
    !       call newdafun('SIN ',s1%j,temp%il)
    !       call newdacop(temp%il,dtant%j)
    !       call newdafun('COS ',s1%j,temp%il)
    !       call newdadiv(dtant%j,temp%il,dtant%j)
    !    endif

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

    datanht=log((1+s1)/sqrt(1-s1))/2.0_dp

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

    ! if(old) then
    call dafun('COS ',s1%i,temp%i)
    call dacop(temp%i,dcost%i)
    !    else
    !       call newdafun('COS ',s1%j,dcost%j)
    !    endif

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
    ! if(old) then
    call dafun('SIN ',s1%i,temp%i)
    call dacop(temp%i,dsint%i)
    !    else
    !       call newdafun('SIN ',s1%j,dsint%j)
    !    endif

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
    ! if(old) then
    call dafun('SINH',s1%i,temp%i)
    call dacop(temp%i,dsinHt%i)
    !    else
    !       call newdafun('SINH',s1%j,dsinHt%j)
    !    endif
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
    ! if(old) then
    call dafun('COSH',s1%i,temp%i)
    call dacop(temp%i,DCOSHT%i)
    !    else
    !       call newdafun('COSH',s1%j,DCOSHT%j)
    !    endif

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
    ! if(old) then

    dtanht=sinh(s1)/cosh(s1)
    !    else
    !       call newdafun('COSH',s1%j,DCOSHT%j)
    !    endif

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
    ! if(old) then
    call dafun('LOG ',s1%i,temp%i)
    call dacop(temp%i,dlogt%i)
    !    else
    !       call newdafun('LOG ',s1%j,dlogt%j)
    !    endif

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

    ! if(old) then
    call dafun('SQRT',s1%i,temp%i)
    call dacop(temp%i,dsqrtt%i)
    !    else
    !       call newdafun('SQRT',s1%j,dsqrtt%j)
    !    endif

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

    ! if(old) then
    call damul(s1%i,s2%i,temp%i)
    call dacop(temp%i,mul%i)
    !    else
    !       call newdamul(s1%j,s2%j,mul%j)
    !    endif

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

    ! if(old) then
    pbbra=0.0_dp
    do i=1,nd
       pbbra=(s1.d.(2*i-1))*(s2.d.(2*i))-(s2.d.(2*i-1))*(s1.d.(2*i))+pbbra
    enddo
    !    call DAPOI(s1%i,s2%i,temp%i,nd)
    !    call dacop(temp%i,pbbra%i)
    !    else
    !       call newDAPOI(s1%j,s2%j,temp%il,nd)
    !       call newdacop(temp%il,pbbra%j)
    !    endif

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

    ! if(old) then
    CALL TAKE(S1%I,S2,temp%i)
    call dacop(temp%i,GETORDER%i)
    !    else
    !       CALL NEWTAKE(S1%J,S2,temp%iL)
    !       call NEWdacop(temp%iL,GETORDER%J)
    !    endif
    master=localmaster

  END FUNCTION GETORDER




  FUNCTION CUTORDER( S1, S2 )
    implicit none
    TYPE (taylor) CUTORDER
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       CUTORDER%i=0
      RETURN
    endif
    localmaster=master

    !    call check(s1)
    call ass(CUTORDER)

    ! if(old) then
    call datrunc(S1%I,s2,CUTORDER%i)
    !    call dacop(S1%I,CUTORDER%i)

    !    DO I=S2,NO
    !       CALL TAKE(CUTORDER%I,I,TEMP)
    !       CALL DASUB(CUTORDER%I,TEMP,CUTORDER%I)
    !    ENDDO
    !    else
    !       call NEWdacop(S1%J,CUTORDER%J)
    !       DO I=S2,NO
    !          CALL NEWTAKE(CUTORDER%J,I,TEMPL)
    !          CALL NEWDASUB(CUTORDER%J,TEMPL,CUTORDER%J)
    !       ENDDO
    !    endif
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
       if(i>nv) then
          if(j(i)>0) then
             dputchar=0.0_dp
             !             call var(dputchar,zero,0)
             return
          endif
       endif
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
    localmaster=master


    call ass(dputint0)

    j=0
    if(s2>nv) then
       dputint0=S1
       return
    endif


    dputint0=0.0_dp
    !    call var(dputint0,zero,s2)

    j(s2)=1
    CALL pok(dputint0,j,s1)

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
    CALL dapek(S1%I,j,r1)
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
    CALL dapek(S1%I,j,r1)
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

    ! if(old) then
    CALL dader(S2,S1%I,temp%i)
    call dacop(temp%i,GETdiff%i)
    !    else
    !       CALL NEWdader(S2,S1%J,temp%iL)
    !       call NEWdacop(temp%iL,GETdiff%J)
    !    endif
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
    master=localmaster

  END FUNCTION GETINTegrate


  FUNCTION GETdatra( S1, S2 )
    implicit none
    TYPE (taylor) GETdatra
    TYPE (taylor), INTENT (IN) :: S1
    INTEGER, INTENT (IN) ::  S2
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       GETdatra%i=0
      RETURN
    endif
    localmaster=master


    !    call check(s1)
    call ass(GETdatra)

    ! if(old) then
    CALL datra(S2,S1%I,temp%i)
    call dacop(temp%i,GETdatra%i)
    !    else
    !       CALL NEWdatra(S2,S1%J,temp%iL)
    !       call NEWdacop(temp%iL,GETdatra%J)
    !    endif
    master=localmaster

  END FUNCTION GETdatra

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

    ! if(old) then
    CALL DACON(temp%i,1.0_dp)

    R22=IABS(R2)
    DO I=1,R22
       CALL DAMUL(temp%i,S1%I,temp%i)
    ENDDO
    IF(R2.LT.0) THEN
       CALL DADIC(temp%i,1.0_dp,temp%i)
    ENDIF
    call dacop(temp%i,POW%i)
 
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

    ! if(old) then
    call dacmu(s1%i,sc,temp%i)
    call dacop(temp%i,dmulsc%i)
    !    else
    !       call newdacmu(s1%j,sc,dmulsc%j)
    !    endif

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

    ! if(old) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,mulsc%i)
    !    else
    !       call newdacmu(s1%j,REAL(sc,kind=DP),mulsc%j)
    !    endif
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


    ! if(old) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,imulsc%i)
    !    else
    !       call newdacmu(s1%j,REAL(sc,kind=DP),imulsc%j)
    !    endif

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

    ! if(old) then
    call dacmu(s1%i,sc,temp%i)
    call dacop(temp%i,dscmul%i)
    !    else
    !       call newdacmu(s1%j,sc,dscmul%j)
    !    endif

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


    ! if(old) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scmul%i)
    !    else
    !       call newdacmu(s1%j,REAL(sc,kind=DP),scmul%j)
    !    endif

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

    ! if(old) then
    call dacmu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscmul%i)
    !    else
    !       call newdacmu(s1%j,REAL(sc,kind=DP),iscmul%j)
    !    endif

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

    ! if(old) then
    call dadiv(s1%i,s2%i,temp%i)
    call dacop(temp%i,div%i)
    !    else
    !       call newdadiv(s1%j,s2%j,templ)
    !       call newdacop(templ,div%j)
    !    endif

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

    ! if(old) then
    call dadic(s1%i,sc,temp%i)
    call dacop(temp%i,dscdiv%i)
    !    else
    !       call newdadic(s1%j,sc,dscdiv%j)
    !    endif

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


    ! if(old) then
    call dadic(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scdiv%i)
    !    else
    !       call newdadic(s1%j,REAL(sc,kind=DP),scdiv%j)
    !    endif

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

    ! if(old) then
    call dadic(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscdiv%i)
    !    else
    !       call newdadic(s1%j,REAL(sc,kind=DP),iscdiv%j)
    !    endif


    master=localmaster
  END FUNCTION iscdiv

  FUNCTION ddivsc( S1, sc )
    implicit none
    TYPE (taylor) ddivsc
    TYPE (taylor), INTENT (IN) :: S1
    real(dp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       ddivsc%i=0
      RETURN
    endif  
    localmaster=master


    !    call check(s1)
    call ass(ddivsc)


    ! if(old) then
    call dacdi(s1%i,sc,temp%i)
    call dacop(temp%i,ddivsc%i)
    !    else
    !       call newdacdi(s1%j,sc,ddivsc%j)
    !    endif
    master=localmaster

  END FUNCTION ddivsc

  FUNCTION divsc( S1, sc )
    implicit none
    TYPE (taylor) divsc
    TYPE (taylor), INTENT (IN) :: S1
    real(sp), INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       divsc%i=0
      RETURN
    endif 
    localmaster=master


    if(real_warning) call real_stop
    !    call check(s1)
    call ass(divsc)

    ! if(old) then
    call dacdi(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,divsc%i)
    !    else
    !       call newdacdi(s1%j,REAL(sc,kind=DP),divsc%j)
    !    endif
    master=localmaster

  END FUNCTION divsc


  FUNCTION idivsc( S1, sc )
    implicit none
    TYPE (taylor) idivsc
    TYPE (taylor), INTENT (IN) :: S1
    integer, INTENT (IN) :: sc
    integer localmaster
    IF(.NOT.C_%STABLE_DA) then
       idivsc%i=0
      RETURN
    endif 
    localmaster=master


    !    call check(s1)
    call ass(idivsc)


    ! if(old) then
    call dacdi(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,idivsc%i)
    !    else
    !       call newdacdi(s1%j,REAL(sc,kind=DP),idivsc%j)
    !    endif
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


    ! if(old) then
    call daadd(s1%i,s2%i,add%i)
    !  call dacop(temp,add%i)
    !  call daadd(s1%i,s2%i,temp)
    !  call dacop(temp,add%i)
    !    else
    !       call newdaadd(s1%j,s2%j,add%j)
    !    endif

    master=localmaster

  END FUNCTION add

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

    ! if(old) then
    call dacad(s1%i,sc,temp%i)
    call dacop(temp%i,daddsc%i)
    !    else
    !       call newdacad(s1%j,sc,daddsc%j)
    !    endif
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


    ! if(old) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,addsc%i)
    !    else
    !       call newdacad(s1%j,REAL(sc,kind=DP),addsc%j)
    !    endif
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

    ! if(old) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iaddsc%i)
    !    else
    !       call newdacad(s1%j,REAL(sc,kind=DP),iaddsc%j)
    !    endif
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

    ! if(old) then
    call dacad(s1%i,sc,temp%i)
    call dacop(temp%i,dscadd%i)
    !    else
    !       call newdacad(s1%j,sc,dscadd%j)
    !    endif
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

    ! if(old) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scadd%i)
    !    else
    !       call newdacad(s1%j,REAL(sc,kind=DP),scadd%j)
    !    endif

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


    ! if(old) then
    call dacad(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscadd%i)
    !    else
    !       call newdacad(s1%j,REAL(sc,kind=DP),iscadd%j)
    !    endif
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



    ! if(old) then
    call dasub(s1%i,s2%i,temp%i)
    call dacop(temp%i,subs%i)
    !    else
    !       call newdasub(s1%j,s2%j,subs%j)
    !    endif
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

    ! if(old) then
    call dacsu(s1%i,sc,temp%i)
    call dacop(temp%i,dsubsc%i)
    !    else
    !       call newdacsu(s1%j,sc,dsubsc%j)
    !    endif

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

    ! if(old) then
    call dacsu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,subsc%i)
    !    else
    !       call newdacsu(s1%j,REAL(sc,kind=DP),subsc%j)
    !    endif
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

    ! if(old) then
    call dacsu(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,isubsc%i)
    !    else
    !       call newdacsu(s1%j,REAL(sc,kind=DP),isubsc%j)
    !    endif
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

    ! if(old) then
    call dasuc(s1%i,sc,temp%i)
    call dacop(temp%i,dscsub%i)
    !    else
    !       call newdasuc(s1%j,sc,dscsub%j)
    !    endif
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

    ! if(old) then
    call dasuc(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,scsub%i)
    !    else
    !       call newdasuc(s1%j,REAL(sc,kind=DP),scsub%j)
    !    endif
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

    ! if(old) then
    call dasuc(s1%i,REAL(sc,kind=DP),temp%i)
    call dacop(temp%i,iscsub%i)
    !    else
    !       call newdasuc(s1%j,REAL(sc,kind=DP),iscsub%j)
    !    endif
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
    IF(.NOT.C_%STABLE_DA) RETURN

    ! if(old) then
    if(s2%i==0) call crap1("shift000  1" )  !call etall1(s2%i)
    CALL DAshift(s1%i,s2%i,s)
    !   else
    !      if(.NOT. ASSOCIATED(s2%j%r))call crap1("shift000  2" )   ! call newetall(s2%j,1)
    !
    !      CALL NEWDAshift(s1%j,s2%j,s)
    !   endif
    !
  END SUBROUTINE shift000


  SUBROUTINE  pek000(S1,J,R1)
    implicit none
    INTEGER,INTENT(IN),dimension(:)::j
    real(dp),INTENT(inOUT)::R1
    type (taylor),INTENT(IN)::S1
 !   integer k
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
    if(s1%i==0) call crap1("pek000  1" )  !call etall1(s1%i)
 !   k=s1%i
!    write(6,*) r1,k
    CALL DApek(s1%i,j,r1)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("pek000  2" ) ! newetall(s1%j,1)
    !
    !      CALL newDApek(s1%j,j,r1)
    !    endif
    !
  END SUBROUTINE pek000

  SUBROUTINE  pok000(S1,J,R1)
    implicit none
    INTEGER,INTENT(in),dimension(:)::j
    real(dp),INTENT(in)::R1
    type (taylor),INTENT(inout)::S1
    IF(.NOT.C_%STABLE_DA) RETURN

    if(check_j(j)/=0) return
    ! if(old) then
    if(s1%i==0) call crap1("pok000 1" )  ! call etall1(s1%i)
    CALL DApok(s1%i,j,r1)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("pok000  2" )  ! call newetall(s1%j,1)
    !
    !       CALL newDApok(s1%j,j,r1)
    !    endif
    !
  END SUBROUTINE pok000

  SUBROUTINE  TAYLOR_ran(S1,r1,R2)
    implicit none
    real(dp),INTENT(in)::R1
    real(dp),INTENT(inout)::R2
    type (taylor),INTENT(inout)::S1
    IF(.NOT.C_%STABLE_DA) RETURN

    !
    !     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
    !     FOR R1 > 0, THE VECTOR IS FILLED WITH REALS,
    !     FOR R1 < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
    !     ABS(R1) IS THE FILLING FACTOR
    ! if(old) then
    if(s1%i==0) call crap1("tAYLOR_ran  1" )  ! call etall1(s1%i)
    call daran(s1%i,r1,R2)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r))call crap1("tAYLOR_ran  2" )  !  call newetall(s1%j,1)
    !
    !       call newdaran(s1%j,r1,R2)
    !    endif
    !
  END SUBROUTINE TAYLOR_ran

  SUBROUTINE  intd_taylor(S1,S2,factor)
    implicit none
    type (taylor),INTENT(inOUT)::S2
    type (taylor),INTENT(IN)::S1(:)
    real(dp),INTENT(IN):: factor
    IF(.NOT.C_%STABLE_DA) RETURN

    ! if(old) then
    if(s1(1)%i==0) call crap1("intd_taylor 1")  !call etall1(s2%h%i)
    CALL intd(S1%i,s2%i,factor)
    !    else
    !       if(.NOT. ASSOCIATED(s1(1)%j%r)) call crap1("intd_taylor 2")  !call etall1(s2%h%i)
    !       CALL newintd(S1%j,s2%j,factor)
    !    endif
  END SUBROUTINE intd_taylor

  SUBROUTINE  DIFd_taylor(S2,S1,factor)
    implicit none
    type (taylor),INTENT(in)::S2
    type (taylor),INTENT(INOUT)::S1(:)
    real(dp),INTENT(IN):: factor
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
    CALL DIFD(S2%i,s1%i,factor)
    !    else
    !       CALL NEWDIFD(S2%j,s1%j,factor)
    !    endif
  END SUBROUTINE DIFd_taylor


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

  SUBROUTINE  CFUR(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    complex(dp) FUN
    EXTERNAL FUN
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
    if(s1%i==0) call crap1("CFUR  1" )  ! call etall1(s1%i)
    CALL DACFUR(s2%i,FUN,s1%i)
    !   else
    !      if(.NOT. ASSOCIATED(s1%j%r))call crap1("CFUR  2" )  !  call newetall(s1%j,1)
    !      CALL NEWDACFUR(s2%J,FUN,s1%J)
    !   endif

  END SUBROUTINE CFUR

  SUBROUTINE  CFUI(S2,FUN,S1)
    implicit none
    type (taylor),INTENT(INOUT)::S1
    type (taylor),INTENT(IN)::S2
    complex(dp) FUN
    EXTERNAL FUN
    IF(.NOT.C_%STABLE_DA) RETURN
    ! if(old) then
    if(s1%i==0)call crap1("CFUI  1" ) ! call etall1(s1%i)
    CALL DACFUI(s2%i,FUN,s1%i)
    !    else
    !       if(.NOT. ASSOCIATED(s1%j%r)) call crap1("CFUI  2" ) !call newetall(s1%j,1)
    !       CALL NEWDACFUI(s2%J,FUN,s1%J)
    !    endif

  END SUBROUTINE CFUI

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
             GETCHARnd2=0.0_dp
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
          ! call !write_e(0)
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
    nd2par=size(s2)
    ndel=0

    !frs get around compiler problem
    !frs    do i=1,len(trim(ADJUSTR (s2)))
    do i=1,nd2par
       Jfil(I)=s2(i)
       if(i>nv) then
          if(Jfil(i)>0) then
             GETintnd2=0.0_dp
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
          ! call !write_e(0)
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
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72),/))'
          w_p%c(1)=" error in GETintnd2t for .part_taylor. "
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

  END SUBROUTINE taylor_cycle


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
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72),/))'
          w_p%fi='(3((1X,i4)))'
          w_p%c(1)= "iass0user(master),scratchda(master)%n,newscheme_max"
          w_p=(/iass0user(master),scratchda(master)%n,newscheme_max/)
          ! call !write_e
          call ndum_warning_user
       endif
       iass0user(master)=0
    case(ndumt+1:)
       w_p=0
       w_p%nc=1
       w_p=(/"Should not be here"/)
       w_p%fc='(1((1X,A72),/))'
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
    do i=1,nd2par
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

  SUBROUTINE  print_for_bmad_parse(S1,MFILE,prec,ind)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(IN)::prec
    integer ,OPTIONAL,INTENT(IN)::ind
    type (TAYLOR),INTENT(IN)::S1

     bmadparser=1
    call pri(S1,MFILE,prec,ind)
     bmadparser=0
  
  end SUBROUTINE  print_for_bmad_parse

  SUBROUTINE  pri(S1,MFILE,prec,ind)
    implicit none
    INTEGER,INTENT(IN)::MFILE
    REAL(DP),OPTIONAL,INTENT(IN)::prec
    integer ,OPTIONAL,INTENT(IN)::ind
    type (TAYLOR),INTENT(IN)::S1
    REAL(DP) PREC1,depst,value
    integer i,j,it,n,indo,k,kt,kl
    integer, allocatable :: jc(:)
    character(255) line,line0

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
        write(mfile,*) line0(1:len_trim(line0))
       endif

    enddo

    deallocate(jc)


   else
   



    ! if(old) then
    if(print77) then
       CALL DAPRI77(s1%i,MFILE)
    else
       CALL DAPRI(s1%i,MFILE)
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
    INTEGER ipresent,k,n,I,illa
    real(dp) value
    INTEGER, allocatable :: j(:)
    call check_snake

    ! if(old) then
    if(s1%i==0)  call crap1("FILL_N 1")
    !    else
    !       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("FILL_N 2")
    !    endif


    IF(ASSOCIATED(S2%N)) S2=-1
    S2=0
    ipresent=1
    call dacycle(S1%I,ipresent,value,n)
    CALL ALLOC_U(S2,N,c_%nv)
    allocate(j(c_%nv))

    do i=1,N
       call dacycle(S1%I,i,value,illa,j)
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
    if(s1%i==0)  call crap1("REFILL_N 1")
    !    else
    !       IF (.NOT. ASSOCIATED(s1%j%r)) call crap1("REFILL_N 2")
    !    endif


    S1=0.0_dp

    IF(.not.ASSOCIATED(S2%N)) THEN
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72),/))'
       w_p%c(1)=" ERROR IN REFILL_N: UNIVERSAL_TAYLOR DOES NOT EXIST"
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
    integer                :: iunit
    integer                :: i,ii

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

    w_p=0
    w_p%nc=2
    w_p%fc='((1X,A72,/),(1X,A72))'
    w_p%c(1)= "ERROR IN :"
    w_p%c(2)= STRING
    ! call !write_e(3478)

  end subroutine crap1

  SUBROUTINE real_stop()
    implicit none
    integer i(1),j

    w_p=0
    w_p%nc=3
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


    w_p=0
    w_p%nc=3
    w_p%fc='(3((1X,A72),/))'
    w_p%c(1)=  " *****************************************************************"
    w_p%c(2)=  " *  Should never be here in New Linked List Scheme               *"
    w_p%c(3)=  " *****************************************************************"
    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A72),/))'
    w_p%c(1)= " do you want a crash? "
    ! call !write_e
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
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72)))'
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
    CALL ETALL1(DUMMY)
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
    CALL DADAL1(DUMMY)
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
       w_p=0
       w_p%nc=1
       w_p=(/" cannot indent anymore "/)
       w_p%fc='(1((1X,A72),/))'
       ! call !write_e(100)
       master=sqrt(-dble(master))
    end select
    !    write(26,*) "   taylor ",master
    call ass0(s1)

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
    ! if(old) then
    s1%i=scratchda(master)%PRESENT%T%i
    !    else
    !       s1%j=scratchda(master)%PRESENT%T%j
    !    endif


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
    ipresent=1
    call dacycle(S1%I,ipresent,value,n)

    allocate(j(c_%nv))

    do i=1,N
       call dacycle(S1%I,i,value,illa,j)
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

  !!! bessel

		FUNCTION I_nr(n,x)
        real(dp) I_nr
		INTEGER(I4B), INTENT(IN) :: n
		REAL(dp), INTENT(IN) :: x
        if(n>=2)  then
         I_nr=bessi_se(n,x)
!         I_nr=bessi(n,x)

        elseif(n==1) then
         I_nr=bessi1_se(x)
        elseif(n==0) then
         I_nr=bessi0_se(x)
        else
         write(6,*) "n<0 in I_n "
         stop
        endif
		END FUNCTION I_nr


		FUNCTION I_nt(n,x)
        type(taylor)  I_nt
		INTEGER(I4B), INTENT(IN) :: n
        type(taylor), INTENT(IN) :: x
        integer localmaster,i,j
        type(taylor)  dx,tx
        real(dp) fac,x0
!        real(dp) :: der(0:2*lno)
!        real(dp) :: dder(0:2*lno),ddert(0:2*lno)
         real(dp), allocatable :: ddert(:),dder(:),der(:)



        localmaster=master
        call ass(I_nt)
             x0=(x.sub.'0')


        if(c_%no==1) then
         I_nt=I_ns(n,x0)+dI_n(n,x0)*(x-x0)
         master=localmaster
         return
        endif
        allocate(der(0:c_%no+n),dder(0:c_%no+n),ddert(0:c_%no+n))

        der=0
        der(0)=I_ns(n,x0)            
 
        do i=n,c_%no+n
         der(i)=I_ns(i,x0)
        enddo
        j=max(0,n-c_%no)

        do i=n-1,j,-1
         der(i)=I_ns(i,x0)
        enddo

       call alloc(dx,tx)

          dder=0.0_dp
          dx=x-x0
          tx=dx
          I_nt=I_ns(n,x0)
          fac=1.0_dp
          dder(n)=1.0_dp

          do i=1,c_%no
            ddert=0.0_dp
          do j=max(n-(i-1),0),n+(i-1)
           ddert(iabs(j-1))=0.5_dp*dder(iabs(j))+ddert(iabs(j-1))
           ddert(iabs(j+1))=0.5_dp*dder(iabs(j))+ddert(iabs(j+1))
          enddo

            ddert=ddert/i
            fac=0
          do j=max(n-i,0),n+i
            fac=fac+ddert(j)*der(j)
          enddo
            I_nt=I_nt+fac*tx
            tx=tx*dx
            dder=ddert
          enddo  
       call kill(dx,tx)
        deallocate(der,dder,ddert)

        master=localmaster

		END FUNCTION I_nt

		FUNCTION In_nt(n,x)
        type(taylor)  in_nt
		INTEGER(I4B), INTENT(IN) :: n
        type(taylor), INTENT(IN) :: x
        integer localmaster,i,j
        type(taylor)  dx,tx
        real(dp) fac,x0
!        real(dp) :: der(0:2*lno)
!        real(dp) :: dder(0:2*lno),ddert(0:2*lno)
         real(dp), allocatable :: ddert(:),dder(:),der(:)



        localmaster=master
        call ass(in_nt)
             x0=(x.sub.'0')


        if(c_%no==1) then
         in_nt=I_n(n,x0)+dI_n(n,x0)*(x-x0)
         master=localmaster
         return
        endif
        allocate(der(0:c_%no+n),dder(0:c_%no+n),ddert(0:c_%no+n))

        der=0
        der(0)=I_n(n,x0)            
 
        do i=n,c_%no+n
         der(i)=I_n(i,x0)
        enddo
        j=max(0,n-c_%no)

        do i=n-1,j,-1
         der(i)=I_n(i,x0)
        enddo

       call alloc(dx,tx)

          dder=0.0_dp
          dx=x-x0
          tx=dx
          in_nt=I_n(n,x0)
          fac=1.0_dp
          dder(n)=1.0_dp

          do i=1,c_%no
            ddert=0.0_dp
          do j=max(n-(i-1),0),n+(i-1)
           ddert(iabs(j-1))=0.5_dp*dder(iabs(j))+ddert(iabs(j-1))
           ddert(iabs(j+1))=0.5_dp*dder(iabs(j))+ddert(iabs(j+1))
          enddo

            ddert=ddert/i
            fac=0
          do j=max(n-i,0),n+i
            fac=fac+ddert(j)*der(j)
          enddo
            in_nt=in_nt+fac*tx
            tx=tx*dx
            dder=ddert
          enddo  
       call kill(dx,tx)
        deallocate(der,dder,ddert)

        master=localmaster

		END FUNCTION In_nt

	FUNCTION bessi_se(n,x)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessi_se
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,m
	REAL(dp) :: bi,bim,bip,tox
	bessi_se=0.0
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
			bessi_se=scale(bessi_se,-IEXP)
			bi=scale(bi,-IEXP)
			bip=scale(bip,-IEXP)
		end if
		if (j == n) bessi_se=bip
	end do
	bessi_se=bessi_se*bessi0_se(x)/bi
	if (x < 0.0 .and. mod(n,2) == 1) bessi_se=-bessi_se
	END FUNCTION bessi_se

	FUNCTION bessi0_se(x)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessi0_se
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
		bessi0_se=poly_e(real((x/3.75_dp)**2,dp),p)
	else  
		bessi0_se=(exp(ax)/sqrt(ax))*poly_e(real(3.75_dp/ax,dp),q)
	end if
	END FUNCTION bessi0_se

	FUNCTION poly_e(x,coeffs)
! poly_rr
INTEGER(I4B), PARAMETER :: NPAR_POLY=8
	REAL(dp), INTENT(IN) :: x
	REAL(dp), DIMENSION(:), INTENT(IN) :: coeffs
	REAL(dp) :: poly_e
	REAL(dp) :: pow
	REAL(dp), DIMENSION(:), ALLOCATABLE :: vec
	INTEGER(I4B) :: i,n,nn
	n=size(coeffs)
	if (n <= 0) then
		poly_e=0.0_dp
	else if (n < NPAR_POLY) then
		poly_e=coeffs(n)
		do i=n-1,1,-1
			poly_e=x*poly_e+coeffs(i)
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
		poly_e=vec(1)
		deallocate(vec)
	end if
	END FUNCTION poly_e

	FUNCTION bessi1_se(x)
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: x
	REAL(dp) :: bessi1_se
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
		bessi1_se=ax*poly_e(real((x/3.75_dp)**2,dp),p)
	else
		bessi1_se=(exp(ax)/sqrt(ax))*poly_e(real(3.75_dp/ax,dp),q)
	end if
	if (x < 0.0) bessi1_se=-bessi1_se
	END FUNCTION bessi1_se

		FUNCTION dI_n(n,x)
        real(dp) dI_n
		INTEGER(I4B), INTENT(IN) :: n
		REAL(dp), INTENT(IN) :: x
        if(n>=1)  then
         dI_n=0.5_dp*(I_n(n+1,x)+I_n(n-1,x))
        elseif(n==0) then
         dI_n=bessi1_se(x)
        else
         write(6,*) "n<0 in dI_n "
         stop
        endif
		END FUNCTION dI_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Takuya OOURA  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bessel I_0(x) function in double precision
!
function dbesi0(x)
 implicit none
integer i,k
 real(dp) a(0 : 64), b(0 : 69), c(0 : 44),w,x,dbesi0,y,t
 data (a(i), i = 0, 12) /8.5246820682016865877d-11,2.5966600546497407288d-9,7.9689994568640180274d-8, &
 1.9906710409667748239d-6,4.0312469446528002532d-5,6.4499871606224265421d-4,7.9012345761930579108d-3, &
 7.1111111109207045212d-2,4.4444444444472490900d-1,1.7777777777777532045d0,4.0000000000000011182d0, &
 3.9999999999999999800d0,1.0000000000000000001d0 / 
 data (a(i), i = 13, 25) /1.1520919130377195927d-10, 2.2287613013610985225d-9,8.1903951930694585113d-8, &
 1.9821560631611544984d-6,4.0335461940910133184d-5, 6.4495330974432203401d-4,7.9013012611467520626d-3, &
 7.1111038160875566622d-2,4.4444450319062699316d-1, 1.7777777439146450067d0,4.0000000132337935071d0, &
 3.9999999968569015366d0,1.0000000003426703174d0 / 
 data (a(i), i = 26, 38) /1.5476870780515238488d-10, 1.2685004214732975355d-9,9.2776861851114223267d-8, &
 1.9063070109379044378d-6,4.0698004389917945832d-5, 6.4370447244298070713d-4,7.9044749458444976958d-3, &
 7.1105052411749363882d-2,4.4445280640924755082d-1, 1.7777694934432109713d0,4.0000055808824003386d0, &
 3.9999977081165740932d0,1.0000004333949319118d0 / 
 data (a(i), i = 39, 51) /2.0675200625006793075d-10,-6.1689554705125681442d-10,1.2436765915401571654d-7, &
 1.5830429403520613423d-6,4.2947227560776583326d-5, 6.3249861665073441312d-4,7.9454472840953930811d-3, &
 7.0994327785661860575d-2,4.4467219586283000332d-1,1.7774588182255374745d0,4.0003038986252717972d0, &
 3.9998233869142057195d0,1.0000472932961288324d0 / 
 data (a(i), i = 52, 64) /2.7475684794982708655d-10, -3.8991472076521332023d-9,1.9730170483976049388d-7, &
 5.9651531561967674521d-7,5.1992971474748995357d-5, 5.7327338675433770752d-4,8.2293143836530412024d-3, &
 6.9990934858728039037d-2,4.4726764292723985087d-1, 1.7726685170014087784d0,4.0062907863712704432d0, &
 3.9952750700487845355d0,1.0016354346654179322d0 / 
 data (b(i), i = 0, 13) /6.7852367144945531383d-8, 4.6266061382821826854d-7,6.9703135812354071774d-6, &
 7.6637663462953234134d-5,7.9113515222612691636d-4, 7.3401204731103808981d-3,6.0677114958668837046d-2, &
 4.3994941411651569622d-1,2.7420017097661750609d0, 14.289661921740860534d0,59.820609640320710779d0, &
 188.78998681199150629d0,399.87313678256011180d0, 427.56411572180478514d0 / 
 data (b(i), i = 14, 27) /1.8042097874891098754d-7, 1.2277164312044637357d-6,1.8484393221474274861d-5, &
 2.0293995900091309208d-4,2.0918539850246207459d-3, 1.9375315654033949297d-2, 1.5985869016767185908d-1, &
 1.1565260527420641724d0,7.1896341224206072113d0, 37.354773811947484532d0,155.80993164266268457d0, &
 489.52113711585409180d0,1030.9147225169564806d0, 1093.5883545113746958d0 / 
 data (b(i), i = 28, 41) /4.8017305613187493564d-7, 3.2613178439123800740d-6,4.9073137508166159639d-5, &
 5.3806506676487583755d-4,5.5387918291051866561d-3, 5.1223717488786549025d-2,4.2190298621367914765d-1, &
 3.0463625987357355872d0,18.895299447327733204d0, 97.915189029455461554d0,407.13940115493494659d0, &
 1274.3088990480582632d0,2670.9883037012547506d0, 2815.7166284662544712d0 / 
 data (b(i), i = 42, 55) /1.2789926338424623394d-6, 8.6718263067604918916d-6,1.3041508821299929489d-4, &
 1.4282247373727478920d-3,1.4684070635768789378d-2, 1.3561403190404185755d-1,1.1152592585977393953d0, &
 8.0387088559465389038d0,49.761318895895479206d0, 257.26842323135291380d0,1066.8543146269566231d0, &
 3328.3874581009636362d0,6948.8586598121634874d0, 7288.4893398212481055d0 / 
 data (b(i), i = 56, 69) /3.4093503681970328930d-6, 2.3079025203103376076d-5,3.4691373283901830239d-4, &
 3.7949949772229085450d-3,3.8974209677945602145d-2, 3.5949483804148783710d-1,2.9522878893539528226d0, &
 21.246564609514287056d0,131.28727387146173141d0, 677.38107093296675421d0,2802.3724744545046518d0, &
 8718.5731420798254081d0,18141.348781638832286d0, 18948.925349296308859d0 / 
 data (c(i), i = 0, 8) /2.5568678676452702768d-15, 3.0393953792305924324d-14,6.3343751991094840009d-13, &
 1.5041298011833009649d-11,4.4569436918556541414d-10, 1.7463930514271679510d-8,1.0059224011079852317d-6, &
 1.0729838945088577089d-4,5.1503226936425277380d-2 / 
 data (c(i), i = 9, 17) /5.2527963991711562216d-15, 7.2021184814210056410d-15,7.2561421229904797156d-13, &
 1.4823121466731042510d-11,4.4602670450376245434d-10, 1.7463600061788679671d-8,1.0059226091322347560d-6, &
 1.0729838937545111487d-4,5.1503226936437300716d-2 / 
 data (c(i), i = 18, 26) /1.3365917359358069908d-14, -1.2932643065888544835d-13,1.7450199447905602915d-12, &
 1.0419051209056979788d-11,4.5804788198059832600d-10, 1.7442405450073548966d-8,1.0059461453281292278d-6, &
 1.0729837434500161228d-4,5.1503226940658446941d-2 / 
 data (c(i), i = 27, 35) /5.3771611477352308649d-14, -1.1396193006413731702d-12,1.2858641335221653409d-11, &
 -5.9802086004570057703d-11,7.3666894305929510222d-10, 1.6731837150730356448d-8,1.0070831435812128922d-6, &
 1.0729733111203704813d-4,5.1503227360726294675d-2 / 
 data (c(i), i = 36, 44) /3.7819492084858931093d-14, -4.8600496888588034879d-13,1.6898350504817224909d-12, &
 4.5884624327524255865d-11,1.2521615963377513729d-10, 1.8959658437754727957d-8,1.0020716710561353622d-6, &
 1.0730371198569275590d-4,5.1503223833002307750d-2 / 
 w = abs(x)
 if (w .lt. 8.5d0) then
          t = w * w * 0.0625d0
          k = 13 * (int(t))
y = (((((((((((a(k) * t + a(k + 1)) * t +a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t +  &
a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t +  &
a(k + 11)) * t + a(k + 12)
      else if (w .lt. 12.5d0) then
          k = int(w)
          t = w - k
          k = 14 * (k - 8)
y = ((((((((((((b(k) * t + b(k + 1)) * t +b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t +  &
b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t +  &
b(k + 11)) * t + b(k + 12)) * t + b(k + 13)
      else
          t = 60 / w
          k = 9 * (int(t))
y = ((((((((c(k) * t + c(k + 1)) * t + c(k + 2)) * t + c(k + 3)) * t + c(k + 4)) * t +  &
c(k + 5)) * t + c(k + 6)) * t + c(k + 7)) * t + c(k + 8)) * sqrt(t) * exp(w)
      end if
      dbesi0 = y
      end function dbesi0



! Bessel I_1(x) function in double precision
!
function dbesi1(x)
 implicit none
 integer i
 real(dp) a(0 : 59), b(0 : 69), c(0 : 44),x,t,k,dbesi1,w,y
data (a(i), i = 0, 11) /1.2787464404046789181d-10, 3.5705860060088241077d-9,9.9611537619347335040d-8, 2.2395070088633043177d-6, &
4.0312466928887462346d-5, 5.6437387840203722356d-4, 5.9259259312934746096d-3, 4.4444444443499008870d-2, 2.2222222222232042719d-1, &
6.6666666666666139867d-1,1.0000000000000001106d0, 4.9999999999999999962d-1 / 
data (a(i), i = 12, 23) /1.7281952384448634449d-10, 3.0647204559976390130d-9,1.0237662138842827028d-7, 2.2299494417341498163d-6, & 
4.0335364374929326943d-5, 5.6433440269141349899d-4,5.9259754885893798654d-3, 4.4444399410880397870d-2, & 
2.2222225112835026730d-1, 6.6666665422146063244d-1,1.0000000032274936821d0, 4.9999999961866867205d-1 / 
data (a(i), i = 24, 35) /2.3216048939948030996d-10, 1.7443372702334489579d-9,1.1596478963485415499d-7, 2.1446755518623035147d-6, & 
4.0697440347437076195d-5, 5.6324394900433192204d-4,5.9283484996093060678d-3, 4.4440673899150997921d-2, & 
2.2222638016852657860d-1, 6.6666358151576732094d-1,1.0000013834029985337d0, 4.9999971643129650249d-1 / 
data (a(i), i = 36, 47) /3.1013758938255172562d-10, -8.4813676145611694984d-10,1.5544980187411802596d-7, 1.7811109378708045726d-6, & 
4.2945322199060856985d-5, 5.5344850176852353639d-4,5.9590327716950614802d-3, 4.4371611097707060659d-2, & 
2.2233578241986401111d-1, 6.6654747300463315310d-1,1.0000756505206705927d0, 4.9997803664415994554d-1 / 
data (a(i), i = 48, 59) /4.1214758313965020365d-10, -5.3613317735347429440d-9,2.4661360807517345161d-7, 6.7144593918926723203d-7, & 
5.1988027944945587571d-5, 5.0165568586065803067d-4,6.1717530047005289953d-3, 4.3745229577317251404d-2, & 
2.2363147971477747996d-1, 6.6475469131117660240d-1, 1.0015686689447547657d0, 4.9941120439785391891d-1 / 
data (b(i), i = 0, 13) /6.6324787943143095845d-8, 4.5125928898466638619d-7, 6.7937793134877246623d-6, 7.4580507871505926302d-5,  &
7.6866382927334005919d-4, 7.1185174803491859307d-3, 5.8721838073486424416d-2, 4.2473949281714196041d-1,  &
2.6396965606282079123d0, 13.710008536637016903d0, 57.158647688180932003d0, 179.46182892089389037d0,  &
377.57997362398478619d0, 399.87313678256009819d0 / 
data (b(i), i = 14, 27) /1.7652713206027939711d-7, 1.1988179244834708057d-6, 1.8037851545747139231d-5, 1.9775785516370314656d-4,  &
2.0354870702829387283d-3, 1.8822164191032253600d-2, 1.5500485219010424263d-1, 1.1190100010560573210d0,  &
6.9391565185406617552d0, 35.948170579648649345d0, 149.41909525103032616d0, 467.42979492780642582d0,  &
979.04227423171290408d0, 1030.9147225169564443d0 / 
data (b(i), i = 28, 41) /4.7022299276154507603d-7, 3.1878571710170115972d-6, 4.7940153875711448496d-5, 5.2496623508411440227d-4,  &
5.3968661134780824779d-3, 4.9837081920693776234d-2, 4.0979593830387765545d-1, 2.9533186922862948404d0,  &
18.278176130722516369d0, 94.476497150189121070d0, 391.66075612645333624d0, 1221.4182034643210345d0,  &
2548.6177980961291004d0, 2670.9883037012546541d0 / 
data (b(i), i = 42, 55) / 1.2535083724002034147d-6, 8.4845871420655708250d-6, 1.2753227372734042108d-4, 1.3950105363562648921d-3,  &
1.4325473993765291906d-2, 1.3212452778932829125d-1, 1.0849287786885151432d0, 7.8068089156260172673d0,  &
48.232254570679165833d0, 248.80659424902394371d0, 1029.0736929484210803d0, 3200.5629438795801652d0,  &
6656.7749162019607914d0, 6948.8586598121632302d0 / 
data (b(i), i = 56, 69) / 3.3439394490599745013d-6, 2.2600596902211837757d-5, 3.3955927589987356838d-4, 3.7105306061050972474d-3,  &
3.8065263634919156421d-2, 3.5068223415665236079d-1, 2.8760027832105027316d0, 20.665999500843274339d0,  &
127.47939148516390205d0, 656.43636874254000885d0, 2709.5242837932479920d0, 8407.1174233600734871d0,  &
17437.146284159740233d0, 18141.348781638831600d0 / 
data (c(i), i = 0, 8) / -2.8849790431465382128d-15, -3.5125350943844774657d-14, &
                        -7.4850867013707419750d-13, -1.8383904048277485153d-11, &
                        -5.7303556446977223342d-10, -2.4449502737311496525d-8, -1.6765373351766929724d-6, &
                        -3.2189516835265773471d-4, 5.1503226936425277377d-2/ 
data (c(i), i = 9, 17) /-5.8674306822281631119d-15, -9.4884898451194085565d-15, &
 -8.5033865136600364340d-13, -1.8142997866945285736d-11,  &
-5.7340238386338193949d-10, -2.4449138101742183665d-8, -1.6765375646678855842d-6, &
-3.2189516826945356325d-4,5.1503226936412017608d-2 / 
 data (c(i), i = 18, 26) /-1.4723362506764340882d-14, 1.3945147385179042899d-13, &
                          -1.9618041857586930923d-12, -1.3343606394065121821d-11, & 
                          -5.8649674606973244159d-10, -2.4426060539669553778d-8, -1.6765631828366988006d-6, &
                          -3.2189515191449587253d-4, 5.1503226931820146445d-2 / 
 data (c(i), i = 27, 35) / -5.8203519372580372987d-14, 1.2266326995309845825d-12, &
                           -1.3921625844526453237d-11, 6.2228025878281625469d-11, & 
                           -8.8636681342142794023d-10, -2.3661241616744818608d-8, &
	       -1.6777870960740520557d-6, -3.2189402882677074318d-4, 5.1503226479551959376d-2 / 
data (c(i), i = 36, 44) /-4.5801527369223291722d-14, 6.7998819697143727209d-13, &
                        -4.1624857909290468421d-12, -3.2849009406112440998d-11, & 
        -3.2478275690431118270d-10, -2.5739209934053714983d-8, -1.6730566573215739195d-6, &
         -3.2190010909008684076d-4, 5.1503229866932077150d-2 / 
      w = abs(x)
      if (w .lt. 8.5d0) then
          t = w * w * 0.0625d0
          k = 12 * (int(t))
y=(((((((((((a(k)*t+a(k+1))*t+a(k+2))*t+a(k+3))*t+a(k+4))*t+a(k+5))*t+a(k+6))*t+a(k+7))*t+a(k+8))*t+a(k+9))*t+a(k+10))*t+a(k+11))*w
elseif(w.lt.12.5d0)then
k=int(w)
t=w-k
k=14*(k-8)
y=((((((((((((b(k)*t+b(k+1))*t+b(k+2))*t+b(k+3))*t+b(k+4))*t+b(k+5))*t+b(k+6))*t+b(k+7))*t+b(k+8))*t+b(k+9))*t+b(k+10))*t+ &
b(k+11))*t+b(k+12))*t+b(k+13)
else
t=60/w
k=9*(int(t))
y=((((((((c(k)*t+c(k+1))*t+c(k+2))*t+c(k+3))*t+c(k+4))*t+c(k+5))*t+c(k+6))*t+c(k+7))*t+c(k + 8)) * sqrt(t) * exp(w)
      end if
      if (x .lt. 0) y = -y
      dbesi1 = y
end function dbesi1
!
function I_ent(n,x,nt)
 implicit none
 integer i,n,nt,m
 real(dp) I_ent,x,xh,xt,I_entb,eps,i0,di0,i1
logical doit



  if(x==0) then
   I_ent=0
   return
  endif

doit=.true.
eps=1.d-6

I_entb=i_ns(n,x)

 xh=x/2
 xt=1
 do i=1,n
  xt=xh*xt/i
 enddo 

 I_ent=xt
 xh=xh**2
 i0=I_ent
 di0=1.d38
 do m=1,nt
  xt=xt*xh/m/(m+n)
  I_ent= I_ent + xt
  i1=abs(i0-I_ent)
  i0=I_ent
!write(6,*) i0,i1,di0,doit
!pause 333
 if(doit) then
   if((abs(I_entb-I_ent)/abs(I_entb)) <eps) then
    doit=.false.
   endif 
  di0=i1
else
    if(i1>=di0) then
      return
     else
      di0=i1
    endif
endif
enddo

Write(6,*) " I_ent did not converge "
stop

end function I_ent


function In_en(n,x)
 implicit none
 integer i,n
 real(dp) In_en,x,i0,i1

 if(x==0) then
  In_en=0
  return
 endif
 if(n==0) then
  In_en=dbesi0(x)
   return
 elseif(n==1) then
  In_en=dbesi1(x)
  return
 endif

i0=dbesi0(x)
i1=dbesi1(x)

do i=2,n
 In_en=i0-2*(i-1)*i1/x
 i0=i1
 i1=In_en
enddo

end function In_en


function In_enz(n,x)
 implicit none
 integer i,n
 real(dp) In_enz,x,p0,p1,p,q,q0,q1

  if(x==0) then
   In_enz=0
   return
  endif

 if(n==0) then
  In_enz=dbesi0(x)
   return
 elseif(n==1) then
  In_enz=dbesi1(x)
  return
 endif

!i0=dbesi0(x)
!i1=dbesi1(x)

p0=1
p1=0
q0=0
q1=1

do i=2,n
 p=p0-2*(i-1)*p1/x
 q=q0-2*(i-1)*q1/x
 p0=p1
 p1=p
 q0=q1
 q1=q
enddo

p0=p*dbesi0(x)
q0=q*dbesi1(x)

q1=abs(p0+q0)/(abs(p0)+abs(q0))


if(q1>switch_bessel) then
  In_enz=p0+q0
else
  In_enz=I_ent(n,x,1000)
endif

end function In_enz

END MODULE  tpsa
