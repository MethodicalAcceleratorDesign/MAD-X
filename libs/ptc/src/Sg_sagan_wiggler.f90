!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN
module sagan_WIGGLER
 ! use S_def_all_kinds
  use S_DEF_KIND
  implicit none
  public
  private INTR,INTP,ZERO_SAGANr,ZERO_SAGANp
  private ALLOC_SAGAN,KILL_SAGAN,POINTERS_SAGANR,POINTERS_SAGANP
  private copy_el_elp ,copy_elp_el ,copy_el_el
  private scale_SAGANR,scale_SAGANP,PRINT_,READ_
  PRIVATE driftsaganR,driftsaganP,driftsagan
  PRIVATE KICKPATHR,KICKPATHP,KICKPATH
  PRIVATE COMPX_R,COMPX_P,COMPY_R,COMPY_P,COMPZ_R,COMPZ_P,b_fieldr,b_fieldp
 ! PRIVATE COMPZ  !,B_FIELD
  PRIVATE KICKR,KICKP,KICK
  PRIVATE KILL_WIGGLER,feval_saganp,feval_saganr,feval_sagan
  PRIVATE ALLOC_WIGGLER
  PRIVATE ZERO_Wr,ZERO_Wp,POINTERS_WP,e_potentialr,e_potentialp
  PRIVATE copy_W_WP,copy_WP_W,copy_W_W,INT_SAGANR,INT_SAGANP
  !  PRIVATE SET_R,SET_P,SET_W
  PRIVATE ADJUST_WIR,ADJUSTP_WI,get_z_wiR,get_z_wiP,kick_integral_r,kick_integral_p
  private ADJUST_like_abellr,ADJUST_like_abellp
  integer :: wiggler_sagan=6
   logical(lp) :: xprime_sagan=.false.,get_out=.false.
integer, parameter :: hyper_y_family_y = 1, hyper_xy_family_y = 2, hyper_x_family_y = 3
integer, parameter :: hyper_y_family_x = 4, hyper_xy_family_x = 5, hyper_x_family_x = 6
integer, parameter :: hyper_y_family_qu = 7, hyper_xy_family_qu = 8, hyper_x_family_qu = 9
integer, parameter :: hyper_y_family_sq = 10, hyper_xy_family_sq = 11, hyper_x_family_sq = 12
private conv_to_xprsagan,conv_to_xppsagan,conv_to_pxrsagan,conv_to_pxpsagan

  integer :: limit_sag(2) =(/4,18/) 
 
  INTERFACE conv_to_xp
     MODULE PROCEDURE conv_to_xprsagan
     MODULE PROCEDURE conv_to_xppsagan
  END INTERFACE

  INTERFACE conv_to_px
     MODULE PROCEDURE conv_to_pxrsagan
     MODULE PROCEDURE conv_to_pxpsagan
  END INTERFACE

  INTERFACE feval_sagan 
     MODULE PROCEDURE feval_saganr
     MODULE PROCEDURE feval_saganp
  END INTERFACE

  INTERFACE rk2_sagan
     MODULE PROCEDURE rk2saganr
     MODULE PROCEDURE rk2saganp
  END INTERFACE

  INTERFACE rk4_sagan
     MODULE PROCEDURE rk4saganr
     MODULE PROCEDURE rk4saganp
  END INTERFACE


  INTERFACE rk6_sagan
     MODULE PROCEDURE rk6saganr
     MODULE PROCEDURE rk6saganp
  END INTERFACE

  INTERFACE get_z_wi 
     MODULE PROCEDURE get_z_wiR
     MODULE PROCEDURE get_z_wip
  END INTERFACE

  INTERFACE driftsagan
     MODULE PROCEDURE driftsaganR
     MODULE PROCEDURE driftsaganP
  END INTERFACE

  INTERFACE KICKPATH
     MODULE PROCEDURE KICKPATHR
     MODULE PROCEDURE KICKPATHP
  END INTERFACE

  INTERFACE KICK
     MODULE PROCEDURE KICKR
     MODULE PROCEDURE KICKP
  END INTERFACE

  !  INTERFACE SET_W
  !     MODULE PROCEDURE SET_R
  !     MODULE PROCEDURE SET_P
  !  END INTERFACE

  INTERFACE COMPX
     MODULE PROCEDURE COMPX_R
     MODULE PROCEDURE COMPX_P
  END INTERFACE

  INTERFACE COMPY
     MODULE PROCEDURE COMPY_R
     MODULE PROCEDURE COMPY_P
  END INTERFACE

  INTERFACE COMPZ
     MODULE PROCEDURE COMPZ_R
     MODULE PROCEDURE COMPZ_P
  END INTERFACE

  INTERFACE B_FIELD
     MODULE PROCEDURE b_fieldr
     MODULE PROCEDURE b_fieldp
  END INTERFACE

  INTERFACE e_field
     MODULE PROCEDURE e_fieldr
     MODULE PROCEDURE e_fieldp
  END INTERFACE

  INTERFACE e_potential
     MODULE PROCEDURE e_potentialr
     MODULE PROCEDURE e_potentialp
  END INTERFACE

  INTERFACE TRACK
     MODULE PROCEDURE INTR
     MODULE PROCEDURE INTP
  END INTERFACE

  INTERFACE ALLOC
     MODULE PROCEDURE ALLOC_SAGAN
     MODULE PROCEDURE ALLOC_WIGGLER
  END INTERFACE

  INTERFACE POINTERS_SAGAN
     MODULE PROCEDURE POINTERS_SAGANR
     MODULE PROCEDURE POINTERS_SAGANP
  END INTERFACE

  INTERFACE POINTERS_W
     MODULE PROCEDURE POINTERS_WR
     MODULE PROCEDURE POINTERS_WP
  END INTERFACE


  INTERFACE KILL
     MODULE PROCEDURE KILL_SAGAN
     MODULE PROCEDURE KILL_WIGGLER
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_el_elp
     MODULE PROCEDURE copy_elp_el
     MODULE PROCEDURE copy_el_el
     MODULE PROCEDURE copy_W_WP
     MODULE PROCEDURE copy_WP_W
     MODULE PROCEDURE copy_W_W
  END INTERFACE

  INTERFACE scale_SAGAN
     MODULE PROCEDURE scale_SAGANR
     MODULE PROCEDURE scale_SAGANP
  END INTERFACE

  INTERFACE PRINT_USER
     MODULE PROCEDURE PRINT_
  END INTERFACE

  INTERFACE READ_USER
     MODULE PROCEDURE READ_
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE ZERO_SAGANr
     MODULE PROCEDURE ZERO_SAGANp
     MODULE PROCEDURE ZERO_Wr
     MODULE PROCEDURE ZERO_Wp
  END INTERFACE

  INTERFACE TRACK_SLICE
     MODULE PROCEDURE INT_SAGANR
     MODULE PROCEDURE INT_SAGANP
  END INTERFACE

  INTERFACE ADJUST_like_abell
     MODULE PROCEDURE ADJUST_like_abellr
     MODULE PROCEDURE ADJUST_like_abellp
  END INTERFACE

  INTERFACE ADJUST_WI
     MODULE PROCEDURE ADJUST_WIR
     MODULE PROCEDURE ADJUSTP_WI
  END INTERFACE

  INTERFACE kick_integral
     MODULE PROCEDURE kick_integral_r
     MODULE PROCEDURE kick_integral_p
  END INTERFACE


contains

  SUBROUTINE ADJUST_WIR(EL,X,k,J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(sagan),INTENT(INOUT):: EL
    TYPE(INTERNAL_STATE),OPTIONAL :: K    

    INTEGER, INTENT(IN) :: J


    IF(J==1.and.el%p%dir==-1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    elseIF(J==2.and.el%p%dir==1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    endif
  END SUBROUTINE ADJUST_WIR

  SUBROUTINE ADJUSTP_WI(EL,X,k,J)
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INOUT) :: X(6)
    TYPE(saganP),INTENT(INOUT):: EL
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    INTEGER, INTENT(IN) :: J

    IF(J==1.and.el%p%dir==-1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    elseIF(J==2.and.el%p%dir==1) then

     X(1)=X(1)-EL%INTERNAL(1)
     X(2)=X(2)-EL%INTERNAL(2)
     X(3)=X(3)-EL%INTERNAL(3)
     X(4)=X(4)-EL%INTERNAL(4)
     X(5)=X(5)-EL%INTERNAL(5)
     if(K%time) then
       X(6)=X(6)-EL%INTERNAL(6)/el%p%beta0
      else
       X(6)=X(6)-EL%INTERNAL(6)
     endif
    endif


  END SUBROUTINE ADJUSTP_WI


subroutine kick_integral_r(el,v,kx,ky,symp)
  implicit none    
  real(dp), INTENT(INOUT) :: v(6),kx,ky
  TYPE(sagan),INTENT(IN):: EL
  real(dp), pointer :: e(:)
  integer symp
  real(dp) Ix1,Ix2,Iy1,Iy2,x,y
  real(dp) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24
  real(dp) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24

! KYMA undulator field integrals on 9/8/2014
! from fitting I1y

   e=>el%w%ex

  a1 =     e(1)
  a2 =     e(2)
  a3 =     e(3)
  a4 =     e(4)
  a5 =     e(5)
  a6 =     e(6)
  a7 =     e(7)
  a8 =     e(8)
  a9 =     e(9)
  a10 =    e(10)
  a11 =    e(11)
  a12 =    e(12)
  a13 =    e(13)
  a14 =    e(14)
  a15 =    e(15)
  a16 =    e(16)
  a17 =    e(17)
  a18 =    e(18)
  a19 =    e(19)
  a20 =    e(20)
  a21 =    e(21)
  a22 =    e(22)
  a23 =    e(23)
  a24 =    e(24)

! from fitting I1x
 
  e=>el%w%ey  

  b1 =     e(1)
  b2 =     e(2)
  b3 =     e(3)
  b4 =     e(4)
  b5 =     e(5)
  b6 =     e(6)
  b7 =     e(7)
  b8 =     e(8)
  b9 =     e(9)
  b10 =    e(10)
  b11 =    e(11)
  b12 =    e(12)
  b13 =    e(13)
  b14 =    e(14)
  b15 =    e(15)
  b16 =    e(16)
  b17 =    e(17)
  b18 =    e(18)
  b19 =    e(19)
  b20 =    e(20)
  b21 =    e(21)
  b22 =    e(22)
  b23 =    e(23)
  b24 =    e(24)


if(symp==1) then
  Iy1=a1*sin(a2*v(1)+a3)*cosh(a2*v(3))+a4*sin(a5*v(1)+a6)*cosh(a5*v(3))+a7*sin(a8*v(1)+a9)*cosh(a8*v(3))+ &
           a10*sin(a11*v(1)+a12)*cosh(a11*v(3))+a13*sin(a14*v(1)+a15)*cosh(a14*v(3))+a16*sin(a17*v(1)+a18)*cosh(a17*v(3))+ &
           a19*sin(a20*v(1)+a21)*cosh(a20*v(3))+a22*sin(a23*v(1)+a24)*cosh(a23*v(3))

  Ix1=a1*cos(a2*v(1)+a3)*sinh(a2*v(3))+a4*cos(a5*v(1)+a6)*sinh(a5*v(3))+a7*cos(a8*v(1)+a9)*sinh(a8*v(3))+ &
           a10*cos(a11*v(1)+a12)*sinh(a11*v(3))+a13*cos(a14*v(1)+a15)*sinh(a14*v(3))+a16*cos(a17*v(1)+a18)*sinh(a17*v(3))+ &
           a19*cos(a20*v(1)+a21)*sinh(a20*v(3))+a22*cos(a23*v(1)+a24)*sinh(a23*v(3))

  Iy1=Iy1*1.0e-4_dp
  Ix1=Ix1*1.0e-4_dp



  Ix2=b1*sin(b2*v(1)+b3)*cosh(b2*v(3))+b4*sin(b5*v(1)+b6)*cosh(b5*v(3))+b7*sin(b8*v(1)+b9)*cosh(b8*v(3))+ &
           b10*sin(b11*v(1)+b12)*cosh(b11*v(3))+b13*sin(b14*v(1)+b15)*cosh(b14*v(3))+b16*sin(b17*v(1)+b18)*cosh(b17*v(3))+ &
           b19*sin(b20*v(1)+b21)*cosh(b20*v(3))+b22*sin(b23*v(1)+b24)*cosh(b23*v(3))
  Iy2=b1*cos(b2*v(1)+b3)*sinh(b2*v(3))+b4*cos(b5*v(1)+b6)*sinh(b5*v(3))+b7*cos(b8*v(1)+b9)*sinh(b8*v(3))+ &
           b10*cos(b11*v(1)+b12)*sinh(b11*v(3))+b13*cos(b14*v(1)+b15)*sinh(b14*v(3))+b16*cos(b17*v(1)+b18)*sinh(b17*v(3))+ &
           b19*cos(b20*v(1)+b21)*sinh(b20*v(3))+b22*cos(b23*v(1)+b24)*sinh(b23*v(3))

  Ix2=Ix2*1.0e-4_dp
  Iy2=Iy2*1.0e-4_dp

  kx=(-Iy1+Iy2)
  ky=(-Ix1+Ix2)
else

x=v(1)
y=v(3)

kx=b13*sinh(b14*y)*cos(b14*x)*cos(b15)-b22*sinh(b23*y)*sin(b23*x)*sin(b24)-b16*sinh(b17*y)*sin(b17*x)*sin(b18)
KX=KX-b19*sinh(b20*y)*sin(b20*x)*sin(b21)+b7*sinh(b8*y)*cos(b8*x)*cos(b9)-b13*sinh(b14*y)*sin(b14*x)*sin(b15)
KX=KX+b16*sinh(b17*y)*cos(b17*x)*cos(b18)+b19*sinh(b20*y)*cos(b20*x)*cos(b21)+b4*sinh(b5*y)*cos(b5*x)*cos(b6)
KX=KX+b1*sinh(b2*y)*cos(b2*x)*cos(b3)-b10*sinh(b11*y)*sin(b11*x)*sin(b12)-b7*sinh(b8*y)*sin(b8*x)*sin(b9)
KX=KX+b22*sinh(b23*y)*cos(b23*x)*cos(b24)+b10*sinh(b11*y)*cos(b11*x)*cos(b12)-b4*sinh(b5*y)*sin(b5*x)*sin(b6)
KX=KX-a10*sin(a11*x)*cos(a12)-a22*sin(a23*x)*cos(a24)-a7*sin(a8*x)*cos(a9)-a16*sin(a17*x)*cos(a18)-a13*sin(a14*x)*cos(a15)
KX=KX-a19*sin(a20*x)*cos(a21)-a1*sin(a2*x)*cos(a3)-a4*sin(a5*x)*cos(a6)-b1*sinh(b2*y)*sin(b2*x)*sin(b3)-a1*cos(a2*x)*sin(a3)
KX=KX-a4*cos(a5*x)*sin(a6)-a7*cos(a8*x)*sin(a9)-a10*cos(a11*x)*sin(a12)-a13*cos(a14*x)*sin(a15)-a16*cos(a17*x)*sin(a18)
KX=KX-a19*cos(a20*x)*sin(a21)-a22*cos(a23*x)*sin(a24)

ky=-a1*cos(a3)*sinh(a2*y)-a4*cos(a6)*sinh(a5*y)-a7*cos(a9)*sinh(a8*y)-a10*cos(a12)*sinh(a11*y)-a13*cos(a15)*sinh(a14*y)
KY=KY-a16*cos(a18)*sinh(a17*y)-a19*cos(a21)*sinh(a20*y)-a22*cos(a24)*sinh(a23*y)+b1*cosh(b2*y)*sin(b2*x)*cos(b3)
KY=KY+b1*cosh(b2*y)*cos(b2*x)*sin(b3)+b4*cosh(b5*y)*sin(b5*x)*cos(b6)+b4*cosh(b5*y)*cos(b5*x)*sin(b6)
KY=KY+b7*cosh(b8*y)*sin(b8*x)*cos(b9)
KY=KY+b7*cosh(b8*y)*cos(b8*x)*sin(b9)+b10*cosh(b11*y)*sin(b11*x)*cos(b12)+b10*cosh(b11*y)*cos(b11*x)*sin(b12)
KY=KY+b13*cosh(b14*y)*sin(b14*x)*cos(b15)+b13*cosh(b14*y)*cos(b14*x)*sin(b15)+b16*cosh(b17*y)*sin(b17*x)*cos(b18)
KY=KY+b16*cosh(b17*y)*cos(b17*x)*sin(b18)+b19*cosh(b20*y)*sin(b20*x)*cos(b21) 
KY=KY+b19*cosh(b20*y)*cos(b20*x)*sin(b21)+b22*cosh(b23*y)*sin(b23*x)*cos(b24)+b22*cosh(b23*y)*cos(b23*x)*sin(b24)

kx=kx*1.e-4_dp
ky=ky*1.e-4_dp
endif

end subroutine kick_integral_r


subroutine kick_integral_p(el,v,kx,ky,symp)
  implicit none    
  type(real_8), INTENT(INOUT) :: v(6),kx,ky
  TYPE(saganp),INTENT(IN):: EL
  real(dp), pointer :: e(:)
  integer symp
  type(real_8) Ix1,Ix2,Iy1,Iy2,x,y
  real(dp) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24
  real(dp) b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24

call ALLOC(Ix1,Ix2,Iy1,Iy2,x,y)

! KYMA undulator field integrals on 9/8/2014
! from fitting I1y

   e=>el%w%ex

  a1 =     e(1)
  a2 =     e(2)
  a3 =     e(3)
  a4 =     e(4)
  a5 =     e(5)
  a6 =     e(6)
  a7 =     e(7)
  a8 =     e(8)
  a9 =     e(9)
  a10 =    e(10)
  a11 =    e(11)
  a12 =    e(12)
  a13 =    e(13)
  a14 =    e(14)
  a15 =    e(15)
  a16 =    e(16)
  a17 =    e(17)
  a18 =    e(18)
  a19 =    e(19)
  a20 =    e(20)
  a21 =    e(21)
  a22 =    e(22)
  a23 =    e(23)
  a24 =    e(24)

! from fitting I1x
 
  e=>el%w%ey  

  b1 =     e(1)
  b2 =     e(2)
  b3 =     e(3)
  b4 =     e(4)
  b5 =     e(5)
  b6 =     e(6)
  b7 =     e(7)
  b8 =     e(8)
  b9 =     e(9)
  b10 =    e(10)
  b11 =    e(11)
  b12 =    e(12)
  b13 =    e(13)
  b14 =    e(14)
  b15 =    e(15)
  b16 =    e(16)
  b17 =    e(17)
  b18 =    e(18)
  b19 =    e(19)
  b20 =    e(20)
  b21 =    e(21)
  b22 =    e(22)
  b23 =    e(23)
  b24 =    e(24)


if(symp==1) then
  Iy1=a1*sin(a2*v(1)+a3)*cosh(a2*v(3))+a4*sin(a5*v(1)+a6)*cosh(a5*v(3))+a7*sin(a8*v(1)+a9)*cosh(a8*v(3))+ &
           a10*sin(a11*v(1)+a12)*cosh(a11*v(3))+a13*sin(a14*v(1)+a15)*cosh(a14*v(3))+a16*sin(a17*v(1)+a18)*cosh(a17*v(3))+ &
           a19*sin(a20*v(1)+a21)*cosh(a20*v(3))+a22*sin(a23*v(1)+a24)*cosh(a23*v(3))

  Ix1=a1*cos(a2*v(1)+a3)*sinh(a2*v(3))+a4*cos(a5*v(1)+a6)*sinh(a5*v(3))+a7*cos(a8*v(1)+a9)*sinh(a8*v(3))+ &
           a10*cos(a11*v(1)+a12)*sinh(a11*v(3))+a13*cos(a14*v(1)+a15)*sinh(a14*v(3))+a16*cos(a17*v(1)+a18)*sinh(a17*v(3))+ &
           a19*cos(a20*v(1)+a21)*sinh(a20*v(3))+a22*cos(a23*v(1)+a24)*sinh(a23*v(3))

  Iy1=Iy1*1.0e-4_dp
  Ix1=Ix1*1.0e-4_dp



  Ix2=b1*sin(b2*v(1)+b3)*cosh(b2*v(3))+b4*sin(b5*v(1)+b6)*cosh(b5*v(3))+b7*sin(b8*v(1)+b9)*cosh(b8*v(3))+ &
           b10*sin(b11*v(1)+b12)*cosh(b11*v(3))+b13*sin(b14*v(1)+b15)*cosh(b14*v(3))+b16*sin(b17*v(1)+b18)*cosh(b17*v(3))+ &
           b19*sin(b20*v(1)+b21)*cosh(b20*v(3))+b22*sin(b23*v(1)+b24)*cosh(b23*v(3))
  Iy2=b1*cos(b2*v(1)+b3)*sinh(b2*v(3))+b4*cos(b5*v(1)+b6)*sinh(b5*v(3))+b7*cos(b8*v(1)+b9)*sinh(b8*v(3))+ &
           b10*cos(b11*v(1)+b12)*sinh(b11*v(3))+b13*cos(b14*v(1)+b15)*sinh(b14*v(3))+b16*cos(b17*v(1)+b18)*sinh(b17*v(3))+ &
           b19*cos(b20*v(1)+b21)*sinh(b20*v(3))+b22*cos(b23*v(1)+b24)*sinh(b23*v(3))

  Ix2=Ix2*1.0e-4_dp
  Iy2=Iy2*1.0e-4_dp

  kx=(-Iy1+Iy2)
  ky=(-Ix1+Ix2)
else

x=v(1)
y=v(3)

!kx=b13*sinh(b14*y)*cos(b14*x)*cos(b15)-b22*sinh(b23*y)*sin(b23*x)*sin(b24)-b16*sinh(b17*y)*sin(b17*x)*sin(b18)-b19*sinh(b20*y) &
!*sin(b20*x)*sin(b21)+b7*sinh(b8*y)*cos(b8*x)*cos(b9)-b13*sinh(b14*y)*sin(b14*x)*sin(b15)+b16*sinh(b17*y)*cos(b17*x)*cos(b18)+ &
!b19*sinh(b20*y)*cos(b20*x)*cos(b21)+b4*sinh(b5*y)*cos(b5*x)*cos(b6)+b1*sinh(b2*y)*cos(b2*x)*cos(b3)-b10*sinh(b11*y)*sin(b11*x)* &
!sin(b12)-b7*sinh(b8*y)*sin(b8*x)*sin(b9)+b22*sinh(b23*y)*cos(b23*x)*cos(b24)+b10*sinh(b11*y)*cos(b11*x)*cos(b12)-b4*sinh(b5*y)* &
!sin(b5*x)*sin(b6)-a10*sin(a11*x)*cos(a12)-a22*sin(a23*x)*cos(a24)-a7*sin(a8*x)*cos(a9)-a16*sin(a17*x)*cos(a18)-a13*sin(a14*x)* &
!cos(a15)-a19*sin(a20*x)*cos(a21)-a1*sin(a2*x)*cos(a3)-a4*sin(a5*x)*cos(a6)-b1*sinh(b2*y)*sin(b2*x)*sin(b3)-a1*cos(a2*x)*sin(a3) &
!-a4*cos(a5*x)*sin(a6)-a7*cos(a8*x)*sin(a9)-a10*cos(a11*x)*sin(a12)-a13*cos(a14*x)*sin(a15)-a16*cos(a17*x)*sin(a18)-a19*cos(a20*x) &
!*sin(a21)-a22*cos(a23*x)*sin(a24)

!ky=-a1*cos(a3)*sinh(a2*y)-a4*cos(a6)*sinh(a5*y)-a7*cos(a9)*sinh(a8*y)-a10*cos(a12)*sinh(a11*y)-a13*cos(a15)*sinh(a14*y)-a16* &
!cos(a18)*sinh(a17*y)-a19*cos(a21)*sinh(a20*y)-a22*cos(a24)*sinh(a23*y)+b1*cosh(b2*y)*sin(b2*x)*cos(b3)+b1*cosh(b2*y)*cos(b2*x)* &
!sin(b3)+b4*cosh(b5*y)*sin(b5*x)*cos(b6)+b4*cosh(b5*y)*cos(b5*x)*sin(b6)+b7*cosh(b8*y)*sin(b8*x)*cos(b9)+b7*cosh(b8*y)*cos(b8*x) &
!*sin(b9)+b10*cosh(b11*y)*sin(b11*x)*cos(b12)+b10*cosh(b11*y)*cos(b11*x)*sin(b12)+b13*cosh(b14*y)*sin(b14*x)*cos(b15)+b13*cosh(b14*y) &
!*cos(b14*x)*sin(b15)+b16*cosh(b17*y)*sin(b17*x)*cos(b18)+b16*cosh(b17*y)*cos(b17*x)*sin(b18)+b19*cosh(b20*y)*sin(b20*x)*cos(b21) &
!+b19*cosh(b20*y)*cos(b20*x)*sin(b21)+b22*cosh(b23*y)*sin(b23*x)*cos(b24)+b22*cosh(b23*y)*cos(b23*x)*sin(b24)

kx=b13*sinh(b14*y)*cos(b14*x)*cos(b15)-b22*sinh(b23*y)*sin(b23*x)*sin(b24)-b16*sinh(b17*y)*sin(b17*x)*sin(b18)
KX=KX-b19*sinh(b20*y)*sin(b20*x)*sin(b21)+b7*sinh(b8*y)*cos(b8*x)*cos(b9)-b13*sinh(b14*y)*sin(b14*x)*sin(b15)
KX=KX+b16*sinh(b17*y)*cos(b17*x)*cos(b18)+b19*sinh(b20*y)*cos(b20*x)*cos(b21)+b4*sinh(b5*y)*cos(b5*x)*cos(b6)
KX=KX+b1*sinh(b2*y)*cos(b2*x)*cos(b3)-b10*sinh(b11*y)*sin(b11*x)*sin(b12)-b7*sinh(b8*y)*sin(b8*x)*sin(b9)
KX=KX+b22*sinh(b23*y)*cos(b23*x)*cos(b24)+b10*sinh(b11*y)*cos(b11*x)*cos(b12)-b4*sinh(b5*y)*sin(b5*x)*sin(b6)
KX=KX-a10*sin(a11*x)*cos(a12)-a22*sin(a23*x)*cos(a24)-a7*sin(a8*x)*cos(a9)-a16*sin(a17*x)*cos(a18)-a13*sin(a14*x)*cos(a15)
KX=KX-a19*sin(a20*x)*cos(a21)-a1*sin(a2*x)*cos(a3)-a4*sin(a5*x)*cos(a6)-b1*sinh(b2*y)*sin(b2*x)*sin(b3)-a1*cos(a2*x)*sin(a3)
KX=KX-a4*cos(a5*x)*sin(a6)-a7*cos(a8*x)*sin(a9)-a10*cos(a11*x)*sin(a12)-a13*cos(a14*x)*sin(a15)-a16*cos(a17*x)*sin(a18)
KX=KX-a19*cos(a20*x)*sin(a21)-a22*cos(a23*x)*sin(a24)

ky=-a1*cos(a3)*sinh(a2*y)-a4*cos(a6)*sinh(a5*y)-a7*cos(a9)*sinh(a8*y)-a10*cos(a12)*sinh(a11*y)-a13*cos(a15)*sinh(a14*y)
KY=KY-a16*cos(a18)*sinh(a17*y)-a19*cos(a21)*sinh(a20*y)-a22*cos(a24)*sinh(a23*y)+b1*cosh(b2*y)*sin(b2*x)*cos(b3)
KY=KY+b1*cosh(b2*y)*cos(b2*x)*sin(b3)+b4*cosh(b5*y)*sin(b5*x)*cos(b6)+b4*cosh(b5*y)*cos(b5*x)*sin(b6)
KY=KY+b7*cosh(b8*y)*sin(b8*x)*cos(b9)
KY=KY+b7*cosh(b8*y)*cos(b8*x)*sin(b9)+b10*cosh(b11*y)*sin(b11*x)*cos(b12)+b10*cosh(b11*y)*cos(b11*x)*sin(b12)
KY=KY+b13*cosh(b14*y)*sin(b14*x)*cos(b15)+b13*cosh(b14*y)*cos(b14*x)*sin(b15)+b16*cosh(b17*y)*sin(b17*x)*cos(b18)
KY=KY+b16*cosh(b17*y)*cos(b17*x)*sin(b18)+b19*cosh(b20*y)*sin(b20*x)*cos(b21) 
KY=KY+b19*cosh(b20*y)*cos(b20*x)*sin(b21)+b22*cosh(b23*y)*sin(b23*x)*cos(b24)+b22*cosh(b23*y)*cos(b23*x)*sin(b24)


kx=kx*1.e-4_dp
ky=ky*1.e-4_dp
endif

call kill(Ix1,Ix2,Iy1,Iy2,x,y)

end subroutine kick_integral_p



  SUBROUTINE INTR(EL,X,k,mid)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(INOUT):: EL
    TYPE(WORM),OPTIONAL,INTENT(INOUT):: mid
    INTEGER I,ENT,EXI
    TYPE(INTERNAL_STATE),OPTIONAL :: K
 
    IF(EL%P%DIR==1) THEN
       ENT=1;EXI=2;
    ELSE
       ENT=2;EXI=1;
    ENDIF

    IF(.NOT.PRESENT(MID))call ADJUST_LIKE_ABELL(EL,X,k,ENT)

    IF(PRESENT(MID)) CALL XMID(MID,X,0)

    DO I=1,EL%P%NST
       call track_slice(el,x,k,i)
       IF(PRESENT(MID)) CALL XMID(MID,X,i)
    ENDDO

   IF(.NOT.PRESENT(MID))call ADJUST_LIKE_ABELL(EL,X,k,EXI)

    call ADJUST_WI(EL,X,k,2)

  END SUBROUTINE INTR

  SUBROUTINE get_z_wir(EL,i,z)
    IMPLICIT NONE
    TYPE(SAGAN),INTENT(INOUT):: EL
    integer i
    real(dp),INTENT(INOUT):: z
    real(dp) d

    D=EL%L/EL%P%NST
    IF(EL%P%DIR==1) THEN
       Z=(i-1)*d
    ELSE
       Z=EL%L-(i-1)*d
    ENDIF

  end SUBROUTINE get_z_wir

  SUBROUTINE get_z_wip(EL,i,z)
    IMPLICIT NONE
    TYPE(SAGANP),INTENT(INOUT):: EL
    integer i
    TYPE(REAL_8),INTENT(INOUT):: z
    TYPE(REAL_8) d

    CALL ALLOC(D)

    D=EL%L/EL%P%NST
    IF(EL%P%DIR==1) THEN
       Z=(i-1)*d
    ELSE
       Z=EL%L-(i-1)*d
    ENDIF

    CALL KILL(D)

  end SUBROUTINE get_z_wip

  SUBROUTINE INT_SAGANR(EL,X,k,i)
    IMPLICIT NONE
    integer ipause, mypause
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(INOUT):: EL
    integer,INTENT(IN):: I
    real(dp) Z
    real(dp) D,DH
    real(dp) D1,D2,DK1,DK2
    real(dp) DF(4),DK(4)
    INTEGER J
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    !    CALL SET_W(EL%W)

    SELECT CASE(EL%P%METHOD)
    CASE(2)



       DH=EL%L/2.0_dp/EL%P%NST
       D=EL%L/EL%P%NST
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       if(el%xprime) then
        call rk2_sagan(z,d,el,x,k)
       else

          Z=Z+EL%P%DIR*DH
          CALL driftsagan(EL,DH,Z,1,X,k)
          CALL driftsagan(EL,DH,Z,2,X,k)
          CALL KICKPATH(EL,DH,X,k)
          CALL KICK(EL,D,Z,X,k)
          CALL KICKPATH(EL,DH,X,k)
          CALL driftsagan(EL,DH,Z,2,X,k)
          CALL driftsagan(EL,DH,Z,1,X,k)
       !       Z=Z+EL%P%DIR*DH
      endif
    CASE(4)
       D=EL%L/EL%P%NST

       DK1=D*FK1
       D1=DK1/2.0_dp
       DK2=D*FK2
       D2=DK2/2.0_dp
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       if(el%xprime) then
        call rk4_sagan(z,d,el,x,k)
       else

       Z=Z+EL%P%DIR*D1
       CALL driftsagan(EL,D1,Z,1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL driftsagan(EL,D1,Z,1,X,k)
       Z=Z+EL%P%DIR*D1+D2
       CALL driftsagan(EL,D2,Z,1,X,k)
       CALL driftsagan(EL,D2,Z,2,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL KICK(EL,DK2,Z,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL driftsagan(EL,D2,Z,2,X,k)
       CALL driftsagan(EL,D2,Z,1,X,k)
       Z=Z+EL%P%DIR*(D1+D2)
       CALL driftsagan(EL,D1,Z,1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL driftsagan(EL,D1,Z,1,X,k)
      endif


    CASE(6)
       DO j =1,4
          DK(j)=EL%L*YOSK(J)/EL%P%NST
          DF(j)=DK(j)/2.0_dp
       ENDDO
       D=EL%L/EL%P%NST
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       if(el%xprime) then
        call rk6_sagan(z,d,el,x,k)
       else


       DO J=4,1,-1
          Z=Z+EL%P%DIR*DF(J)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO
       DO J=2,4
          Z=Z+EL%P%DIR*DF(J)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO

      endif



    CASE DEFAULT
       WRITE(6,*) " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       ipause=mypause(357)
    END SELECT


  END SUBROUTINE INT_SAGANR



  SUBROUTINE INTP(EL,X,k)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(INOUT):: EL
    TYPE(INTERNAL_STATE),OPTIONAL :: K
    INTEGER I,ENT,EXI
 
    IF(EL%P%DIR==1) THEN
       ENT=1;EXI=2;
    ELSE
       ENT=2;EXI=1;
    ENDIF

    call ADJUST_LIKE_ABELL(EL,X,k,ENT)

    DO I=1,EL%P%NST
       call track_slice(el,x,k,i)
    ENDDO

    call ADJUST_LIKE_ABELL(EL,X,k,EXI)

    call ADJUST_WI(EL,X,k,2)

  END SUBROUTINE INTP

  SUBROUTINE INT_SAGANP(EL,X,k,i)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(INOUT):: EL
    integer,INTENT(IN):: I
    TYPE(REAL_8) Z
    TYPE(REAL_8) D,DH
    TYPE(REAL_8) D1,D2,DK1,DK2
    TYPE(REAL_8) DF(4),DK(4)
    INTEGER J
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    CALL ALLOC(Z,D,DH,D1,D2,DK1,DK2)
    CALL ALLOC(DF,4)
    CALL ALLOC(DK,4)

    !    CALL SET_W(EL%W)
    SELECT CASE(EL%P%METHOD)
    CASE(2)
       DH=EL%L/2.0_dp/EL%P%NST
       D=EL%L/EL%P%NST

       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       if(el%xprime) then
        call rk2_sagan(z,d,el,x,k)
       else

          Z=Z+EL%P%DIR*DH
          CALL driftsagan(EL,DH,Z,1,X,k)
          CALL driftsagan(EL,DH,Z,2,X,k)
          CALL KICKPATH(EL,DH,X,k)
          CALL KICK(EL,D,Z,X,k)
          CALL KICKPATH(EL,DH,X,k)
          CALL driftsagan(EL,DH,Z,2,X,k)
          CALL driftsagan(EL,DH,Z,1,X,k)
       !       Z=Z+EL%P%DIR*DH
      endif

    CASE(4)
       D=EL%L/EL%P%NST

       DK1=D*FK1
       D1=DK1/2.0_dp
       DK2=D*FK2
       D2=DK2/2.0_dp

       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       if(el%xprime) then
        call rk4_sagan(z,d,el,x,k)
       else

       Z=Z+EL%P%DIR*D1
       CALL driftsagan(EL,D1,Z,1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL driftsagan(EL,D1,Z,1,X,k)
       Z=Z+EL%P%DIR*D1+D2
       CALL driftsagan(EL,D2,Z,1,X,k)
       CALL driftsagan(EL,D2,Z,2,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL KICK(EL,DK2,Z,X,k)
       CALL KICKPATH(EL,D2,X,k)
       CALL driftsagan(EL,D2,Z,2,X,k)
       CALL driftsagan(EL,D2,Z,1,X,k)
       Z=Z+EL%P%DIR*(D1+D2)
       CALL driftsagan(EL,D1,Z,1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL KICK(EL,DK1,Z,X,k)
       CALL KICKPATH(EL,D1,X,k)
       CALL driftsagan(EL,D1,Z,2,X,k)
       CALL driftsagan(EL,D1,Z,1,X,k)
      endif


    CASE(6)
       DO j =1,4
          DK(j)=EL%L*YOSK(J)/EL%P%NST
          DF(j)=DK(j)/2.0_dp
       ENDDO
       D=EL%L/EL%P%NST
       IF(EL%P%DIR==1) THEN
          Z=(i-1)*d
       ELSE
          Z=EL%L-(i-1)*d
       ENDIF

       if(el%xprime) then
        call rk6_sagan(z,d,el,x,k)
       else


       DO J=4,1,-1
          Z=Z+EL%P%DIR*DF(J)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO
       DO J=2,4
          Z=Z+EL%P%DIR*DF(J)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL KICK(EL,DK(J),Z,X,k)
          CALL KICKPATH(EL,DF(J),X,k)
          CALL driftsagan(EL,DF(J),Z,2,X,k)
          CALL driftsagan(EL,DF(J),Z,1,X,k)
          Z=Z+EL%P%DIR*DF(J)
       ENDDO

      endif

    CASE DEFAULT
       WRITE(6,*) " THE METHOD ",EL%P%METHOD," IS NOT SUPPORTED"
       ipause=mypause(357)
    END SELECT

    CALL KILL(Z,D,DH,D1,D2,DK1,DK2)
    CALL KILL(DF,4)
    CALL KILL(DK,4)



  END SUBROUTINE INT_SAGANP

  SUBROUTINE ZERO_SAGANr(EL,I)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%INTERNAL))  THEN
          DEALLOCATE(EL%INTERNAL)
          DEALLOCATE(EL%n_min)
          DEALLOCATE(EL%xprime)
          EL%W=-1
          DEALLOCATE(EL%W)
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%INTERNAL)
       NULLIFY(EL%n_min)
       NULLIFY(EL%xprime)
       NULLIFY(EL%W)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZERO_SAGANr

  SUBROUTINE ZERO_SAGANp(EL,I)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED KILL AND DEASSOCIATE
       IF(ASSOCIATED(EL%INTERNAL))  THEN
          CALL KILL(EL)              ! FPP DEALLOCATION FIRST OBVIOUSLY
          EL%W=-1
          DEALLOCATE(EL%INTERNAL)
          DEALLOCATE(EL%n_min)
          DEALLOCATE(EL%W)
          DEALLOCATE(EL%xprime)
       ENDIF
    elseif(i==0)       then
       NULLIFY(EL%W)
       NULLIFY(EL%n_min)
       NULLIFY(EL%INTERNAL)  !!! was not there January 2014
       NULLIFY(EL%xprime)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
       !
    endif

  END SUBROUTINE ZERO_SAGANp

  SUBROUTINE ZERO_Wr(EL,I)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED DEASSOCIATE
       IF(ASSOCIATED(EL%K))  THEN
          DEALLOCATE(EL%A)
          DEALLOCATE(EL%F,EL%x0,EL%Y0)
          DEALLOCATE(EL%FORM)
          DEALLOCATE(EL%K)
       ENDIF
       IF(ASSOCIATED(EL%KE))  THEN
          DEALLOCATE(EL%AE)
          DEALLOCATE(EL%FE,EL%x0E,EL%Y0E)
          DEALLOCATE(EL%FORME)
          DEALLOCATE(EL%KE)
       ENDIF
         IF(ASSOCIATED(EL%offset))  DEALLOCATE(EL%offset)
         IF(ASSOCIATED(EL%ex))  DEALLOCATE(EL%ex,EL%ey)
         IF(ASSOCIATED(EL%n))  DEALLOCATE(EL%n)
         IF(ASSOCIATED(EL%ne))  DEALLOCATE(EL%ne)
    elseif(i==0)       then
          NULLIFY(EL%A,EL%AE)
          NULLIFY(EL%F,EL%x0,EL%y0,EL%y0,EL%x0e,EL%y0e)
          NULLIFY(EL%offset)
          NULLIFY(EL%FORM,EL%FORME)
          NULLIFY(EL%K,EL%KE)
          NULLIFY(EL%ex,EL%ey)
          NULLIFY(EL%n,EL%ne)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
    endif

  END SUBROUTINE ZERO_Wr

  SUBROUTINE ZERO_Wp(EL,I)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(inout)::EL
    INTEGER, INTENT(IN)::I
    IF(I==-1) THEN
       !  Set real(dp) variables to zero or whatever  if any
       !  IF POINTED ASSOCIATED KILL AND DEASSOCIATE
          IF(ASSOCIATED(EL%K).OR.ASSOCIATED(EL%KE)) CALL KILL(EL)  
       IF(ASSOCIATED(EL%K))  THEN
          DEALLOCATE(EL%A)
          DEALLOCATE(EL%F,EL%x0,EL%Y0)
          DEALLOCATE(EL%FORM)
          DEALLOCATE(EL%K)
       ENDIF
       IF(ASSOCIATED(EL%KE))  THEN
          DEALLOCATE(EL%AE)
          DEALLOCATE(EL%FE,EL%x0E,EL%Y0E)
          DEALLOCATE(EL%FORME)
          DEALLOCATE(EL%KE)
       ENDIF
         IF(ASSOCIATED(EL%offset))  DEALLOCATE(EL%offset)
         IF(ASSOCIATED(EL%ex))  DEALLOCATE(EL%ex,EL%ey)
         IF(ASSOCIATED(EL%n))  DEALLOCATE(EL%n)
         IF(ASSOCIATED(EL%ne))  DEALLOCATE(EL%ne)

    elseif(i==0)       then
          NULLIFY(EL%A,EL%AE)
          NULLIFY(EL%F,EL%x0,EL%y0,EL%x0e,EL%y0e)
          NULLIFY(EL%offset)
          NULLIFY(EL%FORM,EL%FORME)
          NULLIFY(EL%K,EL%KE)
          NULLIFY(EL%ex,EL%ey)
          NULLIFY(EL%n,EL%ne)
       ! nullifies pointers
       ! And also zeroes for security ordinary variables
       !
    endif

  END SUBROUTINE ZERO_Wp


  SUBROUTINE copy_el_elp(EL,ELP)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(in)::EL
    TYPE(SAGANP), INTENT(inout)::ELP

    ELP%INTERNAL(1)    =EL%INTERNAL(1)
    ELP%INTERNAL(2)    =EL%INTERNAL(2)
    ELP%INTERNAL(3)    =EL%INTERNAL(3)
    ELP%INTERNAL(4)    =EL%INTERNAL(4)
    ELP%INTERNAL(5)    =EL%INTERNAL(5)
    ELP%INTERNAL(6)    =EL%INTERNAL(6)
    ELP%n_min    =EL%n_min
    ELP%xprime   =EL%xprime
    CALL COPY(EL%W,ELP%W)
    !  COPY CODING HERE NO ALLOCATION OF POINTERS OR POLYMORPH NEEDED
    !  IF DONE CORRECTLY

  END SUBROUTINE copy_el_elp

  SUBROUTINE copy_elp_el(EL,ELP)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(in)::EL
    TYPE(SAGAN), INTENT(inout)::ELP

    ELP%INTERNAL(1)    =EL%INTERNAL(1)
    ELP%INTERNAL(2)    =EL%INTERNAL(2)
    ELP%INTERNAL(3)    =EL%INTERNAL(3)
    ELP%INTERNAL(4)    =EL%INTERNAL(4)
    ELP%INTERNAL(5)    =EL%INTERNAL(5)
    ELP%INTERNAL(6)    =EL%INTERNAL(6)
    ELP%n_min    =EL%n_min
    ELP%xprime   =EL%xprime
    CALL COPY(EL%W,ELP%W)

  END SUBROUTINE copy_elp_el

  SUBROUTINE copy_el_el(EL,ELP)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(in)::EL
    TYPE(SAGAN), INTENT(inout)::ELP

    ELP%INTERNAL(1)    =EL%INTERNAL(1)
    ELP%INTERNAL(2)    =EL%INTERNAL(2)
    ELP%INTERNAL(3)    =EL%INTERNAL(3)
    ELP%INTERNAL(4)    =EL%INTERNAL(4)
    ELP%INTERNAL(5)    =EL%INTERNAL(5)
    ELP%INTERNAL(6)    =EL%INTERNAL(6)
    ELP%n_min    =EL%n_min
    ELP%xprime   =EL%xprime
    !  COPY CODING HERE NO ALLOCATION OF POINTERS
    CALL COPY(EL%W,ELP%W)

  END SUBROUTINE copy_el_el



  SUBROUTINE copy_W_W(EL,ELP)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(in)::EL
    TYPE(undu_r), INTENT(inout)::ELP

    INTEGER I,J,n,ne
    n=0
    ne=0
    if(associated(el%K)) n=EL%n
    if(associated(el%Ke)) ne=EL%ne

    if(associated(el%K).or.associated(el%Ke)) CALL POINTERS_W(ELP,n,ne)

    if(associated(el%K)) then
 !      CALL POINTERS_W(ELP,SIZE(EL%A),SIZE(EL%AE))


       DO I=1,3
          DO J=1,SIZE(EL%A)
             ELP%K(I,J)   =EL%K(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%A)
          ELP%A(I)    =EL%A(I)
          ELP%F(I)    =EL%F(I)
          ELP%x0(I)    =EL%x0(I)
          ELP%y0(I)    =EL%y0(I)
          ELP%FORM(I) =EL%FORM(I)
       ENDDO


     ENDIF

if(associated(el%ex))        elp%ex=el%ex
if(associated(el%ey))    elp%ey=el%ey
if(associated(EL%offset))    ELP%offset   =EL%offset
if(associated(EL%n))    ELP%n   =EL%n
if(associated(EL%ne))    ELP%ne   =EL%ne

    if(associated(el%KE)) then
       DO I=1,3
          DO J=1,SIZE(EL%AE)
             ELP%KE(I,J)   =EL%KE(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%AE)
          ELP%AE(I)    =EL%AE(I)
          ELP%FE(I)    =EL%FE(I)
          ELP%x0E(I)    =EL%x0E(I)
          ELP%y0E(I)    =EL%y0E(I)
          ELP%FORME(I) =EL%FORME(I)
       ENDDO

    endif

  END SUBROUTINE copy_W_W

  SUBROUTINE copy_W_WP(EL,ELP)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(in)::EL
    TYPE(undu_p), INTENT(inout)::ELP

    INTEGER I,J,n,ne
    n=0
    ne=0
    if(associated(el%K)) n=EL%n
    if(associated(el%Ke)) ne=EL%ne

    if(associated(el%K).or.associated(el%Ke)) CALL POINTERS_W(ELP,n,ne)

    if(associated(el%K)) then
!       CALL POINTERS_W(ELP,SIZE(EL%A),SIZE(EL%AE))
       DO I=1,3
          DO J=1,SIZE(EL%A)
             ELP%K(I,J)   =EL%K(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%A)
          ELP%A(I)    =EL%A(I)
          ELP%F(I)    =EL%F(I)
          ELP%x0(I)    =EL%x0(I)
          ELP%y0(I)    =EL%y0(I)
          ELP%FORM(I) =EL%FORM(I)
       ENDDO
 ENDIF
if(associated(el%ex))        elp%ex=el%ex
if(associated(el%ey))    elp%ey=el%ey
if(associated(EL%offset))    ELP%offset   =EL%offset
if(associated(EL%n))    ELP%n   =EL%n
if(associated(EL%ne))    ELP%ne   =EL%ne
    if(associated(el%KE)) then
       DO I=1,3
          DO J=1,SIZE(EL%AE)
             ELP%KE(I,J)   =EL%KE(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%AE)
          ELP%AE(I)    =EL%AE(I)
          ELP%FE(I)    =EL%FE(I)
          ELP%x0E(I)    =EL%x0E(I)
          ELP%y0E(I)    =EL%y0E(I)
          ELP%FORME(I) =EL%FORME(I)
       ENDDO

    endif

  END SUBROUTINE copy_W_WP

  SUBROUTINE copy_WP_W(EL,ELP)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(in)::EL
    TYPE(undu_r), INTENT(inout)::ELP

    INTEGER I,J,n,ne
    n=0
    ne=0
    if(associated(el%K)) n=EL%n
    if(associated(el%Ke)) ne=EL%ne

    if(associated(el%K).or.associated(el%Ke)) CALL POINTERS_W(ELP,n,ne)

    if(associated(el%K)) then
!       CALL POINTERS_W(ELP,SIZE(EL%A),SIZE(EL%AE))
       DO I=1,3
          DO J=1,SIZE(EL%A)
             ELP%K(I,J)   =EL%K(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%A)
          ELP%A(I)    =EL%A(I)
          ELP%x0(I)    =EL%x0(I)
          ELP%y0(I)    =EL%y0(I)
          ELP%F(I)    =EL%F(I)
          ELP%FORM(I) =EL%FORM(I)
       ENDDO

ENDIF
if(associated(el%ex))        elp%ex=el%ex
if(associated(el%ey))    elp%ey=el%ey
if(associated(EL%offset))    ELP%offset   =EL%offset
if(associated(EL%n))    ELP%n   =EL%n
if(associated(EL%ne))    ELP%ne   =EL%ne
    if(associated(el%KE)) then
       DO I=1,3
          DO J=1,SIZE(EL%AE)
             ELP%KE(I,J)   =EL%KE(I,J)
          ENDDO
       ENDDO
       DO I=1,SIZE(EL%AE)
          ELP%AE(I)    =EL%AE(I)
          ELP%FE(I)    =EL%FE(I)
          ELP%x0E(I)    =EL%x0E(I)
          ELP%y0E(I)    =EL%y0E(I)
          ELP%FORME(I) =EL%FORME(I)
       ENDDO

    endif

  END SUBROUTINE copy_WP_W

  SUBROUTINE POINTERS_SAGANR(EL)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(INOUT)::EL

    ALLOCATE(EL%INTERNAL(6))
    EL%INTERNAL=0.0_dp
    allocate(EL%n_min,EL%xprime)
     EL%n_min=wiggler_sagan
     EL%xprime   =xprime_sagan
    ALLOCATE(EL%W)
    !CALL POINTERS_W(EL%W)
    el%w=0

    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_SAGANR

  SUBROUTINE POINTERS_SAGANP(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL

    ALLOCATE(EL%INTERNAL(6))
    allocate(EL%n_min,EL%xprime)
     EL%n_min=wiggler_sagan
     EL%xprime   =xprime_sagan
    ALLOCATE(EL%W)
    !CALL POINTERS_W(EL%W)
    el%w=0
    ! ALLOCATE INTERNAL POINTERS IF ANY

  END SUBROUTINE POINTERS_SAGANP

  SUBROUTINE POINTERS_WR(EL,N,NE)
    IMPLICIT NONE
    TYPE(undu_r), INTENT(INOUT)::EL
    INTEGER, INTENT(IN)::N,ne

!    IF(ASSOCIATED(EL%A).OR.ASSOCIATED(EL%AE)) THEN
    IF(ASSOCIATED(EL%offset)) THEN
       EL=-1
    ENDIF
    EL=0

allocate(el%n)
el%n=n
allocate(el%ne)
el%ne=ne

IF(N>0) THEN
    ALLOCATE(EL%A(N))
    ALLOCATE(EL%F(N),EL%x0(N),EL%y0(N))
    ALLOCATE(EL%FORM(N))
    ALLOCATE(EL%K(3,N))

    EL%K=0.0_dp
    EL%x0=0.0_dp
    EL%y0=0.0_dp
    EL%A=0.0_dp
    EL%F=0.0_dp
    EL%FORM=0
ENDIF
    ALLOCATE(EL%offset)
    EL%offset=0.0_dp
IF(NE>0) THEN
    ALLOCATE(EL%AE(NE))
    ALLOCATE(EL%FE(NE),EL%x0E(NE),EL%y0E(NE))
    ALLOCATE(EL%FORME(NE))
    ALLOCATE(EL%KE(3,NE))
    EL%KE=0.0_dp
    EL%x0E=0.0_dp
    EL%y0E=0.0_dp
    EL%AE=0.0_dp
    EL%FE=0.0_dp
    EL%FORME=0
ENDIF
    ALLOCATE(EL%ex(wiggler_suntao),EL%ey(wiggler_suntao))
    EL%ex=0
    EL%ey=0
  END SUBROUTINE POINTERS_WR

  SUBROUTINE POINTERS_WP(EL,N,ne)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER, INTENT(IN)::N,ne

!    IF(ASSOCIATED(EL%A).OR.ASSOCIATED(EL%AE)) THEN
    IF(ASSOCIATED(EL%offset)) THEN
       CALL KILL(EL)
       EL=-1
    ENDIF
    EL=0

allocate(el%n)
el%n=n
allocate(el%ne)
el%ne=ne

IF(N>0) THEN
    ALLOCATE(EL%A(N))
    ALLOCATE(EL%F(N),EL%x0(N),EL%y0(N))
    ALLOCATE(EL%FORM(N))
    ALLOCATE(EL%K(3,N))

ENDIF
    ALLOCATE(EL%offset)
IF(NE>0) THEN
    ALLOCATE(EL%AE(NE))
    ALLOCATE(EL%FE(NE),EL%x0E(NE),EL%y0E(NE))
    ALLOCATE(EL%FORME(NE))
    ALLOCATE(EL%KE(3,NE))
ENDIF

 
!    CALL ALLOC(EL)
    ALLOCATE(EL%ex(wiggler_suntao),EL%ey(wiggler_suntao))
    EL%ex=0
    EL%ey=0
    ! ALLOCATE INTERNAL POINTERS IF ANY
    CALL ALLOC(EL)
  END SUBROUTINE POINTERS_WP



  SUBROUTINE ALLOC_SAGAN(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL
    CALL ALLOC(EL%INTERNAL,6)
 
    ! CALL ALLOC(EL%W)
    ! ALLOC INTERNAL POLYMORPHS IF ANY
  END SUBROUTINE ALLOC_SAGAN

  SUBROUTINE ALLOC_WIGGLER(EL)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER I,J
    ! ALLOC INTERNAL POLYMORPHS IF ANY
    IF(ASSOCIATED(EL%K)) THEN  !DAVID
       DO I=1,3
          DO J=1,SIZE(EL%A)
             CALL ALLOC(EL%K(I,J));
          ENDDO
       ENDDO
       CALL ALLOC(EL%A,SIZE(EL%A));
       CALL ALLOC(EL%F,SIZE(EL%A));
       CALL ALLOC(EL%x0,SIZE(EL%A));
       CALL ALLOC(EL%y0,SIZE(EL%A));
ENDIF
      IF(ASSOCIATED(EL%offset)) CALL ALLOC(EL%offset);
    IF(ASSOCIATED(EL%KE)) THEN  
       DO I=1,3
          DO J=1,SIZE(EL%AE)
             CALL ALLOC(EL%KE(I,J));
          ENDDO
       ENDDO
       CALL ALLOC(EL%AE,SIZE(EL%AE));
       CALL ALLOC(EL%FE,SIZE(EL%AE));
       CALL ALLOC(EL%x0E,SIZE(EL%AE));
       CALL ALLOC(EL%y0E,SIZE(EL%AE));
    ENDIF
  END SUBROUTINE ALLOC_WIGGLER


  SUBROUTINE KILL_SAGAN(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL

    CALL KILL(EL%INTERNAL,6)
    CALL KILL(EL%W)
    ! KILL INTERNAL POLYMORPHS IF ANY

  END SUBROUTINE KILL_SAGAN

  SUBROUTINE KILL_WIGGLER(EL)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER I,J
    ! KILL INTERNAL POLYMORPHS IF ANY
    IF(ASSOCIATED(EL%K)) then! DAVID
       DO I=1,3
          DO J=1,SIZE(EL%A)
             CALL KILL(EL%K(I,J));
          ENDDO
       ENDDO
       CALL KILL(EL%A,SIZE(EL%A));
       CALL KILL(EL%x0,SIZE(EL%A));
       CALL KILL(EL%y0,SIZE(EL%A));
       CALL KILL(EL%F,SIZE(EL%A));
ENDIF

    IF(ASSOCIATED(EL%KE)) then! DAVID
       DO I=1,3
          DO J=1,SIZE(EL%AE)
             CALL KILL(EL%KE(I,J));
          ENDDO
       ENDDO
       CALL KILL(EL%AE,SIZE(EL%AE));
       CALL KILL(EL%x0E,SIZE(EL%AE));
       CALL KILL(EL%y0E,SIZE(EL%AE));
       CALL KILL(EL%FE,SIZE(EL%AE));

    ENDIF
  END SUBROUTINE KILL_WIGGLER

  SUBROUTINE reset_WI(EL)
    IMPLICIT NONE
    TYPE(SAGANP), INTENT(INOUT)::EL

    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS

    CALL resetpoly_R31N(EL%INTERNAL,6)
    CALL reset_WIG(EL%W)

  END SUBROUTINE reset_WI

  SUBROUTINE reset_WIG(EL)
    IMPLICIT NONE
    TYPE(undu_p), INTENT(INOUT)::EL
    INTEGER I,J
    IF(ASSOCIATED(EL%K)) THEN ! DAVID
       DO I=1,3
          DO J=1,SIZE(EL%A)
             CALL resetpoly_R31(EL%K(I,J));
          ENDDO
       ENDDO

       CALL resetpoly_R31N(EL%A,SIZE(EL%A))
       CALL resetpoly_R31N(EL%F,SIZE(EL%A))
       CALL resetpoly_R31N(EL%x0,SIZE(EL%A))
       CALL resetpoly_R31N(EL%y0,SIZE(EL%A))
       CALL resetpoly_R31(EL%offset)
ENDIF
    IF(ASSOCIATED(EL%KE)) THEN ! DAVID
       DO I=1,3
          DO J=1,SIZE(EL%AE)
             CALL resetpoly_R31(EL%KE(I,J));
          ENDDO
       ENDDO

       CALL resetpoly_R31N(EL%AE,SIZE(EL%AE))
       CALL resetpoly_R31N(EL%FE,SIZE(EL%AE))
       CALL resetpoly_R31N(EL%x0E,SIZE(EL%AE))
       CALL resetpoly_R31N(EL%y0E,SIZE(EL%AE))

    ENDIF
  END SUBROUTINE reset_WIG


  SUBROUTINE  ELp_POL_SAGAN(S2,S2R,S1,DONEIT)
    implicit none
    integer ipause,mypause
    type (POL_BLOCK),INTENT(IN):: S1
    TYPE(SAGANp),INTENT(inOUT):: S2
    TYPE(SAGAN),INTENT(inOUT):: S2R
    LOGICAL(lp),INTENT(inOUT)::  DONEIT
    integer i
!!! NOT DONE ETIENNE
    ! ONE CAN LINK INTERNAL POLYMORPHS TO PART OF POL_BLOCK WHICH IS NOT USED
    ! HERE THE VARIABLE "INTERNAL" IS LINKED TO VOLT
    ! We also linked it to the pol_block2


    !   IF(S1%IVOLT>0) THEN
    !      s2%INTERNAL%I=S1%IVOLT+S1%NPARA
    !      s2%INTERNAL%S=S1%SVOLT
    !      s2%INTERNAL%KIND=3
    !!      DONEIT=.TRUE.
    !      IF(S1%SET_TPSAFIT) THEN
    !         s2%INTERNAL%R=s2%INTERNAL%R+s2%INTERNAL%S*s1%TPSAFIT(S1%IVOLT)
    !      ENDIF
    !    ENDIF
    !  or try
    DO I=1,6
       IF(S1%SAGAN%Iinternal(I)>0) THEN
          s2%INTERNAL(I)%I=S1%SAGAN%Iinternal(I)+S1%NPARA
          s2%INTERNAL(I)%S=S1%SAGAN%Sinternal(I)
          s2%INTERNAL(I)%KIND=3
          if(S1%SAGAN%Iinternal(I)>c_%np_pol) c_%np_pol=S1%SAGAN%Iinternal(I)
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN

             s2%INTERNAL(I)%R=s2%INTERNAL(I)%R+s2%INTERNAL(I)%S*s1%TPSAFIT(S1%SAGAN%Iinternal(I))
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2R%INTERNAL(I)=s2%INTERNAL(I)%R
          ENDIF
       ENDIF
    ENDDO
    if(size(s2%w%A)>size(S1%SAGAN%w%ia)) then
       write(6,*) " Pol_block for wiggler must be made bigger ",size(s2%w%A),size(S1%SAGAN%w%ia)
       ipause=mypause(121)
    endif

    do i=1,SIZE(s2%w%A)
       IF(S1%SAGAN%w%ia(i)>0) THEN
          s2%w%a(i)%I=S1%SAGAN%w%ia(i)+S1%NPARA
          s2%w%a(i)%S=S1%SAGAN%w%sa(i)
          s2%w%a(i)%KIND=3
          if(S1%SAGAN%w%ia(i)>c_%np_pol) c_%np_pol=S1%SAGAN%w%ia(i)
          DONEIT=.TRUE.
          IF(S1%SET_TPSAFIT) THEN
             s2%w%a(i)%R=s2%w%a(i)%R+s2%w%a(i)%S*s1%TPSAFIT(S1%SAGAN%w%ia(i))
          ENDIF
          IF(S1%SET_ELEMENT) THEN
             s2R%w%a(i)=s2%w%a(i)%R
          ENDIF
       ENDIF
    enddo

  end SUBROUTINE  ELp_POL_SAGAN

  SUBROUTINE  scale_SAGANR(S2,P0C_OLD,P0C_NEW,power)
    implicit none
    TYPE(SAGAN),INTENT(inOUT):: S2
    real(dp),INTENT(IN)::  P0C_OLD,P0C_NEW
    integer, INTENT(IN):: power
    INTEGER I
    ! EXAMPLE

    !    S2%INTERNAL= S2%INTERNAL*P0C_OLD/P0C_NEW
    if(associated(S2%W%A)) then
     DO I=1,SIZE(S2%W%A)
        S2%W%A(I)=S2%W%A(I)*(P0C_OLD/P0C_NEW)**power
     ENDDO
    endif

   if(associated(S2%W%AE)) then
     DO I=1,SIZE(S2%W%AE)
        S2%W%AE(I)=S2%W%AE(I)*(P0C_OLD/P0C_NEW)**power
     ENDDO
    endif

    S2%W%ex(1:wiggler_suntao)=S2%W%ex(1:wiggler_suntao)*(P0C_OLD/P0C_NEW)**power
    S2%W%ey(1:wiggler_suntao)=S2%W%ey(1:wiggler_suntao)*(P0C_OLD/P0C_NEW)**power

  end SUBROUTINE  scale_SAGANR

  SUBROUTINE  scale_SAGANP(S2,P0C_OLD,P0C_NEW,power)
    implicit none
    TYPE(SAGANp),INTENT(inOUT):: S2
    real(dp),INTENT(IN)::  P0C_OLD,P0C_NEW
    integer, INTENT(IN):: power
    INTEGER I
    ! EXAMPLE

    !    S2%INTERNAL= S2%INTERNAL*P0C_OLD/P0C_NEW
    if(associated(S2%W%A)) then
     DO I=1,SIZE(S2%W%A)
        S2%W%A(I)=S2%W%A(I)*(P0C_OLD/P0C_NEW)**power
     ENDDO
    endif

   if(associated(S2%W%AE)) then
     DO I=1,SIZE(S2%W%AE)
        S2%W%AE(I)=S2%W%AE(I)*(P0C_OLD/P0C_NEW)**power
     ENDDO
    endif

    S2%W%ex(1:wiggler_suntao)=S2%W%ex(1:wiggler_suntao)*(P0C_OLD/P0C_NEW)**power
    S2%W%ey(1:wiggler_suntao)=S2%W%ey(1:wiggler_suntao)*(P0C_OLD/P0C_NEW)**power

  end SUBROUTINE  scale_SAGANP

  ! split driftsagans
  SUBROUTINE driftsaganR(EL,L,Z,PLANE,X,k)
    IMPLICIT NONE
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(INOUT):: X(6)
    real(dp), INTENT(IN):: L,Z
    INTEGER, INTENT(IN)::PLANE
    real(dp) PZ,A,B,AP,BP
    TYPE(INTERNAL_STATE),OPTIONAL :: K


    IF(PLANE==1) THEN
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)-A
       X(4)=X(4)-AP
       if(k%TIME) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(1)=X(1)+L*X(2)/pz
          X(6)=X(6)+((X(2)*X(2))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(1)=X(1)+L*X(2)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(2)*X(2))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)+A
       X(4)=X(4)+AP
    ELSE
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)-BP
       X(4)=X(4)-B
       if(k%TIME) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(3)=X(3)+L*X(4)/pz
          X(6)=X(6)+((X(4)*X(4))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(3)=X(3)+L*X(4)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)+BP
       X(4)=X(4)+B
    ENDIF
  END SUBROUTINE driftsaganR

  SUBROUTINE driftsaganP(EL,L,Z,PLANE,X,k)
    IMPLICIT NONE
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(REAL_8), INTENT(IN):: L,Z
    INTEGER, INTENT(IN)::PLANE
    TYPE(REAL_8) PZ,A,B,AP,BP
    TYPE(INTERNAL_STATE),OPTIONAL :: K

    CALL ALLOC(PZ,A,B,AP,BP)
    IF(PLANE==1) THEN
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)-A
       X(4)=X(4)-AP
       if(k%TIME) then
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(1)=X(1)+L*X(2)/pz
          X(6)=X(6)+((X(2)*X(2))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(1)=X(1)+L*X(2)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(2)*X(2))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPX(EL,Z,X,A,AP)
       X(2)=X(2)+A
       X(4)=X(4)+AP
    ELSE
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)-BP
       X(4)=X(4)-B
       if(k%TIME) then
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)
          X(3)=X(3)+L*X(4)/pz
          X(6)=X(6)+((X(4)*X(4))/2.0_dp/pz**2)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          X(3)=X(3)+L*X(4)/(1.0_dp+X(5))
          X(6)=X(6)+(L/(1.0_dp+X(5)))*(X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL  COMPY(EL,Z,X,B,BP)
       X(2)=X(2)+BP
       X(4)=X(4)+B
    ENDIF

    CALL KILL(PZ,A,B,AP,BP)

  END SUBROUTINE driftsaganP

  SUBROUTINE KICKPATHR(EL,L,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: L
    real(dp) PZ,PZ0,DPZ
    TYPE(INTERNAL_STATE),OPTIONAL :: K
    ! ETIENNE
    IF(EL%P%EXACT) THEN
       if(k%TIME) then
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2-X(2)**2-X(4)**2)
          PZ0=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp/EL%P%BETA0+X(5))/PZ +k%TOTALPATH*L/EL%P%BETA0
          PZ=ROOT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)   ! WRONG NOT SYMPLECTIC 2015.8.11
          X(6)=X(6)-((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          PZ=ROOT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          PZ0=1.0_dp+X(5)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp+X(5))/PZ +k%TOTALPATH*L 
          PZ=ROOT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          X(6)=X(6)-(L/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
    ELSE
       if(k%TIME) then
          X(6)=X(6)+k%TOTALPATH*L/EL%P%BETA0
       else
          X(6)=X(6)+k%TOTALPATH*L
       endif
    ENDIF

  END SUBROUTINE KICKPATHR

  SUBROUTINE KICKPATHP(EL,L,X,k)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: L
    TYPE(REAL_8) PZ,PZ0,DPZ
    TYPE(INTERNAL_STATE),OPTIONAL :: K
    ! ETIENNE
    IF(EL%P%EXACT) THEN
       CALL ALLOC(PZ,PZ0,DPZ)
       if(k%TIME) then
          PZ=sqrt(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2-X(2)**2-X(4)**2)
          PZ0=sqrt(1.0_dp+2.0_dp*X(5)/EL%P%beta0+x(5)**2)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp/EL%P%BETA0+X(5))/PZ +k%TOTALPATH*L/EL%P%BETA0
          PZ=SQRT(1.0_dp+2.0_dp*X(5)/EL%P%BETA0+x(5)**2)   ! WRONG NOT SYMPLECTIC 2015.8.11
          X(6)=X(6)-((X(2)*X(2)+X(4)*X(4))/2.0_dp/pz**2+1.0_dp)*(1.0_dp/EL%P%BETA0+x(5))*L/pz
       else
          PZ=sqrt((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          PZ0=1.0_dp+X(5)
          DPZ=(X(2)**2+X(4)**2)/PZ/PZ0/(PZ+PZ0)   ! = (one/PZ-one/PZ0)
          X(1)=X(1)+L*X(2)*DPZ
          X(3)=X(3)+L*X(4)*DPZ
          X(6)=X(6)+L*(1.0_dp+X(5))/PZ +k%TOTALPATH*L 
          PZ=SQRT((1.0_dp+X(5))**2-X(2)**2-X(4)**2)
          X(6)=X(6)-(L/(1.0_dp+X(5)))*(X(2)*X(2)+X(4)*X(4))/2.0_dp/(1.0_dp+X(5))
       endif
       CALL KILL(PZ,PZ0,DPZ)
    ELSE
       if(k%TIME) then
          X(6)=X(6)+k%TOTALPATH*L/EL%P%BETA0
       else
          X(6)=X(6)+k%TOTALPATH*L
       endif
    ENDIF

  END SUBROUTINE KICKPATHP

  !   X_PLANE
  SUBROUTINE COMPX_R(EL,Z,X,A,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp


   DO I=1,el%w%n  !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =   EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))/2.0_dp)**2*EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2*0.5_dp &
                *COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
!          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          A =  EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))*0.5_dp)**2 *EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2* 0.5_dp &
               *COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          A =  EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
!          B =  EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =  EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
!          A =  EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          A =    EL%W%A(I)*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =   0.5_dp*EL%W%A(I)*sinx_x(0.5d0*EL%W%K(1,i)*(X(1)+EL%W%X0(i)))**2*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(1,i)*EL%W%K(2,i)/EL%W%K(3,i)*(X(1)+EL%W%X0(i))**2 + B
    
       elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          A =    EL%W%A(I)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(3,i)**2 + A
          B =   0.5_dp*EL%W%K(1,i)*EL%W%K(2,i)**2*EL%W%A(I)*sinhx_x(0.5_dp*EL%W%K(1,i)*(X(1)+EL%W%X0(i)))**2   &
                *sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i)**2*(X(1)+EL%W%X0(i))**2 + B

       elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          A =    EL%W%A(I)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(1,i)/EL%W%K(3,i) + A

          B =   -EL%W%A(I)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)**2/EL%W%K(1,i)**2/EL%W%K(3,i) + B

       elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          A =    EL%W%A(I)*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A

          B =   EL%W%A(I)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(3,i)*(X(1)+EL%W%X0(i)) + B
        elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          A =    EL%W%A(I)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(3,i)**2 + A

          B =   EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*(X(1)+EL%W%X0(i))*EL%W%K(2,i)**2/EL%W%K(3,i)**2 + B
        elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          A =    EL%W%A(I)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(1,i)/EL%W%K(3,i) + A

          B =   EL%W%A(I)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)**2/EL%W%K(1,i)**2/EL%W%K(3,i) + B   
           
       endif
     enddo
 
    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 
    
  END SUBROUTINE COMPX_R

  SUBROUTINE COMPX_P(EL,Z,X,A,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp
! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X

   DO I=1,el%w%n ! SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =   EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))/2.0_dp)**2*EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2*0.5_dp &
                *cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
!          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          A =  EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i))*0.5_dp)**2 *EL%W%K(1,i)*(X(1)+EL%W%X0(i))**2* 0.5_dp &
               *cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          A =  EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
!          B =  EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =  EL%W%A(I)*EL%W%K(3,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
!          A =  EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          A =    EL%W%A(I)*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =   0.5_dp*EL%W%A(I)*sinx_x(0.5d0*EL%W%K(1,i)*(X(1)+EL%W%X0(i)))**2*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(1,i)*EL%W%K(2,i)/EL%W%K(3,i)*(X(1)+EL%W%X0(i))**2 + B
       elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          A =    EL%W%A(I)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(3,i)**2 + A
          B =   0.5_dp*EL%W%K(1,i)*EL%W%K(2,i)**2*EL%W%A(I)*sinhx_x(0.5_dp*EL%W%K(1,i)*(X(1)+EL%W%X0(i)))**2* &
                sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i)**2*(X(1)+EL%W%X0(i))**2 + B
       elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          A =    EL%W%A(I)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(1,i)/EL%W%K(3,i) + A

          B =   -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)**2/EL%W%K(1,i)**2/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          A =    EL%W%A(I)*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A

          B =   EL%W%A(I)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(3,i)*(X(1)+EL%W%X0(i)) + B
        elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          A =    EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(3,i)**2 + A

          B =   EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*(X(1)+EL%W%X0(i))*EL%W%K(2,i)**2/EL%W%K(3,i)**2 + B
        elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          A =    EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)/EL%W%K(1,i)/EL%W%K(3,i) + A

          B =   EL%W%A(I)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))*EL%W%K(2,i)**2/EL%W%K(1,i)**2/EL%W%K(3,i) + B      
       endif
     enddo

    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 
  END SUBROUTINE COMPX_P

  !   Y_PLANE
  SUBROUTINE COMPY_R(EL,Z,X,A,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp
 
    DO I=1,el%w%n  !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i))  &
               *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + B
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          A =  -EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2 *EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) & 
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2*EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          A =   -EL%W%A(I)*EL%W%K(1,i)*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(3,i) + A
          B =    EL%W%A(I)*EL%W%K(1,i)**2*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))) &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          A =   -EL%W%A(I)*EL%W%K(1,i)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i)**2 + A
          B =   -0.5_DP*EL%W%A(I)*EL%W%K(1,i)**2*EL%W%K(2,i)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*Sinhx_x(EL%W%K(2,i) &
                *(X(3)+EL%W%Y0(I))*0.5_DP )**2 &
                *sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I))**2/EL%W%K(3,i)**2 + B
       elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          A =   -EL%W%A(I)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =   -EL%W%A(I)*EL%W%K(1,i)*EL%W%K(2,i)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                 *Sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_DP )**2 &
                *sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I))**2*0.5_DP/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          A =   EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i)/EL%W%K(3,i))*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                *coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A

          B =   EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i))**2*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          A =   -EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i)**2)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I)) + A

          B =   -EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))**2*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))* &
                 sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I)) + B
       elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          A =   EL%W%A(I)* sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A

          B =    EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))* &
                sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I)) +B


       endif
ENDDO
 

 
    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 
  END SUBROUTINE COMPY_R

  SUBROUTINE COMPY_P(EL,Z,X,A,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp

    DO I=1,el%w%n   !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i))  &
               *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + B
!          A =  -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(2,i) + A
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          A =  -EL%W%A(I)*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2 *EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  -EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
!          B =  -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + A
          B =  -EL%W%A(I)*EL%W%K(3,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) & 
               *sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_dp)**2*EL%W%K(2,i)*(X(3)+EL%W%Y0(I))**2*0.5_dp &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          A =   -EL%W%A(I)*EL%W%K(1,i)*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(3,i) + A
          B =    EL%W%A(I)*EL%W%K(1,i)**2*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))) &
                *SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          A =   -EL%W%A(I)*EL%W%K(1,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i)**2 + A
          B =   -0.5_DP*EL%W%A(I)*EL%W%K(1,i)**2*EL%W%K(2,i)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))* &
                Sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I))*0.5_DP )**2 &
                *sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I))**2/EL%W%K(3,i)**2 + B
       elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          A =   -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =   -EL%W%A(I)*EL%W%K(1,i)*EL%W%K(2,i)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*Sinx_x(EL%W%K(2,i) &
                 *(X(3)+EL%W%Y0(I))*0.5_DP )**2 &
                *sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I))**2*0.5_DP/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          A =   EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i)/EL%W%K(3,i))*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) &
                *cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + A

          B =   EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i))**2*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          A =   -EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i)**2)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I)) + A
    B =   -EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))**2*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I)) + B
       elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          A =   EL%W%A(I)* sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                sin(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A

          B =    EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))* &
                sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* sin(EL%W%K(3,i)*Z+EL%W%F(I))*(X(3)+EL%W%Y0(I)) +B
       endif
ENDDO


    A=A*EL%P%CHARGE 
    B=B*EL%P%CHARGE 

  END SUBROUTINE COMPY_P

  !   Z_PLANE

  SUBROUTINE COMPZ_R(EL,Z,X,A,B)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: A,B
    INTEGER I
    A=0.0_dp
    B=0.0_dp

    DO I=1,el%w%n   !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_x) THEN
          A =  -EL%W%A(I)*EL%W%K(1,i)**2*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =  EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))   &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  EL%W%A(I)*EL%W%K(1,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(3,i) + A
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) * &
                sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
       ELSEif (EL%W%FORM(I) == hyper_y_family_y) THEN
          A =  -EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          B =  -EL%W%A(I)*EL%W%K(2,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)**2*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
                *SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
 !         B =  -EL%W%A(I)*EL%W%K(2,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
 !               COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =   EL%W%A(I)*EL%W%K(2,i)**2*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
       endif
    ENDDO


    a=a-EL%W%offset
    A=A*EL%P%CHARGE*EL%P%DIR 
    B=B*EL%P%CHARGE*EL%P%DIR 
  END SUBROUTINE COMPZ_R


  SUBROUTINE eval_thin_q(EL,Q,n)
    IMPLICIT NONE
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(OUT):: q
    INTEGER I,n

    q=0.0_dp
    n=el%n_min
    DO I=1,el%w%n   !SIZE(EL%W%A)

        q = q + (EL%W%A(I)**2/EL%W%K(3,i)) * (EL%W%K(1,i)**2+EL%W%K(2,i)**2)/2.d0

    ENDDO

  END SUBROUTINE eval_thin_q

  SUBROUTINE COMPZ_P(EL,Z,X,A,B)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: A,B
    INTEGER I
 
    A=0.0_dp
    B=0.0_dp
    DO I=1,el%w%n   !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_x) THEN
          A =  -EL%W%A(I)*EL%W%K(1,i)**2*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)**2 + A
          B =  EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))   &
               *sinhx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  EL%W%A(I)*EL%W%K(1,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i)/EL%W%K(3,i) + A
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          A =  EL%W%A(I)*EL%W%K(1,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i))) * &
                sinx_x(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))*(X(3)+EL%W%Y0(I)) * &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B
!          A =  EL%W%A(I)*EL%W%K(1,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + A
       ELSEif (EL%W%FORM(I) == hyper_y_family_y) THEN
          A =  -EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)*sinx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
               *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))  + B
!          B =  -EL%W%A(I)*EL%W%K(2,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
!                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + A
          B =  -EL%W%A(I)*EL%W%K(2,i)**2*sinhx_x(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*(X(1)+EL%W%X0(i)) &
                *SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B
 !         B =  -EL%W%A(I)*EL%W%K(2,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
 !               COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)/EL%W%K(3,i) + B
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          A =  -EL%W%A(I)*EL%W%K(2,i)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + A
          B =   EL%W%A(I)*EL%W%K(2,i)**2*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i)**2 + B
       endif
    ENDDO
    a=a-EL%W%offset
    A=A*EL%P%CHARGE*EL%P%DIR 
    B=B*EL%P%CHARGE*EL%P%DIR 
 
  END SUBROUTINE COMPZ_P

  SUBROUTINE e_potentialr(EL,Z,X,v)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: v
    INTEGER I

! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X
v=0.0_dp
 
    DO I=1,el%w%ne   !SIZE(EL%W%AE)
       if (EL%W%FORME(I) == hyper_y_family_x) THEN
          V = EL%W%AE(I)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + V
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          V = EL%W%AE(I)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + V
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          V = EL%W%AE(I)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V
       ELSEif (EL%W%FORM(I) == hyper_y_family_y) THEN
          V = EL%W%AE(I)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + V
  
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          V = EL%W%AE(I)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + V
 
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          V = EL%W%AE(I)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V

        elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(2,i))*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

        elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(3,i))*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

        elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          V = EL%W%AE(I)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V
        elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(2,i))*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

    elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(3,i))*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

        elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          V = -EL%W%AE(I)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V
       else
          print *, 'ERROR IN e_fieldp: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
    ENDDO

    V=-V*volt_c/EL%P%P0C
 
  END SUBROUTINE e_potentialr


  SUBROUTINE e_potentialp(EL,Z,X,v)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    TYPE(SAGANp),INTENT(IN):: EL
    type(real_8),INTENT(IN):: Z
    type(real_8),INTENT(INOUT):: v
    INTEGER I

! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X
v=0.0_dp
 
    DO I=1,el%w%ne  !SIZE(EL%W%AE)
       if (EL%W%FORME(I) == hyper_y_family_x) THEN
          V = EL%W%AE(I)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + V
       elseif (EL%W%forme(I) == hyper_xy_family_x) THEN
          V = EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + V
       elseif (EL%W%forme(I) == hyper_x_family_x) THEN
          V = EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V
       ELSEif (EL%W%forme(I) == hyper_y_family_y) THEN
          V = EL%W%AE(I)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + V
  
       elseif (EL%W%forme(I) == hyper_xy_family_y) THEN
          V = EL%W%AE(I)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + V
 
       elseif (EL%W%forme(I) == hyper_x_family_y) THEN
          V = EL%W%AE(I)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V

        elseif (EL%W%forme(I) == hyper_y_family_qu) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(2,i))*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

        elseif (EL%W%forme(I) == hyper_xy_family_qu) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(3,i))*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

        elseif (EL%W%forme(I) == hyper_x_family_qu) THEN
          V = EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V
        elseif (EL%W%forme(I) == hyper_y_family_sq) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(2,i))*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

    elseif (EL%W%forme(I) == hyper_xy_family_sq) THEN
          V = EL%W%AE(I)*(1.0_DP/EL%W%KE(3,i))*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + V

        elseif (EL%W%forme(I) == hyper_x_family_sq) THEN
          V = -EL%W%AE(I)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + V
       else
          print *, 'ERROR IN e_fieldr: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
    ENDDO

    V=-V*volt_c/EL%P%P0C
 
  END SUBROUTINE e_potentialp


  SUBROUTINE e_fieldr(EL,Z,X,e,kick)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: E(3)
    logical, optional :: kick
    INTEGER I,dire(3)
    e=0.0_dp

! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X

 
    DO I=1,el%w%ne  !SIZE(EL%W%AE)

  if (EL%W%FORME(I) == hyper_y_family_x) THEN


     E(1) = EL%W%AE(I)*EL%W%KE(1,i)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(1)
     E(2) = EL%W%AE(I)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(3)
  elseif (EL%W%forme(I) == hyper_xy_family_x) THEN
     E(1) = EL%W%AE(I)*EL%W%KE(1,i)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(1)
     E(2) = EL%W%AE(I)*EL%W%KE(2,i)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(2) 
     E(3) = -EL%W%AE(I)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
  elseif (EL%W%forme(I) == hyper_x_family_x) THEN
     E(1) = EL%W%AE(I)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = -EL%W%AE(I)*EL%W%KE(2,i)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(2) 
     E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(3)

  ELSEif (EL%W%forme(I) == hyper_y_family_y) THEN


     E(1) = -EL%W%AE(I)*EL%W%KE(1,i)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(1)
     E(2) = EL%W%AE(I)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(3)
  elseif (EL%W%forme(I) == hyper_xy_family_y) THEN
     E(1) = EL%W%AE(I)*EL%W%KE(1,i)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(1)
     E(2) = EL%W%AE(I)*EL%W%KE(2,i)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COSEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(2) 
     E(3) = -EL%W%AE(I)*COSH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SINEH(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
  elseif (EL%W%forme(I) == hyper_x_family_y) THEN
     E(1) = EL%W%AE(I)*SINEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*EL%W%KE(2,i)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(2) 
     E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*COSEH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(3)
   elseif (EL%W%forme(I) == hyper_y_family_qu) THEN
     E(1) = EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(2,i))*cos(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*sin(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(2,i))*sin(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
   elseif (EL%W%forme(I) == hyper_xy_family_qu) THEN
     E(1) = EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(3,i))*coseh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(3,i))*sineh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*sineh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
   elseif (EL%W%forme(I) == hyper_x_family_qu) THEN
     E(1) = EL%W%AE(I)*coseh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(1,i))*sineh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(1,i))*sineh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
   elseif (EL%W%forme(I) == hyper_y_family_sq) THEN
     E(1) = -EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(2,i))*sin(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*cos(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(2,i))*cos(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
    elseif (EL%W%forme(I) == hyper_xy_family_sq) THEN
     E(1) = EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(3,i))*sineh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(3,i))*coseh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sineh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = -EL%W%AE(I)*coseh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*coseh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
   elseif (EL%W%forme(I) == hyper_x_family_sq) THEN
     E(1) = -EL%W%AE(I)*sineh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
     E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(1,i))*coseh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
     E(3) = EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(1,i))*coseh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
           SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
  else
     print *, 'ERROR IN e_fieldr: UNKNOWN FORM FOR WIGGLER TERM!'
     stop
  endif

    ENDDO

      do i=1,3
        e(i)=e(i)*volt_c/EL%P%P0C
      enddo

      if(present(kick)) then 
     if(kick) then
      DIRE=EL%P%DIR; DIRE(1:2)=1;
      do i=1,3
        e(i)=DIRE(I)*el%p%charge*e(i) 
      enddo
     endif
     endif
 
  END SUBROUTINE e_fieldr

 SUBROUTINE e_fieldP(EL,Z,X,e,kick)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    TYPE(SAGANp),INTENT(IN):: EL
    type(real_8),INTENT(IN):: Z
    type(real_8),INTENT(INOUT):: E(3)
    logical, optional :: kick
    INTEGER I,dire(3)

    E(1)=0.0_dp;E(2)=0.0_dp;E(3)=0.0_dp;
! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X

 
    DO I=1,el%w%ne   !SIZE(EL%W%AE)
       if (EL%W%FORME(I) == hyper_y_family_x) THEN


          E(1) = EL%W%AE(I)*EL%W%KE(1,i)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(1)
          E(2) = EL%W%AE(I)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                 COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
          E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(3)
       elseif (EL%W%forme(I) == hyper_xy_family_x) THEN
          E(1) = EL%W%AE(I)*EL%W%KE(1,i)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(1)
          E(2) = EL%W%AE(I)*EL%W%KE(2,i)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                 COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(2) 
          E(3) = -EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
       elseif (EL%W%forme(I) == hyper_x_family_x) THEN
          E(1) = EL%W%AE(I)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
          E(2) = -EL%W%AE(I)*EL%W%KE(2,i)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                 COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(2) 
          E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(3)
 
       ELSEif (EL%W%forme(I) == hyper_y_family_y) THEN


          E(1) = -EL%W%AE(I)*EL%W%KE(1,i)*SIN(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(1)
          E(2) = EL%W%AE(I)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                 COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
          E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*COS(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(2,i) + E(3)
       elseif (EL%W%forme(I) == hyper_xy_family_y) THEN
          E(1) = EL%W%AE(I)*EL%W%KE(1,i)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(1)
          E(2) = EL%W%AE(I)*EL%W%KE(2,i)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                 COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(3,i) + E(2) 
          E(3) = -EL%W%AE(I)*COSH(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
       elseif (EL%W%forme(I) == hyper_x_family_y) THEN
          E(1) = EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
          E(2) = EL%W%AE(I)*EL%W%KE(2,i)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*COS(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                 COS(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(2) 
          E(3) = -EL%W%AE(I)*EL%W%KE(3,i)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*SIN(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
                SIN(EL%W%KE(3,i)*Z+EL%W%FE(I))/EL%W%KE(1,i) + E(3)
        elseif (EL%W%forme(I) == hyper_y_family_qu) THEN
       E(1) = EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(2,i))*cos(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
       E(2) = EL%W%AE(I)*sin(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
              COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
       E(3) = -EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(2,i))*sin(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
     elseif (EL%W%forme(I) == hyper_xy_family_qu) THEN
       E(1) = EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(3,i))*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
       E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(3,i))*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
              COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
       E(3) = -EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
     elseif (EL%W%forme(I) == hyper_x_family_qu) THEN
       E(1) = EL%W%AE(I)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
       E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(1,i))*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
              COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
       E(3) = -EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(1,i))*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
     elseif (EL%W%forme(I) == hyper_y_family_sq) THEN
       E(1) = -EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(2,i))*sin(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
       E(2) = EL%W%AE(I)*cos(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
              COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
       E(3) = -EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(2,i))*cos(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
 elseif (EL%W%forme(I) == hyper_xy_family_sq) THEN
   E(1) = EL%W%AE(I)*(EL%W%KE(1,i)/EL%W%KE(3,i))*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
         COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
   E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(3,i))*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sinh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
          COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
   E(3) = -EL%W%AE(I)*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cosh(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
         SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
 elseif (EL%W%forme(I) == hyper_x_family_sq) THEN
   E(1) = -EL%W%AE(I)*sinh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(1)
      E(2) = EL%W%AE(I)*(EL%W%KE(2,i)/EL%W%KE(1,i))*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*sin(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
             COS(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(2) 
      E(3) = EL%W%AE(I)*(EL%W%KE(3,i)/EL%W%KE(1,i))*cosh(EL%W%KE(1,i)*(X(1)+EL%W%X0E(i)))*cos(EL%W%KE(2,i)*(X(3)+EL%W%Y0E(I)))* &
            SIN(EL%W%KE(3,i)*Z+EL%W%FE(I)) + E(3)
   else
      print *, 'ERROR IN e_fieldp: UNKNOWN FORM FOR WIGGLER TERM!'
      stop
   endif
    ENDDO

       do i=1,3
        e(i)=e(i)*volt_c/EL%P%P0C
      enddo

      if(present(kick)) then 
     if(kick) then
      DIRE=EL%P%DIR; DIRE(1:2)=1;
      do i=1,3
        e(i)=DIRE(I)*el%p%charge*e(i) 
      enddo
     endif
     endif

  END SUBROUTINE e_fieldP

SUBROUTINE K_normal(p)
    IMPLICIT NONE
    TYPE(fibre),target:: p
   TYPE(undu_R),pointer:: w
   TYPE(undu_p),pointer:: wp
    real(dp) k
    INTEGER I

w=>p%mag%wi%w
wp=>p%magp%wi%w

! result = SIGN (a, b)   sign(b)*|a|

    DO I=1,w%n  !SIZE(w%A)
       if (w%FORM(I) == hyper_y_family_x) THEN
          
          k=sqrt(w%k(1,i)**2+w%k(3,i)**2)
          w%k(2,i)=sign(k,w%k(2,i))
          wp%k(2,i)=w%k(2,i) 

       elseif (w%FORM(I) == hyper_xy_family_x) THEN

          k=sqrt(w%k(1,i)**2+w%k(2,i)**2)
          w%k(3,i)=sign(k,w%k(3,i))
          wp%k(3,i)=w%k(3,i)

       elseif (w%FORM(I) == hyper_x_family_x) THEN

          k=sqrt(w%k(2,i)**2+w%k(3,i)**2)
          w%k(1,i)=sign(k,w%k(1,i))
          wp%k(1,i)=w%k(1,i)
 
       ELSEif (w%FORM(I) == hyper_y_family_y) THEN

          k=sqrt(w%k(1,i)**2+w%k(3,i)**2)
          w%k(2,i)=sign(k,w%k(2,i))
          wp%k(2,i)=w%k(2,i)

       elseif (w%FORM(I) == hyper_xy_family_y) THEN


          k=sqrt(w%k(1,i)**2+w%k(2,i)**2)
          w%k(3,i)=sign(k,w%k(3,i))
          wp%k(3,i)=w%k(3,i)   

       elseif (w%FORM(I) == hyper_x_family_y) THEN

          k=sqrt(w%k(2,i)**2+w%k(3,i)**2)
          w%k(1,i)=sign(k,w%k(1,i))
          wp%k(1,i)=w%k(1,i)

        elseif (w%FORM(I) == hyper_y_family_qu) THEN

          k=sqrt(w%k(1,i)**2+w%k(3,i)**2)
          w%k(2,i)=sign(k,w%k(2,i))
          wp%k(2,i)=w%k(2,i)

        elseif (w%FORM(I) == hyper_xy_family_qu) THEN


          k=sqrt(w%k(1,i)**2+w%k(2,i)**2)
          w%k(3,i)=sign(k,w%k(3,i))
          wp%k(3,i)=w%k(3,i)

        elseif (w%FORM(I) == hyper_x_family_qu) THEN

          k=sqrt(w%k(2,i)**2+w%k(3,i)**2)
          w%k(1,i)=sign(k,w%k(1,i))
          wp%k(1,i)=w%k(1,i)
     
        elseif (w%FORM(I) == hyper_y_family_sq) THEN

          k=sqrt(w%k(1,i)**2+w%k(3,i)**2)
          w%k(2,i)=sign(k,w%k(2,i))
          wp%k(2,i)=w%k(2,i)

    elseif (w%FORM(I) == hyper_xy_family_sq) THEN


          k=sqrt(w%k(1,i)**2+w%k(2,i)**2)
          w%k(3,i)=sign(k,w%k(3,i))
          wp%k(3,i)=w%k(3,i)

        elseif (w%FORM(I) == hyper_x_family_sq) THEN

          k=sqrt(w%k(2,i)**2+w%k(3,i)**2)
          w%k(1,i)=sign(k,w%k(1,i))
          wp%k(1,i)=w%k(1,i)

       else
          print *, 'ERROR IN k_normal: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
    ENDDO

    DO I=1,w%ne   !SIZE(w%AE)
       if (w%FORME(I) == hyper_y_family_x) THEN
          
          k=sqrt(w%KE(1,i)**2+w%KE(3,i)**2)
          w%KE(2,i)=sign(k,w%KE(2,i))
          wp%KE(2,i)=w%KE(2,i) 

       elseif (w%FORME(I) == hyper_xy_family_x) THEN

          k=sqrt(w%KE(1,i)**2+w%KE(2,i)**2)
          w%KE(3,i)=sign(k,w%KE(3,i))
          wp%KE(3,i)=w%KE(3,i)

       elseif (w%FORME(I) == hyper_x_family_x) THEN

          k=sqrt(w%KE(2,i)**2+w%KE(3,i)**2)
          w%KE(1,i)=sign(k,w%KE(1,i))
          wp%KE(1,i)=w%KE(1,i)
 
       ELSEif (w%FORME(I) == hyper_y_family_y) THEN

          k=sqrt(w%KE(1,i)**2+w%KE(3,i)**2)
          w%KE(2,i)=sign(k,w%KE(2,i))
          wp%KE(2,i)=w%KE(2,i)

       elseif (w%FORME(I) == hyper_xy_family_y) THEN


          k=sqrt(w%KE(1,i)**2+w%KE(2,i)**2)
          w%KE(3,i)=sign(k,w%KE(3,i))
          wp%KE(3,i)=w%KE(3,i)   

       elseif (w%FORME(I) == hyper_x_family_y) THEN

          k=sqrt(w%KE(2,i)**2+w%KE(3,i)**2)
          w%KE(1,i)=sign(k,w%KE(1,i))
          wp%KE(1,i)=w%KE(1,i)

        elseif (w%FORME(I) == hyper_y_family_qu) THEN

          k=sqrt(w%KE(1,i)**2+w%KE(3,i)**2)
          w%KE(2,i)=sign(k,w%KE(2,i))
          wp%KE(2,i)=w%KE(2,i)

        elseif (w%FORME(I) == hyper_xy_family_qu) THEN


          k=sqrt(w%KE(1,i)**2+w%KE(2,i)**2)
          w%KE(3,i)=sign(k,w%KE(3,i))
          wp%KE(3,i)=w%KE(3,i)

        elseif (w%FORME(I) == hyper_x_family_qu) THEN

          k=sqrt(w%KE(2,i)**2+w%KE(3,i)**2)
          w%KE(1,i)=sign(k,w%KE(1,i))
          wp%KE(1,i)=w%KE(1,i)
     
        elseif (w%FORME(I) == hyper_y_family_sq) THEN

          k=sqrt(w%KE(1,i)**2+w%KE(3,i)**2)
          w%KE(2,i)=sign(k,w%KE(2,i))
          wp%KE(2,i)=w%KE(2,i)

    elseif (w%FORME(I) == hyper_xy_family_sq) THEN


          k=sqrt(w%KE(1,i)**2+w%KE(2,i)**2)
          w%KE(3,i)=sign(k,w%KE(3,i))
          wp%KE(3,i)=w%KE(3,i)

        elseif (w%FORME(I) == hyper_x_family_sq) THEN

          k=sqrt(w%KE(2,i)**2+w%KE(3,i)**2)
          w%KE(1,i)=sign(k,w%KE(1,i))
          wp%KE(1,i)=w%KE(1,i)

       else
          print *, 'ERROR IN k_normal: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
    ENDDO



 
  END SUBROUTINE K_normal

  SUBROUTINE b_fieldr(EL,Z,X,B,kick)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: Z
    real(dp),INTENT(INOUT):: B(3)
    logical, optional :: kick
    INTEGER I,dir(3)
    B=0.0_dp

! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X

 
    DO I=1,el%w%n  !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_x) THEN


          B(1) = EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          B(1) = EL%W%A(I)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = -EL%W%A(I)*EL%W%K(2,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
 
       ELSEif (EL%W%FORM(I) == hyper_y_family_y) THEN


          B(1) = -EL%W%A(I)*EL%W%K(1,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINEH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          B(1) = EL%W%A(I)*SINEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COSEH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
        elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          B(1) = EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i))*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(2,i))*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          B(1) = EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(3,i))*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          B(1) = EL%W%A(I)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(1,i))*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(1,i))*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          B(1) = -EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i))*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(2,i))*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
    elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          B(1) = EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(3,i))*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sineh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*coseh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          B(1) = -EL%W%A(I)*sineh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(1,i))*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(1,i))*coseh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       else
          print *, 'ERROR IN b_fieldr: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
    ENDDO


       b(2)=b(2)+el%w%offset

     if(present(kick)) then 
     if(kick) then
      DIR=EL%P%DIR; DIR(3)=1;
      do i=1,3
        b(i)=DIR(I)*el%p%charge*b(i)
      enddo
     endif
     endif

 
  END SUBROUTINE b_fieldr

  SUBROUTINE b_fieldp(EL,Z,X,B,kick)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: Z
    TYPE(REAL_8),INTENT(INOUT):: B(3)  
    logical, optional :: kick
    INTEGER I,dir(3)
    B(1)=0.0_dp;B(2)=0.0_dp;B(3)=0.0_dp;

 
! COSEH(X) ! REPLACES COSH(X)
! SINEH(X) ! REPLACES SINH(X)
! SINEHX_X(X) ! REPLACES SINH(X)/X
 
 
    DO I=1,el%w%n  !SIZE(EL%W%A)
       if (EL%W%FORM(I) == hyper_y_family_x) THEN


          B(1) = EL%W%A(I)*EL%W%K(1,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_family_x) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_family_x) THEN
          B(1) = EL%W%A(I)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = -EL%W%A(I)*EL%W%K(2,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
 
 
       ELSEif (EL%W%FORM(I) == hyper_y_family_y) THEN


           B(1) = -EL%W%A(I)*EL%W%K(1,i)*SIN(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(1)
          B(2) = EL%W%A(I)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COS(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(2,i) + B(3)
       elseif (EL%W%FORM(I) == hyper_xy_family_y) THEN
          B(1) = EL%W%A(I)*EL%W%K(1,i)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COSH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(3,i) + B(2) 
          B(3) = -EL%W%A(I)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SINH(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       elseif (EL%W%FORM(I) == hyper_x_family_y) THEN
          B(1) = EL%W%A(I)*SINH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*EL%W%K(2,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*COS(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(2) 
          B(3) = -EL%W%A(I)*EL%W%K(3,i)*COSH(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*SIN(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I))/EL%W%K(1,i) + B(3)
        elseif (EL%W%FORM(I) == hyper_y_family_qu) THEN
          B(1) = EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i))*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(2,i))*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_xy_family_qu) THEN
          B(1) = EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(3,i))*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_x_family_qu) THEN
          B(1) = EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(1,i))*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(1,i))*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_y_family_sq) THEN
          B(1) = -EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(2,i))*sin(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(2,i))*cos(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
    elseif (EL%W%FORM(I) == hyper_xy_family_sq) THEN
          B(1) = EL%W%A(I)*(EL%W%K(1,i)/EL%W%K(3,i))*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(3,i))*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sinh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = -EL%W%A(I)*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cosh(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
        elseif (EL%W%FORM(I) == hyper_x_family_sq) THEN
          B(1) = -EL%W%A(I)*sinh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(1)
          B(2) = EL%W%A(I)*(EL%W%K(2,i)/EL%W%K(1,i))*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*sin(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                 COS(EL%W%K(3,i)*Z+EL%W%F(I)) + B(2) 
          B(3) = EL%W%A(I)*(EL%W%K(3,i)/EL%W%K(1,i))*cosh(EL%W%K(1,i)*(X(1)+EL%W%X0(i)))*cos(EL%W%K(2,i)*(X(3)+EL%W%Y0(I)))* &
                SIN(EL%W%K(3,i)*Z+EL%W%F(I)) + B(3)
       else
          print *, 'ERROR IN b_fieldp: UNKNOWN FORM FOR WIGGLER TERM!'
          stop
       endif
     enddo

 
    b(2)=b(2)+el%w%offset
 
     if(present(kick)) then 
     if(kick) then
      DIR=EL%P%DIR; DIR(3)=1;
      do i=1,3
        b(i)=DIR(I)*el%p%charge*b(i)
      enddo
     endif
     endif


  END SUBROUTINE b_fieldp

  SUBROUTINE KICKR(EL,L,Z,X,k)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(SAGAN),INTENT(IN):: EL
    real(dp),INTENT(IN):: L,Z
    real(dp) A,B,kx,ky
    TYPE(INTERNAL_STATE),OPTIONAL :: K


    CALL COMPZ(EL,Z,X,A,B)
    X(2)=X(2)+L*A
    X(4)=X(4)+L*B

    if(el%p%permfringe>0) then
     call kick_integral(el,x,kx,ky,el%p%permfringe)
     X(2)=X(2)+L*kx
     X(4)=X(4)+L*ky
    endif

  END SUBROUTINE KICKR

  SUBROUTINE KICKP(EL,L,Z,X,k)
    IMPLICIT NONE
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    TYPE(SAGANP),INTENT(IN):: EL
    TYPE(REAL_8),INTENT(IN):: L,Z
    TYPE(REAL_8) A,B,AP,BP,kx,ky
    TYPE(INTERNAL_STATE),OPTIONAL :: K


    call alloc(A,B,AP,BP)


    CALL COMPZ(EL,Z,X,A,B)
    X(2)=X(2)+L*A
    X(4)=X(4)+L*B

    if(el%p%permfringe>0) then
     call alloc(kx,ky)
     call kick_integral(el,x,kx,ky,el%p%permfringe)
     X(2)=X(2)+L*kx
     X(4)=X(4)+L*ky
     call kill(kx,ky)
    endif

    call KILL(A,B,AP,BP)

  END SUBROUTINE KICKP

  SUBROUTINE PRINT_(EL,MF)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(INOUT)::EL
    INTEGER MF,I

    !    WRITE(MF,*) EL%INTERNAL
    WRITE(MF,100) "NUMBER OF TERMS ",SIZE(EL%W%A)
    DO I=1,el%w%n  !SIZE(EL%W%A)
       WRITE(MF,101) " A = ",EL%W%A(I), " K = ",EL%W%K(1,I),EL%W%K(2,I),EL%W%K(3,I), &
            & " PHASE = ",EL%W%F(I)," FORM = ",EL%W%FORM(I)
    ENDDO
    DO I=1,el%w%ne  !SIZE(EL%W%A)
       WRITE(MF,101) " A = ",EL%W%Ae(I), " K = ",EL%W%Ke(1,I),EL%W%K(2,I),EL%W%Ke(3,I), &
            & " PHASE = ",EL%W%Fe(I)," FORM = ",EL%W%FORMe(I)
    ENDDO
    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS
100 FORMAT(A16,(1X,I4))
101 FORMAT(A5,(1X,g21.14),A5,3(1X,g21.14),A9,(1X,g21.14),A11,I3)

  END SUBROUTINE PRINT_

  SUBROUTINE READ_(EL,MF)
    IMPLICIT NONE
    TYPE(SAGAN), INTENT(INOUT)::EL
    INTEGER MF,I

    READ(MF,100) I

    CALL POINTERS_W(EL%W,I,i)
    DO I=1,SIZE(EL%W%A)
       READ(MF,101) EL%W%A(I), EL%W%K(1,I),EL%W%K(2,I),EL%W%K(3,I),EL%W%F(I),EL%W%FORM(I)
    ENDDO

    !   READ(MF,*) EL%INTERNAL
    ! CALL resetpoly_R31 ON ALL THE INTERNAL POLYMORPHS
100 FORMAT(16X,(1X,I4))
101 FORMAT(5X,(1X,g21.14),5X,3(1X,g21.14),9X,(1X,g21.14),11X,I3)

  END SUBROUTINE READ_


 subroutine feval_saganr(Z,X,k,f,EL)   !electric teapot s
    IMPLICIT NONE
    real(dp), INTENT(INout) :: X(6)
    real(dp), INTENT(INOUT) :: F(6)
    real(dp), INTENT(INOUT) :: Z
    REAL(DP) PZ,DEL,H,B(3),E(3),VE,beta0   ! vm MAGNETIC POTENTIAL
    TYPE(sagan),  INTENT(INout) :: EL
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
     
    if(k%time) then
       beta0=el%p%beta0; 
    else
       beta0=1.0_dp; 
    endif


      e=0.0_dp
      ve=0.0_dp

      call B_field(EL,Z,X,B)
      call E_field(EL,Z,X,E)
      call E_potential(EL,Z,X,VE)

           DEL=x(5)+VE
      call fx_new(f,x,k,EL%P%EXACT,EL%P%b0,beta0,b,e,ve) 

     global_e= DEL*el%p%p0c
   END subroutine feval_saganr


subroutine feval_saganp(Z,X,k,f,EL)   !electric teapot s
    IMPLICIT NONE
    TYPE(REAL_8), INTENT(INout) :: X(6)
    TYPE(REAL_8),  INTENT(INOUT) :: F(6)
    TYPE(REAL_8),   INTENT(INOUT) :: Z
    TYPE(REAL_8) PZ,DEL,H,B(3),E(3),VM,VE   ! vm MAGNETIC POTENTIAL
    TYPE(saganp),  INTENT(INout) :: EL
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    integer i
    real(dp) beta0

     call alloc(PZ,DEL,H,VM)
     call alloc(E);call alloc(B);
      CALL alloc(VE)

    if(k%time) then
       beta0=el%p%beta0; 
    else
       beta0=1.0_dp; 
    endif

      do i=1,3
       e(i)=0.0_dp
      enddo
      ve=0.0_dp

      call B_field(EL,Z,X,B)
      call E_field(EL,Z,X,E)
      call E_potential(EL,Z,X,VE)

           DEL=x(5)+VE
      call fx_new(f,x,k,EL%P%EXACT,EL%P%b0,beta0,b,e,ve) 

     global_e= DEL*el%p%p0c

     call kill(PZ,DEL,H,VM)
     call kill(E);call kill(B);

      CALL KILL(VE)
   END subroutine feval_saganp


 subroutine rk2saganr(ti,h,GR,y,k)
    IMPLICIT none

    integer ne
    parameter (ne=6)
    real(dp), INTENT(INOUT)::  y(ne)
    real(dp)  yt(ne),f(ne),a(ne),b(ne)
    real(dp)  tt
    type (sagan) ,INTENT(INOUT)::  GR
    integer j
    real(dp), intent(inout) :: ti,h
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K


    call feval_sagan(tI,y,k,f,gr)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/2.0_dp
    enddo

    tt=tI+h/2.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+b(j)
    enddo
    tI=ti+h

    return
  end  subroutine rk2saganr

 subroutine rk4saganr(ti,h,GR,y,k)
    IMPLICIT none

    integer ne
    parameter (ne=6)
    real(dp), INTENT(INOUT)::  y(ne)
    real(dp)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne)
    type (sagan) ,INTENT(INOUT)::  GR
    integer j
    real(dp), intent(inout) :: h
    real(dp), intent(inout) :: ti
    real(dp) TT
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K


    call feval_sagan(tI,y,k,f,gr)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/2.0_dp
    enddo

    tt=tI+h/2.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + b(j)/2.0_dp
    enddo


    !      tt=tI+1
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       c(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+c(j)
    enddo

    tt=tI+h
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       d(j)=h*f(j)
    enddo


    do  j=1,ne
       y(j) = y(j)+(a(j)+2.0_dp*b(j)+2.0_dp*c(j)+d(j))/6.0_dp
    enddo
    tI=tt

    return
  end  subroutine rk4saganr

 subroutine rk6saganr(ti,h,GR,y,k)
    IMPLICIT none
    !  Written by Rob Ryne, Spring 1986, based on a routine of
    !c  J. Milutinovic.
    !c  For a reference, see page 76 of F. Ceschino and J Kuntzmann,
    !c  Numerical Solution of Initial Value Problems, Prentice Hall 1966.
    !c  This integration routine makes local truncation errors at each
    !c  step of order h**7.
    !c  That is, it is locally correct through terms of order h**6.
    !c  Each step requires 8 function evaluations.

    integer ne
    parameter (ne=6)
    real(dp), INTENT(INOUT)::  y(ne)
    real(dp)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne),e(ne),g(ne),o(ne),p(ne)
    real(dp)  tt
    type (sagan) ,INTENT(INOUT)::  GR
    integer j
    real(dp), intent(inout) :: ti,h
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K


    call feval_sagan(tI,y,k,f,gr)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/9.0_dp
    enddo
    tt=tI+h/9.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + (a(j) + 3.0_dp*b(j))/24.0_dp
    enddo
    tt=tI+h/6.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       c(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j)+(a(j)-3.0_dp*b(j)+4.0_dp*c(j))/6.0_dp
    enddo

    tt=tI+h/3.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       d(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j) + (-5.0_dp*a(j) + 27.0_dp*b(j) - 24.0_dp*c(j) + 6.0_dp*d(j))/8.0_dp
    enddo
    tt=tI+0.5_dp*h
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       e(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j) + (221.0_dp*a(j) - 981.0_dp*b(j) + 867.0_dp*c(j)- 102.0_dp*d(j) + e(j))/9.0_dp
    enddo
    tt = tI+2.0_dp*h/3.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do   j=1,ne
       g(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j) = y(j)+(-183.0_dp*a(j)+678.0_dp*b(j)-472.0_dp*c(j)-66.0_dp*d(j)+80.0_dp*e(j) + 3.0_dp*g(j))/48.0_dp
    enddo
    tt = tI + 5.0_dp*h/6.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       o(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j) = y(j)+(716.0_dp*a(j)-2079.0_dp*b(j)+1002.0_dp*c(j)+834.0_dp*d(j)-454.0_dp*e(j)-9.0_dp*g(j)+72.0_dp*o(j))/82.0_dp
    enddo

    tt = tI + h
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       p(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+(41.0_dp*a(j)+216.0_dp*c(j)+27.0_dp*d(j)+272.0_dp*e(j)+27.0_dp*g(j)+216.0_dp*o(j)+41.0_dp*p(j))/840.0_dp
    enddo
    tI=ti+h

    return
  end  subroutine rk6saganr

  subroutine rk2saganp(ti,h,GR,y,k)
    IMPLICIT none

    integer ne
    parameter (ne=6)
    type (real_8), INTENT(INOUT)::  y(ne)
    type (real_8)  yt(ne),f(ne),a(ne),b(ne)
    type (real_8)  tt
    type (saganp) ,INTENT(INOUT)::  GR
    integer j
    type(real_8), intent(inout) :: ti,h
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    call alloc(yt,ne)
    call alloc(f,ne)
    call alloc(a,ne)
    call alloc(b,ne)

    call alloc(tt)

    call feval_sagan(tI,y,k,f,gr)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/2.0_dp
    enddo

    tt=tI+h/2.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+b(j)
    enddo
    tI=ti+h

    call kill(tt)
    call kill(yt,ne)
    call kill(f,ne)
    call kill(a,ne)
    call kill(b,ne)

    return
  end  subroutine rk2saganp

  subroutine rk4saganp(ti,h,GR,y,k)
    IMPLICIT none

    integer ne
    parameter (ne=6)
    type(real_8), INTENT(INOUT)::  y(ne)
    type (saganp) ,INTENT(INOUT)::  GR
    type(real_8), intent(inout) :: h
    type(real_8), intent(inout) :: ti
    type(real_8)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne)
    type(real_8) TT
    integer j
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    call alloc(tt)
    call alloc(yt)
    call alloc(f)
    call alloc(a)
    call alloc(b)
    call alloc(c)
    call alloc(d)

    call feval_sagan(tI,y,k,f,gr)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/2.0_dp
    enddo

    tt=tI+h/2.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + b(j)/2.0_dp
    enddo


    !      tt=tI+1
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       c(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+c(j)
    enddo

    tt=tI+h
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       d(j)=h*f(j)
    enddo


    do  j=1,ne
       y(j) = y(j)+(a(j)+2.0_dp*b(j)+2.0_dp*c(j)+d(j))/6.0_dp
    enddo
    tI=tt

    call kill(tt)
    call kill(yt)
    call kill(f)
    call kill(a)
    call kill(b)
    call kill(c)
    call kill(d)

    return
  end  subroutine rk4saganp

  ! sixth order Runge
  subroutine rk6saganp(ti,h,GR,y,k)
    IMPLICIT none


    integer ne
    parameter (ne=6)
    type (real_8), INTENT(INOUT)::  y(ne)
    type (real_8)  yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne),e(ne),g(ne),o(ne),p(ne)
    type (real_8)  tt
    type (saganp) ,INTENT(INOUT)::  GR
    integer j
    type(real_8), intent(inout) :: ti,h
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    call alloc(yt,ne)
    call alloc(f,ne)
    call alloc(a,ne)
    call alloc(b,ne)
    call alloc(c,ne)
    call alloc(d,ne)
    call alloc(e,ne)
    call alloc(g,ne)
    call alloc(o,ne)
    call alloc(p,ne)
    call alloc(tt)

    call feval_sagan(tI,y,k,f,gr)
    do  j=1,ne
       a(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j)=y(j)+a(j)/9.0_dp
    enddo
    tt=tI+h/9.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       b(j)=h*f(j)
    enddo
    do   j=1,ne
       yt(j)=y(j) + (a(j) + 3.0_dp*b(j))/24.0_dp
    enddo
    tt=tI+h/6.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       c(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j)+(a(j)-3.0_dp*b(j)+4.0_dp*c(j))/6.0_dp
    enddo

    tt=tI+h/3.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       d(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j) + (-5.0_dp*a(j) + 27.0_dp*b(j) - 24.0_dp*c(j) + 6.0_dp*d(j))/8.0_dp
    enddo
    tt=tI+0.5_dp*h
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       e(j)=h*f(j)
    enddo

    do  j=1,ne
       yt(j)=y(j) + (221.0_dp*a(j) - 981.0_dp*b(j) + 867.0_dp*c(j)- 102.0_dp*d(j) + e(j))/9.0_dp
    enddo
    tt = tI+2.0_dp*h/3.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do   j=1,ne
       g(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j) = y(j)+(-183.0_dp*a(j)+678.0_dp*b(j)-472.0_dp*c(j)-66.0_dp*d(j)+80.0_dp*e(j) + 3.0_dp*g(j))/48.0_dp
    enddo
    tt = tI + 5.0_dp*h/6.0_dp
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       o(j)=h*f(j)
    enddo
    do  j=1,ne
       yt(j) = y(j)+(716.0_dp*a(j)-2079.0_dp*b(j)+1002.0_dp*c(j)+834.0_dp*d(j)-454.0_dp*e(j)-9.0_dp*g(j)+72.0_dp*o(j))/82.0_dp
    enddo

    tt = tI + h
    call feval_sagan(tt,yt,k,f,gr)
    do  j=1,ne
       p(j)=h*f(j)
    enddo

    do  j=1,ne
       y(j) = y(j)+(41.0_dp*a(j)+216.0_dp*c(j)+27.0_dp*d(j)+272.0_dp*e(j)+27.0_dp*g(j)+216.0_dp*o(j)+41.0_dp*p(j))/840.0_dp
    enddo
    tI=ti+h
    call kill(tt)
    call kill(yt,ne)
    call kill(f,ne)
    call kill(a,ne)
    call kill(b,ne)
    call kill(c,ne)
    call kill(d,ne)
    call kill(e,ne)
    call kill(g,ne)
    call kill(o,ne)
    call kill(p,ne)

    return
  end  subroutine rk6saganp

  SUBROUTINE ADJUST_like_abellr(EL,X,k,J)
    IMPLICIT NONE
    real(dp), INTENT(INOUT) :: X(6)
    TYPE(sagan),INTENT(INOUT):: EL
    INTEGER, INTENT(IN) :: J
    real(dp) z
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

     if(.not.el%xprime.or.get_out) return

     IF(J==1) then
        z=0
      IF(EL%P%DIR==1) THEN
        call conv_to_xp(el,x,k,z)
      ELSE
        call conv_to_px(el,x,k,z)
      ENDIF
    else
        z=el%l
      IF(EL%P%DIR==1) THEN
          call conv_to_px(el,x,k,z)
      ELSE
          call conv_to_xp(el,x,k,z)
      ENDIF
    endif

  END SUBROUTINE ADJUST_like_abellr

  SUBROUTINE ADJUST_like_abellp(EL,X,k,J)
    IMPLICIT NONE
    type(real_8), INTENT(INOUT) :: X(6)
    TYPE(saganp),INTENT(INOUT):: EL
    INTEGER, INTENT(IN) :: J
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    type(real_8)   z
    call ALLOC(z)
     if(.not.el%xprime.or.get_out) return

     IF(J==1) then
        z=0.0_dp
      IF(EL%P%DIR==1) THEN
        call conv_to_xp(el,x,k,z)
      ELSE
        call conv_to_px(el,x,k,z)
      ENDIF
    else
        z=el%l
      IF(EL%P%DIR==1) THEN
          call conv_to_px(el,x,k,z)
      ELSE
          call conv_to_xp(el,x,k,z)
      ENDIF
    endif
    call kill(z)
  END SUBROUTINE ADJUST_like_abellp


  SUBROUTINE conv_to_xprsagan(EL,X,k,z)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(sagan),INTENT(INOUT):: EL
    real(dp) ti,ve,z,a(3),beta0
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K


 !   call B_E_FIELD(EL,X,Z,PSIE_IN=VE,A_in=a,kick=.true.)
       CALL  COMPX(EL,Z,X,A(1),a(3))
       CALL  COMPy(EL,Z,X,A(2),a(3))
     call E_potential(EL,Z,X,VE)
      if(k%TIME) then
       beta0=el%p%beta0
      else
       beta0=1.0_dp
      endif
      call gen_conv_to_xp(X,a,ve,el%p%exact,beta0,0.0_dp)

  end SUBROUTINE conv_to_xprsagan

  SUBROUTINE conv_to_xppsagan(EL,X,k,z)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    TYPE(saganp),INTENT(INOUT):: EL
    type(real_8) ti,ve,z,a(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    real(dp) beta0

    call alloc(ti,ve)
    call alloc(a)
 
 !   call B_E_FIELD(EL,X,Z,PSIE_IN=VE,A_in=a,kick=.true.)
       CALL  COMPX(EL,Z,X,A(1),a(3))
       CALL  COMPy(EL,Z,X,A(2),a(3))
     call E_potential(EL,Z,X,VE)
      if(k%TIME) then
       beta0=el%p%beta0
      else
       beta0=1.0_dp
      endif
      call gen_conv_to_xp(X,a,ve,el%p%exact,beta0,0.0_dp)


   call kill(ti,ve)
    call kill(a)

  end SUBROUTINE conv_to_xppsagan

  SUBROUTINE conv_to_pxrsagan(EL,X,k,z)
    IMPLICIT NONE
    real(dp),INTENT(INOUT):: X(6)
    TYPE(sagan),INTENT(INOUT):: EL
    real(dp) ti,ve,z,a(3),beta0
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

 
 !   call B_E_FIELD(EL,X,Z,PSIE_IN=VE,A_in=a,kick=.true.)
       CALL  COMPX(EL,Z,X,A(1),a(3))
       CALL  COMPy(EL,Z,X,A(2),a(3))
     call E_potential(EL,Z,X,VE)
  if(k%TIME) then
   beta0=el%p%beta0
  else
   beta0=1.0_dp
  endif
    call gen_conv_to_px(X,a,ve,el%p%exact,beta0,0.0_dp)


  end SUBROUTINE conv_to_pxrsagan

  SUBROUTINE conv_to_pxpsagan(EL,X,k,z)
    IMPLICIT NONE
    type(real_8),INTENT(INOUT):: X(6)
    TYPE(saganp),INTENT(INOUT):: EL
    type(real_8) ti,ve,z,a(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
    real(dp) beta0
    call alloc(ti,ve)
    call alloc(a)
 
 
 !   call B_E_FIELD(EL,X,Z,PSIE_IN=VE,A_in=a,kick=.true.)
       CALL  COMPX(EL,Z,X,A(1),a(3))
       CALL  COMPy(EL,Z,X,A(2),a(3))
     call E_potential(EL,Z,X,VE)

  if(k%TIME) then
   beta0=el%p%beta0
  else
   beta0=1.0_dp
  endif
    call gen_conv_to_px(X,a,ve,el%p%exact,beta0,0.0_dp)

    call kill(ti,ve)
    call kill(a)
  end SUBROUTINE conv_to_pxpsagan

end module sagan_WIGGLER
