!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis

! This library is free software; you can redistribute it and/or
! modify it under the terms of the The Clarified Artistic License.
! However the source files C_dabnew.f90 and D_Lielib.f90 are derived from
! the LBNL versions of Berz's DA-package and Forest's analysis package.
! Their distribution and commercial usage are governed by the laws of
! the United States of America and more specifically by the USA Department
! of Energy. The above license may be partly or totally void with respect to
! these two files.

! All the other files were created as part of Forest's work as an employee
! of the Japanese Ministry of Culture and Education and are, to best of our
! knowledge, his intellectual property under Japanese Law.

! THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF MERCHANTIBILITY AND
! FITNESS FOR A PARTICULAR PURPOSE.


!  Misalignment module

module rotation_mis
  use dabnew
  use precision_constants
  implicit none
  integer,private,parameter::d=3,max=1000
  private EQUAL,EQUALint,mul,sub,add,uadd,usub
  private mulr,rmul,pow,EQUALd ,factorize_m
  private divr,expmat,norm,expvec,expvec2,texpvec,print_rot

  type matrix_PTC
     real(dp) r(d,d)
     real(dp) t(d)
  end type matrix_PTC

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE EQUAL    !    M2=M1
     MODULE PROCEDURE EQUALINT !    M2=0 (Zeroes M), M=-1,-2,OR -3 (Creates L_1,2,3), M=1 Makes an identity matrix
     MODULE PROCEDURE EQUALD   !    M2=r1 makes a translation of length r in the 3rd direction ( a drift);
  END  INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE MUL     !    M3=M2*M1
     MODULE PROCEDURE MULR    !    M3=r1*M2
     MODULE PROCEDURE RMUL    !    M3=M1*r2
  END INTERFACE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE POW     !  M3=M1**n2 but only n=-1 is accepted (inverse)
  END INTERFACE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE DIVR    ! M3=M1%R/r2  ( Translation part set to zero)
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE SUB     ! M3=M2-M1 (Matrix part only ; Translation part set to zero)
     MODULE PROCEDURE USUB    ! M2=-M1   (Matrix part only ; Translation part set to zero)
  END INTERFACE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE ADD    ! M3=M2+M1 (Matrix part only ; Translation part set to zero)
     MODULE PROCEDURE UADD   ! M2=M1
  END INTERFACE

  INTERFACE EXP
     MODULE PROCEDURE EXPMAT    ! M2=exp(M1)
     MODULE PROCEDURE EXPVEC    ! M2=exp(V1_i L_i)
     MODULE PROCEDURE EXPVEC2   ! M3=exp(-V2_1 L_1) exp(-V2_2 L_2)exp(-V2_3 L_3)exp(V1_i L_i)
     MODULE PROCEDURE TEXPVEC   ! M2= exp(V2_3 L_3) exp(V2_2 L_2) exp(V2_1 L_1)   if IFAC=1
     ! M2= exp(-V2_1 L_1)exp(-V2_2L_2)exp(-V2_3 L_3)   if IFAC=-1
  END INTERFACE

  INTERFACE TEXP
     MODULE PROCEDURE EXPMAT    ! Same routines as above
     MODULE PROCEDURE EXPVEC
     MODULE PROCEDURE EXPVEC2
     MODULE PROCEDURE TEXPVEC
  END INTERFACE


contains

  SUBROUTINE  EQUAL(S2,S1)
    implicit none
    type (matrix_PTC),INTENT(inOUT)::S2
    type (matrix_PTC),INTENT(IN)::S1
    integer i,j

    do i=1,d
       s2%t(i)=s1%t(i)
       do j=1,d
          s2%r(i,j)=s1%r(i,j)
       enddo
    enddo


  END SUBROUTINE EQUAL

  SUBROUTINE  EQUALint(S2,S1)
    implicit none
    type (matrix_PTC),INTENT(inOUT)::S2
    integer,INTENT(IN)::S1
    integer i,j

    do i=1,d
       s2%t(i)=zero
       do j=1,d
          s2%r(i,j)=zero
       enddo
    enddo

    select case(s1)

    case(0)

       do i=1,d
          s2%t(i)=zero
          do j=1,d
             s2%r(i,j)=zero
          enddo
       enddo


    case(1)

       do i=1,d
          s2%r(i,i)=one
       enddo

    case(-1)
       s2%r(2,3)=one
       s2%r(3,2)=-one
    case(-2)
       s2%r(1,3)=-one
       s2%r(3,1)=one
    case(-3)
       s2%r(1,2)=one
       s2%r(2,1)=-one
    case default
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72)))'
       write(w_p%c(1),'(a53,1x,i4)') "Not supported in equalint, rotation_mis module, s1 = ",s1
       CALL WRITE_E(-1)
    end select

  END SUBROUTINE EQUALint

  SUBROUTINE  EQUALd(S2,S1)
    implicit none
    type (matrix_PTC),INTENT(inOUT)::S2
    real(dp),INTENT(IN)::S1

    s2=1
    s2%t(3)=s1

  END SUBROUTINE EQUALd

  function mul(a,b)
    implicit none
    type(matrix_PTC) mul,temp
    type(matrix_PTC), intent(in):: a,b
    integer i,j,k

    temp=0

    do i=1,d
       temp%t(i)=temp%t(i)+a%t(i)
       do j=1,d
          temp%t(i)=temp%t(i)+a%r(i,j)*b%t(j)
          do k=1,d
             temp%r(i,k)=temp%r(i,k)+a%r(i,j)*b%r(j,k)
          enddo
       enddo
    enddo

    mul=temp

  end function mul

  function pow(a,b)
    implicit none
    type(matrix_PTC) pow
    type(matrix_PTC), intent(in):: a
    integer, intent(in):: b
    integer i,j
    real(dp) t(d)

    t=zero
    if(b==-1) then
       do i=1,d
          do j=1,d
             pow%r(i,j)=a%r(j,i)
          enddo
       enddo

       do i=1,d
          do j=1,d
             t(i)=-a%r(j,i)*a%t(j)+t(i)
          enddo
       enddo
       pow%t=t
    else
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72)))'
       w_p%c(1)= "Only inverse accepted in pow(a,b) in rotation module"
       CALL WRITE_E(-1)
    endif

  end function pow

  function rmul(a,b)
    implicit none
    type(matrix_PTC) rmul
    real(dp), intent(in):: a
    type(matrix_PTC), intent(in):: b
    integer i,j

    rmul=0
    do i=1,d
       do j=1,d
          rmul%r(i,j)=a*b%r(i,j)
       enddo
    enddo

  end function rmul

  function mulr(b,a)
    implicit none
    type(matrix_PTC) mulr
    real(dp), intent(in):: a
    type(matrix_PTC), intent(in):: b
    integer i,j

    mulr=0
    do i=1,d
       do j=1,d
          mulr%r(i,j)=a*b%r(i,j)
       enddo
    enddo

  end function mulr

  function divr(b,a)
    implicit none
    type(matrix_PTC) divr
    real(dp), intent(in):: a
    type(matrix_PTC), intent(in):: b
    integer i,j

    divr=0
    do i=1,d
       do j=1,d
          divr%r(i,j)=b%r(i,j)/a
       enddo
    enddo

  end function divr


  function sub(a,b)
    implicit none
    type(matrix_PTC) sub
    type(matrix_PTC), intent(in):: a,b
    integer i,j

    sub=0
    do i=1,d
       do j=1,d
          sub%r(i,j)=a%r(i,j)-b%r(i,j)
       enddo
    enddo

  end function sub


  function add(a,b)
    implicit none
    type(matrix_PTC) add
    type(matrix_PTC), intent(in):: a,b
    integer i,j

    add=0
    do i=1,d
       do j=1,d
          add%r(i,j)=a%r(i,j)+b%r(i,j)
       enddo
    enddo

  end function add

  function uadd(a)
    implicit none
    type(matrix_PTC) uadd
    type(matrix_PTC), intent(in):: a
    uadd=a
  end function uadd

  function usub(a)
    implicit none
    type(matrix_PTC) usub
    type(matrix_PTC), intent(in):: a
    integer i,j

    usub=0
    do i=1,d
       do j=1,d
          usub%r(i,j)=-a%r(i,j)
       enddo
    enddo

  end function usub

  function norm(a)
    implicit none
    real(dp) norm
    type(matrix_PTC), intent(in):: a
    integer i,j

    norm=zero
    do i=1,d
       do j=1,d
          norm=norm+abs(a%r(i,j))
       enddo
    enddo

  end function norm

  FUNCTION expmat( S1)
    implicit none
    TYPE (matrix_PTC) expmat,temp
    TYPE (matrix_PTC), INTENT (IN) :: S1
    integer i
    real(dp) normold,normnew,diffold,diffnew
    logical(lp) leps

    expmat=1
    temp=1
    normold=c_1d10
    normnew=c_1d10
    diffold=c_1d10
    diffnew=c_1d10
    !frs eps_rot_mis1=c_1d_6
    leps=.true.

    do i=1,max
       temp=temp*s1/REAL(i,kind=DP)
       expmat=expmat+temp
       normnew=norm(expmat)
       diffnew=abs(normnew-normold)
       if(leps) then
          if(diffnew<eps_rot_mis1) leps=.false.
       else
          if(diffnew>=diffold) goto 100
       endif
       diffold=diffnew
       normold=normnew
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='(1((1X,A72)))'
    w_p%c(1) =  "never converged in expmat, module rotation_mis "
    CALL WRITE_E(-1)

100 continue

  end FUNCTION expmat

  FUNCTION expvec2( S1,s2)
    implicit none
    TYPE (matrix_PTC) expvec2,L(d),o
    real(dp), INTENT (IN) :: S1(d),s2(d)
    integer i

    o=0
    do i=1,d
       L(i)=-i
       o=L(i)*S1(i)+o
       L(i)=s2(i)*L(i)
    enddo

    expvec2=0
    expvec2=exp(-L(1))*exp(-L(2))*exp(-L(3))*exp(o)

  end FUNCTION expvec2

  FUNCTION expvec( S1)
    implicit none
    TYPE (matrix_PTC) expvec,L(d),o
    real(dp), INTENT (IN) :: S1(d)
    integer i

    o=0
    do i=1,d
       L(i)=-i
       o=L(i)*S1(i)+o
    enddo

    expvec=0
    expvec=exp(o)

  end FUNCTION expvec

  FUNCTION texpvec( S1,ifac)
    implicit none
    TYPE (matrix_PTC) texpvec,o
    real(dp), INTENT (IN) :: S1(d)
    integer, INTENT (IN) :: ifac
    real(dp) t(3)
    integer i

    o=1
    if(ifac>0) then
       do i=1,d
          t=zero;t(i)=s1(i);
          o=exp(t)*o
       enddo
    else
       do i=d,1,-1
          t=zero;t(i)=-s1(i);
          o=exp(t)*o
       enddo
    endif

    texpvec=0
    texpvec=o

  end FUNCTION texpvec

  subroutine FACTORIZE_ROTATION(mis,l,alpha,t_in,beta_in,t_out,beta_out)
    implicit none
    real(dp), INTENT (IN) :: mis(2*d),L,alpha
    real(dp), INTENT (inout) :: t_in(d),beta_in(d), t_out(d),beta_out(d)
    type(matrix_PTC) a,ai

    type(matrix_PTC) Y,T_D,r,ri(d)
    real(dp) t(d),beta(d),an(d)
    integer i,dir
    ! direction of propagation
    DIR=1
    an=zero;an(2)=-alpha/two  ! PTC convention
    do i=1,3
       t(i)=mis(i)
       beta(i)=mis(i+d)
    enddo

    Y=exp(an)
    t_d=L/2

    do i=1,d
       an=zero;an(i)=beta(i)*(-1)**(i+1); ! PTC convention
       ri(i)=exp(an)
    enddo
    r=1
    r=ri(3)*ri(2)*ri(1)
    r%t=t

    !       CALL ROT_XZ(  C%LOCAL%ALPHA/two  ,X,  C%MAG%P%BETA0,OUR,C%MAG%P%TIME)   !MAKE MAGNET  THICK AGAIN   1
    !       CALL TRANSZ(  C%LOCAL%L/two      ,X,  C%MAG%P%BETA0,OUR,C%MAG%P%TIME)    !MAKE MAGNET THICK AGAIN   2
    !       CALL ROT_YZ(C%MAG%r(1),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)                ! rotations
    !       CALL ROT_XZ(C%MAG%r(2),X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)
    !       CALL ROT_XY(C%MAG%r(3),X,OU)
    !       CALL TRANS(C%MAG%d,X,C%MAG%P%BETA0,OU,C%MAG%P%TIME)                       !translation

    !       CALL TRANSZ( -C%LOCAL%L/two      ,X,  C%MAG%P%BETA0,OUR,C%MAG%P%TIME)    !MAKE MAGNET THIN        3
    !      CALL ROT_XZ( -C%LOCAL%ALPHA/two  ,X,  C%MAG%P%BETA0,OUR,C%MAG%P%TIME)       !MAKE MAGNET THIN        4
    a= y**(-1)*t_d**(-1)*R*t_d*y        ! PTC representation
    ai= y*t_d*R**(-1)*t_d**(-1)*y**(-1) ! PTC representation
    if(dir==1) then
       call factorize_m2( a,beta_in)
       t_in=a%t
       call factorize_m2( ai,beta_out)
       t_out=ai%t
       beta_in(2)=-beta_in(2)     ! PTC convention
       beta_out(2)=-beta_out(2)   ! PTC convention
       ! call print_rot(y,6)
    else  ! NOT USED ANYMORE
       a = a**(-1)
       ai= ai**(-1)
       call factorize_m2( ai,beta_in)
       t_in=ai%t
       call factorize_m2( a,beta_out)
       t_out=a%t
       beta_in(2) =-beta_in(2)     ! PTC convention
       beta_out(2)=-beta_out(2)   ! PTC convention
       ! dir=-1 symmetries involving z->-z
       t_in(3)=-t_in(3);t_out(3)=-t_out(3);
       beta_in(2)=-beta_in(2);beta_in(1)=-beta_in(1)
       beta_out(2)=-beta_out(2);beta_out(1)=-beta_out(1)
    endif
  end subroutine FACTORIZE_ROTATION

  subroutine factorize_m( S1,s2)
    implicit none
    TYPE (matrix_PTC) L(d),o,id
    real(dp), INTENT (IN) :: S1(d)
    real(dp), INTENT (out) :: S2(d)
    real(dp) t(d)
    integer i
    real(dp) normold,normnew,diffold,diffnew
    logical(lp) leps
    do i=1,d
       L(i)=-i
    enddo

    leps=.true.
    id=1

    t=s1

    normold=c_1d10
    normnew=c_1d10
    diffold=c_1d10
    diffnew=c_1d10
    !frs eps_rot_mis2=c_1d_9

    do i=1,max

       o=texp(s1,t)

       o=o-id
       t(1)=t(1)+o%r(2,3)
       t(2)=t(2)+o%r(3,1)
       t(3)=t(3)+o%r(1,2)

       normnew=abs(t(1))+abs(t(2))+abs(t(3))
       diffnew=abs(normnew-normold)
       if(leps) then
          if(diffnew<eps_rot_mis2) leps=.false.
       else
          if(diffnew>=diffold) goto 100
       endif
       diffold=diffnew
       normold=normnew
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='(1((1X,A72)))'
    w_p%c(1) =  "never converged in factorize_m, module rotation_mis "
    CALL WRITE_E(-1)

100 continue


    s2=t
  end subroutine factorize_m

  subroutine factorize_m2( s1,s2)
    implicit none
    TYPE (matrix_PTC) L(d),id,o,s1i
    TYPE (matrix_PTC), INTENT (INOUT) :: s1
    real(dp), INTENT (INout) :: S2(d)
    real(dp) t(d),tz(d),tt(d),re(d),m(3,3)
    integer i,j,ier,method
    real(dp) normold,normnew,diffold,diffnew
    logical(lp) leps
    ! THIS ROUTINE FACTORIZE A ROTATION IN THE STANDARD ptc ORDER
    ! NOTICE S2(2)=-A_Y

    do i=1,d
       L(i)=-i
    enddo

    leps=.true.
    id=1

    t=zero
    tz=zero
    s1i=s1**(-1)
    normold=c_1d10
    normnew=c_1d10
    diffold=c_1d10
    diffnew=c_1d10
    !frs eps_rot_mis2=c_1d_9
    !frs epsdif=c_1d_6

    method=1
    do i=1,max

       o=texp(tz,t)*s1

       if(method==0) then
          o=o-id
          re(1)=o%r(2,3)
          re(2)=o%r(3,1)
          re(3)=o%r(1,2)

          do j=1,3
             tt=t;tt(j)=tt(j)+epsdif;
             o=texp(tz,tt)*s1
             o=o-id
             m(1,j)=(o%r(2,3)-re(1))/epsdif
             m(2,j)=(o%r(3,1)-re(2))/epsdif
             m(3,j)=(o%r(1,2)-re(3))/epsdif
          enddo
          call matinv(m,m,3,3,ier)
          do j=1,3
             t(j)=t(j)-m(j,1)*re(1)-m(j,2)*re(2)-m(j,3)*re(3)
          enddo
       else
          o=o-id
          t(1)=t(1)+o%r(2,3)
          t(2)=t(2)+o%r(3,1)
          t(3)=t(3)+o%r(1,2)
       endif
       normnew=abs(t(1))+abs(t(2))+abs(t(3))
       diffnew=abs(normnew-normold)
       if(leps) then
          if(diffnew<eps_rot_mis2) leps=.false.
       else
          if(diffnew>=diffold) goto 100
       endif
       diffold=diffnew
       normold=normnew
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='(1((1X,A72)))'
    w_p%c(1) =  "never converged in expmat, module rotation_mis "
    CALL WRITE_E(-1)

100 continue


    s2=t
  end subroutine factorize_m2

  subroutine print_rot( b,s2)
    implicit none
    type(matrix_PTC),intent(inout) ::  b
    integer,intent (in)::  s2

    integer i
    write(s2,*) "%%%%%%%     %%%%%%%"
    write(s2,*) b%t(1),b%t(2),b%t(3)

    do i=1,d
       write(s2,*) b%r(i,1),b%r(i,2),b%r(i,3)
    enddo

    return
  end subroutine print_rot

end module rotation_mis
