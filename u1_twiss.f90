module twisspara
  use Mad_like
  implicit none
  integer, private, parameter :: ndd=ndim2,mynreso=20
  integer, private, dimension(4) :: iia,icoast
  integer, private :: NO,ND,ND2,NP,NDPT,NV,icount=0
  real(dp), private, dimension(ndim2,5) :: rdd
  real(dp), private, dimension(ndim2) :: dicu
  real(dp), private, dimension(2,ndim2) :: angp
  character(len=4), private, dimension(4), target :: str4 = (/'1000','0100','0010','0001'/)
  character(len=5), private, dimension(5), parameter :: str5 = (/'10000','01000','00100','00010','00001'/)
  character(len=6), private, dimension(6), target :: str6 = (/'100000','010000','001000','000100','000001','000010'/)
  private zerotwiss,equaltwiss,alloctwiss,killtwiss

  type twiss
     type(damap) a1
     type(damap) a_t
     type(damap) junk
     type(damap) junk1
     type(normalform) n
     logical nf
     real(dp), dimension(3,3) ::  beta,alfa,gama
     real(dp), dimension(3)   ::  mu
     real(dp), dimension(6)   ::  disp
     real(dp), dimension(3)   ::  tune
  end type twiss

  interface assignment (=)
     module procedure equaltwiss
     module procedure zerotwiss
     module procedure normalform_normalform
  end interface

  interface alloc
     module procedure alloctwiss
  end interface

  interface kill
     module procedure killtwiss
  end interface

contains

  subroutine zerotwiss(s1,i)
    implicit none
    type(twiss), intent(inout)::s1
    integer, intent(in)::i

    if(i==0) then

       call liepeek(iia,icoast)
       NO=iia(1)
       ND=iia(3)
       ND2=iia(3)*2
       NP=iia(2)-nd2
       NDPT=icoast(4)
       NV=iia(2)

       s1%a1=1
       s1%a_t=1
       s1%nf=.false.
       s1%beta(:,:)=zero
       s1%alfa(:,:)=zero
       s1%gama(:,:)=zero
       s1%mu(:)=zero
       s1%disp(:)=zero
       s1%tune(:)=zero
       dicu(:)=zero
       angp(:,:)=zero
    endif

  end subroutine zerotwiss

  subroutine equaltwiss(s1,s2)
    implicit none
    type(twiss), intent(inout)::s1
    type(real_8), intent(in)::s2(ndd)
    integer i,j,j1,ii,ii1,ii2,iii,i2,i3
    real(dp) au(6,6),aui(2),sx,cx,dphi(3)
    character(len=nd2), dimension(:), pointer :: string

    icount=icount+1
    if(nd2.eq.4) string=>str4
    if(nd2.eq.6) string=>str6
    if(s2(1)%kind.eq.2) then

       if(.not.s1%nf) then
          s1%nf=.true.
          do i=1,nd2
             s1%junk%v(i)=s2(i)
          enddo
          s1%n=s1%junk
          s1%a1=s1%n%a1
          s1%a_t=s1%n%a_t
          s1%tune(:)=s1%n%tune(:)
          s1%junk=s1%n%a_t
          do i=1,nd2
             dicu(i)=s1%n%a1%v(i).sub.str5(5)
          enddo
       else
          do i=1,nd2
             s1%junk%v(i)=s2(i)
          enddo
          s1%junk=s1%junk*s1%a_t
       endif
    else
       print*,'A Taylor Map is needed as input for the Twiss calculation'
       return
    endif

    if(nd2.eq.4.and.np.ge.1) then
       do i=1,nd2
          do j=1,nd2
             rdd(i,j)=s2(i).sub.str5(j)
          enddo
          rdd(i,5)=s2(i).sub.str5(5)
       enddo
       do j1=1,2
          ii=2*j1
          iii=2*j1-1
          s1%disp(iii)=rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2)+rdd(ii-1,3)*dicu(3)+rdd(ii-1,4)*dicu(4)+rdd(ii-1,5)
          s1%disp(ii)=rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2)+rdd(ii,3)*dicu(3)+rdd(ii,4)*dicu(4)+rdd(ii,5)
       enddo
    endif
    if(nd.eq.3) then
       s1%junk1=s1%junk**(-1)
    endif
    do j=1,nd
       ii=2*j
       ii1=ii-1
       ii2=ii
       if(j.eq.1) then
          i2=4
          i3=6
       elseif(j.eq.2) then
          i2=2
          i3=6
       elseif(j.eq.3) then
          i2=2
          i3=4
          ii1=ii
          ii2=ii-1
       endif

       angp(1,ii-1)=s1%junk%v(ii1).sub.string(ii-1)
       au(ii,ii-1)=s1%junk%v(ii2).sub.string(ii-1)
       angp(1,ii)=s1%junk%v(ii1).sub.string(ii)
       au(ii,ii)=s1%junk%v(ii2).sub.string(ii)
       au(i2-1,i2-1)=s1%junk%v(ii1).sub.string(i2-1)
       au(i2,i2-1)=s1%junk%v(ii2).sub.string(i2-1)
       au(i2-1,i2)=s1%junk%v(ii1).sub.string(i2)
       au(i2,i2)=s1%junk%v(ii2).sub.string(i2)

       s1%beta(1,j)=angp(1,ii-1)*angp(1,ii-1)+angp(1,ii)*angp(1,ii)
       s1%beta(2,j)=au(i2-1,i2-1)*au(i2-1,i2-1)+au(i2-1,i2)*au(i2-1,i2)
       s1%alfa(1,j)=-(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))
       s1%alfa(2,j)=-(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2))
       s1%gama(1,j)=au(ii,ii-1)*au(ii,ii-1)+au(ii,ii)*au(ii,ii)
       s1%gama(2,j)=au(i2,i2-1)*au(i2,i2-1)+au(i2,i2)*au(i2,i2)
       if(nd.eq.3) then
          au(i3-1,i3-1)=s1%junk%v(ii1).sub.string(i3-1)
          au(i3,i3-1)=s1%junk%v(ii2).sub.string(i3-1)
          au(i3-1,i3)=s1%junk%v(ii1).sub.string(i3)
          au(i3,i3)=s1%junk%v(ii2).sub.string(i3)
          s1%beta(3,j)=au(i3-1,i3-1)*au(i3-1,i3-1)+au(i3-1,i3)*au(i3-1,i3)
          s1%alfa(3,j)=-(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3))
          s1%gama(3,j)=au(i3,i3-1)*au(i3,i3-1)+au(i3,i3)*au(i3,i3)
          aui(1)=s1%junk1%v(6).sub.string(6)
          aui(2)=s1%junk1%v(5).sub.string(6)
          if(j.lt.3) then
             s1%disp(ii-1)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
             s1%disp(ii)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
          else
             s1%disp(5)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
             s1%disp(6)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
          endif
       endif
       sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
       cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
       if(abs(sx).gt.c_1d_15.or.abs(cx).gt.c_1d_15) then
          dphi(j)=atan2(sx,cx)*twopii
       else
          dphi(j)=zero
       endif
       s1%mu(j)=s1%mu(j)+dphi(j)
    enddo
    if(icount.gt.1) then
       do i=1,ndim2
          angp(2,i)=angp(1,i)
       enddo
    endif

  end subroutine equaltwiss

  subroutine  alloctwiss(s1)
    implicit none
    type (twiss),intent(inout)::s1

    call alloc(s1%a1)
    call alloc(s1%a_t)
    call alloc(s1%n)
    call alloc(s1%junk)
    call alloc(s1%junk1)

    s1=0
  end subroutine alloctwiss

  subroutine  killtwiss(s1)
    implicit none
    type (twiss),intent(inout)::s1

    call kill(s1%a1)
    call kill(s1%a_t)
    call kill(s1%n)
    call kill(s1%junk)
    call kill(s1%junk1)

    s1%nf=.false.
    s1%beta(:,:)=zero
    s1%alfa(:,:)=zero
    s1%gama(:,:)=zero
    s1%mu(:)=zero
    s1%disp(:)=zero
    s1%tune(:)=zero

  end subroutine killtwiss

  subroutine normalform_normalform(s1,s2)
    implicit none
    type (normalform),intent(inout)::s1
    type (normalform),intent(in)::s2
    integer i,j

    s1%a_t=s2%a_t
    s1%a1=s2%a1
    s1%a%constant(:)=s2%a%constant(:)
    s1%a%Linear=s2%a%Linear
    s1%a%nonlinear=s2%a%nonlinear
    s1%a%pb=s2%a%pb
    s1%normal%constant(:)=s2%normal%constant(:)
    s1%normal%Linear=s2%normal%Linear
    s1%normal%nonlinear=s2%normal%nonlinear
    s1%normal%pb=s2%normal%pb
    s1%DHDJ=s2%DHDJ
    do i=1,ndim
       s1%TUNE(i)=s2%TUNE(i)
       s1%damping(i)=s2%damping(i)
       s1%plane(i)=s2%plane(i)
       do j=1,mynreso
          s1%m(i,j)=s2%m(i,j)
       enddo
    enddo
    s1%nord=s2%nord
    s1%jtune=s2%jtune
    s1%nres=s2%nres
    s1%AUTO=s2%AUTO
  end subroutine normalform_normalform

end module twisspara
