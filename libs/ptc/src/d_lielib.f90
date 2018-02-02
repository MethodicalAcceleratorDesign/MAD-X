! The Full Polymorphic Package
! The module in this file is, to the best of our knowledge,
! the property of Lawrence Berkeley National Laboratory
! Its distribution and commercial usage may therefore be governed by the laws of the
! United States of America

module lielib_yang_berz
  use dabnew
  !  use precision_constants
  implicit none
  public
  !  private
  PUBLIC XGAM,XGBM  !,FILTRES
  PUBLIC LIEPEEK,INITPERT,HYPER,MAPFLOL
  PUBLIC ETALL1,TAKE,ETALL,DAPEK0,ETINI,DACLRD,DACOPD,DIFD
  PUBLIC INTD,ETCCT,TRXFLO,TRX,FACFLOD,EXPFLO,DALIND,ETINV
  PUBLIC INPUTRES,MAPNORMF,DHDJFLO,GETTURA   !,SETIDPRIDPRSET,
  PUBLIC FLOFACG,FLOFAC,DACMUD,CTORFLO,RTOCFLO,CTOR,RTOC,ETPIN
  PUBLIC LIEINIT,PERTPEEK,FLOWPARA,COMCFU
  PUBLIC DAPOK0,FACFLO,EXPFLOD,gofix
  public getcct,GETINV,gtrx,eig6
  private DFILT,DLIE,FILT,REXT,respoke
  private etallnom,simil
  private dapokzer,davar0,taked,daread,daprid,daflo,daflod,fexpo,etcom,etpoi
  private exp1d,expnd2,liefact,mapnorm,orderflo,nuanaflo,h2pluflo,rotflo,rotiflo
  private ctord,rtocd,resvec,reelflo,midbflo,mulnd2,movearou,movemul,cpart
  private ctoi,itoc,etrtc,etctr,etcjg,etdiv,sympl3,ety,etyt,ety2  !,flip
  integer,public,parameter::ndim=4,nreso=100
  integer,public::no,nv,nd,nd2,ndpt
  integer, private :: ndc,ndc2,ndt,iref,itu,iflow,jtune,nres !,idpr
  integer, private,dimension(ndim)::nplane,idsta,ista
!  real(dp), private,dimension(0:20)::xintex
  real(dp),dimension(0:20)::xintex

  real(dp), private,dimension(ndim)::dsta,sta,angle,rad,ps,rads
  real(dp), private,dimension(ndim,nreso)::mx
  !real(dp),private::epsplane
  !real(dp),private,dimension(ndim)::xplane
  !integer,public,parameter::ndim2=2*ndim,ntt=40
  integer,private,parameter::ndim2=2*ndim,ntt=lnv    ! joahn 2008
  !  integer,private,parameter::ndim2=2*ndim,ntt=100
  character(120), private :: line
  logical :: frankheader=.true.
  logical :: new_ndpt = .true.
  integer,private :: nt_pos,npt_pos
  integer time_plane
  real(dp), private :: stmem(ndim)
  logical(lp) :: courant_snyder=.true.
  logical(lp) :: check_krein=.true.
  real(dp) :: size_krein=1.e-10_dp

  real(dp),dimension(ndim2)::rr_eigen,ri_eigen
contains




  subroutine lieinit(no1,nv1,nd1,ndpt1,time_pos,da_init)   !,nis
    implicit none
    !! Lieinit initializes AD Package and Lielib
    integer i,nd1,ndc1,ndpt1,no1,nv1      !,nis
    real(dp),dimension(ndim)::ang,ra,st
    integer, optional :: time_pos
    logical, optional :: da_init
    logical dai
    dai=.true.
    if(present(da_init)) dai=da_init
    
    if(present(time_pos)) then
       time_plane=time_pos
    else
       time_plane=0
       if(ndpt1==0.and.(nd1==3.or.nd1==4) ) time_plane=3
    endif

    do i=1,ndim
       nplane(i)=2*i-1
       ang(i)=0.0_dp
       ra(i)=0.0_dp
       st(i)=1.0_dp
    enddo
    no=no1
    nv=nv1
    nd=nd1
    nd2=2*nd1
    !    do i=1,100
    !       is(i)=0
    !    enddo
    if(dai) call daini(no,nv,0)
    !    if(nis.gt.0)call etallnom(is,nis,'$$IS      ')

    if(ndpt1.eq.0) then
       ndpt=0
       ndt=0
       ndc1=0
    else
       if(.not.new_ndpt) then
          ndpt=ndpt1
          ndc1=1
          if(ndpt.eq.nd2) then
             ndt=nd2-1
          else
             ndt=nd2
             if(ndpt.ne.nd2-1) then
                line=' LETHAL ERROR IN LIEINIT'
                write(6,*) line
                stop
             endif
          endif
       else
          if(mod(ndpt1,2)==0) then
             ndpt=nd2
             ndt=nd2-1
          else
             ndpt=nd2-1
             ndt=nd2
          endif
          npt_pos=ndpt1
          ndc1=1
          if(mod(npt_pos,2)==0) then
             nt_pos=npt_pos-1
          else
             nt_pos=npt_pos+1
          endif
          if(npt_pos<3.or.npt_pos>nd2) then
             line=' LETHAL ERROR IN LIEINIT'
             write(6,*) line
             stop
          endif
       endif
    endif
    ndc=ndc1
    ndc2=2*ndc1
    iref=0
    call initpert(st,ang,ra)
    !    iref=iref1
    !    if(iref1.eq.0) then
    itu=0
    !    else
    !       itu=1
    !    endif
    ! if(iref1.eq.0)
    iref=-1

    if(lielib_print(1)==1) then
       write(6,'(a17,4(1x,i4))') ' no,nv,nd,ndpt = ',no1,nv1,nd1,ndpt1
    endif

    do i=0,20
       xintex(i)=0.0_dp
    enddo
    xintex(0) =1.0_dp
    xintex(1) =0.5_dp
    xintex(2) =1.0_dp/12.0_dp
    xintex(4) =-1.0_dp/720e0_dp
    xintex(6) =1.0_dp/30240e0_dp
    xintex(8) =-1.0_dp/1209600e0_dp
    xintex(10)=1.0_dp/21772800e0_dp


    return
  end subroutine lieinit
  subroutine flowpara(ifl,jtu)
    implicit none
    integer ifl,jtu

    iflow=ifl
    jtune=jtu
    return
  end subroutine flowpara
  subroutine pertpeek(st,ang,ra)
    implicit none
    integer i
    real(dp),dimension(ndim)::ang,ra,st

    do i=1,nd
       st(i)=sta(i)
       ang(i)=angle(i)
       ra(i)=rad(i)
    enddo
    return
  end subroutine pertpeek
  subroutine inputres(mx1,nres1)
    implicit none
    integer i,j,nres1
    integer,dimension(ndim,nreso)::mx1

    nres=nres1
    do i=1,nreso
       do j=1,ndim
          mx(j,i)=0
       enddo
    enddo

    do i=1,nres
       do j=1,ndim
          mx(j,i)=mx1(j,i)
       enddo
    enddo
    return
  end subroutine inputres
  subroutine respoke(mres,nre,ire)
    implicit none
    integer i,ire,j,nre
    integer,dimension(ndim,nreso)::mres
    real(dp),dimension(ndim)::ang,ra,st

    iref=ire
    nres=nre
    do j=1,nreso
       do i=1,nd
          mx(i,j)=mres(i,j)
       enddo
    enddo
    call initpert(st,ang,ra)
    return
  end subroutine respoke
  subroutine liepeek(iia,icoast)
    implicit none
    integer,dimension(:)::iia,icoast

    iia(1)=no
    iia(2)=nv
    iia(3)=nd
    iia(4)=nd2
    icoast(1)=ndc
    icoast(2)=ndc2
    icoast(3)=ndt
    icoast(4)=ndpt

    return
  end subroutine liepeek


  subroutine etallnom(x,n)
    implicit none
    ! CREATES A AD-VARIABLE WHICH CAN BE DESTROYED BY DADAL
    ! allocates vector of n polynomials and give it the name NOM=A10
    integer i,n,nd2
    integer,dimension(:)::x
    integer,dimension(4)::i1,i2
    !    character(10) nom

    do i=1,iabs(n)
       x(i)=0
       call daall0(x(i))
    enddo
    !    call daallno(x,iabs(n),nom)
    if(n.lt.0) then
       call liepeek(i1,i2)
       nd2=i1(4)
       do i=nd2+1,-n
          call davar(x(i),0.0_dp,i)
       enddo
    endif
    return
  end subroutine etallnom
  subroutine etall(x,n)
    implicit none
    ! allocates vector of n polynomials
    integer i,n,nd2
    integer,dimension(:)::x
    integer,dimension(4)::i1,i2
    do i=1,iabs(n)
       x(i)=0
    enddo

    if(.not.frankheader) then
       do i=1,iabs(n)
          call daall0(x(i))
       enddo
    else
       do i=1,iabs(n)
          call daall1(x(i),'etall     ',no,nv)
       enddo
    endif
    if(n.lt.0) then
       call liepeek(i1,i2)
       nd2=i1(4)
       do i=nd2+1,-n
          call davar(x(i),0.0_dp,i)
       enddo
    endif
    return
  end subroutine etall
  subroutine etall1(x)
    implicit none
    integer x

    x=0
    if(.not.frankheader) then
       call daall0(x)
    else
       call daall1(x,'etall     ',no,nv)
    endif
    return
  end subroutine etall1

  subroutine etcct(x,y,z)
    implicit none
    !  Z=XoY
    integer i,nt
    integer,dimension(ntt)::ie,iv
    integer,dimension(:)::x,y,z
    if(.not.c_%stable_da) return

    nt=nv-nd2
    if(nt.gt.0) then
       call etallnom(ie,nt) !,'IE        ')
       do i=nd2+1,nv
          call davar(ie(i-nd2),0.0_dp,i)
       enddo
       do i=nd2+1,nv
          iv(i)=ie(i-nd2)
       enddo
    endif
    do i=1,nd2
       iv(i)=y(i)
    enddo
    call dacct(x,nd2,iv,nv,z,nd2)
    if(nt.gt.0) then
       call dadal(ie,nt)
    endif
    return
  end subroutine etcct

  subroutine getcct(x,y,z,n)
    implicit none
    !  Z=XoY
    integer i,nt,n
    integer,dimension(ntt)::ie,iv
    integer,dimension(:)::x,y,z
    if(.not.c_%stable_da) return

    nt=nv-n
    if(nt.gt.0) then
       call etallnom(ie,nt)  !,'IE        ')
       do i=n+1,nv
          call davar(ie(i-n),0.0_dp,i)
       enddo
       do i=n+1,nv
          iv(i)=ie(i-n)
       enddo
    endif
    do i=1,n
       iv(i)=y(i)
    enddo
    call dacct(x,n,iv,nv,z,n)
    if(nt.gt.0) then
       call dadal(ie,nt)
    endif
    return
  end subroutine getcct


  subroutine trx(h,rh,y)
    implicit none
    !  :RH: = Y :H: Y^-1 =  :HoY:
    integer i,nt,h,rh
    integer,dimension(ntt)::ie,iv
    integer,dimension(1)::h1,rh1
    integer,dimension(:)::y
    if(.not.c_%stable_da) return

    nt=nv-nd2
    if(nt.gt.0) then
       call etallnom(ie,nt)  !,'IE        ')
       do i=nd2+1,nv
          call davar(ie(i-nd2),0.0_dp,i)
       enddo
       do i=nd2+1,nv
          iv(i)=ie(i-nd2)
       enddo
    endif
    do i=1,nd2
       iv(i)=y(i)
    enddo
    h1(1)=h
    rh1(1)=rh
    call dacct(h1,1,iv,nv,rh1,1)
    if(nt.gt.0) then
       call dadal(ie,nt)
    endif
    return
  end subroutine trx

  subroutine gtrx(h,rh,y,n)
    implicit none
    !  :RH: = Y :H: Y^-1 =  :HoY:
    integer i,nt,h,rh,n
    integer,dimension(ntt)::ie,iv
    integer,dimension(1)::h1,rh1
    integer,dimension(:)::y
    if(.not.c_%stable_da) return

    nt=nv-n
    if(nt.gt.0) then
       call etallnom(ie,nt)  !,'IE        ')
       do i=n+1,nv
          call davar(ie(i-n),0.0_dp,i)
       enddo
       do i=n+1,nv
          iv(i)=ie(i-n)
       enddo
    endif
    do i=1,n
       iv(i)=y(i)
    enddo
    h1(1)=h
    rh1(1)=rh
    call dacct(h1,1,iv,nv,rh1,1)
    if(nt.gt.0) then
       call dadal(ie,nt)
    endif
    return
  end subroutine gtrx

  subroutine trxflo(h,rh,y)
    implicit none
    !  *RH* = Y *H* Y^-1  CHANGE OF A VECTOR FLOW OPERATOR
    integer j,k,b1,b2
    integer,dimension(:)::h,rh,y
    integer,dimension(ndim2)::yi,ht
    if(.not.c_%stable_da) return

    call etallnom(yi,nd2)  !  ,'YI        ')
    call etallnom(ht,nd2)  !  ,'HT        ')
    call etall1(b1 )
    call etall1(b2 )

    call etinv(y,yi)
    !----- HT= H o Y
    call etcct(h,y,ht)
    !----
    call daclrd(rh)
    do j=1,nd2
       do k=1,nd2
          call dader(k,yi(j),b1)
          call trx(b1,b2,y)
          call damul(b2,ht(k),b1)
          call daadd(b1,rh(j),b2)
          call dacop(b2,rh(j))
       enddo
    enddo

    call dadal1(b2)
    call dadal1(b1)
    call dadal(ht,nd2)
    call dadal(yi,nd2)
    return
  end subroutine trxflo

  subroutine simil(a,x,ai,y)
    implicit none
    !  Y= AoXoAI
    integer,dimension(:)::x,y,a,ai
    integer,dimension(ndim2)::w,v
    if(.not.c_%stable_da) return

    call etallnom(w,nd2) ! ,'W         ')
    call etallnom(v,nd2) ! ,'V         ')

    call etcct(a,x,w)
    call etcct(w,ai,v)

    call dacopd(v,y)

    call dadal(v,nd2)
    call dadal(w,nd2)
    return
  end subroutine simil

  subroutine etini(x)
    implicit none
    !  X=IDENTITY
    integer i
    integer,dimension(:)::x
    if(.not.c_%stable_da) return

    do i=1,nd2
       call davar(x(i),0.0_dp,i)
    enddo
    return
  end subroutine etini

  subroutine etinv(x,y)
    implicit none
    ! Y=X^-1
    integer i,nt
    integer,dimension(ntt)::ie1,ie2,iv1,iv2
    integer,dimension(:)::x,y
    if(.not.c_%stable_da) return
    nt=nv-nd2
    if(nt.gt.0) then
       do i=1,nt
          ie1(i)=0
          ie2(i)=0
       enddo
       call etallnom(ie1,nt) !,'IE1       ')
       call etallnom(ie2,nt) !,'IE2       ')
       do i=nd2+1,nv
          call davar(ie1(i-nd2),0.0_dp,i)
       enddo
       do i=nd2+1,nv
          iv1(i)=ie1(i-nd2)
          iv2(i)=ie2(i-nd2)
       enddo
    endif
    do i=1,nd2
       iv1(i)=x(i)
       iv2(i)=y(i)
    enddo

    call dainv(iv1,nv,iv2,nv)
    if(nt.gt.0) then
       call dadal(ie2,nt)
       call dadal(ie1,nt)
    endif
    return
  end subroutine etinv

  subroutine etpin(x,y,jj)
    implicit none

    integer i,nt
    integer,dimension(ntt)::ie1,ie2,iv1,iv2
    integer,dimension(:)::x,y,jj
    if(.not.c_%stable_da) return

    nt=nv-nd2
    if(nt.gt.0) then
       do i=1,nt
          ie1(i)=0
          ie2(i)=0
       enddo
       call etallnom(ie1,nt) !,'IE1       ')
       call etallnom(ie2,nt) !,'IE2       ')
       do i=nd2+1,nv
          call davar(ie1(i-nd2),0.0_dp,i)
       enddo
       do i=nd2+1,nv
          iv1(i)=ie1(i-nd2)
          iv2(i)=ie2(i-nd2)
       enddo
    endif
    do i=1,nd2
       iv1(i)=x(i)
       iv2(i)=y(i)
    enddo

    call dapin(iv1,nv,iv2,nv,jj)
    if(nt.gt.0) then
       call dadal(ie2,nt)
       call dadal(ie1,nt)
    endif
    return
  end subroutine etpin

  subroutine getinv(x,y,n)
    implicit none
    ! Y=X^-1
    integer i,nt,n
    integer,dimension(ntt)::ie1,ie2,iv1,iv2
    integer,dimension(:)::x,y
    if(.not.c_%stable_da) return


    nt=nv-n
    if(nt.gt.0) then
       do i=1,nt
          ie1(i)=0
          ie2(i)=0
       enddo
       call etallnom(ie1,nt) !,'IE1       ')
       call etallnom(ie2,nt) !,'IE2       ')
       do i=n+1,nv
          call davar(ie1(i-n),0.0_dp,i)
       enddo
       do i=n+1,nv
          iv1(i)=ie1(i-n)
          iv2(i)=ie2(i-n)
       enddo
    endif
    do i=1,n
       iv1(i)=x(i)
       iv2(i)=y(i)
    enddo

    call dainv(iv1,nv,iv2,nv)
    if(nt.gt.0) then
       call dadal(ie2,nt)
       call dadal(ie1,nt)
    endif
    return
  end subroutine getinv

  subroutine dapek0(v,x,jj)
    implicit none

    integer i,jj
    integer,dimension(ntt)::jd
    integer,dimension(:)::v
    real(dp),dimension(:)::x
    if(.not.c_%stable_da) return

    do i=1,ntt
       jd(i)=0
    enddo
    do i=1,jj
       call dapek(v(i),jd,x(i))
    enddo
    return
  end subroutine dapek0

  subroutine dapok0(v,x,jj)
    implicit none
    integer i,jj
    integer,dimension(ntt)::jd
    integer,dimension(:)::v
    real(dp),dimension(:)::x
    if(.not.c_%stable_da) return

    do i=1,ntt
       jd(i)=0
    enddo
    do i=1,jj
       call dapok(v(i),jd,x(i))
    enddo
    return
  end subroutine dapok0

  subroutine dapokzer(v,jj)
    implicit none
    integer i,jj
    integer,dimension(:)::v
    integer,dimension(ntt)::jd
    if(.not.c_%stable_da) return

    do i=1,ntt
       jd(i)=0
    enddo
    do i=1,jj
       call dapok(v(i),jd,0.0_dp)
    enddo
    return
  end subroutine dapokzer

  subroutine davar0(v,x,jj)
    implicit none
    integer i,jj
    integer,dimension(:)::v
    real(dp),dimension(:)::x
    if(.not.c_%stable_da) return

    do i=1,jj
       call davar(v(i),x(i),i)
    enddo
    return
  end subroutine davar0

  subroutine comcfu(b,f1,f2,c)
    implicit none
    ! Complex dacfu
    integer,dimension(2)::b,c
    integer,dimension(4)::t
    real(dp),external::f1,f2

    if(.not.c_%stable_da) return



    call etall(t,4)

    call dacfu(b(1),f1,t(1))
    call dacfu(b(1),f2,t(2))
    call dacfu(b(2),f1,t(3))
    call dacfu(b(2),f2,t(4))

    call dasub(t(1),t(4),c(1))
    call daadd(t(2),t(3),c(2))
    call dadal(t,4)

    return
  end subroutine comcfu

  subroutine take(h,m,ht)
    implicit none
    !  HT= H_M  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
    integer i,m,h,ht,b1,b2,b3
    integer,dimension(ntt)::j
    real(dp) r
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etall1(b2)
    call etall1(b3)

    if(no.ge.2) then
       if(m.eq.0) then
          do i=1,ntt
             j(i)=0
          enddo
          call dapek(h,j,r)
          call dacon(ht,r)
       else
          !          call danot(m)
          !          call dacop(h,b1)
          call datrunc(h,m+1,b1)
          !          call danot(m-1)
          !          call dacop(b1,b2)
          call datrunc(b1,m,b2)
          !          call danot(no)
          call dasub(b1,b2,b3)
          call dacop(b3,ht)
       endif
    else
       do i=1,ntt
          j(i)=0
       enddo
       if(m.eq.0) then
          call dapek(h,j,r)
          call dacon(ht,r)
       elseif(m.eq.1)  then
          do i=1,nv
             j(i)=1
             call dapek(h,j,r)
             call dapok(b3,j,r)
             j(i)=0
          enddo
          call dacop(b3,ht)
       else
          call daclr(ht)
       endif
    endif

    call dadal1(b3)
    call dadal1(b2)
    call dadal1(b1)
    return
  end subroutine take

  subroutine taked(h,m,ht)
    implicit none
    !  \VEC{HT}= \VEC{H_M}  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
    integer i,m,b1,b2
    integer,dimension(:)::h,ht
    integer,dimension(ntt)::j
    integer,dimension(ndim2)::x
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etall1(b2)
    call etallnom(x,nd2) !  ,'X         ')


    do i=1,ntt
       j(i)=0
    enddo

    do   i=1,nd2
       call take(h(i),m,ht(i))
    enddo
    call dadal(x,nd2)
    call dadal1(b2)
    call dadal1(b1)
    return
  end subroutine taked

  subroutine daclrd(h)
    implicit none
    ! clear a map : a vector of nd2 polynomials
    integer i
    integer,dimension(:)::h
    if(.not.c_%stable_da) return

    do i=1,nd2
       call daclr(h(i))
    enddo
    return
  end subroutine daclrd

  subroutine dacopd(h,ht)
    implicit none
    !    H goes into HT  (nd2 array)
    integer i
    integer,dimension(:)::h,ht
    if(.not.c_%stable_da) return

    do i=1,nd2
       call dacop(h(i),ht(i))
    enddo
    return
  end subroutine dacopd

  subroutine datruncd(h,io,ht)
    implicit none
    !    H goes into HT  (nd2 array)
    integer i
    integer,dimension(:)::h,ht
    integer io
    if(.not.c_%stable_da) return

    do i=1,nd2
       call datrunc(h(i),io,ht(i))
    enddo
    return
  end subroutine datruncd

  subroutine dacmud(h,sca,ht)
    implicit none
    integer i
    integer,dimension(:)::h,ht
    real(dp) sca
    if(.not.c_%stable_da) return

    do i=1,nd2
       call dacmu(h(i),sca,ht(i))
    enddo
    return
  end subroutine dacmud
  subroutine dalind(h,rh,ht,rt,hr)
    implicit none
    integer i
    integer,dimension(:)::h,ht,hr
    integer,dimension(ndim2)::b
    real(dp) rh,rt
    if(.not.c_%stable_da) return

    call etallnom(b,nd2) !  ,'B         ')

    do i=1,nd2
       call dalin(h(i),rh,ht(i),rt,b(i))
    enddo
    call dacopd(b,hr)
    call dadal(b,nd2)
    return
  end subroutine dalind
  subroutine daread(h,nd1,mfile,xipo)
    implicit none
    !  read a map
    integer i,mfile,nd1
    integer,dimension(:)::h
    integer,dimension(ntt)::j
    real(dp) rx,xipo
    if(.not.c_%stable_da) return

    do i=1,ntt
       j(i)=0
    enddo
    do i=1,nd1
       call darea(h(i),mfile)
       call dapek(h(i),j,rx)
       rx=rx*xipo
       call dapok(h(i),j,rx)
    enddo
    return
  end subroutine daread
  subroutine daprid(h,n1,n2,mfile)
    implicit none
    !  print a map
    integer i,mfile,n1,n2
    integer,dimension(:)::h
    if(.not.c_%stable_da) return

    if(mfile.le.0) return
    do i=n1,n2
       call dapri(h(i),mfile)
    enddo
    return
  end subroutine daprid

  subroutine daflo(h,x,y)
    implicit none
    ! LIE EXPONENT ROUTINES WITH FLOW OPERATORS
    !
    !     \VEC{H}.GRAD X =Y
    integer i,x,y,b1,b2,b3
    integer,dimension(:)::h
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etall1(b2)
    call etall1(b3)

    call daclr(b1)
    call daclr(b2)
    do i=1,nd2
       call dader(i,x,b2)
       call damul(b2,h(i),b3)
       call daadd(b3,b1,b2)
       call dacop(b2,b1)
    enddo
    call dacop(b1,y)
    call dadal1(b3)
    call dadal1(b2)
    call dadal1(b1)
    return
  end subroutine daflo
  subroutine daflod(h,x,y)
    implicit none
    integer i
    integer,dimension(:)::h,x,y
    integer,dimension(ndim2)::b1,b2
    if(.not.c_%stable_da) return

    call etall(b1,nd2)
    call etall(b2,nd2)

    call dacopd(h,b1)
    call dacopd(x,b2)

    do i=1,nd2
       call daflo(b1,b2(i),y(i))
    enddo

    call dadal(b1,nd2)
    call dadal(b2,nd2)
    return
  end subroutine daflod

  subroutine intd(v,h,sca)
    implicit none
    ! IF SCA=-one
    !     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:
    !
    ! IF SCA=one
    !     \VEC{V}.GRAD  = GRAD H . GRAD
    integer i,h,b1,b2,b3,b4
    integer,dimension(:)::v
    integer,dimension(ndim2)::x
    real(dp) sca
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etall1(b2)
    call etall1(b3)
    call etall1(b4)
    call etallnom(x,nd2) !  ,'X         ')

    call daclr(b4)
    call daclr(h)
    call etini(x)
    do i=1,nd
       call dacfu(v(2*i-1),dlie,b3)
       call dacfu(v(2*i),dlie,b1)
       call damul(b1,x(2*i-1),b2)
       call damul(b3,x(2*i),b1)
       call dalin(b2,1.0_dp,b1,sca,b3)
       call daadd(b3,b4,b2)
       call dacop(b2,b4)
    enddo
    call dacop(b4,h)
    call dadal(x,nd2)
    call dadal1(b4)
    call dadal1(b3)
    call dadal1(b2)
    call dadal1(b1)
    return
  end subroutine intd

  subroutine difd(h1,v,sca)
    implicit none
    ! INVERSE OF INTD ROUTINE
    integer i,h1,b1,h
    integer,dimension(:)::v
    real(dp) sca
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etall1(h)
    call dacop(h1,h)
    do i=1,nd
       call dader(2*i-1,h,v(2*i))
       call dader(2*i,h,b1)
       call   dacmu(b1,sca,v(2*i-1))
    enddo
    call dadal1(h)
    call dadal1(b1)
    return
  end subroutine difd

  subroutine expflo(h,x,y,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) X = Y
    logical(lp) more
    integer i,nrmax,x,y,b1,b2,b3,b4
    integer,dimension(:)::h
    real(dp) coe,eps,r,rbefore
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etall1(b2)
    call etall1(b3)
    call etall1(b4)

    call dacop(x,b4)
    call dacop(x,b1)
    more=.true.
    rbefore=1e30_dp
    do i=1,nrmax
       coe=1.0_dp/REAL(i,kind=DP)
       call dacmu(b1,coe,b2)
       call daflo(h,b2,b1)
       call daadd(b4,b1,b3)
       call daabs(b1,r)
       if(more) then
          if(r.gt.eps) then
             rbefore=r
             goto 100
          else
             rbefore=r
             more=.false.
          endif
       else
          if(r.ge.rbefore) then
             call dacop(b3,y)
             call dadal1(b4)
             call dadal1(b3)
             call dadal1(b2)
             call dadal1(b1)
             return
          endif
          rbefore=r
       endif
100    continue
       call dacop(b3,b4)
    enddo
    if(lielib_print(2)==1) then
       write(6,'(a6,1x,G21.14,1x,a25)') ' NORM ',eps,' NEVER REACHED IN EXPFLO '
    endif
    call dacop(b3,y)
    call dadal1(b4)
    call dadal1(b3)
    call dadal1(b2)
    call dadal1(b1)
    return
  end subroutine expflo

  subroutine expflod(h,x,w,eps,nrmax)
    implicit none
    ! DOES EXP( \VEC{H} ) \VEC{X} = \VEC{Y}
    integer j,nrmax,b0
    integer,dimension(:)::x,w,h
    integer,dimension(ndim2)::v
    real(dp) eps
    if(.not.c_%stable_da) return

    call etall1(b0 )
    call etallnom(v,nd2) !  ,'V         ')

    call dacopd(x,v)
    do j=1,nd2
       call expflo(h,v(j),b0,eps,nrmax)
       call dacop(b0,v(j))
    enddo
    call dacopd(v,w)
    call dadal(v,nd2)
    call dadal1(b0)
    return
  end subroutine expflod
  subroutine facflo(h,x,w,nrmin,nrmax,sca,ifac)
    implicit none
    ! IFAC=1
    ! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX ) X= Y
    ! IFAC=-1
    ! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) X= Y
    integer i,ifac,nmax,nrmax,nrmin,x,w,v
    integer,dimension(:)::h
    integer,dimension(ndim2)::bm,b0
    real(dp) eps,sca
    if(.not.c_%stable_da) return

    call etallnom(bm,nd2) !  ,'BM        ')
    call etallnom(b0,nd2) !  ,'B0        ')
    call etall1(v)

    call dacop(x,v)

    !    eps=-one
    !    call daeps(eps)
    eps=epsflo
    nmax=100
    !
    ! IFAC =1 ---> V = EXP(:SCA*H(NRMAX):)...EXP(:SCA*H(NRMIN):)X
    if(ifac.eq.1) then
       do i=nrmax,nrmin,-1
          call taked(h,i,b0)
          call dacmud(b0,sca,bm)

          call expflo(bm,v,b0(1),eps,nmax)
          call dacop(b0(1),v)
       enddo
    else
       ! IFAC =-1 ---> V = EXP(:SCA*H(NRMIN):)...EXP(:SCA*H(NRMAX):)X
       do i=nrmin,nrmax
          call taked(h,i,b0)
          call dacmud(b0,sca,bm)

          call expflo(bm,v,b0(1),eps,nmax)
          call dacop(b0(1),v)
       enddo
    endif
    call dacop(v,w)
    call dadal1(v)
    call dadal(b0,nd2)
    call dadal(bm,nd2)
    return
  end subroutine facflo
  subroutine facflod(h,x,w,nrmin,nrmax,sca,ifac)
    implicit none
    ! IFAC=1
    ! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX )  \VEC{X}= \VEC{Y}
    ! IFAC=-1
    ! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) \VEC{X}= \VEC{Y}
    integer i,ifac,nrmax,nrmin
    integer,dimension(:)::x,w,h
    real(dp) sca
    if(.not.c_%stable_da) return

    do i=1,nd2
       call facflo(h,x(i),w(i),nrmin,nrmax,sca,ifac)
    enddo

    return
  end subroutine facflod
  subroutine fexpo(h,x,w,nrmin,nrmax,sca,ifac)
    implicit none
    !   WRAPPED ROUTINES FOR THE OPERATOR  \VEC{H}=:H:
    ! WRAPPING FACFLOD
    integer ifac,nrma,nrmax,nrmi,nrmin,h
    integer,dimension(:)::x,w
    integer,dimension(ndim2)::v
    real(dp) sca
    if(.not.c_%stable_da) return

    nrmi=nrmin-1
    nrma=nrmax-1
    call etall(v,nd2)
    call difd(h,v,-1.0_dp)
    call facflod(v,x,w,nrmi,nrma,sca,ifac)

    call dadal(v,nd2)

    return
  end subroutine fexpo
  subroutine etcom(x,y,h)
    implicit none
    ! ETCOM TAKES THE BRACKET OF TWO VECTOR FIELDS.
    integer i,j,t1,t2
    integer,dimension(:)::h,x,y
    integer,dimension(ndim2)::t3
    if(.not.c_%stable_da) return

    call etall1(t1)
    call etall1(t2)
    call etall(t3,nd2)

    do j=1,nd2
       do i=1,nd2

          call dader(i,x(j),t1)
          call dader(i,y(j),t2)
          call damul(x(i),t2,t2)
          call damul(y(i),t1,t1)
          call dalin(t2,1.0_dp,t1,-1.0_dp,t1)
          call daadd(t1,t3(j),t3(j))

       enddo
    enddo

    call dacopd(t3,h)

    call dadal1(t1)
    call dadal1(t2)
    call dadal(t3,nd2)
    return
  end subroutine etcom
  subroutine etpoi(x,y,h)
    implicit none
    ! ETPOI TAKES THE POISSON BRACKET OF TWO FUNCTIONS
    integer i,h,x,y,t1,t2,t3
    if(.not.c_%stable_da) return

    call etall1(t1)
    call etall1(t2)
    call etall1(t3)

    do i=1,nd

       call dader(2*i-1,x,t1)
       call dader(2*i,y,t2)
       call damul(t1,t2,t1)

       call dalin(t1,1.0_dp,t3,1.0_dp,t3)
       call dader(2*i-1,y,t1)
       call dader(2*i,x,t2)
       call damul(t1,t2,t1)

       call dalin(t1,-1.0_dp,t3,1.0_dp,t3)

    enddo

    call dacop(t3,h)

    call dadal1(t1)
    call dadal1(t2)
    call dadal1(t3)
    return
  end subroutine etpoi
  subroutine exp1d(h,x,y,eps,non)
    implicit none
    ! WRAPPING EXPFLO
    integer non,h,x,y
    integer,dimension(ndim2)::v
    real(dp) eps
    if(.not.c_%stable_da) return

    call etall(v,nd2)
    call difd(h,v,-1.0_dp)
    call expflo(v,x,y,eps,non)

    call dadal(v,nd2)

    return
  end subroutine exp1d
  subroutine expnd2(h,x,w,eps,nrmax)
    implicit none
    ! WRAPPING EXPFLOD USING EXP1D
    integer j,nrmax,b0,h
    integer,dimension(:)::x,w
    integer,dimension(ndim2)::v
    real(dp) eps
    if(.not.c_%stable_da) return

    call etall1(b0)
    call etallnom(v,nd2) !  ,'V         ')

    call dacopd(x,v)
    do j=1,nd2
       call exp1d(h,v(j),b0,eps,nrmax)
       call dacop(b0,v(j))
    enddo
    call dacopd(v,w)
    call dadal(v,nd2)
    call dadal1(b0)
    return
  end subroutine expnd2

  subroutine flofacg(xy,h,epsone)
    implicit none
    ! GENERAL ONE EXPONENT FACTORIZATION
    logical(lp) more
    integer i,k,kk,nrmax
    integer,dimension(:)::xy,h
    integer,dimension(ndim2)::x,v,w,t,z
    integer,dimension(ntt)::jj
    real(dp) eps,epsone,r,xn,xnbefore,xnorm,xnorm1,xx
    if(.not.c_%stable_da) return

    jj(1)=1
    !
    call etallnom(v,nd2) !  ,'V         ')
    call etallnom(w,nd2) !  ,'W         ')
    call etallnom(t,nd2) !  ,'T         ')
    call etallnom(x,nd2) !  ,'Z         ')
    call etallnom(z,nd2) !  ,'Z         ')

    call etini(v)
    call daclrd(w)
    xnorm1=0.0_dp
    do i=1,nd2
       call daabs(xy(i),r)
       xnorm1=xnorm1+r
    enddo
    xnbefore=1e36_dp
    more=.false.
    eps=1e-5_dp
    nrmax=1000
    xn=1e4_dp

    if(epsone>0.0_dp) then  !epsone>zero
       do k=1,nrmax
          call dacmud(h,-1.0_dp,t)
          call expflod(t,xy,x,eps,nrmax)
          call dalind(x,1.0_dp,v,-1.0_dp,t)
          ! write(20,*) "$$$$$$$$$$$$$$",k,"$$$$$$$$$$$$$$$$$$$$"
          ! call daprid(t,1,1,20)
          if(xn.lt.epsone) then
             if(lielib_print(3)==1) then

                write(6,'(a14,g21.14)') " xn quadratic ",xn

                ! CALL !WRITE_a
             endif
             call daflod(t,t,w)
             call dalind(t,1.0_dp,w,-0.5_dp,t)
             call dacopd(t,z)
             call dacopd(t,w)
             !  second order in W
             call etcom(h,w,x)
             call etcom(x,w,x)
             !  END OF  order in W

             do kk=1,3   !10
                call etcom(h,w,w)
                call dalind(z,1.0_dp,w,xintex(kk),z)
             enddo
             call dacopd(z,t)
             xx=1.0_dp/12.0_dp
             call dalind(x,xx,h,1.0_dp,h)
          endif

          call dalind(t,1.0_dp,h,1.0_dp,h)
          xnorm=0.0_dp
          do i=1,nd2
             call daabs(t(i),r)
             xnorm=xnorm+r
          enddo
          xn=xnorm/xnorm1
          if(xn.ge.epsone.and.(lielib_print(3)==1)) then
             write(6,'(a11,g21.14)') " xn linear ",xn
          endif
          if(xn.lt.eps.or.more) then
             more=.true.
             if(xn.ge.xnbefore) goto 1000
             xnbefore=xn
          endif
       enddo
1000   continue
    else  !epsone>zero
       do k=1,nint(abs(epsone))-1
          call dacmud(h,-1.0_dp,t)
          call expflod(t,xy,x,eps,nrmax)
          call dalind(x,1.0_dp,v,-1.0_dp,t)
          ! write(20,*) "$$$$$$$$$$$$$$",k,"$$$$$$$$$$$$$$$$$$$$"
          ! call daprid(t,1,1,20)

          call dalind(t,1.0_dp,h,1.0_dp,h)
       enddo
    endif
    if(lielib_print(3)==1) WRITE(6,*) " K ", K,epsone
    if(lielib_print(3)==1) then
       write(6,'(a11,i4)') " iteration " , k
    endif
    !  if(lielib_print(3)==1) CALL WRITE_a
    call dadal(x,nd2)
    call dadal(w,nd2)
    call dadal(v,nd2)
    call dadal(t,nd2)
    call dadal(z,nd2)
    return
  end subroutine flofacg

  subroutine flofac(xy,x,h)
    implicit none
    ! GENERAL DRAGT-FINN FACTORIZATION
    integer k
    integer,dimension(:)::xy,x,h
    integer,dimension(ndim2)::v,w
    if(.not.c_%stable_da) return

    call etallnom(v,nd2) !  ,'V         ')
    call etallnom(w,nd2) !  ,'W         ')

    call dacopd(xy,x)
    !    call dacopd(x,v)
    call datruncd(x,2,v)
    call daclrd(w)
    !    call danot(1)
    call etinv(v,w)
    !    call danot(no)
    call etcct(x,w,v)
    !    call danot(1)
    !    call dacopd(xy,x)
    !    call datruncd(xy,1,x)  ! lethal error
    call datruncd(xy,2,x)


    !    call danot(no)
    call dacopd(v,w)
    call daclrd(h)
    do k=2,no
       call taked(w,k,v)
       call dalind(v,1.0_dp,h,1.0_dp,h)
       call facflod(h,w,v,k,k,-1.0_dp,-1)
       call dacopd(v,w)
    enddo
    call dadal(w,nd2)
    call dadal(v,nd2)
    return
  end subroutine flofac
  subroutine liefact(xy,x,h)
    implicit none
    ! SYMPLECTIC DRAGT-FINN FACTORIZATION WRAPPING FLOFAC
    integer h
    integer,dimension(:)::xy,x
    integer,dimension(ndim2)::v
    if(.not.c_%stable_da) return

    call etall(v,nd2)

    call flofac(xy,x,v)
    call intd(v,h,-1.0_dp)
    !
    call dadal(v,nd2)

    return
  end subroutine liefact
  subroutine mapnorm(x,ft,a2,a1,xy,h,nord)
    implicit none
    !--NORMALIZATION ROUTINES OF LIELIB
    !- WRAPPING MAPNORMF
    integer isi,nord,ft,h
    integer,dimension(:)::x,a1,a2,xy
    integer,dimension(ndim2)::hf,ftf
    if(.not.c_%stable_da) return

    call etall(ftf,nd2)
    call etall(hf,nd2)
    isi=0
    call mapnormf(x,ftf,a2,a1,xy,hf,nord,isi)
    call intd(hf,h,-1.0_dp)
    call intd(ftf,ft,-1.0_dp)
    call dadal(ftf,nd2)
    call dadal(hf,nd2)

    return
  end subroutine mapnorm

  subroutine gettura(psq,radsq)
    implicit none
    integer ik,j1,j2
    real(dp),dimension(ndim)::psq,radsq,st
    if(.not.c_%stable_da) return

    if(new_ndpt) then
       if(ndpt/=0) then
          if(mod(ndpt,2)==0) then
             j1=ndpt/2
             j2=npt_pos/2
          else
             j1=(ndpt+1)/2
             j2=(npt_pos+1)/2
          endif
       endif
       do ik=1,nd
          psq(ik)=ps(ik)
          radsq(ik)=rads(ik)
          st(ik)=stmem(ik)
       enddo
       if(ndpt/=0) then
          psq(j2)   = ps(j1)
          psq(j1)   = ps(j2)
          radsq(j2) = rads(j1)
          radsq(j1) = rads(j2)
          st(j2) = stmem(j1)
          st(j1) = stmem(j2)
       endif
       do ik=1,nd
          if(ik/=time_plane) then
             if(st(ik)+1e-3_dp.gt.1.0_dp.and.psq(ik).lt.0.0_dp) psq(ik)=psq(ik)+1.0_dp
          endif
       enddo
       if(time_plane>0) then
          if(st(time_plane)+1e-3_dp.gt.1.0_dp.and.psq(time_plane).lt.-0.5_dp) psq(time_plane)=psq(time_plane)+1.0_dp
       endif
    else
       do ik=1,nd
          psq(ik)=ps(ik)
          radsq(ik)=rads(ik)
       enddo
    endif
    return
  end subroutine gettura

  subroutine setidpr(nplan)
    implicit none
    integer ik
    integer,dimension(ndim)::nplan
    if(.not.c_%stable_da) return

    do ik=1,nd
       nplane(ik)=nplan(ik)
    enddo
    return
  end subroutine setidpr
  !  subroutine idprset(idprint)
  !    implicit none
  !    integer idprint
  !    if(.not.c_%stable_da) return!
  !
  !    idpr=idprint
  !
  !    return
  !  end subroutine idprset
  !   CALL MAPNORMF(junk%V%i,s2%a%nonlinear%v%i,s2%a%linear%v%i,&
  !        & s2%A1%v%i,s2%normal%linear%v%i,s2%normal%nonlinear%v%i,s2%NORD,s2%jtune)
  subroutine mapnormf(x,ft,a2,a1,xy,h,nord,isi)
    implicit none
    integer ij,isi,nord
    integer,dimension(ndim2)::a1i,a2i
    integer,dimension(:)::x,a1,a2,ft,xy,h
    real(dp),dimension(ndim)::angle,rad,st,p
    if(.not.c_%stable_da) return

    call etallnom(a1i,nd2) !  ,'A1I       ')
    call etallnom(a2i,nd2) !  ,'A2I       ')
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
    call dacopd(x,xy)
    ! go to fix point in the parameters + pt to order nord>=1
    if(nv>nd2.or.ndc==1) then
       call gofix(xy,a1,a1i,nord)
       call simil(a1i,xy,a1,xy)
    else    !  this "if" was added to remove crashes when y-=plane is nearly identity
       call etini(a1)   ! in stochastic kick calculations
       call etini(a1i)  ! 2002.10.20
    endif
    ! linear part
    call midbflo(xy,a2,a2i,angle,rad,st)
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
    call simil(a2i,xy,a2,xy)
    call dacopd(xy,a2i)
    call orderflo(h,ft,xy,angle,rad)
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
    call h2pluflo(h,p,rad)
    !      CALL TAKED(A2I,1,XY)
    call taked(a2i,1,a1i)
    call etcct(xy,a1i,xy)


    call dadal(a2i,nd2)
    call dadal(a1i,nd2)
    return
  end subroutine mapnormf

  subroutine get_flip_info(nt_pos1,npt_pos1)
    implicit none
    integer nt_pos1,npt_pos1
    nt_pos1=nt_pos
    npt_pos1=npt_pos
  end   subroutine get_flip_info

  subroutine flip(xy,xyf)
    implicit none

    integer,dimension(:):: xy,xyf
    integer,dimension(ndim2)::x,xi
    if(.not.c_%stable_da) return

    if(nt_pos>=nd2-1) return


    call etallnom(x,nd2)
    call etallnom(xi,nd2)

    call etini(x)

    call davar(x(npt_pos),0.0_dp,ndpt)
    call davar(x(nt_pos),0.0_dp,ndt)
    call davar(x(ndpt),0.0_dp,npt_pos)
    call davar(x(ndt),0.0_dp,nt_pos)
    call etinv(x,xi)

    call simil(x,xy,xi,xyf)

    !       call davar(x(ndpt),zero,ndpt)


    call dadal(xi,nd2)
    call dadal(x,nd2)

  end subroutine flip

  subroutine flip_real_array(xy,xyf,i)
    implicit none
    integer i
    real(dp),dimension(:):: xy,xyf
    real(dp),dimension(ndim)::x
    if(.not.c_%stable_da) return

    if(nt_pos>=nd2-1) return
    x=0.0_dp
    x(1:nd) = xy(1:nd)
    xyf(1:nd) = xy(1:nd)
    if(mod(ndpt,2)==0) then
       if(i==1) then
          xyf(ndpt/2)=x(npt_pos/2)
          xyf(npt_pos/2)=x(ndpt/2)
       else
          xyf(npt_pos/2)=x(ndpt/2)
          xyf(ndpt/2)=x(npt_pos/2)
       endif
    else
       if(i==1) then
          xyf((ndpt+1)/2)=x((npt_pos+1)/2)
          xyf((npt_pos+1)/2)=x((ndpt+1)/2)
       else
          xyf((npt_pos+1)/2)=x((ndpt+1)/2)
          xyf((ndpt+1)/2)=x((npt_pos+1)/2)
       endif
    endif


  end subroutine flip_real_array

  subroutine flip_resonance(xy,xyf,i)
    implicit none
    integer i,j
    integer,dimension(:,:):: xy,xyf
    integer,dimension(NDIM,NRESO) ::x
    if(.not.c_%stable_da) return

    if(nt_pos>=nd2-1) return
    x=0
    x = xy
    xyf  = xy
    if(mod(ndpt,2)==0) then
       do j=1,nreso
          if(i==1) then
             xyf(ndpt/2,j)=x(npt_pos/2,j)
             xyf(npt_pos/2,j)=x(ndpt/2,j)
          else
             xyf(npt_pos/2,j)=x(ndpt/2,j)
             xyf(ndpt/2,j)=x(npt_pos/2,j)
          endif
       enddo
    else
       do j=1,nreso
          if(i==1) then
             xyf((ndpt+1)/2,j)=x((npt_pos+1)/2,j)
             xyf((npt_pos+1)/2,j)=x((ndpt+1)/2,j)
          else
             xyf((npt_pos+1)/2,j)=x((ndpt+1)/2,j)
             xyf((ndpt+1)/2,j)=x((npt_pos+1)/2,j)
          endif
       enddo
    endif


  end subroutine flip_resonance

  subroutine flipflo(xy,xyf,i)
    implicit none
    integer i
    integer,dimension(:):: xy,xyf
    integer,dimension(ndim2)::x,xi
    if(.not.c_%stable_da) return

    if(nt_pos>=nd2-1) return


    call etallnom(x,nd2)
    call etallnom(xi,nd2)

    call etini(x)

    call davar(x(npt_pos),0.0_dp,ndpt)
    call davar(x(nt_pos),0.0_dp,ndt)
    call davar(x(ndpt),0.0_dp,npt_pos)
    call davar(x(ndt),0.0_dp,nt_pos)
    call etinv(x,xi)

    if(i==1) then
       call trxflo(xy,xyf,xi)
    else
       call trxflo(xy,xyf,x)
    endif


    call dadal(xi,nd2)
    call dadal(x,nd2)

  end subroutine flipflo

  subroutine flip_i(xy,xyf,i)
    implicit none
    integer i
    integer  xy,xyf
    integer,dimension(ndim2)::x,xi
    if(.not.c_%stable_da) return

    if(nt_pos>=nd2-1) return


    call etallnom(x,nd2)
    call etallnom(xi,nd2)

    call etini(x)

    call davar(x(npt_pos),0.0_dp,ndpt)
    call davar(x(nt_pos),0.0_dp,ndt)
    call davar(x(ndpt),0.0_dp,npt_pos)
    call davar(x(ndt),0.0_dp,nt_pos)
    call etinv(x,xi)

    if(i==1) then
       call trx(xy,xyf,xi)
    else
       call trx(xy,xyf,x)
    endif


    call dadal(xi,nd2)
    call dadal(x,nd2)

  end subroutine flip_i


  subroutine gofix(xy,a1,a1i,nord)
    implicit none
    ! GETTING TO THE FIXED POINT AND CHANGING TIME APPROPRIATELY IN THE
    ! COASTING BEAM CASE
    !****************************************************************
    ! X = A1 XY A1I WHERE X IS TO THE FIXED POINT TO ORDER NORD
    ! for ndpt not zero, works in all cases. (coasting beam: eigenvalue
    !1 in Jordan form)
    !****************************************************************
    integer i,nord
    integer,dimension(:)::xy,a1,a1i
    integer,dimension(ndim2)::x,w,v,rel
    real(dp) xic
    if(.not.c_%stable_da) return



    call etallnom(x,nd2) !  ,  'X         ')
    call etallnom(w,nd2) !  ,  'W         ')
    call etallnom(v,nd2) !  ,  'V         ')
    call etallnom(rel,nd2) !  ,'REL       ')


    ! COMPUTATION OF A1 AND A1I USING DAINV



    call etini(rel)

    !    call danot(nord)

    call etini(v)

    do i=1,nd2-ndc2
       !       call dacop(xy(i),x(i))
       call datrunc(xy(i),nord+1,x(i))
       call dalin(x(i),1.0_dp,rel(i),-1.0_dp,v(i))
    enddo
    call etinv(v,w)
    call datruncd(w,nord+1,w)
    call daclrd(x)
    if(ndc.eq.1) then
       call davar(x(ndpt),0.0_dp,ndpt)
    endif
    call etcct(w,x,v)
    if(ndc.eq.1) then
       call daclr(v(nd2))
       call daclr(v(nd2-ndc))
    endif
    call dalind(rel,1.0_dp,v,1.0_dp,a1)
    call dalind(rel,1.0_dp,v,-1.0_dp,a1i)

    if(ndpt.ne.0) then

       !  CORRECTIONS
       call daclrd(w)
       call daclrd(v)
       call daclrd(x)

       do i=1,nd2-ndc2
          call dalin(a1(i),1.0_dp,rel(i),-1.0_dp,w(i))
       enddo

       !      COMPUTE Deta/Ddelta
       call dacopd(w,a1)

       do i=1,nd2-ndc2
          call dader(ndpt,w(i),w(i))
       enddo
       !      COMPUTE J*Deta/dDELTA

       do i=1,nd-ndc
          call dacmu(w(2*i),1.0_dp,v(2*i-1) )
          call dacmu(w(2*i-1),-1.0_dp,v(2*i) )
       enddo

       xic=(-1)**(ndt)

       do i=1,nd2-ndc2
          call damul(v(i),rel(i),x(1))
          call daadd(x(1),w(ndt),w(ndt))
          call dacop(a1(i),w(i))
       enddo
       call dacmu(w(ndt),xic,w(ndt))

       call expflod(w,rel,a1,1e-7_dp,10000)
       ! END OF  CORRECTIONS

       call datruncd(a1,nord+1,a1)
       call etinv(a1,a1i)
       call datruncd(a1i,nord+1,a1i)
    endif


    !    call danot(no)

    call dadal(rel,nd2)
    call dadal(v,nd2)
    call dadal(w,nd2)
    call dadal(x,nd2)
    return
  end subroutine gofix


  subroutine orderflo(h,ft,x,ang,ra)
    implicit none
    !-   NONLINEAR NORMALIZATION PIECE OF MAPNORMF
    integer k
    integer,dimension(ndim2)::w,v,rel,roi,b1,b5,b6,b9
    integer,dimension(:)::x,h,ft
    real(dp),dimension(ndim)::ang,ra
    if(.not.c_%stable_da) return

    call etallnom(w,nd2) !  ,'W         ')
    call etallnom(v,nd2) !  ,'V         ')
    call etallnom(rel,nd2) !  ,'REL       ')
    call etallnom(roi,nd2) !  ,'ROI       ')
    call etallnom(b1,nd2) !  ,'B1        ')
    call etallnom(b5,nd2) !  ,'B5        ')
    call etallnom(b6,nd2) !  ,'B6        ')
    call etallnom(b9,nd2) !  ,'B9        ')
    call rotiflo(roi,ang,ra)
    call etini(rel)
    call daclrd(h)
    call daclrd(ft)
    call etcct(x,roi,x)
    do k=2,no
       ! IF K>2 V = H(K)^-1 X(K)
       call facflod(h,x,v,2,k-1,-1.0_dp,-1)
       ! EXTRACTING K TH DEGREE OF V ----> W
       call taked(v,k,w)
       !  write(16,*) "$$$$$$$$  K  $$$$$$$$$$", k
       ! W = EXP(B5) + ...
       call dacopd(w,b5)
       !      CALL INTD(W,B5,-one)
       ! B5 ON EXIT IS THE NEW CONTRIBUTION TO H
       ! B6 IS THE NEW CONTRIBUTION TO FT
       call nuanaflo(b5,b6)
       call dalind(b5,1.0_dp,h,1.0_dp,b1)
       call dacopd(b1,h)
       ! EXP(B9) = EXP( : ROTI B6 :)
       call trxflo(b6,b9,roi)

       ! V = EXP(-B6) REL
       call facflod(b6,rel,v,k,k,-1.0_dp,1)
       ! W = V o X
       call etcct(v,x,w)
       if(lielib_print(5)==1) then

          write(6,'(a13,i4)') ' ORDERFLO K= ', k

       endif
       ! X = EXP(B9) W
       call facflod(b9,w,x,k,k,1.0_dp,1)
       ! B6 IS THE NEW CONTRIBUTION TO FT
       call dalind(b6,1.0_dp,ft,1.0_dp,b1)
       call dacopd(b1,ft)
    enddo
    call dadal(b9,nd2)
    call dadal(b6,nd2)
    call dadal(b5,nd2)
    call dadal(b1,nd2)
    call dadal(roi,nd2)
    call dadal(rel,nd2)
    call dadal(v,nd2)
    call dadal(w,nd2)
    return
  end subroutine orderflo
  subroutine nuanaflo(h,ft)
    implicit none
    ! RESONANCE DENOMINATOR OPERATOR (1-R^-1)^-1
    integer i
    integer,dimension(:)::h,ft
    integer,dimension(ndim2)::br,bi,c,ci
    integer,dimension(2)::t1,t2
    if(.not.c_%stable_da) return

    call etall(br,nd2)
    call etall(bi,nd2)
    call etall(c,nd2)
    call etall(ci,nd2)

    call ctorflo(h,br,bi)

    ! FILTERING RESONANCES AND TUNE SHIFTS
    ! ASSUMING REALITY I.E. B(2*I-1)=CMPCJG(B(2*I))

    do i=1,nd2
       iflow=i
       call dacfu(br(i),filt,c(i))
       call dacfu(bi(i),filt,ci(i))
    enddo
    call rtocflo(c,ci,h)

    do i=1,nd2

       iflow=i
       call dacfu(br(i),dfilt,br(i))
       call dacfu(bi(i),dfilt,bi(i))
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
    enddo

    call rtocflo(c,ci,ft)

    call dadal(br,nd2)
    call dadal(bi,nd2)
    call dadal(c,nd2)
    call dadal(ci,nd2)

    return
  end subroutine nuanaflo

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

  subroutine dhdjflo(h,t)
    implicit none
    ! CONVENIENT TUNE SHIFT FINDED FOR SYMPLECTIC CASE (NU,DL)(H)=T
    integer i,bb1,bb2,j1,j2
    integer,dimension(:)::h,t
    integer,dimension(ndim2)::b1,b2,temp
 
    if(.not.c_%stable_da) return



    call etall(b1,nd2)
    call etall(b2,nd2)
    call etall1(bb1)
    call etall1(bb2)

    call ctorflo(h,b1,b2)

    do i=1,nd-ndc
       call datra(2*i,b2(2*i),bb1)
       call dacmu(bb1,twopii,t(i+nd))
       call dacop(t(i+nd),bb1)
       call daclr(bb2)
       call rtoc(bb1,bb2,bb1)
       call dacop(bb1,t(i))
    enddo

    if(ndpt.ne.0) then
       call dacop(h(ndt),t(nd))
       call dacop(b1(ndt),t(nd2))
    endif

    call dadal1(bb2)
    call dadal1(bb1)
    call dadal(b2,nd2)
    call dadal(b1,nd2)
    return
  end subroutine dhdjflo



  subroutine h2pluflo(h,ang,ra)
    implicit none
    ! POKES IN \VEC{H}  ANGLES AND DAMPING COEFFFICIENTS
    !
    integer i
    integer,dimension(ntt)::j
    integer,dimension(:)::h
    real(dp) r1,r2
    real(dp),dimension(ndim)::ang,ra,st
    if(.not.c_%stable_da) return

    do i=1,nd
       st(i)=2.0_dp*sta(i)-1.0_dp
    enddo

    do i=1,ntt
       j(i)=0
    enddo

    do i=1,nd-ndc
       j(2*i-1)=1
       r1=-ang(i)
       !-----
       call dapok(h(2*i),j,r1)

       r2=ra(i)
       call dapok(h(2*i-1),j,r2)
       j(2*i-1)=0

       j(2*i)=1
       r1=ang(i)*st(i)
       call dapok(h(2*i-1),j,r1)
       call dapok(h(2*i),j,r2)
       j(2*i)=0

    enddo

    if(ndpt.eq.nd2-1) then
       j(ndpt)=1
       call dapok(h(ndt),j,ang(nd))
    elseif(ndpt.eq.nd2) then
       j(ndpt)=1
       call dapok(h(ndt),j,ang(nd))     !!!! correct 2014.2.12 with David Sagan and Chris Mayes
    endif
    return
  end subroutine h2pluflo

  subroutine rotflo(ro,ang,ra)
    implicit none
    ! CREATES R AND R^-1 USING THE EXISTING ANGLES AND DAMPING
    ! COULD BE REPLACED BY A CALL H2PLUFLO FOLLOWED BY EXPFLOD
    ! CREATES R
    integer i
    integer,dimension(ntt)::j
    integer,dimension(:)::ro
    real(dp) ch,sh,sim,xx
    real(dp),dimension(ndim)::co,si,ang,ra
    if(.not.c_%stable_da) return

    call daclrd(ro)
    do i=1,nd-ndc
       xx=EXP(ra(i))
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
       call dapok(ro(2*i-1),j,co(i))
       call dapok(ro(2*i),j,sim)
       j(2*i-1)=0
       j(2*i)=1
       call dapok(ro(2*i),j,co(i))
       call dapok(ro(2*i-1),j,si(i))
       j(2*i)=0
    enddo

    if(ndc.eq.1) then
       j(ndt)=1
       call dapok(ro(ndt),j,1.0_dp)
       call dapok(ro(ndpt),j,0.0_dp)
       j(ndt)=0
       j(ndpt)=1
       call dapok(ro(ndt),j,ang(nd))
       call dapok(ro(ndpt),j,1.0_dp)
       j(ndpt)=0
    endif

    return
  end subroutine rotflo
  subroutine rotiflo(roi,ang,ra)
    implicit none
    ! CREATES  R^-1
    integer i
    real(dp) ch,sh,sim,simv,xx
    integer,dimension(ntt)::j
    integer,dimension(:)::roi
    real(dp),dimension(ndim)::co,si,ang,ra
    if(.not.c_%stable_da) return

    !    do i=1,10
    j=0
    !    enddo

    call daclrd(roi)
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
       call dapok(roi(2*i-1),j,co(i))
       simv=-sim
       call dapok(roi(2*i),j,simv)
       j(2*i-1)=0
       j(2*i)=1
       simv=-si(i)
       call dapok(roi(2*i),j,co(i))
       call dapok(roi(2*i-1),j,simv)
       j(2*i)=0
    enddo

    if(ndc.eq.1) then
       j(ndt)=1
       call dapok(roi(ndt),j,1.0_dp)
       call dapok(roi(ndpt),j,0.0_dp)
       j(ndt)=0
       j(ndpt)=1
       call dapok(roi(ndt),j,-ang(nd))
       call dapok(roi(ndpt),j,1.0_dp)
       j(ndpt)=0
    endif

    return
  end subroutine rotiflo

  subroutine hyper(a,ch,sh)
    implicit none
    real(dp) a,ch,sh,x,xi
    if(.not.c_%stable_da) return
    !   USED IN ROTIFLO AND ROTFLO
    x=EXP(a)
    xi=1.0_dp/x
    ch=(x+xi)/2.0_dp
    sh=(x-xi)/2.0_dp
    return
  end subroutine hyper

  subroutine ctor(c1,r2,i2)
    implicit none
    ! CHANGES OF BASIS
    !   C1------> R2+I R1
    integer c1,r2,i2,b1,b2
    integer,dimension(ndim2)::x

    if(.not.c_%stable_da) return



    call etall1(b1)
    call etall1(b2)
    call etallnom(x,nd2) !  ,'X         ')


    call ctoi(c1,b1)
    call etcjg(x)
    call trx(b1,b2,x)
    call dalin(b1,0.5_dp,b2,0.5_dp,r2)
    call dalin(b1,0.5_dp,b2,-0.5_dp,i2)



    call dadal(x,nd2)
    call dadal1(b2)
    call dadal1(b1)
    return
  end subroutine ctor
  subroutine rtoc(r1,i1,c2)
    implicit none
    !  INVERSE OF CTOR
    integer c2,r1,i1,b1

    if(.not.c_%stable_da) return



    call etall1(b1)

    call daadd(r1,i1,b1)
    call itoc(b1,c2)
    call dadal1(b1)




    return
  end subroutine rtoc
  subroutine ctorflo(c,dr,di)
    implicit none
    ! FLOW CTOR
    integer,dimension(:)::dr,di,c
    if(.not.c_%stable_da) return

    call ctord(c,dr,di)
    call resvec(dr,di)

    return
  end subroutine ctorflo
  subroutine rtocflo(dr,di,c)
    implicit none
    ! FLOW RTOC
    integer,dimension(:)::dr,di,c
    integer,dimension(ndim2)::er,ei
    if(.not.c_%stable_da) return

    call etall(er,nd2)
    call etall(ei,nd2)

    call reelflo(dr,di,er,ei)
    call rtocd(er,ei,c)

    call dadal(er,nd2)
    call dadal(ei,nd2)

    return
  end subroutine rtocflo
  subroutine ctord(c,cr,ci)
    implicit none
    ! ROUTINES USED IN THE INTERMEDIATE STEPS OF CTORFLO AND RTOCFLO
    ! SAME AS CTOR  OVER ARRAYS CONTAINING ND2 COMPONENTS
    ! ROUTINE USEFUL IN INTERMEDIATE FLOW CHANGE OF BASIS
    integer i
    integer,dimension(:)::c,ci,cr
    if(.not.c_%stable_da) return

    do i=1,nd2
       call ctor(c(i),cr(i),ci(i))
    enddo
    return
  end subroutine ctord
  subroutine rtocd(cr,ci,c)
    implicit none
    !  INVERSE OF CTORD
    integer i
    integer,dimension(:)::c,ci,cr
    if(.not.c_%stable_da) return

    do i=1,nd2
       call rtoc(cr(i),ci(i),c(i))
    enddo
    return
  end subroutine rtocd
  subroutine resvec(cr,ci)
    implicit none
    ! DOES THE SPINOR PART IN CTORFLO
    integer i
    integer,dimension(:)::ci,cr
    integer,dimension(2)::tr,ti

    if(.not.c_%stable_da) return


    call etall(tr,2)
    call etall(ti,2)

    do i=1,nd-ndc
       if(ista(i).eq.1) then
          call dasub(cr(2*i-1),ci(2*i),tr(1))
          call daadd(ci(2*i-1),cr(2*i),ti(1))
          call daadd(cr(2*i-1),ci(2*i),tr(2))
          call dasub(ci(2*i-1),cr(2*i),ti(2))
          call dacop(tr(1),cr(2*i-1))
          call dacop(tr(2),cr(2*i))
          call dacop(ti(1),ci(2*i-1))
          call dacop(ti(2),ci(2*i))
       else
          call daadd(cr(2*i-1),cr(2*i),tr(1))
          call daadd(ci(2*i-1),ci(2*i),ti(1))
          call dasub(cr(2*i-1),cr(2*i),tr(2))
          call dasub(ci(2*i-1),ci(2*i),ti(2))
          call dacop(tr(1),cr(2*i-1))
          call dacop(tr(2),cr(2*i))
          call dacop(ti(1),ci(2*i-1))
          call dacop(ti(2),ci(2*i))
       endif
    enddo

    !    do i=nd2-ndc2+1,nd2
    !       call dacop(cr(i),dr(i))
    !       call dacop(ci(i),di(i))
    !    enddo


    call dadal(tr,2)
    call dadal(ti,2)
    return
  end subroutine resvec

  subroutine reelflo(c,ci,f,fi)
    implicit none
    ! DOES THE SPINOR PART IN RTOCFLO
    integer i
    integer,dimension(:)::c,ci,f,fi
    integer,dimension(ndim2)::e,ei
 
    if(.not.c_%stable_da) return

 

    call etall(e,nd2)
    call etall(ei,nd2)

    do i=1,nd-ndc
       call dalin(c(2*i-1),0.5_dp,c(2*i),0.5_dp,e(2*i-1))
       call dalin(ci(2*i-1),0.5_dp,ci(2*i),0.5_dp,ei(2*i-1))
       if(ista(i).eq.1) then
          call dalin(ci(2*i-1),0.5_dp,ci(2*i),-0.5_dp,e(2*i))
          call dalin(c(2*i-1),-0.5_dp,c(2*i),0.5_dp,ei(2*i))
       else
          call dalin(ci(2*i-1),0.5_dp,ci(2*i),-0.5_dp,ei(2*i))
          call dalin(c(2*i-1),0.5_dp,c(2*i),-0.5_dp,e(2*i))
       endif
    enddo

    do i=nd2-ndc2+1,nd2
       call dacop(c(i),e(i))
       call dacop(ci(i),ei(i))
    enddo

    call dacopd(e,f)
    call dacopd(ei,fi)

    call dadal(e,nd2)
    call dadal(ei,nd2)

 

    return
  end subroutine reelflo
  
  subroutine midbflo(c,a2,a2i,q,a,st)
    implicit none
    ! LINEAR EXACT NORMALIZATION USING EIGENVALUE PACKAGE OF NERI
    !
    integer i,j
    integer,dimension(ntt)::jx
    integer,dimension(:)::c,a2,a2i
    real(dp) ch,r,shm
    real(dp),dimension(ndim2,ndim2)::cr,sa,sai,cm
    real(dp),dimension(ndim)::st,q,a
    if(.not.c_%stable_da) return

    do i=1,ntt
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
          call  dapek(c(i),jx,r)
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

    call daclrd(a2)
    call daclrd(a2i)

    do i=1,nd2
       do j=1,nd2
          jx(j)=1
          r=sa(i,j)
          if(r.ne.0.0_dp)call  dapok(a2(i),jx,r)
          jx(j)=1
          r=sai(i,j)
          if(r.ne.0.0_dp)call  dapok(a2i(i),jx,r)
          jx(j)=0
       enddo
    enddo

    return
  end subroutine midbflo

  subroutine mapflol(sa,sai,cr,cm,st)
    implicit none
    !---- FROM TRACKING CODE
    ! ---------------------
    integer i,ier,iunst,j,l,n1
    integer,dimension(ndim)::n
    real(dp) ap,ax,rd,rd1,xd,xsu,xsus
    real(dp),dimension(ndim2,ndim2)::cr,xj,sa,sai,cm,w,vr,vi,s1
    real(dp),dimension(ndim)::x,xx,st
    real(dp),dimension(ndim2)::rr,ri,p
    logical hyp
    if(.not.c_%stable_da) return

    n1=0
    !     frank/etienne
    do i=1,ndim2
       do j=1,ndim2
          cr(j,i)=cm(i,j)
          xj(i,j)=0.0_dp
          s1(i,j)=0.0_dp
       enddo
    enddo

    !     frank/etienne
    do i=1,ndim
       n(i)=0
       xj(2*i-1,2*i)=1.0_dp
       xj(2*i,2*i-1)=-1.0_dp
    enddo
    !     frank/etienne
    do i=1,ndim2
       do j=1,ndim2
          sai(i,j)=0.0_dp
          w(i,j)=cm(i,j)
       enddo
    enddo
    if(ndc.eq.1) then
       s1(nd2-ndc,nd2-ndc)=1.0_dp
       s1(nd2,nd2)=1.0_dp
       sai(nd2-ndc,nd2-ndc)=1.0_dp
       sai(nd2,nd2)=1.0_dp
    endif
    call mulnd2(xj,w)
    call mulnd2(cr,w)


       xsu=0.0_dp
       do i=1,nd2


          do j=1,nd2
             xsu=xsu+abs(w(i,j)-XJ(I,J))
          enddo
       enddo
xsus=(xsu)/ND2
  if(lielib_print(6)==1) then
       !     write(w_p%c(1),'(a29,g23.16,a2)') 'Deviation from symplecticity ',c_100*(xsu)/ND2, ' %'
       write(6,'(a29,g23.16,a2)') 'Deviation from symplecticity ',100.0_dp*(xsu)/ND2, ' %'
       !CALL !WRITE_a
    endif
    call eig6(cr,rr,ri,vr,vi)
    rr_eigen=0.0_dp
    ri_eigen=0.0_dp
    rr_eigen=rr
    ri_eigen=ri
       hyp=.false.
       !      write(6,*) " checking no_hyperbolic_in_normal_form "
       do i=1,nd2-ndc2
          if(ri(i)==0.0_dp) then
             hyp=.true.
             c_%stable_da=.false.
             c_%check_stable=.false.
             messagelost="d_lielib.f90 mapflol : one of ri components is 0"
          endif
       enddo       

      
     if(hyp) then
        if(no_hyperbolic_in_normal_form) then  !no_hyperbolic_in_normal_form
         write(6,*) " Eigenvalues are "
          do i=1,nd2-ndc2
             write(6,*) i,rr(i),ri(i)
          enddo
          write(6,*) " HYPERBOLIC NORMAL FORM DETECTED "
          write(6,*) " All TPSA/DA/LIE CALCULATIONS INTERRUPTED AT YOUR REQUEST "
          write(6,*) " PLEASE RESET STABLE FLAGS "
          else
           write(6,*) " HYPERBOLIC NORMAL FORM DETECTED "
           write(6,*) " HOPE YOU KNOW WHAT YOU ARE DOING "
         endif
    endif  ! no_hyperbolic_in_normal_form
    
    !  checking for Krein    
    if(check_krein.and.(.not.hyp)) then
   
     if(.not.hyp.and.nd2>2) then
       xsu=0.0_dp
       xd=0.0_dp
 
       do i=1,4
        xsu=log(rr(i)**2+ri(i)**2)+xsu
        xd=abs(log(rr(i)**2+ri(i)**2))+xd
       enddo
 
       if(xsu>=0.and.xd>size_krein.and.xsus<=size_krein) then
         write(6,*) " A Krein collision seemed to have happened "
         write(6,*) " All calculations interrupted "
         
         do i=1,nd2-ndc
           write(6,*)"damping ", log(rr(i)**2+ri(i)**2) 
         enddo
         
         do i=1,nd2-ndc
           write(6,*)"tunes ",  atan2(ri(i),rr(i))/twopi
         enddo
         
         do i=1,nd2
           write(6,*)"eigenvalues ",  rr(i),ri(i)
         enddo
         c_%stable_da=.false.
         c_%check_stable=.false.
          
         messagelost="d_lielib.f90: mapflol: A Krein collision seemed to have happened" 
       endif
       
     endif  
endif     
           
    if(lielib_print(7)==-1) then

       do i=1,nd-ndc
          rd1=SQRT(rr(2*i-1)**2+ri(2*i-1)**2)
          rd=SQRT(rr(2*i)**2+ri(2*i)**2)
   !       write(6,*) "modulus ",rd1,rd

          write(6,'(i4,2(1x,g21.14))') 2*i-1,rr(2*i-1),ASIN(ri(2*i-1)/rd1)*twopii
          write(6,'(i4,2(1x,g21.14))') 2*i,rr(2*i),ASIN(ri(2*i)/rd)*twopii
          write(6,'(a8,g21.14)') ' alphas ', LOG(SQRT(rd*rd1))
          !CALL !WRITE_a
       enddo

       !      write(w_p%c(1),'(a8,i4,a40)') ' select ',nd-ndc,' eigenplanes (odd integers <0 real axis)'
       write(6,'(a8,i4,a40)') ' select ',nd-ndc,' eigenplanes (odd integers <0 real axis)'

       call read(n,nd-ndc)
    elseif(lielib_print(8)==-1) then
       do i=1,nd-ndc
          n(i)=nplane(i)
       enddo
       !    elseif(idpr.eq.-101.or.idpr.eq.-102) then
    else  ! new line
       do i=1,nd-ndc
          if(ri(2*i).ne.0.0_dp) then
             n(i)=2*i-1
          else
             n(i)=-2*i+1
          endif
       enddo
       !    else
       !       do i=1,nd-ndc
       !          n(i)=2*i-1
       !       enddo
    endif
    iunst=0
    do i=1,nd-ndc                  ! Frank NDC  kept
       if(n(i).lt.0) then
          n(i)=-n(i)
          st(i)=0.0_dp
          iunst=1
       else
          st(i)=1.0_dp
       endif
       x(i)=0.0_dp
       xx(i)=1.0_dp
       do j=1,nd-ndc
          x(i)=vr(2*j-1,n(i))*vi(2*j,n(i))-vr(2*j,n(i))*vi(2*j-1,n(i))+x(i)
       enddo
    enddo

    do i=1,nd-ndc
       if(x(i).lt.0.0_dp) xx(i)=-1.0_dp
       x(i)=SQRT(abs(x(i)))
       if(.not.courant_snyder) x(i)=1.0_dp
    enddo
    do i=1,nd2-ndc2
       do j=1,nd-ndc
          if(st(j)+1e-3_dp.gt.1.0_dp) then
             sai(2*j-1,i)=vr(i,n(j))*xx(j)/x(j)
             sai(2*j,i)=vi(i,n(j))/x(j)
          else
             ax=vr(i,n(j))*xx(j)/x(j)
             ap=vi(i,n(j))/x(j)
             sai(2*j-1,i)=(ax+ap)/SQRT(2.0_dp)
             sai(2*j,i)=(ap-ax)/SQRT(2.0_dp)
          endif
       enddo
    enddo
    !    if(idpr.eq.-101.or.idpr.eq.-102) then
    if(lielib_print(7)/=-1)  call movearou(sai)
    !    endif
    ! adjust sa such that sa(1,2)=0 and sa(3,4)=zero (courant-snyder-edwards-teng
    ! phase advances)
    if(iunst.ne.1) then
       do i=1,nd-ndc
          p(i)=ATAN(-sai(2*i-1,2*i)/sai(2*i,2*i))
          s1(2*i-1,2*i-1)=COS(p(i))
          s1(2*i,2*i)=COS(p(i))
          s1(2*i-1,2*i)=SIN(p(i))
          s1(2*i,2*i-1)=-SIN(p(i))
       enddo
       call mulnd2(s1,sai)
       ! adjust sa to have sa(1,1)>0 and sa(3,3)>0 rotate by pi if necessary.
       do i=1,nd-ndc
          xd=1.0_dp
          if(sai(2*i-1,2*i-1).lt.0.0_dp) xd=-1.0_dp
          s1(2*i-1,2*i-1)=xd
          s1(2*i-1,2*i)=0.0_dp
          s1(2*i,2*i-1)=0.0_dp
          s1(2*i,2*i)=xd
       enddo
       if(courant_snyder) call mulnd2(s1,sai)
       ! sa is now uniquely and unambigeously determined.
    endif
    !  do i=1,nd2
    !     do l=1,nd2
    !        sa(i,l)=sai(i,l)
    !     enddo
    !  enddo
    call matinv(sai,sa,nd2,ndim2,ier)

    call mulnd2(sai,cm)
    do i=1,nd2
       do j=1,nd2
          cr(i,j)=sa(i,j)
       enddo
    enddo

    call mulnd2(cm,cr)

    !!!  new in case ac modulation: correcting the time 
    if(ndpt/=0) then
     l=ndpt+1
     if(mod(ndpt,2)==0) l=l-2
!     write(6,*) "ndpt,l ",ndpt,l
     vr=0.0_dp
     vi=0.0_dp
     do i=1,nd2
      vr(i,i)=1.0_dp
      vi(i,i)=1.0_dp
     enddo
      ax= cr(l,nd2-3)
      ap= cr(l,nd2-2)
      xsu = cr(nd2-3,nd2-2)/(1.0_dp-cr(nd2-3,nd2-3))
!      write(6,*) ax,ap,xsu
      vr(l,nd2-3)=-0.5_dp*(ax-xsu*ap)
      vr(l,nd2-2)=-0.5_dp*(ap+xsu*ax)
      vi(l,nd2-3)=-vr(l,nd2-3)
      vi(l,nd2-2)=-vr(l,nd2-2)
       
      sa(1:nd2,1:nd2)=matmul(sa(1:nd2,1:nd2),vr(1:nd2,1:nd2))
      sai(1:nd2,1:nd2)=matmul(vi(1:nd2,1:nd2),sai(1:nd2,1:nd2))

!     do i=1,nd2
!           write(6,'(8(1x,e14.7))') vr(i,1:nd2) 
!     enddo
!     do i=1,nd2
!           write(6,'(8(1x,e14.7))') cr(i,1:nd2) 
!     enddo
     call mulnd2(cr,vr)
     call mulnd2(vi,vr)
     cr=vr
 
!     do i=1,nd2
!           write(6,'(8(1x,e14.7))') cr(i,1:nd2) 
!     enddo
    endif
 
    return
  end subroutine mapflol

  subroutine mulnd2(rt,r)
    implicit none
    integer i,ia,j
    real(dp),dimension(ndim2,ndim2)::rt,r,rtt
    if(.not.c_%stable_da) return

    do i=1,nd2
       do j=1,nd2
          rtt(i,j)=0.0_dp
       enddo
    enddo
    do i=1,nd2
       do j=1,nd2
          do ia=1,nd2
             rtt(i,ia)=rt(i,j)*r(j,ia)+rtt(i,ia)
          enddo
       enddo
    enddo

    do i=1,nd2
       do j=1,nd2
          r(i,j)=rtt(i,j)
       enddo
    enddo
    return
  end subroutine mulnd2

  subroutine movearou(rt)
    implicit none
    !      integer ipause, mypause
    integer i,ic,j
    real(dp) xr,xrold
    real(dp),dimension(ndim2,ndim2)::rt,rto,s
    real(dp),dimension(ndim2,ndim2):: xt,yt,zt,xy,xz,yz
    real(dp),dimension(ndim2,ndim2):: xyz,xzy,xyt,yxt,yzt,zyt,xzt,zxt
    real(dp),dimension(ndim2,ndim2):: xyzt,xytz,xzyt,xzty,xtzy,xtyz

    if(.not.c_%stable_da) return

    do i=1,nd2
       do j=1,nd2
          s(i,j)=0.0_dp
          s(i,i)=1.0_dp
       enddo
    enddo
    xt=0.0_dp;yt=0.0_dp;zt=0.0_dp;xy=0.0_dp;xz=0.0_dp;yz=0.0_dp;
    xyzt=0.0_dp;xytz=0.0_dp;xzyt=0.0_dp;xzty=0.0_dp;xtzy=0.0_dp;xtyz=0.0_dp;
    xyz=0.0_dp;xzy=0.0_dp;xyt=0.0_dp;yxt=0.0_dp;yzt=0.0_dp;zyt=0.0_dp;xzt=0.0_dp;zxt=0.0_dp;

    do i=0,1

       xy(1+i,3+i)=1.0_dp
       xy(3+i,1+i)=1.0_dp
       xy(5+i,5+i)=1.0_dp
       xy(7+i,7+i)=1.0_dp

       xz(1+i,5+i)=1.0_dp
       xz(5+i,1+i)=1.0_dp
       xz(3+i,3+i)=1.0_dp
       xz(7+i,7+i)=1.0_dp

       xt(1+i,7+i)=1.0_dp
       xt(7+i,1+i)=1.0_dp
       xt(3+i,3+i)=1.0_dp
       xt(5+i,5+i)=1.0_dp

       yz(3+i,5+i)=1.0_dp
       yz(5+i,3+i)=1.0_dp
       yz(1+i,1+i)=1.0_dp
       yz(7+i,7+i)=1.0_dp

       yt(3+i,7+i)=1.0_dp
       yt(7+i,3+i)=1.0_dp
       yt(1+i,1+i)=1.0_dp
       yt(5+i,5+i)=1.0_dp

       zt(5+i,7+i)=1.0_dp
       zt(7+i,5+i)=1.0_dp
       zt(1+i,1+i)=1.0_dp
       zt(3+i,3+i)=1.0_dp

       xyz(1+i,3+i)=1.0_dp
       xyz(3+i,5+i)=1.0_dp
       xyz(5+i,1+i)=1.0_dp
       xyz(7+i,7+i)=1.0_dp

       xyz(1+i,3+i)=1.0_dp
       xyz(3+i,5+i)=1.0_dp
       xyz(5+i,1+i)=1.0_dp
       xyz(7+i,7+i)=1.0_dp

       xzy(1+i,5+i)=1.0_dp
       xzy(5+i,3+i)=1.0_dp
       xzy(3+i,1+i)=1.0_dp
       xzy(7+i,7+i)=1.0_dp

       xyt(1+i,3+i)=1.0_dp
       xyt(3+i,7+i)=1.0_dp
       xyt(7+i,1+i)=1.0_dp
       xyt(5+i,5+i)=1.0_dp

       yxt(3+i,1+i)=1.0_dp
       yxt(1+i,7+i)=1.0_dp
       yxt(7+i,3+i)=1.0_dp
       yxt(5+i,5+i)=1.0_dp

       yzt(3+i,5+i)=1.0_dp
       yzt(5+i,7+i)=1.0_dp
       yzt(7+i,3+i)=1.0_dp
       yzt(1+i,1+i)=1.0_dp

       zyt(5+i,3+i)=1.0_dp
       zyt(3+i,7+i)=1.0_dp
       zyt(7+i,5+i)=1.0_dp
       zyt(1+i,1+i)=1.0_dp

       xzt(1+i,5+i)=1.0_dp
       xzt(5+i,7+i)=1.0_dp
       xzt(7+i,1+i)=1.0_dp
       xzt(3+i,3+i)=1.0_dp

       zxt(5+i,1+i)=1.0_dp
       zxt(1+i,7+i)=1.0_dp
       zxt(7+i,5+i)=1.0_dp
       zxt(3+i,3+i)=1.0_dp

       xyzt(1+i,3+i)=1.0_dp
       xyzt(3+i,5+i)=1.0_dp
       xyzt(5+i,7+i)=1.0_dp
       xyzt(7+i,1+i)=1.0_dp

       xytz(1+i,3+i)=1.0_dp
       xytz(3+i,7+i)=1.0_dp
       xytz(7+i,5+i)=1.0_dp
       xytz(5+i,1+i)=1.0_dp

       xzyt(1+i,5+i)=1.0_dp
       xzyt(5+i,3+i)=1.0_dp
       xzyt(3+i,7+i)=1.0_dp
       xzyt(7+i,1+i)=1.0_dp

       xzty(1+i,5+i)=1.0_dp
       xzty(5+i,7+i)=1.0_dp
       xzty(7+i,3+i)=1.0_dp
       xzty(3+i,1+i)=1.0_dp

       xtzy(1+i,7+i)=1.0_dp
       xtzy(7+i,5+i)=1.0_dp
       xtzy(5+i,3+i)=1.0_dp
       xtzy(3+i,1+i)=1.0_dp

       xtyz(1+i,7+i)=1.0_dp
       xtyz(7+i,3+i)=1.0_dp
       xtyz(3+i,5+i)=1.0_dp
       xtyz(5+i,1+i)=1.0_dp
    enddo

    !do i=1,8
    !write(6,100) (rt(i,j),j=1,8)
    !enddo
    !write(6,*) " "
    !100  FORMAT(8(1x, E12.6))

    ic=0
    xrold=1e9_dp
    call movemul(rt,s,rto,xr)

    if(xr.lt.xrold) then
       xrold=xr
    endif

    if(nd>=2) then
       call movemul(rt,xy,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=1
       endif
    endif

    if(nd>=3) then
       call movemul(rt,xz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=2
       endif
       call movemul(rt,yz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=3
       endif
       call movemul(rt,xyz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=4
       endif
       call movemul(rt,xzy,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=5
       endif
    endif

    if(nd.eq.4) then
       call movemul(rt,xt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=6
       endif
       call movemul(rt,yt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=7
       endif
       call movemul(rt,zt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=8
       endif

       call movemul(rt,xyt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=9
       endif
       call movemul(rt,yxt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=10
       endif
       call movemul(rt,yzt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=11

       endif
       call movemul(rt,zyt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=12
       endif
       call movemul(rt,xzt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=13
       endif
       call movemul(rt,zxt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=14
       endif
       call movemul(rt,xyzt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=15
       endif
       call movemul(rt,xytz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=16
       endif
       call movemul(rt,xzyt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=17
       endif
       call movemul(rt,xzty,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=18
       endif
       call movemul(rt,xtzy,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=19
       endif
       call movemul(rt,xtyz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=20
       endif
    endif


    i=0
    if(ic.eq.0) then
       call movemul(rt,s,rto,xr)
       i=i+1
!       w_p%c(i)=  " no exchanged"
    elseif(ic.eq.1) then
       call movemul(rt,xy,rto,xr)
       i=i+1
!       w_p%c(i)=  " x-y exchanged"
    elseif(ic.eq.2) then
       call movemul(rt,xz,rto,xr)
       i=i+1
!       w_p%c(i)=  " x-z exchanged"
    elseif(ic.eq.3) then
       call movemul(rt,yz,rto,xr)
       i=i+1
!       w_p%c(i)= " y-z exchanged"
    elseif(ic.eq.4) then
       call movemul(rt,xyz,rto,xr)
       i=i+1
 !      w_p%c(i)=  " x-y-z permuted"
    elseif(ic.eq.5) then
       call movemul(rt,xzy,rto,xr)
       i=i+1
 !      w_p%c(i)=  " x-z-y permuted"
    elseif(ic.eq.6) then
       call movemul(rt,xt,rto,xr)
    elseif(ic.eq.7) then
       call movemul(rt,yt,rto,xr)
    elseif(ic.eq.8) then
       call movemul(rt,zt,rto,xr)
    elseif(ic.eq.9) then
       call movemul(rt,xyt,rto,xr)
    elseif(ic.eq.10) then
       call movemul(rt,yxt,rto,xr)
    elseif(ic.eq.11) then
       call movemul(rt,yzt,rto,xr)
    elseif(ic.eq.12) then
       call movemul(rt,zyt,rto,xr)
    elseif(ic.eq.13) then
       call movemul(rt,xzt,rto,xr)
    elseif(ic.eq.14) then
       call movemul(rt,zxt,rto,xr)
    elseif(ic.eq.15) then
       call movemul(rt,xyzt,rto,xr)
    elseif(ic.eq.16) then
       call movemul(rt,xytz,rto,xr)
    elseif(ic.eq.17) then
       call movemul(rt,xzyt,rto,xr)
    elseif(ic.eq.18) then
       call movemul(rt,xzty,rto,xr)
    elseif(ic.eq.19) then
       call movemul(rt,xtzy,rto,xr)
    elseif(ic.eq.20) then
       call movemul(rt,xtyz,rto,xr)
    endif





    do i=1,nd2
       do j=1,nd2
          rt(i,j)=rto(i,j)
       enddo
    enddo

    return
  end subroutine movearou

  subroutine movemul(rt,xy,rto,xr)
    implicit none
    integer i,j,k
    real(dp) xr
    real(dp),dimension(ndim2,ndim2)::rt,xy,rto
    if(.not.c_%stable_da) return

    do i=1,nd2
       do j=1,nd2
          rto(i,j)=0.0_dp
       enddo
    enddo

    do  i=1,nd2
       do  j=1,nd2
          do  k=1,nd2
             rto(i,k)=xy(i,j)*rt(j,k)+rto(i,k)
          enddo
       enddo
    enddo

    xr=0.0_dp
    do i=1,nd2
       do j=1,nd2
          xr=xr+abs(rto(i,j))
       enddo
    enddo
    do i=1,nd
       xr=xr-abs(rto(2*i-1,2*i-1))
       xr=xr-abs(rto(2*i-1,2*i))
       xr=xr-abs(rto(2*i,2*i))
       xr=xr-abs(rto(2*i,2*i-1))
    enddo
    return
  end subroutine movemul

  subroutine initpert(st,ang,ra)
    implicit none
    !   X-RATED
    !- SETS UP ALL THE !COMMON BLOCKS RELEVANT TO NORMAL FORM AND THE BASIS
    !- CHANGES INSIDE  MAPNORMF
    integer i,nn,ipause,mypauses
    real(dp),dimension(ndim)::ang,ra,st
    if(.not.c_%stable_da) return

    if(iref.gt.0) then
 
       read(iref,*) nres
       if(nres.ge.nreso) then
          line= ' NRESO IN LIELIB TOO SMALL '
          ipause=mypauses(-999,line)
       endif
    elseif(iref.eq.0) then
       nres=0
    endif
    if(nres.ne.0.and.global_verbose) then

    endif
    if(iref.gt.0) then
       do i=1,nres
          read(iref,*) (mx(nn,i),nn=1,nd-ndc)
       enddo
    endif
    do i=nres+1,nreso
       do nn=1,ndim
          mx(nn,i)=0
       enddo
    enddo
    !      frank/Etienne
    do i=1,ndim
       angle(i)=0.0_dp
       rad(i)=0.0_dp
       sta(i)=0.0_dp
       dsta(i)=1.0_dp-sta(i)
       ista(i)=0
       idsta(i)=0
    enddo
    do i=1,nd        !  Frank          -ndc
       angle(i)=ang(i)
       rad(i)=ra(i)
       sta(i)=st(i)
       dsta(i)=1.0_dp-sta(i)
    enddo
    do i=1,nd
       ista(i)=int(sta(i)+1e-2_dp)
       idsta(i)=int(dsta(i)+1e-2_dp)
    enddo
    return
  end subroutine initpert

  real(dp) function dlie(j)
    implicit none
    integer i
    !      INTEGER J(NTT)
    integer,dimension(:)::j
    if(.not.c_%stable_da) return

    dlie=0.0_dp
    do i=1,nd
       dlie=REAL(j(2*i-1)+j(2*i),kind=DP)+dlie
    enddo
    dlie=dlie+1.0_dp
    dlie=1.0_dp/dlie
    return
  end function dlie

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
    !frs      goto(11,12,13,14),mo
    ! 11   rext = one
    !      return
    ! 12   rext = -one
    !      return
    ! 13   rext = -one
    !      return
    ! 14   rext = one
    !      return
    !frs
  end function rext
  subroutine cpart(h,ch)  ! no flip needed
    implicit none
    integer h,ch
    if(.not.c_%stable_da) return

    call dacfu(h,rext,ch)
    return
  end subroutine cpart

  subroutine ctoi(f1,f2) ! no flip needed
    implicit none
    integer f1,f2,b1
    integer,dimension(ndim2)::x
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etallnom(x,nd2) !  ,'X         ')

    call cpart(f1,b1)
    call etctr(x)
    call trx(b1,f2,x)
    call dadal(x,nd2)
    call dadal1(b1)
    return
  end subroutine ctoi

  subroutine itoc(f1,f2) ! no flip needed
    implicit none
    integer f1,f2,b1
    integer,dimension(ndim2)::x
    if(.not.c_%stable_da) return

    call etall1(b1)
    call etallnom(x,nd2) !  ,'X         ')

    call etrtc(x)
    call trx(f1,b1,x)
    call cpart(b1,f2)
    call dadal(x,nd2)
    call dadal1(b1)
    return
  end subroutine itoc

  subroutine etrtc(x) ! no flip needed
    implicit none
    integer i
    integer,dimension(:)::x
    integer,dimension(ndim2)::rel
    if(.not.c_%stable_da) return

    call etallnom(rel,nd2) !  ,'REL       ')

    call etini(rel)
    call etini(x)
    do i=1,nd-ndc
       call daadd(rel(2*i-1),rel(2*i),x(2*i-1))
       call dasub(rel(2*i-1),rel(2*i),x(2*i))
    enddo
    call dadal(rel,nd2)
    return
  end subroutine etrtc

  subroutine etctr(x) ! no flip needed
    implicit none
    integer i
    integer,dimension(:)::x
    integer,dimension(ndim2)::rel
    if(.not.c_%stable_da) return

    call etallnom(rel,nd2) !  ,'REL       ')

    call etini(rel)
    call etini(x)
    do i=1,nd-ndc
       call dalin(rel(2*i-1),0.5_dp,rel(2*i),0.5_dp,x(2*i-1))
       call dalin(rel(2*i-1),0.5_dp,rel(2*i),-0.5_dp,x(2*i))
    enddo
    call dadal(rel,nd2)
    return
  end subroutine etctr

  subroutine etcjg(x) ! no flip needed
    implicit none
    integer i
    integer,dimension(:)::x
    integer,dimension(ndim2)::rel
    if(.not.c_%stable_da) return

    call etallnom(rel,nd2) !  ,'REL       ')

    call etini(rel)
    call etini(x)
    do i=1,nd-ndc
       if(ista(i).eq.1) then
          call dacop(rel(2*i-1),x(2*i))
          call dacop(rel(2*i),x(2*i-1))
       else
          call dacop(rel(2*i-1),x(2*i-1))
          call dacop(rel(2*i),x(2*i))
       endif
    enddo
    call dadal(rel,nd2)
    return
  end subroutine etcjg

  ! Neri's Routine below

  subroutine eig6(fm,reval,aieval,revec,aievec)
    implicit none
    !**************************************************************************

    !  Diagonalization routines of NERI

    !ccccccccccccccccc
    !
    !  this routine finds the eigenvalues and eigenvectors
    !  of the full matrix fm.
    !  the eigenvectors are normalized so that the real and
    !  imaginary part of vectors 1, 3, and 5 have +1 antisymmetric
    !  product:
    !      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1 ;
    !      revec5 J aivec5 = one
    !  the eigenvectors 2 ,4, and 6 have the opposite normalization.
    !  written by F. Neri, Feb 26 1986.
    !
    integer jet,nn,i,i1,ilo,ihi,mdim,info
    real(dp),dimension(ndim2)::reval,aieval,ort
    real(dp),dimension(ndim2,ndim2)::revec,aievec,fm,aa,vv
    INTEGER IPAUSE,MYPAUSES
    if(.not.c_%stable_da) return

    !  copy matrix to temporary storage (the matrix aa is destroyed)
    do i=1,nd2-ndc2
       do i1=1,nd2-ndc2
          aa(i1,i) = fm(i1,i)
       enddo
    enddo
    ilo = 1
    ihi = nd2-ndc2
    mdim = ndim2
    nn = nd2-ndc2
    !  compute eigenvalues and eigenvectors using double
    !  precision Eispack routines:
    call ety(mdim,nn,ilo,ihi,aa,ort)
    call etyt(mdim,nn,ilo,ihi,aa,ort,vv)
    call ety2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)
    if ( info .ne. 0 ) then
       LINE= '  ERROR IN EIG6'
       IPAUSE=MYPAUSES(0,LINE)
    endif
    !      call neigv(vv,pbkt)
    do i=1,nd-ndc
       do jet=1,nd2-ndc2
          revec(jet,2*i-1)=vv(jet,2*i-1)
          revec(jet,2*i)=vv(jet,2*i-1)
          aievec(jet,2*i-1)=vv(jet,2*i)
          aievec(jet,2*i)=-vv(jet,2*i)
       enddo
    enddo
    do i=1,nd2-ndc2
       if(abs(reval(i)**2+aieval(i)**2 -1.0_dp).gt.1e-10_dp) then
 
          if(lielib_print(4)==1) then
 
             write(6,*) sqrt(reval(i)**2+aieval(i)**2)
          endif
       endif
    enddo
    return
  end subroutine eig6

  subroutine ety(nm,n,low,igh,a,ort)
    implicit none
    !
    !     this subroutine is a translation of the algol procedure orthes,
    !     num. math. 12, 349-368(1968) by martin and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
    !
    !     given a real general matrix, this subroutine
    !     reduces a submatrix situated in rows and columns
    !     low through igh to upper hessenberg form by
    !     orthogonal similarity transformations.
    !
    !     on input-
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement,
    !
    !        n is the order of the matrix,
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been used,
    !          set low=1, igh=n,
    !
    !        a contains the input matrix.
    !
    !     on output-
    !
    !        a contains the hessenberg matrix.  information about
    !          the orthogonal transformations used in the reduction
    !          is stored in the remaining triangle under the
    !          hessenberg matrix,
    !
    !        ort contains further information about the transformations.
    !          only elements low through igh are used.
    !
    !     fortran routine by b. s. garbow
    !     modified by filippo neri.
    !
    !
    integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
    real(dp),dimension(nm,n)::a
    real(dp),dimension(igh)::ort
    real(dp) f,g,h,scale
    if(.not.c_%stable_da) return

    la = igh - 1
    kp1 = low + 1
    if (la .lt. kp1) go to 200
    !
    do m = kp1, la
       h = 0.0_dp
       ort(m) = 0.0_dp
       scale = 0.0_dp
       !     ********** scale column (algol tol then not needed) **********
       do i = m, igh
          scale = scale + abs(a(i,m-1))
       enddo
       !
       if (scale .eq. 0.0_dp) go to 180
       mp = m + igh
       !     ********** for i=igh step -1 until m do -- **********
       do ii = m, igh
          i = mp - ii
          ort(i) = a(i,m-1) / scale
          h = h + ort(i) * ort(i)
       enddo
       !
       g = -sign(SQRT(h),ort(m))
       h = h - ort(m) * g
       ort(m) = ort(m) - g
       !     ********** form (i-(u*ut)/h) * a **********
       do j = m, n
          f = 0.0_dp
          !     ********** for i=igh step -1 until m do -- **********
          do ii = m, igh
             i = mp - ii
             f = f + ort(i) * a(i,j)
          enddo
          !
          f = f / h
          !
          do i = m, igh
             a(i,j) = a(i,j) - f * ort(i)
          enddo
          !
       enddo
       !     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
       do i = 1, igh
          f = 0.0_dp
          !     ********** for j=igh step -1 until m do -- **********
          do jj = m, igh
             j = mp - jj
             f = f + ort(j) * a(i,j)
          enddo
          !
          f = f / h
          !
          do j = m, igh
             a(i,j) = a(i,j) - f * ort(j)
          enddo
          !
       enddo
       !
       ort(m) = scale * ort(m)
       a(m,m-1) = scale * g
180    continue
    enddo
    !
200 return
    !     ********** last card of ety **********
  end subroutine ety
  subroutine etyt(nm,n,low,igh,a,ort,z)
    implicit none
    !
    !     this subroutine is a translation of the algol procedure ortrans,
    !     num. math. 16, 181-204(1970) by peters and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
    !
    !     this subroutine accumulates the orthogonal similarity
    !     transformations used in the reduction of a real general
    !     matrix to upper hessenberg form by  ety.
    !
    !     on input-
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement,
    !
    !        n is the order of the matrix,
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been used,
    !          set low=1, igh=n,
    !
    !        a contains information about the orthogonal trans-
    !          formations used in the reduction by  orthes
    !          in its strict lower triangle,
    !
    !          ort contains further information about the trans-
    !          formations used in the reduction by  ety.
    !          only elements low through igh are used.
    !
    !     on output-
    !
    !        z contains the transformation matrix produced in the
    !          reduction by  ety,
    !
    !        ort has been altered.
    !
    !     fortran routine by b. s. garbow.
    !     modified by f. neri.
    !
    !
    integer i,j,n,kl,mm,mp,nm,igh,low,mp1
    real(dp) g
    real(dp),dimension(igh)::ort
    real(dp),dimension(nm,igh)::a
    real(dp),dimension(nm,n)::z
    if(.not.c_%stable_da) return

    !     ********** initialize z to identity matrix **********
    do i = 1, n
       !
       do j = 1, n
          z(i,j) = 0.0_dp
       enddo
       !
       z(i,i) = 1.0_dp
    enddo
    !
    kl = igh - low - 1
    if (kl .lt. 1) go to 200
    !     ********** for mp=igh-1 step -1 until low+1 do -- **********
    do mm = 1, kl
       mp = igh - mm
       if (a(mp,mp-1) .eq. 0.0_dp) go to 140
       mp1 = mp + 1
       !
       do i = mp1, igh
          ort(i) = a(i,mp-1)
       enddo
       !
       do j = mp, igh
          g = 0.0_dp
          !
          do i = mp, igh
             g = g + ort(i) * z(i,j)
          enddo
          !     ********** divisor below is negative of h formed in orthes.
          !                double division avoids possible underflow **********
          g = (g / ort(mp)) / a(mp,mp-1)
          !
          do i = mp, igh
             z(i,j) = z(i,j) + g * ort(i)
          enddo
          !
       enddo
       !
140    continue
    enddo
    !
200 return
    !     ********** last card of etyt **********
  end subroutine etyt
  subroutine ety2(nm,n,low,igh,h,wr,wi,z,ierr)
    implicit none
    !
    !
    !
    !     this subroutine is a translation of the algol procedure hqr2,
    !     num. math. 16, 181-204(1970) by peters and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
    !
    !     this subroutine finds the eigenvalues and eigenvectors
    !     of a real upper hessenberg matrix by the qr method.  the
    !     eigenvectors of a real general matrix can also be found
    !     if  elmhes  and  eltran  or  orthes  and  ortran  have
    !     been used to reduce this general matrix to hessenberg form
    !     and to accumulate the similarity transformations.
    !
    !     on input-
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement,
    !
    !        n is the order of the matrix,
    !
    !        low and igh are integers determined by the balancing
    !          subroutine  balanc.  if  balanc  has not been used,
    !          set low=1, igh=n,
    !
    !        h contains the upper hessenberg matrix,
    !
    !        z contains the transformation matrix produced by  eltran
    !          after the reduction by  elmhes, or by  ortran  after the
    !          reduction by  orthes, if performed.  if the eigenvectors
    !          of the hessenberg matrix are desired, z must contain the
    !          identity matrix.
    !
    !     on output-
    !
    !        h has been destroyed,
    !
    !        wr and wi contain the real and imaginary parts,
    !          respectively, of the eigenvalues.  the eigenvalues
    !          are unordered except that complex conjugate pairs
    !          of values appear consecutively with the eigenvalue
    !          having the positive imaginary part first.  if an
    !          error exit is made, the eigenvalues should be correct
    !          for indices ierr+1,...,n,
    !
    !        z contains the real and imaginary parts of the eigenvectors.
    !          if the i-th eigenvalue is real, the i-th column of z
    !          contains its eigenvector.  if the i-th eigenvalue is complex
    !          with positive imaginary part, the i-th and (i+1)-th
    !          columns of z contain the real and imaginary parts of its
    !          eigenvector.  the eigenvectors are unnormalized.  if an
    !          error exit is made, none of the eigenvectors has been found,
    !
    !        ierr is set to
    !          zero       for normal return,
    !          j          if the j-th eigenvalue has not been
    !                     determined after 200 iterations.
    !
    !     arithmetic is real(dp). complex division
    !     is simulated by routin etdiv.
    !
    !     fortran routine by b. s. garbow.
    !     modified by f. neri.
    !
    !
    logical(lp) notlas
    integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,igh,its,low,mp2,enm2,ierr
    real(dp) p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,z3r,z3i
    real(dp),dimension(n)::wr,wi
    real(dp),dimension(nm,n)::h,z
    if(.not.c_%stable_da) return

    !     ********** machep is a machine dependent parameter specifying
    !                the relative precision of floating point arithmetic.
    !
    !                **********
    !     machep = r1mach(4)
    !
    ierr = 0
    norm = 0.0_dp
    k = 1
    !     ********** store roots isolated by balanc
    !                and compute matrix norm **********
    do i = 1, n
       !
       do j = k, n
          norm = norm + abs(h(i,j))
       enddo
       !
       k = i
       if (i .ge. low .and. i .le. igh) go to 50
       wr(i) = h(i,i)
       wi(i) = 0.0_dp
50     continue
    enddo
    !
    en = igh
    t = 0.0_dp
    !     ********** search for next eigenvalues **********
60  if (en .lt. low) go to 340
    its = 0
    na = en - 1
    enm2 = na - 1
    !     ********** look for single small sub-diagonal element
    !                for l=en step -1 until low do -- **********
70  do ll = low, en
       l = en + low - ll
       if (l .eq. low) go to 100
       s = abs(h(l-1,l-1)) + abs(h(l,l))
       if (s .eq. 0.0_dp) s = norm
       if (abs(h(l,l-1)) .le. machep * s) go to 100
    enddo
    !     ********** form shift **********
100 x = h(en,en)
    if (l .eq. en) go to 270
    y = h(na,na)
    w = h(en,na) * h(na,en)
    if (l .eq. na) go to 280
    if (its .eq. 200) go to 1000
    if (its .ne. 10 .and. its .ne. 20) go to 130
    !     ********** form exceptional shift **********
    t = t + x
    !
    do i = low, en
       h(i,i) = h(i,i) - x
    enddo
    !
    s = abs(h(en,na)) + abs(h(na,enm2))
    x = 0.75_dp * s
    y = x
    w = -0.4375_dp * s * s
130 its = its + 1
    !     ********** look for two consecutive small
    !                sub-diagonal elements.
    !                for m=en-2 step -1 until l do -- **********
    do mm = l, enm2
       m = enm2 + l - mm
       zz = h(m,m)
       r = x - zz
       s = y - zz
       p = (r * s - w) / h(m+1,m) + h(m,m+1)
       q = h(m+1,m+1) - zz - r - s
       r = h(m+2,m+1)
       s = abs(p) + abs(q) + abs(r)
       p = p / s
       q = q / s
       r = r / s
       if (m .eq. l) go to 150
       if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. machep * abs(p) * (abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))) go to 150
    enddo
    !
150 mp2 = m + 2
    !
    do i = mp2, en
       h(i,i-2) = 0.0_dp
       if (i .eq. mp2) go to 160
       h(i,i-3) = 0.0_dp
160    continue
    enddo
    !     ********** double qr step involving rows l to en and
    !                columns m to en **********
    do k = m, na
       notlas = k .ne. na
       if (k .eq. m) go to 170
       p = h(k,k-1)
       q = h(k+1,k-1)
       r = 0.0_dp
       if (notlas) r = h(k+2,k-1)
       x = abs(p) + abs(q) + abs(r)
       if (x .eq. 0.0_dp) go to 260
       p = p / x
       q = q / x
       r = r / x
170    s = sign(SQRT(p*p+q*q+r*r),p)
       if (k .eq. m) go to 180
       h(k,k-1) = -s * x
       go to 190
180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
190    p = p + s
       x = p / s
       y = q / s
       zz = r / s
       q = q / p
       r = r / p
       !     ********** row modification **********
       do j = k, n
          p = h(k,j) + q * h(k+1,j)
          if (.not. notlas) go to 200
          p = p + r * h(k+2,j)
          h(k+2,j) = h(k+2,j) - p * zz
200       h(k+1,j) = h(k+1,j) - p * y
          h(k,j) = h(k,j) - p * x
       enddo
       !
       j = min0(en,k+3)
       !     ********** column modification **********
       do i = 1, j
          p = x * h(i,k) + y * h(i,k+1)
          if (.not. notlas) go to 220
          p = p + zz * h(i,k+2)
          h(i,k+2) = h(i,k+2) - p * r
220       h(i,k+1) = h(i,k+1) - p * q
          h(i,k) = h(i,k) - p
       enddo
       !     ********** accumulate transformations **********
       do i = low, igh
          p = x * z(i,k) + y * z(i,k+1)
          if (.not. notlas) go to 240
          p = p + zz * z(i,k+2)
          z(i,k+2) = z(i,k+2) - p * r
240       z(i,k+1) = z(i,k+1) - p * q
          z(i,k) = z(i,k) - p
       enddo
       !
260    continue
    enddo
    !
    go to 70
    !     ********** one root found **********
270 h(en,en) = x + t
    wr(en) = h(en,en)
    wi(en) = 0.0_dp
    en = na
    go to 60
    !     ********** two roots found **********
280 p = (y - x) / 2.0_dp
    q = p * p + w
    zz = SQRT(abs(q))
    h(en,en) = x + t
    x = h(en,en)
    h(na,na) = y + t
    if (q .lt. 0.0_dp) go to 320
    !     ********** real pair **********
    zz = p + sign(zz,p)
    wr(na) = x + zz
    wr(en) = wr(na)
    if (zz .ne. 0.0_dp) wr(en) = x - w / zz
    wi(na) = 0.0_dp
    wi(en) = 0.0_dp
    x = h(en,na)
    s = abs(x) + abs(zz)
    p = x / s
    q = zz / s
    r = SQRT(p*p+q*q)
    p = p / r
    q = q / r
    !     ********** row modification **********
    do j = na, n
       zz = h(na,j)
       h(na,j) = q * zz + p * h(en,j)
       h(en,j) = q * h(en,j) - p * zz
    enddo
    !     ********** column modification **********
    do i = 1, en
       zz = h(i,na)
       h(i,na) = q * zz + p * h(i,en)
       h(i,en) = q * h(i,en) - p * zz
    enddo
    !     ********** accumulate transformations **********
    do i = low, igh
       zz = z(i,na)
       z(i,na) = q * zz + p * z(i,en)
       z(i,en) = q * z(i,en) - p * zz
    enddo
    !
    go to 330
    !     ********** complex pair **********
320 wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz
330 en = enm2
    go to 60
    !     ********** all roots found.  backsubstitute to find
    !                vectors of upper triangular form **********
340 if (norm .eq. 0.0_dp) go to 1001
    !     ********** for en=n step -1 until 1 do -- **********
    do nn = 1, n
       en = n + 1 - nn
       p = wr(en)
       q = wi(en)
       na = en - 1
       if (q.lt.0) goto 710
       if (q.eq.0) goto 600
       if (q.gt.0) goto 800
       !     ********** real vector **********
600    m = en
       h(en,en) = 1.0_dp
       if (na .eq. 0) go to 800
       !     ********** for i=en-1 step -1 until 1 do -- **********
       do ii = 1, na
          i = en - ii
          w = h(i,i) - p
          r = h(i,en)
          if (m .gt. na) go to 620
          !
          do j = m, na
             r = r + h(i,j) * h(j,en)
          enddo
          !
620       if (wi(i) .ge. 0.0_dp) go to 630
          zz = w
          s = r
          go to 700
630       m = i
          if (wi(i) .ne. 0.0_dp) go to 640
          t = w
          if (w .eq. 0.0_dp) t = machep * norm
          h(i,en) = -r / t
          go to 700
          !     ********** solve real equations **********
640       x = h(i,i+1)
          y = h(i+1,i)
          q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
          t = (x * s - zz * r) / q
          h(i,en) = t
          if (abs(x) .le. abs(zz)) go to 650
          h(i+1,en) = (-r - w * t) / x
          go to 700
650       h(i+1,en) = (-s - y * t) / zz
700       continue
       enddo
       !     ********** end real vector **********
       go to 800
       !     ********** complex vector **********
710    m = na
       !     ********** last vector component chosen imaginary so that
       !                eigenvector matrix is triangular **********
       if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
       h(na,na) = q / h(en,na)
       h(na,en) = -(h(en,en) - p) / h(en,na)
       go to 730
       ! 720    z3 = cmplx(zero,-h(na,en)) / cmplx(h(na,na)-p,q)
       !        h(na,na) = real(z3,kind=dp)
       !        h(na,en) = aimag(z3)
720    call etdiv(z3r,z3i,0.0_dp,-h(na,en),h(na,na)-p,q)
       h(na,na) = z3r
       h(na,en) = z3i
730    h(en,na) = 0.0_dp
       h(en,en) = 1.0_dp
       enm2 = na - 1
       if (enm2 .eq. 0) go to 800
       !     ********** for i=en-2 step -1 until 1 do -- **********
       do ii = 1, enm2
          i = na - ii
          w = h(i,i) - p
          ra = 0.0_dp
          sa = h(i,en)
          !
          do j = m, na
             ra = ra + h(i,j) * h(j,na)
             sa = sa + h(i,j) * h(j,en)
          enddo
          !
          if (wi(i) .ge. 0.0_dp) go to 770
          zz = w
          r = ra
          s = sa
          go to 790
770       m = i
          if (wi(i) .ne. 0.0_dp) go to 780
          !           z3 = cmplx(-ra,-sa) / cmplx(w,q)
          !           h(i,na) = real(z3,kind=dp)
          !           h(i,en) = aimag(z3)
          call etdiv(z3r,z3i,-ra,-sa,w,q)
          h(i,na) = z3r
          h(i,en) = z3i
          go to 790
          !     ********** solve complex equations **********
780       x = h(i,i+1)
          y = h(i+1,i)
          vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
          vi = (wr(i) - p) * 2.0_dp * q
          if (vr .eq. 0.0_dp .and. vi .eq. 0.0_dp) vr = machep * norm  * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
          !           z3 = cmplx(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra) / cmplx(vr,vi)
          !           h(i,na) = real(z3,kind=dp)
          !           h(i,en) = aimag(z3)
          call etdiv(z3r,z3i,x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi)
          h(i,na) = z3r
          h(i,en) = z3i
          if (abs(x) .le. abs(zz) + abs(q)) go to 785
          h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
          h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
          go to 790
          ! 785       z3 = cmplx(-r-y*h(i,na),-s-y*h(i,en)) / cmplx(zz,q)
          !           h(i+1,na) = real(z3,kind=dp)
          !           h(i+1,en) = aimag(z3)
785       call etdiv(z3r,z3i,-r-y*h(i,na),-s-y*h(i,en),zz,q)
          h(i+1,na) = z3r
          h(i+1,en) = z3i
790       continue
       enddo
       !     ********** end complex vector **********
800    continue
    enddo
    !     ********** end back substitution.
    !                vectors of isolated roots **********
    do i = 1, n
       if (i .ge. low .and. i .le. igh) go to 840
       !
       do j = i, n
          z(i,j) = h(i,j)
       enddo
       !
840    continue
    enddo
    !     ********** multiply by transformation matrix to give
    !                vectors of original full matrix.
    !                for j=n step -1 until low do -- **********
    do jj = low, n
       j = n + low - jj
       m = min0(j,igh)
       !
       do i = low, igh
          zz = 0.0_dp
          !
          do k = low, m
             zz = zz + z(i,k) * h(k,j)
          enddo
          !
          z(i,j) = zz
       enddo
    enddo
    !
    go to 1001
    !     ********** set error -- no convergence to an
    !                eigenvalue after 200 iterations **********
1000 ierr = en
1001 return
    !     ********** last card of ety2 **********
  end subroutine ety2
  subroutine etdiv(a,b,c,d,e,f)
    implicit none
    !   computes the complex division
    !     a + ib = (c + id)/(e + if)
    !  very slow, but tries to be as accurate as
    !  possible by changing the order of the
    !  operations, so to avoid under(over)flow
    !  problems.
    !  Written by F. Neri Feb. 12 1986
    !
    integer flip
    real(dp) a,b,c,d,e,f,s,t,cc,dd,ee,ff,temp
    if(.not.c_%stable_da) return

    flip = 0
    cc = c
    dd = d
    ee = e
    ff = f
    if( abs(f).ge.abs(e) ) then
       ee = f
       ff = e
       cc = d
       dd = c
       flip = 1
    endif
    s = 1.0_dp/ee
    t = 1.0_dp/(ee+ ff*(ff*s))
    if ( abs(ff) .ge. abs(s) ) then
       temp = ff
       ff = s
       s = temp
    endif
    if( abs(dd) .ge. abs(s) ) then
       a = t*(cc + s*(dd*ff))
    else if ( abs(dd) .ge. abs(ff) ) then
       a = t*(cc + dd*(s*ff))
    else
       a = t*(cc + ff*(s*dd))
    endif
    if ( abs(cc) .ge. abs(s)) then
       b = t*(dd - s*(cc*ff))
    else if ( abs(cc) .ge. abs(ff)) then
       b = t*(dd - cc*(s*ff))
    else
       b = t*(dd - ff*(s*cc))
    endif
    if (flip.ne.0 ) then
       b = -b
    endif
    return
  end subroutine etdiv
  subroutine sympl3(m)
    implicit none
    !**********************************************************
    !
    !    SYMPL3
    !
    !
    !   On return ,the matrix m(*,*), supposed to be almost
    !   symplectic on entry is made exactly symplectic by
    !   using a non iterative, constructive method.
    !
    !**********************************************************
    !
    !  Written by F. Neri  Feb 7 1986
    !
    integer,parameter::n=3
    integer kp,kq,lp,lq,jp,jq,i
    real(dp) qq,pq,qp,pp
    real(dp),dimension(2*n,2*n)::m
    if(.not.c_%stable_da) return

    !
    do kp=2,2*n,2
       kq = kp-1
       do lp=2,kp-2,2
          lq = lp-1
          qq = 0.0_dp
          pq = 0.0_dp
          qp = 0.0_dp
          pp = 0.0_dp
          do jp=2,2*n,2
             jq = jp-1
             qq = qq + m(lq,jq)*m(kq,jp) - m(lq,jp)*m(kq,jq)
             pq = pq + m(lp,jq)*m(kq,jp) - m(lp,jp)*m(kq,jq)
             qp = qp + m(lq,jq)*m(kp,jp) - m(lq,jp)*m(kp,jq)
             pp = pp + m(lp,jq)*m(kp,jp) - m(lp,jp)*m(kp,jq)
          enddo

          do i=1,2*n
             m(kq,i) = m(kq,i) - qq*m(lp,i) + pq*m(lq,i)
             m(kp,i) = m(kp,i) - qp*m(lp,i) + pp*m(lq,i)
          enddo
       enddo
       qp = 0.0_dp
       do jp=2,2*n,2
          jq = jp-1
          qp = qp + m(kq,jq)*m(kp,jp) - m(kq,jp)*m(kp,jq)
       enddo
       do i=1,2*n
          m(kp,i) = m(kp,i)/qp
       enddo
    enddo
    return
  end subroutine sympl3

  subroutine diagonalise_envelope_a(b,br,a,ai,kick)
    implicit none


    integer i
    real(dp) a(6,6),ai(6,6),b(6,6)

    real(dp) xj(6,6),xn,jb(6,6),kick(3),br(6,6)



    xj=0.0_dp

    do i=1,3
       xj(2*i,2*i-1)=-1.0_dp
       xj(2*i-1,2*i)=1.0_dp
    enddo


    jb=matmul(xj,b)
    xn=30.0_dp*mat_norm(jb)
    jb=jb/xn

    !    mj=0.0_dp
    !    xj=0.0_dp
    !    do i=1,6
    !     xj(i,i)=one
    !    enddo

    !    mj=xj+matmul(xj,jb)

    call mapflol6s(a,ai,br,jb)


    do i=1,3
       kick(i)=sqrt(abs(br(2*i-1,2*i)*xn))
    enddo



  end subroutine diagonalise_envelope_a

  subroutine mapflol6s(sa,sai,cr,cm)
    implicit none
    !---- FROM TRACKING CODE
    ! ---------------------
    integer, parameter :: ndimt=3,ndimt2=6
    integer i,ier,iunst,j,n1,n(ndimt)
    real(dp),dimension(ndimt2,ndimt2)::cr,xj,sa,sai,cm,w,vr,vi,s1
    real(dp),dimension(ndimt)::x,xx,st
    real(dp),dimension(ndimt2)::rr,ri

    if(.not.c_%stable_da) return

    n1=0
    !     frank/etienne
    do i=1,ndimt2
       do j=1,ndimt2
          cr(j,i)=cm(i,j)
       enddo
    enddo
    xj=0.0_dp
    s1=0.0_dp

    !     frank/etienne
    do i=1,ndimt
       n(i)=0
       xj(2*i-1,2*i)=1.0_dp
       xj(2*i,2*i-1)=-1.0_dp
    enddo
    !     frank/etienne


    sai=0.0_dp
    w=cm

    w=matmul(xj,w)
    w=matmul(cr,w)

    !    call mulnd2(xj,w)
    !    call mulnd2(cr,w)

    call eig6s(cr,rr,ri,vr,vi)

    do i=1,6
       write(6,*) rr(i),ri(i)
    enddo


    do i=1,ndimt
       n(i)=2*i-1
       st(i)=1.0_dp
    enddo
    !    elseif(idpr.eq.-101.or.idpr.eq.-102) then

    iunst=0
    do i=1,ndimt                 ! Frank NDC  kept
       x(i)=0.0_dp
       xx(i)=1.0_dp
       do j=1,ndimt
          x(i)=vr(2*j-1,n(i))*vi(2*j,n(i))-vr(2*j,n(i))*vi(2*j-1,n(i))+x(i)
       enddo
    enddo

    do i=1,ndimt
       if(x(i).lt.0.0_dp) xx(i)=-1.0_dp
       x(i)=SQRT(abs(x(i)))
    enddo
    do i=1,ndimt2
       do j=1,ndimt
          sai(2*j-1,i)=vr(i,n(j))*xx(j)/x(j)
          sai(2*j,i)=vi(i,n(j))/x(j)
       enddo
    enddo
    !    if(idpr.eq.-101.or.idpr.eq.-102) then
    call movearous(sai)
    !    endif
    ! adjust sa such that sa(1,2)=0 and sa(3,4)=zero (courant-snyder-edwards-teng
    ! phase advances)
    ! sa=sai


    call matinv(sai,sa,ndimt2,ndimt2,ier)

    !    call mulnd2(sai,cm)
    cm=matmul(sai,cm)
    cr=sa


    cr=matmul(cm,cr)

    !  call mulnd2(cm,cr)

    return
  end subroutine mapflol6s

  subroutine eig6s(fm,reval,aieval,revec,aievec)
    implicit none
    !**************************************************************************

    !  Diagonalization routines of NERI

    !ccccccccccccccccc
    !
    !  this routine finds the eigenvalues and eigenvectors
    !  of the full matrix fm.
    !  the eigenvectors are normalized so that the real and
    !  imaginary part of vectors 1, 3, and 5 have +1 antisymmetric
    !  product:
    !      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1 ;
    !      revec5 J aivec5 = one
    !  the eigenvectors 2 ,4, and 6 have the opposite normalization.
    !  written by F. Neri, Feb 26 1986.
    !
    integer, parameter :: ndimt=3,ndimt2=6
    integer jet,nn,i,i1,ilo,ihi,mdim,info
    real(dp),dimension(ndimt2)::reval,aieval,ort
    real(dp),dimension(ndimt2,ndimt2)::revec,aievec,fm,aa,vv

    if(.not.c_%stable_da) return

    !  copy matrix to temporary storage (the matrix aa is destroyed)
    do i=1,ndimt2
       do i1=1,ndimt2
          aa(i1,i) = fm(i1,i)
       enddo
    enddo
    ilo = 1
    ihi = ndimt2
    mdim = ndimt2
    nn = ndimt2
    !  compute eigenvalues and eigenvectors using double
    !  precision Eispack routines:
    call ety(mdim,nn,ilo,ihi,aa,ort)
    call etyt(mdim,nn,ilo,ihi,aa,ort,vv)
    call ety2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)


    !      call neigv(vv,pbkt)
    do i=1,ndimt
       do jet=1,ndimt2
          revec(jet,2*i-1)=vv(jet,2*i-1)
          revec(jet,2*i)=vv(jet,2*i-1)
          aievec(jet,2*i-1)=vv(jet,2*i)
          aievec(jet,2*i)=-vv(jet,2*i)
       enddo
    enddo
    do i=1,ndimt2
       if(abs(reval(i)**2+aieval(i)**2 -1.0_dp).gt.1e-10_dp) then
 
          if(lielib_print(4)==1) then
 
             write(6,*) sqrt(reval(i)**2+aieval(i)**2)
          endif
       endif
    enddo
    return
  end subroutine eig6s

  subroutine movearous(rt)
    implicit none
    !      integer ipause, mypause
    integer, parameter :: ndimt=3,ndimt2=6
    integer i,ic,j
    real(dp) xr,xrold
    real(dp),dimension(ndimt2,ndimt2)::rt,rto,s
    real(dp),dimension(ndimt2,ndimt2):: xt,yt,zt,xy,xz,yz
    real(dp),dimension(ndimt2,ndimt2):: xyz,xzy,xyt,yxt,yzt,zyt,xzt,zxt
    real(dp),dimension(ndimt2,ndimt2):: xyzt,xytz,xzyt,xzty,xtzy,xtyz

    if(.not.c_%stable_da) return

    do i=1,ndimt2
       do j=1,ndimt2
          s(i,j)=0.0_dp
          s(i,i)=1.0_dp
       enddo
    enddo
    xt=0.0_dp;yt=0.0_dp;zt=0.0_dp;xy=0.0_dp;xz=0.0_dp;yz=0.0_dp;
    xyzt=0.0_dp;xytz=0.0_dp;xzyt=0.0_dp;xzty=0.0_dp;xtzy=0.0_dp;xtyz=0.0_dp;
    xyz=0.0_dp;xzy=0.0_dp;xyt=0.0_dp;yxt=0.0_dp;yzt=0.0_dp;zyt=0.0_dp;xzt=0.0_dp;zxt=0.0_dp;

    do i=0,1

       xy(1+i,3+i)=1.0_dp
       xy(3+i,1+i)=1.0_dp
       xy(5+i,5+i)=1.0_dp
       !  xy(7+i,7+i)=one

       xz(1+i,5+i)=1.0_dp
       xz(5+i,1+i)=1.0_dp
       xz(3+i,3+i)=1.0_dp
       !  xz(7+i,7+i)=one

       !   xt(1+i,7+i)=one
       !   xt(7+i,1+i)=one
       xt(3+i,3+i)=1.0_dp
       xt(5+i,5+i)=1.0_dp

       yz(3+i,5+i)=1.0_dp
       yz(5+i,3+i)=1.0_dp
       yz(1+i,1+i)=1.0_dp
       !   yz(7+i,7+i)=one

       !   yt(3+i,7+i)=one
       !   yt(7+i,3+i)=one
       yt(1+i,1+i)=1.0_dp
       yt(5+i,5+i)=1.0_dp

       !   zt(5+i,7+i)=one
       !   zt(7+i,5+i)=one
       zt(1+i,1+i)=1.0_dp
       zt(3+i,3+i)=1.0_dp

       xyz(1+i,3+i)=1.0_dp
       xyz(3+i,5+i)=1.0_dp
       xyz(5+i,1+i)=1.0_dp
       !   xyz(7+i,7+i)=one

       xyz(1+i,3+i)=1.0_dp
       xyz(3+i,5+i)=1.0_dp
       xyz(5+i,1+i)=1.0_dp
       !   xyz(7+i,7+i)=one

       xzy(1+i,5+i)=1.0_dp
       xzy(5+i,3+i)=1.0_dp
       xzy(3+i,1+i)=1.0_dp
       !   xzy(7+i,7+i)=one

       xyt(1+i,3+i)=1.0_dp
       !  xyt(3+i,7+i)=one
       !   xyt(7+i,1+i)=one
       xyt(5+i,5+i)=1.0_dp

       yxt(3+i,1+i)=1.0_dp
       !   yxt(1+i,7+i)=one
       !   yxt(7+i,3+i)=one
       yxt(5+i,5+i)=1.0_dp

       yzt(3+i,5+i)=1.0_dp
       !   yzt(5+i,7+i)=one
       !   yzt(7+i,3+i)=one
       yzt(1+i,1+i)=1.0_dp

       zyt(5+i,3+i)=1.0_dp
       !   zyt(3+i,7+i)=one
       !   zyt(7+i,5+i)=one
       zyt(1+i,1+i)=1.0_dp

       xzt(1+i,5+i)=1.0_dp
       !    xzt(5+i,7+i)=one
       !    xzt(7+i,1+i)=one
       xzt(3+i,3+i)=1.0_dp

       zxt(5+i,1+i)=1.0_dp
       !    zxt(1+i,7+i)=one
       !    zxt(7+i,5+i)=one
       zxt(3+i,3+i)=1.0_dp

       xyzt(1+i,3+i)=1.0_dp
       xyzt(3+i,5+i)=1.0_dp
       !   xyzt(5+i,7+i)=one
       !   xyzt(7+i,1+i)=one

       xytz(1+i,3+i)=1.0_dp
       !   xytz(3+i,7+i)=one
       !   xytz(7+i,5+i)=one
       xytz(5+i,1+i)=1.0_dp

       xzyt(1+i,5+i)=1.0_dp
       xzyt(5+i,3+i)=1.0_dp
       !   xzyt(3+i,7+i)=one
       !   xzyt(7+i,1+i)=one

       xzty(1+i,5+i)=1.0_dp
       !   xzty(5+i,7+i)=one
       !   xzty(7+i,3+i)=one
       xzty(3+i,1+i)=1.0_dp

       !  xtzy(1+i,7+i)=one
       !  xtzy(7+i,5+i)=one
       xtzy(5+i,3+i)=1.0_dp
       xtzy(3+i,1+i)=1.0_dp

       !  xtyz(1+i,7+i)=one
       !  xtyz(7+i,3+i)=one
       xtyz(3+i,5+i)=1.0_dp
       xtyz(5+i,1+i)=1.0_dp
    enddo

    !do i=1,8
    !write(6,100) (rt(i,j),j=1,8)
    !enddo
    !write(6,*) " "
    !100  FORMAT(8(1x, E12.6))

    ic=0
    xrold=1e9_dp
    call movemuls(rt,s,rto,xr)

    if(xr.lt.xrold) then
       xrold=xr
    endif

    if(ndimt>=2) then
       call movemuls(rt,xy,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=1
       endif
    endif

    if(ndimt>=3) then
       call movemuls(rt,xz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=2
       endif
       call movemuls(rt,yz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=3
       endif
       call movemuls(rt,xyz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=4
       endif
       call movemuls(rt,xzy,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=5
       endif
    endif

    if(ndimt.eq.4) then
       call movemuls(rt,xt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=6
       endif
       call movemuls(rt,yt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=7
       endif
       call movemuls(rt,zt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=8
       endif

       call movemuls(rt,xyt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=9
       endif
       call movemuls(rt,yxt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=10
       endif
       call movemuls(rt,yzt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=11

       endif
       call movemuls(rt,zyt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=12
       endif
       call movemuls(rt,xzt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=13
       endif
       call movemuls(rt,zxt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=14
       endif
       call movemuls(rt,xyzt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=15
       endif
       call movemuls(rt,xytz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=16
       endif
       call movemuls(rt,xzyt,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=17
       endif
       call movemuls(rt,xzty,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=18
       endif
       call movemuls(rt,xtzy,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=19
       endif
       call movemuls(rt,xtyz,rto,xr)
       if(xr.lt.xrold) then
          xrold=xr
          ic=20
       endif
    endif


!    w_p=0
    i=0
    if(ic.eq.0) then
       call movemuls(rt,s,rto,xr)
       i=i+1
 !      w_p%c(i)=  " no exchanged"
    elseif(ic.eq.1) then
       call movemuls(rt,xy,rto,xr)
       i=i+1
!       w_p%c(i)=  " x-y exchanged"
    elseif(ic.eq.2) then
       call movemuls(rt,xz,rto,xr)
       i=i+1
 !      w_p%c(i)=  " x-z exchanged"
    elseif(ic.eq.3) then
       call movemuls(rt,yz,rto,xr)
       i=i+1
!       w_p%c(i)= " y-z exchanged"
    elseif(ic.eq.4) then
       call movemuls(rt,xyz,rto,xr)
       i=i+1
!       w_p%c(i)=  " x-y-z permuted"
    elseif(ic.eq.5) then
       call movemuls(rt,xzy,rto,xr)
       i=i+1
!       w_p%c(i)=  " x-z-y permuted"
    elseif(ic.eq.6) then
       call movemuls(rt,xt,rto,xr)
    elseif(ic.eq.7) then
       call movemuls(rt,yt,rto,xr)
    elseif(ic.eq.8) then
       call movemuls(rt,zt,rto,xr)
    elseif(ic.eq.9) then
       call movemuls(rt,xyt,rto,xr)
    elseif(ic.eq.10) then
       call movemuls(rt,yxt,rto,xr)
    elseif(ic.eq.11) then
       call movemuls(rt,yzt,rto,xr)
    elseif(ic.eq.12) then
       call movemuls(rt,zyt,rto,xr)
    elseif(ic.eq.13) then
       call movemuls(rt,xzt,rto,xr)
    elseif(ic.eq.14) then
       call movemuls(rt,zxt,rto,xr)
    elseif(ic.eq.15) then
       call movemuls(rt,xyzt,rto,xr)
    elseif(ic.eq.16) then
       call movemuls(rt,xytz,rto,xr)
    elseif(ic.eq.17) then
       call movemuls(rt,xzyt,rto,xr)
    elseif(ic.eq.18) then
       call movemuls(rt,xzty,rto,xr)
    elseif(ic.eq.19) then
       call movemuls(rt,xtzy,rto,xr)
    elseif(ic.eq.20) then
       call movemuls(rt,xtyz,rto,xr)
    endif





    do i=1,ndimt2
       do j=1,ndimt2
          rt(i,j)=rto(i,j)
       enddo
    enddo

    return
  end subroutine movearous

  subroutine movemuls(rt,xy,rto,xr)
    implicit none
    integer, parameter :: ndimt=3,ndimt2=6
    integer i,j,k
    real(dp) xr
    real(dp),dimension(ndimt2,ndimt2)::rt,xy,rto
    if(.not.c_%stable_da) return

    do i=1,ndimt2
       do j=1,ndimt2
          rto(i,j)=0.0_dp
       enddo
    enddo

    do  i=1,ndimt2
       do  j=1,ndimt2
          do  k=1,ndimt2
             rto(i,k)=xy(i,j)*rt(j,k)+rto(i,k)
          enddo
       enddo
    enddo

    xr=0.0_dp
    do i=1,ndimt2
       do j=1,ndimt2
          xr=xr+abs(rto(i,j))
       enddo
    enddo
    do i=1,ndimt
       xr=xr-abs(rto(2*i-1,2*i-1))
       xr=xr-abs(rto(2*i-1,2*i))
       xr=xr-abs(rto(2*i,2*i))
       xr=xr-abs(rto(2*i,2*i-1))
    enddo
    return
  end subroutine movemuls

end module lielib_yang_berz
