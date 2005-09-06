!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size
MODULE LIELIB_ETIENNE
  USE define_newda
  USE NEWDA
  USE LIELIB_BERZ
  implicit none
  private
  PUBLIC NEWTAKE,NEWDAPEK0,NEWDAPOK0,NEWETINI,NEWDACLRD,NEWDACOPD
  PUBLIC NEWDIFD,NEWINTD,NEWETCCT,NEWTRXFLO,NEWTRX,NEWFACFLOD,NEWFACFLO,NEWRTOCFLO
  PUBLIC NEWEXPFLO,NEWDALIND,NEWETINV,NEWMAPNORMF,NEWDHDJFLO,NEWFLOFACG,NEWEXPFLOD
  PUBLIC NEWFLOFAC,NEWDACMUD,NEWCTORFLO,NEWCTOR,NEWRTOC,NEWETPIN,NEWLIEINIT,NEWCOMCFU
  public newgetcct,NEWGETINV,newgtrx
  integer,PARAMETER::NDIM2=2*ndim,NTT=100
  character(120) line
CONTAINS

  SUBROUTINE newLIEINIT(NO1,NV1,ND1,NDPT1,IREF1)
    IMPLICIT none
    real(dp) ST(NDIM),ANG(NDIM),RA(NDIM)
    integer ndc1,i
    integer NO1,NV1,ND1,NDPT1,IREF1
    integer ipause,mypauses
    !     call newdaexter
    DO I=1,TOTALdol
       TABLE(I)=.FALSE.
    ENDDO
    numnewda=0
    DO I=1,NDIM
       ANG(I)=zero
       RA(I)=zero
       ST(I)=one
    enddo
    NO=NO1
    NV=NV1
    ND=ND1
    ND2=2*ND1
    call initialize_da(NO,NV)
    IF(NDPT1.EQ.0) THEN
       NDPT=0
       NDT=0
       NDC1=0
    ELSE
       NDPT=NDPT1
       NDC1=1
       IF(NDPT.EQ.ND2) THEN
          NDT=ND2-1
       ELSE
          NDT=ND2
          IF(NDPT.NE.ND2-1) THEN
             line=' LETHAL ERROR IN newLIEINIT'
             ipause=mypauses(-1,line)
          ENDIF
       ENDIF
    ENDIF
    NDC=NDC1
    NDC2=2*NDC1
    IREF=0
    CALL INITPERT(ST,ANG,RA)
    IREF=IREF1
    IF(IREF1.EQ.0) THEN
       ITU=0
    ELSE
       ITU=1
    ENDIF
    IF(IREF1.EQ.0) IREF=-1

    if(idpr.eq.1)  then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,A72))'
       write(w_p%c(1),'(a6,i4,a20)') ' NO = ',NO,' IN DA-CALCULATIONS '
       CALL WRITE_i
    endif
    do i=0,20
       xintex(i)=zero
    enddo
    xintex(0) =one
    xintex(1) =half
    xintex(2) =one/twelve
    xintex(4) =-one/c_720
    xintex(6) =one/c_30240
    xintex(8) =-one/c_1209600
    xintex(10)=one/c_21772800

    RETURN
  END subroutine newLIEINIT

  !  Z=XoY

  SUBROUTINE newETCCT(X,Y,Z)
    IMPLICIT none
    integer i
    type (taylorlow) iv(ntt)
    type (taylorlow),dimension(:)::X,Y,Z

    CALL newETALL(Iv,nv)
    DO  I=ND2+1,NV
       CALL newDAVAR(Iv(I),zero,I)
    enddo
    do  i=1,nd2
       call newdacop(y(i),iv(i))
    enddo
    CALL newDACCT(X,ND2,IV,NV,Z,ND2)

    CALL newdadal(Iv,Nv)
    RETURN
  END subroutine  newETCCT

  SUBROUTINE newgetcct(X,Y,Z,n)
    IMPLICIT none
    integer i,n
    type (taylorlow) iv(ntt)
    type (taylorlow),dimension(:)::X,Y,Z

    CALL newETALL(Iv,nv)
    DO  I=N+1,NV
       CALL newDAVAR(Iv(I),zero,I)
    enddo
    do  i=1,n
       call newdacop(y(i),iv(i))
    enddo
    CALL newDACCT(X,N,IV,NV,Z,N)

    CALL newdadal(Iv,Nv)
    RETURN
  END subroutine  newgetcct

  !  :RH: = Y :H: Y^-1 =  :HoY:

  subroutine newetmtree(y,x)
    implicit none
    integer i,nt
    type (taylorlow),dimension(:)::x,y

    nt=nv-nd2
    if(nt.gt.0) then
       do  i=nd2+1,nv
          call newdavar(x(i),zero,i)
       enddo
    endif

    do i=1,nd2
       call newdacop(y(i),x(i))
    enddo

  end   subroutine newetmtree

  subroutine newetppushlnv(x,xi,xff)
    implicit none
    integer i
    integer j(ntt)
    type (taylorlow) xt(ntt)
    type (taylorlow),dimension(:)::x
    real(dp),dimension(:)::xi,xff

    call newetall(xt,nv)
    do i=1,nv
       j(i)=0
       call newdacon(xt(i),xi(i))
    enddo

    call newdacct(x,nv,xt,nv,xt,nv)

    do i=1,nv
       call  newdapek(xt(i),j,xff(i))
    enddo

    call newdadal(xt,nv)
  end  subroutine newetppushlnv


  SUBROUTINE NEWTRX(H,RH,Y)
    IMPLICIT none
    integer i
    type (taylorlow) h,RH,iv(ntt),HH(1),RHH(1)
    type (taylorlow),dimension(:)::Y

    CALL newETALL(HH,1)
    CALL newETALL(RHH,1)
    CALL NEWDACOP(H,HH(1))
    CALL newETALL(Iv,nv)
    DO  I=ND2+1,NV
       CALL newDAVAR(Iv(I),zero,I)
    enddo
    do  i=1,nd2
       call newdacop(y(i),iv(i))
    enddo
    CALL newDACCT(HH,1,IV,NV,RHH,1)
    CALL NEWDACOP(RHH(1),RH)

    CALL newdadal(Iv,Nv)
    CALL newdadal(RHH,1)
    CALL newdadal(HH,1)

    RETURN
  END subroutine  NEWTRX

  SUBROUTINE NEWgTRX(H,RH,Y,n)
    IMPLICIT none
    integer i,n
    type (taylorlow) h,RH,iv(ntt),HH(1),RHH(1)
    type (taylorlow),dimension(:)::Y

    CALL newETALL(HH,1)
    CALL newETALL(RHH,1)
    CALL NEWDACOP(H,HH(1))
    CALL newETALL(Iv,nv)
    DO  I=n+1,NV
       CALL newDAVAR(Iv(I),zero,I)
    enddo
    do  i=1,n
       call newdacop(y(i),iv(i))
    enddo
    CALL newDACCT(HH,1,IV,NV,RHH,1)
    CALL NEWDACOP(RHH(1),RH)

    CALL newdadal(Iv,Nv)
    CALL newdadal(RHH,1)
    CALL newdadal(HH,1)

    RETURN
  END subroutine  NEWgTRX

  !  *RH! = Y *H! Y^-1  CHANGE OF A VECTOR FLOW OPERATOR

  SUBROUTINE NEWTRXFLO(H,RH,Y)
    IMPLICIT none
    integer J,K
    type (taylorlow) YI(NDIM2),HT(NDIM2),B1,B2
    type (taylorlow),dimension(:)::h,RH,Y
    !
    CALL NEWETALL(YI,ND2)
    CALL NEWETALL(HT,ND2)
    CALL NEWETALL(B1,1)
    CALL NEWETALL(B2,1)

    CALL NEWETINV(Y,YI)
    !----- HT= H o Y
    CALL NEWETCCT(H,Y,HT)
    !----
    CALL NEWDACLRD(RH)
    DO J=1,ND2
       DO K=1,ND2
          CALL NEWDADER(K,YI(J),B1)
          CALL NEWTRX(B1,B2,Y)
          CALL NEWDAMUL(B2,HT(K),B1)
          CALL NEWDAADD(B1,RH(J),B2)
          CALL NEWDACOP(B2,RH(J))
       enddo
    enddo

    CALL NEWDADAL(B2,1)
    CALL NEWDADAL(B1,1)
    CALL NEWDADAL(HT,ND2)
    CALL NEWDADAL(YI,ND2)
    RETURN
  END        SUBROUTINE NEWTRXFLO

  !  Y= AoXoAI

  SUBROUTINE NEWSIMIL(A,X,AI,Y)
    IMPLICIT NONE
    type (taylorlow),dimension(:)::X,Y,A,AI
    type (taylorlow)  W(NDIM2),V(NDIM2)
    !
    CALL NEWETALL(W,ND2)
    CALL NEWETALL(V,ND2)


    CALL NEWETCCT(A,X,W)
    CALL NEWETCCT(W,AI,V)

    CALL NEWDACOPD(V,Y)

    CALL NEWDADAL(V,ND2)
    CALL NEWDADAL(W,ND2)
    RETURN
  END   SUBROUTINE NEWSIMIL

  !  X=IDENTITY

  SUBROUTINE NEWETINI(X)
    IMPLICIT NONE
    INTEGER I
    type (taylorlow),dimension(:)::X
    !*DAEXT(NO,NV) X(NDIM2)
    DO I=1,ND2
       CALL NEWDAVAR(X(I),zero,I)
    enddo
    RETURN
  END   SUBROUTINE NEWETINI
  ! Y=X^-1

  SUBROUTINE NEWETINV(X,Y)
    IMPLICIT NONE
    integer I
    type (taylorlow) IV1(ntt),IV2(ntt)

    type (taylorlow),dimension(:)::X,Y

    CALL NEWETALL(IV1,NV)
    CALL NEWETALL(IV2,NV)
    DO I=ND2+1,NV
       CALL NEWDAVAR(IV1(I),zero,I)
    enddo
    do  i=1,nd2
       CALL NEWDACOP(X(I),IV1(I))
    ENDDO

    CALL NEWDAINV(IV1,NV,IV2,NV)
    do  i=1,nd2
       CALL NEWDACOP(IV2(I),Y(I))
    ENDDO

    CALL NEWDADAL(IV2,NV)
    CALL NEWDADAL(IV1,NV)
    RETURN
  END  SUBROUTINE NEWETINV
  !  Y=PARTIAL INVERSION OF X SEE BERZ'S PACKAGE

  SUBROUTINE NEWETPIN(X,Y,JJ)
    IMPLICIT NONE
    integer I
    integer,dimension(:)::JJ
    type (taylorlow) IV1(ntt),IV2(ntt)

    type (taylorlow),dimension(:)::X,Y


    CALL NEWETALL(IV1,NV)
    CALL NEWETALL(IV2,NV)
    DO I=ND2+1,NV
       CALL NEWDAVAR(IV1(I),zero,I)
    enddo
    do  i=1,nd2
       CALL NEWDACOP(X(I),IV1(I))
    ENDDO

    CALL NEWDAPIN(IV1,NV,IV2,NV,JJ)

    do  i=1,nd2
       CALL NEWDACOP(IV2(I),Y(I))
    ENDDO

    CALL NEWDADAL(IV2,NV)
    CALL NEWDADAL(IV1,NV)

    RETURN
  END  SUBROUTINE NEWETPIN
  !- MORE EXTENSIONS OF BASIC BERZ'S PACKAGE

  SUBROUTINE NEWgETINV(X,Y,n)
    IMPLICIT NONE
    integer I,n
    type (taylorlow) IV1(ntt),IV2(ntt)

    type (taylorlow),dimension(:)::X,Y

    CALL NEWETALL(IV1,NV)
    CALL NEWETALL(IV2,NV)
    DO I=n+1,NV
       CALL NEWDAVAR(IV1(I),zero,I)
    enddo
    do  i=1,n
       CALL NEWDACOP(X(I),IV1(I))
    ENDDO

    CALL NEWDAINV(IV1,NV,IV2,NV)
    do  i=1,n
       CALL NEWDACOP(IV2(I),Y(I))
    ENDDO

    CALL NEWDADAL(IV2,NV)
    CALL NEWDADAL(IV1,NV)
    RETURN
  END  SUBROUTINE NEWgETINV

  SUBROUTINE NEWDAPEK0(V,X,JJ)
    IMPLICIT NONE
    integer I,JJ
    type (taylorlow),dimension(:)::V
    INTEGER JD(NTT)
    real(dp),dimension(:)::X

    DO I=1,NTT
       JD(I)=0
    enddo
    DO I=1,JJ
       CALL NEWDAPEK(V(I),JD,X(I))
    enddo
    RETURN
  END  SUBROUTINE NEWDAPEK0

  SUBROUTINE NEWDAPOK0(V,X,JJ)
    IMPLICIT NONE
    integer I,JJ
    type (taylorlow),dimension(:)::V
    INTEGER JD(NTT)
    real(dp),dimension(:)::X

    DO I=1,NTT
       JD(I)=0
    enddo
    DO I=1,JJ
       CALL NEWDAPOK(V(I),JD,X(I))
    enddo
    RETURN
  END  SUBROUTINE NEWDAPOK0


  SUBROUTINE NEWDAPOKzero(V,JJ)
    IMPLICIT NONE
    integer I,JJ
    type (taylorlow),dimension(:)::V
    INTEGER JD(NTT)

    DO I=1,NTT
       JD(I)=0
    enddo
    DO I=1,JJ
       CALL NEWDAPOK(V(I),JD,zero)
    enddo
    RETURN
  END  SUBROUTINE NEWDAPOKzero

  SUBROUTINE NEWDAVAR0(V,X,JJ)
    IMPLICIT NONE
    integer I,JJ

    type (taylorlow),dimension(:)::V
    real(dp),dimension(:)::X

    DO I=1,JJ
       CALL NEWDAVAR(V(I),X(I),I)
    enddo
    RETURN
  END  SUBROUTINE NEWDAVAR0

  ! Complex dacfu

  SUBROUTINE NEWCOMCFU(B,F1,F2,C)
    IMPLICIT NONE
    real(dp) F1,F2
    EXTERNAL F1,F2
    TYPE (TAYLORLOW) T(4)
    TYPE (TAYLORLOW),dimension(:)::B,C

    CALL NEWETALL(T,4)

    CALL NEWDACFU(B(1),F1,T(1))
    CALL NEWDACFU(B(1),F2,T(2))
    CALL NEWDACFU(B(2),F1,T(3))
    CALL NEWDACFU(B(2),F2,T(4))

    CALL NEWDAsub(T(1),T(4),C(1))
    CALL NEWDAADD(T(2),T(3),C(2))
    CALL NEWDADAL(T,4)
    RETURN
  END SUBROUTINE NEWCOMCFU

  !  HT= H_M  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)

  SUBROUTINE NEWTAKE(H,M,HT)
    IMPLICIT NONE
    INTEGER I,M
    TYPE (TAYLORLOW)H,HT,B1

    CALL NEWETALL(B1,1)

    CALL NEWDACLR(B1)

    IF(M<=NO) THEN
       DO I=LOCS(M),LOCE(M)
          B1%R(I)=H%R(I)
       ENDDO
    ENDIF

    CALL NEWDACOP(B1,HT)
    CALL NEWDADAL(B1,1)
    RETURN
  END SUBROUTINE NEWTAKE

  !  \VEC{HT}= \VEC{H_M}  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)

  SUBROUTINE NEWTAKED(H,M,HT)
    IMPLICIT NONE
    INTEGER I,M
    TYPE (TAYLORLOW),dimension(:)::H,HT


    DO  I=1,ND2
       CALL NEWTAKE(H(I),M,HT(I))
    ENDDO
    RETURN
  END  SUBROUTINE NEWTAKED

  ! clear a map : a vector of nd2 polynomials
  SUBROUTINE NEWDACLRD(H)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H
    DO I=1,ND2
       CALL NEWDACLR(H(I))
    enddo
    RETURN
  END  SUBROUTINE NEWDACLRD
  !      H goes into HT  (nd2 array)
  SUBROUTINE NEWDACOPD(H,HT)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H,HT

    DO I=1,ND2
       CALL NEWDACOP(H(I),HT(I))
    enddo
    RETURN
  END SUBROUTINE NEWDACOPD

  !
  SUBROUTINE NEWDACMUD(H,SCA,HT)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H,HT
    real(dp) SCA

    DO I=1,ND2
       CALL NEWDACMU(H(I),SCA,HT(I))
    enddo
    RETURN
  END  SUBROUTINE NEWDACMUD

  SUBROUTINE NEWDALIND(H,RH,HT,RT,HR)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H,HT,HR
    real(dp) RH,RT


    DO  I=1,ND2
       CALL NEWDALIN(H(I),RH,HT(I),RT,HR(I))
    ENDDO

    RETURN
  END  SUBROUTINE NEWDALIND

  !  read a map
  SUBROUTINE NEWDAREAD(H,ND1,MFILE,XIPO)
    IMPLICIT NONE
    INTEGER I,J(NTT),MFILE,ND1
    TYPE (TAYLORLOW),dimension(:)::H
    real(dp) XIPO,RX

    DO I=1,ntt
       J(I)=0
    enddo
    do i=1,nd1
       CALL NEWDAREA(H(i),MFILE)
       CALL NEWDAPEK(H(i),j,rx)
       rx=rx*xipo
       CALL NEWDAPOK(H(i),j,rx)
    enddo
    RETURN
  END SUBROUTINE NEWDAREAD
  !  read a map
  SUBROUTINE OLDDAREAD(H,ND1,MFILE,XIPO)
    IMPLICIT NONE
    INTEGER I,J(NTT),MFILE,ND1
    TYPE (TAYLORLOW),dimension(:)::H
    real(dp) XIPO,RX

    DO I=1,ntt
       J(I)=0
    enddo
    do i=1,nd1
       CALL OLDDAREA(H(i),MFILE)
       CALL NEWDAPEK(H(i),j,rx)
       rx=rx*xipo
       CALL NEWDAPOK(H(i),j,rx)
    enddo
    RETURN
  END SUBROUTINE OLDDAREAD
  !  print a map
  SUBROUTINE NEWDAPRID(H,N1,N2,MFILE)
    IMPLICIT NONE
    INTEGER I
    INTEGER N1,N2,MFILE
    TYPE (TAYLORLOW),dimension(:)::H
    if(mfile.le.0) return
    DO I=N1,N2
       CALL NEWDAPRI(H(I),MFILE)
    enddo
    RETURN
  END  SUBROUTINE NEWDAPRID
  !  print a map
  SUBROUTINE OLDDAPRID(H,N1,N2,MFILE)
    IMPLICIT NONE
    INTEGER I
    INTEGER N1,N2,MFILE
    TYPE (TAYLORLOW),dimension(:)::H
    if(mfile.le.0) return
    DO I=N1,N2
       CALL OLDDAPRI(H(I),MFILE)
    enddo
    RETURN
  END  SUBROUTINE OLDDAPRID
  !  print a map
  SUBROUTINE NEWDAPRIresflo(H,eps,MFILE)
    IMPLICIT NONE
    INTEGER I
    INTEGER MFILE
    TYPE (TAYLORLOW) B(NDIM2),C(NDIM2)
    TYPE (TAYLORLOW),dimension(:)::H

    real(dp) EPS,dEPS

    call NEWetall(b,ND2)
    call NEWetall(c,ND2)
    call NEWdacopd(h,c)
    do i=1,nd2
       ifilt=(-1)**i
       call  NEWDACFU(c(i),filtres,h(i))
    enddo

    dEPS=-one
    CALL NEWDAEPS(dEPS)
    CALL NEWDAEPS(EPS)

    call NEWdacopd(c,h)
    call NEWdaeps(deps)
    call  NEWdadal(c,ND2)
    call  NEWdadal(b,ND2)
    RETURN
  END  SUBROUTINE NEWDAPRIresflo

  !     \VEC{H}.GRAD X =Y

  SUBROUTINE NEWDAFLO(H,X,Y)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H
    TYPE (TAYLORLOW) B1,B2,B3,X,Y
    !
    CALL NEWETALL(B1,1)
    CALL NEWETALL(B2,1)
    CALL NEWETALL(B3,1)

    CALL NEWDACLR(B1)
    CALL NEWDACLR(B2)
    DO I=1,ND2
       CALL NEWDADER(I,X,B2)
       CALL NEWDAMUL(B2,H(I),B3)
       CALL NEWDAADD(B3,B1,B2)
       CALL NEWDACOP(B2,B1)
    enddo
    CALL NEWDACOP(B1,Y)

    CALL NEWDADAL(B3,1)
    CALL NEWDADAL(B2,1)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWDAFLO

  SUBROUTINE NEWDAFLOD(H,X,Y)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H,X,Y
    TYPE (TAYLORLOW) B1(NDIM2),B2(NDIM2)
    !
    CALL NEWETALL(B1,ND2)
    CALL NEWETALL(B2,ND2)

    CALL NEWDACOPD(H,B1)
    CALL NEWDACOPD(X,B2)

    DO I=1,ND2
       CALL NEWDAFLO(B1,B2(I),Y(I))
    enddo

    CALL NEWDADAL(B1,ND2)
    CALL NEWDADAL(B2,ND2)
    RETURN
  END   SUBROUTINE NEWDAFLOD
  !     \VEC{V}.GRAD  = :H:
  ! IF SCA=1-.e0_dp
  !     \VEC{V}.GRAD  = GRAD H . GRAD
  ! IF SCA=one

  SUBROUTINE NEWINTD(V,H,SCA)
    IMPLICIT NONE
    INTEGER I
    real(dp) SCA
    TYPE (TAYLORLOW) H,B1,B2,B3,B4,X(NDIM2)
    TYPE (TAYLORLOW),dimension(:)::V

    CALL NEWETALL(B1,1)
    CALL NEWETALL(B2,1)
    CALL NEWETALL(B3,1)
    CALL NEWETALL(B4,1)
    CALL NEWETALL(X,ND2)

    CALL NEWDACLR(B4)
    CALL NEWDACLR(H)
    CALL NEWETINI(X)
    DO I=1,ND
       CALL NEWDACFU(V(2*I-1),DLIE,B3)
       CALL NEWDACFU(V(2*I),DLIE,B1)
       CALL NEWDAMUL(B1,X(2*I-1),B2)
       CALL NEWDAMUL(B3,X(2*I),B1)
       CALL NEWDALIN(B2,one,B1,SCA,B3)
       CALL NEWDAADD(B3,B4,B2)
       CALL NEWDACOP(B2,B4)
    enddo
    CALL NEWDACOP(B4,H)
    CALL NEWDADAL(X,ND2)
    CALL NEWDADAL(B4,1)
    CALL NEWDADAL(B3,1)
    CALL NEWDADAL(B2,1)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWINTD

  ! INVERSE OF INTD ROUTINE

  SUBROUTINE NEWDIFD(H1,V,SCA)
    IMPLICIT NONE
    INTEGER I
    real(dp) SCA

    TYPE (TAYLORLOW),dimension(:)::V
    TYPE (TAYLORLOW) H1,B1,h

    call NEWetall(b1,1)
    call NEWetall(h,1)
    call NEWdacop(h1,h)
    DO I=1,ND
       CALL NEWDADER(2*I-1,H,V(2*I))
       CALL NEWDADER(2*I,H,B1)
       CALL   NEWDACMU(B1,SCA,V(2*I-1))
    enddo
    CALL NEWDADAL(h,1)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWDIFD
  ! DOES EXP( \VEC{H} ) X = Y


  SUBROUTINE NEWEXPFLO(H,X,Y,EPS,NRMAX)
    IMPLICIT NONE
    INTEGER I,NRMAX
    real(dp) EPS,rbefore,COE,R
    TYPE (TAYLORLOW),dimension(:)::H
    TYPE (TAYLORLOW) B1,B2,B3,B4,X,Y
    logical(lp) more
    !
    CALL NEWETALL(B1,1)
    CALL NEWETALL(B2,1)
    CALL NEWETALL(B3,1)
    CALL NEWETALL(B4,1)

    CALL NEWDACOP(X,B4)
    CALL NEWDACOP(X,B1)
    more=.true.
    rbefore=c_1d30
    DO I=1,NRMAX
       COE=one/REAL(I,kind=DP)
       CALL NEWDACMU(B1,COE,B2)
       CALL NEWDAFLO(H,B2,B1)
       CALL NEWDAADD(B4,B1,B3)
       CALL NEWDAABS(B1,R)
       if(more) then
          if(r.gt.eps) then
             rbefore=r
             goto 100
          else
             rbefore=r
             more=.false.
          endif
       else
          IF(R.ge.rbefore) THEN
             CALL NEWDACOP(B3,Y)
             CALL NEWDADAL(B4,1)
             CALL NEWDADAL(B3,1)
             CALL NEWDADAL(B2,1)
             CALL NEWDADAL(B1,1)
             RETURN
          ENDIF
          rbefore=r
       endif
100    continue
       CALL NEWDACOP(B3,B4)
    enddo
    IF(IDPR.GE.0) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,A120))'
       write(w_p%c(1),'(a6,g20.14,a25)') ' NORM ',EPS,' NEVER REACHED IN EXPFLO '
       w_p%c(2) =  'NEW IDPR '
       CALL WRITE_i
       call read(idpr)
    ENDIF
    CALL NEWDACOP(B3,Y)
    CALL NEWDADAL(B4,1)
    CALL NEWDADAL(B3,1)
    CALL NEWDADAL(B2,1)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWEXPFLO

  ! DOES EXP( \VEC{H} ) \VEC{X} = \VEC{Y}

  SUBROUTINE NEWEXPFLOD(H,X,W,EPS,NRMAX)
    IMPLICIT NONE
    INTEGER J,NRMAX
    real(dp) EPS
    TYPE (TAYLORLOW),dimension(:)::X,W,H
    TYPE (TAYLORLOW) B0,V(NDIM2)

    CALL NEWETALL(B0,1)
    CALL NEWETALL(V,ND2)

    CALL NEWDACOPD(X,V)
    DO J=1,ND2
       CALL NEWEXPFLO(H,V(J),B0,EPS,NRMAX)
       CALL NEWDACOP(B0,V(J))
    enddo
    CALL NEWDACOPD(V,W)
    CALL NEWDADAL(V,ND2)
    CALL NEWDADAL(B0,1)
    RETURN
  END   SUBROUTINE NEWEXPFLOD

  ! IFAC=1
  ! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX ) X= Y
  ! IFAC=-1
  ! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) X= Y


  SUBROUTINE NEWFACFLO(H,X,W,NRMIN,NRMAX,SCA,IFAC)
    IMPLICIT NONE
    INTEGER NRMIN,NRMAX,IFAC,NMAX,I
    real(dp) SCA,EPS
    TYPE (TAYLORLOW),dimension(:)::H
    TYPE (TAYLORLOW) X,W,BM(NDIM2),B0(NDIM2),V
    !
    CALL NEWETALL(BM,ND2)
    CALL NEWETALL(B0,ND2)
    CALL NEWETALL(V,1)

    CALL NEWDACOP(X,V)

    EPS=-one
    CALL NEWDAEPS(EPS)
    NMAX=100
    !
    ! IFA! =1 ---> V = EXP(:SCA*H(NRMAX):)...EXP(:SCA*H(NRMIN):)X
    IF(IFAC.EQ.1) THEN
       DO I=NRMAX,NRMIN,-1
          CALL NEWTAKED(H,I,B0)
          CALL NEWDACMUD(B0,SCA,BM)

          CALL NEWEXPFLO(BM,V,B0(1),EPS,NMAX)
          CALL NEWDACOP(B0(1),V)
       enddo
    ELSE
       ! IFA! =-1 ---> V = EXP(:SCA*H(NRMIN):)...EXP(:SCA*H(NRMAX):)X
       DO I=NRMIN,NRMAX
          CALL NEWTAKED(H,I,B0)
          CALL NEWDACMUD(B0,SCA,BM)

          CALL NEWEXPFLO(BM,V,B0(1),EPS,NMAX)
          CALL NEWDACOP(B0(1),V)
       enddo
    ENDIF
    CALL NEWDACOP(V,W)
    CALL NEWDADAL(V,1)
    CALL NEWDADAL(B0,ND2)
    CALL NEWDADAL(BM,ND2)
    RETURN
  END SUBROUTINE NEWFACFLO

  ! IFAC=1
  ! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX )  \VEC{X}= \VEC{Y}
  ! IFAC=-1
  ! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) \VEC{X}= \VEC{Y}

  SUBROUTINE NEWFACFLOD(H,X,W,NRMIN,NRMAX,SCA,IFAC)
    IMPLICIT NONE
    INTEGER NRMIN,NRMAX,IFAC,I
    real(dp) SCA
    TYPE (TAYLORLOW),dimension(:)::X,W,H

    DO I=1,ND2
       CALL NEWFACFLO(H,X(I),W(I),NRMIN,NRMAX,SCA,IFAC)
    enddo

    RETURN
  END  SUBROUTINE NEWFACFLOD
  !   WRAPPED ROUTINES FOR THE OPERATOR  \VEC{H}=:H:

  ! WRAPPING FACFLOD

  SUBROUTINE NEWFEXPO(H,X,W,NRMIN,NRMAX,SCA,IFAC)
    IMPLICIT NONE
    INTEGER NRMIN,NRMAX,IFAC,NRMI,NRMA
    real(dp) SCA
    TYPE (TAYLORLOW),dimension(:)::X,W
    TYPE (TAYLORLOW) H,V(NDIM2)


    NRMI=NRMIN-1
    NRMA=NRMAX-1
    CALL NEWETALL(V,ND2)
    CALL NEWDIFD(H,V,-one)
    CALL NEWFACFLOd(V,X,W,NRMI,NRMA,SCA,IFAC)

    CALL NEWDADAL(V,ND2)

    RETURN
  END  SUBROUTINE NEWFEXPO

  ! ETCOM TAKES THE BRACKET OF TWO VECTOR FIELDS.

  SUBROUTINE NEWETCOM(X,Y,H)
    IMPLICIT NONE
    INTEGER I,J
    TYPE (TAYLORLOW),dimension(:)::H,X,Y
    TYPE (TAYLORLOW) T1,T2,T3(NDIM2)
    TYPE (TAYLORLOW) TA,TB

    CALL NEWETALL(TA,1)
    CALL NEWETALL(TB,1)
    CALL NEWETALL(T1,1)
    CALL NEWETALL(T2,1)
    CALL NEWETALL(T3,ND2)


    DO  J=1,ND2
       DO  I=1,ND2

          CALL NEWDADER(I,X(J),T1)
          CALL NEWDADER(I,Y(J),T2)
          CALL NEWDAMUL(X(I),T2,TA)
          CALL NEWDAMUL(Y(I),T1,TB)
          CALL NEWDASUB(TA,TB,T1)
          CALL NEWDAADD(T1,T3(J),TB)
          CALL NEWDACOP(TB,T3(J))

       ENDDO
    ENDDO

    CALL NEWDACOPD(T3,H)

    CALL NEWDADAL(TA,1)
    CALL NEWDADAL(TB,1)
    CALL NEWDADAL(T1,1)
    CALL NEWDADAL(T2,1)
    CALL NEWDADAL(T3,ND2)
    RETURN
  END SUBROUTINE NEWETCOM
  ! ETPOI TAKES THE POISSON BRACKET OF TWO FUNCTIONS

  SUBROUTINE NEWETPOI(X,Y,H)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW)H,X,Y,T1,T2,T3


    CALL NEWETALL(T1,1)
    CALL NEWETALL(T2,1)
    CALL NEWETALL(T3,1)


    DO I=1,ND

       CALL NEWDADER(2*I-1,X,T1)
       CALL NEWDADER(2*I,Y,T2)
       CALL NEWDAMUL(T1,T2,T1)

       CALL NEWDALIN(T1,one,T3,one,T3)
       CALL NEWDADER(2*I-1,Y,T1)
       CALL NEWDADER(2*I,X,T2)
       CALL NEWDAMUL(T1,T2,T1)

       CALL NEWDALIN(T1,-one,T3,one,T3)


    enddo

    CALL NEWDACOP(T3,H)

    CALL NEWDADAL(T1,1)
    CALL NEWDADAL(T2,1)
    CALL NEWDADAL(T3,1)
    RETURN
  END  SUBROUTINE NEWETPOI

  ! WRAPPING EXPFLO

  SUBROUTINE NEWEXP1D(H,X,Y,EPS,NON)
    IMPLICIT NONE
    real(dp) EPS
    INTEGER NON
    TYPE (TAYLORLOW)H,X,Y
    TYPE (TAYLORLOW)V(NDIM2)

    CALL NEWETALL(V,ND2)
    CALL NEWDIFD(H,V,-one)
    CALL NEWEXPFLO(V,X,Y,EPS,NON)

    CALL NEWDADAL(V,ND2)


    RETURN
  END  SUBROUTINE NEWEXP1D

  ! WRAPPING EXPFLOD USING EXP1D

  SUBROUTINE NEWEXPND2(H,X,W,EPS,NRMAX)
    IMPLICIT NONE
    INTEGER J
    real(dp) EPS
    INTEGER NRMAX
    TYPE (TAYLORLOW),dimension(:)::X,W
    TYPE (TAYLORLOW) H,B0,V(NDIM2)

    CALL NEWETALL(B0,1)
    CALL NEWETALL(V,ND2)


    CALL NEWDACOPD(X,V)
    DO J=1,ND2
       CALL NEWEXP1D(H,V(J),B0,EPS,NRMAX)
       CALL NEWDACOP(B0,V(J))
    enddo
    CALL NEWDACOPD(V,W)
    CALL NEWDADAL(V,ND2)
    CALL NEWDADAL(B0,1)
    RETURN
  END   SUBROUTINE NEWEXPND2

  ! GENERAL DRAGT-FINN FACTORIZATION

  SUBROUTINE NEWFLOFAC(XY,X,H)
    IMPLICIT NONE
    INTEGER K
    TYPE (TAYLORLOW),dimension(:)::XY,X,H
    TYPE (TAYLORLOW) V(NDIM2),W(NDIM2)

    CALL NEWETALL(V,ND2)
    CALL NEWETALL(W,ND2)

    CALL NEWDACOPD(XY,X)
    CALL NEWDACOPD(X,V)
    CALL NEWDACLRD(W)
    !      CALL danot(1)
    CALL NEWTAKED(V,1,V)
    CALL NEWETINV(V,W)
    !      CALL danot(NO)
    CALL NEWETCCT(X,W,V)
    !      CALL danot(1)
    CALL NEWTAKED(XY,1,X)
    !      CALL danot(NO)
    CALL NEWDACOPD(V,W)
    CALL NEWDACLRD(H)
    DO K=2,NO
       CALL NEWTAKED(W,K,V)
       CALL NEWDALIND(V,one,H,one,H)
       CALL NEWFACFLOD(H,W,V,K,K,-one,-1)
       CALL NEWDACOPD(V,W)
    enddo
    CALL NEWDADAL(W,ND2)
    CALL NEWDADAL(V,ND2)
    RETURN
  END   SUBROUTINE NEWFLOFAC

  ! GENERAL ONE EXPONENT FACTORIZATION

  SUBROUTINE NEWFLOFACG(XY,H,epsone)
    IMPLICIT NONE
    INTEGER K,I,NRMAX
    logical(lp) more
    real(dp) XNORM1,R,XNBEFORE,EPS,EPSONE,XX,XNORM,XN!,xintex
    TYPE (TAYLORLOW),dimension(:)::XY,H
    TYPE (TAYLORLOW) V(NDIM2),W(NDIM2),T(NDIM2), Z(NDIM2),X(NDIM2)
    integer jj(100)
    jj(1)=1
    !
    CALL NEWETALL(V,ND2 )
    CALL NEWETALL(W,ND2 )
    CALL NEWETALL(T,ND2 )
    CALL NEWETALL(Z,ND2 )
    CALL NEWETALL(X,ND2 )

    CALL NEWETINI(V)
    CALL NEWDACLRD(W)
    XNORM1=zero
    DO I=1,ND2
       CALL NEWDAABS(xy(I),R)
       XNORM1=XNORM1+R
    ENDDO
    xnbefore=c_1d30
    more=.false.
    EPS=c_1d_9
    NRMAX=1000
    XN=c_1d4
    DO K=1,NRMAX
       CALL NEWDACMUD(H,-one,T)
       CALL NEWEXPFLOD(T,XY,X,EPS,NRMAX)
       CALL NEWDALIND(X,one,V,-one,T)
       !     call olddaprid(t,1,1,20)
       IF(XN.LT.epsone) THEN
          if(idpr.ge.0) then
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,A120))'
             write(w_p%c(1),'(a14,g20.14)') " XN QUADRATIC ",XN
             CALL WRITE_i
          endif
          CALL NEWDAFLOD(T,T,W)
          CALL NEWDALIND(T,one,W,-half,T)
          call NEWdacopd(t,z)
          call NEWdacopd(t,w)
          !  second order in W
          CALL NEWETCOM(H,W,X)
          CALL NEWETCOM(X,W,X)
          !  END OF  order in W

          do I=1,3
             CALL NEWETCOM(H,W,W)
             !     call newdaprid(w,1,1,8)
             !     call olddaprid(w,1,1,8)
             CALL newDALIND(Z,one,W,XINTEX(i),Z)
          enddo
          call NEWdacopd(z,t)
          xx=one/twelve
          CALL NEWDALIND(X,xx,H,one,H)
       ENDIF

       CALL NEWDALIND(T,one,H,one,H)
       XNORM=zero
       DO I=1,ND2
          CALL NEWDAABS(T(I),R)
          XNORM=XNORM+R
       ENDDO
       XN=XNORM/XNORM1
       if(xn.ge.epsone.and.(idpr.ge.0)) then
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A120))'
          write(w_p%c(1),'(a11,g20.14)') " XN Linear ",xn
          CALL WRITE_i
       endif
       IF(XN.LT.EPS.or.more) THEN
          more=.true.
          if(xn.ge.xnbefore) goto 1000
          xnbefore=xn
       ENDIF
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A120))'
    write(w_p%c(1),'(a11,i4)') " ITERATION " , K
    CALL WRITE_i
1000 continue
    CALL NEWDADAL(W,ND2)
    CALL NEWDADAL(V,ND2)
    CALL NEWDADAL(T,ND2)
    CALL NEWDADAL(Z,ND2)
    CALL NEWDADAL(X,ND2)
    RETURN
  END  SUBROUTINE NEWFLOFACG

  ! SYMPLECTI! DRAGT-FINN FACTORIZATION WRAPPING FLOFAC

  SUBROUTINE NEWLIEFACT(XY,X,H)
    IMPLICIT NONE
    TYPE (TAYLORLOW),dimension(:)::XY,X
    TYPE (TAYLORLOW) V(NDIM2),H

    CALL NEWETALL(V,ND2)

    CALL NEWFLOFAC(XY,X,V)
    CALL NEWINTD(V,H,-one)
    !
    CALL NEWDADAL(V,ND2)

    RETURN
  END   SUBROUTINE NEWLIEFACT

  !--NORMALIZATION ROUTINES OF LIELIB


  !- WRAPPING MAPNORMF

  SUBROUTINE NEWMAPNORM(X,FT,A2,A1,XY,H,NORD)
    IMPLICIT NONE
    INTEGER ISI,NORD
    TYPE (TAYLORLOW),dimension(:)::X,A1,A2,XY
    TYPE (TAYLORLOW)H,HF(NDIM2),FTF(NDIM2),FT

    CALL NEWETALL(FTF,ND2)
    CALL NEWETALL(HF,ND2)
    ISI=0
    CALL NEWMAPNORMF(X,FTF,A2,A1,XY,HF,NORD,ISI)
    CALL NEWINTD(HF,H,-one)
    CALL NEWINTD(FTF,FT,-one)
    CALL NEWDADAL(FTF,ND2)
    CALL NEWDADAL(HF,ND2)

    RETURN
  END   SUBROUTINE NEWMAPNORM

  SUBROUTINE NEWMAPNORMF(X,FT,A2,A1,XY,H,NORD,ISI)
    IMPLICIT NONE
    real(dp) ANGLE(NDIM),ST(NDIM),p(ndim),RAD(NDIM)
    INTEGER ISI,NORD
    TYPE (TAYLORLOW),dimension(:)::X,A1,A2,FT,XY,H
    TYPE (TAYLORLOW)A1I(NDIM2),A2I(NDIM2)
    INTEGER IJ

    CALL NEWETALL(A1I,ND2)
    CALL NEWETALL(A2I,ND2)

    JTUNE=ISI
    CALL NEWDACOPD(X,XY)
    if(nv>nd2.or.ndc==1) then
       CALL NEWGOFIX(XY,A1,A1I,NORD)
       CALL NEWSIMIL(A1I,xy,A1,XY)
    else    !  this "if" was added to remove crashes when y-=plane is nearly identity
       call newetini(a1)   ! in stochastic kick calculations
       call newetini(a1i)  ! 2002.10.20
    endif

    CALL NEWMIDBFLO(XY,A2,A2I,ANGLE,RAD,ST)
    do ij=1,nd-NDC
       p(ij)=angle(ij)*(st(ij)*(twopii-one)+one)
    enddo
    IF(NDC.EQ.1) P(ND)=ANGLE(ND)
    if(idpr.ge.0) then
       w_p=1
       w_p%nc=1
       w_p%nr=2
       w_p%c(1)='tune    '
       do ij=1,nd
          w_p%r(ij)=p(ij)
       enddo
       w_p%fc='((1X,A8))'
       w_p%fr='(3(1x,g20.14))'
       CALL WRITE_i
       w_p=1
       w_p%nc=1
       w_p%nr=2
       w_p%c(1)='damping '
       do ij=1,nd
          w_p%r(ij)=rad(ij)
       enddo
       w_p%fc='((1X,A8))'
       w_p%fr='(3(1x,g20.14))'
       CALL WRITE_i
    endif
    do ij=1,nd
       !frs       pdol(ij)=p(ij)
       !frs       RADdol(ij)=RAD(ij)
       ps(ij)=p(ij)
       RADs(ij)=RAD(ij)
    enddo
    CALL INITPERT(ST,ANGLE,RAD)
    CALL NEWSIMIL(A2I,xy,A2,XY)
    CALL NEWDACOPD(XY,A2I)
    CALL NEWORDERFLO(H,FT,XY,ANGLE,RAD)
    do ij=1,nd-ndc
       p(ij)=angle(ij)
       if(angle(ij).gt.pi.and.st(ij).gt.zero.and.itu.eq.1)then
          p(ij)=angle(ij)-twopi
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(w_p%c(1),'(i4,a27,g20.14)') ij,' TH TUNE MODIFIED IN H2 TO ',p(ij)*twopii
          CALL WRITE_i
       endif
    enddo
    CALL NEWH2PLUFLO(H,P,RAD)
    !      CALL TAKED(A2I,1,XY)
    CALL NEWTAKED(A2I,1,A1I)
    CALL NEWETCCT(XY,A1I,XY)

    CALL NEWDADAL(A2I,ND2)
    CALL NEWDADAL(A1I,ND2)
    RETURN
  END  SUBROUTINE NEWMAPNORMF

  ! GETTING TO THE FIXED POINT AND CHANGING TIME APPROPRIATELY IN THE
  ! COASTING BEAM CASE


  !****************************************************************
  ! X = A1 XY A1I WHERE X IS TO THE FIXED POINT TO ORDER NORD
  ! for ndpt not zero, works in all cases. (coasting beam: eigenvalue
  !1 in Jordan form)
  !****************************************************************
  SUBROUTINE NEWGOFIX(XY,A1,A1I,NORD)
    IMPLICIT NONE
    INTEGER I,J,K,NORD
    TYPE (TAYLORLOW),dimension(:)::XY,A1,A1I
    TYPE (TAYLORLOW)X(NDIM2),W(NDIM2),V(NDIM2),REL(NDIM2),F1
    real(dp) XIC
    !
    CALL NEWETALL(X,ND2)
    CALL NEWETALL(W,ND2)
    CALL NEWETALL(V,ND2)
    CALL NEWETALL(REL,ND2)
    CALL NEWETALL(F1,1)

    ! COMPUTATION OF A1 AND A1I USING DAINV
    CALL NEWETINI(REL)

    !      CALL danot(NORD)

    CALL newETINI(V)

    DO I=1,ND2-NDC2
       CALL NEWDACOP(XY(I),X(I))
       CALL NEWDALIN(X(I),one,REL(I),-one,V(I))
    enddo

    CALL NEWETINV(V,W)
    CALL NEWDACLRD(X)

    IF(NDC.EQ.1) THEN
       CALL NEWDAVAR(X(NDPT),zero,NDPT)
    ENDIF

    CALL NEWETCCT(W,X,V)

    IF(NDC.EQ.1) THEN
       CALL NEWDACLR(V(ND2))
       CALL NEWDACLR(V(ND2-NDC))
    ENDIF

    CALL NEWDALIND(REL,one,V,one,A1)
    CALL NEWDALIND(REL,one,V,-one,A1I)


    IF(NDPT.NE.0) THEN

       !      CORRECTIONS
       call NEWDACLRD(W)
       call NEWDACLRD(V)
       call NEWDACLRD(X)

       DO I=1,ND2-NDC2
          call NEWdalin(a1(I),one,REL(I),-one,w(I))
       ENDDO

       !      COMPUTE Deta/Ddelta
       CALL NEWDACOPD(W,A1)

       DO I=1,ND2-NDC2
          CALL NEWDADER(NDPT,W(I),W(I))
       ENDDO
       !      COMPUTE J*Deta/dDELTA

       DO I=1,ND-NDC
          CALL NEWDACMU(W(2*I),one,V(2*I-1) )
          CALL NEWDACMU(W(2*I-1),-one,V(2*I) )
       ENDDO


       XIC=(-1)**(NDT)

       DO I=1,ND2-NDC2
          CALL NEWDAMUL(V(I),REL(I),X(1))
          CALL NEWDAADD(X(1),W(NDT),W(NDT))
          CALL NEWDACOP(A1(I),W(I))
       ENDDO
       CALL NEWDACMU(W(NDT),XIC,W(NDT))

       CALL NEWEXPFLOD(W,REL,A1,c_1d_7,10000)
       !     END OF  CORRECTIONS



       CALL NEWETINV(A1,A1I)
    ENDIF
    !      CALL danot(NO)
    DO I=NORD+1,NO
       DO J=LOCS(I),LOCE(I)
          DO K=1,ND2
             A1(K)%R(J)=zero
             A1I(K)%R(J)=zero
          ENDDO
       ENDDO
    ENDDO

    CALL NEWDADAL(F1,1)
    CALL NEWDADAL(REL,ND2)
    CALL NEWDADAL(V,ND2)
    CALL NEWDADAL(W,ND2)
    CALL NEWDADAL(X,ND2)
    RETURN
  END     SUBROUTINE NEWGOFIX

  !-   NONLINEAR NORMALIZATION PIECE OF MAPNORMF

  SUBROUTINE NEWORDERFLO(H,FT,X,ANG,RA)
    IMPLICIT NONE
    real(dp),dimension(:)::ANG,RA
    INTEGER K
    TYPE (TAYLORLOW),dimension(:)::X,H,FT
    TYPE (TAYLORLOW)W(NDIM2),V(NDIM2),REL(NDIM2)
    TYPE (TAYLORLOW)ROI(NDIM2)
    TYPE (TAYLORLOW)B1(NDIM2),B5(NDIM2),B6(NDIM2),B9(NDIM2)
    !
    CALL NEWETALL(W,ND2)
    CALL NEWETALL(V,ND2)
    CALL NEWETALL(REL,ND2)
    CALL NEWETALL(ROI,ND2)
    CALL NEWETALL(B1,ND2)
    CALL NEWETALL(B5,ND2)
    CALL NEWETALL(B6,ND2)
    CALL NEWETALL(B9,ND2)

    CALL NEWROTIFLO(ROI,ANG,RA)
    CALL NEWETINI(REL)
    CALL NEWDACLRD(H)
    CALL NEWDACLRD(FT)
    CALL NEWETCCT(X,ROI,X)
    DO K=2,NO
       ! IF K>2 V = H(K)^-1 X(K)
       CALL NEWFACFLOD(H,X,V,2,K-1,-one,-1)
       ! EXTRACTING K TH DEGREE OF V ----> W
       CALL NEWTAKED(V,K,W)
       !      write(16,*) "$$$$$$$$  K  $$$$$$$$$$", k

       ! W = EXP(B5) + ...
       CALL NEWDACOPD(W,B5)
       !      CALL INTD(W,B5,-one)
       ! B5 ON EXIT IS THE NEW CONTRIBUTION TO H
       ! B6 IS THE NEW CONTRIBUTION TO FT
       call NEWnuanaflo(b5,b6)
       CALL NEWDALIND(B5,one,H,one,B1)
       CALL NEWDACOPD(B1,H)
       ! EXP(B9) = EXP( : ROTI B6 :)
       CALL NEWTRXFLO(B6,B9,ROI)

       ! V = EXP(-B6) REL
       CALL NEWFACFLOD(B6,REL,V,K,K,-one,1)
       ! W = V o X
       CALL NEWETCCT(V,X,W)
       if(idpr.ge.0) then
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,A72),/))'
          write(w_p%c(1),'(a13,i4)') ' ORDERFLO K= ', k
          CALL WRITE_i
       endif
       ! X = EXP(B9) W
       CALL NEWFACFLOD(B9,W,X,K,K,one,1)
       ! B6 IS THE NEW CONTRIBUTION TO FT
       CALL NEWDALIND(B6,one,FT,one,B1)
       CALL NEWDACOPD(B1,FT)
    enddo
34  CALL NEWDADAL(B9,ND2)
    CALL NEWDADAL(B6,ND2)
    CALL NEWDADAL(B5,ND2)
    CALL NEWDADAL(B1,ND2)
    CALL NEWDADAL(ROI,ND2)
    CALL NEWDADAL(REL,ND2)
    CALL NEWDADAL(V,ND2)
    CALL NEWDADAL(W,ND2)
    RETURN
  END   SUBROUTINE NEWORDERFLO

  ! RESONANCE DENOMINATOR OPERATOR (1-R^-1)^-1

  SUBROUTINE NEWNUANAFLO(H,FT)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H,FT
    TYPE (TAYLORLOW)T1(2),T2(2),Br(NDIM2),BI(NDIM2),C(NDIM2),CI(NDIM2)

    CALL NEWETALL(BR,ND2)
    CALL NEWETALL(BI,ND2)
    CALL NEWETALL(C,ND2)
    CALL NEWETALL(CI,ND2)
    CALL NEWETALL(T1,2)
    CALL NEWETALL(T2,2)

    CALL NEWCTORFLO(H,BR,BI)
    ! FILTERING RESONANCES AND TUNE SHIFTS
    ! ASSUMING REALITY I.E. B(2*I-1)=CMPCJG(B(2*I))


    DO I=1,ND2
       IFLOW=I
       CALL NEWDACFU(BR(I),FILT,C(I))
       CALL NEWDACFU(BI(I),FILT,CI(I))
    enddo

    CALL NEWRTOCFLO(C,CI,H)

    DO I=1,ND2

       IFLOW=I
       CALL NEWDACFU(BR(i),DFILT,BR(i))
       CALL NEWDACFU(BI(i),DFILT,BI(i))
    enddo
    !  NOW WE MUST REORDER ! AND CI TO SEPARATE THE REAL AND IMAGINARY PART
    ! THIS IS NOT NECESSARY WITH :H: OPERATORS


    DO I=1,ND2
       CALL NEWDACOP(BR(I),T1(1))
       CALL NEWDACOP(BI(I),T1(2))
       CALL NEWDACOP(C(I),T2(1))
       CALL NEWDACOP(CI(I),T2(2))
       IFLOW=I
       CALL newCOMCFU(T1,XGAM,XGBM,T2)
       CALL NEWDACOP(T2(1),C(I))
       CALL NEWDACOP(T2(2),CI(I))
    enddo

    CALL NEWRTOCFLO(C,CI,FT)


    CALL NEWDADAL(BR,ND2)
    CALL NEWDADAL(BI,ND2)
    CALL NEWDADAL(C,ND2)
    CALL NEWDADAL(CI,ND2)
    CALL NEWDADAL(T1,2)
    CALL NEWDADAL(T2,2)

    RETURN
  END  SUBROUTINE NEWNUANAFLO

  !  SYMPLECTIC SHIFTS OR SINE TERMS!
  SUBROUTINE NEWDHDJFLO(H,T)
    implicit none
    integer i
    TYPE (TAYLORLOW),dimension(:)::H,T

    TYPE (TAYLORLOW) B1(NDIM2),B2(NDIM2),BB1,BB2

    CALL NEWETALL(B1,nd2)
    CALL NEWETALL(B2,nd2)
    CALL NEWETALL(Bb1,1)
    CALL NEWETALL(Bb2,1)

    CALL NEWCTORFLO(H,B1,B2)

    DO I=1,ND-NDC
       CALL NEWDATRA(2*I,B2(2*I),BB1)
       CALL NEWDACMU(BB1,twopii,T(I+ND))
       CALL NEWDACOP(T(I+ND),BB1)
       CALL NEWDACLR(BB2)
       call NEWrtoc(BB1,BB2,BB1)
       CALL NEWDACOP(BB1,T(I))
    enddo

    IF(NDPT.ne.0) THEN
       CALL NEWDACOP(H(NDT),T(ND))
       CALL NEWDACOP(B1(NDT),T(ND2))
    endif

    CALL NEWDADAL(BB2,1)
    CALL NEWDADAL(BB1,1)
    CALL NEWDADAL(B2,ND2)
    CALL NEWDADAL(B1,ND2)
    RETURN
  END SUBROUTINE NEWDHDJFLO

  ! POKES IN \VEC{H}  ANGLES AND DAMPING COEFFFICIENTS

  SUBROUTINE NEWH2PLUFLO(H,ANG,RA)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::H
    real(dp) R1,R2,ST(NDIM)
    real(dp),dimension(:)::ANG,RA
    INTEGER J(NTT)

    DO I=1,ND
       ST(I)=two*STA(I)-one
    enddo

    DO I=1,NTT
       J(I)=0
    enddo

    DO I=1,ND-NDC
       J(2*I-1)=1
       R1=-ANG(I)
       !-----
       CALL NEWDAPOK(H(2*I),J,R1)

       R2=RA(I)
       CALL NEWDAPOK(H(2*I-1),J,R2)
       J(2*I-1)=0

       J(2*I)=1
       R1=ANG(I)*ST(I)
       CALL NEWDAPOK(H(2*I-1),J,R1)
       CALL NEWDAPOK(H(2*I),J,R2)
       J(2*I)=0

    enddo

    IF(NDPT.EQ.ND2-1) THEN
       J(NDPT)=1
       CALL NEWDAPOK(H(NDT),J,ANG(ND))
    ELSEif(ndpt.eq.nd2) then
       J(NDPT)=1
       CALL NEWDAPOK(H(NDT),J,-ANG(ND))
    ENDIF
    RETURN
  END   SUBROUTINE NEWH2PLUFLO

  ! CREATES R AND R^-1 USING THE EXISTING ANGLES AND DAMPING
  ! COULD BE REPLACED BY A CALL H2PLUFLO FOLLOWED BY EXPFLOD

  ! CREATES R

  SUBROUTINE NEWROTFLO(RO,ANG,RA)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::RO
    real(dp) CO(NDIM),SI(NDIM),XX,CH,SH,SIM
    real(dp),dimension(:)::ANG,RA
    INTEGER J(NTT)

    do i=1,ntt
       j(i)=0
    enddo

    CALL NEWDACLRD(RO)
    DO I=1,ND-NDC
       XX=EXP(RA(I))
       IF(ISTA(I).EQ.0) THEN
          CALL HYPER(ANG(I),CH,SH)
          CO(I)=CH*XX
          SI(I)=-SH*XX
       ELSE
          CO(I)=COS(ANG(I))*XX
          SI(I)=SIN(ANG(I))*XX
       ENDIF
    enddo
    DO I=1,ND-NDC
       IF(ISTA(I).EQ.0)THEN
          SIM=SI(I)
       ELSE
          SIM=-SI(I)
       ENDIF
       J(2*I-1)=1
       CALL NEWDAPOK(RO(2*I-1),J,CO(I))
       CALL NEWDAPOK(RO(2*I),J,SIM)
       J(2*I-1)=0
       J(2*I)=1
       CALL NEWDAPOK(RO(2*I),J,CO(I))
       CALL NEWDAPOK(RO(2*I-1),J,SI(I))
       J(2*I)=0
    enddo

    IF(NDC.EQ.1) THEN
       J(NDT)=1
       CALL NEWDAPOK(RO(NDT),J,one)
       CALL NEWDAPOK(RO(NDPT),J,zero)
       J(NDT)=0
       J(NDPT)=1
       CALL NEWDAPOK(RO(NDT),J,ANG(ND))
       CALL NEWDAPOK(RO(NDPT),J,one)
       J(NDPT)=0
    ENDIF

    RETURN
  END  SUBROUTINE NEWROTFLO

  ! CREATES  R^-1

  SUBROUTINE NEWROTIFLO(ROI,ANG,RA)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::ROI
    real(dp) CO(NDIM),SI(NDIM),XX,CH,SH,SIM,SIMV
    real(dp),dimension(:)::ANG,RA
    INTEGER J(NTT)

    do i=1,ntt
       j(i)=0
    enddo


    CALL NEWDACLRD(ROI)
    DO I=1,ND-NDC
       XX=EXP(-RA(I))
       IF(ISTA(I).EQ.0) THEN
          CALL HYPER(ANG(I),CH,SH)
          CO(I)=CH*XX
          SI(I)=-SH*XX
       ELSE
          CO(I)=COS(ANG(I))*XX
          SI(I)=SIN(ANG(I))*XX
       ENDIF
    enddo
    DO I=1,ND-NDC
       IF(ISTA(I).EQ.0)THEN
          SIM=SI(I)
       ELSE
          SIM=-SI(I)
       ENDIF
       J(2*I-1)=1
       CALL NEWDAPOK(ROI(2*I-1),J,CO(I))
       SIMV=-SIM
       CALL NEWDAPOK(ROI(2*I),J,SIMV)
       J(2*I-1)=0
       J(2*I)=1
       SIMV=-SI(I)

       CALL NEWDAPOK(ROI(2*I),J,CO(I))
       CALL NEWDAPOK(ROI(2*I-1),J,SIMV)
       J(2*I)=0
    enddo

    IF(NDC.EQ.1) THEN
       J(NDT)=1
       CALL NEWDAPOK(ROI(NDT),J,one)
       CALL NEWDAPOK(ROI(NDPT),J,zero)
       J(NDT)=0
       J(NDPT)=1
       CALL NEWDAPOK(ROI(NDT),J,-ANG(ND))
       CALL NEWDAPOK(ROI(NDPT),J,one)
       J(NDPT)=0
    ENDIF

    RETURN
  END   SUBROUTINE NEWROTIFLO

  ! CHANGES OF BASIS


  !   C1------> R2+I R1

  SUBROUTINE NEWCTOR(C1,R2,I2)
    IMPLICIT NONE
    TYPE (TAYLORLOW)C1,R2,I2
    TYPE (TAYLORLOW)B1,B2,X(NDIM2)
    !
    CALL NEWETALL(B1,1)
    CALL NEWETALL(B2,1)
    CALL NEWETALL(X,ND2)
    CALL NEWCTOI(C1,B1)
    CALL NEWETCJG(X)
    CALL NEWTRX(B1,B2,X)
    CALL NEWDALIN(B1,half,B2,half,R2)
    CALL NEWDALIN(B1,half,B2,-half,I2)
    CALL NEWDADAL(X,ND2)
    CALL NEWDADAL(B2,1)
    CALL NEWDADAL(B1,1)
    RETURN
  END   SUBROUTINE NEWCTOR

  !  INVERSE OF CTOR

  SUBROUTINE NEWRTOC(R1,I1,C2)
    IMPLICIT NONE
    TYPE (TAYLORLOW)C2,R1,I1
    TYPE (TAYLORLOW)B1
    !
    CALL NEWETALL(B1,1)

    CALL NEWDAADD(R1,I1,B1)
    CALL NEWITOC(B1,C2)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWRTOC

  ! FLOW CTOR

  SUBROUTINE NEWCTORFLO(C,DR,DI)
    IMPLICIT NONE
    TYPE (TAYLORLOW),dimension(:)::DR,DI,C
    CALL NEWCTORD(C,DR,DI)
    CALL NEWRESVEC(DR,DI,DR,DI)
    RETURN
  END  SUBROUTINE NEWCTORFLO

  ! FLOW RTOC


  SUBROUTINE NEWRTOCFLO(DR,DI,C)
    IMPLICIT NONE
    TYPE (TAYLORLOW),dimension(:)::DR,DI,C
    TYPE (TAYLORLOW) ER(NDIM2),EI(NDIM2)

    CALL NEWETALL(ER,ND2)
    CALL NEWETALL(EI,ND2)

    CALL NEWREELFLO(DR,DI,ER,EI)
    CALL NEWRTOCD(ER,EI,C)
    CALL NEWDADAL(ER,ND2)
    CALL NEWDADAL(EI,ND2)

    RETURN
  END  SUBROUTINE NEWRTOCFLO

  ! ROUTINES USED IN THE INTERMEDIATE STEPS OF CTORFLO AND RTOCFLO

  ! SAME AS CTOR  OVER ARRAYS CONTAINING ND2 COMPONENTS
  ! ROUTINE USEFUL IN INTERMEDIATE FLOW CHANGE OF BASIS

  SUBROUTINE NEWCTORD(C,CR,CI)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::C,CI,CR

    DO I=1,ND2
       CALL NEWCTOR(C(I),CR(I),CI(I))
    enddo
    RETURN
  END   SUBROUTINE NEWCTORD

  !  INVERSE OF CTORD

  SUBROUTINE NEWRTOCD(CR,CI,C)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::C,CI,CR
    DO I=1,ND2
       CALL newRTOC(CR(I),CI(I),C(I))
    enddo
    RETURN
  END  SUBROUTINE NEWRTOCD

  ! DOES THE SPINOR PART IN CTORFLO

  SUBROUTINE NEWRESVEC(CR,CI,DR,DI)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW) TR(2),TI(2)
    TYPE (TAYLORLOW),dimension(:)::DR,DI,CI,CR

    CALL NEWETALL(TR,2)
    CALL NEWETALL(TI,2)

    DO I=1,ND-ndc
       if(ista(i).eq.1) then
          CALL NEWDASUB(CR(2*I-1),CI(2*I),TR(1))
          CALL NEWDAADD(CI(2*I-1),CR(2*I),TI(1))
          CALL NEWDAADD(CR(2*I-1),CI(2*I),TR(2))
          CALL NEWDASUB(CI(2*I-1),CR(2*I),TI(2))
          CALL NEWDACOP(TR(1),DR(2*I-1))
          CALL NEWDACOP(TR(2),DR(2*I))
          CALL NEWDACOP(TI(1),DI(2*I-1))
          CALL NEWDACOP(TI(2),DI(2*I))
       else
          CALL NEWDAadd(CR(2*I-1),Cr(2*I),TR(1))
          CALL NEWDAADD(CI(2*I-1),Ci(2*I),TI(1))
          CALL NEWDAsub(CR(2*I-1),Cr(2*I),TR(2))
          CALL NEWDASUB(CI(2*I-1),Ci(2*I),TI(2))
          CALL NEWDACOP(TR(1),DR(2*I-1))
          CALL NEWDACOP(TR(2),DR(2*I))
          CALL NEWDACOP(TI(1),DI(2*I-1))
          CALL NEWDACOP(TI(2),DI(2*I))
       endif
    enddo

    do i=nd2-ndc2+1,nd2
       call NEWdacop(cr(i),dr(i))
       call NEWdacop(ci(i),di(i))
    enddo

    CALL NEWDADAL(TR,2)
    CALL NEWDADAL(TI,2)
    RETURN
  END   SUBROUTINE NEWRESVEC

  ! DOES THE SPINOR PART IN RTOCFLO

  SUBROUTINE NEWREELFLO(C,CI,F,FI)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW) E(NDIM2),EI(NDIM2)
    TYPE (TAYLORLOW),dimension(:)::C,CI,F,FI

    CALL NEWETALL(E,ND2)
    CALL NEWETALL(EI,ND2)


    DO I=1,ND-NDC
       CALL NEWDALIN(C(2*i-1),half,C(2*i),half,E(2*i-1))
       CALL NEWDALIN(CI(2*i-1),half,CI(2*I),half,EI(2*i-1))
       if(ista(i).eq.1) then
          CALL NEWDALIN(CI(2*i-1),half,CI(2*i),-half,E(2*i))
          CALL NEWDALIN(C(2*i-1),-half,C(2*I),half,EI(2*i))
       ELSE
          CALL NEWDALIN(CI(2*i-1),half,CI(2*i),-half,EI(2*i))
          CALL NEWDALIN(C(2*i-1),half,C(2*I),-half,E(2*i))
       ENDIF
    enddo

    do i=nd2-ndc2+1,nd2
       call NEWdacop(c(i),e(i))
       call NEWdacop(ci(i),ei(i))
    enddo

    CALL NEWDACOPD(E,F)
    CALL NEWDACOPD(EI,FI)



    CALL NEWDADAL(E,ND2)
    CALL NEWDADAL(EI,ND2)
    RETURN
  END   SUBROUTINE NEWREELFLO

  ! TAKES THE COMPLEX CONJUGATE IN RESONANCE BASIS OF A POLYNOMIAL


  SUBROUTINE NEWCOMPCJG(CR,CI,DR,DI)
    IMPLICIT NONE
    TYPE (TAYLORLOW)DR,DI,CI,CR,X(NDIM2)

    CALL NEWETALL(X,ND2)


    CALL NEWETCJG(X)
    CALL NEWTRX(CR,DR,X)
    CALL NEWTRX(CI,DI,X)
    CALL NEWDACMU(DI,-one,DI)

    CALL NEWDADAL(X,ND2)
    RETURN
  END  SUBROUTINE NEWCOMPCJG

  subroutine NEWmidBFLO(c,a2,a2i,q,a,ST)
    IMPLICIT NONE
    real(dp) R,CH,SHM
    INTEGER jx(ntt)
    INTEGER I,J
    real(dp) cr(ndim2,ndim2),ST(NDIM),q(ndim),a(ndim)
    real(dp) sa(ndim2,ndim2),sai(ndim2,ndim2),cm(ndim2,ndim2)
    TYPE (TAYLORLOW),dimension(:)::C,A2,A2I
    !*DAEXT(NO,NV) C(NDIM2),A2(NDIM2),A2I(NDIM2)

    do i=1,NTT
       jx(i)=0
    enddo

    do i=1,nd2
       do j=1,nd2
          sai(i,j)=zero
          sa(i,j)=zero
          CM(I,J)=zero
          CR(I,J)=zero
       enddo
    enddo

    do i=1,nd2
       do j=1,nd2
          jx(j)=1
          call  NEWdapek(c(i),jx,r)
          jx(j)=0
          cm(i,j)=r
       enddo
    enddo


    call MAPFLOL(SA,SAI,CR,CM,st)

    do i=1,nd-ndc
       if(st(i)+c_1d_3.gt.one) then
          a(i)=SQRT(cr(2*i-1,2*i-1)**2+cr(2*i-1,2*i)**2)
          q(i)=ACOS(cr(2*i-1,2*i-1)/a(i))
          A(I)=LOG(A(I))
          if(cr(2*i-1,2*i).lt.zero) q(i)=twopi-q(i)
       ELSE
          a(i)=SQRT(cr(2*i-1,2*i-1)**2-cr(2*i-1,2*i)**2)
          CH=cr(2*i-1,2*i-1)/a(i)
          sHm=cr(2*i-1,2*i)/a(i)
          !       CH=CH+SQRT(CH**2-one)
          !       q(i)=LOG(CH)
          q(i)=-LOG(ch+shm)
          !       IF(cr(2*i-1,2*i).gt.zero) Q(I)=-Q(I)
          A(I)=LOG(A(I))
       endif
    enddo


    IF(NDC.EQ.0) THEN
       IF(st(3)+c_1d_3.gt.one.and.ND.EQ.3.AND.Q(nd).GT.half) Q(3)=Q(3)-twopi
    ELSE
       Q(ND)=CR(NDT,NDPT)
    ENDIF

    CALL NEWDACLRD(A2)
    CALL NEWDACLRD(A2I)

    do i=1,nd2
       do j=1,nd2
          jx(j)=1
          r=sa(i,j)
          if(r.ne.zero)call  NEWdapok(a2(i),jx,r)
          jx(j)=1
          r=sai(i,j)
          if(r.ne.zero)call  NEWdapok(a2i(i),jx,r)
          jx(j)=0
       enddo
    enddo

    return
  end    subroutine NEWmidBFLO

  SUBROUTINE NEWCPART(H,CH)
    IMPLICIT NONE
    TYPE (TAYLORLOW) H,CH

    CALL NEWDACFU(H,REXT,CH)
    RETURN
  END  SUBROUTINE NEWCPART

  SUBROUTINE NEWCTOI(F1,F2)
    IMPLICIT NONE
    TYPE (TAYLORLOW) F1,F2
    TYPE (TAYLORLOW) B1,X(NDIM2)
    !
    CALL NEWETALL(B1,1)
    CALL NEWETALL(X,ND2)

    CALL NEWCPART(F1,B1)
    CALL NEWETCTR(X)
    CALL NEWTRX(B1,F2,X)
    CALL NEWDADAL(X,ND2)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWCTOI

  SUBROUTINE NEWITOC(F1,F2)
    IMPLICIT NONE
    TYPE (TAYLORLOW) F1,F2
    TYPE (TAYLORLOW) B1,X(NDIM2)
    !
    CALL NEWETALL(B1,1)
    CALL NEWETALL(X,ND2)

    CALL NEWETRTC(X)
    CALL NEWTRX(F1,B1,X)
    CALL NEWCPART(B1,F2)

    CALL NEWDADAL(X,ND2)
    CALL NEWDADAL(B1,1)
    RETURN
  END  SUBROUTINE NEWITOC

  SUBROUTINE NEWETRTC(X)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::X
    TYPE (TAYLORLOW) REL(NDIM2)
    !
    CALL NEWETALL(REL,ND2)

    CALL NEWETINI(REL)
    CALL NEWETINI(X)
    DO I=1,ND-NDC
       CALL NEWDAADD(REL(2*I-1),REL(2*I),X(2*I-1))
       CALL NEWDASUB(REL(2*I-1),REL(2*I),X(2*I))
    enddo
    CALL NEWDADAL(REL,ND2)
    RETURN
  END SUBROUTINE NEWETRTC

  SUBROUTINE NEWETCTR(X)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::X
    TYPE (TAYLORLOW) REL(NDIM2)
    !
    CALL NEWETALL(REL,ND2)

    CALL NEWETINI(REL)
    CALL NEWETINI(X)
    DO I=1,ND-NDC
       CALL NEWDALIN(REL(2*I-1),half,REL(2*I),half,X(2*I-1))
       CALL NEWDALIN(REL(2*I-1),half,REL(2*I),-half,X(2*I))
    enddo
    CALL NEWDADAL(REL,ND2)
    RETURN
  END   SUBROUTINE NEWETCTR

  SUBROUTINE NEWETCJG(X)
    IMPLICIT NONE
    INTEGER I
    TYPE (TAYLORLOW),dimension(:)::X
    TYPE (TAYLORLOW) REL(NDIM2)
    !
    CALL NEWETALL(REL,ND2)

    CALL NEWETINI(REL)
    CALL NEWETINI(X)
    DO I=1,ND-NDC
       IF(ISTA(I).EQ.1) THEN
          CALL NEWDACOP(REL(2*I-1),X(2*I))
          CALL NEWDACOP(REL(2*I),X(2*I-1))
       ELSE
          CALL NEWDACOP(REL(2*I-1),X(2*I-1))
          CALL NEWDACOP(REL(2*I),X(2*I))
       ENDIF
    enddo
    CALL NEWDADAL(REL,ND2)
    RETURN
  END   SUBROUTINE NEWETCJG



  SUBROUTINE NEWd(X,r)
    IMPLICIT NONE
    INTEGER I
    real(dp) r,r1
    TYPE (TAYLORLOW),dimension(:)::X

    !
    !
    r=zero
    DO  I=1,ND2
       call newdaAbs(x(i),r1)
       r=r1+r
    enddo

    RETURN
  END   SUBROUTINE NEWd


END MODULE LIELIB_ETIENNE
