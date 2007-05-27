!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_extend_poly
  USE tree_element_MODULE
  IMPLICIT NONE
  public
  integer,private,parameter::ndd=6
  private real_8REAL6,REAL6real_8,env_8map,real_8REAL_8,env_8benv
  private print6,PRINTenv
  private ALLOCenv6,killenv6,REAL6env_8,env_8t,tenv_8,scdadd,daddsc
  logical(lp), target :: ALWAYS_knobs=.false.

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE REAL_8REAL6
     MODULE PROCEDURE REAL6REAL_8
     MODULE PROCEDURE ENV_8MAP
     MODULE PROCEDURE REAL6env_8
     MODULE PROCEDURE real_8REAL_8
     MODULE PROCEDURE env_8t
     MODULE PROCEDURE tenv_8
     MODULE PROCEDURE env_8benv
  END  INTERFACE

  INTERFACE operator (+)
     MODULE PROCEDURE scdadd
     MODULE PROCEDURE daddsc
  END  INTERFACE

  INTERFACE PRINT
     MODULE PROCEDURE PRINT6
     MODULE PROCEDURE PRINTenv
  END  INTERFACE

  INTERFACE DAPRINT
     MODULE PROCEDURE PRINT6
     MODULE PROCEDURE PRINTenv
  END  INTERFACE



CONTAINS



  SUBROUTINE ANALYSE_APERTURE_FLAG(I,R)
    IMPLICIT NONE
    INTEGER I,B,K
    INTEGER :: R(:)

    K=I
    B=1
    r=-1
    DO WHILE (K>0.AND.B<=SIZE(R))
       R(B)=MOD(K,2)
       IF(MOD(K,2)==1) THEN
          K=(K-1)/2
       ELSE
          K=K/2
       ENDIF
       B=B+1
    ENDDO

  END   SUBROUTINE ANALYSE_APERTURE_FLAG

  SUBROUTINE PRODUCE_APERTURE_FLAG(I)
    IMPLICIT NONE
    INTEGER I,B
    I=0
    B=1
    IF(.NOT.c_%CHECK_STABLE) THEN
       I=B+I
    ENDIF
    B=B*2
    IF(.NOT.c_%CHECK_MADX_APERTURE) THEN
       I=B+I
    ENDIF
    B=B*2
    !frs 10.03.2006    IF(.NOT.c_%check_iteration) THEN
    !frs 10.03.2006       I=B+I
    !frs 10.03.2006    ENDIF
    B=B*2
    !frs 10.03.2006    IF(.NOT.c_%check_interpolate_x) THEN
    !frs 10.03.2006       I=B+I
    !frs 10.03.2006    ENDIF
    B=B*2
    !frs 10.03.2006    IF(.NOT.c_%check_interpolate_y) THEN
    !frs 10.03.2006       I=B+I
    !frs 10.03.2006    ENDIF
    B=B*2
    IF(.NOT.c_%check_x_min) THEN
       I=B+I
    ENDIF
    B=B*2
    IF(.NOT.c_%check_x_max) THEN
       I=B+I
    ENDIF
    B=B*2
    IF(.NOT.c_%check_y_min) THEN
       I=B+I
    ENDIF
    B=B*2
    IF(.NOT.c_%check_y_max) THEN
       I=B+I
    ENDIF
    B=B*2
    IF(.not.c_%stable_da) THEN
       I=B+I
    ENDIF
    B=B*2

    ALLOW_TRACKING=.TRUE.

  END   SUBROUTINE PRODUCE_APERTURE_FLAG





  REAL(DP) FUNCTION  SINEHX_X(X) ! REPLACES SINH(X)/X
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.CHECK_STABLE) return

    IF((ABS(X)>hyperbolic_aperture).AND.ROOT_CHECK) THEN
       SINEHX_X=ZERO
       CHECK_STABLE=.FALSE.
    ELSEIF(ABS(X)<=hyperbolic_aperture) THEN
       sinehx_x = sinhx_x(x)
    ELSE      !  IF X IS NOT A NUMBER
       sinehx_x=ZERO
       CHECK_STABLE=.FALSE.
    ENDIF

  END FUNCTION SINEHX_X



  ! Some polymorphism

  SUBROUTINE  print6(S1,mf)
    implicit none
    type (real_8),INTENT(INout)::S1(ndd)
    integer        mf,i

    do i=1,ndd
       call print(s1(i),mf)
    enddo

  END SUBROUTINE print6

  SUBROUTINE  printenv(S1,mf)
    implicit none
    type (env_8),INTENT(INout)::S1(ndd)
    integer        mf,i,j
    do i=1,ndd
       do j=1,ndd
          write(mf,*) "Sigma0 ",i,j
          call print(s1(i)%sigma0(j),mf)
       enddo
    enddo

    do i=1,ndd
       do j=1,ndd
          write(mf,*) "Sigmaf ",i,j
          call print(s1(i)%sigmaf(j),mf)
       enddo
    enddo

    write(mf,*) "Map "
    do i=1,ndd
       call print(s1(i)%v,mf)
    enddo

    do i=1,ndd
       do j=1,ndd
          write(mf,*) "dB ",i,j
          call print(s1(i)%e(j),mf)
       enddo
    enddo


  END SUBROUTINE printenv

  SUBROUTINE  allocenv6(S1)
    implicit none
    type (env_8),INTENT(INout)::S1(ndd)

    call alloc(s1,ndd)

  END SUBROUTINE allocenv6

  SUBROUTINE  killenv6(S1)

    implicit none
    type (env_8),INTENT(INout)::S1(ndd)

    call kill(s1,ndd)

  END SUBROUTINE killenv6






  SUBROUTINE  REAL6real_8(S2,S1)
    implicit none
    real(dp),INTENT(inOUT)::S2(ndd)
    type (real_8),INTENT(IN)::S1(ndd)
    integer i


    do i=1,ndd
       s2(i)=s1(i)          !%t
    enddo
  END SUBROUTINE REAL6real_8

  SUBROUTINE  real_8REAL_8(S1,S2)
    implicit none
    type (real_8),INTENT(in)::S2(ndd)
    type (real_8),INTENT(inOUT)::S1(ndd)
    integer i


    do i=1,ndd
       s1(i)=s2(i)
    enddo
  END SUBROUTINE real_8REAL_8

  SUBROUTINE  real_8REAL6(S1,S2)
    implicit none
    real(dp),INTENT(in)::S2(ndd)
    type (real_8),INTENT(inOUT)::S1(ndd)
    integer i


    do i=1,ndd
       s1(i)=s2(i)
    enddo
  END SUBROUTINE real_8REAL6


  SUBROUTINE  REAL6env_8(S2,S1)
    implicit none
    real(dp),INTENT(inout)::S2(ndd)
    type (env_8),INTENT(in)::S1(ndd)
    integer i


    do i=1,ndd
       s2(i)=s1(i)%v
    enddo
  END SUBROUTINE REAL6env_8

  SUBROUTINE  env_8map(S1,S2)
    implicit none
    type (damap),INTENT(in)::S2
    type (env_8),INTENT(inOUT)::S1(ndd)
    integer i


    do i=1,ndd
       s1(i)%v=s2%v(i)
    enddo
  END SUBROUTINE env_8map

  SUBROUTINE  env_8benv(S1,S2)
    implicit none
    type (beamenvelope),INTENT(in)::S2
    type (env_8),INTENT(inOUT)::S1(ndd)

    s1=s2%sij0

  END SUBROUTINE env_8benv

  SUBROUTINE  env_8t(S1,S2)
    implicit none
    type (taylor),INTENT(in)::S2
    type (env_8),INTENT(inOUT)::S1(ndd)
    integer i,j
    character(6) ind,ind0
    real(dp) fac

    ind0='000000'
    do i=1,ndd
       do j=1,ndd
          ind  = ind0
          if(i==j) then
             ind(i:i)='2'
             fac=one
          else
             fac=half
             ind(i:i)='1'
             ind(j:j)='1'
          endif
          s1(i)%sigma0(j)=(s2.par.ind)*fac
          !        s1(i)%e(j)=(s2.par.ind)*fac
       enddo
    enddo
  END SUBROUTINE env_8t

  SUBROUTINE  tenv_8(S2,S1)
    implicit none
    type (taylor),INTENT(inout)::S2
    type (env_8),INTENT(in)::S1(ndd)
    type(damap) id
    type(real_8) x(6),s
    integer i,j

    call alloc(id)
    call alloc(x,6)
    call alloc(s)
    id=1
    x=id
    do i=1,ndd
       do j=1,ndd
          s=s1(i)%sigmaf(j)*x(i)*x(j)+s
       enddo
    enddo
    s2=s
    call kill(id)
    call kill(s)
    call kill(x,6)
  END SUBROUTINE tenv_8

  ! End of Some polymorphism

  FUNCTION daddsc( S1, S2 )
    implicit none
    TYPE (real_8) daddsc(ndd)
    TYPE (damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    integer localmaster,iia(4),ico(4),nd2,i
    call liepeek(iia,ico)
    nd2=iia(4)


    do i=1,nd2
       localmaster=master
       call ass(daddsc(i))
       daddsc(i)=s1%v(i)+s2(i)
       master=localmaster
    enddo
    do i=nd2+1,ndd
       localmaster=master
       call ass(daddsc(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5) then
          daddsc(i)=s2(i)+(one.mono.'00001')
       else
          daddsc(i)=s2(i)
       endif
       master=localmaster
    enddo

  END FUNCTION daddsc

  FUNCTION scdadd( S2,S1  )
    implicit none
    TYPE (real_8) scdadd(ndd)
    TYPE (damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    integer localmaster,iia(4),ico(4),nd2,i
    call liepeek(iia,ico)
    nd2=iia(4)

    do i=1,nd2
       localmaster=master
       call ass(scdadd(i))
       scdadd(i)=s1%v(i)+s2(i)
       master=localmaster
    enddo
    do i=nd2+1,ndd
       localmaster=master
       call ass(scdadd(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5) then
          scdadd(i)=s2(i)+(one.mono.'00001')
       else
          scdadd(i)=s2(i)
       endif
       master=localmaster
    enddo


  END FUNCTION scdadd

end module S_extend_poly



! MODULE ANBN SOLVES MAXWELL'S EQUATION FOR IDEAL SECTOR BEND OF TYPE TEAPOT

Module anbn
  use precision_constants
  implicit none
  public
  !frs real(dp),private::eps_extend_poly=c_1e_10

  TYPE B_CYL
     integer firsttime
     integer, POINTER ::  nmul,n_mono
     integer, DIMENSION(:), POINTER   :: i,j
     real(dp), DIMENSION(:,:), POINTER   :: a_x,a_y,b_x,b_y
  END  TYPE B_CYL


contains

  subroutine print_b(b,nmul,in,mf)
    implicit none
    type(B_CYL) b
    integer i,nmul,mf,in


    do i=1,b%n_mono
       if(in==1.and.abs(b%a_x(nmul,i))>eps_extend_poly)write(mf,*) b%i(i),b%j(i),b%a_x(nmul,i)
       if(in==2.and.abs(b%a_y(nmul,i))>eps_extend_poly)write(mf,*) b%i(i),b%j(i),b%a_y(nmul,i)
       if(in==3.and.abs(b%b_x(nmul,i))>eps_extend_poly)write(mf,*) b%i(i),b%j(i),b%b_x(nmul,i)
       if(in==4.and.abs(b%b_y(nmul,i))>eps_extend_poly)write(mf,*) b%i(i),b%j(i),b%b_y(nmul,i)
    enddo

  end subroutine   print_b

  subroutine make_coef(b,no)
    implicit none
    integer no
    integer i,j,k,a,m , ic

    type(B_CYL) b

    ic=1
    b%firsttime=-100
    allocate(b%nmul)
    allocate(b%n_mono)
    b%nmul=no
    b%n_mono=((no+2-ic)*(no+1-ic))/2
    allocate(b%i(b%n_mono),b%j(b%n_mono))
    allocate(b%a_x(no,b%n_mono),b%a_y(no,b%n_mono))
    allocate(b%b_x(no,b%n_mono),b%b_y(no,b%n_mono))

    do i=1,no
       do j=1,b%n_mono
          b%a_x(i,j)=zero
          b%a_y(i,j)=zero
          b%b_x(i,j)=zero
          b%b_y(i,j)=zero
       enddo
    enddo


    k=0
    m=no-1
    do a=m,1,-1
       do j=m-a,1,-1
          k=k+1
          b%i(k)=a
          b%j(k)=j
       enddo
       k=k+1
       b%i(k)=a
       b%j(k)=0
    enddo
    do j=m,1,-1
       k=k+1
       b%i(k)=0
       b%j(k)=j
    enddo
    k=k+1
    b%i(k)=0
    b%j(k)=0

  end subroutine make_coef

  subroutine nul_coef(b)
    implicit none
    type(B_CYL) b
    if(b%firsttime/=-100) then
       nullify(b%nmul)
       nullify(b%n_mono)
       nullify(b%i)
       nullify(b%j)
       nullify(b%a_x)
       nullify(b%a_y)
       nullify(b%b_x)
       nullify(b%b_y)
       b%firsttime=-100
    else
       deallocate(b%nmul)
       deallocate(b%n_mono)
       deallocate(b%i)
       deallocate(b%j)
       deallocate(b%a_x)
       deallocate(b%a_y)
       deallocate(b%b_x)
       deallocate(b%b_y)

    endif


  end subroutine nul_coef

  subroutine curvebend(b_sol,no)
    use polymorphic_complextaylor
    implicit none
    LOGICAL(lp) :: doneitt=.TRUE.
    integer result,error
    integer no,i,j,jd(lnv)
    type (damap) y0
    type (taylor) t
    type (taylorresonance) tr
    type (complextaylor) z,zb
    type (complextaylor), ALLOCATABLE,dimension(:)::fs
    type (complextaylor), ALLOCATABLE,dimension(:,:)::F
    type(taylor), ALLOCATABLE,dimension(:)::a,b
    type(taylor), ALLOCATABLE,dimension(:,:)::da,db
    real(dp), ALLOCATABLE,dimension(:,:)::MA,MB
    complex(dp) z0
    type(b_cyl) b_sol


    allocate(MA(no,no))
    allocate(MB(no,no))
    allocate(F(no,no))
    allocate(Fs(no))
    allocate(a(no))
    allocate(b(no))
    allocate(da(no,2))
    allocate(db(no,2))

    call init(no,1,0,0,doneitt)
    !call daeps(c_1e_10)

    call alloc(y0)
    call alloc(t)
    call alloc(tr)
    call alloc(z)
    call alloc(zb)
    call alloc(fs,no)
    call alloc(a,no)
    call alloc(b,no)
    do i=1,no
       call alloc(da(i,1))
       call alloc(da(i,2))
       call alloc(db(i,1))
       call alloc(db(i,2))
    enddo
    do i=1,no
       do j=1,no
          ma(i,j)=zero
          mb(i,j)=zero
          call alloc(F(i,j))
       enddo
    enddo



    y0=1
    y0%v(2)=zero

    z0=zero
    !    call var(z,z0,1,2)
    z=z0.var.(/1,2/)
    !z%i=-z%i




    do i=1,no

       F(i,i)=-(z**i)  !/i    ! Harmonic guess

       fs(i)=fs(i)+f(i,i)

       do j=i,no-1
          call s_op(f(i,j),f(i,j+1))
          fs(i)=fs(i)+f(i,j+1)
       enddo

       a(i)=aimag(fs(i))
       b(i)=-dble(fs(i))
       a(i)=(a(i).d.2)/(one+(one.mono.'1'))
       b(i)=(b(i).d.1)/(one+(one.mono.'1'))
    enddo


    jd(:)=0
    do j=1,no

       a(j)=a(j).o.y0
       b(j)=b(j).o.y0

       do i=1,no
          jd(1)=i-1
          call pek(a(j),jd,ma(i,j))
       enddo


       do i=1,no
          jd(1)=i-1
          call pek(b(j),jd,mb(i,j))
       enddo

    enddo

    call matinv(MA,MA,no,no,result)
    if(result/=0) then
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72)))'
       w_p%c(1)= " failed inversion in curvebend part 1"
       CALL WRITE_E(result)
    endif
    call matinv(MB,MB,no,no,result)
    if(result/=0) then
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,A72)))'
       w_p%c(1)= " failed inversion in curvebend part 2"
       CALL WRITE_E(result)
    endif

    do i=1,no
       a(i)=zero
       b(i)=zero
    enddo


    do i=1,no
       do j=1,no
          a(i)=MA(j,i)*aimag(fs(j))+a(i)
          b(i)=MB(j,i)*dble(fs(j))+b(i)
       enddo
       da(i,1)=a(i).d.1
       da(i,2)=a(i).d.2
       db(i,1)=b(i).d.1
       db(i,2)=b(i).d.2
    enddo

    do i=1,no
       do j=1,b_sol%n_mono
          jd(1)=b_sol%i(j)
          jd(2)=b_sol%j(j)
          if(jd(1)+jd(2)>no-1) then
             w_p=1
             w_p%nc=1
             w_p%fc='(1((1X,a20)))'
             w_p=(/jd(1),jd(2)/);w_p%fi= '(2((1X,i4)),/)'
             w_p%c(1) = " curvebend error 1"
             CALL WRITE_E(1)
          endif
          call pek(da(i,1),jd,b_sol%a_x(i,j))
          if(abs(b_sol%a_x(i,j))<eps_extend_poly) b_sol%a_x(i,j)=zero
       enddo
       do j=1,b_sol%n_mono
          jd(1)=b_sol%i(j)
          jd(2)=b_sol%j(j)
          if(jd(1)+jd(2)>no-1) then
             w_p=1
             w_p%nc=1
             w_p%fc='(1((1X,a20)))'
             w_p=(/jd(1),jd(2)/);w_p%fi= '(2((1X,i4)),/)'
             w_p%c(1) = " curvebend error 2"
             CALL WRITE_E(2)
          endif
          call pek(da(i,2),jd,b_sol%a_y(i,j))
          if(abs(b_sol%a_y(i,j))<eps_extend_poly) b_sol%a_y(i,j)=zero
       enddo
       do j=1,b_sol%n_mono
          jd(1)=b_sol%i(j)
          jd(2)=b_sol%j(j)
          if(jd(1)+jd(2)>no-1) then
             w_p=1
             w_p%nc=1
             w_p%fc='(1((1X,a20)))'
             w_p=(/jd(1),jd(2)/);w_p%fi= '(2((1X,i4)),/)'
             w_p%c(1) = " curvebend error 3"
             CALL WRITE_E(3)
          endif
          call pek(db(i,1),jd,b_sol%b_x(i,j))
          if(abs(b_sol%b_x(i,j))<eps_extend_poly) b_sol%b_x(i,j)=zero
       enddo
       do j=1,b_sol%n_mono
          jd(1)=b_sol%i(j)
          jd(2)=b_sol%j(j)
          if(jd(1)+jd(2)>no-1) then
             w_p=1
             w_p%nc=1
             w_p%fc='(1((1X,a20)))'
             w_p=(/jd(1),jd(2)/);w_p%fi= '(2((1X,i4)),/)'
             w_p%c(1) = " curvebend error 4"
             CALL WRITE_E(4)
          endif
          call pek(db(i,2),jd,b_sol%b_y(i,j))
          if(abs(b_sol%b_y(i,j))<eps_extend_poly) b_sol%b_y(i,j)=zero
       enddo
    enddo

    call kill(y0)
    call kill(t)
    call kill(tr)
    call kill(z)
    call kill(zb)
    call kill(fs,no)
    call kill(a,no)
    call kill(b,no)
    do i=1,no
       call kill(da(i,1))
       call kill(da(i,2))
       call kill(db(i,1))
       call kill(db(i,2))
    enddo
    do i=1,no
       do j=1,no
          call kill(F(i,j))
       enddo
    enddo

    IF (ALLOCATED(MA)) THEN
       DEALLOCATE (MA, STAT = error)
       IF(ERROR==0) THEN
          ! WRITE(6,*) " MA ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) = " MA ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(100)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) = " MA ARRAY not DEALLOCATED : PROBLEMS"
       CALL WRITE_E(101)
    ENDIF
    IF (ALLOCATED(MB)) THEN
       DEALLOCATE (MB, STAT = error)
       IF(ERROR==0) THEN
          ! WRITE(6,*) " MB ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" MB ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(102)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) ="MB DID NOT EXIST (curvebend)"
       CALL WRITE_E(103)
    ENDIF
    IF (ALLOCATED(F)) THEN
       DEALLOCATE (F, STAT = error)
       IF(ERROR==0) THEN
          !  WRITE(6,*) " F ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" F ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(104)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) =" F DID NOT EXIST (curvebend)"
       CALL WRITE_E(105)
    ENDIF
    IF (ALLOCATED(FS)) THEN
       DEALLOCATE (FS, STAT = error)
       IF(ERROR==0) THEN
          !  WRITE(6,*) " FS ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" FS ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(106)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) =" FS DID NOT EXIST (curvebend)"
       CALL WRITE_E(106)
    ENDIF
    IF (ALLOCATED(A)) THEN
       DEALLOCATE (A, STAT = error)
       IF(ERROR==0) THEN
          !  WRITE(6,*) " A ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" A ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(107)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) =" A DID NOT EXIST (curvebend)"
       CALL WRITE_E(107)
    ENDIF
    IF (ALLOCATED(B)) THEN
       DEALLOCATE (B, STAT = error)
       IF(ERROR==0) THEN
          !  WRITE(6,*) " B ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" B ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(108)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) =" B DID NOT EXIST (curvebend)"
       CALL WRITE_E(109)
    ENDIF
    IF (ALLOCATED(DA)) THEN
       DEALLOCATE (DA, STAT = error)
       IF(ERROR==0) THEN
          !  WRITE(6,*) " DA ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" DA ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(110)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) =" DA DID NOT EXIST (curvebend)"
       CALL WRITE_E(111)
    ENDIF
    IF (ALLOCATED(DB)) THEN
       DEALLOCATE (DB, STAT = error)
       IF(ERROR==0) THEN
          !  WRITE(6,*) " DB ARRAY DEALLOCATED "
       ELSE
          w_p=0
          w_p%nc=1
          w_p%fc='(1((1X,a72)))'
          w_p%c(1) =" DB ARRAY not DEALLOCATED : PROBLEMS"
          CALL WRITE_E(111)
       ENDIF
    ELSE
       w_p=0
       w_p%nc=1
       w_p%fc='(1((1X,a72)))'
       w_p%c(1) =" DB DID NOT EXIST (curvebend)"
       CALL WRITE_E(111)
    ENDIF

  end subroutine  curvebend





  subroutine s_op(f,g)
    use polymorphic_complextaylor
    implicit none
    type (complextaylor) f,g
    type (taylor) t(2)
    type (taylorresonance) tr(2)

    call alloc(t,2)
    call alloc(tr(1))
    call alloc(tr(2))

    !    g=( ((f.d.1).d.1) +((f.d.2).d.2 ))
    !    g=(f.d.1)-(one.mono.'1')*g

    !    work around a possible compiler bug on the SUN => no logical explanation!
    g%r=((f%r.d.1).d.1) + ((f%r.d.2).d.2)
    g%i=((f%i.d.1).d.1) + ((f%i.d.2).d.2)
    g%r=(f%r.d.1)-(one.mono.'1')*g%r
    g%i=(f%i.d.1)-(one.mono.'1')*g%i

    t(1)=g%r
    t(2)=g%i
    tr(1)=t(1)
    tr(2)=t(2)

    call cfu(tr(1)%cos,op,tr(1)%cos )
    call cfu(tr(2)%cos,op,tr(2)%cos )
    call cfu(tr(1)%sin,op,tr(1)%sin )
    call cfu(tr(2)%sin,op,tr(2)%sin )

    tr(1)%cos=tr(1)%cos*(one.mono.'10')*(one.mono.'01')/four
    tr(1)%sin=tr(1)%sin*(one.mono.'10')*(one.mono.'01')/four
    tr(2)%cos=tr(2)%cos*(one.mono.'10')*(one.mono.'01')/four
    tr(2)%sin=tr(2)%sin*(one.mono.'10')*(one.mono.'01')/four



    t(1)=tr(1)
    t(2)=tr(2)

    g%r=t(1)
    g%i=t(2)

    call kill(t,2)
    call kill(tr(1))
    call kill(tr(2))


  end subroutine s_op

  real(dp) function op(J)
    implicit none
    integer,dimension(:)::J

    op=one/(REAL(j(1),kind=DP)+one)/(REAL(j(2),kind=DP)+one)

  end  function op

end Module anbn
