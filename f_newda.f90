!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size
module newda
  use define_newda
  USE DABNEW
  use lielib_berz
  implicit none
  PRIVATE
  !private location,expon,deldol,del2dol,linear
  !private no,nv, LNV,EPSdolmac,EPSdol,flag !,numnewda ,nm locs,loce,
  INTEGER, ALLOCATABLE,dimension(:)::pow
  INTEGER, ALLOCATABLE,dimension(:,:)::expon
  INTEGER, ALLOCATABLE,dimension(:,:,:)::location
  !  INTEGER, ALLOCATABLE,dimension(:,:)::linear
  integer nm,flag(3),numnewda
  logical(lp) table(totaldol)
  logical(lp)  non1
  type (taylorlow) deldol,del2dol,IPOW,INON,ISCR
  CHARACTER(120) LINE

  integer,PARAMETER::LNV=100
  INTEGER LOCS(-1:4),LOCE(-1:4)
  PUBLIC numnewda,TABLE,locs,loce
  PUBLIC INITIALIZE_DA,NEWETALL,NEWDAVAR,NEWDACOP,NEWDACCT,NEWDADAL
  PUBLIC NEWDACON,NEWDAPEK,NEWDADER,NEWDAMUL,NEWDAADD,NEWDAINV
  PUBLIC NEWDAPIN,NEWDAPOK,NEWDACFU,NEWDASUB,NEWDACLR,NEWDACMU,NEWDALIN
  PUBLIC NEWDAREA,OLDDAREA,NEWDAPRI,OLDDAPRI,NEWDAEPS,NEWDAABS,NEWDATRA
  PUBLIC ALLOCNEWDA,KILLNEWDAS,NEWDAFUN,NEWDADIV,NEWDAPOI,NEWDADIC,NEWDACDI,NEWDACAD
  PUBLIC NEWDACSU,NEWDASUC,NEWDASHIFT,NEWDARAN,NEWDACFUR,NEWDACFUI,OLDDAPRI77,OLDDAREA77
  PUBLIC NEWDANCD,NEWPPUSH1,DE_INITIALIZE_DA,nullnewda


  INTERFACE NEWDADAL
     MODULE PROCEDURE NEWDADALn
     MODULE PROCEDURE NEWDADAL0
  END INTERFACE

  INTERFACE NEWETALL
     MODULE PROCEDURE NEWETALLN
     MODULE PROCEDURE NEWETALL0
  END INTERFACE

  INTERFACE NEWDACCT
     MODULE PROCEDURE NEWDACCTN
     MODULE PROCEDURE NEWDACCT0
  END INTERFACE

contains

  subroutine DE_initialize_da
    implicit none
    integer error,IPAUSE,MYPAUSES
    call KILLNEWDAS(DELdol)
    call KILLNEWDAS(DEL2dol)
    CALL KILLNEWDAS(IPOW)
    CALL KILLNEWDAS(INON)
    CALL KILLNEWDAS(ISCR)

    IF (ALLOCATED(location)) THEN
       DEALLOCATE (location, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          LINE= " LOCATION ARRAY not DEALLOCATED : PROBLEMS"
          IPAUSE=MYPAUSES(-1,LINE)
       ENDIF
    ELSE
       if(warnda) THEN
          LINE=  "location WAS NOT  ALLOCATED "
          IPAUSE=MYPAUSES(1,LINE)
       ENDIF
    ENDIF

    IF (ALLOCATED(expon)) THEN

       DEALLOCATE (expon, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          LINE= " EXPON ARRAY not DEALLOCATED : PROBLEMS"
          IPAUSE=MYPAUSES(-2,LINE)
       ENDIF

    ELSE
       if(warnda) THEN
          LINE=  "expon WAS NOT ALLOCATED "
          IPAUSE=MYPAUSES(2,LINE)
       ENDIF
    ENDIF

    IF (ALLOCATED(pow)) THEN

       DEALLOCATE (pow, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          LINE= " pow ARRAY not DEALLOCATED : PROBLEMS"
          IPAUSE=MYPAUSES(-2,LINE)
       ENDIF

    ELSE
       if(warnda.AND.PACKING) THEN
          LINE=  " pow NOT ALLOCATED "
          IPAUSE=MYPAUSES(3,LINE)
       ENDIF
    ENDIF



  END subroutine DE_initialize_da

  subroutine NEWcheck(A,k)
    implicit none
    integer ipause, mypauses,MYPAUSE
    type (taylorlow) A
    integer i,J ,k

    return


    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " NEWcheck"
          IPAUSE=MYPAUSES(-997,LINE)
       ENDIF
    endif
    if(.not.packing) return
    J=0
    do i=1,nm
       if(a%yes(i)) then
          j=j+1
       endif
    enddo
    if(a%m/=j) then
       w_p=0
       w_p%nc=2
       write(w_p%c(1),'(a30,i8)') "problems in  NEWcheck, a%id = ",a%id
       write(w_p%c(2),'(a9,3(1x,i8))') "a%m ,j,k ",  a%m ,j,k
       w_p%fc='((1X,A72,/),(1X,A72))'
       CALL WRITE_I
       w_p=0
       w_p%nc=1
       do i=1,nm
          if(a%yes(i)) then
             write(w_p%c(1),'(2(1x,i8))') i,a%r(i)
             CALL WRITE_i
          endif
       enddo
       ipause=mypause(9999)
    endif

  end subroutine NEWcheck


  subroutine initialize_da(no1,nv1)
    implicit none
    integer no1,nv1 ,i,j,k,ic,ik,ij, nof(3),ind(3)

    no=no1
    nv=nv1
    nm=0
    if(no==1) then
       non1=.false.
    else
       non1=.true.
    endif
    do ik=-1,4
       LOCS(IK)=0
       LOCE(IK)=0
    enddo
    do i=1,3
       ind(i)=0
       nof(i)=0
       flag(i)=0
    enddo
    do i=1,no
       flag(i)=1
       nof(i)=nv
    enddo

    allocate(location(0:nof(1),0:nof(2),0:nof(3)))

    do i=0,nof(1)
       do j=0,nof(2)
          do k=0,nof(3)
             location(i,j,k)=0
          enddo
       enddo
    enddo
    do ik=0,no  !3
       LOCS(IK)=NM+1
       LOCE(IK-1)=NM
       do i=0,nv  !nof(1)
          do j=0,i*flag(2)
             do k=0,j*flag(3)

                ic=0
                if(i>0) ic=ic+1
                if(j>0) ic=ic+1
                if(k>0) ic=ic+1
                if(ic.eq.ik) then
                   nm=nm+1
                   location(i,j,k)=nm
                endif
             enddo
          enddo
       enddo
    enddo
    LOCE(NO)=NM


    allocate(expon(no,nm))
    do i=1,no
       do j=1,nm
          expon(i,j)=0
       enddo
    enddo

    do ik=0,no  !3
       do i=0,nv
          do j=0,i*flag(2)
             do k=0,j*flag(3)
                ind(1)=i
                ind(2)=j
                ind(3)=k
                ic=0
                if(i>0) ic=ic+1
                if(j>0) ic=ic+1
                if(k>0) ic=ic+1
                if(ic.eq.ik) then
                   do ij=1,no  !3
                      expon(ij,location(i,j,k))=ind(ij)
                   enddo
                endif
             enddo
          enddo
       enddo
    enddo

    do i=0,nv
       do j=0,i*flag(2)
          do k=0,j*flag(3)
             if(no>1) location(j*flag(1),i*flag(2),k*flag(3))=location(i,j,k)
             if(no>2) then
                location(k*flag(1),i*flag(2),j*flag(3))=location(i,j,k)
                location(k*flag(1),j*flag(2),i*flag(3))=location(i,j,k)
                location(j*flag(1),k*flag(2),i*flag(3))=location(i,j,k)
                location(i*flag(1),k*flag(2),j*flag(3))=location(i,j,k)
             endif
          enddo
       enddo
    enddo

    if(packing) then
       allocate(pow(nm))
       DO I=0,NO
          DO J=LOCS(I),LOCE(I)
             pow(j)=i
          enddo
       enddo
    endif





    call allocnewda(deldol)
    call allocnewda(del2dol)
    CALL allocnewda(IPOW)
    CALL allocnewda(INON)
    CALL allocnewda(ISCR)

  end subroutine initialize_da

  subroutine nullnewda(A)
    implicit none
    type (taylorlow) A
    nullify(A%R)
    nullify(A%id)
    nullify(A%NZ)
    nullify(A%yes)
    nullify(A%m)
  end subroutine nullnewda

  subroutine allocnewda(A)
    implicit none
    integer ind
    !logical(lp) incnda
    type (taylorlow) A


    if(newdadynamical) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          ind=1

20        CONTINUE
          IF(TABLE(IND)) THEN
             IND=IND+1
             GOTO 20
          ENDIF
          TABLE(IND)=.TRUE.

          allocate(A%R(NM))
          allocate(A%ID)
          a%id=ind
          numnewda=numnewda+1
          if(packing) then
             allocate(A%nz(NM))
             allocate(A%yes(NM))
             allocate(A%M)
             do ind=1,nm
                A%r(ind)=zero
                A%yes(ind)=.false.
                A%nz(ind)=zero
                A%r(ind)=zero
             enddo
             a%m=0
          else

             call newdaclr(a)

          endif
       else
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(A16,I8)')  " already exists " , a%id
          CALL WRITE_E(-998)
       ENDIF

    endif

  end subroutine allocnewda





  subroutine newdadalN(A,n)
    implicit none
    type (taylorlow),dimension(:)::A
    integer n,i
    do i=1,n
       call KILLNEWDAS(A(i))
    enddo
  end subroutine newdadalN

  subroutine newdadal0(A,n)
    implicit none
    integer n
    type (taylorlow) A
    call KILLNEWDAS(A)
  end subroutine newdadal0

  ! allocates vector of n polynomials
  SUBROUTINE newETALLN(X,N)
    IMPLICIT none
    INTEGER I1(4),I2(4),i,n
    type (taylorlow),dimension(:)::x
    DO  I=1,iabs(n)
       call allocnewda(x(i))
    enddo
    if(N.lt.0) then
       CALL LIEPEEK(I1,I2)
       ND2=I1(4)
       do  i=nd2+1,-n
          call newdavar(x(i),zero,i)
       ENDDO
    endif

  END SUBROUTINE newETALLN
  ! allocates vector of n polynomials
  SUBROUTINE newETALL0(X,N)
    IMPLICIT none
    INTEGER n
    type (taylorlow) x
    call allocnewda(x)
  END SUBROUTINE newETALL0


  subroutine KILLNEWDAS(A)
    implicit none
    integer ipause, mypauseS,MYPAUSE
    integer ind
    type (taylorlow) A

    IF(newdadynamical) THEN
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE=" DID NOT EXIST (KILLNEWDAS)"
          ipause=mypauseS(535,LINE)
       else
          ind=a%id
          if(ind>0.and.ind<=totaldol) then
             TABLE(IND)=.FALSE.
          else
             WRITE(LINE,'(A6,I8)') " a%id ",  a%id
             ipause=mypauseS(345,LINE)
          endif
          numnewda=numnewda-1
          a%id=0
          IF (ASSOCIATED(A%R)) then
             DEallocate(A%R)
          else
             ipause=mypause(1231)
          endif
          IF (ASSOCIATED(A%id)) then
             DEallocate(A%ID)
          else
             ipause=mypause(1232)
          endif
          if(packing) then
             a%m=0
             IF (ASSOCIATED(A%nz)) then
                DEallocate(A%nz)
             else
                ipause=mypause(1233)
             endif
             IF (ASSOCIATED(A%yes)) then
                DEallocate(A%yes)
             else
                ipause=mypause(1234)
             endif
             IF (ASSOCIATED(A%m)) then
                DEallocate(A%m)
             else
                ipause=mypause(1235)
             endif
          endif
          if (ASSOCIATED(A%R)) then
             WRITE(LINE,'(A36,I8)') "error ind= ",ind
             ipause=mypauseS(0,LINE)
          endif
       ENDIF
    ENDIF

  end subroutine KILLNEWDAS



  subroutine newdacmu(A,RA,C)
    !     DACMU(A,RA,C):      PERFORMS C = A * RA
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,C
    integer i
    real(dp) RA
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " NEWDACMU"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(C%r)) THEN
          LINE= " NEWDACMU"
          IPAUSE=MYPAUSES(-999,LINE)
       ENDIF
    endif

    IF(PACKING)  THEN
       ! call newcheck(a,-2)
       ! call newcheck(c,-21)

       if(c%id/=a%id) then
          call newdaclr(c)
          c%m=a%m
          do i=1,c%m
             C%nz(I)=A%nz(I)
             C%yes(A%nz(I))=A%yes(A%nz(I))
             C%R(A%nz(I))=RA*A%R(A%nz(I))
          enddo

       else
          do i=1,a%m
             C%R(A%nz(I))=RA*A%R(A%nz(I))
          enddo

       endif
       ! call newcheck(c,-3)
       ! call newcheck(a,-20)

    else

       do i=1,nm
          C%R(I)=RA*A%R(I)
       enddo

    endif



  end subroutine newdacmu

  subroutine newDACDI(A,RA,C)
    !     DACDI(A,RA,C):      PERFORMS C = A / RA
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,C
    integer i
    real(dp) RA
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " newDACDI"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(C%r)) THEN
          LINE= " newDACDI"
          IPAUSE=MYPAUSES(-999,LINE)
       ENDIF
    endif


    IF(PACKING)  THEN
       ! call newcheck(a,-4)
       ! call newcheck(c,-40)
       if(c%id/=a%id) then
          call newdaclr(c)
          c%m=a%m
          do i=1,c%m
             C%nz(I)=A%nz(I)
             C%yes(A%nz(I))=A%yes(A%nz(I))
             C%R(A%nz(I))=A%R(A%nz(I))/ra
          enddo
       else
          do i=1,a%m
             C%R(A%nz(I))=A%R(A%nz(I))/ra
          enddo
       endif
       ! call newcheck(a,-40)
       ! call newcheck(c,-5)
    else


       do i=1,nm
          C%R(I)=A%R(I)/ra
       enddo


    ENDIF



  end subroutine newDACDI



  subroutine newDACSU(A,RA,C)
    !     DACSU(A,RA,C):      PERFORMS C = A - RA
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,C
    integer i
    real(dp) RA
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " newDACSU"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          LINE= " newDACSU"
          IPAUSE=MYPAUSES(-999,LINE)
       ENDIF
    endif

    IF(PACKING)  THEN
       ! call newcheck(a,-6)
       ! call newcheck(c,-66)
       if(c%id/=a%id) then  !!
          call newdaclr(c)
          c%m=a%m
          do i=1,c%m
             C%nz(I)=A%nz(I)
             C%r(A%nz(I))=A%r(A%nz(I))
             C%yes(A%nz(I))=A%yes(A%nz(I))
          enddo
          if(a%yes(1)) then
             C%R(1)=C%R(1)-ra
          else
             c%m=c%m+1
             c%yes(1)=.true.
             C%nz(c%m)=1
             C%R(1)=-ra
          endif
       else  !!
          if(a%yes(1)) then
             C%R(1)=C%R(1)-ra
          else
             c%m=c%m+1
             c%yes(1)=.true.
             C%nz(c%m)=1
             C%R(1)=-ra
          endif
       endif !!
       ! call newcheck(c,-7)
       ! call newcheck(a,-77)

    else

       !if(a%id/=c%id) then
       do i=1,nm
          C%R(I)=A%R(I)
       enddo
       !endif
       C%R(1)=C%R(1)-ra

    ENDIF


  end subroutine newDACSU

  subroutine newDACAD(A,RA,C)
    !     DACAD(A,RA,C):      PERFORMS C = A + RA
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,C
    integer i
    real(dp) RA
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " newDACAD"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          LINE= " newDACAD"
          IPAUSE=MYPAUSES(-999,LINE)
       ENDIF
    endif



    IF(PACKING)  THEN
       ! call newcheck(a,-8)
       ! call newcheck(c,-888)
       if(c%id/=a%id) then   !!!
          call newdaclr(c)
          c%m=a%m
          do i=1,c%m
             C%nz(I)=A%nz(I)
             C%r(A%nz(I))=A%r(A%nz(I))
             C%yes(A%nz(I))=A%yes(A%nz(I))
          enddo

          if(a%yes(1)) then
             C%R(1)=C%R(1)+ra
          else
             c%m=c%m+1
             c%yes(1)=.true.
             C%nz(c%m)=1
             C%R(1)=ra
          endif
       else !!!
          if(a%yes(1)) then
             C%R(1)=C%R(1)+ra
          else
             c%m=c%m+1
             c%yes(1)=.true.
             C%nz(c%m)=1
             C%R(1)=ra
          endif
       endif      !!!
       ! call newcheck(a,-99)
       ! call newcheck(c,-9)

    else

       !if(a%id/=c%id) then
       do i=1,nm
          C%R(I)=A%R(I)
       enddo
       !endif
       C%R(1)=C%R(1)+ra

    ENDIF

  end subroutine newDACAD

  subroutine NEWDASUC(A,RA,C)
    !     DASUC(A,RA,C):      PERFORMS C = RA - A
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,C
    integer i
    real(dp) RA
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " NEWDASUC"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          LINE= " NEWDASUC"
          IPAUSE=MYPAUSES(-999,LINE)
       ENDIF
    endif

    IF(PACKING)  THEN
       ! call newcheck(a,-101)
       ! call newcheck(c,-10)
       if(c%id/=a%id) then    !!
          call newdaclr(c)
          c%m=a%m
          do i=1,c%m
             C%nz(I)=A%nz(I)
             C%r(A%nz(I))=-A%r(A%nz(I))
             C%yes(A%nz(I))=A%yes(A%nz(I))
          enddo
          if(a%yes(1)) then
             C%R(1)=RA-A%R(1)
          else
             c%m=c%m+1
             c%yes(1)=.true.
             C%nz(c%m)=1
             C%R(1)=ra
          endif
       else  !!!
          if(a%yes(1)) then
             C%R(1)=RA-A%R(1)
          else
             c%m=c%m+1
             c%yes(1)=.true.
             C%nz(c%m)=1
             C%R(1)=ra
          endif
       endif

       ! call newcheck(c,-11)
       ! call newcheck(a,-111)

    else

       C%R(1)=RA-A%R(1)
       do i=2,nm
          C%R(I)=-A%R(I)
       enddo

    ENDIF





  end subroutine NEWDASUC

  subroutine NEWDACMA(A,B,RB,C)
    !     DACMA(A,B,RB,C):    PERFORMS C = A + RB*B
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,B,C
    integer i
    real(dp) RB
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " NEWDACMA"
          IPAUSE=MYPAUSES(-997,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          LINE= " NEWDACMA"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          LINE= " NEWDACMA"
          IPAUSE=MYPAUSES(-999,LINE)
       ENDIF
    endif

    IF(PACKING)  THEN
       ! call newcheck(a,-12)
       ! call newcheck(b,-13)
       ! call newcheck(c,-14)

       call newdaclr(deldol)
       ! call newcheck(deldol,-1411)
       do i=1,A%m
          deldol%r(A%nz(I))=a%r(A%nz(I))
          deldol%yes(A%nz(I))=.true.
          deldol%nz(i)=A%nz(I)
       enddo
       deldol%m=A%m
       do i=1,b%m
          deldol%r(b%nz(I))=rb*b%r(b%nz(I))+deldol%r(b%nz(I))
          if(.not.deldol%yes(b%nz(I))) then
             deldol%m=deldol%m+1
             deldol%nz(deldol%m)=b%nz(I)
             deldol%yes(b%nz(I))=.true.
          endif
       enddo

       ! call newcheck(a,-15)
       ! call newcheck(b,-16)
       ! call newcheck(c,-17)
       ! call newcheck(deldol,-18)

       call newdacop(deldol,c)
    else

       do i=1,nm
          C%R(I)=A%R(I)+RB*B%R(I)
       enddo

    ENDIF




  end subroutine NEWDACMA

  subroutine NEWDACLR(A)
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " NEWDACLR"
          IPAUSE=MYPAUSES(-997,LINE)
       ENDIF
    endif

    if(packing) then
       do i=1,a%m
          A%yes(A%nz(I))=.false.
          A%R(A%nz(I))=zero
       enddo
       A%M=0


    else

       do i=1,nm
          A%R(I)=zero
       enddo

    endif

  end subroutine NEWDACLR

  subroutine NEWDACOP(A,B)
    implicit none
    integer ipause, mypauseS
    type (taylorlow) A,B
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          LINE= " NEWDACOP"
          IPAUSE=MYPAUSES(-997,LINE)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          LINE= " NEWDACOP"
          IPAUSE=MYPAUSES(-998,LINE)
       ENDIF
    endif

    IF(PACKING)  THEN
       ! call newcheck(a,0)
       ! call newcheck(b,-1000)

       if(b%id/=a%id) then  !!
          call newdaclr(b)
          b%m=a%m
          do i=1,B%m
             B%nz(I)=A%nz(I)
             B%yes(B%nz(I))=A%yes(B%nz(I))
             B%R(B%nz(I))=A%R(B%nz(I))
          enddo
       endif
       ! call newcheck(a, 1000)
       ! call newcheck(b, 1001)
    else

       do i=1,nm
          B%R(I)=A%R(I)
       enddo

    ENDIF

  end subroutine NEWDACOP

  subroutine NEWDAPACK(A)
    implicit none
    integer ipause, mypauses
    type (taylorlow) A
    integer i,J
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAPACK"
          ipause=mypauses(-997,line)
       ENDIF
    endif
    if(.not.packing) then
       line= " NOT packing"
       ipause=mypauses(-284,line)
    endif
    J=0
    do i=1,nm
       A%yes(i)=.false.
       IF(ABS(A%R(I))>PUNY) THEN
          J=J+1
          A%NZ(J)=I
          A%yes(i)=.true.
       ENDIF
    enddo
    A%M=J
  end subroutine NEWDAPACK




  subroutine NEWDAMUL(A,B,C)
    !     DAMUL(A,B,C):       PERFORMS C = A * B
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,B,C
    integer i1,I2,J1,J2,TEST,IND(3) ,itot,k
    real(dp) r
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAMUL"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDAMUL"
          ipause=mypauses(-998,line)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          line= " NEWDAMUL"
          ipause=mypauses(-999,line)
       ENDIF
    endif

    if(packing) then  !222

       ! call newcheck(a,11)
       ! call newcheck(b,12)
       ! call newcheck(c,13)
       ! call newcheck(Deldol,140)
       CALL NEWDACLR(Deldol)   !  PACKING
       ! call newcheck(Deldol,14)
       k=0
       DO I1=1,A%M
          DO I2=1,B%M
             J1=A%NZ(I1)
             J2=B%NZ(I2)
             itot=pow(j1)+pow(j2)
             r=A%R(J1)*B%R(J2)
             if(abs(r)>puny) then
                IF(itot<=NO) THEN

                   DO TEST=1,3   !no
                      IND(TEST)=0
                   ENDDO

                   DO TEST=1,pow(j1)
                      IND(TEST)=expon(TEST,J1)
                   ENDDO


                   DO TEST=pow(j1)+1,itot
                      IND(TEST)=expon(TEST-pow(j1),J2)
                   ENDDO


                   Deldol%R(LOCATION(IND(1),IND(2),IND(3)))=Deldol%R(LOCATION(IND(1),IND(2),IND(3)))+A%R(J1)*B%R(J2)
                   if(.not.deldol%yes(LOCATION(IND(1),IND(2),IND(3)))) then
                      k=k+1
                      deldol%nz(k)=LOCATION(IND(1),IND(2),IND(3))
                      deldol%yes(LOCATION(IND(1),IND(2),IND(3)))=.true.
                   endif

                ENDIF
             endif
          ENDDO
       ENDDO

       deldol%m=k

       ! call newcheck(a,1100)
       ! call newcheck(b,1200)
       ! call newcheck(c,1300)
       ! call newcheck(Deldol,15)

       CALL NEWDACOP(Deldol,C)


    else !222
       if(non1)then
          CALL NEWDACLR(Deldol)
          DO I1=0,NO
             DO I2=0,NO
                itot=I1+I2
                IF(itot<=NO)  THEN  ! 1

                   DO J1=LOCS(I1),LOCE(I1)
                      if(A%R(J1)/=zero) then
                         DO J2=LOCS(I2),LOCE(I2)
                            if(B%R(J2)/=zero) then
                               DO TEST=1,3   !no
                                  IND(TEST)=0
                               ENDDO

                               DO TEST=1,i1
                                  IND(TEST)=expon(TEST,J1)
                               ENDDO


                               DO TEST=i1+1,itot
                                  IND(TEST)=expon(TEST-I1,J2)
                               ENDDO


                               Deldol%R(LOCATION(IND(1),IND(2),IND(3)))=Deldol%R(LOCATION(IND(1),IND(2),IND(3)))+A%R(J1)*B%R(J2)
                            endif
                         ENDDO
                      endif
                   ENDDO
                ENDIF  ! 1
             ENDDO
          ENDDO

          CALL NEWDACOP(Deldol,C)


       else

          do i1=2,nm
             C%R(i1)=A%R(i1)*b%R(1)+B%R(i1)*A%R(1)
          enddo

          c%R(1)=A%R(1)*B%R(1)

       endif

    endif !22

  end subroutine NEWDAMUL

  subroutine NEWppush1(A,xi,xf)

    implicit none
    integer ipause, mypauses
    type (taylorlow) A
    integer i1,J1,TEST
    real(dp) xf,xt(0:lnv),r
    real(dp),dimension(:)::xi
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWpush1"
          ipause=mypauses(-997,line)
       ENDIF
    ENDIF

    xf=zero

    if(packing) then !22222

       do i1=1,nv
          xt(0)=one
          xt(i1)=xi(i1)
       enddo

       DO i1=1,a%m

          j1=a%nz(i1)
          if(A%R(J1)/=zero) then

             r=one

             DO TEST=1,no
                r=r*xt(expon(TEST,J1))
             ENDDO




             xf=A%R(J1)*r+xf

          endif

       ENDDO

    else  !22222


       if(non1)then
          do i1=1,nv
             xt(0)=one
             xt(i1)=xi(i1)
          enddo

          DO J1=1,nm

             if(A%R(J1)/=zero) then

                r=one

                DO TEST=1,no
                   r=r*xt(expon(TEST,J1))
                ENDDO




                xf=A%R(J1)*r+xf

             endif

          ENDDO


       else

          do i1=2,nm
             xf=xf+A%R(i1)*xi(i1-1)
          enddo

          xf=xf+A%R(i1)
       endif

    endif !22222

  end subroutine NEWppush1


  subroutine NEWDALIN(A,RA,B,RB,C)
    !     DALIN(A,RA,B,RB,C): PERFORMS C = A*RA + B*RB
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,B,C
    integer i
    real(dp) RA,RB
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDALIN"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDALIN"
          ipause=mypauses(-998,line)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          line= " NEWDALIN"
          ipause=mypauses(-999,line)
       ENDIF
    endif

    IF(PACKING)  THEN
       call newdaclr(deldol)
       do i=1,A%m
          deldol%r(A%nz(I))=ra*a%r(A%nz(I))
          deldol%yes(A%nz(I))=.true.
          deldol%nz(i)=A%nz(I)
       enddo
       deldol%m=A%m
       do i=1,b%m
          deldol%r(b%nz(I))=rb*b%r(b%nz(I))+deldol%r(b%nz(I))
          if(.not.deldol%yes(b%nz(I))) then
             deldol%m=deldol%m+1
             deldol%nz(deldol%m)=b%nz(I)
             deldol%yes(b%nz(I))=.true.
          endif
       enddo
       call newdacop(deldol,c)
    else

       DO I=1,NM
          C%R(I)=RA*A%R(I)+RB*B%R(I)

       ENDDO

    ENDIF

  end subroutine NEWDALIN

  subroutine NEWDAADD(A,B,C)
    !     NEWDAADD(A,B,C):       PERFORMS C = A + B
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,B,C
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAADD"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDAADD"
          ipause=mypauses(-998,line)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          line= " NEWDAADD"
          ipause=mypauses(-999,line)
       ENDIF
    endif

    IF(PACKING)  THEN
       call newdaclr(deldol)
       do i=1,A%m
          deldol%r(A%nz(I))=a%r(A%nz(I))
          deldol%yes(A%nz(I))=.true.
          deldol%nz(i)=A%nz(I)
       enddo
       deldol%m=A%m
       do i=1,b%m
          deldol%r(b%nz(I))=b%r(b%nz(I))+deldol%r(b%nz(I))
          if(.not.deldol%yes(b%nz(I))) then
             deldol%m=deldol%m+1
             deldol%nz(deldol%m)=b%nz(I)
             deldol%yes(b%nz(I))=.true.
          endif
       enddo
       call newdacop(deldol,c)
    else

       DO I=1,NM
          C%R(I)=A%R(I)+B%R(I)
       ENDDO

    ENDIF

  end subroutine NEWDAADD

  subroutine NEWDASUB(A,B,C)
    !     NEWDASUB(A,B,C):       PERFORMS C = A - B
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,B,C
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDASUB"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDASUB"
          ipause=mypauses(-998,line)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          line= " NEWDASUB"
          ipause=mypauses(-999,line)
       ENDIF
    endif

    IF(PACKING)  THEN
       call newdaclr(deldol)
       ! call newcheck(a,-1998)
       ! call newcheck(b,-1999)
       ! call newcheck(c,-1997)

       ! call newcheck(deldol,1)

       do i=1,A%m
          deldol%r(A%nz(I))=a%r(A%nz(I))
          deldol%yes(A%nz(I))=.true.
          deldol%nz(i)=A%nz(I)
       enddo
       deldol%m=A%m

       do i=1,b%m
          deldol%r(b%nz(I))=-b%r(b%nz(I))+deldol%r(b%nz(I))
          if(.not.deldol%yes(b%nz(I))) then
             deldol%m=deldol%m+1
             deldol%nz(deldol%m)=b%nz(I)
             deldol%yes(b%nz(I))=.true.
          endif

       enddo
       ! call newcheck(a,1998)
       ! call newcheck(b,1999)
       ! call newcheck(c,1997)
       ! call newcheck(deldol,2)
       call newdacop(deldol,c)
       ! call newcheck(c,19970)
    else

       DO I=1,NM
          C%R(I)=A%R(I)-B%R(I)
       ENDDO

    ENDIF

  end subroutine NEWDASUB

  subroutine NEWDAshift(A,B,s)
    !     shift index and trashes lower ones
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,B,d
    integer i,s,test,ind(3) ,ind0(3)
    logical(lp) shift
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAshift"
          ipause=mypauses(-998,line)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          line= " NEWDAshift"
          ipause=mypauses(-999,line)
       ENDIF
    endif

    if(packing) then
       CALL allocnewda(D)

       do i=1,nm
          D%r(i)=zero
          D%yes(i)=.false.
       enddo

       do i=1,NM
          shift=.true.
          do test=1,3
             IND(TEST)=0
          enddo
          DO TEST=1,NO
             IND(TEST)=expon(TEST,i)
             IND0(TEST)=IND(TEST)
          ENDDO
          do test=1,3
             if(ind(test)>0.and.ind(test)<=s) then
                shift=.false.
             else
                if(ind(test)/=0) ind(test)=ind(test)-s
             endif
          enddo
          if(shift) then
             if(A%yes(LOCATION(IND0(1),IND0(2),IND0(3)))) then
                D%R(LOCATION(IND(1),IND(2),IND(3)))=A%R(LOCATION(IND0(1),IND0(2),IND0(3)))
             endif
          endif
       enddo

       CALL NEWDAPACK(D)

       CALL NEWDACOP(D,B)

       CALL KILLNEWDAS(D)

    else   !2222

       CALL allocnewda(D)

       do i=1,NM
          shift=.true.
          do test=1,3
             IND(TEST)=0
          enddo
          DO TEST=1,NO
             IND(TEST)=expon(TEST,i)
             IND0(TEST)=IND(TEST)
          ENDDO
          do test=1,3
             if(ind(test)>0.and.ind(test)<=s) then
                shift=.false.
             else
                if(ind(test)/=0) ind(test)=ind(test)-s
             endif
          enddo
          if(shift) then
             D%R(LOCATION(IND(1),IND(2),IND(3)))=A%R(LOCATION(IND0(1),IND0(2),IND0(3)))
          endif
       enddo

       CALL NEWDACOP(D,B)

       CALL KILLNEWDAS(D)

    endif

  end subroutine NEWDAshift


  subroutine NEWDASQR(A,C)
    !     DASQR(A,C):         PERFORMS C = A^2           (SQUARE OF A)
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,C
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDASQR"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDASQR"
          ipause=mypauses(-998,line)
       ENDIF
    endif

    call newdamul(A,A,C)
  end subroutine NEWDASQR

  subroutine NEWDACON(A,ra)
    !     DACON(A,R):         SETS A TO CONSTANT R
    implicit none
    integer ipause, mypauses
    type (taylorlow) A
    real(dp) ra
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDACON"
          ipause=mypauses(-997,line)
       ENDIF
    endif



    IF(PACKING)  THEN
       call newdaclr(A)
       if(abs(ra)>puny) then
          A%M=1
          A%NZ(1)=1
          A%yes(1)=.true.
          a%r(1)=ra
       endif

    else
       call newdaclr(A)
       a%r(1)=ra

    endif



  end subroutine NEWDACON

  subroutine NEWDAABS(A,ra)
    !     DAABS(A,R):         PERFORMS R = |A|           (NORM OF A)
    implicit none
    integer ipause, mypauses
    type (taylorlow) A
    real(dp) ra
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAABS"
          ipause=mypauses(-997,line)
       ENDIF
    endif


    IF(PACKING)  THEN
       ra=zero
       do i=1,a%m
          ra=abs(A%r(A%nz(I)))+ra
       enddo
    else

       ra=zero
       do i=1,nm
          ra=   abs(a%r(i))+ra
       enddo

    endif


  end subroutine NEWDAABS

  subroutine NEWDAPOS(A,C)
    !     DAPOS(A,C):         PERFORMS C(I) = |A(I)|     (MAKE SIGNS POSITIVE)
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,C
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAPOS"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDAPOS"
          ipause=mypauses(-998,line)
       ENDIF
    endif

    IF(PACKING)  THEN
       if(c%id/=a%id) then
          call newdaclr(c)
       endif
       c%m=a%m
       do i=1,c%m
          C%nz(I)=A%nz(I)
          C%yes(A%nz(I))=A%yes(A%nz(I))
          c%r(A%nz(I))= abs(a%r(A%nz(I)))
       enddo
    else
       do i=1,nm
          c%r(i)= abs(a%r(i))
       enddo
    ENDIF

  end subroutine NEWDAPOS

  subroutine NEWDACOM(A,C,ra)
    !     DACOM(A,C,R):       PERFORMS R = |A-C|         (NORM OF A-C)
    implicit none
    integer ipause, mypauses
    type (taylorlow) A,C
    real(dp) ra,r1,r2
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDACOM"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDACOM"
          ipause=mypauses(-998,line)
       ENDIF
    endif


    IF(PACKING)  THEN
       ra=zero

       do i=1,nm

          r1=zero
          r2=zero
          if(a%yes(i)) then
             r1=a%r(i)
          endif
          if(c%yes(i)) then
             r2=c%r(i)
          endif
          ra=   abs(r1-r2)+ra
       enddo

    else
       ra=zero
       do i=1,nm
          ra=   abs(a%r(i)-c%r(i))+ra
       enddo

    ENDIF





  end subroutine NEWDACOM

  subroutine NEWDADER(I,A,C)
    !     DADER(I,A,C):       PERFORMS C = DA/DI (DERIV. WITH RESPECT TO VARIABLE I)
    implicit none
    integer ipause, mypauses,mypause
    type (taylorlow) A,C,D
    integer i,J1,TEST,IND(3)
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDADER"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDADER"
          ipause=mypauses(-998,line)
       ENDIF
    endif

    if(packing) then
       CALL allocnewda(D)

       DO I=1,NM
          D%R(I)=zero
          D%yes(i)=.false.
       ENDDO



       ! position 1 is a constant
       DO J1=2,nm

          do test=1,3
             IND(TEST)=0
          enddo
          DO TEST=1,NO
             IND(TEST)=expon(TEST,J1)
          ENDDO

          DO TEST=1,NO
             if(ind(test)==I) THEN
                IND(TEST)=0
                IF(LOCATION(IND(1),IND(2),IND(3))>NM) then
                   ipause=mypause(-888)
                endif
                IF(LOCATION(IND(1),IND(2),IND(3))<1) then
                   ipause=mypause(-889)
                endif
                if(A%yes(J1)) then
                   D%R(LOCATION(IND(1),IND(2),IND(3)))=D%R(LOCATION(IND(1),IND(2),IND(3))) + A%R(J1)
                endif
                IND(TEST)=I
             ENDIF
          ENDDO




       ENDDO

       CALL NEWDAPACK(D)

       CALL NEWDACOP(D,C)

       CALL KILLNEWDAS(D)


    else
       CALL allocnewda(D)

       ! position 1 is a constant
       DO J1=2,nm

          do test=1,3
             IND(TEST)=0
          enddo
          DO TEST=1,NO
             IND(TEST)=expon(TEST,J1)
          ENDDO

          DO TEST=1,NO
             if(ind(test)==I) THEN
                IND(TEST)=0
                IF(LOCATION(IND(1),IND(2),IND(3))>NM) then
                   ipause=mypause(-888)
                endif
                IF(LOCATION(IND(1),IND(2),IND(3))<1) then
                   ipause=mypause(-889)
                endif
                D%R(LOCATION(IND(1),IND(2),IND(3)))=D%R(LOCATION(IND(1),IND(2),IND(3))) + A%R(J1)
                IND(TEST)=I
             ENDIF
          ENDDO




       ENDDO


       CALL NEWDACOP(D,C)

       CALL KILLNEWDAS(D)

    endif
  end subroutine NEWDADER

  subroutine NEWDATRA(I,A,C)
    !     NEWDATRA(I,A,C):
    implicit none
    integer ipause, mypauses,mypause
    real(dp) RK
    type (taylorlow) A,C,D
    integer i,J1,TEST,IND(3) ,LASTTEST
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDATRA"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDATRA"
          ipause=mypauses(-998,line)
       ENDIF
    endif


    if(packing) then
       CALL allocnewda(D)

       DO I=1,NM
          D%R(I)=zero
          D%yes(i)=.false.
       ENDDO

       ! position 1 is a constant
       DO J1=2,nm
          do test=1,3
             IND(TEST)=0
          enddo
          DO TEST=1,NO
             IND(TEST)=expon(TEST,J1)
          ENDDO

          RK=zero
          LASTTEST=0
          DO TEST=1,NO
             if(ind(test)==I) THEN
                LASTTEST=test
             ENDIF
          ENDDO

          IF(LASTTEST>0) THEN
             IND(LASTTEST)=0
             IF(LOCATION(IND(1),IND(2),IND(3))>NM) then
                ipause=mypause(-888)
             endif
             IF(LOCATION(IND(1),IND(2),IND(3))<1)  then
                ipause=mypause(-889)
             endif

             if(A%yes(J1)) then
                D%R(LOCATION(IND(1),IND(2),IND(3)))=D%R(LOCATION(IND(1),IND(2),IND(3))) + A%R(J1)
             endif

             IND(LASTTEST)=I
          ENDIF


       ENDDO


       CALL NEWDAPACK(D)

       CALL NEWDACOP(D,C)

       CALL KILLNEWDAS(D)


    else

       CALL allocnewda(D)

       ! position 1 is a constant
       DO J1=2,nm
          do test=1,3
             IND(TEST)=0
          enddo
          DO TEST=1,NO
             IND(TEST)=expon(TEST,J1)
          ENDDO

          RK=zero
          LASTTEST=0
          DO TEST=1,NO
             if(ind(test)==I) THEN
                LASTTEST=test
             ENDIF
          ENDDO

          IF(LASTTEST>0) THEN
             IND(LASTTEST)=0
             IF(LOCATION(IND(1),IND(2),IND(3))>NM) then
                ipause=mypause(-888)
             endif
             IF(LOCATION(IND(1),IND(2),IND(3))<1)  then
                ipause=mypause(-889)
             endif
             D%R(LOCATION(IND(1),IND(2),IND(3)))=D%R(LOCATION(IND(1),IND(2),IND(3))) + A%R(J1)
             IND(LASTTEST)=I
          ENDIF


       ENDDO


       CALL NEWDACOP(D,C)

       CALL KILLNEWDAS(D)

    endif

  end subroutine NEWDATRA

  subroutine NEWDAPOI(A,B,C,I)
    !     DAPOI(A,B,C,I):     PERFORMS C = [A,B] (POISSON BRACKET, 2*I: # PHASEVARS
    implicit none
    integer ipause, mypauses
    INTEGER I,K
    type (taylorlow) A,C,B,D,E,F
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAPOI"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(b%R)) THEN
          line= " NEWDAPOI"
          ipause=mypauses(-999,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDAPOI"
          ipause=mypauses(-998,line)
       ENDIF
    endif

    CALL allocnewda(D)

    CALL allocnewda(E)

    CALL allocnewda(F)


    ! position 1 is a constant
    DO K=1,I
       CALL NEWDADER(2*k-1,A,D)
       CALL NEWDADER(2*k,B,E)
       CALL NEWDAMUL(D,E,E)
       CALL NEWDAADD(E,F,F)
       CALL NEWDADER(2*k,A,D)
       CALL NEWDADER(2*k-1,B,E)
       CALL NEWDAMUL(D,E,E)
       CALL NEWDASUB(F,E,F)
    ENDDO

    CALL NEWDACOP(F,C)

    CALL KILLNEWDAS(D)
    CALL KILLNEWDAS(E)
    CALL KILLNEWDAS(F)


  end subroutine NEWDAPOI

  subroutine NEWDAVAR(A,R,I)
    !     DAVAR(A,R,I):       MAKES A INDEPENDENT VARIABLE # I WITH INITIAL VALUE R
    implicit none
    integer ipause, mypauses
    INTEGER I
    real(dp) R
    type (taylorlow) A
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAVAR"
          ipause=mypauses(-997,line)
       ENDIF
    endif

    IF((i<=0).OR.(I>NV)) THEN
       line= " NEWDAVAR"
       ipause=mypauses(-998,line)
    ENDIF


    IF(PACKING)  THEN
       CALL   NEWDACLR(A)
       if(abs(r)>puny) then
          A%M=2
          A%NZ(1)=1
          A%NZ(2)=I+1
          A%yes(A%NZ(1))=.true.
          A%yes(A%NZ(2))=.true.
          A%R(I+1)=one
          A%R(1)=R
       else
          A%M=1
          A%NZ(1)=I+1
          A%yes(A%NZ(1))=.true.
          A%R(1)=zero
          A%R(I+1)=one
       endif
       ! call newcheck(a,-1)
    ELSE
       CALL   NEWDACLR(A)
       A%R(I+1)=one
       A%R(1)=R
    ENDIF

  end subroutine NEWDAVAR

  subroutine NEWDAMULIN(A,B,RA,C,D,RB,E)
    !DAMULIN(A,B,RA,C,D,RB,E):    PERFORMS E = A*B*RA + C*D*RB
    implicit none
    integer ipause, mypauses
    real(dp) RA,RB
    type (taylorlow) A,C,B,D,E,F,G
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDAMULIN"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(b%R)) THEN
          line= " NEWDAMULIN"
          ipause=mypauses(-999,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDAMULIN"
          ipause=mypauses(-998,line)
       ENDIF
       IF (.NOT. ASSOCIATED(D%R)) THEN
          line= " NEWDAMULIN"
          ipause=mypauses(-990,line)
       ENDIF
       IF (.NOT. ASSOCIATED(E%R)) THEN
          line= " NEWDAMULIN"
          ipause=mypauses(-991,line)
       ENDIF
    endif

    CALL allocnewda(F)
    CALL allocnewda(G)

    call newdamul(a,b,f)
    call newdamul(c,d,g)
    call newdalin(f,ra,g,rb,e)

    CALL KILLNEWDAS(F)
    CALL KILLNEWDAS(G)


  end subroutine NEWDAMULIN

  subroutine NEWDADIV(A,B,C)
    implicit none
    integer ipause, mypauses
    integer i1 ,i

    type (taylorlow) A,C,B
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDADIV"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(B%R)) THEN
          line= " NEWDADIV"
          ipause=mypauses(-999,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDADIV"
          ipause=mypauses(-998,line)
       ENDIF
    endif

    if(packing) then   !2222
       if(non1) then ! for nonlinear

          CALL NEWDAFUN('INV ',B,Del2dol)
          CALL NEWDAMUL(A,Del2dol,C)

       else

          call newdaclr(deldol)
          do i=1,A%m
             deldol%r(A%nz(I))=a%r(A%nz(I))/b%R(1)
             deldol%yes(A%nz(I))=.true.
             deldol%nz(i)=A%nz(I)
          enddo
          deldol%m=A%m

          if(   a%yes(1)) then
             deldol%r(1)=deldol%r(1)+a%r(1)/b%R(1) !! doing twice but subtracting
             do i=1,b%m
                deldol%r(b%nz(I))=-b%r(b%nz(I))*A%R(1)/b%R(1)/b%R(1)+deldol%r(b%nz(I))

                if(.not.deldol%yes(b%nz(I))) then
                   deldol%m=deldol%m+1
                   deldol%nz(deldol%m)=b%nz(I)
                   deldol%yes(b%nz(I))=.true.
                endif
             enddo
          endif


          call newdacop(deldol,c)
       endif    ! for linear

    else   !2222

       if(non1) then

          CALL NEWDAFUN('INV ',B,Del2dol)
          CALL NEWDAMUL(A,Del2dol,C)

       else
          do i1=2,nm
             C%R(i1)=A%R(i1)/b%R(1)-B%R(i1)*A%R(1)/b%R(1)/b%R(1)
          enddo

          c%R(1)=A%R(1)/b%R(1)

       endif

    endif !2222


  end subroutine NEWDADIV


  subroutine NEWDADIC(A,RA,C)
    !     DADIC(A,RA,C):      PERFORMS C = RA / A
    implicit none
    integer ipause, mypauses
    INTEGER I1 ,i
    type (taylorlow) A,C
    real(dp) RA
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " NEWDADIC"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(C%R)) THEN
          line= " NEWDADIC"
          ipause=mypauses(-998,line)
       ENDIF
    endif


    IF(PACKING)  THEN !222


       if(non1) then
          CALL NEWDAFUN('INV ',A,Del2dol)
          CALL NEWDACMU(Del2dol,RA,C)
       else

          !do i1=2,nm
          !C%R(i1)=-A%R(i1)*RA/A%R(1)/A%R(1)
          !enddo
          !c%R(1)=RA/A%R(1)

          if(c%id/=a%id) then
             call newdaclr(c)
          endif

          c%m=a%m
          do i=1,c%m
             C%nz(I)=A%nz(I)
             C%yes(A%nz(I))=A%yes(A%nz(I))
             c%r(A%nz(I))=-A%R(A%nz(I))*RA/A%R(1)/A%R(1)
          enddo
          c%R(1)=-c%R(1)     !  factor of - because sum carried over i1=1,nm
       endif



    else  !222

       if(non1) then

          CALL NEWDAFUN('INV ',A,Del2dol)
          CALL NEWDACMU(Del2dol,RA,C)

       else

          do i1=2,nm
             C%R(i1)=-A%R(i1)*RA/A%R(1)/A%R(1)
          enddo

          c%R(1)=RA/A%R(1)




       endif


    endif !222

  end subroutine NEWDADIC


  SUBROUTINE NEWDAFUN(CF,INA,INC)
    implicit none
    integer ipause, mypauses
    CHARACTER(4) CF
    real(dp) XF(0:5),A0, RA,ERA,EA,SA,CA,SCR,E1,A1,A2,A3,A4,A5,P,T,E2
    INTEGER I
    type (taylorlow) INA,INC  !,IPOW,INON,ISCR
    do i=0,5
       xf(i)=zero
    enddo
    if(checkass) then
       IF (.NOT. ASSOCIATED(INA%R)) THEN
          line= " NEWDAFUN"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(INC%R)) THEN
          line= " NEWDAFUN"
          ipause=mypauses(-999,line)
       ENDIF
    endif


    IF(CF.EQ.'SQR ') THEN
       CALL NEWDASQR(INA,INC)
       RETURN
    ENDIF
    !
    !     ALLOCATE VARIABLES, PICK ZEROTH ORDER TERM
    !     ******************************************
    !
    !
    !CALL allocnewda(IPOW)
    CALL NEWDACLR(IPOW)
    !CALL allocnewda(INON)
    CALL NEWDACLR(INON)
    !CALL allocnewda(ISCR)
    CALL NEWDACLR(ISCR)

    !      CALL DAPEK(INA,JJ,A0)
    if(packing) then
       if(ina%yes(1)) then
          A0=INA%R(1)
       else
          A0=zero
       endif
    else
       A0=INA%R(1)
    endif
    !      NO = MIN(NOCUT,INOA,INOC)
    !
    !     BRANCHING TO DIFFERENT FUNCTIONS
    !     ********************************
    !
    IF(CF.EQ.'INV ') THEN
       !        1/(A0+P) = 1/A0*(1-(P/A0)+(P/A0)**2-...)
       IF(A0.EQ.zero) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-1)
       ENDIF
       XF(0) = one/A0
       DO I=1,NO
          XF(I) = -XF(I-1)/A0
       enddo
       !
    ELSEIF(CF.EQ.'SQRT') THEN
       !        SQRT(A0+P) = SQRT(A0)*(1+1/2(P/A0)-1/8*(P/A0)**2+...)
       IF(A0.LE.0) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-2)
       ENDIF
       RA = SQRT(A0)
       XF(0) = RA
       DO I=1,NO
          XF(I) = -XF(I-1)/A0/REAL(2*I,kind=DP)*REAL(2*I-3,kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'ISRT') THEN
       !        1/SQRT(A0+P) = 1/SQRT(A0)*(1-1/2(P/A0)+3/8*(P/A0)**2-...)
       IF(A0.LE.0) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-3)
       ENDIF
       ERA = one/SQRT(A0)
       XF(0) = ERA
       DO I=1,NO
          XF(I) = -XF(I-1)/A0/REAL(2*I,kind=DP)*REAL(2*I-1,kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'EXP ') THEN
       !        EXP(A0+P) = EXP(A0)*(1+P+P**2/2!+...)
       EA  = EXP(A0)
       XF(0) = EA
       DO I=1,NO
          XF(I) = XF(I-1)/REAL(I,kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'LOG ') THEN
       !        LOG(A0+P) = LOG(A0) + (P/A0) - 1/2*(P/A0)**2 + 1/3*(P/A0)**3 - ...)
       IF(A0.LE.0) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-4)
       ENDIF
       EA  = LOG(A0)
       XF(0) = EA
       XF(1) = one/A0
       DO I=2,NO
          XF(I) = -XF(I-1)/A0/REAL(I,kind=DP)*REAL(I-1,kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'SIN ') THEN
       !        SIN(A0+P) = SIN(A0)*(1-P**2/2!+P**4/4!) + COS(A0)*(P-P**3/3!+P**5/5!)
       SA  = SIN(A0)
       CA  = COS(A0)
       XF(0) = SA
       XF(1) = CA
       DO I=2,NO
          XF(I) = -XF(I-2)/REAL(I*(I-1),kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'COS ') THEN
       !        COS(A0+P) = COS(A0)*(1-P**2/2!+P**4/4!) - SIN(A0)*(P-P**3/3!+P**5/5!)
       SA  = SIN(A0)
       CA  = COS(A0)
       XF(0) = CA
       XF(1) = -SA
       DO I=2,NO
          XF(I) = -XF(I-2)/REAL(I*(I-1),kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'SIRX') THEN
       !        SIN(SQRT(P))/SQRT(P) = 1 - P/3! + P**2/5! - P**3/7! + ...
       IF(A0.NE.0) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-5)
       ENDIF
       XF(0)=one
       DO I=1,NO
          XF(I) = -XF(I-1)/REAL(2*I*(2*I+1),kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'CORX') THEN
       !        COS(SQRT(P)) = 1 - P/2! + P**2/4! - P**3/6! + ...
       IF(A0.NE.0) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-6)
       ENDIF
       XF(0)=one
       DO I=1,NO
          XF(I) = -XF(I-1)/REAL(2*I*(2*I-1),kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'SIDX') THEN
       !        SIN(P)/P = 1 - P**2/3! + P**4/5! - P**6/7! + ...
       IF(A0.NE.0) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-7)
       ENDIF
       XF(0)=one
       XF(1)=zero
       DO I=2,NO
          XF(I) = -XF(I-2)/REAL(I*(I+1),kind=DP)
       enddo
       !
    ELSEIF(CF.EQ.'TAN ') THEN
       IF(ABS(COS(A0)).LT.EPSdolMAC) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-8)
       ENDIF
       SA  = SIN(A0)
       CA  = COS(A0)
       XF(0) = SA/CA
       XF(1) = one/CA/CA
       XF(2) = two*SA/CA/CA/CA/two
       XF(3) = (two*CA*CA+six*SA*SA)/CA/CA/CA/CA/six
       XF(4) = (16*SA+eight*SA*SA*SA)/CA/CA/CA/CA/CA/c_24
       XF(5) = (c_16*CA*CA+c_24*CA*CA*SA*SA+c_80*SA*SA+c_40*SA*SA*SA*SA)
       XF(5) = XF(5)/CA/CA/CA/CA/CA/CA/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'COT ') THEN
       IF(ABS(SIN(A0)).LT.EPSdolMAC) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-16)
       ENDIF
       SA  = SIN(A0)
       CA  = COS(A0)
       XF(0) = CA/SA
       XF(1) = -one/SA/SA
       XF(2) = two*CA/SA/SA/SA/two
       XF(3) = -(two*SA*SA+six*CA*CA)/SA/SA/SA/SA/six
       XF(4) = (16*CA+eight*CA*CA*CA)/SA/SA/SA/SA/SA/c_24
       XF(5) = -(c_16*SA*SA+c_24*SA*SA*CA*CA+c_80*CA*CA+c_40*CA*CA*CA*CA)
       XF(5)=XF(5)/SA/SA/SA/SA/SA/SA/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ASIN') THEN
       IF((one-ABS(A0)).LT.zero) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-9)
       ENDIF
       XF(0) = ASIN(A0)
       XF(1) = (one-A0*A0)**(-half)
       XF(2) = A0*XF(1)**three/two
       XF(3) = (1+two*A0*A0)*XF(1)**five/six
       XF(4) = (nine*A0+six*A0*A0*A0)*XF(1)**seven/c_24
       XF(5) = (nine+c_72*A0*A0+c_24*A0*A0*A0*A0)*XF(1)**nine/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ACOS')THEN
       IF((one-ABS(A0)).LT.zero) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-10)
       ENDIF
       XF(0) =  ACOS(A0)
       SCR =  (one-A0*A0)**(-half)
       XF(1) =  -SCR
       XF(2) = -A0*SCR**three/two
       XF(3) = -(1+two*A0*A0)*SCR**five/six
       XF(4) = -(nine*A0+six*A0*A0*A0)*SCR**seven/c_24
       XF(5) = -(nine+c_72*A0*A0+c_24*A0*A0*A0*A0)*SCR**nine/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ATAN') THEN
       !        ATAN(A0+P) = ATAN(A0)+1/(1+A0**2)*P-A0/(1+A0**2)**2*P**2+....)
       XF(0) = ATAN(A0)
       XF(1) = one/(one+A0*A0)
       XF(2) = -A0*(XF(1)*XF(1))
       XF(3) = (A0*A0-one/three)*XF(1)**3
       XF(4) = (A0-A0*A0*A0)*XF(1)**4
       XF(5) = (one/five+A0**4-two*A0*A0)*XF(1)**5
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ACOT') THEN
       XF(0) = two*ATAN(one)-ATAN(A0)
       SCR = one/(one+A0*A0)
       XF(1) = -SCR
       XF(2) = A0*(SCR*SCR)
       XF(3) = -(A0*A0-one/three)*SCR**3
       XF(4) = -(A0-A0*A0*A0)*SCR**4
       XF(5) = -(one/five+A0**4-two*A0*A0)*SCR**5
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    elseif(cf.eq.'SINH') then
       sa  = SINH(a0)
       ca  = COSH(a0)
       xf(0) = sa
       xf(1) = ca
       do i=2,no
          xf(i) = xf(i-2)/REAL(i*(i-1),kind=DP)
       enddo
    elseif(cf.eq.'COSH') then
       sa  = SINH(a0)
       ca  = COSH(a0)
       xf(0) = ca
       xf(1) = sa
       xf(0) = ca
       xf(1) = sa
       do i=2,no
          xf(i) = xf(i-2)/REAL(i*(i-1),kind=DP)
       enddo
    ELSEIF(CF.EQ.'TANH') THEN
       SA  = SINH(A0)
       CA  = COSH(A0)
       XF(0) = SA/CA
       XF(1) = one/CA/CA
       XF(2) = -two*SA/CA/CA/CA/two
       XF(3) = (-two*CA*CA+six*SA*SA)/CA/CA/CA/CA/six
       XF(4) = (16*SA-eight*SA*SA*SA)/CA/CA/CA/CA/CA/c_24
       XF(5) = (c_16*CA*CA-c_24*CA*CA*SA*SA-c_80*SA*SA+c_40*SA*SA*SA*SA)
       XF(5)=XF(5)/CA/CA/CA/CA/CA/CA/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'COTH') THEN
       IF(ABS(SINH(A0)).LT.EPSdolMAC) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-11)
       ENDIF
       SA  = SINH(A0)
       CA  = COSH(A0)
       XF(0) = CA/SA
       XF(1) = -one/SA/SA
       XF(2) =  two*CA/SA/SA/SA/two
       XF(3) = (two*SA*SA-six*CA*CA)/SA/SA/SA/SA/six
       XF(4) = (16*CA+eight*CA*CA*CA)/SA/SA/SA/SA/SA/c_24
       XF(5) = (c_16*SA*SA+c_24*SA*SA*CA*CA-c_80*CA*CA-c_40*CA*CA*CA*CA)
       XF(5)=XF(5)/SA/SA/SA/SA/SA/SA/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ASNH') THEN
       XF(0) = LOG(A0+SQRT(A0*A0+one))
       XF(1) = (one+A0*A0)**(-half)
       XF(2) = -A0*XF(1)**three/two
       XF(3) = (two*A0*A0-one)*XF(1)**five/six
       XF(4) = (nine*A0-six*A0*A0*A0)*XF(1)**seven/c_24
       XF(5) = (nine-c_72*A0*A0+c_24*A0*A0*A0*A0)*XF(1)**nine/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ACSH') THEN
       IF((one-A0).GE.zero) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-12)
       ENDIF
       XF(0) = LOG(A0+SQRT(A0*A0-one))
       XF(1) = (A0*A0-one)**(-half)
       XF(2) = -A0*XF(1)**three/two
       XF(3) = (two*A0*A0+one)*XF(1)**five/six
       XF(4) = (-nine*A0-six*A0*A0*A0)*XF(1)**seven/c_24
       XF(5) = (nine+c_72*A0*A0+c_24*A0*A0*A0*A0)*XF(1)**nine/c_120
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ATNH') THEN
       IF((ABS(A0)-one).GE.zero) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-13)
       ENDIF
       XF(0) =  half*LOG((1+A0)/(1-A0))
       XF(1) =  one/(one-A0*A0)
       XF(2) =  A0*(XF(1)*XF(1))
       XF(3) = (A0*A0+one/three)*XF(1)**3
       XF(4) = (A0+A0*A0*A0)*XF(1)**4
       XF(5) = (one/five+A0**4+two*A0*A0)*XF(1)**5
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ACTH') THEN
       IF(one-ABS(A0).GE.zero) THEN
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,A72))'
          write(W_P%C(1),'(a4,1x,g20.14)')  CF,A0
          CALL WRITE_E(-14)
       ENDIF
       XF(0) =  half*LOG((A0+1)/(A0-1))
       SCR =  one/(-one+A0*A0)
       XF(1) = -SCR
       XF(2) =  A0*(SCR*SCR)
       XF(3) = (-A0*A0-one/three)*SCR**three
       XF(4) = (A0+A0*A0*A0)*SCR**four
       XF(5) = (-one/five-A0**4-two*A0*A0)*SCR**five
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSEIF(CF.EQ.'ERF ') THEN
       !
       !    ERF(X) STELLT DAS INTEGRAL VON 0 BIS X VON [ 2/SQRT(PI) * EXP(-X*X) ]
       !    DAR
       !
       E1 = EXP(-A0*A0)
       A1 = c_0_254829592
       A2 = -c_0_284496736
       A3 = c_1_421413741
       A4 = -c_1_453152027
       A5 = c_1_061405429
       P  = c_0_3275911
       T  = one/(one+P*A0)
       E2 = one-T*(A1+T*(A2+T*(A3+T*(A4+T*A5))))*E1
       XF(0)= E2
       rpi4=sqrt(pi)/two
       XF(1) = E1/RPI4
       XF(2) = -A0*E1/RPI4
       XF(3) = (-two+four*A0*A0)/six*E1/RPI4
       XF(4) = (twelve*A0-eight*A0*A0*A0)/c_24*E1/RPI4
       XF(5) = (c_16*A0*A0*A0*A0-c_48*A0*A0+twelve)/c_120*E1/RPI4
       IF(NO.GT.5) THEN
          line='ERROR IN DAFUN, '//CF//' ONLY UP TO NO = 5'
          ipause=mypauses(-1,line)
       ENDIF
    ELSE
       line='ERROR, UNSOPPORTED FUNCTION '//cf
       ipause=mypauses(-1,line)
    ENDIF
    !
    CALL NEWDACON(INC,XF(0))
    CALL NEWDACOP(INA,INON)
    INON%R(1)=zero
    !       IPOW%R(1)=one
    call newdacon(ipow,one)
    !      CALL DAPOK(INON,JJ,zero)
    !      CALL DACON(IPOW,one)
    !
    DO I=1,  NO   !MIN(NO,NOCUT)
       !
       CALL NEWDAMUL(INON,IPOW,ISCR)
       CALL NEWDACOP(ISCR,IPOW)
       CALL NEWDACMA(INC,IPOW,XF(I),INC)
    enddo

    !CALL KILLNEWDAS(ISCR)
    !CALL KILLNEWDAS(INON)
    !CALL KILLNEWDAS(IPOW)

  END SUBROUTINE  NEWDAFUN

  subroutine newdancd(id,j)
    implicit none
    integer i,id
    integer,dimension(:)::J

    do i=1,lnv
       j(i)=0
    enddo

    do i=1,no
       if(expon(i,id)/=0) J(expon(i,id))=J(expon(i,id))+1
    enddo

  end  subroutine newdancd


  SUBROUTINE newDACFU(INA,FUN,INC)
    !     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
    implicit none
    integer ipause, mypauses
    integer J(LNV)
    real(dp) cfac
    INTEGER I
    !    real(dp) fun
    !    external  fun
    type (taylorlow) INA,INC

    interface
       !       real(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         real(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface


    if(checkass) then
       IF (.NOT. ASSOCIATED(INA%R)) THEN
          line= " newDACFU"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(INC%R)) THEN
          line= " newDACFU"
          ipause=mypauses(-999,line)
       ENDIF
    endif


    IF(PACKING)  THEN  !222
       ! call newcheck(ina,907)
       ! call newcheck(inc,908)
       del2dol%m=0
       do i=1,nm
          del2dol%r(i)=zero
          del2dol%nz(i)=0
          del2dol%yes(i)=.false.
       enddo

       DO  I=1,nm
          CALL newDANCD(i,J)
          if(abs(ina%r(i))>EPSdol) then
             CFAC = FUN(J)
             del2dol%r(i) = ina%r(i) *CFAC
          else
             del2dol%r(i) = zero
          endif
       enddo

       call newdapack(del2dol)
       ! call newcheck(del2dol,909)
       call newdacop(del2dol,inc)
       ! call newcheck(inc,910)

    else   !222
       DO  I=1,nm
          CALL newDANCD(i,J)
          if(abs(ina%r(i))>EPSdol) then
             CFAC = FUN(J)
             inc%r(i) = ina%r(i) *CFAC
          else
             inc%r(i) = zero
          endif
       enddo

    endif !222
  END    SUBROUTINE newDACFU

  SUBROUTINE newDACFUR(INA,FUN,INC)
    !     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
    implicit none
    integer ipause, mypauses
    integer J(LNV)
    real(dp) cfac
    INTEGER I
    !    complex(dp) fun
    !    external  fun
    type (taylorlow) INA,INC

    interface
       !       complex(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface

    if(checkass) then
       IF (.NOT. ASSOCIATED(INA%R)) THEN
          line= " newDACFUR"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(INC%R)) THEN
          line= " newDACFUR"
          ipause=mypauses(-999,line)
       ENDIF
    endif


    IF(PACKING)  THEN  !222
       ! call newcheck(ina,907)
       ! call newcheck(inc,908)
       del2dol%m=0
       do i=1,nm
          del2dol%r(i)=zero
          del2dol%nz(i)=0
          del2dol%yes(i)=.false.
       enddo

       DO  I=1,nm
          CALL newDANCD(i,J)
          if(abs(ina%r(i))>EPSdol) then
             CFAC = REAL(FUN(J),kind=DP)
             del2dol%r(i) = ina%r(i) *CFAC
          else
             del2dol%r(i) = zero
          endif
       enddo

       call newdapack(del2dol)
       ! call newcheck(del2dol,909)
       call newdacop(del2dol,inc)
       ! call newcheck(inc,910)

    else   !222
       DO  I=1,nm
          CALL newDANCD(i,J)
          if(abs(ina%r(i))>EPSdol) then
             CFAC = REAL(FUN(J),kind=DP)
             inc%r(i) = ina%r(i) *CFAC
          else
             inc%r(i) = zero
          endif
       enddo

    endif !222
  END    SUBROUTINE newDACFUR

  SUBROUTINE newDACFUI(INA,FUN,INC)
    !     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
    implicit none
    integer ipause, mypauses
    integer J(LNV)
    real(dp) cfac
    INTEGER I
    !    complex(dp) fun
    !    external  fun
    type (taylorlow) INA,INC
    interface
       !       complex(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface

    if(checkass) then
       IF (.NOT. ASSOCIATED(INA%R)) THEN
          line= " newDACFUI"
          ipause=mypauses(-997,line)
       ENDIF
       IF (.NOT. ASSOCIATED(INC%R)) THEN
          line= " newDACFUI"
          ipause=mypauses(-999,line)
       ENDIF
    endif


    IF(PACKING)  THEN  !222
       ! call newcheck(ina,907)
       ! call newcheck(inc,908)
       del2dol%m=0
       do i=1,nm
          del2dol%r(i)=zero
          del2dol%nz(i)=0
          del2dol%yes(i)=.false.
       enddo

       DO  I=1,nm
          CALL newDANCD(i,J)
          if(abs(ina%r(i))>EPSdol) then
             CFAC = aimag(FUN(J))
             del2dol%r(i) = ina%r(i) *CFAC
          else
             del2dol%r(i) = zero
          endif
       enddo

       call newdapack(del2dol)
       ! call newcheck(del2dol,909)
       call newdacop(del2dol,inc)
       ! call newcheck(inc,910)

    else   !222
       DO  I=1,nm
          CALL newDANCD(i,J)
          if(abs(ina%r(i))>EPSdol) then
             CFAC = aimag(FUN(J))
             inc%r(i) = ina%r(i) *CFAC
          else
             inc%r(i) = zero
          endif
       enddo

    endif !222
  END    SUBROUTINE newDACFUI

  SUBROUTINE newDACCTN(X,Ic,Y,Jc,Z,Kc)
    !     DACCT(X,I,Y,J,Z,K)  CONCATENATES Z = X o Y;   I,J,K: # OF VECTORS IN X,Y,Z
    implicit none
    integer ipause, mypauses
    INTEGER Ic,jc,kc
    INTEGER I,j,k,n
    type (taylorlow) d(3),temp(LNV)
    type (taylorlow),dimension(:)::x,y,z

    if(checkass) then
       IF (.NOT. ASSOCIATED(x(1)%R)) THEN
          line= " newDACCT"
          ipause=mypauses(-990,line)
       ENDIF
       IF (.NOT. ASSOCIATED(y(1)%R)) THEN
          line= " newDACCT"
          ipause=mypauses(-991,line)
       ENDIF
       IF (.NOT. ASSOCIATED(z(1)%R)) THEN
          line= " newDACCT"
          ipause=mypauses(-992,line)
       ENDIF
       if(jc.ne.nv) then
          line= " newDACCT"
          ipause=mypauses(-993,line)
       endif
       if(ic.ne.kc) then
          line= " newDACCT"
          ipause=mypauses(-994,line)
       endif
    endif

    if(packing) then !222
       DO I=1,NO
          CALL allocnewda(D(I))
       ENDDO
       do i=1,kc
          CALL allocnewda(temp(i))
       enddo
       do i=0,nv  !nof(1)
          CALL NEWDAcon(D(1),one)
          if(i.NE.0) call newdamul(y(i),d(1),D(1))
          do j=0,i*flag(2)
             IF(NO>1) CALL NEWDACOP(d(1),d(2))
             if(j.NE.0) call newdamul(y(j),d(2),d(2))
             do k=0,j*flag(3)
                IF(NO>2) CALL NEWDACOP(d(2),d(3))
                if(k.NE.0) call newdamul(y(k),d(3),d(3))
                do n=1,kc
                   if(X(N)%yes(location(i,j,k))) then
                      CALL NEWDALIN(TEMP(N),one,D(NO),X(N)%R(location(i,j,k)),TEMP(N))
                   endif
                enddo
             enddo
          enddo
       enddo
       do n=1,kc
          CALL NEWDACOP(TEMP(N),Z(N))
       enddo
       DO I=1,NO
          CALL KILLNEWDAS(D(I))
       ENDDO
       DO I=1,KC
          CALL KILLNEWDAS(TEMP(I))
       ENDDO
    else  !222
       DO I=1,NO
          CALL allocnewda(D(I))
       ENDDO
       do i=1,kc
          CALL allocnewda(temp(i))
       enddo
       do i=0,nv  !nof(1)
          CALL NEWDAcon(D(1),one)
          if(i.NE.0) call newdamul(y(i),d(1),D(1))
          do j=0,i*flag(2)
             IF(NO>1) CALL NEWDACOP(d(1),d(2))
             if(j.NE.0) call newdamul(y(j),d(2),d(2))
             do k=0,j*flag(3)
                IF(NO>2) CALL NEWDACOP(d(2),d(3))
                if(k.NE.0) call newdamul(y(k),d(3),d(3))
                do n=1,kc

                   CALL NEWDALIN(TEMP(N),one,D(NO),X(N)%R(location(i,j,k)),TEMP(N))
                enddo
             enddo
          enddo
       enddo
       do n=1,kc
          CALL NEWDACOP(TEMP(N),Z(N))
       enddo
       DO I=1,NO
          CALL KILLNEWDAS(D(I))
       ENDDO
       DO I=1,KC
          CALL KILLNEWDAS(TEMP(I))
       ENDDO
    endif  !222

  END    SUBROUTINE newDACCTN


  SUBROUTINE newDACCT0(X,Ic,Y,Jc,Z,Kc)
    !     DACCT(X,I,Y,J,Z,K)  CONCATENATES Z = X o Y;   I,J,K: # OF VECTORS IN X,Y,Z
    implicit none
    integer ipause, mypauses
    INTEGER Ic,jc,kc
    INTEGER I,j,k,n
    type (taylorlow) x,z,d(3),temp(LNV)
    type (taylorlow),dimension(:)::y

    if(checkass) then
       IF (.NOT. ASSOCIATED(x%R)) THEN
          line= " newDACCT0"
          ipause=mypauses(-990,line)
       ENDIF
       IF (.NOT. ASSOCIATED(y(1)%R)) THEN
          line= " newDACCT0"
          ipause=mypauses(-991,line)
       ENDIF
       IF (.NOT. ASSOCIATED(z%R)) THEN
          line= " newDACCT0"
          ipause=mypauses(-992,line)
       ENDIF
       if(jc.ne.nv) then
          line= " newDACCT0"
          ipause=mypauses(-993,line)
       endif
       if(ic.ne.kc) then
          line= " newDACCT0"
          ipause=mypauses(-994,line)
       endif
    endif

    if(packing) then !222
       DO I=1,NO
          CALL allocnewda(D(I))
       ENDDO
       do i=1,kc
          CALL allocnewda(temp(i))
       enddo
       do i=0,nv  !nof(1)
          CALL NEWDAcon(D(1),one)
          if(i.NE.0) call newdamul(y(i),d(1),D(1))
          do j=0,i*flag(2)
             IF(NO>1) CALL NEWDACOP(d(1),d(2))
             if(j.NE.0) call newdamul(y(j),d(2),d(2))
             do k=0,j*flag(3)
                IF(NO>2) CALL NEWDACOP(d(2),d(3))
                if(k.NE.0) call newdamul(y(k),d(3),d(3))
                do n=1,kc
                   if(X%yes(location(i,j,k))) then
                      CALL NEWDALIN(TEMP(N),one,D(NO),X%R(location(i,j,k)),TEMP(N))
                   endif
                enddo
             enddo
          enddo
       enddo
       do n=1,kc
          CALL NEWDACOP(TEMP(N),Z)
       enddo
       DO I=1,NO
          CALL KILLNEWDAS(D(I))
       ENDDO
       DO I=1,KC
          CALL KILLNEWDAS(TEMP(I))
       ENDDO
    else !222
       DO I=1,NO
          CALL allocnewda(D(I))
       ENDDO
       do i=1,kc
          CALL allocnewda(temp(i))
       enddo
       do i=0,nv  !nof(1)
          CALL NEWDAcon(D(1),one)
          if(i.NE.0) call newdamul(y(i),d(1),D(1))
          do j=0,i*flag(2)
             IF(NO>1) CALL NEWDACOP(d(1),d(2))
             if(j.NE.0) call newdamul(y(j),d(2),d(2))
             do k=0,j*flag(3)
                IF(NO>2) CALL NEWDACOP(d(2),d(3))
                if(k.NE.0) call newdamul(y(k),d(3),d(3))
                do n=1,kc
                   CALL NEWDALIN(TEMP(N),one,D(NO),X%R(location(i,j,k)),TEMP(N))
                enddo
             enddo
          enddo
       enddo
       do n=1,kc
          CALL NEWDACOP(TEMP(N),Z)
       enddo
       DO I=1,NO
          CALL KILLNEWDAS(D(I))
       ENDDO
       DO I=1,KC
          CALL KILLNEWDAS(TEMP(I))
       ENDDO
    endif !222

  END    SUBROUTINE newDACCT0

  SUBROUTINE newdapri(INA,mf)
    !     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
    implicit none
    integer ipause, mypauses
    INTEGER I,mf,k,iout
    type (taylorlow) INA

    if(checkass) then
       IF (.NOT. ASSOCIATED(INA%R)) THEN
          line= " newdapri"
          ipause=mypauses(-997,line)
       ENDIF
    endif


    iout=0
    write(mf,*) no,nv
    DO I=1,NM
       IF(abs(ina%r(I)).gt.EPSdol) then
          iout=iout+1
          IF(PACKING) THEN
             IF(INA%YES(I)) write(mf,*) i, (expon(k,i),k=1,no),ina%r(I)

          ELSE
             write(mf,*) i, (expon(k,i),k=1,no),ina%r(I)

          ENDIF

       endif
    ENDDO
    if(iout.eq.0) iout=1
    write(mf,*)-iout, (expon(k,1),k=1,no),zero


  END SUBROUTINE NEWDAPRI


  SUBROUTINE newdarea(INA,mf)
    !     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
    implicit none
    integer ipause, mypauses
    INTEGER I,mf,k,iout,no1,nv1,ind(3)
    real(dp) c
    logical(lp) check
    type (taylorlow) INA

    if(checkass) then
       IF (.NOT. ASSOCIATED(INA%R)) THEN
          line= " newdarea"
          ipause=mypauses(-997,line)
       ENDIF
    endif

    iout=0
    do i=1,3
       ind(i)=0
    enddo
    DO I=1,NM
       ina%r(I)=zero
       ina%yes(i)=.false.
    enddo

    read(mf,*) no1,nv1

    do while(iout>=0)
       read(mf,*) iout, (ind(k),k=1,no1),c
       k=0
       check=.true.
       do i=1,no1
          if(ind(i)>0) k=k+1
          if(ind(i)>nv) check=.false.
       enddo
       if((k<=no).and.check.and.(iout>=0)) ina%r(location(ind(1),ind(2),ind(3)))=c

    ENDDO

    IF(PACKING)  THEN
       CALL NEWDAPACK(INA )
    ENDIF
  END SUBROUTINE newdarea

  SUBROUTINE newDAEPS(DEPSdol)
    implicit none
    real(dp) dEPSdol
    if(dEPSdol.ge.zero) then
       EPSdol = DEPSdol
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,A72))'
       write(W_P%C(1),'(a8,1x,g20.14)')   " EPSdol ",EPSdol
       CALL WRITE_i
    else
       dEPSdol=EPSdol
    endif
  END   subroutine newDAEPS

  SUBROUTINE newDAINV(MA,IA,MBO,IB)
    implicit none
    integer ipause, mypauses
    real(dp) AA(LNV,LNV),AI(LNV,LNV),PROD,R1,R2
    type (taylorlow) MB(LNV),MS(LNV),ML(LNV)
    type (taylorlow),dimension(:)::MA,MBO
    integer i,IE,J,K,IA,IB,IER
    if(checkass) then
       DO I=1,IA
          IF (.NOT. ASSOCIATED(MA(I)%R)) THEN
             write(line,'(a10,i8)')  " newDAINV ",I
             ipause=mypauses(-990,line)
          ENDIF
          IF (.NOT. ASSOCIATED(MBO(I)%R)) THEN
             write(line,'(a10,i8)')  " newDAINV ",I
             ipause=mypauses(-991,line)
          ENDIF
       ENDDO
    endif

    if(packing) then   !2222
       do  ie=1,ib
          CALL allocnewda(MB(IE))
          CALL allocnewda(MS(IE))
          CALL allocnewda(ML(IE))
          call newdaVAR(MS(ie),zero,IE)
       ENDDO

       !     DO  I=1,IB
       !       v(i)= MA(I)%R(1)
       !        MA(I)%R(1)=zero
       !      enddo

       IF(IA.NE.IB) THEN
          line= 'ERROR IN DAINV, IA .NE. IB'
          ipause=mypauses(-1,line)
          STOP 889
       ENDIF
       !
       !     EXTRACTING LINEAR MATRIX, GENERATING NONLINEAR PART OF A
       !     ********************************************************
       !
       DO  I=1,IB
          DO  J=1,IB
             if( MA(I)%yes(1+J) ) then
                AA(I,J) = MA(I)%R(1+J)
             else
                AA(I,J) = zero
             endif
          enddo
       enddo
       !
       !     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
       !     **********************************************************
       !
       CALL MATINV(AA,AI,IA,LNV,IER)

       IF(IER.EQ.132) THEN
          line= 'ERROR IN ROUTINE DAINV'
          ipause=mypauses(-888,line)
       ENDIF
       !
       IER = 0
       DO I=1,IB
          DO J=1,IB
             PROD = zero
             DO K=1,IB
                PROD = PROD + AA(I,K)*AI(K,J)
             enddo
             IF(I.EQ.J) PROD = PROD - one
             IF(abs(PROD).GT.c_100*EPSdolMAC) THEN
                write(line,'(a50,2(1x,i4),3(1x,g12.6))')'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD = ' &
                     & ,I,J,PROD,EPSdolmac,EPSdol
                IER = 1
                ipause=mypauses(-555,line)
             ENDIF
          enddo
       enddo
       DO  I=1,IB
          DO  J=1,IB
             if(MB(I)%yes(1+J)) then
                MB(I)%R(1+J)=AI(I,J)
             else
                MB(I)%yes(1+J)=.true.
                MB(I)%R(1+J)=AI(I,J)
                MB(I)%M=MB(I)%M+1
                MB(I)%NZ(MB(I)%M)=1+J
             endif
             if(ML(I)%yes(1+J)) then
                ML(I)%R(1+J)=AI(I,J)
             else
                ML(I)%yes(1+J)=.true.
                ML(I)%R(1+J)=AI(I,J)
                ML(I)%M=ML(I)%M+1
                ML(I)%NZ(ML(I)%M)=1+J
             endif

          enddo
       enddo

       !     ITERATIVELY COMPUTING DIFFERENT PARTS OF THE INVERSE
       !     ****************************************************

       CALL newDACCT(MA,IA,MB,IB,ML,IA)

       DO I=1,IB
          CALL NEWDACOP(ML(I),MB(I))
       ENDDO

       !

       DO I=2,NO

          DO K=1,IB !k
             DO J=LOCS(I),LOCE(I)

                R1=zero
                R2=zero
                IF(MB(K)%YES(J)) R1=MB(K)%R(J)
                IF(MS(K)%YES(J)) R2=MS(K)%R(J)

                IF(MS(K)%YES(J)) THEN
                   MS(K)%R(J)=R2-R1
                ELSE
                   if(abs(r2-r1)>puny) then
                      MS(K)%YES(J)=.TRUE.
                      MS(K)%M=MS(K)%M+1
                      MS(K)%NZ(MS(K)%M)=J
                      MS(K)%R(J)=R2-R1
                   endif
                ENDIF
                !MS(K)%R(J)=MS(K)%R(J)-MB(K)%R(J)
             ENDDO
          ENDDO  !k

          IF(i<NO) CALL newDACCT(ML,IA,MS,IB,MB,IA)
       ENDDO

       DO  I=1,IB
          CALL NEWDACLR(MB(I))
          DO  J=1,IB
             if(abs(AI(I,J))>puny) then
                MB(I)%yes(1+J)=.true.
                MB(I)%M=MB(I)%M+1
                MB(I)%NZ(MB(I)%M)=1+J
                MB(I)%R(1+J)=AI(I,J)
             endif
          enddo
       enddo

       CALL newDACCT(MB,IB,MS,IA,MBO,IA)

       DO I=1,IA
          CALL KILLNEWDAS(MB(I))
          CALL KILLNEWDAS(ML(I))
          CALL KILLNEWDAS(MS(I))
       ENDDO
       !     DO  I=1,IB
       !        MBO(I)%R(1)= v(i)
       !        Ma(I)%R(1)= v(i)
       !      enddo

    else     !22222
       do  ie=1,ib
          CALL allocnewda(MB(IE))
          CALL allocnewda(MS(IE))
          CALL allocnewda(ML(IE))
          call newdaVAR(MS(ie),zero,IE)
       ENDDO

       !     DO  I=1,IB
       !       v(i)= MA(I)%R(1)
       !        MA(I)%R(1)=zero
       !      enddo

       IF(IA.NE.IB) THEN
          line='ERROR IN DAINV, IA .NE. IB'
          ipause=mypauses(-889,line)
       ENDIF
       !
       !     EXTRACTING LINEAR MATRIX, GENERATING NONLINEAR PART OF A
       !     ********************************************************
       !
       DO  I=1,IB
          DO  J=1,IB
             AA(I,J) = MA(I)%R(1+J)
          enddo
       enddo
       !
       !     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
       !     **********************************************************
       !
       CALL MATINV(AA,AI,IA,LNV,IER)

       IF(IER.EQ.132) THEN
          line='ERROR IN ROUTINE DAINV'
          ipause=mypauses(-888,line)
       ENDIF
       !
       IER = 0
       DO I=1,IB
          DO J=1,IB
             PROD = zero
             DO K=1,IB
                PROD = PROD + AA(I,K)*AI(K,J)
             enddo
             IF(I.EQ.J) PROD = PROD - one
             IF(abs(PROD).GT.c_100*EPSdolMAC) THEN
                write(line,'(a50,2(1x,i4),3(1x,g12.6))')'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD = ' &
                     & ,I,J,PROD,EPSdolmac,EPSdol
                IER = 1
                ipause=mypauses(-555,line)
             ENDIF
          enddo
       enddo
       DO  I=1,IB
          DO  J=1,IB
             MB(I)%R(1+J)=AI(I,J)
             ML(I)%R(1+J)=AI(I,J)
          enddo
       enddo

       !     ITERATIVELY COMPUTING DIFFERENT PARTS OF THE INVERSE
       !     ****************************************************

       CALL newDACCT(MA,IA,MB,IB,ML,IA)

       DO I=1,IB
          CALL NEWDACOP(ML(I),MB(I))
       ENDDO

       !

       DO I=2,NO

          DO K=1,IB !k
             DO J=LOCS(I),LOCE(I)
                MS(K)%R(J)=MS(K)%R(J)-MB(K)%R(J)
             ENDDO
          ENDDO  !k

          IF(i<NO) CALL newDACCT(ML,IA,MS,IB,MB,IA)
       ENDDO

       DO  I=1,IB
          CALL NEWDACLR(MB(I))
          DO  J=1,IB
             MB(I)%R(1+J)=AI(I,J)
          enddo
       enddo

       CALL newDACCT(MB,IB,MS,IA,MBO,IA)

       DO I=1,IA
          CALL KILLNEWDAS(MB(I))
          CALL KILLNEWDAS(ML(I))
          CALL KILLNEWDAS(MS(I))
       ENDDO
       !     DO  I=1,IB
       !        MBO(I)%R(1)= v(i)
       !        Ma(I)%R(1)= v(i)
       !      enddo

    endif !2222

  END    SUBROUTINE newDAINV

  SUBROUTINE newDAPIN(MA,IA,MB,IB,JIND)
    implicit none
    integer ipause, mypauses
    !     DAPIN(X,I,Z,K,JJ)   PARTIALLY INVERTS Z = X^-1; I,J: # OF VECTORS IN X,Y,
    !
    !     THIS SUBROUTINE PERFORMS A PARTIAL INVERSION OF THE ROWS MARKED WITH
    !     NONZERO ENTRIES IN JJ OF THE MATRIX A. THE RESULT IS STORED IN B.
    !
    INTEGER IE,I,IA,IB
    INTEGER,dimension(:)::JIND
    type (taylorlow) MN(LNV),MI(LNV),ME(LNV)
    type (taylorlow),dimension(:)::MA,MB

    IF(IA/=IB) THEN
       line= " IA NOT EQUAL IB IN NEWDAPIN  "
       ipause=mypauses(-900,line)
    ENDIF
    if(checkass) then
       DO I=1,IA
          IF (.NOT. ASSOCIATED(MA(I)%R)) THEN
             line= " newDAPIN"
             ipause=mypauses(-990,line)
          ENDIF
          IF (.NOT. ASSOCIATED(MB(I)%R)) THEN
             line= " newDAPIN"
             ipause=mypauses(-990,line)
          ENDIF
       ENDDO
    endif

    IF(PACKING) THEN !22222
       do  ie=1,ib
          CALL allocnewda(MN(IE))
          CALL allocnewda(MI(IE))
          CALL allocnewda(ME(IE))
       ENDDO


       DO  I=1,IA
          if(MI(I)%yes(1+I)) then
             MI(I)%R(1+I)=one
          else
             MI(I)%yes(1+I)=.true.
             MI(I)%R(1+I)=one
             MI(I)%M=MI(I)%M+1
             MB(I)%NZ(MI(I)%M)=1+I
          endif
       ENDDO

       DO I=1,IA
          CALL NEWDACOP(MA(I),MN(I))
          IF(JIND(I).EQ.0) CALL NEWDACOP(ME(I),MN(I))
       enddo

       CALL NEWDAINV(MN,IA,MI,IA)

       DO I=1,IA
          IF(JIND(I).EQ.0) CALL NEWDACOP(MA(I),ME(I))
       enddo

       CALL NEWDACCT(ME,IA,MI,IA,MB,IB)

       DO I=1,IA
          CALL KILLNEWDAS(ME(I))
          CALL KILLNEWDAS(MI(I))
          CALL KILLNEWDAS(MN(I))
       ENDDO
    ELSE  !2222
       do  ie=1,ib
          CALL allocnewda(MN(IE))
          CALL allocnewda(MI(IE))
          CALL allocnewda(ME(IE))
       ENDDO


       DO  I=1,IA
          ME(I)%R(1+I)=one
       ENDDO

       DO I=1,IA
          CALL NEWDACOP(MA(I),MN(I))
          IF(JIND(I).EQ.0) CALL NEWDACOP(ME(I),MN(I))
       enddo

       CALL NEWDAINV(MN,IA,MI,IA)

       DO I=1,IA
          IF(JIND(I).EQ.0) CALL NEWDACOP(MA(I),ME(I))
       enddo

       CALL NEWDACCT(ME,IA,MI,IA,MB,IB)

       DO I=1,IA
          CALL KILLNEWDAS(ME(I))
          CALL KILLNEWDAS(MI(I))
          CALL KILLNEWDAS(MN(I))
       ENDDO
    ENDIF  !222
  END  SUBROUTINE  newDAPIN


  SUBROUTINE newDADCD(JJ,IH)
    implicit none
    INTEGER I,J,K, jj(lnv),ind(3),IH


    DO I=1,3
       IND(I)=0
    ENDDO
    do i=1,lnv
       do k=nv+1,lnv
          if(jj(k)/=0) THEN
             IH=-1
             RETURN
          ENDIF
       enddo
       if(jj(i)/=0) then
          do j=1,no
             if(ind(j)==0) then
                do k=1,jj(i)
                   ind(j-1+k)=i
                enddo
                GOTO 100
             endif

          enddo
       endif
100    CONTINUE
    enddo

    IH=LOCATION(IND(1),IND(2),IND(3))
    if((ih<1).or.(ih>nm)) then
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,A72))'
       write(W_P%C(1),'(a6,1x,i8)')   " ih=  ",ih
       write(W_P%C(2),'(3(1x,i8))')   ind(1),ind(2),ind(3)
       CALL WRITE_E(14)
    endif


  END SUBROUTINE newDADCD

  SUBROUTINE newDAPEK(A,JJ,R)
    !     DAPEK(A,JJ,R):      RETURNS COEF R OF MONOMIAL WITH EXPONENTS JJ OF A
    implicit none
    integer ipause, mypauses
    INTEGER I,J,K,ind(3)
    INTEGER,dimension(:)::jj
    real(dp) R
    type (taylorlow) A

    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " newDAPEK"
          ipause=mypauses(-990,line)
       ENDIF
    endif

    DO I=1,3
       IND(I)=0
    ENDDO
    do i=1,lnv
       do k=nv+1,lnv
          if(jj(k)/=0) THEN
             R=zero
             RETURN
          ENDIF
       enddo
       if(jj(i)/=0) then
          do j=1,no
             if(ind(j)==0) then
                do k=1,jj(i)
                   ind(j-1+k)=i
                enddo
                GOTO 100
             endif

          enddo
       endif
100    CONTINUE
    enddo

    IF(PACKING) THEN
       IF(A%YES(LOCATION(IND(1),IND(2),IND(3)))) THEN
          R=A%R(LOCATION(IND(1),IND(2),IND(3)) )
       ELSE
          R=zero
       ENDIF
    ELSE
       R=A%R(LOCATION(IND(1),IND(2),IND(3)) )
    ENDIF


  END SUBROUTINE newDAPEK

  SUBROUTINE newDAPOK(A,JJ,R)
    !     DAPOK(A,JJ,R):      SETS COEF OF MONOMIAL WITH EXPONENTS JJ OF A TO R
    implicit none
    integer ipause, mypauses
    INTEGER I,J,K,ind(3)
    INTEGER,dimension(:)::jj
    real(dp) R
    type (taylorlow) A

    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " newDAPOK"
          ipause=mypauses(-990,line)
       ENDIF
    endif

    DO I=1,3
       IND(I)=0
    ENDDO
    do i=1,lnv
       do k=nv+1,lnv
          if(jj(k)/=0) THEN
             RETURN
          ENDIF
       enddo
       if(jj(i)/=0) then
          do j=1,no
             if(ind(j)==0) then
                do k=1,jj(i)
                   ind(j-1+k)=i
                enddo
                GOTO 100
             endif

          enddo
       endif
100    CONTINUE
    enddo

    IF(PACKING) THEN
       IF(ABS(R)>PUNY) THEN
          IF(.NOT.A%YES(LOCATION(IND(1),IND(2),IND(3)))) THEN
             A%YES(LOCATION(IND(1),IND(2),IND(3)))=.TRUE.
             A%M=A%M+1
             A%NZ(A%M)=LOCATION(IND(1),IND(2),IND(3))
          ENDIF
          A%R(LOCATION(IND(1),IND(2),IND(3)) )=R
       else
          IF(A%YES(LOCATION(IND(1),IND(2),IND(3)))) THEN
             A%R(LOCATION(IND(1),IND(2),IND(3)) )=zero
          ENDIF
       endif

    ELSE


       A%R(LOCATION(IND(1),IND(2),IND(3)) )=R
    ENDIF


  END SUBROUTINE newDAPOK

  SUBROUTINE newDARAN(A,CM,xran)
    !     DARAN(A,R,seed):         FILLS A WITH RANDOM NUMBERS. R: FILLFACTOR
    !     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
    !     FOR CM > 0, THE VECTOR IS FILLED WITH REALS,
    !     FOR CM < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
    !     ABS(CM) IS THE FILLING FACTOR
    !
    implicit none
    integer ipause, mypauses
    type (taylorlow) a
    real(dp) xran,cm
    integer i
    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " newDARAN"
          ipause=mypauses(-990,line)
       ENDIF
    endif
    IF(PACKING)  THEN
       CALL NEWDACLR(A )
       DO I=1,NM
          A%R(I)=zero
          a%yes(i)=.false.
       ENDDO
    ENDIF
    DO I=1,NM     !MAX
       IF(CM.GT.zero) THEN
          a%r(i) = BRAN(xRAN)
          IF(a%r(i).GT.CM) a%r(i) = zero
       ELSEIF(CM.LT.zero) THEN
          a%r(i) = INT(1+10*BRAN(xRAN))
          IF(a%r(i).GT.-ten*CM) a%r(i) = zero
       ELSE
          line= 'ERROR IN ROUTINE DARAN'
          ipause=mypauses(-888,line)
       ENDIF
    enddo
    IF(PACKING)  THEN
       CALL NEWDAPACK(A )
    ENDIF
  END  subroutine newDARAN

  SUBROUTINE oldDAPRI(a,IUNIT)
    implicit none
    integer ipause, mypauses
    !     ***************************
    !     Frank
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    integer J(LNV) ,IUNIT,iout,i1,j1 ,i


    !CHARACTER C10*10,K10*10
    type (taylorlow) A

    if(iunit.eq.0) return

    if(checkass) then
       IF (.NOT.ASSOCIATED(A%R)) THEN
          line= " OLDDAPRI"
          ipause=mypauses(-997,line)
       endif
    endif

    WRITE(IUNIT,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)') &
         '          ',', NO =',no,', NV =',nv,', INA =',a%id,   &
         '***********'//'**********************************'

    IOUT = 0

    WRITE(IUNIT,'(A)')'    I  COEFFICIENT          ORDER   EXPONENTS'

    DO I1=0,NO
       DO J1=LOCS(I1),LOCE(I1)
          CALL NEWDANCD(j1,J)

          if(abs(a%r(j1)).gt.EPSdol) then
             IF(PACKING) THEN
                IF(A%YES(J1)) THEN
                   IOUT = IOUT+1
                   WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))') iout,a%r(j1),i1,(J(I),I=1,nv)
                   WRITE(IUNIT,*) a%r(j1)
                ENDIF
             ELSE
                IOUT = IOUT+1
                WRITE(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))') iout,a%r(j1),i1,(J(I),I=1,nv)
                WRITE(IUNIT,*) a%r(j1)
             ENDIF

          endif
       ENDDO
    ENDDO


    WRITE(IUNIT,'(A)') '                                      '

  END  subroutine olddapri

  SUBROUTINE OLDDAPRI77(A,IUNIT)
    implicit none
    integer ipause, mypauses
    INTEGER J(LNV) ,IUNIT,iout,i1,j1 ,i


    CHARACTER(10) C10,K10
    type (taylorlow) A

    if(iunit.eq.0) return


    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " OLDDAPRI77"
          ipause=mypauses(0,line)
       ENDIF
    endif

    WRITE(IUNIT,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)')   &
         "          ",', NO =',no,', NV =',nv,', INA =',a%id,   &
         '***********'//'**********************************'
    WRITE(IUNIT,'(A)')'    I  COEFFICIENT          ORDER   EXPONENTS'
    C10='      NO ='
    K10='      NV ='

    WRITE(IUNIT,'(A10,I6,A10,I6)') C10,no,K10,nv

    IOUT = 0

    DO I1=0,NO
       DO J1=LOCS(I1),LOCE(I1)
          CALL NEWDANCD(j1,J)

          if(abs(a%r(j1)).gt.EPSdol) then
             IF(PACKING) THEN
                IF(A%YES(j1)) THEN
                   WRITE(IUNIT,501) i1,a%r(j1),(J(I),I=1,nv)
                   IOUT = IOUT+1
                ENDIF
             ELSE
                WRITE(IUNIT,501) i1,a%r(j1),(J(I),I=1,nv)
                IOUT = IOUT+1
             ENDIF
          endif
       ENDDO
    ENDDO


    DO I=1,LNV
       J(I)=0
    enddo

    IF(IOUT.EQ.0) IOUT=1

    WRITE(IUNIT,502) -IOUT,zero,(J(I),I=1,nv)
501 FORMAT(' ', I3,1X,G23.16,1X,100(1X,I2))
503 FORMAT(' ', I3,1X,G23.16,1X,100(1X,I2))
502 FORMAT(' ', I5,1X,G23.16,1X,100(1X,I2))

  END subroutine OLDDAPRI77


  SUBROUTINE olddarea(A,IUNIT)
    implicit none
    integer ipause, mypause, mypauses
    !     Frank
    CHARACTER(10) C10
    real(dp) C
    integer J(LNV),NNO,nvjoh,IO1,IUNIT,i,II,iin,IO
    !
    type (taylorlow) A

    if(iunit.eq.0) return

    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " OLDDAREA"
          ipause=mypauses(0,line)
       ENDIF
    endif

    DO I=1,LNV
       J(I) = 0
    enddo

    CALL newDACLR(A)
    IF(PACKING)  THEN
       DO I=1,NM
          A%R(I)=zero
          a%yes(i)=.false.
       ENDDO
    ENDIF

    !
    !
    nvjoh=nv
    !
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(18X,I4)') NNO
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10)') C10


    IIN = 0

10  CONTINUE
    IIN = IIN + 1
    READ(IUNIT,'(I6,2X,G20.14,I5,4X,18(2I2,1X))')II,C,IO,(J(I),I=1,nvjoh)

    IF(II.EQ.0) GOTO 20
    !ETIENNE
    READ(IUNIT,*) C
    !ETIENNE

    IO1 = 0
    DO I=1,nv
       IO1 = IO1 + J(I)
    enddo
    !
    IF(IO1.NE.IO) THEN
       GOTO 10
    ENDIF
    IF(IO.GT.no) THEN
       GOTO 10
    ENDIF

    CALL newDADCD(J,i)
    IF(I.LT.0) then
       ipause=mypause(-666)
    endif
    a%r(I) = C
    GOTO 10


20  CONTINUE


    RETURN
    IF(PACKING)  THEN
       CALL NEWDAPACK(A )
    ENDIF
  END  subroutine  olddarea

  SUBROUTINE OLDDAREA77(A,IUNIT)
    implicit none
    integer ipause, mypauses, mypause
    CHARACTER(10) C10,K10
    integer  J(LNV),NOJOH,NVJOH,iin, iunit,i,ii,k,iche
    real(dp) c

    type (taylorlow) A

    if(iunit.eq.0) return

    if(checkass) then
       IF (.NOT. ASSOCIATED(A%R)) THEN
          line= " OLDDAREA77"
          ipause=mypauses(-997,line)
       ENDIF
    endif

    DO I=1,LNV
       J(I) = 0
    enddo

    CALL newDACLR(A)
    IF(PACKING)  THEN
       DO I=1,NM
          A%R(I)=zero
          a%yes(i)=.false.
       ENDDO
    ENDIF

    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10)') C10
    READ(IUNIT,'(A10,I6,A10,I6)') C10,NOJOH,K10,NVJOH

    IIN = 0

10  CONTINUE
    IIN = IIN + 1
    READ(IUNIT,*) II,C,(J(K),K=1,NVJOH)

    IF(II.LT.0) GOTO 20


    DO I=nv+1,NVJOH
       IF(J(I).NE.0) GOTO 10
    enddo
    iche=0
    DO I=1,nv
       iche=iche+j(i)
    enddo
    if(iche.gt.no) goto 10
    CALL newDADCD(J,i)
    IF(I.LT.0) then
       ipause=mypause(-666)
    endif
    a%r(I) = C
    GOTO 10
20  CONTINUE
    IF(PACKING)  THEN
       CALL NEWDAPACK(A )
    ENDIF
  END  subroutine  OLDDAREA77

end module newda
