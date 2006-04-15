!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size

module da_arrays
  use precision_constants
  use scratch_size
  implicit none
  public
  integer lda,lea,lia,lst
  integer,parameter::lno=200,lnv=100,lnomax=8,lnvmax=9,lstmax=800500,ldamax=16000,leamax=5000,liamax=50000
  logical(lp) :: reallocate = .true.
  logical(lp) :: notallocated = .true.
  logical(lp),parameter::etiennefix=.true.


  real(dp), allocatable,dimension(:)::cc                               ! lst
  integer, allocatable,dimension(:)::i1,i2                             ! lst
  integer, allocatable,dimension(:)::ie1,ie2,ieo                       !lea
  integer, allocatable,dimension(:)::ia1,ia2                           ! 0:lia
  integer, allocatable,dimension(:)::idano,idanv,idapo,idalm,idall     ! lda
  character(10), allocatable,dimension(:)::daname                      ! daname(lda)*10
  logical(lp), allocatable,dimension(:)::allvec                            !lda

  integer nda
  integer :: ndamaxi=0
  real(dp),TARGET :: total_da_size = c_300
  integer nst,nomax,nvmax,nmmax,nocut,lfi
  real(dp) facint(0:lno)
  integer nhole
  integer,TARGET :: lda_used =1000
  logical(lp),TARGET :: stable_da =.true.
  logical,TARGET :: check_da =.false.
  real(dp),target ::  da_absolute_aperture=c_1d6

contains

  subroutine alloc_all(no,nv,nd2t)
    implicit none
    integer no,nv,nd2t
    if(reallocate) then
       call dealloc_all
       call danum0(no,nv,nd2t)
       call alloc_
       notallocated=.false.
    endif
    if(notallocated) then
       lia=liamax
       lst=lstmax
       lda=ldamax
       lea=leamax
       notallocated=.false.
       call alloc_
    endif
  end subroutine alloc_all

  subroutine alloc_
    implicit none
    integer i
    allocate(cc(lst))
    allocate(i1(lst));allocate(i2(lst));
    allocate(ie1(lea));allocate(ie2(lea));allocate(ieo(lea));
    allocate(ia1(0:lia));allocate(ia2(0:lia));
    allocate(idano(lda));allocate(idanv(lda));allocate(idapo(lda));
    allocate(idalm(lda));allocate(idall(lda));
    allocate(daname(lda));
    allocate(allvec(lda));
    do i=1,lst
       cc(I)=zero  ! ADDED BY ETIENNE
       i1(i)=0
       i2(i)=0
    enddo
    do i=1,lea
       ie1(i)=0
       ie2(i)=0
       ieo(i)=0
    enddo
    do i=0,lia
       ia1(i)=0
       ia2(i)=0
    enddo
    do i=1,lda
       idano(i)=0
       idanv(i)=0
       idapo(i)=0
       idalm(i)=0
       idall(i)=0
    enddo
  end subroutine alloc_


  subroutine dealloc_all
    implicit none
    integer error,ipause,mypauses
    IF (ALLOCATED(cc)) THEN
       DEALLOCATE (cc, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(100," cc ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(i1)) THEN
       DEALLOCATE (i1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(101," i1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(i2)) THEN
       DEALLOCATE (i2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(102," i2 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(ie1)) THEN
       DEALLOCATE (ie1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(103," ie1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(ie2)) THEN
       DEALLOCATE (ie2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(104," ie2 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(ieo)) THEN
       DEALLOCATE (ieo, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(105," ieo ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(ia1)) THEN
       DEALLOCATE (ia1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(106," ia1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(ia2)) THEN
       DEALLOCATE (ia2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(107," ia2 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(idano)) THEN
       DEALLOCATE (idano, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(108," idano ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(idanv)) THEN
       DEALLOCATE (idanv, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(109," idanv ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(idapo)) THEN
       DEALLOCATE (idapo, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(110," idapo ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(idalm)) THEN
       DEALLOCATE (idalm, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(111," idalm ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(idall)) THEN
       DEALLOCATE (idall, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(112," idall ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(daname)) THEN
       DEALLOCATE (daname, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(112," daname ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(allvec)) THEN
       DEALLOCATE (allvec, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(112," allvec ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
  end subroutine dealloc_all

  subroutine danum(no,nv,numda)
    implicit none
    integer i,mm,no,numda,nv
    !     *****************************
    !
    !     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS OF
    !     ORDER NO AND NUMBER OF VARIABLES NV
    !
    numda = 1
    mm = max(nv,no)
    !
    do i=1,min(nv,no)
       numda = (numda*(mm+i))/i
    enddo
  end  subroutine danum

  subroutine danum0(no,nv,nd2t)
    use scratch_size
    use file_handler
    implicit none
    integer i,mm,no,nv,ldamin,nd2t ,nvt
    integer mf
    real(dp) size
    !     *****************************
    !
    !     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS OF
    !     ORDER NO AND NUMBER OF VARIABLES NV
    !
    lea = 1
    if(etiennefix) then
       nvt=nv
    else
       if(nv==1) then   ! still unknown feature of daini which must be taken into account
          nvt=2           ! perhaps related to Ying Wu modification in daini
       else
          nvt=nv
       endif
    endif
    mm = max(nvt,no)
    do i=1,min(nvt,no)
       lea = (lea*(mm+i))/i
    enddo
    ldamin=2+no + 2+ 4  !+ nd2t*(2+ndum)      ! (2+ no)daini  + varf +varc + assignmap(tpsalie.f90)


    ldamin=ldamin+2        !  2 + sum of ndumuser  in assign tpsa.f90


    lia=(no+1)**((nv+mod(nv,2))/2)
    lda=lda_used+ldamin
    lst=lda*lea

    size=(three*REAL(lea,kind=DP)+two*(REAL(lia,kind=DP)+one)+five*REAL(lda,kind=DP))*four/c_1024**2
    size=size+ten*REAL(lda,kind=DP)/c_1024**2+REAL(lda,kind=DP)/eight/c_1024**2
    size=size+REAL(lst,kind=DP)*(eight/c_1024**2+two*four/c_1024**2)
    if(size>total_da_size) then
       w_p=0
       w_p%nc=13
       w_p%fc='(12(1X,A72,/),(1X,A72))'
       write(w_p%c(1),'(a10,1x,i4,1x,i4)') " no,nv  = ",no,nv
       write(w_p%c(2),'(a10,1x,i8)') "    LEA = ",lea
       write(w_p%c(3),'(a24,1x,i8)') " ldamin (with nd2=6)  = ",ldamin
       write(w_p%c(4),'(a8,1x,i8)') " lia  = ",lia
       write(w_p%c(5),'(a8,1x,i8)') " lda  = ",lda
       write(w_p%c(6),'(a8,1x,i8)') " lst  = ",lst
       write(w_p%c(7),'(a14,1x,i8)') " ndamaxi    = ",ndamaxi
       write(w_p%c(8),'(a18,1x,g20.14)') " size in Mbytes = ",size
       write(w_p%c(9),'(a25,1x,g20.14)') " Total_da_size Allowed = ",Total_da_size
       w_p%c(10)=" "
       w_p%c(11)="************************************"
       w_p%c(12)="* Execution Continues Nevertheless *"
       w_p%c(13)="************************************"
       call write_e(1000)
       call kanalnummer(mf)
       open(unit=mf,file='too_big_da.txt')
       write(mf,*) "no,nv  = ",no,nv
       write(mf,*) "LEA = ",lea
       write(mf,*) "ldamin (with nd2=6)  = ",ldamin
       write(mf,*) "lia  = ",lia
       write(mf,*) "lda  = ",lda
       write(mf,*) "lst  = ",lst
       write(mf,*) "ndamaxi    = ",ndamaxi
       write(mf,*) "size in Mbytes = ",size
       write(mf,*) "Total_da_size Allowed = ",Total_da_size
       close(mf)
    endif


  end  subroutine danum0


end  module da_arrays
