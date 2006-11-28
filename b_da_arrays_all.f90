!The Full Polymorphic Package
!Copyright (C) Etienne Forest
! Based on an original Fortran77 prototype
! developed

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
  character(120) line

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

  subroutine daallno(ic,l,ccc)
    implicit none
    !     ********************************
    !
    !     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
    !     ORDER NOmax AND NUMBER OF VARIABLES NVmax
    !
    !-----------------------------------------------------------------------------
    !
    logical(lp) incnda
    integer i,ind,l,ndanum,no,nv,ipause,mypauses
    integer,dimension(:)::ic
    real(dp) x
    character(10) c,ccc
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in daallno ", sqrt(crash)
       endif
       return
    endif
    !
    no=nomax
    nv=nvmax
    ind = 1
    do i=1,l
       if(ic(i).gt.0.and.ic(i).le.nda) then
          !         DANAME(IC(I)) = C
          !         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
       else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
             write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
                  &' NOMAX, NVMAX = ',nomax,nvmax
             ipause=mypauses(5,line)
             call dadeb(31,'ERR DAALL ',1)
          endif
          !
          if(nhole.gt.0) then
             ind=nda
20           if (allvec(ind)) then
                ind = ind - 1
                goto 20
             endif
             incnda = .false.
             nhole=nhole-1
          else
             incnda = .true.
             nda = nda + 1
             ind=nda
             if(nda.gt.lda) then
                write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
                ipause=mypauses(6,line)
                call dadeb(31,'ERR DAALL ',1)
             endif
          endif

          allvec(ind) = .true.
          ic(i) = ind
          !
          if(nv.ne.0) then
             call danum(no,nv,ndanum)
          else
             ndanum = no
          endif

          c = ccc
          if(l.ne.1) write(c(6:10),'(I5)') i
          daname(ind) = c

          if (incnda) then
             if(ind.gt.nomax+2) then
                idano(ind) = nomax
                idanv(ind) = nvmax
                idapo(ind) = nst + 1
                idalm(ind) = nmmax
                idall(ind) = 0
                nst = nst + nmmax
             else
                idano(ind) = no
                idanv(ind) = nv
                idapo(ind) = nst + 1
                idalm(ind) = ndanum
                idall(ind) = 0
                nst = nst + ndanum
             endif
          endif
          !
          if(nst.gt.lst) then
             x=-one
             w_p=0
             w_p%nc=5
             w_p%fc='(4(1X,a72,/),(1X,a72))'
             w_p%c(1)= 'ERROR IN DAALL, STACK EXHAUSTED '
             w_p%c(2)=  ' NST,LST '
             write(w_p%c(3),'(i8,1x,i8)') NST,LST
             w_p%c(4)=  ' NDA,NDANUM,NDA*NDANUM '
             write(w_p%c(5),'(i8,1x,i8,1x,i8)') nda,ndanum,nda*ndanum
             CALL WRITE_E(124)
             call dadeb(31,'ERR DAALL ',1)
          endif
          !
          if(nv.eq.0.or.nomax.eq.1) then
             call daclr(ic(i))
             idall(ic(i)) = idalm(ic(i))
          endif
       endif
    enddo
    !
    if(nda.gt.ndamaxi) ndamaxi=nda

    return
  end subroutine daallno

  subroutine daallno1(ic,ccc)
    implicit none
    !     ********************************
    !
    !     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
    !     ORDER NOmax AND NUMBER OF VARIABLES NVmax
    !
    !-----------------------------------------------------------------------------
    !
    logical(lp) incnda
    integer ind,ndanum,no,nv,ic,ipause,mypauses
    character(10) c,ccc
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in daallno1 ", sqrt(crash)
       endif
       return
    endif
    !
    no=nomax
    nv=nvmax
    ind = 1
    if(ic.gt.0.and.ic.le.nda) then
       !         DANAME(IC(I)) = C
       !         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' NOMAX, NVMAX = ',nomax,nvmax
          ipause=mypauses(7,line)
          call dadeb(31,'ERR DAALL ',1)
       endif
       !
       if(nhole.gt.0) then
          ind=nda
20        if (allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          nhole=nhole-1
       else
          incnda = .true.
          nda = nda + 1
          ind=nda
          if(nda.gt.lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(8,line)
             call dadeb(31,'ERR DAALL ',1)
          endif
       endif

       allvec(ind) = .true.
       ic = ind
       !
       if(nv.ne.0) then
          call danum(no,nv,ndanum)
       else
          ndanum = no
       endif

       c = ccc
       write(c(6:10),'(I5)') 1
       daname(ind) = c

       if (incnda) then
          if(ind.gt.nomax+2) then
             idano(ind) = nomax
             idanv(ind) = nvmax
             idapo(ind) = nst + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst = nst + nmmax
          else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst = nst + ndanum
          endif
       endif
       !
       if(nst.gt.lst) then
          w_p=0
          w_p%nc=5
          w_p%fc='(4(1X,a72,/),(1X,a72))'
          w_p%c(1)= 'ERROR IN DAALL, STACK EXHAUSTED '
          w_p%c(2)=  ' NST,LST '
          write(w_p%c(3),'(i8,1x,i8)') NST,LST
          w_p%c(4)=  ' NDA,NDANUM,NDA*NDANUM '
          write(w_p%c(5),'(i8,1x,i8,1x,i8)') nda,ndanum,nda*ndanum
          CALL WRITE_E(125)
          call dadeb(31,'ERR DAALL ',1)
       endif
       !
       if(nv.eq.0.or.nomax.eq.1) then
          call daclr(ic)
          idall(ic) = idalm(ic)
       endif
    endif

    !
    if(nda.gt.ndamaxi) ndamaxi=nda

    return
  end subroutine daallno1

  subroutine daall(ic,l,ccc,no,nv)
    implicit none
    !     ********************************
    !
    !     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
    !     ORDER NO AND NUMBER OF VARIABLES NV
    !
    !-----------------------------------------------------------------------------
    !
    logical(lp) incnda
    integer i,ind,l,ndanum,no,nv,ipause,mypauses
    integer,dimension(:)::ic
    character(10) c,ccc
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in daall ", sqrt(crash)
       endif
       return
    endif
    !
    ind = 1

    do i=1,l
       if(ic(i).gt.0.and.ic(i).le.nda) then
          !         DANAME(IC(I)) = C
          !         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
       else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
             write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
                  &' NOMAX, NVMAX = ',nomax,nvmax
             ipause=mypauses(9,line)
             call dadeb(31,'ERR DAALL ',1)
          endif
          !
          if(nhole.gt.0) then
             ind=nda
20           if (allvec(ind)) then
                ind = ind - 1
                goto 20
             endif
             incnda = .false.
             nhole=nhole-1
          else
             incnda = .true.
             nda = nda + 1
             ind=nda
             if(nda.gt.lda) then
                write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
                ipause=mypauses(10,line)
                call dadeb(31,'ERR DAALL ',1)
             endif
          endif

          allvec(ind) = .true.

          ic(i) = ind
          !
          if(nv.ne.0) then
             call danum(no,nv,ndanum)
          else
             ndanum = no
          endif

          c = ccc
          if(l.ne.1) write(c(6:10),'(I5)') i

          daname(ind) = c

          if (incnda) then
             if(ind.gt.nomax+2) then
                idano(ind) = nomax
                idanv(ind) = nvmax
                idapo(ind) = nst + 1
                idalm(ind) = nmmax
                idall(ind) = 0
                nst = nst + nmmax
             else
                idano(ind) = no
                idanv(ind) = nv
                idapo(ind) = nst + 1
                idalm(ind) = ndanum
                idall(ind) = 0
                nst = nst + ndanum
             endif
          endif
          !
          if(nst.gt.lst) then
             w_p=0
             w_p%nc=5
             w_p%fc='(4(1X,a72,/),(1X,a72))'
             w_p%c(1)= 'ERROR IN DAALL, STACK EXHAUSTED '
             w_p%c(2)=  ' NST,LST '
             write(w_p%c(3),'(i8,1x,i8)') NST,LST
             w_p%c(4)=  ' NDA,NDANUM,NDA*NDANUM '
             write(w_p%c(5),'(i8,1x,i8,1x,i8)') nda,ndanum,nda*ndanum
             CALL WRITE_E(126)
             call dadeb(31,'ERR DAALL ',1)
          endif
          !
          !          IF(NV.EQ.0) THEN
          if(nv.eq.0.or.nomax.eq.1) then
             call daclr(ic(i))
             idall(ic(i)) = idalm(ic(i))
          endif
       endif
    enddo
    !
    if(nda.gt.ndamaxi) ndamaxi=nda

    return
  end subroutine daall


  subroutine daall1(ic,ccc,no,nv)
    implicit none
    !     ********************************
    !
    !     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
    !     ORDER NO AND NUMBER OF VARIABLES NV
    !
    !-----------------------------------------------------------------------------
    !
    logical(lp) incnda
    integer ic,ind,ndanum,no,nv,ipause,mypauses
    character(10) c,ccc
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in daall1 ", sqrt(crash)
       endif
       return
    endif
    !
    ind = 1

    if(ic.gt.0.and.ic.le.nda) then
       !         DANAME(ic) = C
       !         IF(IDANO(ic).EQ.NO.AND.IDANV(ic).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' NOMAX, NVMAX = ',nomax,nvmax
          ipause=mypauses(11,line)
          call dadeb(31,'ERR DAALL ',1)
       endif
       !
       if(nhole.gt.0) then
          ind=nda
20        if (allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          nhole=nhole-1
       else
          incnda = .true.
          nda = nda + 1
          ind=nda
          if(nda.gt.lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(12,line)
             call dadeb(31,'ERR DAALL ',1)
          endif
       endif

       allvec(ind) = .true.

       ic = ind
       !
       if(nv.ne.0) then
          call danum(no,nv,ndanum)
       else
          ndanum = no
       endif

       c = ccc
       write(c(6:10),'(I5)') 1

       daname(ind) = c

       if (incnda) then
          if(ind.gt.nomax+2) then
             idano(ind) = nomax
             idanv(ind) = nvmax
             idapo(ind) = nst + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst = nst + nmmax
          else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst = nst + ndanum
          endif
       endif
       !
       if(nst.gt.lst) then
          w_p=0
          w_p%nc=5
          w_p%fc='(4(1X,a72,/),(1X,a72))'
          w_p%c(1)= 'ERROR IN DAALL, STACK EXHAUSTED '
          w_p%c(2)=  ' NST,LST '
          write(w_p%c(3),'(i8,1x,i8)') NST,LST
          w_p%c(4)=  ' NDA,NDANUM,NDA*NDANUM '
          write(w_p%c(5),'(i8,1x,i8,1x,i8)') nda,ndanum,nda*ndanum
          CALL WRITE_E(127)
          call dadeb(31,'ERR DAALL ',1)
       endif
       !
       !          IF(NV.EQ.0) THEN
       if(nv.eq.0.or.nomax.eq.1) then
          call daclr(ic)
          idall(ic) = idalm(ic)
       endif
    endif
    !
    if(nda.gt.ndamaxi) ndamaxi=nda

    return
  end subroutine daall1
  !
  subroutine dadal(idal,l)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
    !
    !-----------------------------------------------------------------------------
    !
    integer i,l,ipause,mypauses
    integer,dimension(:)::idal
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dadal ", sqrt(crash)
       endif
       return
    endif
    do i=l,1,-1
       if(idal(i).le.nomax+2.or.idal(i).gt.nda) then
          write(line,'(a38,i8,1x,i8)') 'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),nda
          ipause=mypauses(13,line)
          call dadeb(31,'ERR DADAL ',1)
       endif
       if(idal(i).eq.nda) then
          !       deallocate
          nst = idapo(nda) - 1
          nda = nda - 1
       else
          nhole=nhole+1
       endif

       allvec(idal(i)) = .false.

       !        IDANO(IDAL(I)) = 0
       !        IDANV(IDAL(I)) = 0
       !        IDAPO(IDAL(I)) = 0
       !        IDALM(IDAL(I)) = 0
       idall(idal(i)) = 0

       idal(i) = 0
    enddo

    return
  end subroutine dadal

  subroutine dadal1(idal)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
    !
    !-----------------------------------------------------------------------------
    !
    integer idal,ipause,mypauses
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dadal1 ", sqrt(crash)
       endif
       return
    endif
    if(idal.le.nomax+2.or.idal.gt.nda) then
       write(line,'(a35,i8,1x,i8)') 'ERROR IN ROUTINE DADAL, IDAL,NDA = ',idal,nda
       ipause=mypauses(14,line)
       call dadeb(31,'ERR DADAL ',1)
    endif
    if(idal.eq.nda) then
       !       deallocate
       nst = idapo(nda) - 1
       nda = nda - 1
    else
       nhole=nhole+1
    endif

    allvec(idal) = .false.

    !        IDANO(IDAL(I)) = 0
    !        IDANV(IDAL(I)) = 0
    !        IDAPO(IDAL(I)) = 0
    !        IDALM(IDAL(I)) = 0
    idall(idal) = 0

    idal = 0

    return
  end subroutine dadal1

  subroutine count_da(n)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE counts allocate da
    !
    !-----------------------------------------------------------------------------
    !
    integer i,n
    !
    n=0
    do i=1,lda
       if(allvec(i)) n=n+1
    enddo
    return
  end subroutine count_da

  subroutine dadeb(iunit,c,istop)
    implicit none
    !     *******************************
    !
    !     THIS SUBROUTINE SERVES AS A DEBUGGING TOOL. IT PRINTS ALL
    !     NONZERO INFORMATION IN THE COMMON BLOCKS AND ALL DA  VECTORS.
    !
    !-----------------------------------------------------------------------------
    !
    integer istop,iunit,I
    !integer,dimension(0:1)::i8
    character(10) c
    C_%STABLE_DA=.false.
    return
  end subroutine dadeb

  subroutine daclr(inc)
    implicit none
    !     *********************
    !
    !     THIS SUBROUTINE SETS ALL THE STACK SPACE RESERVED FOR VARIABLE
    !     C TO ZERO
    !
    !-----------------------------------------------------------------------------
    !
    integer i,illc,ilmc,inc,inoc,invc,ipoc
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in daclr ", sqrt(crash)
       endif
       return
    endif
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in daclr ", sqrt(crash)
       endif
       return
    endif
    !
    do i=ipoc,ipoc+ilmc-1
       !
       cc(i) = zero
       !
    enddo
    !
    return
  end subroutine daclr
  !
  subroutine dainf(inc,inoc,invc,ipoc,ilmc,illc)
    implicit none
    !     **********************************************
    !
    !     THIS SUBROUTINE SEARCHES THE NUMBER OF DA VECTOR C
    !     AND RETURS THE INFORMATION IN COMMON DA
    !
    !-----------------------------------------------------------------------------
    !
    integer illc,ilmc,inc,inoc,invc,ipoc,ipause,mypauses
    !
    if(inc.ge.1.and.inc.le.nda) then
       inoc = idano(inc)
       invc = idanv(inc)
       ipoc = idapo(inc)
       ilmc = idalm(inc)
       illc = idall(inc)
       return
    endif
    !
    write(line,'(a26,1x,i8,1x,a11)')  'ERROR IN DAINF, DA VECTOR ',inc,' NOT FOUND '
    ipause=mypauses(35,line)
    call dadeb(31,'ERR DAINF ',1)
    !
    return
  end subroutine dainf
end  module da_arrays
