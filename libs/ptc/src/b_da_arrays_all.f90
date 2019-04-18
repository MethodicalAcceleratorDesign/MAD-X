!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size

module da_arrays
  use precision_constants
  use my_own_1D_TPSA
  use scratch_size
  implicit none
  public
  integer lda,lea,lia,lst
  integer,private,parameter::nmax=400
  real(dp),private,parameter::tiny=1e-20_dp
  private ludcmp_nr,lubksb_nr
  ! johan
  ! integer,parameter::lno=2,lnv=6,lnomax=8,lnvmax=9,lstmax=800500,ldamax=16000,leamax=5000,liamax=50000
  !
  integer,parameter::lno=200,lnv=100,lnomax=8,lnvmax=9,lstmax=800500,ldamax=16000,leamax=5000,liamax=50000
  logical(lp) :: reallocate = .true.
  logical(lp) :: notallocated = .true.
  logical(lp),parameter::etiennefix=.true.


  real(dp), allocatable,dimension(:)::cc                               ! lst
  integer, allocatable,dimension(:)::i_1,i_2                             ! lst
  integer, allocatable,dimension(:)::ie1,ie2,ieo                       !lea
  integer, allocatable,dimension(:)::ia1,ia2                           ! 0:lia
  integer, allocatable,dimension(:)::idano,idanv,idapo,idalm,idall     ! lda
  character(10), allocatable,dimension(:)::daname                      ! daname(lda)*10
  logical(lp), allocatable,dimension(:)::allvec                            !lda

  integer nda_dab
  integer :: ndamaxi=0
  real(dp),TARGET :: total_da_size = 1.D38  !c_300
  integer nst0,nomax,nvmax,nmmax,nocut,lfi
  real(dp) facint(0:lno)
  integer nhole
  integer,TARGET :: lda_used =15000

contains

  subroutine alloc_all(no,nv)
    implicit none
    integer no,nv
    if(reallocate) then
       call dealloc_all
       call danum0(no,nv)
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
allocate(cc(lst))
allocate(i_1(lst));allocate(i_2(lst));
allocate(ie1(lea));allocate(ie2(lea));allocate(ieo(lea));
allocate(ia1(0:lia));allocate(ia2(0:lia));
allocate(idano(lda));allocate(idanv(lda));allocate(idapo(lda));
allocate(idalm(lda));allocate(idall(lda));
allocate(daname(lda));
allocate(allvec(lda));

cc=0.0_dp ! ADDED BY ETIENNE
i_1=0
i_2=0
ie1=0
ie2=0
ieo=0
ia1=0
ia2=0
idano=0
idanv=0
idapo=0
idalm=0
idall=0

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
    IF (ALLOCATED(i_1)) THEN
       DEALLOCATE (i_1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(101," i_1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(i_2)) THEN
       DEALLOCATE (i_2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(102," i_2 ARRAY not DEALLOCATED : PROBLEMS")
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

  subroutine danum0(no,nv)
    use scratch_size
    use file_handler
    implicit none
    integer i,mm,no,nv,ldamin ,nvt
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

    if(lingyun_yang) then
       size=REAL(lst,kind=DP)*(8.0_dp/1024.0_dp**2+2.0_dp*4.0_dp/1024.0_dp**2)
    else
       size=(3.0_dp*REAL(lea,kind=DP)+2.0_dp*(REAL(lia,kind=DP)+1.0_dp)+5.0_dp*REAL(lda,kind=DP))*4.0_dp/1024.0_dp**2
       size=size+10.0_dp*REAL(lda,kind=DP)/1024.0_dp**2+REAL(lda,kind=DP)/8.0_dp/1024.0_dp**2
       size=size+REAL(lst,kind=DP)*(8.0_dp/1024.0_dp**2+2.0_dp*4.0_dp/1024.0_dp**2)
    endif
    if(size>total_da_size.or.printdainfo) then
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
       write(6,*) "no,nv  = ",no,nv
       write(6,*) "LEA = ",lea
       write(6,*) "ldamin (with nd2=6)  = ",ldamin
       write(6,*) "lia  = ",lia
       write(6,*) "lda  = ",lda
       write(6,*) "lst  = ",lst
       write(6,*) "ndamaxi    = ",ndamaxi
       write(6,*) "size in Mbytes = ",size
       write(6,*) "Total_da_size Allowed = ",Total_da_size
       close(mf)
    endif


  end  subroutine danum0

!!!!! some da stuff

 
  subroutine matinv(a,ai,n,nmx,ier)

    implicit none

    integer i,ier,j,n,nmx
    integer,dimension(nmax)::indx
    real(dp) d
    real(dp),dimension(nmx,nmx)::a,ai
    real(dp),dimension(nmax,nmax)::aw
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif

    aw(1:n,1:n) = a(1:n,1:n)

    call ludcmp_nr(aw,n,nmax,indx,d,ier)
    if (ier .eq. 132) return

    ai(1:n,1:n) = 0.0_dp
    !    forall (i = 1:n) ai(i,i) = one
    do i=1,n
       ai(i,i) = 1.0_dp
    enddo

    do j=1,n
       call lubksb_nr(aw,n,nmax,indx,ai(1,j),nmx)
    enddo

  end subroutine matinv


  !
  subroutine ludcmp_nr(a,n,np,indx,d,ier)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE DECOMPOSES A MATRIX INTO LU FORMAT
    !     INPUT A: NXN MATRIX - WILL BE OVERWRITTEN BY THE LU DECOMP.
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 35
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ier,imax,j,k,n,np
    integer,dimension(np)::indx
    real(dp) aamax,d,dum,sum
    real(dp),dimension(np,np)::a
    real(dp),dimension(nmax)::vv
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ier=0
    d=1.0_dp
    do i=1,n
       aamax=0.0_dp
       do j=1,n
          if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       enddo
       if(aamax.eq.0.0_dp) then
          ier=132
          return
       endif
       vv(i)=1.0_dp/aamax
    enddo
    do j=1,n
       if(j.gt.1) then
          do i=1,j-1
             sum=a(i,j)
             if(i.gt.1) then
                do k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                enddo
                a(i,j)=sum
             endif
          enddo
       endif
       aamax=0.0_dp
       do i=j,n
          sum=a(i,j)
          if (j.gt.1) then
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if(dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j.ne.imax) then
          do k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if(j.ne.n) then
          if(a(j,j).eq.0.0_dp) a(j,j)=tiny
          dum=1.0_dp/a(j,j)
          do i=j+1,n
             a(i,j)=a(i,j)*dum
          enddo
       endif
    enddo
    if(a(n,n).eq.0.0_dp) a(n,n)=tiny
    return
  end subroutine ludcmp_nr
  !
  subroutine lubksb_nr(a,n,np,indx,b,nmx)
    implicit none
    !     ************************************
    !
    !     THIS SUBROUTINE SOLVES SET OF LINEAR EQUATIONS AX=B,
    !     INPUT A: NXN MATRIX IN lu FORM GIVEN BY ludcmp_nr
    !           NP: PHYSICAL DIMENSION OF A
    !           INDX: ROW PERMUTATION VECTOR
    !           D: EVEN OR ODD ROW INTERCHANGES
    !           B: RHS OF LINEAR EQUATION - WILL BE OVERWRITTEN BY X
    !
    !     REFERENCE: NUMERICAL RECIPIES BY PRESS ET AL (CAMBRIDGE) PG. 36
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,j,ll,n,nmx,np
    integer,dimension(np)::indx
    real(dp) sum
    real(dp),dimension(np,np)::a
    real(dp),dimension(nmx)::b
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    ii = 0
    do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if(ii.ne.0) then
          do j=ii,i-1
             sum = sum-a(i,j)*b(j)
          enddo
       else if (sum.ne.0.0_dp) then
          ii = i
       endif
       b(i)=sum
    enddo
    do i=n,1,-1
       sum=b(i)
       if(i.lt.n) then
          do j=i+1,n
             sum = sum-a(i,j)*b(j)
          enddo
       endif

       b(i)=sum/a(i,i)

    enddo
    return
  end subroutine lubksb_nr
  !


end  module da_arrays
