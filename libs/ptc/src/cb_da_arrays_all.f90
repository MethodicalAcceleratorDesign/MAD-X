!The Full Polymorphic Package
!Copyright (C) Etienne Forest and Frank Schmidt
! See file a_scratch_size

module c_da_arrays
  use precision_constants
  use my_own_1D_TPSA
  use scratch_size
  implicit none
  public
  private ludcmp_nr,lubksb_nr
  integer c_lda,c_lea,c_lia,c_lst
  integer,private,parameter::nmax=400
  real(dp),private,parameter::tiny=1e-20_dp
 
  ! johan
  ! integer,parameter::c_lno=2,c_lnv=6,c_lc_nomax=8,c_lc_nvmax=9,c_lstmax=800500,c_ldamax=16000,c_leamax=5000,c_liamax=50000
  !
  integer,parameter::c_lno=200,c_lnv=100,c_lc_nomax=8,c_lc_nvmax=9,c_lstmax=800500,c_ldamax=16000,c_leamax=5000,c_liamax=50000
  logical(lp) :: c_reallocate = .true.
  logical(lp) :: c_notallocated = .true.
  logical(lp),parameter::c_etiennefix=.true.


  complex(dp), allocatable,dimension(:)::c_cc                               ! c_lst
  integer, allocatable,dimension(:)::c_i_1,c_i_2                             ! c_lst
  integer, allocatable,dimension(:)::c_ie1,c_ie2,c_ieo                       !c_lea
  integer, allocatable,dimension(:)::c_ia1,c_ia2                           ! 0:c_lia
  integer, allocatable,dimension(:)::c_idano,c_idanv,c_idapo,c_idalm,c_idall     ! c_lda
  character(10), allocatable,dimension(:)::c_daname                      ! c_daname(c_lda)*10
  logical(lp), allocatable,dimension(:)::c_allvec                            !c_lda

  integer c_nda_dab
  integer :: c_ndamaxi=0
  real(dp),TARGET :: c_total_da_size = 1.D38  !c_300
  integer c_nst0,c_nomax,c_nvmax,c_nmmax,c_nocut,c_lfi
 ! real(dp) c_facint(0:c_lno)
  integer c_nhole
  integer,TARGET :: c_lda_used =15000

contains

  subroutine alloc_all_c(no,nv)
    implicit none
    integer no,nv
    if(c_reallocate) then
       call dealloc_all_c
       call danum0_c(no,nv)
       call alloc_c
       c_notallocated=.false.
    endif
    if(c_notallocated) then
       c_lia=c_liamax
       c_lst=c_lstmax
       c_lda=c_ldamax
       c_lea=c_leamax
       c_notallocated=.false.
       call alloc_c
    endif
  end subroutine alloc_all_c

  subroutine alloc_c
    implicit none
    allocate(c_cc(c_lst))
    allocate(c_i_1(c_lst));allocate(c_i_2(c_lst));
    allocate(c_ie1(c_lea));allocate(c_ie2(c_lea));allocate(c_ieo(c_lea));
    allocate(c_ia1(0:c_lia));allocate(c_ia2(0:c_lia));
    allocate(c_idano(c_lda));allocate(c_idanv(c_lda));allocate(c_idapo(c_lda));
    allocate(c_idalm(c_lda));allocate(c_idall(c_lda));
    allocate(c_daname(c_lda));
    allocate(c_allvec(c_lda));
!    do i=1,c_lst
!       c_cc(I)=0.0_dp  ! ADDED BY ETIENNE
!       c_i_1(i)=0
!       c_i_2(i)=0
!    enddo
!    do i=1,c_lea
!       c_ie1(i)=0
!       c_ie2(i)=0
!       c_ieo(i)=0
!    enddo
!    do i=0,c_lia
!       c_ia1(i)=0
!       c_ia2(i)=0
!    enddo
!    do i=1,c_lda
!       c_idano(i)=0
!       c_idanv(i)=0
!       c_idapo(i)=0
!       c_idalm(i)=0
!       c_idall(i)=0
!    enddo
c_cc=0.0_dp ! ADDED BY ETIENNE
c_i_1=0
c_i_2=0
c_ie1=0
c_ie2=0
c_ieo=0
c_ia1=0
c_ia2=0
c_idano=0
c_idanv=0
c_idapo=0
c_idalm=0
c_idall=0
  end subroutine alloc_c


  subroutine dealloc_all_c
    implicit none
    integer error,ipause,mypauses
    IF (ALLOCATED(c_cc)) THEN
       DEALLOCATE (c_cc, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(100," c_cc ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_i_1)) THEN
       DEALLOCATE (c_i_1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(101," c_i_1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_i_2)) THEN
       DEALLOCATE (c_i_2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(102," c_i_2 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_ie1)) THEN
       DEALLOCATE (c_ie1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(103," c_ie1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_ie2)) THEN
       DEALLOCATE (c_ie2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(104," c_ie2 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_ieo)) THEN
       DEALLOCATE (c_ieo, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(105," c_ieo ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_ia1)) THEN
       DEALLOCATE (c_ia1, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(106," c_ia1 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_ia2)) THEN
       DEALLOCATE (c_ia2, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(107," c_ia2 ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_idano)) THEN
       DEALLOCATE (c_idano, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(108," c_idano ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_idanv)) THEN
       DEALLOCATE (c_idanv, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(109," c_idanv ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_idapo)) THEN
       DEALLOCATE (c_idapo, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(110," c_idapo ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_idalm)) THEN
       DEALLOCATE (c_idalm, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(111," c_idalm ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_idall)) THEN
       DEALLOCATE (c_idall, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(112," c_idall ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_daname)) THEN
       DEALLOCATE (c_daname, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(112," c_daname ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
    IF (ALLOCATED(c_allvec)) THEN
       DEALLOCATE (c_allvec, STAT = error)
       IF(ERROR==0) THEN
       ELSE
          ipause=mypauses(112," c_allvec ARRAY not DEALLOCATED : PROBLEMS")
       ENDIF
    ENDIF
  end subroutine dealloc_all_c

  subroutine danum_c(no,nv,numda)
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
  end  subroutine danum_c

  subroutine danum0_C(no,nv)
    use scratch_size
    use file_handler
    implicit none
    integer i,mm,no,nv,c_ldamin ,nvt
    integer mf
    real(dp) size
    !     *****************************
    !
    !     THIS SUBROUTINE COMPUTES THE NUMBER OF MONOMIALS OF
    !     ORDER NO AND NUMBER OF VARIABLES NV
    !
    c_lea = 1
    if(c_etiennefix) then
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
       c_lea = (c_lea*(mm+i))/i
    enddo
    c_ldamin=2+no + 2+ 4  !+ nd2t*(2+ndum)      ! (2+ no)daini  + varf +varc + assignmap(tpsalie.f90)


    c_ldamin=c_ldamin+2        !  2 + sum of ndumuser  in assign tpsa.f90


    c_lia=(no+1)**((nv+mod(nv,2))/2)
    c_lda=c_lda_used+c_ldamin
    c_lst=c_lda*c_lea


       size=(3.0_dp*REAL(c_lea,kind=DP)+2.0_dp*(REAL(c_lia,kind=DP)+1.0_dp)+5.0_dp*REAL(c_lda,kind=DP))*4.0_dp/1024.0_dp**2
       size=size+10.0_dp*REAL(c_lda,kind=DP)/1024.0_dp**2+REAL(c_lda,kind=DP)/8.0_dp/1024.0_dp**2
       size=size+REAL(c_lst,kind=DP)*(8.0_dp/1024.0_dp**2+2.0_dp*4.0_dp/1024.0_dp**2)

    if(size>c_total_da_size.or.printdainfo) then
       call kanalnummer(mf)
       open(unit=mf,file='too_big_da.txt')
       write(mf,*) "no,nv  = ",no,nv
       write(mf,*) "c_lea = ",c_lea
       write(mf,*) "c_ldamin (with nd2=6)  = ",c_ldamin
       write(mf,*) "c_lia  = ",c_lia
       write(mf,*) "c_lda  = ",c_lda
       write(mf,*) "c_lst  = ",c_lst
       write(mf,*) "c_ndamaxi    = ",c_ndamaxi
       write(mf,*) "size in Mbytes = ",size
       write(mf,*) "c_total_da_size Allowed = ",c_total_da_size
       write(6,*) "no,nv  = ",no,nv
       write(6,*) "c_lea = ",c_lea
       write(6,*) "c_ldamin (with nd2=6)  = ",c_ldamin
       write(6,*) "c_lia  = ",c_lia
       write(6,*) "c_lda  = ",c_lda
       write(6,*) "c_lst  = ",c_lst
       write(6,*) "c_ndamaxi    = ",c_ndamaxi
       write(6,*) "size in Mbytes = ",size
       write(6,*) "c_total_da_size Allowed = ",c_total_da_size
       close(mf)
    endif


  end  subroutine danum0_C



  subroutine c_matinv(a,ai,n,nmx,ier)

    implicit none

    integer i,ier,j,n,nmx
    integer,dimension(nmax)::indx
    complex(dp) d
    complex(dp),dimension(nmx,nmx)::a,ai
    complex(dp),dimension(nmax,nmax)::aw
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

  end subroutine c_matinv


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
    complex(dp) d,dum,sum
    real(dp) aamax
    complex(dp),dimension(np,np)::a
    complex(dp),dimension(nmax)::vv
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
          if(abs(dum).ge.aamax) then
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
    complex(dp) sum
    complex(dp),dimension(np,np)::a
    complex(dp),dimension(nmx)::b
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
       else if (abs(sum).ne.0.0_dp) then
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


  




end  module c_da_arrays
