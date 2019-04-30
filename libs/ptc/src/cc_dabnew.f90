! The Full Polymorphic Package
! The module in this file is, to the best of our knowledge,
! the property of Lawrence Berkeley National Laboratory
! Its distribution and commercial usage may therefore be governed by the laws of the
! United States of America
module c_dabnew
  use c_da_arrays,eeps=>eps
  implicit none
  private
!  public
  private daalc_lno1,daall,damult,dasqrt,dacmut,dacma,DALINt,dacctt
  private dainvt,dapint,dadert,dacfut
  private dadeb,dapac,dachk,damch,dadcd,dancd,hash,dehash
  public C_STABLE_DA,C_watch_user,c_daini,c_daall0,c_davar,c_damul,c_dacycle
  public c_dapok,c_dapek,c_etall1,c_DAshift,c_dacop,c_dacon,c_daadd,c_dacad,c_dasub
  public c_dacsu,c_dasuc,c_dacmu,c_dadal1,c_dapri,c_dapri77,c_darea,c_darea77,c_daeps
  public c_dacdi,c_dadic,c_count_da,c_mtree,c_dafun,c_DAABS,c_dadiv,c_take,c_datrunc,c_dader,c_datra
  public c_daran,c_dacfu,c_matinv,c_dapek0,c_dapok0,c_dacct,c_dainv,c_etcom,c_danot
  public c_print_eps
  integer,private,parameter:: lsw=1
  integer :: c_lda_max_used=0
  ! integer,private,parameter::nmax=400,lsw=1
  ! real(dp),private,parameter::tiny=c_1d_20
  character(120),private :: line
  logical(lp) C_STABLE_DA,C_watch_user,C_check_stable
  real(dp), private :: eps=1.d-38,epsprint=1.d-38
!real(dp),public :: eps_clean=0
complex(dp), private :: i_=(0.0_dp,1.0_dp)
contains
  !******************************************************************************
  !                                                                             *
  !                                                                             *
  !                                                                             *
  !               DIFFERENTIAL ALGEBRA PACKAGE OF M. BERZ                       *
  !                     ****************************                            *
  !                               holy3 complex                                 *
  !                                                                             *
  !                                                                             *
  !                                                                             *
  !         VERSION FOR MACHINE IN LINE THAT IS NOT COMMENTED OFF               *
  !        TO CREATE DIFFERENT VERSIONS, USE THE PROGRAM 'VERSION'              *
  !                                                                             *
  !                                                                             *
  !                                                                             *
  !                                                                             *
  !        THIS PACKAGE WAS INITIALLY WRITTEN BY PROF. M. BERZ WHILE AT         *
  !        THE LAWRENCE BERKELEY LABORATORY.                                    *
  !        IT HAS BEEN EXTENSIVELY MODIFIED BY THE MEMBERS OF THE ESG GROUP.    *
  !        THEREFORE PROF. BERZ SHOULD NOT BE HELD RESPONSIBLE FOR ANY BUGS.    *
  !                                                                             *
  !                  NEW RULES OF THE GAME (EXHAUSTIVE)                         *
  !                 **********************************                          *
  !                         THERE ARE NONE                                      *
  !                                                                             *
  !******************************************************************************
  !
  !
  !     THIS FILE CONTAINS ROUTINES TO PERFORM DIFFERENTIAL ALGEBRA (DA)
  !     AS AN OPTION, ALSO COMPONENTWISE ALGEBRA (CA) CAN BE PERFORMED.
  !     A DESCRIPTION OF THE INTERNAL ARRAYS USED BY THE ROUTINES CAN
  !     BE FOUND IN BLOCKDATA DABLD.
  !
  !
  !     SHORT REFERENCE CHART
  !     *********************
  !
  !     THE PARAMETERS USED BELOW HAVE THE FOLLOWING MEANING:
  !
  !     A,B:                NAME OF INPUT DA VECTORS   (INTEGER)
  !     C:                  NAME OF OUTPUT DA VECTOR   (INTEGER)
  !     X,Y:                NAME OF INPUT DA MATRIX    (INTEGER(...))
  !     Z:                  NAME OF OUTPUT DA MATRIX   (INTEGER(...))
  !
  !     F:                  NAME OF A DA FUNCTION      (CHARACTER(4))
  !     G:                  NAME OF EXTERNAL FUNCTION  (real(dp))
  !     JJ:                 ARRAY OF EXPONENTS         (INTEGER(20))
  !     O:                  ORDER                      (INTEGER)
  !     N:                  NUMBER OF VARIABLES        (INTEGER)
  !     I,J,K:              INTEGER NUMBER             (INTEGER
  !     R,RA,RB:            REAL NUMBERS               (real(dp))
  !     H:                  ARRAY OF LENGTH LH         (real(dp))
  !     U:                  OUTPUT UNIT NUMBER         (INTEGER)
  !     T:                  COMMENT TEXT               (CHARACTER(10))
  !
  !
  !               SUBROUTINES AND THEIR CALLING PARAMETERS
  !               ****************************************
  !
  !     c_daini(O,N,U):       INITIALIZES CONTROL ARRAYS AND SETS MAX. ORDER O AND
  !                         MAX. NUMBER OF VARIABLES N. MUST BE CALLED BEFORE ANY
  !                         OTHER DA ROUTINE CAN BE USED.
  !
  !     DAALL(A,I,T,O,N):   ALLOCATES SPACE FOR I VECTORS A. T: CHARACTER NAME
  !     DADAL(A,I):         DEALLOCATES THE I VECTORS A.
  !!     c_davar(A,R,I):       MAKES A INDEPENDENT VARIABLE # I WITH INITIAL VALUE R
  !!     c_dacon(A,R):         SETS A TO CONSTANT R
  !     DANOT(O):           SETS NEW TRUNCATION ORDER O FOR DA OPERATIONS
  !     DAEPS(R):           SETS NEW ZERO TOLERANCE EPSILON
  !
  !!     c_dapek(A,JJ,R):      RETURNS COEF R OF MONOMIAL WITH EXPONENTS JJ OF A
  !!     c_dapok(A,JJ,R):      SETS COEF OF MONOMIAL WITH EXPONENTS JJ OF A TO R
  !
  !!     c_dacop(A,C):         PERFORMS C = A
  !!     c_daadd(A,B,C):       PERFORMS C = A + B
  !!    c_c_dasub(A,B,C):       PERFORMS C = A - B
  !!     c_damul(A,B,C):       PERFORMS C = A * B
  !!     c_dadiv(A,B,C):       PERFORMS C = A / B
  !!     DASQR(A,C):         PERFORMS C = A^2           (SQUARE OF A)
  !
  !!     c_dacad(A,RA,C):      PERFORMS C = A + RA
  !!     c_dacsu(A,RA,C):      PERFORMS C = A - RA
  !!     c_dasuc(A,RA,C):      PERFORMS C = RA - A
  !!     c_dacmu(A,RA,C):      PERFORMS C = A * RA
  !!    c_dacdi(A,RA,C):      PERFORMS C = A / RA
  !!     c_dadic(A,RA,C):      PERFORMS C = RA / A
  !!     DACMA(A,B,RB,C):    PERFORMS C = A + RB*B
  !!c_damulIN(A,B,RA,C,D,RB,C):    PERFORMS C = A*B*RA + C*D*RB
  !!     DALIN(A,RA,B,RB,C): PERFORMS C = A*RA + B*RB
  !!     c_dafun(F,A,C):       PERFORMS C = F(A)          (DA FUNCTION)
  !
  !!     c_daabs((A,R):         PERFORMS R = |A|           (NORM OF A)
  !!     DACOM(A,B,R):       PERFORMS R = |A-B|         (NORM OF A-B)
  !!     DAPOS(A,C):         PERFORMS C(I) = |A(I)|     (MAKE SIGNS POSITIVE)
  !
  !!     DACCT(X,I,Y,J,Z,K)  CONCATENATES Z = X o Y;   I,J,K: # OF VECTORS IN X,Y,
  !!     DAINV(X,I,Z,K)      INVERTS Z = X^-1;           I,J: # OF VECTORS IN X,Y
  !!     DAPIN(X,I,Z,K,JJ)   PARTIALLY INVERTS Z = X^-1; I,J: # OF VECTORS IN X,Y,
  !                         JJ: ARRAY; NONZERO ENTRIES DENOTE TO BE INVERTED LINES
  !
  !!     DADER(I,A,C):       PERFORMS C = DA/DI (DERIV. WITH RESPECT TO VARIABLE I
  !!     DAPOI(A,B,C,I):     PERFORMS C = [A,B] (POISSON BRACKET, 2*I: # PHASEVARS
  !!     c_dacfu(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
  !
  !     c_dapri(A,U):         PRINTS DA VECTOR A TO UNIT U
  !     DAREA(A,U):         READS DA VECTOR A FROM UNIT U
  !     DADEB(U,T,I):       DEBUGGER, DUMPS TO U. T: MEMO, I=0: RETURN, I=1:STOP
  !!     DARAN(A,R,seed):         FILLS A WITH RANDOM NUMBERS. R: FILLFACTOR
  !     DANUM(O,N,I):       COMPUTES NUMBER OF MONOMIALS IN N VAR THROUGH ORDER O
  !
  !
  !     ADDITIONAL ROUTINES THE USER DOES NOT NEED TO CALL:
  !
  !     DAINF: RETURNS INFOS ABOUT A DA VECTOR PREVIOUSLY DECLARED
  !     DAPAC: PACKS DA VECTORS
  !     DACHK: CHECKS IF DA VECTORS HAVE COMPATIBLE ATTRIBUTES
  !     DCODE: TRANSFORMS DIGITS IN A CERTAIN BASE TO A DECIMAL INTEGER
  !     NCODE: EXTRACTS DIGITS IN A CERTAIN BASE FROM A DECIMAL INTEGER
  !
  !
  !     FURTHER WISHES
  !     **************
  !
  !     - CHECK DAREA AND c_dapri FOR CA VECTORS
  !     - MAKE DAFUN USE DASQR
  !
  !
  !      BLOCKDATA DABLD
  !     ***************
  !
  !
  !     PARAMETERS:
  !
  !     c_lda: MAXIMUM NUMBER OF DA-VECTORS;    CAN BE CHANGED QUITE ARBITRARILY
  !     c_lst: LENGTH OF MAIN STORAGE STACK;    CAN BE CHANGED QUITE ARBITRARILY
  !     c_lea: MAXIMUM NUMBER OF MONOMIALS;     CAN BE INCREASED FOR LARGE NO,NV
  !     c_lia: DIMENSION OF c_ia1,c_ia2;            CAN BE INCREASED FOR LARGE NO,NV
  !     c_lno: MAXIMUM ORDER;                   CAN BE INCREASED TO ABOUT 1000
  !     c_lnv: MAXIMUM NUMBER OF VARIABLES;     CAN BE INCREASED TO ABOUT 1000
  !
  !     ALL THE CHANGES IN THE VALUES OF PARAMETERS HAVE TO BE MADE BY GLOBAL
  !     SUBSTITUTIONS IN ALL SUBROUTINES.
  !
  !     c_daname:   NAME OF DA VECTOR
  !
  !     CC:       STACK OF DOUBLE PRECISON COEFFICIENTS
  !     c_i_1:       FIRST CHARACTERISTIC INTEGER (CF c_daini)
  !     c_i_2:       SECOND CHARACTERISTIC INTEGER (CF c_daini)
  !
  !     c_ie1:      CHARACTERISTIC INTEGER 1 OF UNPACKED REPRESENTATION (CF c_daini)
  !     c_ie2:      CHARACTERISTIC INTEGER 2 OF UNPACKED REPRESENTATION (CF c_daini)
  !     c_ieo:      ORDER OF ENTRY IN UNPACKED REPRESENTATION
  !     c_ia1:      REVERSE TO c_ie1 (CF c_daini)
  !     c_ia2:      REVERSE TO c_ie2 (CF c_daini)
  !
  !     c_idano:    ORDER OF DA VECTOR; IN CA, NUMBER OF COMPONENTS
  !     c_idanv:    NUMBER OF VARIABLES; IF 0, INDICATES CA VECTOR
  !     c_idapo:    FIRST ADDRESS IN STACK
  !     c_idalm:    NUMBER OF RESERVED STACK POSITIONS
  !     c_idall:    NUMBER OF MOMENTARILY REQUIRED STACK POSITIONS
  !
  !     c_nda_dab:      NUMBER OF DA VECTORS MOMENTARILY DEFINED
  !     NST:      NUMBER OF STACK POSITIONS MOMENTARILY ALLOCATED
  !     c_nomax:    MAXIMUM REQUESTED ORDER  (CF c_daini)
  !     c_nvmax:    MAXIMUM REQUESTED NUMBER OF VARIABLES (CF c_daini)
  !     c_nmmax:    MAXIMUM NUMBER OF MONOMIALS FOR c_nomax, c_nvmax (CF c_daini)
  !     c_nocut:    MOMENTARY TRUNCATION ORDER
  !     EPS:      TRUNCATION ACCURACY (CAN BE SET BY USER)
  !     EPSMAC:   MANTISSA LENGTH OF MACHINE (PESSIMISTIC ESTIMATE)
  !
  !-----------------------------------------------------------------------------
  !
 
  subroutine c_daini(no,nv,iunit)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE SETS UP THE MAJOR ORDERING AND ADDRESSING ARRAYS IN
    !     COMMON BLOCK c_daini. IF IUNIT > 0, THE ARRAYS WILL BE PRINTED TO UNIT
    !     NUMBER IUNIT. AN EXAMPLE FOR THE ARRAYS GENERATED BY c_daini CAN BE
    !     FOUND AFTER THE ROUTINE.
    !
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    !      COMMON / DASCR /  IS(20), RS(20)
    !integer idao,is,iscrri
    !real(dp) rs
    !common/dascr/is(100),rs(100),iscrri(100),idao
    !-----------------------------------------------------------------------------
    !
    integer i,iall,ibase,ic1,ic2,icmax,io1,io2,iout,iunit,jd,jjj,jjjj,jl,js,&
         nn,no,nv,ipause,mypauses
    integer,dimension(c_lnv+1)::n
    integer,dimension(0:c_lnv)::k
    integer,dimension(c_lnv)::j,jj
    character(10) aa
    !
    !frs if(eps.le.zero) eps=c_1d_38
    !      if(EPS.le.zero) eps=c_1d_90
    !frs epsmac=c_1d_7
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    c_last_tpsa=1
    if(nv.eq.0) return
    call alloc_all_c(no,nv)
    c_ndamaxi=0
    !
    do i=1, c_lda
       c_allvec(i) = .false.
    enddo
    c_nhole=0
    !*****************************************
    !     INITIALIZING VARIABLES IN COMMON / DA /
    !     ***************************************
    !
    c_nda_dab   = 0
    c_nst0   = 0
    c_nomax = no
    c_nvmax = nv
    call danum_c(no,nv,c_nmmax)
    c_nocut = no
    c_lfi   = 0
    !
    do i=0,c_lia
       c_ia1(i) = 0
       c_ia2(i) = 0
    enddo
    !
    !    do i=1,100
    !       is(i) = 0
    !    enddo
    !
    if(nv.gt.c_lnv.or.no.gt.c_lno) then
       write(6,*) 'ERROR IN SUBROUTINE c_daini, NO, NV = ',no,nv
       ipause=mypauses(1,line)
       call dadeb !(31,'ERR c_daini ',1)
    endif
    !
    ibase = no+1
    js    = nv/2
    if(float(ibase)**((nv+1)/2).gt.float(c_lia)) then
       write(line,'(a12,i4,a7,i4,a21,i4)') 'ERROR, NO = ',no,', NV = ',nv,' TOO LARGE FOR c_lia = ',c_lia
       ipause=mypauses(2,line)
       call dadeb !(31,'ERR c_daini ',1)
    endif
    !
    icmax = 0
    nn    = 0
    k(0)  = 0
    !
    do io2=0,no
       !     ***************
       !
       n(1)  = io2
       jl    = 0
       jd    = 1
       !
50     jl    = jl + jd
       !
       !old
       !      IF(JL.EQ.0) THEN
       !old
       !     modified according to Wu Ying
       !
       if(jl.le.0) then
          !
          goto 100
       elseif(jd.eq.1) then
          j(jl) = 0
       else
          j(jl) = j(jl) + 1
       endif
       !
       k(jl)    = k(jl-1)*ibase + j(jl)
       n(jl+1)  = n(jl) - j(jl)
       !
       if(j(jl).gt.n(jl)) then
          jd    = -1
          goto 50
       elseif(jl.lt.js) then
          jd    = 1
          goto 50
       else
          j(jl) = n(jl)
          k(jl) = k(jl-1)*ibase + j(jl)
          ic2   = k(jl)
          icmax = max(icmax,ic2)
          k(jl) = 0
          !
          c_ia2(ic2) = nn
          !
          do io1=0,no-io2
             !        ******************
             !
             n(js+1) = io1
             jd      = 1
             !
70           jl      = jl + jd
             !
             if(jl.eq.js) then
                goto 80
             elseif(jd.eq.1) then
                j(jl) = 0
             else
                j(jl) = j(jl) + 1
             endif
             !
             k(jl)    = k(jl-1)*ibase + j(jl)
             n(jl+1)  = n(jl) - j(jl)
             !
             if(j(jl).gt.n(jl)) then
                jd    = -1
                goto 70
             elseif(jl.lt.nv) then
                jd    = 1
                goto 70
             else
                jd    = -1
                j(jl) = n(jl)
                k(jl) = k(jl-1)*ibase + j(jl)
                ic1   = k(jl)
                icmax = max(icmax,ic1)
                nn = nn + 1
                !
                if(c_etiennefix) then
                   if(nn<=c_lea) then   ! Etienne
                      c_ie1(nn) = ic1
                      c_ie2(nn) = ic2
                   endif             ! Etienne
                   if(nn<=c_lst) then   ! Etienne
                      c_i_1 (nn) = ic1
                      c_i_2 (nn) = ic2
                   endif
                   if(ic2.eq.0) c_ia1(ic1) = nn
                   if(nn<=c_lea) then   ! Etienne
                      c_ieo(nn) = io1 + io2
                   endif             ! Etienne
                   !
                else
                   c_ie1(nn) = ic1
                   c_ie2(nn) = ic2
                   c_i_1 (nn) = ic1
                   c_i_2 (nn) = ic2
                   if(ic2.eq.0) c_ia1(ic1) = nn
                   c_ieo(nn) = io1 + io2
                   !
                endif
                goto 70
             endif
             !
80           continue
          enddo
          !
          jd = -1
          goto 50
       endif
       !
100    continue
    enddo
    !
    if(nn.gt.c_lea.and.(.not.c_etiennefix)) then
       write(line,'(a21,i4,a12)') 'ERROR IN c_daini, NN = ',nn,' EXCEEDS c_lea'
       ipause=mypauses(3,line)
       call dadeb !(31,'ERR c_daini ',1)
    endif
    !
    !     ALLOCATING SCRATCH VARIABLES
    !     ****************************
    !
    iall = 0
    call daall1(iall,'$$UNPACK$$',c_nomax,c_nvmax)
    !
    do i=0,c_nomax
       aa = '$$MUL   $$'
       write(aa(6:10),'(I5)') i
       iall = 0
       !      CALL DAALL(IALL,1,AA,I,c_nvmax)
       call daall1(iall,aa,c_nomax,c_nvmax)
    enddo
    !
    c_idall(1) = c_nmmax
    !
    !     DOUBLE CHECKING ARRAYS c_ie1,c_ie2,c_ia1,c_ia2
    !     **************************************
    !
    do i=1,c_nmmax
       !
       jjj = c_ia1(c_ie1(i)) + c_ia2(c_ie2(i))
       if(jjj.ne.i) then
          write(line,'(a48,i4)') 'ERROR IN c_daini IN ARRAYS c_ie1,c_ie2,c_ia1,c_ia2 AT I = ',i
          ipause=mypauses(4,line)
          call dadeb !(31,'ERR c_daini ',1)
       endif
       !
    enddo
    !
    if(iunit.eq.0) return
    !
    write(line,'(a32)') 'ARRAY SETUP DONE, BEGIN PRINTING'
    ipause=mypauses(5,line)
    !
    iout = 32
    open(iout,file='c_daini.DAT',status='NEW')
    !CRAY OPEN(IOUT,FILE='c_daini',STATUS='UNKNOWN',FORM='FORMATTED')          *CRAY
    !CRAY REWIND IOUT                                                        *CRAY
    !
    write(iout,'(/A/A/)') ' ARRAYS c_i_1 THROUGH c_i_20, c_ie1,c_ie2,c_ieo **********************************'
    do i=1,c_nmmax
       call dancd(c_ie1(i),c_ie2(i),jj)
       write(iout,'(1X,I5,2X,4(5i2,1X),3I6)') i,(jj(jjjj),jjjj=1,c_lnv),c_ie1(i),c_ie2(i),c_ieo(i)
    enddo
    !
    write(iout,'(/A/A/)') ' ARRAYS c_ia1,c_ia2 **************'
    do i=0,icmax
       write(iout,'(3i10)') i,c_ia1(i),c_ia2(i)
    enddo
    !
    return
  end subroutine c_daini
  !  subroutine daexter
  !    implicit none
  !     *****************************
  !
  !-----------------------------------------------------------------------------
  !     integer c_i_1,c_i_2,c_ia1,c_ia2,c_idall,c_idalm,c_idano,c_idanv,c_idapo,c_ie1,c_ie2,c_ieo,ifi,c_lfi,nda,c_ndamaxi,c_nmmax,c_nocut,c_nomax,nst,c_nvmax
  !-----------------------------------------------------------------------------
  !
  !    integer i
  !
  !    if(.not.c_notallocated) then
  !       do i=1, c_lda
  !          c_allvec(i)=.false.
  !       enddo
  !    endif
  !    return
  !  end subroutine daexter
  subroutine dalc_lsta(c_ldanow)
    implicit none
    !     *****************************
    !
    !-----------------------------------------------------------------------------
    !
    integer i,c_ldanow
    !
    c_ldanow=0
    do i=1, c_lda
       if(c_allvec(i)) c_ldanow=c_ldanow+1
    enddo
    return
  end subroutine dalc_lsta
  !
  ! EXAMPLE: ARRAYS c_i_1 THROUGH c_i_20, c_ie1,c_ie2,c_ieo (c_nomax=3,c_nvmax=4)
  ! *************************************************************
  !     I   c_i_1               THROUGH               c_i_20     c_ie1   c_ie2   c_ieo
  !     1   0 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      0     0     0
  !     2   1 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      1     0     1
  !     3   0 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      4     0     1
  !     4   2 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      2     0     2
  !     5   1 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      5     0     2
  !     6   0 2 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      8     0     2
  !     7   3 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      3     0     3
  !     8   2 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      6     0     3
  !     9   1 2 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0      9     0     3
  !    10   0 3 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0     12     0     3
  !    11   0 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      0     1     1
  !    12   1 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      1     1     2
  !    13   0 1 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      4     1     2
  !    14   2 0 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      2     1     3
  !    15   1 1 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      5     1     3
  !    16   0 2 0 0 0  0 0 0 0 0  1 0 0 0 0  0 0 0 0 0      8     1     3
  !    17   0 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      0     4     1
  !    18   1 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      1     4     2
  !    19   0 1 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      4     4     2
  !    20   2 0 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      2     4     3
  !    21   1 1 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      5     4     3
  !    22   0 2 0 0 0  0 0 0 0 0  0 1 0 0 0  0 0 0 0 0      8     4     3
  !    23   0 0 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      0     2     2
  !    24   1 0 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      1     2     3
  !    25   0 1 0 0 0  0 0 0 0 0  2 0 0 0 0  0 0 0 0 0      4     2     3
  !    26   0 0 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      0     5     2
  !    27   1 0 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      1     5     3
  !    28   0 1 0 0 0  0 0 0 0 0  1 1 0 0 0  0 0 0 0 0      4     5     3
  !    29   0 0 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      0     8     2
  !    30   1 0 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      1     8     3
  !    31   0 1 0 0 0  0 0 0 0 0  0 2 0 0 0  0 0 0 0 0      4     8     3
  !    32   0 0 0 0 0  0 0 0 0 0  3 0 0 0 0  0 0 0 0 0      0     3     3
  !    33   0 0 0 0 0  0 0 0 0 0  2 1 0 0 0  0 0 0 0 0      0     6     3
  !    34   0 0 0 0 0  0 0 0 0 0  1 2 0 0 0  0 0 0 0 0      0     9     3
  !    35   0 0 0 0 0  0 0 0 0 0  0 3 0 0 0  0 0 0 0 0      0    12     3
  !
  !    ARRAYS c_ia1,c_ia2
  !    **************
  !    I        c_ia1       c_ia2
  !    0         1         0   c_ie1,c_ie2 AND c_ia1,c_ia2 ALLOW THE EASY COMPUTATION
  !    1         2        10   OF THE ADDRESS OF THE PRODUCT OF TWO MONOMIALS.
  !    2         4        22   LET IX AND IY BE THE POSITIONS OF THE TWO
  !    3         7        31   FACTORS. THEN THE POSITION IZ OF THE PRODUCT OF
  !    4         3        16   THE TWO FACTORS IS GIVEN BY
  !    5         5        25
  !    6         8        32   IZ = c_ia1(c_ie1(IX)+c_ie1(IY)) + c_ia2(c_ie2(IX)+c_ie2(IY))
  !    7         0         0
  !    8         6        28
  !    9         9        33   THE OTHER VARIABLES SET BY c_daini WOULD HAVE THE
  !   10         0         0   VALUES
  !   11         0         0
  !   12        10        34   c_nomax = 3,  c_nvmax = 4, c_nmmax = 35
  !
  subroutine daalc_lno1(ic,ccc)
    implicit none
    !     ********************************
    !
    !     THIS SUBROUTINE ALLOCATES STORAGE FOR A DA VECTOR WITH
    !     ORDER c_nomax AND NUMBER OF VARIABLES c_nvmax
    !
    !-----------------------------------------------------------------------------
    !
    logical(lp) incnda
    integer ind,ndanum,no,nv,ic,ipause,mypauses
    character(10) c,ccc
    !    if((.not.C_STABLE_DA)) then
    !       if(C_watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    no=c_nomax
    nv=c_nvmax
    ind = 1
    if(ic.gt.0.and.ic.le.c_nda_dab) then
       !         c_daname(IC(I)) = C
       !         IF(c_idano(IC(I)).EQ.NO.AND.c_idanv(IC(I)).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.c_nomax.or.nv.gt.c_nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' c_nomax, c_nvmax = ',c_nomax,c_nvmax
          ipause=mypauses(7,line)
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(c_nhole.gt.0) then
          ind=c_nda_dab
20        if (c_allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          c_nhole=c_nhole-1
       else
          incnda = .true.
          c_nda_dab = c_nda_dab + 1
          ind=c_nda_dab
          if(c_nda_dab.gt.c_lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(8,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
       endif
       if(ind>c_lda_max_used) c_lda_max_used=ind
       if(ind>c_lda) then
          write(6,*) "ind>c_lda ",c_lda,ind
          print*, 'ERROR IN DAALc_lno1, MAX NUMBER OF DA VECTORS EXHAUSTED: c_lda = ',c_lda
          stop
       endif
       c_allvec(ind) = .true.
       ic = ind
       !
       if(nv.ne.0) then
          call danum_c(no,nv,ndanum)
       else
          ndanum = no
       endif
       c = ccc
       write(c(6:10),'(I5)') 1
       c_daname(ind) = c
       if (incnda) then
          if(ind.gt.c_nomax+2) then
             c_idano(ind) = c_nomax
             c_idanv(ind) = c_nvmax
             c_idapo(ind) = c_nst0 + 1
             c_idalm(ind) = c_nmmax
             c_idall(ind) = 0
             c_nst0 = c_nst0 + c_nmmax
          else
             c_idano(ind) = no
             c_idanv(ind) = nv
             c_idapo(ind) = c_nst0 + 1
             c_idalm(ind) = ndanum
             c_idall(ind) = 0
             c_nst0 = c_nst0 + ndanum
          endif
       endif
       !
       if(c_nst0.gt.c_lst) then
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(nv.eq.0.or.c_nomax.eq.1) then
          call daclr(ic)
          c_idall(ic) = c_idalm(ic)
       endif
    endif
    !
    if(c_nda_dab.gt.c_ndamaxi) c_ndamaxi=c_nda_dab
    return
  end subroutine daalc_lno1

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
    !    if((.not.C_STABLE_DA)) then
    !       if(C_watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    ind = 1
    do i=1,l
       if(ic(i).gt.0.and.ic(i).le.c_nda_dab) then
          !         c_daname(IC(I)) = C
          !         IF(c_idano(IC(I)).EQ.NO.AND.c_idanv(IC(I)).EQ.NV) THEN
       else
          if(nv.ne.0.and.(no.gt.c_nomax.or.nv.gt.c_nvmax)) then
             write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
                  &' c_nomax, c_nvmax = ',c_nomax,c_nvmax
             ipause=mypauses(9,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
          !
          if(c_nhole.gt.0) then
             ind=c_nda_dab
20           if (c_allvec(ind)) then
                ind = ind - 1
                goto 20
             endif
             incnda = .false.
             c_nhole=c_nhole-1
          else
             incnda = .true.
             c_nda_dab = c_nda_dab + 1
             ind=c_nda_dab
             if(c_nda_dab.gt.c_lda) then
                write(6,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
                !    ipause=mypauses(10,line)
                call dadeb !(31,'ERR DAALL ',1)
                stop 111
             endif
          endif
          !write(30,*) no,ind,c_lda,size(c_allvec)
          if(ind>c_lda_max_used) c_lda_max_used=ind
          if(ind>c_lda) then
             write(6,*) "ind>c_lda ",c_lda,ind
             print*, 'ERROR IN DAALc_lno1, MAX NUMBER OF DA VECTORS EXHAUSTED: c_lda = ',c_lda
             stop
          endif
          c_allvec(ind) = .true.
          ic(i) = ind
          !
          if(nv.ne.0) then
             call danum_c(no,nv,ndanum)
          else
             ndanum = no
          endif
          c = ccc
          if(l.ne.1) write(c(6:10),'(I5)') i
          c_daname(ind) = c
          if (incnda) then
             if(ind.gt.c_nomax+2) then
                c_idano(ind) = c_nomax
                c_idanv(ind) = c_nvmax
                c_idapo(ind) = c_nst0 + 1
                c_idalm(ind) = c_nmmax
                c_idall(ind) = 0
                c_nst0 = c_nst0 + c_nmmax
             else
                c_idano(ind) = no
                c_idanv(ind) = nv
                c_idapo(ind) = c_nst0 + 1
                c_idalm(ind) = ndanum
                c_idall(ind) = 0
                c_nst0 = c_nst0 + ndanum
             endif
          endif
          !
          if(c_nst0.gt.c_lst) then

             call dadeb !(31,'ERR DAALL ',1)
          endif
          !
          !          IF(NV.EQ.0) THEN
          if(nv.eq.0.or.c_nomax.eq.1) then
             call daclr(ic(i))
             c_idall(ic(i)) = c_idalm(ic(i))
          endif
       endif
    enddo
    !
    if(c_nda_dab.gt.c_ndamaxi) c_ndamaxi=c_nda_dab
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
    !    if((.not.C_STABLE_DA)) then
    !       if(C_watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    ind = 1
    if(ic.gt.0.and.ic.le.c_nda_dab) then
       !         c_daname(ic) = C
       !         IF(c_idano(ic).EQ.NO.AND.c_idanv(ic).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.c_nomax.or.nv.gt.c_nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' c_nomax, c_nvmax = ',c_nomax,c_nvmax
          ipause=mypauses(11,line)
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(c_nhole.gt.0) then
          ind=c_nda_dab
20        if (c_allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          c_nhole=c_nhole-1
       else
          incnda = .true.
          c_nda_dab = c_nda_dab + 1
          ind=c_nda_dab
          if(c_nda_dab.gt.c_lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(12,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
       endif
       if(ind>c_lda_max_used) c_lda_max_used=ind
       if(ind>c_lda) then
          write(6,*) "ind>c_lda ",c_lda,ind
          print*, 'ERROR IN DAALc_lno1, MAX NUMBER OF DA VECTORS EXHAUSTED: c_lda = ',c_lda
          stop
       endif
       c_allvec(ind) = .true.
       ic = ind
       !
       if(nv.ne.0) then
          call danum_c(no,nv,ndanum)
       else
          ndanum = no
       endif
       c = ccc
       write(c(6:10),'(I5)') 1
       c_daname(ind) = c
       if (incnda) then
          if(ind.gt.c_nomax+2) then
             c_idano(ind) = c_nomax
             c_idanv(ind) = c_nvmax
             c_idapo(ind) = c_nst0 + 1
             c_idalm(ind) = c_nmmax
             c_idall(ind) = 0
             c_nst0 = c_nst0 + c_nmmax
          else
             c_idano(ind) = no
             c_idanv(ind) = nv
             c_idapo(ind) = c_nst0 + 1
             c_idalm(ind) = ndanum
             c_idall(ind) = 0
             c_nst0 = c_nst0 + ndanum
          endif
       endif
       !
       if(c_nst0.gt.c_lst) then

          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       !          IF(NV.EQ.0) THEN
       if(nv.eq.0.or.c_nomax.eq.1) then
          call daclr(ic)
          c_idall(ic) = c_idalm(ic)
       endif
    endif
    !
    if(c_nda_dab.gt.c_ndamaxi) c_ndamaxi=c_nda_dab
    return
  end subroutine daall1

  subroutine c_etall1(ic)
    implicit none
    integer ic
    ic=0
    call c_daall0(ic)
  end subroutine c_etall1

  subroutine c_daall0(ic)
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
    ccc='         '
    no=c_nomax
    nv=c_nvmax
    !    if((.not.C_STABLE_DA)) then
    !       if(C_watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    ind = 1
    if(ic.gt.0.and.ic.le.c_nda_dab) then
       !         c_daname(ic) = C
       !         IF(c_idano(ic).EQ.NO.AND.c_idanv(ic).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.c_nomax.or.nv.gt.c_nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' c_nomax, c_nvmax = ',c_nomax,c_nvmax
          ipause=mypauses(11,line)
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(c_nhole.gt.0) then
          ind=c_nda_dab
20        if (c_allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          c_nhole=c_nhole-1
       else
          incnda = .true.
          c_nda_dab = c_nda_dab + 1
          ind=c_nda_dab
          if(c_nda_dab.gt.c_lda) then
             write(6,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             write(6,*) c_nda_dab,c_lda
             call dadeb !(31,'ERR DAALL ',1)
          endif
       endif
       if(ind>c_lda_max_used) c_lda_max_used=ind
       if(ind>c_lda) then
          write(6,*) "ind>c_lda ",c_lda,ind
          print*, 'ERROR IN DAALc_lno1, MAX NUMBER OF DA VECTORS EXHAUSTED: c_lda = ',c_lda
          stop
       endif
       c_allvec(ind) = .true.
       ic = ind
       !
       if(nv.ne.0) then
          call danum_c(no,nv,ndanum)
       else
          ndanum = no
       endif
       c = ccc
       write(c(6:10),'(I5)') 1
       c_daname(ind) = c
       if (incnda) then
          if(ind.gt.c_nomax+2) then
             c_idano(ind) = c_nomax
             c_idanv(ind) = c_nvmax
             c_idapo(ind) = c_nst0 + 1
             c_idalm(ind) = c_nmmax
             c_idall(ind) = 0
             c_nst0 = c_nst0 + c_nmmax
          else
             c_idano(ind) = no
             c_idanv(ind) = nv
             c_idapo(ind) = c_nst0 + 1
             c_idalm(ind) = ndanum
             c_idall(ind) = 0
             c_nst0 = c_nst0 + ndanum
          endif
       endif
       !
       if(c_nst0.gt.c_lst) then

          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       !          IF(NV.EQ.0) THEN
       if(nv.eq.0.or.c_nomax.eq.1) then
          call daclr(ic)
          c_idall(ic) = c_idalm(ic)
       endif
    endif
    !
    if(c_nda_dab.gt.c_ndamaxi) c_ndamaxi=c_nda_dab
    return
  end subroutine c_daall0
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
    !    if((.not.C_STABLE_DA)) then
    !       if(C_watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    do i=l,1,-1
       if(idal(i).le.c_nomax+2.or.idal(i).gt.c_nda_dab) then
          write(line,'(a38,i8,1x,i8)') 'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),c_nda_dab
          !ipause=mypauses(13,line)
          C_%STABLE_DA = .false.
          l = 1
          return
          call dadeb !(31,'ERR DADAL ',1)
       endif
       if(idal(i).eq.c_nda_dab) then
          !       deallocate
          c_nst0 = c_idapo(c_nda_dab) - 1
          c_nda_dab = c_nda_dab - 1
       else
          c_nhole=c_nhole+1
       endif
       c_allvec(idal(i)) = .false.
       !        c_idano(IDAL(I)) = 0
       !        c_idanv(IDAL(I)) = 0
       !        c_idapo(IDAL(I)) = 0
       !        c_idalm(IDAL(I)) = 0
       c_idall(idal(i)) = 0
       idal(i) = 0
    enddo
    return
  end subroutine dadal

  subroutine c_dadal1(idal)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
    !
    !-----------------------------------------------------------------------------
    !
    integer idal
    !
    !    if((.not.C_STABLE_DA)) then
    !       if(C_watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    if(idal.le.c_nomax+2.or.idal.gt.c_nda_dab) then
       write(6,'(a35,i8,1x,i8)') 'ERROR IN ROUTINE DADAL, IDAL,NDA = ',idal,c_nda_dab
       call dadeb !(31,'ERR DADAL ',1)
    endif
    if(idal.eq.c_nda_dab) then
       !       deallocate
       c_nst0 = c_idapo(c_nda_dab) - 1
       c_nda_dab = c_nda_dab - 1
    else
       c_nhole=c_nhole+1
    endif
    c_allvec(idal) = .false.
    !        c_idano(IDAL(I)) = 0
    !        c_idanv(IDAL(I)) = 0
    !        c_idapo(IDAL(I)) = 0
    !        c_idalm(IDAL(I)) = 0
    c_idall(idal) = 0
    idal = 0
    return
  end subroutine c_dadal1

  subroutine c_count_da(n)
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
    do i=1,c_lda
       if(c_allvec(i)) n=n+1
    enddo
    return
  end subroutine c_count_da

  subroutine c_davar(ina,ckon,i)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE DECLARES THE DA VECTOR
    !     AS THE INDEPENDENT VARIABLE NUMBER I.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic1,ic2,illa,ilma,ina,inoa,inva,ipoa,ipause,mypauses
    complex(dp) ckon
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(i.gt.inva) then
       write(line,'(a20,i8,a16,i8)') 'ERROR IN c_davar, I = ',i,' EXCEEDS INVA = ',inva
       ipause=mypauses(14,line)
       call dadeb !(31,'ERR c_davar ',1)
    endif
    !
    if(c_nomax.eq.1) then
       if(i.gt.inva) then
          print*,'ERROR IN c_davar, I = ',i,' EXCEEDS INVA = ',inva
          !           call dadeb !(31,'ERR c_davar3',1)
       endif
       call daclr(ina)
       c_cc(ipoa) = ckon
       c_cc(ipoa+i) = 1.0_dp
       return
    endif
    ibase = c_nomax+1
    !
    if(i.gt.(c_nvmax+1)/2) then
       ic1 = 0
       ic2 = ibase**(i-(c_nvmax+1)/2-1)
    else
       ic1 = ibase**(i-1)
       ic2 = 0
    endif
    !
    if(abs(ckon).gt.eps) then
       c_idall(ina) = 2
       c_cc(ipoa) = ckon
       c_i_1(ipoa) = 0
       c_i_2(ipoa) = 0
       !
       c_cc(ipoa+1) = 1.0_dp
       c_i_1(ipoa+1) = ic1
       c_i_2(ipoa+1) = ic2
    else
       c_idall(ina) = 1
       c_cc(ipoa) = 1.0_dp
       c_i_1(ipoa) = ic1
       c_i_2(ipoa) = ic2
    endif
    !
    return
  end subroutine c_davar
  !
  subroutine c_dacon(ina,ckon)
    implicit none
    !     **************************
    !
    !     THIS SUBROUTINE SETS THE VECTOR C TO THE CONSTANT CKON
    !
    !-----------------------------------------------------------------------------
    !
    integer illa,ilma,ina,inoa,inva,ipoa
    complex(dp) ckon
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(c_nomax.eq.1) then
       call daclr(ina)
       c_cc(ipoa) = ckon
       return
    endif
    c_idall(ina) = 1
    c_cc(ipoa) = ckon
    c_i_1(ipoa) = 0
    c_i_2(ipoa) = 0
    if(abs(ckon).lt.eps) c_idall(ina) = 0
    !
    return
  end subroutine c_dacon
  !
  subroutine c_danot(not)
    implicit none
    !     *********************
    !
    !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER c_nocut TO A NEW VALUE
    !
    !-----------------------------------------------------------------------------
    !
    integer not,ipause,mypauses
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(not.gt.c_nomax) then
       write(line,'(a15,i8,a17,i8)') 'ERROR, c_nocut = ',c_nocut,' EXCEEDS c_nomax = ',c_nomax
       ipause=mypauses(15,line)
       call dadeb !(31,'ERR DANOT ',1)
    endif
    !
    c_nocut = not
    !
    return
  end subroutine c_danot
  !  subroutine getdanot(not)
  !    implicit none
  !     *********************
  !
  !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER c_nocut TO A NEW VALUE
  !
  !-----------------------------------------------------------------------------
  !
  !    integer not,ipause,mypauses
  !
  !    if(not.gt.c_nomax) then
  !       write(line,'(a15,i8,a17,i8)') 'ERROR, c_nocut = ',c_nocut,' EXCEEDS c_nomax = ',c_nomax
  !       ipause=mypauses(15,line)
  !       call dadeb !(31,'ERR DANOT ',1)
  !    endif
  !
  !    not=c_nocut
  !
  !    return
  !  end subroutine getdanot
  subroutine c_daeps(deps)
    implicit none
    !     **********************
    !
    !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER c_nocut TO A NEW VALUE
    !
    !-----------------------------------------------------------------------------
    !
    real(dp) deps
    !
    if(deps.ge.0.0_dp) then
       eps = deps
    else
       deps=eps
    endif
    !
    return
  end subroutine c_daeps

  subroutine c_print_eps(deps)
    implicit none
    !     **********************
    !
    !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER c_nocut TO A NEW VALUE
    !
    !-----------------------------------------------------------------------------
    !
    real(dp) deps
    !
    if(deps.ge.0.0_dp) then
       epsprint = deps
    else
       deps=epsprint
    endif
    !
    return
  end subroutine c_print_eps


  !
  subroutine c_dapek(ina,jv,cjj)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE ARRAY
    !     OF EXPONENTS JJ AND RETURNS IT IN CJJ
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ic1,ic2,icu,icz,ic_i_1,illa,ilma,ina,inoa,inva,ipek,ipoa,&
         iu,iz,jj1,mchk,ipause,mypauses
    integer,dimension(:)::jv     ! 2002.12.4
    integer,dimension(c_lnv)::jj
    complex(dp) cjj
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    jj=0
    do i=1,size(jv)
       jj(i)=jv(i)
    enddo
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(illa.eq.0) then   ! Etienne shit
       cjj = 0
       return
    endif
    jj1 = 1
    if(inva.eq.0.or.c_nomax.eq.1) then
       if(inva.ne.0.and.c_nomax.eq.1) then
          if(illa.ge.2) then
             do i=1,illa - 1
                if(jj(i).eq.1) jj1 = i + 1
             enddo
          else
             jj1 = jj(1) + 1
          endif
       else
          jj1 = jj(1)
       endif
       if(jj1.lt.1.or.jj1.gt.illa) then
          print*,'ERROR IN c_dapek, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
          !           call dadeb !(31,'ERR c_dapek1',1)
       endif
       ipek = ipoa + jj1 - 1
       cjj = c_cc(ipek)
       return
    endif
    ic_i_1 = (c_nvmax+1)/2
    ibase = c_nomax+1
    !
    !     DETERMINE INDEX TO BE SEARCHED FOR
    !     **********************************
    !
    call dadcd(jj,ic1,ic2)
    !
    !ETIENNE
    if(ic1.gt.c_lia.or.ic2.gt.c_lia) then
       write(line,'(a24,i8)') 'DISASTER IN c_dapek, INA= ',ina
       ipause=mypauses(16,line)
    endif
    !ETIENNE
    ic = c_ia1(ic1) + c_ia2(ic2)
    !
    !     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,c_nocut
    !     ****************************************************************
    !
    !      IF(ICO.GT.INOA.OR.ICV.GT.INVA.OR.ICO.GT.c_nocut) THEN
    !         CJJ = 0
    !         RETURN
    !      ENDIF
    !
    !     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
    !     *************************************************************
    !
    iu = ipoa
    iz = ipoa + illa - 1
    icu = c_ia1(c_i_1(iu))+c_ia2(c_i_2(iu))
    icz = c_ia1(c_i_1(iz))+c_ia2(c_i_2(iz))
    !
    if(illa.eq.0) then
       cjj = 0
       return
    elseif(ic.eq.icu) then
       cjj = c_cc(iu)
       return
    elseif(ic.eq.icz) then
       cjj = c_cc(iz)
       return
    elseif(ic.lt.icu.or.ic.gt.icz) then
       cjj = 0
       return
    endif
    !
    !     SEARCHING PROPER MONOMIAL
    !     *************************
    !
10  continue
    if(iz-iu.le.1) then
       cjj = 0
       return
    endif
    i = (iu+iz)/2
    !
    !     if(c_ia1(c_i_1(i))+c_ia2(c_i_2(i)) - ic) 20,30,40
    mchk=c_ia1(c_i_1(i))+c_ia2(c_i_2(i)) - ic
    if(mchk.lt.0) goto 20
    if(mchk.eq.0) goto 30
    if(mchk.gt.0) goto 40
20  iu = i
    goto 10
30  cjj = c_cc(i)
    return
40  iz = i
    goto 10
    !
  end subroutine c_dapek
  !
  subroutine c_dapok(ina,jv,cjj)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE SETS THE COEFFICIENT OF THE ARRAY
    !     OF EXPONENTS JJ TO THE VALUE CJJ
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ic,ic1,ic2,icu,icz,ii,illa,ilma,ina,inoa,inva,ipoa,ipok,&
         iu,iz,jj1,mchk,ipause,mypauses
    integer,dimension(:)::jv     ! 2002.12.4
    integer,dimension(c_lnv)::jj
    complex(dp) cjj
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    jj=0
    do i=1,size(jv)
       jj(i)=jv(i)
    enddo
    !
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    jj1 = 1
    if(inva.eq.0.or.c_nomax.eq.1) then
       if(inva.ne.0.and.c_nomax.eq.1) then
          if(illa.ge.2) then
             do i=1,illa - 1
                if(jj(i).eq.1) jj1 = i + 1
             enddo
          else
             jj1 = jj(1) + 1
          endif
       else
          jj1 = jj(1)
       endif
       if(jj1.lt.1.or.jj1.gt.illa) then
          print*,'ERROR IN c_dapok, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
          !           call dadeb !(31,'ERR c_dapok1',1)
       endif
       ipok = ipoa + jj1 - 1
       c_cc(ipok) = cjj
       return
    endif
    !     DETERMINE INDEX TO BE SEARCHED FOR
    !     **********************************
    !
    call dadcd(jj,ic1,ic2)
    !
    ic = c_ia1(ic1) + c_ia2(ic2)
    !
    !     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,c_nocut
    !     ****************************************************************
    !
    !
    if(illa.ne.0) then ! etienne shit
       iu = ipoa
       iz = ipoa + illa - 1
       !
       !     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
       !     *************************************************************
       !
       icu = c_ia1(c_i_1(iu))+c_ia2(c_i_2(iu))
       icz = c_ia1(c_i_1(iz))+c_ia2(c_i_2(iz))
    endif
    if(illa.eq.0) then
       i = ipoa
       goto 100
    elseif(ic.eq.icu) then
       c_cc(iu) = cjj
       i = iu
       goto 200
    elseif(ic.eq.icz) then
       c_cc(iz) = cjj
       i = iz
       goto 200
    elseif(ic.lt.icu) then
       i = iu
       goto 100
    elseif(ic.gt.icz) then
       i = iz + 1
       goto 100
    endif
    !
    !
    !     SEARCHING PLACE TO POKE INTO OR BEFORE WHICH TO POKE
    !     ****************************************************
    !
    iu = ipoa
    iz = ipoa + illa
    !
10  continue
    if(iz-iu.le.1) then
       i = iz
       goto 100
    endif
    i = (iu+iz)/2
    !
    !      if(c_ia1(c_i_1(i))+c_ia2(c_i_2(i)) - ic) 20,30,40
    mchk=c_ia1(c_i_1(i))+c_ia2(c_i_2(i)) - ic
    if(mchk.lt.0) goto 20
    if(mchk.eq.0) goto 30
    if(mchk.gt.0) goto 40
20  iu = i
    goto 10
30  c_cc(i) = cjj
    goto 200
40  iz = i
    goto 10
    !
    !     INSERTING THE MONOMIAL, MOVING THE REST
    !     ***************************************
    !
100 continue
    !
    if(abs(cjj).lt.eps) return
    !
    do ii=ipoa+illa,i+1,-1
       c_cc(ii) = c_cc(ii-1)
       c_i_1(ii) = c_i_1(ii-1)
       c_i_2(ii) = c_i_2(ii-1)
    enddo
    !
    c_cc(i) = cjj
    c_i_1(i) = ic1
    c_i_2(i) = ic2
    !
    c_idall(ina) = illa + 1
    if(c_idall(ina).gt.c_idalm(ina)) then
       write(line,'(a15)') 'ERROR IN c_dapok '
       ipause=mypauses(17,line)
       call dadeb !(31,'ERR c_dapok ',1)
    endif
    !
    return
    !
    !     CASE OF CJJ = 0 WHICH MEANS MOVING THE REST
    !     *********************************************
    !
200 continue
    if(abs(cjj).lt.eps) then
       do ii=i,ipoa+illa-2
          c_cc(ii) = c_cc(ii+1)
          c_i_1(ii) = c_i_1(ii+1)
          c_i_2(ii) = c_i_2(ii+1)
       enddo
       c_idall(ina) = illa - 1
    endif
    return
    !
  end subroutine c_dapok
  !
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
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=ipoc,ipoc+ilmc-1
       !
       c_cc(i) = 0.0_dp
       !
    enddo
    !
    return
  end subroutine daclr
  !
  subroutine c_dacop(ina,inb)
    implicit none
    !     *************************
    !
    !     THIS SUBROUTINE COPIES THE DA VECTOR A TO THE DA VECTOR B
    !
    !-----------------------------------------------------------------------------
    !
    !      call dainf(ina,inoa,inva,ipoa,ilma,illa)
    !      call dainf(inb,inob,invb,ipob,ilmb,illb)
    !
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
    !-----------------------------------------------------------------------------
    !
    integer ia,ib,illa,ina,inb,ipoa,ipob
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ipob = c_idapo(inb)
    ipoa = c_idapo(ina)
    illa = c_idall(ina)
    ib = ipob - 1
    !
    !      iif = 0
    !      if(c_nomax.eq.1.or.inva.eq.0) iif = 1
    do ia = ipoa,ipoa+illa-1
       !
       if(c_nomax.gt.1) then
          if(c_ieo(c_ia1(c_i_1(ia))+c_ia2(c_i_2(ia))).gt.c_nocut) goto 100
       endif
       ib = ib + 1
       c_cc(ib) = c_cc(ia)
       c_i_1(ib) = c_i_1(ia)
       c_i_2(ib) = c_i_2(ia)
       !
100    continue
    enddo
    !
    c_idall(inb) = ib - ipob + 1
    return
  end subroutine c_dacop
  subroutine c_daadd(ina,inb,inc)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE PERFORMS A DA ADDITION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ina,ipoa
    integer idaadd,inb,inc,ipoc
    integer ipob
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,c_nvmax
          c_cc(ipoc+i) = c_cc(ipoa+i)   + c_cc(ipob+i)
       enddo
       return
    endif
    if(ina.ne.inc.and.inb.ne.inc) then
       call dalin(ina,(1.0_dp,0.0_dp),inb,(1.0_dp,0.0_dp),inc)
    else
       idaadd = 0
       call daall1(idaadd,'$$DAADD $$',c_nomax,c_nvmax)
       call dalin(ina,(1.0_dp,0.0_dp),inb,(1.0_dp,0.0_dp),idaadd)
       call c_dacop(idaadd,inc)
       call c_dadal1(idaadd)
    endif
    !
    return
  end subroutine c_daadd
  !
  subroutine c_datrunc(ina,io,inb)
    implicit none
    integer ina,io,inb,nt,ipoca,ipocb,i
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    nt=c_nocut
    if(io>c_nomax) then
       if(ina/=inb) call c_dacop(ina,inb)
       return
    endif
    c_nocut=io-1
      
     if(c_nomax==1.and.io<=1) then
       ipoca=c_idapo(ina)
       ipocb=c_idapo(inb)
       do i=1,c_nvmax
          c_cc(ipocb+i) =0.0_dp
       enddo
         c_cc(ipocb)=c_cc(ipoca)*io
    else
     call c_dacop(ina,inb)
    endif
    c_nocut = nt
  end subroutine c_datrunc

  subroutine c_dasub(ina,inb,inc)
    implicit none
    !     THIS SUBROUTINE PERFORMS A DA SUBTRACTION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C.
    !
    !-----------------------------------------------------------------------------
    !
    integer idasub
    integer i,ina,ipoa
    integer inc,ipoc,inb
    integer ipob
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,c_nvmax
          c_cc(ipoc+i) = c_cc(ipoa+i)   - c_cc(ipob+i)
       enddo
       return
    endif
    if(ina.ne.inc.and.inb.ne.inc) then
       call dalin(ina,(1.0_dp,0.0_dp),inb,-(1.0_dp,0.0_dp),inc)
    else
       idasub = -1
       !         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       call daall1(idasub,'$$DASUB $$',c_nomax,c_nvmax)
       call dalin(ina,(1.0_dp,0.0_dp),inb,-(1.0_dp,0.0_dp),idasub)
       call c_dacop(idasub,inc)
       call c_dadal1(idasub)
    endif
    !
    return
  end subroutine c_dasub

  subroutine c_damul(ina,inb,inc)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
    !     OF THE (c_nomax+2) SCRATCH VARIABLES ALLOCATED BY c_daini IS USED.
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inb,inc,incc,ipoc,ipoa,ipob,i
    complex(dp) ccipoa,ccipob
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoa=c_idapo(ina)
       ipob=c_idapo(inb)
       ipoc=c_idapo(inc)
       !         minv = min(inva,invb,invc)
       ccipoa = c_cc(ipoa)
       ccipob = c_cc(ipob)
       c_cc(ipoc) = ccipoa*ccipob
       do i=1,c_nvmax
          c_cc(ipoc+i) = ccipoa*c_cc(ipob+i) + ccipob*c_cc(ipoa+i)
       enddo
       !         do 30 i=ipoc+minv+1,ipoc+invc
       !  30     cc(i) = zero
       return
    endif
    if(ina.eq.inc.or.inb.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1(incc,'$$DAJUNK$$',c_nomax,c_nvmax)
       call damult(ina,inb,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call damult(ina,inb,inc)
    endif
    return
  end subroutine c_damul

  subroutine damult(ina,inb,inc)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
    !     OF THE (c_nomax+2) SCRATCH VARIABLES ALLOCATED BY c_daini IS USED.
    !
    !-----------------------------------------------------------------------------
    !
    !
    !
    !      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)
    !
    !     CASE OF FIRST ORDER ONLY
    !     ************************
    !
    !-----------------------------------------------------------------------------
    !
    integer i,c_i_1ia,c_i_2ia,ia,ib,ic,illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa,inob,inoc,&
         inva,invb,invc,ioffb,ipoa,ipob,ipoc,ipos,noib,nom
    integer,dimension(0:c_lno)::ipno,noff
    complex(dp) ccia,ccipoa,ccipob
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoa=c_idapo(ina)
       ipob=c_idapo(inb)
       ipoc=c_idapo(inc)
       !         minv = min(inva,invb,invc)
       ccipoa = c_cc(ipoa)
       ccipob = c_cc(ipob)
       c_cc(ipoc) = ccipoa*ccipob
       do i=1,c_nvmax
          c_cc(ipoc+i) = ccipoa*c_cc(ipob+i) + ccipob*c_cc(ipoa+i)
       enddo
       !         do 30 i=ipoc+minv+1,ipoc+invc
       !  30     c_cc(i) = zero
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inb,inob,invb,ipob,ilmb,illb)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !     GENERAL CASE
    !     ************
    !
    do i=0,c_nomax
       noff(i) = c_idapo(i+2)
       ipno(i) = 0
    enddo
    !
    call daclr(1)
    !
    !     RE-SORTING THE VECTOR B INTO PIECES THAT ARE OF ONLY ONE ORDER
    !     *************************************************************
    !
    do ib=ipob,ipob+illb-1
       !
       noib = c_ieo(c_ia1(c_i_1(ib))+c_ia2(c_i_2(ib)))
       ipos = ipno(noib) + 1
       ipno(noib) = ipos
       inob = noff(noib) + ipos
       !
       c_cc(inob) = c_cc(ib)
       c_i_1(inob) = c_i_1(ib)
       c_i_2(inob) = c_i_2(ib)
       !
    enddo
    !
    do i=0,c_nomax
       c_idall(i+2) = ipno(i)
    enddo
    !
    !     PERFORMING ACTUAL MULTIPLICATION
    !     ********************************
    !
    nom = min(c_nocut,inoc)
    !
    do ia=ipoa,ipoa+illa-1
       !
       c_i_1ia = c_i_1(ia)
       c_i_2ia = c_i_2(ia)
       ccia = c_cc(ia)
       !
       do noib = 0,nom-c_ieo(c_ia1(c_i_1(ia))+c_ia2(c_i_2(ia)))
          !
          ioffb = noff(noib)
          !
          do ib = ioffb+1,ioffb+ipno(noib)
             !
             ic = c_ia2(c_i_2ia+c_i_2(ib)) + c_ia1(c_i_1ia + c_i_1(ib))
             ! Georg says maybe needs if(ic/=0)
             if(ic/=0) then
                c_cc(ic) = c_cc(ic) + ccia*c_cc(ib)
             else
                write(6,*) " Georg warn me about ic could be 0.0_dp"
                stop 999
             endif
             !
          enddo
       enddo
    enddo
    !
    call dapac(inc)
    !
    return
  end subroutine damult
  !
  subroutine c_dadiv(ina,inb,inc)
    implicit none
    !     *************************
    !
    !     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
    !
    !-----------------------------------------------------------------------------
    !
    !     THIS SUBROUTINE PERFORMS A DA DIVISION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C.
    !
    !-----------------------------------------------------------------------------
    !
    integer idadiv,inb,ina,inc,ipoc,ipoa,ipob,i
    complex(dp) ck,ck1
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       ipoc = c_idapo(inc)
       ck=1.0_dp/c_cc(ipob)
       ck1=c_cc(ipoa)*ck
       do i=1,c_nvmax
          c_cc(ipoc+i) = (c_cc(ipoa+i)-c_cc(ipob+i)*ck1)*ck
       enddo
       c_cc(ipoc)=ck1
       return
    endif
    idadiv = 0
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    call daall1(idadiv,'$$DADIV $$',c_nomax,c_nvmax)
    call c_dafun('INV ',inb,idadiv)
    call c_damul(ina,idadiv,inc)
    call c_dadal1(idadiv)
    !
    return
  end subroutine c_dadiv
  !
  subroutine dasqr(ina,inc)
    implicit none
    !     *************************
    !
    !     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inc,incc,ipoc,i,ipoa
    complex(dp) ccipoa
    !
    !     CASE OF FIRST ORDER ONLY
    !     ************************
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       !         minv = min(inva,invc)
       ccipoa = c_cc(ipoa)
       c_cc(ipoc) = ccipoa*ccipoa
       do i=1,c_nvmax
          c_cc(ipoc+i) = 2.0_dp*ccipoa*c_cc(ipoa+i)
       enddo
       return
    endif
    if(ina.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1(incc,'$$DAJUNK$$',c_nomax,c_nvmax)
       call dasqrt(ina,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call dasqrt(ina,inc)
    endif
    return
  end subroutine dasqr
  subroutine dasqrt(ina,inc)
    implicit none
    !     *************************
    !
    !     THIS SUBROUTINE SQUARES THE VECTOR A AND STORES THE RESULT IN C.
    !
    !-----------------------------------------------------------------------------
    !
    !
    !      CALL DACHK(INA,INOA,INVA,'          ',-1,-1,INC,INOC,INVC)
    !
    !
    !     CASE OF FIRST ORDER ONLY
    !     ************************
    !
    !-----------------------------------------------------------------------------
    !
    integer i,c_i_1ia,c_i_2ia,ia,ib,ib1,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,&
         inva,invc,ioffa,ioffb,ipoa,ipoc,ipos,noia,noib,nom
    integer,dimension(0:c_lno)::ipno,noff
    complex(dp) ccia,ccipoa
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       !         minv = min(inva,invc)
       ccipoa = c_cc(ipoa)
       c_cc(ipoc) = ccipoa*ccipoa
       do i=1,c_nvmax
          c_cc(ipoc+i) = 2.0_dp*ccipoa*c_cc(ipoa+i)
       enddo
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !     GENERAL CASE
    !     ************
    !
    do i=0,c_nomax
       noff(i) = c_idapo(i+2)
       ipno(i) = 0
    enddo
    !
    call daclr(1)
    !
    !     RESORTING THE VECTOR A INTO PIECES THAT ARE OF ONLY ONE ORDER
    !     *************************************************************
    !
    do ia=ipoa,ipoa+illa-1
       !
       noia = c_ieo(c_ia1(c_i_1(ia))+c_ia2(c_i_2(ia)))
       ipos = ipno(noia) + 1
       ipno(noia) = ipos
       inoa = noff(noia) + ipos
       !
       c_cc(inoa) = c_cc(ia)
       c_i_1(inoa) = c_i_1(ia)
       c_i_2(inoa) = c_i_2(ia)
       !
    enddo
    !
    do i=0,c_nomax
       c_idall(i+2) = ipno(i)
    enddo
    !
    !     PERFORMING ACTUAL MULTIPLICATION
    !     ********************************
    !
    nom = min(c_nocut,inoc)
    !
    do noia = 0,nom/2
       !
       ioffa = noff(noia)
       !
       do ia=ioffa+1,ioffa+ipno(noia)
          !
          c_i_1ia = c_i_1(ia)
          c_i_2ia = c_i_2(ia)
          ccia = c_cc(ia)
          !
          ic = c_ia2(c_i_2ia+c_i_2ia) + c_ia1(c_i_1ia+c_i_1ia)
          c_cc(ic) = c_cc(ic) + ccia*ccia
          ccia = ccia + ccia
          !
          do noib = noia,nom-noia
             !
             ioffb = noff(noib)
             if(noib.eq.noia) then
                ib1 = ia + 1
             else
                ib1 = ioffb + 1
             endif
             !
             do ib = ib1,ioffb+ipno(noib)
                !
                ic = c_ia2(c_i_2ia+c_i_2(ib)) + c_ia1(c_i_1ia + c_i_1(ib))
                c_cc(ic) = c_cc(ic) + ccia*c_cc(ib)
                !
             enddo
          enddo
       enddo
    enddo
    !
    call dapac(inc)
    !
    return
  end subroutine dasqrt
  !
  subroutine c_dacad(ina,ckon,inb)
    !  use da_arrays
    implicit none
    !    integer,dimension(c_lnv)::jjy
    !    data jjy / c_lnv*0 /  ! flat zero here
    !     ******************************
    !
    !     THIS SUBROUTINE ADDS THE CONSTANT CKON TO THE VECTOR A
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inb
    integer,parameter,dimension(c_lnv)::jjx=0
    complex(dp) ckon,const
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call c_dacop(ina,inb)
    if(c_nomax.eq.1) then
       c_cc(c_idapo(inb)) = c_cc(c_idapo(inb)) + ckon
       return
    endif
    !
    call c_dapek(inb,jjx,const)
    call c_dapok(inb,jjx,const+ckon)
    !
    return
  end subroutine c_dacad
  !
  subroutine c_dacsu(ina,ckon,inb)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE SUBTRACTS THE CONSTANT CKON FROM THE VECTOR A
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inb
    integer,parameter,dimension(c_lnv)::jjx=0
    complex(dp) ckon,const
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call c_dacop(ina,inb)
    !
    if(c_nomax.eq.1) then
       c_cc(c_idapo(inb)) = c_cc(c_idapo(inb)) - ckon
       return
    endif
    !
    call c_dapek(inb,jjx,const)
    call c_dapok(inb,jjx,const-ckon)
    !
    return
  end subroutine c_dacsu
  !
  subroutine c_dasuc(ina,ckon,inb)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE SUBTRACTS THE VECTOR INA FROM THE CONSTANT CKON
    !
    !-----------------------------------------------------------------------------
    !
    !      call dainf(ina,inoa,inva,ipoa,ilma,illa)
    !      call dainf(inb,inob,invb,ipob,ilmb,illb)
    !-----------------------------------------------------------------------------
    !
    integer i,ina,inb,ipoa,ipob
    complex(dp) ckon
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ipob=c_idapo(inb)
    ipoa=c_idapo(ina)
    if(c_nomax.eq.1) then
       c_cc(ipob) = ckon - c_cc(ipoa)
       do i=1,c_nvmax
          c_cc(ipob+i) =-c_cc(ipoa+i)
       enddo
       return
    endif
    call c_dacsu(ina,ckon,inb)
    call c_dacmu(inb,-(1.0_dp,0.0_dp),inb)
    !
    return
  end subroutine c_dasuc
  !
  subroutine c_dacmu(ina,ckon,inc)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
    !     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
    !     THE DA VECTOR DENOTED WITH THE INTEGER E.
    !
    !-----------------------------------------------------------------------------
    !
    integer ipoa,i,ina,inc,incc,ipoc
    complex(dp) ckon
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = c_idapo(ina)
       ipoc = c_idapo(inc)
       do i=0,c_nvmax
          c_cc(ipoc+i) = c_cc(ipoa+i) * ckon
       enddo
       return
    endif
    if(ina.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daalc_lno1(incc,'$$DAJUNK$$')
       call dacmut(ina,ckon,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call dacmut(ina,ckon,inc)
    endif
    return
  end subroutine c_dacmu

  subroutine dacmut(ina,ckon,inb)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE MULTIPLIES THE DA VECTOR DENOTED BY THE
    !     THE INTEGER A WITH THE CONSTANT C AND STORES THE RESULT IN
    !     THE DA VECTOR DENOTED WITH THE INTEGER E.
    !
    !-----------------------------------------------------------------------------
    !
    !
    !
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INB,INOB,INVB)
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,illa,illb,ilma,ilmb,ina,inb,inoa,inob,inva,invb,ipoa,&
         ipob,ipause,mypauses
    complex(dp) ckon
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       do i=0,c_nvmax
          c_cc(ipob+i) = c_cc(ipoa+i) * ckon
       enddo
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inb,inob,invb,ipob,ilmb,illb)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(abs(ckon).lt.eps) then
       c_idall(inb) = 0
       return
    endif
    !
    ib = ipob - 1
    !
    do ia=ipoa,ipoa+illa-1
       !
       if(c_ieo(c_ia1(c_i_1(ia))+c_ia2(c_i_2(ia))).gt.c_nocut) goto 100
       ib = ib + 1
       c_cc(ib) = c_cc(ia)*ckon
       c_i_1(ib) = c_i_1(ia)
       c_i_2(ib) = c_i_2(ia)
       !
100    continue
    enddo
    !
    c_idall(inb) = ib-ipob+1
    if(c_idall(inb).gt.c_idalm(inb)) then
       write(line,'(a17)') 'ERROR IN DACMU '
       ipause=mypauses(18,line)
       call dadeb !(31,'ERR DACMU ',1)
    endif
    !
    return
  end subroutine dacmut
  !
  subroutine c_dacdi(ina,ckon,inb)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE DIVIDES THE VECTOR INA BY THE CONSTANT CKON
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ina,inb,ipoa,ipob,ipause,mypauses
    complex(dp) ckon
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(abs(ckon)==0.0_dp) then
       if(check_da) then
          C_STABLE_DA=.false.
          messagelost='constant part 0.0_dp in dacdi'
          return
       else
          write(line,'(a38)')  'ERROR IN DACDI  CKON IS 0.0_dp'
          ipause=mypauses(25,line)
       endif
    endif
    if(c_nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       do i=0,c_nvmax
          c_cc(ipob+i) = c_cc(ipoa+i)/ ckon
       enddo
       return
    endif
    !
    call c_dacmu(ina,1.0_dp/ckon,inb)
    !
    return
  end subroutine c_dacdi
  !
  !
  subroutine c_dadic(ina,ckon,inc)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE DIVIDES THE CONSTANT CKON BY THE VECTOR INA
    !
    !-----------------------------------------------------------------------------
    !
    integer i,idadic,ina,inc,ipoa,ipoc
    complex(dp) ckon,ck
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ipoa = c_idapo(ina)
    if(abs(c_cc(ipoa))==0.0_dp) then
       if(check_da) C_STABLE_DA=.false.
       messagelost='constant part 0.0_dp in dadic'
       return
    endif
    if(c_nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoc = c_idapo(inc)
       c_cc(ipoc)=ckon/c_cc(ipoa)
       ck=c_cc(ipoc)/c_cc(ipoa)
       do i=1,c_nvmax
          c_cc(ipoc+i) = -c_cc(ipoa+i)* ck
       enddo
       return
    endif
    !    if(abs(ckon).lt.eps) then    !2002.11.28
    !       call c_dacon(inc,zero)
    !       return
    !    endif
    idadic = 0
    call daall1(idadic,'$$DADIC $$',c_nomax,c_nvmax)
    !    if(ckon.eq.zero) then
    !       write(line,'(a18)') 'ERROR IN DACDI and DADIC, CKON IS ZERO' !2002.11.28
    !       ipause=mypauses(19,line)
    !       call dadeb !(31,'ERR DACDI ',1)
    !    endif
    !    call c_dacdi(ina,ckon,idadic)
    !    call c_dafun('INV ',idadic,inc)
    call c_dafun('INV ',ina,idadic)
    call c_dacmu(idadic,ckon,inc)
    call c_dadal1(idadic)
    !
    return
  end subroutine c_dadic
  !
  subroutine dacma(ina,inb,bfac,inc)
    implicit none
    !     **********************************
    !
    !     THIS SUBROUTINE PERFORMS THE OPERATIONS C = A + B*BFAC, WHERE A,B,C ARE
    !     DA VECTORS AND BFAC IS A complex(dp). A AND C CAN BE IDENTICAL.
    !     CAN LATER BE REPLACED BY SOMETHING LIKE DAADD WITH MINOR CHANGES.
    !
    !-----------------------------------------------------------------------------
    !
    integer idacma,ina,inb,inc,ipoc,ipob,ipoa,i
    complex(dp) bfac
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,c_nvmax
          c_cc(ipoc+i) = c_cc(ipoa+i)   + c_cc(ipob+i) * bfac
       enddo
       !         do 8 i=ipoc+minv+1,ipoc+invc
       ! 8       c_cc(i) = zero
       return
    endif
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    idacma = 0
    call daall1(idacma,'$$DACMA $$',c_nomax,c_nvmax)
    call dalin(ina,(1.0_dp,0.0_dp),inb,bfac,idacma)
    call c_dacop(idacma,inc)
    call c_dadal1(idacma)
    !
    return
  end subroutine dacma
  !
  subroutine dalin(ina,afac,inb,bfac,inc)
    implicit none
    integer ina,inb,inc,incc,ipoc
    complex(dp) afac,bfac
    !     ***************************************
    !
    !     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
    !     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
    !
    !-----------------------------------------------------------------------------
    !
    integer  ipob,ipoa,i
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,c_nvmax
          c_cc(ipoc+i) = c_cc(ipoa+i) * afac + c_cc(ipob+i) * bfac
       enddo
       !         do 8 i=ipoc+minv+1,ipoc+invc
       ! 8       c_cc(i) = zero
       return
    endif
    if(ina.eq.inc.or.inb.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1(incc,'$$DAJUNK$$',c_nomax,c_nvmax)
       call dalint(ina,afac,inb,bfac,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call dalint(ina,afac,inb,bfac,inc)
    endif
    return
  end subroutine dalin

  subroutine dalint(ina,afac,inb,bfac,inc)
    implicit none
    !     ***************************************
    !
    !     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
    !     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
    !
    !-----------------------------------------------------------------------------
    !
    !      CALL DACHK(INA,INOA,INVA, INB,INOB,INVB, INC,INOC,INVC)
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,iamax,ib,ibmax,ic,icmax,illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa,&
         inob,inoc,inva,invb,invc,ipoa,ipob,ipoc,is,ismax,ismin,ja,jb,mchk,&
         ipause,mypauses
    complex(dp) afac,bfac,ccc,copf
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(c_nomax.eq.1) then
       ipoc = c_idapo(inc)
       ipoa = c_idapo(ina)
       ipob = c_idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,c_nvmax
          c_cc(ipoc+i) = c_cc(ipoa+i) * afac + c_cc(ipob+i) * bfac
       enddo
       !         do 8 i=ipoc+minv+1,ipoc+invc
       ! 8       c_cc(i) = zero
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inb,inob,invb,ipob,ilmb,illb)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ia = ipoa
    ib = ipob
    ic = ipoc - 1
    iamax = ipoa+illa-1
    ibmax = ipob+illb-1
    icmax = ipoc+ilmc-1
    ja = c_ia1(c_i_1(ia)) + c_ia2(c_i_2(ia))
    jb = c_ia1(c_i_1(ib)) + c_ia2(c_i_2(ib))
    !
    if(ia.gt.iamax) then
       ismin = ib
       ismax = ibmax
       copf  = bfac
       goto 50
    endif
    if(ib.gt.ibmax) then
       ismin = ia
       ismax = iamax
       copf  = afac
       goto 50
    endif
    !
    !     COMPARING
    !     *********
    !
10  continue
    !      if(ja-jb) 30,20,40
    mchk=ja-jb
    if(mchk.lt.0) goto 30
    if(mchk.eq.0) goto 20
    if(mchk.gt.0) goto 40
    !
    !     ADDING TWO TERMS
    !     ****************
    !
20  continue
    ccc = c_cc(ia)*afac + c_cc(ib)*bfac
    if(abs(ccc).lt.eps) goto 25
    if(c_ieo(c_ia1(c_i_1(ia))+c_ia2(c_i_2(ia))).gt.c_nocut) goto 25
    ic = ic + 1
    c_cc(ic) = ccc
    c_i_1(ic) = c_i_1(ia)
    c_i_2(ic) = c_i_2(ia)
25  continue
    ia = ia + 1
    ib = ib + 1
    if(ia.gt.iamax) then
       ismin = ib
       ismax = ibmax
       copf  = bfac
       goto 50
    endif
    if(ib.gt.ibmax) then
       ismin = ia
       ismax = iamax
       copf  = afac
       goto 50
    endif
    ja = c_ia1(c_i_1(ia)) + c_ia2(c_i_2(ia))
    jb = c_ia1(c_i_1(ib)) + c_ia2(c_i_2(ib))
    goto 10
    !
    !     STORING TERM A
    !     **************
    !
30  continue
    if(c_ieo(c_ia1(c_i_1(ia))+c_ia2(c_i_2(ia))).gt.c_nocut) goto 35
    ccc = c_cc(ia)*afac
    if(abs(ccc).lt.eps) goto 35
    ic = ic + 1
    c_cc(ic) = ccc
    c_i_1(ic) = c_i_1(ia)
    c_i_2(ic) = c_i_2(ia)
35  continue
    ia = ia + 1
    if(ia.gt.iamax) then
       ismin = ib
       ismax = ibmax
       copf  = bfac
       goto 50
    endif
    ja = c_ia1(c_i_1(ia)) + c_ia2(c_i_2(ia))
    goto 10
    !
    !     STORING TERM B
    !     **************
    !
40  continue
    if(c_ieo(c_ia1(c_i_1(ib))+c_ia2(c_i_2(ib))).gt.c_nocut) goto 45
    ccc = c_cc(ib)*bfac
    if(abs(ccc).lt.eps) goto 45
    ic = ic + 1
    c_cc(ic) = ccc
    c_i_1(ic) = c_i_1(ib)
    c_i_2(ic) = c_i_2(ib)
45  continue
    ib = ib + 1
    if(ib.gt.ibmax) then
       ismin = ia
       ismax = iamax
       copf  = afac
       goto 50
    endif
    jb = c_ia1(c_i_1(ib)) + c_ia2(c_i_2(ib))
    goto 10
    !
    !     COPYING THE REST
    !     ****************
    !
50  continue
    do is=ismin,ismax
       if(c_ieo(c_ia1(c_i_1(is))+c_ia2(c_i_2(is))).gt.c_nocut) goto 60
       ccc = c_cc(is)*copf
       if(abs(ccc).lt.eps) goto 60
       ic = ic + 1
       c_cc(ic) = ccc
       c_i_1(ic) = c_i_1(is)
       c_i_2(ic) = c_i_2(is)
60     continue
    enddo
    !
    c_idall(inc) = ic - ipoc + 1
    !
    if(c_idall(inc).gt.c_idalm(inc)) then
       write(line,'(a40)')  'ERROR IN dalint c_idall(inc).gt.c_idalm(inc)'
       ipause=mypauses(21,line)
       call dadeb !(31,'ERR DALIN ',1)
    endif
    !
    return
  end subroutine dalint
  !

  !
  subroutine c_daabs(ina,anorm)
    implicit none
    !     ***************************
    !
    !     THIS SUBROUTINE COMPUTES THE NORM OF THE DA VECTOR A
    !
    !-----------------------------------------------------------------------------
    !
    integer i,illa,ilma,ina,inoa,inva,ipoa
    real(dp) anorm
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    anorm = 0.0_dp
    do i=ipoa,ipoa+illa-1
       anorm = anorm + abs(c_cc(i))
    enddo
    !
    return
  end subroutine c_DAABS
  !
  !
  subroutine c_dacct(ma,ia,mb,ib,mc,ic)
    implicit none
    !     ***********************************
    !
    !     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
    !     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
    !     DA VECTORS, RESPECTIVELY.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ic,ij,illc,ilmc,inoc,invc,ipoc
    integer,dimension(c_lnv)::monx
    integer,dimension(:)::ma,mb,mc
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(ma(1).eq.mc(1).or.mb(1).eq.mc(1)) then
       call dainf(mc(1),inoc,invc,ipoc,ilmc,illc)
       if((.not.C_STABLE_DA)) then
          if(C_watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       do ij=1,ic
          monx(ij)=0
       enddo
       call daall(monx,ic,'$$DAJUNK$$',inoc,invc)
       call dacctt(ma,ia,mb,ib,monx,ic)
       do i=1,ic
          call c_dacop(monx(i),mc(i))
       enddo
       call dadal(monx,ic)
    else
       call dacctt(ma,ia,mb,ib,mc,ic)
    endif
    return
  end subroutine c_dacct

  subroutine dacctt(mb,ib,mc,ic,ma,ia)
    implicit none
    !     ***********************************
    !
    !     THIS SUBROUTINE PERFORMS A CONCATENATION MA = MB o MC
    !     WHERE MA, MB AND MC ARE MATRICES CONSISTING OF IA, IB AND IC
    !     DA VECTORS, RESPECTIVELY.
    !
    !-----------------------------------------------------------------------------
    !
    !      INTEGER MON(c_lno+1),Ic_cc(c_lnv)
    !      INTEGER,dimension(:)::MB,MC,MA
    !ETIENNE
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ic,iia,iib,iic,illa,illb,illc,ilma,ilmb,ilmc,inoa,inob,inoc,inva,invb,invc,&
         ipoa,ipob,ipoc,iv,jl,jv,ipause,mypauses
    integer,dimension(c_lno+1)::mon
    !    integer,dimension(c_lno)::icc
    integer,dimension(c_lnv)::icc  !johan 2008 March
    integer,dimension(:)::ma,mb,mc
    complex(dp) ccf
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !ETIENNE
    !
    !     CONSISTENCY CHECKS
    !     ******************
    !
    iia = ma(1)
    iib = mb(1)
    iic = mc(1)
    call dainf(iia,inoa,inva,ipoa,ilma,illa)
    call dainf(iib,inob,invb,ipob,ilmb,illb)
    call dainf(iic,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    call damch(ma,ia)
    call damch(mb,ib)
    !
    if(ia.ne.ib) then
       write(line,'(a26)')  'ERROR IN DACCT, IA .NE. IB'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DACCT1',1)
    elseif(ic.ne.invb) then
       write(line,'(a26)')  'ERROR IN DACCT, IC.NE.INVB'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DACCT2',1)
    endif
    !
    !     ALLOCATING LOCAL VECTORS AND CALLING MTREE
    !     ******************************************
    !
    do i=1,ib
       icc(i) = 0
    enddo
    !
    do i=1,c_nomax+1
       mon(i) = 0
    enddo
    !
    call daall(icc,ib,'$$DACCT $$',c_nomax,c_nvmax)
    call daall(mon,c_nomax+1,'$$DAMON $$',inoc,invc)
    !
    call c_mtree(mb,ib,icc,ib)
    !
    !     PERFORMING CONCATENATION
    !     ************************
    !
    do i=1,ia
       call c_dacon(ma(i),c_cc(c_idapo(icc(i))))
    enddo
    !
    call c_dacon(mon(1),(1.0_dp,0.0_dp))
    !
    do i=1,c_idall(icc(1))-1
       !
       jl = c_i_1(c_idapo(icc(1))+i)
       jv = c_i_2(c_idapo(icc(1))+i)
       !
       call c_damul(mon(jl),mc(jv),mon(jl+1))
       !
       do iv=1,ia
          !
          ccf = c_cc(c_idapo(icc(iv))+i)
          if(abs(ccf).gt.eps) call dacma(ma(iv),mon(jl+1),ccf,ma(iv))
          !
       enddo
    enddo
    !
    call dadal(mon,c_nomax+1)
    call dadal(icc,ib)
    !
    return
  end subroutine dacctt
  !
  subroutine c_mtree(mb,ib,mc,ic)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE IS USED FOR CONCATENATION AND TRACKING OF VECTORS
    !     THROUGH A DA MAP. IT COMPUTES THE TREE THAT HAS TO BE TRANSVERSED
    !     MB IS THE DA MATRIX WITH IA TERMS. THE OUTPUT MC IS A CA MATRIX WHICH
    !     CONTAINS COEFFICIENTS AND CONTROL INTEGERS USED FOR THE TRAVERSAL.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ib,ib1,ibi,ic,ic1,ic2,iccx,ichk,ii,iib,iic,illb,illc,ilmb,ilmc,&
         inob,inoc,invb,invc,ipob,ipoc,j,jl,jnon,nterm,ntermf,ipause,mypauses
    integer,dimension(c_lnv)::jjy
    integer,dimension(0:c_lno)::jv
    integer,dimension(:)::mb,mc
    complex(dp) apek,bbijj,chkjj
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !     CONSISTENCY CHECKS
    !     ******************
    !
    iib = mb(1)
    iic = mc(1)
    call dainf(iib,inob,invb,ipob,ilmb,illb)
    call dainf(iic,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    call damch(mb,ib)
    call damch(mc,ic)
    !
    if(ib.ne.ic) then
       write(line,'(a26)')  'ERROR IN MTREE, IB .NE. IC'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR MTREE1',1)
    endif
    !
    !     ALLOCATING LOCAL VECTORS
    !     ************************
    !
    ichk = 0
    call daall1(ichk,'$$MTREE $$',c_nomax,c_nvmax)
    !
    !     FIND ALL THE ENTRIES TO BE LOOKED FOR
    !     *************************************
    !
    call daclr(1)
    !
    c_cc(1) = 1.0_dp
    !
    do i=1,ib
       if(c_nomax.eq.1) then
          do ib1 =2,invc+1    ! 2,7    ! Etienne
             c_cc(ib1) = 1.0_dp
          enddo
       else
          do ibi = c_idapo(mb(i)),c_idapo(mb(i))+c_idall(mb(i))-1
             iccx = c_ia1(c_i_1(ibi)) + c_ia2(c_i_2(ibi))
             if(c_ieo(iccx).gt.inob) goto 90
             c_cc(iccx) = 1.0_dp
90           continue
          enddo
       endif
    enddo
    !
    do ii=1,inob
       !
       !     SEARCHING FOR FATHER FOR EACH TERM
       !
       do i=1,c_nmmax
          if(real(c_cc(i)).lt.0.5_dp) goto 140
          !
          jnon = 0
          call dancd(c_i_1(i),c_i_2(i),jjy)
          do j=1,invb
             if(jjy(j).eq.0) goto 130
             jnon = j
             jjy(j) = jjy(j) - 1
             call dadcd(jjy,ic1,ic2)
             apek = c_cc(c_ia1(ic1)+c_ia2(ic2))
             jjy(j) = jjy(j) + 1
             if(real(apek).ge.0.5_dp) goto 140
130          continue
          enddo
          !
          if(jnon.eq.0) goto 140
          !
          !     TERM IS AN ORPHAN, SO CREATE FOSTER FATHER
          !
          jjy(jnon) = jjy(jnon) - 1
          call dadcd(jjy,ic1,ic2)
          c_cc(c_ia1(ic1)+c_ia2(ic2)) = 1.0_dp
          !
140       continue
       enddo
    enddo
    !
    call dapac(ichk)
    !ETIENNE      CALL DAPRI(ICHK,32)
    !
    !     SETTING UP TREE STRUCTURE
    !     *************************
    !
    ntermf = c_idall(ichk)
    !
    !     ZEROTH ORDER TERMS
    !     ******************
    !
    do i=1,c_lnv
       jjy(i) = 0
    enddo
    !
    do i=1,ib
       call c_dapek(mb(i),jjy,bbijj)
       c_i_1(c_idapo(mc(i))) = 0
       c_i_2(c_idapo(mc(i))) = 0
       c_cc(c_idapo(mc(i))) = bbijj
    enddo
    !
    call c_dapek(ichk,jjy,chkjj)
    if(real(chkjj).gt.0.5_dp) then
       call c_dapok(ichk,jjy,(-1.0_dp,0.0_dp))
    else
       write(line,'(a49)')  'ERROR IN MTREE, ZEROTH ORDER TERM OF ICHK IS 0.0_dp'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR MTREE2',1)
    endif
    !
    nterm = 1
    !
    !     HIGHER ORDER TERMS
    !     ******************
    !
    do jl=1,inob
       jv(jl) = 0
    enddo
    !
    jl = 0
    chkjj = 1.0_dp
    !
200 continue
    if(jl.eq.0.and.(real(chkjj).le.0.5_dp)) goto 250
    if(jl.lt.inob.and.real(chkjj).gt.0.5_dp) then
       jl = jl + 1
       jjy(1) = jjy(1) + 1
       jv(jl) = 1
    elseif(jv(jl).eq.invb) then
       jjy(jv(jl)) = jjy(jv(jl)) - 1
       jv(jl) = 0
       jl = jl - 1
       chkjj = 0.0_dp
       goto 200
    else
       jjy(jv(jl)) = jjy(jv(jl)) - 1
       jv(jl) = jv(jl) + 1
       jjy(jv(jl)) = jjy(jv(jl)) + 1
    endif
    !
    call c_dapek(ichk,jjy,chkjj)
    !
    if(real(chkjj).le.0.5_dp) goto 200
    !
    nterm = nterm + 1
    if(nterm.gt.c_idalm(mc(1))) then
       write(line,'(a31)')  'ERROR IN MTREE, NTERM TOO LARGE'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR MTREE3',1)
    endif
    !
    call c_dapok(ichk,jjy,(-1.0_dp,0.0_dp))
    !
    do i=1,ib
       call c_dapek(mb(i),jjy,bbijj)
       c_i_1(c_idapo(mc(i))+nterm-1) = jl
       c_i_2(c_idapo(mc(i))+nterm-1) = jv(jl)
       c_cc(c_idapo(mc(i))+nterm-1) = bbijj
    enddo
    !
    goto 200
    !
250 continue
    !
    do i=1,ib
       c_idall(mc(i)) = nterm
    enddo
    !
    !     PERFORMING CROSS CHECKS
    !     ***********************
    !
    if(nterm.ne.ntermf.or.nterm.ne.c_idall(ichk)) then
       write(line,'(a46,1x,i8,1x,i8,1x,i8)') 'ERROR IN MTREE, NTERM, NTERMF, c_idall(ICHK) =  ',nterm,ntermf,c_idall(ichk)
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR MTREE4',1)
    endif
    !
    do i=c_idapo(ichk),c_idapo(ichk)+nterm-1
       if(abs(c_cc(i)+1.0_dp).gt.epsmac) then
          write(line,'(a44)')  'ERROR IN MTREE, NOT ALL TERMS IN ICHK ARE -1'
          ipause=mypauses(35,line)
          call dadeb !(31,'ERR MTREE5',1)
       endif
    enddo
    !
    call c_dadal1(ichk)
    !
    return
  end subroutine c_mtree

  !
  subroutine ppushprint(mc,ic,mf,jc,line)
    implicit none
    !
    integer i,ic,iv,jc,jl,jv,mf
    integer,dimension(:)::mc
    character(20) line
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(mf.le.0) return
    write(mf,*) 0,0,jc+1,0,line
    do i=1,ic
       jc=1+jc
       write(mf,*) jc,jl,jv,c_cc(c_idapo(mc(i)))
    enddo
    !     xf(i) = c_cc(c_idapo(mc(i)))
    !      xm(1) = one
    do i=1,c_idall(mc(1))-1
       jl = c_i_1(c_idapo(mc(1))+i)
       jv = c_i_2(c_idapo(mc(1))+i)
       !      xx = xm(jl)*xi(jv)
       !      xm(jl+1) = xx
       do iv=1,ic
          jc=1+jc
          write(mf,*) jc,jl,jv,c_cc(c_idapo(mc(iv))+i)
          !        xf(iv) = xf(iv) + c_cc(c_idapo(mc(iv))+i) * xx
       enddo
    enddo
    return
  end subroutine ppushprint
  !
  subroutine ppushstore(mc,nd2,coef,ml,mv)
    implicit none
    !
    integer i,ic,iv,jc,jl,jv,ntot,nd2
    integer,dimension(:), intent(in)::mc
    integer,dimension(:), intent(out)::ml,mv
    complex(dp),dimension(:),intent(out)::coef
    !
    jc=0
    jl=0;jv=0;
    ic=nd2
    ntot=c_idall(mc(1))*ic
    do i=1,ic
       jc=1+jc
       ml(jc)=0
       mv(jc)=0
       coef(jc)=c_cc(c_idapo(mc(i)))
       !       write(mf,*) jc,jl,jv,c_cc(c_idapo(mc(i)))
    enddo
    do i=1,c_idall(mc(1))-1
       jl = c_i_1(c_idapo(mc(1))+i)
       jv = c_i_2(c_idapo(mc(1))+i)
       !      xx = xm(jl)*xi(jv)
       !      xm(jl+1) = xx
       do iv=1,ic
          jc=1+jc
          ml(jc)=jl
          mv(jc)=jv
          coef(jc)=c_cc(c_idapo(mc(iv))+i)
          !         write(mf,*) jc,jl,jv,c_cc(c_idapo(mc(iv))+i)
          !        xf(iv) = xf(iv) + c_cc(c_idapo(mc(iv))+i) * xx
       enddo
    enddo
    return
  end subroutine ppushstore
  subroutine ppushGETN(mc,ND2,ntot)
    implicit none
    !
    integer ntot,ND2
    integer,dimension(:), intent(inout)::mc
    !
    ntot=c_idall(mc(1))*nd2
  end subroutine ppushGETN
  !
  subroutine ppush(mc,ic,xi,xf)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
    !     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ic,iv,jl,jv
    integer,dimension(:)::mc
    complex(dp) xx
    complex(dp),dimension(c_lno+1)::xm
    complex(dp),dimension(c_lno)::xt
    complex(dp),dimension(:)::xf,xi
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,ic
       xt(i)=xi(i)
    enddo
    do i=1,ic
       xf(i) = c_cc(c_idapo(mc(i)))
    enddo
    !
    xm(1) = 1.0_dp
    !
    do i=1,c_idall(mc(1))-1
       !
       jl = c_i_1(c_idapo(mc(1))+i)
       jv = c_i_2(c_idapo(mc(1))+i)
       xx = xm(jl)*xt(jv)
       xm(jl+1) = xx
       !
       do iv=1,ic
          xf(iv) = xf(iv) + c_cc(c_idapo(mc(iv))+i) * xx
       enddo
    enddo
    do i=1,c_nvmax
       if(abs(xf(i))>da_absolute_aperture.and.check_da) then
          C_STABLE_DA=.false.
          write(6,*) "unstable in ppush ",i,xf(i)
       endif
    enddo
    !
    return
  end subroutine ppush
  subroutine ppush1(mc,xi,xf)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
    !     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
    !
    !-----------------------------------------------------------------------------
    !
    integer i,jl,jv,mc
    complex(dp) xf,xx
    complex(dp),dimension(:)::xi
    complex(dp),dimension(c_lno)::xt
    complex(dp),dimension(c_lno+1)::xm
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,c_nvmax
       xt(i)=xi(i)
    enddo
    xf = c_cc(c_idapo(mc))
    !
    xm(1) = 1.0_dp
    !
    do i=1,c_idall(mc)-1
       !
       jl = c_i_1(c_idapo(mc)+i)
       jv = c_i_2(c_idapo(mc)+i)
       xx = xm(jl)*xt(jv)
       xm(jl+1) = xx
       !
       xf = xf + c_cc(c_idapo(mc)+i) * xx
    enddo
    if(abs(xf)>da_absolute_aperture.and.check_da) then
       C_STABLE_DA=.false.
       write(6,*) "unstable in ppush1 ", xf
       write(6,*) xi(1:c_nvmax)
    endif
    return
  end subroutine ppush1

  subroutine c_dainv(ma,ia,mb,ib)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
    !     STORES THE RESULT IN MI
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob
    integer,dimension(c_lnv)::jj,ml
    integer,dimension(:)::ma,mb
    complex(dp),dimension(c_lnv)::x
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,c_lnv
       jj(i)=0
    enddo
    if(ma(1).eq.mb(1)) then
       call dainf(mb(1),inob,invb,ipob,ilmb,illb)
       if((.not.C_STABLE_DA)) then
          if(C_watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       do i=1,ia
          call c_dapok(ma(i),jj,(0.0_dp,0.0_dp))
       enddo
       do ij=1,ib
          ml(ij)=0
       enddo
       call daall(ml,ib,'$$DAJUNK$$',inob,invb)
       call dainvt(ma,ia,ml,ib)
       do i=1,ib
          call c_dacop(ml(i),mb(i))
       enddo
       call dadal(ml,ib)
    else
       do i=1,ia
          call c_dapek(ma(i),jj,x(i))
          call c_dapok(ma(i),jj,(0.0_dp,0.0_dp))
       enddo
       call dainvt(ma,ia,mb,ib)
       do i=1,ia
          call c_dapok(ma(i),jj,x(i))
       enddo
    endif
    return
  end subroutine c_dainv

  subroutine dainvt(ma,ia,mb,ib)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
    !     STORES THE RESULT IN MI
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ie,ier,illa,illb,ilma,ilmb,inoa,inob,inva,invb,&
         ipoa,ipob,j,k,c_nocut0,ipause,mypauses
    integer,dimension(c_lnv)::jj,ms,ml
    integer,dimension(:)::ma,mb
    complex(dp),dimension(c_lnv,c_lnv)::aa,ai
    complex(dp) amjj,amsjj,prod
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
    call dainf(mb(1),inob,invb,ipob,ilmb,illb)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !     CONSISTENCY CHECKS
    !     ******************
    !
    call damch(ma,ia)
    call damch(mb,ib)
    !etienne
    do ie=1,ib
       call c_dacon(mb(ie),(0.0_dp,0.0_dp))
    enddo
    !etienne
    !
    if(ia.ne.ib) then
       write(line,'(a26)')  'ERROR IN DAINV, IA .NE. IB'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DAINV1',1)
    elseif(ia.ne.inva.or.ib.ne.invb) then
       write(line,'(a40)')  'ERROR IN DAINV, IA.NE.INVA.OR.IB.NE.INVB'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DAINV2',1)
    endif
    !
    !     ALLOCATING LOCAL VECTORS
    !     ************************
    !
    do i=1,ia
       ms(i) = 0
       ml(i) = 0
    enddo
    !
    call daall(ms,ia,'$$INV   $$',inoa,inva)
    call daall(ml,ia,'$$INVL  $$',inoa,inva)
    !
    !     EXTRACTING LINEAR MATRIX, GENERATING NONLINEAR PART OF A
    !     ********************************************************
    !
    do i=1,ib
       do j=1,ib
          do k=1,ib
             jj(k) = 0
          enddo
          jj(j) = 1
          call c_dapek(ma(i),jj,amjj)
          if(abs(amjj).gt.eps) call c_dapok(ma(i),jj,(0.0_dp,0.0_dp))
          aa(i,j) = amjj
       enddo
       call c_dacmu(ma(i),-(1.0_dp,0.0_dp),ma(i))
    enddo
    !
    !     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
    !     **********************************************************
    !
    call c_matinv(aa,ai,ia,c_lnv,ier)
    !
    if(ier.eq.132) then
       if(check_da) then
          C_STABLE_DA=.false.
          C_check_stable=.false.
          messagelost='ERROR IN ROUTINE DAINV, ier=132 in matinv'
          write(6,*) messagelost
          call dadal(ml,ia)
          call dadal(ms,ia)
          return
       else
          write(line,'(a22)')  'ERROR IN ROUTINE DAINV'
          ipause=mypauses(35,line)
          call dadeb !(31,'ERR DAINV3',1)
       endif
    endif
    !
    ier = 0
    do i=1,ib
       do j=1,ib
          prod = 0.0_dp
          do k=1,ib
             jj(k) = 0
             prod = prod + aa(i,k)*ai(k,j)
          enddo
          if(i.eq.j) prod = prod - 1.0_dp
          if(abs(prod).gt.100.0_dp*epsmac) then
             write(6,*) " abs(prod) > 100.0_dp*epsmac in dainvt",abs(prod), 100.0_dp*epsmac
             if(check_da) then
                C_STABLE_DA=.false.
                messagelost='ERROR IN ROUTINE DAINV, abs(prod).gt.100.0_dp*epsmac '
                call dadal(ml,ia)
                call dadal(ms,ia)
                return
             else
                write(line,'(a50,2(1x,i4),3(1x,g13.6))')  'ERROR IN DAINV, INVERSION DID NOT WORK,I,J,PROD = ' &
                     &  ,i,j,prod,epsmac,eps
                ipause=mypauses(35,line)
                ier = 1
                !ETIENNE
                return
                !ETIENNE
             endif
          endif
          jj(j) = 1
          call c_dapok(mb(i),jj,ai(i,j))
          call c_dapok(ml(i),jj,ai(i,j))
       enddo
    enddo
    !
    if(ier.eq.1) call dadeb !(31,'ERR DAINV4',1)
    !
    !     ITERATIVELY COMPUTING DIFFERENT PARTS OF THE INVERSE
    !     ****************************************************
    !
    !     MB (OF ORDER I) = A1^-1 o [ E - ANL (NONLINEAR) o MB (OF ORDER I) ]
    !
    c_nocut0 = c_nocut
    !
    do i=2,c_nocut
       !
       c_nocut = i
       !
       call c_dacct(ma,ia,mb,ib,ms,ia)
       do j=1,ib
          do k=1,ib
             jj(k) = 0
          enddo
          jj(j) = 1
          call c_dapek(ms(j),jj,amsjj)
          call c_dapok(ms(j),jj,amsjj+1.0_dp)
       enddo
       !
       call c_dacct(ml,ia,ms,ia,mb,ib)
       !
    enddo
    !
    c_nocut = c_nocut0
    !
    !     FLIPPING BACK SIGN OF A, FILLING UP FIRST ORDER PART AGAIN
    !     **********************************************************
    !
    do i=1,ib
       call c_dacmu(ma(i),(-1.0_dp,0.0_dp),ma(i))
       do j=1,ib
          do k=1,ib
             jj(k) = 0
          enddo
          jj(j) = 1
          call c_dapok(ma(i),jj,aa(i,j))
       enddo
    enddo
    !
    call dadal(ml,ia)
    call dadal(ms,ia)
    !
    return
  end subroutine dainvt
  !
  subroutine dapin(ma,ia,mb,ib,jx)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
    !     STORES THE RESULT IN MI
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob
    integer,dimension(c_lnv)::jj,ml
    integer,dimension(:)::ma,mb,jx
    complex(dp),dimension(c_lnv)::x
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    do i=1,c_lnv
       jj(i)=0
    enddo
    if(ma(1).eq.mb(1)) then
       call dainf(mb(1),inob,invb,ipob,ilmb,illb)
       if((.not.C_STABLE_DA)) then
          if(C_watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       do i=1,ia
          call c_dapok(ma(i),jj,(0.0_dp,0.0_dp))
       enddo
       do ij=1,ib
          ml(ij)=0
       enddo
       call daall(ml,ib,'$$DAJUNK$$',inob,invb)
       call dapint(ma,ia,ml,ib,jx)
       do i=1,ib
          call c_dacop(ml(i),mb(i))
       enddo
       call dadal(ml,ib)
    else
       do i=1,ia
          call c_dapek(ma(i),jj,x(i))
          call c_dapok(ma(i),jj,(0.0_dp,0.0_dp))
       enddo
       call dapint(ma,ia,mb,ib,jx)
       do i=1,ia
          call c_dapok(ma(i),jj,x(i))
       enddo
    endif
    return
  end subroutine dapin

  subroutine dapint(ma,ia,mb,ib,jind)
    implicit none
    !     **********************************
    !
    !     THIS SUBROUTINE PERFORMS A PARTIAL INVERSION OF THE ROWS MARKED WITH
    !     NONZERO ENTRIES IN JJ OF THE MATRIX A. THE RESULT IS STORED IN B.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,illa,ilma,inoa,inva,ipoa,k
    integer,dimension(c_lnv)::jj,mn,mi,me
    integer,dimension(:)::ma,mb,jind
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,ia
       mn(i) = 0
       mi(i) = 0
       me(i) = 0
    enddo
    !
    call daall(mn,ia,'$$PIN1  $$',inoa,inva)
    call daall(mi,ia,'$$PIN2  $$',inoa,inva)
    call daall(me,ia,'$$PIN3  $$',inoa,inva)
    !
    do i=1,ia
       do k=1,c_nvmax
          jj(k) = 0
       enddo
       jj(i) = 1
       call c_dapok(me(i),jj,(1.0_dp,0.0_dp))
    enddo
    !
    do i=1,ia
       call c_dacop(ma(i),mn(i))
       if(jind(i).eq.0) call c_dacop(me(i),mn(i))
    enddo
    !
    call c_dainv(mn,ia,mi,ia)
    !
    do i=1,ia
       if(jind(i).eq.0) call c_dacop(ma(i),me(i))
    enddo
    !
    call c_dacct(me,ia,mi,ia,mb,ib)
    !
    call dadal(me,ia)
    call dadal(mi,ia)
    call dadal(mn,ia)
    !
    return
  end subroutine dapint
  !
  subroutine c_dader(idif,ina,inc)
    !  subroutine dader(idif,ina,inc)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
    !     OF THE VECTOR A AND STORES THE RESULT IN C.
    !
    !-----------------------------------------------------------------------------
    !
    integer idif,illc,ilmc,ina,inc,incc,inoc,invc,ipoc
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       if((.not.C_STABLE_DA)) then
          if(C_watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       incc=0
       call daall1(incc,'$$DAJUNK$$',inoc,invc)
       call dadert(idif,ina,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call dadert(idif,ina,inc)
    endif
    return
  end subroutine c_dader
  subroutine dadert(idif,ina,inc)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE COMPUTES THE DERIVATIVE WITH RESPECT TO VARIABLE I
    !     OF THE VECTOR A AND STORES THE RESULT IN C.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,illc,&
         ilma,ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj,ipause,mypauses
    integer,dimension(c_lnv)::jd
    complex(dp) rr,x,xdivi
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(c_nomax.eq.1) then
       !         PRINT*,'ERROR, DADER CALLED WITH c_nomax = 1'
       !        call dadeb !(31,'ERR DADER1',1)
       !         stop 666
       do i=1,c_lnv
          jd(i)=0
       enddo
       jd(idif)=1
       call c_dapek(ina,jd,rr)
       call c_dacon(inc,rr)
       return
    endif
    !
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !
    ibase = c_nomax + 1
    !
    if(idif.gt.(c_nvmax+1)/2) then
       ider1  = 0
       ider1s = 0
       ider2  = idif-(c_nvmax+1)/2
       ider2s = 1
       do jj=1,ider2-1
          ider2s = ider2s*ibase
       enddo
       xdivi  = ider2s*ibase
    else
       ider1  = idif
       ider1s = 1
       do jj=1,ider1-1
          ider1s = ider1s*ibase
       enddo
       ider2  = 0
       ider2s = 0
       xdivi  = ider1s*ibase
    endif
    !
    ibase = c_nomax+1
    !
    ic = ipoc-1
    !
    do i=ipoa,ipoa+illa-1
       !
       if(ider1.eq.0) then
          iee = c_i_2(i)
       else
          iee = c_i_1(i)
       endif
       !
       x = iee/xdivi
       ifac = int(ibase*(x-int(x+epsmac)+epsmac))
       !
       if(ifac.eq.0) goto 100
       !
       ic = ic + 1
       c_cc(ic) = c_cc(i)*ifac
       c_i_1(ic) = c_i_1(i) - ider1s
       c_i_2(ic) = c_i_2(i) - ider2s
       !
100    continue
    enddo
    !
    c_idall(inc) = ic - ipoc + 1
    if(c_idall(inc).gt.c_idalm(inc)) then
       write(line,'(a15)') 'ERROR IN DADER '
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DADER2',1)
    endif
    !
    return
  end subroutine dadert
  !
  !
    !
  subroutine c_dacfu(ina,fun,inc)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL complex(dp) FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc
    complex(dp),external::fun
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       if((.not.C_STABLE_DA)) then
          if(C_watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       incc=0
       call daall1(incc,'$$DAJUNK$$',inoc,invc)
       call dacfut(ina,fun,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call dacfut(ina,fun,inc)
    endif
    return
  end subroutine c_dacfu

    !
  subroutine dacfut(ina,fun,inc)
    implicit none
    ! external fun
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL complex(dp) FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,&
         invc,ipoa,ipoc,ipause,mypauses
    integer,dimension(c_lnv)::j
    complex(dp) cfac,rr
    !
    interface
       !       real(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(c_nomax.eq.1) then
       do i=1,c_lnv
          j(i)=0
       enddo
       call c_dapek(ina,j,rr)
       cfac = fun(j)
       rr=cfac*rr
       call c_dapok(inc,j,rr)
       do i=1,c_lnv
          j(i)=1
          call c_dapek(ina,j,rr)
          cfac = fun(j)
          rr=cfac*rr
          call c_dapok(inc,j,rr)
          j(i)=0
       enddo
       !         PRINT*,'ERROR, DACFU CALLED WITH c_nomax = 1'
       !        call dadeb !(31,'ERR DACFU ',1)
       !         stop 667
       return
    endif
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !
    ic = ipoc - 1
    !
    do ia=ipoa,ipoa+illa-1
       !
       call dancd(c_i_1(ia),c_i_2(ia),j)
       cfac = fun(j)
       !      IF(abs(CFAC).LT.EPS) GOTO 100
       !      IF(abs(CFAC*c_cc(IA)).LT.EPS) GOTO 100
       if(abs(cfac*c_cc(ia)).lt.eps.or.abs(c_cc(ia)).lt.eps) goto 100
       !
       ic = ic + 1
       c_cc(ic) = c_cc(ia)*cfac
       c_i_1(ic) = c_i_1(ia)
       c_i_2(ic) = c_i_2(ia)
       !
100    continue
    enddo
    !
    c_idall(inc) = ic - ipoc + 1
    if(c_idall(inc).gt.c_idalm(inc)) then
       write(line,'(a15)') 'ERROR IN DACFU '
       ipause=mypauses(38,line)
       call dadeb !(31,'ERR DACFU ',1)
    endif
    !
    return
  end subroutine dacfut
  !  subroutine GET_C_J(ina,I,C,J)
  !  subroutine GET_C_J(ina,I,C,J)
  !    implicit none
  !
  !    INTEGER I,ina
  !    integer, dimension(c_lnv)::j
  !    complex(dp) C
  !
  !    C=c_cc(I)
  !    call dancd(c_i_1(I),c_i_2(I),J)!
  !  END subroutine GET_C_J
  !  END subroutine GET_C_J
  subroutine c_dapri(ina,iunit)
    implicit none
    !     ***************************
    !       Frank
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,iii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,k
    integer,dimension(c_lnv)::j
    real(dp) a,b
    complex(dp) ccc
    logical  long
     long=longprint
      if(iunit/=6) longprint=.true.
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       print*,'ERROR IN c_dapri, INA = ',ina
       C_STABLE_DA=.false.
    endif
    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    !
    write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)') c_daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,&
         '*********************************************'
    !
    iout = 0
    ioa = 0
    if(inva.eq.0) then
       write(iunit,'(A)') '    I  VALUE  '
       do i = ipoa,ipoa+illa-1
          write(iunit,'(I6,2X,G20.13)') i-ipoa, c_clean_complex(c_cc(i))
       enddo
    elseif(c_nomax.eq.1) then
       if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
       if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS 0.0_dp '
       do i=1,illa
          do k=1,inva
             j(k)=0
          enddo
          iout=iout+1
          if(i.ne.1) then
             j(i-1)=1
             ioa=1
          endif
          write(iunit,'(I6,2X,G20.13,1x,G20.13,I5,4X,18(2i2,1X))') iout,c_clean_complex(c_cc(ipoa+i-1)),ioa,(j(iii),iii=1,c_nvmax)
          write(iunit,*) c_clean_complex(c_cc(ipoa+i-1))
       enddo
    else
       if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
       if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS 0.0_dp '
       do ioa = 0,inoa
          do ii=ipoa,ipoa+illa-1
             if(c_ieo(c_ia1(c_i_1(ii))+c_ia2(c_i_2(ii))).ne.ioa) goto 100
             call dancd(c_i_1(ii),c_i_2(ii),j)
             !ETIENNE
             if(abs(c_cc(ii)).gt.epsprint) then
             a=0; b=0;
             if(abs(real(c_cc(ii)))> epsprint) a=c_cc(ii)
             if(abs(aimag(c_cc(ii)))> epsprint) b=aimag(c_cc(ii))
             ccc=a+(0.0_dp,1.0_dp)*b

                !ETIENNE
                
                iout = iout+1
                write(iunit,'(I6,2X,G20.13,1x,G20.13,I5,4X,18(2i2,1X))') iout,ccc,ioa,(j(iii),iii=1,c_nvmax)
                !ETIENNE
                write(iunit,*) c_cc(ii)
             endif
             !ETIENNE
             !
100          continue
          enddo
       enddo
       !
    endif
    write(iunit,'(A)') '                                      '
    !
    return
longprint=long
  end subroutine c_dapri

function c_clean_complex(c)
implicit none 
complex(dp) c_clean_complex,c
real(dp) cr,ci

cr=c
if(abs(cr)<epsprint) cr=0
ci=-i_*c
if(abs(ci)<epsprint) ci=0
c_clean_complex=cr+i_*ci

end function c_clean_complex 


  subroutine c_dapri77(ina,iunit)
    implicit none
    !     ***************************
    !       Etienne
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit
    integer,dimension(c_lnv)::j
    character(10) c10,k10
    real(dp) a,b
    complex(dp) ccc
    logical some,imprime
    logical  long
     long=longprint
      if(iunit/=6) longprint=.true.
    some=.false.
    !
    if(iunit.eq.0) return
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN c_dapri, INA = ',ina
       C_STABLE_DA=.false.
    endif
    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    !

    if(longprint) then
       write(iunit,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)') c_daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,&
         '*********************************************'
    else
        write(iunit,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)') "Properties",', NO =',inoa,', NV =',inva,', INA =',ina,&
         '*********************************************'
    endif
    !
    if(illa.ne.0.and.longprint) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
    if(illa.eq.0.and.longprint) write(iunit,'(A)') '   ALL COMPONENTS 0.0_dp '
    !
    c10='      NO ='
    k10='      NV ='
    if(longprint)write(iunit,'(A10,I6,A10,I6)') c10,inoa,k10,inva
    iout = 0
    !
    !      DO 100 IOA = 0,INOA
    do ioa = 0,c_nocut
       do ii=ipoa,ipoa+illa-1
          if(c_nomax.ne.1) then
             if(c_ieo(c_ia1(c_i_1(ii))+c_ia2(c_i_2(ii))).ne.ioa) goto 100
          endif
          !ETIENNE
          if(abs(c_cc(ii)).gt.epsprint) then
             !ETIENNE
             a=0; b=0; imprime=.false.
             if(abs(real(c_cc(ii)))> epsprint) then
                 a=real(c_cc(ii))
                 imprime=.true.
            endif
             if(abs(aimag(c_cc(ii)))> epsprint) then
               b=aimag(c_cc(ii))
               imprime=.true.
             endif 
             ccc=a+(0.0_dp,1.0_dp)*b
!             ccc=c_cc(ii)
             if(c_nomax.ne.1) then
                call dancd(c_i_1(ii),c_i_2(ii),j)
                iout = iout+1
             else
                if(ii.eq.ipoa.and.ioa.eq.1) goto 100
                if(ii.gt.ipoa.and.ioa.eq.0) goto 100
                do i=1,c_lnv
                   j(i)=0
                enddo
                if(ii.ne.ipoa) j(ii-ipoa)=1
                iout = iout+1
             endif
             !
             !      WRITE(IUNIT,*) IOA,c_cc(II),(J(I),I=1,INVA)
             if(imprime) then
                 some=.true.
                if(epsprint.gt.1e-37_dp) then
                   write(iunit,501) ioa,ccc,(j(i),i=1,inva)
                else
                   write(iunit,503) ioa,ccc,(j(i),i=1,inva)
                endif
             endif
501          format(' ', i3,1x,g23.16,1x,g23.16,1x,100(1x,i2))
503          format(' ', i3,1x,g23.16,1x,g23.16,1x,100(1x,i2))
502          format(' ', i5,1x,g23.16,1x,g23.16,1x,100(1x,i2))
          endif
          !ETIENNE
          !
100       continue
       enddo
    enddo
    !
    do i=1,c_lnv
       j(i)=0
    enddo
    if(iout.eq.0) iout=1
if(longprint) write(iunit,502) -iout,0.0_dp,0.0_dp,(j(i),i=1,inva)
    if((.not.longprint).and.(.not.some)) write(iunit,*) " Complex Polynomial is zero "
if(.not.longprint) write(6,*) " "
    !
    return
longprint=long
  end subroutine c_dapri77

  subroutine c_dashift(ina,inc,ishift)
    implicit none
    !      complex(dp) c
    !       Frank
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,inb,ishift,ich,&
         ik,inc,k
    integer,dimension(c_lnv)::j,jd
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    inb=0
    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN dashift, INA = ',ina
       C_STABLE_DA=.false.
    endif
    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    call daall1(inb,'$$DAJUNK$$',inoa,inva)
    iout = 0
    !
    !      DO 100 IOA = 0,INOA
    do ioa = 0,c_nocut
       do ii=ipoa,ipoa+illa-1
          if(c_nomax.ne.1) then
             if(c_ieo(c_ia1(c_i_1(ii))+c_ia2(c_i_2(ii))).ne.ioa) goto 100
          endif
          !ETIENNE
          if(abs(c_cc(ii)).gt.eps) then
             !ETIENNE
             if(c_nomax.ne.1) then
                call dancd(c_i_1(ii),c_i_2(ii),j)
                iout = iout+1
             else
                if(ii.eq.ipoa.and.ioa.eq.1) goto 100
                if(ii.gt.ipoa.and.ioa.eq.0) goto 100
                do i=1,c_lnv
                   j(i)=0
                enddo
                if(ii.ne.ipoa) j(ii-ipoa)=1
                iout = iout+1
             endif
             do k=1,ishift   ! put 2004 may
                if(j(k)>0  ) then
                   write(6,*) " trouble in dashift "
                   stop 888
                endif
             enddo
             !
             !      WRITE(IUNIT,*) IOA,c_cc(II),(J(I),I=1,INVA)
             if(abs(c_cc(ii)).gt.eps) then
                if(eps.gt.1e-37_dp) then
                   !       write(iunit,501) ioa,c_cc(ii),(j(i),i=1,inva)
                   ich=1
                   do ik=1,ishift
                      if(j(ik).ne.0) ich=0
                   enddo
                   if(ich.eq.1) then
                      do ik=1,c_lnv
                         jd(ik)=0
                      enddo
                      do ik=ishift+1,c_lnv
                         jd(ik-ishift)=j(ik)  !%%%%%%Etienne
                      enddo
                   endif
                   call c_dapok(inb,jd,c_cc(ii))
                else
                   !       write(iunit,503) ioa,c_cc(ii),(j(i),i=1,inva)
                   ich=1
                   do ik=1,ishift
                      if(j(ik).ne.0) ich=0
                   enddo
                   if(ich.eq.1) then
                      do ik=1,c_lnv
                         jd(ik)=0
                      enddo
                      do ik=ishift+1,c_lnv
                         jd(ik-ishift)=j(ik)  !%%%%%%Etienne
                      enddo
                   endif
                   call c_dapok(inb,jd,c_cc(ii))
                endif
             endif
501          format(' ', i3,1x,g23.16,1x,100(1x,i2))
503          format(' ', i3,1x,g23.16,1x,100(1x,i2))
502          format(' ', i5,1x,g23.16,1x,100(1x,i2))
          endif
          !ETIENNE
          !
100       continue
       enddo
    enddo
    !
    do i=1,c_lnv
       j(i)=0
    enddo
    call c_dacop(inb,inc)
    call c_dadal1(inb)
    !
    return
  end subroutine c_dashift

  subroutine c_darea(ina,iunit)
    implicit none
    !       Frank
    !-----------------------------------------------------------------------------
    !
    integer i,ic,iche,ii,ic_i_1,ic_i_2,iin,illa,ilma,ina,inoa,inva,io,io1,ipoa,iunit,&
         iwarin,iwarno,iwarnv,nno
    integer,dimension(c_lnv)::j
    complex(dp) c
    character(10) c10
 
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       C_STABLE_DA=.false.
       print*,'ERROR IN DAREA, INA = ',ina
       !        X = SQRT(-ONE)
       !        PRINT*,X
    endif

    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    !
    do i=1,c_lnv
       j(i) = 0
    enddo
    !
    call daclr(1)
    call daclr(ina)   ! etienne 2008
    !
    ic = 0
    !
    iwarno = 0
    iwarnv = 0
    iwarin = 0
    !
    read(iunit,'(A10)') c10
    read(iunit,'(18X,I4)') nno
    read(iunit,'(A10)') c10
    read(iunit,'(A10)') c10
    read(iunit,'(A10)') c10
     
    !
    !
    iin = 0
    !
10  continue
    iin = iin + 1
!    read(iunit,'(I6,2X,G20.13,I5,4X,18(2i2,1X))') ii,c,io,(j(i),i=1,inva)
    read(iunit,*) ii,c,io,(j(i),i=1,inva)
    !
    if(ii.eq.0) goto 20
    !ETIENNE
    read(iunit,*) c
    !ETIENNE
    if(ii.ne.iin) then
       iwarin = 1
    endif
    io1 = 0
    do i=1,inva
       io1 = io1 + j(i)
    enddo
    !
    if(io1.ne.io) then
       iwarnv = 1
       goto 10
    endif
    if(io.gt.inoa) then
       !        IF(IWARNO.EQ.0) PRINT*,'WARNING IN DAREA, FILE ',
       !    *              'CONTAINS HIGHER ORDERS THAN VECTOR '
       iwarno = 1
       goto 10
    endif
    !
    if(c_nomax.ne.1) then
       ic = ic + 1
       call dadcd(j,ic_i_1,ic_i_2)
       ic = c_ia1(ic_i_1) + c_ia2(ic_i_2)
       c_cc(ic) = c
       goto 10
    else
       iche=0
       do i=1,inva
          if(j(i).eq.1) iche=i
       enddo
       c_cc(ipoa+iche)=c
       goto 10
    endif
    !
20  continue
    !
    if(c_nomax.ne.1) call dapac(ina)
    !
    return
 
  end subroutine c_darea
  !FF
  !
  subroutine c_darea77(ina,iunit)
    implicit none
    !     ***************************
    !     Etienne
    !     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ic,iche,ii,ic_i_1,ic_i_2,iin,illa,ilma,ina,inoa,inva,ipoa,iunit,&
         k,nojoh,nvjoh
    integer,dimension(c_lnv)::j
    real(dp) cr,ci
    character(10) c10,k10
    complex ik
 
    ik=( 0.0_dp,1.0_dp )
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN darea77, INA = ',ina
       C_STABLE_DA=.false.
    endif
    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    !
    do i=1,c_lnv
       j(i) = 0
    enddo
    !
    call daclr(1)
    call daclr(ina)   ! etienne 2008
    !
    !

    read(iunit,'(A10)') c10
    read(iunit,'(A10)') c10
    read(iunit,'(A10)') c10
    read(iunit,'(A10)') c10

    read(iunit,'(A10)') c10
     read(iunit,'(A10,I6,A10,I6)') c10,nojoh,k10,nvjoh
    !
    iin = 0
    !
10  continue
    iin = iin + 1
    read(iunit,*) ii,cr,ci,(j(k),k=1,nvjoh)
    if(ii.lt.0) goto 20
    do i=inva+1,nvjoh
       if(j(i).ne.0) goto 10
    enddo
    iche=0
    do i=1,inva
       iche=iche+j(i)
    enddo
    if(iche.gt.c_nomax) goto 10
    if(c_nomax.ne.1) then
       call dadcd(j,ic_i_1,ic_i_2)
       ic = c_ia1(ic_i_1) + c_ia2(ic_i_2)
       c_cc(ic) = cr+ik*ci
    else
       iche=0
       do i=1,inva
          if(j(i).eq.1) iche=i
       enddo
       c_cc(ipoa+iche)=cr+ik*ci
    endif
    goto 10
    !
20  continue
    !
    if(c_nomax.ne.1) call dapac(ina)
 
    return
  end subroutine c_darea77

  subroutine dadeb   !(iunit,c,istop)
    implicit none
    !     *******************************
    !
    !     THIS SUBROUTINE SERVES AS A DEBUGGING TOOL. IT PRINTS ALL
    !     NONZERO INFORMATION IN THE COMMON BLOCKS AND ALL DA  VECTORS.
    !
    !-----------------------------------------------------------------------------
    !
    !integer,dimension(0:1)::i8
    !
    !etienne
    !    if(check_da) then
    C_STABLE_DA=.false.
 
    write(6,*) "big problem in complex dadeb ", sqrt(crash)
    return
    !    endif
    !READ(*,*) I
    !I=SQRT(DBLE(I))
    !    stop
  end subroutine dadeb
  !
  !
  !
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
    if(inc.ge.1.and.inc.le.c_nda_dab) then
       inoc = c_idano(inc)
       invc = c_idanv(inc)
       ipoc = c_idapo(inc)
       ilmc = c_idalm(inc)
       illc = c_idall(inc)
       return
    endif
    !
    write(line,'(a26,1x,i8,1x,a11)')  'ERROR IN DAINF, DA VECTOR ',inc,' NOT FOUND '
    ipause=mypauses(35,line)
    call dadeb !(31,'ERR DAINF ',1)
    !
    return
  end subroutine dainf
  !
  subroutine dapac(inc)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE PACKS THE INFORMATION IN THE SCRATCH VECTOR 1
    !     INTO THE VECTOR INC. IF LF = 1, THE FILTERING (CF DAMUF) IS
    !     PERFORMED.
    !     INVERSE IS DAUNP.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ic,inc,ipoc,ipause,mypauses
    complex(dp) ccc
    !
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    ipoc = c_idapo(inc)
    !
    ic = ipoc - 1
    !
    do i=1,c_nmmax
       ccc = c_cc(i)
       if(abs(ccc).lt.eps) goto 100
       ic = ic + 1
       c_cc(ic) = ccc
       c_i_1(ic) = c_ie1(i)
       c_i_2(ic) = c_ie2(i)
100    continue
    enddo
    !
    c_idall(inc) = ic - ipoc + 1
    if(c_idall(inc).gt.c_idalm(inc)) then
       write(line,'(a15)')  'ERROR IN DAPAC '
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DAPAC ',1)
    endif
    !
    return
  end subroutine dapac
  !
  !
  subroutine dachk(ina,inoa,inva, inb,inob,invb, inc,inoc,invc)
    implicit none
    !     *************************************************************
    !
    !     THIS SUBROUTINE CHECKS IF THE VECTORS A, B AND C
    !     HAVE COMPATIBLE ATTRIBUTES
    !
    !-----------------------------------------------------------------------------
    !
    integer ierr,ina,inb,inc,inoa,inob,inoc,inva,invb,invc,&
         invsum,ipause,mypauses
    !
    if(lsw.eq.1) return
    !
    ierr = 0
    !
    !     CASE OF A UNARY OPERATION
    !     *************************
    !
    if(inob.eq.-1.and.invb.eq.-1) then
       invsum = inva + invc
       if(invsum.eq.0) then
          if(inoa.gt.inoc) ierr = 1
       elseif(invsum.eq.1) then
          ierr = 1
       else
          if(inoa.gt.inoc.or.inva.gt.invc) ierr = 1
       endif
       if(ierr.eq.1) then
          write(line,'(a16,i8,a5,i8,a17,4(1x,i8))')  'ERROR IN DACHK, ',ina,' AND ',inc, &
               & ' ARE INCOMPATIBLE',inoa,inva,inoc,invc
          ipause=mypauses(35,line)
          call dadeb !(31,'ERR DACHK1',1)
       endif
       !
       !     CASE OF A BINARY OPERATION
       !     **************************
       !
    else
       invsum = inva + invb + invc
       if(invsum.eq.0) then
          if(inoa.gt.inoc.or.inob.gt.inoc) ierr = 1
       elseif(invsum.eq.1.or.invsum.eq.2) then
          ierr = 1
       else
          if(inoa.gt.inoc.or.inob.gt.inoc.or.inva.gt.invc.or.invb.gt.invc) ierr = 1
       endif
       if(ierr.eq.1) then
          write(line,'(a16,i8,a1,i8,a5,i8,a17)')  'ERROR IN DACHK, ',ina,',',inb &
               & ,' AND ',inc,' ARE INCOMPATIBLE'
          ipause=mypauses(35,line)
          call dadeb !(31,'ERR DACHK2',1)
       endif
    endif
    !
    return
  end subroutine dachk
  !
  subroutine damch(iaa,ia)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE CHECKS IF THE IA VECTORS IN THE MATRIX IA HAVE
    !     IDENTICAL ATTRIBUTES.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,illa,ilma,ino1,inoi,inv1,invi,ipoa,ipause,mypauses
    integer,dimension(:)::iaa
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(iaa(1),ino1,inv1,ipoa,ilma,illa)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=2,ia
       call dainf(iaa(i),inoi,invi,ipoa,ilma,illa)
       if((.not.C_STABLE_DA)) then
          if(C_watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       if(ino1.ne.inoi.or.inv1.ne.invi) then
          write(line,'(a24,i8,a5,i8,a18)')  'ERROR IN DAMCH, VECTORS ',iaa(1), &
               & ' AND ',iaa(i),' ARE INCOMPATIBLE '
          ipause=mypauses(35,line)
          stop
       endif
    enddo
    !
    return
  end subroutine damch
  !
  subroutine dadcd(jj,ic1,ic2)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES c_i_1,c_i_2.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic1,ic2,isplit
    integer,dimension(c_lnv)::jj
    !
    ibase = c_nomax + 1
    isplit = (c_nvmax+1)/2
    ic1 = 0
    ic2 = 0
    !
    do i=c_nvmax,isplit+1,-1
       ic2 = ic2*ibase + jj(i)
    enddo
    !
    do i=isplit,1,-1
       ic1 = ic1*ibase + jj(i)
    enddo
    !
    return
  end subroutine dadcd
  !
  subroutine dancd(ic1,ic2,jj)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES c_i_1,c_i_2.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ic1,ic2,isplit
    integer,dimension(:)::jj
    real(dp) x
    !
    ibase = c_nomax + 1
    isplit = (c_nvmax+1)/2
    !
    ic = ic1
    do i=1,isplit
       x  = REAL(ic,kind=DP)/REAL(ibase,kind=DP)
       ic = int(x+epsmac)
       jj(i) = nint(ibase*(x-ic))
    enddo
    !
    ic = ic2
    do i=isplit+1,c_nvmax
       x  = REAL(ic,kind=DP)/REAL(ibase,kind=DP)
       ic = int(x+epsmac)
       jj(i) = nint(ibase*(x-ic))
    enddo
    !
    do i=c_nvmax+1,size(jj)    ! 2002.12.2
       jj(i) = 0
    enddo
    !
    return
  end subroutine dancd
  !ETIENNE
  subroutine c_datra(idif,ina,inc)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE COMPUTES THE PSEUDO DERIVATIVE WITH RESPECT TO VARIABLE I
    !     OF THE VECTOR A AND STORES THE RESULT IN C.
    !
    !     dx^n/dx= x^(n-1)
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ider1,ider1s,ider2,ider2s,idif,iee,ifac,illa,illc,ilma,&
         ilmc,ina,inc,inoa,inoc,inva,invc,ipoa,ipoc,jj,ipause,mypauses
    complex(dp) x,xdivi
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    !       CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !
    if(c_nomax.eq.1) then
       call c_dader(idif,ina,inc)
       return
    endif
    ibase = c_nomax + 1
    !
    if(idif.gt.(c_nvmax+1)/2) then
       ider1  = 0
       ider1s = 0
       ider2  = idif-(c_nvmax+1)/2
       ider2s = 1
       do jj=1,ider2-1
          ider2s = ider2s*ibase
       enddo
       xdivi  = ider2s*ibase
    else
       ider1  = idif
       ider1s = 1
       do jj=1,ider1-1
          ider1s = ider1s*ibase
       enddo
       ider2  = 0
       ider2s = 0
       xdivi  = ider1s*ibase
    endif
    !
    ibase = c_nomax+1
    !
    ic = ipoc-1
    !
    do i=ipoa,ipoa+illa-1
       !
       if(ider1.eq.0) then
          iee = c_i_2(i)
       else
          iee = c_i_1(i)
       endif
       !
       x = iee/xdivi
       ifac = int(ibase*(x-int(x+epsmac)+epsmac))
       !
       if(ifac.eq.0) goto 100
       !
       !etienne      IFAC = INT(IBASE*(X-INT(X)+c_1d_8))
       !
       !etienne      IF(IFAC.EQ.0) GOTO 100
       !
       ic = ic + 1
       c_cc(ic) = c_cc(i)
       c_i_1(ic) = c_i_1(i) - ider1s
       c_i_2(ic) = c_i_2(i) - ider2s
       !
100    continue
    enddo
    !
    c_idall(inc) = ic - ipoc + 1
    if(c_idall(inc).gt.c_idalm(inc)) then
       write(line,'(a16)')  'ERROR IN DADTRA '
       ipause=mypauses(35,line)
       call dadeb !(111,'ERR DADTRA',1)
    endif
    !
    return
  end subroutine c_datra

  subroutine hash(no1,nv1,jj,ic1,ic2)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES c_i_1,c_i_2.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic1,ic2,isplit,no1,nv1
    integer,dimension(:)::jj
    !
    ibase = no1 + 1
    isplit = (nv1+1)/2
    ic1 = 0
    ic2 = 0
    !
    do i=nv1,isplit+1,-1
       ic2 = ic2*ibase + jj(i)
    enddo
    !
    do i=isplit,1,-1
       ic1 = ic1*ibase + jj(i)
    enddo
    !
    return
  end subroutine hash
  !
  subroutine dehash(no1,nv1,ic1,ic2,jj)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES c_i_1,c_i_2.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ic1,ic2,isplit,no1,nv1
    integer,dimension(:)::jj
    real(dp) x
    !
    !frs epsmac=c_1d_7
    ibase = no1 + 1
    isplit = (nv1+1)/2
    !
    ic = ic1
    do i=1,isplit
       x  = REAL(ic,kind=DP)/REAL(ibase,kind=DP)
       ic = int(x+epsmac)
       jj(i) = nint(ibase*(x-ic))
    enddo
    !
    ic = ic2
    do i=isplit+1,nv1
       x  = REAL(ic,kind=DP)/REAL(ibase,kind=DP)
       ic = int(x+epsmac)
       jj(i) = nint(ibase*(x-ic))
    enddo
    !
    return
  end subroutine dehash


  !
  subroutine c_dacycle(ina,ipresent,value,illa,j)
    implicit none
    integer ipause, mypauses
    !     ***************************
    !
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer ii,illa,ilma,ina,inoa,inva,iout,ipoa,ipresent
    integer,optional,dimension(:)::j
    complex(dp) value
    !
    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN dacycle, INA = ',ina
       ipause=mypauses(39,line)
       write(6,*) ina
       stop
    endif
    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    if(.not.present(j)) return
    iout = 0
    if(ipresent.gt.illa.or.ipresent<1) then
       write(6,*) ipresent,illa
       write(6,*) " error in dacycle "
       read(5,*) ipresent
       stop 101
    endif
    ii=ipresent+ipoa-1
    call dancd(c_i_1(ii),c_i_2(ii),j)
    value=c_cc(ii)
    if(c_nomax==1) then
       j=0
       if(ipresent/=1) j(ipresent-1)=1
    endif
    !    ipresent=1+ipresent
    return
  end subroutine c_dacycle
!!!! new stuff lingyun
  subroutine dacc_lean(ina)
    implicit none
    integer ipause, mypauses
    !     ***************************
    !
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer ii,illa,ilma,ina,inoa,inva,iout,ipoa
    real(dp) a,b
 

    if(ina.lt.1.or.ina.gt.c_nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN dacycle, INA = ',ina
       ipause=mypauses(39,line)
       stop
    endif
    !
    inoa = c_idano(ina)
    inva = c_idanv(ina)
    ipoa = c_idapo(ina)
    ilma = c_idalm(ina)
    illa = c_idall(ina)
    iout = 0
    do ii=ipoa,ipoa+ilma-1
  !     if(abs(c_cc(ii))<value) c_cc(ii)=0.0_dp
             a=0; b=0;
             if(abs(real(c_cc(ii)))> eps) a=real(c_cc(ii))
             if(abs(aimag(c_cc(ii)))> eps) b=aimag(c_cc(ii))
             c_cc(ii)=a+(0.0_dp,1.0_dp)*b
    !   
    enddo
    !    call dapac(ina)
    return
  end subroutine dacc_lean


subroutine c_dafun(cf,ina,inc)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
    !     AND STORES THE RESULT IN C.
    !     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
    !     THIS HAS TO BE FIXED IN THE FUTURE.
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inc,incc
    character(4) cf
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       !       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1(incc,'$$DAJUNK$$',c_nomax,c_nvmax)
       call dafunt(cf,ina,incc)
       call c_dacop(incc,inc)
       call c_dadal1(incc)
    else
       call dafunt(cf,ina,inc)
    endif
    return
  end subroutine c_dafun

  subroutine dafunt(cf,ina,inc)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE COMPUTES THE FUNCTION CF OF THE DA VECTOR A
    !     AND STORES THE RESULT IN C.
    !     AT PRESENT, SOME FUNCTIONS CAN BE COMPUTED ONLY TO FIFTH ORDER.
    !     THIS HAS TO BE FIXED IN THE FUTURE.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ina,inc,ind,inon,ipow,iscr,lfun,no,ipause,mypauses
    integer,parameter,dimension(c_lnv)::jjx=0
    complex(dp) a0,ca,ea,ra,sa
    real(dp) ang
    complex(dp),dimension(0:c_lno)::xf
    character(4) cf,cfh
    character(26) abcs,abcc
    !
    data abcs /'abcdefghijklmnopqrstuvwxyz'/
    data abcc /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
    !
    if((.not.C_STABLE_DA)) then
       if(C_watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(cf(1:1).eq.' ') then
       cfh(1:3) = cf(2:4)
       cfh(1:4) = ' '
       cf = cfh
    endif
    !
    do i=1,4
       ind = index(abcs,cf(i:i))
       if(ind.ne.0) cf(i:i) = abcc(ind:ind)
    enddo
    !      call dainf(ina,inoa,inva,ipoa,ilma,illa)
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    !
    !     CASE OF NV = 0 WHICH MEANS COORDINATEWISE OPERATION
    !     ***************************************************
    !
    !     CASE OF NV > 0 WHICH MEANS DIFFERENTIAL ALGEBRAIC OPERATION
    !     ***********************************************************
    !
    if(cf.eq.'SQR ') then
       call dasqr(ina,inc)
       return
    endif
    !
    !     ALLOCATE VARIABLES, PICK ZEROTH ORDER TERM
    !     ******************************************
    !
    ipow = 0
    inon = 0
    iscr = 0
    !
    call daall1(ipow,'$$DAFUN1$$',c_nomax,c_nvmax)
    call daall1(inon,'$$DAFUN2$$',c_nomax,c_nvmax)
    call daall1(iscr,'$$DAFUN3$$',c_nomax,c_nvmax)
    !
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !
    call c_dapek(ina,jjx,a0)
    !
    !      no = min(c_nocut,inoa,inoc)
    no = min(c_nocut,c_nomax)
    !
    !     BRANCHING TO DIFFERENT FUNCTIONS
    !     ********************************
    !
    select case(cf)
    case('INV ')
       !    if(cf.eq.'INV ') then
       !        1/(A0+P) = 1/A0*(1-(P/A0)+(P/A0)**2-...)
       if(abs(a0).eq.0) then
          if(check_da) then
             messagelost="a0.eq.0 for INV in dafun"
             C_STABLE_DA=.false.
             C_check_stable=.false.
             call c_dadal1(iscr)
             call c_dadal1(inon)
             call c_dadal1(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       xf(0) = 1.0_dp/a0
       do i=1,no
          xf(i) = -xf(i-1)/a0
       enddo
       !
        !
    case('LOG ')
       !    elseif(cf.eq.'LOG ') then
       !        LOG(A0+P) = LOG(A0) + (P/A0) - 1/2*(P/A0)**2 + 1/3*(P/A0)**3 - ...)

       ea  = LOG(abs(a0)) + (0.0_dp,1.0_dp)*atan2(aimag(a0),real(a0))

       xf(0) = ea
       xf(1) = 1.0_dp/a0
       do i=2,no
          xf(i) = -xf(i-1)/a0/REAL(i,kind=DP)*REAL(i-1,kind=DP)
       enddo
       !

    case('SIN ')
       !    elseif(cf.eq.'SIN ') then
       !        SIN(A0+P) = SIN(A0)*(1-P**2/2!+P**4/4!) + COS(A0)*(P-P**3/3!+P**5/5!)
       sa  = SIN(a0)
       ca  = COS(a0)
       xf(0) = sa
       xf(1) = ca
       do i=2,no
          xf(i) = -xf(i-2)/REAL(i*(i-1),kind=DP)
       enddo
       !
    case('COS ')
       !    elseif(cf.eq.'COS ') then
       !        COS(A0+P) = COS(A0)*(1-P**2/2!+P**4/4!) - SIN(A0)*(P-P**3/3!+P**5/5!)
       sa  = SIN(a0)
       ca  = COS(a0)
       xf(0) = ca
       xf(1) = -sa
       do i=2,no
          xf(i) = -xf(i-2)/REAL(i*(i-1),kind=DP)
       enddo
       !

    case('EXP ')
       !    elseif(cf.eq.'EXP ') then
       !        EXP(A0+P) = EXP(A0)*(1+P+P**2/2!+...)
       if(abs(a0)>hyperbolic_aperture) then
          if(check_da) then
             messagelost="a0>hyperbolic_aperture for EXP in dafun"
             C_STABLE_DA=.false.
             C_check_stable=.false.
             call c_dadal1(iscr)
             call c_dadal1(inon)
             call c_dadal1(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       ea  = exp(a0)
       xf(0) = ea
       do i=1,no
          xf(i) = xf(i-1)/REAL(i,kind=DP)
       enddo
    case('SQRT')
       !    elseif(cf.eq.'SQRT') then
       !        SQRT(A0+P) = SQRT(A0)*(1+1/2(P/A0)-1/8*(P/A0)**2+...)

!       ra = SQRT(a0)

        ang=atan2(aimag(a0),real(a0))/2.0_dp
        ra=(cos(ang)+(0.0_dp,1.0_dp)*sin(ang))*sqrt(abs(a0)) 
       xf(0) = ra
       do i=1,no
          xf(i) = -xf(i-1)/a0/REAL(2*i,kind=DP)*REAL(2*i-3,kind=DP)
       enddo
       !
       !

    case default
       !    else
       write(line,'(a28,1x,a4)')  'ERROR, UNSOPPORTED FUNCTION ',cf
       ipause=mypauses(35,line)
       !    endif
    end select
    !
    call c_dacon(inc,xf(0))
    call c_dacop(ina,inon)
    call c_dapok(inon,jjx,(0.0_dp,0.0_dp))
    call c_dacon(ipow,(1.0_dp,0.0_dp))
    !
    do i=1,min(no,c_nocut)
       !
       call c_damul(inon,ipow,iscr)
       call c_dacop(iscr,ipow)
       call dacma(inc,ipow,xf(i),inc)
       !
    enddo
    !
1000 format('ERROR IN DAFUN, ',a4,' DOES NOT EXIST FOR VECTOR ',i10,'CONST TERM  = ',e12.5)
    !
    call c_dadal1(iscr)
    call c_dadal1(inon)
    call c_dadal1(ipow)
    !
    return
  end subroutine dafunt

  subroutine c_take(h,m,ht)
    implicit none
    !  HT= H_M  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
    integer i,m,h,ht,b1,b2,b3,no,nv
    integer,dimension(c_lnv)::j
    complex(dp) r
    if(.not.c_stable_da) return
    
    no = c_nomax
    nv = c_nvmax
    call c_etall1(b1)
    call c_etall1(b2)
    call c_etall1(b3)

    if(no.ge.2) then
       if(m.eq.0) then
          do i=1,c_lnv
             j(i)=0
          enddo
          call c_dapek(h,j,r)
          call c_dacon(ht,r)
       else
          !          call danot(m)
          !          call dacop(h,b1)
          call c_datrunc(h,m+1,b1)
          !          call danot(m-1)
          !          call dacop(b1,b2)
          call c_datrunc(b1,m,b2)
          !          call danot(no)
          call c_dasub(b1,b2,b3)
          call c_dacop(b3,ht)
       endif
    else
       do i=1,c_lnv
          j(i)=0
       enddo
       if(m.eq.0) then
          call c_dapek(h,j,r)
          call c_dacon(ht,r)
       elseif(m.eq.1)  then
          do i=1,nv
             j(i)=1
             call c_dapek(h,j,r)
             call c_dapok(b3,j,r)
             j(i)=0
          enddo
          call c_dacop(b3,ht)
       else
          call daclr(ht)
       endif
    endif

    call c_dadal1(b3)
    call c_dadal1(b2)
    call c_dadal1(b1)
    return
  end subroutine c_take

 subroutine c_daran(ina,cm,xran)
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE FILLS THE DA VECTOR A WITH RANDOM ENTRIES.
    !     FOR CM > 0, THE VECTOR IS FILLED WITH REALS,
    !     FOR CM < 0, THE VECTOR IS FILLED WITH SINGLE DIGIT INTEGERS
    !     ABS(CM) IS THE FILLING FACTOR
    !
    !-----------------------------------------------------------------------------
    !
    integer i,illa,ilma,ina,inoa,inva,ipoa,ipause,mypauses
    real(dp) cm,xran
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(inva.eq.0.or.c_nomax.eq.1) then
       do i=ipoa,ipoa+ilma-1
          if(cm.gt.0.0_dp) then
             c_cc(i) = bran(xran)+(0.0_dp,1.0_dp)*bran(xran)
             c_cc(i) =c_cc(i)/abs(c_cc(i)) 
             if(abs(c_cc(i)).gt.cm) c_cc(i) = 0.0_dp
          elseif(cm.lt.0.0_dp) then
             c_cc(i) = int(1+10*bran(xran))
             if(abs(c_cc(i)).gt.-10.0_dp*cm) c_cc(i) = 0.0_dp
          endif
       enddo
       c_idall(ina) = c_idalm(ina)
       return
    endif
    !
    if(inoa.ne.c_nomax.or.inva.ne.c_nvmax) then
       line='ERROR IN DARAN, ONLY VECTORS WITH NO = c_nomax AND NV = c_nvmax ALLOWED'
       ipause=mypauses(31,line)
       call dadeb !(31,'ERR DARAN1',1)
    endif
    !
    call daclr(1)
    !
    do i=1,c_nmmax
       if(cm.gt.0.0_dp) then
          c_cc(i) =  bran(xran)+(0.0_dp,1.0_dp)*bran(xran)
             c_cc(i) =c_cc(i)/abs(c_cc(i)) 
          if(abs(c_cc(i)).gt.cm) c_cc(i) = 0.0_dp
       elseif(cm.lt.0.0_dp) then
          c_cc(i) = int(1+10*bran(xran))
          if(abs(c_cc(i)).gt.-10.0_dp*cm) c_cc(i) = 0.0_dp
       else
          line='ERROR IN ROUTINE DARAN'
          ipause=mypauses(31,line)
          call dadeb !(31,'ERR DARAN2',1)
       endif
    enddo
    !
    call dapac(ina)
    !
    return
  end subroutine c_daran

  subroutine c_dapek0(v,x,jj)
    implicit none

    integer i,jj
    integer,dimension(c_lnv)::jd
    integer,dimension(:)::v
    complex(dp),dimension(:)::x
    if(.not.c_stable_da) return

    do i=1,c_lnv
       jd(i)=0
    enddo
    do i=1,jj
       call c_dapek(v(i),jd,x(i))
    enddo
    return
  end subroutine c_dapek0

  subroutine c_dapok0(v,x,jj)
    implicit none
    integer i,jj
    integer,dimension(c_lnv)::jd
    integer,dimension(:)::v
    complex(dp),dimension(:)::x
    if(.not.c_stable_da) return

    do i=1,c_lnv
       jd(i)=0
    enddo
    do i=1,jj
       call c_dapok(v(i),jd,x(i))
    enddo
    return
  end subroutine c_dapok0

  subroutine c_etcom(x,y,h,nd2)
    implicit none
    ! ETCOM TAKES THE BRACKET OF TWO VECTOR FIELDS.
    integer i,j,t1,t2,nd2
    integer,dimension(:)::h,x,y
    integer,dimension(c_lnv)::t3
    if(.not.c_stable_da) return
    
    call c_etall1(t1)
    call c_etall1(t2)
    do j=1,nd2
     call c_etall1(t3(j))
    enddo


    do j=1,nd2
       do i=1,nd2

          call c_dader(i,x(j),t1)
          call c_dader(i,y(j),t2)
          call c_damul(x(i),t2,t2)
          call c_damul(y(i),t1,t1)
          call dalin(t2,(1.0_dp,0.0_dp),t1,(-1.0_dp,0.0_dp),t1)
          call c_daadd(t1,t3(j),t3(j))

       enddo
    enddo


    do j=1,nd2
     call c_dacop(t3(j),h(j))
    enddo

    call c_dadal1(t1)
    call c_dadal1(t2)
    do j=1,nd2
     call c_dadal1(t3(j))
    enddo

    return
  end subroutine c_etcom

end module c_dabnew
