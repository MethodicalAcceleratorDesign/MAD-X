! The Full Polymorphic Package
! The module in this file is, to the best of our knowledge,
! the property of Lawrence Berkeley National Laboratory
! Its distribution and commercial usage may therefore be governed by the laws of the
! United States of America

!$$$$ module dabnew
module dabnew_b  !$$$$
  use da_arrays
  implicit none
  !  private
  public
  private daallno1,daall,damult,dasqrt,dacmut,dacma,DALINt,dafunt,dacctt
  private dainvt,dapint,dadert,dacfuRt,dacfuIt,dacfut
  private dadeb,dapac,dachk,damch,dadcd,dancd,hash,dehash
  !  public ppushstore,ppushprint,ppushGETN
  !  public DEALLOC_ALL,dainf
  !  PUBLIC DAINI,daall1
  !  PUBLIC DANOT,DAVAR,PPUSH,MTREE,DADAL,DACCT,DADER,DASQR
  !  PUBLIC DAMUL,DAADD,DACOP,DAINV,DAPIN,DAPEK,DAPOK,DACFU,DASUB,DACON,DADAL1
  !  PUBLIC DACLR,DACMU,DALIN,DAREA,DAPRI,DAABS,DAEPS,DATRA,MATINV,DAFUN
  !  PUBLIC DADIV,DADIC,DACDI,DACAD,DACSU,DASUC,DASHIFT,DARAN,DACFUR
  !  PUBLIC DACFUI,DAPRI77,DAREA77,GET_C_J,PPUSH1,DALLSTA,dacycle
  !  public count_da
  integer,private,parameter:: lsw=1
  integer :: lda_max_used=0

  ! integer,private,parameter::nmax=400,lsw=1
  ! real(dp),private,parameter::tiny=c_1d_20
  character(120),private :: line
  private dacsu_b !$$$$
contains
  !******************************************************************************
  !                                                                             *
  !                                                                             *
  !                                                                             *
  !               DIFFERENTIAL ALGEBRA PACKAGE OF M. BERZ                       *
  !                     ****************************                            *
  !                               holy3                                         *
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
  !     DAINI(O,N,U):       INITIALIZES CONTROL ARRAYS AND SETS MAX. ORDER O AND
  !                         MAX. NUMBER OF VARIABLES N. MUST BE CALLED BEFORE ANY
  !                         OTHER DA ROUTINE CAN BE USED.
  !
  !     DAALL(A,I,T,O,N):   ALLOCATES SPACE FOR I VECTORS A. T: CHARACTER NAME
  !     DADAL(A,I):         DEALLOCATES THE I VECTORS A.
  !!     DAVAR(A,R,I):       MAKES A INDEPENDENT VARIABLE # I WITH INITIAL VALUE R
  !!     DACON(A,R):         SETS A TO CONSTANT R
  !     DANOT(O):           SETS NEW TRUNCATION ORDER O FOR DA OPERATIONS
  !     DAEPS(R):           SETS NEW ZERO TOLERANCE EPSILON
  !
  !!     DAPEK(A,JJ,R):      RETURNS COEF R OF MONOMIAL WITH EXPONENTS JJ OF A
  !!     DAPOK(A,JJ,R):      SETS COEF OF MONOMIAL WITH EXPONENTS JJ OF A TO R
  !
  !!     DACOP(A,C):         PERFORMS C = A
  !!     DAADD(A,B,C):       PERFORMS C = A + B
  !!    DASUB(A,B,C):       PERFORMS C = A - B
  !!     DAMUL(A,B,C):       PERFORMS C = A * B
  !!     DADIV(A,B,C):       PERFORMS C = A / B
  !!     DASQR(A,C):         PERFORMS C = A^2           (SQUARE OF A)
  !
  !!     DACAD(A,RA,C):      PERFORMS C = A + RA
  !!     DACSU(A,RA,C):      PERFORMS C = A - RA
  !!     DASUC(A,RA,C):      PERFORMS C = RA - A
  !!     DACMU(A,RA,C):      PERFORMS C = A * RA
  !!    DACDI(A,RA,C):      PERFORMS C = A / RA
  !!     DADIC(A,RA,C):      PERFORMS C = RA / A
  !!     DACMA(A,B,RB,C):    PERFORMS C = A + RB*B
  !!DAMULIN(A,B,RA,C,D,RB,C):    PERFORMS C = A*B*RA + C*D*RB
  !!     DALIN(A,RA,B,RB,C): PERFORMS C = A*RA + B*RB
  !!     DAFUN(F,A,C):       PERFORMS C = F(A)          (DA FUNCTION)
  !
  !!     DAABS(A,R):         PERFORMS R = |A|           (NORM OF A)
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
  !!     DACFU(A,G,C):       MULTIPLIES COEFFICIENTS WITH FUNCTION G(JJ)
  !
  !     DAPRI(A,U):         PRINTS DA VECTOR A TO UNIT U
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
  !     - CHECK DAREA AND DAPRI FOR CA VECTORS
  !     - MAKE DAFUN USE DASQR
  !
  !
  !      BLOCKDATA DABLD
  !     ***************
  !
  !
  !     PARAMETERS:
  !
  !     LDA: MAXIMUM NUMBER OF DA-VECTORS;    CAN BE CHANGED QUITE ARBITRARILY
  !     LST: LENGTH OF MAIN STORAGE STACK;    CAN BE CHANGED QUITE ARBITRARILY
  !     LEA: MAXIMUM NUMBER OF MONOMIALS;     CAN BE INCREASED FOR LARGE NO,NV
  !     LIA: DIMENSION OF IA1,IA2;            CAN BE INCREASED FOR LARGE NO,NV
  !     LNO: MAXIMUM ORDER;                   CAN BE INCREASED TO ABOUT 1000
  !     LNV: MAXIMUM NUMBER OF VARIABLES;     CAN BE INCREASED TO ABOUT 1000
  !
  !     ALL THE CHANGES IN THE VALUES OF PARAMETERS HAVE TO BE MADE BY GLOBAL
  !     SUBSTITUTIONS IN ALL SUBROUTINES.
  !
  !     DANAME:   NAME OF DA VECTOR
  !
  !     CC:       STACK OF DOUBLE PRECISON COEFFICIENTS
  !     i_1:       FIRST CHARACTERISTIC INTEGER (CF DAINI)
  !     i_2:       SECOND CHARACTERISTIC INTEGER (CF DAINI)
  !
  !     IE1:      CHARACTERISTIC INTEGER 1 OF UNPACKED REPRESENTATION (CF DAINI)
  !     IE2:      CHARACTERISTIC INTEGER 2 OF UNPACKED REPRESENTATION (CF DAINI)
  !     IEO:      ORDER OF ENTRY IN UNPACKED REPRESENTATION
  !     IA1:      REVERSE TO IE1 (CF DAINI)
  !     IA2:      REVERSE TO IE2 (CF DAINI)
  !
  !     IDANO:    ORDER OF DA VECTOR; IN CA, NUMBER OF COMPONENTS
  !     IDANV:    NUMBER OF VARIABLES; IF 0, INDICATES CA VECTOR
  !     IDAPO:    FIRST ADDRESS IN STACK
  !     IDALM:    NUMBER OF RESERVED STACK POSITIONS
  !     IDALL:    NUMBER OF MOMENTARILY REQUIRED STACK POSITIONS
  !
  !     nda_dab:      NUMBER OF DA VECTORS MOMENTARILY DEFINED
  !     NST:      NUMBER OF STACK POSITIONS MOMENTARILY ALLOCATED
  !     NOMAX:    MAXIMUM REQUESTED ORDER  (CF DAINI)
  !     NVMAX:    MAXIMUM REQUESTED NUMBER OF VARIABLES (CF DAINI)
  !     NMMAX:    MAXIMUM NUMBER OF MONOMIALS FOR NOMAX, NVMAX (CF DAINI)
  !     NOCUT:    MOMENTARY TRUNCATION ORDER
  !     EPS:      TRUNCATION ACCURACY (CAN BE SET BY USER)
  !     EPSMAC:   MANTISSA LENGTH OF MACHINE (PESSIMISTIC ESTIMATE)
  !
  !-----------------------------------------------------------------------------
  !
  subroutine change_package(i)
    implicit none
    integer, intent(in) :: i
    if(i==2) then
       lingyun_yang=.false.
    elseif(i==1) then
       lingyun_yang=.true.
    else
       write(6,*) " i = 1 or 2 "
       write(6,*) " INPUT IGNORED "
    endif
  end  subroutine change_package

  !$$$$  subroutine daini(no,nv,iunit)
  subroutine daini_b(no,nv,iunit)  !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE SETS UP THE MAJOR ORDERING AND ADDRESSING ARRAYS IN
    !     COMMON BLOCK DAINI. IF IUNIT > 0, THE ARRAYS WILL BE PRINTED TO UNIT
    !     NUMBER IUNIT. AN EXAMPLE FOR THE ARRAYS GENERATED BY DAINI CAN BE
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
    integer,dimension(lnv+1)::n
    integer,dimension(0:lnv)::k
    integer,dimension(lnv)::j,jj
    character(10) aa
    !
    !frs if(eps.le.zero) eps=c_1d_38
    !      if(EPS.le.zero) eps=c_1d_90
    !frs epsmac=c_1d_7
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    last_tpsa=2
    if(nv.eq.0) return

    ! LD: trace
    ! write(6,*) "no=", no, "nv=", nv

    call alloc_all(no,nv)
    ndamaxi=0
    !
    do i=1, lda
       allvec(i) = .false.
    enddo
    nhole=0
    !*****************************************

    !     INITIALIZING VARIABLES IN COMMON / DA /
    !     ***************************************
    !
    nda_dab   = 0
    nst0   = 0
    nomax = no
    nvmax = nv
    call danum(no,nv,nmmax)
    nocut = no
    lfi   = 0
    !
    do i=0,lia
       ia1(i) = 0
       ia2(i) = 0
    enddo
    !
    !    do i=1,100
    !       is(i) = 0
    !    enddo
    !
    if(nv.gt.lnv.or.no.gt.lno) then
       write(6,*) 'ERROR IN SUBROUTINE DAINI, NO, NV = ',no,nv
       ipause=mypauses(1,line)
       call dadeb !(31,'ERR DAINI ',1)
    endif
    !
    ibase = no+1
    js    = nv/2
    if(float(ibase)**((nv+1)/2).gt.float(lia)) then
       write(line,'(a12,i4,a7,i4,a21,i4)') 'ERROR, NO = ',no,', NV = ',nv,' TOO LARGE FOR LIA = ',lia
       ipause=mypauses(2,line)
       call dadeb !(31,'ERR DAINI ',1)
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
          ia2(ic2) = nn
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
                if(etiennefix) then
                   if(nn<=lea) then   ! Etienne
                      ie1(nn) = ic1
                      ie2(nn) = ic2
                   endif             ! Etienne
                   if(nn<=lst) then   ! Etienne
                      i_1 (nn) = ic1
                      i_2 (nn) = ic2
                   endif
                   if(ic2.eq.0) ia1(ic1) = nn
                   if(nn<=lea) then   ! Etienne
                      ieo(nn) = io1 + io2
                   endif             ! Etienne
                   !
                else
                   ie1(nn) = ic1
                   ie2(nn) = ic2
                   i_1 (nn) = ic1
                   i_2 (nn) = ic2
                   if(ic2.eq.0) ia1(ic1) = nn
                   ieo(nn) = io1 + io2
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
    if(nn.gt.lea.and.(.not.etiennefix)) then
       write(line,'(a21,i4,a12)') 'ERROR IN DAINI, NN = ',nn,' EXCEEDS LEA'
       ipause=mypauses(3,line)
       call dadeb !(31,'ERR DAINI ',1)
    endif
    !
    !     ALLOCATING SCRATCH VARIABLES
    !     ****************************
    !
    iall = 0
    call daall1_b(iall,'$$UNPACK$$',nomax,nvmax)
    !
    do i=0,nomax
       aa = '$$MUL   $$'
       write(aa(6:10),'(I5)') i
       iall = 0
       !      CALL DAALL(IALL,1,AA,I,NVMAX)
       call daall1_b(iall,aa,nomax,nvmax)
    enddo
    !
    idall(1) = nmmax
    !
    !     DOUBLE CHECKING ARRAYS IE1,IE2,IA1,IA2
    !     **************************************
    !
    do i=1,nmmax
       !
       jjj = ia1(ie1(i)) + ia2(ie2(i))
       if(jjj.ne.i) then
          write(line,'(a48,i4)') 'ERROR IN DAINI IN ARRAYS IE1,IE2,IA1,IA2 AT I = ',i
          ipause=mypauses(4,line)
          call dadeb !(31,'ERR DAINI ',1)
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
    open(iout,file='DAINI.DAT',status='NEW')
    !CRAY OPEN(IOUT,FILE='DAINI',STATUS='UNKNOWN',FORM='FORMATTED')          *CRAY
    !CRAY REWIND IOUT                                                        *CRAY
    !
    write(iout,'(/A/A/)') ' ARRAYS i_1 THROUGH i_20, IE1,IE2,IEO **********************************'
    do i=1,nmmax
       call dancd(ie1(i),ie2(i),jj)
       write(iout,'(1X,I5,2X,4(5i2,1X),3I6)') i,(jj(jjjj),jjjj=1,lnv),ie1(i),ie2(i),ieo(i)
    enddo
    !
    write(iout,'(/A/A/)') ' ARRAYS IA1,IA2 **************'
    do i=0,icmax
       write(iout,'(3i10)') i,ia1(i),ia2(i)
    enddo
    !
    return
  end subroutine daini_b !$$$$
  !$$$$  end subroutine daini

  !  subroutine daexter
  !    implicit none
  !     *****************************
  !
  !-----------------------------------------------------------------------------
  !     integer i_1,i_2,ia1,ia2,idall,idalm,idano,idanv,idapo,ie1,ie2,ieo,ifi,lfi,nda,ndamaxi,nmmax,nocut,nomax,nst,nvmax
  !-----------------------------------------------------------------------------
  !
  !    integer i
  !
  !    if(.not.notallocated) then
  !       do i=1, lda
  !          allvec(i)=.false.
  !       enddo
  !    endif

  !    return
  !  end subroutine daexter

  subroutine dallsta(ldanow)
    implicit none
    !     *****************************
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ldanow
    !
    ldanow=0
    do i=1, lda
       if(allvec(i)) ldanow=ldanow+1
    enddo

    return
  end subroutine dallsta

  !
  ! EXAMPLE: ARRAYS i_1 THROUGH i_20, IE1,IE2,IEO (NOMAX=3,NVMAX=4)
  ! *************************************************************
  !     I   i_1               THROUGH               i_20     IE1   IE2   IEO
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
  !    ARRAYS IA1,IA2
  !    **************
  !    I        IA1       IA2
  !    0         1         0   IE1,IE2 AND IA1,IA2 ALLOW THE EASY COMPUTATION
  !    1         2        10   OF THE ADDRESS OF THE PRODUCT OF TWO MONOMIALS.
  !    2         4        22   LET IX AND IY BE THE POSITIONS OF THE TWO
  !    3         7        31   FACTORS. THEN THE POSITION IZ OF THE PRODUCT OF
  !    4         3        16   THE TWO FACTORS IS GIVEN BY
  !    5         5        25
  !    6         8        32   IZ = IA1(IE1(IX)+IE1(IY)) + IA2(IE2(IX)+IE2(IY))
  !    7         0         0
  !    8         6        28
  !    9         9        33   THE OTHER VARIABLES SET BY DAINI WOULD HAVE THE
  !   10         0         0   VALUES
  !   11         0         0
  !   12        10        34   NOMAX = 3,  NVMAX = 4, NMMAX = 35
  !

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
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    no=nomax
    nv=nvmax
    ind = 1
    if(ic.gt.0.and.ic.le.nda_dab) then
       !         DANAME(IC(I)) = C
       !         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' NOMAX, NVMAX = ',nomax,nvmax
          ipause=mypauses(7,line)
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(nhole.gt.0) then
          ind=nda_dab
20        if (allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          nhole=nhole-1
       else
          incnda = .true.
          nda_dab = nda_dab + 1
          ind=nda_dab
          if(nda_dab.gt.lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(8,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
       endif

       if(ind>lda_max_used) lda_max_used=ind
       if(ind>lda) then
          write(6,*) "ind>lda ",lda,ind
          print*, 'ERROR IN DAALLNO1, MAX NUMBER OF DA VECTORS EXHAUSTED: LDA = ',LDA
          stop
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
             idapo(ind) = nst0 + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst0 = nst0 + nmmax
          else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst0 + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst0 = nst0 + ndanum
          endif
       endif
       !
       if(nst0.gt.lst) then
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(nv.eq.0.or.nomax.eq.1) then
          call daclr_b(ic)
          idall(ic) = idalm(ic)
       endif
    endif

    !
    if(nda_dab.gt.ndamaxi) ndamaxi=nda_dab

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
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    ind = 1

    do i=1,l
       if(ic(i).gt.0.and.ic(i).le.nda_dab) then
          !         DANAME(IC(I)) = C
          !         IF(IDANO(IC(I)).EQ.NO.AND.IDANV(IC(I)).EQ.NV) THEN
       else
          if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
             write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
                  &' NOMAX, NVMAX = ',nomax,nvmax
             ipause=mypauses(9,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
          !
          if(nhole.gt.0) then
             ind=nda_dab
20           if (allvec(ind)) then
                ind = ind - 1
                goto 20
             endif
             incnda = .false.
             nhole=nhole-1
          else
             incnda = .true.
             nda_dab = nda_dab + 1
             ind=nda_dab
             if(nda_dab.gt.lda) then
                write(6,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
                !    ipause=mypauses(10,line)
                call dadeb !(31,'ERR DAALL ',1)
                stop 111
             endif
          endif
          !write(30,*) no,ind,lda,size(allvec)
          if(ind>lda_max_used) lda_max_used=ind
          if(ind>lda) then
             write(6,*) "ind>lda ",lda,ind
             print*, 'ERROR IN DAALLNO1, MAX NUMBER OF DA VECTORS EXHAUSTED: LDA = ',LDA
             stop
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
                idapo(ind) = nst0 + 1
                idalm(ind) = nmmax
                idall(ind) = 0
                nst0 = nst0 + nmmax
             else
                idano(ind) = no
                idanv(ind) = nv
                idapo(ind) = nst0 + 1
                idalm(ind) = ndanum
                idall(ind) = 0
                nst0 = nst0 + ndanum
             endif
          endif
          !
          if(nst0.gt.lst) then
             call dadeb !(31,'ERR DAALL ',1)
          endif
          !
          !          IF(NV.EQ.0) THEN
          if(nv.eq.0.or.nomax.eq.1) then
             call daclr_b(ic(i))
             idall(ic(i)) = idalm(ic(i))
          endif
       endif
    enddo
    !
    if(nda_dab.gt.ndamaxi) ndamaxi=nda_dab

    return
  end subroutine daall

  subroutine daall1_b(ic,ccc,no,nv) !$$$$
    !$$$$  subroutine daall1(ic,ccc,no,nv)
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
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    ind = 1

    if(ic.gt.0.and.ic.le.nda_dab) then
       !         DANAME(ic) = C
       !         IF(IDANO(ic).EQ.NO.AND.IDANV(ic).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' NOMAX, NVMAX = ',nomax,nvmax
          ipause=mypauses(11,line)
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(nhole.gt.0) then
          ind=nda_dab
20        if (allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          nhole=nhole-1
       else
          incnda = .true.
          nda_dab = nda_dab + 1
          ind=nda_dab
          if(nda_dab.gt.lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(12,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
       endif

       if(ind>lda_max_used) lda_max_used=ind
       if(ind>lda) then
          write(6,*) "ind>lda ",lda,ind
          print*, 'ERROR IN DAALLNO1, MAX NUMBER OF DA VECTORS EXHAUSTED: LDA = ',LDA
          stop
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
             idapo(ind) = nst0 + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst0 = nst0 + nmmax
          else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst0 + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst0 = nst0 + ndanum
          endif
       endif
       !
       if(nst0.gt.lst) then
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       !          IF(NV.EQ.0) THEN
       if(nv.eq.0.or.nomax.eq.1) then
          call daclr_b(ic)
          idall(ic) = idalm(ic)
       endif
    endif
    !
    if(nda_dab.gt.ndamaxi) ndamaxi=nda_dab

    return
    !$$$$  end subroutine daall1
  end subroutine daall1_b  !$$$$

  !$$$$  subroutine daall0(ic)
  subroutine daall0_b(ic) !$$$$
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
    no=nomax
    nv=nvmax
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    !
    ind = 1

    if(ic.gt.0.and.ic.le.nda_dab) then
       !         DANAME(ic) = C
       !         IF(IDANO(ic).EQ.NO.AND.IDANV(ic).EQ.NV) THEN
    else
       if(nv.ne.0.and.(no.gt.nomax.or.nv.gt.nvmax)) then
          write(line,'(a23,i4,a14,i4,1x,i4,a16,i4,1x,i4)') 'ERROR IN DAALL, VECTOR ',c,' HAS NO, NV = ',no,nv, &
               &' NOMAX, NVMAX = ',nomax,nvmax
          ipause=mypauses(11,line)
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       if(nhole.gt.0) then
          ind=nda_dab
20        if (allvec(ind)) then
             ind = ind - 1
             goto 20
          endif
          incnda = .false.
          nhole=nhole-1
       else
          incnda = .true.
          nda_dab = nda_dab + 1
          ind=nda_dab
          if(nda_dab.gt.lda) then
             write(line,'(a50)') 'ERROR IN DAALL, MAX NUMBER OF DA VECTORS EXHAUSTED'
             ipause=mypauses(12,line)
             call dadeb !(31,'ERR DAALL ',1)
          endif
       endif

       if(ind>lda_max_used) lda_max_used=ind
       if(ind>lda) then
          write(6,*) "ind>lda ",lda,ind
          print*, 'ERROR IN DAALLNO1, MAX NUMBER OF DA VECTORS EXHAUSTED: LDA = ',LDA
          stop
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
             idapo(ind) = nst0 + 1
             idalm(ind) = nmmax
             idall(ind) = 0
             nst0 = nst0 + nmmax
          else
             idano(ind) = no
             idanv(ind) = nv
             idapo(ind) = nst0 + 1
             idalm(ind) = ndanum
             idall(ind) = 0
             nst0 = nst0 + ndanum
          endif
       endif
       !
       if(nst0.gt.lst) then
          call dadeb !(31,'ERR DAALL ',1)
       endif
       !
       !          IF(NV.EQ.0) THEN
       if(nv.eq.0.or.nomax.eq.1) then
          call daclr_b(ic)
          idall(ic) = idalm(ic)
       endif
    endif
    !
    if(nda_dab.gt.ndamaxi) ndamaxi=nda_dab

    return
    !$$$$  end subroutine daall0
  end subroutine daall0_b !$$$$

  !
  !$$$$  subroutine dadal(idal,l)
  subroutine dadal_b(idal,l) !$$$$
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
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    do i=l,1,-1
       if(idal(i).le.nomax+2.or.idal(i).gt.nda_dab) then
          write(line,'(a38,i8,1x,i8)') 'ERROR IN ROUTINE DADAL, IDAL(I),NDA = ',idal(i),nda_dab
          !ipause=mypauses(13,line)
          C_%STABLE_DA = .false.
          l = 1
          return
          
          call dadeb !(31,'ERR DADAL ',1)
       endif
       if(idal(i).eq.nda_dab) then
          !       deallocate
          nst0 = idapo(nda_dab) - 1
          nda_dab = nda_dab - 1
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
    !$$$$  end subroutine dadal
  end subroutine dadal_b !$$$$

  !$$$$  subroutine dadal1(idal)
  subroutine dadal1_b(idal) !$$$$
    implicit none
    !     ************************
    !
    !     THIS SUBROUTINE DEALLOCATES THE VECTORS IDAL
    !
    !-----------------------------------------------------------------------------
    !
    integer idal,ipause,mypauses
    !
    !    if((.not.C_%STABLE_DA)) then
    !       if(c_%watch_user) then
    !          write(6,*) "big problem in dabnew ", sqrt(crash)
    !       endif
    !       return
    !    endif
    if(idal.le.nomax+2.or.idal.gt.nda_dab) then
       write(line,'(a35,i8,1x,i8)') 'ERROR IN ROUTINE DADAL, IDAL,NDA = ',idal,nda_dab
       ipause=mypauses(14,line)
       call dadeb !(31,'ERR DADAL ',1)
    endif
    if(idal.eq.nda_dab) then
       !       deallocate
       nst0 = idapo(nda_dab) - 1
       nda_dab = nda_dab - 1
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
    !$$$$  end subroutine dadal1
  end subroutine dadal1_b !$$$$

  !$$$$  subroutine count_da(n)
  subroutine count_da_b(n) !$$$$
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
    !$$$$  end subroutine count_da
  end subroutine count_da_b !$$$$

  !$$$$  subroutine davar(ina,ckon,i)
  subroutine davar_b(ina,ckon,i) !$$$$
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE DECLARES THE DA VECTOR
    !     AS THE INDEPENDENT VARIABLE NUMBER I.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic1,ic2,illa,ilma,ina,inoa,inva,ipoa,ipause,mypauses
    real(dp) ckon
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
    !
    if(i.gt.inva) then
       write(line,'(a20,i8,a16,i8)') 'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
       ipause=mypauses(14,line)
       call dadeb !(31,'ERR DAVAR ',1)
    endif
    !
    if(nomax.eq.1) then
       if(i.gt.inva) then
          print*,'ERROR IN DAVAR, I = ',i,' EXCEEDS INVA = ',inva
          !           call dadeb !(31,'ERR DAVAR3',1)
       endif
       call daclr_b(ina)
       cc(ipoa) = ckon
       cc(ipoa+i) = one
       return
    endif

    ibase = nomax+1
    !
    if(i.gt.(nvmax+1)/2) then
       ic1 = 0
       ic2 = ibase**(i-(nvmax+1)/2-1)
    else
       ic1 = ibase**(i-1)
       ic2 = 0
    endif
    !
    if(abs(ckon).gt.eps) then
       idall(ina) = 2
       cc(ipoa) = ckon
       i_1(ipoa) = 0
       i_2(ipoa) = 0
       !
       cc(ipoa+1) = one
       i_1(ipoa+1) = ic1
       i_2(ipoa+1) = ic2
    else
       idall(ina) = 1
       cc(ipoa) = one
       i_1(ipoa) = ic1
       i_2(ipoa) = ic2
    endif
    !
    return
    !$$$$  end subroutine davar
  end subroutine davar_b !$$$$
  !
  !$$$$  subroutine dacon(ina,ckon)
  subroutine dacon_b(ina,ckon) !$$$$
    implicit none
    !     **************************
    !
    !     THIS SUBROUTINE SETS THE VECTOR C TO THE CONSTANT CKON
    !
    !-----------------------------------------------------------------------------
    !
    integer illa,ilma,ina,inoa,inva,ipoa
    real(dp) ckon
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
    !
    if(nomax.eq.1) then
       call daclr_b(ina)
       cc(ipoa) = ckon
       return
    endif

    idall(ina) = 1
    cc(ipoa) = ckon
    i_1(ipoa) = 0
    i_2(ipoa) = 0
    if(abs(ckon).lt.eps) idall(ina) = 0
    !
    return
    !$$$$  end subroutine dacon
  end subroutine dacon_b !$$$$
  !
  !$$$$  subroutine danot(not)
  subroutine danot_b(not) !$$$$
    implicit none
    !     *********************
    !
    !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
    !
    !-----------------------------------------------------------------------------
    !
    integer not,ipause,mypauses
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(not.gt.nomax) then
       write(line,'(a15,i8,a17,i8)') 'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
       ipause=mypauses(15,line)
       call dadeb !(31,'ERR DANOT ',1)
    endif
    !
    nocut = not
    !
    return
  end subroutine danot_b !$$$$
  !$$$$  end subroutine danot

  !  subroutine getdanot(not)
  !    implicit none
  !     *********************
  !
  !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
  !
  !-----------------------------------------------------------------------------
  !
  !    integer not,ipause,mypauses
  !
  !    if(not.gt.nomax) then
  !       write(line,'(a15,i8,a17,i8)') 'ERROR, NOCUT = ',nocut,' EXCEEDS NOMAX = ',nomax
  !       ipause=mypauses(15,line)
  !       call dadeb !(31,'ERR DANOT ',1)
  !    endif
  !
  !    not=nocut
  !
  !    return
  !  end subroutine getdanot

  !$$$$  subroutine daeps(deps)
  subroutine daeps_b(deps) !$$$$
    implicit none
    !     **********************
    !
    !     THIS SUBROUTINE RESETS THE TRUNCATION ORDER NOCUT TO A NEW VALUE
    !
    !-----------------------------------------------------------------------------
    !
    real(dp) deps
    !
    if(deps.ge.zero) then
       eps = deps
    else
       deps=eps
    endif
    !
    return
    !$$$$  end subroutine daeps
  end subroutine daeps_b !$$$$
  !
  !$$$$  subroutine dapek(ina,jv,cjj)
  subroutine dapek_b(ina,jv,cjj) !$$$$
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE DETERMINES THE COEFFICIENT OF THE ARRAY
    !     OF EXPONENTS JJ AND RETURNS IT IN CJJ
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ic1,ic2,icu,icz,ii_1,illa,ilma,ina,inoa,inva,ipek,ipoa,&
         iu,iz,jj1,mchk,ipause,mypauses
    integer,dimension(:)::jv     ! 2002.12.4
    integer,dimension(lnv)::jj
    real(dp) cjj
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    jj=0
    do i=1,size(jv)
       jj(i)=jv(i)
    enddo


    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
    if(inva.eq.0.or.nomax.eq.1) then
       if(inva.ne.0.and.nomax.eq.1) then
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
          print*,'ERROR IN DAPEK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
          !           call dadeb !(31,'ERR DAPEK1',1)
       endif
       ipek = ipoa + jj1 - 1
       cjj = cc(ipek)
       return
    endif

    ii_1 = (nvmax+1)/2
    ibase = nomax+1
    !
    !     DETERMINE INDEX TO BE SEARCHED FOR
    !     **********************************
    !
    call dadcd(jj,ic1,ic2)
    !
    !ETIENNE
    if(ic1.gt.lia.or.ic2.gt.lia) then
       write(line,'(a24,i8)') 'DISASTER IN DAPEK, INA= ',ina
       ipause=mypauses(16,line)
    endif
    !ETIENNE
    ic = ia1(ic1) + ia2(ic2)
    !
    !     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
    !     ****************************************************************
    !
    !      IF(ICO.GT.INOA.OR.ICV.GT.INVA.OR.ICO.GT.NOCUT) THEN
    !         CJJ = 0
    !         RETURN
    !      ENDIF
    !
    !     DETERMINE IF MONOMIAL IS INSIDE FIRST AND LAST MONOMIALS OF A
    !     *************************************************************
    !
    iu = ipoa
    iz = ipoa + illa - 1
    icu = ia1(i_1(iu))+ia2(i_2(iu))
    icz = ia1(i_1(iz))+ia2(i_2(iz))
    !
    if(illa.eq.0) then
       cjj = 0
       return
    elseif(ic.eq.icu) then
       cjj = cc(iu)
       return
    elseif(ic.eq.icz) then
       cjj = cc(iz)
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
    !     if(ia1(i_1(i))+ia2(i_2(i)) - ic) 20,30,40
    mchk=ia1(i_1(i))+ia2(i_2(i)) - ic
    if(mchk.lt.0) goto 20
    if(mchk.eq.0) goto 30
    if(mchk.gt.0) goto 40
20  iu = i
    goto 10
30  cjj = cc(i)
    return
40  iz = i
    goto 10
    !
    !$$$$  end subroutine dapek
  end subroutine dapek_b !$$$$
  !
  !$$$$  subroutine dapok(ina,jv,cjj)
  subroutine dapok_b(ina,jv,cjj) !$$$$
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
    integer,dimension(lnv)::jj
    real(dp) cjj
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    jj1 = 1
    if(inva.eq.0.or.nomax.eq.1) then
       if(inva.ne.0.and.nomax.eq.1) then
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
          print*,'ERROR IN DAPOK, INDEX OUTSIDE RANGE, JJ(1) = ',jj1
          !           call dadeb !(31,'ERR DAPOK1',1)
       endif
       ipok = ipoa + jj1 - 1
       cc(ipok) = cjj
       return
    endif

    !     DETERMINE INDEX TO BE SEARCHED FOR
    !     **********************************
    !
    call dadcd(jj,ic1,ic2)
    !
    ic = ia1(ic1) + ia2(ic2)
    !
    !     DETERMINE IF MONOMIAL TO BE POKED CONFORMS WITH INOA, INVA,NOCUT
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
       icu = ia1(i_1(iu))+ia2(i_2(iu))
       icz = ia1(i_1(iz))+ia2(i_2(iz))
    endif

    if(illa.eq.0) then
       i = ipoa
       goto 100
    elseif(ic.eq.icu) then
       cc(iu) = cjj
       i = iu
       goto 200
    elseif(ic.eq.icz) then
       cc(iz) = cjj
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
    !      if(ia1(i_1(i))+ia2(i_2(i)) - ic) 20,30,40
    mchk=ia1(i_1(i))+ia2(i_2(i)) - ic
    if(mchk.lt.0) goto 20
    if(mchk.eq.0) goto 30
    if(mchk.gt.0) goto 40
20  iu = i
    goto 10
30  cc(i) = cjj
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
       cc(ii) = cc(ii-1)
       i_1(ii) = i_1(ii-1)
       i_2(ii) = i_2(ii-1)
    enddo
    !
    cc(i) = cjj
    i_1(i) = ic1
    i_2(i) = ic2
    !
    idall(ina) = illa + 1
    if(idall(ina).gt.idalm(ina)) then
       write(line,'(a15)') 'ERROR IN DAPOK '
       ipause=mypauses(17,line)
       call dadeb !(31,'ERR DAPOK ',1)
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
          cc(ii) = cc(ii+1)
          i_1(ii) = i_1(ii+1)
          i_2(ii) = i_2(ii+1)
       enddo
       idall(ina) = illa - 1
    endif
    return
    !
    !$$$$  end subroutine dapok
  end subroutine dapok_b !$$$$
  !
  !$$$$  subroutine daclr(inc)
  subroutine daclr_b(inc) !$$$$
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
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
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
    !$$$$  end subroutine daclr
  end subroutine daclr_b !$$$$
  !
  !$$$$  subroutine dacop(ina,inb)
  subroutine dacop_b(ina,inb) !$$$$
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ipob = idapo(inb)
    ipoa = idapo(ina)
    illa = idall(ina)
    ib = ipob - 1
    !
    !      iif = 0
    !      if(nomax.eq.1.or.inva.eq.0) iif = 1

    do ia = ipoa,ipoa+illa-1
       !
       if(nomax.gt.1) then
          if(ieo(ia1(i_1(ia))+ia2(i_2(ia))).gt.nocut) goto 100
       endif
       ib = ib + 1
       cc(ib) = cc(ia)
       i_1(ib) = i_1(ia)
       i_2(ib) = i_2(ia)
       !
100    continue
    enddo
    !
    idall(inb) = ib - ipob + 1
    return
    !$$$$  end subroutine dacop
  end subroutine dacop_b !$$$$


  !$$$$  subroutine daadd(ina,inb,inc)
  subroutine daadd_b(ina,inb,inc) !$$$$
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i)   + cc(ipob+i)
       enddo
       return
    endif

    if(ina.ne.inc.and.inb.ne.inc) then
       call dalin_b(ina,+one,inb,+one,inc)
    else
       idaadd = 0
       call daall1_b(idaadd,'$$DAADD $$',nomax,nvmax)
       call dalin_b(ina,+one,inb,+one,idaadd)
       call dacop_b(idaadd,inc)
       call dadal1_b(idaadd)
    endif
    !
    return
    !$$$$  end subroutine daadd
  end subroutine daadd_b !$$$$
  !
  !$$$$    subroutine datrunc(ina,io,inb)
  subroutine datrunc_b(ina,io,inb) !$$$$
    implicit none
    integer ina,io,inb,nt
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    nt=nocut
    if(io>nomax) then
       if(ina/=inb) call dacop_b(ina,inb)
       return
    endif
    nocut=io-1

    call dacop_b(ina,inb)

    nocut = nt


    !$$$$    end subroutine datrunc
  end subroutine datrunc_b !$$$$

  !$$$$  subroutine dasub(ina,inb,inc)
  subroutine dasub_b(ina,inb,inc) !$$$$
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i)   - cc(ipob+i)
       enddo
       return
    endif
    if(ina.ne.inc.and.inb.ne.inc) then
       call dalin_b(ina,+one,inb,-one,inc)
    else
       idasub = -1
       !         call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       call daall1_b(idasub,'$$DASUB $$',nomax,nvmax)
       call dalin_b(ina,+one,inb,-one,idasub)
       call dacop_b(idasub,inc)
       call dadal1_b(idasub)
    endif
    !
    return
    !$$$$  end subroutine dasub
  end subroutine dasub_b !$$$$

  !$$$$  subroutine damul(ina,inb,inc)
  subroutine damul_b(ina,inb,inc) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
    !     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inb,inc,incc,ipoc,ipoa,ipob,i
    real(dp) ccipoa,ccipob
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif

    ! LD: count number of calls
    ! call poly_mul_incr

    if(nomax.eq.1) then
       ipoa=idapo(ina)
       ipob=idapo(inb)
       ipoc=idapo(inc)
       !         minv = min(inva,invb,invc)
       ccipoa = cc(ipoa)
       ccipob = cc(ipob)
       cc(ipoc) = ccipoa*ccipob
       do i=1,nvmax
          cc(ipoc+i) = ccipoa*cc(ipob+i) + ccipob*cc(ipoa+i)
       enddo
       !         do 30 i=ipoc+minv+1,ipoc+invc
       !  30     cc(i) = zero
       return
    endif

    if(ina.eq.inc.or.inb.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',nomax,nvmax)
       call damult(ina,inb,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call damult(ina,inb,inc)
    endif

    return
    !$$$$  end subroutine damul
  end subroutine damul_b !$$$$

  subroutine damult(ina,inb,inc)
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE PERFORMS A DA MULTIPLICATION OF THE DA VECTORS A AND B.
    !     THE RESULT IS STORED IN C. AS TEMPORARY STORAGE, THE STACK SPACE
    !     OF THE (NOMAX+2) SCRATCH VARIABLES ALLOCATED BY DAINI IS USED.
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
    integer i,i_1ia,i_2ia,ia,ib,ic,illa,illb,illc,ilma,ilmb,ilmc,ina,inb,inc,inoa,inob,inoc,&
         inva,invb,invc,ioffb,ipoa,ipob,ipoc,ipos,noib,nom
    integer,dimension(0:lno)::ipno,noff
    real(dp) ccia,ccipoa,ccipob
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoa=idapo(ina)
       ipob=idapo(inb)
       ipoc=idapo(inc)
       !         minv = min(inva,invb,invc)
       ccipoa = cc(ipoa)
       ccipob = cc(ipob)
       cc(ipoc) = ccipoa*ccipob
       do i=1,nvmax
          cc(ipoc+i) = ccipoa*cc(ipob+i) + ccipob*cc(ipoa+i)
       enddo
       !         do 30 i=ipoc+minv+1,ipoc+invc
       !  30     cc(i) = zero
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inb,inob,invb,ipob,ilmb,illb)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif

    !     GENERAL CASE
    !     ************
    !
    do i=0,nomax
       noff(i) = idapo(i+2)
       ipno(i) = 0
    enddo
    !
    call daclr_b(1)
    !
    !     RE-SORTING THE VECTOR B INTO PIECES THAT ARE OF ONLY ONE ORDER
    !     *************************************************************
    !
    do ib=ipob,ipob+illb-1
       !
       noib = ieo(ia1(i_1(ib))+ia2(i_2(ib)))
       ipos = ipno(noib) + 1
       ipno(noib) = ipos
       inob = noff(noib) + ipos
       !
       cc(inob) = cc(ib)
       i_1(inob) = i_1(ib)
       i_2(inob) = i_2(ib)
       !
    enddo
    !
    do i=0,nomax
       idall(i+2) = ipno(i)
    enddo
    !
    !     PERFORMING ACTUAL MULTIPLICATION
    !     ********************************
    !
    nom = min(nocut,inoc)
    !
    do ia=ipoa,ipoa+illa-1
       !
       i_1ia = i_1(ia)
       i_2ia = i_2(ia)
       ccia = cc(ia)
       !
       do noib = 0,nom-ieo(ia1(i_1(ia))+ia2(i_2(ia)))
          !
          ioffb = noff(noib)
          !
          do ib = ioffb+1,ioffb+ipno(noib)
             !
             ic = ia2(i_2ia+i_2(ib)) + ia1(i_1ia + i_1(ib))
             ! Georg says maybe needs if(ic/=0)
             if(ic/=0) then
                cc(ic) = cc(ic) + ccia*cc(ib)
             else
                write(6,*) " Georg warn me about ic could be zero"
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
  !$$$$  subroutine dadiv(ina,inb,inc)
  subroutine dadiv_b(ina,inb,inc) !$$$$
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
    real(dp) ck,ck1
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       ipoc = idapo(inc)
       ck=one/cc(ipob)
       ck1=cc(ipoa)*ck
       do i=1,nvmax
          cc(ipoc+i) = (cc(ipoa+i)-cc(ipob+i)*ck1)*ck
       enddo
       cc(ipoc)=ck1
       return
    endif

    idadiv = 0
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    call daall1_b(idadiv,'$$DADIV $$',nomax,nvmax)
    call dafun_b('INV ',inb,idadiv)
    call damul_b(ina,idadiv,inc)
    call dadal1_b(idadiv)
    !
    return
    !$$$$  end subroutine dadiv
  end subroutine dadiv_b !$$$$

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
    real(dp) ccipoa
    !
    !     CASE OF FIRST ORDER ONLY
    !     ************************
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       !         minv = min(inva,invc)
       ccipoa = cc(ipoa)
       cc(ipoc) = ccipoa*ccipoa
       do i=1,nvmax
          cc(ipoc+i) = two*ccipoa*cc(ipoa+i)
       enddo
       return
    endif

    if(ina.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',nomax,nvmax)
       call dasqrt(ina,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
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
    integer i,i_1ia,i_2ia,ia,ib,ib1,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,&
         inva,invc,ioffa,ioffb,ipoa,ipoc,ipos,noia,noib,nom
    integer,dimension(0:lno)::ipno,noff
    real(dp) ccia,ccipoa
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       !         minv = min(inva,invc)
       ccipoa = cc(ipoa)
       cc(ipoc) = ccipoa*ccipoa
       do i=1,nvmax
          cc(ipoc+i) = two*ccipoa*cc(ipoa+i)
       enddo
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !     GENERAL CASE
    !     ************
    !
    do i=0,nomax
       noff(i) = idapo(i+2)
       ipno(i) = 0
    enddo
    !
    call daclr_b(1)
    !
    !     RESORTING THE VECTOR A INTO PIECES THAT ARE OF ONLY ONE ORDER
    !     *************************************************************
    !
    do ia=ipoa,ipoa+illa-1
       !
       noia = ieo(ia1(i_1(ia))+ia2(i_2(ia)))
       ipos = ipno(noia) + 1
       ipno(noia) = ipos
       inoa = noff(noia) + ipos
       !
       cc(inoa) = cc(ia)
       i_1(inoa) = i_1(ia)
       i_2(inoa) = i_2(ia)
       !
    enddo
    !
    do i=0,nomax
       idall(i+2) = ipno(i)
    enddo
    !
    !     PERFORMING ACTUAL MULTIPLICATION
    !     ********************************
    !
    nom = min(nocut,inoc)
    !
    do noia = 0,nom/2
       !
       ioffa = noff(noia)
       !
       do ia=ioffa+1,ioffa+ipno(noia)
          !
          i_1ia = i_1(ia)
          i_2ia = i_2(ia)
          ccia = cc(ia)
          !
          ic = ia2(i_2ia+i_2ia) + ia1(i_1ia+i_1ia)
          cc(ic) = cc(ic) + ccia*ccia
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
                ic = ia2(i_2ia+i_2(ib)) + ia1(i_1ia + i_1(ib))
                cc(ic) = cc(ic) + ccia*cc(ib)
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
  !$$$$  subroutine dacad(ina,ckon,inb)
  subroutine dacad_b(ina,ckon,inb) !$$$$
    !  use da_arrays
    implicit none
    !    integer,dimension(lnv)::jjy
    !    data jjy / lnv*0 /  ! flat zero here

    !     ******************************
    !
    !     THIS SUBROUTINE ADDS THE CONSTANT CKON TO THE VECTOR A
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inb
    integer,parameter,dimension(lnv)::jjx=0
    real(dp) ckon,const
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dacop_b(ina,inb)
    if(nomax.eq.1) then
       cc(idapo(inb)) = cc(idapo(inb)) + ckon
       return
    endif
    !
    call dapek_b(inb,jjx,const)
    call dapok_b(inb,jjx,const+ckon)
    !
    return
    !$$$$  end subroutine dacad
  end subroutine dacad_b !$$$$
  !
  subroutine dacsu_b(ina,ckon,inb) !$$$$
    !$$$$  subroutine dacsu(ina,ckon,inb)
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE SUBTRACTS THE CONSTANT CKON FROM THE VECTOR A
    !
    !-----------------------------------------------------------------------------
    !
    integer ina,inb
    integer,parameter,dimension(lnv)::jjx=0
    real(dp) ckon,const
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dacop_b(ina,inb)
    !
    if(nomax.eq.1) then
       cc(idapo(inb)) = cc(idapo(inb)) - ckon
       return
    endif
    !
    call dapek_b(inb,jjx,const)
    call dapok_b(inb,jjx,const-ckon)
    !
    return
    !$$$$  end subroutine dacsu
  end subroutine dacsu_b !$$$$
  !
  !$$$$  subroutine dasuc(ina,ckon,inb)
  subroutine dasuc_b(ina,ckon,inb) !$$$$
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
    real(dp) ckon
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ipob=idapo(inb)
    ipoa=idapo(ina)
    if(nomax.eq.1) then
       cc(ipob) = ckon - cc(ipoa)
       do i=1,nvmax
          cc(ipob+i) =-cc(ipoa+i)
       enddo
       return
    endif

    call dacsu_b(ina,ckon,inb)
    call dacmu_b(inb,-one,inb)
    !
    return
    !$$$$  end subroutine dasuc
  end subroutine dasuc_b !$$$$
  !
  !$$$$  subroutine dacmu(ina,ckon,inc)
  subroutine dacmu_b(ina,ckon,inc) !$$$$
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
    real(dp) ckon
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = idapo(ina)
       ipoc = idapo(inc)
       do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i) * ckon
       enddo
       return
    endif

    if(ina.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daallno1(incc,'$$DAJUNK$$')
       call dacmut(ina,ckon,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dacmut(ina,ckon,inc)
    endif

    return
    !$$$$  end subroutine dacmu
  end subroutine dacmu_b !$$$$

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
    real(dp) ckon
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       do i=0,nvmax
          cc(ipob+i) = cc(ipoa+i) * ckon
       enddo
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inb,inob,invb,ipob,ilmb,illb)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(abs(ckon).lt.eps) then
       idall(inb) = 0
       return
    endif
    !
    ib = ipob - 1
    !
    do ia=ipoa,ipoa+illa-1
       !
       if(ieo(ia1(i_1(ia))+ia2(i_2(ia))).gt.nocut) goto 100
       ib = ib + 1
       cc(ib) = cc(ia)*ckon
       i_1(ib) = i_1(ia)
       i_2(ib) = i_2(ia)
       !
100    continue
    enddo
    !
    idall(inb) = ib-ipob+1
    if(idall(inb).gt.idalm(inb)) then
       write(line,'(a17)') 'ERROR IN DACMU '
       ipause=mypauses(18,line)
       call dadeb !(31,'ERR DACMU ',1)
    endif
    !
    return
  end subroutine dacmut
  !
  !$$$$  subroutine dacdi(ina,ckon,inb)
  subroutine dacdi_b(ina,ckon,inb) !$$$$
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE DIVIDES THE VECTOR INA BY THE CONSTANT CKON
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ina,inb,ipoa,ipob,ipause,mypauses
    real(dp) ckon
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ckon==zero) then
       if(check_da) then
          C_%STABLE_DA=.false.
          messagelost='constant part zero in dacdi'
          return
       else
          write(line,'(a38)')  'ERROR IN DACDI  CKON IS ZERO'
          ipause=mypauses(25,line)
       endif
    endif

    if(nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       do i=0,nvmax
          cc(ipob+i) = cc(ipoa+i)/ ckon
       enddo
       return
    endif
    !
    call dacmu_b(ina,one/ckon,inb)
    !
    return
    !$$$$  end subroutine dacdi
  end subroutine dacdi_b !$$$$
  !
  !
  !$$$$  subroutine dadic(ina,ckon,inc)
  subroutine dadic_b(ina,ckon,inc) !$$$$
    implicit none
    !     ******************************
    !
    !     THIS SUBROUTINE DIVIDES THE CONSTANT CKON BY THE VECTOR INA
    !
    !-----------------------------------------------------------------------------
    !
    integer i,idadic,ina,inc,ipoa,ipoc
    real(dp) ckon,ck
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    ipoa = idapo(ina)
    if(cc(ipoa)==zero) then
       if(check_da) C_%STABLE_DA=.false.
       messagelost='constant part zero in dadic'
       return
    endif


    if(nomax.eq.1) then
       !         minv = min(inva,invb)
       ipoc = idapo(inc)
       cc(ipoc)=ckon/cc(ipoa)
       ck=cc(ipoc)/cc(ipoa)
       do i=1,nvmax
          cc(ipoc+i) = -cc(ipoa+i)* ck
       enddo
       return
    endif


    !    if(abs(ckon).lt.eps) then    !2002.11.28
    !       call dacon_b(inc,zero)
    !       return
    !    endif

    idadic = 0
    call daall1_b(idadic,'$$DADIC $$',nomax,nvmax)

    !    if(ckon.eq.zero) then
    !       write(line,'(a18)') 'ERROR IN DACDI and DADIC, CKON IS ZERO' !2002.11.28
    !       ipause=mypauses(19,line)
    !       call dadeb !(31,'ERR DACDI ',1)
    !    endif
    !    call dacdi(ina,ckon,idadic)
    !    call dafun_b('INV ',idadic,inc)
    call dafun_b('INV ',ina,idadic)
    call dacmu_b(idadic,ckon,inc)


    call dadal1_b(idadic)
    !
    return
    !$$$$  end subroutine dadic
  end subroutine dadic_b !$$$$
  !
  subroutine dacma(ina,inb,bfac,inc)
    implicit none
    !     **********************************
    !
    !     THIS SUBROUTINE PERFORMS THE OPERATIONS C = A + B*BFAC, WHERE A,B,C ARE
    !     DA VECTORS AND BFAC IS A real(dp). A AND C CAN BE IDENTICAL.
    !     CAN LATER BE REPLACED BY SOMETHING LIKE DAADD WITH MINOR CHANGES.
    !
    !-----------------------------------------------------------------------------
    !
    integer idacma,ina,inb,inc,ipoc,ipob,ipoa,i
    real(dp) bfac
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i)   + cc(ipob+i) * bfac
       enddo
       !         do 8 i=ipoc+minv+1,ipoc+invc
       ! 8       cc(i) = zero
       return
    endif
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    idacma = 0
    call daall1_b(idacma,'$$DACMA $$',nomax,nvmax)
    call dalin_b(ina,+one,inb,bfac,idacma)
    call dacop_b(idacma,inc)
    call dadal1_b(idacma)
    !
    return
  end subroutine dacma
  !
  subroutine dalin_b(ina,afac,inb,bfac,inc) !$$$$
    !$$$$  subroutine dalin(ina,afac,inb,bfac,inc)
    implicit none
    integer ina,inb,inc,incc,ipoc

    real(dp) afac,bfac
    !     ***************************************
    !
    !     THIS SUBROUTINE COMPUTES THE LINEAR COMBINATION
    !     C = AFAC*A + BFAC*B. IT IS ALSO USED TO ADD AND SUBTRACT.
    !
    !-----------------------------------------------------------------------------
    !
    integer  ipob,ipoa,i
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i) * afac + cc(ipob+i) * bfac
       enddo
       !         do 8 i=ipoc+minv+1,ipoc+invc
       ! 8       cc(i) = zero
       return
    endif

    if(ina.eq.inc.or.inb.eq.inc) then
       !        call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',nomax,nvmax)
       call dalint(ina,afac,inb,bfac,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dalint(ina,afac,inb,bfac,inc)
    endif

    return
  end subroutine dalin_b !$$$$
  !$$$$  end subroutine dalin


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
    real(dp) afac,bfac,ccc,copf
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(nomax.eq.1) then
       ipoc = idapo(inc)
       ipoa = idapo(ina)
       ipob = idapo(inb)
       !         minv = min(inva,invb,invc)
       do i=0,nvmax
          cc(ipoc+i) = cc(ipoa+i) * afac + cc(ipob+i) * bfac
       enddo
       !         do 8 i=ipoc+minv+1,ipoc+invc
       ! 8       cc(i) = zero
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inb,inob,invb,ipob,ilmb,illb)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
    ja = ia1(i_1(ia)) + ia2(i_2(ia))
    jb = ia1(i_1(ib)) + ia2(i_2(ib))
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
    ccc = cc(ia)*afac + cc(ib)*bfac
    if(abs(ccc).lt.eps) goto 25
    if(ieo(ia1(i_1(ia))+ia2(i_2(ia))).gt.nocut) goto 25
    ic = ic + 1
    cc(ic) = ccc
    i_1(ic) = i_1(ia)
    i_2(ic) = i_2(ia)
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
    ja = ia1(i_1(ia)) + ia2(i_2(ia))
    jb = ia1(i_1(ib)) + ia2(i_2(ib))
    goto 10
    !
    !     STORING TERM A
    !     **************
    !
30  continue
    if(ieo(ia1(i_1(ia))+ia2(i_2(ia))).gt.nocut) goto 35
    ccc = cc(ia)*afac
    if(abs(ccc).lt.eps) goto 35
    ic = ic + 1
    cc(ic) = ccc
    i_1(ic) = i_1(ia)
    i_2(ic) = i_2(ia)
35  continue
    ia = ia + 1
    if(ia.gt.iamax) then
       ismin = ib
       ismax = ibmax
       copf  = bfac
       goto 50
    endif
    ja = ia1(i_1(ia)) + ia2(i_2(ia))
    goto 10
    !
    !     STORING TERM B
    !     **************
    !
40  continue
    if(ieo(ia1(i_1(ib))+ia2(i_2(ib))).gt.nocut) goto 45
    ccc = cc(ib)*bfac
    if(abs(ccc).lt.eps) goto 45
    ic = ic + 1
    cc(ic) = ccc
    i_1(ic) = i_1(ib)
    i_2(ic) = i_2(ib)
45  continue
    ib = ib + 1
    if(ib.gt.ibmax) then
       ismin = ia
       ismax = iamax
       copf  = afac
       goto 50
    endif
    jb = ia1(i_1(ib)) + ia2(i_2(ib))
    goto 10
    !
    !     COPYING THE REST
    !     ****************
    !
50  continue
    do is=ismin,ismax
       if(ieo(ia1(i_1(is))+ia2(i_2(is))).gt.nocut) goto 60
       ccc = cc(is)*copf
       if(abs(ccc).lt.eps) goto 60
       ic = ic + 1
       cc(ic) = ccc
       i_1(ic) = i_1(is)
       i_2(ic) = i_2(is)
60     continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    !
    if(idall(inc).gt.idalm(inc)) then
       write(line,'(a40)')  'ERROR IN dalint idall(inc).gt.idalm(inc)'
       ipause=mypauses(21,line)
       call dadeb !(31,'ERR DALIN ',1)
    endif
    !
    return
  end subroutine dalint
  !
  subroutine dafun_b(cf,ina,inc) !$$$$
    !$$$$  subroutine dafun(cf,ina,inc)
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       !       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',nomax,nvmax)
       call dafunt(cf,ina,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dafunt(cf,ina,inc)
    endif

    return
  end subroutine dafun_b !$$$$
  !$$$$  end subroutine dafun

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
    integer,parameter,dimension(lnv)::jjx=0
    real(dp) a0,ca,ea,ra,sa
    real(dp),dimension(0:lno)::xf
    character(4) cf,cfh
    character(26) abcs,abcc
    !
    data abcs /'abcdefghijklmnopqrstuvwxyz'/
    data abcc /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
    call daall1_b(ipow,'$$DAFUN1$$',nomax,nvmax)
    call daall1_b(inon,'$$DAFUN2$$',nomax,nvmax)
    call daall1_b(iscr,'$$DAFUN3$$',nomax,nvmax)
    !
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !
    call dapek_b(ina,jjx,a0)
    !
    !      no = min(nocut,inoa,inoc)
    no = min(nocut,nomax)
    !
    !     BRANCHING TO DIFFERENT FUNCTIONS
    !     ********************************
    !
    select case(cf)
    case('INV ')
       !    if(cf.eq.'INV ') then
       !        1/(A0+P) = 1/A0*(1-(P/A0)+(P/A0)**2-...)
       if(a0.eq.0) then
          if(check_da) then
             messagelost="a0.eq.0 for INV in dafun"
             C_%STABLE_DA=.false.
             C_%check_stable=.false.
             call dadal1_b(iscr)
             call dadal1_b(inon)
             call dadal1_b(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       xf(0) = one/a0
       do i=1,no
          xf(i) = -xf(i-1)/a0
       enddo
       !
    case('SQRT')
       !    elseif(cf.eq.'SQRT') then
       !        SQRT(A0+P) = SQRT(A0)*(1+1/2(P/A0)-1/8*(P/A0)**2+...)
       if(a0.le.0) then
          if(check_da) then
             messagelost="a0.le.0 for SQRT in dafun"
             C_%STABLE_DA=.false.
             C_%check_stable=.false.
             call dadal1_b(iscr)
             call dadal1_b(inon)
             call dadal1_b(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       ra = SQRT(a0)
       xf(0) = ra
       do i=1,no
          xf(i) = -xf(i-1)/a0/REAL(2*i,kind=DP)*REAL(2*i-3,kind=DP)
       enddo
       !
       !
    case('EXP ')
       !    elseif(cf.eq.'EXP ') then
       !        EXP(A0+P) = EXP(A0)*(1+P+P**2/2!+...)
       if(a0>hyperbolic_aperture) then
          if(check_da) then
             messagelost="a0>hyperbolic_aperture for EXP in dafun"
             C_%STABLE_DA=.false.
             C_%check_stable=.false.
             call dadal1_b(iscr)
             call dadal1_b(inon)
             call dadal1_b(ipow)
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
       !
    case('LOG ')
       !    elseif(cf.eq.'LOG ') then
       !        LOG(A0+P) = LOG(A0) + (P/A0) - 1/2*(P/A0)**2 + 1/3*(P/A0)**3 - ...)
       if(a0.le.0) then
          if(check_da) then
             messagelost="a0.le.0 for LOG in dafun"
             C_%STABLE_DA=.false.
             C_%check_stable=.false.
             call dadal1_b(iscr)
             call dadal1_b(inon)
             call dadal1_b(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       ea  = LOG(a0)
       xf(0) = ea
       xf(1) = one/a0
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
    case('SINH')
       !    elseif(cf.eq.'SINH') then
       if(a0>hyperbolic_aperture) then
          if(check_da) then
             messagelost="a0>hyperbolic_aperture for SINH in dafun"
             C_%STABLE_DA=.false.
             C_%check_stable=.false.
             call dadal1_b(iscr)
             call dadal1_b(inon)
             call dadal1_b(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       sa  = SINH(a0)
       ca  = COSH(a0)
       xf(0) = sa
       xf(1) = ca
       do i=2,no
          xf(i) = xf(i-2)/REAL(i*(i-1),kind=DP)
       enddo
    case('COSH')
       !    elseif(cf.eq.'COSH') then
       if(a0>hyperbolic_aperture) then
          if(check_da) then
             messagelost="a0>hyperbolic_aperture for COSH in dafun"
             C_%STABLE_DA=.false.
             C_%check_stable=.false.
             call dadal1_b(iscr)
             call dadal1_b(inon)
             call dadal1_b(ipow)
             return
          else
             write(*,1000) cf,ina,a0
             call dadeb !(31,'ERR DAFUN ',1)
             lfun = 0
             return
          endif
       endif
       sa  = SINH(a0)
       ca  = COSH(a0)
       xf(0) = ca
       xf(1) = sa
       xf(0) = ca
       xf(1) = sa
       do i=2,no
          xf(i) = xf(i-2)/REAL(i*(i-1),kind=DP)
       enddo
    case default
       !    else
       write(line,'(a28,1x,a4)')  'ERROR, UNSOPPORTED FUNCTION ',cf
       ipause=mypauses(35,line)
       !    endif
    end select
    !
    call dacon_b(inc,xf(0))
    call dacop_b(ina,inon)
    call dapok_b(inon,jjx,zero)
    call dacon_b(ipow,one)
    !
    do i=1,min(no,nocut)
       !
       call damul_b(inon,ipow,iscr)
       call dacop_b(iscr,ipow)
       call dacma(inc,ipow,xf(i),inc)
       !
    enddo
    !
1000 format('ERROR IN DAFUN, ',a4,' DOES NOT EXIST FOR VECTOR ',i10,'CONST TERM  = ',e12.5)
    !
    call dadal1_b(iscr)
    call dadal1_b(inon)
    call dadal1_b(ipow)
    !
    return
  end subroutine dafunt
  !

  !$$$$  subroutine daabs(ina,anorm)
  subroutine daabs_b(ina,anorm) !$$$$
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
    anorm = zero
    do i=ipoa,ipoa+illa-1
       anorm = anorm + abs(cc(i))
    enddo
    !
    return
    !$$$$  end subroutine daabs
  end subroutine daabs_b !$$$$
  !

  !
  subroutine dacct_b(ma,ia,mb,ib,mc,ic) !$$$$
    !$$$$  subroutine dacct(ma,ia,mb,ib,mc,ic)
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
    integer,dimension(lnv)::monx
    integer,dimension(:)::ma,mb,mc
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(ma(1).eq.mc(1).or.mb(1).eq.mc(1)) then
       call dainf(mc(1),inoc,invc,ipoc,ilmc,illc)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
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
          call dacop_b(monx(i),mc(i))
       enddo
       call dadal_b(monx,ic)
    else
       call dacctt(ma,ia,mb,ib,mc,ic)
    endif

    return
  end subroutine dacct_b !$$$$
  !$$$$  end subroutine dacct

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
    !      INTEGER MON(LNO+1),ICC(LNV)
    !      INTEGER,dimension(:)::MB,MC,MA
    !ETIENNE
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ic,iia,iib,iic,illa,illb,illc,ilma,ilmb,ilmc,inoa,inob,inoc,inva,invb,invc,&
         ipoa,ipob,ipoc,iv,jl,jv,ipause,mypauses
    integer,dimension(lno+1)::mon
    !    integer,dimension(lno)::icc
    integer,dimension(lnv)::icc  !johan 2008 March
    integer,dimension(:)::ma,mb,mc
    real(dp) ccf
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
    do i=1,nomax+1
       mon(i) = 0
    enddo
    !
    call daall(icc,ib,'$$DACCT $$',nomax,nvmax)
    call daall(mon,nomax+1,'$$DAMON $$',inoc,invc)
    !
    call mtree_b(mb,ib,icc,ib)
    !
    !     PERFORMING CONCATENATION
    !     ************************
    !
    do i=1,ia
       call dacon_b(ma(i),cc(idapo(icc(i))))
    enddo
    !
    call dacon_b(mon(1),one)
    !
    do i=1,idall(icc(1))-1
       !
       jl = i_1(idapo(icc(1))+i)
       jv = i_2(idapo(icc(1))+i)
       !
       call damul_b(mon(jl),mc(jv),mon(jl+1))
       !
       do iv=1,ia
          !
          ccf = cc(idapo(icc(iv))+i)
          if(abs(ccf).gt.eps) call dacma(ma(iv),mon(jl+1),ccf,ma(iv))
          !
       enddo
    enddo
    !
    call dadal_b(mon,nomax+1)
    call dadal_b(icc,ib)
    !
    return
  end subroutine dacctt
  !
  !$$$$  subroutine mtree(mb,ib,mc,ic)
  subroutine mtree_b(mb,ib,mc,ic) !$$$$
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
    integer,dimension(lnv)::jjy
    integer,dimension(0:lno)::jv
    integer,dimension(:)::mb,mc
    real(dp) apek,bbijj,chkjj
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
    call daall1_b(ichk,'$$MTREE $$',nomax,nvmax)
    !
    !     FIND ALL THE ENTRIES TO BE LOOKED FOR
    !     *************************************
    !
    call daclr_b(1)
    !
    cc(1) = one
    !
    do i=1,ib
       if(nomax.eq.1) then
          do ib1 =2,invc+1    ! 2,7    ! Etienne
             cc(ib1) = one
          enddo
       else
          do ibi = idapo(mb(i)),idapo(mb(i))+idall(mb(i))-1
             iccx = ia1(i_1(ibi)) + ia2(i_2(ibi))
             if(ieo(iccx).gt.inob) goto 90
             cc(iccx) = one
90           continue
          enddo
       endif
    enddo
    !
    do ii=1,inob
       !
       !     SEARCHING FOR FATHER FOR EACH TERM
       !
       do i=1,nmmax
          if(cc(i).lt.half) goto 140
          !
          jnon = 0
          call dancd(i_1(i),i_2(i),jjy)
          do j=1,invb
             if(jjy(j).eq.0) goto 130
             jnon = j
             jjy(j) = jjy(j) - 1
             call dadcd(jjy,ic1,ic2)
             apek = cc(ia1(ic1)+ia2(ic2))
             jjy(j) = jjy(j) + 1
             if(apek.ge.half) goto 140
130          continue
          enddo
          !
          if(jnon.eq.0) goto 140
          !
          !     TERM IS AN ORPHAN, SO CREATE FOSTER FATHER
          !
          jjy(jnon) = jjy(jnon) - 1
          call dadcd(jjy,ic1,ic2)
          cc(ia1(ic1)+ia2(ic2)) = one
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
    ntermf = idall(ichk)
    !
    !     ZEROTH ORDER TERMS
    !     ******************
    !
    do i=1,lnv
       jjy(i) = 0
    enddo
    !
    do i=1,ib
       call dapek_b(mb(i),jjy,bbijj)
       i_1(idapo(mc(i))) = 0
       i_2(idapo(mc(i))) = 0
       cc(idapo(mc(i))) = bbijj
    enddo
    !
    call dapek_b(ichk,jjy,chkjj)
    if(chkjj.gt.half) then
       call dapok_b(ichk,jjy,-one)
    else
       write(line,'(a49)')  'ERROR IN MTREE, ZEROTH ORDER TERM OF ICHK IS ZERO'
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
    chkjj = one
    !
200 continue
    if(jl.eq.0.and.chkjj.le.half) goto 250
    if(jl.lt.inob.and.chkjj.gt.half) then
       jl = jl + 1
       jjy(1) = jjy(1) + 1
       jv(jl) = 1
    elseif(jv(jl).eq.invb) then
       jjy(jv(jl)) = jjy(jv(jl)) - 1
       jv(jl) = 0
       jl = jl - 1
       chkjj = zero
       goto 200
    else
       jjy(jv(jl)) = jjy(jv(jl)) - 1
       jv(jl) = jv(jl) + 1
       jjy(jv(jl)) = jjy(jv(jl)) + 1
    endif
    !
    call dapek_b(ichk,jjy,chkjj)
    !
    if(chkjj.le.half) goto 200
    !
    nterm = nterm + 1
    if(nterm.gt.idalm(mc(1))) then
       write(line,'(a31)')  'ERROR IN MTREE, NTERM TOO LARGE'
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR MTREE3',1)
    endif
    !
    call dapok_b(ichk,jjy,-one)
    !
    do i=1,ib
       call dapek_b(mb(i),jjy,bbijj)
       i_1(idapo(mc(i))+nterm-1) = jl
       i_2(idapo(mc(i))+nterm-1) = jv(jl)
       cc(idapo(mc(i))+nterm-1) = bbijj
    enddo
    !
    goto 200
    !
250 continue
    !
    do i=1,ib
       idall(mc(i)) = nterm
    enddo
    !
    !     PERFORMING CROSS CHECKS
    !     ***********************
    !
    if(nterm.ne.ntermf.or.nterm.ne.idall(ichk)) then
       write(line,'(a46,1x,i8,1x,i8,1x,i8)') 'ERROR IN MTREE, NTERM, NTERMF, IDALL(ICHK) =  ',nterm,ntermf,idall(ichk)
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR MTREE4',1)
    endif
    !
    do i=idapo(ichk),idapo(ichk)+nterm-1
       if(abs(cc(i)+one).gt.epsmac) then
          write(line,'(a44)')  'ERROR IN MTREE, NOT ALL TERMS IN ICHK ARE -1'
          ipause=mypauses(35,line)
          call dadeb !(31,'ERR MTREE5',1)
       endif
    enddo
    !
    call dadal1_b(ichk)
    !
    return
  end subroutine mtree_b !$$$$
  !$$$$  end subroutine mtree
  !
  subroutine ppushprint(mc,ic,mf,jc,line)
    implicit none
    !
    integer i,ic,iv,jc,jl,jv,mf
    integer,dimension(:)::mc
    character(20) line
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(mf.le.0) return
    write(mf,*) 0,0,jc+1,0,line
    do i=1,ic
       jc=1+jc
       write(mf,*) jc,jl,jv,cc(idapo(mc(i)))
    enddo
    !     xf(i) = cc(idapo(mc(i)))
    !      xm(1) = one
    do i=1,idall(mc(1))-1
       jl = i_1(idapo(mc(1))+i)
       jv = i_2(idapo(mc(1))+i)
       !      xx = xm(jl)*xi(jv)
       !      xm(jl+1) = xx
       do iv=1,ic
          jc=1+jc
          write(mf,*) jc,jl,jv,cc(idapo(mc(iv))+i)
          !        xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
       enddo
    enddo
    return
  end subroutine ppushprint
  !
  !$$$$  subroutine ppushstore(mc,nd2,coef,ml,mv)
  subroutine ppushstore_b(mc,nd2,coef,ml,mv) !$$$$
    implicit none
    !
    integer i,ic,iv,jc,jl,jv,ntot,nd2
    integer,dimension(:), intent(in)::mc
    integer,dimension(:), intent(out)::ml,mv
    real(dp),dimension(:),intent(out)::coef
    !
    jc=0
    jl=0;jv=0;
    ic=nd2
    ntot=idall(mc(1))*ic
    do i=1,ic
       jc=1+jc
       ml(jc)=0
       mv(jc)=0
       coef(jc)=cc(idapo(mc(i)))
       !       write(mf,*) jc,jl,jv,cc(idapo(mc(i)))
    enddo

    do i=1,idall(mc(1))-1
       jl = i_1(idapo(mc(1))+i)
       jv = i_2(idapo(mc(1))+i)
       !      xx = xm(jl)*xi(jv)
       !      xm(jl+1) = xx
       do iv=1,ic
          jc=1+jc
          ml(jc)=jl
          mv(jc)=jv
          coef(jc)=cc(idapo(mc(iv))+i)

          !         write(mf,*) jc,jl,jv,cc(idapo(mc(iv))+i)
          !        xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
       enddo
    enddo
    return
    !$$$$  end subroutine ppushstore
  end subroutine ppushstore_b !$$$$

  !$$$$  subroutine ppushGETN(mc,ND2,ntot)
  subroutine ppushGETN_b(mc,ND2,ntot) !$$$$
    implicit none
    !
    integer ntot,ND2
    integer,dimension(:), intent(inout)::mc
    !
    ntot=idall(mc(1))*nd2
    !$$$$  end subroutine ppushGETN
  end subroutine ppushGETN_b !$$$$
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
    real(dp) xx
    real(dp),dimension(lno+1)::xm
    real(dp),dimension(lno)::xt
    real(dp),dimension(:)::xf,xi
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,ic
       xt(i)=xi(i)
    enddo
    do i=1,ic
       xf(i) = cc(idapo(mc(i)))
    enddo
    !
    xm(1) = one
    !
    do i=1,idall(mc(1))-1
       !
       jl = i_1(idapo(mc(1))+i)
       jv = i_2(idapo(mc(1))+i)
       xx = xm(jl)*xt(jv)
       xm(jl+1) = xx
       !
       do iv=1,ic
          xf(iv) = xf(iv) + cc(idapo(mc(iv))+i) * xx
       enddo
    enddo
    do i=1,nvmax
       if(abs(xf(i))>da_absolute_aperture.and.check_da) then
          C_%STABLE_DA=.false.
          write(6,*) "unstable in ppush ",i,xf(i)
       endif
    enddo
    !
    return
  end subroutine ppush

  !$$$$  subroutine ppush1(mc,xi,xf)
  subroutine ppush1_b(mc,xi,xf) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE MATRIX WHOSE TREE IS STORED IN CA VECTOR MC
    !     TO THE COORDINATES IN XI AND STORES THE RESULT IN XF
    !
    !-----------------------------------------------------------------------------
    !
    integer i,jl,jv,mc
    real(dp) xf,xx
    real(dp),dimension(:)::xi
    real(dp),dimension(lno)::xt
    real(dp),dimension(lno+1)::xm
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,nvmax
       xt(i)=xi(i)
    enddo

    xf = cc(idapo(mc))
    !
    xm(1) = one
    !
    do i=1,idall(mc)-1
       !
       jl = i_1(idapo(mc)+i)
       jv = i_2(idapo(mc)+i)
       xx = xm(jl)*xt(jv)
       xm(jl+1) = xx
       !
       xf = xf + cc(idapo(mc)+i) * xx
    enddo


    if(abs(xf)>da_absolute_aperture.and.check_da) then
       C_%STABLE_DA=.false.
       write(6,*) "unstable in ppush1 ", xf
       write(6,*) xi(1:nvmax)
    endif
    return
    !$$$$  end subroutine ppush1
  end subroutine ppush1_b !$$$$

  !$$$$  subroutine dainv(ma,ia,mb,ib)
  subroutine dainv_b(ma,ia,mb,ib) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
    !     STORES THE RESULT IN MI
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob
    integer,dimension(lnv)::jj,ml
    integer,dimension(:)::ma,mb
    real(dp),dimension(lnv)::x
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=1,lnv
       jj(i)=0
    enddo

    if(ma(1).eq.mb(1)) then
       call dainf(mb(1),inob,invb,ipob,ilmb,illb)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       do i=1,ia
          call dapok_b(ma(i),jj,zero)
       enddo
       do ij=1,ib
          ml(ij)=0
       enddo
       call daall(ml,ib,'$$DAJUNK$$',inob,invb)
       call dainvt(ma,ia,ml,ib)
       do i=1,ib
          call dacop_b(ml(i),mb(i))
       enddo
       call dadal_b(ml,ib)
    else
       do i=1,ia
          call dapek_b(ma(i),jj,x(i))
          call dapok_b(ma(i),jj,zero)
       enddo
       call dainvt(ma,ia,mb,ib)
       do i=1,ia
          call dapok_b(ma(i),jj,x(i))
       enddo
    endif

    return
    !$$$$  end subroutine dainv
  end subroutine dainv_b !$$$$

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
         ipoa,ipob,j,k,nocut0,ipause,mypauses
    integer,dimension(lnv)::jj,ms,ml
    integer,dimension(:)::ma,mb
    real(dp),dimension(lnv,lnv)::aa,ai
    real(dp) amjj,amsjj,prod
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
    call dainf(mb(1),inob,invb,ipob,ilmb,illb)
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
    call damch(ma,ia)
    call damch(mb,ib)
    !etienne
    do ie=1,ib
       call dacon_b(mb(ie),zero)
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
          call dapek_b(ma(i),jj,amjj)
          if(abs(amjj).gt.eps) call dapok_b(ma(i),jj,zero)
          aa(i,j) = amjj
       enddo
       call dacmu_b(ma(i),-one,ma(i))
    enddo
    !
    !     INVERTING LINEAR MATRIX, CHECKING RESULT AND STORING IN ML
    !     **********************************************************
    !
    call matinv(aa,ai,ia,lnv,ier)
    !
    if(ier.eq.132) then
       if(check_da) then
          C_%STABLE_DA=.false.
          messagelost='ERROR IN ROUTINE DAINV, ier=132 in matinv'
          call dadal_b(ml,ia)
          call dadal_b(ms,ia)
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
          prod = zero
          do k=1,ib
             jj(k) = 0
             prod = prod + aa(i,k)*ai(k,j)
          enddo
          if(i.eq.j) prod = prod - one
          if(abs(prod).gt.c_100*epsmac) then
             write(6,*) " abs(prod) > c_100*epsmac in dainvt",abs(prod), c_100*epsmac
             if(check_da) then
                C_%STABLE_DA=.false.
                messagelost='ERROR IN ROUTINE DAINV, abs(prod).gt.c_100*epsmac '
                call dadal_b(ml,ia)
                call dadal_b(ms,ia)
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
          call dapok_b(mb(i),jj,ai(i,j))
          call dapok_b(ml(i),jj,ai(i,j))
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
    nocut0 = nocut
    !
    do i=2,nocut
       !
       nocut = i
       !
       call dacct_b(ma,ia,mb,ib,ms,ia)
       do j=1,ib
          do k=1,ib
             jj(k) = 0
          enddo
          jj(j) = 1
          call dapek_b(ms(j),jj,amsjj)
          call dapok_b(ms(j),jj,amsjj+one)
       enddo
       !
       call dacct_b(ml,ia,ms,ia,mb,ib)
       !
    enddo
    !
    nocut = nocut0
    !
    !     FLIPPING BACK SIGN OF A, FILLING UP FIRST ORDER PART AGAIN
    !     **********************************************************
    !
    do i=1,ib
       call dacmu_b(ma(i),-one,ma(i))
       do j=1,ib
          do k=1,ib
             jj(k) = 0
          enddo
          jj(j) = 1
          call dapok_b(ma(i),jj,aa(i,j))
       enddo
    enddo
    !
    call dadal_b(ml,ia)
    call dadal_b(ms,ia)
    !
    return
  end subroutine dainvt
  !
  !$$$$  subroutine dapin(ma,ia,mb,ib,jx)
  subroutine dapin_b(ma,ia,mb,ib,jx) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE INVERTS THE MATRIX MA WITH IA DA VECTORS AND
    !     STORES THE RESULT IN MI
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ib,ij,illb,ilmb,inob,invb,ipob
    integer,dimension(lnv)::jj,ml
    integer,dimension(:)::ma,mb,jx
    real(dp),dimension(lnv)::x
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    do i=1,lnv
       jj(i)=0
    enddo

    if(ma(1).eq.mb(1)) then
       call dainf(mb(1),inob,invb,ipob,ilmb,illb)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       do i=1,ia
          call dapok_b(ma(i),jj,zero)
       enddo
       do ij=1,ib
          ml(ij)=0
       enddo
       call daall(ml,ib,'$$DAJUNK$$',inob,invb)
       call dapint(ma,ia,ml,ib,jx)
       do i=1,ib
          call dacop_b(ml(i),mb(i))
       enddo
       call dadal_b(ml,ib)
    else
       do i=1,ia
          call dapek_b(ma(i),jj,x(i))
          call dapok_b(ma(i),jj,zero)
       enddo
       call dapint(ma,ia,mb,ib,jx)
       do i=1,ia
          call dapok_b(ma(i),jj,x(i))
       enddo
    endif

    return
    !$$$$  end subroutine dapin
  end subroutine dapin_b !$$$$

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
    integer,dimension(lnv)::jj,mn,mi,me
    integer,dimension(:)::ma,mb,jind
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ma(1),inoa,inva,ipoa,ilma,illa)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
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
       do k=1,nvmax
          jj(k) = 0
       enddo
       jj(i) = 1
       call dapok_b(me(i),jj,one)
    enddo
    !
    do i=1,ia
       call dacop_b(ma(i),mn(i))
       if(jind(i).eq.0) call dacop_b(me(i),mn(i))
    enddo
    !
    call dainv_b(mn,ia,mi,ia)
    !
    do i=1,ia
       if(jind(i).eq.0) call dacop_b(ma(i),me(i))
    enddo
    !
    call dacct_b(me,ia,mi,ia,mb,ib)
    !
    call dadal_b(me,ia)
    call dadal_b(mi,ia)
    call dadal_b(mn,ia)
    !
    return
  end subroutine dapint
  !
  !$$$$  subroutine dader(idif,ina,inc)
  subroutine dader_b(idif,ina,inc) !$$$$
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',inoc,invc)
       call dadert(idif,ina,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dadert(idif,ina,inc)
    endif

    return
  end subroutine dader_b !$$$$
  !$$$$  end subroutine dader

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
    integer,dimension(lnv)::jd
    real(dp) rr,x,xdivi
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(nomax.eq.1) then
       !         PRINT*,'ERROR, DADER CALLED WITH NOMAX = 1'
       !        call dadeb !(31,'ERR DADER1',1)
       !         stop 666
       do i=1,lnv
          jd(i)=0
       enddo
       jd(idif)=1
       call dapek_b(ina,jd,rr)
       call dacon_b(inc,rr)
       return
    endif
    !
    !      CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !
    ibase = nomax + 1
    !
    if(idif.gt.(nvmax+1)/2) then
       ider1  = 0
       ider1s = 0
       ider2  = idif-(nvmax+1)/2
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
    ibase = nomax+1
    !
    ic = ipoc-1
    !
    do i=ipoa,ipoa+illa-1
       !
       if(ider1.eq.0) then
          iee = i_2(i)
       else
          iee = i_1(i)
       endif
       !
       x = iee/xdivi
       ifac = int(ibase*(x-int(x+epsmac)+epsmac))
       !
       if(ifac.eq.0) goto 100
       !
       ic = ic + 1
       cc(ic) = cc(i)*ifac
       i_1(ic) = i_1(i) - ider1s
       i_2(ic) = i_2(i) - ider2s
       !
100    continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    if(idall(inc).gt.idalm(inc)) then
       write(line,'(a15)') 'ERROR IN DADER '
       ipause=mypauses(35,line)
       call dadeb !(31,'ERR DADER2',1)
    endif
    !
    return
  end subroutine dadert
  !

  !
  !$$$$  subroutine dacfuR(ina,fun,inc)
  subroutine dacfuR_b(ina,fun,inc) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc
    complex(dp),external::fun
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',inoc,invc)
       call dacfuRt(ina,fun,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dacfuRt(ina,fun,inc)
    endif

    return
    !$$$$  end subroutine dacfuR
  end subroutine dacfuR_b !$$$$

  subroutine dacfuRt(ina,fun,inc)
    implicit none
    ! external fun
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,&
         invc,ipoa,ipoc,ipause,mypauses
    integer,dimension(lnv)::j
    real(dp) cfac,rr
    !
    interface
       !       complex(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(nomax.eq.1) then
       do i=1,lnv
          j(i)=0
       enddo
       call dapek_b(ina,j,rr)
       cfac = REAL(fun(j),kind=DP)
       rr=cfac*rr
       call dapok_b(inc,j,rr)
       do i=1,lnv
          j(i)=1
          call dapek_b(ina,j,rr)
          cfac = REAL(fun(j),kind=DP)
          rr=cfac*rr
          call dapok_b(inc,j,rr)
          j(i)=0
       enddo
       !         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
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
       call dancd(i_1(ia),i_2(ia),j)
       cfac = REAL(fun(j),kind=DP)
       !      IF(abs(CFAC).LT.EPS) GOTO 100
       !      IF(abs(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100
       !
       ic = ic + 1
       cc(ic) = cc(ia)*cfac
       i_1(ic) = i_1(ia)
       i_2(ic) = i_2(ia)
       !
100    continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    if(idall(inc).gt.idalm(inc)) then
       write(line,'(a15)') 'ERROR IN DACFU '
       ipause=mypauses(36,line)
       call dadeb !(31,'ERR DACFU ',1)
    endif
    !
    return
  end subroutine dacfuRt
  !
  !$$$$  subroutine dacfu(ina,fun,inc)
  subroutine dacfu_b(ina,fun,inc) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL real(dp) FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc
    real(dp),external::fun
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',inoc,invc)
       call dacfut(ina,fun,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dacfut(ina,fun,inc)
    endif

    return
    !$$$$  end subroutine dacfu
  end subroutine dacfu_b !$$$$

  !$$$$  subroutine dacfuI(ina,fun,inc)
  subroutine dacfuI_b(ina,fun,inc) !$$$$
    implicit none
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer illc,ilmc,ina,inc,incc,inoc,invc,ipoc
    complex(dp),external::fun
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.eq.inc) then
       call dainf(inc,inoc,invc,ipoc,ilmc,illc)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
             write(6,*) "big problem in dabnew ", sqrt(crash)
          endif
          return
       endif
       incc=0
       call daall1_b(incc,'$$DAJUNK$$',inoc,invc)
       call dacfuIt(ina,fun,incc)
       call dacop_b(incc,inc)
       call dadal1_b(incc)
    else
       call dacfuIt(ina,fun,inc)
    endif

    return
    !$$$$  end subroutine dacfuI
  end subroutine dacfuI_b !$$$$

  subroutine dacfuIt(ina,fun,inc)
    implicit none
    !  external fun
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL DOUBLE COMPLEX FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,&
         invc,ipoa,ipoc,ipause,mypauses
    integer,dimension(lnv)::j
    real(dp) cfac,rr
    !
    interface
       !       complex(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         complex(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(nomax.eq.1) then
       do i=1,lnv
          j(i)=0
       enddo
       call dapek_b(ina,j,rr)
       cfac = AIMAG(fun(j))
       rr=cfac*rr
       call dapok_b(inc,j,rr)
       do i=1,lnv
          j(i)=1
          call dapek_b(ina,j,rr)
          cfac = AIMAG(fun(j))
          rr=cfac*rr
          call dapok_b(inc,j,rr)
          j(i)=0
       enddo
       !         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
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
       call dancd(i_1(ia),i_2(ia),j)
       cfac = aimag(fun(j))
       !      IF(abs(CFAC).LT.EPS) GOTO 100
       !      IF(abs(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100
       !
       ic = ic + 1
       cc(ic) = cc(ia)*cfac
       i_1(ic) = i_1(ia)
       i_2(ic) = i_2(ia)
       !
100    continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    if(idall(inc).gt.idalm(inc)) then
       write(line,'(a15)') 'ERROR IN DACFU '
       ipause=mypauses(37,line)
       call dadeb !(31,'ERR DACFU ',1)
    endif
    !
    return
  end subroutine dacfuIt
  !

  subroutine dacfut(ina,fun,inc)
    implicit none
    ! external fun
    !     *****************************
    !
    !     THIS SUBROUTINE APPLIES THE EXTERNAL real(dp) FUNCTION
    !     OF THE EXPONENTS FUN TO EACH COEFFICIENT OF A AND STORES THE
    !     RESULT IN C
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ia,ic,illa,illc,ilma,ilmc,ina,inc,inoa,inoc,inva,&
         invc,ipoa,ipoc,ipause,mypauses
    integer,dimension(lnv)::j
    real(dp) cfac,rr
    !
    interface
       !       real(kind(one)) function fun(abc)
       function fun(abc)
         use precision_constants
         implicit none
         real(dp) fun
         integer,dimension(:)::abc
       end function fun
    end interface
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    if(nomax.eq.1) then
       do i=1,lnv
          j(i)=0
       enddo
       call dapek_b(ina,j,rr)
       cfac = fun(j)
       rr=cfac*rr
       call dapok_b(inc,j,rr)
       do i=1,lnv
          j(i)=1
          call dapek_b(ina,j,rr)
          cfac = fun(j)
          rr=cfac*rr
          call dapok_b(inc,j,rr)
          j(i)=0
       enddo
       !         PRINT*,'ERROR, DACFU CALLED WITH NOMAX = 1'
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
       call dancd(i_1(ia),i_2(ia),j)
       cfac = fun(j)
       !      IF(abs(CFAC).LT.EPS) GOTO 100
       !      IF(abs(CFAC*CC(IA)).LT.EPS) GOTO 100
       if(abs(cfac*cc(ia)).lt.eps.or.abs(cc(ia)).lt.eps) goto 100
       !
       ic = ic + 1
       cc(ic) = cc(ia)*cfac
       i_1(ic) = i_1(ia)
       i_2(ic) = i_2(ia)
       !
100    continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    if(idall(inc).gt.idalm(inc)) then
       write(line,'(a15)') 'ERROR IN DACFU '
       ipause=mypauses(38,line)
       call dadeb !(31,'ERR DACFU ',1)
    endif
    !
    return
  end subroutine dacfut

  !  subroutine GET_C_J_b(ina,I,C,J)
!!$$$$  subroutine GET_C_J(ina,I,C,J)
  !    implicit none
  !
  !    INTEGER I,ina
  !    integer, dimension(lnv)::j
  !    real(dp) C
  !
  !    C=CC(I)
  !    call dancd(i_1(I),i_2(I),J)!

  !  END subroutine GET_C_J_b
!!$$$$  END subroutine GET_C_J

  subroutine dapri_b(ina,iunit) !$$$$
    !$$$$  subroutine dapri(ina,iunit)
    implicit none
    !     ***************************
    !       Frank
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,iii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit,k
    integer,dimension(lnv)::j
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.lt.1.or.ina.gt.nda_dab) then
       print*,'ERROR IN DAPRI, INA = ',ina
       C_%STABLE_DA=.false.
    endif
    !
    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    !
    write(iunit,'(/1X,A,A,I5,A,I5,A,I5/1X,A/)') daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,&
         '*********************************************'
    !
    iout = 0
    ioa = 0

    if(inva.eq.0) then
       write(iunit,'(A)') '    I  VALUE  '
       do i = ipoa,ipoa+illa-1
          write(iunit,'(I6,2X,g20.13)') i-ipoa, cc(i)
       enddo
    elseif(nomax.eq.1) then
       if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
       if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
       do i=1,illa
          do k=1,inva
             j(k)=0
          enddo
          iout=iout+1
          if(i.ne.1) then
             j(i-1)=1
             ioa=1
          endif
          write(iunit,'(I6,2X,g20.13,I5,4X,18(2i2,1X))') iout,cc(ipoa+i-1),ioa,(j(iii),iii=1,nvmax)
          write(iunit,*) cc(ipoa+i-1)
       enddo
    else
       if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
       if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
       do ioa = 0,inoa
          do ii=ipoa,ipoa+illa-1
             if(ieo(ia1(i_1(ii))+ia2(i_2(ii))).ne.ioa) goto 100
             call dancd(i_1(ii),i_2(ii),j)
             !ETIENNE
             if(abs(cc(ii)).gt.eps) then
                !ETIENNE
                iout = iout+1
                write(iunit,'(I6,2X,g20.13,I5,4X,18(2i2,1X))') iout,cc(ii),ioa,(j(iii),iii=1,nvmax)
                !ETIENNE
                write(iunit,*) cc(ii)
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
    !$$$$  end subroutine dapri
  end subroutine dapri_b !$$$$

  !$$$$  subroutine dapri77(ina,iunit)
  subroutine dapri77_b(ina,iunit) !$$$$
    implicit none
    !     ***************************
    !       Etienne
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,iunit
    integer,dimension(lnv)::j
    character(10) c10,k10
    !
    if(iunit.eq.0) return

    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.lt.1.or.ina.gt.nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN DAPRI, INA = ',ina
       C_%STABLE_DA=.false.
    endif
    !
    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    !
    write(iunit,'(/1X,A10,A6,I5,A6,I5,A7,I5/1X,A/)') daname(ina),', NO =',inoa,', NV =',inva,', INA =',ina,&
         '*********************************************'
    !
    if(illa.ne.0) write(iunit,'(A)') '    I  COEFFICIENT          ORDER   EXPONENTS'
    if(illa.eq.0) write(iunit,'(A)') '   ALL COMPONENTS ZERO '
    !
    c10='      NO ='
    k10='      NV ='

    write(iunit,'(A10,I6,A10,I6)') c10,inoa,k10,inva

    iout = 0
    !
    !      DO 100 IOA = 0,INOA
    do ioa = 0,nocut
       do ii=ipoa,ipoa+illa-1
          if(nomax.ne.1) then
             if(ieo(ia1(i_1(ii))+ia2(i_2(ii))).ne.ioa) goto 100
          endif
          !ETIENNE
          if(abs(cc(ii)).gt.eps) then
             !ETIENNE

             if(nomax.ne.1) then
                call dancd(i_1(ii),i_2(ii),j)
                iout = iout+1
             else
                if(ii.eq.ipoa.and.ioa.eq.1) goto 100
                if(ii.gt.ipoa.and.ioa.eq.0) goto 100
                do i=1,lnv
                   j(i)=0
                enddo
                if(ii.ne.ipoa) j(ii-ipoa)=1
                iout = iout+1
             endif
             !

             !      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
             if(abs(cc(ii)).gt.eps) then
                if(eps.gt.c_1d_37) then
                   write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
                else
                   write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
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
    do i=1,lnv
       j(i)=0
    enddo

    if(iout.eq.0) iout=1

    write(iunit,502) -iout,zero,(j(i),i=1,inva)
    !
    return
    !$$$$  end subroutine dapri77
  end subroutine dapri77_b !$$$$

  !$$$$  subroutine dashift(ina,inc,ishift)
  subroutine dashift_b(ina,inc,ishift) !$$$$
    implicit none
    !      real(dp) c
    !       Frank
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ii,illa,ilma,ina,inoa,inva,ioa,iout,ipoa,inb,ishift,ich,&
         ik,inc,k
    integer,dimension(lnv)::j,jd
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    inb=0
    if(ina.lt.1.or.ina.gt.nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN dashift, INA = ',ina
       C_%STABLE_DA=.false.
    endif
    !
    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    call daall1_b(inb,'$$DAJUNK$$',inoa,inva)

    iout = 0
    !
    !      DO 100 IOA = 0,INOA
    do ioa = 0,nocut
       do ii=ipoa,ipoa+illa-1
          if(nomax.ne.1) then
             if(ieo(ia1(i_1(ii))+ia2(i_2(ii))).ne.ioa) goto 100
          endif
          !ETIENNE
          if(abs(cc(ii)).gt.eps) then
             !ETIENNE



             if(nomax.ne.1) then
                call dancd(i_1(ii),i_2(ii),j)

                iout = iout+1
             else
                if(ii.eq.ipoa.and.ioa.eq.1) goto 100
                if(ii.gt.ipoa.and.ioa.eq.0) goto 100
                do i=1,lnv
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

             !      WRITE(IUNIT,*) IOA,CC(II),(J(I),I=1,INVA)
             if(abs(cc(ii)).gt.eps) then
                if(eps.gt.c_1d_37) then
                   !       write(iunit,501) ioa,cc(ii),(j(i),i=1,inva)
                   ich=1
                   do ik=1,ishift
                      if(j(ik).ne.0) ich=0
                   enddo
                   if(ich.eq.1) then
                      do ik=1,lnv
                         jd(ik)=0
                      enddo
                      do ik=ishift+1,lnv
                         jd(ik-ishift)=j(ik)  !%%%%%%Etienne
                      enddo
                   endif

                   call dapok_b(inb,jd,cc(ii))
                else
                   !       write(iunit,503) ioa,cc(ii),(j(i),i=1,inva)
                   ich=1
                   do ik=1,ishift
                      if(j(ik).ne.0) ich=0
                   enddo
                   if(ich.eq.1) then
                      do ik=1,lnv
                         jd(ik)=0
                      enddo
                      do ik=ishift+1,lnv
                         jd(ik-ishift)=j(ik)  !%%%%%%Etienne
                      enddo
                   endif
                   call dapok_b(inb,jd,cc(ii))
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
    do i=1,lnv
       j(i)=0
    enddo

    call dacop_b(inb,inc)
    call dadal1_b(inb)
    !
    return
    !$$$$  end subroutine dashift
  end subroutine dashift_b !$$$$

  !$$$$  subroutine darea(ina,iunit)
  subroutine darea_b(ina,iunit) !$$$$
    implicit none
    !       Frank
    !-----------------------------------------------------------------------------
    !
    integer i,ic,iche,ii,ii_1,ii_2,iin,illa,ilma,ina,inoa,inva,io,io1,ipoa,iunit,&
         iwarin,iwarno,iwarnv,nno
    integer,dimension(lnv)::j
    real(dp) c
    character(10) c10
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    if(ina.lt.1.or.ina.gt.nda_dab) then
       C_%STABLE_DA=.false.
       print*,'ERROR IN DAREA, INA = ',ina
       !        X = SQRT(-ONE)
       !        PRINT*,X
    endif
    !
    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    !
    do i=1,lnv
       j(i) = 0
    enddo
    !
    call daclr_b(1)
    call daclr_b(ina)   ! etienne 2008
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
    read(iunit,'(I6,2X,g20.13,I5,4X,18(2i2,1X))') ii,c,io,(j(i),i=1,inva)
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
    if(nomax.ne.1) then
       ic = ic + 1
       call dadcd(j,ii_1,ii_2)
       ic = ia1(ii_1) + ia2(ii_2)
       cc(ic) = c
       goto 10
    else
       iche=0
       do i=1,inva
          if(j(i).eq.1) iche=i
       enddo
       cc(ipoa+iche)=c
       goto 10
    endif

    !
20  continue
    !
    if(nomax.ne.1) call dapac(ina)
    !
    return
    !$$$$  end subroutine darea
  end subroutine darea_b !$$$$
  !FF
  !
  !$$$$  subroutine darea77(ina,iunit)
  subroutine darea77_b(ina,iunit) !$$$$
    implicit none
    !     ***************************
    !     Etienne
    !     THIS SUBROUTINE READS THE DA VECTOR INA FROM UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ic,iche,ii,ii_1,ii_2,iin,illa,ilma,ina,inoa,inva,ipoa,iunit,&
         k,nojoh,nvjoh
    integer,dimension(lnv)::j
    real(dp) c
    character(10) c10,k10
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    if(ina.lt.1.or.ina.gt.nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN darea77, INA = ',ina
       C_%STABLE_DA=.false.
    endif
    !
    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    !
    do i=1,lnv
       j(i) = 0
    enddo
    !
    call daclr_b(1)
    call daclr_b(ina)   ! etienne 2008
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
    read(iunit,*) ii,c,(j(k),k=1,nvjoh)
    if(ii.lt.0) goto 20

    do i=inva+1,nvjoh
       if(j(i).ne.0) goto 10
    enddo
    iche=0
    do i=1,inva
       iche=iche+j(i)
    enddo
    if(iche.gt.nomax) goto 10
    if(nomax.ne.1) then
       call dadcd(j,ii_1,ii_2)
       ic = ia1(ii_1) + ia2(ii_2)
       cc(ic) = c
    else
       iche=0
       do i=1,inva
          if(j(i).eq.1) iche=i
       enddo
       cc(ipoa+iche)=c
    endif
    goto 10
    !
20  continue
    !
    if(nomax.ne.1) call dapac(ina)
    !
    return
    !$$$$  end subroutine darea77
  end subroutine darea77_b !$$$$

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
    C_%STABLE_DA=.false.
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
    if(inc.ge.1.and.inc.le.nda_dab) then
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
    real(dp) ccc
    !
    !      call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    ipoc = idapo(inc)

    !
    ic = ipoc - 1
    !
    do i=1,nmmax
       ccc = cc(i)
       if(abs(ccc).lt.eps) goto 100
       ic = ic + 1
       cc(ic) = ccc
       i_1(ic) = ie1(i)
       i_2(ic) = ie2(i)
100    continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    if(idall(inc).gt.idalm(inc)) then
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
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(iaa(1),ino1,inv1,ipoa,ilma,illa)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    do i=2,ia
       call dainf(iaa(i),inoi,invi,ipoa,ilma,illa)
       if((.not.C_%STABLE_DA)) then
          if(c_%watch_user) then
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
    !     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES i_1,i_2.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic1,ic2,isplit
    integer,dimension(lnv)::jj
    !
    ibase = nomax + 1
    isplit = (nvmax+1)/2
    ic1 = 0
    ic2 = 0
    !
    do i=nvmax,isplit+1,-1
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
    !     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES i_1,i_2.
    !
    !-----------------------------------------------------------------------------
    !
    integer i,ibase,ic,ic1,ic2,isplit
    integer,dimension(:)::jj
    real(dp) x
    !
    ibase = nomax + 1
    isplit = (nvmax+1)/2
    !
    ic = ic1
    do i=1,isplit
       x  = REAL(ic,kind=DP)/REAL(ibase,kind=DP)
       ic = int(x+epsmac)
       jj(i) = nint(ibase*(x-ic))
    enddo
    !
    ic = ic2
    do i=isplit+1,nvmax
       x  = REAL(ic,kind=DP)/REAL(ibase,kind=DP)
       ic = int(x+epsmac)
       jj(i) = nint(ibase*(x-ic))
    enddo
    !
    do i=nvmax+1,size(jj)    ! 2002.12.2
       jj(i) = 0
    enddo
    !
    return
  end subroutine dancd

  !ETIENNE
  !$$$$  subroutine datra(idif,ina,inc)
  subroutine datra_b(idif,ina,inc) !$$$$
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
    real(dp) x,xdivi
    !
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    call dainf(ina,inoa,inva,ipoa,ilma,illa)
    call dainf(inc,inoc,invc,ipoc,ilmc,illc)
    if((.not.C_%STABLE_DA)) then
       if(c_%watch_user) then
          write(6,*) "big problem in dabnew ", sqrt(crash)
       endif
       return
    endif
    !
    !
    !       CALL DACHK(INA,INOA,INVA, '          ',-1,-1,INC,INOC,INVC)
    !

    if(nomax.eq.1) then
       call dader_b(idif,ina,inc)
       return
    endif
    ibase = nomax + 1
    !
    if(idif.gt.(nvmax+1)/2) then
       ider1  = 0
       ider1s = 0
       ider2  = idif-(nvmax+1)/2
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
    ibase = nomax+1
    !
    ic = ipoc-1
    !
    do i=ipoa,ipoa+illa-1
       !
       if(ider1.eq.0) then
          iee = i_2(i)
       else
          iee = i_1(i)
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
       cc(ic) = cc(i)
       i_1(ic) = i_1(i) - ider1s
       i_2(ic) = i_2(i) - ider2s
       !
100    continue
    enddo
    !
    idall(inc) = ic - ipoc + 1
    if(idall(inc).gt.idalm(inc)) then
       write(line,'(a16)')  'ERROR IN DADTRA '
       ipause=mypauses(35,line)
       call dadeb !(111,'ERR DADTRA',1)
    endif
    !
    return
    !$$$$  end subroutine datra
  end subroutine datra_b !$$$$


  subroutine hash(no1,nv1,jj,ic1,ic2)
    implicit none
    !     ****************************
    !
    !     THIS SUBROUTINE CODES THE EXPONENTS IN JJ INTO THEIR DA CODES i_1,i_2.
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
    !     THIS SUBROUTINE ENCODES THE EXPONENTS IN JJ FROM THEIR DA CODES i_1,i_2.
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


  !$$$$  subroutine daran(ina,cm,xran)
  subroutine daran_b(ina,cm,xran) !$$$$
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
    if(inva.eq.0.or.nomax.eq.1) then
       do i=ipoa,ipoa+ilma-1
          if(cm.gt.zero) then
             cc(i) = bran(xran)
             if(cc(i).gt.cm) cc(i) = zero
          elseif(cm.lt.zero) then
             cc(i) = int(1+10*bran(xran))
             if(cc(i).gt.-ten*cm) cc(i) = zero
          endif
       enddo
       idall(ina) = idalm(ina)
       return
    endif
    !
    if(inoa.ne.nomax.or.inva.ne.nvmax) then
       line='ERROR IN DARAN, ONLY VECTORS WITH NO = NOMAX AND NV = NVMAX ALLOWED'
       ipause=mypauses(31,line)
       call dadeb !(31,'ERR DARAN1',1)
    endif
    !
    call daclr_b(1)
    !
    do i=1,nmmax
       if(cm.gt.zero) then
          cc(i) = bran(xran)
          if(cc(i).gt.cm) cc(i) = zero
       elseif(cm.lt.zero) then
          cc(i) = int(1+10*bran(xran))
          if(cc(i).gt.-ten*cm) cc(i) = zero
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
    !$$$$  end subroutine daran
  end subroutine daran_b !$$$$
  !


  !$$$$  subroutine dacycle(ina,ipresent,value,illa,j)
  subroutine dacycle_b(ina,ipresent,value,illa,j) !$$$$
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
    real(dp) value
    !
    if(ina.lt.1.or.ina.gt.nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN dacycle, INA = ',ina
       ipause=mypauses(39,line)
       stop
    endif
    !

    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    if(.not.present(j)) return
    iout = 0
    if(ipresent.gt.illa.or.ipresent<1) then
       write(6,*) ipresent,illa
       write(6,*) " error in dacycle "
       stop 101
    endif
    ii=ipresent+ipoa-1
    call dancd(i_1(ii),i_2(ii),j)
    value=cc(ii)

    if(nomax==1) then
       j=0
       if(ipresent/=1) j(ipresent-1)=1
    endif
    !    ipresent=1+ipresent



    return

    !$$$$  end subroutine dacycle
  end subroutine dacycle_b !$$$$

!!!! new stuff lingyun

  !$$$$  subroutine daclean(ina,value)
  subroutine daclean_b(ina,value) !$$$$
    implicit none
    integer ipause, mypauses
    !     ***************************
    !
    !     THIS SUBROUTINE PRINTS THE DA VECTOR INA TO UNIT IUNIT.
    !
    !-----------------------------------------------------------------------------
    !
    integer ii,illa,ilma,ina,inoa,inva,iout,ipoa
    real(dp) value
    !
    if(ina.lt.1.or.ina.gt.nda_dab) then
       write(line,'(a22,i8)') 'ERROR IN dacycle, INA = ',ina
       ipause=mypauses(39,line)
       stop
    endif
    !

    inoa = idano(ina)
    inva = idanv(ina)
    ipoa = idapo(ina)
    ilma = idalm(ina)
    illa = idall(ina)
    iout = 0


    do ii=ipoa,ipoa+ilma-1
       if(abs(cc(ii))<value) cc(ii)=zero
    enddo


    !    call dapac(ina)

    return

    !$$$$  end subroutine daclean
  end subroutine daclean_b !$$$$


  !$$$$ end module dabnew
end module dabnew_b !$$$$
