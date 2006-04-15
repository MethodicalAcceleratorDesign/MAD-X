!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
module S_fitting
  USE MAD_LIKE
  IMPLICIT NONE
  public
  PRIVATE FIND_ORBIT_LAYOUT,FIND_ENV_LAYOUT, FIND_ORBIT_LAYOUT_noda
  logical(lp), PRIVATE :: VERBOSE = .false.
  integer :: max_fit_iter=20, ierror_fit=0


  INTERFACE FIND_ORBIT
     ! LINKED
     ! no use of TPSA
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
     !RETURN QUADRATIC YS
     MODULE PROCEDURE FIND_ENV_LAYOUT

  END INTERFACE

contains
  SUBROUTINE lattice_GET_CHROM(R,my_state,CHROM)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    REAL(DP) CHROM(:)
    integer i,IB
    TYPE(internal_state) state
    real(dp) closed(6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)

    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

    closed=zero
    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED
    CALL INIT(STATE,2,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y
    CHROM(1)=norm%DHDJ%V(1).SUB.'00001'
    CHROM(2)=norm%DHDJ%V(2).SUB.'00001'
    WRITE(6,*) "Fractional Tunes = ",norm%tune(1:2)
    WRITE(6,*) "CHROMATICITIES = ",CHROM
    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)

  end SUBROUTINE lattice_GET_CHROM

  SUBROUTINE FILL_BETA(r,my_state,BETA,IB,DBETA,tune,tune2,a,ai,mat,clos)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    REAL(DP), ALLOCATABLE :: BETA(:,:,:)
    REAL(DP)DBETA,tune(:),tune2(:)
    type(fibre),pointer :: p
    integer i,IB
    TYPE(internal_state) state
    real(dp) closed(6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)
    real(dp) dbetamax,db1,db2
    real(dp),optional :: a(6,6),ai(6,6),mat(6,6),clos(6)

    if(.not.allocated(beta))   ALLOCATE(BETA(2,2,R%N))


    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    if(present(clos)) then
       closed=clos
    else
       closed=zero
    endif
    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED

    DBETA=ZERO
    dbetamax=zero
    CALL INIT(STATE,1,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y
    if(present(a)) then
       a=norm%a_t
       ai=norm%a_t**(-1)
       id=y
       mat=id
    endif
    if(ib==1) then
       tune(1:2)=norm%tune(1:2)
       Write(6,*) " Tunes ", tune(1:2)
    endif
    tune2(1:2)=norm%tune(1:2)

    y=closed+norm%a_t
    p=>r%start
    do i=1,r%n

       CALL TRACK(R,Y,i,i+1,STATE)
       beta(IB,1,i)=(y(1).sub.'1')**2   + (y(1).sub.'01')**2
       beta(iB,2,i)=(y(3).sub.'001')**2 + (y(3).sub.'0001')**2

       IF(IB==2) THEN
          db1=ABS(beta(2,1,i)-beta(1,1,i))/beta(1,1,i)
          db2=ABS(beta(2,2,i)-beta(1,2,i))/beta(1,2,i)
          DBETA=(db1+db2)/TWO+dbeta
          if( db1>dbetamax) dbetamax=db1
          if( db2>dbetamax) dbetamax=db2
       ENDIF
       p=>p%next
    enddo
    DBETA=DBETA/R%N

    IF(IB==2) WRITE(6,*) "<DBETA/BETA> = ",DBETA
    IF(IB==2) WRITE(6,*) "MAXIMUM OF DBETA/BETA = ",dbetamax

    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)
    if(present(clos)) clos=closed

  end SUBROUTINE FILL_BETA

  SUBROUTINE comp_linear2(r,my_state,a,ai,mat,closed)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    TYPE(internal_state) state
    real(dp) closed(6),a(6,6),ai(6,6),mat(6,6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)


    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)

    CALL INIT(STATE,1,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y

    Write(6,*) " Tunes ",norm%tune(1:2)
    A=ZERO
    AI=ZERO
    MAT=ZERO
    a=norm%a_t
    AI=norm%a_t**(-1)
    id=y
    mat=id
    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)

  end SUBROUTINE comp_linear2


  subroutine lattice_fit_TUNE_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,more
    TYPE(TAYLOR), allocatable:: EQ(:)
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer :: neq=2, no=2,nt,j,it
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow
    !    EPSF=.0001
    epsr=abs(epsf)

    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    CALL INIT(STATE,no,NP,BERZ)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO
    CLOSED(:)=zero
    it=0
100 continue
    it=it+1

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED


    CALL INIT(STATE,no,NP,BERZ)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2)

    eq(1)=       ((NORM%dhdj%v(1)).par.'0000')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'0000')-targ(2)
    epsnow=abs(eq(1))+abs(eq(2))
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile

    do i=1,neq
       eq(i)=eq(i)<=c_%npara
    enddo
    do i=1,neq
       call daprint(eq(i),scratchfile)
    enddo
    close(SCRATCHFILE)
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(id)
    CALL KILL(EQ)



    CALL INIT(1,nt)
    call alloc(g,nt)
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=np+1,nt
       call read(g%v(i),scratchfile)
    enddo
    close(SCRATCHFILE)

    call alloc(t)
    do i=1,np
       g%v(i)=one.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(one.mono.j)*t
       enddo
    enddo
    CALL KILL(t)

    g=g.oo.(-1)
    tpsafit(1:nt)=g

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO
    SET_TPSAFIT=.false.

    CALL ELP_TO_EL(R)

    !    write(6,*) " more "
    !    read(5,*) more
    if(it>=max_fit_iter) goto 101
    if(epsnow<=epsr) goto 102
    GOTO 100

101 continue
    write(6,*) " warning did not converge "

102 continue
    CALL KILL_PARA(R)
    deallocate(eq)

  end subroutine lattice_fit_TUNE_gmap

  subroutine lattice_fit_CHROM_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,more
    TYPE(TAYLOR), allocatable:: EQ(:)
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer :: neq=2, no=3,nt,j,it
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,CHROM(2)
    !    EPSF=.0001
    epsr=abs(epsf)

    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

    CALL INIT(STATE,no,NP,BERZ)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO
    CLOSED(:)=zero
    it=0
100 continue
    it=it+1

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED


    CALL INIT(STATE,no,NP,BERZ)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2)
    CHROM(1)=(NORM%dhdj%v(1)).SUB.'00001'
    CHROM(2)=(NORM%dhdj%v(2)).SUB.'00001'
    write(6,*) " CHROM ",CHROM

    eq(1)=       ((NORM%dhdj%v(1)).par.'00001')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'00001')-targ(2)
    epsnow=abs(eq(1))+abs(eq(2))
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile

    do i=1,neq
       eq(i)=eq(i)<=c_%npara
    enddo
    do i=1,neq
       call daprint(eq(i),scratchfile)
    enddo
    close(SCRATCHFILE)
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(id)
    CALL KILL(EQ)



    CALL INIT(1,nt)
    call alloc(g,nt)
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=np+1,nt
       call read(g%v(i),scratchfile)
    enddo
    close(SCRATCHFILE)

    call alloc(t)
    do i=1,np
       g%v(i)=one.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(one.mono.j)*t
       enddo
    enddo
    CALL KILL(t)

    g=g.oo.(-1)
    tpsafit(1:nt)=g

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO
    SET_TPSAFIT=.false.

    CALL ELP_TO_EL(R)

    !    write(6,*) " more "
    !    read(5,*) more
    if(it>=max_fit_iter) goto 101
    if(epsnow<=epsr) goto 102
    GOTO 100

101 continue
    write(6,*) " warning did not converge "

102 continue
    CALL KILL_PARA(R)
    deallocate(eq)

  end subroutine lattice_fit_CHROM_gmap

  subroutine lattice_fit_SEXT_RES_from_a_gmap_vec(R,my_state,EPSF,MRES,POLY,NPOLY,TARG,NP,neq)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,more
    TYPE(TAYLOR), allocatable:: EQ(:)
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer :: neq, no=3,nt,j,it,MRES(3)
    type(damap) id,B,N,NC
    type(gmap) g
    TYPE(TAYLOR)t,Q1,P1,Q2,P2
    TYPE(PBFIELD)H
    TYPE(PBresonance)Hr
    TYPE(VECRESONANCE)HRES_a,hres_n
    TYPE(ONELIEEXPONENT)HC
    real(dp) epsf,epsr,epsnow,CHROM(2),ARES(2),TR,PREC
    INTEGER, ALLOCATABLE :: EXPON(:),JEXP(:)
    integer i_look
    epsr=abs(epsf)
    PREC=c_1d_10
    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

    CALL INIT(STATE,no,NP,BERZ)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO
    CLOSED(:)=zero
    it=0
    EPSNOW=mybig
100 continue
    it=it+1

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED


    CALL INIT(STATE,no,NP,BERZ)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id,B,N,NC)    ! LOOK AT PAGE 143 OF THE BOOK
    CALL ALLOC(H)
    CALL ALLOC(HC)
    CALL ALLOC(HRES_a)
    CALL ALLOC(hres_n)
    ALLOCATE(EXPON(C_%NV))
    ALLOCATE(JEXP(C_%NPARA))
    call alloc(hr)
    WRITE(6,*) " NPARA = ",C_%NPARA
    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,+STATE)

    ! ETIENNE

    NORM=Y


    HRES_a=NORM%a%nonlinear
    hr=NORM%a%pb

    CALL PRINT(HR%COS,6,PREC)
    CALL PRINT(HR%SIN,6,PREC)
    if(neq==4) then
       HRES_n=NORM%normal%nonlinear
       hr=NORM%normal%pb

       CALL PRINT(HR%COS,6,PREC)
    endif
    EXPON=0
    IF(MRES(1)==0) THEN
       expon(3)=iabs(mres(2))-1
       i_look=4
    ELSEIF(MRES(2)==0) THEN
       expon(1)=iabs(mres(1))-1
       i_look=2
    ELSEIF(MRES(1)*MRES(2)>0) THEN
       expon(1)=iabs(mres(1))-1
       i_look=2
       expon(3)= iabs(mres(2))
    ELSE
       expon(1)=iabs(mres(1))-1
       i_look=2
       expon(4)= iabs(mres(2))
    ENDIF

    ARES(1)=(HRES_a%sin%v(i_look).SUB.EXPON)/(-mres(i_look/2))/two
    ARES(2)=(HRES_a%cos%v(i_look).SUB.EXPON)/(mres(i_look/2))/two
    JEXP=0
    JEXP(1:C_%NPARA)=EXPON(1:C_%NPARA)
    CHROM(1)=(HRES_n%sin%v(2).SUB.'01001')/(PI)/two
    CHROM(2)=(HRES_n%sin%v(4).SUB.'00011')/(PI)/two

    !    write(6,*) " CHROM ",CHROM
    WRITE(6,*) " STRENGTH OF RESONANCE = ",ARES
    WRITE(6,*) " CHROMATICITIES  = ",CHROM


    eq(1)=       ( (HRES_a%sin%v(i_look).par.jexp)/(-mres(i_look/2))/two)-targ(1)
    eq(2)=       ((HRES_a%cos%v(i_look).par.jexp)/(mres(i_look/2))/two)-targ(2)
    if(neq==4) then
       eq(3)=       (HRES_n%sin%v(2).par.'01001')/(PI)/two-targ(3)
       eq(4)=       (HRES_n%sin%v(4).par.'00011')/(PI)/two-targ(4)
    endif
    epsnow=zero
    do i=1,neq
       epsnow=abs(eq(i))+epsnow
    enddo

    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile

    do i=1,neq
       eq(i)=eq(i)<=c_%npara
    enddo
    do i=1,neq
       call daprint(eq(i),scratchfile)
    enddo
    close(SCRATCHFILE)
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(HC)
    CALL KILL(id,B,N,NC)
    CALL KILL(EQ)
    CALL KILL(H)
    CALL KILL(HRES_a)
    CALL KILL(hres_n)

    call KILL(hr)

    DEALLOCATE(EXPON)
    DEALLOCATE(JEXP)



    CALL INIT(1,nt)
    call alloc(g,nt)
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=np+1,nt
       call read(g%v(i),scratchfile)
    enddo
    close(SCRATCHFILE)

    call alloc(t)
    do i=1,np
       g%v(i)=one.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(one.mono.j)*t
       enddo
    enddo
    CALL KILL(t)

    g=g.oo.(-1)
    tpsafit(1:nt)=g

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO
    SET_TPSAFIT=.false.

    CALL ELP_TO_EL(R)

    !    write(6,*) " more "
    !    read(5,*) more
    if(it>=max_fit_iter) goto 101
    if(epsnow<=epsr) goto 102
    GOTO 100

101 continue
    write(6,*) " warning did not converge "

102 continue
    CALL KILL_PARA(R)
    deallocate(eq)

  end subroutine lattice_fit_SEXT_RES_from_a_gmap_vec

  subroutine lattice_PRINT_RES_FROM_A(R,my_state,NO,EMIT0,MRES,FILENAME)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    type(damap) id
    TYPE(PBresonance)Hr
    CHARACTER(*) FILENAME
    REAL(DP) PREC,EMIT0(2),STR
    INTEGER NO,MF,MRES(4)

    PREC=c_1d_10

    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)



    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED


    CALL INIT(STATE,no,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(hr)
    call alloc(id)    ! LOOK AT PAGE 143 OF THE BOOK
    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)

    ! ETIENNE

    NORM=Y


    hr=NORM%a%pb

    CALL KANALNUMMER(MF)

    OPEN(UNIT=MF,FILE=FILENAME)

    CALL PRINT(HR%COS,6,PREC)
    CALL PRINT(HR%SIN,6,PREC)
    STR=(HR%COS%H.SUB.MRES)**2+(HR%SIN%H.SUB.MRES)**2; STR=SQRT(STR)
    WRITE(6,*) "RESONANCE = ",MRES,STR
    WRITE(MF,*) "RESONANCE = ",MRES,STR
    CALL PRINT(HR%COS,MF,PREC)
    CALL PRINT(HR%SIN,MF,PREC)
    WRITE(MF,*) " SCALE AT EMIT = ",EMIT0
    ID=1
    ID%V(1)=ID%V(1)*SQRT(EMIT0(1))
    ID%V(2)=ID%V(2)*SQRT(EMIT0(1))
    ID%V(3)=ID%V(3)*SQRT(EMIT0(2))
    ID%V(4)=ID%V(4)*SQRT(EMIT0(2))
    HR%COS%H=HR%COS%H*ID
    HR%SIN%H=HR%SIN%H*ID
    CALL PRINT(HR%COS,MF,PREC)
    CALL PRINT(HR%SIN,MF,PREC)


    CLOSE(MF)
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(hr)

  end subroutine lattice_PRINT_RES_FROM_A

  subroutine lattice_random_error(R,nom,iseed,cut,n,addi,integrated,cn,cns)
    use gauss_dis
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    integer iseed,n,addi,ic,i
    character(nlp) nom
    type(fibre), pointer :: p
    logical(lp) integrated
    real(dp) x,bn,cn,cns,cut


    call gaussian_seed(iseed)

    call context(nom)
    ic=0
    p=>r%start
    do i=1,r%n
       if(p%mag%name==nom) then
          call GRNF(X,cut)
          bn=cn*x
          if(integrated.and.p%mag%p%ld/=zero) then
             bn=(cns+bn)/p%mag%l
          endif
          call add(p,n,addi,bn)
          ic=ic+1
       endif
       p=>P%next
    enddo

    write(6,*) ic," Magnets modified "
  end  subroutine lattice_random_error

  subroutine toggle_verbose
    implicit none
    verbose=.not.verbose
  end   subroutine toggle_verbose





  ! linked

  SUBROUTINE FIND_ORBIT_LAYOUT(RING,FIX,LOC,STATE,TURNS)  ! Finds orbit with TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: RING
    INTEGER, OPTIONAL:: TURNS
    real(dp)  FIX(6),DIX(6),xdix,xdix0,tiny,freq
    TYPE(REAL_8) X(6)
    TYPE(DAMAP) MX,SX,SXI,IS
    integer NO1,ND2,I,IU,LOC,ITE,npara
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat
    TYPE (fibre), POINTER :: C
    logical(lp) APERTURE,c_da
    INTEGER TURNS0
    c_%stable_da=.true.
    c_da=c_%check_da
    c_%check_da=.true.
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.
    TURNS0=1
    freq=zero
    IF(PRESENT(TURNS)) TURNS0=TURNS
    Nullify(C);
    if(.not.ring%closed) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)= " This line is not ring : FIND_ORBIT_LAYOUT "
       call write_E(100)
    endif
    dix(:)=zero
    tiny=c_1d_40
    xdix0=c_1d4*DEPS_tracking
    NO1=1

    if(.not.present(STATE)) then
       IF(default%NOCAVITY) THEN
          !    ND1=2
          stat=default+only_4d
       ELSE
          !   ND1=3
          STAT=default
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(111)
       ENDIF
    else
       IF(STATE%NOCAVITY) THEN
          !    ND1=2
          STAT=STATE+only_4d
       ELSE
          !   ND1=3
          STAT=STATE
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(112)
       ENDIF
    endif
101 continue
    !    ND2=2*ND1
    if(stat%totalpath.and.(.not.stat%nocavity)) then
       C=>RING%START
       freq=zero
       i=1
       xdix=zero
       do while(i<=RING%n)
          if(associated(c%magp%freq)) then
             IF(FREQ==ZERO) THEN
                freq=c%magp%freq
             ELSEIF(c%magp%freq/=ZERO.AND.c%magp%freq<FREQ) THEN
                freq=c%magp%freq
             ENDIF
          endif
          XDIX=XDIX+c%mag%P%LD/c%mag%P%BETA0
          c=>c%next
          i=i+1
       enddo
       if(freq==zero) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(113)
       endif
       IF(RING%HARMONIC_NUMBER>0) THEN
          FREQ=RING%HARMONIC_NUMBER*CLIGHT/FREQ
       ELSE
          XDIX=XDIX*FREQ/CLIGHT
          FREQ=NINT(XDIX)*CLIGHT/FREQ
       ENDIF
    endif
    CALL INIT(STAT,NO1,0,BERZ,ND2,NPARA)

    !  call init(NO1,ND1,0,0,.TRUE.)


    CALL ALLOC(X,6)
    CALL ALLOC(MX)
    CALL ALLOC(SX)
    CALL ALLOC(SXI)
    CALL ALLOC(IS)


3   continue
    X=NPARA
    DO I=1,6
       X(I)=FIX(I)
    ENDDO

    DO I=1,TURNS0
       CALL TRACK(RING,X,LOC,STAT)
       if(.not.check_stable.or.(.not.c_%stable_da)) then
          CALL KILL(X,6)
          CALL KILL(MX)
          CALL KILL(SX)
          CALL KILL(SXI)
          CALL KILL(IS)
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",1
          call write_i

          return
       endif
    ENDDO
    x(6)=x(6)-TURNS0*freq

    IS=1
    MX=X
    SX=MX-IS
    DIX=SX
    DO I=1,ND2
       DIX(I)=FIX(I)-DIX(I)
    enddo
    IS=0
    IS=DIX

    SXI=SX**(-1)
    SX=SXI.o.IS
    DIX=SX
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=zero
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,a72))'
    write(w_p%c(1),'(a22,g20.14)') " Convergence Factor = ",xdix
    if(verbose) call write_i
    if(xdix.gt.deps_tracking) then
       ite=1
    else
       if(xdix.ge.xdix0.or.xdix<=tiny) then
          ite=0
       else
          ite=1
          xdix0=xdix
       endif
    endif
    if(ite.eq.1)  then
       GOTO 3

    endif
    CALL KILL(X,6)
    CALL KILL(MX)
    CALL KILL(SX)
    CALL KILL(SXI)
    CALL KILL(IS)
    c_%check_da=c_da
    c_%APERTURE_FLAG=APERTURE

  END SUBROUTINE FIND_ORBIT_LAYOUT

  integer function FIND_ORBIT_flag(RING,FIX,LOC,STATE,eps,TURNS) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: RING
    real(dp) , intent(inOUT) :: FIX(6)
    INTEGER , intent(in) :: LOC
    INTEGER, OPTIONAL::TURNS
    real(dp) , optional, intent(in) :: eps
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE

    call find_orbit(RING,FIX,LOC,STATE,eps,TURNS)

    call PRODUCE_APERTURE_FLAG(FIND_ORBIT_flag)
    if(.not.c_%stable_da) then
       c_%stable_da=.true.
    endif
    !   resets Da on its own here only

  END function  FIND_ORBIT_flag





  SUBROUTINE FIND_ORBIT_LAYOUT_noda(RING,FIX,LOC,STATE,eps,TURNS) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: RING
    real(dp) , intent(inOUT) :: FIX(6)
    INTEGER , intent(in) :: LOC
    INTEGER, OPTIONAL::TURNS
    real(dp) , optional,intent(in) :: eps
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat

    real(dp)  DIX(6),xdix,xdix0,tiny,freq
    real(dp) X(6),Y(6),MX(6,6),sxi(6,6),SX(6,6)
    integer NO1,ND2,I,IU,ITE,ier,j
    TYPE (fibre), POINTER :: C
    logical(lp) APERTURE
    INTEGER TURNS0,trackflag
    TURNS0=1
    trackflag=0
    IF(PRESENT(TURNS)) TURNS0=TURNS
    freq=zero
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    if(.not.present(eps)) then
       if(.not.present(STATE)) then
          call FIND_ORBIT_LAYOUT(RING,FIX,LOC,TURNS=TURNS0)
       else
          call FIND_ORBIT_LAYOUT(RING,FIX,LOC,STATE,TURNS=TURNS0)
       endif
       return
    endif


    Nullify(C);

    if(.not.ring%closed) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=" This line is not ring : FIND_ORBIT_LAYOUT_noda "
       call write_e(100)
    endif
    dix(:)=zero
    tiny=c_1d_40
    xdix0=c_1d4*DEPS_tracking
    NO1=1
    if(.not.present(STATE)) then
       IF(default%NOCAVITY) THEN
          !    ND1=2
          stat=default+only_4d
       ELSE
          !   ND1=3
          STAT=default
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72))'
          w_p%c(1)=" No Cavity in the Line "
          w_p%c(2)=" FIND_ORBIT_LAYOUT will crash "
          call write_e(101)
       ENDIF
    else
       IF(STATE%NOCAVITY) THEN
          ND2=4
          STAT=STATE+only_4d
       ELSE
          ND2=6
          STAT=STATE
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72))'
          w_p%c(1)=" No Cavity in the Line "
          w_p%c(2)=" FIND_ORBIT_LAYOUT will crash "
          call write_e(112)
       ENDIF
    endif
101 continue


    if(stat%totalpath.and.(.not.stat%nocavity)) then
       C=>RING%START
       freq=zero
       i=1
       xdix=zero
       do while(i<=RING%n)
          if(associated(c%magp%freq)) then
             IF(FREQ==ZERO) THEN
                freq=c%magp%freq
             ELSEIF(c%magp%freq/=ZERO.AND.c%magp%freq<FREQ) THEN
                freq=c%magp%freq
             ENDIF
          endif
          XDIX=XDIX+c%mag%P%LD/c%mag%P%BETA0
          c=>c%next
          i=i+1
       enddo
       if(freq==zero) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(113)
       endif
       IF(RING%HARMONIC_NUMBER>0) THEN
          FREQ=RING%HARMONIC_NUMBER*CLIGHT/FREQ
       ELSE
          XDIX=XDIX*FREQ/CLIGHT
          FREQ=NINT(XDIX)*CLIGHT/FREQ
       ENDIF
    endif




3   continue

    X=FIX

    DO I=1,TURNS0
       !       CALL TRACK(RING,X,LOC,STAT)
       trackflag=TRACK_flag(RING,X,LOC,STAT)
       if(trackflag/=0) then
          CALL RESET_APERTURE_FLAG
          c_%APERTURE_FLAG=APERTURE
          write(6,*) " Unstable in find_orbit without TPSA"

          return
       endif
       if(.not.check_stable) then
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",2
          call write_i

          return
       endif

    ENDDO

    x(6)=x(6)-freq*turns0

    mx=zero
    DO J=1,ND2

       Y=FIX
       Y(J)=FIX(J)+EPS
       DO I=1,TURNS0
          CALL TRACK(RING,Y,LOC,STAT)
          if(.not.check_stable) then
             w_p=0
             w_p%nc=1
             w_p%fc='((1X,a72))'
             write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",3
             call write_i

             return
          endif
       ENDDO
       y(6)=y(6)-freq*turns0

       do i=1,ND2
          MX(I,J)=(Y(i)-X(i))/eps
       enddo

    ENDDO


    SX=MX;
    DO I=1,nd2   !  6 before
       SX(I,I)=MX(I,I)-one
    ENDDO

    DO I=1,ND2
       DIX(I)=FIX(I)-X(I)
    enddo

    CALL matinv(SX,SXI,ND2,6,ier)
    IF(IER==132)  then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=" Inversion failed in FIND_ORBIT_LAYOUT_noda"
       call write_e(333)
       return
    endif

    x=zero
    do i=1,nd2
       do j=1,nd2
          x(i)=sxi(i,j)*dix(j)+x(i)
       enddo
    enddo
    dix=x
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=zero
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,a72))'
    write(w_p%c(1),'(a22,g20.14)') " Convergence Factor = ",xdix
    if(verbose) call write_i
    if(xdix.gt.deps_tracking) then
       ite=1
    else
       if(xdix.ge.xdix0.or.xdix<=tiny) then
          ite=0
       else
          ite=1
          xdix0=xdix
       endif
    endif
    if(ite.eq.1)  then
       GOTO 3

    endif
    c_%APERTURE_FLAG=APERTURE

  END SUBROUTINE FIND_ORBIT_LAYOUT_noda

  SUBROUTINE FIND_ENV_LAYOUT(RING,YS,FIX,LOC,STATE,emit) ! Finds envelope with TPSA in State or compatible state
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    TYPE(layout),INTENT(inOUT):: RING
    real(dp)  FIX(6),DIX(6),xdix0,mat(6,6),flu(6,6)
    TYPE(REAL_8) X(6)
    TYPE(ENV_8),INTENT(INOUT)::YS(6)
    TYPE(DAMAP) MX,SX,SXI,IS
    type(beamenvelope) env
    integer NO1,ND2,ND1 ,I,LOC ,J
    INTEGER  JJ(LNV)
    real(dp),optional, intent(inout)::emit(3)
    type(internal_state),optional, intent(in)::STATE
    type(internal_state) sss
    TYPE (fibre), POINTER :: C
    logical(lp) APERTURE
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    Nullify(C);

    if(.not.ring%closed) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1) = " This line is not ring : FIND_ENV_LAYOUT "
       call write_e(100)
    endif
    dix(:)=zero

    if(present(state)) then
       if(state%nocavity.or.(.not.state%radiation)) then
          sss=(STATE-nocavity0-only_4d0-delta0)+radiation0
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72))'
          w_p%c(1) = " Inputed State modified temporarily "
          w_p%c(2) = " to be compatible with Beam Envelope"
          if(verbose) call write_i
       else
          sss=STATE
       endif
    else
       if(default%nocavity.or.(.not.default%radiation)) then
          sss=(default-nocavity0-only_4d0-delta0)+radiation0
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72))'
          w_p%c(1) = " Inputed State modified temporarily "
          w_p%c(2) = " to be compatible with Beam Envelope"
          if(verbose) call write_i
       else
          sss=default
       endif
    endif

    C=>RING%START
    do i=1,RING%n
       if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
       C=>C%NEXT
    enddo
    w_p=0
    w_p%nc=2
    w_p%fc='((1X,a72))'
    w_p%c(1) = " No Cavity in the Line "
    w_p%c(2) = " FIND_ENV_LAYOUT will crash "
    call write_e(111)
101 continue

    xdix0=c_1d4*DEPS_tracking

    NO1=1
    ND1=3
    ND2=2*ND1
    DO I=1,LNV
       JJ(I)=0
    ENDDO
    call init(NO1,ND1,0,0,doneitt)

    CALL ALLOC(X,6)
    CALL ALLOC(MX)
    CALL ALLOC(SX)
    CALL ALLOC(SXI)
    CALL ALLOC(IS)
    CALL ALLOC(YS,6)
    CALL FIND_ORBIT(RING,FIX,LOC,SSS,c_1d_8)
    X=ND2
    DO I=1,6
       X(I)=FIX(I)
    ENDDO

    YS=X

    CALL TRACK(RING,YS,LOC,SSS)

    DO I=1,6
       DO J=1,6
          JJ(J)=1
          CALL PEK(YS(i)%V%T,JJ,mat(i,j))
          JJ(J)=0
       ENDDO
    ENDDO

    DO I=1,6
       DO J=1,6
          flu(i,j)=YS(i)%e(j)
       ENDDO
    ENDDO

    CALL KILL(X,6)
    CALL KILL(MX)
    CALL KILL(SX)
    CALL KILL(SXI)
    CALL KILL(IS)
    CALL KILL(YS,6)

    call init(2,ND1,0,0,doneitt)

    CALL alloc(YS)
    CALL alloc(MX)
    call alloc(env)

    DO I=1,6
       mx%v(i)=mx%v(i)+fix(i)
       DO J=1,6
          JJ(J)=1
          CALL POK(mx%v(i),JJ,mat(i,j))
          JJ(J)=0
       ENDDO
    ENDDO
    !ys%v=mx
    !do i=1,6
    !ys(i)%v=mx%v(i)       ! this works
    !enddo
    ys=mx       ! implemented in extend_poly.f90
    DO I=1,6
       DO J=1,6
          YS(i)%e(j)=flu(i,j)
       ENDDO
    ENDDO
    env=ys
    ys=env
    DO I=1,6
       DO J=1,6
          JJ(J)=1
          flu(i,j)=ys(i)%sigma0(j)
          JJ(J)=0
       ENDDO
    ENDDO
    if(present(emit)) emit=env%emittance
    CALL kill(env)
    CALL kill(MX)
    CALL kill(YS)
    call init(1,ND1,0,0,doneitt)
    CALL alloc(YS)
    CALL alloc(MX)
    mx=1
    DO I=1,6
       mx%v(i)=mx%v(i)+fix(i)
    ENDDO
    ys=mx       ! implemented in extend_poly.f90
    DO I=1,6
       DO J=1,6
          ys(i)%sigma0(j)=flu(i,j)
       ENDDO
    ENDDO
    CALL kill(MX)


    c_%APERTURE_FLAG=APERTURE

  END SUBROUTINE FIND_ENV_LAYOUT

  SUBROUTINE find_ENVELOPE(RING,YS,A1,FIX,LOC,STATE)
    ! Finds Envelope with TPSA in state supplied by user which may include parameters
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: RING
    TYPE(real_8),INTENT(INOUT)::A1(6)
    TYPE(ENV_8),INTENT(INOUT)::YS(6)
    real(dp),INTENT(INOUT)::FIX(6)
    INTEGER,INTENT(IN) :: LOC
    TYPE(INTERNAL_STATE),INTENT(IN) :: STATE
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    TYPE(NORMALFORM) NORMAL
    logical(lp) APERTURE
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.


    CALL ALLOC(Y)
    CALL ALLOC(ID)
    CALL ALLOC(NORMAL)

    y=6
    y=FIX

    CALL TRACK(RING,Y,LOC,STATE)
    if(.not.check_stable) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       write(w_p%c(1),'(a16)') " find_ENVELOPE 1"
       call write_i
       CALL KILL(NORMAL)
       CALL KILL(ID)
       CALL KILL(Y)

       return
    endif
    normal= y
    y=normal%a1+FIX

    a1=y
    ys=y
    CALL TRACK(RING,Ys,loc,STATE)
    if(.not.check_stable) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       write(w_p%c(1),'(a16)') " find_ENVELOPE 2"
       call write_i
       CALL KILL(NORMAL)
       CALL KILL(ID)
       CALL KILL(Y)

       return
    endif

    y=YS
    id=y
    id=(normal%a1**(-1))*id
    y=id+FIX
    ys=y


    CALL KILL(NORMAL)
    CALL KILL(ID)
    CALL KILL(Y)
    c_%APERTURE_FLAG=APERTURE


  END SUBROUTINE find_ENVELOPE

  SUBROUTINE fit_all_bends(r,state)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: state
    type(fibre),pointer :: p
    integer i

    p=>r%start

    do i=1,r%n
       if(p%mag%p%b0/=zero) call fit_bare_bend(p,state,r%charge)
       p=>p%next
    enddo

  end SUBROUTINE fit_all_bends

  SUBROUTINE fit_bare_bend(f,state,charge)
    IMPLICIT NONE
    TYPE(fibre),INTENT(INOUT):: f
    TYPE(real_8) y(6)
    TYPE(internal_state), intent(in):: state
    integer,optional,target :: charge
    real(dp) kf,x(6),xdix,xdix0,tiny
    integer ite
    tiny=c_1d_40
    xdix0=c_1d4*DEPS_tracking

    KF=ZERO   ;
    F%MAGP%BN(1)%KIND=3
    F%MAGP%BN(1)%I=1

    CALL INIT(1,1,BERZ)

    CALL ALLOC(Y)

3   continue
    X=ZERO
    Y=X
    CALL TRACK(f,Y,+state,CHARGE)
    x=y
    !    write(6,'(A10,6(1X,G14.7))') " ORBIT IS ",x
    kf=-(y(2).sub.'0')/(y(2).sub.'1')
    xdix=abs(y(2).sub.'0')
    f%MAGP%BN(1) = f%MAGP%BN(1)+KF
    f%MAG%BN(1) = f%MAG%BN(1)+KF
    CALL ADD(f,1,1,zero)     !etienne

    if(xdix.gt.deps_tracking) then
       ite=1
    else
       if(xdix.ge.xdix0.or.xdix<=tiny) then
          ite=0
       else
          ite=1
          xdix0=xdix
       endif
    endif
    if(ite.eq.1) GOTO 3

    F%MAGP%BN(1)%KIND=1
    F%MAGP%BN(1)%I=0
    !    write(6,'(A10,1(1X,g14.7))') " BN(1) IS ",    f%MAG%BN(1)


    CALL KILL(Y)


  end SUBROUTINE fit_bare_bend

  SUBROUTINE  track_aperture(r,my_state,beta,dbeta,tuneold,ib,ITMAX,emit0,aper,pos,nturn,FILENAME,FILEtune,FILESMEAR,resmax)
    IMPLICIT NONE
    INTEGER NTE
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    TYPE(layout), intent(inout) :: R
    REAL(DP), ALLOCATABLE :: BETA(:,:,:)
    integer pos,nturn,i,flag,ib,MF,mft,j,resmax,it,I1
    real(dp) closed(6),MAT(6,6),AI(6,6),A(6,6),emit(2),emit0(5),aper(2),x(6),xn(6),dbeta,tuneold(:)
    real(dp) ra(2),tunenew(2)
    CHARACTER(*) FILENAME,FILEtune,FILESMEAR
    real(dp), allocatable :: dat(:,:),dats(:,:),SMEAR(:,:)
    REAL(DP) JMin(2),JMAX(2), tune1(2),tune2(2),tot_tune(2),epsi,scas,scau,scat
    integer itmax
    epsi=one/nturn
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)
    allocate(dat(0:nturn,6),dats(0:nturn,6))
    allocate(SMEAR(ITMAX,8))
    CLOSED=ZERO
    call FILL_BETA(r,my_state,BETA,IB,DBETA,tuneold,tunenew,a,ai,mat,closed)
    write(6,*) " *****************************************************************"
    write(6,*) "        Tracking with Normalized Aperture "
    write(6,*) "        Tunes = ",tunenew(1:2)

    scau=one
    scas=zero
    dats=zero
    SMEAR=ZERO
    it=0
1001 continue
    it=it+1
    write(6,*) " iteration ",it
    IF(IT==ITMAX+1) GOTO 1002

    scat=(emit0(1)+ it*emit0(3))/aper(1)
    !    scat=(scau+scas)/two    ! etienne
    dat=zero

    xn=zero
    JMAX=ZERO
    JMIN=mybig
    emit(1:2)=scat*aper(1:2)
    write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    write(6,*) " Initial emit = ", emit(1:2)," scale = ",scat
    write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    xn(2)=sqrt(emit(1))
    xn(4)=sqrt(emit(2))
    X=zero
    x=(a*xn)+closed
    WRITE(6,*) " INITIAL X"
    WRITE(6,200) X
    !    x=zero
    !    x(1)=xn(2)/sqrt(a(2,1)**2+a(2,2)**2)
    !    x(3)=xn(4)/sqrt(a(4,3)**2+a(4,4)**2)
    ! Etienne

    dat(0,1:4)=xn(1:4)
    ra(1)=xn(1)**2+xn(2)**2
    ra(2)=xn(3)**2+xn(4)**2
    dat(0,5:6)=ra(1:2)

    flag=0
    do i=1,nturn

       flag=track_flag(r,x,pos,state)
       if(flag/=0) exit
       xn=(ai*x)-closed
       ra(1)=xn(1)**2+xn(2)**2
       ra(2)=xn(3)**2+xn(4)**2
       IF(RA(1)>JMAX(1)) JMAX(1)=RA(1)
       IF(RA(2)>JMAX(2)) JMAX(2)=RA(2)
       IF(RA(1)<JMIN(1)) JMIN(1)=RA(1)
       IF(RA(2)<JMIN(2)) JMIN(2)=RA(2)
       if(ra(1)>aper(1)) then
          flag=101
          exit
       endif
       if(ra(2)>aper(2)) then
          flag=102
          exit
       endif
       dat(i,1:4)=xn(1:4)
       dat(i,5:6)=ra(1:2)
    enddo
    WRITE(6,*) "     MAXIMUM RADIUS SQUARE = ",JMAX
    WRITE(6,*) "      MAXIMUM/INITIAL = ",JMAX(1:2)/emit(1:2)
    WRITE(6,*) "     MINIMUM RADIUS SQUARE = ",JMIN
    WRITE(6,*) "      MINIMUM/INITIAL = ",JMIN(1:2)/emit(1:2)
    WRITE(6,*) "     SMEAR = ",TWO*(JMAX-JMIN)/(JMAX+JMIN)
    SMEAR(IT,1:2)=EMIT(1:2)
    SMEAR(IT,3:4)=JMIN(1:2)
    SMEAR(IT,5:6)=JMAX(1:2)
    SMEAR(IT,7:8)=TWO*(JMAX-JMIN)/(JMAX+JMIN)
    if(flag/=0)then
       scau=scat
       IF(fLAG==101) THEN
          write(6,*)  "          UNSTABLE AT X-NORMALIZED APERTURE "
          !   write(mf,*) "UNSTABLE AT X-NORMALIZED APERTURE "
       ELSEIF(FLAG==102) THEN
          !   write(mf,*) "UNSTABLE AT Y-NORMALIZED APERTURE "
          write(6,*) "          UNSTABLE AT Y-NORMALIZED APERTURE "
       ELSE
          !   write(mf,*) "UNSTABLE: DYNAMIC APERTURE "
          write(6,*) "          UNSTABLE: DYNAMIC APERTURE "
       ENDIF
       goto 1002  ! etienne
       if(abs(scau-scas)<=emit0(3)) then
          goto 1002
       else
          goto 1001
       endif
    ELSE
       write(6,*) "          STABLE "

       write(6,*) "tunes of the ray " ,tot_tune(1:2)
       WRITE(6,201) EMIT,APER, TUNEnew(1:2),DBETA

       scas=scat
       dats=dat

       !  RESONANCE
       tot_tune=zero
       xn(1:4)=dats(0,1:4)
       tune1(1)=atan2(-xn(2),xn(1))/twopi
       tune1(2)=atan2(-xn(4),xn(3))/twopi
       if(tune1(1)<zero)  tune1(1)=tune1(1)+one
       if(tune1(2)<zero)  tune1(2)=tune1(2)+one
       DO I1=0,NTURN
          xn(1:4)=dats(i1,1:4)
          tune2(1)=atan2(-xn(2),xn(1))/twopi
          tune2(2)=atan2(-xn(4),xn(3))/twopi
          if(tune2(1)<zero)  tune2(1)=tune2(1)+one
          if(tune2(2)<zero)  tune2(2)=tune2(2)+one
          tune1=tune2-tune1
          if(tune1(1)<zero)  tune1(1)=tune1(1)+one
          if(tune1(2)<zero)  tune1(2)=tune1(2)+one
          tot_tune =tot_tune+tune1
          tune1=tune2
       ENDDO
       tot_tune=tot_tune/nturn
       do i1=0,resmax
          do j=-resmax,resmax
             dbeta=i1*tot_tune(1)+j*tot_tune(2)
             if(abs(dbeta-nint(dbeta))<epsi ) then
                if(i1+j/=0) write(6,*) i1,j,dbeta," <--- here "
             else
                if(i1+j/=0) write(6,*)i1,j,dbeta
             endif
          enddo
       enddo
       !     PAUSE 100
       ! END  RESONANCE



       if(abs(scau-scas)<=emit0(3)) then
          goto 1002
       else
          goto 1001
       endif

    ENDIF


1002 continue
    CALL KANALNUMMER(MF)

    OPEN(UNIT=MF,FILE=FILENAME)

    WRITE(MF,201) EMIT,APER, TUNEnew(1:2),DBETA
    tot_tune=zero
    xn(1:4)=dats(0,1:4)
    tune1(1)=atan2(-xn(2),xn(1))/twopi
    tune1(2)=atan2(-xn(4),xn(3))/twopi
    if(tune1(1)<zero)  tune1(1)=tune1(1)+one
    if(tune1(2)<zero)  tune1(2)=tune1(2)+one

    DO I=0,NTURN
       WRITE(MF,200)DATs(I,1:6)
       xn(1:4)=dats(i,1:4)
       tune2(1)=atan2(-xn(2),xn(1))/twopi
       tune2(2)=atan2(-xn(4),xn(3))/twopi
       if(tune2(1)<zero)  tune2(1)=tune2(1)+one
       if(tune2(2)<zero)  tune2(2)=tune2(2)+one
       tune1=tune2-tune1
       if(tune1(1)<zero)  tune1(1)=tune1(1)+one
       if(tune1(2)<zero)  tune1(2)=tune1(2)+one
       tot_tune =tot_tune+tune1
       tune1=tune2
    ENDDO
    tot_tune=tot_tune/nturn
    CLOSE(MF)
    CALL KANALNUMMER(MFt)
    OPEN(UNIT=MFt,FILE=FILEtune)
    write(mft,*) " tunes = ",tot_tune(1:2)
    do i=0,resmax
       do j=-resmax,resmax
          dbeta=i*tot_tune(1)+j*tot_tune(2)
          if(abs(dbeta-nint(dbeta))<epsi ) then
             if(i+j/=0) write(mft,*) i,j,dbeta," <--- here "
          else
             if(i+j/=0) write(mft,*)i,j,dbeta
          endif
       enddo
    enddo
    CLOSE(mft)

    CALL KANALNUMMER(MF)
    OPEN(UNIT=MF,FILE=FILESMEAR)
    WRITE(MF,*) " ITERATION   EMIT0(1:2)  JMIN(1:2) JMAX(1:2) TWO*(JMAX-JMIN)/(JMAX+JMIN)"
    DO I=1,ITMAX
       WRITE(MF,202)I, SMEAR(I,1:8)
    ENDDO

    CLOSE(mf)

    emit0(1:2)=scas*aper(1:2)
    emit0(4:5)=tunenew(1:2)

202 FORMAT(1X,I4,8(1X,D18.11))
201 FORMAT(9(1X,D18.11))
200 FORMAT(6(1X,D18.11))
    deallocate(dat,dats,SMEAR)
    WRITE(6,*) "     MAXIMUM RADIUS SQUARE = ",JMAX
    WRITE(6,*) "      MAXIMUM/INITIAL = ",JMAX(1:2)/emit(1:2)
    WRITE(6,*) "     MINIMUM RADIUS SQUARE = ",JMIN
    WRITE(6,*) "      MINIMUM/INITIAL = ",JMIN(1:2)/emit(1:2)
    WRITE(6,*) "     SMEAR = ",TWO*(JMAX-JMIN)/(JMAX+JMIN)
    write(6,*) " *****************************************************************"

  end SUBROUTINE  track_aperture


  SUBROUTINE  THIN_LENS_resplit(R,THIN,lim) ! A re-splitting routine
    IMPLICIT NONE
    INTEGER NTE
    TYPE(layout), intent(inout) :: R
    real(dp), OPTIONAL, intent(inout) :: THIN
    real(dp) gg,RHOI,XL,QUAD,THI
    INTEGER M1,M2,M3, MK1,MK2,MK3,limit(2)  !,limit0(2)
    integer, optional :: lim(2)
    logical(lp) MANUAL,eject,doit
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    nullify(C)

    CALL LINE_L(R,doneit)

    MANUAL=.FALSE.
    eject=.FALSE.

    THI=R%THIN
    IF(PRESENT(THIN)) THI=THIN

    IF(THI<=0) MANUAL=.TRUE.


    IF(MANUAL) THEN
       write(6,*) "thi: thin lens factor (THI<0 TO STOP) "
       read(5,*) thi
       IF(THI<0) eject=.true.
    ENDIF

1001 CONTINUE

    limit(1)=4
    limit(2)=18
    if(present(lim)) limit=lim
    !    limit0(1)=limit(1)
    !    limit0(2)=limit(2)

    M1=0
    M2=0
    M3=0
    MK1=0
    MK2=0
    MK3=0
    r%NTHIN=0

    C=>R%START
    do   WHILE(ASSOCIATED(C))

       doit=(C%MAG%KIND==kind2.or.C%MAG%KIND==kind5)
       doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
       DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
       DOIT=DOIT.OR.(C%MAG%KIND==kind17)

       if(doit) then
          select case(C%MAG%P%METHOD)
          CASE(2)
             M1=M1+1
             MK1=MK1+C%MAG%P%NST
          CASE(4)
             M2=M2+1
             MK2=MK2+3*C%MAG%P%NST
          CASE(6)
             M3=M3+1
             MK3=MK3+7*C%MAG%P%NST
          END SELECT
          r%NTHIN=r%NTHIN+1   !C%MAG%NST
       endif

       C=>C%NEXT

    enddo
    write(6,*) "Present of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3

    if(eject) then
       !      limit(1)=limit0(1)
       !      limit(2)=limit0(2)
       return
    endif
    M1=0
    M2=0
    M3=0
    MK1=0
    MK2=0
    MK3=0


    r%NTHIN=0
    r%THIN=THI

    C=>R%START
    do   WHILE(ASSOCIATED(C))

       doit=(C%MAG%KIND==kind2.or.C%MAG%KIND==kind5)
       doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
       DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
       DOIT=DOIT.OR.(C%MAG%KIND==kind17)

       if(doit)  then


          xl=C%MAG%L
          RHOI=C%MAG%P%B0
          IF(C%MAG%P%NMUL>=2) THEN
             QUAD=SQRT(C%MAG%BN(2)**2+C%MAG%AN(2)**2)
          ELSE
             QUAD=zero
          ENDIF
          if(C%MAG%KIND==kind5.or.C%MAG%KIND==kind17) then
             quad=quad+(C%MAG%b_sol)**2/four
          endif

          GG=XL*(RHOI**2+ABS(QUAD))
          GG=GG/THI
          NTE=INT(GG)
          IF(NTE.LT.limit(1)) THEN
             M1=M1+1
             C%MAG%P%METHOD=2
             IF(NTE.EQ.0) NTE=1
             C%MAG%P%NST=NTE
             MK1=MK1+NTE
          ELSEIF(NTE.GE.limit(1).AND.NTE.LT.limit(2)) THEN
             M2=M2+1
             C%MAG%P%METHOD=4
             NTE=NTE/3
             IF(NTE.EQ.0) NTE=1
             C%MAG%P%NST=NTE
             MK2=MK2+NTE*3
          ELSEIF(NTE.GE.limit(2)) THEN
             M3=M3+1
             C%MAG%P%METHOD=6
             NTE=NTE/7
             IF(NTE.EQ.0) NTE=1
             C%MAG%P%NST=NTE
             MK3=MK3+NTE*7
          ENDIF
          r%NTHIN=r%NTHIN+1  !C%MAG%NST


          call add(C%MAG,C%MAG%P%nmul,1,zero)
          call COPY(C%MAG,C%MAGP)
       ENDIF

       C=>C%NEXT
    enddo


    write(6,*) "Present of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3



    IF(MANUAL) THEN
       write(6,*) "thi: thin lens factor (THI<0 TO STOP) "
       read(5,*) thi
       IF(THI<0) THEN
          THI=R%THIN
          !          limit(1)=limit0(1)
          !          limit(2)=limit0(2)
          RETURN
       ELSE
          GOTO 1001
       ENDIF
    ENDIF


    !    limit(1)=limit0(1)
    !    limit(2)=limit0(2)

    CALL RING_L(R,doneit)

  END SUBROUTINE  THIN_LENS_resplit



  SUBROUTINE  print_bn_an(r,n,title,filename)
    implicit none
    type(layout),intent(inout) ::r
    character(*) filename
    type(fibre),pointer ::p
    integer n,i,mf,j,ntot
    character(*) title

    ntot=0
    call kanalnummer(mf)
    open(unit=mf,file=filename)
    p=>r%start
    write(mf,'(a120)') title
    write(mf,*) n
    do i=1,r%n

       if(associated(p%mag%an)) then
          ntot=ntot+1
          write(mf,*) min(n,p%mag%p%nmul),p%mag%name
          do j=1,min(n,p%mag%p%nmul)
             write(mf,*)j,p%mag%bn(j),p%mag%an(j)
          enddo
       endif
       p=>p%next
    enddo


    close(mf)

    write(6,*) ntot," magnets settings saved to maximum order ",n

  end   SUBROUTINE  print_bn_an

  SUBROUTINE  read_bn_an(r,filename)
    implicit none
    type(layout),intent(inout) ::r
    character(*) filename
    type(fibre),pointer ::p
    integer n,i,mf,j,nt,jt,ntot
    character(nlp) nom
    character*120 title
    real(dp), allocatable :: an(:),bn(:)

    ntot=0
    call kanalnummer(mf)
    open(unit=mf,file=filename)

    p=>r%start
    read(mf,'(a120)') title
    write(6,'(a120)') title
    read(mf,*) n
    allocate(an(n),bn(n))
    an=zero;bn=zero;

    do i=1,r%n

       read(mf,*,end=100) nt,nom
       call context(nom)

       do j=1,nt
          read(mf,*)jt,bn(j),an(j)
       enddo

       do jt=1,r%n
          if(nom==p%mag%name) then
             ntot=ntot+1
             do j=nt,1,-1
                call ADD(p,j,0,bn(j))
                call ADD(p,-j,0,an(j))
             enddo
          endif
          if(nom==p%mag%name) exit
          p=>p%next
       enddo

    enddo

100 continue
    write(6,*) ntot," magnets settings read"

    close(mf)
    deallocate(an,bn)

  end   SUBROUTINE  read_bn_an


end module S_fitting
