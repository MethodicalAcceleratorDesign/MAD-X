!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_fitting
  USE MAD_LIKE
  IMPLICIT NONE
  public
  PRIVATE FIND_ORBIT_LAYOUT,FIND_ENV_LAYOUT, FIND_ORBIT_LAYOUT_noda
  logical(lp), PRIVATE :: VERBOSE = .false.
  integer :: max_fit_iter=20, ierror_fit=0, max_fiND_iter=40
  real(dp) :: fuzzy_split=one
  real(dp) :: max_ds=zero
  integer :: resplit_cutting = 0    ! 0 just magnets , 1 magnets as before / drifts separately
  ! 2  space charge algorithm

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


  SUBROUTINE lattice_GET_tune(R,my_state)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    integer i,IB
    TYPE(internal_state) state
    real(dp) closed(6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)

    STATE=my_state

    closed=zero
    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED
    CALL INIT(STATE,1,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y
    WRITE(6,'(a19,3(1x,g20.14))') "Fractional Tunes = ",norm%tune(1:3)
    if(norm%tune(3)/=zero) &
         WRITE(6,'(a20,(1x,g20.14))') "Synchrotron period = ",1.d0/abs(norm%tune(3))
    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)

  end SUBROUTINE lattice_GET_tune


  SUBROUTINE compute_A_4d(r,my_state,filename,pos,del,no,MY_A)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    integer pos,no,imod,ic,I
    TYPE(internal_state) state
    real(dp) closed(6),del
    type(DAMAP) ID
    TYPE(REAL_8) Y(6)
    CHARACTER(*) FILENAME
    TYPE(FIBRE), POINTER :: P
    TYPE(DAMAP) MY_A
    TYPE(NORMALFORM) NORM


    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    closed=zero
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED

    CALL INIT(STATE,no,0,BERZ)

    CALL ALLOC(Y); CALL ALLOC(MY_A); CALL ALLOC(NORM)
    call alloc(id)
    imod=r%n/10
    id=1
    Y=CLOSED+id
    ic=0
    p=>r%start
    do i=1,r%n

       CALL TRACK(R,Y,i,i+1,STATE)
       if(mod(i,imod)==0) then
          ic=ic+1
          write(6,*) ic*10," % done "
       endif
       p=>p%next
    enddo
    id=y

    call kanalnummer(ic)
    open(unit=ic,file=FILENAME)


    call print(ID,IC)
    NORM=ID
    MY_A=NORM%A_T
    WRITE(6,*) " TUNES ",NORM%TUNE(1:2)
    call print(MY_A,IC)
    CLOSE(IC)


    CALL kill(Y)
    call kill(id)
    call kill(NORM)

  end SUBROUTINE compute_A_4d



  SUBROUTINE compute_map_general(r,my_state,filename,pos,del,no)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    TYPE(internal_state) state
    integer pos,no,is,i,mf
    real(dp) closed(6),del,s
    type(DAMAP) ID
    TYPE(REAL_8) Y(6)
    CHARACTER(NLP) NAME
    CHARACTER(*) FILENAME
    TYPE(FIBRE), POINTER :: P
    type(normalform) norm
    type(taylor) betax,betax2
    integer, allocatable :: expo1(:),expo2(:)
    integer, allocatable :: expo3(:),expo4(:)

    call kanalnummer(mf)
    open(unit=mf,file=filename)

    state=my_state   !+nocavity0

    if(state%nocavity) then
       allocate(expo1(4),expo2(4))
       allocate(expo3(4),expo4(4))
    else
       allocate(expo1(6),expo2(6))
       allocate(expo3(6),expo4(6))
    endif

    expo1=zero;expo2=zero;
    expo3=zero;expo4=zero;

    closed=zero
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,my_state,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED

    call init(state,no,c_%np_pol,berz)
    call alloc(id); call alloc(norm);call alloc(y);call alloc(betax,betax2);

    id=1; y=closed+id;
    is=track_flag(r,y,pos,+state)

    if(is/=0) then
       write(6,*) "HELP"
       stop
    endif


    call kill(id); call kill(norm);call kill(y);call kill(betax,betax2);
    deallocate(expo1,expo2)
    deallocate(expo3,expo4)

  END SUBROUTINE compute_map_general


  SUBROUTINE compute_map_4d(r,my_state,filename,pos,del,no)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    integer pos,no,imod,ic,I
    TYPE(internal_state) state
    real(dp) closed(6),del
    type(DAMAP) ID
    TYPE(REAL_8) Y(6)
    CHARACTER(*) FILENAME
    TYPE(FIBRE), POINTER :: P



    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    closed=zero
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED

    CALL INIT(STATE,no,0,BERZ)

    CALL ALLOC(Y)
    call alloc(id)
    imod=r%n/10
    id=1
    Y=CLOSED+id
    ic=0
    p=>r%start
    do i=1,r%n

       CALL TRACK(R,Y,i,i+1,STATE)
       if(mod(i,imod)==0) then
          ic=ic+1
          write(6,*) ic*10," % done "
       endif
       p=>p%next
    enddo
    id=y

    call kanalnummer(ic)
    open(unit=ic,file=FILENAME)


    call print(ID,IC)

    CLOSE(IC)


    CALL kill(Y)
    call kill(id)

  end SUBROUTINE compute_map_4d

  SUBROUTINE FILL_BETA(r,my_state,pos,BETA,IB,DBETA,tune,tune2,a,ai,mat,clos)
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    REAL(DP), ALLOCATABLE :: BETA(:,:,:)
    REAL(DP)DBETA,tune(:),tune2(:)
    type(fibre),pointer :: p
    integer i,IB,pos,mf
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
    if(pos/=0) then
       CALL FIND_ORBIT(R,CLOSED,pos,STATE,c_1d_5)
       write(6,*) "closed orbit "
       write(6,*) CLOSED
    else
       write(6,*) " Using a map "
    endif
    DBETA=ZERO
    dbetamax=zero
    CALL INIT(STATE,1,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)
    if(pos/=0) then

       id=1
       Y=CLOSED+id

       CALL TRACK(R,Y,1,STATE)
       write(6,*) " stability ", c_%check_stable
       id=y
    else
       call kanalnummer(mf)
       open(unit=mf,file='map.dat')
       call read(id,mf)
       close(mf)
    endif
    NORM=id
    if(present(a)) then
       a=zero
       ai=zero
       a=norm%a_t
       ai=norm%a_t**(-1)
       id=y
       mat=id
       if(pos==0)  Write(6,*) " Tunes ",norm%tune(1:2)
    endif
    if(pos/=0) then
       if(ib==1) then
          tune(1:2)=norm%tune(1:2)
       endif
       tune2(1:2)=norm%tune(1:2)
       Write(6,*) " Tunes ", norm%tune(1:2)

       y=closed+norm%a_t
       p=>r%start
       do i=pos,pos+r%n-1

          CALL TRACK(R,Y,i,i+1,STATE)
          beta(IB,1,i-pos+1)=(y(1).sub.'1')**2   + (y(1).sub.'01')**2
          beta(iB,2,i-pos+1)=(y(3).sub.'001')**2 + (y(3).sub.'0001')**2

          IF(IB==2) THEN
             db1=ABS(beta(2,1,i-pos+1)-beta(1,1,i-pos+1))/beta(1,1,i-pos+1)
             db2=ABS(beta(2,2,i-pos+1)-beta(1,2,i-pos+1))/beta(1,2,i-pos+1)
             DBETA=(db1+db2)/TWO+dbeta
             if( db1>dbetamax) dbetamax=db1
             if( db2>dbetamax) dbetamax=db2
          ENDIF
          p=>p%next
       enddo
       DBETA=DBETA/R%N

       IF(IB==2) WRITE(6,*) "<DBETA/BETA> = ",DBETA
       IF(IB==2) WRITE(6,*) "MAXIMUM OF DBETA/BETA = ",dbetamax
    endif
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
    write(6,*) "closed orbit ", CHECK_STABLE
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
    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2), CHECK_STABLE
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

  subroutine lattice_fit_tune_CHROM_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
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
    integer :: neq=4, no=3,nt,j,it
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,tune(2),CHROM(2)
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
    write(6,*) "closed orbit ", CHECK_STABLE
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
    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2), CHECK_STABLE
    tune(1)=(NORM%dhdj%v(1)).SUB.'0000'
    tune(2)=(NORM%dhdj%v(2)).SUB.'0000'
    CHROM(1)=(NORM%dhdj%v(1)).SUB.'00001'
    CHROM(2)=(NORM%dhdj%v(2)).SUB.'00001'
    write(6,*) " CHROM ",CHROM

    eq(1)=       ((NORM%dhdj%v(1)).par.'00000')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'00000')-targ(2)
    eq(3)=       ((NORM%dhdj%v(1)).par.'00001')-targ(3)
    eq(4)=       ((NORM%dhdj%v(2)).par.'00001')-targ(4)
    epsnow=abs(eq(1))+abs(eq(2))+abs(eq(3))+abs(eq(4))
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

  end subroutine lattice_fit_tune_CHROM_gmap


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

  !          call lattice_random_error(my_ring,name,i1,i2,cut,n,addi,integrated,cn,cns,sc)

  subroutine lattice_random_error(R,nom,i1,i2,cut,n,addi,integrated,cn,cns,per)
    use gauss_dis
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    integer n,addi,ic,i,i1,i2,j
    character(nlp) nom
    type(fibre), pointer :: p
    logical(lp) integrated,f1,f2
    real(dp) x,bn,cn,cns,cut,per

    if(i1>i2) then
       Write(6,*) " error i1 > i2 ",i1,i2
       return
    elseif(i2>nlp)then
       Write(6,*) " error i2 > nlp ",i2,nlp
       return
    endif

    call context(nom)

    ic=0
    p=>r%start
    do i=1,r%n
       f1=.false.
       f2=.false.
       if(i1>=0) then
          f1=(p%mag%name(i1:i2)==nom(i1:i2))
       else
          f1=(p%mag%name ==nom )
       endif


       if(f1) then
          call GRNF(X,cut)
          bn=cns+cn*x
          if(integrated.and.p%mag%p%ld/=zero) then
             bn=bn/p%mag%l
          endif
          call add(p,n,addi,bn)
          f2=.true.
       endif

       if(f1.and.per/=zero ) then
          do j=p%mag%p%nmul,1,-1
             call GRNF(X,cut)
             if(n>0) then
                bn=p%mag%bn(j)
                bn=bn*(one+x*per)
             else
                bn=p%mag%an(j)
                bn=bn*(one+x*per)
             endif
             call add(p,n,0,bn)

          enddo
          f2=.true.
       endif
       if(f2) ic=ic+1
       p=>P%next
    enddo

    write(6,*) ic," Magnets modified "

  end  subroutine lattice_random_error

  subroutine toggle_verbose
    implicit none
    verbose=.not.verbose
  end   subroutine toggle_verbose

  subroutine special_alex_main_ring(r,n_name,targ,sc)
    implicit none
    TYPE(layout), target, intent(inout):: R
    integer  i1,i2,I3,I4,it1,it2,it3,it4
    INTEGER I,N,NU,N2,NP2,mf,nt,NP,J,n_name
    type(fibre), pointer :: p1,p2
    TYPE(POL_BLOCK) QC(11)
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    REAL(DP) X(6),targ(2),xx  ,tas(6)
    TYPE(INTERNAL_STATE) STATE
    LOGICAL(LP) U
    type(normalform) nf
    type(gmap) g
    type(TAYLOR) T,eq(5)

    REAL(DP) ALX,ALY,BEX,BEY,NUX,NUY,TA(2),sc

    !targ(1)=22.43d0
    !targ(2)=20.82d0

    N=10
    NU=11
    if(.not.associated(r%t)) then
       write(6,*) " thin lens lattice not made "
       stop 300
    endif
    p1=>r%start
    do i=1,r%n
       if(p1%mag%name(1:3)=='QDX') then
          i1=i
          exit
       endif
       p1=>p1%next
    enddo
    p2=>p1%next
    do i=i1+1,r%n
       if(p2%mag%name(1:3)=='QDX') then
          i2=i
          exit
       endif
       p2=>p2%next
    enddo
    !    call move_to(r,p1,i1)
    !    call move_to(r,p2,i2)

    write(6,*) p1%mag%name,p1%mag%p%nst
    write(6,*) p2%mag%name,p2%mag%p%nst

    IT1=p1%T1%POS+2 + (p1%mag%p%nst/2 )
    IT2=p2%T1%POS+2 + (p2%mag%p%nst/2 )

    write(6,*) i1,IT1,i2,IT2

    if(mod(p1%mag%p%nst,2)/=0.or.mod(p2%mag%p%nst,2)/=0) then
       write(6,*) " Even number of split needed for fitting in Alex_special "
       stop 100
    endif
    !    call move_to(r,p1,i3)
    !    call move_to(r,p2,i4)
    p1=>p2%next
    do i=i2+1,r%n
       if(p1%mag%name(1:3)=='QFS') then
          i1=i
          exit
       endif
       p1=>p1%next
    enddo

    do i=i1,1,-1
       if(p1%mag%name(1:3)=='QDX') then
          i1=i
          exit
       endif
       p1=>p1%previous
    enddo


    p2=>p1%next
    do i=i1+1,r%n
       if(p2%mag%name(1:3)=='QDX') then
          i2=i
          exit
       endif
       p2=>p2%next
    enddo

    write(6,*) p1%mag%name,p1%mag%p%nst
    write(6,*) p2%mag%name,p2%mag%p%nst

    IT3=p1%T1%POS+2 + (p1%mag%p%nst/2 )
    IT4=p2%T1%POS+2 + (p2%mag%p%nst/2 )

    write(6,*) i1,IT3,i2,IT4
    if(mod(p1%mag%p%nst,2)/=0.or.mod(p2%mag%p%nst,2)/=0) then
       write(6,*) " Even number of split needed for fitting in Alex_special "
       stop 101
    endif

    DO I=1,NU
       QC(I)=0
       QC(I)%n_name=n_name
    ENDDO
    QC(1)%NAME='QFX'
    QC(2)%NAME='QDX'
    QC(3)%NAME='QFN'
    QC(4)%NAME='QDN'
    QC(5)%NAME='QFS'
    QC(6)%NAME='QDS'
    QC(7)%NAME='QFT'
    QC(8)%NAME='QFP'
    QC(9)%NAME='QDT'
    QC(10)%NAME='QFR'
    QC(11)%NAME='QDR'

    DO I=1,NU
       QC(I)%IBN(2)=I
    ENDDO

    DO I=1,NU
       R=QC(I)
    ENDDO

111 STATE=DEFAULT0+ONLY_4D0

    X=0.D0
    CALL INIT(STATE,2,NU,BERZ,N2,NP2)
    CALL ALLOC(ID); call alloc(nf);call alloc(Y);call alloc(eq,5);
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    CALL TRACK(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    nf=y
    TA(1)=0.75d0
    TA(2)=0.68d0
    write(6,*) " arc tunes ",nf%tune(1:2)
    eq(1)=nf%dhdj%v(1) !-TA(1)
    eq(2)=nf%dhdj%v(2) !-TA(2)
    eq(3)=nf%dhdj%v(1)

    !   call print(nf%dhdj%v(1),6)
    !   call print(nf%dhdj%v(2),6)

    ! nf%dhdj%v(1)=(nf%dhdj%v(1)<=4)
    ! nf%dhdj%v(2)=(nf%dhdj%v(2)<=4)

    call kanalnummer(mf)
    open(unit=mf,file='eq.txt')
    eq(4)=(nf%A_T%V(1).par.'1000')**2+(nf%A_T%V(1).par.'0100')**2
    eq(5)=(nf%A_T%V(3).par.'0010')**2+(nf%A_T%V(3).par.'0001')**2
    !call print(eq(4),6)
    !call print(eq(5),6)
    ALX=-(nf%A_T%V(1).SUB.'1')*(nf%A_T%V(2).SUB.'1')-(nf%A_T%V(1).SUB.'01')*(nf%A_T%V(2).SUB.'01')
    ALY=-(nf%A_T%V(3).SUB.'001')*(nf%A_T%V(4).SUB.'001')-(nf%A_T%V(3).SUB.'0001')*(nf%A_T%V(4).SUB.'0001')

    !WRITE(6,*) BEX,BEY
    !WRITE(6,*) ALX,ALY
    x=0.d0
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    CALL TRACK(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
    nf=y
    eq(1)=(eq(1)*8 + (1.d0+nf%dhdj%v(1)))*3-targ(1)
    eq(2)=(eq(2)*8 + (1.d0+nf%dhdj%v(2)))*3-targ(2)
    eq(3)=eq(3)-ta(1)

    eq(4)=eq(4)-(nf%A_T%V(1).par.'1000')**2-(nf%A_T%V(1).par.'0100')**2
    eq(5)=eq(5)-(nf%A_T%V(3).par.'0010')**2-(nf%A_T%V(3).par.'0001')**2
    !call print(eq(4),6)
    !call print(eq(5),6)

    do i=1,5
       eq(i)=eq(i)<=4
       xx=eq(i)
       write(6,*) i,xx
    enddo
    do i=1,5
       call print(eq(i),mf)
    enddo

    close(mf)
    CALL kill(ID); call kill(nf);call kill(Y);call kill(eq,5);


    NP=nu
    nt=NP+5

    CALL INIT(1,nt)
    call alloc(g,nt)
    call alloc(T)

    call kanalnummer(mf)
    open(unit=mf,file='eq.txt')
    do i=np+1,nt
       call read(g%v(i),mf)
       g%v(i)=g%v(i)-(1.d0-sc)*(g%v(i).sub.'0')
    enddo

    call alloc(t)
    do i=1,np
       g%v(i)=one.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(one.mono.j)*t
       enddo
    enddo
    CALL KILL(t)
    write(6,*) "stable ", check_stable
    g=g.oo.(-1)
    write(6,*) "stable ", check_stable,nu
    tpsafit(1:nt)=g
    write(6,*) tpsafit(1:nu)
    CALL KILL(G)

    SET_TPSAFIT=.true.
    DO I=1,NU
       R=QC(I)
    ENDDO
    SET_TPSAFIT=.false.
    close(mf)

    WRITE(6,*) " MORE "
    READ(5,*) MF
    IF(MF==1) GOTO 111


    CALL KILL_PARA(R)

    CALL ELP_TO_EL(R)

  end   subroutine special_alex_main_ring


  subroutine special_alex_main_ring1(r,i1,i2,I3,I4,targ,sc)
    implicit none
    TYPE(layout), target, intent(inout):: R
    integer  i1,i2,I3,I4,it1,it2,it3,it4
    INTEGER I,N,NU,N2,NP2,mf,nt,NP,J
    type(fibre), pointer :: p1,p2
    TYPE(POL_BLOCK) QC(10)
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    REAL(DP) X(6),targ(2),xx  ,tas(6)
    TYPE(INTERNAL_STATE) STATE
    LOGICAL(LP) U
    type(normalform) nf
    type(gmap) g
    type(TAYLOR) T,eq(6)

    REAL(DP) ALX,ALY,BEX,BEY,NUX,NUY,TA(2),sc

    !targ(1)=22.43d0
    !targ(2)=20.82d0

    N=10
    NU=4
    if(.not.associated(r%t)) then
       write(6,*) " thin lens lattice not made "
       stop 300
    endif

    call move_to(r,p1,i1)
    call move_to(r,p2,i2)

    write(6,*) p1%mag%name,p1%mag%p%nst
    write(6,*) p2%mag%name,p2%mag%p%nst

    IT1=p1%T1%POS+2 + (p1%mag%p%nst/2 )
    IT2=p2%T1%POS+2 + (p2%mag%p%nst/2 )

    write(6,*) IT1,IT2

    call move_to(r,p1,i3)
    call move_to(r,p2,i4)

    write(6,*) p1%mag%name,p1%mag%p%nst
    write(6,*) p2%mag%name,p2%mag%p%nst

    IT3=p1%T1%POS+2 + (p1%mag%p%nst/2 )
    IT4=p2%T1%POS+2 + (p2%mag%p%nst/2 )

    write(6,*) IT3,IT4

    DO I=1,NU
       QC(I)=0
    ENDDO
    QC(1)%NAME='QFX'
    QC(2)%NAME='QDX'
    QC(3)%NAME='QFN'
    QC(4)%NAME='QDN'

    DO I=1,NU
       QC(I)%IBN(2)=I
    ENDDO

    DO I=1,NU
       R=QC(I)
    ENDDO

111 STATE=DEFAULT0+ONLY_4D0

    X=0.D0
    CALL INIT(STATE,2,NU,BERZ,N2,NP2)
    CALL ALLOC(ID); call alloc(nf);call alloc(Y);
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    CALL TRACK(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    nf=y
    TA(1)=0.75d0
    TA(2)=0.68d0
    write(6,*) " arc tunes ",nf%tune(1:2)
    nf%dhdj%v(1)=nf%dhdj%v(1)-TA(1)
    nf%dhdj%v(2)=nf%dhdj%v(2)-TA(2)
    call print(nf%dhdj%v(1),6)
    call print(nf%dhdj%v(2),6)

    nf%dhdj%v(1)=(nf%dhdj%v(1)<=4)
    nf%dhdj%v(2)=(nf%dhdj%v(2)<=4)

    call kanalnummer(mf)
    open(unit=mf,file='eq.txt')
    call print(nf%dhdj%v(1),mf)
    call print(nf%dhdj%v(2),mf)
    BEX=(nf%A_T%V(1).SUB.'1')**2+(nf%A_T%V(1).SUB.'01')**2
    BEY=(nf%A_T%V(3).SUB.'001')**2+(nf%A_T%V(3).SUB.'0001')**2
    ALX=-(nf%A_T%V(1).SUB.'1')*(nf%A_T%V(2).SUB.'1')-(nf%A_T%V(1).SUB.'01')*(nf%A_T%V(2).SUB.'01')
    ALY=-(nf%A_T%V(3).SUB.'001')*(nf%A_T%V(4).SUB.'001')-(nf%A_T%V(3).SUB.'0001')*(nf%A_T%V(4).SUB.'0001')

    WRITE(6,*) BEX,BEY
    WRITE(6,*) ALX,ALY
    tas(1:2)=TARG(1:2)
    tas(3)=BEX
    tas(4)=BEY
    tas(5)=ALX
    tas(6)=ALY
    close(mf)
    CALL kill(ID); call kill(nf);call kill(Y);

    NP=nu
    nt=NP+2

    CALL INIT(1,nt)
    call alloc(g,nt)
    call alloc(T)

    call kanalnummer(mf)
    open(unit=mf,file='eq.txt')
    do i=np+1,nt
       call read(g%v(i),mf)
       g%v(i)=g%v(i)-(1.d0-sc)*(g%v(i).sub.'0')
    enddo


    call alloc(t)
    do i=1,np
       g%v(i)=one.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(one.mono.j)*t
       enddo
    enddo
    CALL KILL(t)

    write(6,*) "stable ", check_stable
    g=g.oo.(-1)
    write(6,*) "stable ", check_stable,nu
    tpsafit(1:nt)=g
    CALL KILL(G)
    SET_TPSAFIT=.true.
    DO I=1,NU
       R=QC(I)
    ENDDO
    SET_TPSAFIT=.false.
    close(mf)

    WRITE(6,*) " MORE "
    READ(5,*) MF
    IF(MF==1) GOTO 111
    CALL KILL_PARA(R)

    NU=7
    DO I=1,NU
       QC(I)=0
    ENDDO
    QC(1)%NAME='QFS'
    QC(2)%NAME='QDS'
    QC(3)%NAME='QFT'
    QC(4)%NAME='QFP'
    QC(5)%NAME='QDT'
    QC(6)%NAME='QFR'
    QC(7)%NAME='QDR'

    DO I=1,NU
       QC(I)%IBN(2)=I
    ENDDO

    DO I=1,NU
       R=QC(I)
    ENDDO

112 continue
    X=0.D0
    CALL INIT(STATE,2,NU,BERZ,N2,NP2)
    CALL ALLOC(ID); call alloc(nf);call alloc(Y);
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    CALL TRACK(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
    write(6,*) "before stable nf=y", check_stable
    !    global_verbose=.true.
    nf=y
    !   nf%dhdj%v(1)=nf%dhdj%v(1)-0.75_DP
    !   nf%dhdj%v(2)=nf%dhdj%v(2)-0.68_DP
    write(6,*) "stable nf=y", check_stable



    call alloc(eq,4)

    eq(1)=(TA(1)*8 + (1.d0+nf%dhdj%v(1)))*3
    eq(2)=(TA(2)*8 + (1.d0+nf%dhdj%v(2)))*3

    nux=eq(1)
    nuy=eq(2)

    eq(1)=eq(1)-targ(1)
    eq(2)=eq(2)-targ(2)

    WRITE(6,*) NUX,NUY
    WRITE(6,*) NUX-targ(1),NUY-targ(2)

    eq(3)=(nf%A_T%V(1).par.'1000')**2+(nf%A_T%V(1).par.'0100')**2-bex
    eq(4)=(nf%A_T%V(3).par.'0010')**2+(nf%A_T%V(3).par.'0001')**2-bey
    !eq(5)=-alx-(nf%A_T%V(1).par.'1000')*(nf%A_T%V(2).par.'1')-(nf%A_T%V(1).par.'0100')*(nf%A_T%V(2).par.'0100')
    !eq(6)=-aly-(nf%A_T%V(3).par.'0010')*(nf%A_T%V(4).par.'0010')-(nf%A_T%V(3).par.'0001')*(nf%A_T%V(4).par.'0001')

    do i=1,4
       eq(i)=(eq(i)<=4)
       xx=eq(i)
       write(6,*) i,xx,TAS(I)
    enddo
    call kanalnummer(mf)
    open(unit=mf,file='eq.txt')
    do i=1,4
       call print(eq(i),mf)
    enddo
    close(mf)
    CALL kill(ID); call kill(nf);call kill(Y);
    call kill(eq,4)

    NP=nu
    nt=NP+4

    CALL INIT(1,nt)
    call alloc(g,nt)
    call alloc(T)

    call kanalnummer(mf)
    open(unit=mf,file='eq.txt')
    do i=np+1,nt
       call read(g%v(i),mf)
       g%v(i)=g%v(i)-(1.d0-sc)*(g%v(i).sub.'0')
    enddo

    call alloc(t)
    do i=1,np
       g%v(i)=one.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(one.mono.j)*t
       enddo
    enddo
    CALL KILL(t)
    write(6,*) "stable ", check_stable
    g=g.oo.(-1)
    write(6,*) "stable ", check_stable,nu
    tpsafit(1:nt)=g
    write(6,*) tpsafit(1:nu)
    CALL KILL(G)

    SET_TPSAFIT=.true.
    DO I=1,NU
       R=QC(I)
    ENDDO
    SET_TPSAFIT=.false.
    close(mf)

    WRITE(6,*) " MORE "
    READ(5,*) MF
    IF(MF==1) GOTO 112


    CALL ELP_TO_EL(R)

    CALL KILL_PARA(R)

  end   subroutine special_alex_main_ring1




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

    CALL RESET_APERTURE_FLAG
    !    if(.not.c_%stable_da) then
    !       c_%stable_da=.true.
    !    endif
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
    integer NO1,ND2,I,IU,ITE,ier,j,ITEM
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
          stat=default+    only_4d
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



    ITEM=0
3   continue
    ITEM=ITEM+1
    X=FIX

    DO I=1,TURNS0
       !       CALL TRACK(RING,X,LOC,STAT)
       trackflag=TRACK_flag(RING,X,LOC,STAT)
       if(trackflag/=0) then
          !          CALL RESET_APERTURE_FLAG
          !          c_%APERTURE_FLAG=APERTURE
          !          write(6,*) " Unstable in find_orbit without TPSA"
          !          return
          ITEM=MAX_FIND_ITER+100
       endif
       !       if(.not.check_stable) then
       !          w_p=0
       !          w_p%nc=1
       !          w_p%fc='((1X,a72))'
       !          write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",2
       !          call write_i

       !          return
       !       endif

    ENDDO
    !    write(6,*) x
    x(6)=x(6)-freq*turns0

    mx=zero
    DO J=1,ND2

       Y=FIX
       Y(J)=FIX(J)+EPS
       DO I=1,TURNS0
          CALL TRACK(RING,Y,LOC,STAT)
          !          if(.not.check_stable) then
          !             w_p=0
          !             w_p%nc=1
          !             w_p%fc='((1X,a72))'
          !             write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",3
          !             call write_i

          !             return
          !          endif
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
    !    write(6,*) " Convergence Factor = ",nd2,xdix,deps_tracking
    !    pause 123321
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

    if(iteM>=MAX_FIND_ITER)  then
       C_%stable_da=.FALSE.
       IF(iteM==MAX_FIND_ITER+100) THEN
          write(6,*) " Unstable in find_orbit without TPSA"
       ELSE
          write(6,*) " maximum number of iterations in find_orbit without TPSA"
       ENDIF
       ITE=0
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

  SUBROUTINE fit_bare_bend(f,state,charge,next)
    IMPLICIT NONE
    TYPE(fibre),INTENT(INOUT):: f
    TYPE(real_8) y(6)
    TYPE(internal_state), intent(in):: state
    integer,optional,target :: charge
    real(dp) kf,x(6),xdix,xdix0,tiny
    integer ite
    logical(lp), optional :: next
    logical(lp) nex
    nex=my_false
    if(present(next)) nex=next
    tiny=c_1d_40
    xdix0=c_1d4*DEPS_tracking

    KF=ZERO   ;
    F%MAGP%BN(1)%KIND=3
    F%MAGP%BN(1)%I=1
    if(nex) then
       F%next%MAGP%BN(1)%KIND=3
       F%next%MAGP%BN(1)%I=1
    endif

    CALL INIT(1,1,BERZ)

    CALL ALLOC(Y)

3   continue
    X=ZERO
    Y=X
    CALL TRACK(f,Y,+state,CHARGE)
    if(nex) CALL TRACK(f%next,Y,+state,CHARGE)
    x=y
    !    write(6,'(A10,6(1X,G14.7))') " ORBIT IS ",x
    kf=-(y(2).sub.'0')/(y(2).sub.'1')
    xdix=abs(y(2).sub.'0')
    f%MAGP%BN(1) = f%MAGP%BN(1)+KF
    f%MAG%BN(1) = f%MAG%BN(1)+KF
    if(nex) then
       f%next%MAGP%BN(1) = f%next%MAGP%BN(1)+KF
       f%next%MAG%BN(1) = f%next%MAG%BN(1)+KF
    endif

    CALL ADD(f,1,1,zero)     !etienne
    if(nex) CALL ADD(f%next,1,1,zero)     !etienne

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
    if(nex) then
       F%next%MAGP%BN(1)%KIND=1
       F%next%MAGP%BN(1)%I=0
    endif
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
    integer pos,nturn,i,flag,ib,MF,mft,j,resmax,it,I1,no
    real(dp) closed(6),MAT(6,6),AI(6,6),A(6,6),emit(2),emit0(6),aper(2),x(6),xn(6),dbeta,tuneold(:)
    real(dp) ra(2),tunenew(2),xda(lnv)
    CHARACTER(*) FILENAME,FILEtune,FILESMEAR
    real(dp), allocatable :: dat(:,:),dats(:,:),SMEAR(:,:)
    REAL(DP) JMin(2),JMAX(2), tune1(2),tune2(2),tot_tune(2),epsi,scas(2),scau,scat(2)
    integer itmax
    type(damap) id
    type(tree) monkey



    epsi=one/nturn
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)
    allocate(dat(0:nturn,6),dats(0:nturn,6))
    allocate(SMEAR(ITMAX,8))
    CLOSED=ZERO

    call FILL_BETA(r,my_state,pos,BETA,IB,DBETA,tuneold,tunenew,a,ai,mat,closed)
    write(6,*) " *****************************************************************"
    write(6,*) "        Tracking with Normalized Aperture "
    write(6,*) "        Tunes = ",tunenew(1:2)
    if(pos==0) then
       write(6,*) " give no "
       read(5,*) no
       call init(no,2,0,0)
       call alloc(id); call alloc(monkey)
       call kanalnummer(mf)
       open(unit=mf,file='map.dat')
       call read(id,mf)
       id=closed
       monkey=id
       call kill(id)
       close(mf)
       xda=zero
    endif
    scau=one
    scas=zero
    dats=zero
    SMEAR=ZERO
    it=0
    CALL KANALNUMMER(MFt)
    OPEN(UNIT=MFt,FILE=FILEtune)
1001 continue
    it=it+1
    write(6,*) " iteration ",it
    IF(IT==ITMAX+1) GOTO 1002

    scat(1)=(emit0(1)+ it*emit0(3))/aper(1)
    scat(2)=(emit0(2)+ it*emit0(6))/aper(2)
    !    scat=(scau+scas)/two    ! etienne
    dat=zero

    xn=zero
    JMAX=ZERO
    JMIN=mybig
    emit=scat*aper
    write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    write(6,*) " Initial emit = ", emit(1:2)," scale = ",scat
    write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    xn(2)=sqrt(emit(1))
    xn(4)=sqrt(emit(2))
    X=zero

    x=matmul(a,xn)+closed

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

       if(pos/=0) then
          flag=track_flag(r,x,pos,state)
       else
          flag=0
          xda(1:4)=x(1:4)
          xda=monkey*xda
          x(1:4)=xda(1:4)
       endif
       write(80,*) x(1:4)
       if(flag/=0) exit
       xn=matmul(ai,(x-closed))
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
       ! scau=scat
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
       if(abs(scau-scas(1))<=emit0(3)) then
          goto 1002
       else
          goto 1001
       endif
    ELSE
       write(6,*) "          STABLE "

       !       write(6,*) "tunes of the ray " ,tot_tune(1:2)
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
       write(mft,*) " emit = ",emit(1:2)
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
       !     PAUSE 100
       ! END  RESONANCE



       if(abs(scau-scas(1))<=emit0(3)) then
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

    CLOSE(mft)

    CALL KANALNUMMER(MF)
    OPEN(UNIT=MF,FILE=FILESMEAR)
    WRITE(MF,*) " ITERATION   EMIT0(1:2)  JMIN(1:2) JMAX(1:2) TWO*(JMAX-JMIN)/(JMAX+JMIN)"
    DO I=1,ITMAX
       WRITE(MF,202)I, SMEAR(I,1:8)
    ENDDO

    CLOSE(mf)

    emit0(1:2) =scas*aper
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
    if(pos==0) call kill(monkey)
  end SUBROUTINE  track_aperture


  SUBROUTINE  THIN_LENS_resplit(R,THIN,even,lim,lmax) ! A re-splitting routine
    IMPLICIT NONE
    INTEGER NTE
    TYPE(layout), intent(inout) :: R
    real(dp), OPTIONAL, intent(inout) :: THIN,lmax
    real(dp) gg,RHOI,XL,QUAD,THI,lm,dl
    INTEGER M1,M2,M3, MK1,MK2,MK3,limit(2),parity,inc,nst_tot,ntec  !,limit0(2)
    integer, optional :: lim(2)
    logical(lp) MANUAL,eject,doit
    TYPE (fibre), POINTER :: C
    logical(lp),optional :: even
    logical(lp) doneit
    nullify(C)
    parity=0
    inc=0
    lm=1.0e38_dp
    ntec=0
    max_ds=zero
    if(present(lmax)) lm=abs(lmax)
    if(present(even)) then
       inc=1
       if(even) then
          parity=0
       else
          parity=1
       endIf
    endif
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
    nst_tot=0

    C=>R%START
    do   WHILE(ASSOCIATED(C))

       doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind5)
       doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
       DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
       DOIT=DOIT.OR.(C%MAG%KIND==kind17)
       doit=doit.and.C%MAG%recut

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
       NST_tot=NST_tot+C%MAG%P%nst

       C=>C%NEXT

    enddo
    write(6,*) "Previous of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3
    write(6,*)   "Total NST ", NST_tot

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
    ! CAVITY FOCUSING
    ! TEAPOT SPLITTING....

    r%NTHIN=0
    r%THIN=THI

    nst_tot=0
    C=>R%START
    do   WHILE(ASSOCIATED(C))   !




       !       if(doit)  then

       select case(resplit_cutting)

       case(0)

          doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
          doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
          DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
          DOIT=DOIT.OR.(C%MAG%KIND==kind17)
          doit=doit.and.C%MAG%recut

          if(doit) then
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
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=2
                MK1=MK1+NTE
             ELSEIF(NTE.GE.limit(1).AND.NTE.LT.limit(2)) THEN
                M2=M2+1
                NTE=NTE/3
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=4
                MK2=MK2+NTE*3
             ELSEIF(NTE.GE.limit(2)) THEN
                M3=M3+1
                NTE=NTE/7
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=6
                MK3=MK3+NTE*7
             ENDIF

             r%NTHIN=r%NTHIN+1  !C%MAG%NST
             !         write(6,*)"nte>ntec", nte,ntec
             call add(C%MAG,C%MAG%P%nmul,1,zero)
             call COPY(C%MAG,C%MAGP)
             if(gg>zero) then
                if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst
             endif

          endif  ! doit

       case(1)

          doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
          doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
          DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
          DOIT=DOIT.OR.(C%MAG%KIND==kind17)
          doit=doit.and.C%MAG%recut

          if(doit) then
             xl=C%MAG%L
             RHOI=C%MAG%P%B0
             IF(C%MAG%P%NMUL>=2) THEN !
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
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=2
                MK1=MK1+NTE
             ELSEIF(NTE.GE.limit(1).AND.NTE.LT.limit(2)) THEN
                M2=M2+1
                NTE=NTE/3
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=4
                MK2=MK2+NTE*3
             ELSEIF(NTE.GE.limit(2)) THEN
                M3=M3+1
                NTE=NTE/7
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=6
                MK3=MK3+NTE*7
             ENDIF


             r%NTHIN=r%NTHIN+1  !C%MAG%NST

             if(present(lmax).and.c%mag%kind==kind1) then
                dl=(C%MAG%P%ld/C%MAG%P%nst)
                if(dl>lm*fuzzy_split) then
                   ntec=int(C%MAG%P%ld/lm)+1
                   if(mod(nte,2)/=parity) ntec=ntec+inc
                   C%MAG%P%NST=ntec
                endif
             endif

             call add(C%MAG,C%MAG%P%nmul,1,zero)
             call COPY(C%MAG,C%MAGP)
             !             if(gg>zero) then
             if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst
             !             endif
          endif

       case(2)


          doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
          doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
          DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
          DOIT=DOIT.OR.(C%MAG%KIND==kind17)
          doit=doit.and.C%MAG%recut

          if(doit) then
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
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                if(nte>ntec.or.(.not.present(lmax)) ) then
                   C%MAG%P%NST=NTE
                   C%MAG%P%METHOD=2
                endif
                MK1=MK1+NTE
             ELSEIF(NTE.GE.limit(1).AND.NTE.LT.limit(2)) THEN
                M2=M2+1
                NTE=NTE/3
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                if(nte>ntec.or.(.not.present(lmax)) ) then
                   C%MAG%P%NST=NTE
                   C%MAG%P%METHOD=4
                endif
                MK2=MK2+NTE*3
             ELSEIF(NTE.GE.limit(2)) THEN
                M3=M3+1
                NTE=NTE/7
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                if(nte>ntec.or.(.not.present(lmax)) ) then
                   C%MAG%P%NST=NTE
                   C%MAG%P%METHOD=6
                endif
                MK3=MK3+NTE*7
             ENDIF


             r%NTHIN=r%NTHIN+1  !C%MAG%NST
             !         write(6,*)"nte>ntec", nte,ntec
             if(nte>ntec.or.(.not.present(lmax)) ) then
                call add(C%MAG,C%MAG%P%nmul,1,zero)
                call COPY(C%MAG,C%MAGP)
             endif
             !            if(gg>zero) then
             !               if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst
             !            endif

             if(present(lmax)) then
                dl=(C%MAG%P%ld/C%MAG%P%nst)
                if(dl>lm*fuzzy_split.and.C%MAG%KIND/=kindpa) then
                   nte=int(C%MAG%P%ld/lm)+1
                   if(mod(nte,2)/=parity) nte=nte+inc
                   if(nte > C%MAG%P%NST ) then
                      C%MAG%P%NST=nte
                      call add(C%MAG,C%MAG%P%nmul,1,zero)
                      call COPY(C%MAG,C%MAGP)
                   endif

                elseif(dl>lm.and.C%MAG%KIND==kindpa) then
                   write(6,*) " Pancake cannot be recut "
                endif
             endif
             if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst


          endif

       case default
          stop 988
       end select


       !      endif
       NST_tot=NST_tot+C%MAG%P%nst
       C=>C%NEXT
    enddo   !   end of do   WHILE


    write(6,*) "Present of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3
    write(6,*)   "Total NST ", NST_tot
    write(6,*)   "Biggest ds ", max_ds



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

  SUBROUTINE  THIN_LENS_restart(R) ! A re-splitting routine
    IMPLICIT NONE
    INTEGER NTE
    TYPE(layout), intent(inout) :: R
    real(dp) gg,RHOI,XL,QUAD,THI,lm,dl
    INTEGER M1,M2,M3, MK1,MK2,MK3,limit(2),parity,inc,nst_tot,ntec  !,limit0(2)
    logical(lp) doit
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    nullify(C)

    CALL LINE_L(R,doneit)






    M1=0
    M2=0
    M3=0
    MK1=0
    MK2=0
    MK3=0


    r%NTHIN=0

    nst_tot=0
    C=>R%START
    do   WHILE(ASSOCIATED(C))


       doit=(C%MAG%KIND==kind1.and.C%MAG%KIND==kind2.or.C%MAG%KIND==kind5.or.C%MAG%KIND==kind4)
       doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
       DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
       DOIT=DOIT.OR.(C%MAG%KIND==kind17)


       if(doit)  then

          M1=M1+1
          NTE=1
          C%MAG%P%NST=NTE
          MK1=MK1+NTE
          C%MAG%P%METHOD=2

          r%NTHIN=r%NTHIN+1  !C%MAG%NST

          call add(C%MAG,C%MAG%P%nmul,1,zero)
          call COPY(C%MAG,C%MAGP)
       else
          if(C%MAG%KIND/=kindpa) then
             C%MAG%P%NST=1
             if(associated(C%MAG%bn))call add(C%MAG,C%MAG%P%nmul,1,zero)
             call COPY(C%MAG,C%MAGP)
          endif
       endif

       NST_tot=NST_tot+C%MAG%P%nst
       C=>C%NEXT
    enddo


    write(6,*) "Present of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3
    write(6,*)   "Total NST ", NST_tot



    CALL RING_L(R,doneit)

  END SUBROUTINE  THIN_LENS_restart


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

       if(associated(p%mag%an)) then
          read(mf,*) nt,nom
          call context(nom)

          do j=1,nt
             read(mf,*)jt,bn(j),an(j)
          enddo

          ntot=ntot+1
          do j=nt,1,-1
             call ADD(p,j,0,bn(j))
             call ADD(p,-j,0,an(j))
          enddo
       endif  ! associated
       p=>p%next
    enddo

    write(6,*) ntot," magnets settings read"

    close(mf)
    deallocate(an,bn)

  end   SUBROUTINE  read_bn_an

  ! THIN LENS EXAMPLE

  SUBROUTINE assign_one_aperture(L,pos,kindaper,R,X,Y)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: L
    integer pos,kindaper
    REAL(DP) R,X,Y
    type(fibre), pointer :: P

    call move_to(L,p,pos)

    if(.NOT.ASSOCIATED(P%MAG%p%aperture)) THEN
       call alloc(P%MAG%p%aperture)
       call alloc(P%MAGP%p%aperture)
    ENDIF
    if(kindaper/=0) then
       P%MAG%p%aperture%kind = kindaper
       P%MAGP%p%aperture%kind = kindaper
       P%MAG%p%aperture%r    = R
       P%MAG%p%aperture%x    = X
       P%MAG%p%aperture%y    = y
       P%MAGP%p%aperture%r    = R
       P%MAGP%p%aperture%x    = X
       P%MAGP%p%aperture%y    = y
    endif

  end SUBROUTINE assign_one_aperture

  SUBROUTINE TURN_OFF_ONE_aperture(R,pos)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer pos
    type(fibre), pointer :: P

    call move_to(r,p,pos)

    if(ASSOCIATED(P%MAG%p%aperture)) THEN
       P%MAG%p%aperture%kind = -P%MAG%p%aperture%kind
       P%MAGP%p%aperture%kind = P%MAG%p%aperture%kind
    ENDIF

  end SUBROUTINE TURN_OFF_ONE_aperture

  SUBROUTINE MESS_UP_ALIGNMENT(R,SIG,cut)
    use gauss_dis
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer J,I
    type(fibre), pointer :: P
    REAL(DP) SIG(:),X,MIS(6),cut


    p=>r%start
    do i=1,r%n

       IF(P%MAG%KIND/=KIND0.AND.P%MAG%KIND/=KIND1) THEN
          DO J=1,6
             call GRNF(X,cut)
             MIS(J)=X*SIG(J)
          ENDDO
          call MISALIGN_FIBRE(p,mis)
       ENDIF
       P=>P%NEXT
    ENDDO
  end SUBROUTINE MESS_UP_ALIGNMENT


  !          CALL MESS_UP_ALIGNMENT_name(my_ring,name,i1,i2,SIG,cut)

  subroutine MESS_UP_ALIGNMENT_name(R,nom,i1,i2,sig,cut)
    use gauss_dis
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    integer i1,i2,j,ic,i
    character(nlp) nom
    type(fibre), pointer :: p
    logical(lp) integrated,f1
    real(dp) cut,sig(6),mis(6),x

    if(i1>i2) then
       Write(6,*) " error i1 > i2 ",i1,i2
       return
    elseif(i2>nlp) then
       Write(6,*) " error i2 > nlp ",i2,nlp
       return
    endif

    call context(nom)

    ic=0


    p=>r%start
    do i=1,r%n

       IF(P%MAG%KIND/=KIND0.AND.P%MAG%KIND/=KIND1) THEN
          f1=.false.
          if(i1>=0) then
             f1=(p%mag%name(i1:i2)==nom(i1:i2))
          else
             f1=(p%mag%name ==nom )
          endif
          if(f1) then
             ic=ic+1
             DO J=1,6
                call GRNF(X,cut)
                MIS(J)=X*SIG(J)
             ENDDO
             call MISALIGN_FIBRE(p,mis)
          endif
       ENDIF
       P=>P%NEXT
    ENDDO

    write(6,*) ic," Magnets misalgned "

  end  subroutine MESS_UP_ALIGNMENT_name





  SUBROUTINE dyn_aper(L,x_in,n_in,ang_in,ang_out,del_in,dlam,pos,nturn,ite,state,mf)
    IMPLICIT NONE
    type(layout), intent(inout) :: L
    real(dp) x(6)
    REAL(DP) x_in,del_in,closed(6),r(6),rt(6)
    REAL(DP) lamT,lams,lamu,dlam,DLAMT,DX,ang,ang_in,ang_out
    integer pos,nturn,i,st,ite,ic,mf,J,n_in,j_in
    TYPE(INTERNAL_STATE) STATE
    TYPE(FIBRE), POINTER :: P
    !
    !    TYPE(REAL_8) Y(6)
    !    TYPE(DAMAP) ID
    !    TYPE(NORMALFORM) NORM

    closed=zero
    !    STATE=STATE+NOCAVITY0
    if(state%nocavity) closed(5)=del_in

    CALL FIND_ORBIT(L,CLOSED,pos,STATE,c_1d_5)
    write(6,*) "closed orbit "
    write(6,*) CLOSED
    write(mf,201) closed
    ang= (ang_out-ang_in)/n_in
    lamt=one
    do j_in=0,n_in

       x=zero
       x(1)=x_in*cos(j_in*ang+ang_in)
       x(3)=x_in*sin(j_in*ang+ang_in)
       x(5)=del_in


       dx=0.3_dp

       r=zero;rt=zero;
       lams=zero
       lamu=ZERO

       DLAMT=DX

       !    lamt=ONE
       ic=0
       do while(DLAMT>dlam.and.ic<ite)

          ic=ic+1
          R=ZERO;
          r(1:4)=lamt*x(1:4)
          if(state%nocavity) then
             rt=r+closed
          else
             rt=r+closed
             rt(5)=rt(5)+x(5)
          endif


          do i=1,nturn
             st=track_flag(L,rt,pos,state)
             if(st/=0) exit
          enddo

          if(st/=0) then
             lamu=lamt
             lamt=(lams+lamt)/two
          else
             lams=lamt
             IF(LAMU<DX) THEN
                lamt=DX+lamt
             ELSE
                lamt=(lamu+lamt)/two
             ENDIF
          endif
          DLAMT=sqrt(x(1)**2+x(3)**2)*ABS(LAMU-LAMS)
       enddo
       write(6,*) ic,(j_in*ang+ang_in)/twopi,lamS*x(1),lamS*x(3)

       write(mf,202) lamS*x(1),lamS*x(3),lamS*x(1)+closed(1),lamS*x(3)+closed(3),DLAMT
       lamt=lamt*0.8_dp
    enddo
201 FORMAT(6(1X,D18.11))
202 FORMAT(5(1X,D18.11))

  end SUBROUTINE dyn_aper



  SUBROUTINE THIN_EXAMPLE(R,B,I1,I2,IN_STATE,MF)
    IMPLICIT NONE
    TYPE(BEAM), INTENT(INOUT) :: B(:)
    TYPE(LAYOUT),TARGET :: R
    integer i1,i2,MF,i,IT1,IT2,INSIDE_POS1,INSIDE_POS2
    type(fibre), pointer :: FIBRE1,FIBRE2
    TYPE(INTEGRATION_NODE),POINTER :: SLICE1,SLICE2
    TYPE(INTERNAL_STATE) IN_STATE
    TYPE(INTERNAL_STATE) STATE
    real(dp) X(6),S_INITIAL,S_FINAL,TOTAL_LENGTH

    STATE=IN_STATE

    CALL COPY_BEAM(B(1),B(2))





    !  TRACKING AS IN PTC PROPER

    !   We locate the I1th fibre and the I2th


    call move_to(r,FIBRE1,i1)
    call move_to(r,FIBRE2,i2)


    Write(MF,*)"Magnets ",i1," AND ",i2," are ",FIBRE1%mag%name(1:LEN_TRIM(FIBRE1%mag%name))," AND ",&
         FIBRE2%mag%name(1:LEN_TRIM(FIBRE2%mag%name))

    ! WE TRACK FROM ENTRANCE OF FIBRE1 TO ENTRANCE OF FIBRE2
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,P1=FIBRE1,P2=FIBRE2 ) "
    WRITE(MF,*) "                                            "
    CALL TRACK(R,B(2),STATE,P1=FIBRE1,P2=FIBRE2 )   ! LOOK AT ROUTINE TRACK_LAYOUT_12  COMMENTS BELOW



    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),6)

    ! Regular PTC tracking: it should agree
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "___________________ ORDINARY PTC RESULTS    ____________________________"

    DO I=1,B(2)%N
       X=B(1)%X(I,1:6)      ! B(2) WAS SAVED IN B(1)
       CALL TRACK(R,X,I1,I2,STATE)
       WRITE(MF,*) "_________________________________________________________________________"
       WRITE(MF,*) " PARTICLE # ",I
       WRITE(MF,*) " TIME  = ",X(6)
       WRITE(MF,*) " X,Y = ", X(1),X(3)
       WRITE(MF,*) " PX,PY = ",X(2),X(4)
       WRITE(MF,*) " ENERGY VARIABLE = ",X(5)

    ENDDO


!!!!!!!
    ! THESE TWO FIBRES START AT THE SLICE NUMBER IT1 AND IT2 RESPECTIVELY
    !
    CALL COPY_BEAM(B(1),B(2))

    IT1=FIBRE1%T1%POS
    IT2=FIBRE2%T1%POS

    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,POS1=IT1,POS2=IT2 ) "
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  INDEX OF THE FIRST SLICE   =",IT1
    WRITE(MF,*) "  INDEX OF THE FINAL SLICE   =",IT2
    CALL TRACK(R,B(2),STATE,POS1=IT1,POS2=IT2 ) ! LOOK AT ROUTINE TRACK_LAYOUT_12  COMMENTS BELOW
    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),MF)

    ! TRACKING USING THE POINTERS SLICE1 AND SLICE2 RESPECTIVELY
    !
    CALL COPY_BEAM(B(1),B(2))

    SLICE1=>FIBRE1%T1
    SLICE2=>FIBRE2%T1

    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,T1=SLICE1,T2=SLICE2 ) "
    WRITE(MF,*) "                                            "
    CALL TRACK(R,B(2),STATE,T1=SLICE1,T2=SLICE2 ) ! LOOK AT ROUTINE TRACK_LAYOUT_12  COMMENTS BELOW
    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),MF)

    ! TRACKING USING THE POINTERS S_INITIAL AND S_FINAL RESPECTIVELY
    !
    CALL COPY_BEAM(B(1),B(2))

    S_INITIAL = SLICE1%S(3)
    S_FINAL   = SLICE2%S(3)
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,S1=S_INITIAL,S2=S_FINAL ) "
    WRITE(MF,*) "                                            "
    WRITE(MF,*) " INITIAL POSITION IN METRES =",S_INITIAL
    WRITE(MF,*) " FINAL POSITION IN METRES   =",S_FINAL
    CALL TRACK_LAYOUT_S12(R,B(2),STATE,S1=S_INITIAL,S2=S_FINAL ) ! LOOK AT ROUTINE TRACK_LAYOUT_S12  COMMENTS BELOW
    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),MF)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$  going inside a magnet  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    WRITE(MF,*) "                                            "
    write(mf,*) fibre1%mag%name," has ",fibre1%mag%p%nst," steps "
    write(mf,*) fibre2%mag%name," has ",fibre2%mag%p%nst," steps "
    WRITE(MF,*) "                                            "
    write(mf,*) "$$$$$$$$$$$$$$$$ STARTING AND ENDING IN THE MIDDLE   $$$$$$$$$$$$$$$$$$$$$$$$"

    INSIDE_POS1=3+fibre1%mag%p%nst/2
    INSIDE_POS2=3+fibre2%mag%p%nst/2
    CALL COPY_BEAM(B(1),B(2))
    ! WE TRACK FROM ENTRANCE OF FIBRE1,INSIDE_POS1 TO ENTRANCE OF FIBRE2,INSIDE_POS2
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,P1=FIBRE1,IN_P1=INSIDE_POS1,P2=FIBRE2,IN_P2=INSIDE_POS2 )  "
    WRITE(MF,*) "                                            "
    CALL TRACK(R,B(2),STATE,P1=FIBRE1,IN_P1=INSIDE_POS1,P2=FIBRE2,IN_P2=INSIDE_POS2 )   ! LOOK AT ROUTINE TRACK_LAYOUT_12  COMMENTS BELOW
    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),MF)

    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$  going inside a magnet using s $$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"

    S_INITIAL = SLICE1%S(3)+ FIBRE1%MAG%L/TWO
    S_FINAL   = SLICE2%S(3)+ FIBRE2%MAG%L/TWO
    CALL COPY_BEAM(B(1),B(2))
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,S1=S_INITIAL,S2=S_FINAL ) "
    WRITE(MF,*) "                                            "
    WRITE(MF,*) " INITIAL POSITION IN METRES =",S_INITIAL
    WRITE(MF,*) " FINAL POSITION IN METRES   =",S_FINAL
    CALL TRACK_LAYOUT_S12(R,B(2),STATE,S1=S_INITIAL,S2=S_FINAL ) ! LOOK AT ROUTINE TRACK_LAYOUT_S12  COMMENTS BELOW
    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),MF)

    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$  NOW WE LOOK AT TIME TRACKING  $$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$    THE AGREEMENT IS PERFECT    $$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$ IF WE LOOK IN A DRIFT SECTION  $$$$$$$$$$$$$$$$$$$$$$"
    write(mf,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    WRITE(MF,*) " WE USE TOTAL TOTAL TIME BECAUSE THIS IS THE ONLY THING IMPLEMENTED WITH TIME TRACKING"
    WRITE(MF,*) " STATE=IN_STATE+TOTALPATH0 "
    STATE=IN_STATE+TOTALPATH0

    CALL COPY_BEAM(B(1),B(2))

    ! FIRST MOVE THE BEAM 1.0 METRE INSIDE D1
    S_INITIAL = zero
    S_FINAL   = ONE
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(R,B(2),STATE,S1=S_INITIAL,S2=S_FINAL ) "
    WRITE(MF,*) "                                            "
    WRITE(MF,*) " INITIAL POSITION IN METRES =",S_INITIAL
    WRITE(MF,*) " FINAL POSITION IN METRES   =",S_FINAL
    CALL TRACK_LAYOUT_S12(R,B(2),STATE,S1=S_INITIAL,S2=S_FINAL ) ! LOOK AT ROUTINE TRACK_LAYOUT_S12  COMMENTS BELOW
    ! WE PRINT THE RESULTS
    CALL PRINT_beam(B(2),MF,I=1)
    CALL COPY_BEAM(B(2),B(1))

    ! THEN TRACK FOR A TIME EQUAL TO THE TOTAL LENGTH OF THE MACHINE: THIS PUTS US BACK SOMEWHERE IN D1
    TOTAL_LENGTH=R%T%END%S(3)
    WRITE(MF,*) "                                            "
    WRITE(MF,*) "  TRACK(B(2),TOTAL_LENGTH,STATE) "
    WRITE(MF,*) "                                            "

    CALL TRACK(B(2),TOTAL_LENGTH,STATE)       !  TRACK_THIN_T(B,DT,K)

    CALL PRINT_beam(B(2),MF,I=1)

    WRITE(MF,*) " NOW WE REPRODUCE THIS WITH 'S' TRACKING TO SEE IF WE HAVE CONSISTENCY"
    WRITE(MF,*) " WE DO ONLY THE FIRST PARTICLE: SAME CHECK WOULD APPLY FOR PARTICLE #2"

    WRITE(MF,*) "  S_INITIAL = ONE "
    WRITE(MF,*) "  S_FINAL   = TOTAL_LENGTH+B(2)%POS(1)%THINLENS%S(3)+B(2)%X(1,7)"
    WRITE(MF,*) "  CALL TRACK(R,B(1),STATE,S1=S_INITIAL,S2=S_FINAL )"


    S_INITIAL = ONE
    S_FINAL   = TOTAL_LENGTH+B(2)%POS(1)%NODE%S(3)+B(2)%X(1,7)

    CALL TRACK_LAYOUT_S12(R,B(1),STATE,S1=S_INITIAL,S2=S_FINAL ) ! LOOK AT ROUTINE TRACK_LAYOUT_S12  COMMENTS BELOW

    CALL PRINT_beam(B(1),MF,I=1)


  END SUBROUTINE THIN_EXAMPLE

  !  INTERFACE TRACK
  !     MODULE PROCEDURE TRACK_LAYOUT_12
  !     MODULE PROCEDURE TRACK_LAYOUT_S12
  !     MODULE PROCEDURE TRACK_THIN_T
  !  END INTERFACE

  !  SUBROUTINE TRACK_LAYOUT_12( R,B,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2 ) OF SMA_MULTIPARTICLE.F90

  ! Tracks through the thin lens structure R%T of the layout R if it exists.
  ! Several posibilities:
  !1) Pos1 and Pos2 are given :  tracks from thin lens position pos1 to thin lens position pos2
  !2) Thin lens points T1 and T2 are given: Same as above pos1 and pos2 are derived from T1 and T2
  !3) P1 and P2 are fibres: Tracks from the first thin lens of P1 to the first of P2. Results should agree
  !   with plain PTC.
  !4) P1, IN_P1 and P2 , IN_P2 are give: movies to IN_P1 thin lens in P1 and tracks to IN_P2 position in p2
  !  For example, if the input is (P1,1) and (P2,1) the results is same as item #3 (same as plain PTC).

  ! In all the above cases one can elect to give only the first input. Then it tracks one turn around as in
  ! plain PTC.
  ! interfaced as TRACK_LAYOUT_USING_THIN_S

  !  SUBROUTINE TRACK_LAYOUT_S12( R,B,K,S1,S2 )

  ! Tracks through the thin lens structure from position S1 to position S2 (defined as the S(3) variables
  ! of the thin lens.
  ! The final position is stored as in time tracking but is obviously the same for all the particles.
  ! interfaced as TRACK_LAYOUT_USING_THIN_S

  !  SUBROUTINE TRACK_THIN_T(B,DT,K)
  ! Tracks to full beam for a time DT
  ! All the particles are at different locations
  ! Notice that the layout is hidden: this is consistant with time tracking
  ! Magnets are not ontological objects




end module S_fitting
