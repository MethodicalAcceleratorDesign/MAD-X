!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and Frank Schmidt
! See file Sa_rotation_mis
module S_fitting
  USE MAD_LIKE
  IMPLICIT NONE
  PRIVATE lattice_fit_TUNE_L,lattice_fit_L   !,LAGRANGE
  PRIVATE THINLENS_L_B,THINLENS_L_L,THINLENS_L_2
  PRIVATE FIND_ORBIT_LAYOUT,FIND_ORBIT_M_LAYOUT,FIND_ENV_LAYOUT, FIND_ORBIT_LAYOUT_noda
  logical(lp), PRIVATE :: VERBOSE = .false.


  INTERFACE lattice_fit_TUNE
     ! LINKED
     MODULE PROCEDURE lattice_fit_TUNE_L
  END INTERFACE

  INTERFACE lattice_fit
     MODULE PROCEDURE lattice_fit_L
  END INTERFACE

  INTERFACE THINLENS
     !LINKED
     MODULE PROCEDURE THINLENS_L_L
     MODULE PROCEDURE THINLENS_L_2
  END INTERFACE

  INTERFACE FIND_ORBIT
     ! LINKED
     ! no use of TPSA
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
     ! returns linear matrix
     MODULE PROCEDURE FIND_ORBIT_M_LAYOUT
     !RETURN QUADRATIC YS
     MODULE PROCEDURE FIND_ENV_LAYOUT

  END INTERFACE

contains

  subroutine toggle_verbose
    implicit none
    verbose=.not.verbose
  end   subroutine toggle_verbose



  subroutine Lagrange(nfam,nfit,targ)
    implicit none
    logical(lp) :: doneitt=.true.
    INTEGER resultI,I,K,j
    INTEGER SCRATCHFILE
    INTEGER,PARAMETER::NEQ=8*4

    INTEGER, intent(in) :: nfam,nfit
    real(dp), intent(in),dimension(:)::TARG
    real(dp)  DELT ,DELTw
    real(dp) result  ,La(neq,neq),Lar(neq,neq),VLa(neq),vlb(neq)
    type (taylor) eq(NEQ),scrap,TEST,DF(NEQ)

    tpsafit(:)=zero
    VLa(:)=zero
    vlb(:)=zero
    LA(:,:)=zero
    Lar(:,:)=zero
    DELT=zero
    DELTw=zero

    call init(2,NFAM,doneitt)

    call alloc(eq,neq)
    call alloc(scrap)
    call alloc(test)
    call alloc(df,neq)


    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=1,NFIT
       call dainput(eq(i),scratchfile)
       eq(i)=eq(i)-targ(i)
    enddo
    SCRATCHFILE=CLOSEFILE



    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A120))'
    w_p%c(1)= " norms"
    call write_i
    DELT=zero
    do i=1,NFIT
       DELTw=eq(i).sub.'0'
       w_p=0
       w_p%nr=1
       w_p%fr='(1(1X,g20.14))'
       w_p%r(1)= DELTw
       call write_i
       DELT=DELT+abs(DELTw)
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A72))'
    write(w_p%c(1),'(1(1X,A6,1x,g20.14))')  " norm ",DELT
    call write_i


    DELT=one

    do i=1,NFIT
       eq(i)=eq(i)+(DELT-one)*(eq(i).sub.'0')
    enddo


    scrap=zero




    do i=1,neq
       vlb(i)=zero
       vla(i)=zero
       do k=1,neq
          la(i,k)= zero
       enddo
    enddo
    do i=NFAM+nfit+1,neq              ! initializing beyond useful part of array
       la(i,i)= one
    enddo

    do i=1,nfam
       la(i,i)= one      ! equal weight
    enddo


    do i=1,nfit
       result=eq(i).sub.'0'
       vla(i+nfam)=-result

       do k=1,NFAM
          scrap=eq(i).d.k

          result=       scrap.sub.'0'
          La(i+nfam,k)=result
          La(k,i+nfam)=result
       enddo

    enddo


    CALL MATINV(La,Lar,neq,neq,resultI)


    IF(resultI/=0) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='(1(1X,A72))'
       write(w_p%c(1),'(1(1X,A22,1x,i4))')  "MATRIX INVERSION FLAG ",resultI
       w_p%c(2)="PROBLEMS IN INVERSION MATINV"
       call write_e(1000)
    ENDIF

    do i=1,neq
       do k=1,neq
          vlb(i)=lar(i,k)*vla(k)+vlb(i)
       enddo
    enddo
    w_p=0
    w_p%nc=neq/8+2
    !    w_p%fc='(1(1X,A22,1x,i4))'
    w_p%c(1)='********    AJUSTMENTS   *********'
    do i=1,neq/8
       j=i*8-7
       write(w_p%c(i+1),'(8(2x,1pe12.5))') (vLb(k),k=j,j+7)
    enddo
    w_p%c(i+1) ='**********************************'
    call write_i



    do i=1,NFAM
       tpsafit(i)=vlb(i)   !vlb(i+1)
    enddo






    call kill(eq,neq)
    call kill(test)
    call kill(scrap)
    call kill(df,neq)


  end subroutine Lagrange

  ! linked
  subroutine lattice_fit_TUNE_L(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(2)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+ONLY_4D)-RADIATION0)
    CALL INIT(STATE,2,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,2,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,2)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(1)).par.'0000')
    eq(2)=       ((NORM%dhdj%v(2)).par.'0000')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=2
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,2)
    call Lagrange(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_TUNE_L


  subroutine lattice_fit_L(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,more
    INTEGER,parameter::nt=5
    TYPE(TAYLOR) EQ(nt)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer ipause, mypause

    SCRATCHFILE=90

    STATE=((DEFAULT+ONLY_4D+delta)-RADIATION0)
    CALL INIT(STATE,2,NP,BERZ,ND2,NPARA)
    ipause=mypause(888)
    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,nt)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(1)).par.'00000')
    eq(2)=       ((NORM%dhdj%v(2)).par.'00000')
    eq(3)=       (y(1).par.'00001')
    eq(4)=       ((NORM%dhdj%v(1)).par.'00001')
    eq(5)=       ((NORM%dhdj%v(2)).par.'00001')
    !  eq(6)=       (y(1).par.'00000')


    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,nt)
    call Lagrange(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1) = " More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100

    CALL ELP_TO_EL(R)

    CALL KILL_PARA(R)
  end subroutine lattice_fit_L

  SUBROUTINE  THINLENS_L_B(R,THI,kind_) ! A re-splitting routine
    IMPLICIT NONE
    INTEGER NTE
    TYPE(layout), intent(inout) :: R
    real(dp),intent(inout) :: THI
    real(dp) gg,RHOI,XL,QUAD
    INTEGER M1,M2,M3, MK1,MK2,MK3,limit(2),limit0(2)
    logical(lp) MANUAL,eject,doit
    integer,intent(in) :: kind_
    TYPE (fibre), POINTER :: C
    logical(lp) doneit
    nullify(C)

    CALL LINE_L(R,doneit)

    MANUAL=.FALSE.
    eject=.FALSE.
    if(kind_/=kind2.and.kind_/=kind7.and.kind_/=kind6.and.kind_/=kind10.and.kind_/=kind0.and.kind_/=kind16) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)=" THINLENS_L_B "
       call write_e(451)
    endif

    IF(THI<0) MANUAL=.TRUE.


    IF(MANUAL) THEN
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)="thi: thin lens factor (THI<0 TO STOP) "
       call write_I
       call read(thi)
       IF(THI<0) eject=.true.
    ENDIF

1001 CONTINUE

    limit(1)=3
    limit(2)=14
    if(kind_==kind6) then
       limit(1)=100000
       limit(2)=1000000
    endif
    limit0(1)=limit(1)
    limit0(2)=limit(2)

    M1=0
    M2=0
    M3=0
    MK1=0
    MK2=0
    MK3=0
    r%NTHIN=0

    C=>R%START
    do   WHILE(ASSOCIATED(C))

       if(kind_==kind0) then
          doit=C%MAG%KIND==kind2.or.C%MAG%KIND==kind10.or.C%MAG%KIND==kind16
       else
          doit=C%MAG%KIND==kind_
       endif
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
    w_p=0
    w_p%nc=6
    w_p%fc='(5(1X,a72,/),(1X,a72))'
    WRITE(w_p%c(1),'(1x,A13,1x,A24)') " Magnet type ",MYTYPE(kind_)
    WRITE(w_p%c(2),'(1x,A23,1x,I4)') "Present of thin lenses ",r%NTHIN
    WRITE(w_p%c(3),'(1x,A9,2(1x,I4))') "METHOD 2 ",M1,MK1
    WRITE(w_p%c(4),'(1x,A9,2(1x,I4))') "METHOD 4 ",M2,MK2
    WRITE(w_p%c(5),'(1x,A9,2(1x,I4))') "METHOD 6 ",M3,MK3
    WRITE(w_p%c(6),'(1x,A16,(1x,I4))') "number of KICKS ", MK1+MK2+MK3
    call write_I
    if(eject) then
       limit(1)=limit0(1)
       limit(2)=limit0(2)
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

       if(kind_==kind0) then
          doit=C%MAG%KIND==kind2.or.C%MAG%KIND==kind10.or.C%MAG%KIND==kind16
       else
          doit=C%MAG%KIND==kind_
       endif
       if(doit)  then


          xl=C%MAG%L
          RHOI=C%MAG%P%B0
          IF(C%MAG%P%NMUL>=2) THEN
             QUAD=SQRT(C%MAG%BN(2)**2+C%MAG%AN(2)**2)
          ELSE
             QUAD=zero
          ENDIF

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


    w_p=0
    w_p%nc=6
    w_p%fc='(5(1X,a72,/),(1X,a72))'
    WRITE(w_p%c(1),'(1x,A13,1x,A24)') " Magnet type ",MYTYPE(kind_)
    WRITE(w_p%c(2),'(1x,A23,1x,I4)') "Present of thin lenses ",r%NTHIN
    WRITE(w_p%c(3),'(1x,A9,2(1x,I4))') "METHOD 2 ",M1,MK1
    WRITE(w_p%c(4),'(1x,A9,2(1x,I4))') "METHOD 4 ",M2,MK2
    WRITE(w_p%c(5),'(1x,A9,2(1x,I4))') "METHOD 6 ",M3,MK3
    WRITE(w_p%c(6),'(1x,A16,(1x,I4))') "number of KICKS ", MK1+MK2+MK3
    call write_I


    IF(MANUAL) THEN
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1)= "thi: thin lens factor (THI<0 TO STOP) "
       call write_I
       call read(thi)
       IF(THI<0) THEN
          THI=R%THIN
          limit(1)=limit0(1)
          limit(2)=limit0(2)
          RETURN
       ELSE
          GOTO 1001
       ENDIF

    ENDIF

    limit(1)=limit0(1)
    limit(2)=limit0(2)

    CALL RING_L(R,doneit)

  END SUBROUTINE  THINLENS_L_B


  SUBROUTINE  THINLENS_L_L(R,THI,kind_) ! Interface to THINLENS_L_B
    IMPLICIT NONE
    TYPE(layout), intent(inout) :: R
    real(dp),intent(inout) :: THI
    integer,intent(in) :: kind_

    call THINLENS_L_B(R,THI,kind_)

  end SUBROUTINE THINLENS_L_L

  SUBROUTINE  THINLENS_L_2(R,THI) ! Interface to THINLENS_L_B
    IMPLICIT NONE
    TYPE(layout), intent(inout) :: R
    real(dp),intent(inout) :: THI
    integer :: kind_ !frs adding ::

    kind_=kind0
    call THINLENS_L_B(R,THI,kind_)

  end SUBROUTINE THINLENS_L_2

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
    logical(lp) APERTURE
    INTEGER TURNS0
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.
    TURNS0=1
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
             if(C%magp%kind==kind4) goto 101
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
             if(C%magp%kind==kind4) goto 101
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
       freq=0.d0
       i=1
       do while(i<=RING%n.and.freq==0.d0)
          c=>c%next
          if(associated(c%magp%freq)) then
             freq=c%magp%freq
          endif
          i=i+1
       enddo
       if(freq==0.d0) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(113)
       endif
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
       if(.not.check_stable) then
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
    c_%APERTURE_FLAG=APERTURE

  END SUBROUTINE FIND_ORBIT_LAYOUT


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
    INTEGER TURNS0
    TURNS0=1
    IF(PRESENT(TURNS)) TURNS0=TURNS

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
             if(C%magp%kind==kind4) goto 101
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
             if(C%magp%kind==kind4) goto 101
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
       freq=0.d0
       i=1
       do while(i<=RING%n.and.freq==0.d0)
          c=>c%next
          if(associated(c%magp%freq)) then
             freq=c%magp%freq
          endif
          i=i+1
       enddo
       if(freq==0.d0) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(113)
       endif
    endif




3   continue

    X=FIX

    DO I=1,TURNS0
       CALL TRACK(RING,X,LOC,STAT)
       if(.not.check_stable) then
          w_p=0
          w_p%nc=1
          w_p%fc='((1X,a72))'
          write(w_p%c(1),'(a30,i4)') " Lost in Fixed Point Searcher ",2
          call write_i

          return
       endif

    ENDDO



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

       do i=1,ND2
          MX(I,J)=(Y(i)-X(i))/eps
       enddo

    ENDDO


    SX=MX;
    DO I=1,6
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

  SUBROUTINE FIND_ORBIT_polymorph_noda(RING,FIX,LOC,STATE,eps,TURNS) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: RING
    real(dp) , intent(inOUT) :: FIX(6)
    INTEGER , intent(in) :: LOC
    INTEGER, OPTIONAL::TURNS
    real(dp) , optional,intent(in) :: eps
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat

    real(dp)  DIX(6),xdix,xdix0,tiny,freq
    real(dp) MX(6,6),sxi(6,6),SX(6,6)
    type(real_8) X(6),Y(6)
    integer NO1,ND2,I,IU,ITE,ier,j
    TYPE (fibre), POINTER :: C
    logical(lp) APERTURE
    INTEGER TURNS0
    TURNS0=1
    IF(PRESENT(TURNS)) TURNS0=TURNS

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
             if(C%magp%kind==kind4) goto 101
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
             if(C%magp%kind==kind4) goto 101
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
       freq=0.d0
       i=1
       do while(i<=RING%n.and.freq==0.d0)
          c=>c%next
          if(associated(c%magp%freq)) then
             freq=c%magp%freq
          endif
          i=i+1
       enddo
       if(freq==0.d0) then
          w_p=0
          w_p%nc=2
          w_p%fc='((1X,a72,/),(1X,a72))'
          w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          call write_E(113)
       endif
    endif


    call alloc(x)
    call alloc(y)
3   continue

    X=FIX

    DO I=1,TURNS0
       CALL TRACK(RING,X,LOC,STAT)
    ENDDO



    DO J=1,ND2

       Y=FIX
       Y(J)=FIX(J)+EPS
       DO I=1,TURNS0
          CALL TRACK(RING,Y,LOC,STAT)
       ENDDO
       do i=1,ND2
          MX(I,J)=(Y(i)-X(i))/eps
       enddo

    ENDDO


    SX=MX;
    DO I=1,6
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

    do i=1,6
       x(i)=zero
    enddo
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
    call kill(x)
    call kill(y)

  END SUBROUTINE FIND_ORBIT_polymorph_noda




  SUBROUTINE FIND_ORBIT_M_LAYOUT(RING,X,LOC,STATE,TURNS) ! Finds orbit and linear Map with TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),INTENT(INOUT):: RING
    TYPE(REAL_8),INTENT(INOUT):: X(6)
    real(dp)  FIX(6),DIX(6),xdix,xdix0,tiny
    TYPE(DAMAP) MX,SX,SXI,IS
    integer NO1,ND2,I,IU,LOC,ITE,npara
    TYPE(INTERNAL_STATE) , optional, intent(in) :: STATE      !_in
    INTEGER,OPTIONAL :: TURNS
    TYPE(INTERNAL_STATE)  STAT
    !   TYPE(INTERNAL_STATE) state
    TYPE(fibre), POINTER :: C
    logical(lp) APERTURE
    INTEGER TURNS0
    TURNS0=1
    IF(PRESENT(TURNS)) TURNS0=TURNS
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.
    Nullify(C);

    if(.not.ring%closed) then
       w_p=0
       w_p%nc=1
       w_p%fc='((1X,a72))'
       w_p%c(1) = " This line is not ring : FIND_ORBIT_M_LAYOUT "
       call write_e(100)
    endif

    if(present(state)) then
       stat=state
    else
       stat=default
    endif


    dix(:)=zero
    tiny=c_1d_40
    NO1=1
    xdix0=c_1d4*DEPS_tracking
    IF(stat%NOCAVITY) THEN
       !    ND1=2
       !    state=STATE_in+only_4d
    ELSE
       !    state=STATE_in
       !   ND1=3
       C=>RING%START
       do i=1,RING%n
          if(C%magp%kind==kind4) goto 101
          C=>C%NEXT
       enddo
       w_p=0
       w_p%nc=2
       w_p%fc='((1X,a72))'
       w_p%c(1) = " No Cavity in the Line "
       w_p%c(2) = " FIND_ORBIT_LAYOUT will crash "
       call write_e(111)
    ENDIF
101 continue
    !    ND2=2*ND1
    CALL INIT(stat,NO1,0,BERZ,ND2,NPARA)


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
    ENDDO

    IS=1
    MX=X
    SX=MX-IS
    DIX=SX
    if(nd2==6.and.stat%NOCAVITY)  then
       SX%v(5)=one.mono.'000010'
       SX%v(6)=one.mono.'000001'
    endif
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
    CALL KILL(MX)
    CALL KILL(SX)
    CALL KILL(SXI)
    CALL KILL(IS)
    c_%APERTURE_FLAG=aperture

  END SUBROUTINE FIND_ORBIT_M_LAYOUT

  SUBROUTINE FIND_ENV_LAYOUT(RING,YS,FIX,LOC,STATE) ! Finds envelope with TPSA in State or compatible state
    IMPLICIT NONE
    logical(lp) :: doneitt=.true.
    TYPE(layout),INTENT(INOUT):: RING
    real(dp)  FIX(6),DIX(6),xdix,xdix0,mat(6,6),flu(6,6)
    TYPE(REAL_8) X(6)
    TYPE(ENV_8),INTENT(INOUT)::YS(6)
    TYPE(DAMAP) MX,SX,SXI,IS
    type(beamenvelope) env
    integer NO1,ND2,ND1 ,I,IU,LOC,ITE ,J
    INTEGER  JJ(LNV)
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
       if(C%magp%kind==kind4) goto 101
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


3   continue
    X=ND2
    DO I=1,6
       X(I)=FIX(I)
    ENDDO



    CALL TRACK(RING,X,LOC,SSS)

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
       if(xdix.ge.xdix0) then
          ite=0
       else
          ite=1
          xdix0=xdix
       endif
    endif
    if(ite.eq.1) then
       GOTO 3
    endif

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

!!!!!   Fernando !!!!!!!!!!
  subroutine lattice_fit_CHROM_ONLY(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(2)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+ONLY_4D+DELTA)-RADIATION0)
    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,2)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(1)).par.'00001')
    eq(2)=       ((NORM%dhdj%v(2)).par.'00001')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=2
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,2)
    call Lagrange(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_CHROM_ONLY

  subroutine lattice_fit_CHROM_alpha2(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(3)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+nocavity0)-RADIATION0)
    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,3)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(1)).par.'000010')
    eq(2)=       ((NORM%dhdj%v(2)).par.'000010')
    eq(3)=       ((NORM%dhdj%v(3)).par.'000020')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=3
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,3)
    call Lagrange(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_CHROM_alpha2

  subroutine lattice_fit_alpha23(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(2)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+nocavity0)-RADIATION0)
    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,5,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,3)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(3)).par.'000020')
    eq(2)=       ((NORM%dhdj%v(3)).par.'000030')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=2
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,2)
    call Lagrange_t(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_alpha23

  subroutine lattice_fit_CHROM_alpha23(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(4)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+nocavity0)-RADIATION0)
    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,4,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,4)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(3)).par.'000020')
    eq(2)=       ((NORM%dhdj%v(3)).par.'000030')
    eq(3)=       ((NORM%dhdj%v(1)).par.'000010')
    eq(4)=       ((NORM%dhdj%v(2)).par.'000010')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=4
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,4)
    call Lagrange_t(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_CHROM_alpha23

  subroutine lattice_fit_alpha234(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(3)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+nocavity0)-RADIATION0)
    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,5,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,3)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(3)).par.'000020')
    eq(2)=       ((NORM%dhdj%v(3)).par.'000030')
    eq(3)=       ((NORM%dhdj%v(3)).par.'000040')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=3
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,3)
    call Lagrange_t(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_alpha234

  subroutine latticefitchromxoryalpha234(R,plane,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more,plane
    TYPE(TAYLOR) EQ(4)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+nocavity0)-RADIATION0)
    CALL INIT(STATE,3,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,5,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,4)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(3)).par.'000020')
    eq(2)=       ((NORM%dhdj%v(3)).par.'000030')
    eq(3)=       ((NORM%dhdj%v(3)).par.'000040')
    eq(4)=       ((NORM%dhdj%v(plane)).par.'000010')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=4
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,4)
    call Lagrange_t(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine latticefitchromxoryalpha234

  subroutine lattice_fit_TUNE_alpha(R,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE,ND2,NPARA,nt,more
    TYPE(TAYLOR) EQ(3)

    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM

    SCRATCHFILE=90

    STATE=((DEFAULT+nocavity0)-RADIATION0)
    CALL INIT(STATE,2,NP,BERZ,ND2,NPARA)

    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       POLY(i)%NPARA=NPARA

       R=POLY(i)
    ENDDO



    CLOSED(:)=zero
100 continue

    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_7)

    CALL INIT(STATE,2,NP,BERZ,ND2,NPARA)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ,3)

    Y=NPARA
    Y=CLOSED
    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    write(w_p%c(1),'(8(2x,g20.14))') NORM%TUNE(1), NORM%TUNE(2)
    call write_i

    eq(1)=       ((NORM%dhdj%v(1)).par.'000000')
    eq(2)=       ((NORM%dhdj%v(2)).par.'000000')
    eq(3)=       ((NORM%dhdj%v(3)).par.'000010')

    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    nt=3
    do i=1,nt
       call shiftda(eq(i),eq(i),NPARA)
       call daprint(eq(i),scratchfile)
    enddo
    SCRATCHFILE=CLOSEFILE
    CALL KILL(NORM)
    CALL KILL(Y)
    CALL KILL(EQ,3)
    call Lagrange(np,nt,targ)

    SET_TPSAFIT=.true.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    w_p=0
    w_p%nc=1
    w_p%fc='((1X,A72))'
    w_p%c(1)=" More =>  yes=1"
    call write_i
    call read(more)
    if(more==1) goto 100
    CALL ELP_TO_EL(R)
    CALL KILL_PARA(R)
  end subroutine lattice_fit_TUNE_alpha


  subroutine Lagrange_t(nfam,nfit,targ)
    implicit none
    logical(lp) :: doneitt=.true.
    INTEGER resultI,I,K,j
    INTEGER SCRATCHFILE
    INTEGER,PARAMETER::NEQ=8*4

    INTEGER, intent(in) :: nfam,nfit
    real(dp), intent(in),dimension(:)::TARG
    real(dp)  DELT ,DELTw,eq0
    real(dp) result ,La(neq,neq),Lar(neq,neq),VLa(neq),vlb(neq)
    type (taylor) eq(NEQ),scrap,TEST,DF(NEQ)

    tpsafit(:)=zero
    VLa(:)=zero
    vlb(:)=zero
    LA(:,:)=zero
    Lar(:,:)=zero
    DELT=zero
    DELTw=zero

    call init(2,NFAM,doneitt)

    call alloc(eq,neq)
    call alloc(scrap)
    call alloc(test)
    call alloc(df,neq)


    SCRATCHFILE=NEWFILE
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=1,NFIT
       call dainput(eq(i),scratchfile)
       eq0=abs(eq(i))
       write(6,*) eq0,targ(i)
       eq(i)=eq(i)-targ(i)

    enddo
    SCRATCHFILE=CLOSEFILE



    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A120))'
    w_p%c(1)= " norms"
    call write_i
    DELT=zero
    do i=1,NFIT
       DELTw=eq(i).sub.'0'
       w_p=0
       w_p%nr=1
       w_p%fr='(1(1X,g20.14))'
       w_p%r(1)= DELTw
       call write_i
       DELT=DELT+abs(DELTw)
    enddo
    w_p=0
    w_p%nc=1
    w_p%fc='(1(1X,A72))'
    write(w_p%c(1),'(1(1X,A6,1x,g20.14))')  " norm ",DELT
    call write_i


    DELT=one
    write(6,*) "delt "
    read(5,*) delt

    do i=1,NFIT
       eq(i)=eq(i)+(DELT-one)*(eq(i).sub.'0')
    enddo


    scrap=zero




    do i=1,neq
       vlb(i)=zero
       vla(i)=zero
       do k=1,neq
          la(i,k)= zero
       enddo
    enddo
    do i=NFAM+nfit+1,neq              ! initializing beyond useful part of array
       la(i,i)= one
    enddo

    do i=1,nfam
       la(i,i)= one      ! equal weight
    enddo


    do i=1,nfit
       result=eq(i).sub.'0'
       vla(i+nfam)=-result

       do k=1,NFAM
          scrap=eq(i).d.k

          result=       scrap.sub.'0'
          La(i+nfam,k)=result
          La(k,i+nfam)=result
       enddo

    enddo

    CALL MATINV(La,Lar,neq,neq,resultI)

    IF(resultI/=0) THEN
       w_p=0
       w_p%nc=2
       w_p%fc='(1(1X,A72))'
       write(w_p%c(1),'(1(1X,A22,1x,i4))')  "MATRIX INVERSION FLAG ",resultI
       w_p%c(2)="PROBLEMS IN INVERSION MATINV"
       call write_e(1000)
    ENDIF

    do i=1,neq
       do k=1,neq
          vlb(i)=lar(i,k)*vla(k)+vlb(i)
       enddo
    enddo
    w_p=0
    w_p%nc=neq/8+2
    !    w_p%fc='(1(1X,A22,1x,i4))'
    w_p%c(1)='********    AJUSTMENTS   *********'
    do i=1,neq/8
       j=i*8-7
       write(w_p%c(i+1),'(8(2x,1pe12.5))') (vLb(k),k=j,j+7)
    enddo
    w_p%c(i+1) ='**********************************'
    call write_i



    do i=1,NFAM
       tpsafit(i)=vlb(i)   !vlb(i+1)
    enddo






    call kill(eq,neq)
    call kill(test)
    call kill(scrap)
    call kill(df,neq)


  end subroutine Lagrange_t


end module S_fitting
