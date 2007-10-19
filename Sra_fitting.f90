!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_fitting_new
  USE ptc_spin
  IMPLICIT NONE
  public

contains

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
    !   CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    CALL TRACK_probe_x(R,Y,+STATE,U,node1=IT1,node2=IT2)
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
    CALL TRACK_probe_x(R,Y,+STATE,U,node1=IT3,node2=IT4)
    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
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
    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    CALL TRACK_probe_x(R,Y,+STATE,U,node1=IT1,node2=IT2)

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
    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
    CALL TRACK_probe_x(R,Y,+STATE,U,node1=IT3,node2=IT4)
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

  subroutine lattice_linear_res_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
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
    integer :: neq=6, no=2,nt,j,it
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,tune(2),co(4)
    !    EPSF=.0001
    epsr=abs(epsf)

    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0))+only_4d0)-RADIATION0)

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
    co(1)=(y(1).SUB.'0010')
    co(2)=(y(1).SUB.'0001')
    co(3)=(y(2).SUB.'0010')
    co(4)=(y(2).SUB.'0001')
    id=y
    id=NORM%a1**(-1)*id*NORM%a1

    eq(1)=       ((NORM%dhdj%v(1)).par.'00000')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'00000')-targ(2)
    eq(3)=       (id%v(1).par.'0010') !-targ(3)
    eq(4)=       (id%v(1).par.'0001') !-targ(4)
    eq(5)=       (id%v(2).par.'0010') !-targ(5)
    eq(6)=       (id%v(2).par.'0001') !-targ(6)
    epsnow=abs(eq(1))+abs(eq(2))+abs(eq(3))+abs(eq(4))+abs(eq(5))+abs(eq(6))
    write(6,*) " closeness ",epsnow
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

  end subroutine lattice_linear_res_gmap


end module S_fitting_new
