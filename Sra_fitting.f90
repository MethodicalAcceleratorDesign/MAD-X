!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_fitting_new
  USE ptc_spin
  IMPLICIT NONE
  public
  integer m_turn
  integer :: with_c=1
  TYPE fibre_monitor_data
     type(fibre), pointer :: p    ! fibre location
     integer, pointer ::  turn,kind  ! kind=1 x, kind = 2 y
     real(dp), pointer :: r(:,:)  ! store fake experiment from alex_track_monitors
     real(dp), pointer :: xf(:,:)  ! real data put here
     real(dp), pointer :: mom(:,:)
     real(dp), pointer :: A(:,:)
     real(dp), pointer :: At(:,:)
  END TYPE fibre_monitor_data

  TYPE(fibre_monitor_data), allocatable :: monitors(:)

contains

  subroutine special_alex_main_ring(r,n_name,targ,sc)
    implicit none
    TYPE(layout), target, intent(inout):: R
    integer  i1,i2,it1,it2,it3,it4
    INTEGER I,N,NU,N2,NP2,mf,nt,NP,J,n_name
    type(fibre), pointer :: p1,p2
    TYPE(POL_BLOCK) QC(11)
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    REAL(DP) X(6),targ(2),xx
    TYPE(INTERNAL_STATE) STATE
    LOGICAL(LP) U
    type(normalform) nf
    type(gmap) g
    type(TAYLOR) T,eq(5)

    REAL(DP) ALX,ALY,TA(2),sc

    !targ(1)=22.43d0
    !targ(2)=20.82d0

    N=10
    NU=11
    if(.not.associated(r%t)) then
       write(6,*) " thin lens lattice not made "
       CALL MAKE_node_LAYOUT(r)
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
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT1,node2=IT2)
    U=.NOT.CHECK_STABLE
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
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT3,node2=IT4)
    U=.NOT.CHECK_STABLE

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

  subroutine special_alex_main_ring_auto(r,n_name,targ,sc,epsa)
    implicit none
    TYPE(layout), target, intent(inout):: R
    integer  i1,i2,it1,it2,it3,it4
    INTEGER I,N,NU,N2,NP2,mf,nt,NP,J,n_name
    type(fibre), pointer :: p1,p2
    TYPE(POL_BLOCK) QC(11)
    TYPE(REAL_8) Y(6)
    TYPE(DAMAP) ID
    REAL(DP) X(6),targ(2),xx,epsa,nxx
    TYPE(INTERNAL_STATE) STATE
    LOGICAL(LP) U
    type(normalform) nf
    type(gmap) g
    type(TAYLOR) T,eq(5)

    REAL(DP) ALX,ALY,TA(2),sc

    !targ(1)=22.43d0
    !targ(2)=20.82d0

    N=10
    NU=11
    if(.not.associated(r%t)) then
       write(6,*) " thin lens lattice not made "
       CALL MAKE_node_LAYOUT(r)
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
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT1,node2=IT2)
    U=.NOT.CHECK_STABLE

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
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT3,node2=IT4)
    U=.NOT.CHECK_STABLE

    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
    nf=y
    eq(1)=(eq(1)*8 + (1.d0+nf%dhdj%v(1)))*3-targ(1)
    eq(2)=(eq(2)*8 + (1.d0+nf%dhdj%v(2)))*3-targ(2)
    eq(3)=eq(3)-ta(1)

    eq(4)=eq(4)-(nf%A_T%V(1).par.'1000')**2-(nf%A_T%V(1).par.'0100')**2
    eq(5)=eq(5)-(nf%A_T%V(3).par.'0010')**2-(nf%A_T%V(3).par.'0001')**2
    !call print(eq(4),6)
    !call print(eq(5),6)
    nxx=zero
    do i=1,5
       eq(i)=eq(i)<=4
       xx=eq(i)
       write(6,*) i,xx
       nxx=nxx+abs(xx)
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
    mf=1
    if(  nxx<=epsa) mf=0
    !    WRITE(6,*) " MORE "
    !    READ(5,*) MF
    IF(MF==1) GOTO 111


    CALL KILL_PARA(R)

    CALL ELP_TO_EL(R)

  end   subroutine special_alex_main_ring_auto


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
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT1,node2=IT2)
    U=.NOT.CHECK_STABLE

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
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT3,node2=IT4)
    U=.NOT.CHECK_STABLE

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
    INTEGER I,SCRATCHFILE
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


  subroutine special_alex_main_ring_removal(main)
    IMPLICIT NONE

!!!!!!!!!
    type(layout),pointer :: main
    type(fibre),pointer :: p,p1,p2,pp
    real(dp) an,bn
    logical doit
    integer ij,mf,i,k,jk

    call kanalnummer(mf,"diag.txt")


    if(associated(main%t)) call kill(main%t)

    ij=0
    jk=0
    p=>main%start
    do i=1,main%n
       p1=>p%next
       p2=>p1%next
       doit=p%mag%name(1:len_trim(p%mag%name)-1)==p1%mag%name(1:len_trim(p1%mag%name)-1)
       doit=doit.and.(p%mag%name(1:len_trim(p%mag%name)-1)==p2%mag%name(1:len_trim(p1%mag%name)-1))
       if(doit) then
          ij=ij+1
          do k=p%mag%p%nmul,1,-1
             an=p%mag%an(k)/p1%mag%l*2.d0
             bn=p%mag%bn(k)/p1%mag%l*2.d0
             call add(p1,-k,1,an)
             call add(p1,k,1,bn)
             !   call add(p,-k,0,0.d0)
             !   call add(p,k,0,0.d0)
             !   call add(p2,-k,0,0.d0)
             !   call add(p2,k,0,0.d0)
          enddo
          ! el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle
          if(p%mag%kind==kind3) then
             bn=p%mag%k3%thin_h_angle/p1%mag%l*2.d0
             call add(p1,1,1,bn)
             an=p%mag%k3%thin_v_angle/p1%mag%l*2.d0
             call add(p1,-1,1,an)
             doit=(p%mag%k3%thin_h_foc/=zero.or.p%mag%k3%thin_v_foc/=zero)  !.or.p%mag%k3%thin_v_angle/=zero)
             if(doit) then
                write(6,*) " cannot handle the stuff in kind3 related to fake thinlens tracking of MAD8 "
                stop 8
             endif
          endif

          write(mf,*) p%mag%name,p1%mag%name,p2%mag%name
          !  goto 200
          pp=>p%previous
          pp%next=>p1
          p1%previous=>pp
          pp=>p2%next
          pp%previous=>p1
          p1%next=>pp

          call super_kill(p)
          call super_kill(p2)
          p=>p1
          !   200 continue
          !   p=>p2
       else
          if(p%mag%kind==kind3) then
             jk=jk+1
             bn=p%mag%k3%thin_h_angle
             call add(p,1,1,bn)
             an=p%mag%k3%thin_v_angle
             call add(p,-1,1,an)
             p%mag%k3%thin_v_angle=zero
             p%mag%k3%thin_h_angle=zero
             p%magp%k3%thin_v_angle=zero
             p%magp%k3%thin_h_angle=zero
             doit=(p%mag%k3%thin_h_foc/=zero.or.p%mag%k3%thin_v_foc/=zero)
             if(doit) then
                write(6,*) " cannot handle the stuff in kind3 related to fake thinlens tracking of MAD8 "
                stop 9
             endif
          endif
       endif
       p=>p%next
       if(associated(p,main%start)) exit


    enddo

    write(mf,*) ij, " magnet sandwiches done"
    write(mf,*) jk, " single magnets done"

    ij=0
    p=>main%start
    do i=1,main%n
       ij=ij+1
       p%pos=i
       if(p%mag%l==zero.and.p%mag%kind>kind1) then
          write(mf,*) i,p%mag%name,p%mag%kind,p%next%mag%name
       endif
       p=>p%next
       if(associated(p,main%start)) exit
    enddo
    write(6,*) ij,main%n
    main%n=ij

    close(mf)
  end subroutine special_alex_main_ring_removal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!  real stuff   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  alex_track_monitors(r,x,state) ! one particle tracking m_turn
    implicit none
    type(layout), target :: r
    integer nturn,i,j,nm
    real(dp) x(6),y(6)
    type(fibre),pointer::p1
    type(internal_state) state

    p1=>r%start

    nturn=SIZE (monitors(1)%r,dim=2)
    nm=size(monitors)
    y=x
    CALL FIND_ORBIT(R,y,1,STATE,c_1d_5)
    write(6,*)" closed orbit "
    write(6,*) y(1:3)
    write(6,*) y(4:6)
    x(1:4)=x(1:4)+y(1:4)
    do i=1,nturn
       call track(x,state,p1,monitors(1)%p)

       call track(y,state,p1,monitors(1)%p)
       monitors(1)%r(1:4,i)=x(1:4)-y(1:4)

       do j=1,nm-1
          call track(x,state,monitors(j)%p,monitors(j+1)%p)
          call track(y,state,monitors(j)%p,monitors(j+1)%p)
          monitors(j+1)%r(1:4,i)=x(1:4)-y(1:4)

       enddo
       call track(x,state,monitors(nm)%p,p1)
       call track(y,state,monitors(nm)%p,p1)

    enddo





  end subroutine  alex_track_monitors


  subroutine  alex_mom_real_monitors(ring,jm,jn,x,state_in)
    implicit none
    integer ipause, mypause
    integer jm,nm,i,jinv(4),i1,i11,i2,i22,jn(4)
    type(layout), target :: ring
    real(dp), pointer :: mom(:,:),r(:,:),a(:,:),at(:,:)
    real(dp)ar,x(6),z(4),zz(lnv),zx(4),z0(4),zf(4)
    type(damap) id,m11,m2,m22,m1,ex,map1,fin
    type(pbfield) h
    type(normalform) norm
    type(fibre),pointer::p1,p11,p2,p22
    type(real_8) y(6)
    type(internal_state) state,state_in
    type(gmap) g,gi
!!! will produce the real data at jm

    z0=zero

    nm=size(monitors)

    state=state_in+only_4d0

    CALL INIT(STATE,1,0)

    call alloc(y)
    call alloc(id,m11,m2,m22,m1,ex,fin)
    call alloc(map1)
    call alloc(g,gi)

    CALL FIND_ORBIT(ring,x,1,STATE,c_1d_5)
    write(6,*) monitors(jm)%kind
    jn(1)=jm
    if(monitors(jm)%kind==1) then

       p1=>monitors(jm)%p
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1) then
             p11=>monitors(mod_n(i,nm))%p
             jn(2)=mod_n(i,nm)
             i11=i
             exit
          endif
       enddo
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2) then
             p2=>monitors(mod_n(i,nm))%p
             jn(3)=mod_n(i,nm)
             i2=i
             exit
          endif
       enddo
       do i=i2+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2) then
             p22=>monitors(mod_n(i,nm))%p
             jn(4)=mod_n(i,nm)
             i22=i
             exit
          endif
       enddo
       write(6,*)
       write(6,*) p1%mag%name,p11%mag%name
       write(6,*) p2%mag%name,p22%mag%name
       write(6,*) p1%pos,p11%pos
       write(6,*) p2%pos,p22%pos

       ! track from ip to p1 : jm x monitor

       call track(x,state,ring%start,p1)
       !
       ! x is closed orbit at p1



       id=1

       y=x+id
       call track(y,state,p1,p11)
       m11=y
       ! m11 is map from p1 to p11  (second x monitor)
       m11=z0  ! removed zeroth order

       y=x+id
       call track(y,state,p1,p2)
       m2=y
       m2=z0
       ! m2 is map from p1 to p2  (first y monitor)


       y=x+id
       call track(y,state,p1,p22)
       m22=y
       m22=z0
       ! m22 is map from p1 to p22  (second y monitor)

!!!!!!!!!!!!   testing and computing m11 !!!!!!!!!!!!!!!!!
       fin=1
       g=1

       zF=0.d0

       z=0.d0
       z(1)=0.00001d0
       z(2)=0.00002d0
       z(3)=0.00003d0
       z(4)=0.00004d0
       zz=0.d0
       zf(1)=z(1)
       zz=m11*z
       write(6,*) "testing tracking from jm to p1"
       ! testing tracking from jm to p1
       WRITE(6,'(4(1x,E15.8))')z
       WRITE(6,'(4(1x,E15.8))')zz(1:4)

       ipause=mypause(777)

       jinv=0
       jinv(1)=1
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=m11*ex
       m11=m11**jinv
       !m11=m11**(-1)
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=ex*m11*ex

       ! m11 is now the map from x,x11,y,py

       fin%v(2)=m11%v(2)   !   px_0(x_0,x_p11,y_0,py_0)
       g%v(2)=m11%v(2)     !   px_0(x_0,x_p11,y_0,py_0)
       zf(2)=zz(1)

       zx(1)=z(1)
       zx(2)=zz(1)
       zx(3)=z(3)
       zx(4)=z(4)

       ex=0
       ex=zx



       map1=m11.o.ex

       zz(1:4)=map1

       WRITE(6,'(4(1x,E15.8))')zz(1:4)
       call print(map1%v(2),6)
       ipause=mypause(778)
!!!!!!!!!!!!   testing and computing m22  !!!!!!!!!!!!!!!!!
       z=0.d0
       z(1)=0.00001d0
       z(2)=0.00002d0
       z(3)=0.00003d0
       z(4)=0.00004d0
       zz=0.d0

       zz=m22*z
       write(6,*) "testing tracking from jm to p22"
       ! testing tracking from jm to p2

       WRITE(6,'(4(1x,E15.8))')z
       WRITE(6,'(4(1x,E15.8))')zz(1:4)

       ipause=mypause(777)

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=m22*ex
       m22=m22**jinv
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=ex*m22*ex
       ! m22 is now the map from x,px,y,y22


       fin%v(4)=m22%v(4)    !   py_0(x_0,px_0,y_0,y_p22)
       g%v(4)=m22%v(4)

       zf(4)=zz(3)

       zx(1)=z(1)
       zx(2)=z(2)
       zx(3)=z(3)
       zx(4)=zz(3)

       ex=0
       ex=zx



       map1=m22.o.ex

       zz(1:4)=map1
       WRITE(6,'(4(1x,E15.8))')zz(1:4)
       call print(map1%v(4),6)
       ipause=mypause(778)
!!!!!!!!!!!!   testing and computing m2  !!!!!!!!!!!!!!!!!
       z=0.d0
       z(1)=0.00001d0
       z(2)=0.00002d0
       z(3)=0.00003d0
       z(4)=0.00004d0
       zz=0.d0

       zz=m2*z

       WRITE(6,'(4(1x,E15.8))')z
       WRITE(6,'(4(1x,E15.8))')zz(1:4)

       ipause=mypause(777)

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m2=m2*ex
       m2=m2**jinv
       !m11=m11**(-1)
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m2=ex*m2*ex

       ! m2 is now the map from x,px,y2,py

       fin%v(3)=m2%v(3)      !   y_0(x_0,px_0,y_p2,py_0)
       g%v(3)=m2%v(3)
       zf(3)=zz(3)

       zx(1)=z(1)   !x_1
       zx(2)=z(2)    ! x_0
       zx(3)=zz(3)    !
       zx(4)=z(4)

       ex=0
       ex=zx



       map1=m2.o.ex

       zz(1:4)=map1
       WRITE(6,'(4(1x,E15.8))')zz(1:4)
       call print(map1%v(3),6)
       ipause=mypause(778)

       call print(fin,6)



       gi=g



    else  !!!! if(monitors(jm)%kind==2)

       write(6,*) " crotte "
       ipause=mypause(0)
       p2=>monitors(jm)%p
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2) then
             p22=>monitors(mod_n(i,nm))%p
             jn(2)=mod_n(i,nm)
             i22=i
             exit
          endif
       enddo
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1) then
             p1=>monitors(mod_n(i,nm))%p
             jn(3)=mod_n(i,nm)
             i1=i
             exit
          endif
       enddo
       do i=i22+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1) then
             p11=>monitors(mod_n(i,nm))%p
             jn(4)=mod_n(i,nm)
             i11=i
             exit
          endif
       enddo
       write(6,*)
       write(6,*) p1%mag%name,p11%mag%name
       write(6,*) p2%mag%name,p22%mag%name
       write(6,*) p1%pos,p11%pos
       write(6,*) p2%pos,p22%pos

       ! track from ip to p1 : jm x monitor

       call track(x,state,ring%start,p2)
       !
       ! x is closed orbit at p2



       id=1

       y=x+id
       call track(y,state,p2,p22)
       m22=y
       ! m22 is map from p2 to p22  (second y monitor)
       m22=z0  ! removed zeroth order

       y=x+id
       call track(y,state,p2,p1)
       m1=y
       m1=z0
       ! m1 is map from p2 to p1  (first x monitor)


       y=x+id
       call track(y,state,p2,p11)
       m11=y
       m11=z0
       ! m11 is map from p2 to p11  (second x monitor)

!!!!!!!!!!!!   testing and computing m22 !!!!!!!!!!!!!!!!!  asdad
       fin=1
       g=1

       zF=0.d0

       z=0.d0
       z(1)=0.00001d0
       z(2)=0.00002d0
       z(3)=0.00003d0
       z(4)=0.00004d0
       zz=0.d0
       zf(3)=z(3)
       zz=m22*z
       write(6,*) "testing tracking from jm to p22"
       ! testing tracking from jm to p1
       WRITE(6,'(4(1x,E15.8))')z
       WRITE(6,'(4(1x,E15.8))')zz(1:4)

       ipause=mypause(777)

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=m22*ex
       m22=m22**jinv
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=ex*m22*ex

       ! m22 is now the map from x,px,y,y22

       fin%v(4)=m22%v(4)   !   py_0(x_0,p_x0,y_0,y_p22)
       g%v(4)=m22%v(4)     !   py_0(x_0,p_x0,y_0,y_p22)
       zf(4)=zz(3)

       zx(1)=z(1)
       zx(2)=z(2)
       zx(3)=z(3)
       zx(4)=zz(3)

       ex=0
       ex=zx



       map1=m22.o.ex

       zz(1:4)=map1

       WRITE(6,'(4(1x,E15.8))')zz(1:4)
       call print(map1%v(4),6)
       ipause=mypause(666)
!!!!!!!!!!!!   testing and computing m1  !!!!!!!!!!!!!!!!!
       z=0.d0
       z(1)=0.00001d0
       z(2)=0.00002d0
       z(3)=0.00003d0
       z(4)=0.00004d0
       zz=0.d0

       zz=m1*z
       write(6,*) "testing tracking from jm to p1"
       ! testing tracking from jm to p1

       WRITE(6,'(4(1x,E15.8))')z
       WRITE(6,'(4(1x,E15.8))')zz(1:4)

       ipause=mypause(777)

       jinv=0
       jinv(1)=1
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m1=m1*ex
       m1=m1**jinv
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m1=ex*m1*ex
       ! m1 is now the map from x,x1,y,py


       fin%v(2)=m1%v(2)    !   px_0(x_0,x1,y_0,py0)
       g%v(2)=m1%v(2)

       zf(2)=zz(1)

       zx(1)=z(1)
       zx(2)=zz(1)
       zx(3)=z(3)
       zx(4)=z(4)

       ex=0
       ex=zx



       map1=m1.o.ex

       zz(1:4)=map1
       WRITE(6,'(4(1x,E15.8))')zz(1:4)
       call print(map1%v(2),6)
       ipause=mypause(667)
!!!!!!!!!!!!   testing and computing m11  !!!!!!!!!!!!!!!!!
       z=0.d0
       z(1)=0.00001d0
       z(2)=0.00002d0
       z(3)=0.00003d0
       z(4)=0.00004d0
       zz=0.d0

       zz=m11*z

       WRITE(6,'(4(1x,E15.8))')z
       WRITE(6,'(4(1x,E15.8))')zz(1:4)

       ipause=mypause(777)

       jinv=0
       jinv(1)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=m11*ex
       m11=m11**jinv
       !m11=m11**(-1)
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=ex*m11*ex

       ! m11 is now the map from x2,px,y,py

       fin%v(1)=m11%v(1)      !   x_0(x_0,px_0,y_p2,py_0)
       g%v(1)=m11%v(1)
       zf(1)=zz(1)

       zx(3)=z(3)   !x_1
       zx(2)=z(2)    ! x_0
       zx(1)=zz(1)    !
       zx(4)=z(4)

       ex=0
       ex=zx



       map1=m11.o.ex

       zz(1:4)=map1
       WRITE(6,'(4(1x,E15.8))')zz(1:4)
       call print(map1%v(1),6)
       ipause=mypause(668)

       call print(fin,6)



       gi=g


    endif !!!! if(monitors(jm)%kind==2)

    call invert_monitors(monitors(jm)%kind,zf,gi)
    write(6,*) zf(1:4)
    ipause=mypause(555)



    call kill(y)
    call kill(id,m11,m2,m22,m1,ex,fin)
    call kill(map1)
    call kill(g,gi)



  end subroutine  alex_mom_real_monitors

  subroutine invert_monitors(kind,zf0,g)
    implicit none
    type(damap) ex
    type(gmap) gi,g
    real(dp) zf(4),zx(6),zf0(4)
    integer kind,i

    zf=zf0

    call alloc(ex)
    call alloc(gi)
    ! zf contains measuments xf(x_1,x_11,y_2,y_22)
    gi=g

    if(kind==1) then
       !write(6,*) zf(1:4)
       zx=zero
       zx(1)=zf(1)
       zx(2)=zf(2)   !!!
       ex=0
       ex=zx
       ex%v(3)=ex%v(3)+(1.d0.mono.3)
       ex%v(4)=ex%v(4)+(1.d0.mono.4)
       gi%v(2)=(gi%v(2).o.ex)-(1.d0.mono.2)
       zx=zero
       zx(1)=zf(1)
       zx(3)=zf(3)   !!!
       ex=0
       ex=zx
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       ex%v(4)=ex%v(4)+(1.d0.mono.4)

       gi%v(3)=(gi%v(3).o.ex)-(1.d0.mono.3)
       zx=zero
       zx(1)=zf(1)
       zx(4)=zf(4)   !!!
       ex=0
       ex=zx
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       ex%v(3)=ex%v(3)+(1.d0.mono.3)
       gi%v(4)=(gi%v(4).o.ex)-(1.d0.mono.4)
       !call print(gi,6)
       gi%v(1)=gi%v(1)-zf(1)

    else
       !write(6,*) zf(1:4)

       zx=zero
       zx(3)=zf(3)
       zx(4)=zf(4)   !!!
       ex=0
       ex=zx
       ex%v(1)=ex%v(1)+(1.d0.mono.1)
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       gi%v(4)=(gi%v(4).o.ex)-(1.d0.mono.4)
       zx=zero
       zx(1)=zf(1)
       zx(3)=zf(3)   !!!
       ex=0
       ex=zx
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       ex%v(4)=ex%v(4)+(1.d0.mono.4)

       gi%v(1)=(gi%v(1).o.ex)-(1.d0.mono.1)
       zx=zero
       zx(3)=zf(3)
       zx(2)=zf(2)   !!!
       ex=0
       ex=zx
       ex%v(4)=ex%v(4)+(1.d0.mono.4)
       ex%v(1)=ex%v(1)+(1.d0.mono.1)
       gi%v(2)=(gi%v(2).o.ex)-(1.d0.mono.2)
       !call print(gi,6)

       gi%v(3)=gi%v(3)-zf(3)
    endif

    gi=gi.oo.(-1)


    zf0=gi

    call kill(ex)
    call kill(gi)

  end subroutine invert_monitors


  subroutine  alex_mom_monitors
    implicit none
    integer nturn,i,j,nm,k,jm
    real(dp), pointer :: mom(:,:),r(:,:),a(:,:),at(:,:)
    real(dp)ar
    type(damap) id,m12,ex
    type(pbfield) h
    type(normalform) norm

    CALL INIT(2,2,0,0)    ! fpp
    call alloc(id,m12,ex)
    call alloc(norm)
    call alloc(h)


    nturn=SIZE (monitors(1)%r,dim=2)
    nm=size(monitors)

    do jm=1,nm

       mom=>monitors(jm)%mom
       r=>monitors(jm)%r
       a=>monitors(jm)%a
       at=>monitors(jm)%at
       mom=0.d0
       do i=1,4
          do j=i,4
             do k=1,nturn
                mom(i,j)=r(i,k)*r(j,k)+mom(i,j)
             enddo
          enddo
       enddo

       mom=mom/nturn
       do i=1,4
          do j=i,4
             mom(j,i)=mom(i,j)
          enddo
       enddo


       h%h=zero

       do i=1,2
          do j=1,2
             if(with_c==0.and.i/=j) cycle
             h%h=h%h + mom(2*i-1,2*j-1)*(1.d0.mono.(2*i))*(1.d0.mono.(2*j))
             h%h=h%h + mom(2*i,2*j)*(1.d0.mono.(2*i-1))*(1.d0.mono.(2*j-1))
             h%h=h%h - mom(2*i,2*j-1)*(1.d0.mono.(2*i-1))*(1.d0.mono.(2*j))
             h%h=h%h - mom(2*i-1,2*j)*(1.d0.mono.(2*i))*(1.d0.mono.(2*j-1))
          enddo
       enddo

       h%h=-h%h
       ar=full_abs(h%h)
       ar=ar*20.d0
       h%h=h%h/ar

       id=1
       ex=texp(h,id)
       norm=ex

       a=norm%a_t
       at=matmul(at,a)
    enddo   ! jm=1,nm


    call kill(id,m12,ex)
    call kill(norm)
    call kill(h)

  end subroutine  alex_mom_monitors

  subroutine  alex_apply_A_on_data
    implicit none
    integer jm,ier,nm
    real(dp), pointer :: mom(:,:),r(:,:),a(:,:)
    real(dp) ai(4,4)



    nm=size(monitors)

    do jm=1,nm

       r=>monitors(jm)%r
       a=>monitors(jm)%a
       call matinv(a,ai,4,4,ier)
       if(ier/=0) then
          write(6,*) " error in matinv ",jm
          stop 999
       endif

       r=matmul(ai,r)



    enddo   ! jm=1,nm



  end subroutine  alex_apply_A_on_data

  subroutine  alex_count_monitors(r,n) !locate monitors and create array monitors(n)
    implicit none
    type(layout), target :: r
    type(fibre), pointer :: p
    integer i,n,nturn
!!! find qf and qds
    !qf are horizontal monitors
    !qd are vertical monitors
    !!

    p=>r%start
    n=0
    do i=1,r%n

       if(p%mag%name(1:2)=="QF") then
          n=n+1
          !write(6,*) n,p%mag%name, "     x"
       endif
       if(p%mag%name(1:2)=="QD") then
          n=n+1
          !write(6,*) n,p%mag%name, "     y"
       endif

       p=>p%next
    enddo

    allocate( monitors(n))

    p=>r%start
    n=0
    do i=1,r%n

       if(p%mag%name(1:2)=="QF") then
          n=n+1
          call alloc_fibre_monitor_data(monitors(n),m_turn,p)  ! m_)turn is number of turns
          monitors(n)%kind=1
          ! write(6,*) n,p%mag%name
       endif
       if(p%mag%name(1:2)=="QD") then
          n=n+1
          call alloc_fibre_monitor_data(monitors(n),m_turn,p)
          monitors(n)%kind=2
          ! write(6,*) n,p%mag%name
       endif

       p=>p%next
    enddo

  end subroutine alex_count_monitors

  subroutine alloc_fibre_monitor_data(c,turn,p)
    implicit none
    TYPE(fibre_monitor_data), target :: c
    TYPE(fibre), target :: p
    integer turn,i

    nullify(c%p)
    c%p=>p
    allocate(c%turn)
    allocate(c%kind)
    allocate(c%r(4,turn))
    allocate(c%xf(4,turn))
    allocate(c%mom(4,4))
    allocate(c%A(4,4))
    allocate(c%At(4,4))
    c%turn=0
    c%kind=0
    c%r=zero
    c%xf=zero
    c%mom=zero
    c%A=zero
    do i=1,4
       c%a(i,i)=one
    enddo
    c%At=zero
    do i=1,4
       c%at(i,i)=one
    enddo
  end subroutine alloc_fibre_monitor_data

  subroutine kill_fibre_monitor_data(c)
    implicit none
    TYPE(fibre_monitor_data), pointer :: c
    integer turn

    deallocate(c%turn)
    deallocate(c%kind)
    deallocate(c%r)
    deallocate(c%xf)
    deallocate(c%mom)
    deallocate(c%a)
    deallocate(c%at)
    nullify(c%turn)
    nullify(c%kind)
    nullify(c%r)
    nullify(c%xf)
    nullify(c%mom)
    nullify(c%a)
    nullify(c%at)

  end subroutine kill_fibre_monitor_data



  subroutine  alex_mom_real_monitors_old(ring,jm,jn,x,state_in)
    implicit none
    integer ipause, mypause
    integer jm,nm,i,jinv(4),i1,i11,i2,i22,jn(4)
    type(layout), target :: ring
    real(dp), pointer :: mom(:,:),r(:,:),a(:,:),at(:,:)
    real(dp)ar,x(6),z(4),zz(lnv),zx(4),z0(4),zf(4)
    type(damap) id,m11,m2,m22,m1,ex,map1,fin
    type(pbfield) h
    type(normalform) norm
    type(fibre),pointer::p1,p11,p2,p22
    type(real_8) y(6)
    type(internal_state) state,state_in
    type(gmap) g
!!! will produce the real data at jm

    z0=zero

    nm=size(monitors)

    state=state_in+only_4d0

    CALL INIT(STATE,1,0)

    call alloc(y)
    call alloc(id,m11,m2,m22,m1,ex)
    call alloc(map1)
    call alloc(g)

    CALL FIND_ORBIT(ring,x,1,STATE,c_1d_5)
    write(6,*) monitors(jm)%kind
    if(monitors(jm)%kind==1) then
       jn(1)=jm

       p1=>monitors(jm)%p
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1) then
             p11=>monitors(mod_n(i,nm))%p
             jn(2)=mod_n(i,nm)
             i11=i
             exit
          endif
       enddo
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2) then
             p2=>monitors(mod_n(i,nm))%p
             jn(3)=mod_n(i,nm)
             i2=i
             exit
          endif
       enddo
       do i=i2+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2) then
             p22=>monitors(mod_n(i,nm))%p
             jn(4)=mod_n(i,nm)
             i22=i
             exit
          endif
       enddo
       write(6,*)
       write(6,*) p1%mag%name,p11%mag%name
       write(6,*) p2%mag%name,p22%mag%name
       write(6,*) p1%pos,p11%pos
       write(6,*) p2%pos,p22%pos

       ! track from ip to p1 : jm x monitor

       call track(x,state,ring%start,p1)
       !
       ! x is closed orbit at p1



       id=1

       y=x+id
       call track(y,state,p1,p11)
       m11=y
       ! m11 is map from p1 to p11  (second x monitor)
       m11=z0  ! removed zeroth order

       y=x+id
       call track(y,state,p1,p2)
       m2=y
       m2=z0
       ! m2 is map from p1 to p2  (first y monitor)


       y=x+id
       call track(y,state,p1,p22)
       m22=y
       m22=z0
       ! m22 is map from p1 to p22  (second y monitor)

!!!!!!!!!!!!   testing and computing m11 !!!!!!!!!!!!!!!!!
       g=1


       jinv=0
       jinv(1)=1
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=m11*ex
       m11=m11**jinv
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=ex*m11*ex

       ! m11 is now the map from x,x11,y,py

       g%v(2)=m11%v(2)     !   px_0(x_0,x_p11,y_0,py_0)

       ! testing tracking from jm to p2

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=m22*ex
       m22=m22**jinv
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=ex*m22*ex
       ! m22 is now the map from x,px,y,y22


       g%v(4)=m22%v(4)

!!!!!!!!!!!!   testing and computing m2  !!!!!!!!!!!!!!!!!

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m2=m2*ex
       m2=m2**jinv
       !m11=m11**(-1)
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m2=ex*m2*ex

       ! m2 is now the map from x,px,y2,py

       g%v(3)=m2%v(3)


    else  !!!! if(monitors(jm)%kind==2)
       jn(3)=jm

       write(6,*) " crotte "
       ipause=mypause(0)
       p2=>monitors(jm)%p
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2) then
             p22=>monitors(mod_n(i,nm))%p
             jn(4)=mod_n(i,nm)
             i22=i
             exit
          endif
       enddo
       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1) then
             p1=>monitors(mod_n(i,nm))%p
             jn(2)=mod_n(i,nm)
             i1=i
             exit
          endif
       enddo
       do i=i22+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1) then
             p11=>monitors(mod_n(i,nm))%p
             jn(1)=mod_n(i,nm)
             i11=i
             exit
          endif
       enddo
       write(6,*)
       write(6,*) p1%mag%name,p11%mag%name
       write(6,*) p2%mag%name,p22%mag%name
       write(6,*) p1%pos,p11%pos
       write(6,*) p2%pos,p22%pos

       ! track from ip to p1 : jm x monitor

       call track(x,state,ring%start,p2)
       !
       ! x is closed orbit at p2



       id=1

       y=x+id
       call track(y,state,p2,p22)
       m22=y
       ! m22 is map from p2 to p22  (second y monitor)
       m22=z0  ! removed zeroth order

       y=x+id
       call track(y,state,p2,p1)
       m1=y
       m1=z0
       ! m1 is map from p2 to p1  (first x monitor)


       y=x+id
       call track(y,state,p2,p11)
       m11=y
       m11=z0
       ! m11 is map from p2 to p11  (second x monitor)

!!!!!!!!!!!!   testing and computing m22 !!!!!!!!!!!!!!!!!  asdad
       g=1


       jinv=0
       jinv(3)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=m22*ex
       m22=m22**jinv
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.4
       ex%v(4)=1.0d0.mono.3
       m22=ex*m22*ex

       g%v(4)=m22%v(4)     !   py_0(x_0,p_x0,y_0,y_p22)

       ! m22 is now the map from x,px,y,y22

!!!!!!!!!!!!   testing and computing m1  !!!!!!!!!!!!!!!!!

       jinv=0
       jinv(1)=1
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m1=m1*ex
       m1=m1**jinv
       ex%v(1)=1.0d0.mono.2
       ex%v(2)=1.0d0.mono.1
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m1=ex*m1*ex
       ! m1 is now the map from x,x1,y,py

       g%v(2)=m1%v(2)


!!!!!!!!!!!!   testing and computing m11  !!!!!!!!!!!!!!!!!
       jinv=0
       jinv(1)=1
       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=m11*ex
       m11=m11**jinv

       ex%v(1)=1.0d0.mono.1
       ex%v(2)=1.0d0.mono.2
       ex%v(3)=1.0d0.mono.3
       ex%v(4)=1.0d0.mono.4
       m11=ex*m11*ex

       ! m11 is now the map from x2,px,y,py

       g%v(1)=m11%v(1)


    endif !!!! if(monitors(jm)%kind==2)


    zf=0.d0
    zf(1)=monitors(jn(1))%r(1,1)
    zf(2)=monitors(jn(2))%r(1,1)
    zf(3)=monitors(jn(3))%r(3,1)
    zf(4)=monitors(jn(4))%r(3,1)
    write(6,*) zf(1:4)
    write(6,*) monitors(jm)%r(1:4,1)
    ipause=mypause(200)
    call invert_monitors(monitors(jm)%kind,zf,g)
    !call invert_monitors_old(monitors(jm)%kind,g)
    !do i=1,4
    !zz(i)=g%v(i)*zf
    !enddo
    !write(6,*) zz(1:4)
    write(6,*) zf(1:4)


    ipause=mypause(201)





    call kill(y)
    call kill(id,m11,m2,m22,m1,ex)
    call kill(map1)
    call kill(g)



  end subroutine  alex_mom_real_monitors_old

  subroutine invert_monitors_old(kind,g)
    implicit none
    type(damap) ex
    type(gmap) gi,g
    integer kind


    call alloc(ex)
    call alloc(gi)
    ! zf contains measuments xf(x_1,x_11,y_2,y_22)
    gi=g

    if(kind==1) then
       !write(6,*) zf(1:4)
       ex=0
       ex%v(3)=ex%v(3)+(1.d0.mono.3)
       ex%v(4)=ex%v(4)+(1.d0.mono.4)
       gi%v(2)=(gi%v(2).o.ex)-(1.d0.mono.2)
       ex=0
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       ex%v(4)=ex%v(4)+(1.d0.mono.4)
       gi%v(3)=(gi%v(3).o.ex)-(1.d0.mono.3)
       ex=0
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       ex%v(3)=ex%v(3)+(1.d0.mono.3)
       gi%v(4)=(gi%v(4).o.ex)-(1.d0.mono.4)

       gi%v(1)=gi%v(1)

    else

       ex=0
       ex%v(1)=ex%v(1)+(1.d0.mono.1)
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       gi%v(4)=(gi%v(4).o.ex)-(1.d0.mono.4)
       ex=0
       ex%v(2)=ex%v(2)+(1.d0.mono.2)
       ex%v(4)=ex%v(4)+(1.d0.mono.4)
       gi%v(1)=(gi%v(1).o.ex)-(1.d0.mono.1)
       ex=0
       ex%v(4)=ex%v(4)+(1.d0.mono.4)
       ex%v(1)=ex%v(1)+(1.d0.mono.1)
       gi%v(2)=(gi%v(2).o.ex)-(1.d0.mono.2)
       gi%v(3)=gi%v(3)
    endif

    gi=gi.oo.(-1)

    g=gi

    call kill(ex)
    call kill(gi)

  end subroutine invert_monitors_old


end module S_fitting_new
