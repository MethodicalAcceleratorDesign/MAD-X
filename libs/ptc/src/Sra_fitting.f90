!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_fitting_new
  USE ptc_spin
  IMPLICIT NONE
  public
  integer:: m_turn,m_skip=0
  integer :: with_c=1
  TYPE fibre_monitor_data
     type(fibre), pointer :: p    ! fibre location
     integer, pointer ::  turn,kind  ! kind=1 x, kind = 2 y
     real(dp), pointer :: bpm(:,:)  ! store fake experiment from alex_track_monitors
     real(dp), pointer :: r(:,:)  ! store fake experiment from alex_track_monitors
     real(dp), pointer :: xf(:,:)  ! real data put here
     real(dp), pointer :: xn(:,:)  ! real data put here
     real(dp), pointer :: mom(:,:)
     real(dp), pointer :: A(:,:)
     real(dp), pointer :: At(:,:)
     logical full
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

    X=0.0_dp
    CALL INIT(STATE,2,NU,BERZ,N2,NP2)
    CALL ALLOC(ID); call alloc(nf);call alloc(Y);call alloc(eq,5);
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    !   CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT1,node2=IT2)
    U=.NOT.CHECK_STABLE
    nf=y
    TA(1)=0.75_dp
    TA(2)=0.68_dp
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
    x=0.0_dp
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT3,node2=IT4)
    U=.NOT.CHECK_STABLE

    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
    nf=y
    eq(1)=(eq(1)*8 + (1.0_dp+nf%dhdj%v(1)))*3-targ(1)
    eq(2)=(eq(2)*8 + (1.0_dp+nf%dhdj%v(2)))*3-targ(2)
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
       g%v(i)=g%v(i)-(1.0_dp-sc)*(g%v(i).sub.'0')
    enddo

    call alloc(t)
    do i=1,np
       g%v(i)=1.0_dp.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(1.0_dp.mono.j)*t
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

    X=0.0_dp
    CALL INIT(STATE,2,NU,BERZ,N2,NP2)
    CALL ALLOC(ID); call alloc(nf);call alloc(Y);call alloc(eq,5);
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    !   CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT1,node2=IT2)
    U=.NOT.CHECK_STABLE

    nf=y
    TA(1)=0.75_dp
    TA(2)=0.68_dp
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
    x=0.0_dp
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT3,node2=IT4)
    U=.NOT.CHECK_STABLE

    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT3,POS2=IT4)
    nf=y
    eq(1)=(eq(1)*8 + (1.0_dp+nf%dhdj%v(1)))*3-targ(1)
    eq(2)=(eq(2)*8 + (1.0_dp+nf%dhdj%v(2)))*3-targ(2)
    eq(3)=eq(3)-ta(1)

    eq(4)=eq(4)-(nf%A_T%V(1).par.'1000')**2-(nf%A_T%V(1).par.'0100')**2
    eq(5)=eq(5)-(nf%A_T%V(3).par.'0010')**2-(nf%A_T%V(3).par.'0001')**2
    !call print(eq(4),6)
    !call print(eq(5),6)
    nxx=0.0_dp
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
       g%v(i)=g%v(i)-(1.0_dp-sc)*(g%v(i).sub.'0')
    enddo

    call alloc(t)
    do i=1,np
       g%v(i)=1.0_dp.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(1.0_dp.mono.j)*t
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

    X=0.0_dp
    CALL INIT(STATE,2,NU,BERZ,N2,NP2)
    CALL ALLOC(ID); call alloc(nf);call alloc(Y);
    ID=1
    Y=X+ID
    !( R,X,U,K,POS1,POS2,T1,T2,P1,P2,IN_P1,IN_P2,POS1_FIBRE,POS2_FIBRE)
    !    CALL TRACK_BEAM_x(R,Y,U,+STATE,POS1=IT1,POS2=IT2)
    CALL TRACK_probe_x(R,Y,+STATE,node1=IT1,node2=IT2)
    U=.NOT.CHECK_STABLE

    nf=y
    TA(1)=0.75_dp
    TA(2)=0.68_dp
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
       g%v(i)=g%v(i)-(1.0_dp-sc)*(g%v(i).sub.'0')
    enddo


    call alloc(t)
    do i=1,np
       g%v(i)=1.0_dp.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(1.0_dp.mono.j)*t
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
    X=0.0_dp
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

    eq(1)=(TA(1)*8 + (1.0_dp+nf%dhdj%v(1)))*3
    eq(2)=(TA(2)*8 + (1.0_dp+nf%dhdj%v(2)))*3

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
       g%v(i)=g%v(i)-(1.0_dp-sc)*(g%v(i).sub.'0')
    enddo

    call alloc(t)
    do i=1,np
       g%v(i)=1.0_dp.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(1.0_dp.mono.j)*t
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
    CLOSED(:)=0.0_dp
    it=0
100 continue
    it=it+1

    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)
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
       g%v(i)=1.0_dp.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(1.0_dp.mono.j)*t
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
             an=p%mag%an(k)/p1%mag%l*2.0_dp
             bn=p%mag%bn(k)/p1%mag%l*2.0_dp
             call add(p1,-k,1,an)
             call add(p1,k,1,bn)
             !   call add(p,-k,0,zero)
             !   call add(p,k,0,zero)
             !   call add(p2,-k,0,zero)
             !   call add(p2,k,0,zero)
          enddo
          ! el%k3%thin_h_foc,el%k3%thin_v_foc,el%k3%thin_h_angle,el%k3%thin_v_angle
          if(p%mag%kind==kind3) then
             bn=p%mag%k3%thin_h_angle/p1%mag%l*2.0_dp
             call add(p1,1,1,bn)
             an=p%mag%k3%thin_v_angle/p1%mag%l*2.0_dp
             call add(p1,-1,1,an)
             doit=(p%mag%k3%thin_h_foc/=0.0_dp.or.p%mag%k3%thin_v_foc/=0.0_dp)  !.or.p%mag%k3%thin_v_angle/=zero)
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
             p%mag%k3%thin_v_angle=0.0_dp
             p%mag%k3%thin_h_angle=0.0_dp
             p%magp%k3%thin_v_angle=0.0_dp
             p%magp%k3%thin_h_angle=0.0_dp
             doit=(p%mag%k3%thin_h_foc/=0.0_dp.or.p%mag%k3%thin_v_foc/=0.0_dp)
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
       if(p%mag%l==0.0_dp.and.p%mag%kind>kind1) then
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
    CALL FIND_ORBIT(R,y,1,STATE,1e-5_dp)
    write(6,*)" closed orbit "
    write(6,*) y(1:3)
    write(6,*) y(4:6)
    x(1:4)=x(1:4)+y(1:4)
    do i=1,nturn
       call track(x,state,p1,monitors(1)%p)

       call track(y,state,p1,monitors(1)%p)
       monitors(1)%r(1:4,i)=x(1:4)-y(1:4)
       monitors(1)%bpm(1,i)=monitors(1)%r(1,i)
       monitors(1)%bpm(2,i)=monitors(1)%r(3,i)

       do j=1,nm-1
          call track(x,state,monitors(j)%p,monitors(j+1)%p)
          call track(y,state,monitors(j)%p,monitors(j+1)%p)
          monitors(j+1)%r(1:4,i)=x(1:4)-y(1:4)
          monitors(j+1)%bpm(1,i)=monitors(j)%r(1,i)
          monitors(j+1)%bpm(2,i)=monitors(j)%r(3,i)

       enddo
       call track(x,state,monitors(nm)%p,p1)
       call track(y,state,monitors(nm)%p,p1)

    enddo





  end subroutine  alex_track_monitors



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
       zx=0.0_dp
       zx(1)=zf(1)
       zx(2)=zf(2)   !!!
       ex=0
       ex=zx
       ex%v(3)=ex%v(3)+(1.0_dp.mono.3)
       ex%v(4)=ex%v(4)+(1.0_dp.mono.4)
       gi%v(2)=(gi%v(2).o.ex)-(1.0_dp.mono.2)
       zx=0.0_dp
       zx(1)=zf(1)
       zx(3)=zf(3)   !!!
       ex=0
       ex=zx
       ex%v(2)=ex%v(2)+(1.0_dp.mono.2)
       ex%v(4)=ex%v(4)+(1.0_dp.mono.4)

       gi%v(3)=(gi%v(3).o.ex)-(1.0_dp.mono.3)
       zx=0.0_dp
       zx(1)=zf(1)
       zx(4)=zf(4)   !!!
       ex=0
       ex=zx
       ex%v(2)=ex%v(2)+(1.0_dp.mono.2)
       ex%v(3)=ex%v(3)+(1.0_dp.mono.3)
       gi%v(4)=(gi%v(4).o.ex)-(1.0_dp.mono.4)
       !call print(gi,6)
       gi%v(1)=gi%v(1)-zf(1)

    else
       !write(6,*) zf(1:4)

       zx=0.0_dp
       zx(3)=zf(3)
       zx(4)=zf(4)   !!!
       ex=0
       ex=zx
       ex%v(1)=ex%v(1)+(1.0_dp.mono.1)
       ex%v(2)=ex%v(2)+(1.0_dp.mono.2)
       gi%v(4)=(gi%v(4).o.ex)-(1.0_dp.mono.4)
       zx=0.0_dp
       zx(1)=zf(1)
       zx(3)=zf(3)   !!!
       ex=0
       ex=zx
       ex%v(2)=ex%v(2)+(1.0_dp.mono.2)
       ex%v(4)=ex%v(4)+(1.0_dp.mono.4)

       gi%v(1)=(gi%v(1).o.ex)-(1.0_dp.mono.1)
       zx=0.0_dp
       zx(3)=zf(3)
       zx(2)=zf(2)   !!!
       ex=0
       ex=zx
       ex%v(4)=ex%v(4)+(1.0_dp.mono.4)
       ex%v(1)=ex%v(1)+(1.0_dp.mono.1)
       gi%v(2)=(gi%v(2).o.ex)-(1.0_dp.mono.2)
       !call print(gi,6)

       gi%v(3)=gi%v(3)-zf(3)
    endif

    gi=gi.oo.(-1)


    zf0=gi

    call kill(ex)
    call kill(gi)

  end subroutine invert_monitors


  subroutine  alex_mom_monitors(jm,nt)
    implicit none
    integer nturn,i,j,nm,k,jm,nt
    real(dp), pointer :: mom(:,:),r(:,:),a(:,:),at(:,:)
    real(dp)ar, momt(6,6),kicke(3),a6(6,6),ai6(6,6),br(6,6)
    type(damap) id,m12,ex
    type(pbfield) h
    type(normalform) norm

    CALL INIT(2,2,0,0)    ! fpp
    call alloc(id,m12,ex)
    call alloc(norm)
    call alloc(h)


    nturn=SIZE (monitors(1)%r,dim=2)
    nm=size(monitors)

    !do jm=1,nm

    mom=>monitors(jm)%mom
    r=>monitors(jm)%xf
    a=>monitors(jm)%a
    at=>monitors(jm)%at
    mom=0.0_dp
    do i=1,4
       do j=i,4
          do k=1,nt
             mom(i,j)=r(i,k)*r(j,k)+mom(i,j)
          enddo
       enddo
    enddo

    mom=mom/nt
    do i=1,4
       do j=i,4
          mom(j,i)=mom(i,j)
       enddo
    enddo

    i=0

    if(i==0) then

       momt=0.0_dp
       do i=1,2
          do j=1,2
             if(with_c==0.and.i/=j) cycle
             ! momt(2*i-1,2*j-1)=mom(2*i-1,2*j-1)
             ! momt(2*i-1,2*j)=mom(2*i-1,2*j)
             ! momt(2*i,2*j-1)=mom(2*i,2*j-1)
             ! momt(2*i,2*j)=mom(2*i,2*j)
             momt(2*i,2*j)=mom(2*i-1,2*j-1)
             momt(2*i,2*j-1)=-mom(2*i-1,2*j)
             momt(2*i-1,2*j)=-mom(2*i,2*j-1)
             momt(2*i-1,2*j-1)=mom(2*i,2*j)
          enddo
       enddo

       momt(5,5)=mom(1,1)
       momt(6,6)=mom(1,1)

       call diagonalise_envelope_a(momt,br,a6,ai6,kicke)
       a=a6(1:4,1:4)
       at=matmul(at,a)
    else
       h%h=0.0_dp

       do i=1,2
          do j=1,2
             if(with_c==0.and.i/=j) cycle
             h%h=h%h + mom(2*i-1,2*j-1)*(1.0_dp.mono.(2*i))*(1.0_dp.mono.(2*j))
             h%h=h%h + mom(2*i,2*j)*(1.0_dp.mono.(2*i-1))*(1.0_dp.mono.(2*j-1))
             h%h=h%h - mom(2*i,2*j-1)*(1.0_dp.mono.(2*i-1))*(1.0_dp.mono.(2*j))
             h%h=h%h - mom(2*i-1,2*j)*(1.0_dp.mono.(2*i))*(1.0_dp.mono.(2*j-1))
          enddo
       enddo

       h%h=-h%h
       ar=full_abs(h%h)
       ar=ar*20.0_dp
       h%h=h%h/ar

       id=1
       ex=texp(h,id)
       norm=ex

       a=norm%a_t
       at=matmul(at,a)
       !enddo   ! jm=1,nm

    endif

    !write(6,*) "betax", a(1,1)**2+a(1,2)**2, sqrt(a(1,1)**2+a(1,2)**2)
    !write(6,*) "alphax", -a(1,1)*a(2,1)-a(1,2)*a(2,2)
    !write(6,*) "betay", a(3,3)**2+a(3,4)**2, sqrt(a(3,3)**2+a(3,4)**2)
    !write(6,*) "alphay", -a(3,3)*a(4,3)-a(3,4)*a(4,4)
    !write(6,*) "betaxy", a(1,3)**2+a(1,4)**2
    !write(6,*) "betayx", a(3,1)**2+a(3,2)**2


    call kill(id,m12,ex)
    call kill(norm)
    call kill(h)

  end subroutine  alex_mom_monitors

  subroutine  alex_apply_A_on_data(jm)
    implicit none
    integer jm,ier,nm
    real(dp), pointer :: mom(:,:),r(:,:),a(:,:),xn(:,:)
    real(dp) ai(4,4)



    nm=size(monitors)

    !do jm=1,nm

    r=>monitors(jm)%xf
    xn=>monitors(jm)%xn
    a=>monitors(jm)%a
    call matinv(a,ai,4,4,ier)
    if(ier/=0) then
       write(6,*) " error in matinv ",jm
       stop 999
    endif

    xn=matmul(ai,r)



    !enddo   ! jm=1,nm



  end subroutine  alex_apply_A_on_data

  subroutine  alex_count_monitors(r,n,kind) !locate monitors and create array monitors(n)
    implicit none
    type(layout), target :: r
    type(fibre), pointer :: p
    integer i,n,nturn
    integer, optional :: kind
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
          if(present(kind)) then
             monitors(n)%kind=kind
          else
             monitors(n)%kind=1
          endif
          ! write(6,*) n,p%mag%name
       endif
       if(p%mag%name(1:2)=="QD") then
          n=n+1
          call alloc_fibre_monitor_data(monitors(n),m_turn,p)
          if(present(kind)) then
             monitors(n)%kind=kind
          else
             monitors(n)%kind=2
          endif
          ! write(6,*) n,p%mag%name
       endif

       p=>p%next
    enddo

  end subroutine alex_count_monitors

  subroutine  alex_count_jparc_monitors(r,n,kind) !locate monitors and create array monitors(n)
    implicit none
    type(layout), target :: r
    type(fibre), pointer :: p
    integer i,n,nturn,mf,i1,i2,j,i3,kind
    integer, allocatable :: ind(:)
    character(255) line
!!! find qf and qds
    !qf are horizontal monitors
    !qd are vertical monitors
    !!

    call kanalnummer(mf,"J_parc_monitor_id.txt")

    n=0
    do i=1,r%n
       read(mf,*,end=100) i1,i2

       n=n+1
    enddo
100 write(6,*) n
    allocate(ind(n))
    ind=0

    rewind mf

    do i=1,n
       read(mf,*) i1,i2
       ind(i1)=i2
    enddo


    close(mf)


    allocate( monitors(n))

    p=>r%start
    do i=1,r%n

       if(p%mag%name(1:2)=="QF") then


          i1=0
          line=p%mag%name(4:4)
          read(line,*) i2
          i1=i2*100+i1
          line=p%mag%name(5:5)
          read(line,*) i2
          i1=i2*10+i1
          line=p%mag%name(6:6)
          read(line,*) i2
          i1=i2+i1
          i3=i1
          if(i1>n) then
             i3=n
          endif
          i2=0
          do j=i3,1,-1
             if(ind(j)==i1) then
                i2=ind(j)
                exit
             endif
          enddo
          if(i2/=0) then
             call alloc_fibre_monitor_data(monitors(j),m_turn,p)
             monitors(j)%kind=kind
          endif
       endif

       if(p%mag%name(1:2)=="QD") then

          i1=0
          line=p%mag%name(4:4)
          read(line,*) i2
          i1=i2*100+i1
          line=p%mag%name(5:5)
          read(line,*) i2
          i1=i2*10+i1
          line=p%mag%name(6:6)
          read(line,*) i2
          i1=i2+i1

          i3=i1
          if(i1>n) then
             i3=n
          endif
          i2=0
          do j=i3,1,-1
             if(ind(j)==i1) then
                i2=ind(j)
                exit
             endif
          enddo

          if(i2/=0) then
             call alloc_fibre_monitor_data(monitors(j),m_turn,p)
             monitors(j)%kind=kind
          endif


       endif


       p=>p%next
    enddo



    deallocate(ind)

  end subroutine alex_count_jparc_monitors

  subroutine  alex_read_r_jparc(filename1,filename2,n_max,nav)
    implicit none
    integer jm,mf1,mf2,n,i,na,i1,i2
    CHARACTER(*) filename1,filename2
    integer, optional :: n_max,nav
    real(dp) x,y,jj(4),xa(2),ya(2),norm(2)
    real(dp), allocatable :: bpmx(:), bpmy(:)

    allocate(bpmx(size(monitors)),bpmy(size(monitors)))

    call kanalnummer(mf1,filename1)
    call kanalnummer(mf2,filename2)

    n=m_turn
    na=0
    if(present(n_max)) n=n_max
    if(present(nav)) na=nav


    x=0.0_dp
    y=0.0_dp
    do i=1,m_turn+m_skip
       read(mf1,*)i1,i2,bpmx
       read(mf2,*)i1,i2,bpmy
       if(i>m_skip) then
          do jm=1,size(monitors)
             monitors(jm)%r(1,i-m_skip)=bpmx(jm)/1.d3
             monitors(jm)%r(3,i-m_skip)=bpmy(jm)/1.d3
             monitors(jm)%bpm(1,i-m_skip)=bpmx(jm)/1.d3
             monitors(jm)%bpm(2,i-m_skip)=bpmy(jm)/1.d3
          enddo
       endif
    enddo

write(6,*) size(monitors), "monitors "
    do jm=1,size(monitors)
       norm=0.0_dp
       do i=1,m_turn
          norm(1)=norm(1)+abs(monitors(jm)%r(1,i))
          norm(2)=norm(2)+abs(monitors(jm)%r(3,i))
       enddo
       if(norm(1)<1.e-10.or. norm(2)<1.e-10) then
          write(6,*) " monitor ",jm, " has no data "
          monitors(jm)%full=.false.
       endif
    enddo

    close(mf1)
    close(mf2)
    deallocate(bpmx,bpmy)
    ! 171683 59 68 1 2.6794 -1.1248

  end   subroutine  alex_read_r_jparc


  subroutine scale_bpm
    implicit none
    integer i,n,ib
    real x,y,x0,y0,sx,sy,ax,ay

    if(.not.allocated(monitors)) return

    do ib=1,size(monitors)
       n=m_turn/10


       do i=1,n
          x0=abs(monitors(ib)%r(1,i))+x0
          y0=abs(monitors(ib)%r(3,i))+y0
       enddo

       do i=m_turn-n+1,m_turn
          x=abs(monitors(ib)%r(1,i))+x
          y=abs(monitors(ib)%r(3,i))+y
       enddo

       sx=(x0/x-1)/(m_turn-1)
       ax=1-sx
       sy=(y0/y-1)/(m_turn-1)
       ay=1-sy

       do i=1,m_turn
          monitors(ib)%r(1,i)=monitors(ib)%r(1,i)*(i*sx+ax)
          monitors(ib)%r(3,i)=monitors(ib)%r(3,i)*(i*sy+ay)
       enddo

    enddo

  end subroutine scale_bpm



  subroutine  alex_average_r_jparc
    implicit none
    integer jm,mf1,mf2,n,i,na,i1,i2
    real(dp) x,y,jj(4),xa(2),ya(2)


    n=m_turn

    do jm=1,size(monitors)
       monitors(jm)%r(1,:)=monitors(jm)%bpm(1,:)
       monitors(jm)%r(3,:)=monitors(jm)%bpm(2,:)
    enddo
    do jm=1,size(monitors)

       x=0.0_dp
       y=0.0_dp
       do i=1,n
          x=monitors(jm)%r(1,i)+x
          y=monitors(jm)%r(3,i)+y
       enddo

       x=x/n
       y=y/n
       write(6,*) "Averages ",x,y
       do i=1,n
          monitors(jm)%r(1,i)=(monitors(jm)%r(1,i)-x)
          monitors(jm)%r(3,i)=(monitors(jm)%r(3,i)-y)
       enddo

    enddo


    ! 171683 59 68 1 2.6794 -1.1248

  end   subroutine  alex_average_r_jparc



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
    allocate(c%bpm(2,turn))
    allocate(c%xf(4,turn))
    allocate(c%xn(4,turn))
    allocate(c%mom(4,4))
    allocate(c%A(4,4))
    allocate(c%At(4,4))
    c%turn=0
    c%kind=0
    c%r=0.0_dp
    c%bpm=0.0_dp
    c%xf=0.0_dp
    c%xn=0.0_dp
    c%mom=0.0_dp
    c%A=0.0_dp
    c%full=.true.
    do i=1,4
       c%a(i,i)=1.0_dp
    enddo
    c%At=0.0_dp
    do i=1,4
       c%at(i,i)=1.0_dp
    enddo
  end subroutine alloc_fibre_monitor_data

  subroutine kill_fibre_monitor_data(c)
    implicit none
    TYPE(fibre_monitor_data), pointer :: c
    integer turn

    deallocate(c%turn)
    deallocate(c%kind)
    deallocate(c%r)
    deallocate(c%bpm)
    deallocate(c%xf)
    deallocate(c%xn)
    deallocate(c%mom)
    deallocate(c%a)
    deallocate(c%at)
    nullify(c%turn)
    nullify(c%kind)
    nullify(c%r)
    nullify(c%xf)
    nullify(c%xn)
    nullify(c%mom)
    nullify(c%a)
    nullify(c%at)

  end subroutine kill_fibre_monitor_data


!!!1 this is the correct one
  subroutine  alex_mom_real_monitors(ring,jm,jn,x,state_in,no,nt)
    implicit none
    integer jm,nm,i,jinv(4),i1,i11,i2,i22,jn(4),mf1,mf2,no,nt
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

    z0=0.0_dp

    nm=size(monitors)

    state=state_in+only_4d0

    CALL INIT(STATE,no,0)

    call alloc(y)
    call alloc(id,m11,m2,m22,m1,ex)
    call alloc(map1)
    call alloc(g)

    CALL FIND_ORBIT(ring,x,1,STATE,1e-5_dp)
    write(6,*) monitors(jm)%kind
    if(monitors(jm)%kind==1.or.monitors(jm)%kind==3) then
       jn(1)=jm

       p1=>monitors(jm)%p
       if(monitors(jm)%kind==3) then
          p2=>monitors(jm)%p
          i2=jm
          jn(3)=jm
       else
          do i=jm+1,jm+nm
             if(monitors(mod_n(i,nm))%kind==2) then
                p2=>monitors(mod_n(i,nm))%p
                jn(3)=mod_n(i,nm)
                i2=i
                exit
             endif
          enddo

       endif

       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1.or.monitors(mod_n(i,nm))%kind==3) then
             p11=>monitors(mod_n(i,nm))%p
             jn(2)=mod_n(i,nm)
             i11=i
             exit
          endif
       enddo

       if(monitors(mod_n(jn(2),nm))%kind==3) then
          p22=>monitors(mod_n(jn(2),nm))%p
          jn(4)=jn(2)
          i22=i2+1
       else
          do i=i2+1,jm+nm
             if(monitors(mod_n(i,nm))%kind==2) then
                p22=>monitors(mod_n(i,nm))%p
                jn(4)=mod_n(i,nm)
                i22=i
                exit
             endif
          enddo
       endif

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
       ex%v(1)=1.0_dp.mono.2
       ex%v(2)=1.0_dp.mono.1
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m11=m11*ex
       m11=m11**jinv
       ex%v(1)=1.0_dp.mono.2
       ex%v(2)=1.0_dp.mono.1
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m11=ex*m11*ex

       ! m11 is now the map from x,x11,y,py

       g%v(2)=m11%v(2)     !   px_0(x_0,x_p11,y_0,py_0)

       ! testing tracking from jm to p2

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.4
       ex%v(4)=1.0_dp.mono.3
       m22=m22*ex
       m22=m22**jinv
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.4
       ex%v(4)=1.0_dp.mono.3
       m22=ex*m22*ex
       ! m22 is now the map from x,px,y,y22


       g%v(4)=m22%v(4)

!!!!!!!!!!!!   testing and computing m2  !!!!!!!!!!!!!!!!!

       jinv=0
       jinv(3)=1
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m2=m2*ex
       m2=m2**jinv
       !m11=m11**(-1)
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m2=ex*m2*ex

       ! m2 is now the map from x,px,y2,py

       g%v(3)=m2%v(3)


    else  !!!! if(monitors(jm)%kind==2)
       jn(3)=jm

       !write(6,*) " crotte "
       !pause
       p2=>monitors(jm)%p

       do i=jm+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==2.or.monitors(mod_n(i,nm))%kind==3) then
             p22=>monitors(mod_n(i,nm))%p
             jn(4)=mod_n(i,nm)
             i22=i
             exit
          endif
       enddo

       if(monitors(mod_n(jn(4),nm))%kind==3) then

          p1=>monitors(mod_n(jn(4),nm))%p
          i1=i22
          jn(2)=jn(4)
       else
          do i=jm+1,jm+nm
             if(monitors(mod_n(i,nm))%kind==1) then
                p1=>monitors(mod_n(i,nm))%p
                jn(2)=mod_n(i,nm)
                i1=i
                exit
             endif
          enddo

       endif


       do i=i22+1,jm+nm
          if(monitors(mod_n(i,nm))%kind==1.or.monitors(mod_n(i,nm))%kind==3) then
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
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.4
       ex%v(4)=1.0_dp.mono.3
       m22=m22*ex
       m22=m22**jinv
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.4
       ex%v(4)=1.0_dp.mono.3
       m22=ex*m22*ex

       g%v(4)=m22%v(4)     !   py_0(x_0,p_x0,y_0,y_p22)

       ! m22 is now the map from x,px,y,y22

!!!!!!!!!!!!   testing and computing m1  !!!!!!!!!!!!!!!!!

       jinv=0
       jinv(1)=1
       ex%v(1)=1.0_dp.mono.2
       ex%v(2)=1.0_dp.mono.1
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m1=m1*ex
       m1=m1**jinv
       ex%v(1)=1.0_dp.mono.2
       ex%v(2)=1.0_dp.mono.1
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m1=ex*m1*ex
       ! m1 is now the map from x,x1,y,py

       g%v(2)=m1%v(2)


!!!!!!!!!!!!   testing and computing m11  !!!!!!!!!!!!!!!!!
       jinv=0
       jinv(1)=1
       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m11=m11*ex
       m11=m11**jinv

       ex%v(1)=1.0_dp.mono.1
       ex%v(2)=1.0_dp.mono.2
       ex%v(3)=1.0_dp.mono.3
       ex%v(4)=1.0_dp.mono.4
       m11=ex*m11*ex

       ! m11 is now the map from x2,px,y,py

       g%v(1)=m11%v(1)


    endif !!!! if(monitors(jm)%kind==2)

    !call kanalnummer(mf1,"original.dat")
    !call kanalnummer(mf2,"new.dat")
    zf=0.0_dp
    zf(1)=monitors(jn(1))%r(1,1)
    zf(2)=monitors(jn(2))%r(1,1)
    zf(3)=monitors(jn(3))%r(3,1)
    zf(4)=monitors(jn(4))%r(3,1)
    write(6,*) zf(1:4)
    write(6,*) monitors(jm)%r(1:4,1)
    !pause 200
    call invert_monitors(monitors(jm)%kind,zf,g)
    write(6,*) zf(1:4)
    !pause 201

    write(6,*) " M_turn used in Eikonal",m_turn
    do i=1,nt
       zf=0.0_dp
       zf(1)=monitors(jn(1))%r(1,i)
       zf(2)=monitors(jn(2))%r(1,i)
       zf(3)=monitors(jn(3))%r(3,i)
       zf(4)=monitors(jn(4))%r(3,i)
       call invert_monitors(monitors(jm)%kind,zf,g)
       ! WRITE(mf1,'(4(1x,E15.8))')monitors(jm)%r(1:4,i)
       ! WRITE(mf2,'(4(1x,E15.8))')zf(1:4)
       monitors(jm)%xf(1:4,i)=zf(1:4)
    enddo

    !close(mf1)
    !close(mf2)

    call kill(y)
    call kill(id,m11,m2,m22,m1,ex)
    call kill(map1)
    call kill(g)



  end subroutine  alex_mom_real_monitors

  subroutine  alex_print_xf(jm,filename,n_max)
    implicit none
    integer jm,mf1,n,i
    CHARACTER(*) FILENAME
    integer, optional :: n_max
    call kanalnummer(mf1,filename)

    n=m_turn
    if(present(n_max)) n=n_max

    do i=1,n
       WRITE(mf1,'(4(1x,E15.8))')monitors(jm)%xf(1:4,i)
    enddo

    close(mf1)

  end   subroutine  alex_print_xf

  subroutine  alex_read_r(sc,jm,filename,n_max,nav)
    implicit none
    integer jm,mf1,n,i,na
    CHARACTER(*), optional ::  FILENAME
    integer, optional :: n_max,nav
    real(dp) x,y,jj(4),xa(2),ya(2),sc
    call kanalnummer(mf1,filename)

    n=m_turn
    na=0
    if(present(n_max)) n=n_max
    if(present(nav)) na=nav

    if(present(filename)) then
       x=0.0_dp
       y=0.0_dp
       do i=1,n
          read(mf1,*)jj,monitors(jm)%r(1,i),monitors(jm)%r(3,i)
          monitors(jm)%r(3,i)=sc*monitors(jm)%r(3,i)
          x=monitors(jm)%r(1,i)+x
          y=monitors(jm)%r(3,i)+y
       enddo

       x=x/n
       y=y/n
       write(6,*) "Averages ",x,y
       do i=1,n
          monitors(jm)%r(1,i)=(monitors(jm)%r(1,i)-x)/1000.0_dp
          monitors(jm)%r(3,i)=(monitors(jm)%r(3,i)-y)/1000.0_dp
       enddo

    endif  ! filename


    if(na>0) then
       xa=0.0_dp
       ya=0.0_dp

       do i=1,na
          xa(1)=monitors(jm)%r(1,i)**2+xa(1)
          ya(1)=monitors(jm)%r(3,i)**2+ya(1)
          xa(2)=monitors(jm)%r(1,n+1-i)**2+xa(2)
          ya(2)=monitors(jm)%r(3,n+1-i)**2+ya(2)
       enddo
       xa(2)=sqrt(xa(1)/xa(2))
       ya(2)=sqrt(ya(1)/ya(2))

       write(6,*) " scales ",xa(2),ya(2)
       xa(1)=(-1.0_dp+xa(2) )/(n-1)
       ya(1)=(-1.0_dp+ya(2) )/(n-1)
       xa(2)=(n-xa(2))/(n-1)
       ya(2)=(n-ya(2))/(n-1)

       do i=1,n
          monitors(jm)%r(1,i)=monitors(jm)%r(1,i) * (i*xa(1)+xa(2))
          monitors(jm)%r(3,i)=monitors(jm)%r(3,i) * (i*ya(1)+ya(2))
       enddo

       x=0.0_dp
       y=0.0_dp
       do i=1,n
          x=monitors(jm)%r(1,i)+x
          y=monitors(jm)%r(3,i)+y
       enddo

       x=x/n
       y=y/n
       write(6,*) "Averages ",x,y
       do i=1,n
          monitors(jm)%r(1,i)=(monitors(jm)%r(1,i)-x)
          monitors(jm)%r(3,i)=(monitors(jm)%r(3,i)-y)
       enddo

    endif

    close(mf1)

    ! 171683 59 68 1 2.6794 -1.1248


  end   subroutine  alex_read_r

  SUBROUTINE compute_twiss(r,my_state,filename,pos,del,no,thinlens,name,teng,notusingteng)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    integer pos,no,I,mf,ind,j
    TYPE(internal_state) state
    real(dp) closed(6),del,bx,by,ax,ay,s,px,py,ex,exp,dpx,dpy,phi,E1(6,6),E2(6,6)
    type(DAMAP) ID,a_f,a_l,a_nl,DR,R_TE,cs_te
    logical(lp) COSLIKE,thinlens, doname, printres
    TYPE(REAL_8) Y(6)
    CHARACTER(*) FILENAME
    TYPE(FIBRE), POINTER :: P
    TYPE(integration_node), POINTER :: t
    TYPE(NORMALFORM) NORM
    character(nlp) name
    character(9) typec
    logical(lp) printmap,teng,notusingteng
    printmap=my_false
    if(.not.associated(r%t)) then
       call make_node_layout(r)
    endif
    printres=my_false
    doname=my_true
    call context(name)
    if(index(name,'###')/=0) then
       doname=my_false
       printres=my_true
    endif
    typec=' '
    call kanalnummer(mf,filename)

    STATE=my_state+nocavity0
    call print(state,6)
   ! pause 777
    closed=0.0_dp
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,STATE,1e-5_dp)
    write(6,*) "closed orbit "
    write(6,*) CLOSED

    CALL INIT(STATE,no,0)
    CALL ALLOC(Y);   CALL ALLOC(NORM)
    call alloc(id,a_f,a_l,a_nl,DR,R_TE,cs_te)
    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,pos,STATE)
    id=y
    if(index(filename,'MAP')/=0.or.index(filename,'map')/=0) then
      printmap=my_true
    endif

    NORM=ID
    write(mf,*) " fractional tunes ",norm%tune(1:2)
    p=>r%start
    s=0.0_dp
    do i=1,pos-1
       s=s+p%mag%p%ld
       p=>p%next
    enddo
    px=0.0_dp
    py=0.0_dp
    if(teng) then
     call factor(norm%a_t,a_f,a_l,a_nl,DR=dr,r_te=r_te,cs_te=cs_te,COSLIKE=COSLIKE)
    else
     if(notusingteng) then
      call factor_am_special(norm%a_t,a_f,a_l,a_nl,DR=dr)
       else
      call factor(norm%a_t,a_f,a_l,a_nl,DR=dr)
     endif
    endif
    
    Y=CLOSED+norm%a_t
    !    bx,by,ay,ay,s,px,py
    bx=(norm%a_t%v(1).sub.'1')**2+(norm%a_t%v(1).sub.'01')**2
    by=(norm%a_t%v(3).sub.'001')**2+(norm%a_t%v(3).sub.'0001')**2
    ax=-(norm%a_t%v(1).sub.'1')*(norm%a_t%v(2).sub.'1')-(norm%a_t%v(1).sub.'01')*(norm%a_t%v(2).sub.'01')
    ay=-(norm%a_t%v(3).sub.'001')*(norm%a_t%v(4).sub.'001')-(norm%a_t%v(3).sub.'0001')*(norm%a_t%v(4).sub.'0001')
    ex=a_f%v(1).sub.'00001'
    exp=a_f%v(2).sub.'00001'
    if(teng) then
     phi=(r_te%v(1).sub.'1')-1.0_dp   !acos
    else
     phi=(norm%a_t%v(1).sub.'001')**2+(norm%a_t%v(1).sub.'0001')**2
    endif
    if(abs(phi)<epsflo) phi=0.0_dp
    if(teng) then
     if(COSLIKE) then
       typec=' = COS-1 '
     else
       typec=' = COSH-1'
     endif
    endif
    if(doname) then
       if(index(p%mag%name,name(1:len_trim(name)))/=0) then
          printres=my_true
       endif
    endif
    if(teng) then
    write(mf,'(6x,a4,6x,(6x,a1,6x),2(4x,a5,4x),2(4x,a6,3x),2(4x,a6,3x),2(3x,a7,3x),a25)') &
         "name", "s", "betax","betay","alphax","alphay","eta_x ","etap_x","Phase x","Phase y","Teng-Edwards Cos or Cosh "
    else
    write(mf,'(6x,a4,6x,(6x,a1,6x),2(4x,a5,4x),2(4x,a6,3x),2(4x,a6,3x),2(3x,a7,3x),a25)') &
         "name", "s", "betax","betay","alphax","alphay","eta_x ","etap_x","Phase x","Phase y","  d<x**2>/de_y (Ripken)  "
    endif
    if(printres) write(mf,'(a16,10(1x,g12.5),1x,a9)') p%mag%name(1:16),s,bx,by,ax,ay,ex,exp,px,py,phi,typec
    if(doname) printres=my_false

    do i=pos,pos+r%n-1
       t=>p%t1
       do while(.not.associated(t,p%t2%next))
          if(thinlens) then
             call track_probe_x(y,STATE, node1=t,node2=t%next)
          else
             call track_probe_x(y,STATE, node1=t,node2=p%t2%next)
          endif
          closed=y
          norm%a_t=y
        if(teng) then
         call factor(norm%a_t,a_f,a_l,a_nl,DR=dr,r_te=r_te,cs_te=cs_te,COSLIKE=COSLIKE)
        else
          if(notusingteng) then
           call factor_am_special(norm%a_t,a_f,a_l,a_nl,DR=dr)
            else
           call factor(norm%a_t,a_f,a_l,a_nl,DR=dr)
          endif        
        endif

!          call factor(norm%a_t,a_f,a_l,a_nl,DR=dr,r_te=r_te,cs_te=cs_te,COSLIKE=COSLIKE)
          px=px+asin(dr%v(1).sub.'01')/twopi
          py=py+asin(dr%v(3).sub.'0001')/twopi
          Y=CLOSED+norm%a_t

          if(thinlens) then
             t=>t%next
          else
             t=>p%t2%next
          endif
       enddo
       s=s+p%mag%p%ld
       p=>p%next

       bx=(norm%a_t%v(1).sub.'1')**2+(norm%a_t%v(1).sub.'01')**2
       by=(norm%a_t%v(3).sub.'001')**2+(norm%a_t%v(3).sub.'0001')**2
       ax=-(norm%a_t%v(1).sub.'1')*(norm%a_t%v(2).sub.'1')-(norm%a_t%v(1).sub.'01')*(norm%a_t%v(2).sub.'01')
       ay=-(norm%a_t%v(3).sub.'001')*(norm%a_t%v(4).sub.'001')-(norm%a_t%v(3).sub.'0001')*(norm%a_t%v(4).sub.'0001')
       ex=a_f%v(1).sub.'00001'
       exp=a_f%v(2).sub.'00001'
       if(teng) then
        phi=(r_te%v(1).sub.'1')-1.0_dp   !acos
       else
        phi=(norm%a_t%v(1).sub.'001')**2+(norm%a_t%v(1).sub.'0001')**2
       endif
       if(abs(phi)<epsflo) phi=0.0_dp
        if(teng) then
         if(COSLIKE) then
           typec=' = COS-1 '
         else
           typec=' = COSH-1'
         endif
        endif
       if(doname) then
          if(index(p%mag%name,name(1:len_trim(name)))/=0) then
             printres=my_true
          endif
       endif
       if(printres) write(mf,'(a16,10(1x,g12.5),1x,a9)') p%mag%name(1:16),s,bx,by,ax,ay,ex,exp,px,py,phi,typec
       if(doname) printres=my_false
    enddo
    if(printmap) then
        write(6,*) "stable_da, check_stable",stable_da, check_stable
        stable_da=my_true
        check_stable=my_true
        call print(id,mf)
    endif

    a_l=norm%a_t
    a_f=a_l**(-1)
    a_nl=0
    a_nl%v(1)=1.0_dp.mono.2
    a_nl%v(2)=-(1.0_dp.mono.1)
    a_nl%v(3)=1.0_dp.mono.4
    a_nl%v(4)=-(1.0_dp.mono.3)
    
    id=1
    id%v(3)=0
    id%v(4)=0    
    DR=A_l*a_nl*id*a_f*a_nl
    e1=dr
    e1=-e1
    id=1
    id%v(1)=0
    id%v(2)=0    
    DR=A_l*a_nl*id*a_f*a_nl
    e2=dr
    e2=-e2
    
  !! write(mf,'(6x,a4,6x,(6x,a1,6x),2(4x,a5,4x),2(4x,a6,3x),2(4x,a6,3x),2(3x,a7,3x),a25)')
          write(mf,*)  "   "
          write(mf,*)  "  All 4x4 Ripken Lattice Functions  "
          write(mf,*)  "   "
          write(mf,*)  " i , j, d<x_i x_j>/dJ_1  , d<x_i x_j>/dJ_2  "
    do i=1,4
    do j=i,4
        write(mf,'(1x,i2,2x,i2,2x,(1x,g12.5),(9x,g12.5))')  i,j,e1(i,j),e2(i,j)
    enddo
    enddo
          write(mf,*)  "   "
          write(mf,*)  "   4x4 Dispersions  "
          write(mf,*)  "   "
     do i=1,4
         e1(1,1)=a_l%v(i).sub.'00001'
        write(mf,'(1x,i2,(1x,g12.5))')  i,e1(1,1)
    enddo
    
    CLOSE(mf)
 
    CALL kill(Y)
    call kill(NORM)
    call kill(id,a_f,a_l,a_nl,DR,R_TE,cs_te)

  end SUBROUTINE compute_twiss
  
!!!! special rcs

  subroutine lattice_fit_bump_rcs(R,EPSF)
    IMPLICIT NONE
    TYPE(layout), target,intent(inout):: R
    integer, parameter :: neq=6,np=7
    real(dp)  TARG(neq)
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE, MF
    TYPE(TAYLOR) EQ(neq)
    TYPE(REAL_8) Y(6)
    integer ::  no=2,nt,j,it,pos1,pos2
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,gam(2)
    type(fibre), pointer:: p
    TYPE(POL_BLOCK) poly(np)
    
do i=1,np
 poly(i)=0
enddo

poly(1)%name='QDX2701'
poly(2)%name='QFL0101'
poly(3)%name='QDL0101'
poly(4)%name='QFM0201'
poly(5)%name='QDL0201'
poly(6)%name='QFL0301'
poly(7)%name='QDX0301'

do i=1,np
 poly(i)%ibn(2)=i
 r=poly(i)
enddo
!r%lastpos=1
!r%last=>R%start
call move_to( R,p,poly(1)%name,pos1,reset=my_true)
call move_to( R,p,poly(7)%name,pos2)
pos2=pos2+1

SET_TPSAFIT=.FALSE.

 
targ(1)=-1.216091585893243e0_dp     !  1 1
targ(2)=-6.061686871135153e0_dp     !  1 2
targ(3)=-0.7900090444471683E-01_dp  ! 2 1 
targ(4)=2.120414509984705e0_dp      ! 3 3
targ(5)=-19.02200463334060e0_dp     ! 3 4
targ(6)=-0.1837954391003473e0_dp    ! 4 3





    epsr=abs(epsf)

    nt=neq+np
    STATE=only_4d0

    CALL INIT(STATE,no,NP)


    it=0
100 continue
    it=it+1

    CLOSED=0.0_dp


    CALL INIT(STATE,no,NP)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,pos1,pos2,+STATE)

    !stop

 !   write(6,*) "c_%no,c_%nv,c_%nd,c_%nd2"
 !   write(6,*) c_%no,c_%nv,c_%nd,c_%nd2
 !   write(6,*) "c_%ndpt,c_%npara,c_%npara,c_%np_pol"
 !   write(6,*)  c_%ndpt,c_%npara,c_%npara,c_%np_pol

    eq(1)=(y(1)%t.par.'1000')-targ(1)
    eq(2)=(y(1)%t.par.'0100')-targ(2)
    eq(3)=(y(2)%t.par.'1000')-targ(3)
    eq(4)=(y(3)%t.par.'0010')-targ(4)
    eq(5)=(y(3)%t.par.'0001')-targ(5)
    eq(6)=(y(4)%t.par.'0010')-targ(6)

    epsnow=0.0_dp
    do i=1,neq
     epsnow=epsnow+abs(eq(i)) 
    enddo
     write(6,*) " deviation ",epsnow
     
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
       g%v(i)=1.0_dp.mono.i
       do j=np+1,nt
          t=g%v(j).d.i
          g%v(i)=g%v(i)+(1.0_dp.mono.j)*t
       enddo
    enddo
    CALL KILL(t)

    g=g.oo.(-1)
    tpsafit=0.0_dp
    tpsafit(1:nt)=g

    SET_TPSAFIT=.true.

    do i=1,7
     r=poly(i)
    enddo


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


  end subroutine lattice_fit_bump_rcs
  
    subroutine lattice_fit_bump_min_rcs(R0,R1,EPSF,poly,np,sca)
    IMPLICIT NONE
    TYPE(layout), target,intent(inout):: R0,R1
    integer, parameter :: neq=2
    real(dp)  TARG(neq),bety0m,bety1m,bety1,bety0,dbx,dby,betx1ma,bety1ma
    real(dp) X0(6),X1(6),CLOSED(6),betx0m,betx1m,betx1,betx0,dbyr,dbyrm,dbxr,dbxrm
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE, MF
    TYPE(TAYLOR) EQ(neq)
    TYPE(REAL_8) Y(6),Y1(6),Y0(6)
    integer ::  no=2,nt,j,it,pos1,pos2,np 
    type(damap) id1,ID0
    type(gmap) g
    TYPE(TAYLOR)t,bp(2),v
    type(taylor), allocatable :: dr(:)
    real(dp) epsf,epsr,epsnow,gam(2),sca
    type(fibre), pointer:: p,px,py
    TYPE(POL_BLOCK) poly(:)
    type(normalform) norm
    real(dp) s0,ds
 !   write(6,*)" correct ? "
 !   read(5,*) j
 !   if(j==1)    call lattice_fit_bump_rcs(R1,EPSF)
!do i=1,np
! poly(i)=0
!enddo

!poly(1)%name='QDX2701'
!poly(2)%name='QFL0101'
!poly(3)%name='QDL0101'
!poly(4)%name='QFM0201'
!poly(5)%name='QDL0201'
!poly(6)%name='QFL0301'
!poly(7)%name='QDX0301'

!do i=1,np
! poly(i)%ibn(2)=i
! r1=poly(i)
!enddo

    do i=1,np
     r1=poly(i)
    enddo

SET_TPSAFIT=.FALSE.
    STATE=only_4d0


    CALL FIND_ORBIT(R0,X0,1,STATE,1e-5_dp)
    write(6,*) "closed orbit ", CHECK_STABLE
    write(6,*) X0
    CALL FIND_ORBIT(R1,X1,1,STATE,1e-5_dp)
    write(6,*) "closed orbit ", CHECK_STABLE
    write(6,*) X1

   CALL INIT(STATE,1,0)
   
   CALL ALLOC(Y0)
   CALL ALLOC(Y1)
   CALL ALLOC(norm)
   CALL ALLOC(ID1,ID0)
   id1=1
   id0=1
   
   y0=x0+id0
   y1=x1+id1

    CALL TRACK(R0,Y0,1,STATE)
    CALL TRACK(R1,Y1,1,STATE)
    
    norm=y0
    id0=norm%a_t
    targ(1:2) = norm%tune(1:2)
    norm=y1
    id1=norm%a_t
    y0=x0+id0
    y1=x1+id1
  
    betx0m=0
    betx1m=0 
    bety0m=0
    bety1m=0 
    dbxrm=0
    dbyrm=0
   
    dbxr=0
    dbyr=0
   
    dbx=0
    dby=0
    s0=0
    ds=0
    
    p=>r0%start
    
    do i=1,r0%n
    
    ds=p%mag%p%ld
    
      call TRACK(R0,Y0,i,i+1,STATE)      
      betx0=(y0(1)%t.sub.'1')**2+(y0(1)%t.sub.'01')**2
      bety0=(y0(3)%t.sub.'001')**2+(y0(3)%t.sub.'0001')**2
        call TRACK(R1,Y1,i,i+1,STATE)          
      betx1=(y1(1)%t.sub.'1')**2+(y1(1)%t.sub.'01')**2
      bety1=(y1(3)%t.sub.'001')**2+(y1(3)%t.sub.'0001')**2
     if(betx0>betx0m) betx0m=betx0
     if(betx1>betx1m) betx1m=betx1   
     if(bety0>bety0m) bety0m=bety0
     if(bety1>bety1m) bety1m=bety1  
     dbxr= abs(betx1-betx0)  !/bety0
     dbyr= abs(bety1-bety0)  !/bety0
     if(dbxr>dbxrm) then
      dbxrm=dbxr
     ! write(6,*)"x", betx0,betx1
      px=>p
     endif
     if(dbyr>dbyrm) then
      dbyrm=dbyr
      !write(6,*)"y", bety0,bety1
      py=>p
     endif
     dbx=ds*(betx1-betx0)**2+dbx
     dby=ds*(bety1-bety0)**2+dby
     p=>p%next
     s0=s0+ds
    enddo
    
     dbx=dbx/s0
     dby=dby/s0
    write(6,*) " maximum "
    write(6,*) sqrt(betx0m),sqrt(betx1m)
    write(6,*) sqrt(bety0m),sqrt(bety1m)
 !   write(6,*) sqrt(dbx),sqrt(dby),s0
 !   write(6,*) px%mag%name,py%mag%name
 !   write(6,*) dbxr,dbyr
        
    
   CALL kill(Y0)
   CALL kill(Y1)
   CALL kill(norm)
   CALL kill(ID1,ID0)



    epsr=abs(epsf)

    nt=neq+np

    CALL INIT(STATE,no,NP)


    it=0
100 continue
    it=it+1

    x1=0.0_dp


    CALL INIT(STATE,no,NP)
    
    CALL ALLOC(Y0)
    CALL ALLOC(Y1)
    CALL ALLOC(EQ)
    call alloc(id0)
    call alloc(id1)
    call alloc(bp)
    call alloc(norm)
    call alloc(t,v)

    CALL FIND_ORBIT(R1,X1,1,STATE,1e-5_dp)

    id0=1
    Y0=x0+id0
    id1=1
    Y1=x1+id1

    CALL TRACK(R0,Y0,1,STATE)
    CALL TRACK(R1,Y1,1,+STATE)

    norm=y0
    id0=norm%a_t
    targ(1:2) = norm%tune(1:2)
    targ(1)=targ(1)+0.02d0
    targ(2)=targ(2)+0.01d0
    norm=y1
    id1=norm%a_t
     
    y0=x0+id0
    y1=x1+id1


    betx1ma=0
    bety1ma=0 

    s0=0
    p=>r0%start
    
    do i=1,r0%n
    
    ds=p%mag%p%ld
    
      call TRACK(R0,Y0,i,i+1,STATE)      
      betx0=(y0(1)%t.sub.'1')**2+(y0(1)%t.sub.'01')**2
      bety0=(y0(3)%t.sub.'001')**2+(y0(3)%t.sub.'0001')**2


      call TRACK(R1,Y1,i,i+1,+STATE)          
      bp(1)=ds*((y1(1)%t.par.'1000')**2+(y1(1)%t.par.'0100')**2-betx0)**2+bp(1)
      bp(2)=ds*((y1(3)%t.par.'0010')**2+(y1(3)%t.par.'0001')**2-bety0)**2+bp(2)
      
      betx1=(y1(1)%t.sub.'1000')**2+(y1(1)%t.sub.'0100')**2
      bety1=(y1(3)%t.sub.'0010')**2+(y1(3)%t.sub.'0001')**2
      if(betx1>betx1ma) betx1ma=betx1
      if(bety1>bety1ma) bety1ma=bety1



     p=>p%next
     s0=s0+ds
    enddo

    write(6,*) " maximum "
    write(6,*) sqrt(betx0m),sqrt(betx1m)
    write(6,*) sqrt(bety0m),sqrt(bety1m)
    write(6,*) " maximum now"
    write(6,*) sqrt(betx0m),sqrt(betx1ma)
    write(6,*) sqrt(bety0m),sqrt(bety1ma)

    !stop

 !   write(6,*) "c_%no,c_%nv,c_%nd,c_%nd2"
 !   write(6,*) c_%no,c_%nv,c_%nd,c_%nd2
 !   write(6,*) "c_%ndpt,c_%npara,c_%npara,c_%np_pol"
 !   write(6,*)  c_%ndpt,c_%npara,c_%npara,c_%np_pol

    eq(1)=       ((NORM%dhdj%v(1)).par.'0000')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'0000')-targ(2)
    bp(1)=bp(1)/s0
    bp(2)=bp(2)/s0
    epsnow=0.0_dp
    do i=1,neq
     epsnow=epsnow+abs(eq(i)) 
    enddo
    
     write(6,*) " deviation ",epsnow
!write(6,*) " scale "
!read(5,*) sca
    do i=1,neq
     eq(i)=eq(i)-(1.0_dp-sca)*(eq(i).sub.'0')
    enddo
    epsnow=0.0_dp
    do i=1,neq
     epsnow=epsnow+abs(eq(i)) 
    enddo
      write(6,*) " deviation ",epsnow
    
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    
     t=(bp(1)+bp(2))
         ds=t
     write(6,*) " Merit function  ",ds

    ! t=0
    ! do i=1,np
    !  t=t+half*(one.mono.(i+c_%nd2))**2
    ! enddo
     
    do i=1,np
     v=t.d.(i+c_%nd2)
     v=v<=c_%npara
       call daprint(v,scratchfile)  
    enddo
 
    do i=1,neq
        eq(i)=eq(i)<=c_%npara
      call daprint(eq(i),scratchfile)
    enddo
    close(SCRATCHFILE)

    CALL KILL(Y0)
    CALL KILL(Y1)
    call KILL(id0)
    CALL KILL(id1)
    CALL KILL(EQ)
    call KILL(bp)
    call KILL(norm)
    call KILL(t,v)


     
    CALL INIT(1,nt)
    
    call alloc(g,nt)
    call alloc(t,v)
    allocate(dr(nt))
    call alloc(dr,nt)
    
    call kanalnummer(SCRATCHFILE)
    OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
    rewind scratchfile
    do i=1,nt
       call read(dr(i),scratchfile)
    enddo

    do i=1,np
      t=0.0_dp
      do j=1,neq
       v=(dr(np+j).d.i)
       t=(1.0_dp.mono.(np+j))*v+t
      enddo 
       g%v(i)=dr(i)+t
    enddo
        do j=1,neq
         g%v(np+j)=dr(np+j)
       enddo  
    close(SCRATCHFILE)


    g=g.oo.(-1)
    tpsafit=0.0_dp
    tpsafit(1:nt)=g

    SET_TPSAFIT=.true.

    do i=1,np
     r1=poly(i)
    enddo


    SET_TPSAFIT=.false.

    CALL ELP_TO_EL(R1)
    CALL KILL(t,v)
    CALL KILL(g)
    call KILL(dr,nt)
    deallocate(dr)

     !   write(6,*) " more "
     !   read(5,*) i
     !   if(i==0) goto 102
    if(it>=max_fit_iter/sca**2) goto 101
    if(epsnow<=epsr) goto 102
    GOTO 100

101 continue
    write(6,*) " warning did not converge "

102 continue
    CALL KILL_PARA(R1)




  end subroutine lattice_fit_bump_min_rcs
  

  
end module S_fitting_new
