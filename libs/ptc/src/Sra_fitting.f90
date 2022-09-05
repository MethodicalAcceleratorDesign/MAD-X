!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_fitting_new
  USE ptc_spin
  IMPLICIT NONE
  public
  private FIND_ORBIT_LAYOUT_noda,FIND_ORBIT_LAYOUT_noda_object
  private FIND_ORBIT_LAYOUT_da,FIND_ORBIT_LAYOUT_da_object
  integer:: m_turn,m_skip=0
  integer :: with_c=1
  integer ::  other_fix=0
  logical :: check_longitudinal=.false.,piotr_fix=.true.

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

    INTERFACE FIND_ORBIT_TPSA_x
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_da
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_da_object
  END INTERFACE

  INTERFACE FIND_ORBIT_probe_tpsa_x
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_da
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_da_object
  END INTERFACE

  INTERFACE FIND_ORBIT_x
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda_object
  END INTERFACE

  INTERFACE FIND_ORBIT_probe_x
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda_object
  END INTERFACE
  
  
contains

  subroutine lattice_fit_TUNE_gmap_rad(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout), target,intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,NP
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE, MF
    TYPE(TAYLOR), allocatable:: EQ(:)
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer :: neq=2, no=2,nt,j,it
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,gam(2)
    !    EPSF=.0001
    epsr=abs(epsf)

    allocate(eq(neq))

    nt=neq+np
    STATE=time0+radiation0

    CALL INIT(STATE,no,NP)

    SET_TPSAFIT=.FALSE.


    DO I=1,NPOLY
       R=POLY(i)
    ENDDO

    CLOSED(:)=0.0_dp
    it=0
100 continue
    it=it+1
      call FIND_ORBIT_x(r,CLOSED,state,1.0e-7_dp,fibre1=1)
 
    write(6,*) "closed orbit "
    write(6,*) CLOSED


    CALL INIT(STATE,no,NP,BERZ)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK_probe_x(R,Y,+STATE,fibre1=1)
    write(6,*) "c_%no,c_%nv,c_%nd,c_%nd2"
    write(6,*) c_%no,c_%nv,c_%nd,c_%nd2
    write(6,*) "c_%ndpt,c_%npara,c_%npara,c_%np_pol"
    write(6,*)  c_%ndpt,c_%npara,c_%npara,c_%np_pol

     id=y
    NORM=id
    gam(1)=(norm%a_t%v(2).sub.'1')**2+(norm%a_t%v(2).sub.'01')**2
    gam(2)=(norm%a_t%v(4).sub.'001')**2+(norm%a_t%v(4).sub.'0001')**2
    write(6,*) "  Gamma= ",GAM
    !      CALL KANALNUMMER(MF)
   ! OPEN(UNIT=1111,FILE='GAMMA.TXT')
   ! WRITE(1111,*) "  Gamma= ",GAM

    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2)

    eq(1)=       ((NORM%dhdj%v(1)).par.'000000')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'000000')-targ(2)
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

  end subroutine lattice_fit_TUNE_gmap_rad
 !!! my_default default0 (my units)
subroutine find_time_patch(kekb,my_default,emax,bmadpatch,wipeout,kf,kb)
! time patch for SAD
implicit none
type(layout), pointer :: kekb
real(dp) closed_orbit(6),ee
integer, optional :: kf,kb
integer kc,ke,i,cc,jc
type(internal_state) my_default,state
type(fibre), pointer :: f,f1,f2
logical, optional :: bmadpatch,wipeout
logical bm,wi
real(dp), optional :: emax

ee=1.d40
if(present(emax)) ee=emax

bm=.true.
wi=.false.
if(present(bmadpatch)) bm =bmadpatch
if(present(wipeout)) wi =wipeout
ke=0
kc=0
if(wi) then
 f=> kekb%start
 do i=1,kekb%n
  f%patch%time=0
  f%patch%a_t=0
  f%patch%b_t=0
  f%patch%a_L=0
  f%patch%b_L=0
 f=>f%next
 enddo
endif


state=(my_default+nocavity0)-radiation0-spin0-time0
closed_orbit=0.d0

call find_orbit_x(kekb,closed_orbit,STATE,1.e-8_dp,fibre1=1)  

cc=0
!call propagate(kekb,closed_orbit,state,fibre1=1)
closed_orbit(6)=0.d0

f1=> kekb%start
if(bm) then
f=> kekb%start
 do i=1,kekb%n
   if(f%mag%kind==kind4) then 
     cc=cc+1
     if(cc==1) then 
      f1=>f
      jc=i
     endif
     f2=>f
   endif
   f=>f%next
 enddo
endif


f=> kekb%start
do i=1,kekb%n

call propagate(kekb,closed_orbit,state,fibre1=i,fibre2=i+1)


if(bm) then

 if(f%next%mag%kind==kind4) then
  f%patch%time=2
  f%patch%B_L=closed_orbit(6)+f%patch%B_L
  f%patch%B_T=closed_orbit(6)/F%beta0+f%patch%B_T
 ke=ke+1
 elseif(f%mag%kind==kind4) then
  if(associated(f,f2)) then
   f%next%patch%time=1
   f%patch%A_L=closed_orbit(6)+f%patch%A_L
   f%patch%A_T=closed_orbit(6)/F%beta0+f%patch%A_T
  kc=kc+1
  endif
 closed_orbit(6)=0.d0
 endif

else

if(abs(closed_orbit(6))>ee.or.f%next%mag%kind==kind4.or.f%mag%kind==kind4) then

 if(f%next%mag%kind==kind4) then
  f%next%patch%time=1
  f%next%patch%A_L=closed_orbit(6)+f%next%patch%A_L
  f%next%patch%A_T=closed_orbit(6)/F%beta0+f%next%patch%A_T
 kc=kc+1
 elseif(f%mag%kind==kind4) then
  f%patch%time=3
  f%patch%B_L=closed_orbit(6)+f%patch%B_L
  f%patch%B_T=closed_orbit(6)/F%beta0+f%patch%B_T
  ke=ke+1
 closed_orbit(6)=0.d0
 else
  f%patch%time=2
  f%patch%B_L=closed_orbit(6)+f%patch%B_L
  f%patch%B_T=closed_orbit(6)/F%beta0+f%patch%B_T
  ke=ke+1
  closed_orbit(6)=0.d0
 endif

endif

endif
 

f=>f%next
enddo


if(bm) then
  f2%next%patch%A_L=closed_orbit(6)+f2%next%patch%A_L
  f2%next%patch%A_T=closed_orbit(6)/F%beta0+f2%next%patch%A_T
else
f=> kekb%end
  f%patch%time=2
  f%patch%B_L=closed_orbit(6)+f%patch%B_L
  f%patch%B_T=closed_orbit(6)/F%beta0+f%patch%B_T
  ke=ke+1
endif
 
if(present(kb)) kb=ke
if(present(kf)) kf=kc

end subroutine find_time_patch

  subroutine compute_linear_one_magnet_maps(f,state,del)
    implicit none
    TYPE(fibre), pointer, intent(inout):: f
    TYPE(layout), pointer :: als
    type(internal_state) state
    integer i,no
    logical rad
    real(dp), optional :: del
    real(dp) closed(6),m(6,6) 
    type(c_damap) c_map,d_map,id_s
    type(probe) xs0
    type(probe_8) xs
    TYPE(fibre), pointer :: f1
    rad=state%radiation
    closed=0.d0
    no=1
    als=>f%parent_layout

    if(present(del)) closed(5+ndpt_bmad)=del
    call find_orbit_x(CLOSED,STATE,1.e-8_dp,fibre1=f)       
    call init_all(STATE,no,0)
call alloc(c_map,d_map)
call alloc(id_s)
 
   m=0.0_dp
 
!!!! Polymorphic probe is created in the usual manner 





! Copy probe_8 into a complex damap 


f1=>f
do i=1,als%n
   XS0=CLOSED    
   ID_S=1        
   XS=XS0+ID_S 
f1%i%fix0=xs%x
 CALL propagate(XS,STATE,FIBRE1=f1,fibre2=f1%next)
f1%i%fix=xs%x
d_map=xs
f1%i%m=d_map
CLOSED=xs%x

f1=>f1%next

enddo



call kill(c_map,d_map)
call kill(id_s)

end subroutine compute_linear_one_magnet_maps

  subroutine compute_linear_one_turn_maps(f,state,del)
    implicit none
    TYPE(fibre), pointer, intent(inout):: f
    TYPE(layout), pointer :: als
    type(internal_state) state
    integer i,no
    logical rad
    real(dp), optional :: del
    real(dp) closed(6),m(6,6),s(6,6)
    type(c_normal_form) c_n
    type(c_damap) c_map,d_map,id_s
    type(probe) xs0
    type(probe_8) xs
    TYPE(fibre), pointer :: f1
    rad=state%radiation
    closed=0.d0
    no=1
    als=>f%parent_layout

    if(present(del)) closed(5+ndpt_bmad)=del
    call find_orbit_x(CLOSED,STATE,1.e-8_dp,fibre1=f)    
    call init_all(STATE,no,0)
call alloc(c_map,d_map)
call alloc(c_n)
call alloc(id_s)
   s=0.0_dp
   m=0.0_dp
do i=1,3
 S(2*i-1,2*i) = 1.d0; S(2*i,2*i-1) = -1.d0;
enddo
!!!! Polymorphic probe is created in the usual manner 
   XS0=CLOSED    
   ID_S=1        
   XS=XS0+ID_S 

!!!! get spin polymorphic probe after one turn   
CALL propagate(XS,STATE,FIBRE1=f)  ! (4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Copy probe_8 into a complex damap 
c_map=XS ! (5)
m=c_map
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
!call c_normal(c_map,c_n,dospin=state%spin)  ! (6)
! write(6,'(4(1x,g21.14))') c_n%tune(1:3), c_n%spin_tune
f%i%m=m
f%i%fix0=xs%x

   XS0=CLOSED    
   ID_S=1        
   XS=XS0+ID_S 

f1=>f
do i=1,als%n
 
 CALL propagate(XS,STATE,FIBRE1=f1,fibre2=f1%next)
f1=>f1%next
d_map=xs
f1%i%fix0=xs%x
if(.not.state%radiation) then
f1%i%m=d_map
f1%i%m=matmul(f1%i%m, matmul(m,matmul(S,matmul(transpose(f1%i%m),transpose(S)))))
else
f1%i%m=d_map*c_map*d_map**(-1)
endif

enddo



call kill(c_map,d_map)
call kill(c_n)
call kill(id_s)

end subroutine compute_linear_one_turn_maps

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
  
      SUBROUTINE FIND_ORBIT_LAYOUT_da(RING,FIX,STATE,TURNS,fibre1,node1,total) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    real(dp) , intent(inOUT) :: FIX(6)
    real(dp), optional :: total
    INTEGER , optional,intent(in) :: TURNS,node1,fibre1
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    type(fibre), pointer :: object_fibre1
    type(integration_node), pointer :: object_node1
    integer i

    if(present(fibre1)) then
       object_fibre1=>ring%start
      do i=1,fibre1-1
         object_fibre1=>object_fibre1%next
      enddo   
      call FIND_ORBIT_LAYOUT_da_object(FIX,STATE,TURNS,fibre1=object_fibre1,total=total)
     else
       object_node1=>ring%t%start
      do i=1,node1-1
         object_node1=>object_node1%next
      enddo 
      call FIND_ORBIT_LAYOUT_da_object(FIX,STATE,TURNS,node1=object_node1,total=total)
     endif


     end  SUBROUTINE FIND_ORBIT_LAYOUT_da
 
    SUBROUTINE FIND_ORBIT_LAYOUT_da_object(FIX0,STATE,TURNS,fibre1,node1,total) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    real(dp) , intent(inOUT) :: FIX0(6)
    real(dp), optional :: total
    INTEGER , optional,intent(in) :: TURNS 
    type(fibre), optional, pointer :: fibre1
    type(integration_node), optional, pointer :: node1
    real(dp)  eps,TOT,freq
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat
    real(dp)  DIX(6),xdix,xdix0,tiny,beta1,freqmin,t6
    real(dp) X(6),Y(6),MX(6,6),sxi(6,6),SX(6,6),fix(6)
    integer NO1,ND2,I,IU,ITE,ier,j,ITEM
    TYPE (fibre), POINTER :: C
    TYPE (integration_node), POINTER :: t
    logical(lp) APERTURE,use_bmad_units_temp
    INTEGER TURNS0,trackflag
    type(damap) id
    type(real_8) yy(6)
    type(layout), pointer :: ring
    logical isStableFixPoint

    fix=fix0
    tot=0
    if(present(fibre1)) then
     ring=>fibre1%parent_layout
    else
     ring=>node1%parent_fibre%parent_layout
    endif

    if (.not.STATE%NOCAVITY.and.check_longitudinal) then
        freqmin=1.d38
        t6=0
        c=>ring%start    
      do i=1,ring%n
        if(c%mag%kind==kind4.or.c%mag%kind==kind21) then
          if(abs(c%mag%freq)>0) then
            if( abs(c%mag%freq)<=freqmin) freqmin=c%mag%freq
          endif
        endif
         t6=t6+c%mag%p%ld
       c=>c%next
      enddo 
      if(state%time) t6=t6/c%beta0
     if(global_verbose) write(6,*) " Harmonic # = ", freqmin*t6/clight
     if(global_verbose) write(6,*) " freqmin , dt", freqmin,clight/freqmin/2
    endif


    if(.not.associated(RING%t)) call MAKE_NODE_LAYOUT(ring)
    !!    xs%x=zero
    !!    xs%s%x=zero
    use_bmad_units_temp=use_bmad_units
 
    if(use_bmad_units_temp) then 
          if(present(fibre1)) then
           beta1=fibre1%mag%p%beta0
          else
           beta1=node1%parent_fibre%mag%p%beta0
          endif
      call convert_bmad_to_ptc(fix,beta1,STATE%TIME)
    endif

    use_bmad_units=.false.
    TURNS0=1
    trackflag=0
    tot=0
    IF(PRESENT(TURNS)) TURNS0=TURNS
 
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    !!    call move_to(ring,c,loc)
    !!    loct=c%t1%pos


    Nullify(C);

    if(.not.ring%closed) then
       !w_p=0
       !w_p%nc=1
       !w_p%fc='((1X,a72))'
       !w_p%c(1)=" This line is not ring : FIND_ORBIT_LAYOUT_noda "
       ! call !write_e(100)
    endif
    dix(:)=0.0_dp
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking
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
          
          messagelost= " No Cavity in the Line "
          check_stable=.false.
          return
 
       ENDIF
    else   ! (.not.present(STATE)) t
       IF(STATE%NOCAVITY) THEN
          ND2=4
          STAT=STATE+only_4d0
          if(state%radiation) then
             check_stable=.false.

             messagelost= " Cavity needed when radiation present "
             return
          endif
       ELSE
          ND2=6
          STAT=STATE
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          check_stable=.false.
          messagelost= " State present; no cavity: FIND_ORBIT_LAYOUT will crash => exiting"
         return

       ENDIF
    endif
101 continue



    if((stat%totalpath==1).and.(.not.stat%nocavity)) then
       C=>RING%START
       freq=0.0_dp
       i=1
       xdix=0.0_dp
       do while(i<=RING%n)
          if(associated(c%mag%freq)) then
             IF(FREQ==0.0_dp) THEN
                freq=c%mag%freq
             ELSEIF(c%mag%freq/=0.0_dp.AND.c%mag%freq<FREQ) THEN
                freq=c%mag%freq
             ENDIF
          endif
          c=>c%next
          i=i+1
       enddo
       
       if(ring%harmonic_number==0) then
          t=>ring%t%start
          tot=0
          
          do i=1,ring%t%n
             tot= tot + t%ds_ac
             t=>t%next
          enddo
       
          if(state%time) tot=tot/c%beta0

          if(freq_redefine) then
             ring%harmonic_number=tot*freq/twopi
          else
             ring%harmonic_number=tot*freq/CLIGHT
          endif

       else

          if(freq_redefine) then
             tot=RING%HARMONIC_NUMBER*twopi/FREQ
          else
             tot=RING%HARMONIC_NUMBER*CLIGHT/FREQ
          endif
       endif
      
   endif
    



call init(stat,1,0)
call alloc(yy); call alloc(id);

1111 continue

    ITEM=0
3   continue
    ITEM=ITEM+1
!    X=FIX
    id=1
    yy=fix+id
    DO I=1,TURNS0
       !       CALL TRACK(RING,X,LOC,STAT)
       !       trackflag=TRACK_flag(RING,X,LOC,STAT)
       !!       xs%x=x

       call TRACK_probe_X(yy,stat,fibre1=fibre1,node1=node1)
 
       if(.not.check_stable) then
          messagelost(len_trim(messagelost)+1:255)=" -> Unstable tracking guessed orbit "
          c_%APERTURE_FLAG=APERTURE
          return
       endif
       !     write(6,*) item,check_stable
       !!       call TRACK_PROBE(Ring,xs,loct,loct+ring%t%n,stat)
       !!       x=xs%x
       !       if(trackflag/=0) then
       !         ITEM=MAX_FIND_ITER+100
       !       endif

    ENDDO
    !    write(6,*) x
  !  x(6)=x(6)-freq*turns0
    id=yy
    mx=0.0_dp
    mx=id
    do i=1,nd2
     x(i)=yy(i)
    enddo
 

    SX=MX;
    DO I=1,nd2   !  6 before
       SX(I,I)=MX(I,I)-1.0_dp
    ENDDO

    DO I=1,ND2
       DIX(I)=FIX(I)-X(I)
      if(i==6) dix(6)=dix(6)+tot
    enddo
    
    CALL matinv(SX,SXI,ND2,6,ier)
    IF(IER==132)  then
       messagelost= " Inversion failed in FIND_ORBIT_LAYOUT_da"
        check_stable=.false.
       return
    endif

    x=0.0_dp
    do i=1,nd2
       do j=1,nd2
          x(i)=sxi(i,j)*dix(j)+x(i)
       enddo
    enddo
    dix=x
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=0.0_dp
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo
    !    write(6,*) " Convergence Factor = ",nd2,xdix,deps_tracking
    !    pause 123321
  !  if(verbose) write(6,*) " Convergence Factor = ",xdix
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
       !   C_%stable_da=.FALSE.
       !      IF(iteM==MAX_FIND_ITER+100) THEN
       !        write(6,*) " Unstable in find_orbit without TPSA"
       messagelost= "Maximum number of iterations in find_orbit with TPSA"
       xlost=fix
       check_stable=my_false
       !     ENDIF
       ITE=0
    endif
    !   write(6,*) item,xdix,xdix0,tiny

    if(ite.eq.1.or.item<item_min)  then

       GOTO 3

    endif


    if (ND2 == 6.and.check_longitudinal) then
       eps=1e-8_dp
      isStableFixPoint = is_ORBIT_STABLE(FIX,EPS,STAT,fibre1,node1)
      if (isStableFixPoint .eqv. .false.) then
        
     if(global_verbose) print*,"Orbit seemed to be unstable in longitudinal"
        fix = fix0
        fix(6)= fix(6)+ clight/freqmin/2

        goto 1111
        
      endif
    endif

call kill(yy); call kill(id);

    if(use_bmad_units_temp) then 
 
      call convert_ptc_to_bmad(fix,beta1,STATE%TIME)
    endif
   use_bmad_units=use_bmad_units_temp
    !    FIX(6)=FIX(6)+freq*turns0
    c_%APERTURE_FLAG=APERTURE
    fix0=fix
    if(present(total)) total=tot
  END SUBROUTINE FIND_ORBIT_LAYOUT_da_object



  SUBROUTINE FIND_ORBIT_LAYOUT_noda(RING,FIX,STATE,eps,TURNS,fibre1,node1,total) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    real(dp) , intent(inOUT) :: FIX(6)
    real(dp), optional :: total
    INTEGER , optional,intent(in) :: TURNS,node1,fibre1
    real(dp)  eps,TOT,freq,t6
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    type(fibre), pointer :: object_fibre1
    type(integration_node), pointer :: object_node1
    integer i

    if(present(fibre1)) then
       object_fibre1=>ring%start
      do i=1,fibre1-1
         object_fibre1=>object_fibre1%next
      enddo   
      if(piotr_fix) then
      call FIND_ORBIT_LAYOUT_noda_object(FIX,STATE,eps,TURNS,fibre1=object_fibre1,total=total)
       else
      call FIND_ORBIT_LAYOUT_noda_object_orig(FIX,STATE,eps,TURNS,fibre1=object_fibre1,total=total)
      endif
     else
       object_node1=>ring%t%start
      do i=1,node1-1
         object_node1=>object_node1%next
      enddo 
      if(piotr_fix) then
       call FIND_ORBIT_LAYOUT_noda_object(FIX,STATE,eps,TURNS,node1=object_node1,total=total)
       else
       call FIND_ORBIT_LAYOUT_noda_object_orig(FIX,STATE,eps,TURNS,node1=object_node1,total=total)
      endif
     endif

  end SUBROUTINE FIND_ORBIT_LAYOUT_noda


  SUBROUTINE FIND_ORBIT_LAYOUT_noda_object(FIX0,STATE,eps,TURNS,fibre1,node1,total) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),pointer :: RING
    real(dp) , intent(inOUT) :: FIX0(6)
    real(dp), optional :: total
    INTEGER , optional,intent(in) :: TURNS
    type(fibre), optional, pointer :: fibre1
    type(integration_node), optional, pointer :: node1
    real(dp)  eps,freq,tot
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat

    real(dp)  DIX(6),xdix,xdix0,tiny,beta1,t6,freqmin
    real(dp) X(6),Y(6),MX(6,6),sxi(6,6),SX(6,6),fix(6), Y6start
    integer NO1,ND2,I,IU,ITE,ier,j,ITEM !,try
    TYPE (fibre), POINTER :: C
    TYPE (integration_node), POINTER :: t
    logical(lp) APERTURE,use_bmad_units_temp
    logical isStableFixPoint, didPhaseJump
    INTEGER TURNS0,trackflag


    fix=fix0

    didPhaseJump = .false.

    tot=0
    if(present(fibre1)) then
     ring=>fibre1%parent_layout
    else
     ring=>node1%parent_fibre%parent_layout
    endif


    if (.not.STATE%NOCAVITY.and.check_longitudinal) then
        freqmin=1.d38
        t6=0
        c=>ring%start    
      do i=1,ring%n
        if(c%mag%kind==kind4.or.c%mag%kind==kind21) then
          if(abs(c%mag%freq)>0) then
            if( abs(c%mag%freq)<=freqmin) freqmin=c%mag%freq
          endif
        endif
         t6=t6+c%mag%p%ld
       c=>c%next
      enddo 
      if(state%time) t6=t6/c%beta0
     if(global_verbose) write(6,*) " Harmonic # = ", freqmin*t6/clight
     if(global_verbose) write(6,*) " freqmin , dt", freqmin,clight/freqmin/2
    endif

    !    fixed_found=my_true
    !!    type(probe) xs
    if(.not.associated(RING%t)) call MAKE_NODE_LAYOUT(ring)
    !!    xs%x=zero
    !!    xs%s%x=zero
    use_bmad_units_temp=use_bmad_units
    if(use_bmad_units_temp) then 
          if(present(fibre1)) then
           beta1=fibre1%beta0
  !         beta1=fibre1%mag%p%beta0

          else
!           beta1=node1%parent_fibre%mag%p%beta0
           beta1=node1%parent_fibre%beta0

          endif
      call convert_bmad_to_ptc(fix,beta1,STATE%TIME)
    endif
    
    use_bmad_units=.false.

    TURNS0=1
    trackflag=0
    IF(PRESENT(TURNS)) TURNS0=TURNS
 
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    !!    call move_to(ring,c,loc)
    !!    loct=c%t1%pos


    Nullify(C);

 
    dix(:)=0.0_dp
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking
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
          
          messagelost= " No Cavity in the Line "
          check_stable=.false.
            use_bmad_units=use_bmad_units_temp
          return
 
       ENDIF
    else   ! (.not.present(STATE)) t
       IF(STATE%NOCAVITY) THEN
          ND2=4
          STAT=STATE+only_4d0
          if(state%radiation) then
             check_stable=.false.

             messagelost= " Cavity needed when radiation present "
            use_bmad_units=use_bmad_units_temp
             return
          endif
       ELSE
          ND2=6
          STAT=STATE
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          check_stable=.false.
          messagelost= " State present; no cavity: FIND_ORBIT_LAYOUT will crash => exiting"
            use_bmad_units=use_bmad_units_temp
         return

       ENDIF
    endif
101 continue


    if((stat%totalpath==1).and.(.not.stat%nocavity)) then
       C=>RING%START
       freq=0.0_dp
       i=1
       xdix=0.0_dp
       do while(i<=RING%n)
          if(associated(c%magp%freq)) then
             IF(FREQ==0.0_dp) THEN
                freq=c%magp%freq
             ELSEIF(c%magp%freq/=0.0_dp.AND.c%magp%freq<FREQ) THEN
                freq=c%magp%freq
             ENDIF
          endif
          c=>c%next
          i=i+1
       enddo
       
       
       if(global_verbose) write(6,*) " Using frequency ", freq

       
       if(ring%harmonic_number==0) then
          t=>ring%t%start
          tot=0
          do i=1,ring%t%n
             tot= tot + t%ds_ac
             t=>t%next
          enddo

          if(state%time) tot=tot/c%beta0

          if(global_verbose) write(6,*) " Integrated tot of the machine ", tot

          if(freq_redefine) then
             ring%harmonic_number=tot*freq/twopi
          else
             ring%harmonic_number=tot*freq/CLIGHT
          endif

          if(global_verbose) write(6,*) " Corresponding harmonic number ", ring%harmonic_number

       else

         if(freq_redefine) then
           tot=RING%HARMONIC_NUMBER*twopi/FREQ
         else
           tot=RING%HARMONIC_NUMBER*CLIGHT/FREQ
         endif

         if(global_verbose) write(6,*) " Harmonic number already defined", ring%harmonic_number
         if(global_verbose) write(6,*) " Corresponding tot ", tot
         
      endif

   endif
    



    
1111 continue


    ITEM=0
3   continue
    ITEM=ITEM+1
    X=FIX

    DO I=1,TURNS0
       call TRACK_probe_X(x,stat,fibre1=fibre1,node1=node1)

       if(.not.check_stable) then
          messagelost(len_trim(messagelost)+1:255)=" -> Unstable tracking guessed orbit "
          c_%APERTURE_FLAG=APERTURE
 !                if(try>0) goto 1111
            use_bmad_units=use_bmad_units_temp
          return
       endif

    ENDDO

    if(global_verbose) then
      write(6,*) "#############################################"
      write(6,*) "ITERATION ", ITEM
      write(6,*) ""
      write(6,*) "FIX start :", FIX(:)
      write(6,*) "FIX end   :", X(:)
      write(6,*) "DIFF      :", FIX - X 
    endif
 


    mx=0.0_dp
    DO J=1,ND2
       Y=FIX
       Y(J)=FIX(J)+EPS
       Y6start=Y(6)
       DO I=1,TURNS0
          call TRACK_probe_X(Y,stat,fibre1=fibre1,node1=node1)

          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
 !                   if(try>0) goto 1111
            use_bmad_units=use_bmad_units_temp
             return
          endif

          
       ENDDO

       if(stat%totalpath==1) then
         y(6)=y(6)-TURNS0*gettot((y(6)/TURNS0) - Y6start , freq)
       endif
       
       do i=1,ND2
            MX(I,J)=Y(i)/2/eps+MX(I,J)   
       enddo
       
       Y=FIX
       Y(J)=FIX(J)-EPS
       Y6start=Y(6)
       DO I=1,TURNS0
          !          CALL TRACK(RING,Y,LOC,STAT)
          !!       xs%x=y
          call TRACK_probe_X(Y,stat,fibre1=fibre1,node1=node1)

          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
 !                   if(try>0) goto 1111
            use_bmad_units=use_bmad_units_temp
             return
          endif
 
       ENDDO
 
      if(stat%totalpath==1) then
         y(6)=y(6)-TURNS0*gettot(y(6)/TURNS0 - Y6start, freq)
       endif

       do i=1,ND2
            MX(I,J)=-Y(i)/2/eps+MX(I,J)   
       enddo

       
    ENDDO

    if(global_verbose) then
      write(6,*) ""
      write(6,*) "MX"
      do i=1,ND2
         write(6,*) "    ",   MX(I,:)
      enddo
    endif

    SX=MX;
    DO I=1,nd2   !  6 before
       SX(I,I)=MX(I,I)-1.0_dp
    ENDDO

    if(global_verbose) then
      write(6,*) ""
      write(6,*) "SX"
      do i=1,ND2
         write(6,*) "    ",   SX(I,:)
      enddo
    endif

    DO I=1,ND2
       DIX(I)=FIX(I)-X(I)
       if(i==6 .and. stat%totalpath==1) dix(6)=dix(6)+gettot(X(6) - FIX(6),freq)
    enddo
    
    CALL matinv(SX,SXI,ND2,6,ier)
    IF(IER==132)  then
       messagelost= " Inversion failed in FIND_ORBIT_LAYOUT_noda"
        check_stable=.false.
            use_bmad_units=use_bmad_units_temp
       return
    endif

    x=0.0_dp
    do i=1,nd2
       do j=1,nd2
          x(i)=sxi(i,j)*dix(j)+x(i)
       enddo
    enddo
    dix=x
    
    
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=0.0_dp
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo

    if(global_verbose) then
      write(6,*) ""
      write(6,*) "DIX (corr): ", DIX
      write(6,*) "PENALTY F : ", XDIX
    endif
    
    
    !    write(6,*) " Convergence Factor = ",nd2,xdix,deps_tracking
    !    pause 123321
  !  if(verbose) write(6,*) " Convergence Factor = ",xdix
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
       !   C_%stable_da=.FALSE.
       !      IF(iteM==MAX_FIND_ITER+100) THEN
       !        write(6,*) " Unstable in find_orbit without TPSA"
       messagelost= "Maximum number of iterations in find_orbit without TPSA"
       xlost=fix
       check_stable=my_false
       !     ENDIF
       ITE=0
   !           if(try>0) goto 1111
            use_bmad_units=use_bmad_units_temp
      return
    endif
    !   write(6,*) item,xdix,xdix0,tiny

    if(ite.eq.1.or.item<item_min)  then

       GOTO 3

    endif

    if (ND2 == 6.and.check_longitudinal) then
      isStableFixPoint = is_ORBIT_STABLE(FIX,EPS,STAT,fibre1,node1)

      if (isStableFixPoint .eqv. .false.) then
        if (didPhaseJump ) then

          messagelost= "Found unstable fixed point"
          xlost=fix
          check_stable=my_false

        else
        
         if(global_verbose) print*,"Orbit seemed to be unstable in longitudinal"
   
         fix = fix0
         fix(6)= fix(6)+ clight/freqmin/2
         
         didPhaseJump = .true.  ! it is protection against infinite lopp
         
         goto 1111
        endif
       endif
    endif


    if(use_bmad_units_temp) then 
 
      call convert_ptc_to_bmad(fix,beta1,STATE%TIME)
    endif
   use_bmad_units=use_bmad_units_temp
    !    FIX(6)=FIX(6)+freq*turns0
    c_%APERTURE_FLAG=APERTURE
    fix0=fix
    if(present(total)) total=tot
  END SUBROUTINE FIND_ORBIT_LAYOUT_noda_object




  SUBROUTINE FIND_ORBIT_LAYOUT_noda_object_orig(FIX0,STATE,eps,TURNS,fibre1,node1,total) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),pointer :: RING
    real(dp) , intent(inOUT) :: FIX0(6)
    real(dp), optional :: total
    INTEGER , optional,intent(in) :: TURNS
    type(fibre), optional, pointer :: fibre1
    type(integration_node), optional, pointer :: node1
    real(dp)  eps,freq,tot
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat

    real(dp)  DIX(6),xdix,xdix0,tiny,beta1,t6,freqmin
    real(dp) X(6),Y(6),MX(6,6),sxi(6,6),SX(6,6),dt,dl,fix(6)
    integer NO1,ND2,I,IU,ITE,ier,j,ITEM !,try
    TYPE (fibre), POINTER :: C
    TYPE (integration_node), POINTER :: t
    logical(lp) APERTURE,use_bmad_units_temp
    logical isStableFixPoint, didPhaseJump
    INTEGER TURNS0,trackflag


    fix=fix0

    didPhaseJump = .false.

    tot=0
    if(present(fibre1)) then
     ring=>fibre1%parent_layout
    else
     ring=>node1%parent_fibre%parent_layout
    endif


    if (.not.STATE%NOCAVITY.and.check_longitudinal) then
        freqmin=1.d38
        t6=0
        c=>ring%start    
      do i=1,ring%n
        if(c%mag%kind==kind4.or.c%mag%kind==kind21) then
          if(abs(c%mag%freq)>0) then
            if( abs(c%mag%freq)<=freqmin) freqmin=c%mag%freq
          endif
        endif
         t6=t6+c%mag%p%ld
       c=>c%next
      enddo 
      if(state%time) t6=t6/c%beta0
     if(global_verbose) write(6,*) " Harmonic # = ", freqmin*t6/clight
     if(global_verbose) write(6,*) " freqmin , dt", freqmin,clight/freqmin/2
    endif

    !    fixed_found=my_true
    !!    type(probe) xs
    if(.not.associated(RING%t)) call MAKE_NODE_LAYOUT(ring)
    !!    xs%x=zero
    !!    xs%s%x=zero
    use_bmad_units_temp=use_bmad_units
    if(use_bmad_units_temp) then 
          if(present(fibre1)) then
           beta1=fibre1%mag%p%beta0
          else
           beta1=node1%parent_fibre%mag%p%beta0
          endif
      call convert_bmad_to_ptc(fix,beta1,STATE%TIME)
    endif
    
    use_bmad_units=.false.

    TURNS0=1
    trackflag=0
    IF(PRESENT(TURNS)) TURNS0=TURNS
 
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    !!    call move_to(ring,c,loc)
    !!    loct=c%t1%pos


    Nullify(C);

 
    dix(:)=0.0_dp
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking
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
          
          messagelost= " No Cavity in the Line "
          check_stable=.false.
          return
 
       ENDIF
    else   ! (.not.present(STATE)) t
       IF(STATE%NOCAVITY) THEN
          ND2=4
          STAT=STATE+only_4d0
          if(state%radiation) then
             check_stable=.false.

             messagelost= " Cavity needed when radiation present "
             return
          endif
       ELSE
          ND2=6
          STAT=STATE
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          check_stable=.false.
          messagelost= " State present; no cavity: FIND_ORBIT_LAYOUT will crash => exiting"
         return

       ENDIF
    endif
101 continue


    if((stat%totalpath==1).and.(.not.stat%nocavity)) then
       C=>RING%START
       freq=0.0_dp
       i=1
       xdix=0.0_dp
       do while(i<=RING%n)
          if(associated(c%magp%freq)) then
             IF(FREQ==0.0_dp) THEN
                freq=c%magp%freq
             ELSEIF(c%magp%freq/=0.0_dp.AND.c%magp%freq<FREQ) THEN
                freq=c%magp%freq
             ENDIF
          endif
          c=>c%next
          i=i+1
       enddo
       
       
       if(global_verbose) write(6,*) " Using frequency ", freq

       
       if(ring%harmonic_number==0) then
          t=>ring%t%start
          tot=0
          do i=1,ring%t%n
             tot= tot + t%ds_ac
             t=>t%next
          enddo

          if(state%time) tot=tot/c%beta0

          if(global_verbose) write(6,*) " Integrated tot of the machine ", tot

          if(freq_redefine) then
             ring%harmonic_number=tot*freq/twopi
          else
             ring%harmonic_number=tot*freq/CLIGHT
          endif

          if(global_verbose) write(6,*) " Corresponding harmonic number ", ring%harmonic_number

       else

         if(freq_redefine) then
           tot=RING%HARMONIC_NUMBER*twopi/FREQ
         else
           tot=RING%HARMONIC_NUMBER*CLIGHT/FREQ
         endif

         if(global_verbose) write(6,*) " Harmonic number already defined", ring%harmonic_number
         if(global_verbose) write(6,*) " Corresponding tot ", tot
         
      endif

   endif
    



    
1111 continue


    ITEM=0
3   continue
    ITEM=ITEM+1
    X=FIX

    DO I=1,TURNS0
       call TRACK_probe_X(x,stat,fibre1=fibre1,node1=node1)

       if(.not.check_stable) then
          messagelost(len_trim(messagelost)+1:255)=" -> Unstable tracking guessed orbit "
          c_%APERTURE_FLAG=APERTURE
 !                if(try>0) goto 1111
          return
       endif
 

    ENDDO

    mx=0.0_dp
    DO J=1,ND2
       Y=FIX
       Y(J)=FIX(J)+EPS
       DO I=1,TURNS0
          call TRACK_probe_X(Y,stat,fibre1=fibre1,node1=node1)

          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
 !                   if(try>0) goto 1111
             return
          endif

          
       ENDDO
 
       if(stat%totalpath==1) then
         y(6)=y(6)-TURNS0*tot
       endif
       do i=1,ND2
            MX(I,J)=Y(i)/2/eps+MX(I,J)   
       enddo
       Y=FIX
       Y(J)=FIX(J)-EPS
       DO I=1,TURNS0
          !          CALL TRACK(RING,Y,LOC,STAT)
          !!       xs%x=y
          call TRACK_probe_X(Y,stat,fibre1=fibre1,node1=node1)

          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
 !                   if(try>0) goto 1111
             return
          endif
 
       ENDDO
 
      if(stat%totalpath==1) then
         y(6)=y(6)-TURNS0*tot
       endif

       do i=1,ND2
            MX(I,J)=-Y(i)/2/eps+MX(I,J)   
       enddo

       
    ENDDO

    SX=MX;
    DO I=1,nd2   !  6 before
       SX(I,I)=MX(I,I)-1.0_dp
    ENDDO

    DO I=1,ND2
       DIX(I)=FIX(I)-X(I)
      if(i==6) dix(6)=dix(6)+tot
    enddo
    
    CALL matinv(SX,SXI,ND2,6,ier)
    IF(IER==132)  then
       messagelost= " Inversion failed in FIND_ORBIT_LAYOUT_noda"
        check_stable=.false.
       return
    endif

    x=0.0_dp
    do i=1,nd2
       do j=1,nd2
          x(i)=sxi(i,j)*dix(j)+x(i)
       enddo
    enddo
    dix=x
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=0.0_dp
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo
    !    write(6,*) " Convergence Factor = ",nd2,xdix,deps_tracking
    !    pause 123321
  !  if(verbose) write(6,*) " Convergence Factor = ",xdix
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
       !   C_%stable_da=.FALSE.
       !      IF(iteM==MAX_FIND_ITER+100) THEN
       !        write(6,*) " Unstable in find_orbit without TPSA"
       messagelost= "Maximum number of iterations in find_orbit without TPSA"
       xlost=fix
       check_stable=my_false
       !     ENDIF
       ITE=0
   !           if(try>0) goto 1111
      return
    endif
    !   write(6,*) item,xdix,xdix0,tiny

    if(ite.eq.1.or.item<item_min)  then

       GOTO 3

    endif

    if (ND2 == 6.and.check_longitudinal) then
      isStableFixPoint = is_ORBIT_STABLE(FIX,EPS,STAT,fibre1,node1)

      if (isStableFixPoint .eqv. .false.) then
        if (didPhaseJump ) then

          messagelost= "Found unstable fixed point"
          xlost=fix
          check_stable=my_false

        else
        
         if(global_verbose) print*,"Orbit seemed to be unstable in longitudinal"
   
         fix = fix0
         fix(6)= fix(6)+ clight/freqmin/2
         
         didPhaseJump = .true.  ! it is protection against infinite lopp
         
         goto 1111
        endif
       endif
    endif


    if(use_bmad_units_temp) then 
 
      call convert_ptc_to_bmad(fix,beta1,STATE%TIME)
    endif
   use_bmad_units=use_bmad_units_temp
    !    FIX(6)=FIX(6)+freq*turns0
    c_%APERTURE_FLAG=APERTURE
    fix0=fix
    if(present(total)) total=tot
  END SUBROUTINE FIND_ORBIT_LAYOUT_noda_object_orig


 
    !!!!!!!!!!!!   about radiation and tapering  !!!!!!!!!!!!!!!!!

   SUBROUTINE FIND_ORBIT_tapering(FIX,eps,stat,f1) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    real(dp) , intent(inOUT) :: FIX(6)
    real(dp)  eps
    TYPE(INTERNAL_STATE) stat
    real(dp)  DIX(6),xdix,xdix0,tiny 
    real(dp) X(6),Y(6),MX(6,6),sxi(6,6),SX(6,6)
    integer NO1,ND2,I,IU,ITE,ier,j,ITEM,k
    TYPE (fibre), POINTER :: p,f1
    TYPE (integration_node), POINTER :: t

    logical(lp) APERTURE
    INTEGER trackflag
    type(work) w,we
    

    if(.not.associated(f1%parent_layout%t)) call MAKE_NODE_LAYOUT(f1%parent_layout)

 
          ND2=6
 
 
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.

    dix(:)=0.0_dp
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking
    
we=1


    ITEM=0
3   continue
    ITEM=ITEM+1
 !   goto 111
    X=FIX
    
       p=>f1
      do k=1,f1%parent_layout%n
          if(p%mag%kind/=kind4.and.p%mag%kind/=kind0.and.p%mag%kind/=kind1) then
              w=p
              we=1
             call find_energy(we,energy=p%mag%ene)
            p=we

            p%mag%ene= w%energy + x(5)*w%p0c
            we=-1

            call find_energy(we,energy=p%mag%ene)
            p=we
 
        endif
       call TRACK_probe_X(x,stat,fibre1=p,fibre2=p%next)
        if(.not.check_stable) then
           radfac=1
         write(6,*) "tapering failed in FIND_ORBIT_tapering 1",item
           return  
         endif

       p=>p%next
      enddo

    mx=0.0_dp
    DO J=1,ND2
       Y=FIX
       Y(J)=FIX(J)+EPS


       p=>f1
      do k=1,f1%parent_layout%n
          if(p%mag%kind/=kind4.and.p%mag%kind/=kind0.and.p%mag%kind/=kind1) then
              w=p
              we=1
            call find_energy(we,energy=p%mag%ene)

            p=we
            p%mag%ene= w%energy + y(5)*w%p0c
            we=-1
            call find_energy(we,energy=p%mag%ene)
            p=we
        endif
       call TRACK_probe_X(Y,stat,fibre1=p,fibre2=p%next)
        if(.not.check_stable) then
           radfac=1
         write(6,*) "tapering failed in FIND_ORBIT_tapering 2",item,j
           return  
         endif
       p=>p%next
      enddo

          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
             return
          endif

     
       do i=1,ND2

           MX(I,J)=Y(i)/2/eps+MX(I,J)
       enddo

    ENDDO

        DO J=1,ND2
       Y=FIX
       Y(J)=FIX(J)-EPS

       p=>f1
      do k=1,f1%parent_layout%n
          if(p%mag%kind/=kind4.and.p%mag%kind/=kind0.and.p%mag%kind/=kind1) then
              w=p
              we=1

           call find_energy(we,energy=p%mag%ene)
            p=we

            p%mag%ene= w%energy + y(5)*w%p0c
            we=-1

            call find_energy(we,energy=p%mag%ene)
            p=we
        endif
       call TRACK_probe_X(Y,stat,fibre1=p,fibre2=p%next)
        if(.not.check_stable) then
           radfac=1
         write(6,*) "tapering failed in FIND_ORBIT_tapering 3",item, j
           return  
         endif       
       p=>p%next
      enddo

          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
             return
          endif

     
       do i=1,ND2

           MX(I,J)=-Y(i)/2/eps+MX(I,J)
       enddo

    ENDDO
    
    SX=MX;
    DO I=1,nd2   !  6 before
       SX(I,I)=MX(I,I)-1.0_dp
    ENDDO

    DO I=1,ND2
       DIX(I)=FIX(I)-X(I)
    enddo
    
    CALL matinv(SX,SXI,ND2,6,ier)
    IF(IER==132)  then
       messagelost= " Inversion failed in FIND_ORBIT_LAYOUT_noda"
        check_stable=.false.
       return
    endif

    x=0.0_dp
    do i=1,nd2
       do j=1,nd2
          x(i)=sxi(i,j)*dix(j)+x(i)
       enddo
    enddo
    dix=x
    DO  I=1,ND2
       FIX(I)=FIX(I)+DIX(I)
    ENDDO

    xdix=0.0_dp
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo
    
    !    write(6,*) " Convergence Factor = ",nd2,xdix,deps_tracking
    !    pause 123321
  !  if(verbose) write(6,*) " Convergence Factor = ",xdix
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
       !   C_%stable_da=.FALSE.
       !      IF(iteM==MAX_FIND_ITER+100) THEN
       !        write(6,*) " Unstable in find_orbit without TPSA"
       write(6,*) "Maximum number of iterations in find_orbit without TPSA"
        check_Stable=.false.
        radfac=1
       return
       
    endif
    !   write(6,*) item,xdix,xdix0,tiny

    if(ite.eq.1.or.item<item_min)  then

       GOTO 3

    endif
    
    X=FIX
    
       p=>f1
      do k=1,f1%parent_layout%n
          if(p%mag%kind/=kind4.and.p%mag%kind/=kind0.and.p%mag%kind/=kind1) then
              w=p
              we=1
              
            call find_energy(we,energy=p%mag%ene)
            p=we
            p%mag%ene= w%energy + x(5)*w%p0c            
            we=-1

            call find_energy(we,energy=p%mag%ene)            
            p=we
        endif
       call TRACK_probe_X(x,stat,fibre1=p,fibre2=p%next)
        if(.not.check_stable) then
           radfac=1
         write(6,*) "tapering failed in FIND_ORBIT_tapering 4" 
           return  
         endif      
       p=>p%next
      enddo
fix=x
    c_%APERTURE_FLAG=APERTURE

   END SUBROUTINE FIND_ORBIT_tapering

subroutine taper(f1,fix,nsf,state,eps,file)  

implicit none
!!!!!!  PTC stuff
real(dp) x(6),fix(6),se,eps 
type(layout), pointer :: ring
type(internal_state) state
integer i,mf,k,nsf
type(fibre), pointer ::p,f1
type(work) w,we
character(*), optional :: file
type(c_damap) id
type(probe_8) rayp
type(probe) ray
type(c_normal_form) nf

ring=>f1%parent_layout

!!! fill in design energy

     w=0
        p=>ring%start
      do k=1,ring%n
        w=p
        p%mag%ene=w%energy
       p=>p%next
      enddo 
      
    
x=fix

do k=0,nsf
radfac=k
radfac=radfac/nsf

write(6,*) "iteration ",k

call FIND_ORBIT_tapering(x,eps,STATE,f1) 

call init(state,1,0)

call alloc(id)
call alloc(nf)
call alloc(rayp)
id=1
ray=x
rayp=id+ray

call propagate(rayp,state,fibre1=f1)
id=rayp
call c_normal(id,nf)

write(6,*) " tunes in tapering "
write(6,format3) nf%tune(1:3)

ray=rayp
call print(id%v(5))
call print(id%v(6))
write(6,format6) x
write(6,format6) ray%x
call kill(id)
call kill(nf)
call kill(rayp)

if(.not.check_stable) then
radfac=1
write(6,*) "tapering failed "
 return  
endif
enddo

write(6,*) "tapering apparently succesful "

if(present(file).and.file(1:6)/="nofile") then
write(6,*) " printing tapered closed orbit on file ", file(1:len_trim(file))
call kanalnummer(mf,file)


write(mf,'(7(1x,g12.5))') se,x(1:6)
write(6,'(6(1x,g12.5))') x(1:6)
p=>f1
do i=1,ring%n



call track_probe_x(x,state,fibre1=p,fibre2=p%next)
se=p%next%t1%s(1)

write(mf,'(7(1x,g12.5))') se,x(1:6)
p=>p%next
enddo
write(6,'(6(1x,g12.5))') x(1:6)
  close(mf)
endif

fix=x


end subroutine taper

subroutine untaper(f1)  

implicit none
!!!!!!  PTC stuff

type(layout), pointer :: ring
integer k 
type(fibre), pointer ::p,f1
type(work) w,we

ring=>f1%parent_layout

!!! fill in design energy

     w=0
        p=>ring%start
      do k=1,ring%n
          if(p%mag%kind/=kind4.and.p%mag%kind/=kind0.and.p%mag%kind/=kind1) then
              w=p
              we=1
              
            call find_energy(we,energy=p%mag%ene)
            p=we
            p%mag%ene= w%energy      
            we=-1
            call find_energy(we,energy=p%mag%ene)            
            p=we
        endif
       p=>p%next
      enddo 
      


end subroutine untaper

!!!! From Piotr Skowronski at CERN

  FUNCTION gettot(t6,freq)
    IMPLICIT NONE
    real(dp) gettot
    real(dp) t6,freq
    real(dp) harmon
    integer  harmoni
    
    if(freq_redefine) then
       harmon=t6*freq/twopi
    else
       harmon=t6*freq/CLIGHT
    endif
    
    harmoni = NINT(harmon)
    
    
    if(freq_redefine) then
      gettot=harmoni*twopi/FREQ
    else
      gettot=harmoni*CLIGHT/FREQ
    endif
    
    if(global_verbose) write(6,*) "cT=",t6, " Rounding Harmonic number ", harmon, " to ", harmoni, " corrsponding tot ", gettot
    
    
  END FUNCTION gettot
    

  FUNCTION is_ORBIT_STABLE(FIX,EPS,STAT,fibre1,node1)
    IMPLICIT NONE
    logical is_ORBIT_STABLE
    real(dp) , intent(inOUT) :: FIX(6)
    real(dp)  eps
    TYPE(INTERNAL_STATE) stat
    type(fibre), optional, pointer :: fibre1
    type(integration_node), optional, pointer :: node1
    integer i, nfailed
    real(dp)  :: Yi(6), Yf(6), Ydiff(6)
    
    is_ORBIT_STABLE = .true.
     
    nfailed = 0
   !  X=FIX
   !  call TRACK_probe_X(x,stat)
  
     !! TEST 1  
     Yi = FIX
     Yi(6)=FIX(6)+EPS
     Yf = Yi
     
     do i=1,3
       call TRACK_probe_X(Yf,stat,fibre1=fibre1,node1=node1)
     enddo
     
     Ydiff = Yf - Yi
     if(global_verbose) then
      write(*,'(2(f24.18,1x))') FIX(5:6)
      write(*,'(2(f24.18,1x))') Yi(5:6)
      write(*,'(2(f24.18,1x))') Yf(5:6)
      write(*,'(2(f24.18,1x))') Ydiff(5:6)
     endif
     
     
     if (Yf(6) .gt. Yi(6) ) then
     if(global_verbose) print*,"Test 1 failed"
       nfailed = nfailed + 1
     else
     if(global_verbose) print*,"Test 1 OK"
     endif
      

     !! TEST 2
     Yi = FIX 
     Yi(6)=FIX(6)-EPS
     Yf = Yi

     do i=1,3
       call TRACK_probe_X(Yf,stat,fibre1=fibre1,node1=node1)
     enddo
     
     Ydiff = Yf - Yi
     if(global_verbose) then
      write(*,'(2(f24.18,1x))') FIX(5:6)
      write(*,'(2(f24.18,1x))') Yi(5:6)
      write(*,'(2(f24.18,1x))') Yf(5:6)
      write(*,'(2(f24.18,1x))') Ydiff(5:6)
    endif
     if (Yf(6) .lt. Yi(6) ) then
     if(global_verbose)print*,"Test 2 failed"
       nfailed = nfailed + 1
     else
     if(global_verbose) print*,"Test 2 OK"
     endif

     !! TEST 3 
     Yi = FIX
     Yi(5)=FIX(5)+EPS
     Yf = Yi

     do i=1,3
       call TRACK_probe_X(Yf,stat,fibre1=fibre1,node1=node1)
     enddo
     
     Ydiff = Yf - Yi
     if(global_verbose) then
      write(*,'(2(f24.18,1x))') FIX(5:6)
      write(*,'(2(f24.18,1x))') Yi(5:6)
      write(*,'(2(f24.18,1x))') Yf(5:6)
      write(*,'(2(f24.18,1x))') Ydiff(5:6)
     endif
     if (Yf(5) .gt. Yi(5) ) then
     if(global_verbose) print*,"Test 3 failed"
       nfailed = nfailed + 1
     else
     if(global_verbose) print*,"Test 3 OK"
     endif

     !! TEST 4
     Yi = FIX 
     Yi(5)=FIX(5)-EPS
     Yf = Yi

     do i=1,3
       call TRACK_probe_X(Yf,stat,fibre1=fibre1,node1=node1)
     enddo
     
     Ydiff = Yf - Yi
     if(global_verbose) then
      write(*,'(2(f24.18,1x))') FIX(5:6)
      write(*,'(2(f24.18,1x))') Yi(5:6)
      write(*,'(2(f24.18,1x))') Yf(5:6)
      write(*,'(2(f24.18,1x))') Ydiff(5:6)
     endif
     if (Yf(5) .lt. Yi(5) ) then
     if(global_verbose) print*,"Test 4 failed"
       nfailed = nfailed + 1
     else
     if(global_verbose) print*,"Test 4 OK"
     endif
      
     if (nfailed .gt. 3) then
        is_ORBIT_STABLE = .FALSE.
     endif
  
  END FUNCTION is_ORBIT_STABLE

subroutine SMALL_CODE_TWISS(ring,no,general)
implicit none
integer no,np,mf,MFM,i,k
real(dp) fix(6),s(3,3,0:6)
type(internal_state), target  :: state,state_graph
type(layout), target  :: ring
type(c_damap) id,m,a_cs,disp,A_L,A_NL,a_spin,eval
type(probe) r0
type(probe_8) r
type(c_normal_form) normal
type(c_taylor) phase(3),phase_spin,betaxx,betaxy,x,c_spin_tune,d_CS,h,hb,pb_field
type(fibre),pointer :: f,sf,ff
integer, allocatable :: j1(:)
real(dp) phaser(3),spin_tune(2),damping(3),betx(2),dnu_dko,phx,a(6,6),mux,betasfx
type(c_linear_map)  q_cs,q_as,q_rot
logical general
type(c_spinor) e_y,isf
type(c_vector_field) field
complex(dp)  f4,f4green,f2,f2green
!   MAD-X LATTICE
!   L : drift, L= 0.2;
!   alpha= .3141592653589793238462643383279502;
! QF : SBEND,L= 1.0, ANGLE=ALPHA,k1=1.0;   
! QD : SBEND,L= 1.0, ANGLE=ALPHA,k1=-1.0;   
!   Oct : octupole,  K3= 0.0;
!lattice : LINE=  (QF,Oct,L,QD,Oct,L);

use_quaternion=.true.
!write(6,*) " General algorithm "
!read(5,*) general

!write(6,*) " give e1: 0 (normal quadrupole) or WM.INC secret value (0.25 for example)"
!read(5,*) e1_cas



MF=16
MFM=17
if(general) then
 call kanalnummer(mf,"TWISS.TXT")
 call kanalnummer(mf,"MAPS.TXT")
else
 call kanalnummer(mf,"TWISS_FAST.TXT")
 call kanalnummer(mf,"MAPS_FAST.TXT")
endif
 WRITE(MF,'(A4,6X,A7,7X,3X,A15,3X,4X,A14,9X,A10,9X,A11)')  &
"NAME", "PHASE_X", "DPHASE_X/DDELTA","DPHASE_X/DR2_X","   BETA   ","   ALPHA   "
 WRITE(MF,*) " "
 
phaser=0

 


np=0

 state=nocavity0+spin0
 

call init_all(state,no,np)

allocate(j1(c_%nd2t))

 eval%n=c_%nv
call alloc(r)
call alloc(id,m,a_cs,disp,A_L,A_NL,a_spin,eval)
call alloc(normal)
call alloc(phase_spin,betaxx,betaxy,x,d_CS,h,hb,pb_field)
call alloc(phase); call alloc(field)  
call alloc(e_y);call alloc(isf);

FIX=0.0_DP  ! FIXED POINT
fix(5)=0.00d0
call find_orbit(ring,fix(1:6),1,state,1.e-5_dp)   
 write(6,*)
 write(6,'(a12,5(1x,g12.5))')"Closed Orbit",  FIX(1:5)
 write(6,*)
! INITIALIZE THE RAY AS  --> 
!RAY = FIXED POINT + IDENTITY (TAYLOR MAP)
r0=fix

ID=1
R=r0+ID 

!  COMPUTING A ONE-TURN MAP TO ORDER MY_ORDER

call propagate(ring,r,state,fibre1=1)


M=R    

 CALL c_normal(M,NORMAL,phase=phase,nu_spin=phase_spin,dospin=state%spin)   
 write(6,'(a12,4(1x,g12.5))') "Linear tune ", NORMAL%tune(1:c_%nd),normal%spin_tune
 mux=NORMAL%tune(1)*twopi



 write(mf,*) "TOTAL TUNES "
 call print(phase(1:c_%nd),mf)


 if(state%spin) then
    call clean(phase_spin,phase_spin,prec=1.e-10_dp)
    write(mf,*) " Spin tune "
   call print(phase_spin,mf)
 endif
 

if(general) then
 call c_full_canonise(NORMAL%Atot,a_cs,a_spin,disp,A_L,A_NL)
else
 call c_fast_canonise(NORMAL%Atot,a_cs,dospin=state%spin)
endif
 

phase=0.0_dp
phaser=0
spin_tune=0
damping=0
phase_spin=0.0_dp
dnu_dko=0
d_CS=0.0_dp
h =(1.0_dp.cmono.1)+i_*(1.0_dp.cmono.2)
hb=(1.0_dp.cmono.1)-i_*(1.0_dp.cmono.2)
f=>ring%start
r=a_cs+r0
eval=1
eval%v(3)=0.0e0_dp
eval%v(4)=0.0e0_dp
eval%v(5)=0.0e0_dp

d_CS=h*hb
d_cs=d_cs*NORMAL%Atot**(-1)
call print(d_cs)

 do i=1,ring%n

call propagate(ring,r,+state,fibre1=i,fibre2=i+1)

 write(6,*) I,f%mag%name

a_cs=r
if(general) then
 call c_full_canonise(a_cs,a_cs,a_spin, &
  disp,A_L,A_NL,phase=phase,nu_spin=phase_spin)
else
 call c_fast_canonise(a_cs,a_cs,phaser,damping,q_cs=q_cs,q_as=q_as, &
 q_rot=q_rot,spin_tune=spin_tune ,dospin=state%spin)
endif

r0=r
r=a_cs+r0


!!!!!!!!!!!!! Do something !!!!!!!!!!!!!!!!!!!!!
write(mf,*) f%mag%name
if(general) then
 j1=0
 j1(1)=1
 betaxx=(A_L%v(1).par.j1)**2
 j1=0
 j1(2)=1
 betaxx=betaxx + (A_L%v(1).par.j1)**2
 write(mf,*) " Betax_1 "
 call print(betaxx,mf)
 j1=0
 j1(3)=1
 betaxy=(A_L%v(1).par.j1)**2
 j1=0
 j1(4)=1
 betaxy=betaxy + (A_L%v(1).par.j1)**2
 write(mf,*) " Betax_2 "
 call print(betaxy,mf)
write(mf,*) " ISF "
 e_y=2
 call makeso3(a_spin)
 isf=a_spin%s*e_y
 call clean(isf,isf,prec=1.0e-10_dp)
 call print(isf,mf)

 betx(1)=betaxx
 write(mf,*) " phases  "
 call clean(phase(1:c_%nd),phase(1:c_%nd),prec=1.0e-10_dp)
 call print(phase(1:c_%nd),mf)
 if(state%spin) then
  write(mf,*) " spin tune "
  call clean(phase_spin,phase_spin,prec=1.0e-10_dp)
  call print(phase_spin,mf)
 endif
else
 betx(1)=q_cs%mat(1,1)**2+q_cs%mat(1,2)**2
 betx(2)=q_cs%mat(1,3)**2+q_cs%mat(1,3)**2
 write(mf,*) " Betax_1 , Betax_2"
 write(mf,*) betx
 write(mf,*) " Phases "
 write(mf,*) phaser

 if(state%spin) then
  write(mf,*) " spin tune "
  write(mf,*) spin_tune
  write(mf,*) " ISF "
  call MAKESO3(q_as,s)

  write(mf,'(10x,a4,15x,a7,16x,a9,12x,a7,15x,a9,11x,a11)') 'n0_x','dn_x/dx','dn_x/dp_x','dn_x/dy','dn_x/dp_y','dn_x/ddelta' 
  write(mf,'(6(1x,G21.14))') s(1,2,0:5)
  write(mf,'(10x,a4,15x,a7,16x,a9,12x,a7,15x,a9,11x,a11)') 'n0_y','dn_y/dx','dn_y/dp_x','dn_y/dy','dn_y/dp_y','dn_y/ddelta' 
  write(mf,'(6(1x,G21.14))') s(2,2,0:5)
  write(mf,'(10x,a4,15x,a7,16x,a9,12x,a7,15x,a9,11x,a11)') 'n0_z','dn_z/dx','dn_z/dp_x','dn_z/dy','dn_z/dp_y','dn_z/ddelta' 
  write(mf,'(6(1x,G21.14))') s(3,2,0:5)
 endif


 endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 

!!!!!!!!!!!!  Analytical  !!!!!!!!!!!!!!!!!!!!



Write(mf,*) " change in invariant due to octupole : TPSA "
d_CS=(h*hb)*a_cs**(-1)
d_cs=d_cs-(d_cs.cut.5)
d_cs=(d_cs.o.eval).d.c_%nv
call clean(d_cs,d_cs,prec=1.0e-10_dp)
call print(d_cs,mf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



f=>f%next
enddo




mux=phase(1)
mux=mux*twopi
Write(mf,*) " phase, total phase and beta "
 write(mf,*) phx,mux,betx(1)

!phx=phx*twopi
!d_CS=8*i_*exp(-i_*4*phx)/(1.0_dp-exp(-i_*4.0_dp*mux))*h**4+16*i_*exp(-i_*2*phx)/(1.0_dp-exp(-i_*2.0_dp*mux))*h**3*hb
!d_cs=d_cs*(a_cs.cut.2)**(-1)
!d_cs=2.0_dp*real(d_cs)*betasfx**2/64.0_dp
!d_cs=d_cs.o.eval
!call clean(d_cs,d_cs,prec=1.d-10)
!Write(mf,*) " change in invariant due to octupole : analytical "
!call print(d_cs,mf)


write(6,*) "dnu_dko = ",dnu_dko
 write(mf,*) " <x^2> "
   x=2.0_dp*(1.0e0_dp.cmono.1)**2
   call average(x,a_l,x)
   call print(x,mf)
 
 WRITE(MF,'(A4,6X,A7,7X,3X,A15,3X,4X,A14,9X,A10,9X,A11)')  &
"NAME", "PHASE_X", "DPHASE_X/DDELTA","DPHASE_X/DR2_X","   BETA   ","   ALPHA   "

!!!!!! compute guignard  !!!!!!



 
223 CLOSE(MF)
CLOSE(MFM) 

END subroutine SMALL_CODE_TWISS



end module S_fitting_new
