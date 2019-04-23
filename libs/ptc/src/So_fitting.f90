!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_fitting 
  USE MAD_LIKE
  IMPLICIT NONE
  public
  PRIVATE FIND_ORBIT_LAYOUT, FIND_ORBIT_LAYOUT_noda
  logical(lp), PRIVATE :: VERBOSE = .false.
  integer :: max_fit_iter=20, ierror_fit=0, max_fiND_iter=40
  real(dp) :: fuzzy_split=1.0_dp
  real(dp) :: max_ds=0.0_dp
  integer :: resplit_cutting = 0    ! 0 just magnets , 1 magnets as before / drifts separately

  logical :: sagan_even=my_true
  ! 2  space charge algorithm
  logical(lp) :: radiation_bend_split=my_false
  type(mad_universe),private, pointer :: m_u=>null()
  type(mad_universe),private, pointer :: m_t=>null()

  INTERFACE FIND_ORBIT
     ! LINKED
     ! no use of TPSA
     MODULE PROCEDURE FIND_ORBIT_LAYOUT_noda
  END INTERFACE


contains
  SUBROUTINE lattice_GET_CHROM(R,my_state,CHROM)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    REAL(DP) CHROM(:)
    TYPE(internal_state) state
    real(dp) closed(6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)

    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

    closed=0.0_dp
    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)
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
    if(.not.check_stable) then
     CALL RESET_APERTURE_FLAG
     write(6,*) " Flags were reset in lattice_GET_CHROM"
    endif

  end SUBROUTINE lattice_GET_CHROM


  SUBROUTINE lattice_GET_tune(R,my_state,mf,targ)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    TYPE(internal_state) state
    integer mf
    real(dp) closed(6),targ(3)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)
    
    

    STATE=my_state
  
    closed=0.0_dp
    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)
    write(6,*) "closed orbit "
    WRITE(6,'(6(1x,g21.14))') CLOSED
    CALL INIT(STATE,1,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y
    closed=y
    targ(1:3)=NORM%tune(1:3)
    WRITE(6,'(6(1x,g21.14),a24)') CLOSED," <-- should be identical"
    if(mf==6) then
     WRITE(6,'(a19,3(1x,g21.14))') "Fractional Tunes = ",norm%tune(1:3)
     if(norm%tune(3)/=0.0_dp.and.c_%ndpt==0) &
     WRITE(6,'(a20,(1x,g21.14))') "Synchrotron period = ",1.d0/abs(norm%tune(3))
     else
     if(norm%tune(3)/=0.0_dp.and.c_%ndpt==0) then
       WRITE(mf,'(4(1x,g21.14))') xsm0t/clight,norm%tune(1:3)
     else
       WRITE(mf,'(3(1x,g21.14))') xsm0t/clight,norm%tune(1:2)
     endif
    endif
    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)
    if(.not.check_stable) then
     CALL RESET_APERTURE_FLAG
     write(6,*) " Flags were reset lattice_GET_tune"
    endif
    
  end SUBROUTINE lattice_GET_tune



  SUBROUTINE compute_A_4d(r,my_state,filename,pos,del,no,MY_A)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
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

    closed=0.0_dp
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,STATE,1e-5_dp)
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
    if(.not.check_stable) then
     CALL RESET_APERTURE_FLAG
     write(6,*) " Flags were reset in compute_A_4d"
    endif

  end SUBROUTINE compute_A_4d



  SUBROUTINE compute_map_general(r,my_state,filename,pos,del,no)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    TYPE(internal_state) state
    integer pos,no,is,mf
    real(dp) closed(6),del
    type(DAMAP) ID
    TYPE(REAL_8) Y(6)
    CHARACTER(*) FILENAME
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

    expo1=0.0_dp;expo2=0.0_dp;
    expo3=0.0_dp;expo4=0.0_dp;

    closed=0.0_dp
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,my_state,1e-5_dp)
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
    id=y

    call print(id,mf)
    close(mf)
    call kill(id); call kill(norm);call kill(y);call kill(betax,betax2);
    deallocate(expo1,expo2)
    deallocate(expo3,expo4)

  END SUBROUTINE compute_map_general


  SUBROUTINE compute_map_4d(r,my_state,filename,pos,del,no)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    integer pos,no,imod,ic,I
    TYPE(internal_state) state
    real(dp) closed(6),del
    type(DAMAP) ID
    TYPE(REAL_8) Y(6)
    CHARACTER(*) FILENAME
    TYPE(FIBRE), POINTER :: P



    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    closed=0.0_dp
    closed(5)=del
    CALL FIND_ORBIT(R,CLOSED,pos,STATE,1e-5_dp)
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
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    REAL(DP), pointer :: BETA(:,:,:)
    REAL(DP)DBETA,tune(:),tune2(:)
    type(fibre),pointer :: p
    integer i,IB,pos,mf
    TYPE(internal_state) state
    real(dp) closed(6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)
    real(dp) dbetamax,db1,db2,db31,db32,eta1,eta2
    real(dp),optional :: a(6,6),ai(6,6),mat(6,6),clos(6)

    if(.not.associated(beta))   ALLOCATE(BETA(2,3,R%N))

eta1=0.0_dp
eta2=0.0_dp
    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)
    state=state+delta0

    if(present(clos)) then
       closed=clos
    else
       closed=0.0_dp
    endif
    if(pos/=0) then
       CALL FIND_ORBIT(R,CLOSED,pos,STATE,1e-5_dp)
       write(6,*) "closed orbit "
       write(6,*) CLOSED
    else
       write(6,*) " Using a map "
    endif
    DBETA=0.0_dp
    dbetamax=0.0_dp
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
       a=0.0_dp
       ai=0.0_dp
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
          beta(IB,3,i-pos+1)=y(1).sub.'00001'
          
          IF(IB==2) THEN
             db1=ABS(beta(2,1,i-pos+1)-beta(1,1,i-pos+1))/beta(1,1,i-pos+1)
             db2=ABS(beta(2,2,i-pos+1)-beta(1,2,i-pos+1))/beta(1,2,i-pos+1)
             db31=ABS(beta(1,3,i-pos+1))
             db32=ABS(beta(2,3,i-pos+1)-beta(1,3,i-pos+1))
             DBETA=(db1+db2)/2.0_dp+dbeta
             eta1=db31+eta1
             eta2=db32+eta2
             if( db1>dbetamax) dbetamax=db1
             if( db2>dbetamax) dbetamax=db2
          ENDIF
          p=>p%next
       enddo
       DBETA=DBETA/R%N

       IF(IB==2) WRITE(6,*) "<DBETA/BETA> = ",DBETA
       IF(IB==2) WRITE(6,*) "MAXIMUM OF DBETA/BETA = ",dbetamax
       IF(IB==2) WRITE(6,*) "<DETA/ETA> = ",eta2/eta1
 
    endif
    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)
    if(present(clos)) clos=closed

  end SUBROUTINE FILL_BETA

  SUBROUTINE comp_longitudinal_accel(r,my_state,no,h,filename)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state) state
    real(dp) closed(6),l1,l2,p0c
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)
    integer no,mf,i,h
    CHARACTER(*) FILENAME
    integer, allocatable :: j(:)
    TYPE(internal_state), intent(in):: my_state
    TYPE(fibre), pointer :: p
    logical first
    type(work) w
    p=> r%start

    l1=0.d0
    l2=0.d0
    first=.true.
    w=p
    p0c=w%p0c

    call kanalnummer(mf,filename)

    write(6,*)h,w%mass,w%p0c,w%beta0
    write(mf,*) h
    write(mf,*) w%mass,w%p0c,w%beta0

    do i=1,r%n
       if(p%mag%kind/=kind4) then
          if(first) then
             l1=l1+p%mag%p%ld
          else
             l2=l2+p%mag%p%ld
          endif
       else
          if(.not.first) stop 111
          l1=l1+p%mag%p%ld/2
          l2=p%mag%p%ld/2
          first=.false.
       endif
       p=>p%next
    enddo
    write(6,*) l1,l2
    write(mf,*) l1,l2


    STATE=my_state+nocavity0

    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)

    CALL INIT(STATE,no,0)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y

    allocate(j(c_%nv))
    j=0
    do i=1,no
       j(5)=i
       write(6,*) norm%dhdj%v(3).sub.j
       write(mf,*) norm%dhdj%v(3).sub.j
    enddo

    deallocate(j)

    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)
    close(mf)

  end SUBROUTINE comp_longitudinal_accel


  SUBROUTINE comp_linear2(r,my_state,a,ai,mat,closed)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: my_state
    TYPE(internal_state) state
    real(dp) closed(6),a(6,6),ai(6,6),mat(6,6)
    type(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)


    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)

    CALL INIT(STATE,1,0,BERZ)

    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,STATE)
    NORM=Y

    Write(6,*) " Tunes ",norm%tune(1:2)
    A=0.0_dp
    AI=0.0_dp
    MAT=0.0_dp
    a=norm%a_t
    AI=norm%a_t**(-1)
    id=y
    mat=id
    CALL kill(NORM)
    CALL kill(Y)
    call kill(id)

  end SUBROUTINE comp_linear2


  subroutine lattice_fit_TUNE_gmap_auto(R,my_state,EPSF,TARG,name)
    IMPLICIT NONE
    TYPE(layout), target,intent(inout):: R
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6)
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE, MF
    TYPE(TAYLOR), allocatable:: EQ(:)
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer :: neq=2, no=2,nt,j,it,NP
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    character(nlp) name(2)
    real(dp) epsf,epsr,epsnow,gam(2)
    type(fibre), pointer:: p
    logical dname(2)
    TYPE(POL_BLOCK) poly(2)
    call context(name(1))
    call context(name(2))
    poly(1)=0
    poly(2)=0
    poly(1)%ibn(2)=1
    poly(2)%ibn(2)=2
    !    EPSF=.0001
    SET_TPSAFIT=.FALSE.
    dname=.false.
    p=>r%start
    do i=1,r%n
       if(index (p%mag%name,name(1)(1:len_trim(name(1))) )   /=0) then
          call  EL_POL_force(p,poly(1))
          dname(1)=.true.
       elseif(index(p%mag%name,name(2)(1:len_trim(name(2))))/=0) then
          call  EL_POL_force(p,poly(2))
          dname(2)=.true.
       endif

       p=>p%next
    enddo

    if(.not.(dname(1).and.dname(2))) then
       CALL ELP_TO_EL(R)
       write(6,*) " lattice_fit_TUNE_gmap_auto ---> FAILED"
       return
    endif


    epsr=abs(epsf)

    np=2

    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

    CALL INIT(STATE,no,NP)


    !   DO I=1,NPOLY
    !      R=POLY(i)
    !   ENDDO


    CLOSED(:)=0.0_dp
    it=0
100 continue
    it=it+1

    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)
    write(6,*) "closed orbit "
    write(6,*) CLOSED


    CALL INIT(STATE,no,NP)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,+STATE)
    write(6,*) "c_%no,c_%nv,c_%nd,c_%nd2"
    write(6,*) c_%no,c_%nv,c_%nd,c_%nd2
    write(6,*) "c_%ndpt,c_%npara,c_%npara,c_%np_pol"
    write(6,*)  c_%ndpt,c_%npara,c_%npara,c_%np_pol
    NORM=Y
    gam(1)=(norm%a_t%v(2).sub.'1')**2+(norm%a_t%v(2).sub.'01')**2
    gam(2)=(norm%a_t%v(4).sub.'001')**2+(norm%a_t%v(4).sub.'0001')**2
    write(6,*) "  Gamma= ",GAM
    !      CALL KANALNUMMER(MF)
  !  OPEN(UNIT=1111,FILE='GAMMA.TXT')
  !  WRITE(1111,*) "  Gamma= ",GAM

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

    p=>r%start
    do i=1,r%n
       if(index (p%mag%name,name(1)(1:len_trim(name(1))) )   /=0) then
          call  EL_POL_force(p,poly(1))
       elseif(index(p%mag%name,name(2)(1:len_trim(name(2))))/=0) then
          call  EL_POL_force(p,poly(2))
       endif

       p=>p%next
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
    deallocate(eq)

  end subroutine lattice_fit_TUNE_gmap_auto

  subroutine lattice_fit_TUNE_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
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
    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)

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
    write(6,*) "c_%no,c_%nv,c_%nd,c_%nd2"
    write(6,*) c_%no,c_%nv,c_%nd,c_%nd2
    write(6,*) "c_%ndpt,c_%npara,c_%npara,c_%np_pol"
    write(6,*)  c_%ndpt,c_%npara,c_%npara,c_%np_pol


    NORM=Y
    gam(1)=(norm%a_t%v(2).sub.'1')**2+(norm%a_t%v(2).sub.'01')**2
    gam(2)=(norm%a_t%v(4).sub.'001')**2+(norm%a_t%v(4).sub.'0001')**2
    write(6,*) "  Gamma= ",GAM
    !      CALL KANALNUMMER(MF)
   ! OPEN(UNIT=1111,FILE='GAMMA.TXT')
   ! WRITE(1111,*) "  Gamma= ",GAM

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

  end subroutine lattice_fit_TUNE_gmap

  subroutine lattice_fit_CHROM_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout),target, intent(inout):: R
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

  end subroutine lattice_fit_CHROM_gmap

  subroutine lattice_fit_CHROM_gmap2(R,my_state,EPSF,POLY,NPOLY,TARG,np,n_extra,mf)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(layout),target, intent(inout):: R
    TYPE(POL_BLOCK), intent(inout),dimension(:)::POLY
    INTEGER, intent(in):: NPOLY,np
    real(dp) , intent(IN),dimension(:)::TARG
    real(dp) CLOSED(6),co
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    INTEGER I,SCRATCHFILE
    TYPE(TAYLOR), allocatable:: EQ(:)
    TYPE(REAL_8) Y(6)
    TYPE(NORMALFORM) NORM
    integer ::  neq,no=3,nt,j,it,n_extra,mf
    type(damap) id
    type(vecresonance) vr
    type(pbresonance) fr
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,CHROM(2)
    integer, allocatable:: res(:,:),jt(:),je(:)
    real(dp), allocatable :: mat(:,:),v0(:)
    integer kb,kbmax,kpos,npi,ier


    neq=2
    allocate(res(n_extra,4))
    allocate(jt(5))
    res=0
    do i=1,n_extra
       read(mf,*) res(i,:)
    enddo
    !    EPSF=.0001
    epsr=abs(epsf)
    neq=neq+2*n_extra

    allocate(eq(neq))

    nt=neq+np
    allocate(mat(nt,nt),v0(nt))

    kbmax=0

    DO I=1,NPOLY
       if(POLY(i)%nb>kbmax)  kbmax= POLY(i)%nb
    ENDDO
    SET_TPSAFIT=.FALSE.

    DO I=1,NPOLY
       R=POLY(i)
    ENDDO


    it=0
100 continue
    it=it+1
    mat=0.0_dp
    v0=0.0_dp
    do i=1,np
       mat(i,i)=1.0_dp
    enddo

    do kb=1,kbmax
       STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

       DO I=1,NPOLY
          if(POLY(i)%nb==kb) then
             npi=POLY(i)%np
             kpos=POLY(i)%g-1
             exit
          endif
       ENDDO

       write(6,*) " np in batch ",kb," = ",npi
       !         pause 1

       CLOSED(:)=0.0_dp

       CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)
       write(6,*) "closed orbit ", CHECK_STABLE
       write(6,*) CLOSED


       CALL INIT(STATE,no,NPi,BERZ)
       nb_=kb
       CALL ALLOC(NORM)
       CALL ALLOC(Y)
       CALL ALLOC(EQ)
       call alloc(id)
       call alloc(vr)
       call alloc(fr)

       id=1
       Y=CLOSED+id

       CALL TRACK(R,Y,1,+STATE)
       NORM=Y
       vr=norm%a%nonlinear
       fr=norm%a%pb
       !    call print(vr%cos%v(2),6)
       !    call print(vr%sin%v(2),6)
       !    pause 1
       ! call print(fr%cos,6)
       ! call print(fr%sin,6)
       ! pause 2


       write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2), CHECK_STABLE
       CHROM(1)=(NORM%dhdj%v(1)).SUB.'00001'
       CHROM(2)=(NORM%dhdj%v(2)).SUB.'00001'
       write(6,*) " CHROM ",CHROM

       eq(1)=       ((NORM%dhdj%v(1)).par.'00001')-targ(1)
       eq(2)=       ((NORM%dhdj%v(2)).par.'00001')-targ(2)
       do i=1,n_extra
          jt=0
          jt(1:4)=res(i,:)
          jt(1)=jt(1)-1
          eq(2+2*i-1)=       ((vr%cos%v(2)).par.jt)
          eq(2+2*i)=       ((vr%sin%v(2)).par.jt)
       enddo

       epsnow=0.0_dp
       do i=1,neq
          epsnow=abs(eq(i))+epsnow
          if(kb==1) then
             co=eq(i)
             if(i<=2) then
                co=co+targ(i)
             endif
             write(6,*) i,co
          endif
       enddo
       write(6,*) "epsnow ", epsnow
       ipause=mypause(123)
       call kanalnummer(SCRATCHFILE)
       OPEN(UNIT=SCRATCHFILE,FILE='EQUATION.TXT')
       rewind scratchfile
       allocate(je(c_%nv))
       !    write(6,*) nt,neq,np,kpos,npi
       je=0
       do i=1,neq
          eq(i)=eq(i)<=c_%npara
          v0(i+np)=-(eq(i).sub.'0')
          do j=1,npi
             je(j)=1
             co=eq(i).sub.je
             mat(np+i,j+kpos)=co
             mat(j+kpos,np+i)=co
             je(j)=0
          enddo
       enddo
       deallocate(je)
       do i=1,neq
          call daprint(eq(i),scratchfile)
       enddo
       close(SCRATCHFILE)
       CALL KILL(NORM)
       CALL KILL(Y)
       CALL KILL(id)
       CALL KILL(vr)
       CALL KILL(fr)
       CALL KILL(EQ)
       !    pause 888
    enddo ! kbmax
    write(6,*) "Iteration # ",it

    !    pause 2000
    !    do i=1,nt
    !    do j=1,nt
    !    if(mat(i,j)/=zero) write(6,*) i,j,mat(i,j)
    !    enddo
    !    enddo
    call  matinv(mat,mat,nt,nt,ier)
    if(ier/=0 ) then
       write(6,*) ier
       write(6,*) " inversion error "

       stop
    endif
    !    pause 2001

    v0=matmul(mat,v0)
    tpsafit=0.0_dp
    tpsafit(1:nt)=v0

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
    deallocate(mat)
    deallocate(eq)
    deallocate(res)
    deallocate(jt)

  end subroutine lattice_fit_CHROM_gmap2

  subroutine lattice_fit_CHROM_gmap1(R,my_state,EPSF,POLY,NPOLY,TARG,NP,n_extra,mf)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(layout),target, intent(inout):: R
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
    integer :: neq=2, no=3,nt,j,it,n_extra,mf
    type(damap) id
    type(vecresonance) vr
    type(pbresonance) fr
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,CHROM(2)
    integer, allocatable:: res(:,:),jt(:)

    neq=2

    allocate(res(n_extra,4))
    allocate(jt(5))
    res=0
    do i=1,n_extra
       read(mf,*) res(i,:)
    enddo
    !    EPSF=.0001
    epsr=abs(epsf)
    neq=neq+2*n_extra
    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

    !    CALL INIT(STATE,no,NP,BERZ)

    !    SET_TPSAFIT=.FALSE.

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
    call alloc(vr)
    call alloc(fr)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,+STATE)
    NORM=Y
    vr=norm%a%nonlinear
    fr=norm%a%pb
 !       call print(vr%cos%v(2),6)
 !       call print(vr%sin%v(2),6)
 !       pause 1
    call print(fr%cos,6)
    call print(fr%sin,6)
!write(6,*) " pausing -> type 1 "
!read(5,*) i


    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2), CHECK_STABLE
    CHROM(1)=(NORM%dhdj%v(1)).SUB.'00001'
    CHROM(2)=(NORM%dhdj%v(2)).SUB.'00001'
    write(6,*) " CHROM ",CHROM

    eq(1)=       ((NORM%dhdj%v(1)).par.'00001')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'00001')-targ(2)
    do i=1,n_extra
       jt=0
       jt(1:4)=res(i,:)
       jt(1)=jt(1)-1
       eq(2+2*i-1)=       ((vr%cos%v(2)).par.jt)
       eq(2+2*i)=       ((vr%sin%v(2)).par.jt)
    enddo

    epsnow=0.0_dp
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
    CALL KILL(id)
    CALL KILL(vr)
    CALL KILL(fr)
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
    deallocate(res)
    deallocate(jt)

  end subroutine lattice_fit_CHROM_gmap1

 subroutine c_lattice_fit_CHROM_gmap1(R,my_state,EPSF,POLY,NPOLY,TARG,NP,n_extra,mf)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(layout),target, intent(inout):: R
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
    integer :: neq=2, no=3,nt,j,it,n_extra,mf
    type(damap) id
    type(vecresonance) vr
    type(pbresonance) fr
    type(gmap) g
    TYPE(TAYLOR)t
    type(c_normal_form) cn
    type(c_damap) cmap
    type(c_vector_field) cvec

    real(dp) epsf,epsr,epsnow,CHROM(2)
    integer, allocatable:: res(:,:),jt(:)

    neq=2

    allocate(res(n_extra,4))
    allocate(jt(5))
    res=0
    do i=1,n_extra
       read(mf,*) res(i,:)
    enddo
    !    EPSF=.0001
    epsr=abs(epsf)
    neq=neq+2*n_extra
    allocate(eq(neq))

    nt=neq+np
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)

    !    CALL INIT(STATE,no,NP,BERZ)

    !    SET_TPSAFIT=.FALSE.

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


    CALL INIT_all(STATE,no,NP,BERZ)
    CALL ALLOC(NORM)
    CALL ALLOC(Y)
    CALL ALLOC(EQ)
    call alloc(id)
    call alloc(vr)
    call alloc(fr)

    call alloc(cn);call alloc(cmap); 
    call alloc(cvec)

    id=1
    Y=CLOSED+id

    CALL TRACK(R,Y,1,+STATE)
    norm=y
    cmap=y
    call c_normal(cmap,cn,dospin=.false.)
    call flatten_c_factored_lie(cn%ker,cvec)
    cvec=transform_vector_field_by_map(cvec,to_phasor())
    cmap=(cn%a_t)**(-1)
    cvec=transform_vector_field_by_map(cvec,cmap) 
    cmap=(cn%a_t).sub.1
    cvec=transform_vector_field_by_map(cvec,cmap)
    cvec=transform_vector_field_by_map(cvec,from_phasor())
  vr%cos%v(2)=real(cvec%v(2))
  vr%sin%v(2)=aimag(cvec%v(2))
!    NORM=Y
!    vr=norm%a%nonlinear
!    fr=norm%a%pb
 !       call print(vr%cos%v(2),6)
 !       call print(vr%sin%v(2),6)
 !       pause 1
!call print(cvec%v(2),6)
!    call print(vr%cos%v(2),6)
!    call print(vr%sin%v(2),6)
!pause 10
!write(6,*) " pausing -> type 1 "
!read(5,*) i


    write(6,*) " tunes ",NORM%TUNE(1), NORM%TUNE(2), CHECK_STABLE
    CHROM(1)=(NORM%dhdj%v(1)).SUB.'00001'
    CHROM(2)=(NORM%dhdj%v(2)).SUB.'00001'
    write(6,*) " CHROM ",CHROM

    eq(1)=       ((NORM%dhdj%v(1)).par.'00001')-targ(1)
    eq(2)=       ((NORM%dhdj%v(2)).par.'00001')-targ(2)
    do i=1,n_extra
       jt=0
       jt(1:4)=res(i,:)
       jt(1)=jt(1)-1
       eq(2+2*i-1)=       ((vr%cos%v(2)).par.jt)
       eq(2+2*i)=       ((vr%sin%v(2)).par.jt)
    enddo

    epsnow=0.0_dp
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
    CALL KILL(id)
    CALL KILL(vr)
    CALL KILL(fr)
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
    deallocate(res)
    deallocate(jt)

  end subroutine c_lattice_fit_CHROM_gmap1


  subroutine lattice_fit_tune_CHROM_gmap(R,my_state,EPSF,POLY,NPOLY,TARG,NP)
    IMPLICIT NONE
    TYPE(layout),target, intent(inout):: R
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
    integer :: neq=4, no=3,nt,j,it
    type(damap) id
    type(gmap) g
    TYPE(TAYLOR)t
    real(dp) epsf,epsr,epsnow,tune(2),CHROM(2),gam(2)
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
    gam(1)=(norm%a_t%v(2).sub.'1')**2+(norm%a_t%v(2).sub.'01')**2
    gam(2)=(norm%a_t%v(4).sub.'001')**2+(norm%a_t%v(4).sub.'0001')**2
    write(6,*) "  Gamma= ",GAM
    !      CALL KANALNUMMER(MF)
  !  OPEN(UNIT=1111,FILE='GAMMA.TXT')
  !  WRITE(1111,*) "  Gamma= ",GAM

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

  end subroutine lattice_fit_tune_CHROM_gmap


  subroutine lattice_PRINT_RES_FROM_A(R,my_state,NO,EMIT0,MRES,FILENAME)
    IMPLICIT NONE
    TYPE(layout),target, intent(inout):: R
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

    PREC=1e-10_dp

    STATE=((((my_state+nocavity0)-delta0)+only_4d0)-RADIATION0)



    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)
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



  subroutine lattice_random_error_new(R,nom,full,iseed,cut,n,addi,integrated,cn,cns,per,pr)
    use gauss_dis
    IMPLICIT NONE
    TYPE(layout),target, intent(inout):: R
    integer n,addi,ic,i,iseed,j
    character(nlp) nom
    type(fibre), pointer :: p
    logical(lp) integrated,f1,f2,full,pr
    real(dp) x,bn,cn,cns,cut,per

    if(iseed/=0) call gaussian_seed(iseed)


    call context(nom)

    ic=0
    p=>r%start
    do i=1,r%n
       f1=.false.
       f2=.false.
       if(full) then
          f1=(p%mag%name ==nom )
       else
          f1=(   index(p%mag%name,nom(1:len_trim(nom)))   >0)
       endif


       if(f1.and.per<=0) then
          call GRNF(X,cut)
          bn=cns+cn*x
          if(integrated.and.p%mag%p%ld/=0.0_dp) then
             bn=bn/p%mag%l
          endif
          if(bn/=0.0_dp) call add(p,n,addi,bn)
          f2=.true.
       endif

       if(f1.and.per/=0.0_dp) then
          call GRNF(X,cut)
          do j=p%mag%p%nmul,1,-1
             !        if(n>0) then
             bn=p%mag%bn(j)
             bn=bn*(1.0_dp+x*per)
             call add(p,j,0,bn)
             !        else
             bn=p%mag%an(j)
             bn=bn*(1.0_dp+x*per)
             !        endif
             call add(p,-j,0,bn)

          enddo
          f2=.true.
       endif
       if(f2) then 
         ic=ic+1
        if(pr) write(6,*) p%mag%name
       endif
       p=>P%next
    enddo

   if(pr) write(6,*) ic," Magnets modified "

  end  subroutine lattice_random_error_new


  subroutine toggle_verbose
    implicit none
    verbose=.not.verbose
  end   subroutine toggle_verbose





  ! linked

  SUBROUTINE FIND_ORBIT_LAYOUT(RING,FIX,LOC,STATE,TURNS)  ! Finds orbit with TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
    INTEGER, OPTIONAL:: TURNS
    real(dp)  FIX(6),DIX(6),xdix,xdix0,tiny,freq
    TYPE(REAL_8) X(6)
    TYPE(DAMAP) MX,SX,SXI,IS
    integer NO1,ND2,I,IU,LOC,ITE,npara
    TYPE(INTERNAL_STATE),optional, intent(in) :: STATE
    TYPE(INTERNAL_STATE) stat
    TYPE (fibre), POINTER :: C
    logical(lp)  c_da,s_da
    INTEGER TURNS0
    s_da=c_%stable_da
    c_da=c_%check_da
    !   APERTURE=c_%APERTURE_FLAG

    !  c_%APERTURE_FLAG=.false.
    c_%stable_da=.true.
    c_%check_da=.true.
    TURNS0=1
    freq=0.0_dp
    IF(PRESENT(TURNS)) TURNS0=TURNS
    Nullify(C);
    if(.not.ring%closed) then
       !w_p=0
       !w_p%nc=1
       !w_p%fc='((1X,a72))'
       !w_p%c(1)= " This line is not ring : FIND_ORBIT_LAYOUT "
       ! call !write_e(100)
    endif
    dix(:)=0.0_dp
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking
    NO1=1

    if(.not.present(STATE)) then
       IF(default%NOCAVITY) THEN
          !    ND1=2
          stat=default+only_4d-spin0
       ELSE
          !   ND1=3
          STAT=default-spin0
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,a72,/),(1X,a72))'
          !w_p%c(1)=  " No Cavity in the Line "
          !w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          ! call !write_e(111)
       ENDIF
    else
       IF(STATE%NOCAVITY) THEN
          !    ND1=2
          STAT=STATE+only_4d0-spin0
       ELSE
          !   ND1=3
          STAT=STATE-spin0
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,a72,/),(1X,a72))'
          !w_p%c(1)=  " No Cavity in the Line "
          !w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          ! call !write_e(112)
       ENDIF
    endif
101 continue
    !    ND2=2*ND1
    if(stat%totalpath==1.and.(.not.stat%nocavity)) then
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
          IF(stat%TIME) THEN
             XDIX=XDIX+c%mag%P%LD/c%BETA0
          ELSE
             XDIX=XDIX+c%mag%P%LD
          ENDIF
          c=>c%next
          i=i+1
       enddo
       if(freq==0.0_dp) then
          !w_p=0
          !w_p%nc=2
          !w_p%fc='((1X,a72,/),(1X,a72))'
          !w_p%c(1)=  " No Cavity in the Line or Frequency = 0 "
          !w_p%c(2)=  " FIND_ORBIT_LAYOUT will crash "
          ! call !write_e(113)
       endif
       IF(RING%HARMONIC_NUMBER>0) THEN
          FREQ=RING%HARMONIC_NUMBER*CLIGHT/FREQ
          stop 475
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
          !w_p=0
          !w_p%nc=1
          !w_p%fc='((1X,a72))'
          write(6,'(a30,i4)') " Lost in Fixed Point Searcher ",1
          messagelost(len_trim(messagelost)+1:255)=" -> Lost in Fixed Point Searcher "
          ! call ! WRITE_I

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

    xdix=0.0_dp
    do iu=1,ND2
       xdix=abs(dix(iu))+xdix
    enddo

    if(verbose) write(6,'(a22,g21.14)') " Convergence Factor = ",xdix
    !    if(verbose) ! call ! WRITE_I
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
    !  c_%APERTURE_FLAG=APERTURE
    c_%stable_da=s_da

    !  write(6,*) " 2 ",APERTURE,APERTURE_FLAG

  END SUBROUTINE FIND_ORBIT_LAYOUT


  SUBROUTINE FIND_ORBIT_LAYOUT_noda(RING,FIX,LOC,STATE,eps,TURNS) ! Finds orbit without TPSA in State or compatible state
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: RING
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
    freq=0.0_dp
    APERTURE=c_%APERTURE_FLAG
    c_%APERTURE_FLAG=.false.
    messagelost=' Orbit most likely found'
    if(state%radiation) then
     write(6,*) "You have radiation : use find_orbit_x "
     stop
    endif

    if(.not.present(eps)) then
       if(.not.present(STATE)) then
          call FIND_ORBIT_LAYOUT(RING,FIX,LOC,TURNS=TURNS0)
       else
          call FIND_ORBIT_LAYOUT(RING,FIX,LOC,STATE,TURNS=TURNS0)
       endif
       c_%APERTURE_FLAG=APERTURE
       return
    endif


    Nullify(C);

    if(.not.ring%closed) then
       write(6,*) " This line is not ring : FIND_ORBIT_LAYOUT_noda "
        check_stable=.false.
       ! call !write_e(100)
    endif
    dix(:)=0.0_dp
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking
    NO1=1
    if(.not.present(STATE)) then
       IF(default%NOCAVITY) THEN
          !    ND1=2
          stat=default+    only_4d-spin0
       ELSE
          !   ND1=3
          STAT=default-spin0
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          messagelost= " FIND_ORBIT_LAYOUT will crash : exiting"
         check_stable=.false.
          return
       ENDIF
    else
       IF(STATE%NOCAVITY) THEN
          ND2=4
          STAT=STATE+only_4d0-spin0
       ELSE
          ND2=6
          STAT=STATE-spin0
          C=>RING%START
          do i=1,RING%n
             if(C%magp%kind==kind4.OR.C%magp%kind==kind21) goto 101
             C=>C%NEXT
          enddo
          messagelost=" State present; no cavity: FIND_ORBIT_LAYOUT will crash => exiting"
         check_stable=.false.
         return
       ENDIF
    endif
101 continue


    if(stat%totalpath==1.and.(.not.stat%nocavity)) then
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
          XDIX=XDIX+c%mag%P%LD/c%BETA0
          c=>c%next
          i=i+1
       enddo
       if(freq==0.0_dp) then
       
          messagelost= " No Cavity in the Line or Frequency = 0 (totalpath==1)"
         check_stable=.false.
         return
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
       call TRACK(RING,X,LOC,STAT)
       if(.not.check_stable) then
          messagelost(len_trim(messagelost)+1:255)=" -> Unstable tracking: closed orbit not found"
          c_%APERTURE_FLAG=APERTURE
          return
       endif
       !       if(.not.check_stable) then
       !          !w_p=0
       !          !w_p%nc=1
       !          !w_p%fc='((1X,a72))'
       !          write(6,'(a30,i4)') " Lost in Fixed Point Searcher ",2
       !          ! call ! WRITE_I

       !          return
       !       endif

    ENDDO
    !    write(6,*) x
    x(6)=x(6)-freq*turns0

    mx=0.0_dp
    DO J=1,ND2

       Y=FIX
       Y(J)=FIX(J)+EPS
       DO I=1,TURNS0
          CALL TRACK(RING,Y,LOC,STAT)
          if(.not.check_stable) then
             messagelost(len_trim(messagelost)+1:255)=" -> Unstable while tracking small rays around the guessed orbit "
             !   fixed_found=my_false
             c_%APERTURE_FLAG=APERTURE
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

    if(verbose) write(6,*) " Convergence Factor = ",xdix
    !  if(verbose) ! call ! WRITE_I
 !  if(item<10)        then

 !   write(6,*) xdix,xdix0
 !   pause 777
 !    GOTO 3 
 !  endif
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
    endif

    if(ite>=1)  then
       GOTO 3

    endif
    c_%APERTURE_FLAG=APERTURE

  END SUBROUTINE FIND_ORBIT_LAYOUT_noda


  SUBROUTINE fit_all_bends(r,state)
    IMPLICIT NONE
    TYPE(layout),target,INTENT(INOUT):: r
    TYPE(internal_state), intent(in):: state
    type(fibre),pointer :: p
    integer i

    p=>r%start

    do i=1,r%n
       if(p%mag%p%b0/=0.0_dp) call fit_bare_bend(p,state)
       p=>p%next
    enddo

  end SUBROUTINE fit_all_bends

  SUBROUTINE fit_bare_bend(f,state,next)
    IMPLICIT NONE
    TYPE(fibre),INTENT(INOUT):: f
    TYPE(real_8) y(6)
    TYPE(internal_state), intent(in):: state
    !    integer,optional,target :: charge
    real(dp) kf,x(6),xdix,xdix0,tiny
    integer ite
    logical(lp), optional :: next
    logical(lp) nex
    nex=my_false
    if(present(next)) nex=next
    tiny=1e-40_dp
    xdix0=1e4_dp*DEPS_tracking

    KF=0.0_dp   ;
    F%MAGP%BN(1)%KIND=3
    F%MAGP%BN(1)%I=1
    if(nex) then
       F%next%MAGP%BN(1)%KIND=3
       F%next%MAGP%BN(1)%I=1
    endif

    CALL INIT(1,1)

    CALL ALLOC(Y)

3   continue
    X=0.0_dp
    Y=X
    CALL TRACK(f,Y,+state)  !,CHARGE)
    if(nex) CALL TRACK(f%next,Y,+state)  !,CHARGE)
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

    CALL ADD(f,1,1,0.0_dp)     !etienne
    if(nex) CALL ADD(f%next,1,1,0.0_dp)     !etienne

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
    TYPE(INTERNAL_STATE), intent(IN):: my_STATE
    TYPE(INTERNAL_STATE) STATE
    TYPE(layout),target, intent(inout) :: R
    REAL(DP), pointer :: BETA(:,:,:)
    integer pos,nturn,i,flag,ib,MF,mft,j,resmax,it,I1,no
    real(dp) closed(6),MAT(6,6),AI(6,6),A(6,6),emit(2),emit0(6),aper(2),x(6),xn(6),dbeta,tuneold(:)
    real(dp) ra(2),tunenew(2),xda(lnv)
    CHARACTER(*) FILENAME,FILEtune,FILESMEAR
    real(dp), allocatable :: dat(:,:),dats(:,:),SMEAR(:,:)
    REAL(DP) JMin(2),JMAX(2), tune1(2),tune2(2),tot_tune(2),epsi,scas(2),scau,scat(2)
    integer itmax
    type(damap) id
    type(tree) monkey



    epsi=1.0_dp/nturn
    STATE=((((my_state+nocavity0)+delta0)+only_4d0)-RADIATION0)
    allocate(dat(0:nturn,6),dats(0:nturn,6))
    allocate(SMEAR(ITMAX,8))
    CLOSED=0.0_dp

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
       xda=0.0_dp
    endif
    scau=1.0_dp
    scas=0.0_dp
    dats=0.0_dp
    SMEAR=0.0_dp
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
    dat=0.0_dp

    xn=0.0_dp
    JMAX=0.0_dp
    JMIN=mybig
    emit=scat*aper
    write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    write(6,*) " Initial emit = ", emit(1:2)," scale = ",scat
    write(6,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    xn(2)=sqrt(emit(1))
    xn(4)=sqrt(emit(2))
    X=0.0_dp

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
    WRITE(6,*) "     SMEAR = ",2.0_dp*(JMAX-JMIN)/(JMAX+JMIN)
    SMEAR(IT,1:2)=EMIT(1:2)
    SMEAR(IT,3:4)=JMIN(1:2)
    SMEAR(IT,5:6)=JMAX(1:2)
    SMEAR(IT,7:8)=2.0_dp*(JMAX-JMIN)/(JMAX+JMIN)
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
       tot_tune=0.0_dp
       xn(1:4)=dats(0,1:4)
       tune1(1)=atan2(-xn(2),xn(1))/twopi
       tune1(2)=atan2(-xn(4),xn(3))/twopi
       if(tune1(1)<0.0_dp)  tune1(1)=tune1(1)+1.0_dp
       if(tune1(2)<0.0_dp)  tune1(2)=tune1(2)+1.0_dp
       DO I1=0,NTURN
          xn(1:4)=dats(i1,1:4)
          tune2(1)=atan2(-xn(2),xn(1))/twopi
          tune2(2)=atan2(-xn(4),xn(3))/twopi
          if(tune2(1)<0.0_dp)  tune2(1)=tune2(1)+1.0_dp
          if(tune2(2)<0.0_dp)  tune2(2)=tune2(2)+1.0_dp
          tune1=tune2-tune1
          if(tune1(1)<0.0_dp)  tune1(1)=tune1(1)+1.0_dp
          if(tune1(2)<0.0_dp)  tune1(2)=tune1(2)+1.0_dp
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
    tot_tune=0.0_dp
    xn(1:4)=dats(0,1:4)
    tune1(1)=atan2(-xn(2),xn(1))/twopi
    tune1(2)=atan2(-xn(4),xn(3))/twopi
    if(tune1(1)<0.0_dp)  tune1(1)=tune1(1)+1.0_dp
    if(tune1(2)<0.0_dp)  tune1(2)=tune1(2)+1.0_dp

    DO I=0,NTURN
       WRITE(MF,200)DATs(I,1:6)
       xn(1:4)=dats(i,1:4)
       tune2(1)=atan2(-xn(2),xn(1))/twopi
       tune2(2)=atan2(-xn(4),xn(3))/twopi
       if(tune2(1)<0.0_dp)  tune2(1)=tune2(1)+1.0_dp
       if(tune2(2)<0.0_dp)  tune2(2)=tune2(2)+1.0_dp
       tune1=tune2-tune1
       if(tune1(1)<0.0_dp)  tune1(1)=tune1(1)+1.0_dp
       if(tune1(2)<0.0_dp)  tune1(2)=tune1(2)+1.0_dp
       tot_tune =tot_tune+tune1
       tune1=tune2
    ENDDO
    tot_tune=tot_tune/nturn
    CLOSE(MF)

    CLOSE(mft)

    CALL KANALNUMMER(MF)
    OPEN(UNIT=MF,FILE=FILESMEAR)
    WRITE(MF,*) " ITERATION   EMIT0(1:2)  JMIN(1:2) JMAX(1:2) 2.0_dp*(JMAX-JMIN)/(JMAX+JMIN)"
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
    WRITE(6,*) "     SMEAR = ",2.0_dp*(JMAX-JMIN)/(JMAX+JMIN)
    write(6,*) " *****************************************************************"
    if(pos==0) call kill(monkey)
  end SUBROUTINE  track_aperture

 subroutine point_m_u(m1,m2)
 implicit none
 type(mad_universe), target :: m1,m2

  m_u=>m1
  m_t=>m2

 end  subroutine point_m_u

  SUBROUTINE  THIN_LENS_resplit(R,THIN,even,lim,limit_wiggler,lmax0,xbend,sexr,fib,useknob,universe) ! A re-splitting routine
    IMPLICIT NONE
    INTEGER NTE
    TYPE(layout),target, intent(inout) :: R
    real(dp), OPTIONAL, intent(inout) :: THIN
    real(dp), OPTIONAL, intent(in) :: lmax0
    real(dp), OPTIONAL, intent(in) ::xbend
    real(dp), OPTIONAL, intent(in) ::sexr
    type(fibre), OPTIONAL, target :: fib
    logical(lp), OPTIONAL :: useknob,universe
    integer, optional :: limit_wiggler(2),lim(2)
    real(dp) gg,RHOI,XL,QUAD,THI,lm,dl,ggbt,xbend1,gf(7),sexr0,quad0,dq
    INTEGER M1,M2,M3, MK1,MK2,MK3,limit(2),parity,inc,nst_tot,ntec,ii,metb,sexk
    integer incold ,parityold, nt,nsag,lim0(2),lims(2),kkk

    logical(lp) MANUAL,eject,doit,DOBEND
    TYPE (fibre), POINTER :: C
    logical(lp),optional :: even
    type(layout), pointer :: L
    logical f1,f2,m_t_pres,uni
    
!!!! gymnastic for M_u and M_t (these are local m_u and m_t)
     kkk=lielib_print(12) 
     lielib_print(12)=0

m_t_pres=my_false
if(associated(m_u)) then
 if(associated(m_t%start)) then
  m_t_pres=associated(m_t%start%start) 
 endif
endif
    uni=.false.
    if(present(universe)) uni=universe

    if(uni) then
      if(r%parent_universe%nf==0) call TIE_MAD_UNIVERSE(r%parent_universe)
      nt=r%parent_universe%nf
    else
      nt=r%n
    endif

    sexr0=0.0_dp


    if(present(sexr)) sexr0=sexr

    f1=.false.
    f2=.false.

if(m_t_pres) then
   l=>m_t%start
       do ii=1,m_t%n
          call kill(l%t)
          l=>l%next
       enddo
   l=>m_u%start
       do ii=1,m_u%n
          call kill(l%t)
          l=>l%next
       enddo
else
    if(associated(r%parent_universe)) then
       f1=.true.
       l=>r%parent_universe%start
       do ii=1,r%parent_universe%n
          call kill(l%t)
          l=>l%next
       enddo
    elseif(associated(r%t)) then
       call kill(r%t)
       f2=.true.
    endif
endif
    !    logical(lp) doneit
    nullify(C)
    parity=0
    inc=0
    lm=1.0e38_dp
    ntec=0
    max_ds=0.0_dp
    xbend1=-1.0_dp

    if(present(xbend)) xbend1=xbend
    if(present(lmax0)) lm=abs(lmax0)
    if(present(even)) then
       inc=1
       if(even) then
          parity=0
       else
          parity=1
       endIf
    endif
    if(sagan_even.and.parity==1) then
     inc=0
     parity=0
    endif
    parityold=parity
    incold=inc
    !   CALL LINE_L(R,doneit)

    MANUAL=.FALSE.
    eject=.FALSE.

    THI=R%THIN
    IF(PRESENT(THIN)) THI=THIN

    IF(THI<=0) MANUAL=.TRUE.


    IF(MANUAL) THEN
       write(6,*) "thi: thin lens factor (THI<0 TO STOP), sextupole factor and Bend factor "
       read(5,*) thi,sexr0,xbend1
       IF(THI<0) eject=.true.
    ENDIF

1001 CONTINUE

    limit= limit_int0 
 
    if(present(lim)) limit=lim
    lims=limit
    if(sixtrack_compatible) then
       limit(1)=1000000
       limit(2)=1000001
    endif
    lim0=limit
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
    sexk=0
    if(uni) then
     c=>r%parent_universe%start%start
    else
     C=>R%START
    endif
    do  ii=1,nt    ! WHILE(ASSOCIATED(C))
       doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
       doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
       DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
       DOIT=DOIT.OR.(C%MAG%KIND==kind17.or.C%MAG%KIND==kindwiggler.or.C%MAG%KIND==KINDhel)
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

       if(uni) then
        c=>c%n
       else
        C=>C%NEXT
       endif
    enddo
if(lielib_print(14)==1) then
    write(6,*) "Previous of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3
    write(6,*)   "Total NST ", NST_tot
endif
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
    ggbt=0.0_dp
    r%NTHIN=0
    r%THIN=THI

    nst_tot=0
    if(uni) then
     c=>r%parent_universe%start%start
    else
     C=>R%START
    endif
    do  ii=1,nt   ! WHILE(ASSOCIATED(C))
       doit=.true.
       if(present(useknob)) then
          if(useknob) then
             doit=c%magp%knob
          else
             doit=.not.c%magp%knob
          endif
       endif
       if(present(fib)) then
          doit=associated(c,fib)
       endif
       if(.not.doit) then
          NST_tot=NST_tot+C%MAG%P%nst
          if(c%mag%even) then
             parity=parityold
             inc=incold
          endif
             if(uni) then
              c=>c%n
             else
              C=>C%NEXT
             endif
          cycle
       endif

       if(c%mag%even) then
          parity=0
          inc=1
       endif
             if(C%MAG%KIND==kindwiggler)  then
                if(present(limit_wiggler)) then
                 limit=limit_wiggler
                else
                 limit=limit_sag
                endif
             else
               limit=lims
             endif
       !       if(doit)  then

       select case(resplit_cutting)

       case(0)

          doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
          doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
          DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
          DOIT=DOIT.OR.(C%MAG%KIND==kind17.or.C%MAG%KIND==kindwiggler.or.C%MAG%KIND==KINDhel)
          doit=doit.and.C%MAG%recut
          if(doit) then
             xl=C%MAG%L
             RHOI=0.0_dp
             QUAD=0.0_dp
             QUAD0=0.0_dp
             IF(C%MAG%P%NMUL>=1) THEN
                !               RHOI=C%MAG%P%B0
                RHOI=abs(C%MAG%bn(1))+abs(C%MAG%an(1))
             endif
             IF(C%MAG%P%NMUL>=2) THEN
                QUAD=SQRT(C%MAG%BN(2)**2+C%MAG%AN(2)**2)
                IF(C%MAG%P%NMUL>=3) THEN
                 if(sexr0>0) quad0=SQRT(C%MAG%BN(3)**2+C%MAG%AN(3)**2)*sexr0
                 QUAD=QUAD+quad0
                endif
             ELSE
                QUAD=0.0_dp
             ENDIF
             if(C%MAG%KIND==kind5.or.C%MAG%KIND==kind17) then
                quad=quad+(C%MAG%b_sol)**2/4.0_dp !+abs(C%MAG%b_sol/2.0_dp)
             endif
             if(C%MAG%KIND==kind10) then
               if(associated(c%mag%tp10%ae).and.associated(c%mag%tp10%be)) then
                quad=quad+SQRT(c%mag%tp10%ae(2)**2+c%mag%tp10%be(2)**2)*volt_c/c%mag%p%p0c  
              endif
             endif
             if(C%MAG%KIND==kindwiggler) then
               call eval_thin_q(C%MAG%wi,dQ,nsag)
               quad=quad+dq
     !          write(6,*) " wiggler detected ",dq,quad,nsag
             else
               nsag=1
             endif
             DOBEND=MY_FALSE
             IF(xbend1>0.0_dp) THEN
                IF(C%MAG%KIND==kind10) THEN
                   IF(C%MAG%TP10%DRIFTKICK) THEN
                      DOBEND=MY_TRUE
                   ENDIF
                ENDIF
                IF(C%MAG%KIND==kind16.OR.C%MAG%KIND==kind20) THEN
                   IF(C%MAG%K16%DRIFTKICK) THEN
                      DOBEND=MY_TRUE
                   ENDIF
                ENDIF
                if(rhoi/=0.0_dp.and.radiation_bend_split)DOBEND=MY_TRUE
             ENDIF
             !  ETIENNE
             GG=XL*(RHOI**2+ABS(QUAD))
             GG=GG/THI
             NTE=INT(GG)

             sexk=sexk+xl*ABS(QUAD0)/thi
             metb=0
             if(dobend) then
                call check_bend(xl,gg,rhoi,xbend1,gf,metb)
                if(gf(metb)>gg) then
                   gg=gf(metb)
                   NTE=INT(GG)
                   ggbt=ggbt+NTE
                else
                   metb=0
                endif
             endif

             if(C%MAG%KIND==kindwiggler) then
              if(nsag>nte) nte=nsag
             endif 

             IF(NTE.LT.limit(1).or.metb==2) THEN
                M1=M1+1
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=2
                MK1=MK1+NTE
             ELSEIF((NTE.GE.limit(1).AND.NTE.LT.limit(2)).or.metb==4) THEN
                M2=M2+1
                NTE=NTE/3
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=4
                MK2=MK2+NTE*3
             ELSEIF(NTE.GE.limit(2).or.metb==6) THEN
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
             call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
             call COPY(C%MAG,C%MAGP)
             if(gg>0.0_dp) then
                if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst
             endif

          endif  ! doit

       case(1)

          doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
          doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
          DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
          DOIT=DOIT.OR.(C%MAG%KIND==kind17.or.C%MAG%KIND==kindwiggler.or.C%MAG%KIND==KINDhel)
          doit=doit.and.C%MAG%recut

          if(doit) then
             xl=C%MAG%L
             RHOI=0.0_dp
             QUAD=0.0_dp
             QUAD0=0.0_dp
             IF(C%MAG%P%NMUL>=1) THEN
                !               RHOI=C%MAG%P%B0
                RHOI=abs(C%MAG%bn(1))+abs(C%MAG%an(1))
             endif
             IF(C%MAG%P%NMUL>=2) THEN
                QUAD=SQRT(C%MAG%BN(2)**2+C%MAG%AN(2)**2)
                IF(C%MAG%P%NMUL>=3) THEN
                 if(sexr0>0) quad0=SQRT(C%MAG%BN(3)**2+C%MAG%AN(3)**2)*sexr0
                 QUAD=QUAD+quad0
                endif
             ELSE
                QUAD=0.0_dp
             ENDIF
             if(C%MAG%KIND==kind5.or.C%MAG%KIND==kind17) then
                quad=quad+(C%MAG%b_sol)**2/4.0_dp   !+abs(C%MAG%b_sol/2.0_dp)
             endif
             if(C%MAG%KIND==kind10) then
               if(associated(c%mag%tp10%ae).and.associated(c%mag%tp10%be)) then
                quad=quad+SQRT(c%mag%tp10%ae(2)**2+c%mag%tp10%be(2)**2)*volt_c/c%mag%p%p0c  
              endif
             endif
             if(C%MAG%KIND==kindwiggler) then
               call eval_thin_q(C%MAG%wi,dQ,nsag)
               quad=quad+dq
        !       write(6,*) " wiggler detected ",dq,quad,nsag
             else
               nsag=1
             endif
             DOBEND=MY_FALSE
             IF(xbend1>0.0_dp) THEN
                IF(C%MAG%KIND==kind10) THEN
                   IF(C%MAG%TP10%DRIFTKICK) THEN
                      DOBEND=MY_TRUE
                   ENDIF
                ENDIF
                IF(C%MAG%KIND==kind16.OR.C%MAG%KIND==kind20) THEN
                   IF(C%MAG%K16%DRIFTKICK) THEN
                      DOBEND=MY_TRUE
                   ENDIF
                ENDIF
                if(rhoi/=0.0_dp.and.radiation_bend_split)DOBEND=MY_TRUE
             ENDIF
             !  ETIENNE
             GG=XL*(RHOI**2+ABS(QUAD))
             GG=GG/THI
             sexk=sexk+xl*ABS(QUAD0)/thi
             NTE=INT(GG)
             metb=0
             if(dobend) then
                call check_bend(xl,gg,rhoi,xbend1,gf,metb)
                if(gf(metb)>gg) then
                   gg=gf(metb)
                   NTE=INT(GG)
                   ggbt=ggbt+NTE
                else
                   metb=0
                endif
             endif

             if(C%MAG%KIND==kindwiggler) then
              if(nsag>nte) nte=nsag
             endif 

             IF(NTE.LT.limit(1).or.metb==2) THEN
                M1=M1+1
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=2
                MK1=MK1+NTE
             ELSEIF((NTE.GE.limit(1).AND.NTE.LT.limit(2)).or.metb==4) THEN
                M2=M2+1
                NTE=NTE/3
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=4
                MK2=MK2+NTE*3
             ELSEIF(NTE.GE.limit(2).or.metb==6) THEN
                M3=M3+1
                NTE=NTE/7
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=6
                MK3=MK3+NTE*7
             ENDIF


             r%NTHIN=r%NTHIN+1  !C%MAG%NST

             if(present(lmax0).and.c%mag%kind==kind1) then
                dl=(C%MAG%P%ld/C%MAG%P%nst)
                if(dl>lm*fuzzy_split) then
                   ntec=int(C%MAG%P%ld/lm)+1
                   if(mod(ntec,2)/=parity) ntec=ntec+inc
                   C%MAG%P%NST=ntec
                endif
             endif

             call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
             call COPY(C%MAG,C%MAGP)
             !             if(gg>zero) then
             if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst
             !             endif
          endif

       case(2,-2)

          doit=.not.(resplit_cutting==-2.and.C%MAG%KIND==kind1)

          doit=(C%MAG%KIND==kind1.or.C%MAG%KIND==kind2.or.C%MAG%KIND==kind4.or.C%MAG%KIND==kind5)
          doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
          DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
          DOIT=DOIT.OR.(C%MAG%KIND==kind17.or.C%MAG%KIND==kindwiggler.or.C%MAG%KIND==KINDhel)
          doit=doit.and.C%MAG%recut

          if(doit) then
             xl=C%MAG%L
             RHOI=0.0_dp
             QUAD=0.0_dp
             QUAD0=0.0_dp
             IF(C%MAG%P%NMUL>=1) THEN
                !               RHOI=C%MAG%P%B0
                RHOI=abs(C%MAG%bn(1))+abs(C%MAG%an(1))
             endif
             IF(C%MAG%P%NMUL>=2) THEN
                QUAD=SQRT(C%MAG%BN(2)**2+C%MAG%AN(2)**2)
                IF(C%MAG%P%NMUL>=3) THEN
                 if(sexr0>0) quad0=SQRT(C%MAG%BN(3)**2+C%MAG%AN(3)**2)*sexr0
                 QUAD=QUAD+quad0
                endif
             ELSE
                QUAD=0.0_dp
             ENDIF
             if(C%MAG%KIND==kind5.or.C%MAG%KIND==kind17) then
                quad=quad+(C%MAG%b_sol)**2/4.0_dp   !+abs(C%MAG%b_sol/2.0_dp)
             endif
             if(C%MAG%KIND==kind10) then
               if(associated(c%mag%tp10%ae).and.associated(c%mag%tp10%be)) then
                quad=quad+SQRT(c%mag%tp10%ae(2)**2+c%mag%tp10%be(2)**2)*volt_c/c%mag%p%p0c  
              endif
             endif
             if(C%MAG%KIND==kindwiggler) then
               call eval_thin_q(C%MAG%wi,dQ,nsag)
               quad=quad+dq
          !     write(6,*) " wiggler detected ",dq,quad,nsag
             else
               nsag=1
             endif

             DOBEND=MY_FALSE
             IF(xbend1>0.0_dp) THEN
                IF(C%MAG%KIND==kind10) THEN
                   IF(C%MAG%TP10%DRIFTKICK) THEN
                      DOBEND=MY_TRUE
                   ENDIF
                ENDIF
                IF(C%MAG%KIND==kind16.OR.C%MAG%KIND==kind20) THEN
                   IF(C%MAG%K16%DRIFTKICK) THEN
                      DOBEND=MY_TRUE
                   ENDIF
                ENDIF
                if(rhoi/=0.0_dp.and.radiation_bend_split)DOBEND=MY_TRUE
             ENDIF
             !  ETIENNE
             GG=XL*(RHOI**2+ABS(QUAD))
             GG=GG/THI
             NTE=INT(GG)
             sexk=sexk+xl*ABS(QUAD0)/thi
             metb=0
             if(dobend) then
                call check_bend(xl,gg,rhoi,xbend1,gf,metb)
                if(gf(metb)>gg) then
                   gg=gf(metb)
                   NTE=INT(GG)
                   ggbt=ggbt+NTE
                else
                   metb=0
                endif
             endif

             if(C%MAG%KIND==kindwiggler) then
              if(nsag>nte) nte=nsag
             endif 

             IF(NTE.LT.limit(1).or.metb==2) THEN
                M1=M1+1
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=2
                MK1=MK1+NTE
             ELSEIF((NTE.GE.limit(1).AND.NTE.LT.limit(2)).or.metb==4) THEN
                M2=M2+1
                NTE=NTE/3
                IF(NTE.EQ.0) NTE=1
                if(mod(nte,2)/=parity) nte=nte+inc
                C%MAG%P%NST=NTE
                C%MAG%P%METHOD=4
                MK2=MK2+NTE*3
             ELSEIF(NTE.GE.limit(2).or.metb==6) THEN
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
             if(nte>ntec.or.(.not.present(lmax0)) ) then
                call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
                call COPY(C%MAG,C%MAGP)
             endif
             !            if(gg>zero) then
             !               if(c%mag%l/c%mag%p%nst>max_ds) max_ds=c%mag%l/c%mag%p%nst
             !            endif

             if(present(lmax0)) then
                dl=(C%MAG%P%ld/C%MAG%P%nst)
                if(dl>lm*fuzzy_split.and.C%MAG%KIND/=kindpa) then
                   nte=int(C%MAG%P%ld/lm)+1
                   if(mod(nte,2)/=parity) nte=nte+inc
                   if(nte > C%MAG%P%NST ) then
                      C%MAG%P%NST=nte
                      call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
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
 
             if(C%MAG%KIND==kindwiggler)  then
                if(present(limit_wiggler)) then
                 limit=limit_wiggler
                else
                 limit=limit_sag
                endif
             else
               limit=lims
             endif             

       !      endif
       NST_tot=NST_tot+C%MAG%P%nst
       if(c%mag%even) then
          parity=parityold
          inc=incold
       endif
             if(uni) then
              c=>c%n
             else
              C=>C%NEXT
             endif

    enddo   !   end of do   WHILE

if(lielib_print(14)==1) then
    write(6,*) "Present of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3
    write(6,*)   "Total NST ", NST_tot
endif
    if(radiation_bend_split) then
if(lielib_print(14)==1) then
       write(6,*)   "Total NST due to Bend Closed Orbit ", int(ggbt)
       write(6,*)   "Restricted to method=2 for radiation or spin "
endif
    else
if(lielib_print(14)==1) then
       write(6,*)   "Total NST due to Bend Closed Orbit ", int(ggbt)
endif
    endif
if(lielib_print(14)==1) then
    write(6,*)   "Total NST due to Sextupoles ", sexk
    write(6,*)   "Biggest ds ", max_ds
endif


    IF(MANUAL) THEN
       write(6,*) "thi: thin lens factor (THI<0 TO STOP), sextupole factor and Bend factor "
       read(5,*) thi,sexr0, xbend1
       IF(THI<0) THEN
          THI=R%THIN
          !          limit(1)=limit0(1)
          !          limit(2)=limit0(2)
          if(f1) then
             l=>r%parent_universe%start
             do ii=1,r%parent_universe%n
                call make_node_layout(l)
                l=>l%next
             enddo
          elseif(f2) then
             call make_node_layout(r)  !!! bug (l) was wrong
          endif
          RETURN
       ELSE
          GOTO 1001
       ENDIF
    else
          if(f1) then
             l=>r%parent_universe%start
             do ii=1,r%parent_universe%n
                call make_node_layout(l)
                l=>l%next
             enddo
          elseif(f2) then
             call make_node_layout(r)  !!! bug (l) was wrong
          endif

    ENDIF

     lielib_print(12)=kkk
    !    limit(1)=limit0(1)
    !    limit(2)=limit0(2)

    !    CALL RING_L(R,doneit)

  END SUBROUTINE  THIN_LENS_resplit



  SUBROUTINE  RECUT_KIND7_one(C,lmax0,drift,ido,idc) ! A re-splitting routine
    IMPLICIT NONE
    TYPE(layout),POINTER :: L
    real(dp) lmax0
    TYPE(FIBRE),target :: C
    INTEGER I,f0,ido,idc
    logical(lp) drift,doit
   
    if(associated(c%parent_layout)) then
       if(associated(c%parent_layout%parent_universe)) then
          l=>c%parent_layout%parent_universe%start
          do i=1,c%parent_layout%parent_universe%n
             call kill(l%t)
             l=>l%next
          enddo
       else
          call kill(c%parent_layout%t)
       endif
    endif
    if(drift.and.C%MAG%KIND==KIND1) then
       f0=nint(C%MAG%l/lmax0)
       if(f0==0) f0=1
       C%MAG%p%nst=f0
       C%MAGp%p%nst=C%MAG%p%nst
       call COPY(C%MAG,C%MAGP)
    endif
    if(mod(C%MAG%p%method,2)==0) then  !!!
       doit=C%MAG%KIND==KIND7.or.(C%MAG%KIND==KIND2.and.C%MAG%p%method==2)
       if(associated(C%MAG%K16)) then
          doit=(C%MAG%K16%DRIFTKICK.and.C%MAG%p%method==2)
       endif
       if(associated(C%MAG%TP10)) then
          doit=(C%MAG%TP10%DRIFTKICK.and.C%MAG%p%method==2)
       endif
       IF(doit) then
          ido=ido+1
          !         f0=nint(C%MAG%l/lmax0)
          f0=nint(C%MAG%l/lmax0/C%MAG%p%nst/2)
          if(C%MAG%p%method==6) f0=nint(C%MAG%l/lmax0/C%MAG%p%nst/4)
          if(f0==0) f0=1
          if(C%MAG%p%method==2) then
             C%MAG%p%nst=C%MAG%p%nst*f0*2
             C%MAGp%p%nst=C%MAG%p%nst
             C%MAG%p%method=1
             C%MAGp%p%method=1
          elseif(C%MAG%p%method==4) then
             C%MAG%p%nst=C%MAG%p%nst*f0*2
             C%MAGp%p%nst=C%MAG%p%nst
             C%MAG%p%method=3
             C%MAGp%p%method=3
          elseif(C%MAG%p%method==6) then
             C%MAG%p%nst=C%MAG%p%nst*f0*4
             C%MAGp%p%nst=C%MAG%p%nst
             C%MAG%p%method=5
             C%MAGp%p%method=5
          endif
          call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
          call COPY(C%MAG,C%MAGP)
          if(C%MAG%KIND==KIND7) then
             C%MAG%t7%f=f0
             C%MAGp%t7%f=f0
          elseif(associated(C%MAG%K16)) then
             C%MAG%K16%f=f0
             C%MAGp%K16%f=f0
          elseif(associated(C%MAG%TP10)) then
             C%MAG%TP10%f=f0
             C%MAGp%TP10%f=f0
          else
             C%MAG%k2%f=f0
             C%MAGp%k2%f=f0
          endif
       ENDIF
    else !!!
       idc=idc+1
       f0=nint(C%MAG%l/lmax0/C%MAG%p%nst)
       if(f0>=1) then
          C%MAG%p%nst=C%MAG%p%nst*f0
          C%MAGp%p%nst=C%MAG%p%nst
          call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
          call COPY(C%MAG,C%MAGP)
          if(C%MAG%KIND==KIND7) then
             C%MAG%t7%f=f0*C%MAG%t7%f
             C%MAGp%t7%f=C%MAG%t7%f
          elseif(associated(C%MAG%K16)) then
             C%MAG%K16%f=f0*C%MAG%K16%f
             C%MAGp%K16%f=C%MAG%K16%f
          elseif(associated(C%MAG%TP10)) then
             C%MAG%TP10%f=f0*C%MAG%TP10%f
             C%MAGp%TP10%f=C%MAG%TP10%f
          else
             C%MAG%k2%f=f0*C%MAG%k2%f
             C%MAGp%k2%f=C%MAG%k2%f
          endif
       endif
    endif !!!
  END SUBROUTINE  RECUT_KIND7_one


  SUBROUTINE  RECUT_KIND7(R,lmax0,drift) ! A re-splitting routine
    IMPLICIT NONE
    TYPE(layout),target, intent(inout) :: R
    real(dp) lmax0
    TYPE(FIBRE),POINTER :: C
    INTEGER I,ido,idc
    logical(lp) drift

    ido=0
    idc=0


    C=>R%START
    DO I=1,R%N
       call RECUT_KIND7_one(c,lmax0,drift,ido,idc)
       C=>C%NEXT
    ENDDO
    write(6,*) ido," elements changed to odd methods " 
    write(6,*) idc," elements only " 
  END SUBROUTINE  RECUT_KIND7



subroutine dipole_check_step(m,n,L,b0,k1,x)
implicit none
integer m,n,i,j
real(dp) L,b0,k1,h,x(2)
real(dp) DH,D,D1,D2,DK1,DK2,DF(4),DK(4)

x=0
 
h=b0

    SELECT CASE(m)
 
    CASE(2)
       DH=L/2.0_dp/n
       D=L/n
do i=1,n
       call prot_drift(dh,X,h)
       CALL SKICKt(x,D,b0,k1,h)
       call prot_drift(dh,X,h)
enddo
    CASE(4)

       D1=L*FD1/n
       D2=L*FD2/n
       DK1=L*FK1/n
       DK2=L*FK2/n
do i=1,n
       call prot_drift(D1,X,h)
       CALL SKICKt (x,DK1,b0,k1,h)
       call prot_drift(D2,X,h)
       CALL SKICKt (x,DK2,b0,k1,h)

       call prot_drift(D2,X,h)
 
       CALL SKICKt (x,DK1,b0,k1,h)
       call prot_drift(D1,X,h)
enddo

    CASE(6)
       DO I =1,4
          DF(I)=L*YOSD(I)/n
          DK(I)=L*YOSK(I)/n
       ENDDO
do i=1,n
       DO J=4,2,-1
       call prot_drift(Df(j),X,h)
       CALL SKICKt (x,DK(J),b0,k1,h)
       ENDDO
 
       call prot_drift(Df(1),X,h)
       CALL SKICKt (x,DK(1),b0,k1,h)

       call prot_drift(Df(1),X,h)
       DO J=2,4
          CALL SKICKt (x,DK(J),b0,k1,h)
       call prot_drift(Df(j),X,h)
       ENDDO
enddo


    CASE DEFAULT
       !w_p=0
       !w_p%nc=1
       !w_p%fc='(1(1X,A72))'
         write(6,'(a12,1x,i4,1x,a17)') " THE METHOD ",m," IS NOT SUPPORTED"
       ! call !write_e(357)
    END SELECT

end subroutine dipole_check_step

subroutine SKICKt(x,yL,b0,k1,h)
implicit none
real(dp) yL,b0,k1,h,x(2)


x(2)=x(2)-yL*(b0*(1.0_dp+h*x(1))+k1*x(1))

end subroutine SKICKt

  SUBROUTINE prot_drift(YL,X,h)
    IMPLICIT NONE
    real(dp) X(2)
    real(dp) YL,h
 
    real(dp) XN(2),PZ,pt
    real(dp)  A,R
           PZ=sqrt(1.0_dp-X(2)**2)
    if(h/=0.0_dp) then
       A=YL*h
       R=1.0_dp/h
 
 

          PT=1.0_dp-X(2)*TAN(A)/PZ
          !       XN(1)=(X(1)+R)/COS(A)/PT-R
          XN(1)=(X(1)+R*(2.0_dp*sin(a/2.0_dp)**2+X(2)*sin(A)/PZ))/COS(A)/PT
          XN(2)=X(2)*COS(A)+SIN(A)*PZ
 
    
       X(1)=XN(1)
       X(2)=XN(2)
    else
     x(1)=x(1)+x(2)*yl/pz
    endif

  END SUBROUTINE prot_drift



  SUBROUTINE  check_bend(xl,ggi,rhoi,xbend1,gf,met) ! A re-splitting routine
    IMPLICIT NONE
    real(dp) xl,gg,ggi,rhoi,ggb,ar,co(7),xbend1,gf(7)
    integer i,met

    gg=int(ggi)
    if(gg==0.d0) gg=1.d0
    co(2)=1.d0/12.0_dp
    co(4)=  0.17e0_dp
    co(6)=0.17e-1_dp

    gf=0.0_dp
 
    do i=3,7,2
           ar=i-1
       gf(i)= (co(i-1)/xbend1)**(1.d0/ar)*abs(rhoi)*xl 
       if(gf(i)<gg) gf(i)=gg
    enddo
    
    gf(2)=gf(3)
    gf(4)=gf(5)*3
    gf(6)=gf(7)*7
    

    met=2

    if(gf(4)<gf(2)) met=4
    if(gf(6)<gf(4).and.gf(6)<gf(2)) met=6
    if(radiation_bend_split) met=2
    if(sixtrack_compatible) met=2
    

  end SUBROUTINE  check_bend
! routine which computed check_bends co(i)
!call init(15,2,2,0,.true.)
!call alloc(h1,h2,h)
!call alloc(id,idh,idi)
!call alloc(onel)
!call alloc(dl,del,bb,t1,t2)
!del=1.d0.mono.3
!dl=1.d0.mono.6
!bb=1.d0.mono.5
!id=1
!t1=(-(1.d0+bb*id%v(1))*sqrt((1.d0+del)**2-id%v(2)**2)+del)*dl
!t2=+bb*(id%v(1)+id%v(1)**2/2.d0*bb)*dl
!nn=1
!h%h=t1+t2

!h1%h=t1/2.d0
!h2%h=t2
!idi=exp(h1,id)
!idi=exp(h2,idi)
!idi=exp(h1,idi)
!onel=idi
!
!idh=exp(h,id)
!
!idi%v(nn)=idi%v(nn)-idh%v(nn)
!call print(idi%v(nn),6)
!pause 1
!
!h%h=t1+t2
!id=1
!h1%h=fd1*t1
!idi=exp(h1,id) 
!h2%h=fk1*t2
!idi=exp(h2,idi) 
!h1%h=fd2*t1    
!idi=exp(h1,idi) 
!h2%h=fk2*t2 
!idi=exp(h2,idi) 
! 
!h1%h=fd2*t1    
!idi=exp(h1,idi) 
!h2%h=fk1*t2
!idi=exp(h2,idi) 
!h1%h=fd1*t1
!idi=exp(h1,idi) 
!onel=idi 
! 
!idh=exp(h,id)

!idi%v(nn)=idi%v(nn)-idh%v(nn)
!call print(idi%v(nn),6)
!pause 2
!h%h=t1+t2
!id=1
!
!call MAKE_YOSHIDA
!
!       DO i=4,2,-1
!        h1%h=YOSD(I)*t1
!       if(i==4) then
!        idi=exp(h1,id) 
!       else
!        idi=exp(h1,idi) 
!       endif
!       h2%h=YOSK(I)*t2
!       idi=exp(h2,idi) 
!       ENDDO
!        h1%h=YOSD(1)*t1
!       idi=exp(h1,idi) 
!       h2%h=YOSK(1)*t2
!       idi=exp(h2,idi) 
!        h1%h=YOSD(1)*t1
!       idi=exp(h1,idi) 
!       DO i=2,4
!       h2%h=YOSK(I)*t2
!       idi=exp(h2,idi) 
!        h1%h=YOSD(I)*t1
!        idi=exp(h1,idi) 
!       ENDDO
!onel=idi 
! 
!idh=exp(h,id)
!
!idi%v(nn)=idi%v(nn)-idh%v(nn)
!call print(idi%v(nn),6)

  SUBROUTINE  THIN_LENS_restart(R,fib,useknob,universe,ignore_recut) ! A re-splitting routine
    IMPLICIT NONE
    INTEGER NTE
    TYPE(layout),target, intent(inout) :: R
    INTEGER M1,M2,M3, MK1,MK2,MK3,nst_tot,ii,nt  !,limit0(2)
    type(fibre), OPTIONAL, target :: fib
    logical(lp), OPTIONAL :: useknob,universe,ignore_recut
    logical(lp) doit,uni,m_t_pres
    TYPE (fibre), POINTER :: C
    TYPE (layout), POINTER :: l
     logical(lp) ignore 
    nullify(C)
ignore=.false.
if(present(ignore_recut)) ignore=ignore_recut
    !    CALL LINE_L(R,doneit)
m_t_pres=my_false
if(associated(m_u)) then
 if(associated(m_t%start)) then
  m_t_pres=associated(m_t%start%start) 
 endif
endif
  
!  if(associated(r%parent_universe)) then
!       l=>r%parent_universe%start
!       do ii=1,r%parent_universe%n
!          call kill(l%t)
!          l=>l%next
!       enddo
!    else
!       call kill(r%t)
!    endif





if(m_t_pres) then
   l=>m_t%start
       do ii=1,m_t%n
          call kill(l%t)
          l=>l%next
       enddo
   l=>m_u%start
       do ii=1,m_u%n
          call kill(l%t)
          l=>l%next
       enddo
else
    if(associated(r%parent_universe)) then
       l=>r%parent_universe%start
       do ii=1,r%parent_universe%n
          call kill(l%t)
          l=>l%next
       enddo
    elseif(associated(r%t)) then
       call kill(r%t)
    endif
endif

    uni=.false.
    if(present(universe)) uni=universe

    if(uni) then
      if(r%parent_universe%nf==0) call TIE_MAD_UNIVERSE(r%parent_universe)
      nt=r%parent_universe%nf
    else
      nt=r%n
    endif





    M1=0
    M2=0
    M3=0
    MK1=0
    MK2=0
    MK3=0


    r%NTHIN=0

    nst_tot=0
    if(uni) then
     c=>r%parent_universe%start%start
    else
     C=>R%START
    endif
    do  ii=1,nt    ! WHILE(ASSOCIATED(C))

       doit=c%mag%recut.or.(ignore)

       if(present(useknob)) then
          if(useknob) then
             doit=c%magp%knob
          else
             doit=.not.c%magp%knob
          endif
       endif
       if(present(fib)) then
          doit=associated(c,fib)
       endif
       if(.not.doit) then
          NST_tot=NST_tot+C%MAG%P%nst
             if(uni) then
              c=>c%n
             else
              C=>C%NEXT
             endif
          cycle
       endif


       doit=(C%MAG%KIND==kind1.and.C%MAG%KIND==kind2.or.C%MAG%KIND==kind5.or.C%MAG%KIND==kind4)
       doit=DOIT.OR.(C%MAG%KIND==kind6.or.C%MAG%KIND==kind7)
       DOIT=DOIT.OR.(C%MAG%KIND==kind10.or.C%MAG%KIND==kind16)
       DOIT=DOIT.OR.(C%MAG%KIND==kind17.or.C%MAG%KIND==kindwiggler)


       if(doit)  then

          M1=M1+1
          NTE=1
          C%MAG%P%NST=NTE
          MK1=MK1+NTE
          C%MAG%P%METHOD=2

          r%NTHIN=r%NTHIN+1  !C%MAG%NST

          call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
          call COPY(C%MAG,C%MAGP)
       else
          if(C%MAG%KIND/=kindpa) then
             C%MAG%P%NST=1
             if(associated(C%MAG%bn))call add(C%MAG,C%MAG%P%nmul,1,0.0_dp)
             call COPY(C%MAG,C%MAGP)
          endif
       endif

       NST_tot=NST_tot+C%MAG%P%nst
             if(uni) then
              c=>c%n
             else
              C=>C%NEXT
             endif
    enddo

if(lielib_print(14)==1) then
    write(6,*) "Present of cutable Elements ",r%NTHIN
    write(6,*) "METHOD 2 ",M1,MK1
    write(6,*) "METHOD 4 ",M2,MK2
    write(6,*) "METHOD 6 ",M3,MK3
    write(6,*)   "number of Slices ", MK1+MK2+MK3
    write(6,*)   "Total NST ", NST_tot
endif


    !    CALL RING_L(R,doneit)

  END SUBROUTINE  THIN_LENS_restart


  SUBROUTINE  print_bn_an(r,n,title,filename)
    implicit none
    type(layout),target,intent(inout) ::r
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
    type(layout),target,intent(inout) ::r
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
    an=0.0_dp;bn=0.0_dp;

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


  SUBROUTINE REVERSE_BEAM_LINE(R, changeanbn )
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer J,I
    type(fibre), pointer :: P
    logical(lp), optional:: changeanbn
    logical(lp) changeanbn0

    changeanbn0=my_true
    if(present(changeanbn)) changeanbn0=changeanbn

    p=>r%start
    do i=1,r%n
       p%dir=-1
       if(changeanbn0) then
          if(associated(p%mag%an)) then
             do j=1,p%mag%p%nmul
                p%mag%bn(j)=-p%magp%bn(j)
                p%mag%an(j)=-p%magp%an(j)
                p%magp%bn(j)=-p%magp%bn(j)
                p%magp%an(j)=-p%magp%an(j)
             enddo
             if(abs(abs(p%mag%bn(1))-abs(p%mag%p%b0)).gt.1e-11_dp.or. &
                  abs(p%mag%p%b0).lt.1e-11_dp) then
                p%mag%bn(1)=-p%magp%bn(1)
                p%mag%an(1)=-p%magp%an(1)
                p%magp%bn(1)=-p%magp%bn(1)
                p%magp%an(1)=-p%magp%an(1)
             endif
             if(p%mag%p%nmul>0) call add(p,1,1,0.0_dp)
          endif
          if(associated(p%mag%volt)) p%mag%volt=-p%mag%volt
          if(associated(p%magp%volt)) p%magp%volt=-p%magp%volt
       endif
       P=>P%next
    ENDDO
  end SUBROUTINE REVERSE_BEAM_LINE

  SUBROUTINE PUTFRINGE(R, changeanbn )
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer I,perm
    type(fibre), pointer :: P
    logical(lp) changeanbn

     perm=0
    if(changeanbn) perm=1 
    p=>r%start
    do i=1,r%n
       p%mag%p%PERMFRINGE =perm
       p%magP%p%PERMFRINGE=perm
       P=>P%next
    ENDDO
  end SUBROUTINE PUTFRINGE

  SUBROUTINE PUTbend_fringe(R, changeanbn )
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    integer I
    type(fibre), pointer :: P
    logical(lp) changeanbn

    p=>r%start
    do i=1,r%n
    if(p%mag%p%b0/=0.0_dp) then
       p%mag%p%bend_fringe =changeanbn
       p%magp%p%bend_fringe =changeanbn
       write(6,*) P%mag%name, " changed to ",changeanbn
    endif
       P=>P%next
    ENDDO
  end SUBROUTINE PUTbend_fringe


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

  subroutine MESS_UP_ALIGNMENT_name(R,nom,iseed,full,sig,cut,pr)
    use gauss_dis
    IMPLICIT NONE
    TYPE(layout),target, intent(inout):: R
    integer iseed,j,ic,i
    character(nlp) nom
    type(fibre), pointer :: p
    logical(lp) f1,full,pr
    real(dp) cut,sig(6),mis(6),x,taxi(6)

    if(iseed/=0) call gaussian_seed(iseed)

    taxi=0.d0
    call context(nom)

    ic=0

    p=>r%start
    do i=1,r%n
       f1=.false.

       if(full) then
          f1=(p%mag%name ==nom )
       else
          f1=(   index(p%mag%name,nom(1:len_trim(nom)))   >0)
       endif


          if(f1) then
             ic=ic+1
             DO J=1,6
                call GRNF(X,cut)
                MIS(J)=X*SIG(J)
                taxi(j)=taxi(j)+abs(MIS(J))
             ENDDO
             call MISALIGN_FIBRE(p,mis)
          endif

       if(f1) then 
         ic=ic+1
       if(pr)   write(6,*) p%mag%name
       endif

       p=>P%next
    enddo

  if(pr) then  
   write(6,*) ic," Magnets modified "



    write(6,*) ic," Magnets misaligned "
    taxi=taxi/ic

    write(6,'(a21,3(1x,E15.8))') " <|displacements|> = ",taxi(1:3)
    write(6,'(a21,3(1x,E15.8))') " <|rotations|>     = ",taxi(4:6)
endif

  end  subroutine MESS_UP_ALIGNMENT_name

  subroutine Sigma_of_alignment(r,sig,ave)
    IMPLICIT NONE
    TYPE(layout),target, intent(inout):: R
    integer i,j,is(6)
    type(fibre), pointer :: p
    real(dp) sig(6),ave(6)
    character*23 lab(6)
    lab(1)=" Dx Average and Sigma  "
    lab(2)=" Dy Average and Sigma  "
    lab(3)=" Dz Average and Sigma  "
    lab(4)=" Dax Average and Sigma "
    lab(5)=" Day Average and Sigma "
    lab(6)=" Daz Average and Sigma "

    AVE=0.0_dp
    SIG=0.0_dp
    p=>r%start
    is=0
    do i=1,r%n
       do j=1,3
          if(p%chart%D_IN(j)/=0.0_dp) then
             is(j)=is(j)+1
             ave(j)=p%chart%D_IN(j)+ave(j)
             sig(j)=p%chart%D_IN(j)**2+sig(j)
          endif
       enddo
       do j=4,6
          if(p%chart%ANG_IN(j)/=0.0_dp) then
             is(j)=is(j)+1
             ave(j)=p%chart%ANG_IN(j)+ave(j)
             sig(j)=p%chart%ANG_IN(j)**2+sig(j)
          endif
       enddo
       p=>p%next
    enddo

    do i=1,6
       if(is(i)/=0) then
          ave(i)=ave(i)/is(i)
          sig(i)=sig(i)/is(i)
          sig(i)=sqrt(sig(i)-ave(i)**2)
          write(6,*) is(i), " Magnets misaligned "
          write(6,"(1x,a23,2(1x,E15.8))") lab(i),ave(i),sig(i)
       endif
    enddo

  end subroutine Sigma_of_alignment




  SUBROUTINE dyn_aper(L,x_in,n_in,ang_in,ang_out,del_in,dlam,pos,nturn,ite,state,mf)
    IMPLICIT NONE
    type(layout),target, intent(inout) :: L
    real(dp) x(6)
    REAL(DP) x_in,del_in,closed(6),r(6),rt(6)
    REAL(DP) lamT,lams,lamu,dlam,DLAMT,DX,ang,ang_in,ang_out
    integer pos,nturn,i,st,ite,ic,mf,n_in,j_in
    TYPE(INTERNAL_STATE) STATE
    !
    !    TYPE(REAL_8) Y(6)
    !    TYPE(DAMAP) ID
    !    TYPE(NORMALFORM) NORM

    closed=0.0_dp
    !    STATE=STATE+NOCAVITY0
    if(state%nocavity) closed(5)=del_in

    CALL FIND_ORBIT(L,CLOSED,pos,STATE,1e-5_dp)
    write(6,*) "closed orbit "
    write(6,*) CLOSED
    write(mf,201) closed
    ang= (ang_out-ang_in)/n_in
    lamt=1.0_dp
    do j_in=0,n_in

       x=0.0_dp
       x(1)=x_in*cos(j_in*ang+ang_in)
       x(3)=x_in*sin(j_in*ang+ang_in)
       x(5)=del_in


       dx=0.3_dp

       r=0.0_dp;rt=0.0_dp;
       lams=0.0_dp
       lamu=0.0_dp

       DLAMT=DX

       !    lamt=ONE
       ic=0
       do while(DLAMT>dlam.and.ic<ite)

          ic=ic+1
          R=0.0_dp;
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
             lamt=(lams+lamt)/2.0_dp
          else
             lams=lamt
             IF(LAMU<DX) THEN
                lamt=DX+lamt
             ELSE
                lamt=(lamu+lamt)/2.0_dp
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



  SUBROUTINE dyn_aperalex(L,x_in,del_in,dx,dlam,pos,nturn,ite,state,mf,targ_tune,fixp)
    IMPLICIT NONE
    type(layout),target, intent(inout) :: L
    real(dp) x(6)
    REAL(DP) x_in,del_in,closed(6),r(6),rt(6),targ_tune(2)
    REAL(DP) lamT,lams,lamu,dlam,DLAMT,DX,ang,ang_in,ang_out
    integer pos,nturn,i,st,ite,ic,mf,n_in,j_in
    TYPE(INTERNAL_STATE) STATE
    logical(lp) fixp
    !
    !    TYPE(REAL_8) Y(6)
    !    TYPE(DAMAP) ID
    !    TYPE(NORMALFORM) NORM

    closed=0.0_dp
    !    STATE=STATE+NOCAVITY0

    !    if(state%nocavity)
    closed(5)=del_in
    if(fixp) then
       CALL FIND_ORBIT(L,CLOSED,pos,STATE,1e-5_dp)
       write(6,*) "closed orbit "
       write(6,*) CLOSED
    else
       closed(1:4)=0.0_dp
       closed(6)=0.0_dp
    endif
    !    write(mf,201) closed
    n_in=1
    ang_in=pi/4.0_dp
    ang_out=pi/4.0_dp
    ang= (ang_out-ang_in)/n_in
    lamt=1.0_dp
    j_in=0
    !  do j_in=0,n_in

    x=0.0_dp
    x(1)=x_in*cos(j_in*ang+ang_in)
    x(3)=x_in*sin(j_in*ang+ang_in)
    !       x(5)=del_in



    r=0.0_dp;rt=0.0_dp;
    lams=0.0_dp
    lamu=0.0_dp

    DLAMT=DX

    !    lamt=ONE
    ic=0
    do while(DLAMT>dlam.and.ic<ite)

       ic=ic+1
       R=0.0_dp;
       r(1:4)=lamt*x(1:4)
       if(state%nocavity) then
          if(fixp) then
             rt=r+closed
          else
             rt(1:4)=r(1:4)
             rt(6)=0.0_dp
             rt(5)=del_in
          endif
       else
          if(fixp) then
             rt=r+closed
             rt(5)=rt(5)+del_in
          else
             rt=r
             rt(5)=rt(5)+del_in
          endif
       endif

       !   write(6,*) rt
       do i=1,nturn
          st=track_flag(L,rt,pos,state)
          if(st/=0) exit
       enddo
       !   write(6,*) i,check_stable
       !   write(6,*) rt
       !   pause
       if(st/=0) then
          lamu=lamt
          lamt=(lams+lamt)/2.0_dp
          !  write(mf,*) "unstable ",lamt

       else
          lams=lamt
          IF(LAMU<DX) THEN
             lamt=DX+lamt
          ELSE
             lamt=(lamu+lamt)/2.0_dp
          ENDIF
          !  write(mf,*) "stable ",lamt

       endif
       DLAMT=sqrt(x(1)**2+x(3)**2)*ABS(LAMU-LAMS)
       !              write(mf,*) "dlamt ",dlamt,sqrt(x(1)**2+x(3)**2)

    enddo
    write(6,*) ic,lamS*x(1)   !,lamS*x(3),(j_in*ang+ang_in)/twopi

    write(mf,202) targ_tune,lamS*x(1)      !,lamS*x(3),lamS*x(1)+closed(1),lamS*x(3)+closed(3),DLAMT
    lamt=lamt*0.8_dp
    !    enddo
201 FORMAT(6(1X,D18.11))
202 FORMAT(3(1X,D18.11))

  end SUBROUTINE dyn_aperalex




end module S_fitting
