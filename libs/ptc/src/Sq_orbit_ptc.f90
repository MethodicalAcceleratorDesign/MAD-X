module orbit_ptc
  use madx_keywords
  !  use pointer_lattice
  !  use madx_ptc_module
  implicit none
  public
  ! TYPE(INTERNAL_STATE),POINTER :: MY_ORBIT_STATE

  PRIVATE ORBIT_TRACK_NODEP,ORBIT_TRACK_NODE_Standard_R

  REAL(dp)  X_ORBIT(6)
  REAL(DP) :: XBIG = 1.D10
  type(work)  :: w1_orbit,w2_orbit
  type(fibre), pointer :: p_orbit
  integer :: ptc_node_old=-1
  integer :: mdebug = 0
  integer :: mbugplot = 0
  logical :: include_patch =.false.
  logical :: fill_patch =.false.
  integer :: n_patch=0
  integer :: n_fill_patch=0
  integer :: n_used_patch=0, extra_node=0  
  real(dp) :: t0_main=0.0_dp
  character(nlp), allocatable :: orbitname(:)
  !   integer mfff
  INTERFACE ORBIT_TRACK_NODE
     !LINKED
     MODULE PROCEDURE ORBIT_TRACK_NODE_Standard_R
     ! for speed
     !          MODULE PROCEDURE ORBIT_TRACK_NODER
     MODULE PROCEDURE ORBIT_TRACK_NODEP
  END INTERFACE

  INTERFACE ORBIT_TRACK_NODE_Standard
     !LINKED
     ! MODULE PROCEDURE ORBIT_TRACK_NODER
     ! MODULE PROCEDURE ORBIT_TRACK_NODEP
     MODULE PROCEDURE ORBIT_TRACK_NODE_Standard_R
     MODULE PROCEDURE ORBIT_TRACK_NODEP
  END INTERFACE

  INTERFACE ORBIT_MAKE_NODE_LAYOUT
     !LINKED
     MODULE PROCEDURE ORBIT_MAKE_NODE_LAYOUT_accel
  END INTERFACE

  TYPE(ORBIT_LATTICE), pointer :: my_ORBIT_LATTICE

contains



  SUBROUTINE PUT_RAY(X1,X2,X3,X4,X5,X6)
    IMPLICIT NONE
    REAL(DP), INTENT(IN):: X1,X2,X3,X4,X5,X6
    X_ORBIT(1)=X1
    X_ORBIT(2)=X2
    X_ORBIT(3)=X3
    X_ORBIT(4)=X4
    X_ORBIT(5)=X5
    X_ORBIT(6)=X6
  END SUBROUTINE PUT_RAY

  SUBROUTINE PUT_RAY_track(node,X1,X2,X3,X4,X5,X6)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT):: X1,X2,X3,X4,X5,X6
    integer, intent(in) :: node
    real(dp) X(6)
    X(1)=X1
    X(2)=X2
    X(3)=X3
    X(4)=X4
    X(5)=X5
    X(6)=X6
    CALL ORBIT_TRACK_NODE(node,X)
    X1=X(1)
    X2=X(2)
    X3=X(3)
    X4=X(4)
    X5=X(5)
    X6=X(6)

  END SUBROUTINE PUT_RAY_track



  SUBROUTINE GET_RAY(X1,X2,X3,X4,X5,X6)
    IMPLICIT NONE
    REAL(DP), INTENT(OUT):: X1,X2,X3,X4,X5,X6
    X1=X_ORBIT(1)
    X2=X_ORBIT(2)
    X3=X_ORBIT(3)
    X4=X_ORBIT(4)
    X5=X_ORBIT(5)
    X6=X_ORBIT(6)
  END SUBROUTINE GET_RAY

  !===========================================================
  ! This subroutine should be called before particle tracking.
  !  It specifies the type of the task that will be performed
  !  in ORBIT before particle tracking for particular node.
  !  i_task = 0 - do not do anything
  !  i_task = 1 - energy of the sync. particle changed
  !===========================================================
  SUBROUTINE GET_task(i_node,i_task)
    IMPLICIT NONE
    INTEGER  i_node,i_task

    i_task    =my_ORBIT_LATTICE%ORBIT_NODES(i_node)%entering_task

  END SUBROUTINE GET_task

  SUBROUTINE GET_info(i_node,dl,betax,betay,alphax,alphay,etax,etaxp)
    IMPLICIT NONE
    INTEGER  i_node
    real(dp) dl,betax,betay,alphax,alphay,etax,etaxp

    dl    =my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(1)
    betax =my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(2)
    betay =my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(3)
    alphax=my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(4)
    alphay=my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(5)
    etax  =my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(6)
    etaxp =my_ORBIT_LATTICE%ORBIT_NODES(i_node)%lattice(7)

  END SUBROUTINE GET_info


  SUBROUTINE GET_deltae(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_deltae
  END SUBROUTINE GET_deltae

  !  SUBROUTINE GET_dppfac(X)
  !    IMPLICIT NONE
  !    REAL(DP) X
  !    X=my_ORBIT_LATTICE%ORBIT_dppfac
  !  END SUBROUTINE GET_dppfac

  SUBROUTINE GET_gamma(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_energy
  END SUBROUTINE GET_gamma

  SUBROUTINE GET_total_energy(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_energy
  END SUBROUTINE GET_total_energy

  SUBROUTINE GET_brho(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_brho
  END SUBROUTINE GET_brho

  SUBROUTINE GET_BETA0(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_BETA0
  END SUBROUTINE GET_BETA0

  SUBROUTINE GET_omega(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_OMEGA
  END SUBROUTINE GET_omega

  SUBROUTINE GET_HARMONIC(X)
    IMPLICIT NONE
    INTEGER  X
    X=NINT(my_ORBIT_LATTICE%orbit_HARMONIC)
  END SUBROUTINE GET_HARMONIC

  SUBROUTINE GET_kinetic(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%orbit_kinetic
  END SUBROUTINE GET_kinetic

  SUBROUTINE GET_CHARGE(X)
    IMPLICIT NONE
    integer X
 !   real(dp) X
!    X=my_ORBIT_LATTICE%ORBIT_CHARGE
   X=nint(my_ORBIT_LATTICE%ORBIT_CHARGE)
  END SUBROUTINE GET_CHARGE

  SUBROUTINE GET_MASS_AMU(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_mass_in_amu
  END SUBROUTINE GET_MASS_AMU

  SUBROUTINE GET_GAMMAT(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_gammat
  END SUBROUTINE GET_GAMMAT

  SUBROUTINE GET_CIRCUMFERENCE(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_L
  END SUBROUTINE GET_CIRCUMFERENCE

  SUBROUTINE GET_LMAX(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_LMAX
  END SUBROUTINE GET_LMAX

  SUBROUTINE GET_P0C(X)
    IMPLICIT NONE
    REAL(DP) X
    X=my_ORBIT_LATTICE%ORBIT_P0C
  END SUBROUTINE GET_P0C

  SUBROUTINE GET_N_NODE(X)
    IMPLICIT NONE
    INTEGER X
    X=my_ORBIT_LATTICE%ORBIT_N_NODE
  END SUBROUTINE GET_N_NODE

  SUBROUTINE TRACK_ONE_NODE(K)
    IMPLICIT NONE
    INTEGER K
    CALL ORBIT_TRACK_NODE(K,X_ORBIT)
  END SUBROUTINE TRACK_ONE_NODE

  SUBROUTINE TRACK_ONE_TURN(K)
    IMPLICIT NONE
    INTEGER K
    CALL ORBIT_TRACK_ONE_TURN(K,X_ORBIT)
  END SUBROUTINE TRACK_ONE_TURN


  !  TRACKING AND SETTING UP
  SUBROUTINE ORBIT_TRACK_ONE_TURN(K,X)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    INTEGER K,I
    DO I=K,my_ORBIT_LATTICE%ORBIT_N_NODE
       CALL ORBIT_TRACK_NODE(I,X)
    ENDDO
    DO I=1,K-1
       CALL ORBIT_TRACK_NODE(I,X)
    ENDDO
  END SUBROUTINE ORBIT_TRACK_ONE_TURN


  SUBROUTINE Locate_orbit_start(n,m)
    IMPLICIT NONE
    INTEGER I,n,m,k,k1
    TYPE(INTEGRATION_NODE), POINTER  :: T


    k1=0

    do k=1,my_ORBIT_LATTICE%ORBIT_N_NODE

       T=>my_ORBIT_LATTICE%ORBIT_NODES(K)%NODE
       DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(K)%dpos
          if(associated(t,t%parent_fibre%t1)) k1=k1+1
          if(k1==n) then
             m=k
             goto 100
          endif
          T=>T%NEXT
       ENDDO

    enddo

100 continue
    T=>my_ORBIT_LATTICE%ORBIT_NODES(m)%NODE
    write(6,*) " Fibre position ",n,t%parent_fibre%mag%name
    write(6,*) "  position in fibre ",t%pos_in_fibre,t%parent_fibre%mag%p%nst
    write(6,*) " Orbit node ",m

  end SUBROUTINE Locate_orbit_start





  SUBROUTINE TRACK_NODE_fake_totalpath_half(T,X,STATE,after)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), target :: STATE
    logical(lp) after
    real(dp) beta0
    type(work) w
    integer nst,PATCHT

    beta0=t%PARENT_FIBRE%beta0
    nst=t%pos_in_fibre-2

    if(t%PARENT_FIBRE%mag%kind==kind4.and.accelerate) then

       if(after.and.t%cas==case0.and.associated(t%parent_fibre%mag%c4%acc)) then
          w=t%parent_fibre
          x(5)= x(5) + t%parent_fibre%mag%c4%acc%de(nst)/w%p0c
       endif !  kind4

    endif

    !  if(state%time) then
    x(6)=x(6)+T%ds_AC*0.5_dp*(1.0_dp/beta0+x(5))/root(1.0_dp+2*x(5)/beta0+x(5)**2)
    ! else
    !    x(6)=x(6)+T%ds_AC*half

    ! endif


    !    if(T%CAS==CASEP1.and.include_patch) then
    !       PATCHT=t%parent_fibre%patch%time
    !       IF(PATCHT/=0.AND.PATCHT/=2) THEN
    !         if(state%time) then
    !          X(6)=X(6)+t%parent_fibre%PATCH%a_T*half/t%parent_fibre%beta0
    !         else
    !          X(6)=X(6)+t%parent_fibre%PATCH%a_T*half
    !         endif
    !       ENDIF
    !    endif

    !    if(T%CAS==CASEP2.and.include_patch) then
    !       PATCHT=t%parent_fibre%patch%time
    !       IF(PATCHT/=0.AND.PATCHT/=1) THEN
    !         if(state%time) then
    !          X(6)=X(6)+t%parent_fibre%PATCH%b_T*half/t%parent_fibre%beta0
    !         else
    !          X(6)=X(6)+t%parent_fibre%PATCH%b_T*half
    !         endif
    !       ENDIF
    !    endif

  end SUBROUTINE TRACK_NODE_fake_totalpath_half

  SUBROUTINE TRACK_NODE_fake_totalpath_half_plain(T,X,STATE)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), target :: STATE
    integer PATCHT
    !    if(state%time) then
    x(6)=x(6)+T%ds_AC*0.5_dp/t%PARENT_FIBRE%beta0
    !    else
    !       x(6)=x(6)+T%ds_AC*half
    !    endif

    !    if(T%CAS==CASEP1.and.include_patch) then
    !       PATCHT=t%parent_fibre%patch%time
    !       IF(PATCHT/=0.AND.PATCHT/=2) THEN
    !         if(state%time) then
    !          X(6)=X(6)+t%parent_fibre%PATCH%a_T*half/t%parent_fibre%beta0
    !         else
    !          X(6)=X(6)+t%parent_fibre%PATCH%a_T*half
    !         endif
    !       ENDIF
    !    endif

    !    if(T%CAS==CASEP2.and.include_patch) then
    !       PATCHT=t%parent_fibre%patch%time
    !       IF(PATCHT/=0.AND.PATCHT/=1) THEN
    !         if(state%time) then
    !          X(6)=X(6)+t%parent_fibre%PATCH%b_T*half/t%parent_fibre%beta0
    !         else
    !          X(6)=X(6)+t%parent_fibre%PATCH%b_T*half
    !         endif
    !       ENDIF
    !    endif

  end SUBROUTINE TRACK_NODE_fake_totalpath_half_plain

  subroutine get_ideal_harmonic(R,xharm,dt,state)
    IMPLICIT NONE
    REAL(DP)  X(6),xharm,dt,om,om0
    TYPE(LAYOUT), target  :: R
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), target :: STATE
    integer i
    logical found

    found=.false.
    if(.not.associated(R%T))then
       call make_node_layout(r)
    endif
    x=0.0_dp
    t=>R%t%start
    om=1.e38_dp
    do i=1,R%T%N

       if(t%parent_fibre%mag%kind==kind4) then
          om0=t%parent_fibre%mag%freq
          if(om0/=0.0_dp) then
             if(om0<om) om=om0
             found=.true.
          endif
       endif

       call TRACK_NODE_fake_totalpath_half_plain(T,X,STATE)
       call TRACK_NODE_fake_totalpath_half_plain(T,X,STATE)
       t=>t%next
    enddo

    if(found) then
       xharm=om*x(6)/clight
       dt=clight/om
    else
       write(6,*) "No cavities found"
       xharm=0.0_dp
       dt=0.0_dp
    endif

  end subroutine get_ideal_harmonic

  SUBROUTINE ORBIT_TRACK_NODE_Standard_R(K,X,STATE)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    INTEGER K,I
    LOGICAL(LP) U,cav
    REAL(DP) X5,dt,dt_orbit_sync
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), target, OPTIONAL :: STATE
    TYPE(INTERNAL_STATE), pointer :: STATE0

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       x(1:4)=x(1:4)*1.e-3_dp
       X5=X(5)
       X(5)=X(6)/my_ORBIT_LATTICE%ORBIT_P0C
       X(6)=X5/my_ORBIT_LATTICE%ORBIT_OMEGA
    ENDIF

    u=my_false

    T=>my_ORBIT_LATTICE%ORBIT_NODES(K)%NODE

    if(present(state)) then
       state0=>state
    else
       state0=>my_ORBIT_LATTICE%STATE
    endif


    IF(accelerate.OR.RAMP) THEN !accelerate
       if(K/=ptc_node_old) then !accelerate
          ptc_node_old=k !accelerate
          first_particle=.true.
          if(k==1) n_used_patch=n_used_patch+1
       endif !accelerate
    ENDIF !accelerate




    DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(K)%dpos



       if(first_particle.and.accelerate) then !accelerate
          if(t%previous%parent_fibre%mag%kind==kind4.and.associated(t%previous%parent_fibre%mag%c4%acc)) then
             !  if(t%previous%cas==case0) t%previous%parent_fibre=t%previous%parent_fibre%mag%c4%acc%w2
             if(t%previous%cas==case0) then
                call ORBIT_up_grade_mag(t%previous)
             endif
          endif
          if(t%parent_fibre%mag%kind==kind4.and.accelerate) then !accelerate
             cav=(t%cas==case0)
             if(associated(t%parent_fibre%mag%c4%acc).and.cav) then !accelerate
                if(t%parent_fibre%mag%c4%acc%pos==1.and.t%pos_in_fibre==3) call find_all_energies(t,state0) !accelerate
             endif !accelerate
             if(associated(t%parent_fibre%mag%c4%acc).and.cav) call set_cavity(t,state0,dt_orbit_sync)
             !          if(cav) call set_cavity(t,state0,dt_orbit_sync)
             p_orbit=>t%parent_fibre   !p_orbit%mag%c4%t
          endif !accelerate



       endif !accelerate

       if(first_particle.and.(accelerate.or.ramp)) then !accelerate
          call TRACK_NODE_fake_totalpath_half(T,x_orbit_sync,STATE0,my_true)  ! accelerate
       endif

       if(RAMP.and.first_particle) then !modulate
          if(t%parent_fibre%mag%slow_ac.and.t%cas==CASEP1) CALL do_ramping_r(T,x_orbit_sync(6),STATE0) !modulate
       endif !modulate

       if(u) exit
       CALL TRACK_NODE_SINGLE(T,X,STATE0) !,my_ORBIT_LATTICE%ORBIT_CHARGE

       if(.not.CHECK_STABLE) then
          CALL RESET_APERTURE_FLAG
          u=my_true
! EB - comment by Elena Benedetto - Only the 4 transverse coordinates are multiplied by xbig, 
! to avoid loosing the information on the phase (when the mod 2pi is taken)
         ! x(1)=XBIG
          x(1:4) =x(1:4)*XBIG
          if(wherelost==1) then
             t%lost=t%lost+1
          endif
          exit
       endif



       if(first_particle.and.(accelerate.or.ramp)) then !accelerate
          call TRACK_NODE_fake_totalpath_half(T,x_orbit_sync,STATE0,my_true)  ! accelerate
       endif


       if(t%parent_fibre%mag%kind==kind4.and.accelerate.and.associated(t%parent_fibre%mag%c4%acc)) then
          !    if(t%parent_fibre%mag%kind==kind4.and.accelerate) then
          if(t%cas==case0) then
             x(6)=x(6)-dt_orbit_sync
             call ORBIT_up_grade_x(x,t)
             if(first_particle) call ORBIT_up_grade_x(x_orbit_sync,t)
          endif
       endif

       if(mbugplot/=0) then
          write(mbugplot,*) t%pos,t%cas,t%parent_fibre%mag%name
          write(mbugplot,'(4(1X,D18.11))') x(1:2),x(5:6)  !, dt_orbit_sync
       endif
       if(fill_patch.and.associated(t,my_ORBIT_LATTICE%tp)) then
          n_fill_patch=n_fill_patch+1
          my_ORBIT_LATTICE%dt(n_fill_patch)=x(6)
          x(6)=0.d0
       endif

       if(associated(t,my_ORBIT_LATTICE%tp).and.n_used_patch<=n_patch.and.include_patch) then
          if(state0%time) then
             X(6)=X(6)-my_ORBIT_LATTICE%dt(n_used_patch)/t%parent_fibre%beta0
             !  write(30,*)n_used_patch, my_ORBIT_LATTICE%dt(n_used_patch),t%parent_fibre%beta0
          else
             X(6)=X(6)-my_ORBIT_LATTICE%dt(n_used_patch)
          endif
       endif

       T=>T%NEXT
    ENDDO
    first_particle=.false.

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       x(1:4)=x(1:4)*1.e3_dp
       X5=X(5)
       X(5)=X(6)*my_ORBIT_LATTICE%ORBIT_OMEGA
       X(6)=X5*my_ORBIT_LATTICE%ORBIT_P0C
    ENDIF

  end SUBROUTINE ORBIT_TRACK_NODE_Standard_R


  subroutine energize_ORBIT_lattice(t)
    implicit none
    real(dp), optional :: t
    type(fibre), pointer :: p
    type(work) werk,travail
    real(dp) e_in,t0,freqs
    integer i
    logical found

    found=.false.
    freqs=1.e38_dp
    werk=0
    travail=0
    t0=x_orbit_sync(6) ! read from the initial_settings file and multiplied by clight
    if(present(t)) t0=t
    if(fill_patch) then
       write(6,*) " filling patches with t= x0 from main program "
       t0=t0_main
    endif
    write(6,*) "energize at time ", t0,t0/clight
    write(6,*) "Initial Frequency of First Cavity", paccfirst%mag%c4%freq

    call find_acc_energy(paccfirst,t0,e_in,my_true) ! new
    !      call find_acc_energy(paccfirst,x_orbit_sync(6),e_in,.false.)
    call find_energy(werk,kinetic=e_in)



    p=>paccfirst%parent_layout%start
    do i=1,paccfirst%parent_layout%n
       travail=p
       if(p%mag%kind==kind4) then
          p%mag%c4%freq  = p%mag%c4%freq*werk%beta0/travail%beta0
          p%magp%c4%freq = p%mag%c4%freq
          if(p%mag%c4%freq<freqs)then
             freqs=p%mag%c4%freq
             found=.true.
          endif
          call find_acc_energy(paccfirst,t0,e_in,my_true)  ! new
       endif
       p=werk
       p=>p%next
    enddo
    if(freqs/=0.and.found) then
       my_ORBIT_LATTICE%ORBIT_OMEGA=twopi*freqs/CLIGHT

! EB  beginning - added by Elena Benedetto - to have correct energy already at the beginning
       my_ORBIT_LATTICE%ORBIT_P0C=werk%P0C
       my_ORBIT_LATTICE%ORBIT_kinetic=werk%kinetic
       my_ORBIT_LATTICE%orbit_beta0=werk%beta0
       my_ORBIT_LATTICE%orbit_brho=werk%brho
       my_ORBIT_LATTICE%orbit_energy=werk%energy
       my_ORBIT_LATTICE%orbit_gamma=1.0_dp/werk%gamma0i
! EB end - added by Elena Benedetto

    else
       write(6,*) " cavity with frequency problems ", freqs,found
       stop
    endif
    write(6,*) "Final Frequency of First Cavity", paccfirst%mag%c4%freq
    write(6,*) "Initial and Final beta0 ",travail%beta0,werk%beta0
    write(6,*) "Starting time of simulations =",t0/clight ," and kinetic energy =",my_ORBIT_LATTICE%ORBIT_kinetic

  end subroutine energize_ORBIT_lattice

  subroutine ORBIT_up_grade_mag(t)
    implicit none
    type(fibre), pointer :: p
    integer i
    TYPE(INTEGRATION_NODE), POINTER  :: T
    type(acceleration), pointer :: a
    real(dp) freqs
    logical found

    found=.false.
    freqs=1.e38_dp

    a=>t%parent_fibre%mag%c4%acc


    p=>t%parent_fibre

    do i=1,p%parent_layout%n  !-1
       if(p%mag%kind==kind4) then
          p%mag%c4%freq  = p%mag%c4%freq*a%w2%beta0/a%w1%beta0
          p%magp%c4%freq = p%mag%c4%freq
          if(p%mag%c4%freq<freqs.and.associated(p%mag%c4%acc))then
             freqs=p%mag%c4%freq
             found=.true.
          endif
       endif
       p=a%w2
       p=>p%next
    enddo

    if(freqs/=0.and.found) then
      my_ORBIT_LATTICE%ORBIT_OMEGA=twopi*freqs/CLIGHT

! EB beginning - added by Elena Benedetto - to pass the correct energy from PTC RF table to Orbit
       my_ORBIT_LATTICE%ORBIT_P0C=a%w2%P0C
       my_ORBIT_LATTICE%ORBIT_kinetic=a%w2%kinetic
       my_ORBIT_LATTICE%orbit_beta0=a%w2%beta0
       my_ORBIT_LATTICE%orbit_brho=a%w2%brho
       my_ORBIT_LATTICE%orbit_energy=a%w2%energy
       my_ORBIT_LATTICE%orbit_gamma=1.0_dp/a%w2%gamma0i       
! EB end - added by Elena Benedetto

    else
       write(6,*) " ORBIT_up_grade_mag ", freqs,found
       write(6,*) " cavity with frequency problems ", freqs,found
       stop
    endif

  end subroutine ORBIT_up_grade_mag

  subroutine ORBIT_up_grade_x(x,t)
    implicit none
    real(dp),intent(inout) :: x(6)
    TYPE(INTEGRATION_NODE), POINTER  :: T
    type(acceleration), pointer :: a

    a=>t%parent_fibre%mag%c4%acc

    X(2)=X(2)*a%w1%P0C/a%w2%P0C
    X(4)=X(4)*a%w1%P0C/a%w2%P0C
    X(5)=root(1.0_dp+2.0_dp*X(5)/a%w1%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
    X(5)=X(5)*a%w1%P0C/a%w2%P0C-1.0_dp !X(5) = DP/P0C_NEW
    X(5)=(2.0_dp*X(5)+X(5)**2)/(root(1.0_dp/a%w2%BETA0**2+2.0_dp*X(5)+X(5)**2) &
         +1.0_dp/a%w2%BETA0)
    !if(enforce_zero_x5) x(5)=zero
    !    x(6)=x(6)*my_ORBIT_LATTICE%orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA

    !write(6,*) orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA
  END subroutine ORBIT_up_grade_x

  SUBROUTINE find_all_energies(t,state0)
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE), POINTER  :: T,c
    TYPE(internal_state), target  :: state0
    integer i,ipause,mypause
    real(dp)e_fin,beta0,e_in,de,e_in0,xsync(6),t_fin
    !    type(probe) xs
    type(acceleration), pointer :: a
    type(work) w

    call find_one_turn_final_energy(t,state0,e_fin)
    !    xs%x=zero
    !    xs%ac%om=zero
    !    xs%ac%t=xsm0%ac%t
    w=t%parent_fibre%previous

    e_in0=w%kinetic
    e_in=e_in0

    xsync=x_orbit_sync

    i=0
    c=>t

    do while(.true.)

       beta0=C%PARENT_FIBRE%beta0
       call find_energy(w,kinetic=e_in)
       C%PARENT_FIBRE%beta0=w%beta0
       call TRACK_NODE_fake_totalpath_half_plain(c,xsync,STATE0)
       C%PARENT_FIBRE%beta0=beta0

       if(c%parent_fibre%mag%kind==kind4.and.c%cas==case0) then  !.and.c%pos_in_fibre==3) then

          a=>c%parent_fibre%mag%c4%acc

          de=(e_fin-e_in0)*a%r/a%nst

          a%de(c%pos_in_fibre-2)=de
          a%e_in(c%pos_in_fibre-2)=e_in

          a=>c%parent_fibre%magp%c4%acc
          a%de(c%pos_in_fibre-2)=de
          a%e_in(c%pos_in_fibre-2)=e_in



          e_in=e_in+de
       endif


       if(c%parent_fibre%mag%kind==kind4) then
          a=>c%parent_fibre%mag%c4%acc
          if(associated(c,c%parent_fibre%t2)) then
             if(associated(c%parent_fibre%mag%c4%acc%next%t1%next%next,t)) then

                e_in=e_in-de      ! remove last

                t_fin=xsync(6)   ! to insure ridiculous self-consistancy on bare ideal no patch lattice
                !        call find_acc_energy(c%parent_fibre,t_fin,e_fin,.false.)  ! new
                call find_acc_energy(c%parent_fibre,t_fin,e_fin,my_true)  ! to insure ridiculous self-consistancy

                de=e_fin-e_in

                a=>c%parent_fibre%mag%c4%acc
                a%de(c%pos_in_fibre-4)=de
                a%e_in(c%pos_in_fibre-4)=e_in
                a=>c%parent_fibre%magp%c4%acc
                a%de(c%pos_in_fibre-4)=de
                a%e_in(c%pos_in_fibre-4)=e_in
                !           write(6,*) i,de,e_in


                exit
             endif
          endif
       endif

       beta0=C%PARENT_FIBRE%beta0
       call find_energy(w,kinetic=e_in)
       C%PARENT_FIBRE%beta0=w%beta0
       call TRACK_NODE_fake_totalpath_half_plain(c,xsync,STATE0)
       C%PARENT_FIBRE%beta0=beta0

       c=>c%next


       i=i+1
       if(i>c%parent_fibre%parent_layout%t%n) then
          write(6,*) " error acceleration loop find_final_energy "
          ipause=mypause( 1)
          stop 1950
       endif


    enddo
  end SUBROUTINE find_all_energies

  SUBROUTINE find_all_tc_for_restarting(r,tc,nc)
    IMPLICIT NONE
    TYPE(layout), POINTER  :: r
    TYPE(fibre), POINTER  :: f
    real(dp), allocatable :: tc(:)
    integer i,n,nc



    f=>r%start
    n=0
    do i=1,r%n

       if(f%mag%kind==kind4) then
          n=n+1
       endif
       f=>f%next
    enddo
    nc=n
    allocate(tc(nc))
    n=0
    do i=1,r%n

       if(f%mag%kind==kind4) then
          n=n+1
          tc(n)=f%mag%c4%t
       endif
       f=>f%next
    enddo

  end SUBROUTINE find_all_tc_for_restarting

  SUBROUTINE set_all_tc_for_restarting(r,tc,nc)
    IMPLICIT NONE
    TYPE(layout), POINTER  :: r
    TYPE(fibre), POINTER  :: f
    real(dp), allocatable :: tc(:)
    integer i,n,nc,ne



    f=>r%start

    n=0
    ne=0
    do i=1,r%n

       if(f%mag%kind==kind4) then
          if(associated(f%mag%c4%acc)) then
             n=n+1
             f%mag%c4%t=tc(n)
          else
             ne=ne+1
          endif
       endif
       f=>f%next
    enddo
    write(6,*)ne," cavities have no table "

  end SUBROUTINE set_all_tc_for_restarting

  SUBROUTINE find_one_turn_final_energy(t,state0,e_fin)
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE), POINTER  :: T,c
    TYPE(internal_state), target  :: state0
    integer i,ipause,mypause
    real(dp)e_fin,t_fin,beta0,e_in,xsync(6)
    !    type(probe) xs
    type(acceleration), pointer :: a
    type(work) w

    !    xs%x=zero
    !    xs%ac%om=zero
    !    xs%ac%t=xsm0%ac%t
    w=t%parent_fibre%previous   !  correspounds to last fibre's energy

    xsync=x_orbit_sync

    e_in=w%kinetic

    i=0
    c=>t

    do while(.true.)

       beta0=C%PARENT_FIBRE%beta0
       call find_energy(w,kinetic=e_in)
       C%PARENT_FIBRE%beta0=w%beta0
       call TRACK_NODE_fake_totalpath_half_plain(c,xsync,STATE0)
       C%PARENT_FIBRE%beta0=beta0

       if(c%parent_fibre%mag%kind==kind4.and.c%cas==case0) then
          t_fin=xsync(6)
          !       t_fin=xs%ac%t
          call find_acc_energy(c%parent_fibre,t_fin,e_fin,my_false)
          e_in=e_fin
       endif

       beta0=C%PARENT_FIBRE%beta0
       call find_energy(w,kinetic=e_in)
       C%PARENT_FIBRE%beta0=w%beta0
       call TRACK_NODE_fake_totalpath_half_plain(c,xsync,STATE0)
       C%PARENT_FIBRE%beta0=beta0

       if(c%parent_fibre%mag%kind==kind4) then
          a=>c%parent_fibre%mag%c4%acc
          if(associated(c,c%parent_fibre%t2)) then
             if(associated(c%parent_fibre%mag%c4%acc%next%t1%next%next,t)) then
                exit
             endif
          endif
       endif
       c=>c%next


       i=i+1
       if(i>c%parent_fibre%parent_layout%t%n) then
          write(6,*) " error acceleration loop find_final_energy "
          ipause=mypause( 1)
          stop 1951
       endif


    enddo

    !  write(6,*) " one turn cav x(6) ",xsync(6)
  end SUBROUTINE find_one_turn_final_energy



  SUBROUTINE find_acc_energy(p,t_fin,e_fin,insert)
    IMPLICIT NONE
    TYPE(fibre), POINTER  :: p
    integer i,it
    real(dp)  dtot,e_fin ,ti,t_fin,rat
    type(work) w
    type(acceleration), pointer :: a
    type(element), pointer :: el
    type(elementp), pointer :: elp
    logical(lp) insert

    el=>p%mag
    elp=>p%magp
    a=>el%c4%acc


    ti=t_fin/clight    ! time in milliseconds
    if(ti>a%tableau(a%n)%temps.or.ti<a%tableau(1)%temps) then
       ! write(6,*) " Acceleration must stop time =",ti
       if(ti>a%tableau(a%n)%temps) then
          if(insert) then
             do i=1,el%c4%nf
                el%c4%f(i)= a%tableau(a%n)%volt(i)
                el%c4%ph(i)= a%tableau(a%n)%phase(i)
                elp%c4%f(i)=a%tableau(a%n)%volt(i)
                elp%c4%ph(i)= a%tableau(a%n)%phase(i)
             enddo
          endif
          e_fin=a%tableau(a%n)%energie
       else
          if(insert) then
             do i=1,el%c4%nf
                el%c4%f(i)= a%tableau(1)%volt(i)
                el%c4%ph(i)= a%tableau(1)%phase(i)
                elp%c4%f(i)=a%tableau(1)%volt(i)
                elp%c4%ph(i)= a%tableau(1)%phase(i)
             enddo
          endif
          e_fin=a%tableau(1)%energie
       endif
    else  ! time inside table


       dtot=(a%tableau(a%n)%temps-a%tableau(1)%temps)/(a%n-1)

       ti=(ti-a%tableau(1)%temps)/dtot+1

       !    it=idint(ti)
       it=int(ti)

       rat=(ti-it)
       if(insert) then
          do i=1,el%c4%nf
             el%c4%f(i)=(a%tableau(it+1)%volt(i)-a%tableau(it)%volt(i))*rat + a%tableau(it)%volt(i)
             el%c4%ph(i)=(a%tableau(it+1)%phase(i)-a%tableau(it)%phase(i))*rat + a%tableau(it)%phase(i)
             elp%c4%f(i)=(a%tableau(it+1)%volt(i)-a%tableau(it)%volt(i))*rat + a%tableau(it)%volt(i)
             elp%c4%ph(i)=(a%tableau(it+1)%phase(i)-a%tableau(it)%phase(i))*rat + a%tableau(it)%phase(i)
          enddo
       endif
       e_fin=(a%tableau(it+1)%energie-a%tableau(it)%energie)*rat + a%tableau(it)%energie

    endif

  end SUBROUTINE find_acc_energy

  SUBROUTINE set_cavity(tIN,state0,dt)
    IMPLICIT NONE
    TYPE(INTEGRATION_NODE), POINTER  :: TIN,T,t0
    TYPE(internal_state), target  :: state0
    integer n,it,i,nit,count
    real(dp) r,ti,rat,tot,dtot,x(6),tc,dep,dtc,en0,dtc0,e_fin,dt
    real(dp) small_tc, energy0
    type(work) w
    type(acceleration), pointer :: a
    type(element), pointer :: el
    type(elementp), pointer :: elp
    integer :: hh=0
    hh=hh+1
    t=>tin
    nit=1000

    w=t%parent_fibre
    energy0=w%kinetic   !+x_orbit_sync(5)*w%p0c

    el=>t%parent_fibre%mag
    elp=>t%parent_fibre%magp
    a=>el%c4%acc


    n=t%pos_in_fibre-2

    t0=>t     !!!!!!!!!!!!!!!

    a%w1=0
    a%w2=0


    small_tc=(twopi*EL%freq/CLIGHT)**(-1)*1.e-7_dp
    dep=small_tc

    !  write(6,*) small_tc
    ! write(6,*) energy0

    !     energy0=a%e_in(n)
    call find_energy(a%w1,kinetic=energy0)
    energy0=energy0+a%de(n)   ! final energy after that cavity

    call find_energy(a%w2,kinetic=energy0)

    if(a%de(n)/=0.0_dp ) then
       if(mdebug/=0) then
          write(mdebug,*) hh,a%de(n)
          tc=el%c4%t
          dtc=(EL%freq/CLIGHT)**(-1)/50
          do count=-50,50
             x=0.d0;
             el%c4%t=count*dtc
             CALL TRACK_NODE_SINGLE(T,X,STATE0)
             write(mdebug,*) el%c4%t, x(5)*w%p0c,a%de(n)
          enddo
          el%c4%t=tc
       endif
       dtc0=1.e38_dp
       do i=1,nit
          tc=el%c4%t
          x=0.d0;
          CALL TRACK_NODE_SINGLE(T,X,STATE0)
          en0=x(5)*w%p0c

          x=0.d0;
          el%c4%t=tc+dep
          CALL TRACK_NODE_SINGLE(T,X,STATE0)
          !  !  write(6,*) x(5)*w%p0c
          dtc=(x(5)*w%p0c-en0)/dep
          dtc=(a%de(n)-en0)/dtc
          el%c4%t=tc+dtc

          if(i>100) then
             if(abs(dtc)<small_tc.and.abs(dtc)>=dtc0) exit
             dtc0=abs(dtc)
             !    pause 123
          endif
          !    write(6,*) " more ",i,dtc,small_tc
       enddo
       el%c4%t=tc+dtc
       !     elp%c4%t=tc+dtc
       x=0.d0;
       CALL TRACK_NODE_SINGLE(T,X,STATE0)
       dt=x(6)

       if(mdebug/=0) write(mdebug,*) "final tc = ", el%c4%t


       if(i>nit-1) then
          write(6,*) " NO convergence in set_cavity "
          stop 1939
       endif
    endif  ! de=zero
  end SUBROUTINE set_cavity

  SUBROUTINE ORBIT_restore_ANBN
    IMPLICIT NONE
    INTEGER K,I
    TYPE(INTEGRATION_NODE), POINTER  :: T


    DO K=1,my_ORBIT_LATTICE%ORBIT_N_NODE

       T=>my_ORBIT_LATTICE%ORBIT_NODES(K)%NODE


       DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(K)%dpos
          IF(T%PARENT_FIBRE%MAG%slow_ac) THEN
             !         if(old_mod) then
             !          CALL restore_ANBN(T%PARENT_FIBRE%MAG,T%PARENT_FIBRE%MAGP,1)
             !         else
             CALL restore_ANBN_SINGLE(T%PARENT_FIBRE%MAG,T%PARENT_FIBRE%MAGP)
             !        endif
          ENDIF

          T=>T%NEXT
       ENDDO

    enddo
  end SUBROUTINE ORBIT_restore_ANBN





  SUBROUTINE ORBIT_TRACK_NODEP(K,X,STATE)
    IMPLICIT NONE
    type(real_8),  INTENT(INOUT) :: X(6)
    INTEGER K,I,j
    LOGICAL(LP) U
    type(real_8) X5
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), OPTIONAL :: STATE

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       call alloc(x5)
       do i=1,4
          x(i)=x(i)*1.e-3_dp
       enddo
       X5=X(5)
       X(5)=X(6)/my_ORBIT_LATTICE%ORBIT_P0C
       X(6)=X5/my_ORBIT_LATTICE%ORBIT_OMEGA
    ENDIF


    u=my_false

    T=>MY_ORBIT_LATTICE%ORBIT_NODES(K)%NODE
    !    IF(T%USE_TPSA_MAP) THEN !1
    !       DO I=1,6
    !          X(I)=X(I)-T%ORBIT(I)
    !       ENDDO
    !       CALL TRACK(T%TPSA_MAP,X)
    !       if(.not.CHECK_STABLE) then
    !          CALL RESET_APERTURE_FLAG
    !          u=my_true
    !          x(1)=XBIG
    !       endif
    !    ELSE
    DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(K)%dpos
       if(u) exit
       IF(PRESENT(STATE)) THEN
          CALL TRACK_NODE_SINGLE(T,X,STATE)  !,my_ORBIT_LATTICE%ORBIT_CHARGE
       ELSE
          CALL TRACK_NODE_SINGLE(T,X,my_ORBIT_LATTICE%STATE) !,my_ORBIT_LATTICE%ORBIT_CHARGE
       ENDIF
       if(.not.CHECK_STABLE) then
          !          CALL RESET_APERTURE_FLAG
          u=my_true
! EB - comment by Elena Benedetto - Only the 4 transverse coordinates are multiplied by xbig, 
! to avoid loosing the information on the phase (when the mod 2pi is taken)
        !  x(1)=XBIG
          do j=1,4
            x(j)=x(j)*XBIG
          enddo
          exit
       endif
       T=>T%NEXT
    ENDDO
    !    ENDIF

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       do i=1,4
          x(i)=x(i)*1.e3_dp
       enddo
       X5=X(5)
       X(5)=X(6)*my_ORBIT_LATTICE%ORBIT_OMEGA
       X(6)=X5*my_ORBIT_LATTICE%ORBIT_P0C
       call kill(x5)

    ENDIF

  end SUBROUTINE ORBIT_TRACK_NODEP




  SUBROUTINE orbit_to_ptc(x)
    implicit none
    real(dp), intent(inout) :: x(6)
    real(dp) x5

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       x(1:4)=x(1:4)*1.e-3_dp
       X5=X(5)
       X(5)=X(6)/my_ORBIT_LATTICE%ORBIT_P0C
       X(6)=X5/my_ORBIT_LATTICE%ORBIT_OMEGA
    ENDIF

  END SUBROUTINE orbit_to_ptc

  SUBROUTINE ptc_to_orbit(x)
    implicit none
    real(dp), intent(inout) :: x(6)
    real(dp) x5

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       x(1:4)=x(1:4)*1.e3_dp
       X5=X(5)
       X(5)=X(6)*my_ORBIT_LATTICE%ORBIT_OMEGA
       X(6)=X5*my_ORBIT_LATTICE%ORBIT_P0C
    ENDIF

  END SUBROUTINE ptc_to_orbit

  SUBROUTINE ORBIT_MAKE_NODE_LAYOUT_accel(R,no_end_mag)
    IMPLICIT NONE
    TYPE(LAYOUT),TARGET :: R
    TYPE(FIBRE),POINTER :: P
    TYPE(INTEGRATION_NODE),POINTER :: T
    REAL(DP) DL,L,DLMAX,FREQ,CLOSED(6)
    LOGICAL(LP) no_end_mag
    LOGICAL(LP) END_MAG
    INTEGER I,K,NL
    TYPE(INTERNAL_STATE) STATE
    TYPE(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)
    REAL(DP) BET(2),ALF(2),ETA,ETAP
    TYPE(ORBIT_NODE), pointer :: ORBIT_NODES(:)
    logical(lp) doit,cav
     integer mf,i1

    END_MAG=.not.no_end_mag
    ALLOCATE(R%T%ORBIT_LATTICE)
    CALL Set_Up_ORBIT_LATTICE(R%T%ORBIT_LATTICE,0,MY_TRUE)
    my_ORBIT_LATTICE=>R%T%ORBIT_LATTICE
    my_ORBIT_LATTICE%parent_layout=>R

    my_ORBIT_LATTICE%ORBIT_WARNING=0
    CALL GET_FREQ(R,FREQ)
    IF(FREQ/=0.0_dp) THEN
       my_ORBIT_LATTICE%ORBIT_OMEGA=twopi*FREQ/CLIGHT
    ELSE
       WRITE(6,*) " COULD NOT LOCALIZE NON 0.0_dp CAVITY FREQUENCY "
       my_ORBIT_LATTICE%ORBIT_OMEGA=1.0_dp
       my_ORBIT_LATTICE%ORBIT_WARNING=1
    ENDIF
    my_ORBIT_LATTICE%ORBIT_P0C=R%START%mag%P%P0C
    !    my_ORBIT_LATTICE%ORBIT_BETA0=R%START%mag%P%BETA0
    my_ORBIT_LATTICE%ORBIT_BETA0=R%START%BETA0

    NL=7
    K=1
    DLMAX=0.D0
    L=0.0_dp
    DL=0.0_dp
     if(extra_node==1) then  
             k=k+1
     endif 
    T=>R%T%START
    DO I=1,R%T%N-1
       cav=.false.
       !       if(t%parent_fibre%mag%kind==kind4) then
       !          doit=.false.
       !       else
       !          doit=((T%CAS==CASEP1.AND.END_MAG).OR.T%NEXT%S(1)-L>=LMAX)
       !       endif
!!!!!

       doit=((T%CAS==CASEP1.AND.END_MAG).OR.T%NEXT%S(1)-L>=LMAX)
       if(allocated(orbitname)) then
        do i1=1,size(orbitname)
         doit=(t%parent_fibre%mag%name==orbitname(i1).and.T%CAS==CASEP1).or.doit
        enddo
       endif
       if(t%parent_fibre%mag%kind==kind4.or.T%previous%parent_fibre%MAG%KIND==KIND4) then
          DOIT=.TRUE.
       ENDIF
!!!!!
       !       if(T%CAS==CASEP1.AND.t%parent_fibre%mag%kind==kind4) cav=.true.
       !       if(T%CAS==CASEP1.AND.t%previous%parent_fibre%mag%kind==kind4) cav=.true.
       doit=doit.or.cav
       doit=doit.and.I/=1
       IF(doit) THEN    !
          K=K+1
          DL=T%S(1)-L
          !          IF(DL>DLMAX.and.t%previous%parent_fibre%mag%kind/=kind4) DLMAX=DL
          IF(DL>DLMAX) DLMAX=DL
          L=T%S(1)
       ENDIF
       T=>T%NEXT
    ENDDO

    DL=T%S(1)-L
    IF(DL>DLMAX) DLMAX=DL

    K=K+1
    !   if(ORBIT_N_NODE/=0) DEALLOCATE(ORBIT_NODES)
    CALL Set_Up_ORBIT_LATTICE(my_ORBIT_LATTICE,K,MY_TRUE)
    nullify(ORBIT_NODES)
    ORBIT_NODES=>my_ORBIT_LATTICE%ORBIT_NODES

    K=1
    cav=.false.
    L=0.0_dp
    DL=0.0_dp
     if(extra_node==1) then
             ORBIT_NODES(K)%NODE=>T
             CALL ALLOC_ORBIT_NODE1(ORBIT_NODES(K),NL)    
             k=k+1
     endif    
    T=>R%T%START
    ORBIT_NODES(K)%NODE=>T

    !    CALL ALLOC_ORBIT_NODE(ORBIT_NODES,K,NL)
    CALL ALLOC_ORBIT_NODE1(ORBIT_NODES(K),NL)
    IF(t%parent_fibre%mag%kind==kind4) ORBIT_NODES(K)%CAVITY=MY_TRUE;
    DO I=1,R%T%N-1
       cav=.false.
       !       if(t%parent_fibre%mag%kind==kind4) then
       !          doit=.false.
       !       else
       !          doit=((T%CAS==CASEP1.AND.END_MAG).OR.T%NEXT%S(1)-L>=LMAX)
       !       endif
!!!!!

       doit=((T%CAS==CASEP1.AND.END_MAG).OR.T%NEXT%S(1)-L>=LMAX)
       if(allocated(orbitname)) then
        do i1=1,size(orbitname)
         doit=(t%parent_fibre%mag%name==orbitname(i1).and.T%CAS==CASEP1).or.doit
        enddo
       endif
       if(t%parent_fibre%mag%kind==kind4.or.T%previous%parent_fibre%MAG%KIND==KIND4) then
          DOIT=.TRUE.
       ENDIF
!!!!!
       !       if(T%CAS==CASEP1.AND.t%parent_fibre%mag%kind==kind4) cav=.true.
       !       if(T%CAS==CASEP1.AND.t%previous%parent_fibre%mag%kind==kind4) cav=.true.
       doit=doit.or.cav
       doit=doit.and.I/=1
       IF(doit) THEN    !
          K=K+1
          DL=T%S(1)-L
          L=T%S(1)
          ORBIT_NODES(K)%NODE=>T
          ORBIT_NODES(K-1)%LATTICE(1)=DL;
          CALL ALLOC_ORBIT_NODE1(ORBIT_NODES(K),NL)
          !          CALL ALLOC_ORBIT_NODE(ORBIT_NODES,K,NL)
          IF(t%parent_fibre%mag%kind==kind4) ORBIT_NODES(K)%CAVITY=MY_TRUE;
       ENDIF
       T=>T%NEXT
    ENDDO

    K=K+1
    DL=T%S(1)-L
    ORBIT_NODES(K)%NODE=>R%T%START
    ORBIT_NODES(K-1)%LATTICE(1)=DL;
    CALL ALLOC_ORBIT_NODE1(ORBIT_NODES(K),NL)
    !    CALL ALLOC_ORBIT_NODE(ORBIT_NODES,K,NL)

    !   write(6,*) k,DLMAX,LMAX
    !      MY_ORBIT_STATE=>DEFAULT
    my_ORBIT_LATTICE%ORBIT_CHARGE=R%start%CHARGE
    my_ORBIT_LATTICE%ORBIT_N_NODE=k-1
    
    !    write(6,*) size(ORBIT_NODES),my_ORBIT_LATTICE%ORBIT_N_NODE,k

    IF(DLMAX>LMAX) THEN
       WRITE(6,*) " DLMAX > LMAX ",DLMAX
       WRITE(6,*)  "CONSIDER RESPLITTING LATTICE "
       my_ORBIT_LATTICE%ORBIT_WARNING=10+my_ORBIT_LATTICE%ORBIT_WARNING
    ENDIF

    DO K=1,my_ORBIT_LATTICE%ORBIT_N_NODE
       ORBIT_NODES(K)%DPOS=-ORBIT_NODES(K)%NODE%POS+ORBIT_NODES(K+1)%NODE%POS
       IF(ORBIT_NODES(K)%DPOS<=0)ORBIT_NODES(K)%DPOS=ORBIT_NODES(K)%DPOS+r%t%n
       if(ORBIT_NODES(K)%node%parent_fibre%mag%p%p0c/=my_ORBIT_LATTICE%ORBIT_P0C) then
          write(6,*) " element ",ORBIT_NODES(K)%node%parent_fibre%pos,ORBIT_NODES(K)%node%parent_fibre%mag%name
          write(6,*) " has different energy ",ORBIT_NODES(K)%node%parent_fibre%mag%p%p0c,my_ORBIT_LATTICE%ORBIT_P0C
          write(6,*) " fatal "
          stop
       endif
    ENDDO
   if(extra_node==1) ORBIT_NODES(1)%DPOS=0
    P=>R%START
    DO K=1,R%N
       IF(P%PATCH%PATCH/=0) THEN
          IF(ABS(P%PATCH%A_D(3))>my_ORBIT_LATTICE%ORBIT_MAX_PATCH_TZ) THEN
             my_ORBIT_LATTICE%ORBIT_WARNING=100+my_ORBIT_LATTICE%ORBIT_WARNING
             WRITE(6,*)  "LARGE FRONT PATCH "
          ENDIF
          IF(ABS(P%PATCH%B_D(3))>my_ORBIT_LATTICE%ORBIT_MAX_PATCH_TZ) THEN
             my_ORBIT_LATTICE%ORBIT_WARNING=1000+my_ORBIT_LATTICE%ORBIT_WARNING
             WRITE(6,*)  "LARGE BACK PATCH "
          ENDIF
       ENDIF
       P=>P%NEXT
    ENDDO
    !   write(6,*)ORBIT_WARNING, DLMAX,LMAX

    !  COMPUTE LATTICE FUNCTIONS
    my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=MY_FALSE
 !   STATE=(DEFAULT0+NOCAVITY0)
    STATE=(DEFAULT0+NOCAVITY0)  !+time0

    CALL INIT(STATE,1,0,BERZ)
    CALL ALLOC(ID);CALL ALLOC(Y);CALL ALLOC(NORM);

    closed=0.0_dp

 !   CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)
closed(1)=0.001d0
call kanalnummer(mf,"junk.txt")
p=>r%start
do i=1,r%n
    CALL TRACK(R,closed,i,i+1,STATE)
write(mf,*) i,p%mag%name
write(mf,*) closed(1:2)
p=>p%next
enddo
!write(6,*) closed
close(mf)
!pause 123

    ID=1
    Y=CLOSED+ID
    CALL TRACK(R,Y,1,STATE)

    NORM=Y

    Y=CLOSED+NORM%A_T
    BET(1)=(Y(1).SUB.'1')**2+(Y(1).SUB.'01')**2
    BET(2)=(Y(3).SUB.'001')**2+(Y(3).SUB.'0001')**2
    ALF(1)=-((Y(1).SUB.'1')*(Y(2).SUB.'1')+(Y(1).SUB.'01')*(Y(2).SUB.'01'))
    ALF(2)=-((Y(3).SUB.'001')*(Y(4).SUB.'001')+(Y(3).SUB.'0001')*(Y(4).SUB.'0001'))
    ETA=Y(1).SUB.'00001'
    ETAP=Y(2).SUB.'00001'

    WRITE(6,*) "TWISS PARAMETERS AT THE ENTRANCE"
    WRITE(6,*) "BETAS ", BET
    WRITE(6,*) "ALPHAS ",ALF
    WRITE(6,*) "ETA, ETAP ",ETA,ETAP
    DO K=1,my_ORBIT_LATTICE%ORBIT_N_NODE
       CALL ORBIT_TRACK_NODE(K,Y,STATE)
       BET(1)=(Y(1).SUB.'1')**2+(Y(1).SUB.'01')**2
       BET(2)=(Y(3).SUB.'001')**2+(Y(3).SUB.'0001')**2
       ALF(1)=-((Y(1).SUB.'1')*(Y(2).SUB.'1')+(Y(1).SUB.'01')*(Y(2).SUB.'01'))
       ALF(2)=-((Y(3).SUB.'001')*(Y(4).SUB.'001')+(Y(3).SUB.'0001')*(Y(4).SUB.'0001'))
       ETA=Y(1).SUB.'00001'
       ETAP=Y(2).SUB.'00001'
       ORBIT_NODES(K)%LATTICE(2:3)=BET
       ORBIT_NODES(K)%LATTICE(4:5)=ALF
       ORBIT_NODES(K)%LATTICE(6)=ETA
       ORBIT_NODES(K)%LATTICE(7)=ETAP
    ENDDO
    WRITE(6,*) "TWISS PARAMETERS AT THE EXIT"
    WRITE(6,*) "BETAS ", BET
    WRITE(6,*) "ALPHAS ",ALF
    WRITE(6,*) "ETA, ETAP ",ETA,ETAP


    my_ORBIT_LATTICE%ORBIT_L=r%t%end%s(1)
    my_ORBIT_LATTICE%ORBIT_gammat=norm%tune(3)/my_ORBIT_LATTICE%ORBIT_L
    if(my_ORBIT_LATTICE%ORBIT_gammat<0) then
       my_ORBIT_LATTICE%ORBIT_gammat=-sqrt(-1.0_dp/my_ORBIT_LATTICE%ORBIT_gammat)
    else
       my_ORBIT_LATTICE%ORBIT_gammat=sqrt(1.0_dp/my_ORBIT_LATTICE%ORBIT_gammat)
    endif

    my_ORBIT_LATTICE%ORBIT_mass_in_amu=r%start%mass*pmae_amu/pmae
    my_ORBIT_LATTICE%orbit_kinetic=sqrt(r%start%mass**2+my_ORBIT_LATTICE%ORBIT_P0C**2) &
         -r%start%mass
    my_ORBIT_LATTICE%orbit_harmonic=my_ORBIT_LATTICE%ORBIT_L*my_ORBIT_LATTICE%ORBIT_OMEGA/TWOPI/my_ORBIT_LATTICE%ORBIT_BETA0
    my_ORBIT_LATTICE%ORBIT_LMAX=DLMAX
    !my_ORBIT_LATTICE%ORBIT_BETA0
    my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=MY_TRUE
    CALL KILL(ID);CALL KILL(Y);CALL KILL(NORM);

    my_ORBIT_LATTICE%orbit_harmonic=nint(my_ORBIT_LATTICE%orbit_harmonic)
    my_ORBIT_LATTICE%ORBIT_OMEGA=my_ORBIT_LATTICE%orbit_harmonic/my_ORBIT_LATTICE%ORBIT_L*twopi*my_ORBIT_LATTICE%ORBIT_BETA0
    my_ORBIT_LATTICE%ORBIT_OMEGA_after=my_ORBIT_LATTICE%ORBIT_OMEGA

    w1_orbit=r%start
    my_ORBIT_LATTICE%orbit_brho=w1_orbit%brho
    my_ORBIT_LATTICE%orbit_energy=w1_orbit%energy
    my_ORBIT_LATTICE%orbit_gamma=1.0_dp/w1_orbit%gamma0i
    !    my_ORBIT_LATTICE%orbit_dppfac=one/sqrt(w1_orbit%beta0)/w1_orbit%energy
    my_ORBIT_LATTICE%orbit_deltae=0.0_dp;

    write(6,*) my_ORBIT_LATTICE%ORBIT_L
    write(6,*) my_ORBIT_LATTICE%ORBIT_OMEGA
    write(6,*) my_ORBIT_LATTICE%orbit_beta0
    write(6,*) my_ORBIT_LATTICE%orbit_brho
    write(6,*) my_ORBIT_LATTICE%orbit_p0c
    write(6,*) my_ORBIT_LATTICE%orbit_energy
    write(6,*) my_ORBIT_LATTICE%orbit_kinetic
    write(6,*) my_ORBIT_LATTICE%orbit_gamma
    write(6,*) r%start%mass,w1_orbit%mass
    my_ORBIT_LATTICE%state=default

    if(allocated(orbitname)) deallocate(orbitname)
  END SUBROUTINE ORBIT_MAKE_NODE_LAYOUT_accel


  SUBROUTINE update_twiss_for_orbit 
    IMPLICIT NONE
    type(layout), pointer :: R
    REAL(DP) CLOSED(6)

    INTEGER K
    TYPE(INTERNAL_STATE) STATE
    TYPE(DAMAP) ID
    TYPE(NORMALFORM) NORM
    TYPE(REAL_8) Y(6)
    REAL(DP) BET(2),ALF(2),ETA,ETAP
    TYPE(ORBIT_NODE), pointer :: ORBIT_NODES(:)

    

      !  COMPUTE LATTICE FUNCTIONS
    r=>my_ORBIT_LATTICE%parent_layout
    ORBIT_NODES=>my_ORBIT_LATTICE%ORBIT_NODES
    my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=MY_FALSE
 
    STATE=(my_ORBIT_LATTICE%state-time0)+delta0   

    CALL INIT(STATE,1,0,BERZ)
    CALL ALLOC(ID);CALL ALLOC(Y);CALL ALLOC(NORM);

    closed=0.0_dp
    CALL FIND_ORBIT(R,CLOSED,1,STATE,1e-5_dp)


    ID=1
    Y=CLOSED+ID
    CALL TRACK(R,Y,1,STATE)

    NORM=Y

    Y=CLOSED+NORM%A_T
    BET(1)=(Y(1).SUB.'1')**2+(Y(1).SUB.'01')**2
    BET(2)=(Y(3).SUB.'001')**2+(Y(3).SUB.'0001')**2
    ALF(1)=-((Y(1).SUB.'1')*(Y(2).SUB.'1')+(Y(1).SUB.'01')*(Y(2).SUB.'01'))
    ALF(2)=-((Y(3).SUB.'001')*(Y(4).SUB.'001')+(Y(3).SUB.'0001')*(Y(4).SUB.'0001'))
    ETA=Y(1).SUB.'00001'
    ETAP=Y(2).SUB.'00001'

    WRITE(6,*) "TWISS PARAMETERS AT THE ENTRANCE"
    WRITE(6,*) "BETAS ", BET
    WRITE(6,*) "ALPHAS ",ALF
    WRITE(6,*) "ETA, ETAP ",ETA,ETAP
    DO K=1,my_ORBIT_LATTICE%ORBIT_N_NODE
       CALL ORBIT_TRACK_NODE(K,Y,STATE)
       BET(1)=(Y(1).SUB.'1')**2+(Y(1).SUB.'01')**2
       BET(2)=(Y(3).SUB.'001')**2+(Y(3).SUB.'0001')**2
       ALF(1)=-((Y(1).SUB.'1')*(Y(2).SUB.'1')+(Y(1).SUB.'01')*(Y(2).SUB.'01'))
       ALF(2)=-((Y(3).SUB.'001')*(Y(4).SUB.'001')+(Y(3).SUB.'0001')*(Y(4).SUB.'0001'))
       ETA=Y(1).SUB.'00001'
       ETAP=Y(2).SUB.'00001'
       ORBIT_NODES(K)%LATTICE(2:3)=BET
       ORBIT_NODES(K)%LATTICE(4:5)=ALF
       ORBIT_NODES(K)%LATTICE(6)=ETA
       ORBIT_NODES(K)%LATTICE(7)=ETAP
    ENDDO
    WRITE(6,*) "TWISS PARAMETERS AT THE EXIT"
    WRITE(6,*) "BETAS ", BET
    WRITE(6,*) "ALPHAS ",ALF
    WRITE(6,*) "ETA, ETAP ",ETA,ETAP

    my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=my_true
 
    end SUBROUTINE update_twiss_for_orbit 

end module orbit_ptc

subroutine ptc_track_particle(node_index, x,xp,y,yp,phi,dE)

  USE orbit_ptc
  IMPLICIT NONE
  REAL(DP) x,xp,y,yp,phi,dE
  INTEGER node_index
  INTEGER i

  i = node_index + 1


  call PUT_RAY(x,xp,y,yp,phi,dE)

  call TRACK_ONE_NODE(i)

  call GET_RAY(x,xp,y,yp,phi,dE)

  IF(I==1.AND.MF_HERD/=0) THEN
     WRITE(MF_HERD,'(4(1X,E15.8))') PHI,DE,my_ORBIT_LATTICE%orbit_p0c, &
          x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*1e3_dp
     !       WRITE(MF_HERD,'(6(1X,E15.8))') PHI,DE,X_ORBIT(6),X_ORBIT(5), &
     !x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*1000.d0,my_ORBIT_LATTICE%ORBIT_OMEGA
     !        ,my_ORBIT_LATTICE%ORBIT_P0C
  ENDIF

  return
end subroutine ptc_track_particle
