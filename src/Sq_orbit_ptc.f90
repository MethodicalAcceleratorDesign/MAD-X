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
  REAL(dp) :: x_orbit_sync(6)= zero,dt_orbit_sync=zero
  type(work)  :: w1_orbit,w2_orbit
  type(fibre), pointer :: p_orbit
  integer :: ptc_node_old=-1
  logical(lp) :: accelerate=my_false, first_particle=my_false
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
    INTEGER X
    X=my_ORBIT_LATTICE%ORBIT_CHARGE
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

  SUBROUTINE ORBIT_TRACK_NODE_fake(K,X,STATE)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    INTEGER K,I
    REAL(DP) X5
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), OPTIONAL :: STATE

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       x(1:4)=x(1:4)*1.e-3_dp
       X5=X(5)
       X(5)=X(6)/my_ORBIT_LATTICE%ORBIT_P0C
       X(6)=X5/my_ORBIT_LATTICE%ORBIT_OMEGA
    ENDIF




    T=>my_ORBIT_LATTICE%ORBIT_NODES(K)%NODE
    DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(K)%dpos
       IF(PRESENT(STATE)) THEN
          !          CALL TRACK_NODE_SINGLE(T,X,STATE) !,my_ORBIT_LATTICE%ORBIT_CHARGE
          if(state%time) then
             x(6)=x(6)+STATE%totalpath*T%previous%s(5)/t%parent_fibre%beta0
          else
             x(6)=x(6)+STATE%totalpath*T%previous%s(5)
          endif
       ELSE
          if(my_ORBIT_LATTICE%STATE%time) then
             x(6)=x(6)+my_ORBIT_LATTICE%STATE%totalpath*T%previous%s(5)/t%parent_fibre%beta0
          else
             x(6)=x(6)+my_ORBIT_LATTICE%STATE%totalpath*T%previous%s(5)
          endif
          !          CALL TRACK_NODE_SINGLE(T,X,my_ORBIT_LATTICE%STATE) !,my_ORBIT_LATTICE%ORBIT_CHARGE
       ENDIF
       T=>T%NEXT
    ENDDO

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       x(1:4)=x(1:4)*1.e3_dp
       X5=X(5)
       X(5)=X(6)*my_ORBIT_LATTICE%ORBIT_OMEGA
       X(6)=X5*my_ORBIT_LATTICE%ORBIT_P0C
    ENDIF



  end SUBROUTINE ORBIT_TRACK_NODE_fake


  SUBROUTINE TRACK_NODE_fake_totalpath_half(T,X,STATE,after)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), target :: STATE
    logical(lp) after
    real(dp) beta0
    type(work) w
    integer nst,nc

    beta0=t%PARENT_FIBRE%beta0
    nst=t%pos_in_fibre-2

    if(t%PARENT_FIBRE%mag%kind==kind4) then

       if(after.and.t%cas==case0) then
          w=t%parent_fibre
          x(5)= x(5) + t%parent_fibre%mag%c4%acc%de(nst)/w%p0c
       endif !  kind4

    endif

    if(state%time) then
       x(6)=x(6)+T%ds_AC*half*(one/beta0+x(5))/root(one+2*x(5)/beta0+x(5)**2)
    else
       x(6)=x(6)+T%ds_AC*half

    endif

  end SUBROUTINE TRACK_NODE_fake_totalpath_half

  SUBROUTINE TRACK_NODE_fake_totalpath_half_plain(T,X,STATE)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    TYPE(INTEGRATION_NODE), POINTER  :: T
    TYPE(INTERNAL_STATE), target :: STATE

    if(state%time) then
       x(6)=x(6)+T%ds_AC*half/t%PARENT_FIBRE%beta0
    else
       x(6)=x(6)+T%ds_AC*half
    endif

  end SUBROUTINE TRACK_NODE_fake_totalpath_half_plain

  SUBROUTINE ORBIT_TRACK_NODE_Standard_R(K,X,STATE)
    IMPLICIT NONE
    REAL(DP),  INTENT(INOUT) :: X(6)
    INTEGER K,I
    LOGICAL(LP) U,cav
    REAL(DP) X5,dt
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
    !    IF(T%USE_TPSA_MAP) THEN  ! 2
    !       X=X-T%ORBIT
    !       CALL TRACK(T%TPSA_MAP,X)
    !       if(.not.CHECK_STABLE) then
    !          CALL RESET_APERTURE_FLAG
    !          u=my_true
    !          x(1)=XBIG
    !       endif
    !    ELSE

    if(present(state)) then
       state0=>state
    else
       state0=>my_ORBIT_LATTICE%STATE
    endif

    IF(STATE0%MODULATION) THEN !modulate
       if(K/=ptc_node_old) then !modulate
          ptc_node_old=k !modulate
          xsm=xsm0 !modulate
          first_particle=.true.
       endif !modulate
       xsm0=xsm !modulate
    ENDIF !modulate

    IF(accelerate) THEN !accelerate
       if(K/=ptc_node_old) then !accelerate
          ptc_node_old=k !accelerate
          first_particle=.true.
       endif !accelerate
    ENDIF !accelerate




    DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(K)%dpos
       if(STATE0%MODULATION) then !modulate
          if(t%parent_fibre%mag%slow_ac) CALL MODULATE(T,XSM0,STATE0) !modulate
          CALL TRACK_MODULATION(T,XSM0,STATE0) !modulate
       endif !modulate


       if(first_particle.and.accelerate) then !accelerate
          if(t%previous%parent_fibre%mag%kind==kind4) then
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
             if(cav) call set_cavity(t,state0,dt_orbit_sync)
          endif !accelerate

          call TRACK_NODE_fake_totalpath_half(T,x_orbit_sync,STATE0,my_false)  ! accelerate
          !            write(6,*)t%parent_fibre%mag%name,t%parent_fibre%pos,t%pos, x_orbit_sync(6)
          !            pause 777
       endif !accelerate

       if(u) exit
       CALL TRACK_NODE_SINGLE(T,X,STATE0) !,my_ORBIT_LATTICE%ORBIT_CHARGE
       !   WRITE(MFFF,*) t%parent_fibre%mag%name,t%pos_in_fibre
       !   WRITE(MFFF,*) xsm0%ac%t
       !   WRITE(MFFF,*) X(1),X(2)
       if(.not.CHECK_STABLE) then
          CALL RESET_APERTURE_FLAG
          u=my_true
          x(1)=XBIG
          if(wherelost==1) then
             t%lost=t%lost+1
          endif
          exit
       endif
       if(STATE0%MODULATION) then !modulate
          IF(T%PARENT_FIBRE%MAG%slow_ac) THEN !modulate
             !        if(old_mod) then
             !          CALL restore_ANBN(T%PARENT_FIBRE%MAG,T%PARENT_FIBRE%MAGP,1) !modulate
             !        else
             CALL restore_ANBN_SINGLE(T%PARENT_FIBRE%MAG,T%PARENT_FIBRE%MAGP)
             !        endif
          ENDIF !modulate
       endif  !modulate

       if(first_particle.and.accelerate) then !accelerate
          call TRACK_NODE_fake_totalpath_half(T,x_orbit_sync,STATE0,my_true)  ! accelerate
       endif

       if(t%parent_fibre%mag%kind==kind4.and.accelerate) then
          if(t%cas==case0) then
             x(6)=x(6)-dt_orbit_sync
             call ORBIT_up_grade_x(x,t)
             if(first_particle) call ORBIT_up_grade_x(x_orbit_sync,t)
             !       t%parent_fibre=t%parent_fibre%mag%c4%acc%w2
          endif

          !   if(associated(t%parent_fibre%t2,t)) then
          !     call ORBIT_up_grade_mag(t)
          !   endif
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


  subroutine energize_ORBIT_lattice
    implicit none
    type(fibre), pointer :: p
    type(work) werk,travail
    real(dp) e_in
    integer i


    werk=0
    travail=0

    write(6,*) "energize ", paccfirst%mag%c4%freq

    call find_acc_energy(paccfirst,x_orbit_sync(6),e_in,my_false)
    call find_energy(werk,kinetic=e_in)



    p=>paccfirst%parent_layout%start
    do i=1,paccfirst%parent_layout%n
       travail=p
       if(p%mag%kind==kind4) then
          p%mag%c4%freq  = p%mag%c4%freq*werk%beta0/travail%beta0
          p%magp%c4%freq = p%mag%c4%freq
       endif
       p=werk
       p=>p%next
    enddo

    write(6,*) "energize 2", paccfirst%mag%c4%freq
    write(6,*) werk%beta0,travail%beta0

  end subroutine energize_ORBIT_lattice

  subroutine ORBIT_up_grade_mag(t)
    implicit none
    type(fibre), pointer :: p
    integer i
    TYPE(INTEGRATION_NODE), POINTER  :: T
    type(acceleration), pointer :: a

    a=>t%parent_fibre%mag%c4%acc


    p=>t%parent_fibre

    do i=1,p%parent_layout%n  !-1
       if(p%mag%kind==kind4) then
          p%mag%c4%freq  = p%mag%c4%freq*a%w2%beta0/a%w1%beta0
          p%magp%c4%freq = p%mag%c4%freq
       endif
       p=a%w2
       p=>p%next
    enddo

  end subroutine ORBIT_up_grade_mag

  subroutine ORBIT_up_grade_x(x,t)
    implicit none
    real(dp),intent(inout) :: x(6)
    TYPE(INTEGRATION_NODE), POINTER  :: T
    type(acceleration), pointer :: a

    a=>t%parent_fibre%mag%c4%acc

    X(2)=X(2)*a%w1%P0C/a%w2%P0C
    X(4)=X(4)*a%w1%P0C/a%w2%P0C
    X(5)=root(one+two*X(5)/a%w1%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
    X(5)=X(5)*a%w1%P0C/a%w2%P0C-one !X(5) = DP/P0C_NEW
    X(5)=(two*X(5)+X(5)**2)/(root(one/a%w2%BETA0**2+two*X(5)+X(5)**2) &
         +one/a%w2%BETA0)
    !if(enforce_zero_x5) x(5)=zero
    !    x(6)=x(6)*my_ORBIT_LATTICE%orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA

    !write(6,*) orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA
  END subroutine ORBIT_up_grade_x

  SUBROUTINE find_all_energies(t,state0)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(INTEGRATION_NODE), POINTER  :: T,c
    TYPE(internal_state), target  :: state0
    integer i
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
                call find_acc_energy(c%parent_fibre,t_fin,e_fin,my_false)  ! to insure ridiculous self-consistancy

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
          ipause=mypause(1)
          stop 1950
       endif


    enddo
  end SUBROUTINE find_all_energies

  SUBROUTINE find_one_turn_final_energy(t,state0,e_fin)
    IMPLICIT NONE
    integer ipause, mypause
    TYPE(INTEGRATION_NODE), POINTER  :: T,c
    TYPE(internal_state), target  :: state0
    integer i
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
          ipause=mypause(1)
          stop 1951
       endif


    enddo

    !  write(6,*) " one turn cav x(6) ",xsync(6)
  end SUBROUTINE find_one_turn_final_energy



  SUBROUTINE find_acc_energy(p,t_fin,e_fin,insert)
    IMPLICIT NONE
    integer ipause, mypause
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


    ti=t_fin/clight/a%unit_time    ! time in milliseconds

    if(ti>a%tableau(a%n)%temps.or.ti<a%tableau(1)%temps) then
       write(6,*) " Tracking must stop time =",ti
       ipause=mypause(1947)
       stop 1947
    endif

    dtot=(a%tableau(a%n)%temps-a%tableau(1)%temps)/(a%n-1)

    ti=(ti-a%tableau(1)%temps)/dtot+1

    it=idint(ti)

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
    !   write(6,*) " e_fin as a function of time only ", e_fin
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

    if(a%de(n)/=zero ) then
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

          if(i>10) then
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
    INTEGER K,I
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
          x(1)=XBIG
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


    END_MAG=.not.no_end_mag
    ALLOCATE(R%T%ORBIT_LATTICE)
    CALL Set_Up_ORBIT_LATTICE(R%T%ORBIT_LATTICE,0,MY_TRUE)
    my_ORBIT_LATTICE=>R%T%ORBIT_LATTICE

    my_ORBIT_LATTICE%ORBIT_WARNING=0
    CALL GET_FREQ(R,FREQ)
    IF(FREQ/=ZERO) THEN
       my_ORBIT_LATTICE%ORBIT_OMEGA=twopi*FREQ/CLIGHT
    ELSE
       WRITE(6,*) " COULD NOT LOCALIZE NON ZERO CAVITY FREQUENCY "
       my_ORBIT_LATTICE%ORBIT_OMEGA=ONE
       my_ORBIT_LATTICE%ORBIT_WARNING=1
    ENDIF
    my_ORBIT_LATTICE%ORBIT_P0C=R%START%mag%P%P0C
    !    my_ORBIT_LATTICE%ORBIT_BETA0=R%START%mag%P%BETA0
    my_ORBIT_LATTICE%ORBIT_BETA0=R%START%BETA0

    NL=7
    K=1
    DLMAX=0.D0
    L=ZERO
    DL=ZERO
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
    L=ZERO
    DL=ZERO
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
    STATE=(DEFAULT0+NOCAVITY0)

    CALL INIT(STATE,1,0,BERZ)
    CALL ALLOC(ID);CALL ALLOC(Y);CALL ALLOC(NORM);

    closed=zero
    CALL FIND_ORBIT(R,CLOSED,1,STATE,c_1d_5)


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
       my_ORBIT_LATTICE%ORBIT_gammat=-sqrt(-one/my_ORBIT_LATTICE%ORBIT_gammat)
    else
       my_ORBIT_LATTICE%ORBIT_gammat=sqrt(one/my_ORBIT_LATTICE%ORBIT_gammat)
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
    my_ORBIT_LATTICE%orbit_gamma=one/w1_orbit%gamma0i
    !    my_ORBIT_LATTICE%orbit_dppfac=one/sqrt(w1_orbit%beta0)/w1_orbit%energy
    my_ORBIT_LATTICE%orbit_deltae=zero;

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


  END SUBROUTINE ORBIT_MAKE_NODE_LAYOUT_accel



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
          x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*c_1d3
     !       WRITE(MF_HERD,'(6(1X,E15.8))') PHI,DE,X_ORBIT(6),X_ORBIT(5), &
     !x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*1000.d0,my_ORBIT_LATTICE%ORBIT_OMEGA
     !        ,my_ORBIT_LATTICE%ORBIT_P0C
  ENDIF

  return
end subroutine ptc_track_particle
