module accel_ptc
  use beam_beam_ptc
  ! use orbit_ptc
  implicit none
  public
  integer n_cav_accel
  real(dp) :: vrf=210.d-3   !mev
  real(dp),parameter :: dphsa_final=30.0_dp*deg_to_rad_  ! 30.0_dp*deg_to_rad_
  real(dp) :: dphsa_by_dt=dphsa_final/(100.e-3_dp*CLIGHT)
  real(dp) :: harmon=9.d0
  real(dp) :: freq=0.d0
  real(dp) :: circum=0.d0
  real(dp) :: total_cav_L=0.d0
  logical(lp) :: enforce_zero_x5=my_false
  logical(lp) :: energy_reached=my_false
  real(dp) :: present_kinetic_energy=0.d0
  real(dp) :: final_kinetic_energy=40.d0  !3.1d0

  !   Other mode of acceleration
  real(dp) :: T_parabolic=100.e-3_dp*CLIGHT
  real(dp) :: T_rise=13.e-3_dp*CLIGHT,T_END_LINEAR=630091002.783574d0,T_END=660070248.583574d0
  real(dp) :: p0c0
  real(dp) :: p0c1 , kinetic1
  real(dp) :: p01_ratio = .181_dp/.143_dp
  real(dp) :: p0c2,kinetic2,p0cf
  real(dp) :: as1,as2,slope12
  real(dp) :: vrf0=40.d-3,vrff=280.d-3, slopev    !mev
  integer :: count_cav
  TYPE fibre_array
     type(fibre), pointer :: p
     integer, pointer :: pos
  END TYPE fibre_array
  type(fibre_array), allocatable ::cavity_location(:)
  TYPE accel_data
     real(dp), pointer :: p0c
     real(dp), pointer :: volt(:)
     real(dp), pointer :: phase(:)
  END TYPE accel_data
  type(accel_data), allocatable :: accel(:)

  ! accel_ptc_locate_accel_cavity
  ! Locates the cavities and store then in cavity_location: Pointer to fibre and integer position
  !
  ! accel_ptc_design_no_accel
  !
  !  Sets the phase at the cavity and also c_%phase0 so that VRF>0 => acceleration and stability
  !
  !
  !
  !

contains

  !!  common to both methods
  !call make_table(w1%p0c,"ACCWAVE_40kV_280kV_350ms.DAT")

  subroutine make_table(filename)
    implicit none
    character*(*) filename
    real(dp)  x
    real(dp) p0c,b0
    integer k,i,mf,mode,t

    p0c=my_ORBIT_LATTICE%orbit_p0c
    my_ORBIT_LATTICE%accel=my_true
    call kanalnummer(mf)

    open(unit=mf,file=filename)

    k=0
    do i=1,1000000
       read(mf,*,end=1) x
       k=k+1
    enddo
1   continue
    write(6,*)  k ," lines"

    close(mf)

    allocate(accel(0:k-2))

    open(unit=mf,file=filename)

    read(mf,*) mode
    mode=1
    do i=0,k-2
       allocate(accel(i)%p0c)
       allocate(accel(i)%volt(mode))
       allocate(accel(i)%phase(mode))
       read(mf,* ) t,accel(i)%p0c,accel(i)%volt(1),accel(i)%phase(1)
       if(i==0) b0=accel(i)%p0c
       accel(i)%p0c=accel(i)%p0c/b0*p0c
       accel(i)%volt(1)=accel(i)%volt(1)*1.e-3_dp
    enddo

    close(mf)

  end subroutine make_table

  subroutine get_from_table_volt(time,p0,volt,phase)
    implicit none
    integer ti
    real(dp) tr,volt,phase,time,p0

    !write(6,*) size(accel)

    tr=time/clight*1000
    ti=int(tr)
    !write(6,*) tr,ti,accel(ti)%p0c

    if(ti>size(accel)) stop

    volt= (tr-ti)*(accel(ti+1)%volt(1)-accel(ti)%volt(1)) + accel(ti)%volt(1)
    phase= (tr-ti)*(accel(ti+1)%phase(1)-accel(ti)%phase(1)) + accel(ti)%phase(1)
    p0= (tr-ti)*(accel(ti+1)%p0c -accel(ti)%p0c) + accel(ti)%p0c





  end subroutine get_from_table_volt


  subroutine accel_ORBIT_up_grade_x(x)
    implicit none
    real(dp),intent(inout) :: x(6)

    call orbit_to_ptc(x)

    X(2)=X(2)*w1_ORBIT%P0C/w2_ORBIT%P0C
    X(4)=X(4)*w1_ORBIT%P0C/w2_ORBIT%P0C
    X(5)=root(one+two*X(5)/w1_ORBIT%BETA0+X(5)**2)  !X(5) = 1+DP/P0C_OLD
    X(5)=X(5)*w1_ORBIT%P0C/w2_ORBIT%P0C-one !X(5) = DP/P0C_NEW
    X(5)=(two*X(5)+X(5)**2)/(root(one/w2_ORBIT%BETA0**2+two*X(5)+X(5)**2) &
         +one/w2_ORBIT%BETA0)
    !if(enforce_zero_x5) x(5)=zero
    if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) x(6)=x(6)*my_ORBIT_LATTICE%orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA
    call ptc_to_orbit(x)
    !write(6,*) orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA
  END subroutine accel_ORBIT_up_grade_x

  subroutine accel_ORBIT_up_grade_mag
    implicit none
    type(fibre), pointer :: p
    integer i


    p=>P_orbit%next

    do i=1,p%parent_layout%n-1
       p=w2_ORBIT
       p=>p%next
    enddo

  end subroutine accel_ORBIT_up_grade_mag

  subroutine accel_ORBIT_up_grade_mag_all
    implicit none
    type(fibre), pointer :: p
    type(LAYOUT), pointer :: L
    integer i


    L=>my_ORBIT_LATTICE%ORBIT_NODES(1)%node%parent_fibre%parent_layout
    P=>L%START
    do i=1,L%n

       ! if(p%mag%kind==kind4) then
       !   p%mag%freq= my_ORBIT_LATTICE%freqa
       !   p%mag%volt= my_ORBIT_LATTICE%volta
       !   p%mag%phas= my_ORBIT_LATTICE%phasa
       !   p%magp%freq=my_ORBIT_LATTICE%freqa
       !   p%magp%volt=my_ORBIT_LATTICE%volta
       !   p%magp%phas=my_ORBIT_LATTICE%phasa
       ! endif

       p=w2_ORBIT

       p=>p%next
    enddo


  end subroutine accel_ORBIT_up_grade_mag_all




!!!!!!!!!!! fake orbit !!!!!!!!!!!!!!!!!!

  subroutine accel_orbit_beam(ring)
    !use accel_ptc
    implicit none
    integer n_turn,i,k,npart,j,kk
    type(layout), target :: ring
    TYPE(BEAM),target :: RAYS
    real(dp) sig0(6),x(6)
    integer mf
    type(fibre), pointer :: p
    npart=9
    sig0=1.e-6

    CALL create_beam(RAYS,npart,2.D0,SIG0)

    my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=.true.
    open(unit=mf, file='wp2_029_300kW_J54_INJALL_TRM_210kVM_Bf016.dat')
    if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) then
       do i=1,RAYS%n
          read(mf,*) RAYS%x(i,1:4), RAYS%x(i,5), RAYS%x(i,6)
       enddo
    else
       stop 999
       do i=1,RAYS%n
          read(mf,*) RAYS%x(i,1:4), RAYS%x(i,6), RAYS%x(i,5)
          RAYS%x(i,1:4)=RAYS%x(i,1:4)/1.0e3_dp
          RAYS%x(i,5)=RAYS%x(i,5)/my_ORBIT_LATTICE%orbit_p0c
          RAYS%x(i,6)=RAYS%x(i,6)/my_ORBIT_LATTICE%ORBIT_omega
       enddo
    endif

    close(mf)

    !call make_table("NOACC_ACC_210.DAT")
    call make_table("RF_Pattern_210kV_INJ120mc_ACC350ms.DAT")
    !call make_table("ACCWAVE_40kV_280kV_350ms.DAT")
    !call make_table("ACCWAVE_210KVH9_350ms.DAT")
    !call make_table("noaccel.DAT")

    write(6,*) my_ORBIT_LATTICE%ORBIT_harmonic
    write(6,*) my_ORBIT_LATTICE%ORBIT_omega
    write(6,*) c_%phase0
    default=my_ORBIT_LATTICE%state

    IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       call kanalnummer(mf)
       open(unit=mf,file='crazy_orbit_unit.dat')
    else
       stop 998
       call kanalnummer(mf)
       open(unit=mf,file='crazy_ptc_unit.dat')
    endif
    WRITE(6,*) " n_TURN "
    READ(5,*) N_TURN
    !n_turn=20

    !rays%x(2,1:6)=0.d0
    !rays%x(2,6)=0.00001d0

    call ptc_synchronous_set(-1)
    write(6,*) " reading rays after some turns ?"
    read(5,*) i
    if(i==1) then
       call ptc_synchronous_set(-2)
       if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) then
          do i=1,RAYS%n
             read(mf,*) RAYS%x(i,1:4), RAYS%x(i,5), RAYS%x(i,6)
          enddo
       else
          stop 997
          do i=1,RAYS%n
             read(mf,*) RAYS%x(i,1:4), RAYS%x(i,6), RAYS%x(i,5)
             RAYS%x(i,1:4)=RAYS%x(i,1:4)/1.0e3_dp
             RAYS%x(i,5)=RAYS%x(i,5)/my_ORBIT_LATTICE%orbit_p0c
             RAYS%x(i,6)=RAYS%x(i,6)/my_ORBIT_LATTICE%ORBIT_omega
          enddo
       endif

       close(mf)
       IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
          call kanalnummer(mf)
          open(unit=mf,file='crazy_orbit_unit2.dat')
       else
          stop 996
          call kanalnummer(mf)
          open(unit=mf,file='crazy_ptc_unit2.dat')
       endif

    endif

    !       write(6,*)" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
    !       write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
    !       write(6,*)my_ORBIT_LATTICE%freqa,my_ORBIT_LATTICE%freqb
    !!       write(6,*)my_ORBIT_LATTICE%volta,my_ORBIT_LATTICE%voltb
    !      write(6,*)my_ORBIT_LATTICE%phasa,my_ORBIT_LATTICE%phasb
    !      write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
    !      write(6,*)" $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ "
    call print(my_ORBIT_LATTICE%state,6)

    do j=1,n_turn
       !if(j==2) my_ORBIT_LATTICE%accel=.false.
       if(mod(j,100)==0) then
          write(6,*) " turn ",j,rays%x(1,6),x_orbit_sync(6)
       endif
       do i=0,my_ORBIT_LATTICE%ORBIT_N_NODE-1

          call ptc_synchronous_set(i)
          do k=1,rays%n

             IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN

                call ptc_track_particle(i,rays%x(k,1),rays%x(k,2),rays%x(k,3),rays%x(k,4),&
                     rays%x(k,5),rays%x(k,6))
             else
                stop 990
                call ptc_track_particle(i,rays%x(k,1),rays%x(k,2),rays%x(k,3),rays%x(k,4),&
                     rays%x(k,5),rays%x(k,6))

             endif

          enddo

          call ptc_synchronous_after(i)

       enddo
       write(25,*)  rays%x(1,5),rays%x(1,6)
       kk=2
       !IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
       ! write(mf,'(E15.8,1x,E15.8,1x)') rays%x(kk,5)/my_ORBIT_LATTICE%ORBIT_OMEGA,rays%x(kk,6)/my_ORBIT_LATTICE%ORBIT_P0C
    enddo
    call ptc_synchronous_after(-2)

    DO KK=1,RAYS%N
       !write(mf,'(6(E15.8,1x))') rays%x(kk,1:6)
       write(mf,*) rays%x(kk,1:6)
    ENDDO
    !do kk=1,rays%n
    !IF(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) THEN
    ! write(mf,'(6(E15.8,1x))') rays%x(kk,1:6)
    ! write(mf,'(E15.8,1x,E15.8,1x)') rays%x(kk,5)/my_ORBIT_LATTICE%ORBIT_OMEGA,rays%x(kk,6)/my_ORBIT_LATTICE%ORBIT_P0C
    !else
    ! write(mf,'(E15.8,1x,E15.8,1x)') rays%x(kk,6),rays%x(kk,5)
    !endif
    !enddo
    close(mf)
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(6,*) x_orbit_sync
    write(6,*) rays%x(1,1:6)
    !write(6,*) rays%x(2,1:6)
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    write(6,*) my_ORBIT_LATTICE%ORBIT_harmonic
    write(6,*) my_ORBIT_LATTICE%ORBIT_omega
    write(6,*) c_%phase0
    write(6,*) "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    !stop
  end subroutine accel_orbit_beam



end module accel_ptc




!===========================================================
! This subroutine should be called before particle tracking.
! It tells PTC to do something related to acceleration.
!
!
!
!===========================================================
SUBROUTINE ptc_synchronous_set(i_node)

  USE accel_ptc     !,vrff=>vrf,freqf=>freq
  IMPLICIT NONE
  INTEGER  i_node
  INTEGER  i_node1,i,mf,j
  type(internal_state) state
  real(dp) p0,vrfx,dphat,freqf,dphase0,dphase
  TYPE(INTEGRATION_NODE), POINTER  :: T


  i_node1 = i_node + 1

  if(i_node==-1) then
     x_orbit_sync=zero
     w1_orbit=0
     w2_orbit=0
     my_ORBIT_LATTICE%orbit_omega_after=my_ORBIT_LATTICE%ORBIT_omega
     my_ORBIT_LATTICE%state=my_ORBIT_LATTICE%state+time0
     default=my_ORBIT_LATTICE%state
     my_ORBIT_LATTICE%first=my_true
     my_ORBIT_LATTICE%xs6=-1.e38_dp  ! very big number
     !for speed
     !     call PUT_state(default,my_ORBIT_LATTICE%ORBIT_NODES(1)%node%PARENT_FIBRE%PARENT_LAYOUT)
     !  end for speed
     return
  elseif(i_node>=0) then

     if(.not.my_ORBIT_LATTICE%accel) return

     state=my_ORBIT_LATTICE%state
     my_ORBIT_LATTICE%state=my_ORBIT_LATTICE%state+totalpath0
     ! for speed
     !       t=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE
     !       DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%dpos
     !          t%parent_fibre%mag=my_ORBIT_LATTICE%state
     !          T=>T%NEXT
     !       ENDDO
     ! end for speed
     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE%parent_fibre
        if(.not.my_ORBIT_LATTICE%first) then
           !  write(6,*)" ############################################## "
           !  write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
           !  write(6,*)my_ORBIT_LATTICE%freqa,my_ORBIT_LATTICE%freqb
           !  write(6,*)my_ORBIT_LATTICE%volta,my_ORBIT_LATTICE%voltb
           !  write(6,*)my_ORBIT_LATTICE%phasa,my_ORBIT_LATTICE%phasb
           !  write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
           !  write(6,*)" ############################################## "
           my_ORBIT_LATTICE%freqb=p_orbit%mag%freq
           my_ORBIT_LATTICE%voltb=p_orbit%mag%volt
           my_ORBIT_LATTICE%phasb=p_orbit%mag%phas
           !  write(6,*)" ############################################## "
           !  write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
           !  write(6,*)my_ORBIT_LATTICE%freqa,my_ORBIT_LATTICE%freqb
           !  write(6,*)my_ORBIT_LATTICE%volta,my_ORBIT_LATTICE%voltb
           !  write(6,*)my_ORBIT_LATTICE%phasa,my_ORBIT_LATTICE%phasb
           !  write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
           !  write(6,*)" ############################################## "
        endif
        if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) then
           my_ORBIT_LATTICE%dxs6=x_orbit_sync(5)-my_ORBIT_LATTICE%xs6
           my_ORBIT_LATTICE%xs6=x_orbit_sync(5)
        else
           my_ORBIT_LATTICE%dxs6=x_orbit_sync(6)-my_ORBIT_LATTICE%xs6
           my_ORBIT_LATTICE%xs6=x_orbit_sync(6)
        endif

        call orbit_to_ptc(x_orbit_sync)
        !       p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE%parent_fibre
        my_ORBIT_LATTICE%orbit_deltae=zero
        w1_orbit=p_orbit
        w2_orbit=0
        call get_from_table_volt(x_orbit_sync(6),p0,vrfx,dphat)

        call find_energy(w2_orbit,p0c=p0)
        my_ORBIT_LATTICE%orbit_deltae=(w2_orbit%energy-w1_orbit%energy)
        freqf=w1_orbit%beta0*my_ORBIT_LATTICE%ORBIT_harmonic*clight/my_ORBIT_LATTICE%ORBIT_L
        p_orbit%mag%freq=freqf
        if(p_orbit%mag%l/=zero) then
           p_orbit%mag%volt=-vrfx/p_orbit%mag%l
           dphase0=asin(-my_ORBIT_LATTICE%orbit_deltae/p_orbit%DIR/p_orbit%parent_laYOUT%CHARGE &
                /p_orbit%mag%volt/c_1d_3/p_orbit%mag%L)-c_%phase0
        else
           p_orbit%mag%volt=-vrfx
           dphase0=asin(-my_ORBIT_LATTICE%orbit_deltae/p_orbit%DIR/p_orbit%parent_laYOUT%CHARGE/p_orbit%mag%volt/c_1d_3)-c_%phase0
        endif
        dphase=dphase0-twopi*p_orbit%mag%freq*x_orbit_sync(6)/CLIGHT
        p_orbit%mag%phas=dphase
        call ptc_to_orbit(x_orbit_sync)

     endif
     ! for speed
     CALL ORBIT_TRACK_NODE(i_node1,x_orbit_sync)

     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        if(.not.state%totalpath==1) then
           p_orbit%mag%phas=dphase0
           p_orbit%magp%phas=dphase0
        endif
        if(my_ORBIT_LATTICE%first) then
           my_ORBIT_LATTICE%freqb=p_orbit%mag%freq
           my_ORBIT_LATTICE%freqa=p_orbit%mag%freq
           my_ORBIT_LATTICE%voltb=p_orbit%mag%volt
           my_ORBIT_LATTICE%volta=p_orbit%mag%volt
           my_ORBIT_LATTICE%phasa=p_orbit%mag%phas
           my_ORBIT_LATTICE%phasb=p_orbit%mag%phas
           my_ORBIT_LATTICE%first=my_false
        else
           my_ORBIT_LATTICE%freqa=p_orbit%mag%freq
           my_ORBIT_LATTICE%volta=p_orbit%mag%volt
           my_ORBIT_LATTICE%phasa=p_orbit%mag%phas
        endif

        p_orbit%magp%freq=p_orbit%mag%freq
        p_orbit%magp%volt=p_orbit%mag%volt
        p_orbit%magp%phas=p_orbit%mag%phas
        if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) then
           dphat=x_orbit_sync(6)/my_ORBIT_LATTICE%ORBIT_P0C
        else
           dphat=x_orbit_sync(5)
        endif
        call find_energy(w2_orbit,ENERGY=dphat*w1_orbit%p0c+w1_orbit%energy)
        !  Changes the px,py,delta using the new p0C i.e. w2%p0c
        my_ORBIT_LATTICE%orbit_omega_after=twopi*p_orbit%mag%FREQ/CLIGHT
        CALL accel_ORBIT_up_grade_x(x_orbit_sync)
        !  upgrades all the reference energy from after cavity 1 to cavity 2 included
        call accel_ORBIT_up_grade_mag

     endif


     my_ORBIT_LATTICE%state=state
     ! for speed
     !       t=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE
     !       DO I=1,my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%dpos
     !          t%parent_fibre%mag=my_ORBIT_LATTICE%state
     !          T=>T%NEXT
     !       ENDDO
     ! end for speed
  elseif(i_node==-2) then
     if(.not.my_ORBIT_LATTICE%accel) return
     call kanalnummer(mf)
     open(unit=mf,file=initial_setting)
     read(MF,*)   my_ORBIT_LATTICE%freqb
     read(MF,*)   my_ORBIT_LATTICE%freqa
     read(MF,*)   my_ORBIT_LATTICE%voltb
     read(MF,*)   my_ORBIT_LATTICE%volta
     read(MF,*)   my_ORBIT_LATTICE%phasa
     read(MF,*)   my_ORBIT_LATTICE%phasb
     read(MF,*)   my_ORBIT_LATTICE%xs6
     read(MF,*)   my_ORBIT_LATTICE%dxs6
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_omega
     read(MF,*)   my_ORBIT_LATTICE%orbit_omega_after
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_gamma
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_P0C
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_BETA0
     read(MF,*)   my_ORBIT_LATTICE%orbit_kinetic
     read(MF,*)   my_ORBIT_LATTICE%orbit_energy
     read(MF,*)   my_ORBIT_LATTICE%orbit_brho
     read(MF,*)   my_ORBIT_LATTICE%orbit_dppfac
     do i=1,my_ORBIT_LATTICE%ORBIT_N_NODE
        if(my_ORBIT_LATTICE%ORBIT_NODES(i)%CAVITY) then
           p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i)%NODE%parent_fibre
           exit
        endif
     enddo
     read(MF,*)  p_orbit%mag%freq
     read(MF,*)  p_orbit%mag%volt
     read(MF,*)  p_orbit%mag%phas
     p_orbit%magp%freq=p_orbit%mag%freq
     p_orbit%magp%volt=p_orbit%mag%volt
     p_orbit%magp%phas=p_orbit%mag%phas
     read(MF,*) x_orbit_sync(1:2)
     read(MF,*) x_orbit_sync(3:4)
     read(MF,*) x_orbit_sync(5:6)
     read(MF,*) w1_ORBIT
     read(MF,*) w2_ORBIT
     read(MF,*) my_ORBIT_LATTICE%state
     my_ORBIT_LATTICE%first=.false.
     CLOSE(MF)
     call accel_ORBIT_up_grade_mag_all
     return
  elseif(i_node==-3) then
     stop 553
     call PUT_state(default,my_ORBIT_LATTICE%ORBIT_NODES(1)%node%PARENT_FIBRE%PARENT_LAYOUT)
     return
  elseif(i_node==-4) then
     stop 554
     return
  elseif(i_node==-5) then
     !    my_ORBIT_LATTICE%state=my_estate
     !    call print(my_ORBIT_LATTICE%state,6)
     stop 555
  endif


END SUBROUTINE  ptc_synchronous_set

SUBROUTINE ptc_synchronous_after(i_node)

  USE accel_ptc    !,vrff=>vrf,freqf=>freq
  IMPLICIT NONE
  INTEGER  i_node
  INTEGER  i_node1,mf,i
  if(.not.my_ORBIT_LATTICE%accel) return
  if(i_node>=0) then
     i_node1 = i_node + 1


     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE%parent_fibre
        p_orbit=w2_ORBIT

        my_ORBIT_LATTICE%ORBIT_OMEGA=my_ORBIT_LATTICE%orbit_omega_after
        my_ORBIT_LATTICE%ORBIT_gamma=1.d0/w2_ORBIT%gamma0i
        my_ORBIT_LATTICE%ORBIT_P0C=p_orbit%mag%P%P0C
        my_ORBIT_LATTICE%ORBIT_BETA0=p_orbit%mag%P%BETA0
        my_ORBIT_LATTICE%orbit_kinetic=w2_ORBIT%kinetic
        my_ORBIT_LATTICE%orbit_energy=w2_ORBIT%energy
        my_ORBIT_LATTICE%orbit_brho=w2_ORBIT%brho
        my_ORBIT_LATTICE%orbit_dppfac=one/sqrt(w2_orbit%beta0)/w2_orbit%energy
     ENDIF

  elseif(i_node==-2) then
     call kanalnummer(mf)
     open(unit=mf,file=final_setting)
     WRITE(MF,*)   my_ORBIT_LATTICE%freqb, "freqb"
     WRITE(MF,*)   my_ORBIT_LATTICE%freqa, "freqa"
     WRITE(MF,*)   my_ORBIT_LATTICE%voltb, "voltb"
     WRITE(MF,*)   my_ORBIT_LATTICE%volta, "volta"
     WRITE(MF,*)   my_ORBIT_LATTICE%phasa, "phasa"
     WRITE(MF,*)   my_ORBIT_LATTICE%phasb, "phasb"
     WRITE(MF,*)   my_ORBIT_LATTICE%xs6, "xs6"
     WRITE(MF,*)   my_ORBIT_LATTICE%dxs6, "dxs6"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_omega, "ORBIT_omega"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_omega_after, "orbit_omega_after"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_gamma, "ORBIT_gamma"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_P0C, "ORBIT_P0C"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_BETA0, "ORBIT_BETA0"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_kinetic, "orbit_kinetic"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_energy, "orbit_energy"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_brho, "orbit_brho"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_dppfac, "orbit_dppfac"
     WRITE(MF,*)  p_orbit%mag%freq
     WRITE(MF,*)  p_orbit%mag%volt
     WRITE(MF,*)  p_orbit%mag%phas
     WRITE(MF,*) x_orbit_sync(1:2)
     WRITE(MF,*) x_orbit_sync(3:4)
     WRITE(MF,*) x_orbit_sync(5:6)
     WRITE(MF,*) w1_ORBIT
     WRITE(MF,*) w2_ORBIT
     WRITE(MF,*) my_ORBIT_LATTICE%state
     CLOSE(MF)
     return
  endif

END SUBROUTINE  ptc_synchronous_after

subroutine ptc_track_particle(node_index, x,xp,y,yp,phi,dE)

  USE accel_ptc
  IMPLICIT NONE
  REAL(DP) x,xp,y,yp,phi,dE,x6
  INTEGER node_index
  INTEGER i

  i = node_index + 1

  call PUT_RAY(x,xp,y,yp,phi,dE)

  if(my_ORBIT_LATTICE%accel) then
     !  goto 111
     if(my_ORBIT_LATTICE%ORBIT_NODES(i)%CAVITY) then
        if(my_ORBIT_LATTICE%state%totalpath==1) then
           if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) then
              x6=(x_orbit(5)- my_ORBIT_LATTICE%xs6)/my_ORBIT_LATTICE%dxs6
           else
              x6=(x_orbit(6)-my_ORBIT_LATTICE%xs6)/my_ORBIT_LATTICE%dxs6
           endif
        else  ! relative path
           if(my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS) then
              x6=x_orbit(5)/my_ORBIT_LATTICE%dxs6
           else
              x6=x_orbit(6)/my_ORBIT_LATTICE%dxs6
           endif
           !  write(6,*)" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% "
           !  write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
           !  write(6,*)my_ORBIT_LATTICE%freqa,my_ORBIT_LATTICE%freqb
           !  write(6,*)my_ORBIT_LATTICE%volta,my_ORBIT_LATTICE%voltb
           !  write(6,*)my_ORBIT_LATTICE%phasa,my_ORBIT_LATTICE%phasb
           !  write(6,*)x_orbit(5),my_ORBIT_LATTICE%dxs6
           !  write(6,*)" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% "
           !  pause 12
        endif
        p_orbit%mag%freq=my_ORBIT_LATTICE%freqa+x6*(my_ORBIT_LATTICE%freqa-my_ORBIT_LATTICE%freqb)
        p_orbit%mag%volt=my_ORBIT_LATTICE%volta+x6*(my_ORBIT_LATTICE%volta-my_ORBIT_LATTICE%voltb)
        p_orbit%mag%phas=my_ORBIT_LATTICE%phasa+x6*(my_ORBIT_LATTICE%phasa-my_ORBIT_LATTICE%phasb)
        p_orbit%magp%freq=p_orbit%mag%freq
        p_orbit%magp%volt=p_orbit%mag%volt
        p_orbit%magp%phas=p_orbit%mag%phas

     endif  ! if cavity
111  continue
     call TRACK_ONE_NODE(i)

     if(my_ORBIT_LATTICE%ORBIT_NODES(i)%CAVITY) then
        CALL accel_ORBIT_up_grade_x(x_orbit)
     endif

     call GET_RAY(x,xp,y,yp,phi,dE)

  else
     call PUT_RAY(x,xp,y,yp,phi,dE)

     call TRACK_ONE_NODE(i)

     call GET_RAY(x,xp,y,yp,phi,dE)
  endif
  return
end subroutine ptc_track_particle
