module accel_ptc
  !  use beam_beam_ptc
  use orbit_ptc

  ! use orbit_ptc
  implicit none
  public
  integer n_cav_accel
  integer :: nvolt=20, NHARMON=1
  !  real(dp) :: vrf=210.d-3   !mev
  real(dp),parameter :: dphsa_final=30.0_dp*deg_to_rad_  ! 30.0_dp*deg_to_rad_
  real(dp) :: dphsa_by_dt=dphsa_final/(100.e-3_dp*CLIGHT)
  real(dp) :: harmon=nine
  real(dp) :: freq=zero
  real(dp) :: circum=zero
  real(dp) :: total_cav_L=zero
  logical(lp) :: enforce_zero_x5=my_false
  logical(lp) :: energy_reached=my_false
  real(dp) :: present_kinetic_energy=zero
  real(dp) :: final_kinetic_energy=40.0_dp  !3.1d0

  !   Other mode of acceleration
  real(dp) :: T_parabolic=100.e-3_dp*CLIGHT
  real(dp) :: T_rise=13.e-3_dp*CLIGHT,T_END_LINEAR=630091002.783574_dp,T_END=660070248.583574_dp
  real(dp) :: p0c0
  real(dp) :: p0c1 , kinetic1
  real(dp) :: p01_ratio = .181_dp/.143_dp
  real(dp) :: p0c2,kinetic2,p0cf
  real(dp) :: as1,as2
  real(dp) :: vrf0=40.d-3,vrff=280.d-3, slopev    !mev
  integer :: count_cav
  logical(lp) :: oldway= .false.
  integer :: slope_sign=1,slope_flip=1
  logicaL :: autoflip=.False.
  real(dp) :: maximum_phase=one
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
    integer k,i,mf,mode,t,j,imax
    integer, allocatable :: imode(:)
    character*120 line

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
    read(mf,'(a120)') line
    read(line,*) mode
    allocate(imode(mode))
    imode=0
    read(line,*) mode,imode
    imax=mode
    do i=1,mode
       if(imode(i)>imax) imax=imode(i)
    enddo

    NHARMON=IMAX
    !mode=1
    do i=0,k-2
       allocate(accel(i)%p0c)
       allocate(accel(i)%volt(imax))
       allocate(accel(i)%phase(imax))
       accel(i)%volt=zero
       accel(i)%phase=zero
       read(mf,* ) t,accel(i)%p0c,( accel(i)%volt(imode(j)),accel(i)%phase(imode(j)),j=1,mode)

       if(i==0) b0=accel(i)%p0c
       accel(i)%p0c=accel(i)%p0c/b0*p0c
       do j=1,  imax
          accel(i)%volt(j)=accel(i)%volt(j)*1.e-3_dp
       enddo

    enddo

    close(mf)
    deallocate(imode)
  end subroutine make_table

  subroutine get_from_table_volt(time,p0,vo,ph)
    implicit none
    integer ti,i
    real(dp) tr,time,p0,vo(nvolt),ph(nvolt)

    vo=zero
    ph=zero
    !write(6,*) size(accel)

    tr=time/clight*1000
    ti=int(tr)
    !write(6,*) tr,ti,accel(ti)%p0c

    if(ti>size(accel)) stop

    vo(1)= (tr-ti)*(accel(ti+1)%volt(1)-accel(ti)%volt(1)) + accel(ti)%volt(1)
    ph(1)= (tr-ti)*(accel(ti+1)%phase(1)-accel(ti)%phase(1)) + accel(ti)%phase(1)
    p0= (tr-ti)*(accel(ti+1)%p0c -accel(ti)%p0c) + accel(ti)%p0c


    do i=2,size(accel(1)%volt(:))
       vo(i)= (tr-ti)*(accel(ti+1)%volt(i)-accel(ti)%volt(i)) + accel(ti)%volt(i)
       ph(i)= (tr-ti)*(accel(ti+1)%phase(i)-accel(ti)%phase(i)) + accel(ti)%phase(i)
    enddo


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
    x(6)=x(6)*my_ORBIT_LATTICE%orbit_omega_after/my_ORBIT_LATTICE%ORBIT_OMEGA
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

  subroutine accel_orbit_beam(ring,filetable,npart)
    implicit none
    integer ipause, mypause
    integer n_turn,i,k,npart,j,kk
    type(layout), target :: ring
    TYPE(BEAM),target :: RAYS
    real(dp) sig0(6),x(6)
    integer mf
    type(fibre), pointer :: p
    character(*) filetable
    logical ooo

    sig0=1.e-6
    write(6,*) " old"
    read (5,*) ooo

    CALL create_beam(RAYS,npart,two,SIG0)

    my_ORBIT_LATTICE%ORBIT_USE_ORBIT_UNITS=.true.
    open(unit=mf, file='init.dat')
    !open(unit=mf, file='wp2_029_300kW_J54_INJALL_TRM_210kVM_Bf016.dat')
    do i=1,RAYS%n
       read(mf,*) RAYS%x(i,1:4), RAYS%x(i,5), RAYS%x(i,6)
    enddo
    close(mf)

    call make_table(filetable)

    !call make_table("NOACC_ACC_210.DAT")
    !call make_table("RF_Pattern_210kV_INJ120mc_ACC350ms.DAT")
    !call make_table("ACCWAVE_40kV_280kV_350ms.DAT")
    !call make_table("ACCWAVE_210KVH9_350ms.DAT")
    !call make_table("noaccel.DAT")

    write(6,*) "Read RF_table ... from RF_file: ", filetable(1:len_trim(filetable))
    ipause=mypause(321)

    write(6,*) my_ORBIT_LATTICE%ORBIT_harmonic
    write(6,*) my_ORBIT_LATTICE%ORBIT_omega
    default=my_ORBIT_LATTICE%state

    call kanalnummer(mf)
    open(unit=mf,file='outnew.dat')

    WRITE(6,*) " n_TURN ... pseudo ORBIT debugging"
    READ(5,*) N_TURN
    !n_turn=20

    !rays%x(2,1:6)=0.d0
    !rays%x(2,6)=0.00001d0

    if(ooo) then
       call ptc_synchronous_set_old(-1)
    else
       call ptc_synchronous_set(-1)
    endif
    write(6,*) " reading rays after some turns ?"
    read(5,*) i

    if(i==1) then

       call ptc_synchronous_set(-2)

       do i=1,RAYS%n
          read(mf,*) RAYS%x(i,1:4), RAYS%x(i,5), RAYS%x(i,6)
       enddo


       close(mf)

       call kanalnummer(mf)
       open(unit=mf,file='out2.dat')


    endif


    !call print(my_ORBIT_LATTICE%state,6)

    do j=1,n_turn

       do i=0,my_ORBIT_LATTICE%ORBIT_N_NODE-1

          if(ooo) then
             call ptc_synchronous_set_old(i)
          else
             call ptc_synchronous_set(i)
          endif
          do k=1,rays%n


             call ptc_track_particle(i,rays%x(k,1),rays%x(k,2),rays%x(k,3),rays%x(k,4),&
                  rays%x(k,5),rays%x(k,6))

          enddo

          call ptc_synchronous_after(i)

       enddo

       DO KK=1,RAYS%N
          X(5)=rays%x(kk,6)+my_ORBIT_LATTICE%ORBIT_P0C
          X(6)=x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*1000
          write(mf,'(1x,E25.17,1x,i8,1x,1(1x,E25.17))') x(6),j,x(5)     !rays%x(kk,1:6)
       ENDDO
    enddo  ! turn
    call ptc_synchronous_after(-2)

    !   DO KK=1,RAYS%N
    !      write(mf,'(6(1x,E25.17))') rays%x(kk,1:6)
    !   ENDDO
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


SUBROUTINE ptc_synchronous_set(i_node)

  USE accel_ptc     !,vrff=>vrf,freqf=>freq
  IMPLICIT NONE
  INTEGER  i_node
  INTEGER  i_node1,i,mf,j,nf
  type(internal_state) state
  real(dp) p0,vrfx,dphat,freqf,dt0,dt,x6,vo(nvolt),ph(nvolt)
  TYPE(INTEGRATION_NODE), POINTER  :: T


  i_node1 = i_node + 1

  if(i_node==-1) then
     x_orbit_sync=zero
     w1_orbit=0
     w2_orbit=0
     my_ORBIT_LATTICE%orbit_omega_after=my_ORBIT_LATTICE%ORBIT_omega
     my_ORBIT_LATTICE%state=my_ORBIT_LATTICE%state+time0
     default=my_ORBIT_LATTICE%state

     write(6,*) " Orbit set for acceleration "
     return
  elseif(i_node>=0) then

     if(.not.my_ORBIT_LATTICE%accel) return

     state=my_ORBIT_LATTICE%state
     my_ORBIT_LATTICE%state=my_ORBIT_LATTICE%state+totalpath0

     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE%parent_fibre


        call orbit_to_ptc(x_orbit_sync)

        my_ORBIT_LATTICE%orbit_deltae=zero
        w1_orbit=p_orbit
        w2_orbit=0
        !         write(6,*) "x_orbit_sync(6) ",x_orbit_sync(6)
        call get_from_table_volt(x_orbit_sync(6),p0,vo,ph)
        call find_energy(w2_orbit,p0c=p0)
        my_ORBIT_LATTICE%orbit_deltae=(w2_orbit%energy-w1_orbit%energy)
        freqf=w1_orbit%beta0*my_ORBIT_LATTICE%ORBIT_harmonic*clight/my_ORBIT_LATTICE%ORBIT_L
        p_orbit%mag%freq=freqf

        call compute_phase(x_orbit_sync,my_ORBIT_LATTICE%state,vo,ph,dt0)
        x6=x_orbit_sync(6)

        x_orbit_sync(5)=x_orbit_sync(5)+my_ORBIT_LATTICE%orbit_deltae/p_orbit%mag%p%p0c

        call ptc_to_orbit(x_orbit_sync)

     endif
     ! for speed
     call ORBIT_TRACK_NODE_fake(i_node1,x_orbit_sync)

     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        call adjust_phase(x6,state,dt0)

        p_orbit%magp%freq=p_orbit%mag%freq
        p_orbit%magp%volt=p_orbit%mag%volt
        p_orbit%magp%phas=p_orbit%mag%phas
        dphat=x_orbit_sync(6)/my_ORBIT_LATTICE%ORBIT_P0C

        call find_energy(w2_orbit,ENERGY=dphat*w1_orbit%p0c+w1_orbit%energy)
        !  Changes the px,py,delta using the new p0C i.e. w2%p0c
        my_ORBIT_LATTICE%orbit_omega_after=twopi*p_orbit%mag%FREQ/CLIGHT
        CALL accel_ORBIT_up_grade_x(x_orbit_sync)
        !  upgrades all the reference energy from after cavity 1 to cavity 2 included
        call accel_ORBIT_up_grade_mag

     endif


     my_ORBIT_LATTICE%state=state

  elseif(i_node==-2) then

     if(.not.my_ORBIT_LATTICE%accel) return
     call kanalnummer(mf)
     open(unit=mf,file=initial_setting)
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_omega
     read(MF,*)   my_ORBIT_LATTICE%orbit_omega_after
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_gamma
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_P0C
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_BETA0
     read(MF,*)   my_ORBIT_LATTICE%orbit_kinetic
     read(MF,*)   my_ORBIT_LATTICE%orbit_energy
     read(MF,*)   my_ORBIT_LATTICE%orbit_brho

     do i=1,my_ORBIT_LATTICE%ORBIT_N_NODE
        if(my_ORBIT_LATTICE%ORBIT_NODES(i)%CAVITY) then
           p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i)%NODE%parent_fibre
           exit
        endif
     enddo
     read(MF,*)  p_orbit%mag%freq
     read(MF,*)  p_orbit%mag%volt
     read(MF,*)  p_orbit%mag%phas
     read(MF,*)  p_orbit%mag%c4%CAVITY_TOTALPATH
     read(MF,*)  slope_sign
     read(MF,*)  slope_flip
     read(MF,*)  maximum_phase
     p_orbit%magp%freq=p_orbit%mag%freq
     p_orbit%magp%volt=p_orbit%mag%volt
     p_orbit%magp%phas=p_orbit%mag%phas
     p_orbit%magp%c4%CAVITY_TOTALPATH=p_orbit%mag%c4%CAVITY_TOTALPATH
     read(mf,*) nf
     if(size(p_orbit%mag%c4%ph)<nf) then
        write(6,*) " error not a big enough size for modes "
        stop 476
     endif
     do i=1,nf
        read(mf,*) j,p_orbit%mag%c4%f(i),p_orbit%mag%c4%ph(i)
        p_orbit%magp%c4%f(i)  = p_orbit%mag%c4%f(i)
        p_orbit%magp%c4%ph(i) = p_orbit%mag%c4%ph(i)
     enddo
     read(MF,*) x_orbit_sync(1:2)
     read(MF,*) x_orbit_sync(3:4)
     read(MF,*) x_orbit_sync(5:6)
     read(MF,*) w1_ORBIT
     read(MF,*) w2_ORBIT
     read(MF,*) my_ORBIT_LATTICE%state
     CLOSE(MF)
     call accel_ORBIT_up_grade_mag_all
     return
  elseif(i_node==-3) then
     stop 553
     !     call PUT_state(default,my_ORBIT_LATTICE%ORBIT_NODES(1)%node%PARENT_FIBRE%PARENT_LAYOUT)
     !     return
  elseif(i_node==-4) then
     stop 554
     return
  elseif(i_node==-5) then
     !    my_ORBIT_LATTICE%state=my_estate
     !    call print(my_ORBIT_LATTICE%state,6)
     stop 555
  endif


END SUBROUTINE  ptc_synchronous_set



!===========================================================
! This subroutine should be called before particle tracking.
! It tells PTC to do something related to acceleration.
!
!
!
!===========================================================
SUBROUTINE ptc_synchronous_set_old(i_node)

  USE accel_ptc     !,vrff=>vrf,freqf=>freq
  IMPLICIT NONE
  INTEGER  i_node
  INTEGER  i_node1,i,mf,j,nf
  type(internal_state) state
  real(dp) p0,vrfx,dphat,freqf,dt0,dt,x6,vo(nvolt),ph(nvolt)
  TYPE(INTEGRATION_NODE), POINTER  :: T


  i_node1 = i_node + 1

  if(i_node==-1) then
     x_orbit_sync=zero
     w1_orbit=0
     w2_orbit=0
     my_ORBIT_LATTICE%orbit_omega_after=my_ORBIT_LATTICE%ORBIT_omega
     my_ORBIT_LATTICE%state=my_ORBIT_LATTICE%state+time0
     default=my_ORBIT_LATTICE%state

     write(6,*) " Orbit set for acceleration "
     return
  elseif(i_node>=0) then

     if(.not.my_ORBIT_LATTICE%accel) return

     state=my_ORBIT_LATTICE%state
     my_ORBIT_LATTICE%state=my_ORBIT_LATTICE%state+totalpath0

     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%NODE%parent_fibre


        call orbit_to_ptc(x_orbit_sync)

        my_ORBIT_LATTICE%orbit_deltae=zero
        w1_orbit=p_orbit
        w2_orbit=0
        !         write(6,*) "x_orbit_sync(6) ",x_orbit_sync(6)
        call get_from_table_volt(x_orbit_sync(6),p0,vo,ph)
        call find_energy(w2_orbit,p0c=p0)
        my_ORBIT_LATTICE%orbit_deltae=(w2_orbit%energy-w1_orbit%energy)
        freqf=w1_orbit%beta0*my_ORBIT_LATTICE%ORBIT_harmonic*clight/my_ORBIT_LATTICE%ORBIT_L
        p_orbit%mag%freq=freqf

        call compute_phase(x_orbit_sync,my_ORBIT_LATTICE%state,vo,ph,dt0)
        x6=x_orbit_sync(6)


        call ptc_to_orbit(x_orbit_sync)

     endif
     ! for speed
     CALL ORBIT_TRACK_NODE(i_node1,x_orbit_sync)

     if(my_ORBIT_LATTICE%ORBIT_NODES(i_node1)%CAVITY) then
        call adjust_phase(x6,state,dt0)

        p_orbit%magp%freq=p_orbit%mag%freq
        p_orbit%magp%volt=p_orbit%mag%volt
        p_orbit%magp%phas=p_orbit%mag%phas
        dphat=x_orbit_sync(6)/my_ORBIT_LATTICE%ORBIT_P0C

        call find_energy(w2_orbit,ENERGY=dphat*w1_orbit%p0c+w1_orbit%energy)
        !  Changes the px,py,delta using the new p0C i.e. w2%p0c
        my_ORBIT_LATTICE%orbit_omega_after=twopi*p_orbit%mag%FREQ/CLIGHT
        CALL accel_ORBIT_up_grade_x(x_orbit_sync)
        !  upgrades all the reference energy from after cavity 1 to cavity 2 included
        call accel_ORBIT_up_grade_mag

     endif


     my_ORBIT_LATTICE%state=state

  elseif(i_node==-2) then

     if(.not.my_ORBIT_LATTICE%accel) return
     call kanalnummer(mf)
     open(unit=mf,file=initial_setting)
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_omega
     read(MF,*)   my_ORBIT_LATTICE%orbit_omega_after
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_gamma
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_P0C
     read(MF,*)   my_ORBIT_LATTICE%ORBIT_BETA0
     read(MF,*)   my_ORBIT_LATTICE%orbit_kinetic
     read(MF,*)   my_ORBIT_LATTICE%orbit_energy
     read(MF,*)   my_ORBIT_LATTICE%orbit_brho

     do i=1,my_ORBIT_LATTICE%ORBIT_N_NODE
        if(my_ORBIT_LATTICE%ORBIT_NODES(i)%CAVITY) then
           p_orbit=>my_ORBIT_LATTICE%ORBIT_NODES(i)%NODE%parent_fibre
           exit
        endif
     enddo
     read(MF,*)  p_orbit%mag%freq
     read(MF,*)  p_orbit%mag%volt
     read(MF,*)  p_orbit%mag%phas
     read(MF,*)  p_orbit%mag%c4%CAVITY_TOTALPATH
     read(MF,*)  slope_sign
     read(MF,*)  slope_flip
     read(MF,*)  maximum_phase
     p_orbit%magp%freq=p_orbit%mag%freq
     p_orbit%magp%volt=p_orbit%mag%volt
     p_orbit%magp%phas=p_orbit%mag%phas
     p_orbit%magp%c4%CAVITY_TOTALPATH=p_orbit%mag%c4%CAVITY_TOTALPATH
     read(mf,*) nf
     if(size(p_orbit%mag%c4%ph)<nf) then
        write(6,*) " error not a big enough size for modes "
        stop 476
     endif
     do i=1,nf
        read(mf,*) j,p_orbit%mag%c4%f(i),p_orbit%mag%c4%ph(i)
        p_orbit%magp%c4%f(i)  = p_orbit%mag%c4%f(i)
        p_orbit%magp%c4%ph(i) = p_orbit%mag%c4%ph(i)
     enddo
     read(MF,*) x_orbit_sync(1:2)
     read(MF,*) x_orbit_sync(3:4)
     read(MF,*) x_orbit_sync(5:6)
     read(MF,*) w1_ORBIT
     read(MF,*) w2_ORBIT
     read(MF,*) my_ORBIT_LATTICE%state
     CLOSE(MF)
     call accel_ORBIT_up_grade_mag_all
     return
  elseif(i_node==-3) then
     stop 553
     !     call PUT_state(default,my_ORBIT_LATTICE%ORBIT_NODES(1)%node%PARENT_FIBRE%PARENT_LAYOUT)
     !     return
  elseif(i_node==-4) then
     stop 554
     return
  elseif(i_node==-5) then
     !    my_ORBIT_LATTICE%state=my_estate
     !    call print(my_ORBIT_LATTICE%state,6)
     stop 555
  endif


END SUBROUTINE  ptc_synchronous_set_old



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
        my_ORBIT_LATTICE%ORBIT_gamma=one/w2_ORBIT%gamma0i
        my_ORBIT_LATTICE%ORBIT_P0C=p_orbit%mag%P%P0C
        my_ORBIT_LATTICE%ORBIT_BETA0=p_orbit%mag%P%BETA0
        my_ORBIT_LATTICE%orbit_kinetic=w2_ORBIT%kinetic
        my_ORBIT_LATTICE%orbit_energy=w2_ORBIT%energy
        my_ORBIT_LATTICE%orbit_brho=w2_ORBIT%brho
     ENDIF

  elseif(i_node==-2) then
     call kanalnummer(mf)
     open(unit=mf,file=final_setting)

     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_omega, "ORBIT_omega"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_omega_after, "orbit_omega_after"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_gamma, "ORBIT_gamma"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_P0C, "ORBIT_P0C"
     WRITE(MF,*)   my_ORBIT_LATTICE%ORBIT_BETA0, "ORBIT_BETA0"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_kinetic, "orbit_kinetic"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_energy, "orbit_energy"
     WRITE(MF,*)   my_ORBIT_LATTICE%orbit_brho, "orbit_brho"
     WRITE(MF,*)  p_orbit%mag%freq
     WRITE(MF,*)  p_orbit%mag%volt
     WRITE(MF,*)  p_orbit%mag%phas
     WRITE(MF,*)  p_orbit%mag%c4%CAVITY_TOTALPATH
     WRITE(MF,*)  slope_sign
     WRITE(MF,*)  slope_flip
     WRITE(MF,*)  maximum_phase
     write(mf,*) p_orbit%mag%c4%nf, " Modes "
     do i=1,p_orbit%mag%c4%nf
        write(mf,*) i,p_orbit%mag%c4%f(i),p_orbit%mag%c4%ph(i)
     enddo
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


SUBROUTINE compute_phase(x,state,v,ph,dt0)
  USE accel_ptc    !,vrff=>vrf,freqf=>freq
  IMPLICIT NONE
  real(dp) a,x(6),o,dl,driv,dt0,dtp0
  real(dp) v(nvolt),ph(nvolt),normb,norm
  integer n,i,k,j
  TYPE(CAV4),pointer :: EL
  TYPE(CAV4p),pointer :: ELp
  type(real_8) y(6)
  type(internal_state) state,local_state
  logical doit,conti
  local_state=state-totalpath0

  el=>p_orbit%mag%c4
  elp=>p_orbit%magp%c4
  n=size(p_orbit%mag%c4%f)
  if(n>nvolt) then
     write(6,*) " stopped in compute_phase "
     stop 546
  endif
  if(n<NHARMON) then
     write(6,*) "  NOT ENOUGH HARMONICS IN FLAT FILE"
     stop 547
  endif
  if(p_orbit%mag%l/=zero) then
     p_orbit%mag%volt=v(1)/p_orbit%mag%l*slope_flip
     p_orbit%magp%volt=v(1)/p_orbit%mag%l*slope_flip
     el%f(1)=one
     elp%f(1)=one
     el%ph(1)=ph(1)
     elp%ph(1)=ph(1)
     dl=p_orbit%mag%l
  else
     p_orbit%mag%volt=v(1)*slope_flip
     p_orbit%magp%volt=v(1)*slope_flip
     el%f(1)=one
     elp%f(1)=one
     el%ph(1)=ph(1)
     elp%ph(1)=ph(1)
     dl=one
  endif

  do i=2,n
     el%f(i)=v(i)/v(1)
     elp%f(i)=v(i)/v(1)
     el%ph(i)=ph(i)
     elp%ph(i)=ph(i)
  enddo

  O=twopi*p_orbit%mag%freq/CLIGHT
  el%phase0=zero
  elp%phase0=zero
  el%phas=zero
  elp%phas=zero
  el%t=zero
  elp%t=zero
  dt0=zero
  dtp0=zero
  norm=1.e38_dp
  normb=1.e38_dp
  doit=.true.
  conti=.true.
  k=0
  call init(1,1)
  call alloc(y)

  do while(conti.and.k<100)
     k=k+1
     !   a=zero
     !   driv=zero


     !     do i=1,n

     !      a=a-dl*el%f(i)*p_orbit%dir*p_orbit%charge*EL%volt*c_1d_3*sin(i*o*dt0+EL%PH(i)+EL%phas+EL%phase0)
     !      driv=driv-i*o*dl*el%f(i)*p_orbit%dir*p_orbit%charge*EL%volt*c_1d_3*cos(i*o*dt0+EL%PH(i)+EL%phas+EL%phase0)

     !   enddo
     !write(6,*) my_ORBIT_LATTICE%orbit_deltae,driv,o*dt0

     do j=1,5; y(j)=x(j); enddo;

        y(6)=dt0+(one.mono.1)
        call track(p_orbit,y,local_state)

        driv=y(5).sub.'1'
        !     write(6,*) (y(5).sub.'0'),driv,my_ORBIT_LATTICE%orbit_deltae,p_orbit%mag%p%p0c
        !      pause 777
        if(slope_sign*driv<0.and.autoflip) then  !!!!
           slope_flip=-slope_flip
           EL%volt=-EL%volt
           ELp%volt=-ELp%volt
           write(6,*) " flipped voltage "
        else  !!!!

           normb=abs(o*(dt0-dtp0))
           dtp0=dt0
           dt0=(my_ORBIT_LATTICE%orbit_deltae/p_orbit%mag%p%p0c-(y(5).sub.'0'))/(y(5).sub.'1')+dt0
           norm=abs(o*(dt0-dtp0))
           if(norm<1.e-9_dp.and.doit) then
              doit=.false.
           elseif(.not.doit)then
              if(normb<=norm) conti=.false.
           endif

        endif

     enddo

     !        write(6,*) "phase shift is  ",abs(o*(dt0))*rad_to_deg_ ," degrees ",my_ORBIT_LATTICE%orbit_deltae/p_orbit%mag%p%p0c
     call kill(y)
     if(k>100) then
        write(6,*) " did not converge in compute_phase "
        write(6,*) "phase shift is  ",abs(o*(dt0))*rad_to_deg_ ," degrees "
        norm=my_ORBIT_LATTICE%orbit_deltae/p_orbit%mag%p%p0c-(y(5).sub.'0')
        write(6,*) " Norm of DE ",norm
        write(6,*) " Normb  ",normb
        stop 100
     endif
     if(abs(o*(dt0))>maximum_phase) then
        write(6,*) "phase shift is huge ",abs(o*(dt0))*rad_to_deg_ ," degrees "
        norm=my_ORBIT_LATTICE%orbit_deltae/p_orbit%mag%p%p0c-(y(5).sub.'0')
        write(6,*) " Norm of DE ",norm
        write(6,*) " Normb  ",normb
        stop 101
     endif
     !      pause 888

     el%t=x(6)-dt0
     elp%t=x(6)-dt0


   end SUBROUTINE compute_phase


   SUBROUTINE adjust_phase(x6,state,dt0)
     USE accel_ptc    !,vrff=>vrf,freqf=>freq
     IMPLICIT NONE
     real(dp) x6,dt0
     integer n,i
     TYPE(CAV4),pointer :: EL
     TYPE(CAV4p),pointer :: ELp
     type(internal_state) state

     if(.not.state%totalpath==1) then
        p_orbit%mag%c4%t=-dt0
        p_orbit%magp%c4%t=-dt0
     endif
   end SUBROUTINE adjust_phase
   subroutine ptc_track_particle(node_index, x,xp,y,yp,phi,dE)

     USE accel_ptc
     IMPLICIT NONE
     REAL(DP) x,xp,y,yp,phi,dE,x6
     INTEGER node_index
     INTEGER i

     i = node_index + 1

     call PUT_RAY(x,xp,y,yp,phi,dE)

     if(my_ORBIT_LATTICE%accel) then

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
     IF(I==1.AND.MF_HERD/=0) THEN
        WRITE(MF_HERD,'(4(1X,E15.8))') PHI,DE,my_ORBIT_LATTICE%orbit_p0c, x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*c_1d3
        !       WRITE(MF_HERD,'(6(1X,E15.8))') PHI,DE,X_ORBIT(6),X_ORBIT(5), &
        !x_orbit_sync(5)/my_ORBIT_LATTICE%ORBIT_OMEGA/clight*1000.d0,my_ORBIT_LATTICE%ORBIT_OMEGA
        !        ,my_ORBIT_LATTICE%ORBIT_P0C
     ENDIF

     return
   end subroutine ptc_track_particle
